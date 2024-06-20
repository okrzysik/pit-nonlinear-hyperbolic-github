% Implement the "step" function for MGRIT for a Godunov discretization of
% the advection equation:
%
%   w_t \pm c0*w_x = 0.
%
% Where it is asserted that c0 > 0. The wave-speed sign (as in the \pm in 
% the equation) is expected as a field in disc_pa as disc_pa.c0_sign
% A vector of c0 > 0 discretized on the mesh must also appear as a field in
% disc_pa as disc_pa.c_cc
function [w1, MGRIT_object] = step_advection_godunov_MGRIT(w0, step_status, MGRIT_object, mesh_pa, disc_pa)

    level  = step_status.level;
    t0_idx = step_status.t0_idx;
    dt     = step_status.dt;
    h      = mesh_pa.h;
    nx     = mesh_pa.nx;

    if ~isfield(MGRIT_object, 'hierarchy')
        MGRIT_object.hierarchy = struct();
    end
    
        
    % Does the time-stepping operator for this level need to be assembled?
    % Note that on each level these are time independent, so they're just
    % built once per level.
    if ~isfield(MGRIT_object.hierarchy(level), 'Phi') || isempty(MGRIT_object.hierarchy(level).Phi)
        
        %% Assemble Phi for level 1.
        % Phi is a FE+U1 discretization of w_t + c0*w_x = 0
        if level == 1
                        
            % Build and the spatial discretization operators.
            I  = speye(nx);

            % Build D1: First-order upwind derivative
            derivative = 1;
            bias = 'U';
            approx_order = 1;
            [p_nodes, p_weights] = FD_stencil(derivative, bias, approx_order);
            p_weights = p_weights/(h^derivative);
            spatial_problem.stencil  = {{p_weights; p_nodes}};             
            spatial_problem.nDOFs = nx;
            spatial_problem.BCs.type = 'periodic';
            L = get_toeplitz_spatial_disc(spatial_problem); 
            
            % Assemble discretization as a matrix.
            if disc_pa.c0_sign == -1
                Phi_matrix = I - dt * (-disc_pa.c_cc) .* (-L.');
            elseif disc_pa.c0_sign == 1
                Phi_matrix = I - dt * (+disc_pa.c_cc) .* (+L);
            end
            
            MGRIT_object.hierarchy(level).Phi = @(b) Phi_matrix*b;
            
            
            %% Compute additional information required for stepping on coarse levels
            % Truncation error. This is the same for Phi_1 and Phi_2
            qRK = 1;
            eRK = -0.5;
            eFV =  0.5;
            FV_error = dt * h * eFV * disc_pa.c_cc;
            RK_error = eRK *(dt * disc_pa.c_cc).^(qRK+1) ;                          
            MGRIT_object.hierarchy(level).error = -(RK_error + FV_error);
            
            % Departure points             
            arrive_points = mesh_pa.x_centers;
            if disc_pa.c0_sign == -1
                MGRIT_object.hierarchy(level).depart_points = arrive_points - dt * -disc_pa.c_cc;
            elseif disc_pa.c0_sign == 1
                MGRIT_object.hierarchy(level).depart_points = arrive_points - dt * +disc_pa.c_cc;
            end

        % End of building Phi for level == 1
        
        %% Assemble Phi for all levels >= 2
        % Phi is a modified SL discretization of w_t \pm c0*w_x = 0
        elseif level >= 2
            
            t0_idx_fine = MGRIT_object.hierarchy(level-1).cpoint_indptr(t0_idx);
            t1_idx_fine = MGRIT_object.hierarchy(level-1).cpoint_indptr(t0_idx+1);
            m           = t1_idx_fine - t0_idx_fine;
            
            % Compute depature points
            MGRIT_object.hierarchy(level).depart_points = depart_point_lin_interp_steady(m, mesh_pa.x_centers, MGRIT_object.hierarchy(level-1).depart_points);
            
            poly_interp_degree = 1;
            xmin = mesh_pa.x_centers(1);
            
            % Build matrix that interpolates at the feet of characteristics with slope \pm c0
            temp = (MGRIT_object.hierarchy(level).depart_points - xmin)/h;
            depart_east_neighbour = ceil(temp)+1;
            depart_east_epsilon = depart_east_neighbour - temp-1;
            [SL_matrix, SL_error] = semiLagrangianMatrix_new(poly_interp_degree, nx, depart_east_neighbour, depart_east_epsilon);  
            
            
            %% Build truncation error correction matrices.
            I  = speye(nx);

            % D2: Second-order second derivative
            derivative = 2;
            bias = 'C';
            approx_order = 2;
            [p_nodes, p_weights] = FD_stencil(derivative, bias, approx_order);
            p_weights = p_weights/(h^derivative);
            spatial_problem.stencil  = {{p_weights; p_nodes}};             
            spatial_problem.nDOFs = nx;
            spatial_problem.BCs.type = 'periodic';
            D2 = get_toeplitz_spatial_disc(spatial_problem); 
            
            % Accumulate error from m steps on the fine grid
            ideal_error = m * MGRIT_object.hierarchy(level-1).error;
            
            % Error from one step on current grid
            direct_error = SL_error / factorial(poly_interp_degree+1) * h^(poly_interp_degree+1);
            
            % Combine to get total error on current grid
            gamma = -ideal_error + direct_error;
            
            MGRIT_object.hierarchy(level).error = -gamma;

            BE_matrix = (I - gamma .* D2);
            

            % How are the BE matrices inverted?
            % Factor B
            if strcmp(MGRIT_object.BE_coarse_solve.id, 'LU')
                [BE_matrix_lower, BE_matrix_upper] = lu(BE_matrix); 
                BE_matrix_inv = @(b) BE_matrix_upper\(BE_matrix_lower \ b);
                
            % Iteratively solve BE system with GMRES
            elseif strcmp(MGRIT_object.BE_coarse_solve.id, 'GMRES')
                
                restart = []; % This means no restarting, i.e., full-memory GMRES
                BE_matrix_inv = @(b) krylov_withoutoutput(@gmres, {BE_matrix, b, ...
                        restart, MGRIT_object.BE_coarse_solve.rtol, MGRIT_object.BE_coarse_solve.maxiter});

            % Don't solve BE system at all. Just take an SL step (of course we
            % don't need all the code above to create B then, but just to keep
            % the code readable, just leave it, especially since this is an
            % edge case).
            elseif strcmp(MGRIT_object.BE_coarse_solve.id, 'NONE')
                
                BE_matrix_inv = @(b) b;
            end
            
            MGRIT_object.hierarchy(level).Phi = @(b) BE_matrix_inv(SL_matrix * b);
            
        end
        % End of building Phi for level >= 2
        
    end
    % End of building Phi

    % Apply the step: Evolve w0 into w1.
    w1 = MGRIT_object.hierarchy(level).Phi(w0);

end

% This is just a wrapper function that allows us to call a Krylov solver
% and supress the output from printing to the console.
function x = krylov_withoutoutput(fun_handle, fun_inputs)
    [x, ~] = fun_handle(fun_inputs{:});
end