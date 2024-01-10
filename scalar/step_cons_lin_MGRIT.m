% The step function for MGRIT applied to the linear conservation law
%   e_t + (alpha(x, t)*e)_x = 0.
%
% This routine is called by "cons_lin_st_MGRIT.m"
%
% All of the details of the PDE are contained in the object "my_cons_law"
function [u1, MGRIT_object] = step_cons_lin_MGRIT(u0, step_status, MGRIT_object, my_cons_law)

    level = step_status.level;

    if level == 1
        [u1, MGRIT_object] = step_MOL(u0, step_status, MGRIT_object, my_cons_law);
        
    else
        [u1, MGRIT_object] = step_SL(u0, step_status, MGRIT_object, my_cons_law);
        
    end
end


% Level 1. This just takes a MOL step for the linear conservation law, and
% it includes some extra bits and pieces that we need to run MGRIT on the
% problem.
function [u1, MGRIT_object] = step_MOL(u0, step_status, MGRIT_object, my_cons_law)

    level  = step_status.level;
    t0_idx = step_status.t0_idx;
    t0     = step_status.t0;
    dt     = step_status.dt;

    % Take step
    u1 = my_cons_law.step(t0_idx, u0);
    
    %% Compute error terms to be used in truncation error correction if not done previously.
    if ~isfield(MGRIT_object.hierarchy(level), 'stepped') || isempty(MGRIT_object.hierarchy(level).stepped)
        MGRIT_object.hierarchy(level).stepped  = false(MGRIT_object.hierarchy(level).nt-1, 1);
    end


    if ~MGRIT_object.hierarchy(level).stepped(t0_idx)
    
        % Unpack parameters 
        k  = my_cons_law.reconstruction.k;
        p  = my_cons_law.reconstruction.p;
        q  = my_cons_law.num_RK_stages;
        nx = my_cons_law.mesh_pa.nx;
        h  = my_cons_law.mesh_pa.h;
        
        % Note we ignore the 1st interface for evaluating the wave-speed,
        % the LF dissipation, and the departure points. 
        x_interfaces = my_cons_law.mesh_pa.x_interfaces(2:nx+1);

        % Error constant
        %eRK = 0 - 1/factorial(qRK+1);
        if q == 1
            eRK = -0.500000000000000;
        elseif q == 3
            eRK = -0.041666666666667;
        end
        
        % Mesh-independent constant that appears out the front of the FV error
        %eFV = -(-1)^k * factorial(k) * factorial(k-1) / factorial(2*k);
        if k == 1
            eFV = 0.5;
        elseif k == 2
            eFV = -0.083333333333333;
        elseif k == 3
            eFV = 0.016666666666667;
        end

        % Choose a time to evaluate the wave-speed function at. I think
        % that perhaps by a Taylor series argument the collocation time in
        % the first stage is most important. (When I have fiddled with this
        % in the past it doesn't seem to have much of an impact which time
        % I choose though).
        %
        % The dissipation coefficient nu enters into the truncation
        % estimate as a difference. Specifically, for cell i, we have
        % a difference of nu_{i+1/2} - nu_{i-1/2}, and the way that this
        % difference is implemented is through the application of a
        % 1st-order upwind (wind blowing left -> right) finite-difference
        % matrix applied to a vector nu. So, to get the desired behaviour 
        % should have nu be a nx-dimensional vector, where its first entry
        % is the dissipation coefficient corresponding to the second
        % interface, so e.g., (D*nu)_1 = nu_2 - nu_{nx+1} = nu_2 - nu_1
        % because the (nx+1)-st interface is the same as the first
        % interface. So, this means we evaluate the wave-speed a at every
        % right interface (i.e., we skip the first one).
        %
        % Note that the same situation holds for where the wave-speed
        % enters into the error estimate: In each cell, we take a
        % difference between the wave-speed at the right interface and the
        % left interface. 
        
        alpha = my_cons_law.wave_speed(x_interfaces, t0);

        % Compute the numerical dissipation coefficient used in the 
        % Lax--Friedrichs flux.
        % Compute numerical dissipation
        if strcmp(my_cons_law.disc_pa.num_flux_id, 'GLF')
            nu = my_cons_law.f0_prime_max;
        elseif strcmp(my_cons_law.disc_pa.num_flux_id, 'LLF')
            nu = abs(alpha);
        end

        
        % The error terms associated with the FV disc and the RK disc.
        FV_error = dt * h^p * eFV * nu; 
        RK_error = eRK *(-dt * alpha).^(q+1);                          

        MGRIT_object.hierarchy(level).error_FV{t0_idx} = FV_error;
        MGRIT_object.hierarchy(level).error_RK{t0_idx} = RK_error;


        %% If computing coarse-grid departure points by using fine-grid 
        % information then estimate fine-grid departure points by using a
        % forward Euler step.
        if strcmp(MGRIT_object.depart_coarse_strategy, 'backtrack+interp')
            % Note we eval the wave-speed at all but the first interface, 
            % consistent with how departure points are computed on the coarse grid.
            arrive_points = x_interfaces;
            alpha_arrive  = my_cons_law.wave_speed(arrive_points, t0 + dt); 
            depart_points = arrive_points - dt*alpha_arrive;

            % Save departure points
            MGRIT_object.hierarchy(level).depart_points{t0_idx} = depart_points;
        end
    end
    % End of computation of truncation error terms
    
    
end
% step_MOL_linear_wrapper


% Phi on coarser levels is a SL discretization with a dissipative
% correction to account for the dissipation of the fine-grid scheme.
function [u1, MGRIT_object] = step_SL(u0, step_status, MGRIT_object, my_cons_law)

    level = step_status.level;
    t0_idx = step_status.t0_idx;
    t0     = step_status.t0;
    dt     = step_status.dt;


    if ~isfield(MGRIT_object.hierarchy(level), 'SL_opers') || isempty(MGRIT_object.hierarchy(level).SL_opers)
        MGRIT_object.hierarchy(level).SL_opers = cell(MGRIT_object.hierarchy(level).nt-1, 1);
        MGRIT_object.hierarchy(level).stepped  = false(MGRIT_object.hierarchy(level).nt-1, 1);
    end

    %% Haven't stepped here before so construct the truncation error 
    %% correction and the coarse-grid SL method
    if ~MGRIT_object.hierarchy(level).stepped(t0_idx)
        MGRIT_object.hierarchy(level).stepped(t0_idx) = true;

        t0_idx_fine = MGRIT_object.hierarchy(level-1).cpoint_indptr(t0_idx);
        t1_idx_fine = MGRIT_object.hierarchy(level-1).cpoint_indptr(t0_idx+1);
        m           = t1_idx_fine - t0_idx_fine;
        
        nx = my_cons_law.mesh_pa.nx;
        h  = my_cons_law.mesh_pa.h;
        poly_interp_degree = my_cons_law.reconstruction.p; % Note naming inconsistency here between p and k for SL problem. 

        % Compute departure points associated with east boundaries only
        % (i.e., the first component in our numerical flux vector will be
        % associated with the second interface). This sets up the error
        % vector correctly below such that when its hit with the
        % first-order difference matrix D1 we get the first entry being the
        % difference of fluxes at the second and the first interfaces.
        arrive_points = my_cons_law.mesh_pa.x_interfaces(2:nx+1);

        depart_points = depart_point_lin_interp_vector(t0_idx_fine, m, arrive_points, MGRIT_object, level);

        if strcmp(MGRIT_object.depart_coarse_strategy, 'backtrack+interp')            
            depart_points = depart_point_lin_interp_vector(t0_idx_fine, m, arrive_points, MGRIT_object, level);

        elseif strcmp(MGRIT_object.depart_coarse_strategy, 'finest-step-ERK')
            tspan = [t0+dt:-dt/m^(level-1):t0]; % Take m^(level-1) steps
            depart_points = MGRIT_object.depart_solver_ERK(tspan, arrive_points);

        elseif strcmp(MGRIT_object.depart_coarse_strategy, 'coarse-step-ERK')
            tspan = [t0+dt:-dt:t0]; % Take one step
            depart_points = MGRIT_object.depart_solver_ERK(tspan, arrive_points);

        end
        
        % Save departure points
        MGRIT_object.hierarchy(level).depart_points{t0_idx} = depart_points;

        % Decompose departure point as required for SL step.
        [integral_translation, epsilon] = SL_decompose_FV_cell_evolution(arrive_points, depart_points, h);

        % Get stencil and discretization error for each departure point.
        [d, nodes, error_SL_coarse] = SL_get_stencil_and_error(epsilon, nx, poly_interp_degree);

        error_SL_coarse = error_SL_coarse / factorial(poly_interp_degree+1) * h^(poly_interp_degree+1);
        
        % Set handle to take SL step with this stencil (called below).
        my_SL_step = @(u0) SL_apply_stencil(u0, d, nodes, integral_translation, nx, poly_interp_degree);


        %% Assemble dissipative backward Euler matrix 

        % Sum fine-grid truncation error coefficients across the given 
        % coarse interval. Note: We separate the FV and RK
        % contrubutions because they potentially are associated with
        % different difference matrices. 
        ideal_error_RK = 0;
        ideal_error_FV = 0;

        % To make the code more intelligible, the error arrays on
        % different levels are labelled slightly differently.
        % On first coarse level, accumulate the finest-grid errors.
        if level == 2
            for i = t0_idx_fine:t1_idx_fine-1
                ideal_error_FV = ideal_error_FV ...
                    + MGRIT_object.hierarchy(level-1).error_FV{i};

                ideal_error_RK = ideal_error_RK ...
                    + MGRIT_object.hierarchy(level-1).error_RK{i};
            end

        % On all coarser levels, accumulate the errors on the level
        % above us (which themselves are accumulations of errors above
        % them).
        else
            for i = t0_idx_fine:t1_idx_fine-1
                ideal_error_FV = ideal_error_FV ...
                    + MGRIT_object.hierarchy(level-1).ideal_error_FV{i};

                ideal_error_RK = ideal_error_RK ...
                    + MGRIT_object.hierarchy(level-1).ideal_error_RK{i};
            end
        end

        % Store accumlated errors such that they can be used on the
        % next level.
        MGRIT_object.hierarchy(level).ideal_error_FV{t0_idx} = ideal_error_FV;
        MGRIT_object.hierarchy(level).ideal_error_RK{t0_idx} = ideal_error_RK;


        % If an SL method was used on the fine level, accumulate the
        % associated truncation error coefficients so they can be added
        % into the correction for the rediscretized SL method.
        ideal_error_SL = 0;
        if level > 2
            for i = t0_idx_fine:t1_idx_fine-1
                ideal_error_SL = ideal_error_SL ...
                    + MGRIT_object.hierarchy(level-1).error_SL{i};
            end
        end

        % Compute the total SL error that needs to be included in the
        % correction (i.e., the coarse-grid SL error minus the ideal SL 
        % error). (note the signs are reversed here due to using a BE 
        % disc such that this means we add in the error from the ideal 
        % SL method, and we subtract out the coarse-grid SL method's error). 
        error_SL = error_SL_coarse - ideal_error_SL;
        % To make it more obvious what's happening here, we could
        % instead include in the BE matrix below a -ideal_error_SL and 
        % a +error_SL_coarse, but then on the next level remember that
        % we have to accumulate all of the terms in the m BE matrices 
        % from this level, so it's just easier to add these two things
        % now and store this.

        % Store this so it can be accumulated on next level
        MGRIT_object.hierarchy(level).error_SL{t0_idx} = error_SL;

        % Build and save the spatial discretization operators.
        if ~isfield(MGRIT_object, 'I') || size(MGRIT_object.I, 1) ~= nx
            MGRIT_object.I  = speye(nx);

            % D1: First-order derivative
            derivative = 1;
            bias = 'U';
            approx_order = 1;
            [p_nodes, p_weights] = FD_stencil(derivative, bias, approx_order);
            p_weights = p_weights/(h^derivative);
            spatial_problem.stencil  = {{p_weights; p_nodes}};             
            spatial_problem.nDOFs = nx;
            spatial_problem.BCs.type = 'periodic';
            MGRIT_object.D1 = get_toeplitz_spatial_disc(spatial_problem); 

            % D_space_order
            derivative = my_cons_law.disc_pa.spatial_order;
            bias = 'U';
            approx_order = 1;
            [p_nodes, p_weights] = FD_stencil(derivative, bias, approx_order);
            p_weights = p_weights/(h^derivative);
            spatial_problem.stencil  = {{p_weights; p_nodes}};             
            MGRIT_object.D_space = get_toeplitz_spatial_disc(spatial_problem); 

            % D_time_order
            derivative = my_cons_law.num_RK_stages;
            bias = 'U';
            approx_order = 1;
            [p_nodes, p_weights] = FD_stencil(derivative, bias, approx_order);
            p_weights = p_weights/(h^derivative);
            spatial_problem.stencil  = {{p_weights; p_nodes}};             
            MGRIT_object.D_time = get_toeplitz_spatial_disc(spatial_problem); 

        end

        % Subtract out error contributions from coarse-grid rediscretized
        % coarse-grid SL method and add in error contributions from m
        % steps of the fine level method (note the signs are opposite
        % here due to using a BE disc).
        B = MGRIT_object.I ...
            + MGRIT_object.D1*(...
            ( ideal_error_FV + error_SL ) .* (MGRIT_object.D_space') +...
            ( ideal_error_RK            ) .* (MGRIT_object.D_time')...
            );
             
        
        %% Assemble solver for error correction matrix B.
        
        % Set handle for computing the action of Phi: MGRIT_object.hierarchy(level).Phi is a cell array
        % Use an LU factorization to solve BE system
        if strcmp(MGRIT_object.BE_coarse_solve.id, 'LU')
            [B_lower, B_upper] = lu(B); % Factor B
    
            MGRIT_object.hierarchy(level).SL_opers{t0_idx} = @(u0) B_upper\(B_lower\(my_SL_step(u0)));
        
        % Iteratively solve BE system with GMRES
        elseif strcmp(MGRIT_object.BE_coarse_solve.id, 'GMRES')
        
            restart = []; % This means no restarting, i.e., full-memory GMRES
            MGRIT_object.hierarchy(level).SL_opers{t0_idx} = @(u0) ...
                krylov_withoutoutput(@gmres, {B, my_SL_step(u0), ...
                    restart, MGRIT_object.BE_coarse_solve.rtol, MGRIT_object.BE_coarse_solve.maxiter});
        
        % Don't solve BE system at all. Just take an SL step (of course we
        % don't need all the code above to create B then, but just to keep
        % the code readable, just leave it, especially since this is an
        % edge case).
        elseif strcmp(MGRIT_object.BE_coarse_solve.id, 'NONE')
            
            MGRIT_object.hierarchy(level).SL_opers{t0_idx} = @(u0) my_SL_step(u0);
        end

    end
    
    % Take the step!
    u1 = MGRIT_object.hierarchy(level).SL_opers{t0_idx}(u0);
end

% This is just a wrapper function that allows us to call a Krylov solver
% and supress the output from printing to the console.
function x = krylov_withoutoutput(fun_handle, fun_inputs)
    [x, ~] = fun_handle(fun_inputs{:});
end