% The step function for MGRIT applied to the system P*e = r, where P
% represents the discretization of a linearized conservation law
% w_t + (alpha*w)_x = 0, where w is one of the characteristic variables of 
% an underlying nonlinear hyperbolic PDE system, i.e., the one with
% wave-speed alpha.
%
% This routine is called by "cons_law_st.m" when MGRIT is used as the
% solver for the linearized characteristic systems.
%
% All of the details of the PDE are contained in the object "my_cons_law_system"
% "char_idx" is index of the particular characteristic variable.
%
% This script is based on a more general script written for the MGRIT
% solution of nonlinear scalar conservation laws using high-order WENO
% discretizations.
%
%
% NOTE: This implementation is not particularly efficient in terms of
% storing things between outer iterations and for each of the different
% characteristic variables individually, etc.
function [w1, MGRIT_object] = step_cons_law_lin_char_MGRIT(w0, step_status, MGRIT_object, my_cons_law_system, char_idx)

    assert(strcmp(my_cons_law_system.pde_pa.bcs, 'periodic'), 'MGRIT solves only implemented for periodic BCs')

    level = step_status.level;

    if level == 1
        [w1, MGRIT_object] = step_MOL_linearized_wrapper(w0, step_status, MGRIT_object, my_cons_law_system, char_idx);
        
    elseif level > 1
        
        [w1, MGRIT_object] = step_SL_modified(w0, step_status, MGRIT_object, my_cons_law_system);    
    end
end


% This advances the linearized problem over the same interval of points
% that the fine-grid problem is advanced over during its relaxation. I.e.,
% this is the ideal coarse-grid operator of the linearized problem over
% whatever interval of points was coarsened on the fine grid.
function [w1, MGRIT_object] = step_MOL_linearized_wrapper(w0, step_status, MGRIT_object, my_cons_law_system, char_idx)

    assert(step_status.level == 1, 'This linearized problem only makes sense on level 1!');
    level = step_status.level;

    %% Step in the linearized problem.
    t0_idx = step_status.t0_idx;
    w1 = my_cons_law_system.step_linearized_characteristic_variable(t0_idx, w0, char_idx);
    
    
    %% Compute and store information relating to truncation error correction.
    % Unpack parameters 
    nx = my_cons_law_system.mesh_pa.nx;
    h  = my_cons_law_system.mesh_pa.h;
    dt = my_cons_law_system.mesh_pa.dt;
    
    %% Approximate departure points
    % Computing coarse-grid departure points by using fine-grid 
    % information then estimate fine-grid departure points by using a
    % forward Euler step.
    % Note we eval the wave-speed at all but the first interface, 
    % consistent with how departure points are computed on the coarse grid.
    arrive_points = my_cons_law_system.mesh_pa.x_interfaces(2:nx+1);

    t1_idx = t0_idx + 1;
    % The solution is not stored at the final time point, so if this is
    % the point where we need the wave-speed, then just use constant
    % interpolation back to the second last time point. 
    if t1_idx == my_cons_law_system.mesh_pa.nt; t1_idx = t0_idx; end


    %% get alpha
    my_cons_law_system.linearization_data.t0_idx = t1_idx;
    [alpha_arrive_neg, alpha_arrive_pos] = my_cons_law_system.linearized_wave_speed(char_idx);
    alpha_arrive_av  = 0.5*( alpha_arrive_neg + alpha_arrive_pos );
    alpha_arrive_av  = alpha_arrive_av(2:nx+1);
    depart_points    = arrive_points - dt*alpha_arrive_av;

    % Save departure points
    MGRIT_object.hierarchy(level).depart_points{t0_idx} = depart_points;
    
    
    %% Truncation error
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
    % interface. So, this means we evaluate the wave-speed a all interfaces
    % but skip the first one.
    %
    % Note that the same situation holds for where the wave-speed
    % enters into the error estimate: In each cell, we take a
    % difference between the wave-speed at the right interface and the
    % left interface.     
    
    % Unpack linearized wave-speed and remove first interface
    [alpha_neg, alpha_pos] = my_cons_law_system.linearized_wave_speed(char_idx);
    alpha_av = 0.5*( alpha_neg + alpha_pos );
    alpha_av = alpha_av(2:nx+1);

    % Linearized dissipation---this is not really what's used, but we just
    % approximate it with this.
    nu = abs(alpha_av);

    % Discretization is only ever 1st order
    % Order of ERK method
    qRK = 1;
    % Error constant
    eRK = -0.500000000000000;
    k = 1;
    % Mesh-independent constant that appears out the front of the FV error
    %eFV = -(-1)^k * factorial(k) * factorial(k-1) / factorial(2*k);
    eFV = 0.5;
    % The error terms associated with the FV disc and the RK disc.
    FV_error = -dt * h^(2*k-1) * eFV * nu    ;
    RK_error = -eRK *(dt * alpha_av).^(qRK+1);                          
    MGRIT_object.hierarchy(level).error{t0_idx} = RK_error + FV_error;
end


% Phi on coarser levels is a SL discretization with a dissipative
% correction to account for the dissipation of the fine-grid scheme.
function [w1, MGRIT_object] = step_SL_modified(w0, step_status, MGRIT_object, my_cons_law_system)

    level = step_status.level;
    
    assert(level > 1, 'SL step only makes sense on levels > 1!')
    
    t0_idx = step_status.t0_idx;
    t0     = step_status.t0;
    dt     = step_status.dt;

    t0_idx_fine = MGRIT_object.hierarchy(level-1).cpoint_indptr(t0_idx);
    t1_idx_fine = MGRIT_object.hierarchy(level-1).cpoint_indptr(t0_idx+1);
    m  = t1_idx_fine - t0_idx_fine;

    nx = my_cons_law_system.mesh_pa.nx;
    h  = my_cons_law_system.mesh_pa.h;
    poly_interp_degree = 1; % Discretization is only ever 1st order
    

    % Compute departure points associated with east boundaries only
    % (i.e., the first component in our numerical flux vector will be
    % associated with the second interface). This sets up the error
    % vector correctly below such that when its hit with the
    % first-order difference matrix D1 we get the first entry being the
    % difference of fluxes at the second and the first interfaces.
    arrive_points = my_cons_law_system.mesh_pa.x_interfaces(2:nx+1);

    depart_points = depart_point_lin_interp_vector(t0_idx_fine, m, arrive_points, MGRIT_object, level);

    % Save departure points
    MGRIT_object.hierarchy(level).depart_points{t0_idx} = depart_points;

    % Decompose departure point as required for SL step.
    [integral_translation, epsilon] = SL_decompose_FV_cell_evolution(arrive_points, depart_points, h);

    % Get stencil and discretization error for each departure point.
    [d, nodes, error_SL] = SL_get_stencil_and_error(epsilon, nx, poly_interp_degree);

    error_SL = error_SL / factorial(poly_interp_degree+1) * h^(poly_interp_degree+1);

    % Set handle to take SL step with this stencil (called below).
    my_SL_step = @(u0) SL_apply_stencil(u0, d, nodes, integral_translation, nx, poly_interp_degree);

    
    % Don't apply any BE correction, just take an SL step
    if strcmp(MGRIT_object.BE_coarse_solve.id, 'NONE')
        w1 = my_SL_step(w0);
        return
    end
    

    %% Assemble dissipative backward Euler matrix and it's LU factorization.

    % Store current disc error
    MGRIT_object.hierarchy(level).error{t0_idx} = error_SL;

    % Sum fine-grid error
    ideal_error = 0;
    for i = t0_idx_fine:t1_idx_fine-1
        ideal_error = ideal_error + MGRIT_object.hierarchy(level-1).error{i};
    end

    % Two-level coefficient
    gamma = error_SL - ideal_error;
    
    % Multilevel coefficient
    if level == 2
        nu = gamma;
    elseif level > 2
        nu = gamma;
        for i = t0_idx_fine:t1_idx_fine-1
            nu = nu + MGRIT_object.hierarchy(level-1).nu{i};
        end
    end

    % Store this so it can be accumulated on next level
    MGRIT_object.hierarchy(level).nu{t0_idx} = nu;


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
        %derivative = my_cons_law_system.disc_pa.spatial_order;
        derivative = 1;
        bias = 'U';
        approx_order = 1;
        [p_nodes, p_weights] = FD_stencil(derivative, bias, approx_order);
        p_weights = p_weights/(h^derivative);
        spatial_problem.stencil  = {{p_weights; p_nodes}};             
        MGRIT_object.D_space_T = get_toeplitz_spatial_disc(spatial_problem)'; % Take the transpose here 

        % D_time_order
        %derivative = my_cons_law_system.num_RK_stages;
        derivative = 1;
        bias = 'U';
        approx_order = 1;
        [p_nodes, p_weights] = FD_stencil(derivative, bias, approx_order);
        p_weights = p_weights/(h^derivative);
        spatial_problem.stencil  = {{p_weights; p_nodes}};             
        MGRIT_object.D_time = get_toeplitz_spatial_disc(spatial_problem); 

    end

    B = MGRIT_object.I ...
        + MGRIT_object.D1*(( nu ).*(MGRIT_object.D_space_T) );        
    % End of construction of backward Euler matrix.
    
    
    % Set handle for computing the action of Phi: MGRIT_object.hierarchy(level).Phi is a cell array
    % Use an LU factorization to solve BE system
    if strcmp(MGRIT_object.BE_coarse_solve.id, 'LU')
        %[B_lower, B_upper] = lu(B); % Factor B

        %MGRIT_object.hierarchy(level).SL_opers{t0_idx} = @(u0) B_upper\(B_lower\(my_SL_step(u0)));

        %u1 = B_upper\(B_lower\(my_SL_step(u0)));
        
        w1 = B \ (my_SL_step(w0));
        
    % Iteratively solve BE system with GMRES
    elseif strcmp(MGRIT_object.BE_coarse_solve.id, 'GMRES')

        restart = []; % This means no restarting, i.e., full-memory GMRES
%         MGRIT_object.hierarchy(level).SL_opers{t0_idx} = @(u0) ...
%             krylov_withoutoutput(@gmres, {B, my_SL_step(u0), ...
%                 restart, MGRIT_object.BE_coarse_solve.rtol, MGRIT_object.BE_coarse_solve.maxiter});

%         figure(1)
%         plot(abs(fft(my_SL_step(u0))))

        w1 = krylov_withoutoutput(@gmres, {B, my_SL_step(w0), ...
                 restart, MGRIT_object.BE_coarse_solve.rtol, MGRIT_object.BE_coarse_solve.maxiter});
    end

    
%     % Take the step!
%     u1 = MGRIT_object.hierarchy(level).SL_opers{t0_idx}(u0);
end

% This is just a wrapper function that allows us to call a Krylov solver
% and supress the output from printing to the console.
function x = krylov_withoutoutput(fun_handle, fun_inputs)
    [x, ~] = fun_handle(fun_inputs{:});
end