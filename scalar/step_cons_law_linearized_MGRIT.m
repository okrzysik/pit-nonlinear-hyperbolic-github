% The step function for MGRIT applied to the system P*e = r, where P
% represents the discretization of a linearized conservation law.
%
% This routine is called by "cons_law_st.m" when MGRIT is used as the
% solver for the linearized systems.
%
% All of the details of the PDE are contained in the object "my_cons_law"
function [u1, MGRIT_object] = step_cons_law_linearized_MGRIT(u0, step_status, MGRIT_object, my_cons_law)

    level = step_status.level;

    if level == 1
        [u1, MGRIT_object] = step_MOL_linearized_wrapper(u0, step_status, MGRIT_object, my_cons_law);
        
    elseif level > 1
        
        [u1, MGRIT_object] = step_SL_modified(u0, step_status, MGRIT_object, my_cons_law);    
    end
end


%% Options for applying trunc. correction
function correction = truncation_correction_options(k)
    
    correction.use_linear_weights_in_time = ~true; 
    correction.use_linear_weights = ~true; % If using a WENO disc, approx. trunc error using linear weights or WENO weights...

    % For 1st order
    if k == 1
        correction.diss_idx   = 0;
        correction.wave_idx   = 0;
        correction.weight_idx = 0;
        correction.fudge_FV   = 1;
        correction.fudge_RK   = 1;
        
        correction.weighted_stencil = ~true; % Approximate truncation error on a big stencil (i.e. a non-weighted stencil).

    % For 3rd order
    elseif k == 2
        correction.diss_idx   = 1;
        correction.wave_idx   = 1;
        correction.weight_idx = 1;
        correction.fudge_FV   = 4/3;
        correction.fudge_RK   = 4*4/3;
        
        correction.weighted_stencil = true; % Approximate truncation error on a weighted stencil.
    end
end




% This advances the linearized problem over the same interval of points
% that the fine-grid problem is advanced over during its relaxation. I.e.,
% this is the ideal coarse-grid operator of the linearized problem over
% whatever interval of points was coarsened on the fine grid.
function [u1, MGRIT_object] = step_MOL_linearized_wrapper(u0, step_status, MGRIT_object, my_cons_law)

    assert(step_status.level == 1, 'This linearized problem only makes sense on level 1!');
    level = step_status.level;

    %% Step in the linearized problem.
    t0_idx = step_status.t0_idx;
    u1 = my_cons_law.step_linearized(t0_idx, u0); 
    
    
    %% Compute and store information relating to truncation error correction.
    % Unpack parameters 
    nx = my_cons_law.mesh_pa.nx;
    h  = my_cons_law.mesh_pa.h;
    dt = my_cons_law.mesh_pa.dt;
    k  = my_cons_law.reconstruction.k; 
    correction = truncation_correction_options(k);
    
    %% Approximate departure points
    % Computing coarse-grid departure points by using fine-grid 
    % information then estimate fine-grid departure points by using a
    % forward Euler step.
    % Note we eval the wave-speed at all but the first interface, 
    % consistent with how departure points are computed on the coarse grid.
    arrive_points = my_cons_law.mesh_pa.x_interfaces(2:nx+1);

    t1_idx = t0_idx + 1;
    % The solution is not stored at the final time point, so if this is
    % the point where we need the wave-speed, then just use constant
    % interpolation back to the second last time point. 
    if t1_idx == my_cons_law.mesh_pa.nt; t1_idx = t0_idx; end

    stage_idx = 1;
    alpha_arrive_neg = my_cons_law.linearization_data.wave_speed{t1_idx}{stage_idx}.neg;
    alpha_arrive_pos = my_cons_law.linearization_data.wave_speed{t1_idx}{stage_idx}.pos;
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
    switch correction.diss_idx
        case 0
            diss_idx = t0_idx; 
        case 1
            diss_idx = t1_idx; 
    end
    
    switch correction.wave_idx
        case 0
            wave_idx = t0_idx; 
        case 1
            wave_idx = t1_idx; 
    end
    
    
    % Unpack linearized dissipation and remove first interface
    stage_idx = 1;
    nu = my_cons_law.linearization_data.dissipation{diss_idx}{stage_idx};
    if numel(nu) > 1; nu = nu(2:nx+1); end % For GLF, nu is just a scalar.

    % Unpack linearized wave-speed and remove first interface
    alpha_neg = my_cons_law.linearization_data.wave_speed{wave_idx}{stage_idx}.neg;
    alpha_pos = my_cons_law.linearization_data.wave_speed{wave_idx}{stage_idx}.pos;
    alpha_av = 0.5*( alpha_neg + alpha_pos );
    %alpha_av = max( alpha_neg, alpha_pos );
    alpha_av = alpha_av(2:nx+1);

    % Order of ERK method
    qRK = my_cons_law.num_RK_stages;
    % Error constant
    %eRK = 0 - 1/factorial(qRK+1);
    if qRK == 1
        eRK = -0.500000000000000;
    elseif qRK == 3
        eRK = -0.041666666666667;
    else
        error('Implementation assumes number of RK stages is 1 or 3')
    end
    
    % Use a single big stencil for the correction
    if ~correction.weighted_stencil
        % Mesh-independent constant that appears out the front of the FV error
        %eFV = -(-1)^k * factorial(k) * factorial(k-1) / factorial(2*k);
        if k == 1
            eFV = 0.5;
        elseif k == 2
            eFV = -0.083333333333333; % == -1/12
        elseif k == 3
            eFV = 0.016666666666667;
        end

        % The error terms associated with the FV disc and the RK disc.
        FV_error = -dt * h^(2*k-1) * eFV * nu     * correction.fudge_FV;
        RK_error = -eRK *(dt * alpha_av).^(qRK+1) * correction.fudge_RK;                          
        
        assert(2*k-1 == qRK, 'implementation assumes these are equal');
        MGRIT_object.hierarchy(level).error{t0_idx} = RK_error + FV_error;
        
        
    % Estimate truncation error using a weighted stencil.
    elseif correction.weighted_stencil
        
        switch correction.weight_idx 
            case 0
                weight_idx = t0_idx; 
            case 1
                weight_idx = t1_idx; 
        end
        
        
        if k ~= 2; error('this correction assumes k == 2!'); end
        
        % Coefficients in weighted error estimates
        zeta0       = -1;
        zeta1       =  2;
        zeta0_tilde =  2;
        zeta1_tilde = -1;
    
        % Use linear weights for approximating truncation error.
        if strcmp(my_cons_law.disc_pa.reconstruction_id, 'linear') || ...
            correction.use_linear_weights
            
            % The optimal linear weights.
            b0          = 2/3; % ell = 0, RHS interface reconstruction
            b1          = 1/3; % ell = 1, RHS interface reconstruction 
            b0_tilde    = 1/3; % ell = 0, LHS interface reconstruction
            b1_tilde    = 2/3; % ell = 1, LHS interface reconstruction
        
        % Use WENO weights for approximating truncation error
        else

            stage_idx = 1;            
            % u_{i+1/2}^{-} WENO weights
            b0       = my_cons_law.linearization_data.reconstruction_weights{weight_idx}{stage_idx}.right(:, 1); % ell = 0
            b1       = my_cons_law.linearization_data.reconstruction_weights{weight_idx}{stage_idx}.right(:, 2); % ell = 1

            % u_{i+1/2}^{+} WENO weights
            b0_tilde = my_cons_law.linearization_data.reconstruction_weights{weight_idx}{stage_idx}.left(:, 1); % ell = 0
            b1_tilde = my_cons_law.linearization_data.reconstruction_weights{weight_idx}{stage_idx}.left(:, 2); % ell = 1

            % Shift so the first entries are the LHS weights of the second
            % cell, corresponding to what's fed into the numerical flux for the
            % second interface. This allows us to more easily implement the
            % required w_{i+1} - w_{i}
            b0_tilde = [b0_tilde(2:nx); b0_tilde(1)];
            b1_tilde = [b1_tilde(2:nx); b1_tilde(1)];

        end
        % End of picking weights b, t_tilde.
        
        % Constant. There is dt * 0.5 * nu, then the 3/4 * (h^2 / 6) is the
        % constant out the front of estimate of the reconstruction
        % differences.
        FV_error = dt * 0.5 * nu;
        c = FV_error * 3/4 * (h^2 / 6);

        c = c * correction.fudge_FV;
        
        MGRIT_object.hierarchy(level).error0{t0_idx} = c .* (b0_tilde * zeta0_tilde - b0 * zeta0);
        MGRIT_object.hierarchy(level).error1{t0_idx} = c .* (b1_tilde * zeta1_tilde - b1 * zeta1);

        
        % For powers of the spatial disc, do we want to use a weighted,
        % WENO-aware stencil or not? Selecting the linear weights below
        % will just mean its the same as doing it on a big stencil.
        %
        % Burgers Riemann problem: LLF: Doesn't make any difference. (this also seems true for Newton iteration.)  
        %
        % Burgers Riemann problem: GLF: Marginally better (maybe) at
        % larger mesh resolutions if WENO weights are used here rather than 
        % linear weights. (doesn't help for Newton)
        %
        % BL(3): GLF: Seems like it's possibly marginally better to use
        % WENO weights rather than linear weights.
        %
        % BL(3): LLF: Seems like also possibly marginally better to use
        % WENO weights rather than linear weights.
        
        % Choosing these weights means the RK error is the same as if
        % discretized on a big stencil without knowing about the WENO
        % weights.
        if correction.use_linear_weights_in_time
            b0          = 2/3;
            b1          = 1/3;
            b0_tilde    = 1/3;
            b1_tilde    = 2/3;     
        end
    
        RK_error = -eRK *(dt * alpha_av).^(qRK+1);       
        c = RK_error * 3/4 / h; % Modify constant so it can be split across a weighted stencil.

        c = c * correction.fudge_RK;
        
        MGRIT_object.hierarchy(level).error0{t0_idx} = MGRIT_object.hierarchy(level).error0{t0_idx} + c .* (b0_tilde * zeta0_tilde - b0 * zeta0);
        MGRIT_object.hierarchy(level).error1{t0_idx} = MGRIT_object.hierarchy(level).error1{t0_idx} + c .* (b1_tilde * zeta1_tilde - b1 * zeta1);
        
    end
    % End of computing truncation error on weighted stencils

end


% Phi on coarser levels is a SL discretization with a dissipative
% correction to account for the dissipation of the fine-grid scheme.
function [u1, MGRIT_object] = step_SL_modified(u0, step_status, MGRIT_object, my_cons_law)

    level = step_status.level;
    
    assert(level > 1, 'SL step only makes sense on levels > 1!')
    
    t0_idx = step_status.t0_idx;
    t0     = step_status.t0;
    dt     = step_status.dt;

    t0_idx_fine = MGRIT_object.hierarchy(level-1).cpoint_indptr(t0_idx);
    t1_idx_fine = MGRIT_object.hierarchy(level-1).cpoint_indptr(t0_idx+1);
    m  = t1_idx_fine - t0_idx_fine;

    nx = my_cons_law.mesh_pa.nx;
    h  = my_cons_law.mesh_pa.h;
    poly_interp_degree = my_cons_law.reconstruction.p; % Note naming inconsistency here between p and k for SL problem. 
    
    k  = my_cons_law.reconstruction.k; 
    correction = truncation_correction_options(k);

    % Compute departure points associated with east boundaries only
    % (i.e., the first component in our numerical flux vector will be
    % associated with the second interface). This sets up the error
    % vector correctly below such that when its hit with the
    % first-order difference matrix D1 we get the first entry being the
    % difference of fluxes at the second and the first interfaces.
    arrive_points = my_cons_law.mesh_pa.x_interfaces(2:nx+1);

    depart_points = depart_point_lin_interp_vector(t0_idx_fine, m, arrive_points, MGRIT_object, level);

    % Save departure points
    MGRIT_object.hierarchy(level).depart_points{t0_idx} = depart_points;

    % Decompose departure point as required for SL step.
    [integral_translation, epsilon] = SL_decompose_FV_cell_evolution(arrive_points, depart_points, h);

    % Get stencil and discretization error for each departure point.
    [d, nodes, error_SL_coarse] = SL_get_stencil_and_error(epsilon, nx, poly_interp_degree);

    error_SL_coarse = error_SL_coarse / factorial(poly_interp_degree+1) * h^(poly_interp_degree+1);

    % Set handle to take SL step with this stencil (called below).
    my_SL_step = @(u0) SL_apply_stencil(u0, d, nodes, integral_translation, nx, poly_interp_degree);

    
    % Don't apply any BE correction, just take an SL step
    if strcmp(MGRIT_object.BE_coarse_solve.id, 'NONE')
        u1 = my_SL_step(u0);
        return
    end
    

    %% Assemble dissipative backward Euler matrix and it's LU factorization.
    if ~correction.weighted_stencil
        % Sum fine-grid truncation error coefficients across the given 
        % coarse interval. Note: We separate the FV and RK
        % contrubutions because they potentially are associated with
        % different difference matrices. 

        assert(my_cons_law.disc_pa.spatial_order == my_cons_law.num_RK_stages, 'implementation assumes these are equal');
        
        ideal_error = 0;
        for i = t0_idx_fine:t1_idx_fine-1
                ideal_error = ideal_error + MGRIT_object.hierarchy(level-1).error{i};
        end
        
        
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
            MGRIT_object.D_space_T = get_toeplitz_spatial_disc(spatial_problem)'; % Take the transpose here 

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
%         B = MGRIT_object.I ...
%             + MGRIT_object.D1*((-ideal_error_FV + error_SL).*(MGRIT_object.D_space') ...
%             - ideal_error_RK.*(MGRIT_object.D_time'));        
        % End of assembling error correction matrix B.
        % End of construction of backward Euler matrix.
        
        B = MGRIT_object.I ...
            + MGRIT_object.D1*((-ideal_error + error_SL).*(MGRIT_object.D_space_T) );        
        

    %%
    elseif correction.weighted_stencil
        
        % Sum fine-grid truncation error coefficients across the given 
        % coarse interval. Note: We separate the FV and RK
        % contrubutions because they potentially are associated with
        % different difference matrices. 
%         ideal_error_RK = zeros(nx, 2);
%         ideal_error_FV = zeros(nx, 2);
%         
%         %ideal_error_FV_wave = zeros(nx, 1);
% 
%         % To make the code more intelligible, the error arrays on
%         % different levels are labelled slightly differently.
%         % On first coarse level, accumulate the finest-grid errors.
%         if level == 2
%             for i = t0_idx_fine:t1_idx_fine-1
%                 ideal_error_FV(:, 1) = ideal_error_FV(:, 1) ...
%                     + MGRIT_object.hierarchy(level-1).error_FV{i}{1};
%                 
%                 ideal_error_FV(:, 2) = ideal_error_FV(:, 2) ...
%                     + MGRIT_object.hierarchy(level-1).error_FV{i}{2};
%                 
% %                 ideal_error_FV_wave = ideal_error_FV_wave ...
% %                     + MGRIT_object.hierarchy(level-1).error_FV_wave{i};
% 
%                 ideal_error_RK(:, 1) = ideal_error_RK(:, 1) ...
%                     + MGRIT_object.hierarchy(level-1).error_RK{i}{1};
%                 
%                 ideal_error_RK(:, 2) = ideal_error_RK(:, 2) ...
%                     + MGRIT_object.hierarchy(level-1).error_RK{i}{2};
%             end
% 
%         % On all coarser levels, accumulate the errors on the level
%         % above us (which themselves are accumulations of errors above
%         % them).
%         else
%             error('multilevel is not implemented...')
%         end
        

        ideal_error_0 = zeros(nx, 1);
        ideal_error_1 = zeros(nx, 1);
        for i = t0_idx_fine:t1_idx_fine-1
            ideal_error_0 = ideal_error_0 ...
                + MGRIT_object.hierarchy(level-1).error0{i};

            ideal_error_1 = ideal_error_1 ...
                + MGRIT_object.hierarchy(level-1).error1{i};
         end


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

            % The left shifts here are relative to the stencil [u_{i-1}, u_{i}, u_{i+1}, u_{i+2}]
            % A left shift of zero is biased to the right. A left
            % shift of one is biased to the left.
            % D_2_0 <-- left shift of zero.
            p_nodes   = [0; 1; 2];
            p_weights = [1; -2; 1] / (h^2);
            spatial_problem.stencil  = {{p_weights; p_nodes}};             
            MGRIT_object.D_2_0 = get_toeplitz_spatial_disc(spatial_problem); 
            
            % D_2_1 <-- left shift of one.
            p_nodes   = [-1; 0; 1];
            p_weights = [1; -2; 1] / (h^2);
            spatial_problem.stencil  = {{p_weights; p_nodes}};             
            MGRIT_object.D_2_1 = get_toeplitz_spatial_disc(spatial_problem); 

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
        + MGRIT_object.D1*( ...
        + ideal_error_0.*MGRIT_object.D_2_0 ... % ell = 0
        + ideal_error_1.*MGRIT_object.D_2_1 ... % ell = 1
        + error_SL.*(MGRIT_object.D_space') ); 


    end
%     rng(1)
%     b = rand(nx, 1);
%     z = B \ b;
%     z(1:5)
%     Bw = B;
    
    
    
    
    % End of construction of backward Euler matrix.
    
    
    % Set handle for computing the action of Phi: MGRIT_object.hierarchy(level).Phi is a cell array
    % Use an LU factorization to solve BE system
    if strcmp(MGRIT_object.BE_coarse_solve.id, 'LU')
        %[B_lower, B_upper] = lu(B); % Factor B

        %MGRIT_object.hierarchy(level).SL_opers{t0_idx} = @(u0) B_upper\(B_lower\(my_SL_step(u0)));

        %u1 = B_upper\(B_lower\(my_SL_step(u0)));
        
        u1 = B \ (my_SL_step(u0));
        
    % Iteratively solve BE system with GMRES
    elseif strcmp(MGRIT_object.BE_coarse_solve.id, 'GMRES')

        restart = []; % This means no restarting, i.e., full-memory GMRES
%         MGRIT_object.hierarchy(level).SL_opers{t0_idx} = @(u0) ...
%             krylov_withoutoutput(@gmres, {B, my_SL_step(u0), ...
%                 restart, MGRIT_object.BE_coarse_solve.rtol, MGRIT_object.BE_coarse_solve.maxiter});

%         figure(1)
%         plot(abs(fft(my_SL_step(u0))))

        u1 = krylov_withoutoutput(@gmres, {B, my_SL_step(u0), ...
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