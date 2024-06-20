% Preconditioned residual iteration iteration to solve the system of 
% nonlinear algebraic equations arising from the space-time discretizations
% of a nonlinear conservation law. Examples include SWE and Euler.
%
% The basics of this script come from cons_law_time_stepping.m which uses
% sequential time-stepping to solve these PDEs. 
% 
% Select the PDE+domain+initial condition by ensuring the associated code
% is uncommented below.
%
% There are various settings that can be changed both for the outer
% nonlinear iteration and for the inner iteration.
%
% "linearized_solver": this option determines how the linearized systems
% are solved at each nonlinear iteration. Best-case convergence is achieved
% with "direct" which solves the systems with sequential time-stepping. To
% use characteristic-variable block-preconditioning choose "block-prec"
% A number of options are then specified inside the "prec" object:
%   "blocks":
%       "exact" means the preconditioner uses the true blocks from the
%           linearized time-stepping operator in characteristic variables
%       "approx" means the preconditioner is based on blocks that
%           approximate those in the above-mentioned time-stepping operator.
%   "type": Specified either diagonal or triangular preconditioning.
%   "diag_solves":
%       "direct" means diagonal blocks are inverted directly with
%           time-stepping
%       "MGRIT" means diagonal blocks are inverted inexactly using MGRIT
%   "maxit": Number of preconditioner applications.
%
% NOTE: MGRIT implementation only available for problems with periodic
% boundaries in space. 


tmp = matlab.desktop.editor.getActive;
cd(fileparts(tmp.Filename));

clc
clear
close all

%% How are the linearized problems solved at each nonlinear iteration?

linearized_solver = 'direct';     % Directly with time-stepping. Essentially a Newton method. 
linearized_solver = 'block-prec'; % Approximately using a block preconditioner

% Further options for block preconditioner. 
prec = struct();
if strcmp(linearized_solver, 'block-prec')
    prec.blocks = 'exact';
    prec.blocks = 'approx';
    
    prec.type = 'lower'; % Block lower triangular
    prec.type = 'diag';  % Block diagonal
    
    % If using approximate diagonal blocks in the linearized problem, how should they be inverted?
    prec.diag_solves = 'direct';
    prec.diag_solves = 'MGRIT';
end


%% PDE and problem parameters

nx_array = 2.^(6:11);

%disc_pa.num_flux_id = 'LLF';
disc_pa.num_flux_id = 'ROE'; disc_pa.delta_smoothing_parameter = 1e-6;


%% PDE, domain and initial condition parameters
%% SWE: Initial depth perturbation with Gaussian of size epsilon
pde_pa.pde_id = 'shallow-water';
pde_pa.ic_id = 'idp1'; 
pde_pa.bcs = 'periodic';

%pde_pa.ic_epsilon = 0.05; 

pde_pa.ic_epsilon = 0.1; % Only weak shocks form here.
prec.maxit = 1;

% pde_pa.ic_epsilon = 0.6; % Shocks form.
% prec.maxit = 1;
% prec.maxit = 2;
% % prec.maxit = 3;

mesh_pa.xmin = -5;
mesh_pa.xmax =  5;
disc_pa.CFL_number = 0.8; % max-wave-speed * dt/h
mesh_pa.tmax = 10;


%% SWE: Initial cosine pert to h
% pde_pa.pde_id = 'shallow-water';
% pde_pa.ic_id = 'idp2'; 
% 
% mesh_pa.xmin = -5;
% mesh_pa.xmax =  5;
% disc_pa.CFL_number = 0.8; % max-wave-speed * dt/h
% mesh_pa.tmax = 10;
% 
% prec.maxit = 1;
% pde_pa.ic_epsilon = 0.2; % Shocks don't form, but there is steepening
% 
% prec.maxit = 3;
% pde_pa.ic_epsilon = 0.6; % 3 iters with best prec is almost scalable.


%% SWE: Dam break problem: h0 has a jump of height eps.
% eps = 2 is the usual dam break problem from LeVeque p. 259
% pde_pa.pde_id = 'shallow-water';
% pde_pa.ic_id  = 'dam-break'; 
% pde_pa.bcs    = 'constant';
% 
% pde_pa.ic_epsilon = 0.1;
% prec.maxit = 1;
% 
% % pde_pa.ic_epsilon = 2;
% % prec.maxit = 1;
% % %prec.maxit = 2;
% % prec.maxit = 3;
% % prec.maxit = 4;
% 
% 
% mesh_pa.xmin = -10;
% mesh_pa.xmax =  10;
% disc_pa.CFL_number = 0.7; % max-wave-speed * dt/h
% mesh_pa.tmax = 5;


%% Euler: Smooth initial perturbation
% pde_pa.pde_id = 'euler'; 
% pde_pa.ic_id = 'idp1'; 
% 
% pde_pa.ic_epsilon = 0.2;
% prec.maxit = 1;
% 
% pde_pa.ic_epsilon = 1.2;
% prec.maxit = 1;
% % prec.maxit = 2;
% % prec.maxit = 3;
% 
% mesh_pa.xmin = -5;
% mesh_pa.xmax =  5;
% disc_pa.CFL_number = 0.7; % max-wave-speed * dt/h
% mesh_pa.tmax = 10;
% pde_pa.bcs = 'periodic';


%% Euler: Sod problem
% pde_pa.pde_id = 'euler'; 
% pde_pa.ic_id = 'sod'; 
% mesh_pa.xmin = 0;
% mesh_pa.xmax = 1;
% disc_pa.CFL_number = 0.45; % max-wave-speed * dt/h
% mesh_pa.tmax = 0.25;
% pde_pa.bcs = 'constant';
% 
% 
% pde_pa.ic_epsilon = 0.125; % Weakly nonlinear shock tube problem
% prec.maxit = 1;
% 
% 
% pde_pa.ic_epsilon = 0.875; % This is the original Sod problem
% prec.maxit = 2;
% % prec.maxit = 3;
% prec.maxit = 4;
% prec.maxit = 5;


%% Nonlinear iteration and linearization parameters

rtol_outer  = 1e-10;
maxit_outer = 15;
diverge_tol_outer = 1e3;

% Add nonlinear relaxation.
outer_solve_pa.cf = 8; % Coarsening factor used to define FC splitting. 
outer_solve_pa.relax_scheme = 'F'; 
%outer_solve_pa.relax_scheme = 'FCF'; 

resnorm = 2;
errnorm = 2;

% Linearize about the exact solution and run the block-preconditioned
% iteration on it.
run_inner_tests = ~true;
plot_inner_residual_space_time = ~true;


%% MGRIT parameters for space-time linear advection solves
if strcmp(linearized_solver, 'block-prec') && strcmp(prec.diag_solves, 'MGRIT')

    MGRIT_maxiter         = 1; % Maximum number of MGRIT iterations per richardson iteration.
    MGRIT_cf              = outer_solve_pa.cf; % MGRIT coarsening factor. Must be the same as for the outer nonlinear iteration!
    MGRIT_maxlevels       = 12; % MGRIT max levels.
    
    MGRIT_relax  = 'F-FCF'; % F-relax on level 1, and FCF on all other levels.
    
    MGRIT_res_halt_tol = 0; % MGRIT relative residual halting tol. We want to apply a single MGRIT iteration. We don't care what the residual is.
    MGRIT_verbose      = false;
    
    % Choose how coarse BE systems are solved
    BE_coarse_solver = 'GMRES'; % Further parameters are set below.
    %BE_coarse_solver = 'LU';
    %BE_coarse_solver = 'NONE';
    
    % Should the F-relax be applied at the end of an MGRIT V-cycle? If the
    % outer iteration begins with an F-relaxation, then no. Recall that the
    % outer iteration does u <- u + e, with e the MGRIT solution. Then
    % if nonlinear relaxation is done that begins with F, all F-points of u are
    % immediately overwritten. Thus, there's no point in ensuring that e is up
    % to date at F-points. 
    MGRIT_final_F_relax = ~true;

    % Package MGRIT parameters
    myMGRIT_solver_params = struct();
    myMGRIT_solver_params.cf            = MGRIT_cf;
    myMGRIT_solver_params.maxlevels     = MGRIT_maxlevels;
    myMGRIT_solver_params.maxiter       = MGRIT_maxiter;
    myMGRIT_solver_params.final_F_relax = MGRIT_final_F_relax;
    myMGRIT_solver_params.res_halt_tol  = MGRIT_res_halt_tol;
    myMGRIT_solver_params.min_coarse_nt = 2;
    myMGRIT_solver_params.verbose       = MGRIT_verbose;

    if strcmp(MGRIT_relax, 'F')
        myMGRIT_solver_params.pre_relax = 'F';
    elseif strcmp(MGRIT_relax, 'FCF')
        myMGRIT_solver_params.pre_relax ='FCF';
    elseif strcmp(MGRIT_relax, 'F-FCF')
        myMGRIT_solver_params.pre_relax = @(level) F_then_FCF_relax(level);
    end

end
% End of setting MGRIT options


%% Plotting options
% If res fig is saved then the other figures produced are also saved. 
% Further plotting settings are also determined later when saving the figure.
save_res_fig = ~true; 
res_fig_dir = './figures/paper/';

plot_pa.ncon_lvl = 11;
plot_pa.con_tol  = 1e-7;
plot_pa.invisible_col_bar = ~true;

plot_pa.plot_solution_contours         = true;
plot_pa.plot_solution_cross_sections   = true;
plot_pa.plot_wave_speed_contours       = ~true;
plot_pa.plot_wave_speed_cross_sections = ~true;

figno = 20;
show_error_history = ~true;

% Line style 
if strcmp(linearized_solver, 'direct')
    plot_pa.ls = ':'; % Direct solves uses dotted lines

elseif strcmp(linearized_solver, 'block-prec')
    if strcmp(prec.type, 'diag') % P_D uses solid lines
        plot_pa.ls = '-';
    elseif strcmp(prec.type, 'lower') % P_L uses dashed lines
        plot_pa.ls = '--';
    end

    if strcmp(prec.diag_solves, 'MGRIT') 
        plot_pa.ls = '-.';
    end
end



%% Loop over different mesh resolutions
for nx_idx = 1:numel(nx_array)
    
    % Re-name some quantities needed for nested iteration
    if nx_idx > 1
        mesh_pa_coarse = mesh_pa;
        q_coarse = q;
    end
    
    %% Setup problem at current resolution
    nx = nx_array(nx_idx);
    mesh_pa.nx = nx;
    
    %% Set up mesh
    % Create spatial mesh
    [mesh_pa.x_interfaces, mesh_pa.x_centers, mesh_pa.h] = space_mesh(mesh_pa.xmin, mesh_pa.xmax, mesh_pa.nx);

    % Get initial conditon.
    if strcmp(pde_pa.pde_id, 'shallow-water')
        my_cons_law = swe_system(pde_pa);
        [h0, u0, ic_description] = my_cons_law.initial_condition(mesh_pa.x_centers);
        q0 = [h0; h0.*u0];

    elseif strcmp(pde_pa.pde_id, 'euler')
        my_cons_law = euler_system(pde_pa);
        [rho0, u0, E0, ic_description] = my_cons_law.initial_condition(mesh_pa.x_centers);
        q0 = [rho0; rho0.*u0; E0];

    end

    % Compute fastest wave-speed associated with initial condition
    my_cons_law.mesh_pa = mesh_pa; % The wave-speed function requires this to be set
    lambda0 = my_cons_law.wave_speeds(q0);
    abs_fprime0_max = max(abs(lambda0(:)));

    % Get temporal mesh.
    [mesh_pa.t, mesh_pa.dt] = time_mesh(mesh_pa.h, disc_pa.CFL_number, abs_fprime0_max, 1, mesh_pa.tmax);
    mesh_pa.nt = numel(mesh_pa.t);

    fprintf('\nnx=%d, nt=%d\n', mesh_pa.nx, mesh_pa.nt);
    
    %% Finalize construction of PDE class
    my_cons_law.disc_pa = disc_pa;
    my_cons_law.mesh_pa = mesh_pa;
    
    
    %% Build all-at-once nonlinear and linearized systems
    % Create data structure for holding linearization data. Note this is
    % updated as if it's a pointer.
    if strcmp(pde_pa.pde_id, 'shallow-water')
        linearization_data = swe_linearization_data(mesh_pa.nt);
    elseif strcmp(pde_pa.pde_id, 'euler')
        linearization_data = euler_linearization_data(mesh_pa.nt);
    end
    % Store in the conservation law so it can access it as required.
    my_cons_law.linearization_data = linearization_data;
    
    disc_pa.step          = @(t0idx, q) my_cons_law.step(t0idx, q); % Step in the nonlinear problem.
    my_cons_law_st_system = one_step_st_system(disc_pa.step, my_cons_law.m * nx, mesh_pa.nt, outer_solve_pa);     

    disc_pa_linearized = disc_pa;
    disc_pa_linearized.step          = @(t0idx, q) my_cons_law.step_linearized(t0idx, q); % Step in the linearized problem.
    my_linearized_cons_law_st_system = one_step_st_system(disc_pa_linearized.step, my_cons_law.m * nx, mesh_pa.nt, outer_solve_pa);     
        
    % Build all-at-once RHS vector
    b = zeros(my_cons_law.m * mesh_pa.nx, mesh_pa.nt);
    b(:, 1) = q0;
    b = b(:);
    
    % Get exact solution on current mesh by time-stepping
    q_time_stepping = my_cons_law_st_system.forward_solve(b);

    
    %% Nested iteration
    % Get exact solution with time-stepping and then skip to next
    % resolution.
    if nx_idx == 1
    
        q = q_time_stepping;
        continue
    
    %% Interpolate the coarse-grid solution
    elseif nx_idx > 1

        % Get space-time interpolation matrix.
        nx_coarse = mesh_pa_coarse.nx;
        nt_coarse = mesh_pa_coarse.nt;
        q_coarse = reshape(q_coarse, [my_cons_law.m * nx_coarse, nt_coarse]);
        t_coarse = mesh_pa_coarse.t;
        t        = mesh_pa.t;
        nt       = mesh_pa.nt;
        I        = bilinear_space_time_interpolation_FV(nx, t_coarse, t, strcmp(pde_pa.bcs, 'periodic'));

        % Unpack all components of q, interpolating them one at a time.
        q_init = zeros(my_cons_law.m * nx, mesh_pa.nt);
        for k = 1:my_cons_law.m
            qk_coarse = q_coarse( (k-1)*nx_coarse+1:k*nx_coarse, :);

            qk_init   = I*qk_coarse(:); 
            qk_init   = reshape(qk_init, [nx, nt]);

            % qk_init = qk_init + 0.01*randn(size(qk_init)); Adds random
            % noise to initial iterate... Useful sometimes.

            % Insert the correct fine-grid initial condition
            qk_init(:, 1) = q0( (k-1)*nx+1:k*nx );

            q_init( (k-1)*nx+1:k*nx, :) = qk_init;
        end
        q = q_init;
    end
    
    % Start from the exact solution, linearizing about it
    if run_inner_tests 
        q = q_time_stepping;
    end
    
    %% Compute initial residual. 
    [q, r] = my_cons_law_st_system.relaxation(q, b); 
    resnorm_array = norm(r, resnorm);
    errnorm_array = norm(q(:) - q_time_stepping(:), resnorm);
    % Print residual norm.
    fprintf('it %d: ||r_0|| = %.2e\n', 0, resnorm_array(1))

    
    %% Setup MGRIT object
    if strcmp(linearized_solver, 'block-prec') && strcmp(prec.diag_solves, 'MGRIT')
        myMGRIT_object = struct();
        myMGRIT_object.BE_coarse_solve.id      = BE_coarse_solver;
        myMGRIT_object.BE_coarse_solve.maxiter = 10;
        myMGRIT_object.BE_coarse_solve.rtol    = 0.01;
    
        myMGRIT_object.t = mesh_pa.t;
        myMGRIT_object.block_size = mesh_pa.nx;
    
        % Store MGRIT object in prec object
        prec.myMGRIT_object        = myMGRIT_object;
        prec.myMGRIT_solver_params = myMGRIT_solver_params;
    end

    %% Iterate Richardson 
    for iter_outer = 1:maxit_outer
        
        % Check if halting tolerance reached
        if resnorm_array(end)/resnorm_array(1) < rtol_outer
            break
        % Check if diverged
        elseif resnorm_array(end)/resnorm_array(1) > diverge_tol_outer 
            error('rich iteration has diverged...')
            break
        end

        % Solve linearized problem A*e = r directly.
        if strcmp(linearized_solver, 'direct')
            e = my_linearized_cons_law_st_system.forward_solve(r);
            
        % Approximately solve linearized problem A*e = r using a block preconditioner
        elseif strcmp(linearized_solver, 'block-prec')
            % Reshape to facillitate eigenvector transformations...
            q = reshape(q, [my_cons_law.m * mesh_pa.nx, mesh_pa.nt]);
            
            % Create handles for mapping back and forth between variables; these
            % are parsed to the space-time system class
            R_handle    = @(t0idx, b) my_cons_law.right_eigenvector_map(q(:, t0idx), b);
            Rinv_handle = @(t0idx, a) my_cons_law.left_eigenvector_map( q(:, t0idx), a);
            
            
            if ~run_inner_tests
                %e = zeros(size(r)); % Guess e = 0
                e = r; % Guess e = r
            else
                e = randn(size(r));
            end
            
            r_lin_norm_array = [];
            
            % One iteration of the preconditioner has an F-relax+C-point 
            % residual evaluation, followed by a solve of a block linear
            % system.
            for iter = 1:prec.maxit
                % Initial relaxation
                [e, r_lin] = my_linearized_cons_law_st_system.relaxation(e, r);
                r_lin_norm_array = [r_lin_norm_array; norm(r_lin, resnorm)];
    
                % Map linearized residual from primitive space to characteristic
                r_lin_hat = my_linearized_cons_law_st_system.all_point_transformation(r_lin, Rinv_handle);

                % Preconditioner uses exact blocks
                if strcmp(prec.blocks, 'exact')
                    e_lin_hat = block_prec_exact_blocks(mesh_pa, my_cons_law, r_lin_hat, Rinv_handle, R_handle, prec);

                % Preconditioner uses approximate blocks
                elseif strcmp(prec.blocks, 'approx')
                    e_lin_hat = block_prec_approx_blocks(mesh_pa, my_cons_law, r_lin_hat, Rinv_handle, R_handle, prec);

                end

                % Map linearized error correction from characteristic space to primitive space
                e_lin = my_linearized_cons_law_st_system.all_point_transformation(e_lin_hat, R_handle);

                % Add correction
                e = e + e_lin;
                
                % Plot the current algebraic error in space-time.
                if plot_inner_residual_space_time
                    fh = plot_space_time_error(figno, r_lin, pde_pa, mesh_pa, iter); figno = figno+1;
                    pause
                end
            end
            
            q = q(:);
            
            if run_inner_tests; break; end
        end

        % Add linearized correction to update q.
        q = q + e;

        % Nonlinear relaxation
        [q, r] = my_cons_law_st_system.relaxation(q, b); 
        
        % Compute residual and error norms
        resnorm_array = [resnorm_array; norm(r, resnorm)];
        e_alg = q(:) - q_time_stepping(:);
        errnorm_array = [errnorm_array; norm(e_alg, errnorm)];

%         fh = plot_space_time_error(figno+2+iter_outer, e_alg, pde_pa, mesh_pa, iter_outer);
%         pause
        
        % Print residual norm.
        fprintf('it %d: ||r_%d||/||r_0|| = %.2e, ||r_%d||/||r_%d|| = %.2f\n', ...
            iter_outer, iter_outer, resnorm_array(end)/resnorm_array(1), ...
            iter_outer, iter_outer-1, resnorm_array(end)/resnorm_array(end-1))
    end
    % End of Rich iter loop
    
    nameval = lineprops(nx_idx);
    nameval{4} = plot_pa.ls;
    
    %% Plots: Richardson convergence
    if run_inner_tests
        resnorm_array = r_lin_norm_array;
    end
    
    % Residuals
    res_fig = figure(figno);
    if resnorm == 1
        resnorm_array = resnorm_array * (mesh_pa.h*mesh_pa.dt);
    elseif resnorm == 2
        resnorm_array = resnorm_array * sqrt(mesh_pa.h*mesh_pa.dt);
    end
    displayname = sprintf('$(n_x,n_t)=(%d,%d)$', mesh_pa.nx, mesh_pa.nt);
    semilogy(0:numel(resnorm_array)-1, resnorm_array/resnorm_array(1), nameval{:}, ...
        'LineWidth', 2, ...
        'DisplayName', displayname)
    hold on
    
    % Error
    if show_error_history
        figure(figno+1)
        if errnorm == 1
            errnorm_array = errnorm_array * (mesh_pa.h*mesh_pa.dt);
        elseif errnorm == 2
            errnorm_array = errnorm_array * sqrt(mesh_pa.h*mesh_pa.dt);
        end
        semilogy(0:numel(errnorm_array)-1, errnorm_array, nameval{:}, ...
            'LineWidth', 2, ...
            'DisplayName', displayname)
        hold on
    end

    drawnow
    
    
end
% End of looping over different mesh resolutions.

%% Pretty Richardson conv. plots
figure(figno)

axis tight
ylim([1e-10, 10])
xlim([0 maxit_outer])
yticks((10).^(-10:2:0))
grid minor
grid on
grid minor
set(gca,'YMinorTick', 'Off')
%grid minor

lh = legend();
lh.set('Location', 'best')
xlabel('$k$')
ylabel('$\Vert \mathbf{r}_k \Vert / \Vert \mathbf{r}_0 \Vert$')

% Title includes linearized solver description
if strcmp(linearized_solver, 'direct')
    title('residual history: direct solve')
else
    if strcmp(linearized_solver, 'block-prec')
        % block-prec-exact-blocks  -> \wh{\cal P}
        if strcmp(prec.blocks, 'exact')
            res_fig_title = strcat('residual history: $\widehat{\mathcal{P}}$', sprintf('(%d)', prec.maxit));
        % block-prec-approx-blocks -> \wt{\cal P}
        elseif strcmp(prec.blocks, 'approx')
            res_fig_title = strcat('residual history: $\widetilde{\mathcal{P}}$', sprintf('(%d)', prec.maxit));
        end
    end
    title(res_fig_title);
end



% Save plot of residual history
if save_res_fig
    
    % Create string describing problem
    res_fig_name = 'res-';

    % PDE description
    if strcmp(pde_pa.pde_id, 'shallow-water')
        res_fig_name = strcat(res_fig_name, 'SWE');
        res_fig_dir = strcat(res_fig_dir, '/SWE/', pde_pa.ic_id);

    elseif strcmp(pde_pa.pde_id, 'euler')
        res_fig_name = strcat(res_fig_name, 'euler');
        res_fig_dir = strcat(res_fig_dir, '/euler/', pde_pa.ic_id);
    end

    % Initial condition description
    res_fig_name = strcat(res_fig_name, '-', ic_description);

    % Linearized solver description
    if strcmp(linearized_solver, 'direct')
        res_fig_name = strcat(res_fig_name, '-direct');

    % Inexact solvers have more parameters
    elseif strcmp(linearized_solver, 'block-prec')
        % block-prec-exact-blocks  -> bpeb
        if strcmp(prec.blocks, 'exact')
            res_fig_name = strcat(res_fig_name, '-bpeb');

        % block-prec-approx-blocks -> bpab
        elseif strcmp(prec.blocks, 'approx')
            res_fig_name = strcat(res_fig_name, '-bpab');
        end

        % Add MGRIT into name if using it. 
        if strcmp(prec.diag_solves, 'MGRIT')
            res_fig_name = strcat(res_fig_name, '-MGRIT');
        end

        res_fig_name = strcat(res_fig_name, '-', prec.type, '-it', num2str(prec.maxit));
    end

    % Save the figure
    figure(res_fig);
    figure_saver(gcf, sprintf('%s/%s', res_fig_dir, res_fig_name), true); 
end
figno = figno+1;



if show_error_history && ~run_inner_tests
    figure(figno)
    xlabel('$k$')
    title('error history')
    ylabel('$\Vert \mathbf{e}_k \Vert$')
    xlim([0 maxit_outer])
    ylim([rtol_outer, 10])
    axis tight
    box on
    lh = legend();
    lh.set('Location', 'Best')
    figno = figno+1;
end




%% Plot solution, etc.
figno = figno+1;

if save_res_fig
    plot_pa.save_figs = true;
    plot_pa.save_dir = res_fig_dir;
    plot_pa.save_str = ic_description;
end

if contains(ic_description, 'sod')
    plot_pa.title_str = sprintf('$\\varepsilon = %.3f$', pde_pa.ic_epsilon); % Sod problem needs 3 digits visible
else
    plot_pa.title_str = sprintf('$\\varepsilon = %.1f$', pde_pa.ic_epsilon);
end
[fh, figno] = plot_cons_law_system_solutions(q, my_cons_law, pde_pa, disc_pa, mesh_pa, figno, plot_pa);



% Plot error and residual in space time
function fh = plot_space_time_error(figno, r, pde_pa, mesh_pa, iteration_idx)
    fh = figure(figno);
    nx = mesh_pa.nx;
    nt = mesh_pa.nt;
    [X, T] = meshgrid(mesh_pa.x_centers, mesh_pa.t);
    r = reshape(r, [pde_pa.m*nx, nt]);
    %subplot(1, 2, 1)
    r_logged = log10(abs(r(1:nx, :)))';
    r_logged(r_logged < -8) = NaN;
    mesh(X, T, r_logged)
    %subplot(1, 2, 2)
    %r_reshaped = reshape(r, [2*nx nt]);
    %mesh(X, T, log10(abs(r_reshaped(1:nx, :)))')
    title(sprintf('$n_x = %d$. it=%d. error $= p - p_k$', nx, iteration_idx))
    xlim([mesh_pa.x_centers(1), mesh_pa.x_centers(end)])
    ylim([mesh_pa.t(1), mesh_pa.t(end)])
    zlim([-8, 0])
    xlabel('$x$')
    ylabel('$t$')
    zlabel('$\log_{10}| r_k |$')
end


% Solve P*e = r where P is a block m x m preconditioner.
% Diagonal blocks of P are approximate and are inverted exactly.
function [e, inner_solve_pa] = block_prec_approx_blocks(mesh_pa, my_cons_law, r, L_handle, R_handle, prec)

    %assert(~strcmp(prec_id, 'lower'), 'lower prec not implemented for approx blocks. Only diag')

    nx = mesh_pa.nx;
    nt = mesh_pa.nt;

    r  = reshape(r, [my_cons_law.m * nx, nt]);
    e  = zeros(size(r)); 
    
    % Solve for kth char variable
    for char_idx = 1:my_cons_law.m

        char_inds = (char_idx-1)*nx+1:char_idx*nx; % Spatial indices of kth char var

        rk = r(char_inds, :);      
        ek = zeros(size(rk));
       

        %% Direct inversion of diagonal block via sequential time-stepping 
        if strcmp(prec.diag_solves, 'direct')
   
            ek(:, 1) = rk(:, 1); % Initial condition
            for n = 1:mesh_pa.nt-1
    
                ek(:, n+1) = my_cons_law.step_linearized_characteristic_variable(n, ek(:, n), char_idx) + rk(:, n+1);
    
                % Add in coupling to first k-1 char variables at time n if using lower triangular prec.
                if strcmp(prec.type, 'lower')
                    if char_idx > 1
                        ek_pad = zeros(my_cons_law.m * nx, 1); % 
                        ek_pad(1:char_inds(1)-1) = e(1:char_inds(1)-1, n); 
                    
                        % Compute a vector temp that is Phi^{hat, n} * ek_padded
                        % temp <- R*e
                        temp = R_handle(n, ek_pad);
                        % temp <- Phi*temp
                        temp = my_cons_law.step_linearized(n, temp);
                        % temp <- L*temp. Note L is the left eigenvectors at time n+1, not n.
                        temp = L_handle(n+1, temp);
                        ek(:, n+1) = ek(:, n+1) + temp(char_inds);
                    end
                end
            end

        %% Approximate inversion of diagonal block using MGRIT
        elseif strcmp(prec.diag_solves, 'MGRIT')

            % Triangular prec would require us to update the global RHS
            % vector that's parsed to MGRIT. 
            assert(strcmp(prec.type, 'diag'), 'Inexact inversion of diagonal blocks only implemented for block diag prec.')

            %step_cons_law_lin_char_MGRIT(u0, step_status, MGRIT_object, my_cons_law_system, char_idx)

            % Faltten space-time matrix into vectors
            rk = rk(:);
            ek = rk; % Initial iterate is the right-hand side vector
            [ek, ~, prec.myMGRIT_object] = ...
                mgrit(ek, rk, @(a, b, c) step_cons_law_lin_char_MGRIT(a, b, c, my_cons_law, char_idx), ...
                prec.myMGRIT_object, prec.myMGRIT_solver_params);

            % Reshape into space-time matrix.
            ek  = reshape(ek, [nx, nt]);
        end

        % Insert ek into the global e matrix.
        e(char_inds, :) = ek;
    end
    e = e(:);
        
end



% Solve P*e^{hat} = r^{hat} where P is a block m x mx preconditioner and 
% e_hat and r_hat are error and residual in characteristic space.
% The true system reads e^{hat, n+1} = Phi^{hat, n} * e^{hat, n} + r^{hat, n+1}
% Phi^{hat,n} = L^{n+1} * Phi^{n} * R^{n} is an m-dimensional block matrix 
% The R maps the char var e^{hat} into a cons var, then it's hit with Phi, 
% and then it's mapped back into a char var by L.
%
% Phi^{hat} is an m x m block matrix that is responsible for coupling the
% characteristic variables together. It's not implemented directly here,
% but it's implemented using the actions of L, Phi, and R.
% 
% For the block diagonal preconditioner, Phi^{hat} is approximated by its 
% block diagonal 
% For the block lower triangular preconditioner, Phi^{hat} is approximated 
% by its block lower triangular part.
%
% Say we want to just compute the action of the block diagonal of Phi on a 
% block vector e, then (Phi*ek)_k is the action of the kth diagonal block
% on the kth block of e if ek is a m-dimensional block vector with zero
% blocks except its kth block equal to that of e.
% Similarly, to compute the kth block of the action of the block lower 
% triangular part of Phi on e is (Phi*ek)_k where now ek is the same as e,
% except its blocks j > k are zero'd out.
%
% So this is a really inefficient way to implement the preconditioner
% because it requires k applications of Phi^{hat} to compute the action of
% its block diagonal or block lower triangular. Anyway, it's the only way I
% can think to do it.
%
function e = block_prec_exact_blocks(mesh_pa, my_cons_law, r, L_handle, R_handle, prec)

    assert(strcmp(prec.diag_solves, 'direct'), 'Only direct solves implemented for exact-block prec')

    nx = mesh_pa.nx;
    nt = mesh_pa.nt;
    r  = reshape(r, [my_cons_law.m * nx, nt]);
    e  = zeros(size(r)); 
    
    % Global temporal solve for kth characteristic variable. For the block diagonal 
    % preconditioner all of these solves are independent. For the block 
    % lower triangular preconditioner, the kth characteristic variable
    % depends on all the variables that were solved before it.
    for char_idx = 1:my_cons_law.m
        char_inds = (char_idx-1)*nx+1:char_idx*nx; % Spatial indices of kth char var

        rk = r(char_inds, :);      
        ek = zeros(size(rk));

        % Create m*nx - dimensional vector. This is a place holder for the
        % block components from e(:, n) required to compute the kth block
        % of either the block diagonal of Phi^{hat,n} applied to en or the
        % block lower triangular of Phi^{hat,n} applied to en.
        ek_pad = zeros(my_cons_law.m * nx, 1); % 

        ek(:, 1) = rk(:, 1); % Initial condition

        %error_max = 0.0;

        for n = 1:mesh_pa.nt-1

            % Pack ek^n into bigger vector with zeros everywhere else.
            ek_pad(char_inds) = ek(:, n); 

            % Add in coupling to first k-1 char variables if using lower triangular prec.
            if strcmp(prec.type, 'lower')
                if char_idx > 1
                    ek_pad(1:char_inds(1)-1) = e(1:char_inds(1)-1, n); 
                end
            end

            % Compute a vector temp that is Phi^{hat, n} * ek_padded
            % temp <- R*e
            temp = R_handle(n, ek_pad);
            % temp <- Phi*temp
            temp = my_cons_law.step_linearized(n, temp);
            % temp <- L*temp. Note L is the left eigenvectors at time n+1, not n.
            temp = L_handle(n+1, temp);

            % This code snippet can be used to compare the action of the
            % true kth diagonal block vs. that resulting from the
            % approximate block implementation.
            % temp_true = temp(k_inds);
            %test_vec = ek(:, n);
            % temp_test = my_cons_law.step_linearized_characteristic_variable(n, test_vec, k);
            % 
            % %[temp_true, temp_test, abs(temp_true-temp_test) ./ norm(temp_true, inf)]
            % error_max_temp = norm(temp_true - temp_test, inf) ./ norm(temp_true, inf);
            % error_max = max(error_max, error_max_temp);


            % Compute ek^{n+1} by extracting the kth block from temp.
            ek(:, n+1) = temp(char_inds) + rk(:, n+1);
        end
        %error_max
        e(char_inds, :) = ek;
    end
    e = e(:);
        
end

% Implement level-dependent F-FCF MGRIT relaxation
function nu = F_then_FCF_relax(level)
    if level == 1 
        nu = 'F';
    else
        nu = 'FCF';
    end
end