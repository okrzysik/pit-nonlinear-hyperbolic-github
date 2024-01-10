% Applies a preconditioned residual correction iteration to solve the
% discretized space-time system A(u) = b.
%
% Linearized systems may be solved DIRECTLY with sequential time-stepping
% or may be solved APPROXIMATELY using MGRIT. 

tmp = matlab.desktop.editor.getActive;
cd(fileparts(tmp.Filename));

clear
clc
close all



%% Nonlinear iteration parameters
rtol_outer  = 1e-10;
maxit_outer = 20;
diverge_tol_outer = 1e3;

% Add nonlinear relaxation.
outer_solve_pa.cf = 8; % Coarsening factor used to define FC splitting. 
outer_solve_pa.relax_scheme = ''; % No relax.
outer_solve_pa.relax_scheme = 'F'; 
%outer_solve_pa.relax_scheme = 'FC';
%outer_solve_pa.relax_scheme = 'FCF';

% In which norm do we measure conv. of iterates?
resnorm = 2;
errnorm = 2;

initial_iterate = 'nested-iteration-it'; % Interpolate from iterate solution, unless solver diverged, then use ts solution
%initial_iterate = 'nested-iteration-ts'; % Interpolate from time-stepped solution
%initial_iterate = 'rand-pert-of-exact-sol'; rand_pert_mag = 1e-2;

%% linearization parameters
weno_linearization = 'picard'; % Freezes weights 
%weno_linearization = 'newton'; % Uses gradient 
weno_linearization = 'newton-FD-approx'; % Approximates gradient with FD


%% Linear solver parameters
inner_solve_pa.linear_solve = 'direct';
inner_solve_pa.linear_solve = 'MGRIT';

if strcmp(inner_solve_pa.linear_solve, 'MGRIT')

MGRIT_maxiter         = 1; % Maximum number of MGRIT iterations per richardson iteration.
MGRIT_cf              = outer_solve_pa.cf; % MGRIT coarsening factor. Must be the same as for the outer nonlinear iteration!
MGRIT_maxlevels       = 2; % MGRIT max levels.

%MGRIT_relax  = 'FCF';
MGRIT_relax  = 'F';
%MGRIT_relax  = 'CF';

MGRIT_res_halt_tol = 0; % MGRIT relative residual halting tol. We want to apply a single MGRIT iteration. We don't care what the residual is.
MGRIT_verbose      = false;

% Choose how coarse BE systems are solved
%BE_coarse_solver = 'GMRES'; % Further parameters are set below.
BE_coarse_solver = 'LU';
%BE_coarse_solver = 'NONE';

% Should the F-relax be applied at the end of an MGRIT V-cycle? If the
% outer iteration begins with an F-relaxation, then no. Recall that the
% outer iteration does u <- u + e, with e the MGRIT solution. Then
% if nonlinear relaxation is done that begins with F, all F-points of u are
% immediately overwritten. Thus, there's no point in ensuring that e is up
% to date at F-points. 
MGRIT_final_F_relax = ~true;

end


%% Discretization parameters
spatial_order = 1;
spatial_order = 3;

reconstruction_id = 'linear';
reconstruction_id = 'WENO';

num_flux_id = 'GLF'; % The most dissipative option, but this is differentiable.
num_flux_id = 'LLF'; % The usual local LF flux. NOT differentiable.


%nx_array = 2^7;
%if spatial_order == 1; nx_array = 2.^(5:9); end
if spatial_order == 1; nx_array = 2.^(5:9); end
if spatial_order == 3; nx_array = 2.^(5:7); end
%if spatial_order == 3; nx_array = 2.^(7); end

CFL_number = 0.8;

% Only limit reconstructions for LLF solution of Buckley--Leverett
limit_reconstructions = false;


%% PDE parameters

pde_id = 'burgers'; 
%u0_id = 1; tmax = 4; 
%u0_id = 2; tmax = 4; 
u0_id = 3; tmax = 4; 

% pde_id = 'buckley-leverett'; if strcmp(num_flux_id, 'LLF'); limit_reconstructions = true; end
% u0_id = 3; tmax = 2; % Riemann problem. Two compound waves.

xmin = -1;
xmax = 1;


%% Plotting
figno = 108;
save_res_fig    = ~true; % Save plot of residual history
save_matlab_fig = ~true; % Save the native MATLAB figure too.
plot_pa.fig_dir = './figures/paper/inexact/';


show_st_alg_error = ~true; % Make a plot of st error at each Rich iter
show_iter_error   = ~true; % Plot error norms.

% Line style 
if strcmp(inner_solve_pa.linear_solve, 'direct')
    if strcmp(outer_solve_pa.relax_scheme, '') 
        plot_pa.ls = '-'; % Direct solve without relax is solid line
    else
        plot_pa.ls = '--'; % Direct solve with relax is dashed line
    end
elseif strcmp(inner_solve_pa.linear_solve, 'MGRIT')
    plot_pa.ls = '-'; % MGRIT solve is solid line.
end  

plot_pa.PDE_ic  = false; % Show IC in title
plot_pa.relax   = false; % Show relax scheme in title

%plot_pa.lh_loc = 'SouthWest';
plot_pa.lh_loc = 'NorthEast';

plot_pa.ncon_lvl = 10;
plot_pa.con_tol  = 1e-3;
plot_pa.invisible_col_bar = ~true;

plot_pa.plot_solution_contours       = true;
plot_pa.plot_solution_cross_sections = true;


%% Create nonlinear PDE object.
pde_pa.ic_id = u0_id;
if strcmp(pde_id, 'burgers')
    my_cons_law = burgers(pde_pa);
    
elseif strcmp(pde_id, 'buckley-leverett')
    my_cons_law = buckley_leverett(pde_pa);
    
end

% Get max initial wave-speed
my_cons_law.compute_and_set_zero_extremal_values(xmin, xmax);


%% Package mesh-resolution-independent parameters
% Package mesh params
mesh_pa.xmin = xmin;
mesh_pa.xmax = xmax;
mesh_pa.tmax = tmax;

% Package disc params
disc_pa.spatial_order         = spatial_order;
disc_pa.limit_reconstructions = limit_reconstructions;
disc_pa.CFL_number            = CFL_number;
disc_pa.reconstruction_id     = reconstruction_id;
disc_pa.num_flux_id           = num_flux_id;

% Create linearization pa
linearization_pa = scalar_linearization_pa(weno_linearization);

%% Package MGRIT parameters
if strcmp(inner_solve_pa.linear_solve, 'MGRIT')
myMGRIT_solver_params = struct();
myMGRIT_solver_params.cf            = MGRIT_cf;
myMGRIT_solver_params.maxlevels     = MGRIT_maxlevels;
myMGRIT_solver_params.maxiter       = MGRIT_maxiter;
myMGRIT_solver_params.final_F_relax = MGRIT_final_F_relax;
myMGRIT_solver_params.pre_relax     = MGRIT_relax;
myMGRIT_solver_params.res_halt_tol  = MGRIT_res_halt_tol;
myMGRIT_solver_params.min_coarse_nt = 2;
myMGRIT_solver_params.verbose       = MGRIT_verbose;
end


%% Loop over different mesh resolutions
tstart_nx_loop = tic;
solve_diverged        = ~true; % Flag for current iterate having diverged
coarse_solve_diverged = ~true; % Flag for iterate on coarser mesh having diverged
for nx_idx = 1:numel(nx_array)
    
    
    % Re-name some quantities needed for nested iteration
    if nx_idx > 1
        mesh_pa_coarse = mesh_pa;
        u_coarse = u;
        if solve_diverged; coarse_solve_diverged = true; end
        solve_diverged = false; 
    end
    
    %% Compute and setup mesh-resolution-dependent parameters
    % Get spatial mesh.
    mesh_pa.nx = nx_array(nx_idx);
    [mesh_pa.x_interfaces, mesh_pa.x_centers, mesh_pa.h] = space_mesh(xmin, xmax, mesh_pa.nx);

    % Get temporal mesh.
    [mesh_pa.t, mesh_pa.dt] = time_mesh(mesh_pa.h, CFL_number, my_cons_law.f0_prime_max, disc_pa.spatial_order, mesh_pa.tmax);
    mesh_pa.nt = numel(mesh_pa.t);
    fprintf('nx=%d, nt=%d\n', mesh_pa.nx, mesh_pa.nt)

    % Initialize class for spatial reconstructions
    reconstruction = weighted_reconstruction(disc_pa.spatial_order, mesh_pa.nx, disc_pa.reconstruction_id);

    %% Finalize construction of cons_law_scalar
    my_cons_law.disc_pa          = disc_pa;
    my_cons_law.mesh_pa          = mesh_pa;
    my_cons_law.linearization_pa = linearization_pa;
    my_cons_law.reconstruction   = reconstruction;
    
    % Create data structure for holding linearization data. Note this is
    % updated as if it's a pointer.
    linearization_data = scalar_linearization_data(mesh_pa.nt);
    % Store in the conservation law so it can access it as required.
    my_cons_law.linearization_data = linearization_data;
    
    
    %% Build all-at-once systems for nonlinear and linearized problems
    step                  = @(t0idx, u) my_cons_law.step(t0idx, u); % Step in the nonlinear problem.
    my_cons_law_st_system = one_step_st_system(step, mesh_pa.nx, mesh_pa.nt, outer_solve_pa);     

    step_linearized    = @(t0idx, u) my_cons_law.step_linearized(t0idx, u); % Step in the linearized problem.
    my_linearized_cons_law_st_system = one_step_st_system(step_linearized, mesh_pa.nx, mesh_pa.nt, outer_solve_pa);     
    
    % Build all-at-once RHS vector for nonlinear problem
    u0 = cell_average(@(x) my_cons_law.initial_condition(x), mesh_pa.x_interfaces); % Initial condition.
    b = zeros(mesh_pa.nx, mesh_pa.nt);
    b(:, 1) = u0;
    b = b(:);
    
    % Get exact solution on current mesh by time-stepping
    tstart_seq_ts = tic;
    u_time_stepping = my_cons_law_st_system.forward_solve(b);
    fprintf('Seq. ts timer: %.2f\n', toc(tstart_seq_ts))

    
    %% Initialize solution as random perturbation of exact solution
    if strcmp(initial_iterate, 'rand-pert-of-exact-sol')
        u_init = u_time_stepping + rand_pert_mag*randn(size(u_time_stepping));

        % Initialize iterate with initial condition evaluated on current grid!
        u_init(:, 1) = u0;
    end
    
    %% Initialize solution using nested iteration
    % Get exact solution with time-stepping and then skip to next
    % resolution.
    if nx_idx == 1
        u = u_time_stepping;
        continue
        
    % Take the initial iterate as a small perturbation of the exact solution. 
    elseif nx_idx > 1 
        
        % Get space-time interpolation matrix.
        I = bilinear_space_time_interpolation_FV(mesh_pa.nx, mesh_pa_coarse.t, mesh_pa.t);
        
        if strcmp(initial_iterate, 'nested-iteration-ts') || coarse_solve_diverged
            u_init = I*u_time_stepping(:);
        elseif strcmp(initial_iterate, 'nested-iteration-it')
            u_init = I*u_coarse(:);
        end
        u_init = reshape(u_init, [mesh_pa.nx, mesh_pa.nt]);
        
        % Initialize iterate with initial condition evaluated on current grid!
        u_init(:, 1) = u0; 
        
        u = u_init(:);
    end
    
    %% Compute initial residual. 
    [u, r] = my_cons_law_st_system.relaxation(u, b); 
    resnorm_array = norm(r, resnorm);
    errnorm_array = norm(u(:) - u_time_stepping(:), resnorm);
    % Print residual norm.
    fprintf('it %d: ||r_0|| = %.2e\n', 0, resnorm_array(1))
    
    
    %% Setup MGRIT object
    if strcmp(inner_solve_pa.linear_solve, 'MGRIT')
    myMGRIT_object = struct();
    myMGRIT_object.BE_coarse_solve.id      = BE_coarse_solver;
    myMGRIT_object.BE_coarse_solve.maxiter = 10;
    myMGRIT_object.BE_coarse_solve.rtol    = 0.000001;

    myMGRIT_object.t = mesh_pa.t;
    myMGRIT_object.block_size = mesh_pa.nx;
    end

    
    %% Iterate Richardson 
    for iter_outer = 1:maxit_outer

        % Check if halting tolerance reached
        if resnorm_array(end)/resnorm_array(1) < rtol_outer
            break
        % Check if diverged
        elseif resnorm_array(end)/resnorm_array(1) > diverge_tol_outer 
            solve_diverged = true;
            error('rich iteration has diverged...')
            break
        end

        % Solve linearized problem: P*e = r.
        tstart_lin_solve = tic;
        if strcmp(inner_solve_pa.linear_solve, 'direct')
            e = my_linearized_cons_law_st_system.forward_solve(r);

        elseif strcmp(inner_solve_pa.linear_solve, 'MGRIT')
            %e = zeros(size(r)); % Initial iterate is zero vector
            e = r; % Initial iterate is the right-hand side vector
            myMGRIT_object.e = reshape(e, [mesh_pa.nx, mesh_pa.nt]);
            [e, mgrit_resnorm_array, myMGRIT_object] = ...
                mgrit(e, r, @(a, b, c) step_cons_law_linearized_MGRIT(a, b, c, my_cons_law), myMGRIT_object, myMGRIT_solver_params);
        end
        fprintf('lin. solve timer: %.2f\n', toc(tstart_lin_solve))

        % Add linearized correction to update q.
        u = u + e;

        % Nonlinear relaxation
        [u, r] = my_cons_law_st_system.relaxation(u, b); 
        

        % Compute residual and error norms
        resnorm_array = [resnorm_array; norm(r, resnorm)];
        errnorm_array = [errnorm_array; norm(u(:) - u_time_stepping(:), resnorm)];

        % Print residual norm.
        fprintf('it %d: ||r_%d||/||r_0|| = %.2e, ||r_%d||/||r_%d|| = %.2f\n', ...
            iter_outer, iter_outer, resnorm_array(end)/resnorm_array(1), ...
            iter_outer, iter_outer-1, resnorm_array(end)/resnorm_array(end-1))
        
        
        % Plot space-time error
        if show_st_alg_error
            if iter_richardson == 0; caxis0 = []; cb_limits0 = []; end
            [caxis0, cb_limits0] = plot_st_alg_error(u_bar, u_time_stepping, mesh_pa, iter_outer, figno, caxis0, cb_limits0);
            pause
        end
    end
    % End of Rich iter loop

    nameval = lineprops(nx_idx-1);
    nameval{4} = plot_pa.ls;

    %% Plots: Richardson convergence
    % Residuals
    figure(figno+4)
    if resnorm == 1
        resnorm_array = resnorm_array * (mesh_pa.h*mesh_pa.dt);
    elseif resnorm == 2
        resnorm_array = resnorm_array * sqrt(mesh_pa.h*mesh_pa.dt);
    end
    displayname = sprintf('$n_x=%d$', mesh_pa.nx); 
    semilogy(0:numel(resnorm_array)-1, resnorm_array/resnorm_array(1), nameval{:}, ...
        'LineWidth', 2, ...
        'DisplayName', displayname)
    hold on
    
    % Error
    if show_iter_error
        figure(figno+5)
        if errnorm == 1
            errnorm_array = errnorm_array * (mesh_pa.h*mesh_pa.dt);
        elseif errnorm == 2
            errnorm_array = errnorm_array * sqrt(mesh_pa.h*mesh_pa.dt);
        end
        semilogy(0:numel(errnorm_array)-1, errnorm_array, nameval{:}, ...
            'LineWidth', 2, ...
            'DisplayName', sprintf('$n_x=%d$', mesh_pa.nx))
        hold on
    end
    
end
toc(tstart_nx_loop)
% End of looping over different mesh resolutions.

%% Pretty Richardson conv. plots
figure(figno+4)
[plot_title, fig_name] = rich_strings(my_cons_law, disc_pa, pde_pa, outer_solve_pa, linearization_pa, plot_pa);
title(plot_title)
lh = legend();
lh.set('Location', plot_pa.lh_loc)
xlabel('$k$')
ylabel('$\Vert \mathbf{r}_k \Vert / \Vert \mathbf{r}_0 \Vert$')
%ylim([10^-10, 2])
%axis([0 20 1e-16 1])
axis tight
%xlim([0 maxit_richardson])
ylim([rtol_outer, 1e1])
box on

% Save plot of residual history
if save_res_fig
    figure(figno+4);
    figure_saver(gcf, sprintf('%s/%s', plot_pa.fig_dir, fig_name), save_matlab_fig); 
end

if show_iter_error
    figure(figno+5)
    title(plot_title)
    lh = legend();
    lh.set('Location', 'SouthWest')
    xlabel('$k$')
    ylabel('$\Vert \mathbf{e}_k \Vert$')
    % ylim([10^-10, 2])
    % axis([0 20 1e-16 1])
    xlim([0 maxit_outer])
    ylim([rtol_outer, 1])
    axis tight
    box on
end

% Plot the solution
[con_fh, cs_hf, figno] = plot_cons_law_scalar_solutions(u, my_cons_law, pde_pa, disc_pa, mesh_pa, figno, plot_pa);