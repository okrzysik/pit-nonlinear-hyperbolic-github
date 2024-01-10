% Use MGRIT to solve the conservative advection problem
%   u_t + (alpha(x, t)*u)_x = 0.
%
% Fine grid discretization is conservative (linear) MOL. 
% Coarse grid is modified conservative SL. 

tmp = matlab.desktop.editor.getActive;
cd(fileparts(tmp.Filename));

clear
clc


%% PDE parameters

wave_speed_id = 1;
%wave_speed_id = 2;
wave_speed_id = 3;
wave_speed_id = 4;
wave_speed_id = 5;

xmin = -1;
xmax = 1;
tmax = 4;

%% Discretization parameters

nx_array = 2.^(6:8);

num_flux_id = 'GLF'; 
%num_flux_id = 'LLF'; 

spatial_order = 1;
spatial_order = 3;
%spatial_order = 5;

CFL_number  = 0.8; 

reconstruction_id = 'linear';
pde_id = 'linear'; 
u0_id = 1; 


%% MGRIT parameters
%m = 2; % Coarsening factor
m = 8; % Coarsening factor
%m = 16; % Coarsening factor
% m = 32; % Coarsening factor

% Choose how coarse BE systems are solved
BE_coarse_solver = 'GMRES'; % Further parameters are set below.
%BE_coarse_solver = 'LU';
%BE_coarse_solver = 'NONE'; % This means no correction to the SL scheme will be added.

% Maximum number of levels in MGRIT hierarchy
%MGRIT_maxlevels = 2;
MGRIT_maxlevels = 10;

% Maximum number of MGRIT iterations
MGRIT_maxiter   = 20; 


%% Plotting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
save_res_fig  = ~true; % Save plot of residual history
save_sol_fig  = ~true; % Save plot of sol
save_wave_fig = ~true; % Save plot of wave-speed

plot_pa.fig_dir = './figures/paper/linear/';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figno = 18;

plot_pa.show_legend = true; 
plot_pa.y_label     = true;  

% Line style
if strcmp(num_flux_id, 'LLF'); plot_pa.ls = '--'; end
if strcmp(num_flux_id, 'GLF'); plot_pa.ls = '-'; end


%% Set up coarse-grid SL method
% Order of the ODE integrator used to locate departure points. 
SL_ERK_order = spatial_order; 
% How depart points on coarse levels are located. 
depart_coarse_strategy = 'backtrack+interp'; % Note that if this is the case then depart points on the fine level are estimated with a forward Euler step
%depart_coarse_strategy = 'finest-step-ERK';
% depart_coarse_strategy = 'coarse-step-ERK';
switch SL_ERK_order
    case 1
        butcher = butcher_table('FE');
    case 3
        butcher = butcher_table('SSPERK3');
    case 5
        butcher = butcher_table('RK5');  
    otherwise
        error('No ERK method implemented for ODE_order=%d', SL_ERK_order)
end

% This handle is parsed to the time-stepping operators. 
depart_solver_ERK = @(tspan, arrive) ERK_solver(@(t, x) wave_speed(x, t), tspan, arrive, butcher, 'final');

    
%% Package everything into a struct to pass to mgrit --- %
myMGRIT_solver_params = struct();
myMGRIT_solver_params.cf = m;
myMGRIT_solver_params.maxlevels = MGRIT_maxlevels;
myMGRIT_solver_params.maxiter = MGRIT_maxiter;
myMGRIT_solver_params.final_F_relax = ~true;
myMGRIT_solver_params.pre_relax = 'FCF';
myMGRIT_solver_params.res_halt_tol = 1e-10;
myMGRIT_solver_params.min_coarse_nt = 2;

myMGRIT_object = struct();
myMGRIT_object.depart_solver_ERK = depart_solver_ERK;
myMGRIT_object.depart_coarse_strategy = depart_coarse_strategy;

myMGRIT_object.BE_coarse_solve.id      = BE_coarse_solver;
myMGRIT_object.BE_coarse_solve.maxiter = 10;
myMGRIT_object.BE_coarse_solve.rtol    = 0.01;


%% Package mesh-resolution-independent parameters
% Package mesh params
mesh_pa.xmin = xmin;
mesh_pa.xmax = xmax;
mesh_pa.tmax = tmax;

% Package disc params
disc_pa.spatial_order         = spatial_order;
disc_pa.CFL_number            = CFL_number;
disc_pa.reconstruction_id     = reconstruction_id;
disc_pa.num_flux_id           = num_flux_id;

% Create PDE object.
pde_pa.ic_id = u0_id;
my_cons_law = cons_lin_scalar(pde_pa, wave_speed_id);

% Get max initial wave-speed
my_cons_law.compute_and_set_zero_extremal_values(xmin, xmax);


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% --- Solve problem at different spatial resolutions --- %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(figno)
%hold on
for nx_idx = 1:numel(nx_array)
    tic
    
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
    my_cons_law.disc_pa        = disc_pa;
    my_cons_law.mesh_pa        = mesh_pa;
    my_cons_law.reconstruction = reconstruction;

    %% Setup MGRIT
    myMGRIT_object.t          = mesh_pa.t;
    myMGRIT_object.block_size = mesh_pa.nx;
    
    % ST initial condition and RHS vector
    u       = randn(mesh_pa.nx, mesh_pa.nt); 
    u(:, 1) = cell_average(@(x) my_cons_law.initial_condition(x), mesh_pa.x_interfaces); % Initial condition.
    g       = zeros(mesh_pa.nx, mesh_pa.nt); 
    g(:, 1) = u(:, 1);
    u       = u(:);
    g       = g(:);
    
    %% MGRIT Solve 
    [u, rnorm, myMGRIT_object] = mgrit(u, g, @(a, b, c) step_cons_lin_MGRIT(a, b, c, my_cons_law), myMGRIT_object, myMGRIT_solver_params);

    
    %% Plot residuals from this solve
    nameval = lineprops(nx_idx);
    nameval{4} = plot_pa.ls;

    %% Plots: Richardson convergence
    % Residuals
    figure(figno)
    if plot_pa.show_legend
        displayname = sprintf('$(n_x,n_t)=(%d,%d)$', mesh_pa.nx, mesh_pa.nt); 
        semilogy(0:numel(rnorm)-1, rnorm/rnorm(1), nameval{:}, ...
            'LineWidth', 2, ...
            'DisplayName', displayname)
    else
        semilogy(0:numel(rnorm)-1, rnorm/rnorm(1), nameval{:}, ...
            'LineWidth', 2, ...
            'HandleVisibility', 'Off')
    end
    hold on
  
toc    
end
% End of looping over different spatial resolutions.



%%
%%%%%%%%%%%%%%%%%
% --- Plots --- %
%%%%%%%%%%%%%%%%%

%% MGRIT residual plots
if plot_pa.show_legend
    lh = legend();
    lh.set('Location', 'best')
end
if plot_pa.y_label
    ylabel('$\Vert \mathbf{r}_k \Vert / \Vert \mathbf{r}_0 \Vert$')
else
    %yticklabels([])
    set(gca,'Yticklabel',[])
end
xlabel('$k$')
axis tight
ylim([1e-10, 1])

if strcmp(BE_coarse_solver, 'NONE')
   ylim([1e-10, 1e5]) 
end

% Save plot of residual history
if save_res_fig
    figure(figno);
    
    % Basic name of fig. 
    base_fig_name = sprintf('%s/cons_lin_MGRIT_MOL-a%d-p%d-q%d-RES', ...
                            plot_pa.fig_dir, wave_speed_id, spatial_order, MOL_RK_order);

    % MGRIT parameters
    mgrit_details = sprintf('m%d-%s-%d', myMGRIT_solver_params.cf, ...
                                         myMGRIT_object.BE_coarse_solve.id, ...
                                         myMGRIT_solver_params.maxlevels);

    fig_name = sprintf('%s-%s', base_fig_name, mgrit_details);
    figure_saver(gcf, fig_name); 
end

%% Solution plot
if save_sol_fig
    figure(figno+1)
    [X, T] = meshgrid(mesh_pa.x_centers, mesh_pa.t);
    U = reshape(u, [mesh_pa.nx, mesh_pa.nt])';
    mesh(X, T, U); view(2)
    axis tight
    box on
    c = colorbar();
    
    title('$e(x, t)$')
    xlabel('$x$')
    ylabel('$t$')
    
    if save_sol_fig
        fig_name = sprintf('%s/cons_lin_MGRIT_MOL-a%d-p%d-q%d-SOL', ...
                            plot_pa.fig_dir, wave_speed_id, spatial_order, MOL_RK_order);
        figure_saver(gcf, fig_name);
    end
end

%% Wave-speed plot
if save_wave_fig
    figure(figno+1)
    
    % If wave-speed varies in time then plot it as such
    if norm(wave_speed(mesh_pa.x_centers, rand(1)) - wave_speed(mesh_pa.x_centers, rand(1)), inf) > 1e-10        
        [X, T] = meshgrid(mesh_pa.x_centers, mesh_pa.t);
        ALPHA = wave_speed(X, T);
        mesh(X, T, ALPHA); view(2)
        axis tight
        box on
        c = colorbar();
        
        title('$\alpha(x, t)$')
        xlabel('$x$')
        ylabel('$t$')
    
    % If it varies in space only then just do a 1D plot.
    else
        
        alpha = wave_speed(mesh_pa.x_centers, []);
        plot(mesh_pa.x_centers, alpha, 'b', 'LineWidth', 2)
        axis tight
        box on        
        title('$\alpha(x)$')
        xlabel('$x$')
    end
    
    fig_name = sprintf('%s/cons_lin_MGRIT_MOL-a%d-ALPHA', ...
                            plot_pa.fig_dir, wave_speed_id);
    figure_saver(gcf, fig_name);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% ----- HELPER FUNCTIONS ----- %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Save figure
function figure_saver(fig, fig_name)
    fig.PaperPositionMode = 'auto';
    fig_pos = fig.PaperPosition;
    fig.PaperSize = [fig_pos(3) fig_pos(4)];    
    set(gcf, 'Color', 'w'); % Otherwise saved fig will have grey background
    export_fig(strcat(fig_name, '.png'), '-m4')
end

