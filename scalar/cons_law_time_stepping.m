% Does time-stepping to solve the 1D scalar PDE u_t + f(u)_x = 0.
%
% This includes the linear conservation law where f(u) = alpha(x, t)*u for
% some prescribed function alpha.
%
% Upon completing, a space-time contour of the numerical solution will be
% shown and cross-sections of the solution will be shown at several time
% points.

tmp = matlab.desktop.editor.getActive;
cd(fileparts(tmp.Filename));
clear
clc


%% Plotting stuff
%%%%%%%%%%%%%%%%%%%%
figno = 90;
save_sol_figs = ~true;
fig_dir = './figures/paper/sol/';

plot_pa.ncon_lvl = 10;
plot_pa.con_tol  = 1e-3;
plot_pa.invisible_col_bar = ~true;

plot_pa.plot_solution_contours       = true;
plot_pa.plot_solution_cross_sections = true;


%% Discretization parameters
nx = 2^9;

spatial_order = 1;
spatial_order = 3;
%spatial_order = 5;

%reconstruction_id = 'linear';
reconstruction_id = 'WENO';

num_flux_id = 'GLF'; 
num_flux_id = 'LLF'; 

limit_reconstructions = false;
CFL_number  = 0.8;

%% Linear conservation law
% pde_id = 'linear'; 
% u0_id = 1; tmax = 4; 
% wave_speed_id = 1;
% %wave_speed_id = 2;
% wave_speed_id = 3;
% wave_speed_id = 4;
% %wave_speed_id = 5;


%% Burgers
pde_id = 'burgers'; 
%u0_id = 1; tmax = 4; 
%u0_id = 2; tmax = 4; 
u0_id = 3; tmax = 4; 

%% Buckley--Leverett
% pde_id = 'buckley-leverett'; if strcmp(num_flux_id, 'LLF'); limit_reconstructions = true; end
% u0_id = 3; tmax = 2; % Riemann problem. Two compound waves.


%% PDE parameters
pde_pa.ic_id = u0_id;
xmin = -1;
xmax = 1;

% Create PDE object.
if strcmp(pde_id, 'linear')
    my_cons_law = cons_lin_scalar(pde_pa, wave_speed_id);

elseif strcmp(pde_id, 'burgers')
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

%% Compute and setup mesh-resolution-dependent parameters
% Get spatial mesh.
mesh_pa.nx = nx;
[mesh_pa.x_interfaces, mesh_pa.x_centers, mesh_pa.h] = space_mesh(xmin, xmax, mesh_pa.nx);

% Get temporal mesh.
[mesh_pa.t, mesh_pa.dt] = time_mesh(mesh_pa.h, CFL_number, my_cons_law.f0_prime_max, disc_pa.spatial_order, mesh_pa.tmax);
mesh_pa.nt = numel(mesh_pa.t);
fprintf('nx=%d, nt=%d\n', mesh_pa.nx, mesh_pa.nt)

% Initialize class for spatial reconstructions
reconstruction = weighted_reconstruction(spatial_order, nx, reconstruction_id);

%% Finalize construction of cons_law_scalar
my_cons_law.disc_pa        = disc_pa;
my_cons_law.mesh_pa        = mesh_pa;
my_cons_law.reconstruction = reconstruction;


%% Initialize solution
% Step u from t(tidx) -> t(tidx+1)
tstart_seq_ts = tic;
u        = randn(mesh_pa.nx, mesh_pa.nt); 
u(:, 1)  = cell_average(@(x) my_cons_law.initial_condition(x), mesh_pa.x_interfaces); % Initial condition.
for tidx = 1:mesh_pa.nt-1
    u(:, tidx+1) = my_cons_law.step(tidx, u(:, tidx));
end
fprintf('Seq. ts timer: %.2f\n', toc(tstart_seq_ts))

% Plot the solution
[con_fh, cs_hf, figno] = plot_cons_law_scalar_solutions(u, my_cons_law, pde_pa, disc_pa, mesh_pa, figno, plot_pa);