% Time-stepping on SWE or Euler equations.
%
% Just uncomment below the PDE+domain+initial condition combination below 
% that you want to solve and then run the script.
%
% Upon completion there's the option to plot various quantities relating to
% the solution (see the plotting options below).

tmp = matlab.desktop.editor.getActive;
cd(fileparts(tmp.Filename));

clc
clear
close all

%% Plotting options
figno = 2;

fig_dir = 'figures/';
PDE_name = 'swe';
save_sol_figs = ~true;

plot_pa.ncon_lvl = 10;
plot_pa.con_tol  = 1e-3;
plot_pa.invisible_col_bar = ~true;

plot_pa.plot_solution_contours         = true;
plot_pa.plot_solution_cross_sections   = true;
plot_pa.plot_wave_speed_contours       = ~true;
plot_pa.plot_wave_speed_cross_sections = ~true;


%% Discretization parameters
nx = 2^9;

disc_pa.num_flux_id = 'LLF';
disc_pa.num_flux_id = 'ROE'; disc_pa.delta_smoothing_parameter = 1e-6;


%% PDE, domain and initial condition parameters
%% SWE: Initial depth perturbation with Gaussian of size epsilon
% pde_pa.pde_id = 'shallow-water';
% pde_pa.ic_id = 'idp1'; 
% pde_pa.bcs = 'periodic';

% % pde_pa.ic_epsilon = 0.1; % Only weak shocks form here.
% pde_pa.ic_epsilon = 0.6; % Shocks form.
% 
% mesh_pa.xmin = -5;
% mesh_pa.xmax =  5;
% disc_pa.CFL_number = 0.8; % max-wave-speed * dt/h
% mesh_pa.tmax = 10;


%% SWE: Initial cosine pert to h
% pde_pa.pde_id = 'shallow-water';
% pde_pa.ic_id = 'idp2'; 
% 
% mesh_pa.xmin = -5;
% mesh_pa.xmax =  5;
% disc_pa.CFL_number = 0.8; % max-wave-speed * dt/h
% mesh_pa.tmax = 10;
% 
% pde_pa.ic_epsilon = 0.2; % Shocks don't form, but there is steepening
% pde_pa.ic_epsilon = 0.6; 


%% SWE: Dam break problem: h0 has a jump of height eps.
% eps = 2 is the usual dam break problem from LeVeque p. 259
% pde_pa.pde_id = 'shallow-water';
% pde_pa.ic_id  = 'dam-break'; 
% pde_pa.bcs    = 'constant';
% 
% pde_pa.ic_epsilon = 0.1;
% pde_pa.ic_epsilon = 2;
% 
% mesh_pa.xmin = -10;
% mesh_pa.xmax =  10;
% disc_pa.CFL_number = 0.7; % max-wave-speed * dt/h
% mesh_pa.tmax = 5;


%% Euler: Smooth initial density and pressure perturbation
% pde_pa.pde_id = 'euler'; 
% pde_pa.ic_id = 'idp1'; 
% 
% pde_pa.ic_epsilon = 0.2;
% pde_pa.ic_epsilon = 1.2;
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
% pde_pa.ic_epsilon = 0.875; % This is the original Sod problem


    
%% Set up mesh
% Create spatial mesh
mesh_pa.nx = nx;
[mesh_pa.x_interfaces, mesh_pa.x_centers, mesh_pa.h] = space_mesh(mesh_pa.xmin, mesh_pa.xmax, mesh_pa.nx);

% Get initial conditon.
if strcmp(pde_pa.pde_id, 'shallow-water')
    my_cons_law = swe_system(pde_pa);
    [h0, u0] = my_cons_law.initial_condition(mesh_pa.x_centers);
    q0 = [h0; h0.*u0];
    
elseif strcmp(pde_pa.pde_id, 'euler')
    my_cons_law = euler_system(pde_pa);
    [rho0, u0, E0] = my_cons_law.initial_condition(mesh_pa.x_centers);
    q0 = [rho0; rho0.*u0; E0];
    
end

% Compute fastest wave-speed associated with initial condition
my_cons_law.mesh_pa = mesh_pa; % The wave-speed function requires this to be set
lambda0 = my_cons_law.wave_speeds(q0);
abs_fprime0_max = max(abs(lambda0(:)));

% Get temporal mesh.
[mesh_pa.t, mesh_pa.dt] = time_mesh(mesh_pa.h, disc_pa.CFL_number, abs_fprime0_max, 1, mesh_pa.tmax);
mesh_pa.nt = numel(mesh_pa.t);

%% Finalize construction of PDE class
my_cons_law.disc_pa = disc_pa; % Set this again so that it's updated.
my_cons_law.mesh_pa = mesh_pa;

nx = mesh_pa.nx;
nt = mesh_pa.nt;

fprintf('\nnx=%d, nt=%d\n', mesh_pa.nx, mesh_pa.nt);

% Initialize solution vector
q = zeros(my_cons_law.m * nx, nt);

% Insert initial condition
q(:, 1) = q0;

%% Time-step!
for n = 1:mesh_pa.nt-1
    q(:, n+1) = my_cons_law.step(n, q(:, n));
end
fprintf('Time-stepping done...\n')

%% Plot solution, etc.
[fh, figno] = plot_cons_law_system_solutions(q, my_cons_law, pde_pa, disc_pa, mesh_pa, figno, plot_pa);