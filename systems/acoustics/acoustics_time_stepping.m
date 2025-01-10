% Solves the acoustic equations:
%   [p,u]_t + A(x)*[p,u]_x = 0, A = [0, K(x); 1/rho(x), 0].
%
% A 1st-order accurate Godunov discretization is implemented. 
%
% This code generates plots of the material parameters, cross-sections
% of the solution, and space-time contours of the solution. 
%


tmp = matlab.desktop.editor.getActive;
cd(fileparts(tmp.Filename));

clc
clear 
%close all

%% Plotting options
figno = 19;
save_figs = ~true;
fig_dir  = 'figures/paper/time-stepping/'; 

plot_contours       = true;
plot_cross_sec      = true; pa.num_cross_sec = 3;
plot_material_param = true;

%% PDE and disc parameters
nx_array = 2.^(11);

disc_pa.high_res = true; % Apply high-res corrections to disc or not.

%pde_pa.mat_param_id = 0; % c = Z = 1.
%pde_pa.mat_param_id = 1; % Bale et al. example 1.
%pde_pa.mat_param_id = 2; % Bale et al. example 2.
%pde_pa.mat_param_id = 3; % Bale et al. example 3.
%pde_pa.mat_param_id = 4; % Z and c are 1, and jump to 2 and 0.5, respectively.
%pde_pa.mat_param_id = 5; % Periodically layered medium
%pde_pa.mat_param_id = 6; % Randomly layered medium

% Extra parameters for layered media:
pde_pa.mat_param_num_layers = 16; 
%pde_pa.mat_param_homogenize = ~true;

mesh_pa.tmax       = 1;
% % tmax's for examples in Bale et al.
% mesh_pa.tmax       = 0.4; % For ex. 1
% mesh_pa.tmax       = 0.4; % For ex. 2
% mesh_pa.tmax       = 0.5; % For ex. 3

disc_pa.CFL_number    = 0.85;

% Spatial domain is set to [0,1] in Bale et al.
mesh_pa.xmin = 0;
mesh_pa.xmax = 1;


% Loop over different spatial resolutions
for nx_idx = 1:numel(nx_array)
    mesh_pa.nx = nx_array(nx_idx);
    
    %% Set up discretization
    % Create spatial mesh
    [mesh_pa.x_interfaces, mesh_pa.x_centers, mesh_pa.h] = space_mesh(mesh_pa.xmin, mesh_pa.xmax, mesh_pa.nx);
    
    % Get material parameters. Spatial disc expects these stored like this.
    [disc_pa.c_cc, disc_pa.Z_cc] = acoustics_material_parameters(mesh_pa.x_centers, pde_pa);
    pde_pa.abs_fprime_max = max(disc_pa.c_cc);
    
    % Get temporal mesh.
    [mesh_pa.t, mesh_pa.dt] = time_mesh(mesh_pa.h, disc_pa.CFL_number, pde_pa.abs_fprime_max, 1, mesh_pa.tmax);
    mesh_pa.nt = numel(mesh_pa.t);
    
    nx = mesh_pa.nx;
    nt = mesh_pa.nt;
    fprintf('\nnx=%d, nt=%d\n', nx, nt);
    
    % PDE initial condition.
    q               = randn(2*nx, nt); 
    q(1:nx, 1)      = initial_condition_p(mesh_pa.x_centers); 
    q(nx+1:2*nx, 1) = initial_condition_u(mesh_pa.x_centers); 
    tstart_seq_ts = tic;
    %% Time-step
    for tidx = 1:nt-1
        q(:, tidx+1) = step_acoustics(q(:, tidx), mesh_pa, disc_pa);
    end
    fprintf('Seq. ts timer: %.2f\n', toc(tstart_seq_ts))
    
    
end
% End of looping over different spatial resolutions.


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- Plots of solution etc --- %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if plot_material_param
    mat_pa_fig = acoustics_plot_material_parameters(figno, mesh_pa, disc_pa, true); figno = figno+1;
end


[mesh_pa.X, mesh_pa.T] = meshgrid(mesh_pa.x_centers, mesh_pa.t);
pa.disc_pa = disc_pa;
pa.mesh_pa = mesh_pa;
pa.ncon_lvl = 10;
pa.con_tol  = 1e-3;
pa.invisible_col_bar = ~true;

P = q(1:nx, :);
U = q(nx+1:2*nx, :);

mesh_nx_tol = 2*1024;
if mesh_pa.nx > mesh_nx_tol
    pa.meshing_down_sample_factor = mesh_pa.nx / mesh_nx_tol;
end

if plot_contours
    p_con_fig = contour_plot(figno, P, '$p(x, t)$', pa); figno = figno+1;
    u_con_fig = contour_plot(figno, U, '$u(x, t)$', pa); figno = figno+1;
end

if plot_cross_sec
    p_cs_fig = cross_sec_plot(figno, P, '$p(x, t)$', pa); figno = figno+1;
    u_cs_fig = cross_sec_plot(figno, U, '$u(x, t)$', pa); figno = figno+1;
end


if save_figs 
    % p
    fig_name = sprintf('acoustics-ex%d-p-nx%d', pde_pa.mat_param_id, mesh_pa.nx);
    figure_saver(figure(p_con_fig), strcat(fig_dir, fig_name), ~true);
    % c and Z
    fig_name = sprintf('acoustics-ex%d-mat-pa', pde_pa.mat_param_id);
    figure_saver(figure(mat_pa_fig), strcat(fig_dir, fig_name), ~true);
end

% See Bale et al. eq. (5.5)
function p0 = initial_condition_p(x)
    p0    = ones(size(x));
    I     = find(x > 0.4 & x < 0.6);
    p0(I) = 7/4 - 3/4*cos( 10*pi*x(I) - 4*pi );
end

function u0 = initial_condition_u(x)
    u0    = zeros(size(x));
end


function q1 = step_acoustics(q0, mesh_pa, disc_pa)
    
    q1 = q0 + mesh_pa.dt * acoustics_spatial_discretization_godunov(q0, mesh_pa, disc_pa);
    
    if norm(q1, inf) > 10*norm(q0, inf)
        error('step: q1 > 10*q0')
    end
end