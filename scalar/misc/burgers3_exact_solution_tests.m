% Make a few plots to inspect the exact solutions are implemented properly,
% and also look at the difference between the point-wise solution evaluated
% in cell centers and the cell-averaged solution. For a smooth problem
% these things agree to second order. In this problem, where the solution
% is smooth they agree, but where they are not smooth they do not. 

nx = 2^5;
tmax = 3.9;

xmin = -1;
xmax =  1;
CFL_number = 0.8;

% Package mesh params
mesh_pa.nx = nx;
mesh_pa.xmin = xmin;
mesh_pa.xmax = xmax;
mesh_pa.tmax = tmax;
disc_pa.CFL_number = CFL_number;
pde_pa.abs_fprime_max = 1;
disc_pa.spatial_order = 1;

% Create spatial mesh
[mesh_pa.x_interfaces, mesh_pa.x_centers, mesh_pa.h] = space_mesh(mesh_pa.xmin, mesh_pa.xmax, mesh_pa.nx);

% Get temporal mesh.
[mesh_pa.t, mesh_pa.dt] = time_mesh(mesh_pa.h, disc_pa.CFL_number, pde_pa.abs_fprime_max, disc_pa.spatial_order, mesh_pa.tmax);
mesh_pa.nt = numel(mesh_pa.t);


%% Compute point-wise solution
u_cell_centers = zeros(mesh_pa.nx, mesh_pa.nt);
for n = 1:mesh_pa.nt
    u_cell_centers(:, n) = burgers3_exact_sol_pt_wise(mesh_pa.x_centers, mesh_pa.t(n));
end
[X, T] = meshgrid(mesh_pa.x_centers, mesh_pa.t);
figure(1)
mesh(X, T, u_cell_centers')
axis tight

%% Compute cell-averaged solution
u_cell_avg = zeros(mesh_pa.nx, mesh_pa.nt);
for n = 1:mesh_pa.nt
    u_cell_avg(:, n) = burgers3_exact_sol_cell_avg(mesh_pa.x_centers, mesh_pa.t(n));
end
[X, T] = meshgrid(mesh_pa.x_centers, mesh_pa.t);
figure(2)
mesh(X, T, u_cell_avg')
axis tight

%% Look at the difference
figure(3)
mesh(X, T, (u_cell_avg - u_cell_centers)')
axis tight