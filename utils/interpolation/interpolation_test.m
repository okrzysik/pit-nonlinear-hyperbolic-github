clc
clear

% Here's a code snippet demonstrating how the interpolation routines should 
% be used. It's useful
% for checking the correctness, because the solid line that MATLAB plots
% between the values of the coarse function is just a linear interpolant,
% and so our interpolated cell averages and points should lie exactly on 
% top of these line segments. 

xmin = -1;
xmax =  1;
spatial_order = 3;
%spatial_order = 5;
tmax = 2;

ut_fun = @(t) sin(pi*t);
ux_fun = @(x) cos(pi*x).^2;


nx_coarse = 2^3;
nx_fine   = 2*nx_coarse;

%% Space
[x_interfaces_coarse, x_centers_coarse, h_coarse] = space_mesh(xmin, xmax, nx_coarse);
[x_interfaces_fine,   x_centers_fine,   h_fine]   = space_mesh(xmin, xmax, nx_fine);
Px = linear_interpolation_cell_average(nx_fine);

ux_bar_coarse = cell_average(ux_fun, x_interfaces_coarse);
ux_bar_fine   = cell_average(ux_fun, x_interfaces_fine);
ux_bar_interp = Px*ux_bar_coarse;

figure()
plot(x_centers_coarse, ux_bar_coarse, '-bo', 'MarkerSize', 12, 'MarkerFaceColor', 'b', 'DisplayName', 'u bar coarse')
hold on
plot(x_centers_fine, ux_bar_fine, '--r>', 'DisplayName', 'u bar fine')
plot(x_centers_fine, ux_bar_interp, '*k', 'DisplayName', 'interpolant')
lh = legend();
title('Cell-average interpolation in space')
lh.set('Location', 'Best')

%% Time
[t_coarse, dt_coarse] = time_mesh(h_coarse, tmax, spatial_order);
[t_fine, dt_fine]     = time_mesh(h_fine, tmax, spatial_order);
Pt = linear_interpolation_non_nested(t_coarse, t_fine);
ut_coarse = ut_fun(t_coarse);
ut_fine   = ut_fun(t_fine);
ut_interp = Pt*ut_coarse;

figure()
plot(t_coarse, ut_coarse, '-bo', 'MarkerSize', 12, 'MarkerFaceColor', 'b', 'DisplayName', 'u coarse')
hold on
plot(t_fine, ut_fine, '--r>', 'DisplayName', 'u fine')
plot(t_fine, ut_interp, '*k', 'DisplayName', 'interpolant')
title('Point-wise interpolation in time')
lh = legend();
lh.set('Location', 'Best')


%% Assorted functions
% Compute the cell average of a functiot(10n by integration of its function
% handle
function u_bar = cell_average(u_fun, x_interfaces)
    nx = numel(x_interfaces) - 1;
    u_bar = zeros(nx, 1);
    for i = 1:nx
        xl = x_interfaces(i);
        xr = x_interfaces(i+1);
        h = (xr - xl);
        %u0(i) = integral(@(x) u0_fun(x) / h, xl, xr, 'AbsTol', 1.e-15, 'RelTol', 1.e-15);
        u_bar(i) = gauss_legendre_quad(u_fun, xl, xr, 7) / h;
    end
end

function [x_interfaces, x_centers, h] = space_mesh(xmin, xmax, nx)
    x_interfaces = linspace(xmin, xmax, nx+1)'; % There are nx+1 interfaces 
    h            = (xmax - xmin)/nx; % Set spatial step according to number of spatial DOFs
    x_centers    = x_interfaces(1:end-1) + h/2; % There are nx FV cells
end

function [t, dt] = time_mesh(h, tmax, spatial_order)
    
    if spatial_order == 5
        dt = 0.8 * h^(5/3); 
    else
        dt = 0.8 * h; 
    end
    nt = floor(tmax/dt)+1; % Ensure integer number of time steps (from dt = tmax/(nt+1))
    tmax_true = dt*nt;

    % Temporal grid
    t  = linspace(0, tmax_true, nt+1)'; % There are nt+1 points in between 0 and tmax_true; 
end