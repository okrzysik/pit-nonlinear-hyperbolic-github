% Compute the cell average of a function by integration of its function
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