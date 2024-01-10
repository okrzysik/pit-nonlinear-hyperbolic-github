function [x_interfaces, x_centers, h] = space_mesh(xmin, xmax, nx)
    x_interfaces = linspace(xmin, xmax, nx+1)'; % There are nx+1 interfaces 
    h            = (xmax - xmin)/nx; % Set spatial step according to number of spatial DOFs
    x_centers    = x_interfaces(1:end-1) + h/2; % There are nx FV cells
end




