% Apply the above stencil to advance u0 into u1.
function u1 = SL_apply_stencil(u0, d, nodes, depart_east_neighbour_trans, nx, k)

    i_star = (1:nx)' + depart_east_neighbour_trans;

    % Evaluate sum that is the integral of the reconstruction.
    j = 1;    
    flux_east     =             d(:, j).*u0(mod(i_star + nodes(j) + nx-1, nx)+1);
    for j = 2:k    
        flux_east = flux_east + d(:, j).*u0(mod(i_star + nodes(j) + nx-1, nx)+1);
    end

    % The flux also has a component that is the integral from the
    % east interface of the departure point up to the arrival
    % point. This is just the sum of the cell avergaes across all
    % these cells, but some care is needed to do this. Note that if
    % the arrival point is the east interface of the departure
    % point, i.e., i_star = i, then the intergral is over an empty
    % interval.
    mesh_inds = (1:nx)';
    I_pos     = find(i_star < mesh_inds)';
    I_neg     = find(i_star > mesh_inds)';
    for i = I_pos
        integral_cells = (i_star(i)+1):i;
        % Ensure indices are on the mesh.
        integral_cells = mod(integral_cells + nx-1, nx)+1;
        flux_east(i) = flux_east(i) + sum(u0(integral_cells)); 
    end

    for i = I_neg
        integral_cells = (i+1:i_star(i));
        % Ensure indices are on the mesh.
        integral_cells = mod(integral_cells + nx-1, nx)+1;
        flux_east(i) = flux_east(i) - sum(u0(integral_cells)); 
    end           

    % Numerical update of the solution. By the periodic boundaries, the
    % flux on the west side of the first cell is the east side flux on
    % the last cell. 
    u1 = u0 - [flux_east(1)-flux_east(end); flux_east(2:end)-flux_east(1:end-1)];

end