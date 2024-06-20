function [Phi, interp_error] = semiLagrangianMatrix_new(interp_degree, nx, depart_east_neighbour, depart_east_epsilon)
% Build the SL interpolation matrix after having located departure points
%
%INPUT:
%   INTERP_DEGREE: Degree of interpolating polynomial. Non-negative
%   integer.
%   NX: Number of spatial points to interpolate. Positive integer.
%   DEPART_DELTA_EAST: Distance from depature point to arrival point
%   divided by the mesh length. Array, with value for all NX points. Or if
%   constant for all points (as in spatially-independent wave-speed case)
%
% epsilon are the distances from the departure point to its right
% neighbour.
%
% error are the errors associated with performing the polynomial
% interpolation, but have not been dived by the (interp_order+1)! term.


% if numel(depart_delta_east) ~= nx && numel(depart_delta_east) ~= 1
%     error('Difference from departure to arrival point must be provided for every point unless it is the same for all points.')
% end

N = interp_degree + 1; % Short-hand to save computing this all the time; interp_degree interpolation uses interp_degree+1 points

% Each DOF couples to interp_degree+1 DOFs, including itself
nz_per_dof = N;

% Allocate CSR structure
nz_per_row = nz_per_dof;
rows = zeros(nz_per_row * nx, 1);
cols = zeros(nz_per_row * nx, 1);
data = zeros(nz_per_row * nx, 1);



% % The departure arrival difference is the same for all points. Save some
% % time by doing interpolation once only.
% if numel(depart_delta_east) == 1
%     interp_error = zeros(1, 1);
%     [weights, nodes, interp_error(1)] = poly_interpolation(epsilon_east, interp_degree);
%     
%     % Assemble row of Phi corresponding to DOF i. Note DOF i has index i.
%     for i = 1:nx
%         % Add offset from node i to node i itself, and periodically wrap the nodes if necessary
%         rows((i-1)*nz_per_row+1:i*nz_per_row) = i;
%         cols((i-1)*nz_per_row+1:i*nz_per_row) = mod(east_distance_local(i)+i + nodes+nx-1, nx)+1; % one's based index here
%         data((i-1)*nz_per_row+1:i*nz_per_row) = weights; 
%     end
%     
%     
% % The departure arrival difference is different for every point
% else
    
    interp_error = zeros(nx, 1);
    % Assemble row of Phi corresponding to DOF i. Note DOF i has index i.
    for i = 1:nx
        [weights, nodes, interp_error(i)] = poly_interpolation(depart_east_epsilon(i), ...
                                                                interp_degree);

        % Add offset from node i to node i itself, and periodically wrap 
        % the nodes if necessary
        nodes = mod(depart_east_neighbour(i) + nodes+nx-1, nx)+1; % one's based index here

        rows((i-1)*nz_per_row+1:i*nz_per_row) = i;
        cols((i-1)*nz_per_row+1:i*nz_per_row) = nodes;
        data((i-1)*nz_per_row+1:i*nz_per_row) = weights; 
    end
%end
        
Phi = sparse(rows, cols, data, nx, nx, nz_per_row * nx);
end


