function Px = linear_interpolation_cell_average(nx_fine, periodic_BCs)
% In space, P interpolates a cell average vector with nx/2 cells to a cell
% average vector of nx cells, where each of the nx/2 cells is split into
% two new cells. 
% BOUNDARIES: By default, periodic boundary conditions are used that assume 
% the left most interface is the same as the right most interface. However,
% if the "PERIODIC_BCS" is set to false, then the left and right boundary
% interpolants are set by constant interpolation.
%
% The interpolation is based on assuming the cell average agrees with the
% value at the cell center and the cell center values are interpolated from
% the coarse mesh to the fine mesh. Note that none of the cell centers on
% the two meshes are co-located, and hence why the stencil looks a bit
% weird compared to the FD-based interpolation stencil in which the fine
% mesh comes from just adding a point in between each two coarse mesh
% points. 

% Make periodic BCs the default
if nargin == 1
    periodic_BCs = true;
end

assert(mod(nx_fine, 2) == 0, 'nx_fine must be divisible by 2')
nx_coarse = nx_fine/2;


% Top left-hand side corner of matrix is special due to periodic BC
Px = sparse(nx_fine, nx_coarse);
Px(1, 1) = 3;
Px(1, nx_coarse) = 1;
Px(2, 1:2) = [3, 1];
stencil = [1, 3, 0; 0, 3, 1];

% Build the matrix by looping over pairs of consecutive rows and inserting
% the above stencil into the corresponding three columns (I guess there are
% many ways to build the matrix, this was the first that came to mind).
for block_row = 2:nx_coarse-1
    rows = [2*block_row-1, 2*block_row];
    cols = [block_row - 1, block_row, block_row + 1];
    Px(rows, cols) = stencil;
end

% Bottom right-hand corner of the matrix is special due to periodic BC
Px(nx_fine-1, [nx_coarse-1, nx_coarse]) = [1, 3];
Px(nx_fine, nx_fine/2) = 3;
Px(nx_fine, 1) = 1;

Px = 0.25*Px;


% Patch up matrix by re-setting entries using periodicity.
if ~periodic_BCs
    % Top LHS corner
    Px(1, 1) = 1;   % This was 3/4
    Px(1, end) = 0; % This was 1/4
    % Bottom RHS corner
    Px(end, end) = 1; % This was 3/4
    Px(end, 1) = 0;   % This was 1/4
end

end