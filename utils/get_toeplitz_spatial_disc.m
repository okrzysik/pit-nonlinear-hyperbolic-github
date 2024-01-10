function L = get_toeplitz_spatial_disc(spatial_info)
%GET_TOEPLITZ_SPATIAL_DISC get Toeplitz spatial discretization matrix L.
%
%INPUT:
% spatial_info  :   STRUCT. Holds information about the spatial
%                           discretization. Must have the fields:
%                           1. nDOFs        :  ARRAY. Number of spatial degrees of 
%                                               freedom in each spatial dimension. 
%                                               1D problem is just an int, 2D problem 
%                                               is nDOFs in x-direction, then nDOFs 
%                                               in y-direction.
%                           2. h            :  ARRAY. Grid spacing in each direction.
%                                               Same format as above.
%                           3. stencils     :  CELL. Cell array, with each
%                                               entry pointing to two cell arrays containing
%                                               diagonal entries of the given direction and
%                                               then diagonal indices for that direction. 
%                           4. BCs          :  STRUCT. Must have a "type" field, if type == "periodic", 
%                                               L is cicrculant in 1D (BCCB in 2D). Otherwise, it's
%                                               just Toeplitz in 1D (BTTB in 2D).
%                                               
%
%OUTPUT:
%   L   :     MATRIX. Sparse matrix (stored in sparse format). Is Toeplitz
%                   in 1D, or BTTB in 2D.
%
%
%NOTES:
%   -Row-wise lexicographical ordering is used for spatial unknowns, i.e.,
%       x first then y in 2D.

% Number of spatial dimensions
spatial_dim = numel(spatial_info.nDOFs);

% 1 spatial dimension
if spatial_dim == 1

    nx = spatial_info.nDOFs(1);
    
    if strcmp(spatial_info.BCs.type, 'periodic') == 1
        periodic = 1;
        % Pass stencil entries and diagonal indices
        L = stencil2matrix(spatial_info.stencil{1}{1}, spatial_info.stencil{1}{2}, nx, periodic);
    else
        periodic = 0;
        % Pass stencil entries and diagonal indices
        L = stencil2matrix(spatial_info.stencil{1}{1}, spatial_info.stencil{1}{2}, nx, periodic);
    end
        

    
% 2 spatial dimensions. Use Kronecker sum formulation.
elseif spatial_dim == 2
    
    nx = spatial_info.nDOFs(1);
    ny = spatial_info.nDOFs(2);
    
    if strcmp(spatial_info.BCs.type, 'periodic') == 1
        periodic = 1;
        % Pass stencil entries and diagonal indices
        Lx = stencil2matrix(spatial_info.stencil{1}{1}, spatial_info.stencil{1}{2}, nx, periodic);
        Ly = stencil2matrix(spatial_info.stencil{2}{1}, spatial_info.stencil{2}{2}, ny, periodic);
    else
        periodic = 0;
        % Pass stencil entries and diagonal indices
        Lx = stencil2matrix(spatial_info.stencil{1}{1}, spatial_info.stencil{1}{2}, nx, periodic);
        Ly = stencil2matrix(spatial_info.stencil{2}{1}, spatial_info.stencil{2}{2}, ny, periodic);
    end
    
    % Form Kronecker sum of Lx and Ly, row-wise lexicographic ordering.
    L = kron(speye(ny), Lx) + kron(Ly, speye(nx));
    
else
    error('Only 1D and 2D spatial discretizations implemented')
end

end