function Px = linear_interpolation(nx_fine)

% P will interpolate a vector on a space grid with nx/2 
% points to a grid with nx points
% The fine mesh has h that is exactly half of the coarse mesh. The
% refinement happens by the addition of a fine point to the right of a
% coarse point. I.e., the LHS spatial boundary
% remain constant as the mesh is refined.
%
% Periodic boundaries are used in space. 
%
%% TODO: Be careful before blindly using this code. You need to check how
% exactly the fine grid relates to the coarse grid and if this is how you
% intend of it to be set up
%%

e_x = 0.5*ones(nx_fine,1); % Dummy matrix elements.
% Build a square matrix to extract every second column from.
Px = spdiags([e_x 2*e_x e_x], -1:1, nx_fine, nx_fine); 
% Periodic boundaries in space.
Px(1, end) = 0.5; 
Px(end, 1) = 0.5;
Px = Px(:, 1:2:nx_fine-1); % The 1-dimensional interpolation operator in space.

end