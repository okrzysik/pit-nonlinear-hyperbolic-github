function Pxt = bilinear_space_time_interpolation_FV(nx_fine, t_coarse, t_fine, periodic_BCs)

if nargin == 3
    periodic_BCs = true;
end

Px = linear_interpolation_cell_average(nx_fine, periodic_BCs);
Pt = linear_interpolation_non_nested(t_coarse, t_fine);

% The 2-dimensional operator is the kronecker product 
Pxt = kron(Pt, Px); 

end