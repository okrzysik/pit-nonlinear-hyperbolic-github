% Map from primitive variables q to characteristic variables w. This
% is achieved by applying R^-1 to q.
%
% q is expected to be a spatial vector of length 2nx.
function w = left_eigenvector_map_acoustics(q, mesh_pa, disc_pa)
    nx           = mesh_pa.nx;
    w            = zeros(2*nx, 1);
    w(1:nx)      = 0.5*(-q(1:nx) ./ disc_pa.Z_cc + q(nx+1:2*nx));
    w(nx+1:2*nx) = 0.5*( q(1:nx) ./ disc_pa.Z_cc + q(nx+1:2*nx));
end


