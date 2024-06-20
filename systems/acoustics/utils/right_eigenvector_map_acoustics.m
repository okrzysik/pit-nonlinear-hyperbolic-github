% Map from characteristic variables w to primitive variables q. This
% is achieved by applying R to w.
%
% q is expected to be a spatial vector of length 2nx.
function q = right_eigenvector_map_acoustics(w, mesh_pa, disc_pa)
    nx           = mesh_pa.nx;
    q            = zeros(2*nx, 1);
    q(1:nx)      = disc_pa.Z_cc .* (-w(1:nx) + w(nx+1:2*nx));
    q(nx+1:2*nx) =                   w(1:nx) + w(nx+1:2*nx);
end