% Returns the 4 time-stepping blocks needed for time-stepping the solution
% of the accoustics system from (p,u) at time t to (p,u) at time t+dt.
%
% The Godunov discretization is covered in Section 9.11 of LeVeque (2004).
%
% linearized_Z if true returns a discretization in which the stencil for 
% the ith DOF ith-row has been linearized about Z_i. I.e., the update for
% (u,p)_i only involves Z_i, and not Z_{i-1} and Z_{i+1}.
% This is a discretization I came up with when considering preconditioning 
% ideas etc. I don't know how much sense it makes. So just be careful with 
% it... (it does seem to do a pretty good job though, and it seems to be an
% excellent preconditioner for the true discretization, even when the
% impedance is wildly discontinuous!)
function [Phi_pp, Phi_pu, Phi_up, Phi_uu] = Phi_blocks_acoustics_godunov(dt, dx, c, Z, nx, linearized_Z)

    if nargin <= 5
        linearized_Z = ~true;
    end

    % Get all stencil coefficients used in the discretization.
   if ~linearized_Z
        [s_pp, s_pu, s_up, s_uu] = stencil_coefficients(dt, dx, c, Z, nx);
    else
        [s_pp, s_pu, s_up, s_uu] = stencil_coefficients_linearized_Z(dt, dx, c, Z, nx);
    end

    Phi_pp = periodic_tridiagonal_matrix(s_pp);
    Phi_pu = periodic_tridiagonal_matrix(s_pu);
    Phi_up = periodic_tridiagonal_matrix(s_up);
    Phi_uu = periodic_tridiagonal_matrix(s_uu);
    
end

% Based on how the stencils are constructed, and how spdiags builds
% matrices we need to re-organize the sub- and super-diag stencils.
% This is a quick example working out how to re-shuffle the sub- and
% super-diag stencil so that I end up with the correct thing which is: The
% first element in the sud-diagonal stencil should go in the top right of
% the matrix, and the last element in the super-diagonal stencil should go
% in the bottom left of the matrix. 
% n = 5;
% d_neg = (1:n)';
% d_zero = (n+1:2*n)';
% d_pos = (2*n+1:3*n)';
% A = full(spdiags([d_neg d_zero d_pos], -1:1, n, n))
% 
% d_neg = [d_neg(2:end); d_neg(1)];
% d_pos = [d_pos(end);   d_pos(1:end-1)];
% A = full(spdiags([d_neg d_zero d_pos], -1:1, n, n));
% A(1, end) = d_neg(end);
% A(end, 1) = d_pos(1);
% A
function A = periodic_tridiagonal_matrix(stencil)

    n = size(stencil, 1);

    d_neg  = stencil(:, 1);
    d_zero = stencil(:, 2);
    d_pos  = stencil(:, 3);

    d_neg = [d_neg(2:end); d_neg(1)];
    d_pos = [d_pos(end);   d_pos(1:end-1)];
    A = spdiags([d_neg d_zero d_pos], -1:1, n, n);
    
    % Patch up boundaries:
    A(1, end) = d_neg(end); % Top right entry
    A(end, 1) = d_pos(1);   % Bottom left entry
end


% The -1, 0, +1 connections are in the 1st, 2nd, and 3rd columns of each 
% stencil, respectively. By periodicity, the 1st entry the -1 columns
% should connect to the last variable, and the last entry in the +1 columns
% should connect to the first variable.
function [spp, spu, sup, suu] = stencil_coefficients(dt, h, c, Z, nx)
    
    spp = zeros(nx, 3);
    spu = zeros(nx, 3);
    sup = zeros(nx, 3);
    suu = zeros(nx, 3);

    for i = 1:nx

        left  = i - 1; % Index of neighbour to the LEFT of i
        right = i + 1; % Index of neighbour to the RIGHT of i

        % Fix these up if at a boundary.
        if i == 1
            left  = nx;    
        elseif i == nx
            right = 1;
        end

        % spp: p connections to p
        % -1 coefficient
        spp(i, 1) =   - dt/h * c(i) .* (  -1./(Z(left) + Z(i))  ) .* Z(i);
        % 0 coefficient
        spp(i, 2) = 1 - dt/h * c(i) .* (   1./(Z(left) + Z(i)) + 1./(Z(i) + Z(right))  ) .* Z(i);
        % 1 coefficient
        spp(i, 3) =   - dt/h * c(i) .* (  -1./(Z(i) + Z(right))  ) .* Z(i);

        % spu: p connections to u
        % -1 coefficient
        spu(i, 1) =   - dt/h * c(i) .* (  -Z(left)./(Z(left) + Z(i))  ) .* Z(i);
        % 0 coefficient
        spu(i, 2) =   - dt/h * c(i) .* (   Z(left)./(Z(left) + Z(i)) - Z(right)./(Z(i) + Z(right))  ) .* Z(i);
        % 1 coefficient
        spu(i, 3) =   - dt/h * c(i) .* (   Z(right)./(Z(i) + Z(right))  ) .* Z(i);


        % sup: u connections to p
        % -1 coefficient
        sup(i, 1) =   - dt/h * c(i) .* (  -1./(Z(left) + Z(i))  );
        % 0 coefficient
        sup(i, 2) =   - dt/h * c(i) .* (   1./(Z(left) + Z(i)) - 1./(Z(i) + Z(right))  );
        % 1 coefficient
        sup(i, 3) =   - dt/h * c(i) .* (   1./(Z(i) + Z(right))  );

        % suu: u connections to u
        % -1 coefficient
        suu(i, 1) =   - dt/h * c(i) .* (  -Z(left)./(Z(left) + Z(i))  );
        % 0 coefficient
        suu(i, 2) = 1 - dt/h * c(i) .* (   Z(left)./(Z(left) + Z(i)) + Z(right)./(Z(i) + Z(right))  );
        % 1 coefficient
        suu(i, 3) =   - dt/h * c(i) .* (  -Z(right)./(Z(i) + Z(right))  );
    end
    
end


% The -1, 0, +1 connections are in the 1st, 2nd, and 3rd columns of each 
% stencil, respectively. By periodicity, the 1st entry the -1 columns
% should connect to the last variable, and the last entry in the +1 columns
% should connect to the first variable.
function [spp, spu, sup, suu] = stencil_coefficients_linearized_Z(dt, h, c, Z, nx)
    
    spp = zeros(nx, 3);
    spu = zeros(nx, 3);
    sup = zeros(nx, 3);
    suu = zeros(nx, 3);

    for i = 1:nx

        left  = i - 1; % Index of neighbour to the LEFT of i
        right = i + 1; % Index of neighbour to the RIGHT of i

        % Fix these up if at a boundary.
        if i == 1
            left  = nx;    
        elseif i == nx
            right = 1;
        end

        
        % spp: p connections to p
        % -1 coefficient
        spp(i, 1) =     dt/h * c(i) / 2;
        % 0 coefficient
        spp(i, 2) = 1 - dt/h * c(i);
        % 1 coefficient
        spp(i, 3) =     dt/h * c(i) / 2;

        % spu: p connections to u
        % -1 coefficient
        spu(i, 1) =     dt/h * c(i) .* Z(i) / 2;
        % 0 coefficient
        spu(i, 2) = 0;
        % 1 coefficient
        spu(i, 3) =   - dt/h * c(i) .* Z(i) / 2;


        % sup: u connections to p
        % -1 coefficient
        sup(i, 1) =     dt/h * c(i) / 2 ./ Z(i);
        % 0 coefficient
        sup(i, 2) = 0;
        % 1 coefficient
        sup(i, 3) =   - dt/h * c(i) / 2 ./ Z(i);

        % suu: u connections to u
        % -1 coefficient
        suu(i, 1) =     dt/h * c(i) / 2;
        % 0 coefficient
        suu(i, 2) = 1 - dt/h * c(i);
        % 1 coefficient
        suu(i, 3) =     dt/h * c(i) / 2;
    end
    
end        