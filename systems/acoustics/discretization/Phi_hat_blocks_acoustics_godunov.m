% Computes each of the blocks in the 2x2 time-stepping operator Phi_hat
% that corresponds to the Godunov discretization in characteristic
% variables.
%
function [Phi_hat_11, Phi_hat_12, Phi_hat_21, Phi_hat_22] = Phi_hat_blocks_acoustics_godunov(dt, dx, c, Z, nx)


    % Get all stencil coefficients used in the discretization.
    [s_11, s_12, s_21, s_22] = stencil_coefficients(dt, dx, c, Z, nx);
    
    Phi_hat_11 = periodic_tridiagonal_matrix(s_11);
    Phi_hat_12 = periodic_tridiagonal_matrix(s_12);
    Phi_hat_21 = periodic_tridiagonal_matrix(s_21);
    Phi_hat_22 = periodic_tridiagonal_matrix(s_22);
    
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
function [s11, s12, s21, s22] = stencil_coefficients(dt, h, c, Z, nx)
    
    s11 = zeros(nx, 3);
    s12 = zeros(nx, 3);
    s21 = zeros(nx, 3);
    s22 = zeros(nx, 3);

    for i = 1:nx

        left  = i - 1; % Index of neighbour to the LEFT of i
        right = i + 1; % Index of neighbour to the RIGHT of i

        % Fix these up if at a boundary.
        if i == 1
            left  = nx;    
        elseif i == nx
            right = 1;
        end

        % 11: 1 connections to 1
        % -1 coefficient
        s11(i, 1) =   0;
        % 0 coefficient
        s11(i, 2) = 1 - dt/h * c(i) .* 1;
        % 1 coefficient
        s11(i, 3) =   - dt/h * c(i) .* -2*(Z(right)./(Z(i) + Z(right)));

        % 12: 1 connections to 2
        % -1 coefficient
        s12(i, 1) =   0;
        % 0 coefficient
        s12(i, 2) =   dt/h * c(i) .* ( Z(i) - Z(right) )./( Z(i)+Z(right) );
        % 1 coefficient
        s12(i, 3) =   0;


        % 21: 2 connections to 1
        % -1 coefficient
        s21(i, 1) =   0;
        % 0 coefficient
        s21(i, 2) =   - dt/h * c(i) .* ( Z(left)-Z(i) )./( Z(left)+Z(i) );
        % 1 coefficient
        s21(i, 3) =   0;

        % 22: 2 connections to 2
        % -1 coefficient
        s22(i, 1) =   - dt/h * c(i) .* ( -2*Z(left)./(Z(left) + Z(i))  );
        % 0 coefficient
        s22(i, 2) = 1 - dt/h * c(i) .* 1;
        % 1 coefficient
        s22(i, 3) =   0;
    end
    
end
