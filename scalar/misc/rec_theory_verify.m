% Check the theoretical error estimates for shifted linear reconstructions,
% each using k cells, and the weighted combination of them that uses 2k-1
% cells. 
%
% Everything seems in order.
%
% Note: For the k=2 case, it seemed like when I evaluated the derivative at
% the interface value I was getting an extra order of accuracy than
% expected. 

clc
clear
close all

recon_type = 'linear';
%recon_type = 'WENO';

make_recoc_plots = ~true;


p = 1;
p = 3;
%p = 5;

k = (p+1)/2;

[zeta_lhs, zeta_rhs] = reconstruction_error_weights(k);

u_handle     = my_fun(0);
u_derivative = my_fun(k); % The estimate involves the kth derivative of u.
u_weighted_derivative = my_fun(2*k-1); % The estimate involves the 2k-1st derivative of u.

nx_array = 2.^(4:7);

% There are k reconstructions at each mesh resolution nx.
error_norm_lhs = zeros(numel(nx_array), k);
error_norm_rhs = zeros(numel(nx_array), k);

% Error in the error estimate
error_error_est_norm_lhs = zeros(numel(nx_array), k);
error_error_est_norm_rhs = zeros(numel(nx_array), k);

weighted_error_norm_lhs = zeros(numel(nx_array), 1);
weighted_error_norm_rhs = zeros(numel(nx_array), 1);

error_weighted_error_est_norm_lhs = zeros(numel(nx_array), 1);
error_weighted_error_est_norm_rhs = zeros(numel(nx_array), 1);


for nx_idx = 1:numel(nx_array)
    
    nx = nx_array(nx_idx);

    % Mesh
    h                 = 2/nx;
    x_interfaces     = linspace(-1, 1, nx+1)';
    x_interfaces_lhs = x_interfaces(1:end-1);
    x_interfaces_rhs = x_interfaces(2:end);
    x_centers        = x_interfaces(1:end-1) + h/2;

    u_exact_lhs = u_handle(x_interfaces_lhs);
    u_exact_rhs = u_handle(x_interfaces_rhs);

    % Compute cell-averaged initial condition.
    u_bar = zeros(nx, 1);
    for i = 1:nx
        xl    = x_centers(i)-h/2;
        xr    = x_centers(i)+h/2;
        u_bar(i) = integral(@(x) u_handle(x) / h, xl, xr, 'AbsTol', 1.e-15, 'RelTol', 1.e-15);
        %u0(i) = gauss_legendre_quad(u0_fun, xl, xr, 7) / h;
    end

    my_recon = weighted_reconstruction(p, nx, recon_type);

    % Plot the reconstructions
    if make_recoc_plots
        figure(1)
        nameval = lineprops(1);
        plot(x_interfaces_lhs, u_exact_lhs, nameval{:}, 'DisplayName', '$u$')
        hold on

        figure(2)
        nameval = lineprops(1);
        plot(x_interfaces_rhs, u_exact_rhs, nameval{:}, 'DisplayName', '$u$')
        hold on
    end


    % Do the k different reconstructions, each with a different left-shift
    % ell.
    for ell = 0:k-1
        ell_idx = ell+1;

        u_reconstruct_lhs = my_recon.reconstruct_matrices_left{ell_idx} * u_bar;
        u_reconstruct_rhs = my_recon.reconstruct_matrices_right{ell_idx} * u_bar;

        error_lhs = u_exact_lhs - u_reconstruct_lhs;
        error_rhs = u_exact_rhs - u_reconstruct_rhs;

        error_norm_lhs(nx_idx, ell_idx) = norm(error_lhs, inf);
        error_norm_rhs(nx_idx, ell_idx) = norm(error_rhs, inf);

        % The estimate is for the derivative anywhere in the interpolation
        % interval, so just pick the middle of the cell.
        error_est_lhs = (-h)^k / factorial(k+1) * zeta_lhs(ell_idx) * u_derivative(x_interfaces_lhs + h/2);
        error_est_rhs = (-h)^k / factorial(k+1) * zeta_rhs(ell_idx) * u_derivative(x_interfaces_rhs - h/2);

        error_error_est_norm_lhs(nx_idx, ell_idx) = norm(error_lhs - error_est_lhs, inf);
        error_error_est_norm_rhs(nx_idx, ell_idx) = norm(error_rhs - error_est_rhs, inf);
        
        % Plot the reconstructions
        if make_recoc_plots
            nameval = lineprops(ell_idx+1);
            figure(1)
            plot(x_interfaces_lhs, u_reconstruct_lhs, nameval{:}, 'DisplayName', sprintf('$\\ell = %d$', ell))
            figure(2)
            plot(x_interfaces_rhs, u_reconstruct_rhs, nameval{:}, 'DisplayName', sprintf('$\\ell = %d$', ell))

            % Plot the error and the estimate of the error
            figure(3)
            hold on
            plot(x_interfaces_lhs, error_lhs, nameval{:}, 'DisplayName', sprintf('error: $\\ell = %d$', ell))
            plot(x_interfaces_lhs, error_est_lhs, nameval{:}, 'LineStyle', 'none', 'DisplayName', sprintf('error est: $\\ell = %d$', ell))

            figure(4)
            hold on
            plot(x_interfaces_rhs, error_rhs, nameval{:}, 'DisplayName', sprintf('error: $\\ell = %d$', ell))
            plot(x_interfaces_rhs, error_est_rhs, nameval{:}, 'LineStyle', 'none', 'DisplayName', sprintf('error est: $\\ell = %d$', ell))
        end
    end
    
    % Reconstructions on big stencil with 2k-1 cells
    [u_reconstruct_lhs, u_reconstruct_rhs] = my_recon.interface_reconstruction(u_bar);

    error_lhs = u_exact_lhs - u_reconstruct_lhs;
    error_rhs = u_exact_rhs - u_reconstruct_rhs;
    
    weighted_error_norm_lhs(nx_idx) = norm(error_lhs, inf);
    weighted_error_norm_rhs(nx_idx) = norm(error_rhs, inf);

    
    % The estimate is for the derivative anywhere in the interpolation
    % interval, so just pick the middle of the cell.
    error_est_lhs = (-1)^(k) * h^(2*k-1) * factorial(k) * factorial(k-1) / factorial(2*k) * u_weighted_derivative(x_interfaces_lhs + h/2);
    error_est_rhs = (-1)^(k+1) * h^(2*k-1) * factorial(k) * factorial(k-1) / factorial(2*k) * u_weighted_derivative(x_interfaces_rhs - h/2);

    error_weighted_error_est_norm_lhs(nx_idx) = norm(error_lhs - error_est_lhs, inf); 
    error_weighted_error_est_norm_rhs(nx_idx) = norm(error_rhs - error_est_rhs, inf); 

    if make_recoc_plots
        figure(1)
        lh = legend();
        title('LHS')
        xlabel('$x$')

        figure(2)
        lh = legend();
        title('RHS')
        xlabel('$x$')

        figure(3)
        lh = legend();
        title('LHS')
        xlabel('$x$')

        figure(3)
        lh = legend();
        title('RHS')
        xlabel('$x$')
    end

end





%% Print stats about convergence rates
% fprintf('Log10 error inf norms:\n')
% fprintf('LHS: ell = %s \n', num2str(0:k-1))
% disp(log10(error_norm_lhs))
% fprintf('RHS: ell = %s \n', num2str(0:k-1))
% disp(log10(error_norm_rhs))
fprintf('Rates for error:\n')
fprintf('LHS: ell = %s \n', num2str(0:k-1))
disp([log2(error_norm_lhs(1:end-1, :) ./ error_norm_lhs(2:end, :))]);
fprintf('RHS: ell = %s \n', num2str(0:k-1))
disp([log2(error_norm_rhs(1:end-1, :) ./ error_norm_rhs(2:end, :))]);


% The estimate should be one order higher
fprintf('Rates for error estimate\n')
fprintf('LHS: ell = %s \n', num2str(0:k-1))
disp([log2(error_error_est_norm_lhs(1:end-1, :) ./ error_error_est_norm_lhs(2:end, :))]);
fprintf('RHS: ell = %s \n', num2str(0:k-1))
disp([log2(error_error_est_norm_rhs(1:end-1, :) ./ error_error_est_norm_rhs(2:end, :))]);


fprintf('\n\nRates for weighted error:\n')
fprintf('LHS:\n')
disp([log2(weighted_error_norm_lhs(1:end-1) ./ weighted_error_norm_lhs(2:end))]);
fprintf('RHS:\n')
disp([log2(weighted_error_norm_rhs(1:end-1) ./ weighted_error_norm_rhs(2:end))]);

fprintf('Rates for estimate of weighted error:\n')
fprintf('LHS:\n')
disp([log2(error_weighted_error_est_norm_lhs(1:end-1) ./ error_weighted_error_est_norm_lhs(2:end))]);
fprintf('RHS:\n')
disp([log2(error_weighted_error_est_norm_rhs(1:end-1) ./ error_weighted_error_est_norm_rhs(2:end))]);

% 
% fprintf('Log10 weighted error inf norms:\n')
% fprintf('LHS: ell = %s \n', num2str(0:k-1))
% disp(log10(weighted_error_norm_lhs))
% fprintf('RHS: ell = %s \n', num2str(0:k-1))
% disp(log10(weighted_error_norm_rhs))
% fprintf('\n\nRates for weighted error:\n')
% 
% % Do non-weighted, linear reconstruction on a big stencil. Does this
% % produce exactly the same error as when we take a weighted combination of 
% % reconstructions on small stencils?
% my_recon = weighted_reconstruction(5, nx, 'linear');
% ell_idx = 2;
% u_reconstruct_lhs = my_recon.reconstruct_matrices_left{ell_idx}  * u_bar;
% u_reconstruct_rhs = my_recon.reconstruct_matrices_right{ell_idx} * u_bar;
% 
% error_lhs = u_exact_lhs - u_reconstruct_lhs;
% error_rhs = u_exact_rhs - u_reconstruct_rhs;
% 
% error_norm_lhs = norm(error_lhs, inf);
% error_norm_rhs = norm(error_rhs, inf);
% fprintf('Log10 ... error inf norms:\n')
% fprintf('LHS: ell = %s \n', num2str(0:k-1))
% disp(log10(error_norm_lhs))
% fprintf('RHS: ell = %s \n', num2str(0:k-1))
% disp(log10(error_norm_rhs))


%% Helper functions
function [zeta_left, zeta_right] = reconstruction_error_weights(k)

    zeta_left = [];
    zeta_right = [];

    for ell = 0:k-1
        zeta_left  = [zeta_left;  (-1)^(ell)*factorial(ell)*factorial(k-ell)      ];
        zeta_right = [zeta_right; (-1)^(ell+1)*factorial(ell+1)*factorial(k-ell-1)];
    end

end

% The function used in the test is just a sin. The k>=0 here is the degree
% of the derivative (this is required to evaluate the estimate).
function u_handle = my_fun(k)
    a = 1;
    if mod(k, 2) == 0
        u_handle = @(x) (-1)^(k/2)*(a*pi)^k * sin(a*pi*x);
    else 
        u_handle = @(x) (-1)^((k-1)/2)*(a*pi)^k * cos(a*pi*x);
    end
end