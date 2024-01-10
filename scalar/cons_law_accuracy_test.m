%% Measures discretization accuracy for solution computed by time-stepping
% The PDEs considered are
%   u_t + f(u)_x = 0,
% with
%   f(u) = u^2 / 2       -- Burgers equation
%   f(u) = alpha(x, t)*u -- Linear conservation law
%
% The exact solutions are implemented using something along the lines of
% the method-of-characteristics. 
%
% The code that computes the exact solution to Burgers equation will work
% only so long as the solution is sufficiently smooth (once characteristics
% intersect or almost intersect it breaks down).
%
% Tests can be done for 1st-, 3rd-, and 5th-order discretizations, using
% either linear or WENO reconstructions in the FV procedure. 
%
% The results indicate that roughly the correct order of convergence is 
% reached, but it can sometimes be a bit messy. Sometimes better than 
% expected, sometimes worse.
% You also have to be careful on the norm you measure in (e.g., classical 
% WENO degrades to second order at critical points, so if you measure in 
% the infinity norm you won't see this).
% Convergence rates are definitely much cleaner for the linear
% reconstructions relative to WENO reconstructions.
%
%
%% Burger's equation with square-wave initial condition.
% For Burgers there's also an exact solution implemented for the square
% wave initial condition, where the solution has a rarefaction and a shock
% wave, and these intersect at t = 1. 
% Some of the literature on Burgers equation says that convergence rate
% should be O(h), but much of this considers only problems with (possibly
% interacting) shocks, and not rarefaction waves. Considering the last section of
% the paper: "Rarefactions and Large Time Behavior for Parabolic Equations
% and Monotone Schemes" by Eduard Harabetian, the convergence rate for
% monotone schemes in the presence of rarefaction waves is O(h*log(h)).
% Numerically, if using a 1st-order reconstruction, I indeed measure a
% convergence rate of O(h*log(h)); however, using a 3rd-order WENO
% reconstruction I observe this improve to O(h). Interestingly, to get
% O(h) convergence rates with the 3rd-order WENO reconstruction, it appears 
% necessary to compare the numerical solution (which approximates a cell 
% average, since we're using FV methods) to the cell-average of the exact
% solution and not just the exact solution evaluated in cell centers!
% 
% Another interesting observation is that if a 3rd-order linear
% discretization is used then the L1-rate still seems to be O(h). The
% numerical solution has spurious oscillations along the shock, but the
% L1-norm doesn't seem to care enough about these for the rate to be
% reduced (probably because in space-time the numerical shock has a narrow 
% width that's fixed independent of h, and so errors at the shock probably
% do not overall contribute significantly to the total L1 error given that
% the point-wise error everwhere else in the domain is bigger than zero,
% and the total contribution of the error from these points increases as
% the mesh is refined since there are more and more points added).



clear
clc
figno = 18;

%% Linear conservation law
% pde_id = 'linear'; tmax = 1;
% 
% wave_speed_id = 1;
% wave_speed_id = 3;
% wave_speed_id = 4;
% %wave_speed_id = 5;
% u0_id = 1;  


%% Burgers equation
pde_id = 'burgers'; 

% Initial condition: Just use basic sine wave, and for Burgers we need to
% make sure that the time is small enough such that a shock has not yet
% developed. 
u0_id = 1; tmax = 0.3; 

% The exact solution for Burgers with a square-wave initial condition is
% implemented also.
%u0_id = 3; tmax = 3.9;
% Choose which exact solution of Burgers equation is used. The one based on
% evaluating the point-wise solution at cell centers, or properly computing
% the cell average solution. 
exact_burgers3_solution = @burgers3_exact_sol_pt_wise;
exact_burgers3_solution = @burgers3_exact_sol_cell_avg;
s

%% Discretization parameters

spatial_order = 1;
spatial_order = 3;
%spatial_order = 5;

reconstruction_id = 'linear';
%reconstruction_id = 'WENO';

num_flux_id = 'GLF'; 
num_flux_id = 'LLF'; 

limit_reconstructions = false;
CFL_number  = 0.8;

%% Mesh parameters
% Solve the problem over a sequence of meshes with increasing resolution.
nx_array = 2.^(5:10);

xmin = -1;
xmax = 1;
pde_pa.ic_id = u0_id;

%% Plotting stuff
%%%%%%%%%%%%%%%%%%%%
figno = 90;

plot_pa.ncon_lvl = 10;
plot_pa.con_tol  = 1e-3;
plot_pa.invisible_col_bar = ~true;

plot_pa.plot_solution_contours       = true;
plot_pa.plot_solution_cross_sections = true;

%% Create PDE object.
if strcmp(pde_id, 'linear')
    my_cons_law = cons_lin_scalar(pde_pa, wave_speed_id);

elseif strcmp(pde_id, 'burgers')
    my_cons_law = burgers(pde_pa);
    
end

% Get max initial wave-speed
my_cons_law.compute_and_set_zero_extremal_values(xmin, xmax);


%% Package mesh-resolution-independent parameters
% Package mesh params
mesh_pa.xmin = xmin;
mesh_pa.xmax = xmax;
mesh_pa.tmax = tmax;

% Package disc params
disc_pa.spatial_order         = spatial_order;
disc_pa.limit_reconstructions = limit_reconstructions;
disc_pa.CFL_number            = CFL_number;
disc_pa.reconstruction_id     = reconstruction_id;
disc_pa.num_flux_id           = num_flux_id;


%% Solve problems at different spatial resolutions. 
error_array_1norm   = [];
error_array_2norm   = [];
error_array_infnorm = [];
for nx_idx = 1:numel(nx_array)
    
    %% Compute and setup mesh-resolution-dependent parameters
    % Get spatial mesh.
    mesh_pa.nx = nx_array(nx_idx);
    [mesh_pa.x_interfaces, mesh_pa.x_centers, mesh_pa.h] = space_mesh(xmin, xmax, mesh_pa.nx);

    % Get temporal mesh.
    [mesh_pa.t, mesh_pa.dt] = time_mesh(mesh_pa.h, CFL_number, my_cons_law.f0_prime_max, disc_pa.spatial_order, mesh_pa.tmax);
    mesh_pa.nt = numel(mesh_pa.t);
    fprintf('nx=%d, nt=%d\n', mesh_pa.nx, mesh_pa.nt)

    % Initialize class for spatial reconstructions
    reconstruction = weighted_reconstruction(spatial_order, mesh_pa.nx, reconstruction_id);

    %% Finalize construction of cons_law_scalar
    my_cons_law.disc_pa        = disc_pa;
    my_cons_law.mesh_pa        = mesh_pa;
    my_cons_law.reconstruction = reconstruction;


    %% Initialize solution
    % Step u from t(tidx) -> t(tidx+1)
    tstart_seq_ts = tic;
    u        = randn(mesh_pa.nx, mesh_pa.nt); 
    u(:, 1)  = cell_average(@(x) my_cons_law.initial_condition(x), mesh_pa.x_interfaces); % Initial condition.
    for tidx = 1:mesh_pa.nt-1
        u(:, tidx+1) = my_cons_law.step(tidx, u(:, tidx));
    end
    fprintf('Seq. ts timer: %.2f\n', toc(tstart_seq_ts))
    
    if strcmp(pde_id, 'linear')
        u_cell_avg_exact = cons_lin_exact_sol_cell_avg(my_cons_law);
    
    elseif strcmp(pde_id, 'burgers')
        if u0_id == 1
            u_cell_avg_exact = burgers_exact_sol_cell_avg(@(x) my_cons_law.initial_condition(x), mesh_pa.x_interfaces, mesh_pa.t(end));
        elseif u0_id == 3
            u_cell_avg_exact = exact_burgers3_solution(mesh_pa.x_centers, mesh_pa.t(end));
        end
    end
    
    % Compute discrete error norm.
    error_array_1norm   = [error_array_1norm;   norm(u(:, end) - u_cell_avg_exact, 1  )*mesh_pa.h];
    error_array_2norm   = [error_array_2norm;   norm(u(:, end) - u_cell_avg_exact, 2  )*sqrt(mesh_pa.h)];
    error_array_infnorm = [error_array_infnorm; norm(u(:, end) - u_cell_avg_exact, inf)];
    

end
% End of looping over all meshes

% Plot the solution
[con_fh, cs_hf, figno] = plot_cons_law_scalar_solutions(u, my_cons_law, pde_pa, disc_pa, mesh_pa, figno, plot_pa);
xlabel('$x$')

% Compute exact solution of PDE at final time and plot this over the top.
if strcmp(pde_id, 'linear')
    u_cell_avg_exact = cons_lin_exact_sol_cell_avg(my_cons_law);
elseif strcmp(pde_id, 'burgers')
    if u0_id == 1
        u_cell_avg_exact = burgers_exact_sol_cell_avg(@(x) my_cons_law.initial_condition(x), mesh_pa.x_interfaces, mesh_pa.t(end));
    elseif u0_id == 3
        u_cell_avg_exact = exact_burgers3_solution(mesh_pa.x_centers, mesh_pa.t(end));
    end
end
plot(mesh_pa.x_centers, u_cell_avg_exact, ...
    'ro', 'MarkerSize', 10, ...
    'DisplayName', 'exact')
lh = legend();
lh.set('Location', 'Best')
box on


%% Print stats about the error.
fprintf('Log10 discretization error in discrete 1, 2, and inf norms:\n')
disp([log10(error_array_1norm), log10(error_array_2norm), log10(error_array_infnorm)])
fprintf('Discretization order in discrete 1, 2, and inf norms:\n')
disp([log2(error_array_1norm(1:end-1)  ./error_array_1norm(2:end)), ...
      log2(error_array_2norm(1:end-1)  ./error_array_2norm(2:end)), ...
      log2(error_array_infnorm(1:end-1)./error_array_infnorm(2:end))])
  
  
%% Error plots
figure(figno)
if strcmp(disc_pa.reconstruction_id, 'WENO')
    mk = '--o';
elseif strcmp(disc_pa.reconstruction_id, 'linear')
    mk = '-x';
end

mycols  = {'k', 'r', 'b', [0, 0.5, 0], [0.75, 0, 0.95], [0.9290, 0.6940, 0.1250], [0.3010 0.7450 0.9330]};


%% For smooth initial condition plot expected rate and one order less
if u0_id == 1

    semilogy(log2(nx_array), error_array_1norm, mk, ...
        'Color', mycols{1}, 'LineWidth', 2, 'DisplayName', sprintf('$\\ell = 1$, %s', disc_pa.reconstruction_id))
    hold on
    % semilogy(log2(nx_array), error_array_2norm, ls, 'Marker', mk, ...
    %     'Color', mycols{2}, 'LineWidth', 2, 'DisplayName', sprintf('$\\ell = 2$, %s', reconstruction_type))
    semilogy(log2(nx_array), error_array_infnorm, mk,  ...
        'Color', mycols{3}, 'LineWidth', 2, 'DisplayName', sprintf('$\\ell = \\infty$, %s', disc_pa.reconstruction_id))

    semilogy(log2(nx_array), 0.7*error_array_infnorm(end-1)*(nx_array/nx_array(end-1)).^-disc_pa.spatial_order, ...
        'Color', mycols{4}, 'DisplayName', sprintf('$h^{%d}$', disc_pa.spatial_order));
    semilogy(log2(nx_array), 0.7*error_array_infnorm(end-1)*(nx_array/nx_array(end-1)).^-(disc_pa.spatial_order-1), ...
        'Color', mycols{5}, 'DisplayName', sprintf('$h^{%d}$', disc_pa.spatial_order-1));

%% For non-smooth problem, plot rates h and h*log(h)
else
    semilogy(log2(nx_array), error_array_1norm, mk, ...
        'Color', mycols{1}, 'LineWidth', 2, 'DisplayName', sprintf('$\\ell = 1$, %s', disc_pa.reconstruction_id))
    hold on
    
    h_array = 1./nx_array;
    
    semilogy(log2(nx_array), 0.7*error_array_1norm(end-1) * h_array./h_array(end-1), ...
        'Color', mycols{4}, 'DisplayName', '$h$');
    semilogy(log2(nx_array), 1.5*error_array_1norm(end-1) * h_array./h_array(end-1) .* abs(log(h_array)) / abs(log(h_array(end-1))), ...
        'Color', mycols{5}, 'DisplayName', '$h |\log h|$');
    
end

lh = legend();
xlabel('$\log_2(n_x)$')
ylabel('$\Vert \mathbf{e} \Vert_{\ell}$')
ax = gca;
ax.XTick = unique( round(ax.XTick) );
axis tight
box on

if strcmp(pde_id, 'linear')
    title_str = sprintf('%s ($\\alpha=%d$), $p=%d$', pde_id, wave_speed_id, spatial_order);
else
    title_str = sprintf('%s, $p=%d$', pde_id, spatial_order);
end
title(title_str)



%% Space-time discretization error plot for non-smooth Burgers problem
if mesh_pa.nx < 1025 && u0_id == 3 && strcmp(pde_id, 'burgers')
    error_mask_tol = 10^-8; % Errors < tol are masked to zero to help with visualization

    u_cell_avg_exact = zeros(size(u));
    for n = 1:mesh_pa.nt
        u_cell_avg_exact(:, n) = exact_burgers3_solution(mesh_pa.x_centers, mesh_pa.t(n));
    end
    disc_error = u - u_cell_avg_exact;

    logged_error = log10(abs(disc_error));
    logged_error(logged_error < log10(error_mask_tol)) = NaN;
    figure(figno+5)
    [X, T] = meshgrid(mesh_pa.x_centers, mesh_pa.t);
    mesh(X, T, logged_error'); view(2); view(-25, 45)
    cb = colorbar();
    cb.Limits = [log10(error_mask_tol), 0];
    
    % Make title look nice, including error norms.
    [man1, exp1] = get_scientific_decomposition(error_array_1norm(end));   e1_str   = sprintf('%.0f \\times 10^{%d}', man1, exp1);
    [man2, exp2] = get_scientific_decomposition(error_array_infnorm(end)); einf_str = sprintf('%.0f \\times 10^{%d}', man2, exp2);
    e_label_str  = '\Vert \bar{\mathbf{e}}_{\fontsize{8}{0}\selectfont\textrm{disc}} \Vert(1,\infty)';
    
    title(sprintf('$%s=(%s,%s)$', e_label_str, e1_str, einf_str))
    xlabel('$x$')
    ylabel('$t$')
    zlabel('$\log_{10} |\bar{e}_{\fontsize{8}{0}\selectfont\textrm{disc}}|$')
    axis tight
    ylim([0 max(mesh_pa.t)]) % Since function is zero at t=0, that tick label disappears 
    box on
    
    zticks([log10(error_mask_tol):2:0])
    zlim([log10(error_mask_tol), 0])
end


function u_bar = burgers_exact_sol_cell_avg(u0, x_interfaces, t)

    % Handle for computing exact solution of PDE at final time.
    u_exact_ptwise = @(x) burgers_exact_sol_pt_wise(u0, x, t(end));

    % Compute cell average of exact solution.
    u_bar = cell_average(u_exact_ptwise, x_interfaces); % Initial condition.
end

function u = burgers_exact_sol_pt_wise(u0, x, t)
    % Use MATLAB's fsolve to solve the nonlinear equation. Take as the
    % initial guess the initial condition evaluated at time zero.
    options = optimoptions('fsolve', 'Display', 'off');
    options = optimoptions(options, 'OptimalityTolerance', 1e-12);
    [u, ~, exit_flag] = fsolve(@(u) u - u0(x - u*t), u0(x), options);    
    if exit_flag~= 1
        error('fsolve did not converge...')
    end
end

% Compute exact solution using the fact that cell average at new time is
% average over spatial region at old time with lateral boundaries given by
% the characteristics. 
function u_bar = cons_lin_exact_sol_cell_avg(my_cons_law)
    h             = my_cons_law.mesh_pa.h;
    nx            = my_cons_law.mesh_pa.nx;

    tspan         = [my_cons_law.mesh_pa.t(end); 0];
    x_arrival     = my_cons_law.mesh_pa.x_interfaces;

    % Compute characteristics that bound the evolution of FV cells at final
    % time. 
    options       = odeset('RelTol', 1e-12, 'AbsTol', 1e-12);
    [~, x_depart] = ode45(@(t, x) my_cons_law.wave_speed(x, t), tspan, x_arrival, options);
    x_depart      = x_depart(end, :);
    u_bar         = zeros(nx, 1);
    
    u0 = @(x) my_cons_law.initial_condition(x);
    % Note that we have to scale the answer by 1/h since the numerical
    % solution I compare to is the cell average at the final time which is
    % the intergral divided by the mesh width.
    for i = 1:nx
        xl       = x_depart(i);
        xr       = x_depart(i+1);
        u_bar(i) = gauss_legendre_quad(u0, xl, xr, 7) / h; 
        %u_bar(i) = integral(@(x) u0(x) / h, xl, xr, 'AbsTol', 1.e-15, 'RelTol', 1.e-15);
    end    
end