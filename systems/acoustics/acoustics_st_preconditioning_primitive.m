% Solve space-time formulation of acoustic equations using a space-time 
% block-preconditioning strategy to the primitive variables of the system.
%
% There are not results for this in the main paper, but they are in the
% Supplementary Materials section, shown in Figure SM1.
%
% This file is a modification of "acoustics_st_preconditioning.m," which 
% applies block preconditioning on the space-time system after converting 
% to characteristic variables.
%
% The results produced by this script show that applying block
% preconditioning in primitive variables does not work.
%
% To generate Figure SM1:
%   Left:  Use pde_pa.mat_param_id = 1
%   Right: Use pde_pa.mat_param_id = 2



tmp = matlab.desktop.editor.getActive;
cd(fileparts(tmp.Filename));


%close all
clear
clc
rng(1)


%% 2x2 block preconditioning parameters
%prec_id = 'upper';
prec_id = 'lower';
%prec_id = 'diag';


%% Outer solver parameters
rtol_outer    = 1e-10; % Residual reduction tol to meet
maxit_outer   = 20;  % Max number of iterations
div_tol_outer = 1e40; % Tolerance for halting 


%% PDE and disc parameters
nx_array = 2.^(7:10);

disc_pa.high_res = ~true; % Apply high-res corrections to disc or not.

%pde_pa.mat_param_id = 0; % c = Z = 1.
pde_pa.mat_param_id = 1; % Bale et al. example 1.
%pde_pa.mat_param_id = 2; % Bale et al. example 2.
%pde_pa.mat_param_id = 3; % Bale et al. example 3.
% % pde_pa.mat_param_id = 4; % Z and c are 1, and jump to 2 and 0.5, respectively.
%pde_pa.mat_param_id = 5; % Periodically layered medium
%pde_pa.mat_param_id = 6; % Randomly layered medium

% Extra parameters for layered media:
pde_pa.mat_param_num_layers = 16; 
%pde_pa.mat_param_homogenize = ~true;

mesh_pa.tmax       = 1;
disc_pa.CFL_number = 0.85;

% Spatial domain is set to [0,1] in Bale et al.
mesh_pa.xmin = 0;
mesh_pa.xmax = 1;


%% Norm to measure residual in.
resnorm = 2;

%% Plotting options
figno = 19;
save_res_fig = ~true;
res_fig_dir = './figures/paper/block-prec-prim/';

plot_contours       = ~true;
plot_cross_sec      = ~true;
plot_material_param = true;
plot_iterative_space_time_error = true;

if strcmp(prec_id, 'diag'); ls = '--'; else; ls = '-'; end % DIAG prec uses dashed lines; LOWER uses solid lines
if strcmp(prec_id, 'diag'); hv = false; else; hv = true; end


%% Package parameters

% Package PDE params
pde_pa.m = 2; % Number of variables in PDE

% Solver params
solve_pa = '';


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% --- Solve problem at different spatial resolutions --- %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for nx_idx = 1:numel(nx_array)
    
    mesh_pa.nx = nx_array(nx_idx);
    
    %% Set up discretization
    % Create spatial mesh
    [mesh_pa.x_interfaces, mesh_pa.x_centers, mesh_pa.h] = space_mesh(mesh_pa.xmin, mesh_pa.xmax, mesh_pa.nx);
    
    % Get material parameters. Spatial disc expects these stored like this.
    [disc_pa.c_cc, disc_pa.Z_cc] = acoustics_material_parameters(mesh_pa.x_centers, pde_pa);
    pde_pa.abs_fprime_max = max(disc_pa.c_cc);
    
    % Get temporal mesh.
    [mesh_pa.t, mesh_pa.dt] = time_mesh(mesh_pa.h, disc_pa.CFL_number, pde_pa.abs_fprime_max, 1, mesh_pa.tmax);
    mesh_pa.nt = numel(mesh_pa.t);

    nx = mesh_pa.nx;
    nt = mesh_pa.nt;
    fprintf('\nnx=%d, nt=%d\n', nx, nt);

    % PDE initial condition.
    p0 = initial_condition_p(mesh_pa.x_centers); 
    u0 = initial_condition_u(mesh_pa.x_centers); 

    % Initialize solution vector. This is all random except for the ICS.
    q = randn(2*nx, nt);
    q(1:2*nx, 1) = [p0; u0];
    q = q(:);
    
    % Initialize RHS vector of system. This is all zeros except for the ICs.
    b = zeros(2*nx, nt);
    b(1:2*nx, 1) = q(1:2*nx, 1);
    b = b(:);

    % Get blocks of 2x2 time-stepping operator
    [disc_pa.Phi_pp, disc_pa.Phi_pu, disc_pa.Phi_up, disc_pa.Phi_uu] = Phi_blocks_acoustics_godunov(mesh_pa.dt, mesh_pa.h, disc_pa.c_cc, disc_pa.Z_cc, mesh_pa.nx, false);
        

    % Get acoustic system
    disc_pa.step = @(t0_idx, q0) step_acoustics(q0, mesh_pa, disc_pa);
    my_acoustic_st_system = one_step_st_system(disc_pa.step, 2*nx, nt, solve_pa);
    % Get exact solution by time-stepping
    q_time_stepping = my_acoustic_st_system.forward_solve(b);
    %q = q_time_stepping;

    %% Compute initial residual. 
    r = space_time_residual(q, b, mesh_pa, disc_pa);
    resnorm_array = norm(r, resnorm);
    % Print residual norm.
    fprintf('it %d: ||r_0|| = %.2e\n', 0, resnorm_array(1))

    
    %% Iterate
    for iteration_idx = 1:maxit_outer
    
        % Check if halting tolerance reached
        if resnorm_array(end)/resnorm_array(1) < rtol_outer
            break
        % Check if diverged
        elseif resnorm_array(end)/resnorm_array(1) > div_tol_outer
            solve_diverged = true;
            break
        end
        
        % Apply preconditioner to get approximate error
        e = block_preconditioner(r, prec_id, mesh_pa, disc_pa);
        
        % Apply error correction
        q = q + e;
        
        % Compute residual norm.
        r = space_time_residual(q, b, mesh_pa, disc_pa);
        resnorm_array = [resnorm_array; norm(r, resnorm)];

        % % Just a sanity check that we're computing the residual correctly
        % % in "space_time_residual"
        % norm(my_acoustic_st_system.residual(q, b) - space_time_residual(q, b, mesh_pa, disc_pa))/norm(q)

        % Print residual norm.
        fprintf('it %d: ||r_%d||/||r_0|| = %.2e, ||r_%d||/||r_%d|| = %.2f\n', ...
            iteration_idx, iteration_idx,   resnorm_array(end)/resnorm_array(1), ...
            iteration_idx, iteration_idx-1, resnorm_array(end)/resnorm_array(end-1))
        
    end
    % End of outer iterations

    %% Plot residuals from this solve
    nameval = lineprops(nx_idx);

    %% Plots: Richardson convergence
    % Residuals
    res_fig = figure(figno);
    % Show handle for nx x nt or not. 
    if hv
        %displayname = sprintf('$(n_x,n_t)=(%d,%d)$', nx, nt); 
        displayname = sprintf('$%d \\times %d$', nx, nt); 
        display = {'DisplayName', displayname};
    else
        %displayname = sprintf('$(n_x,n_t)=(%d,%d)$', nx, nt); 
        displayname = sprintf('$%d \\times %d$', nx, nt); 
        display = {'HandleVisibility', 'Off'};
    end
    semilogy(0:numel(resnorm_array)-1, resnorm_array/resnorm_array(1), nameval{:}, ...
        'LineWidth', 2, ...
        'LineStyle', ls, ...
        display{:})
    hold on
        
    
end
% End of looping over different spatial resolutions.

%%
%%%%%%%%%%%%%%%%%
% --- Plots --- %
%%%%%%%%%%%%%%%%%
%% Residual plot
% Pretty up the plot
figure(res_fig); 


%title(sprintf('residual history: %s', prec_label_title))
axis tight
%ylim([1e-10, 1])
%xlim([0 maxit_outer])
%
%yticks((10).^(-10:2:0))
grid minor
grid on
grid minor
set(gca,'YMinorTick','Off')
%grid minor

lh = legend();
lh.set('Location', 'best')
xlabel('$k$')
ylabel('$\Vert \mathbf{r}_k \Vert / \Vert \mathbf{r}_0 \Vert$')

% Save plot of residual history
if save_res_fig
    res_fig_name = sprintf('acoustics-block-prec-prim-ex%d-T%.1f-P%s', ...
        pde_pa.mat_param_id, mesh_pa.tmax, prec_id);
    
    figure(res_fig);
    figure_saver(gcf, sprintf('%s/%s', res_fig_dir, res_fig_name), true); 
end
figno = figno+1;

%% Material parameters plot
if plot_material_param
    mat_pa_fig = acoustics_plot_material_parameters(figno, mesh_pa, disc_pa); figno = figno+1;
    
    if save_res_fig
        mat_pa_fig_name = sprintf('acoustics-block-prec-prim-ex%d-mat-pa', ...
        pde_pa.mat_param_id);
        figure(mat_pa_fig);
        figure_saver(gcf, sprintf('%s/%s', res_fig_dir, mat_pa_fig_name), ~true);
    end
end


%% Solution plot
[mesh_pa.X, mesh_pa.T] = meshgrid(mesh_pa.x_centers, mesh_pa.t);
pa.disc_pa = disc_pa;
pa.mesh_pa = mesh_pa;
pa.ncon_lvl = 10;
pa.con_tol  = 1e-3;
pa.invisible_col_bar = ~true;

q = reshape(q, [2*mesh_pa.nx, mesh_pa.nt]);
P = q(1:mesh_pa.nx, :);
U = q(mesh_pa.nx+1:2*mesh_pa.nx, :);

mesh_nx_tol = 2*1024;
if mesh_pa.nx > mesh_nx_tol
    pa.meshing_down_sample_factor = mesh_pa.nx / mesh_nx_tol;
end

if plot_contours
    p_con_fig = contour_plot(figno, P, '$p(x, t)$', pa); figno = figno+1;
    u_con_fig = contour_plot(figno, U, '$u(x, t)$', pa); figno = figno+1;
end

if plot_cross_sec
    p_cs_fig = cross_sec_plot(figno, P, '$p(x, t)$', pa); figno = figno+1;
    u_cs_fig = cross_sec_plot(figno, U, '$u(x, t)$', pa); figno = figno+1;
end
        
 
% Compute residual of 2x2 space-time system. q and b are in the standard
% ordering where p and u are blocked together at each time point. To compute 
% this residual we unpack them into p and u at all time points and evaluate the
% residual in block 2x2 form. We then repackage them into a vector with the
% same ordering as q.
function r = space_time_residual(q, b, mesh_pa, disc_pa)
    nx = mesh_pa.nx;
    nt = mesh_pa.nt;

    % Unpack b into p and u components at all time points
    b   = reshape(b, [2*nx, nt]);
    b_p = b(1:nx, :);      b_p = b_p(:);
    b_u = b(nx+1:2*nx, :); b_u = b_u(:);

    % Unpack q into p and u components at all time points
    q = reshape(q, [2*nx, nt]);
    q_p = q(1:nx, :);      q_p = q_p(:);
    q_u = q(nx+1:2*nx, :); q_u = q_u(:);

    % 
    r_p = b_p - (A_kk_action(q_p, disc_pa.Phi_pp, mesh_pa) + A_kj_action(q_u, disc_pa.Phi_pu, mesh_pa));
    r_u = b_u - (A_kj_action(q_p, disc_pa.Phi_up, mesh_pa) + A_kk_action(q_u, disc_pa.Phi_uu, mesh_pa));

    % Package r so it has the same ordering as q.
    r = [reshape(r_p, [nx, nt]); reshape(r_u, [nx, nt])];
    r = r(:);

end

% Solve P*e = r where P is a block 2x2 preconditioner based on primitive
% variables.
% r is ordered such that p and u and blocked together at each time point.
% Diagonal blocks of P are inverted exactly by time-stepping
function e = block_preconditioner(r, prec_id, mesh_pa, disc_pa)
    nx = mesh_pa.nx;
    nt = mesh_pa.nt;

    r  = reshape(r, [2*nx, nt]);
    r_p = r(1:nx, :);      r_p = r_p(:); % Residual associated with p
    r_u = r(nx+1:2*nx, :); r_u = r_u(:); % Residual associated with u
    
    % Block diagonal preconditioner
    if strcmp(prec_id, 'diag')

        % Invert 11 block
        e_p = A_kk_time_stepping_solve(r_p, disc_pa.Phi_pp, mesh_pa);
        % Invert 22 block
        e_u = A_kk_time_stepping_solve(r_u, disc_pa.Phi_uu, mesh_pa);

    % Block lower triangular preconditioner
    elseif strcmp(prec_id, 'lower')

        % Invert 11 block
        e_p = A_kk_time_stepping_solve(r_p, disc_pa.Phi_pp, mesh_pa);
        
        % Update RHS of first variable by adding solution from second
        % variable mutliplied by the (2,1) block
        r_u = r_u - A_kj_action(e_p, disc_pa.Phi_up, mesh_pa);
        
        % Invert 22 block
        e_u = A_kk_time_stepping_solve(r_u, disc_pa.Phi_uu, mesh_pa);
        
        
    % Block upper triangular preconditioner
    elseif strcmp(prec_id, 'upper')
        
        % Invert 22 block
        e_u = A_kk_time_stepping_solve(r_u, disc_pa.Phi_hat_uu, mesh_pa);
        
        % Update RHS of first variable by adding solution from second
        % variable mutliplied by the (1,2) block
        r_p = r_p - A_kj_action(e_u, disc_pa.Phi_pu, mesh_pa);
        
        % Invert 11 block
        e_p = A_kk_time_stepping_solve(r_p, disc_pa.Phi_pp, mesh_pa);
        
    else
        error('prec_id = %s not recognised', prec_id)
    end
    
    e = [reshape(e_p, [nx, nt]); reshape(e_u, [nx, nt])];
    e = e(:);
        
end

% Compute the action of an diagonal block of A on the vector y.
% Both off-diagonal blocks are block lower bi-diagonal matrices with
% identities on the diagonal and -Phi_kk blocks on the subdiagonal.
function x = A_kk_action(y, Phi_kk, mesh_pa)
    y = reshape(y, [mesh_pa.nx, mesh_pa.nt]);
    
    x = zeros(size(y));

    x(:, 1) = y(:, 1); % First block row
    % All later block rows
    for n = 2:mesh_pa.nt
        x(:, n) = y(:, n) - (Phi_kk * y(:, n-1));
    end
    
    x = x(:);
end

% Compute the action of an off-diagonal block of A on the vector y.
% Both off-diagonal blocks are all zero matrices except for a non-zero
% subdiagonal which has the block parsed to this function.
% So, to compute this action, x is zero in the first block, and then
% all other blocks are lagged one behind y, and scaled by -Phi_kj
function x = A_kj_action(y, Phi_kj, mesh_pa)
    y = reshape(y, [mesh_pa.nx, mesh_pa.nt]);
    
    x = zeros(size(y));
    for n = 2:mesh_pa.nt
        x(:, n) = -(Phi_kj * y(:, n-1));
    end
    
    x = x(:);
end

% Invert a diagonal block of A to solve the system A * x = y.
% The diagonal block is a lower bi-diagonal matrix with identities on the
% diagonal and -Phi_kk on the subdiagonal.
function x = A_kk_time_stepping_solve(y, Phi_hat_kk, mesh_pa)

    y = reshape(y, [mesh_pa.nx, mesh_pa.nt]);
    
    x = zeros(size(y));
    x(:, 1) = y(:, 1);
    for n = 1:mesh_pa.nt-1
        x(:, n+1) = Phi_hat_kk * x(:, n) + y(:, n+1);
    end
    
    x = x(:);
end


% See Bale et al. eq. (5.5)
function p0 = initial_condition_p(x)
    p0    = ones(size(x));
    I     = find(x > 0.4 & x < 0.6);
    p0(I) = 7/4 - 3/4*cos( 10*pi*x(I) - 4*pi );
end

function u0 = initial_condition_u(x)
    u0    = zeros(size(x));
end

function q1 = step_acoustics(q0, mesh_pa, disc_pa)
    
    q1 = q0 + mesh_pa.dt * acoustics_spatial_discretization_godunov(q0, mesh_pa, disc_pa);
    
    if norm(q1, inf) > 10*norm(q0, inf)
        error('step: q1 > 10*q0')
    end
end
