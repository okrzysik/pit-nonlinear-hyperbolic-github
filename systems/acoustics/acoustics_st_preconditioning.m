% Solve space-time formulation of acoustic equations using a space-time 
% block-preconditioning strategy based on transforming to characteristic
% variables.
%
%
% Plots from Figures 1 and 2 in the paper can be generated as follows:
% --------------------------------------------------------------------
%
% Plots in Figure 1: Shows results for applying the preconditioners exactly,
% i.e., using SinT solves for the diagonal blocks.
%   prec_application = 'exact';
%       middle column: approximate_diag_blocks = ~true; 
%       right  column: approximate_diag_blocks =  true; 
% Plots in Figure 2: Shows results for applying the preconditioners 
% inexactly, i.e., using approximate PinT solves for the diagonal blocks.
%       prec_application = 'inexact';
%       approximate_diag_blocks =  true; 
%
% NOTES:
%   -Many plots in the paper combine results for both diagonal and 
%   triangular preconditioners. These can be generated by running the code 
%   for the diagonal preconditioner, leaving the figure open, and then 
%   running the code for the triangular preconditioner.



% Comments unrelated to content in the paper:
% Jan 11, 2023: Just before the main fixed-point loop, there's some code 
% that (if uncommented will solve the space-time system with GMRES using the
% characteristic-variable preconditioner as its preconditioner. Note that
% the GMRES approach does not apply any relaxation to the system, since in
% the original fixed-point iteration the relaxation essentially comes for
% free during residual computation.
% I only tested GMRES in the context that the diagonal blocks are
% inverted exactly (i.e., not with MGRIT), and it did almost nothing (maybe
% it saved one iteration here or there). This is not really surprising
% though since it's not like the inadequecies of the preconditioner live in
% some low-dimensional space that such that GMRES could easily mop them up
% (or at least this is my intuition, anyway).
%
% Jan 17, 2023: I added some code that implements better approximations to
% the Schur complement. Specifically, I included the second subdiagonal and
% third subdiagonal. This seemed to make practially no difference at all
% to the convergence rates (when using lower triangular prec)... What
% gives?



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

% Replace diagonal blocks with approximations based on Godunov advection discretizations.
approximate_diag_blocks = true; 

%prec_application = 'exact'; % Invert diag blocks with time-stepping
prec_application = 'inexact'; % Invert diag blocks with MGRIT


%% Outer solver parameters
rtol_outer    = 1e-10; % Residual reduction tol to meet
maxit_outer   = 20;  % Max number of iterations
div_tol_outer = 1e3; % Tolerance for halting 

% Use FC-based relaxation on outer iteration. 
relax_cf_OUTER     = 8; 
relax_scheme_OUTER = '';
relax_scheme_OUTER = 'F';
%relax_scheme_OUTER = 'FC';
% relax_scheme_OUTER = 'C';


%% PDE and disc parameters
nx_array = 2.^(7:11);

disc_pa.high_res = ~true; % Apply high-res corrections to disc or not.

%pde_pa.mat_param_id = 0; % c = Z = 1.
pde_pa.mat_param_id = 1; % Bale et al. example 1.
%pde_pa.mat_param_id = 2; % Bale et al. example 2.
%pde_pa.mat_param_id = 3; % Bale et al. example 3.
%pde_pa.mat_param_id = 4; % Z and c are 1, and jump to 2 and 0.5, respectively.
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


%% MGRIT parameters
cf_MGRIT = relax_cf_OUTER;
%cf_MGRIT = 8;
%maxlevels_MGRIT = 1;  % Sequential time-stepping 
%maxlevels_MGRIT = 2;  % Two-level MGRIT 
maxlevels_MGRIT = 50; % Multilevel MGRIT 

% Relaxation strategy
%relax_scheme_MGRIT = 'FCF';
relax_scheme_MGRIT = 'F-FCF'; % F-relax on level 1, FCF on all others

% Choose how coarse BE systems are solved
BE_coarse_solver = 'GMRES'; % Further parameters are set below.
%BE_coarse_solver = 'LU';
%BE_coarse_solver = 'NONE'; % This means no correction to the SL scheme will be added


% --- Package everything into a struct to pass to driver_mgrit --- %
myMGRIT_solver_params = struct();
myMGRIT_solver_params.cf = cf_MGRIT;
myMGRIT_solver_params.maxlevels = maxlevels_MGRIT;
myMGRIT_solver_params.maxiter = 1; % Ensure only a single iteration
myMGRIT_solver_params.res_halt_tol = 0;
myMGRIT_solver_params.min_coarse_nt = 2;
myMGRIT_solver_params.verbose = ~true;
myMGRIT_solver_params.verbose_assembly = true;
myMGRIT_solver_params.warm_restart = true; % Ensure same MGRIT method used over all outer iterations.

if strcmp(relax_scheme_MGRIT, 'FCF')
    myMGRIT_solver_params.pre_relax = @(level) 'FCF';
elseif strcmp(relax_scheme_MGRIT, 'F-FCF')
    myMGRIT_solver_params.pre_relax = @(level) F_then_FCF_relax(level);
end

% If outer iteration begins with an F-relaxation, there's no need to do 
% the final F-relax in MGRIT to update the error correction at F-points
% because any update to the outer iteration vector at F-points will be
% immediately overwritten on the next iteration.
if strcmp(relax_scheme_OUTER, '') 
    myMGRIT_solver_params.final_F_relax = true;
elseif strcmp(relax_scheme_OUTER(1), 'F') 
    myMGRIT_solver_params.final_F_relax = ~true;
end

myMGRIT_object = struct();



%% Norm to measure residual in.
resnorm = 2;
%resnorm = inf;


%% Plotting options
figno = 19;
save_res_fig = ~true;
res_fig_dir = './figures/paper/block-prec/';

plot_contours       = ~true;
plot_cross_sec      = ~true;
plot_material_param = true;
plot_iterative_space_time_error = ~true;

if strcmp(prec_id, 'diag'); ls = '--'; else; ls = '-'; end % DIAG prec uses dashed lines; LOWER uses solid lines
if strcmp(prec_id, 'diag'); hv = false; else; hv = true; end


%% Package parameters

% Package PDE params
pde_pa.m = 2; % Number of variables in PDE

% Solver params
solve_pa.cf           = relax_cf_OUTER;
solve_pa.relax_scheme = relax_scheme_OUTER;


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
    

    %% Setup instance of space-time system class
    % The class expects a step function in disc pa. 
    disc_pa.step = @(t0_idx, q0) step_acoustics(q0, mesh_pa, disc_pa);

    % Get acoustic system
    my_acoustic_st_system = one_step_st_system(disc_pa.step, 2*nx, nt, solve_pa);

    % Initialize RHS vector of system. This is all zeros except for the ICs.
    b = zeros(2*nx, nt);
    b(1:2*nx, 1) = q(1:2*nx, 1);
    b = b(:);
    
    % Get exact solution by time-stepping
    q_time_stepping = my_acoustic_st_system.forward_solve(b);
    %q = q_time_stepping;
        
    % Create handles for mapping back and forth between variables; these
    % are parsed to the space-time system class
    R_handle    = @(t0_idx, w0) right_eigenvector_map_acoustics(w0, mesh_pa, disc_pa);
    Rinv_handle = @(t0_idx, q0) left_eigenvector_map_acoustics(q0, mesh_pa, disc_pa);

    % Set up discretization in characteristic variables for preconditoner.
    [disc_pa.Phi_hat_11, disc_pa.Phi_hat_12, disc_pa.Phi_hat_21, disc_pa.Phi_hat_22] = ...
        Phi_hat_blocks_acoustics_godunov(mesh_pa.dt, mesh_pa.h, disc_pa.c_cc, disc_pa.Z_cc, mesh_pa.nx);
    
    % This replaces the exact diagonal blocks with inexact blocks
    % representing discretizations of advection equations that arise from 
    % ignoring changes in the impedance
    if approximate_diag_blocks
        [disc_pa.Phi_hat_11, disc_pa.Phi_hat_22] = Phi_blocks_advection_godunov(mesh_pa, disc_pa);
    else
       if disc_pa.high_res 
           error('Exact diagonal blocks not implemented for high-res discretization')
       end
    end
        

    %% Compute initial residual. 
    [q, r] = my_acoustic_st_system.relaxation(q, b); 
    resnorm_array = norm(r, resnorm);
    % Print residual norm.
    fprintf('it %d: ||r_0|| = %.2e\n', 0, resnorm_array(1))

    
    %% Re-fresh MGRIT object
    myMGRIT_object = struct();
    myMGRIT_object.t          = mesh_pa.t;
    myMGRIT_object.block_size = nx;
    
    myMGRIT_object.BE_coarse_solve.id      = BE_coarse_solver;
    myMGRIT_object.BE_coarse_solve.maxiter = 10;
    myMGRIT_object.BE_coarse_solve.rtol    = 0.01;
    
    % Create a copy for both pos and neg advection equations
    myMGRIT_object_pos = myMGRIT_object;
    myMGRIT_object_neg = myMGRIT_object;
    
    % Step function for the Godunov advection problems expects the sign of
    % c0 to be specified:
    disc_pa_neg = disc_pa; disc_pa_neg.c0_sign = -1;
    disc_pa_pos = disc_pa; disc_pa_pos.c0_sign = +1;
    
    % Make handles for advection steps that to be called by MGRIT
    step_advection_neg = @(w0, step_status, MGRIT_object) step_advection_godunov_MGRIT(w0, step_status, MGRIT_object, mesh_pa, disc_pa_neg);
    step_advection_pos = @(w0, step_status, MGRIT_object) step_advection_godunov_MGRIT(w0, step_status, MGRIT_object, mesh_pa, disc_pa_pos);
                                        
    % Create handles for solving A*e = r for some space-time matrix A.
    % These functions will output [e, rnorm, MGRIT_object]
    mgrit_solve_neg = @(e_init, r, MGRIT_object) mgrit(e_init, r, step_advection_neg, MGRIT_object, myMGRIT_solver_params);
    mgrit_solve_pos = @(e_init, r, MGRIT_object) mgrit(e_init, r, step_advection_pos, MGRIT_object, myMGRIT_solver_params);
    
    
%     %% Setup GMRES: If you want to use this, then comment out the main
%     loop below which just applied fixed-point iteration.
%     % Handle for applying matrix to a vector
%     A_handle = @(b) my_acoustic_st_system.system_action(b);
%     % Handle for applying the preconditioner (which is just the block prec, it doesn't involve aspects of the relaxation).
%     M_handle = @(b) apply_block_preconditioner_for_GMRES(b, prec_id, mesh_pa, disc_pa, Rinv_handle, R_handle, my_acoustic_st_system);
%     
%     restart = [];
%     tol     = rtol_outer;
%     maxit   = maxit_outer;
%     %x0 = b;
%     %x0 = randn(size(b)); x0(1:2*mesh_pa.nx) = b(1:2*mesh_pa.nx);
%     x0 = zeros(size(b));
%     % Solve with GMRES
%     [X,FLAG,RELRES,ITER,RESVEC] = gmres(A_handle, b, restart, tol, maxit, M_handle, [], x0);
%     % Note this is the preconditioned residual.
%     resnorm_array = RESVEC;
    
    
    %% Loop outer iterations
    for iteration_idx = 1:maxit_outer
    
        % Check if halting tolerance reached
        if resnorm_array(end)/resnorm_array(1) < rtol_outer
            break
        % Check if diverged
        elseif resnorm_array(end)/resnorm_array(1) > div_tol_outer
            solve_diverged = true;
            break
        end
        
        %% Compute characteristic residual from primitive residual
        % If no relaxation, then residual at all points is non-zero, so
        % transform the whole vector.
        if strcmp(relax_scheme_OUTER, '')
            r_hat = my_acoustic_st_system.all_point_transformation(r, Rinv_handle);
        
        % If ended on an F-relax, then only C-point residuals are non-zero,
        % so just transform these.
        elseif strcmp(relax_scheme_OUTER(end), 'F')
            r_hat = my_acoustic_st_system.C_point_transformation(r, Rinv_handle);
            
        % If ended on a C-relax, then only F-point residuals are non-zero,
        % so just transform these.
        elseif strcmp(relax_scheme_OUTER(end), 'C')
            r_hat = my_acoustic_st_system.F_point_transformation(r, Rinv_handle);

        end
        
        e_hat = r_hat;
        
        %% Solve for characteristic error.
        if strcmp(prec_application, 'exact')
            e_hat = block_preconditioner_EXACT(prec_id, mesh_pa, disc_pa, r_hat);
            
        elseif strcmp(prec_application, 'inexact')
            assert(approximate_diag_blocks, 'INEXACT prec via MGRIT requires the use of approximate diagonal blocks')
            [e_hat, myMGRIT_object_neg, myMGRIT_object_pos] = ...
                block_preconditioner_INEXACT(prec_id, mesh_pa, disc_pa, r_hat, ...
                                                mgrit_solve_neg, mgrit_solve_pos, ...
                                                myMGRIT_object_neg, myMGRIT_object_pos);
        end

        %% Transform characteristic error to primitive error.
        % If no relaxation, then error at all points is non-zero, so
        % transform the whole vector.
        if strcmp(relax_scheme_OUTER, '')
            e = my_acoustic_st_system.all_point_transformation(e_hat, R_handle);
        
        % If outer relaxation begins with an F-relaxation, then there's no
        % point transforming preconditioned error at F-points because the
        % correction will be immediately overwritten
        elseif strcmp(relax_scheme_OUTER(1), 'F')
            e = my_acoustic_st_system.C_point_transformation(e_hat, R_handle);
            
        % If outer relaxation begins with an C-relaxation, then there's no
        % point transforming preconditioned error at C-points because the
        % correction will be immediately overwritten
        elseif strcmp(relax_scheme_OUTER(1), 'C')
            e = my_acoustic_st_system.F_point_transformation(e_hat, R_handle);
        end
        
        %exact = my_acoustic_st_system.forward_solve(r);
                
       
        %% Apply error correction
        q = q + e;
        
        %% Relax and compute new residual.
        [q, r] = my_acoustic_st_system.relaxation(q, b);

        % Plot the current algebraic error in space-time.
        if plot_iterative_space_time_error
            fh = plot_space_time_error(figno, q, q_time_stepping, mesh_pa, iteration_idx); figno = figno+1;
            pause
        end
        
        % Compute error and residual norms.
        resnorm_array = [resnorm_array; norm(r, resnorm)];

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

prec_description = '';
prec_label_title = '';
if approximate_diag_blocks
    prec_label_title = '$\widetilde{\mathcal{P}}$'; % Inexact diag blocks
    prec_description = 'approx-diag';
else
    prec_label_title = '$\widehat{\mathcal{P}}$'; % Exact diag blocks
    prec_description = 'true-diag';
end

if strcmp(prec_application, 'exact')
    prec_label_title = strcat(prec_label_title, ' (SinT)');
    prec_description = strcat(prec_description, '-SinT');
elseif strcmp(prec_application, 'inexact')
    prec_label_title = strcat(prec_label_title, ' (PinT)');
    prec_description = strcat(prec_description, '-PinT');
end

title(sprintf('residual history: %s', prec_label_title))
axis tight
ylim([1e-10, 1])
xlim([0 maxit_outer])
%
yticks((10).^(-10:2:0))
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
    res_fig_name = sprintf('acoustics-block-prec-ex%d-m%d-T%.1f-P%s', ...
        pde_pa.mat_param_id, relax_cf_OUTER, mesh_pa.tmax, prec_description);
    
    figure(res_fig);
    figure_saver(gcf, sprintf('%s/%s', res_fig_dir, res_fig_name), true); 
end
figno = figno+1;

%% Material parameters plot
if plot_material_param
    mat_pa_fig = acoustics_plot_material_parameters(figno, mesh_pa, disc_pa); figno = figno+1;
    
    if save_res_fig
        mat_pa_fig_name = sprintf('acoustics-block-prec-ex%d-mat-pa', ...
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

        
% Solve P*e = r where P is a block 2x2 preconditioner.
% Diagonal blocks of P are approximarely inverted with MGRIT.
function [e, myMGRIT_object_neg, myMGRIT_object_pos] = ...
                block_preconditioner_INEXACT(prec_id, mesh_pa, disc_pa, r, ...
                                                mgrit_solve_neg, mgrit_solve_pos, ...
                                                myMGRIT_object_neg, myMGRIT_object_pos)
    nx = mesh_pa.nx;
    nt = mesh_pa.nt;

    r  = reshape(r, [2*nx, nt]);
    r1 = r(1:nx, :);      r1 = r1(:);
    r2 = r(nx+1:2*nx, :); r2 = r2(:);
    
    % Block diagonal preconditioner
    if strcmp(prec_id, 'diag')

        % Invert 11 block
        [e1, ~, myMGRIT_object_neg] = mgrit_solve_neg(r1, r1, myMGRIT_object_neg);        
        % Invert 22 block
        [e2, ~, myMGRIT_object_pos] = mgrit_solve_pos(r2, r2, myMGRIT_object_pos);

    % Block lower triangular preconditioner
    elseif strcmp(prec_id, 'lower')

        % Invert 11 block
        [e1, ~, myMGRIT_object_neg] = mgrit_solve_neg(r1, r1, myMGRIT_object_neg);        
        
        % Update RHS of first variable by adding solution from second
        % variable mutliplied by the (2,1) block
        r2 = r2 - A_hat_kj_action(e1, disc_pa.Phi_hat_21, mesh_pa);
        
        % Invert 22 block
        [e2, ~, myMGRIT_object_pos] = mgrit_solve_pos(r2, r2, myMGRIT_object_pos);
        
    % Block upper triangular preconditioner
    elseif strcmp(prec_id, 'upper')
        
        % Invert 22 block
        [e2, ~, myMGRIT_object_pos] = mgrit_solve_pos(r2, r2, myMGRIT_object_pos);
        
        % Update RHS of first variable by adding solution from second
        % variable mutliplied by the (1,2) block
        r1 = r1 - A_hat_kj_action(e2, disc_pa.Phi_hat_12, mesh_pa);
        
        % Invert 11 block
        [e1, ~, myMGRIT_object_neg] = mgrit_solve_neg(r1, r1, myMGRIT_object_neg);  
        
    else
        error('prec_id = %s not recognised', prec_id)
    end
    
    e = [reshape(e1, [nx, nt]); reshape(e2, [nx, nt])];
    e = e(:);
        
end
        
        
% Solve P*e = r where P is a block 2x2 preconditioner.
% Diagonal blocks of P are inverted exactly.
function e = block_preconditioner_EXACT(prec_id, mesh_pa, disc_pa, r)
    nx = mesh_pa.nx;
    nt = mesh_pa.nt;

    r  = reshape(r, [2*nx, nt]);
    r1 = r(1:nx, :);      r1 = r1(:);
    r2 = r(nx+1:2*nx, :); r2 = r2(:);
    
    % Block diagonal preconditioner
    if strcmp(prec_id, 'diag')

        % Invert 11 block
        e1 = A_hat_kk_time_stepping_solve(r1, disc_pa.Phi_hat_11, mesh_pa);
        % Invert 22 block
        e2 = A_hat_kk_time_stepping_solve(r2, disc_pa.Phi_hat_22, mesh_pa);

    % Block lower triangular preconditioner
    elseif strcmp(prec_id, 'lower')

        % Invert 11 block
        e1 = A_hat_kk_time_stepping_solve(r1, disc_pa.Phi_hat_11, mesh_pa);
        
        % Update RHS of first variable by adding solution from second
        % variable mutliplied by the (2,1) block
        r2 = r2 - A_hat_kj_action(e1, disc_pa.Phi_hat_21, mesh_pa);
        
        % Invert 22 block
        e2 = A_hat_kk_time_stepping_solve(r2, disc_pa.Phi_hat_22, mesh_pa);
        
        %e2 = schur_comp_approx_2_solve(r2, disc_pa.Phi_hat_22, disc_pa.Phi_hat_12, disc_pa.Phi_hat_21, mesh_pa);
        %e2 = schur_comp_approx_3_solve(r2, disc_pa.Phi_hat_22, disc_pa.Phi_hat_11, disc_pa.Phi_hat_12, disc_pa.Phi_hat_21, mesh_pa);
        
    % Block upper triangular preconditioner
    elseif strcmp(prec_id, 'upper')
        
        % Invert 22 block
        e2 = A_hat_kk_time_stepping_solve(r2, disc_pa.Phi_hat_22, mesh_pa);
        
        % Update RHS of first variable by adding solution from second
        % variable mutliplied by the (1,2) block
        r1 = r1 - A_hat_kj_action(e2, disc_pa.Phi_hat_12, mesh_pa);
        
        % Invert 11 block
        e1 = A_hat_kk_time_stepping_solve(r1, disc_pa.Phi_hat_11, mesh_pa);
        
    else
        error('prec_id = %s not recognised', prec_id)
    end
    
    e = [reshape(e1, [nx, nt]); reshape(e2, [nx, nt])];
    e = e(:);
        
end

% Compute the action of an off-diagonal block of A_hat on the vector y.
% Both off-diagonal blocks are all zero matrices except for a non-zero
% subdiagonal which has the block parsed to this function.
% So, to compute this action, x is zero in the first block, and then
% all other blocks are lagged one behind y, and scaled by -Phi_hat_kj
function x = A_hat_kj_action(y, Phi_hat_kj, mesh_pa)
    y = reshape(y, [mesh_pa.nx, mesh_pa.nt]);
    
    x = zeros(size(y));
    for n = 2:mesh_pa.nt
        x(:, n) = -(Phi_hat_kj * y(:, n-1));
    end
    
    x = x(:);
end

% Invert a diagonal block of A_hat to solve the system A_hat * x = y.
% The diagonal block is a lower bi-diagonal matrix with identities on the
% diagonal and -Phi_hat_kk on the subdiagonal.
function x = A_hat_kk_time_stepping_solve(y, Phi_hat_kk, mesh_pa)

    y = reshape(y, [mesh_pa.nx, mesh_pa.nt]);
    
    x = zeros(size(y));
    x(:, 1) = y(:, 1);
    for n = 1:mesh_pa.nt-1
        x(:, n+1) = Phi_hat_kk * x(:, n) + y(:, n+1);
    end
    
    x = x(:);
end

% Do a solve of the 22 block when it is not just A_22, but it also includes 
% the second subdiagonal from the true Schur complement.
function x = schur_comp_approx_2_solve(y, Phi_hat_22, Phi_hat_12, Phi_hat_21, mesh_pa)

    % gamma = 1 is the actualy value this should be, but I was wondering if
    % changing it to some other constant might work better. In any event,
    % it doesn't seem to do much...
    gamma = 1; 

    y = reshape(y, [mesh_pa.nx, mesh_pa.nt]);
    
    x = zeros(size(y));
    x(:, 1) = y(:, 1);
    x(:, 2) = Phi_hat_22 * x(:, 1) + y(:, 2);
    for n = 2:mesh_pa.nt-1
        x(:, n+1) = Phi_hat_22 * x(:, n) ...
                    + gamma * (Phi_hat_21 * (Phi_hat_12 * x(:, n-1))) ...
                    + y(:, n+1);
    end
    
    x = x(:);
end

% Do a solve of the 22 block when it is not just A_22, but it also includes 
% the second and third subdiagonal from the true Schur complement.
function x = schur_comp_approx_3_solve(y, Phi_hat_22, Phi_hat_11, Phi_hat_12, Phi_hat_21, mesh_pa)

    y = reshape(y, [mesh_pa.nx, mesh_pa.nt]);
    
    x = zeros(size(y));
    x(:, 1) = y(:, 1);
    x(:, 2) = Phi_hat_22 * x(:, 1) + y(:, 2);
    x(:, 3) = Phi_hat_22 * x(:, 2) + y(:, 3) + Phi_hat_21 * (Phi_hat_12 * x(:, 1));
    for n = 3:mesh_pa.nt-1
        x(:, n+1) = Phi_hat_22 * x(:, n) ...
                    + Phi_hat_21 *                (Phi_hat_12 * x(:, n-1)) ...
                    + Phi_hat_21 * ( Phi_hat_11 * (Phi_hat_12 * x(:, n-2))) ...
                    + y(:, n+1);
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


% Plot error and residual in space time
function fh = plot_space_time_error(figno, q, q_time_stepping, mesh_pa, iteration_idx)
    fh = figure(figno);
    nx = mesh_pa.nx;
    nt = mesh_pa.nt;
    [X, T] = meshgrid(mesh_pa.x_centers, mesh_pa.t);
    error_algebraic = (reshape(q, [2*nx, nt]) - reshape(q_time_stepping, [2*nx, nt]));
    %subplot(1, 2, 1)
    mesh(X, T, log10(abs(error_algebraic(1:nx, :)))')
    %subplot(1, 2, 2)
    %r_reshaped = reshape(r, [2*nx nt]);
    %mesh(X, T, log10(abs(r_reshaped(1:nx, :)))')
    title(sprintf('$n_x = %d$. it=%d. error $= p - p_k$', nx, iteration_idx))
    xlabel('$x$')
    ylabel('$t$')
    zlabel('$\log_{10}| p - p_k |$')
end


% Implement level-dependent F-FCF MGRIT relaxation
function nu = F_then_FCF_relax(level)
    if level == 1 
        nu = 'F';
    else
        nu = 'FCF';
    end
end

% Apply block preconditioner to the vector r. I.e., transform to char var,
% apply the 2x2 prec in char var and map back.
function e = apply_block_preconditioner_for_GMRES(r, prec_id, mesh_pa, disc_pa, Rinv_handle, R_handle, my_acoustic_st_system)

    r_hat = my_acoustic_st_system.all_point_transformation(r, Rinv_handle);
    e_hat = block_preconditioner_EXACT(prec_id, mesh_pa, disc_pa, r_hat);
    e     = my_acoustic_st_system.all_point_transformation(e_hat, R_handle);

end



            
%             if n > 2
%                 e2(:, n+1) = e2(:, n+1) + disc_pa.LUMP .* e2(:, n-1);
%             end

    
           % trying to do something with the largest term in the Schur complement that's ignored...     
%     Z = disc_pa.Z_cc;
%     Z_left  = [Z(end); Z(1:end-1)];
%     Z_right = [Z(2:end); Z(1)];         
%     LUMP = (mesh_pa.dt / mesh_pa.h * disc_pa.c_cc).^2 .* (Z - Z_right)./(Z + Z_right) .* (Z_left - Z)./(Z_left + Z);
%     disc_pa.LUMP = LUMP;
        

%         [Phi_pp, Phi_pu, Phi_up, Phi_uu] = Godunov_blocks(mesh_pa.dt, mesh_pa.h, disc_pa.c_cc, disc_pa.Z_cc, mesh_pa.nx, true);
%         Phi      = [Phi_pp, Phi_pu; Phi_up, Phi_uu];
%         Z_mat     = diag(disc_pa.Z_cc);
%         Z_mat_inv = diag(1./disc_pa.Z_cc);
%         I    = speye(mesh_pa.nx);
%         R    =     [-Z_mat, Z_mat; I, I];
%         Rinv = 0.5*[-Z_mat_inv, I; Z_mat_inv, I];
%         Phi_11_char_true = 0.5*( (Phi_uu + Z_mat_inv * Phi_pp * Z_mat) - (Phi_up*Z_mat + Z_mat_inv*Phi_pu) );
%         Phi_22_char_true = 0.5*( (Phi_uu + Z_mat_inv * Phi_pp * Z_mat) + (Phi_up*Z_mat + Z_mat_inv*Phi_pu) );
% 
%         disc_pa.Phi_hat_11 = Phi_11_char_true;
%         disc_pa.Phi_hat_22 = Phi_22_char_true;
    
    
%     % For the periodic medium, rather than using Z-linearized diagonal
%     % blocks, we could use homogenized blocks. When I tested this it didn't
%     % seem to work in the slightest :( It seems like a nice idea though... 
%     if homogenized_diag_blocks
%         temp = pde_pa.mat_param_id_homogenize;
%         pde_pa.mat_param_id_homogenize = true;
%         [disc_pa.c_cc, disc_pa.Z_cc] = material_parameters(mesh_pa.x_centers, pde_pa);
%         [disc_pa.Phi_hat_11, disc_pa.Phi_hat_22] = godunov_advection_blocks(mesh_pa, disc_pa);
%         
%         % Restore original material parameters, whatever they were.
%         pde_pa.mat_param_id_homogenize = temp;
%         [disc_pa.c_cc, disc_pa.Z_cc] = material_parameters(mesh_pa.x_centers, pde_pa);
%     end