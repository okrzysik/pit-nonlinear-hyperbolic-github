% This is experimenting numerically to understand what the applicability is
% of the condition number bound I derived (it's not very much).
%
% Basically it's applicable to very weakly varying impedance problems on short
% time intervals.

tmp = matlab.desktop.editor.getActive;
cd(fileparts(tmp.Filename));


%close all
clear
clc
rng(1)


%% PDE and disc parameters
nx_array = 2.^(5:6);

%disc_pa.high_res = true; % Apply high-res corrections to disc or not.

pde_pa.mat_param_id = 0; % c = Z = 1.
pde_pa.mat_param_id = 1; % Bale et al. example 1.
pde_pa.mat_param_id = 2; % Bale et al. example 2.
pde_pa.mat_param_id = 3; % Bale et al. example 3.
% % % pde_pa.mat_param_id = 4; % Z and c are 1, and jump to 2 and 0.5, respectively.
%pde_pa.mat_param_id = 5; % Periodically layered medium
%pde_pa.mat_param_id = 6; % Randomly layered medium

% Extra parameters for layered media:
pde_pa.mat_param_num_layers = 16; 
%pde_pa.mat_param_homogenize = ~true;

mesh_pa.tmax       = 20;

disc_pa.CFL_number = 0.85;

% Spatial domain is set to [0,1] in Bale et al.
mesh_pa.xmin = 0;
mesh_pa.xmax = 1;




%% Package parameters

solve_pa = [];

% Package PDE params
pde_pa.m = 2; % Number of variables in PDE

sigma_11 = [];
sigma_12 = [];
sigma_21 = [];
sigma_22 = [];

Phi_hat_11_powers = [];

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% --- Solve problem at different spatial resolutions --- %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
temp = [];
temp2 = [];
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

    
    % Set up discretization in characteristic variables for preconditoner.
    [Phi_hat_11, Phi_hat_12, Phi_hat_21, Phi_hat_22] = ...
        Phi_hat_blocks_acoustics_godunov(mesh_pa.dt, mesh_pa.h, disc_pa.c_cc, disc_pa.Z_cc, mesh_pa.nx);
    Phi = Phi_hat_11;
    
    
    T = @(theta) speye(nx) - exp(1i*theta)*Phi;
    
    N = 32;
    theta = linspace(0, pi, N); % symmetry means we only need to consider [0, pi]
    for n = 1:N   
        
        %Phi_hat_11_phi = [Phi_hat_11_phi; svds( T, 1, 'smallest' )];
        
        T_min_sing(n) = sqrt(eigs( full(T(theta(n))*T(theta(n))'), 1, 'smallestreal' ));
        
        
    end
    
    temp = [temp; 1/sqrt(min(T_min_sing))];
    
    figure(2) 
    plot(theta, T_min_sing, 'o-', 'LineWidth', 2, 'DisplayName', sprintf('$n_t = %d$', nt))
    hold on
    xlabel('$\theta$')
    lh = legend();


%     
L = Phi - speye(nx);
%     
%     F = @(theta) -Phi*exp(-1i*theta) + (Phi * Phi' + speye(nx)) - Phi'*exp(1i*theta);
%     %F = @(theta) -Phi*exp(-1i*theta) - Phi'*exp(1i*theta);
%     
%     %F = @(theta) + (Phi * Phi' + speye(nx));
%     
%     %F = @(theta) (1-exp(-1i*theta))*L + (1-exp(+1i*theta))*L';
%     
%     N = 32;
%     theta = linspace(0, pi, N); % symmetry means we only need to consider [0, pi]
%     F_eigs = [];
%     for n = 1:N  
%         F_eigs(n) = eigs( full(F(theta(n))), 1, 'smallestreal' );
%         
%     end
%     temp = [temp; 1/sqrt(min(F_eigs))];
%     
%     %
%     %temp2 = [temp2; 1/sqrt( eigs( full( L * L' + 4*speye(nx) + 2*(L + L') ), 1, 'smallestabs' ) )];
    temp2 = [temp2; 1/sqrt( eigs( full( L * L' ), 1, 'smallestabs' ) )];
    
%     figure(2) 
%     plot(theta, F_eigs, 'o-', 'LineWidth', 2, 'DisplayName', sprintf('$n_t = %d$', nt))
%     hold on
%     %plot(theta, F_eigs + 2*(1 - cos(theta)), '>-', 'LineWidth', 2, 'DisplayName', sprintf('$n_t = %d$', nt))
%     xlabel('$\theta$')
%     F_eigs(1)
        
    
    
%     % Looking at norms of powers of the blocks. The sum of these things can
%     % be used to bound the 2-norm of the inverse space-time matrix. 
%     Phi_hat_11_powers = [];
%     N = 600;
%     for n = 1:N   
%        Phi_hat_11_powers = [Phi_hat_11_powers; norm( full(Phi_hat_11^n) )];
%     end
%     figure(1) 
%     plot(1:N, Phi_hat_11_powers, 'o-', 'LineWidth', 2, 'DisplayName', sprintf('$n_t = %d$', nt))
%     hold on
%     lh = legend();

       %Phi_hat_11_powers = [Phi_hat_11_powers; norm( full(Phi_hat_11^nt) )];

       
   %temp = [temp; (norm( full(Phi_hat_11) )^nt - 1) / (norm(full(Phi_hat_11)) - 1)];
       
    
    %% Setup instance of space-time system class
    % The class expects a step function in disc pa. 
    Phi_11_step = @(t0_idx, q0) Phi_hat_11 * q0;
    Phi_12_step = @(t0_idx, q0) Phi_hat_12 * q0;
    Phi_21_step = @(t0_idx, q0) Phi_hat_21 * q0;
    Phi_22_step = @(t0_idx, q0) Phi_hat_22 * q0;
    Phi_11_step_transp = @(t0_idx, q0) Phi_hat_11.' * q0;
    Phi_12_step_transp = @(t0_idx, q0) Phi_hat_12.' * q0;
    Phi_21_step_transp = @(t0_idx, q0) Phi_hat_21.' * q0;
    Phi_22_step_transp = @(t0_idx, q0) Phi_hat_22.' * q0;

    st_system_11 = one_step_st_system(Phi_11_step, nx, nt, solve_pa);
    st_system_11.step_transp = Phi_11_step_transp;
    
    st_system_12 = one_step_st_system(Phi_12_step, nx, nt, solve_pa);
    st_system_12.step_transp = Phi_12_step_transp;
    
    st_system_21 = one_step_st_system(Phi_21_step, nx, nt, solve_pa);
    st_system_21.step_transp = Phi_21_step_transp;
    
    st_system_22 = one_step_st_system(Phi_22_step, nx, nt, solve_pa);
    st_system_22.step_transp = Phi_22_step_transp;
%     
%     
%     sigma_max = svds(@(b, label) Delta_action(st_system_11, st_system_12, st_system_21, st_system_22, b, label), ...
%         [nx*nt, nx*nt], 1, 'largest')
    
%     sigma_min = svds(@(b, label) Delta_action(st_system_11, st_system_12, st_system_21, st_system_22, b, label), ...
%         [nx*nt, nx*nt], 1, 'smallest');
    
    sigma    = svds(@(b, label) A_diag_block_inv_action(st_system_11, b, label), [nx*nt, nx*nt], 1);
    sigma_11 = [sigma_11; sigma];
%     
%     
%     sigma    = svds(@(b, label) A_diag_block_inv_action(st_system_22, b, label), [nx*nt, nx*nt], 1);
%     sigma_22 = [sigma_22; sigma];
%     
%     sigma_12 = [sigma_12; norm(Phi_hat_12, inf)];
%     sigma_21 = [sigma_21; norm(Phi_hat_21, inf)];
    
end
% End of looping over different spatial resolutions.

%[temp, (1 + Phi_hat_11_powers).*temp]
%temp(2:end)./temp(1:end-1)


figure(3)
hold on
%plot(nx_array, temp * sigma_11(end) / temp(end), '-c*')
plot(nx_array, temp, '-c*')
plot(nx_array, temp2, '--r>')
hold on
%figure
plot(nx_array, sigma_11, '-bo')

%sigma = sigma_22 .* sigma_21 .* sigma_11 .* sigma_12;

%kappa = (1 + sigma) ./ (1 - sigma);

%figure(2)
%semilogx(nx_array, kappa, '-o' ,'DisplayName', sprintf('$T = %.2f$', mesh_pa.tmax))
xlabel('$n_x$')
%title('Cond number bound')
lh = legend();

% P = I - Delta, 
% Delta = (A_{22}^-1 * A_{21}) * (A_{11}^-1 * A_{12})
% Delta^T = (A_{11}^-1 * A_{12})^T * (A_{22}^-1 * A_{21})^T
function x = Delta_action(st_system_11, st_system_12, st_system_21, st_system_22, b, label)
    if strcmp(label, 'notransp')
        % Comptute x = Delta * b.
        x = b;
        x = A_off_diag_block_action(st_system_12, x, label);
        x = A_diag_block_inv_action(st_system_11, x, label);
        
        x = A_off_diag_block_action(st_system_21, x, label);
        x = A_diag_block_inv_action(st_system_22, x, label);
        
        %x = b - x;
        
    elseif strcmp(label, 'transp')
        
        % Comptute x = Delta' * b.
        x = b;
        x = A_off_diag_block_action(st_system_22, x, label);
        x = A_diag_block_inv_action(st_system_21, x, label);
        
        x = A_diag_block_inv_action(st_system_11, x, label);
        x = A_off_diag_block_action(st_system_12, x, label);

        %x = b - x;
    end
end

function x = A_diag_block_inv_action(st_system, b, label)
    if strcmp(label, 'notransp')
        x = st_system.forward_solve(b);
    elseif strcmp(label, 'transp')
        x = st_system.backward_solve(b);
    end
end

% The off-diag blocks A_{12} and A_{21} have a block matrix on their block
% subdiagonal. So to compute their action we can just apply the action and
% then subtract out the vecctor to account for the identity diagonal in 
% the action code. 
function x = A_off_diag_block_action(st_system, b, label)
    if strcmp(label, 'notransp')
        x = st_system.system_action(b);
        x = x - b;
    elseif strcmp(label, 'transp')
        x = st_system.system_transp_action(b);
        x = x - b;
    end
end
