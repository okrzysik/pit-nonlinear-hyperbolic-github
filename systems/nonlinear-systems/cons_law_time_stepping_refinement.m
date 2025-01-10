% Time-stepping on SWE or Euler equations.
%
% Just uncomment below the PDE+domain+initial condition combination below 
% that you want to solve and then run the script.
%
% This function does largely the same thing as "cons_law_time_stepping.m"
% except that it plots cross-sections of the solution at two different mesh
% resolutions. The coarser resolution solutions are shown as colored lines,
% while the higher resolution are shown as thin black lines. This is to
% provide some sense of accuracy of the discretization.
%
% NOTE: This code doesn't have the full functionality of
% "cons_law_time_stepping.m" in terms of plotting and saving solution
% plots. It only draws cross-sections of the plots and it uses a local
% function to plot the cross-sections.

tmp = matlab.desktop.editor.getActive;
cd(fileparts(tmp.Filename));

clc
clear
close all

%% Plotting options
figno = 2;

fig_dir = 'figures/paper/';
save_sol_figs = ~true;

%% Discretization parameters

%nx_array = [2.^(11); 2.^(13)]; % Just two values in this array. 
nx_array = [2.^(13); 2.^(11)]; % Just two values in this array. 

disc_pa.num_flux_id = 'LLF';
disc_pa.num_flux_id = 'ROE'; disc_pa.delta_smoothing_parameter = 1e-6;


%% PDE, domain and initial condition parameters
%% SWE: Initial depth perturbation with Gaussian of size epsilon
%% Figures 3 and 4
% pde_pa.pde_id = 'shallow-water';
% pde_pa.ic_id = 'idp1'; 
% pde_pa.bcs = 'periodic';
% 
% Figure 3
% %pde_pa.ic_epsilon = 0.1; % Only weak shocks form here.
% Figure 4
% pde_pa.ic_epsilon = 0.6; % Shocks form.
% 
% mesh_pa.xmin = -5;
% mesh_pa.xmax =  5;
% disc_pa.CFL_number = 0.8; % max-wave-speed * dt/h
% mesh_pa.tmax = 10;


%% Euler: Smooth initial density and pressure perturbation
%% Figures 5 and 6
pde_pa.pde_id = 'euler'; 
pde_pa.ic_id = 'idp1'; 

% Figure 5
pde_pa.ic_epsilon = 0.2;
% Figure 6
%pde_pa.ic_epsilon = 1.2;

mesh_pa.xmin = -5;
mesh_pa.xmax =  5;
disc_pa.CFL_number = 0.7; % max-wave-speed * dt/h
mesh_pa.tmax = 10;
pde_pa.bcs = 'periodic';


%% Loop over different mesh resolutions
for nx_idx = 1:numel(nx_array)

    nx = nx_array(nx_idx);

    
    % Set the same max wave-speed for all resolutions by finding that for
    % the highest resolution mesh. Since the time-steps are proportional to
    % this, this allows us to have a nested temporal grid across levels,
    % and hence overlay solutions are two different mesh resolutions
    if nx_idx == 1
        nx_max = max(nx_array);
        % Create spatial mesh
        mesh_pa.nx = nx_max;
        [mesh_pa.x_interfaces, mesh_pa.x_centers, mesh_pa.h] = space_mesh(mesh_pa.xmin, mesh_pa.xmax, mesh_pa.nx);
        
        % Get initial conditon.
        if strcmp(pde_pa.pde_id, 'shallow-water')
            my_cons_law = swe_system(pde_pa);
            [h0, u0, ic_description] = my_cons_law.initial_condition(mesh_pa.x_centers);
            q0 = [h0; h0.*u0];
            
        elseif strcmp(pde_pa.pde_id, 'euler')
            my_cons_law = euler_system(pde_pa);
            [rho0, u0, E0, ic_description] = my_cons_law.initial_condition(mesh_pa.x_centers);
            q0 = [rho0; rho0.*u0; E0];
            
        end
        
        % Compute fastest wave-speed associated with initial condition
        my_cons_law.mesh_pa = mesh_pa; % The wave-speed function requires this to be set
        lambda0 = my_cons_law.wave_speeds(q0);
        abs_fprime0_max = max(abs(lambda0(:)));

    end

    %% Set up mesh
    % Create spatial mesh
    mesh_pa.nx = nx;
    [mesh_pa.x_interfaces, mesh_pa.x_centers, mesh_pa.h] = space_mesh(mesh_pa.xmin, mesh_pa.xmax, mesh_pa.nx);
    
    % Get initial conditon.
    if strcmp(pde_pa.pde_id, 'shallow-water')
        my_cons_law = swe_system(pde_pa);
        [h0, u0] = my_cons_law.initial_condition(mesh_pa.x_centers);
        q0 = [h0; h0.*u0];
        
    elseif strcmp(pde_pa.pde_id, 'euler')
        my_cons_law = euler_system(pde_pa);
        [rho0, u0, E0] = my_cons_law.initial_condition(mesh_pa.x_centers);
        q0 = [rho0; rho0.*u0; E0];
        
    end

    
    % Get temporal mesh.
    [mesh_pa.t, mesh_pa.dt] = time_mesh(mesh_pa.h, disc_pa.CFL_number, abs_fprime0_max, 1, mesh_pa.tmax);
    mesh_pa.nt = numel(mesh_pa.t);
    
    %% Finalize construction of PDE class
    my_cons_law.disc_pa = disc_pa; % Set this again so that it's updated.
    my_cons_law.mesh_pa = mesh_pa;
    
    nx = mesh_pa.nx;
    nt = mesh_pa.nt;
    
    fprintf('\nnx=%d, nt=%d\n', mesh_pa.nx, mesh_pa.nt);
    
    % Initialize solution vector
    q = zeros(my_cons_law.m * nx, nt);
    
    % Insert initial condition
    q(:, 1) = q0;
    
    %% Time-step!
    for n = 1:mesh_pa.nt-1
        q(:, n+1) = my_cons_law.step(n, q(:, n));
    end
    fprintf('Time-stepping done...\n')
    
    % At the coarsest mesh resolution pick out num_cross_sec times at which to
    % plot the solution.
    t    = mesh_pa.t;
    nt   = mesh_pa.nt;
    tmax = t(end);
    %t_plot_times = [0; 2.5; 5; 7.5; 10];
    t_plot_times = [2.5; 5; 7.5; 10]; % Ignore t=0 since this doesn't change under mesh refinement and just obscures the details of the figure
    
    % Find indices of time points in t_plot_times
    t_plot_inds = [];
    for ind = 1:numel(t_plot_times)
        [~, t_plot_ind] = min(abs(t - t_plot_times(ind)));
        t_plot_inds = [t_plot_inds; t_plot_ind];
    end
    
    
    % Unpack first two variables
    if strcmp(pde_pa.pde_id, 'shallow-water')
            H  = q(1:nx, :);
            HU = q(nx+1:2*nx, :);
    
            Q1 = H;
            Q2 = HU;
    
    elseif strcmp(pde_pa.pde_id, 'euler')
        RHO = q(1:nx, :);
        U   = q(nx+1:2*nx, :) ./ RHO;
        E   = q(2*nx+1:3*nx, :);
    
        Q1 = RHO;
        Q2 = U;
    end

    % Create figures
    if nx_idx == 1
        fh1 = figure();
        fh2 = figure();
    end
    my_cross_sec_plot(fh1, t_plot_inds, Q1, nx == max(nx_array), mesh_pa);
    my_cross_sec_plot(fh2, t_plot_inds, Q2, nx == max(nx_array), mesh_pa);
    
end

% Pretty up plots


plot_pa.title_str = sprintf('$\\varepsilon = %.1f$', pde_pa.ic_epsilon);
% Get appropriate title strings and saving names
if strcmp(pde_pa.pde_id, 'shallow-water')
    fh1_title_str  = sprintf('$h(x, t)$: %s', plot_pa.title_str);
    fh2_title_str = sprintf('$hu(x, t)$: %s', plot_pa.title_str);

    plot_pa.save_dir = strcat(fig_dir, '/SWE/', pde_pa.ic_id);
    fh1_save_name  = sprintf('%s/SWE-h-cs-eps%.1f-nx%d-%d-ref', plot_pa.save_dir, pde_pa.ic_epsilon, min(nx_array), max(nx_array));
    fh2_save_name = sprintf('%s/SWE-hu-cs-eps%.1f-nx%d-%d-ref', plot_pa.save_dir, pde_pa.ic_epsilon, min(nx_array), max(nx_array));
elseif strcmp(pde_pa.pde_id, 'euler')
    fh1_title_str  = sprintf('$\\rho(x, t)$: %s', plot_pa.title_str);
    fh2_title_str = sprintf('$u(x, t)$: %s', plot_pa.title_str);

    plot_pa.save_dir = strcat(fig_dir, '/euler/', pde_pa.ic_id);
    fh1_save_name  = sprintf('%s/euler-rho-cs-eps%.1f-nx%d-%d-ref', plot_pa.save_dir, pde_pa.ic_epsilon, min(nx_array), max(nx_array));
    fh2_save_name = sprintf('%s/euler-u-cs-eps%.1f-nx%d-%d-ref', plot_pa.save_dir, pde_pa.ic_epsilon, min(nx_array), max(nx_array));
end

% Add legends, axis labels, and axis ticks.
figure(fh1)
title(fh1_title_str)
pretty_cs_figure(fh1, Q1(:, t_plot_inds), mesh_pa)
figure(fh2)
title(fh2_title_str)
pretty_cs_figure(fh2, Q2(:, t_plot_inds), mesh_pa)


if save_sol_figs
    figure_saver(figure(fh1), fh1_save_name, false);
    figure_saver(figure(fh2), fh2_save_name, false);
end

function my_cross_sec_plot(fh, t_plot_inds, Q, nx_max, mesh_pa)
% Plot cross sections of Q and the indices provided on figure fh. If nx ==
% nx_max the cross-sections are grey, while for nx not equal to this they
% are colored

    t = mesh_pa.t;
    ls = {"-"; "--"; ":"; "-."};
    figure(fh);
    hold on
    for ind = 1:numel(t_plot_inds)
        t_plot_ind = t_plot_inds(ind);

        ls_ind = mod(ind, 4)+1;

        if nx_max
            plot(mesh_pa.x_centers, Q(:, t_plot_ind), 'HandleVisibility', 'Off', ...
                'LineStyle', '-', ...
                'LineWidth', 2, ...
                'Color', 'k')
        else
            nameval = lineprops(ind+1);
            plot(mesh_pa.x_centers, Q(:, t_plot_ind), ...
            nameval{:}, ...
            'LineStyle', ls{ls_ind}, ...
            'Marker', 'None', ...
            'LineWidth', 3, ...
            'DisplayName', sprintf('$t = %.2f$', t(t_plot_ind)))
        end

    end

end


function pretty_cs_figure(fh, plotted_data, mesh_pa)

    figure(fh)
    lh = legend();
    lh.set('Location', 'Best', 'NumColumns', 2) 
    box on
    xlabel('$x$')
    
    % Adjust bounding box of figure.
    data_min = min(plotted_data(:));
    data_max = max(plotted_data(:));
    height   = data_max - data_min;
    tol      = 0.05;
    ylim([data_min - tol*height, data_max + tol*height])
    xmin     = mesh_pa.xmin;
    xmax     = mesh_pa.xmax;
    xlim([xmin, xmax])
    xticks([xmin:(xmax-xmin)/10:xmax])
    xticks([xmin:2*(xmax-xmin)/10:xmax])
end
