% Plot solution and solution-related quantities for SWE or Euler
% 
function [fh, figno] = plot_cons_law_system_solutions(q, cons_law, pde_pa, disc_pa, mesh_pa, figno, plot_pa)


    %% Plots of solution, etc 
    [mesh_pa.X, mesh_pa.T] = meshgrid(mesh_pa.x_centers, mesh_pa.t);
    plot_pa.disc_pa = disc_pa;
    plot_pa.mesh_pa = mesh_pa;
    
    % Reshape data if needed
    nx = mesh_pa.nx;
    [~, ncols] = size(q);
    if ncols == 1
        q = reshape(q, [cons_law.m * nx, mesh_pa.nt]);
    end
    
    fh = [];

    %% Plotting customized for SWE
    if strcmp(pde_pa.pde_id, 'shallow-water')
        H  = q(1:nx, :);
        HU = q(nx+1:2*nx, :);
        LAMBDA = zeros(size(q));
        for n = 1:mesh_pa.nt
            LAMBDA(:, n) = cons_law.wave_speeds(q(:, n));
        end
        LAMBDA1 = LAMBDA(1:nx, :);
        LAMBDA2 = LAMBDA(nx+1:2*nx, :);

        mesh_nx_tol = 2*1024;
        if mesh_pa.nx > mesh_nx_tol
            plot_pa.meshing_down_sample_factor = mesh_pa.nx / mesh_nx_tol;
        end

        mesh_nx_tol = 512;
        if mesh_pa.nx > mesh_nx_tol
            plot_pa.meshing_down_sample_factor = mesh_pa.nx / mesh_nx_tol;
        end
        
        if isfield(plot_pa, 'title_str')
            h_title_str = sprintf('$h(x, t)$: %s', plot_pa.title_str);
            hu_title_str = sprintf('$hu(x, t)$: %s', plot_pa.title_str);
        else
            h_title_str = sprintf('$h(x, t)$');
            hu_title_str = sprintf('$hu(x, t)$');
        end
            

        if plot_pa.plot_solution_contours
            h_con_fig  = contour_plot(figno, H, h_title_str, plot_pa); figno = figno+1;    
            hu_con_fig = contour_plot(figno, HU, hu_title_str, plot_pa); figno = figno+1;
            
            if isfield(plot_pa, 'save_figs') && plot_pa.save_figs
                h_save_name  = sprintf('%s/SWE-h-con-%s-nx%d', plot_pa.save_dir, plot_pa.save_str, nx);
                figure_saver(figure(h_con_fig), h_save_name, false);

                hu_save_name = sprintf('%s/SWE-hu-con-%s-nx%d', plot_pa.save_dir, plot_pa.save_str, nx);
                figure_saver(figure(hu_con_fig), hu_save_name, false);
            end
        end

        if plot_pa.plot_solution_cross_sections
            h_cs_fig  = cross_sec_plot(figno, H, h_title_str, plot_pa); figno = figno+1;
            hu_cs_fig = cross_sec_plot(figno, HU, hu_title_str, plot_pa); figno = figno+1;
            
            if isfield(plot_pa, 'save_figs') && plot_pa.save_figs
                h_save_name = sprintf('%s/SWE-h-cs-%s-nx%d', plot_pa.save_dir, plot_pa.save_str, nx);
                figure_saver(figure(h_cs_fig), h_save_name, false);

                hu_save_name = sprintf('%s/SWE-hu-cs-%s-nx%d', plot_pa.save_dir, plot_pa.save_str, nx);
                figure_saver(figure(hu_cs_fig), hu_save_name, false);
            end
        end

        if plot_pa.plot_wave_speed_contours
            lambda1_con_fig = contour_plot(figno, LAMBDA1, '$\lambda^1(x, t)$', plot_pa); figno = figno+1;
            lambda2_con_fig = contour_plot(figno, LAMBDA2, '$\lambda^2(x, t)$', plot_pa); figno = figno+1;
        end
        
        if plot_pa.plot_wave_speed_cross_sections
            lambda1_cs_fig = cross_sec_plot(figno, LAMBDA1, '$\lambda^1(x, t)$', plot_pa); figno = figno+1;
            lambda2_cs_fig = cross_sec_plot(figno, LAMBDA2, '$\lambda^2(x, t)$', plot_pa); figno = figno+1;
        end
        
    %% Plotting customized for Euler equations
    elseif strcmp(pde_pa.pde_id, 'euler')
        
        RHO = q(1:nx, :);
        U   = q(nx+1:2*nx, :) ./ RHO;
        E   = q(2*nx+1:3*nx, :);
        P   = cons_law.pressure(RHO, U, E);
        LAMBDA = zeros(size(q));
        for n = 1:mesh_pa.nt
            LAMBDA(:, n) = cons_law.wave_speeds(q(:, n));
        end
        LAMBDA1 = LAMBDA(1:nx, :);
        LAMBDA2 = LAMBDA(nx+1:2*nx, :);
        LAMBDA3 = LAMBDA(2*nx+1:3*nx, :);

        mesh_nx_tol = 2*1024;
        if mesh_pa.nx > mesh_nx_tol
            plot_pa.meshing_down_sample_factor = mesh_pa.nx / mesh_nx_tol;
        end

        mesh_nx_tol = 2*1024;
        if mesh_pa.nx > mesh_nx_tol
            plot_pa.meshing_down_sample_factor = mesh_pa.nx / mesh_nx_tol;
        end

        if isfield(plot_pa, 'title_str')
            rho_title_str = sprintf('$\\rho(x, t)$: %s', plot_pa.title_str);
            u_title_str   = sprintf('$u(x, t)$: %s', plot_pa.title_str);
            p_title_str   = sprintf('$p(x, t)$: %s', plot_pa.title_str);
        else
            rho_title_str = sprintf('$\\rho(x, t)$');
            u_title_str   = sprintf('$u(x, t)$');
            p_title_str   = sprintf('$p(x, t)$');
        end

        if plot_pa.plot_solution_contours
            rho_con_fig = contour_plot(figno, RHO, rho_title_str, plot_pa);  figno = figno+1;
            u_con_fig   = contour_plot(figno, U,   u_title_str,    plot_pa); figno = figno+1;
            p_con_fig   = contour_plot(figno, P,   p_title_str,    plot_pa); figno = figno+1;

            if isfield(plot_pa, 'save_figs') && plot_pa.save_figs
                rho_save_name = sprintf('%s/euler-rho-con-%s-nx%d', plot_pa.save_dir, plot_pa.save_str, nx);
                figure_saver(figure(rho_con_fig), rho_save_name, false);

                u_save_name = sprintf('%s/euler-u-con-%s-nx%d', plot_pa.save_dir, plot_pa.save_str, nx);
                figure_saver(figure(u_con_fig), u_save_name, false);

                p_save_name = sprintf('%s/euler-p-con-%s-nx%d', plot_pa.save_dir, plot_pa.save_str, nx);
                figure_saver(figure(p_con_fig), p_save_name, false);
            end
        end

        if plot_pa.plot_solution_cross_sections
            rho_cs_fig = cross_sec_plot(figno, RHO, rho_title_str, plot_pa); figno = figno+1;
            u_cs_fig   = cross_sec_plot(figno, U,   u_title_str,   plot_pa); figno = figno+1;
            p_cs_fig   = cross_sec_plot(figno, P,   p_title_str,   plot_pa); figno = figno+1;

            if isfield(plot_pa, 'save_figs') && plot_pa.save_figs
                rho_save_name = sprintf('%s/euler-rho-cs-%s-nx%d', plot_pa.save_dir, plot_pa.save_str, nx);
                figure_saver(figure(rho_cs_fig), rho_save_name, false);

                u_save_name = sprintf('%s/euler-u-cs-%s-nx%d', plot_pa.save_dir, plot_pa.save_str, nx);
                figure_saver(figure(u_cs_fig), u_save_name, false);

                p_save_name = sprintf('%s/euler-p-cs-%s-nx%d', plot_pa.save_dir, plot_pa.save_str, nx);
                figure_saver(figure(p_cs_fig), p_save_name, false);
            end
        end

        if plot_pa.plot_wave_speed_contours
            lambda1_con_fig = contour_plot(figno, LAMBDA1, '$\lambda^1(x, t)$', plot_pa); figno = figno+1;
            lambda2_con_fig = contour_plot(figno, LAMBDA2, '$\lambda^2(x, t)$', plot_pa); figno = figno+1;
            lambda3_con_fig = contour_plot(figno, LAMBDA3, '$\lambda^3(x, t)$', plot_pa); figno = figno+1;
        end
        
        if plot_pa.plot_wave_speed_cross_sections
            lambda1_cs_fig = cross_sec_plot(figno, LAMBDA1, '$\lambda^1(x, t)$', plot_pa); figno = figno+1;
            lambda2_cs_fig = cross_sec_plot(figno, LAMBDA2, '$\lambda^2(x, t)$', plot_pa); figno = figno+1;
            lambda3_cs_fig = cross_sec_plot(figno, LAMBDA3, '$\lambda^3(x, t)$', plot_pa); figno = figno+1;
        end
        
    end

end