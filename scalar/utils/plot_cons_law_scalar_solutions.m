% Plot solution and solution-related quantities for SWE or Euler
% 
function [con_fig, cs_fig, figno] = plot_cons_law_scalar_solutions(u, cons_law, pde_pa, disc_pa, mesh_pa, figno, plot_pa)


    %% Plots of solution, etc 
    [mesh_pa.X, mesh_pa.T] = meshgrid(mesh_pa.x_centers, mesh_pa.t);
    plot_pa.disc_pa = disc_pa;
    plot_pa.mesh_pa = mesh_pa;
    
    % Reshape data if needed
    nx = mesh_pa.nx;
    [~, ncols] = size(u);
    if ncols == 1
        u = reshape(u, [nx, mesh_pa.nt]);
    end
    
    mesh_nx_tol = 512;
    if mesh_pa.nx > mesh_nx_tol
        plot_pa.meshing_down_sample_factor = mesh_pa.nx / mesh_nx_tol;
    end

    if isfield(plot_pa, 'ic_title_str')
        u_title_str = sprintf('%s: %s', cons_law.id_abbreviation, plot_pa.ic_title_str);
    else
        u_title_str = sprintf('%s', cons_law.id_abbreviation);
    end

    if plot_pa.plot_solution_contours
        con_fig  = contour_plot(figno, u, u_title_str, plot_pa); figno = figno+1;    

        if isfield(plot_pa, 'save_figs') && plot_pa.save_figs
            u_save_name = sprintf('%s/%s-con-%s-nx%d', plot_pa.save_dir, cons_law.id, plot_pa.ic_save_str, nx);
            figure_saver(figure(con_fig), u_save_name, false);
        end
    end

    if plot_pa.plot_solution_cross_sections
        cs_fig  = cross_sec_plot(figno, u, u_title_str, plot_pa); figno = figno+1;

        if isfield(plot_pa, 'save_figs') && plot_pa.save_figs
            u_save_name = sprintf('%s/%s-cs-%s-nx%d', plot_pa.save_dir, cons_law.id, plot_pa.ic_save_str, nx);
            figure_saver(figure(cs_fig), u_save_name, false);
        end
    end

        

end