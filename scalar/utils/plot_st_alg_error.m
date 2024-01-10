% Make a plot of the algebraic space-time error associated with the current
% Richardson iterate.
function [caxis0, cb_limits0] = plot_st_alg_error(u_bar, u_bar_exact, mesh_pa, iter_richardson, figno, caxis0, cb_limits0)
    error_mask_tol = 10^-8; % Errors < tol are masked to zero to help with visualization
    
    u_bar = reshape(u_bar, [mesh_pa.nx, mesh_pa.nt]);
    % Algebraic error
    alg_error = u_bar - u_bar_exact; 
    % Compute its norms
    alg_error_infnorm = norm(alg_error(:), inf);
    alg_error_1norm   = norm(alg_error(:), 1) * mesh_pa.h * mesh_pa.dt;
    
    [X, T] = meshgrid(mesh_pa.x_centers, mesh_pa.t);
    figure(figno+10+iter_richardson)
    %mesh(X, T, error'); view(2)
    logged_alg_error = log10(abs(alg_error));
    logged_alg_error(logged_alg_error < log10(error_mask_tol)) = NaN;
    mesh(X, T, logged_alg_error'); view(2); view(-25, 45)
    box on
    cb = colorbar();
    cb.Limits = [log10(error_mask_tol), 0];
    
    if iter_richardson == 0
        cb_limits0 = cb.Limits;
        caxis0     = caxis;
    else
        cb.Limits = cb_limits0;
        caxis(caxis0);
    end
    
    % Make title look nice, including error norms.
    [m1, e1] = get_scientific_decomposition(alg_error_1norm);   e1_str   = sprintf('%.0f \\times 10^{%d}', m1, e1);
    [m2, e2] = get_scientific_decomposition(alg_error_infnorm); einf_str = sprintf('%.0f \\times 10^{%d}', m2, e2);
    e_label_str   = '\Vert \bar{\mathbf{e}}_{\fontsize{8}{0}\selectfont\textrm{alg}} \Vert({1},\infty)';
    title(sprintf('$k=%d:%s=(%s,%s)$', iter_richardson, e_label_str, e1_str, einf_str))
    xlabel('$x$')
    ylabel('$t$')
    zlabel('$\log_{10} |\bar{e}_{\fontsize{8}{0}\selectfont\textrm{alg}}|$')
    axis tight
    ylim([0 max(mesh_pa.t)]) % Since function is zero at t=0, that tick label disappears 
    box on
    
    zticks([log10(error_mask_tol):2:0])
    zlim([log10(error_mask_tol), 0])
end