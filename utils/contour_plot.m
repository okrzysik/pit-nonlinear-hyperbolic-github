function con_fig = contour_plot(figno, data, label, pa)

    con_tol  = pa.con_tol;
    ncon_lvl = pa.ncon_lvl;

    if isfield(pa, 'meshing_down_sample_factor')
        ds_fac = pa.meshing_down_sample_factor;
    else
        ds_fac = 1;
    end

    data_min  = min(data(:));
    data_max  = max(data(:));

    con_fig = figure(figno);
    mesh(pa.mesh_pa.X(1:ds_fac:end, 1:ds_fac:end), pa.mesh_pa.T(1:ds_fac:end, 1:ds_fac:end), data(1:ds_fac:end, 1:ds_fac:end)'); view(2);
    axis tight
    xlabel('$x$')
    ylabel('$t$')
    title(label)
    con0 = data_min + con_tol; 
    con1 = data_max - con_tol;
    ch = (con1 - con0)/ncon_lvl;
    hold on
    contour3(pa.mesh_pa.X(1:ds_fac:end, 1:ds_fac:end), pa.mesh_pa.T(1:ds_fac:end, 1:ds_fac:end), data(1:ds_fac:end, 1:ds_fac:end)', [con0:ch:con1] ,'linew', 2, 'color', 'k')
    cb = colorbar();
    xmin = pa.mesh_pa.xmin;
    xmax = pa.mesh_pa.xmax;
    xlim([xmin, xmax])
    xticks([xmin:(xmax-xmin)/10:xmax])
    xticks([xmin:2*(xmax-xmin)/10:xmax])
end