function cs_fig = cross_sec_plot(figno, data, label, pa)
    
    t = pa.mesh_pa.t;
    tmax = t(end);
    
    % Number of cross sections to plot
    if isfield(pa, 'num_cross_sec')
        num_cross_sec = pa.num_cross_sec; 
    else
        num_cross_sec = 5; 
    end 
    if num_cross_sec == 1
        t_plot_times = tmax;
    else
        t_plot_times = [0:tmax/(num_cross_sec - 1):tmax];
    end

    % Line styles
    ls = {"-"; "--"; ":"; "-."};

    plotted_data = [];
    for t_plot_idx = 1:numel(t_plot_times)

        t_plot_approx = t_plot_times(t_plot_idx);

        [~, t_idx] = min(abs(t - t_plot_approx));
        t_plot = t(t_idx);
        %[t_plot, t_plot_approx]

        nameval = lineprops(t_plot_idx);
        % h
        cs_fig = figure(figno);
        plot(pa.mesh_pa.x_centers, data(:, t_idx), ...
            nameval{:}, ...
            'Marker', 'None', ...
            'LineWidth', 2, ...
            'LineStyle', ls{mod(t_plot_idx-1, numel(ls))+1}, ...
            'DisplayName', sprintf('$t = %.2f$', t_plot))
        hold on
        plotted_data = [plotted_data; data(:, t_idx)];
    end
    
    % Adjust bounding box of figure.
    data_min = min(plotted_data(:));
    data_max = max(plotted_data(:));
    height   = data_max - data_min;
    tol      = 0.05;
    ylim([data_min - tol*height, data_max + tol*height])
    xmin     = pa.mesh_pa.xmin;
    xmax     = pa.mesh_pa.xmax;
    xlim([xmin, xmax])
    xticks([xmin:(xmax-xmin)/10:xmax])
    xticks([xmin:2*(xmax-xmin)/10:xmax])
    
    lh = legend();
    lh.set('Location', 'Best') 
    box on
    hold on
    %ylabel(label)
    title(label)
    xlabel('$x$')

    if pa.invisible_col_bar
        cb = colorbar(); scale = 0.998; cb.Color = scale*[1 1 1]; myColorMap = scale*ones(1, 3); colormap(gcf, myColorMap); 
    end

end