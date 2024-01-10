% Used for saving figures to nice high-resolution pngs.
function figure_saver(fig, fig_name, save_native_fig)
    fig.PaperPositionMode = 'auto';
    fig_pos = fig.PaperPosition;
    fig.PaperSize = [fig_pos(3) fig_pos(4)];
    set(gcf, 'Color', 'w'); % Otherwise saved fig will have grey background
    export_fig(strcat(fig_name, '.png'), '-m4')
    
    % Store MATLAB figure too?
    if nargin == 3 
        if save_native_fig; saveas(gcf, strcat(fig_name, '.fig')); end
    end
end