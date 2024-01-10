% I just saved the over-solving figures as .figs then here I create
% high-res png versions of them.

tmp = matlab.desktop.editor.getActive;
cd(fileparts(tmp.Filename));

clc
clear 
close all 



open_fig_dir = '../figures/paper/over-solving/';
fig_names = {
    'rich-inexact-conv-k0', ...
    'rich-inexact-conv-k2', ...
    'rich-inexact-conv-k4', ...
    'rich-inexact-conv-k6', ...
    'time-stepping-sol-disc-error',  ...
    'time-stepping-sol'};

for idx = 1:numel(fig_names)
    figure_name = sprintf('%s%s', open_fig_dir, fig_names{idx}); 
    fh = openfig(figure_name);
    figure_saver(fh, figure_name)
    close(fh)
end




% Save the figure
function figure_saver(fig, fig_name)
    fig.PaperPositionMode = 'auto';
    fig_pos = fig.PaperPosition;
    fig.PaperSize = [fig_pos(3) fig_pos(4)];
    set(gcf, 'Color', 'w'); % Otherwise saved fig will have grey background
    export_fig(strcat(fig_name, '.png'), '-m4')
    %saveas(gcf, strcat(fig_name, '.fig'));
end