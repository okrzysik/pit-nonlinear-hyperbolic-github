% This file is for taking two (or potentially 3) native MATLAB figures,
% merging their contents into a single figure and then saving the result.

tmp = matlab.desktop.editor.getActive;
cd(fileparts(tmp.Filename));

clc
clear 
close all 

save_fig = ~true;

upper_xlim = 20;

% Merge plots for DIRECT solution of linear problems using F- and NO-relaxation
suf_new_fig = '-DIRECT';
suf_fig1 = '-F(8)-DIRECT';
suf_fig2 = '-DIRECT';
open_fig_dir = '../figures/paper/exact/';

% % Merge plots for DIRECT and APPROXIMATE solutions of linear problems
% suf_new_fig = '-F(8)-DIRECT-MGRIT';
% suf_fig2 = '-F(8)-MGRIT-F(8)-LU-2';
% suf_fig1 = '-F(8)-DIRECT';
% %suf_fig3 = '-F(8)-MGRIT-F(8)-NONE-2'; % An example where there's no truncation correction
% open_fig_dir = '../figures/paper/inexact/';

% Manyally choose legend location
leg_loc = 'SouthWest';
leg_loc = 'NorthEast';

pref = 'B(3)-GLF-lin(1)';  show_y_label = true;  show_legend = true; 
%pref = 'B(3)-LLF-lin(1)';  show_y_label = true;  show_legend = false;
pref = 'BL(3)-GLF-lin(1)'; show_y_label = false; show_legend = false;
pref = 'BL(3)-LLF-lin(1)'; show_y_label = false; show_legend = false;

pref = 'B(3)-LLF-WENO(3,P)';  show_y_label = true; show_legend = true; 
pref = 'B(3)-LLF-WENO(3,N)';  show_y_label = true; show_legend = false;
pref = 'BL(3)-LLF-WENO(3,P)'; show_y_label = false; show_legend = false;
pref = 'BL(3)-LLF-WENO(3,N)'; show_y_label = false; show_legend = false;

pref = 'B(3)-GLF-WENO(3,P)';  show_y_label = true; show_legend = true; 
pref = 'B(3)-GLF-WENO(3,N)';  show_y_label = true; show_legend = false;
pref = 'BL(3)-GLF-WENO(3,P)'; show_y_label = false; show_legend = false;
pref = 'BL(3)-GLF-WENO(3,N)'; show_y_label = false; show_legend = false;


save_fig_name = sprintf('%s/merged/%s%s', open_fig_dir, pref, suf_new_fig);

fig1_name = sprintf('%s%s%s', open_fig_dir, pref, suf_fig1);
fh1 = openfig(fig1_name);
l1  = findall(get(fh1, 'Number'), 'type','legend');
delete(l1)
%if ~isempty(l1);

fig2_name = sprintf('%s%s%s.fig', open_fig_dir, pref, suf_fig2);
fh2 = openfig(fig2_name);
l2  = findall(get(fh2, 'Number'), 'type','legend');
%delete(l2)


% FINDALL and not FINDOBJ because after saving the line properties of the
% figures are empty. 
copyobj(findall(get(fh1, 'Number'), 'type', 'line'), findall(get(fh2, 'Number'), 'type', 'axes'))

figure(fh2)
if ~show_y_label
    ylabel('')
    yticklabels([])
end

if ~show_legend
    l2 = findall(get(fh2, 'Number'), 'type','legend');
    delete(l2)
else
   l2.String(1:numel(l2.String)/2) = '';
   set(l2, 'Location', leg_loc)
end


% In some special cases merge in a 3rd set of lines too.
if exist('suf_fig3','var')
    xl = xlim();
    fig3_name = sprintf('%s%s%s.fig', open_fig_dir, pref, suf_fig3);
    fh3 = openfig(fig3_name);
    l3  = findall(get(fh3, 'Number'), 'type','legend');
    delete(l3)
    copyobj(findall(get(fh3, 'Number'), 'type', 'line'), findall(get(fh2, 'Number'), 'type', 'axes'))
    %figure(fh2)
    xlim(xl)
end


if exist('upper_xlim','var')
   xlim([0 upper_xlim]) 
end

grid minor
grid on
grid minor
%set(gca,'YMinorTick','Off')
%grid minor
    
if save_fig
    figure_saver(figure(fh2), save_fig_name)
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