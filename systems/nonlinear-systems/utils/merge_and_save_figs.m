% This file is for taking two native MATLAB figures,
% merging their contents into a single figure and then saving the result.

tmp = matlab.desktop.editor.getActive;
cd(fileparts(tmp.Filename));

clc
clear 
close all 

save_fig = ~true;

% Merge plots that live in the following directory:
%% SWE plots
%open_fig_dir = '../figures/paper/SWE/idp1/';
%open_fig_dir = '../figures/paper/SWE/dam-break/';
%% Euler plots
%open_fig_dir = '../figures/paper/euler/idp1/';
%open_fig_dir = '../figures/paper/euler/sod/';


% Manually choose legend location
%leg_loc = 'SouthWest';
leg_loc = 'NorthEast';

%% SWE, IDP small eps. Figure 3.
% Fig 3, bottom left
% suf_new_fig = '-bpeb-diag-it1';
% suf_fig1    = '-direct';
% suf_fig2    = '-bpeb-diag-it1';
% pref        = 'res-SWE-idp1-eps0.10'; show_y_label = true;  show_legend = true; 

% % Fig 3, bottom right
% suf_new_fig = '-bpap-diag-it1';
% suf_fig1    = '-bpab-MGRIT-diag-it1';
% suf_fig2    = '-bpab-diag-it1';
% pref        = 'res-SWE-idp1-eps0.10'; show_y_label = true;  show_legend = true; 


%% SWE, IDP big eps. Figure 4.
% % Fig 4, bottom left
% suf_new_fig = '-bpeb-diag-it1';
% suf_fig1    = '-direct';
% suf_fig2    = '-bpeb-diag-it1';
% pref        = 'res-SWE-idp1-eps0.60'; show_y_label = true;  show_legend = true; 

% % Fig 4, bottom middle
% suf_new_fig = '-bpeb-diag-it2';
% suf_fig1 = '-direct';
% suf_fig2 = '-bpeb-diag-it2';
% pref = 'res-SWE-idp1-eps0.60'; show_y_label = true;  show_legend = true; 

% Fig 4, bottom right
% suf_new_fig = '-bpap-diag-it2';
% suf_fig1 = '-bpab-MGRIT-diag-it2';
% suf_fig2 = '-bpab-diag-it2';
% pref = 'res-SWE-idp1-eps0.60'; show_y_label = true;  show_legend = true; 


%% Euler, IDP small eps. Figure 5
% % Fig 5, bottom left
% suf_new_fig = '-bpeb-diag-it1';
% suf_fig1 = '-direct';
% suf_fig2 = '-bpeb-diag-it1';
% pref = 'res-euler-idp1-eps0.20'; show_y_label = true;  show_legend = true; 

% % Fig 5, bottom right
% suf_new_fig = '-bpab-diag-it1';
% suf_fig1 = '-bpab-MGRIT-diag-it1';
% suf_fig2 = '-bpab-diag-it1';
% pref = 'res-euler-idp1-eps0.20'; show_y_label = true;  show_legend = true; 


%% Euler, IDP big eps. Figure 6
% % Fig 5, bottom left
% suf_new_fig = '-bpeb-diag-it2';
% suf_fig1 = '-direct';
% suf_fig2 = '-bpeb-diag-it2';
% pref = 'res-euler-idp1-eps1.20'; show_y_label = true;  show_legend = true; 

% % Fig 5, bottom middle
% suf_new_fig = '-bpeb-diag-it3';
% suf_fig1 = '-direct';
% suf_fig2 = '-bpeb-diag-it3';
% pref = 'res-euler-idp1-eps1.20'; show_y_label = true;  show_legend = true; 

% % Fig 5, bottom right
% suf_new_fig = '-bpab-diag-it3';
% suf_fig1 = '-bpab-MGRIT-diag-it3';
% suf_fig2 = '-bpab-diag-it3';
% pref = 'res-euler-idp1-eps1.20'; show_y_label = true;  show_legend = true; 


%% SWE, DB small eps. SM Figure 2
% % Fig SM 2, bottom left
% suf_new_fig = '-bpeb-diag-it1';
% suf_fig1 = '-direct';
% suf_fig2 = '-bpeb-diag-it1';
% pref = 'res-SWE-dam-break-eps0.10'; show_y_label = true;  show_legend = true; 

% % Fig SM 2, bottom right
% % This figure doesn't merge anything.

%% SWE, DB big eps. SM Figure 3
% % Fig SM 3, bottom left
% suf_new_fig = '-bpeb-diag-it1';
% suf_fig1 = '-direct';
% suf_fig2 = '-bpeb-diag-it1';
% pref = 'res-SWE-dam-break-eps2.00'; show_y_label = true;  show_legend = true; 
% 
% % Fig SM 2, bottom middle
% suf_new_fig = '-bpeb-diag-it3';
% suf_fig1 = '-direct';
% suf_fig2 = '-bpeb-diag-it3';
% pref = 'res-SWE-dam-break-eps2.00'; show_y_label = true;  show_legend = true;  

% % Fig SM 2, bottom right
% suf_new_fig = '-bpab-diag-it4';
% suf_fig1 = '-direct';
% suf_fig2 = '-bpab-diag-it4';
% pref = 'res-SWE-dam-break-eps2.00'; show_y_label = true;  show_legend = true; 


%% Euler, Sod small eps. SM Figure 4.
% % SM Fig 4, bottom left.
% suf_new_fig = '-bpeb-diag-it1';
% suf_fig1 = '-direct';
% suf_fig2 = '-bpeb-diag-it1';
% pref = 'res-euler-sod-eps0.12'; show_y_label = true;  show_legend = true; 

% % SM Fig 4, bottom right.
% suf_new_fig = '-bpab-diag-it1';
% suf_fig1 = '-direct';
% suf_fig2 = '-bpab-diag-it1';
% pref = 'res-euler-sod-eps0.12'; show_y_label = true;  show_legend = true; 

%% Euler, Sod big eps. SM Figure 5.
% SM Fig 5, bottom left.
suf_new_fig = '-bpeb-diag-it2';
suf_fig1 = '-direct';
suf_fig2 = '-bpeb-diag-it2';
pref = 'res-euler-sod-eps0.88'; show_y_label = true;  show_legend = true; 

% % SM Fig 5, bottom left.
% suf_new_fig = '-bpeb-diag-it4';
% suf_fig1 = '-direct';
% suf_fig2 = '-bpeb-diag-it4';
% pref = 'res-euler-sod-eps0.88'; show_y_label = true;  show_legend = true; 



save_fig_name = sprintf('%s/merged/%s%s', open_fig_dir, pref, suf_new_fig);

fig1_name = sprintf('%s%s%s', open_fig_dir, pref, suf_fig1, '.fig');
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
    export_fig(strcat(fig_name, '.png'), '-m2')
    %saveas(gcf, strcat(fig_name, '.fig'));
end