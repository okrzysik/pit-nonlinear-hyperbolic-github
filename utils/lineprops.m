function nameval = lineprops(idx)
%   LINEPROPS(idx) will return a 1x6 cell array of name-value pairs that
%   specify a marker, linestyle, and color combination identifies by idx
%   which is a positive integer.  These name-value pairs can be applied 
%   to to the inputs of many plotting function using, 
%       plot(___, nameval{:}); 
%
%   If you're only varying linestyle and color, set the LineStyleOrder 
%   and ColorOrder properties of the axes, though be aware of some 
%   compatibility issues prior to Matlab r2019b (see footnotes [1,2]).  
%
%   LINEPROPS('count') returns the number of property combinations. 
%   combinations
%
%   Examples
%   ------------
%   figure()
%   hold on
%   for i = 1:10
%       nameval = lineprops(i); 
%       plot(0:.1:1,rand(1,11)+i, nameval{:}, 'LineWidth', 1.5)
%   end
%   legend('Location','BestOutside')
%
% lineprops() is available at
% https://www.mathworks.com/matlabcentral/answers/697655#answer_579560
% Author: Adam Danz
% Copyright (c) 2020  All rights reserved

% Version history
% vs 1.0 18-Dec-2020 created and uploaded to https://www.mathworks.com/matlabcentral/answers/697655#answer_579560
% vs 1.1 18-Dec-2020 Added URL source to help section

%% Input validation
narginchk(1,1)
if (ischar(idx) || isa(idx,'string')) && all(strcmpi(idx,'count'))
    returnCount = true;
else
    validateattributes(idx,{'numeric'},{'scalar','numel',1})
    returnCount = false;
end

%% line properties
% List a bunch of markers; they will be selected in
% order and then the selection will start again if
% there are more lines than markers.
markers = {'o','+','*','x','s','d','v','^','>','p','v','<','h'};

% List a bunch of colors; like the markers, they
% will be selected circularly.
%colors = {'r' 'g' 'b' 'c' 'm' 'k'}; % not use 'w' or 'y' due to visibility on white axes
colors    = {'k', 'r', 'b', [0, 0.5, 0], [0.75, 0, 0.95], [0.9290, 0.6940, 0.1250], [0.3010 0.7450 0.9330], 'm', 'g', 'c'};

% Same with line styles
%linestyles = {'-','--',':','-.'};
linestyles = {'-','--','-.'};

if returnCount
    nameval = lcm(numel(markers),lcm(numel(colors),numel(linestyles)));
else
    first = @(v)v{1};
    prop = @(options, idx)first(circshift(options,-idx+1));
    nameval = {'Marker', prop(markers,idx), 'LineStyle', prop(linestyles,idx), 'Color', prop(colors,idx)};
end

%% Footnotes
% [1] https://www.mathworks.com/help/matlab/creating_plots/defining-the-color-of-lines-for-plotting.html#mw_bfc11d57-558f-4191-a39f-b19be8d2e4bd
% [2] https://www.mathworks.com/help/matlab/ref/matlab.graphics.axis.axes-properties.html#budumk7_sep_shared-mw_4ec3600f-1440-4d69-aeff-313e363abcda    
