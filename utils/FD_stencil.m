function [nodes, weights] = FD_stencil(derivative, bias, approx_order)
%[nodes, weights] = FD_stencil(derivative, bias, approx_order)
%
%Reads in the FD nodes and weights from text files in the current directory
%and returns them as MATLAB arrays.
%
%INPUT:
%   derivative:     INT. The derivative to approximate
%   bias:           STRING/CHAR. The bias of the stencil. 'upwind', or 'U', or
%                           'central', or 'C'
%   approx_order:   INT. The order of the approximation.


filename = strcat('D', num2str(derivative), '_', upper(bias(1)), num2str(approx_order), '.txt');

fileID = fopen(filename, 'r');
if fileID == -1
    error('The file ''%s'' does not exist', filename);
end

n = fscanf(fileID, '%d', 1);
nodes = fscanf(fileID, '%d', n);
weights = fscanf(fileID, '%f', n);

fclose(fileID);

end