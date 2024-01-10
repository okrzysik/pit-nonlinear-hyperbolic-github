function [weights, nodes, error] = poly_interpolation_odd_vector(epsilon, interp_order)

% Basically a vectorized implementation of polynomial interpolation, but it
% holds only for odd-degree interpolation polynomials (the stencil for
% these is independent of epsilon, but for even-degree polynomials the
% stencil shifts based on epsilon, so it makes the calculation a bit hard
% to vectorize).

p = interp_order;

assert(mod(interp_order, 2) == 1, 'interpolation holds only for odd degree polynomials')


% Left-most and right-most nodes are set centrally about right neighbour of
% interpolation point.
xL = -(p+1)/2;
xR = (p-1)/2;

nodes = (xL : xR)';    
nx    = numel(epsilon);
error = ones(nx, 1);


% Compute the `error` coefficient. This appears in every coefficient when
% we divide out the corresponding value of (node(i) + epsilon)
for i = 1:numel(nodes)
    error = error .* (nodes(i) + epsilon);
end

% Divide out this value.
weights = error ./ (nodes' + epsilon);

for i = 1:numel(nodes)
    temp = 1;
    for j = 1:numel(nodes)
        if j ~= i
            temp = temp * (nodes(j) - nodes(i));
        end
    end

    % Do one division at end since this is faster than many divisions
    weights(:, i) = weights(:, i) / temp;
end

        % The general formula breaks below when the interpolation point coincides 
% with an interpolation node (when epsilon == 0 it is its right neighbour, 
% and when epsilon = 1 it is its left neighbour). So manually fix-up these
% cases
ii = find(epsilon == 0);
weights(ii, :) = 0;
weights(ii, find(nodes == 0)) = 1;
error(ii) = 0;

ii = find(epsilon == 1);
weights(ii, :) = 0;
weights(ii, find(nodes == -1)) = 1;
error(ii) = 0;


end