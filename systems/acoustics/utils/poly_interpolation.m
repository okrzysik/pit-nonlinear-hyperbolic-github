function [weights, nodes, error] = poly_interpolation(epsilon, interp_order)
% interp_order is the degree of the interpolation polynomial.
%
% epsilon is the distance from the interpolation point to its right neighbour. It should
% always be in the interval [0, 1).
%
% The p+1 nodes are chosen such that they are closest nodes to the
% interpolation point. 
%
% The node immediately to the right of the interpolation point has index 0.
%
% error is the error coefficient arising from the interpolation, but it
% isn't divided by the (interp_order+1)! term.

p = interp_order;
% Polynomial is even degree, use an odd number of nodes
if mod(interp_order, 2) == 0
    % Left-most and right-most nodes are set centrally about right neighbour of
    % interpolation point.
    xL = -p/2;
    xR = p/2;
    % If interpolation point is closer to left neighbour than right neighbour,
    % then shift stencil one to the left.
    if epsilon > 0.5
        xL = xL - 1;
        xR = xR - 1;
    end
% Polynomial is odd degree, use an even number of nodes
else
    % Left-most and right-most nodes are set centrally about right neighbour of
    % interpolation point.
    xL = -(p+1)/2;
    xR = (p-1)/2;
end

nodes = (xL : xR)';
    

% The general formula breaks below when the interpolation point coincides 
% with an interpolation node (when epsilon == 0 it is its right neighbour, 
% and when epsilon = 1 it is its left neighbour). So manually fix-up these
% cases
if epsilon == 0
    weights = zeros(size(nodes));
    weights(find(nodes == 0)) = 1;
    error = 0;
elseif epsilon == 1
    weights = zeros(size(nodes));
    weights(find(nodes == -1)) = 1;
    error = 0;

% General formula when interpolation point does not coincide with an
% interpolation node.
else

    % Compute the `error` coefficient. This appears in every coefficient when
    % we divide out the corresponding value of (node(i) + epsilon)
    error = 1;
    for i = 1:numel(nodes)
        error = error*(nodes(i) + epsilon);
    end

    % Divide out this value.
    weights = error ./ (nodes + epsilon);

    for i = 1:numel(nodes)
        temp = 1;
        for j = 1:numel(nodes)
            if j ~= i
                temp = temp * (nodes(j) - nodes(i));
            end
        end
        % Do one division at end since this is faster than many divisions
        weights(i) = weights(i) / temp;
    end

end

%%%%%%%%%%%%%% Older implementation %%%%%%%%%%%%%%
%%% This implementation uses internal MATLAB commands to execute products,
%%% but it's actually slower upon timing. Maybe because there are so few
%%% numbers in the products that the overhead dominates. 
% error = prod(nodes + epsilon);
% weights = error ./ (nodes + epsilon);
%
% for i = 1:numel(nodes)
%     temp = nodes(i);
%     nodes(i) = nodes(i) + 1;
%     weights(i) = weights(i) / prod( nodes - temp );
%     nodes(i) = nodes(i) - 1;
% end

%%%%%%%%%%%%%% Older implementation %%%%%%%%%%%%%%
% weights = ones(p + 1, 1);
% % Compute Lagrangian interpolation polynomial coefficients
% for i = 1:interp_order+1
%     j = xL+i-1;
%     for m = xL:xR
%         if (m ~= j)
%             weights(i) = weights(i) * (m + epsilon)/(m - j);
%         end
%     end
% end
% 
% % Compute error coefficient
% if nargout > 2
%     error_coefficient = 1;
%     for m = xL:xR
%         error_coefficient = error_coefficient * (m + epsilon);
%     end
% end

end