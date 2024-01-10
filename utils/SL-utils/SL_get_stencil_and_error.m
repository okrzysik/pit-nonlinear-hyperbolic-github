% Compute the stencil and error associated with a single step. This strips
% out the (main) part of the SL step that we can save and store. 
function [d, nodes, interp_error] = SL_get_stencil_and_error(depart_east_epsilon, nx, k)

    % A is an upper bidiagonal difference matrix. 
    e = ones(k, 1);
    A = diag(-e, 0) + diag(e(1:end-1), 1);

    [weights_primitive, nodes_primitive, interp_error] = poly_interpolation_odd_vector(depart_east_epsilon, k);
    weights_primitive = weights_primitive';
    nodes = nodes_primitive(2:end);

    chat = weights_primitive(2:end, :);
    if k == 1
        chat(1, :) = chat(1, :) - 1;
    elseif k == 3
        chat(2, :) = chat(2, :) - 1;
    elseif k == 5
        chat(3, :) = chat(3, :) - 1;
    else
        error('Coefficients not set up for other k...') % Note that for even k we'd need to actually check where the zero node is.
    end

    % Solve for interplolation weights for all points.
    d = A \ chat;

    d = d'; % Re-orient for corrected usage.
end