% The exact solution to Burgers equation u_t + (u^2/2)_x = 0 when the
% initial data is u(x, 0) = 1 for -0.5 < x < 0, and zero otherwise in (-1,
% 1), and periodic BCs are used.
%
% The solution is valid only up to time t = 4, since it doesn't account for
% the shockwave running into the base of the rarefaction wave.
%
% NOTE: This is the exact point-wise solution, not the cell-averaged
% solution.

function u = B3_exact_sol_pt_wise(x, t)
    assert(numel(t) == 1, 'Solution accepts only scalar t')
    assert(t <= 4, 'Solution not implemented for t > 4')
    
    u = zeros(size(x));
    for i = 1:numel(x)
        u(i) = exact_sol(x(i), t);
    end
end


function u0 = initial_condition(x)
    u0 = zeros(size(x));
    I = find(x > -1/2 & x < 0);
    u0(I) = 1;
end

% Note: This function is not vectorized, it can only accept scalar x and t.
function u = exact_sol(x, t)
    % Ensure periodicity by mapping (-1,-1/2) into (1,3/2). 
    if x < -0.5; x = x + 2; end 

    % t = 0
    if t == 0
        u = initial_condition(x);
    % 0 < t < 1
    elseif t < 1
        xs  = shock_position(t);
        xf = front_position(t);
        
        if x < -1/2 
            u = 0;
        elseif x < xf
            u = (x + 0.5)/t;
        elseif x < xs
            u = 1;
        else
            u = 0;
        end
        
    % t >= 1
    else
        xs  = shock_position(t);
        
        if x < -1/2 
            u = 0;
        elseif x < xs
            u = (x + 0.5)/t;
        else
            u = 0;
        end
    end
end

% Front of rarefaction wave
function xf = front_position(t)
    if t < 1
        xf = t - 0.5;
    else
        xf = shock_position(t);
    end
end

% Shock position
function xs = shock_position(t)
    if t == 0
        xs = 0;
    elseif t < 1
        xs = t/2;
    else
        xs = sqrt(t) - 0.5;
    end
end