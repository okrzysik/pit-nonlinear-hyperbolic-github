% The exact solution to Burgers equation u_t + (u^2/2)_x = 0 when the
% initial data is u(x, 0) = 1 for -0.5 < x < 0, and zero otherwise in (-1,
% 1), and periodic BCs are used.
%
% The solution is valid only up to time t = 4, since it doesn't account for
% the shockwave running into the base of the rarefaction wave.
%
% NOTE: This is the exact cell-averaged solution, not the point-wise
% solution. The cell-averaged solution is computed from integrating the
% exact point-wise solution over a given FV cell. This is a tedious task
% because there are many cases to consider, including ones where the FV
% cell is split over an interval where the point-wise solution changes
% piecewise. For a fully smooth function one would normally just
% compute the cell averages numerically, but in this that doesn't really
% give accurate results due to the non-smoothness of the point-wise
% solution.
%
% NOTE: The results here probably only make sense in cases where h < 1/2.

% For a given time t, the vector of cell averages is evaluated, where the
% cells are centered on x_centers, and are assumed to be of equal width.
function u_bar = B3_exact_sol_cell_avg(x_centers, t)
    assert(numel(t) == 1, 'Solution accepts only scalar t')
    assert(t <= 4, 'Solution not implemented for t > 4')
    
    x = x_centers;
    nx = numel(x_centers);
    u_bar = zeros(size(x));
    h = x(2) - x(1);
    
    % Ensure periodicity by mapping (-1,-1/2) into (1,3/2). 
    I = find(x < -0.5);
    x(I) = x(I) + 2;
    %if x < -0.5; x = x + 2; end 
    
    % t = 0
    if t == 0
        %u = initial_condition(xl, xr, h);
        for i = 1:nx
           u_bar(i) = initial_condition(x(i) - h/2, x(i) + h/2, h);
        end
        
        
    % 0 < t < 1
    elseif t < 1
        xs  = shock_position(t);
        xf = front_position(t);
        
        %u_bar = pre_merge(xl, xr, h, t, xs, xf);
        for i = 1:nx
           u_bar(i) = pre_merge(x(i) - h/2, x(i) + h/2, h, t, xs, xf);
        end
        
    % t >= 1
    else
        
        xs  = shock_position(t);
        %u_bar = post_merge(xl, xr, h, t, xs);
        for i = 1:nx
           u_bar(i) = post_merge(x(i) - h/2, x(i) + h/2, h, t, xs);
        end
        
    end
end



function u_bar = initial_condition(xl, xr, h)
    
    % Case 1
    if xr <= -0.5
        u_bar = 0;
        
    % Case 2
    elseif xl <= -0.5 && xr > -0.5
        u_bar = (xr + 1/2)/h;
        
    % Case 3
    elseif xl > -0.5 && xr < 0
        u_bar = 1;
        
    % Case 4
    elseif xl < 0 && xr >= 0
        u_bar = -xl/h;
       
    % Case 5
    elseif xr > 0
        u_bar = 0;
    
    else
        error('I missed a case...')
    end

end

function u_bar = pre_merge(xl, xr, h, t, xs, xf)

    % Case 1
    if xr <= - 0.5
        u_bar = 0;
        
        
    % Case 2a: xl is LHS of base of rarefaction, and RHS is in rarefaction
    elseif xl < -0.5 && -0.5 < x_r && xr <= xf
        u_bar = ( xr*(xr+1) + 0.25 ) / h / (2*t);
        
        
    % Case 2b: xl is LHS of base of rarefaction, and RHS is past the rarefaction front
    elseif xl < -0.5 && xr > xf && xr <= xs
        u_bar = ( xf*(xf+1) + 0.25 ) / h / (2*t) + (xr - xf)/h;
        
    % Case 3
    elseif xr > -0.5 && xr <= xf
        u_bar = ( xr*(xr+1) - xl*(xl+1) ) / h / (2*t); 
        
    % Case 4a
    elseif xl <= xf && xr <= xs
        u_bar = ( xf*(xf+1) - xl*(xl+1) ) / h / (2*t) + (xr - xf)/h; 
       
    % Case 4b
    elseif xl <= xf && xr > xs
        u_bar = ( xf*(xf+1) - xl*(xl+1) ) / h / (2*t) + (xs - xf)/h; 
        
    % Case 5a
    elseif xl > xf && xr < xs
        u_bar = 1;
        
    % Case 5b
    elseif xf < xl && xl < xs && xr >= xs
        u_bar = (xs - xl)/h;
        
    % Case 6
    elseif xl >= xs
        u_bar = 0;

    else
        
%         [xl, xr, xs, xf]
%         x = [xl; xr];
%         plot(x, B3_exact_sol(x, t), '-ro')
%         hold on
%         x = linspace(-1, 1, 1e3);
%         plot(x, B3_exact_sol(x, t), 'k')
        
        error('I missed a case...')
    end

end



function u_bar = post_merge(xl, xr, h, t, xs)

    % Case 1
    if xr <= - 0.5
        u_bar = 0;
        
    % Case 2
    elseif xl < -0.5 && xr >= -0.5 
        u_bar = ( xr*(xr+1) + 0.25 ) / h / (2*t);
        
    % Case 3
    elseif xr > -0.5 && xr <= xs
        u_bar = ( xr*(xr+1) - xl*(xl+1) ) / h / (2*t); 
        
    % Case 4
    elseif xl < xs && xr > xs
        u_bar = ( xs*(xs+1) - xl*(xl+1) ) / h / (2*t); 

    % Case 5
    elseif xl >= xs
        u_bar = 0;
        
    else
        
        error('I missed a case...')
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