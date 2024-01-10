function P = linear_interpolation_non_nested(t_coarse, t_fine)
%
% The individual grids are structured with equal spacing, but they are such
% that fine time points don't necessarily coincide with coarse time points.
% I.e., the coarse grid is not a subset of the fine time grid (they are
% non-nested).
% The below implementation makes an assumption that there's at least 3
% fine-grid time points in the inbetween two coarse-grid time points. 
% (with the possible exception that the fine grid extends out past the
% end of the coarse grid, and these points just use constant interpolation
% to the last coarse point rather than linear).


nt_coarse = numel(t_coarse);
nt_fine   = numel(t_fine);

P = sparse(nt_fine, nt_coarse);

a = 1; % Index of coarse time point that is below current fine time point
b = 2; % Index of coarse time point that is above current fine time point

% The solution at current fine time point t is estimated by linearly
% interpolating to t_coarse(a) and t_coarse(b).

% Sequentially step through fine time points since this makes it's easy to
% find the two nearest coarse time points that are next to us.
for n = 1:nt_fine
    
    t_current = t_fine(n);

    ta = t_coarse(a);
    tb = t_coarse(b);
    
    % Is t_current actually in [t_coarse(a), t_coarse(b)]? If not then we 
    % increment to the next interval [t_coarse(a+1), t_coarse(b+1)] which
    % it must be in due to the structuring of the grid (every coarse-grid
    % interval [t_coarse(a), t_coarse(b)] contains at least 3 fine-grid
    % time points in it).
    if t_current > tb
        % If not at the last interval in the coarse grid then increment to
        % the next interval. If the fine-grid point extends past the last 
        % coarse interval then don't update a or b, and this will result in
        % the below interpolation formula being extrapolated from the last
        % coarse interval (this seems more accurate than just constant
        % interpolation).
        if b < nt_coarse
            b = b + 1; tb = t_coarse(b);
            a = a + 1; ta = t_coarse(a);
        end
    end
        
    wa = (t_current - tb)/(ta - tb); % Weight for ta
    wb = (t_current - ta)/(tb - ta); % Weight for tb

    P(n, a) = wa;
    P(n, b) = wb;
end