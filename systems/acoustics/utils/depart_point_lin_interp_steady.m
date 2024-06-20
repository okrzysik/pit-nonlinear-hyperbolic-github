% UPDATED: This computes departure points when the wave-speed is
% independent of time, so you can just pass the departure points for a
% single fine-grid step. I've commented out below the couple of lines where
% in the original code you'd extract the time-dependent fine-grid departure
% points (and I simplified the input too).
%
function depart_points_coarse = depart_point_lin_interp_steady(m, x, depart_points_fine)
%
% Uses a linear interpolation strategy to advance a coarse-grid
% characterstic from the arrival time to the departure time.
%
% The diagram below shows how this works on the first coarse level. Its the
% same on coarser levels, just with dt -> m^(level-1)*dt
%
% C, k = m,    t = t_n + m \delta t
% F, k = m-1,  t = t_n + (m-1) \delta t
% F, k = m-2,  t = t_n + (m-2) \delta t
% .
% .
% .
% F, k = 2,   t = t_n + \delta t
% C, k = 1,   t = t_n


% Some spatial domain knowledge we need.
nx   = numel(x);
dx   = x(2) - x(1);
xmin = x(1);
xmax = x(end) + dx;
domain_length = xmax - xmin;



% At the last F-point in CF-interval, the coarse- and fine-grid 
% characteristics are the same, i.e., they intersect the x-axis at the same 
% location location.
%k = m-1;
%depart_points_coarse = MGRIT_object.hierarchy(level-1).depart_points{t0idx_fine + k};
depart_points_coarse = depart_points_fine;

% Step backwards m-1 times to get to evolve the coarse-grid cgaracteristic 
% from the last F-point in the CF-interval to the C-point at the beginning 
% of the interval.
%
% We're estimating the coare-grid departure point at time t_n + (k-1)*dt
% using the nearest neighbours of the fine-grid departure points at that
% time on either side of us.
for k = m-2:-1:0    

    % Departure points for stepping from the jth F-point to its
    % right neighbour
    %depart_points_fine = MGRIT_object.hierarchy(level-1).depart_points{t0idx_fine + k};
    % depart_points_fine = depart_points_fine;

    % Index of mesh point to right of the current value of the 
    % coarse-grid characteristic.
    % NOTE: This is not necessarily on the grid since the
    % characteristic is free to exit the boundaries.
    % An index in [1,nx] is on the mesh, otherwise it exceeds it.
    rn_idx = ceil((depart_points_coarse - xmin) / dx)+1; 

    % Distance the characteristic is from the mesh point to its right
    rn_epsilon = depart_points_coarse - dx * (rn_idx-1) - xmin;

    %     Periodically map the index onto the set of index points for the 
    %     physical mesh. 
    %     That is, we want to translate q by some interger multiple of nx
    %     such that it ends up in the interval [1,nx].
    %     Decompose q as
    %         q = a*nx + b, 
    %     with abs(b) < nx, and for a,b integers. We have that 
    %     a = floor(q/nx), b = mod(q,nx).    
    %     Then define
    %         p == b      = q - a*nx,      if b > 0,
    %         p == b + nx = q - a*nx + nx, if b <= 0,
    %     clearly p is the desired translate of q.                
    b = mod(rn_idx, nx);
    b(b <= 0) = b(b <= 0) + nx;
    rn_idx_physical = b;
        
    % The below is how I did the above previously, but the below makes more
    % sense, I think.
    %rn_idx_physical = mod(rn_idx + nx-1, nx)+1;
%     temp = norm(mod(rn_idx + nx-1, nx)+1 - b);
%     if temp > 0
%         [mod(rn_idx + nx-1, nx)+1, b]
%         error('oops')
%     end
    
    

    % By how many lengths of the domain did the characteristic get 
    % shifted when mapping the index.
    period_shifts = (rn_idx - rn_idx_physical)/nx;

    % Loop over all coare-grid characteristics, extending them by an
    % amount dt using linear interpolation of the fine-grid
    % characteristics on either side of the starting point.

    % Compute the difference between the arrival points of the
    % fine-grid characteristics that arrive at either side of me. 
    % The index of the departure point to the left is always 1 less
    % than the one to the right, unless my right neighbour is the
    % first mesh point, then by periodicity my left neighbour has
    % index nx. AND, in  this case need to account for the fact
    % that there fine depart point is offset the length of the
    % domain.

    
    % Find any troublesome indices.
    I = find(rn_idx_physical == 1);
    temp = rn_idx_physical(I);
    rn_idx_physical(I) = 2; % Set to a dummy value to allow vectorized calc.
    
    % Do vectorized calculation for all points.
    depart_fine_difference = depart_points_fine(rn_idx_physical) - depart_points_fine(rn_idx_physical-1);
    
    % Restore true values.
    rn_idx_physical(I) = temp;
    % Do correct calculation for the few troublesome points.
    depart_fine_difference(I) = depart_points_fine(rn_idx_physical(I)) - depart_points_fine(nx) + domain_length;

    % Do linear interpolation
    depart_points_coarse = depart_fine_difference/dx .* rn_epsilon + ...
                                             depart_points_fine(rn_idx_physical) + period_shifts*domain_length;
end

end