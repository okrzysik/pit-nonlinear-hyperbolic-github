% Create time mesh. This will ensure dt is set according to the specified 
% CFL number, but it will do so at the expense of not necessarily stopping 
% at tmax_approx (since we require an integer number of time steps).
function [t, dt] = time_mesh(h, CFL_number, abs_fprime_max, spatial_order, tmax_approx)
    
    % Set time step according to CFL number and ensure we maintain a
    % temporal accuracy that matches the spatial accuracy.
    if spatial_order == 1 || spatial_order == 3
        dt = CFL_number * h / abs_fprime_max; 

    elseif spatial_order == 5
        dt = CFL_number * h^(5/3) / abs_fprime_max; 
    end

    nt   = floor(tmax_approx/dt)+1; % Ensure integer number of time steps (from dt = tmax/(nt+1))
    % nt = 2^floor(log2(nt)); % Ensure number of time steps is a power of 2, but maintain the same CFL number.
    tmax_true = dt*nt;

    % Temporal grid
    t  = linspace(0, tmax_true, nt+1)'; % There are nt+1 points in between 0 and tmax_true; 
  
end