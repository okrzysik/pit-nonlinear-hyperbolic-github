function u = ERK_solver(odefun, t, u0, butcher_table, storage)
%ERK_SOLVER apply ERK to the system du/dt = odefun(t,u), u(0) = u0.
%
%INPUT:
%   odefun          :   HANDLE. odefun(t, u) must return the RHS of the ODEs evaulated at time t
%                               and solution u. 
%                               dt is the size of the current time step. 
%                               i is between 1 and s, and is the index of the current stage being solved for. 
%   t               :   ARRAY. Discrete times to approximate the solution at. This must
%                               include the initial time where u == u0.
%   u0              :   ARRAY. Initial condition.
%   butcher_table   :   STRUCT. The Butcher table of the scheme. Must have
%                               the three fields 'A', 'b', 'c'.
%   storage         :   STRING. How much of the solution to store
%                           'final'     :   u contains the solution at t(end) only.
%                           'all'       :   u constains the solution at all times in t (including t(1)).
%
%OUTPUT:
%   u   :   ARRAY. ERK solution of the system. The amount of the solution contained in u 
%                       dependeds on 'storage_option.'
%
%NOTES:

if nargin(odefun) > 2
    error('odefun should only have 2 arguments; an old version accepted 4')
end

m = size(u0, 1); % Dimension of solution.

% Unapck Butcher table.
A = butcher_table.A;
b = butcher_table.b;
c = butcher_table.c;

if sum(sum(abs(triu(A)))) ~= 0
    error('Stricly lower trinagular Butcher table necessary for ERK')
end

s = size(A, 1); % Number of stages
k = zeros(m, s); % Allocate space for stage vectors (recycled at every step).
nt = numel(t)-1; % Number of time steps to take.

% Allocate storage for solution at all times t.
if strcmp(storage, 'all')
    u = zeros(m, nt+1);
    u(:, 1) = u0; % Include initial solution.
else
    if ~strcmp(storage, 'final')
        error('solution storage option ''%s'' is invalid', storage)
    end
end


% Time step!
step = 1; % == number of step we're taking
while step <= nt
    dt = t(step+1) - t(step); % Time step to take (difference between where we are and where we're headed).
    
    u1 = u0; % Initialize solution we're targeting.
    % Compute the s stage vectors k_i    
    for i = 1:s
        k(:, i) = u0; % Initialize stage vector.
        xi = t(step) + dt*c(i); % Node that the stage is evaluated at.
        for j = 1:i-1
            if A(i,j) ~= 0
                k(:, i) = k(:, i) + dt*A(i,j)*k(:, j);
            end
        end
        k(:, i) = odefun(xi, k(:, i)); % Evaluate ODEs at this point to get stage vector.
    end

    % Build solution at new time (remembering we initialized with u1 = u0).
    for i = 1:s
        if b(i) ~= 0
            u1 = u1 + dt*b(i)*k(:, i);
        end
    end
    
    % Add current solution into global solution matrix.
    if strcmp(storage, 'all')
        u(:, step+1) = u1;
    end
    
    % Increment counter.
    step = step + 1;  
    if step <= nt
        u0 = u1; % Initial condition for next step.
    end
end

% Return the solution at final time only.
if strcmp(storage, 'final')
    u = u1;
end

end