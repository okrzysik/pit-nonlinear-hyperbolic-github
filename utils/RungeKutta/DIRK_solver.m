function u = DIRK_solver(odesolver, t, u0, butcher_table, storage)
%DIRK_SOLVER apply DIRK to the system du/dt = odefun(t,u), u(0) = u0.
%
%INPUT:
%   odesolver       :   HANDLE. odesolver(xi, dt, i, q, p0) must return the solution p of
%                               the nonlinear system: 
%                                   p - dt*A(i,i)*odefun(xi, p) - q == 0,
%                               where: 
%                                       1. xi is the current time (xi == t + c(i)*dt),
%                                       2. dt is the current time step,
%                                       3. i is the index of the stage we're in the process of computing,
%                                       4. A(i,i) is the ith diagonal entry from the Butcher table matrix A,
%                                       5. p0 is a good initial guess at p (or so I think...).
%                                       
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

m = size(u0, 1); % Dimension of solution.

% Unapck Butcher table.
A = butcher_table.A;
b = butcher_table.b;
c = butcher_table.c;

if sum(sum(abs(triu(A, 1)))) ~= 0
    error('Lower trinagular Butcher table necessary for DIRK')
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
        q = u0; % Initialize stage vector.
        xi = t(step) + dt*c(i); % Node that the stage is evaluated at.
        for j = 1:i-1
            if A(i,j) ~= 0
                q = q + dt*A(i,j)*k(:, j);
            end
        end
        p = odesolver(xi, dt, i, q, q); % Solve the m-dimensional nonlinear problem. Use q as initial guess for p.
        k(:, i) = (p - q)/(dt*A(i,i));
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