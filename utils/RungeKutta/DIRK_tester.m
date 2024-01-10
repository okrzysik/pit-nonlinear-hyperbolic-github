% Test DIRK implementation is correct, by comparing numerical solutions with
% exact solutions.
%
% Seems like all schemes are converging with the right order!

clc
clear
close all

% Various RHS of ODE systems. This ODE system is stiff. It seems like
% adaptive time stepping probably needed...
myodes = @(t, u) [u(2); 10*(1-u(1)^2)*u(2)-u(1)]; % vdp1000, as in ode15s examples

mytol = 1e-12; % Ramp up this tol to get 'exact solution'
options = odeset('RelTol', mytol, 'AbsTol', mytol);

tspan = [0; 4];
y0 = [2; 0];   
[t,y] = ode15s(myodes, tspan, y0, options);

plot(t,y(:,1), '-ro');
hold on
plot(t,y(:,2), '-rx');


%butcher = butcher_table('BE');
butcher = butcher_table('SDIRK2');
%butcher = butcher_table('SDIRK3');

% Set up nonlinear solver for computing stage vectors. See DIRK_solver.m
A = butcher.A;
options = optimset('TolFun', 1e-12, 'Display', 'None');
H = @(xi, dt, i, q, p) p - dt*A(i,i)*myodes(xi, p) - q; % We want to solve the equation H(p) == 0.
mysolver = @(xi, dt, i, q, p0) fsolve(@(p) H(xi, dt, i, q, p), p0, options); % Here is the solver to solve H(p) == 0.


u = DIRK_solver(mysolver, t, y0, butcher, 'all');
plot(t, u(1, :), '-bo')
plot(t, u(2, :), '-bx')

% % Just get the solution at the end of the time interval and compare this
% % with what MATLAB's ODE solvers get. If we make them solve the ODEs really
% % accurately, then it's like comparing against exact solution.
% u = DIRK_solver(mysolver, t, y0, butcher, 'final');
% plot(t(end), u(1), '-bo')
% plot(t(end), u(2), '-bx')
% 
% Refine time step by factor of half, check convergence rate. Should be 2^p
dt = 0.1;
for k = 2:6
    nt = tspan(2)/dt;
    t = linspace(0, tspan(2), nt);
    u = DIRK_solver(mysolver, t, y0, butcher, 'final');
    error(k) = norm(u - y(end, :)', inf);
    dt = dt/2;
end
[error(1:end-1)./error(2:end)]'
