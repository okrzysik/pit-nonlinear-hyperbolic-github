% Test ERK implementation is correct, by comparing numerical solutions with
% exact solutions.
%
% Seems like all schemes are converging with the right order!

clc
clear
close all

% Various RHS of ODE systems.
myodes = @(t, u) [cos(4*t); sin(t)];
%myodes = @(t, u) vdp1(t, u); % As in the ODE 45 solver examples
%myodes = @(t, u) vdp1(t, u) + [cos(t); sin(t)];


mytol = 1e-14; % Ramp up this tol to get 'exact solution'
options = odeset('RelTol', mytol, 'AbsTol', mytol);

tspan = [0; 20];
y0 = [0.8; 0.4];   
[t,y]=ode45(myodes, tspan, y0, options);

plot(t,y(:,1), '-ro');
hold on
plot(t,y(:,2), '-rx');


%butcher = butcher_table('FE');
%butcher = butcher_table('Ralston2');
%butcher = butcher_table('Heun2');
%butcher = butcher_table('Heun3');
butcher = butcher_table('SSPERK3');
%butcher = butcher_table('RK4classic');
%butcher = butcher_table('RK5');
%butcher = butcher_table('RK6');

u = ERK_solver(myodes, t, y0, butcher, 'all');
plot(t, u(1, :), '-bo')
plot(t, u(2, :), '-bx')

% Just get the solution at the end of the time interval and compare this
% with what MATLAB's ODE solvers get. If we make them solve the ODEs really
% accurately, then it's like comparing against exact solution.

u = ERK_solver(myodes, t, y0, butcher, 'final');
plot(t(end), u(1), '-bo')
plot(t(end), u(2), '-bx')

% Refine time step by factor of half, check convergence rate. Should be 2^p
dt = 0.5;
for k = 1:8
    nt = tspan(2)/dt;
    t = linspace(0, tspan(2), nt);
    u = ERK_solver(myodes, t, y0, butcher, 'final');
    u
    y(end, :)'
    error(k) = norm(u - y(end, :)', inf);
    dt = dt/2;
end
[error(1:end-1)./error(2:end)]'
