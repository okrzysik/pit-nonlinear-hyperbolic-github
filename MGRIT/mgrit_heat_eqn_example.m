% Apply MGRIT to solve the 1D heat equation: 
%    u_t=eta*u_xx
% An implicit Euler discretization is used in time, and second-order finite
% differences are used in space. 

clc
clear

eta = 0.025;        % Diffusion coefficient.
nx = 100;           % Number of spatial points
h  = 1/nx;          % Spatial step size
dt = 10*h^2;        % Time step size
nt = 1/dt;          % Number of time points 
nt = round(nt);
m  = 8;             % Coarsening factor

x = linspace(-1,1,nx); 
t = linspace(0,1,nt);
u0 = sin(pi*x);       % Initial condition

% Package everything into a struct to pass to mgrit solver
myMGRIT_object = struct();

% Tell mgrit solver what the dimension of Phi is.
myMGRIT_object.block_size = nx;

% Need to give mgrit solver the fine time grid.
myMGRIT_object.t = t;

% Get default solver parameters (see function below).
solver_params = my_heat_eqn_mgrit_solver_params();
% Now update the solver_params struct with the parameters we want to alter from the default values
solver_params.cf = m;

% Set up fine-level linear system
% Assemble RHS of global fine-level linear system. This only holds initial condition at t=0.
g       = zeros(nt*nx, 1);
g(1:nx) = u0; % Update first block row of g.
% Solution vector
u       = rand(nx*nt, 1);
u(1:nx) = u0;

% MGRIT Solve
tic
[u, rnorm, myMGRIT_object] = mgrit(u, g, @(a, b, c) my_heat_eqn_step(a, b, c, eta), myMGRIT_object, solver_params);
fprintf('MGRIT timer = %.2fsecs\n', toc)

%% Residuals plot
figure
semilogy((0:numel(rnorm)-1), rnorm/rnorm(1), ...
    'g>-', 'LineWidth', 2)
fs = {'Interpreter', 'latex', 'FontSize', 24};
xlabel('$k$', fs{:})
ylabel('$\Vert \mathbf{r}_k \Vert / \Vert \mathbf{r}_0 \Vert$', fs{:})
title('MGRIT residual history')

%% Solution plot
figure
[X, T] = meshgrid(x, t);
U = reshape(u, [nx nt]).';
mesh(X, T, U)
box on
xlabel('$x$')
ylabel('$t$')
title('Heat eqn. sol. via MGRIT')


%% Helper functions

% This is a function that has to evolve the solution u at time t0 into the
% solution at time t1. Information about the current level, the time-step
% size, etc. is saved in step_status. MGRIT_object holds a bunch of
% information that the solver requires. Notice that MGRIT_object is also
% returned by this function, so you can store some things in it (which is
% what we do here with the handle for applying the LU factorization of
% Phi).
function [u1, MGRIT_object] = my_heat_eqn_step(u0, step_status, MGRIT_object, eta)
    % Extract basic step information from step_status
    level  = step_status.level;
    dt     = step_status.dt;
    
    % Compute the LU factorization of the matrix once and get handle for 
    % applying it. We then store this LU factorization handle in the 
    % MGRIT_object so that we can use it whenever this my_stepper function 
    % is called again.
    % Check if LU has already been computed, and if not then compute it
    if ~isfield(MGRIT_object.hierarchy(level), 'step_handle') || ...
            isempty(MGRIT_object.hierarchy(level).step_handle)
        
        nx = MGRIT_object.block_size;
        
        A = laplacian_1d(nx);
        B = speye(nx) - dt*eta*A;
        [L, U] = lu(B);
        
        % Store the handle in the MGRIT_object.hierarchy.
        MGRIT_object.hierarchy(level).step_handle = @(b) U \ (L \ b);
    end
    
    % Apply Phi to u0 to get u1
    u1 = MGRIT_object.hierarchy(level).step_handle(u0);
end

function A = laplacian_1d(nx)
    e = ones(nx, 1);
    A = spdiags([e -2*e e], -1:1, nx, nx) * nx^2;
end


function solver_params = my_heat_eqn_mgrit_solver_params()
    solver_params                 = struct();
    solver_params.cf              = 8;
    solver_params.res_halt_tol    = 1e-10; 
    solver_params.res_reduction   = 1; 
    solver_params.maxlevels       = 10; 
    solver_params.maxiter         = 10; 
    solver_params.min_coarse_nt   = 2; 
    solver_params.pre_relax       = 'FCF';
    solver_params.final_F_relax   = false;
end