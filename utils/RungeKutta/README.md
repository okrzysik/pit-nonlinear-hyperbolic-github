Here live some RK solvers for linear/nonlinear ODE systems of the form
	du/dt = f(t,u), u(0) == u0.

At the moment, I currently have the following implemented
1. ERK schemes, where the user specifies the function f.
2. DIRK schemes, where the user has to pass a solver for inverting the nonlinear system that's required to be solved at each time step.

The idea is that the user can just pass a Butcher table for their favourite method, and then use it to solve the system. See the Butcher table m file in the current directory.

The solution is computed at fixed times specified by the user, i.e., adaptive time-stepping is not used (I really don't want this capability for the problems I'm solving anyways).

There are also a couple of scripts here that compare the solutions generated with MATLAB's ODE solvers, just to check they're converging at the correct rate.
