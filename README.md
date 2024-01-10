This project is about the parallel-in-time solution of (potentially nonlinear) hyperbolic PDEs. Everything is limited to one spatial dimension, and the problems use periodic boundary conditions in space. 

* `startup.m`: Run this to put everything in your MATLAB path (NOTE: this will also modify your MGRIT plotting settings, at least until you re-open MATLAB, so comment out that part of the code if you don't want this to happen)

* `MGRIT/`: Implementation of the multigrid reduction-in-time algorithm.

* `scalar/`: Considers the scalar PDE u_t + f(u)_x = 0 where f(u) is the (potentially nonlinear) flux function. There are several examples, including:
    1. Linear conservation law:         f(u) = alpha(x,t)
    2. The Burgers equation:            f(u) = u^2/2
    3. The Buckley--Leverett equation:  f(u) = 4u^2 / (4u^2 + (1-u)^2)

* `utils/`: Various pieces of code for doing helpful things.