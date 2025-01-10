This project is about the parallel-in-time solution of hyperbolic PDEs. Everything is limited to one spatial dimension, and the problems use periodic boundary conditions in space (at least in most cases).

Presently, there are two papers associated with this repository:
1. `Parallel-in-time solution of scalar nonlinear conservation laws` by H. De Sterck, R. D. Falgout, O. A. Krzysik, J. B. Schroder
    See arXiv:2401.04936
    The code associated with this paper is in `scalar/`, as discussed below
2. `Parallel-in-time solution of hyperbolic PDE systems via characteristic-variable block preconditioning` by H. De Sterck, R. D. Falgout, O. A. Krzysik, J. B. Schroder
    See arXiv:2407.03873
    The code associated with this paper is in `systems/`, as discussed below 
  

Before running any of the examples in the directories described below, run `startup.m` to put everything in your MATLAB path and to configure your plotting settings. 

The two main directories in this repository are:
* `scalar/`: Considers the scalar PDE u_t + f(u)_x = 0 where f(u) is the (potentially nonlinear) flux function. There are several examples, including:
    1. Linear conservation law:         f(u) = alpha(x,t)
    2. The Burgers equation:            f(u) = u^2/2
    3. The Buckley--Leverett equation:  f(u) = 4u^2 / (4u^2 + (1-u)^2)
    
* `systems/`: Considers PDE systems including:
    1. The acoustic equations q_t + A(x)q_x = 0 -- a linear, non-conservative system in 2 variables
    2. The general nonlinear conservation law q_t +  f(q)_x = 0 -- examples include the shallow water equations and the Euler equations of gas dynamics

Additional information:
* `MGRIT/`: Implementation of the multigrid reduction-in-time algorithm.
* `utils/`: Various pieces of code for doing helpful things.
