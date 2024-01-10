In this directory we consider the solution of the scalar conservation law u_t + f(u)_x = 0, with possibly nonlinear flux function f.

* `cons_law_scalar.m`: Abstract class implementing PDE and its discretization. Specific PDEs implemented using this class are: 
    1. The linear conservation law:   f(u) = alpha(x, t)*u for some prescribed function alpha. See `cons_lin_scalar.m`
    2. The Burgers equation:          f(u) = u^2/2. See `burgers.m`
    3. The Buckley--Levertt equation: f(u) = 4u^2 ./ (4u^2 + (1-u)^2). See `buckley_leverett.m`

* `cons_law_time_stepping.m`: Solves the PDE on some time interval using time-stepping

* `cons_law_st.m`: Solve the PDE on some time interval using a preconditioned residual correction scheme applied to the whole space-time system. The linearized problems at each iteration may be solved directly via sequential time-stepping or approximately with MGRIT. When MGRIT is used as the linear solver, the associated MGRIT stepping methodology is implemented in `step_cons_law_linearized_MGRIT.m`

* `cons_law_accuracy_test.m`: Measures discretization error for a linear and Burgers equations.

* `cons_lin_st_MGRIT.m`: Uses MGRIT to solve the linear conservation law e_t + (alpha(x, t)*e)_x = 0, where alpha(x, t) is some prescribed function, and the PDE is discretized with a standard linear MOL discretization. The coarse-grid operator is a modified, conservative semi-Lagrangian method. All the stepping functionality for MGRIT is implemented in `step_cons_lin_MGRIT.m`




