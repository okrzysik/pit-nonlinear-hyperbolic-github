MGRIT implementation

* `mgrit.m` Implements the multigrid reduction-in-time algorithm to solve the linear system A*u = g, where A is a lower bi-diagonal matrix with some (potentially) time-dependent, linear time-stepping operator.

Note that the implementation is serial. 

* `mgrit_heat_example.m` Shows an example using the code to solve a one-dimensional heat equation.

