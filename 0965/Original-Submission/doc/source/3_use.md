# Using the RIDC Library                                                                

To utilize the RIDC library, a main program should specify the ODE
class, which specifies the system of ODEs to be solved, the ODE
parameters (including the number of equations, number of time steps,
size of the time step, and the initial and final time), as well as the
step function: an Euler integerator that advances the solution from
time t(n) to t(n+1).  This step routine may be complicated requiring
large scale linear algebra provided by external external libraries or
possibly a nonlinear solve.  The `examples/brusselator_gsl` directory
contains such an example.  This example uses a backward Euler step for
a nonlinear system of ODEs.  The step function uses a Newton iteration
(see the functions `newt` and `jac`) and the GNU scientific library
(GSL) to solve for the Newton step.  The functions `newt` and `jac`
required by step are defined and declared in `brusselator.h`.
Finally, the solution is integrated using a call to the `ridc_fe`
or `ridc_be` functions.


To link against the RIDC library, include the
following arguments to the compilation command:

`-L/home/user/opt/libridc/lib -I/home/user/opt/libridc/include -lridc`

