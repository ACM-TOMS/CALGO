# Running the Examples                                                                  {#running}

The directory `examples/` includes five examples of utilizing the RIDC
library, and one example, `examples/brusselator_radau_mkl` that
implements a three stage, fifth-order Radau method to provide a basis
of comparison with the RIDC integrators.  Depending on the options
selected in the `./configure` step, some or all of these examples are
built and run during during the `make check` process.  Alternatively,
a user can compile and run an example seperately after the
`./configure` step.  For example, the subdirectory
`examples/explicit/` contains the code to solve a system of ODES using
RIDC with an explicit Euler step function.  To compile this specific
example, move into the `examples/explicit` subdirectory and type `make
explicit`.  The executable `explicit` takes as input the order
required and the number of time steps. For example `./explicit.exe 4
100` solves the system of ODEs using fourth order RIDC with 100 time
steps.  A shell script `run.sh` is provided to run the RIDC integrator
with different numbers of time steps for a convergence study.  A
simple matlab or octave script `convergence.m` is included in that
subdirectory to test the order of convergence.  `octave
convergence.m` gives the slope and intercept for the linear fit of log
of the error versus log of the time step.  In this example we obtain a
slope of -4.0630 indicating the we indeed have an order 4 method.

