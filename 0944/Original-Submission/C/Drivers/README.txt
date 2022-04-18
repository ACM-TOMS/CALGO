README - DRIVER

This subfolder contains driver programs; each of them is related to a particular
implementation of Talbot Suite:

    - OMP_main for OMP_Talbot1 and OMP_Talbot2 functions

    - MPI_main for MPI_Talbot1 and MPI_Talbot2 functions

    - HYB_main for HYB_Talbot3 function

Talbot Suite's functions are located in SRC subfolder of Talbot_Suite.


There are also:

*   an example of Makefile (to build the executable by means of the gcc compiler);

*   an example of shell script (run-me.sh) which, interactively, manages the
    compiling/linking phase for a particular Talbot Suite function and launches
    the executable setting the required parameters.
    The user is asked to choose which implementation has to be used (OMP, MPI, HYB)
    and, if required, the number of parallel processes and the Talbot Suite function
    to be used.
    Make the run-me.sh shell script executable.
