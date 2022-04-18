README - EXAMPLE

This subfolder contains sample programs for Talbot Suite's functions.
Each main program calls a single user-level function in Talbot Suite.
Talbot Suite's functions are located in SRC subfolder of Talbot_Suite.

In particular

    - OMP implementation:
      * OMP_main1 for OMP_Talbot1 function
      * OMP_main2 for OMP_Talbot2 function

    - MPI implementation:
      * MPI_main1 for MPI_Talbot1 function
      * MPI_main2 for MPI_Talbot2 function

    - HYB implementation:
      * HYB_main for HYB_Talbot3 function


The subfolder also contains:

*   an example of Makefile (to build the executable by means of the gcc compiler);

*   an example of shell script (run-me.sh) which, interactively, manages the
    compiling/linking phase for a particular Talbot Suite function and launches
    the executable setting the required parameters.
    The user is asked to choose which implementation has to be used (OMP, MPI, HYB)
    and, if required, the number of parallel processes and the Talbot Suite function
    to be used.
    Make the run-me.sh shell script executable.
