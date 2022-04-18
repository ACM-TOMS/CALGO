README - Talbot Suite

This folder contains four subfolders:


    1) ./Src
        It contains all the C functions and header files of the Talbot Suite
        implementations. In particular:

        - COM_Talbot_pack.c [COM_Talbot_pack.h] for the shared utility functions;

        - OMP_Talbot_pack.c [OMP_Talbot_pack.h] for the OpenMP implementation;

        - MPI_Talbot_pack.c [MPI_Talbot_pack.h] for the MPI implementation;

        - HYB_Talbot_pack.c [HYB_Talbot_pack.h] for the MPI/OMP hybrid implementation;

        - main_include.h: header file used by all the implementations.

        To use the OMP implementation of Talbot Suite only the following files are
        required: COM_Talbot_pack (.c/.h), OMP_Talbot_pack (.c/.h) and main_include.h.
        The same holds for the others.


    2) ./Examples
        It contains sample main programs. Each of them is related to just a single function
        in Talbot Suite: for example, MPI_main2.c calls MPI_Talbot2, which implements
        the fine grain parallel algorithm in the MPI-based implementation. 
        A Makefile (for the gcc compiler) and a shell script are provided: they are intended
        as examples to build executables and to run them.


    3) ./Drivers
        It contains driver programs. These codes are addressed to the referee for testing
        the code, because, unlike the sample programs, they also compute absolute
        and relative errors in numerical approximations and elapsed times.
        The test function is the same as in the paper.
        A Makefile (for the gcc compiler) and a shell script are provided: they are intended
        as examples to build executables and to run them.


    4) ./{Drivers,Examples}/Results
        It contains examples of outputs.


Make the run-me.sh shell script executable before launching it.


Main differences between sample and driver programs consist of:

    1)  Each sample main calls a single function in Talbot Suite, while each
        driver program is designed for a particular implementation, switching
        if necessary between the two functions (such as OMP_Talbot1/2 or MPI_Talbot1/2).

    2)  Driver programs compute elapsed times, while sample programs do not.

    3)  Driver programs compare the exact inverse LT function values with
        their numerical approximations computing absolute and relative errors,
        while sample programs do not.

    4)  The Laplace Transform function is different between sample and driver
        programs. In driver programs it is the same as in the paper.
