For the GNU gfortran compiler

In order to use OpenMP directives in FORTRAN code, we need to abilitate
pre-processing and conditional compilation in gfortran.
To do this, we changed the ".f" extension of FORTRAN code into ".F" and
we use #ifdef-#else-#endif as usual.
