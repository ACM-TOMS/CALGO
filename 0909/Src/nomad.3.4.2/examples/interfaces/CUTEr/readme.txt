
This example illustrates how to use NOMAD with a CUTEr test problem,
given in a SIF files.

You need to have CUTEr installed, including a SIF decoder.
Everything is available at http://hsl.rl.ac.uk/cuter-www/.

Once the CUTEr installation is complete, define these environment variables:

setenv CUTER      /home/user_name/CUTEr
setenv SIFDEC     $CUTER/sifdec
setenv MYSIFDEC   $SIFDEC

Then use the sifdecode program on the problem file:
$SIFDEC/double/bin/sifdecode PROBLEM.SIF

Compile with the C wrapper 'bb.c' by running the script 'compile'
(it uses gfortran and gcc, but you can use other C and Fortran compilers).

The black-box executable 'bb.exe' should be created.

You can test it with the command 'bb.exe x0.txt', and run nomad with the command
'nomad parameters.txt'.