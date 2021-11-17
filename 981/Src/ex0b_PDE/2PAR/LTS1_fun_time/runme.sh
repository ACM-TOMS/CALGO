#!/bin/sh

echo ""
echo "runme.sh: shell script to run a driver program"
echo "           for OpenMP parallel Talbot Suite DE"
echo "           Example 0b [Duffy] - LT samples by a function"
echo "                  Test on execution times"
echo "     Comparison between Talbot Suite and Talbot Suite DE"
echo ""


# BUILD THE EXECUTABLE
echo "\n-------------------------------------------------\n"
echo "\tCompile and link with Makefile: make time\n"
make time
echo "\n-------------------------------------------------\n"


# RUN THE PROGRAM
echo "\n-------------------------------------------------\n"
echo "\tRun the executable as: ./ex0.exe 1e-12 1 120 120\n"
echo "\twhere  tol  = 1e-12"
echo "\t      jFUN  = 1 (OMP_Talbot11 function)"
echo "\t      NTval = 120 (number of t-values in [0.5, 20])"
echo "\t      NXval = 120 (number of x-values in [0, 1])"
echo ""
./ex0.exe 1e-12 1 120 120
echo "\n-------------------------------------------------\n"


# REMOVE EXECUTABLE AND OBJECT FILES
echo "\n-------------------------------------------------\n"
echo "\tClean executable and object files: make clean\n"
make clean
echo "\n-------------------------------------------------\n"
