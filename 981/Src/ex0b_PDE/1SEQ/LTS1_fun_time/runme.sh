#!/bin/sh

clear
echo ""
echo "runme.sh: shell script to build and run a driver program"
echo "          for Talbot Suite DE"
echo "          Example 0b [Duffy] - LT samples by a function"
echo "                  Test on execution times"
echo "     Comparison between Talbot Suite and Talbot Suite DE"
echo ""


# BUILD THE EXECUTABLE
echo "\n-------------------------------------------------\n"
echo "     Compile and link with Makefile: make time"
make time
echo "\n-------------------------------------------------\n"


# RUN THE PROGRAM
echo "\n-------------------------------------------------\n"
echo "     Run the executable: ./ex1_time.exe 1e-6"
echo "\n-------------------------------------------------\n"
./ex0_time.exe 1e-6


# REMOVE EXECUTABLE AND OBJECT FILES
echo "\n-------------------------------------------------\n"
echo "3) clean executable and object files: make clean\n"
make clean
echo "\n-------------------------------------------------\n"
