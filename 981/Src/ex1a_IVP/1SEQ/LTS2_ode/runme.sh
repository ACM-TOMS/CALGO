#!/bin/sh

echo ""
echo "runme.sh: shell script to run a driver program"
echo "           for Talbot Suite DE"
echo "           Example 1a - LT samples by ode.c"
echo ""
echo "\nWhich test?\n"
echo "        1)  On accuracy"
echo "        2)  On execution times"
echo -n "Enter a selection : "
read choice

# CHECK IF choice IS CORRECT, OTHERWISE EXIT
if [ $choice -lt 1 ] || [ $choice -gt 2 ]; then
	echo "\n\n***   wrong choice ==> exit \n"
	exit
fi # NOW choice IS CORRECT
echo "choice = " $choice

case $choice in
1)	# 1) BUILD THE EXECUTABLE
	echo "\n-------------------------------------------------\n"
	echo "     Compile and link with Makefile: make acc"
	echo "\n-------------------------------------------------\n"
	make acc
	# 2) RUN THE PROGRAM
	echo "\n-------------------------------------------------\n"
	echo "     Run the executable: ./ex1_acc.exe 1e-8"
	echo "\n-------------------------------------------------\n"
	./ex1_acc.exe 1e-8
	;;
2)	# 1) BUILD THE EXECUTABLE
	echo "\n-------------------------------------------------\n"
	echo "    Compile and link with Makefile: make time"
	echo "\n-------------------------------------------------\n"
	make time
	# 2) RUN THE PROGRAM
	echo "\n-------------------------------------------------\n"
	echo "     Run the executable: ./ex1_time.exe 1e-8"
	echo "\n-------------------------------------------------\n"
	./ex1_time.exe 1e-8
	;;
esac


# REMOVE EXECUTABLE AND OBJECT FILES
echo "\n-------------------------------------------------\n"
echo "     Clean executable and object files: make clean"
echo "\n-------------------------------------------------\n"
make clean

