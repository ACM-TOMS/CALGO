#!/bin/sh

echo ""
echo "runme.sh: shell script to run a driver program"
echo "           for Talbot Suite DE"
echo "           Example 1b - LT samples by ode.c"
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
	echo "1) compile and link with Makefile: make acc\n"
	make acc
	echo "\n-------------------------------------------------\n"
	# 2) RUN THE PROGRAM
	echo "\n-------------------------------------------------\n"
	echo "2) run the executable: ./ex1_acc.exe tol\n"
	echo "                       where tol is a tolerance"
	echo "\n-------------------------------------------------\n"
	./ex1_acc.exe 1e-8
	;;
2)	# 1) BUILD THE EXECUTABLE
	echo "\n-------------------------------------------------\n"
	echo "1) compile and link with Makefile: make time\n"
	make time
	echo "\n-------------------------------------------------\n"
	# 2) RUN THE PROGRAM
	echo "\n-------------------------------------------------\n"
	echo "2) run the executable: ./ex1_time.exe tol\n"
	echo "                       where tol is a tolerance"
	echo "\n-------------------------------------------------\n"
	./ex1_time.exe 1e-8
	;;
esac


# REMOVE EXECUTABLE AND OBJECT FILES
echo "\n-------------------------------------------------\n"
echo "3) clean executable and object files: make clean\n"
make clean
echo "\n-------------------------------------------------\n"
