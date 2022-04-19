#!/bin/sh

clear
echo "\n-------------------------------------------------\n"
echo ""
echo "runme.sh: shell script to build and run a sample program"
echo "          for Talbot Suite DE implementations"
echo ""
echo "\n-------------------------------------------------\n"


# 1) BUILD THE EXECUTABLE
echo "\n-------------------------------------------------\n"
echo "     Compile and link with Makefile: make SEQ"
echo "\n-------------------------------------------------\n"
make SEQ
echo ""


# 2) CHOOSE THE TEST
echo "\n-------------------------------------------------\n"
echo "     Choose the test"
echo "\n-------------------------------------------------\n"
echo ""
echo "\nWhich problem?\n"
echo "        1)      20 t in [1000, 3000]"
echo "        2)     120 t in [1000, 3000]"
echo "        3)     120 t in [ 100,  500]"
echo "        4)      20 t in [ 100,  500]"
echo -n "Enter a selection : "
read choice

# CHECK IF choice IS CORRECT, OTHERWISE EXIT
if [ $choice -lt 1 ] || [ $choice -gt 4 ]; then
	echo "\n\n***   wrong choice ==> exit \n"
	exit
fi # NOW choice IS CORRECT


# 3) RUN THE PROGRAM
echo "\n-------------------------------------------------\n"
echo "     Run the program: ./ex0.exe $choice"
echo "\n-------------------------------------------------\n"
echo ""

case $choice in
	1)	echo "Inverting at 20 t in [1000, 3000],  tol=1e-12"
		./ex0.exe 1
		;;
	2)	echo "Inverting at 120 t in [1000, 3000],  tol=1e-12"
		./ex0.exe 2
		;;
	3)	echo "Inverting at 120 t in [100, 500],  tol=1e-12"
		./ex0.exe 3
		;;
	4)	echo "Inverting at 20 t in [100, 500],  tol=1e-12"
		./ex0.exe 4
		;;
esac
echo ""
echo "\n-------------------------------------------------\n"


# 3) REMOVE EXECUTABLE AND OBJECT FILES
echo "\n-------------------------------------------------\n"
echo "3) clean executable and object files: make clean\n"
make clean
echo "\n-------------------------------------------------\n"
