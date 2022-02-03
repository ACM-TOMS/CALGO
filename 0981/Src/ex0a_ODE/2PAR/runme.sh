#!/bin/sh

clear
echo "\n-------------------------------------------------\n"
echo ""
echo "runme.sh: shell script to build and run a sample program"
echo "           for Talbot Suite DE implementations"
echo ""
echo "\n-------------------------------------------------\n"


# 1) BUILD THE EXECUTABLE
echo "\n-------------------------------------------------\n"
echo "1) compile and link with Makefile: make OMP\n"
make OMP
echo "\n-------------------------------------------------\n"


# 2) RUN THE PROGRAM
echo "\n-------------------------------------------------\n"
echo ""
echo "RUN THE PROGRAM\n"
echo "\nWhich problem?\n"
echo "        1)      20 t in [1000, 3000]"
echo "        2)     120 t in [1000, 3000]"
echo "        3)     120 t in [ 100,  500]"
echo -n "Enter a selection : "
read choice

# CHECK IF choice IS CORRECT, OTHERWISE EXIT
if [ $choice -lt 1 ] || [ $choice -gt 3 ]; then
	echo "\n\n***   wrong choice ==> exit \n"
	exit
fi # NOW choice IS CORRECT

echo "\n-------------------------------------------------\n"
echo ""
echo "\nWhich Talbot Suite function?\n"
echo "        0 : SEQ_Talbot1_DE"
echo "        1 : OMP_Talbot11_DE"
echo "        2 : OMP_Talbot12_DE"
echo "        3 : OMP_Talbot13_DE"
echo "        4 : SEQ_Talbot2_DE"
echo -n "Enter a selection : "
read jFUN

case $choice in
	1)	echo "Inverting at 20 t in [1000, 3000],  tol=1e-12"
		./ex0.exe 1 $jFUN
		;;
	2)	echo "Inverting at 120 t in [1000, 3000],  tol=1e-12"
		./ex0.exe 2 $jFUN
		;;
	3)	echo "Inverting at 120 t in [100, 500],  tol=1e-12"
		./ex0.exe 3 $jFUN
		;;
esac
echo ""
echo "\n-------------------------------------------------\n"


# 3) REMOVE EXECUTABLE AND OBJECT FILES
echo "\n-------------------------------------------------\n"
echo "3) clean executable and object files: make clean\n"
make clean
echo "\n-------------------------------------------------\n"
