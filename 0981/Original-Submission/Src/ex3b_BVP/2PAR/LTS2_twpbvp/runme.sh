#!/bin/sh

echo ""
echo "runme.sh: shell script to run a driver program"
echo "           for OpenMP parallel Talbot Suite DE"
echo "           Example 3b - LT samples by twpbvp.F"
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
1)	# 1) BUILD THE EXECUTABLE (ACCURACY)
	echo "\n-------------------------------------------------\n"
	echo "1) compile and link with Makefile: make acc\n"
	make acc
	echo "\n-------------------------------------------------\n"
	# 2) RUN THE PROGRAM
	echo "\n-------------------------------------------------\n"
	echo "2) run the executable:\n"
	echo "                ./ex3_acc.exe tol thrds1\n"
	echo "          or"
	echo "                ./ex3_acc.exe tol thrds1 thrds2\n"
	echo "          where\n"
	echo "             tol: tolerance"
	echo "             thrds1,thrds2 are the number of parallel threads"
	echo "\n-------------------------------------------------\n"
	./ex3_acc.exe 1e-12 4 2
	;;
2)	# 1) BUILD THE EXECUTABLE (TIMES)
	echo "\n-------------------------------------------------\n"
	echo "1) compile and link with Makefile: make time\n"
	make time
	echo "\n-------------------------------------------------\n"
	# 2) RUN THE PROGRAM
	echo "\n-------------------------------------------------\n"
	echo "2) run the executable:\n"
	echo "                ./ex3_time.exe tol jFUN NTval NXval\n"
	echo "          where"
	echo "                tol : tolerance"
	echo "                jFUN: 1 for OMP_Talbot11 (coarse-grained parallelism)"
	echo "                      2 for OMP_Talbot12 ( fine-grained  parallelism)"
	echo "                      3 for OMP_Talbot13 (    nested     parallelism)"
	echo "                NTval: number of t-values"
	echo "                NXval: number of x-values"
	echo "\n-------------------------------------------------\n"
	./ex3_time.exe 1e-12 1 120 120
	;;
esac



# REMOVE EXECUTABLE AND OBJECT FILES
echo "\n-------------------------------------------------\n"
echo "3) clean executable and object files: make clean\n"
make clean
echo "\n-------------------------------------------------\n"
