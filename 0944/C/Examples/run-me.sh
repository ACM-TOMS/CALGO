#!/bin/sh

# INITIALIZE
nprocs=1		# number of MPI processes
threads=1		# number of OMP threads

echo ""
echo "run_me.sh: shell script to run sample main programs"
echo "           for Talbot Suite implementations"
echo ""
echo "\nWhich implementation of Talbot Suite?\n"
echo "        1)  OMP, run as SEQUENTIAL"
echo "        2)  OMP"
echo "        3)  MPI"
echo "        4)  HYB"
echo -n "Enter a selection : "
read choice

# CHECK IF choice IS CORRECT, OTHERWISE EXIT
if [ $choice -lt 1 ] || [ $choice -gt 4 ]; then
	echo "\n\n***   wrong choice ==> exit \n"
	exit
fi # NOW choice IS CORRECT

# SET THE NUMBER OF PARALLEL PROCESSES
if [ $choice -gt 1 ]; then
	echo "\n\nHow many parallel processes?\n"
	case $choice in
	2)	echo "Talbot Suite's OMP implementation"
		echo -n "enter the number of OMP threads [1,2,...]: "
		read threads
		Tver='OMP'
		;;
	3)	echo "Talbot Suite's MPI implementation"
		echo -n "enter the number of MPI procs [1,2,...]: "
		read nprocs
		Tver='MPI'
		;;
	4)	echo "Talbot Suite's HYB implementation"
		echo -n "number of MPI procs [1,2,...]: "
		read nprocs
		echo -n "number of OMP threads [1,2,...]: "
		read threads
		Tver='HYB'
		;;
	esac
fi

# SELECT THE Talbot_Suite FUNCTION
if [ $choice -eq 4 ]
then
	# HYB
	Tfun=3
else
	if [ $choice -eq 1 ]
	then
		# OMP AS SEQUENTIAL
		echo "\n\nWhich Talbot's method?\n"
		echo "        1) Modified Talbot's method"
		echo "        2) Classical Talbot's method"
		Tver='OMP'
	else
		# PARALLEL (MPI OR OMP)
		echo "\n\nChoose the level of parallelism\n"
		echo "        1) Coarse grain parallelism"
		echo "        2) Fine grain parallelism"
	fi
	echo -n "Enter a selection : "
	read Tfun	# switch between _Talbot1 and _Talbot2
fi
echo "Selected function:  ${Tver}_Talbot$Tfun()"

if [ $choice -eq 1 ]; then
	Tver='SEQ'
fi


# 1) BUILD THE EXECUTABLE
echo "\n-------------------------------------------------\n"
echo "1) compile and link with Makefile: make ${Tver}$Tfun\n"
make ${Tver}$Tfun
echo "\n-------------------------------------------------\n"


# 2) RUN THE PROGRAM
echo "\n-------------------------------------------------\n"
case $choice in

#   1) OMP AS SEQUENTIAL
1)	echo "2) run the executable: ./${Tver}_talbot.exe"
	./${Tver}_talbot.exe
	;;

#   2) OMP
2)	echo "2) run the executable: ./OMP_talbot.exe $threads"
	./OMP_talbot.exe $threads
	;;

#   3) MPI
#   WITHOUT -machinefile OPTION TO mpiexec
3)	echo "2) run the executable: mpiexec -np $nprocs ./MPI_talbot.exe"
	mpiexec  -np $nprocs ./MPI_talbot.exe

#   OR, ALTERNATIVELY, WITH -machinefile OPTION TO mpiexec
#3)	echo "2) run the executable: mpiexec -machinefile hostlist -np $nprocs ./MPI_talbot.exe"
#	mpiexec -machinefile hostlist -np $nprocs ./MPI_talbot.exe
	;;

#   4) HYB
#   WITHOUT -machinefile OPTION TO mpiexec
4)	echo "2) run the executable: mpiexec -np $nprocs ./HYB_talbot.exe $threads"
	mpiexec  -np $nprocs ./HYB_talbot.exe $threads

#   OR, ALTERNATIVELY, WITH -machinefile OPTION TO mpiexec
#4)	echo "2) run the executable: mpiexec -machinefile hostlist -np $nprocs ./HYB_talbot.exe $threads"
#	mpiexec -machinefile hostlist -np $nprocs ./HYB_talbot.exe $threads
	;;

esac
echo "\n-------------------------------------------------\n"


# REMOVE EXECUTABLE AND OBJECT FILES
echo "\n-------------------------------------------------\n"
echo "3) clean executable and object files: make clean\n"
make clean
echo "\n-------------------------------------------------\n"
