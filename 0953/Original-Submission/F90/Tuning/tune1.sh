#!/bin/sh
#
# Auto-tuning script.
#
# This file needs to be modified according to the computational environment.
#

MPIRUN=mpirun
NP=-np

rm -f tune1.out

make tune1

$MPIRUN $NP 1 ./test.exe < tune1.in >> tune1.out
$MPIRUN $NP 4 ./test.exe < tune1.in >> tune1.out
$MPIRUN $NP 16 ./test.exe < tune1.in >> tune1.out
$MPIRUN $NP 64 ./test.exe < tune1.in >> tune1.out

grep "%SOLVER,N,NB,NPROW,NPCOL,TIME" tune1.out > summary1.txt
./tune1.exe < summary1.txt | tee suggestion1.txt
