#!/bin/sh
#
# Auto-tuning script.
#
# This file needs to be modified according to the computational environment.
#

MPIRUN=mpirun
NP=-np

rm -f tune2.out

make tune2

$MPIRUN $NP 4 ./test10.exe < tune2_1.in >> tune2.out
$MPIRUN $NP 4 ./test20.exe < tune2_1.in >> tune2.out
$MPIRUN $NP 4 ./test30.exe < tune2_1.in >> tune2.out
$MPIRUN $NP 4 ./test40.exe < tune2_1.in >> tune2.out
$MPIRUN $NP 4 ./test50.exe < tune2_1.in >> tune2.out
$MPIRUN $NP 4 ./test60.exe < tune2_1.in >> tune2.out
$MPIRUN $NP 4 ./test70.exe < tune2_1.in >> tune2.out
$MPIRUN $NP 4 ./test80.exe < tune2_1.in >> tune2.out
$MPIRUN $NP 4 ./test90.exe < tune2_1.in >> tune2.out
$MPIRUN $NP 16 ./test10.exe < tune2_2.in >> tune2.out
$MPIRUN $NP 16 ./test20.exe < tune2_2.in >> tune2.out
$MPIRUN $NP 16 ./test30.exe < tune2_2.in >> tune2.out
$MPIRUN $NP 16 ./test40.exe < tune2_2.in >> tune2.out
$MPIRUN $NP 16 ./test50.exe < tune2_2.in >> tune2.out
$MPIRUN $NP 16 ./test60.exe < tune2_2.in >> tune2.out
$MPIRUN $NP 16 ./test70.exe < tune2_2.in >> tune2.out
$MPIRUN $NP 16 ./test80.exe < tune2_2.in >> tune2.out
$MPIRUN $NP 16 ./test90.exe < tune2_2.in >> tune2.out
$MPIRUN $NP 64 ./test10.exe < tune2_3.in >> tune2.out
$MPIRUN $NP 64 ./test20.exe < tune2_3.in >> tune2.out
$MPIRUN $NP 64 ./test30.exe < tune2_3.in >> tune2.out
$MPIRUN $NP 64 ./test40.exe < tune2_3.in >> tune2.out
$MPIRUN $NP 64 ./test50.exe < tune2_3.in >> tune2.out
$MPIRUN $NP 64 ./test60.exe < tune2_3.in >> tune2.out
$MPIRUN $NP 64 ./test70.exe < tune2_3.in >> tune2.out
$MPIRUN $NP 64 ./test80.exe < tune2_3.in >> tune2.out
$MPIRUN $NP 64 ./test90.exe < tune2_3.in >> tune2.out

grep "%SOLVER,NIBBLE,N,NB,NPROW,NPCOL,TIME" tune2.out > summary2.txt
./tune2.exe < summary2.txt | tee suggestion2.txt
