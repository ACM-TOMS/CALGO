#!/bin/sh
#
# This file needs to be modified according to the computational environment.
#

MPIRUN=mpirun
NP=-np

rm -f quick.out

make

$MPIRUN $NP 1 ./test_random.exe < quick1.in > quick.out
$MPIRUN $NP 1 ./test_bench.exe < quick2.in >> quick.out

grep "out of" quick.out
