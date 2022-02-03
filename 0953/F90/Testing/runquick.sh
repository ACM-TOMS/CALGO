#!/bin/sh
#
# This file needs to be modified according to the computational environment.
#

MPIRUN=mpirun
NP=-np

rm -f quick.out

make

$MPIRUN $NP 1 ./test_random.exe < quick1.in >> quick.out
$MPIRUN $NP 1 ./test_bench.exe < quick2.in >> quick.out
$MPIRUN $NP 4 ./test_random.exe < quick1.in >> quick.out
$MPIRUN $NP 4 ./test_bench.exe < quick2.in >> quick.out
$MPIRUN $NP 9 ./test_random.exe < quick1.in >> quick.out
$MPIRUN $NP 9 ./test_bench.exe < quick2.in >> quick.out
$MPIRUN $NP 16 ./test_random.exe < quick1.in >> quick.out
$MPIRUN $NP 16 ./test_bench.exe < quick2.in >> quick.out
$MPIRUN $NP 36 ./test_random.exe < quick1.in >> quick.out
$MPIRUN $NP 36 ./test_bench.exe < quick2.in >> quick.out
$MPIRUN $NP 2 ./test_random.exe < quick1.in >> quick.out
$MPIRUN $NP 2 ./test_bench.exe < quick2.in >> quick.out
$MPIRUN $NP 3 ./test_random.exe < quick1.in >> quick.out
$MPIRUN $NP 3 ./test_bench.exe < quick2.in >> quick.out
$MPIRUN $NP 6 ./test_random.exe < quick1.in >> quick.out
$MPIRUN $NP 6 ./test_bench.exe < quick2.in >> quick.out
$MPIRUN $NP 12 ./test_random.exe < quick1.in >> quick.out
$MPIRUN $NP 12 ./test_bench.exe < quick2.in >> quick.out
$MPIRUN $NP 24 ./test_random.exe < quick1.in >> quick.out
$MPIRUN $NP 24 ./test_bench.exe < quick2.in >> quick.out

grep "out of" quick.out
