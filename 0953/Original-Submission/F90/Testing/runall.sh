#!/bin/sh
#
# This file needs to be modified according to the computational environment.
#
# WARNING: This test will take a VERY LONG time!
#

MPIRUN=mpirun
NP=-np

rm -f rand*.out bench*out rect*out big*out

make

#
# Make sure the following four benchmark matrices have been downloaded.
#
wget -c ftp://math.nist.gov/pub/MatrixMarket2/NEP/olmstead/olm5000.mtx.gz
wget -c ftp://math.nist.gov/pub/MatrixMarket2/NEP/dwave/dw8192.mtx.gz
wget -c ftp://math.nist.gov/pub/MatrixMarket2/NEP/crystal/cry10000.mtx.gz
wget -c ftp://math.nist.gov/pub/MatrixMarket2/NEP/airfoil/af23560.mtx.gz
gunzip olm5000.mtx.gz
gunzip dw8192.mtx.gz
gunzip cry10000.mtx.gz
gunzip af23560.mtx.gz

# 1x1 processor grid
$MPIRUN $NP 1 ./test_random.exe < rand1.in > rand_1x1.out
$MPIRUN $NP 1 ./test_bench.exe < bench1.in > bench_1x1.out

# 2x2 processor grid
$MPIRUN $NP 4 ./test_random.exe < rand1.in > rand_2x2.out
$MPIRUN $NP 4 ./test_bench.exe < bench2.in > bench_2x2.out

# 3x3 processor grid
$MPIRUN $NP 9 ./test_bench.exe < bench3.in > bench_3x3.out

# 4x4 processor grid
$MPIRUN $NP 16 ./test_random.exe < rand2.in > rand_4x4.out
$MPIRUN $NP 16 ./test_bench.exe < bench4.in > bench_4x4.out

# 6x6 processor grid
$MPIRUN $NP 36 ./test_random.exe < rand3.in > rand_6x6.out
$MPIRUN $NP 36 ./test_bench.exe < bench4.in > bench_6x6.out

# 8x8 processor grid
$MPIRUN $NP 64 ./test_random.exe < rand3.in > rand_8x8.out
$MPIRUN $NP 64 ./test_bench.exe < bench4.in > bench_8x8.out

# 10x10 processor grid
$MPIRUN $NP 100 ./test_random.exe < rand3.in > rand_10x10.out
$MPIRUN $NP 100 ./test_bench.exe < bench4.in > bench_10x10.out

# 12x12 processor grid
$MPIRUN $NP 144 ./test_bench.exe < bench12.in > bench_12x12.out

# rectangle processor grids
$MPIRUN $NP 2 ./test_random.exe < quick1.in >> rect.out
$MPIRUN $NP 2 ./test_bench.exe < quick2.in >> rect.out
$MPIRUN $NP 6 ./test_random.exe < quick1.in >> rect.out
$MPIRUN $NP 6 ./test_bench.exe < quick2.in >> rect.out
$MPIRUN $NP 12 ./test_random.exe < quick1.in >> rect.out
$MPIRUN $NP 12 ./test_bench.exe < quick2.in >> rect.out
$MPIRUN $NP 24 ./test_random.exe < quick1.in >> rect.out
$MPIRUN $NP 24 ./test_bench.exe < quick2.in >> rect.out

# Very big problems
#$MPIRUN $NP 256 ./test_random.exe < rand4.in >> big.out
#$MPIRUN $NP 576 ./test_random.exe < rand4.in >> big.out
#$MPIRUN $NP 1024 ./test_random.exe < rand4.in >> big.out
#$MPIRUN $NP 1600 ./test_random.exe < rand4.in >> big.out

grep "out of" *.out | tee summary.txt
