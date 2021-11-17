#!/bin/bash

# script file to obtain timing results for 4th order ridc using 4 cores

ORDER=4
EXE=brusselator_mkl
Nx=1000
export OMP_THREADS=${ORDER}

# convergence study
STEPS=100
echo "running order ${ORDER} ridc with nt = ${STEPS}" >> ${EXE}.log
time ./${EXE} ${ORDER} ${STEPS} ${Nx} > n1.dat

STEPS=200
echo "running order ${ORDER} ridc with nt = ${STEPS}" >> ${EXE}.log
time ./${EXE} ${ORDER} ${STEPS} ${Nx} > n2.dat

STEPS=400
echo "running order ${ORDER} ridc with nt = ${STEPS}" >> ${EXE}.log
time ./${EXE} ${ORDER} ${STEPS} ${Nx} > n3.dat

STEPS=800
echo "running order ${ORDER} ridc with nt = ${STEPS}" >> ${EXE}.log
time ./${EXE} ${ORDER} ${STEPS} ${Nx} > n4.dat

STEPS=1600
echo "running order ${ORDER} ridc with nt = ${STEPS}" >> ${EXE}.log
time ./${EXE} ${ORDER} ${STEPS} ${Nx} > n5.dat


