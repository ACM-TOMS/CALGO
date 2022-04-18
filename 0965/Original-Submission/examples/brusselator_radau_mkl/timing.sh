#!/bin/bash

EXE=brusselator_radau

# convergence study
STEPS=100
Nx=1000

echo "running radau with with nt = ${STEPS}, neq=${Nx}" >> ${EXE}.log
time ./${EXE} ${STEPS} ${Nx} > n1.dat

STEPS=200
echo "running radau with with nt = ${STEPS}, neq=${Nx}" >> ${EXE}.log
time ./${EXE} ${STEPS} ${Nx} > n2.dat

STEPS=400
echo "running radau with with nt = ${STEPS}, neq=${Nx}" >> ${EXE}.log
time ./${EXE} ${STEPS} ${Nx} > n3.dat

STEPS=800
echo "running radau with with nt = ${STEPS}, neq=${Nx}" >> ${EXE}.log
time ./${EXE} ${STEPS} ${Nx} > n4.dat

STEPS=1600
echo "running radau with with nt = ${STEPS}, neq=${Nx}" >> ${EXE}.log
time ./${EXE} ${STEPS} ${Nx} > n5.dat

