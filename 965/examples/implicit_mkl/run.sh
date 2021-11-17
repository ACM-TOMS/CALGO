#!/bin/bash

ORDER=4
EXE=implicit_mkl
export OMP_THREADS=${ORDER}

# convergence study
STEPS=10
echo "running order ${ORDER} ridc with nt = ${STEPS}" >> ${EXE}.log
./${EXE} ${ORDER} ${STEPS} > n1.dat
STEPS=20
echo "running order ${ORDER} ridc with nt = ${STEPS}" >> ${EXE}.log
./${EXE} ${ORDER} ${STEPS} > n2.dat
STEPS=40
echo "running order ${ORDER} ridc with nt = ${STEPS}" >> ${EXE}.log
./${EXE} ${ORDER} ${STEPS} > n3.dat
STEPS=80
echo "running order ${ORDER} ridc with nt = ${STEPS}" >> ${EXE}.log
./${EXE} ${ORDER} ${STEPS} > n4.dat
STEPS=160
echo "running order ${ORDER} ridc with nt = ${STEPS}" >> ${EXE}.log
./${EXE} ${ORDER} ${STEPS} > n5.dat


rm -f diff.log

# compare with reference solution
diff n1.dat ref/n1.ref >> diff.log
diff n2.dat ref/n2.ref >> diff.log
diff n3.dat ref/n3.ref >> diff.log
diff n4.dat ref/n4.ref >> diff.log
diff n5.dat ref/n5.ref >> diff.log


if [ -s diff.log ]
then
	echo ":test-result: FAIL" >> ${EXE}.trs
else
	echo ":test-result: PASS" >> ${EXE}.trs
fi
