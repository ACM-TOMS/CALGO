#
# Makefile for benchmark program to triangularize 2n x n random matrices
# using Givens transformations.
#
# This makefile builds the binaries for the single and double precision
# versions of the benchmark code.
# 
# To run a particular benchmark run
#           ./dtimer < data
# for the double precision benchmark (use stime for the single precision
# version)
#
# The data file should contain the following lines of data
#    line 1: the problem size N (this tests for the efficiency of the
#            transformation using vectors of size 1 .. N-1)
#    line 2: the number of randomly generated problems to be run, NRAND
#    line 3: the compiler flags used (printed in the results file)
#    line 4: the number of times the complete test should be run, NTIMES
#            (NRAND*NTIMES problems are run in total)
#    line 5: name of the output file for the results
#
# For example: 
#           ./dtimer < data_dp
#           ./stimer < data_sp
#
# will generate two results files, dbenchmark.res and sbenchmark.res
#
.SUFFIXES: .f90

# F90 should be set to the name of the Fortran 90 compiler to be used
F90         = f95
F90         = /opt/NAGWare/bin/f95

# F90OPTS should be set to compiler flags to be used
F90OPTS     = -C -g -u -xcheck=stkovf -ansi -XlistE
F90OPTS     = -g90 -gline -dcfuns -nan -save -C

# F90LINK and F90LINKOPTS should be set to the relevant linker and
# linker options. For most systems these are the same as for the compiler
F90LINK     = $(F90)
F90LINKOPTS = $(F90OPTS)


.f90.o:
	$(F90) $(F90OPTS) -c $*.f90
.f.o:
	$(F90) $(F90OPTS) -c $*.f

STIMEROBJS=Import_Kinds.o Givens_Rotations.o SBenchmark.o srotmg.o blas.o
DTIMEROBJS=Import_Kinds.o Givens_Rotations.o DBenchmark.o drotmg.o blas.o

all: dtimer stimer


stimer: $(STIMEROBJS)
	$(F90LINK) -o stimer $(F90LINKOPTS) $(STIMEROBJS)

dtimer: $(DTIMEROBJS)
	$(F90LINK) -o dtimer $(F90LINKOPTS) $(DTIMEROBJS)

Import_Kinds.o: ../Src/Import_Kinds.f90
	$(F90) $(F90OPTS) -c ../Src/Import_Kinds.f90

Givens_Rotations.o: Import_Kinds.o Givens_Rotations.f90
	$(F90) $(F90OPTS) -c Givens_Rotations.f90

DBenchmark.o: DBenchmark.f90
	$(F90) $(F90OPTS) -c DBenchmark.f90

SBenchmark.o: SBenchmark.f90
	$(F90) $(F90OPTS) -c SBenchmark.f90

blas.o: ../../Remark/Drivers/blas.f
	$(F90) $(F90OPTS) -c ../../Remark/Drivers/blas.f

drotmg.o: ../../Remark/Src/drotmg.f
	$(F90) $(F90OPTS) -c ../../Remark/Src/drotmg.f

srotmg.o: ../../Remark/Src/srotmg.f
	$(F90) $(F90OPTS) -c ../../Remark/Src/srotmg.f

clean:
	rm -rf ./stimer ./dtimer
	rm -rf *.o *.g90 *.d *.mod ./a.out core
