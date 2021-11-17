#--------------------------------------------------------------------
#
#    file to compile and link drivers for the package PREQN:
#  
#    driver1: is the small test driver
#    driver2: is the single driver
#
#--------------------------------------------------------------------
#
# 
#

$FC = epcf90
$FFLAGS = -sloppy -C -d5 -g -temp=/tmp -u

all: res1 res2

driver1: driver1.o extras.o blas.o src.o
	$(FC) $(FFLAGS) driver1.o extras.o blas.o src.o -o driver1

driver2: driver2.o extras.o blas.o src.o
	$(FC) $(FFLAGS) driver2.o extras.o blas.o src.o -o driver2

res1: driver1 data1
	driver1 < data1 > res1
res2: driver2 data2 data2.1 data2.2 data2.3 data2.4 data2.5 data2.6
	driver2 < data2 > res2
