# ###################################################################
#
#  This is the makefile for SRRIT - A FORTRAN subroutine to calculate 
#  the Dominant Invariant Subspace of a Nonsymmetric Matrix developed 
#  by Z. Bai at University of Kentucky, and G. W. Stewart at University
#  of Maryland.  
#
#  If you have installed BLAS and LAPACK on your local machine, you may 
#  replace lapacksub.o and blas.o by your BLAS and LAPACK libraries.
#
#  For any comments, please send to 
#          bai@ms.uky.edu or stewart@cs.umd.edu
#
# ###################################################################
#
.SUFFIXES: .f .o
FFLAGS = -sloppy -C -d1 -g -temp=/tmp -u
F77 = epcf90

.f.o:
	$(F77) -c $(FFLAGS) $*.f
BLAS = 
LAPACK =
#
SRRIT = srrit.o sblas.o slapack.o
DSRRIT = dsrrit.o dblas.o dlapack.o
MAKE = make single double
#

#all:	
#	( $(MAKE); \
#	  subspiter < data; \
#	  dsubspiter < data )
#
single: driver.o $(SRRIT) ;\
	$(F77) -g -o subspiter driver.o $(SRRIT) $(LAPACK) $(BLAS)

double: driver.o $(DSRRIT) ;\
	$(F77) -g -o dsubspiter driver.o $(DSRRIT) $(LAPACK) $(BLAS)

clean:
	rm -f *.o
