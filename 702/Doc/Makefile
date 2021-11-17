#----------------------------------------------------------
#       Makefile for TNPACK tests
#---------------------------------------------------------

F77    =  epcf90
FFLAGS =  epcf90 -sloppy -C -d5 -g -temp=/tmp -u

.SUFFIXES:      .f .o
.f.o:;$(F77) $(FFLAGS) -c .f


test1:	tn1.o tnpack.o blas1.o
	$(F77) $(LFLAGS) -o   tn1.o tnpack.o blas1.o

test2: tn2.o tnpack.o blas1.o
	$(F77) $(LFLAGS) -o   tn2.o tnpack.o blas1.o

test3:  tn3.o tnpack.o alg566.o blas1.o
	$(F77) $(LFLAGS) -o   tn3.o tnpack.o alg566.o blas1.o
