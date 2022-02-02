AR= ar
F90C= f95
FFLAGS= -c -u
INC1= lapack
INC2= condnum

all : mylapack mycondnum babdcr clean

# compile BABDCR files and create the executable file babdcr
babdcr : 
	$(F90C) $(FFLAGS) precision.f90 babd_matmult.f90 babdcr.f90 \
	babdcr_solvt.f90 babdcr_factsolv.f90 babdcr_condnum.f90
	$(F90C) -o babdcr main.f90 *.o -L$(INC1) -lmylapack -L$(INC2) -lmycondnum

# compile  all blas/lapack routines used in the package
# and save in mylapack.a
mylapack : 
	make libmylapack.a

libmylapack.a :
	cd $(INC1); $(F90C) $(FFLAGS) *.f; $(AR) -r $@ *.o; rm *.o

# compile all the routines used to compute the condition number 
# of the BABD matrix and save in mycondnum.a
mycondnum : 
	make libmycondnum.a

libmycondnum.a :
	cd $(INC2); $(F90C) -c *.f; $(AR) -r $@ *.o; rm *.o

#delete all files except babdcr
clean : 
	rm -f *.mod *.o $(INC1)/*.a $(INC2)/*.a

# delete all files
cleanall : 
	rm -f babdcr *.mod *.o $(INC1)/*.a $(INC2)/*.a
