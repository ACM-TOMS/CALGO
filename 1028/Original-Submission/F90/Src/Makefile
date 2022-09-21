# Set the compiler and flags here
F90 = gfortran -std=f2008 -fopenmp
F77 = gfortran -std=legacy -fopenmp
# Set the build flags here
COMP = -c
#LIBS = -llapack -lblas
# List of object files
OBJS = vtmop.o vtmop_func.o linear_shepard.o svtdirect.o bvtdirect.o \
	shared_modules.o delsparse.o qnstop.o slatec.o lapack.o blas.o

# Execute the entire project
all : samples samplep
	./samples
	export OMP_NESTED=TRUE && export OMP_NUM_THREADS=2,2 && ./samplep

samples : samples.f90 $(OBJS)
	$(F90) samples.f90 $(OBJS) $(LIBS) -o samples

samplep : samplep.f90 $(OBJS)
	$(F90) samplep.f90 $(OBJS) $(LIBS) -o samplep

vtmop.o : vtmop.f90 delsparse.o linear_shepard.o svtdirect.o bvtdirect.o \
	  qnstop.o
	$(F90) $(COMP) vtmop.f90 -o vtmop.o

vtmop_func.o: vtmop_func.f90 shared_modules.o
	$(F90) $(COMP) vtmop_func.f90 -o vtmop_func.o

delsparse.o : delsparse.f90 shared_modules.o
	$(F90) $(COMP) delsparse.f90 -o delsparse.o

linear_shepard.o : linear_shepard.f90 shared_modules.o
	$(F90) $(COMP) linear_shepard.f90 -o linear_shepard.o

qnstop.o : qnstop.f90 shared_modules.o
	$(F90) $(COMP) qnstop.f90 -o qnstop.o

svtdirect.o : sVTdirect.f90 shared_modules.o
	$(F90) $(COMP) sVTdirect.f90 -o svtdirect.o

bvtdirect.o : bVTdirect.f90 shared_modules.o
	$(F90) $(COMP) bVTdirect.f90 -o bvtdirect.o

shared_modules.o : shared_modules.f90
	$(F90) $(COMP) shared_modules.f90 -o shared_modules.o

slatec.o : slatec.f
	$(F77) $(COMP) slatec.f -o slatec.o

lapack.o : lapack.f
	$(F77) $(COMP) lapack.f -o lapack.o

blas.o : blas.f
	$(F77) $(COMP) blas.f -o blas.o

clean :
	rm -f *.o *.mod samples samplep vtmop.dat vtmop.chkpt
