FORTRAN   = f77
FOPTS =  -u -C
LOPTS =
BLAS   = -L/export/home/giraud/Lib/Blas -lcblas -lf77blas -latlas -lf2c
LAPACK = /usr/lib/liblapack.a
OPTS    = -O
#
# Mac OS-X F77
FORTRAN = mpif77
CC = cc
OPTS    = -u
BLAS   = -L/sw/lib/ -lcblas -lf77blas -latlas -lptf77blas -lg2c
LAPACK = -llapack
OTHER  = -lmpi
#
# Mac OS-X F90
FORTRAN = mpif90
FOPTS =  
CC = mpicc
OPTS    = 
BLAS   = -L/sw/lib/ -lcblas -lf77blas -latlas -lg2c
LAPACK = -llapack
MPI  = 
# Laptop - PGI
FORTRAN   = mpif90
FOPTS = 
LOPTS =
BLAS   =  -lblas
LAPACK =  -llapack 
OTHER =
MPI  = 
# SUN
#FORTRAN   = f77
#FOPTS =  -u -O
#LOPTS = -dalign -O
#BLAS   = /usr/local/lib/libblas.a
#LAPACK = /usr/local/lib/liblapack.a
#OTHER =
#
# IBM
#FORTRAN   = xlf
#FOPTS =  -u
#LOPTS = 
#BLAS   = /usr/local/lib/libblas.a
#LAPACK = /usr/local/lib/liblapack.a
#OTHER =
#
# DEC
#FORTRAN   = f77
#FOPTS =  -u
#LOPTS = 
#BLAS   = /usr/lib/libblas.a
#LAPACK = /usr/lib/liblapack.a
#OTHER = -lfor -lm -lots -lm_4sqrt
#
# CRAY - T3D
#FORTRAN = TARGET=cray-t3d f90
#FOPTS = -O2 -dp -I/usr/include/mpp
#LOPTS = -O2 -dp -I/usr/include/mpp
#BLAS   = /opt/ctl/craylibs_m/craylibs_m/libsci.a
#LAPACK = 
#
# SGI
#FORTRAN = f77
#FOPTS =  -u
#LOPTS = 
#BLAS   = 
#LAPACK =
#OTHER  = -lcomplib.sgimath
#
# user's object files
#
GMRES = dPackgmres.o dPackfgmres.o
#
.SUFFIXES: .f .o
.SUFFIXES: .F .f
#
all: dTestfgmres dDriver dParDriver
#
#
dDriver: dDriver.o $(GMRES) 
	$(FORTRAN) $(LOPTS) -o $@ $@.o $(GMRES) $(LAPACK) $(BLAS) 
#
dParDriver: dParDriver.o $(GMRES) 
	$(FORTRAN) $(LOPTS) -o $@ $@.o $(GMRES) $(LAPACK) $(BLAS) $(MPI)
#
dTestfgmres: dTestfgmres.o $(GMRES) 
	$(FORTRAN) $(LOPTS) -o $@ $@.o $(GMRES) $(LAPACK) $(BLAS) 
#
.f.o :
	$(FORTRAN) $(FOPTS) -c $*.f
clean :
	/bin/rm -f *.o dTestfgmres fort.* sol_Testfgmres dDriver dParDriver dTestfgmres
