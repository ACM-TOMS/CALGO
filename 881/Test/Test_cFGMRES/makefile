# Laptop - PGI
FORTRAN   = mpif90
FOPTS =
LOPTS =
BLAS   = -lblas
LAPACK = -llapack
OTHER = 
#
#
# SUN
#FORTRAN   = f77
#FOPTS =  -u -g
#LOPTS = -dalign -g
#BLAS   = /usr/local/lib/libblas.a
#LAPACK = /usr/local/lib/liblapack.a
#OTHER = blas_rot.o
#
# IBM
#FORTRAN   = xlf
#FOPTS =  -u
#LOPTS = 
#BLAS   = /usr/local/lib/libblas.a
#LAPACK = /usr/local/lib/liblapack.a
#OTHER = blas_rot.o
#
# SGI
#FORTRAN = f77
#FOPTS =  -u
#LOPTS = 
#BLAS   = 
#LAPACK = -lcomplib.sgimath
#OTHER  = blas_rot.o
#
# user's object files
#
GMRES = cPackgmres.o cPackfgmres.o blas_rot.o 
#
.SUFFIXES: .f .o
.SUFFIXES: .F .f
#
all: cTestfgmres cDriver cParDriver
#
#
cParDriver: cParDriver.o $(GMRES) 
	$(FORTRAN) $(LOPTS) -o $@ $@.o $(GMRES) $(LAPACK) $(BLAS) $(OTHER)
#
cDriver: cDriver.o $(GMRES)  $(OTHER)
	$(FORTRAN) $(LOPTS) -o $@ $@.o $(GMRES) $(OTHER) $(LAPACK) $(BLAS)
#
cTestfgmres: cTestfgmres.o $(GMRES) $(OTHER)
	$(FORTRAN) $(LOPTS) -o $@ $@.o $(GMRES) $(OTHER) $(LAPACK) $(BLAS) 
#
.f.o :
	$(FORTRAN) $(FOPTS) -c $*.f
clean :
	/bin/rm -f *.o cTestfgmres fort.* sol_Testfgmres cDriver cParDriver
