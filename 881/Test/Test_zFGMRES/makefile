# LAPTOP DELL - GNU
FORTRAN   = f77
FOPTS =  -u -C
LOPTS =
BLAS   = -L/export/home/giraud/Lib/Blas -lcblas -lf77blas -latlas -lf2c
LAPACK = /usr/lib/liblapack.a
OPTS    = -O
# Laptop - PGI
FORTRAN   = mpif90
FOPTS =
LOPTS =
BLAS   = -lblas
LAPACK = -llapack
OTHER =
#
# SUN
#FORTRAN   = f77
#FOPTS =  -u -g
#LOPTS = -dalign -g
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
GMRES = zPackgmres.o zPackfgmres.o blas_rot.o
#
.SUFFIXES: .f .o
.SUFFIXES: .F .f
#
all: zTestfgmres zDriver zParDriver
#
#
zDriver: zDriver.o $(GMRES) 
	$(FORTRAN) $(LOPTS) -o $@ $@.o $(GMRES) $(LAPACK) $(BLAS) $(OTHER)
#
zParDriver: zParDriver.o $(GMRES) 
	$(FORTRAN) $(LOPTS) -o $@ $@.o $(GMRES) $(LAPACK) $(BLAS) $(OTHER)
#
zTestfgmres: zTestfgmres.o $(GMRES) 
	$(FORTRAN) $(LOPTS) -o $@ $@.o $(GMRES) $(LAPACK) $(BLAS) $(OTHER)
#
.f.o :
	$(FORTRAN) $(FOPTS) -c $*.f
clean :
	/bin/rm -f *.o zTestfgmres fort.20 fort.30 sol_Testfgmres zDriver zParDriver
