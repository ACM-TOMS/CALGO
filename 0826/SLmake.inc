FORTRAN = f90
OPTS = -O3
DRVOPTS  = $(OPTS)
NOOPT    =
LOADER   = f90
LOADOPTS =

CC = cc
CCFLAGS = $(OPTS)
CDEFS = -DAdd_ $(USEMPI)

ARCH     = ar
ARCHFLAGS= cr
RANLIB   = echo

MACH = origin
PLAT = _irix65

SCALAPACK = /usr/local/usp/PETtools/lib/libscalapackn32.a
BLACSMPI = /usr/local/usp/PETtools/lib/libmpiblacsn32.a

BLASLIB     = -lblas
LAPACKLIB    = $(HOME)/LAPACK/lapack$(PLAT).a
MATGENLIB    = $(HOME)/LAPACK/tmglib$(PLAT).a
PXLAHQRLIB    = $(HOME)/$(MACH)/ScaLAPACK/pxlahqrlib.a

SERLIBS = $(MATGENLIB) $(LAPACKLIB) $(BLASLIB)
PARLIBS = $(PXLAHQRLIB) $(SCALAPACK) $(BLACSMPI) $(SERLIBS) -lscs -lmpi

