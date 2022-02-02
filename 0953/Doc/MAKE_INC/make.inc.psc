#
# make.inc
#
# This file needs to be modified according to your computational environment.
#

# Settings for the Pathscale compiler
#   module add psc
#   module add openmpi/psc

#
# Compiler and linker
#
FC      = mpif90
OPENMP  =
#OPENMP  = -mp -DUSE_OMP
FFLAGS  = -O3 -fPIC -C -cpp $(OPENMP)
LOADER  = $(FC)
LDFLAGS = $(FFLAGS)
AR      = ar
ARFLAGS = cr
RANLIB  = ranlib

#
# External libraries
#
BLASLIB = $(GOTO2_LDFLAGS) $(GOTO2_LIBS)
LAPACKLIB = $(LAPACK_LDFLAGS) $(LAPACK_LIBS)
SCALAPACKLIB = $(HOME)/pfs/Akka/lib/libscalapack.a
LIBS    = $(SCALAPACKLIB) $(LAPACKLIB) $(BLASLIB)
