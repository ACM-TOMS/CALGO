########################################################################
# Options for Solaris 10, SUN4                                         #
########################################################################

# The machine (platform) identifier appended to the library names ######

PLAT      = _SOL10_SUN4

# Fortran compiler #####################################################

# Use "psrinfo -v" to get processor information

# Libraries defined by FCLDFLGS (below):
# -lf77compat : required when mixing f77 and f90 
# -lm         : common elementary mathematical functions
# -lsunmath   : Sun value-added math library
# -lmvec      : Vector math function library
# -lcpc       : CPU performance counters, required by PAPI

PROF      = 
FC        = f90
FCFLAGS   = -O5 -xarch=v8plusb -dalign -ftrap=%none -u $(PROF)
FCFLAGS   = -g -xarch=v8plusb -dalign -ftrap=%none -u $(PROF)
FCMODLS   = -M./LIB -moddir=./LIB
FCFRMT    = -free
FCLDFLGS  = -lm -lsunmath -lcpc -lmvec -xarch=v8plusb -dalign $(PROF) 
FCNOOPT   = -O0 -u -ftrap=%none $(PROF)

# LAPACK and BLAS libraries to be linked to ############################

LAPACKDIR = $(HOME)/LIB/LAPACK-3.1.1
LAPACKLIB = $(LAPACKDIR)/lapack$(PLAT).a
TMGLIB    = $(LAPACKDIR)/tmglib$(PLAT).a

BLASLIB   = -xlic_lib=sunperf # Sun BLAS library #

# Location of the STETESTER library ####################################

STETESTERDIR = $(HOME)/Projects/LAPACK/stetester
STETESTERLIB = $(STETESTERDIR)/LIB/stetesterlib$(PLAT).a

# Interface for the timing function and flops counting #################

# F77T: uses DSECND/SECOND (Fortran 77)
# F95T: uses CPU_TIME (Fortran 90/95)
# PAPI: uses PAPI's PAPIF_GET_REAL_USEC

TIMER   = F95T

# NONE: flop counting is disabled
# PAPI: uses PAPI's PAPIF_FLOPS

FLOPS   = NONE

# If PAPI is used set PAPIINC and PAPILIB accordingly ##################

PAPIINC =
PAPILIB =

# The archiver and flag(s) to use when building archive (library) ######

ARCH      = ar
ARCHFLAGS = cr
RANLIB    = touch
