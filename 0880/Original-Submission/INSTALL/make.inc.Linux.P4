########################################################################
# Options for Linux, Pentium 4 (Intel Fortran compiler)                #
########################################################################

# The machine (platform) identifier appended to the library names ######

PLAT      = _LINUX_IFORT

# Fortran compiler #####################################################

FC        = ifort
FCFLAGS   = -O3 -tpp7 -arch pn4 -mp -module ./LIB
FCFRMT    = -FR 
FCLDFLGS  = -O3 -tpp7 -arch pn4 -mp
FCNOOPT   = -O0 -fp_port

# LAPACK and BLAS libraries to be linked to ############################

LAPACKDIR = $(HOME)/LIB/LAPACK
LAPACKLIB = $(LAPACKDIR)/lapack$(PLAT).a
TMGLIB    = $(LAPACKDIR)/tmglib$(PLAT).a

BLASLIB   = -L$(HOME)/LIB/intel/mkl/8.0.2/lib/32 -lmkl_ia32 -lguide -lpthread 

# Location of the STETESTER library ####################################

STETESTERDIR = $(HOME)/Projects/LAPACK/stetester
STETESTERLIB = $(STETESTERDIR)/LIB/stetesterlib$(PLAT).a

# Interface for the timing function and flops counting #################

# F77T: uses DSECND/SECOND (Fortran 77)
# F95T: uses CPU_TIME (Fortran 90/95)
# PAPI: uses PAPI's PAPIF_GET_REAL_USEC

TIMER   = F77T

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
