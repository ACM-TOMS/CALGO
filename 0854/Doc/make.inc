#######################################################################
#
# The user must set this options:
#
#  Modify the FORTRAN and OPTS definitions to refer to the
#  compiler and desired compiler options for your machine. 
#
FORTRAN  = f77
OPTS     = -O4
LOADER   = f77
LOADOPTS = $(HAPACKLIB) $(LAPACKLIB) $(BLASLIB)
#
#  The archiver and the flag(s) to use when building archive (library)
#  If your system has no ranlib, set RANLIB = echo.
#
ARCH      = ar
ARCHFLAGS = cr
RANLIB    = ranlib
#
#  The location of the libraries to which you will link.  (The 
#  machine-specific, optimized BLAS library should be used whenever
#  possible.)
#
BLASLIB      = -L/usr/site-local/lib -lblas
LAPACKLIB    = -L/usr/site-local/lib -llapack
HAPACKLIB    = ../hapack.a
