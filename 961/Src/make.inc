####################################################################
#  SHHEIG make include file.                                       #
#  SHHEIG, Release 1.0                          ./shheig/make.inc  #
#  Matthias Voigt, MPI Magdeburg                                   #
#  July 23, 2013.                                                  #
####################################################################
#
#  The machine (platform) identifier to append to the library names
#
#  PLAT = _sun4
#  
#  Modify the FORTRAN and OPTS definitions to refer to the
#  compiler and desired compiler options for your machine.  NOOPT
#  refers to the compiler options desired when NO OPTIMIZATION is
#  selected.  Define LOADER and LOADOPTS to refer to the loader and 
#  desired load options for your machine.
#
#  IMPORTANT: Use the options -fPIC -fdefault-integer-8 for 64bit
#             architectures.
#
FORTRAN  = gfortran
OPTS     = -O2 -fPIC -fdefault-integer-8
NOOPT    = 
LOADER   = gfortran
LOADOPTS = $(SHHEIGLIB) $(SLICOTLIB) $(LPKAUXLIB) $(LAPACKLIB) \
	   $(BLASLIB)
#
#  The archiver and the flag(s) to use when building archive (library)
#  If your system has no ranlib, set RANLIB = echo.
#
ARCH     = ar
ARCHFLAGS= r
# ARCHFLAGS= cr
# RANLIB   = ranlib
#
#  The location of the libraries to which you will link.  (The 
#  machine-specific, optimized BLAS library should be used whenever
#  possible.)
#
BLASLIB      = ../libs/blas64.a
LAPACKLIB    = ../libs/lapack64.a
SLICOTLIB    = ../libs/slicot64.a
LPKAUXLIB    = ../libs/lpkaux64.a
SHHEIGLIB    = ../shheig64.a
