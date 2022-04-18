##########################################################################
#
#  Modify the FORTRAN and _OPT definitions to refer to the compiler
#  and desired compiler options for your machine. Define LOADER and
#  LOADOPTS to refer to the loader and desired load options for your
#  machine. The compiler options are specified as follows for the
#  different programs and libraries:
#
#     GBL3BOPT : the GEMM-based level 3 BLAS routines
#     GBBENOPT : the GEMM-based performance benchmark programs
#     TSTOPT   : the test programs for correctness tests
#     UNDEROPT : the underlying BLAS routines
#     TMGOPT   : the timing routines
#     AUXOPT   : various auxiliary routines
#
FORTRAN   = f77
GBL3BOPT  = -O2
GBBENOPT  = -O2
TSTOPT    = -O2
UNDEROPT  = -O2
TMGOPT    = -O2
AUXOPT    = -O2
LOADER    = f77
LOADOPT   =
#
#  The archiver and the flag(s) to use when building archive (library)
#  If you system has no ranlib, set RANLIB = echo.
#
ARCH       = ar
ARCHFLAGS  = cr
RANLIB     = ranlib
#
#  Libraries to be created. If the libraries are not specified, the
#  libraries are not created. The libraries are specified relative
#  to the subdirectories where they are produced.
#
#     CREGBLIB  : the GEMM-based level 3 BLAS library
#     CREULIB   : a library for the underlying BLAS routines
#     CRETMG    : a library for the timing routines
#
CREGBLIB  = ../../libgbl3b.a
CREULIB   = ../libunderlib.a
CRETMG    = ../../libtmglib.a
#
#  Libraries used by the testing and timing programs.
#
#     USELIB    : user specified level 3 BLAS library to be tested
#                 and/or timed.
#     USEULIB   : a library for the underlying BLAS routines
#     USETMG    : a library for the timing routines
#
USELIB    = ../../libgbl3b.a
USEULIB   = ../../libunderlib.a
USETMG    = ../../libtmglib.a
