# Makefile for F90 version of MINRESQLP for symmetric systems.
# To run this file:
#    make
# Or to remove object files and executable:
#    make clean
#
# Authors:
#     Sou-Cheng Choi <sctchoi@uchicago.edu>
#     Computation Institute (CI)
#     University of Chicago
#     Chicago, IL 60637, USA
#
#     Michael Saunders <saunders@stanford.edu>
#     Systems Optimization Laboratory (SOL)
#     Stanford University
#     Stanford, CA 94305-4026, USA
#
# 26 Aug 2012: First release version.
# 15 Oct 2007: Don't use compiler option -r8 to get double precision
#              because it is nonstandard.
#              Use minresqlpDataModule to define real(kind=dp).
# 11 Oct 2007: First version for compiling minresqlpTestProgram.f90
#              and associated modules.
#              All files listed are .f90 source code.
#              C is not currently used.

# Set exactly one of these to yes
# (Extend the list if necessary)

  USEg95        = no
  USEgfortran   = yes
  USEgeneric90  = no
  USEgeneric    = no
  USEnag        = no

ifeq ($(USEg95),yes)
  FC      =  g95
  FFLAGS1 = -g -O0 -pedantic -Wall -Wextra -fbounds-check -ftrace=full
  FFLAGS2 = -g -O
endif

ifeq ($(USEgfortran),yes)
  FC      =  gfortran
  FFLAGS1 = -g -O0 -pedantic -Wall -W -fbounds-check
  FFLAGS2 = -g -O
endif

ifeq ($(USEgeneric),yes)
  FC      =  f95
  FFLAGS1 = -g -O0
  FFLAGS2 = -g -O
endif

ifeq ($(USEgeneric90),yes)
  FC      =  f90
  FFLAGS1 = -g -O0
  FFLAGS2 = -g -O
endif

ifeq ($(USEnag),yes)       
  FC      =  nagfor
  FFLAGS1 = -C=all -C=undefined -nan -gline -f2003  -g90 -u -kind=byte
  FFLAGS2 = -nan -gline -f2003  -g90 -u -kind=byte
endif

# Select one of these
#FFLAGS  = ${FFLAGS1}    # for development
 FFLAGS  = ${FFLAGS2}    # for release


  CC      =  gcc           # Not currently used
  CFLAGS  = -g -O

# Clear suffix list, then define the ones we want
  .SUFFIXES:
  .SUFFIXES: .c .f .f90 .o

  .f90.o:; ${FC} ${FFLAGS} -c -o $@ $<
  .f.o:;   ${FC} ${FFLAGS} -c -o $@ $<
  .c.o:;   $(CC) $(CFLAGS) -c -o $@ $<

  files = minresqlpDataModule.o minresqlpBlasModule.o minresqlpModule.o mm_ioModule.o minresqlpReadMtxModule.o minresqlpTestModule.o minresqlpTestProgram.o

minresqlptest: ${files}
	${FC} ${FFLAGS} -o $@ ${files}

clean:
	\rm -f *.o *.mod *.exe minresqlptest minresqlp.mexmaci64 
