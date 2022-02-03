# Time-stamp:  <2016-09-01 11:44:06 m>
# # # This make file by Fred T. Krogh, benefiting from lots of testing help from
# # # Philip W. Sharp and Richard J. Hanson.  At the suggestion of a referee
# # # parts of the old Fortran and C makefiles have been put here.

# # # For a new installation, after unpacking the zip file in a new
# # # directory, change this make file as appropriate for you environment.
# # # Or just type the command
# # # make FC=<your Fortran compiler> CC=<Your C compiler>
# # # Running make will compile all the codes, run them and compare results
# # # make clean will clean up things
# # # make replace_orig may simplify comparisons if you change NUMT.

# For examaple set FC/CC= gfortran/gcc on linux, ifort/icc for Intel on linux,
# ftn/cc on Cray, nagfor/gcc for the NAG compiler, ifc/icc for Intel on Windows,
# and pgf90/pgcc for the Portland Group compiler.

#See notes for the Cray ftn compiler below.

#SPTEST=1 # Uncommented may enable some special testing.

# If you want to define your compiler in this make file do it just below.
# Leave CC=undefined If you don't want the C support compiled.
FC=
CC=

#OPT=-O2 # Commented gives extra checks, Uncomment and set compiler flag for the
#          optimization level you want for regular use.

# Set NUMT = 2 to test threads, set to 1 for no thread support.  (Gives a
#    slightly smaller object files for messy.)
NUMT=2
# Precisions to generate for Fortran and C: single, double, quad
PLETS=s d q
CPLETS=s d
#NOQ=1 # Uncomment this line if PLETS does not contain q.
# Subdirectory containing source files (from Fortran90 or C directories)
SRC=Src
# Subdirectory containing library files (from Fortran90 or C directories)
LIBS=../Libs
# Subdirectory with test drivers (from Fortran90 or C directories)
DRV=Drivers
# Subdirecotry where documentation is stored (from Fortran90 or C directories)
DOC=../UserManuals
# Subdirecotry for correct results are stored (from Fortran90 or C directories)
RO=$(DRV)/Results
# Subdirecotry for new results are stored (from Fortran90 or C directories)
R=$(DRV)/NewResults
# Compiler options for threads
T=-fopenmp -lgomp -Dnumt_=$(NUMT)
# NOTE: Spcifications below may be altered later for different compiler choices.
# Specify compiler options for Fortran and C libraries
LIBF= -L$(LIBS) -lmathf
LIBC= -L$(LIBS) -lmathc -lm
# Specify include directories for the main code and additional for the drivers.
INCL= -I$(LIBS)
INCLD= -I$(DRV)
# Specify module directories for the main code and additional for the drivers.
MODU=-J$(LIBS)
MODUD=-J$(DRV)

ifeq ($(FC),myfor) # Set up for Krogh's latest Fortran compiler
  std_flags=-march=native -funroll-loops -pipe \
   -ftree-vectorize -floop-strip-mine  -floop-block -fimplicit-none \
   -finit-real=NAN -Wall -Wextra -Wno-compare-reals -Wno-c-binding-type
# -ffpe-trap=zero,overflow,invalid
  ifdef OPT
     FFLAGS=$(OPT) $(std_flags)
     CFLAGS=$(OPT)
  else
     FFLAGS=-g $(std_flags) -fbounds-check
     CFLAGS=-g -fbounds-check
  endif
endif

ifeq ($(FC),gfortran)
  std_flags=-march=native -funroll-loops -pipe \
   -ftree-vectorize -floop-strip-mine  -floop-block -fimplicit-none \
   -finit-real=NAN -Wall -Wextra -Wno-compare-reals -Wno-c-binding-type
# -ffpe-trap=zero,overflow,invalid
  ifdef OPT
     FFLAGS=$(OPT) $(std_flags)
     CFLAGS=$(OPT)
  else
     FFLAGS=-g $(std_flags) -fbounds-check
     CFLAGS=-g -fbounds-check
  endif
endif
ifeq ($(FC),ifort)
  MODU=-module $(LIBS)
  MODUD=-module $(DRV)
  T=-fopenmp -lomp_db -liomp5 -Dnumt_=$(NUMT)
  ifdef OPT
     FFLAGS=$(OPT) -nofor-main
     CFLAGS=$(OPT)
  else
     FFLAGS=-g -nofor-main
     CFLAGS=-g
  endif
endif
ifeq ($(FC),itc)
  ifdef OPT
     FFLAGS=$(OPT)
     CFLAGS=$(OPT)
  else
     FFLAGS=-g -fpp
     CFLAGS=-g -fpp
  endif
  MODU=-module $(LIBS)
  MODUD=-module $(DRV)
endif

ifeq ($(FC),nagfor)
  ifdef OPT
     FFLAGS=$(OPT) -DNAG
     CFLAGS=$(OPT)
  else
     FFLAGS=-g -C=all -nan -DNAG
     CFLAGS=-g
  endif
  MODU=-mdir $(LIBS)
  MODUD=-mdir $(DRV)
#  NUMT=1 # Version 6.0 of nagfor does not support threadprivate arrays.
  T=-openmp -Dnumt_=$(NUMT)
  ifeq ($(CC),gcc)
    MODUC=-J
  endif
endif
ifndef MODUC
  MODUC=$(MODU)
endif
RUN=
RUNC=

ifeq ($(FC),pgf90)
# No support for q here.  Worked with no threads once, but no longer.
  NOQ=1
  PLETS=s d
  NUMT=1
  T=-mp -Dnumt_=$(NUMT)
  MODU=-module $(LIBS)
  MODUD=-module $(DRV)
  ENVOUT=export NO_STOP_MESSAGE=yes
  RUN:=$(ENVOUT);
endif
ifdef NOQ
  FFLAGS+= -DNOQ
endif

ifeq ($(FC),ftn)
# At one time we had a working version for the Cray. What is here has not be
# successfully tested, and likely needs some work.  In the past we need to
# precede execution with a command like
#   qsub -I -l mppwidth=2 -q small -l mppnppn=2
  MODU=
  MODID=
  ifdef OPT
    FFLAGS=$(OPT)
    CFLAGS=$(OPT)
  else
# Only a single thread can be used with -g
    FFLAGS=-O2
    CFLAGS=-O2
  endif
  ifndef OPT # Cray can only do one thread if no optimization.
    NUMT = 1
  endif
  T=-Dnumt_=$(NUMT)
  RUN=aprun -q -n1 -d1
endif
FFLAGS+= $(T)
CFLAGS+= $(T)

VERS=Fortran
ifdef CC
VERS+= CCODE
endif

all : $(VERS)

.PHONY : all
.EXPORT_ALL_VARIABLES:
Fortran:
	$(MAKE) -C Fortran90 -e

CCODE:
	$(MAKE) -C C -e

clean: # Used prior to zipping or to check after major changes.
	$(MAKE) -C Fortran90 -e clean;\
  $(MAKE) -C C -e clean;

replace_orig:
	$(MAKE) -C Fortran90 -e replace_orig;\
  $(MAKE) -C C -e replace_orig
