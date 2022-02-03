# Suffix for Fortran 90 programs
.SUFFIXES:
.SUFFIXES : .f90 $(SUFFIXES)
#-----------------------------------------------------------------------

# Name of Fortran 90 compiler
FC = f95

#  Select your compiler flags
#  Example for SUN's f95
#CFLAGS = -M../Src -fast
#  Example for NAGWare f95
CFLAGS = -ieee=full -nan -C=all -I../Src
#CFLAGS = -I../Src

# Suffix for Modules (created by compiler)
MSUFF = mod

#-----------------------------------------------------------------------

# Select main program (do not give a suffix!)
# With this distribution the following test programs are distributed:
EXAMPLES= details ex_qag ex_qags ex_triex ex_decuhr2d ex_cutet ex_decuhr3d simplexpapertest
# details    : prints info on current distribution
# ex_qag     : 1-dimensional integration, tests of QAG from QUADPACK
# ex_qags    : 1-dimensional integration, tests of QAGS from QUADPACK
# ex_triex   : 2-dimensional integration (triangle), tests of TRIEX
# ex_decuhr2d: 2-dimensional integration (square), test of DECUHR
# ex_cutet   : 3-dimensional integration (tetrahedron), test of DCUTET
# ex_decuhr3d: 3-dimensional integration (cube), tests of DECUHR
# simplexpapertest : 5-dimension integration (simplex)

MAIN = simplexpapertest

#-----------------------------------------------------------------------

# The target (executable) will have the same name as the main
# program, without suffix.
TARGET = $(MAIN)

#-----------------------------------------------------------------------
.f90.o:
	$(FC) -c $(CFLAGS) $<

#-----------------------------------------------------------------------

all: $(MAIN).o
	$(FC) -o $(TARGET) $(MAIN).o $(libcui)

#-----------------------------------------------------------------------


# Place library in current directory (or change the -L../Src)
libcui = -L../Src -lcubpack


clean:
	/bin/rm -f *.o *.$(MSUFF) $(EXAMPLES)
