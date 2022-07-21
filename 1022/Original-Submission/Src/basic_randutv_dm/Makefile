###############################################################################
#
# Makefile
#
###############################################################################

#
# Define the Fortran90 compiler and loader, and their flags.
#
F90           = mpif90
F90FLAGS      = -O3 -Wall

F90LOADER     = $(F90)
F90LOADFLAGS  =

#
# Define the libraries.
#
MKLROOT       = /opt/intel/compilers_and_libraries_2018.2.199/linux/mkl
LIBS          = -L${MKLROOT}/lib/intel64 \
                -Wl,--no-as-needed \
                -lmkl_scalapack_lp64 \
                -lmkl_gf_lp64 -lmkl_gnu_thread -lmkl_core \
                -lmkl_blacs_intelmpi_lp64 -lgomp -lpthread -lm -ldl

#
# Define objects and dependencies.
#
OBJECTS = pdttrandutv.o \
          pdgeutv.o

all : pdttrandutv.x

pdttrandutv.x : $(OBJECTS)
	$(F90LOADER) $(F90LOADFLAGS) -o $@  \
            $(OBJECTS) $(LIBS)

pdttrandutv.o : pdttrandutv.f90
	$(F90) -c $(F90FLAGS) $<

pdgeutv.o : pdgeutv.f90
	$(F90) -c $(F90FLAGS) $<

clean :
	- rm -f *.o *.x core

.f90.o : 
	$(F90) -c $(F90FLAGS) $*.f

