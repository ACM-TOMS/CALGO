F95 = gfortran

F77OPTS = -O 
F95OPTS = $(F77OPTS)

# Switch off all suffix rules
.SUFFIXES:

VPATH = ../Src

# Generate .o files from .f95 or .f90
%.o: %.f90
	$(F95) $(F95OPTS) -o $@ -c $<

%.o: %.f
	$(F95) $(F77OPTS) -o $@ -c $<

# Turn off the stupid default rules for Modula 2
# as some compilers generate .mod files
%.o: %.mod

# generate executable depending on .o files
%: %.o
	$(F95) $(F95LINKOPTS) -o $@ $^ $(LIBS)

# Precompiled libraries for linpack and the blas could be used here
OBJS = routine.o linpack.o blas.o timer.o

driver1 driver2 driver3 driver4 driver5 driver6: $(OBJS)

all :  driver1 driver2 driver3 driver4 driver5 driver6
