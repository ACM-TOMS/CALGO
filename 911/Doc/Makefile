# Makefile.inc sets the compiler and flags to be used
# in all other Makefiles
#
# F95 needs to be set to the Fortran compiler you wish to use
# F95FLAGS needs to be set to the required compilation flags
#          We recommend you switch on as much checking as possible
#          when you first try this out

F95=sunf95
F95FLAGS=

F95=nagfor
F95FLAGS=-C=all -C=undefined -nan -gline -f2003  -g90 -u
F95FLAGS=-O3

F95=gfortran
F95FLAGS=

# LINLIB is the Linpack library
LINLIB=

%.o : %.f95
	$(F95) -c $(F95FLAGS) $< -o $@

# assume the compiler treats .f files as fixed 
%.o : %.f
	$(F95) -c $(F95FLAGS) $< -o $@

# turn off the rule that attempts (on some ancient systems) to
# compile .mod files as if they contained Modula code!
%.o : %.mod

all: TestFM SampleFM

TestFM:  fmsave.o  FM.o  FMZM90.o  TestFM.o
	$(F95) $(F95FLAGS) -o $@ $^

SampleFM:  fmsave.o  FM.o  FMZM90.o  SampleFM.o
	$(F95) $(F95FLAGS) -o $@ $^

fmsave.o : fmsave.f95

FM.o: FM.f95 fmsave.o

FMZM90.o: FMZM90.f95 FM.o fmsave.o

clean:
	rm -rf *.o *.mod *.g90 TestFM SampleFM TEMPFM FMERRMSG.OUT
