PROG =	testgamstar

SRCS =	setprecision.f90 cons.f90  gammaError.f90 incgamNEG.f90 testgamstar.f90

OBJS =	setprecision.o cons.o  gammaError.o incgamNEG.o testgamstar.o

FFLAGS = -O
FC = gfortran
LDFLAGS =

all: $(PROG)

$(PROG): $(OBJS)
	$(FC) $(LDFLAGS) -o $@ $(OBJS)

clean:
	rm -f $(PROG) $(OBJS) *.mod

.SUFFIXES: $(SUFFIXES) .f90

%.o:%.f90
	$(FC) $(FFLAGS) -c $<

setprecision.o: setprecision.f90
cons.o: setprecision.o constants cons.f90
gammaError.o: setprecision.o cons.o gammaError.f90
incgamNEG.o: setprecision.o cons.o incgamNEG.f90
testgamstar.o: setprecision.o cons.o  gammaError.o incgamNEG.o testgamstar.f90
