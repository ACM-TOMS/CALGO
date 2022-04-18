PROG =	testmarq

SRCS =	cons.f90 incgam.f90 marcumq.f90 testmarq.f90

OBJS =	cons.o incgam.o marcumq.o testmarq.o

FFLAGS = -O
FC = g95
LDFLAGS =

all: $(PROG)

$(PROG): $(OBJS)
	$(FC) $(LDFLAGS) -o $@ $(OBJS)

clean:
	rm -f $(PROG) $(OBJS) *.mod

.SUFFIXES: $(SUFFIXES) .f90

%.o:%.f90
	$(FC) $(FFLAGS) -c $<

cons.o: constants cons.f90
incgam.o: cons.o incgam.f90
marcumq.o: cons.o incgam.o marcumq.f90
testmarq.o: cons.o marcumq.o testmarq.f90
