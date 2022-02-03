include ../make.gbinc

### Timing library #####################################################
all: $(CTMG)

TIMELIB = $(CTMG)

CSECS = second.f dsecnd.f
CSEC = second.o dsecnd.o

$(CTMG):	$(CSEC)
	$(ARCH) $(ARCHFLAGS) $(TIMELIB) $(CSEC)
	$(RANLIB) $(TIMELIB)

$(CSEC):	$(CSECS)
	$(FORTRAN) -c $(TMGOPT) $(CSECS)

clean:
	rm -f *.o
