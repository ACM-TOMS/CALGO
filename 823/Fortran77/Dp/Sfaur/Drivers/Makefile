include makefile.inc

all: Res1

SRCLIBS= $(LINPACK) $(EISPACK_DP) $(BLAS) $(PORT)

src.o: $(INCSRC)/src.f
	$(F77) $(F77OPTS) -c $(INCSRC)/src.f

extra.o: $(INCSRC)/extra.f
	$(F77) $(F77OPTS) -c $(INCSRC)/extra.f

DRIVERS= driver
RESULTS= Res1
SRCLIBS = $(TIMER)

Objs1= driver.o src.o extra.o
driver: $(Objs1)
	$(F77LINK) $(F77LINKOPTS) -o driver $(Objs1) $(SRCLIBS)
Res1: driver 
	./driver >Res1

diffres: Res1 res1
	echo "Differences in results from driver"
	$(DIFF)  Res1 res1

clean: 
	rm -rf *.o $(DRIVERS) $(CLEANUP) $(RESULTS)
