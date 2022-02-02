include makefile.inc

all: Res

src.o: $(INCSRC)/src.f
	$(F77) $(F77OPTS) -c $(INCSRC)/src.f

DRIVERS= driver1
RESULTS= Res

Libs1= $(PORT)
Objs1= driver.o src.o
driver1: $(Objs1)
	$(F77LINK) $(F77LINKOPTS) -o driver1 $(Objs1) $(Libs1)
Res: driver1 data 
	./driver1 <data >Res

diffres:Res res
	echo "Differences in results from driver"
	$(DIFF) Res res

clean: 
	rm -rf *.o $(DRIVERS) $(CLEANUP) $(RESULTS)
