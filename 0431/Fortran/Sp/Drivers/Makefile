all: Res

src.o: src.f
	$(F77) $(F77OPTS) -c src.f

DRIVERS= driver
RESULTS= Res

Objs1= src.o
driver: $(Objs1)
	$(F77) $(F77OPTS) -o driver $(Objs1) $(SRCLIBS)
Res: driver  data
	./driver < data >Res

diffres:Res res
	echo "Differences in results from driver"
	$(DIFF) Res res

clean: 
	rm -rf *.o $(DRIVERS) $(CLEANUP) $(RESULTS)
