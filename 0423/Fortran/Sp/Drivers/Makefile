all: Res

src.o: src.f
	$(F77) $(F77OPTS) -c src.f

driver.o: driver.f
	$(F77) $(F77OPTS) -c driver.f

DRIVERS= driver1
RESULTS= Res

Objs1= driver.o src.o
driver1: $(Objs1)
	$(F77) $(F77OPTS) -o driver1 $(Objs1) $(SRCLIBS)
Res: driver1 
	./driver1 >Res

diffres:Res res
	echo "Differences in results from driver"
	$(DIFF) Res res

clean: 
	rm -rf *.o $(DRIVERS) $(CLEANUP) $(RESULTS)
