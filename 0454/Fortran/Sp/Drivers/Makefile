all: Res

src.o: src.f
	$(F77) $(F77OPTS) -c src.f

driver1.o: driver1.f
	$(F77) $(F77OPTS) -c driver1.f

driver2.o: driver2.f
	$(F77) $(F77OPTS) -c driver2.f

driver3.o: driver3.f
	$(F77) $(F77OPTS) -c driver3.f

DRIVERS= driver1 driver2 driver3
RESULTS= Res1 Res2 Res3

Objs1= driver1.o src.o
Objs2= driver2.o src.o
Objs3= driver3.o src.o

driver1: $(Objs1)
	$(F77) $(F77OPTS) -o driver1 $(Objs1) $(SRCLIBS)

driver2: $(Objs2)
	$(F77) $(F77OPTS) -o driver2 $(Objs2) $(SRCLIBS)

driver3: $(Objs3)
	$(F77) $(F77OPTS) -o driver3 $(Objs3) $(SRCLIBS)

Res: driver1 driver2 driver3 
	./driver1 >Res1
	./driver2 >Res2
	./driver3 >Res3

diffres:Res1 Res2 Res3 res1 res2 res3
	echo "Differences in results from driver"
	$(DIFF) Res1 res1
	$(DIFF) Res2 res2
	$(DIFF) Res3 res3

clean: 
	rm -rf *.o $(DRIVERS) $(CLEANUP) $(RESULTS)
