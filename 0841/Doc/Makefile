INCSRC = ../Src

src.o: $(INCSRC)/src.f
	$(F77) $(F77OPTS) -c $(INCSRC)/src.f
driver2.o: driver2.f
	$(F77) $(F77OPTS) -c driver2.f
driver1.o: driver1.f
	$(F77) $(F77OPTS) -c driver1.f

DRIVERS= driver1 driver2

Libs1= $(BLAS)

Objs1= driver1.o src.o
Objs2= driver2.o src.o

all: res1 res2

driver1: $(Objs1)
	$(F77) $(F77OPTS) -o driver1 $(Objs1) $(Libs1)
res1: driver1 data
	driver1 < data > res1

driver2: $(Objs2)
	$(F77) $(F77OPTS) -o driver2 $(Objs2) $(Libs1)
res2: driver2 data
	driver2 < data > res2

diffres:Res1 res1 Res2 res2
	echo "Differences in results from driver1"
	$(DIFF) Res1 res1
	echo "Differences in results from driver2"
	$(DIFF) Res2 res2

clean: 
	rm -rf *.o $(DRIVERS) $(CLEANUP) $(RESULTS)
