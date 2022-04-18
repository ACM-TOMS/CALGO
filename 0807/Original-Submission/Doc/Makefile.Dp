include makefile.inc
all: res1 res2 res3 res4
SRCLIBS= $(PORT) $(BLAS)
DRIVERLIBS = $(TIMER) $(LAPACK) 
Objs1= driver.o src.o

driver1: $(Objs1)
	$(F77LINK) $(F77LINKOPTS) -o driver1 $(Objs1) $(DRIVERLIBS) $(SRCLIBS)
res1: driver1 data1 
	driver1 <data1 >res1
res2: driver1 data2 
	driver1 <data2 >res2
res3: driver1 data3 
	driver1 <data3 >res3
res4: driver1 data4 
	driver1 <data4 >res4
