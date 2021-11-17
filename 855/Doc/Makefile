include makefile.inc

OBJS = driver.o mcnr.o mcnf.o

all: res

res: $(OBJS) 
	$(CPP) $(CPPOPTS) -o driver $(OBJS)
	driver  > res

driver.o: driver.cpp mcnf.cpp mcnr.cpp
	$(CPP) $(CPPOPTS) -c  driver.cpp

mcnr.o: mcnr.cpp
	$(CPP) $(CPPOPTS) -c mcnr.cpp

mcnf.o: mcnf.cpp
	$(CPP) $(CPPOPTS) -c mcnf.cpp
