FFLAGS= -sloppy -C -d5 -g -temp=/tmp -u
FC=epcf90

OBJS = d1mach.o lasys.o rfft_pack.o

neuman : $(OBJS) neuman.o driver1.o
	$(FC) $(FFLAGS) $(OBJS) neuman.o driver1.o -o driver1
	driver1 < data1
drchlt : $(OBJS) drchlt.o driver2.o
	$(FC) $(FFLAGS) $(OBJS) drchlt.o driver2.o -o driver2
	driver2 < data2 
