FFLAGS= -sloppy -C -d5 -g -temp=/tmp -u
FC=epcf90
VPATH = .

OBJS = driver.o gmpsg.o

planar : $(OBJS)
	$(FC) $(FFLAGS) $(OBJS) -o $(VPATH)/planar
