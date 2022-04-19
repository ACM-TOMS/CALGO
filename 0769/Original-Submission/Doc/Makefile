FFLAGS= -sloppy -C -d5 -g -temp=/tmp -u
FC=epcf90
VPATH = .

OBJS = driver.o src.o

gqaps : $(OBJS)
	$(FC) $(FFLAGS) $(OBJS) -o $(VPATH)/gqaps
