include ../../SLmake.inc

ZQRSRC = pzlahqr.o pzlaconsb.o pzlasmsub.o pzlacp3.o pzlawil.o pzrot.o \
         zlamsh.o zlaref.o zlanv2.o zlahqr.o
ZVCSRC = pztrevc.o 
#pzlatrs.o pzaxpy.o pzdotc.o pzdotu.o

all: complex16

complex16: $(ZQRSRC) $(ZVCSRC)
	$(ARCH) $(ARCHFLAGS) $(PXLAHQRLIB) $(ZQRSRC) $(ZVCSRC)
	$(RANLIB) $(PXLAHQRLIB)

$(ZQRSRC): $(FRC)
$(ZVCSRC): $(FRC)

FRC:
	@FRC=$(FRC)

clean:
	rm -f *.o

.c.o: 
	$(CC) $(CDEFS) $(CCFLAGS) -c $*.c

.f.o: 
	$(FORTRAN) $(OPTS) -c $*.f

