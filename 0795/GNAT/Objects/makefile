# compiler environment :

CC=gcc
BIND=gnatbl
INCLULIBS = -I../System -I../Math_Lib/Numbers -I../Math_Lib/Matrices -I../Math_Lib/Polynomials -I../Math_Lib/Supports -I../Homotopy -I../Continuation -I../Root_Counts/Product -I../Root_Counts/Implift -I../Root_Counts/Stalift -I../Root_Counts/Dynlift -I../Root_Counts/Symmetry -I../Main
GNATFLAGS = -gnatv -O3 -gnatp
# GNATFLAGS = -O3 -gnatp
# GNATFLAGS = -gnatv
.SUFFIXES: .adb .ads .ali

.ads.o:
	$(CC) -c $(GNATFLAGS) $<
.adb.o:
	$(CC) -c $(GNATFLAGS) $<
.c.o:
	$(CC) -c $(CFLAGS) $<

# target routine :

phc:
	gnatmake -c $(INCLULIBS) $(GNATFLAGS) dispatch.adb
	gnatbl -o ../../bin/phc dispatch.ali

# example of using PHCPACK :

use_phc:
	gnatmake -c $(INCLULIBS) $(GNATFLAGS) use_phc.adb
	gnatbl -o /tmp/use_phc use_phc.ali

# cleaning up the object files :

clean: force
	/bin/rm -f  *.o *.ali

force:
