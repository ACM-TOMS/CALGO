
#  --------------------------------------------------------------
#  makefile (SRC) of ADOL-C version 1.6 as of January 1,   1995
#  --------------------------------------------------------------
#  makefile for the library libad.a of ADOL-C
#  written for  GNU's  gcc and g++ compilers.
#

CFLAGS =  -O 
CC = g++
MCC = gcc 
LIBS = -lad -lm
libad: adouble.o avector.o taputil1.o taputil2.o taputil3.o hos_forward.o hov_forward.o fov_forward.o hos_reverse.o fos_reverse.o hov_reverse.o fov_reverse.o tayutil.o drivers.o driversc.o utils.o
	ranlib libad.a
	@echo 'Library created'
adouble.o: adouble.c adouble.h avector.h oplate.h taputil1.h 
	$(CC) -c $(CFLAGS) $(LIB) adouble.c
	ar rcv libad.a adouble.o
avector.o: avector.c adouble.h avector.h oplate.h taputil1.h
	$(CC) -c $(CFLAGS) $(LIB) avector.c
	ar rcv libad.a avector.o
taputil1.o: taputil1.c oplate.h taputil2.h tayutil.h usrparms.h
	$(MCC) -c $(CFLAGS) $(LIB) taputil1.c
	ar rcv libad.a taputil1.o
taputil2.o: taputil2.c dvlparms.h oplate.h taputil3.h usrparms.h
	$(MCC) -c $(CFLAGS) $(LIB) taputil2.c
	ar rcv libad.a taputil2.o
taputil3.o: taputil3.c dvlparms.h oplate.h tayutil.h taputil1.h taputil2.h usrparms.h
	$(MCC) -c $(CFLAGS) $(LIB) taputil3.c
	ar rcv libad.a taputil3.o
hos_forward.o: hos_forward.c dvlparms.h taputil1.h taputil2.h taputil3.h tayutil.h oplate.h usrparms.h
	$(MCC) -c $(CFLAGS) $(LIB) hos_forward.c
	ar rcv libad.a hos_forward.o
hov_forward.o: hov_forward.c dvlparms.h taputil1.h taputil2.h taputil3.h tayutil.h oplate.h usrparms.h
	$(MCC) -c $(CFLAGS) $(LIB) hov_forward.c
	ar rcv libad.a hov_forward.o
fov_forward.o: fov_forward.c dvlparms.h taputil1.h taputil2.h taputil3.h tayutil.h oplate.h usrparms.h
	$(MCC) -c $(CFLAGS) $(LIB) fov_forward.c
	ar rcv libad.a fov_forward.o
hos_reverse.o: hos_reverse.c dvlparms.h taputil1.h taputil2.h taputil3.h tayutil.h oplate.h usrparms.h
	$(MCC) -c $(CFLAGS) $(LIB) hos_reverse.c
	ar rcv libad.a hos_reverse.o
fos_reverse.o: fos_reverse.c dvlparms.h taputil1.h taputil2.h taputil3.h tayutil.h oplate.h usrparms.h
	$(MCC) -c $(CFLAGS) $(LIB) fos_reverse.c
	ar rcv libad.a fos_reverse.o
hov_reverse.o: hov_reverse.c dvlparms.h taputil1.h taputil2.h taputil3.h tayutil.h oplate.h usrparms.h
	$(MCC) -c $(CFLAGS) $(LIB) hov_reverse.c
	ar rcv libad.a hov_reverse.o
fov_reverse.o: fov_reverse.c dvlparms.h taputil1.h taputil2.h taputil3.h tayutil.h oplate.h usrparms.h
	$(MCC) -c $(CFLAGS) $(LIB) fov_reverse.c
	ar rcv libad.a fov_reverse.o
tayutil.o: tayutil.c dvlparms.h usrparms.h tayutil.h
	$(MCC) -c $(CFLAGS) $(LIB) tayutil.c
	ar rcv libad.a tayutil.o
utils.o: utils.c tayutil.h taputil3.h 
	$(CC) -c $(CFLAGS) $(LIB) utils.c
	ar rcv libad.a utils.o
driversc.o: driversc.c dvlparms.h adutilsc.h
	$(MCC) -c $(CFLAGS) $(LIB) driversc.c
	ar rcv libad.a driversc.o
drivers.o: drivers.c dvlparms.h adutils.h
	$(CC) -c $(CFLAGS) $(LIB) drivers.c
	ar rcv libad.a drivers.o
clean:
	/bin/rm *.o

