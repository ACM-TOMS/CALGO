include ../make.inc

LIB  =  ../libpdhseqr.a

all: tune1 tune2

test1: test.o
	$(LOADER) $(LDFLAGS) -o test.exe test.o $(LIB) $(LIBS)

tune1: test1 tune1.o
	$(LOADER) $(LDFLAGS) -o tune1.exe tune1.o $(LIB) $(LIBS)

test2:
	$(FC) -c -cpp -DSTNB=10 piparmq.f -o nibble10.o
	$(FC) -c -cpp -DSTNB=10 test.f -o test10.o
	$(LOADER) $(LDFLAGS) -o test10.exe test10.o nibble10.o $(LIB) $(LIBS)
	$(FC) -c -cpp -DSTNB=20 piparmq.f -o nibble20.o
	$(FC) -c -cpp -DSTNB=20 test.f -o test20.o
	$(LOADER) $(LDFLAGS) -o test20.exe test20.o nibble20.o $(LIB) $(LIBS)
	$(FC) -c -cpp -DSTNB=30 piparmq.f -o nibble30.o
	$(FC) -c -cpp -DSTNB=30 test.f -o test30.o
	$(LOADER) $(LDFLAGS) -o test30.exe test30.o nibble30.o $(LIB) $(LIBS)
	$(FC) -c -cpp -DSTNB=40 piparmq.f -o nibble40.o
	$(FC) -c -cpp -DSTNB=40 test.f -o test40.o
	$(LOADER) $(LDFLAGS) -o test40.exe test40.o nibble40.o $(LIB) $(LIBS)
	$(FC) -c -cpp -DSTNB=50 piparmq.f -o nibble50.o
	$(FC) -c -cpp -DSTNB=50 test.f -o test50.o
	$(LOADER) $(LDFLAGS) -o test50.exe test50.o nibble50.o $(LIB) $(LIBS)
	$(FC) -c -cpp -DSTNB=60 piparmq.f -o nibble60.o
	$(FC) -c -cpp -DSTNB=60 test.f -o test60.o
	$(LOADER) $(LDFLAGS) -o test60.exe test60.o nibble60.o $(LIB) $(LIBS)
	$(FC) -c -cpp -DSTNB=70 piparmq.f -o nibble70.o
	$(FC) -c -cpp -DSTNB=70 test.f -o test70.o
	$(LOADER) $(LDFLAGS) -o test70.exe test70.o nibble70.o $(LIB) $(LIBS)
	$(FC) -c -cpp -DSTNB=80 piparmq.f -o nibble80.o
	$(FC) -c -cpp -DSTNB=80 test.f -o test80.o
	$(LOADER) $(LDFLAGS) -o test80.exe test80.o nibble80.o $(LIB) $(LIBS)
	$(FC) -c -cpp -DSTNB=90 piparmq.f -o nibble90.o
	$(FC) -c -cpp -DSTNB=90 test.f -o test90.o
	$(LOADER) $(LDFLAGS) -o test90.exe test90.o nibble90.o $(LIB) $(LIBS)

tune2: test2 tune2.o
	$(LOADER) $(LDFLAGS) -o tune2.exe tune2.o nibble10.o $(LIB) $(LIBS)

clean:
	rm -f test*.o test*.exe tune?.o tune?.exe nibble*.o
