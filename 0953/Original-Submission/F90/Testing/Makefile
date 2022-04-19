include ../make.inc

LIB  =  ../libpdhseqr.a

all: test_random test_bench

test_random: test_random.o
	$(LOADER) $(LDFLAGS) -o test_random.exe test_random.o $(LIB) $(LIBS)

test_bench: test_bench.o
	$(LOADER) $(LDFLAGS) -o test_bench.exe test_bench.o $(LIB) $(LIBS)

quicktest: test_random test_bench
	./test_random.exe < quick1.in > quick.out
	./test_bench.exe < quick2.in >> quick.out
	grep "out of" quick.out

clean:
	rm -f test_random.o test_random.exe test_bench.o test_bench.exe quick.out
