include ../make.inc

LIB  =  ../libpdhseqr.a

all: example1 example2

example1: example1.o
	$(LOADER) $(LDFLAGS) -o example1.exe example1.o $(LIB) $(LIBS)

example2: example2.o
	$(LOADER) $(LDFLAGS) -o example2.exe example2.o $(LIB) $(LIBS)

clean:
	rm -f example1.o example1.exe example2.o example2.exe
