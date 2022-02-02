include ../make.inc

OBJS =  matpde.o matran.o pdlaprnt.o pdmatgen2.o pmatgeninc.o pqrbench.o \
	pqrhlps.o pqrrmmm.o
LIB = ../libpdhseqr.a

all: $(OBJS)
	$(AR) $(ARFLAGS) $(LIB) $(OBJS)
	$(RANLIB) $(LIB)

clean:
	rm -f $(OBJS)
