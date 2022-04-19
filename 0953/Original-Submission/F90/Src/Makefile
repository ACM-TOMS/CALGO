include ../make.inc

# pdlahqr.f is not up to date for testing purpose
OBJS =  pdrot.o pdgebal.o \
	pdhseqr.o pdlaqr0.o \
	pdlaqr1.o pdlaqr2.o pdlaqr4.o \
	pdlaqr3.o pdlaqr5.o dlaqr6.o \
	pdtrord.o pdtrsen.o \
	bdlaapp.o bdlaexc.o bdtrexc.o \
	pilaenvx.o piparmq.o \
	pdlamve.o dgemx.o dtrmx.o dlaref.o pdlahqr.o
LIB = ../libpdhseqr.a

all: lib
lib: $(OBJS)
	$(AR) $(ARFLAGS) $(LIB) $(OBJS)
	$(RANLIB) $(LIB)

clean:
	rm -f $(OBJS) $(LIB)
