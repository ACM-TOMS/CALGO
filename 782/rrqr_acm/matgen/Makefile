#######################################################################
#
#  Makefile for generating Matrix Generator of LAPACK 2.0
#
#######################################################################

FORTRAN   = f77
OPTS      = -u -O -dalign

ARCH      = ar
ARCHFLAGS = cr
RANLIB    = ranlib

#
#######################################################################

# No more changes are required beyond this line.

#######################################################################
#
#######################################################################
 
TMGLIB   = ../tmglib.a



SCATGEN = slatm1.o slaran.o slarnd.o

SMATGEN = slatms.o slatme.o slatmr.o \
   slagge.o slagsy.o slarge.o slaror.o slarot.o slatm2.o slatm3.o

CMATGEN = clatms.o clatme.o clatmr.o \
   clagge.o claghe.o clagsy.o clarge.o claror.o clarot.o clatm1.o \
   clarnd.o clatm2.o clatm3.o

DZATGEN = dlatm1.o dlaran.o dlarnd.o

DMATGEN = dlatms.o dlatme.o dlatmr.o \
   dlagge.o dlagsy.o dlarge.o dlaror.o dlarot.o dlatm2.o dlatm3.o

ZMATGEN = zlatms.o zlatme.o zlatmr.o \
   zlagge.o zlaghe.o zlagsy.o zlarge.o zlaror.o zlarot.o zlatm1.o \
   zlarnd.o zlatm2.o zlatm3.o

all: single complex double complex16

single: $(SMATGEN) $(SCATGEN)
	$(ARCH) $(ARCHFLAGS) $(TMGLIB) $(SMATGEN) $(SCATGEN)
	$(RANLIB) $(TMGLIB)

complex: $(CMATGEN) $(SCATGEN)
	$(ARCH) $(ARCHFLAGS) $(TMGLIB) $(CMATGEN) $(SCATGEN)
	$(RANLIB) $(TMGLIB)

double: $(DMATGEN) $(DZATGEN)
	$(ARCH) $(ARCHFLAGS) $(TMGLIB) $(DMATGEN) $(DZATGEN)
	$(RANLIB) $(TMGLIB)

complex16: $(ZMATGEN) $(DZATGEN)
	$(ARCH) $(ARCHFLAGS) $(TMGLIB) $(ZMATGEN) $(DZATGEN)
	$(RANLIB) $(TMGLIB)

$(SCATGEN): $(FRC)
$(SMATGEN): $(FRC)
$(CMATGEN): $(FRC)
$(DZATGEN): $(FRC)
$(DMATGEN): $(FRC)
$(ZMATGEN): $(FRC)
 
FRC:
	@FRC=$(FRC)

clean: ; \
        rm -f *.o

.f.o: ; $(FORTRAN) -c $(OPTS) $<
