include ../make.inc

SRC = \
    dcroot.o dgesqb.o dgesqr.o dgesub.o dgesuv.o \
    dhabal.o dhabak.o dhaesu.o dhaex2.o dhaord.o \
    dhapvl.o dhapvb.o dhasub.o dhgpqr.o dlabmx.o \
    dladfb.o dlaesb.o dlaest.o dlanha.o dlaors.o \
    dlapqr.o dlapv2.o dlapvb.o dlasub.o dosgpv.o \
    dosgsb.o dosgsq.o dosgsu.o dosmsb.o dosmpv.o \
    dosmsq.o dshbal.o dshes.o  dshpvb.o dshpvl.o \
    dskmv.o  dskr2.o  dskr2k.o dskrkb.o dskupd.o \
    dtgpx2.o dtgpy2.o dtrqml.o ilahap.o

all: $(SRC)
	$(ARCH) $(ARCHFLAGS) $(HAPACKLIB) $(SRC)

clean:
	rm -f *.o

.f.o:
	$(FORTRAN) $(OPTS) -c $<
