include ../make.gbinc

### Compiler options #####################################################
#
#  The compiler options are specified as follows for the different
#  programs and libraries:
#
#     ZGBFLG  : the GEMM-based level 3 BLAS routines
#     ZAXFLG  : GEMM-based specific auxiliary routines
#     AXOPT   : other auxiliary routines
#
ZGBFLG = $(GBOPT)
ZAXFLG = $(GBOPT)
AXFLG  = $(AXOPT)
#
### GEMM-based Level 3 BLAS ##############################################

ZGBS =	 zsymm.f zhemm.f zsyrk.f zherk.f zsyr2k.f zher2k.f \
	 ztrmm.f ztrsm.f

ZGB =    zsymm.o zhemm.o zsyrk.o zherk.o zsyr2k.o zher2k.o \
	 ztrmm.o ztrsm.o

ZAUXS =  zbigp.f zcld.f

ZAUX =   zbigp.o zcld.o

AUXS =	 lsame.f xerbla.f

AUX =	 lsame.o xerbla.o

ZGPMS =  zsgpm.f getwrd.f eoln.f

ZGPM  =  zsgpm.o getwrd.o eoln.o

########################################################################
all:	$(GBL3B) zsgpm

$(GBL3B):	$(ZGB) $(ZAUX) $(AUX)
	$(ARCH) $(ARCHFLAGS) $(GBL3B) $(ZGB) $(ZAUX) $(AUX)
	$(RANLIB) $(GBL3B)

zsgpm:		$(ZGPM)
	$(LOADER) $(LOADOPT) -o zsgpm $(ZGPM)

$(ZGB):		$(ZGBS)
	$(FORTRAN) -c $(ZGBFLG) $(ZGBS)

$(ZAUX):        $(ZAUXS)
	$(FORTRAN) -c $(ZAXFLG) $(ZAUXS)

$(AUX):		$(AUXS)
	$(FORTRAN) -c $(AXFLG) $(AUXS)

$(ZGPM):	$(ZGPMS)
	$(FORTRAN) -c $(AXFLG) $(ZGPMS)

clean:
	rm -f *.o zsgpm
