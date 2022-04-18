include ../make.gbinc

### Compiler options #####################################################
#
#  The compiler options are specified as follows for the different
#  programs and libraries:
#
#     CGBFLG  : the GEMM-based level 3 BLAS routines
#     CAXFLG  : GEMM-based specific auxiliary routines
#     AXOPT   : other auxiliary routines
#
CGBFLG = $(GBOPT)
CAXFLG = $(GBOPT)
AXFLG  = $(AXOPT)
#
### GEMM-based Level 3 BLAS ##############################################

CGBS =	 csymm.f chemm.f csyrk.f cherk.f csyr2k.f cher2k.f \
	 ctrmm.f ctrsm.f

CGB =    csymm.o chemm.o csyrk.o cherk.o csyr2k.o cher2k.o \
	 ctrmm.o ctrsm.o

CAUXS =  cbigp.f ccld.f

CAUX =   cbigp.o ccld.o

AUXS =	 lsame.f xerbla.f

AUX =	 lsame.o xerbla.o

CGPMS =  csgpm.f getwrd.f eoln.f

CGPM  =  csgpm.o getwrd.o eoln.o

########################################################################
all:	$(GBL3B) csgpm

$(GBL3B):	$(CGB) $(CAUX) $(AUX)
	$(ARCH) $(ARCHFLAGS) $(GBL3B) $(CGB) $(CAUX) $(AUX)
	$(RANLIB) $(GBL3B)

csgpm:		$(CGPM)
	$(LOADER) $(LOADOPT) -o csgpm $(CGPM)

$(CGB):		$(CGBS)
	$(FORTRAN) -c $(CGBFLG) $(CGBS)

$(CAUX):        $(CAUXS)
	$(FORTRAN) -c $(CAXFLG) $(CAUXS)

$(AUX):		$(AUXS)
	$(FORTRAN) -c $(AXFLG) $(AUXS)

$(CGPM):	$(CGPMS)
	$(FORTRAN) -c $(AXFLG) $(CGPMS)

clean:
	rm -f *.o csgpm
