include ../make.gbinc

### Compiler options #####################################################
#
#  The compiler options are specified as follows for the different
#  programs and libraries:
#
#     DGBFLG  : the GEMM-based level 3 BLAS routines
#     DAXFLG  : GEMM-based specific auxiliary routines
#     AXOPT   : other auxiliary routines
#
DGBFLG = $(GBOPT)
DAXFLG = $(GBOPT)
AXFLG  = $(AXOPT)
#
### GEMM-based Level 3 BLAS ##############################################

AUXS =	 lsame.f xerbla.f

AUX =	 lsame.o xerbla.o

DGBS =	 dsymm.f dsyrk.f dsyr2k.f dtrmm.f dtrsm.f

DGB =	 dsymm.o dsyrk.o dsyr2k.o dtrmm.o dtrsm.o

DAUXS =	 dbigp.f dcld.f

DAUX =	 dbigp.o dcld.o

DGPMS =  dsgpm.f getwrd.f eoln.f

DGPM  =  dsgpm.o getwrd.o eoln.o

########################################################################
all:	$(GBL3B) dsgpm

$(GBL3B):	$(DGB) $(DAUX) $(AUX)
	$(ARCH) $(ARCHFLAGS) $(GBL3B) $(DGB) $(DAUX) $(AUX)
	$(RANLIB) $(GBL3B)

dsgpm:		$(DGPM)
	$(LOADER) $(LOADOPT) -o dsgpm $(DGPM)

$(DGB):		$(DGBS)
	$(FORTRAN) -c $(DGBFLG) $(DGBS)

$(DAUX):	$(DAUXS)
	$(FORTRAN) -c $(DAXFLG) $(DAUXS)

$(AUX):		$(AUXS)
	$(FORTRAN) -c $(AXFLG) $(AUXS)

$(DGPM):	$(DGPMS)
	$(FORTRAN) -c $(AXFLG) $(DGPMS)

clean:
	rm -f *.o dsgpm
