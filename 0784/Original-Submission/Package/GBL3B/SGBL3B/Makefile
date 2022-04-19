include ../make.gbinc

### Compiler options #####################################################
#
#  The compiler options are specified as follows for the different
#  programs and libraries:
#
#     SGBFLG  : the GEMM-based level 3 BLAS routines
#     SAXFLG  : GEMM-based specific auxiliary routines
#     AXOPT   : other auxiliary routines
#
SGBFLG = $(GBOPT)
SAXFLG = $(GBOPT)
AXFLG  = $(AXOPT)
#
### GEMM-based Level 3 BLAS ##############################################

AUXS =	 lsame.f xerbla.f

AUX =	 lsame.o xerbla.o

SGBS =	 ssymm.f ssyrk.f ssyr2k.f strmm.f strsm.f

SGB =	 ssymm.o ssyrk.o ssyr2k.o strmm.o strsm.o

SAUXS =	 sbigp.f scld.f

SAUX =	 sbigp.o scld.o

SGPMS =  ssgpm.f getwrd.f eoln.f

SGPM  =  ssgpm.o getwrd.o eoln.o

########################################################################
all:	$(GBL3B) ssgpm

$(GBL3B):	$(SGB) $(SAUX) $(AUX)
	$(ARCH) $(ARCHFLAGS) $(GBL3B) $(SGB) $(SAUX) $(AUX)
	$(RANLIB) $(GBL3B)

ssgpm:		$(SGPM)	
	$(LOADER) $(LOADOPT) -o ssgpm $(SGPM)

$(SGB):		$(SGBS)
	$(FORTRAN) -c $(SGBFLG) $(SGBS)

$(SAUX):	$(SAUXS)
	$(FORTRAN) -c $(SAXFLG) $(SAUXS)

$(AUX):		$(AUXS)
	$(FORTRAN) -c $(AXFLG) $(AUXS)

$(SGPM):	$(SGPMS)
	$(FORTRAN) -c $(AXFLG) $(SGPMS)

clean:
	rm -f *.o ssgpm
