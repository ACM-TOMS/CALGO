include ../make.gbinc

### GEMM-Based Level 3 BLAS Benchmark ####################################
#
#  The following libraries are specified, a user specified level 3 BLAS
#  library to be timed (LIB3B), a library with underlying BLAS routines
#  (LIB12B) where the underlying BLAS routine SGEMM may be specified
#  separately, and the library with the timing functions SECOND and
#  DSECND (SSEC).
#
LIB3B  = $(ULIB)
SGEMM  = $(UULIB)
LIB12B = $(UULIB)
SSEC   = $(UTMG)
#
#  LIB specifies the order in which the libraries are linked with the
#  benchmark programs. Notice that the built-in GEMM-based routines
#  will be linked the first SGEMM, level 1 and 2 BLAS routines found
#  as underlying routines. Change the order in which the libraries are
#  linked as desired.
#
LIB = $(SSEC) $(LIB3B) $(SGEMM) $(LIB12B)
#
### Compiler options #####################################################
#
#  The compiler options are specified as follows for the different
#  programs and libraries:
#
#     SBMFLG  : the GEMM-based performance benchmark programs
#     SPRFLG  : routines that print the output results
#     SGBFLG  : the built-in GEMM-based level 3 BLAS routines
#     SAXFLG  : GEMM-based specific auxiliary routines
#     AXOPT   : other auxiliary routines
#
SBMFLG = $(GBBOPT)
SPRFLG = $(AXOPT)
SGBFLG = $(GBOPT)
SAXFLG = $(GBOPT)
AXFLG  = $(AXOPT)
#
########################################################################

STIMS =	sgbtim.f

STIM =	sgbtim.o

SBMS =	sgbt01.f sgbt02.f

SBM =	sgbt01.o sgbt02.o

SPRS =  sgbtp1.f sgbtp2.f

SPR =   sgbtp1.o sgbtp2.o

SGBS =	sgb02.f sgb04.f sgb06.f sgb08.f sgb09.f

SGB =	sgb02.o sgb04.o sgb06.o sgb08.o sgb09.o

SAUXS =	sgb90.f sgb91.f

SAUX =	sgb90.o sgb91.o

AUXS =  lsame.f xerbla.f

AUX =   lsame.o xerbla.o

SPRMS = ssbpm.f

SPRM =  ssbpm.o

AUXS2 =	getwrd.f eoln.f

AUX2 =	getwrd.o eoln.o

OBJ1 = $(STIM) $(SBM) $(SGB) $(SAUX) $(AUX) $(AUX2) $(SPR)

OBJ2 = $(SPRM) $(AUX2)

########################################################################
all:	sgbtim ssbpm

sgbtim:		$(OBJ1)
	$(LOADER) $(LOADOPT) -o sgbtim $(OBJ1) $(LIB)

ssbpm:		$(OBJ2)
	$(LOADER) $(LOADOPT) -o ssbpm $(OBJ2)

$(STIM):	$(STIMS)
	$(FORTRAN) -c $(SBMFLG) $(STIMS)

$(SBM):		$(SBMS)
	$(FORTRAN) -c $(SBMFLG) $(SBMS)

$(SPR):		$(SPRS)
	$(FORTRAN) -c $(SPRFLG) $(SPRS)

$(SGB):		$(SGBS)
	$(FORTRAN) -c $(SGBFLG) $(SGBS)

$(SAUX):	$(SAUXS)
	$(FORTRAN) -c $(SAXFLG) $(SAUXS)

$(AUX):		$(AUXS)
	$(FORTRAN) -c $(AXFLG) $(AUXS)

$(SPRM):	$(SPRMS)
	$(FORTRAN) -c $(AXFLG) $(SPRMS)

$(AUX2):	$(AUXS2)
	$(FORTRAN) -c $(AXFLG) $(AUXS2)

clean:
	rm -f *.o sgbtim ssbpm
