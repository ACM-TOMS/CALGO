include ../make.gbinc

### GEMM-Based Level 3 BLAS Benchmark ####################################
#
#  The following libraries are specified, a user specified level 3 BLAS
#  library to be timed (LIB3B), a library with underlying BLAS routines
#  (LIB12B) where the underlying BLAS routine CGEMM may be specified
#  separately, and the library with the timing functions SECOND and
#  DSECND (CSEC).
#
LIB3B  = $(ULIB)
CGEMM  = $(UULIB)
LIB12B = $(UULIB)
CSEC   = $(UTMG)
#
#  LIB specifies the order in which the libraries are linked with the
#  benchmark programs. Notice that the built-in GEMM-based routines
#  will be linked the first CGEMM, level 1 and 2 BLAS routines found
#  as underlying routines. Change the order in which the libraries are
#  linked as desired.
#
LIB = $(CSEC) $(LIB3B) $(CGEMM) $(LIB12B)
#
### Compiler options #####################################################
#
#  The compiler options are specified as follows for the different
#  programs and libraries:
#
#     CBMFLG  : the GEMM-based performance benchmark programs
#     CPRFLG  : routines that print the output results
#     CGBFLG  : the built-in GEMM-based level 3 BLAS routines
#     CAXFLG  : GEMM-based specific auxiliary routines
#     AXOPT   : other auxiliary routines
#
CBMFLG = $(GBBOPT)
CPRFLG = $(AXOPT)
CGBFLG = $(GBOPT)
CAXFLG = $(GBOPT)
AXFLG  = $(AXOPT)
#
########################################################################

CTIMS =	cgbtim.f

CTIM =	cgbtim.o

CBMS =	cgbt01.f cgbt02.f

CBM =	cgbt01.o cgbt02.o

CPRS =  cgbtp1.f cgbtp2.f

CPR =   cgbtp1.o cgbtp2.o

CGBS =	cgb02.f cgb03.f cgb04.f cgb05.f cgb06.f cgb07.f cgb08.f cgb09.f

CGB =	cgb02.o cgb03.o cgb04.o cgb05.o cgb06.o cgb07.o cgb08.o cgb09.o

CAUXS =	cgb90.f cgb91.f

CAUX =	cgb90.o cgb91.o

AUXS =  lsame.f xerbla.f

AUX =   lsame.o xerbla.o

CPRMS = csbpm.f

CPRM =  csbpm.o

AUXS2 =	getwrd.f eoln.f

AUX2 =	getwrd.o eoln.o

OBJ1 = $(CTIM) $(CBM) $(CGB) $(CAUX) $(AUX) $(AUX2) $(CPR)

OBJ2 = $(CPRM) $(AUX2)

########################################################################
all:	cgbtim csbpm

cgbtim:		$(OBJ1)
	$(LOADER) $(LOADOPT) -o cgbtim $(OBJ1) $(LIB)

csbpm:		$(OBJ2)
	$(LOADER) $(LOADOPT) -o csbpm $(OBJ2)

$(CTIM):	$(CTIMS)
	$(FORTRAN) -c $(CBMFLG) $(CTIMS)

$(CBM):		$(CBMS)
	$(FORTRAN) -c $(CBMFLG) $(CBMS)

$(CPR):		$(CPRS)
	$(FORTRAN) -c $(CPRFLG) $(CPRS)

$(CGB):		$(CGBS)
	$(FORTRAN) -c $(CGBFLG) $(CGBS)

$(CAUX):	$(CAUXS)
	$(FORTRAN) -c $(CAXFLG) $(CAUXS)

$(AUX):		$(AUXS)
	$(FORTRAN) -c $(AXFLG) $(AUXS)

$(CPRM):	$(CPRMS)
	$(FORTRAN) -c $(AXFLG) $(CPRMS)

$(AUX2):	$(AUXS2)
	$(FORTRAN) -c $(AXFLG) $(AUXS2)

clean:
	rm -f *.o cgbtim csbpm
