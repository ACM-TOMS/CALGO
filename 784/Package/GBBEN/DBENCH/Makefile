include ../make.gbinc

### GEMM-Based Level 3 BLAS Benchmark ####################################
#
#  The following libraries are specified, a user specified level 3 BLAS
#  library to be timed (LIB3B), a library with underlying BLAS routines
#  (LIB12B) where the underlying BLAS routine DGEMM may be specified
#  separately, and the library with the timing functions SECOND and
#  DSECND (DSEC).
#
LIB3B  = $(ULIB)
DGEMM  = $(UULIB)
LIB12B = $(UULIB)
DSEC   = $(UTMG)
#
#  LIB specifies the order in which the libraries are linked with the
#  benchmark programs. Notice that the built-in GEMM-based routines
#  will be linked the first DGEMM, level 1 and 2 BLAS routines found
#  as underlying routines. Change the order in which the libraries are
#  linked as desired.
#
LIB = $(DSEC) $(LIB3B) $(DGEMM) $(LIB12B)
#
### Compiler options #####################################################
#
#  The compiler options are specified as follows for the different
#  programs and libraries:
#
#     DBMFLG  : the GEMM-based performance benchmark programs
#     DPRFLG  : routines that print the output results
#     DGBFLG  : the built-in GEMM-based level 3 BLAS routines
#     DAXFLG  : GEMM-based specific auxiliary routines
#     AXOPT   : other auxiliary routines
#
DBMFLG = $(GBBOPT)
DPRFLG = $(AXOPT)
DGBFLG = $(GBOPT)
DAXFLG = $(GBOPT)
AXFLG  = $(AXOPT)
#
########################################################################

DTIMS =	dgbtim.f

DTIM =	dgbtim.o

DBMS =	dgbt01.f dgbt02.f

DBM =	dgbt01.o dgbt02.o

DPRS =  dgbtp1.f dgbtp2.f

DPR =   dgbtp1.o dgbtp2.o

DGBS =	dgb02.f dgb04.f dgb06.f dgb08.f dgb09.f

DGB =	dgb02.o dgb04.o dgb06.o dgb08.o dgb09.o

DAUXS =	dgb90.f dgb91.f

DAUX =	dgb90.o dgb91.o

AUXS =  lsame.f xerbla.f

AUX =   lsame.o xerbla.o

DPRMS = dsbpm.f

DPRM =  dsbpm.o

AUXS2 =	getwrd.f eoln.f

AUX2 =	getwrd.o eoln.o

OBJ1 = $(DTIM) $(DBM) $(DGB) $(DAUX) $(AUX) $(AUX2) $(DPR)

OBJ2 = $(DPRM) $(AUX2)

########################################################################
all:	dgbtim dsbpm

dgbtim:		$(OBJ1)
	$(LOADER) $(LOADOPT) -o dgbtim $(OBJ1) $(LIB)

dsbpm:		$(OBJ2)
	$(LOADER) $(LOADOPT) -o dsbpm $(OBJ2)

$(DTIM):	$(DTIMS)
	$(FORTRAN) -c $(DBMFLG) $(DTIMS)

$(DBM):		$(DBMS)
	$(FORTRAN) -c $(DBMFLG) $(DBMS)

$(DPR):		$(DPRS)
	$(FORTRAN) -c $(DPRFLG) $(DPRS)

$(DGB):		$(DGBS)
	$(FORTRAN) -c $(DGBFLG) $(DGBS)

$(DAUX):	$(DAUXS)
	$(FORTRAN) -c $(DAXFLG) $(DAUXS)

$(AUX):		$(AUXS)
	$(FORTRAN) -c $(AXFLG) $(AUXS)

$(DPRM):	$(DPRMS)
	$(FORTRAN) -c $(AXFLG) $(DPRMS)

$(AUX2):	$(AUXS2)
	$(FORTRAN) -c $(AXFLG) $(AUXS2)

clean:
	rm -f *.o dgbtim dsbpm
