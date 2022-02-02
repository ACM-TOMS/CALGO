include ../make.gbinc

### GEMM-Based Level 3 BLAS Benchmark ####################################
#
#  The following libraries are specified, a user specified level 3 BLAS
#  library to be timed (LIB3B), a library with underlying BLAS routines
#  (LIB12B) where the underlying BLAS routine ZGEMM may be specified
#  separately, and the library with the timing functions SECOND and
#  DSECND (ZSEC).
#
LIB3B  = $(ULIB)
ZGEMM  = $(UULIB)
LIB12B = $(UULIB)
ZSEC   = $(UTMG)
#
#  LIB specifies the order in which the libraries are linked with the
#  benchmark programs. Notice that the built-in GEMM-based routines
#  will be linked the first ZGEMM, level 1 and 2 BLAS routines found
#  as underlying routines. Change the order in which the libraries are
#  linked as desired.
#
LIB = $(ZSEC) $(LIB3B) $(ZGEMM) $(LIB12B)
#
### Compiler options #####################################################
#
#  The compiler options are specified as follows for the different
#  programs and libraries:
#
#     ZBMFLG  : the GEMM-based performance benchmark programs
#     ZPRFLG  : routines that print the output results
#     ZGBFLG  : the built-in GEMM-based level 3 BLAS routines
#     ZAXFLG  : GEMM-based specific auxiliary routines
#     AXOPT   : other auxiliary routines
#
ZBMFLG = $(GBBOPT)
ZPRFLG = $(AXOPT)
ZGBFLG = $(GBOPT)
ZAXFLG = $(GBOPT)
AXFLG  = $(AXOPT)
#
########################################################################

ZTIMS =	zgbtim.f

ZTIM =	zgbtim.o

ZBMS =	zgbt01.f zgbt02.f

ZBM =	zgbt01.o zgbt02.o

ZPRS =  zgbtp1.f zgbtp2.f

ZPR =   zgbtp1.o zgbtp2.o

ZGBS =	zgb02.f zgb03.f zgb04.f zgb05.f zgb06.f zgb07.f zgb08.f zgb09.f

ZGB =	zgb02.o zgb03.o zgb04.o zgb05.o zgb06.o zgb07.o zgb08.o zgb09.o

ZAUXS =	zgb90.f zgb91.f

ZAUX =	zgb90.o zgb91.o

AUXS =  lsame.f xerbla.f

AUX =   lsame.o xerbla.o

ZPRMS = zsbpm.f

ZPRM =  zsbpm.o

AUXS2 =	getwrd.f eoln.f

AUX2 =	getwrd.o eoln.o

OBJ1 = $(ZTIM) $(ZBM) $(ZGB) $(ZAUX) $(AUX) $(AUX2) $(ZPR)

OBJ2 = $(ZPRM) $(AUX2)

########################################################################
all:	zgbtim zsbpm

zgbtim:		$(OBJ1)
	$(LOADER) $(LOADOPT) -o zgbtim $(OBJ1) $(LIB)

zsbpm:		$(OBJ2)
	$(LOADER) $(LOADOPT) -o zsbpm $(OBJ2)

$(ZTIM):	$(ZTIMS)
	$(FORTRAN) -c $(ZBMFLG) $(ZTIMS)

$(ZBM):		$(ZBMS)
	$(FORTRAN) -c $(ZBMFLG) $(ZBMS)

$(ZPR):		$(ZPRS)
	$(FORTRAN) -c $(ZPRFLG) $(ZPRS)

$(ZGB):		$(ZGBS)
	$(FORTRAN) -c $(ZGBFLG) $(ZGBS)

$(ZAUX):	$(ZAUXS)
	$(FORTRAN) -c $(ZAXFLG) $(ZAUXS)

$(AUX):		$(AUXS)
	$(FORTRAN) -c $(AXFLG) $(AUXS)

$(ZPRM):	$(ZPRMS)
	$(FORTRAN) -c $(AXFLG) $(ZPRMS)

$(AUX2):	$(AUXS2)
	$(FORTRAN) -c $(AXFLG) $(AUXS2)

clean:
	rm -f *.o zgbtim zsbpm

