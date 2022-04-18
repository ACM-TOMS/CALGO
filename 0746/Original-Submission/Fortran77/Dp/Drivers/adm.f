C*********************************************************************
C
C          PROGRAM FOR TESTING AUTOMATIC DIFFERENTIATION ROUTINE
C
C                            P C O M P
C
C
C     The program reads a PCOMP input file and some data, to test
C     various options of the automatic differentiation algorithm.
C     Functions and derivatives are evaluated in the following way:
C
C     1. Symbolic interpretation and evaluation of input functions
C     2. Symbolic interpretation of input functions and evaluation
C        of gradients by automatic differentiation in forward mode
C     3. Symbolic interpretation of input functions and numerical
C        approximation of gradient by forward difference formula
C     4. Generation and execution of Fortran code for function
C        evaluation
C     5. Generation and execution of Fortran code for function
C        and gradient evaluation, in the latter case by automatic
C        differentiation in reverse mode
C     6. Generation and execution of Fortran code for function
C        evaluation and approximation of gradients by forward difference
C        formula
C
C     The numerical performance is evaluated by execution time in
C     seconds by mean values w.r.t. a user-provided number of repeated
C     test runs. For each type of gradient evaluation, work ratios are
C     computed, i.e. the relative time of one function and gradient
C     calculation divided by the time required for one function
C     calculation.
C
C     The test program ADM requires a timing routine, i.e. a subroutine
C     called EVLTIM(T) to fix the current system time. This subroutine
C     must be adapted if the compiler is changed.
C
C     To execute ADM, nonlinear functions must be defined in the PCOMP
C     language in a file called <example>.FUN. The syntax is described
C     in
C       Dobmann M., Liepelt M., Schittkowski K. (1995):
C       Algorithm 746: PCOMP: A FORTRAN code for automatic
C       differentiation, ACM Transactions on Mathematical Software,
C       Vol. 21, No. 3, 233-266
C     and some extensions of Version 5.3 in forthcoming remarks.
C     To evaluate functions and gradients, we need some data in a
C     file <example>.INP in the following form:
C     line 1: arbitrary name of the example (A50)
C     line 2: REPEAT, i.e. number of repeated function and gradient
C             evaluations to get a more reliable measure for execution
C             time (10X,I5)
C     line 3: N, i.e. number of variables (10X,I5)
C     line 4: M, i.e. number of functions (10X,I5)
C     line 5: I, i.e. index of function to be evaluated, 0<I<=M (10X,I5)
C     line 6: X(1), i.e. value of variable no. 1 (10X,G20.5)
C             ...
C     line 5+N: X(N), i.e. value of variable no. N (10X,G20.5)
C
C
C     Author:    K. Schittkowski
C                Department of Mathematics
C                University of Bayreuth
C                D - 98440 Bayreuth
C
C                +921 553278
C                klaus.schittkowski@uni-bayreuth.de
C                http://www.uni-bayreuth.de/math/dmv/~kschittkowski/
C                       home.htm
C
C     Copyright: Schittkowski, April 1997
C
C*********************************************************************
C
      INTEGER NMAX,MMAX,LRSYM,LISYM
      PARAMETER (NMAX=100,MMAX=50,LRSYM=20000,LISYM=4000)
C
      DOUBLE PRECISION X(NMAX),F(MMAX),FX(MMAX),DF(MMAX,NMAX),
     /  DFNUM(MMAX,NMAX),DFX(MMAX,NMAX),DFXNUM(MMAX,NMAX),FEPS(MMAX),
     /  FXEPS(MMAX),RSYM(LRSYM)
      INTEGER ISYM(LISYM),IDFX(NMAX)
      LOGICAL ACT(MMAX)
      CHARACTER NAME*50
      DOUBLE PRECISION ZERO,ONE,EPS,T1,T2,TSF,TSG,XEPS,XEPSI,TNU,ERRDF,
     /  TXF,TXG,TNX,ERRDFX,WRSYM,WRGEN,WRNGEN,WRNSYM
      INTEGER I,J,MAXFUN,N,M,NOFUNC,LROW,LARSYM,LAISYM,IERR,NV,NF,IEVL
C
C   OPEN FILES
C
      OPEN(1,FILE='ADM.INP',STATUS='UNKNOWN')
      OPEN(2,FILE='ADM.FUN',STATUS='UNKNOWN')
C      OPEN(3,FILE='ADM.SYM',STATUS='UNKNOWN')
      OPEN(10,FILE='ADM.RES',STATUS='UNKNOWN')
C
C   READ INPUT DATA
C
      READ(1,1000,ERR=850,END=9999) NAME
      READ(1,1100,ERR=850,END=9999) MAXFUN
      READ(1,1100,ERR=850,END=9999) N
      READ(1,1100,ERR=850,END=9999) M
      READ(1,1100,ERR=850,END=9999) NOFUNC
      DO 10 I=1,N
      IDFX(I)=I
        READ(1,1200,ERR=850,END=9999) X(I)
   10 CONTINUE
C
C   SET SOME CONSTANT DATA
C
      ZERO=0.D+0
      ONE=1.D+0
      EPS=1.D-7
      LROW=0
      DO 1 I=1,M
        ACT(I)=.TRUE.
    1 CONTINUE
C
C   PARSE FUNCTION INPUT FILE
C
C      CALL SYMPRP(3,RSYM,LRSYM,ISYM,LISYM,LARSYM,LAISYM,IERR,1,NV,NF)
      CALL SYMINP(2,0,RSYM,LRSYM,ISYM,LISYM,LARSYM,LAISYM,IERR,LROW,
     /            1,NV,NF,50,0)
      IF (IERR.GT.0) GOTO 900
      IF ((N.GT.NMAX).OR.(N.NE.NV)) GOTO 800
      IF ((M.GT.MMAX).OR.(M.NE.NF)) GOTO 810
C     IF (IERR.GT.0) GOTO 900
C
C   SYMBOLIC FUNCTION EVALUATION
C
      CALL EVLTIM(T1)
      DO 100 IEVL=1,N*MAXFUN
        CALL SYMFUN(X,N,F,M,ACT,RSYM,LARSYM,ISYM,LAISYM,IDFX,N,IERR)
  100 CONTINUE
      CALL EVLTIM(T2)
      IF (IERR.GT.0) GOTO 900
      TSF=(T2 - T1)/DBLE(N*MAXFUN)
C
C   SYMBOLIC GRADIENT EVALUATION
C
      CALL EVLTIM(T1)
      DO 200 IEVL=1,MAXFUN
        CALL SYMGRA(X,N,F,M,DF,MMAX,ACT,RSYM,LARSYM,ISYM,LAISYM,
     /            IDFX,N,IERR)
  200 CONTINUE
      CALL EVLTIM(T2)
      IF (IERR.GT.0) GOTO 900
      TSG=(T2 - T1)/DBLE(MAXFUN)
C
C   NUMERICAL GRADIENT EVALUATION OF SYMBOLIC FUNCTION INPUT
C
      CALL EVLTIM(T1)
      DO 300 IEVL=1,MAXFUN
      DO 301 I=1,N
      XEPS=EPS*DMAX1(ONE,DABS(X(I)))
      XEPSI=ONE/XEPS
      X(I)=X(I) + XEPS
      CALL SYMFUN(X,N,FEPS,M,ACT,RSYM,LARSYM,ISYM,LAISYM,
     /            IDFX,N,IERR)
      DO 302 J=1,M
      IF (.NOT.ACT(J)) GOTO 302
      DFNUM(J,I)=(FEPS(J) - F(J))*XEPSI
  302 CONTINUE
      X(I)=X(I) - XEPS
  301 CONTINUE
  300 CONTINUE
      CALL EVLTIM(T2)
      IF (IERR.GT.0) GOTO 900
      TNU=TSF + (T2 - T1)/DBLE(MAXFUN)
      ERRDF=ZERO
      DO 320 I=1,N
      DO 310 J=1,M
      ERRDF=ERRDF + DABS(DF(J,I)-DFNUM(J,I))
  310 CONTINUE
  320 CONTINUE
C
C   FUNCTION EVALUATION OF GENERATED FORTRAN CODE
C
      CALL EVLTIM(T1)
      DO 400 IEVL=1,5*N*MAXFUN
      CALL XFUN(X,N,FX,M,ACT,IERR)
  400 CONTINUE
      CALL EVLTIM(T2)
      IF (IERR.GT.0) GOTO 900
      TXF=(T2 - T1)/DBLE(5*N*MAXFUN)
C
C   GRADIENT EVALUATION OF GENERATED FORTRAN CODE
C
      CALL EVLTIM(T1)
      DO 500 IEVL=1,5*MAXFUN
      CALL XGRA(X,N,FX,M,DFX,MMAX,ACT,IERR)
  500 CONTINUE
      CALL EVLTIM(T2)
      IF (IERR.GT.0) GOTO 900
      TXG=(T2 - T1)/DBLE(5*MAXFUN)
C
C   NUMERICAL GRADIENT EVALUATION OF GENERATED FORTRAN CODE
C
      CALL EVLTIM(T1)
      DO 600 IEVL=1,5*MAXFUN
      DO 601 I=1,N
      XEPS=EPS*DMAX1(ONE,DABS(X(I)))
      XEPSI=ONE/XEPS
      X(I)=X(I) + XEPS
      CALL XFUN(X,N,FXEPS,M,ACT,IERR)
      DO 602 J=1,M
      IF (.NOT.ACT(J)) GOTO 602
      DFXNUM(J,I)=(FXEPS(J) - FX(J))*XEPSI
  602 CONTINUE
      X(I)=X(I) - XEPS
  601 CONTINUE
  600 CONTINUE
      CALL EVLTIM(T2)
      IF (IERR.GT.0) GOTO 900
      TNX=TXF + (T2 - T1)/DBLE(5*MAXFUN)
      ERRDFX=ZERO
      DO 620 I=1,N
      DO 610 J=1,M
C      WRITE(10,*) DF(J,I), DFNUM(J,I)
      ERRDFX=ERRDFX + DABS(DF(J,I)-DFXNUM(J,I))
  610 CONTINUE
  620 CONTINUE
C
C   OUTPUT ON RESULT
C
      WRITE(10,3000) NAME,MAXFUN,N,M,NOFUNC
      WRITE(10,3100) F(NOFUNC),FX(NOFUNC),ERRDF,ERRDFX
      IF (TSF.EQ.0) TSF=1.0
      IF (TXF.EQ.0) TXF=1.0
      WRSYM=TSG/TSF
      WRNSYM=TNU/TSF
      WRGEN=TXG/TXF
      WRNGEN=TNX/TXF
      TSF=100.0*TSF
      TXF=100.0*TXF
      TSG=100.0*TSG
      TNU=100.0*TNU
      TXG=100.0*TXG
      TNX=100.0*TNX
C
      WRITE(10,3200) TSF,TSG,TNU,WRSYM,WRNSYM,
     /               TXF,TXG,TNX,WRGEN,WRNGEN
      GOTO 9999
C
C   ERROR SITUATIONS
C
  800 WRITE(*,8000)
      GOTO 9999
  810 WRITE(*,8100)
      GOTO 9999
  850 WRITE(*,8500)
      GOTO 9999
  900 CALL SYMERR(LROW,IERR)
 9999 CONTINUE
C
C   CLOSE FILES
C
      CLOSE(1)
      CLOSE(2)
C      CLOSE(3)
      CLOSE(10)
C
C   FORMATTING STATEMENTS
C
 1000 FORMAT(A50)
 1100 FORMAT(10X,I5)
 1200 FORMAT(10X,G20.5)
 3000 FORMAT(/,' *** ',A50,/,
     /        /'     Function repeat:            ',I5,
     /        /'     Number of variables:        ',I5,
     /        /'     Number of functions:        ',I5,
     /        /'     Index of selected function: ',I5)
 3100 FORMAT(  '     Symbolic function value:    ',D22.14,
     /       /,'     Generated function value:   ',D22.14,
     /       /,'     Symbolic derivative error:  ',D10.3,
     /       /,'     Generated derivative error: ',D10.3,
     /       ////,'     Execution times:')
 3200 FORMAT(
     /    /,4X,'   code      function    gradient  num gradient ',
     /          'work-ratio  num work-ratio ',
     /    /,4X,75('-'),
     /    /,4X,'symbolic ',3(F12.6),2(F12.2),
     /    /,4X,'generated',3(F12.6),2(F12.2),/)
 8000 FORMAT(' *** ERROR: Wrong number of variables!')
 8100 FORMAT(' *** ERROR: Wrong number of functions!')
 8500 FORMAT(' *** ERROR: Input error in data file!')
C
C   END OF MAIN PROGRAMM
C
      STOP
      END
C
C-----------------------------------------------------------------------
C
      SUBROUTINE EVLTIM(T)
      DOUBLE PRECISION T
C     INTEGER*2 H,MIN,SEC,HSEC
      real second
C     CHARACTER*11 S
C
C  MS-FORTRAN:
C
C     CALL GETTIM(H,MIN,SEC,HSEC)
C
C  LAHEY:
C
C     CALL TIME(S)
C     READ(S,100) H
C 100 FORMAT(I2)
C     READ(S,200) MIN
C 200 FORMAT(3X,I2)
C     READ(S,300) SEC
C 300 FORMAT(6X,I2)
C     READ(S,400) HSEC
C 400 FORMAT(9X,I2)
C
C     T=0.01*DBLE(H*360000+MIN*6000+SEC*100+HSEC)
C     T=0.0D0
      T = dble(second())
      RETURN
      END
