* ----------------------------------------------------------
*
* SEGMENT 3:
* 3. Driver and Routines for Testing TNPACK-Minimization
*    for the Functions defined in Alg. 566
*
************************************************************************
      PROGRAM MINFUC566

* -----------------------------------------------------------------
* TNPACK-Minimization Driver for the Functions in Alg.566
* -----------------------------------------------------------------



C     .. Parameters ..
      INTEGER MAXM,MAXLH
      PARAMETER (MAXM=5000,MAXLH=MAXM* (MAXM-1)/2)
      INTEGER NN,NZ,MP,LW,LIW
      PARAMETER (NN=1000,NZ=3*NN,MP=6,LW=10*NN+5*NZ,LIW=7*NN+5*NZ+1)
C     ..
C     .. Scalars in Common ..
      INTEGER NPROB
C     ..
C     .. Arrays in Common ..
      DOUBLE PRECISION HESD(MAXM),HESL(MAXLH)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION F,FACTOR
      INTEGER INFORM,N
C     ..
C     .. Local Arrays ..
      DOUBLE PRECISION G(NN),PLIST(20),W(LW),X(NN)
      INTEGER IW(LIW),OPLIST(20)
C     ..
C     .. External Subroutines ..
      EXTERNAL INITPT,PROBSET,SETLIS,TNMIN,TRGFUN,TRGHDP,TRGPAT
C     ..
C     .. Common blocks ..
      COMMON NPROB
      COMMON /HESIAN/HESL,HESD
C     ..
      NPROB = 0
   10 CONTINUE
      NPROB = NPROB + 1

      WRITE (MP,FMT=9010)
      WRITE (MP,FMT=9000) NPROB
      WRITE (MP,FMT=9010)

*
*--- Enter the value of N
*
      CALL PROBSET(NPROB,N)

      IF (N.GT.MAXM) THEN
          WRITE (MP,FMT=9020) MAXM
          GO TO 20

      END IF
*
*--- Initial guess
*
      FACTOR = 1.0
      CALL INITPT(N,X,NPROB,FACTOR)

      CALL SETLIS(N,OPLIST,PLIST,INFORM)

      OPLIST(1) = 1
      OPLIST(4) = 10000
      OPLIST(5) = 100
      OPLIST(6) = 4000
      OPLIST(10) = 0
      OPLIST(13) = 0
      PLIST(2) = 1.D-8
      PLIST(3) = 1.0D0

*-- xie
      OPLIST(15) = 2
      OPLIST(16) = 2
      OPLIST(17) = 2
      PLIST(8) = 10.D0
*------------------------------

      CALL TNMIN(N,X,F,G,OPLIST,PLIST,INFORM,NZ,W,LW,IW,LIW,TRGFUN,
     +           TRGPAT,TRGHDP)

   20 CONTINUE

      IF (NPROB.LT.18) GO TO 10


      STOP

 9000 FORMAT (/,20X,'Problem ',I2)
 9010 FORMAT (/,14X,'=======================')
 9020 FORMAT (/,5X,' N > MAXM,     MAXM  = ',I6,/,6X,
     +       'reset PARAMETER statement for MAXM',/)
      END
************************************************************************
      SUBROUTINE TRGFUN(N,X,F,G,A,IA,JA,NOUT)

* --------------------------------------------------------
*  ASSEMBLE F,G,H,M
* --------------------------------------------------------

C     .. Parameters ..
      INTEGER MAXM,MAXLH
      PARAMETER (MAXM=5000,MAXLH=MAXM* (MAXM-1)/2)
      INTEGER LBAND
      PARAMETER (LBAND=1)
C     ..
C     .. Scalar Arguments ..
      DOUBLE PRECISION F
      INTEGER N,NOUT
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION A(*),G(N),X(N)
      INTEGER IA(*),JA(*)
C     ..
C     .. Scalars in Common ..
      INTEGER NPROB
C     ..
C     .. Arrays in Common ..
      DOUBLE PRECISION HESD(MAXM),HESL(MAXLH)
C     ..
C     .. Local Scalars ..
      INTEGER I,II,J,JJ,K,M
C     ..
C     .. External Subroutines ..
      EXTERNAL GRDFCN,HESFCN,OBJFCN
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC MIN
C     ..
C     .. Common blocks ..
      COMMON NPROB
      COMMON /HESIAN/HESL,HESD
C     ..
C     .. Statement Functions ..
      INTEGER IX
C     ..
C     .. Save statement ..
      SAVE /HESIAN/
C     ..
C     .. Statement Function definitions ..

      IX(II,JJ) = ((II-1)* (II-2)/2) + JJ
C     ..

*---  THIS SUBROUTINE DEFINES THE OBJECTIVE FUNCTIONS OF EIGHTEEN
*     NONLINEAR UNCONSTRAINED MINIMIZATION PROBLEMS.
*     THE SUBROUTINE STATEMENT IS
*
*       SUBROUTINE OBJFCN(N,X,F,NPROB)
*
*     WHERE
*
*       N IS A POSITIVE INTEGER INPUT VARIABLE.
*
*       X IS AN INPUT ARRAY OF LENGTH N.
*
*       F IS AN OUTPUT VARIABLE WHICH CONTAINS THE VALUE OF
*         THE NPROB OBJECTIVE FUNCTION EVALUATED AT X.
*
*       NPROB IS A POSITIVE INTEGER INPUT VARIABLE WHICH DEFINES THE
*         NUMBER OF THE PROBLEM. NPROB MUST NOT EXCEED 18.

      CALL OBJFCN(N,X,F,NPROB)

      IF (NOUT.EQ.0) RETURN

*---  THIS SUBROUTINE DEFINES THE GRADIENT VECTORS OF EIGHTEEN
*     NONLINEAR UNCONSTRAINED MINIMIZATION PROBLEMS.
*
*       N IS A POSITIVE INTEGER INPUT VARIABLE.
*
*       X IS AN INPUT ARRAY OF LENGTH N.
*
*       G IS AN OUTPUT ARRAY OF LENGTH N WHICH CONTAINS THE COMPONENTS
*         OF THE GRADIENT VECTOR OF THE NPROB OBJECTIVE FUNCTION
*         EVALUATED AT X.
*
*       NPROB IS A POSITIVE INTEGER INPUT VARIABLE WHICH DEFINES THE
*         NUMBER OF THE PROBLEM. NPROB MUST NOT EXCEED 18.
*-----------------------------------------------------------------------

      CALL GRDFCN(N,X,G,NPROB)

      IF (NOUT.EQ.1) RETURN

*--- THIS SUBROUTINE DEFINES THE HESSIAN MATRICES OF 18
*     NONLINEAR UNCONSTRAINED MINIMIZATION PROBLEMS.
*
*       N IS A POSITIVE INTEGER INPUT VARIABLE.
*
*       X IS AN INPUT ARRAY OF LENGTH N.
*
*       HESD IS AN OUTPUT ARRAY OF LENGTH N WHICH CONTAINS THE
*         DIAGONAL COMPONENTS OF THE HESSIAN MATRIX OF THE NPROB
*         OBJECTIVE FUNCTION EVALUATED AT X.
*
*       HESL IS AN OUTPUT ARRAY OF LENGTH N*(N-1)/2 WHICH CONTAINS
*         THE LOWER TRIANGULAR PART OF THE HESSIAN MATRIX OF THE
*         NPROB OBJECTIVE FUNCTION EVALUATED AT X.
*
*       NPROB IS A POSITIVE INTEGER INPUT VARIABLE WHICH DEFINES THE
*         NUMBER OF THE PROBLEM.  NPROB MUST NOT EXCEED 18.
*--------------------------------------------------------------------

      CALL HESFCN(N,X,HESD,HESL,NPROB)

      IF (NOUT.EQ.2) RETURN

      K = 1
      DO 20 I = 1,N
          A(K) = HESD(I)
          M = MIN(LBAND/2,N-I)
          IF (M.GT.0) THEN
              DO 10 J = 1,M
                  A(K+J) = HESL(IX(I+J,I))
   10         CONTINUE
          END IF

          K = K + M + 1
   20 CONTINUE

      RETURN

      END
************************************************************************
      SUBROUTINE TRGHDP(N,D,HD,XC,GC)

* ------------------------------------------------------
*  COMPUTE HD
* ------------------------------------------------------

C     .. Parameters ..
      INTEGER MAXM,MAXLH
      PARAMETER (MAXM=5000,MAXLH=MAXM* (MAXM-1)/2)
C     ..
C     .. Scalar Arguments ..
      INTEGER N
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION D(N),GC(N),HD(N),XC(N)
C     ..
C     .. Arrays in Common ..
      DOUBLE PRECISION HESD(MAXM),HESL(MAXLH)
C     ..
C     .. Local Scalars ..
      INTEGER I,II,J,JJ
C     ..
C     .. Common blocks ..
      COMMON /HESIAN/HESL,HESD
C     ..
C     .. Statement Functions ..
      INTEGER IX
C     ..
C     .. Save statement ..
      SAVE /HESIAN/
C     ..
C     .. Statement Function definitions ..

      IX(II,JJ) = ((II-1)* (II-2)/2) + JJ
C     ..

      DO 10 I = 1,N
          HD(I) = HESD(I)*D(I)
   10 CONTINUE
      DO 40 I = 1,N
          DO 20 J = 1,I - 1
              HD(I) = HD(I) + HESL(IX(I,J))*D(J)
   20     CONTINUE
          DO 30 J = I + 1,N
              HD(I) = HD(I) + HESL(IX(J,I))*D(J)
   30     CONTINUE
   40 CONTINUE

      RETURN

      END
************************************************************************
      SUBROUTINE TRGPAT(N,X,A,IA,JA)

* --------------------------------------------------------------
* DETERMINE PATTERN OF M
* --------------------------------------------------------------


C     .. Parameters ..
      INTEGER LBAND
      PARAMETER (LBAND=1)
C     ..
C     .. Scalar Arguments ..
      INTEGER N
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION A(*),X(*)
      INTEGER IA(N+1),JA(*)
C     ..
C     .. Local Scalars ..
      INTEGER I,J,K,M
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC MIN
C     ..
      K = 1
      DO 20 I = 1,N
          IA(I) = K
          M = MIN(LBAND/2,N-I)
          DO 10 J = 0,M
              JA(K+J) = I + J
   10     CONTINUE
          K = K + M + 1
   20 CONTINUE
      IA(N+1) = K

      RETURN

      END

      SUBROUTINE PROBSET(NPROB,N)
C     .. Scalar Arguments ..
      INTEGER N,NPROB
C     ..
C     .. Local Arrays ..
      INTEGER NVARS(18)
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC MAX,MIN,MOD
C     ..
C     .. Data statements ..
      DATA NVARS/3,6,3,2,3,0,0,0,0,2,4,3,0,0,0,2,4,0/
C     ..

      IF (NPROB.EQ.0) STOP
      IF (NPROB.LT.1 .OR. NPROB.GT.18) THEN
          WRITE (*,FMT=9000)
          STOP

      END IF

      IF (NVARS(NPROB).NE.0) N = NVARS(NPROB)
      IF (NPROB.EQ.7) N = MAX(2,MIN(N,31))
      IF (NPROB.EQ.14) N = MAX(2,N-MOD(N,2))
      IF (NPROB.EQ.15) N = MAX(4,N-MOD(N,4))
      IF (NPROB.EQ.18) N = MAX(2,MIN(N,50))

      RETURN

 9000 FORMAT (/,4X,'ERROR IN INPUT FILE',/)
      END
