* --------------------------------------------------------------------
*
* SEGMENT 2:
* 2. Driver and Routines for Testing TNPACK-Minimization
*    for the Trigonometric Function.
*
*    To test the re-ordering option in TNPACK, a preconditioner
*    M = ( m_{ij} ), a matrix of order n, is defined by
*
*                m_{i,j} = H_{ii}  if  i = j
*                          0.1  if  i = 1, j = n-1
*                         -0.1  if  i = 1, j = n
*                          0.1  if  i = n-1, j = 1
*                         -0.1  if  i = n, j = 1
*                          0.0  otherwise.
*
*     where H_{ii} is the diagonal element of Hessian matrix.
*
************************************************************************

      PROGRAM MINTRG

* -----------------------------------------------------------------
* TNPACK-MinimizatiON dRIVer for the Trigonometric Function
* -----------------------------------------------------------------


C     .. Parameters ..
      INTEGER MAXM,MAXLH
      PARAMETER (MAXM=5000,MAXLH=MAXM* (MAXM-1)/2)
      INTEGER N,NZ,MP,LW,LIW
      PARAMETER (N=1000,NZ=3*N,MP=6,LW=9*N+5*NZ+7*N,LIW=7*N+5*NZ+1)
C     ..
C     .. Scalars in Common ..
      INTEGER NPROB
C     ..
C     .. Arrays in Common ..
      DOUBLE PRECISION HESD(MAXM),HESL(MAXLH)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION F,RI,XZERO
      REAL TIMEED,TIMEST,TIMEUR
      INTEGER I,INFORM
C     ..
C     .. Local Arrays ..
      DOUBLE PRECISION G(N),PLIST(20),W(LW),X(N)
      REAL TARRAY(2)
      INTEGER IW(LIW),OPLIST(20)
C     ..
C     .. External Subroutines ..
      EXTERNAL SETLIS,TNMIN,TRGFUN,TRGHDP,TRGPAT
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC COS,DBLE
C     ..
C     .. Common blocks ..
      COMMON NPROB
      COMMON /HESIAN/HESL,HESD
C     ..

C     .. External Functions ..
      REAL ETIME
      EXTERNAL ETIME
C     ..
      NPROB = 1

      IF (N.GT.MAXM) THEN
          WRITE (MP,FMT=9000) MAXM
          GO TO 20

      END IF

      XZERO = 1.D0/DBLE(N)
      DO 10 I = 1,N
          RI = I
          X(I) = XZERO + 0.2D0*COS(RI)
   10 CONTINUE

      CALL SETLIS(N,OPLIST,PLIST,INFORM)

* ---------------------------------------------
* Modify OPLIST elements as desired:
*    OPLIST(1)  - Preconditioning option
*    OPLIST(4)  - Max. Newton itns.
*    OPLIST(5)  - Max. PCG itns.
*    OPLIST(6)  - Max. F&G evals.
*    OPLIST(7)  - 1: reoder; 0: do not reorder
*    OPLIST(10) - Hd mult. option
*  Four new elements:
*    OPLIST(15) - 1: Test 1A'   2: Test 2A
*                 Default value is 1
*    OPLIST(16) - Option selector for line search stopping rule
*                1: Criterion 1
*                2: Criterion 2 (more lenient stopping rule)
*                Default value is 1
*    OPLIST(17) - Option selector for modified Cholesky factorization
*                1: the MC by Gill and Murray
*                2: the UMC by Xie and Schlick
*                Default value is 1
*    PLIST(8)   - UMC parameter tau
* ---------------------------------------------

      OPLIST(1) = 1
      OPLIST(5) = 40
      OPLIST(6) = 10*N
      OPLIST(7) = 0
      OPLIST(10) = 0
      OPLIST(12) = 1


c--- xie
      OPLIST(15) = 2
      OPLIST(16) = 2
      OPLIST(17) = 2
      PLIST(8) = 0.5D0
c-------------------------


      TIMEST = ETIME(TARRAY)
      TIMEUR = TARRAY(1)

      CALL TNMIN(N,X,F,G,OPLIST,PLIST,INFORM,NZ,W,LW,IW,LIW,TRGFUN,
     +           TRGPAT,TRGHDP)


      TIMEED = ETIME(TARRAY) - TIMEST
      TIMEUR = TARRAY(1) - TIMEUR

      WRITE (*,FMT=*) ' The CPU time in seconds ',TIMEUR

   20 CONTINUE


      STOP

 9000 FORMAT (/,5X,' N > MAXM,     MAXM  = ',I6,/,6X,
     +       'reset PARAMETER statement for MAXM',/)
      END
************************************************************************
      SUBROUTINE TRGFUN(N,X,F,G,A,IA,JA,NOUT)

* --------------------------------------------------------
* TRIGONOMETRIC FUNCTION OF DIMENSION N - ASSEMBLE F,G,H,M
* --------------------------------------------------------

C     .. Parameters ..
      INTEGER MAXM,MAXLH
      PARAMETER (MAXM=5000,MAXLH=MAXM* (MAXM-1)/2)
C     ..
C     .. Scalar Arguments ..
      DOUBLE PRECISION F
      INTEGER N,NOUT
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION A(*),G(N),X(N)
      INTEGER IA(*),JA(*)
C     ..
C     .. Arrays in Common ..
      DOUBLE PRECISION HESD(MAXM),HESL(MAXLH)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION SUMCOS,SUMTRM,TWO,ZERO
      INTEGER I,II,J,JJ,K
C     ..
C     .. Local Arrays ..
      DOUBLE PRECISION COSINE(MAXM),SINE(MAXM),TERM(MAXM)
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC COS,DBLE,SIN
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

      ZERO = 0.D0
      TWO = 2.D0

      SUMCOS = ZERO
      DO 10 I = 1,N
          SINE(I) = SIN(X(I))
          COSINE(I) = COS(X(I))
          SUMCOS = SUMCOS + COSINE(I)
   10 CONTINUE
      F = ZERO
      DO 20 I = 1,N
          TERM(I) = DBLE(N+I) - SUMCOS - DBLE(I)*COSINE(I) - SINE(I)
          F = F + TERM(I)**2
   20 CONTINUE
      IF (NOUT.EQ.0) RETURN

      SUMTRM = ZERO
      DO 30 I = 1,N
          SUMTRM = SUMTRM + TERM(I)
   30 CONTINUE
      DO 40 I = 1,N
          G(I) = SINE(I)*SUMTRM + (DBLE(I)*SINE(I)-COSINE(I))*TERM(I)
          G(I) = TWO*G(I)
   40 CONTINUE
      IF (NOUT.EQ.1) RETURN

      DO 50 I = 1,N
          HESD(I) = DBLE(I* (I+2)+N)*SINE(I)**2 +
     +              COSINE(I)* (SUMTRM+COSINE(I)-
     +              (DBLE(2*I+2)*SINE(I))) +
     +              TERM(I)* (DBLE(I)*COSINE(I)+SINE(I))
          HESD(I) = TWO*HESD(I)
   50 CONTINUE
      DO 70 I = 2,N
          DO 60 J = 1,I - 1
              K = IX(I,J)
              HESL(K) = DBLE(N+I+J)*SINE(I)*SINE(J) -
     +                  SINE(I)*COSINE(J) - SINE(J)*COSINE(I)
              HESL(K) = TWO*HESL(K)
   60     CONTINUE
   70 CONTINUE
      IF (NOUT.EQ.2) RETURN

      DO I = 4,N + 2
          A(I) = HESD(I-2)
      END DO
      A(1) = HESD(1)
      A(2) = 0.1D0
      A(3) = -0.1D0

      RETURN

      END
************************************************************************
      SUBROUTINE TRGHDP(N,D,HD,XC,GC)

* ------------------------------------------------------
* TRIGONOMETRIC FUNCTION OF DIMENSION N - COMPUTE HD
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
* TRIGONOMETRIC FUNCTION OF DIMENSION N - DETERMINE PATTERN OF M
* --------------------------------------------------------------


C     .. Scalar Arguments ..
      INTEGER N
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION A(*),X(*)
      INTEGER IA(N+1),JA(*)
C     ..
C     .. Local Scalars ..
      INTEGER I
C     ..
      DO 10 I = 1,N
          IA(I) = I
   10 CONTINUE

      DO I = 4,N + 2
          JA(I) = I - 2
      END DO
      JA(1) = 1
      JA(2) = N - 1
      JA(3) = N

      DO I = 2,N
          IA(I) = IA(I) + 2
      END DO

      IA(N+1) = IA(N) + 1

      RETURN

      END
