C ----------------------------------------------------------
C SEGMENT 1:
C 1. Driver and Routines for Testing TNPACK-Minimization
C    for the Extended Rosenbrock Function
C
C*********************************************************************
      PROGRAM MINROS

C -----------------------------------------------------------------
C TNPACK-Minimization Driver for the Extended Rosenbrock Function
C of dimension N (N even)
C -----------------------------------------------------------------


C     .. Parameters ..
      INTEGER MAXN
      PARAMETER (MAXN=1000)
      INTEGER N,NZ,MP,LW,LIW
      PARAMETER (N=1000,NZ=N,MP=6,LW=10*N+5*NZ,LIW=7*N+5*NZ+1)
C     ..
C     .. Scalars in Common ..
      INTEGER NPROB
C     ..
C     .. Arrays in Common ..
      DOUBLE PRECISION HESD(MAXN),HESOD(MAXN)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION F,RJ
      REAL TIMEED,TIMEST,TIMEUR
      INTEGER INFORM,J
C     ..
C     .. Local Arrays ..
      DOUBLE PRECISION G(N),PLIST(20),W(LW),X(N)
      REAL TARRAY(2)
      INTEGER IW(LIW),OPLIST(20)
C     ..
C     .. External Subroutines ..
      EXTERNAL ROSFUN,ROSHDP,ROSPAT,SETLIS,TNMIN
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC COS
C     ..
C     .. Common blocks ..
C     .. CPU time
      COMMON NPROB
      COMMON /HESIAN/HESD,HESOD
C     ..
C     .. External Functions ..
      REAL ETIME
      EXTERNAL ETIME
C     ..
      NPROB = 1

      IF (N.GT.MAXN) THEN
          WRITE (MP,FMT=9000) MAXN
          GO TO 20

      END IF
c
C--- initial guess
c
      DO 10 J = 1,N - 1,2
          RJ = J
          X(J) = -1.2D0 - COS(RJ)
          X(J+1) = 1.D0 + COS(RJ)
   10 CONTINUE

c
c--- Sets sample parameters for minimization
c
      CALL SETLIS(N,OPLIST,PLIST,INFORM)

C ---------------------------------------------
C Modify OPLIST and PLIST elements as desired:
C    OPLIST(1)  - Preconditioning option
C    OPLIST(5)  - Max. PCG itns.
C    OPLIST(10) - Hd mult. option
C    PLIST(2)   - Convergence parameter EPSG
C    PLIST(3)   - Residual truncation parameter
c Four new elements:
c    OPLIST(15) - Option selector for termination rule of PCG
C                 1: the modified negative curvature test;
C                 2: the descent direction test
C                 Default value is 1
C    OPLIST(16) - Option selector for line search stopping rule
C                1: Criterion 1
C                2: Criterion 2 (more lenient stopping rule)
C                Default value is 1
C    OPLIST(17) -   option selector for modified Cholesky factorization
C                1: the MC by Gill and Murray
C                2: the UMC by Xie and Schlick
C                Default value is 1
C    PLIST(8)   - Parameter tau of UMC. Default value is 10.d0
C ---------------------------------------------

      OPLIST(1) = 1
      OPLIST(5) = 100
      OPLIST(10) = 0
      OPLIST(12) = 1
      PLIST(2) = 1.0D-08
      PLIST(3) = 0.1

c--- xie
      OPLIST(15) = 1
      OPLIST(16) = 2
      OPLIST(17) = 2
      PLIST(8) = 10.0D0
c-------------------------

c
c--- Truncated Newton Interface Routine
c

      TIMEST = ETIME(TARRAY)
      TIMEUR = TARRAY(1)

      CALL TNMIN(N,X,F,G,OPLIST,PLIST,INFORM,NZ,W,LW,IW,LIW,ROSFUN,
     +           ROSPAT,ROSHDP)

      TIMEED = ETIME(TARRAY) - TIMEST
      TIMEUR = TARRAY(1) - TIMEUR

      WRITE (*,FMT=*) ' The CPU time in seconds ',TIMEUR

   20 CONTINUE

      STOP

 9000 FORMAT (/,5X,' N > MAXN,     MAXN  = ',I6,/,6X,
     +       'reset PARAMETER statement for MAXN',/)
      END

C*********************************************************************
      SUBROUTINE ROSFUN(N,X,F,G,A,IA,JA,NOUT)

C -------------------------------------------------------
C ROSENBROCK'S FUNCTION OF DIMENSION N - ASSEMBLE F,G,H,M
C -------------------------------------------------------

c
c--- F: ROSENBROCK'S FUNCTION
c
C     .. Parameters ..
      INTEGER MAXN
      PARAMETER (MAXN=1000)
C     ..
C     .. Scalar Arguments ..
      DOUBLE PRECISION F
      INTEGER N,NOUT
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION A(*),G(N),X(N)
      INTEGER IA(N+1),JA(*)
C     ..
C     .. Arrays in Common ..
      DOUBLE PRECISION HESD(MAXN),HESOD(MAXN)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION T1,T2
      INTEGER J
C     ..
C     .. Common blocks ..
      COMMON /HESIAN/HESD,HESOD
C     ..
      F = 0.D0
      DO 10 J = 1,N - 1,2
          T1 = 1.D0 - X(J)
          T2 = 10.D0* (X(J+1)-X(J)*X(J))
          F = F + T1*T1 + T2*T2
   10 CONTINUE

      IF (NOUT.EQ.0) RETURN
c
c--- Computing the componets G(j) of the gradient of F
c
      DO 20 J = 1,N - 1,2
          T1 = 1.D0 - X(J)
          G(J+1) = 200.D0* (X(J+1)-X(J)*X(J))
          G(J) = -2.D0* (X(J)*G(J+1)+T1)
   20 CONTINUE

      IF (NOUT.EQ.1) RETURN

C ------------------------------------------------------
c      Computing Hessian matrix
C H is stored in 2 arrays: HESD for diagonals, HESOD for
C off-diagonals, stored by rows for the upper triangle.
C HESOD is zero every even entry since H has the 2x2
C block-diagonal pattern:
C
C                    * *
C                    * *
C                        * *
C                        * *
C                                Etc.
C
C M is set to the diagonal of H.
c
c  H_jj = HESD( J )   H_{j,j+1} = HESOD(J)
c                     H_{j+1,j+1} = HESD(J+1)
C ------------------------------------------------------

      DO 30 J = 1,N - 1,2
          HESD(J+1) = 200.D0
          HESD(J) = 2.D0 - 400.D0*X(J+1) + 1200.D0*X(J)*X(J)
          HESOD(J) = -400.D0*X(J)
   30 CONTINUE

      IF (NOUT.EQ.2) RETURN
c
c--- Computing the preconditor M: the diagonal of H
c
      DO 40 J = 1,N - 1,2
          A(J) = -2.D0* (200.D0*X(J+1)-1.D0-600.D0* (X(J)**2))
          A(J+1) = 200.0D0
   40 CONTINUE

      RETURN

      END
C***********************************************************************

      SUBROUTINE ROSHDP(N,D,HD,Y1,Y2)

C -------------------------------------------------
C COMPUTE
c               HD =  H*D
c where H is the Hessian matrix and D is a vector.
C -------------------------------------------------


C     .. Parameters ..
      INTEGER MAXN
      PARAMETER (MAXN=1000)
C     ..
C     .. Scalar Arguments ..
      INTEGER N
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION D(N),HD(N),Y1(N),Y2(N)
C     ..
C     .. Arrays in Common ..
      DOUBLE PRECISION HESD(MAXN),HESOD(MAXN)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION OFDIAG
      INTEGER J
C     ..
C     .. Common blocks ..
      COMMON /HESIAN/HESD,HESOD
C     ..
      DO 10 J = 1,N - 1,2
          OFDIAG = HESOD(J)
          HD(J) = HESD(J)*D(J) + OFDIAG*D(J+1)
          HD(J+1) = OFDIAG*D(J) + HESD(J+1)*D(J+1)
   10 CONTINUE

      RETURN

      END

C***********************************************************************
      SUBROUTINE ROSPAT(N,X,A,IA,JA)

C -------------------------------------------------------------
C      DETERMINE The PATTERN OF M
C -------------------------------------------------------------


C     .. Scalar Arguments ..
      INTEGER N
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION A(*),X(N)
      INTEGER IA(N+1),JA(*)
C     ..
C     .. Local Scalars ..
      INTEGER I
C     ..
      DO 10 I = 1,N
          IA(I) = I
          JA(I) = I
   10 CONTINUE
      IA(N+1) = N + 1

      RETURN

      END
