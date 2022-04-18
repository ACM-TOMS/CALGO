      MODULE REAL_PRECISION
! This is for 64-bit arithmetic.
      INTEGER, PARAMETER:: R8=SELECTED_REAL_KIND(13)
      END MODULE REAL_PRECISION


!  This module provides global allocatable arrays used for the sparse
!  matrix data structures, and by the polynomial system solver.  The
!  MODULE HOMOTOPY uses this module.
!
      MODULE HOMPACK90_GLOBAL
      USE REAL_PRECISION
      INTEGER, DIMENSION(:), ALLOCATABLE:: COLPOS, IPAR, ROWPOS
      REAL (KIND=R8), DIMENSION(:), ALLOCATABLE:: PAR, PP, QRSPARSE
      END MODULE HOMPACK90_GLOBAL


      MODULE HOMOTOPY       ! Interfaces for user written subroutines.
      USE REAL_PRECISION, ONLY : R8
      USE HOMPACK90_GLOBAL
!
C Interface for subroutine that evaluates F(X) and returns it in the vector V.
      INTERFACE
         SUBROUTINE F(X,V)
         USE REAL_PRECISION
         REAL (KIND=R8), DIMENSION(:), INTENT(IN):: X
         REAL (KIND=R8), DIMENSION(:), INTENT(OUT):: V
         END SUBROUTINE F
      END INTERFACE
!
C Interface for subroutine that returns in V the K-th column of the Jacobian 
C matrix of F(X) evaluated at X. 
      INTERFACE
         SUBROUTINE FJAC(X,V,K)
         USE REAL_PRECISION
         REAL (KIND=R8), DIMENSION(:), INTENT(IN):: X
         REAL (KIND=R8), DIMENSION(:), INTENT(OUT):: V
         INTEGER, INTENT(IN):: K
         END SUBROUTINE FJAC
      END INTERFACE
!
C Interface for subroutine that evaluates RHO(A,LAMBDA,X) and returns it 
C in the vector V.
      INTERFACE
         SUBROUTINE RHO(A,LAMBDA,X,V)
         USE REAL_PRECISION
         REAL (KIND=R8), INTENT(IN):: A(:),X(:)
         REAL (KIND=R8), INTENT(IN OUT):: LAMBDA
         REAL (KIND=R8), INTENT(OUT):: V(:)
         END SUBROUTINE RHO
      END INTERFACE
C The following code is specifically for the polynomial system driver
C POLSYS1H, and should be used verbatim with POLSYS1H in the external 
C subroutine RHO.  
!     USE HOMPACK90_GLOBAL, ONLY: IPAR, PAR  ! FOR POLSYS1H ONLY.
!     INTERFACE
!       SUBROUTINE HFUNP(N,A,LAMBDA,X)
!       USE REAL_PRECISION
!       INTEGER, INTENT(IN):: N
!       REAL (KIND=R8), INTENT(IN):: A(2*N),LAMBDA,X(2*N)
!       END SUBROUTINE HFUNP
!     END INTERFACE
!     INTEGER:: J,NPOL
C FORCE PREDICTED POINT TO HAVE  LAMBDA .GE. 0  .
!     IF (LAMBDA .LT. 0.0) LAMBDA=0.0
!     NPOL=IPAR(1)
!     CALL HFUNP(NPOL,A,LAMBDA,X)
!     DO J=1,2*NPOL
!       V(J)=PAR(IPAR(3 + (4-1)) + (J-1))
!     END DO
!     RETURN
C If calling FIXP?? or STEP?? directly, supply appropriate replacement
C code in the external subroutine RHO.
!
C Interface for subroutine that calculates and returns in A the vector
C Z such that RHO(Z,LAMBDA,X) = 0 .
      INTERFACE
         SUBROUTINE RHOA(A,LAMBDA,X)
         USE REAL_PRECISION
         REAL (KIND=R8), DIMENSION(:), INTENT(OUT):: A
         REAL (KIND=R8), INTENT(IN):: LAMBDA,X(:)
         END SUBROUTINE RHOA
      END INTERFACE
!
C Interface for subroutine that returns in the vector V the Kth column
C of the Jacobian matrix [D RHO/D LAMBDA, D RHO/DX] evaluated at the
C point (A, LAMBDA, X).
      INTERFACE
         SUBROUTINE RHOJAC(A,LAMBDA,X,V,K)
         USE REAL_PRECISION
         REAL (KIND=R8), INTENT(IN):: A(:),X(:)
         REAL (KIND=R8), INTENT(IN OUT):: LAMBDA
         REAL (KIND=R8), INTENT(OUT):: V(:)
         INTEGER, INTENT(IN):: K
         END SUBROUTINE RHOJAC
      END INTERFACE
C The following code is specifically for the polynomial system driver
C POLSYS1H, and should be used verbatim with POLSYS1H in the external 
C subroutine RHOJAC.  
!     USE HOMPACK90_GLOBAL, ONLY: IPAR, PAR  ! FOR POLSYS1H ONLY.
!     INTERFACE
!       SUBROUTINE HFUNP(N,A,LAMBDA,X)
!       USE REAL_PRECISION
!       INTEGER, INTENT(IN):: N
!       REAL (KIND=R8), INTENT(IN):: A(2*N),LAMBDA,X(2*N)
!       END SUBROUTINE HFUNP
!     END INTERFACE
!     INTEGER:: J,NPOL,N2
!     NPOL=IPAR(1)
!     N2=2*NPOL
!     IF (K .EQ. 1) THEN
C FORCE PREDICTED POINT TO HAVE  LAMBDA .GE. 0  .
!       IF (LAMBDA .LT. 0.0) LAMBDA=0.0
!       CALL HFUNP(NPOL,A,LAMBDA,X)
!       DO J=1,N2
!         V(J)=PAR(IPAR(3 + (6-1)) + (J-1))
!       END DO
!       RETURN
!     ELSE
!       DO J=1,N2
!         V(J)=PAR(IPAR(3 + (5-1)) + (J-1) + N2*(K-2))
!       END DO
!     ENDIF
C
!     RETURN
C If calling FIXP?? or STEP?? directly, supply appropriate replacement
C code in the external subroutine RHOJAC.
!
!
C Interface for subroutine that evaluates a sparse Jacobian matrix of
C F(X) at X, and operates as follows:
C
C If MODE = 1,
C evaluate the N x N symmetric Jacobian matrix of F(X) at X, and return
C the result in packed skyline storage format in QRSPARSE.  LENQR is the
C length of QRSPARSE, and ROWPOS contains the indices of the diagonal
C elements of the Jacobian matrix within QRSPARSE.  ROWPOS(N+1) and
C ROWPOS(N+2) are set by subroutine FODEDS.  The allocatable array COLPOS
C is not used by this storage format.
C
C If MODE = 2,
C evaluate the N x N Jacobian matrix of F(X) at X, and return the result
C in sparse row storage format in QRSPARSE.  LENQR is the length of
C QRSPARSE, ROWPOS contains the indices of where each row begins within
C QRSPARSE, and COLPOS (of length LENQR) contains the column indices of
C the corresponding elements in QRSPARSE.  Even if zero, the diagonal
C elements of the Jacobian matrix must be stored in QRSPARSE.
      INTERFACE
         SUBROUTINE FJACS(X)
         USE REAL_PRECISION
         USE HOMPACK90_GLOBAL, ONLY: QRSPARSE, ROWPOS, COLPOS
         REAL (KIND=R8), DIMENSION(:), INTENT(IN):: X
         END SUBROUTINE FJACS
      END INTERFACE
!
!
C Interface for subroutine that evaluates a sparse Jacobian matrix of
C RHO(A,X,LAMBDA) at (A,X,LAMBDA), and operates as follows:
C
C If MODE = 1,
C evaluate the N X N symmetric Jacobian matrix [D RHO/DX] at
C (A,X,LAMBDA), and return the result in packed skyline storage format in
C QRSPARSE.  LENQR is the length of QRSPARSE, and ROWPOS contains the
C indices of the diagonal elements of [D RHO/DX] within QRSPARSE.  PP
C contains -[D RHO/D LAMBDA] evaluated at (A,X,LAMBDA).  Note the minus
C sign in the definition of PP.  The allocatable array COLPOS is not used
C in this storage format.
C
C If MODE = 2,
C evaluate the N X (N+1) Jacobian matrix [D RHO/DX, D RHO/DLAMBDA] at
C (A,X,LAMBDA), and return the result in sparse row storage format in
C QRSPARSE.  LENQR is the length of QRSPARSE, ROWPOS contains the indices
C of where each row begins within QRSPARSE, and COLPOS (of length LENQR)
C contains the column indices of the corresponding elements in QRSPARSE.
C Even if zero, the diagonal elements of the Jacobian matrix must be
C stored in QRSPARSE.  The allocatable array PP is not used in this
C storage format.
C
      INTERFACE
         SUBROUTINE RHOJS(A,LAMBDA,X)
         USE REAL_PRECISION
         USE HOMPACK90_GLOBAL, ONLY: QRSPARSE, ROWPOS, COLPOS
         REAL (KIND=R8), INTENT(IN):: A(:),LAMBDA,X(:)
         END SUBROUTINE RHOJS
      END INTERFACE
      END MODULE HOMOTOPY


      MODULE HOMPACK90
!
!  This MODULE is an encapsulation of the HOMPACK90 drivers, and uses
!  the modules REAL_PRECISION (defines real precision for all
!  routines), HOMPACK90_GLOBAL (allocatable global data structures for
!  sparse matrices), and HOMOTOPY (interfaces to user written routines
!  defining the problem).
!
!  The intended usage is that the calling program would include a
!  statement like, for example,
!     USE HOMPACK90, ONLY : FIXPNF
!
      USE REAL_PRECISION, ONLY : R8    ! Kind for all reals.
      USE HOMPACK90_GLOBAL             ! Allocated data structures.
      CONTAINS
!
      SUBROUTINE FIXPDF(N,Y,IFLAG,ARCTOL,EPS,TRACE,A,NDIMA,
     &   NFE,ARCLEN)
C
C Subroutine  FIXPDF  finds a fixed point or zero of the
C N-dimensional vector function F(X), or tracks a zero curve
C of a general homotopy map RHO(A,LAMBDA,X).  For the fixed 
C point problem F(X) is assumed to be a C2 map of some ball 
C into itself.  The equation  X = F(X)  is solved by
C following the zero curve of the homotopy map
C
C  LAMBDA*(X - F(X)) + (1 - LAMBDA)*(X - A)  ,
C
C starting from LAMBDA = 0, X = A.  The curve is parameterized
C by arc length S, and is followed by solving the ordinary
C differential equation  D(HOMOTOPY MAP)/DS = 0  for
C Y(S) = (LAMBDA(S), X(S)).
C
C For the zero finding problem F(X) is assumed to be a C2 map
C such that for some R > 0,  X*F(X) >= 0  whenever NORM(X) = R.
C The equation  F(X) = 0  is solved by following the zero curve
C of the homotopy map
C
C   LAMBDA*F(X) + (1 - LAMBDA)*(X - A)
C
C emanating from LAMBDA = 0, X = A.
C
C  A  must be an interior point of the above mentioned balls.
C
C For the curve tracking problem RHO(A,LAMBDA,X) is assumed to
C be a C2 map from E**M X [0,1) X E**N into E**N, which for
C almost all parameter vectors A in some nonempty open subset
C of E**M satisfies
C
C  rank [D RHO(A,LAMBDA,X)/D LAMBDA , D RHO(A,LAMBDA,X)/DX] = N
C
C for all points (LAMBDA,X) such that RHO(A,LAMBDA,X)=0.  It is
C further assumed that
C
C           rank [ D RHO(A,0,X0)/DX ] = N  .
C
C With A fixed, the zero curve of RHO(A,LAMBDA,X) emanating
C from  LAMBDA = 0, X = X0  is tracked until  LAMBDA = 1  by
C solving the ordinary differential equation
C D RHO(A,LAMBDA(S),X(S))/DS = 0  for  Y(S) = (LAMBDA(S), X(S)),
C where S is arc length along the zero curve.  Also the homotopy
C map RHO(A,LAMBDA,X) is assumed to be constructed such that
C
C              D LAMBDA(0)/DS > 0  .
C
C This code is based on the algorithm in L. T. Watson, A
C globally convergent algorithm for computing fixed points of
C C2 maps, Appl. Math. Comput., 5 (1979) 297-311.
C
C
C For the fixed point and zero finding problems, the user
C must supply a subroutine  F(X,V)  which evaluates F(X) at X
C and returns the vector F(X) in V, and a subroutine  FJAC(X,V,K)
C which returns in V the Kth column of the Jacobian matrix of 
C F(X) evaluated at X.  For the curve tracking problem, the user must
C supply a subroutine  RHOA(V,LAMBDA,X)  which given 
C (LAMBDA,X) returns a parameter vector A in V such that 
C RHO(A,LAMBDA,X)=0, and a subroutine  RHOJAC(A,LAMBDA,X,V,K)
C which returns in V the Kth column of the N X (N+1) Jacobian 
C matrix [D RHO/D LAMBDA, D RHO/DX] evaluated at (A,LAMBDA,X).
C Whichever of the routines  F,  FJAC,  RHOA,  RHOJAC  are required
C should be supplied as external subroutines, conforming with the
C interfaces in the module  HOMOTOPY.
C FIXPDF  directly or indirectly uses the subroutines
C   F (or  RHOA ),  FJAC (or  RHOJAC ),  FODE,   ROOT,
C   SINTRP,  STEPS,  the LAPACK routine  DGEQPF,  auxiliary LAPACK 
C routines, and the BLAS functions  DCOPY,  DDOT,  DGEMV,  DGER,  
C   DNRM2,  DSCAL,  DSWAP,  IDAMAX.  
C The module  REAL_PRECISION  specifies 64-bit real arithmetic, which
C the user may want to change.
C
C ***Warning:  this subroutine is generally more robust than  FIXPNF
C and  FIXPQF, but may be slower than those subroutines by a
C factor of two.
C
C
C ON INPUT:
C
C N  is the dimension of X, F(X), and RHO(A,LAMBDA,X).
C
C Y  is an array of length  N + 1.  (Y(2),...,Y(N+1)) = A  is the
C    starting point for the zero curve for the fixed point and 
C    zero finding problems.  (Y(2),...,Y(N+1)) = X0  for the curve
C    tracking problem.
C
C IFLAG  can be -2, -1, 0, 2, or 3.  IFLAG  should be 0 on the 
C    first call to  FIXPDF  for the problem  X=F(X), -1 for the
C    problem  F(X)=0, and -2 for the problem  RHO(A,LAMBDA,X)=0.
C    In certain situations  IFLAG  is set to 2 or 3 by  FIXPDF,
C    and  FIXPDF  can be called again without changing  IFLAG.
C
C ARCTOL  is the local error allowed the ODE solver when
C    following the zero curve.  If  ARCTOL .LE. 0.0  on input
C    it is reset to  .5*SQRT(EPS).  Normally  ARCTOL  should
C    be considerably larger than  EPS.
C
C EPS  is the local error allowed the ODE solver when very
C    near the fixed point(zero).  EPS  is approximately the
C    mixed absolute and relative error in the computed fixed 
C    point(zero).
C
C TRACE  is an integer specifying the logical I/O unit for
C    intermediate output.  If  TRACE .GT. 0  the points computed on
C    the zero curve are written to I/O unit  TRACE .
C
C A(1:NDIMA) contains the parameter vector  A.  For the fixed point
C    and zero finding problems, A  need not be initialized by the
C    user, and is assumed to have length  N.  For the curve
C    tracking problem, A  has length  NDIMA  and must be initialized
C    by the user.
C
C NDIMA  is the dimension of  A, used for the curve tracking problem,
C    and must be N for the fixed point and zero finding problems.
C
C Y, ARCTOL, EPS, ARCLEN, NFE, and IFLAG should all be
C variables in the calling program.
C
C
C ON OUTPUT:
C
C N  and  TRACE  are unchanged.
C
C Y(1) = LAMBDA, (Y(2),...,Y(N+1)) = X, and Y is an approximate
C    zero of the homotopy map.  Normally LAMBDA = 1 and X is a
C    fixed point(zero) of F(X).  In abnormal situations LAMBDA
C    may only be near 1 and X is near a fixed point(zero).
C
C IFLAG =
C  -2   causes  FIXPDF  to initialize everything for the problem
C       RHO(A,LAMBDA,X) = 0 (use on first call).
C
C  -1   causes  FIXPDF  to initialize everything for the problem
C       F(X) = 0 (use on first call).
C
C   0   causes  FIXPDF  to initialize everything for the problem
C       X = F(X) (use on first call).
C
C   1   Normal return.
C
C   2   Specified error tolerance cannot be met.  EPS has been
C       increased to a suitable value.  To continue, just call
C       FIXPDF  again without changing any parameters.
C
C   3   STEPS  has been called 1000 times.  To continue, call
C       FIXPDF  again without changing any parameters.
C
C   4   Jacobian matrix does not have full rank.  The algorithm
C       has failed (the zero curve of the homotopy map cannot be
C       followed any further).
C
C   5   EPS  (or  ARCTOL) is too large.  The problem should be
C       restarted by calling  FIXPDF  with a smaller  EPS  (or
C       ARCTOL) and  IFLAG = 0 (-1, -2).
C
C   6   I - DF(X)  is nearly singular at the fixed point (DF(X) is
C       nearly singular at the zero, or  D RHO(A,LAMBDA,X)/DX  is
C       nearly singular at  LAMBDA = 1 ).  Answer may not be
C       accurate.
C
C   7   Illegal input parameters, a fatal error.
C
C   8   Memory allocation error, fatal.
C
C ARCTOL = EPS after a normal return (IFLAG = 1).
C
C EPS  is unchanged after a normal return (IFLAG = 1).  It is
C    increased to an appropriate value on the return IFLAG = 2.
C
C A  will (normally) have been modified.
C
C NFE  is the number of function evaluations (= number of
C    Jacobian matrix evaluations).
C
C ARCLEN  is the length of the path followed.
C
C
C Automatic work arrays:
C
C YP(1:N+1) is a work array containing the current tangent
C    vector to the zero curve.
C
C YPOLD(1:N+1) is a work array containing the previous tangent
C    vector to the zero curve.
C
C QR(1:N,1:N+1), ALPHA(1:3*N+3), TZ(1:N+1), and PIVOT(1:N+1) are 
C    all work arrays used by  FODE  to calculate the tangent
C    vector YP.
C
C WT(1:N+1), PHI(1:N+1,1:16), and P(1:N+1) are all work arrays
C    used by the ODE subroutine  STEPS  .
C
      USE HOMOTOPY
      USE REAL_PRECISION
C
      INTEGER, INTENT(IN)::N,NDIMA,TRACE
      REAL (KIND=R8), DIMENSION(:), INTENT(IN OUT)::A,Y
      INTEGER, INTENT(IN OUT)::IFLAG
      REAL (KIND=R8), INTENT(IN OUT)::ARCTOL,EPS
      INTEGER, INTENT(OUT)::NFE
      REAL (KIND=R8), INTENT(OUT)::ARCLEN
C
C LOCAL VARIABLES.
      REAL (KIND=R8), SAVE:: CURSW,CURTOL,EPSSTP,EPST,H,HOLD,
     &  S,S99,SA,SB,SOUT,SQNP1,XOLD,Y1SOUT
      INTEGER, SAVE:: IFLAGC,ITER,IVC,JW,K,KGI,KOLD,
     &  KSTEPS,LCODE,LIMIT,NFEC,NP1
      LOGICAL, SAVE:: CRASH,START,ST99
C
C *****  ARRAY DECLARATIONS.  *****
C
C ARRAYS NEEDED BY THE ODE SUBROUTINE  STEPS .
      REAL (KIND=R8), ALLOCATABLE, SAVE:: P(:),PHI(:,:),WT(:),YP(:)
      REAL (KIND=R8), SAVE:: ALPHAS(12),G(13),GI(11),W(12)
      INTEGER, SAVE:: IV(10)
C
C ARRAYS NEEDED BY  FIXPDF , FODE , AND LAPACK ROUTINES.
      REAL (KIND=R8), DIMENSION(:), ALLOCATABLE, SAVE:: YPOLD
      REAL (KIND=R8):: ALPHA(3*N+3),AOLD(NDIMA),QR(N,N+1),TZ(N+1)
      INTEGER:: PIVOT(N+1)
C
C *****  END OF DIMENSIONAL INFORMATION.  *****
C
C LIMITD  IS AN UPPER BOUND ON THE NUMBER OF STEPS.  IT MAY BE
C CHANGED BY CHANGING THE FOLLOWING PARAMETER STATEMENT:
      INTEGER, PARAMETER:: LIMITD=1000
C
      INTERFACE
        SUBROUTINE FODE(S,Y,YP,YPOLD,A,QR,ALPHA,TZ,PIVOT,NFE,N,IFLAG)
        USE REAL_PRECISION
        REAL (KIND=R8):: S
        INTEGER:: IFLAG,N,NFE
        REAL (KIND=R8):: A(:),Y(:),YP(N+1),YPOLD(N+1)
        REAL (KIND=R8):: ALPHA(3*N+3),QR(N,N+1),TZ(N+1)
        INTEGER, DIMENSION(N+1):: PIVOT
        END SUBROUTINE FODE
        SUBROUTINE STEPS(F,NEQN,Y,X,H,EPS,WT,START,HOLD,K,KOLD,CRASH,
     &  PHI,P,YP,ALPHA,W,G,KSTEPS,XOLD,IVC,IV,KGI,GI,  FPWA1,FPWA2,
     &  FPWA3,FPWA4,FPWA5,IFPWA1,IFPC1,IFPC2)
        USE REAL_PRECISION
        EXTERNAL F
        REAL (KIND=R8):: ALPHA,EPS,FPWA1,FPWA2,FPWA3,FPWA4,FPWA5,
     &  G,GI,H,HOLD,P,PHI,W,WT,X,XOLD,Y,YP
        INTEGER:: IFPC1,IFPC2,IFPWA1,IV,IVC,K,KGI,KOLD,KSTEPS,NEQN
        LOGICAL:: CRASH,START
        DIMENSION Y(:),WT(NEQN),PHI(NEQN,16),P(NEQN),YP(NEQN),
     &  ALPHA(12),W(12),G(13),GI(11),IV(10),
     &  FPWA1(NEQN),FPWA2(:),FPWA3(NEQN-1,NEQN),FPWA4(3*NEQN),
     &  FPWA5(NEQN),IFPWA1(NEQN)
        END SUBROUTINE STEPS
      END INTERFACE
C
C
C :  :  :  :  :  :  :  :  :  :  :  :  :  :  :  :  :  :  :  :  :
      IF (N .LE. 0  .OR.  EPS .LE. 0.0  .OR.  N+1 .NE. SIZE(Y)
     &  .OR.  NDIMA .NE. SIZE(A)  .OR.
     &  ((IFLAG .EQ. -1  .OR.  IFLAG .EQ. 0) .AND.  N .NE. SIZE(A)))
     &  IFLAG=7
      IF (IFLAG .GE. -2  .AND.  IFLAG .LE. 0) GO TO 10
      IF (IFLAG .EQ. 2) GO TO 35
      IF (IFLAG .EQ. 3) GO TO 30
C ONLY VALID INPUT FOR  IFLAG  IS -2, -1, 0, 2, 3.
      IFLAG=7
      RETURN
C
C *****  INITIALIZATION BLOCK.  *****
C
10    ARCLEN=0.0
      S=0.0
      IF (ARCTOL .LE. 0.0) ARCTOL=.5*SQRT(EPS)
      NFEC=0
      IFLAGC=IFLAG
      NP1=N+1
      SQNP1=SQRT(DBLE(NP1))
      IF (ALLOCATED(P)) DEALLOCATE(P)
      IF (ALLOCATED(PHI)) DEALLOCATE(PHI)
      IF (ALLOCATED(WT)) DEALLOCATE(WT)
      IF (ALLOCATED(YP)) DEALLOCATE(YP)
      IF (ALLOCATED(YPOLD)) DEALLOCATE(YPOLD)
      ALLOCATE(P(NP1),PHI(NP1,16),WT(NP1),YP(NP1),YPOLD(NP1),
     &  STAT=JW)
      IF (JW /= 0) THEN
        IFLAG=8
        RETURN
      END IF
C
C SWITCH FROM THE TOLERANCE  ARCTOL  TO THE (FINER) TOLERANCE  EPS  IF
C THE CURVATURE OF ANY COMPONENT OF  Y  EXCEEDS  CURSW.
C
      CURSW=10.0
C
      ST99=.FALSE.
      START=.TRUE.
      CRASH=.FALSE.
      HOLD=1.0
      H=.1
      EPSSTP=ARCTOL
      KSTEPS=0
C SET INITIAL CONDITIONS FOR ORDINARY DIFFERENTIAL EQUATION.
      YPOLD(1)=1.0
      YP(1)=1.0
      Y(1)=0.0
      YPOLD(2:NP1)=0.0
      YP(2:NP1)=0.0
C LOAD  A  FOR THE FIXED POINT AND ZERO FINDING PROBLEMS.
      IF (IFLAGC .GE. -1) THEN
        A=Y(2:NP1)
      ENDIF
30    LIMIT=LIMITD
C
C *****  END OF INITIALIZATION BLOCK.  *****
C
35    MAIN_LOOP: DO ITER=1,LIMIT  ! *****  MAIN LOOP.  *****
      IF (Y(1) .LT. 0.0) THEN
        ARCLEN=ARCLEN+S
        IFLAG=5
        CALL CLEANUP ; RETURN
      ENDIF
      IF (S .LE. 7.0*SQNP1) GO TO 80
C ARC LENGTH IS GETTING TOO LONG, THE PROBLEM WILL BE
C RESTARTED WITH A DIFFERENT  A  VECTOR.
      ARCLEN=ARCLEN+S
      S=0.0
60    START=.TRUE.
      CRASH=.FALSE.
C COMPUTE A NEW  A  VECTOR.
      IF (IFLAGC .EQ. -2) THEN
        AOLD=A
        CALL RHOA(A,Y(1),Y(2:NP1))
C IF NEW AND OLD  A  DIFFER BY TOO MUCH, TRACKING SHOULD NOT CONTINUE.
        IF (ANY(ABS(A-AOLD) .GT. 1.0+ABS(AOLD))) THEN
          ARCLEN=ARCLEN+S
          IFLAG=5
          CALL CLEANUP ; RETURN
        ENDIF
      ELSE
        CALL F(Y(2:NP1),YP(1:N))
        AOLD=A
        IF (IFLAGC .EQ. -1) THEN
          A=Y(1)*YP(1:N)/(1.0 - Y(1)) + Y(2:NP1)
        ELSE
          A=(Y(2:NP1) - Y(1)*YP(1:N))/(1.0 - Y(1))
        ENDIF
C IF NEW AND OLD  A  DIFFER BY TOO MUCH, TRACKING SHOULD NOT CONTINUE.
        IF (ANY(ABS(A-AOLD) .GT. 1.0+ABS(AOLD))) THEN
          ARCLEN=ARCLEN+S
          IFLAG=5
          CALL CLEANUP ; RETURN
        ENDIF
      ENDIF
      GO TO 100
80    IF (Y(1) .LE. .99  .OR. ST99) GO TO 100
C WHEN LAMBDA REACHES .99, THE PROBLEM WILL BE RESTARTED WITH
C A NEW  A  VECTOR.
90    ST99=.TRUE.
      EPSSTP=EPS
      ARCTOL=EPS
      GO TO 60
C
C SET DIFFERENT ERROR TOLERANCE FOR HIGH CURVATURE COMPONENTS OF THE
C TRAJECTORY Y(S).
100   CURTOL=CURSW*HOLD
      EPST=EPS/EPSSTP
      WHERE (ABS(YP-YPOLD) .LE. CURTOL)
        WT=(ABS(Y)+1.0)
      ELSEWHERE
        WT=(ABS(Y)+1.0)*EPST
      END WHERE
C
C TAKE A STEP ALONG THE CURVE.
      CALL STEPS(FODE,NP1,Y,S,H,EPSSTP,WT,START,HOLD,K,KOLD,CRASH,
     &     PHI,P,YP,ALPHAS,W,G,KSTEPS,XOLD,IVC,IV,KGI,GI,
     &     YPOLD,A,QR,ALPHA,TZ,PIVOT,NFEC,IFLAGC)
C PRINT LATEST POINT ON CURVE IF REQUESTED.
      IF (TRACE .GT. 0) THEN
        WRITE (TRACE,117) ITER,NFEC,S,Y(1),(Y(JW),JW=2,NP1)
117     FORMAT(/' STEP',I5,3X,'NFE =',I5,3X,'ARC LENGTH =',F9.4,3X,
     &  'LAMBDA =',F7.4,5X,'X VECTOR:'/(1X,6ES12.4))
      ENDIF
      NFE=NFEC
C CHECK IF THE STEP WAS SUCCESSFUL.
      IF (IFLAGC .EQ. 4) THEN
        ARCLEN=ARCLEN+S
        IFLAG=4
        CALL CLEANUP ; RETURN
      ENDIF
      IF (CRASH) THEN
C RETURN CODE FOR ERROR TOLERANCE TOO SMALL.
        IFLAG=2
C CHANGE ERROR TOLERANCES.
        EPS=EPSSTP
        IF (ARCTOL .LT. EPS) ARCTOL=EPS
C CHANGE LIMIT ON NUMBER OF ITERATIONS.
        LIMIT=LIMIT-ITER
        RETURN
      ENDIF
C
      IF (Y(1) .GE. 1.0) THEN
        IF (ST99) GO TO 160
C
C IF LAMBDA .GE. 1.0 BUT THE PROBLEM HAS NOT BEEN RESTARTED
C WITH A NEW  A  VECTOR, BACK UP AND RESTART.
C
        S99=S-.5*HOLD
C GET AN APPROXIMATE ZERO Y(S) WITH  Y(1)=LAMBDA .LT. 1.0  .
135     CALL SINTRP(S,Y,S99,WT,YP,NP1,KOLD,PHI,IVC,IV,KGI,GI,
     &     ALPHAS,G,W,XOLD,P)
        IF (WT(1) .LT. 1.0) GO TO 140
        S99=.5*(S-HOLD+S99)
        GO TO 135
C
140     Y=WT
        YPOLD=YP
        S=S99
        GO TO 90
      ENDIF
C
      END DO MAIN_LOOP  ! *****  END OF MAIN LOOP.  *****
C
C LAMBDA HAS NOT REACHED 1 IN 1000 STEPS.
      IFLAG=3
      RETURN
C
C
C USE INVERSE INTERPOLATION TO GET THE ANSWER AT LAMBDA = 1.0 .
C
160   SA=S-HOLD
      SB=S
      LCODE=1
170   CALL ROOT(SOUT,Y1SOUT,SA,SB,EPS,EPS,LCODE)
C ROOT  FINDS S SUCH THAT Y(1)(S) = LAMBDA = 1 .
      IF (LCODE .GT. 0) GO TO 190
      CALL SINTRP(S,Y,SOUT,WT,YP,NP1,KOLD,PHI,IVC,IV,KGI,GI,
     &     ALPHAS,G,W,XOLD,P)
      Y1SOUT=WT(1)-1.0
      GO TO 170
190   IFLAG=1
C SET IFLAG = 6 IF  ROOT  COULD NOT GET  LAMBDA = 1.0  .
      IF (LCODE .GT. 2) IFLAG=6
      ARCLEN=ARCLEN+SA
C LAMBDA(SA) = 1.0 .
      CALL SINTRP(S,Y,SA,WT,YP,NP1,KOLD,PHI,IVC,IV,KGI,GI,
     &     ALPHAS,G,W,XOLD,P)
C
      Y=WT
      CALL CLEANUP ; RETURN
C
      CONTAINS
        SUBROUTINE CLEANUP
        IF (ALLOCATED(P)) DEALLOCATE(P)
        IF (ALLOCATED(PHI)) DEALLOCATE(PHI)
        IF (ALLOCATED(WT)) DEALLOCATE(WT)
        IF (ALLOCATED(YP)) DEALLOCATE(YP)
        IF (ALLOCATED(YPOLD)) DEALLOCATE(YPOLD)
        END SUBROUTINE CLEANUP
      END SUBROUTINE FIXPDF
!
      SUBROUTINE FIXPDS(N,Y,IFLAG,ARCTOL,EPS,TRACE,A,NDIMA,NFE,
     &     ARCLEN,MODE,LENQR)
C
C Subroutine  FIXPDS  finds a fixed point or zero of the
C N-dimensional vector function F(X), or tracks a zero curve
C of a general homotopy map RHO(A,X,LAMBDA).  For the fixed 
C point problem F(X) is assumed to be a C2 map of some ball 
C into itself.  The equation  X = F(X)  is solved by
C following the zero curve of the homotopy map
C
C  LAMBDA*(X - F(X)) + (1 - LAMBDA)*(X - A)  ,
C
C starting from LAMBDA = 0, X = A.  The curve is parameterized
C by arc length S, and is followed by solving the ordinary
C differential equation  D(HOMOTOPY MAP)/DS = 0  for
C Y(S) = (X(S), LAMBDA(S)).
C
C For the zero finding problem F(X) is assumed to be a C2 map
C such that for some R > 0,  X*F(X) >= 0  whenever NORM(X) = R.
C The equation  F(X) = 0  is solved by following the zero curve
C of the homotopy map
C
C   LAMBDA*F(X) + (1 - LAMBDA)*(X - A)
C
C emanating from LAMBDA = 0, X = A.
C
C  A  must be an interior point of the above mentioned balls.
C
C For the curve tracking problem RHO(A,X,LAMBDA) is assumed to
C be a C2 map from E**M X E**N X [0,1) into E**N, which for
C almost all parameter vectors A in some nonempty open subset
C of E**M satisfies
C
C  rank [D RHO(A,X,LAMBDA)/D LAMBDA , D RHO(A,X,LAMBDA)/DX] = N
C
C for all points (X,LAMBDA) such that RHO(A,X,LAMBDA)=0.  It is
C further assumed that
C
C           rank [ D RHO(A,X0,0)/DX ] = N  .
C
C With A fixed, the zero curve of RHO(A,X,LAMBDA) emanating
C from  LAMBDA = 0, X = X0  is tracked until  LAMBDA = 1  by
C solving the ordinary differential equation
C D RHO(A,X(S),LAMBDA(S))/DS = 0  for  Y(S) = (X(S), LAMBDA(S)),
C where S is arc length along the zero curve.  Also the homotopy
C map RHO(A,X,LAMBDA) is assumed to be constructed such that
C
C              D LAMBDA(0)/DS > 0  .
C
C This code is based on the algorithm in L. T. Watson, A
C globally convergent algorithm for computing fixed points of
C C2 maps, Appl. Math. Comput., 5 (1979) 297-311.
C
C
C For the fixed point and zero finding problems, the user
C must supply a subroutine  F(X,V)  which evaluates F(X) at X
C and returns the vector F(X) in V, and a subroutine
C  FJACS(X)  which evaluates, if
C MODE = 1,
C   the (symmetric) Jacobian matrix of F(X) at X, and returns the
C   symmetric Jacobian matrix in packed skyline storage format in
C   QR, or if
C MODE = 2,
C   returns the (nonsymmetric) Jacobian matrix in sparse row format
C   in QR.  The MODE 1 format is defined by QR, LENQR, ROWPOS; the
C   MODE 2 format is defined by QR, LENQR, ROWPOS, COLPOS.
C
C For the curve tracking problem, the user must supply a subroutine
C  RHOA(V,LAMBDA,X)  which given (X,LAMBDA) returns a
C parameter vector A in V such that RHO(A,X,LAMBDA)=0, and a 
C subroutine  RHOJS(A,LAMBDA,X)  which, if
C MODE = 1,
C   returns in QR the symmetric N X N Jacobian matrix [D RHO/DX] 
C   evaluated at (A,X,LAMBDA) and stored in packed skyline format, 
C   and returns in PP the vector -(D RHO/D LAMBDA) evaluated at 
C   (A,X,LAMBDA).  This data structure is described by QR, LENQR,
C   ROWPOS, PP.  *** Note the minus sign in the definition of PP. ***  If
C MODE = 2,
C   the (nonsymmetric) N X (N+1) Jacobian matrix [D RHO/DX, D RHO/DLAMBDA]
C   evaluated at (A,X,LAMBDA) is returned in a data structure described
C   by QR, LENQR, ROWPOS, COLPOS.
C
C Whichever of the routines  F,  FJACS,  RHOA,  RHOJS  are required
C should be supplied as external subroutines, conforming with the
C interfaces in the module  HOMOTOPY.
C
C
C Functions and subroutines directly or indirectly called by FIXPDS:
C DLAIC1  and  DLAMCH (LAPACK), F (or  RHOA ), FJACS (or  RHOJS ),
C FODEDS , GMFADS , GMRES , GMRILUDS , ILUFDS , ILUSOLVDS , MULTDS ,
C MULT2DS , PCGDS , ROOT , SINTRP , SOLVDS , STEPDS , and the BLAS
C functions  DDOT , DNRM2.  The module  REAL_PRECISION  specifies 64-bit
C real arithmetic, which the user may want to change.
C
C ***Warning:  this subroutine is generally more robust than  FIXPNS
C and  FIXPQS, but may be slower than those subroutines by a
C factor of two.
C
C
C ON INPUT:
C
C N  is the dimension of X, F(X), and RHO(A,X,LAMBDA).
C
C Y  is an array of length  N + 1.  (Y(1),...,Y(N)) = A  is the
C    starting point for the zero curve for the fixed point and 
C    zero finding problems.  (Y(1),...,Y(N)) = X0  for the curve
C    tracking problem.
C
C IFLAG  can be -2, -1, 0, 2, or 3.  IFLAG  should be 0 on the 
C    first call to  FIXPDS  for the problem  X=F(X), -1 for the
C    problem  F(X)=0, and -2 for the problem  RHO(A,X,LAMBDA)=0.
C    In certain situations  IFLAG  is set to 2 or 3 by  FIXPDS,
C    and  FIXPDS  can be called again without changing  IFLAG.
C
C ARCTOL  is the local error allowed the ODE solver when
C    following the zero curve.  If  ARCTOL .LE. 0.0  on input
C    it is reset to  .5*SQRT(EPS).  Normally  ARCTOL  should
C    be considerably larger than  EPS.
C
C EPS  is the local error allowed the ODE solver when very
C    near the fixed point(zero).  EPS  is approximately the
C    mixed absolute and relative error in the computed fixed 
C    point(zero).
C
C TRACE  is an integer specifying the logical I/O unit for
C    intermediate output.  If  TRACE .GT. 0  the points computed on
C    the zero curve are written to I/O unit  TRACE .
C
C A(1:NDIMA) contains the parameter vector  A .  For the fixed point
C    and zero finding problems, A  need not be initialized by the
C    user, and is assumed to have length  N.  For the curve
C    tracking problem, A  has length  NDIMA  and must be initialized
C    by the user.
C
C NDIMA  is the dimension of  A, used for the curve tracking problem,
C    and must be N for the fixed point and zero finding problems.
C
C MODE = 1 if the Jacobian matrix is symmetric and stored in a packed
C          skyline format;
C      = 2 if the Jacobian matrix is stored in a sparse row format.
C
C LENQR  is the number of nonzero entries in the sparse Jacobian
C    matrices, used to determine the sparse matrix data structures.
C
C A, Y, ARCTOL, EPS, ARCLEN, NFE, and IFLAG should all be
C variables in the calling program.
C
C
C ON OUTPUT:
C
C N  and  TRACE  are unchanged.
C
C (Y(1),...,Y(N)) = X, Y(N+1) = LAMBDA, and Y is an approximate
C    zero of the homotopy map.  Normally LAMBDA = 1 and X is a
C    fixed point(zero) of F(X).  In abnormal situations LAMBDA
C    may only be near 1 and X is near a fixed point(zero).
C
C IFLAG =
C  -2   causes  FIXPDS  to initialize everything for the problem
C       RHO(A,X,LAMBDA) = 0 (use on first call).
C
C  -1   causes  FIXPDS  to initialize everything for the problem
C       F(X) = 0 (use on first call).
C
C   0   causes  FIXPDS  to initialize everything for the problem
C       X = F(X) (use on first call).
C
C   1   Normal return.
C
C   2   Specified error tolerance cannot be met.  EPS has been
C       increased to a suitable value.  To continue, just call
C       FIXPDS  again without changing any parameters.
C
C   3   STEPDS  has been called 1000 times.  To continue, call
C       FIXPDS  again without changing any parameters.
C
C   4   Jacobian matrix does not have full rank or has a zero on the
C       diagonal, and/or the conjugate gradient iteration for the
C       kernel of the Jacobian matrix failed to converge.  The
C       algorithm has failed (the zero curve of the homotopy map
C       cannot be followed any further).
C
C   5   EPS  (or  ARCTOL ) is too large.  The problem should be
C       restarted by calling  FIXPDS  with a smaller  EPS  (or
C       ARCTOL ) and  IFLAG = 0 (-1, -2).
C
C   6   I - DF(X)  is nearly singular at the fixed point (DF(X) is
C       nearly singular at the zero, or  D RHO(A,X,LAMBDA)/DX  is
C       nearly singular at  LAMBDA = 1 ).  Answer may not be
C       accurate.
C
C   7   Illegal input parameters, a fatal error.
C
C ARCTOL = EPS after a normal return (IFLAG = 1).
C
C EPS  is unchanged after a normal return (IFLAG = 1).  It is
C    increased to an appropriate value on the return IFLAG = 2.
C
C A  will (normally) have been modified.
C
C NFE  is the number of function evaluations (= number of
C    Jacobian evaluations).
C
C ARCLEN  is the length of the path followed.
C
C
C Allocatable and automatic work arrays:
C
C YP(1:N+1) is a work array containing the current tangent
C    vector to the zero curve.
C
C YPOLD(1:N+1) is a work array containing the previous tangent
C    vector to the zero curve.
C
C QR(1:LENQR), PP(1:N), ROWPOS(1:N+2), COLPOS(1:LENQR) are all work
C    arrays used to define the sparse Jacobian matrices, allocated
C    here, and distributed via the module  HOMOTOPY .
C
C WT(1:N+1), PHI(1:N+1,1:16), and P(1:N+1) are all work arrays
C    used by the ODE subroutine  STEPDS  .
C
      USE HOMOTOPY, QR => QRSPARSE
      USE REAL_PRECISION
      INTEGER, INTENT(IN)::LENQR,MODE,N,NDIMA,TRACE
      REAL (KIND=R8), DIMENSION(:), INTENT(IN OUT)::A,Y
      INTEGER, INTENT(IN OUT)::IFLAG
      REAL (KIND=R8), INTENT(IN OUT)::ARCTOL,EPS
      INTEGER, INTENT(OUT)::NFE
      REAL (KIND=R8), INTENT(OUT)::ARCLEN
C
C *****  LOCAL VARIABLES.  *****
C
      REAL (KIND=R8), SAVE:: CURSW,CURTOL,EPSSTP,EPST,
     &  H,HOLD,S,S99,SA,SB,SOUT,SQNP1,XOLD,Y1SOUT
      INTEGER, SAVE:: IFLAGC,ITER,IVC,JW,K,KGI,KOLD,
     &  KSTEPS,LCODE,LIMIT,NFEC,NP1
      LOGICAL, SAVE:: CRASH,START,ST99
C
C ARRAYS NEEDED BY THE ODE SUBROUTINE  STEPDS .
      REAL (KIND=R8), SAVE:: ALPHAS(12),G(13),GI(11),W(12)
      REAL (KIND=R8), ALLOCATABLE, SAVE:: P(:),PHI(:,:),WT(:),YP(:)
      INTEGER, SAVE:: IV(10)
C
C ARRAYS NEEDED BY  FIXPDS , FODEDS , AND  PCGDS .
      REAL (KIND=R8), ALLOCATABLE, DIMENSION(:), SAVE:: AOLD,YPOLD
C
C *****  END OF DIMENSIONAL INFORMATION.  *****
C
      INTERFACE
        SUBROUTINE FODEDS(S,Y,YP,N,IFLAG,YPOLD,A,NDIMA,LENQR,MODE,NFE)
        USE HOMOTOPY, QR => QRSPARSE
        USE REAL_PRECISION
        INTEGER:: IFLAG,LENQR,MODE,N,NDIMA,NFE
        REAL (KIND=R8):: A(NDIMA),S,Y(N+1),YP(N+1),YPOLD(N+1)
        END SUBROUTINE FODEDS
      END INTERFACE
C
C LIMITD  IS AN UPPER BOUND ON THE NUMBER OF STEPS.  IT MAY BE
C CHANGED BY CHANGING THE FOLLOWING PARAMETER STATEMENT:
      INTEGER, PARAMETER:: LIMITD=1000
C
C
C :  :  :  :  :  :  :  :  :  :  :  :  :  :  :  :  :  :  :  :  :
      IF (N .LE. 0  .OR.  EPS .LE. 0.0  .OR.  N+1 .NE. SIZE(Y)
     &  .OR.  NDIMA .NE. SIZE(A)  .OR.
     &  ((IFLAG .EQ. -1  .OR.  IFLAG .EQ. 0) .AND.  N .NE. SIZE(A))
     &  .OR.  MODE .LE. 0  .OR.  MODE .GE. 3)
     &  IFLAG=7
      IF (IFLAG .GE. -2  .AND.  IFLAG .LE. 0) GO TO 10
      IF (IFLAG .EQ. 2) GO TO 35
      IF (IFLAG .EQ. 3) GO TO 30
C ONLY VALID INPUT FOR  IFLAG  IS -2, -1, 0, 2, 3.
      IFLAG=7
      RETURN
C
C *****  INITIALIZATION BLOCK.  *****
C
10    ARCLEN=0.0
      S=0.0
      IF (ARCTOL .LE. 0.0) ARCTOL=.5*SQRT(EPS)
      NFEC=0
      IFLAGC=IFLAG
      NP1=N+1
      SQNP1=SQRT(REAL(NP1,KIND=R8))
C
C SWITCH FROM THE TOLERANCE  ARCTOL  TO THE (FINER) TOLERANCE  EPS  IF
C THE CURVATURE OF ANY COMPONENT OF  Y  EXCEEDS  CURSW.
C
      CURSW=10.0
C
      ST99=.FALSE.
      START=.TRUE.
      CRASH=.FALSE.
      HOLD=1.0
      H=.1
      EPSSTP=ARCTOL
      KSTEPS=0
C ALLOCATE SAVED WORK ARRAYS.
      IF (ALLOCATED(AOLD)) DEALLOCATE(AOLD)
      IF (ALLOCATED(P)) DEALLOCATE(P)
      IF (ALLOCATED(PHI)) DEALLOCATE(PHI)
      IF (ALLOCATED(WT)) DEALLOCATE(WT)
      IF (ALLOCATED(YP)) DEALLOCATE(YP)
      IF (ALLOCATED(YPOLD)) DEALLOCATE(YPOLD)
      ALLOCATE(AOLD(NDIMA),P(NP1),PHI(NP1,16),WT(NP1),YP(NP1),
     &  YPOLD(NP1))
C SET INITIAL CONDITIONS FOR ORDINARY DIFFERENTIAL EQUATION.
      YPOLD(NP1)=1.0
      YP(NP1)=1.0
      Y(NP1)=0.0
      YPOLD(1:N)=0.0
      YP(1:N)=0.0
C LOAD  A  FOR THE FIXED POINT AND ZERO FINDING PROBLEMS.
      IF (IFLAGC .GE. -1) THEN
        A=Y(1:N)
      ENDIF
30    LIMIT=LIMITD
C ALLOCATE ARRAYS FOR SPARSE JACOBIAN MATRIX DATA STRUCTURE.
35    SELECT CASE (MODE)
        CASE (1)
          IF (.NOT. ALLOCATED(QR)) ALLOCATE(QR(LENQR))
          IF (.NOT. ALLOCATED(ROWPOS)) ALLOCATE(ROWPOS(N+2))
          IF (.NOT. ALLOCATED(PP)) ALLOCATE(PP(N))
        CASE (2)
          IF (.NOT. ALLOCATED(QR)) ALLOCATE(QR(LENQR))
          IF (.NOT. ALLOCATED(ROWPOS)) ALLOCATE(ROWPOS(N+2))
          IF (.NOT. ALLOCATED(COLPOS)) ALLOCATE(COLPOS(LENQR))
          IF ((.NOT. ALLOCATED(PP)) .AND. (IFLAGC .GE. -1))
     &      ALLOCATE(PP(N))
      END SELECT
C
C *****  END OF INITIALIZATION BLOCK.  *****
C
      MAIN_LOOP: DO ITER=1,LIMIT  ! *****  MAIN LOOP.  *****
      IF (Y(NP1) .LT. 0.0) THEN
        ARCLEN=ARCLEN+S
        IFLAG=5
        CALL CLEANUPALL
        RETURN
      ENDIF
      IF (S .LE. 7.0*SQNP1) GO TO 80
C ARC LENGTH IS GETTING TOO LONG, THE PROBLEM WILL BE
C RESTARTED WITH A DIFFERENT  A  VECTOR.
      ARCLEN=ARCLEN+S
      S=0.0
60    START=.TRUE.
      CRASH=.FALSE.
C COMPUTE A NEW  A  VECTOR.
      IF (IFLAGC .EQ. -2) THEN
        AOLD=A
        CALL RHOA(A,Y(NP1),Y(1:N))
C IF NEW AND OLD  A  DIFFER BY TOO MUCH, TRACKING SHOULD NOT CONTINUE.
        IF (ANY(ABS(A-AOLD) .GT. 1.0+ABS(AOLD))) THEN
            ARCLEN=ARCLEN+S
            IFLAG=5
            CALL CLEANUPALL
            RETURN
        ENDIF
      ELSE
        CALL F(Y(1:N),YP(1:N))
        AOLD=A
        IF (IFLAGC .EQ. -1) THEN
          A=Y(NP1)*YP(1:N)/(1.0 - Y(NP1)) + Y(1:N)
        ELSE
          A=(Y(1:N) - Y(NP1)*YP(1:N))/(1.0 - Y(NP1))
        ENDIF
C IF NEW AND OLD  A  DIFFER BY TOO MUCH, TRACKING SHOULD NOT CONTINUE.
        IF (ANY(ABS(A-AOLD) .GT. 1.0+ABS(AOLD))) THEN
            ARCLEN=ARCLEN+S
            IFLAG=5
            CALL CLEANUPALL
            RETURN
        ENDIF
      ENDIF
      GO TO 100
80    IF (Y(NP1) .LE. .99  .OR. ST99) GO TO 100
C WHEN LAMBDA REACHES .99, THE PROBLEM WILL BE RESTARTED WITH
C A NEW  A  VECTOR.
90    ST99=.TRUE.
      EPSSTP=EPS
      ARCTOL=EPS
      GO TO 60
C
C SET DIFFERENT ERROR TOLERANCE FOR HIGH CURVATURE COMPONENTS OF THE
C TRAJECTORY Y(S).
100   CURTOL=CURSW*HOLD
      EPST=EPS/EPSSTP
      WHERE (ABS(YP-YPOLD) .LE. CURTOL)
        WT=(ABS(Y)+1.0)
      ELSEWHERE
        WT=(ABS(Y)+1.0)*EPST
      END WHERE
C
C TAKE A STEP ALONG THE CURVE.
      CALL STEPDS(FODEDS,NP1,Y,S,H,EPSSTP,WT,START,HOLD,K,KOLD,CRASH,
     &     PHI,P,YP,ALPHAS,W,G,KSTEPS,XOLD,IVC,IV,KGI,GI,
     &     IFLAGC,YPOLD,A,NDIMA,LENQR,MODE,NFEC)
C PRINT LATEST POINT ON CURVE IF REQUESTED.
      IF (TRACE .GT. 0) THEN
        WRITE (TRACE,117) ITER,NFEC,S,Y(NP1),(Y(JW),JW=1,N)
117     FORMAT(/' STEP',I5,3X,'NFE =',I5,3X,'ARC LENGTH =',F9.4,3X,
     &  'LAMBDA =',F7.4,5X,'X vector:'/(1X,6ES12.4))
      ENDIF
      NFE=NFEC
C CHECK IF THE STEP WAS SUCCESSFUL.
      IF (IFLAGC .EQ. 4) THEN
        ARCLEN=ARCLEN+S
        IFLAG=4
        CALL CLEANUPALL
        RETURN
      ENDIF
      IF (CRASH) THEN
C RETURN CODE FOR ERROR TOLERANCE TOO SMALL.
        IFLAG=2
C CHANGE ERROR TOLERANCES.
        EPS=EPSSTP
        IF (ARCTOL .LT. EPS) ARCTOL=EPS
C CHANGE LIMIT ON NUMBER OF ITERATIONS.
        LIMIT=LIMIT-ITER
        CALL CLEANUP
        RETURN
      ENDIF
C
      IF (Y(NP1) .GE. 1.0) THEN
        IF (ST99) GO TO 160
C
C IF LAMBDA .GE. 1.0 BUT THE PROBLEM HAS NOT BEEN RESTARTED
C WITH A NEW  A  VECTOR, BACK UP AND RESTART.
C
        S99=S-.5*HOLD
C GET AN APPROXIMATE ZERO Y(S) WITH  Y(NP1)=LAMBDA .LT. 1.0  .
135     CALL SINTRP(S,Y,S99,WT,YP,NP1,KOLD,PHI,IVC,IV,KGI,GI,
     &     ALPHAS,G,W,XOLD,P)
        IF (WT(NP1) .LT. 1.0) GO TO 140
        S99=.5*(S-HOLD+S99)
        GO TO 135
C
140     Y=WT
        YPOLD=YP
        S=S99
        GO TO 90
      ENDIF
C
      END DO MAIN_LOOP ! *****  END OF MAIN LOOP.  *****
C
C LAMBDA HAS NOT REACHED 1 IN 1000 STEPS.
      IFLAG=3
      CALL CLEANUP
      RETURN
C
C USE INVERSE INTERPOLATION TO GET THE ANSWER AT LAMBDA = 1.0 .
C
160   SA=S-HOLD
      SB=S
      LCODE=1
170   CALL ROOT(SOUT,Y1SOUT,SA,SB,EPS,EPS,LCODE)
C ROOT  FINDS S SUCH THAT Y(NP1)(S) = LAMBDA = 1 .
      IF (LCODE .GT. 0) GO TO 190
      CALL SINTRP(S,Y,SOUT,WT,YP,NP1,KOLD,PHI,IVC,IV,KGI,GI,
     &     ALPHAS,G,W,XOLD,P)
      Y1SOUT=WT(NP1)-1.0
      GO TO 170
190   IFLAG=1
C SET IFLAG = 6 IF  ROOT  COULD NOT GET  LAMBDA = 1.0  .
      IF (LCODE .GT. 2) IFLAG=6
      ARCLEN=ARCLEN+SA
C LAMBDA(SA) = 1.0 .
      CALL SINTRP(S,Y,SA,WT,YP,NP1,KOLD,PHI,IVC,IV,KGI,GI,
     &     ALPHAS,G,W,XOLD,P)
C
      Y=WT
      CALL CLEANUPALL
      RETURN
C
      CONTAINS
      SUBROUTINE CLEANUPALL
      IF (ALLOCATED(AOLD)) DEALLOCATE(AOLD)
      IF (ALLOCATED(P)) DEALLOCATE(P)
      IF (ALLOCATED(PHI)) DEALLOCATE(PHI)
      IF (ALLOCATED(WT)) DEALLOCATE(WT)
      IF (ALLOCATED(YP)) DEALLOCATE(YP)
      IF (ALLOCATED(YPOLD)) DEALLOCATE(YPOLD)
      CALL CLEANUP
      RETURN
      END SUBROUTINE CLEANUPALL
      SUBROUTINE CLEANUP
      IF (ALLOCATED(QR)) DEALLOCATE(QR)
      IF (ALLOCATED(ROWPOS)) DEALLOCATE(ROWPOS)
      IF (ALLOCATED(COLPOS)) DEALLOCATE(COLPOS)
      IF (ALLOCATED(PP)) DEALLOCATE(PP)
      IF (ALLOCATED(PAR)) DEALLOCATE(PAR)
      IF (ALLOCATED(IPAR)) DEALLOCATE(IPAR)
      RETURN
      END SUBROUTINE CLEANUP
      END SUBROUTINE FIXPDS
!
      SUBROUTINE FIXPNF(N,Y,IFLAG,ARCRE,ARCAE,ANSRE,ANSAE,TRACE,A,
     &   SSPAR,NFE,ARCLEN,  POLY_SWITCH)
C
C Subroutine  FIXPNF  finds a fixed point or zero of the
C N-dimensional vector function F(X), or tracks a zero curve
C of a general homotopy map RHO(A,LAMBDA,X).  For the fixed 
C point problem F(X) is assumed to be a C2 map of some ball 
C into itself.  The equation  X = F(X)  is solved by
C following the zero curve of the homotopy map
C
C  LAMBDA*(X - F(X)) + (1 - LAMBDA)*(X - A)  ,
C
C starting from LAMBDA = 0, X = A.  The curve is parameterized
C by arc length S, and is followed by solving the ordinary
C differential equation  D(HOMOTOPY MAP)/DS = 0  for
C Y(S) = (LAMBDA(S), X(S)) using a Hermite cubic predictor and a
C corrector which returns to the zero curve along the flow normal
C to the Davidenko flow (which consists of the integral curves of
C D(HOMOTOPY MAP)/DS ).
C
C For the zero finding problem F(X) is assumed to be a C2 map
C such that for some R > 0,  X*F(X) >= 0  whenever NORM(X) = R.
C The equation  F(X) = 0  is solved by following the zero curve
C of the homotopy map
C
C   LAMBDA*F(X) + (1 - LAMBDA)*(X - A)
C
C emanating from LAMBDA = 0, X = A.
C
C  A  must be an interior point of the above mentioned balls.
C
C For the curve tracking problem RHO(A,LAMBDA,X) is assumed to
C be a C2 map from E**M X [0,1) X E**N into E**N, which for
C almost all parameter vectors A in some nonempty open subset
C of E**M satisfies
C
C  rank [D RHO(A,LAMBDA,X)/D LAMBDA , D RHO(A,LAMBDA,X)/DX] = N
C
C for all points (LAMBDA,X) such that RHO(A,LAMBDA,X)=0.  It is
C further assumed that
C
C           rank [ D RHO(A,0,X0)/DX ] = N  .
C
C With A fixed, the zero curve of RHO(A,LAMBDA,X) emanating
C from  LAMBDA = 0, X = X0  is tracked until  LAMBDA = 1  by
C solving the ordinary differential equation
C D RHO(A,LAMBDA(S),X(S))/DS = 0  for  Y(S) = (LAMBDA(S), X(S)),
C where S is arc length along the zero curve.  Also the homotopy
C map RHO(A,LAMBDA,X) is assumed to be constructed such that
C
C              D LAMBDA(0)/DS > 0  .
C
C
C For the fixed point and zero finding problems, the user must supply 
C a subroutine  F(X,V)  which evaluates F(X) at X and returns the 
C vector F(X) in V, and a subroutine  FJAC(X,V,K)  which returns in V 
C the Kth column of the Jacobian matrix of F(X) evaluated at X.  For 
C the curve tracking problem, the user must supply a subroutine  
C  RHO(A,LAMBDA,X,V)  which evaluates the homotopy map RHO at 
C (A,LAMBDA,X) and returns the vector RHO(A,LAMBDA,X) in V, and a
C subroutine  RHOJAC(A,LAMBDA,X,V,K)  which returns in V the Kth
C column of the N X (N+1) Jacobian matrix [D RHO/D LAMBDA, D RHO/DX]
C evaluated at (A,LAMBDA,X).  FIXPNF  directly or indirectly uses
C the subroutines  F (or  RHO ),  FJAC (or  RHOJAC ), 
C   ROOT,  ROOTNF,  STEPNF,  the LAPACK routines  DGEQPF,  DORMQR,  
C their auxiliary routines, and the BLAS routines  DCOPY,
C   DDOT,  DGEMM,  DGEMV,  DGER,  DNRM2,  DSCAL,  DSWAP,  DTRMM,  DTRMV, 
C   IDAMAX.  The module  REAL_PRECISION  specifies 64-bit
C real arithmetic, which the user may want to change.
C
C
C ON INPUT:
C
C N  is the dimension of X, F(X), and RHO(A,LAMBDA,X).
C
C Y(:)  is an array of length  N + 1.  (Y(2),...,Y(N+1)) = A  is the
C    starting point for the zero curve for the fixed point and 
C    zero finding problems.  (Y(2),...,Y(N+1)) = X0  for the curve
C    tracking problem.
C
C IFLAG  can be -2, -1, 0, 2, or 3.  IFLAG  should be 0 on the 
C    first call to  FIXPNF  for the problem  X=F(X), -1 for the
C    problem  F(X)=0, and -2 for the problem  RHO(A,LAMBDA,X)=0.
C    In certain situations  IFLAG  is set to 2 or 3 by  FIXPNF,
C    and  FIXPNF  can be called again without changing  IFLAG.
C
C ARCRE , ARCAE  are the relative and absolute errors, respectively,
C    allowed the normal flow iteration along the zero curve.  If
C    ARC?E .LE. 0.0  on input it is reset to  .5*SQRT(ANS?E) .
C    Normally  ARC?E should be considerably larger than  ANS?E .
C
C ANSRE , ANSAE  are the relative and absolute error values used for
C    the answer at LAMBDA = 1.  The accepted answer  Y = (LAMBDA, X)
C    satisfies
C
C       |Y(1) - 1|  .LE.  ANSRE + ANSAE           .AND.
C
C       ||Z||  .LE.  ANSRE*||X|| + ANSAE          where
C
C    (.,Z) is the Newton step to Y.
C
C TRACE  is an integer specifying the logical I/O unit for
C    intermediate output.  If  TRACE .GT. 0  the points computed on
C    the zero curve are written to I/O unit  TRACE .
C
C A(:)  contains the parameter vector  A .  For the fixed point
C    and zero finding problems, A  need not be initialized by the
C    user, and is assumed to have length  N.  For the curve
C    tracking problem, A  must be initialized by the user.
C
C SSPAR(1:8) = (LIDEAL, RIDEAL, DIDEAL, HMIN, HMAX, BMIN, BMAX, P)  is
C    a vector of parameters used for the optimal step size estimation.
C    If  SSPAR(J) .LE. 0.0  on input, it is reset to a default value
C    by  FIXPNF .  Otherwise the input value of  SSPAR(J)  is used.
C    See the comments below and in  STEPNF  for more information about
C    these constants.
C
C POLY_SWITCH  is an optional logical variable used only by the driver
C    POLSYS1H  for polynomial systems.
C
C
C ON OUTPUT:
C
C N , TRACE , A  are unchanged.
C
C Y(1) = LAMBDA, (Y(2),...,Y(N+1)) = X, and Y is an approximate
C    zero of the homotopy map.  Normally LAMBDA = 1 and X is a
C    fixed point(zero) of F(X).  In abnormal situations LAMBDA
C    may only be near 1 and X is near a fixed point(zero).
C
C IFLAG =
C  -2   causes  FIXPNF  to initialize everything for the problem
C       RHO(A,LAMBDA,X) = 0 (use on first call).
C
C  -1   causes  FIXPNF  to initialize everything for the problem
C       F(X) = 0 (use on first call).
C
C   0   causes  FIXPNF  to initialize everything for the problem
C       X = F(X) (use on first call).
C
C   1   Normal return.
C
C   2   Specified error tolerance cannot be met.  Some or all of
C       ARCRE , ARCAE , ANSRE , ANSAE  have been increased to 
C       suitable values.  To continue, just call  FIXPNF  again 
C       without changing any parameters.
C
C   3   STEPNF  has been called 1000 times.  To continue, call
C       FIXPNF  again without changing any parameters.
C
C   4   Jacobian matrix does not have full rank.  The algorithm
C       has failed (the zero curve of the homotopy map cannot be
C       followed any further).
C
C   5   The tracking algorithm has lost the zero curve of the
C       homotopy map and is not making progress.  The error tolerances
C       ARC?E  and  ANS?E  were too lenient.  The problem should be
C       restarted by calling  FIXPNF  with smaller error tolerances
C       and  IFLAG = 0 (-1, -2).
C
C   6   The normal flow Newton iteration in  STEPNF  or  ROOTNF
C       failed to converge.  The error tolerances  ANS?E  may be too
C       stringent.
C
C   7   Illegal input parameters, a fatal error.
C
C   8   Memory allocation error, fatal.
C
C ARCRE , ARCAE , ANSRE , ANSAE  are unchanged after a normal return 
C    (IFLAG = 1).  They are increased to appropriate values on the 
C    return  IFLAG = 2 .
C
C NFE  is the number of function evaluations (= number of
C    Jacobian matrix evaluations).
C
C ARCLEN  is the length of the path followed.
C
C
C Allocatable and automatic work arrays:
C
C YP(1:N+1)  is a work array containing the tangent vector to 
C    the zero curve at the current point  Y .
C
C YOLD(1:N+1)  is a work array containing the previous point found
C    on the zero curve.
C
C YPOLD(1:N+1)  is a work array containing the tangent vector to 
C    the zero curve at  YOLD .
C
C QR(1:N,1:N+2), ALPHA(1:3*N+3), TZ(1:N+1), PIVOT(1:N+1), W(1:N+1),
C    WP(1:N+1), Z0(1:N+1), Z1(1:N+1)  are all work arrays used by
C    STEPNF  to calculate the tangent vectors and Newton steps.
C
C
      USE REAL_PRECISION
      INTEGER, INTENT(IN)::N,TRACE
      REAL (KIND=R8), DIMENSION(:), INTENT(IN OUT)::A,Y
      INTEGER, INTENT(IN OUT)::IFLAG
      REAL (KIND=R8), INTENT(IN OUT)::ANSAE,ANSRE,ARCAE,ARCRE,
     &    SSPAR(8)
      INTEGER, INTENT(OUT)::NFE
      REAL (KIND=R8), INTENT(OUT)::ARCLEN
      LOGICAL, INTENT(IN), OPTIONAL::POLY_SWITCH
C
C LOCAL VARIABLES.
      REAL (KIND=R8), SAVE:: ABSERR,CURTOL,H,HOLD,RELERR,S
      INTEGER, SAVE:: IFLAGC,ITER,JW,LIMIT,NC,NFEC,NP1
      LOGICAL, SAVE:: CRASH,POLSYS,START
      INTERFACE
        FUNCTION DNRM2(N,X,STRIDE)
        USE REAL_PRECISION
        INTEGER:: N,STRIDE
        REAL (KIND=R8):: DNRM2,X(N)
        END FUNCTION DNRM2
      END INTERFACE
C
C ALLOCATABLE AND AUTOMATIC ARRAYS.
      REAL (KIND=R8), DIMENSION(:), ALLOCATABLE, SAVE:: YOLD,YP,YPOLD
      REAL (KIND=R8):: ALPHA(3*N+3),QR(N,N+2),TZ(N+1),
     &  W(N+1),WP(N+1),Z0(N+1),Z1(N+1)
      INTEGER:: PIVOT(N+1)
C
C ***** END OF DIMENSIONAL INFORMATION. *****
C
C LIMITD  IS AN UPPER BOUND ON THE NUMBER OF STEPS.  IT MAY BE
C CHANGED BY CHANGING THE FOLLOWING PARAMETER STATEMENT:
      INTEGER, PARAMETER:: LIMITD=1000
C
C SWITCH FROM THE TOLERANCE  ARC?E  TO THE (FINER) TOLERANCE  ANS?E  IF
C THE CURVATURE OF ANY COMPONENT OF  Y  EXCEEDS  CURSW.
      REAL (KIND=R8), PARAMETER:: CURSW=10.0
C
      INTERFACE
        SUBROUTINE STEPNF(N,NFE,IFLAG,START,CRASH,HOLD,H,RELERR,
     &    ABSERR,S,Y,YP,YOLD,YPOLD,A,QR,ALPHA,TZ,PIVOT,W,WP,
     &    Z0,Z1,SSPAR)
        USE REAL_PRECISION
        REAL (KIND=R8):: ABSERR,H,HOLD,RELERR,S
        INTEGER:: IFLAG,N,NFE
        LOGICAL:: CRASH,START
        REAL (KIND=R8):: A(:),ALPHA(3*N+3),QR(N,N+2),SSPAR(8),TZ(N+1),
     &    W(N+1),WP(N+1),Y(:),YOLD(N+1),YP(N+1),YPOLD(N+1),
     &    Z0(N+1),Z1(N+1)
        INTEGER:: PIVOT(N+1)
        END SUBROUTINE STEPNF
        SUBROUTINE ROOTNF(N,NFE,IFLAG,RELERR,ABSERR,Y,YP,YOLD,
     &    YPOLD,A,QR,ALPHA,TZ,PIVOT,W,WP)
        USE REAL_PRECISION
        REAL (KIND=R8):: ABSERR,RELERR
        INTEGER:: IFLAG,N,NFE
        REAL (KIND=R8):: A(:),ALPHA(3*N+3),QR(N,N+2),TZ(N+1),W(N+1),
     &    WP(N+1),Y(:),YOLD(N+1),YP(N+1),YPOLD(N+1)
        INTEGER:: PIVOT(N+1)
        END SUBROUTINE ROOTNF
      END INTERFACE
C
C
C :  :  :  :  :  :  :  :  :  :  :  :  :  :  :  :  :  :  :  :  :  :  :  :
C TEST LOGICAL SWITCH TO REFLECT INTENDED USAGE OF FIXPNF.
      IF (PRESENT(POLY_SWITCH)) THEN
        POLSYS=.TRUE.
      ELSE
        POLSYS=.FALSE.
      ENDIF
C
      IF (N .LE. 0  .OR.  ANSRE .LE. 0.0  .OR.  ANSAE .LT. 0.0
     &  .OR.  N+1 .NE. SIZE(Y)  .OR.
     &  ((IFLAG .EQ. -1  .OR.  IFLAG .EQ. 0) .AND.  N .NE. SIZE(A)))
     &  IFLAG=7
      IF (IFLAG .GE. -2  .AND.  IFLAG .LE. 0) GO TO 20
      IF (IFLAG .EQ. 2) GO TO 120
      IF (IFLAG .EQ. 3) GO TO 90
C ONLY VALID INPUT FOR  IFLAG  IS -2, -1, 0, 2, 3.
      IFLAG=7
      RETURN
C
C *****  INITIALIZATION BLOCK.  *****
C
20    ARCLEN=0.0
      IF (ARCRE .LE. 0.0) ARCRE=.5*SQRT(ANSRE)
      IF (ARCAE .LE. 0.0) ARCAE=.5*SQRT(ANSAE)
      NC=N
      NFEC=0
      IFLAGC=IFLAG
      NP1=N+1
      IF (ALLOCATED(YP)) DEALLOCATE(YP)
      IF (ALLOCATED(YOLD)) DEALLOCATE(YOLD)
      IF (ALLOCATED(YPOLD)) DEALLOCATE(YPOLD)
      ALLOCATE(YP(NP1),YOLD(NP1),YPOLD(NP1),STAT=JW)
      IF (JW /= 0) THEN
        IFLAG=8
        RETURN
      END IF
C SET INITIAL CONDITIONS FOR FIRST CALL TO  STEPNF .
      START=.TRUE.
      CRASH=.FALSE.
      HOLD=1.0
      H=.1
      S=0.0
      YPOLD(1)=1.0
      YP(1)=1.0
      Y(1)=0.0
      YPOLD(2:NP1)=0.0
      YP(2:NP1)=0.0
C SET OPTIMAL STEP SIZE ESTIMATION PARAMETERS.
C LET Z[K] DENOTE THE NEWTON ITERATES ALONG THE FLOW NORMAL TO THE
C DAVIDENKO FLOW AND Y THEIR LIMIT.
C IDEAL CONTRACTION FACTOR:  ||Z[2] - Z[1]|| / ||Z[1] - Z[0]||
      IF (SSPAR(1) .LE. 0.0) SSPAR(1)= .5
C IDEAL RESIDUAL FACTOR:  ||RHO(A, Z[1])|| / ||RHO(A, Z[0])||
      IF (SSPAR(2) .LE. 0.0) SSPAR(2)= .01
C IDEAL DISTANCE FACTOR:  ||Z[1] - Y|| / ||Z[0] - Y||
      IF (SSPAR(3) .LE. 0.0) SSPAR(3)= .5
C MINIMUM STEP SIZE  HMIN .
      IF (SSPAR(4) .LE. 0.0) SSPAR(4)=(SQRT(N+1.0)+4.0)*EPSILON(1.0_R8)
C MAXIMUM STEP SIZE  HMAX .
      IF (SSPAR(5) .LE. 0.0) SSPAR(5)= 1.0
C MINIMUM STEP SIZE REDUCTION FACTOR  BMIN .
      IF (SSPAR(6) .LE. 0.0) SSPAR(6)= .1_R8
C MAXIMUM STEP SIZE EXPANSION FACTOR  BMAX .
      IF (SSPAR(7) .LE. 0.0) SSPAR(7)= 3.0
C ASSUMED OPERATING ORDER  P .
      IF (SSPAR(8) .LE. 0.0) SSPAR(8)= 2.0
C
C LOAD  A  FOR THE FIXED POINT AND ZERO FINDING PROBLEMS.
      IF (IFLAGC .GE. -1) THEN
        A=Y(2:NP1)
      ENDIF
90    LIMIT=LIMITD
C
C *****  END OF INITIALIZATION BLOCK.  *****
C
120   MAIN_LOOP: DO ITER=1,LIMIT  ! *****  MAIN LOOP.  *****
      IF (Y(1) .LT. 0.0) THEN
        ARCLEN=S
        IFLAG=5
        CALL CLEANUP ; RETURN
      ENDIF
C
C SET DIFFERENT ERROR TOLERANCE IF THE TRAJECTORY Y(S) HAS ANY HIGH 
C CURVATURE COMPONENTS.
      CURTOL=CURSW*HOLD
      RELERR=ARCRE
      ABSERR=ARCAE
      IF (ANY(ABS(YP-YPOLD) .GT. CURTOL)) THEN
        RELERR=ANSRE
        ABSERR=ANSAE
      ENDIF
C
C TAKE A STEP ALONG THE CURVE.
      CALL STEPNF(NC,NFEC,IFLAGC,START,CRASH,HOLD,H,RELERR,ABSERR,
     &     S,Y,YP,YOLD,YPOLD,A,QR,ALPHA,TZ,PIVOT,W,WP,Z0,Z1,SSPAR)
C PRINT LATEST POINT ON CURVE IF REQUESTED.
      IF (TRACE .GT. 0) THEN
        WRITE (TRACE,217) ITER,NFEC,S,Y(1),(Y(JW),JW=2,NP1)
217     FORMAT(/' STEP',I5,3X,'NFE =',I5,3X,'ARC LENGTH =',F9.4,3X,
     &  'LAMBDA =',F7.4,5X,'X VECTOR:'/(1X,6ES12.4))
      ENDIF
      NFE=NFEC
C CHECK IF THE STEP WAS SUCCESSFUL.
      IF (IFLAGC .GT. 0) THEN
        ARCLEN=S
        IFLAG=IFLAGC
        CALL CLEANUP ; RETURN
      ENDIF
      IF (CRASH) THEN
C RETURN CODE FOR ERROR TOLERANCE TOO SMALL.
        IFLAG=2
C CHANGE ERROR TOLERANCES.
        IF (ARCRE .LT. RELERR) ARCRE=RELERR
        IF (ANSRE .LT. RELERR) ANSRE=RELERR
        IF (ARCAE .LT. ABSERR) ARCAE=ABSERR
        IF (ANSAE .LT. ABSERR) ANSAE=ABSERR
C CHANGE LIMIT ON NUMBER OF ITERATIONS.
        LIMIT=LIMIT-ITER
        RETURN
      ENDIF
C
      IF (Y(1) .GE. 1.0) THEN
C
C USE HERMITE CUBIC INTERPOLATION AND NEWTON ITERATION TO GET THE 
C ANSWER AT LAMBDA = 1.0 .
C
C SAVE  YOLD  FOR ARC LENGTH CALCULATION LATER.
        Z0=YOLD
        CALL ROOTNF(NC,NFEC,IFLAGC,ANSRE,ANSAE,Y,YP,YOLD,YPOLD,
     &              A,QR,ALPHA,TZ,PIVOT,W,WP)
C
        NFE=NFEC
        IFLAG=1
C SET ERROR FLAG IF  ROOTNF  COULD NOT GET THE POINT ON THE ZERO
C CURVE AT  LAMBDA = 1.0  .
        IF (IFLAGC .GT. 0) IFLAG=IFLAGC
C CALCULATE FINAL ARC LENGTH.
        W=Y-Z0
        ARCLEN=S - HOLD + DNRM2(NP1,W,1)
        CALL CLEANUP ; RETURN
      ENDIF
C
C FOR POLYNOMIAL SYSTEMS AND THE  POLSYS1H  HOMOTOPY MAP,
C D LAMBDA/DS .GE. 0 NECESSARILY.  THIS CONDITION IS FORCED HERE IF
C THE  POLY_SWITCH  VARIABLE IS PRESENT.
C
      IF (POLSYS) THEN
        IF (YP(1) .LT. 0.0) THEN
C REVERSE TANGENT DIRECTION SO D LAMBDA/DS = YP(1) > 0 .
          YP=-YP
          YPOLD=YP
C FORCE  STEPNF  TO USE THE LINEAR PREDICTOR FOR THE NEXT STEP ONLY.
          START=.TRUE.
        ENDIF
      ENDIF
C
      END DO MAIN_LOOP   ! *****  END OF MAIN LOOP.  *****
C
C LAMBDA HAS NOT REACHED 1 IN 1000 STEPS.
      IFLAG=3
      ARCLEN=S
      RETURN
C
      CONTAINS
        SUBROUTINE CLEANUP
        IF (ALLOCATED(YP)) DEALLOCATE(YP)
        IF (ALLOCATED(YOLD)) DEALLOCATE(YOLD)
        IF (ALLOCATED(YPOLD)) DEALLOCATE(YPOLD)
        END SUBROUTINE CLEANUP
      END SUBROUTINE FIXPNF
!
      SUBROUTINE FIXPNS(N,Y,IFLAG,ARCRE,ARCAE,ANSRE,ANSAE,TRACE,A,
     &   NFE,ARCLEN,MODE,LENQR,SSPAR)
C
C Subroutine  FIXPNS  finds a fixed point or zero of the
C N-dimensional vector function F(X), or tracks a zero curve
C of a general homotopy map RHO(A,X,LAMBDA).  For the fixed 
C point problem F(X) is assumed to be a C2 map of some ball 
C into itself.  The equation  X = F(X)  is solved by
C following the zero curve of the homotopy map
C
C  LAMBDA*(X - F(X)) + (1 - LAMBDA)*(X - A)  ,
C
C starting from LAMBDA = 0, X = A.  The curve is parameterized
C by arc length S, and is followed by solving the ordinary
C differential equation  D(HOMOTOPY MAP)/DS = 0  for
C Y(S) = (X(S), LAMBDA(S)) using a Hermite cubic predictor and a
C corrector which returns to the zero curve along the flow normal
C to the Davidenko flow (which consists of the integral curves of
C D(HOMOTOPY MAP)/DS ).
C
C For the zero finding problem F(X) is assumed to be a C2 map
C such that for some R > 0,  X*F(X) >= 0  whenever NORM(X) = R.
C The equation  F(X) = 0  is solved by following the zero curve
C of the homotopy map
C
C   LAMBDA*F(X) + (1 - LAMBDA)*(X - A)
C
C emanating from LAMBDA = 0, X = A.
C
C  A  must be an interior point of the above mentioned balls.
C
C For the curve tracking problem RHO(A,X,LAMBDA) is assumed to
C be a C2 map from E**M X E**N X [0,1) into E**N, which for
C almost all parameter vectors A in some nonempty open subset
C of E**M satisfies
C
C  rank [D RHO(A,X,LAMBDA)/DX , D RHO(A,X,LAMBDA)/D LAMBDA] = N
C
C for all points (X,LAMBDA) such that RHO(A,X,LAMBDA)=0.  It is
C further assumed that
C
C           rank [ D RHO(A,X0,0)/DX ] = N  .
C
C With A fixed, the zero curve of RHO(A,X,LAMBDA) emanating
C from  LAMBDA = 0, X = X0  is tracked until  LAMBDA = 1  by
C solving the ordinary differential equation
C D RHO(A,X(S),LAMBDA(S))/DS = 0  for  Y(S) = (X(S), LAMBDA(S)),
C where S is arc length along the zero curve.  Also the homotopy
C map RHO(A,X,LAMBDA) is assumed to be constructed such that
C
C              D LAMBDA(0)/DS > 0  .
C
C
C For the fixed point and zero finding problems, the user
C must supply a subroutine  F(X,V)  which evaluates F(X) at X
C and returns the vector F(X) in V, and a subroutine
C  FJACS(X)  which evaluates, if
C MODE = 1,
C   the (symmetric) Jacobian matrix of F(X) at X, and returns the
C   symmetric Jacobian matrix in packed skyline storage format in
C   QR, or if
C MODE = 2,
C   returns the (nonsymmetric) Jacobian matrix in sparse row format
C   in QR.  The MODE 1 format is defined by QR, LENQR, ROWPOS; the
C   MODE 2 format is defined by QR, LENQR, ROWPOS, COLPOS.
C
C For the curve tracking problem, the user must supply a subroutine
C  RHO(A,LAMBDA,X,V)  which evaluates the homotopy map RHO 
C at (A,X,LAMBDA) and returns the vector RHO(A,X,LAMBDA) in V, and 
C a subroutine  RHOJS(A,LAMBDA,X)  which, if
C MODE = 1,
C   returns in QR the symmetric N X N Jacobian matrix [D RHO/DX] 
C   evaluated at (A,X,LAMBDA) and stored in packed skyline format, 
C   and returns in PP the vector -(D RHO/D LAMBDA) evaluated at 
C   (A,X,LAMBDA).  This data structure is described by QR, LENQR,
C   ROWPOS, PP.  *** Note the minus sign in the definition of PP. ***  If
C MODE = 2,
C   the (nonsymmetric) N X (N+1) Jacobian matrix [D RHO/DX, D RHO/DLAMBDA]
C   evaluated at (A,X,LAMBDA) is returned in a data structure described
C   by QR, LENQR, ROWPOS, COLPOS.
C
C Whichever of the routines  F,  FJACS,  RHO,  RHOJS  are required
C should be supplied as external subroutines, conforming with the
C interfaces in the module  HOMOTOPY.
C
C
C Functions and subroutines directly or indirectly called by FIXPNS:
C F (or  RHO ), FJACS (or  RHOJS ), GMFADS , GMRES , GMRILUDS ,
C ILUFDS , ILUSOLVDS , MULTDS , MULT2DS , PCGDS , ROOT , ROOTNS ,
C SOLVDS , STEPNS , TANGNS , and the BLAS functions  DDOT , DLAIC1 ,
C DLAMCH , DNRM2 .  The module  REAL_PRECISION  specifies 64-bit
C real arithmetic, which the user may want to change.
C 
C
C ON INPUT:
C
C N  is the dimension of X, F(X), and RHO(A,X,LAMBDA).
C
C Y  is an array of length  N + 1.  (Y(1),...,Y(N)) = A  is the
C    starting point for the zero curve for the fixed point and 
C    zero finding problems.  (Y(1),...,Y(N)) = X0  for the curve
C    tracking problem.
C
C IFLAG  can be -2, -1, 0, 2, or 3.  IFLAG  should be 0 on the 
C    first call to  FIXPNS  for the problem  X=F(X), -1 for the
C    problem  F(X)=0, and -2 for the problem  RHO(A,X,LAMBDA)=0.
C    In certain situations  IFLAG  is set to 2 or 3 by  FIXPNS,
C    and  FIXPNS  can be called again without changing  IFLAG.
C
C ARCRE , ARCAE  are the relative and absolute errors, respectively,
C    allowed the normal flow iteration along the zero curve.  If
C    ARC?E .LE. 0.0  on input it is reset to  .5*SQRT(ANS?E) .
C    Normally  ARC?E should be considerably larger than  ANS?E .
C
C ANSRE , ANSAE  are the relative and absolute error values used for
C    the answer at LAMBDA = 1.  The accepted answer  Y = (X, LAMBDA)
C    satisfies
C
C       |Y(NP1) - 1|  .LE.  ANSRE + ANSAE           .AND.
C
C       ||Z||  .LE.  ANSRE*||X|| + ANSAE          where
C
C    (Z,.) is the Newton step to Y.
C
C TRACE  is an integer specifying the logical I/O unit for
C    intermediate output.  If  TRACE .GT. 0  the points computed on
C    the zero curve are written to I/O unit  TRACE .
C
C A(:)  contains the parameter vector  A.  For the fixed point
C    and zero finding problems, A  need not be initialized by the
C    user, and is assumed to have length  N.  For the curve
C    tracking problem, A  must be initialized by the user.
C
C MODE = 1 if the Jacobian matrix is symmetric and stored in a packed
C          skyline format;
C      = 2 if the Jacobian matrix is stored in a sparse row format.
C
C LENQR  is the number of nonzero entries in the sparse Jacobian
C    matrices, used to determine the sparse matrix data structures.
C
C SSPAR(1:8) = (LIDEAL, RIDEAL, DIDEAL, HMIN, HMAX, BMIN, BMAX, P)  is
C    a vector of parameters used for the optimal step size estimation.
C    If  SSPAR(J) .LE. 0.0  on input, it is reset to a default value
C    by  FIXPNS .  Otherwise the input value of  SSPAR(J)  is used.
C    See the comments below and in  STEPNS  for more information about
C    these constants.
C
C
C ON OUTPUT:
C
C N , TRACE , A  are unchanged.
C
C (Y(1),...,Y(N)) = X, Y(NP1) = LAMBDA, and Y is an approximate
C    zero of the homotopy map.  Normally LAMBDA = 1 and X is a
C    fixed point(zero) of F(X).  In abnormal situations LAMBDA
C    may only be near 1 and X is near a fixed point(zero).
C
C IFLAG =
C  -2   causes  FIXPNS  to initialize everything for the problem
C       RHO(A,X,LAMBDA) = 0 (use on first call).
C
C  -1   causes  FIXPNS  to initialize everything for the problem
C       F(X) = 0 (use on first call).
C
C   0   causes  FIXPNS  to initialize everything for the problem
C       X = F(X) (use on first call).
C
C   1   Normal return.
C
C   2   Specified error tolerance cannot be met.  Some or all of
C       ARCRE , ARCAE , ANSRE , ANSAE  have been increased to 
C       suitable values.  To continue, just call  FIXPNS  again 
C       without changing any parameters.
C
C   3   STEPNS  has been called 1000 times.  To continue, call
C       FIXPNS  again without changing any parameters.
C
C   4   The preconditioned conjugate gradient iteration failed to
C       converge, or the Jacobian matrix does not have full rank
C       or has a zero on the diagonal.  The algorithm has failed
C       (the zero curve of the homotopy map cannot be followed any
C       further).
C
C   5   The tracking algorithm has lost the zero curve of the
C       homotopy map and is not making progress.  The error tolerances
C       ARC?E  and  ANS?E  were too lenient.  The problem should be
C       restarted by calling  FIXPNS  with smaller error tolerances
C       and  IFLAG = 0 (-1, -2).
C
C   6   The normal flow Newton iteration in  STEPNS  or  ROOTNS
c       failed to converge.  The error tolerances  ANS?E  may be too
C       stringent.
C
C   7   Illegal input parameters, a fatal error.
C
C ARCRE , ARCAE , ANSRE , ANSAE  are unchanged after a normal return 
C    (IFLAG = 1).  They are increased to appropriate values on the 
C    return  IFLAG = 2 .
C
C NFE  is the number of function evaluations (= number of
C    Jacobian evaluations).
C
C ARCLEN  is the length of the path followed.
C
C
C Allocatable and automatic work arrays:
C
C YP(1:N+1)  is a work array containing the tangent vector to 
C    the zero curve at the current point  Y .
C
C YOLD(1:N+1)  is a work array containing the previous point found
C    on the zero curve.
C
C YPOLD(1:N+1)  is a work array containing the tangent vector to 
C    the zero curve at  YOLD .
C
C QR(1:LENQR), PP(1:N), ROWPOS(1:N+2), COLPOS(1:LENQR) are all work
C    arrays used to define the sparse Jacobian matrices, allocated
C    here, and distributed via the module  HOMOTOPY.
C
C
C
      USE HOMOTOPY, QR => QRSPARSE
      USE REAL_PRECISION
      INTEGER, INTENT(IN)::LENQR,MODE,N,TRACE
      REAL (KIND=R8), DIMENSION(:), INTENT(IN OUT)::A,Y
      INTEGER, INTENT(IN OUT)::IFLAG
      REAL (KIND=R8), INTENT(IN OUT)::ANSAE,ANSRE,ARCAE,ARCRE,SSPAR(8)
      INTEGER, INTENT(OUT)::NFE
      REAL (KIND=R8), INTENT(OUT)::ARCLEN
C
C *****  LOCAL VARIABLES.  *****
C
      REAL (KIND=R8), SAVE:: ABSERR,CURTOL,H,HOLD,RELERR,S
      INTEGER, SAVE:: IFLAGC,ITER,JW,LIMIT,NC,NFEC,NP1
      LOGICAL, SAVE:: CRASH,START
      INTERFACE
        FUNCTION DNRM2(N,X,STRIDE)
        USE REAL_PRECISION
        INTEGER:: N,STRIDE
        REAL (KIND=R8):: DNRM2,X(N)
        END FUNCTION DNRM2
      END INTERFACE
C ***** WORK ARRAYS. *****
      REAL (KIND=R8), ALLOCATABLE, DIMENSION(:), SAVE:: YP,YOLD,YPOLD
      REAL (KIND=R8):: TZ(N+1),W(N+1),WP(N+1),Z0(N+1),Z1(N+1)
C
C LIMITD  IS AN UPPER BOUND ON THE NUMBER OF STEPS.  IT MAY BE
C CHANGED BY CHANGING THE FOLLOWING PARAMETER STATEMENT:
      INTEGER, PARAMETER:: LIMITD=1000
C
C SWITCH FROM THE TOLERANCE  ARC?E  TO THE (FINER) TOLERANCE  ANS?E  IF
C THE CURVATURE OF ANY COMPONENT OF  Y  EXCEEDS  CURSW.
      REAL (KIND=R8), PARAMETER:: CURSW=10.0
C
      INTERFACE
        SUBROUTINE ROOTNS(NC,NFEC,IFLAGC,ANSRE,ANSAE,Y,YP,YOLD,YPOLD,
     &     A,MODE,LENQR)
        USE HOMOTOPY, QR => QRSPARSE
        USE REAL_PRECISION
        INTEGER, INTENT(IN):: LENQR,MODE,NC
        INTEGER, INTENT(IN OUT):: IFLAGC,NFEC
        REAL (KIND=R8), INTENT(IN):: A(:)
        REAL (KIND=R8), INTENT(IN):: ANSAE,ANSRE
        REAL (KIND=R8), DIMENSION(:), INTENT(IN OUT):: Y,YOLD,YP,YPOLD
        END SUBROUTINE ROOTNS
C
        SUBROUTINE STEPNS(NC,NFEC,IFLAGC,START,CRASH,HOLD,H,RELERR,
     &     ABSERR,S,Y,YP,YOLD,YPOLD,A,MODE,LENQR,SSPAR,TZ,W,WP,Z0,Z1)
        USE HOMOTOPY, QR => QRSPARSE
        USE REAL_PRECISION
        INTEGER, INTENT(IN):: LENQR,MODE,NC
        INTEGER, INTENT(IN OUT):: IFLAGC,NFEC
        LOGICAL, INTENT(IN OUT):: CRASH,START
        REAL (KIND=R8), INTENT(IN):: A(:),SSPAR(8)
        REAL (KIND=R8), INTENT(IN OUT):: ABSERR,H,HOLD,RELERR,S,
     &    Y(:),YOLD(:),YP(:),YPOLD(:)
        REAL (KIND=R8), INTENT(OUT), DIMENSION(:):: TZ,W,WP,Z0,Z1
        END SUBROUTINE STEPNS
      END INTERFACE
C
C ***** END OF SPECIFICATION INFORMATION. *****
C
C :  :  :  :  :  :  :  :  :  :  :  :  :  :  :  :  :  :  :  :  :  :  :  :
      IF (N .LE. 0  .OR.  ANSRE .LE. 0.0  .OR.  ANSAE .LT. 0.0
     &  .OR.  N+1 .NE. SIZE(Y) .OR.
     &  ((IFLAG .EQ. -1  .OR.  IFLAG .EQ. 0) .AND.  N .NE. SIZE(A))
     &  .OR.  MODE .LE. 0  .OR.  MODE .GE. 3)
     &                                                     IFLAG=7
      IF (IFLAG .GE. -2  .AND.  IFLAG .LE. 0) GO TO 20
      IF (IFLAG .EQ. 2) GO TO 120
      IF (IFLAG .EQ. 3) GO TO 90
C ONLY VALID INPUT FOR  IFLAG  IS -2, -1, 0, 2, 3.
      IFLAG=7
      RETURN
C
C *****  INITIALIZATION BLOCK.  *****
C
20    ARCLEN=0.0
      IF (ARCRE .LE. 0.0) ARCRE=.5*SQRT(ANSRE)
      IF (ARCAE .LE. 0.0) ARCAE=.5*SQRT(ANSAE)
      NC=N
      NFEC=0
      IFLAGC=IFLAG
      NP1=N+1
C ALLOCATE SAVED WORK ARRAYS.
      IF (ALLOCATED(YOLD)) DEALLOCATE(YOLD)
      IF (ALLOCATED(YP)) DEALLOCATE(YP)
      IF (ALLOCATED(YPOLD)) DEALLOCATE(YPOLD)
      ALLOCATE(YOLD(NP1),YP(NP1),YPOLD(NP1))
C SET INITIAL CONDITIONS FOR FIRST CALL TO  STEPNS .
      START=.TRUE.
      CRASH=.FALSE.
      HOLD=1.0
      H=.1
      S=0.0
      YPOLD(NP1)=1.0
      YP(NP1)=1.0
      Y(NP1)=0.0
      YPOLD(1:N)=0.0
      YP(1:N)=0.0
C SET OPTIMAL STEP SIZE ESTIMATION PARAMETERS.
C LET Z[K] DENOTE THE NEWTON ITERATES ALONG THE FLOW NORMAL TO THE
C DAVIDENKO FLOW AND Y THEIR LIMIT.
C IDEAL CONTRACTION FACTOR:  ||Z[2] - Z[1]|| / ||Z[1] - Z[0]||
      IF (SSPAR(1) .LE. 0.0) SSPAR(1)= .5_R8
C IDEAL RESIDUAL FACTOR:  ||RHO(A, Z[1])|| / ||RHO(A, Z[0])||
      IF (SSPAR(2) .LE. 0.0) SSPAR(2)= .01_R8
C IDEAL DISTANCE FACTOR:  ||Z[1] - Y|| / ||Z[0] - Y||
      IF (SSPAR(3) .LE. 0.0) SSPAR(3)= .5_R8
C MINIMUM STEP SIZE  HMIN .
      IF (SSPAR(4) .LE. 0.0) SSPAR(4)=(SQRT(N+1.0)+4.0)*EPSILON(1.0_R8)
C MAXIMUM STEP SIZE  HMAX .
      IF (SSPAR(5) .LE. 0.0) SSPAR(5)= 1.0
C MINIMUM STEP SIZE REDUCTION FACTOR  BMIN .
      IF (SSPAR(6) .LE. 0.0) SSPAR(6)= .1_R8
C MAXIMUM STEP SIZE EXPANSION FACTOR  BMAX .
      IF (SSPAR(7) .LE. 0.0) SSPAR(7)= 3.0
C ASSUMED OPERATING ORDER  P .
      IF (SSPAR(8) .LE. 0.0) SSPAR(8)= 2.0
C
C LOAD  A  FOR THE FIXED POINT AND ZERO FINDING PROBLEMS.
      IF (IFLAGC .GE. -1) A(1:N) = Y(1:N)
90    LIMIT=LIMITD
C ALLOCATE ARRAYS FOR SPARSE JACOBIAN MATRIX DATA STRUCTURE.
120   SELECT CASE (MODE)
        CASE (1)
          IF (.NOT. ALLOCATED(QR)) ALLOCATE(QR(LENQR))
          IF (.NOT. ALLOCATED(ROWPOS)) ALLOCATE(ROWPOS(N+2))
          IF (.NOT. ALLOCATED(PP)) ALLOCATE(PP(N))
        CASE (2)
          IF (.NOT. ALLOCATED(QR)) ALLOCATE(QR(LENQR))
          IF (.NOT. ALLOCATED(ROWPOS)) ALLOCATE(ROWPOS(N+2))
          IF (.NOT. ALLOCATED(COLPOS)) ALLOCATE(COLPOS(LENQR))
          IF ((.NOT. ALLOCATED(PP)) .AND. (IFLAGC .GE. -1))
     &      ALLOCATE(PP(N))
      END SELECT
C
C *****  END OF INITIALIZATION BLOCK.  *****
C
      MAIN_LOOP: DO ITER=1,LIMIT  ! *****  MAIN LOOP.  *****
      IF (Y(NP1) .LT. 0.0) THEN
        ARCLEN=S
        IFLAG=5
        CALL CLEANUPALL
        RETURN
      ENDIF
C
C SET DIFFERENT ERROR TOLERANCE IF THE TRAJECTORY Y(S) HAS ANY HIGH 
C CURVATURE COMPONENTS.
      CURTOL=CURSW*HOLD
      RELERR=ARCRE
      ABSERR=ARCAE
        IF (ANY(ABS(YP-YPOLD) .GT. CURTOL)) THEN
          RELERR=ANSRE
          ABSERR=ANSAE
        ENDIF
C
C TAKE A STEP ALONG THE CURVE.
      CALL STEPNS(NC,NFEC,IFLAGC,START,CRASH,HOLD,H,RELERR,
     &     ABSERR,S,Y,YP,YOLD,YPOLD,A,MODE,LENQR,SSPAR,TZ,W,WP,Z0,Z1)
C PRINT LATEST POINT ON CURVE IF REQUESTED.
      IF (TRACE .GT. 0) THEN
        WRITE (TRACE,217) ITER,NFEC,S,Y(NP1),(Y(JW),JW=1,NC)
217     FORMAT(/' STEP',I5,3X,'NFE =',I5,3X,'ARC LENGTH =',F9.4,3X,
     &  'LAMBDA =',F7.4,5X,'X vector:'/(1X,6ES12.4))
      ENDIF
      NFE=NFEC
C CHECK IF THE STEP WAS SUCCESSFUL.
      IF (IFLAGC .GT. 0) THEN
        ARCLEN=S
        IFLAG=IFLAGC
        CALL CLEANUPALL
        RETURN
      ENDIF
      IF (CRASH) THEN
C RETURN CODE FOR ERROR TOLERANCE TOO SMALL.
        IFLAG=2
C CHANGE ERROR TOLERANCES.
        IF (ARCRE .LT. RELERR) ARCRE=RELERR
        IF (ANSRE .LT. RELERR) ANSRE=RELERR
        IF (ARCAE .LT. ABSERR) ARCAE=ABSERR
        IF (ANSAE .LT. ABSERR) ANSAE=ABSERR
C CHANGE LIMIT ON NUMBER OF ITERATIONS.
        LIMIT=LIMIT-ITER
        CALL CLEANUP
        RETURN
      ENDIF
C
      IF (Y(NP1) .GE. 1.0) THEN
C
C USE HERMITE CUBIC INTERPOLATION AND NEWTON ITERATION TO GET THE 
C ANSWER AT LAMBDA = 1.0 .
C
C SAVE  YOLD  FOR ARC LENGTH CALCULATION LATER.
        W=YOLD
C
        CALL ROOTNS(NC,NFEC,IFLAGC,ANSRE,ANSAE,Y,YP,YOLD,YPOLD,
     &              A,MODE,LENQR)
C
        NFE=NFEC
        IFLAG=1
C SET ERROR FLAG IF  ROOTNS  COULD NOT GET THE POINT ON THE ZERO
C CURVE AT  LAMBDA = 1.0  .
        IF (IFLAGC .GT. 0) IFLAG=IFLAGC
C CALCULATE FINAL ARC LENGTH.
        W = Y - W
        ARCLEN = S - HOLD + DNRM2(NP1,W,1)
        CALL CLEANUPALL
        RETURN
      ENDIF
C
      END DO MAIN_LOOP  !  *****  END OF MAIN LOOP.  *****
C
C LAMBDA HAS NOT REACHED 1 IN 1000 STEPS.
      IFLAG=3
      ARCLEN=S
      CALL CLEANUP
      RETURN
C
      CONTAINS
      SUBROUTINE CLEANUPALL
      IF (ALLOCATED(YOLD)) DEALLOCATE(YOLD)
      IF (ALLOCATED(YP)) DEALLOCATE(YP)
      IF (ALLOCATED(YPOLD)) DEALLOCATE(YPOLD)
      CALL CLEANUP
      RETURN
      END SUBROUTINE CLEANUPALL
      SUBROUTINE CLEANUP
      IF (ALLOCATED(QR)) DEALLOCATE(QR)
      IF (ALLOCATED(ROWPOS)) DEALLOCATE(ROWPOS)
      IF (ALLOCATED(COLPOS)) DEALLOCATE(COLPOS)
      IF (ALLOCATED(PP)) DEALLOCATE(PP)
      IF (ALLOCATED(PAR)) DEALLOCATE(PAR)
      IF (ALLOCATED(IPAR)) DEALLOCATE(IPAR)
      RETURN
      END SUBROUTINE CLEANUP
      END SUBROUTINE FIXPNS
!
      SUBROUTINE FIXPQF(N,Y,IFLAG,ARCRE,ARCAE,ANSRE,ANSAE,TRACE,A,
     &     SSPAR,NFE,ARCLEN)
C
C Subroutine  FIXPQF  finds a fixed point or zero of the 
C N-dimensional vector function  F(X), or tracks a zero curve of a 
C general homotopy map  RHO(A,LAMBDA,X).  For the fixed point problem
C F(X) is assumed to be a C2 map of some ball into itself.  The 
C equation  X=F(X)  is solved by following the zero curve of the 
C homotopy map
C
C  LAMBDA*(X - F(X)) + (1 - LAMBDA)*(X - A) ,
C
C starting from  LAMBDA = 0, X = A.   The curve is parameterized
C by arc length  S, and is followed by solving the ordinary 
C differential equation  D(HOMOTOPY MAP)/DS = 0  for  
C Y(S) = (LAMBDA(S), X(S)).  This is done by using a Hermite cubic 
C predictor and a corrector which returns to the zero curve in a 
C hyperplane perpendicular to the tangent to the zero curve at the 
C most recent point.
C
C For the zero finding problem  F(X)  is assumed to be a C2 map
C such that for some  R > 0,  X*F(X) >= 0  whenever  NORM(X) = R.
C The equation  F(X) = 0  is solved by following the zero curve of
C the homotopy map
C
C  LAMBDA*F(X) + (1 - LAMBDA)*(X - A)
C
C emanating from  LAMBDA = 0, X = A.
C
C A  must be an interior point of the above mentioned balls.
C
C For the curve tracking problem RHO(A,LAMBDA,X) is assumed to 
C be a C2 map from  E**M X [0,1) X E**N  into  E**N, which for 
C almost all parameter vectors  A  in some nonempty open subset
C of E**M satisfies
C
C  rank [D RHO(A,LAMBDA,X)/D LAMBDA, D RHO(A,LAMBDA,X)/DX] = N
C
C for all points  (LAMBDA,X)  such that  RHO(A,LAMBDA,X) = 0.  It is
C further assumed that
C
C         rank [ D RHO(A,0,X0)/DX ] = N.
C
C With  A  fixed, the zero curve of  RHO(A,LAMBDA,X)  emanating from
C LAMBDA = 0, X = X0  is tracked until  LAMBDA = 1  by solving the 
C ordinary differential equation    D RHO(A,LAMBDA(S),X(S))/DS = 0
C for  Y(S) = (LAMBDA(S), X(S)), where  S  is arc length along the
C zero curve.  Also the homotopy map  RHO(A,LAMBDA,X)  is assumed to
C be constructed such that
C
C         D LAMBDA(0)/DS > 0.
C
C For the fixed point and zero finding problems, the user must supply
C a subroutine  F(X,V)  which evaluates  F(X)  at  X  and returns the
C vector F(X) in  V, and a subroutine  FJAC(X,V,K)  which returns in  V
C the Kth column of the Jacobian matrix of F(X) evaluated at X.  For
C the curve tracking problem, the user must supply a subroutine
C RHO(A,LAMBDA,X,V)  which evaluates the homotopy map  RHO at
C (A,LAMBDA,X)  and returns the vector  RHO(A,LAMBDA,X)  in  V, and
C a subroutine  RHOJAC(A,LAMBDA,X,V,K)  which returns in  V
C the Kth column of the  N X (N+1)  Jacobian matrix  
C [D RHO/D LAMBDA, D RHO/DX]  evaluated at  (A,LAMBDA,X).  FIXPQF
C directly or indirectly uses the subroutines  F (or RHO), 
C   FJAC (or RHOJAC),  ROOT,  ROOTQF,  STEPQF,  TANGQF,  UPQRQF,
C the LAPACK routines  DGEQRF,  DORGQR, their auxiliary routines,
C and the BLAS routines  DCOPY,  DDOT,  DGEMM,  DGEMV,
C   DGER,  DNRM2,  DSCAL,  DTPMV,  DTPSV,  DTRMM,  and   DTRMV.
C The module  REAL_PRECISION  specifies 64-bit real arithmetic,
C which the user may want to change.
C
C
C ON INPUT:
C
C N  is the dimension of X, F(X), and RHO(A,LAMBDA,X).
C
C Y(1:N+1)  contains the starting point for tracking the homotopy map.
C    (Y(2),...,Y(N+1)) = A  for the fixed point and zero finding 
C    problems.  (Y(2),...,Y(N+1)) = X0  for the curve tracking problem.
C    Y(1)  need not be defined by the user.
C
C IFLAG can be -2, -1, 0, 2, or 3.  IFLAG should be 0 on the first
C    call to  FIXPQF  for the problem  X=F(X), -1 for the problem
C    F(X)=0, and -2 for the problem  RHO(A,LAMBDA,X)=0.   In certain
C    situations  IFLAG  is set to 2 or 3 by  FIXPQF, and  FIXPQF  can
C    be called again without changing  IFLAG.
C
C ARCRE, ARCAE  are the relative and absolute errors, respectively,
C    allowed the quasi-Newton iteration along the zero curve.  If
C    ARC?E .LE. 0.0  on input, it is reset to  .5*SQRT(ANS?E).
C    Normally  ARC?E  should be considerably larger than  ANS?E.
C
C ANSRE, ANSAE  are the relative and absolute error values used for 
C    the answer at  LAMBDA = 1.  The accepted answer  Y = (LAMBDA, X)
C    satisfies
C
C      |Y(1) - 1| .LE. ANSRE + ANSAE      .AND.
C  
C      ||DZ|| .LE. ANSRE*||Y|| + ANSAE      where
C
C      DZ is the quasi-Newton step to Y.
C
C TRACE  is an integer specifying the logical I/O unit for
C    intermediate output.  If  TRACE .GT. 0  the points computed on
C    the zero curve are written to I/O unit  TRACE .
C
C A(:)  contains the parameter vector  A.  For the fixed point
C    and zero finding problems,  A  need not be initialized by the 
C    user, and is assumed to have length  N.  For the curve
C    tracking problem,  A  must be initialized by the user.
C
C SSPAR(1:4) =  (HMIN, HMAX, BMIN, BMAX)  is a vector of parameters 
C    used for the optimal step size estimation.  A default value
C    can be specified for any of these four parameters by setting it
C    .LE. 0.0  on input.  See the comments in  STEPQF  for more
C    information about these parameters.
C
C
C ON OUTPUT:
C
C N , TRACE , A  are unchanged.
C
C Y(1) = LAMBDA, (Y(2),...,Y(N+1)) = X, and  Y  is an approximate
C    zero of the homotopy map.  Normally  LAMBDA = 1  and  X  is a
C    fixed point or zero of  F(X).   In abnormal situations,  LAMBDA
C    may only be near 1 and  X  near a fixed point or zero.
C
C IFLAG =
C
C   1   Normal return.
C
C   2   Specified error tolerance cannot be met.  Some or all of
C       ARCRE, ARCAE, ANSRE, ANSAE  have been increased to 
C       suitable values.  To continue, just call  FIXPQF  again
C       without changing any parameters.
C
C   3   STEPQF  has been called 1000 times.  To continue, call
C       FIXPQF  again without changing any parameters.
C
C   4   Jacobian matrix does not have full rank.  The algorithm
C       has failed (the zero curve of the homotopy map cannot be
C       followed any further).
C
C   5   The tracking algorithm has lost the zero curve of the 
C       homotopy map and is not making progress.  The error 
C       tolerances  ARC?E  and  ANS?E  were too lenient.  The problem 
C       should be restrarted by calling  FIXPQF  with smaller error 
C       tolerances and  IFLAG = 0 (-1, -2).
C
C   6   The quasi-Newton iteration in  STEPQF  or  ROOTQF  failed to
C       converge.  The error tolerances  ANS?E  may be too stringent.
C
C   7   Illegal input parameters, a fatal error.
C
C   8   Memory allocation error, fatal.
C
C ARCRE, ARCAE, ANSRE, ANSAE  are unchanged after a normal return
C    (IFLAG = 1).  They are increased to appropriate values on the
C    return  IFLAG = 2.
C
C NFE  is the number of Jacobian evaluations.
C
C ARCLEN  is the approximate length of the zero curve.  
C
C
C Allocatable and automatic work arrays:
C
C YP(1:N+1)  is a work array containing the tangent vector to the
C    zero curve at the current point  Y.
C
C YOLD(1:N+1) is a work array containing the previous point found
C    on the zero curve.
C
C YPOLD(1:N+1) is a work array containing the tangent vector to
C    the zero curve at  YOLD.
C
C Q(1:N+1,1:N+1), R((N+1)*(N+2)/2), F0(1:N+1), F1(1:N+1), Z0(1:N+1),
C    DZ(1:N+1), W(1:N+1), T(1:N+1), YSAV(1:N+1)  are all work arrays 
C    used by  STEPQF, TANGQF and ROOTQF to calculate the tangent 
C    vectors and quasi-Newton steps.
C
C
C ***** DECLARATIONS *****
      USE REAL_PRECISION
C
C     FUNCTION DECLARATIONS 
C 
      INTERFACE
        FUNCTION DNRM2(N,X,STRIDE)
        USE REAL_PRECISION
        INTEGER:: N,STRIDE
        REAL (KIND=R8):: DNRM2,X(N)
        END FUNCTION DNRM2
      END INTERFACE
C
C     LOCAL VARIABLES 
C
      REAL (KIND=R8), SAVE:: ABSERR, H, HOLD, RELERR, S, WK 
      INTEGER, SAVE:: IFLAGC, ITER, JW, LIMIT, NP1
      LOGICAL, SAVE:: CRASH, START       
C
C     SCALAR ARGUMENTS 
C
      REAL (KIND=R8):: ARCRE, ARCAE, ANSRE, ANSAE, ARCLEN
      INTEGER:: N,IFLAG,TRACE,NFE
C
C     ARRAY DECLARATIONS 
C
      REAL (KIND=R8), DIMENSION(:), ALLOCATABLE, SAVE::
     &  R,YOLD,YP,YPOLD
      REAL (KIND=R8), DIMENSION(:,:), ALLOCATABLE, SAVE:: Q
      REAL (KIND=R8):: A(:), DZ(N+1), F0(N+1), F1(N+1),
     &    SSPAR(4), T(N+1), W(N+1), Y(:), YSAV(N+1), Z0(N+1)
C 
C ***** END OF DECLARATIONS *****
      INTERFACE
        SUBROUTINE STEPQF(N,NFE,IFLAG,START,CRASH,HOLD,H,
     &    WK,RELERR,ABSERR,S,Y,YP,YOLD,YPOLD,A,Q,R,
     &    F0,F1,Z0,DZ,W,T,SSPAR)
        USE REAL_PRECISION
        INTEGER:: N, NFE, IFLAG
        LOGICAL:: START, CRASH
        REAL (KIND=R8):: HOLD, H, WK, RELERR, ABSERR, S
        REAL (KIND=R8):: A(:), DZ(N+1), F0(N+1), F1(N+1), 
     &    Q(N+1,N+1), R((N+1)*(N+2)/2), SSPAR(4), T(N+1), W(N+1),
     &    Y(:), YOLD(N+1), YP(N+1), YPOLD(N+1), Z0(N+1)
        END SUBROUTINE STEPQF
        SUBROUTINE ROOTQF(N,NFE,IFLAG,RELERR,ABSERR,Y,YP,YOLD,
     &    YPOLD,A,Q,R,DZ,Z,W,T,F0,F1)
        USE REAL_PRECISION
        REAL (KIND=R8):: RELERR, ABSERR
        INTEGER:: N, NFE, IFLAG
        REAL (KIND=R8):: A(:), DZ(N+1), F0(N+1), F1(N+1), 
     &    Q(N+1,N+1), R((N+1)*(N+2)/2), T(N+1), W(N+1),
     &    Y(:), YOLD(N+1), YP(N+1), YPOLD(N+1), Z(N+1)
        END SUBROUTINE ROOTQF
      END INTERFACE
C
C LIMITD IS AN UPPER BOUND ON THE NUMBER OF STEPS.  IT MAY BE 
C CHANGED BY CHANGING THE FOLLOWING PARAMETER STATEMENT:
      INTEGER, PARAMETER:: LIMITD = 1000
C
C ***** FIRST EXECUTABLE STATEMENT *****
C
C CHECK IFLAG
C
      IF (N .LE. 0  .OR.  ANSRE .LE. 0.0  .OR.  ANSAE .LT. 0.0
     &  .OR.  N+1 .NE. SIZE(Y)  .OR.
     &  ((IFLAG .EQ. -1  .OR.  IFLAG .EQ. 0) .AND.  N .NE. SIZE(A)))
     &  IFLAG=7
      IF (IFLAG .GE. -2 .AND. IFLAG .LE. 0) GO TO 10
      IF (IFLAG .EQ. 2) GO TO 50
      IF (IFLAG .EQ. 3) GO TO 40
C
C ONLY VALID INPUT FOR IFLAG IS -2, -1, 0, 2, 3.
C
      IFLAG = 7
      RETURN
C
C ***** INITIALIZATION BLOCK  *****
C
 10   ARCLEN = 0.0
      IF (ARCRE .LE. 0.0) ARCRE = .5*SQRT(ANSRE)
      IF (ARCAE .LE. 0.0) ARCAE = .5*SQRT(ANSAE)
      NFE=0
      IFLAGC = IFLAG
      NP1=N+1
      IF (ALLOCATED(Q)) DEALLOCATE(Q)
      IF (ALLOCATED(R)) DEALLOCATE(R)
      IF (ALLOCATED(YOLD)) DEALLOCATE(YOLD)
      IF (ALLOCATED(YP)) DEALLOCATE(YP)
      IF (ALLOCATED(YPOLD)) DEALLOCATE(YPOLD)
      ALLOCATE(Q(NP1,NP1),R(NP1*(N+2)/2),YOLD(NP1),YP(NP1),YPOLD(NP1),
     &  STAT=JW)
      IF (JW /= 0) THEN
        IFLAG=8
        RETURN
      END IF
C 
C SET INITIAL CONDITIONS FOR FIRST CALL TO STEPQF.
C
      START=.TRUE.
      CRASH=.FALSE.
      RELERR = ARCRE
      ABSERR = ARCAE
      HOLD=1.0
      H=0.1
      S=0.0
      YPOLD(1) = 1.0
      Y(1) = 0.0
      YPOLD(2:NP1)=0.0
C
C SET OPTIMAL STEP SIZE ESTIMATION PARAMETERS.
C
C     MINIMUM STEP SIZE HMIN
      IF (SSPAR(1) .LE. 0.0) SSPAR(1)=(SQRT(N+1.0)+4.0)*EPSILON(1.0_R8)
C     MAXIMUM STEP SIZE HMAX
      IF (SSPAR(2) .LE. 0.0) SSPAR(2)= 1.0
C     MINIMUM STEP REDUCTION FACTOR BMIN
      IF (SSPAR(3) .LE. 0.0) SSPAR(3)= 0.1_R8
C     MAXIMUM STEP EXPANSION FACTOR BMAX
      IF (SSPAR(4) .LE. 0.0) SSPAR(4)= 7.0
C
C LOAD  A  FOR THE FIXED POINT AND ZERO FINDING PROBLEMS.
C
      IF (IFLAGC .GE. -1) THEN
        A=Y(2:NP1)
      ENDIF
C
40    LIMIT=LIMITD
C
C ***** END OF INITIALIZATION BLOCK. *****
C
50    DO ITER=1,LIMIT   ! ***** MAIN LOOP. *****
      IF (Y(1) .LT. 0.0) THEN
        ARCLEN = S
        IFLAG = 5
        CALL CLEANUP ; RETURN
      END IF
C
C TAKE A STEP ALONG THE CURVE.
C
      CALL STEPQF(N,NFE,IFLAGC,START,CRASH,HOLD,H,WK,
     &    RELERR,ABSERR,S,Y,YP,YOLD,YPOLD,A,Q,R,F0,F1,Z0,DZ,
     &    W,T,SSPAR) 
C
C PRINT LATEST POINT ON CURVE IF REQUESTED.
C
      IF (TRACE .GT. 0) THEN
         WRITE (TRACE,217) ITER,NFE,S,Y(1),(Y(JW),JW=2,NP1)
 217     FORMAT(/' STEP',I5,3X,'NFE =',I5,3X,'ARC LENGTH =',F9.4,3X,
     &   'LAMBDA =',F7.4,5X,'X VECTOR:'/(1X,6ES12.4))
      ENDIF
C
C CHECK IF THE STEP WAS SUCCESSFUL.
C
      IF (IFLAGC .GT. 0) THEN
        ARCLEN=S
        IFLAG=IFLAGC
        CALL CLEANUP ; RETURN
      END IF
C
      IF (CRASH) THEN
C
C         RETURN CODE FOR ERROR TOLERANCE TOO SMALL.
C      
        IFLAG=2
C
C         CHANGE ERROR TOLERANCES.
C
        IF (ARCRE .LT. RELERR) THEN
          ARCRE=RELERR
          ANSRE=RELERR
        END IF
        IF (ARCAE .LT. ABSERR) ARCAE=ABSERR
C
C         CHANGE LIMIT ON NUMBER OF ITERATIONS.
C
        LIMIT = LIMIT - ITER
        RETURN
      END IF
C
C IF LAMBDA >= 1.0, USE ROOTQF TO FIND SOLUTION.
C
      IF (Y(1) .GE. 1.0) GOTO 500
C
      END DO   ! ***** END OF MAIN LOOP *****
C
C DID NOT CONVERGE IN  LIMIT  ITERATIONS, SET  IFLAG  AND RETURN.
C
      ARCLEN = S
      IFLAG = 3
      RETURN
C
C ***** FINAL STEP -- FIND SOLUTION AT LAMBDA=1 *****
C
C SAVE  YOLD  FOR ARC LENGTH CALCULATION LATER.
C
 500  YSAV=YOLD
C
C FIND SOLUTION.
C
      CALL ROOTQF(N,NFE,IFLAGC,ANSRE,ANSAE,Y,YP,YOLD,
     &    YPOLD,A,Q,R,DZ,Z0,W,T,F0,F1)
C
C CHECK IF SOLUTION WAS FOUND AND SET  IFLAG  ACCORDINGLY.
C
      IFLAG=1
C
C     SET ERROR FLAG IF ROOTQF COULD NOT GET THE POINT ON THE ZERO
C     CURVE AT  LAMBDA = 1.0.
C
      IF (IFLAGC .GT. 0) IFLAG=IFLAGC
C
C CALCULATE FINAL ARC LENGTH.
C
      DZ = Y - YSAV
      ARCLEN = S - HOLD + DNRM2(NP1,DZ,1)
C
C ***** END OF FINAL STEP *****
C
      CALL CLEANUP ; RETURN
C
      CONTAINS
        SUBROUTINE CLEANUP
        IF (ALLOCATED(Q)) DEALLOCATE(Q)
        IF (ALLOCATED(R)) DEALLOCATE(R)
        IF (ALLOCATED(YOLD)) DEALLOCATE(YOLD)
        IF (ALLOCATED(YP)) DEALLOCATE(YP)
        IF (ALLOCATED(YPOLD)) DEALLOCATE(YPOLD)
        END SUBROUTINE CLEANUP
      END SUBROUTINE FIXPQF
!
      SUBROUTINE FIXPQS(N,Y,IFLAG,ARCRE,ARCAE,ANSRE,ANSAE,TRACE,A,
     &     NFE,ARCLEN,MODE,LENQR,SSPAR)
C
C Subroutine  FIXPQS  finds a fixed point or zero of the 
C N-dimensional vector function  F(X), or tracks a zero curve of a 
C general homotopy map  RHO(A,X,LAMBDA).  For the fixed point problem
C F(X) is assumed to be a C2 map of some ball into itself.  The 
C equation  X=F(X)  is solved by following the zero curve of the 
C homotopy map
C
C  LAMBDA*(X - F(X)) + (1 - LAMBDA)*(X - A),
C
C starting from  LAMBDA = 0, X = A.   The curve is parameterized
C by arc length  S, and is followed by solving the ordinary 
C differential equation  D(HOMOTOPY MAP)/DS = 0  for  
C Y(S) = (X(S),LAMBDA(S)).  This is done by using a Hermite cubic 
C predictor and a corrector which returns to the zero curve in a 
C hyperplane perpendicular to the tangent to the zero curve at the 
C most recent point.
C
C For the zero finding problem  F(X)  is assumed to be a C2 map such
C that for some  R > 0,  X*F(X) >= 0  whenever  NORM(X) = R.
C The equation  F(X) = 0  is solved by following the zero curve of
C the homotopy map
C
C  LAMBDA*F(X) + (1 - LAMBDA)*(X - A)
C
C emanating from  LAMBDA = 0, X = A.
C
C A  must be an interior point of the above mentioned balls.
C
C For the curve tracking problem RHO(A,X,LAMBDA) is assumed to 
C be a C2 map from  E**M X [0,1) X E**N  into  E**N, which for 
C almost all parameter vectors  A  in some nonempty open subset
C of E**M satisfies
C
C  rank [D RHO(A,X,LAMBDA)/D LAMBDA, D RHO(A,X,LAMBDA)/DX] = N
C
C for all points  (X,LAMBDA)  such that  RHO(A,X,LAMBDA) = 0.  It is
C further assumed that
C
C         rank [ D RHO(A,X0,0)/DX ] = N.
C
C With  A  fixed, the zero curve of  RHO(A,X,LAMBDA)  emanating from
C LAMBDA = 0, X = X0  is tracked until  LAMBDA = 1  by solving the 
C ordinary differential equation  D RHO(A,X(S),LAMBDA(S))/DS = 0
C for  Y(S) = (X(S),LAMBDA(S)), where  S  is arc length along the
C zero curve.  Also the homotopy map  RHO(A,X,LAMBDA)  is assumed to
C be constructed such that
C
C         D LAMBDA(0)/DS > 0.
C
C For the fixed point and zero finding problems, the user
C must supply a subroutine  F(X,V)  which evaluates F(X) at X
C and returns the vector F(X) in V, and a subroutine
C  FJACS(X)  which evaluates, if
C MODE = 1,
C   the (symmetric) Jacobian matrix of F(X) at X, and returns the
C   symmetric Jacobian matrix in packed skyline storage format in
C   QR, or if
C MODE = 2,
C   returns the (nonsymmetric) Jacobian matrix in sparse row format
C   in QR.  The MODE 1 format is defined by QR, LENQR, ROWPOS; the
C   MODE 2 format is defined by QR, LENQR, ROWPOS, COLPOS.
C
C For the curve tracking problem, the user must supply a subroutine
C  RHO(A,LAMBDA,X,V)  which evaluates the homotopy map RHO 
C at (A,X,LAMBDA) and returns the vector RHO(A,X,LAMBDA) in V, and 
C a subroutine  RHOJS(A,LAMBDA,X)  which, if
C MODE = 1,
C   returns in QR the symmetric N X N Jacobian matrix [D RHO/DX] 
C   evaluated at (A,X,LAMBDA) and stored in packed skyline format, 
C   and returns in PP the vector -(D RHO/D LAMBDA) evaluated at 
C   (A,X,LAMBDA).  This data structure is described by QR, LENQR,
C   ROWPOS, PP.  *** Note the minus sign in the definition of PP. ***  If
C MODE = 2,
C   the (nonsymmetric) N X (N+1) Jacobian matrix [D RHO/DX, D RHO/DLAMBDA]
C   evaluated at (A,X,LAMBDA) is returned in a data structure described
C   by QR, LENQR, ROWPOS, COLPOS.
C
C Whichever of the routines  F,  FJACS,  RHO,  RHOJS  are required
C should be supplied as external subroutines, conforming with the
C interfaces in the module  HOMOTOPY.
C
C
C FIXPQS directly or indirectly uses the subroutines  
C F (or  RHO ), FJACS (or  RHOJS ), GMFADS , GMRES , GMRILUDS ,
C ILUFDS , ILUSOLVDS , MULTDS , MULT2DS , PCGDS , ROOT , ROOTNS ,
C SOLVDS , STEPNS , TANGNS , and the BLAS functions  DDOT , DLAIC1 ,
C DLAMCH , DNRM2 .  The module  REAL_PRECISION  specifies 64-bit
C real arithmetic, which the user may want to change.
C 
C
C ON INPUT:
C
C N  is the dimension of X, F(X), and RHO(A,X,LAMBDA).
C
C Y(1:N+1)  contains the starting point for tracking the homotopy map.
C    (Y(1),...,Y(N)) = A  for the fixed point and zero finding 
C    problems.  (Y(1),...,Y(N)) = X0  for the curve tracking problem.
C    Y(N+1)  need not be defined by the user.
C
C IFLAG  can be -2, -1, 0, 2, or 3.  IFLAG  should be 0 on the first
C    call to  FIXPQS  for the problem  X=F(X), -1 for the problem
C    F(X)=0, and -2 for the problem  RHO(A,X,LAMBDA)=0.   In certain
C    situations  IFLAG  is set to 2 or 3 by  FIXPQS, and  FIXPQS  can
C    be called again without changing  IFLAG.
C
C ARCRE, ARCAE  are the relative and absolute errors, respectively,
C    allowed the iteration along the zero curve.  If
C    ARC?E .LE. 0.0  on input, it is reset to  .5*SQRT(ANS?E).
C    Normally  ARC?E  should be considerably larger than  ANS?E.
C
C ANSRE, ANSAE  are the relative and absolute error values used for 
C    the answer at  LAMBDA = 1.  The accepted answer  Y = (X,LAMBDA)
C    satisfies
C
C      |Y(1) - 1| .LE. ANSRE + ANSAE      .AND.
C  
C      ||DZ|| .LE. ANSRE*||Y|| + ANSAE      where
C
C      DZ is the Newton step to Y.
C
C TRACE  is an integer specifying the logical I/O unit for
C    intermediate output.  If  TRACE .GT. 0  the points computed on
C    the zero curve are written to I/O unit  TRACE .
C
C A(:)  contains the parameter vector  A.  For the fixed point
C    and zero finding problems,  A  need not be initialized by the 
C    user, and is assumed to have length  N.  For the curve
C    tracking problem,  A  must be initialized by the user.
C
C MODE = 1 if the Jacobian matrix is symmetric and stored in a packed
C          skyline format;
C      = 2 if the Jacobian matrix is stored in a sparse row format.
C
C LENQR  is the number of nonzero entries in the sparse Jacobian
C    matrices, used to determine the sparse matrix data structures.
C
C SSPAR(1:4) =  (HMIN, HMAX, BMIN, BMAX)  is a vector of parameters 
C    used for the optimal step size estimation.  A default value
C    can be specified for any of these four parameters by setting it
C    .LE. 0.0  on input.  See the comments in  STEPQS  for more
C    information about these parameters.
C
C
C ON OUTPUT:
C
C N , TRACE , A , LENQR  are unchanged.
C
C Y(N+1) = LAMBDA, (Y(1),...,Y(N)) = X, and  Y  is an approximate
C    zero of the homotopy map.  Normally  LAMBDA = 1  and  X  is a
C    fixed point or zero of  F(X).   In abnormal situations,  LAMBDA
C    may only be near 1 and  X  near a fixed point or zero.
C
C IFLAG =
C
C   1   Normal return.
C
C   2   Specified error tolerance cannot be met.  Some or all of
C       ARCRE, ARCAE, ANSRE, ANSAE  have been increased to 
C       suitable values.  To continue, just call  FIXPQS  again
C       without changing any parameters.
C
C   3   STEPQS  has been called 1000 times.  To continue, call
C       FIXPQS  again without changing any parameters.
C
C   4   Jacobian matrix does not have full rank.  The algorithm
C       has failed (the zero curve of the homotopy map cannot be
C       followed any further).
C
C   5   The tracking algorithm has lost the zero curve of the 
C       homotopy map and is not making progress.  The error 
C       tolerances  ARC?E  and  ANS?E  were too lenient.  The problem 
C       should be restrarted by calling  FIXPQS  with smaller error 
C       tolerances and  IFLAG = 0 (-1, -2).
C
C   6   The Newton iteration in  STEPQS  or  ROOTNS  failed to
C       converge.  The error tolerances  ANS?E  may be too stringent.
C
C   7   Illegal input parameters, a fatal error.
C
C ARCRE, ARCAE, ANSRE, ANSAE  are unchanged after a normal return
C    (IFLAG = 1).  They are increased to appropriate values on the
C    return  IFLAG = 2.
C
C NFE  is the number of Jacobian evaluations.
C
C ARCLEN  is the approximate length of the zero curve.  
C
C
C Allocatable and automatic work arrays:
C
C YP(1:N+1)  is a work array containing the tangent vector to 
C    the zero curve at the current point  Y .
C
C YOLD(1:N+1)  is a work array containing the previous point found
C    on the zero curve.
C
C YPOLD(1:N+1)  is a work array containing the tangent vector to 
C    the zero curve at  YOLD .
C
C QR(1:LENQR), PP(1:N), ROWPOS(1:N+2), COLPOS(1:LENQR) are all work
C    arrays used to define the sparse Jacobian matrices, allocated
C    here, and distributed via the module  HOMOTOPY .
C
C
C
      USE HOMOTOPY, QR => QRSPARSE
      USE REAL_PRECISION
      INTEGER, INTENT(IN)::LENQR,MODE,N,TRACE
      REAL (KIND=R8), DIMENSION(:), INTENT(IN OUT)::A,Y
      INTEGER, INTENT(IN OUT)::IFLAG
      REAL (KIND=R8), INTENT(IN OUT)::ANSAE,ANSRE,ARCAE,ARCRE,SSPAR(4)
      INTEGER, INTENT(OUT)::NFE
      REAL (KIND=R8), INTENT(OUT)::ARCLEN
C
C     FUNCTION DECLARATIONS 
C 
      INTERFACE
        FUNCTION DNRM2(N,X,STRIDE)
        USE REAL_PRECISION
        INTEGER:: N,STRIDE
        REAL (KIND=R8):: DNRM2,X(N)
        END FUNCTION DNRM2
      END INTERFACE
C
C     LOCAL VARIABLES 
C
      REAL (KIND=R8), SAVE:: ABSERR, H, HOLD, RELERR, S, WK 
      INTEGER, SAVE:: IFLAGC, ITER, JW, LIMIT, NP1
      LOGICAL, SAVE:: CRASH, START       
C
C     WORK ARRAYS 
C
      REAL (KIND=R8), ALLOCATABLE, DIMENSION(:), SAVE:: YP,YOLD,YPOLD
      REAL (KIND=R8):: DZ(N+1),T(N+1),Z0(N+1) 
C
C LIMITD  IS AN UPPER BOUND ON THE NUMBER OF STEPS.  IT MAY BE 
C CHANGED BY CHANGING THE FOLLOWING PARAMETER STATEMENT:
      INTEGER, PARAMETER:: LIMITD = 1000
C
      INTERFACE
        SUBROUTINE ROOTNS(N,NFE,IFLAGC,ANSRE,ANSAE,Y,YP,
     &    YOLD,YPOLD,A,MODE,LENQR)
        USE HOMOTOPY, QR => QRSPARSE
        USE REAL_PRECISION
        INTEGER, INTENT(IN):: LENQR,MODE,N
        INTEGER, INTENT(IN OUT):: IFLAGC,NFE
        REAL (KIND=R8), INTENT(IN):: A(:)
        REAL (KIND=R8), INTENT(IN):: ANSAE,ANSRE
        REAL (KIND=R8), DIMENSION(:), INTENT(IN OUT):: Y,YOLD,YP,YPOLD
        END SUBROUTINE ROOTNS
C
        SUBROUTINE STEPQS(N,NFE,IFLAGC,MODE,LENQR,START,CRASH,HOLD,H,
     &    WK,RELERR,ABSERR,S,Y,YP,YOLD,YPOLD,A,Z0,DZ,T,SSPAR)
        USE HOMOTOPY, QR => QRSPARSE
        USE REAL_PRECISION
        INTEGER, INTENT(IN):: LENQR,MODE,N
        INTEGER, INTENT(IN OUT):: IFLAGC,NFE
        LOGICAL, INTENT(IN OUT):: CRASH,START
        REAL (KIND=R8), INTENT(IN):: A(:),SSPAR(4)
        REAL (KIND=R8), INTENT(IN OUT):: ABSERR,H,HOLD,RELERR,S,WK,
     &    Y(:),YOLD(:),YP(:),YPOLD(:)
        REAL (KIND=R8), INTENT(OUT), DIMENSION(:):: DZ,T,Z0
        END SUBROUTINE STEPQS
      END INTERFACE
C
C ***** FIRST EXECUTABLE STATEMENT *****
C
C CHECK IFLAG
C
      IF (N .LE. 0  .OR.  ANSRE .LE. 0.0  .OR.  ANSAE .LT. 0.0
     &  .OR.  N+1 .NE. SIZE(Y) .OR.
     &  ((IFLAG .EQ. -1  .OR.  IFLAG .EQ. 0) .AND.  N .NE. SIZE(A))
     &  .OR.  MODE .LE. 0  .OR.  MODE .GE. 3)
     &                                                     IFLAG=7
      IF (IFLAG .GE. -2 .AND. IFLAG .LE. 0) GO TO 10
      IF (IFLAG .EQ. 2) GO TO 50
      IF (IFLAG .EQ. 3) GO TO 40
C
C ONLY VALID INPUT FOR IFLAG IS -2, -1, 0, 2, 3.
C
      IFLAG = 7
      RETURN
C
C ***** INITIALIZATION BLOCK  *****
C
 10   ARCLEN = 0.0
      IF (ARCRE .LE. 0.0) ARCRE = .5*SQRT(ANSRE)
      IF (ARCAE .LE. 0.0) ARCAE = .5*SQRT(ANSAE)
      NFE=0
      IFLAGC = IFLAG
      NP1=N+1
C ALLOCATE SAVED WORK ARRAYS.
      IF (ALLOCATED(YOLD)) DEALLOCATE(YOLD)
      IF (ALLOCATED(YP)) DEALLOCATE(YP)
      IF (ALLOCATED(YPOLD)) DEALLOCATE(YPOLD)
      ALLOCATE(YOLD(NP1),YP(NP1),YPOLD(NP1))
C 
C SET INITIAL CONDITIONS FOR FIRST CALL TO STEPQS.
C
        START=.TRUE.
        CRASH=.FALSE.
        RELERR = ARCRE
        ABSERR = ARCAE
        HOLD=1.0
        H=0.1
        S=0.0
        YPOLD(NP1) = 1.0
        Y(NP1) = 0.0
        YPOLD(1:N)=0.0
C
C SET OPTIMAL STEP SIZE ESTIMATION PARAMETERS.
C
C     MINIMUM STEP SIZE HMIN
      IF (SSPAR(1) .LE. 0.0) SSPAR(1)=(SQRT(N+1.0)+4.0)*EPSILON(1.0_R8)
C     MAXIMUM STEP SIZE HMAX
      IF (SSPAR(2) .LE. 0.0) SSPAR(2)= 1.0
C     MINIMUM STEP REDUCTION FACTOR BMIN
      IF (SSPAR(3) .LE. 0.0) SSPAR(3)= 0.1_R8
C     MAXIMUM STEP EXPANSION FACTOR BMAX
      IF (SSPAR(4) .LE. 0.0) SSPAR(4)= 7.0
C
C LOAD  A  FOR THE FIXED POINT AND ZERO FINDING PROBLEMS.
C
      IF (IFLAGC .GE. -1) A(1:N) = Y(1:N)
C
 40   LIMIT=LIMITD
C ALLOCATE ARRAYS FOR SPARSE JACOBIAN MATRIX DATA STRUCTURE.
 50   SELECT CASE (MODE)
        CASE (1)
          IF (.NOT. ALLOCATED(QR)) ALLOCATE(QR(LENQR))
          IF (.NOT. ALLOCATED(ROWPOS)) ALLOCATE(ROWPOS(N+2))
          IF (.NOT. ALLOCATED(PP)) ALLOCATE(PP(N))
        CASE (2)
          IF (.NOT. ALLOCATED(QR)) ALLOCATE(QR(LENQR))
          IF (.NOT. ALLOCATED(ROWPOS)) ALLOCATE(ROWPOS(N+2))
          IF (.NOT. ALLOCATED(COLPOS)) ALLOCATE(COLPOS(LENQR))
          IF ((.NOT. ALLOCATED(PP)) .AND. (IFLAGC .GE. -1))
     &      ALLOCATE(PP(N))
      END SELECT
C
C ***** END OF INITIALIZATION BLOCK. *****
C
      MAIN_LOOP: DO ITER=1,LIMIT ! ***** MAIN LOOP. *****
        IF (Y(NP1) .LT. 0.0) THEN
          ARCLEN = S
          IFLAG = 5
          CALL CLEANUPALL
          RETURN
        END IF
C
C TAKE A STEP ALONG THE CURVE.
C
        CALL STEPQS(N,NFE,IFLAGC,MODE,LENQR,START,CRASH,HOLD,H,
     &    WK,RELERR,ABSERR,S,Y,YP,YOLD,YPOLD,A,Z0,DZ,T,SSPAR)
C
C PRINT LATEST POINT ON CURVE IF REQUESTED.
C
        IF (TRACE .GT. 0) THEN
          WRITE (TRACE,217) ITER,NFE,S,Y(NP1),(Y(JW),JW=1,N)
217       FORMAT(/' STEP',I5,3X,'NFE =',I5,3X,'ARC LENGTH =',F9.4,3X,
     &    'LAMBDA =',F7.4,5X,'X vector:'/(1X,6ES12.4))
        ENDIF
C
C CHECK IF THE STEP WAS SUCCESSFUL.
C
        IF (IFLAGC .GT. 0) THEN
          ARCLEN=S
          IFLAG=IFLAGC
          CALL CLEANUPALL
          RETURN
        END IF
C
        IF (CRASH) THEN
C
C         RETURN CODE FOR ERROR TOLERANCE TOO SMALL.
C      
          IFLAG=2
C
C         CHANGE ERROR TOLERANCES.
C
          IF (ARCRE .LT. RELERR) THEN
            ARCRE=RELERR
            ANSRE=RELERR
          ENDIF
          IF (ARCAE .LT. ABSERR) ARCAE=ABSERR
C
C         CHANGE LIMIT ON NUMBER OF ITERATIONS.
C
          LIMIT = LIMIT - ITER
          CALL CLEANUP
          RETURN
        END IF
C
C IF  LAMBDA >= 1.0,  USE  ROOTNS  TO FIND SOLUTION.
C
        IF (Y(NP1) .GE. 1.0) GO TO 500
C
      END DO MAIN_LOOP   ! ***** END OF MAIN LOOP *****
C
C DID NOT CONVERGE IN  LIMIT  ITERATIONS, SET  IFLAG  AND RETURN.
C
      ARCLEN = S
      IFLAG = 3
      CALL CLEANUP
      RETURN
C
C ***** FINAL STEP -- FIND SOLUTION AT LAMBDA=1 *****
C
C SAVE  YOLD  FOR ARC LENGTH CALCULATION LATER.
C
 500  T = YOLD
C
C FIND SOLUTION.
C
      CALL ROOTNS(N,NFE,IFLAGC,ANSRE,ANSAE,Y,YP,
     &    YOLD,YPOLD,A,MODE,LENQR)
C
C CHECK IF SOLUTION WAS FOUND AND SET  IFLAG  ACCORDINGLY.
C
      IFLAG=1
C
C     SET ERROR FLAG IF ROOTNS COULD NOT GET THE POINT ON THE ZERO
C     CURVE AT  LAMBDA = 1.0 .
C
      IF (IFLAGC .GT. 0) IFLAG=IFLAGC
C
C CALCULATE FINAL ARC LENGTH.
C
      DZ = Y - T
      ARCLEN = S - HOLD + DNRM2(NP1,DZ,1)
C
C ***** END OF FINAL STEP *****
C
      CALL CLEANUPALL
      RETURN
C
      CONTAINS
        SUBROUTINE CLEANUPALL
        IF (ALLOCATED(YOLD)) DEALLOCATE(YOLD)
        IF (ALLOCATED(YP)) DEALLOCATE(YP)
        IF (ALLOCATED(YPOLD)) DEALLOCATE(YPOLD)
        CALL CLEANUP
        RETURN
        END SUBROUTINE CLEANUPALL
        SUBROUTINE CLEANUP
        IF (ALLOCATED(QR)) DEALLOCATE(QR)
        IF (ALLOCATED(ROWPOS)) DEALLOCATE(ROWPOS)
        IF (ALLOCATED(COLPOS)) DEALLOCATE(COLPOS)
        IF (ALLOCATED(PP)) DEALLOCATE(PP)
        IF (ALLOCATED(PAR)) DEALLOCATE(PAR)
        IF (ALLOCATED(IPAR)) DEALLOCATE(IPAR)
        RETURN
        END SUBROUTINE CLEANUP
      END SUBROUTINE FIXPQS
!
      SUBROUTINE POLSYS1H(N,NUMT,COEF,KDEG,IFLG1,IFLG2,EPSBIG,EPSSML,
     &     SSPAR,NUMRR,LAMBDA,ROOTS,ARCLEN,NFE)
C
C POLSYS1H finds all (complex) solutions to a system
C F(X)=0 of N polynomial equations in N unknowns
C with real coefficients. If IFLG=10 or IFLG=11, POLSYS1H
C returns the solutions at infinity also.
C
C The system F(X)=0 is described via the coefficents,
C "COEF", and the parameters "N, NUMT, KDEG", as follows.
C
C
C       NUMT(J)
C
C F(J) = SUM  COEF(J,K) * X(1)**KDEG(J,1,K)...X(N)**KDEG(J,N,K)
C
C        K=1
C
C FOR J=1, ..., N.
C
C
C POLSYS1H has two main run options:  automatic scaling and
C the projective transformation.  These are evoked via the
C flag "IFLG1", as described below.  The other input
C parameters are the same whether one or both of these options
C are specified, and the output is always returned unscaled
C and untransformed.
C
C If automatic scaling is specified, then the input
C coefficients are modified by subroutine  SCLGNP . The problem
C is solved with the scaled coefficients and scaled variables.
C The coefficients are returned scaled.
C
C If the projective transformation is specified, then
C essentially the system is reformulated in homogeneous
C coordinates, Z(1), ..., Z(N+1), and solved in complex
C projective space.  The resulting solutions are
C untransformed via
C
C X(J) = Z(J)/Z(N+1)   J=1, ..., N.
C
C On return,
C
C ROOTS(1,J,M) = real part of X(J) for the Mth path,
C
C ROOTS(2,J,M) = imaginary part of X(J) for the Mth path,
C
C for J=1, ..., N, and
C
C ROOTS(1,N+1,M) = real part of Z(N+1) for the Mth path,
C
C ROOTS(2,N+1,M) = imaginary part of Z(N+1) for the Mth path.
C
C If ROOTS(*,N+1,M) is small, then the associated solution
C should be regarded as being "near infinity".  Note that,
C when the projective transformation has been specified, the
C ROOTS values have been untransformed -- that is, divided
C through by Z(N+1) -- unless such division would have caused
C overflow.  In this latter case, the affected components of
C ROOTS are set to the largest floating point number (machine
C infinity).
C
C The code can be modified easily to solve systems with complex
C coefficients,  COEF .  Only the subroutines  INITP  and  FFUNP
C need be changed.
C
C The FORTRAN COMPLEX declaration is not used in POLSYS1H.
C Complex variables are represented by real arrays with first
C index dimensioned 2 and complex operations are evoked by
C subroutine calls.
C
C The total number of paths that will be tracked (if
C IFLG2(M)=-2 for all M) is equal to the "total degree" of the
C system, TOTDG.   TOTDG is equal to the products of the
C degrees of all the equations in the system.  The degree of
C an equation is the maximum of the degrees of its terms.  The
C degree of a term is the sum of the degrees of the variables.
C Thus, TOTDG = IDEG(1) * ... * IDEG(N) where IDEG(J) =
C MAX {JDEG(J,K) | K=1,...,NUMT(J)} where JDEG(J,K) = KDEG(J,1,K) +
C ... + KDEG(J,N,K).
C
C IFLG1  determines whether the system is to be automatically
C scaled by  POLSYS1H  and whether the projective transformation
C of the system is to be automatically evoked by POLSYS1H.  See
c "ON INPUT" below.
C
C IFLG2, EPSBIG, EPSSML, and  SSPAR  tell the path tracker
C FIXPNF  which paths to track and set parameters for the path
C tracker.
C
C NUMRR  tells  POLSYS1H  how many multiples of 1000 steps to try
C before abandoning a path.
C
C The output consists of  IFLG1, and of  LAMBDA, ROOTS, ARCLEN, and
C NFE  for each path.  IFLG1  returns input data error information.
C ROOTS  gives the solutions themselves, while  LAMBDA, ARCLEN,
C and  NFE  give information about the associated paths.
C
C
C The following subroutines are used directly or indirectly by
C POLSYS1H: 
C         Special for POLSYS1H:
C           INITP , STRPTP , OTPUTP , RHO , RHOJAC ,
C           HFUNP , HFUN1P , GFUNP , FFUNP ,
C           MULP , POWP , DIVP , SCLGNP .
C         From the general HOMPACK routines:
C           FIXPNF , ROOT , ROOTNF , STEPNF , TANGNF .
C         From LAPACK routines:
C           DGEQPF , DGEQRF , DORMQR .
C         From BLAS routines:
C           DCOPY ,  DDOT ,  DGEMM ,  DGEMV ,  DGER ,  
C           DNRM2 ,  DSCAL ,  DSWAP ,  DTRMM ,  DTRMV , DTRSV ,
C           IDAMAX ,  LSAME , XERBLA . 
C
C ON INPUT:
C
C N  is the number of equations and variables.
C
C NUMT(1:N)  is an integer array.  NUMT(J)  is the number of terms
C   in the Jth equation for J=1 to N.
C
C COEF(1:N,1:)  is a real array.  COEF(J,K)  is 
C   the Kth coefficient of the Jth equation for J=1 to N,
C   K=1 to NUMT(J).  The second dimension must be greater than or equal
C   to the maximum number of terms in each equation.  In other words,
C   SIZE(COEF,DIM=2) .GE. MAXT = MAX {NUMT(J) | J=1, ..., N} .
C
C KDEG(1:N,1:N+1,1:)  is an integer array.  
C   KDEG(J,L,K)  is the degree of the Lth variable in the Kth
C   term of the Jth equation for  J=1 to N, L=1 to N, K=1 to NUMT(J).
C   SIZE(KDEG,DIM=3) .GE. MAXT = MAX {NUMT(J) | J=1, ..., N} .
C
C IFLG1 =
C   00  if the problem is to be solved without
C       calling POLSYS1H' scaling routine, SCLGNP, and
C       without using the projective transformtion.
C
C   01  if scaling but no projective transformation is to be used.
C
C   10  if no scaling but projective transformation is to be used.
C
C   11  if both scaling and projective transformation are to be used.
C
C IFLG2(1:TOTDG)  is an integer array.  If IFLG2(M) = -2, then the 
C   Mth path is tracked.  Otherwise the Mth path is skipped.
C   Thus, to find all solutions set IFLG2(M) = -2 for M=1,...,TOTDG.
C   Selected paths can be rerun by setting IFLG2(M)=-2 for
C   the paths to be rerun and IFLG2(M).NE.-2 for the others.
C
C EPSBIG  is the local error tolerance allowed the path tracker along
C   the path.  ARCRE and ARCAE (in  FIXPNF ) are set to  EPSBIG.
C
C EPSSML  is the accuracy desired for the final solution.  ANSRE and
C   ANSAE (in  FIXPNF ) are set to  EPSSML.
C
C SSPAR(1:8) = (LIDEAL, RIDEAL, DIDEAL, HMIN, HMAX, BMIN, BMAX, P)  is
C    a vector of parameters used for the optimal step size estimation.
C    If  SSPAR(J) .LE. 0.0  on input, it is reset to a default value
C    by  FIXPNF .  Otherwise the input value of  SSPAR(J)  is used.
C    See the comments in  FIXPNF  and in  STEPNF  for more information
C    about these constants.
C
C NUMRR  is the number of multiples of 1000 steps that will be tried
C   before abandoning a path.
C
C
C ON OUTPUT:
C
C N, NUMT, COEF, KDEG, EPSBIG, EPSSML, and NUMRR are unchanged.
C
C IFLG1=
C   -1  if  NUMT  is incorrectly dimensioned or invalid.
C   -2  if  COEF  is incorrectly dimensioned.
C   -3  if  KDEG  is incorrectly dimensioned or invalid.
C   -4  if any of  IFLG2, LAMBDA, ROOTS, ARCLEN, or  NFE  are
C       incorrectly dimensioned.
C   -5  if the global work arrays  IPAR  and  PAR  could not be
C       allocated.
C   -6  if  IFLG1  on input is not 00 or 01 or 10 or 11.
C   Unchanged otherwise.
C
C IFLG2(1:TOTDG)  gives information about how the Mth path terminated:
C IFLG2(M) =
C   1   Normal return.
C
C   2   Specified error tolerance cannot be met.  Increase  EPSBIG
C       and  EPSSML  and rerun.
C
C   3   Maximum number of steps exceeded.  To track the path further,
C       increase  NUMRR  and rerun the path.  However, the path may
C       be diverging, if the  LAMBDA  value is near 1 and the  ROOTS 
C       values are large.
C
C   4   Jacobian matrix does not have full rank.  The algorithm
C       has failed (the zero curve of the homotopy map cannot be
C       followed any further).
C
C   5   The tracking algorithm has lost the zero curve of the
C       homotopy map and is not making progress.  The error tolerances
C       EPSBIG  and  EPSSML  were too lenient.  The problem should be
C       restarted with smaller error tolerances.
C
C   6   The normal flow Newton iteration in  STEPNF  or  ROOTNF
C       failed to converge.  The error tolerances  EPSBIG  or  EPSSML
C       may be too stringent.
C
C   7   Illegal input parameters, a fatal error.
C
C LAMBDA(M)  is the final LAMBDA value for the Mth path, M = 1, ...,
C   TOTDG, where LAMBDA is the continuation parameter.
C
C ROOTS(1,J,M), ROOTS(2,J,M)  are the real and imaginary parts
C   of the Jth variable respectively, for J = 1,...,N, for
C   the Mth path, for M = 1,...,TOTDG.  If  IFLG1 = 10 or 11, then
C   ROOTS(1,N+1,M)  and  ROOTS(2,N+1,M)  are the real and
C   imaginary parts respectively of the projective
C   coordinate of the solution.
C
C ARCLEN(M)  is the arc length of the Mth path for M = 1, ..., TOTDG.
C
C NFE(M)  is the number of Jacobian matrix evaluations required to 
C   track the Mth path for M =1, ..., TOTDG.
C
C ----------------------------------------------------------------------
      USE HOMOTOPY
      USE REAL_PRECISION
      INTERFACE
        SUBROUTINE INITP(IFLG1,N,NUMT,KDEG,COEF,
     &                              IDEG,FACV,CL,PDG,QDG,R)
        USE HOMOTOPY
        USE REAL_PRECISION
        INTEGER, INTENT(IN):: IFLG1,N,NUMT(:)
        INTEGER, INTENT(IN OUT):: KDEG(:,:,:)
        REAL (KIND=R8), INTENT(IN OUT):: COEF(:,:)
        INTEGER, INTENT(OUT):: IDEG(N)
        REAL (KIND=R8), INTENT(OUT):: 
     &    FACV(N),CL(2,N+1),PDG(2,N),QDG(2*N),R(2,N)
        END SUBROUTINE INITP
C
        SUBROUTINE STRPTP(N,ICOUNT,IDEG,R,X)
        USE REAL_PRECISION
        INTEGER:: N,ICOUNT(N),IDEG(N)
        REAL (KIND=R8):: R(2,N),X(2*N)
        END SUBROUTINE STRPTP
C
!       SUBROUTINE FIXPNF(N,Y,IFLAG,ARCRE,ARCAE,ANSRE,ANSAE,TRACE,A,
!    &    SSPAR,NFE,ARCLEN,POLY_SWITCH)
!       USE REAL_PRECISION
!       INTEGER, INTENT(IN)::N,TRACE
!       REAL (KIND=R8), DIMENSION(:), INTENT(IN OUT)::A,Y
!       INTEGER, INTENT(IN OUT)::IFLAG
!       REAL (KIND=R8), INTENT(IN OUT)::ANSAE,ANSRE,ARCAE,ARCRE,
!    &    SSPAR(8)
!       INTEGER, INTENT(OUT)::NFE
!       REAL (KIND=R8), INTENT(OUT)::ARCLEN
!       LOGICAL, INTENT(IN), OPTIONAL::POLY_SWITCH
!       END SUBROUTINE FIXPNF
C
        SUBROUTINE OTPUTP(N,NUMPAT,CL,FACV,CLX,X,XNP1)
        USE REAL_PRECISION
        INTEGER, INTENT(IN):: N,NUMPAT
        REAL (KIND=R8), INTENT(IN):: CL(2,N+1),FACV(N)
        REAL (KIND=R8), INTENT(IN OUT):: CLX(2,N),X(2,N),XNP1(2)
        END SUBROUTINE OTPUTP
      END INTERFACE
C
C TYPE DECLARATIONS FOR INPUT AND OUTPUT
C
      INTEGER, INTENT(IN):: N,NUMT(:),NUMRR
      REAL (KIND=R8), INTENT(IN OUT):: COEF(:,:),SSPAR(8)
      INTEGER, INTENT(IN OUT):: KDEG(:,:,:),IFLG1,IFLG2(:)
      REAL (KIND=R8), INTENT(IN):: EPSBIG,EPSSML
      REAL (KIND=R8), INTENT(OUT):: LAMBDA(:),ROOTS(:,:,:),ARCLEN(:)
      INTEGER, INTENT(OUT):: NFE(:)
C
C TYPE DECLARATIONS FOR LOCAL VARIABLES
C
      INTEGER:: I,ICOUNT(N),IDEG(N),IDUMMY,IFLAG,IJ,
     &  IPROFF(15),J,LIPAR(15),LPAR(25),MAXT,N2,N2P1,
     &  NNFE,NP1,NUMPAT,PROFF(25),TOTDG,TRACE
      REAL (KIND=R8):: AARCLN,ANSAE,ANSRE,ARCAE,ARCRE,CL(2,N+1),
     &  FACV(N),PDG(2,N),QDG(2*N),R(2,N),XNP1(2),Y(2*N+1)
C
C ----------------------------------------------------------------------
      N2=2*N
      NP1=N+1
      N2P1=N2+1
C
C CHECK THAT DIMENSIONS ARE VALID.
C
      IF ((SIZE(NUMT) /= N) .OR. ANY(NUMT .LE. 0)) THEN
        IFLG1=-1
        RETURN
      END IF
      MAXT = MAXVAL(NUMT)
      IF ((SIZE(COEF,DIM=1) /= N) .OR. (SIZE(COEF,DIM=2) < MAXT)) THEN
        IFLG1=-2
        RETURN
      END IF
      KDEG = ABS(KDEG)
      IF ((SIZE(KDEG,DIM=1) /= N) .OR. (SIZE(KDEG,DIM=2) /= NP1) .OR.
     &  (SIZE(KDEG,DIM=3) < MAXT) ) THEN
        IFLG1=-3
        RETURN
      END IF
      DO J=1,N
        IDEG(J)=MAXVAL(SUM(KDEG(J,1:N,1:NUMT(J)),DIM=1))
      END DO
      TOTDG = PRODUCT(IDEG)
      IF ((SIZE(IFLG2) < TOTDG) .OR. (SIZE(LAMBDA) < TOTDG) .OR.
     &  (SIZE(ROOTS,DIM=3) < TOTDG) .OR. (SIZE(ARCLEN) < TOTDG) .OR.
     &  (SIZE(NFE) < TOTDG) .OR. 
     &  (IFLG1 <= 1 .AND. SIZE(ROOTS,DIM=2) /= N) .OR.
     &  (IFLG1 >= 10 .AND. SIZE(ROOTS,DIM=2) /= NP1)) THEN
        IFLG1=-4
        RETURN
      END IF
      IF (IFLG1 /= 0 .AND. IFLG1 /= 1 .AND.
     &  IFLG1 /= 10 .AND. IFLG1 /= 11) THEN
        IFLG1=-6
        RETURN
      END IF
C
C ALLOCATE THE GLOBAL WORK ARRAYS  IPAR  AND  PAR, USED TO COMMUNICATE
C DATA BETWEEN SUBROUTINES VIA THE MODULE HOMOTOPY.
C
      ALLOCATE(IPAR(42 + 2*N + N*(N+1)*MAXT),
     &  PAR(2 + 28*N + 6*N**2 + 7*N*MAXT + 4*N**2*MAXT),STAT=IJ)
      IF (IJ .NE. 0) THEN
        IFLG1=-5
        RETURN
      END IF
C      
C INITIALIZATION
C
      CALL INITP(IFLG1,N,NUMT,KDEG,COEF,
     &                              IDEG,FACV,CL,PDG,QDG,R)
C
C INTEGER VARIABLES AND ARRAYS TO BE PASSED IN IPAR:
C
C    IPAR INDEX     VARIABLE NAME       LENGTH
C    ----------     -------------    -----------------
C          1                N               1
C          2             MAXT               1
C          3            PROFF               25
C          4           IPROFF               15
C          5             IDEG               N
C          6             NUMT               N
C          7             KDEG               N*(N+1)*MAXT
C
C
C DOUBLE PRECISION VARIABLES AND ARRAYS TO BE PASSED IN PAR:
C
C     PAR INDEX     VARIABLE NAME       LENGTH
C    ----------     -------------    -----------------
C          1              PDG               2*N
C          2               CL               2*(N+1)
C          3             COEF               N*MAXT
C          4                H               N2
C          5              DHX               N2*N2
C          6              DHT               N2
C          7            XDGM1               2*N
C          8              XDG               2*N
C          9              G                 2*N
C         10             DG                 2*N
C         11           PXDGM1               2*N
C         12             PXDG               2*N
C         13               F                2*N
C         14              DF                2*N*(N+1)
C         15               XX               2*N*(N+1)*MAXT
C         16              TRM               2*N*MAXT
C         17             DTRM               2*N*(N+1)*MAXT
C         18              CLX               2*N
C         19            DXNP1               2*N
C
C SET LENGTHS OF VARIABLES
      LIPAR(1)=1
      LIPAR(2)=1
      LIPAR(3)=25
      LIPAR(4)=15
      LIPAR(5)=N
      LIPAR(6)=N
      LIPAR(7)=N*(N+1)*MAXT
      LPAR( 1)=2*N
      LPAR( 2)=2*NP1
      LPAR( 3)=N*MAXT
      LPAR( 4)=N2
      LPAR( 5)=N2*N2
      LPAR( 6)=N2
      LPAR( 7)=2*N
      LPAR( 8)=2*N
      LPAR( 9)=2*N
      LPAR(10)=2*N
      LPAR(11)=2*N
      LPAR(12)=2*N
      LPAR(13)=2*N
      LPAR(14)=2*N*NP1
      LPAR(15)=2*N*NP1*MAXT
      LPAR(16)=2*N*MAXT
      LPAR(17)=2*N*NP1*MAXT
      LPAR(18)=2*N
      LPAR(19)=2*N
C
C PROFF AND IPROFF ARE OFFSETS THAT DEFINE THE VARIABLES LISTED ABOVE
      PROFF(1)=1
      DO I=2,19
          PROFF(I)=PROFF(I-1)+LPAR(I-1)
      END DO
      IPROFF(1)=1
      DO I=2,7
          IPROFF(I)=IPROFF(I-1)+LIPAR(I-1)
      END DO
C
C DEFINE VARIABLES
      IPAR(1)=N
      IPAR(2)=MAXT
      IPAR(IPROFF(3):IPROFF(3)+18) = PROFF(1:19)
      IPAR(IPROFF(4):IPROFF(4)+ 6) = IPROFF(1:7)
      IPAR(IPROFF(5):IPROFF(5)+N-1) = IDEG(1:N)
      IPAR(IPROFF(6):IPROFF(6)+N-1) = NUMT(1:N)
      IPAR(IPROFF(7):IPROFF(7)+LIPAR(7)-1) =
     &  PACK(KDEG(:,:,1:MAXT),.TRUE.)
      PAR(PROFF(1):PROFF(1)+LPAR(1)-1) = PACK(PDG,.TRUE.)
      PAR(PROFF(2):PROFF(2)+LPAR(2)-1) = PACK(CL,.TRUE.)
      PAR(PROFF(3):PROFF(3)+LPAR(3)-1) = PACK(COEF(:,1:MAXT),.TRUE.)
C
C ICOUNT IS A COUNTER USED BY "STRPTP"
      ICOUNT(1)=0
      ICOUNT(2:N)=1
C
C PATHS LOOP -- ITERATE THROUGH PATHS
C
      PATHS: DO NUMPAT = 1,TOTDG
C         GET A START POINT, Y, FOR THE PATH.
          Y(1) = 0.0
          CALL STRPTP(N,ICOUNT,IDEG,R,Y(2:N2P1))
C         CHECK WHETHER PATH IS TO BE FOLLOWED.
          IFLAG = IFLG2(NUMPAT)
          IF (IFLAG .NE. -2) CYCLE PATHS
          ARCRE = EPSBIG
          ARCAE = ARCRE
          ANSRE = EPSSML
          ANSAE = ANSRE
          TRACE = 0
C         TRACK A HOMOTOPY PATH.
          DO IDUMMY=1,MAX(NUMRR,1)
            CALL FIXPNF(N2,Y,IFLAG,ARCRE,ARCAE,ANSRE,ANSAE,TRACE,
     &        QDG,SSPAR,NNFE,AARCLN, POLY_SWITCH=.TRUE.)
            IF (IFLAG .NE. 2 .AND. IFLAG .NE. 3) EXIT
          END DO
C         UNSCALE AND UNTRANSFORM COMPUTED SOLUTION.
          CALL OTPUTP(N,NUMPAT,CL,FACV,
     &      PAR(PROFF(18):PROFF(18)+LPAR(18)-1),Y(2:N2P1),XNP1)
          LAMBDA(NUMPAT) = Y(1)
          ROOTS(1,1:N,NUMPAT) = Y(2:N2P1:2)
          ROOTS(2,1:N,NUMPAT) = Y(3:N2P1:2)
          ROOTS(1:2,NP1,NUMPAT) = XNP1
C
          ARCLEN(NUMPAT)= AARCLN
          NFE(NUMPAT)   = NNFE
          IFLG2(NUMPAT) = IFLAG
      END DO PATHS
C CLEAN UP WORK SPACE.
      IF (ALLOCATED(IPAR)) DEALLOCATE(IPAR)
      IF (ALLOCATED(PAR))  DEALLOCATE(PAR)
      RETURN
      END SUBROUTINE POLSYS1H
      END MODULE HOMPACK90


      SUBROUTINE DIVP(XXXX,YYYY,ZZZZ,IERR)
C
C THIS SUBROUTINE PERFORMS DIVISION  OF COMPLEX NUMBERS:
C ZZZZ = XXXX/YYYY
C
C ON INPUT:
C
C XXXX  IS AN ARRAY OF LENGTH TWO REPRESENTING THE FIRST COMPLEX
C       NUMBER, WHERE XXXX(1) = REAL PART OF XXXX AND XXXX(2) =
C       IMAGINARY PART OF XXXX.
C
C YYYY  IS AN ARRAY OF LENGTH TWO REPRESENTING THE SECOND COMPLEX
C       NUMBER, WHERE YYYY(1) = REAL PART OF YYYY AND YYYY(2) =
C       IMAGINARY PART OF YYYY.
C
C ON OUTPUT:
C
C ZZZZ  IS AN ARRAY OF LENGTH TWO REPRESENTING THE RESULT OF
C       THE DIVISION, ZZZZ = XXXX/YYYY, WHERE ZZZZ(1) =
C       REAL PART OF ZZZZ AND ZZZZ(2) = IMAGINARY PART OF ZZZZ.
C
C IERR =
C  1   IF DIVISION WOULD HAVE CAUSED OVERFLOW.  IN THIS CASE, THE
C      APPROPRIATE PARTS OF ZZZZ ARE SET EQUAL TO THE LARGEST
C      FLOATING POINT NUMBER, AS GIVEN BY FUNCTION  HUGE .
C
C  0   IF DIVISION DOES NOT CAUSE OVERFLOW.
C
C DECLARATION OF INPUT
      USE REAL_PRECISION
      REAL (KIND=R8), DIMENSION(2), INTENT(IN):: XXXX,YYYY
C
C DECLARATION OF OUTPUT
      INTEGER, INTENT(OUT):: IERR
      REAL (KIND=R8), DIMENSION(2), INTENT(OUT):: ZZZZ
C
C DECLARATION OF VARIABLES
      REAL (KIND=R8):: DENOM,XNUM
C
      IERR = 0
      DENOM = YYYY(1)*YYYY(1) + YYYY(2)*YYYY(2)
      XNUM    =   XXXX(1)*YYYY(1) + XXXX(2)*YYYY(2)
      IF (ABS(DENOM) .GE. 1.0  .OR.  ( ABS(DENOM) .LT. 1.0   .AND.
     & ABS(XNUM)/HUGE(1.0_R8) .LT. ABS(DENOM) ) ) THEN
            ZZZZ(1) = XNUM/DENOM
          ELSE
            ZZZZ(1) = HUGE(1.0_R8)
            IERR =1
          END IF
      XNUM    =   XXXX(2)*YYYY(1) - XXXX(1)*YYYY(2)
      IF (ABS(DENOM) .GE. 1.0  .OR.  ( ABS(DENOM) .LT. 1.0   .AND.
     & ABS(XNUM)/HUGE(1.0_R8) .LT. ABS(DENOM) ) ) THEN
            ZZZZ(2) = XNUM/DENOM
          ELSE
            ZZZZ(2) = HUGE(1.0_R8)
            IERR =1
          END IF
      RETURN
      END SUBROUTINE DIVP
      SUBROUTINE FFUNP(N,NUMT,MAXT,KDEG,COEF,CL,X,
     &  XX,TRM,DTRM,CLX,DXNP1,F,DF)
C
C FFUNP  EVALUATES THE SYSTEM "F(X)=0" AND ITS PARTIAL
C DERIVATIVES, USING THE "TABLEAU" INPUT: N,NUMT,KDEG,COEF.
C
C FFUNP  CAN BE MADE MORE EFFICIENT BY CUSTOMIZING IT TO
C PARTICULAR SYSTEM TYPES.  FOR EXAMPLE,        
C IF X(1)**2 AND X(1)**3 ARE USED IN SEVERAL
C EQUATIONS, THE CURRENT  FFUNP  RECOMPUTES BOTH OF THESE FOR
C EACH EQUATION.  BUT (OF COURSE) WE CAN COMPUTE
C X1SQ=X(1)**2 AND X1CU=XSQ(1)*X(1), AND 
C USE THESE IN EACH OF THE EQUATIONS.
C
C THE PART OF THE CODE BELOW LABELED "BLOCK A" CAN BE
C CUSTOMIZED IN THIS WAY.   (THE CODE OUTSIDE OF
C BLOCK A CONCERNS THE PROJECTIVE TRANSFORMATION AND NEED NOT
C BE CHANGED.)  HOWEVER, BLOCK A REQUIRES THE HOMOGENEOUS FORM
C OF THE POLYNOMIALS RATHER THAN THE STANDARD FORM.  FURTHER,
C THE PARTIAL DERIVATIVES WITH RESPECT TO ALL N+1 PROJECTIVE
C VARIABLES MUST BE COMPUTED.  MORE EXPLICITLY,
C THE ORIGINAL SYSTEM, F(X)=0, IS GIVEN IN "NON-HOMOGENEOUS FORM" AS
C DESCRIBED IN SUBROUTINE POLSYS.  F(X)  IS 
C REPRESENTED IN "HOMOGENEOUS FORM" AS FOLLOWS:
C
C              NUMT(J)
C
C    F(J) =     SUM   TRM(J,K)
C
C               K=1
C
C WHERE  TRM(J,K)=COEF(J,K) * XX(J,1,K)*XX(J,2,K)* ... *XX(J,N+1,K)
C
C WITH XX(J,L,K) = X(L)**KDEG(J,L,K) FOR J=1 TO N, L=1 TO N, AND
C K=1 TO NUMT(J), AND WITH XX(J,N+1,K) = XNP1**KDEG(J,N+1,K) FOR J=1 TO
C N AND K=1 TO NUMT(J), WHERE  XNP1  IS THE "HOMOGENEOUS COORDINATE,"
C KDEG(J,N+1,K)=IDEG(J)-(KDEG(J,1,K)+ ... + KDEG(J,N,K)),
C AND IDEG(J) THE DEGREE OF THE J-TH EQUATION.   XNP1  IS GENERATED
C FROM  X  AND  CL  BEFORE BLOCK A.
C
C IN THIS DISCUSSION WE HAVE OMITTED, FOR SIMPLICITY OF 
C EXPOSITION, THE LEADING INDEX, WHICH DIFFERENTIATES THE 
C REAL AND IMAGINARY PARTS.  HOWEVER, THIS INDEX MUST NOT BE 
C OMITTED IN THE CODE.  
C
C WE COMPLETE THE EXPOSITION OF "REPLACING BLOCK A WITH MORE EFFICIENT
C CODE" WITH AN EXPLICIT EXAMPLE.  FIRST, THE SYSTEM IS DESCRIBED.
C THEN THE CODE THAT SHOULD BE USED IS GIVEN (COMMENTED OUT).
C IN TESTS  POLSYS  WITH THE MORE EFFICIENT  FFUNP  RAN ABOUT TWICE AS
C FAST AS WITH THE GENERIC  FFUNP .
C
C HERE IS THE SYSTEM TO BE SOLVED:
C         
C     F(1) = COEF(1,1) * X(1)**4
C    &     + COEF(1,2) * X(1)**3 * X(2) 
C    &     + COEF(1,3) * X(1)**3
C    &     + COEF(1,4) * X(1)
C    &     + COEF(1,5)
C     F(2) = COEF(2,1) * X(1)     * X(2)**2
C    &     + COEF(2,2)              X(2)**2
C    &     + COEF(2,3) 
C
C THE REPLACEMENT CODE REQUIRES THE FOLLOWING DECLARATIONS:
C     DOUBLE PRECISION X1SQ,X1CU,X2SQ,X3SQ,X3CU,
C    &  TEMPA,TEMPB,TEMPC,TEMPD,TEMPE,TEMPF
C     DIMENSION X1SQ(2),X1CU(2),X2SQ(2),X3SQ(2),X3CU(2),
C    &  TEMPA(2),TEMPB(2),TEMPC(2),TEMPD(2),TEMPE(2),TEMPF(2)
C
C HERE IS CODE TO REPLACE BLOCK A:
C
C******************  BEGIN BLOCK A  *******************
C
C     CALL MULP(X(1,1),X(1,1),X1SQ)
C     CALL MULP(X1SQ  ,X(1,1),X1CU)
C     CALL MULP(X(1,2),X(1,2),X2SQ)
C     CALL MULP(XNP1,  XNP1,  X3SQ)
C     CALL MULP(X3SQ  ,XNP1,  X3CU)
C     
C     DO 1 I=1,2
C       TEMPA(I)=   COEF(1,1) * X(I,1)
C    &            + COEF(1,2) * X(I,2) 
C    &            + COEF(1,3) * XNP1(I)
C       TEMPB(I)=   COEF(1,4) * X(I,1)
C    &            + COEF(1,5) * XNP1(I) 
C 1   CONTINUE
C
C     CALL MULP(X1SQ,  TEMPA,TEMPC)
C     CALL MULP(X(1,1),TEMPC,TEMPD)
C     CALL MULP(X3SQ,  TEMPB,TEMPE)
C     CALL MULP(XNP1,  TEMPE,TEMPF)
C     
C     DO 2 I=1,2
C       F(I,1)=TEMPD(I) + TEMPF(I)
C       DF(I,1,1)= 3. *TEMPC(I) + COEF(1,1)*X1CU(I) + COEF(1,4)*X3CU(I)
C       DF(I,1,2)= COEF(1,2) * X1CU(I)
C       DF(I,1,3)= COEF(1,3)*X1CU(I) + 3. *TEMPE(I) + COEF(1,5)*X3CU(I) 
C
C       TEMPA(I) = COEF(2,1) * X(I,1) + COEF(2,2) * XNP1(I)
C  2  CONTINUE
C
C     CALL MULP(TEMPA,X(1,2),TEMPB)
C     CALL MULP(TEMPB,X(1,2),TEMPC)
C
C     DO 3 I=1,2
C       F(I,2) = TEMPC(I) + COEF(2,3) * X3CU(I)
C       DF(I,2,1) = COEF(2,1) * X2SQ(I)  
C       DF(I,2,2) = 2. * TEMPB(I)
C       DF(I,2,3) = COEF(2,2) * X2SQ(I) + COEF(2,3) * 3. * X3SQ(I)      
C  3  CONTINUE
C******************  END OF BLOCK A  *******************
C
C ON INPUT:
C
C N  IS THE NUMBER OF EQUATIONS AND VARIABLES.
C
C NUMT(J)  IS THE NUMBER OF TERMS IN THE JTH EQUATION.
C
C MAXT  IS AN UPPER BOUND ON NUMT(J) FOR J=1 TO N.
C
C KDEG(J,L,K)  IS THE DEGREE OF THE L-TH VARIABLE IN THE K-TH TERM
C   OF THE J-TH EQUATION.
C
C COEF(J,K)  IS THE K-TH COEFFICIENT OF THE J-TH EQUATION.
C
C CL  IS USED TO DEFINE THE PROJECTIVE TRANSFORMATION.  IF 
C   THE PROJECTIVE TRANSFORMATION IS NOT SPECIFIED, THEN  CL
C   CONTAINS DUMMY VALUES.
C
C X(1,J), X(2,J)  ARE THE REAL AND IMAGINARY PARTS RESPECTIVELY OF
C   THE J-TH INDEPENDENT VARIABLE.
C
C XX, TRM, DTRM, CLX, DXNP1  ARE WORKSPACE VARIABLES.  
C
C ON OUTPUT:
C
C F(1,J), F(2,J)  ARE THE REAL AND IMAGINARY PARTS RESPECTIVELY OF 
C   THE J-TH EQUATION.
C
C DF(1,J,K), DF(2,J,K)  ARE THE REAL AND IMAGINARY PARTS RESPECTIVELY
C   OF THE K-TH PARTIAL DERIVATIVE OF THE J-TH EQUATION.
C
C
C VARIABLES: XNP1,TEMP1,TEMP2.
C 
C NOTE:  XNP1(1), XNP1(2)  ARE THE REAL AND IMAGINARY PARTS, 
C   RESPECTIVELY, OF THE PROJECTIVE VARIABLE.  XNP1  IS UNITY 
C   IF THE PROJECTIVE TRANSFORMATION IS NOT SPECIFIED.
C
C  SUBROUTINES: MULP,POWP,DIVP.
C
      USE REAL_PRECISION
C
      INTERFACE
        SUBROUTINE DIVP(XXXX,YYYY,ZZZZ,IERR)
        USE REAL_PRECISION
        REAL (KIND=R8), DIMENSION(2), INTENT(IN):: XXXX,YYYY
        INTEGER, INTENT(OUT):: IERR
        REAL (KIND=R8), DIMENSION(2), INTENT(OUT):: ZZZZ
        END SUBROUTINE DIVP
        SUBROUTINE MULP(XXXX,YYYY,ZZZZ)
        USE REAL_PRECISION
        REAL (KIND=R8), DIMENSION(2), INTENT(IN):: XXXX,YYYY
        REAL (KIND=R8), DIMENSION(2), INTENT(OUT):: ZZZZ
        END SUBROUTINE MULP
        SUBROUTINE POWP(NNNN,XXXX,YYYY)
        USE REAL_PRECISION
        INTEGER, INTENT(IN):: NNNN
        REAL (KIND=R8), DIMENSION(2), INTENT(IN):: XXXX
        REAL (KIND=R8), DIMENSION(2), INTENT(IN OUT):: YYYY
        END SUBROUTINE POWP
      END INTERFACE
C
C DECLARATION OF INPUT AND OUTPUT:
      INTEGER, INTENT(IN):: N,NUMT(N),MAXT,KDEG(N,N+1,MAXT)
      REAL (KIND=R8), INTENT(IN):: COEF(N,MAXT),CL(2,N+1),X(2,N)
      REAL (KIND=R8), INTENT(IN OUT):: 
     &  XX(2,N,N+1,MAXT),TRM(2,N,MAXT),DTRM(2,N,N+1,MAXT)
      REAL (KIND=R8), INTENT(OUT):: 
     &  CLX(2,N),DXNP1(2,N),F(2,N),DF(2,N,N+1)
C
C DECLARATION OF LOCAL VARIABLES:
      INTEGER:: IERR,J,K,L,M,NNNN,NP1
      REAL (KIND=R8), DIMENSION(2):: TEMP1,TEMP2,XNP1
C
      NP1=N+1
C
C GENERATE XNP1, THE PROJECTIVE COORDINATE, AND ITS DERIVATIVES.
      DO J=1,N
        CALL MULP(CL(1,J),X(1,J),CLX(1,J))
      END DO
C
      XNP1(1:2)=CL(1:2,NP1) + SUM(CLX(1:2,1:N),DIM=2)
      DXNP1(1:2,1:N)=CL(1:2,1:N)
C
C******************  BEGIN BLOCK A  *******************
C
C "BLOCK A" TAKES  X  AND  XNP1  AS INPUT AND RETURNS  F 
C AND  DF  AS OUTPUT.   F  IS THE HOMOGENEOUS FORM OF THE
C ORIGINAL  F, AND  DF  CONSISTS OF THE PARTIAL 
C DERIVATIVES OF THE HOMOGENEOUS FORM OF  F  WITH RESPECT 
C TO THE N+1 VARIABLES X(1), ... ,X(N), XNP1.
C
C BEGIN "COMPUTE F"
C
      DO J=1,N
        DO K=1,NUMT(J)
          CALL POWP(KDEG(J,NP1,K),XNP1, XX(1,J,NP1,K))
          DO L=1,N
            CALL POWP(KDEG(J, L,K),X(1,L),XX(1,J,  L,K))
          END DO
        END DO
      END DO
      TRM = 0.0
      DO J=1,N
        DO K=1,NUMT(J)
          TRM(1,J,K)=COEF(J,K)
          DO L=1,NP1
            CALL MULP(XX(1,J,L,K), TRM(1,J,K),TEMP1)
            TRM(1:2,J,K ) = TEMP1(1:2)
          END DO
        END DO
      END DO
      F(1:2,1:N) = SUM(TRM(1:2,1:N,:),DIM=3)
C
C END OF "COMPUTE F"
C
C BEGIN "COMPUTE DF"
C
      J_LOOP: DO J=1,N
        K_LOOP: DO K=1,NUMT(J)
        M_LOOP: DO M=1,NP1
C
C IF TERM DOES NOT INCLUDE X(M), SET PARTIAL DERIVATIVE OF TERM
C   EQUAL TO ZERO.
          IF(KDEG(J,M,K) .EQ. 0) THEN
            DTRM(1:2,J,M,K) = 0.0
          ELSE
C
C IF TERM DOES INCLUDE X(M), TRY COMPUTING THE PARTIAL BY DIVIDING
C   THE TERM BY X(M).
            IF(M.LE.N) CALL DIVP(TRM(1,J,K),X(1,M),DTRM(1,J,M,K),IERR)
            IF(M.EQ.NP1) CALL DIVP(TRM(1,J,K),XNP1,DTRM(1,J,M,K),IERR)
            IF (IERR .EQ. 0) THEN
              DTRM(1:2,J,M,K)=KDEG(J,M,K)*DTRM(1:2,J,M,K)
            ELSE
C
C IF DIVISION WOULD CAUSE OVERFLOW, GENERATE THE PARTIAL BY
C   THE POLYNOMIAL FORMULA.
              DTRM(1,J,M,K)=COEF(J,K)
              DTRM(2,J,M,K)=0.0
              DO L=1,NP1
                IF (L .EQ. M) CYCLE
                CALL MULP(XX(1,J,L,K),DTRM(1,J,M,K),TEMP1)
                DTRM(1:2,J,M,K)=TEMP1(1:2)
              END DO
              NNNN=KDEG(J,M,K)-1
              IF (M .LE. N) CALL POWP(NNNN,X(1,M),TEMP2)
              IF (M .EQ. NP1) CALL POWP(NNNN,XNP1 ,TEMP2)
              CALL MULP(TEMP2,TEMP1,DTRM(1,J,M,K))
              DTRM(1:2,J,M,K)=KDEG(J,M,K)*DTRM(1:2,J,M,K)
            END IF
          END IF
        END DO M_LOOP
        END DO K_LOOP
      END DO J_LOOP
      DO J=1,N
        DF(1:2,J,1:NP1) = SUM(DTRM(1:2,J,1:NP1,1:NUMT(J)), DIM=3)
      END DO
C
C END OF "COMPUTE DF"
C*******************  END BLOCK A  ********************
C
C CONVERT  DF  TO BE PARTIALS WITH RESPECT TO  X(1), ... ,X(N),
C BY APPLYING THE CHAIN RULE WITH  XNP1  CONSIDERED A FUNCTION OF 
C OF  X(1), ... ,X(N).
C
      DO J=1,N
        DO K=1,N
          CALL MULP(DF(1,J,NP1),DXNP1(1,K),TEMP1)
          DF(1:2,J,K)=DF(1:2,J,K)+TEMP1(1:2)
        END DO
      END DO
      RETURN
      END SUBROUTINE FFUNP
      SUBROUTINE FODE(S,Y,YP,YPOLD,A,QR,ALPHA,TZ,PIVOT,NFE,N,IFLAG)
C
C SUBROUTINE  FODE  IS USED BY SUBROUTINE  STEPS  TO SPECIFY THE
C ORDINARY DIFFERENTIAL EQUATION  DY/DS = G(S,Y) , WHOSE SOLUTION
C IS THE ZERO CURVE OF THE HOMOTOPY MAP.  S = ARC LENGTH,
C YP = DY/DS, AND  Y(S) = (LAMBDA(S), X(S)) .
C
C CALLS  DGEQPF , DNRM2 .
C
      USE HOMOTOPY
      USE REAL_PRECISION
      REAL (KIND=R8):: DNRM2,S,YPNORM
      INTEGER:: I,IFLAG,IK,K,KP1,LW,N,NFE,NP1
C
C *****  ARRAY DECLARATIONS.  *****
C
      REAL (KIND=R8):: A(:),Y(:),YP(N+1),YPOLD(N+1)
C
C ARRAYS FOR COMPUTING THE JACOBIAN MATRIX AND ITS KERNEL.
      REAL (KIND=R8):: ALPHA(3*N+3),QR(N,N+1),TZ(N+1)
      INTEGER, DIMENSION(N+1):: PIVOT
C
C *****  END OF DIMENSIONAL INFORMATION.  *****
C
C
      NP1=N+1
      NFE=NFE+1
C NFE CONTAINS THE NUMBER OF JACOBIAN EVALUATIONS.
C   *   *   *   *   *   *   *   *   *   *   *   *   *   *   *   *   *
C
C COMPUTE THE JACOBIAN MATRIX, STORE IT IN QR.
C
      IF (IFLAG .EQ. -2) THEN
C
C  QR = ( D RHO(A,LAMBDA,X)/D LAMBDA , D RHO(A,LAMBDA,X)/DX )  .
C
        DO K=1,NP1
          CALL RHOJAC(A,Y(1),Y(2:NP1),QR(:,K),K)
        END DO
      ELSE
        CALL F(Y(2:NP1),TZ(1:N))
        QR(:,1)=A
        IF (IFLAG .EQ. 0) THEN
C
C      QR = ( A - F(X), I - LAMBDA*DF(X) )  .
C
          QR(:,1)=QR(:,1)-TZ(1:N)
          DO K=1,N
            CALL FJAC(Y(2:NP1),TZ(1:N),K)
            KP1=K+1
            QR(:,KP1)=-Y(1)*TZ(1:N)
            QR(K,KP1)=1.0+QR(K,KP1)
          END DO
        ELSE
C
C   QR = ( F(X) - X + A, LAMBDA*DF(X) + (1 - LAMBDA)*I ) .
C
          QR(:,1)=QR(:,1)+TZ(1:N)-Y(2:NP1)
          DO K=1,N
            CALL FJAC(Y(2:NP1),TZ(1:N),K)
            KP1=K+1
            QR(:,KP1)=Y(1)*TZ(1:N)
            QR(K,KP1)=1.0-Y(1)+QR(K,KP1)
          END DO
        ENDIF
      ENDIF
C
C   *   *   *   *   *   *   *   *   *   *   *   *   *   *   *   *   *
C REDUCE THE JACOBIAN MATRIX TO UPPER TRIANGULAR FORM.
C
        PIVOT = 0
C
      CALL DGEQPF(N,NP1,QR,N,PIVOT,YP,ALPHA,K)
C
      IF (ABS(QR(N,N)) .LE. ABS(QR(1,1))*EPSILON(1.0_R8)) THEN 
        IFLAG=4
        RETURN
      ENDIF 
C COMPUTE KERNEL OF JACOBIAN, WHICH SPECIFIES YP=DY/DS.
      TZ(NP1)=1.0
      DO LW=1,N
        I=NP1-LW
        IK=I+1
        TZ(I)=-DOT_PRODUCT(QR(I,IK:NP1),TZ(IK:NP1))/QR(I,I)
      END DO
      YPNORM=DNRM2(NP1,TZ,1)
      YP(PIVOT)=TZ/YPNORM
      IF (DOT_PRODUCT(YP,YPOLD) .LT. 0.0) YP=-YP
C
C SAVE CURRENT DERIVATIVE (= TANGENT VECTOR) IN  YPOLD .
      YPOLD=YP
      RETURN
      END SUBROUTINE FODE
      SUBROUTINE FODEDS(S,Y,YP,N,IFLAG,YPOLD,A,NDIMA,LENQR,MODE,NFE)
C
C SUBROUTINE  FODEDS  IS USED BY SUBROUTINE  STEPDS  TO SPECIFY THE
C ORDINARY DIFFERENTIAL EQUATION  DY/DS = G(S,Y) , WHOSE SOLUTION
C IS THE ZERO CURVE OF THE HOMOTOPY MAP.  S = ARC LENGTH,
C YP = DY/DS, AND  Y(S) = (X(S), LAMBDA(S)) .
C
      USE HOMOTOPY, QR => QRSPARSE
      USE REAL_PRECISION
      REAL (KIND=R8):: LAMBDA,S,YPNORM
      INTEGER:: IFLAG,J,JPOS,LENQR,MODE,N,NDIMA,NFE,NP1
      REAL (KIND=R8):: A(NDIMA),Y(N+1),YP(N+1),YPOLD(N+1)
      INTERFACE
        SUBROUTINE PCGDS(N,LENQR,IFLAG,YP,RHO)
          USE REAL_PRECISION
          INTEGER, INTENT(IN):: LENQR,N
          INTEGER, INTENT(IN OUT):: IFLAG
          REAL (KIND=R8), INTENT(IN OUT):: YP(N+1)
          REAL (KIND=R8), OPTIONAL, INTENT(IN):: RHO(N)
        END SUBROUTINE PCGDS
        SUBROUTINE GMRILUDS(N,LENQR,IFLAG,YP,RHO)
          USE REAL_PRECISION
          INTEGER, INTENT(IN):: LENQR,N
          INTEGER, INTENT(IN OUT):: IFLAG
          REAL (KIND=R8), INTENT(IN OUT):: YP(N+1)
          REAL (KIND=R8), OPTIONAL, INTENT(IN):: RHO(N)
        END SUBROUTINE GMRILUDS
        FUNCTION DNRM2(N,X,STRIDE)
          USE REAL_PRECISION
          INTEGER:: N,STRIDE
          REAL (KIND=R8):: DNRM2,X(N)
        END FUNCTION DNRM2
      END INTERFACE
C
C *****  END OF SPECIFICATION INFORMATION.  *****
C
      NP1=N+1
      NFE=NFE+1
C NFE CONTAINS THE NUMBER OF JACOBIAN EVALUATIONS.
      LAMBDA=Y(NP1)
      ROWPOS(NP1)=LENQR+1
C   *   *   *   *   *   *   *   *   *   *   *   *   *   *   *   *   *
C MODE = 1 STORAGE FORMAT.
C   *   *   *   *   *   *   *   *   *   *   *   *   *   *   *   *   *
C
      IF (MODE .EQ. 1) THEN
C COMPUTE THE JACOBIAN MATRIX, STORE IT IN  [QR | -PP] .
C
      IF (IFLAG .EQ. -2) THEN
C
C  [QR | -PP] = [ D RHO(A,X,LAMBDA)/DX | D RHO(A,X,LAMBDA)/D LAMBDA ]  .
C
C  PP = - (D RHO(A,X,LAMBDA)/D LAMBDA) .
        CALL RHOJS(A,LAMBDA,Y(1:N))
C
      ELSE
        CALL F(Y(1:N),PP)
        IF (IFLAG .EQ. 0) THEN
C
C      [QR | -PP] = [ I - LAMBDA*DF(X) | A - F(X) ]  .
C
          PP = PP - A(1:N)
          CALL FJACS(Y(1:N))
          QR = (-LAMBDA)*QR
          QR(ROWPOS(1:N)) = QR(ROWPOS(1:N)) + 1.0
        ELSE
C
C   [QR | -PP] = [ LAMBDA*DF(X) + (1 - LAMBDA)*I | F(X) - X + A ] .
C
          PP = Y(1:N) - A(1:N) - PP
          CALL FJACS(Y(1:N))
          QR = LAMBDA*QR
          QR(ROWPOS(1:N)) = QR(ROWPOS(1:N)) + 1.0 - LAMBDA
        ENDIF
      ENDIF
      ELSE
C   *   *   *   *   *   *   *   *   *   *   *   *   *   *   *   *   *
C MODE = 2 STORAGE FORMAT.
C   *   *   *   *   *   *   *   *   *   *   *   *   *   *   *   *   *
C
      IF (IFLAG .EQ. -2) THEN
C
C  [QR] = [ D RHO(A,X,LAMBDA)/DX , D RHO(A,X,LAMBDA)/D LAMBDA ]  .
C
        CALL RHOJS(A,LAMBDA,Y(1:N))
C
      ELSE
        CALL F(Y(1:N),PP)
        IF (IFLAG .EQ. 0) THEN
C
C      [QR | -PP] = [ I - LAMBDA*DF(X) | A - F(X) ]  .
C
          PP = PP - A(1:N)
          CALL FJACS(Y(1:N))
          QR = (-LAMBDA)*QR
C FIND INDEX JPOS OF DIAGONAL ELEMENT IN JTH ROW OF QR.
          DO J=1,N
            JPOS=ROWPOS(J)
            DO
              IF (COLPOS(JPOS) .EQ. J) EXIT
              JPOS=JPOS+1
              IF (JPOS < ROWPOS(J+1)) CYCLE
              IFLAG=4
              RETURN
            END DO
            QR(JPOS) = QR(JPOS) + 1.0
          END DO
        ELSE
C
C   [QR | -PP] = [ LAMBDA*DF(X) + (1 - LAMBDA)*I | F(X) - X + A ] .
C
          PP = Y(1:N) - A(1:N) - PP
          CALL FJACS(Y(1:N))
          QR = LAMBDA*QR
C FIND INDEX JPOS OF DIAGONAL ELEMENT IN JTH ROW OF QR.
          DO J=1,N
            JPOS=ROWPOS(J)
            DO
              IF (COLPOS(JPOS) .EQ. J) EXIT
              JPOS=JPOS+1
              IF (JPOS < ROWPOS(J+1)) CYCLE
              IFLAG=4
              RETURN
            END DO
            QR(JPOS) = QR(JPOS) + 1.0 - LAMBDA
          END DO
        ENDIF
      ENDIF
      ENDIF
C
C   *   *   *   *   *   *   *   *   *   *   *   *   *   *   *   *   *
      YP=YPOLD
C COMPUTE KERNEL OF JACOBIAN, WHICH SPECIFIES YP=DY/DS, USING A
C PRECONDITIONED CONJUGATE GRADIENT ALGORITHM.
      SELECT CASE (MODE)
        CASE (1)
        CALL PCGDS(N,LENQR,IFLAG,YP)
        CASE (2)
        CALL GMRILUDS(N,LENQR,IFLAG,YP)
      END SELECT
      IF (IFLAG .GT. 0) RETURN
C
C NORMALIZE TANGENT VECTOR YP.
      YPNORM=DNRM2(NP1,YP,1)
      YP = (1.0/YPNORM)*YP
C
C CHOOSE UNIT TANGENT VECTOR DIRECTION TO MAINTAIN CONTINUITY.
      IF (DOT_PRODUCT(YP,YPOLD) .LT. 0.0) YP = -YP
C
C SAVE CURRENT DERIVATIVE (= TANGENT VECTOR) IN  YPOLD .
      YPOLD = YP
C
      RETURN
      END SUBROUTINE FODEDS
      SUBROUTINE GFUNP(N,IDEG,PDG,QDG,X,XDGM1,XDG,PXDGM1,PXDG,G,DG)
C
C GFUNP  EVALUATES THE START EQUATION "G".
C
C ON INPUT:
C
C N  IS THE NUMBER OF VARIABLES.
C
C IDEG(J)  IS THE DEGREE OF THE J-TH EQUATION.
C
C PDG(1,J), PDG(2,J)  ARE THE REAL AND IMAGINARY PARTS
C   OF THE POWERS OF P USED TO DEFINE G.
C
C QDG(1,J), QDG(2,J)  ARE THE REAL AND IMAGINARY PARTS
C   OF THE POWERS OF Q USED TO DEFINE G.
C
C X(1,J), X(2,J)  ARE THE REAL AND IMAGINARY PARTS OF THE
C   J-TH INDEPENDENT VARIABLE.
C
C XDGM1,XDG,PXDGM1,PXDG ARE WORKSPACE ARRAYS.
C
C ON OUTPUT:
C
C N,IDEG,PDG,QDG, AND X  ARE UNCHANGED. 
C
C G(1,J),G(2,J)  ARE THE REAL AND IMAGINARY PARTS OF THE
C   J-TH START EQUATION.
C
C DG(1,J),DG(2,J)  ARE THE REAL AND IMAGINARY PARTS OF THE
C   PARTIAL DERIVATIVES OF THE J-TH START EQUATION WITH RESPECT TO THE
C   J-TH INDEPENDENT VARIABLE.
C
C SUBROUTINES:  MULP, POWP.
C
      USE REAL_PRECISION
C
      INTERFACE
        SUBROUTINE DIVP(XXXX,YYYY,ZZZZ,IERR)
        USE REAL_PRECISION
        REAL (KIND=R8), DIMENSION(2), INTENT(IN):: XXXX,YYYY
        INTEGER, INTENT(OUT):: IERR
        REAL (KIND=R8), DIMENSION(2), INTENT(OUT):: ZZZZ
        END SUBROUTINE DIVP
        SUBROUTINE MULP(XXXX,YYYY,ZZZZ)
        USE REAL_PRECISION
        REAL (KIND=R8), DIMENSION(2), INTENT(IN):: XXXX,YYYY
        REAL (KIND=R8), DIMENSION(2), INTENT(OUT):: ZZZZ
        END SUBROUTINE MULP
        SUBROUTINE POWP(NNNN,XXXX,YYYY)
        USE REAL_PRECISION
        INTEGER, INTENT(IN):: NNNN
        REAL (KIND=R8), DIMENSION(2), INTENT(IN):: XXXX
        REAL (KIND=R8), DIMENSION(2), INTENT(IN OUT):: YYYY
        END SUBROUTINE POWP
      END INTERFACE
C
C DECLARATION OF INPUT AND OUTPUT:
      INTEGER, INTENT(IN):: N,IDEG(N)
      REAL (KIND=R8), INTENT(IN):: PDG(2,N),QDG(2,N),X(2,N)
      REAL (KIND=R8), INTENT(IN OUT):: XDGM1(2,N),XDG(2,N),PXDGM1(2,N),
     &  PXDG(2,N)
      REAL (KIND=R8), INTENT(OUT):: G(2,N),DG(2,N)
C
C DECLARATION LOCAL OF VARIABLES
      INTEGER:: I,J
C
C COMPUTE THE (IDEG-1)-TH AND IDEG-TH POWER OF X
      DO J=1,N
        CALL POWP(IDEG(J)-1,X(1,J), XDGM1(1,J))
        CALL MULP(X(1,J),XDGM1(1,J), XDG(1,J))
      END DO
C
C COMPUTE THE PRODUCT OF PDG AND XDGM1
      DO J=1,N
          CALL MULP( PDG(1,J), XDGM1(1,J), PXDGM1(1,J) )
      END DO
C
C COMPUTE THE PRODUCT OF PDG AND XDG
      DO J=1,N
          CALL MULP( PDG(1,J), XDG(1,J), PXDG(1,J) )
      END DO
      G = PXDG - QDG
      DO J=1,N
        DG(1:2,J)= IDEG(J)*PXDGM1(1:2,J)
      END DO
      RETURN
      END SUBROUTINE GFUNP
      SUBROUTINE GMFADS(NN,A,NWK,MAXA)
C
C     This subroutine computes the LDU decomposition of a symmetric positive
C     definite matrix B where only the upper triangular skyline structure
C     is stored.  The decomposition is done by the Gill-Murray
C     strategy from P.E. Gill and W. Murray, Newton type Methods
C     for Unconstrained and Linearly Constrained Optimization,
C     Mathematical Programming, 7, 311-350 (1974) and gives an
C     approximate decomposition in the case of a nonpositive
C     definite or ill-conditioned matrix.
C
C     Input variables:
C
C        NN -- dimension of B.
C
C        A -- one dimensional real array containing the upper 
C             triangular skyline portion of a symmetric matrix B in
C             packed skyline storage format.
C
C        NWK -- number of elements in A.
C
C        MAXA -- an integer array of dimension NN+1 containing the 
C                locations of the diagonal elements of B in A.  
C                By convention, MAXA(NN+1)=NWK+1.  
C
C     Output variables:
C
C        A -- the upper triangular skyline portion of the LDU 
C             decomposition of the symmetric matrix B (or B + E if B
C             was not sufficiently positive definite).
C
C
C     No working storage is required by this routine.
C
C
      USE REAL_PRECISION
      INTEGER, INTENT(IN):: NN,NWK,MAXA(NN+1)
      REAL (KIND=R8), INTENT(IN OUT):: A(NWK)
      INTEGER:: I,I0,I1,I2,I3,I4,J,J1,K,K1,K2,KH,KL,KN,KU,KZ,L,L1,
     &   L2,L3,M,M1,N1,NNN
      REAL (KIND=R8):: BET,DEL,DJ,G,GAM,GAM1,PHI,
     &   THE,THE1,XT1,XT2,ZET,ZET1
      G=0.0
      GAM=0.0
      DO I=1,NN
         K=MAXA(I)
         G=G+A(K)*A(K)
         GAM1=ABS(A(K))
         IF(GAM1.GT.GAM)GAM=GAM1
      END DO
      ZET=0.0
      DO I=1,NN
         K=MAXA(I)
         K1=MAXA(I+1)-1
         K2=K1-K
         IF (K2.EQ.0) CYCLE
         L=K+1
         DO J=L,K1
            G=G+2.0*A(J)*A(J)
            ZET1=ABS(A(J))
            IF(ZET1.GT.ZET)ZET=ZET1
         END DO
      END DO
      ZET=ZET/NN
      DEL=EPSILON(1.0_R8)
      BET=DEL
      IF (ZET .GT. BET) BET=ZET
      IF (GAM .GT. BET) BET=GAM
      G=SQRT(G)
      IF (G .GT. 1.0) DEL=DEL*G
      DO I=1,NN
         N1=I-1
         KN=MAXA(I)
         KL=KN+1
         KU=MAXA(I+1)-1
         KH=KU-KL
         PHI=A(KN)
         IF (KH .GE. 0) THEN
           K1=KN+1
           K2=I
           DO J=K1,KU
              K2=K2-1
              KZ=MAXA(K2)
              PHI=PHI-A(J)*A(J)*A(KZ)
           END DO
         END IF
         PHI=ABS(PHI)
         L=I+1
         THE=0.0
         NNN=NN+1
         IF (L .NE. NNN) THEN
           DO J=L,NN
              L1=MAXA(J)
              L2=MAXA(J+1)
              L3=L2-L1-1
              M=J-I
              IF (L3 .LT. M) CYCLE
              M1=L1+M
              IF (N1 .NE. 0) THEN
                DO J1=1,N1
                  I0=MAXA(J1)
                  I1=MAXA(L)
                  I2=I-J1
                  I3=I1-KN-1
                  I4=J-J1
                  IF (I3 .LT. I2) CYCLE
                  IF (L3 .LT. I4) CYCLE
                  XT1=A(KN+I2)
                  XT2=A(L1+I4)
                  A(M1)=A(M1)-XT1*XT2*A(I0)
                END DO
              END IF
            THE1=ABS(A(M1))
            IF (THE .LT. THE1) THE=THE1
            END DO
         END IF
         THE=THE*THE/BET
         DJ=DEL
         IF (PHI .GT. DJ) DJ=PHI
         IF (THE .GT. DJ) DJ=THE
         A(KN)=DJ
         IF (L .EQ. NNN) CYCLE
         DO J=L,NN
            L1=MAXA(J)
            L2=MAXA(J+1)
            L3=L2-L1-1
            M=J-I
            IF (L3 .LT. M) CYCLE
            M1=L1+M
            A(M1)=A(M1)/A(KN)
         END DO
      END DO
      RETURN
      END SUBROUTINE GMFADS
      SUBROUTINE GMRES(N, KDMAX, ITMAX, RHS, X, KVAL,
     &                PRECON, IFLAG, ROWPOSP, COLPOSP)
C
C THIS ROUTINE IS AN EXTENSION OF THE STANDARD RESTARTED GMRES METHOD OF
C Y. SAAD AND M. SCHULTZ, "GMRES: A GENERALIZED MINIMAL RESIDUAL METHOD
C FOR SOLVING NONSYMMETRIC LINEAR SYSTEMS", SIAM J. SCI. STAT. COMPUT., 7
C (1986), PP. 856-869, THAT ALLOWS THE MAXIMUM KRYLOV SUBSPACE DIMENSION
C TO INCREASE (UP TO A MAXIMUM SPECIFIED PARAMETER VALUE) IF STAGNATION
C OCCURS.  THE ARNOLDI BASIS VECTORS ARE GENERATED USING ORTHOGONALIZATION
C WITH HOUSEHOLDER TRANSFORMATIONS.
C ON RESTARTS, RESIDUAL VECTORS ARE COMPUTED BY DIRECT EVALUATION.
C CONDITIONING OF THE GMRES LEAST-SQUARES PROBLEM IS MONITORED USING 
C LAPACK ROUTINE DLAIC1 (SEE P. N. BROWN AND H. F. WALKER, "GMRES ON 
C (NEARLY) SINGULAR SYSTEMS", UTAH STATE UNIV.  C MATH. STAT. DEPT. RES.
C REPORT 2/94/73, FEBRUARY, 1994).
C
C
C INPUT VARIABLES:
C
C N  IS THE DIMENSION OF THE LINEAR SYSTEM AA*X = RHS.  THE SPARSE
C   MATRIX DATA STRUCTURE FOR AA IS STORED IN THE MODULE HOMOTOPY.
C
C KDMAX  IS THE MAXIMUM KRYLOV SUBSPACE DIMENSION ALLOWED BEFORE
C   STAGNATION DETECTION BEGINS.
C
C ITMAX  IS THE MAXIMUM OVERALL NUMBER OF GMRES ITERATIONS ALLOWED.
C
C RHS  IS THE RIGHT HAND SIDE OF THE LINEAR SYSTEM.
C
C X  IS THE INITIAL APPROXIMATE SOLUTION OF THE LINEAR SYSTEM.
C
C KVAL  IS AN INDEX INTO X REQUIRED FOR MATRIX-VECTOR MULTIPLICATION.
C  
C PRECON  IS A ONE-DIMENSIONAL REAL ARRAY CONTAINING A PRECONDITIONING
C   MATRIX Q FOR AA.   
C
C IFLAG  INDICATES WHAT DATA IS USED IN MATRIX-VECTOR MULTIPLICATION.
C
C ROWPOSP  IS AN OPTIONAL INTEGER ARRAY USED IN THE DATA STRUCTURE
C   DESCRIBING THE SPARSE PRECONDITIONING MATRIX Q. 
C
C COLPOSP  IS AN OPTIONAL INTEGER ARRAY USED IN THE DATA STRUCTURE
C   DESCRIBING THE SPARSE PRECONDITIONING MATRIX Q. 
C
C OUTPUT VARIABLES:
C
C KDMAX  IS THE INCREASED MAXIMUM KRYLOV SUBSPACE DIMENSION ON
C   TERMINATION.
C
C IFLAG  IS UNCHANGED ON NORMAL TERMINATION (DESIRED TOLERANCE MET);
C   = 4, IF GMRES ITERATION FAILS TO CONVERGE BECAUSE
C      - DESIRED TOLERANCE NOT MET AFTER ITMAX ITERATIONS;
C      - DANGEROUS NEAR-SINGULARITY DETECTED;
C      - LIMITS OF NUMERICAL ACCURACY REACHED;
C      - STAGNATION OCCURS.
C
C ITMAX  IS THE TOTAL NUMBER OF GMRES ITERATIONS PERFORMED.
C
C X  IS THE FINAL APPROXIMATE SOLUTION.
C   
C
C  CALLS DNRM2, ILUSOLVDS, MULTDS, MULT2DS, SOLVDS, AND INTERNAL
C  SUBROUTINES MULM1 AND MULM2 AND INTERNAL FUNCTIONS APPL_HOUSE, HOUSE.
C
      USE HOMOTOPY, AA=>QRSPARSE
      USE REAL_PRECISION
      INTEGER, INTENT(IN):: KVAL, N
      INTEGER, INTENT(IN OUT):: IFLAG, ITMAX, KDMAX
      REAL (KIND=R8), INTENT(IN):: PRECON(:), RHS(N)  
      REAL (KIND=R8), INTENT(IN OUT):: X(N)
      INTEGER, INTENT(IN), OPTIONAL:: COLPOSP(:), ROWPOSP(:)
C
C IPRINT  IS A PARAMETER FOR THE FORTRAN UNIT NUMBER, WHICH IF POSITIVE
C   ACTIVATES TRACING OF THE GMRES ITERATIONS.
C M  IS A PARAMETER USED IN THE ADAPTIVE STRATEGY FOR INCREMENTING THE 
C   KRYLOV SUBSPACE DIMENSION (VALUE KDMAX).  
C STRONG_VERSION  IS A LOGICAL PARAMETER, WHICH IF TRUE CAUSES GMRES TO 
C   STOP IN TWO CASES: WHEN THE ESTIMATED CONDITION NUMBER IS GREATER 
C   THAN 1.0/(50*EPSILON); 
C   WHEN THE TRUE RESIDUAL NORM IS LARGE AND IT INCREASED DURING THE LAST 
C   RESTART CYCLE.
C
      INTEGER, PARAMETER:: IPRINT=0, M=4 
      LOGICAL, PARAMETER:: STRONG_VERSION=.TRUE.
C
C LOCAL VARIABLES.
C
      INTEGER:: I, IFLAGC, IFLAGI, IQUIT, ITNO, KD, 
     &  KDP1, KDLIMIT, LENAA 
      REAL (KIND=R8), ALLOCATABLE:: C(:), R(:,:), S(:), SVBIG(:), 
     &  SVSML(:), V(:,:), W(:)
      REAL (KIND=R8):: VTEMP(N), VTEMP2(N)
      REAL (KIND=R8):: BIG, BIGCND, CC, CNDMAX,
     &  RSN, RSNOLD, SESTPR, SMALL, SS, TEMP, TOL 
      INTERFACE
        FUNCTION DNRM2(N,X,STRIDE)
        USE REAL_PRECISION
        INTEGER:: N,STRIDE
        REAL (KIND=R8):: DNRM2,X(N)
        END FUNCTION DNRM2
      END INTERFACE
C
C DEFINE MAXIMUM ALLOWED KRYLOV SUBSPACE DIMENSION.
C
      KDLIMIT = MIN(MAX(200, INT(SQRT(REAL(N)))), (N+1)/2)
      KDMAX = MIN(KDMAX,KDLIMIT)
C
C ALLOCATE LOCAL ARRAYS.
C
      IF(.NOT. ALLOCATED(C)) ALLOCATE(C(KDLIMIT))
      IF(.NOT. ALLOCATED(R)) ALLOCATE(R(KDLIMIT,KDLIMIT))
      IF(.NOT. ALLOCATED(S)) ALLOCATE(S(KDLIMIT))
      IF(.NOT. ALLOCATED(SVBIG)) ALLOCATE(SVBIG(KDLIMIT))
      IF(.NOT. ALLOCATED(SVSML)) ALLOCATE(SVSML(KDLIMIT))
      IF(.NOT. ALLOCATED(V)) ALLOCATE(V(N,KDLIMIT+1))
      IF(.NOT. ALLOCATED(W)) ALLOCATE(W(KDLIMIT+1))
C
C PERFORM INITIALIZATIONS.
C
      ITNO = 0
      IFLAGC = 0
      IFLAGI = 1
      LENAA = ROWPOS(N)-1
      TOL = MAX(100.0, 1.01*REAL(LENAA)/REAL(N))*EPSILON(1.0_R8)
      VTEMP = X
      IF (PRESENT(ROWPOSP) .AND. PRESENT( COLPOSP)) THEN
        CALL MULM2(V(:,1),VTEMP)            ! MODE=2 
      ELSE 
        CALL MULM1(V(:,1),VTEMP)            ! MODE=1 
      END IF
      V(:,1) = RHS - V(:,1)
      RSN = DNRM2(N, V(:,1), 1)
      TEMP = DNRM2(N, RHS, 1)
      TOL = MAX(RSN,TEMP)*TOL
      IF (RSN .LE. TOL) THEN 
         ITMAX = ITNO
         CALL CLEANUP
         RETURN
      ELSE
         IF (ITMAX .EQ. 0) THEN 
            ITMAX = ITNO
            IFLAG = 4
            CALL CLEANUP
            RETURN
         ENDIF
      ENDIF
      CNDMAX = 1.0/(50.0*EPSILON(1.0_R8))
      IQUIT = 0
      BIGCND = 0.0
C
C FOR PRINTING:
C
      IF (IPRINT .GT. 0) THEN 
         WRITE(IPRINT, 19) ITNO, RSN
19       FORMAT(//,9X,'GMRES ITERATIONS:',//,T4,'ITERATION',TR4,
     &   'RESIDUAL NORM',TR5,'CONDITION NUMBER',/,T4,I6,TR7,ES16.8)
      ENDIF
C
C BEGIN THE OUTER LOOP. 
C
      OUTER: DO
C
C FOR PRINTING:
C
      IF (IPRINT .GT. 0 ) THEN
        IF (IFLAGI .NE. 0) THEN
          WRITE (IPRINT, 21) KDMAX
21        FORMAT(/' RESTART WITH KRYLOV SUBSPACE DIMENSION =',I2)
        ELSE
           WRITE(IPRINT, 22) 
22         FORMAT(/' RESTART')
        ENDIF
      ENDIF
C
      KD = 0
      RSNOLD = RSN
      IQUIT = 0
      IFLAGI = 0
C
C FIND HOUSEHOLDER VECTOR V(:,1).
C
      W(1)=V(1,1)
      CALL HREFG(N,V(2:N,1),V(1,1),W(1))
C
C BEGIN THE INNER LOOP.
C
 200  INNER: DO
      KD = KD + 1
      KDP1 = KD + 1
      ITNO = ITNO + 1
      IQUIT = 0 
C
C FIND THE KD-TH ORTHOGONAL VECTOR.
C
      VTEMP = 0.0 
      VTEMP(KD) = 1.0
      DO I=KD,1,-1
        TEMP=V(I,I)
        V(I,I)=1.0
        CALL HREFX(N-I+1,V(I:N,I),TEMP,VTEMP(I:N))
        V(I,I)=TEMP
      END DO
C
C EVALUATE AA*(KD-TH KRYLOV SUBSPACE BASIS VECTOR). 
C
        IF (PRESENT(ROWPOSP) .AND. PRESENT(COLPOSP)) THEN  ! MODE = 2
          CALL ILUSOLVDS(N, PRECON(1:ROWPOSP(N+1)-1), 
     &      ROWPOSP(N+1)-1, ROWPOSP(1:N+1), 
     &      COLPOSP(1:ROWPOSP(N+1)-1), VTEMP(1:N))
          CALL MULM2(VTEMP2, VTEMP(1:N))
        ELSE                                               ! MODE = 1 
          CALL SOLVDS(N, PRECON, ROWPOS(N+1)-1, ROWPOS(1:N+1),
     &      VTEMP(1:N))
          CALL MULM1(VTEMP2, VTEMP(1:N))
        ENDIF
C
C FIND HOUSEHOLDER VECTOR V(:,KDP1).
C
      VTEMP(1:N)=VTEMP2(1:N)
      DO I=1,KD
        TEMP=V(I,I)
        V(I,I)=1.0
        CALL HREFX(N-I+1,V(I:N,I),TEMP,VTEMP(I:N))
        V(I,I)=TEMP
      END DO 
      IF ( MAXVAL(ABS(VTEMP(KDP1:N))) .NE. 0.0) THEN
        V(KDP1+1:N, KDP1) = VTEMP(KDP1+1:N)
        CALL HREFG(N-KDP1+1,V(KDP1+1:N, KDP1),V(KDP1,KDP1), 
     &   VTEMP(KDP1))
      ENDIF
C
C IF KD .GT. 1, APPLY THE PREVIOUS GIVENS ROTATIONS. 
C
      DO I = 1, KD-1
        TEMP = VTEMP(I) 
        VTEMP(I) = C(I)*TEMP + S(I)*VTEMP(I+1)
        VTEMP(I+1) = -S(I)*TEMP + C(I)*VTEMP(I+1)
      END DO   
C
C DETERMINE AND APPLY THE NEXT ROTATION. 
C
      IF (VTEMP(KDP1) .NE. 0.0) THEN
        TEMP = VTEMP(KD)
        R(KD,KD) = SQRT(VTEMP(KD)**2 + VTEMP(KDP1)**2) 
        C(KD) = TEMP/R(KD,KD)
        S(KD) = VTEMP(KDP1)/R(KD,KD)
        VTEMP(KD) = C(KD)*TEMP + S(KD)*VTEMP(KDP1)
      END IF 
      R(1:KD, KD) = VTEMP(1:KD)
C
C COMPUTE AND TEST INCREMENTAL CONDITION NUMBER.
C
      IF (KD .EQ. 1) THEN 
        BIG = R(KD,KD) 
        SMALL = BIG 
        SVBIG(1) = 1.0
        SVSML(1) = 1.0
      ELSE
        I = 1
        CALL DLAIC1(I, KD-1, SVBIG, BIG, R(1,KD), R(KD,KD),
     &    SESTPR, SS, CC)
        BIG = SESTPR
        SVBIG(1:KD-1) = SS*SVBIG(1:KD-1)
        SVBIG(KD) = CC
        I = 2
        CALL DLAIC1(I, KD-1, SVSML, SMALL, R(1,KD), R(KD,KD),
     &    SESTPR, SS, CC)
        SMALL = SESTPR
        SVSML(1:KD-1) = SS*SVSML(1:KD-1)
        SVSML(KD) = CC
      ENDIF
      IF (STRONG_VERSION) THEN
        IF (BIG .GE. SMALL*CNDMAX) THEN
          IF (IPRINT .GT. 0) THEN
            WRITE (IPRINT, 230) CNDMAX
 230        FORMAT(/,4X, 'IN GMRES CONDITION NUMBER .GE.',
     &      ES16.8)
          ENDIF
          ITNO = ITNO - 1
          IFLAGC = 4
          EXIT OUTER
        ENDIF
      ENDIF  
      TEMP = BIG/SMALL
      IF (TEMP .GT. BIGCND) BIGCND = TEMP
C
C UPDATE W AND THE RESIDUAL NORM. 
C
      IF (VTEMP(KDP1) .NE. 0.0) THEN
        W(KDP1)= -S(KD)*W(KD)
        W(KD) = C(KD)*W(KD)
      ELSE
        W(KDP1) = 0.0
      ENDIF 
      RSN = ABS(W(KDP1))
C
C FOR PRINTING:
C
      IF (IPRINT .GT. 0 ) THEN
        WRITE (IPRINT, 240) ITNO, RSN, BIG/SMALL                             
 240    FORMAT(T4,I6, TR7,2ES16.8)
      ENDIF
C
C TEST FOR TERMINATION OF THE INNER LOOP. 
C
      IF (RSN .LE. TOL .OR. KD .EQ. KDMAX .OR. ITNO .EQ. ITMAX)
     &  EXIT INNER
      END DO INNER
C
C IF TERMINATING THE INNER LOOP, TEST FOR TERMINATION OF THE OUTER LOOP, 
C COMPUTE THE CORRECTION, AND UPDATE THE APPROXIMATE SOLUTION. 
C
C TEST FOR TERMINATION OF THE OUTER LOOP. 
C
      IF (RSN .LE. TOL .OR. ITNO .EQ. ITMAX ) THEN
        IQUIT = 1
C
C TEST FOR STAGNATION.
C
      ELSE 
        TEMP = KD*LOG(TOL/RSN)/
     &    LOG(RSN/((1.0 + 10.0*EPSILON(1.0_R8))*RSNOLD))
        IF (TEMP .GE. 40.0*(ITMAX - ITNO)) THEN
          IF (KDMAX .LE. KDLIMIT-M) THEN
            IFLAGI = 3
            KDMAX = KDMAX + M ! INCREASE KDMAX BY M
            GO TO 200
          ELSE IF (KDMAX .NE. KDLIMIT) THEN
            IFLAGI = 3 
            KDMAX = KDLIMIT
            GO TO 200 
          END IF
        ENDIF 
      ENDIF             
C
C COMPUTE THE CORRECTION IN V(:,1).
C
      DO I = KD, 1, -1
        W(I) = W(I)/R(I,I)
        IF (I .GT. 1) W(1:I-1) = W(1:I-1)-W(I)*R(1:I-1,I) 
      END DO
C
C COMPUTE KD ORTHOGONAL VECTORS FROM HOUSEHOLDER VECTORS.
C
      VTEMP = 0.0
      VTEMP(1:KD)=W(1:KD)
      DO I=KD,1,-1
        TEMP=V(I,I)
        V(I,I)=1.0
        CALL HREFX(N-I+1,V(I:N,I), TEMP, VTEMP(I:N))
        V(I,I)=TEMP
      END DO
C
C ADD THE CORRECTION TO THE APPROXIMATE SOLUTION. 
C
      IF (PRESENT(ROWPOSP) .AND. PRESENT(COLPOSP)) THEN    ! MODE = 2
        CALL ILUSOLVDS(N, PRECON(1:ROWPOSP(N+1)-1), 
     &     ROWPOSP(N+1)-1, ROWPOSP(1:N+1), 
     &     COLPOSP(1:ROWPOSP(N+1)-1), VTEMP(1:N))
      ELSE                                                 ! MODE = 1
        CALL SOLVDS(N, PRECON, ROWPOS(N+1)-1, ROWPOS(1:N+1),
     &     VTEMP(1:N))
      ENDIF
      X = X+VTEMP
C
C UPDATE THE RESIDUAL VECTOR BY DIRECT EVALUATION IN (V:,1) AND RECOMPUTE 
C THE RESIDUAL NORM.
C
      IF (IQUIT .EQ. 0 ) THEN
        VTEMP=X
        IF (PRESENT(ROWPOSP) .AND. PRESENT( COLPOSP)) THEN
          CALL MULM2(V(:,1),VTEMP)            ! MODE=2
        ELSE
          CALL MULM1(V(:,1),VTEMP)            ! MODE=1
        END IF
        V(:,1) = RHS - V(:,1)
        RSN = DNRM2(N, V(:,1), 1)
      ENDIF
C
C TERMINATE, OR GO TO THE TOP OF THE OUTER LOOP.
C
      IF (RSN .LE. TOL) THEN 
        IFLAGC = 0
        EXIT OUTER 
      ENDIF
      IF (ITNO .EQ. ITMAX) THEN 
        IF (IPRINT .GT. 0) THEN
          WRITE (IPRINT, 250) ITMAX 
 250      FORMAT(/,4X, 'MAXIMUM NUMBER OF ITERATIONS REACHED:',I7)
        ENDIF
        IFLAGC = 4 
        EXIT OUTER
      END IF 
      IF (IQUIT .EQ. 0 ) THEN
        IF (RSN .GT. RSNOLD ) THEN
          IF (IPRINT .GT. 0) THEN
            WRITE (IPRINT, 260) RSN     
 260        FORMAT(/,4X, '(TRUE) RESIDUAL NORM INCREASED TO',ES16.8)
          ENDIF
          IF (RSN .LE. TOL**(2.0/3.0))  THEN
            IFLAGC = 0
          ELSE  IF (STRONG_VERSION) THEN
            IF (IPRINT .GT. 0) THEN
              WRITE (IPRINT, 270) TOL**(2.0/3.0)                
 270          FORMAT(/,4X,'(TRUE) RESIDUAL NORM IS LARGER THAN',ES16.8)
            ENDIF
            IFLAGC = 4
          END IF
          EXIT OUTER
        END IF
C
C TEST FOR STAGNATION USING TRUE RESIDUAL NORM.
C
        TEMP = KD*LOG(TOL/RSN)/
     &     LOG(RSN/((1.0 + 10.0*EPSILON(1.0_R8))*RSNOLD))
        IF (TEMP .GE. 70.0*(ITMAX - ITNO)) THEN
          IFLAGI = 3
          IF ( KDMAX .LE. KDLIMIT-M ) THEN
            KDMAX = KDMAX + M            ! INCREASE KDMAX
          ELSE IF (KDMAX .NE. KDLIMIT) THEN
            KDMAX = KDLIMIT
          END IF
        ELSE  IF (TEMP .GE. 1000.0*(ITMAX - ITNO)) THEN
          IFLAGC = 4
          EXIT OUTER
        END IF
      END IF
      END DO OUTER
C
C RETURN. 
C
      ITMAX = ITNO
      IF (IFLAGC .NE. 0) IFLAG = IFLAGC
      CALL CLEANUP
      RETURN
C
C INTERNAL SUBROUTINES.
C
      CONTAINS
      SUBROUTINE MULM1(Y,X)
C
C MATRIX-VECTOR MULTIPLY FOR MODE = 1.
C 
        REAL (KIND=R8), INTENT(IN):: X(N) 
        REAL (KIND=R8), INTENT(OUT):: Y(N) 
        Y = 0.0
        CALL MULTDS(Y(1:N-1), AA, X(1:N-1), ROWPOS(1:N),
     &    N-1, ROWPOS(N)-1)       
C
C       RESULT MODIFIED ACCORDING TO KVAL.
C
        Y(KVAL) = Y(KVAL)+X(N)
        IF (KVAL .LT. N) 
     &    Y(N) = X(KVAL) + X(N)*(1.0/ABS(AA(ROWPOS(KVAL)))+1.0)
        RETURN
      END SUBROUTINE MULM1
C       
      SUBROUTINE MULM2(Y,X)
C
C MATRIX-VECTOR MULTIPLY FOR MODE = 2.
C
        REAL (KIND=R8), INTENT(IN):: X(N)
        REAL (KIND=R8), INTENT(OUT):: Y(N) 
        INTERFACE
          SUBROUTINE MULT2DS(Y, B, X, ROWPOS, COLPOS, N, LENB)
            USE REAL_PRECISION
            INTEGER, INTENT (IN):: LENB, N, ROWPOS(N+1),   
     &        COLPOS(LENB)
            REAL (KIND=R8), INTENT(IN):: X(:), B(LENB)
            REAL (KIND=R8), INTENT (OUT):: Y(N)
          END SUBROUTINE MULT2DS
        END INTERFACE
        Y = 0.0
        IF (IFLAG .NE. -2) THEN 
          CALL MULT2DS(Y(1:N-1), AA, X(1:N-1), ROWPOS(1:N), 
     &      COLPOS(1:LENAA), N-1, LENAA)
          Y(1:N-1) = Y(1:N-1)-X(N)*PP(1:N-1)
        ELSE
          CALL MULT2DS(Y(1:N-1), AA, X(1:N), ROWPOS(1:N),
     &      COLPOS(1:LENAA), N-1, LENAA)
        END IF
C
C       RESULT MODIFIED ACCORDING TO KVAL.
C
        Y(N) = X(KVAL)
        RETURN
      END SUBROUTINE MULM2
      SUBROUTINE CLEANUP    ! DEALLOCATE TEMPORARY LOCAL STORAGE.
        IF(ALLOCATED(C)) DEALLOCATE(C)
        IF(ALLOCATED(R)) DEALLOCATE(R)
        IF(ALLOCATED(S)) DEALLOCATE(S)
        IF(ALLOCATED(SVBIG)) DEALLOCATE(SVBIG)
        IF(ALLOCATED(SVSML)) DEALLOCATE(SVSML)
        IF(ALLOCATED(V)) DEALLOCATE(V)
        IF(ALLOCATED(W)) DEALLOCATE(W)
      END SUBROUTINE CLEANUP
C 
      SUBROUTINE HREFG(N, X, TAU, ALPHA)
C
C HOUSEHOLDER VECTOR CONSTRUCTION. HREFG IS  THE LAPACK AUXILIARY 
C ROUTINE DLARFG WITH TRIVIAL MODIFICATIONS FOR COMPATIBILITY WITH
C HOMPACK90.
C
      INTEGER, INTENT(IN):: N
      REAL (KIND=R8), INTENT(OUT):: ALPHA, TAU
      REAL (KIND=R8), INTENT(IN OUT):: X(N-1)
!
!  Purpose
!  =======
!
!  HREFG generates a real elementary reflector H of order N, such
!  that
!
!        H * ( ALPHA ) = ( BETA ),   H' * H = I.
!            (   X   )   (   0  )
!
!  where ALPHA and BETA are scalars, and X is an (N-1)-element real
!  vector. H is represented in the form
!
!        H = I - TAU * ( 1 ) * ( 1 V' ) ,
!                      ( V )
!
!  where TAU is a real scalar and V is a real (N-1)-element
!  vector.
!
!  If the elements of X are all zero, then TAU = 0 and H is taken to be
!  the unit matrix.
!
!  Otherwise  1 <= TAU <= 2.  
!
!  Arguments
!  =========
!
!  N       (input)
!          The order of the elementary reflector.
!
!
!  X       (input/output) 
!          On entry, the vector X.
!          On exit, it is overwritten with the vector V.
!
!  TAU     (output)
!          The value TAU.
! 
!  ALPHA   (INPUT/OUTPUT)
!          On entry, the value ALPHA.
!          On exit, it is overwritten with the value BETA.
!
!  =====================================================================
!
!     .. Parameters ..
      REAL (KIND=R8), PARAMETER:: ONE = 1.0_R8, ZERO = 0.0_R8
!     .. Local Scalars ..
      INTEGER::          J, KNT
      REAL (KIND=R8)::   BETA, RSAFMN, SAFMIN, XNORM
!       
!     .. External Functions ..
      INTERFACE
        FUNCTION DLAMCH(CMACH)
          USE REAL_PRECISION
          CHARACTER (LEN=1):: CMACH
          REAL (KIND=R8):: DLAMCH
        END FUNCTION DLAMCH
        FUNCTION DLAPY2(X,Y)
          USE REAL_PRECISION
          REAL (KIND=R8):: DLAPY2,X,Y
        END FUNCTION DLAPY2
        FUNCTION DNRM2(N,X,STRIDE)
          USE REAL_PRECISION
          INTEGER:: N,STRIDE
          REAL (KIND=R8):: DNRM2,X(N)
        END FUNCTION DNRM2
      END INTERFACE
!      
!     .. Executable Statements ..
!
      IF( N.LE.1 ) THEN
         TAU = ZERO
         RETURN
      END IF
      XNORM = DNRM2( N-1, X, 1)
      IF( XNORM.EQ.ZERO ) THEN
!
!        H  =  I
!
         TAU = ZERO
      ELSE
!
!        General case.
!
         BETA = -SIGN( DLAPY2( ALPHA, XNORM ), ALPHA )
         SAFMIN = DLAMCH( 'S' )
         IF( ABS( BETA ).LT.SAFMIN ) THEN
!
!           XNORM, BETA may be inaccurate; scale X and recompute them.
!
            RSAFMN = ONE / SAFMIN
            KNT = 0
            DO WHILE (ABS(BETA) < SAFMIN)
              KNT = KNT + 1
              X=RSAFMN*X
              BETA = BETA*RSAFMN
              ALPHA = ALPHA*RSAFMN
            END DO
!
!           New BETA is at most 1, at least SAFMIN.
!
            XNORM = DNRM2( N-1, X, 1 )
            BETA = -SIGN( DLAPY2( ALPHA, XNORM ), ALPHA )
            TAU = ( BETA-ALPHA ) / BETA
            X = X/( ALPHA-BETA )
!
!           If ALPHA is subnormal, it may lose relative accuracy.
!
            ALPHA = BETA
            DO J = 1, KNT
               ALPHA = ALPHA*SAFMIN
            END DO
         ELSE
            TAU = ( BETA-ALPHA) / BETA
            X = X/( ALPHA-BETA )
            ALPHA = BETA
         END IF
      END IF
      RETURN
      END SUBROUTINE HREFG
C
      SUBROUTINE HREFX(M, V, TAU, C)
C
C HOUSEHOLDER VECTOR APPLICATION. HREFX IS A FORTRAN 90 VERSION OF
C THE LAPACK ROUTINE DLARFX WITH MODIFICATIONS FOR
C COMPATIBILITY WITH HOMPACK 90.
C
      INTEGER, INTENT(IN):: M
      REAL (KIND=R8), INTENT(IN):: TAU, V(M)
      REAL (KIND=R8), INTENT(IN OUT):: C(M) 
!
!  Purpose
!  =======
!
!  HREFX applies a real elementary reflector H to a real M-element
!  vector C. H is represented in the form
!
!        H = I - TAU * V * V'
!
!  where TAU is a real scalar and V is a real vector.
!
!  If TAU = 0, then H is taken to be the unit matrix.
!
!  This version uses inline code if H has order < 11.
!
!  Arguments
!  =========
!
!  M       (input)
!          The order of the vector C.
!
!  V       (input) 
!          The vector V in the representation of H.
!
!  TAU     (input) 
!          The value TAU in the representation of H.
!
!  C       (input/output) 
!          On entry, the M-element vector C.
!          On exit, C is overwritten by H * C.
!
!  =====================================================================
!
!     ..  Parameters ..
      REAL (KIND=R8), PARAMETER:: ONE=1.0_R8, ZERO=0.0_R8
!
!     .. Local Scalars ..
      REAL (KIND=R8):: SUM, T1, T10, T2, T3, T4, T5, T6, T7, T8, T9,
     &                   V1, V10, V2, V3, V4, V5, V6, V7, V8, V9
!
!     .. Executable Statements ..  
!
      IF( TAU.EQ.ZERO ) RETURN
!
!        Form  H * C, where H has order M.
!
      SELECT CASE(M)
!
!     Code for general M.
!
      CASE (11:)
!
!     C := C - TAU * V * V' * C
!
           C = C - (TAU * DOT_PRODUCT(V,C)) * V
!
!        Special code for 1 x 1 Householder. 
!
      CASE(1)
           T1 = ONE - TAU*V( 1 )*V( 1 )
           C(1) = T1*C(1)
!
!        Special code for 2 x 2 Householder.
!
      CASE(2)
           V1 = V( 1 )
           T1 = TAU*V1
           V2 = V( 2 )
           T2 = TAU*V2
           SUM = V1*C(1) + V2*C(2)
           C(1) = C(1) - SUM*T1
           C(2) = C(2) - SUM*T2
!
!        Special code for 3 x 3 Householder.
!
      CASE(3)
           V1 = V( 1 )
           T1 = TAU*V1
           V2 = V( 2 )
           T2 = TAU*V2
           V3 = V( 3 )
           T3 = TAU*V3
           SUM = V1*C(1) + V2*C(2) + V3*C(3)
           C(1) = C(1) - SUM*T1
           C(2) = C(2) - SUM*T2
           C(3) = C(3) - SUM*T3
!
!        Special code for 4 x 4 Householder.
!
      CASE(4)
           V1 = V( 1 )
           T1 = TAU*V1
           V2 = V( 2 )
           T2 = TAU*V2
           V3 = V( 3 )
           T3 = TAU*V3
           V4 = V( 4 )
           T4 = TAU*V4
           SUM = V1*C(1) + V2*C(2) + V3*C(3) + V4*C(4)
           C( 1 ) = C( 1 ) - SUM*T1
           C( 2 ) = C( 2 ) - SUM*T2
           C( 3 ) = C( 3 ) - SUM*T3
           C( 4 ) = C( 4 ) - SUM*T4
!
!        Special code for 5 x 5 Householder.
!
      CASE(5)
           V1 = V( 1 )
           T1 = TAU*V1
           V2 = V( 2 )
           T2 = TAU*V2
           V3 = V( 3 )
           T3 = TAU*V3
           V4 = V( 4 )
           T4 = TAU*V4
           V5 = V( 5 )
           T5 = TAU*V5
           SUM = V1*C( 1 ) + V2*C( 2 ) + V3*C( 3 ) +
     &            V4*C( 4 ) + V5*C( 5 )
           C( 1 ) = C( 1 ) - SUM*T1
           C( 2 ) = C( 2 ) - SUM*T2
           C( 3 ) = C( 3 ) - SUM*T3
           C( 4 ) = C( 4 ) - SUM*T4
           C( 5 ) = C( 5 ) - SUM*T5
!
!        Special code for 6 x 6 Householder.
!
      CASE(6)
           V1 = V( 1 )
           T1 = TAU*V1
           V2 = V( 2 )
           T2 = TAU*V2
           V3 = V( 3 )
           T3 = TAU*V3
           V4 = V( 4 )
           T4 = TAU*V4
           V5 = V( 5 )
           T5 = TAU*V5
           V6 = V( 6 )
           T6 = TAU*V6
           SUM = V1*C( 1 ) + V2*C( 2 ) + V3*C( 3 ) +
     &            V4*C( 4 ) + V5*C( 5 ) + V6*C( 6 )
           C( 1 ) = C( 1 ) - SUM*T1
           C( 2 ) = C( 2 ) - SUM*T2
           C( 3 ) = C( 3 ) - SUM*T3
           C( 4 ) = C( 4 ) - SUM*T4
           C( 5 ) = C( 5 ) - SUM*T5
           C( 6 ) = C( 6 ) - SUM*T6
!
!        Special code for 7 x 7 Householder.
!
      CASE(7)
           V1 = V( 1 )
           T1 = TAU*V1
           V2 = V( 2 )
           T2 = TAU*V2
           V3 = V( 3 )
           T3 = TAU*V3
           V4 = V( 4 )
           T4 = TAU*V4
           V5 = V( 5 )
           T5 = TAU*V5
           V6 = V( 6 )
           T6 = TAU*V6
           V7 = V( 7 )
           T7 = TAU*V7
           SUM = V1*C( 1 ) + V2*C( 2 ) + V3*C( 3 ) +
     &            V4*C( 4 ) + V5*C( 5 ) + V6*C( 6 ) +
     &            V7*C( 7 )
           C( 1 ) = C( 1 ) - SUM*T1
           C( 2 ) = C( 2 ) - SUM*T2
           C( 3 ) = C( 3 ) - SUM*T3
           C( 4 ) = C( 4 ) - SUM*T4
           C( 5 ) = C( 5 ) - SUM*T5
           C( 6 ) = C( 6 ) - SUM*T6
           C( 7 ) = C( 7 ) - SUM*T7
!
!        Special code for 8 x 8 Householder.
!
      CASE(8) 
           V1 = V( 1 )
           T1 = TAU*V1
           V2 = V( 2 )
           T2 = TAU*V2
           V3 = V( 3 )
           T3 = TAU*V3
           V4 = V( 4 )
           T4 = TAU*V4
           V5 = V( 5 )
           T5 = TAU*V5
           V6 = V( 6 )
           T6 = TAU*V6
           V7 = V( 7 )
           T7 = TAU*V7
           V8 = V( 8 )
           T8 = TAU*V8
           SUM = V1*C( 1 ) + V2*C( 2 ) + V3*C( 3 ) +
     &            V4*C( 4 ) + V5*C( 5 ) + V6*C( 6 ) +
     &            V7*C( 7 ) + V8*C( 8 )
           C( 1 ) = C( 1 ) - SUM*T1
           C( 2 ) = C( 2 ) - SUM*T2
           C( 3 ) = C( 3 ) - SUM*T3
           C( 4 ) = C( 4 ) - SUM*T4
           C( 5 ) = C( 5 ) - SUM*T5
           C( 6 ) = C( 6 ) - SUM*T6
           C( 7 ) = C( 7 ) - SUM*T7
           C( 8 ) = C( 8 ) - SUM*T8
!
!        Special code for 9 x 9 Householder.
!
      CASE(9)
           V1 = V( 1 )
           T1 = TAU*V1
           V2 = V( 2 )
           T2 = TAU*V2
           V3 = V( 3 )
           T3 = TAU*V3
           V4 = V( 4 )
           T4 = TAU*V4
           V5 = V( 5 )
           T5 = TAU*V5
           V6 = V( 6 )
           T6 = TAU*V6
           V7 = V( 7 )
           T7 = TAU*V7
           V8 = V( 8 )
           T8 = TAU*V8
           V9 = V( 9 )
           T9 = TAU*V9
           SUM = V1*C( 1 ) + V2*C( 2 ) + V3*C( 3 ) +
     &            V4*C( 4 ) + V5*C( 5 ) + V6*C( 6 ) +
     &            V7*C( 7 ) + V8*C( 8 ) + V9*C( 9 )
           C( 1 ) = C( 1 ) - SUM*T1
           C( 2 ) = C( 2 ) - SUM*T2
           C( 3 ) = C( 3 ) - SUM*T3
           C( 4 ) = C( 4 ) - SUM*T4
           C( 5 ) = C( 5 ) - SUM*T5
           C( 6 ) = C( 6 ) - SUM*T6
           C( 7 ) = C( 7 ) - SUM*T7
           C( 8 ) = C( 8 ) - SUM*T8
           C( 9 ) = C( 9 ) - SUM*T9
!
!        Special code for 10 x 10 Householder.
!
      CASE (10)
           V1 = V( 1 )
           T1 = TAU*V1
           V2 = V( 2 )
           T2 = TAU*V2
           V3 = V( 3 )
           T3 = TAU*V3
           V4 = V( 4 )
           T4 = TAU*V4
           V5 = V( 5 )
           T5 = TAU*V5
           V6 = V( 6 )
           T6 = TAU*V6
           V7 = V( 7 )
           T7 = TAU*V7
           V8 = V( 8 )
           T8 = TAU*V8
           V9 = V( 9 )
           T9 = TAU*V9
           V10 = V( 10 )
           T10 = TAU*V10
           SUM = V1*C( 1 ) + V2*C( 2 ) + V3*C( 3 ) +
     &            V4*C( 4 ) + V5*C( 5 ) + V6*C( 6 ) +
     &            V7*C( 7 ) + V8*C( 8 ) + V9*C( 9 ) +
     &            V10*C( 10 )
           C( 1 ) = C( 1 ) - SUM*T1
           C( 2 ) = C( 2 ) - SUM*T2
           C( 3 ) = C( 3 ) - SUM*T3
           C( 4 ) = C( 4 ) - SUM*T4
           C( 5 ) = C( 5 ) - SUM*T5
           C( 6 ) = C( 6 ) - SUM*T6
           C( 7 ) = C( 7 ) - SUM*T7
           C( 8 ) = C( 8 ) - SUM*T8
           C( 9 ) = C( 9 ) - SUM*T9
           C( 10 ) = C( 10 ) - SUM*T10
      END SELECT   
      RETURN
      END SUBROUTINE HREFX
      END SUBROUTINE GMRES
      SUBROUTINE GMRILUDS(NN,LENAA,IFLAG,START,RHS)
C
C     This subroutine solves a system of equations using a
C     preconditioned adaptive GMRES(k).
C     GMRILUDS is the MODE=2 equivalent of subroutine PCGDS, which
C     is only for MODE=1 storage.
C
C     The linear system to be solved is  Mx=b.
C     If IFLAG = -1 or 0,
C 
C        +--          --+ 
C        |        |     |
C        |   AA   | -PP |
C    M = |        |     | ,
C        +--------+-----+ 
C        |    E(k)**t   | 
C        +--          --+
C
C        where AA is an (NN x NN) matrix stored in compressed-row format,
C        PP is an (NN x 1) vector.
C        It is assumed that rank [AA,-PP]=NN and M is invertible.
C      If IFLAG = -2,
C
C        +--          --+ 
C        |              |
C        |     AA       | 
C    M = |              | , 
C        +--------+-----+ 
C        |    E(k)**t   | 
C        +--          --+ 
C
C        where AA is an (NN x (NN+1)) matrix stored in compressed-row 
C        format. It is assumed that rank [AA]=NN and M is invertible.
C
C        +-   -+          +-    -+
C        |  0  |          |      |
C        | ... |          | -RHS |
C    b = |  0  |,  or b = |      |, 
C        +-----+          +------+
C        |  T  |          |  T   |
C        +-   -+          +-    -+
C
C        T = START(k), where |START(k)|= max|START(i)|
C                                       1<=i<=NN+1
C
C        b is of length NN+1 and E(k)**t is the ( 1 x (NN+1) ) vector
C        consisting of all zeros, except for a '1' in its k-th position.
C
C  Input variables:
C
C      NN -- row dimension of the matrix packed in AA.
C
C      LENAA -- number of elements in the packed array AA.
C
C      IFLAG -- indicator of how M is assembled.  
C
C      START -- vector of length NN+1, normally the solution to the
C               previous linear system; used to determine the index k .
C
C      RHS -- optional vector of length NN, used to define right hand
C             side for normal flow correction calculation.  It is
C             assumed that GMRILUDS is called without RHS present before
C             it is called with RHS present.  An ILU factorization based 
C             preconditioner is computed only when RHS is not present.
C
C  Input variables defined in module HOMOTOPY:
C
C      AA -- one dimensional real array containing a submatrix of M
C            in compressed-row storage form.  The logical dimensions
C            of AA depend on IFLAG.
C
C      ROWPOS -- integer array used for specifying information about AA.
C                Using compressed-row storage, it has length NN+2, and
C                stores the indices of where each row begins within AA.
C                ROWPOS(NN+1) = LENAA + 1 and ROWPOS(NN+2) = LENAA + 2.
C                (NOTE:  The value of ROWPOS(NN+2) is set by this
C                subroutine when the preconditioning matrix Q is
C                initialized.)
C
C      COLPOS -- integer array used for specifying information about AA.
C                Using compressed-row storage, it has length LENAA,
C                and contains the column indices of the corresponding
C                elements in AA.
C
C             For example, using the compressed-row storage scheme 
C             with IFLAG = -2, a 5 x 6 matrix of the form
C
C             +--               --+
C             |  1  3  0  0  0 10 |
C             |  3  2  0  7  0 11 |
C             |  0  0  4  6  0 12 |
C             |  0  7  6  5  9 13 |
C             |  0  0  0  9  8 14 |
C             +--               --+
C 
C             would result in NN=5, LENAA=18, ROWPOS=(1,4,8,11,16,19,*),
C             AA=(1,3,10,3,2,7,11,4,6,12,7,6,5,9,13,9,8,14),
C             COLPOS=(1,2,6,1,2,6,3,4,6,2,3,4,5,6,4,5,6).
C
C      PP -- vector of length NN, used for (NN+1)st column of
C            augmented matrix M when IFLAG = -1 or 0 .
C
C  Output variables:
C
C      START -- solution vector x of  M x = b  (defined above).
C
C      IFLAG -- normally unchanged on output.  If the GMRES
C               iteration fails to converge in 10*(NN+1) iterations (most
C               likely due to a singular Jacobian matrix), GMRILUDS returns
C               with  IFLAG = 4 , and does not compute x.
C
C  Calls subroutines ILUFDS and GMRES.
C
      USE HOMOTOPY, AA => QRSPARSE, WORK => PAR, IWORK => IPAR
      USE REAL_PRECISION
C
      INTEGER, INTENT(IN):: NN,LENAA
      INTEGER, INTENT(IN OUT):: IFLAG
      REAL (KIND=R8), INTENT(IN OUT):: START(NN+1)
      REAL (KIND=R8), INTENT(IN), OPTIONAL :: RHS(NN)
C
C LOCAL VARIABLES.
C
      INTEGER:: I, ITMAX, K, KD, NP1, QIND, ZBIND  
      INTEGER:: CIND, RIND, ROWL, STRT 
      REAL (KIND=R8):: STARTK
      REAL (KIND=R8):: RHSC(NN+1)
C
C GMRES PARAMETER.
C
      INTEGER, PARAMETER:: SUBSPACE=8         ! KRYLOV SUBSPACE VALUE.
C
      INTERFACE
        SUBROUTINE GMRES(N, KDMAX, ITMAX, RHSC, X, KVAL,
     &                Q, IFLAG, ROWPOSP, COLPOSP)
          USE HOMOTOPY
          USE REAL_PRECISION
          INTEGER, INTENT(IN):: KVAL, N
          INTEGER, INTENT(IN OUT):: IFLAG, ITMAX, KDMAX
          REAL (KIND=R8), INTENT(IN):: Q(:), RHSC(N)
          REAL (KIND=R8), INTENT(IN OUT):: X(N)
          INTEGER, INTENT(IN), OPTIONAL:: COLPOSP(:), ROWPOSP(:)
        END SUBROUTINE GMRES
        SUBROUTINE ILUFDS(NN, Q, LENQ, ROWPOS, COLPOS)
          USE REAL_PRECISION
          INTEGER, INTENT(IN):: LENQ, NN, COLPOS(LENQ), ROWPOS(NN+1) 
          REAL (KIND=R8), INTENT(IN OUT):: Q(LENQ) 
        END SUBROUTINE ILUFDS
      END INTERFACE 
C 
      NP1=NN+1 
C
C     INITIALIZE START POSITIONS WITHIN WORK AND IWORK.
C
      ZBIND = 1
      QIND = NP1+1
      RIND = 1
      CIND = NP1+2
C
      IF (.NOT. ALLOCATED(WORK)) THEN
        IF (IFLAG .EQ. -2) THEN
          ALLOCATE(WORK(NP1+LENAA+2)) 
        ELSE
          ALLOCATE(WORK(NP1+LENAA+NN+2))
        END IF
        WORK(1:NP1) = 0.0
      END IF
      IF (.NOT. ALLOCATED(IWORK)) THEN
        IF (IFLAG .EQ. -2) THEN
          ALLOCATE(IWORK(NP1+1+LENAA+2))
        ELSE
          ALLOCATE(IWORK(NP1+1+LENAA+NN+2))
        END IF
      END IF
C
C     FIND THE ELEMENT OF LARGEST MAGNITUDE IN THE INITIAL VECTOR, AND
C     RECORD ITS POSITION IN K.
C
      K = MAXVAL(MAXLOC(ABS(START)))
      STARTK = START(K)
C
C     SET VALUES OF ROWPOS(NN+1) AND ROWPOS(NN+2), AND
C     COMPUTE THE PRECONDITIONER Q FOR M.
C
      IF (.NOT. PRESENT(RHS)) THEN
        ROWPOS(NP1) = LENAA+1
        ROWPOS(NN+2) = LENAA+2
        IF (IFLAG .EQ. -2) THEN
          WORK(QIND:QIND+LENAA-1) = AA(1:LENAA)
          IWORK(RIND:RIND+NP1) = ROWPOS
          IWORK(CIND:CIND+LENAA-1) = COLPOS
        ELSE
C       MERGE AA AND -PP INTO Q ONLY FOR IFLAG >= -1. 
          IWORK(RIND) = 1
          DO I=1,NN
            STRT = IWORK(RIND+I-1)
            ROWL = ROWPOS(I+1)-ROWPOS(I)
            WORK(QIND+STRT-1:QIND+STRT+ROWL-2) = 
     &        AA(ROWPOS(I):ROWPOS(I+1)-1)
            WORK(QIND+STRT+ROWL-1) = -PP(I)
            IWORK(CIND+STRT-1:CIND+STRT+ROWL-2) =
     &        COLPOS(ROWPOS(I):ROWPOS(I+1)-1)
            IWORK(CIND+STRT+ROWL-1) = NN+1
            IWORK(RIND+I) = ROWPOS(I+1)+I
          END DO 
          IWORK(RIND+NP1) = ROWPOS(NN+2)+NN
        END IF
        WORK(QIND+IWORK(RIND+NN)-1) = 1.0
        IWORK(CIND+IWORK(RIND+NN)-1) = K
        IF (K. LT. NP1) THEN
          WORK(QIND+IWORK(RIND+NN)) = 0.0
          IWORK(CIND+IWORK(RIND+NN)) = NP1
          IWORK(RIND+NP1) = IWORK(RIND+NP1)+1
        END IF
        CALL ILUFDS(NP1, WORK(QIND:QIND+IWORK(RIND+NP1)-2),
     &    IWORK(RIND+NP1)-1, IWORK(RIND:RIND+NP1),
     &    IWORK(CIND:CIND+IWORK(RIND+NP1)-2))
      END IF
C
C     COMPUTE RIGHT HAND SIDE FOR Mx=b.
C
      RHSC(NP1) = STARTK
      IF (PRESENT(RHS)) THEN 
        RHSC(1:NN) = -RHS
      ELSE
        RHSC(1:NN) = 0.0
      END IF
C
C CALL GMRES FOR THE Mx=b SYSTEM.
C
      ITMAX=15*NP1
      KD=SUBSPACE          
      CALL GMRES(NP1, KD, ITMAX, RHSC, WORK(ZBIND:ZBIND+NN), 
     &  K, WORK(QIND:QIND+IWORK(RIND+NP1)-2), IFLAG,
     &  ROWPOSP=IWORK(RIND:RIND+NP1),
     &  COLPOSP=IWORK(CIND:CIND+IWORK(RIND+NP1)-2))
      IF ( IFLAG .GT. 0) RETURN
      START(1:NP1) = WORK(ZBIND:ZBIND+NN)
      RETURN
      END SUBROUTINE GMRILUDS
      SUBROUTINE HFUN1P(QDG,LAMBDA,X,
     & PDG,CL,COEF,RHO,
     & DRHOX,DRHOL,XDGM1,XDG,
     & G,DG,PXDGM1,PXDG,
     & F,DF,XX,TRM,
     & DTRM,CLX,DXNP1,
     & N,MAXT,IDEG,
     & NUMT,KDEG)
C
C  HFUN1P  EVALUATES THE CONTINUATION EQUATION "RHO".
C
C  NOTE THAT:
C    DRHOX IS THE "REALIFICATION" OF DCRHOX, WHERE
C    DCRHOX DENOTES THE (COMPLEX) PARTIAL
C    DERIVATIVE MATRIX OF THE CONTINUATION SYSTEM
C    WITH RESPECT TO X,  AND
C    DRHOL IS THE "REALIFICATION" OF DCRHOL, WHERE
C    DCRHOL DENOTES THE (COMPLEX) PARTIAL
C    DERIVATIVE MATRIX OF THE CONTINUATION SYSTEM
C    WITH RESPECT TO LAMBDA. THUS
C      DRHOX(2J-1,2K-1) = DCRHOX(1,J,K)
C      DRHOX(2J  ,2K  ) = DCRHOX(1,J,K)
C      DRHOX(2J-1,2K  ) =-DCRHOX(2,J,K)
C      DRHOX(2J  ,2K-1) = DCRHOX(2,J,K)
C      DRHOL(2J-1,N2P1) = DCRHOL(1,J)
C      DRHOL(2J  ,N2P1) = DCRHOL(2,J)
C       RHO(2J-1)      = CRHO(1,J)
C       RHO(2J  )      = CRHO(2,J)
C    WHERE CRHO DENOTES THE (COMPLEX) CONTINUATION SYSTEM,
C    THE INITIAL "1" OR "2" DENOTES REAL OR IMAGINARY PARTS,
C    RESPECTIVELY, "J" INDEXES THE EQUATION, "K" INDEXES THE PARTIAL
C    DERIVATIVE, AND NEITHER DCRHOX NOR DCRHOL ARE PROGRAM VARIABLES.
C
C  ON INPUT:
C
C    QDG  IS THE "RANDOM" PARAMETER "A".
C
C    LAMBDA  IS THE CONTINUATION PARAMETER.
C
C    X    IS THE INDEPENDENT VARIABLE.
C
C    PDG  IS ONE OF THE PARAMETERS THAT DEFINES G (SEE SUBROUTINE
C         GFUNP).
C
C    CL   IS ONE OF THE PARAMETERS THAT DEFINES F (SEE SUBROUTINE
C         FFUNP).
C
C    COEF  IS ONE OF THE PARAMETERS THAT DEFINES F (SEE SUBROUTINE
C         FFUNP).
C
C  ON OUTPUT:
C
C    RHO    IS THE HOMOTOPY.
C
C    DRHOX  CONTAINS THE PARTIAL DERIVATIVES OF RHO WITH RESPECT
C         TO X.
C
C    DRHOL  CONTAINS THE PARTIAL DERIVATIVES OF RHO WITH RESPECT
C         TO LAMBDA.
C
C  THE FOLLOWING ARE VARIABLES WHOSE WORKSPACE IS PASSED FROM HFUNP:
C    XDGM1
C    XDG
C    G
C    DG
C    PXDGM1
C    PXDG
C    F
C    DF
C    XX
C    TRM
C    DTRM
C    CLX
C    DXNP1
C    N
C    MAXT
C    IDEG
C    NUMT
C    KDEG
C
C  OTHER VARIABLES:
C    ONEML
C
C  SUBROUTINES:  GFUNP, FFUNP.
C
      USE REAL_PRECISION
C DECLARATION OF INPUT, WORKSPACE, AND OUTPUT:
      INTEGER, INTENT(IN):: N,MAXT,IDEG(N),NUMT(N),KDEG(N,N+1,MAXT)
      REAL (KIND=R8), INTENT(IN):: QDG(2,N),LAMBDA,X(2,N),
     &  PDG(2,N),CL(2,N+1),COEF(N,MAXT)
      REAL (KIND=R8), INTENT(OUT):: RHO(2*N),DRHOX(2*N,2*N),DRHOL(2*N)
      REAL (KIND=R8), INTENT(IN OUT):: XDGM1(2,N),XDG(2,N),
     &  G(2,N),DG(2,N),PXDGM1(2,N),PXDG(2,N),
     &  F(2,N), DF(2,N,N+1),XX(2,N,N+1,MAXT),TRM(2,N,MAXT),
     &  DTRM(2,N,N+1,MAXT),CLX(2,N),DXNP1(2,N)
C
C DECLARATION OF LOCAL VARIABLES:
      INTEGER:: J,J2,J2M1,K,K2,K2M1
      REAL (KIND=R8):: ONEML
C
      CALL GFUNP(N,IDEG,PDG,QDG,X,XDGM1,XDG,PXDGM1,PXDG,G,DG)
      CALL FFUNP(N,NUMT,MAXT,KDEG,COEF,CL,X,XX,TRM,DTRM,CLX,DXNP1,F,DF)
      ONEML=1.0 - LAMBDA
      DO J=1,N
          J2=2*J
          J2M1=J2-1
          DO K=1,N
              K2=2*K
              K2M1=K2-1
              DRHOX(J2M1,K2M1)= LAMBDA*DF(1,J,K)
              DRHOX(J2  ,K2  )= DRHOX(J2M1,K2M1)
              DRHOX(J2  ,K2M1)= LAMBDA*DF(2,J,K)
              DRHOX(J2M1,K2  )=-DRHOX(J2  ,K2M1)
          END DO
          DRHOX(J2M1,J2M1)= DRHOX(J2M1,J2M1) + ONEML*DG(1,J)
          DRHOX(J2  ,J2  )= DRHOX(J2M1,J2M1)
          DRHOX(J2  ,J2M1)= DRHOX(J2  ,J2M1) + ONEML*DG(2,J)
          DRHOX(J2M1,J2  )=-DRHOX(J2  ,J2M1)
          DRHOL(J2M1)     =   F(1,J)      -        G(1,J)
          DRHOL(J2)       =   F(2,J)      -        G(2,J)
          RHO(J2M1)      = LAMBDA*F(1,J) + ONEML* G(1,J)
          RHO(J2  )      = LAMBDA*F(2,J) + ONEML* G(2,J)
      END DO
      RETURN
      END SUBROUTINE HFUN1P
      SUBROUTINE HFUNP(N,QDG,LAMBDA,X)
C
C HFUNP ALLOCATES STORAGE FOR SUBROUTINE HFUN1P FROM THE WORK ARRAYS
C PAR AND IPAR, AS FOLLOWS:
C
C DOUBLE PRECISION VARIABLES AND ARRAYS PASSED IN PAR
C
C     PAR INDEX     VARIABLE NAME       LENGTH
C    ----------     -------------    -----------------
C          1              PDG               2*N
C          2               CL               2*(N+1)
C          3             COEF               N*MAXT
C          4              RHO               N2
C          5              DRHOX             N2*N2
C          6              DRHOL             N2
C          7            XDGM1               2*N
C          8              XDG               2*N
C          9              G                 2*N
C         10             DG                 2*N
C         11           PXDGM1               2*N
C         12             PXDG               2*N
C         13               F                2*N
C         14              DF                2*N*(N+1)
C         15               XX               2*N*(N+1)*MAXT
C         16              TRM               2*N*MAXT
C         17             DTRM               2*N*(N+1)*MAXT
C         18              CLX               2*N
C         19            DXNP1               2*N
C
C INTEGER VARIABLES AND ARRAYS PASSED IN IPAR
C
C    IPAR INDEX     VARIABLE NAME       LENGTH            OFFSET
C    ----------     -------------    -----------------
C          1                N               1               1
C          2             MAXT               1               2
C          3            PROFF               25              3
C          4           IPROFF               15              28
C          5             IDEG               N               43
C          6             NUMT               N               43+N
C          7             KDEG               N*(N+1)*MAXT   43+N2
C
C ON INPUT:
C
C N  IS THE NUMBER OF EQUATIONS AND VARIABLES.
C
C QDG  IS THE "RANDOM" VECTOR DENOTED  "A"  IN HOMPACK DOCUMENTATION.
C
C LAMBDA  IS THE CONTINUATION PARAMETER.
C
C X  IS THE INDEPENDENT VARIABLE.
C
C ON OUTPUT:
C
C THE GLOBAL WORK ARRAYS PAR AND IPAR HAVE BEEN UPDATED.
C
C SUBROUTINES:  HFUN1P.
C
      USE HOMOTOPY
      USE REAL_PRECISION
      INTEGER, INTENT(IN):: N
      REAL (KIND=R8), INTENT(IN):: QDG(2,N),LAMBDA,X(2,N)
C
      CALL HFUN1P(QDG,LAMBDA,X,
     & PAR( IPAR(3 + ( 1-1))), PAR( IPAR(3 + ( 2-1))),
     & PAR( IPAR(3 + ( 3-1))), PAR( IPAR(3 + ( 4-1))),
     & PAR( IPAR(3 + ( 5-1))), PAR( IPAR(3 + ( 6-1))),
     & PAR( IPAR(3 + ( 7-1))), PAR( IPAR(3 + ( 8-1))),
     & PAR( IPAR(3 + ( 9-1))), PAR( IPAR(3 + (10-1))),
     & PAR( IPAR(3 + (11-1))), PAR( IPAR(3 + (12-1))),
     & PAR( IPAR(3 + (13-1))), PAR( IPAR(3 + (14-1))),
     & PAR( IPAR(3 + (15-1))), PAR( IPAR(3 + (16-1))),
     & PAR( IPAR(3 + (17-1))), PAR( IPAR(3 + (18-1))),
     & PAR( IPAR(3 + (19-1))),
     &IPAR( IPAR(28+ ( 1-1))),IPAR( IPAR(28+ ( 2-1))),
     &IPAR( IPAR(28+ ( 5-1))),IPAR( IPAR(28+ ( 6-1))),
     &IPAR( IPAR(28+ ( 7-1))) )
C
      RETURN
      END SUBROUTINE HFUNP
      SUBROUTINE ILUFDS(NP1, B, LENB, ROWPOSP, COLPOSP)
c
C     Computes the incomplete LU factorization of the matrix B,
C     where B is NP1 x NP1.  B is assumed to be stored in the general
C     sparse row scheme.
C
C     The method used is that found in TR 89-41, Department of Computer
C     Science, VPI&SU, Blacksburg, VA, 1989: 'Preconditioned
C     conjugate gradient algorithms for homotopy curve tracking',
C     page 10.
C---------------------------------------------------------------------
C
C     Input variables:
C       B       matrix to be factored.
C       ROWPOSP  indices of row-starts within B.
C       COLPOSP  column indices for matrix B stored in the general
C                sparse row scheme.
C       LENB    number of entries in B.
C       NP1     the dimension of B.
C
C     Output variables:
C       B       the ILU factors of input matrix B.
C-------------------------------------------------------------
C
      USE REAL_PRECISION
      INTEGER, INTENT(IN):: LENB, NP1, ROWPOSP(NP1+1), COLPOSP(LENB)
      REAL (KIND=R8), INTENT(IN OUT):: B(LENB)
C
C LOCAL VARIABLES
      REAL (KIND=R8):: SIJ, LIT, LII
      INTEGER:: I, J, COUNT, ISTRT, IFIN, TMAX, K, T, M
C
C
      DO I = 1, NP1
         ISTRT = ROWPOSP(I)
         IFIN  = ROWPOSP(I+1) - 1
C---------------------------------------  FOR EACH ELEMENT IN ROW I,
C                                         COMPUTE THE COLUMN NUMBER
C                                         AND T.
         DO COUNT = ISTRT, IFIN
            J = COLPOSP(COUNT)
            TMAX = MIN0(I,J) - 1
            SIJ = B(COUNT)
C---------------------------------------  COMPUTE THE CORRESPONDING
C                                         SUM OF PRODUCTS OF ELEMENTS
C                                         OF L AND U.
            K = ISTRT
 42         T = COLPOSP(K)
            IF (T .LE. TMAX) THEN
              LIT = B(K)
              M = ROWPOSP(T)
C---------------------------------------  FIND VALUE OF U_{TJ}.
 20           IF (COLPOSP(M) .LT. J) THEN
                 M = M + 1
                 GOTO 20
              ENDIF
              IF (COLPOSP(M) .EQ. J) SIJ = SIJ - LIT*B(M)
              K = K + 1
              GOTO 42
           ENDIF
C---------------------------------------  END OF 'T' LOOP.
           K = ISTRT
C---------------------------------------  FIND VALUE OF L_{II}.
 30        IF (COLPOSP(K) .LT. I) THEN
              K = K+1
              GOTO 30
           ENDIF
           IF (COLPOSP(K) .EQ. I) THEN
              LII = B(K)
              IF (DABS(LII) .EQ. 0.0) THEN
                 LII = 0.00001
                 B(K) = 0.00001
              ENDIF
           ELSE
              LII = 0.00001
           ENDIF
C---------------------------------------  UPDATE L OR U, AS NEEDED.
           IF (I .GE. J) THEN
              B(COUNT) = SIJ
           ELSE
              B(COUNT) = SIJ/LII
           ENDIF
        END DO
C---------------------------------------  END OF 'COUNT' LOOP.
      END DO
C
      RETURN
      END SUBROUTINE ILUFDS
      SUBROUTINE ILUSOLVDS(NN, Q, LENQ, ROWPOSP, COLPOSP, B)
C
C     Computes Q^{-1}*B -- returns result as B.
C---------------------------------------------------------------------
C
C     Input variables:
C
C       Q        triangular factors of preconditioning matrix, stored
C                in the general sparse row scheme.
C       ROWPOSP  indices of row-starts within B.
C       COLPOSP  column indices for matrix B stored in the general
C                sparse row scheme.
C       NN       logical row dimension of Q.
C       LENQ     number of data entries in Q.
C       B        right hand side -- should have dimension NN.
C
C     Output variables:
C
C       B        solution of Qx = B.
C
C---------------------------------------------------------------------
C
      USE REAL_PRECISION
      INTEGER, INTENT(IN) ::  NN, LENQ, ROWPOSP(NN+1), COLPOSP(LENQ)
      REAL (KIND=R8), INTENT(IN) :: Q(LENQ) 
      REAL (KIND=R8), INTENT(IN OUT) :: B(NN)
C 
C LOCAL VARIABLES
      INTEGER:: DIAG(NN), I, K, J
C------------------------------------------------ COMPUTE B = INV(L)*B
      B(1) = B(1)/Q(1)
      DIAG(1) = 1
      DO  I = 2, NN
        K = ROWPOSP(I)
 42     J = COLPOSP(K)
        IF (J .LT. I) THEN
          B(I) = B(I) - Q(K)*B(J)
          K = K + 1
          GOTO 42
        ELSE 
          DIAG(I) = K
          B(I) = B(I)/Q(K)
        ENDIF
      END DO
C------------------------------------------------ COMPUTE B = INV(U)*B
      DO  I = NN-1, 1, -1
        DO  K = DIAG(I)+1, ROWPOSP(I+1)-1
          B(I) = B(I) - Q(K)*B(COLPOSP(K))
        END DO
      END DO
      RETURN
      END SUBROUTINE ILUSOLVDS
      SUBROUTINE INITP(IFLG1,N,NUMT,KDEG,COEF,
     &                              IDEG,FACV,CL,PDG,QDG,R)
C
C INITP  INITIALIZES THE CONSTANTS THAT DEFINE THE POLSYS HOMOTOPY,
C INITIALIZES THE CONSTANTS THAT DEFINE THE PROJECTIVE TRANSFORMATION,
C AND SCALES THE COEFFICIENTS (IF SCALING IS SPECIFIED).
C
C ON INPUT:
C
C IFLG1  IS A FLAG THAT SPECIFIES WHETHER THE COEFFICIENTS ARE TO
C   BE SCALED OR NOT AND WHETHER THE PROJECTIVE TRANSFORMATION IS TO
C   BE USED OR NOT.  IFLG1=A*10+B.  SCALING IS SPECIFIED WHEN B=1.  THE 
C   PROJECTIVE TRANSFORMATION IS SPECIFIED WHEN A=1.  OTHERWISE, A AND/OR 
C   B =0.  SCALING IS EVOKED BY A CALL TO THE SUBROUTINE  SCLGNP.  THE 
C   PROJECTIVE TRANSFORMATION IS EVOKED BY SETTING THE  CL  ARRAY EQUAL
C   TO RANDOM COMPLEX NUMBERS.  OTHERWISE,  CL  IS SET TO NOMINAL VALUES.
C
C N  IS THE NUMBER OF EQUATIONS AND VARIABLES.
C
C NUMT(J)  IS THE NUMBER OF TERMS IN EQUATION J, FOR J=1 TO N.
C
C KDEG(J,L,K)  IS THE DEGREE OF THE L-TH VARIABLE, X(L), IN THE K-TH
C  TERM OF THE J-TH EQUATION, WHERE J=1 TO N, L=1 TO N+1, AND K=1 TO 
C  NUMT(J).  THE CASE "L=N+1" IS SPECIAL, AND  KDEG  IS NOT AN INPUT
C  VALUE TO  POLSYS , BUT RATHER IS COMPUTED IN THIS SUBROUTINE.  
C 
C COEF(J,K)  IS THE COEFFICIENT OF THE K-TH TERM FOR THE J-TH
C   EQUATION, WHERE J=1 TO N AND K=1 TO NUMT(J).
C
C
C ON OUTPUT:
C
C IDEG(J)  IS THE DEGREE OF THE J-TH EQUATION FOR J=1 TO N.
C
C FACV(J)  IS THE SCALE FACTOR FOR THE J-TH VARIABLE.
C
C CL(2,1:N+1)  IS AN ARRAY USED TO DEFINE THE PROJECTIVE
C   TRANSFORMATION.  IT IS USED IN SUBROUTINES  FFUNP  AND  OTPUTP
C   TO DEFINE THE PROJECTIVE COORDINATE, XNP1.    
C
C PDG  IS USED IN SUBROUTINE  GFUNP  TO DEFINE THE INITIAL SYSTEM,
C   G(X)=0.
C
C QDG  IS USED IN SUBROUTINE  GFUNP  TO DEFINE THE INITIAL SYSTEM,
C   G(X)=0.
C
C R  IS USED IN SUBROUTINE  STRPTP  TO GENERATE SOLUTIONS TO G(X)=0.
C
C
      USE HOMOTOPY
      USE REAL_PRECISION
C DECLARATIONS OF INPUT AND OUTPUT:
      INTEGER, INTENT(IN):: IFLG1,N,NUMT(:)
      INTEGER, INTENT(IN OUT):: KDEG(:,:,:)
      REAL (KIND=R8), INTENT(IN OUT):: COEF(:,:)
      INTEGER, INTENT(OUT):: IDEG(N)
      REAL (KIND=R8), INTENT(OUT):: 
     &  FACV(N),CL(2,N+1),PDG(2,N),QDG(2,N),R(2,N)
C
C DECLARATIONS OF LOCAL VARIABLES:
      INTEGER:: IERR,J,JJ,MAXT,N2,NP1
      REAL (KIND=R8):: CCL(2,11),P(2,10),Q(2,10),ZERO
C
      INTERFACE
        SUBROUTINE DIVP(XXXX,YYYY,ZZZZ,IERR)
        USE REAL_PRECISION
        REAL (KIND=R8), DIMENSION(2), INTENT(IN):: XXXX,YYYY
        INTEGER, INTENT(OUT):: IERR
        REAL (KIND=R8), DIMENSION(2), INTENT(OUT):: ZZZZ
        END SUBROUTINE DIVP
        SUBROUTINE MULP(XXXX,YYYY,ZZZZ)
        USE REAL_PRECISION
        REAL (KIND=R8), DIMENSION(2), INTENT(IN):: XXXX,YYYY
        REAL (KIND=R8), DIMENSION(2), INTENT(OUT):: ZZZZ
        END SUBROUTINE MULP
        SUBROUTINE POWP(NNNN,XXXX,YYYY)
        USE REAL_PRECISION
        INTEGER, INTENT(IN):: NNNN
        REAL (KIND=R8), DIMENSION(2), INTENT(IN):: XXXX
        REAL (KIND=R8), DIMENSION(2), INTENT(IN OUT):: YYYY
        END SUBROUTINE POWP
        SUBROUTINE SCLGNP(N,MAXT,NUMT,DEG,MODE,EPS0,COEF,
     &    NNUMT,DDEG,CCOEF,ALPHA,BETA,RWORK,XWORK,
     &    FACV,FACE,COESCL,IERR)
        USE REAL_PRECISION
        INTEGER, INTENT(IN):: N,MAXT,NUMT(:),DEG(:,:,:),MODE
        REAL (KIND=R8), INTENT(IN):: EPS0,COEF(:,:)
        INTEGER, INTENT(IN OUT):: NNUMT(N),DDEG(N,N+1,MAXT)
        REAL (KIND=R8), INTENT(IN OUT):: CCOEF(N,MAXT),ALPHA(2*N,2*N),
     &    BETA(2*N),RWORK(N*(2*N+1)),XWORK(2*N)
        REAL (KIND=R8), INTENT(OUT):: FACV(N),FACE(N),COESCL(N,MAXT)
        INTEGER, INTENT(OUT):: IERR
        END SUBROUTINE SCLGNP
      END INTERFACE
C
      MAXT = MAXVAL(NUMT)
      N2 =2*N
      NP1=N+1
      ZERO=0.0
      DO J=1,N
        IDEG(J)=MAXVAL(SUM(KDEG(J,1:N,1:NUMT(J)),DIM=1))
      END DO
      DO J=1,N
        KDEG(J,NP1,1:NUMT(J))=IDEG(J)-SUM(KDEG(J,1:N,1:NUMT(J)),DIM=1)
      END DO
      IF ( IFLG1 .EQ. 10  .OR.  IFLG1 .EQ. 00) THEN
C
C       DON'T SCALE THE COEFFICIENTS.  SET  FACV  EQUAL TO NOMINAL 
C       VALUES.
C
        FACV = 0.0
      ELSE
C
C SET UP THE WORKSPACE FOR SUBROUTINE  SCLGNP  AND CALL  SCLGNP  TO
C SCALE THE COEFFICIENTS.
C
C*****************************************************************
C VARIABLES THAT ARE PASSED IN ARRAY PAR.
C
C    VARIABLE NAME   LENGTH        OFFSET
C
C    1   CCOEF       N*MAXT        1
C    2   ALPHA       4*N**2        1+N*MAXT
C    3   BETA        2*N           1+N*MAXT+4*N**2
C    4   RWORK       N*(2*N+1)     1+N*MAXT+4*N**2+2*N
C    5   XWORK       2*N           1+N*MAXT+4*N**2+2*N+N*(2*N+1)
C    6   FACE        N             1+N*MAXT+4*N**2+4*N+N*(2*N+1)
C    7   COESCL      N*MAXT        1+N*MAXT+4*N**2+5*N+N*(2*N+1)
C
C*****************************************************************
C VARIABLES THAT ARE PASSED IN ARRAY IPAR.
C
C    VARIABLE NAME       LENGTH               OFFSET
C
C    1   NNUMT             N                  1
C    2   KKDEG             N*(N+1)*MAXT      1+N
C
C*****************************************************************
C
        CALL SCLGNP(N,MAXT,NUMT,KDEG,0,ZERO,COEF,
     &    IPAR(1:N),
     &    IPAR(1+N:N+N*(N+1)*MAXT),
     &    PAR(1:N*MAXT),
     &    PAR(1+N*MAXT:N*MAXT+4*N**2),
     &    PAR(1+N*MAXT+4*N**2:N*MAXT+4*N**2+2*N),
     &    PAR(1+N*MAXT+4*N**2+2*N:N*MAXT+4*N**2+2*N+N*(2*N+1)),
     &    PAR(1+N*MAXT+4*N**2+2*N+N*(2*N+1):
     &        N*MAXT+4*N**2+4*N+N*(2*N+1)),
     &    FACV,
     &    PAR(1+N*MAXT+4*N**2+4*N+N*(2*N+1):
     &        N*MAXT+4*N**2+5*N+N*(2*N+1)),
     &    PAR(1+N*MAXT+4*N**2+5*N+N*(2*N+1):
     &        2*N*MAXT+4*N**2+5*N+N*(2*N+1)),
     &    IERR)
C
C       SET COEF EQUAL TO THE SCALED COEFFICIENTS
C
        IF (IERR .EQ. 0) THEN
          COEF(:,1:MAXT) = RESHAPE(PAR(1+N*MAXT+4*N**2+5*N+N*(2*N+1):
     &      2*N*MAXT+4*N**2+5*N+N*(2*N+1)), (/N,MAXT/) )
        END IF
      END IF
C
      P(1, 1)= .12324754231_R8
          P(2, 1)= .76253746298_R8
      P(1, 2)= .93857838950_R8
          P(2, 2)=-.99375892810_R8
      P(1, 3)=-.23467908356_R8
          P(2, 3)= .39383930009_R8
      P(1, 4)= .83542556622_R8
          P(2, 4)=-.10192888288_R8
      P(1, 5)=-.55763522521_R8
          P(2, 5)=-.83729899911_R8
      P(1, 6)=-.78348738738_R8
          P(2, 6)=-.10578234903_R8
      P(1, 7)= .03938347346_R8
          P(2, 7)= .04825184716_R8
      P(1, 8)=-.43428734331_R8
          P(2, 8)= .93836289418_R8
      P(1, 9)=-.99383729993_R8
          P(2, 9)=-.40947822291_R8
      P(1,10)= .09383736736_R8
          P(2,10)= .26459172298_R8
C
      Q(1, 1)= .58720452864_R8
          Q(2, 1)= .01321964722_R8
      Q(1, 2)= .97884134700_R8
          Q(2, 2)=-.14433009712_R8
      Q(1, 3)= .39383737289_R8
          Q(2, 3)= .41543223411_R8
      Q(1, 4)=-.03938376373_R8
          Q(2, 4)=-.61253112318_R8
      Q(1, 5)= .39383737388_R8
          Q(2, 5)=-.26454678861_R8
      Q(1, 6)=-.00938376766_R8
          Q(2, 6)= .34447867861_R8
      Q(1, 7)=-.04837366632_R8
          Q(2, 7)= .48252736790_R8
      Q(1, 8)= .93725237347_R8
          Q(2, 8)=-.54356527623_R8
      Q(1, 9)= .39373957747_R8
          Q(2, 9)= .65573434564_R8
      Q(1,10)=-.39380038371_R8
          Q(2,10)= .98903450052_R8
C
      CCL(1, 1)=-.03485644332_R8
          CCL(2, 1)= .28554634336_R8
      CCL(1, 2)= .91453454766_R8
          CCL(2, 2)= .35354566613_R8
      CCL(1, 3)=-.36568737635_R8
          CCL(2, 3)= .45634642477_R8
      CCL(1, 4)=-.89089767544_R8
          CCL(2, 4)= .34524523544_R8
      CCL(1, 5)= .13523462465_R8
          CCL(2, 5)= .43534535555_R8
      CCL(1, 6)=-.34523544445_R8
          CCL(2, 6)= .00734522256_R8
      CCL(1, 7)=-.80004678763_R8
          CCL(2, 7)=-.009387123644_R8
      CCL(1, 8)=-.875432124245_R8
          CCL(2, 8)= .00045687651_R8
      CCL(1, 9)= .65256352333_R8
          CCL(2, 9)=-.12356777452_R8
      CCL(1,10)= .09986798321548_R8
          CCL(2,10)=-.56753456577_R8
      CCL(1,11)= .29674947394739_R8
          CCL(2,11)= .93274302173_R8
C
C IF THE PROJECTIVE TRANSFORMATION IS TO BE USED, THEN  CL  IS
C SET EQUAL TO THE  CCL  VALUES.  OTHERWISE,  CL  IS SET
C EQUAL TO NOMINAL VALUES.
C
      IF (IFLG1 .EQ. 01  .OR.  IFLG1 .EQ. 00) THEN 
        CL(1:2,1:N)=0.0
        CL(1,NP1)=1.0
        CL(2,NP1)=0.0
      ELSE
        DO J=1,NP1
          JJ=MOD(J-1,11)+1
          CL(1:2,J)=CCL(1:2,JJ)
        END DO
      END IF
C
C COMPUTE POWERS OF P AND Q, AND R=Q/P
      DO J=1,N
        JJ=MOD(J-1,10)+1
        CALL POWP(IDEG(J),P(1,JJ),PDG(1,J))
        CALL POWP(IDEG(J),Q(1,JJ),QDG(1,J))
        CALL DIVP(Q(1,JJ),P(1,JJ),R(1,J),IERR)
      END DO
      RETURN
      END SUBROUTINE INITP
      SUBROUTINE MULP(XXXX,YYYY,ZZZZ)
C
C THIS SUBROUTINE PERFORMS MULTIPLICATION OF COMPLEX NUMBERS:
C ZZZZ = XXXX*YYYY
C
C NOTE:  IN THE CALLING ROUTINE, ZZZZ SHOULD NOT BE THE SAME
C AS XXXX OR YYYY.  HOWEVER, XXXX MAY BE THE SAME AS YYYY.
C THUS, "CALL MULP(X,X,Z)" IS OK, BUT "CALL MULP(X,Y,X)" IS NOT.
C
C ON INPUT:
C
C XXXX  IS AN ARRAY OF LENGTH TWO REPRESENTING THE FIRST COMPLEX
C       NUMBER, WHERE XXXX(1) = REAL PART OF XXXX AND XXXX(2) =
C       IMAGINARY PART OF XXXX.
C
C YYYY  IS AN ARRAY OF LENGTH TWO REPRESENTING THE SECOND COMPLEX
C       NUMBER, WHERE YYYY(1) = REAL PART OF YYYY AND YYYY(2) =
C       IMAGINARY PART OF YYYY.
C
C ON OUTPUT:
C
C ZZZZ  IS AN ARRAY OF LENGTH TWO REPRESENTING THE RESULT OF
C       THE MULTIPLICATION, ZZZZ = XXXX*YYYY, WHERE ZZZZ(1) =
C       REAL PART OF ZZZZ AND ZZZZ(2) = IMAGINARY PART OF ZZZZ.
C
C DECLARATION OF INPUT
      USE REAL_PRECISION
      REAL (KIND=R8), DIMENSION(2), INTENT(IN):: XXXX,YYYY
C
C DECLARATION OF OUTPUT
      REAL (KIND=R8), DIMENSION(2), INTENT(OUT):: ZZZZ
C
      ZZZZ(1) = XXXX(1)*YYYY(1) - XXXX(2)*YYYY(2)
      ZZZZ(2) = XXXX(1)*YYYY(2) + XXXX(2)*YYYY(1)
      RETURN
      END SUBROUTINE MULP
      SUBROUTINE MULT2DS(Y, B, X, ROWPOS, COLPOS, N, LENB)
C
C     Returns B*X as Y.
C---------------------------------------------------------------------
C
C     Input variables:
C
C       B       matrix stored in the sparse row scheme.
C       ROWPOSP  indices of row-starts within B.
C       COLPOSP  column indices for matrix B stored in the general
C                sparse row scheme.
C       N       logical row dimension of B
C       LENB    number of data entries in B.
C       X       source vector -- should be compatible with B.
C       Y       target vector -- should have dimension N.
C
C     Output variables:
C
C       Y       value of B*X.
C
C---------------------------------------------------------------------
C
      USE REAL_PRECISION
      INTEGER, INTENT (IN):: LENB, N, ROWPOS(N+1), COLPOS(LENB)
      REAL (KIND=R8), INTENT(IN):: X(:), B(LENB)
      REAL (KIND=R8), INTENT (OUT) :: Y(N)
C
C LOCAL VARIABLES.
C
      INTEGER:: I, FIN, STRT, K
      REAL (KIND=R8):: TMP
C
        DO I = 1, N
         STRT = ROWPOS(I)
         FIN = ROWPOS(I+1)-1
         TMP = 0.0
         DO K = STRT, FIN
            TMP = TMP + B(K)*X(COLPOS(K))
         END DO 
         Y(I) = TMP
        END DO
      RETURN
      END SUBROUTINE MULT2DS
      SUBROUTINE MULTDS(Y,AA,X,MAXA,NN,LENAA)
C
C     This subroutine accepts a matrix, AA, in packed skyline storage form and
C       a vector, x, and returns the product AA*x in y.
C
C     Input Variables:
C
C       AA -- one dimensional real array containing the NN x NN matrix in 
C             packed skyline storage form.
C
C       x -- real vector of length NN to be multiplied by AA.
C
C       MAXA -- integer array used for specifying information about AA.
C               MAXA has length NN+1, and stores the indices of the 
C               diagonal elements of the matrix packed in AA.  By 
C               convention, MAXA(NN+1) = LENAA + 1 .
C
C       NN -- dimension of the matrix packed in AA .
C
C       LENAA -- number of elements in AA.
C
C
C     Output Variables:
C
C       y -- real vector of length NN containing the product  AA*x .
C
C
C
      USE REAL_PRECISION
      INTEGER, INTENT(IN):: LENAA,NN,MAXA(NN+1)
      REAL (KIND=R8), INTENT(IN):: AA(LENAA),X(NN)
      REAL (KIND=R8), INTENT(OUT):: Y(NN)
      INTEGER:: I,II,KK,KL,KU
      REAL (KIND=R8):: B,CC
      IF (LENAA .LE. NN) THEN
        DO I=1,NN
          Y(I)=AA(I)*X(I)
        END DO
        RETURN
      END IF
      DO I=1,NN
        Y(I)=0.00
      END DO
      DO I=1,NN
        KL=MAXA(I)
        KU=MAXA(I+1)-1
        II=I+1
        CC=X(I)
        DO KK=KL,KU
          II=II-1
          Y(II)=Y(II)+AA(KK)*CC
        END DO
      END DO
      IF (NN .EQ. 1) RETURN
      DO I=2,NN
        KL=MAXA(I)+1
        KU=MAXA(I+1)-1
        IF (KU-KL .LT. 0) CYCLE
        II=I
        B=0.00
        DO KK=KL,KU
          II=II-1
          B=B+AA(KK)*X(II)
        END DO
        Y(I)=Y(I)+B
      END DO
      RETURN
      END SUBROUTINE MULTDS
      SUBROUTINE OTPUTP(N,NUMPAT,CL,FACV,CLX,X,XNP1)
C
C OTPUTP  POSTPROCESSES THE ENDPOINTS OF THE PATHS, UNTRANSFORMING 
C AND UNSCALING THEM.
C
C ON INPUT:
C
C N  IS THE NUMBER OF EQUATIONS AND VARIABLES.
C
C NUMPAT  IS THE CURRENT PATH NUMBER.         
C
C CL  IS THE ARRAY THAT DEFINES THE PROJECTIVE TRANSFORMATION.
C
C FACV  CONTAINS THE VARIABLE SCALING FACTORS.
C
C X  IS THE ENDPOINT OF THE PATH, POSSIBLY TRANSFORMED AND/OR SCALED 
C   DEPENDING ON THE  POLSYS  INPUT FLAG  IFLG1.
C
C CLX  IS WORKSPACE.
C
C ON OUTPUT:
C
C N, NUMPAT, CL, AND  FACV  ARE UNCHANGED.
C
C X  IS THE UNTRANSFORMED AND UNSCALED VERSION OF X.
C
C XNP1  IS THE PROJECTIVE COORDINATE "X(N+1)".  XNP1  EQUALS UNITY IF
C   THE PROJECTIVE TRANSFORMATION IS NOT SPECIFIED.
C
      USE REAL_PRECISION
C
      INTERFACE
        SUBROUTINE DIVP(XXXX,YYYY,ZZZZ,IERR)
        USE REAL_PRECISION
        REAL (KIND=R8), DIMENSION(2), INTENT(IN):: XXXX,YYYY
        INTEGER, INTENT(OUT):: IERR
        REAL (KIND=R8), DIMENSION(2), INTENT(OUT):: ZZZZ
        END SUBROUTINE DIVP
        SUBROUTINE MULP(XXXX,YYYY,ZZZZ)
        USE REAL_PRECISION
        REAL (KIND=R8), DIMENSION(2), INTENT(IN):: XXXX,YYYY
        REAL (KIND=R8), DIMENSION(2), INTENT(OUT):: ZZZZ
        END SUBROUTINE MULP
        SUBROUTINE POWP(NNNN,XXXX,YYYY)
        USE REAL_PRECISION
        INTEGER, INTENT(IN):: NNNN
        REAL (KIND=R8), DIMENSION(2), INTENT(IN):: XXXX
        REAL (KIND=R8), DIMENSION(2), INTENT(IN OUT):: YYYY
        END SUBROUTINE POWP
      END INTERFACE
C
C DECLARATIONS OF INPUT, WORKSPACE, AND OUTPUT:
      INTEGER, INTENT(IN):: N,NUMPAT
      REAL (KIND=R8), INTENT(IN):: CL(2,N+1),FACV(N)
      REAL (KIND=R8), INTENT(IN OUT):: CLX(2,N),X(2,N),XNP1(2)
C
C DECLARATION OF LOCAL VARIABLES
      INTEGER:: I,IERR,J,NP1 
      REAL (KIND=R8):: FAC,TEMP(2)
C
      NP1=N+1
C COMPUTE XNP1
      DO J=1,N
        CALL MULP(CL(1,J),X(1,J),CLX(1,J))
      END DO
      XNP1 = CL(:,NP1) + SUM(CLX,DIM=2)
C UNTRANSFORM VARIABLES
      DO J=1,N
        CALL DIVP(X(1,J),XNP1,TEMP,IERR)
        X(1,J)=TEMP(1)
        X(2,J)=TEMP(2)
      END DO
C UNSCALE VARIABLES
      TEMP(1) = HUGE(1.0_R8)
      DO J=1,N
        FAC=10.**FACV(J)
        DO I=1,2
          IF( (ABS(X(I,J))/TEMP(1))*FAC .LT. 1.0 ) X(I,J)=FAC*X(I,J)
        END DO
      END DO
      RETURN
      END SUBROUTINE OTPUTP
      SUBROUTINE PCGDS(NN,LENAA,IFLAG,START,RHS)
C
C     PCGDS computes a tangent vector or normal flow correction using
C     a preconditioned conjugate gradient method (adaptive GMRES(k)).
C
C     The system to be solved is in the form Bx=b, where
C
C        +--          --+        +-   -+    +-    -+
C        |        |     |        |  0  |    |      | 
C        |   AA   | -PP |        | ... |    | -RHS |
C    B = |        |     | ,  b = |  0  | or |      |, 
C        +--------+-----+        +-----+    +------+ 
C        |    E(k)**t   |        |  T  |    |  T   | 
C        +--          --+        +-   -+    +-    -+
C
C        T = START(k), where |START(k)|=     max    |START(i)|.
C                                        1<=i<=NN+1
C                           
C        AA is an (NN x NN) symmetric matrix, PP is an (NN x 1) vector,
C        b is of length NN+1 and E(k)**t is the ( 1 x (NN+1) ) vector
C        consisting of all zeros, except for a '1' in its k-th position.
C        It is assumed that rank [AA,-PP]=NN and B is invertible.
C
C   The system is solved by splitting B into two matrices M and L, where
C
C       +-        -+                                +-     -+
C       |      |   |                                |       |
C       |  AA  | c |                                | -PP-c |
C   M = |      |   |  ,  L = u * [E(NN+1)**t],  u = |       | ,
C       +------+---+                                +-------+
C       |  c   | d |                                |  d'   |
C       +-        -+                                +-     -+
C
C   d = 1 and d' = 0 if k = NN+1, otherwise d = -d' = 1 + 1/M(k,k).
C   E(NN+1) is the (NN+1) x 1 vector consisting of all zeros except for
C   a '1' in its last position, and x**t is the transpose of x.
C
C    The final solution vector, x, is given by
C
C            +-                                    -+
C            |           [sol(u)]*[E(NN+1)**t]      |
C       x =  | I  -  -----------------------------  | * sol(b)
C            |        {[(sol(u))**t]*E(NN+1)}+1.0   |
C            +-                                    -+
C
C     where sol(a)=[M**(-1)]*a.  The two systems (Mz=u, Mz=b) are solved
C     by a preconditioned GMRES algorithm.
C
C  Input variables:
C
C        NN -- dimension of the matrix packed in AA.
C
C        LENAA -- number of elements in the packed array AA.
C
C        START -- vector of length NN+1, normally the solution to the
C                 previous linear system; used to determine the index k .
C
C        RHS -- optional vector of length NN, used to define right hand
C               side for normal flow correction calculation.  It is
C               assumed that PCGDS is called without RHS present before
C               it is called with RHS present.  A Gill-Murray LL^t
C               factorization based preconditioner is computed only when
C               RHS is not present.
C
C  Input variables defined in module HOMOTOPY:
C
C        AA -- one dimensional real array containing the leading NN x NN
C              submatrix of B in packed skyline storage form.
C
C        ROWPOS -- integer array used for specifying information about AA.
C                  Using packed skyline storage, it has length NN+2, and
C                  stores the indices of the diagonal elements within AA.
C                  ROWPOS(NN+1) = LENAA + 1 and ROWPOS(NN+2) = 
C                  LENAA + NN + 3 - k (k as defined above) by convention.
C                  (NOTE:  The value of ROWPOS(NN+2) is set by this
C                  subroutine when the preconditioning matrix Q is
C                  initialized.)
C
C                For example, using the packed storage scheme,
C                a symmetric 5 x 5 matrix of the form
C
C                +--             --+
C                |  1  3  0  0  0  |
C                |  3  2  0  7  0  |
C                |  0  0  4  6  0  |
C                |  0  7  6  5  9  |
C                |  0  0  0  9  8  |
C                +--             --+
C
C                would result in NN=5, LENAA=9, ROWPOS=(1,2,4,5,8,10,*),
C                and AA=(1,2,3,4,5,6,7,8,9).
C
C        PP -- vector of length NN, used for (NN+1)st column of
C              augmented matrix B .
C
C  Output variables:
C
C        START -- solution vector x of  B x = b  (defined above).
C
C        IFLAG -- normally unchanged on output.  If the GMRES
C                 iteration fails to converge in 10*(NN+1) iterations (most
C                 likely due to a singular Jacobian matrix), PCGDS returns
C                 with  IFLAG = 4 , and does not compute x.
C
C    Calls subroutines GMFADS and GMRES.
C
      USE HOMOTOPY, AA => QRSPARSE, WORK => PAR
      USE REAL_PRECISION
      INTEGER, INTENT(IN):: LENAA, NN
      INTEGER, INTENT(IN OUT):: IFLAG
      REAL (KIND=R8), INTENT(IN OUT):: START(NN+1)
      REAL (KIND=R8), INTENT(IN), OPTIONAL:: RHS(NN)
C
C LOCAL VARIABLES.
C
      INTEGER:: ITMAX, K, KD, NP1, QIND, ZBIND, ZUIND
      REAL (KIND=R8):: STARTK
      REAL (KIND=R8):: RHSC(NN+1)            ! RIGHT-HAND SIDE FOR GMRES.
C
C GMRES PARAMETERS.
C 
      INTEGER, PARAMETER:: SUBSPACE=8        ! KRYLOV SUBSPACE VALUE.
C
      INTERFACE 
        SUBROUTINE GMFADS(NN,A,NWK,MAXA)
          USE REAL_PRECISION
          INTEGER, INTENT(IN):: NN,NWK,MAXA(NN+1)
          REAL (KIND=R8), INTENT(IN OUT):: A(NWK)
        END SUBROUTINE GMFADS
        SUBROUTINE GMRES(N, KDMAX, ITMAX, RHSC, X, KVAL,
     &                Q, IFLAG, ROWPOSP, COLPOSP)
          USE HOMOTOPY
          USE REAL_PRECISION
          INTEGER, INTENT(IN):: KVAL, N
          INTEGER, INTENT(IN OUT):: IFLAG, ITMAX, KDMAX
          REAL (KIND=R8), INTENT(IN):: Q(:), RHSC(N)
          REAL (KIND=R8), INTENT(IN OUT):: X(N)
          INTEGER, INTENT(IN), OPTIONAL:: COLPOSP(:), ROWPOSP(:)
        END SUBROUTINE GMRES
      END INTERFACE
C    
      NP1 = NN+1
C
C     INITIALIZE START POSITIONS WITHIN WORK.
C
      ZBIND = 1
      ZUIND = NP1+1
      QIND = 2*NP1+1
C
      IF (.NOT. ALLOCATED(WORK)) THEN
        ALLOCATE(WORK(2*NP1+LENAA+NN+1))
        WORK(1:2*NP1) = 0.0
      END IF
C
C     FIND THE ELEMENT OF LARGEST MAGNITUDE IN THE INITIAL VECTOR, AND
C     RECORD ITS POSITION IN K.
C 
      K = MAXVAL(MAXLOC(ABS(START)))
      STARTK = START(K)
C
C     SET VALUES OF ROWPOS(NN+1) AND ROWPOS(NN+2), AND
C     COMPUTE THE PRECONDITIONER Q FOR M.
C
      IF (.NOT. PRESENT(RHS)) THEN
        WORK(QIND:QIND+LENAA-1) = AA(1:LENAA)
        ROWPOS(NP1) = LENAA+1
        ROWPOS(NN+2) = LENAA+NN+3-K
        WORK(QIND+LENAA+NN+1-K) = 1.0
        IF (K .LT. NP1) THEN
          WORK(QIND+LENAA) = 1.0 + 1.0/ABS(AA(ROWPOS(K)))
          WORK(QIND+LENAA+1:QIND+LENAA+NN-K) = 0.0
        END IF
        CALL GMFADS(NP1, WORK(QIND:QIND+LENAA+NP1-K),
     &            ROWPOS(NN+2)-1, ROWPOS(1:NN+2))
      END IF
C
C     COMPUTE RIGHT HAND SIDE FOR MZ=B.
C
      RHSC(NP1) = STARTK
      IF (PRESENT(RHS))  THEN
        RHSC(1:NN) = -RHS
      ELSE
        RHSC(1:NN) = 0.0
      END IF
C
C CALL TO GMRES (MZ=B SYSTEM). 
C
      ITMAX = 30*NP1
      KD = SUBSPACE          
      CALL GMRES(NP1, KD, ITMAX, RHSC, WORK(ZBIND:ZBIND+NN),  
     &           K, WORK(QIND:QIND+LENAA+NP1-K), IFLAG) 
      IF ( IFLAG .GT. 0) RETURN
C
C COMPUTE RIGHT HAND SIDE FOR  MZ=U.
C
      RHSC(1:NN) = -PP(1:NN)
      IF (K .LT. NP1) THEN 
        RHSC(K) = RHSC(K)-1.0
        RHSC(NP1) = -(1.0+1.0/ABS(AA(ROWPOS(K))))
      ELSE
        RHSC(NP1) = 0.0
      END IF
C
C CALL TO GMRES (MZ=U SYSTEM).      
C
      ITMAX = 30*NP1-ITMAX
      KD = SUBSPACE
      CALL GMRES(NP1, KD, ITMAX, RHSC, WORK(ZUIND:ZUIND+NN),  
     &           K, WORK(QIND:QIND+LENAA+NP1-K), IFLAG) 
      IF ( IFLAG .GT. 0) RETURN
C
C COMPUTE THE FINAL SOLUTION BY SHERMAN-MORRISON FORMULA.
C
      STARTK = -WORK(ZBIND+NN)/(1.0+WORK(ZUIND+NN))
      START(1:NP1) = WORK(ZBIND:ZBIND+NN) + STARTK*WORK(ZUIND:ZUIND+NN)
      RETURN
      END SUBROUTINE PCGDS
      SUBROUTINE POWP(NNNN,XXXX,YYYY)
C
C THIS SUBROUTINE TAKES A NON-NEGATIVE POWER OF A COMPLEX NUMBER:
C YYYY = XXXX**NNNN USING DE MOIVRE'S FORMULA:
C
C     YYYY = R**NNNN * (COS(NNNN*THETA),SIN(NNNN*THETA)),
C
C WHERE R=DNRM2(2,XXXX,1) AND THETA=ATAN2(XXXX(2),XXXX(1)).
C
C NOTE: POWP SETS 0**0 EQUAL TO 1.
C
C ON INPUT:
C
C NNNN  IS A NON-NEGATIVE INTEGER.
C
C XXXX  IS AN ARRAY OF LENGTH TWO REPRESENTING A COMPLEX
C       NUMBER, WHERE XXXX(1) = REAL PART OF XXXX AND XXXX(2) =
C       IMAGINARY PART OF XXXX.
C
C ON OUTPUT:
C
C YYYY  IS AN ARRAY OF LENGTH TWO REPRESENTING THE RESULT OF
C       THE POWER, YYYY = XXXX**NNNN, WHERE YYYY(1) =
C       REAL PART OF YYYY AND YYYY(2) = IMAGINARY PART OF YYYY.
C
C SUBROUTINES: COS, SIN, ATAN2, DNRM2
C
C DECLARATION OF INPUT
      USE REAL_PRECISION
      INTEGER, INTENT(IN):: NNNN
      REAL (KIND=R8), DIMENSION(2), INTENT(IN):: XXXX
C
C DECLARATION OF OUTPUT
      REAL (KIND=R8), DIMENSION(2), INTENT(IN OUT):: YYYY
C
C DECLARATION OF VARIABLES
      REAL (KIND=R8):: R,RR,T,TT
C
C DECLARATION OF FUNCTIONS
      REAL (KIND=R8)::  DNRM2
C
      IF (NNNN .EQ. 0) THEN
          YYYY(1)=1.
          YYYY(2)=0.
          RETURN
      ENDIF
      IF (NNNN .EQ. 1) THEN
          YYYY(1)=XXXX(1)
          YYYY(2)=XXXX(2)
          RETURN
      ENDIF
      R = DNRM2(2,XXXX,1)
      IF (R .EQ. 0.0) THEN
          YYYY(1)=0.0
          YYYY(2)=0.0
          RETURN
      END IF
      RR= R**NNNN
      T = ATAN2(XXXX(2),XXXX(1))
      TT= NNNN*T
      YYYY(1) = RR*COS(TT)
      YYYY(2) = RR*SIN(TT)
      RETURN
      END SUBROUTINE POWP
      SUBROUTINE ROOT(T,FT,B,C,RELERR,ABSERR,IFLAG)
C
C  ROOT COMPUTES A ROOT OF THE NONLINEAR EQUATION F(X)=0
C  WHERE F(X) IS A CONTINOUS REAL FUNCTION OF A SINGLE REAL
C  VARIABLE X.  THE METHOD USED IS A COMBINATION OF BISECTION
C  AND THE SECANT RULE.
C
C  NORMAL INPUT CONSISTS OF A CONTINUOS FUNCTION F AND AN
C  INTERVAL (B,C) SUCH THAT F(B)*F(C).LE.0.0.  EACH ITERATION
C  FINDS NEW VALUES OF B AND C SUCH THAT THE INTERVAL(B,C) IS
C  SHRUNK AND F(B)*F(C).LE.0.0.  THE STOPPING CRITERION IS
C
C          DABS(B-C).LE.2.0*(RELERR*DABS(B)+ABSERR)
C
C  WHERE RELERR=RELATIVE ERROR AND ABSERR=ABSOLUTE ERROR ARE
C  INPUT QUANTITIES.  SET THE FLAG, IFLAG, POSITIVE TO INITIALIZE
C  THE COMPUTATION.  AS B,C AND IFLAG ARE USED FOR BOTH INPUT AND
C  OUTPUT, THEY MUST BE VARIABLES IN THE CALLING PROGRAM.
C
C  IF 0 IS A POSSIBLE ROOT, ONE SHOULD NOT CHOOSE ABSERR=0.0.
C
C  THE OUTPUT VALUE OF B IS THE BETTER APPROXIMATION TO A ROOT
C  AS B AND C ARE ALWAYS REDEFINED SO THAT DABS(F(B)).LE.DABS(F(C)).
C
C  TO SOLVE THE EQUATION, ROOT MUST EVALUATE F(X) REPEATEDLY. THIS
C  IS DONE IN THE CALLING PROGRAM.  WHEN AN EVALUATION OF F IS
C  NEEDED AT T, ROOT RETURNS WITH IFLAG NEGATIVE.  EVALUATE FT=F(T)
C  AND CALL ROOT AGAIN.  DO NOT ALTER IFLAG.
C
C  WHEN THE COMPUTATION IS COMPLETE, ROOT RETURNS TO THE CALLING
C  PROGRAM WITH IFLAG POSITIVE=
C
C     IFLAG=1  IF F(B)*F(C).LT.0 AND THE STOPPING CRITERION IS MET.
C
C          =2  IF A VALUE B IS FOUND SUCH THAT THE COMPUTED VALUE
C              F(B) IS EXACTLY ZERO.  THE INTERVAL (B,C) MAY NOT
C              SATISFY THE STOPPING CRITERION.
C
C          =3  IF DABS(F(B)) EXCEEDS THE INPUT VALUES DABS(F(B)),
C              DABS(F(C)).  IN THIS CASE IT IS LIKELY THAT B IS CLOSE
C              TO A POLE OF F.
C
C          =4  IF NO ODD ORDER ROOT WAS FOUND IN THE INTERVAL.  A
C              LOCAL MINIMUM MAY HAVE BEEN OBTAINED.
C
C          =5  IF TOO MANY FUNCTION EVALUATIONS WERE MADE.
C              (AS PROGRAMMED, 500 ARE ALLOWED.)
C
C  THIS CODE IS A MODIFICATION OF THE CODE ZEROIN WHICH IS COMPLETELY
C  EXPLAINED AND DOCUMENTED IN THE TEXT  NUMERICAL COMPUTING:  AN
C  INTRODUCTION,  BY L. F. SHAMPINE AND R. C. ALLEN.
C
C
      USE REAL_PRECISION
      REAL (KIND=R8):: A,ABSERR,ACBS,ACMB,AE,B,C,CMB,FA,FB,
     &  FC,FT,FX,P,Q,RE,RELERR,T,TOL,U
      INTEGER IC,IFLAG,KOUNT
      SAVE
C
      IF(IFLAG.GE.0) GO TO 100
      IFLAG=ABS(IFLAG)
      GO TO (200,300,400), IFLAG
  100 U=EPSILON(1.0_R8)
      RE=MAX(RELERR,U)
      AE=MAX(ABSERR,0.0_R8)
      IC=0
      ACBS=ABS(B-C)
      A=C
      T=A
      IFLAG=-1
      RETURN
  200 FA=FT
      T=B
      IFLAG=-2
      RETURN
  300 FB=FT
      FC=FA
      KOUNT=2
      FX=MAX(ABS(FB),ABS(FC))
    1 IF(ABS(FC).GE.ABS(FB))GO TO 2
C
C  INTERCHANGE B AND C SO THAT ABS(F(B)).LE.ABS(F(C)).
C
      A=B
      FA=FB
      B=C
      FB=FC
      C=A
      FC=FA
    2 CMB=0.5*(C-B)
      ACMB=ABS(CMB)
      TOL=RE*ABS(B)+AE
C
C  TEST STOPPING CRITERION AND FUNCTION COUNT.
C
      IF(ACMB.LE.TOL)GO TO 8
      IF(KOUNT.GE.500)GO TO 12
C
C  CALCULATE NEW ITERATE EXPLICITLY AS B+P/Q
C  WHERE WE ARRANGE P.GE.0.  THE IMPLICIT
C  FORM IS USED TO PREVENT OVERFLOW.
C
      P=(B-A)*FB
      Q=FA-FB
      IF(P.GE.0.0)GO TO 3
      P=-P
      Q=-Q
C
C  UPDATE A, CHECK IF REDUCTION IN THE SIZE OF BRACKETING
C  INTERVAL IS SATISFACTORY.  IF NOT BISECT UNTIL IT IS.
C
    3 A=B
      FA=FB
      IC=IC+1
      IF(IC.LT.4)GO TO 4
      IF(8.0*ACMB.GE.ACBS)GO TO 6
      IC=0
      ACBS=ACMB
C
C  TEST FOR TOO SMALL A CHANGE.
C
    4 IF(P.GT.ABS(Q)*TOL)GO TO 5
C
C  INCREMENT BY TOLERANCE
C
      B=B+SIGN(TOL,CMB)
      GO TO 7
C
C  ROOT OUGHT TO BE BETWEEN B AND (C+B)/2
C
    5 IF(P.GE.CMB*Q)GO TO 6
C
C  USE SECANT RULE.
C
      B=B+P/Q
      GO TO 7
C
C  USE BISECTION.
C
    6 B=0.5*(C+B)
C
C  HAVE COMPLETED COMPUTATION FOR NEW ITERATE B.
C
    7 T=B
      IFLAG=-3
      RETURN
  400 FB=FT
      IF(FB.EQ.0.0)GO TO 9
      KOUNT=KOUNT+1
      IF(SIGN(1.0_R8,FB).NE.SIGN(1.0_R8,FC))GO TO 1
      C=A
      FC=FA
      GO TO 1
C
C FINISHED.  SET IFLAG.
C
    8 IF(SIGN(1.0_R8,FB).EQ.SIGN(1.0_R8,FC))GO TO 11
      IF(ABS(FB).GT.FX)GO TO 10
      IFLAG=1
      RETURN
    9 IFLAG=2
      RETURN
   10 IFLAG=3
      RETURN
   11 IFLAG=4
      RETURN
   12 IFLAG=5
      RETURN
      END SUBROUTINE ROOT
      SUBROUTINE ROOTNF(N,NFE,IFLAG,RELERR,ABSERR,Y,YP,YOLD,
     &   YPOLD,A,QR,ALPHA,TZ,PIVOT,W,WP)
C
C ROOTNF  FINDS THE POINT  YBAR = (1, XBAR)  ON THE ZERO CURVE OF THE
C HOMOTOPY MAP.  IT STARTS WITH TWO POINTS YOLD=(LAMBDAOLD,XOLD) AND
C Y=(LAMBDA,X) SUCH THAT  LAMBDAOLD < 1 <= LAMBDA , AND ALTERNATES
C BETWEEN SECANT ESTIMATES OF  YBAR  AND NEWTON ITERATION UNTIL
C CONVERGENCE.
C
C THE FOLLOWING INTERFACE BLOCK SHOULD BE INCLUDED IN THE CALLING
C PROGRAM:
C
C     INTERFACE
C       SUBROUTINE ROOTNF(N,NFE,IFLAG,RELERR,ABSERR,Y,YP,YOLD,
C    &    YPOLD,A,QR,ALPHA,TZ,PIVOT,W,WP)
C       USE REAL_PRECISION
C       REAL (KIND=R8):: ABSERR,RELERR
C       INTEGER:: IFLAG,N,NFE
C       REAL (KIND=R8):: A(:),ALPHA(3*N+3),QR(N,N+2),TZ(N+1),W(N+1),
C    &    WP(N+1),Y(:),YOLD(N+1),YP(N+1),YPOLD(N+1)
C       INTEGER:: PIVOT(N+1)
C       END SUBROUTINE ROOTNF
C     END INTERFACE
C
C
C ON INPUT:
C
C N = DIMENSION OF X AND THE HOMOTOPY MAP.
C
C NFE = NUMBER OF JACOBIAN MATRIX EVALUATIONS.
C
C IFLAG = -2, -1, OR 0, INDICATING THE PROBLEM TYPE.
C
C RELERR, ABSERR = RELATIVE AND ABSOLUTE ERROR VALUES.  THE ITERATION IS
C    CONSIDERED TO HAVE CONVERGED WHEN A POINT Y=(LAMBDA,X) IS FOUND 
C    SUCH THAT
C
C    |Y(1) - 1| <= RELERR + ABSERR              AND
C
C    ||Z|| <= RELERR*||X|| + ABSERR  ,          WHERE
C
C    (?,Z) IS THE NEWTON STEP TO Y=(LAMBDA,X).
C
C Y(1:N+1) = POINT (LAMBDA(S), X(S)) ON ZERO CURVE OF HOMOTOPY MAP.
C
C YP(1:N+1) = UNIT TANGENT VECTOR TO THE ZERO CURVE OF THE HOMOTOPY MAP
C    AT  Y .
C
C YOLD(1:N+1) = A POINT DIFFERENT FROM  Y  ON THE ZERO CURVE.
C
C YPOLD(1:N+1) = UNIT TANGENT VECTOR TO THE ZERO CURVE OF THE HOMOTOPY
C    MAP AT  YOLD .
C
C A(:) = PARAMETER VECTOR IN THE HOMOTOPY MAP.
C
C QR(1:N,1:N+2), ALPHA(1:3*N+3), TZ(1:N+1), PIVOT(1:N+1), W(1:N+1), 
C    WP(1:N+1)  ARE WORK ARRAYS USED FOR THE QR FACTORIZATION (IN THE
C    NEWTON STEP CALCULATION) AND THE INTERPOLATION.
C
C
C ON OUTPUT:
C
C N , RELERR , ABSERR , A  ARE UNCHANGED.
C
C NFE  HAS BEEN UPDATED.
C
C IFLAG 
C    = -2, -1, OR 0 (UNCHANGED) ON A NORMAL RETURN.
C
C    = 4 IF A JACOBIAN MATRIX WITH RANK < N HAS OCCURRED.  THE
C        ITERATION WAS NOT COMPLETED.
C
C    = 6 IF THE ITERATION FAILED TO CONVERGE.  Y  AND  YOLD  CONTAIN
C        THE LAST TWO POINTS FOUND ON THE ZERO CURVE.
C
C Y  IS THE POINT ON THE ZERO CURVE OF THE HOMOTOPY MAP AT  LAMBDA = 1 .
C
C
C CALLS  DNRM2 , ROOT , TANGNF .
C
      USE REAL_PRECISION
      REAL (KIND=R8):: ABSERR,AERR,
     &   DD001,DD0011,DD01,DD011,DELS,F0,F1,FP0,FP1,
     &   QOFS,QSOUT,RELERR,RERR,S,SA,SB,SOUT,U
      INTEGER:: IFLAG,JUDY,JW,LCODE,LIMIT,N,NFE,NP1
      LOGICAL:: BRACK
C
C ***** ARRAY DECLARATIONS. *****
C
      REAL (KIND=R8):: A(:),ALPHA(3*N+3),QR(N,N+2),TZ(N+1),W(N+1),
     &   WP(N+1),Y(:),YOLD(N+1),YP(N+1),YPOLD(N+1)
      INTEGER:: PIVOT(N+1)
C
C ***** END OF DIMENSIONAL INFORMATION. *****
C
      INTERFACE
        FUNCTION DNRM2(N,X,STRIDE)
        USE REAL_PRECISION
        INTEGER:: N,STRIDE
        REAL (KIND=R8):: DNRM2,X(N)
        END FUNCTION DNRM2
        SUBROUTINE TANGNF(RHOLEN,Y,YP,YPOLD,A,QR,ALPHA,TZ,PIVOT,
     &    NFE,N,IFLAG)
        USE REAL_PRECISION
        REAL (KIND=R8):: RHOLEN
        INTEGER:: IFLAG,N,NFE
        REAL (KIND=R8):: A(:),Y(:),YP(N+1),YPOLD(N+1)
        REAL (KIND=R8):: ALPHA(3*N+3),QR(N,N+2),TZ(N+1)
        INTEGER:: PIVOT(N+1)
        END SUBROUTINE TANGNF
      END INTERFACE
C
C DEFINITION OF HERMITE CUBIC INTERPOLANT VIA DIVIDED DIFFERENCES.
C
      DD01(F0,F1,DELS)=(F1-F0)/DELS
      DD001(F0,FP0,F1,DELS)=(DD01(F0,F1,DELS)-FP0)/DELS
      DD011(F0,F1,FP1,DELS)=(FP1-DD01(F0,F1,DELS))/DELS
      DD0011(F0,FP0,F1,FP1,DELS)=(DD011(F0,F1,FP1,DELS) - 
     &                            DD001(F0,FP0,F1,DELS))/DELS
      QOFS(F0,FP0,F1,FP1,DELS,S)=((DD0011(F0,FP0,F1,FP1,DELS)*(S-DELS) +
     &   DD001(F0,FP0,F1,DELS))*S + FP0)*S + F0
C
C
      U=EPSILON(1.0_R8)
      RERR=MAX(RELERR,U)
      AERR=MAX(ABSERR,0.0_R8)
      NP1=N+1
C
C THE LIMIT ON THE NUMBER OF ITERATIONS ALLOWED MAY BE CHANGED BY
C CHANGING THE FOLLOWING PARAMETER STATEMENT:
      LIMIT=2*(INT(ABS(LOG10(AERR+RERR)))+1)
C
      TZ=Y - YOLD
      DELS=DNRM2(NP1,TZ,1)
C
C USING TWO POINTS AND TANGENTS ON THE HOMOTOPY ZERO CURVE, CONSTRUCT 
C THE HERMITE CUBIC INTERPOLANT Q(S).  THEN USE  ROOT  TO FIND THE S
C CORRESPONDING TO  LAMBDA = 1 .  THE TWO POINTS ON THE ZERO CURVE ARE
C ALWAYS CHOSEN TO BRACKET LAMBDA=1, WITH THE BRACKETING INTERVAL
C ALWAYS BEING [0, DELS].
C
      SA=0.0
      SB=DELS
      LCODE=1
130   CALL ROOT(SOUT,QSOUT,SA,SB,RERR,AERR,LCODE)
      IF (LCODE .GT. 0) GO TO 140
      QSOUT=QOFS(YOLD(1),YPOLD(1),Y(1),YP(1),DELS,SOUT) - 1.0
      GO TO 130
C IF LAMBDA = 1 WERE BRACKETED,  ROOT  CANNOT FAIL.
140   IF (LCODE .GT. 2) THEN
        IFLAG=6
        RETURN
      ENDIF
C
C CALCULATE Q(SA) AS THE INITIAL POINT FOR A NEWTON ITERATION.
      DO JW=1,NP1
        W(JW)=QOFS(YOLD(JW),YPOLD(JW),Y(JW),YP(JW),DELS,SA)
      END DO
C
C ***** END OF CALCULATION OF CUBIC INTERPOLANT. *****
C
C TANGENT INFORMATION  YP  IS NO LONGER NEEDED.  HEREAFTER,  YP
C REPRESENTS THE MOST RECENT POINT WHICH IS ON THE OPPOSITE SIDE OF
C THE HYPERPLANE  LAMBDA = 1 FROM  Y.
C
C    PREPARE FOR MAIN LOOP.
C
      YP=YOLD
C
C INITIALIZE  BRACK  TO INDICATE THAT THE POINTS  Y  AND  YOLD  BRACKET
C LAMBDA = 1,  THUS  YOLD = YP .
C
      BRACK = .TRUE.
C
C ***** MAIN LOOP. *****
      DO JUDY = 1,LIMIT
C CALCULATE NEWTON STEP AT CURRENT ESTIMATE  W .
      CALL TANGNF(SA,W,WP,YPOLD,A,QR,ALPHA,TZ,PIVOT,NFE,N,IFLAG)
      IF (IFLAG .GT. 0) RETURN
C
C NEXT POINT = CURRENT POINT + NEWTON STEP.
C
      W = W + TZ
C
C CHECK FOR CONVERGENCE.
C
      IF ((ABS(W(1)-1.0) .LE. RERR+AERR) .AND.
     &    (DNRM2(NP1,TZ,1) .LE. RERR*DNRM2(N,W(2:NP1),1)+AERR)) THEN
        Y = W
        RETURN
      ENDIF
C
C PREPARE FOR NEXT ITERATION.
C
      IF (ABS(W(1)-1.0) .LE. RERR+AERR) THEN
         YPOLD=WP
         CYCLE
      ENDIF
C
C    UPDATE  Y  AND  YOLD .
C
      YOLD=Y
      Y=W
C
C    UPDATE  YP  SUCH THAT  YP  IS THE MOST RECENT POINT
C    OPPOSITE OF  LAMBDA = 1  FROM  Y .  SET  BRACK = .TRUE.  IFF
C    Y  AND  YOLD  BRACKET  LAMBDA = 1  SO THAT  YP = YOLD .
C
          IF ((Y(1)-1.0)*(YOLD(1)-1.0) .GT. 0) THEN
            BRACK = .FALSE.
          ELSE
            BRACK = .TRUE.
            YP=YOLD
          END IF
C
C    COMPUTE  DELS = ||Y-YP||.
C
          TZ=Y - YP
          DELS=DNRM2(NP1,TZ,1)
C
C       COMPUTE  TZ  FOR THE LINEAR PREDICTOR   W = Y + TZ,
C           WHERE  TZ = SA*(YOLD-Y).
C
          SA = (1.0-Y(1))/(YOLD(1)-Y(1))
          TZ = SA*(YOLD - Y)
C
C       TO INSURE STABILITY, THE LINEAR PREDICTION MUST BE NO FARTHER
C       FROM  Y  THAN  YP  IS.  THIS IS GUARANTEED IF  BRACK = .TRUE.
C       IF LINEAR PREDICTION IS TOO FAR AWAY, USE BRACKETING POINTS
C       TO COMPUTE LINEAR PREDICTION.
C
          IF (.NOT. BRACK) THEN
            IF (DNRM2(NP1,TZ,1) .GT. DELS) THEN
C
C             COMPUTE  TZ = SA*(YP-Y).
C
              SA = (1.0-Y(1))/(YP(1)-Y(1))
              TZ = SA*(YP - Y)
            END IF
          END IF
C
C       COMPUTE ESTIMATE  W = Y + TZ  AND SAVE OLD TANGENT VECTOR.
C
           W = W + TZ
           YPOLD = WP
       END DO
C
C ***** END OF MAIN LOOP. *****
C
C THE ALTERNATING SECANT ESTIMATION AND NEWTON ITERATION
C HAS NOT CONVERGED IN  LIMIT  STEPS.  ERROR RETURN.
      IFLAG=6
      RETURN
      END SUBROUTINE ROOTNF
      SUBROUTINE ROOTNS(N,NFE,IFLAG,RELERR,ABSERR,Y,YP,YOLD,YPOLD,
     &   A,MODE,LENQR)
C
C ROOTNS  FINDS THE POINT  YBAR = (XBAR, 1)  ON THE ZERO CURVE OF THE
C HOMOTOPY MAP.  IT STARTS WITH TWO POINTS YOLD=(XOLD,LAMBDAOLD) AND
C Y=(X,LAMBDA) SUCH THAT  LAMBDAOLD < 1 <= LAMBDA , AND ALTERNATES
C BETWEEN USING A SECANT METHOD TO FIND A PREDICTED POINT ON THE 
C HYPERPLANE  LAMBDA=1, AND TAKING A NEWTON STEP TO RETURN TO THE 
C ZERO CURVE OF THE HOMOTOPY MAP.
C
C THE FOLLOWING INTERFACE BLOCK SHOULD BE INCLUDED IN THE CALLING
C PROGRAM:
C
C     INTERFACE
C       SUBROUTINE ROOTNS(NC,NFEC,IFLAGC,ANSRE,ANSAE,Y,YP,YOLD,YPOLD,
C    &     A,MODE,LENQR)
C       USE REAL_PRECISION
C       INTEGER, INTENT(IN):: LENQR,MODE,NC
C       INTEGER, INTENT(IN OUT):: IFLAGC,NFEC
C       REAL (KIND=R8), INTENT(IN):: A(:)
C       REAL (KIND=R8), INTENT(IN):: ANSAE,ANSRE
C       REAL (KIND=R8), DIMENSION(:), INTENT(IN OUT):: Y,YOLD,YP,YPOLD
C       END SUBROUTINE ROOTNS
C     END INTERFACE
C
C
C ON INPUT:
C
C N = DIMENSION OF X AND THE HOMOTOPY MAP.
C
C NFE = NUMBER OF JACOBIAN MATRIX EVALUATIONS.
C
C IFLAG = -2, -1, OR 0, INDICATING THE PROBLEM TYPE.
C
C RELERR, ABSERR = RELATIVE AND ABSOLUTE ERROR VALUES.  THE ITERATION IS
C    CONSIDERED TO HAVE CONVERGED WHEN A POINT Y=(X,LAMBDA) IS FOUND 
C    SUCH THAT
C
C    |Y(NP1) - 1| <= RELERR + ABSERR              AND
C
C    ||Z|| <= RELERR*||X|| + ABSERR  ,          WHERE
C
C    (Z,?) IS THE NEWTON STEP TO Y=(X,LAMBDA).
C
C Y(1:N+1) = POINT (X(S), LAMBDA(S)) ON ZERO CURVE OF HOMOTOPY MAP.
C
C YP(1:N+1) = UNIT TANGENT VECTOR TO THE ZERO CURVE OF THE HOMOTOPY MAP
C    AT  Y .
C
C YOLD(1:N+1) = A POINT DIFFERENT FROM  Y  ON THE ZERO CURVE.
C
C YPOLD(1:N+1) = UNIT TANGENT VECTOR TO THE ZERO CURVE OF THE HOMOTOPY
C    MAP AT  YOLD .
C
C A(:) = PARAMETER VECTOR IN THE HOMOTOPY MAP.
C
C MODE = 1 IF THE JACOBIAN MATRIX IS SYMMETRIC AND STORED IN A PACKED
C          SKYLINE FORMAT;
C      = 2 IF THE JACOBIAN MATRIX IS STORED IN A SPARSE ROW FORMAT.
C
C LENQR  IS THE NUMBER OF NONZERO ENTRIES IN THE SPARSE JACOBIAN
C    MATRICES, USED TO DETERMINE THE SPARSE MATRIX DATA STRUCTURES.
C
C
C ON OUTPUT:
C
C N , RELERR , ABSERR , A  ARE UNCHANGED.
C
C NFE  HAS BEEN UPDATED.
C
C IFLAG 
C    = -2, -1, OR 0 (UNCHANGED) ON A NORMAL RETURN.
C
C    = 4 IF THE PRECONDITIONED CONJUGATE GRADIENT ITERATION FAILED TO
C        CONVERGE (MOST LIKELY DUE TO A JACOBIAN MATRIX WITH RANK < N).
C        THE ITERATION WAS NOT COMPLETED.
C
C    = 6 IF THE INTERPOLATION/NEWTON ITERATION FAILED TO CONVERGE.  
C        Y  AND  YOLD  CONTAIN THE LAST TWO POINTS FOUND ON THE 
C        ZERO CURVE.
C
C Y  IS THE POINT ON THE ZERO CURVE OF THE HOMOTOPY MAP AT  LAMBDA = 1 .
C
C
C CALLS  DNRM2 , ROOT , TANGNS .
C
      USE REAL_PRECISION
      INTEGER, INTENT(IN):: LENQR,MODE,N
      INTEGER, INTENT(IN OUT):: IFLAG,NFE
      REAL (KIND=R8), INTENT(IN):: A(:)
      REAL (KIND=R8), INTENT(IN):: ABSERR,RELERR
      REAL (KIND=R8), DIMENSION(:), INTENT(IN OUT):: Y,YOLD,YP,YPOLD
C
C ***** LOCAL VARIABLES. *****
C
      REAL (KIND=R8):: AERR,DD001,DD0011,DD01,DD011,DELS,
     &   F0,F1,FP0,FP1,QOFS,QSOUT,RERR,S,SA,SB,SOUT,U
      INTEGER:: JUDY,JW,LCODE,NP1
      LOGICAL:: BRACK
      INTERFACE
        FUNCTION DNRM2(N,X,STRIDE)
        USE REAL_PRECISION
        INTEGER:: N,STRIDE
        REAL (KIND=R8):: DNRM2,X(N)
        END FUNCTION DNRM2
      END INTERFACE
C
C ***** AUTOMATIC WORK ARRAYS. *****
C
      REAL (KIND=R8):: TZ(N+1),W(N+1),WP(N+1)
C
C THE LIMIT ON THE NUMBER OF ITERATIONS ALLOWED MAY BE CHANGED BY
C CHANGING THE FOLLOWING PARAMETER STATEMENT:
      INTEGER, PARAMETER:: LIMIT=20
C
      INTERFACE
        SUBROUTINE TANGNS(RHOLEN,Y,YP,TZ,YPOLD,A,MODE,LENQR,
     &    NFE,N,IFLAG)
        USE REAL_PRECISION
        REAL (KIND=R8), INTENT(IN), DIMENSION(:):: A,Y,YPOLD
        REAL (KIND=R8), INTENT(IN OUT):: RHOLEN
        REAL (KIND=R8), INTENT(OUT), DIMENSION(:):: TZ,YP
        INTEGER, INTENT(IN):: LENQR,MODE,N
        INTEGER, INTENT(IN OUT):: IFLAG,NFE
        END SUBROUTINE TANGNS
      END INTERFACE
C
C DEFINITION OF HERMITE CUBIC INTERPOLANT VIA DIVIDED DIFFERENCES.
C
      DD01(F0,F1,DELS)=(F1-F0)/DELS
      DD001(F0,FP0,F1,DELS)=(DD01(F0,F1,DELS)-FP0)/DELS
      DD011(F0,F1,FP1,DELS)=(FP1-DD01(F0,F1,DELS))/DELS
      DD0011(F0,FP0,F1,FP1,DELS)=(DD011(F0,F1,FP1,DELS) - 
     &                            DD001(F0,FP0,F1,DELS))/DELS
      QOFS(F0,FP0,F1,FP1,DELS,S)=((DD0011(F0,FP0,F1,FP1,DELS)*(S-DELS) +
     &   DD001(F0,FP0,F1,DELS))*S + FP0)*S + F0
C
C ***** END OF SPECIFICATION INFORMATION. *****
C
      U=EPSILON(1.0_R8)
      RERR=MAX(RELERR,U)
      AERR=MAX(ABSERR,0.0_R8)
      NP1=N+1
      TZ=Y - YOLD
      DELS=DNRM2(NP1,TZ,1)
C
C USING TWO POINTS AND TANGENTS ON THE HOMOTOPY ZERO CURVE, CONSTRUCT 
C THE HERMITE CUBIC INTERPOLANT Q(S).  THEN USE  ROOT  TO FIND THE S
C CORRESPONDING TO  LAMBDA = 1 .  THE TWO POINTS ON THE ZERO CURVE ARE
C ALWAYS CHOSEN TO BRACKET LAMBDA=1, WITH THE BRACKETING INTERVAL
C ALWAYS BEING [0, DELS].
C
      SA=0.0
      SB=DELS
      LCODE=1
130   CALL ROOT(SOUT,QSOUT,SA,SB,RERR,AERR,LCODE)
      IF (LCODE .GT. 0) GO TO 140
      QSOUT=QOFS(YOLD(NP1),YPOLD(NP1),Y(NP1),YP(NP1),DELS,SOUT) - 1.0
      GO TO 130
C IF LAMBDA = 1 WERE BRACKETED,  ROOT  CANNOT FAIL.
140   IF (LCODE .GT. 2) THEN
        IFLAG=6
        RETURN
      ENDIF
C
C CALCULATE Q(SA) AS THE INITIAL POINT FOR A NEWTON ITERATION.
      DO JW=1,NP1
        W(JW)=QOFS(YOLD(JW),YPOLD(JW),Y(JW),YP(JW),DELS,SA)
      END DO
C
C ***** END OF CALCULATION OF CUBIC INTERPOLANT. *****
C
C TANGENT INFORMATION  YP  IS NO LONGER NEEDED.  HEREAFTER,  YP
C REPRESENTS THE MOST RECENT POINT WHICH IS ON THE OPPOSITE SIDE OF
C THE HYPERPLANE  LAMBDA = 1 FROM  Y.
C
C    PREPARE FOR MAIN LOOP.
C
      YP=YOLD
C
C INITIALIZE  BRACK  TO INDICATE THAT THE POINTS  Y  AND  YOLD  BRACKET
C LAMBDA = 1,  THUS  YOLD = YP .
C
      BRACK = .TRUE.
C
      DO JUDY = 1,LIMIT         ! ***** MAIN LOOP. *****
C CALCULATE NEWTON STEP AT CURRENT ESTIMATE  W .
      SA = -1.0
      CALL TANGNS(SA,W,WP,TZ,YPOLD,A,MODE,LENQR,NFE,N,IFLAG)
      IF (IFLAG .GT. 0) RETURN
C
C NEXT POINT = CURRENT POINT + NEWTON STEP.
C
      W = W + TZ
C
C CHECK FOR CONVERGENCE.
C
      IF ((ABS(W(NP1)-1.0) .LE. RERR+AERR) .AND.
     &    (DNRM2(NP1,TZ,1) .LE. RERR*DNRM2(N,W(1:N),1)+AERR)) THEN
        Y = W
        RETURN
      ENDIF
C
C PREPARE FOR NEXT ITERATION.
C
      IF (ABS(W(NP1)-1.0) .LE. RERR+AERR) THEN
         YPOLD=WP
         CYCLE
      ENDIF
C
C    UPDATE  Y  AND  YOLD .
C
      YOLD=Y
      Y=W
C
C    UPDATE  YP  SUCH THAT  YP  IS THE MOST RECENT POINT
C    OPPOSITE OF  LAMBDA = 1  FROM  Y .  SET  BRACK = .TRUE.  IFF
C    Y  AND  YOLD  BRACKET  LAMBDA = 1  SO THAT  YP = YOLD .
C
          IF ((Y(NP1)-1.0)*(YOLD(NP1)-1.0) .GT. 0) THEN
            BRACK = .FALSE.
          ELSE
            BRACK = .TRUE.
            YP=YOLD
          END IF
C
C    COMPUTE  DELS = ||Y-YP||.
C
          TZ=Y - YP
          DELS=DNRM2(NP1,TZ,1)
C
C       COMPUTE  TZ  FOR THE LINEAR PREDICTOR   W = Y + TZ,
C           WHERE  TZ = SA*(YOLD-Y).
C
          SA = (1.0-Y(NP1))/(YOLD(NP1)-Y(NP1))
          TZ = SA*(YOLD - Y)
C
C       TO INSURE STABILITY, THE LINEAR PREDICTION MUST BE NO FARTHER
C       FROM  Y  THAN  YP  IS.  THIS IS GUARANTEED IF  BRACK = .TRUE.
C       IF LINEAR PREDICTION IS TOO FAR AWAY, USE BRACKETING POINTS
C       TO COMPUTE LINEAR PREDICTION.
C
          IF (.NOT. BRACK) THEN
            IF (DNRM2(NP1,TZ,1) .GT. DELS) THEN
C
C             COMPUTE  TZ = SA*(YP-Y).
C
              SA = (1.0-Y(NP1))/(YP(NP1)-Y(NP1))
              TZ = SA*(YP - Y)
            END IF
          END IF
C
C       COMPUTE ESTIMATE  W = Y + TZ  AND SAVE OLD TANGENT VECTOR.
C
           W = W + TZ
           YPOLD = WP
       END DO          ! ***** END OF MAIN LOOP. *****
C
C THE ALTERNATING SECANT ESTIMATION AND NEWTON ITERATION
C HAS NOT CONVERGED IN  LIMIT  STEPS.  ERROR RETURN.
      IFLAG=6
      RETURN
      END SUBROUTINE ROOTNS
      SUBROUTINE ROOTNX(N,NFE,IFLAG,RELERR,ABSERR,Y,YP,YOLD,
     &   YPOLD,A,GOFW,TZ,W,WP)
C
C  ROOTNX  is an expert user version of ROOTN(F|S), written using the
C reverse call protocol.  All matrix data structures and numerical linear
C algebra are the responsibility of the calling program.  ROOTNX
C indicates to the calling program, via flags, at which points
C RHO(A,LAMBDA,X)  and [ D RHO(A,LAMBDA,X)/D LAMBDA, D RHO(A,LAMBDA,X)/DX ] 
C must be evaluated, and what linear algebra must be done with these
C functions.  ROOTNX  solves the following problem:  given
C YOLD = (LAMBDA(S1),X(S1)), YPOLD = (LAMBDA'(S1),X'(S1)),
C Y = (LAMBDA(S2),X(S2)), YP = (LAMBDA'(S2),X'(S2)), S1 < S2, and a
C continuous function G(Y) = G(LAMBDA,X) such that  G(YOLD) G(Y) <= 0,
C find the point Y(S) = (LAMBDA(S),X(S)), S1 <= S <= S2, such that
C G(Y(S)) = 0.  ROOTNX  alternates between secant estimates of  Y(S)
C and Newton iteration until convergence.
C
C The following interface block should be inserted in the calling
C program:
C
C     INTERFACE
C       SUBROUTINE ROOTNX(N,NFE,IFLAG,RELERR,ABSERR,Y,YP,YOLD,
C    &   YPOLD,A,GOFW,TZ,W,WP)
C       USE HOMOTOPY
C       USE REAL_PRECISION
C       INTEGER, INTENT(IN):: N
C       INTEGER, INTENT(IN OUT):: NFE,IFLAG
C       REAL (KIND=R8), INTENT(IN):: RELERR,ABSERR
C       REAL (KIND=R8), DIMENSION(:), INTENT(IN):: A
C       REAL (KIND=R8), INTENT(IN OUT):: GOFW
C       REAL (KIND=R8), DIMENSION(:), INTENT(IN OUT):: Y,YP,YOLD,YPOLD,
C    &    TZ,W,WP
C       END SUBROUTINE ROOTNX
C     END INTERFACE
C
C
C ON INPUT:
C
C N = dimension of X and the homotopy map.
C
C NFE = number of Jacobian matrix evaluations.
C
C IFLAG = -2, -1, or 0, indicating the problem type, on the first
C         call to  ROOTNX .  ROOTNX  does not distinguish between
C         these values, but they are permitted for consistency with
C         the rest of HOMPACK.
C
C       = 0-10*R, -1-10*R, OR -2-10*R, R = 1,...,5, indicate to  ROOTNX
C         where to resume after a reverse call.  The calling program
C         must not modify  IFLAG  after a reverse call.
C
C RELERR, ABSERR = relative and absolute error values.  The iteration is
C    considered to have converged when a point Y=(LAMBDA,X) is found 
C    such that
C
C    | G(Y(S)) | <= RELERR*||Y|| + ABSERR       and
C
C    ||Z|| <= RELERR*||Y|| + ABSERR  ,          where
C
C    Z is the Newton step to Y=(LAMBDA,X).
C
C Y(1:N+1) = point (LAMBDA(S), X(S)) on zero curve of homotopy map.
C
C YP(1:N+1) = unit tangent vector to the zero curve of the homotopy map
C    at  Y .
C
C YOLD(1:N+1) = a point different from  Y  on the zero curve.
C
C YPOLD(1:N+1) = unit tangent vector to the zero curve of the homotopy
C    map at  YOLD .
C
C A(:) = parameter vector in the homotopy map.
C
C GOFW = G(W), the value requested by some reverse calls.
C
C TZ(1:N+1), W(1:N+1), and WP(1:N+1)  are work arrays used for the
C    Newton step calculation and the interpolation.  On reentry after
C    a reverse call,  WP  and  TZ  contain the tangent vector and
C    Newton step, respectively, at the point  W .  Precisely,
C    D RHO(A,W)/DW WP = 0,  WP^T YPOLD > 0,  ||WP|| = 1,
C    and  TZ  is the minimum norm solution of
C    D RHO(A,W)/DW TZ = - RHO(A,W).
C
C
C ON OUTPUT:
C
C N , RELERR , ABSERR , A  are unchanged.
C
C NFE  has been updated.
C
C IFLAG 
C    = -52, -51, or -50 requests the calling program to
C      return the unit tangent vector in  WP  and the normal flow Newton
C      step in  TZ , all evaluated at the point  W .
C
C    = 0-10*R, -1-10*R, or -2-10*R, R = 1,...,4, requests the calling
C      program to return in  GOFW  the scalar function  G(Y)  evaluated
C      at  W .  The calling program must not modify  IFLAG  after a
C      reverse call.
C
C    = -2, -1, or 0 (unchanged) on a normal return.
C
C    = 4 if a Jacobian matrix with rank < N has occurred.  The
C        iteration was not completed.
C
C    = 6 if the iteration failed to converge.  Y  and  YOLD  contain
C        the last two points found on the zero curve.
C
C    = 7 if input arguments or array sizes are invalid, or  IFLAG  was
C        changed during a reverse call.
C
C Y  is the point on the zero curve of the homotopy map such that
C    G(Y) = 0.
C
C
C Calls  DNRM2 , ROOT .
C
      USE HOMOTOPY
      USE REAL_PRECISION
      INTEGER, INTENT(IN):: N
      INTEGER, INTENT(IN OUT):: NFE,IFLAG
      REAL (KIND=R8), INTENT(IN):: RELERR,ABSERR
      REAL (KIND=R8), DIMENSION(:), INTENT(IN):: A
      REAL (KIND=R8), INTENT(IN OUT):: GOFW
      REAL (KIND=R8), DIMENSION(:), INTENT(IN OUT):: Y,YP,YOLD,YPOLD,
     &    TZ,W,WP
C
C ***** LOCAL VARIABLES. *****
C
      REAL (KIND=R8), SAVE:: AERR,DELS,F0,F1,FP0,FP1,GOFY,
     &  GOFYOLD,GOFYP,RERR,S,SA,SB,SOUT,U
      REAL (KIND=R8):: DD001,DD0011,DD01,DD011,DNRM2,QOFS
      INTEGER, SAVE:: IFLAGC,JUDY,JW,LCODE,LIMIT,NP1
      LOGICAL, SAVE:: BRACK
C
C ***** END OF SPECIFICATION INFORMATION. *****
C
C DEFINITION OF HERMITE CUBIC INTERPOLANT VIA DIVIDED DIFFERENCES.
C
      DD01(F0,F1,DELS)=(F1-F0)/DELS
      DD001(F0,FP0,F1,DELS)=(DD01(F0,F1,DELS)-FP0)/DELS
      DD011(F0,F1,FP1,DELS)=(FP1-DD01(F0,F1,DELS))/DELS
      DD0011(F0,FP0,F1,FP1,DELS)=(DD011(F0,F1,FP1,DELS) - 
     &                            DD001(F0,FP0,F1,DELS))/DELS
      QOFS(F0,FP0,F1,FP1,DELS,S)=((DD0011(F0,FP0,F1,FP1,DELS)*(S-DELS) +
     &   DD001(F0,FP0,F1,DELS))*S + FP0)*S + F0
C
C
      NP1=N+1
      IF (IFLAG > 0) RETURN
      IF ((MOD(-IFLAG,10) > 2) .OR. SIZE(Y) /= NP1 .OR.
     &  SIZE(YP) /= NP1 .OR. SIZE(YOLD) /= NP1 .OR.
     &  SIZE(YPOLD) /= NP1 .OR. SIZE(TZ) /= NP1 .OR.
     &  SIZE(W) /= NP1 .OR. SIZE(WP) /= NP1 .OR.
     &  (IFLAG < -2 .AND. IFLAG /= IFLAGC)) THEN
        IFLAG=7
        RETURN
      ENDIF
      IFLAGC=-MOD(-IFLAG,10)
C
C PICK UP EXECUTION WEHRE IT LEFT OFF AFTER A REVERSE CALL.
C
      IF (IFLAG < -2) THEN
        GO TO (100,110,130,210,200), ABS(IFLAG)/10
      ENDIF
      U=EPSILON(1.0_R8)
      RERR=MAX(RELERR,U)
      AERR=MAX(ABSERR,0.0_R8)
C
C THE LIMIT ON THE NUMBER OF ITERATIONS ALLOWED MAY BE CHANGED BY
C CHANGING THE FOLLOWING STATEMENT:
      LIMIT=2*(INT(ABS(LOG10(AERR+RERR)))+1)
C
      TZ=Y - YOLD
      DELS=DNRM2(NP1,TZ,1)
C EVALUATE  G  AT  YOLD  AND  Y .
      W = YOLD
      IFLAG = IFLAGC - 10
      IFLAGC = IFLAG
      RETURN
 100  GOFYOLD = GOFW
      W = Y
      IFLAG = IFLAGC - 20
      IFLAGC = IFLAG
      RETURN
 110  GOFY = GOFW
C
C USING TWO POINTS AND TANGENTS ON THE HOMOTOPY ZERO CURVE, CONSTRUCT 
C THE HERMITE CUBIC INTERPOLANT Q(S).  THEN USE  ROOT  TO FIND THE S
C CORRESPONDING TO  G(Q(S)) = 0 .  THE TWO POINTS ON THE ZERO CURVE ARE
C ALWAYS CHOSEN TO BRACKET  G(Y(S)) = 0, WITH THE BRACKETING INTERVAL
C ALWAYS BEING [0, DELS].
C
      SA=0.0
      SB=DELS
      LCODE=1
130   DO
        CALL ROOT(SOUT,GOFW,SA,SB,RERR,AERR,LCODE)
        IF (LCODE .GT. 0) EXIT
        DO JW=1,NP1
          W(JW)=QOFS(YOLD(JW),YPOLD(JW),Y(JW),YP(JW),DELS,SOUT)
        END DO
C REQUEST VALUE  G(Q(SOUT))  BY REVERSE CALL.
        IFLAG = IFLAGC - 30
        IFLAGC = IFLAG
        RETURN
      END DO
C IF G(Y(S)) = 0 WERE BRACKETED,  ROOT  CANNOT FAIL.
      IF (LCODE .GT. 2) THEN
        IFLAG=6
        RETURN
      ENDIF
C
C CALCULATE Q(SA) AS THE INITIAL POINT FOR A NEWTON ITERATION.
      DO JW=1,NP1
        W(JW)=QOFS(YOLD(JW),YPOLD(JW),Y(JW),YP(JW),DELS,SA)
      END DO
C
C ***** END OF CALCULATION OF CUBIC INTERPOLANT. *****
C
C TANGENT INFORMATION  YP  IS NO LONGER NEEDED.  HEREAFTER,  YP
C REPRESENTS THE MOST RECENT POINT WHICH IS ON THE OPPOSITE SIDE OF
C THE SOLUTION TO  G(Y(S)) = 0 FROM  Y.
C
C    PREPARE FOR MAIN LOOP.
C
      YP=YOLD
      GOFYP=GOFYOLD
C
C INITIALIZE  BRACK  TO INDICATE THAT THE POINTS  Y  AND  YOLD  BRACKET
C G(Y(S)) = 0,  THUS  YOLD = YP .
C
      BRACK = .TRUE.
C
C ***** MAIN LOOP. *****
      JUDY=1                                  ! DO JUDY = 1,LIMIT
 170  IF (JUDY > LIMIT) GO TO 600
C CALCULATE NEWTON STEP  TZ  AT CURRENT ESTIMATE  W .
      IFLAG = IFLAGC - 50
      IFLAGC = IFLAG
      NFE = NFE+1
      RETURN
C
C NEXT POINT = CURRENT POINT + NEWTON STEP.
C
 200  W = W + TZ
C
C GET FUNCTION VALUE  G(W) .
C
      IFLAG = IFLAGC - 40
      IFLAGC = IFLAG
      RETURN
C
C CHECK FOR CONVERGENCE.
C
 210  SA = RERR*DNRM2(NP1,W,1)+AERR
      IF ((ABS(GOFW) .LE. SA) .AND. (DNRM2(NP1,TZ,1) .LE. SA)) THEN
        Y = W
        IFLAG = IFLAGC
        RETURN
      ENDIF
C
C PREPARE FOR NEXT ITERATION.
C
      IF (ABS(GOFW) .LE. SA) THEN
         YPOLD=WP
         GO TO 590    ! CYCLE
      ENDIF
C
C    UPDATE  Y  AND  YOLD .
C
      YOLD=Y
      Y=W
      GOFYOLD=GOFY
      GOFY=GOFW
C
C    UPDATE  YP  SUCH THAT  YP  IS THE MOST RECENT POINT
C    OPPOSITE OF  G(Y(S)) = 0  FROM  Y .  SET  BRACK = .TRUE.  IFF
C    Y  AND  YOLD  BRACKET  G(Y(S)) = 0  SO THAT  YP = YOLD .
C
      IF ( GOFY * GOFYOLD .GT. 0) THEN
        BRACK = .FALSE.
      ELSE
        BRACK = .TRUE.
        YP=YOLD
        GOFYP=GOFYOLD
      END IF
C
C    COMPUTE  DELS = ||Y-YP||.
C
      TZ=Y - YP
      DELS=DNRM2(NP1,TZ,1)
C
C     COMPUTE  TZ  FOR THE LINEAR PREDICTOR   W = Y + TZ,
C     WHERE  TZ = SA*(YOLD-Y).
C
      S = ABS(GOFY - GOFYOLD)
      IF (S .GE. 1.0) THEN
        SA = GOFY/(GOFY - GOFYOLD)
        TZ = SA*(YOLD - Y)
      ELSE IF (ANY(ABS(GOFY*(YOLD-Y)) .GE. S*HUGE(1.0_R8))) THEN
        TZ = DELS
      ELSE
        SA = GOFY/(GOFY - GOFYOLD)
        TZ = SA*(YOLD - Y)
      END IF
C
C     TO INSURE STABILITY, THE LINEAR PREDICTION MUST BE NO FARTHER
C     FROM  Y  THAN  YP  IS.  THIS IS GUARANTEED IF  BRACK = .TRUE.
C     IF LINEAR PREDICTION IS TOO FAR AWAY, USE BRACKETING POINTS
C     TO COMPUTE LINEAR PREDICTION.
C
      IF (.NOT. BRACK) THEN
        IF (DNRM2(NP1,TZ,1) .GT. DELS) THEN
C
C         COMPUTE  TZ = SA*(YP-Y).
C
          SA = GOFY/(GOFY - GOFYP)
          TZ = SA*(YP - Y)
        END IF
      END IF
C
C     COMPUTE ESTIMATE  W = Y + TZ  AND SAVE OLD TANGENT VECTOR.
C
      W = W + TZ
      YPOLD = WP
 590  JUDY=JUDY+1
      GO TO 170                          ! END DO
C
C ***** END OF MAIN LOOP. *****
C
C THE ALTERNATING SECANT ESTIMATION AND NEWTON ITERATION
C HAS NOT CONVERGED IN  LIMIT  STEPS.  ERROR RETURN.
 600  IFLAG=6
      RETURN
      END SUBROUTINE ROOTNX
      SUBROUTINE ROOTQF(N,NFE,IFLAG,RELERR,ABSERR,Y,YP,YOLD,
     &   YPOLD,A,Q,R,DZ,Z,W,T,F0,F1)
C
C ROOTQF  FINDS THE POINT  YBAR = (1, XBAR)  ON THE ZERO CURVE OF THE
C HOMOTOPY MAP.  IT STARTS WITH TWO POINTS  YOLD=(LAMBDAOLD,XOLD)  AND
C Y=(LAMBDA,X)  SUCH THAT  LAMBDAOLD < 1 <= LAMBDA, AND ALTERNATES
C BETWEEN USING A SECANT METHOD TO FIND A PREDICTED POINT ON THE 
C HYPERPLANE  LAMBDA=1, AND TAKING A QUASI-NEWTON STEP TO RETURN TO THE 
C ZERO CURVE OF THE HOMOTOPY MAP.
C
C THE FOLLOWING INTERFACE BLOCK SHOULD BE INCLUDED IN THE CALLING
C PROGRAM:
C
C     INTERFACE
C       SUBROUTINE ROOTQF(N,NFE,IFLAG,RELERR,ABSERR,Y,YP,YOLD,
C    &    YPOLD,A,Q,R,DZ,Z,W,T,F0,F1)
C       USE REAL_PRECISION
C       REAL (KIND=R8):: RELERR, ABSERR
C       INTEGER:: N, NFE, IFLAG
C       REAL (KIND=R8):: A(:), DZ(N+1), F0(N+1), F1(N+1), 
C    &    Q(N+1,N+1), R((N+1)*(N+2)/2), T(N+1), W(N+1),
C    &    Y(:), YOLD(N+1), YP(N+1), YPOLD(N+1), Z(N+1)
C       END SUBROUTINE ROOTQF
C     END INTERFACE
C
C
C ON INPUT:
C
C N = DIMENSION OF X.
C
C NFE = NUMBER OF JACOBIAN MATRIX EVALUATIONS.
C
C IFLAG = -2, -1, OR 0, INDICATING THE PROBLEM TYPE.
C
C RELERR, ABSERR = RELATIVE AND ABSOLUTE ERROR VALUES.  THE ITERATION IS
C    CONSIDERED TO HAVE CONVERGED WHEN A POINT Y=(LAMBDA,X) IS FOUND 
C    SUCH THAT
C
C    |Y(1) - 1| <= RELERR + ABSERR              AND
C
C    ||DZ|| <= RELERR*||Y|| + ABSERR,           WHERE
C
C    DZ  IS THE QUASI-NEWTON STEP TO Y.
C
C Y(1:N+1) = POINT (LAMBDA(S), X(S)) ON ZERO CURVE OF HOMOTOPY MAP.
C
C YP(1:N+1) = UNIT TANGENT VECTOR TO THE ZERO CURVE OF THE HOMOTOPY MAP
C    AT  Y.
C
C YOLD(1:N+1) = A POINT DIFFERENT FROM  Y  ON THE ZERO CURVE.
C
C YPOLD(1:N+1) = UNIT TANGENT VECTOR TO THE ZERO CURVE OF THE HOMOTOPY
C    MAP AT  YOLD.
C
C A(:) = PARAMETER VECTOR IN THE HOMOTOPY MAP.
C
C Q(1:N+1,1:N+1) CONTAINS  Q  OF THE QR FACTORIZATION OF
C    THE AUGMENTED JACOBIAN MATRIX EVALUATED AT THE POINT Y.
C
C R((N+1)*(N+2)/2) CONTAINS THE UPPER TRIANGLE OF THE R PART OF
C    THE QR FACTORIZATION, STORED BY COLUMNS.
C
C DZ(1:N+1), Z(1:N+1), W(1:N+1), T(1:N+1), F0(1:N+1), F1(1:N+1)
C    ARE WORK ARRAYS USED FOR THE QUASI-NEWTON STEP AND THE SECANT
C    STEP.
C
C
C ON OUTPUT:
C
C N, RELERR, ABSERR, AND A  ARE UNCHANGED.
C
C NFE  HAS BEEN UPDATED.
C
C IFLAG 
C    = -2, -1, OR 0 (UNCHANGED) ON A NORMAL RETURN.
C
C    = 4 IF A SINGULAR JACOBIAN MATRIX OCCURRED.  THE
C        ITERATION WAS NOT COMPLETED.
C
C    = 6 IF THE ITERATION FAILED TO CONVERGE.  Y  AND  YOLD  CONTAIN
C        THE LAST TWO POINTS OBTAINED BY QUASI-NEWTON STEPS, AND  YP
C        CONTAINS A POINT OPPOSITE OF THE HYPERPLANE  LAMBDA=1  FROM
C        Y.
C
C Y  IS THE POINT ON THE ZERO CURVE OF THE HOMOTOPY MAP AT  LAMBDA = 1.
C
C YP  AND  YOLD  CONTAIN POINTS NEAR THE SOLUTION.
C
C CALLS  DGEMV, DNRM2, DTPSV, F (OR RHO), ROOT, UPQRQF.
C
C ***** DECLARATIONS ***** 
      USE HOMOTOPY
      USE REAL_PRECISION
C
C FUNCTION DECLARATIONS 
C
      INTERFACE
        FUNCTION DNRM2(N,X,STRIDE)
        USE REAL_PRECISION
        INTEGER:: N,STRIDE
        REAL (KIND=R8):: DNRM2,X(N)
        END FUNCTION DNRM2
      END INTERFACE
      REAL (KIND=R8):: QOFS
C
C LOCAL VARIABLES 
C
      REAL (KIND=R8):: AERR, DD001, DD0011, DD01, DD011, DELS, ETA, 
     &   ONE, P0, P1, PP0, PP1, QSOUT, RERR, S, SA, SB, SOUT,
     &   U, ZERO
      INTEGER:: ISTEP, I, LCODE, LIMIT,NP1
      LOGICAL:: BRACK
C
C SCALAR ARGUMENTS 
C
      REAL (KIND=R8):: RELERR, ABSERR
      INTEGER:: N, NFE, IFLAG
C
C ARRAY DECLARATIONS 
C
      REAL (KIND=R8):: A(:), DZ(N+1), F0(N+1), F1(N+1), 
     &   Q(N+1,N+1), R((N+1)*(N+2)/2), T(N+1), W(N+1),
     &   Y(:), YOLD(N+1), YP(N+1), YPOLD(N+1), Z(N+1)
C
C ***** END OF DECLARATIONS *****
C
C
C DEFINITION OF HERMITE CUBIC INTERPOLANT VIA DIVIDED DIFFERENCES.
C
      DD01(P0,P1,DELS)=(P1-P0)/DELS
      DD001(P0,PP0,P1,DELS)=(DD01(P0,P1,DELS)-PP0)/DELS
      DD011(P0,P1,PP1,DELS)=(PP1-DD01(P0,P1,DELS))/DELS
      DD0011(P0,PP0,P1,PP1,DELS)=(DD011(P0,P1,PP1,DELS) - 
     &                              DD001(P0,PP0,P1,DELS))/DELS
      QOFS(P0,PP0,P1,PP1,DELS,S)=((DD0011(P0,PP0,P1,PP1,DELS)*
     &     (S-DELS) + DD001(P0,PP0,P1,DELS))*S + PP0)*S + P0
C
C ***** FIRST EXECUTABLE STATEMENT *****
C
C ***** INITIALIZATION *****
C
C ETA = PARAMETER FOR BROYDEN'S UPDATE.
C LIMIT = MAXIMUM NUMBER OF ITERATIONS ALLOWED.
C
      ONE=1.0
      ZERO=0.0
      U=EPSILON(1.0_R8)
      RERR=MAX(RELERR,U)
      AERR=MAX(ABSERR,ZERO)
      NP1=N+1
      ETA = 100.0*U
      LIMIT = 2*(INT(ABS(LOG10(AERR+RERR)))+1)
C
C F0 = (RHO(Y), YP*Y) TRANSPOSE.
C
      IF (IFLAG .EQ. -2) THEN
C
C CURVE TRACKING PROBLEM.
C
        CALL RHO(A,Y(1),Y(2:NP1),F0(1:N))
      ELSE IF (IFLAG .EQ. -1) THEN
C
C ZERO FINDING PROBLEM.
C
        CALL F(Y(2:NP1),F0(1:N))
        F0(1:N) = Y(1)*F0(1:N) + (1.0-Y(1))*(Y(2:NP1)-A(1:N))
      ELSE
C
C FIXED POINT PROBLEM.
C
        CALL F(Y(2:NP1),F0(1:N))
        F0(1:N) = Y(1)*(A(1:N)-F0(1:N))+Y(2:NP1)-A(1:N)
      END IF
      F0(NP1) = DOT_PRODUCT(YP,Y)
C
C ***** END OF INITIALIZATION BLOCK *****
C
C ***** COMPUTE FIRST INTERPOLANT WITH A HERMITE CUBIC *****
C
C FIND DISTANCE BETWEEN Y AND YOLD.  DZ=||Y-YOLD||.
C
      DZ = Y - YOLD
      DELS=DNRM2(NP1,DZ,1)
C
C USING TWO POINTS AND TANGENTS ON THE HOMOTOPY ZERO CURVE, CONSTRUCT 
C THE HERMITE CUBIC INTERPOLANT Q(S).  THEN USE  ROOT  TO FIND THE S
C CORRESPONDING TO  LAMBDA = 1.  THE TWO POINTS ON THE ZERO CURVE ARE
C ALWAYS CHOSEN TO BRACKET  LAMBDA=1, WITH THE BRACKETING INTERVAL
C ALWAYS BEING [0, DELS].
C
      SA=0.0
      SB=DELS
      LCODE=1
 40     CALL ROOT(SOUT,QSOUT,SA,SB,RERR,AERR,LCODE)
        IF (LCODE .GT. 0) GO TO 50
        QSOUT=QOFS(YOLD(1),YPOLD(1),Y(1),YP(1),DELS,SOUT) - 1.0
      GO TO 40
C
C IF  LAMBDA = 1  WERE BRACKETED,  ROOT  CANNOT FAIL.
C
 50   IF (LCODE .GT. 2) THEN
        IFLAG=6
        RETURN
      ENDIF
C
C CALCULATE  Q(SA)  AS THE INITIAL POINT FOR A NEWTON ITERATION.
C
      DO I=1,NP1
        Z(I)=QOFS(YOLD(I),YPOLD(I),Y(I),YP(I),DELS,SA)
      END DO
C
C CALCULATE DZ = Z-Y.
C
      DZ = Z - Y
C
C ***** END OF CALCULATION OF CUBIC INTERPOLANT *****
C
C TANGENT INFORMATION  YPOLD  IS NO LONGER NEEDED.  HEREAFTER,  YPOLD
C REPRESENTS THE MOST RECENT POINT WHICH IS ON THE OPPOSITE SIDE OF
C LAMBDA=1  FROM  Y.
C
C ***** PREPARE FOR MAIN LOOP *****
C
      YPOLD = YOLD
C
C INITIALIZE  BRACK  TO INDICATE THAT THE POINTS  Y  AND  YOLD  BRACKET
C LAMBDA=1,  THUS YOLD = YPOLD.
C
      BRACK = .TRUE.
C
      DO ISTEP=1,LIMIT   ! ***** MAIN LOOP *****
C
C UPDATE JACOBIAN MATRIX.
C  
C       F1=(RHO(Z), YP*Z) TRANSPOSE.
C
        IF (IFLAG .EQ. -2) THEN
          CALL RHO(A,Z(1),Z(2:NP1),F1(1:N))
        ELSE IF (IFLAG .EQ. -1) THEN
          CALL F(Z(2:NP1),F1(1:N))
          F1(1:N) = Z(1)*F1(1:N) + (1.0-Z(1))*(Z(2:NP1)-A(1:N))
        ELSE
          CALL F(Z(2:NP1),F1(1:N))
          F1(1:N) = Z(1)*(A(1:N)-F1(1:N))+Z(2:NP1)-A(1:N)
        END IF
        F1(NP1) = DOT_PRODUCT(YP,Z)
C
C
C PERFORM BROYDEN UPDATE.
C
        CALL UPQRQF(NP1,ETA,DZ,F0,F1,Q,R,W,T)
C
C QUASI-NEWTON STEP.
C
C COMPUTE NEWTON STEP.
C
        T(1:N) = -F1(1:N)
        T(NP1) = 0.0
        CALL DGEMV('T',NP1,NP1,ONE,Q,NP1,T,1,ZERO,DZ,1)
        CALL DTPSV('U', 'N', 'N', NP1, R, DZ, 1)
C
C TAKE NEWTON STEP.
C
        W = Z
        Z = Z + DZ
C
C CHECK FOR CONVERGENCE. 
C
        IF ((ABS(Z(1)-1.0) .LE. RERR+AERR) .AND. 
     &        (DNRM2(NP1,DZ,1) .LE. RERR*DNRM2(N,Z(2),1)+AERR)) THEN
           Y = Z
           RETURN
        END IF
C
C PREPARE FOR NEXT ITERATION.
C
        F0 = F1  
C
C IF  Z(1) = 1.0  THEN PERFORM QUASI-NEWTON ITERATION AGAIN
C WITHOUT COMPUTING A NEW PREDICTOR.
C
        IF (ABS(Z(1)-1.0) .LE. RERR+AERR) THEN
           DZ = Z - W
           CYCLE
        END IF
C
C       UPDATE  Y  AND  YOLD.
C
        YOLD = Y
        Y = Z
C
C UPDATE  YPOLD  SUCH THAT  YPOLD  IS THE MOST RECENT POINT 
C OPPOSITE OF  LAMBDA=1  FROM  Y.  SET  BRACK = .TRUE.  IFF  
C Y & YOLD  BRACKET  LAMBDA=1  SO THAT  YPOLD=YOLD. 
C
        IF ((Y(1)-1.0)*(YOLD(1)-1.0) .GT. 0) THEN
          BRACK = .FALSE.
        ELSE
          BRACK = .TRUE.
          YPOLD = YOLD
        END IF
C
C COMPUTE DELS = ||Y-YPOLD||.
C              
        DZ = Y - YPOLD
        DELS=DNRM2(NP1,DZ,1)
C
C COMPUTE  DZ  FOR THE LINEAR PREDICTOR   Z = Y + DZ,
C           WHERE  DZ = SA*(YOLD-Y).
C
        SA = (1.0-Y(1))/(YOLD(1)-Y(1))
        DZ = SA*(YOLD - Y)
C
C TO INSURE STABILITY, THE LINEAR PREDICTION MUST BE NO FARTHER 
C FROM  Y  THAN  YPOLD  IS.  THIS IS GUARANTEED IF  BRACK = .TRUE.
C IF LINEAR PREDICTION IS TOO FAR AWAY, USE BRACKETING POINTS
C COMPUTE LINEAR PREDICTION.
C
        IF (.NOT. BRACK) THEN
          IF (DNRM2(NP1,DZ,1) .GT. DELS) THEN
C
C COMPUTE  DZ = SA*(YPOLD-Y).
C          
            SA = (1.0-Y(1))/(YPOLD(1)-Y(1))
            DZ = SA*(YPOLD - Y)
          END IF
        END IF
C
C COMPUTE PREDICTOR Z = Y+DZ, AND DZ = NEW Z  - OLD Z (USED FOR
C QUASI-NEWTON UPDATE).
C
        Z = Z + DZ
        DZ = Z - W
      END DO  ! ***** END OF MAIN LOOP. *****
C
C THE ALTERNATING OSCULATORY LINEAR PREDICTION AND QUASI-NEWTON 
C CORRECTION HAS NOT CONVERGED IN  LIMIT  STEPS.  ERROR RETURN.
      IFLAG=6
      RETURN
C
      END SUBROUTINE ROOTQF
      SUBROUTINE SCLGNP(N,MAXT,NUMT,DEG,MODE,EPS0,COEF,
     &  NNUMT,DDEG,CCOEF,ALPHA,BETA,RWORK,XWORK,
     &  FACV,FACE,COESCL,IERR)
C
C SCLGNP  SCALES THE COEFFICIENTS OF A POLYNOMIAL SYSTEM OF N
C EQUATIONS IN N UNKNOWNS, F(X)=0, WHERE THE JTH TERM OF
C THE ITH EQUATION LOOKS LIKE:
C
C    COEF(I,J) * X(1)**DEG(I,1,J) ... X(N)**DEG(I,N,J)
C
C THE ITH EQUATION IS SCALED BY 10**FACE(I).  THE KTH
C VARIABLE IS SCALED BY 10**FACV(K).  IN OTHER WORDS, X(K) =
C 10**FACV(K) * Y(K), WHERE Y SOLVES THE SCALED EQUATION.
C THE SCALED EQUATION HAS THE SAME FORM AS THE ORIGINAL
C EQUATION, EXCEPT THAT COESCL(I,J) REPLACES COEF(I,J), WHERE
C
C COESCL(I,J)=COEF(I,J)* 10**( FACE(I) + FACV(1)*DEG(I,1,J)+ ...
C                                       +FACV(N)*DEG(I,N,J) )
C
C THE CRITERION FOR GENERATING FACE AND FACV IS THAT OF
C MINIMIZING THE SUM OF SQUARES OF THE EXPONENTS OF THE SCALED
C COEFFICIENTS.  IT TURNS OUT THAT THIS CRITERION REDUCES TO
C SOLVING A SINGLE LINEAR SYSTEM, ALPHA*X = BETA, AS DEFINED
C IN THE CODE BELOW.  FURTHER, THE FORM OF THE POLYNOMIAL
C SYSTEM ALONE DETERMINES THE MATRIX ALPHA.  THUS, IN CASES
C IN WHICH MANY SYSTEMS OF THE SAME FORM, BUT WITH DIFFERENT
C COEFFICIENTS, ARE TO BE SCALED, THE MATRIX ALPHA IS
C UNCHANGED AND MAY BE FACTORED ONLY ONCE (BY  DGEQRF).  WHEN
C SCLGNP  IS CALLED WITH MODE=1,  SCLGNP  DOES NOT RECOMPUTE OR
C REFACTOR THE MATRIX ALPHA.  SEE MEINTJES AND MORGAN "A
C METHODOLOGY FOR SOLVING CHEMICAL EQUILIBRIUM SYSTEMS"
C (GENERAL MOTORS RESEARCH LABORATORIES TECHNICAL REPORT
C GMR-4971).
C
C CALLS DIRECTLY: THE LAPACK ROUTINES  DGEQRF,  DORMQR,  THE BLAS 
C ROUTINE  DTRSV.
C
C N  IS THE NUMBER OF EQUATIONS AND THE NUMBER OF VARIABLES.
C
C MAXT  IS THE LEAST UPPER BOUND OF THE SET NUMT(I), I=1 TO N.
C
C NUMT(I)  IS THE NUMBER OF TERMS IN THE I-TH EQUATION FOR I=1 TO N.
C
C DEG(I,K,J)  IS THE DEGREE OF THE K-TH VARIABLE IN THE
C   J-TH TERM OF THE I-TH EQUATION FOR I=1 TO N, J=1 TO NUMT(I), AND
C   K=1 TO N.
C
C MODE  
C  =1  THIS IS NOT THE FIRST CALL TO  SCLGNP, AND THE FORM OF THE
C      SYSTEM HAS NOT CHANGED.
C  =0  THIS IS THE FIRST CALL TO  SCLGNP.
C
C EPS0  ZERO-EPSILON FOR TERMS (TERMS LESS THAN  EPS0  IN MAGNITUDE
C   ARE TREATED AS ZERO BY THE SCALING ALGORITHM).
C
C COEF(I,J)  IS THE COEFFICIENT OF THE JTH TERM OF THE ITH EQUATION
C   FOR I=1 TO N AND J=1 TO NUMT(N).  (COEF(I,J) MAY BE ZERO.)
C
C NNUMT, DDEG, CCOEF, ALPHA, BETA, RWORK, AND  XWORK  ARE WORKSPACES.
C
C ON OUTPUT:
C
C N, NUMT, DEG, MODE, EPS0, AND  COEF  ARE UNCHANGED.
C
C FACV(I)  IS THE VARIABLE SCALE FACTOR FOR THE I-TH VARIABLE, FOR
C   I=1 TO N.
C
C FACE(I)  IS THE EQUATION SCALE FACTOR FOR THE I-TH EQUATION, FOR
C   I=1 TO N.
C
C COESCL(I,J)  IS THE SCALED VERSION OF COEFFICIENT  COEF(I,J), FOR
C   I=1 TO N, J=1 TO NUMT(I), UNLESS IERR=1.
C
C IERR
C   =0  IF SCALING MATRIX, ALPHA, IS WELL CONDITIONED.  
C   =1  OTHERWISE.  IN THIS CASE, ALPHA IS "REPAIRED" AND A
C         SCALING IS COMPUTED.
C
      USE REAL_PRECISION
C
C DECLARATION OF INPUT
      INTEGER, INTENT(IN):: N,MAXT,NUMT(:),DEG(:,:,:),MODE
      REAL (KIND=R8), INTENT(IN):: EPS0,COEF(:,:)
C
C DECLARATION OF WORKSPACE
      INTEGER, INTENT(IN OUT):: NNUMT(N),DDEG(N,N+1,MAXT)
      REAL (KIND=R8), INTENT(IN OUT):: CCOEF(N,MAXT),ALPHA(2*N,2*N),
     &  BETA(2*N),RWORK(N*(2*N+1)),XWORK(2*N)
C
C DECLARATION OF OUTPUT
      REAL (KIND=R8), INTENT(OUT):: FACV(N),FACE(N),COESCL(N,MAXT)
      INTEGER, INTENT(OUT):: IERR
C
C DECLARATION OF LOCAL VARIABLES
      INTEGER:: I,IDAMAX,ICMAX,IRMAX,J,JJ,K,LENR,N2,S
      REAL (KIND=R8):: DUM,LMFPN,NTUR,RTOL,TUM
C
      SAVE
C
      IERR=0
      N2=2*N
      LMFPN=HUGE(1.0_R8)
      NTUR=EPSILON(1.0_R8)*N
      LENR=N*(N+1)/2
C
C  DELETE NEAR ZERO TERMS
      NNUMT = 0
      DO I=1,N
        JJ=0
        DO J=1,NUMT(I)
          IF (ABS(COEF(I,J)) .GT. EPS0) THEN
            JJ=JJ+1
            NNUMT(I)=NNUMT(I)+1
            CCOEF(I,JJ)=COEF(I,J)
            DDEG(I,1:N,JJ)=DEG(I,1:N,J)
          END IF
        END DO
      END DO
      DO I=1,N
        COESCL(I,1:NNUMT(I)) = LOG10(ABS(CCOEF(I,1:NNUMT(I))))
      END DO
C
C SKIP OVER THE GENERATION AND DECOMPOSITON OF MATRIX ALPHA IF MODE=1
      MODE0: IF (MODE .EQ. 0) THEN
C
C GENERATE THE MATRIX ALPHA
      ALPHA(1:N,1:N) = 0.0
      DO S=1,N
        ALPHA(S,S)=NNUMT(S)
      END DO
      DO I=1,N
        ALPHA(N+1:2*N,I) = SUM(DDEG(I,1:N,1:NNUMT(I)),DIM=2)
      END DO
      DO S=1,N
        DO K=1,N
          TUM=0
          DO I=1,N
            DO J=1,NNUMT(I)
              TUM=TUM+DDEG(I,S,J)*DDEG(I,K,J)
            END DO
          END DO
          ALPHA(N+S,N+K)=TUM
        END DO
      END DO
      DO S=1,N
        ALPHA(S,N+1:2*N) = SUM(DDEG(S,1:N,1:NNUMT(S)),DIM=2)
      END DO
C
C COMPUTE QR FACTORIZATION OF MATRIX ALPHA
      CALL DGEQRF(2*N,2*N,ALPHA,2*N,XWORK,BETA,2*N,I)
C
C REPAIR ILL-CONDITIONED SCALING MATRIX
      IRMAX=1        
      ICMAX=1
      DO J=2,N
        I=IDAMAX(J,ALPHA(1,J),1)
        IF (ABS(ALPHA(I,J)) .GT. ABS(ALPHA(IRMAX,ICMAX))) THEN 
          IRMAX=I
          ICMAX=J
        ENDIF
      END DO
      RTOL=ABS(ALPHA(IRMAX,ICMAX))*NTUR
      DO I=N,1,-1
        IF (ABS(ALPHA(I,I)) .LT. RTOL) THEN
          ALPHA(I,I)=LMFPN
          IERR=1
        ENDIF
      END DO
C      
      ENDIF MODE0
C
C CONTROL PASSES HERE IF MODE=1
C
C
C GENERATE THE COLUMN BETA
      DO S=1,N
        BETA(S)=-SUM(COESCL(S,1:NNUMT(S)))
      END DO
      DO S=1,N
        TUM=0
        DO I=1,N
          TUM = TUM + DOT_PRODUCT(COESCL(I,1:NNUMT(I)),
     &      DDEG(I,S,1:NNUMT(I)))
        END DO
        BETA(N+S)=-TUM
      END DO
C
C SOLVE THE LINEAR SYSTEM ALPHA * X = BETA
      CALL DORMQR('L','T',2*N,1,2*N-1,ALPHA,2*N,XWORK,BETA,2*N,RWORK,
     &  N*(2*N+1),I) 
      CALL DTRSV('U','N','N',2*N,ALPHA,2*N,BETA,1) 
C
C GENERATE FACE, FACV, AND THE MATRIX COESCL
      FACE(1:N)=BETA(1:N)
      FACV(1:N)=BETA(N+1:2*N)
      DO I=1,N
        DO J=1,NUMT(I)
          DUM = ABS(COEF(I,J))
          IF (DUM .EQ. 0.0) THEN
            COESCL(I,J) = 0.0
          ELSE
            TUM = FACE(I) + LOG10( DUM ) + DOT_PRODUCT(FACV(1:N),
     &        DEG(I,1:N,J))
            COESCL(I,J) = SIGN(10.0**(TUM), COEF(I,J))
          ENDIF
        END DO
      END DO
      RETURN
      END SUBROUTINE SCLGNP
      SUBROUTINE SINTRP(X,Y,XOUT,YOUT,YPOUT,NEQN,KOLD,PHI,IVC,IV,KGI,GI,
     &                                                ALPHA,G,W,XOLD,P)
C 
C***BEGIN PROLOGUE  SINTRP
C***DATE WRITTEN   740101   (YYMMDD)
C***REVISION DATE  840201   (YYMMDD)
C***CATEGORY NO.  D2A2
C***KEYWORDS  INITIAL VALUE ORDINARY DIFFERENTIAL EQUATIONS,
C             VARIABLE ORDER ADAMS METHODS, SMOOTH INTERPOLANT FOR
C             DEABM IN THE DEPAC PACKAGE
C***AUTHOR  SHAMPINE, L.F.,  SNLA 
C           GORDON, M.K.
C             MODIFIED BY H.A. WATTS
C***PURPOSE  APPROXIMATES THE SOLUTION AT XOUT BY EVALUATING THE
C            POLYNOMIAL COMPUTED IN STEPS AT XOUT.  MUST BE USED IN 
C            CONJUNCTION WITH STEPS.
C***DESCRIPTION 
C 
C   WRITTEN BY L. F. SHAMPINE AND M. K. GORDON
C 
C   ABSTRACT
C 
C 
C   THE METHODS IN SUBROUTINE  STEPS  APPROXIMATE THE SOLUTION NEAR  X
C   BY A POLYNOMIAL.  SUBROUTINE  SINTRP  APPROXIMATES THE SOLUTION AT
C   XOUT  BY EVALUATING THE POLYNOMIAL THERE.  INFORMATION DEFINING THIS
C   POLYNOMIAL IS PASSED FROM  STEPS  SO  SINTRP  CANNOT BE USED ALONE. 
C 
C   THIS CODE IS COMPLETELY EXPLAINED AND DOCUMENTED IN THE TEXT, 
C   COMPUTER SOLUTION OF ORDINARY DIFFERENTIAL EQUATIONS, THE INITIAL 
C   VALUE PROBLEM  BY L. F. SHAMPINE AND M. K. GORDON.
C   FURTHER DETAILS ON USE OF THIS CODE ARE AVAILABLE IN *SOLVING 
C   ORDINARY DIFFERENTIAL EQUATIONS WITH ODE, STEP, AND INTRP*, 
C   BY L. F. SHAMPINE AND M. K. GORDON, SLA-73-1060.
C 
C   INPUT TO SINTRP --
C 
C   THE USER PROVIDES STORAGE IN THE CALLING PROGRAM FOR THE ARRAYS IN
C   THE CALL LIST 
C      DIMENSION Y(NEQN),YOUT(NEQN),YPOUT(NEQN),PHI(NEQN,16),P(NEQN),
C                ALPHA(12),G(13),W(12),GI(11),IV(10)
C   AND DEFINES 
C      XOUT -- POINT AT WHICH SOLUTION IS DESIRED.
C   THE REMAINING PARAMETERS ARE DEFINED IN  STEPS  AND PASSED TO 
C   SINTRP  FROM THAT SUBROUTINE.
C 
C   OUTPUT FROM  SINTRP --
C 
C      YOUT(*) -- SOLUTION AT  XOUT 
C      YPOUT(*) -- DERIVATIVE OF SOLUTION AT  XOUT
C 
C   THE REMAINING PARAMETERS ARE RETURNED UNALTERED FROM THEIR INPUT
C   VALUES.  INTEGRATION WITH  STEPS  MAY BE CONTINUED. 
C 
C***REFERENCES  SHAMPINE L.F., GORDON M.K., *SOLVING ORDINARY 
C                 DIFFERENTIAL EQUATIONS WITH ODE, STEP, AND INTRP*,
C                 SLA-73-1060, SANDIA LABORATORIES, 1973. 
C               WATTS H.A., SHAMPINE L.F., *A SMOOTHER INTERPOLANT FOR
C                 DE/STEP,INTRP : II*, SAND84-0293, SANDIA LABORATORIES,
C                 1984. 
C***ROUTINES CALLED  (NONE) 
C***END PROLOGUE  SINTRP
C 
      USE REAL_PRECISION
      REAL (KIND=R8):: ALP,ALPHA,C,G,GAMMA,GDI,GDIF,GI,GTEMP,
     &  H,HI,HMU,P,PHI,RMU,SIGMA,TEMP1,TEMP2,TEMP3,W,WTEMP,
     &  X,XI,XIM1,XIQ,XOLD,XOUT,Y,YOUT,YPOUT
      INTEGER I,IQ,IV,IVC,IW,J,JQ,KGI,KOLD,KP1,KP2,L,M,NEQN
C
      DIMENSION Y(NEQN),YOUT(NEQN),YPOUT(NEQN),PHI(NEQN,16),P(NEQN)
      DIMENSION GTEMP(13),C(13),WTEMP(13),G(13),W(12),ALPHA(12),
     &          GI(11),IV(10) 
C 
C***FIRST EXECUTABLE STATEMENT
      KP1 = KOLD + 1
      KP2 = KOLD + 2
C 
      HI = XOUT - XOLD
      H = X - XOLD
      XI = HI/H 
      XIM1 = XI - 1.
C 
C   INITIALIZE WTEMP(*) FOR COMPUTING GTEMP(*)
C 
      XIQ = XI
      DO IQ = 1,KP1
        XIQ = XI*XIQ
        TEMP1 = IQ*(IQ+1) 
        WTEMP(IQ) = XIQ/TEMP1 
      END DO
C 
C   COMPUTE THE DOUBLE INTEGRAL TERM GDI
C 
      IF (KOLD .LE. KGI) THEN
        GDI = GI(KOLD)
      ELSE
        IF (IVC .GT. 0) THEN
          IW = IV(IVC)
          GDI = W(IW)
          M = KOLD - IW + 3 
        ELSE
          GDI = 1.0/TEMP1 
          M = 2 
        END IF
        IF (M .LE. KOLD) THEN
          DO I = M,KOLD
            GDI = W(KP2-I) - ALPHA(I)*GDI
          END DO
        END IF
      END IF
C 
C   COMPUTE GTEMP(*) AND C(*) 
C 
      GTEMP(1) = XI 
      GTEMP(2) = 0.5*XI*XI
      C(1) = 1.0
      C(2) = XI 
      IF (KOLD .GE. 2) THEN
        DO I = 2,KOLD
          ALP = ALPHA(I)
          GAMMA = 1.0 + XIM1*ALP
          L = KP2 - I 
          DO JQ = 1,L
            WTEMP(JQ) = GAMMA*WTEMP(JQ) - ALP*WTEMP(JQ+1) 
          END DO
          GTEMP(I+1) = WTEMP(1) 
          C(I+1) = GAMMA*C(I) 
        END DO
      END IF
C 
C   DEFINE INTERPOLATION PARAMETERS 
C 
      SIGMA = (WTEMP(2) - XIM1*WTEMP(1))/GDI
      RMU = XIM1*C(KP1)/GDI 
      HMU = RMU/H 
C 
C   INTERPOLATE FOR THE SOLUTION -- YOUT
C   AND FOR THE DERIVATIVE OF THE SOLUTION -- YPOUT 
C 
      YOUT = 0.0 
      YPOUT = 0.0
      DO J = 1,KOLD 
        I = KP2 - J 
        GDIF = G(I) - G(I-1)
        TEMP2 = (GTEMP(I) - GTEMP(I-1)) - SIGMA*GDIF
        TEMP3 = (C(I) - C(I-1)) + RMU*GDIF
        YOUT = YOUT + TEMP2*PHI(:,I)
        YPOUT = YPOUT + TEMP3*PHI(:,I)
      END DO
      YOUT = ((1.0 - SIGMA)*P + SIGMA*Y) +
     &             H*(YOUT + (GTEMP(1) - SIGMA*G(1))*PHI(:,1))
      YPOUT = HMU*(P - Y) + (YPOUT + (C(1) + RMU*G(1))*PHI(:,1))
C 
      RETURN
      END SUBROUTINE SINTRP
      SUBROUTINE SOLVDS(NN,A,NWK,MAXA,V)
C
C     This subroutine solves a system of linear equations Bx=b, where
C     B is symmetric, and is represented by its LDU factorization.
C
C     Input variables:
C
C        NN  -- dimension of B.
C
C        A -- one dimensional real array containing the upper
C             triangular skyline portion of the LDU decomposition 
C             of the symmetric matrix B.  
C
C        NWK  -- number of elements in A.
C
C        MAXA -- an integer array of length NN+1 which contains the
C                location in A of the diagonal elements of B.  
C                By convention, MAXA(NN+1) = NWK+1 .
C
C        V -- real array of length NN containing the vector b.
C
C 
C     Output variables:
C
C        V -- solution of the system of equations B x = b .
C
C
C     No working storage is required by this routine.
C
      USE REAL_PRECISION
      INTEGER, INTENT(IN):: NN,MAXA(NN+1),NWK
      REAL (KIND=R8), INTENT(IN):: A(NWK)
      REAL (KIND=R8), INTENT(IN OUT):: V(NN)
C local variables.
      INTEGER:: K,KK,KL,KU,L,N
      REAL (KIND=R8):: C
      DO N=1,NN
         KL=MAXA(N)+1
         KU=MAXA(N+1)-1
         IF (KU-KL < 0) CYCLE
         K=N
         C=0.0
         DO KK=KL,KU
            K=K-1
            C=C+A(KK)*V(K)
         END DO
         V(N)=V(N)-C
      END DO
      DO N=1,NN
         K=MAXA(N)
         V(N)=V(N)/A(K)
      END DO
      IF (NN.EQ.1) RETURN
      N=NN
      DO L=2,NN
         KL=MAXA(N) + 1
         KU=MAXA(N+1) - 1
         IF (KU-KL .GE. 0) THEN
           K=N
           DO KK=KL,KU
             K=K - 1
             V(K)=V(K) - A(KK)*V(N)
           END DO
         END IF
         N = N - 1
      END DO
      RETURN
      END SUBROUTINE SOLVDS
      SUBROUTINE STEPDS(F,NEQN,Y,X,H,EPS,WT,START,HOLD,K,KOLD,
     &   CRASH,PHI,P,YP,ALPHA,W,G,KSTEPS,XOLD,IVC,IV,KGI,GI,  
     &   IFLAGC,YPOLD,A,NDIMA,LENQR,MODE,NFEC)
C 
C   STEPDS  IS A MODIFIED FORTRAN 90 VERSION OF  STEPS
C   WRITTEN BY L. F. SHAMPINE AND M. K. GORDON.
C 
C   ABSTRACT
C 
C   SUBROUTINE  STEPS  IS NORMALLY USED INDIRECTLY THROUGH SUBROUTINE 
C   DEABM .  BECAUSE  DEABM  SUFFICES FOR MOST PROBLEMS AND IS MUCH 
C   EASIER TO USE, USING IT SHOULD BE CONSIDERED BEFORE USING  STEPS
C   ALONE.
C 
C   SUBROUTINE STEPS INTEGRATES A SYSTEM OF  NEQN  FIRST ORDER ORDINARY 
C   DIFFERENTIAL EQUATIONS ONE STEP, NORMALLY FROM X TO X+H, USING A
C   MODIFIED DIVIDED DIFFERENCE FORM OF THE ADAMS PECE FORMULAS.  LOCAL 
C   EXTRAPOLATION IS USED TO IMPROVE ABSOLUTE STABILITY AND ACCURACY. 
C   THE CODE ADJUSTS ITS ORDER AND STEP SIZE TO CONTROL THE LOCAL ERROR 
C   PER UNIT STEP IN A GENERALIZED SENSE.  SPECIAL DEVICES ARE INCLUDED 
C   TO CONTROL ROUNDOFF ERROR AND TO DETECT WHEN THE USER IS REQUESTING 
C   TOO MUCH ACCURACY.
C 
C   THIS CODE IS COMPLETELY EXPLAINED AND DOCUMENTED IN THE TEXT, 
C   COMPUTER SOLUTION OF ORDINARY DIFFERENTIAL EQUATIONS, THE INITIAL 
C   VALUE PROBLEM  BY L. F. SHAMPINE AND M. K. GORDON.
C   FURTHER DETAILS ON USE OF THIS CODE ARE AVAILABLE IN *SOLVING 
C   ORDINARY DIFFERENTIAL EQUATIONS WITH ODE, STEP, AND INTRP*, 
C   BY L. F. SHAMPINE AND M. K. GORDON, SLA-73-1060.
C 
C 
C   THE PARAMETERS REPRESENT -- 
C      F -- SUBROUTINE TO EVALUATE DERIVATIVES
C      NEQN -- NUMBER OF EQUATIONS TO BE INTEGRATED 
C      Y(*) -- SOLUTION VECTOR AT X 
C      X -- INDEPENDENT VARIABLE
C      H -- APPROPRIATE STEP SIZE FOR NEXT STEP.  NORMALLY DETERMINED BY
C           CODE
C      EPS -- LOCAL ERROR TOLERANCE 
C      WT(*) -- VECTOR OF WEIGHTS FOR ERROR CRITERION 
C      START -- LOGICAL VARIABLE SET .TRUE. FOR FIRST STEP,  .FALSE.
C           OTHERWISE 
C      HOLD -- STEP SIZE USED FOR LAST SUCCESSFUL STEP
C      K -- APPROPRIATE ORDER FOR NEXT STEP (DETERMINED BY CODE)
C      KOLD -- ORDER USED FOR LAST SUCCESSFUL STEP
C      CRASH -- LOGICAL VARIABLE SET .TRUE. WHEN NO STEP CAN BE TAKEN,
C           .FALSE. OTHERWISE.
C      YP(*) -- DERIVATIVE OF SOLUTION VECTOR AT  X  AFTER SUCCESSFUL 
C           STEP
C      KSTEPS -- COUNTER ON ATTEMPTED STEPS 
C
C   THE VARIABLES X,XOLD,KOLD,KGI AND IVC AND THE ARRAYS Y,PHI,ALPHA,G, 
C   W,P,IV AND GI ARE REQUIRED FOR THE INTERPOLATION SUBROUTINE SINTRP. 
C   THE ARRAYS  YPOLD  AND  A  AND INTEGER CONSTANTS IFLAGC, NDIMA,
C   LENQR, MODE, NFEC ARE WORKING STORAGE PASSED DIRECTLY THROUGH TO
C   FODEDS.
C 
C   INPUT TO STEPS
C 
C      FIRST CALL --
C 
C   THE USER MUST PROVIDE STORAGE IN HIS CALLING PROGRAM FOR ALL ARRAYS 
C   IN THE CALL LIST, NAMELY
C 
C     DIMENSION Y(NEQN),WT(NEQN),PHI(NEQN,16),P(NEQN),YP(NEQN), 
C    &  ALPHA(12),W(12),G(13),GI(11),IV(10), YPOLD(NEQN),A(NDIMA)
C
C                              --                --    **NOTE** 
C 
C   THE USER MUST ALSO DECLARE  START  AND  CRASH 
C   LOGICAL VARIABLES AND  F  AN EXTERNAL SUBROUTINE, SUPPLY THE
C   SUBROUTINE  F(X,Y,YP,NEQN-1,IFLAGC,YPOLD,A,NDIMA,LENQR,MODE,NFEC)
C   TO EVALUATE
C      DY(I)/DX = YP(I) = F(X,Y(1),Y(2),...,Y(NEQN))
C   AND INITIALIZE ONLY THE FOLLOWING PARAMETERS. 
C      NEQN -- NUMBER OF EQUATIONS TO BE INTEGRATED 
C      Y(*) -- VECTOR OF INITIAL VALUES OF DEPENDENT VARIABLES
C      X -- INITIAL VALUE OF THE INDEPENDENT VARIABLE 
C      H -- NOMINAL STEP SIZE INDICATING DIRECTION OF INTEGRATION 
C           AND MAXIMUM SIZE OF STEP.  MUST BE VARIABLE 
C      EPS -- LOCAL ERROR TOLERANCE PER STEP.  MUST BE VARIABLE 
C      WT(*) -- VECTOR OF NON-ZERO WEIGHTS FOR ERROR CRITERION
C      START -- .TRUE.
C      KSTEPS -- SET KSTEPS TO ZERO 
C   DEFINE U TO BE THE MACHINE UNIT ROUNDOFF QUANTITY BY CALLING
C   THE FORTRAN 90 INTRINSIC FUNCTION EPSILON.  U IS THE SMALLEST
C   POSITIVE NUMBER SUCH THAT 1.0+U .GT. 1.0.
C 
C   STEPS  REQUIRES THAT THE L2 NORM OF THE VECTOR WITH COMPONENTS
C   LOCAL ERROR(L)/WT(L)  BE LESS THAN  EPS  FOR A SUCCESSFUL STEP.  THE
C   ARRAY  WT  ALLOWS THE USER TO SPECIFY AN ERROR TEST APPROPRIATE 
C   FOR HIS PROBLEM.  FOR EXAMPLE,
C      WT(L) = 1.0  SPECIFIES ABSOLUTE ERROR, 
C            = ABS(Y(L))  ERROR RELATIVE TO THE MOST RECENT VALUE OF THE
C                 L-TH COMPONENT OF THE SOLUTION, 
C            = ABS(YP(L))  ERROR RELATIVE TO THE MOST RECENT VALUE OF 
C                 THE L-TH COMPONENT OF THE DERIVATIVE, 
C            = MAX(WT(L),ABS(Y(L)))  ERROR RELATIVE TO THE LARGEST
C                 MAGNITUDE OF L-TH COMPONENT OBTAINED SO FAR,
C            = ABS(Y(L))*RELERR/EPS + ABSERR/EPS  SPECIFIES A MIXED 
C                 RELATIVE-ABSOLUTE TEST WHERE  RELERR  IS RELATIVE 
C                 ERROR,  ABSERR  IS ABSOLUTE ERROR AND  EPS =
C                 MAX(RELERR,ABSERR) .
C 
C      SUBSEQUENT CALLS --
C 
C   SUBROUTINE  STEPS  IS DESIGNED SO THAT ALL INFORMATION NEEDED TO
C   CONTINUE THE INTEGRATION, INCLUDING THE STEP SIZE  H  AND THE ORDER 
C   K , IS RETURNED WITH EACH STEP.  WITH THE EXCEPTION OF THE STEP 
C   SIZE, THE ERROR TOLERANCE, AND THE WEIGHTS, NONE OF THE PARAMETERS
C   SHOULD BE ALTERED.  THE ARRAY  WT  MUST BE UPDATED AFTER EACH STEP
C   TO MAINTAIN RELATIVE ERROR TESTS LIKE THOSE ABOVE.  NORMALLY THE
C   INTEGRATION IS CONTINUED JUST BEYOND THE DESIRED ENDPOINT AND THE 
C   SOLUTION INTERPOLATED THERE WITH SUBROUTINE  SINTRP .  IF IT IS 
C   IMPOSSIBLE TO INTEGRATE BEYOND THE ENDPOINT, THE STEP SIZE MAY BE 
C   REDUCED TO HIT THE ENDPOINT SINCE THE CODE WILL NOT TAKE A STEP 
C   LARGER THAN THE  H  INPUT.  CHANGING THE DIRECTION OF INTEGRATION,
C   I.E., THE SIGN OF  H , REQUIRES THE USER SET  START = .TRUE. BEFORE 
C   CALLING  STEPS  AGAIN.  THIS IS THE ONLY SITUATION IN WHICH  START
C   SHOULD BE ALTERED.
C 
C   OUTPUT FROM STEPS 
C 
C      SUCCESSFUL STEP -- 
C 
C   THE SUBROUTINE RETURNS AFTER EACH SUCCESSFUL STEP WITH  START  AND
C   CRASH  SET .FALSE. .  X  REPRESENTS THE INDEPENDENT VARIABLE
C   ADVANCED ONE STEP OF LENGTH  HOLD  FROM ITS VALUE ON INPUT AND  Y 
C   THE SOLUTION VECTOR AT THE NEW VALUE OF  X .  ALL OTHER PARAMETERS
C   REPRESENT INFORMATION CORRESPONDING TO THE NEW  X  NEEDED TO
C   CONTINUE THE INTEGRATION. 
C 
C      UNSUCCESSFUL STEP -- 
C 
C   WHEN THE ERROR TOLERANCE IS TOO SMALL FOR THE MACHINE PRECISION,
C   THE SUBROUTINE RETURNS WITHOUT TAKING A STEP AND  CRASH = .TRUE. .
C   AN APPROPRIATE STEP SIZE AND ERROR TOLERANCE FOR CONTINUING ARE 
C   ESTIMATED AND ALL OTHER INFORMATION IS RESTORED AS UPON INPUT 
C   BEFORE RETURNING.  TO CONTINUE WITH THE LARGER TOLERANCE, THE USER
C   JUST CALLS THE CODE AGAIN.  A RESTART IS NEITHER REQUIRED NOR 
C   DESIRABLE.
C***REFERENCES  SHAMPINE L.F., GORDON M.K., *SOLVING ORDINARY 
C                 DIFFERENTIAL EQUATIONS WITH ODE, STEP, AND INTRP*,
C                 SLA-73-1060, SANDIA LABORATORIES, 1973. 
C 
      USE REAL_PRECISION
      REAL (KIND=R8):: ABSH,EPS,ERK,ERKM1,ERKM2,ERKP1,ERR,FOURU,H,
     &  HNEW,HOLD,P5EPS,R,REALI,REALNS,RHO,ROUND,SUM,TAU,
     &  TEMP1,TEMP2,TEMP3,TEMP4,TEMP5,TEMP6,TWOU,X,XOLD
      INTEGER:: I,IFAIL,IFLAGC,IM1,IP1,IQ,IV(10),IVC,J,JV,K,KGI,
     &  KM1,KM2,KNEW,KOLD,KP1,KP2,KPREV,KSTEPS,L,LENQR,LIMIT1,LIMIT2,
     &  MODE,NDIMA,NEQN,NFEC,NS,NSM2,NSP1,NSP2
      LOGICAL:: CRASH,NORND,PHASE1,START
C
      REAL (KIND=R8):: A(NDIMA),ALPHA(12),BETA(12),G(13),GI(11),
     &  P(NEQN),PHI(NEQN,16),PSI(12),SIG(13),V(12),W(12),
     &  WT(NEQN),Y(NEQN),YP(NEQN),YPOLD(NEQN)
C
C   ALL LOCAL VARIABLES ARE SAVED, RATHER THAN PASSED, IN THIS
C   SPECIALIZED VERSION OF STEPS.
C
      SAVE
C
      INTERFACE
        SUBROUTINE F(S,Y,YP,N,IFLAG,YPOLD,A,NDIMA,LENQR,MODE,NFE)
        USE REAL_PRECISION
        INTEGER:: IFLAG,LENQR,MODE,N,NDIMA,NFE
        REAL (KIND=R8):: A(NDIMA),S,Y(N+1),YP(N+1),YPOLD(N+1)
        END SUBROUTINE F
      END INTERFACE
C 
      REAL (KIND=R8), DIMENSION(13)::
     &  TWO=(/2.0,4.0,8.0,16.0,32.0,64.0,128.0,256.0,512.0,1024.0,
     &  2048.0,4096.0,8192.0/),
     &  GSTR=(/0.500,0.0833,0.0417,0.0264,0.0188,0.0143,0.0114,0.00936,
     &  0.00789,0.00679,0.00592,0.00524,0.00468/)
C 
C 
C       ***     BEGIN BLOCK 0     *** 
C   CHECK IF STEP SIZE OR ERROR TOLERANCE IS TOO SMALL FOR MACHINE
C   PRECISION.  IF FIRST STEP, INITIALIZE PHI ARRAY AND ESTIMATE A
C   STARTING STEP SIZE. 
C                   *** 
C 
C   IF STEP SIZE IS TOO SMALL, DETERMINE AN ACCEPTABLE ONE
C 
C***FIRST EXECUTABLE STATEMENT
      TWOU = 2.0 * EPSILON(1.0_R8)
      FOURU = TWOU + TWOU
      CRASH = .TRUE.
      IF(ABS(H) .GE. FOURU*ABS(X)) GO TO 5
      H = SIGN(FOURU*ABS(X),H)
      RETURN
 5    P5EPS = 0.5*EPS 
C 
C   IF ERROR TOLERANCE IS TOO SMALL, INCREASE IT TO AN ACCEPTABLE VALUE 
C 
      ROUND = 0.0 
      DO L = 1,NEQN
        ROUND = ROUND + (Y(L)/WT(L))**2 
      END DO
      ROUND = TWOU*SQRT(ROUND)
      IF(P5EPS .GE. ROUND) GO TO 15 
      EPS = 2.0*ROUND*(1.0 + FOURU) 
      RETURN
 15   CRASH = .FALSE. 
      G(1) = 1.0
      G(2) = 0.5
      SIG(1) = 1.0
      IF (.NOT.START) GO TO 99 
C 
C   INITIALIZE.  COMPUTE APPROPRIATE STEP SIZE FOR FIRST STEP 
C 
      CALL F(X,Y,YP,NEQN-1,IFLAGC,YPOLD,A,NDIMA,LENQR,MODE,NFEC)
      IF (IFLAGC .GT. 0) RETURN
      SUM = 0.0 
      DO L = 1,NEQN
        PHI(L,1) = YP(L)
        PHI(L,2) = 0.0
        SUM = SUM + (YP(L)/WT(L))**2
      END DO
      SUM = SQRT(SUM) 
      ABSH = ABS(H) 
      IF(EPS .LT. 16.0*SUM*H*H) ABSH = 0.25*SQRT(EPS/SUM) 
      H = SIGN(MAX(ABSH,FOURU*ABS(X)),H)
      HOLD = 0.0
      K = 1 
      KOLD = 0
      KPREV = 0 
      START = .FALSE. 
      PHASE1 = .TRUE. 
      NORND = .TRUE.
      IF(P5EPS .GT. 100.0*ROUND) GO TO 99 
      NORND = .FALSE. 
      PHI(1:NEQN,15) = 0.0 
 99   IFAIL = 0 
C       ***     END BLOCK 0     *** 
C 
C       ***     BEGIN BLOCK 1     *** 
C   COMPUTE COEFFICIENTS OF FORMULAS FOR THIS STEP.  AVOID COMPUTING
C   THOSE QUANTITIES NOT CHANGED WHEN STEP SIZE IS NOT CHANGED. 
C                   *** 
C 
 100  KP1 = K+1 
      KP2 = K+2 
      KM1 = K-1 
      KM2 = K-2 
C 
C   NS IS THE NUMBER OF STEPS TAKEN WITH SIZE H, INCLUDING THE CURRENT
C   ONE.  WHEN K.LT.NS, NO COEFFICIENTS CHANGE
C 
      IF(H .NE. HOLD) NS = 0
      IF (NS.LE.KOLD) NS = NS+1 
      NSP1 = NS+1 
      IF (K .LT. NS) GO TO 199
C 
C   COMPUTE THOSE COMPONENTS OF ALPHA(*),BETA(*),PSI(*),SIG(*) WHICH
C   ARE CHANGED 
C 
      BETA(NS) = 1.0
      REALNS = NS 
      ALPHA(NS) = 1.0/REALNS
      TEMP1 = H*REALNS
      SIG(NSP1) = 1.0 
      IF(K .LT. NSP1) GO TO 110 
      DO I = NSP1,K 
        IM1 = I-1 
        TEMP2 = PSI(IM1)
        PSI(IM1) = TEMP1
        BETA(I) = BETA(IM1)*PSI(IM1)/TEMP2
        TEMP1 = TEMP2 + H 
        ALPHA(I) = H/TEMP1
        REALI = I 
        SIG(I+1) = REALI*ALPHA(I)*SIG(I)
      END DO
 110  PSI(K) = TEMP1
C 
C   COMPUTE COEFFICIENTS G(*) 
C 
C   INITIALIZE V(*) AND SET W(*). 
C 
      IF(NS .GT. 1) GO TO 120 
      DO IQ = 1,K 
        TEMP3 = IQ*(IQ+1) 
        V(IQ) = 1.0/TEMP3 
        W(IQ) = V(IQ) 
      END DO
      IVC = 0 
      KGI = 0 
      IF (K .EQ. 1) GO TO 140 
      KGI = 1 
      GI(1) = W(2)
      GO TO 140 
C 
C   IF ORDER WAS RAISED, UPDATE DIAGONAL PART OF V(*) 
C 
 120  IF (K .LE. KPREV) GO TO 130
      IF (IVC .EQ. 0) GO TO 122 
      JV = KP1 - IV(IVC)
      IVC = IVC - 1 
      GO TO 123 
 122  JV = 1
      TEMP4 = K*KP1 
      V(K) = 1.0/TEMP4
      W(K) = V(K) 
      IF (K .NE. 2) GO TO 123 
      KGI = 1 
      GI(1) = W(2)
 123  NSM2 = NS-2 
      IF (NSM2 .LT. JV) GO TO 130
      DO J = JV,NSM2
        I = K-J 
        V(I) = V(I) - ALPHA(J+1)*V(I+1) 
        W(I) = V(I) 
      END DO
      IF (I .NE. 2) GO TO 130 
      KGI = NS - 1
      GI(KGI) = W(2)
C 
C   UPDATE V(*) AND SET W(*)
C 
 130  LIMIT1 = KP1 - NS 
      TEMP5 = ALPHA(NS) 
      DO IQ = 1,LIMIT1
        V(IQ) = V(IQ) - TEMP5*V(IQ+1) 
        W(IQ) = V(IQ) 
      END DO
      G(NSP1) = W(1)
      IF (LIMIT1 .EQ. 1) GO TO 137
      KGI = NS
      GI(KGI) = W(2)
 137  W(LIMIT1+1) = V(LIMIT1+1) 
      IF (K .GE. KOLD) GO TO 140
      IVC = IVC + 1 
      IV(IVC) = LIMIT1 + 2
C 
C   COMPUTE THE G(*) IN THE WORK VECTOR W(*)
C 
 140  NSP2 = NS + 2 
      KPREV = K 
      IF (KP1 .GE. NSP2) THEN
        DO I = NSP2,KP1 
          LIMIT2 = KP2 - I
          TEMP6 = ALPHA(I-1)
          DO IQ = 1,LIMIT2
            W(IQ) = W(IQ) - TEMP6*W(IQ+1) 
          END DO
          G(I) = W(1) 
        END DO
      END IF
 199  CONTINUE
C       ***     END BLOCK 1     *** 
C 
C       ***     BEGIN BLOCK 2     *** 
C   PREDICT A SOLUTION P(*), EVALUATE DERIVATIVES USING PREDICTED 
C   SOLUTION, ESTIMATE LOCAL ERROR AT ORDER K AND ERRORS AT ORDERS K, 
C   K-1, K-2 AS IF CONSTANT STEP SIZE WERE USED.
C                   *** 
C 
C   INCREMENT COUNTER ON ATTEMPTED STEPS
C 
      KSTEPS = KSTEPS + 1 
C 
C   CHANGE PHI TO PHI STAR
C 
      IF (K .LT. NSP1) GO TO 215 
      DO I = NSP1,K 
        PHI(1:NEQN,I) = BETA(I)*PHI(1:NEQN,I) 
      END DO
C 
C   PREDICT SOLUTION AND DIFFERENCES
C 
 215  PHI(1:NEQN,KP2) = PHI(1:NEQN,KP1) 
      PHI(1:NEQN,KP1) = 0.0
      P(1:NEQN) = 0.0
      DO J = 1,K
        I = KP1 - J 
        IP1 = I+1 
        P(1:NEQN) = P(1:NEQN) + G(I)*PHI(1:NEQN,I)
        PHI(1:NEQN,I) = PHI(1:NEQN,I) + PHI(1:NEQN,IP1)
      END DO
      IF (NORND) THEN
        P(1:NEQN) = Y(1:NEQN) + H*P(1:NEQN)
      ELSE
        DO L = 1,NEQN 
          TAU = H*P(L) - PHI(L,15)
          P(L) = Y(L) + TAU 
          PHI(L,16) = (P(L) - Y(L)) - TAU 
        END DO
      END IF
      XOLD = X
      X = X + H 
      ABSH = ABS(H) 
      CALL F(X,P,YP,NEQN-1,IFLAGC,YPOLD,A,NDIMA,LENQR,MODE,NFEC)
      IF (IFLAGC .GT. 0) RETURN
C 
C   ESTIMATE ERRORS AT ORDERS K,K-1,K-2 
C 
      ERKM2 = 0.0 
      ERKM1 = 0.0 
      ERK = 0.0 
      DO L = 1,NEQN 
        TEMP3 = 1.0/WT(L) 
        TEMP4 = YP(L) - PHI(L,1)
        IF (KM2 > 0) ERKM2 = ERKM2 + ((PHI(L,KM1)+TEMP4)*TEMP3)**2 
        IF (KM2 .GE. 0) ERKM1 = ERKM1 + ((PHI(L,K)+TEMP4)*TEMP3)**2 
        ERK = ERK + (TEMP4*TEMP3)**2
      END DO
      IF (KM2 > 0) ERKM2 = ABSH*SIG(KM1)*GSTR(KM2)*SQRT(ERKM2) 
      IF (KM2 .GE. 0) ERKM1 = ABSH*SIG(K)*GSTR(KM1)*SQRT(ERKM1) 
      TEMP5 = ABSH*SQRT(ERK)
      ERR = TEMP5*(G(K)-G(KP1)) 
      ERK = TEMP5*SIG(KP1)*GSTR(K)
      KNEW = K
C 
C   TEST IF ORDER SHOULD BE LOWERED 
C 
      IF (KM2 > 0) THEN
        IF(MAX(ERKM1,ERKM2) .LE. ERK) KNEW = KM1
      ELSE IF (KM2 .EQ. 0) THEN
        IF(ERKM1 .LE. 0.5*ERK) KNEW = KM1 
      END IF
C 
C   TEST IF STEP SUCCESSFUL 
C 
      IF(ERR .LE. EPS) GO TO 400
C       ***     END BLOCK 2     *** 
C 
C       ***     BEGIN BLOCK 3     *** 
C   THE STEP IS UNSUCCESSFUL.  RESTORE  X, PHI(*,*), PSI(*) . 
C   IF THIRD CONSECUTIVE FAILURE, SET ORDER TO ONE.  IF STEP FAILS MORE 
C   THAN THREE TIMES, CONSIDER AN OPTIMAL STEP SIZE.  DOUBLE ERROR
C   TOLERANCE AND RETURN IF ESTIMATED STEP SIZE IS TOO SMALL FOR MACHINE
C   PRECISION.
C                   *** 
C 
C   RESTORE X, PHI(*,*) AND PSI(*)
C 
      PHASE1 = .FALSE.
      X = XOLD
      DO I = 1,K
        TEMP1 = 1.0/BETA(I) 
        IP1 = I+1 
        PHI(1:NEQN,I) = TEMP1*(PHI(1:NEQN,I) - PHI(1:NEQN,IP1))
      END DO
      IF (K .GE. 2) THEN
        DO I = 2,K
          PSI(I-1) = PSI(I) - H 
        END DO
      END IF
C 
C   ON THIRD FAILURE, SET ORDER TO ONE.  THEREAFTER, USE OPTIMAL STEP 
C   SIZE
C 
      IFAIL = IFAIL + 1 
      TEMP2 = 0.5 
      IF (IFAIL > 3) THEN
        IF (P5EPS .LT. 0.25*ERK) TEMP2 = SQRT(P5EPS/ERK) 
      ENDIF
      IF (IFAIL .GE. 3) KNEW = 1
      H = TEMP2*H 
      K = KNEW
      NS = 0
      IF(ABS(H) .GE. FOURU*ABS(X)) GO TO 340
      CRASH = .TRUE.
      H = SIGN(FOURU*ABS(X),H)
      EPS = EPS + EPS 
      RETURN
 340  GO TO 100 
C       ***     END BLOCK 3     *** 
C 
C       ***     BEGIN BLOCK 4     *** 
C   THE STEP IS SUCCESSFUL.  CORRECT THE PREDICTED SOLUTION, EVALUATE 
C   THE DERIVATIVES USING THE CORRECTED SOLUTION AND UPDATE THE 
C   DIFFERENCES.  DETERMINE BEST ORDER AND STEP SIZE FOR NEXT STEP. 
C                   *** 
 400  KOLD = K
      HOLD = H
C 
C   CORRECT AND EVALUATE
C 
      TEMP1 = H*G(KP1)
      IF (NORND) THEN
        DO L = 1,NEQN 
          TEMP3 = Y(L)
          Y(L) = P(L) + TEMP1*(YP(L) - PHI(L,1))
          P(L) = TEMP3
        END DO
      ELSE
        DO L = 1,NEQN 
          TEMP3 = Y(L)
          RHO = TEMP1*(YP(L) - PHI(L,1)) - PHI(L,16)
          Y(L) = P(L) + RHO 
          PHI(L,15) = (Y(L) - P(L)) - RHO 
          P(L) = TEMP3
        END DO
      END IF
      CALL F(X,Y,YP,NEQN-1,IFLAGC,YPOLD,A,NDIMA,LENQR,MODE,NFEC)
      IF (IFLAGC .GT. 0) RETURN
C 
C   UPDATE DIFFERENCES FOR NEXT STEP
C 
      PHI(1:NEQN,KP1) = YP(1:NEQN) - PHI(1:NEQN,1) 
      PHI(1:NEQN,KP2) = PHI(1:NEQN,KP1) - PHI(1:NEQN,KP2)
      DO I = 1,K
        PHI(1:NEQN,I) = PHI(1:NEQN,I) + PHI(1:NEQN,KP1)
      END DO
C 
C   ESTIMATE ERROR AT ORDER K+1 UNLESS: 
C     IN FIRST PHASE WHEN ALWAYS RAISE ORDER, 
C     ALREADY DECIDED TO LOWER ORDER, 
C     STEP SIZE NOT CONSTANT SO ESTIMATE UNRELIABLE 
C 
      ERKP1 = 0.0 
      IF(KNEW .EQ. KM1  .OR.  K .EQ. 12) PHASE1 = .FALSE. 
      IF(PHASE1) GO TO 450
      IF(KNEW .EQ. KM1) GO TO 455 
      IF(KP1 .GT. NS) GO TO 460 
      DO L = 1,NEQN 
        ERKP1 = ERKP1 + (PHI(L,KP2)/WT(L))**2 
      END DO
      ERKP1 = ABSH*GSTR(KP1)*SQRT(ERKP1)
C 
C   USING ESTIMATED ERROR AT ORDER K+1, DETERMINE APPROPRIATE ORDER 
C   FOR NEXT STEP 
C 
      IF(K .GT. 1) GO TO 445
      IF(ERKP1 .GE. 0.5*ERK) GO TO 460
      GO TO 450 
 445  IF(ERKM1 .LE. MIN(ERK,ERKP1)) GO TO 455 
      IF(ERKP1 .GE. ERK  .OR.  K .EQ. 12) GO TO 460 
C 
C   HERE ERKP1 .LT. ERK .LT. MAX(ERKM1,ERKM2) ELSE ORDER WOULD HAVE 
C   BEEN LOWERED IN BLOCK 2.  THUS ORDER IS TO BE RAISED
C 
C   RAISE ORDER 
C 
 450  K = KP1 
      ERK = ERKP1 
      GO TO 460 
C 
C   LOWER ORDER 
C 
 455  K = KM1 
      ERK = ERKM1 
C 
C   WITH NEW ORDER DETERMINE APPROPRIATE STEP SIZE FOR NEXT STEP
C 
 460  HNEW = H + H
      IF(PHASE1) GO TO 465
      IF(P5EPS .GE. ERK*TWO(K+1)) GO TO 465 
      HNEW = H
      IF(P5EPS .GE. ERK) GO TO 465
      TEMP2 = K+1 
      R = (P5EPS/ERK)**(1.0/TEMP2)
      HNEW = ABSH*MAX(0.5_R8,MIN(0.9_R8,R)) 
      HNEW = SIGN(MAX(HNEW,FOURU*ABS(X)),H) 
 465  H = HNEW
      RETURN
C       ***     END BLOCK 4     *** 
      END SUBROUTINE STEPDS
      SUBROUTINE STEPNF(N,NFE,IFLAG,START,CRASH,HOLD,H,RELERR,
     &   ABSERR,S,Y,YP,YOLD,YPOLD,A,QR,ALPHA,TZ,PIVOT,W,WP,
     &   Z0,Z1,SSPAR)
C
C  STEPNF  TAKES ONE STEP ALONG THE ZERO CURVE OF THE HOMOTOPY MAP
C USING A PREDICTOR-CORRECTOR ALGORITHM.  THE PREDICTOR USES A HERMITE
C CUBIC INTERPOLANT, AND THE CORRECTOR RETURNS TO THE ZERO CURVE ALONG
C THE FLOW NORMAL TO THE DAVIDENKO FLOW.  STEPNF  ALSO ESTIMATES A
C STEP SIZE H FOR THE NEXT STEP ALONG THE ZERO CURVE.  NORMALLY
C  STEPNF  IS USED INDIRECTLY THROUGH  FIXPNF , AND SHOULD BE CALLED
C DIRECTLY ONLY IF IT IS NECESSARY TO MODIFY THE STEPPING ALGORITHM'S
C PARAMETERS.
C
C THE FOLLOWING INTERFACE BLOCK SHOULD BE INCLUDED IN THE CALLING
C PROGRAM:
C
C     INTERFACE
C       SUBROUTINE STEPNF(N,NFE,IFLAG,START,CRASH,HOLD,H,RELERR,
C    &    ABSERR,S,Y,YP,YOLD,YPOLD,A,QR,ALPHA,TZ,PIVOT,W,WP,
C    &    Z0,Z1,SSPAR)
C       USE REAL_PRECISION
C       REAL (KIND=R8):: ABSERR,H,HOLD,RELERR,S
C       INTEGER:: IFLAG,N,NFE
C       LOGICAL:: CRASH,START
C       REAL (KIND=R8):: A(:),ALPHA(3*N+3),QR(N,N+2),SSPAR(8),TZ(N+1),
C    &    W(N+1),WP(N+1),Y(:),YOLD(N+1),YP(N+1),YPOLD(N+1),
C    &    Z0(N+1),Z1(N+1)
C       INTEGER:: PIVOT(N+1)
C       END SUBROUTINE STEPNF
C     END INTERFACE
C
C ON INPUT:
C
C N = DIMENSION OF X AND THE HOMOTOPY MAP.
C
C NFE = NUMBER OF JACOBIAN MATRIX EVALUATIONS.
C
C IFLAG = -2, -1, OR 0, INDICATING THE PROBLEM TYPE.
C
C START = .TRUE. ON FIRST CALL TO  STEPNF , .FALSE. OTHERWISE.
C
C HOLD = ||Y - YOLD||; SHOULD NOT BE MODIFIED BY THE USER.
C
C H = UPPER LIMIT ON LENGTH OF STEP THAT WILL BE ATTEMPTED.  H  MUST BE
C    SET TO A POSITIVE NUMBER ON THE FIRST CALL TO  STEPNF .
C    THEREAFTER  STEPNF  CALCULATES AN OPTIMAL VALUE FOR  H , AND  H
C    SHOULD NOT BE MODIFIED BY THE USER.
C
C RELERR, ABSERR = RELATIVE AND ABSOLUTE ERROR VALUES.  THE ITERATION IS
C    CONSIDERED TO HAVE CONVERGED WHEN A POINT W=(LAMBDA,X) IS FOUND 
C    SUCH THAT
C
C    ||Z|| <= RELERR*||W|| + ABSERR  ,          WHERE
C
C    Z IS THE NEWTON STEP TO W=(LAMBDA,X).
C
C S = (APPROXIMATE) ARC LENGTH ALONG THE HOMOTOPY ZERO CURVE UP TO
C    Y(S) = (LAMBDA(S), X(S)).
C
C Y(1:N+1) = PREVIOUS POINT (LAMBDA(S), X(S)) FOUND ON THE ZERO CURVE OF 
C    THE HOMOTOPY MAP.
C
C YP(1:N+1) = UNIT TANGENT VECTOR TO THE ZERO CURVE OF THE HOMOTOPY MAP
C    AT  Y .
C
C YOLD(1:N+1) = A POINT BEFORE  Y  ON THE ZERO CURVE OF THE HOMOTOPY MAP.
C
C YPOLD(1:N+1) = UNIT TANGENT VECTOR TO THE ZERO CURVE OF THE HOMOTOPY
C    MAP AT  YOLD .
C
C A(:) = PARAMETER VECTOR IN THE HOMOTOPY MAP.
C
C QR(1:N,1:N+2), ALPHA(1:3*N+3), TZ(1:N+1), PIVOT(1:N+1), W(1:N+1), 
C    WP(1:N+1)  ARE WORK ARRAYS USED FOR THE QR FACTORIZATION (IN THE
C    NEWTON STEP CALCULATION) AND THE INTERPOLATION.
C
C Z0(1:N+1), Z1(1:N+1)  ARE WORK ARRAYS USED FOR THE ESTIMATION OF THE
C    NEXT STEP SIZE  H .
C
C SSPAR(1:8) = (LIDEAL, RIDEAL, DIDEAL, HMIN, HMAX, BMIN, BMAX, P)  IS
C    A VECTOR OF PARAMETERS USED FOR THE OPTIMAL STEP SIZE ESTIMATION.
C
C
C ON OUTPUT:
C
C N , A , SSPAR  ARE UNCHANGED.
C
C NFE  HAS BEEN UPDATED.
C
C IFLAG  
C    = -2, -1, OR 0 (UNCHANGED) ON A NORMAL RETURN.
C
C    = 4 IF A JACOBIAN MATRIX WITH RANK < N HAS OCCURRED.  THE
C        ITERATION WAS NOT COMPLETED.
C
C    = 6 IF THE ITERATION FAILED TO CONVERGE.  W  CONTAINS THE LAST
C        NEWTON ITERATE.
C
C START = .FALSE. ON A NORMAL RETURN.
C
C CRASH 
C    = .FALSE. ON A NORMAL RETURN.
C
C    = .TRUE. IF THE STEP SIZE  H  WAS TOO SMALL.  H  HAS BEEN
C      INCREASED TO AN ACCEPTABLE VALUE, WITH WHICH  STEPNF  MAY BE
C      CALLED AGAIN.
C
C    = .TRUE. IF  RELERR  AND/OR  ABSERR  WERE TOO SMALL.  THEY HAVE
C      BEEN INCREASED TO ACCEPTABLE VALUES, WITH WHICH  STEPNF  MAY
C      BE CALLED AGAIN.
C
C HOLD = ||Y - YOLD||.
C
C H = OPTIMAL VALUE FOR NEXT STEP TO BE ATTEMPTED.  NORMALLY  H  SHOULD
C    NOT BE MODIFIED BY THE USER.
C
C RELERR, ABSERR  ARE UNCHANGED ON A NORMAL RETURN.
C
C S = (APPROXIMATE) ARC LENGTH ALONG THE ZERO CURVE OF THE HOMOTOPY MAP 
C    UP TO THE LATEST POINT FOUND, WHICH IS RETURNED IN  Y .
C
C Y, YP, YOLD, YPOLD  CONTAIN THE TWO MOST RECENT POINTS AND TANGENT
C    VECTORS FOUND ON THE ZERO CURVE OF THE HOMOTOPY MAP.
C
C
C CALLS  DNRM2 , TANGNF .
C
      USE REAL_PRECISION
      REAL (KIND=R8):: ABSERR,DCALC,DD001,DD0011,DD01,
     &   DD011,DELS,F0,F1,FOURU,FP0,FP1,H,HFAIL,HOLD,HT,
     &   LCALC,QOFS,RCALC,RELERR,RHOLEN,S,TEMP,TWOU
      INTEGER:: IFLAG,ITNUM,J,JUDY,N,NFE,NP1
      LOGICAL:: CRASH,FAIL,START
C
C ***** ARRAY DECLARATIONS. *****
C
      REAL (KIND=R8):: A(:),ALPHA(3*N+3),QR(N,N+2),SSPAR(8),TZ(N+1),
     &  W(N+1),WP(N+1),Y(:),YOLD(N+1),YP(N+1),YPOLD(N+1),
     &  Z0(N+1),Z1(N+1)
      INTEGER:: PIVOT(N+1)
C
C ***** END OF DIMENSIONAL INFORMATION. *****
C
      INTERFACE
        FUNCTION DNRM2(N,X,STRIDE)
        USE REAL_PRECISION
        INTEGER:: N,STRIDE
        REAL (KIND=R8):: DNRM2,X(N)
        END FUNCTION DNRM2
        SUBROUTINE TANGNF(RHOLEN,Y,YP,YPOLD,A,QR,ALPHA,TZ,PIVOT,
     &    NFE,N,IFLAG)
        USE REAL_PRECISION
        REAL (KIND=R8):: RHOLEN
        INTEGER:: IFLAG,N,NFE
        REAL (KIND=R8):: A(:),Y(:),YP(N+1),YPOLD(N+1)
        REAL (KIND=R8):: ALPHA(3*N+3),QR(N,N+2),TZ(N+1)
        INTEGER:: PIVOT(N+1)
        END SUBROUTINE TANGNF
      END INTERFACE
C
C THE LIMIT ON THE NUMBER OF NEWTON ITERATIONS ALLOWED BEFORE REDUCING
C THE STEP SIZE  H  MAY BE CHANGED BY CHANGING THE FOLLOWING PARAMETER 
C STATEMENT:
      INTEGER, PARAMETER:: LITFH=4
C
C DEFINITION OF HERMITE CUBIC INTERPOLANT VIA DIVIDED DIFFERENCES.
C
      DD01(F0,F1,DELS)=(F1-F0)/DELS
      DD001(F0,FP0,F1,DELS)=(DD01(F0,F1,DELS)-FP0)/DELS
      DD011(F0,F1,FP1,DELS)=(FP1-DD01(F0,F1,DELS))/DELS
      DD0011(F0,FP0,F1,FP1,DELS)=(DD011(F0,F1,FP1,DELS) - 
     &                            DD001(F0,FP0,F1,DELS))/DELS
      QOFS(F0,FP0,F1,FP1,DELS,S)=((DD0011(F0,FP0,F1,FP1,DELS)*(S-DELS) +
     &   DD001(F0,FP0,F1,DELS))*S + FP0)*S + F0
C
C
      TWOU=2.0*EPSILON(1.0_R8)
      FOURU=TWOU+TWOU
      NP1=N+1
      CRASH=.TRUE.
C THE ARCLENGTH  S  MUST BE NONNEGATIVE.
      IF (S .LT. 0.0) RETURN
C IF STEP SIZE IS TOO SMALL, DETERMINE AN ACCEPTABLE ONE.
      IF (H .LT. FOURU*(1.0+S)) THEN
        H=FOURU*(1.0+S)
        RETURN
      ENDIF
C IF ERROR TOLERANCES ARE TOO SMALL, INCREASE THEM TO ACCEPTABLE VALUES.
      TEMP=DNRM2(NP1,Y,1)+1.0
      IF (.5*(RELERR*TEMP+ABSERR) .LT. TWOU*TEMP) THEN
        IF (RELERR .NE. 0.0) THEN
          RELERR=FOURU*(1.0+FOURU)
          ABSERR=MAX(ABSERR,0.0_R8)
        ELSE
          ABSERR=FOURU*TEMP
        ENDIF
        RETURN
      ENDIF
      CRASH=.FALSE.
      STARTUP: IF (START) THEN
C
C *****  STARTUP SECTION (FIRST STEP ALONG ZERO CURVE).  *****
C
      FAIL=.FALSE.
      START=.FALSE.
C DETERMINE SUITABLE INITIAL STEP SIZE.
      H=MIN(H, .10_R8, SQRT(SQRT(RELERR*TEMP+ABSERR)))
C USE LINEAR PREDICTOR ALONG TANGENT DIRECTION TO START NEWTON ITERATION.
      YPOLD(1)=1.0
      YPOLD(2:NP1)=0.0
      CALL TANGNF(S,Y,YP,YPOLD,A,QR,ALPHA,TZ,PIVOT,NFE,N,IFLAG)
      IF (IFLAG .GT. 0) RETURN
      LP: DO
      W=Y + H*YP
      Z0=W
      DO JUDY=1,LITFH
        RHOLEN=-1.0
C CALCULATE THE NEWTON STEP  TZ  AT THE CURRENT POINT  W .
        CALL TANGNF(RHOLEN,W,WP,YPOLD,A,QR,ALPHA,TZ,PIVOT,NFE,N,IFLAG)
        IF (IFLAG .GT. 0) RETURN
C
C TAKE NEWTON STEP AND CHECK CONVERGENCE.
        W=W + TZ
        ITNUM=JUDY
C COMPUTE QUANTITIES USED FOR OPTIMAL STEP SIZE ESTIMATION.
        IF (JUDY .EQ. 1) THEN
          LCALC=DNRM2(NP1,TZ,1)
          RCALC=RHOLEN
          Z1=W
        ELSE IF (JUDY .EQ. 2) THEN
          LCALC=DNRM2(NP1,TZ,1)/LCALC
          RCALC=RHOLEN/RCALC
        ENDIF
C GO TO MOP-UP SECTION AFTER CONVERGENCE.
        IF (DNRM2(NP1,TZ,1) .LE. RELERR*DNRM2(NP1,W,1)+ABSERR)
     &                                                 GO TO 600
C
      END DO
C
C NO CONVERGENCE IN  LITFH  ITERATIONS.  REDUCE  H  AND TRY AGAIN.
      IF (H .LE. FOURU*(1.0 + S)) THEN
        IFLAG=6
        RETURN
      ENDIF
      H=.5 * H
      END DO LP
      END IF STARTUP
C
C ***** END OF STARTUP SECTION. *****
C
C ***** PREDICTOR SECTION. *****
C
      FAIL=.FALSE.
C COMPUTE POINT PREDICTED BY HERMITE INTERPOLANT.  USE STEP SIZE  H
C COMPUTED ON LAST CALL TO  STEPNF .
      HP: DO
      DO J=1,NP1
        W(J)=QOFS(YOLD(J),YPOLD(J),Y(J),YP(J),HOLD,HOLD+H)
      END DO
      Z0=W 
C
C ***** END OF PREDICTOR SECTION. *****
C
C ***** CORRECTOR SECTION. *****
C
      CORRECTOR: DO JUDY=1,LITFH
        RHOLEN=-1.0
C CALCULATE THE NEWTON STEP  TZ  AT THE CURRENT POINT  W .
        CALL TANGNF(RHOLEN,W,WP,YP,A,QR,ALPHA,TZ,PIVOT,NFE,N,IFLAG)
        IF (IFLAG .GT. 0) RETURN
C
C TAKE NEWTON STEP AND CHECK CONVERGENCE.
        W=W + TZ
        ITNUM=JUDY
C COMPUTE QUANTITIES USED FOR OPTIMAL STEP SIZE ESTIMATION.
        IF (JUDY .EQ. 1) THEN
          LCALC=DNRM2(NP1,TZ,1)
          RCALC=RHOLEN
          Z1=W
        ELSE IF (JUDY .EQ. 2) THEN
          LCALC=DNRM2(NP1,TZ,1)/LCALC
          RCALC=RHOLEN/RCALC
        ENDIF
C GO TO MOP-UP SECTION AFTER CONVERGENCE.
        IF (DNRM2(NP1,TZ,1) .LE. RELERR*DNRM2(NP1,W,1)+ABSERR)
     &                                                 GO TO 600
C
      END DO CORRECTOR
C
C NO CONVERGENCE IN  LITFH  ITERATIONS.  RECORD FAILURE AT CALCULATED  H , 
C SAVE THIS STEP SIZE, REDUCE  H  AND TRY AGAIN.
      FAIL=.TRUE.
      HFAIL=H
      IF (H .LE. FOURU*(1.0 + S)) THEN
        IFLAG=6
        RETURN
      ENDIF
      H=.5 * H
      END DO HP
C
C ***** END OF CORRECTOR SECTION. *****
C
C ***** MOP-UP SECTION. *****
C
C YOLD  AND  Y  ALWAYS CONTAIN THE LAST TWO POINTS FOUND ON THE ZERO
C CURVE OF THE HOMOTOPY MAP.  YPOLD  AND  YP  CONTAIN THE TANGENT
C VECTORS TO THE ZERO CURVE AT  YOLD  AND  Y , RESPECTIVELY.
C
600   YPOLD=YP
      YOLD=Y
      Y=W
      YP=WP
      W=Y - YOLD
C UPDATE ARC LENGTH.
      HOLD=DNRM2(NP1,W,1)
      S=S+HOLD
C
C ***** END OF MOP-UP SECTION. *****
C
C ***** OPTIMAL STEP SIZE ESTIMATION SECTION. *****
C
C CALCULATE THE DISTANCE FACTOR  DCALC .
      TZ=Z0 - Y
      W=Z1 - Y
      DCALC=DNRM2(NP1,TZ,1)
      IF (DCALC .NE. 0.0) DCALC=DNRM2(NP1,W,1)/DCALC
C
C THE OPTIMAL STEP SIZE HBAR IS DEFINED BY
C
C   HT=HOLD * [MIN(LIDEAL/LCALC, RIDEAL/RCALC, DIDEAL/DCALC)]**(1/P)
C
C     HBAR = MIN [ MAX(HT, BMIN*HOLD, HMIN), BMAX*HOLD, HMAX ]
C
C IF CONVERGENCE HAD OCCURRED AFTER 1 ITERATION, SET THE CONTRACTION
C FACTOR  LCALC  TO ZERO.
      IF (ITNUM .EQ. 1) LCALC = 0.0
C FORMULA FOR OPTIMAL STEP SIZE.
      IF (LCALC+RCALC+DCALC .EQ. 0.0) THEN
        HT = SSPAR(7) * HOLD
      ELSE 
        HT = (1.0/MAX(LCALC/SSPAR(1), RCALC/SSPAR(2), DCALC/SSPAR(3)))
     &       **(1.0/SSPAR(8)) * HOLD
      ENDIF
C  HT  CONTAINS THE ESTIMATED OPTIMAL STEP SIZE.  NOW PUT IT WITHIN
C REASONABLE BOUNDS.
      H=MIN(MAX(HT,SSPAR(6)*HOLD,SSPAR(4)), SSPAR(7)*HOLD, SSPAR(5))
      IF (ITNUM .EQ. 1) THEN
C IF CONVERGENCE HAD OCCURRED AFTER 1 ITERATION, DON'T DECREASE  H .
        H=MAX(H,HOLD)
      ELSE IF (ITNUM .EQ. LITFH) THEN
C IF CONVERGENCE REQUIRED THE MAXIMUM  LITFH  ITERATIONS, DON'T
C INCREASE  H .
        H=MIN(H,HOLD)
      ENDIF
C IF CONVERGENCE DID NOT OCCUR IN  LITFH  ITERATIONS FOR A PARTICULAR
C H = HFAIL , DON'T CHOOSE THE NEW STEP SIZE LARGER THAN  HFAIL .
      IF (FAIL) H=MIN(H,HFAIL)
C
C
      RETURN
      END SUBROUTINE STEPNF
      SUBROUTINE STEPNS(N,NFE,IFLAG,START,CRASH,HOLD,H,RELERR,
     &   ABSERR,S,Y,YP,YOLD,YPOLD,A,MODE,LENQR,SSPAR,TZ,W,WP,Z0,Z1)
C
C  STEPNS  TAKES ONE STEP ALONG THE ZERO CURVE OF THE HOMOTOPY MAP
C USING A PREDICTOR-CORRECTOR ALGORITHM.  THE PREDICTOR USES A HERMITE
C CUBIC INTERPOLANT, AND THE CORRECTOR RETURNS TO THE ZERO CURVE ALONG
C THE FLOW NORMAL TO THE DAVIDENKO FLOW.  STEPNS  ALSO ESTIMATES A
C STEP SIZE H FOR THE NEXT STEP ALONG THE ZERO CURVE.  NORMALLY
C  STEPNS  IS USED INDIRECTLY THROUGH  FIXPNS , AND SHOULD BE CALLED
C DIRECTLY ONLY IF IT IS NECESSARY TO MODIFY THE STEPPING ALGORITHM'S
C PARAMETERS.  SEE ALSO THE REVERSE CALL ROUTINE  STEPNX .
C
C THE CALLING PROGRAM MUST INCLUDE THE FOLLOWING INTERFACE BLOCK:
C
C     INTERFACE
C       SUBROUTINE STEPNS(N,NFE,IFLAG,START,CRASH,HOLD,H,RELERR,
C    &    ABSERR,S,Y,YP,YOLD,YPOLD,A,MODE,LENQR,SSPAR,TZ,W,WP,Z0,Z1)
C       USE REAL_PRECISION
C       INTEGER, INTENT(IN):: LENQR,MODE,N
C       INTEGER, INTENT(IN OUT):: IFLAG,NFE
C       LOGICAL, INTENT(IN OUT):: CRASH,START
C       REAL (KIND=R8), INTENT(IN):: A(:),SSPAR(8)
C       REAL (KIND=R8), INTENT(IN OUT):: ABSERR,H,HOLD,RELERR,S,
C    &    Y(:),YOLD(:),YP(:),YPOLD(:)
C       REAL (KIND=R8), INTENT(OUT), DIMENSION(:):: TZ,W,WP,Z0,Z1
C       END SUBROUTINE STEPNS
C     END INTERFACE
C
C
C ON INPUT:
C
C N = DIMENSION OF X AND THE HOMOTOPY MAP.
C
C NFE = NUMBER OF JACOBIAN MATRIX EVALUATIONS.
C
C IFLAG = -2, -1, OR 0, INDICATING THE PROBLEM TYPE.
C
C START = .TRUE. ON FIRST CALL TO  STEPNS , .FALSE. OTHERWISE.
C
C HOLD = ||Y - YOLD||; SHOULD NOT BE MODIFIED BY THE USER.
C
C H = UPPER LIMIT ON LENGTH OF STEP THAT WILL BE ATTEMPTED.  H  MUST BE
C    SET TO A POSITIVE NUMBER ON THE FIRST CALL TO  STEPNS .
C    THEREAFTER  STEPNS  CALCULATES AN OPTIMAL VALUE FOR  H , AND  H
C    SHOULD NOT BE MODIFIED BY THE USER.
C
C RELERR, ABSERR = RELATIVE AND ABSOLUTE ERROR VALUES.  THE ITERATION IS
C    CONSIDERED TO HAVE CONVERGED WHEN A POINT W=(X,LAMBDA) IS FOUND 
C    SUCH THAT
C
C    ||Z|| <= RELERR*||W|| + ABSERR  ,          WHERE
C
C    Z IS THE NEWTON STEP TO W=(X,LAMBDA).
C
C S = (APPROXIMATE) ARC LENGTH ALONG THE HOMOTOPY ZERO CURVE UP TO
C    Y(S) = (X(S), LAMBDA(S)).
C
C Y(1:N+1) = PREVIOUS POINT (X(S), LAMBDA(S)) FOUND ON THE ZERO CURVE OF 
C    THE HOMOTOPY MAP.
C
C YP(1:N+1) = UNIT TANGENT VECTOR TO THE ZERO CURVE OF THE HOMOTOPY MAP
C    AT  Y .
C
C YOLD(1:N+1) = A POINT BEFORE  Y  ON THE ZERO CURVE OF THE HOMOTOPY MAP.
C
C YPOLD(1:N+1) = UNIT TANGENT VECTOR TO THE ZERO CURVE OF THE HOMOTOPY
C    MAP AT  YOLD .
C
C A(:) = PARAMETER VECTOR IN THE HOMOTOPY MAP.
C
C MODE = 1 IF THE JACOBIAN MATRIX IS SYMMETRIC AND STORED IN A PACKED
C          SKYLINE FORMAT;
C      = 2 IF THE JACOBIAN MATRIX IS STORED IN A SPARSE ROW FORMAT.
C
C LENQR  IS THE NUMBER OF NONZERO ENTRIES IN THE SPARSE JACOBIAN
C    MATRICES, USED TO DETERMINE THE SPARSE MATRIX DATA STRUCTURES.
C
C SSPAR(1:8) = (LIDEAL, RIDEAL, DIDEAL, HMIN, HMAX, BMIN, BMAX, P)  IS
C    A VECTOR OF PARAMETERS USED FOR THE OPTIMAL STEP SIZE ESTIMATION.
C
C TZ(1:N+1), W(1:N+1), WP(1:N+1), Z0(1:N+1), AND  Z1(1:N+1)  ARE WORK
C    ARRAYS USED FOR THE CALCULATION OF THE JACOBIAN MATRIX KERNEL, THE
C    NEWTON STEP, INTERPOLATION, AND THE ESTIMATION OF THE NEXT STEP
C    SIZE  H .
C
C
C ON OUTPUT:
C
C N , A , SSPAR  ARE UNCHANGED.
C
C NFE  HAS BEEN UPDATED.
C
C IFLAG  
C    = -2, -1, OR 0 (UNCHANGED) ON A NORMAL RETURN.
C
C    = 4 IF THE CONJUGATE GRADIENT ITERATION FAILED TO CONVERGE
C        (MOST LIKELY DUE TO A JACOBIAN MATRIX WITH RANK < N).  THE
C        ITERATION WAS NOT COMPLETED.
C
C    = 6 IF THE NEWTON ITERATION FAILED TO CONVERGE.  W  CONTAINS 
C        THE LAST NEWTON ITERATE.
C
C START = .FALSE. ON A NORMAL RETURN.
C
C CRASH 
C    = .FALSE. ON A NORMAL RETURN.
C
C    = .TRUE. IF THE STEP SIZE  H  WAS TOO SMALL.  H  HAS BEEN
C      INCREASED TO AN ACCEPTABLE VALUE, WITH WHICH  STEPNS  MAY BE
C      CALLED AGAIN.
C
C    = .TRUE. IF  RELERR  AND/OR  ABSERR  WERE TOO SMALL.  THEY HAVE
C      BEEN INCREASED TO ACCEPTABLE VALUES, WITH WHICH  STEPNS  MAY
C      BE CALLED AGAIN.
C
C HOLD = ||Y - YOLD||.
C
C H = OPTIMAL VALUE FOR NEXT STEP TO BE ATTEMPTED.  NORMALLY  H  SHOULD
C    NOT BE MODIFIED BY THE USER.
C
C RELERR, ABSERR  ARE UNCHANGED ON A NORMAL RETURN.
C
C S = (APPROXIMATE) ARC LENGTH ALONG THE ZERO CURVE OF THE HOMOTOPY MAP 
C    UP TO THE LATEST POINT FOUND, WHICH IS RETURNED IN  Y .
C
C Y, YP, YOLD, YPOLD  CONTAIN THE TWO MOST RECENT POINTS AND TANGENT
C    VECTORS FOUND ON THE ZERO CURVE OF THE HOMOTOPY MAP.
C
C
C CALLS  DNRM2 , TANGNS .
C
      USE REAL_PRECISION
      INTEGER, INTENT(IN):: LENQR,MODE,N
      INTEGER, INTENT(IN OUT):: IFLAG,NFE
      LOGICAL, INTENT(IN OUT):: CRASH,START
      REAL (KIND=R8), INTENT(IN):: A(:),SSPAR(8)
      REAL (KIND=R8), INTENT(IN OUT):: ABSERR,H,HOLD,RELERR,S,
     &  Y(:),YOLD(:),YP(:),YPOLD(:)
      REAL (KIND=R8), INTENT(OUT), DIMENSION(:):: TZ,W,WP,Z0,Z1
C
C *****  LOCAL VARIABLES.  *****
C
      REAL (KIND=R8):: DCALC,DD001,DD0011,DD01,DD011,DELS,F0,F1,
     &   FOURU,FP0,FP1,HFAIL,HT,LCALC,QOFS,RCALC,RHOLEN,TEMP,TWOU
      INTEGER:: ITNUM,J,JUDY,NP1
      LOGICAL:: FAIL
C
      INTERFACE
        FUNCTION DNRM2(N,X,STRIDE)
        USE REAL_PRECISION
        INTEGER:: N,STRIDE
        REAL (KIND=R8):: DNRM2,X(N)
        END FUNCTION DNRM2
        SUBROUTINE TANGNS(RHOLEN,Y,YP,TZ,YPOLD,A,MODE,LENQR,
     &    NFE,N,IFLAG)
        USE REAL_PRECISION
        REAL (KIND=R8), INTENT(IN), DIMENSION(:):: A,Y,YPOLD
        REAL (KIND=R8), INTENT(IN OUT):: RHOLEN
        REAL (KIND=R8), INTENT(OUT), DIMENSION(:):: TZ,YP
        INTEGER, INTENT(IN):: LENQR,MODE,N
        INTEGER, INTENT(IN OUT):: IFLAG,NFE
        END SUBROUTINE TANGNS
      END INTERFACE
C
C THE LIMIT ON THE NUMBER OF NEWTON ITERATIONS ALLOWED BEFORE REDUCING
C THE STEP SIZE  H  MAY BE CHANGED BY CHANGING THE FOLLOWING PARAMETER 
C STATEMENT:
      INTEGER, PARAMETER:: LITFH=4
C
C DEFINITION OF HERMITE CUBIC INTERPOLANT VIA DIVIDED DIFFERENCES.
C
      DD01(F0,F1,DELS)=(F1-F0)/DELS
      DD001(F0,FP0,F1,DELS)=(DD01(F0,F1,DELS)-FP0)/DELS
      DD011(F0,F1,FP1,DELS)=(FP1-DD01(F0,F1,DELS))/DELS
      DD0011(F0,FP0,F1,FP1,DELS)=(DD011(F0,F1,FP1,DELS) - 
     &                            DD001(F0,FP0,F1,DELS))/DELS
      QOFS(F0,FP0,F1,FP1,DELS,S)=((DD0011(F0,FP0,F1,FP1,DELS)*(S-DELS) +
     &   DD001(F0,FP0,F1,DELS))*S + FP0)*S + F0
C
C ***** END OF SPECIFICATION INFORMATION. *****
C
C
      TWOU=2.0*EPSILON(1.0_R8)
      FOURU=TWOU+TWOU
      NP1=N+1
      CRASH=.TRUE.
C THE ARCLENGTH  S  MUST BE NONNEGATIVE.
      IF (S .LT. 0.0) RETURN
C IF STEP SIZE IS TOO SMALL, DETERMINE AN ACCEPTABLE ONE.
      IF (H .LT. FOURU*(1.0+S)) THEN
        H=FOURU*(1.0+S)
        RETURN
      ENDIF
C IF ERROR TOLERANCES ARE TOO SMALL, INCREASE THEM TO ACCEPTABLE VALUES.
      TEMP=DNRM2(NP1,Y,1)+1.0
      IF (.5*(RELERR*TEMP+ABSERR) .LT. TWOU*TEMP) THEN
        IF (RELERR .NE. 0.0) THEN
          RELERR=FOURU*(1.0+FOURU)
          ABSERR=MAX(ABSERR,0.0_R8)
        ELSE
          ABSERR=FOURU*TEMP
        ENDIF
        RETURN
      ENDIF
      CRASH=.FALSE.
      STARTUP: IF (START) THEN
C
C *****  STARTUP SECTION (FIRST STEP ALONG ZERO CURVE).  *****
C
      FAIL=.FALSE.
      START=.FALSE.
C DETERMINE SUITABLE INITIAL STEP SIZE.
      H=MIN(H, .10_R8, SQRT(SQRT(RELERR*TEMP+ABSERR)))
C USE LINEAR PREDICTOR ALONG TANGENT DIRECTION TO START NEWTON ITERATION.
      YPOLD(NP1)=1.0
      YPOLD(1:N)=0.0
      CALL TANGNS(S,Y,YP,TZ,YPOLD,A,MODE,LENQR,NFE,N,IFLAG)
      IF (IFLAG .GT. 0) RETURN
      LP: DO
      W=Y + H*YP
      Z0=W
      DO JUDY=1,LITFH
        RHOLEN=-1.0
C CALCULATE THE NEWTON STEP  TZ  AT THE CURRENT POINT  W .
        CALL TANGNS(RHOLEN,W,WP,TZ,YPOLD,A,MODE,LENQR,NFE,N,IFLAG)
        IF (IFLAG .GT. 0) RETURN
C
C TAKE NEWTON STEP AND CHECK CONVERGENCE.
        W=W + TZ
        ITNUM=JUDY
C COMPUTE QUANTITIES USED FOR OPTIMAL STEP SIZE ESTIMATION.
        IF (JUDY .EQ. 1) THEN
          LCALC=DNRM2(NP1,TZ,1)
          RCALC=RHOLEN
          Z1=W
        ELSE IF (JUDY .EQ. 2) THEN
          LCALC=DNRM2(NP1,TZ,1)/LCALC
          RCALC=RHOLEN/RCALC
        ENDIF
C GO TO MOP-UP SECTION AFTER CONVERGENCE.
        IF (DNRM2(NP1,TZ,1) .LE. RELERR*DNRM2(NP1,W,1)+ABSERR)
     &                                                 GO TO 600
C
      END DO
C
C NO CONVERGENCE IN  LITFH  ITERATIONS.  REDUCE  H  AND TRY AGAIN.
      IF (H .LE. FOURU*(1.0 + S)) THEN
        IFLAG=6
        RETURN
      ENDIF
      H=.5 * H
      END DO LP
      END IF STARTUP
C
C ***** END OF STARTUP SECTION. *****
C
C ***** PREDICTOR SECTION. *****
C
      FAIL=.FALSE.
C COMPUTE POINT PREDICTED BY HERMITE INTERPOLANT.  USE STEP SIZE  H
C COMPUTED ON LAST CALL TO  STEPNF .
      HP: DO
      DO J=1,NP1
        W(J)=QOFS(YOLD(J),YPOLD(J),Y(J),YP(J),HOLD,HOLD+H)
      END DO
      Z0=W 
C
C ***** END OF PREDICTOR SECTION. *****
C
C ***** CORRECTOR SECTION. *****
C
      CORRECTOR: DO JUDY=1,LITFH
        RHOLEN=-1.0
C CALCULATE THE NEWTON STEP  TZ  AT THE CURRENT POINT  W .
        CALL TANGNS(RHOLEN,W,WP,TZ,YP,A,MODE,LENQR,NFE,N,IFLAG)
        IF (IFLAG .GT. 0) RETURN
C
C TAKE NEWTON STEP AND CHECK CONVERGENCE.
        W=W + TZ
        ITNUM=JUDY
C COMPUTE QUANTITIES USED FOR OPTIMAL STEP SIZE ESTIMATION.
        IF (JUDY .EQ. 1) THEN
          LCALC=DNRM2(NP1,TZ,1)
          RCALC=RHOLEN
          Z1=W
        ELSE IF (JUDY .EQ. 2) THEN
          LCALC=DNRM2(NP1,TZ,1)/LCALC
          RCALC=RHOLEN/RCALC
        ENDIF
C GO TO MOP-UP SECTION AFTER CONVERGENCE.
        IF (DNRM2(NP1,TZ,1) .LE. RELERR*DNRM2(NP1,W,1)+ABSERR)
     &                                                 GO TO 600
C
      END DO CORRECTOR
C
C NO CONVERGENCE IN  LITFH  ITERATIONS.  RECORD FAILURE AT CALCULATED  H , 
C SAVE THIS STEP SIZE, REDUCE  H  AND TRY AGAIN.
      FAIL=.TRUE.
      HFAIL=H
      IF (H .LE. FOURU*(1.0 + S)) THEN
        IFLAG=6
        RETURN
      ENDIF
      H=.5 * H
      END DO HP
C
C ***** END OF CORRECTOR SECTION. *****
C
C ***** MOP-UP SECTION. *****
C
C YOLD  AND  Y  ALWAYS CONTAIN THE LAST TWO POINTS FOUND ON THE ZERO
C CURVE OF THE HOMOTOPY MAP.  YPOLD  AND  YP  CONTAIN THE TANGENT
C VECTORS TO THE ZERO CURVE AT  YOLD  AND  Y , RESPECTIVELY.
C
600   YPOLD=YP
      YOLD=Y
      Y=W
      YP=WP
      W=Y - YOLD
C UPDATE ARC LENGTH.
      HOLD=DNRM2(NP1,W,1)
      S=S+HOLD
C
C ***** END OF MOP-UP SECTION. *****
C
C ***** OPTIMAL STEP SIZE ESTIMATION SECTION. *****
C
C CALCULATE THE DISTANCE FACTOR  DCALC .
      TZ=Z0 - Y
      W=Z1 - Y
      DCALC=DNRM2(NP1,TZ,1)
      IF (DCALC .NE. 0.0) DCALC=DNRM2(NP1,W,1)/DCALC
C
C THE OPTIMAL STEP SIZE HBAR IS DEFINED BY
C
C   HT=HOLD * [MIN(LIDEAL/LCALC, RIDEAL/RCALC, DIDEAL/DCALC)]**(1/P)
C
C     HBAR = MIN [ MAX(HT, BMIN*HOLD, HMIN), BMAX*HOLD, HMAX ]
C
C IF CONVERGENCE HAD OCCURRED AFTER 1 ITERATION, SET THE CONTRACTION
C FACTOR  LCALC  TO ZERO.
      IF (ITNUM .EQ. 1) LCALC = 0.0
C FORMULA FOR OPTIMAL STEP SIZE.
      IF (LCALC+RCALC+DCALC .EQ. 0.0) THEN
        HT = SSPAR(7) * HOLD
      ELSE 
        HT = (1.0/MAX(LCALC/SSPAR(1), RCALC/SSPAR(2), DCALC/SSPAR(3)))
     &       **(1.0/SSPAR(8)) * HOLD
      ENDIF
C  HT  CONTAINS THE ESTIMATED OPTIMAL STEP SIZE.  NOW PUT IT WITHIN
C REASONABLE BOUNDS.
      H=MIN(MAX(HT,SSPAR(6)*HOLD,SSPAR(4)), SSPAR(7)*HOLD, SSPAR(5))
      IF (ITNUM .EQ. 1) THEN
C IF CONVERGENCE HAD OCCURRED AFTER 1 ITERATION, DON'T DECREASE  H .
        H=MAX(H,HOLD)
      ELSE IF (ITNUM .EQ. LITFH) THEN
C IF CONVERGENCE REQUIRED THE MAXIMUM  LITFH  ITERATIONS, DON'T
C INCREASE  H .
        H=MIN(H,HOLD)
      ENDIF
C IF CONVERGENCE DID NOT OCCUR IN  LITFH  ITERATIONS FOR A PARTICULAR
C H = HFAIL , DON'T CHOOSE THE NEW STEP SIZE LARGER THAN  HFAIL .
      IF (FAIL) H=MIN(H,HFAIL)
C
C
      RETURN
      END SUBROUTINE STEPNS
      SUBROUTINE STEPNX(N,NFE,IFLAG,START,CRASH,HOLD,H,RELERR,
     &   ABSERR,S,Y,YP,YOLD,YPOLD,A,TZ,W,WP,RHOLEN,SSPAR)
C
C  STEPNX  takes one step along the zero curve of the homotopy map
C using a predictor-corrector algorithm.  The predictor uses a Hermite
C cubic interpolant, and the corrector returns to the zero curve along
C the flow normal to the Davidenko flow.  STEPNX  also estimates a
C step size H for the next step along the zero curve.  STEPNX  is an
C expert user version of STEPN(F|S), written using the reverse call
C protocol.  All matrix data structures and numerical linear algebra
C are the responsibility of the calling program.  STEPNX  indicates to
C the calling program, via flags, at which points  RHO(A,LAMBDA,X)  and
C [ D RHO(A,LAMBDA,X)/D LAMBDA, D RHO(A,LAMBDA,X)/DX ]  must be
C evaluated, and what linear algebra must be done with these functions.
C Out of range arguments can also be signaled to  STEPNX , which will
C attempt to modify its steplength algorithm to reflect this
C information.
C
C The following interface block should be inserted in the calling
C program:
C
C     INTERFACE
C       SUBROUTINE STEPNX(N,NFE,IFLAG,START,CRASH,HOLD,H,RELERR,
C    &    ABSERR,S,Y,YP,YOLD,YPOLD,A,TZ,W,WP,RHOLEN,SSPAR)
C       USE HOMOTOPY
C       USE REAL_PRECISION
C       INTEGER, INTENT(IN):: N
C       INTEGER, INTENT(IN OUT):: NFE,IFLAG
C       LOGICAL, INTENT(IN OUT):: START,CRASH
C       REAL (KIND=R8), INTENT(IN OUT):: HOLD,H,RELERR,ABSERR,S,RHOLEN,
C    &    SSPAR(8)
C       REAL (KIND=R8), DIMENSION(:), INTENT(IN):: A
C       REAL (KIND=R8), DIMENSION(:), INTENT(IN OUT):: Y,YP,YOLD,YPOLD,
C    &    TZ,W,WP
C       REAL (KIND=R8), DIMENSION(:), ALLOCATABLE, SAVE:: Z0,Z1
C       END SUBROUTINE STEPNX
C     END INTERFACE
C
C ON INPUT:
C
C N = dimension of X and the homotopy map.
C
C NFE = number of Jacobian matrix evaluations.
C
C IFLAG = -2, -1, or 0, indicating the problem type, on the first
C         call to  STEPNX .  STEPNX  does not distinguish between
C         these values, but they are permitted for consistency with
C         the rest of HOMPACK.
C
C       = 0-10*R, -1-10*R, or -2-10*R, R = 1,2,3, indicate to  STEPNX
C         where to resume after a reverse call.  The calling program
C         must not modify  IFLAG  after a reverse call, except as
C         noted next.
C
C       = -40, -41, or -42, used for a final call to deallocate working
C         storage, after all path tracking is finished.  START  and
C         IFLAG  are reset on return.
C
C       = -100-10*R, -101-10*R, -102-10*R, R = 1,2,3, indicate to
C         STEPNX  where to resume after a reverse call, and that the
C         requested evaluation point was out of range.  STEPNX  will
C         reduce  H  and try again.
C
C START = .TRUE. on first call to  STEPNX , .FALSE. otherwise.
C
C HOLD = ||Y - YOLD||; should not be modified by the user.
C
C H = upper limit on length of step that will be attempted.  H  must be
C    set to a positive number on the first call to  STEPNX .
C    Thereafter  STEPNX  calculates an optimal value for  H , and  H
C    should not be modified by the user.
C
C RELERR, ABSERR = relative and absolute error values.  The iteration is
C    considered to have converged when a point W=(LAMBDA,X) is found 
C    such that
C
C    ||Z|| <= RELERR*||W|| + ABSERR  ,          where
C
C    Z is the Newton step to W=(LAMBDA,X).
C
C S = (approximate) arc length along the homotopy zero curve up to
C    Y(S) = (LAMBDA(S), X(S)).
C
C Y(1:N+1) = previous point (LAMBDA(S), X(S)) found on the zero curve of 
C    the homotopy map.
C
C YP(1:N+1) = unit tangent vector to the zero curve of the homotopy map
C    at  Y .
C
C YOLD(1:N+1) = a point before  Y  on the zero curve of the homotopy map.
C
C YPOLD(1:N+1) = unit tangent vector to the zero curve of the homotopy
C    map at  YOLD .
C
C A(:) = parameter vector in the homotopy map.
C
C TZ(1:N+1), W(1:N+1), and WP(1:N+1)  are work arrays used for the
C    Newton step calculation and the interpolation.  On reentry after
C    a reverse call,  WP  and  TZ  contain the tangent vector and
C    Newton step, respectively, at the point  W .  Precisely,
C    D RHO(A,W)/DW WP = 0,  WP^T YP > 0,  ||WP|| = 1,
C    and  TZ  is the minimum norm solution of
C    D RHO(A,W)/DW TZ = - RHO(A,W).
C
C RHOLEN = ||RHO(A,W)||_2 is required by some reverse calls.
C
C SSPAR(1:8) = (LIDEAL, RIDEAL, DIDEAL, HMIN, HMAX, BMIN, BMAX, P)  is
C    a vector of parameters used for the optimal step size estimation.
C    If  SSPAR(J) .LE. 0.0  on input, it is reset to a default value
C    by  STEPNX .  Otherwise the input value of  SSPAR(J)  is used.
C    See the comments below in  STEPNX  for more information about
C    these constants.
C
C
C ON OUTPUT:
C
C N  and  A  are unchanged.
C
C NFE  has been updated.
C
C IFLAG  
C    = -22, -21, -20, -32, -31, or -30 requests the calling program to
C      return the unit tangent vector in  WP , the normal flow Newton
C      step in  TZ , and the 2-norm of the homotopy map in  RHOLEN ,
C      all evaluated at the point  W .
C
C    = -12, -11, or -10 requests the calling program to return in  WP
C      the unit tangent vector at  W .
C
C    = -2, -1, or 0 (unchanged) on a normal return after a successful
C      step.
C
C    = 4 if a Jacobian matrix with rank < N has occurred.  The
C        iteration was not completed.
C
C    = 6 if the iteration failed to converge.  W  contains the last
C        Newton iterate.
C
C    = 7 if input arguments or array sizes are invalid, or  IFLAG  was
C        changed during a reverse call.
C
C START = .FALSE. on a normal return.
C
C CRASH 
C    = .FALSE. on a normal return.
C
C    = .TRUE. if the step size  H  was too small.  H  has been
C      increased to an acceptable value, with which  STEPNX  may be
C      called again.
C
C    = .TRUE. if  RELERR  and/or  ABSERR  were too small.  They have
C      been increased to acceptable values, with which  STEPNX  may
C      be called again.
C
C HOLD = ||Y - YOLD||.
C
C H = optimal value for next step to be attempted.  Normally  H  should
C    not be modified by the user.
C
C RELERR, ABSERR  are unchanged on a normal return.
C
C S = (approximate) arc length along the zero curve of the homotopy map 
C    up to the latest point found, which is returned in  Y .
C
C Y, YP, YOLD, YPOLD  contain the two most recent points and tangent
C    vectors found on the zero curve of the homotopy map.
C
C SSPAR  may have been changed to default values.
C
C
C Z0(1:N+1), Z1(1:N+1)  are allocatable work arrays used for the
C    estimation of the next step size  H .
C
C Calls  DNRM2 .
C
      USE HOMOTOPY
      USE REAL_PRECISION
      INTEGER, INTENT(IN):: N
      INTEGER, INTENT(IN OUT):: NFE,IFLAG
      LOGICAL, INTENT(IN OUT):: START,CRASH
      REAL (KIND=R8), INTENT(IN OUT):: HOLD,H,RELERR,ABSERR,S,RHOLEN,
     &  SSPAR(8)
      REAL (KIND=R8), DIMENSION(:), INTENT(IN):: A
      REAL (KIND=R8), DIMENSION(:), INTENT(IN OUT):: Y,YP,YOLD,YPOLD,
     &  TZ,W,WP
      REAL (KIND=R8), DIMENSION(:), ALLOCATABLE, SAVE:: Z0,Z1
C
C ***** LOCAL VARIABLES. *****
C
      REAL (KIND=R8), SAVE:: DCALC,DELS,F0,F1,FOURU,FP0,FP1,
     &  HFAIL,HT,LCALC,RCALC,TEMP,TWOU
      INTEGER, SAVE:: IFLAGC,ITNUM,J,JUDY,NP1
      LOGICAL, SAVE:: FAIL
C
C ***** END OF SPECIFICATION INFORMATION. *****
C
C THE LIMIT ON THE NUMBER OF NEWTON ITERATIONS ALLOWED BEFORE REDUCING
C THE STEP SIZE  H  MAY BE CHANGED BY CHANGING THE FOLLOWING PARAMETER 
C STATEMENT:
      INTEGER, PARAMETER:: LITFH=4
C
C DEFINITION OF HERMITE CUBIC INTERPOLANT VIA DIVIDED DIFFERENCES.
C
      REAL (KIND=R8):: DD001,DD0011,DD01,DD011,DNRM2,QOFS
      DD01(F0,F1,DELS)=(F1-F0)/DELS
      DD001(F0,FP0,F1,DELS)=(DD01(F0,F1,DELS)-FP0)/DELS
      DD011(F0,F1,FP1,DELS)=(FP1-DD01(F0,F1,DELS))/DELS
      DD0011(F0,FP0,F1,FP1,DELS)=(DD011(F0,F1,FP1,DELS) - 
     &                            DD001(F0,FP0,F1,DELS))/DELS
      QOFS(F0,FP0,F1,FP1,DELS,S)=((DD0011(F0,FP0,F1,FP1,DELS)*(S-DELS) +
     &   DD001(F0,FP0,F1,DELS))*S + FP0)*S + F0
C
C
      NP1=N+1
      IF (IFLAG > 0) RETURN
      IF ((START .AND. IFLAG < -2) .OR. SIZE(Y) /= NP1 .OR.
     &  SIZE(YP) /= NP1 .OR. SIZE(YOLD) /= NP1 .OR.
     &  SIZE(YPOLD) /= NP1 .OR. SIZE(TZ) /= NP1 .OR.
     &  SIZE(W) /= NP1 .OR. SIZE(WP) /= NP1 .OR.
     &  (.NOT. START .AND. -MOD(-IFLAG,100) /= IFLAGC .AND.
     &  ABS(IFLAG)/10 /= 4)) THEN
        IFLAG=7
        RETURN
      ENDIF
      IFLAGC=-MOD(-IFLAG,10)
C
C PICK UP EXECUTION WEHRE IT LEFT OFF AFTER A REVERSE CALL.
C
      IF (IFLAG < -2) THEN
        GO TO (50,100,400,700), MOD(ABS(IFLAG),100)/10
      ENDIF
      TWOU=2.0*EPSILON(1.0_R8)
      FOURU=TWOU+TWOU
      CRASH=.TRUE.
C THE ARCLENGTH  S  MUST BE NONNEGATIVE.
      IF (S .LT. 0.0) RETURN
C IF STEP SIZE IS TOO SMALL, DETERMINE AN ACCEPTABLE ONE.
      IF (H .LT. FOURU*(1.0+S)) THEN
        H=FOURU*(1.0+S)
        RETURN
      ENDIF
C IF ERROR TOLERANCES ARE TOO SMALL, INCREASE THEM TO ACCEPTABLE VALUES.
      TEMP=DNRM2(NP1,Y,1)+1.0
      IF (.5*(RELERR*TEMP+ABSERR) .LT. TWOU*TEMP) THEN
        IF (RELERR .NE. 0.0) THEN
          RELERR=FOURU*(1.0+FOURU)
          ABSERR=MAX(ABSERR,0.0_R8)
        ELSE
          ABSERR=FOURU*TEMP
        ENDIF
        RETURN
      ENDIF
      CRASH=.FALSE.
      IF (.NOT. START) GO TO 300
C
C *****  STARTUP SECTION (FIRST STEP ALONG ZERO CURVE).  *****
C
      FAIL=.FALSE.
      START=.FALSE.
      IF (ALLOCATED(Z0)) DEALLOCATE(Z0)
      IF (ALLOCATED(Z1)) DEALLOCATE(Z1)
      ALLOCATE(Z0(NP1),Z1(NP1))
C
C SET OPTIMAL STEP SIZE ESTIMATION PARAMETERS.
C LET Z[K] DENOTE THE NEWTON ITERATES ALONG THE FLOW NORMAL TO THE
C DAVIDENKO FLOW AND Y THEIR LIMIT.
C IDEAL CONTRACTION FACTOR:  ||Z[2] - Z[1]|| / ||Z[1] - Z[0]||
      IF (SSPAR(1) .LE. 0.0) SSPAR(1)= .5
C IDEAL RESIDUAL FACTOR:  ||RHO(A, Z[1])|| / ||RHO(A, Z[0])||
      IF (SSPAR(2) .LE. 0.0) SSPAR(2)= .01
C IDEAL DISTANCE FACTOR:  ||Z[1] - Y|| / ||Z[0] - Y||
      IF (SSPAR(3) .LE. 0.0) SSPAR(3)= .5
C MINIMUM STEP SIZE  HMIN .
      IF (SSPAR(4) .LE. 0.0) SSPAR(4)=(SQRT(N+1.0)+4.0)*EPSILON(1.0_R8)
C MAXIMUM STEP SIZE  HMAX .
      IF (SSPAR(5) .LE. 0.0) SSPAR(5)= 1.0
C MINIMUM STEP SIZE REDUCTION FACTOR  BMIN .
      IF (SSPAR(6) .LE. 0.0) SSPAR(6)= .1_R8
C MAXIMUM STEP SIZE EXPANSION FACTOR  BMAX .
      IF (SSPAR(7) .LE. 0.0) SSPAR(7)= 3.0
C ASSUMED OPERATING ORDER  P .
      IF (SSPAR(8) .LE. 0.0) SSPAR(8)= 2.0
C
C DETERMINE SUITABLE INITIAL STEP SIZE.
      H=MIN(H, .10_R8, SQRT(SQRT(RELERR*TEMP+ABSERR)))
C USE LINEAR PREDICTOR ALONG TANGENT DIRECTION TO START NEWTON ITERATION.
      YPOLD(1)=1.0
      YPOLD(2:NP1)=0.0
C REQUEST TANGENT VECTOR AT Y VIA REVERSE CALL.
      W=Y
      YP=YPOLD
      IFLAG=IFLAGC-10
      IFLAGC=IFLAG
      NFE=NFE+1
      RETURN
 50   YP=WP
C IF THE STARTING POINT IS OUT OF RANGE, GIVE UP.
      IF (IFLAG .LE. -100) THEN
        IFLAG=6
        RETURN
      ENDIF
 70   W=Y + H*YP
      Z0=W
      JUDY=1                                    ! DO JUDY=1,LITFH
 80   IF (JUDY > LITFH) GO TO 200
C REQUEST THE CALCULATION OF THE NEWTON STEP  TZ  AT THE CURRENT
C POINT  W  VIA REVERSE CALL.
        IFLAG=IFLAGC-20
        IFLAGC=IFLAG
        NFE=NFE+1
        RETURN
100     IF (IFLAG .LE. -100) GO TO 200
C
C TAKE NEWTON STEP AND CHECK CONVERGENCE.
        W=W + TZ
        ITNUM=JUDY
C COMPUTE QUANTITIES USED FOR OPTIMAL STEP SIZE ESTIMATION.
        IF (JUDY .EQ. 1) THEN
          LCALC=DNRM2(NP1,TZ,1)
          RCALC=RHOLEN
          Z1=W
        ELSE IF (JUDY .EQ. 2) THEN
          LCALC=DNRM2(NP1,TZ,1)/LCALC
          RCALC=RHOLEN/RCALC
        ENDIF
C GO TO MOP-UP SECTION AFTER CONVERGENCE.
        IF (DNRM2(NP1,TZ,1) .LE. RELERR*DNRM2(NP1,W,1)+ABSERR)
     &                                                 GO TO 600
C
        JUDY=JUDY+1
      GO TO 80                                   ! END DO
C
C NO CONVERGENCE IN  LITFH  ITERATIONS.  REDUCE  H  AND TRY AGAIN.
200   IF (H .LE. FOURU*(1.0 + S)) THEN
        IFLAG=6
        RETURN
      ENDIF
      H=.5 * H
      GO TO 70
C
C ***** END OF STARTUP SECTION. *****
C
C ***** PREDICTOR SECTION. *****
C
300   FAIL=.FALSE.
C COMPUTE POINT PREDICTED BY HERMITE INTERPOLANT.  USE STEP SIZE  H
C COMPUTED ON LAST CALL TO  STEPNX .
320   DO J=1,NP1
        W(J)=QOFS(YOLD(J),YPOLD(J),Y(J),YP(J),HOLD,HOLD+H)
      END DO
      Z0=W 
C
C ***** END OF PREDICTOR SECTION. *****
C
C ***** CORRECTOR SECTION. *****
C
      JUDY=1                          ! CORRECTOR: DO JUDY=1,LITFH
350   IF (JUDY > LITFH) GO TO 500
C REQUEST THE CALCULATION OF THE NEWTON STEP  TZ  AT THE CURRENT
C POINT  W  VIA REVERSE CALL.
        IFLAG=IFLAGC-30
        IFLAGC=IFLAG
        NFE=NFE+1
        RETURN
400     IF (IFLAG .LE. -100) GO TO 500
C
C TAKE NEWTON STEP AND CHECK CONVERGENCE.
        W=W + TZ
        ITNUM=JUDY
C COMPUTE QUANTITIES USED FOR OPTIMAL STEP SIZE ESTIMATION.
        IF (JUDY .EQ. 1) THEN
          LCALC=DNRM2(NP1,TZ,1)
          RCALC=RHOLEN
          Z1=W
        ELSE IF (JUDY .EQ. 2) THEN
          LCALC=DNRM2(NP1,TZ,1)/LCALC
          RCALC=RHOLEN/RCALC
        ENDIF
C GO TO MOP-UP SECTION AFTER CONVERGENCE.
        IF (DNRM2(NP1,TZ,1) .LE. RELERR*DNRM2(NP1,W,1)+ABSERR)
     &                                                 GO TO 600
C
        JUDY=JUDY+1
      GO TO 350                              ! END DO CORRECTOR
C
C NO CONVERGENCE IN  LITFH  ITERATIONS.  RECORD FAILURE AT CALCULATED  H , 
C SAVE THIS STEP SIZE, REDUCE  H  AND TRY AGAIN.
500   FAIL=.TRUE.
      HFAIL=H
      IF (H .LE. FOURU*(1.0 + S)) THEN
        IFLAG=6
        RETURN
      ENDIF
      H=.5 * H
      GO TO 320
C
C ***** END OF CORRECTOR SECTION. *****
C
C ***** MOP-UP SECTION. *****
C
C YOLD  AND  Y  ALWAYS CONTAIN THE LAST TWO POINTS FOUND ON THE ZERO
C CURVE OF THE HOMOTOPY MAP.  YPOLD  AND  YP  CONTAIN THE TANGENT
C VECTORS TO THE ZERO CURVE AT  YOLD  AND  Y , RESPECTIVELY.
C
600   YPOLD=YP
      YOLD=Y
      Y=W
      YP=WP
      W=Y - YOLD
C UPDATE ARC LENGTH.
      HOLD=DNRM2(NP1,W,1)
      S=S+HOLD
C
C ***** END OF MOP-UP SECTION. *****
C
C ***** OPTIMAL STEP SIZE ESTIMATION SECTION. *****
C
C CALCULATE THE DISTANCE FACTOR  DCALC .
      TZ=Z0 - Y
      W=Z1 - Y
      DCALC=DNRM2(NP1,TZ,1)
      IF (DCALC .NE. 0.0) DCALC=DNRM2(NP1,W,1)/DCALC
C
C THE OPTIMAL STEP SIZE HBAR IS DEFINED BY
C
C   HT=HOLD * [MIN(LIDEAL/LCALC, RIDEAL/RCALC, DIDEAL/DCALC)]**(1/P)
C
C     HBAR = MIN [ MAX(HT, BMIN*HOLD, HMIN), BMAX*HOLD, HMAX ]
C
C IF CONVERGENCE HAD OCCURRED AFTER 1 ITERATION, SET THE CONTRACTION
C FACTOR  LCALC  TO ZERO.
      IF (ITNUM .EQ. 1) LCALC = 0.0
C FORMULA FOR OPTIMAL STEP SIZE.
      IF (LCALC+RCALC+DCALC .EQ. 0.0) THEN
        HT = SSPAR(7) * HOLD
      ELSE 
        HT = (1.0/MAX(LCALC/SSPAR(1), RCALC/SSPAR(2), DCALC/SSPAR(3)))
     &       **(1.0/SSPAR(8)) * HOLD
      ENDIF
C  HT  CONTAINS THE ESTIMATED OPTIMAL STEP SIZE.  NOW PUT IT WITHIN
C REASONABLE BOUNDS.
      H=MIN(MAX(HT,SSPAR(6)*HOLD,SSPAR(4)), SSPAR(7)*HOLD, SSPAR(5))
      IF (ITNUM .EQ. 1) THEN
C IF CONVERGENCE HAD OCCURRED AFTER 1 ITERATION, DON'T DECREASE  H .
        H=MAX(H,HOLD)
      ELSE IF (ITNUM .EQ. LITFH) THEN
C IF CONVERGENCE REQUIRED THE MAXIMUM  LITFH  ITERATIONS, DON'T
C INCREASE  H .
        H=MIN(H,HOLD)
      ENDIF
C IF CONVERGENCE DID NOT OCCUR IN  LITFH  ITERATIONS FOR A PARTICULAR
C H = HFAIL , DON'T CHOOSE THE NEW STEP SIZE LARGER THAN  HFAIL .
      IF (FAIL) H=MIN(H,HFAIL)
C
C
      IFLAG=IFLAGC
      RETURN
C CLEAN UP ALLOCATED WORKING STORAGE.
 700  START=.TRUE.
      IFLAG=IFLAGC
      IF (ALLOCATED(Z0)) DEALLOCATE(Z0)
      IF (ALLOCATED(Z1)) DEALLOCATE(Z1)
      RETURN
      END SUBROUTINE STEPNX
      SUBROUTINE STEPQF(N,NFE,IFLAG,START,CRASH,HOLD,H,
     &   WK,RELERR,ABSERR,S,Y,YP,YOLD,YPOLD,A,Q,R,
     &   F0,F1,Z0,DZ,W,T,SSPAR)
C
C SUBROUTINE  STEPQF  TAKES ONE STEP ALONG THE ZERO CURVE OF THE 
C HOMOTOPY MAP  RHO(LAMBDA,X)  USING A PREDICTOR-CORRECTOR ALGORITHM.
C THE PREDICTOR USES A HERMITE CUBIC INTERPOLANT, AND THE CORRECTOR 
C RETURNS TO THE ZERO CURVE USING A QUASI-NEWTON ALGORITHM, REMAINING
C IN A HYPERPLANE PERPENDICULAR TO THE MOST RECENT TANGENT VECTOR.
C  STEPQF  ALSO ESTIMATES A STEP SIZE  H  FOR THE NEXT STEP ALONG THE 
C ZERO CURVE.
C
C THE FOLLOWING INTERFACE BLOCK SHOULD BE INCLUDED IN THE CALLING
C PROGRAM:
C
C     INTERFACE
C       SUBROUTINE STEPQF(N,NFE,IFLAG,START,CRASH,HOLD,H,
C    &    WK,RELERR,ABSERR,S,Y,YP,YOLD,YPOLD,A,Q,R,
C    &    F0,F1,Z0,DZ,W,T,SSPAR)
C       USE REAL_PRECISION
C       INTEGER:: N, NFE, IFLAG
C       LOGICAL:: START, CRASH
C       REAL (KIND=R8):: HOLD, H, WK, RELERR, ABSERR, S
C       REAL (KIND=R8):: A(:), DZ(N+1), F0(N+1), F1(N+1), 
C    &    Q(N+1,N+1), R((N+1)*(N+2)/2), SSPAR(4), T(N+1), W(N+1),
C    &    Y(:), YOLD(N+1), YP(N+1), YPOLD(N+1), Z0(N+1)
C       END SUBROUTINE STEPQF
C     END INTERFACE
C
C
C ON INPUT:
C 
C N = DIMENSION OF  X. 
C
C NFE = NUMBER OF JACOBIAN MATRIX EVALUATIONS.
C
C IFLAG = -2, -1, OR 0, INDICATING THE PROBLEM TYPE.
C
C START = .TRUE. ON FIRST CALL TO  STEPQF, .FALSE. OTHERWISE.
C         SHOULD NOT BE MODIFIED BY THE USER AFTER THE FIRST CALL.
C
C HOLD = ||Y - YOLD|| ; SHOULD NOT BE MODIFIED BY THE USER.
C
C H = UPPER LIMIT ON LENGTH OF STEP THAT WILL BE ATTEMPTED.  H  MUST
C    BE SET TO A POSITIVE NUMBER ON THE FIRST CALL TO  STEPQF.
C    THEREAFTER,  STEPQF  CALCULATES AN OPTIMAL VALUE FOR  H, AND  H
C    SHOULD NOT BE MODIFIED BY THE USER.
C
C WK = APPROXIMATE CURVATURE FOR THE LAST STEP (COMPUTED BY PREVIOUS
C    CALL TO  STEPQF).  UNDEFINED ON FIRST CALL.  SHOULD NOT BE
C    MODIFIED BY THE USER.
C  
C RELERR, ABSERR = RELATIVE AND ABSOLUTE ERROR VALUES.  THE ITERATION
C    IS CONSIDERED TO HAVE CONVERGED WHEN A POINT  Z=(LAMBDA,X)  IS 
C    FOUND SUCH THAT
C       ||DZ|| .LE. RELERR*||Z|| + ABSERR,
C    WHERE  DZ  IS THE LAST QUASI-NEWTON STEP.
C
C S  = (APPROXIMATE) ARC LENGTH ALONG THE HOMOTOPY ZERO CURVE UP TO
C    Y(S) = (LAMBDA(S), X(S)).
C
C Y(1:N+1) = PREVIOUS POINT (LAMBDA(S),X(S)) FOUND ON THE ZERO CURVE
C    OF THE HOMOTOPY MAP.
C
C YP(1:N+1) = UNIT TANGENT VECTOR TO THE ZERO CURVE OF THE HOMOTOPY
C    MAP AT  Y.  INPUT IN THIS VECTOR IS NOT USED ON THE FIRST CALL
C    TO  STEPQF.
C
C YOLD(1:N+1) = A POINT BEFORE  Y  ON THE ZERO CURVE OF THE HOMOTOPY
C    MAP.  INPUT IN THIS VECTOR IS NOT USED ON THE FIRST CALL TO 
C    STEPQF.
C
C YPOLD(1:N+1) = UNIT TANGENT VECTOR TO THE ZERO CURVE OF THE 
C    HOMOTOPY MAP AT  YOLD.
C
C A(:) = PARAMETER VECTOR IN THE HOMOTOPY MAP.
C
C Q(1:N+1,1:N+1) =  Q  OF THE QR FACTORIZATION OF
C    THE AUGMENTED JACOBIAN MATRIX AT  Y.
C
C R((N+1)*(N+2)/2) = THE UPPER TRIANGLE  R  OF THE QR 
C    FACTORIZATION, STORED BY COLUMNS.
C
C F0(1:N+1), F1(1:N+1), Z0(1:N+1), DZ(1:N+1), W(1:N+1), T(1:N+1) ARE
C    WORK ARRAYS.  
C 
C SSPAR(1:4) = PARAMETERS USED FOR COMPUTATION OF THE OPTIMAL STEP SIZE.  
C    SSPAR(1) = HMIN, SSPAR(2) = HMAX, SSPAR(3) = BMIN, SSPAR(4) = BMAX.  
C    THE OPTIMAL STEP  H  IS RESTRICTED SUCH THAT 
C       HMIN .LE. H .LE. HMAX, AND  BMIN*HOLD .LE. H .LE. BMAX*HOLD.
C
C
C ON OUTPUT:
C
C NFE HAS BEEN UPDATED.
C
C IFLAG
C
C    = -2, -1, OR 0 (UNCHANGED) ON A NORMAL RETURN.
C
C    = 4 IF A JACOBIAN MATRIX WITH RANK <  N  HAS OCCURRED.  THE
C        ITERATION WAS NOT COMPLETED.
C
C    = 6 IF THE ITERATION FAILED TO CONVERGE. 
C
C START = .FALSE. ON A NORMAL RETURN.
C
C CRASH 
C
C    = .FALSE. ON A NORMAL RETURN.
C
C    = .TRUE. IF THE STEP SIZE  H  WAS TOO SMALL.  H  HAS BEEN
C      INCREASED TO AN ACCEPTABLE VALUE, WITH WHICH  STEPQF  MAY BE
C      CALLED AGAIN.
C
C    = .TRUE. IF  RELERR  AND/OR  ABSERR  WERE TOO SMALL.  THEY HAVE
C      BEEN INCREASED TO ACCEPTABLE VALUES, WITH WHICH  STEPQF  MAY
C      BE CALLED AGAIN.
C
C HOLD = ||Y-YOLD||.
C
C H = OPTIMAL VALUE FOR NEXT STEP TO BE ATTEMPTED.  NORMALLY  H  SHOULD
C     NOT BE MODIFIED BY THE USER.
C
C WK = APPROXIMATE CURVATURE FOR THE STEP TAKEN BY  STEPQF.
C
C S = (APPROXIMATE) ARC LENGTH ALONG THE ZERO CURVE OF THE HOMOTOPY 
C     MAP UP TO THE LATEST POINT FOUND, WHICH IS RETURNED IN  Y.
C
C RELERR, ABSERR  ARE UNCHANGED ON A NORMAL RETURN.  THEY ARE POSSIBLY
C     CHANGED IF  CRASH  = .TRUE. (SEE DESCRIPTION OF  CRASH  ABOVE).
C
C Y, YP, YOLD, YPOLD  CONTAIN THE TWO MOST RECENT POINTS AND TANGENT
C     VECTORS FOUND ON THE ZERO CURVE OF THE HOMOTOPY MAP.
C
C Q, R  STORE THE QR FACTORIZATION OF THE AUGMENTED JACOBIAN MATRIX 
C     EVALUATED AT  Y.
C
C
C CALLS  DGEMV, DGEQRF, DNRM2, DORGQR, DTPSV, F (OR RHO),
C     FJAC (OR RHOJAC), TANGQF, UPQRQF.
C
C ***** DECLARATIONS *****
      USE HOMOTOPY
      USE REAL_PRECISION
C
C     FUNCTION DECLARATIONS  
C
      INTERFACE
        FUNCTION DNRM2(N,X,STRIDE)
        USE REAL_PRECISION
        INTEGER:: N,STRIDE
        REAL (KIND=R8):: DNRM2,X(N)
        END FUNCTION DNRM2
      END INTERFACE
      REAL (KIND=R8):: DD001, DD0011, DD01, DD011, QOFS
C
C     LOCAL VARIABLES
C
      REAL (KIND=R8):: ALPHA, DELS, ETA, FOURU, GAMMA, HFAIL, HTEMP,
     &  IDLERR, ONE, P0, P1, PP0, PP1, TEMP, TWOU, WKOLD, ZERO         
      INTEGER:: I, ITCNT, LITFH, J, JP1, NP1
      LOGICAL:: FAILED
C
C     SCALAR ARGUMENTS 
C
      INTEGER:: N, NFE, IFLAG
      LOGICAL:: START, CRASH
      REAL (KIND=R8):: HOLD, H, WK, RELERR, ABSERR, S
C
C     ARRAY DECLARATIONS
C
      REAL (KIND=R8):: A(:), DZ(N+1), F0(N+1), F1(N+1), 
     &   Q(N+1,N+1), R((N+1)*(N+2)/2), SSPAR(4), T(N+1), W(N+1),
     &   Y(:), YOLD(N+1), YP(N+1), YPOLD(N+1), Z0(N+1)
C
      SAVE
C
C ***** END OF DECLARATIONS *****
C
      INTERFACE
        SUBROUTINE TANGQF(Y,YP,YPOLD,A,Q,R,W,S,T,N,IFLAG,NFE)
        USE HOMOTOPY
        USE REAL_PRECISION
        INTEGER:: N, IFLAG, NFE
        REAL (KIND=R8):: A(:), Q(N+1,N+1), R((N+1)*(N+2)/2),
     &    S(N+1), T(N+1), W(N+1), Y(:), YP(N+1), YPOLD(N+1)
        END SUBROUTINE TANGQF
      END INTERFACE
C
C DEFINITION OF HERMITE CUBIC INTERPOLANT VIA DIVIDED DIFFERENCES.
C
      DD01(P0,P1,DELS) = (P1-P0)/DELS
      DD001(P0,PP0,P1,DELS) = (DD01(P0,P1,DELS)-PP0)/DELS
      DD011(P0,P1,PP1,DELS) = (PP1-DD01(P0,P1,DELS))/DELS
      DD0011(P0,PP0,P1,PP1,DELS) = (DD011(P0,P1,PP1,DELS) -
     &                                DD001(P0,PP0,P1,DELS))/DELS
      QOFS(P0,PP0,P1,PP1,DELS,S) = ((DD0011(P0,PP0,P1,PP1,DELS)*
     &    (S-DELS) + DD001(P0,PP0,P1,DELS))*S + PP0)*S + P0
C
C ***** FIRST EXECUTABLE STATEMENT *****
C
C
C ***** INITIALIZATION *****
C
C ETA = PARAMETER FOR BROYDEN'S UPDATE.
C LITFH = MAXIMUM NUMBER OF QUASI-NEWTON ITERATIONS ALLOWED.
C
      ONE = 1.0
      ZERO = 0.0
      TWOU = 2.0*EPSILON(1.0_R8)
      FOURU = TWOU + TWOU
      NP1 = N+1
      FAILED = .FALSE.
      CRASH = .TRUE.
      ETA = 50.0*TWOU
      LITFH = 2*(INT(ABS(LOG10(ABSERR+RELERR)))+1)
C
C CHECK THAT ALL INPUT PARAMETERS ARE CORRECT.
C
C     THE ARCLENGTH  S MUST BE NONNEGATIVE.
C
      IF (S .LT. 0.0) RETURN
C
C     IF STEP SIZE IS TOO SMALL, DETERMINE AN ACCEPTABLE ONE.
C   
      IF (H .LT. FOURU*(1.0+S)) THEN
        H=FOURU*(1.0 + S)
        RETURN
      END IF
C
C     IF ERROR TOLERANCES ARE TOO SMALL, INCREASE THEM TO ACCEPTABLE 
C     VALUES.
C
      TEMP=DNRM2(NP1,Y,1) + 1.0
      IF (.5*(RELERR*TEMP+ABSERR) .LT. TWOU*TEMP) THEN
        IF (RELERR .NE. 0.0) THEN
          RELERR = FOURU*(1.0+FOURU)
          TEMP = 0.0
          ABSERR = MAX(ABSERR,TEMP)
        ELSE
          ABSERR=FOURU*TEMP
        END IF
        RETURN
      END IF
C
C     INPUT PARAMETERS WERE ALL ACCEPTABLE.
C
      CRASH = .FALSE.
C
C COMPUTE  YP  ON FIRST CALL.
C NOTE:  DZ  IS USED SIMPLY AS A WORK ARRAY HERE.
C
      IF (START) THEN
        CALL TANGQF(Y,YP,YPOLD,A,Q,R,W,DZ,T,N,IFLAG,NFE)
        IF (IFLAG .GT. 0) RETURN
      END IF
C
C F0 = (RHO(Y), YP*Y) TRANSPOSE (DIFFERENT FOR EACH PROBLEM TYPE).
C
      IF (IFLAG .EQ. -2) THEN
C
C CURVE TRACKING PROBLEM.
C
        CALL RHO(A,Y(1),Y(2:NP1),F0(1:N))
      ELSE IF (IFLAG .EQ. -1) THEN
C
C ZERO FINDING PROBLEM.
C
        CALL F(Y(2:NP1),F0(1:N))
        F0(1:N) = Y(1)*F0(1:N) + (1.0-Y(1))*(Y(2:NP1)-A(1:N))
      ELSE
C
C FIXED POINT PROBLEM.
C
        CALL F(Y(2:NP1),F0(1:N))
        F0(1:N) = Y(1)*(A(1:N)-F0(1:N))+Y(2:NP1)-A(1:N)
      END IF
C
C DEFINE LAST ROW OF F0  =  YP*Y.
C
       F0(NP1) = DOT_PRODUCT(YP,Y)
C
C ***** END OF INITIALIZATION *****
C
C ***** COMPUTE PREDICTOR POINT Z0 *****
C
 20   IF (START) THEN
C           
C COMPUTE Z0 WITH LINEAR PREDICTOR USING Y, YP --
C Z0 = Y+H*YP.
C
        Z0 = Y + H*YP
C         
      ELSE
C
C COMPUTE Z0 WITH CUBIC PREDICTOR.
C
        DO I=1,NP1
          Z0(I) = QOFS(YOLD(I),YPOLD(I),Y(I),YP(I),HOLD,HOLD+H) 
        END DO
C
      END IF
C
C F1 = (RHO(Z0), YP*Z0) TRANSPOSE.
C
      IF (IFLAG .EQ. -2) THEN
        CALL RHO(A,Z0(1),Z0(2:NP1),F1(1:N))
      ELSE IF (IFLAG .EQ. -1) THEN
        CALL F(Z0(2:NP1),F1(1:N))
        F1(1:N) = Z0(1)*F1(1:N) + (1.0-Z0(1))*(Z0(2:NP1)-A(1:N))
      ELSE
        CALL F(Z0(2:NP1),F1(1:N))
        F1(1:N) = Z0(1)*(A(1:N)-F1(1:N))+Z0(2:NP1)-A(1:N)
      END IF
      F1(NP1) = DOT_PRODUCT(YP,Z0)
C
C ***** END OF PREDICTOR SECTION *****
C
C ***** SET-UP FOR QUASI-NEWTON ITERATION *****
C
      IF (FAILED) THEN
C        
C GENERATE Q = AUGMENTED JACOBIAN MATRIX FOR POINT Z0=(LAMBDA,X).
C        
        IF (IFLAG .EQ. -2) THEN
C
C CURVE TRACKING PROBLEM:
C D(RHO) = (D RHO(A,LAMBDA,X)/D LAMBDA, D RHO(A,LAMBDA,X)/DX).
C
          DO J = 1,NP1
            CALL RHOJAC(A,Z0(1),Z0(2:NP1),Q(1:N,J),J)
          END DO
        ELSE IF (IFLAG .EQ. -1) THEN
C
C ZERO FINDING PROBLEM:
C D(RHO) = (F(X) - X + A, LAMBDA*DF(X) + (1-LAMBDA)*I).
C
          CALL F(Z0(2:NP1),Q(1:N,1))
          Q(1:N,1) = A - Z0(2:NP1) + Q(1:N,1)
          DO J= 1,N
            JP1 = J+1
            CALL FJAC(Z0(2:NP1),Q(1:N,JP1),J)
            Q(1:N,JP1) = Z0(1)*Q(1:N,JP1)
            Q(J,JP1) = 1.0 - Z0(1) + Q(J,JP1)
          END DO
        ELSE 
C 
C FIXED POINT PROBLEM:
C D(RHO) = (A - F(X), I - LAMBDA*DF(X)).
C
          CALL F(Z0(2:NP1),Q(1:N,1))
          Q(1:N,1) = A - Q(1:N,1)
          DO J=1,N
            JP1 = J+1
            CALL FJAC(Z0(2:NP1),Q(1:N,JP1),J)
            Q(1:N,JP1) = -Z0(1)*Q(1:N,JP1)
            Q(J,JP1) = 1.0 + Q(J,JP1)
          END DO
        END IF
C
C DEFINE LAST ROW OF Q = YP.
C
        Q(NP1,:) = YP
C
C COUNT JACOBIAN EVALUATION.
C
        NFE = NFE+1
C
C DO FIRST QUASI-NEWTON STEP.
C
C FACTOR AUG.
C
        CALL DGEQRF(NP1,NP1,Q,NP1,T,W,NP1,I)
C
C PACK UPPER TRIANGLE INTO ARRAY R.
C
        DO I=1,NP1
          R((I*(I-1))/2 + 1:(I*(I-1))/2 + I) = Q(1:I,I)
        END DO
C
C CHECK FOR SINGULARITY.
C
        J = 1
        DO I = 1, N
          IF( R(J+I-1) .EQ. ZERO ) THEN
            IFLAG = 4
            RETURN
          END IF
          J = J + I
        END DO
C
C EXPAND HOUSEHOLDER REFLECTIONS INTO FULL MATRIX Q .
C
        CALL DORGQR(NP1, NP1, N, Q, NP1, T, W, NP1, I)
C
C COMPUTE NEWTON STEP.
C
        T(1:N) = -F1(1:N)
        T(NP1) = 0.0
        CALL DGEMV('T',NP1,NP1,ONE,Q,NP1,T,1,ZERO,DZ,1)
        CALL DTPSV('U', 'N', 'N', NP1, R, DZ, 1)
C
C TAKE STEP AND SET F0 = F1.
C
        Z0 = Z0 + DZ
        F0 = F1
C
C F1 = (RHO(Z0), YP*Z0) TRANSPOSE.
C
        IF (IFLAG .EQ. -2) THEN
          CALL RHO(A,Z0(1),Z0(2:NP1),F1(1:N))
        ELSE IF (IFLAG .EQ. -1) THEN
          CALL F(Z0(2:NP1),F1(1:N))
          F1(1:N) = Z0(1)*F1(1:N) + (1.0-Z0(1))*(Z0(2:NP1)-A(1:N))
        ELSE
          CALL F(Z0(2:NP1),F1(1:N))
          F1(1:N) = Z0(1)*(A(1:N)-F1(1:N))+Z0(2:NP1)-A(1:N)
        END IF
        F1(NP1) = DOT_PRODUCT(YP,Z0)
C
      ELSE
C
C IF NOT FAILED THEN DEFINE  DZ=Z0-Y  PRIOR TO MAIN LOOP.
C
        DZ = Z0 - Y
      END IF
C
C ***** END OF PREPARATION FOR QUASI-NEWTON ITERATION *****
C
      DO ITCNT = 1,LITFH  ! ***** QUASI-NEWTON ITERATION *****
C
C PERFORM UPDATE FOR NEWTON STEP JUST TAKEN.
C
        CALL UPQRQF(NP1,ETA,DZ,F0,F1,Q,R,W,T)
C
C COMPUTE NEXT NEWTON STEP.
C
        T(1:N) = -F1(1:N)
        T(NP1) = 0.0
        CALL DGEMV('T',NP1,NP1,ONE,Q,NP1,T,1,ZERO,DZ,1)
        CALL DTPSV('U', 'N', 'N', NP1, R, DZ, 1)
C
C TAKE STEP.
C
        Z0 = Z0 + DZ
C
C CHECK FOR CONVERGENCE.
C
        IF (DNRM2(NP1,DZ,1) .LE. RELERR*DNRM2(NP1,Z0,1)+ABSERR) THEN
           GO TO 180
        END IF
C
C IF NOT CONVERGED, PREPARE FOR NEXT ITERATION.
C
        F0 = F1
C
C F1 = (RHO(Z0), YP*Z0) TRANSPOSE.
C
        IF (IFLAG .EQ. -2) THEN
          CALL RHO(A,Z0(1),Z0(2:NP1),F1(1:N))
        ELSE IF (IFLAG .EQ. -1) THEN
          CALL F(Z0(2:NP1),F1(1:N))
          F1(1:N) = Z0(1)*F1(1:N) + (1.0-Z0(1))*(Z0(2:NP1)-A(1:N))
        ELSE
          CALL F(Z0(2:NP1),F1(1:N))
          F1(1:N) = Z0(1)*(A(1:N)-F1(1:N))+Z0(2:NP1)-A(1:N)
        END IF
        F1(NP1) = DOT_PRODUCT(YP,Z0)
C
      END DO  ! ***** END OF QUASI-NEWTON LOOP *****
C
C ***** DIDN'T CONVERGE OR TANGENT AT NEW POINT DID NOT MAKE
C       AN ACUTE ANGLE WITH YPOLD -- TRY AGAIN WITH A SMALLER H *****
C      
 170  FAILED = .TRUE.
      HFAIL = H
      IF (H .LE. FOURU*(1.0 + S)) THEN
        IFLAG = 6
        RETURN
      ELSE
        H = .5 * H
      END IF
      GO TO 20
C
C ***** END OF CONVERGENCE FAILURE SECTION *****
C
C ***** CONVERGED -- MOP UP AND RETURN *****
C
C COMPUTE TANGENT & AUGMENTED JACOBIAN AT  Z0.
C NOTE:  DZ  AND  F1  ARE USED SIMPLY AS WORK ARRAYS HERE.
C
 180  CALL TANGQF(Z0,T,YP,A,Q,R,W,DZ,F1,N,IFLAG,NFE)
      IF (IFLAG .GT. 0) RETURN
C
C CHECK THAT COMPUTED TANGENT  T  MAKES AN ANGLE NO LARGER THAN
C 60 DEGREES WITH CURRENT TANGENT  YP.  (I.E. COS OF ANGLE < .5)
C IF NOT, STEP SIZE WAS TOO LARGE, SO THROW AWAY Z0, AND TRY
C AGAIN WITH A SMALLER STEP.
C
      ALPHA = DOT_PRODUCT(T,YP)
      IF (ALPHA .LT. 0.5) GOTO 170
      ALPHA = ACOS(ALPHA)
C
C SET UP VARIABLES FOR NEXT CALL.
C
      YOLD = Y
      Y = Z0
      YPOLD = YP
      YP = T
C
C UPDATE ARCLENGTH   S = S + ||Y-YOLD||.
C
      HTEMP = HOLD
      Z0 = Z0 - YOLD
      HOLD = DNRM2(NP1,Z0,1)
      S = S+HOLD
C
C COMPUTE OPTIMAL STEP SIZE. 
C   IDLERR = DESIRED ERROR FOR NEXT PREDICTOR STEP.
C   WK = APPROXIMATE CURVATURE = 2*SIN(ALPHA/2)/HOLD  WHERE 
C        ALPHA = ARCCOS(YP*YPOLD).
C   GAMMA = EXPECTED CURVATURE FOR NEXT STEP, COMPUTED BY 
C        EXTRAPOLATING FROM CURRENT CURVATURE  WK, AND LAST 
C        CURVATURE  WKOLD.  GAMMA IS FURTHER REQUIRED TO BE 
C        POSITIVE.
C
      IF (.NOT. START) WKOLD = WK
      IDLERR = SQRT(SQRT(ABSERR + RELERR*DNRM2(NP1,Y,1)))
C
C     IDLERR SHOULD BE NO BIGGER THAN 1/2 PREVIOUS STEP.
C
      IDLERR = MIN(.5*HOLD,IDLERR)
      WK = 2.0*ABS(SIN(.5*ALPHA))/HOLD
      IF (START) THEN
         GAMMA = WK
      ELSE 
         GAMMA = WK + HOLD/(HOLD+HTEMP)*(WK-WKOLD)
      END IF
      GAMMA = MAX(GAMMA, 0.01*ONE)
      H = SQRT(2.0*IDLERR/GAMMA)
C
C     ENFORCE RESTRICTIONS ON STEP SIZE SO AS TO ENSURE STABILITY.
C        HMIN <= H <= HMAX, BMIN*HOLD <= H <= BMAX*HOLD.
C
      H = MIN(MAX(SSPAR(1),SSPAR(3)*HOLD,H),SSPAR(4)*HOLD,SSPAR(2))
      IF (FAILED) H = MIN(HFAIL,H)
      START = .FALSE.
C
C ***** END OF MOP UP SECTION *****
C
      RETURN
C
      END SUBROUTINE STEPQF
      SUBROUTINE STEPQS(N,NFE,IFLAG,MODE,LENQR,START,CRASH,HOLD,H,
     &  WK,RELERR,ABSERR,S,Y,YP,YOLD,YPOLD,A,Z0,DZ,T,SSPAR)
C
C SUBROUTINE  STEPQS  TAKES ONE STEP ALONG THE ZERO CURVE OF THE 
C HOMOTOPY MAP  RHO(X,LAMBDA)  USING A PREDICTOR-CORRECTOR ALGORITHM.
C THE PREDICTOR USES A HERMITE CUBIC INTERPOLANT, AND THE CORRECTOR 
C RETURNS TO THE ZERO CURVE USING A NEWTON ITERATION, REMAINING
C IN A HYPERPLANE PERPENDICULAR TO THE MOST RECENT TANGENT VECTOR.
C  STEPQS  ALSO ESTIMATES A STEP SIZE  H  FOR THE NEXT STEP ALONG THE 
C ZERO CURVE.  SEE ALSO THE REVERSE CALL ROUTINE  STEPNX .
C 
C THE CALLING PROGRAM MUST CONTAIN THE FOLLOWING INTERFACE BLOCK:
C
C     INTERFACE
C       SUBROUTINE STEPQS(N,NFE,IFLAG,MODE,LENQR,START,CRASH,HOLD,H,
C    &    WK,RELERR,ABSERR,S,Y,YP,YOLD,YPOLD,A,Z0,DZ,T,SSPAR)
C       USE REAL_PRECISION
C       INTEGER, INTENT(IN):: LENQR,MODE,N
C       INTEGER, INTENT(IN OUT):: IFLAG,NFE
C       LOGICAL, INTENT(IN OUT):: CRASH,START
C       REAL (KIND=R8), INTENT(IN):: A(:),SSPAR(4)
C       REAL (KIND=R8), INTENT(IN OUT):: ABSERR,H,HOLD,RELERR,S,WK,
C    &    Y(:),YOLD(:),YP(:),YPOLD(:)
C       REAL (KIND=R8), INTENT(OUT), DIMENSION(:):: DZ,T,Z0
C       END SUBROUTINE STEPQS
C     END INTERFACE
C
C
C ON INPUT:
C 
C N = DIMENSION OF  X. 
C
C NFE = NUMBER OF JACOBIAN MATRIX EVALUATIONS.
C
C IFLAG = -2, -1, OR 0, INDICATING THE PROBLEM TYPE.
C
C MODE = 1 IF THE JACOBIAN MATRIX IS SYMMETRIC AND STORED IN A PACKED
C          SKYLINE FORMAT;
C      = 2 IF THE JACOBIAN MATRIX IS STORED IN A SPARSE ROW FORMAT.
C
C LENQR  IS THE NUMBER OF NONZERO ENTRIES IN THE SPARSE JACOBIAN
C    MATRICES, USED TO DETERMINE THE SPARSE MATRIX DATA STRUCTURES.
C
C START = .TRUE. ON FIRST CALL TO  STEPQS, .FALSE. OTHERWISE.
C         SHOULD NOT BE MODIFIED BY THE USER AFTER THE FIRST CALL.
C
C HOLD = ||Y - YOLD|| ; SHOULD NOT BE MODIFIED BY THE USER.
C
C H = UPPER LIMIT ON LENGTH OF STEP THAT WILL BE ATTEMPTED.  H  MUST
C    BE SET TO A POSITIVE NUMBER ON THE FIRST CALL TO  STEPQS.
C    THEREAFTER,  STEPQS  CALCULATES AN OPTIMAL VALUE FOR  H, AND  H
C    SHOULD NOT BE MODIFIED BY THE USER.
C
C WK = APPROXIMATE CURVATURE FOR THE LAST STEP (COMPUTED BY PREVIOUS
C    CALL TO  STEPQS).  UNDEFINED ON FIRST CALL.  SHOULD NOT BE
C    MODIFIED BY THE USER.
C  
C RELERR, ABSERR = RELATIVE AND ABSOLUTE ERROR VALUES.  THE ITERATION
C    IS CONSIDERED TO HAVE CONVERGED WHEN A POINT  Z=(X,LAMBDA)  IS 
C    FOUND SUCH THAT
C       ||DZ|| .LE. RELERR*||Z|| + ABSERR,
C    WHERE  DZ  IS THE LAST NEWTON STEP.
C
C S  = (APPROXIMATE) ARC LENGTH ALONG THE HOMOTOPY ZERO CURVE UP TO
C    Y(S) = (X(S),LAMBDA(S)).
C
C Y(1:N+1) = PREVIOUS POINT (X(S),LAMBDA(S)) FOUND ON THE ZERO CURVE
C    OF THE HOMOTOPY MAP.
C
C YP(1:N+1) = UNIT TANGENT VECTOR TO THE ZERO CURVE OF THE HOMOTOPY
C    MAP AT  Y.  INPUT IN THIS VECTOR IS NOT USED ON THE FIRST CALL
C    TO  STEPQS.
C
C YOLD(1:N+1) = A POINT BEFORE  Y  ON THE ZERO CURVE OF THE HOMOTOPY
C    MAP.  INPUT IN THIS VECTOR IS NOT USED ON THE FIRST CALL TO 
C    STEPQS.
C
C YPOLD(1:N+1) = UNIT TANGENT VECTOR TO THE ZERO CURVE OF THE 
C    HOMOTOPY MAP AT  YOLD.
C
C A(:) = PARAMETER VECTOR IN THE HOMOTOPY MAP.
C
C QR(1:LENQR), PP(1:N), ROWPOS(1:N+2), COLPOS(1:LENQR) ARE ALL WORK
C    ARRAYS USED TO DEFINE THE SPARSE JACOBIAN MATRICES, ALLOCATED
C    IN FIXPQS, AND DISTRIBUTED VIA THE MODULE  HOMOTOPY .
C
C Z0(1:N+1), DZ(1:N+1), T(1:N+1)  ARE ALL WORK ARRAYS USED TO
C    CALCULATE THE TANGENT VECTORS AND NEWTON STEPS.
C    
C SSPAR(1:4) = PARAMETERS USED FOR COMPUTATION OF THE OPTIMAL STEP SIZE.
C    SSPAR(1) = HMIN, SSPAR(2) = HMAX, SSPAR(3) = BMIN, SSPAR(4) = BMAX.
C    THE OPTIMAL STEP  H  IS RESTRICTED SUCH THAT 
C       HMIN .LE. H .LE. HMAX, AND  BMIN*HOLD .LE. H .LE. BMAX*HOLD.
C
C
C ON OUTPUT:
C
C N, LENQR, A  ARE UNCHANGED.
C
C NFE HAS BEEN UPDATED.
C
C IFLAG
C
C    = -2, -1, OR 0 (UNCHANGED) ON A NORMAL RETURN.
C
C    = 4 IF A JACOBIAN MATRIX WITH RANK <  N  HAS OCCURRED.  THE
C        ITERATION WAS NOT COMPLETED.
C
C    = 6 IF THE ITERATION FAILED TO CONVERGE. 
C
C START = .FALSE. ON A NORMAL RETURN.
C
C CRASH 
C
C    = .FALSE. ON A NORMAL RETURN.
C
C    = .TRUE. IF THE STEP SIZE  H  WAS TOO SMALL.  H  HAS BEEN
C      INCREASED TO AN ACCEPTABLE VALUE, WITH WHICH  STEPQS  MAY BE
C      CALLED AGAIN.
C
C    = .TRUE. IF  RELERR  AND/OR  ABSERR  WERE TOO SMALL.  THEY HAVE
C      BEEN INCREASED TO ACCEPTABLE VALUES, WITH WHICH  STEPQS  MAY
C      BE CALLED AGAIN.
C
C HOLD = ||Y-YOLD||.
C
C H = OPTIMAL VALUE FOR NEXT STEP TO BE ATTEMPTED.  NORMALLY  H  SHOULD
C     NOT BE MODIFIED BY THE USER.
C
C WK = APPROXIMATE CURVATURE FOR THE STEP TAKEN BY  STEPQS.
C
C S = (APPROXIMATE) ARC LENGTH ALONG THE ZERO CURVE OF THE HOMOTOPY 
C     MAP UP TO THE LATEST POINT FOUND, WHICH IS RETURNED IN  Y.
C
C RELERR, ABSERR  ARE UNCHANGED ON A NORMAL RETURN.  THEY ARE POSSIBLY
C     CHANGED IF  CRASH  = .TRUE. (SEE DESCRIPTION OF  CRASH  ABOVE).
C
C Y, YP, YOLD, YPOLD  CONTAIN THE TWO MOST RECENT POINTS AND TANGENT
C     VECTORS FOUND ON THE ZERO CURVE OF THE HOMOTOPY MAP.
C
C
C CALLS  DNRM2, TANGNS.
C
      USE REAL_PRECISION
      INTEGER, INTENT(IN):: LENQR,MODE,N
      INTEGER, INTENT(IN OUT):: IFLAG,NFE
      LOGICAL, INTENT(IN OUT):: CRASH,START
      REAL (KIND=R8), INTENT(IN):: A(:),SSPAR(4)
      REAL (KIND=R8), INTENT(IN OUT):: ABSERR,H,HOLD,RELERR,S,WK,
     &    Y(:),YOLD(:),YP(:),YPOLD(:)
      REAL (KIND=R8), INTENT(OUT), DIMENSION(:):: DZ,T,Z0
C
C     FUNCTION DECLARATIONS.  
C
      INTERFACE
        FUNCTION DNRM2(N,X,STRIDE)
        USE REAL_PRECISION
        INTEGER:: N,STRIDE
        REAL (KIND=R8):: DNRM2,X(N)
        END FUNCTION DNRM2
      END INTERFACE
      REAL (KIND=R8):: DD001,DD0011,DD01,DD011,QOFS
C
C     LOCAL VARIABLES.
C
      REAL (KIND=R8), SAVE:: ACOF(12), ALPHA, CORDIS, DELS, FOURU,
     &  GAMMA, HFAIL, HTEMP, IDLERR, OMEGA, P0, P1, PP0, PP1, 
     &  SIGMA, TEMP, THETA, TWOU, WKOLD, WRGE(8), XSTEP
      INTEGER:: I, ITCNT, LK, LST, NP1
      LOGICAL:: FAILED
      DATA WRGE  /
     &   .8735115E+00_R8, .1531947E+00_R8, .3191815E-01_R8,
     &   .3339946E-10_R8, .4677788E+00_R8, .6970123E-03_R8,
     &   .1980863E-05_R8, .1122789E-08_R8/
      DATA ACOF  /
     &   .9043128E+00_R8, -.7075675E+00_R8, -.4667383E+01_R8,
     &  -.3677482E+01_R8,  .8516099E+00_R8, -.1953119E+00_R8,
     &  -.4830636E+01_R8, -.9770528E+00_R8,  .1040061E+01_R8,
     &   .3793395E-01_R8,  .1042177E+01_R8,  .4450706E-01_R8/
C
      INTERFACE
        SUBROUTINE TANGNS(RHOLEN,Y,YP,TZ,YPOLD,A,MODE,LENQR,
     &    NFE,N,IFLAG)
        USE REAL_PRECISION
        REAL (KIND=R8), INTENT(IN), DIMENSION(:):: A,Y,YPOLD
        REAL (KIND=R8), INTENT(IN OUT):: RHOLEN
        REAL (KIND=R8), INTENT(OUT), DIMENSION(:):: TZ,YP
        INTEGER, INTENT(IN):: LENQR,MODE,N
        INTEGER, INTENT(IN OUT):: IFLAG,NFE
        END SUBROUTINE TANGNS
      END INTERFACE
C
C THE LIMIT ON THE NUMBER OF NEWTON ITERATIONS ALLOWED BEFORE REDUCING
C THE STEP SIZE  H  MAY BE CHANGED BY CHANGING THE FOLLOWING PARAMETER 
C STATEMENT:
      INTEGER, PARAMETER:: LITFH = 10
C
C DEFINITION OF HERMITE CUBIC INTERPOLANT VIA DIVIDED DIFFERENCES.
C
      DD01(P0,P1,DELS) = (P1-P0)/DELS
      DD001(P0,PP0,P1,DELS) = (DD01(P0,P1,DELS)-PP0)/DELS
      DD011(P0,P1,PP1,DELS) = (PP1-DD01(P0,P1,DELS))/DELS
      DD0011(P0,PP0,P1,PP1,DELS) = (DD011(P0,P1,PP1,DELS) -
     &  DD001(P0,PP0,P1,DELS))/DELS
      QOFS(P0,PP0,P1,PP1,DELS,S) = ((DD0011(P0,PP0,P1,PP1,DELS)*
     &  (S-DELS) + DD001(P0,PP0,P1,DELS))*S + PP0)*S + P0
C
C ***** END OF SPECIFICATION SECTION. *****
C
C ***** INITIALIZATION. *****
C
      TWOU = 2.0*EPSILON(1.0_R8)
      FOURU = TWOU + TWOU
      NP1 = N+1
      FAILED = .FALSE.
      CRASH = .TRUE.
C 
C CHECK THAT ALL INPUT PARAMETERS ARE CORRECT.
C
C     THE ARCLENGTH  S  MUST BE NONNEGATIVE.
C
      IF (S .LT. 0.0) RETURN
C
C     IF STEP SIZE IS TOO SMALL, DETERMINE AN ACCEPTABLE ONE.
C   
      IF (H .LT. FOURU*(1.0+S)) THEN
          H=FOURU*(1.0 + S)
          RETURN
      END IF
C
C     IF ERROR TOLERANCES ARE TOO SMALL, INCREASE THEM TO ACCEPTABLE 
C     VALUES.
C
      TEMP=DNRM2(NP1,Y,1) + 1.0
      IF (.5*(RELERR*TEMP+ABSERR) .LT. TWOU*TEMP) THEN
          IF (RELERR .NE. 0.0) THEN
            RELERR = FOURU*(1.0+FOURU)
            TEMP = 0.0
            ABSERR = MAX(ABSERR,TEMP)
          ELSE
            ABSERR=FOURU*TEMP
          END IF
          RETURN
      END IF
C
C     INPUT PARAMETERS WERE ALL ACCEPTABLE.
C
      CRASH = .FALSE.
C
C COMPUTE  YP  ON FIRST CALL.
C
      IF (START) THEN
C
C         INITIALIZE THE IDEAL ERROR USED FOR STEP SIZE ESTIMATION.
C
          IDLERR=SQRT(SQRT(ABSERR))
C
          CALL TANGNS(S,Y,YP,DZ,YPOLD,A,MODE,LENQR,NFE,N,IFLAG)
          IF (IFLAG .GT. 0) RETURN
      END IF
C
      CONV: DO
C
C ***** COMPUTE PREDICTOR POINT Z0. *****
C
        IF (START) THEN
C           
C         COMPUTE Z0 WITH LINEAR PREDICTOR USING Y, YP --
C         
          Z0 = Y + H*YP
        ELSE
C
C         COMPUTE Z0 WITH CUBIC PREDICTOR.
C
          DO I=1,NP1
            Z0(I) = QOFS(YOLD(I),YPOLD(I),Y(I),YP(I),HOLD,HOLD+H) 
          END DO
        END IF
C
C ***** END OF PREDICTOR SECTION. *****
C
        NEWTON: DO ITCNT = 1,LITFH   ! ***** NEWTON ITERATION. *****
C
C COMPUTE TANGENT  T  AND MINIMUM NORM NEWTON STEP  DZ  AT THE
C CURRENT POINT  Z0 .
C
          TEMP = -1.0
          CALL TANGNS(TEMP,Z0,T,DZ,YP,A,MODE,LENQR,NFE,N,IFLAG)
          IF (IFLAG .GT. 0) RETURN
C
C CHECK THAT COMPUTED TANGENT  T  MAKES AN ANGLE NO LARGER THAN
C 60 DEGREES WITH CURRENT TANGENT  YP.  (I.E., COS OF ANGLE < .5)
C IF NOT, STEP SIZE WAS TOO LARGE, SO THROW AWAY Z0, AND TRY
C AGAIN WITH A SMALLER STEP.
C
          ALPHA = DOT_PRODUCT(T,YP)
          IF (ALPHA < 0.5) EXIT NEWTON
C
C MAKE  DZ  ORTHOGONAL TO TANGENT DIRECTION  YP .
C
          SIGMA = -DOT_PRODUCT(DZ,YP)/DOT_PRODUCT(T,YP)
          DZ = DZ + SIGMA*T
C
C TAKE NEWTON STEP.
C
          Z0 = Z0 + DZ
C
C CHECK FOR CONVERGENCE.
C
          XSTEP=DNRM2(NP1,DZ,1)
          IF (XSTEP .LE. RELERR*DNRM2(NP1,Z0,1)+ABSERR) EXIT CONV
C
        END DO NEWTON   ! ***** END OF NEWTON LOOP. *****
C
C DIDN'T CONVERGE OR TANGENT AT NEW POINT DID NOT MAKE
C AN ANGLE SMALLER THAN 60 DEGREES WITH  YPOLD -- 
C TRY AGAIN WITH A SMALLER H.
C      
        FAILED = .TRUE.
        HFAIL = H
        IF (H .LE. FOURU*(1.0 + S)) THEN
          IFLAG = 6
          RETURN
        ELSE
          H = .5 * H
        END IF
C
C END OF CONVERGENCE FAILURE SECTION.
C
      END DO CONV
C
C ***** CONVERGED -- MOP UP AND RETURN. *****
C
C COMPUTE TANGENT  T  AT  Z0 .
C
      CALL TANGNS(S,Z0,T,DZ,YP,A,MODE,LENQR,NFE,N,IFLAG)
      IF (IFLAG .GT. 0) RETURN
      ALPHA = DOT_PRODUCT(T,YP)
      ALPHA = ACOS(ALPHA)
C
C COMPUTE CORRECTOR DISTANCE.
C
      IF (START) THEN
        DZ = Y + H*YP
      ELSE
        DO I=1,NP1
          DZ(I)=QOFS(YOLD(I),YPOLD(I),Y(I),YP(I),HOLD,HOLD+H)
        END DO
      ENDIF
      DZ = DZ - Z0
      CORDIS = DNRM2(NP1,DZ,1)
C
C SET UP VARIABLES FOR NEXT CALL.
C
      YOLD = Y
      Y = Z0
      YPOLD = YP
      YP = T
C
C UPDATE ARCLENGTH   S = S + ||Y-YOLD||.
C
      HTEMP = HOLD
      Z0 = Z0 - YOLD
      HOLD = DNRM2(NP1,Z0,1)
      S = S+HOLD
C
C COMPUTE IDEAL ERROR FOR STEP SIZE ESTIMATION.
C
      IF (ITCNT .LE. 1) THEN
          THETA = 8.0
      ELSE IF (ITCNT .EQ. 4) THEN
          THETA = 1.0
      ELSE
          OMEGA=XSTEP/CORDIS
          IF (ITCNT .LT. 4) THEN
            LK = 4*ITCNT-7
            IF (OMEGA .GE. WRGE(LK)) THEN
              THETA = 1.0
            ELSE IF (OMEGA .GE. WRGE(LK+1)) THEN
              THETA = ACOF(LK) + ACOF(LK+1)*LOG(OMEGA)
            ELSE IF (OMEGA .GE. WRGE(LK+2)) THEN
              THETA = ACOF(LK+2) + ACOF(LK+3)*LOG(OMEGA)
            ELSE 
              THETA = 8.0
            END IF
          ELSE IF (ITCNT .GE. 7) THEN
            THETA = 0.125
          ELSE
            LK = 4*ITCNT - 16
            IF (OMEGA .GT. WRGE(LK)) THEN
              LST = 2*ITCNT - 1
              THETA = ACOF(LST) + ACOF(LST+1)*LOG(OMEGA)
            ELSE
              THETA = 0.125
            END IF
          END IF
      END IF
      IDLERR=THETA*IDLERR
C
C IDLERR SHOULD BE NO BIGGER THAN 1/2 PREVIOUS STEP.
C
      IDLERR = MIN(.5*HOLD,IDLERR)
C
C COMPUTE OPTIMAL STEP SIZE. 
C   WK = APPROXIMATE CURVATURE = 2*SIN(ALPHA/2)/HOLD  WHERE 
C        ALPHA = ARCCOS(YP*YPOLD).
C   GAMMA = EXPECTED CURVATURE FOR NEXT STEP, COMPUTED BY 
C        EXTRAPOLATING FROM CURRENT CURVATURE  WK, AND LAST 
C        CURVATURE  WKOLD.  GAMMA  IS FURTHER REQUIRED TO BE 
C        POSITIVE.
C
      IF (.NOT. START) WKOLD = WK
      WK = 2.0*ABS(SIN(.5*ALPHA))/HOLD
      IF (START) THEN
        GAMMA = WK
      ELSE 
        GAMMA = WK + HOLD/(HOLD+HTEMP)*(WK-WKOLD)
      END IF
      GAMMA = MAX(GAMMA, 0.01_R8)
      H = SQRT(2.0*IDLERR/GAMMA)
C
C     ENFORCE RESTRICTIONS ON STEP SIZE SO AS TO ENSURE STABILITY.
C        HMIN <= H <= HMAX, BMIN*HOLD <= H <= BMAX*HOLD.
C
      H = MIN(MAX(SSPAR(1),SSPAR(3)*HOLD,H),SSPAR(4)*HOLD,SSPAR(2))
      IF (FAILED) H = MIN(HFAIL,H)
      START = .FALSE.
C
C ***** END OF MOP UP SECTION. *****
C
      RETURN
      END SUBROUTINE STEPQS
      SUBROUTINE STEPS(FODE,NEQN,Y,X,H,EPS,WT,START,HOLD,K,KOLD,CRASH,
     &  PHI,P,YP,ALPHA,W,G,KSTEPS,XOLD,IVC,IV,KGI,GI, 
     &  FPWA1,FPWA2,FPWA3,FPWA4,FPWA5,IFPWA1,IFPC1,IFPC2)
C 
C   FORTRAN 90 MODIFICATION OF THE SUBROUTINE  STEP  WRITTEN BY
C   L. F. SHAMPINE AND M. K. GORDON
C 
C   ABSTRACT
C 
C   SUBROUTINE  STEPS  IS NORMALLY USED INDIRECTLY THROUGH SUBROUTINE 
C   DEABM .  BECAUSE  DEABM  SUFFICES FOR MOST PROBLEMS AND IS MUCH 
C   EASIER TO USE, USING IT SHOULD BE CONSIDERED BEFORE USING  STEPS
C   ALONE.
C 
C   SUBROUTINE STEPS INTEGRATES A SYSTEM OF  NEQN  FIRST ORDER ORDINARY 
C   DIFFERENTIAL EQUATIONS ONE STEP, NORMALLY FROM X TO X+H, USING A
C   MODIFIED DIVIDED DIFFERENCE FORM OF THE ADAMS PECE FORMULAS.  LOCAL 
C   EXTRAPOLATION IS USED TO IMPROVE ABSOLUTE STABILITY AND ACCURACY. 
C   THE CODE ADJUSTS ITS ORDER AND STEP SIZE TO CONTROL THE LOCAL ERROR 
C   PER UNIT STEP IN A GENERALIZED SENSE.  SPECIAL DEVICES ARE INCLUDED 
C   TO CONTROL ROUNDOFF ERROR AND TO DETECT WHEN THE USER IS REQUESTING 
C   TOO MUCH ACCURACY.
C 
C   THIS CODE IS COMPLETELY EXPLAINED AND DOCUMENTED IN THE TEXT, 
C   COMPUTER SOLUTION OF ORDINARY DIFFERENTIAL EQUATIONS, THE INITIAL 
C   VALUE PROBLEM  BY L. F. SHAMPINE AND M. K. GORDON.
C   FURTHER DETAILS ON USE OF THIS CODE ARE AVAILABLE IN *SOLVING 
C   ORDINARY DIFFERENTIAL EQUATIONS WITH ODE, STEP, AND INTRP*, 
C   BY L. F. SHAMPINE AND M. K. GORDON, SLA-73-1060.
C 
C 
C   THE PARAMETERS REPRESENT -- 
C      FODE -- SUBROUTINE TO EVALUATE DERIVATIVES
C      NEQN -- NUMBER OF EQUATIONS TO BE INTEGRATED 
C      Y(*) -- SOLUTION VECTOR AT X 
C      X -- INDEPENDENT VARIABLE
C      H -- APPROPRIATE STEP SIZE FOR NEXT STEP.  NORMALLY DETERMINED BY
C           CODE
C      EPS -- LOCAL ERROR TOLERANCE 
C      WT(*) -- VECTOR OF WEIGHTS FOR ERROR CRITERION 
C      START -- LOGICAL VARIABLE SET .TRUE. FOR FIRST STEP,  .FALSE.
C           OTHERWISE 
C      HOLD -- STEP SIZE USED FOR LAST SUCCESSFUL STEP
C      K -- APPROPRIATE ORDER FOR NEXT STEP (DETERMINED BY CODE)
C      KOLD -- ORDER USED FOR LAST SUCCESSFUL STEP
C      CRASH -- LOGICAL VARIABLE SET .TRUE. WHEN NO STEP CAN BE TAKEN,
C           .FALSE. OTHERWISE.
C      YP(*) -- DERIVATIVE OF SOLUTION VECTOR AT  X  AFTER SUCCESSFUL 
C           STEP
C      KSTEPS -- COUNTER ON ATTEMPTED STEPS 
C
C   THE VARIABLES X,XOLD,KOLD,KGI AND IVC AND THE ARRAYS Y,PHI,ALPHA,G, 
C   W,P,IV AND GI ARE REQUIRED FOR THE INTERPOLATION SUBROUTINE SINTRP. 
C   THE ARRAYS FPWA* AND IFPWA1 AND INTEGER CONSTANTS IFPC* ARE
C   WORKING STORAGE PASSED DIRECTLY THROUGH TO  FODE.
C 
C   INPUT TO STEPS
C 
C      FIRST CALL --
C 
C   THE USER MUST PROVIDE STORAGE IN HIS CALLING PROGRAM FOR ALL ARRAYS 
C   IN THE CALL LIST, NAMELY
C 
C     DIMENSION Y(NEQN),WT(NEQN),PHI(NEQN,16),P(NEQN),YP(NEQN), 
C    1  ALPHA(12),W(12),G(13),GI(11),IV(10),   FPWA1(NEQN),
C    2  FPWA2(NEQN-1),FPWA3(NEQN-1,NEQN),FPWA4(3*NEQN),
C    3  FPWA5(NEQN),IFPWA1(NEQN)
C                              --                --    **NOTE** 
C 
C   THE USER MUST ALSO DECLARE  START  AND  CRASH 
C   LOGICAL VARIABLES AND  FODE  AN EXTERNAL SUBROUTINE, SUPPLY THE
C   SUBROUTINE  FODE(X,Y,YP,FPWA1,FPWA2,FPWA3,FPWA4,FPWA5,IFPWA1,IFPC1,
C                 NEQN-1,IFPC2) TO EVALUATE
C      DY(I)/DX = YP(I) = F(X,Y(1),Y(2),...,Y(NEQN))
C   AND INITIALIZE ONLY THE FOLLOWING PARAMETERS. 
C      NEQN -- NUMBER OF EQUATIONS TO BE INTEGRATED 
C      Y(*) -- VECTOR OF INITIAL VALUES OF DEPENDENT VARIABLES
C      X -- INITIAL VALUE OF THE INDEPENDENT VARIABLE 
C      H -- NOMINAL STEP SIZE INDICATING DIRECTION OF INTEGRATION 
C           AND MAXIMUM SIZE OF STEP.  MUST BE VARIABLE 
C      EPS -- LOCAL ERROR TOLERANCE PER STEP.  MUST BE VARIABLE 
C      WT(*) -- VECTOR OF NON-ZERO WEIGHTS FOR ERROR CRITERION
C      START -- .TRUE.
C      KSTEPS -- SET KSTEPS TO ZERO 
C   DEFINE U TO BE THE MACHINE UNIT ROUNDOFF QUANTITY BY CALLING
C   THE INTRINSIC FUNCTION  EPSILON(1.0_R8), OR BY 
C   COMPUTING U SO THAT U IS THE SMALLEST POSITIVE NUMBER SUCH
C   THAT 1.0+U .GT. 1.0.
C 
C   STEPS  REQUIRES THAT THE L2 NORM OF THE VECTOR WITH COMPONENTS
C   LOCAL ERROR(L)/WT(L)  BE LESS THAN  EPS  FOR A SUCCESSFUL STEP.  THE
C   ARRAY  WT  ALLOWS THE USER TO SPECIFY AN ERROR TEST APPROPRIATE 
C   FOR HIS PROBLEM.  FOR EXAMPLE,
C      WT(L) = 1.0  SPECIFIES ABSOLUTE ERROR, 
C            = ABS(Y(L))  ERROR RELATIVE TO THE MOST RECENT VALUE OF THE
C                 L-TH COMPONENT OF THE SOLUTION, 
C            = ABS(YP(L))  ERROR RELATIVE TO THE MOST RECENT VALUE OF 
C                 THE L-TH COMPONENT OF THE DERIVATIVE, 
C            = MAX(WT(L),ABS(Y(L)))  ERROR RELATIVE TO THE LARGEST
C                 MAGNITUDE OF L-TH COMPONENT OBTAINED SO FAR,
C            = ABS(Y(L))*RELERR/EPS + ABSERR/EPS  SPECIFIES A MIXED 
C                 RELATIVE-ABSOLUTE TEST WHERE  RELERR  IS RELATIVE 
C                 ERROR,  ABSERR  IS ABSOLUTE ERROR AND  EPS =
C                 MAX(RELERR,ABSERR) .
C 
C      SUBSEQUENT CALLS --
C 
C   SUBROUTINE  STEPS  IS DESIGNED SO THAT ALL INFORMATION NEEDED TO
C   CONTINUE THE INTEGRATION, INCLUDING THE STEP SIZE  H  AND THE ORDER 
C   K , IS RETURNED WITH EACH STEP.  WITH THE EXCEPTION OF THE STEP 
C   SIZE, THE ERROR TOLERANCE, AND THE WEIGHTS, NONE OF THE PARAMETERS
C   SHOULD BE ALTERED.  THE ARRAY  WT  MUST BE UPDATED AFTER EACH STEP
C   TO MAINTAIN RELATIVE ERROR TESTS LIKE THOSE ABOVE.  NORMALLY THE
C   INTEGRATION IS CONTINUED JUST BEYOND THE DESIRED ENDPOINT AND THE 
C   SOLUTION INTERPOLATED THERE WITH SUBROUTINE  SINTRP .  IF IT IS 
C   IMPOSSIBLE TO INTEGRATE BEYOND THE ENDPOINT, THE STEP SIZE MAY BE 
C   REDUCED TO HIT THE ENDPOINT SINCE THE CODE WILL NOT TAKE A STEP 
C   LARGER THAN THE  H  INPUT.  CHANGING THE DIRECTION OF INTEGRATION,
C   I.E., THE SIGN OF  H , REQUIRES THE USER SET  START = .TRUE. BEFORE 
C   CALLING  STEPS  AGAIN.  THIS IS THE ONLY SITUATION IN WHICH  START
C   SHOULD BE ALTERED.
C 
C   OUTPUT FROM STEPS 
C 
C      SUCCESSFUL STEP -- 
C 
C   THE SUBROUTINE RETURNS AFTER EACH SUCCESSFUL STEP WITH  START  AND
C   CRASH  SET .FALSE. .  X  REPRESENTS THE INDEPENDENT VARIABLE
C   ADVANCED ONE STEP OF LENGTH  HOLD  FROM ITS VALUE ON INPUT AND  Y 
C   THE SOLUTION VECTOR AT THE NEW VALUE OF  X .  ALL OTHER PARAMETERS
C   REPRESENT INFORMATION CORRESPONDING TO THE NEW  X  NEEDED TO
C   CONTINUE THE INTEGRATION. 
C 
C      UNSUCCESSFUL STEP -- 
C 
C   WHEN THE ERROR TOLERANCE IS TOO SMALL FOR THE MACHINE PRECISION,
C   THE SUBROUTINE RETURNS WITHOUT TAKING A STEP AND  CRASH = .TRUE. .
C   AN APPROPRIATE STEP SIZE AND ERROR TOLERANCE FOR CONTINUING ARE 
C   ESTIMATED AND ALL OTHER INFORMATION IS RESTORED AS UPON INPUT 
C   BEFORE RETURNING.  TO CONTINUE WITH THE LARGER TOLERANCE, THE USER
C   JUST CALLS THE CODE AGAIN.  A RESTART IS NEITHER REQUIRED NOR 
C   DESIRABLE.
C***REFERENCES  SHAMPINE L.F., GORDON M.K., *SOLVING ORDINARY 
C                 DIFFERENTIAL EQUATIONS WITH ODE, STEP, AND INTRP*,
C                 SLA-73-1060, SANDIA LABORATORIES, 1973. 
C 
      USE REAL_PRECISION
      REAL (KIND=R8):: ABSH,ALPHA,BETA,EPS,ERK,ERKM1,ERKM2,
     &  ERKP1,ERR,FOURU,FPWA1,FPWA2,FPWA3,FPWA4,FPWA5,G,GI,GSTR,H,
     &  HNEW,HOLD,P,P5EPS,PHI,PSI,R,REALI,REALNS,RHO,ROUND,SIG,
     &  SUM,TAU,TEMP1,TEMP2,TEMP3,TEMP4,TEMP5,TEMP6,TWO,TWOU,V,
     &  W,WT,X,XOLD,Y,YP
      INTEGER I,IFAIL,IFPC1,IFPC2,IFPWA1,IM1,IP1,IQ,IV,IVC,
     &  J,JV,K,KGI,KM1,KM2,KNEW,KOLD,KP1,KP2,KPREV,KSTEPS,
     &  L,LIMIT1,LIMIT2,NEQN,NS,NSM2,NSP1,NSP2
      LOGICAL START,CRASH,PHASE1,NORND
C
      DIMENSION Y(:),WT(NEQN),PHI(NEQN,16),P(NEQN),YP(NEQN),PSI(12), 
     &  ALPHA(12),BETA(12),SIG(13),V(12),W(12),G(13),GI(11),IV(10), 
     &  FPWA1(NEQN),FPWA2(:),FPWA3(NEQN-1,NEQN),FPWA4(3*NEQN),
     &  FPWA5(NEQN),IFPWA1(NEQN)
      DIMENSION TWO(13),GSTR(13)
C
C   ALL LOCAL VARIABLES ARE SAVED, RATHER THAN PASSED, IN THIS
C   SPECIALIZED VERSION OF STEPS.
C
      SAVE
C
      INTERFACE
        SUBROUTINE FODE(S,Y,YP,YPOLD,A,QR,ALPHA,TZ,PIVOT,NFE,N,IFLAG)
        USE REAL_PRECISION
        REAL (KIND=R8):: S
        INTEGER:: IFLAG,N,NFE
        REAL (KIND=R8):: A(:),Y(:),YP(N+1),YPOLD(N+1)
        REAL (KIND=R8):: ALPHA(3*N+3),QR(N,N+1),TZ(N+1)
        INTEGER, DIMENSION(N+1):: PIVOT
        END SUBROUTINE FODE
      END INTERFACE
C
      DATA TWO /2.0_R8, 4.0_R8, 8.0_R8, 16.0_R8, 32.0_R8, 64.0_R8,
     &  128.0_R8, 256.0_R8, 512.0_R8, 1024.0_R8, 2048.0_R8,
     &  4096.0_R8, 8192.0_R8/
      DATA GSTR /0.500_R8, 0.0833_R8, 0.0417_R8, 0.0264_R8,
     &  0.0188_R8, 0.0143_R8, 0.0114_R8, 0.00936_R8, 0.00789_R8,
     &  0.00679_R8, 0.00592_R8, 0.00524_R8, 0.00468_R8/
C 
C 
C       ***     BEGIN BLOCK 0     *** 
C   CHECK IF STEP SIZE OR ERROR TOLERANCE IS TOO SMALL FOR MACHINE
C   PRECISION.  IF FIRST STEP, INITIALIZE PHI ARRAY AND ESTIMATE A
C   STARTING STEP SIZE. 
C                   *** 
C 
C   IF STEP SIZE IS TOO SMALL, DETERMINE AN ACCEPTABLE ONE
C 
C***FIRST EXECUTABLE STATEMENT
      TWOU = 2.0 * EPSILON(1.0_R8)
      FOURU = TWOU + TWOU
      CRASH = .TRUE.
      IF(ABS(H) .GE. FOURU*ABS(X)) GO TO 5
      H = SIGN(FOURU*ABS(X),H)
      RETURN
 5    P5EPS = 0.5*EPS 
C 
C   IF ERROR TOLERANCE IS TOO SMALL, INCREASE IT TO AN ACCEPTABLE VALUE 
C 
      ROUND = 0.0 
      DO L = 1,NEQN
        ROUND = ROUND + (Y(L)/WT(L))**2 
      END DO
      ROUND = TWOU*SQRT(ROUND)
      IF(P5EPS .GE. ROUND) GO TO 15 
      EPS = 2.0*ROUND*(1.0 + FOURU) 
      RETURN
 15   CRASH = .FALSE. 
      G(1) = 1.0
      G(2) = 0.5
      SIG(1) = 1.0
      IF(.NOT.START) GO TO 99 
C 
C   INITIALIZE.  COMPUTE APPROPRIATE STEP SIZE FOR FIRST STEP 
C 
      CALL FODE(X,Y,YP,FPWA1,FPWA2,FPWA3,FPWA4,FPWA5,IFPWA1,
     &       IFPC1,NEQN-1,IFPC2)
      IF (IFPC2 .GT. 0) RETURN
      SUM = 0.0 
      DO L = 1,NEQN
        PHI(L,1) = YP(L)
        PHI(L,2) = 0.0
        SUM = SUM + (YP(L)/WT(L))**2
      END DO
      SUM = SQRT(SUM) 
      ABSH = ABS(H) 
      IF(EPS .LT. 16.0*SUM*H*H) ABSH = 0.25*SQRT(EPS/SUM) 
      H = SIGN(MAX(ABSH,FOURU*ABS(X)),H)
C 
C*      U = D1MACH(3) 
C*      BIG = SQRT(D1MACH(2)) 
C*      CALL HSTART (F,NEQN,X,X+H,Y,YP,WT,1,U,BIG,
C*     1             PHI(1,3),PHI(1,4),PHI(1,5),PHI(1,6),RPAR,IPAR,H) 
C 
      HOLD = 0.0
      K = 1 
      KOLD = 0
      KPREV = 0 
      START = .FALSE. 
      PHASE1 = .TRUE. 
      NORND = .TRUE.
      IF(P5EPS .GT. 100.0*ROUND) GO TO 99 
      NORND = .FALSE. 
      DO L = 1,NEQN
        PHI(L,15) = 0.0 
      END DO
 99   IFAIL = 0 
C       ***     END BLOCK 0     *** 
C 
C       ***     BEGIN BLOCK 1     *** 
C   COMPUTE COEFFICIENTS OF FORMULAS FOR THIS STEP.  AVOID COMPUTING
C   THOSE QUANTITIES NOT CHANGED WHEN STEP SIZE IS NOT CHANGED. 
C                   *** 
C 
 100  KP1 = K+1 
      KP2 = K+2 
      KM1 = K-1 
      KM2 = K-2 
C 
C   NS IS THE NUMBER OF STEPS TAKEN WITH SIZE H, INCLUDING THE CURRENT
C   ONE.  WHEN K.LT.NS, NO COEFFICIENTS CHANGE
C 
      IF(H .NE. HOLD) NS = 0
      IF (NS.LE.KOLD) NS = NS+1 
      NSP1 = NS+1 
      IF (K .LT. NS) GO TO 199
C 
C   COMPUTE THOSE COMPONENTS OF ALPHA(*),BETA(*),PSI(*),SIG(*) WHICH
C   ARE CHANGED 
C 
      BETA(NS) = 1.0
      REALNS = NS 
      ALPHA(NS) = 1.0/REALNS
      TEMP1 = H*REALNS
      SIG(NSP1) = 1.0 
      IF(K .LT. NSP1) GO TO 110 
      DO I = NSP1,K 
        IM1 = I-1 
        TEMP2 = PSI(IM1)
        PSI(IM1) = TEMP1
        BETA(I) = BETA(IM1)*PSI(IM1)/TEMP2
        TEMP1 = TEMP2 + H 
        ALPHA(I) = H/TEMP1
        REALI = I 
        SIG(I+1) = REALI*ALPHA(I)*SIG(I)
      END DO
 110  PSI(K) = TEMP1
C 
C   COMPUTE COEFFICIENTS G(*) 
C 
C   INITIALIZE V(*) AND SET W(*). 
C 
      IF(NS .GT. 1) GO TO 120 
      DO IQ = 1,K 
        TEMP3 = IQ*(IQ+1) 
        V(IQ) = 1.0/TEMP3 
        W(IQ) = V(IQ) 
      END DO
      IVC = 0 
      KGI = 0 
      IF (K .EQ. 1) GO TO 140 
      KGI = 1 
      GI(1) = W(2)
      GO TO 140 
C 
C   IF ORDER WAS RAISED, UPDATE DIAGONAL PART OF V(*) 
C 
 120  IF(K .LE. KPREV) GO TO 130
      IF (IVC .EQ. 0) GO TO 122 
      JV = KP1 - IV(IVC)
      IVC = IVC - 1 
      GO TO 123 
 122  JV = 1
      TEMP4 = K*KP1 
      V(K) = 1.0/TEMP4
      W(K) = V(K) 
      IF (K .NE. 2) GO TO 123 
      KGI = 1 
      GI(1) = W(2)
 123  NSM2 = NS-2 
      IF(NSM2 .LT. JV) GO TO 130
      DO J = JV,NSM2
        I = K-J 
        V(I) = V(I) - ALPHA(J+1)*V(I+1) 
        W(I) = V(I) 
      END DO
      IF (I .NE. 2) GO TO 130 
      KGI = NS - 1
      GI(KGI) = W(2)
C 
C   UPDATE V(*) AND SET W(*)
C 
 130  LIMIT1 = KP1 - NS 
      TEMP5 = ALPHA(NS) 
      DO IQ = 1,LIMIT1
        V(IQ) = V(IQ) - TEMP5*V(IQ+1) 
        W(IQ) = V(IQ) 
      END DO
      G(NSP1) = W(1)
      IF (LIMIT1 .EQ. 1) GO TO 137
      KGI = NS
      GI(KGI) = W(2)
 137  W(LIMIT1+1) = V(LIMIT1+1) 
      IF (K .GE. KOLD) GO TO 140
      IVC = IVC + 1 
      IV(IVC) = LIMIT1 + 2
C 
C   COMPUTE THE G(*) IN THE WORK VECTOR W(*)
C 
 140  NSP2 = NS + 2 
      KPREV = K 
      IF(KP1 .LT. NSP2) GO TO 199 
      DO I = NSP2,KP1 
        LIMIT2 = KP2 - I
        TEMP6 = ALPHA(I-1)
        DO IQ = 1,LIMIT2
          W(IQ) = W(IQ) - TEMP6*W(IQ+1) 
        END DO
        G(I) = W(1) 
      END DO
 199  CONTINUE
C       ***     END BLOCK 1     *** 
C 
C       ***     BEGIN BLOCK 2     *** 
C   PREDICT A SOLUTION P(*), EVALUATE DERIVATIVES USING PREDICTED 
C   SOLUTION, ESTIMATE LOCAL ERROR AT ORDER K AND ERRORS AT ORDERS K, 
C   K-1, K-2 AS IF CONSTANT STEP SIZE WERE USED.
C                   *** 
C 
C   INCREMENT COUNTER ON ATTEMPTED STEPS
C 
      KSTEPS = KSTEPS + 1 
C 
C   CHANGE PHI TO PHI STAR
C 
      IF(K .LT. NSP1) GO TO 215 
      DO I = NSP1,K 
        TEMP1 = BETA(I) 
        DO L = 1,NEQN 
          PHI(L,I) = TEMP1*PHI(L,I) 
        END DO
      END DO
C 
C   PREDICT SOLUTION AND DIFFERENCES
C 
 215  DO L = 1,NEQN 
        PHI(L,KP2) = PHI(L,KP1) 
        PHI(L,KP1) = 0.0
        P(L) = 0.0
      END DO
      DO J = 1,K
        I = KP1 - J 
        IP1 = I+1 
        TEMP2 = G(I)
        DO L = 1,NEQN 
          P(L) = P(L) + TEMP2*PHI(L,I)
          PHI(L,I) = PHI(L,I) + PHI(L,IP1)
        END DO
      END DO
      IF(NORND) GO TO 240 
      DO L = 1,NEQN 
        TAU = H*P(L) - PHI(L,15)
        P(L) = Y(L) + TAU 
        PHI(L,16) = (P(L) - Y(L)) - TAU 
      END DO
      GO TO 250 
 240  DO L = 1,NEQN 
        P(L) = Y(L) + H*P(L)
      END DO
 250  XOLD = X
      X = X + H 
      ABSH = ABS(H) 
      CALL FODE(X,P,YP,FPWA1,FPWA2,FPWA3,FPWA4,FPWA5,IFPWA1,
     &       IFPC1,NEQN-1,IFPC2)
      IF (IFPC2 .GT. 0) RETURN
C 
C   ESTIMATE ERRORS AT ORDERS K,K-1,K-2 
C 
      ERKM2 = 0.0 
      ERKM1 = 0.0 
      ERK = 0.0 
      DO L = 1,NEQN 
        TEMP3 = 1.0/WT(L) 
        TEMP4 = YP(L) - PHI(L,1)
        IF (KM2 .GT. 0) ERKM2 = ERKM2 + ((PHI(L,KM1)+TEMP4)*TEMP3)**2
        IF (KM2 .GE. 0) ERKM1 = ERKM1 + ((PHI(L,K)+TEMP4)*TEMP3)**2 
        ERK = ERK + (TEMP4*TEMP3)**2
      END DO
      IF (KM2 .GT. 0) ERKM2 = ABSH*SIG(KM1)*GSTR(KM2)*SQRT(ERKM2) 
      IF (KM2 .GE. 0) ERKM1 = ABSH*SIG(K)*GSTR(KM1)*SQRT(ERKM1) 
      TEMP5 = ABSH*SQRT(ERK)
      ERR = TEMP5*(G(K)-G(KP1)) 
      ERK = TEMP5*SIG(KP1)*GSTR(K)
      KNEW = K
C 
C   TEST IF ORDER SHOULD BE LOWERED 
C 
      IF (KM2 > 0) THEN
        IF(MAX(ERKM1,ERKM2) .LE. ERK) KNEW = KM1
      ELSE IF (KM2 == 0) THEN
        IF(ERKM1 .LE. 0.5*ERK) KNEW = KM1 
      ENDIF
C 
C   TEST IF STEP SUCCESSFUL 
C 
      IF(ERR .LE. EPS) GO TO 400
C       ***     END BLOCK 2     *** 
C 
C       ***     BEGIN BLOCK 3     *** 
C   THE STEP IS UNSUCCESSFUL.  RESTORE  X, PHI(*,*), PSI(*) . 
C   IF THIRD CONSECUTIVE FAILURE, SET ORDER TO ONE.  IF STEP FAILS MORE 
C   THAN THREE TIMES, CONSIDER AN OPTIMAL STEP SIZE.  DOUBLE ERROR
C   TOLERANCE AND RETURN IF ESTIMATED STEP SIZE IS TOO SMALL FOR MACHINE
C   PRECISION.
C                   *** 
C 
C   RESTORE X, PHI(*,*) AND PSI(*)
C 
      PHASE1 = .FALSE.
      X = XOLD
      DO I = 1,K
        TEMP1 = 1.0/BETA(I) 
        IP1 = I+1 
        DO L = 1,NEQN 
          PHI(L,I) = TEMP1*(PHI(L,I) - PHI(L,IP1))
        END DO
      END DO
      IF(K .LT. 2) GO TO 320
      DO I = 2,K
        PSI(I-1) = PSI(I) - H 
      END DO
C 
C   ON THIRD FAILURE, SET ORDER TO ONE.  THEREAFTER, USE OPTIMAL STEP 
C   SIZE
C 
 320  IFAIL = IFAIL + 1 
      TEMP2 = 0.5 
      IF (IFAIL-3 .GT. 0) THEN
        IF(P5EPS .LT. 0.25*ERK) TEMP2 = SQRT(P5EPS/ERK) 
      ENDIF
      IF (IFAIL-3 .GE. 0) KNEW = 1
      H = TEMP2*H 
      K = KNEW
      NS = 0
      IF(ABS(H) .GE. FOURU*ABS(X)) GO TO 340
      CRASH = .TRUE.
      H = SIGN(FOURU*ABS(X),H)
      EPS = EPS + EPS 
      RETURN
 340  GO TO 100 
C       ***     END BLOCK 3     *** 
C 
C       ***     BEGIN BLOCK 4     *** 
C   THE STEP IS SUCCESSFUL.  CORRECT THE PREDICTED SOLUTION, EVALUATE 
C   THE DERIVATIVES USING THE CORRECTED SOLUTION AND UPDATE THE 
C   DIFFERENCES.  DETERMINE BEST ORDER AND STEP SIZE FOR NEXT STEP. 
C                   *** 
 400  KOLD = K
      HOLD = H
C 
C   CORRECT AND EVALUATE
C 
      TEMP1 = H*G(KP1)
      IF (NORND) THEN
        DO L = 1,NEQN 
          TEMP3 = Y(L)
          Y(L) = P(L) + TEMP1*(YP(L) - PHI(L,1))
          P(L) = TEMP3
        END DO
      ELSE
        DO L = 1,NEQN 
          TEMP3 = Y(L)
          RHO = TEMP1*(YP(L) - PHI(L,1)) - PHI(L,16)
          Y(L) = P(L) + RHO 
          PHI(L,15) = (Y(L) - P(L)) - RHO 
          P(L) = TEMP3
        END DO
      ENDIF
      CALL FODE(X,Y,YP,FPWA1,FPWA2,FPWA3,FPWA4,FPWA5,IFPWA1,
     &       IFPC1,NEQN-1,IFPC2)
      IF (IFPC2 .GT. 0) RETURN
C 
C   UPDATE DIFFERENCES FOR NEXT STEP
C 
      DO L = 1,NEQN 
        PHI(L,KP1) = YP(L) - PHI(L,1) 
        PHI(L,KP2) = PHI(L,KP1) - PHI(L,KP2)
      END DO
      DO I = 1,K
        DO L = 1,NEQN 
          PHI(L,I) = PHI(L,I) + PHI(L,KP1)
        END DO
      END DO
C 
C   ESTIMATE ERROR AT ORDER K+1 UNLESS: 
C     IN FIRST PHASE WHEN ALWAYS RAISE ORDER, 
C     ALREADY DECIDED TO LOWER ORDER, 
C     STEP SIZE NOT CONSTANT SO ESTIMATE UNRELIABLE 
C 
      ERKP1 = 0.0 
      IF(KNEW .EQ. KM1  .OR.  K .EQ. 12) PHASE1 = .FALSE. 
      IF(PHASE1) GO TO 450
      IF(KNEW .EQ. KM1) GO TO 455 
      IF(KP1 .GT. NS) GO TO 460 
      DO L = 1,NEQN 
        ERKP1 = ERKP1 + (PHI(L,KP2)/WT(L))**2 
      END DO
      ERKP1 = ABSH*GSTR(KP1)*SQRT(ERKP1)
C 
C   USING ESTIMATED ERROR AT ORDER K+1, DETERMINE APPROPRIATE ORDER 
C   FOR NEXT STEP 
C 
      IF(K .GT. 1) GO TO 445
      IF(ERKP1 .GE. 0.5*ERK) GO TO 460
      GO TO 450 
 445  IF(ERKM1 .LE. MIN(ERK,ERKP1)) GO TO 455 
      IF(ERKP1 .GE. ERK  .OR.  K .EQ. 12) GO TO 460 
C 
C   HERE ERKP1 .LT. ERK .LT. MAX(ERKM1,ERKM2) ELSE ORDER WOULD HAVE 
C   BEEN LOWERED IN BLOCK 2.  THUS ORDER IS TO BE RAISED
C 
C   RAISE ORDER 
C 
 450  K = KP1 
      ERK = ERKP1 
      GO TO 460 
C 
C   LOWER ORDER 
C 
 455  K = KM1 
      ERK = ERKM1 
C 
C   WITH NEW ORDER DETERMINE APPROPRIATE STEP SIZE FOR NEXT STEP
C 
 460  HNEW = H + H
      IF(PHASE1) GO TO 465
      IF(P5EPS .GE. ERK*TWO(K+1)) GO TO 465 
      HNEW = H
      IF(P5EPS .GE. ERK) GO TO 465
      TEMP2 = K+1 
      R = (P5EPS/ERK)**(1.0/TEMP2)
      HNEW = ABSH*MAX(0.5_R8,MIN(0.9_R8,R)) 
      HNEW = SIGN(MAX(HNEW,FOURU*ABS(X)),H) 
 465  H = HNEW
      RETURN
C       ***     END BLOCK 4     *** 
      END SUBROUTINE STEPS
      SUBROUTINE STRPTP(N,ICOUNT,IDEG,R,X)
C
C COMPUTES INITIAL POINTS FOR PATHS.
C
C ON INPUT:
C
C N  IS THE NUMBER OF (COMPLEX) VARIABLES.
C
C ICOUNT  IS A COUNTER USED TO INCREMENT EACH
C   VARIABLE AROUND THE UNIT CIRCLE SO THAT EVERY
C   COMBINATION OF START VALUES IS CHOSEN.  ICOUNT  IS
C   INITIALIZED IN  POLSYS1H .
C
C IDEG(J)  IS THE DEGREE OF THE J-TH EQUATION.
C
C R(I,J)  IS A (COMPLEX) ARRAY GENERATED BY SUBROUTINE  INITP.
C   R(1,J), AND R(2,J) ARE THE REAL AND IMAGINARY PARTS, RESPECTIVELY.
C
C ON OUTPUT:
C
C X(1:2*N)  IS INITIALIZED TO THE START VALUES FOR THE CURRENT PATH,
C   WITH X(2*J-1) AND X(2*J) THE REAL AND IMAGINARY PARTS OF THE 
C   J-TH VARIABLE, RESPECTIVELY.
C
C FUNCTIONS USED:  ATAN, COS, SIN.
C
      USE REAL_PRECISION
C DECLARATION OF INPUT AND OUTPUT:
      INTEGER:: N,ICOUNT(N),IDEG(N)
      REAL (KIND=R8):: R(2,N),X(2*N)
C
C DECLARATION OF LOCAL VARIABLES:
      INTEGER:: J
      REAL (KIND=R8):: ANGLE,TWOPI
      COMPLEX (KIND=R8):: XXXX
C
      DO J=1,N
        IF (ICOUNT(J) .GE. IDEG(J)) THEN
          ICOUNT(J)=1
        ELSE
          ICOUNT(J)=ICOUNT(J)+1
          EXIT
        END IF
      END DO
      TWOPI = 8.0_R8*ATAN(1.0_R8)
      DO J=1,N
        ANGLE = ( TWOPI/IDEG(J) )*ICOUNT(J)
        XXXX = CMPLX(COS(ANGLE),SIN(ANGLE),KIND=R8)*
     &        CMPLX(R(1,J),R(2,J),KIND=R8)
        X(2*J-1) = REAL(XXXX)
        X(2*J) = AIMAG(XXXX)
      END DO
      RETURN
      END SUBROUTINE STRPTP
      SUBROUTINE TANGNF(RHOLEN,Y,YP,YPOLD,A,QR,ALPHA,TZ,PIVOT,
     &   NFE,N,IFLAG)
C
C THIS SUBROUTINE BUILDS THE JACOBIAN MATRIX OF THE HOMOTOPY MAP,
C COMPUTES A QR DECOMPOSITION OF THAT MATRIX, AND THEN CALCULATES THE
C (UNIT) TANGENT VECTOR AND THE NEWTON STEP.
C
C THE FOLLOWING INTERFACE BLOCK SHOULD BE INCLUDED IN THE CALLING
C PROGRAM:
C
C     INTERFACE
C       SUBROUTINE TANGNF(RHOLEN,Y,YP,YPOLD,A,QR,ALPHA,TZ,PIVOT,
C    &    NFE,N,IFLAG)
C       USE REAL_PRECISION
C       REAL (KIND=R8):: RHOLEN
C       INTEGER:: IFLAG,N,NFE
C       REAL (KIND=R8):: A(:),Y(:),YP(N+1),YPOLD(N+1)
C       REAL (KIND=R8):: ALPHA(3*N+3),QR(N,N+2),TZ(N+1)
C       INTEGER:: PIVOT(N+1)
C       END SUBROUTINE TANGNF
C     END INTERFACE
C
C
C ON INPUT:
C
C RHOLEN < 0 IF THE NORM OF THE HOMOTOPY MAP EVALUATED AT
C    (A, LAMBDA, X) IS TO BE COMPUTED.  IF  RHOLEN >= 0  THE NORM IS NOT
C    COMPUTED AND  RHOLEN  IS NOT CHANGED.
C
C Y(1:N+1) = CURRENT POINT (LAMBDA(S), X(S)).
C
C YPOLD(1:N+1) = UNIT TANGENT VECTOR AT PREVIOUS POINT ON THE ZERO
C    CURVE OF THE HOMOTOPY MAP.
C
C A(:) = PARAMETER VECTOR IN THE HOMOTOPY MAP.
C
C QR(1:N,1:N+2), ALPHA(1:3*N+3), TZ(1:N+1), PIVOT(1:N+1)  ARE WORK
C    ARRAYS USED FOR THE QR FACTORIZATION.
C
C NFE = NUMBER OF JACOBIAN MATRIX EVALUATIONS = NUMBER OF HOMOTOPY
C    FUNCTION EVALUATIONS.
C
C N = DIMENSION OF X.
C
C IFLAG = -2, -1, OR 0, INDICATING THE PROBLEM TYPE.
C
C
C ON OUTPUT:
C
C RHOLEN = ||RHO(A, LAMBDA(S), X(S)|| IF  RHOLEN < 0  ON INPUT.
C    OTHERWISE  RHOLEN  IS UNCHANGED.
C
C Y, YPOLD, A, N  ARE UNCHANGED.
C
C YP(1:N+1) = DY/DS = UNIT TANGENT VECTOR TO INTEGRAL CURVE OF
C    D(HOMOTOPY MAP)/DS = 0  AT  Y(S) = (LAMBDA(S), X(S)) .
C
C TZ = THE NEWTON STEP = -(PSEUDO INVERSE OF  (D RHO(A,Y(S))/D LAMBDA ,
C    D RHO(A,Y(S))/DX)) * RHO(A,Y(S)) .
C
C NFE  HAS BEEN INCRMENTED BY 1.
C
C IFLAG  IS UNCHANGED, UNLESS THE QR FACTORIZATION DETECTS A RANK < N,
C    IN WHICH CASE THE TANGENT AND NEWTON STEP VECTORS ARE NOT COMPUTED
C    AND  TANGNF  RETURNS WITH  IFLAG = 4 .
C
C
C CALLS  DGEQPF , DNRM2 , DORMQR , F (OR  RHO ), FJAC (OR  RHOJAC ).
C
      USE HOMOTOPY
      USE REAL_PRECISION
      REAL (KIND=R8):: LAMBDA,RHOLEN,SIGMA,YPNORM
      INTEGER:: I,IFLAG,J,K,KP1,N,NFE,NP1,NP2
      INTERFACE
        FUNCTION DNRM2(N,X,STRIDE)
        USE REAL_PRECISION
        INTEGER:: N,STRIDE
        REAL (KIND=R8):: DNRM2,X(N)
        END FUNCTION DNRM2
      END INTERFACE
C
C *****  ARRAY DECLARATIONS.  *****
C
      REAL (KIND=R8):: A(:),Y(:),YP(N+1),YPOLD(N+1)
C
C ARRAYS AND FLAG FOR COMPUTING THE JACOBIAN MATRIX AND ITS KERNEL.
      REAL (KIND=R8):: ALPHA(3*N+3),QR(N,N+2),TZ(N+1)
      INTEGER:: PIVOT(N+1)
C
C *****  END OF DIMENSIONAL INFORMATION.  *****
C
C
      LAMBDA=Y(1)
      NP1=N+1
      NP2=N+2
      NFE=NFE+1
C NFE CONTAINS THE NUMBER OF JACOBIAN EVALUATIONS.
C   *   *   *   *   *   *   *   *   *   *   *   *   *   *   *   *   *
C
C COMPUTE THE JACOBIAN MATRIX, STORE IT AND HOMOTOPY MAP IN QR.
C
      IF (IFLAG .EQ. -2) THEN
C
C  QR = ( D RHO(A,LAMBDA,X)/D LAMBDA , D RHO(A,LAMBDA,X)/DX ,
C                                              RHO(A,LAMBDA,X) )  .
C
        DO K=1,NP1
          CALL RHOJAC(A,LAMBDA,Y(2:NP1),QR(:,K),K)
        END DO
        CALL RHO(A,LAMBDA,Y(2:NP1),QR(:,NP2))
      ELSE
        CALL F(Y(2:NP1),TZ(1:N))
        IF (IFLAG .EQ. 0) THEN
C
C      QR = ( A - F(X), I - LAMBDA*DF(X) ,
C                                 X - A + LAMBDA*(A - F(X)) )  .
C
          QR(:,1)=A - TZ(1:N)
          QR(:,NP2)=Y(2:NP1) - A + LAMBDA*QR(:,1)
          DO K=1,N
            CALL FJAC(Y(2:NP1),TZ(1:N),K)
            KP1=K+1
            QR(:,KP1)=-LAMBDA*TZ(1:N)
            QR(K,KP1)=1.0+QR(K,KP1)
          END DO
        ELSE
C
C   QR = ( F(X) - X + A, LAMBDA*DF(X) + (1 - LAMBDA)*I ,
C                                  X - A + LAMBDA*(F(X) - X + A) )  .
C
          QR(:,1)=TZ(1:N) - Y(2:NP1) + A
          QR(:,NP2)=Y(2:NP1) - A + LAMBDA*QR(:,1)
          DO K=1,N
            CALL FJAC(Y(2:NP1),TZ(1:N),K)
            KP1=K+1
            QR(:,KP1)=LAMBDA*TZ(1:N)
            QR(K,KP1)=1.0-LAMBDA+QR(K,KP1)
          END DO
        ENDIF
      ENDIF
C
C   *   *   *   *   *   *   *   *   *   *   *   *   *   *   *   *   *
C COMPUTE THE NORM OF THE HOMOTOPY MAP IF IT WAS REQUESTED.
      IF (RHOLEN .LT. 0.0) RHOLEN=DNRM2(N,QR(:,NP2),1)
C
C REDUCE THE JACOBIAN MATRIX TO UPPER TRIANGULAR FORM.
C
      PIVOT = 0
C
      CALL DGEQPF(N,NP1,QR,N,PIVOT,YP,ALPHA,K)
C
      IF (ABS(QR(N,N)) .LE. ABS(QR(1,1))*EPSILON(1.0_R8)) THEN 
        IFLAG=4
        RETURN
      ENDIF
C
      CALL DORMQR('L','T',N,1,N,QR,N,YP,QR(:,NP2),N,
     &           ALPHA,3*N+3,K)
C
      DO I=1,N
        ALPHA(I)=QR(I,I)
      END DO
C
C COMPUTE KERNEL OF JACOBIAN, WHICH SPECIFIES YP=DY/DS.
      TZ(NP1)=1.0
      DO I=N,1,-1
        J=I+1
        TZ(I)=-DOT_PRODUCT(QR(I,J:NP1),TZ(J:NP1))/ALPHA(I)
      END DO
      YPNORM=DNRM2(NP1,TZ,1)
      YP(PIVOT)=TZ/YPNORM
      IF (DOT_PRODUCT(YP,YPOLD) .LT. 0.0) YP = -YP
C YP  IS THE UNIT TANGENT VECTOR IN THE CORRECT DIRECTION.
C
C COMPUTE THE MINIMUM NORM SOLUTION OF [D RHO(Y(S))] V = -RHO(Y(S)).
C V IS GIVEN BY  P - (P,Q)Q  , WHERE P IS ANY SOLUTION OF
C [D RHO] V = -RHO  AND Q IS A UNIT VECTOR IN THE KERNEL OF [D RHO].
C
      ALPHA(NP1)=1.0
      DO I=N,1,-1
        J=I+1
        ALPHA(I)=-(DOT_PRODUCT(QR(I,J:NP1),ALPHA(J:NP1)) + QR(I,NP2))
     &          /ALPHA(I)
      END DO
      TZ(PIVOT)=ALPHA(1:NP1)
C TZ NOW CONTAINS A PARTICULAR SOLUTION P, AND YP CONTAINS A VECTOR Q
C IN THE KERNEL(THE TANGENT).
      SIGMA=DOT_PRODUCT(TZ,YP)
      TZ=TZ - SIGMA*YP
C TZ IS THE NEWTON STEP FROM THE CURRENT POINT Y(S) = (LAMBDA(S), X(S)).
      RETURN
      END SUBROUTINE TANGNF
      SUBROUTINE TANGNS(RHOLEN,Y,YP,TZ,YPOLD,A,MODE,LENQR,NFE,N,IFLAG)
C
C THIS SUBROUTINE BUILDS THE JACOBIAN MATRIX OF THE HOMOTOPY MAP,
C AND THEN CALCULATES THE (UNIT) TANGENT VECTOR AND THE NEWTON STEP
C USING A PRECONDITIONED CONJUGATE GRADIENT ALGORITHM.
C
C THE CALLING PROGRAM MUST CONTAIN THE FOLLOWING INTERFACE BLOCK:
C
C     INTERFACE
C       SUBROUTINE TANGNS(RHOLEN,Y,YP,TZ,YPOLD,A,MODE,LENQR,
C    &    NFE,N,IFLAG)
C       USE HOMOTOPY, QR => QRSPARSE
C       USE REAL_PRECISION
C       REAL (KIND=R8), INTENT(IN), DIMENSION(:):: A,Y,YPOLD
C       REAL (KIND=R8), INTENT(IN OUT):: RHOLEN
C       REAL (KIND=R8), INTENT(OUT), DIMENSION(:):: TZ,YP
C       INTEGER, INTENT(IN):: LENQR,MODE,N
C       INTEGER, INTENT(IN OUT):: IFLAG,NFE
C       END SUBROUTINE TANGNS
C     END INTERFACE
C
C
C ON INPUT:
C
C RHOLEN < 0 IF THE NORM OF THE HOMOTOPY MAP EVALUATED AT
C    (A, X, LAMBDA) IS TO BE COMPUTED.  IF  RHOLEN >= 0  THE NORM IS NOT
C    COMPUTED AND  RHOLEN  IS NOT CHANGED.
C
C Y(1:N+1) = CURRENT POINT (X(S), LAMBDA(S)).
C
C YPOLD(1:N+1) = UNIT TANGENT VECTOR AT PREVIOUS POINT ON THE ZERO
C    CURVE OF THE HOMOTOPY MAP.
C
C A(:) = PARAMETER VECTOR IN THE HOMOTOPY MAP.
C
C MODE = 1 IF THE JACOBIAN MATRIX IS SYMMETRIC AND STORED IN A PACKED
C          SKYLINE FORMAT;
C      = 2 IF THE JACOBIAN MATRIX IS STORED IN A SPARSE ROW FORMAT.
C
C LENQR  IS THE NUMBER OF NONZERO ENTRIES IN THE SPARSE JACOBIAN
C    MATRICES, USED TO DETERMINE THE SPARSE MATRIX DATA STRUCTURES.
C
C NFE = NUMBER OF JACOBIAN MATRIX EVALUATIONS = NUMBER OF HOMOTOPY
C    FUNCTION EVALUATIONS.
C
C N = DIMENSION OF X.
C
C IFLAG = -2, -1, OR 0, INDICATING THE PROBLEM TYPE.
C
C
C ON OUTPUT:
C
C RHOLEN = ||RHO(A, X(S), LAMBDA(S)|| IF  RHOLEN < 0  ON INPUT.
C    OTHERWISE  RHOLEN  IS UNCHANGED.
C
C Y, YPOLD, A, N  ARE UNCHANGED.
C
C YP(1:N+1) = DY/DS = UNIT TANGENT VECTOR TO INTEGRAL CURVE OF
C    D(HOMOTOPY MAP)/DS = 0  AT  Y(S) = (X(S), LAMBDA(S)) .
C
C TZ(1:N+1) = THE NEWTON STEP = -(PSEUDO INVERSE OF  (D RHO(A,Y(S))/DX ,
C    D RHO(A,Y(S))/D LAMBDA)) * RHO(A,Y(S)) .  THE NEWTON STEP IS
C    CALCULATED ONLY IF RHOLEN < 0 ON INPUT.
C
C NFE  HAS BEEN INCRMENTED BY 1.
C
C IFLAG  IS UNCHANGED, UNLESS THE PRECONDITIONED CONJUGATE GRADIENT
C    ITERATION FAILS TO CONVERGE, IN WHICH CASE THE TANGENT AND NEWTON 
C    STEP VECTORS ARE NOT COMPUTED AND  TANGNS  RETURNS WITH  IFLAG = 4 .
C
C
C CALLS  F (OR  RHO ), FJACS (OR  RHOJS ), PCGDS , GMRILUDS , AND THE
C    BLAS ROUTINE  DNRM2 .
C
        USE HOMOTOPY, QR => QRSPARSE
        USE REAL_PRECISION
        REAL (KIND=R8), INTENT(IN), DIMENSION(:):: A,Y,YPOLD
        REAL (KIND=R8), INTENT(IN OUT):: RHOLEN
        REAL (KIND=R8), INTENT(OUT), DIMENSION(:):: TZ,YP
        INTEGER, INTENT(IN):: LENQR,MODE,N
        INTEGER, INTENT(IN OUT):: IFLAG,NFE
C
C ***** LOCAL VARIABLES AND AUTOMATIC WORK ARRAYS. *****
C
      REAL (KIND=R8):: LAMBDA,RHOVEC(N),SIGMA,YPNORM
      INTEGER:: J,NP1,JPOS
C
      INTERFACE
        FUNCTION DNRM2(N,X,STRIDE)
          USE REAL_PRECISION
          INTEGER:: N,STRIDE
          REAL (KIND=R8):: DNRM2,X(N)
        END FUNCTION DNRM2
        SUBROUTINE PCGDS(N,LENQR,IFLAG,YP,RHS)
          USE REAL_PRECISION
          INTEGER, INTENT(IN):: LENQR,N
          INTEGER, INTENT(IN OUT):: IFLAG
          REAL (KIND=R8), INTENT(IN OUT):: YP(N+1)
          REAL (KIND=R8), OPTIONAL, INTENT(IN):: RHS(N)
        END SUBROUTINE PCGDS
        SUBROUTINE GMRILUDS(N,LENQR,IFLAG,YP,RHS)
          USE REAL_PRECISION
          INTEGER, INTENT(IN):: LENQR,N
          INTEGER, INTENT(IN OUT):: IFLAG
          REAL (KIND=R8), INTENT(IN OUT):: YP(N+1)
          REAL (KIND=R8), OPTIONAL, INTENT(IN):: RHS(N)
        END SUBROUTINE GMRILUDS
      END INTERFACE
C
C *****  END OF SPECIFICATION INFORMATION.  *****
C
      NP1=N+1
      NFE=NFE+1
C NFE CONTAINS THE NUMBER OF JACOBIAN EVALUATIONS.
      LAMBDA=Y(NP1)
      ROWPOS(NP1)=LENQR+1
C   *   *   *   *   *   *   *   *   *   *   *   *   *   *   *   *   *
C MODE = 1 STORAGE FORMAT.
C   *   *   *   *   *   *   *   *   *   *   *   *   *   *   *   *   *
C
      IF (MODE .EQ. 1) THEN
C COMPUTE THE JACOBIAN MATRIX, STORE IT IN  [QR | -PP] .
C
      IF (IFLAG .EQ. -2) THEN
C
C  [QR | -PP] = [ D RHO(A,X,LAMBDA)/DX | D RHO(A,X,LAMBDA)/D LAMBDA ]  .
C  RHOVEC = RHO(A,X,LAMBDA) .
C
C  PP = - (D RHO(A,X,LAMBDA)/D LAMBDA) .
        CALL RHOJS(A,LAMBDA,Y(1:N))
        IF (RHOLEN < 0) CALL RHO(A,LAMBDA,Y(1:N),RHOVEC)
C
      ELSE
        CALL F(Y(1:N),PP)
        IF (IFLAG .EQ. 0) THEN
C
C      [QR | -PP] = [ I - LAMBDA*DF(X) | A - F(X) ]  .
C      RHOVEC = X - A + LAMBDA*(A - F(X)) .
C
          PP = PP - A(1:N)
          IF (RHOLEN < 0) RHOVEC = Y(1:N) - A(1:N) - LAMBDA*PP
          CALL FJACS(Y(1:N))
          QR = (-LAMBDA)*QR
          QR(ROWPOS(1:N)) = QR(ROWPOS(1:N)) + 1.0
        ELSE
C
C   [QR | -PP] = [ LAMBDA*DF(X) + (1 - LAMBDA)*I | F(X) - X + A ] .
C   RHOVEC = X - A + LAMBDA*(F(X) - X + A) .
C
          PP = Y(1:N) - A(1:N) - PP
          IF (RHOLEN < 0) RHOVEC = Y(1:N) - A(1:N) - LAMBDA*PP
          CALL FJACS(Y(1:N))
          QR = LAMBDA*QR
          QR(ROWPOS(1:N)) = QR(ROWPOS(1:N)) + 1.0 - LAMBDA
        ENDIF
      ENDIF
      ELSE
C   *   *   *   *   *   *   *   *   *   *   *   *   *   *   *   *   *
C MODE = 2 STORAGE FORMAT.
C   *   *   *   *   *   *   *   *   *   *   *   *   *   *   *   *   *
C
      IF (IFLAG .EQ. -2) THEN
C
C  [QR] = [ D RHO(A,X,LAMBDA)/DX , D RHO(A,X,LAMBDA)/D LAMBDA ]  .
C  RHOVEC = RHO(A,X,LAMBDA) .
C
        CALL RHOJS(A,LAMBDA,Y(1:N))
        IF (RHOLEN < 0) CALL RHO(A,LAMBDA,Y(1:N),RHOVEC)
C
      ELSE
        CALL F(Y(1:N),PP)
        IF (IFLAG .EQ. 0) THEN
C
C      [QR | -PP] = [ I - LAMBDA*DF(X) | A - F(X) ]  .
C      RHOVEC = X - A + LAMBDA*(A - F(X)) .
C
          PP = PP - A(1:N)
          IF (RHOLEN < 0) RHOVEC = Y(1:N) - A(1:N) - LAMBDA*PP
          CALL FJACS(Y(1:N))
          QR = (-LAMBDA)*QR
C FIND INDEX JPOS OF DIAGONAL ELEMENT IN JTH ROW OF QR.
          DO J=1,N
            JPOS=ROWPOS(J)
            DO
              IF (COLPOS(JPOS) .EQ. J) EXIT
              JPOS=JPOS+1
              IF (JPOS < ROWPOS(J+1)) CYCLE
              IFLAG=4
              RETURN
            END DO
            QR(JPOS) = QR(JPOS) + 1.0
          END DO
        ELSE
C
C   [QR | -PP] = [ LAMBDA*DF(X) + (1 - LAMBDA)*I | F(X) - X + A ] .
C   RHOVEC = X - A + LAMBDA*(F(X) - X + A) .
C
          PP = Y(1:N) - A(1:N) - PP
          IF (RHOLEN < 0) RHOVEC = Y(1:N) - A(1:N) - LAMBDA*PP
          CALL FJACS(Y(1:N))
          QR = LAMBDA*QR
C FIND INDEX JPOS OF DIAGONAL ELEMENT IN JTH ROW OF QR.
          DO J=1,N
            JPOS=ROWPOS(J)
            DO
              IF (COLPOS(JPOS) .EQ. J) EXIT
              JPOS=JPOS+1
              IF (JPOS < ROWPOS(J+1)) CYCLE
              IFLAG=4
              RETURN
            END DO
            QR(JPOS) = QR(JPOS) + 1.0 - LAMBDA
          END DO
        ENDIF
      ENDIF
      ENDIF
C
C   *   *   *   *   *   *   *   *   *   *   *   *   *   *   *   *   *
      YP=YPOLD
C COMPUTE KERNEL OF JACOBIAN, WHICH SPECIFIES YP=DY/DS, USING A
C PRECONDITIONED CONJUGATE GRADIENT ALGORITHM.
      SELECT CASE (MODE)
        CASE (1)
        CALL PCGDS(N,LENQR,IFLAG,YP)
        CASE (2)
        CALL GMRILUDS(N,LENQR,IFLAG,YP)
      END SELECT
      IF (IFLAG .GT. 0) RETURN
C
C NORMALIZE TANGENT VECTOR YP.
      YPNORM=DNRM2(NP1,YP,1)
      YP = (1.0/YPNORM)*YP
C
C CHOOSE UNIT TANGENT VECTOR DIRECTION TO MAINTAIN CONTINUITY.
      IF (DOT_PRODUCT(YP,YPOLD) .LT. 0.0) YP = -YP
C
C COMPUTE THE NORM OF THE HOMOTOPY MAP AND THE NEWTON STEP
C IF IT WAS REQUESTED.
      IF (RHOLEN < 0.0) THEN
        RHOLEN=DNRM2(N,RHOVEC,1)
C COMPUTE THE MINIMUM NORM SOLUTION OF [D RHO(Y(S))] V = -RHO(Y(S)).
C V IS GIVEN BY  P - (P,Q)Q  , WHERE P IS ANY SOLUTION OF
C [D RHO] V = -RHO  AND Q IS A UNIT VECTOR IN THE KERNEL OF [D RHO].
        TZ=YPOLD
        SELECT CASE (MODE)
          CASE (1)
          CALL PCGDS(N,LENQR,IFLAG,TZ,RHS=RHOVEC)
          CASE (2)
          CALL GMRILUDS(N,LENQR,IFLAG,TZ,RHS=RHOVEC)
        END SELECT
        IF (IFLAG .GT. 0) RETURN
C TZ NOW CONTAINS A PARTICULAR SOLUTION P, AND YP CONTAINS A VECTOR Q
C IN THE KERNEL(THE TANGENT).
        SIGMA=DOT_PRODUCT(TZ(1:NP1),YP(1:NP1))
        TZ = TZ - SIGMA*YP
C TZ IS THE NEWTON STEP FROM THE CURRENT POINT Y(S) = (X(S), LAMBDA(S)).
      END IF
C
      RETURN
      END SUBROUTINE TANGNS
        SUBROUTINE TANGQF(Y,YP,YPOLD,A,Q,R,W,S,T,N,IFLAG,NFE)
C
C SUBROUTINE  TANGQF  COMPUTES THE UNIT TANGENT VECTOR  YP  TO THE
C ZERO CURVE OF THE HOMOTOPY MAP AT  Y  BY GENERATING THE AUGMENTED 
C JACOBIAN MATRIX  
C
C           --           --
C           |  D(RHO(Y))  |      
C     AUG = |        T    |,   WHERE RHO IS THE HOMOTOPY MAP,
C           |   YPOLD     | 
C           --           --
C
C SOLVING THE SYSTEM
C                                T
C         AUG*YPT = (0,0,...,0,1)    FOR YPT,
C
C AND FINALLY COMPUTING  YP = YPT/||YPT||.
C
C IN ADDITION, THE MATRIX AUG IS UPDATED SO THAT THE LAST ROW IS
C YP  INSTEAD OF  YPOLD  ON RETURN.
C
C THE FOLLOWING INTERFACE BLOCK SHOULD BE INCLUDED IN THE CALLING
C PROGRAM:
C
C     INTERFACE
C       SUBROUTINE TANGQF(Y,YP,YPOLD,A,Q,R,W,S,T,N,IFLAG,NFE)
C       USE HOMOTOPY
C       USE REAL_PRECISION
C       INTEGER:: N, IFLAG, NFE
C       REAL (KIND=R8):: A(:), Q(N+1,N+1), R((N+1)*(N+2)/2),
C    &    S(N+1), T(N+1), W(N+1), Y(:), YP(N+1), YPOLD(N+1)
C       END SUBROUTINE TANGQF
C     END INTERFACE
C
C
C ON INPUT:
C 
C Y(1:N+1) = CURRENT POINT (LAMBDA(S), X(S)).
C
C YP(1:N+1)  IS UNDEFINED ON INPUT.
C 
C YPOLD(1:N+1) = UNIT TANGENT VECTOR AT THE PREVIOUS POINT ON THE 
C    ZERO CURVE OF THE HOMOTOPY MAP.
C
C A(:)  IS THE PARAMETER VECTOR IN THE HOMOTOPY MAP.
C
C W(1:N+1), S(1:N+1), T(1:N+1)  ARE WORK ARRAYS.
C
C N  IS THE DIMENSION OF X, WHERE  Y=(LAMBDA(S),X(S)).
C
C IFLAG  IS -2, -1, OR 0, INDICATING THE PROBLEM TYPE.
C
C NFE  IS THE NUMBER OF JACOBIAN EVALUATIONS.
C
C
C ON OUTPUT:
C
C Y, YPOLD, A, N  ARE UNCHANGED.
C
C YP(1:N+1)  CONTAINS THE NEW UNIT TANGENT VECTOR TO THE ZERO
C    CURVE OF THE HOMOTOPY MAP AT  Y(S) = (LAMBDA(S), X(S)).
C
C Q(1:N+1,1:N+1)  CONTAINS  Q  OF THE QR FACTORIZATION OF
C    THE JACOBIAN MATRIX OF RHO EVALUATED AT  Y  AUGMENTED BY  
C    YP TRANSPOSE.
C
C R(1:(N+1)*(N+2)/2)  CONTAINS THE UPPER TRIANGLE (STORED BY COLUMNS)
C    OF THE  R  PART OF THE QR FACTORIZATION OF THE AUGMENTED JACOBIAN
C    MATRIX.
C
C IFLAG  = -2, -1, OR 0, (UNCHANGED) ON A NORMAL RETURN.
C        = 4 IF THE AUGMENTED JACOBIAN MATRIX HAS RANK LESS THAN N+1.
C 
C NFE  HAS BEEN INCREMENTED BY 1.
C
C
C CALLS  DGEQRF, DNRM2, DORGQR, DTPSV, F (OR RHO IF IFLAG = -2),
C FJAC (OR RHOJAC, IF IFLAG = -2),
C R1UPQF (WHICH IS AN ENTRY POINT OF UPQRQF).
C        
C ***** DECLARATIONS *****
      USE HOMOTOPY
      USE REAL_PRECISION
C
C FUNCTION DECLARATIONS
C
      INTERFACE
        FUNCTION DNRM2(N,X,STRIDE)
        USE REAL_PRECISION
        INTEGER:: N,STRIDE
        REAL (KIND=R8):: DNRM2,X(N)
        END FUNCTION DNRM2
      END INTERFACE
C
C LOCAL VARIABLES
C
      REAL (KIND=R8):: LAMBDA, YPNRM
      INTEGER:: I, J, JP1, NP1
C
C SCALAR ARGUMENTS
C
      INTEGER:: N, IFLAG, NFE
C
C ARRAY DECLARATIONS 
C
      REAL (KIND=R8):: A(:), Q(N+1,N+1), R((N+1)*(N+2)/2),
     &  S(N+1), T(N+1), W(N+1), Y(:), YP(N+1), YPOLD(N+1)
C
C ***** END OF DECLARATIONS *****
C
C ***** FIRST EXECUTABLE STATEMENT *****
C
      NFE = NFE + 1
      NP1 = N + 1
      LAMBDA = Y(1)
C        
C ***** DEFINE THE AUGMENTED JACOBIAN MATRIX *****
C
C Q = AUG.
C
      IF (IFLAG .EQ. -2) THEN
C
C CURVE TRACKING PROBLEM:
C         D(RHO) = (D RHO(A,LAMBDA,X)/D LAMBDA, D RHO(A,LAMBDA,X)/DX).
C
        DO J = 1,NP1
          CALL RHOJAC(A,LAMBDA,Y(2:NP1),Q(1:N,J),J)
        END DO
      ELSE IF (IFLAG .EQ. -1) THEN
C
C ZERO FINDING PROBLEM:
C         D(RHO) = (F(X) - X + A, LAMBDA*DF(X) + (1-LAMBDA)*I)
C
        CALL F(Y(2:NP1),Q(1:N,1))
        Q(1:N,1) = A(1:N) - Y(2:NP1) + Q(1:N,1)
        DO J= 1,N
          JP1 = J+1
          CALL FJAC(Y(2:NP1),Q(1:N,JP1),J)
          Q(1:N,JP1) = LAMBDA*Q(1:N,JP1)
          Q(J,JP1) = 1.0 - LAMBDA + Q(J,JP1)
        END DO
      ELSE 
C
C FIXED POINT PROBLEM:
C         D(RHO) = (A - F(X), I - LAMBDA*DF(X)).
C
        CALL F(Y(2:NP1),Q(1:N,1))
        Q(1:N,1) = A(1:N) - Q(1:N,1)
        DO J=1,N
          JP1 = J+1
          CALL FJAC(Y(2:NP1),Q(1:N,JP1),J)
          Q(1:N,JP1) = -LAMBDA*Q(1:N,JP1)
          Q(J,JP1) = 1.0 + Q(J,JP1)
        END DO
      END IF
C
C DEFINE LAST ROW OF Q = YPOLD.
C
      Q(NP1,:) = YPOLD
C
C ***** END OF DEFINITION OF AUGMENTED JACOBIAN MATRIX *****
C
C                                          T
C ***** SOLVE SYSTEM  AUG*YPT = (0,...,0,1)  *****
C
C FACTOR MATRIX.
C           
      CALL DGEQRF(NP1,NP1,Q,NP1,T,W,NP1,I)
C
C PACK UPPER TRIANGLE INTO ARRAY R .
C
      DO I=1,NP1
        R((I*(I-1))/2 + 1:(I*(I-1))/2 + I) = Q(1:I,I)
      END DO
C
C IF MATRIX IS SINGULAR, THEN RETURN WITH IFLAG = 4,
C ELSE SOLVE SYSTEM  R*YP = QT*(0,...,0,1)  FOR YP.
C
C
C CHECK FOR SINGULARITY.
C
      J = 1
      DO I = 1, N
        IF( R( J+I-1 ).EQ. 0.0 ) THEN
          IFLAG = 4
          RETURN
        END IF
        J = J + I
      END DO
C
C EXPAND HOUSEHOLDER REFLECTIONS INTO FULL MATRIX Q . 
C
      CALL DORGQR(NP1, NP1, N, Q, NP1, T, W, NP1, I)
C
      YP = Q(NP1,:)
      CALL DTPSV('U', 'N', 'N', NP1, R, YP, 1)
C
C COMPUTE UNIT VECTOR.
C
      YPNRM = 1.0/DNRM2(NP1,YP,1)
      YP = YPNRM*YP
C
C ***** SYSTEM SOLVED *****
C
C ***** UPDATE AUGMENTED SYSTEM SO THAT LAST ROW IS YP *****
C                        
C S=YP-YPOLD,  T = E(NP1)T*Q.
C      
      S = YP - YPOLD
      T = Q(NP1,:)
      CALL R1UPQF(NP1,S,T,Q,R,W)        
C
      RETURN
C
      END SUBROUTINE TANGQF
        SUBROUTINE UPQRQF(N,ETA,S,F0,F1,Q,R,W,T)
C
C SUBROUTINE  UPQRQF  PERFORMS A BROYDEN UPDATE ON THE  Q R  
C FACTORIZATION OF A MATRIX  A, (AN APPROXIMATION TO J(X0)), 
C RESULTING IN THE FACTORIZATION  Q+ R+ OF
C
C       A+  =  A  +  (Y - A*S) (ST)/(ST * S),
C
C (AN APPROXIMATION TO J(X1))
C WHERE S = X1 - X0, ST = S TRANSPOSE,  Y = F(X1) - F(X0).
C
C THE ENTRY POINT  R1UPQF  PERFORMS THE RANK ONE UPDATE ON THE QR
C FACTORIZATION OF 
C
C       A+ =  A + Q*(T*ST).
C
C
C ON INPUT:
C
C N  IS THE DIMENSION OF X AND F(X).
C
C ETA  IS A NOISE PARAMETER.  IF (Y-A*S)(I) .LE. ETA*(|F1(I)|+|F0(I)|)
C    FOR 1 .LE. I .LE. N, THEN NO UPDATE IS PERFORMED.
C
C S(1:N) = X1 - X0   (OR S FOR THE ENTRY POINT R1UPQF).
C
C F0(1:N) = F(X0).
C
C F1(1:N) = F(X1).
C
C Q(1:N,1:N)  CONTAINS THE OLD Q , WHERE  A = Q*R .
C
C R(1:N*(N+1)/2)  CONTAINS THE OLD R, STORED BY COLUMNS.
C
C W(1:N), T(1:N)  ARE WORK ARRAYS ( T  CONTAINS THE VECTOR T FOR THE
C    ENTRY POINT  R1UPQF ).
C
C 
C ON OUTPUT:
C
C N  AND  ETA  ARE UNCHANGED.
C
C Q  CONTAINS Q+ .
C
C R   CONTAINS R+, STORED BY COLUMNS.
C
C S, F0, F1, W, AND T  HAVE ALL BEEN CHANGED.
C
C
C CALLS   DGEMV, DNRM2, DTPMV.
C
C ***** DECLARATIONS *****
      USE REAL_PRECISION
C
C FUNCTION DECLARATIONS 
C
      INTERFACE
        FUNCTION DNRM2(N,X,STRIDE)
        USE REAL_PRECISION
        INTEGER:: N,STRIDE
        REAL (KIND=R8):: DNRM2,X(N)
        END FUNCTION DNRM2
      END INTERFACE
C
C LOCAL VARIABLES 
C
      REAL (KIND=R8):: C, DEN, ONE, SS, WW, YY, ZERO
      INTEGER:: I, INDEXC, INDEXD, INDXC2, J, K
      LOGICAL:: SKIPUP
C
C SCALAR ARGUMENTS 
C
      REAL (KIND=R8):: ETA
      INTEGER:: N
C
C ARRAY DECLARATIONS  
C
      REAL (KIND=R8)::  S(N), F0(N), F1(N), Q(N,N), R(N*(N+1)/2),
     &    W(N), T(N), TT(2)
C
C ***** END OF DECLARATIONS *****  
C
C ***** FIRST EXECUTABLE STATEMENT *****
C
      ONE = 1.0
      ZERO = 0.0
      SKIPUP = .TRUE.
C
C ***** DEFINE T AND S SUCH THAT *****
C
C           A+ = Q*(R + T*ST). 
C
C T = R*S.
C
      T = S
      CALL DTPMV('U','N','N',N,R,T,1)
C
C W = Y - Q*T  = Y - A*S.
C
      W = F1 - F0 - MATMUL(Q,T)
C
C IF W(I) IS NOT SMALL, THEN UPDATE MUST BE PERFORMED,
C OTHERWISE SET W(I) TO 0.
C
      WHERE (ABS(W) .LE. ETA*(ABS(F1) + ABS(F0))) W = 0.0
      IF (ANY(ABS(W) .GT. ETA*(ABS(F1) + ABS(F0)))) SKIPUP = .FALSE.
C
C IF NO UPDATE IS NECESSARY, THEN RETURN.
C
      IF (SKIPUP) RETURN
C
C T = QT*W = QT*Y - R*S.
C
      CALL DGEMV('T',N,N,ONE,Q,N,W,1,ZERO,T,1)
C
C S = S/(ST*S).
C
      S = (1.0/DOT_PRODUCT(S,S))*S
C
C ***** END OF COMPUTATION OF  T & S      *****
C       AT THIS POINT,  A+ = Q*(R + T*ST). 
C
      ENTRY R1UPQF(N,S,T,Q,R,W)
C
C ***** COMPUTE THE QR FACTORIZATION Q- R- OF (R + T*S).  THEN,  *****
C       Q+ = Q*Q-,  AND  R+ = R-.
C
C FIND THE LARGEST  K  SUCH THAT  T(K) .NE. 0.
C
      K = N
      DO
        IF (T(K) .NE. 0.0 .OR. K .LE. 1) EXIT
        K=K-1
      END DO
C
C COMPUTE THE INDEX OF R(K-1,K-1).
C         
      INDEXD = (K*(K-1))/2
C
C ***** TRANSFORM R+T*ST INTO AN UPPER HESSENBERG MATRIX *****
C
C DETERMINE JACOBI ROTATIONS WHICH WILL ZERO OUT ROWS 
C N, N-1,...,2  OF THE MATRIX  T*ST,  AND APPLY THESE
C ROTATIONS TO  R.  (THIS IS EQUIVALENT TO APPLYING THE
C SAME ROTATIONS TO  R+T*ST, EXCEPT FOR THE FIRST ROW.
C THUS, AFTER AN ADJUSTMENT FOR THE FIRST ROW, THE 
C RESULT IS AN UPPER HESSENBERG MATRIX.  THE
C SUBDIAGONAL ELEMENTS OF WHICH WILL BE STORED IN  W.
C
C NOTE:  ROWS N,N-1,...,K+1 ARE ALREADY ALL ZERO.
C
      JACOBI: DO I=K-1,1,-1
C
C         DETERMINE THE JACOBI ROTATION WHICH WILL ZERO OUT
C         ROW  I+1  OF THE  T*ST  MATRIX.
C
        IF (T(I) .EQ. 0.0) THEN
          C = 0.0
C         SS = SIGN(-T(I+1))= -T(I+1)/|T(I+1)|
          SS = -SIGN(ONE,T(I+1))
        ELSE
          DEN = DNRM2(2,T(I),1)
          C = T(I) / DEN
          SS = -T(I+1)/DEN
        END IF
C
C         PREMULTIPLY  R  BY THE JACOBI ROTATION.
C
        YY = R(INDEXD)
        WW = 0.0
        R(INDEXD) = C*YY - SS*WW
        W(I+1) = SS*YY + C*WW
        DO J= I+1,N
C           YY = R(I,J)
C           WW = R(I+1,J)
            INDEXC = ((J-1)*J)/2 + I 
            INDXC2 = INDEXC + 1
            YY = R(INDEXC)
            WW = R(INDXC2)
C           R(I,J) = C*YY - SS*WW
C           R(I+1,J) = SS*YY + C*WW
            R(INDEXC) = C*YY - SS*WW
            R(INDXC2) = SS*YY + C*WW
        END DO
C
C         MULTIPLY  Q  BY THE JACOBI ROTATION.
C
        DO J=1,N
          YY = Q(J,I)
          WW = Q(J,I+1)
          Q(J,I) = C*YY - SS*WW
          Q(J,I+1) = SS*YY + C*WW
        END DO
C
C         UPDATE  T(I)  SO THAT  T(I)*ST(J)  IS THE  (I,J)TH  COMPONENT
C         OF  T*ST, PREMULTIPLIED BY ALL OF THE JACOBI ROTATIONS SO
C         FAR.
C
        IF (T(I) .EQ. 0.0) THEN
          T(I) = ABS(T(I+1))
        ELSE
          T(I) = DNRM2(2,T(I),1)
        END IF
C
C         LET INDEXD = THE INDEX OF R(I-1,I-1).
C
        INDEXD = INDEXD - I
C
      END DO JACOBI
C
C UPDATE THE FIRST ROW OF  R  SO THAT  R  HOLDS  (R+T*ST) 
C PREMULTIPLIED BY ALL OF THE ABOVE JACOBI ROTATIONS.
C
      J=1
      DO I=1,N 
        R(J) = T(1)*S(I) + R(J)
        J=I+J
      END DO
C
C ***** END OF TRANSFORMATION TO UPPER HESSENBERG *****
C
C
C ***** TRANSFORM UPPER HESSENBERG MATRIX INTO UPPER *****
C       TRIANGULAR MATRIX. 
C
C       INDEXD = INDEX OF R(I,I).
C        
      INDEXD = 1
      HESSEN: DO I=1,K-1
C
C         DETERMINE APPROPRIATE JACOBI ROTATION TO ZERO OUT
C         R(I+1,I).
C
        IF (R(INDEXD) .EQ. 0.0) THEN
          C = 0.0
          SS = -SIGN(ONE,W(I+1))
        ELSE
          TT(1) = R(INDEXD)
          TT(2) = W(I+1)
          DEN = DNRM2(2,TT,1)
          C = R(INDEXD) / DEN
          SS = -W(I+1)/DEN
        END IF
C
C         PREMULTIPLY  R  BY JACOBI ROTATION.
C
        YY = R(INDEXD)
        WW = W(I+1)
        R(INDEXD) = C*YY - SS*WW
        W(I+1) = 0.0
        DO J= I+1,N
C           YY = R(I,J)
C           WW = R(I+1,J)  
          INDEXC = ((J-1)*J)/2 + I
          INDXC2 = INDEXC + 1 
          YY = R(INDEXC)
          WW = R(INDXC2)
C           R(I,J) = C*YY -SS*WW
C           R(I+1,J) = SS*YY + C*WW
          R(INDEXC) = C*YY - SS*WW
          R(INDXC2) = SS*YY + C*WW
        END DO
        INDEXD = INDEXD + I + 1
C
C         MULTIPLY  Q  BY JACOBI ROTATION.
C
        DO J=1,N
          YY = Q(J,I)
          WW = Q(J,I+1)
          Q(J,I) = C*YY - SS*WW
          Q(J,I+1) = SS*YY + C*WW
        END DO
      END DO HESSEN
C
C ***** END OF TRANSFORMATION TO UPPER TRIANGULAR *****
C
C
C ***** END OF UPDATE *****
C
C
      RETURN
      END SUBROUTINE UPQRQF
