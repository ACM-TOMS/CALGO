C  MAIN PROGRAM TO TEST FIXPQS, FIXPNS, AND FIXPDS;
C
C       THIS PROGRAM TESTS THE HOMPACK ROUTINES FIXPQS, FIXPNS, AND
C       FIXPDS.
C
C       THE OUTPUT FROM THIS ROUTINE SHOULD BE AS FOLLOWS, WITH THE
C       EXECUTION TIMES CORRESPONDING TO A DEC AXP 3000/600.
C
C       TESTING FIXPQS WITH STORAGE MODE = 1
C
C LAMBDA = 1.00000000  FLAG = 1      33 JACOBIAN EVALUATIONS
C EXECUTION TIME(SECS) =     0.119    ARC LENGTH =     1.274
C   4.00864019E-01  2.65454893E-01  8.40421103E-02  4.83042527E-01
C   3.01797132E-01  2.32508994E-01  4.96639853E-01  3.00908894E-01
C
C       TESTING FIXPNS WITH STORAGE MODE = 1
C
C LAMBDA = 1.00000000  FLAG = 1      20 JACOBIAN EVALUATIONS
C EXECUTION TIME(SECS) =     0.013    ARC LENGTH =     1.275
C   4.00864019E-01  2.65454893E-01  8.40421103E-02  4.83042527E-01
C   3.01797132E-01  2.32508994E-01  4.96639853E-01  3.00908894E-01
C
C       TESTING FIXPDS WITH STORAGE MODE = 1
C
C LAMBDA = 1.00000000  FLAG = 1      70 JACOBIAN EVALUATIONS
C EXECUTION TIME(SECS) =     0.031    ARC LENGTH =     1.281
C   4.00864019E-01  2.65454893E-01  8.40421103E-02  4.83042527E-01
C   3.01797132E-01  2.32508994E-01  4.96639853E-01  3.00908894E-01
C
C       TESTING FIXPQS WITH STORAGE MODE = 2
C
C LAMBDA = 1.00000000  FLAG = 1      33 JACOBIAN EVALUATIONS
C EXECUTION TIME(SECS) =     0.015    ARC LENGTH =     1.274
C   4.00864019E-01  2.65454893E-01  8.40421103E-02  4.83042527E-01
C   3.01797132E-01  2.32508994E-01  4.96639853E-01  3.00908894E-01
C
C       TESTING FIXPNS WITH STORAGE MODE = 2
C
C LAMBDA = 1.00000000  FLAG = 1      20 JACOBIAN EVALUATIONS
C EXECUTION TIME(SECS) =     0.011    ARC LENGTH =     1.275
C   4.00864019E-01  2.65454893E-01  8.40421103E-02  4.83042527E-01
C   3.01797132E-01  2.32508994E-01  4.96639853E-01  3.00908894E-01
C
C       TESTING FIXPDS WITH STORAGE MODE = 2
C
C LAMBDA = 1.00000000  FLAG = 1      70 JACOBIAN EVALUATIONS
C EXECUTION TIME(SECS) =     0.022    ARC LENGTH =     1.281
C   4.00864019E-01  2.65454893E-01  8.40421103E-02  4.83042527E-01
C   3.01797132E-01  2.32508994E-01  4.96639853E-01  3.00908894E-01
C 
C
      MODULE SWITCH
C
C  ROWSET  IS USED TO INITIALIZE SPARSE MATRIX DATA STRUCTURES ONLY
C  ONCE, AFTER THEY ARE ALLOCATED.
C
      LOGICAL::  ROWSET
      END MODULE SWITCH
!
      PROGRAM TESTS
      USE SWITCH
      USE REAL_PRECISION, ONLY : R8
      USE HOMPACK90, ONLY : FIXPDS, FIXPNS, FIXPQS
      IMPLICIT NONE
      INTEGER, PARAMETER:: N=8, NDIMA=8
      REAL (KIND=R8):: A(N),ANSAE,ANSRE,ARCAE,ARCRE,
     &  ARCLEN,DTIME,SSPAR(8),Y(N+1)
      INTEGER:: IFLAG,II,J,LENQR,MODE,NFE,NP1,TIMENEW(8),
     &  TIMEOLD(8),TRACE
      CHARACTER (LEN=6):: NAME
! If using a subroutine library of the HOMPACK90 subroutines rather than
! the MODULE HOMPACK90 (as above), then the following INTERFACE
! statements are necessary.
!     INTERFACE
!       SUBROUTINE FIXPDS(N,Y,IFLAG,ARCTOL,EPS,TRACE,A,NDIMA,
!    &    NFE,ARCLEN,MODE,LENQR)
!       USE REAL_PRECISION
!       INTEGER, INTENT(IN)::LENQR,MODE,N,NDIMA,TRACE
!       REAL (KIND=R8), DIMENSION(:), INTENT(IN OUT)::A,Y
!       INTEGER, INTENT(IN OUT)::IFLAG
!       REAL (KIND=R8), INTENT(IN OUT)::ARCTOL,EPS
!       INTEGER, INTENT(OUT)::NFE
!       REAL (KIND=R8), INTENT(OUT)::ARCLEN
!       END SUBROUTINE FIXPDS
C
!       SUBROUTINE FIXPNS(N,Y,IFLAG,ARCRE,ARCAE,ANSRE,ANSAE,TRACE,A,
!    &    NFE,ARCLEN,MODE,LENQR,SSPAR)
!       USE REAL_PRECISION
!       INTEGER, INTENT(IN)::LENQR,MODE,N,TRACE
!       REAL (KIND=R8), DIMENSION(:), INTENT(IN OUT)::A,Y
!       INTEGER, INTENT(IN OUT)::IFLAG
!       REAL (KIND=R8), INTENT(IN OUT)::ANSAE,ANSRE,ARCAE,ARCRE,
!    &    SSPAR(8)
!       INTEGER, INTENT(OUT)::NFE
!       REAL (KIND=R8), INTENT(OUT)::ARCLEN
!       END SUBROUTINE FIXPNS
C
!       SUBROUTINE FIXPQS(N,Y,IFLAG,ARCRE,ARCAE,ANSRE,ANSAE,TRACE,A,
!    &    NFE,ARCLEN,MODE,LENQR,SSPAR)
!       USE REAL_PRECISION
!       INTEGER, INTENT(IN)::LENQR,MODE,N,TRACE
!       REAL (KIND=R8), DIMENSION(:), INTENT(IN OUT)::A,Y
!       INTEGER, INTENT(IN OUT)::IFLAG
!       REAL (KIND=R8), INTENT(IN OUT)::ANSAE,ANSRE,ARCAE,ARCRE,
!    &    SSPAR(4)
!       INTEGER, INTENT(OUT)::NFE
!       REAL (KIND=R8), INTENT(OUT)::ARCLEN
!       END SUBROUTINE FIXPQS
!     END INTERFACE
C
C TEST EACH OF THE THREE ALGORITHMS WITH BOTH STORAGE MODES.
      DO MODE=1,2
        SELECT CASE (MODE)
          CASE(1)
            LENQR=18
          CASE(2)
            LENQR=36
        END SELECT
      DO II=1,3
C
C DEFINE ARGUMENTS FOR CALL TO HOMPACK PROCEDURE.
C
         NP1=N+1
         ARCRE=0.5D-4
         ARCAE=0.5D-4
         ANSRE=1.0D-12
         ANSAE=1.0D-12
         TRACE=0
         SSPAR=0.0
         IFLAG=-MODE
         Y(1:N)=0.5_R8
         IF(IFLAG .EQ. -2) A=Y(1:N)
         ROWSET = .FALSE.
C
C GET CURRENT DATE AND TIME.
C
         CALL DATE_AND_TIME(VALUES=TIMEOLD)
C
C CALL TO HOMPACK ROUTINE.
C
        IF (II .EQ. 1) THEN
          NAME='FIXPQS'
          CALL FIXPQS(N,Y,IFLAG,ARCRE,ARCAE,ANSRE,ANSAE,TRACE,A,
     &       NFE,ARCLEN,MODE,LENQR,SSPAR)
        ELSE IF (II .EQ. 2) THEN
          NAME='FIXPNS'
          CALL FIXPNS(N,Y,IFLAG,ARCRE,ARCAE,ANSRE,ANSAE,TRACE,A,
     &       NFE,ARCLEN,MODE,LENQR,SSPAR)
        ELSE
          NAME='FIXPDS'
          CALL FIXPDS(N,Y,IFLAG,ARCRE,ANSRE,TRACE,A,NDIMA,NFE,
     &      ARCLEN,MODE,LENQR)
        END IF
C
C CALCULATE EXECUTION TIME.
C
        CALL DATE_AND_TIME(VALUES=TIMENEW)
        IF (TIMENEW(8) .LT. TIMEOLD(8)) THEN
          TIMENEW(8)=TIMENEW(8)+1000
          TIMENEW(7)=TIMENEW(7)-1
        ENDIF
        IF (TIMENEW(7) .LT. TIMEOLD(7)) THEN
          TIMENEW(7)=TIMENEW(7)+60
          TIMENEW(6)=TIMENEW(6)-1
        ENDIF
        IF (TIMENEW(6) .LT. TIMEOLD(6)) THEN
          TIMENEW(6)=TIMENEW(6)+60
          TIMENEW(5)=TIMENEW(5)-1
        ENDIF
        IF (TIMENEW(5) .LT. TIMEOLD(5)) TIMENEW(5)=TIMENEW(5)+24
        DTIME=DOT_PRODUCT(TIMENEW(5:8)-TIMEOLD(5:8),
     &    (/3600000,60000,1000,1/) )/1000.0
C
        WRITE (6,45) NAME, MODE
45      FORMAT (//,7X,'TESTING',1X,A6,' WITH STORAGE MODE =',I2)
        WRITE (6,50) Y(NP1),IFLAG,NFE,DTIME,ARCLEN,(Y(J),J=1,N)
50      FORMAT(/' LAMBDA =',F11.8,'  FLAG =',I2,I8,' JACOBIAN ',
     &    'EVALUATIONS',/,1X,'EXECUTION TIME(SECS) =',F10.3,4X,
     &    'ARC LENGTH =',F10.3/(1X,4ES16.8))
      END DO
      END DO
      STOP
      END PROGRAM TESTS
!
! SAMPLE USER WRITTEN HOMOTOPY SUBROUTINES FOR TESTING FIXP*S.
!
      SUBROUTINE F(X,V)
C***********************************************************************
C
C SUBROUTINE F(X,V) COMPUTES F AT THE POINT X, RETURNING THE VALUE IN V.
C
C***********************************************************************
      USE REAL_PRECISION, ONLY : R8
      IMPLICIT NONE
      REAL (KIND=R8), INTENT(IN):: X(:)
      REAL (KIND=R8), INTENT(OUT):: V(:)
       V(1)=X(1)**3+6.0*X(2)*X(3)-1+2.0*X(1)
       V(2)=6.0*X(1)*X(3)+X(2)**4*X(5)-1+3.0*X(2)
       V(3)=6.0*X(1)*X(2)+X(3)*X(5)-1+4.0*X(3)
       V(4)=X(4)**3*X(8)-1+2.0*X(4)
       V(5)=X(2)**5/5.0 + X(3)**2/2.0 + X(8)*X(5)-1+3.0*X(5)
       V(6)=X(6)*X(8)-1+4.0*X(6)
       V(7)=X(7)**2*X(8)**3-1+2.0*X(7)
       V(8)=X(4)**4/4.0 + X(5)**2/2.0 + X(6)**2/2.0 + X(7)**3*
     &    X(8)**2-1+3.0*X(8)
      RETURN
      END SUBROUTINE F
C
      SUBROUTINE FJACS(X)
C******************************************************************
C
C SUBROUTINE FJACS(X) COMPUTES THE JACOBIAN MATRIX OF F AT THE POINT
C X, RETURNING THE JACOBIAN MATRIX IN PACKED SKYLINE FORM (MODE=1)
C IN THE ARRAYS  QRSPARSE  AND  ROWPOS .
C
C*****************************************************************
      USE REAL_PRECISION, ONLY : R8
      USE HOMPACK90_GLOBAL, ONLY: QRSPARSE, ROWPOS, COLPOS
      USE SWITCH
      REAL (KIND=R8), INTENT(IN):: X(:)
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
C
      INTEGER:: N
      N=SIZE(X)
      IF (.NOT. ROWSET) THEN
        ROWSET=.TRUE.
        ROWPOS(1:N+1) = (/ 1,2,4,7,8,12,13,14,19 /)
      END IF
      QRSPARSE(1)=3.0*X(1)**2+2.0
      QRSPARSE(2)=4.0*X(2)**3*X(5)+3.0
      QRSPARSE(3)=6.0*X(3)
      QRSPARSE(4)=X(5)+4.0
      QRSPARSE(5)=6.0*X(1)
      QRSPARSE(6)=6.0*X(2)
      QRSPARSE(7)=3.0*X(4)**2*X(8)+2.0
      QRSPARSE(8)=X(8)+3.0
      QRSPARSE(9)=0.0
      QRSPARSE(10)=X(3)
      QRSPARSE(11)=X(2)**4
      QRSPARSE(12)=X(8)+4.0
      QRSPARSE(13)=2.0*X(7)*X(8)**3+2.0
      QRSPARSE(14)=2.0*X(7)**3*X(8)+3.0
      QRSPARSE(15)=3.0*X(7)**2*X(8)**2
      QRSPARSE(16)=X(6)
      QRSPARSE(17)=X(5)
      QRSPARSE(18)=X(4)**3
      RETURN
      END SUBROUTINE FJACS
C
      SUBROUTINE RHO(A,LAMBDA,X,V)
      USE REAL_PRECISION, ONLY : R8
      REAL (KIND=R8), INTENT(IN):: A(:),X(:)
      REAL (KIND=R8), INTENT(IN OUT):: LAMBDA
      REAL (KIND=R8), INTENT(OUT):: V(:)
      INTERFACE
        SUBROUTINE F(X,V)
        USE REAL_PRECISION
        REAL (KIND=R8), DIMENSION(:), INTENT(IN):: X
        REAL (KIND=R8), DIMENSION(:), INTENT(OUT):: V
        END SUBROUTINE F
      END INTERFACE
C
C   EVALUATE  RHO(A,LAMBDA,X)  AND RETURN IN THE VECTOR  V .
C
      INTEGER:: N
      N=SIZE(X)
      CALL F(X(1:N), V(1:N))
      V(1:N) = LAMBDA*V(1:N) + (1.0 - LAMBDA)*(X(1:N) - A(1:N))
C
      RETURN
      END SUBROUTINE RHO
C
      SUBROUTINE RHOA(A,LAMBDA,X)
      USE REAL_PRECISION, ONLY : R8
      REAL (KIND=R8), INTENT(OUT):: A(:)
      REAL (KIND=R8), INTENT(IN):: LAMBDA,X(:)
      INTERFACE
        SUBROUTINE F(X,V)
        USE REAL_PRECISION
        REAL (KIND=R8), DIMENSION(:), INTENT(IN):: X
        REAL (KIND=R8), DIMENSION(:), INTENT(OUT):: V
        END SUBROUTINE F
      END INTERFACE
C
C   CALCULATE AND RETURN IN  A  THE VECTOR Z SUCH THAT
C   RHO(Z,LAMBDA,X) = 0 .
C   N=NDIMA FOR THIS TEST PROBLEM.
C
      INTEGER:: N
      N=SIZE(X)
      CALL F(X(1:N),A(1:N))
      A(1:N)=LAMBDA*A(1:N)/(1.0 - LAMBDA) + X(1:N)
      RETURN
      END SUBROUTINE RHOA
C
      SUBROUTINE RHOJS(A,LAMBDA,X)
C*****************************************************************
C
C Subroutine RHOJS(A,LAMBDA,X) computes the Jacobian matrix of
C rho(a,x,lambda) = lambda*F(x) + (1 - lambda)*(x - a) at the 
C point (A,X,LAMBDA), returning the Jacobian matrix in sparse row
C storage format (MODE = 2) in the arrays QRSPARSE, ROWPOS, and
C COLPOS.
C
C*****************************************************************
      USE REAL_PRECISION, ONLY : R8
      USE HOMPACK90_GLOBAL, ONLY: QRSPARSE, ROWPOS, COLPOS
      USE SWITCH
      REAL (KIND=R8), INTENT(IN):: A(:),LAMBDA,X(:)
      INTERFACE
        SUBROUTINE F(X,V)
        USE REAL_PRECISION
        REAL (KIND=R8), DIMENSION(:), INTENT(IN):: X
        REAL (KIND=R8), DIMENSION(:), INTENT(OUT):: V
        END SUBROUTINE F
      END INTERFACE
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
C
C---------------------------------------------------------------------
C    [QRSPARSE] = [ LAMBDA*DF(X) + (1 - LAMBDA)*I | F(X) - X + A ] .  
C---------------------------------------------------------------------
C
C  LOCAL VARIABLES
C
      INTEGER, PARAMETER:: N=8
      INTEGER:: J, JPOS, ELEM(N) = (/4,9,14,18,24,27,30,36/)
      REAL (KIND=R8):: DRHODL(N)
C
      IF (.NOT. ROWSET) THEN  
        ROWSET=.TRUE.
        ROWPOS(1:N+1) = (/ 1,5,10,15,19,25,28,31,37 /)
        COLPOS(1:36) = (/1,2,3,9,1,2,3,5,9,1,2,3,5,9,4,5,8,9,
     &                   2,3,4,5,8,9,6,8,9,7,8,9,4,5,6,7,8,9/)                
      END IF
C
      QRSPARSE = 0.0
C   ROW 1.
      QRSPARSE(1)=3.0*X(1)**2+2.0
      QRSPARSE(2)=6.0*X(3)
      QRSPARSE(3)=6.0*X(2)
C   ROW 2.
      QRSPARSE(5)=6.0*X(3)
      QRSPARSE(6)=4.0*X(2)**3*X(5)+3.0
      QRSPARSE(7)=6.0*X(1)
      QRSPARSE(8)=X(2)**4
C   ROW 3.
      QRSPARSE(10)=6.0*X(2)
      QRSPARSE(11)=6.0*X(1)
      QRSPARSE(12)=X(5)+4.0
      QRSPARSE(13)=X(3)
C   ROW 4.
      QRSPARSE(15)=3.0*X(4)**2*X(8)+2.0
      QRSPARSE(16)=0.0
      QRSPARSE(17)=X(4)**3
C   ROW 5.
      QRSPARSE(19)=X(2)**4
      QRSPARSE(20)=X(3)
      QRSPARSE(21)=0.0
      QRSPARSE(22)=X(8)+3.0
      QRSPARSE(23)=X(5)
C   ROW 6.
      QRSPARSE(25)=X(8)+4.0
      QRSPARSE(26)=X(6)
      COLPOS(25)=6
      COLPOS(26)=8
      COLPOS(27)=9
C   ROW 7.
      QRSPARSE(28)=2.0*X(7)*X(8)**3+2.0
      QRSPARSE(29)=3.0*X(7)**2*X(8)**2
C   ROW 8.
      QRSPARSE(31)=X(4)**3
      QRSPARSE(32)=X(5)
      QRSPARSE(33)=X(6)
      QRSPARSE(34)=3.0*X(7)**2*X(8)**2
      QRSPARSE(35)=2.0*X(7)**3*X(8)+3.0
C
      QRSPARSE=LAMBDA*QRSPARSE
C
C   FIND INDEX JPOS OF DIAGONAL ELEMENT IN JTH ROW OF QR.
C
          DO J=1,N
            JPOS=ROWPOS(J)
            DO
              IF (COLPOS(JPOS) .EQ. J) EXIT
              JPOS=JPOS+1
            END DO
            QRSPARSE(JPOS) = QRSPARSE(JPOS) + 1.0 - LAMBDA
          END DO
C
C   INITIALIZE (N+1)ST COLUMN.
C
      CALL F(X(1:N),DRHODL(1:N))
      DRHODL = DRHODL - X(1:N) + A(1:N)
      QRSPARSE(ELEM) = DRHODL(1:8)
C
      RETURN
      END SUBROUTINE RHOJS
C **********************************************************************
C
C THE REST OF THESE SUBROUTINES ARE NOT USED BY PROGRAM TESTS, AND ARE
C INCLUDED HERE SIMPLY FOR COMPLETENESS AND AS TEMPLATES FOR THEIR USE.
C
      SUBROUTINE FJAC(X,V,K)
      USE REAL_PRECISION, ONLY : R8
      REAL (KIND=R8), INTENT(IN):: X(:)
      REAL (KIND=R8), INTENT(OUT):: V(:)
      INTEGER, INTENT(IN):: K
C
C RETURN IN  V  THE KTH COLUMN OF THE JACOBIAN MATRIX OF
C F(X) EVALUATED AT  X .
C
      V(1)=X(1) ! INTENT(OUT) VARIABLE MUST BE DEFINED.
      RETURN
      END SUBROUTINE FJAC
C
      SUBROUTINE RHOJAC(A,LAMBDA,X,V,K)
      USE REAL_PRECISION, ONLY : R8
      USE HOMPACK90_GLOBAL, ONLY: IPAR, PAR
      REAL (KIND=R8), INTENT(IN):: A(:),X(:)
      REAL (KIND=R8), INTENT(IN OUT):: LAMBDA
      REAL (KIND=R8), INTENT(OUT):: V(:)
      INTEGER, INTENT(IN):: K
C
C RETURN IN THE VECTOR  V  THE KTH COLUMN OF THE JACOBIAN
C MATRIX [D RHO/D LAMBDA, D RHO/DX] EVALUATED AT THE POINT
C (A, LAMBDA, X).
C
C THE FOLLOWING CODE IS SPECIFICALLY FOR THE POLYNOMIAL SYSTEM DRIVER
C POLSYS1H , AND SHOULD BE USED VERBATIM WITH  POLSYS1H .  IF THE USER IS
C CALLING  FIXP??  OR   STEP??  DIRECTLY, HE MUST SUPPLY APPROPRIATE
C REPLACEMENT CODE HERE.
      INTERFACE
        SUBROUTINE HFUNP(N,A,LAMBDA,X)
        USE REAL_PRECISION
        INTEGER, INTENT(IN):: N
        REAL (KIND=R8), INTENT(IN):: A(2*N),LAMBDA,X(2*N)
        END SUBROUTINE HFUNP
      END INTERFACE
      INTEGER:: J,NPOL,N2
      NPOL=IPAR(1)
      N2=2*NPOL
      IF (K .EQ. 1) THEN
C FORCE PREDICTED POINT TO HAVE  LAMBDA .GE. 0  .
        IF (LAMBDA .LT. 0.0) LAMBDA=0.0
!       CALL HFUNP(NPOL,A,LAMBDA,X)
        DO J=1,N2
          V(J)=PAR(IPAR(3 + (6-1)) + (J-1))
        END DO
        RETURN
      ELSE
        DO J=1,N2
          V(J)=PAR(IPAR(3 + (5-1)) + (J-1) + N2*(K-2))
        END DO
      ENDIF
C
      RETURN
      END SUBROUTINE RHOJAC
