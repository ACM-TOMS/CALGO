! Template for user written subroutines.
!
!  USE statements for the modules REAL_PRECISION and HOMPACK90_GLOBAL
!  should always be present in the user's subroutines.  The user
!  written subroutines must conform to the interfaces in the module
!  HOMOTOPY.
!
C
C All data and subroutines defining the problem should be coded here.
C
      SUBROUTINE F(X,V)
      USE REAL_PRECISION, ONLY : R8
      REAL (KIND=R8), DIMENSION(:), INTENT(IN):: X(:)
      REAL (KIND=R8), DIMENSION(:), INTENT(OUT):: V(:)
C
C Evaluate  F(X)  and return in the vector  V .
C
      RETURN
      END SUBROUTINE F

      SUBROUTINE FJAC(X,V,K)
      USE REAL_PRECISION, ONLY : R8
      REAL (KIND=R8), DIMENSION(:), INTENT(IN):: X(:)
      REAL (KIND=R8), DIMENSION(:), INTENT(OUT):: V(:)
      INTEGER, INTENT(IN):: K
C
C Return in  V  the Kth column of the Jacobian matrix of
C F(X) evaluated at  X .
C
      RETURN
      END SUBROUTINE FJAC

      SUBROUTINE RHO(A,LAMBDA,X,V)
      USE REAL_PRECISION, ONLY : R8
      USE HOMPACK90_GLOBAL, ONLY: IPAR, PAR  ! for POLSYS1H only.
      REAL (KIND=R8), INTENT(IN):: A(:),X(:)
      REAL (KIND=R8), INTENT(IN OUT):: LAMBDA
      REAL (KIND=R8), INTENT(OUT):: V(:)
C
C Evaluate  RHO(A,LAMBDA,X)  and return in the vector  V .
C
C The following code is specifically for the polynomial system driver
C  POLSYS1H , and should be used verbatim with  POLSYS1H .  If the user is
C calling  FIXP??  or   STEP??  directly, he must supply appropriate
C replacement code here.
      INTERFACE
        SUBROUTINE HFUNP(N,A,LAMBDA,X)
        USE REAL_PRECISION
        INTEGER, INTENT(IN):: N
        REAL (KIND=R8), INTENT(IN):: A(2*N),LAMBDA,X(2*N)
        END SUBROUTINE HFUNP
      END INTERFACE
      INTEGER:: J,NPOL
C Force predicted point to have  LAMBDA .GE. 0  .
      IF (LAMBDA .LT. 0.0) LAMBDA=0.0
      NPOL=IPAR(1)
      CALL HFUNP(NPOL,A,LAMBDA,X)
      DO J=1,2*NPOL
        V(J)=PAR(IPAR(3 + (4-1)) + (J-1))
      END DO
C
      RETURN
      END SUBROUTINE RHO

      SUBROUTINE RHOA(A,LAMBDA,X)
      USE REAL_PRECISION, ONLY : R8
      REAL (KIND=R8), INTENT(OUT):: A(:)
      REAL (KIND=R8), INTENT(IN):: LAMBDA,X(:)
C
C Calculate and return in  A  the vector Z such that
C  RHO(Z,LAMBDA,X) = 0 .
C
      RETURN
      END SUBROUTINE RHOA

      SUBROUTINE RHOJAC(A,LAMBDA,X,V,K)
      USE REAL_PRECISION, ONLY : R8
      USE HOMPACK90_GLOBAL, ONLY: IPAR, PAR  ! for POLSYS1H only.
      REAL (KIND=R8), INTENT(IN):: A(:),X(:)
      REAL (KIND=R8), INTENT(IN OUT):: LAMBDA
      REAL (KIND=R8), INTENT(OUT):: V(:)
      INTEGER, INTENT(IN):: K
C
C Return in the vector  V  the Kth column of the Jacobian
C matrix [D RHO/D LAMBDA, D RHO/DX] evaluated at the point
C (A, LAMBDA, X).
C
C The following code is specifically for the polynomial system driver
C  POLSYS1H , and should be used verbatim with  POLSYS1H .  If the user is
C calling  FIXP??  or   STEP??  directly, he must supply appropriate
C replacement code here.
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
C Force predicted point to have  LAMBDA .GE. 0  .
        IF (LAMBDA .LT. 0.0) LAMBDA=0.0
        CALL HFUNP(NPOL,A,LAMBDA,X)
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

      SUBROUTINE FJACS(X)
      USE REAL_PRECISION, ONLY : R8
      USE HOMPACK90_GLOBAL
      REAL (KIND=R8), DIMENSION(:), INTENT(IN):: X
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
      RETURN
      END SUBROUTINE FJACS

      SUBROUTINE RHOJS(A,LAMBDA,X)
      USE REAL_PRECISION, ONLY : R8
      USE HOMPACK90_GLOBAL
      REAL (KIND=R8), INTENT(IN):: A(:),LAMBDA,X(:)
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
      RETURN
      END SUBROUTINE RHOJS
