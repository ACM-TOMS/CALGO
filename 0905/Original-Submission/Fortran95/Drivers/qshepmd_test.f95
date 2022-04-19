  PROGRAM QSHEPM_TEST
!
!
! Algorithm 661, Collected Algorithms from ACM.
!
! This program tests the scattered data interpolation
! package QSHEPMD by printing the maximum errors associated
! with interpolated values and gradients on a 5 by 5 by 5
! uniform grid in the unit cube. The data set consists
! of 64 nodes with data values taken from a quadratic 
! function for which the method is exact. The ratio of maximum
! interpolation error relative to the machine precision is
! also printed. This should be O(1). The interpolated
! value from QSMVAL is compared with the function for accuracy.
!
    USE SHEPPACK, ONLY: QSHEPM, QSMVAL
    USE REAL_PRECISION     
    IMPLICIT NONE
 
    INTEGER :: I, IER, J, K, L
    INTEGER, PARAMETER :: N = 64, NQ = 17, NR = 3, NW = 32

    REAL(KIND=R8) :: EPS, EQ, EQX, EQY, EQZ, Q1, &
      RMAX, RQ, YL, ZL
    REAL(KIND=R8), DIMENSION(5) :: P            
    REAL(KIND=R8), DIMENSION(N) :: F 
    REAL(KIND=R8), DIMENSION(3,N) ::  X
    REAL(KIND=R8), DIMENSION(3) ::  PX

!   INTERFACE 
!     SUBROUTINE QSHEPM (M, N, X, F, NQ, NW, NR, &
!       RMAX, IER)
!       USE REAL_PRECISION     
!       INTEGER :: IER, N, NQ, NW, NR
!       REAL(KIND=R8)::RMAX
!       REAL(KIND=R8), DIMENSION(N) :: F
!       REAL(KIND=R8), DIMENSION(3,N) :: X
!     END SUBROUTINE QSHEPM
!   END INTERFACE
 

    WRITE ( *, * ) 'QSHEPM'
    WRITE ( *, * ) 'Test for QSHEPM.'
    WRITE ( *, * ) 'Testing 3D case.'
!
! Generate a 4 by 4 by 4 grid of nodes in the unit cube.
!
    L = 0
    DO K = 1,4
      ZL = REAL (K - 1,KIND=R8) / 3.0E+00_R8
      DO J = 1, 4
        YL = REAL (J - 1,KIND=R8) / 3.0E+00_R8
        DO I = 1, 4
          L = L + 1
          X(1,L) = REAL (I - 1,KIND=R8) / 3.0E+00_R8
          X(2,L) = YL
          X(3,L) = ZL
        END DO
      END DO
    END DO
!
! Compute the data values.
!
    DO L = 1, N
      F(L) = FQ (X(1,L), X(2,L), X(3,L))
    END DO
!
! Compute parameters defining the interpolant Q.
!
    CALL QSHEPM (3, N, X, F, NQ, NW, NR, &
      RMAX, IER)

    IF (IER /= 0) THEN
      WRITE ( *, * ) ' '
      WRITE ( *, * ) 'QSHEPM'
      WRITE ( *, * ) 'Error return from QSHEPM, IER = ', IER
      STOP
    END IF
!
! Generate a 5 by 5 by 5 uniform grid of interpolation
! points (p(i),p(j),p(k)) in the unit cube.  The eight
! corners coincide with nodes.
!
    DO I = 1, 5
      P(I) = REAL (I - 1, KIND=R8) / 4.0_R8
    END DO
!
! Compute the machine precision EPS.
!
    EPS = EPSILON (EPS)
!
! Compute interpolation errors and test for agreement in the
! Q value returned by QSMVAL.
!
    EQ = 0.0_R8
    EQX = 0.0_R8
    EQY = 0.0_R8
    EQZ = 0.0_R8
    DO K = 1, 5
      PX(3) = P(K)
      DO J = 1, 5
        PX(2) = P(J)
        DO I = 1, 5
          PX(1) = P(I)
          Q1 = QSMVAL (3, N, PX, X, F, NR, RMAX)
          EQ  = MAX (EQ,  ABS (FQ(PX(1),PX(2),PX(3)) - Q1 ))
        END DO
      END DO
    END DO
!
! Print errors and the ratio eq/eps.
!
    RQ = EQ / EPS
    WRITE (*, *) 'Maximum absolute errors in the interpolant Q '
    WRITE (*, *) ' '
    WRITE (*, *) '  FUNCTION   MAX ERROR  MAX ERROR/EPS'
    WRITE (*, '(''      Q       '',E9.3,''       '',F4.2)') EQ, RQ
    WRITE (*, *) ' '
    WRITE (*, *) 'Normal end of execution.'
    STOP
!
! Quadratic test function and its partial derivatives.
!
    CONTAINS
      FUNCTION FQ(XX, YY, ZZ)
        USE REAL_PRECISION           
        REAL(KIND=R8):: FQ, XX, YY, ZZ
        FQ = ((XX + 2.0_R8 * YY + 3.0_R8 * ZZ) / 6.0_R8)**2
        RETURN
      END FUNCTION FQ
    
      FUNCTION FX(XX, YY, ZZ)
        USE REAL_PRECISION           
        REAL(KIND=R8):: FX, XX, YY, ZZ
        FX = (XX + 2.0_R8 * YY + 3.0_R8 * ZZ) / 18.0_R8
        RETURN
      END FUNCTION FX

      FUNCTION FY(XX, YY, ZZ)
        USE REAL_PRECISION     
        REAL(KIND=R8):: FY, XX, YY, ZZ
        FY = (XX + 2.0_R8 * YY + 3.0_R8 * ZZ) / 9.0_R8
        RETURN
      END FUNCTION FY

      FUNCTION FZ(XX, YY, ZZ)
        REAL(KIND=R8):: FZ, XX, YY, ZZ
        FZ = (XX + 2.0_R8 * YY + 3.0_R8 * ZZ) / 6.0_R8
        RETURN
      END FUNCTION FZ
  END PROGRAM QSHEPM_TEST
