  PROGRAM QSHEP3_TEST
!
! DRIVER_QSHEP3 tests QSHEP3.
!
! Algorithm 661, Collected Algorithms from ACM.
!
! This program tests the scattered data interpolation
! package QSHEP3D by printing the maximum errors associated
! with interpolated values and gradients on a 5 by 5 by 5
! uniform grid in the unit cube. The data set consists
! of 64 nodes with data values taken from a quadratic 
! function for which the method is exact. The ratio of maximum
! interpolation error relative to the machine precision is
! also printed. This should be O(1). The interpolated
! values from QS3VAL and QS3GRD are compared for agreement.
!
    USE SHEPPACK, ONLY: QSHEP3, QS3VAL, QS3GRD
    USE REAL_PRECISION     
    IMPLICIT NONE
 
    INTEGER :: I, IER, J, K, L
    INTEGER, PARAMETER :: N = 64, NQ = 17, NR = 3, NW = 32
    INTEGER, DIMENSION (N) :: LNEXT           
    INTEGER, DIMENSION (NR,NR,NR) :: LCELL    

    REAL(KIND=R8) :: EPS, EQ, EQX, EQY, EQZ, PX, PY, PZ, Q, Q1, &
      QX, QY, QZ, RMAX, RQ, YL, ZL
    REAL(KIND=R8), DIMENSION(5) :: P            
    REAL(KIND=R8), DIMENSION(N) :: F, RSQ, X
    REAL(KIND=R8), DIMENSION(3) :: XYZDEL, XYZMIN
    REAL(KIND=R8), DIMENSION(N) :: Y, Z   
    REAL(KIND=R8), DIMENSION(9, N) :: A           

!   INTERFACE 
!     SUBROUTINE QSHEP3 (N, X, Y, Z, F, NQ, NW, NR, LCELL, LNEXT, &
!       XYZMIN, XYZDEL, RMAX, RSQ, A, IER)
!       USE REAL_PRECISION     
!       INTEGER :: IER, N, NQ, NW, NR
!       INTEGER, DIMENSION (N) :: LNEXT           
!       INTEGER, DIMENSION (NR,NR,NR) :: LCELL    
!       REAL(KIND=R8)::RMAX
!       REAL(KIND=R8), DIMENSION(N) :: F, RSQ, X
!       REAL(KIND=R8), DIMENSION(3) :: XYZDEL, XYZMIN
!       REAL(KIND=R8), DIMENSION(N) :: Y, Z   
!       REAL(KIND=R8), DIMENSION(9, N) :: A
!     END SUBROUTINE QSHEP3
!
!     SUBROUTINE QS3GRD (PX, PY, PZ, N, X, Y, Z, F, NR, LCELL, LNEXT, &
!       XYZMIN, XYZDEL, RMAX, RSQ, A, Q, QX, QY, QZ, IER)
!       USE REAL_PRECISION     
!       INTEGER :: IER, N, NR
!       INTEGER, DIMENSION (N) :: LNEXT           
!       INTEGER, DIMENSION (NR,NR,NR) :: LCELL    
!       REAL(KIND=R8):: PX, PY, PZ, Q, QX, QY, QZ, RMAX
!       REAL(KIND=R8), DIMENSION(N) :: F, RSQ, X
!       REAL(KIND=R8), DIMENSION(3) :: XYZDEL, XYZMIN
!       REAL(KIND=R8), DIMENSION(N) :: Y, Z   
!       REAL(KIND=R8), DIMENSION(9, N) :: A      
!     END SUBROUTINE QS3GRD
!   END INTERFACE

    WRITE ( *, * ) 'QSHEP3'
    WRITE ( *, * ) 'Test for QSHEP3.'
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
          X(L) = REAL (I - 1,KIND=R8) / 3.0E+00_R8
          Y(L) = YL
          Z(L) = ZL
        END DO
      END DO
    END DO
!
! Compute the data values.
!
    DO L = 1, N
      F(L) = FQ (X(L), Y(L), Z(L))
    END DO
!
! Compute parameters defining the interpolant Q.
!
    CALL QSHEP3 (N, X, Y, Z, F, NQ, NW, NR, LCELL, LNEXT, XYZMIN, &
      XYZDEL, RMAX, RSQ, A, IER)

    IF (IER /= 0) THEN
      WRITE ( *, * ) ' '
      WRITE ( *, * ) 'QSHEP3'
      WRITE ( *, * ) 'Error return from QSHEP3, IER = ', IER
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
! Q values returned by QS3VAL and QS3GRD.
!
    EQ = 0.0_R8
    EQX = 0.0_R8
    EQY = 0.0_R8
    EQZ = 0.0_R8
    DO K = 1, 5
      PZ = P(K)
      DO J = 1, 5
        PY = P(J)
        DO I = 1, 5
          PX = P(I)
          Q1 = QS3VAL (PX, PY, PZ, N, X, Y, Z, F, NR, LCELL, LNEXT, &
            XYZMIN, XYZDEL, RMAX, RSQ, A)
          CALL QS3GRD (PX, PY, PZ, N, X, Y, Z, F, NR, LCELL, LNEXT, &
            XYZMIN, XYZDEL, RMAX, RSQ, A, Q, QX, QY, QZ, IER)
          IF (IER /= 0) THEN
            WRITE (*, *) ' '
            WRITE (*, *) 'QSHEP3'
            WRITE (*, *) 'Error return from QS3GRD, IER = ', IER
            STOP
          END IF
          IF (ABS (Q1 - Q) > 3.0_R8 * ABS (Q) * EPS) THEN
            WRITE (*, *) ' '
            WRITE (*, *) 'QSHEP3 - Error.'
            WRITE (*, *) 'The interpolated values Q1 (QS3VAL) and'
            WRITE (*, *) 'Q (QS3GRD) differ.'
            WRITE (*, *) 'Q1 = ', Q1
            WRITE (*, *) 'Q  = ', Q
            STOP
          END IF
          EQ  = MAX (EQ,  ABS (FQ(PX,PY,PZ) - Q ))
          EQX = MAX (EQX, ABS (FX(PX,PY,PZ) - QX))
          EQY = MAX (EQY, ABS (FY(PX,PY,PZ) - QY))
          EQZ = MAX (EQZ, ABS (FZ(PX,PY,PZ) - QZ))
        END DO
      END DO
    END DO
!
! Print errors and the ratio eq/eps.
!
    RQ = EQ / EPS
    WRITE (*, *) 'Maximum absolute errors in the interpolant Q '
    WRITE (*, *) 'and partial derivatives (QX,QY,QZ) relative '
    WRITE (*, *) 'to machine precision EPS.'
    WRITE (*, *) ' '
    WRITE (*, *) '  FUNCTION   MAX ERROR  MAX ERROR/EPS'
    WRITE (*, '(''      Q       '',E9.3,''       '',F4.2)') EQ, RQ
    WRITE (*, '(''      QX      '',E9.3)') EQX
    WRITE (*, '(''      QY      '',E9.3)') EQY
    WRITE (*, '(''      QZ      '',E9.3)') EQZ
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
  END PROGRAM QSHEP3_TEST
