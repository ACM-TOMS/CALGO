  PROGRAM QSHEP2_TEST
!
! QSHEP2_TEST is a test for QSHEP2.
!
! Author:
!   Robert Renka
!   University of North Texas
!
! Modified by:
!   Manjula Iyer
!   Virginia Polytechnic Institute and State University
! 
! Modified:   
!   15 November 2005  
!
    USE REAL_PRECISION
    USE SHEPPACK, ONLY: QSHEP2, QS2VAL, QS2GRD
    IMPLICIT NONE

    INTEGER :: I, IER, J, K
    INTEGER, PARAMETER :: N = 36
    INTEGER, PARAMETER :: NQ = 13
    INTEGER, PARAMETER :: NR = 3
    INTEGER, PARAMETER :: NW = 19
    INTEGER, DIMENSION(N) :: LNEXT
    INTEGER, DIMENSION(NR,NR) :: LCELL
  
    REAL(KIND=R8) :: DX, DY, EPS, EQ, EQX, EQY, PX, PY, Q, Q1, &
      QX, QY, RMAX, RQ, XMIN, YMIN
    REAL(KIND=R8), DIMENSION(10) :: P
    REAL(KIND=R8), DIMENSION(N) :: F, RSQ, X, Y
    REAL(KIND=R8), DIMENSION(5,N) :: A
 
    !INTERFACE
    !  SUBROUTINE QSHEP2 (N, X, Y, F, NQ, NW, NR, LCELL, LNEXT, XMIN, &
    !    YMIN, DX, DY, RMAX, RSQ, A, IER)
    !    USE REAL_PRECISION
    !    INTEGER  :: IER, N, NQ, NR, NW
    !    INTEGER, DIMENSION(N) :: LNEXT
    !    INTEGER, DIMENSION(NR,NR) :: LCELL
    !    REAL(KIND=R8):: DX, DY, RMAX, XMIN, YMIN
    !    REAL(KIND=R8), DIMENSION(N) :: F, RSQ, X, Y
    !    REAL(KIND=R8), DIMENSION(5,N) :: A
    !  END SUBROUTINE
      
  !    SUBROUTINE QS2GRD (PX, PY, N, X, Y, F, NR, LCELL, LNEXT, XMIN, &
  !      YMIN, DX, DY, RMAX, RSQ, A, Q, QX, QY, IER)
  !      USE REAL_PRECISION
  !      INTEGER :: IER, N, NR
  !      INTEGER, DIMENSION(N) :: LNEXT
  !      INTEGER, DIMENSION(NR, NR) :: LCELL
  !      REAL(KIND=R8) :: DX, DY, PX, PY, Q, QX, QY, RMAX, XMIN, YMIN
  !      REAL(KIND=R8), DIMENSION(N) :: F, RSQ, X, Y
  !      REAL(KIND=R8), DIMENSION(5,N) :: A
  !    END SUBROUTINE
  ! END INTERFACE
 
    WRITE (*, *) ' '
    WRITE (*, *) 'Test for QSHEP2.'
    WRITE (*, *) ' '
!
! Generate a 6 by 6 grid of nodes in the unit square with
! the natural ordering.
!
    K = 0
    DO J = 5, 0, -1
      DO I = 0, 5
        K = K + 1
        X(K) = REAL(I, KIND=R8) / 5.0E+00_R8
        Y(K) = REAL(J, KIND=R8) / 5.0E+00_R8
      END DO
    END DO
!
! Compute the data values.
!
    DO K = 1, N
      F(K) = FQ (X(K), Y(K))
    END DO
!
! Call QSHEP2 to define the interpolant Q to the data.
!
    CALL QSHEP2 (N, X, Y, F, NQ, NW, NR, LCELL, LNEXT, XMIN, YMIN, &
      DX, DY, RMAX, RSQ, A, IER)

    IF (IER /= 0) THEN
      WRITE (*, *) 'QSHEP2_TEST - Error!'
      WRITE (*, *) ' Error in QSHEP2, IER = ', IER
      STOP
    END IF
!
! Generate a 10 by 10 uniform grid of interpolation points
! (P(I),P(J)) in the unit square.  
!
    DO I = 1, 10
      P(I) = REAL (I - 1, KIND=R8) / 9.0E+00_R8
    END DO
!
! Compute the machine precision EPS.
!
    EPS = EPSILON (EPS)
!
! Compute the interpolation errors.
!
    EQ = 0.0E+00_R8
    EQX = 0.0E+00_R8
    EQY = 0.0E+00_R8

    DO J = 1, 10
      PY = P(J)
      DO I = 1, 10
        PX = P(I)
        Q1 = QS2VAL (PX, PY, N, X, Y, F, NR, LCELL, LNEXT, XMIN, &
          YMIN, DX, DY, RMAX, RSQ, A)
        CALL QS2GRD (PX, PY, N, X, Y, F, NR, LCELL, LNEXT, XMIN, &
          YMIN, DX, DY, RMAX, RSQ, A, Q, QX, QY, IER)
        IF (IER /= 0) THEN
          WRITE (*, *) 'QSHEP2_TEST - ERROR!'
          WRITE (*, *) 'Error in QS2GRD, IER = ', IER
          STOP
        END IF
        IF (ABS (Q1 - Q) > 3.0E+00_R8 * ABS (Q) * EPS) THEN
          WRITE (*, *) 'QSHEP2_TEST - Error!'
          WRITE (*, *) 'The interpolated values Q1 (QS2VAL)'
          WRITE (*, *) 'and Q (QS2GRD) differ.'
          WRITE (*, *) 'Q1 = ', Q1
          WRITE (*, *) 'Q  = ', Q
          STOP
        END IF
        EQ = MAX (EQ, ABS (FQ (PX, PY) - Q))
        EQX = MAX (EQX, ABS (FX (PX, PY) - QX))
        EQY = MAX (EQY, ABS (FY (PX, PY) - QY))
      END DO
    END DO
!
! Print the maximum errors and the ratio EQ / EPS.
!
    RQ = EQ / EPS
    WRITE (*, *) 'Maximum absolute errors in the interpolant Q and partial'
    WRITE (*, *) 'derivatives QX and QY relative to machine precision EPS.'
    WRITE (*, *) ' '
    WRITE (*, *) 'Function   Max error   Max error/EPS'
    WRITE (*, *) ' ' 
    WRITE (*, '(''  Q         '',E9.3,''     '',F4.2)') EQ, RQ
    WRITE (*, '(''  QX        '',E9.3)' ) EQX
    WRITE (*, '(''  QY        '',E9.3)' ) EQY
    WRITE (*, *) ' '
    WRITE (*, *) 'Normal end of execution.'
    STOP
!
! Quadratic test function and partial derivatives.
!
    CONTAINS 
      FUNCTION FQ(XX,YY)
        REAL(KIND=R8) :: FQ, XX, YY
        FQ  = (( XX + 2.0E+00_R8 * YY) / 3.0E+00_R8)**2
        RETURN
      END FUNCTION FQ
     
      FUNCTION FX(XX,YY)
        REAL(KIND=R8) :: FX, XX, YY
        FX = 2.0E+00_R8 * (XX + 2.0E+00_R8 * YY) / 9.0E+00_R8
        RETURN
      END FUNCTION FX
    
      FUNCTION FY(XX,YY)
        REAL(KIND=R8) :: FY, XX, YY
        FY = 4.0E+00_R8 * (XX + 2.0E+00_R8 * YY) / 9.0E+00_R8
        RETURN
      END FUNCTION FY
      
  END PROGRAM QSHEP2_TEST
