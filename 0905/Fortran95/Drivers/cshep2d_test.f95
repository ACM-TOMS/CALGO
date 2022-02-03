  PROGRAM CS2TEST
!
!                          CS2TEST
!                          02/19/97
!
!
!   This program tests the scattered data interpolation
! package CSHEP2D by printing the maximum errors associated
! with interpolated values and gradients on a 20 by 20
! uniform grid in the unit square.  The data set consists
! of 100 nodes with data values taken from a cubic function
! for which the method is exact.  The ratio of maximum
! interpolation error relative to the machine precision is
! also printed.  This should be O(1).  The interpolated
! values from CS2VAL, CS2GRD, and CS2HES are compared for
! agreement.
!
    USE SHEPPACK, ONLY:CSHEP2, CS2VAL, CS2HES, CS2GRD
    USE REAL_PRECISION
    IMPLICIT NONE
    INTEGER :: I, IER, J, K
    INTEGER, PARAMETER :: N = 100, NC = 17, NR = 6, NW = 30
    INTEGER, DIMENSION(N) :: LNEXT
    INTEGER, DIMENSION (NR, NR) :: LCELL
  
    REAL (KIND=R8) :: C1, C2, C3, CX, CY, CXX, CXY, CYY, DX, DY, EC,&
                      ECX, ECY, ECXX, ECXY, ECYY, EPS, PX, PY, RC, RMAX,&
                      XMIN, YMIN, YK
    REAL (KIND=R8), DIMENSION (N) :: F, RW, X, Y
    REAL (KIND=R8), DIMENSION (20) :: P
    REAL (KIND=R8), DIMENSION (9,N) :: A
    
!
! Generate a 10 by 10 grid of nodes in the unit square with
!   the natural ordering.
!
    K = 0
    DO J = 1,10
      YK = (10-J)/9.0_R8
      DO I = 1,10
        K = K + 1
        X(K) = (I-1)/9.0_R8
        Y(K) = YK
      END DO
    END DO
!
! Compute the data values.
!
    DO K = 1, N
      F(K) = FC(X(K),Y(K))
    END DO
!
! Compute parameters defining the interpolant C.
!
    CALL CSHEP2 (N, X, Y, F, NC, NW, NR, LCELL, LNEXT, XMIN, YMIN, DX,&
                 DY,RMAX,RW,A,IER)
    IF (IER .NE. 0) THEN
      WRITE (*, *) 'Error in CSHEP2, IER = '
      WRITE (*, *) IER
      STOP
    END IF
! Generate a 20 by 20 uniform grid of interpolation points
!   (P(I),P(J)) in the unit square.  The four corners coin-
!   cide with nodes.
!
    DO I = 1,20
      P(I) = (I-1)/19.0_R8
    END DO
        
!
! Compute the machine precision EPS.
!
    EPS = EPSILON (EPS)
!
! Compute interpolation errors and test for agreement in the
!   C values returned by CS2VAL, CS2GRD, and CS2HES.
!
    EC = 0.0_R8
    ECX = 0.0_R8
    ECY = 0.0_R8
    ECXX = 0.0_R8
    ECXY = 0.0_R8
    ECYY = 0.0_R8
    DO J = 1,20
      PY = P(J)
      DO I = 1,20
        PX = P(I)
        C1 = CS2VAL (PX, PY, N, X, Y, F, NR, LCELL, LNEXT, XMIN,&
                     YMIN, DX, DY, RMAX, RW, A)
        CALL CS2GRD (PX, PY, N, X, Y, F, NR, LCELL, LNEXT, XMIN,&
                     YMIN, DX, DY, RMAX, RW, A, C2, CX, CY, IER)
        CALL CS2HES (PX, PY, N, X, Y, F ,NR, LCELL, LNEXT, XMIN,&
                     YMIN, DX, DY, RMAX, RW, A, C3, CX, CY, CXX,&
                     CXY, CYY, IER)
        IF (IER .NE. 0) THEN
          WRITE (*, *) 'Error in CSHEP2, IER = ' 
          WRITE (*, *) IER
          STOP
        END IF
        IF (ABS(C1-C2) .GT. 8.0_R8*ABS(C1)*EPS  .OR. &
            ABS(C1-C3) .GT. 8.0_R8*ABS(C1)*EPS) THEN
! Values returned by CS2VAL, CS2GRD, and CS2HES differ by a
!  relative amount greater than 4*EPS.
!
          WRITE (*, *) ' '
          WRITE (*, *) 'CSHEP2 - Error.'
          WRITE (*, *) 'The interpolated values C1 (CS2VAL), C2(CS2GRD) and'
          WRITE (*, *) 'C3(CS2HES) differ.'
          WRITE (*, *) 'C1 = ', C1
          WRITE (*, *) 'C2 = ', C2
          WRITE (*, *) 'C3 = ', C3
          STOP
        END IF           
        EC = MAX(EC,ABS(FC(PX,PY)-C1))
        ECX = MAX(ECX,ABS(FX(PX,PY)-CX))
        ECY = MAX(ECY,ABS(FY(PX,PY)-CY))
        ECXX = MAX(ECXX,ABS(FXX(PX,PY)-CXX))
        ECXY = MAX(ECXY,ABS(FXY(PX,PY)-CXY))
        ECYY = MAX(ECYY,ABS(FYY(PX,PY)-CYY))
      END DO
    END DO
!
! Print errors and the ratio EC/EPS.
!
    RC = EC/EPS
    WRITE (*, *) ' '
    WRITE (*, *) 'Maximum absolute errors in the interpolant C '
    WRITE (*, *) 'and partial derivatives (CXX, ... , CYY) relative '
    WRITE (*, *) 'to machine precision EPS.' 
    WRITE (*, *) ' '
    WRITE (*, *) '  FUNCTION   MAX ERROR  MAX ERROR/EPS'
    WRITE (*, '(''      C       '',E9.3,''       '',F4.2)') EC, RC
    WRITE (*, '(''      CX      '',E9.3)') ECX
    WRITE (*, '(''      CY      '',E9.3)') ECY
    WRITE (*, '(''      CXX     '',E9.3)')ECXX
    WRITE (*, '(''      CXY     '',E9.3)')ECXY
    WRITE (*, '(''      CYY     '',E9.3)')ECYY
    WRITE (*, *) ' '
    WRITE (*, *) 'Normal end of execution.'
    STOP
!
! Cubic test function and partial derivatives:
!

   CONTAINS
      FUNCTION FC(XX, YY)
        USE REAL_PRECISION           
        REAL(KIND=R8):: FC, XX, YY
        FC = ((XX + 2.0_R8 * YY) / 3.0_R8)**3
        RETURN
      END FUNCTION FC
       
      FUNCTION FX(XX, YY)
        USE REAL_PRECISION           
        REAL(KIND=R8):: FX, XX, YY
        FX = ((XX + 2.0_R8 * YY) / 3.0_R8)**2
        RETURN
      END FUNCTION FX
   
      FUNCTION FY(XX, YY)
        USE REAL_PRECISION     
        REAL(KIND=R8):: FY, XX, YY
        FY = 2.0_R8*((XX + 2.0_R8 * YY) / 3.0_R8)**2
        RETURN
      END FUNCTION FY
      
      FUNCTION FXX(XX, YY)
        REAL(KIND=R8):: FXX, XX, YY
        FXX = 2.0_R8*(XX + 2.0_R8 * YY) / 9.0_R8
        RETURN
      END FUNCTION FXX
   
      FUNCTION FXY(XX, YY)
        REAL(KIND=R8):: FXY, XX, YY
        FXY = 4.0_R8*(XX + 2.0_R8 * YY) / 9.0_R8
        RETURN
      END FUNCTION FXY
    
      FUNCTION FYY(XX, YY)
        REAL(KIND=R8):: FYY, XX, YY
        FYY = 8.0_R8*(XX + 2.0_R8 * YY) / 9.0_R8
        RETURN
      END FUNCTION FYY
 
   END PROGRAM CS2TEST
