C
C                          TS2TEST
C                          02/20/97
C
C
C   This program tests the scattered data interpolation
C package TSHEP2D by printing the maximum errors associated
C with interpolated values and gradients on a 20 by 20
C uniform grid in the unit square.  The data set consists
C of 100 nodes with data values taken from a cosine series
C function for which the method is exact.  The ratio of
C maximum interpolation error relative to the machine
C precision is also printed.  This should be O(1).  The
C interpolated values from TS2VAL, TS2GRD, and TS2HES are
C compared for agreement.
C
      INTEGER N, NR
      PARAMETER (N=100, NR=6)
C
C Array storage:
C
      INTEGER LCELL(NR,NR), LNEXT(N)
      DOUBLE PRECISION X(N), Y(N), F(N), RW(N), A(10,N),
     .                 P(20)
C
      INTEGER I, IER, J, K, LOUT, NC, NW
      DOUBLE PRECISION C1, C2, C3, CX, CY, CXX, CXY, CYY,
     .                 DX, DY, EC, ECX, ECY, ECXX, ECXY,
     .                 ECYY, EP1, EPS, PI, PX, PY, RC, RMAX,
     .                 SX, SY, XMIN, YMIN, YK, XX, YY
      DOUBLE PRECISION STORE, TS2VAL, FC, FX, FY, FXX, FXY,
     .                 FYY
C
C TSHEP2 parameters and logical unit for output:
C
      DATA NC/17/,  NW/30/,  LOUT/1/
C
C Cosine series test function and partial derivatives:
C
      FC(XX,YY) = COS(2.*PI*XX)*COS(PI*YY)
      FX(XX,YY) = -2.*PI*SIN(2.*PI*XX)*COS(PI*YY)
      FY(XX,YY) = -PI*COS(2.*PI*XX)*SIN(PI*YY)
      FXX(XX,YY) = -4.*PI*PI*COS(2.*PI*XX)*COS(PI*YY)
      FXY(XX,YY) = 2.*PI*PI*SIN(2.*PI*XX)*SIN(PI*YY)
      FYY(XX,YY) = -PI*PI*COS(2.*PI*XX)*COS(PI*YY)
C
      PI = ACOS(-1.D0)
C
C Output file:
C
      OPEN (LOUT,FILE='res2')
C
C Generate a 10 by 10 grid of nodes in the unit square with
C   the natural ordering.
C
      K = 0
      DO 2 J = 1,10
        YK = DBLE(10-J)/9.
        DO 1 I = 1,10
          K = K + 1
          X(K) = DBLE(I-1)/9.
          Y(K) = YK
    1     CONTINUE
    2   CONTINUE
C
C Compute the data values.
C
      DO 3 K = 1,N
        F(K) = FC(X(K),Y(K))
    3   CONTINUE
C
C Compute parameters defining the interpolant C.
C
      CALL TSHEP2 (N,X,Y,F,NC,NW,NR, LCELL,LNEXT,XMIN,YMIN,
     .             DX,DY,SX,SY,RMAX,RW,A,IER)
      IF (IER .NE. 0) GO TO 10
C
C Generate a 20 by 20 uniform grid of interpolation points
C   (P(I),P(J)) in the unit square.  The four corners coin-
C   cide with nodes.
C
      DO 4 I = 1,20
        P(I) = DBLE(I-1)/19.
    4   CONTINUE
C
C Compute the machine precision EPS.
C
      EPS = 1.
    5 EPS = EPS/2.
      EP1 = EPS + 1.
      IF (STORE(EP1) .GT. 1.) GO TO 5
      EPS = EPS*2.
C
C Compute interpolation errors and test for agreement in the
C   C values returned by TS2VAL, TS2GRD, and TS2HES.
C
      EC = 0.
      ECX = 0.
      ECY = 0.
      ECXX = 0.
      ECXY = 0.
      ECYY = 0.
      DO 7 J = 1,20
        PY = P(J)
        DO 6 I = 1,20
          PX = P(I)
          C1 = TS2VAL (PX,PY,N,X,Y,F,NR,LCELL,LNEXT,XMIN,
     .                 YMIN,DX,DY,SX,SY,RMAX,RW,A)
          CALL TS2GRD (PX,PY,N,X,Y,F,NR,LCELL,LNEXT,XMIN,
     .                 YMIN,DX,DY,SX,SY,RMAX,RW,A, C2,CX,
     .                 CY,IER)
          CALL TS2HES (PX,PY,N,X,Y,F,NR,LCELL,LNEXT,XMIN,
     .                 YMIN,DX,DY,SX,SY,RMAX,RW,A, C3,CX,CY,
     .                 CXX,CXY,CYY,IER)
          IF (IER .NE. 0) GO TO 11
          IF (ABS(C1-C2) .GT. 5.*ABS(C1)*EPS  .OR.
     .        ABS(C1-C3) .GT. 5.*ABS(C1)*EPS) THEN
C
C Values returned by TS2VAL, TS2GRD, and TS2HES differ by a
C   relative amount greater than 3*EPS.
C
            WRITE (LOUT,90) C1, C2, C3
   90       FORMAT (///1X,'*** Error -- interpolated ',
     .              'values C1 (TS2VAL), C2 (TS2GRD), and'/
     .              5X,'C3 (TS2HES) differ:'//
     .              5X,'C1 = ',D21.14/5X,'C2 = ',D21.14/
     .              5X,'C3 = ',D21.14)
          ENDIF
          EC = MAX(EC,ABS(FC(PX,PY)-C1))
          ECX = MAX(ECX,ABS(FX(PX,PY)-CX))
          ECY = MAX(ECY,ABS(FY(PX,PY)-CY))
          ECXX = MAX(ECXX,ABS(FXX(PX,PY)-CXX))
          ECXY = MAX(ECXY,ABS(FXY(PX,PY)-CXY))
          ECYY = MAX(ECYY,ABS(FYY(PX,PY)-CYY))
    6     CONTINUE
    7   CONTINUE
C
C Print errors and the ratio EC/EPS.
C
      RC = EC/EPS
      WRITE (LOUT,100)
      WRITE (LOUT,110) EC, RC
      WRITE (LOUT,120) ECX
      WRITE (LOUT,130) ECY
      WRITE (LOUT,140) ECXX
      WRITE (LOUT,150) ECXY
      WRITE (LOUT,160) ECYY
      STOP
  100 FORMAT (///1X,'Maximum absolute errors in the ',
     .        'interpolant C and partial'/
     .        1X,'derivatives CX, ..., CYY relative ',
     .        'to machine precision EPS'//
     .        1X,10X,'Function',3X,'Max Error',3X,
     .        'Max Error/EPS'/)
  110 FORMAT (1X,13X,'C',7X,E9.3,7X,F5.2)
  120 FORMAT (1X,13X,'CX',6X,E9.3)
  130 FORMAT (1X,13X,'CY',6X,E9.3)
  140 FORMAT (1X,13X,'CXX',5X,E9.3)
  150 FORMAT (1X,13X,'CXY',5X,E9.3)
  160 FORMAT (1X,13X,'CYY',5X,E9.3)
C
C Error in TSHEP2.
C
   10 WRITE (LOUT,200) IER
      STOP
  200 FORMAT (///1X,'*** Error in TSHEP2 -- IER =',I2,
     .        ' ***')
C
C Error in TS2GRD.
C
   11 WRITE (LOUT,210) IER
      STOP
  210 FORMAT (///1X,'*** Error in TS2GRD -- IER =',I2,
     .        ' ***')
      END
      DOUBLE PRECISION FUNCTION STORE (X)
      DOUBLE PRECISION X
C
C***********************************************************
C
C                                               From TRIPACK
C                                            Robert J. Renka
C                                  Dept. of Computer Science
C                                       Univ. of North Texas
C                                           renka@cs.unt.edu
C                                                   03/18/90
C
C   This function forces its argument X to be stored in a
C memory location, thus providing a means of determining
C floating point number characteristics (such as the machine
C precision) when it is necessary to avoid computation in
C high precision registers.
C
C
C On input:
C
C       X = Value to be stored.
C
C X is not altered by this function.
C
C On output:
C
C       STORE = Value of X after it has been stored and
C               possibly truncated or rounded to the single
C               precision word length.
C
C Modules required by STORE:  None
C
C***********************************************************
C
      DOUBLE PRECISION Y
      COMMON/STCOM/Y
C
      Y = X
      STORE = Y
      RETURN
      END

