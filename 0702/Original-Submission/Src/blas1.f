C***********************************************************************
C***********************************************************************
C     BLAS LEVEL 1, DOUBLE PRECISION
C     FROM NETLIB, FRI OCT 12 17:03:00 EDT 1990
C***********************************************************
      SUBROUTINE DAXPY(N,DA,DX,INCX,DY,INCY)
C --------------------------------------------------------
C  CONSTANT TIMES A VECTOR PLUS A VECTOR.
C  USES UNROLLED LOOPS FOR INCREMENTS EQUAL TO ONE.
C  JACK DONGARRA, LINPACK, 3/11/78.
C --------------------------------------------------------
      DOUBLE PRECISION DX(*),DY(*),DA
      INTEGER I,INCX,INCY,IX,IY,M,MP1,N

      IF (N. LE. 0) RETURN
      IF (DA .EQ. 0.0D0) RETURN
      IF (INCX.EQ.1 .AND. INCY.EQ.1) GO TO 20
C --------------------------------------------------------
C  CODE FOR UNEQUAL INCREMENTS OR EQUAL INCREMENTS
C  NOT EQUAL TO 1
C --------------------------------------------------------
      IX = 1
      IY = 1
      IF (INCX .LT. 0) IX = (-N+1)*INCX + 1
      IF (INCY .LT. 0) IY = (-N+1)*INCY + 1
      DO 10 I = 1,N
        DY(IY) = DY(IY) + DA*DX(IX)
        IX = IX + INCX
        IY = IY + INCY
   10 CONTINUE
      RETURN
C --------------------------------------------------------
C  CODE FOR BOTH INCREMENTS EQUAL TO 1
C  CLEAN-UP LOOP
C --------------------------------------------------------
   20 M = MOD(N,4)
      IF (M .EQ. 0) GO TO 40
      DO 30 I = 1,M
        DY(I) = DY(I) + DA*DX(I)
   30 CONTINUE
      IF (N .LT. 4) RETURN
   40 MP1 = M + 1
      DO 50 I = MP1,N,4
        DY(I) = DY(I) + DA*DX(I)
        DY(I + 1) = DY(I + 1) + DA*DX(I + 1)
        DY(I + 2) = DY(I + 2) + DA*DX(I + 2)
        DY(I + 3) = DY(I + 3) + DA*DX(I + 3)
   50 CONTINUE
      RETURN
      END
C***********************************************************************
C***********************************************************************
      SUBROUTINE  DCOPY(N,DX,INCX,DY,INCY)
C --------------------------------------------------------
C  COPIES A VECTOR, X, TO A VECTOR, Y.
C  USES UNROLLED LOOPS FOR INCREMENTS EQUAL TO ONE.
C  JACK DONGARRA, LINPACK, 3/11/78.
C --------------------------------------------------------
      DOUBLE PRECISION DX(*),DY(*)
      INTEGER I,INCX,INCY,IX,IY,M,MP1,N

      IF (N .LE. 0) RETURN
      IF (INCX.EQ.1 .AND. INCY.EQ.1) GO TO 20
C --------------------------------------------------------
C  CODE FOR UNEQUAL INCREMENTS OR EQUAL INCREMENTS
C  NOT EQUAL TO 1
C --------------------------------------------------------
      IX = 1
      IY = 1
      IF (INCX .LT. 0) IX = (-N+1)*INCX + 1
      IF (INCY .LT. 0) IY = (-N+1)*INCY + 1
      DO 10 I = 1,N
        DY(IY) = DX(IX)
        IX = IX + INCX
        IY = IY + INCY
   10 CONTINUE
      RETURN
C --------------------------------------------------------
C  CODE FOR BOTH INCREMENTS EQUAL TO 1
C  CLEAN-UP LOOP
C --------------------------------------------------------
   20 M = MOD(N,7)
      IF (M .EQ. 0) GO TO 40
      DO 30 I = 1,M
        DY(I) = DX(I)
   30 CONTINUE
      IF (N .LT. 7) RETURN
   40 MP1 = M + 1
      DO 50 I = MP1,N,7
        DY(I) = DX(I)
        DY(I + 1) = DX(I + 1)
        DY(I + 2) = DX(I + 2)
        DY(I + 3) = DX(I + 3)
        DY(I + 4) = DX(I + 4)
        DY(I + 5) = DX(I + 5)
        DY(I + 6) = DX(I + 6)
   50 CONTINUE
      RETURN
      END
C***********************************************************************
C***********************************************************************
      DOUBLE PRECISION FUNCTION DDOT(N,DX,INCX,DY,INCY)
C --------------------------------------------------------
C  FORMS THE DOT PRODUCT OF TWO VECTORS.
C  USES UNROLLED LOOPS FOR INCREMENTS EQUAL TO ONE.
C  JACK DONGARRA, LINPACK, 3/11/78.
C --------------------------------------------------------
      DOUBLE PRECISION DX(*),DY(*),DTEMP
      INTEGER I,INCX,INCY,IX,IY,M,MP1,N

      DDOT = 0.0D0
      DTEMP = 0.0D0
      IF (N .LE. 0) RETURN
      IF (INCX.EQ.1. AND. INCY.EQ.1) GO TO 20
C --------------------------------------------------------
C  CODE FOR UNEQUAL INCREMENTS OR EQUAL INCREMENTS
C  NOT EQUAL TO 1
C --------------------------------------------------------
      IX = 1
      IY = 1
      IF (INCX .LT. 0) IX = (-N+1)*INCX + 1
      IF (INCY .LT. 0) IY = (-N+1)*INCY + 1
      DO 10 I = 1,N
        DTEMP = DTEMP + DX(IX)*DY(IY)
        IX = IX + INCX
        IY = IY + INCY
   10 CONTINUE
      DDOT = DTEMP
      RETURN
C --------------------------------------------------------
C   CODE FOR BOTH INCREMENTS EQUAL TO 1
C   CLEAN-UP LOOP
C --------------------------------------------------------
   20 M = MOD(N,5)
      IF (M .EQ. 0) GO TO 40
      DO 30 I = 1,M
        DTEMP = DTEMP + DX(I)*DY(I)
   30 CONTINUE
      IF (N .LT. 5) GO TO 60
   40 MP1 = M + 1
      DO 50 I = MP1,N,5
        DTEMP = DTEMP + DX(I)*DY(I) + DX(I + 1)*DY(I + 1) +
     *   DX(I + 2)*DY(I + 2) + DX(I + 3)*DY(I + 3) + DX(I + 4)*DY(I + 4)
   50 CONTINUE
   60 DDOT = DTEMP
      RETURN
      END
C***********************************************************************
C***********************************************************************
      DOUBLE PRECISION FUNCTION DNRM2 (N,DX,INCX)
      INTEGER NEXT
      DOUBLE PRECISION DX(*),CUTLO,CUTHI,HITEST,SUM,XMAX,ZERO,ONE
      DATA ZERO,ONE /0.0D0,1.0D0/
C --------------------------------------------------------
C  EUCLIDEAN NORM OF THE N-VECTOR STORED IN DX() WITH STORAGE
C  INCREMENT INCX .
C  IF    N .LE. 0 RETURN WITH RESULT = 0.
C  IF N .GE. 1 THEN INCX MUST BE .GE. 1
C        C.L.LAWSON, 1978 JAN 08
C  FOUR PHASE METHOD     USING TWO BUILT-IN CONSTANTS THAT ARE
C  HOPEFULLY APPLICABLE TO ALL MACHINES.
C      CUTLO = MAXIMUM OF  DSQRT(U/EPS)  OVER ALL KNOWN MACHINES.
C      CUTHI = MINIMUM OF  DSQRT(V)      OVER ALL KNOWN MACHINES.
C  WHERE
C      EPS = SMALLEST NO. SUCH THAT EPS + 1. .GT. 1.
C      U   = SMALLEST POSITIVE NO.   (UNDERFLOW LIMIT)
C      V   = LARGEST  NO.            (OVERFLOW  LIMIT)
C
C  BRIEF OUTLINE OF ALGORITHM..
C  PHASE 1    SCANS ZERO COMPONENTS.
C  MOVE TO PHASE 2 WHEN A COMPONENT IS NONZERO AND .LE. CUTLO
C  MOVE TO PHASE 3 WHEN A COMPONENT IS .GT. CUTLO
C  MOVE TO PHASE 4 WHEN A COMPONENT IS .GE. CUTHI/M
C  WHERE M = N FOR X() REAL AND M = 2*N FOR COMPLEX.
C  VALUES FOR CUTLO AND CUTHI..
C  FROM THE ENVIRONMENTAL PARAMETERS LISTED IN THE IMSL CONVERTER
C  DOCUMENT THE LIMITING VALUES ARE AS FOLLOWS..
C  CUTLO, S.P.   U/EPS = 2**(-102) FOR  HONEYWELL.  CLOSE SECONDS ARE
C                UNIVAC AND DEC AT 2**(-103)
C                THUS CUTLO = 2**(-51) = 4.44089E-16
C  CUTHI, S.P.   V = 2**127 FOR UNIVAC, HONEYWELL, AND DEC.
C                THUS CUTHI = 2**(63.5) = 1.30438E19
C  CUTLO, D.P.   U/EPS = 2**(-67) FOR HONEYWELL AND DEC.
C                THUS CUTLO = 2**(-33.5) = 8.23181D-11
C  CUTHI, D.P.   SAME AS S.P.  CUTHI = 1.30438D19
C  DATA CUTLO, CUTHI / 8.232D-11, 1.304D19 /
C  DATA CUTLO, CUTHI / 4.441E-16, 1.304E19 /

      INTEGER N, NN, INCX, i, J
      DATA CUTLO, CUTHI / 8.232D-11, 1.304D19 /

      IF (N .GT. 0) GO TO 10
         DNRM2  = ZERO
         GO TO 300

   10 ASSIGN 30 TO NEXT
      SUM = ZERO
      NN = N * INCX
C ------------------
C BEGIN MAIN LOOP
C ------------------
      I = 1
   20    GO TO NEXT,(30, 50, 70, 110)
   30 IF (DABS(DX(I)) .GT. CUTLO) GO TO 85
      ASSIGN 50 TO NEXT
      XMAX = ZERO
C ------------------
C PHASE 1.  SUM IS ZERO
C ------------------
   50 IF (DX(I) .EQ. ZERO) GO TO 200
      IF (DABS(DX(I)) .GT. CUTLO) GO TO 85
C ------------------
C PREPARE FOR PHASE 2.
C ------------------
      ASSIGN 70 TO NEXT
      GO TO 105
C ------------------
C PREPARE FOR PHASE 4.
C ------------------
  100 I = J
      ASSIGN 110 TO NEXT
      SUM = (SUM / DX(I)) / DX(I)
  105 XMAX = DABS(DX(I))
      GO TO 115
C ------------------
C SUM IS SMALL.
C SCALE TO AVOID DESTRUCTIVE UNDERFLOW.
C ------------------
   70 IF (DABS(DX(I)) .GT. CUTLO) GO TO 75
C ------------------
C COMMON CODE FOR PHASES 2 AND 4.
C IN PHASE 4 SUM IS LARGE.  SCALE TO AVOID OVERFLOW.
C ------------------
  110 IF (DABS(DX(I)) .LE. XMAX) GO TO 115
         SUM = ONE + SUM * (XMAX / DX(I))**2
         XMAX = DABS(DX(I))
         GO TO 200
  115 SUM = SUM + (DX(I)/XMAX)**2
      GO TO 200
C ------------------
C PREPARE FOR PHASE 3.
C ------------------
   75 SUM = (SUM * XMAX) * XMAX
C ------------------
C FOR REAL OR D.P. SET HITEST = CUTHI/N
C FOR COMPLEX      SET HITEST = CUTHI/(2*N)
C ------------------
   85 HITEST = CUTHI/FLOAT(N)
C ------------------
C PHASE 3.  SUM IS MID-RANGE.  NO SCALING.
C ------------------
      DO 95 J =I,NN,INCX
      IF (DABS(DX(J)) .GE. HITEST) GO TO 100
   95    SUM = SUM + DX(J)**2
      DNRM2 = DSQRT(SUM)
      GO TO 300

  200 CONTINUE
      I = I + INCX
      IF (I .LE. NN) GO TO 20
C ------------------
C END OF MAIN LOOP.
C COMPUTE SQUARE ROOT AND ADJUST FOR SCALING.
C ------------------
      DNRM2 = XMAX * DSQRT(SUM)
  300 CONTINUE
      RETURN
      END
