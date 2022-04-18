      SUBROUTINE  qcopy(n,dx,incx,dy,incy)
************************************************************************
*     (Quadruple-precision vector copy)
*     Copy n elements of a vector dx(*) to a vector dy(*), starting at
*     index 1, and then stepping with increments incx and incy,
*     respectively.
*
*     Use unrolled loops for increments equal to one.
*
*     This is a LINPACK BLAS Level 1 primitive, with revised header
*     comments.
*     (11-Mar-1978)
************************************************************************
*
*     Argument variables
*
      REAL*16    dx(*),       dy(*)
*
      INTEGER             incx,        incy,        n
*
*     Local variables
*
      INTEGER             i,           ix,          iy,          m
      INTEGER             mp1
*
*
*     COPIES A VECTOR, X, TO A VECTOR, Y.
*     USES UNROLLED LOOPS FOR INCREMENTS EQUAL TO ONE.
*     JACK DONGARRA, LINPACK, 3/11/78.
*
*
      IF(N.LE.0)RETURN
      IF(INCX.EQ.1.AND.INCY.EQ.1)GO TO 20
*
*        CODE FOR UNEQUAL INCREMENTS OR EQUAL INCREMENTS
*          NOT EQUAL TO 1
*
      IX = 1
      IY = 1
      IF(INCX.LT.0)IX = (-N+1)*INCX + 1
      IF(INCY.LT.0)IY = (-N+1)*INCY + 1
      DO 10 I = 1,N
        DY(IY) = DX(IX)
        IX = IX + INCX
        IY = IY + INCY
   10 CONTINUE
      RETURN
*
*        CODE FOR BOTH INCREMENTS EQUAL TO 1
*
*
*        CLEAN-UP LOOP
*
   20 M = MOD(N,7)
      IF( M .EQ. 0 ) GO TO 40
      DO 30 I = 1,M
        DY(I) = DX(I)
   30 CONTINUE
      IF( N .LT. 7 ) RETURN
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
