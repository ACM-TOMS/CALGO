      DOUBLE PRECISION FUNCTION DDOT(N,DX,INCX,DY,INCY)
*
*  Purpose
*  =======
*
*  DDOT returns the REAL dot product of the elements of two
*       double precision vectors,
*           i.e., w = sum(x(j)*y(k))
*  where j = 1 + (i-1)*incx and k = 1 + (i-1)*incy for i = 1 to n.
*
*  Uses unrolled loops.
*
*  Parameters
*  ==========
*
*  N    - INTEGER
*         On entry, N specifies the number of elements in the sum.
*         Unchanged on exit.
*  DX   - DOUBLE PRECISION
*         On entry, DX specifies the vector x above.
*         Unchanged on exit.
*  INCX - INTEGER
*         On entry, INCX specifies the increment parameter used to step
*         through the array DX.
*         Unchanged on exit.
*  DY   - DOUBLE PRECISION
*         On entry, DY specifies the vector y above.
*         Unchanged on exit.
*  INCY - INTEGER
*         On entry, INCY specifies the increment parameter used to step
*         through the array DY.
*         Unchanged on exit.
*
*
*  Level 1 Blas routine
*
*  Toms algorithm 539 -- Lawson et al, 1979
*  Fortran 77 version -- Tim Hopkins, 1994
*
*     .. Scalar Arguments ..
      INTEGER          INCX,INCY,N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION DX(*),DY(*)
*     ..
*     .. Local Scalars ..
      DOUBLE PRECISION DTEMP
      INTEGER          I,IX,IY,M,MP1
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC        MOD
*     ..
      DDOT = 0.0d0
      DTEMP = 0.0d0
      IF (N.LE.0) RETURN
      IF (INCX.EQ.1 .AND. INCY.EQ.1) GO TO 20
*
*        code for unequal increments or equal increments
*          not equal to 1
*
      IX = 1
      IY = 1
      IF (INCX.LT.0) IX = (-N+1)*INCX + 1
      IF (INCY.LT.0) IY = (-N+1)*INCY + 1
      DO 10 I = 1,N
          DTEMP = DTEMP + DX(IX)*DY(IY)
          IX = IX + INCX
          IY = IY + INCY
   10 CONTINUE
      DDOT = DTEMP
      RETURN
*
*        code for both increments equal to 1
*
*
*        clean-up loop
*
   20 M = MOD(N,5)
      IF (M.EQ.0) GO TO 40
      DO 30 I = 1,M
          DTEMP = DTEMP + DX(I)*DY(I)
   30 CONTINUE
      IF (N.LT.5) GO TO 60
   40 MP1 = M + 1
      DO 50 I = MP1,N,5
          DTEMP = DTEMP + DX(I)*DY(I) + DX(I+1)*DY(I+1) +
     +            DX(I+2)*DY(I+2) + DX(I+3)*DY(I+3) + DX(I+4)*DY(I+4)
   50 CONTINUE
   60 DDOT = DTEMP
      RETURN

      END
      SUBROUTINE DAXPY(N,DA,DX,INCX,DY,INCY)
*
*  Purpose
*  =======
*
*  DAXPY computes y = a * x + y where x and y are double precision
*        vectors and a is a double precision scalar.
*
*  The elements of the x-vector used are
*
*     1 + (i-1)*incx,   if incx >= 0,
*     1 + (n-i)*|incx|, if incx < 0.
*
*  and similarly for y and incy.
*
*  Uses unrolled loops.
*
*  Parameters
*  ==========
*
*  N    - INTEGER
*         On entry, N specifies the number of elements to be computed.
*         Unchanged on exit.
*  DA   - DOUBLE PRECISION
*         On entry, DA specifies the value of the scalar a above.
*         Unchanged on exit.
*  DX   - DOUBLE PRECISION
*         On entry, DX specifies the vector x above.
*         Unchanged on exit.
*  INCX - INTEGER
*         On entry, INCX specifies the increment parameter used to step
*         through the array DX.
*         Unchanged on exit.
*  DY   - DOUBLE PRECISION
*         On entry, DY specifies the vector y above.
*         On exit, DY contains DA * DX + DY.
*  INCY - INTEGER
*         On entry, INCY specifies the increment parameter used to step
*         through the array DY.
*         Unchanged on exit.
*
*
*  Level 1 Blas routine
*
*  Toms algorithm 539 -- Lawson et al, 1979
*  Fortran 77 version -- Tim Hopkins, 1994
*
*     .. Scalar Arguments ..
      DOUBLE PRECISION DA
      INTEGER          INCX,INCY,N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION DX(*),DY(*)
*     ..
*     .. Local Scalars ..
      INTEGER          I,IX,IY,M,MP1
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC        MOD
*     ..
      IF (N.LE.0) RETURN
      IF (DA.EQ.0.0d0) RETURN
      IF (INCX.EQ.1 .AND. INCY.EQ.1) GO TO 20
*
*        code for unequal increments or equal increments
*          not equal to 1
*
      IX = 1
      IY = 1
      IF (INCX.LT.0) IX = (-N+1)*INCX + 1
      IF (INCY.LT.0) IY = (-N+1)*INCY + 1
      DO 10 I = 1,N
          DY(IY) = DY(IY) + DA*DX(IX)
          IX = IX + INCX
          IY = IY + INCY
   10 CONTINUE
      RETURN
*
*        code for both increments equal to 1
*
*
*        clean-up loop
*
   20 M = MOD(N,4)
      IF (M.EQ.0) GO TO 40
      DO 30 I = 1,M
          DY(I) = DY(I) + DA*DX(I)
   30 CONTINUE
      IF (N.LT.4) RETURN
   40 MP1 = M + 1
      DO 50 I = MP1,N,4
          DY(I) = DY(I) + DA*DX(I)
          DY(I+1) = DY(I+1) + DA*DX(I+1)
          DY(I+2) = DY(I+2) + DA*DX(I+2)
          DY(I+3) = DY(I+3) + DA*DX(I+3)
   50 CONTINUE
      RETURN

      END
      SUBROUTINE DSCAL(N,DA,DX,INCX)
*
*  Purpose
*  =======
*
*  DSCAL computes x = a*x for a double precision scalar a and
*        elements of a double precision vector x.
*
*  The elements of the x-vector used are
*
*     1 + (i-1)*incx,   if incx >= 0.
*
*  Uses unrolled loops.
*
*  Parameters
*  ==========
*
*  N    - INTEGER
*         On entry, N specifies the number of elements to be scaled.
*         Unchanged on exit.
*  DA   - DOUBLE PRECISION
*         On entry, DA specifies the value of the scalar a above.
*         Unchanged on exit.
*  DX   - DOUBLE PRECISION
*         On entry, DX specifies the vector x above.
*         On exit, DX contains the scaled data.
*  INCX - INTEGER
*         On entry, INCX specifies the increment parameter used to step
*         through the array DX.
*         Unchanged on exit.
*
*
*  Level 1 Blas routine
*
*  Toms algorithm 539 -- Lawson et al, 1979
*  Fortran 77 version -- Tim Hopkins, 1994
*
*     .. Scalar Arguments ..
      DOUBLE PRECISION DA
      INTEGER          INCX,N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION DX(*)
*     ..
*     .. Local Scalars ..
      INTEGER          I,M,MP1,NINCX
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC        MOD
*     ..
      IF (N.LE.0 .OR. INCX.LE.0) RETURN
      IF (INCX.EQ.1) GO TO 20
*
*        code for increment not equal to 1
*
      NINCX = N*INCX
      DO 10 I = 1,NINCX,INCX
          DX(I) = DA*DX(I)
   10 CONTINUE
      RETURN
*
*        code for increment equal to 1
*
*
*        clean-up loop
*
   20 M = MOD(N,5)
      IF (M.EQ.0) GO TO 40
      DO 30 I = 1,M
          DX(I) = DA*DX(I)
   30 CONTINUE
      IF (N.LT.5) RETURN
   40 MP1 = M + 1
      DO 50 I = MP1,N,5
          DX(I) = DA*DX(I)
          DX(I+1) = DA*DX(I+1)
          DX(I+2) = DA*DX(I+2)
          DX(I+3) = DA*DX(I+3)
          DX(I+4) = DA*DX(I+4)
   50 CONTINUE
      RETURN

      END
      INTEGER FUNCTION IDAMAX(N,DX,INCX)
*
*  Purpose
*  =======
*
*  IDAMAX determines, for the double precision vector, x the smallest
*         index i such that
*
*     |x_i| = max{ |x_k| }
*
*  where k = 1 + (j-1)*incx, j = 1 to n.
*
*  Parameters
*  ==========
*
*  N    - INTEGER
*         On entry, N specifies the number of elements to be tested.
*         Unchanged on exit.
*  DX   - DOUBLE PRECISION
*         On entry, DX specifies the vector x above.
*         Unchanged on exit.
*  INCX - INTEGER
*         On entry, INCX specifies the increment parameter used to step
*         through the array DX.
*         Unchanged on exit.
*
*
*  Level 1 Blas routine
*
*  Toms algorithm 539 -- Lawson et al, 1979
*  Fortran 77 version -- Tim Hopkins, 1994
*
*     .. Scalar Arguments ..
      INTEGER          INCX,N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION DX(*)
*     ..
*     .. Local Scalars ..
      DOUBLE PRECISION DMAX
      INTEGER          I,IX
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC        DABS
*     ..
      IDAMAX = 0
      IF (N.LT.1 .OR. INCX.LE.0) RETURN
      IDAMAX = 1
      IF (N.EQ.1) RETURN
      IF (INCX.EQ.1) GO TO 30
*
*        code for increment not equal to 1
*
      IX = 1
      DMAX = DABS(DX(1))
      IX = IX + INCX
      DO 20 I = 2,N
          IF (DABS(DX(IX)).LE.DMAX) GO TO 10
          IDAMAX = I
          DMAX = DABS(DX(IX))
   10     IX = IX + INCX
   20 CONTINUE
      RETURN
*
*        code for increment equal to 1
*
   30 DMAX = DABS(DX(1))
      DO 40 I = 2,N
          IF (DABS(DX(I)).LE.DMAX) GO TO 40
          IDAMAX = I
          DMAX = DABS(DX(I))
   40 CONTINUE
      RETURN

      END
      SUBROUTINE DSWAP(N,DX,INCX,DY,INCY)
*
*  Purpose
*  =======
*
*  DSWAP exchanges elements of the double precision vectors x and y.
*
*  The elements of the x-vector used are
*
*     1 + (i-1)*incx,   if incx >= 0,
*     1 + (n-i)*|incx|, if incx < 0.
*
*  and similarly for y and incy.
*
*  Uses unrolled loops.
*
*  Parameters
*  ==========
*
*  N    - INTEGER
*         On entry, N specifies the number of elements to be swapped.
*         Unchanged on exit.
*  DX   - DOUBLE PRECISION
*         On entry, DX specifies the vector x above.
*         On exit, DX contains the data originally in DY.
*  INCX - INTEGER
*         On entry, INCX specifies the increment parameter used to step
*         through the array DX.
*         Unchanged on exit.
*  DY   - DOUBLE PRECISION
*         On entry, DY specifies the vector y above.
*         On exit, DY contains the data originally in DX.
*  INCY - INTEGER
*         On entry, INCY specifies the increment parameter used to step
*         through the array DY.
*         Unchanged on exit.
*
*
*  Level 1 Blas routine
*
*  Toms algorithm 539 -- Lawson et al, 1979
*  Fortran 77 version -- Tim Hopkins, 1994
*
*     .. Scalar Arguments ..
      INTEGER          INCX,INCY,N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION DX(*),DY(*)
*     ..
*     .. Local Scalars ..
      DOUBLE PRECISION DTEMP
      INTEGER          I,IX,IY,M,MP1
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC        MOD
*     ..
      IF (N.LE.0) RETURN
      IF (INCX.EQ.1 .AND. INCY.EQ.1) GO TO 20
*
*       code for unequal increments or equal increments not equal
*         to 1
*
      IX = 1
      IY = 1
      IF (INCX.LT.0) IX = (-N+1)*INCX + 1
      IF (INCY.LT.0) IY = (-N+1)*INCY + 1
      DO 10 I = 1,N
          DTEMP = DX(IX)
          DX(IX) = DY(IY)
          DY(IY) = DTEMP
          IX = IX + INCX
          IY = IY + INCY
   10 CONTINUE
      RETURN
*
*       code for both increments equal to 1
*
*
*       clean-up loop
*
   20 M = MOD(N,3)
      IF (M.EQ.0) GO TO 40
      DO 30 I = 1,M
          DTEMP = DX(I)
          DX(I) = DY(I)
          DY(I) = DTEMP
   30 CONTINUE
      IF (N.LT.3) RETURN
   40 MP1 = M + 1
      DO 50 I = MP1,N,3
          DTEMP = DX(I)
          DX(I) = DY(I)
          DY(I) = DTEMP
          DTEMP = DX(I+1)
          DX(I+1) = DY(I+1)
          DY(I+1) = DTEMP
          DTEMP = DX(I+2)
          DX(I+2) = DY(I+2)
          DY(I+2) = DTEMP
   50 CONTINUE
      RETURN

      END
      SUBROUTINE DROT(N,DX,INCX,DY,INCY,C,S)
*
*  Purpose
*  =======
*
*  DROT applies a plane rotation of the form
*
*         ( x_i )     ( c  s) ( x_i )
*         (     )  =  (     ) (     )
*         ( y_i )     (-s  c) ( y_i )
*
*  where the ith element of the double precision vector x is
*
*     1 + (i-1)*incx,   if incx >= 0,
*     1 + (n-i)*|incx|, if incx < 0.
*
*  and similarly for y and incy.
*
*  Parameters
*  ==========
*
*  N    - INTEGER
*         On entry, N specifies the number of columns.
*         Unchanged on exit.
*  DX   - DOUBLE PRECISION
*         On entry, DX specifies the vector x above.
*         On exit, DX specified the transformed data.
*  INCX - INTEGER
*         On entry, INCX specifies the increment parameter used to step
*         through the array DX.
*         Unchanged on exit.
*  DY   - DOUBLE PRECISION
*         On entry, DY specifies the vector y above.
*         On exit, DY specified the transformed data.
*  INCY - INTEGER
*         On entry, INCY specifies the increment parameter used to step
*         through the array DY.
*         Unchanged on exit.
*  C    - DOUBLE PRECISION
*         On entry, C specifies the value of c above.
*         Unchanged on exit.
*  S    - DOUBLE PRECISION
*         On entry, S specifies the value of s above.
*         Unchanged on exit.
*
*
*  Level 1 Blas routine
*
*  Toms algorithm 539 -- Lawson et al, 1979
*  Fortran 77 version -- Tim Hopkins, 1994
*
*     .. Scalar Arguments ..
      DOUBLE PRECISION C,S
      INTEGER          INCX,INCY,N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION DX(*),DY(*)
*     ..
*     .. Local Scalars ..
      DOUBLE PRECISION DTEMP
      INTEGER          I,IX,IY
*     ..
      IF (N.LE.0) RETURN
      IF (INCX.EQ.1 .AND. INCY.EQ.1) GO TO 20
*
*       code for unequal increments or equal increments not equal
*         to 1
*
      IX = 1
      IY = 1
      IF (INCX.LT.0) IX = (-N+1)*INCX + 1
      IF (INCY.LT.0) IY = (-N+1)*INCY + 1
      DO 10 I = 1,N
          DTEMP = C*DX(IX) + S*DY(IY)
          DY(IY) = C*DY(IY) - S*DX(IX)
          DX(IX) = DTEMP
          IX = IX + INCX
          IY = IY + INCY
   10 CONTINUE
      RETURN
*
*       code for both increments equal to 1
*
   20 DO 30 I = 1,N
          DTEMP = C*DX(I) + S*DY(I)
          DY(I) = C*DY(I) - S*DX(I)
          DX(I) = DTEMP
   30 CONTINUE
      RETURN

      END
      DOUBLE PRECISION FUNCTION DNRM2(N,DX,INCX)
*
*  Purpose
*  =======
*
*  DNRM2 returns the DOUBLE PRECISION value of the Euclidean norm of
*  elements of the double precision vector x. It thus returns
*
*     ||x|| = sqrt(sum(x_i*x_i)) for i = 1 to n
*
*  where the elements of the x-vector used are
*
*     1 + (i-1)*incx
*
*
*  Parameters
*  ==========
*
*  N    - INTEGER
*         On entry, N specifies the number of elements to be used.
*         Unchanged on exit.
*  DX   - DOUBLE PRECISION
*         On entry, DX specifies the vector whose Euclidean norm is
*         to be computed.
*         Unchanged on exit.
*  INCX - INTEGER
*         On entry, INCX specifies the increment parameter used to step
*         through the array DX.
*         Unchanged on exit.
*
*
*  Level 1 Blas routine
*
*  Toms algorithm 539 -- Lawson et al, 1979
*  Fortran 77 version -- Tim Hopkins, 1994
*
*     .. Scalar Arguments ..
      INTEGER          INCX,N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION DX(*)
*     ..
*     .. Local Scalars ..
      DOUBLE PRECISION CUTHI,CUTLO,HITEST,ONE,SUM,XMAX,ZERO
      INTEGER          I,IX,J,NEXT
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC        DABS,DSQRT,FLOAT
*     ..
*     .. Data statements ..
      DATA             ZERO,ONE/0.0d0,1.0d0/
      DATA             CUTLO,CUTHI/8.232d-11,1.304d19/
*
*     four phase method     using two built-in constants that are
*     hopefully applicable to all machines.
*         cutlo = maximum of  dsqrt(u/eps)  over all known machines.
*         cuthi = minimum of  dsqrt(v)      over all known machines.
*     where
*         eps = smallest no. such that eps + 1. .gt. 1.
*         u   = smallest positive no.   (underflow limit)
*         v   = largest  no.            (overflow  limit)
*
*     brief outline of algorithm..
*
*     phase 1    scans zero components.
*     move to phase 2 when a component is nonzero and .le. cutlo
*     move to phase 3 when a component is .gt. cutlo
*     move to phase 4 when a component is .ge. cuthi/m
*     where m = n for x() real and m = 2*n for complex.
*
*     values for cutlo and cuthi..
*     from the environmental parameters listed in the imsl converter
*     document the limiting values are as follows..
*     cutlo, s.p.   u/eps = 2**(-102) for  honeywell.  close seconds are
*                   univac and dec at 2**(-103)
*                   thus cutlo = 2**(-51) = 4.44089e-16
*     cuthi, s.p.   v = 2**127 for univac, honeywell, and dec.
*                   thus cuthi = 2**(63.5) = 1.30438e19
*     cutlo, d.p.   u/eps = 2**(-67) for honeywell and dec.
*                   thus cutlo = 2**(-33.5) = 8.23181d-11
*     cuthi, d.p.   same as s.p.  cuthi = 1.30438d19
*     data cutlo, cuthi / 8.232d-11,  1.304d19 /
*     data cutlo, cuthi / 4.441e-16,  1.304e19 /
*     ..
*
      IF (N.GT.0 .AND. INCX.GT.0) GO TO 10
      DNRM2 = ZERO
      GO TO 140
*
   10 ASSIGN 30 TO NEXT
      SUM = ZERO
      I = 1
      IX = 1
*                                                 begin main loop
   20 GO TO NEXT(30,40,70,80)

   30 IF (DABS(DX(I)).GT.CUTLO) GO TO 110
      ASSIGN 40 TO NEXT
      XMAX = ZERO
*
*                        phase 1.  sum is zero
*
   40 IF (DX(I).EQ.ZERO) GO TO 130
      IF (DABS(DX(I)).GT.CUTLO) GO TO 110
*
*                                prepare for phase 2.
      ASSIGN 70 TO NEXT
      GO TO 60
*
*                                prepare for phase 4.
*
   50 CONTINUE
      IX = J
      ASSIGN 80 TO NEXT
      SUM = (SUM/DX(I))/DX(I)
   60 XMAX = DABS(DX(I))
      GO TO 90
*
*                   phase 2.  sum is small.
*                             scale to avoid destructive underflow.
*
   70 IF (DABS(DX(I)).GT.CUTLO) GO TO 100
*
*                     common code for phases 2 and 4.
*                     in phase 4 sum is large.  scale to avoid overflow.
*
   80 IF (DABS(DX(I)).LE.XMAX) GO TO 90
      SUM = ONE + SUM* (XMAX/DX(I))**2
      XMAX = DABS(DX(I))
      GO TO 130
*
   90 SUM = SUM + (DX(I)/XMAX)**2
      GO TO 130
*
*
*                  prepare for phase 3.
*
  100 SUM = (SUM*XMAX)*XMAX
*
*
*     for real or d.p. set hitest = cuthi/n
*     for complex      set hitest = cuthi/(2*n)
*
  110 HITEST = CUTHI/FLOAT(N)
*
*                   phase 3.  sum is mid-range.  no scaling.
*
      DO 120 J = IX,N
          IF (DABS(DX(I)).GE.HITEST) GO TO 50
          SUM = SUM + DX(I)**2
          I = I + INCX
  120 CONTINUE
      DNRM2 = DSQRT(SUM)
      GO TO 140
*
  130 CONTINUE
      IX = IX + 1
      I = I + INCX
      IF (IX.LE.N) GO TO 20
*
*              end of main loop.
*
*              compute square root and adjust for scaling.
*
      DNRM2 = XMAX*DSQRT(SUM)
  140 CONTINUE
      RETURN

      END
      SUBROUTINE DROTG(DA,DB,C,S)
*
*  Purpose
*  =======
*
*  DROTG constructs a Givens transformation, i.e, it computes
*        c = cos(theta) and s = sin(theta) such that
*
*         ( c  s) ( a )   ( r )
*         (     ) (   ) = (   )
*         (-s  c) ( b )   ( 0 )
*  It also computes the value z defined as
*
*         ( s   if |s| < c or c = 0,
*     z = (
*         ( 1/c if 0 < |c| <= s.
*
*  If the user later wishes to reconstruct c and s from z,
*  it can be done as follows:
*
*      If z = 1    set c = 0 and s = 1,
*      If |z| < 1  set c = sqrt(1-z^2) and s = z,
*      If |z| >= 1 set c = 1/z and s = sqrt(1-c^2).
*
*  Parameters
*  ==========
*
*  DA  - DOUBLE PRECISION
*        On entry, DA is the first element of the 2-vector.
*        On exit, DA is overwritten by the value of r (see above).
*  DB  - DOUBLE PRECISION
*        On entry, DB is the second element of the 2-vector.
*        On exit, DB is overwritten by the value of z (see above).
*  C   - DOUBLE PRECISION
*        On exit, contains the value of c (see above).
*  S   - DOUBLE PRECISION
*        On exit, contains the value of s (see above).
*
*
*  Level 1 Blas routine
*
*  Toms algorithm 539 -- Lawson et al, 1979
*  Fortran 77 version -- Tim Hopkins, 1994
*
*     .. Scalar Arguments ..
      DOUBLE PRECISION C,DA,DB,S
*     ..
*     .. Local Scalars ..
      DOUBLE PRECISION R,ROE,SCALE,Z
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC        DABS,DSIGN,DSQRT
*     ..
      ROE = DB
      IF (DABS(DA).GT.DABS(DB)) ROE = DA
      SCALE = DABS(DA) + DABS(DB)
      IF (SCALE.NE.0.0d0) GO TO 10
      C = 1.0d0
      S = 0.0d0
      R = 0.0d0
      Z = 0.0d0
      GO TO 20

   10 R = SCALE*DSQRT((DA/SCALE)**2+ (DB/SCALE)**2)
      R = DSIGN(1.0d0,ROE)*R
      C = DA/R
      S = DB/R
      Z = 1.0d0
      IF (DABS(DA).GT.DABS(DB)) Z = S
      IF (DABS(DB).GE.DABS(DA) .AND. C.NE.0.0d0) Z = 1.0d0/C
   20 DA = R
      DB = Z
      RETURN

      END
      DOUBLE PRECISION FUNCTION DASUM(N,DX,INCX)
*
*  Purpose
*  =======
*
*  DASUM returns the DOUBLE PRECISION sum of the absolute values
*        of elements of a double precision vector x.
*
*  The elements of the x-vector used are
*
*     1 + (i-1)*incx,   if incx >= 0.
*
*  Uses unrolled loops.
*
*  Parameters
*  ==========
*
*  N    - INTEGER
*         On entry, N specifies the number of elements to be summed.
*         Unchanged on exit.
*  DX   - DOUBLE PRECISION
*         On entry, DX specifies the vector x above.
*         Unchanged on exit.
*  INCX - INTEGER
*         On entry, INCX specifies the increment parameter used to step
*         through the array DX.
*         Unchanged on exit.
*
*
*  Level 1 Blas routine
*
*  Toms algorithm 539 -- Lawson et al, 1979
*  Fortran 77 version -- Tim Hopkins, 1994
*
*     .. Scalar Arguments ..
      INTEGER          INCX,N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION DX(*)
*     ..
*     .. Local Scalars ..
      DOUBLE PRECISION DTEMP
      INTEGER          I,M,MP1,NINCX
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC        DABS,MOD
*     ..
      DASUM = 0.0d0
      DTEMP = 0.0d0
      IF (N.LE.0 .OR. INCX.LE.0) RETURN
      IF (INCX.EQ.1) GO TO 20
*
*        code for increment not equal to 1
*
      NINCX = N*INCX
      DO 10 I = 1,NINCX,INCX
          DTEMP = DTEMP + DABS(DX(I))
   10 CONTINUE
      DASUM = DTEMP
      RETURN
*
*        code for increment equal to 1
*
*
*        clean-up loop
*
   20 M = MOD(N,6)
      IF (M.EQ.0) GO TO 40
      DO 30 I = 1,M
          DTEMP = DTEMP + DABS(DX(I))
   30 CONTINUE
      IF (N.LT.6) GO TO 60
   40 MP1 = M + 1
      DO 50 I = MP1,N,6
          DTEMP = DTEMP + DABS(DX(I)) + DABS(DX(I+1)) + DABS(DX(I+2)) +
     +            DABS(DX(I+3)) + DABS(DX(I+4)) + DABS(DX(I+5))
   50 CONTINUE
   60 DASUM = DTEMP
      RETURN

      END
