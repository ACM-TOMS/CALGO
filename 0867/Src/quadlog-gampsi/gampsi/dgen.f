      SUBROUTINE  dgen(unit, iflag, nrand, fmin, fmax, xmin, xmax,
     X     nmin, nmax, fun)
************************************************************************
*     (Generate double-precision data)
*     Generate lines of logarithmically-distributed random arguments
*     and function values in the file open on unit.
*
*     There are nrand integer values, f, generated in the range
*     fmin..fmax for each power, 2**n, from n = nmin to nmax.
*
*     Function arguments are represented as x = f * 2**n, where f is
*     an integer value exactly representable in a double-precision
*     number on all current arithmetic systems, so that arguments can
*     be reconstructed exactly, avoiding errors from inaccurate
*     decimal<->binary conversion.
*
*     The function arguments are additionally constrained to lie
*     in the interval x = xmin..xmax.
*
*     If iflag is 1, then only +x is generated.  If iflag is -1, then
*     only -x is generated.  If iflag is 2, then both +x and -x are
*     generated.
*
*     [18-Jun-2000]
************************************************************************
*
*     External functions
*
      EXTERNAL            diran,       fun
*
      DOUBLE PRECISION    diran,       fun
*
*     Parameter variables
*
      DOUBLE PRECISION    two
      PARAMETER           (two = 2.0d+00)
*
*     Argument variables
*
      DOUBLE PRECISION    fmax,        fmin,        xmax,        xmin
*
      INTEGER             iflag,       nmax,        nmin,        nrand
      INTEGER             unit
*
*     Local variables
*
      DOUBLE PRECISION    f,           fofx,        twoton,      x
*
      INTEGER             m,           n
*
      IF ((fmin .LT. fmax) .AND. (xmin .LT. xmax)) THEN
           DO 300 n = nmin, nmax
                twoton = two**n
                IF ((fmax*twoton .LT. xmin) .OR.
     X               (fmin*twoton .GT. xmax)) THEN
*
*                   Skip this n value, since the generated interval
*                   would lie outside xmin..xmax.
*
                ELSE
                     DO 200 m = 1,nrand
  100                     f = diran(fmin,fmax)
                          x = f*twoton
                          IF ((xmin .LE. x) .AND. (x .LE. xmax)) THEN
                               IF (iflag .eq. 1) THEN
                                 fofx = fun(x)
                                 WRITE (unit,10000) 2,f,n,0,fofx,0,0
                               ELSE IF (iflag .eq. -1) THEN
                                 fofx = fun(-x)
                                 WRITE (unit,10000) 2,-f,n,0,fofx,0,0
                               ELSE IF (iflag .eq. 2) THEN
                                 fofx = fun(-x)
                                 WRITE (unit,10000) 2,-f,n,0,fofx,0,0
                                 fofx = fun(x)
                                 WRITE (unit,10000) 2,f,n,0,fofx,0,0
                               END IF
                          ELSE
                               GO TO 100
                          END IF
  200                CONTINUE
                END IF
  300      CONTINUE
      END IF
*
10000 FORMAT (I1, 1X, F21.0, 1X, I5, 1X, I1, 1X, 1P, E26.18E3, 1X,
     X     I1, 1X, I1)
      END
