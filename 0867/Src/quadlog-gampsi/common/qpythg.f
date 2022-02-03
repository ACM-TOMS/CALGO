      REAL*16 FUNCTION qpythg(a,b)
************************************************************************
*     (Quadruple-precision Euclidean 2-norm)
*     Find dsqrt(a**2 + b**2) without overflow or destructive underflow,
*     and handle Infinity and NaN arguments correctly.
*
*     The original version of this function (from EISPACK-2) due to
*
*         Cleve Moler and Donald Morrison, ``Replacing Square Roots by
*         Pythagorean Sums'', IBM J.  Research and Development, 27,
*         577--581 (1983)
*
*         Augustin A. Dubrulle, ``A Class of Numerical Methods for the
*         Computation of Pythagorean Sums'', IBM J.  Research and
*         Development, 27, 582--589 (1983)
*
*     did not correctly handle Infinity and NaN arguments: it went
*     into an infinite loop.  This version is more robust, with two
*     extra tests to detect such arguments.
*
*     (01-May-2000)
************************************************************************
*
      REAL*16   qabs,        qmax1,       qmin1
*
*     Argument variables
*
      REAL*16   a,           b
*
*     Local variables
*
      REAL*16   p,           r,           s,           t
      REAL*16   u
*
      p = qmax1(qabs(a), qabs(b))
      IF (p .EQ. 0.0Q0) GO TO 20
*
*     [01-May-2000] Add test for NaN to prevent infinite loop:
*
      IF (p .NE. p) GO TO 20
*
*     [01-May-2000] Add test for Infinity to prevent infinite loop:
*
      IF (p .EQ. (p + p)) GO TO 20
*
      r = (qmin1(qabs(a), qabs(b))/p)**2
   10 CONTINUE
      t = 4.0Q0 + r
      IF (t .EQ. 4.0Q0) GO TO 20
      s = r/t
      u = 1.0Q0 + 2.0Q0*s
      p = u*p
      r = (s/u)**2*r
      GO TO 10
   20 qpythg = p
C
      END
