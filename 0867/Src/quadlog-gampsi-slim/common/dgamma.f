      DOUBLE PRECISION FUNCTION  dgamma(x)
************************************************************************
*     (Double-precision Gamma(x))
*     Compute and return the value of the Gamma(x) function for
*     double-precision x.
*
*     This code correctly handles the case where x is NaN, for which
*     Gamma(NaN) is a NaN. and the case where x is sufficiently large
*     and positive, or takes one of the integer values 0, -1, -2, -3,
*     ..., for which Gamma(x) = +Infinity, a special value in IEEE 754
*     arithmetic.
*
*     This code is substantially similar to that given by
*
*         W. J. Cody
*         Algorithm 715: SPECFUN --- A Portable FORTRAN
*         Package of Special Function Routines and Test Drivers
*         ACM Trans. Math. Software 19(1) 22--32, March 1993.
*
*     but has been augmented for support of IEEE 754 arithmetic.
*     (04-Aug-2000)
************************************************************************
*     FOREWORD: This function is derived from code available on the
*     World-Wide Web at ftp://ftp.netlib.org/toms/715, corresponding to
*     the article
*
*         W. J. Cody
*         Algorithm 715: SPECFUN --- A Portable FORTRAN
*         Package of Special Function Routines and Test Drivers
*         ACM Trans. Math. Software 19(1) 22--32 March 1993.
*
*     Additional code has been added to handle the cases of x = NaN, for
*     which Gamma(NaN) = NaN, and large x, and x = 0, -1, -2, -3, ...,
*     for which Gamma(x) = +Infinity, a special value in IEEE 754
*     arithmetic.
*
*     The code has been prettyprinted to match the conventions of the
*     quadrature package in which it is distributed, and the
*     single-precision code embedded in comments has been discarded.
*
*     The constants eps, xbig, xminin, and xinf below have been given
*     values suitable for IEEE 754 double-precision arithmetic.  This is
*     applicable on all current Intel, Motorola, and RISC CPUs, and even
*     on IBM S/390 mainframes since 1999.  Older Compaq/DEC VAX, Convex,
*     Cray, and IBM S/3x0 systems will require changes in those values.
*
************************************************************************
*
*     This routine calculates the gamma function for a real argument x.
*
*     Computation is based on an algorithm outlined in W. J. Cody, ``An
*     overview of software development for special functions'', Lecture
*     Notes in Mathematics, 506, Numerical Analysis Dundee, 1975,
*     G. A. Watson (Ed.), Springer Verlag, Berlin, 1976.
*
*     The program uses rational functions that approximate the gamma
*     function to at least 20 significant decimal digits.  Coefficients
*     for the approximation over the interval (1,2) are unpublished.
*     Those for the approximation for x .ge. 12 are from Hart, et. al.,
*     Computer Approximations, Wiley and Sons, New York, 1968.  Lower
*     order approximations can be substituted for these on machines with
*     less precise arithmetic.
*
************************************************************************
************************************************************************
*
*     Explanation of machine-dependent constants:
*
*     eps    - the smallest positive floating-point number such that
*              1.0 + eps .gt. 1.0
*     xbig   - the largest argument for which Gamma(x) is
*              representable in the machine, i.e., the solution to the
*              equation
*                      Gamma(xbig) = xinf.
*     xminin - the smallest positive floating-point number such that
*              1/xminin is machine representable.
*     xinf   - the largest machine representable floating-point number.
*
*     Approximate values for some important machines are:
*
*             IBM/195    CDC/7600  UNIVAC/1108   VAX 11/780 (UNIX)
*              (d.p.)  (s.p.,rndg)    (d.p.)     (s.p.)     (d.p.)
*
*     eps     2.221d-16  3.553e-15   1.735d-18   5.961e-08  1.388d-17
*     xbig    57.574     177.802     171.489     34.844     34.844
*     xminin  1.382d-76  3.132e-294  1.113d-308  5.883e-39  5.883d-39
*     xinf    7.23d+75   1.26e+322   8.98d+307   1.70e+38   1.70d+38
*
************************************************************************
************************************************************************
*
*
*     Error returns:
*
*      The program returns the value xinf for singularities or
*      when overflow would occur.  The computation is believed
*      to be free of underflow and overflow.
*
*
*
*     Other subprograms required (single precision version):
*
*         alog,exp,float,ifix,sin
*
*     Other subprograms required (double precision version):
*
*         dble,dexp,dlog,dsin,float,ifix,sngl
*
*
*
*      Author: W. J. Cody
*              Applied Mathematics Division
*              Argonne National Laboratory
*              Argonne, IL 60439
*
*      Latest modification: May 18, 1982
*
************************************************************************
*
*     Intrinsic functions
*
      INTRINSIC           dcos,        dexp,        dint,        dlog
      INTRINSIC           dsin,        idint,       mod
*
*     Built-in functions
*
      DOUBLE PRECISION    dcos,        dexp,        dint,        dlog
      DOUBLE PRECISION    dsin
*
      INTEGER             idint,       mod
*
*     External functions
*
      EXTERNAL            dfloat,      dinf,        dnan,        isdnan
*
      DOUBLE PRECISION    dfloat,      dinf,        dnan
*
      LOGICAL             isdnan
*
*     Parameter variables
*
      DOUBLE PRECISION    four
      PARAMETER           (four = 4.0d+00)
*
      DOUBLE PRECISION    half
      PARAMETER           (half = 0.5d0)
*
      DOUBLE PRECISION    one
      PARAMETER           (one = 1.0d0)
*
*     DOUBLE PRECISION    pi
*     PARAMETER           (pi =
*    X    3.14159265358979323846264338327950288419716939937511d+00)
*
*     piov4 = pi / 4
*
      DOUBLE PRECISION    piov4
      PARAMETER           (piov4 =
     X    7.853981633974483096156608458198757210492923498437764552d-01)
*
      DOUBLE PRECISION    two
      PARAMETER           (two = 2.0d0)
*
      INCLUDE 'dl2pi2.inc'
*
      INCLUDE 'deps.inc'
*
      INCLUDE 'dxinf.inc'
*
      INCLUDE 'dxbig.inc'
*
      DOUBLE PRECISION    cutoff
*     PARAMETER           (cutoff = 12.0d0)
      PARAMETER           (cutoff = xbig + one)
*
      INCLUDE 'dxmin.inc'
*
*     Maple V5.1:
*         Digits := 25;
*         printf("%.25e\n", 2^(-1022)):
*             2.2250738585072013830902330e-308
*         printf("%.25e\n", 2^1022):
*             4.4942328371557897693232630e+307
*
*     Thus, the reciprocal is representable, since it is less than xinf.
*     There are smaller denormalized values whose reciprocals are
*     representable, the smallest being
*
*         printf("%.25e\n", 1/((2 - 2^(-53)) * 2^1023)):
*             5.5626846462680037665166100e-309
*
*     However, some architectures carry severe performance penalties for
*     denormalized values, so we avoid using them.
*
      DOUBLE PRECISION    xminin
      PARAMETER           (xminin = xmin)
*
      DOUBLE PRECISION    zero
      PARAMETER           (zero = 0.0d0)
*
*     Argument variables
*
      DOUBLE PRECISION    x
*
*     Local variables
*
      DOUBLE PRECISION    c(7),        fact,        p(8),        q(8)
      DOUBLE PRECISION    result,      sinval,      sum,         u
      DOUBLE PRECISION    v,           w,           xden,        xnum
      DOUBLE PRECISION    y,           y1,          yfrac,       yint
      DOUBLE PRECISION    ysq,         z
*
      INTEGER             i,           n,           umod8
*
      LOGICAL             parity
*
      DATA p/ -1.71618513886549492533811d+0,
     X    2.47656508055759199108314d+1, -3.79804256470945635097577d+2,
     X    6.29331155312818442661052d+2, 8.66966202790413211295064d+2,
     X    -3.14512729688483675254357d+4, -3.61444134186911729807069d+4,
     X    6.64561438202405440627855d+4/
      DATA q/ -3.08402300119738975254353d+1,
     X    3.15350626979604161529144d+2, -1.01515636749021914166146d+3,
     X    -3.10777167157231109440444d+3, 2.25381184209801510330112d+4,
     X    4.75584627752788110767815d+3, -1.34659959864969306392456d+5,
     X    -1.15132259675553483497211d+5/
      DATA c/ -1.910444077728d-03, 8.4171387781295d-04,
     X    -5.952379913043012d-04, 7.93650793500350248d-04,
     X    -2.777777777777681622553d-03, 8.333333333333333331554247d-02,
     X    5.7083835261d-03/
*
      parity = .FALSE.
      fact = one
      n = 0
      y = x
      IF (isdnan(x)) THEN
*-----------------------------------------------------------------------
*          Argument is a NaN, so Gamma(NaN) is too.
*-----------------------------------------------------------------------
           dgamma = dnan()
           RETURN
      ELSE IF (y .LE. zero) THEN
*-----------------------------------------------------------------------
*         Argument is zero or negative.  Gamma(x) is finite, except at
*         the poles where x = 0, -1, -2, -3, ...
*-----------------------------------------------------------------------
          y = -x
          yint = dint(y)
          yfrac = y - yint
          IF (yfrac .NE. zero) THEN
              parity = (yint .NE. dint(yint*half)*two)
*-----------------------------------------------------------------------
*             We want to set fact = -pi/dsin(pi*yfrac).  However, pi
*             has two leading zero bits in IBM S/360 arithmetic, so
*             compute this more accurately as fact = -4*(pi/4) /
*             sin((pi/4)*(4*yfrac)), and reduce this further by
*             4*yfrac = w = u (integer) + v (fraction), where w, u,
*             and v are all nonnegative and small, producing eight
*             cases involving arguments with v and 1-v (which does not
*             suffer from subtraction loss).
*-----------------------------------------------------------------------
              w = four * yfrac
              u = dint(w)
              v = w - u
              umod8 = mod(idint(u),8)
              IF (umod8 .eq. 0) THEN
                   sinval = dsin(piov4 * v)
              ELSE IF (umod8 .eq. 1) THEN
                   sinval = dcos(piov4 * (one -v))
              ELSE IF (umod8 .eq. 2) THEN
                   sinval = dcos(piov4 * v)
              ELSE IF (umod8 .eq. 3) THEN
                   sinval = dsin(piov4 * (one - v))
              ELSE IF (umod8 .eq. 4) THEN
                   sinval = -dsin(piov4 * v)
              ELSE IF (umod8 .eq. 5) THEN
                   sinval = -dcos(piov4 * (one - v))
              ELSE IF (umod8 .eq. 6) THEN
                   sinval = -dcos(piov4 * v)
              ELSE
                   sinval = -dsin(piov4 * (one - v))
              END IF
              fact = -four * (piov4 / sinval)
*OLD:         fact = (-pi)/dsin(pi*yfrac)
              y = y + one
          ELSE
*-----------------------------------------------------------------------
*             Argument is negative integer, so Gamma(x) = +/-Infinity
*-----------------------------------------------------------------------
              result = dinf()
              GO TO 400
          END IF
      END IF
*-----------------------------------------------------------------------
*     Argument is now positive, with fact and parity recording history
*     of negative argument.
*-----------------------------------------------------------------------
      IF (y .LT. eps) THEN
*-----------------------------------------------------------------------
*         Argument .LT. eps
*-----------------------------------------------------------------------
          IF (y .GE. xminin) THEN
              result = one/y
          ELSE
              result = dinf()
              GO TO 400
          END IF
      ELSE IF (y .LT. cutoff) THEN
          y1 = y
          IF (y .LT. one) THEN
*-----------------------------------------------------------------------
*             0.0 .LT. argument .LT. 1.0
*-----------------------------------------------------------------------
              z = y
              y = y + one
          ELSE
*-----------------------------------------------------------------------
*            1.0 .LT. argument .LT. cutoff, reduce argument if necessary
*-----------------------------------------------------------------------
              n = idint(y) - 1
              y = y - dfloat(n)
              z = y - one
          END IF
*-----------------------------------------------------------------------
*         Evaluate approximation for 1.0 .LT. argument .LT. 2.0
*-----------------------------------------------------------------------
          xnum = zero
          xden = one
          DO 100 i = 1, 8
              xnum = (xnum + p(i))*z
              xden = xden*z + q(i)
  100     CONTINUE
          result = xnum/xden + one
          IF (y1 .LT. y) THEN
*-----------------------------------------------------------------------
*             Adjust result for case  0.0 .LT. argument .LT. 1.0
*-----------------------------------------------------------------------
              result = result/y1
          ELSE IF (y1 .GT. y) THEN
*-----------------------------------------------------------------------
*             Adjust result for case  2.0 .LT. argument .LT. cutoff
*-----------------------------------------------------------------------
              DO 200 i = 1, n
                  result = result*y
                  y = y + one
  200         CONTINUE
          END IF
      ELSE
*-----------------------------------------------------------------------
*         Argument .GE. cutoff
*-----------------------------------------------------------------------
          IF (y .LE. xbig) THEN
*-----------------------------------------------------------------------
*             Use asymptotic series for ln(Gamma(x)).  The four
*             parts are summed in order from smallest to largest,
*             positive parts first.
*-----------------------------------------------------------------------
              ysq = y*y
              sum = c(7)
              DO 300 i = 1, 6
                  sum = sum/ysq + c(i)
  300         CONTINUE
              sum = ((y - half)*dlog(y) +(ln2pi2 + sum/y)) - y
              result = dexp(sum)
          ELSE
              result = dinf()
          END IF
      END IF
*-----------------------------------------------------------------------
*     Final adjustments and return
*-----------------------------------------------------------------------
  400 IF (parity) result = -result
      IF (fact .NE. one) result = fact/result
      dgamma = result
      END
