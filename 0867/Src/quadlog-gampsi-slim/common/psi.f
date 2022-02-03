      REAL FUNCTION   psi(x)
************************************************************************
*     (Single-precision psi(x))
*     Compute and return the value of the psi(x) function for
*     single-precision x.
*
*     The psi(x) function is the logarithmic derivative of the
*     Gamma(x) function:
*
*          psi(x) = d/dx (Gamma(x)) / Gamma(x) = d/dx (ln Gamma(x))
*
*     This code correctly handles the case where x is NaN, for which
*     psi(NaN) is a NaN, and the case where x is sufficiently large
*     and positive, or takes one of the integer values 0, -1, -2, -3,
*     ..., for which psi(x) = +Infinity, a special value in IEEE 754
*     arithmetic.
*
*     This code is derived from code given by
*
*         W. J. Cody
*         Algorithm 715: SPECFUN --- A Portable FORTRAN
*         Package of Special Function Routines and Test Drivers
*         ACM Trans. Math. Software 19(1) 22--32, March 1993.
*
*     but has been augmented for support of IEEE 754 arithmetic, and
*     extended with a new rational Pade approximation for psi(x) with x
*     in [1,2], and new algorithms for argument reduction.
*     (03-Aug-2000)
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
*     The code has been prettyprinted to match the conventions of the
*     quadrature package in which it is distributed, and the
*     single-precision code embedded in comments has been discarded.
*
*     Additional code has been added to handle the cases of x = NaN, for
*     which psi(NaN) = NaN, and large x, and x = 0, -1, -2, -3, ...,
*     for which psi(x) = +Infinity, a special value in IEEE 754
*     arithmetic.
*
*     The constants eps, xbig, xminin, and xinf below have been given
*     values suitable for IEEE 754 single-precision arithmetic.  This
*     is applicable on all current Intel, Motorola, and most RISC CPUs,
*     and even on IBM S/390 mainframes since 1999.  Older Compaq/DEC
*     VAX, Convex, Cray, and IBM S/3x0 systems will require changes in
*     those values.
*
************************************************************************
*
*     This function program evaluates the logarithmic derivative of the
*     gamma function,
*
*          psi(x) = d/dx (Gamma(x)) / Gamma(x) = d/dx (ln Gamma(x))
*
*       for real x, where either
*
*              -xmax1 < x < -xmin (x not a negative integer), or
*                xmin < x.
*
*       The calling sequence for this function is
*
*                      Y = PSI(X)
*
*       The main computation uses rational Chebyshev approximations
*       published in Math. Comp. 27, 123-127 (1973) by Cody, Strecok
*       and Thacher.  This transportable program is patterned after
*       the machine-dependent FUNPACK program PSI(X), but cannot match
*       that version for efficiency or accuracy.  This version uses
*       rational approximations that are theoretically accurate to 20
*       significant decimal digits.  The accuracy achieved depends on
*       the arithmetic system, the compiler, the intrinsic functions,
*       and proper selection of the machine-dependent constants.
*
************************************************************************
************************************************************************
*
*     Explanation of machine-dependent constants
*
*       XINF   = largest positive machine number
*       XMAX1  = beta ** (p-1), where beta is the radix for the
*                floating-point system, and p is the number of base-beta
*                digits in the floating-point significand.  This is an
*                upper bound on non-integral floating-point numbers, and
*                the negative of the lower bound on acceptable negative
*                arguments for PSI.  If rounding is necessary, round
*                this value down.
*       XMIN1  = the smallest in magnitude acceptable argument.  We
*                recommend XMIN1 = MAX(1/XINF,xmin) rounded up, where
*                xmin is the smallest positive floating-point number.
*       XSMALL = absolute argument below which  PI*COTAN(PI*X)  may be
*                represented by 1/X.  We recommend XSMALL < sqrt(3
*                eps)/pi, where eps is the smallest positive number such
*                that 1+eps > 1.
*       XLARGE = argument beyond which PSI(X) may be represented by
*                LOG(X).  The solution to the equation
*                   x*ln(x) = beta ** p
*                is a safe value.
*
*      Approximate values for some important machines are
*
*                            beta  p     eps     xmin       XINF
*
*      CDC 7600      (S.P.)    2  48  7.11E-15  3.13E-294  1.26E+322
*      CRAY-1        (S.P.)    2  48  7.11E-15  4.58E-2467 5.45E+2465
*      IEEE (IBM/XT,
*        SUN, etc.)  (S.P.)    2  24  1.19E-07  1.18E-38   3.40E+38
*      IEEE (IBM/XT,
*        SUN, etc.)  (D.P.)    2  53  1.11D-16  2.23E-308  1.79D+308
*      IBM 3033      (D.P.)   16  14  1.11D-16  5.40D-79   7.23D+75
*      SUN 3/160     (D.P.)    2  53  1.11D-16  2.23D-308  1.79D+308
*      VAX 11/780    (S.P.)    2  24  5.96E-08  2.94E-39   1.70E+38
*                    (D.P.)    2  56  1.39D-17  2.94D-39   1.70D+38
*       (G Format)   (D.P.)    2  53  1.11D-16  5.57D-309  8.98D+307
*
*                             XMIN1      XMAX1     XSMALL    XLARGE
*
*      CDC 7600      (S.P.)  3.13E-294  1.40E+14  4.64E-08  9.42E+12
*      CRAY-1        (S.P.)  1.84E-2466 1.40E+14  4.64E-08  9.42E+12
*      IEEE (IBM/XT,
*        SUN, etc.)  (S.P.)  1.18E-38   8.38E+06  1.90E-04  1.20E+06
*      IEEE (IBM/XT,
*        SUN, etc.)  (D.P.)  2.23D-308  4.50D+15  5.80D-09  2.71D+14
*      IBM 3033      (D.P.)  1.39D-76   4.50D+15  5.80D-09  2.05D+15
*      SUN 3/160     (D.P.)  2.23D-308  4.50D+15  5.80D-09  2.71D+14
*      VAX 11/780    (S.P.)  5.89E-39   8.38E+06  1.35E-04  1.20E+06
*                    (D.P.)  5.89D-39   3.60D+16  2.05D-09  2.05D+15
*       (G Format)   (D.P.)  1.12D-308  4.50D+15  5.80D-09  2.71D+14
*
************************************************************************
*
*     Error Returns
*
*      The program returns XINF for X < -XMAX1, for X zero or a negative
*      integer, or when X lies in (-XMIN1, 0), and returns -XINF when X
*      lies in (0, XMIN1).
*
*     Intrinsic functions required are:
*
*         ABS, AINT, DBLE, INT, LOG, REAL, TAN
*
*
*      Author: W. J. Cody
*              Mathematics and Computer Science Division
*              Argonne National Laboratory
*              Argonne, IL 60439
*
*      Latest modification: June 8, 1988
*
************************************************************************
*
*     External functions
*
      EXTERNAL            ainf,        anan,        isanan
*
      INTEGER             int
*
      LOGICAL             isanan
*
      REAL                abs,         ainf,        alog,        anan
      REAL                float,       tan
*
*     Parameter variables
*
      INTEGER             npq
      PARAMETER           (npq = 4)
*
*     The asymptotic series for psi(x) is used for x in (cutoff,INF).
*     The number of terms computed is nc (see above), which determines
*     the size of the array c(*) below, and the number of precomputed
*     initial values stored in it.
*
*     Subject to that limit, the cutoff can be adjusted upwards to
*     reduce the accuracy loss from the use of alog(), and minimize the
*     run time.
*
*     For IEEE 754 quadruple precision, p = 113, and ulp = 2**(p-1) =
*     1.93e-34.
*
*     For IBM RS/6000 and SGI quadruple precision, p = 105, and ulp =
*     2**(p-1) = 2.47e-32.
*
*     The numbers of terms tabulated below were determined by a Maple
*     program that exited the summation loop when the relative size of
*     the last-added term dropped below the ulp value.  Each tabulated
*     cutoff value represents a change in the required number of terms
*     shown in column 2 or column 3.
*
*      =======     ===============
*       cutoff     number of terms
*                  IEEE 754 32-bit
*      =======     ===============
*           10            3
*           13            2
*          346            1
*      =======     ===============
*
      REAL                cutoff
      PARAMETER           (cutoff = 13.0)
*
*     ncuse is the number of terms actually used in the summation of the
*     asymptotic series; it must be set to match cutoff values in the
*     table above, and must also be in the range 1..nc.
*
      INTEGER             ncuse
      PARAMETER           (ncuse = 2)
*
      REAL                four
      PARAMETER           (four = 4.0)
*
      REAL                fourth
      PARAMETER           (fourth = 0.25)
*
      REAL                half
      PARAMETER           (half = 0.5)
*
      REAL                one
      PARAMETER           (one = 1.0)
*
*     piov4 = pi / 4
*
      REAL                piov4
      PARAMETER           (piov4 =
     X    7.853981633974483096156608458198757210492923498437764552e-01)
*
      REAL                two
      PARAMETER           (two = 2.0)
*
*     Argument variables
*
      REAL                x
*
*     Local variables
*
      INTEGER             i,           n,           nq
*
      REAL                aug,         den,         p(npq),      q(npq)
      REAL                sgn,         sum,         upper,       w
      REAL                xx,          z
*
*     zero of psi(x): x = (x01/x01d) + x02 = 1.4616321449683623412626...
*     from Maple:
*         Digits := 65;
*         x0 := fsolve(Psi(x) = 0, x = 1.46 .. 1.47);
*     1.4616321449683623412626595423257213284681962040064463512959884086
*         printf("%.55e\n", x0 - 187/128):
*     6.9464496836234126265954232572132846819620400644635129599e-04
*
      REAL                x01
      PARAMETER           (x01 = 187.0)
*
      REAL                x01d
      PARAMETER           (x01d = 128.0)
*
      REAL                x02
      PARAMETER           (x02 = 6.9464496836234126266e-04)
*
      INCLUDE 'axinf.inc'
*
      INCLUDE 'axlarg.inc'
*
      INCLUDE 'axmax1.inc'
*
      INCLUDE 'axsmal.inc'
*
      INCLUDE 'axmin.inc'
*
*     The reciprocal, 1/xmin, is representable, since it is less than
*     xinf.
*
      REAL                xmin1
      PARAMETER           (xmin1 = xmin)
*
      REAL                zero
      PARAMETER           (zero = 0.0)
*
      INCLUDE 'cpsi.inc'
*
*-----------------------------------------------------------------------
*     Machine-dependent constants
*-----------------------------------------------------------------------
*      DATA xmax1/4.50D+15/, xsmall/1.q-17/, xlarge/1.q+33/
*-----------------------------------------------------------------------
*     Coefficients for approximation to psi(x) over [1.0, 2.0],
*     computed in Maple from psipade.map, the critical kernel of which
*     is
*
*         Digits := 250:
*         x0 := fsolve(Psi(x) = 0, x, 1.46 .. 1.47):
*         the_approx := pade(Psi('x') / ('x' - x0), 'x' = 1.5, [3,3]):
*
*     Maximum relative error in psi(x)/(x - x0) = 6.34e-08
*     at x = 1.0000 for x on [1,2] with Pade degree [3,3]
*
*-----------------------------------------------------------------------
      DATA p /
     X   1.6929944300435845443528597769186657283702956688262609693e-01,
     X   3.6627786162176466409940360877417124471557238249483727573e-01,
     X   1.0143566483084039804223634825017375139043229231586468940e-01,
     X   1.2175984371663803056684280072967881854178802134866685952e-03 /
*
      DATA q /
     X   3.8103879855189839163863629049518839962260307850953867625e-04,
     X   2.4463838698955786398782486219743574621900023320683755352e-01,
     X   2.3386383096200242630862344789634679380182434827064096801e-01,
     X   3.1546003274846176720658166710095083323784374371049467561e-02 /
*-----------------------------------------------------------------------
      xx = x
      w = abs(xx)
      aug = zero
*-----------------------------------------------------------------------
*     Check for valid arguments, then branch to appropriate algorithm
*-----------------------------------------------------------------------
      IF (isanan(xx)) THEN
*-----------------------------------------------------------------------
*          Argument is a NaN, so psi(NaN) is too.
*-----------------------------------------------------------------------
           psi = anan()
           GO TO 600
      ELSE IF (xx .LT. -xmax1) THEN
           psi = anan()
           GO TO 600
      ELSE IF (w .LT. xmin1) THEN
           IF (xx .LT. zero) THEN
               psi = ainf()
           ELSE
               psi = -ainf()
           END IF
           GO TO 600
      ELSE IF (xx .GE. one) THEN
           GO TO 200
*-----------------------------------------------------------------------
*          xx < 0.0, use reflection formula: psi(1-xx) = psi(xx) + pi *
*          cot(pi*xx).  Use 1/xx for pi*cotan(pi*xx) when xmin1 < |xx|
*          <= xsmall.
*-----------------------------------------------------------------------
      ELSE IF (xx .GT. xsmall) THEN
           aug = (-one)/xx
           xx = xx + one
           GO TO 200
      ELSE IF (w .LE. xsmall) THEN
           aug = (-one)/xx
           GO TO 100
      END IF
*-----------------------------------------------------------------------
*     Argument reduction for cot.  We want to compute tan(pi*xx) without
*     unnecessary bit loss.  pi has two leading zero bits in IBM S/360
*     arithmetic, but pi/4 has none, so we form instead
*     tan((pi/4)*(4x)).  Then split |4x| into integer and fractional
*     parts, u and v, and analytically eliminate the integer parts,
*     producing four separate cases indexed by nq, involving tan()
*     arguments of (pi/4)v and (pi/4)(1-v).  The four cases are
*     collapsed into two by setting w = v or w = 1 - v.
*-----------------------------------------------------------------------
      IF (xx .LT. zero) THEN
           sgn = piov4
      ELSE
*-----------------------------------------------------------------------
*          This assignment should never be reached.  It was present
*          in the original ACM Algorithm 715 code.
*-----------------------------------------------------------------------
           sgn = -piov4
      END IF
      w = w - aint(w)
      nq = int(w*four)
      w = four*(w - float(nq)*fourth)
*-----------------------------------------------------------------------
*     w is now related to the fractional part of 4.00 * xx.  Adjust
*     argument to correspond to values in the first quadrant and
*     determine the sign.
*-----------------------------------------------------------------------
      n = nq/2
      IF ((n + n) .NE. nq) w = one - w
      z = piov4*w
      IF (mod(n, 2) .NE. 0) sgn = -sgn
*-----------------------------------------------------------------------
*     Determine the final value for  -pi * cotan(pi*xx)
*-----------------------------------------------------------------------
      n = (nq + 1)/2
      IF (mod(n, 2) .EQ. 0) THEN
*-----------------------------------------------------------------------
*          Check for singularity.  The infinities at negative integer
*          arguments are positive and negative: we arbitrarily chose the
*          positive one.
*-----------------------------------------------------------------------
           IF (z .EQ. zero) THEN
               psi = ainf()
               GO TO 600
           END IF
           aug = sgn*(four/tan(z))
      ELSE
           aug = sgn*(four*tan(z))
      END IF
*-----------------------------------------------------------------------
*     We arrive here when xx < 0, or |xx| <= xsmall.  The relevant
*     reduction formulas are
*         psi(xx) = psi(1-xx) - 1/xx
*         psi(-xx) = psi(xx+1) + pi/tan(pi*xx)
*-----------------------------------------------------------------------
  100 xx = one - xx
      IF (xx .LT. one) THEN
           aug = aug - one/xx
           xx = xx + one
      END IF
  200 IF (xx .LT. cutoff) THEN
*-----------------------------------------------------------------------
*          When xx > 2, reduce to interval [1,2] by
*               psi(xx+1) = psi(xx) + 1/xx
*-----------------------------------------------------------------------
  300      IF (xx .GT. two) THEN
                xx = xx - one
                aug = aug + one/xx
                GO TO 300
           END IF
*-----------------------------------------------------------------------
*          1.0 <= xx <= 2.0
*          Compute the Pade approximation
*              psi(xx) = (sum(k=1,npq)p(k)xx^(k-1)/
*                        sum(k=1,npq)q(k)xx^(k-1))
*                       * (xx - x0) + aug
*          where (xx - x0) is computed in two steps to avoid accuracy
*          loss from IBM S/360 wobbling precision: x0 = 1.4616... has
*          3 zero leading bits in that system.
*-----------------------------------------------------------------------
           den = q(npq)
           upper = p(npq)
           DO 400 i = npq - 1, 1, -1
                den = den*xx + q(i)
                upper = upper*xx + p(i)
  400      CONTINUE
           psi = ((xx - x01/x01d) - x02) * (upper / den) + aug
           GO TO 600
      END IF
*-----------------------------------------------------------------------
*     cutoff <= xx: use asymptotic series, or if xx is at or beyond
*     xlarge, just log(xx).
*-----------------------------------------------------------------------
      IF (xx .LT. xlarge) THEN
           w = one/(xx*xx)
           sum = c(ncuse)
           DO 500 i = (ncuse - 1), 1, -1
                sum = sum*w + c(i)
  500      CONTINUE
           aug = aug - half/xx - sum*w
      END IF
      psi = aug + alog(xx)
*
  600 RETURN
      END
