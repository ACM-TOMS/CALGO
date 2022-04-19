      DOUBLE PRECISION FUNCTION       dpsi(x)
************************************************************************
*     (Double-precision psi(x))
*     Compute and return the value of the psi(x) function for
*     double-precision x.
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
*     but has been augmented for support of IEEE 754 arithmetic.
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
*     Additional code has been added to handle the cases of x = NaN, for
*     which psi(NaN) = NaN, and large x, and x = 0, -1, -2, -3, ...,
*     for which psi(x) = +Infinity, a special value in IEEE 754
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
*     (03-Aug-2000)
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
      EXTERNAL            dfloat,      dinf,        isdnan,      dnan
*
      DOUBLE PRECISION    dabs,        dfloat,      dinf,        dlog
      DOUBLE PRECISION    dnan,        dtan
*
      INTEGER             int
*
      LOGICAL             isdnan
*
*     Parameter variables
*
      DOUBLE PRECISION    four
      PARAMETER           (four = 4.0d0)
*
      DOUBLE PRECISION    fourth
      PARAMETER           (fourth = 0.25d0)
*
      DOUBLE PRECISION    half
      PARAMETER           (half = 0.5d0)
*
      DOUBLE PRECISION    one
      PARAMETER           (one = 1.0d0)
*
*     piov4 = pi / 4
*
      DOUBLE PRECISION    piov4
      PARAMETER           (piov4 =
     X    7.853981633974483096156608458198757210492923498437764552d-01)
*
      DOUBLE PRECISION    three
      PARAMETER           (three = 3.0d0)
*
*     zero of psi(x): x = (x01/x01d) + x02 = 1.4616321449683623412626...
*
      DOUBLE PRECISION    x01
      PARAMETER           (x01 = 187.0d0)
*
      DOUBLE PRECISION    x01d
      PARAMETER           (x01d = 128.0d0)
*
      DOUBLE PRECISION    x02
      PARAMETER           (x02 = 6.9464496836234126266d-04)
*
      INCLUDE 'dxinf.inc'
*
*     Maple V5.1:
*         Digits := 25;
*         printf("%40.25e", fsolve(x*ln(x) = 2^53, x))
*             2.7102974691294103366505040e+14
*
      DOUBLE PRECISION    xlarge
      PARAMETER           (xlarge = 2.71d+14)
*
*     Maple V5.1:
*         Digits := 25;
*         printf("%40.25e", ^52);
*             4.5035996273704960000000000e+15
*
      DOUBLE PRECISION    xmax1
      PARAMETER           (xmax1 = 4.5035996273704960000000000d+15)
*
      INCLUDE 'dxmin.inc'
*
*     The reciprocal, 1/xmin, is representable, since it is less than
*     xinf.
*
      DOUBLE PRECISION    xmin1
      PARAMETER           (xmin1 = xmin)
*
      DOUBLE PRECISION    xsmall
      PARAMETER           (xsmall = 5.80d-09)
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
      DOUBLE PRECISION    aug,         den,         p1(9),       p2(7)
      DOUBLE PRECISION    q1(8),       q2(6),       sgn,         upper
      DOUBLE PRECISION    w,           xx,          z
*
      INTEGER             i,           n,           nq
*
*     Coefficients for approximation to  psi(x)/(x-x0)  over [0.5, 3.0]
*
      DATA p1/4.5104681245762934160d-03, 5.4932855833000385356d+00,
     X    3.7646693175929276856d+02, 7.9525490849151998065d+03,
     X    7.1451595818951933210d+04, 3.0655976301987365674d+05,
     X    6.3606997788964458797d+05, 5.8041312783537569993d+05,
     X    1.6585695029761022321d+05/
      DATA q1/9.6141654774222358525d+01, 2.6287715790581193330d+03,
     X    2.9862497022250277920d+04, 1.6206566091533671639d+05,
     X    4.3487880712768329037d+05, 5.4256384537269993733d+05,
     X    2.4242185002017985252d+05, 6.4155223783576225996d-08/
*
*     Coefficients for approximation to  psi(x) - ln(x) + 1/(2x)
*      for  x > 3.0
*
      DATA p2/ - 2.7103228277757834192d+00, - 1.5166271776896121383d+01,
     X    - 1.9784554148719218667d+01, - 8.8100958828312219821d+00, -
     X    1.4479614616899842986d+00, - 7.3689600332394549911d-02, -
     X    6.5135387732718171306d-21/
      DATA q2/4.4992760373789365846d+01, 2.0240955312679931159d+02,
     X    2.4736979003315290057d+02, 1.0742543875702278326d+02,
     X    1.7463965060678569906d+01, 8.8427520398873480342d-01/
*-----------------------------------------------------------------------
      xx = x
      w = dabs(xx)
      aug = zero
*-----------------------------------------------------------------------
*     Check for valid arguments, then branch to appropriate algorithm
*-----------------------------------------------------------------------
      IF (isdnan(xx)) THEN
*-----------------------------------------------------------------------
*          Argument is a NaN, so psi(NaN) is too.
*-----------------------------------------------------------------------
          dpsi = dnan()
          GO TO 600
      ELSE IF (xx .LT. -xmax1) THEN
          dpsi = dnan()
          GO TO 600
      ELSE IF (w .LT. xmin1) THEN
          IF (xx .LT. zero) THEN
              dpsi = dinf()
          ELSE
              dpsi = -dinf()
          END IF
          GO TO 600
      ELSE IF (xx .GE. half) THEN
          GO TO 200
*-----------------------------------------------------------------------
*         xx < 0.5, use reflection formula: psi(1-xx) = psi(xx) + pi *
*         cot(pi*xx).  Use 1/xx for pi*cotan(pi*xx) when xmin1 < |xx| <=
*         xsmall.
*-----------------------------------------------------------------------
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
*         This assignment should never be reached.  It was present
*         in the original ACM Algorithm 715 code.
*-----------------------------------------------------------------------
          sgn = -piov4
      END IF
      w = w - aint(w)
      nq = int(w*four)
      w = four*(w - dfloat(nq)*fourth)
*-----------------------------------------------------------------------
*     w is now related to the fractional part of 4.0 * xx.  Adjust
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
*         Check for singularity.  The infinities at negative integer
*         arguments are positive and negative: we arbitrarily choose the
*         positive one.
*-----------------------------------------------------------------------
          IF (z .EQ. zero) THEN
              dpsi = dinf()
              GO TO 600
          END IF
          aug = sgn*(four/dtan(z))
      ELSE
          aug = sgn*(four*dtan(z))
      END IF
*-----------------------------------------------------------------------
*     We arrive here when xx < 0, or |xx| <= xsmall.  The relevant
*     reduction formulas are
*         psi(xx) = psi(1-xx) - 1/xx
*         psi(-xx) = psi(xx+1) + pi/tan(pi*xx)
*-----------------------------------------------------------------------
  100 xx = one - xx
  200 IF (xx .GT. three) GO TO 400
*-----------------------------------------------------------------------
*     0.5 <= xx <= 3.0
*-----------------------------------------------------------------------
      den = xx
      upper = p1(1)*xx
      DO 300 i = 1, 7
          den = (den + q1(i))*xx
          upper = (upper + p1(i + 1))*xx
  300 CONTINUE
      den = (upper + p1(9))/(den + q1(8))
      xx = (xx - x01/x01d) - x02
      dpsi = den*xx + aug
      GO TO 600
*-----------------------------------------------------------------------
*     3.0 < xx
*-----------------------------------------------------------------------
  400 IF (xx .LT. xlarge) THEN
          w = one/(xx*xx)
          den = w
          upper = p2(1)*w
          DO 500 i = 1, 5
              den = (den + q2(i))*w
              upper = (upper + p2(i + 1))*w
  500     CONTINUE
          aug = (upper + p2(7))/(den + q2(6)) - half/xx + aug
      END IF
      dpsi = aug + dlog(xx)
*
  600 RETURN
      END
