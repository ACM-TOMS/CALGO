      DOUBLE PRECISION FUNCTION dpsiln(x)
************************************************************************
*     (Double-precision psi(x) - ln(x))
*     Return the value of psi(x) - ln(x), computed so as to avoid
*     unnecessary subtraction loss.
*     [01-Aug-2000]
************************************************************************
*
*     Intrinsic functions
*
      INTRINSIC           dint,        dlog,        idint
*
*     Built-in functions
*
      INTEGER             idint
*
      DOUBLE PRECISION    dint,        dlog
*
*     External functions
*
      EXTERNAL            dfloat,      dnan,        dpsi,        idceil
      EXTERNAL            isdnan
*
      DOUBLE PRECISION    dfloat,      dnan,        dpsi
*
      INTEGER             idceil
*
      LOGICAL             isdnan
*
*     Parameter variables
*
*     See dpsi.f for a discussion of the setting of ncuse and cutoff.
*
      INTEGER             ncuse
      PARAMETER           (ncuse = 27)
*
      DOUBLE PRECISION    cutoff
      PARAMETER           (cutoff = 13.0d+00)
*
      DOUBLE PRECISION    half
      PARAMETER           (half = 0.5d+00)
*
      DOUBLE PRECISION    one
      PARAMETER           (one = 1.0d+00)
*
      DOUBLE PRECISION    zero
      PARAMETER           (zero = 0.0d+00)
*
*     For x values above x1, psi(x) - ln(x) suffers bit loss:
*
*     Maple V5.1:
*     Digits := 60;
*     x1 := fsolve(Psi(x)/ln(x) = 1/2, x = 1.5 .. 2);
*         1.81953794823878564441886676084334572004393270618128250146815
*
      DOUBLE PRECISION                x1
      PARAMETER (x1 =
     X1.81953794823878564441886676084334572004393270618128250146815d+00)
*
*     Argument variables
*
      DOUBLE PRECISION    x
*
*     Local variables
*
      DOUBLE PRECISION    sum,         xfrac,       xsqinv,      y
      DOUBLE PRECISION    ysqinv
*
      INTEGER             i,           n
*
      INCLUDE 'dcpsi.inc'
*
*     Computation of psi(x) - ln(x) for x >= 0 is handled by three
*     cases:
*
*     (a) 0 < x < x1:
*         Return
*             psi(x) - alog(x)  [no bit loss possible]
*
*     (b) x1 <= x < cutoff:
*         Use downward recursion from
*                 psi(x-1) = psi(x) - 1/(x-1)
*                 psiln(x-1) = psiln(x) - ln((x-1)/x) - 1/(x-1)
*                 ...
*                 psiln(x-n) = psiln(x) - ln((x-n)/x) -
*                              sum(i=1:n)(1/(x-i))
*              or
*                 psiln(x) = psiln(x+n) - ln(x/(x+n)) -
*                              sum(i=1:n)(1/(x+n-i))
*         starting in the asymptotic series region.
*
*         In the last formula, psiln(x) and psiln(x+n) are negative,
*         ln(x/(x+n)) is negative, and the sum is positive, so there is
*         subtraction loss from the ln() term and from the x+n-i terms.
*         Indeed, up to 3 bits are lost in forming
*
*             ln(x/(x+n)) + sum(i=1:n)(1/(x+n-i)).
*
*         We should therefore compute that part in higher precision
*         where feasible.
*
*         The straightforward approach, forming psi(x) - ln(x) directly,
*         has a worst-case bit loss at x = cutoff.  For cutoff = 10, 6
*         bits are lost.
*
*     (c) cutoff <= x <= Infinity
*             asymptotic series
*
*     We should do the computation in quadruple precision to reduce
*     the impact of the bit loss in case (b), but that is not
*     universally available.  The result is that about one decimal
*     digit is lost in the range [x1,cutoff).
*
      IF (isdnan(x)) THEN
          dpsiln = dnan()
      ELSE IF (x .LE. zero) THEN
          dpsiln = dnan()
      ELSE IF (x .LT. x1) THEN
          dpsiln = dpsi(x) - dlog(x)
      ELSE IF (x .LT. cutoff) THEN
          xfrac = x - dint(x)
          i = idceil(cutoff)
          n = i - idint(x)
          y = dfloat(i) + xfrac
*
*         We now have y = x + n: form sum = psiln(y) = psiln(x+n)
*
          ysqinv = one/(y * y)
          sum = c(ncuse)
          DO 10 i  = (ncuse - 1), 1, -1
              sum = sum * ysqinv + c(i)
   10     CONTINUE
          sum = -half/y - sum*ysqinv
*
*         Now subtract ln(x/(x+n)) + sum(i=1:n)(1/(x+n-i))
*
          DO 20 i = 1,n
              sum = sum - one / (y - dfloat(i))
   20     CONTINUE
          dpsiln = sum - dlog(x/y)
      ELSE
          xsqinv = one/(x * x)
          sum = c(ncuse)
          DO 30 i  = (ncuse - 1), 1, -1
              sum = sum * xsqinv + c(i)
   30     CONTINUE
          dpsiln = -half/x - sum*xsqinv
      END IF
*
      END
