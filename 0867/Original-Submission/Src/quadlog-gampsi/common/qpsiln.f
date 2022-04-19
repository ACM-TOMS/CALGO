      REAL*16 FUNCTION qpsiln(x)
************************************************************************
*     (Quadruple-precision psi(x) - ln(x))
*     Return the value of psi(x) - ln(x), computed so as to avoid
*     unnecessary subtraction loss.
*     [01-Aug-2000]
************************************************************************
*
*     Intrinsic functions
*
      INTRINSIC           iqint,      qint,         qlog
*
*     Built-in functions
*
      INTEGER             iqint
*
      REAL*16             qint,        qlog
*
*     External functions
*
      EXTERNAL            iqceil,      isqnan,      qfloat,      qnan
      EXTERNAL            qpsi
*
      INTEGER             iqceil
*
      LOGICAL             isqnan
*
      REAL*16             qfloat,      qnan,        qpsi
*
*     Parameter variables
*
*     See qpsi.f for a discussion of the setting of ncuse and cutoff.
*
      INTEGER             ncuse
      PARAMETER           (ncuse = 27)
*
      REAL*16             cutoff
      PARAMETER           (cutoff = 13.0q+00)
*
      REAL*16             half
      PARAMETER           (half = 0.5q+00)
*
      REAL*16             one
      PARAMETER           (one = 1.0q+00)
*
      REAL*16             zero
      PARAMETER           (zero = 0.0q+00)
*
*     For x values above x1, psi(x) - ln(x) suffers bit loss:
*
*     Maple V5.1:
*     Digits := 60;
*     x1 := fsolve(Psi(x)/ln(x) = 1/2, x = 1.5 .. 2);
*         1.81953794823878564441886676084334572004393270618128250146815
*
      REAL*16             x1
      PARAMETER (x1 =
     X1.81953794823878564441886676084334572004393270618128250146815Q+00)
*
*     Argument variables
*
      REAL*16             x
*
*     Local variables
*
      INTEGER             i,           n
*
      REAL*16             sum,         xfrac,       xsqinv,      y
      REAL*16             ysqinv
*
      INCLUDE 'qcpsi.inc'
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
*     We should do the computation in octuple precision to reduce the
*     impact of the bit loss in case (b), but that is not feasible in
*     any Fortran implementation.  The result is that about one decimal
*     digit is lost in the range [x1,cutoff).
*
      IF (isqnan(x)) THEN
          qpsiln = qnan()
      ELSE IF (x .LE. zero) THEN
          qpsiln = qnan()
      ELSE IF (x .LT. x1) THEN
          qpsiln = qpsi(x) - qlog(x)
      ELSE IF (x .LT. cutoff) THEN
          xfrac = x - qint(x)
          i = iqceil(cutoff)
          n = i - iqint(x)
          y = qfloat(i) + xfrac
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
              sum = sum - one / (y - qfloat(i))
   20     CONTINUE
          qpsiln = sum - qlog(x/y)
      ELSE
          xsqinv = one/(x * x)
          sum = c(ncuse)
          DO 30 i  = (ncuse - 1), 1, -1
              sum = sum * xsqinv + c(i)
   30     CONTINUE
          qpsiln = -half/x - sum*xsqinv
      END IF
*
      END
