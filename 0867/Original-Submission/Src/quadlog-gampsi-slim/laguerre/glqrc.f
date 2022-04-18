      SUBROUTINE  glqrc (a, b, s, t, alpha, nquad, ierr)
************************************************************************
*     (Gauss-Laguerre Logarithmic Quadrature Recursion Coefficients)
*
*     Compute the recursion coefficients and zeroth and first moments
*     of the monic polynomials corresponding to the positive weight
*     function
*
*         w(x,\alpha) = (x - 1 - ln(x)) * exp(-x) * x^\alpha
*
*     with recursion relation (n = 0, 1, 2, ...)
*
*         P_{n+1}^\alpha(x) = (x - B_n^\alpha) * P_n^\alpha(x)
*                             - A_n^\alpha * P_{n-1}^\alpha(x)
*
*     and initial conditions
*
*         P_{-1}^\alpha(x) = 0
*         P_{0}^\alpha(x) = 1
*
*     Except in the weight function, the superscripts indicate
*     dependence on \alpha, NOT exponentiation.
*
*     The required moments are:
*
*         T_n^\alpha = \int_0^\infty w(x,\alpha) (P_n^\alpha(x))^2 dx
*         S_n^\alpha = \int_0^\infty w(x,\alpha) (P_n^\alpha(x))^2 x dx
*
*     From these moments, the recursion coefficents are computed as:
*
*         A_n^\alpha = T_n^\alpha / T_{n-1}^\alpha
*         B_n^\alpha = S_n^\alpha / T_n^\alpha
*
*     On entry:
*
*          alpha           Power of x in the integrand (alpha > -1).
*
*          nquad           Number of quadrature points to compute.  It
*                          must be less than the limit MAXPTS defined
*                          in the header file, maxpts.inc.  The default
*                          value chosen there should be large enough
*                          for any realistic application.
*
*     On return:
*
*          a(0..nquad)     Recursion coefficients: a(n) = A_n^\alpha.
*
*          b(0..nquad)     Recursion coefficients: b(n) = B_n^\alpha.
*
*          s(0..nquad)     First moments: s(n) = S_n^\alpha
*
*          t(0..nquad)     Zeroth moments: t(n) = T_n^\alpha
*
*          ierr            Error indicator:
*                          = 0  (success),
*                            1  (eigensolution could not be obtained),
*                            2  (destructive overflow),
*                            3  (nquad out of range),
*                            4  (alpha out of range).
*
*     [27-Oct-2003] -- new code in separate subroutine to allow access
*                      to recursion coefficients and moments
*     [18-Mar-2000] -- original code embedded in glqfd()
************************************************************************
*
*     External functions
*
      EXTERNAL            dgamma,      dpsi,        dvsum
*
      DOUBLE PRECISION    dgamma,      dpsi,        dvsum
*
*     Parameter variables
*
      DOUBLE PRECISION    ONE
      PARAMETER           (ONE = 1.0d+00)
      DOUBLE PRECISION    TWO
      PARAMETER           (TWO = 2.0d+00)
      DOUBLE PRECISION    ZERO
      PARAMETER           (ZERO = 0.0d+00)
*
      INCLUDE 'ecodes.inc'
      INCLUDE 'maxpts.inc'
*
*     Argument variables
*
      DOUBLE PRECISION    a(0:MAXPTS), alpha,       b(0:MAXPTS)
      DOUBLE PRECISION    s(0:MAXPTS), t(0:MAXPTS)
*
      INTEGER             ierr,        nquad
*
*     Local variables
*
      DOUBLE PRECISION    deltaw(MAXPTS),           deltax(MAXPTS)
      DOUBLE PRECISION    p(MAXPTS),   pm(MAXPTS),  pp(MAXPTS)
      DOUBLE PRECISION    ppm(MAXPTS), sn,          tmp,         tmpp
      DOUBLE PRECISION    tn,          w(MAXPTS),   wxm1dw(MAXPTS)
      DOUBLE PRECISION    x(MAXPTS)
*
      INTEGER             m,           n,           nqp2
*
*     Accurate computation of S(nquad) requires a quadrature order at
*     least TWO higher than requested, since S(n) has an extra x
*     factor, and the weight function has an x factor, giving an
*     integrand of polynomial order 2*n + 2 = 2*(n + 1).
*
      nqp2 = nquad + 2
      CALL glqfd(x, w, deltaw, deltax, alpha, nqp2, ierr)
      IF (ierr .NE. eokay) RETURN
*
*     Compute the quadrature weights needed in glqli() and for the
*     initial moments t(0) and s(0).
*
      DO 100 n = 1, nqp2
          wxm1dw(n) = (x(n) - ONE)*w(n) - deltaw(n)
*         g(n) = wxm1dw(n) * x(n)
  100 CONTINUE
*
*     We could compute the initial moments from (exact to
*     machine-precision) quadrature, but it is more accurate to use
*     the analytic expressions in dgamma() and dpsi() given below.
*
*     t(0) = \int_0^\infty w(x,\alpha) (P_0^\alpha(x))^2 dx
*          = \int_0^\infty w(x,\alpha) dx
*          = sum of weights wxm1dw(*)
*          = gamma(\alpha + 1) * (\alpha - psi(\alpha + 1))
*
*     s(0) = \int_0^\infty w(x,\alpha) (P_0^\alpha(x))^2 x dx
*          = \int_0^\infty w(x,\alpha) x dx
*          = sum(i=1:n)(g(i) - deltax(i))
*          = gamma(\alpha + 2) * (\alpha + 1 - psi(\alpha + 2))
*
*     t(0) = dvsum(wxm1dw,nqp2)
*     g(0) = -dvsum(deltax, nqp2)
*     s(0) = dvsum(g, 1 + nqp2)
      t(0) = dgamma(alpha + ONE) * (alpha - dpsi(alpha + ONE))
      s(0) = dgamma(alpha + TWO) * (alpha + ONE - dpsi(alpha + TWO))
*
      a(0) = t(0)
      b(0) = s(0) / t(0)
*
*     Initialize monic polynomial starting values.
*
      DO 200 n = 1, nqp2
          p(n) = ONE
          pm(n) = ZERO
          ppm(n) = ZERO
          pp(n) = ZERO
  200 CONTINUE
*
*     Because a(0), b(0), s(0), and t(0) have already been determined
*     by accurate analytic formulas, on the n = 0 iteration, we avoid
*     replacing them by slightly less accurate values determined by
*     numerical quadrature.
*
      DO 400 n = 0, nquad
          CALL glqli(p, pp, deltax, wxm1dw, x, nqp2, sn, tn)
          IF (n .GT. 0) THEN
              s(n) = sn
              t(n) = tn
              a(n) = t(n) / t(n - 1)
              b(n) = s(n) / t(n)
          END IF
          DO 300 m = 1, nqp2
              tmp = p(m)
              tmpp = pp(m)
              p(m) = (x(m) - b(n))*tmp - a(n)*pm(m)
              pp(m) = ((x(m) - b(n))*tmpp + tmp) - a(n)*ppm(m)
              pm(m) = tmp
              ppm(m) = tmpp
  300     CONTINUE
  400 CONTINUE
*
      END
