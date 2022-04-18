      SUBROUTINE  gjqrc (a, b, s, t, alpha, beta, nquad, ierr)
************************************************************************
*     (Gauss-Jacobi Logarithmic Quadrature Recursion Coefficients)
*
*     Compute the recursion coefficients and zeroth and first moments
*     of the monic polynomials corresponding to the positive weight
*     function
*
*         w(x,\alpha,\beta) = (1-x)^\alpha (1+x)^\beta (-ln((1+x)/2))
*
*     with recursion relation (n = 0, 1, 2, ...)
*
*         P_{n+1}^{\alpha,\beta}(x) =
*                       (x - B_n^{\alpha,\beta}) * P_n^{\alpha,\beta}(x)
*                       - A_n^{\alpha,\beta} * P_{n-1}^{\alpha,\beta}(x)
*
*     and initial conditions
*
*         P_{-1}^{\alpha,\beta}(x) = 0
*         P_{0}^{\alpha,\beta}(x) = 1
*
*     Except in the weight function, the superscripts indicate
*     dependence on \alpha, NOT exponentiation.
*
*     The required moments are:
*
*         T_n^{\alpha,\beta} = \int_0^\infty w(x,\alpha,\beta) *
*                                           (P_n^{\alpha,\beta}(x))^2 dx
*         S_n^{\alpha,\beta} = \int_0^\infty w(x,\alpha,\beta) *
*                                         (P_n^{\alpha,\beta}(x))^2 x dx
*
*     From these moments, the recursion coefficents are computed as:
*
*         A_n^{\alpha,\beta} =
*                 T_n^{\alpha,\beta} / T_{n-1}^{\alpha,\beta}
*
*         B_n^{\alpha,\beta} = S_n^{\alpha,\beta} / T_n^{\alpha,\beta}
*
*     On entry:
*
*          alpha           Power of (1-x) in the integrand (alpha > -1).
*
*          beta            Power of (1+x) in the integrand (beta > -1).
*
*          nquad           Number of quadrature points to compute.  It
*                          must be less than the limit MAXPTS defined
*                          in the header file, maxpts.inc.  The default
*                          value chosen there should be large enough
*                          for any realistic application.
*
*     On return:
*
*          a(0..nquad)     Recursion coefficients:
*                          a(n) = A_n^{\alpha,\beta}.
*
*          b(0..nquad)     Recursion coefficients:
*                          b(n) = B_n^{\alpha,\beta}.
*
*          s(0..nquad)     First moments: s(n) = S_n^{\alpha,\beta}.
*
*          t(0..nquad)     Zeroth moments: t(n) = T_n^{\alpha,\beta}.
*
*          ierr            Error indicator:
*                          = 0  (success),
*                            1  (eigensolution could not be obtained),
*                            2  (destructive overflow),
*                            3  (nquad out of range),
*                            4  (alpha out of range).
*
*     [30-Oct-2003] -- new code in separate subroutine to allow access
*                      to recursion coefficients and moments
*     [18-Mar-2000] -- original code embedded in gjqfd()
************************************************************************
*
*     External functions
*
      EXTERNAL            dvsum
*
      DOUBLE PRECISION    dvsum
*
*     Parameter variables
*
      DOUBLE PRECISION    ONE
      PARAMETER           (ONE = 1.0d+00)
      DOUBLE PRECISION    ZERO
      PARAMETER           (ZERO = 0.0d+00)
*
      INCLUDE 'dlgtwo.inc'
      INCLUDE 'ecodes.inc'
      INCLUDE 'maxpts.inc'
*
*     Argument variables
*
      DOUBLE PRECISION    a(0:MAXPTS), alpha,       b(0:MAXPTS), beta
      DOUBLE PRECISION    s(0:MAXPTS), t(0:MAXPTS)
*
      INTEGER             ierr,        nquad
*
*     Local variables
*
      DOUBLE PRECISION    deltaw(MAXPTS),           deltax(MAXPTS)
      DOUBLE PRECISION    p(MAXPTS),   pm(MAXPTS),  pp(MAXPTS)
      DOUBLE PRECISION    ppm(MAXPTS), sn,          stmp(3,MAXPTS)
      DOUBLE PRECISION    tmp,         tmpp,        tn
      DOUBLE PRECISION    ttmp(2,MAXPTS),           w(MAXPTS)
      DOUBLE PRECISION    wdw(MAXPTS), x(MAXPTS)
*
      INTEGER             m,           n,           nqp1
*
*     Accurate computation of S(nquad) requires a quadrature order at
*     least one higher than requested, since S(n) has an extra x
*     factor, giving an integrand of polynomial order 2 * nquad + 1.
*
      nqp1 = nquad + 1
      CALL gjqfd(x, w, deltaw, deltax, alpha, beta, nqp1, ierr)
      IF (ierr .NE. eokay) RETURN
*
*     Compute the quadrature weights needed in gjqli() and for the
*     initial moments t(0) and s(0).
*
      DO 100 n = 1, nqp1
          wdw(n) = DLGTWO*w(n) - deltaw(n)
          ttmp(1,n) = DLGTWO * w(n)
          ttmp(2,n) = -deltaw(n)
          stmp(1,n) = DLGTWO * w(n) * x(n)
          stmp(2,n) = -deltaw(n) * x(n)
          stmp(3,n) = -deltax(n)
  100 CONTINUE
*
*     t(0) = \int_0^\infty w(...) (P_0^{\alpha,\beta}(x))^2 dx
*          = \int_0^\infty w(...) (1)^2 dx
*          = \int_0^\infty w(...) dx
*          = sum(n=1:nqp1)[log(2)*w(n) - deltaw(n)]
*
*     s(0) = \int_0^\infty w(...) (P_0^{\alpha,\beta}(x))^2 x dx
*          = \int_0^\infty w(...) (1)^2 x dx
*          = \int_0^\infty w(...) x dx
*          = sum(n=1:nqp1)[log(2)*w(n)*x(n) - deltaw(n)*x(n) -
*                           deltax(n)]
*
      t(0) = dvsum(ttmp, 2*nqp1)
      s(0) = dvsum(stmp, 3*nqp1)
      a(0) = t(0)
      b(0) = s(0) / t(0)
*
*     Initialize starting values of the monic polynomials.
*     pp(*) is the derivative of p(*) with respect to x.
*
      DO 200 n = 1, nqp1
          p(n) = ONE
          pm(n) = ZERO
          ppm(n) = ZERO
          pp(n) = ZERO
  200 CONTINUE
*
*     Because a(0), b(0), s(0), and t(0) have already been determined,
*     on the n = 0 iteration, we avoid replacing them.
*
      DO 400 n = 0, nquad
          CALL gjqli(p, pp, deltax, wdw, x, nqp1, sn, tn)
          IF (n .GT. 0) THEN
              s(n) = sn
              t(n) = tn
              a(n) = t(n) / t(n - 1)
              b(n) = s(n) / t(n)
          END IF
          DO 300 m = 1, nqp1
              tmp = p(m)
              tmpp = pp(m)
              p(m) = (x(m) - b(n))*tmp - a(n)*pm(m)
              pp(m) = (x(m) - b(n))*tmpp + tmp - a(n)*ppm(m)
              pm(m) = tmp
              ppm(m) = tmpp
  300     CONTINUE
  400 CONTINUE
*
      END
