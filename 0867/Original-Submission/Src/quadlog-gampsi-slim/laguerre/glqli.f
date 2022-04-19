      SUBROUTINE  glqli(p, pprime, deltax, wxm1dw, x, nquad, sn, tn)
************************************************************************
*     (Gauss-Laguerre Log Quadrature - Integrals)
*
*     Compute two integrals needed by glqf() for computing the recursion
*     coefficients for monic orthogonal polynomials, P_n^{\alpha}(x),
*     corresponding to the positive weight function (x - 1 - ln(x)) *
*     exp(-x) * x**alpha.
*
*     On entry:
*
*          p(1..nquad)          Array of values p(i) =
*                               P_n^{\alpha}(x(i)).  The caller steps
*                               n recursively on successive calls to
*                               this routine.
*
*          pprime(1..nquad)     Array of values pprime(i) =
*                               dP_n^{\alpha}(x)/dx at x = x(i).
*
*          wxm1dw(1..nquad)     Array of values wxm1dw(i) =
*                               (w(i) * (x(i) - 1) - deltaw(i)), where
*                               w(*), x(*), and deltaw(*) are quadrature
*                               weights and nodes from glqfd().
*
*          deltax(1..nquad)     Array of quadrature weights.
*
*          x(1..nquad)          Array of quadrature nodes.
*
*          nquad                Number of quadrature points to compute.
*                               It must be less than the limit MAXPTS
*                               defined in the header file,
*                               maxpts.inc.  The default value chosen
*                               there should be large enough for any
*                               realistic application.
*
*     On return:
*
*          sn        \int_0^\infty w(x,\alpha) (P_n^\alpha(x))^2 x dx
*
*          tn        \int_0^\infty w(x,\alpha) (P_n^\alpha(x))^2 dx
*
*     [18-Mar-2000]
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
      DOUBLE PRECISION    TWO
      PARAMETER           (TWO = 2.0d+00)
*
      INCLUDE 'maxpts.inc'
*
*     Argument variables
*
      DOUBLE PRECISION    deltax(*),   p(*),        pprime(*),   sn
      DOUBLE PRECISION    tn,          wxm1dw(*),   x(*)
*
      INTEGER             nquad
*
*     Local variables
*
      DOUBLE PRECISION    snterm(3,MAXPTS),         tnterm(2,MAXPTS)
*
      INTEGER             i
*
*     The values are computed by quadrature as:
*
*         S(n) = sum_{i=1}^{n}[(w(i) * (x(i) - 1) - deltaw(i)) * g_n(x)
*                              - deltax(i)* gprime_n(x)]
*         T(n) = sum_{i=1}^{n}[(w(i) * (x(i) - 1) - deltaw(i)) * f_n(x)
*                              - deltax(i)* fprime_n(x)]
*
*     where
*
*         f_n(x) =     (P_n^{\alpha}(x))^2
*         g_n(x) = x * (P_n^{\alpha}(x))^2
*
*     and thus
*
*         fprime_n(x) =     2 * P_n^{\alpha}(x) * Pprime_n^{\alpha}(x)
*         gprime_n(x) = (P_n^{\alpha}(x))^2 +
*                        2 * x * P_n^{\alpha}(x) * Pprime_n^{\alpha}(x)
*
*     The sums for sn and tn here contain terms of the form
*
*         x(i)*p(i)**2*wxm1dw(i)
*         (p(i)**2 - TWO*x(i)*p(i)*pprime(i))*deltax(i)
*         wxm1dw(i)*p(i)**2
*
*     For large n, p(n) >> 1, pprime(n) >> 1, and for modest alpha,
*     wxm1dw(n) << 1, and deltax(n) << 1, so there is the very real
*     danger of premature overflow from p(n)**2 or p(n)*pprime(n).
*
*     We therefore parenthesize subexpressions to keep them of
*     reasonable magnitude, and use dvsum() for accurate vector
*     summation.
*
      DO 100 i = 1, nquad
*
* Original code:
*         fi = x(i)*p(i)**2
*         dfi = p(i)**2 + TWO*x(i)*p(i)*pprime(i)
*         sn = sn + wxm1dw(i)*fi - dfi*deltax(i)
*         tn = tn + wxm1dw(i)*p(i)**2 - TWO*p(i)*pprime(i)*deltax(i)
*
          snterm(1,i) = ((wxm1dw(i) * p(i)) * p(i)) * x(i)
          snterm(2,i) = -(deltax(i) * p(i)) * p(i)
          snterm(3,i) = -TWO * ((deltax(i) * p(i)) * pprime(i)) * x(i)
*
*         tn = tn + p(i)*((wxm1dw(i)*p(i)) - TWO*pprime(i)*deltax(i))
*
          tnterm(1,i) = ((wxm1dw(i) * p(i)) * p(i))
          tnterm(2,i) = -TWO * ((deltax(i) * p(i)) * pprime(i))
  100 CONTINUE
      sn = dvsum(snterm, 3 * nquad)
      tn = dvsum(tnterm, 2 * nquad)
*
      END
