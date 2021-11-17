      SUBROUTINE  gjqli(p, pprime, deltax, wdw, x, nquad, sn, tn)
************************************************************************
*     (Gauss-Jacobi Log Quadrature - Integrals)
*
*     Compute two integrals needed by gjqf() for computing the recursion
*     coefficients for monic orthogonal polynomials, P_n^{\alpha}(x),
*     corresponding to the positive weight function -(1 - x)**alpha *
*     (1 + x)**beta ln((1+x)/2) where alpha > -1 and beta > -1.
*
*     On entry:
*
*
*          p(1..nquad)          Array of values p(i) = P_n^{\alpha}(x).
*                               The caller steps n recursively on
*                               successive calls to this routine.
*
*          pprime(1..nquad)     Array of values pprime(i) =
*                               dP_n^{\alpha}(x)/dx at x = x(i).
*
*          wdw(1..nquad)        Array of values wdx(i) = dlog(2.0)*w(i)
*                               - deltaw(i), where w(*) and deltaw(*)
*                               are quadrature weights and nodes from
*                               gjqfd().
*
*          deltax(1..nquad)     Array of quadrature weights.
*
*          x(1..nquad)          Array of quadrature nodes.
*
*          nquad                Number of quadrature points.
*
*     On return:
*
*          sn    \int_{-1}^{+1} w(x,\alpha,\beta)
*                               (P_n^{\alpha,\beta}(x))^2 x dx
*
*          tn    \int_{-1}^{+1} w(x,\alpha,\beta)
*                               (P_n^{\alpha,\beta}(x))^2 dx
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
      DOUBLE PRECISION    ZERO
      PARAMETER           (ZERO = 0.0d+00)
*
      INCLUDE 'maxpts.inc'
*
*     Argument variables
*
      DOUBLE PRECISION    deltax(*),   p(*),        pprime(*),   sn
      DOUBLE PRECISION    tn,          wdw(*),      x(*)
*
      INTEGER             nquad
*
*     Local variables
*
      DOUBLE PRECISION    snterm(3,MAXPTS),         tnterm(2,MAXPTS)
*
      INTEGER             i
*
      sn = ZERO
      tn = ZERO
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
      DO 100 i = 1, nquad
*
* Original code:
*         fi = x(i)*p(i)**2
*         dfi = p(i)**2 + TWO*x(i)*p(i)*pprime(i)
*         sn = sn + wdw(i)*fi - dfi*deltax(i)
*
          snterm(1,i) = ((wdw(i) * p(i)) * p(i)) * x(i)
          snterm(2,i) = -(deltax(i) * p(i)) * p(i)
          snterm(3,i) = -TWO * ((deltax(i) * p(i)) * pprime(i)) * x(i)
*
*         tn = tn + p(i)*((wdw(i)*p(i)) - TWO*pprime(i)*deltax(i))
*
          tnterm(1,i) = ((wdw(i) * p(i)) * p(i))
          tnterm(2,i) = -TWO * ((deltax(i) * p(i)) * pprime(i))
  100 CONTINUE
      sn = dvsum(snterm, 3 * nquad)
      tn = dvsum(tnterm, 2 * nquad)
*
      END
