      SUBROUTINE  glqfd(x, w, deltaw, deltax, alpha, nquad, ierr)
************************************************************************
*     (Gauss-Laguerre Quadrature with Function and Derivative values)
*
*     Compute the nodes and weights for the evaluation of the integral
*
*         \int_0^\infty x^\alpha e^{-x} \ln(x) f(x) dx
*
*     as the quadrature sum:
*
*         \sum_{i=1}^{N}[\delta W_i(\alpha) f(x_i(\alpha))
*                      + \delta x_i(\alpha) f'(x_i(\alpha))]
*
*     The nonlogarithmic ordinary Gauss-Laguerre integral
*
*         \int_0^\infty x^\alpha e^{-x} f(x) dx
*
*     can be computed from the quadrature sum
*
*         \sum_{i=1}^{N}[\W_i(\alpha) f(x_i(\alpha))]
*
*     The quadrature is exact to machine precision for f(x) of
*     polynomial order less than or equal to 2*nquad - 1.
*
*     This form of the quadrature requires values of the function AND
*     ITS DERIVATIVE at N (== nquad) points.  For a derivative-free
*     quadrature at 2N points, see the companion routine, glqf().
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
*          x(1..nquad)     Nodes of both parts of the quadrature,
*                          denoted x_i(\alpha) above.
*
*          w(1..nquad)     Internal weights of both parts of the
*                          quadrature, denoted W_i(\alpha) above.
*
*          deltaw(1..nquad) Weights of the second part of the
*                          quadrature, denoted \delta W_i(\alpha) above.
*
*          deltax(1..nquad) Weights of the first part of the quadrature,
*                          denoted \delta x_i above.
*
*          ierr            Error indicator:
*                          = 0  (success),
*                            1  (eigensolution could not be obtained),
*                            2  (destructive overflow),
*                            3  (nquad out of range),
*                            4  (alpha out of range).
*
*     The integral can then be computed by code like this:
*
*          sum = 0.0d+00
*          do 10 i = 1,nquad
*              sum = sum + deltaw(i)*f(x(i)) + deltax(i)*fprime(x(i))
*       10 continue
*
*     where fprime(x(i)) is the derivative of the function f(x) with
*     respect to x, evaluated at x == x(i).
*
*     The nonlogarithmic integral can be computed by:
*
*          sum = 0.0d+00
*          do 20 i = 1,nquad
*              sum = sum + w(i)*f(x(i))
*       20 continue
*
*     [20-Mar-2000]
************************************************************************
*
*     External functions
*
      EXTERNAL            dfloat,      dgamma,      dpsi
*
      DOUBLE PRECISION    dfloat,      dgamma,      dpsi,        dsqrt
*
*     Parameter variables
*
      DOUBLE PRECISION    HALF
      PARAMETER           (HALF = 0.5d+00)
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
      DOUBLE PRECISION    alpha,       deltaw(*),   deltax(*),   w(*)
      DOUBLE PRECISION    x(*)
*
      INTEGER             ierr,        nquad
*
*     Local variables
*
      DOUBLE PRECISION    d(MAXPTS),   dd(MAXPTS),  dpsi1
      DOUBLE PRECISION    e(MAXPTS+1), ee(MAXPTS+1),             fj
      DOUBLE PRECISION    fl,          l0alfa,      p0,          p1
      DOUBLE PRECISION    p2,          phi0,        phi1,        phi2
      DOUBLE PRECISION    ps2,         psi0,        qialfa,      rialfa
      DOUBLE PRECISION    s,           sialfa
*
      INTEGER             i,           j,           n
*
*     Sanity checks on arguments:
*
      IF ((nquad .LE. 0) .OR. (nquad .GE. MAXPTS)) THEN
          ierr = enquad
          RETURN
      ELSE IF (alpha .LE. -ONE) THEN
          ierr = ealpha
          RETURN
      ELSE
          ierr = eokay
      END IF
*
*     Setup the tridiagonal Jacobi matrix, T_N{\alpha}, for the
*     generalized Laguerre polynomial.
*
*     Here,
*         d(n  ) == B_{n-1}^\alpha
*         e(n+1) == -sqrt(A_n^\alpha)
*
      DO 100 n = 1, nquad
          fl = dfloat(n)
          d(n) = TWO*fl - ONE + alpha
          e(n + 1) = -dsqrt(fl*(fl + alpha))
  100 CONTINUE
*
*     Save the tridiagonal Jacobi matrix so that we can calculate the
*     eigenvectors later, since it will be destroyed in the call to
*     tql1().
*
      CALL dcopy(nquad,d,1,dd,1)
      CALL dcopy(nquad,e(2),1,ee(2),1)
*
*     tql1() solves the tridiagonal system, and returns the eigenvalues
*     in ascending order in d(*); these are the nodes x_i^\alpha for the
*     generalized Laguerre quadrature.  See glqf() for comments on this
*     particular choice of tridiagonal eigensolver.
*
      CALL tql1(nquad, d, e, ierr)
      IF (ierr .NE. 0) THEN
          ierr = eeigen
          RETURN
      END IF
*
*     Form perturbation matrix d/da(j). Note diagonal element is one.
*     This corresponds to the derivative with respect to \alpha of
*     A_n^\alpha.
*
      DO 200 j = 1, nquad
          fj = dfloat(j)
          e(j) = HALF*ee(j + 1)/(fj + alpha)
  200 CONTINUE
*
*     Form normalized eigenvectors of T_N{\alpha} and dl/dalpha from
*     recurrence formula.
*
*     In this loop, the triple (p0, p1, p2) corresponds to
*     (P_{n-1}^\alpha(x_i^\alpha), P_n^\alpha(x_i^\alpha),
*     P_{n+1}^\alpha(x_i^\alpha)) in the three-term recurrence.
*
*     The triple (phi0, phi1, phi2) corresponds to
*     (\phi_{n-1}^\alpha(x_i^\alpha), \phi_n^\alpha(x_i^\alpha),
*     \phi_{n+1}^\alpha(x_i^\alpha)) in the three-term recurrence.
*
*     deltax(i) accumulates the sum for
*     (dx_i^\alpha/dx)S_i^\alpha
*
*     rialfa is denoted x_i(\alpha) * R_i^{\alpha}.
*
      l0alfa = ONE/dsqrt(dgamma(ONE + alpha))
      dpsi1 = l0alfa/dsqrt(ONE + alpha)
      ps2 = (-HALF)*dpsi(TWO + alpha)
      psi0 = (-HALF)*l0alfa*dpsi(ONE + alpha)
      DO 400 i = 1, nquad
          p0 = l0alfa
          p1 = l0alfa*(d(i) - dd(1))/ee(2)
          phi0 = psi0
          phi1 = p1*ps2 + dpsi1
          sialfa = p0**2
          qialfa = p0*phi0
          rialfa = ZERO
          deltax(i) = (e(1)*p1 + p0)*p0
          DO 300 j = 2, nquad
              p2 = ((d(i) - dd(j))*p1 - ee(j)*p0)/ee(j + 1)
              sialfa = sialfa + p1**2
              qialfa = qialfa + p1*phi1
              rialfa = rialfa + (dfloat(j - 1)*p1 + ee(j)*p0)*p1
              deltax(i) = deltax(i) + (e(j - 1)*p0 + p1 + e(j)*p2)*p1
              s = ((-e(j))*p2 - p1 - e(j - 1)*p0)
              phi2 = ((d(i) - dd(j))*phi1 + s - ee(j)*phi0)/ee(j + 1)
              p0 = p1
              p1 = p2
              phi0 = phi1
              phi1 = phi2
  300     CONTINUE
*
*         Weights for generalized Laguerre quadrature
*
          w(i) = ONE/sialfa
          IF (w(i) .EQ. ZERO) THEN
              ierr = eovflo
              RETURN
          END IF
*
*         Calculate deltax(i) = w_i(\alpha)*(dx_i(\alpha)/d\alpha)
*
          deltax(i) = deltax(i) - e(nquad)*p1*p0
          deltax(i) = deltax(i)*w(i)
*
*         Calculate deltaw(i) = \delta W_i(\alpha)/d\alpha
*
          deltaw(i) = (-(TWO*w(i)**2))*(qialfa + deltax(i)/d(i)*rialfa)
          IF (deltaw(i) .EQ. ZERO) THEN
              ierr = eovflo
              RETURN
          END IF
          deltax(i) = deltax(i)*w(i)
          IF (deltax(i) .EQ. ZERO) THEN
              ierr = eovflo
              RETURN
          END IF
          x(i) = d(i)
  400 CONTINUE
*
      END
