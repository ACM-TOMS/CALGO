      SUBROUTINE  gjqfd(x, w, deltaw, deltax, alpha, beta, nquad, ierr)
************************************************************************
*     (Gauss-Jacobi Quadrature with Function and Derivative values)
*
*     Compute the nodes and weights for the evaluation of the integral
*
*         \int_{-1}^{1} (1-x)^\alpha (1+x)^\beta \ln(1+x) f(x) dx
*                                 (\alpha > -1, \beta > -1)
*
*     as the quadrature sum
*
*          \sum_{i=1}^{N}[\deltaW_i(\alpha,\beta) f(\x_i(\alpha,\beta))
*                      + \deltax_i(\alpha,\beta) f'(\x_i(\alpha,\beta))]
*
*     The nonlogarithmic integral
*
*         \int_{-1}^{1} (1-x)^\alpha (1+x)^\beta f(x) dx
*                                 (\alpha > -1, \beta > -1)
*
*     can be computed from the quadrature sum
*
*          \sum_{i=1}^{N}W_i(\alpha, \beta) f(x_i(\alpha, \beta)).
*
*     The quadrature is exact to machine precision for f(x) of
*     polynomial order less than or equal to 2*nquad - 1.
*
*     This form of the quadrature requires values of the function AND
*     ITS DERIVATIVE at N (== nquad) points.  For a derivative-free
*     quadrature at 2N points, see the companion routine, gjqf().
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
*          x(1..nquad)     Nodes of the Jacobi quadrature, denoted
*                          x_i(\alpha,\beta) above.
*
*          w(1..nquad)     Weights of the Jacobi quadrature, denoted
*                          W_i(\alpha,\beta) above.
*
*          deltaw(1..nquad)    Weights of the quadrature, denoted
*                              \deltaW_i(\alpha,\beta) above.
*
*          deltax(1..nquad)    Weights of the quadrature, denoted
*                              \deltax_i(\alpha,\beta) above.
*
*          ierr            Error indicator:
*                          = 0  (success),
*                            1  (eigensolution could not be obtained),
*                            2  (destructive overflow),
*                            3  (nquad out of range),
*                            4  (alpha out of range),
*                            5  (beta out of range).
*
*     The logarithmic integral can then be computed by code like this:
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
      INCLUDE 'dlgtwo.inc'
      INCLUDE 'ecodes.inc'
      INCLUDE 'maxpts.inc'
*
*     Argument variables
*
      DOUBLE PRECISION    alpha,       beta,        deltaw(*)
      DOUBLE PRECISION    deltax(*),   w(*),        x(*)
*
      INTEGER             ierr,        nquad
*
*     Local variables
*
      DOUBLE PRECISION    an,          ap,          bn,          bp
      DOUBLE PRECISION    dd(MAXPTS),  e(MAXPTS),   ee(MAXPTS),  fab
      DOUBLE PRECISION    ff,          fj,          fl
      DOUBLE PRECISION    g(MAXPTS),   h(MAXPTS),   p0,          p1
      DOUBLE PRECISION    p2,          ph0,         ph1,         ph2
      DOUBLE PRECISION    phz,         s,           sab,         sum
      DOUBLE PRECISION    sum1,        sum2,        u0,          u0s
      DOUBLE PRECISION    y(MAXPTS)
*
      INTEGER             i,           j
*
*     Sanity checks on arguments:
*
      IF ((nquad .LE. 0) .OR. (nquad .GE. MAXPTS)) THEN
          ierr = enquad
          RETURN
      ELSE IF (alpha .LE. -ONE) THEN
          ierr = ealpha
          RETURN
      ELSE IF (beta .LE. -ONE) THEN
          ierr = ebeta
          RETURN
      ELSE
          ierr = eokay
      END IF
*
      ap = ONE + alpha
      bp = ONE + beta
      sab = ap + bp
      u0s = dgamma(sab)/(dgamma(ap)*dgamma(bp)*TWO**(ap + beta))
      u0 = dsqrt(u0s)
      fab = alpha**2 - beta**2
      phz = HALF*u0*(dpsi(ap + bp) - dpsi(bp) - DLGTWO)
      DO 100 i = 1, nquad
          fl = dfloat(i)
          an = fl + alpha
          bn = fl + beta
          IF (i .EQ. 1) THEN
              ff = ONE
          ELSE
              ff = (an + beta)/(an + bn - ONE)
          END IF
          e(i + 1) = TWO*dsqrt(fl*an*bn*ff/(an + bn + ONE))/(an + bn)
          x(i) = (-fab)/(an + bn - TWO)/(an + bn)
          g(i) = (alpha - beta)/(an + bn - TWO)
          y(i) = (-x(i))*(ONE/(an + bn) + ONE/(an + bn - TWO)) +
     X        TWO*beta/(an + bn)/(an + bn - TWO)
          dd(i) = x(i)
          ee(i + 1) = e(i + 1)
          h(i + 1) = e(i + 1)*(an + bn + ONE)
  100 CONTINUE
      x(1) = (-(alpha - beta))/(alpha + beta + TWO)
      g(1) = (-x(1))
      dd(1) = x(1)
      y(1) = TWO*ap/(ap + bp)**2
*
*     Diagonalize Jacobi matrix.  On return, x(*) contains the
*     quadrature nodes in ascending order, and e(*) is trashed.
*
      CALL tql1(nquad, x, e, ierr)
      IF (ierr .NE. 0) THEN
          ierr = eeigen
          RETURN
      END IF
*
*     Form perturbation matrix d/da(j).
*
      e(1) = HALF*ee(2)*( (-TWO)/(ap + bp) + ONE/bp - ONE/(ap + bp +
     X    ONE))
      DO 200 j = 2, nquad
          fj = dfloat(j)
          an = fj + alpha
          bn = fj + beta
          e(j) = HALF*ee(j + 1)*( (-TWO)/(an + bn) + ONE/bn + ONE/(an
     X        + beta) - ONE/(an + bn - ONE) - ONE/(an + bn + ONE))
  200 CONTINUE
*
*     The accumulator variables are:
*         sum:   S^\alpha_i
*         sum1:  Q^\alpha_i
*         sum2:  (1 - x)^2 * R^\alpha_i
*         deltax(*): (d/d\alpha)x_i(\alpha,\beta)*S^\alpha_i
*
      DO 400 i = 1, nquad
          p0 = u0
          p1 = u0*(x(i) - dd(1))/ee(2)
          ph0 = phz
          ph1 = ph0/u0*p1 - (u0*y(1) + p1*e(1))/ee(2)
          sum = p0**2
          sum1 = p0*ph0
          sum2 = ZERO
          deltax(i) = (e(1)*p1 + y(1)*p0)*p0
          DO 300 j = 2, nquad
              p2 = ((x(i) - dd(j))*p1 - ee(j)*p0)/ee(j + 1)
              sum = sum + p1**2
              sum1 = sum1 + p1*ph1
              sum2 = sum2 + ((g(j) - x(i))*dfloat(j - 1)*p1 + h(j)*p0)*
     X            p1
              deltax(i) = deltax(i) + (e(j - 1)*p0 + y(j)*p1 +
     X            e(j)*p2)*p1
              s = ( (-e(j))*p2 - y(j)*p1 - e(j - 1)*p0)
              ph2 = ((x(i) - dd(j))*ph1 + s - ee(j)*ph0)/ee(j + 1)
              p0 = p1
              p1 = p2
              ph0 = ph1
              ph1 = ph2
  300     CONTINUE
*
*         Renormalize deltax(*) to recover (d/d\alpha)x_i(\alpha,\beta)
*
          w(i) = ONE/sum
          IF (w(i) .EQ. ZERO) THEN
              ierr = eovflo
              RETURN
          END IF
          deltax(i) = deltax(i) - e(nquad)*p1*p0
          deltax(i) = deltax(i)*w(i)
          IF (deltax(i) .EQ. ZERO) THEN
              ierr = eovflo
              RETURN
          END IF
*
*         Calculate deltaw(x(i))/da = deltaw(i).  We remove the factor
*         (1 - x^2).
*
          deltaw(i) = (sum1 + sum2*deltax(i)/(ONE - x(i)**2))
          deltaw(i) = (-TWO)*w(i)**2*deltaw(i)
          deltax(i) = w(i)*deltax(i)
  400 CONTINUE
*
      END
