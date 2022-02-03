      SUBROUTINE  gjqf(x, w, y, z, alpha, beta, nquad, ierr)
************************************************************************
*     (Gauss-Jacobi Quadrature with Function values)
*
*     Compute the nodes and weights for the evaluation of the integral
*
*         \int_{-1}^{1} (1-x)^\alpha(1+x)^\beta\ln(1+x)f(x)dx
*                                 (\alpha > -1, \beta > -1)
*
*     as the quadrature sum
*
*         \sum_{i=1}^N[W_i(\alpha,\beta) (\ln 2) f(x_i(\alpha,\beta))
*                      - Z_i(\alpha,\beta)f(y_i(\alpha,\beta))]
*
*     The nonlogarithmic ordinary Gauss-Jacobi integral
*
*         \int_{-1}^{1} (1-x)^\alpha(1+x)^\beta f(x)dx
*                                 (\alpha > -1, \beta > -1)
*
*     can be computed from the quadrature sum
*
*         \sum_{i=1}^N[W_i(\alpha,\beta) f(x_i(\alpha,\beta)]
*
*     The quadrature is exact to machine precision for f(x) of
*     polynomial order less than or equal to 2*nquad - 1.
*
*     This form of the quadrature requires only values of the
*     function, at 2*nquad points.  For a faster, and slightly more
*     accurate, quadrature that requires values of the function and
*     its derivative at nquad points, see the companion routine,
*     gjqfd().
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
*          y(1..nquad)     Nodes of the quadrature for
*                              -(1-x)^\alpha * (1+x)^\beta ln((1+x)/2),
*                          denoted y_i(\alpha,\beta) above.
*
*          z(1..nquad)     Weights of the quadrature for
*                              -(1-x)^\alpha * (1+x)^\beta ln((1+x)/2),
*                          denoted Z_i(\alpha,\beta) above.
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
*          dlgtwo = dlog(2.0d+00)
*          sum = 0.0d+00
*          do 10 i = 1,nquad
*              sum = sum + dlgtwo*w(i)*f(x(i)) - z(i)*f(y(i))
*       10 continue
*
*     The nonlogarithmic integral can be computed by:
*
*          sum = 0.0d+00
*          do 20 i = 1,nquad
*              sum = sum + w(i)*f(x(i))
*       20 continue
*
*     [29-Apr-2000]
************************************************************************
*
*     External functions
*
      EXTERNAL            dvsum
*
      DOUBLE PRECISION    dsqrt,       dvsum
*
*     Parameter variables
*
      DOUBLE PRECISION    ONE
      PARAMETER           (ONE = 1.0D+00)
      DOUBLE PRECISION    ZERO
      PARAMETER           (ZERO = 0.0D+00)
*
      INCLUDE 'dlgtwo.inc'
      INCLUDE 'ecodes.inc'
      INCLUDE 'maxpts.inc'
*
*     Argument variables
*
      DOUBLE PRECISION    alpha,       beta,        w(*),        x(*)
      DOUBLE PRECISION    y(*),        z(*)
*
      INTEGER             ierr,        nquad
*
*     Local variables
*
      DOUBLE PRECISION    a(0:MAXPTS), b(0:MAXPTS), deltaw(MAXPTS)
      DOUBLE PRECISION    deltax(MAXPTS),           e(MAXPTS)
      DOUBLE PRECISION    p(MAXPTS),   pm(MAXPTS),  s(0:MAXPTS)
      DOUBLE PRECISION    t(0:MAXPTS)
*
      INTEGER             m,           n
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
*     Compute recursion coefficients, nodes, and weights.
*
      CALL gjqrc(a, b, s, t, alpha, beta, nquad, ierr)
      IF (ierr .NE. eokay) RETURN
*
*     Setup the tridiagonal Jacobi matrix whose elements are determined
*     by the a(*) and b(*) recursion coefficients.
*
      DO 100 n = 1, nquad
          y(n) = b(n-1)
          e(n) = dsqrt(a(n-1))
  100 CONTINUE
*
*     Diagonalize Jacobi matrix.  On return, y(*) contains the
*     quadrature nodes in ascending order, and e(*) is trashed.
*
      CALL tql1(nquad, y, e, ierr)
      IF (ierr .NE. 0) THEN
          ierr = eeigen
          RETURN
      END IF
*
*     Compute the nodes and weights for order nquad for the first part
*     of the quadrature.
*
      CALL gjqfd(x, w, deltaw, deltax, alpha, beta, nquad, ierr)
      IF (ierr .NE. eokay) RETURN
*
*     Normalize eigenvectors.
*
      DO 300 n = 1, nquad
          p(1) = ONE
          p(2) = (y(n) - b(0))*p(1)
          pm(1) = p(1)*(p(1)/t(0))
          if (nquad .GT. 1) pm(2) = p(2)*(p(2)/t(1))
          DO 200 m = 2, nquad - 1
              p(m + 1) = (y(n) - b(m - 1))*p(m) - a(m - 1)*p(m - 1)
              pm(m + 1) = p(m + 1)*(p(m + 1)/t(m))
  200     CONTINUE
          z(n) = ONE/dvsum(pm, nquad)
          IF (z(n) .EQ. ZERO) THEN
              ierr = eovflo
              RETURN
          END IF
  300 CONTINUE
*
      END
