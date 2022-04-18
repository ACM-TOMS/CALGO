      SUBROUTINE  glqf(x, w, wxm1, y, z, alpha, nquad, ierr)
************************************************************************
*     (Gauss-Laguerre Quadrature with Function values)
*
*     Compute the nodes and weights for the evaluation of the integral
*
*         \int_0^\infty x^\alpha e^{-x} \ln(x) f(x) dx
*
*     as the quadrature sum
*
*         \sum_{i=1}^N[W_i(\alpha)(x_i(\alpha) -1)f(x_i(\alpha)) -
*                      Z_i(\alpha)f(y_i(\alpha))]
*
*     The nonlogarithmic integral
*
*         \int_0^\infty x^\alpha e^{-x} f(x) dx
*
*     can be computed from the quadrature sum
*
*         \sum_{i=1}^N[W_i(\alpha) f(x_i(\alpha))]
*
*     The quadrature is exact to machine precision for f(x) of
*     polynomial order less than or equal to 2*nquad - 2 (logarithmic)
*     or 2*nquad - 1 (nonlogarithmic).
*
*     This form of the quadrature requires only values of the
*     function at 2*nquad points.  For a faster, and slightly more
*     accurate, quadrature that requires values of the function and
*     its derivative at nquad points, see the companion routine,
*     glqfd().
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
*          x(1..nquad)     Nodes of the first part of the quadrature,
*                          denoted x_i(\alpha) above.
*
*          w(1..nquad)     Weights of the first part of the quadrature,
*                          denoted W_i(\alpha) above.
*
*          wxm1(1..nquad)  Scaled weights of the first part of the
*                          quadrature, wxm1(i) = w(i)*(x(i) - 1).
*
*          y(1..nquad)     Nodes of the second part of the quadrature,
*                          denoted y_i(\alpha) above.
*
*          z(1..nquad)     Weights of the second part of the quadrature,
*                          denoted -Z_i(\alpha) above.
*
*          ierr            Error indicator:
*                          = 0  (success),
*                            1  (eigensolution could not be obtained),
*                            2  (destructive overflow),
*                            3  (nquad out of range),
*                            4  (alpha out of range).
*
*     The logarithmic integral can then be computed by code like this:
*
*          sum = 0.0d+00
*          do 10 i = 1,nquad
*              sum = sum + wxm1(i)*f(x(i)) - z(i)*f(y(i))
*       10 continue
*
*     The nonlogarithmic integral can be computed by:
*
*          sum = 0.0d+00
*          do 20 i = 1,nquad
*              sum = sum + w(i)*f(x(i))
*       20 continue
*
*     [18-Mar-2000]
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
      PARAMETER           (ONE = 1.0d+00)
      DOUBLE PRECISION    ZERO
      PARAMETER           (ZERO = 0.0d+00)
*
      INCLUDE 'ecodes.inc'
      INCLUDE 'maxpts.inc'
*
*     Argument variables
*
      DOUBLE PRECISION    alpha,       w(*),        wxm1(*),     x(*)
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
      CALL glqrc(a, b, s, t, alpha, nquad, ierr)
      IF (ierr .NE. eokay) RETURN
*
      DO 100 n = 1, nquad
          y(n) = b(n-1)
          e(n) = dsqrt(a(n-1))
  100 CONTINUE
*
*     Diagonalize the tridiagonal Jacobi matrix using a standard EISPACK
*     routine, tql1().  We could equivalently use another EISPACK
*     routine, imtql1(), with the same calling sequence, or EISPACK
*     routine tqlrat() with a slightly different calling sequence, or
*     the more modern LAPACK routine, dstev(), with a different calling
*     sequence.  We choose the older EISPACK code simply because it is
*     compact, needing only a single external function, while the LAPACK
*     routine calls about a dozen other LAPACK routines; we do not want
*     the user of this package to first have to install either of those
*     large and complex packages.
*
*     On return, the off-diagonal elements, e(*), have been destroyed,
*     and the diagonal elements, y(*), have been replaced by the
*     eigenvalues in ascending order, which are the nodes of the
*     quadrature.
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
      CALL glqfd(x, w, deltaw, deltax, alpha, nquad, ierr)
      IF (ierr .NE. eokay) RETURN
*
*     Compute the weights of the second part of the quadrature.
*
      DO 300 n = 1, nquad
          p(1) = ONE
          p(2) = (y(n) - b(0))*p(1)
*
*         Normalize the eigenvectors.  For large n, p(n) >> 1, pp(n) >>
*         1, and t(n) >> 1, so there is the very real danger of
*         premature overflow from p(n)**2.
*
*         We therefore parenthesize subexpressions to keep them of
*         reasonable magnitude.
*
*         Original code:
*         sum = p(1)**2/t(1) + p(2)**2/t(2)
*
          pm(1) = p(1)*(p(1)/t(0))
          if (nquad .GT. 1) pm(2) = p(2)*(p(2)/t(1))
*
*         The sum here is the value S_n(\alpha).
*
          DO 200 m = 2, nquad - 1
              p(m + 1) = (y(n) - b(m - 1))*p(m) - a(m - 1)*p(m - 1)
*
*             Original code:
*             sum = sum + p(m + 1)**2/t(m + 1)
*
              pm(m + 1) = p(m + 1)*(p(m + 1)/t(m))
  200     CONTINUE
          wxm1(n) = (x(n) - ONE)*w(n)
          z(n) = ONE/dvsum(pm, nquad)
          IF (z(n) .EQ. ZERO) THEN
              ierr = eovflo
              RETURN
          END IF
  300 CONTINUE
*
      END
