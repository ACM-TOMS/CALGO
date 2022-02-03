      PROGRAM tqgjf2
************************************************************************
*     (Test Gauss-Jacobi Quadrature with Function values)
*
*     This is the driver program for subroutine qgjqf() which returns
*     the nodes and weights necessary for the nonlog-Jacobi
*     quadrature.
*
*     The test functions are (1+x)^p, p = 0 to pmax, with weight
*
*         (1-x)^alpha * (1+x)^beta,  (alpha > -1, beta > -1)
*
*     integrated from -1 to +1.
*
*     Normal input is on stdin, and normal output on stdout.
*
*     However, there is also an auxiliary file created called simply
*     "output.dat"; it contains data suitable for comparing against
*     similar files from other programs in an integral-independent
*     way.
*
*     Lines in this file are either comment lines (beginning with a
*     sharp sign), or blank or empty lines, or whitespace-separated
*     integral data lines of the form
*
*     np p(1) p(2) ... p(np) neval result opt-relerr opt-abserr
*
*     Here, np is the number of parameters that follow in the p(*)
*     values, neval is the number of function evaluations (0 if
*     unknown or unavailable), result is the floating-point computed
*     value of the integral, opt-relerr is an optional relative error
*     estimate, and opt-abserr is an optional absolute error estimate.
*     If opt-relerr is omitted, then opt-abserr cannot be specified.
*     A zero value for either implies that the value is unknown.
*
*     Normally, comment lines will document what integral is
*     evaluated, and relate any parameters in the integrand to the
*     array elements p(*).  It is acceptable for np to be 0: no p(*)
*     elements are then provided.
*
*     The availability of data files in this standard format makes it
*     relatively easy to compare results of different integrators.  In
*     particular, high-precision values can be computed in symbolic
*     algebra systems, such as Maple, Mathematica, Axiom, Reduce,
*     muPAD, ..., and used to evaluate the accuracy of results from
*     other integrators.
*
*     (13-May-2000)
************************************************************************
*
*     External functions
*
      EXTERNAL            qeps,        qerbit,      qgamma,      qnan
      EXTERNAL            qvsum
*
      REAL*16             qabs,        qeps,        qerbit,      qgamma
      REAL*16             qnan,        qvsum
*
*     Statement functions
*
      REAL*16             f
*
*     Parameter variables
*
      CHARACTER*(*)       METHOD
      PARAMETER           (METHOD = 'qgjqf()')
*
      REAL*16             ONE
      PARAMETER           (ONE = 1.0q+00)
*
      REAL*16             TWO
      PARAMETER           (TWO = 2.0q+00)
*
      REAL*16             ZERO
      PARAMETER           (ZERO = 0.0q+00)
*
      INCLUDE 'maxpts.inc'
      INCLUDE 'stdio.inc'
*
*     Local variables
*
      REAL*16             alpha,       beta,        exact,       relavg
      REAL*16             relerr,      relmax,      relmin,      relulp
      REAL*16             result,      ulp,         v(MAXPTS)
      REAL*16             w(MAXPTS),   x(MAXPTS),   xarg
      REAL*16             y(MAXPTS),   z(MAXPTS),   zinit(1)
*
      INTEGER             i,           ierr,        karg,        neval
      INTEGER             nqmax,       nqmin,       nquad,       nruns
      INTEGER             p,           pmax
*
      f(xarg,karg) = (ONE + xarg)**karg
*
*     Sun floating-point error trapping:
*
*     INTEGER             ieee_handler
*     EXTERNAL            trapit
*     ierr = ieee_handler ( 'set', 'invalid', trapit )
*     if (ierr .ne. 0) print *,'failed to set trap for invalid'
*
*     Initialize all local floating-point variables to NaN, or at
*     least an approximation thereto.
*
      zinit(1) = qnan()
      CALL qcopy (MAXPTS, zinit, 0, w, 1)
      CALL qcopy (MAXPTS, zinit, 0, x, 1)
      CALL qcopy (MAXPTS, zinit, 0, y, 1)
      CALL qcopy (MAXPTS, zinit, 0, z, 1)
      alpha = zinit(1)
      beta = zinit(1)
      exact = zinit(1)
      relavg = zinit(1)
      relerr = zinit(1)
      relmax = zinit(1)
      relmin = zinit(1)
      result = zinit(1)
      ulp = zinit(1)
*
      OPEN (UNIT=stddat, FILE='output.dat', STATUS='unknown',
     X    FORM='formatted')
*
      WRITE (stddat,'(A)') '### Numerical integration with ' // METHOD
      WRITE (stddat,'(A)') '###'
      WRITE (stddat,'(A)') '### int((1+x)^p * (1-x)^alpha * ' //
     X    '(1+x)^beta, x = -1..+1)'
      WRITE (stddat,'(A)') '###'
      WRITE (stddat,'(A)')
     X    '### Line format: 3 p alpha beta neval result relerr abserr'
      WRITE (stddat,'(A)') '###'
*
      ulp = qeps(ONE)
      nqmax = 0
      nqmin = 0
      nruns = 0
      relavg = ZERO
      relmax = ZERO
      relmin = ONE
*
*     Write a standard output header identifying the integral, and the
*     host precision
*
      CALL qprthd(stdout,
     X    'int(x^p * (1-x)^alpha * (1+x)^beta, x = -1..+1)')
*
  100 READ (stdin, *, END=500, ERR=500) nquad, alpha, beta, pmax
      WRITE (stdout, 50000) alpha, beta, nquad
*
*     Quadrature tests with
*
*         -(1-x)^\alpha * (1+x)^\beta
*
*     where
*
*         f(x) = (1+x)^p
*
*     and p is the loop index below.
*
*     This integral can be done exactly:
*
*         2^(\alpha + \beta + p + 1) * Gamma(\alpha + 1) *
*         Gamma(\beta + p + 1) / Gamma(\alpha + \beta + p + 2)
*
*     In this nonlogarithmic case, there is no quadrature sum term
*     involving (y(*),z(*)).
*
      CALL qgjqf(x, w, y, z, alpha, beta, nquad, ierr)
      IF (ierr .NE. 0) THEN
          WRITE (stderr, 80000) METHOD, ierr
          GO TO 100
      END IF
      WRITE (stdout, 10000)
      DO 200 i = 1, nquad
          WRITE (stdout, 40000) i, x(i), w(i), y(i), z(i)
  200 CONTINUE
      WRITE (stdout, 30000)
      DO 400 p = 0, pmax
          DO 300 i = 1, nquad
              v(i) = w(i) * f(x(i),p)
  300     CONTINUE
          neval = nquad
          result = qvsum(v, nquad)
          exact = TWO**(alpha + beta + qfloat(p) + 1) *
     X        qgamma(alpha + 1) * qgamma(beta + qfloat(p) + 1) /
     X        qgamma(alpha + beta + p + 2)
          relerr = (exact - result)/exact
          relulp = relerr/ulp
          IF (qabs(relulp) .LT. 10000.0q+00) THEN
              WRITE (stdout, 20000) p, result, exact, relerr, relulp,
     X            qerbit(relerr,ulp)
          ELSE
              WRITE (stdout, 25000) p, result, exact, relerr, relulp,
     X            qerbit(relerr,ulp)
          END IF
          WRITE (stddat, 70000) 3, p, alpha, beta, neval, result,
     X        relerr, (exact - result)
          IF (relerr .eq. relerr) THEN
*
*             Track the average, minimum, and maximum relative error as
*             long as it is not a NaN.
*
              relavg = relavg + qabs(relerr)
              nruns = nruns + 1
              IF (qabs(relerr) .GT. qabs(relmax)) THEN
                  relmax = relerr
                  nqmax = nquad
              END IF
              IF (qabs(relerr) .LT. qabs(relmin)) THEN
                  relmin = relerr
                  nqmin = nquad
              END IF
          END IF
  400 CONTINUE
      GO TO 100
*
  500 IF (nruns .GT. 0) THEN
          relavg = relavg / qfloat(nruns)
          WRITE (stdout, '()')
          WRITE (stdout, 90000) 'Maximum', relmax, relmax/ulp,
     X        qerbit(relmax,ulp), 'at nquad =', nqmax
          WRITE (stdout, 90000) 'Minimum', relmin, relmin/ulp,
     X        qerbit(relmin,ulp), 'at nquad =', nqmin
          WRITE (stdout, 90000) 'Average', relavg, relavg/ulp,
     X        qerbit(relavg,ulp)
      END IF
*
      WRITE (stdout, 60000)
      CLOSE (UNIT=stddat)
*
10000 FORMAT (/, 2X, 4X, 'i', 16X, 'x(i)', 20X, 3X, 16X, 'w(i)', 24X,
     X    16X, 'y(i)', 25X, 15X, 'z(i)')
20000 FORMAT (2X, I5, 1P, 1X, E43.34, 1X, E43.34, 1X, E10.2, 2X, 0P,
     X    1X, F9.2, 1X, F9.2)
25000 FORMAT (2X, I5, 1P, 1X, E43.34, 1X, E43.34, 1X, E10.2, 2X,
     X    1X, E9.2, 0P, 1X, F9.2)
30000 FORMAT (/, 2X, 4X, 'p', 15X, 'Quadrature Result', 13X,
     X    14X, 'Exact Integral', 16X, 'Rel. Error', 1X, 'RelE (ULPs)',
     X    2X, 'Err (bits)')
40000 FORMAT (2X, I5, 1X, F40.34, 3X, 1P, 1X, E43.34, 1X, E43.34, 1X,
     X    E43.34)
50000 FORMAT (/, 2X, 'alpha = ', F20.16, 2X, 'beta = ', F20.16,
     X    5X, 'nquad = ', I6)
60000 FORMAT (/, 2X, 'Done')
70000 FORMAT (I1, 1X, I4, 1X, 1P, 1X, E44.35, 1X, E44.35, 1X, I10, 1X,
     X    E44.35, 1X, E9.2, 1X, E9.2)
80000 FORMAT (/, 2X, 'ERROR: ', A, ' returns ierr = ', I10)
90000 FORMAT (2X, A, ' relative error = ', 1P, E10.2, 1X,
     X    0P, F10.2, ' ulps', 1X, F10.2, ' bits', 1X, A, I5)
*
      END
