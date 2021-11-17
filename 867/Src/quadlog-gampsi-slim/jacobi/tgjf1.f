      PROGRAM tgjf1
************************************************************************
*     (Test Gauss-Jacobi Quadrature with Function values)
*
*     This is the driver program for subroutine gjqf() which returns
*     the nodes and weights necessary for the log-Jacobi quadrature.
*
*     The test functions are x^p, p = 0 to pmax, with weight
*
*         (1-x)^alpha * (1+x)^beta * log((1+x)/2)
*                                                (alpha > -1, beta > -1)
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
      EXTERNAL            deps,        derbit,      dfloat,      dgamma
      EXTERNAL            dnan,        dpsi,        dvsum
*
      DOUBLE PRECISION    dabs,        deps,        derbit,      dfloat
      DOUBLE PRECISION    dgamma,      dnan,        dpsi,        dvsum
*
*     Statement functions
*
      DOUBLE PRECISION    f
*
*     Parameter variables
*
      CHARACTER*(*)       METHOD
      PARAMETER           (METHOD = 'gjqf()')
*
      DOUBLE PRECISION    ONE
      PARAMETER           (ONE = 1.0d+00)
*
      DOUBLE PRECISION    TWO
      PARAMETER           (TWO = 2.0d+00)
*
      DOUBLE PRECISION    ZERO
      PARAMETER           (ZERO = 0.0d+00)
*
      INCLUDE 'maxpts.inc'
      INCLUDE 'stdio.inc'
*
*     Local variables
*
      DOUBLE PRECISION    alpha,       beta,        exact,       pp1
      DOUBLE PRECISION    relavg,      relerr,      relmax,      relmin
      DOUBLE PRECISION    relulp,      result,      v(MAXPTS)
      DOUBLE PRECISION    ulp,         w(MAXPTS),   x(MAXPTS),   xarg
      DOUBLE PRECISION    y(MAXPTS),   z(MAXPTS),   zinit(1)
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
      zinit(1) = dnan()
      CALL dcopy (MAXPTS, zinit, 0, w, 1)
      CALL dcopy (MAXPTS, zinit, 0, x, 1)
      CALL dcopy (MAXPTS, zinit, 0, y, 1)
      CALL dcopy (MAXPTS, zinit, 0, z, 1)
      alpha = zinit(1)
      beta = zinit(1)
      exact = zinit(1)
      pp1 = zinit(1)
      relavg = zinit(1)
      relerr = zinit(1)
      relmax = zinit(1)
      relmin = zinit(1)
      result = zinit(1)
      ulp = zinit(1)
      xarg = zinit(1)
*
      OPEN (UNIT=stddat, FILE='output.dat', STATUS='unknown',
     X    FORM='formatted')
*
      WRITE (stddat,'(A)') '### Numerical integration with ' // METHOD
      WRITE (stddat,'(A)') '###'
      WRITE (stddat,'(A)') '### int((1+x)^p * (1-x)^alpha * ' //
     X    '(1+x)^beta * ln((1+x)/2), x = -1..+1)'
      WRITE (stddat,'(A)') '###'
      WRITE (stddat,'(A)')
     X    '### Line format: 3 p alpha beta neval result relerr abserr'
      WRITE (stddat,'(A)') '###'
*
      ulp = deps(ONE)
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
      CALL prthdr(stdout, 'int((1+x)^p * (1-x)^alpha * ' //
     X    '(1+x)^beta * ln((1+x)/2),  x = -1..+1)')
*
  100 READ (stdin, *, END=500, ERR=500) nquad, alpha, beta, pmax
      WRITE (stdout, 50000) alpha, beta, nquad
*
*     Quadrature tests with
*
*         -(1-x)^\alpha * (1+x)^\beta * ln((1 + x)/2) * f(x)
*
*     where
*
*         f(x) = (1+x)^p
*
*     and p is the loop index below.
*
*     This integral can be done exactly; see Eq. (19) in paper 2.
*
*     There is no quadrature sum term involving (x(*),w(*)).
*
      CALL gjqf(x, w, y, z, alpha, beta, nquad, ierr)
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
          pp1 = dfloat(p + 1)
          DO 300 i = 1, nquad
              v(i) =  -z(i) * f(y(i),p)
  300     CONTINUE
          neval = nquad
          result = dvsum(v, nquad)
          exact = TWO**(alpha + beta + pp1)*dgamma(beta + pp1)*
     X        dgamma(alpha + ONE)/dgamma(alpha + beta + ONE + pp1)*
     X        (dpsi(pp1 + beta) - dpsi(alpha + beta + pp1 + ONE))
          relerr = (exact - result)/exact
          relulp = relerr/ulp
          IF (dabs(relulp) .LT. 10000.0d+00) THEN
              WRITE (stdout, 20000) p, result, exact, relerr, relulp,
     X            derbit(relerr,ulp)
          ELSE
              WRITE (stdout, 25000) p, result, exact, relerr, relulp,
     X            derbit(relerr,ulp)
          END IF
          WRITE (stddat, 70000) 3, p, alpha, beta, neval, result,
     X        relerr, (exact - result)
          IF (relerr .eq. relerr) THEN
*
*             Track the average, minimum, and maximum relative error as
*             long as it is not a NaN.
*
              relavg = relavg + dabs(relerr)
              nruns = nruns + 1
              IF (dabs(relerr) .GT. dabs(relmax)) THEN
                  relmax = relerr
                  nqmax = nquad
              END IF
              IF (dabs(relerr) .LT. dabs(relmin)) THEN
                  relmin = relerr
                  nqmin = nquad
              END IF
          END IF
  400 CONTINUE
      GO TO 100
*
  500 IF (nruns .GT. 0) THEN
          relavg = relavg / dfloat(nruns)
          WRITE (stdout, '()')
          WRITE (stdout, 90000) 'Maximum', relmax, relmax/ulp,
     X        derbit(relmax,ulp), 'at nquad =', nqmax
          WRITE (stdout, 90000) 'Minimum', relmin, relmin/ulp,
     X        derbit(relmin,ulp), 'at nquad =', nqmin
          WRITE (stdout, 90000) 'Average', relavg, relavg/ulp,
     X        derbit(relavg,ulp)
      END IF
*
      WRITE (stdout, 60000)
      CLOSE (UNIT=stddat)
*
10000 FORMAT (/, 2X, 3X, 'i', 6X, 'x(i)', 10X, 7X, 'w(i)', 14X,
     X    7X, 'y(i)', 15X, 6X, 'z(i)')
20000 FORMAT (2X, I5, 1P, 1X, E23.15, 1X, E23.15, 1X, E10.2, 2X, 0P,
     X    1X, F9.2, 1X, F9.2)
25000 FORMAT (2X, I5, 1P, 1X, E23.15, 1X, E23.15, 1X, E10.2, 2X,
     X    1X, E9.2, 0P, 1X, F9.2)
30000 FORMAT (/, 2X, 4X, 'p', 5X, 'Quadrature Result', 3X,
     X    4X, 'Exact Integral', 6X, 'Rel. Error', 1X, 'RelE (ULPs)',
     X    2X, 'Err (bits)')
40000 FORMAT (2X, I4, 1X, F20.15, 1P, 1X, E24.15, 1X, E24.15, 1X,
     X    E24.15)
50000 FORMAT (/, 2X, 'alpha = ', F10.6, 2X, 'beta = ', F10.6,
     X    5X, 'nquad = ', I6)
60000 FORMAT (/, 2X, 'Done')
70000 FORMAT (I1, 1X, I4, 1X, 1P, 1X, E26.18, 1X, E26.18, 1X, I10, 1X,
     X    E27.18, 1X, E9.2, 1X, E9.2)
80000 FORMAT (/, 2X, 'ERROR: ', A, ' returns ierr = ', I10)
90000 FORMAT (2X, A, ' relative error = ', 1P, E10.2, 1X,
     X    0P, F10.2, ' ulps', 1X, F10.2, ' bits', 1X, A, I5)
*
      END
