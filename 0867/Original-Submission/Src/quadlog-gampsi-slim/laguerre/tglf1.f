      PROGRAM tglf1
************************************************************************
*     (Test Gauss-Laguerre Quadrature with Function values)
*
*     This is the driver program for subroutine glqf() which returns
*     the nodes and weights necessary for the log-Laguerre quadrature.
*
*     The test functions are x^p, p = 0 to pmax, with weight
*
*             x^alpha * exp(-x) * log(x),      (alpha > -1)
*
*     integrated from 0 to infinity.
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
      EXTERNAL            dpsi,        dnan,        dvsum
*
      DOUBLE PRECISION    dabs,        deps,        derbit,      dfloat
      DOUBLE PRECISION    dgamma,      dpsi,        dnan,        dvsum
*
*     Statement functions
*
      DOUBLE PRECISION    f
*
*     Parameter variables
*
      CHARACTER*(*)       METHOD
      PARAMETER           (METHOD = 'glqf()')
*
      DOUBLE PRECISION    ONE
      PARAMETER           (ONE = 1.0d+00)
*
      DOUBLE PRECISION    ZERO
      PARAMETER           (ZERO = 0.0d+00)
*
      INCLUDE 'stdio.inc'
      INCLUDE 'maxpts.inc'
*
*     Local variables
*
      DOUBLE PRECISION    alpha,       b,           c,           exact
      DOUBLE PRECISION    fn,          relavg,      relerr,      relmax
      DOUBLE PRECISION    relmin,      relulp,      result,      ulp
      DOUBLE PRECISION    v(2,MAXPTS), w(MAXPTS)
      DOUBLE PRECISION    wxm1(MAXPTS),             x(MAXPTS),   xarg
      DOUBLE PRECISION    y(MAXPTS),   z(MAXPTS),   zinit(1)
*
      INTEGER             i,           ierr,        narg,        neval
      INTEGER             nqmax,       nqmin,       nquad,       nruns
      INTEGER             p,           pmax
*
      f(xarg,narg)  = xarg**narg
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
      CALL dcopy(MAXPTS, zinit, 0, w, 1)
      CALL dcopy(MAXPTS, zinit, 0, wxm1, 1)
      CALL dcopy(MAXPTS, zinit, 0, x, 1)
      CALL dcopy(MAXPTS, zinit, 0, y, 1)
      CALL dcopy(MAXPTS, zinit, 0, z, 1)
      alpha = zinit(1)
      b = zinit(1)
      c = zinit(1)
      exact = zinit(1)
      fn = zinit(1)
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
      WRITE (stddat,'(A)') '### int(x^p * x^alpha * exp(-x) * ln(x),' //
     X    ' x = 0..infinity)'
      WRITE (stddat,'(A)') '###'
      WRITE (stddat,'(A)')
     X    '### Line format: 2 p alpha neval result relerr abserr'
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
      CALL prthdr(stdout,
     X     'int(x^p * x^alpha * exp(-x) * ln(x), x = 0..infinity)')
*
*     Loop reading input data:
*
  100 READ (stdin, *, END=500, ERR=500) nquad, alpha, pmax
      WRITE (stdout, 50000) alpha, nquad
*
      CALL glqf(x, w, wxm1, y, z, alpha, nquad, ierr)
      IF (ierr .NE. 0) THEN
          WRITE (stderr, 80000) METHOD, ierr
          GO TO 100
      END IF
*
      WRITE (stdout, 10000)
      DO 200 i = 1, nquad
          WRITE (stdout, 20000) i, x(i), w(i), wxm1(i), y(i), z(i)
  200 CONTINUE
*
*     Inside the loop, we have these invariants:
*
*         b = Gamma(p + alpha + 1)
*         c = psi(p + alpha + 1)
*         exact = b*c
*         f(x) = x**p
*
*     Instead of invoking dgamma() and dpsi() each time, the following
*     recursion relations are used to maintain those loop invariants:
*
*         Gamma(x+1) = x * Gamma(x)
*         psi(x+1) = psi(x) + 1/x
*
*     where x is called fn.  Unfortunately, when x is close to -1, the
*     Gamma and psi values cannot be determined accurately, and the
*     errors then propagate to higher p values.  We therefore revert
*     to direct evaluation.
*
      b = dgamma(ONE + alpha)
      c = dpsi(ONE + alpha)
      WRITE (stdout, 30000)
      DO 400 p = 0, pmax
*
*         Evaluate test integral as quadrature sum:
*
          DO 300 i = 1, nquad
              v(1,i) = wxm1(i) * f(x(i),p)
              v(2,i) = -z(i) * f(y(i),p)
  300     CONTINUE
          neval = 2*nquad
          result = dvsum(v, 2 * nquad)
*
*         Evaluate test integral analytically:
*
          fn = dfloat(p + 1) + alpha
*         exact = b*c
          exact = dgamma(fn) * dpsi(fn)
*
*         Compute the relative error as a raw number, a multiple of
*         one ulp, and as a number of bits lost:
*
          relerr = (exact - result)/exact
          relulp = relerr/ulp
          IF (dabs(relulp) .LT. 10000.0d+00) THEN
              WRITE (stdout, 40000) p, result, exact, relerr, relulp,
     X            derbit(relerr,ulp)
          ELSE
              WRITE (stdout, 45000) p, result, exact, relerr, relulp,
     X            derbit(relerr,ulp)
          END IF
          WRITE (stddat, 70000) 2, p, (fn - ONE), neval, result,
     X        relerr, (exact - result)
          b = b*fn
          c = c + ONE/fn
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
     X    7X, 'wxm1(i)', 13X, 5X, 'y(i)', 14X, 7X, 'z(i)')
20000 FORMAT (2X, I4, 1X, F20.15, 1P, 1X, E24.15, 1X, E24.15, 1X,
     X    E24.15, 1X, E24.15)
30000 FORMAT (/, 2X, 4X, 'p', 5X, 'Quadrature Result', 3X,
     X    4X, 'Exact Integral', 6X, 'Rel. Error', 1X, 'RelE (ULPs)',
     X    2X, 'Err (bits)')
40000 FORMAT (2X, I5, 1P, 1X, E23.15, 1X, E23.15, 1X, E10.2, 2X, 0P, 1X,
     X    F9.2, 1X, F9.2)
45000 FORMAT (2X, I5, 1P, 1X, E23.15, 1X, E23.15, 1X, E10.2, 2X, 1X,
     X    E9.2, 0P, 1X, F9.2)
50000 FORMAT (/, 2X, 'alpha = ', F10.6, '  nquad = ', I6)
60000 FORMAT (/, 2X, 'Done')
70000 FORMAT (I1, 1X, I3, 1X, 1P, E27.18, 1X, I10, 1X, E27.18, 1X, E9.2,
     X    1X, E9.2)
80000 FORMAT (/, 2X, 'ERROR: ', A, ' returns ierr = ', I10)
90000 FORMAT (2X, A, ' relative error = ', 1P, E10.2, 1X,
     X    0P, F10.2, ' ulps', 1X, F10.2, ' bits', 1X, A, I5)
*
      END
