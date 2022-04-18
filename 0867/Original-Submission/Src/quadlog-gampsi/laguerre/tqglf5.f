      PROGRAM tqglf5
************************************************************************
*     (Test Gauss-Laguerre Quadrature with Function values)
*
*     This is a driver program for subroutine qglqf() which returns the
*     nodes and weights necessary for the ordinary Laguerre quadrature.
*
*     The test functions are x^p, p = 0 to pmax, with weight
*
*             x^alpha * exp(-x),      (alpha > -1)
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
      PARAMETER           (METHOD = 'qglqf()')
*
      REAL*16             ONE
      PARAMETER           (ONE = 1.0q+00)
*
      REAL*16             ZERO
      PARAMETER           (ZERO = 0.0q+00)
*
      INCLUDE 'stdio.inc'
      INCLUDE 'maxpts.inc'
*
*     Local variables
*
      REAL*16             alpha,       exact,       relavg,      relerr
      REAL*16             relmax,      relmin,      relulp,      result
      REAL*16             ulp
      REAL*16             v(MAXPTS),   w(MAXPTS),   wxm1(MAXPTS)
      REAL*16             x(MAXPTS),   xarg,        y(MAXPTS)
      REAL*16             z(MAXPTS),   zinit(1)
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
      zinit(1) = qnan()
      CALL qcopy(MAXPTS, zinit, 0, w, 1)
      CALL qcopy(MAXPTS, zinit, 0, wxm1, 1)
      CALL qcopy(MAXPTS, zinit, 0, x, 1)
      CALL qcopy(MAXPTS, zinit, 0, y, 1)
      CALL qcopy(MAXPTS, zinit, 0, z, 1)
      alpha = zinit(1)
      exact = zinit(1)
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
      WRITE (stddat,'(A)') '### int(x^p * x^alpha * exp(-x),' //
     X    ' x = 0..infinity)'
      WRITE (stddat,'(A)') '###'
      WRITE (stddat,'(A)')
     X    '### Line format: 2 p alpha neval result relerr abserr'
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
     X     'int(x^p * x^alpha * exp(-x), x = 0..infinity)')
*
*     Loop reading input data:
*
  100 READ (stdin, *, END=500, ERR=500) nquad, alpha, pmax
      WRITE (stdout, 60000) alpha, nquad
*
      CALL qglqf(x, w, wxm1, y, z, alpha, nquad, ierr)
      IF (ierr .NE. 0) THEN
          WRITE (stderr, 90000) METHOD, ierr
          GO TO 100
      END IF
*
      WRITE (stdout, 10000)
      DO 200 i = 1, nquad
          WRITE (stdout, 20000) i, x(i), w(i), wxm1(i), y(i), z(i)
  200 CONTINUE
*
      WRITE (stdout, 30000)
      DO 400 p = 0, pmax
*
*         Evaluate test integral as quadrature sum:
*
          DO 300 i = 1, nquad
              v(i) = w(i) * f(x(i),p)
  300     CONTINUE
          neval = nquad
          result = qvsum(v, nquad)
*
*         Evaluate test integral analytically:
*
          exact = qgamma(alpha + qfloat(p + 1))
*
*         Compute the relative error as a raw number, a multiple of
*         one ulp, and as a number of bits lost:
*
          relerr = (exact - result)/exact
          relulp = relerr/ulp
          IF (qabs(relulp) .LT. 10000.0q+00) THEN
              WRITE (stdout, 40000) p, result, exact, relerr, relulp,
     X            qerbit(relerr,ulp)
          ELSE
              WRITE (stdout, 50000) p, result, exact, relerr, relulp,
     X            qerbit(relerr,ulp)
          END IF
          WRITE (stddat, 80000) 2, p, alpha, neval, result,
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
          WRITE (stdout, 95000) 'Maximum', relmax, relmax/ulp,
     X        qerbit(relmax,ulp), 'at nquad =', nqmax
          WRITE (stdout, 95000) 'Minimum', relmin, relmin/ulp,
     X        qerbit(relmin,ulp), 'at nquad =', nqmin
          WRITE (stdout, 95000) 'Average', relavg, relavg/ulp,
     X        qerbit(relavg,ulp)
      END IF
*
      WRITE (stdout, 70000)
      CLOSE (UNIT=stddat)
*
10000 FORMAT (/, 2X, 3X, 'i', 6X, 'x(i)', 10X, 8X, 'w(i)', 14X,
     X    8X, 'wxm1(i)', 13X, 6X, 'y(i)', 14X, 8X, 'z(i)')
20000 FORMAT (2X, I4, 1X, F20.15, 1P, 1X, E25.15, 1X, E25.15, 1X,
     X    E25.15, 1X, E25.15)
30000 FORMAT (/, 2X, 4X, 'p', 5X, 'Quadrature Result', 3X,
     X    4X, 'Exact Integral', 6X, 'Rel. Error', 1X, 'RelE (ULPs)',
     X    2X, 'Err (bits)')
40000 FORMAT (2X, I5, 1P, 1X, E23.15, 1X, E23.15, 1X, E10.2, 2X, 0P,
     X    1X, F9.2, 1X, F9.2)
50000 FORMAT (2X, I5, 1P, 1X, E23.15, 1X, E23.15, 1X, E10.2, 2X,
     X    1X, E9.2, 0P, 1X, F9.2)
60000 FORMAT (/, 2X, 'alpha = ', F10.6, '  nquad = ', I6)
70000 FORMAT (/, 2X, 'Done')
80000 FORMAT (I1, 1X, I3, 1X, 1P, E27.18, 1X, I10, 1X, E27.18,
     X    1X, E9.2, 1X, E9.2)
90000 FORMAT (/, 2X, 'ERROR: ', A, ' returns ierr = ', I10)
95000 FORMAT (2X, A, ' relative error = ', 1P, E10.2, 1X,
     X    0P, F10.2, ' ulps', 1X, F10.2, ' bits', 1X, A, I5)
*
      END
