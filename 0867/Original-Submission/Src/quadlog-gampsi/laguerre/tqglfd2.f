      PROGRAM tqglfd2
************************************************************************
*     (Test Gauss-Laguerre Quadrature w/ Function and Derivative values)
*
*     This is the driver program for subroutine qglqfd() which returns
*     the nodes and weights necessary for the log-Laguerre quadrature.
*
*     The test function is f(x) = 1, with weight
*
*             x^alpha * exp(-x) * log(x),      (alpha > -1)
*
*     integrated from 0 to infinity.  The integral can be evaluated
*     exactly as
*
*         \int_0^\infty \exp(-t) \ln(t) dt = -\gamma
*
*     where gamma is the Euler-Mascheroni constant (sometimes called
*     C).  Thus, we should have from eq. (15):
*
*         \sum_{i=1}^{N}[(f'(\xia)\delta x_i + f(\xia)dW_i] = -\gamma
*
*     for any N.  In practice, in IEEE 754 double precision
*     arithmetic, overflows and NaNs appear in the computation
*     for N > 100.
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
      EXTERNAL            qeps,        qerbit,      qnan,        qvsum
*
      REAL*16             qabs,        qeps,        qerbit,      qnan
      REAL*16             qvsum
*
*     Parameter variables
*
      CHARACTER*(*)       METHOD
      PARAMETER           (METHOD = 'qglqfd()')
*
      REAL*16             EULERC
      PARAMETER           (EULERC =
     X    0.57721566490153286060651209008240243104215933593992q+00)
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
      REAL*16             deltax(MAXPTS),           deltaw(MAXPTS)
      REAL*16             exact,       relavg,      relerr,      relmax
      REAL*16             relmin,      relulp,      result,      ulp
      REAL*16             w(MAXPTS),   x(MAXPTS),   zinit(1)
*
      INTEGER             ierr,        neval,       nqmax,       nqmin
      INTEGER             nquad,       nruns
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
      CALL qcopy (MAXPTS, zinit, 0, deltax, 1)
      CALL qcopy (MAXPTS, zinit, 0, deltaw, 1)
      CALL qcopy (MAXPTS, zinit, 0, w, 1)
      CALL qcopy (MAXPTS, zinit, 0, x, 1)
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
      WRITE (stddat,'(A)') '### int(exp(-x) * ln(x),' //
     X    ' x = 0..infinity)'
      WRITE (stddat,'(A)') '###'
      WRITE (stddat,'(A)')
     X    '### Line format: 1 nquad neval result relerr abserr'
      WRITE (stddat,'(A)') '###'
      WRITE (stddat,'(A)') '### where nquad is the quadrature order'
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
     X     'int(exp(-x) * ln(x), x = 0..infinity) = -gamma')
      WRITE (stdout, 10000)
*
*     Loop reading input data:
*
  100 READ (stdin, *, END=300, ERR=300) nquad
*
      CALL qglqfd(x, w, deltaw, deltax, ZERO, nquad, ierr)
      IF (ierr .NE. 0) THEN
          WRITE (stderr, 50000) METHOD, ierr
          GO TO 100
      END IF
*
*     Evaluate test integral as quadrature sum:
*
      neval = 0
      result = qvsum(deltaw, nquad)
*
*     Analytic calculation of test integral
*
      exact = -EULERC
      relerr = (exact - result)/exact
      relulp = relerr/ulp
      IF (qabs(relulp) .LT. 10000.0q+00) THEN
          WRITE (stdout, 20000) nquad, result, exact, relerr, relulp,
     X        qerbit(relerr,ulp)
      ELSE
          WRITE (stdout, 25000) nquad, result, exact, relerr, relulp,
     X        qerbit(relerr,ulp)
      END IF
      WRITE (stddat, 30000) 1, nquad, neval, result, relerr,
     X    (exact - result)
      if (relerr .eq. relerr) THEN
*
*         Track the average, minimum, and maximum relative error as
*         long as it is not a NaN.
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
      GO TO 100
*
  300 IF (nruns .GT. 0) THEN
          relavg = relavg / qfloat(nruns)
          WRITE (stdout, '()')
          WRITE (stdout, 60000) 'Maximum', relmax, relmax/ulp,
     X        qerbit(relmax,ulp), 'at nquad =', nqmax
          WRITE (stdout, 60000) 'Minimum', relmin, relmin/ulp,
     X        qerbit(relmin,ulp), 'at nquad =', nqmin
          WRITE (stdout, 60000) 'Average', relavg, relavg/ulp,
     X        qerbit(relavg,ulp)
      END IF
      WRITE (stdout, 40000)
      CLOSE (UNIT=stddat)
*
10000 FORMAT (/, 2X, 'nquad', 15X, 'Quadrature Result', 13X,
     X    14X, 'Exact Integral', 16X, 'Rel. Error', 1X, 'RelE (ULPs)',
     X    2X, 'Err (bits)')
20000 FORMAT (2X, I5, 1P, 1X, E43.34, 1X, E43.34, 1X, E10.2, 2X, 0P,
     X    1X, F9.2, 1X, F9.2)
25000 FORMAT (2X, I5, 1P, 1X, E43.34, 1X, E43.34, 1X, E10.2, 2X,
     X    1X, E9.2, 0P, 1X, F9.2)
30000 FORMAT (I1, 1X, I3, 1X, I10, 1X, 1P, E43.34, 1X, E9.2, 1X, E9.2)
40000 FORMAT (/, 2X, 'Done')
50000 FORMAT (/, 2X, 'ERROR: ', A, ' returns ierr = ', I10)
60000 FORMAT (2X, A, ' relative error = ', 1P, E10.2, 1X,
     X    0P, F10.2, ' ulps', 1X, F10.2, ' bits', 1X, A, I5)
*
      END
