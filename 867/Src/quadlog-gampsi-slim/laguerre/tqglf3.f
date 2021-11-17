      PROGRAM tqglf3
************************************************************************
*     (Test Gauss-Laguerre Quadrature with Function values)
*
*     This is a driver program for subroutine qglqf() which returns the
*     nodes and weights necessary for the log-Laguerre quadrature.
*
*     The test functions are x^(p - alpha) * exp((1-mu)*x), p = 0 to
*     pmax, with weight
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
      EXTERNAL            qeps,        qerbit,      qnan,        qvsum
*
      REAL*16             qabs,        qeps,        qerbit,      qlog
      REAL*16             qnan,        qvsum
*
*     Parameter variables
*
      CHARACTER*(*)       METHOD
      PARAMETER           (METHOD = 'qglqf()')
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
      INCLUDE 'maxpts.inc'
      INCLUDE 'stdio.inc'
*
*     Local variables
*
      REAL*16             alpha,       exact,       factn,       mu
      REAL*16             mumax,       mumin,       munp2,       qp
      REAL*16             relavg,      relerr,      relmax,      relmin
      REAL*16             relulp,      result,      ulp
      REAL*16             v(2,MAXPTS), w(MAXPTS),   wxm1(MAXPTS)
      REAL*16             x(MAXPTS),   y(MAXPTS),   z(MAXPTS)
      REAL*16             zinit(1)
*
      INTEGER             i,           ierr,        neval,       nqmax
      INTEGER             nqmin,       nquad,       nruns,       p
      INTEGER             pmax,        ppmax,       ppmin
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
      CALL qcopy (MAXPTS, zinit, 0, wxm1, 1)
      CALL qcopy (MAXPTS, zinit, 0, x, 1)
      CALL qcopy (MAXPTS, zinit, 0, y, 1)
      CALL qcopy (MAXPTS, zinit, 0, z, 1)
      alpha = zinit(1)
      exact = zinit(1)
      factn = zinit(1)
      mu = zinit(1)
      mumax = zinit(1)
      mumin = zinit(1)
      munp2 = zinit(1)
      qp = zinit(1)
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
      WRITE (stddat,'(A)') '### int(x^p * exp(-mu*x)' //
     X     ' * ln(x), x = 0..infinity)'
      WRITE (stddat,'(A)') '###'
      WRITE (stddat,'(A)')
     X    '### Line format: 2 p mu neval result relerr abserr'
      WRITE (stddat,'(A)') '###'
*
      ulp = qeps(ONE)
      mumin = 1.0q+75
      mumax = ZERO
      ppmax = 0
      ppmin = 0
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
     X     'int(x^p * exp(-mu*x) * ln(x), x = 0..infinity)')
*
*     Loop reading input data:
*
  100 READ (stdin, *, END=500, ERR=500) nquad, mu, pmax
      IF (mu .LE. ZERO) THEN
            WRITE (stderr, 80000) mu
            GO TO 100
      END IF
*
*     Loop entry invariants:
*         qp = Q_p
*            = \int_0^{\infty} x^p exp(-mu*x) ln(x) dx
*            = (mu)**(-(p+1)) p! [1 + 1/2 + 1/3 + ... + 1/p - EULERC -
*                                 ln(mu)]
*     This can be updated by recursion:
*         Q_(p+1) = mu**(-1)*(p + 1)*Q_p + (mu)**(-(p+2))p!
*
      factn = 1
      munp2 = mu**2
      qp = -(EULERC + log(mu))/mu
      WRITE (stdout, 40000) nquad, mu, pmax
      WRITE (stdout, 10000)
      DO 400 p = 0, pmax
          alpha = qfloat(p)
          CALL qglqf(x, w, wxm1, y, z, alpha, nquad, ierr)
          IF (ierr .NE. 0) THEN
              WRITE (stderr, 70000) METHOD, ierr
              GO TO 400
          END IF
*
*         Evaluate test integral as quadrature sum:
*
          DO 200 i = 1, nquad
              v(1,i) = wxm1(i)
              v(2,i) = -z(i)
  200     CONTINUE
          neval = 0
          result = qvsum(v, 2 * nquad)
          result = (result - qlog(mu)*factn)*mu/munp2
*
*         Analytic calculation of test integral
*
          exact = qp
          relerr = (exact - result)/exact
          relulp = relerr/ulp
          IF (qabs(relulp) .LT. 10000.0q+00) THEN
              WRITE (stdout, 20000) p, result, exact, relerr, relulp,
     X            qerbit(relerr,ulp)
          ELSE
              WRITE (stdout, 30000) p, result, exact, relerr, relulp,
     X            qerbit(relerr,ulp)
          END IF
          WRITE (stddat, 60000) 2, p, mu, neval, result, relerr,
     X        (exact - result)
          if (relerr .eq. relerr) THEN
*
*             Track the average, minimum, and maximum relative error as
*             long as it is not a NaN.
*
              relavg = relavg + qabs(relerr)
              nruns = nruns + 1
              IF (qabs(relerr) .GT. qabs(relmax)) THEN
                  relmax = relerr
                  ppmax = p
                  nqmax = nquad
                  mumax = mu
              END IF
              IF (qabs(relerr) .LT. qabs(relmin)) THEN
                  ppmin = p
                  relmin = relerr
                  nqmin = nquad
                  mumin = mu
              END IF
          END IF
*
*         Update loop variables to maintain loop entry invariant:
*
          qp = qfloat(p + 1)*qp/mu + factn/munp2
          factn = qfloat(p + 1) * factn
          munp2 = mu * munp2
  400 CONTINUE
      GO TO 100
*
  500 IF (nruns .GT. 0) THEN
          relavg = relavg / qfloat(nruns)
          WRITE (stdout, '()')
          WRITE (stdout, 90000) 'Maximum', relmax, relmax/ulp,
     X        qerbit(relmax,ulp), 'at nquad =', nqmax, 'mu =', mumax,
     X        'p = ',ppmax
          WRITE (stdout, 90000) 'Minimum', relmin, relmin/ulp,
     X        qerbit(relmin,ulp), 'at nquad =', nqmin, 'mu =', mumin,
     X        'p = ',ppmin
          WRITE (stdout, 90000) 'Average', relavg, relavg/ulp,
     X        qerbit(relavg,ulp)
      END IF
      WRITE (stdout, 50000)
      CLOSE (UNIT=stddat)
*
10000 FORMAT (/, 2X, 4X, 'p', 5X, 'Quadrature Result', 3X,
     X    4X, 'Exact Integral', 6X, 'Rel. Error', 1X, 'RelE (ULPs)',
     X    2X, 'Err (bits)')
20000 FORMAT (2X, I5, 1P, 1X, E23.15, 1X, E23.15, 1X, E10.2, 2X, 0P,
     X    1X, F9.2, 1X, F9.2)
30000 FORMAT (2X, I5, 1P, 1X, E23.15, 1X, E23.15, 1X, E10.2, 2X,
     X    1X, E9.2, 0P, 1X, F9.2)
40000 FORMAT (/, 2X, 'nquad = ', I6, '  mu = ', F10.6, '  pmax = ', I6)
50000 FORMAT (/, 2X, 'Done')
60000 FORMAT (I1, 1X, I3, 1X, 1P, E27.18, 1X, I10, 1X, E27.18, 1X, E9.2,
     X    1X, E9.2)
70000 FORMAT (/, 2X, 'ERROR: ', A, ' returns ierr = ', I10)
80000 FORMAT (/, 2X, 'Illegal value of mu:', 1P, E11.2, ' <= 0')
90000 FORMAT (2X, A, ' relative error = ', 1P, E10.2, 1X,
     X    0P, F10.2, ' ulps', 1X, F10.2, ' bits', 1X, A, I5, 1X, A,
     X    1P, E15.7, 1X, A, I5)
*
      END
