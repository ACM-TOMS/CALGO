      PROGRAM tqglf4
************************************************************************
*     (Test Gauss-Laguerre Quadrature with Function values)
*
*     This is a driver program for subroutine qglqf() which returns the
*     nodes and weights necessary for the log-Laguerre quadrature.
*
*     The test function is sin(sigma*x) with weight
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
      EXTERNAL            qeps,        qerbit,      qgamma,      qnan
      EXTERNAL            qpsi,        qvsum
*
      REAL*16             qabs,        qatan,       qeps,        qerbit
      REAL*16             qgamma,      qlog,        qpsi,        qsin
      REAL*16             qtan,        qnan,        qvsum
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
      REAL*16             HALF
      PARAMETER           (HALF = 0.5q+00)
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
      REAL*16             alfmax,      alfmin,      alpha,       dsda
      REAL*16             exact,       onepa,       oneps2,      relavg
      REAL*16             relerr,      relmax,      relmin,      relulp
      REAL*16             result,      s,           sigarg,      sigma
      REAL*16             sigmax,      sigmin,      theta,       ulp
      REAL*16             v(2,MAXPTS), w(MAXPTS)
      REAL*16             wxm1(MAXPTS),             x(MAXPTS),   xarg
      REAL*16             y(MAXPTS),   z(MAXPTS),   zinit(1)
*
      INTEGER             i,           ierr,        neval,       nqmax
      INTEGER             nqmin,       nquad,       nruns
*
      f(xarg,sigarg) = qsin(sigarg*xarg)
*     fprime(xarg,sigarg) = sigarg * qcos(sigarg*xarg)
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
      alfmax = zinit(1)
      alfmin = zinit(1)
      alpha = zinit(1)
      exact = zinit(1)
      oneps2 = zinit(1)
      relavg = zinit(1)
      relerr = zinit(1)
      relmax = zinit(1)
      relmin = zinit(1)
      relulp = zinit(1)
      result = zinit(1)
      sigarg = zinit(1)
      sigma = zinit(1)
      sigmax = zinit(1)
      sigmin = zinit(1)
      ulp = zinit(1)
      xarg = zinit(1)
*
      OPEN (UNIT=stddat, FILE='output.dat', STATUS='unknown',
     X    FORM='formatted')
*
      WRITE (stddat,'(A)') '### Numerical integration with ' // METHOD
      WRITE (stddat,'(A)') '###'
      WRITE (stddat,'(A)') '### int(sin(sigma*x) * x^alpha * ' //
     X     'exp(-x) * ln(x), x = 0..infinity)'
      WRITE (stddat,'(A)') '###'
      WRITE (stddat,'(A)')
     X    '### Line format: 2 sigma alpha neval result relerr abserr'
      WRITE (stddat,'(A)') '###'
*
      ulp = qeps(ONE)
      sigmin = 1.0q+75
      sigmax = ZERO
      nqmax = 0
      nqmin = 0
      nruns = 0
      relavg = ZERO
      relmax = ZERO
      relmin = 1.0q+75
*
*     Write a standard output header identifying the integral, and the
*     host precision
*
      CALL qprthd(stdout,
     X     'int(sin(sigma*x) * x^alpha * exp(-x) * ln(x), ' //
     X     'x = 0..infinity)')
      WRITE (stdout, 10000)
*
*     Loop reading input data:
*
  100 READ (stdin, *, END=300, ERR=300) nquad, sigma, alpha
      CALL qglqf(x, w, wxm1, y, z, alpha, nquad, ierr)
      IF (ierr .NE. 0) THEN
          WRITE (stderr, 60000) METHOD, ierr
          GO TO 100
      END IF
*
*     Evaluate test integral as quadrature sum:
*
      DO 200 i = 1, nquad
          v(1,i) = wxm1(i) * f(x(i),sigma)
          v(2,i) = -z(i) * f(y(i),sigma)
  200 CONTINUE
      neval = 2*nquad
      result = qvsum(v, 2 * nquad)
*
*     Analytic calculation of test integral.  We start with the known
*     integral (Gradshteyn and Ryzhik, 4th edition, 3.944.5, p. 490):
*
*     S = \int_{0}^{\infty} dx\, x^\alpha e^{-x} \sin(\sigma x)
*       = \Gamma(1 + \alpha) (1 + \sigma^2)^{-(1 + \alpha)/2}
*         \sin((1 + \alpha)\theta)
*
*     where
*
*     \theta = arctan(\sigma)
*
*     and differentiate with respect to \alpha to obtain
*
*     (dS/d\alpha) = \int_{0}^{\infty} dx\,
*                          x^\alpha e^{-x} \ln(x) \sin(\sigma x)
*                  = S[\psi(1 + \alpha)
*                      + \theta / \tan((1 + \alpha)\theta)
*                      - (1/2)\ln(1 + \sigma^2)]
*
      onepa = ONE + alpha
      oneps2 = ONE + sigma**2
      theta = qatan(qabs(sigma))
      s = qgamma(onepa) * qsin(onepa * theta) / oneps2**(HALF * onepa)
      IF (sigma .LT. ZERO) s = -s
      dsda = s*(qpsi(onepa) + theta/qtan(onepa * theta) -
     X    HALF*qlog(oneps2))
      exact = dsda
      relerr = (exact - result)/exact
      relulp = relerr/ulp
*
*     For low quadrature order, the relative error in ulps can be
*     quite large, so we switch from Fw.d to Ew.d output format.
*
      IF (qabs(relulp) .LT. 10000.0q+00) THEN
          WRITE (stdout, 20000) nquad, sigma, alpha, result, exact,
     X         relerr, relulp, qerbit(relerr,ulp)
      ELSE
          WRITE (stdout, 30000) nquad, sigma, alpha, result, exact,
     X         relerr, relulp, qerbit(relerr,ulp)
      END IF
      WRITE (stddat, 50000) 2, sigma, alpha, neval, result, relerr,
     X        (exact - result)
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
              sigmax = sigma
              alfmax = alpha
          END IF
          IF (qabs(relerr) .LT. qabs(relmin)) THEN
              relmin = relerr
              nqmin = nquad
              sigmin = sigma
              alfmin = alpha
          END IF
      END IF
      GO TO 100
*
  300 IF (nruns .GT. 0) THEN
          relavg = relavg / qfloat(nruns)
          WRITE (stdout, '()')
          WRITE (stdout, 70000) 'Maximum', relmax, relmax/ulp,
     X        qerbit(relmax,ulp), 'at nquad =', nqmax,
     X        'sigma =', sigmax, 'alpha =', alfmax
          WRITE (stdout, 70000) 'Minimum', relmin, relmin/ulp,
     X        qerbit(relmin,ulp), 'at nquad =', nqmin,
     X        'sigma =', sigmin, 'alpha =', alfmin
          WRITE (stdout, 70000) 'Average', relavg, relavg/ulp,
     X        qerbit(relavg,ulp)
      END IF
      WRITE (stdout, 40000)
      CLOSE (UNIT=stddat)
*
10000 FORMAT (/,2X, 1X, 'nquad', 3X, 2X, 'sigma', 3X, 2X, 'alpha', 3X,
     X    2X, 'Quadrature Result', 5X,
     X    3X, 'Exact Integral', 4X,
     X    1X, 'Rel. Error', 1X, 'RelE (ULPs)', 2X, 'Err (bits)')
20000 FORMAT (2X, I6, 1X, F9.4, 1X, F9.4, 1P, 1X, E23.15, 1X, E23.15,
     X    1X, E10.2, 2X, 0P, 1X, F9.2, 1X, F9.2)
30000 FORMAT (2X, I6, 1X, F9.4, 1X, F9.4, 1P, 1X, E23.15, 1X, E23.15,
     X    1X, E10.2, 2X, E10.2, 0P, 1X, F9.2)
40000 FORMAT (/, 2X, 'Done')
50000 FORMAT (I1, 1X, 1X, 1P, 1X, E26.18, 1X, E26.18, 1X, I10, 1X,
     X    E27.18, 1X, E9.2, 1X, E9.2)
60000 FORMAT (/, 2X, 'ERROR: ', A, ' returns ierr = ', I10)
70000 FORMAT (2X, A, ' relative error = ', 1P, E10.2, 1X,
     X    1P, E10.2, ' ulps', 1X, 0P, F10.2, ' bits', 1X, A, I5,
     X    1X, A, F10.4, 1X, A, F10.4)
*
      END
