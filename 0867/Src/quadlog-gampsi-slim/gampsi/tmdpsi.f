      PROGRAM tmdpsi
************************************************************************
*     (Monotonicity check)
*     Do approximate monotonicity checks for dpsi(x) at boundaries
*     where the approximation changes.
*     [12-Jul-2000]
************************************************************************
*
*     External functions
*
      EXTERNAL            dinf,        dpsi
*
      DOUBLE PRECISION    dinf,        dpsi
*
*     Parameter variables
*
      DOUBLE PRECISION    beta
      PARAMETER           (beta = 2.0d+00)
*
      DOUBLE PRECISION    one
      PARAMETER           (one = 1.0d0)
*
      INTEGER             pmono
      PARAMETER           (pmono = 52)
*
      DOUBLE PRECISION    delta
      PARAMETER           (delta = 1024.0d+00 * beta**(-pmono))
*
      DOUBLE PRECISION    smin
      PARAMETER           (smin = one - delta)
*
      DOUBLE PRECISION    smax
      PARAMETER           (smax = one + delta)
*
      INCLUDE 'dxmin.inc'
*
      DOUBLE PRECISION    xmin1
      PARAMETER           (xmin1 = xmin)
*
      DOUBLE PRECISION    xlarge
      PARAMETER           (xlarge = 2.71d+14)
*
      DOUBLE PRECISION    xsmall
      PARAMETER           (xsmall = 5.80d-09)
*
*     CALL monchk (dpsi, 'dpsi', -2.0d+00, -1.0d+00)
*     CALL monchk (dpsi, 'dpsi', 1.4617d+00, dinf())
      CALL monchk (dpsi, 'dpsi', -xmin1*smax, -xmin1*smin, beta,pmono)
      CALL monchk (dpsi, 'dpsi', xmin1*smin, xmin1*smax, beta,pmono)
      CALL monchk (dpsi, 'dpsi', xsmall*smin, xsmall*smax, beta,pmono)
      CALL monchk (dpsi, 'dpsi', 0.5d+00*smin, 0.5d+00*smax, beta,pmono)
      CALL monchk (dpsi, 'dpsi', 3.0d+00*smin, 3.0d+00*smax, beta,pmono)
      CALL monchk (dpsi, 'dpsi', xlarge*smin, xlarge*smax, beta,pmono)
*
      END


      SUBROUTINE monchk (f, name, xmin, xmax, beta, pmono)
************************************************************************
*     (Monotonicity check)
*     Do an approximate monotonicity check for f(x) in (xmin,xmax).  We
*     cannot step by machine epsilon, because that would take too
*     long, so instead, we step x by 1 + beta**(-pmono), where pmono
*     is chosen to give reasonable running times.
*     [12-Jul-2000]
************************************************************************
*
*     External functions
*
      EXTERNAL            f
*
      DOUBLE PRECISION    f
*
      INCLUDE 'stdio.inc'
*
*     Argument variables
*
      CHARACTER*(*)       name
      DOUBLE PRECISION    beta
      INTEGER             pmono
*
      DOUBLE PRECISION    xmax,        xmin
*
*     Parameter variables
*
      DOUBLE PRECISION    one
      PARAMETER           (one = 1.0d+00)
*
      DOUBLE PRECISION    zero
      PARAMETER           (zero = 0.0d+00)
*
*     Local variables
*
      DOUBLE PRECISION    flast,       fthis,       x,           xlast
      DOUBLE PRECISION    xmult
*
      INTEGER             ncheck
*
      ncheck = 0
      IF (xmin .LT. zero) THEN
          xmult = one/(one + beta**(-pmono))
      ELSE
          xmult = one + beta**(-pmono)
      END IF
      x = xmin * xmult
      xlast = x
      flast = f(xlast)
      x = x * xmult
      WRITE (stdout,10000) idint(beta), -pmono
   10 IF (x .LT. xmax) THEN
          ncheck = ncheck + 1
          fthis = f(x)
          IF (fthis .LT. flast) THEN
              WRITE (stdout,20000) name, xlast, flast, name, x, fthis
          END IF
          xlast = x
          flast = fthis
          x = x * xmult
          GO TO 10
      END IF
      WRITE (stdout,30000) name, ncheck, xmin, xmax
10000 FORMAT ('Stepsize = ', i2, '**(',i3,')')
20000 FORMAT ('Monotonicity failure:'/
     X      1p, 4x, a,'(', e27.17e4, ') = ', e27.17e4/
     X      1p, 4x, a,'(', e27.17e4, ') = ', e27.17e4)
30000 FORMAT (a, '(x): completed ', I12, ' monotonicity checks on (',
     X      1p, e27.17e4, ',', e27.17e4, ')')
      END
