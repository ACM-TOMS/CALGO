      PROGRAM tmpsi
************************************************************************
*     (Monotonicity check)
*     Do approximate monotonicity checks for psi(x) at boundaries
*     where the approximation changes.
*     [12-Jul-2000]
************************************************************************
*
*     External functions
*
      EXTERNAL            ainf,        psi
*
      REAL                ainf,        psi
*
*     Parameter variables
*
      INTEGER             pmono
      PARAMETER           (pmono = 23)
*
      REAL                beta
      PARAMETER           (beta = 2.0)
*
      REAL                cutoff
      PARAMETER           (cutoff = 13.0)
*
      REAL                delta
      PARAMETER           (delta = 1024.0 * beta**(-pmono))
*
      REAL                one
      PARAMETER           (one = 1.0)
*
      REAL                smin
      PARAMETER           (smin = one - delta)
*
      REAL                smax
      PARAMETER           (smax = one + delta)
*
      INCLUDE 'axmin.inc'
*
      REAL                xmin1
      PARAMETER           (xmin1 = xmin)
*
      INCLUDE 'axlarg.inc'
*
      INCLUDE 'axsmal.inc'
*
*     CALL monchk (psi, 'psi', -2.0d+00, -1.0d+00)
*     CALL monchk (psi, 'psi', 1.4617d+00, ainf())
      CALL monchk (psi, 'psi', -xmin1*smax, -xmin1*smin, beta, pmono)
      CALL monchk (psi, 'psi', xmin1*smin, xmin1*smax, beta, pmono)
      CALL monchk (psi, 'psi', xsmall*smin, xsmall*smax, beta, pmono)
      CALL monchk (psi, 'psi', 0.5*smin, 0.5*smax, beta, pmono)
      CALL monchk (psi, 'psi', 3.0*smin, 3.0*smax, beta, pmono)
      CALL monchk (psi, 'psi', cutoff*smin, cutoff*smax, beta, pmono)
      CALL monchk (psi, 'psi', xlarge*smin, xlarge*smax, beta, pmono)
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
      INTEGER             int
*
      REAL                f
*
*     Parameter variables
*
      INCLUDE 'stdio.inc'
*
      REAL                one
      PARAMETER           (one = 1.0)
*
      REAL                zero
      PARAMETER           (zero = 0.0)
*
*     Argument variables
*
      CHARACTER*(*)       name
*
      INTEGER             pmono
*
      REAL                beta,        xmax,        xmin
*
*     Local variables
*
      INTEGER             ncheck
*
      REAL                flast,       fthis,       x,           xlast
      REAL                xmult
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
      WRITE (stdout,10000) int(beta), -pmono
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
     X      1p, 4x, a,'(', e17.7e3, ') = ', e17.7e3/
     X      1p, 4x, a,'(', e17.7e3, ') = ', e17.7e3)
30000 FORMAT (a, '(x): completed ', I12, ' monotonicity checks on (',
     X      1p, e17.7e3, ',', e17.7e3, ')')
      END
