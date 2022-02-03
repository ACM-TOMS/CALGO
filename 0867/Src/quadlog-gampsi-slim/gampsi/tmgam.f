      PROGRAM tmgam
************************************************************************
*     (Monotonicity check)
*     Do approximate monotonicity checks for gamma(x) at boundaries
*     where the approximation changes.
*     [12-Jul-2000]
************************************************************************
*
*     External functions
*
      EXTERNAL            gamma
*
      REAL                gamma
*
*     Parameter variables
*
      INTEGER             pmono
      PARAMETER           (pmono = 23)
*
      REAL                beta
      PARAMETER           (beta = 2.0)
*
      REAL                one
      PARAMETER           (one = 1.0)
*
      REAL                two
      PARAMETER           (two = 2.0)
*
      REAL                delta
      PARAMETER           (delta = 1024.0 * beta**(-pmono))
*
      REAL                smin
      PARAMETER           (smin = one - delta)
*
      REAL                smax
      PARAMETER           (smax = one + delta)
*
      INCLUDE 'aeps.inc'
*
      CALL monchk (gamma,'gamma',eps*smin,eps*smax,beta,pmono,.false.)
      CALL monchk (gamma,'gamma',one*smin,one*smax,beta,pmono,.false.)
      CALL monchk (gamma,'gamma',two*smin,two*smax,beta,pmono,.true.)
*
      END


      SUBROUTINE monchk (f, name, xmin, xmax, beta, pmono, rising)
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
      LOGICAL             rising
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
          IF (rising .AND. (fthis .LT. flast)) THEN
              WRITE (stdout,20000) name, xlast, flast, name, x, fthis
          else IF (.NOT.rising .AND. (fthis .GT. flast)) THEN
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
