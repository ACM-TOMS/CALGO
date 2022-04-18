      PROGRAM tmdgam
************************************************************************
*     (Monotonicity check)
*     Do approximate monotonicity checks for dgamma(x) at boundaries
*     where the approximation changes.
*     [12-Jul-2000]
************************************************************************
*
*     External functions
*
      EXTERNAL            dgamma
*
      DOUBLE PRECISION    dgamma
*
*     Parameter variables
*
      DOUBLE PRECISION    beta
      PARAMETER           (beta = 2.0d+00)
*
      DOUBLE PRECISION    one
      PARAMETER           (one = 1.0d0)
*
      DOUBLE PRECISION    two
      PARAMETER           (two = 2.0d0)
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
      INCLUDE 'deps.inc'
*
      CALL monchk (dgamma,'dgamma',eps*smin,eps*smax,beta,pmono,.false.)
      CALL monchk (dgamma,'dgamma',one*smin,one*smax,beta,pmono,.false.)
      CALL monchk (dgamma,'dgamma',two*smin,two*smax,beta,pmono,.true.)
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
      DOUBLE PRECISION    f
*
      INCLUDE 'stdio.inc'
*
*     Argument variables
*
      CHARACTER*(*)       name
      DOUBLE PRECISION    beta
      INTEGER             pmono
      LOGICAL             rising
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
     X      1p, 4x, a,'(', e27.17e4, ') = ', e27.17e4/
     X      1p, 4x, a,'(', e27.17e4, ') = ', e27.17e4)
30000 FORMAT (a, '(x): completed ', I12, ' monotonicity checks on (',
     X      1p, e27.17e4, ',', e27.17e4, ')')
      END
