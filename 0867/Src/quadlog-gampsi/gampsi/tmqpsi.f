      PROGRAM tmqpsi
************************************************************************
*     (Monotonicity check)
*     Do approximate monotonicity checks for qpsi(x) at boundaries
*     where the approximation changes.
*     [12-Jul-2000]
************************************************************************
*
*     External functions
*
      EXTERNAL            qinf,        qpsi
*
      REAL*16             qinf,        qpsi
*
*     Parameter variables
*
      REAL*16             beta
      PARAMETER           (beta = 2.0d+00)
*
      REAL*16             one
      PARAMETER           (one = 1.0d0)
*
      REAL*16             two
      PARAMETER           (two = 2.0d0)
*
*     This value should be p-1: 111 for full IEEE 754 quadruple
*     precision, and 105 for IEEE 754 paired double.
*
      INTEGER             pmono
      PARAMETER           (pmono = 111)
*
      REAL*16             delta
      PARAMETER           (delta = 1024.0q+00 * beta**(-pmono))
*
      REAL*16             smin
      PARAMETER           (smin = one - delta)
*
      REAL*16             smax
      PARAMETER           (smax = one + delta)
*
      INCLUDE 'qxlarg.inc'
*
      INCLUDE 'qxmin.inc'
      REAL*16             xmin1
      PARAMETER           (xmin1 = xmin)
*
      INCLUDE 'qxsmal.inc'
*
      REAL*16             cutoff
      PARAMETER           (cutoff = 20.0q+00)
*
*     CALL monchk (qpsi, 'qpsi', -2.0q+00, -1.0q+00)
*     CALL monchk (qpsi, 'qpsi', 1.4617q+00, qinf())
      CALL monchk (qpsi, 'qpsi', -xmin1*smax, -xmin1*smin, beta,pmono)
      CALL monchk (qpsi, 'qpsi', xmin1*smin, xmin1*smax, beta,pmono)
      CALL monchk (qpsi, 'qpsi', xsmall*smin, xsmall*smax, beta,pmono)
      CALL monchk (qpsi, 'qpsi', one*smin, one*smax, beta,pmono)
      CALL monchk (qpsi, 'qpsi', two*smin, two*smax, beta,pmono)
      CALL monchk (qpsi, 'qpsi', cutoff*smin, cutoff*smax, beta,pmono)
      CALL monchk (qpsi, 'qpsi', xlarge*smin, xlarge*smax, beta,pmono)
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
      REAL*16             f
*
      INCLUDE 'stdio.inc'
*
*     Argument variables
*
      CHARACTER*(*)       name
      REAL*16             beta
      INTEGER             pmono
*
      REAL*16             xmax,        xmin
*
*     Parameter variables
*
      REAL*16             one
      PARAMETER           (one = 1.0q+00)
*
      REAL*16             zero
      PARAMETER           (zero = 0.0q+00)
*
*     Local variables
*
      REAL*16             flast,       fthis,       x,           xlast
      REAL*16             xmult
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
      WRITE (stdout,10000) iqint(beta), -pmono
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
10000 FORMAT ('Stepsize = ', i2, '**(',i4,')')
20000 FORMAT ('Monotonicity failure:'/
     X      1p, 4x, a,'(', e45.35e4, ') = ', e45.35e4/
     X      1p, 4x, a,'(', e45.35e4, ') = ', e45.35e4)
30000 FORMAT (a, '(x): completed ', I12, ' monotonicity checks on (',
     X      1p, e45.35e4, ',', e45.35e4, ')')
      END
