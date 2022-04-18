      PROGRAM tlgam
************************************************************************
*     Test the log-gamma and gamma functions in each precision for
*     modest arguments, where plots show high errors in the log-gamma
*     functions.
*     [04-Aug-2000]
************************************************************************
      INTRINSIC exp, dexp, qexp
      EXTERNAL algam, gamma, dlgam, dgamma, qlgam, qgamma,aeps,deps,qeps
      REAL algam, gamma, exp, aeps
      DOUBLE PRECISION dlgam, dgamma, dexp, deps
      REAL*16 qlgam, qgamma, qexp, qmach, qeps
      REAL slg, sg, selg, se, sx, smach
      DOUBLE PRECISION dlg, dg, delg, de, dx, dmach
      REAL*16 qlg, qg, qelg, qe, qx
      INTEGER n
*
      INCLUDE 'stdio.inc'
*
      smach = aeps(1.0e+00)
      dmach = deps(1.0d+00)
      qmach = qeps(1.0q+00)
   10 READ (stdin,*,end=20) qx

      qlg = qlgam(qx)
      qg  = qgamma(qx)
      qelg = qexp(qlg)
      qe = (qelg - qg)/qg
      n = nint(qe/qmach)
      WRITE (stdout, '(''QP'',f10.4, 3f40.25, 1p, e9.2, i10)')
     X     qx, qlg, qg, qelg, qe, n

      dx = qx
      dlg = dlgam(dx)
      dg  = dgamma(dx)
      delg = dexp(dlg)
      de = (delg - dg)/dg
      n = nint(de/dmach)
      WRITE (stdout, '(''DP'',f10.4, 3(f30.15,10x), 1p, e9.2, i10)')
     X    dx, dlg, dg, delg, de, n

      sx = qx
      slg = algam(sx)
      sg  = gamma(sx)
      selg = exp(slg)
      se = (selg - sg)/sg
      n = nint(se/smach)
      WRITE (stdout, '(''SP'',f10.4, 3(f21.6,19x), 1p, e9.2, i10)')
     X    sx, slg, sg, selg, se, n

      WRITE (stdout,'(72(''-''))')

      GO TO 10
   20 END
