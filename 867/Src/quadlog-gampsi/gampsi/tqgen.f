      PROGRAM tqgen
************************************************************************
*     (Quadruple-precision test data generator)
*     This program servers as a data-driven generator of test data for
*     function testing.
*
*     The input file on stdin has the following format:
*
*         functionname
*         nrand
*         iflag
*         pmin pmax
*         fmin fmax
*         xmin xmax
*         pxmin pxmax
*
*     The supported functionname values are: Gamma cos exp lnGamma log
*     psi psiln sin tan.
*
*     If iflag is nonzero, this will result in the generation of nrand
*     values of functionname for x = (randint) * 2**p, for p from pmin
*     to pmax, ensuring that x always lies in xmin..xmax.  If iflag is
*     1, then only +x is generated.  If iflag is -1, then only -x is
*     generated.  If iflag is 2, then both +x and -x are generated.
*
*     The input values of xmin and xmax are multiplied by 2**pxmin and
*     2**pxmax respectively before use.  This feature allows precise
*     specification of fractional xmin and xmax.
*
*     If iflag is zero, then the next two lines are read, but their
*     data is not used, and then additional lines are read with
*
*         num den p
*
*     Each such line defines x = (num / den) * 2**p, and results in
*     one function evaluation for that x.
*
*     To add support for a new function, add another branch in the IF
*     statement in funhdr(), and another branch in the main IF
*     statement.
*
*     Because of a bug in the Cray Fortran 90 compiler, which is
*     propagated (in slightly different ways) into the derivative SGI
*     f90 compiler and Sun f90 and f95 compilers, with those compilers,
*     it is not possible to compile code in which certain
*     quadruple-precision intrinsic functions are passed as arguments to
*     other routines.  We therefore have to provide private wrappers for
*     the intrinsics that we need to use this way; these wrappers have
*     the same name as the intrinsic, but with a q prefix.
*
*     [21-Jun-2000]
************************************************************************
*
*     External functions
*
      EXTERNAL            qgamma,      qlgam,       qpsi,        qpsiln
      EXTERNAL            qqcos,       qqexp,       qqlog,       qqsin
      EXTERNAL            qqtan
*
      REAL*16             qgamma,      qlgam,       qpsi,        qpsiln
      REAL*16             qqcos,       qqexp,       qqlog,       qqsin
      REAL*16             qqtan
*
*     Parameter variables
*
      REAL*16             two
      PARAMETER           (two = 2.0q+00)
*
      INCLUDE 'stdio.inc'
*
*     Local variables
*
      CHARACTER*7         thefun
*
      INTEGER             iflag,       nrand,       p,           pmax
      INTEGER             pmin,        pxmax,       pxmin
*
      REAL*16             den,         f,           fmax,        fmin
      REAL*16             num,         x,           xmax,        xmin
*
      READ (stdin, '(A)') thefun
      CALL funhdr(thefun)
      READ (stdin, *) nrand
      READ (stdin, *) iflag
      READ (stdin, *) pmin, pmax
      READ (stdin, *) fmin, fmax
      READ (stdin, *) xmin, xmax
      READ (stdin, *) pxmin, pxmax
      xmin = xmin * two**pxmin
      xmax = xmax * two**pxmax
      IF (iflag .eq. 0) THEN
  100     READ (stdin, *, end = 200) num, den, p
          x = (num / den) * two**p
          IF (thefun .EQ. 'Gamma ') THEN
              f = qgamma(x)
          ELSE IF (thefun .EQ. 'cos') THEN
              f = qqcos(x)
          ELSE IF (thefun .EQ. 'exp') THEN
              f = qqexp(x)
          ELSE IF (thefun .EQ. 'lnGamma') THEN
              f = qlgam(x)
          ELSE IF (thefun .EQ. 'log') THEN
              f = qqlog(x)
          ELSE IF (thefun .EQ. 'psi') THEN
              f = qpsi(x)
          ELSE IF (thefun .EQ. 'psiln') THEN
              f = qpsiln(x)
          ELSE IF (thefun .EQ. 'sin') THEN
              f = qqsin(x)
          ELSE IF (thefun .EQ. 'tan') THEN
              f = qqtan(x)
          END IF
          WRITE (stdout,10000) 2, (num / den), p, 0, f, 0, 0
          GO TO 100
  200     CONTINUE
      ELSE
          IF (thefun .EQ. 'Gamma ') THEN
              CALL qgen(stdout, iflag, nrand, fmin, fmax, xmin, xmax,
     X            pmin, pmax, qgamma)
          ELSE IF (thefun .EQ. 'cos') THEN
              CALL qgen(stdout, iflag, nrand, fmin, fmax, xmin, xmax,
     X            pmin, pmax, qqcos)
          ELSE IF (thefun .EQ. 'exp') THEN
              CALL qgen(stdout, iflag, nrand, fmin, fmax, xmin, xmax,
     X            pmin, pmax, qqexp)
          ELSE IF (thefun .EQ. 'lnGamma') THEN
              CALL qgen(stdout, iflag, nrand, fmin, fmax, xmin, xmax,
     X            pmin, pmax, qlgam)
          ELSE IF (thefun .EQ. 'log') THEN
              CALL qgen(stdout, iflag, nrand, fmin, fmax, xmin, xmax,
     X            pmin, pmax, qqlog)
          ELSE IF (thefun .EQ. 'psi') THEN
              CALL qgen(stdout, iflag, nrand, fmin, fmax, xmin, xmax,
     X            pmin, pmax, qpsi)
          ELSE IF (thefun .EQ. 'psiln') THEN
              CALL qgen(stdout, iflag, nrand, fmin, fmax, xmin, xmax,
     X            pmin, pmax, qpsiln)
          ELSE IF (thefun .EQ. 'sin') THEN
              CALL qgen(stdout, iflag, nrand, fmin, fmax, xmin, xmax,
     X            pmin, pmax, qqsin)
          ELSE IF (thefun .EQ. 'tan') THEN
              CALL qgen(stdout, iflag, nrand, fmin, fmax, xmin, xmax,
     X            pmin, pmax, qqtan)
          END IF
      END IF
*
10000 FORMAT (I1, 1X, F22.1, 1X, I7, 1X, I1, 1X, 1P, E45.35E4, 1X, I1,
     X    1X, I1)
*
      END

      REAL*16 FUNCTION qqcos(x)
*
*     Intrinsic functions
*
      INTRINSIC           qcos
*
*     Built-in functions
*
      REAL*16             qcos
*
*     Argument variables
*
      REAL*16             x
*
      qqcos = qcos(x)
*
      END

      REAL*16 FUNCTION qqexp(x)
*
*     Intrinsic functions
*
      INTRINSIC           qexp
*
*     Built-in functions
*
      REAL*16             qexp
*
*     Argument variables
*
      REAL*16             x
*
      qqexp = qexp(x)
*
      END

      REAL*16 FUNCTION qqlog(x)
*
*     Intrinsic functions
*
      INTRINSIC           qlog
*
*     Built-in functions
*
      REAL*16             qlog
*
*     Argument variables
*
      REAL*16             x
*
      qqlog = qlog(x)
*
      END

      REAL*16 FUNCTION qqsin(x)
*
*     Intrinsic functions
*
      INTRINSIC           qsin
*
*     Built-in functions
*
      REAL*16             qsin
*
*     Argument variables
*
      REAL*16             x
*
      qqsin = qsin(x)
*
      END

      REAL*16 FUNCTION qqtan(x)
*
*     Intrinsic functions
*
      INTRINSIC           qtan
*
*     Built-in functions
*
      REAL*16             qtan
*
*     Argument variables
*
      REAL*16             x
*
      qqtan = qtan(x)
*
      END
