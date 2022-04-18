      PROGRAM tdgen
************************************************************************
*     (Double-precision test data generator)
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
*     [18-Jun-2000]
************************************************************************
*
*     Intrinsic functions
*
      INTRINSIC           dcos,        dexp,        dlog,        dsin
      INTRINSIC           dtan
*
*     Built-in functions
*
      DOUBLE PRECISION    dcos,        dexp,        dlog,        dsin
      DOUBLE PRECISION    dtan
*
*     External functions
*
      EXTERNAL            dgamma,      dlgam,       dpsi,        dpsiln
*
      DOUBLE PRECISION    dgamma,      dlgam,       dpsi,        dpsiln
*
*     Parameter variables
*
      DOUBLE PRECISION    two
      PARAMETER           (two = 2.0d+00)
*
      INCLUDE 'stdio.inc'
*
*     Local variables
*
      CHARACTER*7         thefun
*
      DOUBLE PRECISION    den,         f,           fmax,        fmin
      DOUBLE PRECISION    num,         x,           xmax,        xmin
*
      INTEGER             iflag,       nrand,       p,           pmax
      INTEGER             pmin,        pxmax,       pxmin
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
              f = dgamma(x)
          ELSE IF (thefun .EQ. 'cos') THEN
              f = dcos(x)
          ELSE IF (thefun .EQ. 'exp') THEN
              f = dexp(x)
          ELSE IF (thefun .EQ. 'lnGamma') THEN
              f = dlgam(x)
          ELSE IF (thefun .EQ. 'log') THEN
              f = dlog(x)
          ELSE IF (thefun .EQ. 'psi') THEN
              f = dpsi(x)
          ELSE IF (thefun .EQ. 'psiln') THEN
              f = dpsiln(x)
          ELSE IF (thefun .EQ. 'sin') THEN
              f = dsin(x)
          ELSE IF (thefun .EQ. 'tan') THEN
              f = dtan(x)
          END IF
          WRITE (stdout,10000) 2, (num / den), p, 0, f, 0, 0
          GO TO 100
  200     CONTINUE
      ELSE
          IF (thefun .EQ. 'Gamma ') THEN
              CALL dgen(stdout, iflag, nrand, fmin, fmax, xmin, xmax,
     X            pmin, pmax, dgamma)
          ELSE IF (thefun .EQ. 'cos') THEN
              CALL dgen(stdout, iflag, nrand, fmin, fmax, xmin, xmax,
     X            pmin, pmax, dcos)
          ELSE IF (thefun .EQ. 'exp') THEN
              CALL dgen(stdout, iflag, nrand, fmin, fmax, xmin, xmax,
     X            pmin, pmax, dexp)
          ELSE IF (thefun .EQ. 'lnGamma') THEN
              CALL dgen(stdout, iflag, nrand, fmin, fmax, xmin, xmax,
     X            pmin, pmax, dlgam)
          ELSE IF (thefun .EQ. 'log') THEN
              CALL dgen(stdout, iflag, nrand, fmin, fmax, xmin, xmax,
     X            pmin, pmax, dlog)
          ELSE IF (thefun .EQ. 'psi') THEN
              CALL dgen(stdout, iflag, nrand, fmin, fmax, xmin, xmax,
     X            pmin, pmax, dpsi)
          ELSE IF (thefun .EQ. 'psiln') THEN
              CALL dgen(stdout, iflag, nrand, fmin, fmax, xmin, xmax,
     X            pmin, pmax, dpsiln)
          ELSE IF (thefun .EQ. 'sin') THEN
              CALL dgen(stdout, iflag, nrand, fmin, fmax, xmin, xmax,
     X            pmin, pmax, dsin)
          ELSE IF (thefun .EQ. 'tan') THEN
              CALL dgen(stdout, iflag, nrand, fmin, fmax, xmin, xmax,
     X            pmin, pmax, dtan)
          END IF
      END IF
*
10000 FORMAT (I1, 1X, F22.1, 1X, I7, 1X, I1, 1X, 1P, E26.18E3, 1X, I1,
     X    1X, I1)
*
      END
