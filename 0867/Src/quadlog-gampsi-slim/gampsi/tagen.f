      PROGRAM tagen
************************************************************************
*     (Single-precision test data generator)
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
      INTRINSIC           alog,        cos,         exp,         sin
      INTRINSIC           tan
*
*     Built-in functions
*
      REAL                alog,        cos,         exp,         sin
      REAL                tan
*
*     External functions
*
      EXTERNAL            algam,       gamma,       psi,         psiln
*
      REAL                algam,       gamma,       psi,         psiln
*
*     Parameter variables
*
      REAL                two
      PARAMETER           (two = 2.0)
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
      REAL                den,         f,           fmax,        fmin
      REAL                num,         x,           xmax,        xmin
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
          IF (thefun .EQ. 'Gamma') THEN
              f = gamma(x)
          ELSE IF (thefun .EQ. 'cos') THEN
              f = cos(x)
          ELSE IF (thefun .EQ. 'exp') THEN
              f = exp(x)
          ELSE IF (thefun .EQ. 'lnGamma') THEN
              f = algam(x)
          ELSE IF (thefun .EQ. 'log') THEN
              f = alog(x)
          ELSE IF (thefun .EQ. 'psi') THEN
              f = psi(x)
          ELSE IF (thefun .EQ. 'psiln') THEN
              f = psiln(x)
          ELSE IF (thefun .EQ. 'sin') THEN
              f = sin(x)
          ELSE IF (thefun .EQ. 'tan') THEN
              f = tan(x)
          END IF
          WRITE (stdout,10000) 2, (num / den), p, 0, f, 0, 0
          GO TO 100
  200     CONTINUE
      ELSE
          IF (thefun .EQ. 'Gamma ') THEN
              CALL agen(stdout, iflag, nrand, fmin, fmax, xmin, xmax,
     X            pmin, pmax, gamma)
          ELSE IF (thefun .EQ. 'cos') THEN
              CALL agen(stdout, iflag, nrand, fmin, fmax, xmin, xmax,
     X            pmin, pmax, cos)
          ELSE IF (thefun .EQ. 'exp') THEN
              CALL agen(stdout, iflag, nrand, fmin, fmax, xmin, xmax,
     X            pmin, pmax, exp)
          ELSE IF (thefun .EQ. 'lnGamma') THEN
              CALL agen(stdout, iflag, nrand, fmin, fmax, xmin, xmax,
     X            pmin, pmax, algam)
          ELSE IF (thefun .EQ. 'log') THEN
              CALL agen(stdout, iflag, nrand, fmin, fmax, xmin, xmax,
     X            pmin, pmax, alog)
          ELSE IF (thefun .EQ. 'psi') THEN
              CALL agen(stdout, iflag, nrand, fmin, fmax, xmin, xmax,
     X            pmin, pmax, psi)
          ELSE IF (thefun .EQ. 'psiln') THEN
              CALL agen(stdout, iflag, nrand, fmin, fmax, xmin, xmax,
     X            pmin, pmax, psiln)
          ELSE IF (thefun .EQ. 'sin') THEN
              CALL agen(stdout, iflag, nrand, fmin, fmax, xmin, xmax,
     X            pmin, pmax, sin)
          ELSE IF (thefun .EQ. 'tan') THEN
              CALL agen(stdout, iflag, nrand, fmin, fmax, xmin, xmax,
     X            pmin, pmax, tan)
          END IF
      END IF
*
10000 FORMAT (I1, 1X, F22.1, 1X, I7, 1X, I1, 1X, 1P, E16.8, 1X, I1,
     X    1X, I1)
*
      END
