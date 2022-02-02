      SUBROUTINE funhdr(thefun)
************************************************************************
*     (Function header)
*     Generate a data file header corresponding to thefun.
*     [01-Aug-2000]
************************************************************************
*
*     Parameter variables
*
      INCLUDE 'stdio.inc'
*
*     Argument variables
*
      CHARACTER*(*)       thefun
*
      IF (thefun .EQ. 'Gamma ') THEN
           WRITE (stdout,'(''### Numerical evaluation of Gamma(x)'')')
           WRITE (stdout,'(''###'')')
           WRITE (stdout,'(''### where Gamma(n+1) = factorial(n)'')')
      ELSE IF (thefun .EQ. 'cos') THEN
           WRITE (stdout,'(''### Numerical evaluation of cos(x)'')')
      ELSE IF (thefun .EQ. 'exp') THEN
           WRITE (stdout,'(''### Numerical evaluation of exp(x)'')')
      ELSE IF (thefun .EQ. 'lnGamma') THEN
           WRITE (stdout,
     X         '(''### Numerical evaluation of ln(Gamma(x))'')')
      ELSE IF (thefun .EQ. 'log') THEN
           WRITE (stdout,'(''### Numerical evaluation of log(x)'')')
      ELSE IF (thefun .EQ. 'psi') THEN
           WRITE (stdout,'(''### Numerical evaluation of psi(x)'')')
           WRITE (stdout,'(''###'')')
           WRITE (stdout,
     X      '(''### where psi(x) (also called the digamma function)'')')
           WRITE (stdout,
     X  '(''### is the derivative of ln(Gamma(x)) with respect to x'')')
      ELSE IF (thefun .EQ. 'psiln') THEN
           WRITE (stdout,'(''### Numerical evaluation of psiln(x)'')')
           WRITE (stdout,'(''###'')')
           WRITE (stdout,'(''### where psiln(x) = psi(x) - ln(x),'')')
           WRITE (stdout,
     X      '(''### and psi(x) (also called the digamma function)'')')
           WRITE (stdout,
     X  '(''### is the derivative of ln(Gamma(x)) with respect to x'')')
      ELSE IF (thefun .EQ. 'sin') THEN
           WRITE (stdout,'(''### Numerical evaluation of sin(x)'')')
      ELSE IF (thefun .EQ. 'tan') THEN
           WRITE (stdout,'(''### Numerical evaluation of tan(x)'')')
      ELSE
          STOP 'Unsupported function name'
      END IF
      WRITE (stdout,'(''###'')')
      WRITE (stdout,'(''### x = f * 2^p'')')
      WRITE (stdout,'(''###'')')
      WRITE (stdout,
     X    '(''### Line format: 2 f p neval result relerr abserr'')')
      WRITE (stdout,'(''###'')')
*
      END
