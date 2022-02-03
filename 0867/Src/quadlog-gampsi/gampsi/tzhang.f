      PROGRAM tzhang
************************************************************************
*     (Zhang/Jin test)
*     Generate a table of high-precision values of gamma(x) for arguments
*     tabulated in the book
*
*         Zhang, S. and Jin, J. 1996.
*         Computation of special functions.
*         Wiley, New York, NY, USA.
*
*     [16-Jun-2000]
************************************************************************
*
*     External functions
*
      EXTERNAL            qgamma
*
      REAL*16             qgamma
*
*     Local variables
*
      INTEGER             k,           lineno
*
      REAL*16             x
*
      lineno = 0
*
      DO 100 k = 0, 100
          x = k + 1
          CALL show(k, qgamma(x), lineno)
  100 CONTINUE
*
      lineno = 0
      DO 200 k = 100, 240, 10
          x = k + 1
          CALL show(k, qgamma(x), lineno)
  200 CONTINUE
*
      DO 300 k = 250, 450, 50
          x = k + 1
          CALL show(k, qgamma(x), lineno)
  300 CONTINUE
*
      DO 400 k = 500, 3000, 100
          x = k + 1
          CALL show(k, qgamma(x), lineno)
  400 CONTINUE
*
      END

      SUBROUTINE show(n, g, lineno)
*
*     Argument variables
*
      INTEGER             lineno,      n
*
      REAL*16             g
*
*     Local variables
*
      CHARACTER*45        s
*
      INTEGER             k
*
      INCLUDE 'stdio.inc'
*
      WRITE (s,'(1p, e45.35e4)') g
      WRITE (stdout,'(1x,i4,4x,a,10(a,1x))') n,
     X     s(1:4), (s(k:k+4), k = 5, 40, 5)
      lineno = lineno + 1
      IF (mod(lineno, 5) .EQ. 0) write (stdout,'()')
      END
