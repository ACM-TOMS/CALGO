      PROGRAM chkran
************************************************************************
*     (Check ran())
*     Check the return values of ran() by writing a file containing mx
*     sequential values in hexadecimal, for comparison between different
*     compilers and architectures.
*     [31-Jul-2000]
************************************************************************
*
*     External functions
*
      EXTERNAL            ran
*
      REAL                ran
*
*     Parameter variables
*
      INTEGER             mx
      PARAMETER           (mx = 100000)
*
*     Local variables
*
      INTEGER             ixx(1),      k
*
      REAL                xx
*
      EQUIVALENCE         (ixx(1), xx)
*
      INCLUDE 'stdio.inc'
*
      DO 10 k = 1, mx
           xx = ran()
           WRITE (stdout, '(i7,2x,z8.8)') k, ixx
   10 CONTINUE
      END
