      PROGRAM chkiran
************************************************************************
*     (Check airan())
*     Check the return values of airan() by writing a file containing mx
*     sequential values in hexadecimal, for comparison between different
*     compilers and architectures.
*     [31-Jul-2000]
************************************************************************
*
*     External functions
*
      EXTERNAL            airan
*
      REAL                airan
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
           xx = airan(0.0, 8388607.0)
           WRITE (stdout, '(i7,2x,z8.8)') k, ixx
   10 CONTINUE
      END
