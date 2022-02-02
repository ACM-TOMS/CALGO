      REAL FUNCTION  avsum(x,n)
*     (Single-precision vector sum)
*     Return an accurate estimate of the sum x(1) + x(2) + ... + x(n)
*     using double-precision summation.
*     [04-Oct-2003]
*
      DOUBLE PRECISION    dble
*
      REAL                sngl
*
*     Argument variables
*
      INTEGER             n
*
      REAL                x(*)
*
*     Local variables
*
      DOUBLE PRECISION    s
*
      INTEGER             i
*
      s = 0.0d+00
      DO 100 i = 1, n
          s = s + dble(x(i))
  100 CONTINUE
      avsum = sngl(s)
*
      END
