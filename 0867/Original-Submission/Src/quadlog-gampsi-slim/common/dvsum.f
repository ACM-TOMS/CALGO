      DOUBLE PRECISION FUNCTION dvsum(x,n)
*     (Double-precision vector sum)
*     Return an accurate estimate of the sum x(1) + x(2) + ... + x(n)
*     using the Kahan compensated summation.  This is a reasonable
*     compromise that does not require explicit support for higher
*     precision, and avoids the complexity of the more sophisticated
*     (and accurate) summation methods described in Nicholas J. Higham,
*     ``Accuracy and stability of numerical algorithms'', SIAM,
*     2002, ISBN 0-89871-521-0, and I. J. Anderson, ``A Distillation
*     Algorithm for Floating-Point Summation'', SIAM J. Sci. Comp.
*     20(5), 1797--1806, September 1999.
*     [04-Oct-2003]
*
*     Argument variables
*
      DOUBLE PRECISION    x(*)
*
      INTEGER             n
*
*     Local variables
*
      DOUBLE PRECISION    a,           b,           e,           s
*
      INTEGER             i
*
      e = 0.0d+00
      s = e
      DO 100 i = 1, n
          a = s
          b = x(i) + e
          s = a + b
          e = (a - s) + b
  100 CONTINUE
      dvsum = s
      END
