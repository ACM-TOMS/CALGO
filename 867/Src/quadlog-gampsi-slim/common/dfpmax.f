      DOUBLE PRECISION FUNCTION  dfpmax()
************************************************************************
*     (Double-precision floating-point maximum)
*     Return the largest finite representable double-precision
*     floating-point number.
*     (03-Jun-2000)
************************************************************************
*
*     Parameter variables
*
      DOUBLE PRECISION    half
      PARAMETER           (half = 0.5d+00)
*
      DOUBLE PRECISION    one
      PARAMETER           (one = 1.0d+00)
*
*     Local variables
*
      DOUBLE PRECISION    fpmax,       x,           x2,          x3
      DOUBLE PRECISION    y,           z,           zplusy
*
      LOGICAL             first
*
      SAVE                first,       fpmax
*
      DATA first /.TRUE./
*
      IF (first) THEN
           first = .FALSE.
*
*          Compute the largest finite representable power of two by
*          successive doubling.  This code CRUCIALLY DEPENDS on IEEE
*          754 nonstop computing: an overflow must not terminate
*          execution!  The calls to dstore() are critical: they foil
*          machines (e.g., Honeywell, Intel x86, Motorola 68K) in
*          which the floating-point registers have range and/or
*          precision beyond that of floating-point values in memory.
*
           x = one
   10      x2 = x + x
           CALL dstore(x2)
           x3 = x2 + x
           CALL dstore(x3)
           IF (x3 .NE. x2) THEN
                CALL dstore(x)
                x = x + x
                GO TO 10
           END IF
*
*          Compute the large finite representable number by successive
*          appending of a one bit to the significand.  We have to work
*          with (x/2) instead of x, so that the last addition does not
*          overflow.
*
           z = half * x
           y = half * z
   20      zplusy = z + y
           CALL dstore(zplusy)
           IF ((zplusy .NE. z) .AND. (zplusy .LT. x)) THEN
                CALL dstore(z)
                CALL dstore(y)
                z = zplusy
                y = half * y
                GO TO 20
           END IF
*
           fpmax = z + z
*
      END IF
*
      dfpmax = fpmax
*
      END
