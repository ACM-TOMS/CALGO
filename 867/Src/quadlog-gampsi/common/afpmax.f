      REAL FUNCTION  afpmax()
************************************************************************
*     (Single-precision floating-point maximum)
*     Return the largest finite representable single-precision
*     floating-point number.
*     (03-Jun-2000)
************************************************************************
*
*     Parameter variables
*
      REAL                half
      PARAMETER           (half = 0.5)
*
      REAL                one
      PARAMETER           (one = 1.0)
*
*     Local variables
*
      LOGICAL             first
*
      REAL                fpmax,       x,           x2,          x3
      REAL                y,           z,           zplusy
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
*          execution!  The calls to astore() are critical: they foil
*          machines (e.g., Honeywell, Intel x86, Motorola 68K) in
*          which the floating-point registers have range and/or
*          precision beyond that of floating-point values in memory.
*
           x = one
   10      x2 = x + x
           CALL astore(x2)
           x3 = x2 + x
           CALL astore(x3)
           IF (x3 .NE. x2) THEN
                CALL astore(x)
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
           CALL astore(zplusy)
           IF ((zplusy .NE. z) .AND. (zplusy .LT. x)) THEN
                CALL astore(z)
                CALL astore(y)
                z = zplusy
                y = half * y
                GO TO 20
           END IF
*
           fpmax = z + z
*
      END IF
*
      afpmax = fpmax
*
      END
