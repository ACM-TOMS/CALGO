      DOUBLE PRECISION FUNCTION       getnan()
************************************************************************
*     (Get a NaN)
*     Attempt to construct a NaN dynamically, without assuming host
*     IEEE 754 arithmetic.  The code may produce one overflow on the
*     first call, but never a zero divide error.
*     If a NaN cannot be generated, try to make an Infinity instead,
*     and if that fails,
*     (04-May-2000)
************************************************************************
*
*     Parameter variables
*
      DOUBLE PRECISION   four
      PARAMETER (four = 4.0d+00)
*
      DOUBLE PRECISION   one
      PARAMETER (one = 1.0d+00)
*
      DOUBLE PRECISION   three
      PARAMETER (three = 3.0d+00)
*
      DOUBLE PRECISION   zero
      PARAMETER (zero = 0.0d+00)
*
*     Local variables
*
      DOUBLE PRECISION   nan,         x,           y
*
      INTEGER            n
*
      SAVE               nan
*
*     This is near the largest floating-point value on the system with
*     the smallest double-precision exponent range: IBM S/360: 16^63 =
*     7.23700557733226d+75.  Its sign serves as a first-time-through
*     flag.
*
      DATA nan/ -7.237d+75/
*
*     First time through, initialize nan to proper value
*
      IF (nan .LT. zero) THEN
*
*         The next three statements should produce
*         x = 0.111111... (binary):
*
           x = four/three
           y = x - one
           x = y + y + y
*
*         Double x repeatedly until we reach an overflow condition.  n
*         serves as a loop limit, and should be large enough to handle
*         the largest current exponent: 2^8191 on a Cray Y-MP system.
*
           n = 0
   10      IF (x .LT. (x + x)) THEN
                n = n + 1
                x = x + x
                IF (n .LT. 10000) GO TO 10
           END IF
*
*         In IEEE 754 arithmetic, we now have Infinity, and the
*         difference of two such values must be a NaN.  On non-IEEE
*         754 systems, or with compilers that improperly optimize x-x
*         to zero, we may have to make do with Infinity, or a positive
*         large floating-point value.
*
           y = x - x
           IF (y .NE. zero) THEN
                nan = y
           ELSE IF (x .EQ. (x + x)) THEN
                nan = x
           ELSE
                nan = - nan
           END IF
      END IF
*
      getnan = nan
*
      END
