      DOUBLE PRECISION FUNCTION dnan()
************************************************************************
*     (Double-precision NaN)
*     Return a run-time trappable NaN.
*     [16-Jun-2000]
************************************************************************
*
*     External functions
*
      EXTERNAL            dstorf
*
      DOUBLE PRECISION    dstorf
*
*     Parameter variables
*
      DOUBLE PRECISION    zero
      PARAMETER           (zero = 0.0d+00)
*
*     Local variables
*
      DOUBLE PRECISION    nan
*
      nan = dstorf(zero)
      dnan = nan / dstorf(nan)
*
      END
