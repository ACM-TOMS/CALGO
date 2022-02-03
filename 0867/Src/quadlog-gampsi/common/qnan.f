      REAL*16 FUNCTION qnan()
************************************************************************
*     (Quadruple-precision NaN)
*     Return a run-time trappable NaN.
*     [16-Jun-2000]
************************************************************************
*
*     External functions
*
      EXTERNAL            qstorf
*
      REAL*16             qstorf
*
*     Parameter variables
*
      REAL*16             zero
      PARAMETER           (zero = 0.0q+00)
*
*     Local variables
*
      REAL*16             nan
*
      nan = zero
      qnan = nan / qstorf(nan)
*
      END
