      REAL FUNCTION anan()
************************************************************************
*     (Single-precision NaN)
*     Return a run-time trappable NaN.
*     [16-Jun-2000]
************************************************************************
*
*     External functions
*
      EXTERNAL            astorf
*
      REAL                astorf
*
*     Parameter variables
*
      REAL                zero
      PARAMETER           (zero = 0.0)
*
*     Local variables
*
      REAL                nan
*
      nan = zero
      anan = nan / astorf(nan)
*
      END
