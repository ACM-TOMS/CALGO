      INTEGER FUNCTION  iceil(x)
************************************************************************
*     (Ceiling of single-precision x)
*     Return the least integer >= x.
*     (03-Jun-2000)
************************************************************************
*
*     Argument variables
*
      REAL                x
*
*     Local variables
*
      REAL                xint
*
      iceil = x
      xint = iceil
      IF (xint .LT. x) iceil = iceil + 1
*
      END
