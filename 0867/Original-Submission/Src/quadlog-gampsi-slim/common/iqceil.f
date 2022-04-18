      INTEGER FUNCTION  iqceil(x)
************************************************************************
*     (Ceiling of quadruple-precision x)
*     Return the least integer >= x.
*     (03-Jun-2000)
************************************************************************
*
*     Argument variables
*
      REAL*16             x
*
*     Local variables
*
      REAL*16             xint
*
      iqceil = x
      xint = iqceil
      IF (xint .LT. x) iqceil = iqceil + 1
*
      END
