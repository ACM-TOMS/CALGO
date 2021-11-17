      INTEGER FUNCTION  idceil(x)
************************************************************************
*     (Ceiling of double-precision x)
*     Return the least integer >= x.
*     (03-Jun-2000)
************************************************************************
*
*     Argument variables
*
      DOUBLE PRECISION    x
*
*     Local variables
*
      DOUBLE PRECISION    xint
*
      idceil = x
      xint = idceil
      IF (xint .LT. x) idceil = idceil + 1
*
      END
