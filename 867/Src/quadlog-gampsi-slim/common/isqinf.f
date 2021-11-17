      LOGICAL FUNCTION isqinf(x)
************************************************************************
*     (Is quadruple-precision x infinite?)
*     Return .TRUE. if x is infinite, and .FALSE. otherwise.
*     [12-Jun-2000]
************************************************************************
*
*     Parameter variables
*
      REAL*16             zero
      PARAMETER           (zero = 0.0q+00)
*
*     Argument variables
*
      REAL*16             x
*
      isqinf = (x .NE. zero) .AND. ((x + x) .EQ. x)
*
      END
