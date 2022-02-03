      LOGICAL FUNCTION isainf(x)
************************************************************************
*     (Is single-precision x infinite?)
*     Return .TRUE. if x is infinite, and .FALSE. otherwise.
*     [12-Jun-2000]
************************************************************************
*
*     Parameter variables
*
      REAL                zero
      PARAMETER           (zero = 0.0)
*
*     Argument variables
*
      REAL                x
*
      isainf = (x .NE. zero) .AND. ((x + x) .EQ. x)
*
      END
