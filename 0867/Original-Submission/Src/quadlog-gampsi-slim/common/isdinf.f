      LOGICAL FUNCTION isdinf(x)
************************************************************************
*     (Is double-precision x infinite?)
*     Return .TRUE. if x is infinite, and .FALSE. otherwise.
*     [12-Jun-2000]
************************************************************************
*
*     Parameter variables
*
      DOUBLE PRECISION    zero
      PARAMETER           (zero = 0.0d+00)
*
*     Argument variables
*
      DOUBLE PRECISION    x
*
      isdinf = (x .NE. zero) .AND. ((x + x) .EQ. x)
*
      END
