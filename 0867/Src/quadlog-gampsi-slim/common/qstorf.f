      REAL*16 FUNCTION qstorf(x)
************************************************************************
*     (Quadruple-precision store and return x)
*     Store and return x, to force a variable into memory.
*
*     This action is necessary on some architectures to force an
*     expression to be converted to storage precision, when it might
*     otherwise be held in a machine register at higher internal
*     precision.
*     [10-Jun-2000]
************************************************************************
*
*     Argument variables
*
      REAL*16             x
*
      qstorf = x
*
      END
