      SUBROUTINE qstore(x)
************************************************************************
*     (Quadruple-precision store x)
*     Store x, to force a variable into memory.
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
      END
