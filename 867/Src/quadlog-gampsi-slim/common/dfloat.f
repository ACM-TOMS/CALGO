      DOUBLE PRECISION FUNCTION dfloat(n)
************************************************************************
*     (Integer to double-precision)
*     Return the integer n converted to double-precision floating-point.
*
*     Although almost all Fortran 66 and 77 implementations provide
*     this function, it was regrettably omitted from both Fortran 66
*     and 77 language standards.
*     (12-Dec-1999)
************************************************************************
      INTEGER             n
*
      dfloat = n
*
      END
