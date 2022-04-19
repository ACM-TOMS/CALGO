      DOUBLE PRECISION FUNCTION  derbit(relerr,ulp)
************************************************************************
*     (Double-precision error in bits)
*     Return the error in bits corresponding to relerr, given machine
*     precision ulp.
*     (02-May-2000)
************************************************************************
*
*     External functions
*
      EXTERNAL            dlog2
*
      DOUBLE PRECISION    dabs,        dlog2
*
*     Parameter variables
*
      DOUBLE PRECISION    zero
      PARAMETER           (zero = 0.0d+00)
*
*     Argument variables
*
      DOUBLE PRECISION    relerr,      ulp
*
      IF (relerr .EQ. zero) THEN
          derbit = zero
      ELSE
          derbit = dabs(dlog2(dabs(relerr/ulp)))
      END IF
*
      END
