      REAL*16 FUNCTION  qerbit(relerr,ulp)
************************************************************************
*     (Double-precision error in bits)
*     Return the error in bits corresponding to relerr, given machine
*     precision ulp.
*     (02-May-2000)
************************************************************************
*
*     External functions
*
      EXTERNAL            qlog2
*
      REAL*16             qabs,        qlog2
*
*     Parameter variables
*
      REAL*16             ZERO
      PARAMETER           (ZERO = 0.0q+00)
*
*     Argument variables
*
      REAL*16    relerr,      ulp
*
      IF (relerr .EQ. ZERO) THEN
          qerbit = ZERO
      ELSE
          qerbit = qabs(qlog2(qabs(relerr/ulp)))
      END IF
*
      END
