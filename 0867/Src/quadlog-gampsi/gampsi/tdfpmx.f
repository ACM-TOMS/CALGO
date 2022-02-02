      PROGRAM tdfpmx
************************************************************************
*     (Test dfpmax())
*     Print the result returned by dfpmax, in decimal and hexadecimal.
*     [03-Jun-2000]
************************************************************************
*
*     External functions
*
      EXTERNAL            dfpmax
*
      DOUBLE PRECISION    dfpmax
*
      INCLUDE 'stdio.inc'
*
*     Local variables
*
      DOUBLE PRECISION    x
*
      INTEGER             ix(2)
*
*     Some poorly-designed Fortran implementations do not permit
*     writing floating-point data with a hexadecimal format, so we
*     must overlay an integer array on x, and handle big-endian and
*     little-endian cases, sigh...
*
      EQUIVALENCE         (x, ix(1))
*
      x = dfpmax()
      IF (ix(2) .eq. -1) THEN
          WRITE (stdout,'(''Computed: '', 1p, e25.15e4, 2x, 2z9.8)')
     X        x, ix(1), ix(2)
      ELSE
          WRITE (stdout,'(''Computed: '', 1p, e25.15e4, 2x, 2z9.8)')
     X        x, ix(2), ix(1)
      END IF
*
      WRITE (stdout,'(''Expected: '', A)')
     X     '  1.797693134862316E+0308  7FEFFFFF FFFFFFFF'
*
      END
