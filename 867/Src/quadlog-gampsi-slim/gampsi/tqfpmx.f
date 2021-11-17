      PROGRAM tqfpmx
************************************************************************
*     (Test qfpmax())
*     Print the result returned by qfpmax, in decimal and hexadecimal.
*     [03-Jun-2000]
************************************************************************
*
*     External functions
*
      EXTERNAL            qfpmax
*
      REAL*16             qfpmax
*
      INCLUDE 'stdio.inc'
*
*     Local variables
*
      INTEGER             ix(4)
*
      REAL*16             x
*
*     Some poorly-designed Fortran implementations do not permit
*     writing floating-point data with a hexadecimal format, so we
*     must overlay an integer array on x, and handle big-endian and
*     little-endian cases, sigh...
*
      EQUIVALENCE         (x, ix(1))
*
      x = qfpmax()
      IF ((ix(1) .LT. 0) .AND. (ix(4) .GT. 0)) THEN
          WRITE (stdout,'(''Computed: '', 1p, e45.35e4, 2x, 4z9.8)')
     X        x, ix(4), ix(3), ix(2), ix(1)
      ELSE
          WRITE (stdout,'(''Computed: '', 1p, e45.35e4, 2x, 4z9.8)')
     X        x, ix(1), ix(2), ix(3), ix(4)
      END IF
*
      WRITE (stdout,'(''Expected: '', A)')
     X    '  1.18973149535723176508575932662800702E+4932 ' //
     X    '  7FFEFFFF FFFFFFFF FFFFFFFF FFFFFFFF'
*
      END
