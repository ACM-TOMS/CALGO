      PROGRAM tafpmx
************************************************************************
*     (Test afpmax())
*     Print the result returned by afpmax, in decimal and hexadecimal.
*     [03-Jun-2000]
************************************************************************
*
*     External functions
*
      EXTERNAL            afpmax
*
      REAL             afpmax
*
      INCLUDE 'stdio.inc'
*
*     Local variables
*
      INTEGER             ix(1)
*
      REAL             x
*
*     Some poorly-designed Fortran implementations do not permit
*     writing floating-point data with a hexadecimal format, so we
*     must overlay an integer array on x, sigh...
*
      EQUIVALENCE         (x, ix(1))
*
      x = afpmax()
      WRITE (stdout,'(''Computed: '', 1p, e15.7, 2x, z9.8)')
     X    x, ix(1)
*
      WRITE (stdout,'(''Expected: '', A)')
     X    '  3.4028235E+38 ' //
     X    '  7F7FFFFF'
*
      END
