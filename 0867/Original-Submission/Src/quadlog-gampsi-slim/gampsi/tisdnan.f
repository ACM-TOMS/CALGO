      PROGRAM tisdna
      DOUBLE PRECISION nan, zero
      LOGICAL isdnan, result
      INTEGER inan(2)
      EQUIVALENCE (nan, inan(1))
*
      INCLUDE 'stdio.inc'
*
      zero = 0.0d+00
      nan = zero
      nan = nan / zero
*
      WRITE (stdout,*) 'Test of isdnan()'
*
      result =  isdnan(nan)
      WRITE (stdout,*) 'isdnan(', nan, ') = ', result
      WRITE (stdout,'(1x, 2z9.8)') inan
*
      END

      LOGICAL FUNCTION isdnan(x)
      INTEGER             iand,        ishft
      DOUBLE PRECISION    x
      DOUBLE PRECISION    xx
      INTEGER             hi,          ix(2),       lo
      EQUIVALENCE         (ix(1), xx)
      xx = 1.0d+00
      IF (ix(1) .EQ. 0) THEN
          hi = 2
          lo = 1
      ELSE
          hi = 1
          lo = 2
      END IF
      xx = x
      isdnan = ((iand(2047,ishft(ix(hi), -20)) .EQ. 2047) .AND.
     X         ((iand(ix(hi),1048575) .NE. 0) .OR. (ix(lo) .NE. 0)))
      END
