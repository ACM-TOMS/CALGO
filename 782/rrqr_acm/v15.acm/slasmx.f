      REAL             FUNCTION SLASMX( I )
      INTEGER I
*
      REAL             OTHIRD
      PARAMETER ( OTHIRD = 1.0E+0/3.0E+0 )
      INTRINSIC REAL
      SLASMX = REAL( I )**OTHIRD
      RETURN
      END
