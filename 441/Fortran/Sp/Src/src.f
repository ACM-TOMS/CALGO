      REAL FUNCTION DIPOLE ( A, B, SEED )
      EXTERNAL R11
      REAL X, Y, A, B, R11
      INTEGER SEED
10    X = R11 ( SEED )
      Y = R11 ( SEED )
      IF ( 1.0 - X * X - Y * Y ) 10, 10, 20
20    DIPOLE = ( Y + B ) / ( X + A )
      RETURN
      END
