      DOUBLE PRECISION FUNCTION dpsiln(x)
      DOUBLE PRECISION x
      DOUBLE PRECISION s14acf
      EXTERNAL s14acf
      INTEGER ifail
      ifail = 1
      dpsiln = s14acf(x,ifail)
      END
