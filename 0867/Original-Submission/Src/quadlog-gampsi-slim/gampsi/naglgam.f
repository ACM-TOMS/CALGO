      DOUBLE PRECISION FUNCTION dlgam(x)
      DOUBLE PRECISION x
      DOUBLE PRECISION s14abf
      EXTERNAL s14abf
      INTEGER ifail
      ifail = 1
      dlgam = s14abf(x,ifail)
      END
