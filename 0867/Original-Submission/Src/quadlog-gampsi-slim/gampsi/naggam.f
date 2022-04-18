      DOUBLE PRECISION FUNCTION dgamma(x)
      DOUBLE PRECISION x
      DOUBLE PRECISION s14aaf
      EXTERNAL s14aaf
      INTEGER ifail
      ifail = 1
      dgamma = s14aaf(x,ifail)
      END
