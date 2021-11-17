      DOUBLE PRECISION FUNCTION dpsiln(x)
      DOUBLE PRECISION x, ans(1)
      CALL dpsifn(x, 0, 2, 1, ans)
      dpsiln = -ans(1)
      END
