      DOUBLE PRECISION FUNCTION  d2norm(a,b)
************************************************************************
*     (Double-precision Euclidean 2-norm)
*     Find dsqrt(a**2 + b**2) without overflow or destructive underflow,
*     and handle Infinity and NaN arguments correctly.
*     (03-Nov-2003)
************************************************************************
*
      DOUBLE PRECISION   dabs,        dsqrt
*
*     Argument variables
*
      DOUBLE PRECISION   a,           b
*
*     Local variables
*
      DOUBLE PRECISION   aabs,        babs,         r
*
      aabs = dabs(a)
      babs = dabs(b)
      IF (aabs .GT. babs) THEN
          r = babs/aabs
          d2norm = aabs * dsqrt(1.0d+00 + r * r)
      ELSE
          r = aabs/babs
          d2norm = babs * dsqrt(1.0d+00 + r * r)
      END IF
*
      END
