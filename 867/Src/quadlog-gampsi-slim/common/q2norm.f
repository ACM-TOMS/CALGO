      REAL*16 FUNCTION  q2norm(a,b)
************************************************************************
*     (Double-precision Euclidean 2-norm)
*     Find dsqrt(a**2 + b**2) without overflow or destructive underflow,
*     and handle Infinity and NaN arguments correctly.
*     (03-Nov-2003)
************************************************************************
*
      REAL                qabs,        qsqrt
*
*     Argument variables
*
      REAL*16             a,           b
*
*     Local variables
*
      REAL*16             aabs,        babs,        r
*
      aabs = qabs(a)
      babs = qabs(b)
      IF (aabs .GT. babs) THEN
          r = babs/aabs
          q2norm = aabs * qsqrt(1.0q+00 + r * r)
      ELSE
          r = aabs/babs
          q2norm = babs * qsqrt(1.0q+00 + r * r)
      END IF
*
      END
