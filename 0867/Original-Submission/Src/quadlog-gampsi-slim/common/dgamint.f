       DOUBLE PRECISION FUNCTION dgamma(x)
************************************************************************
*      (Double-precision interface to quadruple-precision Gamma())
*      Simple interface to extended precision Gamma() function
*      (31-Jul-2000)
************************************************************************
       DOUBLE PRECISION x
       EXTERNAL qgamma
       REAL*16 qgamma, qx
       qx = x
       dgamma = qgamma(qx)
       END
