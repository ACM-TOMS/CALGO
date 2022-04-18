       DOUBLE PRECISION FUNCTION dpsi(x)
*      (Double-precision interface to quadruple-precision psi(x))
*      Simple interface to extended precision psi() function
*      (31-Jul-2000)
       DOUBLE PRECISION x
       EXTERNAL qpsi
       REAL*16 qpsi, qx
       qx = x
       dpsi = qpsi(qx)
       END
