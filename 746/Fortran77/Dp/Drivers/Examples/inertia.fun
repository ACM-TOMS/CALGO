C
C     GEAR TRAIN OF MINIMUM INERTIA (TP328)
C
*     SET OF INDICES
      IND2 = 1..2
C
*     VARIABLE
      X(I), I IN IND2
C
*     FUNCTION F
      A=(0.1D1+X(2)**2)/X(1)**2
      B=((X(1)*X(2))**2+0.1D3)/(X(1)*X(2))**4
      F=(0.12D2+X(1)**2+A+B)/0.1D2
C
*     END
