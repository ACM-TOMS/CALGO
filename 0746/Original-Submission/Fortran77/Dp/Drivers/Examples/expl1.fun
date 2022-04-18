C-----------------------------------------------------------------------
C
C     Problem:  EXPL1 (Lindstroem no. 1)
C
C     ....  Parameter estimation with explicit model function and
C           two nonlinear equality constraints
C
C-----------------------------------------------------------------------
C
*     VARIABLE
      X1, X2, X3, X4, T
C
*     REAL CONSTANT
      TE = 4.
      YE = 0.1957
      T1 = 0.0625
      Y1 = 0.0246
C
*     FUNCTION F
      F = X1*(T**2 + X2*T)/(T**2 + X3*T + X4)
C
*     FUNCTION G1
      G1 = X1*T1*(T1 + X2)/(T1**2 + X3*T1 + X4) - Y1
C
*     FUNCTION G2
      G2 = X1*TE*(TE + X2)/(TE**2 + X3*TE + X4) - YE
C
*     END
C
