C
C-------------------------------------------------------------
C
C     Problem:  CHEM_REA
C
C               Right-hand side of ordinary differential equation
C               describing chemical reaction
C
C-------------------------------------------------------------
C
*     REAL CONSTANT
      R=1.987
C
*     VARIABLE
      A0, EP, M0, A, B, C, M, T
C
*     SPLINE TEMPERATURE
      0.0      291.0
      720.0    303.0
      960.0    311.0
      1200.0   316.0
      1500.0   322.0
      1800.0   331.0
      1980.0   338.0
      2160.0   346.0
      2460.0   364.0
      2640.0   365.0
      2820.0   364.0
C
*     FUNCTION MP
C
      QD0=1.0E-4*(A*TEMPERATURE(T)**2 + B*TEMPERATURE(T) + C)
      MP=-A0*DEXP(-EP/(R*TEMPERATURE(T)))*M*QD0
C
*     END
C
