C
C   PROBLEM:  EXPONENTIAL PARAMETER ESTIMATION
C
*     PARAMETER
      M=14
C
*     SET OF INDICES
      INDOBS=1..M
C
*     REAL CONSTANT
      EPS=1.0D-12
      T(I)=0.1*I, I IN INDOBS
      T(1)=0.0
C
*     TABLE Y(I), I IN INDOBS
      1    0.238
      2    0.578
      3    0.612
      4    0.650
      5    0.661
      6    0.658
      7    0.652
      8    0.649
      9    0.647
      10   0.645
      11   0.644
      12   0.644
      13   0.643
      14   0.644
C
*     VARIABLE
      X1,X2,X3,X4,X5,X6,TAU
C
*     FUNCTION F(I), I IN INDOBS
      X42=X4 - X2
      X32=X3 - X2
      X52=X5 - X2
      X62=X6 - X2
      X43=X4 - X3
      X53=X5 - X3
      X63=X6 - X3
      X45=X4 - X5
      X65=X6 - X5
      X46=X4 - X6
C
      Z1=X32*X52*X62
      IF (ABS(Z1).LT.EPS) THEN
         Z1=EPS
      ENDIF
      Z2=-X32*X53*X63
      IF (ABS(Z2).LT.EPS) THEN
         Z2=EPS
      ENDIF
      Z3=X52*X53*X65
      IF (ABS(Z3).LT.EPS) THEN
         Z3=EPS
      ENDIF
      Z4=-X62*X63*X65
      IF (ABS(Z4).LT.EPS) THEN
         Z4=EPS
      ENDIF
C
      F(I)=X1*X2*X3*
     /    (X42/Z1*DEXP(-X2*(T(I)-TAU))
     /    +X43/Z2*DEXP(-X3*(T(I)-TAU))
     /    +X45/Z3*DEXP(-X5*(T(I)-TAU))
     /    +X46/Z4*DEXP(-X6*(T(I)-TAU))) - Y(I)
C
*     END
C
