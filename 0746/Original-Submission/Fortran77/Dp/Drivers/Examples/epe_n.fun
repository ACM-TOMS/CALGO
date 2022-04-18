C
C   PROBLEM:  EXPONENTIAL PARAMETER ESTIMATION
C
*     SET OF INDICES
      INDOBS=1..14
C
*     REAL CONSTANT
      EPS=1.0D-12
      OBS(I)=DEXP(I), I IN INDOBS
      TIM(I)=I, I IN INDOBS
C
*     VARIABLE
      X1,X2,X3,X4,X5,X6,X7
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
      Y1=X32*X52*X62
      IF ((Y1.LT.EPS) .AND. (Y1.GT.-EPS)) THEN
        GOTO 100
      ENDIF
      Y2=-X32*X53*X63
      IF ((Y2.LT.EPS) .AND. (Y2.GT.-EPS)) THEN
        GOTO 100
      ENDIF
      Y3=X52*X53*X65
      IF ((Y3.LT.EPS) .AND. (Y3.GT.-EPS)) THEN
        GOTO 100
      ENDIF
      Y4=-X62*X63*X65
      IF ((Y4.LT.EPS) .AND. (Y4.GT.-EPS)) THEN
         GOTO 100
      ENDIF
C
      F(I)=X1*X2*X3*
     /    (X42/Y1*DEXP(-X2*(TIM(I)-X7))
     /    +X43/Y2*DEXP(-X3*(TIM(I)-X7))
     /    +X45/Y3*DEXP(-X5*(TIM(I)-X7))
     /    +X46/Y4*DEXP(-X6*(TIM(I)-X7))) - OBS(I)
      GOTO 200
  100 CONTINUE
      F(I)=0.0
  200 CONTINUE
C
*     END
C
