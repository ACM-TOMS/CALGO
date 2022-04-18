C
C-----------------------------------------------------------------------
C
C     Problem:  HF
C
C               Right-hand side of ordinary differential equation
C               describing a chemical engineering model
C
C-----------------------------------------------------------------------
C
*     PARAMETER
      NTEMP=7
C
*     INDEX
      K
C
*     SET OF INDICES
      IND=1..NTEMP
C
*     REAL CONSTANT
      R=1.987
      E=26920.0
      K0=1.572E-11
      IE=0.56779E-3
      TMEAN=350.0
C
*     INTEGER CONSTANT
      N=3
C
*     TABLE TIME(I), I IN IND
      1     0.0
      2  1200.0
      3  1980.0
      4  2160.0
      5  2280.0
      6  2640.0
      7  3390.0
C
*     LININT TEMP
      0.0      291.0
      1000.0   316.0
      2000.0   338.0
      2200.0   346.0
      2500.0   356.0
      5000.0   363.0
C
*     REAL CONSTANT
      KT(I)=K0*DEXP(E/(R*TEMP(TIME(I)))) , I IN IND
      PT(I)=1/(3*N*KT(I)) , I IN IND
      QT(I)=-IE/(2*N*KT(I)) , I IN IND
      WT(I)=-QT(I)/DSQRT(PT(I))**3, I IN IND
      PHI(I)=DLOG(WT(I) + DSQRT(WT(I)**2 + 1)) , I IN IND
      QD0(I)=2*DSQRT(PT(I))*DSINH(PHI(I)/3), I IN IND
C
*     VARIABLE
      A0, EP, M0, Y, T
C
*     FUNCTION YP
      VAL=0
      K=2
 4    CONTINUE
      IF (T.LE.TIME(K)) THEN
         GOTO 3
      ENDIF
      K=K+1
      IF (K.LE.NTEMP) THEN
         GOTO 4
      ENDIF
 3    CONTINUE
      K=K-1
      IF (K.LT.NTEMP) THEN
         VAL1=(TEMP(K+1)-TEMP(K))*(T-TIME(K))/(TIME(K+1)-TIME(K))
     /         + TEMP(K)
         VAL2=(QD0(K+1)-QD0(K))*(T-TIME(K))/(TIME(K+1)-TIME(K))
     /         + QD0(K)
      ELSE
         VAL1=TEMP(K)
         VAL2=QD0(K)
      ENDIF
      YP=-A0*DEXP(-EP/(R*VAL1) + EP/(R*TMEAN))*Y*VAL2
C
*     END
C
