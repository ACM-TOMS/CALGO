        PROGRAM SCTEST
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCC TEST PROGRAM OF HIZ, GIZ. TEST OF THE ADDITION FORMULAS:
CCC A) UNSCALED FUNCTIONS:
CCC      HI(Z+W)=SUM_K W**K/K!HI^(K)(Z)  (SAME FOR GI(Z))
CCC B) SCALED FUNCTIONS:   
CCC      HI(Z+W)=EXP((2/3)*Z**(3/2))*
CCC              EXP(-(2/3)*(Z+W)**(3/2))*
CCC              SUM_K W**K/K!HI^(K)(Z)  (SAME FOR GI(Z))
CCC
CCC  REMEMBER THAT HI(Z) IS SCALED IN THE SECTOR |ARG(Z)|<=PI/3
CCC  WHILE GI(Z) IS SCALED IN THE SECTOR  PI/3<=ARG(Z)<=PI.
CCC  BOTH THE UNSCALED AND SCALED FUNCTIONS ARE TESTED.
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCC  WARNING:  LOSS OF RELATIVE ACCURACY FOR THE TEST OF THE 
CCC            PHASE IS UNAVOIDABLE WHEN THE REAL OR IMAGINARY
CCC            PART VANISHES.
CCC            LOSS OF RELATIVE ACCURACY FOR THE TEST OF THE
CCC            MODULUS IS UNAVOIDABLE NEAR THE ZEROS OF THE
CCC            FUNCTION. 
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCC  WE TEST THE RELATIONS FOR Z=R*EXP(I*PHI) WHERE
CCC               0<=R<=30,  0<=PHI<=PI 
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
        DOUBLE PRECISION PI,X,Y,R,PH,W
        DOUBLE PRECISION ERRO1M,ERRO1P,ERRO2M,ERRO2P,ERRO3M,
     +  ERRO3P,ERRO4M,ERRO4P
        DOUBLE PRECISION EM1MAX,EM2MAX,EM3MAX,EM4MAX,EP1MAX,
     +  EP2MAX,EP3MAX,EP4MAX
        INTEGER I,J,IFAC 
        PI=3.1415926535897932385D0
        OPEN(UNIT=10,FILE='Res3',STATUS='UNKNOWN') 
        WRITE(10,*)' '
        WRITE(10,*)'TEST OF THE UNSCALED HI(Z)' 
        WRITE(10,*)'CHECK THE FOLLOWING RELATION:'
        WRITE(10,*)' '
        WRITE(10,*)' HI(Z+W)=SUM_K W**K/K!HI^(K)(Z) '  
        WRITE(10,*)' '   
        WRITE(10,30)'X','Y','MODULUS','PHASE'
        OPEN(UNIT=11,FILE='Res4',STATUS='UNKNOWN') 
        WRITE(11,*)' '
        WRITE(11,*)'TEST OF THE SCALED HI(Z)' 
        WRITE(11,*)'CHECK THE FOLLOWING RELATION:'
        WRITE(11,*)' '
        WRITE(11,*)' HI(Z+W)=EXP((2/3)*Z**(3/2))* '
        WRITE(11,*)'         EXP(-(2/3)*(Z+W)**(3/2))*'   
        WRITE(11,*)'         SUM_K W**K/K!HI^(K)(Z) '   
        WRITE(11,30)'X','Y','MODULUS','PHASE'
        OPEN(UNIT=12,FILE='Res1',STATUS='UNKNOWN') 
        WRITE(12,*)' '
        WRITE(12,*)' '
        WRITE(12,*)'TEST OF THE UNSCALED GI(Z)' 
        WRITE(12,*)'CHECK THE FOLLOWING RELATION:'
        WRITE(12,*)' '
        WRITE(12,*)' GI(Z+W)=SUM_K W**K/K!GI^(K)(Z) '  
        WRITE(12,*)' '   
        WRITE(12,30)'X','Y','MODULUS','PHASE'
        OPEN(UNIT=14,FILE='Res2',STATUS='UNKNOWN') 
        WRITE(14,*)' '
        WRITE(14,*)'TEST OF THE SCALED GI(Z)' 
        WRITE(14,*)'CHECK THE FOLLOWING RELATION:'
        WRITE(14,*)' '
        WRITE(14,*)' GI(Z+W)=EXP((2/3)*Z**(3/2))* '
        WRITE(14,*)'         EXP(-(2/3)*(Z+W)**(3/2))*'   
        WRITE(14,*)'         SUM_K W**K/K!GI^(K)(Z) '   
        WRITE(14,30)'X','Y','MODULUS','PHASE'
        EM1MAX=0.D0
        EM2MAX=0.D0
        EM3MAX=0.D0
        EM4MAX=0.D0 
        EP1MAX=0.D0
        EP2MAX=0.D0
        EP3MAX=0.D0
        EP4MAX=0.D0 
        W=0.1D0
        DO 1 I=1,300,6
          R=I*0.1D0
          DO 2 J=10,1800,6
            PH=J*PI/1800.D0
            X=R*COS(PH)
            Y=R*SIN(PH)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC          
CCC    CHECK THE UNSCALED HI(Z)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC 
            IFAC=1   
            CALL  HITEST(IFAC,X,Y,W,ERRO1M,ERRO1P)
            IF (ERRO1M.GT.EM1MAX) THEN
              EM1MAX=ERRO1M
            ENDIF
            IF (ERRO1P.GT.EP1MAX) THEN
              EP1MAX=ERRO1P
            ENDIF
            WRITE(10,33)X,Y,ERRO1M,ERRO1P    
            IF (PH.LE.PI/3.D0) THEN
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC          
CCC    CHECK THE SCALED HI(Z)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
              IFAC=2    
              CALL  HITEST(IFAC,X,Y,W,ERRO2M,ERRO2P)
              IF (ERRO2M.GT.EM2MAX) THEN
                EM2MAX=ERRO2M
              ENDIF
              IF (ERRO2P.GT.EP2MAX) THEN
                EP2MAX=ERRO2P
              ENDIF
              WRITE(11,33)X,Y,ERRO2M,ERRO2P
            ENDIF
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC          
CCC    CHECK THE UNSCALED GI(Z)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
            IFAC=1    
            CALL  GITEST(IFAC,X,Y,W,ERRO3M,ERRO3P)
            IF (ERRO3M.GT.EM3MAX) THEN
              EM3MAX=ERRO3M
            ENDIF
            IF (ERRO3P.GT.EP3MAX) THEN
              EP3MAX=ERRO3P
            ENDIF
            WRITE(12,33)X,Y,ERRO3M,ERRO3P
            IF (R.GT.0.D0) THEN    
              IF ((PH.GT.PI/3.D0).AND.(PH.LT.PI)) THEN
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC          
CCC    CHECK THE SCALED GI(Z)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
                IFAC=2    
                CALL  GITEST(IFAC,X,Y,-W,ERRO4M,ERRO4P)
                IF (ERRO4M.GT.EM4MAX) THEN
                  EM4MAX=ERRO4M
                ENDIF
                IF (ERRO4P.GT.EP4MAX) THEN
                  EP4MAX=ERRO4P
                ENDIF
                WRITE(14,33)X,Y,ERRO4M,ERRO4P
              ENDIF
            ENDIF
    2     CONTINUE
 1      CONTINUE
        WRITE(10,*)' '
        WRITE(10,*)'MAXIMUM VALUE OF THE MODULUS OF THE CHECK =',EM1MAX
        WRITE(10,*)' '
        WRITE(10,*)'MAXIMUM VALUE OF THE PHASE OF THE CHECK =',EP1MAX
        WRITE(11,*)' '
        WRITE(11,*)'MAXIMUM VALUE OF THE MODULUS OF THE CHECK =',EM2MAX
        WRITE(11,*)' '
        WRITE(11,*)'MAXIMUM VALUE OF THE PHASE OF THE CHECK =',EP2MAX
        WRITE(12,*)' '
        WRITE(12,*)'MAXIMUM VALUE OF THE MODULUS OF THE CHECK =',EM3MAX
        WRITE(12,*)' '
        WRITE(12,*)'MAXIMUM VALUE OF THE PHASE OF THE CHECK =',EP3MAX
        WRITE(14,*)' '
        WRITE(14,*)'MAXIMUM VALUE OF THE MODULUS OF THE CHECK =',EM4MAX
        WRITE(14,*)' '
        WRITE(14,*)'MAXIMUM VALUE OF THE PHASE OF THE CHECK =',EP4MAX
 30     FORMAT(6X,A1,18X,A1,12X,A5,11X,A5)                  
C31     FORMAT(6X,A1,18X,A1,12X,A6)
 33     FORMAT(D16.9,1X,D16.9,1X,D16.9,1X,D16.9)
        END   
      SUBROUTINE HITEST(IFAC,X,Y,W,ERROR1,ERROR2)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCC    TEST FOR HI(Z):
CCCC    IFAC=1: 
CCCC        THE TEST IS BASED ON THE ADDITION
CCCC        FORMULA:
CCCC         HI(Z+W)=SUM(K=0,INFTY)W**K/K!HI(K)(Z)
CCCC    IFAC=2: 
CCCC        TEST FOR NORMALIZED FUNCTIONS.
CCCC        THE ADDITION FORMULA READS NOW:
CCCC         HI(Z+W)=EXP((2/3*Z**3/2)*
CCCC                (1-(1+W/Z)**(3/2)))*
CCCC                SUM(K=0,INFTY)W**K/K!HI(K)(Z)   
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      DOUBLE PRECISION X,Y,W
      DOUBLE PRECISION EPS,D1MACH,PI,F23,RH0,RH1,RH2,RH3,IH0,IH1,
     +IH2,IH3,RES,IMS,D,T,P,Q,IMSU
      DOUBLE PRECISION PH,PHASE,R,R2,R32,TH15,TH2,S1,C1,F23R,DF1,
     +DF2,S11,C11,DEX,DRE,DIMA,RZ32,RZ2
      DOUBLE PRECISION W2,W3,RZW32,RDENO,RNUME,DMDENO,DRE2,DIMA2,
     +DRE3,DIMA3,RESU
      DOUBLE PRECISION DMODU1,DMODU2,DPHAS1,DPHAS2,ERROR1,ERROR2
      DOUBLE PRECISION IZ2,IZ32,IZW32,IDENO,INUME
      INTEGER IFAC,K,N,IERROH
      EPS=D1MACH(3)
      PI=ACOS(-1.D0)
      F23=.6666666666666667D0
      IF (IFAC.EQ.2) THEN
         PH=PHASE(X,Y)
         IF (PH.GT.PI/3.D0) THEN
           WRITE(*,*)'WARNING: THE TEST IS MEANINGLESS FOR H'
           WRITE(*,*)'Z IS OUTSIDE THE NORMALIZATION SECTOR'
         ENDIF
      ENDIF
      CALL HIZ(IFAC,X,Y,RH0,IH0,RH1,IH1,IERROH)
      IF (IFAC.EQ.1) THEN
        RH2=X*RH0-Y*IH0+1.D0/PI
        IH2=Y*RH0+X*IH0
      ELSE
        PH=PHASE(X,Y)
        R=SQRT(X*X+Y*Y)
        R2=R*R
        R32=R*SQRT(R)
        TH15=1.5D0*PH
        TH2=2.D0*PH
        S1=SIN(TH15)
        C1=COS(TH15)
        F23R=F23*R32
        DF1=F23R*C1
        DF2=F23R*S1
        S11=SIN(DF2)
        C11=COS(DF2)
        DEX=EXP(-DF1)
        DRE=DEX*C11
        DIMA=-DEX*S11
        RH2=X*RH0-Y*IH0+DRE/PI
        IH2=Y*RH0+X*IH0+DIMA/PI
        RZ32=R32*C1
        IZ32=R32*S1
        RZ2=R2*COS(TH2)
        IZ2=R2*SIN(TH2)
      ENDIF
      T=W
      RES=RH0+T*RH1
      IMS=IH0+T*IH1
      T=T*W/2.D0
      RES=RES+T*RH2
      IMS=IMS+T*IH2
      K=0
      N=3
      D=1.D0
 50   RH3=X*RH1-Y*IH1+(K+1)*RH0
      IH3=Y*RH1+X*IH1+(K+1)*IH0
      T=T*W/N
      N=N+1
      K=K+1
      P=T*RH3
      Q=T*IH3
      D=ABS(P)+ABS(Q)
      RES=RES+P
      IMS=IMS+Q
      RH0=RH1
      IH0=IH1
      RH1=RH2
      IH1=IH2
      RH2=RH3
      IH2=IH3
      IF (RES.NE.0) THEN
        IF (ABS(P/RES).GT.EPS) GOTO 50
      ENDIF
      IF (IMS.NE.0) THEN
        IF (ABS(Q/IMS).GT.EPS) GOTO 50
      ENDIF
      CALL HIZ(IFAC,X+W,Y,RH0,IH0,RH1,IH1,IERROH)
      IF (IFAC.EQ.2) THEN
        W2=W*W
        W3=W2*W
        PH=PHASE(X+W,Y)
        R=SQRT((X+W)*(X+W)+Y*Y)
        R32=R*SQRT(R)
        TH15=1.5D0*PH
        S1=SIN(TH15)
        C1=COS(TH15)
        RZW32=R32*C1
        IZW32=R32*S1
        RDENO=RZ32+RZW32
        IDENO=IZ32+IZW32
        RNUME=-(2.D0*RZ2*W+2.D0*X*W2+2.D0/3.D0*W3)
        INUME=-(2.D0*IZ2*W+2.D0*Y*W2)
        DMDENO=RDENO*RDENO+IDENO*IDENO 
        DRE2=(RNUME*RDENO+INUME*IDENO)/DMDENO
        DIMA2=(INUME*RDENO-RNUME*IDENO)/DMDENO
        DRE3=EXP(DRE2)*COS(DIMA2)
        DIMA3=EXP(DRE2)*SIN(DIMA2)
        RESU=RES*DRE3-IMS*DIMA3
        IMSU=RES*DIMA3+IMS*DRE3
        RES=RESU
        IMS=IMSU
      ENDIF
      DMODU1=ABS(RH0)+ABS(IH0)
      DMODU2=ABS(RES)+ABS(IMS)
      DPHAS1=IH0/RH0
      DPHAS2=IMS/RES
      IF (DMODU2.NE.0.D0) THEN
        ERROR1=ABS(1.D0-DMODU1/DMODU2)
      ENDIF
      IF (DPHAS2.NE.0.D0) THEN
        ERROR2=ABS(1.D0-DPHAS1/DPHAS2)
      ENDIF
      END
      SUBROUTINE GITEST(IFAC,X,Y,W,ERROR1,ERROR2)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCC    TEST FOR GI(Z):
CCCC    IFAC=1:
CCCC        THE TEST IS BASED ON THE ADDITION
CCCC        FORMULA:
CCCC         GI(Z+W)=SUM(K=0,INFTY)W**K/K!GI(K)(Z)
CCCC    IFAC=2:
CCCC        TEST FOR NORMALIZED FUNCTIONS.
CCCC        THE ADDITION FORMULA READS NOW:
CCCC         GI(Z+W)=EXP(-(2/3*Z**3/2)*
CCCC                (1-(1+W/Z)**(3/2)))*
CCCC                SUM(K=0,INFTY)W**K/K!GI(K)(Z) 
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      DOUBLE PRECISION X,Y,W
      DOUBLE PRECISION EPS,D1MACH,PI,F23,RH0,RH1,RH2,RH3,IH0,IH1,
     +IH2,IH3,RES,IMS,D,T,P,Q,IMSU
      DOUBLE PRECISION PH,PHASE,R,R2,R32,TH15,TH2,S1,C1,F23R,DF1,
     +DF2,S11,C11,DEX,DRE,DIMA,RZ32,RZ2
      DOUBLE PRECISION W2,W3,RZW32,RDENO,RNUME,DMDENO,DRE2,DIMA2,
     +DRE3,DIMA3,RESU
      DOUBLE PRECISION DMODU1,DMODU2,DPHAS1,DPHAS2,ERROR1,ERROR2
      DOUBLE PRECISION IZ2,IZ32,IZW32,IDENO,INUME
      INTEGER IFAC,K,N,IERROH
      EPS=D1MACH(3)
      PI=ACOS(-1.D0)
      F23=.6666666666666667D0
      IF (IFAC.EQ.2) THEN
         PH=PHASE(X,Y)
         IF (PH.LT.PI/3.D0) THEN
           WRITE(*,*)'WARNING: THE TEST IS MEANINGLESS FOR G'
           WRITE(*,*)'Z IS OUTSIDE THE NORMALIZATION SECTOR'
         ENDIF
      ENDIF
      CALL GIZ(IFAC,X,Y,RH0,IH0,RH1,IH1,IERROH)
      IF (IFAC.EQ.1) THEN
        RH2=X*RH0-Y*IH0-1.D0/PI
        IH2=Y*RH0+X*IH0
      ELSE
        PH=PHASE(X,Y)
        R=SQRT(X*X+Y*Y)
        R2=R*R
        R32=R*SQRT(R)
        TH15=1.5D0*PH
        TH2=2.D0*PH
        S1=SIN(TH15)
        C1=COS(TH15)
        F23R=F23*R32
        DF1=F23R*C1
        DF2=F23R*S1
        S11=SIN(DF2)
        C11=COS(DF2)
        DEX=EXP(DF1)
        DRE=DEX*C11
        DIMA=DEX*S11
        RH2=X*RH0-Y*IH0-DRE/PI
        IH2=Y*RH0+X*IH0-DIMA/PI
        RZ32=R32*C1
        IZ32=R32*S1
        RZ2=R2*COS(TH2)
        IZ2=R2*SIN(TH2)
      ENDIF
      T=W
      RES=RH0+T*RH1
      IMS=IH0+T*IH1
      T=T*W/2.D0
      RES=RES+T*RH2
      IMS=IMS+T*IH2
      K=0
      N=3
      D=1.D0
 51   RH3=X*RH1-Y*IH1+(K+1)*RH0
      IH3=Y*RH1+X*IH1+(K+1)*IH0
      T=T*W/N
      N=N+1
      K=K+1
      P=T*RH3
      Q=T*IH3
      D=ABS(P)+ABS(Q)
      RES=RES+P
      IMS=IMS+Q
      RH0=RH1
      IH0=IH1
      RH1=RH2
      IH1=IH2
      RH2=RH3
      IH2=IH3
      IF (RES.NE.0) THEN
        IF (ABS(P/RES).GT.EPS) GOTO 51
      ENDIF
      IF (IMS.NE.0) THEN
        IF (ABS(Q/IMS).GT.EPS) GOTO 51
      ENDIF
      CALL GIZ(IFAC,X+W,Y,RH0,IH0,RH1,IH1,IERROH)
      IF (IFAC.EQ.2) THEN
        W2=W*W
        W3=W2*W
        PH=PHASE(X+W,Y)
        R=SQRT((X+W)*(X+W)+Y*Y)
        R32=R*SQRT(R)
        TH15=1.5D0*PH
        S1=SIN(TH15)
        C1=COS(TH15)
        RZW32=R32*C1
        IZW32=R32*S1
        RDENO=RZ32+RZW32
        IDENO=IZ32+IZW32
        RNUME=(2.D0*RZ2*W+2.D0*X*W2+2.D0/3.D0*W3)
        INUME=(2.D0*IZ2*W+2.D0*Y*W2)
        DMDENO=RDENO*RDENO+IDENO*IDENO 
        DRE2=(RNUME*RDENO+INUME*IDENO)/DMDENO
        DIMA2=(INUME*RDENO-RNUME*IDENO)/DMDENO
        DRE3=EXP(DRE2)*COS(DIMA2)
        DIMA3=EXP(DRE2)*SIN(DIMA2)
        RESU=RES*DRE3-IMS*DIMA3
        IMSU=RES*DIMA3+IMS*DRE3
        RES=RESU
        IMS=IMSU
      ENDIF
      DMODU1=ABS(RH0)+ABS(IH0)
      DMODU2=ABS(RES)+ABS(IMS)
      DPHAS1=IH0/RH0
      DPHAS2=IMS/RES
      IF (DMODU2.NE.0.D0) THEN
        ERROR1=ABS(1.D0-DMODU1/DMODU2)
      ENDIF
      IF (DPHAS2.NE.0.D0) THEN
        ERROR2=ABS(1.D0-DPHAS1/DPHAS2)
      ENDIF
      END
