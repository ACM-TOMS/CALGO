C     ALGORITHM 599, COLLECTED ALGORITHMS FROM ACM.
C     ALGORITHM APPEARED IN ACM-TRANS. MATH. SOFTWARE, VOL.9, NO. 2,
C     JUN., 1983, P. 255-257.
C**********************************************************************CSUN   10
C**********************************************************************CSUN   20
C**********************************************************************CSUN   30
C                                                                      CSUN   40
C                                                                      CSUN   50
C                                                                      CSUN   60
C     F O R T R A N  SOFTWARE PACKAGE FOR RANDOM NUMBER GENERATION     CSUN   70
C                                                                      CSUN   80
C                                                                      CSUN   90
C                                                                      CSUN  100
C**********************************************************************CSUN  110
C**********************************************************************CSUN  120
C**********************************************************************CSUN  130
C                                                                       SUN  140
C                                                                       SUN  150
C                                                                       SUN  160
C     CONTENTS:                                                         SUN  170
C                                                                       SUN  180
C     1) SUNIF  -  0,1 -UNIFORM DISTRIBUTION                            SUN  190
C                                                                       SUN  200
C     2) SEXPO  - (STANDARD-) EXPONENTIAL DISTRIBUTION                  SUN  210
C                                                                       SUN  220
C     3) SNORM  - (STANDARD-) NORMAL DISTRIBUTION                       SUN  230
C                                                                       SUN  240
C     4) SGAMMA - (STANDARD-) GAMMA DISTRIBUTION                        SUN  250
C                                                                       SUN  260
C     5) KPOISS - POISSON DISTRIBUTION                                  SUN  270
C                                                                       SUN  280
C                                                                       SUN  290
C     THIS PACKAGE CONSTITUTES A FORTRAN-77 DOCUMENTATION OF A SET OF   SUN  300
C     ASSEMBLER FUNCTIONS FOR SAMPLING FROM THE ABOVE DISTRIBUTIONS.    SUN  310
C     ALL ROUTINES MAKE AMPLE USE OF BINARY REPRESENTATIONS OF NUMBERS, SUN  320
C     THEY ARE AMONG THE MOST ACCURATE AND FAST SAMPLING FUNCTIONS      SUN  330
C     KNOWN. THE FORTRAN PROGRAMS BELOW YIELD THE SAME RANDOM NUMBER    SUN  340
C     SEQUENCES AS THE ONES FROM OUR ASSEMBLER PACKAGE, BUT THEY ARE    SUN  350
C     OF COURSE MUCH SLOWER (BY FACTORS 5-8 ON OUR SIEMENS 7760         SUN  360
C     COMPUTER.)                                                        SUN  370
C     THE SET OF ROUTINES WILL ALSO BE ACCEPTABLE TO FORTRAN IV         SUN  380
C     COMPILERS WHICH ALLOW DATA STATEMENTS FOR ARRAYS WITHOUT          SUN  390
C     IMPLICIT DO-LOOPS.                                                SUN  400
C                                                                       SUN  410
C                                                                       SUN  420
C     REMARKS:                                                          SUN  430
C                                                                       SUN  440
C     -  NO CARE IS TAKEN TO ENSURE THAT THE PARAMETER VALUES LIE       SUN  450
C        IN THE ALLOWED RANGE (E.G. A/MU > 0.0 FOR SGAMMA/KPOISS).      SUN  460
C                                                                       SUN  470
C     -  THE PARAMETER 'IR' MUST BE SET TO SOME  4*K+1 > 0  BEFORE      SUN  480
C        THE FIRST CALL OF ANY OF THE GENERATORS. THEREAFTER IR         SUN  490
C        MUST NOT BE ALTERED UNTIL A NEW INITIALIZATION IS DESIRED.     SUN  500
C                                                                       SUN  510
C     -  THE PACKAGE PROVIDES RANDOM DEVIATES OF 6-7 DIGITS ACCURACY.   SUN  520
C        ON MORE ACCURATE COMPUTERS THE CONSTANTS IN SEXPO, SNORM,      SUN  530
C        SGAMMA AND KPOISS OUGHT TO BE ADJUSTED ACCORDING TO LOCAL      SUN  540
C        COMMENTS OR WITH THE AID OF THE TABLES IN THE LITERATURE       SUN  550
C        QUOTED AT THE BEGINNING OF EACH FUNCTION.                      SUN  560
C                                                                       SUN  570
C                                                                       SUN  580
C**********************************************************************CSUN  590
C**********************************************************************CSUN  600
C                                                                      CSUN  610
C                                                                      CSUN  620
C       0 , 1   - U N I F O R M  DISTRIBUTION                          CSUN  630
C                                                                      CSUN  640
C                                                                      CSUN  650
C**********************************************************************CSUN  660
C**********************************************************************CSUN  670
C                                                                      CSUN  680
C     FOR DETAILS SEE:                                                 CSUN  690
C                                                                      CSUN  700
C               AHRENS, J.H., DIETER, U. AND GRUBE, A.                 CSUN  710
C               PSEUDO-RANDOM NUMBERS:  A NEW PROPOSAL                 CSUN  720
C                     FOR THE CHOICE OF MULTIPLICATORS                 CSUN  730
C               COMPUTING, 6 (1970), 121 - 138                         CSUN  740
C                                                                      CSUN  750
C**********************************************************************CSUN  760
C                                                                       SUN  770
      REAL FUNCTION SUNIF(IR)                                           SUN  780
      DOUBLE PRECISION R,FACTOR,TWO28
      SAVE R
C
C     FACTOR - INTEGER OF THE FORM 8*K+5 AS CLOSE AS POSSIBLE
C              TO  2**26 * (SQRT(5)-1)/2     (GOLDEN SECTION)
C     TWO28  = 2**28  (I.E. 28 SIGNIFICANT BITS FOR DEVIATES)
C
      DATA FACTOR /41475557.0D0/, TWO28 /268435456.0D0/
C
C     RETURNS SAMPLE U FROM THE  0,1 -UNIFORM DISTRIBUTION
C     BY A MULTIPLICATIVE CONGRUENTIAL GENERATOR OF THE FORM
C        R := R * FACTOR (MOD 1) .
C     IN THE FIRST CALL R IS INITIALIZED TO
C        R := IR / 2**28 ,
C     WHERE IR MUST BE OF THE FORM  IR = 4*K+1.
C     THEN R ASSUMES ALL VALUES  0 < (4*K+1)/2**28 < 1 DURING
C     A FULL PERIOD 2**26 OF SUNIF.
C     THE PARAMETER IR IS USED ONLY IN THE FIRST CALL FOR
C     INITIALIZATION OF SUNIF. THEREAFTER (WHEN NEGATIVE)
C     IR BECOMES A DUMMY VARIABLE.
C
      IF (IR .GE. 0) GO TO 1
C
C     STANDARD CASE:  SAMPLING
C
      R=DMOD(R*FACTOR,1.0D0)
      SUNIF=SNGL(R)
      RETURN
C
C     FIRST CALL: INITIALIZATION
C
1     R=DBLE(FLOAT(IR))/TWO28
      R=DMOD(R*FACTOR,1.0D0)
      SUNIF=SNGL(R)
      IR=-1
      RETURN
      END
C                                                                       SEX   10
C**********************************************************************CSEX   20
C**********************************************************************CSEX   30
C                                                                      CSEX   40
C                                                                      CSEX   50
C     (STANDARD-)  E X P O N E N T I A L   DISTRIBUTION                CSEX   60
C                                                                      CSEX   70
C                                                                      CSEX   80
C**********************************************************************CSEX   90
C**********************************************************************CSEX  100
C                                                                      CSEX  110
C     FOR DETAILS SEE:                                                 CSEX  120
C                                                                      CSEX  130
C               AHRENS, J.H. AND DIETER, U.                            CSEX  140
C               COMPUTER METHODS FOR SAMPLING FROM THE                 CSEX  150
C               EXPONENTIAL AND NORMAL DISTRIBUTIONS.                  CSEX  160
C               COMM. ACM, 15,10 (OCT. 1972), 873 - 882.               CSEX  170
C                                                                      CSEX  180
C     ALL STATEMENT NUMBERS CORRESPOND TO THE STEPS OF ALGORITHM       CSEX  190
C     'SA' IN THE ABOVE PAPER (SLIGHTLY MODIFIED IMPLEMENTATION)       CSEX  200
C                                                                      CSEX  210
C**********************************************************************CSEX  220
C                                                                       SEX  230
      REAL FUNCTION SEXPO(IR)                                           SEX  240
      DIMENSION Q(8)
      EQUIVALENCE (Q(1),Q1)
C
C     Q(N) = SUM(ALOG(2.0)**K/K])    K=1,..,N ,      THE HIGHEST N
C     (HERE 8) IS DETERMINED BY Q(N)=1.0 WITHIN STANDARD PRECISION
C
      DATA Q/.6931472,.9333737,.9888778,.9984959,
     ,.9998293,.9999833,.9999986,.9999999/
C
   1  A=0.0
      U=SUNIF(IR)
      GO TO 2
   3  A=A+Q1
   2  U=U+U
      IF (U.LE.1.0) GO TO 3
   4  U=U-1.0
      IF (U.GT.Q1) GO TO 6
   5  SEXPO=A+U
      RETURN
   6  I=1
      USTAR=SUNIF(IR)
      UMIN=USTAR
   7  USTAR=SUNIF(IR)
      IF (USTAR.LT.UMIN) UMIN=USTAR
   8  I=I+1
      IF (U.GT.Q(I)) GO TO 7
   9  SEXPO=A+UMIN*Q1
      RETURN
      END
C                                                                       SNO   10
C**********************************************************************CSNO   20
C**********************************************************************CSNO   30
C                                                                      CSNO   40
C                                                                      CSNO   50
C     (STANDARD-)  N O R M A L  DISTRIBUTION                           CSNO   60
C                                                                      CSNO   70
C                                                                      CSNO   80
C**********************************************************************CSNO   90
C**********************************************************************CSNO  100
C                                                                      CSNO  110
C     FOR DETAILS SEE:                                                 CSNO  120
C                                                                      CSNO  130
C               AHRENS, J.H. AND DIETER, U.                            CSNO  140
C               EXTENSIONS OF FORSYTHE'S METHOD FOR RANDOM             CSNO  150
C               SAMPLING FROM THE NORMAL DISTRIBUTION.                 CSNO  160
C               MATH. COMPUT., 27,124 (OCT. 1973), 927 - 937.          CSNO  170
C                                                                      CSNO  180
C     ALL STATEMENT NUMBERS CORRESPOND TO THE STEPS OF ALGORITHM 'FL'  CSNO  190
C     (M=5) IN THE ABOVE PAPER     (SLIGHTLY MODIFIED IMPLEMENTATION)  CSNO  200
C                                                                      CSNO  210
C**********************************************************************CSNO  220
C                                                                       SNO  230
      REAL FUNCTION SNORM(IR)                                           SNO  240
      DIMENSION A(32),D(31),T(31),H(31)
C
C     THE DEFINITIONS OF THE CONSTANTS A(K), D(K), T(K) AND
C     H(K) ARE ACCORDING TO THE ABOVEMENTIONED ARTICLE
C
      DATA A/0.0,.3917609E-1,.7841241E-1,.1177699,.1573107,
     ,.1970991,.2372021,.2776904,.3186394,.3601299,.4022501,
     ,.4450965,.4887764,.5334097,.5791322,.6260990,.6744898,
     ,.7245144,.7764218,.8305109,.8871466,.9467818,1.009990,
     ,1.077516,1.150349,1.229859,1.318011,1.417797,1.534121,
     ,1.675940,1.862732,2.153875/
      DATA D/5*0.0,.2636843,.2425085,.2255674,.2116342,.1999243,
     ,.1899108,.1812252,.1736014,.1668419,.1607967,.1553497,
     ,.1504094,.1459026,.1417700,.1379632,.1344418,.1311722,
     ,.1281260,.1252791,.1226109,.1201036,.1177417,.1155119,
     ,.1134023,.1114027,.1095039/
      DATA T/.7673828E-3,.2306870E-2,.3860618E-2,.5438454E-2,
     ,.7050699E-2,.8708396E-2,.1042357E-1,.1220953E-1,.1408125E-1,
     ,.1605579E-1,.1815290E-1,.2039573E-1,.2281177E-1,.2543407E-1,
     ,.2830296E-1,.3146822E-1,.3499233E-1,.3895483E-1,.4345878E-1,
     ,.4864035E-1,.5468334E-1,.6184222E-1,.7047983E-1,.8113195E-1,
     ,.9462444E-1,.1123001,.1364980,.1716886,.2276241,.3304980,
     ,.5847031/
      DATA H/.3920617E-1,.3932705E-1,.3950999E-1,.3975703E-1,
     ,.4007093E-1,.4045533E-1,.4091481E-1,.4145507E-1,.4208311E-1,
     ,.4280748E-1,.4363863E-1,.4458932E-1,.4567523E-1,.4691571E-1,
     ,.4833487E-1,.4996298E-1,.5183859E-1,.5401138E-1,.5654656E-1,
     ,.5953130E-1,.6308489E-1,.6737503E-1,.7264544E-1,.7926471E-1,
     ,.8781922E-1,.9930398E-1,.1155599,.1404344,.1836142,.2790016,
     ,.7010474/
C
   1  U=SUNIF(IR)
      S=0.0
      IF (U.GE.0.5) S=1.0
      U=U+U-S
   2  U=32.0*U
      I=INT(U)
      IF (I.EQ.0) GO TO 9
C
C                                START CENTER
C
   3  USTAR=U-FLOAT(I)
      AA=A(I)
   4  IF (USTAR.LE.T(I)) GO TO 5
      W=(USTAR-T(I))*H(I)
C
C                                EXIT   (BOTH CASES)
C
  17  Y=AA+W
      SNORM=Y
      IF (S.EQ.1.0) SNORM=-Y
      RETURN
C
C                                CENTER CONTINUED
C
   5  U=SUNIF(IR)
      W=U*(A(I+1)-AA)
      TT=(0.5*W+AA)*W
      GO TO 6
   8  TT=U
      USTAR=SUNIF(IR)
   6  IF (USTAR.GT.TT) GO TO 17
   7  U=SUNIF(IR)
      IF (USTAR.GE.U) GO TO 8
      USTAR=SUNIF(IR)
      GO TO 4
C
C                                START TAIL
C
   9  I=6
      AA=A(32)
      GO TO 10
  11  AA=AA+D(I)
      I=I+1
  10  U=U+U
      IF (U.LT.1.0) GO TO 11
  12  U=U-1.0
  13  W=U*D(I)
      TT=(0.5*W+AA)*W
      GO TO 14
  16  TT=U
  14  USTAR=SUNIF(IR)
      IF (USTAR.GT.TT) GO TO 17
  15  U=SUNIF(IR)
      IF (USTAR.GE.U) GO TO 16
      U=SUNIF(IR)
      GO TO 13
      END
C                                                                       SGA   10
C**********************************************************************CSGA   20
C**********************************************************************CSGA   30
C                                                                      CSGA   40
C                                                                      CSGA   50
C     (STANDARD-)  G A M M A  DISTRIBUTION                             CSGA   60
C                                                                      CSGA   70
C                                                                      CSGA   80
C**********************************************************************CSGA   90
C**********************************************************************CSGA  100
C                                                                      CSGA  110
C               PARAMETER  A >= 1.0  ]                                 CSGA  120
C                                                                      CSGA  130
C**********************************************************************CSGA  140
C                                                                      CSGA  150
C     FOR DETAILS SEE:                                                 CSGA  160
C                                                                      CSGA  170
C               AHRENS, J.H. AND DIETER, U.                            CSGA  180
C               GENERATING GAMMA VARIATES BY A                         CSGA  190
C               MODIFIED REJECTION TECHNIQUE.                          CSGA  200
C               COMM. ACM, 25,1 (JAN. 1982), 47 - 54.                  CSGA  210
C                                                                      CSGA  220
C     STEP NUMBERS CORRESPOND TO ALGORITHM 'GD' IN THE ABOVE PAPER     CSGA  230
C                                 (STRAIGHTFORWARD IMPLEMENTATION)     CSGA  240
C                                                                      CSGA  250
C**********************************************************************CSGA  260
C                                                                      CSGA  270
C               PARAMETER  0.0 < A < 1.0  ]                            CSGA  280
C                                                                      CSGA  290
C**********************************************************************CSGA  300
C                                                                      CSGA  310
C     FOR DETAILS SEE:                                                 CSGA  320
C                                                                      CSGA  330
C               AHRENS, J.H. AND DIETER, U.                            CSGA  340
C               COMPUTER METHODS FOR SAMPLING FROM GAMMA,              CSGA  350
C               BETA, POISSON AND BINOMIAL DISTRIBUTIONS.              CSGA  360
C               COMPUTING, 12 (1974), 223 - 246.                       CSGA  370
C                                                                      CSGA  380
C     (ADAPTED IMPLEMENTATION OF ALGORITHM 'GS' IN THE ABOVE PAPER)    CSGA  390
C                                                                      CSGA  400
C**********************************************************************CSGA  410
C                                                                       SGA  420
      REAL FUNCTION SGAMMA(IR,A)                                        SGA  430
C
C     INPUT:  IR=CURRENT STATE OF BASIC RANDOM NUMBER GENERATOR
C             A =PARAMETER (MEAN) OF THE STANDARD GAMMA DISTRIBUTION
C     OUTPUT: SGAMMA = SAMPLE FROM THE GAMMA-(A)-DISTRIBUTION
C
C     COEFFICIENTS Q(K) - FOR Q0 = SUM(Q(K)*A**(-K))
C     COEFFICIENTS A(K) - FOR Q = Q0+(T*T/2)*SUM(A(K)*V**K)
C     COEFFICIENTS E(K) - FOR EXP(Q)-1 = SUM(E(K)*Q**K)
C
      DATA Q1,Q2,Q3,Q4,Q5,Q6,Q7 /.04166669,.02083148,
     ,.00801191,.00144121,-.00007388,.00024511,.00024240/
      DATA A1,A2,A3,A4,A5,A6,A7 /.3333333,-.2500030,
     ,.2000062,-.1662921,.1423657,-.1367177,.1233795/
      DATA E1,E2,E3,E4,E5 /1.,.4999897,.1668290,.0407753,.0102930/
C
C     PREVIOUS A PRE-SET TO ZERO - AA IS A', AAA IS A"
C     SQRT32 IS THE SQUAREROOT OF 32 = 5.656854249492380
C
      DATA AA /0.0/, AAA /0.0/, SQRT32 /5.656854/
C
      IF (A .EQ. AA) GO TO 1
      IF (A .LT. 1.0) GO TO 12
C
C     STEP  1:  RECALCULATIONS OF S2,S,D IF A HAS CHANGED
C
      AA=A
      S2=A-0.5
      S=SQRT(S2)
      D=SQRT32-12.0*S
C
C     STEP  2:  T=STANDARD NORMAL DEVIATE,
C               X=(S,1/2)-NORMAL DEVIATE.
C               IMMEDIATE ACCEPTANCE (I)
C
   1  T=SNORM(IR)
      X=S+0.5*T
      SGAMMA=X*X
      IF (T .GE. 0.0) RETURN
C
C     STEP  3:  U= 0,1 -UNIFORM SAMPLE. SQUEEZE ACCEPTANCE (S)
C
      U=SUNIF(IR)
      IF (D*U .LE. T*T*T) RETURN
C
C     STEP  4:  RECALCULATIONS OF Q0,B,SI,C IF NECESSARY
C
      IF (A .EQ. AAA) GO TO 4
      AAA=A
      R=1.0/A
      Q0=((((((Q7*R+Q6)*R+Q5)*R+Q4)*R+Q3)*R+Q2)*R+Q1)*R
C
C               APPROXIMATION DEPENDING ON SIZE OF PARAMETER A
C               THE CONSTANTS IN THE EXPRESSIONS FOR B, SI AND
C               C WERE ESTABLISHED BY NUMERICAL EXPERIMENTS
C
      IF (A .LE. 3.686) GO TO 3
      IF (A .LE. 13.022) GO TO 2
C
C               CASE 3:  A .GT. 13.022
C
      B=1.77
      SI=.75
      C=.1515/S
      GO TO 4
C
C               CASE 2:  3.686 .LT. A .LE. 13.022
C
   2  B=1.654+.0076*S2
      SI=1.68/S+.275
      C=.062/S+.024
      GO TO 4
C
C               CASE 1:  A .LE. 3.686
C
   3  B=.463+S-.178*S2
      SI=1.235
      C=.195/S-.079+.016*S
C
C     STEP  5:  NO QUOTIENT TEST IF X NOT POSITIVE
C
   4  IF (X .LE. 0.0) GO TO 7
C
C     STEP  6:  CALCULATION OF V AND QUOTIENT Q
C
      V=T/(S+S)
      IF (ABS(V) .LE. 0.25) GO TO 5
      Q=Q0-S*T+0.25*T*T+(S2+S2)*ALOG(1.0+V)
      GO TO 6
   5  Q=Q0+0.5*T*T*((((((A7*V+A6)*V+A5)*V+A4)*V+A3)*V+A2)*V+A1)*V
C
C     STEP  7:  QUOTIENT ACCEPTANCE (Q)
C
   6  IF (ALOG(1.0-U) .LE. Q) RETURN
C
C     STEP  8:  E=STANDARD EXPONENTIAL DEVIATE
C               U= 0,1 -UNIFORM DEVIATE
C               T=(B,SI)-DOUBLE EXPONENTIAL (LAPLACE) SAMPLE
C
   7  E=SEXPO(IR)
      U=SUNIF(IR)
      U=U+U-1.0
      T=B+SIGN(SI*E,U)
C
C     STEP  9:  REJECTION IF T .LT. TAU(1) = -.71874483771719
C
      IF (T .LT. (-.7187449)) GO TO 7
C
C     STEP 10:  CALCULATION OF V AND QUOTIENT Q
C
      V=T/(S+S)
      IF (ABS(V) .LE. 0.25) GO TO 8
      Q=Q0-S*T+0.25*T*T+(S2+S2)*ALOG(1.0+V)
      GO TO 9
   8  Q=Q0+0.5*T*T*((((((A7*V+A6)*V+A5)*V+A4)*V+A3)*V+A2)*V+A1)*V
C
C     STEP 11:  HAT ACCEPTANCE (H) (IF Q NOT POSITIVE GO TO STEP 8)
C
   9  IF (Q .LE. 0.0) GO TO 7
      IF (Q .LE. 0.5) GO TO 10
      W=EXP(Q)-1.0
      GO TO 11
  10  W=((((E5*Q+E4)*Q+E3)*Q+E2)*Q+E1)*Q
C
C               IF T IS REJECTED, SAMPLE AGAIN AT STEP 8
C
  11  IF (C*ABS(U) .GT. W*EXP(E-0.5*T*T)) GO TO 7
      X=S+0.5*T
      SGAMMA=X*X
      RETURN
C
C     ALTERNATE METHOD FOR PARAMETERS A BELOW 1  (.3678794=EXP(-1.))
C
  12  AA=0.0
      B=1.0+.3678794*A
  13  P=B*SUNIF(IR)
      IF (P .GE. 1.0) GO TO 14
      SGAMMA=EXP(ALOG(P)/A)
      IF (SEXPO(IR) .LT. SGAMMA) GO TO 13
      RETURN
  14  SGAMMA=-ALOG((B-P)/A)
      IF (SEXPO(IR) .LT. (1.0-A)*ALOG(SGAMMA)) GO TO 13
      RETURN
      END
C                                                                       KPO   10
C**********************************************************************CKPO   20
C**********************************************************************CKPO   30
C                                                                      CKPO   40
C                                                                      CKPO   50
C     P O I S S O N  DISTRIBUTION                                      CKPO   60
C                                                                      CKPO   70
C                                                                      CKPO   80
C**********************************************************************CKPO   90
C**********************************************************************CKPO  100
C                                                                      CKPO  110
C     FOR DETAILS SEE:                                                 CKPO  120
C                                                                      CKPO  130
C               AHRENS, J.H. AND DIETER, U.                            CKPO  140
C               COMPUTER GENERATION OF POISSON DEVIATES                CKPO  150
C               FROM MODIFIED NORMAL DISTRIBUTIONS.                    CKPO  160
C               ACM TRANS. MATH. SOFTWARE, 8,2 (JUNE 1982), 163 - 179. CKPO  170
C                                                                      CKPO  180
C     (SLIGHTLY MODIFIED VERSION OF THE PROGRAM IN THE ABOVE ARTICLE)  CKPO  190
C                                                                      CKPO  200
C**********************************************************************CKPO  210
C                                                                       KPO  220
      INTEGER FUNCTION KPOISS(IR,MU)                                    KPO  230
C
C     INPUT:  IR=CURRENT STATE OF BASIC RANDOM NUMBER GENERATOR
C             MU=MEAN MU OF THE POISSON DISTRIBUTION
C     OUTPUT: KPOISS=SAMPLE FROM THE POISSON-(MU)-DISTRIBUTION
C
      REAL MU, MUPREV, MUOLD
C
C     MUPREV=PREVIOUS MU, MUOLD=MU AT LAST EXECUTION OF STEP P OR B.
C     TABLES: COEFFICIENTS A0-A7 FOR STEP F. FACTORIALS FACT
C     COEFFICIENTS A(K) - FOR PX = FK*V*V*SUM(A(K)*V**K)-DEL
C
      DIMENSION FACT(10), PP(35)
      DATA MUPREV,MUOLD /0.,0./
      DATA A0,A1,A2,A3,A4,A5,A6,A7 /-.5,.3333333,-.2500068,
     ,.2000118,-.1661269,.1421878,-.1384794,.1250060/
      DATA FACT /1.,1.,2.,6.,24.,120.,720.,5040.,40320.,362880./
C
C     SEPARATION OF CASES A AND B
C
      IF (MU .EQ. MUPREV) GO TO 1
      IF (MU .LT. 10.0) GO TO 12
C
C     C A S E  A. (RECALCULATION OF S,D,L IF MU HAS CHANGED)
C
      MUPREV=MU
      S=SQRT(MU)
      D=6.0*MU*MU
C
C             THE POISSON PROBABILITIES PK EXCEED THE DISCRETE NORMAL
C             PROBABILITIES FK WHENEVER K >= M(MU). L=IFIX(MU-1.1484)
C             IS AN UPPER BOUND TO M(MU) FOR ALL MU >= 10 .
C
      L=IFIX(MU-1.1484)
C
C     STEP N. NORMAL SAMPLE - SNORM(IR) FOR STANDARD NORMAL DEVIATE
C
   1  G=MU+S*SNORM(IR)
      IF (G .LT. 0.0) GO TO 2
      KPOISS=IFIX(G)
C
C     STEP I. IMMEDIATE ACCEPTANCE IF KPOISS IS LARGE ENOUGH
C
      IF (KPOISS .GE. L) RETURN
C
C     STEP S. SQUEEZE ACCEPTANCE - SUNIF(IR) FOR (0,1)-SAMPLE U
C
      FK=FLOAT(KPOISS)
      DIFMUK=MU-FK
      U=SUNIF(IR)
      IF (D*U .GE. DIFMUK*DIFMUK*DIFMUK) RETURN
C
C     STEP P. PREPARATIONS FOR STEPS Q AND H.
C             (RECALCULATIONS OF PARAMETERS IF NECESSARY)
C             .3989423=(2*PI)**(-.5)  .416667E-1=1./24.  .1428571=1./7.
C             THE QUANTITIES B1, B2, C3, C2, C1, C0 ARE FOR THE HERMITE
C             APPROXIMATIONS TO THE DISCRETE NORMAL PROBABILITIES FK.
C             C=.1069/MU GUARANTEES MAJORIZATION BY THE 'HAT'-FUNCTION.
C
   2  IF (MU .EQ. MUOLD) GO TO 3
      MUOLD=MU
      OMEGA=.3989423/S
      B1=.4166667E-1/MU
      B2=.3*B1*B1
      C3=.1428571*B1*B2
      C2=B2-15.*C3
      C1=B1-6.*B2+45.*C3
      C0=1.-B1+3.*B2-15.*C3
      C=.1069/MU
   3  IF (G .LT. 0.0) GO TO 5
C
C             'SUBROUTINE' F IS CALLED (KFLAG=0 FOR CORRECT RETURN)
C
      KFLAG=0
      GO TO 7
C
C     STEP Q. QUOTIENT ACCEPTANCE (RARE CASE)
C
   4  IF (FY-U*FY .LE. PY*EXP(PX-FX)) RETURN
C
C     STEP E. EXPONENTIAL SAMPLE - SEXPO(IR) FOR STANDARD EXPONENTIAL
C             DEVIATE E AND SAMPLE T FROM THE LAPLACE 'HAT'
C             (IF T <= -.6744 THEN PK < FK FOR ALL MU >= 10.)
C
   5  E=SEXPO(IR)
      U=SUNIF(IR)
      U=U+U-1.0
      T=1.8+SIGN(E,U)
      IF (T .LE. (-.6744)) GO TO 5
      KPOISS=IFIX(MU+S*T)
      FK=FLOAT(KPOISS)
      DIFMUK=MU-FK
C
C             'SUBROUTINE' F IS CALLED (KFLAG=1 FOR CORRECT RETURN)
C
      KFLAG=1
      GO TO 7
C
C     STEP H. HAT ACCEPTANCE (E IS REPEATED ON REJECTION)
C
   6  IF (C*ABS(U) .GT. PY*EXP(PX+E)-FY*EXP(FX+E)) GO TO 5
      RETURN
C
C     STEP F. 'SUBROUTINE' F. CALCULATION OF PX,PY,FX,FY.
C             CASE KPOISS .LT. 10 USES FACTORIALS FROM TABLE FACT
C
   7  IF (KPOISS .GE. 10) GO TO 8
      PX=-MU
      PY=MU**KPOISS/FACT(KPOISS+1)
      GO TO 11
C
C             CASE KPOISS .GE. 10 USES POLYNOMIAL APPROXIMATION
C             A0-A7 FOR ACCURACY WHEN ADVISABLE
C             .8333333E-1=1./12.  .3989423=(2*PI)**(-.5)
C
   8  DEL=.8333333E-1/FK
      DEL=DEL-4.8*DEL*DEL*DEL
      V=DIFMUK/FK
      IF (ABS(V) .LE. 0.25) GO TO 9
      PX=FK*ALOG(1.0+V)-DIFMUK-DEL
      GO TO 10
   9  PX=FK*V*V*(((((((A7*V+A6)*V+A5)*V+A4)*V+A3)*V+A2)*V+A1)*V+A0)-DEL
  10  PY=.3989423/SQRT(FK)
  11  X=(0.5-DIFMUK)/S
      XX=X*X
      FX=-0.5*XX
      FY=OMEGA*(((C3*XX+C2)*XX+C1)*XX+C0)
      IF (KFLAG) 4,4,6
C
C     C A S E  B. (START NEW TABLE AND CALCULATE P0 IF NECESSARY)
C
  12  MUPREV=0.0
      IF (MU .EQ. MUOLD) GO TO 13
      MUOLD=MU
      M=MAX0(1,IFIX(MU))
      L=0
      P=EXP(-MU)
      Q=P
      P0=P
C
C     STEP U. UNIFORM SAMPLE FOR INVERSION METHOD
C
  13  U=SUNIF(IR)
      KPOISS=0
      IF (U .LE. P0) RETURN
C
C     STEP T. TABLE COMPARISON UNTIL THE END PP(L) OF THE
C             PP-TABLE OF CUMULATIVE POISSON PROBABILITIES
C             (0.458=PP(9) FOR MU=10)
C
      IF (L .EQ. 0) GO TO 15
      J=1
      IF (U .GT. 0.458) J=MIN0(L,M)
      DO 14 K=J,L
      IF (U .LE. PP(K)) GO TO 18
  14  CONTINUE
      IF (L .EQ. 35) GO TO 13
C
C     STEP C. CREATION OF NEW POISSON PROBABILITIES P
C             AND THEIR CUMULATIVES Q=PP(K)
C
  15  L=L+1
      DO 16 K=L,35
      P=P*MU/FLOAT(K)
      Q=Q+P
      PP(K)=Q
      IF (U .LE. Q) GO TO 17
  16  CONTINUE
      L=35
      GO TO 13
  17  L=K
  18  KPOISS=K
      RETURN
      END
C**********************************************************************CMAN   10
C**********************************************************************CMAN   20
C**********************************************************************CMAN   30
C                                                                      CMAN   40
C                                                                      CMAN   50
C                                                                      CMAN   60
C     AUTOMATIC TEST DRIVER                                            CMAN   70
C     FOR                                                              CMAN   80
C     RANDOM NUMBER PACKAGE AHRENS/DIETER/KOHRT                        CMAN   90
C     (THIS DRIVER MAY CAUSE MANY EXP() UNDERFLOWS.)                    MAN  100
C                                                                      CMAN  110
C                                                                      CMAN  120
C                                                                      CMAN  130
C**********************************************************************CMAN  140
C**********************************************************************CMAN  150
C**********************************************************************CMAN  160
C     PROGRAM TEST                                                      MAN  170
      REAL MU                                                           MAN  180
      DIMENSION VPAR4(22),VPAR5(24),SAMPLE(10000),II(5)                 MAN  190
C                                                                       MAN  200
C     VPAR4 - VECTOR OF PARAMETER VALUES FOR CASE 4 : SGAMMA            MAN  210
C                                                                       MAN  220
      DATA VPAR4 /.0001,.25,.5,.75,.9999,1.,1.5,2.,                     MAN  230
     ,3.,4.,5.,7.,10.,15.,20.,30.,50.,100.,1000.,                       MAN  240
     ,10000.,100000.,1000000./                                          MAN  250
C                                                                       MAN  260
C     VPAR5 - VECTOR OF PARAMETER VALUES FOR CASE 5 : KPOISS            MAN  270
C                                                                       MAN  280
      DATA VPAR5 /.0001,1.,2.,5.,9.99,10.,                              MAN  290
     ,12.,15.,20.,25.,30.,40.,50.,75.,100.,150.,                        MAN  300
     ,200.,500.,1000.,2000.,5000.,1.E4,1.E5,1.E6/                       MAN  310
C                                                                       MAN  320
C     II - NUMBER OF RUNS FOR EACH DISTRIBUTION                         MAN  330
C                                                                       MAN  340
      DATA II /1,1,1,22,24/                                             MAN  350
C                                                                       MAN  360
C     FORMAT STATEMENTS                                                 MAN  370
C                                                                       MAN  380
   1  FORMAT(' ',/,' ',/,' LISTING OF TRIAL RUNS',                      MAN  390
     ,' FOR RANDOM NUMBER PACKAGE AHRENS/DIETER/KOHRT',                 MAN  400
C    ,/,' ===============================',                             MAN  410
     ,'====================================')                           MAN  420
   2  FORMAT('   FIRST 100 SAMPLES:',/,                                 MAN  430
     ,'   ..................',/,' ',/,(5E15.6))                         MAN  440
   3  FORMAT(' ',/,' ',/,'   TEST DATA:',                               MAN  450
     ,'     (]BASED ON 10000 SAMPLES])',                                MAN  460
     ,/,'   ..........',/,20X,'MEAN',11X,                               MAN  470
     ,'STD.DEV.',7X,'SKEWNESS',                                         MAN  480
     ,7X,'EXCESS',/,' ',/,                                              MAN  490
     ,'   TRUE VALUES:',4E15.6,/,                                       MAN  500
     ,'   SAMPLE DATA:',4E15.6,/,' ')                                   MAN  510
   4  FORMAT(' ',/,' 1.)   0,1 -UNIFORM DISTRIBUTION:',                 MAN  520
     ,/,' ********************************')                            MAN  530
   5  FORMAT(' ',/,' 2.)  (STANDARD-) EXPONENTIAL DISTRIBUTION:',       MAN  540
     ,/,' ******************************************')                  MAN  550
   6  FORMAT(' ',/,' 3.)  (STANDARD-) NORMAL DISTRIBUTION:',            MAN  560
     ,/,' ******************************************')                  MAN  570
   7  FORMAT(' ',/,' 4.)  (STANDARD-) GAMMA-(A) DISTRIBUTION:',         MAN  580
     ,/,' ****************************************')                    MAN  590
   8  FORMAT(' ',/,' 5.)  POISSON-(MU) DISTRIBUTION:',                  MAN  600
     ,/,' *******************************',/,' ',                       MAN  610
     ,/,'   (INTEGER SAMPLES ARE DISPLAYED AS REALS])',                 MAN  620
     ,/,' ',/,' ')                                                      MAN  630
   9  FORMAT(43X,'    GAMMA-(A):  A =',E13.6,                           MAN  640
     ,/,43X,'    ----------------------------')                         MAN  650
  10  FORMAT(43X,'POISSON-(MU):  MU =',E13.6,                           MAN  660
     ,/,43X,'--------------------------------')                         MAN  670
  11  FORMAT(' ',/,' ')                                                 MAN  680
C                                                                       MAN  690
C     DEFINE OUTPUT UNIT NUMBER                                         MAN  700
C                                                                       MAN  710
      NOUT=10                                                           MAN  720
C                                                                       MAN  730
C     OUTPUT: MAIN HEADING                                              MAN  740
C                                                                       MAN  750
      WRITE (NOUT,1)                                                    MAN  760
C                                                                       MAN  770
C     TRIAL RUNS FOR 5 DIFFERENT CASES:                                 MAN  780
C                                                                       MAN  790
C       NDIS=1 :   0,1 -UNIFORM DISTRIBUTION                            MAN  800
C       NDIS=2 :  (STANDARD-) EXPONENTIAL DISTRIBUTION                  MAN  810
C       NDIS=3 :  (STANDARD-) NORMAL DISTRIBUTION                       MAN  820
C       NDIS=4 :  (STANDARD-) GAMMA-(A) DISTRIBUTION                    MAN  830
C       NDIS=5 :  POISSON-(MU) DISTRIBUTION                             MAN  840
C                                                                       MAN  850
      DO 27 NDIS=1,5                                                    MAN  860
      NRUN=II(NDIS)                                                     MAN  870
C                                                                       MAN  880
C     OUTPUT: CASE HEADING                                              MAN  890
C                                                                       MAN  900
      WRITE (NOUT,11)                                                   MAN  910
      IF (NDIS .EQ. 1) WRITE (NOUT,4)                                   MAN  920
      IF (NDIS .EQ. 2) WRITE (NOUT,5)                                   MAN  930
      IF (NDIS .EQ. 3) WRITE (NOUT,6)                                   MAN  940
      IF (NDIS .EQ. 4) WRITE (NOUT,7)                                   MAN  950
      IF (NDIS .EQ. 5) WRITE (NOUT,8)                                   MAN  960
C                                                                       MAN  970
C     EACH CASE: ONE RUN FOR EVERY PARAMETER VALUE                      MAN  980
C                                                                       MAN  990
      DO 26 NPAR=1,NRUN                                                 MAN 1000
C                                                                       MAN 1010
C     CASE 4 AND 5: SET PARAMETER VALUES ACCORDING TO DATA VECTOR       MAN 1020
C                                                                       MAN 1030
      IF (NDIS .EQ. 4) A=VPAR4(NPAR)                                    MAN 1040
      IF (NDIS .EQ. 5) MU=VPAR5(NPAR)                                   MAN 1050
C                                                                       MAN 1060
C     SEED FOR UNIFORM RANDOM NUMBER GENERATOR IS INITIALIZED TO 4*0+1  MAN 1070
C                                                                       MAN 1080
      IR=1                                                              MAN 1090
C                                                                       MAN 1100
C     EACH CASE SEPARATELY: SAMPLING AND TEST DATA                      MAN 1110
C        T2 - STANDARD DEVIATION (=SQRT(VARIANCE)),                     MAN 1120
C        T1 - MEAN,   T3 - SKEWNESS,   T4 - EXCESS.                     MAN 1130
C                                                                       MAN 1140
      GO TO (12,14,16,18,20) , NDIS                                     MAN 1150
C                                                                       MAN 1160
C     CASE 1 :   0,1 -UNIFORM DISTRIBUTION                              MAN 1170
C                                                                       MAN 1180
  12  DO 13 I=1,10000                                                   MAN 1190
  13  SAMPLE(I)=SUNIF(IR)                                               MAN 1200
      T1=0.5                                                            MAN 1210
      T2=1.0/SQRT(12.0)                                                 MAN 1220
      T3=0.0                                                            MAN 1230
      T4=-1.2                                                           MAN 1240
      GO TO 23                                                          MAN 1250
C                                                                       MAN 1260
C     (STANDARD-) EXPONENTIAL DISTRIBUTION                              MAN 1270
C                                                                       MAN 1280
  14  DO 15 I=1,10000                                                   MAN 1290
  15  SAMPLE(I)=SEXPO(IR)                                               MAN 1300
      T1=1.0                                                            MAN 1310
      T2=1.0                                                            MAN 1320
      T3=2.0                                                            MAN 1330
      T4=6.0                                                            MAN 1340
      GO TO 23                                                          MAN 1350
C                                                                       MAN 1360
C     (STANDARD-) NORMAL DISTRIBUTION                                   MAN 1370
C                                                                       MAN 1380
  16  DO 17 I=1,10000                                                   MAN 1390
  17  SAMPLE(I)=SNORM(IR)                                               MAN 1400
      T1=0.0                                                            MAN 1410
      T2=1.0                                                            MAN 1420
      T3=0.0                                                            MAN 1430
      T4=0.0                                                            MAN 1440
      GO TO 23                                                          MAN 1450
C                                                                       MAN 1460
C     (STANDARD-) GAMMA-(A) DISTRIBUTION                                MAN 1470
C                                                                       MAN 1480
  18  DO 19 I=1,10000                                                   MAN 1490
  19  SAMPLE(I)=SGAMMA(IR,A)                                            MAN 1500
      T1=A                                                              MAN 1510
      T2=SQRT(A)                                                        MAN 1520
      T3=2.0/T2                                                         MAN 1530
      T4=6.0/A                                                          MAN 1540
      GO TO 22                                                          MAN 1550
C                                                                       MAN 1560
C     POISSON-(MU) DISTRIBUTION                                         MAN 1570
C                                                                       MAN 1580
  20  DO 21 I=1,10000                                                   MAN 1590
  21  SAMPLE(I)=FLOAT(KPOISS(IR,MU))                                    MAN 1600
      T1=MU                                                             MAN 1610
      T2=SQRT(MU)                                                       MAN 1620
      T3=1.0/T2                                                         MAN 1630
      T4=1.0/MU                                                         MAN 1640
C                                                                       MAN 1650
C     CASE 4 AND 5:  OUTPUT: PARAMETER VALUE                            MAN 1660
C                                                                       MAN 1670
  22  IF (NPAR .NE. 1) WRITE (NOUT,11)                                  MAN 1680
      IF (NDIS .EQ. 4) WRITE (NOUT, 9) A                                MAN 1690
      IF (NDIS .EQ. 5) WRITE (NOUT,10) MU                               MAN 1700
C                                                                       MAN 1710
  23  CONTINUE                                                          MAN 1720
C                                                                       MAN 1730
C     OUTPUT : FIRST 100 RANDOM DEVIATES FOR EACH RUN                   MAN 1740
C              (INTEGER SAMPLES ARE DISPLAYED AS REALS])                MAN 1750
C                                                                       MAN 1760
      IF (NDIS .LE. 3) WRITE (NOUT,11)                                  MAN 1770
      WRITE (NOUT,2) (SAMPLE(I),I=1,100)                                MAN 1780
C                                                                       MAN 1790
C     EVALUATION OF SAMPLE MEAN:    E1 - (1/N)*SUM(SAMPLE(I))           MAN 1800
C                                                                       MAN 1810
      S1=0.0                                                            MAN 1820
      DO 24 I=1,10000                                                   MAN 1830
  24  S1=S1+SAMPLE(I)                                                   MAN 1840
      E1=S1/10000.0                                                     MAN 1850
C                                                                       MAN 1860
C     EVALUATION OF FURTHER SAMPLE ESTIMATES :                          MAN 1870
C                                                                       MAN 1880
C     SK       - (1/N)*SUM((SAMPLE(I)-E1)**K)                           MAN 1890
C                SAMPLE CENTRAL MOMENTS  (K=2,3,4)                      MAN 1900
C                WITH RESPECT TO SAMPLE MEAN                            MAN 1910
C                                                                       MAN 1920
C     E2       - SQRT(S2)                                               MAN 1930
C                SAMPLE STANDARD DEVIATION                              MAN 1940
C                (=SQRT(SAMPLE VARIANCE))                               MAN 1950
C                                                                       MAN 1960
C     E3       - S3/S2**(3/2)                                           MAN 1970
C                SAMPLE SKEWNESS                                        MAN 1980
C                                                                       MAN 1990
C     E4       - S4/S2**2-3                                             MAN 2000
C                SAMPLE EXCESS                                          MAN 2010
C                                                                       MAN 2020
      S2=0.0                                                            MAN 2030
      S3=0.0                                                            MAN 2040
      S4=0.0                                                            MAN 2050
C                                                                       MAN 2060
      DO 25 I=1,10000                                                   MAN 2070
C                                                                       MAN 2080
      X1=SAMPLE(I)-E1                                                   MAN 2090
      X2=X1*X1                                                          MAN 2100
      X3=X2*X1                                                          MAN 2110
      X4=X2*X2                                                          MAN 2120
C                                                                       MAN 2130
      S2=S2+X2                                                          MAN 2140
      S3=S3+X3                                                          MAN 2150
      S4=S4+X4                                                          MAN 2160
C                                                                       MAN 2170
  25  CONTINUE                                                          MAN 2180
C                                                                       MAN 2190
      S2=S2/10000.0                                                     MAN 2200
      S3=S3/10000.0                                                     MAN 2210
      S4=S4/10000.0                                                     MAN 2220
C                                                                       MAN 2230
      E2=SQRT(S2)                                                       MAN 2240
      E3=S3/SQRT(S2*S2*S2)                                              MAN 2250
      E4=S4/(S2*S2)-3.0                                                 MAN 2260
C                                                                       MAN 2270
C     END OF EVALUATION                                                 MAN 2280
C     OUTPUT: CHARACTERISTIC DATA                                       MAN 2290
C                                                                       MAN 2300
      WRITE (NOUT,3) T1,T2,T3,T4,E1,E2,E3,E4                            MAN 2310
C                                                                       MAN 2320
C     END OF PARAMETER LOOP                                             MAN 2330
C                                                                       MAN 2340
  26  CONTINUE                                                          MAN 2350
C                                                                       MAN 2360
C     END OF DISTRIBUTION LOOP                                          MAN 2370
C                                                                       MAN 2380
  27  CONTINUE                                                          MAN 2390
C                                                                       MAN 2400
C     END OF PROGRAM                                                    MAN 2410
C                                                                       MAN 2420
      STOP                                                              MAN 2430
      END                                                               MAN 2440
