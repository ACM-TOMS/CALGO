      SUBROUTINE GAM(A, X, ACC, G, GSTAR, IFLG, IFLGST)                 GAM   10
C     LET  GAMMA(A)  DENOTE THE GAMMA FUNCTION AND  GAM(A,X)  THE
C (COMPLEMENTARY) INCOMPLETE GAMMA FUNCTION,
C
C    GAM(A,X)=INTEGRAL FROM T=X TO T=INFINITY OF EXP(-T)*T**(A-1).
C
C LET  GAMSTAR(A,X)  DENOTE TRICOMI:S FORM OF THE INCOMPLETE GAMMA
C FUNCTION, WHICH FOR A.GT.0. IS DEFINED BY
C
C  GAMSTAR(A,X)=(X**(-A)/GAMMA(A))*INTEGRAL FROM T=0 TO T=X OF
C                EXP(-T)*T**(A-1),
C
C AND FOR A.LE.0. BY ANALYTIC CONTINUATION. FOR THE PURPOSE OF
C THIS SUBROUTINE, THESE FUNCTIONS ARE NORMALIZED AS FOLLOWS&
C
C             GAM(A,X)/GAMMA(A),  IF A.GT.0.,
C
C     G(A,X)=
C
C             EXP(X)*X**(-A)*GAM(A,X),  IF A.LE.0.,
C
C     GSTAR(A,X)=(X**A)*GAMSTAR(A,X).
C
C THE PROGRAM BELOW ATTEMPTS TO EVALUATE  G(A,X)  AND  GSTAR(A,X),
C BOTH TO AN ACCURACY OF ACC SIGNIFICANT DECIMAL DIGITS, FOR ARBI-
C TRARY REAL VALUES OF  A  AND NONNEGATIVE VALUES OF  X. THE SUB-
C ROUTINE AUTOMATICALLY CHECKS FOR UNDERFLOW AND OVERFLOW CONDI-
C TIONS AND RETURNS APPROPRIATE WARNINGS THROUGH THE OUTPUT PARA-
C METERS  IFLG, IFLGST. A RESULT THAT UNDERFLOWS IS RETURNED WITH
C THE VALUE 0., ONE THAT OVERFLOWS WITH THE VALUE OF THE LARGEST
C MACHINE-REPRESENTABLE NUMBER.
C
C     NEAR LINES IN THE (A,X)-PLANE, A.LT.0., ALONG WHICH  GSTAR
C VANISHES, THE ACCURACY SPECIFIED WILL BE ATTAINED ONLY IN TERMS
C OF ABSOLUTE ERROR, NOT RELATIVE ERROR. THERE ARE OTHER (RARE)
C INSTANCES IN WHICH THE ACCURACY ATTAINED IS SOMEWHAT LESS THAN
C THE ACCURACY SPECIFIED. THE DISCREPANCY, HOWEVER, SHOULD NEVER
C EXCEED ONE OR TWO (DECIMAL) ORDERS OF ACCURACY# NO INDICATION
C OF THIS IS GIVEN THROUGH ERROR FLAGS.
C
C     PARAMETER LIST&
C
C        A - - THE FIRST ARGUMENT OF G AND GSTAR. TYPE REAL.
C        X - - THE SECOND ARGUMENT OF G AND GSTAR. TYPE REAL.
C      ACC - - THE NUMBER OF CORRECT SIGNIFICANT DECIMAL DIGITS
C              DESIRED IN THE RESULTS. TYPE REAL.
C        G - - AN OUTPUT VARIABLE RETURNING THE VALUE OF G(A,X).
C              TYPE REAL.
C    GSTAR - - AN OUTPUT VARIABLE RETURNING THE VALUE OF
C              GSTAR(A,X). TYPE REAL.
C     IFLG - - AN ERROR FLAG INDICATING A NUMBER OF ERROR CONDI-
C              TIONS IN G UPON EXECUTION. TYPE INTEGER.
C   IFLGST - - AN ERROR FLAG INDICATING A NUMBER OF ERROR CONDI-
C              TIONS IN GSTAR UPON EXECUTION. TYPE INTEGER.
C              THE VALUES OF IFLG AND IFLGST HAVE THE FOLLOWING
C              MEANINGS&
C              0 - NO ERROR CONDITION.
C              1 - ILLEGAL NEGATIVE ARGUMENT X. THE ROUTINE EXITS
C                  WITH THE VALUES ZERO FOR G AND GSTAR.
C              2 - INFINITELY LARGE RESULT AT X=0. THE ROUTINE
C                  RETURNS THE LARGEST MACHINE-REPRESENTABLE NUMBER.
C              3 - (ONLY FOR IFLGST) GSTAR IS INDETERMINATE AT
C                  A=0. AND X=0. THE ROUTINE RETURNS THE VALUE 1.,
C                  WHICH IS THE LIMIT VALUE AS X TENDS TO ZERO FOR
C                  FIXED A=0.
C              4 - THE RESULT UNDERFLOWS. IT IS SET EQUAL TO 0.
C              5 - THE RESULT OVERFLOWS. IT IS SET EQUAL TO THE
C                  LARGEST MACHINE-REPRESENTABLE NUMBER, WITH
C                  PROPER SIGN.
C              6 - CONVERGENCE FAILS WITHIN 600 ITERATIONS, EITHER
C                  IN TAYLOR:S SERIES OR IN LEGENDRE:S CONTINUED
C                  FRACTION. REASON UNKNOWN. THE COMPUTATION IS
C                  ABORTED AND THE ROUTINE RETURNS WITH ZERO
C                  VALUES FOR G AND GSTAR.
C
C     ALL MACHINE-DEPENDENT PARAMETERS ARE COLLECTED IN THE FIRST
C DATA DECLARATION. THEY ARE AS FOLLOWS&
C
C     PREC - - THE SINGLE PRECISION ACCURACY, TO BE SET APPROXI-
C              MATELY EQUAL TO BETA*ALOG(2.)/ALOG(10.), WHERE BETA
C              IS THE NUMBER OF BINARY DIGITS AVAILABLE IN THE MAN-
C              TISSA OF THE SINGLE PRECISION FLOATING-POINT WORD.
C              TYPE REAL.
C   TOPEXP - - APPROXIMATELY THE LARGEST POSITIVE NUMBER T SUCH
C              THAT 10.**T IS STILL REPRESENTABLE ON THE COMPUTER
C              IN SINGLE PRECISION FLOATING-POINT ARITHMETIC.
C              TYPE REAL.
C   BOTEXP - - APPROXIMATELY THE SMALLEST NEGATIVE NUMBER T SUCH
C              THAT 10.**T IS STILL REPRESENTABLE ON THE COMPUTER
C              IN SINGLE PRECISION FLOATING-POINT ARITHMETIC.
C              TYPE REAL.
C
C IN THE PROGRAM BELOW THESE PARAMETERS ARE SET TO CORRESPOND TO
C THE MACHINE CHARACTERISTICS OF THE CDC 6500 COMPUTER.
C
C     THE SECOND DATA DECLARATION CONTAINS THE SINGLE PRECISION
C VALUE OF ALOG(10.). THE NEXT DATA DECLARATION CONTAINS THE SUCCES-
C SIVE COEFFICIENTS IN THE MACLAURIN EXPANSION OF (1/GAMMA(A+1))-1.
C THEY ARE GIVEN TO AS MANY DECIMAL PLACES AS IS NECESSARY TO ACHIEVE
C MACHINE PRECISION (ON THE CDC 6500 COMPUTER) IN THE RANGE
C ABS(A).LE..5. MORE ACCURATE VALUES OF THESE COEFFICIENTS (TO
C 31 DECIMAL PLACES) CAN BE FOUND IN TABLE 5 OF J.W.WRENCH,JR.,
C CONCERNING TWO SERIES FOR THE GAMMA FUNCTION, MATH. COMPUT.
C 22, 1968, 617-626.
C
C     THE SUBROUTINE CALLS ON A FUNCTION SUBROUTINE, NAMED  ALGA,
C WHICH IS TO PROVIDE SINGLE PRECISION VALUES OF THE LOGARITHM OF
C THE GAMMA FUNCTION FOR ARGUMENTS LARGER THAN OR EQUAL TO .5.
C A POSSIBLE VERSION OF SUCH A FUNCTION SUBROUTINE IS APPENDED
C TO THE PRESENT SUBROUTINE. IT IS TAYLORED TO THE ACCURACY RE-
C QUIREMENTS OF THE CDC 6500 COMPUTER, AND USES RATIONAL APPROXI-
C MATIONS DUE TO CODY AND HILLSTROM (MATH. COMPUT. 21, 1967, 198-
C 203).
C
C     REFERENCE - W. GAUTSCHI, ::A COMPUTATIONAL PROCEDURE FOR
C INCOMPLETE GAMMA FUNCTIONS, ACM TRANS. MATH. SOFTWARE.
C
      DIMENSION C(18)
      DATA PREC, TOPEXP, BOTEXP /1.445E1,3.22E2,-2.93E2/
      DATA AL10 /2.30258509299405E0/
      DATA C /.577215664901533E0,-.655878071520254E0,
     * -4.2002635034095E-2,.16653861138229E0,-4.219773455554E-2,
     * -9.62197152788E-3,7.21894324666E-3,-1.1651675919E-3,
     * -2.152416741E-4,1.28050282E-4,-2.0134855E-5,-1.25049E-6,
     * 1.13303E-6,-2.0563E-7,6.12E-9,5.00E-9,-1.2E-9,1.E-10/
      G = 0.
      GSTAR = 0.
      IF (X.LT.0.) GO TO 290
C
C INITIALIZATION
C
      IFLG = 0
      IFLGST = 0
      I = 0
      IF (X.GT.0.) ALX = ALOG(X)
      ALPHA = X + .25
      IF (X.LT..25 .AND. X.GT.0.) ALPHA = ALOG(.5)/ALX
      ALPREC = AL10*PREC
      TOP = AL10*TOPEXP
      BOT = AL10*BOTEXP
      AINF = 10.**TOPEXP
      EPS = .5*10.**(-ACC)
      EPS1 = EPS/100.
      SGA = 1.
      IF (A.LT.0.) SGA = -SGA
      AE = A
      AA = ABS(A)
      AP1 = A + 1.
      AEP1 = AP1
      MA = .5 - A
      FMA = FLOAT(MA)
      AEPS = A + FMA
      SGAE = 1.
      IF (AEPS.LT.0.) SGAE = -SGAE
      AAEPS = ABS(AEPS)
      ALGP1 = 0.
C
C EVALUATION OF THE LOGARITHM OF THE ABSOLUTE VALUE OF
C GAMMA(A+1.) AND DETERMINATION OF THE SIGN OF GAMMA(A+1.)
C
      SGGA = 1.
      IF (MA.LE.0) GO TO 10
      IF (AEPS.EQ.0.) GO TO 20
      SGGA = SGAE
      IF (MA.EQ.2*(MA/2)) SGGA = -SGGA
      ALGP1 = ALGA(AEPS+1.) - ALOG(AAEPS)
      IF (MA.EQ.1) GO TO 20
      ALGP1 = ALGP1 + ALGA(1.-AEPS) - ALGA(FMA-AEPS)
      GO TO 20
   10 ALGP1 = ALGA(AP1)
   20 ALGEP1 = ALGP1
      IF (X.GT.0.) GO TO 60
C
C EVALUATION OF GSTAR(A,0.) AND G(A,0.)
C
      IF (A) 30, 40, 50
   30 IFLGST = 2
      GSTAR = AINF
      G = 1./AA
      RETURN
   40 IFLGST = 3
      GSTAR = 1.
      IFLG = 2
      G = AINF
      RETURN
   50 G = 1.
      RETURN
   60 IF (A.GT.ALPHA) GO TO 220
      IF (X.GT.1.5) GO TO 240
      IF (A.LT.-.5) GO TO 170
C
C DIRECT EVALUATION OF G(A,X) AND GSTAR(A,X) FOR X.LE.1.5
C AND -.5.LE.A.LE.ALPHA(X)
C
      GSTAR = 1.
      IF (A.GE..5) GO TO 110
   70 SUM = C(18)
      DO 80 K=1,17
        K1 = 18 - K
        SUM = AE*SUM + C(K1)
   80 CONTINUE
      GA = -SUM/(1.+AE*SUM)
      Y = AE*ALX
      IF (ABS(Y).GE.1.) GO TO 100
      SUM = 1.
      TERM = 1.
      K = 1
   90 K = K + 1
      IF (K.GT.600) GO TO 330
      TERM = Y*TERM/FLOAT(K)
      SUM = SUM + TERM
      IF (ABS(TERM).GT.EPS1*SUM) GO TO 90
      U = GA - SUM*ALX
      GO TO 120
  100 U = GA - (EXP(Y)-1.)/AE
      GO TO 120
  110 U = EXP(ALGA(A)) - (X**A)/A
  120 P = AE*X
      Q = AEP1
      R = AE + 3.
      TERM = 1.
      SUM = 1.
      K = 1
  130 K = K + 1
      IF (K.GT.600) GO TO 330
      P = P + X
      Q = Q + R
      R = R + 2.
      TERM = -P*TERM/Q
      SUM = SUM + TERM
      IF (ABS(TERM).GT.EPS1*SUM) GO TO 130
      V = (X**AEP1)*SUM/AEP1
      G = U + V
      IF (I.EQ.1) GO TO 180
      IF (A) 140, 150, 160
  140 T = EXP(X)*X**(-A)
      G = T*G
      GSTAR = 1. - A*G*EXP(-ALGP1)/T
      RETURN
  150 G = EXP(X)*G
      RETURN
  160 G = A*G*EXP(-ALGP1)
      GSTAR = 1. - G
      RETURN
C
C RECURSIVE EVALUATION OF G(A,X) FOR X.LE.1.5 AND A.LT.-.5
C
  170 I = 1
      AE = AEPS
      AEP1 = AEPS + 1.
      IF (X.LT..25 .AND. AE.GT.ALPHA) GO TO 210
      GO TO 70
  180 G = G*EXP(X)*X**(-AE)
      DO 190 K=1,MA
        G = (1.-X*G)/(FLOAT(K)-AE)
  190 CONTINUE
      ALG = ALOG(G)
C
C EVALUATION OF GSTAR(A,X) IN TERMS OF G(A,X)
C
  200 GSTAR = 1.
      IF (MA.GE.0 .AND. AEPS.EQ.0.) RETURN
      SGT = SGA*SGGA
      T = ALOG(AA) - X + A*ALX + ALG - ALGP1
      IF (T.LT.-ALPREC) RETURN
      IF (T.GE.TOP) GO TO 320
      GSTAR = 1. - SGT*EXP(T)
      RETURN
  210 I = 2
      ALGEP1 = ALGA(AEP1)
C
C EVALUATION OF GSTAR(A,X) FOR A.GT.ALPHA(X) BY TAYLOR
C EXPANSION
C
  220 G = 1.
      TERM = 1.
      SUM = 1.
      K = 0
  230 K = K + 1
      IF (K.GT.600) GO TO 340
      TERM = X*TERM/(AE+FLOAT(K))
      SUM = SUM + TERM
      IF (ABS(TERM).GT.EPS*SUM) GO TO 230
      ALGS = AE*ALX - X + ALOG(SUM) - ALGEP1
      IF (ALGS.LE.BOT) GO TO 310
      GSTAR = EXP(ALGS)
      G = 1. - GSTAR
      IF (I.NE.2) RETURN
      G = G*EXP(ALGEP1)/AE
      GO TO 180
C
C EVALUATION OF G(A,X) FOR X.GT.1.5 AND A.LE.ALPHA(X) BY
C MEANS OF THE LEGENDRE CONTINUED FRACTION
C
  240 GSTAR = 1.
      XPA = X + 1. - A
      XMA = X - 1. - A
      P = 0.
      Q = XPA*XMA
      R = 4.*XPA
      S = -A + 1.
      TERM = 1.
      SUM = 1.
      RHO = 0.
      K = 1
  250 K = K + 1
      IF (K.GT.600) GO TO 330
      P = P + S
      Q = Q + R
      R = R + 8.
      S = S + 2.
      T = P*(1.+RHO)
      RHO = T/(Q-T)
      TERM = RHO*TERM
      SUM = SUM + TERM
      IF (ABS(TERM).GT.EPS*SUM) GO TO 250
      IF (A) 260, 270, 280
  260 G = SUM/XPA
      ALG = ALOG(G)
      GO TO 200
  270 G = SUM/XPA
      RETURN
  280 ALG = A*ALX - X + ALOG(A*SUM/XPA) - ALGP1
      IF (ALG.LE.BOT) GO TO 300
      G = EXP(ALG)
      GSTAR = 1. - G
      RETURN
  290 IFLG = 1
      IFLGST = 1
      RETURN
  300 IFLG = 4
      RETURN
  310 IFLGST = 4
      RETURN
  320 IFLGST = 5
      GSTAR = -SGT*AINF
      RETURN
  330 IFLG = 6
      RETURN
  340 IFLGST = 6
      RETURN
      END
      FUNCTION ALGA(X)                                                  ALG   10
      DIMENSION CNUM(8), CDEN(8)
      DATA CNUM /4.120843185847770,85.68982062831317,243.175243524421,
     * -261.7218583856145,-922.2613728801522,-517.6383498023218,
     * -77.41064071332953,-2.208843997216182/, CDEN
     * /1.,45.64677187585908,377.8372484823942,951.323597679706,
     * 846.0755362020782,262.3083470269460,24.43519662506312,
     * .4097792921092615/
      XI = AINT(X)
      IF (X-XI.GT..5) XI = XI + 1.
      M = IFIX(XI) - 1
      XE = X
      IF (M.EQ.-1) XE = X + 1.
      IF (M.GT.0) XE = X - FLOAT(M)
      SNUM = CNUM(1)
      SDEN = CDEN(1)
      DO 10 K=2,8
        SNUM = XE*SNUM + CNUM(K)
        SDEN = XE*SDEN + CDEN(K)
   10 CONTINUE
      ALGA = (XE-1.)*SNUM/SDEN
      IF (M.GT.-1) GO TO 20
      ALGA = ALGA - ALOG(X)
      RETURN
   20 IF (M.EQ.0) RETURN
      P = XE
      IF (M.EQ.1) GO TO 40
      MM1 = M - 1
C
C THE NEXT STATEMENT IS DESIGNED TO AVOID POSSIBLE OVERFLOW IN THE
C COMPUTATION OF P. THE CONDITION IN THE IF-CLAUSE EXPRESSES THE
C INEQUALITY 1*3*5* ... *(2*M+1)/(2**M).GE.Q, WHERE Q IS THE LARGEST
C MACHINE-REPRESENTABLE NUMBER.
C
      IF (M.GE.176) GO TO 50
      DO 30 K=1,MM1
        P = (XE+FLOAT(K))*P
   30 CONTINUE
   40 ALGA = ALGA + ALOG(P)
      RETURN
   50 ALGA = ALGA + ALOG(XE)
      DO 60 K=1,MM1
        ALGA = ALGA + ALOG(XE+FLOAT(K))
   60 CONTINUE
      RETURN
      END
C                                                                       MAN1  10
C DRIVER1 - ERROR FUNCTIONS                                             MAN1  20
C                                                                       MAN1  30
C ERF X  FOR  X=0(.05)1.5  IN SINGLE AND DOUBLE PRECISION WITH          MAN1  40
C RELATIVE ERROR. CHECK AGAINST TABLE 7.1 IN NBS HANDBOOK.              MAN1  50
C                                                                       MAN1  60
      DOUBLE PRECISION DPI, DC, DX, DXSQ, DG, DGSTAR, DERF, DXMSQ, DERFC  1   70
      PI = 4.*ATAN(1.)                                                  MAN1  80
      DPI = 4.D0*DATAN(1.D0)                                            MAN1  90
      WRITE (6,99999)                                                   MAN1 100
      DO 10 I=1,31                                                      MAN1 110
        X = FLOAT(I-1)*.05                                              MAN1 120
        DX = DBLE(FLOAT(I-1))*5.D-2                                     MAN1 130
        XSQ = X*X                                                       MAN1 140
        DXSQ = DX*DX                                                    MAN1 150
        CALL GAM(.5, XSQ, 8., G, ERF, IFLG, IFLGST)                     MAN1 160
        CALL DGAM(.5D0, DXSQ, 16., DG, DERF, IFGD, IFGSTD)              MAN1 170
        ERROR = SNGL(DABS(DBLE(ERF)-DERF))                              MAN1 180
        IF (DERF.NE.0.D0) ERROR = ABS(ERROR/SNGL(DERF))                 MAN1 190
        WRITE (6,99998) X, ERF, DERF, ERROR, IFLG, IFLGST, IFGD, IFGSTD MAN1 200
   10 CONTINUE                                                          MAN1 210
C                                                                       MAN1 220
C X*EXP(X*X)*ERFC X  FOR X**(-2)=.005(.005).25  IN SINGLE AND           MAN1 230
C DOUBLE PRECISION WITH RELATIVE ERROR. CHECK AGAINST TABLE 7.3         MAN1 240
C IN NBS HANDBOOK.                                                      MAN1 250
C                                                                       MAN1 260
      WRITE (6,99997)                                                   MAN1 270
      DO 30 I=1,50                                                      MAN1 280
        XMSQ = FLOAT(I)*.005                                            MAN1 290
        DXMSQ = DBLE(FLOAT(I))*5.D-3                                    MAN1 300
        XSQ = 1./XMSQ                                                   MAN1 310
        DXSQ = 1.D0/DXMSQ                                               MAN1 320
        CALL GAM(.5, XSQ, 7., G, GSTAR, IFLG, IFLGST)                   MAN1 330
        CALL DGAM(.5D0, DXSQ, 16., DG, DGSTAR, IFGD, IFGSTD)            MAN1 340
        ERFC = 0.                                                       MAN1 350
        DERFC = 0.D0                                                    MAN1 360
        IF (XSQ.GT.740.) GO TO 20                                       MAN1 370
        ERFC = SQRT(XSQ)*EXP(XSQ)*G                                     MAN1 380
        DERFC = DSQRT(DXSQ)*DEXP(DXSQ)*DG                               MAN1 390
   20   ERROR = SNGL(DABS(DBLE(ERFC)-DERFC))                            MAN1 400
        IF (DERFC.NE.0.D0) ERROR = ABS(ERROR/SNGL(DERFC))               MAN1 410
        WRITE (6,99996) XMSQ, ERFC, DERFC, ERROR, IFLG, IFLGST, IFGD,   MAN1 420
     *   IFGSTD                                                         MAN1 430
   30 CONTINUE                                                          MAN1 440
C                                                                       MAN1 450
C ERFC(SQRT(N*PI)) FOR N=1(1)10 IN SINGLE AND DOUBLE PRECISION          MAN1 460
C WITH RELATIVE ERROR. CHECK AGAINST TABLE 7.3 IN NBS HANDBOOK.         MAN1 470
C                                                                       MAN1 480
      WRITE (6,99995)                                                   MAN1 490
      DO 40 N=1,10                                                      MAN1 500
        X = FLOAT(N)*PI                                                 MAN1 510
        DX = DBLE(FLOAT(N))*DPI                                         MAN1 520
        CALL GAM(.5, X, 8., ERFC, GSTAR, IFLG, IFLGST)                  MAN1 530
        CALL DGAM(.5D0, DX, 16., DERFC, DGSTAR, IFGD, IFGSTD)           MAN1 540
        ERROR = ABS(SNGL((DBLE(ERFC)-DERFC)/DERFC))                     MAN1 550
        WRITE (6,99994) N, ERFC, DERFC, ERROR, IFLG, IFLGST, IFGD,      MAN1 560
     *   IFGSTD                                                         MAN1 570
   40 CONTINUE                                                          MAN1 580
      STOP                                                              MAN1 590
99999 FORMAT (/26X, 1HX, 8X, 5HERF X, 14X, 5HERF X, 12X, 5HERROR, 4X,   MAN1 600
     * 4HIFLG, 1X, 6HIFLGST, 2X, 4HIFGD, 1X, 6HIFGSTD/)                 MAN1 610
99998 FORMAT (20X, E10.2, E15.7, D23.15, E10.2, 4I6)                    MAN1 620
99997 FORMAT (//23X, 7HX**(-2), 10X, 17HX*EXP(X*X)*ERFC X, 13X, 5HERROR,  1  630
     * 4X, 4HIFLG, 1X, 6HIFLGST, 2X, 4HIFGD, 1X, 6HIFGSTD/)             MAN1 640
99996 FORMAT (20X, E10.2, E14.6, D23.15, E10.2, 4I6)                    MAN1 650
99995 FORMAT (//27X, 1HN, 13X, 16HERFC(SQRT(N*PI)), 14X, 5HERROR, 4X,   MAN1 660
     * 4HIFLG, 1X, 6HIFLGST, 2X, 4HIFGD, 1X, 6HIFGSTD/)                 MAN1 670
99994 FORMAT (I28, E17.7, D23.15, E10.2, 4I6)                           MAN1 680
      END                                                               MAN1 690
C                                                                       MAN3  10
C DRIVER3 - EXPONENTIAL INTEGRAL                                        MAN3  20
C                                                                       MAN3  30
C ESUBN(X)  FOR  N=0(1)20, X=VAR,  IN SINGLE AND DOUBLE PRECISION       MAN3  40
C                                                                       MAN3  50
C ESUBN(X)  FOR  N=0(1)20, X=VAR,  IN SINGLE AND DOUBLE PRECISION       MAN3  60
C WITH RELATIVE ERROR. CHECK AGAINST TABLES I AND II IN PAGUROVA.       MAN3  70
C                                                                       MAN3  80
      DOUBLE PRECISION DX1, DX2, DX, DA, DG, DGSTAR, DESUBN, DP, DANU,  MAN3  90
     * DESBNU, DLGA                                                     MAN3 100
      DIMENSION X1(10), X2(6), DX1(10), DX2(6)                          MAN3 110
      DATA X1 /0.,.01,.05,.2,.5,1.5,5.1,10.,14.7,19.8/                  MAN3 120
      DATA X2 /.01,.37,1.44,3.02,6.57,20./                              MAN3 130
      DATA DX1 /0.D0,1.D-2,5.D-2,.2D0,.5D0,1.5D0,5.1D0,1.D1,1.47D1,     MAN3 140
     * 1.98D1/                                                          MAN3 150
      DATA DX2 /1.D-2,.37D0,1.44D0,3.02D0,6.57D0,2.D1/                  MAN3 160
      WRITE (6,99999)                                                   MAN3 170
      DO 50 I=1,10                                                      MAN3 180
        X = X1(I)                                                       MAN3 190
        DX = DX1(I)                                                     MAN3 200
        DO 40 J=1,21                                                    MAN3 210
          N = J - 1                                                     MAN3 220
          A = FLOAT(-N+1)                                               MAN3 230
          DA = DBLE(A)                                                  MAN3 240
          CALL GAM(A, X, 8., G, GSTAR, IFLG, IFLGST)                    MAN3 250
          CALL DGAM(DA, DX, 16., DG, DGSTAR, IFGD, IFGSTD)              MAN3 260
          IF (X.GT.0.) GO TO 10                                         MAN3 270
          ESUBN = 0.                                                    MAN3 280
          DESUBN = 0.D0                                                 MAN3 290
          IF (N.LE.1) GO TO 30                                          MAN3 300
          ESUBN = G                                                     MAN3 310
          DESUBN = DG                                                   MAN3 320
          GO TO 30                                                      MAN3 330
   10     IF (N.NE.0) GO TO 20                                          MAN3 340
          ESUBN = G/X                                                   MAN3 350
          DESUBN = DG/DX                                                MAN3 360
          GO TO 30                                                      MAN3 370
   20     ESUBN = EXP(-X)*G                                             MAN3 380
          DESUBN = DEXP(-DX)*DG                                         MAN3 390
   30     ERROR = SNGL(DABS(DBLE(ESUBN)-DESUBN))                        MAN3 400
          IF (DESUBN.NE.0.D0) ERROR = ABS(ERROR/SNGL(DESUBN))           MAN3 410
          WRITE (6,99998) N, A, X, ESUBN, DESUBN, ERROR, IFLG, IFLGST,  MAN3 420
     *     IFGD, IFGSTD                                                 MAN3 430
   40   CONTINUE                                                        MAN3 440
        WRITE (6,99997)                                                 MAN3 450
   50 CONTINUE                                                          MAN3 460
C                                                                       MAN3 470
C ESUBNU(X)  FOR  NU=0(.1)1., X=VAR,  IN SINGLE AND DOUBLE              MAN3 480
C PRECISION WITH RELATIVE ERROR. CHECK AGAINST TABLE III IN             MAN3 490
C PAGUROVA.                                                             MAN3 500
C                                                                       MAN3 510
      WRITE (6,99996)                                                   MAN3 520
      DO 90 I=1,6                                                       MAN3 530
        X = X2(I)                                                       MAN3 540
        DX = DX2(I)                                                     MAN3 550
        DO 80 J=1,11                                                    MAN3 560
          ANU = FLOAT(J-1)*.1                                           MAN3 570
          A = -ANU + 1.                                                 MAN3 580
          DANU = DBLE(FLOAT(J-1))*.1D0                                  MAN3 590
          DA = -DANU + 1.D0                                             MAN3 600
          CALL GAM(A, X, 7., G, GSTAR, IFLG, IFLGST)                    MAN3 610
          CALL DGAM(DA, DX, 16., DG, DGSTAR, IFGD, IFGSTD)              MAN3 620
          IF (J.LT.11) GO TO 60                                         MAN3 630
          ESUBNU = EXP(-X)*G                                            MAN3 640
          DESBNU = DEXP(-DX)*DG                                         MAN3 650
          GO TO 70                                                      MAN3 660
   60     ESUBNU = X**(-A)*EXP(ALGA(A+1.))*G/A                          MAN3 670
          DESBNU = DX**(-DA)*DEXP(DLGA(DA+1.D0))*DG/DA                  MAN3 680
   70     ERROR = SNGL(DABS(DBLE(ESUBNU)-DESBNU))                       MAN3 690
          IF (DESBNU.NE.0.D0) ERROR = ABS(ERROR/SNGL(DESBNU))           MAN3 700
          WRITE (6,99995) ANU, A, X, ESUBNU, DESBNU, ERROR, IFLG,       MAN3 710
     *     IFLGST, IFGD, IFGSTD                                         MAN3 720
   80   CONTINUE                                                        MAN3 730
        WRITE (6,99998)                                                 MAN3 740
   90 CONTINUE                                                          MAN3 750
      STOP                                                              MAN3 760
99999 FORMAT (/7X, 1HN, 5X, 1HA, 8X, 1HX, 8X, 8HESUBN(X), 10X, 6HESUBN(,  1  770
     * 2HX), 10X, 5HERROR, 5X, 4HIFLG, 1X, 6HIFLGST, 2X, 4HIFGD, 1X,    MAN3 780
     * 6HIFGSTD/)                                                       MAN3 790
99998 FORMAT (6X, I2, E10.1, E10.2, E15.7, D22.15, E10.2, 4I6)          MAN3 800
99997 FORMAT (/)                                                        MAN3 810
99996 FORMAT (//5X, 2HNU, 7X, 1HA, 8X, 1HX, 8X, 9HESUBNU(X), 9X,        MAN3 820
     * 9HESUBNU(X), 9X, 5HERROR, 5X, 4HIFLG, 1X, 6HIFLGST, 2X, 4HIFGD,  MAN3 830
     * 1X, 6HIFGSTD/)                                                   MAN3 840
99995 FORMAT (1X, 2E9.1, E10.2, E15.7, D22.15, E10.2, 4I6)              MAN3 850
      END                                                               MAN3 860
C                                                                       MAN5  10
C DRIVER5 - CHISQUARE DISTRIBUTION                                      MAN5  20
C                                                                       MAN5  30
C P(CHISQUARE,NU) AND Q(CHISQUARE,NU) FOR SELECTED VALUES OF            MAN5  40
C CHISQUARE AND NU. CHECK AGAINST TABLE 26.7 IN NBS HANDBOOK.           MAN5  50
C                                                                       MAN5  60
      DIMENSION CCHSQ(9), NUMAX(9)                                      MAN5  70
      DATA CCHSQ /.1,1.,2.,4.,6.,8.,15.,20.,60./                        MAN5  80
      DATA NUMAX /6,12,16,22,27,4*30/                                   MAN5  90
      WRITE (6,99999)                                                   MAN5 100
      DO 20 I=1,9                                                       MAN5 110
        CHSQ = CCHSQ(I)                                                 MAN5 120
        NUI = NUMAX(I)                                                  MAN5 130
        DO 10 NU=1,NUI                                                  MAN5 140
          A = .5*FLOAT(NU)                                              MAN5 150
          X = .5*CHSQ                                                   MAN5 160
          CALL GAM(A, X, 5., Q, P, IFLG, IFLGST)                        MAN5 170
          CALL GAM(A, X, 14., Q0, P0, IFLG0, IFLGS0)                    MAN5 180
          ERRP = ABS(P-P0)/P0                                           MAN5 190
          ERRQ = ABS(Q-Q0)/Q0                                           MAN5 200
          WRITE (6,99998) NU, CHSQ, X, P, Q, ERRP, ERRQ, IFLG, IFLGST   MAN5 210
   10   CONTINUE                                                        MAN5 220
        WRITE (6,99997)                                                 MAN5 230
   20 CONTINUE                                                          MAN5 240
      STOP                                                              MAN5 250
99999 FORMAT (//6X, 2HNU, 4X, 9HCHISQUARE, 7X, 1HX, 14X, 1HP, 14X, 1HQ, MAN5 260
     * 11X, 7HERROR P, 8X, 7HERROR Q, 5X, 4HIFLG, 2X, 6HIFLGST/)        MAN5 270
99998 FORMAT (5X, I3, 2E13.3, 4E15.4, 2I6)                              MAN5 280
99997 FORMAT (/)                                                        MAN5 290
      END                                                               MAN5 300
C                                                                       MAN7  10
C DRIVER7 - MOLECULAR INTEGRALS                                         MAN7  20
C ASUBN(ALPHA) FOR N=0(1)16 AND ALPHA=VAR TO 16 DECIMAL PLACES.         MAN7  30
C CHECK AGAINST TABLES IN MILLER,GERHAUSEN AND MATSEN.                  MAN7  40
C                                                                       MAN7  50
      DOUBLE PRECISION ALPHA, X, P, A, G, GSTAR, DG, DGSTAR, ASUBN,     MAN7  60
     * DASUBN                                                           MAN7  70
      DIMENSION ALPHA(10)                                               MAN7  80
      DATA ALPHA /.125D0,.5D0,1.625D0,4.25D0,7.375D0,9.875D0,12.625D0,  MAN7  90
     * 17.125D0,21.25D0,25.D0/                                          MAN7 100
      ACC = 16.                                                         MAN7 110
      DACC = 25.                                                        MAN7 120
      DO 40 I=1,10                                                      MAN7 130
        X = ALPHA(I)                                                    MAN7 140
        WRITE (6,99999) X                                               MAN7 150
        WRITE (6,99998)                                                 MAN7 160
        DO 30 J=1,17                                                    MAN7 170
          N = J - 1                                                     MAN7 180
          P = 1.D0/X                                                    MAN7 190
          IF (N.EQ.0) GO TO 20                                          MAN7 200
          DO 10 K=1,N                                                   MAN7 210
            P = DBLE(FLOAT(K))*P/X                                      MAN7 220
   10     CONTINUE                                                      MAN7 230
   20     A = DBLE(FLOAT(J))                                            MAN7 240
          CALL DGAM(A, X, ACC, G, GSTAR, IFLG, IFLGST)                  MAN7 250
          CALL DGAM(A, X, DACC, DG, DGSTAR, IFLGD, IFLGSD)              MAN7 260
          ASUBN = P*G                                                   MAN7 270
          DASUBN = P*DG                                                 MAN7 280
          ERR = ABS(SNGL((ASUBN-DASUBN)/DASUBN))                        MAN7 290
          WRITE (6,99997) N, ASUBN, DASUBN, ERR, IFLG, IFLGST, IFLGD,   MAN7 300
     *     IFLGSD                                                       MAN7 310
   30   CONTINUE                                                        MAN7 320
   40 CONTINUE                                                          MAN7 330
      STOP                                                              MAN7 340
99999 FORMAT (/5X, 6HALPHA=, D11.4//)                                   MAN7 350
99998 FORMAT (12X, 1HN, 7X, 12HASUBN(ALPHA), 17X, 12HASUBN(ALPHA), 16X, MAN7 360
     * 5HERROR, 6X, 4HIFLG, 1X, 6HIFLGST, 1X, 5HIFLGD, 1X, 6HIFLGSD/)   MAN7 370
99997 FORMAT (10X, I3, D24.15, D33.24, E15.4, 4I6)                      MAN7 380
      END                                                               MAN7 390
      SUBROUTINE DGAM(A, X, ACC, G, GSTAR, IFLG, IFLGST)                DGA   10
      DOUBLE PRECISION A, X, G, GSTAR, C, AL10, ALX, ALPHA, ALPREC,
     * TOP, BOT, AINF, EPS, EPS1, ES, SGA, AE, AA, AP1, AM1, AEP1,
     * AEM1, FMA, AEPS, SGAE, AAEPS, ALGP1, SGGA, ALGEP1, SGGS, ALGS,
     * ALG, SUM, GA, Y, TERM, U, P, Q, R, V, T, H, SGT, A1X, RHO, XPA,
     * XMA, S, DLGA
      DIMENSION C(29)
      DATA PREC, TOPEXP, BOTEXP /28.8989,322.,-293./
      DATA AL10 /2.3025850929940456840179914547D0/
      DATA C /.57721566490153286060651209008D0,-.65587807152025388107701
     * 951515D0,-4.200263503409523552900393488D-2,.166538611382291489501
     * 7007951D0,-4.21977345555443367482083013D-2,-9.6219715278769735621
     * 149217D-3,7.2189432466630995423950103D-3,-1.165167591859065112113
     * 971D-3,-2.15241674114950972815730D-4,1.2805028238811618615320D-4,
     * -2.013485478078823865569D-5,-1.2504934821426706573D-6,
     * 1.1330272319816958824D-6,-2.056338416977607103D-7,
     * 6.1160951044814158D-9,5.0020076444692229D-9,-1.181274570487020D-9
     * ,1.04342671169110D-10,7.782263439905D-12,-3.696805618642D-12,
     * 5.1003702875D-13,-2.058326054D-14,-5.34812254D-15,1.2267786D-15,
     * -1.181259D-16,1.187D-18,1.412D-18,-2.30D-19,1.7D-20/
      G = 0.D0
      GSTAR = 0.D0
      IF (X.LT.0.D0) GO TO 290
C
C INITIALIZATION
C
      IFLG = 0
      IFLGST = 0
      I = 0
      IF (X.GT.0.D0) ALX = DLOG(X)
      ALPHA = X + .25D0
      IF (X.LT..25D0 .AND. X.GT.0.D0) ALPHA = DLOG(.5D0)/ALX
      ALPREC = AL10*DBLE(PREC)
      TOP = AL10*DBLE(TOPEXP)
      BOT = AL10*DBLE(BOTEXP)
      AINF = 10.D0**TOPEXP
      EPS = .5D0*10.D0**(-ACC)
      EPS1 = EPS/1.D2
      SGA = 1.D0
      IF (A.LT.0.D0) SGA = -SGA
      AE = A
      AA = DABS(A)
      AP1 = A + 1.D0
      AEP1 = AP1
      MA = SNGL(.5D0-A)
      FMA = DBLE(FLOAT(MA))
      AEPS = A + FMA
      SGAE = 1.D0
      IF (AEPS.LT.0.D0) SGAE = -SGAE
      AAEPS = DABS(AEPS)
      ALGP1 = 0.D0
C
C EVALUATION OF THE LOGARITHM OF THE ABSOLUTE VALUE OF
C GAMMA(A+1.) AND DETERMINATION OF THE SIGN OF GAMMA(A+1.)
C
      SGGA = 1.D0
      IF (MA.LE.0) GO TO 10
      IF (AEPS.EQ.0.D0) GO TO 20
      SGGA = SGAE
      IF (MA.EQ.2*(MA/2)) SGGA = -SGGA
      ALGP1 = DLGA(AEPS+1.D0) - DLOG(AAEPS)
      IF (MA.EQ.1) GO TO 20
      ALGP1 = ALGP1 + DLGA(1.D0-AEPS) - DLGA(FMA-AEPS)
      GO TO 20
   10 ALGP1 = DLGA(AP1)
   20 ALGEP1 = ALGP1
      IF (X.GT.0.D0) GO TO 60
C
C EVALUATION OF GSTAR(A,0.) AND G(A,0.)
C
      IF (A) 30, 40, 50
   30 IFLGST = 2
      GSTAR = AINF
      G = 1.D0/AA
      RETURN
   40 IFLGST = 3
      GSTAR = 1.D0
      IFLG = 2
      G = AINF
      RETURN
   50 G = 1.D0
      RETURN
   60 IF (A.GT.ALPHA) GO TO 220
      IF (X.GT.1.5D0) GO TO 240
      IF (A.LT.-.5D0) GO TO 170
C
C DIRECT EVALUATION OF G(A,X) AND GSTAR(A,X) FOR X.LE.1.5
C AND -.5.LE.A.LE.ALPHA(X)
C
      GSTAR = 1.D0
      IF (A.GE..5D0) GO TO 110
   70 SUM = C(29)
      DO 80 K=1,28
        K1 = 29 - K
        SUM = AE*SUM + C(K1)
   80 CONTINUE
      GA = -SUM/(1.D0+AE*SUM)
      Y = AE*ALX
      IF (DABS(Y).GE.1.D0) GO TO 100
      SUM = 1.D0
      TERM = 1.D0
      K = 1
   90 K = K + 1
      IF (K.GT.600) GO TO 330
      TERM = Y*TERM/DBLE(FLOAT(K))
      SUM = SUM + TERM
      IF (DABS(TERM).GT.EPS1*SUM) GO TO 90
      U = GA - SUM*ALX
      GO TO 120
  100 U = GA - (DEXP(Y)-1.D0)/AE
      GO TO 120
  110 U = DEXP(DLGA(A)) - (X**A)/A
  120 P = AE*X
      Q = AEP1
      R = AE + 3.D0
      TERM = 1.D0
      SUM = 1.D0
      K = 1
  130 K = K + 1
      IF (K.GT.600) GO TO 330
      P = P + X
      Q = Q + R
      R = R + 2.D0
      TERM = -P*TERM/Q
      SUM = SUM + TERM
      IF (DABS(TERM).GT.EPS1*SUM) GO TO 130
      V = (X**AEP1)*SUM/AEP1
      G = U + V
      IF (I.EQ.1) GO TO 180
      IF (A) 140, 150, 160
  140 T = DEXP(X)*X**(-A)
      G = T*G
      GSTAR = 1.D0 - A*G*DEXP(-ALGP1)/T
      RETURN
  150 G = DEXP(X)*G
      RETURN
  160 G = A*G*DEXP(-ALGP1)
      GSTAR = 1.D0 - G
      RETURN
C
C RECURSIVE EVALUATION OF G(A,X) FOR X.LE.1.5 AND A.LT.-.5
C
  170 I = 1
      AE = AEPS
      AEP1 = AEPS + 1.D0
      IF (X.LT..25D0 .AND. AE.GT.ALPHA) GO TO 210
      GO TO 70
  180 G = G*DEXP(X)*X**(-AE)
      DO 190 K=1,MA
        G = (1.D0-X*G)/(DBLE(FLOAT(K))-AE)
  190 CONTINUE
      ALG = DLOG(G)
C
C EVALUATION OF GSTAR(A,X) IN TERMS OF G(A,X)
C
  200 GSTAR = 1.D0
      IF (MA.GE.0 .AND. AEPS.EQ.0.D0) RETURN
      SGT = SGA*SGGA
      T = DLOG(AA) - X + A*ALX + ALG - ALGP1
      IF (T.LT.-ALPREC) RETURN
      IF (T.GE.TOP) GO TO 320
      GSTAR = 1.D0 - SGT*DEXP(T)
      RETURN
  210 I = 2
      ALGEP1 = DLGA(AEP1)
C
C EVALUATION OF GSTAR(A,X) FOR A.GT.ALPHA(X) BY TAYLOR
C EXPANSION
C
  220 G = 1.D0
      TERM = 1.D0
      SUM = 1.D0
      K = 0
  230 K = K + 1
      IF (K.GT.600) GO TO 340
      TERM = X*TERM/(AE+DBLE(FLOAT(K)))
      SUM = SUM + TERM
      IF (DABS(TERM).GT.EPS*SUM) GO TO 230
      ALGS = AE*ALX - X + DLOG(SUM) - ALGEP1
      IF (ALGS.LE.BOT) GO TO 310
      GSTAR = DEXP(ALGS)
      G = 1.D0 - GSTAR
      IF (I.NE.2) RETURN
      G = G*DEXP(ALGEP1)/AE
      GO TO 180
C
C EVALUATION OF G(A,X) FOR X.GT.1.5 AND A.LE.ALPHA(X) BY
C MEANS OF THE LEGENDRE CONTINUED FRACTION
C
  240 GSTAR = 1.D0
      XPA = X + 1.D0 - A
      XMA = X - 1.D0 - A
      P = 0.D0
      Q = XPA*XMA
      R = 4.D0*XPA
      S = -A + 1.D0
      TERM = 1.D0
      SUM = 1.D0
      RHO = 0.D0
      K = 1
  250 K = K + 1
      IF (K.GT.600) GO TO 330
      P = P + S
      Q = Q + R
      R = R + 8.D0
      S = S + 2.D0
      T = P*(1.D0+RHO)
      RHO = T/(Q-T)
      TERM = RHO*TERM
      SUM = SUM + TERM
      IF (DABS(TERM).GT.EPS*SUM) GO TO 250
      IF (A) 260, 270, 280
  260 G = SUM/XPA
      ALG = DLOG(G)
      GO TO 200
  270 G = SUM/XPA
      RETURN
  280 ALG = A*ALX - X + DLOG(A*SUM/XPA) - ALGP1
      IF (ALG.LE.BOT) GO TO 300
      G = DEXP(ALG)
      GSTAR = 1.D0 - G
      RETURN
  290 IFLG = 1
      IFLGST = 1
      RETURN
  300 IFLG = 4
      RETURN
  310 IFLGST = 4
      RETURN
  320 IFLGST = 5
      GSTAR = -SGT*AINF
      RETURN
  330 IFLG = 6
      RETURN
  340 IFLGST = 6
      RETURN
      END
      DOUBLE PRECISION FUNCTION DLGA(DX)                                DLG   10
      DOUBLE PRECISION DBNUM, DBDEN, DX, DC, DP, DY, DT, DS
      DIMENSION DBNUM(8), DBDEN(8)
      DATA DBNUM /-3.617D3,1.D0,-6.91D2,1.D0,-1.D0,1.D0,-1.D0,1.D0/,
     * DBDEN /1.224D5,1.56D2,3.6036D5,1.188D3,1.68D3,1.26D3,3.6D2,1.2D1/
      DC = .5D0*DLOG(8.D0*DATAN(1.D0))
      DP = 1.D0
      DY = DX
      Y = SNGL(DY)
C
C THE CONDITIONAL CLAUSE IN THE NEXT STATEMENT EXPRESSES THE
C INEQUALITY Y.GT.EXP(.121189*DPREC+.053905), WHERE DPREC IS THE
C NUMBER OF DECIMAL DIGITS CARRIED IN DOUBLE PRECISION FLOATING-POINT
C ARITHMETIC.
C
   10 IF (Y.GT.35.027) GO TO 20
      DP = DY*DP
      DY = DY + 1.D0
      Y = SNGL(DY)
      GO TO 10
   20 DT = 1.D0/(DY*DY)
      DS = 4.3867D4/2.44188D5
      DO 30 I=1,8
        DS = DT*DS + DBNUM(I)/DBDEN(I)
   30 CONTINUE
      DLGA = (DY-.5D0)*DLOG(DY) - DY + DC + DS/DY - DLOG(DP)
      RETURN
      END
