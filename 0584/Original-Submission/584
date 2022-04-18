C     ALGORITHM 584, COLLECTED ALGORITHMS FROM ACM.
C     ALGORITHM APPEARED IN ACM-TRANS. MATH. SOFTWARE, VOL.8, NO. 2,
C     JUN., 1982, P. 210.
C     PROGRAM KONYN(OUTPUT,TAPE6=OUTPUT)                                MAN   10
C     DRIVER MAIN PROGRAM FOR TESTING SUBROUTINE CUBTRI.  SEVEN         MAN   20
C     INTEGRALS ARE EVALUATED AT TOLERANCES RANGING FROM                MAN   30
C     1E-3 THROUGH 1E-9, TAKING ADVANTAGE OF THE RESTARTABILITY FEATURE.MAN   40
      DIMENSION ACTUAL(7)                                               MAN   50
      REAL ANS, EPS, ERR, T(2,3), W(5000), NINE, RDATA(1)               MAN   60
      INTEGER IDATA(1)                                                  MAN   70
C                                                                       MAN   80
      EXTERNAL F                                                        MAN   90
C                                                                       MAN  100
C     FOR THE FIRST THREE EXAMPLES,T IS THE TRIANGLE CORRESPONDING TO   MAN  110
C     INTEGRAL (0 TO 1)(0 TO X) F(X,Y,IDATA,RDATA) DY DX                MAN  120
C                                                                       MAN  130
      DATA ZERO /0.E0/, ONE /1.E0/, TWO /2.E0/, THREE /3.E0/, SIX /6.E0/MAN  140
      DATA NINE /9.E0/, YTIEN /18.E0/, POINT1 /.1E0/, POINT5 /.5E0/,    MAN  150
     * P002 /1.E-2/                                                     MAN  160
      T(1,1) = ZERO                                                     MAN  170
      T(1,2) = ONE                                                      MAN  180
      T(1,3) = ONE                                                      MAN  190
      T(2,1) = ZERO                                                     MAN  200
      T(2,2) = ZERO                                                     MAN  210
      T(2,3) = ONE                                                      MAN  220
C                                                                       MAN  230
C     SET CORRECT ANSWERS                                               MAN  240
C                                                                       MAN  250
      S = SQRT(THREE)                                                   MAN  260
      ACTUAL(1) = (ALOG(TWO+S))/S                                       MAN  270
      ACTUAL(2) = ACTUAL(1)                                             MAN  280
      ACTUAL(3) = ATAN(ONE)/TWO                                         MAN  290
C       ACTUAL(4) IS THE INTEGRAL (0 TO 1)(0 TO X) OF                   MAN  300
C       EXP(SIN(X)*COS(Y)) DY DX. THIS VALUE, TO 16 DIGITS,IS GIVEN     MAN  310
C       BY A.HAEGEMANS, COMPUTING 19(1977),179-187.                     MAN  320
      ACTUAL(4) = .6918104506612316                                     MAN  330
      ACTUAL(5) = POINT5 - ATAN(ONE) + POINT5*ALOG(TWO)                 MAN  340
      ACTUAL(6) = SIN(THREE)/NINE - SIN(SIX)/YTIEN                      MAN  350
      ACTUAL(7) = TWO/THREE                                             MAN  360
C                                                                       MAN  370
      NFILE = 6                                                         MAN  380
      WRITE (NFILE,99999)                                               MAN  390
      WRITE (NFILE,99998)                                               MAN  400
C                                                                       MAN  410
      DO 40 IEXAMP=1,7                                                  MAN  420
        IDATA(1) = IEXAMP                                               MAN  430
        IF (IEXAMP.NE.4) GO TO 10                                       MAN  440
C                                                                       MAN  450
C       FOR LAST FOUR EXAMPLES,T IS THE TRIANGLE CORRESPONDING TO       MAN  460
C       INTEGRAL (0 TO 1)(0 TO 1-X) F(X,Y,IDATA,RDATA) DY DX            MAN  470
C                                                                       MAN  480
        T(1,2) = ZERO                                                   MAN  490
        T(2,1) = ONE                                                    MAN  500
        T(2,3) = ZERO                                                   MAN  510
   10   EPS = P002                                                      MAN  520
        NCALLS = 0                                                      MAN  530
        MCALLS = 20000                                                  MAN  540
        NW = 5000                                                       MAN  550
        IER = 0                                                         MAN  560
        DO 20 I=1,7                                                     MAN  570
          EPS = EPS*POINT1                                              MAN  580
          IF (IER.GE.4) GO TO 30                                        MAN  590
C                                                                       MAN  600
C          DO NOT GO ON AFTER ROUNDOFF ERRORS HAVE BEEN SIGNALLED       MAN  610
C                                                                       MAN  620
          CALL CUBTRI(F, T, EPS, MCALLS, ANS, ERR, NCALLS, W, NW,       MAN  630
     *     IDATA, RDATA, IER)                                           MAN  640
          ACTU = ABS(ACTUAL(IEXAMP)-ANS)                                MAN  650
          WRITE (NFILE,99997) EPS, ANS, NCALLS, ERR, ACTU, IER          MAN  660
   20   CONTINUE                                                        MAN  670
   30   WRITE (NFILE,99996)                                             MAN  680
   40 CONTINUE                                                          MAN  690
      STOP                                                              MAN  700
C                                                                       MAN  710
99999 FORMAT (1H1, 41X, 6HCUBTRI/42X, 6(1H-)///)                        MAN  720
99998 FORMAT (2X, 8HREQUIRED, 3X, 13HAPPROXIMATION, 5X, 9HNUMBER OF,    MAN  730
     * 5X, 9HESTIMATED, 9X, 6HACTUAL, 12X, 11HTERMINATION/2X, 7HTOLERAN,MAN  740
     * 2HCE, 2X, 11HOF INTEGRAL, 7X, 5HCALLS, 9X, 14HABSOLUTE ERROR,    MAN  750
     * 4X, 14HABSOLUTE ERROR, 4X, 9HINDICATOR/2X, 9(1H-), 2X, 13(1H-),  MAN  760
     * 5X, 9(1H-), 5X, 14(1H-), 4X, 14(1H-), 4X, 11(1H-)/)              MAN  770
99997 FORMAT (3X, 1PE6.0, 3X, 0PF15.12, 3X, I8, 2(10X, 1PE7.1), 5X, I6) MAN  780
99996 FORMAT (//)                                                       MAN  790
      END                                                               MAN  800
      FUNCTION F(X, Y, IDATA, RDATA)                                    F     10
C
      DIMENSION IDATA(1), RDATA(1)
C
      DATA ZERO /0.E0/, ONE /1.E0/, THREE /3.E0/, FOUR /4.E0/, PTWO5
     * /.25E0/
      DATA POINT5 /.5E0/, SIX /6.E0/
      DATA EENE6 /1E-6/, AGTE6 /8E6/, VYFE6 /5E6/
C
      IFUNCT = IDATA(1)
      GO TO (10, 20, 30, 50, 60, 70, 80), IFUNCT
C
   10 F = ONE/SQRT(X*X+THREE*Y*Y)
      RETURN
C
   20 PI = FOUR*ATAN(ONE)
      RAND = EENE6*COS(AMOD(AGTE6*X+VYFE6*Y,PI))
      F = ONE/SQRT(X*X+THREE*Y*Y)*(ONE+RAND)
      RETURN
C
   30 VRAAG = (X-POINT5)**2 + (Y-POINT5)**2
      IF (VRAAG.LE.PTWO5) GO TO 40
      F = ZERO
      RETURN
   40 F = ONE
      RETURN
C
   50 F = EXP(SIN(X)*COS(Y))
      RETURN
C
   60 XX = X*X
      F = XX/(ONE+XX)
      RETURN
C
   70 F = SIN(THREE*X+SIX*Y)
      RETURN
C
   80 F = ONE/SQRT(X+Y)
      RETURN
C
      END
      SUBROUTINE CUBTRI(F, T, EPS, MCALLS, ANS, ERR, NCALLS, W, NW,     CUB   10
     * IDATA, RDATA, IER)
C
C       ADAPTIVE CUBATURE OVER A TRIANGLE
C
C       PARAMETERS
C          F     - USER SUPPLIED EXTERNAL FUNCTION OF THE FORM
C                  F(X,Y,IDATA,RDATA)
C                  WHERE X AND Y ARE THE CARTESIAN COORDINATES OF A
C                  POINT IN THE PLANE, AND IDATA AND RDATA ARE INTEGER
C                  AND REAL VECTORS IN WHICH DATA MAY BE PASSED.
C          T     - ARRAY OF DIMENSION (2,3) WHERE T(1,J) AND T(2,J)
C                  ARE THE X AND Y COORDINATES OF THE J-TH VERTEX OF
C                  THE GIVEN TRIANGLE (INPUT)
C          EPS   - REQUIRED TOLERANCE (INPUT).  IF THE COMPUTED
C                  INTEGRAL IS BETWEEN-1 AND 1, AN ABSOLUTE ERROR
C                  TEST IS USED, ELSE A RELATIVE ERROR TEST IS USED.
C          MCALLS- MAXIMUM PERMITTED NUMBER OF CALLS TO F (INPUT)
C          ANS   - ESTIMATE FOR THE VALUE OF THE INTEGRAL OF F OVER
C                  THE GIVEN TRIANGLE (OUTPUT)
C          ERR   - ESTIMATED ABSOLUTE ERROR IN ANS (OUTPUT)
C          NCALLS- ACTUAL NUMBER OF CALLS TO F (OUTPUT).  THIS
C                  PARAMETER MUST BE INITIALIZED TO 0 ON THE FIRST
C                  CALL TO CUBTRI FOR A GIVEN INTEGRAL (INPUT)
C          W     - WORK SPACE.  MAY NOT BE DESTROYED BETWEEN CALLS TO
C                  CUBTRI IF RESTARTING IS INTENDED
C          NW    - LENGTH OF WORK SPACE (INPUT).
C                  IF NW .GE. 3*(19+3*MCALLS)/38, TERMINATION DUE TO
C                  FULL WORK SPACE WILL NOT OCCUR.
C          IER   - TERMINATION INDICATOR (OUTPUT)
C                  IER=0   NORMAL TERMINATION, TOLERANCE SATISFIED
C                  IER=1   MAXIMUM NUMBER OF CALLS REACHED
C                  IER=2   WORK SPACE FULL
C                  IER=3   FURTHER SUBDIVISION OF TRIANGLES IMPOSSIBLE
C                  IER=4   NO FURTHER IMPROVEMENT IN ACCURACY IS
C                        POSSIBLE DUE TO ROUNDING ERRORS IN FUNCTION
C                        VALUES
C                  IER=5   NO FURTHER IMPROVEMENT IN ACCURACY IS
C                        POSSIBLE BECAUSE SUBDIVISION DOES NOT
C                        CHANGE THE ESTIMATED INTEGRAL. MACHINE
C                        ACCURACY HAS PROBABLY BEEN REACHED BUT
C                        THE ERROR ESTIMATE IS NOT SHARP ENOUGH.
C
C       CUBTRI IS DESIGNED TO BE CALLED REPEATEDLY WITHOUT WASTING
C       EARLIER WORK.  THE PARAMETER NCALLS IS USED TO INDICATE TO
C       CUBTRI AT WHAT POINT TO RESTART, AND MUST BE RE-INITIALIZED
C       TO 0 WHEN A NEW INTEGRAL IS TO BE COMPUTED.  AT LEAST ONE OF
C       THE PARAMETERS EPS, MCALLS AND NW MUST BE CHANGED BETWEEN
C       CALLS TO CUBTRI, ACCORDING TO THE RETURNED VALUE OF IER. NONE
C       OF THE OTHER PARAMETERS MAY BE CHANGED IF RESTARTING IS DONE.
C       IF IER=3 IS ENCOUNTERED, THERE PROBABLY IS A SINGULARITY
C       SOMEWHERE IN THE REGION.  THE ERROR MESSAGE INDICATES THAT
C       FURTHER SUBDIVISION IS IMPOSSIBLE BECAUSE THE VERTICES OF THE
C       SMALLER TRIANGLES PRODUCED WILL BEGIN TO COALESCE TO THE
C       PRECISION OF THE COMPUTER.  THIS SITUATION CAN USUALLY BE
C       RELIEVED BY SPECIFYING THE REGION IN SUCH A WAY THAT THE
C       SINGULARITY IS LOCATED AT THE THIRD VERTEX OF THE TRIANGLE.
C       IF IER=4 IS ENCOUNTERED, THE VALUE OF THE INTEGRAL CANNOT BE
C       IMPROVED ANY FURTHER. THE ONLY EXCEPTION TO THIS OCCURS WHEN A
C       FUNCTION WITH HIGHLY IRREGULAR BEHAVIOUR IS INTEGRATED (E.G.
C       FUNCTIONS WITH JUMP DISCONTINUITIES OR VERY HIGHLY OSCILLATORY
C       FUNCTIONS). IN SUCH A CASE THE USER CAN DISABLE THE ROUNDING
C       ERROR TEST BY REMOVING THE IF STATEMENT SHORTLY AFTER STATEMENT
C       NUMBER 70.
C
      EXTERNAL F
      INTEGER IDATA(1), IER, MCALLS, NCALLS, NW
      REAL ALFA, ANS, ANSKP, AREA, EPS, ERR, ERRMAX, H, Q1, Q2, R1, R2,
     * RDATA(1), D(2,4), S(4), T(2,3), VEC(2,3), W(6,NW), X(2)
C       ACTUAL DIMENSION OF W IS (6,NW/6)
C
      DOUBLE PRECISION TANS, TERR, DZERO
      COMMON /CUBSTA/ TANS, TERR
C       THIS COMMON IS REQUIRED TO PRESERVE TANS AND TERR BETWEEN CALLS
C       AND TO SAVE VARIABLES IN FUNCTION RNDERR
      DATA NFE /19/, S(1), S(2), S(3), S(4) /3*1E0,-1E0/, D(1,1),
     * D(2,1) /0.0,0.0/, D(1,2), D(2,2) /0.0,1.0/, D(1,3), D(2,3)
     * /1.0,0.0/, D(1,4), D(2,4) /1.0,1.0/
C       NFE IS THE NUMBER OF FUNCTION EVALUATIONS PER CALL TO CUBRUL.
      DATA ZERO /0.E0/, ONE /1.E0/, DZERO /0.D0/, POINT5 /.5E0/
C
C      CALCULATE DIRECTION VECTORS, AREA AND MAXIMUM NUMBER
C      OF SUBDIVISIONS THAT MAY BE PERFORMED
      DO 20 I=1,2
        VEC(I,3) = T(I,3)
        DO 10 J=1,2
          VEC(I,J) = T(I,J) - T(I,3)
   10   CONTINUE
   20 CONTINUE
      MAXC = (MCALLS/NFE+3)/4
      IER = 1
      MAXK = MIN0(MAXC,(NW/6+2)/3)
      IF (MAXC.GT.MAXK) IER = 2
      AREA = ABS(VEC(1,1)*VEC(2,2)-VEC(1,2)*VEC(2,1))*POINT5
      K = (NCALLS/NFE+3)/4
      MW = 3*(K-1) + 1
      IF (NCALLS.GT.0) GO TO 30
C
C       TEST FOR TRIVIAL CASES
      TANS = DZERO
      TERR = DZERO
      IF (AREA.EQ.ZERO) GO TO 90
      IF (MCALLS.LT.NFE) GO TO 100
      IF (NW.LT.6) GO TO 110
C
C       INITIALIZE DATA LIST
      K = 1
      MW = 1
      W(1,1) = ZERO
      W(2,1) = ZERO
      W(3,1) = ONE
      CALL CUBRUL(F, VEC, W(1,1), IDATA, RDATA)
      TANS = W(5,1)
      TERR = W(6,1)
      NCALLS = NFE
C
C       TEST TERMINATION CRITERIA
   30 ANS = TANS
      ERR = TERR
      IF (ERR.LT.AMAX1(ONE,ABS(ANS))*EPS) GO TO 90
      IF (K.EQ.MAXK) GO TO 120
C
C       FIND TRIANGLE WITH LARGEST ERROR
      ERRMAX = ZERO
      DO 40 I=1,MW
        IF (W(6,I).LE.ERRMAX) GO TO 40
        ERRMAX = W(6,I)
        J = I
   40 CONTINUE
C
C       SUBDIVIDE TRIANGLE INTO FOUR SUBTRIANGLES AND UPDATE DATA LIST
      DO 50 I=1,2
        X(I) = W(I,J)
   50 CONTINUE
      H = W(3,J)*POINT5
      IF (RNDERR(X(1),H,X(1),H).NE.ZERO) GO TO 130
      IF (RNDERR(X(2),H,X(2),H).NE.ZERO) GO TO 130
      ANSKP = SNGL(TANS)
      TANS = TANS - DBLE(W(5,J))
      TERR = TERR - DBLE(W(6,J))
      R1 = W(4,J)
      R2 = W(5,J)
      JKP = J
      Q1 = ZERO
      Q2 = ZERO
      DO 70 I=1,4
        DO 60 L=1,2
          W(L,J) = X(L) + H*D(L,I)
   60   CONTINUE
        W(3,J) = H*S(I)
        CALL CUBRUL(F, VEC, W(1,J), IDATA, RDATA)
        Q2 = Q2 + W(5,J)
        Q1 = Q1 + W(4,J)
        J = MW + I
   70 CONTINUE
      ALFA = 1E15
      IF (Q2.NE.R2) ALFA = ABS((Q1-R1)/(Q2-R2)-ONE)
      J = JKP
      DO 80 I=1,4
        W(6,J) = W(6,J)/ALFA
        TANS = TANS + W(5,J)
        TERR = TERR + W(6,J)
        J = MW + I
   80 CONTINUE
      MW = MW + 3
      NCALLS = NCALLS + 4*NFE
      K = K + 1
C
C       IF ANSWER IS UNCHANGED, IT CANNOT BE IMPROVED
      IF (ANSKP.EQ.SNGL(TANS)) GO TO 150
C
C       REMOVE THIS IF STATEMENT TO DISABLE ROUNDING ERROR TEST
      IF (K.GT.3 .AND. ABS(Q2-R2).GT.ABS(Q1-R1)) GO TO 140
      GO TO 30
C
C       EXITS FROM SUBROUTINE
   90 IER = 0
      GO TO 120
  100 IER = 1
      GO TO 120
  110 IER = 2
  120 ANS = TANS
      ERR = TERR
      RETURN
  130 IER = 3
      GO TO 120
  140 IER = 4
      GO TO 120
  150 IER = 5
      GO TO 120
      END
      FUNCTION RNDERR(X, A, Y, B)                                       RND   10
C       THIS FUNCTION COMPUTES THE ROUNDING ERROR COMMITTED WHEN THE
C       SUM X+A IS FORMED.  IN THE CALLING PROGRAM, Y MUST BE THE SAME
C       AS X AND B MUST BE THE SAME AS A.  THEY ARE DECLARED AS
C       DISTINCT VARIABLES IN THIS FUNCTION, AND THE INTERMEDIATE
C       VARIABLES S AND T ARE PUT INTO COMMON, IN ORDER TO DEFEND
C       AGAINST THE WELL-MEANING ACTIONS OF SOME OFFICIOUS OPTIMIZING
C       FORTRAN COMPILERS.
      COMMON /CUBATB/ S, T
      S = X + A
      T = S - Y
      RNDERR = T - B
      RETURN
      END
      SUBROUTINE CUBRUL(F, VEC, P, IDATA, RDATA)                        CUB   10
C
C       BASIC CUBATURE RULE PAIR OVER A TRIANGLE
C
C       PARAMETERS
C         F  - EXTERNAL FUNCTION - SEE COMMENTS TO CUBTRI
C         VEC- MATRIX OF BASE VECTORS AND ORIGIN (INPUT)
C         P  - TRIANGLE DESCRIPTION VECTOR OF DIMENSION 6
C               P(1) - TRANSFORMED X COORDINATE OF ORIGIN VERTEX(INPUT)
C               P(2) - TRANSFORMED Y COORDINATE OF ORIGIN VERTEX(INPUT)
C               P(3) - DISTANCE OF OTHER VERTICES IN THE DIRECTIONS
C                     OF THE BASE VECTORS (INPUT)
C               P(4) - LESS ACCURATE ESTIMATED INTEGRAL (OUTPUT)
C               P(5) - MORE ACCURATE ESTIMATED INTEGRAL (OUTPUT)
C               P(6) - ABS(P(5)-P(4))   (OUTPUT)
C
C       CUBRUL EVALUATES A LINEAR COMBINATION OF BASIC INTEGRATION
C       RULES HAVING D3 SYMMETRY.  THE AREAL COORDINATES PERTAINING TO
C       THE J-TH RULE ARE STORED IN W(I,J),I=1,2,3.  THE CORRESPONDING
C       WEIGHTS ARE W(4,J) AND W(5,J), WITH W(5,J) BELONGING TO THE
C       MORE ACCURATE FORMULA.  IF W(1,J).EQ.W(2,J), THE INTEGRATION
C       POINT IS THE CENTROID, ELSE IF W(2,J).EQ.W(3,J), THE EVALUATION
C       POINTS ARE ON THE MEDIANS.  IN BOTH CASES ADVANTAGE IS TAKEN OF
C       SYMMETRY TO AVOID REPEATING FUNCTION EVALUATIONS.
C
C       THE FOLLOWING DOUBLE PRECISION VARIABLES ARE USED TO AVOID
C       UNNECESSARY ROUNDING ERRORS IN FLOATING POINT ADDITION.
C       THEY MAY BE DECLARED SINGLE PRECISION IF DOUBLE PRECISION IS
C       NOT AVAILABLE AND FULL ACCURACY IS NOT NEEDED.
      DOUBLE PRECISION A1, A2, S, SN, DZERO, DONE, DTHREE, DSIX
      REAL AREA, ORIGIN(2), P(6), RDATA(1), TVEC(2,3), VEC(2,3), W(5,6)
      INTEGER IDATA(1)
C
C       W CONTAINS POINTS AND WEIGHTS OF THE INTEGRATION FORMULAE
C       NQUAD - NUMBER OF BASIC RULES USED
C
C       THIS PARTICULAR RULE IS THE 19 POINT EXTENSION (DEGREE 8) OF
C       THE FAMILIAR 7 POINT RULE (DEGREE 5).
C
C     SIGMA=SQRT(7)
C     PHI=SQRT(15)
C     W(1,1),W(2,1),W(3,1) = 1/3
C     W(4,1) = 9/40
C     W(5,1) = 7137/62720 - 45*SIGMA/1568
C     W(1,2) = 3/7 + 2*PHI/21
C     W(2,2),W(3,2) = 2/7 - PHI/21
C     W(4,2) = 31/80 - PHI/400
C     W(5,2) = - 9301697/4695040 - 13517313*PHI/23475200
C            + 764885*SIGMA/939008 + 198763*PHI*SIGMA/939008
C     W(*,3) = W(*,2) WITH PHI REPLACED BY -PHI
C     W(1,5) = 4/9 + PHI/9 + SIGMA/9 - SIGMA*PHI/45
C     W(2,5),W(3,5) = 5/18 - PHI/18 - SIGMA/18 + SIGMA*PHI/90
C     W(4,5) = 0
C     W(5,5) = 102791225/59157504 + 23876225*PHI/59157504
C            - 34500875*SIGMA/59157504 - 9914825*PHI*SIGMA/59157504
C     W(*,4) = W(*,5) WITH PHI REPLACED BY -PHI
C     W(1,6) = 4/9 + SIGMA/9
C     W(2,6) = W(2,4)
C     W(3,6) = W(2,5)
C     W(4,6) = 0
C     W(5,6) = 11075/8064 - 125*SIGMA/288
C
      DATA NQUAD /6/, W(1,1), W(2,1), W(3,1) /3*.33333333333333333333333
     * 33E0/, W(4,1), W(5,1) /.225E0,.3786109120031468330830822E-1/,
     * W(1,2), W(2,2), W(3,2) /.7974269853530873223980253E0,2*
     * .1012865073234563388009874E0/, W(4,2), W(5,2)
     * /.3778175416344814577870518E0,.1128612762395489164329420E0/,
     * W(1,3), W(2,3), W(3,3) /.5971587178976982045911758E-1,2*
     * .4701420641051150897704412E0/, W(4,3), W(5,3)
     * /.3971824583655185422129482E0,.2350720567323520126663380E0/
      DATA W(1,4), W(2,4), W(3,4) /.5357953464498992646629509E0,2*
     * .2321023267750503676685246E0/, W(4,4), W(5,4)
     * /0.E0,.3488144389708976891842461E0/, W(1,5), W(2,5), W(3,5)
     * /.9410382782311208665596304E0,2*.2948086088443956672018481E-1/,
     * W(4,5), W(5,5) /0.E0,.4033280212549620569433320E-1/, W(1,6),
     * W(2,6), W(3,6) /.7384168123405100656112906E0,
     * .2321023267750503676685246E0,.2948086088443956672018481E-1/,
     * W(4,6), W(5,6) /0.E0,.2250583347313904927138324E0/
C
      DATA DZERO /0.D0/, DONE /1.D0/, DTHREE /3.D0/, DSIX /6.D0/,
     * POINT5 /.5E0/
C
C       SCALE BASE VECTORS AND OBTAIN AREA
      DO 20 I=1,2
        ORIGIN(I) = VEC(I,3) + P(1)*VEC(I,1) + P(2)*VEC(I,2)
        DO 10 J=1,2
          TVEC(I,J) = P(3)*VEC(I,J)
   10   CONTINUE
   20 CONTINUE
      AREA = POINT5*ABS(TVEC(1,1)*TVEC(2,2)-TVEC(1,2)*TVEC(2,1))
      A1 = DZERO
      A2 = DZERO
C
C       COMPUTE ESTIMATES FOR INTEGRAL AND ERROR
      DO 40 K=1,NQUAD
        X = ORIGIN(1) + W(1,K)*TVEC(1,1) + W(2,K)*TVEC(1,2)
        Y = ORIGIN(2) + W(1,K)*TVEC(2,1) + W(2,K)*TVEC(2,2)
        S = DBLE(F(X,Y,IDATA,RDATA))
        SN = DONE
        IF (W(1,K).EQ.W(2,K)) GO TO 30
        X = ORIGIN(1) + W(2,K)*TVEC(1,1) + W(1,K)*TVEC(1,2)
        Y = ORIGIN(2) + W(2,K)*TVEC(2,1) + W(1,K)*TVEC(2,2)
        S = S + DBLE(F(X,Y,IDATA,RDATA))
        X = ORIGIN(1) + W(2,K)*TVEC(1,1) + W(3,K)*TVEC(1,2)
        Y = ORIGIN(2) + W(2,K)*TVEC(2,1) + W(3,K)*TVEC(2,2)
        S = S + DBLE(F(X,Y,IDATA,RDATA))
        SN = DTHREE
        IF (W(2,K).EQ.W(3,K)) GO TO 30
        X = ORIGIN(1) + W(1,K)*TVEC(1,1) + W(3,K)*TVEC(1,2)
        Y = ORIGIN(2) + W(1,K)*TVEC(2,1) + W(3,K)*TVEC(2,2)
        S = S + DBLE(F(X,Y,IDATA,RDATA))
        X = ORIGIN(1) + W(3,K)*TVEC(1,1) + W(1,K)*TVEC(1,2)
        Y = ORIGIN(2) + W(3,K)*TVEC(2,1) + W(1,K)*TVEC(2,2)
        S = S + DBLE(F(X,Y,IDATA,RDATA))
        X = ORIGIN(1) + W(3,K)*TVEC(1,1) + W(2,K)*TVEC(1,2)
        Y = ORIGIN(2) + W(3,K)*TVEC(2,1) + W(2,K)*TVEC(2,2)
        S = S + DBLE(F(X,Y,IDATA,RDATA))
        SN = DSIX
   30   S = S/SN
        A1 = A1 + W(4,K)*S
        A2 = A2 + W(5,K)*S
   40 CONTINUE
      P(4) = SNGL(A1)*AREA
      P(5) = SNGL(A2)*AREA
      P(6) = ABS(P(5)-P(4))
      RETURN
      END
