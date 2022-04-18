      SUBROUTINE C05AZF(X,Y,FX,TOLX,IR,C,IND,IFAIL)
C     MARK 8 RE-ISSUE. NAG COPYRIGHT 1979.
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     MARK 12A REVISED. IER-496 (AUG 1986).
C     MARK 13 REVISED. USE OF MARK 12 X02 FUNCTIONS (APR 1988).
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='C05AZF')
C     .. Scalar Arguments ..
      DOUBLE PRECISION  FX, TOLX, X, Y
      INTEGER           IFAIL, IND, IR
C     .. Array Arguments ..
      DOUBLE PRECISION  C(17)
C     .. Local Scalars ..
      DOUBLE PRECISION  AB, DIFF, DIFF1, DIFF2, REL, RMAX, TOL, TOL1
      INTEGER           I
      LOGICAL           T
C     .. Local Arrays ..
      CHARACTER*1       P01REC(1)
C     .. External Functions ..
      DOUBLE PRECISION  X02AJF, X02AKF
      INTEGER           P01ABF
      EXTERNAL          X02AJF, X02AKF, P01ABF
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, MAX, SIGN, DBLE, SQRT, INT
C     .. Executable Statements ..
      I = 0
      IF ((IND.GT.0 .AND. IND.LE.4) .OR. IND.EQ.-1) GO TO 20
C     USER NOT CHECKED IND OR CHANGED IT
      I = 2
      IND = 0
      GO TO 640
   20 IF (TOLX.GT.0.D0 .AND. (IR.EQ.0 .OR. IR.EQ.1 .OR. IR.EQ.2))
     *    GO TO 40
      I = 3
      IND = 0
      GO TO 640
   40 REL = 1.D0
      AB = 1.D0
      IF (IR.EQ.1) REL = 0.D0
      IF (IR.EQ.2) AB = 0.D0
      IF (IND.EQ.-1) GO TO 80
      GO TO (60,100,180,480) IND
   60 C(3) = X
      IND = 2
      RETURN
   80 C(3) = X
  100 IF (FX.NE.0.D0) GO TO 140
  120 Y = X
      IND = 0
      I = 0
      GO TO 640
  140 C(4) = FX
      C(15) = ABS(FX)
      C(16) = 0.D0
      X = Y
      Y = C(3)
      C(2) = C(4)
      C(5) = X
      IF (IND.EQ.-1) GO TO 160
      IND = 3
      RETURN
  160 FX = C(1)
      IND = 3
  180 IF (FX.EQ.0.D0) GO TO 120
      IF (SIGN(1.D0,FX).NE.SIGN(1.D0,C(2))) GO TO 200
      IND = 0
      I = 1
      GO TO 640
  200 C(6) = FX
      C(13) = SQRT(X02AJF())
      C(15) = MAX(C(15),ABS(FX))
      C(14) = X02AKF()
      C(16) = 0.0D0
  220 C(1) = C(5)
      C(2) = C(6)
      C(17) = 0.D0
  240 IF (ABS(C(2)).GE.ABS(C(4))) GO TO 280
      IF (C(1).EQ.C(5)) GO TO 260
      C(7) = C(5)
      C(8) = C(6)
  260 C(5) = C(3)
      C(6) = C(4)
      X = C(1)
      C(3) = X
      C(4) = C(2)
      C(1) = C(5)
      C(2) = C(6)
  280 TOL = 0.5D0*TOLX*MAX(AB,REL*ABS(C(3)))
      TOL1 = 2.0D0*X02AJF()*MAX(AB,REL*ABS(C(3)))
      DIFF2 = 0.5D0*(C(1)-C(3))
      C(12) = DIFF2
      DIFF2 = DIFF2 + C(3)
      IF (C(12).EQ.0.D0) GO TO 340
      IF (ABS(C(12)).LE.TOL) GO TO 580
      IF (ABS(C(12)).LE.TOL1) GO TO 340
      IF (C(17).LT.2.5D0) GO TO 300
      C(11) = C(12)
      GO TO 460
  300 TOL = TOL*SIGN(1.D0,C(12))
      DIFF1 = (C(3)-C(5))*C(4)
      IF (C(17).GT.1.5D0) GO TO 320
      DIFF = C(6) - C(4)
      GO TO 380
  320 IF (C(7).NE.C(3) .AND. C(7).NE.C(5)) GO TO 360
  340 IND = 0
      I = 5
      GO TO 640
  360 C(9) = (C(8)-C(4))/(C(7)-C(3))
      C(10) = (C(8)-C(6))/(C(7)-C(5))
      DIFF1 = C(10)*DIFF1
      DIFF = C(9)*C(6) - C(10)*C(4)
  380 IF (DIFF1.GE.0.D0) GO TO 400
      DIFF1 = -DIFF1
      DIFF = -DIFF
  400 IF (ABS(DIFF1).GT.C(14) .AND. DIFF1.GT.DIFF*TOL) GO TO 420
      C(11) = TOL
      GO TO 460
  420 IF (DIFF1.GE.C(12)*DIFF) GO TO 440
      C(11) = DIFF1/DIFF
      GO TO 460
  440 C(11) = C(12)
  460 C(7) = C(5)
      C(8) = C(6)
      C(5) = C(3)
      C(6) = C(4)
      C(3) = C(3) + C(11)
      X = C(3)
      Y = C(1)
      IND = 4
      RETURN
  480 IF (FX.EQ.0.D0) GO TO 120
      C(4) = FX
      RMAX = ABS(FX)
      IF (C(13)*RMAX.LE.C(15)) GO TO 500
      IF (C(16).EQ.1.D0) C(16) = -1.D0
      IF (C(16).EQ.0.D0) C(16) = 1.D0
      GO TO 520
  500 C(16) = 0.D0
  520 IF (C(2).GE.0.D0) GO TO 540
      T = C(4) .LE. 0.D0
      GO TO 560
  540 T = C(4) .GE. 0.D0
  560 IF (T) GO TO 220
      I = INT(C(17)+0.1D0)
      I = I + 1
      IF (C(11).EQ.C(12)) I = 0
      C(17) = DBLE(I)
      GO TO 240
  580 IF (C(16).GE.0.D0) GO TO 600
      I = 4
      GO TO 620
  600 Y = C(1)
      I = 0
  620 IND = 0
  640 IFAIL = P01ABF(IFAIL,I,SRNAME,0,P01REC)
      RETURN
      END
      SUBROUTINE D02KAY(NIT,IFLAG,ELAM,FINFO)
C     MARK 11.5 RELEASE. NAG COPYRIGHT 1986.
C     DUMMY MONIT ROUTINE FOR D02KAF, D02KDF OR D02KEF
C     .. Scalar Arguments ..
      DOUBLE PRECISION  ELAM
      INTEGER           IFLAG, NIT
C     .. Array Arguments ..
      DOUBLE PRECISION  FINFO(15)
C     .. Executable Statements ..
      RETURN
      END
      DOUBLE PRECISION FUNCTION D02KDS(X)
C     MARK 7 RELEASE. NAG COPYRIGHT 1978.
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     MARK 13 REVISED. USE OF MARK 12 X02 FUNCTIONS (APR 1988).
C     EXP AVOIDING UNDERFLOW ERROR
C     .. Scalar Arguments ..
      DOUBLE PRECISION                 X
C     .. External Functions ..
      DOUBLE PRECISION                 X02AMF
      EXTERNAL                         X02AMF
C     .. Intrinsic Functions ..
      INTRINSIC                        EXP, LOG
C     .. Executable Statements ..
      D02KDS = 0.D0
      IF (X.GE.LOG(X02AMF())) D02KDS = EXP(X)
      RETURN
      END
      SUBROUTINE D02KDT(V,Y,PYP,K,IFAIL)
C     MARK 7 RELEASE. NAG COPYRIGHT 1978.
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     .. Scalar Arguments ..
      DOUBLE PRECISION  PYP, Y
      INTEGER           IFAIL, K
C     .. Array Arguments ..
      DOUBLE PRECISION  V(3)
C     .. Scalars in Common ..
      DOUBLE PRECISION  BP, LAMDA, MINSC, ONE, PI, PSIGN, TWO, ZER
      INTEGER           JINT
C     .. Arrays in Common ..
      DOUBLE PRECISION  YL(3), YR(3)
C     .. Local Scalars ..
      DOUBLE PRECISION  B, PHI, R2
C     .. Intrinsic Functions ..
      INTRINSIC         LOG, ATAN2, DBLE
C     .. Common blocks ..
      COMMON            /AD02KD/ZER, ONE, TWO, PI, LAMDA, PSIGN, MINSC,
     *                  BP, YL, YR, JINT
C     .. Executable Statements ..
      B = V(1)
      IF (B.LE.ZER) GO TO 100
      R2 = Y*Y*B + PYP*PYP/B
      IF (R2.EQ.ZER) GO TO 80
      V(3) = LOG(R2)
      PHI = ATAN2(Y*B,PYP)
      IF (K.GE.0) GO TO 20
C     INITIAL BOUNDARY CONDITION
      IF (PHI.LT.ZER) PHI = PHI + PI
      IF (PHI.GE.PI) PHI = PHI - PI
      GO TO 40
C     FINAL BOUNDARY CONDITION
   20 IF (PHI.LE.ZER) PHI = PHI + PI
      IF (PHI.GT.PI) PHI = PHI - PI
      PHI = PHI + DBLE(K)*PI
   40 V(2) = TWO*PHI
      IFAIL = 0
   60 RETURN
   80 IFAIL = 1
      GO TO 60
  100 IFAIL = 2
      GO TO 60
      END
      DOUBLE PRECISION FUNCTION D02KDU(X,COEFFN)
C     MARK 7 RELEASE. NAG COPYRIGHT 1978.
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     COEFFN
C     .. Scalar Arguments ..
      DOUBLE PRECISION                 X
C     .. Subroutine Arguments ..
      EXTERNAL                         COEFFN
C     .. Scalars in Common ..
      DOUBLE PRECISION                 BP, LAMDA, MINSC, ONE, PI, PSIGN,
     *                                 TWO, ZER
      INTEGER                          JINT
C     .. Arrays in Common ..
      DOUBLE PRECISION                 YL(3), YR(3)
C     .. Local Scalars ..
      DOUBLE PRECISION                 DQDL, P, Q
C     .. Intrinsic Functions ..
      INTRINSIC                        SQRT
C     .. Common blocks ..
      COMMON                           /AD02KD/ZER, ONE, TWO, PI, LAMDA,
     *                                 PSIGN, MINSC, BP, YL, YR, JINT
C     .. Executable Statements ..
      CALL COEFFN(P,Q,DQDL,X,LAMDA,JINT)
      Q = P*Q
      P = MINSC*P
      P = P*P
      D02KDU = SQRT(SQRT(P*P+Q*Q))
      RETURN
      END
      SUBROUTINE D02KDV(Y,BNEW,IFAIL)
C     MARK 7 RELEASE. NAG COPYRIGHT 1978.
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     MARK 12A REVISED. IER-499 (AUG 1986).
C     .. Scalar Arguments ..
      DOUBLE PRECISION  BNEW
      INTEGER           IFAIL
C     .. Array Arguments ..
      DOUBLE PRECISION  Y(3)
C     .. Scalars in Common ..
      DOUBLE PRECISION  BP, LAMDA, MINSC, ONE, PI, PSIGN, TWO, ZER
      INTEGER           JINT
C     .. Arrays in Common ..
      DOUBLE PRECISION  YL(3), YR(3)
C     .. Local Scalars ..
      DOUBLE PRECISION  B, BUP, CPHI, SPHI
C     .. Intrinsic Functions ..
      INTRINSIC         LOG, ATAN2, COS, SIN
C     .. Common blocks ..
      COMMON            /AD02KD/ZER, ONE, TWO, PI, LAMDA, PSIGN, MINSC,
     *                  BP, YL, YR, JINT
C     .. Executable Statements ..
      B = Y(1)
      IF (B.LE.ZER .OR. BNEW.LE.ZER) GO TO 20
      Y(1) = BNEW
      BUP = B/BNEW
      B = BNEW/B
      CPHI = COS(Y(2))
      SPHI = SIN(Y(2))
      Y(2) = Y(2) + TWO*ATAN2((B-ONE)*SPHI,B+ONE-(B-ONE)*CPHI)
      Y(3) = Y(3) + LOG((B+BUP-(B-BUP)*CPHI)/TWO)
      IFAIL = 0
      RETURN
   20 IFAIL = 1
      RETURN
      END
      SUBROUTINE D02KDW(N,X,V,F,COEFFN,COEFF1,M,ARR)
C     MARK 7 RELEASE. NAG COPYRIGHT 1978.
C     MARK 8 REVISED. IER-227 (APR 1980).
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     POLYGONAL SCALING METHOD
C     SEVERAL FORMULAE FOR P,Q OVER RANGE,SELECTED BY JINT
C     COEFF1, COEFFN
C     .. Scalar Arguments ..
      DOUBLE PRECISION  X
      INTEGER           M, N
C     .. Array Arguments ..
      DOUBLE PRECISION  ARR(M), F(N), V(N)
C     .. Subroutine Arguments ..
      EXTERNAL          COEFF1, COEFFN
C     .. Scalars in Common ..
      DOUBLE PRECISION  BP, LAMDA, MINSC, ONE, PI, PSIGN, TWO, ZER
      INTEGER           JINT
C     .. Arrays in Common ..
      DOUBLE PRECISION  YL(3), YR(3)
C     .. Local Scalars ..
      DOUBLE PRECISION  B, C, DQDL, P, Q, S, T1, T2
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, COS, SIN
C     .. Common blocks ..
      COMMON            /AD02KD/ZER, ONE, TWO, PI, LAMDA, PSIGN, MINSC,
     *                  BP, YL, YR, JINT
C     .. Executable Statements ..
      CALL COEFFN(P,Q,DQDL,X,LAMDA,JINT)
C     TEST IF P(X) WAS 0 OR CHANGED SIGN EARLIER (PSIGN=0)
C     OR AT THIS CALL (PSIGN .NE. 0 BUT P*PSIGN .LE. 0)
      IF (P*PSIGN.LE.ZER) GO TO 20
      IF (P.LT.ZER) Q = -Q
      P = ABS(P)
      B = V(1)
      T1 = B/P - Q/B
      T2 = BP/B
      C = COS(V(2))
      S = SIN(V(2))
      F(1) = BP
      F(2) = B/P + Q/B + T1*C + T2*S
      F(3) = -T2*C + T1*S
      RETURN
   20 F(1) = ZER
      F(2) = ZER
      F(3) = ZER
      PSIGN = ZER
      RETURN
      END
      SUBROUTINE D02KDY(X,XEND,N,Y,CIN,TOL,FCN,COMM,CONST,COUT,W,IW,IW1,
     *                  COEFFN,COEFF1,ARR,M,IFAIL)
C     MARK 7 RELEASE. NAG COPYRIGHT 1978.
C     MARK 7F REVISED. IER-209 (OCT 1979)
C     MARK 7G REVISED. IER-216 (FEB 1980)
C     MARK 8 REVISED. IER-227 (APR 1980), IER-246 (MAY 1980).
C     MARK 8A REVISED. IER-251 (AUG 1980).
C     MARK 11 REVISED. IER-419 (FEB 1984).
C     MARK 11D REVISED. IER-469 (NOV.1985).
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     MARK 13 REVISED. USE OF MARK 12 X02 FUNCTIONS (APR 1988).
C     MARK 13 REVISED. IER-633 (APR 1988).
C     MARK 14 REVISED. IER-715 (DEC 1989).
C     MARK 14B REVISED. IER-837 (MAR 1990).
C     INTEGRATES THE N DIFFERENTIAL EQUATIONS DEFINED
C     BY FCN FROM X TO XEND.  THE VALUES AT X
C     OF THE SOLUTION MUST BE GIVEN IN Y
C     AND THE CALCULATED VALUES ARE RETURNED
C     IN THE SAME VECTOR.  THE LOCAL ERROR
C     PER STEP IS CONTROLLED BY TOL. VARIOUS
C     OPTIONS ARE CONTROLLED BY CIN AND
C     THE WORKSPACE W WHICH MUST HAVE
C     FIRST DIMENSION N1 .GE. N AND
C     SECOND DIMENSION .GE.7.  USEFUL
C     OUTPUT IS ALSO RETURNED IN CIN AND
C     W, PERMITTING EFFICIENT INTEGRATION.
C
C     COEFF1, COEFFN, FCN
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='D02KDY')
C     .. Scalar Arguments ..
      DOUBLE PRECISION  TOL, X, XEND
      INTEGER           IFAIL, IW, IW1, M, N
C     .. Array Arguments ..
      DOUBLE PRECISION  ARR(M), CIN(6), COMM(5), CONST(3), COUT(14),
     *                  W(IW,IW1), Y(N)
C     .. Subroutine Arguments ..
      EXTERNAL          COEFF1, COEFFN, FCN
C     .. Local Scalars ..
      DOUBLE PRECISION  COUT12, EXP, EXPI, FAC, FAC1, FAC2, FAC3, FAC4,
     *                  FAC5, FAC6, FAC7, FAC8, FAC9, FNORM, HEST,
     *                  HEST1, RAT, RAT1, S, SMALL, T, TOLEST, XORIG,
     *                  YNORM
      INTEGER           I, I1, IEND, IND, ISIG, ISTEP, J, K, N2
      LOGICAL           CALLS, INIT, INTER, REDUCE, START, TEST
C     .. Local Arrays ..
      CHARACTER*1       P01REC(1)
C     .. External Functions ..
      DOUBLE PRECISION  D02PAY, X02AJF, X02AKF
      INTEGER           P01ABF
      EXTERNAL          D02PAY, X02AJF, X02AKF, P01ABF
C     .. External Subroutines ..
      EXTERNAL          D02KDZ
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, MAX, MIN, SIGN, DBLE, SQRT, INT
C     .. Save statement ..
      SAVE              START
C     .. Data statements ..
      DATA              FAC/10.D0/, FAC1/0.1D0/, FAC3/2.D0/,
     *                  FAC4/0.8D0/, FAC5/0.25D0/, FAC6/0.25D0/,
     *                  FAC8/0.2D0/, FAC9/0.2D0/
C     .. Executable Statements ..
C
C     SET CONSTANTS
C
      XORIG = X
      HEST1 = 0.D0
      INIT = .FALSE.
      EXP = FAC6
      EXPI = FAC5
      IND = 0
      ISIG = 0
      INTER = .FALSE.
      ISTEP = 2
      N2 = 7
      TEST = .FALSE.
      CALLS = .FALSE.
      COUT12 = 0.D0
      IF (CIN(1).NE.7.D0) GO TO 20
      COUT12 = COUT(12)
      CIN(1) = 1.D0
      GO TO 40
   20 COUT(12) = 0.D0
   40 IF (CIN(1).EQ.0.D0) GO TO 100
      DO 60 I = 1, 3
         IF (COMM(I).NE.0.D0) GO TO 100
   60 CONTINUE
      IF (CIN(1).NE.5.D0 .OR. COMM(4).LE.0.D0) GO TO 80
      TEST = .TRUE.
      IF (CIN(2).GE.3.D0) GO TO 360
      GO TO 860
   80 IF (CIN(1).NE.6.D0 .OR. COMM(4).GE.0.D0) GO TO 100
      TEST = .TRUE.
      IF (CIN(2).GE.3.D0) GO TO 360
      GO TO 720
C
C     ARGUMENT TESTS
C
  100 IF (IW.GE.N .AND. N.GT.0 .AND. IW1.GE.7) GO TO 140
C     WORKSPACE WRONGLY DIMENSIONED
  120 IND = 1
      CIN(1) = -1.D0
      GO TO 2620
  140 IF (TOL.GT.0.D0) GO TO 160
C     ERROR TOLERANCE NOT POSITIVE
      IND = 1
      CIN(1) = -2.D0
      GO TO 2620
  160 IF (CIN(1).EQ.0.D0) GO TO 180
      IF (CIN(1).EQ.1.D0) GO TO 240
      IF (CIN(1).GE.2.D0 .AND. CIN(1).LE.6.D0) GO TO 280
C     CIN(1) OUT OF RANGE
      IND = 1
      CIN(1) = -3.D0
      GO TO 2620
C
C     SET INPUT AND OUTPUT PARAMETERS
C
  180 DO 200 I = 2, 5
         CIN(I) = 0.D0
  200 CONTINUE
      DO 220 I = 1, 4
         COMM(I) = 0.D0
  220 CONTINUE
      CONST(1) = 0.D0
      CONST(2) = FAC3
      CONST(3) = FAC4
  240 COUT(11) = X02AJF()
      SMALL = X02AKF()
      COUT(13) = ABS(CIN(3))
      COUT(14) = ABS(CIN(4))
      CIN(6) = 0.D0
      ISTEP = 0
      S = 0.D0
      DO 260 I = 1, N
         S = MAX(S,ABS(Y(I)))
  260 CONTINUE
      INIT = .TRUE.
      COUT(1) = 0.D0
      COUT(2) = 0.D0
      COUT(3) = 0.D0
      COUT(4) = X
      COUT(5) = X
      COUT(6) = S
      COUT(7) = S
      COUT(8) = 0.D0
      COUT(9) = 0.D0
      COUT(10) = 0.D0
      IF (CIN(1).EQ.0.D0) GO TO 800
      GO TO 320
C
C     ARGUMENT TESTS FOR CIN(1).GT.0.
C
  280 IF ((COUT(4).GE.COUT(5) .AND. X.LT.COUT(4)) .OR. (COUT(4)
     *    .LT.COUT(5) .AND. X.GE.COUT(4))) GO TO 300
      IF (X.EQ.XEND) GO TO 320
      IF ((X.GE.COUT(4) .AND. XEND.GE.X) .OR. (X.LT.COUT(4)
     *    .AND. XEND.LT.X)) GO TO 320
  300 CIN(1) = 1.D0
      INTER = .TRUE.
      GO TO 240
  320 IF (CIN(2).EQ.0.D0 .OR. CIN(2).EQ.1.D0 .OR. CIN(2)
     *    .EQ.2.D0 .OR. CIN(2).EQ.3.D0 .OR. CIN(2).EQ.4.D0) GO TO 340
C     CIN(2) OUT OF RANGE
      IND = 1
      CIN(1) = -4.D0
      GO TO 2620
  340 IF (CIN(2).LT.3.D0) GO TO 420
  360 IF (IW1.LT.8) GO TO 120
      N2 = 8
      DO 380 I = 1, N
         IF (W(I,7).GT.0.D0) GO TO 380
C        SCALING OF TOLERANCE NON-POSITIVE
         IND = 1
         CIN(1) = -5.D0
         GO TO 2620
  380 CONTINUE
      IF (CIN(2).EQ.4.D0) GO TO 420
      IF (IW1.LT.9) GO TO 120
      N2 = 9
      DO 400 I = 1, N
         IF (Y(I).NE.0.D0 .OR. W(I,8).GT.0.D0) GO TO 400
C        FLOOR ZERO AND SOLUTION ZERO IN SAME COMPONENT
         IND = 1
         CIN(1) = -6.D0
         GO TO 2620
  400 CONTINUE
  420 IF ( .NOT. TEST) GO TO 440
      IF (COMM(4).LT.0.D0) GO TO 720
      GO TO 860
  440 IF (CONST(1).EQ.0.D0) GO TO 480
      IF (CONST(1).EQ.1.D0) GO TO 460
C     CONST(1) INVALID
      IND = 1
      CIN(1) = -7.D0
      GO TO 2620
  460 EXPI = FAC9
      EXP = FAC8
  480 IF (CONST(2).GE.0.D0 .AND. CONST(3).GE.0.D0) GO TO 520
C     CONST(2) OR CONST(3) INVALID
  500 IND = 1
      CIN(1) = -8.D0
      GO TO 2620
  520 IF (CONST(2).EQ.0.D0) CONST(2) = FAC3
      IF (CONST(3).EQ.0.D0) CONST(3) = FAC4
      IF (CONST(2).LE.1.D0/CONST(3) .OR. CONST(2).LE.1.D0 .OR. CONST(3)
     *    .GE.1.D0) GO TO 500
      IF (COMM(1).GE.0.D0) GO TO 540
C     COMM(1) OUT OF RANGE
      IND = 1
      CIN(1) = -9.D0
      GO TO 2620
  540 IF (COMM(1).GT.0.D0) CALLS = .TRUE.
      IF (COMM(2).EQ.0.D0) GO TO 580
      IF (COMM(2).GT.0.D0) GO TO 560
C     COMM(2) OUT OF RANGE
      IND = 1
      CIN(1) = -10.D0
      GO TO 2620
  560 IF (COUT(6).LT.COMM(2)) GO TO 580
C     INITIAL VECTOR TOO LARGE
      IND = 1
      COMM(2) = 0.D0
      CIN(1) = -11.D0
      GO TO 2620
  580 IF (COMM(3).EQ.0.D0) GO TO 640
      IF (IW1.LT.10) GO TO 120
      N2 = 8
      IF (CIN(2).EQ.3.D0) N2 = 11
      IF (IW1.LT.N2) GO TO 120
      DO 600 I = 1, N
         IF (Y(I).NE.W(I,9)) GO TO 600
C        INITIAL VECTOR ATTAINS GIVEN VALUE
         IND = 1
         COMM(3) = 0.D0
         CIN(1) = -12.D0
         GO TO 2620
  600 CONTINUE
      DO 620 I = 1, N
         W(I,10) = Y(I)
  620 CONTINUE
C
C     TEST INPUT PARAMETERS FOR INTERRUPT ON POSITION
C
  640 IF (COMM(4).GE.0.D0) GO TO 800
      IF (CIN(1).GT.1.D0) GO TO 720
      IF (X.EQ.XEND .AND. COMM(5).EQ.X) GO TO 880
      IF (X.NE.XEND) GO TO 680
C     VALUE OF CIN(1) DOES NOT PERMIT EXTRAPOLATION
  660 IND = 1
      CIN(1) = -13.D0
      GO TO 780
  680 IF ((X.GE.COMM(5) .AND. XEND.LT.COMM(5)) .OR. (X.LT.COMM(5)
     *    .AND. XEND.GE.COMM(5))) GO TO 800
      IF (INTER) GO TO 700
C     VALUE OF CIN(1) DOES NOT PERMIT EXTRAPOLATION
      GO TO 660
C     ORDER OF POINTS X,COUT(4) AND COUT(5)  INCORRECT
  700 IND = 1
      CIN(1) = -14.D0
      GO TO 780
  720 IF (X.EQ.COUT(5) .OR. X.EQ.COUT(4)) GO TO 700
      IF ((X.GE.COUT(4) .AND. COUT(4).LT.COUT(5)) .OR. (X.LT.COUT(4)
     *     .AND. COUT(4).GE.COUT(5))) GO TO 700
      IF (X.EQ.XEND) GO TO 880
      IF (SIGN(1.D0,XEND-X).NE.SIGN(1.D0,X-COUT(5))) GO TO 700
      IF (COMM(5).EQ.COUT(5) .OR. COMM(5).EQ.XEND) GO TO 740
      IF (SIGN(1.D0,COMM(5)-COUT(5)).NE.SIGN(1.D0,COMM(5)-XEND))
     *    GO TO 740
C     INTERRUPT POINT NOT IN RANGE
      IND = 1
      CIN(1) = -15.D0
      GO TO 780
  740 IF (COMM(5).EQ.COUT(5) .OR. COMM(5).EQ.X) GO TO 760
      IF (SIGN(1.D0,COMM(5)-COUT(5)).EQ.SIGN(1.D0,COMM(5)-X))
     *    GO TO 800
  760 CIN(1) = 6.D0
  780 COMM(4) = 0.D0
      GO TO 2620
C
C     SET DEFAULTS
C
  800 IF (TEST) GO TO 860
      IF (ABS(CIN(4)).GT.ABS(CIN(3)) .OR. CIN(4).EQ.0.D0) GO TO 820
C     USER SET HMIN.GT.HMAX
      IND = 1
      CIN(1) = -16.D0
      GO TO 2620
  820 COUT(13) = MAX(2.D0*COUT(11)*ABS(XEND-X),ABS(CIN(3)))
      FAC7 = 0.5D0
      IF (CIN(1).GE.2.D0) FAC7 = 1.D0
      COUT(14) = FAC7*ABS(XEND-X)
      IF (CIN(4).EQ.0.D0) GO TO 840
      COUT(14) = MIN(ABS(COUT(14)),ABS(CIN(4)))
  840 IF (CIN(1).GE.2.D0) GO TO 860
      CALL FCN(N,X,Y,W(1,1),COEFFN,COEFF1,M,ARR)
      IF (CALLS) COMM(1) = COMM(1) - 1.D0
C
C     INITIALISATION
C
  860 IF (X.NE.XEND) GO TO 920
      IF (CIN(1).LE.1.D0) GO TO 900
      IF (CIN(1).NE.2.D0 .OR. CIN(6).NE.0.D0) GO TO 900
C     REPEATED CALL WITH X=XEND
      IND = 1
      CIN(1) = -17.D0
      GO TO 2620
C
C     RETURN WHEN X=XEND
C
  880 COMM(4) = 0.D0
  900 CIN(5) = 0.D0
      CIN(6) = 0.D0
      CIN(1) = 2.D0
      GO TO 2620
C
C     X.NE.XEND INITIALLY
C
  920 IEND = 0
      REDUCE = .FALSE.
      RAT1 = 1.D0
      IF (CIN(1).EQ.0.D0) GO TO 960
      IF (CIN(1).EQ.1.D0) CIN(6) = CIN(5)
      IF (CIN(6).EQ.0.0D0) GO TO 940
      IF (SIGN(1.D0,CIN(6)).NE.SIGN(1.D0,XEND-X)) CIN(6) = 0.D0
  940 IF (CIN(2).EQ.2.D0) FAC2 = SMALL/COUT(11)
      IF (ABS(CIN(6)).GE.COUT(13) .AND. CIN(1).GE.2.D0) GO TO 1500
  960 ISTEP = 0
      START = .TRUE.
      IF (CIN(6).NE.0.D0) GO TO 1160
C
C     ESTIMATE STEP
C
      CIN(1) = 1.D0
      FNORM = 0.D0
      YNORM = 0.D0
      DO 980 I = 1, N
         YNORM = MAX(YNORM,ABS(Y(I)))
         FNORM = MAX(FNORM,ABS(W(I,1)))
  980 CONTINUE
      TOLEST = TOL
      J = INT(CIN(2)+0.1D0) + 1
      GO TO (1000,1140,1020,1040,1080) J
 1000 S = MAX(1.D0,YNORM)
      GO TO 1120
 1020 S = MAX(YNORM,FAC2)
      GO TO 1120
 1040 DO 1060 I = 1, N
         T = W(I,7)*MAX(W(I,8),ABS(Y(I)))
         IF (I.EQ.1) S = T
         S = MIN(S,T)
 1060 CONTINUE
      GO TO 1120
 1080 DO 1100 I = 1, N
         IF (I.EQ.1) S = W(1,7)
         S = MIN(S,W(I,7))
 1100 CONTINUE
 1120 TOLEST = TOLEST*S
 1140 S = SQRT(COUT(11))
      CIN(6) = (XEND-X)*(TOLEST*COUT(11))**EXPI*MAX(S,YNORM)
     *         /MAX(S,FNORM)
 1160 DO 1180 I = 1, N
         W(I,4) = Y(I)
         W(I,5) = W(I,1)
 1180 CONTINUE
      IF (ABS(CIN(6)).LT.COUT(13)) CIN(6) = SIGN(COUT(13),XEND-X)
      IF (ABS(CIN(6)).GT.COUT(14)) CIN(6) = SIGN(COUT(14),CIN(6))
      ISIG = 1
      GO TO 1600
C
C     RETURN FOR INITIAL STEP
C
 1200 IF (COMM(1).GT.0.D0 .OR. .NOT. CALLS) GO TO 1260
C     TOO MANY FCN CALLS TAKEN STARTING
 1220 IND = 7
      DO 1240 I = 1, N
         Y(I) = W(I,4)
 1240 CONTINUE
      COUT(9) = COUT(9) + 1.D0
      GO TO 2620
 1260 IF (ABS(HEST).LT.ABS(CIN(6))) GO TO 1280
C
C     ESTIMATED STEP LARGER THAN INITIAL STEP
C
      IF (ABS(CIN(6))*CONST(2).GT.COUT(14)) GO TO 1420
      IF (ABS(HEST).GT.COUT(14)) HEST = SIGN(COUT(14),HEST)
      IF (REDUCE) GO TO 1420
      IF (ABS(HEST).LT.CONST(2)*CONST(2)*ABS(CIN(6))) GO TO 1420
      CIN(6) = HEST
      GO TO 1460
 1280 IF (ABS(CIN(6)).GT.COUT(13)) GO TO 1340
C     ERROR TOLERANCE TOO SMALL FOR INITIAL STEP
 1300 COUT(9) = COUT(9) + 1.D0
      DO 1320 I = 1, N
         Y(I) = W(I,4)
 1320 CONTINUE
      IND = 4
      GO TO 2620
C
C     ESTIMATED STEP SMALLER THAN INITIAL STEP
C
 1340 T = CIN(6)/CONST(2)
      IF (ABS(HEST).GT.ABS(T)) GO TO 1420
      T = T*FAC1
      IF (ABS(HEST).LT.ABS(T)) HEST = T
      CIN(6) = HEST
      IF (ABS(CIN(6)).LT.COUT(13)) CIN(6) = SIGN(COUT(13),XEND-X)
      GO TO 1440
C
C     INSIGNIFICANT ERROR ESTIMATE ON INITIAL STEP
C
 1360 CONTINUE
      IF (REDUCE) GO TO 1380
      IF (ABS(CIN(6)).EQ.COUT(14)) GO TO 1400
      CIN(6) = CIN(6)*FAC
      IF (ABS(CIN(6)).GT.COUT(14)) CIN(6) = SIGN(COUT(14),CIN(6))
      GO TO 1460
 1380 IF (ISIG.EQ.3) GO TO 1300
      ISIG = 3
      CIN(6) = CIN(6)/FAC3
      GO TO 1440
C
C     INITIAL STEP ACCEPTED
C
 1400 HEST = CIN(6)
 1420 START = .FALSE.
      GO TO 2000
C
C     INITIAL STEP REJECTED TRY AGAIN
C
 1440 REDUCE = .TRUE.
 1460 DO 1480 I = 1, N
         Y(I) = W(I,2)
 1480 CONTINUE
      COUT(9) = COUT(9) + 1.D0
      GO TO 1560
C
C     TAKE A STEP, CHECK AGAINST LENGTH OF RANGE
C
 1500 IF (ISTEP.LT.2) ISTEP = ISTEP + 1
 1520 IF (D02PAY(X,CIN(6),XEND)*SIGN(1.0D0,XEND-X).LT.0.0D0)
     *    GO TO 1560
      IF (ISTEP.EQ.-1) GO TO 1540
      CIN(6) = XEND - X
      RAT1 = ABS(COUT(4)-X)/ABS(CIN(6))
 1540 IEND = 1
 1560 DO 1580 I = 1, N
         W(I,4) = W(I,2)
         W(I,5) = W(I,3)
 1580 CONTINUE
 1600 CALL D02KDZ(X,CIN(6),N,Y,FCN,W,IW,N2,COEFFN,COEFF1,M,ARR)
      IF (CALLS) COMM(1) = COMM(1) - 4.D0
C
C     ESTIMATE NEW STEP
C
      IF (CIN(2).GE.3.D0) GO TO 1800
      S = 0.D0
      T = 0.D0
      J = 0
      DO 1620 I = 1, N
         IF (S.GT.ABS(W(I,6))) GO TO 1620
         J = I
         T = W(I,N2)
         S = ABS(W(I,6))
 1620 CONTINUE
      IF (T.EQ.0.D0) GO TO 1680
 1640 RAT = RAT1
      IF (START) GO TO 1360
      IF (J.NE.INT(COUT(10)+0.1D0)) GO TO 1980
      HEST = CIN(6)
      IF (ABS(HEST).EQ.COUT(14) .OR. IEND.EQ.1) GO TO 1980
      IF (ABS(HEST).LT.ABS(HEST1)) GO TO 1660
      HEST1 = HEST
      RAT = CONST(2)
      GO TO 1980
 1660 IND = 3
      GO TO 2620
 1680 J = 0
      K = INT(CIN(2)+0.1D0) + 1
      GO TO (1720,1700,1760) K
 1700 RAT = (TOL/S)**EXP
      GO TO 1980
 1720 T = 1.D0
      DO 1740 I = 1, N
         T = MAX(T,ABS(Y(I)))
 1740 CONTINUE
      RAT = (TOL*T/S)**EXP
      GO TO 1980
 1760 T = FAC2
      DO 1780 I = 1, N
         T = MAX(T,ABS(Y(I)))
 1780 CONTINUE
      RAT = (TOL*T/S)**EXP
      GO TO 1980
 1800 IF (CIN(2).EQ.4.D0) GO TO 1840
      DO 1820 I = 1, N
C        RELATIVE ERROR FAILURE
         IF (W(I,8).GT.0.0D0) GO TO 1820
         IF ((Y(I).GE.0.0D0 .AND. W(I,2).GE.0.0D0) .OR. (Y(I)
     *       .LT.0.0D0 .AND. W(I,2).LT.0.0D0)) GO TO 1820
         IND = 6
         GO TO 2620
 1820 CONTINUE
 1840 J = -1
      I1 = 0
      DO 1960 I = 1, N
         IF (W(I,6).NE.0.D0) GO TO 1880
         IF (W(I,N2).EQ.0.D0) GO TO 1860
         IF (I1.EQ.0) J = I
         GO TO 1960
 1860    IF (J.LT.0) J = I
         GO TO 1960
 1880    IF (CIN(2).EQ.4.D0) GO TO 1900
         T = MAX(W(I,8),ABS(Y(I)))*W(I,7)/ABS(W(I,6))
         GO TO 1920
 1900    T = W(I,7)/ABS(W(I,6))
 1920    IF (I1.EQ.0) GO TO 1940
         IF (S.LE.T) GO TO 1960
 1940    S = T
         J = 0
         I1 = 1
         IF (W(I,N2).EQ.1.D0) J = I
 1960 CONTINUE
      IF (J.NE.0) GO TO 1640
      RAT = (TOL*S)**EXP
 1980 HEST = RAT*CIN(6)
      IF (START) GO TO 1200
C
C     TEST ESTIMATED STEP
C
 2000 IF (ABS(HEST).GE.ABS(CIN(6))) GO TO 2140
C
C     STEP REJECTED
C
      RAT1 = 1
      IEND = 0
      IF (ISTEP.LT.0) GO TO 2100
      COUT12 = 1.D0
      IF (X.EQ.XORIG) COUT(12) = 1.D0
      HEST = HEST*CONST(3)
      COUT(9) = COUT(9) + 1.D0
      DO 2020 I = 1, N
         Y(I) = W(I,2)
         W(I,3) = W(I,5)
         W(I,2) = W(I,4)
 2020 CONTINUE
      T = CIN(6)/CONST(2)
      IF (ABS(HEST).LT.ABS(T)) HEST = T
      CIN(6) = HEST
      IF (ABS(CIN(6)).GE.COUT(13)) GO TO 2060
C     STEP LENGTH TOO SMALL
 2040 IND = 2
      GO TO 2620
 2060 IF (COMM(1).GT.0.D0 .OR. .NOT. CALLS) GO TO 1520
      IF (ISTEP.LT.2) GO TO 1220
C     TOO MANY FCN CALLS
 2080 IND = 5
      GO TO 2620
C
C     STEP REJECTED AFTER INTERRUPT IN FIRST STEP
C
 2100 IF (COMM(1).LT.0.D0) GO TO 1220
      CIN(6) = (X-COUT(4))*0.5D0
      IF (ABS(CIN(6)).LT.COUT(13)) GO TO 2040
      X = COUT(4)
      DO 2120 I = 1, N
         Y(I) = W(I,4)
         W(I,1) = W(I,5)
 2120 CONTINUE
      COUT(9) = COUT(9) + 2.D0
      INIT = .TRUE.
      ISTEP = 0
      GO TO 1600
C
C     STEP ACCEPTED
C
 2140 COUT(5) = COUT(4)
      COUT(4) = X
      X = X + CIN(6)
      HEST1 = 0.D0
      COUT(10) = DBLE(J)
      COUT(8) = COUT(8) + 1.D0
      IF ( .NOT. INIT) GO TO 2160
      INIT = .FALSE.
      CIN(5) = CIN(6)
      COUT(1) = CIN(6)
      COUT(2) = CIN(6)
      GO TO 2180
 2160 IF (ABS(CIN(6)).LT.ABS(COUT(1))) COUT(1) = CIN(6)
      IF (ABS(CIN(6)).GT.ABS(COUT(2))) COUT(2) = CIN(6)
 2180 IF (HEST.EQ.CIN(6) .AND. ABS(HEST).EQ.COUT(14)) GO TO 2260
      HEST = HEST*CONST(3)
      IF (ABS(CIN(6)).EQ.COUT(14) .OR. IEND.EQ.1) COUT(3) = COUT(3) +
     *    1.D0
      IF (COUT12.EQ.0.D0) GO TO 2200
      COUT12 = 0.D0
      HEST = 0.5D0*(HEST+CIN(6))
 2200 IF (IEND.EQ.1) GO TO 2220
      T = CONST(2)*CIN(6)
      IF (ABS(HEST).GT.ABS(T)) HEST = T
      IF (ABS(HEST).GT.COUT(14)) HEST = SIGN(COUT(14),HEST)
      IF (ABS(HEST).LT.COUT(13)) HEST = SIGN(COUT(13),XEND-X)
      GO TO 2240
 2220 CIN(6) = CIN(6)*RAT1
      IF (ABS(CIN(6)).GE.ABS(HEST)) GO TO 2260
 2240 CIN(6) = HEST
 2260 S = 0.D0
      DO 2280 I = 1, N
         S = MAX(S,ABS(Y(I)))
 2280 CONTINUE
      COUT(6) = MAX(S,COUT(6))
      COUT(7) = MIN(S,COUT(7))
      CALL FCN(N,X,Y,W(1,1),COEFFN,COEFF1,M,ARR)
      IF (CIN(1).EQ.0.D0) GO TO 2580
      IF (COMM(3).EQ.0.D0) GO TO 2320
      DO 2300 I = 1, N
         IF (ABS(Y(I)-W(I,9)).LT.ABS(W(I,10)-W(I,9))) W(I,10) = Y(I)
 2300 CONTINUE
 2320 IF ( .NOT. CALLS) GO TO 2340
      COMM(1) = COMM(1) - 1.D0
      IF (COMM(1).LE.0.D0) GO TO 2080
 2340 IF (IEND.EQ.1) GO TO 2600
      IF (ISTEP.GE.0) GO TO 2360
      ISTEP = -ISTEP
      GO TO (2460,2440,2620,2520) ISTEP
C
C     TEST INTERRUPTS
C
 2360 IF (COMM(3).EQ.0.D0) GO TO 2420
      DO 2400 I = 1, N
         IF (Y(I).EQ.W(I,9)) GO TO 2380
         IF ((W(I,2).GE.W(I,9) .AND. Y(I).GE.W(I,9)) .OR. (W(I,2)
     *       .LT.W(I,9) .AND. Y(I).LT.W(I,9))) GO TO 2400
C        COMPONENTS ACHIEVE GIVEN VALUE
 2380    CIN(1) = 3.D0
         IF (ISTEP.EQ.0) GO TO 2540
         GO TO 2460
 2400 CONTINUE
C     NORM OF SOLUTION TOO LARGE TOO CONTINUE
 2420 IF (COMM(2).EQ.0.D0 .OR. COUT(6).LT.COMM(2)) GO TO 2480
      CIN(1) = 4.D0
      IF (ISTEP.EQ.0) GO TO 2540
 2440 COMM(2) = 0.0D0
      GO TO 2620
 2460 COMM(3) = 0.0D0
      GO TO 2620
 2480 IF (COMM(4).EQ.0.D0) GO TO 2560
      IF (COMM(4).LT.0.D0) GO TO 2500
C     INTERRUPT EVERY STEP
      CIN(1) = 5.D0
      IF (ISTEP.EQ.0) GO TO 2540
      GO TO 2620
 2500 IF ((COMM(5).GE.X .AND. XEND.GE.X) .OR. (COMM(5)
     *    .LT.X .AND. XEND.LT.X)) GO TO 2560
C     INTERRUPT AT SPECIFIED POINT
      CIN(1) = 6.D0
      IF (ISTEP.EQ.0) GO TO 2540
 2520 COMM(4) = 0.0D0
      GO TO 2620
 2540 IF (COMM(1).LT.0.D0) GO TO 1220
      ISTEP = -INT(CIN(1)-1.5D0)
      GO TO 1520
 2560 IF (COMM(1).LT.0.D0) GO TO 2080
 2580 IF (IEND.EQ.0) GO TO 1500
C     NORMAL RETURN
 2600 CIN(1) = 2.D0
      X = XEND
C
C     RETURN TO MAIN PROGRAM
C
 2620 IF (CALLS .AND. COMM(1).EQ.0.D0) COMM(1) = -1.D0
      IFAIL = P01ABF(IFAIL,IND,SRNAME,0,P01REC)
      IF (CIN(1).GE.0.D0 .AND. IFAIL.GT.0) CIN(1) = 8.D0
      RETURN
      END
      SUBROUTINE D02KDZ(X,H,N,Y,FCN,W,IW1,IW2,COEFFN,COEFF1,M,ARR)
C     MARK 7 RELEASE. NAG COPYRIGHT 1978.
C     MARK 8 REVISED. IER-227 (APR 1980).
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     MARK 13 REVISED. USE OF MARK 12 X02 FUNCTIONS (APR 1988).
C     COEFF1, COEFFN, FCN
C     .. Scalar Arguments ..
      DOUBLE PRECISION  H, X
      INTEGER           IW1, IW2, M, N
C     .. Array Arguments ..
      DOUBLE PRECISION  ARR(M), W(IW1,IW2), Y(N)
C     .. Subroutine Arguments ..
      EXTERNAL          COEFF1, COEFFN, FCN
C     .. Local Scalars ..
      DOUBLE PRECISION  EPS, S
      INTEGER           I, N1
C     .. Local Arrays ..
      DOUBLE PRECISION  C(10)
C     .. External Functions ..
      DOUBLE PRECISION  X02AJF
      EXTERNAL          X02AJF
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, MAX
C     .. Executable Statements ..
      EPS = X02AJF()
      C(1) = 1.D0/6.D0
      C(2) = 1.D0/3.D0
      C(3) = 0.125D0
      C(4) = 0.375D0
      C(5) = 0.5D0
      C(6) = 1.5D0
      C(7) = 2.0D0
      C(8) = 2.D0/3.D0
      C(9) = 0.2D0
      C(10) = 4.D0/3.D0
      N1 = 6
      IF (IW2.EQ.4) N1 = 4
      IF (IW2.EQ.6) N1 = 5
      DO 20 I = 1, N
         W(I,3) = Y(I) + C(2)*H*W(I,1)
   20 CONTINUE
      CALL FCN(N,X+C(2)*H,W(1,3),W(1,N1),COEFFN,COEFF1,M,ARR)
      DO 40 I = 1, N
         IF (N1.EQ.5) W(I,4) = 0.5D0*(W(I,N1)-W(I,1))*C(2)*H
         W(I,3) = Y(I) + C(1)*H*(W(I,1)+W(I,N1))
   40 CONTINUE
      CALL FCN(N,X+C(2)*H,W(1,3),W(1,N1),COEFFN,COEFF1,M,ARR)
      DO 60 I = 1, N
         W(I,2) = Y(I) + H*(C(3)*W(I,1)+C(4)*W(I,N1))
   60 CONTINUE
      CALL FCN(N,X+C(5)*H,W(1,2),W(1,3),COEFFN,COEFF1,M,ARR)
      DO 100 I = 1, N
         IF (N1.EQ.4) GO TO 80
         W(I,IW2) = -C(2)*W(I,1) - C(10)*W(I,3) + C(6)*W(I,N1)
   80    W(I,2) = Y(I) + H*(C(5)*W(I,1)-C(6)*W(I,N1)+C(7)*W(I,3))
  100 CONTINUE
      CALL FCN(N,X+H,W(1,2),W(1,N1),COEFFN,COEFF1,M,ARR)
      DO 140 I = 1, N
         W(I,2) = Y(I)
         Y(I) = Y(I) + H*(C(1)*(W(I,1)+W(I,N1))+C(8)*W(I,3))
         IF (N1.EQ.4) GO TO 120
         W(I,3) = W(I,N1)
         W(I,N1) = C(9)*H*(W(I,IW2)+C(1)*W(I,N1))
         S = 0.D0
         IF (ABS(W(I,N1)).LE.30.D0*C(9)*EPS*ABS(H)*MAX(ABS(W(I,IW2))
     *       ,C(1)*ABS(W(I,3)))) S = 1.D0
         W(I,IW2) = S
  120    W(I,3) = W(I,1)
  140 CONTINUE
      RETURN
      END
      SUBROUTINE D02KEF(XPOINT,NXP,IC1,COEFFN,BDYVAL,K,TOL,ELAM,DELAM,
     *                  HMAX,MAXIT,MAXFUN,MONIT,REPORT,IFAIL)
C     MARK 8 RELEASE. NAG COPYRIGHT 1979.
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     MARK 13 REVISED. USE OF MARK 12 X02 FUNCTIONS (APR 1988).
C
C     MEANING OF COMMON VARIABLES:
C     ZER,ONE,TWO,PI - CONSTANTS WITH INDICATED VALUES
C     LAMDA - CURRENT VALUE OF EIGENVALUE ESTIMATE. SET AND USED IN
C     EIGODE. USED IN AUX,OPTSC.
C     PSIGN   - SAMPLE P(X) VALUE. SET IN PRUFER, USED IN AUX.
C     MINSC - MINIMUM ALLOWED SCALEFACTOR. SET IN EIGODE, USED IN
C     OPTSC.
C     BP    - USED BY THE POLYGONAL SCALING METHOD. THE CURRENT
C     VALUE OF THE DERIVATIVE OF B(X). SET IN PRUFER, USED I
C
C     BDYVAL, COEFFN, MONIT, REPORT
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='D02KEF')
C     .. Scalar Arguments ..
      DOUBLE PRECISION  DELAM, ELAM, TOL
      INTEGER           IC1, IFAIL, K, MAXFUN, MAXIT, NXP
C     .. Array Arguments ..
      DOUBLE PRECISION  HMAX(2,NXP), XPOINT(NXP)
C     .. Subroutine Arguments ..
      EXTERNAL          BDYVAL, COEFFN, MONIT, REPORT
C     .. Scalars in Common ..
      DOUBLE PRECISION  BP, LAMDA, MINSC, ONE, PI, PSIGN, TWO, ZER
      INTEGER           JINT
C     .. Arrays in Common ..
      DOUBLE PRECISION  YL(3), YR(3)
C     .. Local Scalars ..
      DOUBLE PRECISION  D, DELAM1, DETOL, DETOLL, DETOLR, DMIN, EL0,
     *                  EL1, F0, F1, GAMMA, SIGNSI, TEMP0, TEMP1, TEMP2,
     *                  TEMP3, TEMP4, TEMP5, TEMP6, TEMP7, TEMP8, TENU,
     *                  TOL2, V3BDYL, V3BDYR, V3L, V3R, X, XL, XOLD, XR,
     *                  ZETA
      INTEGER           I, IBACK, IC, IFAIL1, IFLAG, ILOOP, IND, ISTAT1,
     *                  ISTATE, ITOP, IXP, MAXFN1, NINT, NXP1
C     .. Local Arrays ..
      DOUBLE PRECISION  C(17), F(15), V(21), VL(7)
      CHARACTER*1       P01REC(1)
C     .. External Functions ..
      DOUBLE PRECISION  D02KDS, D02KDU, X01AAF, X02AJF
      INTEGER           P01ABF
      EXTERNAL          D02KDS, D02KDU, X01AAF, X02AJF, P01ABF
C     .. External Subroutines ..
      EXTERNAL          C05AZF, D02KDT, D02KDV, D02KEZ
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, MAX, MIN, SIGN, LOG, DBLE
C     .. Common blocks ..
      COMMON            /AD02KD/ZER, ONE, TWO, PI, LAMDA, PSIGN, MINSC,
     *                  BP, YL, YR, JINT
C     .. Executable Statements ..
      ZER = 0.D0
      ONE = 1.D0
      TWO = 2.D0
      PI = X01AAF(0.0D0)
      TENU = 10.D0*X02AJF()
C
C     PARAMETER CHECKS, DEFAULT VALUES ETC.
      IF (NXP.LT.4) GO TO 320
      IF (K.LT.0) GO TO 340
      IF (TOL.LE.ZER) GO TO 360
      NINT = NXP - 3
      XL = XPOINT(2)
      XR = XPOINT(NXP-1)
C     CHECK XPOINT(I) IN ASC. ORDER AND FIND POSN. OF MATCHING
C     POINT.
      DMIN = TWO*ABS(XR-XL)
      NXP1 = NXP - 1
      XOLD = XPOINT(1)
      DO 40 I = 2, NXP1
         X = XPOINT(I)
         IF (X.LT.XOLD) GO TO 380
         D = ABS((XL-X)+(XR-X))
         IF (D.GT.DMIN) GO TO 20
         IC = I
         DMIN = D
   20    XOLD = X
   40 CONTINUE
      IF ((IC1.LT.2) .OR. (IC1.GT.NXP1)) IC1 = IC
      IC = IC1
      I = NXP
      IF (XPOINT(NXP).LT.XOLD) GO TO 380
C
      MINSC = ONE/TENU
      GAMMA = ZER
      IF (XL.EQ.XR) GO TO 60
      MINSC = 4.D0/ABS(XR-XL)
      GAMMA = (XPOINT(IC)-XL)/(XR-XL)
   60 CONTINUE
C
      DO 80 I = 1, NXP
         HMAX(2,I) = ZER
   80 CONTINUE
C
      TOL2 = TOL/TWO
      DELAM1 = DELAM
      IF (DELAM1.EQ.ZER) DELAM1 = .25D0*MAX(ONE,ABS(ELAM))
C     HAVE REPORT ROUTINE SWITCHED OFF INITIALLY
      V(6) = ZER
C     INITIAL RHO VALUES AT XL,XR TAKEN AS ZERO
C     ON FIRST ITERATION
      V3L = ZER
      V3R = ZER
C
C     INITIAL VALUES OF TOLERANCE PARAMS LOG(ETA) AT XL,XR.
      DETOLL = LOG(1.D-4)
      DETOLR = DETOLL
C     INITIALIZE FAIL-FLAG FOR MISS DISTANCE CODE.
      IFAIL1 = 0
      EL1 = ELAM
      EL0 = ELAM + DELAM1
C     INITIAL EVALUATION, F=F(EL1) .
      ISTATE = 1
      LAMDA = EL1
      IBACK = 1
      GO TO 520
  100 IF (IFAIL1.GT.0) GO TO 460
      F1 = F(1)
C     ENTER BRACKETING LOOP.
      IF (F1.EQ.ZER) GO TO 180
C     THROUGHOUT THE ROOTFINDING, MAINTAIN EL1,EL0 AS TWO MOST
C     RECENT
C     APPROXIMATIONS, WITH ABS(F(EL0)).GE.ABS(F(EL1))
C
      IBACK = 2
      LAMDA = EL0
      ILOOP = 0
  120 ILOOP = ILOOP + 1
      GO TO 520
  140 IF (IFAIL1.GT.0) GO TO 440
      EL0 = LAMDA
      F0 = F(1)
      IF (ABS(F0).GT.ABS(F1)) GO TO 160
      EL0 = EL1
      F0 = F1
      EL1 = LAMDA
      F1 = F(1)
  160 IF (F1.EQ.ZER) GO TO 180
      IF (F0*F1.LT.ZER) GO TO 200
      LAMDA = EL1 - TWO*(EL0-EL1)
      IF (ILOOP.LE.8) GO TO 120
C     HERE HAVE FAILED TO BRACKET ROOT AFTER 10 EVALUATIONS
      IFAIL1 = 5
      GO TO 440
C     ON FINDING FMISS EXACTLY 0 ON
C     1ST EVAL OR WHEN BRACKETING
  180 EL0 = EL1
      GO TO 280
C     ROOT BETWEEN EL0 AND EL1..  START CALLING ROOTFINDER
  200 IFLAG = 1
      ISTATE = 2
      IBACK = 3
      IND = -1
      C(1) = F0
  220 CALL C05AZF(EL1,EL0,F1,TOL2,0,C,IND,IFLAG)
      IF (IND.EQ.0) GO TO 260
      IF (IND.NE.4 .OR. IFLAG.NE.1) GO TO 420
      LAMDA = EL1
      GO TO 520
  240 IF (IFAIL1.GT.0) GO TO 440
      F1 = F(1)
      GO TO 220
C     ON EXIT FROM LOOP TEST FAILURE PARAM OF C05AZF
  260 IF (IFLAG.NE.0) GO TO 420
C     SWITCH REPORT ROUTINE ON AND
C     REPEAT INTEGRATION
  280 V(6) = ONE
      LAMDA = EL1
      IBACK = 4
      GO TO 520
  300 IF (IFAIL1.GT.0) GO TO 460
      IFAIL1 = 0
      GO TO 440
C     ***
C     ERROR PROCESSING FOR MAIN ROUTINE
C     ***
C     PARAMETER ERROR NXP,K OR TOL
  320 HMAX(2,1) = 1.D0
      GO TO 400
  340 HMAX(2,1) = 2.D0
      GO TO 400
  360 HMAX(2,1) = 3.D0
      GO TO 400
  380 HMAX(2,1) = 4.D0
      HMAX(2,2) = DBLE(I)
  400 IFAIL1 = 1
      HMAX(2,2) = ZER
      IF (HMAX(2,1).EQ.4.D0) HMAX(2,2) = DBLE(I)
      GO TO 500
  420 IFAIL1 = 12
      IF (IFLAG.EQ.4) IFAIL1 = 10
      IF (IFLAG.EQ.5) IFAIL1 = 9
  440 ELAM = EL1
      DELAM = ABS(F(2)) + ABS(EL0-EL1)
  460 HMAX(2,1) = ZER
      HMAX(2,2) = ZER
      IF (IFAIL1.EQ.0) GO TO 500
      IF (IFAIL1.EQ.5) DELAM = EL0 - EL1
      IF (IFAIL1.NE.11) GO TO 480
      HMAX(2,1) = TEMP1
      HMAX(2,2) = TEMP2
  480 IF (IFAIL1.NE.12) GO TO 500
      HMAX(2,1) = DBLE(IFLAG)
      HMAX(2,2) = DBLE(IND)
  500 IFAIL = P01ABF(IFAIL,IFAIL1,SRNAME,0,P01REC)
      RETURN
C     ***
C     CODE TO COMPUTE MISS DISTANCE.
C     ***
  520 F(3) = ZER
  540 F(3) = F(3) + ONE
      IF (IFAIL1.LT.0) GO TO 960
      MAXIT = MAXIT - 1
      MAXFN1 = MAXFUN
      F(1) = ZER
      F(2) = ZER
      DO 560 I = 4, 15
         F(I) = ZER
  560 CONTINUE
      HMAX(1,NXP-1) = ZER
      HMAX(1,NXP) = ZER
      YL(3) = ZER
      YR(3) = ZER
      CALL BDYVAL(XL,XR,LAMDA,YL,YR)
C     SET PARAMETERS APPLICABLE TO SHOOTING EITHER DIRECTION.
      V(7) = MAX(ONE,ABS(LAMDA))*TOL2
C     SET PARAMS FOR SHOOTING C FROM XL.
      JINT = 1
      V(1) = D02KDU(XL,COEFFN)
      CALL D02KDT(V,YL(1),YL(2),-1,IFAIL1)
      IF (IFAIL1.NE.0) GO TO 900
      V3BDYL = V(3)
      V(3) = V3L
      V(4) = ZER
      V(5) = DETOLL + V(3)
      X = XL
      IF (V(6).NE.ZER) CALL REPORT(X,V,0)
      IF (IC.LT.3) GO TO 600
      DO 580 IXP = 3, IC
         JINT = IXP - 2
         CALL D02KEZ(X,XPOINT(IXP),V,COEFFN,MAXFN1,HMAX(1,JINT)
     *               ,REPORT,IFAIL1)
         F(5) = F(5) + V(15)
         F(6) = F(6) + V(16)
         F(7) = F(7) + V(10)
         IF (IFAIL1.NE.0) GO TO 760
  580 CONTINUE
C
  600 DO 620 I = 1, 7
         VL(I) = V(I)
  620 CONTINUE
C
C     SET PARAMS FOR SHOOTING C FROM XR.
      JINT = NINT
      V(1) = D02KDU(XR,COEFFN)
      CALL D02KDT(V,YR(1),YR(2),K,IFAIL1)
      IF (IFAIL1.NE.0) GO TO 920
      V3BDYR = V(3)
      V(3) = V3R
      V(4) = ZER
      V(5) = DETOLR + V(3)
      X = XR
      IF (V(6).NE.ZER) CALL REPORT(X,V,NINT+1)
      ITOP = NXP - IC
      IF (2.GT.ITOP) GO TO 660
      DO 640 I = 2, ITOP
         IXP = NXP - I
         JINT = IXP - 1
         CALL D02KEZ(X,XPOINT(IXP),V,COEFFN,MAXFN1,HMAX(1,JINT)
     *               ,REPORT,IFAIL1)
         F(5) = F(5) + V(15)
         F(6) = F(6) + V(16)
         F(7) = F(7) + V(10)
         IF (IFAIL1.NE.0) GO TO 760
  640 CONTINUE
C
  660 CONTINUE
C     HAVE NOW SHOT C FROM XL, RESULTS IN VL, C FROM XR, RESULTS IN
C     V.
C     CONVERT TO SAME SCALE AND COMPUTE MISS-DISTANCE.
      CALL D02KDV(VL,V(1),IFAIL1)
      F(1) = (VL(2)-V(2))/(PI+PI)
C     COMPUTE ERROR ESTIMATE INFO.
      TEMP0 = VL(5) - VL(3)
      TEMP1 = V(5) - V(3)
      TEMP2 = D02KDS(-ABS(TEMP0-TEMP1))
      IF (GAMMA.EQ.ZER .OR. (GAMMA.NE.ONE .AND. TEMP0.LT.TEMP1))
     *    GO TO 680
      TEMP3 = TEMP0
      TEMP4 = ONE
      TEMP5 = TEMP2
      GO TO 700
  680 TEMP3 = TEMP1
      TEMP4 = TEMP2
      TEMP5 = ONE
  700 TEMP6 = TEMP4*GAMMA + TEMP5*(ONE-GAMMA)
      TEMP7 = TEMP4*VL(4) - TEMP5*V(4)
      IF (TEMP7.EQ.0.0D0) TEMP7 = TENU
      TEMP7 = SIGN(MAX(ABS(TEMP7),TENU*ABS(TEMP6)),TEMP7)
      F(2) = -TEMP7/TEMP6
      ZETA = ABS(F(2)*V(7))
      F(2) = ONE/F(2)
      TEMP8 = -TEMP3 - LOG(ABS(TEMP7)/TWO)
      V3L = V3L + TEMP8 - VL(3)
      V3R = V3R + TEMP8 - V(3)
      SIGNSI = SIGN(ONE,TEMP7)
      HMAX(1,NXP-1) = -D02KDS(-V3BDYL+V3L)*SIGNSI
      HMAX(1,NXP) = D02KDS(V3R-V3BDYR)*SIGNSI
      TEMP2 = TEMP3 + LOG(.75D0*TEMP6*MAX(ZETA,1.D-2))
C     DETOL=IDEALISED VALUE OF V(5) FOR NEXT ITERATION)
      DETOL = TEMP2 + TEMP8
C     BUT INTRODUCE CAUTION
      DETOLL = MIN(DETOLL+TWO,DETOL-V3L)
      DETOLR = MIN(DETOLR+TWO,DETOL-V3R)
C
      IFAIL1 = 0
      IF (MAXIT.EQ.0) IFAIL1 = -1
      ISTAT1 = ISTATE
  720 CONTINUE
C
      F(4) = MAXFUN - MAXFN1
      CALL MONIT(MAXIT,ISTAT1,LAMDA,F)
      IF (IFAIL1.GT.0) GO TO 740
      IF (ZETA.LT.ONE) GO TO 540
  740 GO TO (100,140,240,300) IBACK
  760 GO TO (780,800,800,820,840,780,840,860,780)
     *       IFAIL1
  780 IFAIL1 = 11
      TEMP1 = DBLE(IFAIL1)
      TEMP2 = V(18)
      GO TO 880
  800 IFAIL1 = 8
      GO TO 880
  820 IFAIL1 = 7
      GO TO 880
  840 IFAIL1 = 6
      GO TO 880
  860 IFAIL1 = 3
  880 ISTAT1 = -IFAIL1
      F(9) = JINT
      F(10) = X
      F(11) = V(1)
      F(12) = V(2)
      F(13) = V(3)
      F(14) = V(18)
      F(15) = V(11)
      GO TO 720
C     BDYVAL ERROR AT XL
  900 F(9) = ZER
      F(10) = XL
      GO TO 940
C     BDYVAL ERROR AT XR
  920 F(9) = DBLE(NXP-2)
      F(10) = XR
  940 IFAIL1 = 2
      ISTAT1 = -IFAIL1
      GO TO 720
C
C     MAXIT ERROR
  960 IFAIL1 = 4
      GO TO 740
      END
      SUBROUTINE D02KEZ(X,XEND,V,COEFFN,MAXFUN,HINFO,REPORT,IFAIL)
C     MARK 8 RELEASE. NAG COPYRIGHT 1979.
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     MARK 12A REVISED. IER-501 (AUG 1986).
C     COEFFN, REPORT
C     .. Scalar Arguments ..
      DOUBLE PRECISION  X, XEND
      INTEGER           IFAIL, MAXFUN
C     .. Array Arguments ..
      DOUBLE PRECISION  HINFO(2), V(21)
C     .. Subroutine Arguments ..
      EXTERNAL          COEFFN, REPORT
C     .. Scalars in Common ..
      DOUBLE PRECISION  BP, LAMDA, MINSC, ONE, PI, PSIGN, TWO, ZER
      INTEGER           JINT
C     .. Arrays in Common ..
      DOUBLE PRECISION  YL(3), YR(3)
C     .. Local Scalars ..
      DOUBLE PRECISION  BNEW, CACC1, CACC2, DQDL, EPSDE, FAC1, FAC2, P,
     *                  Q, SY, SYOLD, TEMP, V2OLD, X1, XMAXFN, XOLD
      INTEGER           IFAIL1
C     .. Local Arrays ..
      DOUBLE PRECISION  CIN(6), COMM(5), CONST(3), W(3,7)
C     .. External Functions ..
      DOUBLE PRECISION  D02KDS, D02KDU
      EXTERNAL          D02KDS, D02KDU
C     .. External Subroutines ..
      EXTERNAL          D02KDV, D02KDW, D02KDY
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, MAX, SIGN, LOG, COS, DBLE, SIN
C     .. Common blocks ..
      COMMON            /AD02KD/ZER, ONE, TWO, PI, LAMDA, PSIGN, MINSC,
     *                  BP, YL, YR, JINT
C     .. Data statements ..
      DATA              COMM(1), COMM(2), COMM(3), COMM(4),
     *                  COMM(5)/0.D0, 0.D0, 0.D0, 1.D0, 0.D0/, CONST(1),
     *                  CONST(2), CONST(3)/0.D0, 0.D0, 0.D0/,
     *                  CIN(2)/1.D0/, CACC1/4.D0/
C     .. Executable Statements ..
      CACC2 = LOG(CACC1)
      CIN(1) = ONE
      CIN(3) = ZER
      CIN(4) = HINFO(1)
      CIN(5) = HINFO(2)
C     IF MAXFUN.LE.0 ADD ON BIG QUANTITY
      XMAXFN = ZER
      IF (MAXFUN.LE.0) XMAXFN = 32767.D0
      COMM(1) = DBLE(MAXFUN) + XMAXFN
C     GIVE CIN(6) A STARTING VALUE AS IT IS USED IN CALCULATING BP
      CIN(6) = CIN(5)
      IF (CIN(6).EQ.ZER) CIN(6) = (XEND-X)*0.125D0
C     RE-SCALE PRUFER VARIABLES IF NECESSARY
      BNEW = D02KDU(X,COEFFN)
      IF (BNEW.NE.V(1)) CALL D02KDV(V,BNEW,IFAIL)
C     INITIAL EVALUATION OF BP
      IF (X.EQ.XEND) GO TO 120
      X1 = X + CIN(6)
C     INITIAL EVALUATION OF INTEGRAND FOR SENSITIVITY INTEGRAL
      CALL COEFFN(P,Q,DQDL,X,LAMDA,JINT)
      TEMP = MAX(D02KDU(X1,COEFFN)-V(1),-ABS(Q*CIN(6)))
      BP = TEMP/CIN(6)
      SY = D02KDS(V(3)-V(5))*DQDL/V(1)
C     STORE INITIAL P(X) IN COMMON SO AUX CAN CHECK SIGN CHANGE
      IF (P.EQ.0.0D0) GO TO 180
      PSIGN = SIGN(ONE,P)
C
C     MAIN LOOP - ADVANCES INTEGRATION ONE STEP
C
   20 CONTINUE
C     SET SOFT FAIL
      IFAIL1 = 1
C     SET LOCAL ERROR TOLERANCE FOR THIS STEP
      EPSDE = D02KDS(V(5)-V(3))
      EPSDE = EPSDE/(ONE+100.D0*EPSDE)
C     STORE OLD VALUES FOR SENSITIVITY INTEGRAL
      XOLD = X
      SYOLD = SY
      V2OLD = V(2)
C
      CALL D02KDY(X,XEND,3,V,CIN,EPSDE,D02KDW,COMM,CONST,V(8)
     *            ,W,3,7,COEFFN,COEFFN,HINFO,1,IFAIL1)
C
C     PSIGN SET TO 0 SHOWS P(X) ZERO OR CHANGED SIGN
      IF (PSIGN.EQ.ZER) GO TO 180
      IF (IFAIL1.NE.0) GO TO 160
C
C     ADVANCE ESTIMATE OF SENSITIVITY INTEGRAL
      CALL COEFFN(P,Q,DQDL,X,LAMDA,JINT)
      SY = D02KDS(V(3)-V(5))*DQDL/V(1)
      TEMP = V(2) - V2OLD
      IF (ABS(TEMP).LE.0.75D0) GO TO 40
      FAC1 = ONE - (SIN(V(2))-SIN(V2OLD))/TEMP
      FAC2 = FAC1
      GO TO 60
   40 FAC1 = ONE - COS(V2OLD)
      FAC2 = ONE - COS(V(2))
   60 CONTINUE
      V(4) = V(4) + (X-XOLD)*(SYOLD*FAC1+SY*FAC2)/TWO
C     CHECK IF ACCURACY UNNECESSARILY HIGH
   80 IF (ABS(V(4)*V(7)).LE.CACC1) GO TO 100
      V(4) = V(4)/CACC1
      SY = SY/CACC1
      V(5) = V(5) + CACC2
      GO TO 80
  100 CONTINUE
C
C     PREPARE FOR NEXT STEP, AND LOOP BACK
      IF (V(6).NE.ZER) CALL REPORT(X,V,JINT)
      IF (X.NE.X1) CALL D02KDV(V,D02KDU(X,COEFFN),IFAIL)
      IF (X.EQ.XEND) GO TO 120
      IF (CIN(1).NE.5.D0) GO TO 200
      X1 = X + CIN(6)
      IF (ABS(XEND-X).LT.ABS(CIN(6))) GO TO 110
      TEMP = MAX(D02KDU(X1,COEFFN)-V(1),-ABS(Q*CIN(6)))
      BP = TEMP/CIN(6)
  110 CONTINUE
      CALL D02KDW(3,X,V,W(1,1),COEFFN,COEFFN,1,HINFO)
      COMM(1) = COMM(1) - ONE
C     AND LOOP BACK
      IF (COMM(1).GT.ZER) GO TO 20
      IFAIL1 = 5
      GO TO 160
C     BEFORE EXIT, SET OUTPUT VALUES ETC.
  120 IFAIL = 0
  140 MAXFUN = COMM(1) - XMAXFN
      HINFO(2) = CIN(5)
      RETURN
C
C     FAILURE EXITS
  160 IFAIL = IFAIL1
      V(18) = EPSDE
      GO TO 140
  180 IFAIL = 8
      V(18) = EPSDE
      GO TO 140
  200 IFAIL = 9
      V(18) = CIN(1)
      GO TO 140
      END
      DOUBLE PRECISION FUNCTION D02PAY(X,H,XEND)
C     MARK 14 RELEASE.  NAG COPYRIGHT 1989.
C     .. Scalar Arguments ..
      DOUBLE PRECISION                 H, X, XEND
C     .. Local Scalars ..
      DOUBLE PRECISION                 TEMP
C     .. External Functions ..
      DOUBLE PRECISION                 D02PAZ, X02AJF
      EXTERNAL                         D02PAZ, X02AJF
C     .. Intrinsic Functions ..
      INTRINSIC                        ABS, MAX, SIGN
C     .. Executable Statements ..
      TEMP = D02PAZ(X+H) - XEND
      IF (ABS(TEMP).LE.2.0D0*X02AJF()*MAX(ABS(X),ABS(XEND)))
     *    TEMP = X02AJF()*SIGN(1.0D0,XEND-X)
      D02PAY = TEMP
      RETURN
      END
      DOUBLE PRECISION FUNCTION D02PAZ(X)
C     MARK 14 RELEASE.  NAG COPYRIGHT 1989.
C     .. Scalar Arguments ..
      DOUBLE PRECISION                 X
C     .. Executable Statements ..
      D02PAZ = X
      RETURN
      END
      INTEGER FUNCTION P01ABF(IFAIL,IERROR,SRNAME,NREC,REC)
C     MARK 11.5(F77) RELEASE. NAG COPYRIGHT 1986.
C     MARK 13 REVISED. IER-621 (APR 1988).
C     MARK 13B REVISED. IER-668 (AUG 1988).
C
C     P01ABF is the error-handling routine for the NAG Library.
C
C     P01ABF either returns the value of IERROR through the routine
C     name (soft failure), or terminates execution of the program
C     (hard failure). Diagnostic messages may be output.
C
C     If IERROR = 0 (successful exit from the calling routine),
C     the value 0 is returned through the routine name, and no
C     message is output
C
C     If IERROR is non-zero (abnormal exit from the calling routine),
C     the action taken depends on the value of IFAIL.
C
C     IFAIL =  1: soft failure, silent exit (i.e. no messages are
C                 output)
C     IFAIL = -1: soft failure, noisy exit (i.e. messages are output)
C     IFAIL =-13: soft failure, noisy exit but standard messages from
C                 P01ABF are suppressed
C     IFAIL =  0: hard failure, noisy exit
C
C     For compatibility with certain routines included before Mark 12
C     P01ABF also allows an alternative specification of IFAIL in which
C     it is regarded as a decimal integer with least significant digits
C     cba. Then
C
C     a = 0: hard failure  a = 1: soft failure
C     b = 0: silent exit   b = 1: noisy exit
C
C     except that hard failure now always implies a noisy exit.
C
C     S.Hammarling, M.P.Hooper and J.J.du Croz, NAG Central Office.
C
C     .. Scalar Arguments ..
      INTEGER                 IERROR, IFAIL, NREC
      CHARACTER*(*)           SRNAME
C     .. Array Arguments ..
      CHARACTER*(*)           REC(*)
C     .. Local Scalars ..
      INTEGER                 I, NERR
      CHARACTER*72            MESS
C     .. External Subroutines ..
      EXTERNAL                P01ABZ, X04AAF, X04BAF
C     .. Intrinsic Functions ..
      INTRINSIC               ABS, MOD
C     .. Executable Statements ..
      IF (IERROR.NE.0) THEN
C        Abnormal exit from calling routine
         IF (IFAIL.EQ.-1 .OR. IFAIL.EQ.0 .OR. IFAIL.EQ.-13 .OR.
     *       (IFAIL.GT.0 .AND. MOD(IFAIL/10,10).NE.0)) THEN
C           Noisy exit
            CALL X04AAF(0,NERR)
            DO 20 I = 1, NREC
               CALL X04BAF(NERR,REC(I))
   20       CONTINUE
            IF (IFAIL.NE.-13) THEN
               WRITE (MESS,FMT=99999) SRNAME, IERROR
               CALL X04BAF(NERR,MESS)
               IF (ABS(MOD(IFAIL,10)).NE.1) THEN
C                 Hard failure
                  CALL X04BAF(NERR,
     *                     ' ** NAG hard failure - execution terminated'
     *                        )
                  CALL P01ABZ
               ELSE
C                 Soft failure
                  CALL X04BAF(NERR,
     *                        ' ** NAG soft failure - control returned')
               END IF
            END IF
         END IF
      END IF
      P01ABF = IERROR
      RETURN
C
99999 FORMAT (' ** ABNORMAL EXIT from NAG Library routine ',A,': IFAIL',
     *  ' =',I6)
      END
      SUBROUTINE P01ABZ
C     MARK 11.5(F77) RELEASE. NAG COPYRIGHT 1986.
C
C     Terminates execution when a hard failure occurs.
C
C     ******************** IMPLEMENTATION NOTE ********************
C     The following STOP statement may be replaced by a call to an
C     implementation-dependent routine to display a message and/or
C     to abort the program.
C     *************************************************************
C     .. Executable Statements ..
      STOP
      END
      DOUBLE PRECISION FUNCTION X01AAF(X)
C     MARK 8 RE-ISSUE. NAG COPYRIGHT 1980.
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C
C     RETURNS THE VALUE OF THE MATHEMATICAL CONSTANT PI.
C
C     X IS A DUMMY ARGUMENT
C
C     IT MAY BE NECESSARY TO ROUND THE REAL CONSTANT IN THE
C     ASSIGNMENT STATEMENT TO A SMALLER NUMBER OF SIGNIFICANT
C     DIGITS IN ORDER TO AVOID COMPILATION PROBLEMS.  IF SO, THEN
C     THE NUMBER OF DIGITS RETAINED SHOULD NOT BE LESS THAN
C     .     2 + INT(FLOAT(IT)*ALOG10(IB))
C     WHERE  IB  IS THE BASE FOR THE REPRESENTATION OF FLOATING-
C     .             POINT NUMBERS
C     . AND  IT  IS THE NUMBER OF IB-ARY DIGITS IN THE MANTISSA OF
C     .             A FLOATING-POINT NUMBER.
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION                 X
C     .. Executable Statements ..
      X01AAF = 3.14159265358979323846264338328D0
      RETURN
      END
      DOUBLE PRECISION FUNCTION X02AJF()
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C
C     RETURNS  (1/2)*B**(1-P)  IF ROUNDS IS .TRUE.
C     RETURNS  B**(1-P)  OTHERWISE
C
C     .. Executable Statements ..
      x02ajf = 1.110223024625157D-16
      RETURN
      END
      DOUBLE PRECISION FUNCTION X02AKF()
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C
C     RETURNS  B**(EMIN-1)  (THE SMALLEST POSITIVE MODEL NUMBER)
C
C     .. Executable Statements ..
      x02akf = 1d-300
c      x02akf = 2.9387358770557188D-39
      RETURN
      END
      DOUBLE PRECISION FUNCTION X02AMF()
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C
C     RETURNS THE 'SAFE RANGE' PARAMETER
C     I.E. THE SMALLEST POSITIVE MODEL NUMBER Z SUCH THAT
C     FOR ANY X WHICH SATISFIES X.GE.Z AND X.LE.1/Z
C     THE FOLLOWING CAN BE COMPUTED WITHOUT OVERFLOW, UNDERFLOW OR OTHER
C     ERROR
C
C        -X
C        1.0/X
C        SQRT(X)
C        LOG(X)
C        EXP(LOG(X))
C        Y**(LOG(X)/LOG(Y)) FOR ANY Y
C
C     .. Executable Statements ..
      x02amf = 1d-300
c      x02amf = 2.9387358770557188D-39
      RETURN
      END
      SUBROUTINE X04AAF(I,NERR)
C     MARK 7 RELEASE. NAG COPYRIGHT 1978
C     MARK 7C REVISED IER-190 (MAY 1979)
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     MARK 14 REVISED. IER-829 (DEC 1989).
C     IF I = 0, SETS NERR TO CURRENT ERROR MESSAGE UNIT NUMBER
C     (STORED IN NERR1).
C     IF I = 1, CHANGES CURRENT ERROR MESSAGE UNIT NUMBER TO
C     VALUE SPECIFIED BY NERR.
C
C     .. Scalar Arguments ..
      INTEGER           I, NERR
C     .. Local Scalars ..
      INTEGER           NERR1
C     .. Save statement ..
      SAVE              NERR1
C     .. Data statements ..
      DATA              NERR1/6/
C     .. Executable Statements ..
      IF (I.EQ.0) NERR = NERR1
      IF (I.EQ.1) NERR1 = NERR
      RETURN
      END
      SUBROUTINE X04BAF(NOUT,REC)
C     MARK 11.5(F77) RELEASE. NAG COPYRIGHT 1986.
C
C     X04BAF writes the contents of REC to the unit defined by NOUT.
C
C     Trailing blanks are not output, except that if REC is entirely
C     blank, a single blank character is output.
C     If NOUT.lt.0, i.e. if NOUT is not a valid Fortran unit identifier,
C     then no output occurs.
C
C     .. Scalar Arguments ..
      INTEGER           NOUT
      CHARACTER*(*)     REC
C     .. Local Scalars ..
      INTEGER           I
C     .. Intrinsic Functions ..
      INTRINSIC         LEN
C     .. Executable Statements ..
      IF (NOUT.GE.0) THEN
C        Remove trailing blanks
         DO 20 I = LEN(REC), 2, -1
            IF (REC(I:I).NE.' ') GO TO 40
   20    CONTINUE
C        Write record to external file
   40    WRITE (NOUT,FMT=99999) REC(1:I)
      END IF
      RETURN
C
99999 FORMAT (A)
      END
