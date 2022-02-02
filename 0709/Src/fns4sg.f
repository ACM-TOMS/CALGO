      SUBROUTINE ZZFNS4 ( N, X, F, G, WORK, SIZE, FONLY, GONLY,
     -                    FIRST, FARG, FUNCNO, STATUS )

C## A R G U M E N T S:
                 INTEGER            N, SIZE, FUNCNO, STATUS(7)

                 LOGICAL            FIRST,   FONLY,  GONLY

                 REAL               F, X(N), G(N), WORK(SIZE), FARG(*)
C!!!!            DOUBLE PRECISION   F, X(N), G(N), WORK(SIZE), FARG(*)

C## S T A T U S:
C               SINGLE/DOUBLE CONVERSION: NEEDED (SEE CONVRT).
C
C               IGNORE LINES BEGINNING WITH  "C!!!!" .
C
C               THIS VERSION IS IN   S I N G L E   PRECISION.
C!!!!           THIS VERSION IS IN   D O U B L E   PRECISION.
C
C               SYSTEM  DEPENDENCE:                      NONE.
C
C>RCS $HEADER: FNS4.F,V 1.2 91/12/31 14:38:18 BUCKLEY EXP $
C>RCS $LOG:     FNS4.F,V $
C>RCSREVISION 1.2  91/12/31  14:38:18  BUCKLEY
C>RCSFINAL SUBMISSION TO TOMS
C>RCS
C>RCSREVISION 1.1  91/11/20  10:52:53  BUCKLEY
C>RCSFINAL SUBMISSION TO TOMS
C>RCS
C>RCSREVISION 1.0  90/07/31  13:01:57  BUCKLEY
C>RCSINITIAL REVISION
C>RCS
C
C## D E S C R I P T I O N:
C
C         THIS SUBROUTINE EVALUATES ONE OF THE STANDARD TEST FUNCTIONS.
C     THE TEST FUNCTIONS ARE DIVIDED AMOUNG FOUR ROUTINES BECAUSE A
C     SINGLE ROUTINE WOULD BE TOO LARGE FOR SOME COMPILERS.
C     THE ARGUMENTS IN THE CALLING SEQUENCE HAVE PRECISELY THE SAME
C     MEANING AS IN THE ROUTINE ZZEVAL.
C
C         THE VALUE OF THE INTEGER PARAMETER FUNCNO SPECIFIES
C     WHICH OF THE TEST FUNCTIONS IS TO BE  USED; THE FUNCTION
C     IS CHOSEN USING A COMPUTED GOTO.
C
C         THE PARAMETERS FONLY AND GONLY SPECIFY FUNCTION AND
C     GRADIENT EVALUATIONS.  THE PARAMETER FIRST SPECIFIES CODE TO BE
C     EVALUATED ONLY ON THE FIRST CALL TO A PARTICULAR FUNCTION.
C     THE PARAMETER STATUS STORES THE RETURN CODES.
C
C## E N T R Y   P O I N T S:  ZZFNS4    THE NATURAL ENTRY POINT.
C
C## S U B R O U T I N E S:
C
C     PREDEFINED FUNCTIONS : SIN, COS, TAN, ACOS, ATAN, ABS, MAX, NINT
C                            EXP, LOG, MIN, MOD, SIGN, SQRT, REAL(DBLE)
C
C     STATEMENT FUNCTION:      RD
C
C## P A R A M E T E R S:

      REAL              ZERO,       ONE,       TWO,       THREE
C!!!! DOUBLE PRECISION  ZERO,       ONE,       TWO,       THREE
      PARAMETER (       ZERO = 0D0, ONE = 1D0, TWO = 2D0, THREE = 3D0)

      REAL              FOUR,       FIVE,      SIX,       SEVEN
C!!!! DOUBLE PRECISION  FOUR,       FIVE,      SIX,       SEVEN
      PARAMETER (       FOUR = 4D0, FIVE = 5D0,SIX = 6D0, SEVEN = 7D0)

      REAL              EIGHT,         NINE,          TEN
C!!!! DOUBLE PRECISION  EIGHT,         NINE,          TEN
      PARAMETER (       EIGHT = 8D0,   NINE = 9D0,    TEN = 10D0     )

      REAL              TENTH,         FIFTH,        HALF
C!!!! DOUBLE PRECISION  TENTH,         FIFTH,        HALF
      PARAMETER (       TENTH = .1D0,  FIFTH = .2D0, HALF = .5D0      )

      REAL              RPT9,        RPT8,        RD29
C!!!! DOUBLE PRECISION  RPT9,        RPT8,        RD29
      PARAMETER (       RPT9 = .9D0, RPT8 = .8D0, RD29 = 1D0/29D0 )
      REAL              R11,        R12,        R13,       R14
C!!!! DOUBLE PRECISION  R11,        R12,        R13,       R14
      PARAMETER (       R11 = 11D0, R12 = 12D0, R13 = 13D0,R14 = 14D0)

      REAL              R15,        R16,        R17,       R18
C!!!! DOUBLE PRECISION  R15,        R16,        R17,       R18
      PARAMETER (       R15 = 15D0, R16 = 16D0, R17 = 17D0,R18 = 18D0)

      REAL              R19,        R20,        R25,       R29
C!!!! DOUBLE PRECISION  R19,        R20,        R25,       R29
      PARAMETER (       R19 = 19D0, R20 = 20D0, R25 = 25D0,R29 = 29D0)

      REAL              R32,        R36,        R40,       R42
C!!!! DOUBLE PRECISION  R32,        R36,        R40,       R42
      PARAMETER (       R32 = 32D0, R36 = 36D0, R40 = 40D0,R42 = 42D0)

      REAL              R45,        R49
C!!!! DOUBLE PRECISION  R45,        R49
      PARAMETER (       R45 = 45D0, R49 = 49D0 )

      REAL              R50,        R56,        R84,       R90
C!!!! DOUBLE PRECISION  R50,        R56,        R84,       R90
      PARAMETER (       R50 = 50D0, R56 = 56D0, R84 = 84D0,R90 = 90D0)

      REAL              R100,            R180,           R200
C!!!! DOUBLE PRECISION  R100,            R180,           R200
      PARAMETER (       R100 = 100D0,    R180 = 180D0,   R200 = 200D0 )

      REAL              R256,            R360,           R400
C!!!! DOUBLE PRECISION  R256,            R360,           R400
      PARAMETER (       R256 = 256D0,    R360 = 360D0,   R400 = 400D0 )

      REAL              R600,            R681,           R991
C!!!! DOUBLE PRECISION  R600,            R681,           R991
      PARAMETER (       R600 = 600D0,    R681 = 681D0,   R991 = 991D0 )

      REAL              R1162,                 R2324
C!!!! DOUBLE PRECISION  R1162,                 R2324
      PARAMETER (       R1162 = 1162D0,        R2324 = 2324D0         )

      REAL              R10000,                R40000
C!!!! DOUBLE PRECISION  R10000,                R40000
      PARAMETER (       R10000 = 10000D0,      R40000 = 40000D0       )
      REAL              R1PD6,                 R2PDM6
C!!!! DOUBLE PRECISION  R1PD6,                 R2PDM6
      PARAMETER (       R1PD6 = 1D6,           R2PDM6 = 2D-6         )

      REAL              RP04,        RP01,          R1PZ1
C!!!! DOUBLE PRECISION  RP04,        RP01,          R1PZ1
      PARAMETER (       RP04 = 4D-2, RP01 = .01D0,  R1PZ1 = 1.0001D0 )

      REAL              R1P2,         R7P5,         RP1136
C!!!! DOUBLE PRECISION  R1P2,         R7P5,         RP1136
      PARAMETER (       R1P2 = 1.2D0, R7P5 = 7.5D0, RP1136 = 0.1136D0 )

      REAL              R1P5,         R2P5,         R2P625
C!!!! DOUBLE PRECISION  R1P5,         R2P5,         R2P625
      PARAMETER (       R1P5 = 1.5D0, R2P5 = 2.5D0, R2P625 = 2.625D0 )

      REAL              R10P1,         R19P8,         R20P2
C!!!! DOUBLE PRECISION  R10P1,         R19P8,         R20P2
      PARAMETER (       R10P1 = 10.1D0,R19P8 = 19.8D0,R20P2 = 20.2D0 )

      REAL              R2D3,          R4D3,          R7D3
C!!!! DOUBLE PRECISION  R2D3,          R4D3,          R7D3
      PARAMETER (       R2D3 = 2D0/3D0,R4D3 = 4D0/3D0,R7D3 = 7D0/3D0 )

      REAL              R2P25
C!!!! DOUBLE PRECISION  R2P25
      PARAMETER (       R2P25 = 2.25D0 )

C## L O C A L   D E C L:

      INTEGER     OK,  ABORT,  LIMIT,  NOF,  NOG,  NOFG

      INTEGER           F1, I, J, K, J0, J1, J2, J3,  JLO, JHI
      INTEGER           IR, I3, IL, N1, RET, P, M, L, IDUMMY
      INTEGER           RHO1, RHO2, K1, K2, K3, K4, I1, I2

      LOGICAL           GFIRST,   DONE,   PROB13

      REAL              ZZMPAR, HUGE
C!!!! DOUBLE PRECISION  ZZMPAR, HUGE

C--------- VARIABLES FOR THE TEST FUNCTIONS.

      REAL              X1,     X2,     X3,     X4,   S1
C!!!! DOUBLE PRECISION  X1,     X2,     X3,     X4,   S1

      REAL              W1,     W2,     W3,    W4,    W5,    W6
C!!!! DOUBLE PRECISION  W1,     W2,     W3,    W4,    W5,    W6

      REAL              W7,     W8,     W9,    W10
C!!!! DOUBLE PRECISION  W7,     W8,     W9,    W10

      REAL               R,      S,      T,    R1, BIGGST, SMLLST
C!!!! DOUBLE PRECISION   R,      S,      T,    R1, BIGGST, SMLLST

      REAL              R2,     R3,     RI,    TI
C!!!! DOUBLE PRECISION  R2,     R3,     RI,    TI

      REAL              XI,     YI,    PI
C!!!! DOUBLE PRECISION  XI,     YI,    PI

      REAL             XP1,    XM1,    R2P,      RD, TPI, TPIS
C!!!! DOUBLE PRECISION XP1,    XM1,    R2P,      RD, TPI, TPIS

      REAL              KAP1, KAP2, KAP3,  H,  TK, ALPHA
C!!!! DOUBLE PRECISION  KAP1, KAP2, KAP3,  H,  TK, ALPHA


C## S A V E:

      SAVE  GFIRST, PI,  TPI, TPIS, R2P, BIGGST, SMLLST
      SAVE  HUGE, PROB13, DONE, M
      SAVE   OK, ABORT, LIMIT, NOF, NOG, NOFG

C## E Q U I V A L E N C E S: NONE ARE DEFINED.
C## C O M M O N:             NONE IS DEFINED.
C## D A T A:

      DATA  GFIRST/.TRUE./

C##                                                E X E C U T I O N
C##                                                E X E C U T I O N
C---------                     FUNCTION DEFINITION
      RD (IDUMMY)  = REAL (IDUMMY)
C!!!! RD (IDUMMY)  = DBLE (IDUMMY)

C--------- SOME ONE TIME ONLY CONSTANTS.
      IF ( GFIRST ) THEN
         PI   = ACOS(-ONE)
         TPI  = TWO * PI
         TPIS = TPI * PI
         R2P  = ONE / TPI
         HUGE   = ZZMPAR(3)/TEN
         SMLLST = LOG(ZZMPAR(2)*TEN)
         BIGGST = LOG(HUGE)
         OK     = STATUS(1)
         ABORT  = STATUS(2)
         LIMIT  = STATUS(3)
         NOF    = STATUS(4)
         NOG    = STATUS(5)
         NOFG   = STATUS(6)
      ENDIF

C---------                 SET LOGICAL FLAGS AND SELECT FUNCTION.
      RET   = OK
      GOTO(
     -     1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000,
     -     1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000,
     -     1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000,
     -     1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000,
     -     1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000,
     -     1000, 1000, 6300, 6400, 6500, 6600, 6700, 6800, 6900, 7000,
     -     7100, 7200, 7300, 7400, 7500, 7600, 7700, 7800, 7900, 8000,
     -     8100, 8200, 8300, 7700, 1000, 1000, 1000, 1000, 1000, 1000
     -    ) FUNCNO
 1000 STATUS(7) = ABORT
      GOTO 91000


C----                       ARGQDN ( N, X, F, G, IFG, NINT(FARG(1)))
C----                       NINT(FARG(1)) IS M
 6300 F1 = NINT ( FARG ( 1 ) )

      IF ( .NOT. GONLY ) THEN
         F = ZERO
      ENDIF

      S  = ZERO
      W4 = ZERO
      IF ( F1 .GE. N ) THEN
         DO 6310 I = 1,N
            S  = S + X(I)
 6310    CONTINUE

         W2 = -TWO*S/FARG(1) - ONE
         DO 6320 I = 1,N

            RI = X(I) + W2

            IF ( .NOT. GONLY ) THEN
               F = F + RI**2
            ENDIF

            IF ( .NOT. FONLY ) THEN
               W4   =  W4 + RI
               G(I) = TWO * RI
            ENDIF

 6320    CONTINUE

         W1 = RD( F1 - N )
         IF ( .NOT. GONLY ) THEN
            F = F + W1*W2**2
         ENDIF

         IF ( .NOT. FONLY ) THEN
            W4 = FOUR * ( W4 + W1*W2 ) / FARG(1)
            DO 6330 I = 1,N
               G(I) = G(I) - W4
 6330       CONTINUE
         ENDIF

      ELSE IF ( .NOT. FONLY ) THEN
         DO 6340 I = 1,N
            G(I) = ZERO
 6340    CONTINUE
      ENDIF

      GOTO 90000

C--------- MINTEST FUNCTION   ARGQDO ( N, X, F, G, IFG, NINT(FARG(1)))
C---------                     NINT(FARG(1)) IS M
 6400 F1 = NINT( FARG(1) )

      IF ( .NOT. GONLY ) THEN
         F = ZERO
      ENDIF

      S  = ZERO
      W1 = ZERO

      IF ( F1 .GE. N ) THEN
         DO 6410 I = 1,N
            S = S + RD(I)*X(I)
 6410    CONTINUE
         DO 6420 I = 1,F1
            W2 = RD(I)
            RI = W2*S - ONE
            IF ( .NOT. GONLY ) THEN
               F = F + RI*RI
            ENDIF
            IF ( .NOT. FONLY ) THEN
               W1 = W1 + RI*W2
            ENDIF
 6420    CONTINUE
      ENDIF

      IF ( .NOT. FONLY ) THEN
         W1 = TWO * W1
         DO 6430 I = 1,N
            G(I) = RD(I) * W1
 6430    CONTINUE
      ENDIF
      GOTO 90000

C--------- MINTEST FUNCTION   ARGQDZ ( N, X, F, G, IFG, NINT(FARG(1)))
C---------                     NINT(FARG(1)) IS M
 6500 F1 = NINT( FARG(1) )

      IF ( .NOT. GONLY ) THEN
         F = TWO
      ENDIF

      S  = ZERO
      W1 = ZERO

      IF ( F1 .GE. N ) THEN

         DO 6510 I = 2,N-1
            S = S + RD(I)*X(I)
 6510    CONTINUE

         DO 6520 I = 2,F1-1
            W2 = RD(I-1)
            RI = W2 * S - ONE

            IF ( .NOT. GONLY ) THEN
               F = F + RI*RI
            ENDIF

            IF ( .NOT. FONLY ) THEN
               W1 = W1 + RI*W2
            ENDIF

 6520    CONTINUE

      ENDIF

      IF ( .NOT. FONLY ) THEN
         W1 = TWO * W1
         G(1) = ZERO
         DO 6530 I = 2,N-1
            G(I) = RD(I) * W1
 6530    CONTINUE
         G(N) = ZERO
      ENDIF

      GOTO 90000


C--------- MINTEST FUNCTION     MOREBV


 6600 IF ( .NOT. GONLY ) THEN
         F = ZERO
      ENDIF

      W1 = ONE/(RD(N)+ONE)
      W2 = W1*W1 / TWO
      W3 = THREE*W2
      XI = X(1)
      W4 = XI
      W8 = ZERO

      DO 6610 I = 1,N

         TI = RD(I)*W1

         IF ( I .LT. N ) THEN
            W10 = X(I+1)
         ELSE
            W10 = ZERO
         ENDIF

         W6 = W10 - XI
         W7 = XI + TI + ONE
         RI = W4 - W6 + W2*W7*W7*W7

         IF ( .NOT. GONLY ) THEN
            F = F + RI*RI
         ENDIF

         IF ( .NOT. FONLY ) THEN

            IF ( I .NE. 1 ) THEN
               G(I-1) = TWO*(-W5 - RI + W8*W9)
            ENDIF

            W9 = TWO + W3 * W7*W7
            W5 = W8
            W8 = RI

         ENDIF

         XI  = W10
         W4  = W6

 6610 CONTINUE

      IF ( .NOT. FONLY ) THEN
         G(N) = TWO*(-W5 + RI*W9)
      ENDIF
      GOTO 90000

C---                                 BROY7D
 6700 IF ( .NOT. GONLY ) THEN
         F = ZERO
      ENDIF

      XM1 = ZERO
      XI  = X(1)
      W4  = ZERO

      DO 6710 I = 1,60
         IF (I .LT. 60) THEN
            XP1 = X(I+1)
         ELSE
            XP1 = ZERO
         ENDIF

         YI = XM1 - (THREE-XI/TWO)*XI + TWO*XP1 - ONE

         IF ( .NOT. GONLY ) THEN
            F = F + (ABS(YI))**R7D3
         ENDIF

         IF (I .LE. 30) THEN
            W2 = XI + X(I+30)
            IF ( .NOT. GONLY ) THEN
               F = F + (ABS(W2))**R7D3
            ENDIF
         ELSE
            W2 = X(I-30) + XI
         ENDIF

         IF ( .NOT. FONLY ) THEN
            TI = (ABS(YI))**R4D3
            TI = SIGN(TI,YI)

            IF ( I .GT. 1 ) THEN
               G(I-1) = R7D3*(W1 + TI)
            ENDIF

            W3 = (ABS(W2))**R4D3
            W3 = SIGN(W3,W2)
            W1  = TWO*W4 + TI*(XI - THREE) + W3
            W4  = TI

         ENDIF

         XM1 = XI
         XI  = XP1

 6710 CONTINUE

      IF ( .NOT. FONLY ) THEN
         G(60) = R7D3 * W1
      ENDIF

      GOTO 90000


C--------- MINTEST FUNCTION     BRKMCC


 6800 X1 = X(1)
      X2 = X(2)

      W1 = X1 - TWO
      W2 = X2 - ONE
      W3 = X1 - TWO*X2 + ONE
      W4 = -X1*X1/FOUR - X2*X2 + ONE

      IF ( .NOT. GONLY ) THEN
         F = W1*W1 + W2*W2 + RP04/W4 + FIVE*W3*W3
      ENDIF

      IF ( .NOT. FONLY ) THEN
         W4   = W4*W4
         G(1) = TWO*W1 +     X1/(R50*W4) + TEN*W3
         G(2) = TWO*W2 + TWO*X2/(R25*W4) - R20*W3
      ENDIF

      GOTO 90000


C--------- MINTEST FUNCTION     HIMM29


 6900 X1 = X(1)
      X2 = X(2)
      W3 = X1*X1
      W1 = W3 + R12*X2 - ONE
      W2 = R49*W3 + R49*X2*X2 + R84*X1 + R2324*X2 - R681

      IF ( .NOT. GONLY ) THEN
         F = W1*W1 + W2*W2
      ENDIF

      IF ( .NOT. FONLY ) THEN
         G(1) = FOUR * (   X1*W1 + ( R49*X1 + R42   ) * W2  )
         G(2) = FOUR * (  SIX*W1 + ( R49*X2 + R1162 ) * W2  )
      ENDIF

      GOTO 90000


C--------- MINTEST FUNCTION     HIMM33


 7000 X1 = X(1)
      X2 = X(2)
      ALPHA = FARG(1)
      IF ( -(X1+X2) .LE. BIGGST ) THEN
         W1 = EXP(-(X1+X2))
      ELSE
         RET = NOFG
         GOTO 90000
      ENDIF
      W2 = TWO*X1**2 + THREE*X2**2

      IF ( .NOT. GONLY ) THEN
         F = ALPHA*W1*W2
      ENDIF

      IF ( .NOT. FONLY ) THEN
         G(1) = ALPHA * W1 * (FOUR*X1 - W2)
         G(2) = ALPHA * W1 * (SIX*X2  - W2)
      ENDIF
      IF ( ABS(F) .LE. 100. ) THEN
         WRITE(6,*) 'F     = ', F
         WRITE(6,*) 'G(1)  = ', G(1), '   G(2) = ', G(2)
C         WRITE(6,*) 'G(2)  = ', G(2)
         WRITE(6,*) 'X1    = ', X1, '   X2 = ', X2
C         WRITE(6,*) 'X2    = ', X2
         WRITE(6,*) ' '
      END IF

      GOTO 90000


C---------                  HIMM30
 7100 X1 = X(1)
      X2 = X(2)
      X3 = X(3)
      W1 = X1 + X2
      R1 = X3 - W1*W1/FOUR
      R2 = ONE - X1
      R3 = ONE - X2
      W2 = R100*R1

      IF ( .NOT. GONLY ) THEN
         F = W2*R1 + R2*R2 + R3*R3
      ENDIF

      IF ( .NOT. FONLY ) THEN
         G(1) = -W2*W1 - TWO*R2
         G(2) = -W2*W1 - TWO*R3
         G(3) =  W2*TWO
      ENDIF
      GOTO 90000

C----                       GOTTFR
 7200           CONTINUE
      W1 = X(1) - RP1136*(X(1)+THREE*X(2))*(ONE-X(1))
      W2 = X(2) + R7P5*  (TWO*X(1)-X(2))  *(ONE-X(2))
      IF ( .NOT. GONLY ) THEN
         F = W1**2 + W2**2
      ENDIF
      IF ( .NOT. FONLY ) THEN
         G(1) = TWO*( W1*(ONE-RP1136*(-TWO*X(1)+1-THREE*X(2)))
     -               +W2*(R15*(ONE-X(2))))
         G(2) = TWO*( W1*(-RP1136*THREE*(ONE-X(1)))
     -               +W2*(ONE+R7P5*(TWO*X(2)-ONE-TWO*X(1))))
      ENDIF
      GOTO 90000

C---------                    BRYBND ( N, X, F, G, IFG, WORK, SIZE )
 7300 IF ( SIZE .LT. N ) THEN
         F = ZERO
         DO 7308 I = 1, N
            G(I) = ZERO
 7308    CONTINUE

      ELSE
         KAP1 =      FARG(1)
         KAP2 =      FARG(2)
         KAP3 =      FARG(3)
         RHO1 = NINT(FARG(4))
         RHO2 = NINT(FARG(5))
         DO 7309 I = 1, N
            WORK(I) = ZERO
            IF ( .NOT. FONLY ) THEN
               G(I) = ZERO
            ENDIF
 7309    CONTINUE

         DO 7313 I = 1, N
            W1 = X(I) * ( KAP1 + KAP2 * X(I)**2 )

            JLO = MAX (I-RHO1, 1)
            JHI = MIN (I+RHO2, N)
            DO 7311 J = JLO, I-1
               WORK(I) = WORK(I) + X(J) + X(J)**2
 7311       CONTINUE
            DO 7312 J = I+1, JHI
               WORK(I) = WORK(I) + X(J) + X(J)**2
 7312       CONTINUE
            WORK(I) = W1 + ONE - KAP3*WORK(I)
 7313    CONTINUE

         IF ( .NOT. GONLY ) THEN
            F = ZERO
         ENDIF

         DO 7320 I = 1, N
            IF ( .NOT. GONLY ) THEN
               F = F + WORK(I) **2
            ENDIF

            IF ( .NOT. FONLY ) THEN
               W2 = TWO * KAP3 *( ONE + TWO * X(I) )
               W1 = TWO * WORK(I) * ( KAP1 + THREE*KAP2* X(I)**2 )

               JLO = MAX (I-RHO2, 1)
               JHI = MIN (I+RHO1, N)
               DO 7314 J = JLO, I-1
                  G(I) = G(I) + WORK(J)
 7314          CONTINUE
               DO 7316 J = I+1, JHI
                  G(I) = G(I) + WORK(J)
 7316          CONTINUE
            G(I) = W1 - W2 * G(I)
            ENDIF
 7320    CONTINUE
      ENDIF
      GO TO 90000

C-----                       CLUSTR
 7400                  CONTINUE
      W1 = SIN(X(1))
      W2 = SIN(X(2))
      W3 = COS(X(1))
      W4 = COS(X(2))
      R1 = (X(1)-X(2)**2)*(X(1)-W2)
      R2 = (W4-X(1))*(X(2)-W3)
      IF ( .NOT. GONLY ) THEN
         F = R1**2 + R2**2
      ENDIF
      IF ( .NOT. FONLY ) THEN
         G(1) = TWO*(R1*(TWO*X(1)-W2-X(2)**2)
     -             + R2*(-X(2)+W3+(W4-X(1))*W1))
         G(2) = TWO*(R1*(-TWO*X(2)*(X(1)-W2)-W4*(X(1)-X(2)**2))
     -             + R2*(-X(1)+W4+(W3-X(2))*W2))
      ENDIF
      GOTO 90000
C------                      ARTRIG(N,X,F,G,IFG,WORK)
 7500 IF ( .NOT. GONLY ) THEN
         F  = ZERO
      ENDIF

      IF ( .NOT. FONLY ) THEN
         DO 7505 I = 1,N
            G(I)  = ZERO
            W3    = ZERO
 7505    CONTINUE
      ENDIF

      IF ( SIZE .LT. 2*N ) THEN
         RET = NOFG
         GOTO 90000
      ELSE
         W1 = ZERO
         DO 7507 I = 1,N
            WORK(I)   = SIN(X(I))
            WORK(N+I) = COS(X(I))
            W1 = W1 + WORK(N+I)
 7507    CONTINUE
      ENDIF

      DO 7580 I = 1,N
         W2 =  N - W1 + I*(ONE-WORK(N+I)) - WORK(I)
         IF ( .NOT. GONLY ) THEN
            F = F + W2**2
         ENDIF
         IF ( .NOT. FONLY ) THEN
            G(I) = G(I) + W2 * ( (ONE+I) * WORK(I)   - WORK(N+I) )
            W3   = W3 + W2
            WORK(N+I) = W2
         ENDIF
 7580 CONTINUE

      IF ( .NOT. FONLY ) THEN
         DO 7590 I = 1, N
            G(I) = G(I) + WORK(I) * ( W3 - WORK(N+I) )
 7590    CONTINUE
         CALL ZZSCAL ( N, TWO, G, 1 )
      ENDIF

      GOTO 90000

C------------------     FUNCTION SQRTMX(N,X,F,G,IFG,FARG(1))
C     FARG(1) MUST BE SET TO A NON-ZERO VALUE FOR
C             THE LIU-NOCEDAL PROBLEM #13 (=0 OR 1) OR #15 (=2).
 7600 CONTINUE
      IF ( FIRST .AND. FARG(1) .NE. 2 ) THEN
         PROB13 = FARG(1) .NE. ZERO
         M      = NINT( SQRT(RD(N)) )
         FIRST = .FALSE.

         IF ( SIZE .GE. N ) THEN
            DO 7610 I = 1,N
               L  = MOD(I-1,M)
               K  =  (I-1)/M
               W1 = ZERO
               DO 7605 J = 1,M
                  P = M*L + J
                  IF ( PROB13 .AND. P .EQ. 2*M+1) THEN
                     W2 = ZERO
                  ELSE
                     W2 = SIN(RD(P**2))
                  ENDIF
                  P = M*(J-1) + K + 1
                  IF ( PROB13 .AND. P .EQ. 2*M+1) THEN
                     W3 = ZERO
                  ELSE
                     W3 = SIN(RD(P**2))
                  ENDIF
                  W1 = W1 + W2*W3
 7605          CONTINUE
               WORK(I) = W1
 7610       CONTINUE
            DONE  = .FALSE.
         ELSE
            DONE  = .FALSE.
         ENDIF
      ENDIF

      IF ( .NOT. GONLY ) THEN
         F = ZERO
      ENDIF

      IF ( .NOT. FONLY ) THEN
         DO 7620 I = 1,N
            G(I) = ZERO
 7620    CONTINUE
      ENDIF

      IF ( FARG(1) .NE. 2 ) THEN
         DO 7640 I = 1,N
            L  = MOD(I-1,M)
            K  = (I-1)/M
            W1 = ZERO
            DO 7625 J = 1,M
               W1 = W1 + X(M*L+J) * X(M*(J-1)+K+1)
 7625       CONTINUE
            IF ( DONE ) THEN
               W4 = WORK(I)
            ELSE
               W4 = ZERO
               DO 7630 J = 1,M
                  P = M*L + J
                  IF ( PROB13 .AND. P .EQ. 2*M+1) THEN
                     W2 = ZERO
                  ELSE
                     W2 = SIN(RD(P**2))
                  ENDIF
                  P = M*(J-1)+K+1
                  IF ( PROB13 .AND. P .EQ. 2*M+1) THEN
                     W3 = ZERO
                  ELSE
                     W3 = SIN(RD(P**2))
                  ENDIF
                  W4 = W4 + W2*W3
 7630          CONTINUE
            ENDIF

            IF ( .NOT. GONLY ) THEN
               F = F + (W4 - W1)**2
            ENDIF

            IF ( .NOT. FONLY ) THEN
               DO 7635 J = 1,M
                  I1 = M*L + J
                  I2 = M*(J-1) + K + 1
                  G(I1) = G(I1) + (W4-W1) * X(I2)
                  G(I2) = G(I2) + (W4-W1) * X(I1)
 7635          CONTINUE
            ENDIF
 7640    CONTINUE
         IF ( .NOT. FONLY ) THEN
            CALL ZZSCAL ( N, -TWO, G, 1 )
         ENDIF
      ELSE
C----                      SPARSE SQUARE ROOT MATRIX
                  CONTINUE
         M  = (N+2)/3
         K1 = 5*M-6
         I = 1
         J = 0

         DO 7690 K=1,K1
C                         COMPUTE (I,J) FOR A GIVEN K
            J = J+1
            IF ( K .EQ. 4 .OR. K .EQ. 8 .OR. K .EQ. 13) THEN
               I = I+1
               J = 1
            ENDIF
 7660       IF (J .EQ. 1) THEN
               J1 = I-3
               K3 = J1
               IF (J1 .LE. 0) J1 = 0
               J = J+J1
            ENDIF
            IF ( J-K3 .EQ. 6 .OR. J .EQ. M+1) THEN
               I = I+1
               J = 1
               GO TO 7660
            ENDIF

            IL  = 3*(I-1)
            I2  = IL
            IR  = I*3-1
            IF (IR .GT. N) IR = N
            J2  = 3*(J-1)-1
            J3  = J2
            J0  = 3*J
            IF (J0 .GT. N) J0 = N
            I3 =  I-2
            J1 =  J-2
            IF (I .EQ. 1) THEN
               IL  = 1
               I2  = IL
               IR  = 2
               I3  = 0
            ENDIF
            IF (J .EQ. 1) THEN
               J2  = 1
               J3  = J2
               J0  = 3
               J1  = 0
            ENDIF
            IF (I3 .LE. J1) THEN
               N1 = J1-I3
               IL = IL+N1
               S  =  ZERO
               DO 7665 K2 = IL,IR
                  S  = S + X(K2)*X(J2) - SIN(RD(K2**2))*SIN(RD(J2**2))
                  J2 = J2+2
                  IF (J2 .GT. J0) GO TO 7670
 7665          CONTINUE
 7670          S1 =  S
               IF ( .NOT. FONLY ) THEN
                  DO 7675 K2=IL,IR
                     G(K2) = G(K2) + TWO*S*X(J3)
                     G(J3) = G(J3) + TWO*S*X(K2)
                     J3    = J3 + 2
                     IF (J3.GT.J0) GO TO 7690
 7675             CONTINUE
               ENDIF
            ELSE
               N1 = I3-J1
               J2 = J2 + 2*N1
               S  =  ZERO
               DO 7680 K2 = J2,J0,2
                  S  =  S + X(K2)*X(IL)- SIN(RD(K2**2))*SIN(RD(IL**2))
                  IL = IL+1
 7680          CONTINUE
               S1 = S
               IF ( .NOT. FONLY ) THEN
                  DO 7685 K2=J2,J0,2
                     G(K2) = G(K2) + TWO*S*X(I2)
                     G(I2) = G(I2) + TWO*S*X(K2)
                     I2    = I2 + 1
 7685             CONTINUE
               ENDIF
            ENDIF
            IF ( .NOT. GONLY ) THEN
               F =  F+ S1**2
            ENDIF
 7690    CONTINUE
      ENDIF
      GO TO 90000

C------------------     FUNCTIONS MNSRF1 AND MNSRF2(N,X,F,G,IFG,WORK,FARG(1))
C  TOINT'S MINIMUM SURFACE FUNCTIONS:
C      MNSRF1: LINEAR AND NON-LINEAR CASES (TOI83B #11 AND #12).
C              SET FARG(1) = 0.
C      MNSRF2: SHIFTED LINEAR AND SHIFTED NON-LINEAR CASES (TOI83B #13 AND #14).
C              SET FARG(1) = 1.

 7700 CONTINUE
      F1 = NINT ( FARG ( 1 ) )
      T  = SQRT(RD(N))
      M  = NINT(T)
      I1 = (M-1)**2
      W1 = RD(I1)

      IF ( .NOT. GONLY ) THEN
        F = ZERO
      ENDIF

      IF ( .NOT. FONLY ) THEN
         DO 7720 I=1,N
            G(I)= ZERO
 7720    CONTINUE
      ENDIF

      DO 7770 I=1,I1
         S  = RD(I)/(T-ONE)
         I2 = INT(S)*M + MOD(I,M-1)
         IF ( MOD(I,M-1) .EQ. 0 ) I2 = I2 - 1
         W3 = (X(I2)-X(I2+M+1))**2+(X(I2+1)-X(I2+M))**2
         W4 = X(I2)**2 / W1
         W5 = X(I2+M+1)**2 / W1
         IF ( .NOT. GONLY ) THEN
            F  = F + SQRT( ONE + W1*W3/TWO )

            IF ( F1 .EQ. 1 ) THEN
               IF ( (I .GE. 1 .AND. I .LE. M-1) .OR.
     -              (MOD(I,M-1) .EQ. 1) ) THEN
                  F = F - W4
               ELSE IF ( (I .GE. I1-M+1 .AND. I .LE. I1) .OR.
     -                   (MOD(I,M-1) .EQ. 0) ) THEN
                  F = F + W5
               ELSE
                  F = F - W4 + W5
               ENDIF
            ENDIF

         ENDIF
         IF ( .NOT. FONLY ) THEN
            K= I2+M+1
            R= TWO*SQRT(ONE+W1*W3/TWO)
            G(I2)   = G(I2)   + W1*(X(I2)-X(K))    /R
            G(K)    = G(K)    - W1*(X(I2)-X(K))    /R
            G(K-1)  = G(K-1)  - W1*(X(I2+1)-X(K-1))/R
            G(I2+1) = G(I2+1) + W1*(X(I2+1)-X(K-1))/R
            IF ( F1 .EQ. 1 ) THEN
               IF ( (I .GE. 1 .AND. I .LE. M-1) .OR.
     -              (MOD(I,M-1) .EQ. 1) ) THEN
                  G(I2)= G(I2)   -  TWO*X(I2)/W1
               ELSE IF ( (I .GE. I1-M+1 .AND. I .LE. I1) .OR.
     -                   (MOD(I,M-1) .EQ. 0) ) THEN
                  G(K) = G(K)    +  TWO*X(K )/W1
               ELSE
                  G(I2)= G(I2)   -  TWO*X(I2)/W1
                  G(K) = G(K)    +  TWO*X(K )/W1
               ENDIF
            ENDIF
         ENDIF
 7770 CONTINUE

      IF ( .NOT. FONLY ) THEN
         DO 7790 I=1,M
            G(I)         = ZERO
            G((I-1)*M+1) = ZERO
            G(I+M*(M-1)) = ZERO
            G(I*M)       = ZERO
 7790    CONTINUE
      ENDIF
      GOTO 90000

C----                      HYPCIR
 7800             CONTINUE
      W1 = X(1)*X(2) - ONE
      W2 = X(1)**2 + X(2)**2 - FOUR
      IF ( .NOT. GONLY ) THEN
         F = W1**2 + W2**2
      ENDIF
      IF ( .NOT. FONLY ) THEN
         G(1) = TWO*(W1*X(2) + TWO*W2*X(1))
         G(2) = TWO*(W1*X(1) + TWO*W2*X(2))
      ENDIF
      GOTO 90000

C----                      SISSER
 7900             CONTINUE
      IF ( .NOT. GONLY ) THEN
         F = THREE*X(1)**4 - TWO*X(1)**2*X(2)**2 + THREE*X(2)**4
      ENDIF
      IF ( .NOT. FONLY ) THEN
         G(1) = R12*X(1)**3 - FOUR*X(1)*X(2)**2
         G(2) = -FOUR*X(1)**2*X(2) + R12*X(2)**3
      ENDIF
      GOTO 90000

C----                      DIXON 7 DIAGONAL DIX7DG
 8000             CONTINUE
      M = N/3
      IF (3*M .NE. N) THEN
         RET = NOFG
         GOTO 90000
      ENDIF

      IF ( .NOT. GONLY ) THEN
         F = ONE
      ENDIF

      W1 =      FARG(1)
      W2 =      FARG(2)
      W3 =      FARG(3)
      W4 =      FARG(4)
      K1 = NINT(FARG(5))
      K2 = NINT(FARG(6))
      K3 = NINT(FARG(7))
      K4 = NINT(FARG(8))

      DO 8090 I = 1,N
         I1 = I  + M
         I2 = I1 + M
         X1 = W1*(RD(I)/RD(N))**K1
         X2 = W2*(RD(I)/RD(N))**K2
         X3 = W3*(RD(I)/RD(N))**K3
         X4 = W4*(RD(I)/RD(N))**K4
         IF ( .NOT. GONLY ) THEN
            F = F + HALF*X1*X(I)**2
            IF ( I .NE. N ) THEN
               F = F + X2*X(I)**2 * (X(I+1)+X(I+1)**2)**2
            ENDIF
            IF ( I1 .LE. N ) THEN
               F = F + X3*X(I)**2 * X(I1)**4
            ENDIF
            IF ( I2 .LE. N ) THEN
               F = F + X4*X(I) * X(I2)
            ENDIF
         ENDIF
         IF ( .NOT. FONLY ) THEN
            G(I) = X1*X(I)
            IF ( I  .NE. N ) THEN
               G(I) = G(I) + TWO*X2*X(I) * (X(I+1)+X(I+1)**2)**2
            ENDIF
            IF ( I  .NE. 1 ) THEN
               X2 = W2*(RD(I-1)/RD(N))**K2
               G(I) = G(I) + TWO*X2*X(I-1)**2 * (X(I)+X(I)**2)*
     -                                             (ONE+TWO*X(I))
            ENDIF
            IF ( I .LE. 2*M ) THEN
               G(I) = G(I) + TWO*X3*X(I) * X(I1)**4
            ENDIF
            IF ( I  .GT. M ) THEN
               X3 = W3*(RD(I-M)/RD(N))**K3
               G(I) = G(I) + FOUR*X3*X(I-M)**2 * X(I)**3
            ENDIF
            IF ( I .LE. M ) THEN
               G(I) = G(I) + X4*X(I2)
            ENDIF
            IF ( I  .GT. 2*M ) THEN
               X4 = W4*(RD(I-2*M)/RD(N))**K4
               G(I) = G(I) + X4*X(I-2*M)
            ENDIF
         ENDIF
 8090 CONTINUE
      GOTO 90000

C----                       MORCIN  MORCIN, MORE COSNARD INTEGRAL PROBLEM
 8100             CONTINUE
      IF ( SIZE .LT. N ) THEN
         RET = NOFG
         GOTO 90000
      ENDIF
      H = ONE/RD(N+1)
      IF ( .NOT. GONLY ) THEN
         F = ZERO
      ENDIF
      W1 = ZERO
      DO 8120 K = 1,N
         TK = RD(K)*H
         W1 = W1 + (ONE-TK)*(X(K)+TK+ONE)**3
 8120 CONTINUE
      W3 = ZERO
      W5 = ZERO
      W6 = ZERO
      DO 8140 K = 1,N
         TK = RD(K)*H
         W2 = X(K) + TK + ONE
         W3 = W3 + TK*W2**3
         W6 = W6 + (ONE-TK)*W2**3
         W4 = X(K) + H*( (ONE-TK)*W3 + TK*(W1-W6) ) /TWO
         W5 = W5 + W4*(ONE-TK)
         WORK(K) = W4
         IF ( .NOT. GONLY ) THEN
            F = F + W4**2
         ENDIF
         IF ( .NOT. FONLY ) THEN
            G(K) = TWO*W4
         ENDIF
 8140 CONTINUE
      IF ( .NOT. FONLY ) THEN
         W1 = ZERO
         W2 = ZERO
         DO 8160 K = 1,N
            TK = RD(K)*H
            G(K) = G(K) +
     -            H*THREE*(X(K)+TK+ONE)**2*(TK*(W5-W2) + (ONE-TK)*W1)
            W1 = W1 + WORK(K)*TK
            W2 = W2 + WORK(K)*(ONE-TK)
 8160    CONTINUE
      ENDIF
      GOTO 90000
C----                      BOOTH'S QUADRATIC
 8200             CONTINUE
      W1 =     X(1) + TWO*X(2) - SEVEN
      W2 = TWO*X(1) +     X(2) -  FIVE

      IF ( .NOT. GONLY ) THEN
         F = W1**2 + W2**2
      ENDIF
      IF ( .NOT. GONLY ) THEN
         G(1) =  TWO*W1 + FOUR*W2
         G(2) = FOUR*W1 +  TWO*W2
      ENDIF
      GOTO 90000
C----                      FROM POWELL
 8300             CONTINUE
      W1 = X(1) + TENTH
      W2 = (TEN * X(1) / W1) + TWO * X(2)**2
      IF ( .NOT. GONLY ) THEN
         F = X(1)**2 + W2**2
      ENDIF
      IF ( .NOT. FONLY ) THEN
         G(1) = TWO * X(1) + TWO*W2 / W1**2
         G(2) = EIGHT * X(2) * W2
      ENDIF

      GOTO 90000


C## E X I T
90000       STATUS(7) =  RET
91000       GFIRST    = .FALSE.
      RETURN

C## F O R M A T S:  NONE ARE DEFINED.
C##                 E N D         OF ZZFNS4.
                    END
