      SUBROUTINE ZZFNS1 ( N, X, F, G, WORK, SIZE, FONLY, GONLY,
     -                    FIRST, FARG, FUNCNO, STATUS )

C## A R G U M E N T S:
                 INTEGER            N, SIZE, FUNCNO, STATUS(7)

                 LOGICAL            FIRST,   FONLY,  GONLY

                 DOUBLE PRECISION   F, X(N), G(N), WORK(SIZE), FARG(*)
C!!!!            REAL               F, X(N), G(N), WORK(SIZE), FARG(*)

C## S T A T U S:
C               SINGLE/DOUBLE CONVERSION: NEEDED (SEE CONVRT).
C
C               IGNORE LINES BEGINNING WITH  "C!!!!" .
C
C               THIS VERSION IS IN   D O U B L E   PRECISION.
C!!!!           THIS VERSION IS IN   S I N G L E   PRECISION.
C
C               SYSTEM  DEPENDENCE:                      NONE.
C
C>RCS $HEADER: FNS1.F,V 1.2 92/01/07 15:12:07 BUCKLEY EXP $
C>RCS $LOG:     FNS1.F,V $
C>RCSREVISION 1.2  92/01/07  15:12:07  BUCKLEY
C>RCSMINOR FIX
C>RCS
C>RCSREVISION 1.1  91/11/20  10:52:47  BUCKLEY
C>RCSFINAL SUBMISSION TO TOMS
C>RCS
C>RCSREVISION 1.0  90/07/31  13:01:55  BUCKLEY
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
C## E N T R Y   P O I N T S:  ZZFNS1    THE NATURAL ENTRY POINT.
C
C## S U B R O U T I N E S:
C
C     PREDEFINED FUNCTIONS : SIN, COS, TAN, ACOS, ATAN, ABS, MAX, NINT
C                            EXP, LOG, MIN, MOD, SIGN, SQRT, REAL(DBLE)
C
C     STATEMENT FUNCTION:      RD
C
C## P A R A M E T E R S:

      DOUBLE PRECISION  ZERO,       ONE,       TWO,       THREE
C!!!! REAL              ZERO,       ONE,       TWO,       THREE
      PARAMETER (       ZERO = 0D0, ONE = 1D0, TWO = 2D0, THREE = 3D0)

      DOUBLE PRECISION  FOUR,       FIVE,      SIX,       SEVEN
C!!!! REAL              FOUR,       FIVE,      SIX,       SEVEN
      PARAMETER (       FOUR = 4D0, FIVE = 5D0,SIX = 6D0, SEVEN = 7D0)

      DOUBLE PRECISION  EIGHT,         NINE,          TEN
C!!!! REAL              EIGHT,         NINE,          TEN
      PARAMETER (       EIGHT = 8D0,   NINE = 9D0,    TEN = 10D0     )

      DOUBLE PRECISION  TENTH,         FIFTH,        HALF
C!!!! REAL              TENTH,         FIFTH,        HALF
      PARAMETER (       TENTH = .1D0,  FIFTH = .2D0, HALF = .5D0      )

      DOUBLE PRECISION  RPT9,        RPT8,        RD29
C!!!! REAL              RPT9,        RPT8,        RD29
      PARAMETER (       RPT9 = .9D0, RPT8 = .8D0, RD29 = 1D0/29D0 )
      DOUBLE PRECISION  R11,        R12,        R13,       R14
C!!!! REAL              R11,        R12,        R13,       R14
      PARAMETER (       R11 = 11D0, R12 = 12D0, R13 = 13D0,R14 = 14D0)

      DOUBLE PRECISION  R15,        R16,        R17,       R18
C!!!! REAL              R15,        R16,        R17,       R18
      PARAMETER (       R15 = 15D0, R16 = 16D0, R17 = 17D0,R18 = 18D0)

      DOUBLE PRECISION  R19,        R20,        R25,       R29
C!!!! REAL              R19,        R20,        R25,       R29
      PARAMETER (       R19 = 19D0, R20 = 20D0, R25 = 25D0,R29 = 29D0)

      DOUBLE PRECISION  R32,        R36,        R40,       R42
C!!!! REAL              R32,        R36,        R40,       R42
      PARAMETER (       R32 = 32D0, R36 = 36D0, R40 = 40D0,R42 = 42D0)

      DOUBLE PRECISION  R45,        R49
C!!!! REAL              R45,        R49
      PARAMETER (       R45 = 45D0, R49 = 49D0 )

      DOUBLE PRECISION  R50,        R56,        R84,       R90
C!!!! REAL              R50,        R56,        R84,       R90
      PARAMETER (       R50 = 50D0, R56 = 56D0, R84 = 84D0,R90 = 90D0)

      DOUBLE PRECISION  R100,            R180,           R200
C!!!! REAL              R100,            R180,           R200
      PARAMETER (       R100 = 100D0,    R180 = 180D0,   R200 = 200D0 )

      DOUBLE PRECISION  R256,            R360,           R400
C!!!! REAL              R256,            R360,           R400
      PARAMETER (       R256 = 256D0,    R360 = 360D0,   R400 = 400D0 )

      DOUBLE PRECISION  R600,            R681,           R991
C!!!! REAL              R600,            R681,           R991
      PARAMETER (       R600 = 600D0,    R681 = 681D0,   R991 = 991D0 )

      DOUBLE PRECISION  R1162,                 R2324
C!!!! REAL              R1162,                 R2324
      PARAMETER (       R1162 = 1162D0,        R2324 = 2324D0         )

      DOUBLE PRECISION  R10000,                R40000
C!!!! REAL              R10000,                R40000
      PARAMETER (       R10000 = 10000D0,      R40000 = 40000D0       )
      DOUBLE PRECISION  R1PD6,                 R2PDM6
C!!!! REAL              R1PD6,                 R2PDM6
      PARAMETER (       R1PD6 = 1D6,           R2PDM6 = 2D-6         )

      DOUBLE PRECISION  RP04,        RP01,          R1PZ1
C!!!! REAL              RP04,        RP01,          R1PZ1
      PARAMETER (       RP04 = 4D-2, RP01 = .01D0,  R1PZ1 = 1.0001D0 )

      DOUBLE PRECISION  R1P2,         R7P5,         RP1136
C!!!! REAL              R1P2,         R7P5,         RP1136
      PARAMETER (       R1P2 = 1.2D0, R7P5 = 7.5D0, RP1136 = 0.1136D0 )

      DOUBLE PRECISION  R1P5,         R2P5,         R2P625
C!!!! REAL              R1P5,         R2P5,         R2P625
      PARAMETER (       R1P5 = 1.5D0, R2P5 = 2.5D0, R2P625 = 2.625D0 )

      DOUBLE PRECISION  R10P1,         R19P8,         R20P2
C!!!! REAL              R10P1,         R19P8,         R20P2
      PARAMETER (       R10P1 = 10.1D0,R19P8 = 19.8D0,R20P2 = 20.2D0 )

      DOUBLE PRECISION  R2D3,          R4D3,          R7D3
C!!!! REAL              R2D3,          R4D3,          R7D3
      PARAMETER (       R2D3 = 2D0/3D0,R4D3 = 4D0/3D0,R7D3 = 7D0/3D0 )

      DOUBLE PRECISION  R2P25
C!!!! REAL              R2P25
      PARAMETER (       R2P25 = 2.25D0 )

C## L O C A L   D E C L:

      INTEGER     OK,  ABORT,  LIMIT,  NOF,  NOG,  NOFG

      INTEGER     I1, I2, IB, IY, I, IDUMMY, J, K, RET

      LOGICAL           GFIRST

      DOUBLE PRECISION  ZZMPAR, HUGE
C!!!! REAL              ZZMPAR, HUGE

C--------- VARIABLES FOR THE TEST FUNCTIONS.

      DOUBLE PRECISION  X1,     X2,     X3,    X4,    X5,    X6
C!!!! REAL              X1,     X2,     X3,    X4,    X5,    X6

      DOUBLE PRECISION  G1,     G2,     G3,    G4,    G5,    G6
C!!!! REAL              G1,     G2,     G3,    G4,    G5,    G6

      DOUBLE PRECISION  W1,     W2,     W3,    W4,    W5,    W6
C!!!! REAL              W1,     W2,     W3,    W4,    W5,    W6

      DOUBLE PRECISION  W7,     W8,     W9,    W10,   W11,   W12
C!!!! REAL              W7,     W8,     W9,    W10,   W11,   W12

      DOUBLE PRECISION   T,     BIGGST,        SMLLST
C!!!! REAL               T,     BIGGST,        SMLLST

      DOUBLE PRECISION  RI,     TI,     YI,    PI,    AI,    BI
C!!!! REAL              RI,     TI,     YI,    PI,    AI,    BI

      DOUBLE PRECISION R2P,     RD,     TPI,   TPIS
C!!!! REAL             R2P,     RD,     TPI,   TPIS

      DOUBLE PRECISION RF1,     RF2
C!!!! REAL             RF1,     RF2

C--------- DATA ARRAYS FOR FUNCTIONS

      DOUBLE PRECISION  BARD7Y (15)
C!!!! REAL              BARD7Y (15)

C## S A V E:

      SAVE  GFIRST, PI, TPI, TPIS, R2P, BIGGST, SMLLST, HUGE
      SAVE  OK, ABORT, LIMIT, NOF, NOG, NOFG

C--------- SAVE DATA ARRAYS FOR THE TEST FUNCTIONS.

      SAVE  BARD7Y

C## E Q U I V A L E N C E S: NONE ARE DEFINED.
C## C O M M O N:             NONE IS DEFINED.
C## D A T A:

      DATA  GFIRST/.TRUE./

C--------- DATA FOR MINTEST FUNCTION     BARD70

      DATA  BARD7Y( 1), BARD7Y( 2), BARD7Y( 3), BARD7Y( 4), BARD7Y( 5)
     -     /  .14 D0  ,   .18 D0  ,   .22 D0  ,   .25 D0  ,   .29 D0  /

      DATA  BARD7Y( 6), BARD7Y( 7), BARD7Y( 8), BARD7Y( 9), BARD7Y(10)
     -     /  .32 D0  ,   .35 D0  ,   .39 D0  ,   .37 D0  ,   .58 D0  /

      DATA  BARD7Y(11), BARD7Y(12), BARD7Y(13), BARD7Y(14), BARD7Y(15)
     -     /  .73 D0  ,   .96 D0  ,  1.34 D0  ,  2.10 D0  ,  4.39 D0  /

C##                                                E X E C U T I O N
C##                                                E X E C U T I O N
C---------                     FUNCTION DEFINITION
      RD (IDUMMY)  = DBLE (IDUMMY)
C!!!! RD (IDUMMY)  = REAL (IDUMMY)

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
     -     1100, 1000, 1300, 1000, 1000, 1000, 1000, 1000, 1000, 1000,
     -     2100, 1000, 1000, 1000, 1000, 2600, 1000, 1000, 1000, 1000,
     -     1000, 1000, 1000, 3400, 3500, 3600, 3700, 1000, 1000, 1000,
     -     1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 5000,
     -     5100, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000,
     -     1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000,
     -     1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000,
     -     1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000
     -    ) FUNCNO
 1000 STATUS(7) = ABORT
      GOTO 91000


C---------                             ROSENB
C       AI IS ALPHA, TAKEN FROM FARG 1
C       BI IS BETA,  TAKEN FROM FARG 2
C       IY IS GAMMA, TAKEN FROM FARG 3

 1100 CONTINUE
      AI = FARG(1)
      BI = FARG(2)
      IY = NINT(FARG(3))
      W1 = X(1)

      IF ( .NOT. GONLY ) THEN
         F = BI * (ONE-W1)**2
      ENDIF

      IF ( .NOT. FONLY ) THEN
         G(1) = -TWO * BI * (ONE-W1)
      ENDIF

      DO 1110 I=2,N
         W2 = X(I)
         W4 = W1**(IY-1)
         W3 = W2 - W1*W4
         IF ( .NOT. GONLY ) THEN
            F = F + AI * W3**2
         ENDIF
         IF ( .NOT. FONLY ) THEN
            G(I-1) = G(I-1) - TWO * AI * IY * W4 * W3
            G(I) = TWO*AI * W3
         ENDIF
         W1 = W2

 1110 CONTINUE
      GOTO 90000

C--------- MINTEST FUNCTION     PWSING


 1300 IF ( .NOT. GONLY )  THEN
         F = ZERO
      ENDIF

      IF ( 4 * (N/4) .NE. N ) THEN

         IF ( .NOT. FONLY ) THEN
            DO 1310 I = 1,N
               G(I) = ZERO
 1310       CONTINUE
         ENDIF

      ELSE

         DO 1320 I=1,N/4

            J  = 4*I

            W1 = X(J-3)
            W2 = X(J-2)
            W3 = X(J-1)
            W4 = X(J  )

            W5 = W1 + TEN * W2
            W6 = W3 - W4
            W2 = W2 - TWO * W3
            W3 = W2 * W2 *W2
            W1 = W1 - W4
            W4 = W1 * W1 * W1

            IF ( .NOT. GONLY ) THEN
               F = F + W5*W5 + FIVE*W6*W6 + W2*W3 + TEN*W1*W4
            ENDIF

            IF ( .NOT. FONLY ) THEN
               G(J-3) =   TWO  * W5  + R40   * W4
               G(J-2) =   R20  * W5  + FOUR  * W3
               G(J-1) =   TEN  * W6  - EIGHT * W3
               G(J  ) =  -TEN  * W6  - R40   * W4
            ENDIF

 1320    CONTINUE
      ENDIF
      GOTO 90000

C--------- MINTEST FUNCTION   BIGGS ( IFG, N, X, F, G, NINT(FARG(1)))
C--------- NINT(FARG(1)) IS M


 2100 IB = NINT(FARG(1))

      X1 = X(1)
      X2 = X(2)
      IF ( N .EQ. 2 ) THEN
         X3 = ONE
         X4 = FIVE
         X5 = ZERO
         X6 = ZERO
      ELSE IF ( N .EQ. 3 ) THEN
         X3 = ONE
         X4 = X(3)
         X5 = ZERO
         X6 = ZERO
      ELSE IF ( N .EQ. 4 ) THEN
         X3 = X(3)
         X4 = X(4)
         X5 = ZERO
         X6 = ZERO
      ELSE IF ( N .EQ. 5 ) THEN
         X3 = X(3)
         X4 = X(4)
         X5 = X(5)
         X6 = THREE
      ELSE IF ( N .EQ. 6 ) THEN
         X3 = X(3)
         X4 = X(4)
         X5 = X(5)
         X6 = X(6)
      ENDIF

      IF ( .NOT. GONLY ) THEN
         F = ZERO
      ENDIF

      IF ( .NOT. FONLY ) THEN
         G(1) = ZERO
         G(2) = ZERO
         G(3) = ZERO
         G(4) = ZERO
         G(5) = ZERO
         G(6) = ZERO
      ENDIF

      DO 2110 I = 1, IB
         T  = RD(I)
         TI = T/TEN
         IF ( MAX(-T,-TI*FOUR,-TI*X1,-TI*X2,-TI*X5) .LE. BIGGST ) THEN
            IF ( N .LE. 4 ) THEN
               YI = EXP(-TI) - FIVE * EXP(-T)
            ELSE
               YI = EXP(-TI) - FIVE * EXP(-T) + THREE*EXP(-FOUR*TI)
            END IF
            W1 = EXP(-TI*X1)
            W2 = EXP(-TI*X2)
            W5 = EXP(-TI*X5)
         ELSE
            RET = NOFG
            GOTO 90000
         ENDIF
         RI = X3*W1 - X4*W2 + X6*W5 - YI

         IF ( .NOT. GONLY ) THEN
            F = F + RI*RI
         ENDIF

         IF ( .NOT. FONLY ) THEN
            W7   = TI*RI
            G(1) = G(1) - TWO * X3 * W1 * W7
            G(2) = G(2) + TWO * X4 * W2 * W7

            IF ( N .EQ. 3 ) THEN
               G(3) = G(3) - TWO * W2 * RI
            ELSE IF ( N .GE. 4 ) THEN
               G(3) = G(3) + TWO * W1 * RI
               G(4) = G(4) - TWO * W2 * RI

               IF ( N .GE. 5 ) THEN
                  G(5) = G(5) - TWO * X6 * W5 * W7
                  IF ( N .EQ. 6 ) G(6) = G(6) + TWO * W5 * RI
               ENDIF

            ENDIF
         ENDIF

 2110 CONTINUE

      GOTO 90000


C--------- MINTEST FUNCTION   BOX663 ( N, X, F, G, IFG, NINT(FARG(1)))
C--------- NINT(FARG(1)) IS M

 2600 IF ( .NOT. GONLY ) THEN
         F = ZERO
      ENDIF

      IF ( .NOT. FONLY ) THEN
         G1 = ZERO
         G2 = ZERO
         IF ( N .EQ. 3 ) G3 = ZERO
      ENDIF

      DO 2610  I = 1,NINT(FARG(1))
         W2 = RD(I)
         TI = W2/TEN
         IF ( MAX(-W2,-TI,-TI*X(1),-TI*X(2)) .LE. BIGGST ) THEN
            W3 = EXP(-TI * X(1))
            W4 = EXP(-TI * X(2))
            W5 = EXP(-TI) - EXP(-W2)
         ELSE
            RET = NOFG
            GOTO 90000
         ENDIF
         IF ( N .EQ. 3 ) THEN
            RI = W3 - W4 - W5*X(3)
         ELSE
            RI = W3 - W4 - W5
         ENDIF

         IF ( .NOT. GONLY ) THEN
            IF ( ABS(RI) .LE. SQRT(HUGE-MAX(F,ZERO)) ) THEN
               F = F + RI**2
            ELSE
               RET = NOFG
               GOTO 90000
            ENDIF
         ENDIF

         IF ( .NOT. FONLY ) THEN
            W2 = TI*RI
            G1 = G1 - W3*W2
            G2 = G2 + W4*W2
            IF ( N .EQ. 3 ) G3 = G3 - W5*RI
         ENDIF

 2610 CONTINUE

      IF ( .NOT. FONLY ) THEN
         G(1) = TWO * G1
         G(2) = TWO * G2
         IF ( N .EQ. 3 ) G(3) = TWO * G3
      ENDIF
      GOTO 90000


C--------- MINTEST FUNCTION     SCHMVT     EXTENSION DUE TO TOI83B #35.
C---- FARG(1) IS ALPHA.


 3400 CONTINUE
      AI = FARG(1)
      IF ( .NOT. GONLY ) THEN
         F = ZERO
      ENDIF

      IF ( .NOT. FONLY ) THEN
         DO 3410 I = 1, N
            G(I) = ZERO
 3410    CONTINUE
      END IF

      DO 3420 I = 1, N-2
         W1 = X(I) - X(I+1)
         W2 = X(I) + X(I+2)

         W3 = ONE + W1*W1
         W4 = (PI*X(I+1) + X(I+2)) / TWO
         W5 = (W2/X(I+1)) - TWO
         IF ( -W5**2 .LE. BIGGST ) THEN
            W6 =  EXP(-W5*W5)
         ELSE
            RET = NOFG
            GOTO 90000
         ENDIF

         IF ( .NOT. GONLY ) THEN
            F = F + AI - ((ONE/W3) + SIN(W4) + W6 )
         ENDIF

         IF ( .NOT. FONLY ) THEN

            W3 = TWO*W1/(W3*W3)
            W4 = COS(W4)/TWO
            W6 = TWO*W5*W6/X(I+1)

            G(I  ) = G(I  ) + W3 +    W6
            G(I+1) = G(I+1) - W3 - PI*W4 - W6*W2/X(I+1)
            G(I+2) = G(I+2) - W4 +    W6

         ENDIF

 3420 CONTINUE

      GOTO 90000

C--------- MINTEST FUNCTION     ENGVL2


 3500 X1 = X(1)
      X2 = X(2)
      X3 = X(3)

      W1 = X1*X1
      W2 = X1*W1
      W3 = X2*X2
      W4 = X3*X3

      W5 = X3 - TWO
      W6 = FIVE*X3 - X1 + ONE
      W7 = W1 + W3 - ONE

      W8  = W7 + W4
      W9  = W7 + W5*W5
      W10 = X1 + X2 + X3 - ONE
      W11 = X1 + X2 - X3 + ONE
      W12 = W2 + THREE*W3 + W6*W6 - R36

      IF ( .NOT. GONLY ) THEN
         F = W8*W8 + W9*W9 + W10*W10 + W11*W11 + W12*W12
      ENDIF

      IF ( .NOT. FONLY ) THEN
         W10  = W8 + W9
         G(1) = TWO*(TWO*X1*W10 + TWO*(X1+X2) + W12*(THREE*W1-TWO*W6))
         G(2) = TWO*(TWO*X2*W10 + TWO*(X1+X2) + SIX*W12*X2)
         G(3) = TWO*(TWO*(W8*X3+W5*W9) + TWO*X3 - TWO + TEN*W12*W6)
      ENDIF

      GOTO 90000


C--------- MINTEST FUNCTION     BARD70


 3600 X1 = X(1)
      X2 = X(2)
      X3 = X(3)

      IF ( .NOT. GONLY ) THEN
         F = ZERO
      ENDIF

      IF ( .NOT. FONLY ) THEN
         G1 = ZERO
         G2 = ZERO
         G3 = ZERO
      ENDIF

      DO 3610 I=1,15

         W1 = RD(I)
         W2 = RD(16-I)
         W3 = MIN(W1,W2)

         W4 = X2*W2 + X3*W3
         RI = BARD7Y(I) - (X1 + W1/W4)
         W4 = W4*W4

         IF ( .NOT. GONLY ) THEN
            F = F + RI*RI
         ENDIF

         IF ( .NOT. FONLY ) THEN
            W4 = RI*W1/W4
            G1 = G1 - RI
            G2 = G2 + W2*W4
            G3 = G3 + W3*W4
         ENDIF

 3610 CONTINUE

      IF ( .NOT. FONLY ) THEN
         G(1) = G1*TWO
         G(2) = G2*TWO
         G(3) = G3*TWO
      ENDIF

      GOTO 90000

C--------- MINTEST FUNCTION     CRGLVY


 3700 K = (N-2) / 2
      F = ZERO

      IF ( .NOT. GONLY ) THEN
         DO 3720 I = 1, K
            W1 = X(2*I)   - X(2*I+1)
            W2 = X(2*I+1) - X(2*I+2)
            W3 = X(2*I+2) - ONE

            IF ( X(2*I-1) .LE. BIGGST ) THEN
               W4 =  EXP(X(2*I-1))
            ELSE
               RET = NOFG
               GOTO 90000
            ENDIF
            W5 =  W4 - X(2*I)
            W6 =  TAN(W2)

            F = F  +  W5**4 + R100*W1**6 + W6**4 + X(2*I-1)**8 + W3*W3
 3720    CONTINUE
      ENDIF

      IF ( .NOT. FONLY ) THEN

         W1 = EXP(X(1))
         W2 = FOUR*(W1 - X(2))**3
         G(1) =  W1*W2 + EIGHT*X(1)**7
         G(2) = -W2 + R600*(X(2) - X(3))**5

         DO 3740 I = 3, N-2
            IF ( INT(I/2)*2 .NE. I ) THEN
C   ODD I
               W1 = R600 * (X(I-1) - X(I))**5
               W2 = X(I) - X(I+1)
               W3 = FOUR * TAN(W2)**3 / COS(W2)**2
               W4 = EXP(X(I))
               W5 = FOUR*(W4 - X(I+1))**3
               G(I) = -W1 + W3 + W4*W5 + EIGHT*X(I)**7
            ELSE
C   EVEN I
               W6 = R600 * (X(I) - X(I+1))**5
               G(I) = -W3 + TWO*(X(I) - ONE) - W5 + W6
            ENDIF
 3740    CONTINUE

         W2 = X(N-1) - X(N)
         W3 = FOUR * TAN(W2)**3 / COS(W2)**2
         G(N-1) = -R600 * (X(N-2) - X(N-1))**5  +  W3
         G(N  ) = -W3  +  TWO*(X(N) - ONE)
      ENDIF
      GOTO 90000

C--------- MINTEST FUNCTION     PENAL1 ( N, X, F, G, IFG,
C                                         FARG(1), FARG(2)      )
C--------- FARG(1) IS A
C--------- FARG(2) IS B

 5000 RF1 = FARG ( 1 )
      RF2 = FARG ( 2 )

      W1 = - ONE / FOUR
      W2 =  ZERO

      DO 5010 J = 1, N
         W3 = X(J)
         W1 = W1 + W3*W3
         W3 = W3 - ONE
         W2 = W2 + W3*W3
 5010 CONTINUE

      IF ( .NOT. GONLY ) THEN
         F = RF1*W2 + RF2 *W1*W1
      ENDIF

      IF ( .NOT. FONLY ) THEN
         W1 = FOUR*RF2*W1
         W2 =  TWO*RF1
         DO 5020 J = 1, N
            W3   = X(J)
            G(J) = W2 * (W3 - ONE) + W3*W1
 5020    CONTINUE
      ENDIF

      GOTO 90000


C--------- MINTEST FUNCTION     PENAL2 ( N, X, F, G, IFG,
C                                         FARG(1), FARG(2), WORK, SIZE)
C--------- FARG(1) IS A
C--------- FARG(2) IS B

 5100 RF1 = FARG ( 1 )
      RF2 = FARG ( 2 )

      IF ( SIZE .LT. 2 * N ) THEN
         F = ZERO
         DO 5110 K = 1, N
            G(K) = ZERO
 5110    CONTINUE
         GO TO 90000
      ENDIF

      W1 = EXP(TENTH)
      W2 = EXP(-TENTH)
      W3 = ZERO

      I1 = 0
      I2 = N

      DO 5120 K = 1, N
         W4 = X(K)
         W3  = W3 + RD( N - K + 1 ) * W4 * W4
         IF ( TENTH*W4 .LE. BIGGST ) THEN
            W5 = EXP (TENTH * W4)
         ELSE
            RET = NOFG
            GOTO 90000
         ENDIF

         IF ( K .EQ. 1 ) THEN
            W6 = ZERO
            W7 = ONE

         ELSE
            W7  = W9 * W1
            W10 = W5 + W8 - (W7 + W9)
            W11 = W5 - W2

            IF ( .NOT. FONLY ) THEN
               WORK(I1+K) = W10
               WORK(I2+K) = W11
            ENDIF

            IF ( .NOT. GONLY ) THEN
               W6 = W6 + W10*W10 + W11*W11
            ENDIF

         ENDIF

         W8 = W5
         W9 = W7

 5120 CONTINUE

      W1 = X(1) - FIFTH
      W2 = W3   - ONE

      IF ( .NOT. GONLY ) THEN
         F = RF1 * W6  +  RF2* ( W1*W1 + W2*W2 )
      ENDIF

      IF ( .NOT. FONLY ) THEN
         W3 = FIFTH * RF1
         W2 = FOUR  * RF2 * W2

         DO 5130 K = 1, N

C        ---NOTE THAT W8 DOES NOT NEED TO BE PRE-DEFINED WHEN K = 1.

            W4 = X(K)
            IF ( TENTH*W4 .LE. BIGGST ) THEN
               W5 = EXP(TENTH * W4)
            ELSE
               RET = NOFG
               GOTO 90000
            ENDIF
            W6 = W8
            W7 = WORK(I2+K)

            IF ( K .LT. N ) THEN
               W8 = WORK(I1+K+1)
               IF ( K .EQ. 1 ) THEN
                  G(1) = W3 * W5 * (           W8 )
     -                 + W2 * W4 * RD(N)  +  W1 * TWO * RF2

               ELSE
                  G(K) = W3 * W5 * ( W6 + W7 + W8 )
     -                 + W2 * W4 * RD( N - K + 1 )

               ENDIF

            ELSE
                  G(N) = W3 * W5 * ( W6 + W7      )
     -                 + W2 * W4

            ENDIF

 5130    CONTINUE

      ENDIF

      GOTO 90000

C## E X I T
90000       STATUS(7) =  RET
91000       GFIRST    = .FALSE.
      RETURN

C## F O R M A T S:  NONE ARE DEFINED.
C##                 E N D         OF ZZFNS1.
                    END
