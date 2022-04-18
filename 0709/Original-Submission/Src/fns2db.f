      SUBROUTINE ZZFNS2 ( N, X, F, G, WORK, SIZE, FONLY, GONLY,
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
C>RCS $HEADER: FNS2.F,V 1.1 91/11/20 10:52:48 BUCKLEY EXP $
C>RCS $LOG:     FNS2.F,V $
C>RCSREVISION 1.1  91/11/20  10:52:48  BUCKLEY
C>RCSFINAL SUBMISSION TO TOMS
C>RCS
C>RCSREVISION 1.0  90/07/31  13:01:56  BUCKLEY
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
C## E N T R Y   P O I N T S:  ZZFNS2    THE NATURAL ENTRY POINT.
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

      INTEGER           ALPHA,          BETA,          GAMMA
      PARAMETER (       ALPHA = 5,      BETA = 14,     GAMMA = 3     )

C## L O C A L   D E C L:

      INTEGER     OK,  ABORT,  LIMIT,  NOF,  NOG,  NOFG

      INTEGER     I1, IC, IR, IS, IV, IDUMMY, I, J, K, M, RET, ALPHA2

      LOGICAL           GFIRST

      DOUBLE PRECISION  ZZMPAR, HUGE
C!!!! REAL              ZZMPAR, HUGE

C--------- VARIABLES FOR THE TEST FUNCTIONS.

      DOUBLE PRECISION  X1,     X2,     X3,    X4,    G1,    G2
C!!!! REAL              X1,     X2,     X3,    X4,    G1,    G2

      DOUBLE PRECISION  W1,     W2,     W3,    W4,    W5,    W6
C!!!! REAL              W1,     W2,     W3,    W4,    W5,    W6

      DOUBLE PRECISION  W7,     W8,     W9,    W10
C!!!! REAL              W7,     W8,     W9,    W10

      DOUBLE PRECISION  BIGGST, SMLLST
C!!!! REAL              BIGGST, SMLLST

      DOUBLE PRECISION  RI,     XI,     PI
C!!!! REAL              RI,     XI,     PI

      DOUBLE PRECISION R2P,     RD,     TPI,   TPIS
C!!!! REAL             R2P,     RD,     TPI,   TPIS

      DOUBLE PRECISION RF1, RF2, RF3, RF4
C!!!! REAL             RF1, RF2, RF3, RF4

C## S A V E:

      SAVE  GFIRST, PI,  TPI, TPIS, R2P, BIGGST, SMLLST
      SAVE  HUGE, M
      SAVE   OK, ABORT, LIMIT, NOF, NOG, NOFG

C## E Q U I V A L E N C E S: NONE ARE DEFINED.
C## C O M M O N:             NONE IS DEFINED.
C## D A T A:

      DATA  GFIRST/.TRUE./

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
     -     1000, 1200, 1000, 1400, 1500, 1600, 1700, 1800, 1000, 2000,
     -     1000, 2200, 2300, 2400, 2500, 1000, 2700, 2800, 2900, 3000,
     -     3100, 3200, 3300, 1000, 1000, 1000, 1000, 3800, 3900, 1000,
     -     1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000,
     -     1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000,
     -     1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000,
     -     1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000,
     -     1000, 1000, 1000, 1000, 8500, 8600, 8700, 8800, 8900, 1000
     -    ) FUNCNO
 1000 STATUS(7) = ABORT
      GOTO 91000


C---------                                   HELIX
 1200 X1 = X(1)
      X2 = X(2)
      X3 = X(3)

      W4 = X1*X1 + X2*X2
      W2 = SQRT(W4)
      W1 = R2P*ATAN( X2/X1 )

      IF ( X1 .LT. ZERO ) THEN
         W1 = W1 + HALF
      ENDIF

      W3 = X3 - TEN*W1
      W1 = W2 - ONE

      IF ( .NOT. GONLY ) THEN
         F = R100 * ( W3*W3 + W1*W1 ) + X3*X3
      ENDIF

      IF ( .NOT. FONLY ) THEN
         W2   = ONE - ONE/W2
         W4   = TEN * W3 * R2P / W4
         G(1) = R200 * ( X2 * W4 + X1 * W2 )
         G(2) = R200 * (-X1 * W4 + X2 * W2 )
         G(3) = R200*W3 + TWO*X3
      ENDIF
      GOTO 90000

C----                         WOODS
 1400               CONTINUE
      IF ( .NOT. GONLY ) THEN
         F = ZERO
      ENDIF
      M = N/4
      IF (M*4 .NE. N) THEN
         RET = NOFG
         GOTO 90000
      ENDIF
      DO 1480 K = 1,M
          X1 = X(4*K-3)
          X2 = X(4*K-2)
          X3 = X(4*K-1)
          X4 = X(4*K  )

          W1 = X2 - X1*X1
          W2 = ONE - X1
          W3 = X4 - X3*X3
          W4 = ONE - X3
          W5 = X2 - ONE
          W6 = X4 - ONE

          IF (.NOT. GONLY) THEN
             F = F + R100*W1*W1 + W2*W2 + R90*W3*W3 + W4*W4
     -             + R10P1*(W5*W5 + W6*W6) + R19P8*W5*W6
          ENDIF

          IF (.NOT. FONLY) THEN
             G(4*K-3) = -R400*X1*W1 -    TWO*W2
             G(4*K-2) =  R200 *  W1 +  R20P2*W5  +  R19P8*W6
             G(4*K-1) = -R360*X3*W3 -    TWO*W4
             G(4*K  ) =  R180 *  W3 +  R20P2*W6  +  R19P8*W5
          ENDIF

 1480 CONTINUE
      GOTO 90000


C--------- MINTEST FUNCTION     NONDIA


 1500 IF ( .NOT. GONLY ) THEN
         F = ZERO
      ENDIF

      IF ( .NOT. FONLY ) THEN
         G1 = ZERO
      ENDIF

      X1 = X(1)

      DO 1510 I=2,N

         XI = X(I)
         W2 = X1 - XI*XI
         W3 = ONE - XI

         IF ( .NOT. GONLY ) THEN
            F = F + R100* W2*W2 + W3 * W3
         ENDIF

         IF ( .NOT. FONLY ) THEN
            G1   = G1 + W2
            G(I) = - ( R400 * XI * W2 )  -  TWO * W3
         ENDIF

 1510 CONTINUE

      IF ( .NOT. FONLY ) THEN
         G(1) = R200 * G1
      ENDIF

      GOTO 90000


C--------- MINTEST FUNCTION     TRIDIA ( N, X, F, G, IFG,
C                                   FARG(1), FARG(2), FARG(3), FARG(4))
C           RF1 IS ALPHA.
C           RF2 IS BETA.
C           RF3 IS GAMMA.
C           RF4 IS DELTA.


 1600 W1 = X(1)
      W4 = ZERO
      RF1 = FARG(1)
      RF2 = FARG(2)
      RF3 = FARG(3)
      RF4 = FARG(4)

      IF ( .NOT. GONLY ) THEN
            F  = RF3 * (RF4*W1 - ONE)**2
      ENDIF

      DO 1610 I=2,N

         W2 = X(I)
         W3 = RF1*W2 - RF2*W1

         IF ( .NOT. GONLY ) THEN
            F = F + RD(I) * W3**2
         ENDIF

         IF ( .NOT. FONLY ) THEN
            W5     = TWO * RD(I) * W3
            G(I-1) = RF1*W4  - RF2*W5
            W4     = W5
         ENDIF

         W1 = W2

 1610 CONTINUE

      IF ( .NOT. FONLY ) THEN
         G(1) = G(1) + TWO * RF3 * RF4 * ( RF4*X(1) - ONE )
         G(N) = RF1 * W4
      ENDIF

      GOTO 90000


C--------- MINTEST FUNCTION     POWER
C      W3 = FARG(1) IS C.
C      W5 = FARG(2) IS ALPHA.
C      W6 = FARG(3) IS BETA.

 1700 W1 = ZERO
      W3 = FARG(1)
      W5 = FARG(2)
      W6 = FARG(3)

      DO 1710 I=1,N
         W2 = X(I)
         W1 = W1 + W2 * W2 * RD(I)**W6
 1710  CONTINUE

      W4 = W1 + W3
      IF ( .NOT. GONLY ) THEN
         F = W4**W5
      ENDIF

      IF ( .NOT. FONLY ) THEN

         DO 1720 I = 1,N
            G(I) = TWO * W5 * RD(I)**W6 * X(I) * W4**(W5-ONE)
 1720    CONTINUE

      ENDIF

      GOTO 90000


C--------- MINTEST FUNCTION     EXTRSN


 1800 IF ( .NOT. GONLY ) THEN
         F = ZERO
         DO 1810 I = 1,N/2
            J  = 2*I
            W1 = X(J-1)
            W2 = X(J) - W1*W1
            W3 = ONE - W1
            F  = F + R100*W2*W2 + W3*W3
 1810    CONTINUE
      ENDIF

      IF ( .NOT. FONLY ) THEN

         DO 1820 I = 1,N
            IF ( MOD(I,2) .NE. 1 ) THEN
               W1   = X(I-1)
               G(I) = R200 * (X(I) - W1*W1 )
            ELSE
               W1   =   X(I)
               G(I) = - TWO * (ONE - W1 + W1 * (X(I+1)-W1*W1)*R200 )
            ENDIF

 1820    CONTINUE
      ENDIF

      GOTO 90000


C--------- MINTEST FUNCTION     MANCIN ( N, X, F, G, IFG, WORK, SIZE )


 2000 IF ( .NOT. FONLY ) THEN
         DO 2010 I =1,N
            G(I) = ZERO
 2010    CONTINUE
      ENDIF

      IF ( .NOT. GONLY ) THEN
         F = ZERO
      ENDIF

      IF ( (.NOT. FONLY .AND. 4*N .LE. SIZE) .OR.
     -     (      FONLY .AND.   N .LE. SIZE) ) THEN

         IR = 0
         IV = IR + N
         IC = IV + N
         IS = IC + N

         W1     = RD( BETA * N )
         W2     = RD(   N/2    )
         ALPHA2 = ALPHA - 2

         DO 2040 I = 1,N

            W3 = W1*X(I)
            W4 = RD(I)
            W5 = (W4 - W2)**GAMMA
            W6 = ZERO

C           COMPUTE THE SUM OF H(I,J) IN Y.

            DO 2020 J = 1,N
               W7 = X(J)
               IF ( I .NE. J ) THEN
                  W8         = SQRT( W7*W7 + W4/RD(J) )
                  IF ( .NOT. FONLY ) WORK(IV+J) = W8
                  W9         = LOG(W8)
                  W10        = SIN(W9)
                  IF ( .NOT. FONLY ) WORK(IS+J) = W10
                  W9         = COS(W9)
                  IF ( .NOT. FONLY ) WORK(IC+J) = W9
                  W6         = W6 + W8*( W10**ALPHA + W9**ALPHA )
               ENDIF
 2020       CONTINUE

            RI = W6 + W3 + W5

            IF ( .NOT. GONLY ) THEN
               F = F + RI*RI
               WORK(IR+I) = RI
            ENDIF

            IF ( .NOT. FONLY ) THEN
               DO 2030 K = 1,N
                  IF ( K .NE. I ) THEN
                     W5 = WORK(IS+K)
                     W6 = W5 ** ALPHA2
                     W8 = WORK(IC+K)
                     W9 = W8 ** ALPHA2
                     G(K) = G(K) + RI * X(K) * ( W6*W5*W5 + W9*W8*W8 +
     -                       RD(ALPHA) * W5* W8* (W6-W9))/WORK(IV+K)
                  ELSE
                     G(K) = G(K) + RI*W1
                  ENDIF
 2030          CONTINUE
            ENDIF
 2040    CONTINUE

         IF ( .NOT. FONLY ) THEN
            CALL ZZSCAL ( N, TWO, G, 1 )
         ENDIF

      ENDIF

      GOTO 90000

C---------        MINTEST FUNCTION     POWBSC
 2200 CONTINUE

C   CHECK FOR POSSIBLE OVERFLOW

      DO 2220 I = 1, N
         IF ( -X(I) .GT. BIGGST ) THEN
            RET = NOFG
            GOTO 90000
         ENDIF
 2220 CONTINUE

      IF ( .NOT. GONLY ) THEN
         F = ZERO

         W1 = EXP(-X(1))
         DO 2240 I = 1, N-1
            W2 = EXP(-X(I+1))
            W3 = W1 + W2 - R1PZ1
            W4 = R10000 * X(I) * X(I+1)  -  ONE
            W1 = W2

            F = F + W4*W4 + W3*W3
 2240    CONTINUE
      ENDIF

      IF ( .NOT. FONLY ) THEN
         W1   = EXP(-X(1))
         W2   = EXP(-X(2))
         W4   = TWO * R10000
         W5   = R10000*X(1)*X(2) - ONE
         W7   = W1 + W2 - R1PZ1
         G(1) = W4 * X(2) * W5  -  TWO * W1 * W7

         DO 2260 I = 2, N-1
            W3   = EXP(-X(I+1))
            W6   = R10000*X(I)*X(I+1) - ONE
            W8   = W2 + W3 - R1PZ1
            G(I) = W4 * ( X(I-1)*W5 + X(I+1)*W6 ) - TWO * W2 * (W7+W8)
            W2   = W3
            W5   = W6
            W7   = W8
 2260    CONTINUE

         G(N) = W4*X(N-1)*W5  -  TWO*W2*(EXP(-X(N-1)) + W2 - R1PZ1)
      ENDIF

      GOTO 90000

C--------- MINTEST FUNCTION   JENSMP ( N, X, F, G, IFG, NINT(FARG(1)))
C--------- NINT(FARG(1)) IS M1


 2300 X1 = X(1)
      X2 = X(2)

      IF ( .NOT. GONLY ) THEN
         F = ZERO
      ENDIF

      IF ( .NOT. FONLY ) THEN
         G1 = ZERO
         G2 = ZERO
      ENDIF

      DO 2310  I = 1, NINT(FARG(1))
         W3 = RD(I)
         IF ( MAX(W3*X1,W3*X2) .LE. BIGGST ) THEN
            W1 = EXP(W3 * X1)
            W2 = EXP(W3 * X2)
         ELSE
            RET = NOFG
            GOTO 90000
         ENDIF
         RI = TWO + TWO*W3 - W1 - W2

         IF ( .NOT. GONLY ) THEN
             F = F + RI*RI
         ENDIF

         IF ( .NOT. FONLY ) THEN
            W3 = RI*W3
            G1 = G1 + W1*W3
            G2 = G2 + W2*W3
         ENDIF
 2310 CONTINUE

      IF ( .NOT. FONLY ) THEN
         G(1) = -TWO * G1
         G(2) = -TWO * G2
      ENDIF

      GOTO 90000


C--------- MINTEST FUNCTION     FRDRTH


 2400 CONTINUE

      IF ( .NOT. GONLY ) THEN
         F = ZERO
      ENDIF

      W3 = ZERO
      W4 = ZERO

      DO 2450 I = 1,N-1

         W1 = X(I) + X(I+1)*(X(I+1)*(FIVE -  X(I+1)) - TWO) - R13
         W2 = X(I) + X(I+1)*(X(I+1)*( X(I+1)  + ONE) - R14) - R29

         IF ( .NOT. GONLY ) THEN
            F = F + W1*W1 + W2*W2
         ENDIF

         IF ( .NOT. FONLY ) THEN
            G(I) = TWO * ( W1 + W2
     -                   + W3 * (X(I)*(TEN-THREE*X(I))-TWO)
     -                   + W4 * (X(I)*(TWO+THREE*X(I))-R14)   )
            W3 = W1
            W4 = W2
         ENDIF

 2450 CONTINUE

      IF ( .NOT. FONLY ) THEN
         G(N) =  TWO * ( W3 * (X(N)*(TEN-THREE*X(N))-TWO)
     -                 + W4 * (X(N)*(TWO+THREE*X(N))-R14)   )
      ENDIF

      GOTO 90000


C--------- MINTEST FUNCTION     BROWNB


 2500 X1 = X(1)
      X2 = X(2)

      W1 = X1 - R1PD6
      W2 = X2 - R2PDM6
      W3 = X1*X2 - TWO

      IF ( .NOT. GONLY ) THEN
         F = W1*W1 + W2*W2 + W3*W3
      ENDIF

      IF ( .NOT. FONLY ) THEN
         G(1) = TWO * (W1 + X2*W3)
         G(2) = TWO * (W2 + X1*W3)
      ENDIF
      GOTO 90000

C--------- MINTEST FUNCTION     HILBRT ( N, X, F, G, IFG, FARG(1) )
C--------- FARG(1) IS FACTOR A.
C--------- FARG(2) IS DIAGONAL ELEMENT D.

 2700 IF ( .NOT. GONLY ) THEN
         F = ZERO
      ENDIF
      W5 = FARG(1)

      DO 2710 I=1,N
         W2 = ZERO
         I1 = I - 1
         DO 2720 J = 1,N
            IF ( I .EQ. J ) THEN
               W3 = FARG(2)
            ELSE
               W3 = ZERO
            ENDIF
            W2 = W2 + X(J)*( W3 + ONE/RD(I1+J) )
 2720    CONTINUE
         IF ( .NOT. FONLY ) THEN
CCC            G(I) = W2
            G(I) = TWO * W5 * W2
         ENDIF
         IF ( .NOT. GONLY ) THEN
            F = F + X(I)*W2
         ENDIF
 2710 CONTINUE

      IF ( .NOT. GONLY ) THEN
CCC         F= F / TWO
         F= W5 * F
      ENDIF
      GOTO 90000

C     ---------        MINTEST FUNCTION     ZANGW1
 2800 X1 = X(1)
      X2 = X(2)

      IF (.NOT. GONLY) THEN
         F = ( R16*(X1*X1 + X2*X2) - EIGHT*X1*(X2+SEVEN)
     -                               -  R256*X2 +  R991  ) / R15
      ENDIF

      IF (.NOT. FONLY) THEN
         G(1) = (R32*X1 - EIGHT*X2 - R56 ) / R15
         G(2) = (R32*X2 - EIGHT*X1 - R256) / R15
      ENDIF
            GOTO 90000
C      ---------      MINTEST FUNCTION     HIMLN3
 2900 X1 = X(1)
      X2 = X(2)

      W1 = X1 * X1
      W2 = X1 * W1
      W3 = X2 * X2

      IF ( .NOT. GONLY ) THEN
         F = W2 + W3 - THREE*X1 - TWO*X2 + TWO
      ENDIF

      IF ( .NOT. FONLY ) THEN
         G(1) = THREE*W1 - THREE
         G(2) = TWO*X2 - TWO
      ENDIF

      GOTO 90000


C--------- MINTEST FUNCTION     BEAL58


 3000 X1 = X(1)
      X2 = X(2)

      W1 = X2 * X2

      W3 = ONE - X2
      W2 = ONE - W1
      W1 = ONE - W1*X2

      W4 = R1P5   - X1*W3
      W5 = R2P25  - X1*W2
      W6 = R2P625 - X1*W1

      IF ( .NOT. GONLY ) THEN
         F = W4*W4 + W5*W5 + W6*W6
      ENDIF

      IF ( .NOT. FONLY ) THEN
         G(1) = -TWO * ( W4*W3 + W5*W2 + W6*W1 )
         G(2) =  TWO * X1 * ( W4 + (TWO*W5 + THREE*W6*X2) * X2 )
      ENDIF

      GOTO 90000


C--------- MINTEST FUNCTION     ENGVL1


 3100 CONTINUE

      IF ( .NOT. GONLY ) THEN
         F = ZERO
      ENDIF

      W1 = X(1)**2
      W3 = ZERO

      DO 3150 I = 1,N-1

         W2 = X(I+1)**2
         W4 = W1 + W2

         IF ( .NOT. GONLY ) THEN
            F = F + W4**2 - FOUR*X(I) + THREE
         ENDIF

         IF ( .NOT. FONLY ) THEN
            G(I) = FOUR * (X(I) * (W3+W4) - ONE)
            W3 = W4
         ENDIF

         W1 = W2
 3150 CONTINUE

      G(N) = FOUR * (X(N) * W3)
      GOTO 90000

C---------      DIXON
 3200 CONTINUE
      IF ( .NOT. GONLY ) THEN
         F = (ONE-X(1))**2 + (ONE-X(10))**2
         DO 3250 I = 2,9
            F = F + (X(I)-X(I+1))**2
 3250    CONTINUE
      ENDIF

      IF ( .NOT. FONLY ) THEN
         G(1) = -TWO*(ONE -X(1))
         G(2) =  TWO*(X(2)-X(3))
         DO 3275 I = 3,9
            G(I) = TWO*(TWO*X(I) - X(I-1) - X(I+1) )
 3275    CONTINUE
         G(10) = -TWO*(ONE + X(9) - TWO*X(10) )
      ENDIF
      GOTO 90000

C---------      ZANGW2
 3300 X1 = X(1)
      X2 = X(2)
      X3 = X(3)

      W1 =  X1 - X2 + X3
      W2 = -X1 + X2 + X3
      W3 =  X1 + X2 - X3

      IF ( .NOT. GONLY ) THEN
         F = W1*W1 + W2*W2 + W3*W3
      ENDIF

      IF ( .NOT. FONLY ) THEN
         G(1) = TWO * ( THREE*X1 -       X2 -       X3 )
         G(2) = TWO * (      -X1 + THREE*X2 -       X3 )
         G(3) = TWO * (      -X1 -       X2 + THREE*X3 )
      ENDIF

      GOTO 90000


C---------       FUNCTION     HIMM25
 3800 CONTINUE

      W1 = X(1) - FIVE
      W2 = X(2) - SIX

      IF ( .NOT. GONLY ) THEN
         F = FOUR*W1*W1 + W2*W2
      ENDIF

      IF ( .NOT. FONLY ) THEN
        G(1) = EIGHT*W1
        G(2) =   TWO*W2
      ENDIF
      GOTO 90000

C---------                  QUARTC
 3900            CONTINUE
      IF ( .NOT. GONLY ) THEN
         F = ZERO
         DO 3930 I = 1,N
            F = F + (X(I) - RD(I))**4
 3930    CONTINUE
      ENDIF

      IF ( .NOT. FONLY ) THEN
         DO 3960 I = 1,N
            G(I) = FOUR*(X(I) - RD(I))**3
 3960    CONTINUE
      ENDIF

      GOTO 90000

C----                    GENRSN
C     GILL AND MURRAY'S GENERALIZED ROSENBROCK FUNCTION.
C     (GILM79 #7, AND TOI83B #10 AND #48).

 8500 CONTINUE

      IF ( .NOT. GONLY ) THEN
         F = ZERO

         DO 8510 I = 2, N
            F = F  +  R100*( X(I) - X(I-1)**2 )**2  +  ( X(I) - ONE )**2
 8510    CONTINUE
      ENDIF

      IF ( .NOT. FONLY ) THEN
         W1   =  X(2) - X(1)**2
         G(1) = -R400 * X(1) * W1

         DO 8520 I = 2, N-1
            W2   = X(I+1) - X(I)**2
            G(I) = R200 * W1  +  TWO * (X(I) - ONE) - R400 * X(I) * W2
            W1   = W2
 8520    CONTINUE

         G(N) =  R200 * W1  +  TWO * ( X(N) - ONE )
      ENDIF

      GOTO 90000

C----                     GREKAR (FUNCTION TRIDIAG FROM BRE73 PP 142-143.)

 8600 CONTINUE

      IF ( .NOT. GONLY ) THEN
         F = (X(1) - X(2))*X(1) + (-X(N-1) + TWO*X(N))*X(N) - TWO*X(1)

         DO 8610 I = 2, N-1
            F = F  +  ( -X(I-1) + TWO*X(I) - X(I+1) ) * X(I)
 8610    CONTINUE
      ENDIF

      IF ( .NOT. FONLY ) THEN
         G(1) = TWO * ( X(1) - X(2) - ONE )

         DO 8620 I = 2, N-1
            G(I) = -TWO * ( X(I-1) - TWO*X(I) + X(I+1) )
 8620    CONTINUE

         G(N) = -TWO * X(N-1)  +  FOUR * X(N)
      ENDIF

      GOTO 90000

C----                    TDQUAD
C     TOINT'S DIAGONAL QUADRATIC FUNCTION  (PROBLEM #22 IN TOI83B).

 8700 CONTINUE

      IF ( .NOT. GONLY ) THEN
         F = ZERO

         DO 8710 I = 1, N-2
            F = F  +  R100 * ( X(I+1)**2 + X(I+2)**2 )  +  X(I)**2
 8710    CONTINUE
      ENDIF

      IF ( .NOT. FONLY ) THEN
         DO 8720 I = 1, N
            G(I) = ZERO
 8720    CONTINUE

         DO 8730 I = 1, N-2
            G(I)   = G(I)    +  TWO  * X(I)
            G(I+1) = G(I+1)  +  R200 * X(I+1)
            G(I+2) = G(I+2)  +  R200 * X(I+2)
 8730    CONTINUE
      ENDIF

      GOTO 90000

C----                    TOIN2
C     TOINT'S SECOND PROBLEM  (PROBLEM #2 IN TOI83B).

 8800 CONTINUE

      IF ( .NOT. GONLY ) THEN
         F = X(1)**2 + ( (X(1) - X(2))**2 + (X(2) - X(3))**2 ) / TWO
      ENDIF

      IF ( .NOT. FONLY ) THEN
         G(1) = THREE * X(1)  -  X(2)
         G(2) = -X(1)  +  TWO * X(2)  -  X(3)
         G(3) = X(3) - X(2)
      ENDIF

      GOTO 90000

C----                    TOIN4
C     TOINT'S FOURTH PROBLEM  (PROBLEM #4 IN TOI83B).

 8900 CONTINUE

      X1 = X(1)
      X2 = X(2)
      X3 = X(3)
      X4 = X(4)

      IF ( .NOT. GONLY ) THEN
         F = X1**2 + X2**2 + TWO + ( (X1 - X2)**2 + (X3 - X4)**2 ) / TWO
      ENDIF

      IF ( .NOT. FONLY ) THEN
         G(1) = THREE * X1  -  X2
         G(2) = -X1  +  THREE * X2
         G(3) = X3 - X4
         G(4) = X4 - X3
      ENDIF

      GOTO 90000

C## E X I T
90000       STATUS(7) =  RET
91000       GFIRST    = .FALSE.
      RETURN

C## F O R M A T S:  NONE ARE DEFINED.
C##                 E N D         OF ZZFNS2.
                    END
