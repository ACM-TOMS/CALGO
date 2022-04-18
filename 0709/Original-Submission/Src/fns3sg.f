      SUBROUTINE ZZFNS3 ( N, X, F, G, WORK, SIZE, FONLY, GONLY,
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
C>RCS $HEADER: FNS3.F,V 1.1 91/11/20 10:52:50 BUCKLEY EXP $
C>RCS $LOG:     FNS3.F,V $
C>RCSREVISION 1.1  91/11/20  10:52:50  BUCKLEY
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
C## E N T R Y   P O I N T S:  ZZFNS3    THE NATURAL ENTRY POINT.
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

      INTEGER     I1, I2, F1, I, IDUMMY, J, K, NOVER2, IEG, IY, RET, L

      LOGICAL           EVEN,   GFIRST, ERROR

      REAL              ZZMPAR, HUGE
C!!!! DOUBLE PRECISION  ZZMPAR, HUGE

C--------- VARIABLES FOR THE TEST FUNCTIONS.

      REAL              X1, X2, X3, X4, X5, X6, S1, T1, T2, T3, T4
C!!!! DOUBLE PRECISION  X1, X2, X3, X4, X5, X6, S1, T1, T2, T3, T4

      REAL              X7,     X8,     X9,    X10,   X11
C!!!! DOUBLE PRECISION  X7,     X8,     X9,    X10,   X11

      REAL              G1,     G2,     G3,    G4,    G5,    G6
C!!!! DOUBLE PRECISION  G1,     G2,     G3,    G4,    G5,    G6

      REAL              G7,     G8,     G9,    G10,   G11
C!!!! DOUBLE PRECISION  G7,     G8,     G9,    G10,   G11

      REAL              W1,     W2,     W3,    W4,    W5,    W6
C!!!! DOUBLE PRECISION  W1,     W2,     W3,    W4,    W5,    W6

      REAL              W7,     W8,     W9
C!!!! DOUBLE PRECISION  W7,     W8,     W9

      REAL               R,      S,      T,    BIGGST, SMLLST
C!!!! DOUBLE PRECISION   R,      S,      T,    BIGGST, SMLLST

      REAL              RI,    RK ,    SK,    TI
C!!!! DOUBLE PRECISION  RI,    RK ,    SK,    TI

      REAL              XI,     XK,     YI,    PI,     U
C!!!! DOUBLE PRECISION  XI,     XK,     YI,    PI,     U

      REAL             R2P,     RD,    TPI,  TPIS
C!!!! DOUBLE PRECISION R2P,     RD,    TPI,  TPIS

      REAL              KAP1, KAP2,    RF1
C!!!! DOUBLE PRECISION  KAP1, KAP2,    RF1

C--------- DATA ARRAYS FOR FUNCTIONS

      REAL              AL (50),    ARGASY (15)
C!!!! DOUBLE PRECISION  AL (50),    ARGASY (15)

      REAL              HIM32A (7),     HIM32B (7)
C!!!! DOUBLE PRECISION  HIM32A (7),     HIM32B (7)

      REAL              KOWOSU (11),    KOWOSY (11)
C!!!! DOUBLE PRECISION  KOWOSU (11),    KOWOSY (11)

      REAL             ORBETA (33), OD (33), MEY (16)
C!!!! DOUBLE PRECISION ORBETA (33), OD (33), MEY (16)

      REAL             OSB1Y (33), OSB2Y (65)
C!!!! DOUBLE PRECISION OSB1Y (33), OSB2Y (65)

      INTEGER     A (50), B (56)

C## S A V E:

      SAVE  GFIRST, PI,  TPI, TPIS, R2P, BIGGST, SMLLST
      SAVE  HUGE
      SAVE   OK, ABORT, LIMIT, NOF, NOG, NOFG

C--------- SAVE DATA ARRAYS FOR THE TEST FUNCTIONS.

      SAVE  ARGASY, AL, HIM32A, HIM32B, KOWOSU, KOWOSY
      SAVE  ORBETA, MEY   , OD, OSB1Y , OSB2Y  , A     , B

C## E Q U I V A L E N C E S: NONE ARE DEFINED.
C## C O M M O N:             NONE IS DEFINED.
C## D A T A:

      DATA  GFIRST/.TRUE./

C--------- DATA FOR FUNCTION     ARGAUS

      DATA  ARGASY( 1), ARGASY( 2), ARGASY( 3), ARGASY( 4), ARGASY( 5)
     -     / 9.000 D-4,  4.400 D-3,  1.750 D-2,  5.400 D-2,  1.295 D-1 /

      DATA  ARGASY( 6), ARGASY( 7), ARGASY( 8), ARGASY( 9), ARGASY(10)
     -     / 2.420 D-1,  3.521 D-1,  3.989 D-1,  3.521 D-1,  2.420 D-1 /

      DATA  ARGASY(11), ARGASY(12), ARGASY(13), ARGASY(14), ARGASY(15)
     -     / 1.295 D-1,  5.400 D-2,  1.750 D-2,  4.400 D-3,  9.000 D-4 /


C--------- DATA FOR FUNCTION     HIMM32

      DATA HIM32A(1),  HIM32A(2),  HIM32A(3),  HIM32A(4)
     -     /  0.0D0,      4.28D-4,    1.0D-3,     1.61D-3    /

      DATA HIM32A(5),  HIM32A(6),  HIM32A(7)
     -     /  2.09D-3,    3.48D-3,    5.25D-3                /

      DATA HIM32B(1),  HIM32B(2),  HIM32B(3),  HIM32B(4)
     -     /  7.391D0,    1.118D1,    1.644D1,    1.62D1     /

      DATA HIM32B(5),  HIM32B(6),  HIM32B(7)
     -     /  2.22D1,     2.402D1,    3.132D1                /

C--------- DATA FOR FUNCTION     KOWOSB

      DATA KOWOSU(1),  KOWOSU(2),  KOWOSU(3),  KOWOSU(4)
     -     /  4.0D0,      2.0D0,      1.0D0,      0.5D0      /

      DATA KOWOSU(5),  KOWOSU(6),  KOWOSU(7),  KOWOSU(8)
     -     /  0.25D0,     0.167D0,    0.125D0,    0.1D0      /

      DATA KOWOSU(9),  KOWOSU(10),  KOWOSU(11)
     -     /  0.0833D0,   0.0714D0,    0.0625D0              /

      DATA KOWOSY(1),  KOWOSY(2),  KOWOSY(3),  KOWOSY(4)
     -     /  0.1957D0,   0.1947D0,   0.1735D0,   0.1600D0  /

      DATA KOWOSY(5),  KOWOSY(6),  KOWOSY(7),  KOWOSY(8)
     -     /  0.0844D0,   0.0627D0,   0.0456D0,   0.0342D0  /

      DATA KOWOSY(9),  KOWOSY(10),  KOWOSY(11)
     -     /  0.0323D0,   0.0235D0,    0.0246D0             /

C--------- DATA FOR FUNCTION        MEYER

      DATA   MEY(1),   MEY(2),   MEY(3),   MEY(4),  MEY(5),  MEY(6)
     -     /3.478D4,  2.861D4,  2.365D4,  1.963D4, 1.637D4, 1.372D4/

      DATA   MEY(7),   MEY(8),   MEY(9),   MEY(10),  MEY(11), MEY(12)
     -     /1.154D4,  9.744D3,  8.261D3,  7.030D3, 6.005D3, 5.147D3/

      DATA    MEY(13),   MEY(14),   MEY(15),   MEY(16)
     -     /4.427D3,  3.820D3,  3.307D3,  2.872D3/

C--------- DATA FOR FUNCTION         ORTOIT

      DATA ORBETA(1),ORBETA(2),ORBETA(3),ORBETA(4),ORBETA(5)
     -   /1.0D0, 1.5D0, 1.0D0, 0.1D0, 1.5D0/

      DATA ORBETA(6),ORBETA(7),ORBETA(8),ORBETA(9),ORBETA(10)
     -   /2.0D0, 1.0D0, 1.5D0, 3.0D0, 2.0D0/

      DATA ORBETA(11),ORBETA(12),ORBETA(13),ORBETA(14),ORBETA(15)
     -   /1.0D0, 3.0D0, 0.1D0, 1.5D0, 0.15D0/

      DATA ORBETA(16),ORBETA(17),ORBETA(18),ORBETA(19),ORBETA(20)
     -   /2.0D0, 1.0D0, 0.1D0, 3.0D0, 0.1D0/

      DATA ORBETA(21),ORBETA(22),ORBETA(23),ORBETA(24),ORBETA(25)
     -   /1.2D0, 1.0D0, 0.1D0, 2.0D0, 1.2D0/

      DATA ORBETA(26),ORBETA(27),ORBETA(28),ORBETA(29),ORBETA(30)
     -   /3.0D0, 1.5D0, 3.0D0, 2.0D0, 1.0D0/

      DATA ORBETA(31),ORBETA(32),ORBETA(33)
     -   /1.2D0, 2.0D0, 1.0D0/

      DATA OD(1), OD(2), OD(3), OD(4), OD(5), OD(6)
     -   /  5.0D0,5.0D0,5.0D0,2.5D0,6.0D0,6.0D0 /

      DATA OD(7), OD(8), OD(9), OD(10), OD(11), OD(12)
     -   /  5.0D0,6.0D0,10.0D0,6.0D0,5.0D0,9.0D0 /

      DATA OD(13), OD(14), OD(15), OD(16), OD(17), OD(18)
     -   /  2.0D0,7.0D0,2.5D0,6.0D0,5.0D0,2.0D0 /

      DATA OD(19), OD(20), OD(21), OD(22), OD(23), OD(24)
     -   /  9.0D0,2.0D0,5.0D0,5.0D0,2.5D0,5.0D0 /

      DATA OD(25), OD(26), OD(27), OD(28), OD(29), OD(30)
     -   /  6.0D0,10.0D0,7.0D0,10.0D0,6.0D0,5.0D0 /

      DATA OD(31), OD(32), OD(33)
     -   /  4.0D0,4.0D0,4.0D0 /

      DATA A(1),A(2),A(3),A(4),A(5),A(6),A(7),A(8),A(9)
     -   /  -31,-1,-2,-4,-6,-8,-10,-12,+11 /
      DATA A(10),A(11),A(12),A(13),A(14),A(15),A(16),A(17),A(18)
     -   /  +13,-14,-16,+9,-18,+5,+20,-21,-19 /
      DATA A(19),A(20),A(21),A(22),A(23),A(24),A(25),A(26),A(27)
     -   /  -23,+7,-25,-28,-29,-32,+3,-33,-35 /
      DATA A(28),A(29),A(30),A(31),A(32),A(33),A(34),A(35),A(36)
     -   /  -36,+30,-37,+38,-39,-40,-41,-44,-46 /
      DATA A(37),A(38),A(39),A(40),A(41),A(42),A(43),A(44),A(45)
     -   /  +42,+45,+48,-50,+26,+34,-43,+15,+17 /
      DATA A(46),A(47),A(48),A(49),A(50)
     -   /  +24,-47,-49,-22,-27 /

      DATA B(1),B(2),B(3),B(4),B(5),B(6),B(7),B(8),B(9)
     -   /  -1,+2,-3,+4,-5,+6,-7,+8,-9 /
      DATA B(10),B(11),B(12),B(13),B(14),B(15),B(16),B(17),B(18)
     -   /  +10,-11,+12,-13,+14,-15,+16,-17,+18 /
      DATA B(19),B(20),B(21),B(22),B(23),B(24),B(25),B(26),B(27)
     -   /  -19,-20,0,+22,+23,-24,+25,-26,+27 /
      DATA B(28),B(29),B(30),B(31),B(32),B(33),B(34),B(35),B(36)
     -   /  -28,+29,-30,+31,-32,+33,-34,-35,+21 /
      DATA B(37),B(38),B(39),B(40),B(41),B(42),B(43),B(44),B(45),B(46)
     -   /  -36,+37,-38,-39,-40,+41,-42,+43,+44,-50 /
      DATA B(47),B(48),B(49),B(50),B(51),B(52),B(53),B(54),B(55),B(56)
     -   /  +45,+46,-47,-48,-49,0,0,0,0,0 /

C-----  ARRAY AL IS USED BY MINTEST FUNCTIONS ORTOIT AND CHNRSN

      DATA   AL(1), AL(2), AL(3), AL(4), AL(5), AL(6), AL(7), AL(8)
     -   /  1.25D0,1.40D0,2.40D0,1.40D0,1.75D0,1.20D0,2.25D0,1.20D0/

      DATA AL(9) ,AL(10),AL(11),AL(12),AL(13),AL(14),AL(15),AL(16)
     -   /  1.00D0,1.10D0,1.50D0,1.60D0,1.25D0,1.25D0,1.20D0,1.20D0/

      DATA AL(17),AL(18),AL(19),AL(20),AL(21),AL(22),AL(23),AL(24)
     -   /  1.40D0,0.50D0,0.50D0,1.25D0,1.80D0,0.75D0,1.25D0,1.40D0/

      DATA AL(25),AL(26),AL(27),AL(28),AL(29),AL(30)
     -   /  1.60D0,2.00D0,1.00D0,1.60D0,1.25D0,2.75D0/

      DATA AL(31),AL(32),AL(33),AL(34),AL(35),AL(36),AL(37),AL(38)
     -   /  1.25D0,1.25D0,1.25D0,3.00D0,1.50D0,2.00D0,1.25D0,1.40D0/

      DATA AL(39),AL(40),AL(41),AL(42),AL(43),AL(44),AL(45),AL(46)
     -   /  1.80D0,1.50D0,2.20D0,1.40D0,1.50D0,1.25D0,2.00D0,1.50D0/

      DATA AL(47),AL(48),AL(49),AL(50)
     -   /  1.25D0,1.40D0,0.60D0,1.50D0/

C--------- DATA FOR FUNCTION        OSBRN1

      DATA    OSB1Y(1),  OSB1Y(2),  OSB1Y(3),  OSB1Y(4),  OSB1Y(5)
     -/.844D0,  .908D0,  .932D0,  .936D0,  .925D0/

      DATA    OSB1Y(6),  OSB1Y(7),  OSB1Y(8),  OSB1Y(9),  OSB1Y(10)
     -/.908D0,  .881D0,  .850D0,  .818D0,  .784D0/

      DATA    OSB1Y(11),  OSB1Y(12),  OSB1Y(13),  OSB1Y(14),  OSB1Y(15)
     -/.751D0,  .718D0,  .685D0,  .658D0,  .628D0/

      DATA    OSB1Y(16),  OSB1Y(17),  OSB1Y(18),  OSB1Y(19),  OSB1Y(20)
     -/.603D0,  .580D0,  .558D0,  .538D0,  .522D0/

      DATA    OSB1Y(21),  OSB1Y(22),  OSB1Y(23),  OSB1Y(24),  OSB1Y(25)
     -/.506D0,  .490D0,  .478D0,  .467D0,  .457D0/

      DATA    OSB1Y(26), OSB1Y(27),  OSB1Y(28),  OSB1Y(29),  OSB1Y(30)
     -/.448D0,  .438D0,  .431D0,  .424D0,  .420D0/

      DATA    OSB1Y(31),  OSB1Y(32),  OSB1Y(33)
     -/.414D0, .411D0, .406D0/

C--------- DATA FOR FUNCTION        OSBRN2

      DATA    OSB2Y(1),  OSB2Y(2),  OSB2Y(3),  OSB2Y(4),  OSB2Y(5)
     -/1.366D0,  1.191D0,  1.112D0,  1.013D0,  .991D0/

      DATA    OSB2Y(6),  OSB2Y(7),  OSB2Y(8),  OSB2Y(9),  OSB2Y(10)
     -/.885D0,  .831D0,  .847D0, .786D0,  .725D0/

      DATA    OSB2Y(11),  OSB2Y(12),  OSB2Y(13),  OSB2Y(14),  OSB2Y(15)
     -/.746D0,  .679D0,  .608D0,  .655D0,  .616D0/

      DATA    OSB2Y(16),  OSB2Y(17),  OSB2Y(18),  OSB2Y(19),  OSB2Y(20)
     -/.606D0,  .602D0,  .626D0,  .651D0,  .724D0/

      DATA    OSB2Y(21),  OSB2Y(22),  OSB2Y(23),  OSB2Y(24),  OSB2Y(25)
     -/.649D0,  .649D0,  .694D0,  .644D0,  .624D0/

      DATA    OSB2Y(26),  OSB2Y(27),  OSB2Y(28),  OSB2Y(29),  OSB2Y(30)
     -/.661D0,  .612D0,  .558D0,  .533D0,  .495D0/

      DATA    OSB2Y(31),  OSB2Y(32),  OSB2Y(33),  OSB2Y(34),  OSB2Y(35)
     -/.50D0,  .423D0,  .395D0,  .375D0,  .372D0/

      DATA    OSB2Y(36),  OSB2Y(37),  OSB2Y(38),  OSB2Y(39),  OSB2Y(40)
     -/.391D0,  .396D0,  .405D0,  .428D0,  .429D0/

      DATA    OSB2Y(41),  OSB2Y(42),  OSB2Y(43),  OSB2Y(44),  OSB2Y(45)
     -/.523D0,  .562D0,  .607D0,  .653D0,  .672D0/

      DATA    OSB2Y(46),  OSB2Y(47),  OSB2Y(48),  OSB2Y(49),  OSB2Y(50)
     -/.708D0,  .633D0,  .668D0,  .645D0, .632D0/

      DATA    OSB2Y(51),  OSB2Y(52),  OSB2Y(53),  OSB2Y(54),  OSB2Y(55)
     -/.591D0,  .559D0,  .597D0,  .625D0,  .739D0/

      DATA    OSB2Y(56),  OSB2Y(57),  OSB2Y(58),  OSB2Y(59),  OSB2Y(60)
     -/.710D0,  .729D0,  .720D0,  .636D0,  .581D0/

      DATA    OSB2Y(61),  OSB2Y(62),  OSB2Y(63),  OSB2Y(64),  OSB2Y(65)
     -/.428D0,  .292D0,  .162D0,  .098D0,  .054D0/

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
     -     1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1900, 1000,
     -     1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000,
     -     1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 4000,
     -     4100, 4200, 4300, 4400, 4500, 4600, 4700, 4800, 4900, 1000,
     -     1000, 5200, 5300, 5400, 5500, 5600, 5700, 5800, 5900, 6000,
     -     6100, 6200, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000,
     -     1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000,
     -     1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000
     -    ) FUNCNO
 1000 STATUS(7) = ABORT
      GOTO 91000


C--------- MINTEST FUNCTION     CHNRSN


 1900 IF ( .NOT. GONLY ) THEN

         F  = ZERO
         W1 = X(1)

         DO 1910 I=2,N
            W2 = X(I)
            W3 = W1 - W2*W2
            W4 = ONE - W2
            F = F + FOUR * AL(I) * W3*W3 +  W4*W4
            W1 = W2
 1910    CONTINUE

      ENDIF

      IF ( .NOT. FONLY ) THEN

         W2   = X(2)

         W1   = EIGHT * AL(2) * ( X(1) - W2*W2 )
         G(1) = W1

         DO 1920 I = 3,N
            W4     = X(I)
            W3     = EIGHT * AL(I) * (W2 - W4*W4)
            G(I-1) = TWO*(-W1*W2 - (ONE - W2))  +  W3
            W1     = W3
            W2     = W4
 1920    CONTINUE

         G(N) = -TWO * (W1*W4 + (ONE - W4))

      ENDIF

      GOTO 90000


C---------                 TRIGTO
 4000            CONTINUE
      IF ( .NOT. GONLY ) THEN
         F = ZERO
      ENDIF
      IF ( .NOT. FONLY ) THEN
         DO 4010 K=1,N
            G(K)=ZERO
 4010    CONTINUE
      ENDIF

      DO 4070 I=1,N
         S  = RD(I)
         S1 = ONE + S/TEN
         DO 4020 J=1,I-1
             T    = RD(J)
             T1   = ONE + T/TEN
             T2   = RD(5*(1+MOD(I,5)+MOD(J,5)))
             T3   = (S+T)/TEN
             T4   = S1*X(I) + T1*X(J) + T3
             IF ( .NOT. GONLY ) THEN
                F    = F + T2*SIN(T4)
             ENDIF
             IF ( .NOT. FONLY ) THEN
                G(I) = G(I) + S1*T2*COS(T4)
             ENDIF
 4020    CONTINUE
         IF ( .NOT. FONLY ) THEN
           DO 4040 J=I+1,N
             T    = RD(J)
             T1   = ONE + T/TEN
             T2   = RD(5*(1+MOD(I,5)+MOD(J,5)))
             T3   = (S+T)/TEN
             T4   = S1*X(I) + T1*X(J) + T3
             G(I) = G(I) + T2*S1*COS(T4)
 4040      CONTINUE
         ENDIF
 4070 CONTINUE
      GOTO 90000

C---------                  MINTEST FUNCTION     HIMM28
 4100 X1 = X(1)
      X2 = X(2)

      W1 = X1*X1 + X2 - R11
      W2 = X2*X2 + X1 - SEVEN

      IF ( .NOT. GONLY ) THEN
         F = W1*W1 + W2*W2
      ENDIF

      IF ( .NOT. FONLY ) THEN
         G(1) = FOUR*X1*W1 + TWO*W2
         G(2) = FOUR*X2*W2 + TWO*W1
      ENDIF

      GOTO 90000


C--------- MINTEST FUNCTION   ORTOIT ( N, X, F, G, IFG, NINT(FARG(1)),
C                                                  WORK, SIZE          )
C
C--------- NINT(FARG(1)) SELECTS QOR, GOQ OR PSP - TOINT.

 4200 F1 = NINT( FARG( 1 ) )

      IF ( .NOT. GONLY ) THEN
         F = ZERO
      ENDIF

      IF ( SIZE .LT. 83 ) THEN
         IF ( .NOT. FONLY ) THEN
            DO 4201 I = 1,N
               G(I) = ZERO
 4201       CONTINUE
         ENDIF
         GOTO 90000
      ENDIF

C-----SET POINTERS FOR WORK ARRAYS.

      IY  = 0
      IEG = IY + 33

      DO 4202 I = 1, 33
         WORK(IY+I) = OD(I)
 4202 CONTINUE

C *** ADD (SUBTRACT) X(K) TO (FROM) THE PROPER Y(I).

      J = 1
      L = 1
      DO 4206 K = 1, 50

         W2   = X(K)
         W1   = AL(K)

         I    = A(K)
         WORK(IY+J) = WORK(IY+J) - X(ABS(I))
         IF ( I .LT. 0 ) J = J + 1

         I    = B(K)
         IF ( I .NE. 0 ) WORK(IY+L) = WORK(IY+L) + X(ABS(I))
         IF ( I .LE. 0 ) L = L + 1

C        ****** AS APPROPRIATE ADD THE FIRST SUM TO F AND THE FIRST
C        TERM TO G(K).  THE SUM AND TERMS INVOLVE ALPHA(K) AND X(K).

         GO TO (4203, 4204, 4205), F1

C        ****** QORTOI ONLY.

 4203    IF ( .NOT. GONLY ) THEN
            F    = F + W1 * W2 * W2
         ENDIF

         IF ( .NOT. FONLY ) THEN
            G(K) =     W1 * W2
         ENDIF

         GO TO 4206

C        ****** GORTOI ONLY.

 4204    IF ( W2 .GE. ZERO ) THEN
            W3 = ONE + W2

            IF ( .NOT. GONLY ) THEN
               F    = F + W1 * W2 * LOG(W3)
            ENDIF

            IF ( .NOT. FONLY ) THEN
               G(K) =     W1 * ( (W2 / W3) + LOG(W3) )
            ENDIF

         ELSE
            W3 = ONE - W2

            IF ( .NOT. GONLY ) THEN
               F    = F - W1 * W2 * LOG(W3)
            ENDIF

            IF ( .NOT. FONLY ) THEN
               G(K) =     W1 * ( (W2 / W3) - LOG(W3) )
            ENDIF

         ENDIF
         GO TO 4206

C        ****** PSPTOI ONLY.

 4205    W2 = W2 - FIVE

         IF ( .NOT. GONLY ) THEN
            F    = F + W1 * W2 * W2
         ENDIF

         IF ( .NOT. FONLY ) THEN
            G(K) =     W1 * W2 * TWO
         ENDIF

 4206 CONTINUE

C *** FOR THE FUNCTION VALUE, ADD THE SUM INVOLVING Y(I).

      IF ( .NOT. GONLY ) THEN

         DO 4210 I = 1, 33

            W3 = WORK(IY+I)
            W1 = W3 * W3
            W2 = ORBETA(I)

            GO TO (4207, 4208, 4209), F1

C           ********* QORTOI ONLY.

 4207       F = F + W2 * W1
            GO TO 4210

C           ********* GORTOI ONLY.

 4208       IF ( W3 .GE. ZERO ) THEN
               F = F + W2 * W1 * LOG(ONE + W3)
            ELSE
               F = F + W2 * W1
            ENDIF
            GO TO 4210

C           ********* PSPTOI ONLY.

 4209      IF ( W3 .GE. TENTH ) THEN
               F = F + W2 / W3
            ELSE
               F = F + W2 * R20 * (ONE - FIVE * W3)
            ENDIF

 4210    CONTINUE

      ENDIF

C *** FOR THE GRADIENTS ADD THE SUMS INVOLVING Y(I).

      IF ( .NOT. FONLY ) THEN

         J = 1
         L = 1
         DO 4214 K = 1, 50

            GO TO (4211, 4212, 4213), F1

C           ********* QORTOI ONLY.

 4211       I  = A(K)
            W1 =    - ORBETA(J) * WORK(IY+J)
            G(ABS(I)) =   G(ABS(I)) + W1
            IF ( I .LT. 0 ) J = J + 1

            I  = B(K)
            W1 =      ORBETA(L) * WORK(IY+L)
            IF ( I .NE. 0 ) G(ABS(I)) =   G(ABS(I)) + W1
            IF ( I .LE. 0 ) L = L + 1

            GO TO 4214

C           ********* GORTOI ONLY.

 4212       I  = A(K)
            W2 = WORK(IY+J)

            IF ( W2 .GE. ZERO ) THEN
               W1 = -ORBETA(J) * W2 * ( (TWO * LOG(ONE + W2)) +
     -                        (W2 / (ONE + W2)) )
            ELSE
               W1 = -ORBETA(J) * TWO * W2
            ENDIF

            G(ABS(I)) =   G(ABS(I)) + W1
            IF ( I .LT. 0 ) J = J + 1

            I  = B(K)
            IF ( I .NE. 0 ) THEN

               W2 = WORK(IY+L)

               IF ( W2 .GE. ZERO ) THEN
                  W1 = ORBETA(L) * W2 * ( (TWO * LOG(ONE + W2)) +
     -                            (W2 / (ONE + W2)) )
               ELSE
                  W1 = ORBETA(L) * TWO * W2
               ENDIF
               G(ABS(I)) = G(ABS(I)) + W1
            ENDIF
            IF ( I .LE. 0 ) L = L + 1

            GO TO 4214

C           ********* PSPTOI ONLY.  (NOTE THAT THE MINUS SIGNS ON
C              THE SUM AND FUNCTION HAVE CANCELED EACH OTHER.)

 4213       I  = A(K)
            W2 = WORK(IY+J)

            IF ( W2 .GE. TENTH ) THEN
               W1 = ORBETA(J) / ( W2 * W2 )
            ELSE
               W1 = ORBETA(J) * R100
            ENDIF

            G(ABS(I)) =   G(ABS(I)) + W1
            IF ( I .LT. 0 ) J = J + 1

            I  = B(K)
            IF ( I .NE. 0 ) THEN

               W2 = WORK(IY+L)

C              ********* NOTE THAT THE MINUS SIGN IS DUE TO
C                        THE FUNCTIONS USED.

               IF ( W2 .GE. TENTH ) THEN
                  W1 = - ORBETA(L) / ( W2 * W2 )
               ELSE
                  W1 = - ORBETA(L) * R100
               ENDIF
               G(ABS(I)) =   G(ABS(I)) + W1
            ENDIF
            IF ( I .LE. 0 ) L = L + 1

 4214   CONTINUE

        DO 4229 K = 1,50

        GOTO ( 4222, 4224, 4226 ), F1

 4222      G(K) = G(K) * TWO
           GOTO 4229

 4224      CONTINUE
           GOTO 4229

 4226      CONTINUE
           GOTO 4229
 4229   CONTINUE

      ENDIF

      GO TO 90000


C--------- MINTEST FUNCTION   GULF ( N, X, F, G, IFG, NINT(FARG(1)))
C--------- NINT(FARG(1)) IS M


 4300 X1 = X(1)
      X2 = X(2)
      X3 = X(3)

      IF ( X1 .EQ. ZERO ) THEN
         RET = NOFG
         GOTO 90000
      ENDIF

      IF ( .NOT. GONLY ) THEN
         F = ZERO
      ENDIF

      IF ( .NOT. FONLY ) THEN
         G1 = ZERO
         G2 = ZERO
         G3 = ZERO
      ENDIF

      DO 4310  I = 1,NINT(FARG(1))
         W1 = RD(I)/R100
         YI = (-R50*LOG(W1)) ** R2D3 + R25 - X2
         IF ( X3*LOG(ABS(YI)) .LE. BIGGST ) THEN
            W2 = ( ABS(YI)**X3 ) /  X1
         ELSE
            RET = NOFG
            GOTO 90000
         ENDIF
         IF ( -W2 .LE. BIGGST ) THEN
            W3 = EXP(-W2)
         ELSE
            RET = NOFG
            GOTO 90000
         ENDIF
         RI = W3 - W1

         IF ( .NOT. GONLY ) THEN
            F = F + RI*RI
         ENDIF

         IF ( .NOT. FONLY ) THEN
            W1 = RI*W2*W3
            G1 = G1 + W1
            G2 = G2 + W1 / YI
            G3 = G3 - W1 * LOG(ABS(YI))
         ENDIF

 4310 CONTINUE

      IF ( .NOT. FONLY ) THEN
         G(1) = TWO * G1 / X1
         G(2) = TWO * G2 * X3
         G(3) = TWO * G3
      ENDIF

      GOTO 90000


C--------- MINTEST FUNCTION     ARGAUS


 4400 IF ( .NOT. GONLY ) THEN
         F = ZERO
      ENDIF

      X1 = X(1)
      X2 = X(2)
      X3 = X(3)

      G1 = ZERO
      G2 = ZERO
      G3 = ZERO

      DO 4410 I = 1, 15
         W1 = FOUR - RD(I)/TWO - X3
         W2 = -HALF * X2 * W1 * W1
         IF ( W2 .LE. BIGGST ) THEN
            W3 =  EXP(W2)
         ELSE
            RET = NOFG
            GOTO 90000
         ENDIF
         RI = X1 * W3 - ARGASY(I)

         IF ( .NOT. GONLY ) THEN
            F = F + RI * RI
         ENDIF

         IF ( .NOT. FONLY ) THEN
            W2 = W3 * RI
            W3 = W1 * W2
            G1 = G1 + W2
            G2 = G2 - W1 * W3
            G3 = G3 + W3
         ENDIF

 4410 CONTINUE

      IF ( .NOT. FONLY ) THEN
         G(1) = TWO * G1
         G(2) =  X1 * G2
         G(3) = TWO * X1 * X2 * G3
      ENDIF

      GOTO 90000


C--------- MINTEST FUNCTION     MEYER

 4500 X1 = X(1)
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

      DO 4510 I = 1,16
         TI = R45 + RD(I)*FIVE
         W1 = TI + X3
         IF ( X2/W1 .LE. BIGGST ) THEN
            W2 = EXP(X2/W1)
         ELSE
            RET = NOFG
            GOTO 90000
         ENDIF
         RI = X1*W2 - MEY(I)

         IF ( .NOT. GONLY ) THEN
            F = F + RI*RI
         ENDIF

         IF ( .NOT. FONLY ) THEN
            W2 = RI*W2
            G1 = G1 + W2
            G2 = G2 + W2/W1
            G3 = G3 + W2/(W1*W1)
         ENDIF

 4510 CONTINUE

      IF ( .NOT. FONLY ) THEN
         G(1) =    TWO*G1
         G(2) =    TWO*G2*X1
         G(3) =  - TWO*G3*X1*X2
      ENDIF

      GOTO 90000

C--------- MINTEST FUNCTION   BROWND ( N, X, F, G, IFG, NINT(FARG(1)))
C--------- NINT(FARG(1)) IS M


 4600 IF ( .NOT. GONLY ) THEN
         F = ZERO
      ENDIF

      IF ( .NOT. FONLY ) THEN
         G1 = ZERO
         G2 = ZERO
         G3 = ZERO
         G4 = ZERO
      ENDIF

      DO  4610 I = 1,NINT(FARG(1))
         W1 = RD(I)/FIVE
         W2 = SIN(W1)
         IF ( W1 .LE. BIGGST ) THEN
            W3 = X(1) + W1*X(2) - EXP(W1)
         ELSE
            RET = NOFG
            GOTO 90000
         ENDIF
         W4 = X(3) + W2*X(4) - COS(W1)
         RI = W3*W3 + W4*W4

         IF ( .NOT. GONLY ) THEN
            F = F + RI*RI
         ENDIF

         IF ( .NOT. FONLY ) THEN
            W3 = W3*RI
            W4 = W4*RI
            G1 = G1 + W3
            G2 = G2 + W1*W3
            G3 = G3 + W4
            G4 = G4 + W2*W4
         ENDIF

 4610 CONTINUE

      IF ( .NOT. FONLY ) THEN
         G(1) = FOUR * G1
         G(2) = FOUR * G2
         G(3) = FOUR * G3
         G(4) = FOUR * G4
      ENDIF

      GOTO 90000


C--------- MINTEST FUNCTION     KOWOSB


 4700 IF ( .NOT. GONLY ) THEN
         F = ZERO
      ENDIF

      IF ( .NOT. FONLY ) THEN
         G1 = ZERO
         G2 = ZERO
         G3 = ZERO
         G4 = ZERO
      ENDIF

      X1 = X(1)
      X2 = X(2)
      X3 = X(3)
      X4 = X(4)

      DO 4710 I = 1,11
         W1 = KOWOSU(I)
         W2 = W1*W1
         W3 = W2 + W1*X2
         W4 = W2 + W1*X3 + X4
         W5 = W3/W4
         RI = KOWOSY(I) - X1*W5

         IF ( .NOT. GONLY ) THEN
            F = F + RI*RI
         ENDIF

         IF ( .NOT. FONLY ) THEN
            G1 = G1 + RI*W5
            G2 = G2 + RI*W1/W4
            W5 = RI*W5/W4
            G3 = G3 + W5*W1
            G4 = G4 + W5
         ENDIF
 4710 CONTINUE

         IF ( .NOT. FONLY ) THEN
            G(1) = - TWO*G1
            G(2) = - TWO*G2*X1
            G(3) =   TWO*G3*X1
            G(4) =   TWO*G4*X1
         ENDIF

      GOTO 90000


C--------- MINTEST FUNCTION     OSBRN1

 4800 IF ( .NOT. GONLY ) THEN
         F = ZERO
      ENDIF

      IF ( .NOT. FONLY ) THEN
         G1 = ZERO
         G2 = ZERO
         G3 = ZERO
         G4 = ZERO
         G5 = ZERO
      ENDIF

      X2 = X(2)
      X3 = X(3)

      DO 4810 I = 1,33
         TI = TEN*(REAL(I) - ONE)
         IF ( MAX(-TI*X(4),-TI*X(5)) .LE. BIGGST ) THEN
            W1 = EXP(-TI*X(4))
            W2 = EXP(-TI*X(5))
         ELSE
            RET = NOFG
            GOTO 90000
         ENDIF
         RI = X(1) + X2*W1 + X3*W2 - OSB1Y(I)

         IF ( .NOT. GONLY ) THEN
            F = F + RI*RI
         ENDIF

         IF ( .NOT. FONLY ) THEN
            W1 =  W1*RI
            W2 =  W2*RI
            G1 = G1 + RI
            G2 = G2 + W1
            G3 = G3 + W2
            G4 = G4 + W1*TI
            G5 = G5 + W2*TI
         ENDIF

 4810 CONTINUE

      IF ( .NOT. FONLY ) THEN
         G(1) =  TWO *  G1
         G(2) =  TWO *  G2
         G(3) =  TWO *  G3
         G(4) = -TWO*X2*G4
         G(5) = -TWO*X3*G5
      ENDIF

      GOTO 90000

C--------- MINTEST FUNCTION     OSBRN2

 4900 IF ( .NOT. GONLY ) THEN
         F = ZERO
      ENDIF

      IF ( .NOT. FONLY ) THEN
         G1 = ZERO
         G2 = ZERO
         G3 = ZERO
         G4 = ZERO
         G5 = ZERO
         G6 = ZERO
         G7 = ZERO
         G8 = ZERO
         G9 = ZERO
         G10= ZERO
         G11= ZERO
      ENDIF

      X1 = X(1)
      X2 = X(2)
      X3 = X(3)
      X4 = X(4)
      X5 = X(5)
      X6 = X(6)
      X7 = X(7)
      X8 = X(8)

      DO 4910 I = 1,65

         TI = (REAL(I) - ONE)/TEN
         X9 = TI - X(9)
         W6 = X9*X9
         X10 = TI - X(10)
         W7 = X10*X10
         X11 = TI - X(11)
         W8 = X11*X11
         IF ( MAX(-TI*X5,-X6*W6,-X7*W7,-X8*W8) .LE. BIGGST ) THEN
            W1 = EXP(-TI*X5)
            W2 = EXP(-X6*W6)
            W3 = EXP(-X7*W7)
            W4 = EXP(-X8*W8)
         ELSE
            RET = NOFG
            GOTO 90000
         ENDIF

         RI = X1*W1 + X2*W2 + X3*W3 + X4*W4 - OSB2Y(I)

         IF ( .NOT. GONLY ) THEN
            F = F + RI*RI
         ENDIF

         IF ( .NOT. FONLY ) THEN
            W1 = RI*W1
            W2 = RI*W2
            W3 = RI*W3
            W4 = RI*W4
            G1 = G1 + W1
            G2 = G2 + W2
            G3 = G3 + W3
            G4 = G4 + W4
            G5 = G5 + W1 * TI
            G6 = G6 + W2 * W6
            G7 = G7 + W3 * W7
            G8 = G8 + W4 * W8
            G9 = G9 + W2 * X9
            G10 = G10 + W3 * X10
            G11 = G11 + W4 * X11
         ENDIF

 4910 CONTINUE

      IF ( .NOT. FONLY ) THEN
         G(1) =   TWO * G1
         G(2) =   TWO * G2
         G(3) =   TWO * G3
         G(4) =   TWO * G4
         G(5) = - TWO * X1 * G5
         G(6) = - TWO * X2 * G6
         G(7) = - TWO * X3 * G7
         G(8) = - TWO * X4 * G8
         G(9) =  FOUR * X2 * X6 * G9
         G(10) = FOUR * X3 * X7 * G10
         G(11) = FOUR * X4 * X8 * G11
      ENDIF

      GOTO 90000


C---------                              WATSON
 5200 ERROR = N .LT. 2  .OR.  N .GT. 31
      IF ( .NOT. GONLY  .OR.  ERROR ) THEN
         F = ZERO
      ENDIF
      IF ( .NOT. FONLY  .OR.  ERROR ) THEN
         DO 5210 I = 1, N
            G(I) = ZERO
 5210    CONTINUE
      ENDIF

      IF ( .NOT. ERROR ) THEN
         DO 5280 I = 1, 29
            W1 = RD29 * RD(I)
            W2 = ZERO
            W3 = ONE
            DO 5230 J = 2, N
               W2 = W2 + W3 * X(J) * RD( J - 1 )
               W3 = W3 * W1
 5230       CONTINUE

            W4 = ZERO
            W3 = ONE
            DO 5240 J = 1, N
               W4 = W4 + W3 * X(J)
               W3 = W3 * W1
 5240       CONTINUE

            W5 = W2 - W4 * W4 - ONE
            IF ( .NOT. GONLY ) THEN
               F = F + W5 * W5
            ENDIF

            IF ( .NOT. FONLY ) THEN
               W6 = TWO * W1 * W4
               W3 = TWO / W1
               DO 5250 J = 1, N
                  G(J) = G(J) + W3 * W5 * (REAL( J - 1) - W6 )
                  W3   = W3 * W1
 5250          CONTINUE
            ENDIF
 5280    CONTINUE

         W1 = X(1)
         W2 = X(2)
         W3 = W2 - W1 * W1 - ONE
         IF ( .NOT. GONLY ) THEN
            F = F + W1 * W1 + W3 * W3
         ENDIF
         IF ( .NOT. FONLY ) THEN
            G(1) = G(1) + W1 * ( TWO - FOUR * W3 )
            G(2) = G(2) + TWO * W3
         ENDIF

      ENDIF
      GO TO 90000

C--------- MINTEST FUNCTION     PENAL3 ( N, X, F, G, IFG,
C                                         FARG(1), WORK, SIZE )
C--------- FARG(1) IS A

 5300 RF1 = FARG ( 1 )

      IF ( SIZE .LT. 2 * N ) THEN
         F = ZERO
         DO 5310 K = 1, N
            G(K) = ZERO
 5310    CONTINUE
         GO TO 90000
      ENDIF

      I1 = 0
      I2 = N

      W1     = RD(N)
      NOVER2 = N / 2

      W2 = X(N)
      W3 = X(N-1)
      IF ( MAX(W2,W3) .LE. BIGGST ) THEN
         W4  = EXP(W2)
         W5  = EXP(W3)
      ELSE
         RET = NOFG
         GOTO 90000
      ENDIF

      R   = ZERO
      S   = ZERO
      T   = ( W2*W2 ) + ( W3*W3 ) - ( TWO*W1 )

      IF ( .NOT. GONLY ) THEN
         U = ZERO
      ENDIF

      DO 5320  K = (N - 2), 1, -1
         XK = X(K)
         RK =     XK + TWO*W3 + TEN*W2 - ONE
         SK = TWO*XK +     W3          - THREE
         R  = R + RK*RK
         S  = S + SK*SK
         T  = T + XK*XK - W1

         IF ( .NOT. FONLY ) THEN
            WORK(I1+K) = RK
            WORK(I2+K) = SK
         ENDIF

         IF ( .NOT. GONLY .AND. K .LE. NOVER2 ) THEN
            RK = XK - ONE
            U  = U + RK*RK
         ENDIF

         W2 = W3
         W3 = XK

 5320 CONTINUE

      IF ( .NOT. GONLY ) THEN
         F = RF1 * ( ONE + W4*R + W5*S + R*S ) + T*T + U
      ENDIF

      IF ( .NOT. FONLY ) THEN
         W2 = TWO * ( W4 + S )
         W1 = TWO * ( W5 + R )
         W3 = FOUR*T

         W6 = WORK(I1+1)
         W7 = WORK(I1+2)
         W8 = WORK(I2+2)

         XK   = X(1)
         G(1) =  RF1 * ( W2 * ( W6                       )
     -        +        W1 * ( TWO*WORK(I2+1)           ) )
     -        +  W3*XK + TWO*(XK - ONE)

         XK   = X(2)
         G(2) =  RF1 * ( W2 * (     W7 + TWO*W6          )
     -        +        W1 * ( TWO*W8 + WORK(I2+1)      ) )
     -        +  W3*XK + TWO*(XK - ONE)

         DO 5330 K = 3, (N - 2)

            RK = WORK(I1+K)
            SK = WORK(I2+K)
            XK = X(K)

            W9 = RF1 * ( W2 * (     RK + TWO*W7 + TEN*W6 )
     -         +       W1 * ( TWO*SK +     W8          ) )
     -         + W3*XK

            IF ( K .LE. NOVER2 ) THEN
               W9 = W9 + TWO * ( XK - ONE )
            ENDIF

            G(K) = W9

            W6 = W7
            W7 = RK
            W8 = SK

 5330    CONTINUE

C     ---NOTE THAT BELOW W7 = RK(N-2), W6 = RK(N-3)
C     ---AND W8 = SK(N-2).

         G(N-1) = RF1 * ( W2 * (          TWO*W7 + TEN*W6 )
     -          +       W1 * (              W8          )
     -          +       S  * W5                           )
     -          + W3*X(N-1)

         G(N)   = RF1 * ( W2 * (                   TEN*W7 )
     -          +       R  * W4                           )
     -          + W3*X(N)

      ENDIF

      GOTO 90000


C--------- MINTEST FUNCTION     VAROSB ( N, X, F, G, IFG, FARG(1) )
C--------- FARG(1) IS LAMBDA


 5400 W2 = ONE / ( RD(N) + ONE )
      W5 = TWO / W2
      W2 = TWO * FARG(1) * W2

      W6 = ZERO
      W9 =  ONE

      W1 = ZERO
      TI = ZERO

      DO 5410 I = 1,N

         W7 = X(I)
         IF ( W7 .LE. BIGGST ) THEN
            W8 = EXP(W7)
         ELSE
            RET = NOFG
            GOTO 90000
         ENDIF
         W3 = W7 - W6

         IF ( .NOT. GONLY ) THEN

            W1 = W1 + W7*W7 - W6*W7

            IF ( W3 .NE. ZERO ) THEN
               TI = TI + (W8 - W9)/ W3
            ELSE
               TI = TI +     W8
            ENDIF

         ENDIF

         IF ( .NOT. FONLY ) THEN

            IF ( W3 .NE. ZERO ) THEN
               W4 = ( W9 + W8 * ( W3 - ONE ) ) / ( W3*W3 )
            ELSE
               W4 = W9 / TWO
            ENDIF

            G(I) = W5 * W3   +  W2 * W4

            IF ( I .GT. 1 ) THEN

               IF ( W3 .NE. ZERO ) THEN
                  W4 = (W8 - W9*(ONE + W3))/(W3*W3)
               ELSE
                  W4 =  W8 / TWO
               ENDIF

               G(I-1) = G(I-1) - W5 * W3  +  W2 * W4

            ENDIF

         ENDIF

         W6  = W7
         W9  = W8

 5410 CONTINUE

      IF ( .NOT. FONLY ) THEN

         IF ( W7 .NE. ZERO ) THEN
            W4 = ( ONE + ( W7 - ONE ) * W8 ) / ( W7*W7)
         ELSE
            W4 = W8 / TWO
         ENDIF

         G(N) = G(N) + W5 * W7  +  W2 * W4

      ENDIF

      IF ( .NOT. GONLY ) THEN

         IF ( W6 .NE. ZERO ) THEN
            W3 = ( W9 - ONE ) / W6
         ELSE
            W3 =   W9
         ENDIF

         F = W5 * W1 + W2 * ( TI + W3 )

      ENDIF

      GOTO 90000


C---------      VARDIM
 5500 W1 = ZERO
      W2 = ZERO

      DO 5510 J = 1, N
         W3 = X(J) - ONE
         W1 = W1 + RD(J)*W3
         W2 = W2 + W3*W3
 5510 CONTINUE

      W4 = W1*W1

      IF ( .NOT. GONLY ) THEN
         F = W2 + W4 * (ONE + W4)
      ENDIF

      IF ( .NOT. FONLY ) THEN

         W3 = W1 * (ONE + TWO * W4)

         DO 5520 J = 1, N
            G(J) = TWO * (X(J) - ONE + RD(J) * W3)
 5520    CONTINUE

      ENDIF
      GOTO 90000

C---------                  RECIPE
 5600               CONTINUE
      W2 = (X(2) - X(1))**2
      W3 = (X(2) - X(1))**3

      IF ( .NOT. GONLY ) THEN
         F = (X(1)-FIVE)**2 + X(2)**2 + X(3)**2/W2
      ENDIF

      IF ( .NOT. FONLY ) THEN
         G(1) = TWO*(X(1) - FIVE) + TWO*X(3)**2/W3
         G(2) = TWO*X(2) - TWO*X(3)**2/W3
         G(3) = TWO*X(3)/W2
      ENDIF
      GOTO 90000

C---------                  CLIFF
 5700             CONTINUE
      W1 = (X(1) - THREE)/R100
      W2 = EXP(R20*(X(1)-X(2)))
      IF ( .NOT. GONLY ) THEN
         F = W1**2 - X(1) + X(2) + W2
      ENDIF
      IF ( .NOT. FONLY ) THEN
         G(1) = W1/R50 - ONE + R20*W2
         G(2) = ONE - R20*W2
      ENDIF
      GOTO 90000

C---------                  CHEBYQ ( N, X, F, G, IFG, WORK, SIZE )
 5800 IF ( SIZE .LT. N ) THEN
         F = ZERO
         DO 5810 I = 1, N
            G(I) = ZERO
 5810    CONTINUE

      ELSE
         DO 5820 I = 1, N
            WORK(I) = ZERO
 5820    CONTINUE

         DO 5830 J = 1, N
            W1 = ONE
            W2 = TWO * X(J) - ONE
            W3 = TWO * W2
            DO 5840 I = 1, N
               WORK(I) = WORK(I) + W2
               W4      = W3 * W2 - W1
               W1      = W2
               W2      = W4
 5840       CONTINUE
 5830    CONTINUE

         IF ( .NOT. GONLY ) THEN
            F = ZERO
         ENDIF

         W7 = ONE / RD(N)
         EVEN = .FALSE.
         DO 5850 I = 1, N
            W2 = W7 * WORK(I)
            IF ( EVEN ) THEN
               W2 = W2 + ONE / ( RD(I)**2 - ONE )
            ENDIF
            EVEN = .NOT. EVEN
            IF ( .NOT. GONLY ) THEN
               F = F + W2 * W2
            ENDIF
            IF ( .NOT. FONLY ) THEN
               WORK(I) = W2
            ENDIF
 5850    CONTINUE

         IF ( .NOT. FONLY ) THEN
            DO 5860 J = 1, N
               G(J) = ZERO
               W1   = ONE
               W2   = TWO * X(J) - ONE
               W3   = TWO * W2
               W4   = ZERO
               W5   = TWO

               DO 5870 I = 1, N
                  G(J) = G(J) + WORK(I) * W5
                  W6   = FOUR * W2 + W3 * W5 - W4
                  W4   = W5
                  W5   = W6
                  W6   = W3 * W2 - W1
                  W1   = W2
                  W2   = W6
 5870          CONTINUE
 5860       CONTINUE

            W7 = W7 * TWO
            CALL ZZSCAL ( N, W7, G, 1 )

         ENDIF

      ENDIF
      GOTO 90000


C--------- MINTEST FUNCTION     HIMM32


 5900 IF ( .NOT. GONLY ) THEN
         F = ZERO
      ENDIF
      IF ( .NOT. FONLY ) THEN
         G1 = ZERO
         G2 = ZERO
         G3 = ZERO
         G4 = ZERO
      ENDIF
      X1 = X(1)
      X2 = X(2)
      X3 = X(3)
      X4 = X(4)
      W1 = X1*X1
      W2 = X2*X2
      W3 = X3*X3
      W4 = X4*X4

      DO 5910 I = 1,7
         W8 = HIM32A(I)
         W9 = HIM32B(I)
         W5 = W1 + W2*W8 + W3*W8*W8
         W6 = W9*(ONE + W4*W8)
         RI = W5/W6 - ONE
         IF ( .NOT. GONLY ) THEN
            F = F + RI*RI
         ENDIF
         IF ( .NOT. FONLY ) THEN
            W7 = RI/W6
            G1 = G1 + W7
            W7 = W7*W8
            G2 = G2 + W7
            G3 = G3 + W7*W8
            G4 = G4 + W7*W5*W9/W6
         ENDIF
 5910 CONTINUE
      IF ( .NOT. GONLY ) THEN
         F = F*R10000
      ENDIF
      IF ( .NOT. FONLY ) THEN
         G(1) =  R40000*X1*G1
         G(2) =  R40000*X2*G2
         G(3) =  R40000*X3*G3
         G(4) = -R40000*X4*G4
      ENDIF

      GOTO 90000


C--------- MINTEST FUNCTION     HIMM27


 6000 X1 = X(1)
      X2 = X(2)
      W1 = X1*X2
      W2 = ONE - X1
      W3 = W2*W2*W2*W2
      W4 = ONE - X2 - X1*W2*W3
      W5 = W1*W2*W4

      IF ( .NOT. GONLY ) THEN
         F = W5*W5
      ENDIF

      IF ( .NOT. FONLY ) THEN
         W5   = TWO*W5
         G(1) = W5 * (  X2*W2*W4 - W1*W4 + W1*W2*W3*( SIX*X1 - ONE )  )
         G(2) = W5*W2 * ( X1*W4 - W1 )
      ENDIF

      GOTO 90000


C---------                       BRYTRI
 6100 IF ( .NOT. GONLY ) THEN
         F = ZERO
      ENDIF

      KAP1 = FARG(1)
      KAP2 = FARG(2)
      W1 = X(1)
      W2 = ZERO
      W3 = ZERO

      DO 6110 I = 1,N

         IF ( I .LT. N ) THEN
            W4 = X(I+1)
         ELSE
            W4 = ZERO
         ENDIF
         RI = (THREE - KAP1*W1)*W1 - W2 - TWO*W4 + KAP2

         IF ( .NOT. GONLY ) THEN
            F = F + RI*RI
         ENDIF

         IF ( .NOT. FONLY ) THEN
            IF ( I .GT. 1 ) THEN
               G(I-1) = TWO*(W5 - RI)
            ENDIF
            W5 = -TWO*W3 + RI*(THREE - TWO*KAP1*W1)
            W3 = RI
         ENDIF

         W2 = W1
         W1 = W4
 6110 CONTINUE
      G(N) = TWO*W5
      GOTO 90000

C---------           BRWNAL
 6200 S  = ZERO
      T  = ONE
      W1 = ZERO

      DO 6210 I = 1,N
         XI = X(I)
         S  = S + XI
         T  = T * XI
 6210 CONTINUE

      W2 = T - ONE
      IF ( .NOT. GONLY ) THEN
         F = W2 * W2
      ENDIF

      W2 = W2*T
      W3 = RD(N+1)
      DO 6220 I = 1,N-1
         XI = X(I)
         RI = XI + S - W3
         IF ( .NOT. GONLY ) THEN
            F = F + RI*RI
         ENDIF
         IF ( .NOT. FONLY ) THEN
            G(I) =  RI + W2/XI
            W1 = W1 + RI
         ENDIF
 6220 CONTINUE

      DO 6230 I = 1,N-1
         G(I) = TWO * (W1 + G(I))
 6230 CONTINUE
      G(N) = TWO * (W1 + W2/X(N))

      GOTO 90000


C## E X I T
90000       STATUS(7) =  RET
91000       GFIRST    = .FALSE.
      RETURN

C## F O R M A T S:  NONE ARE DEFINED.
C##                 E N D         OF ZZFNS3.
                    END
