      SUBROUTINE ZZVECT ( ACTION, N, U, V, ARGS )

C## A R G U M E N T S:
                     INTEGER           ACTION,    N
                     REAL              U(N), V(N), ARGS(*)
C!!!!                DOUBLE PRECISION  U(N), V(N), ARGS(*)
C## S T A T U S:
C               IGNORE LINES BEGINNING WITH  "C!!!!" .
C
C               SYSTEM  DEPENDENCE:                      NONE.
C               SINGLE/DOUBLE CONVERSION: NEEDED (SEE CONVRT).
C
C               THIS VERSION IS IN   S I N G L E   PRECISION.
C!!!!           THIS VERSION IS IN   D O U B L E   PRECISION.

C>RCS $HEADER: VECT.F,V 2.2 91/11/20 10:53:12 BUCKLEY EXP $
C>RCS $LOG:     VECT.F,V $
C>RCS REVISION 2.2  91/11/20  10:53:12  BUCKLEY
C>RCS FINAL SUBMISSION TO TOMS
C>RCS
C>RCS REVISION 2.1  90/10/04  14:44:32  BUCKLEY
C>RCS FIXED MISSING GOTO 10000 FOR CASE 9.
C>RCS
C>RCS REVISION 2.0  90/07/31  11:38:24  BUCKLEY
C>RCS MINOR MOD TO USE BLAS
C>RCS
C>RCS REVISION 1.9  89/06/30  13:39:55  BUCKLEY
C>RCS PREPARING SUBMITTED VERSION OF MT
C>RCS
C>RCS REVISION 1.3.1.1  89/05/20  16:43:03  BUCKLEY
C>RCS TEMP. TEST OF MT BEFORE SUBMITTING
C>RCS
C>RCS REVISION 1.1  89/01/17  16:50:14  BUCKLEY
C>RCS INITIAL REVISION
C>RCS

C## D E S C R I P T I O N:

C    COMPUTE STANDARD VECTORS AND PUT INTO V. THE VALUE PUT INTO
C    V DEPENDS ON THE VALUE OF ACTION:
C
C      IF ACTION = 0, IT IS EQUIVALENT TO SPECIFYING ACTION = N.
C      IF ACTION > 0, A CYCLE IS DEFINED AND MUST BE EXPANDED. IN
C                       THIS CASE, THE FIRST "ACTION" COMPONENTS
C                       OF U MUST DEFINED. THEY ARE CYCLICALLY
C                       COPIED INTO V. IN THIS CASE, U AND V MAY
C                       REFER TO THE SAME VECTORS IN THE CALL.
C      IF ACTION < 0  A SELECTED FORMULA IS USED TO COMPUTE VALUES
C                       FOR THE COMPONENTS OF V.
C
C## E N T R Y   P O I N T S: THE NATURAL ENTRY ZZVECT.

C## S U B R O U T I N E S:
C                         REAL(DBLE),ABS,MOD  INTRINSIC
C                         RD                  STATEMENT FUNCTION

C## P A R A M E T E R S:

      REAL              TENTH,         FIFTH,        HALF
C!!!! DOUBLE PRECISION  TENTH,         FIFTH,        HALF
      PARAMETER (       TENTH = .1D0,  FIFTH = .2D0, HALF = .5D0      )

      REAL              RPT9,        RPT8,        RD29
C!!!! DOUBLE PRECISION  RPT9,        RPT8,        RD29
      PARAMETER (       RPT9 = .9D0, RPT8 = .8D0, RD29 = 1D0/29D0 )
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
      REAL              ZERO,       ONE,       TWO,       THREE
C!!!! DOUBLE PRECISION  ZERO,       ONE,       TWO,       THREE
      PARAMETER (       ZERO = 0D0, ONE = 1D0, TWO = 2D0, THREE = 3D0)

      REAL              FOUR,       FIVE,      SIX,       SEVEN
C!!!! DOUBLE PRECISION  FOUR,       FIVE,      SIX,       SEVEN
      PARAMETER (       FOUR = 4D0, FIVE = 5D0,SIX = 6D0, SEVEN = 7D0)

      REAL              EIGHT,         NINE,          TEN
C!!!! DOUBLE PRECISION  EIGHT,         NINE,          TEN
      PARAMETER (       EIGHT = 8D0,   NINE = 9D0,    TEN = 10D0     )

      INTEGER     JUSTF,     BOTH,     JUSTG,      NOOP
      PARAMETER ( JUSTF = 1, BOTH = 0, JUSTG = -1, NOOP = 2 )

C## L O C A L   D E C L:
                         INTEGER           I,  K
                         REAL              R,  Z ,RD
C!!!!                    DOUBLE PRECISION  R,  Z ,RD

C## S A V E:                 NONE SELECTED.
C## E Q U I V A L E N C E S: NONE ARE DEFINED.
C## C O M M O N:             NONE IS DEFINED.
C## D A T A:                 NO DATA VALUES SET.

C##                                             E X E C U T I O N
C##                                             E X E C U T I O N
C------               FUNCTION STATEMENT.
      RD(I) = REAL(I)
C!!!! RD(I) = DBLE(I)
C-----                NOW EXECUTE.

      IF ( ACTION .LT. 0 ) THEN

        GOTO(
     -    100,  200,  300,  400,  500,  600,  700,  700,  900,  900,
     -   1100, 1200, 1300, 1400, 1500, 1600, 1700, 1800, 1900, 2000 )
     -         , ABS( ACTION )

C        ---VECTOR 1 : V(I) = I/(N+1)
  100       R = ONE / RD(N+1)
            DO 150 I = 1, N
               V(I) = R * RD(I)
  150       CONTINUE
         GOTO 10000

C        ---VECTOR 2 : V(I) = .1*I*(1-I)/(N+1)
  200       R = TENTH / RD(N+1)
            DO 250 I = 1, N
               V(I) = R * RD(I) * RD(1-I)
  250       CONTINUE
         GOTO 10000

C        ---VECTOR 3 : V = -1.2, 1, 1, 1, ... , 1
  300       V(1) = - R1P2
            DO 350 I = 2, N
               V(I) = ONE
  350       CONTINUE
         GOTO 10000

C        ---VECTOR 4 : V = 1, 2, 3, ..., N
  400       DO 450 I = 1, N
               V(I) = RD(I)
  450       CONTINUE
         GOTO 10000

C        ---VECTOR 5 : V(I) = 1 - I/N
  500       R = ONE / RD(N)
            DO 550 I = 1, N
               V(I) = ONE - R * RD(I)
  550       CONTINUE
         GOTO 10000

C        ---VECTOR 6 : V(I) = Z*(Z - 1) WHERE Z = I/(N+1)
  600       R = ONE / RD(N+1)
            DO 650 I = 1, N
               Z     = R * RD(I)
               V(I) = Z * ( Z - ONE )
  650       CONTINUE
         GOTO 10000

C        ---VECTORS 7,8 : V(I) = 0.2 * SIN((REAL(I**2))
C                       IN VECTOR 8, WHEN I = 2M+1, THE CONSTANT 0.2 IS
C                       REPLACED BY -0.8. NOTE THAT N = M^2.

  700       DO 730 I = 1,N
               V(I) = SIN(RD(I**2)) * FIFTH
  730       CONTINUE

            IF ( ABS(ACTION) .EQ. 8 ) THEN
               I = 2 * NINT(SQRT(RD(N))) + 1
               V(I) = - RPT8 * SIN(RD(I**2))
            ENDIF
         GOTO 10000

C        ---VECTORS 9,10 : V(I) DOCUMENTED IN TOI83B.
C              FOR LINEAR (9) AND NONLINEAR (10) MINIMUM SURFACE PROBLEMS,
C              FUNCTIONS MNSRF1 AND MNSRF2.
  900       CONTINUE
            DO 950 I=1,N
               V(I)= ZERO
  950       CONTINUE

            K = NINT( SQRT(RD(N)) )

            IF ( ABS(ACTION) .EQ. 9 ) THEN
               DO 952 I=1,K
                  Z  = (RD(I)-ONE)/(RD(K)-ONE)
                  V(I)           =  ONE + FOUR*Z
                  V((I-1)*K + 1) =  ONE + EIGHT*Z
                  V(I + K*(K-1)) = NINE + FOUR*Z
                  V(I*K)         = FIVE + EIGHT*Z
  952          CONTINUE
            ELSE IF ( ABS(ACTION) .EQ. 10 ) THEN
               DO 954 I=1,K
                  Z  = (RD(I)-ONE)/(RD(K)-ONE)
                  V(I)           =  ONE + FOUR*Z  + TEN*(ONE+Z)**2
                  V((I-1)*K + 1) =  ONE + EIGHT*Z + TEN*(ONE-Z)**2
                  V(I + K*(K-1)) = NINE + FOUR*Z  + TEN*Z**2
                  V(I*K)         = FIVE + EIGHT*Z + TEN*(TWO-Z)**2
  954          CONTINUE
            ENDIF
         GOTO 10000

C        ---VECTOR 11 : V = 1/N.  USED IN FUNCTION ARTRIG.
 1100    CONTINUE
           DO 1150 I = 1, N
               V(I) = ONE/RD(N)
 1150      CONTINUE
         GOTO 10000

C        ---VECTOR 12.  USED IN FUNCTION MANCIN.
 1200    CONTINUE
           DO 1250 I = 1, N
               U(I) = ZERO
 1250      CONTINUE
           CALL ZZFNS ( JUSTF, N, U, Z, U, K, U, V )
           Z = -RD(14*N) / ( 196*N**2 - (36*(N-1)**2) )
           CALL ZZSCAL ( N, Z, V, 1 )
         GOTO 10000

C        ---VECTOR 13 : V = 4/I.  USED IN FUNCTION HILBRT.
 1300    CONTINUE
           DO 1350 I = 1, N
               V(I) = FOUR/RD(I)
 1350      CONTINUE
         GOTO 10000

C        ---VECTOR 14 : V = -1, 1, 1, 1, ... , 1, 1, 1
C                     USED IN FUNCTION GENRSN.
 1400       V(1) = - ONE
            DO 1450 I = 2, N
               V(I) = ONE
 1450       CONTINUE
         GOTO 10000

C        ---VECTOR 15 : V = (0, 0,...,0, N+1).   USED IN FUNCTION BRWNAL.
 1500      DO 1550 I = 1, N-1
               V(I) = ZERO
 1550      CONTINUE
           V(N) = RD(N+1)
         GOTO 10000

C        ---VECTOR 16 : V = (T,T,...,T,T^{1-N}).  USED IN FUNCTION BRWNAL
C        THE SOLUTION COMPONENT T IS PASSED INTO ZZVECT FROM ZZDSOL VIA
C        VARIABLE ARGS.  IF THE VECTOR NUMBER IS CHANGED, IT MAY BE
C        NECESSARY TO CHANGE THE VALUE OF ARGS IN ZZDSOL.
 1600    CONTINUE
         DO 1650 I = 1, N-1
            V(I) = ARGS(1)
 1650    CONTINUE
         V(N) = ARGS(1)**(1-N)
         GOTO 10000

C        ---VECTOR 17 : V = (T,T,...,T).  UNUSED AT PRESENT.
C        THE SOLUTION COMPONENT T IS PASSED INTO ZZVECT FROM ZZDSOL VIA
C        VARIABLE ARGS.  IF THE VECTOR NUMBER IS CHANGED, IT MAY BE
C        NECESSARY TO CHANGE THE VALUE OF ARGS IN ZZDSOL.
 1700    CONTINUE
         DO 1750 I = 1, N
            V(I) = ARGS(1)
 1750    CONTINUE
         GOTO 10000

C        ---VECTOR 18 : V = (T,T,0).  USED IN FUNCTION BOX66.
C        THE SOLUTION COMPONENT T IS PASSED INTO ZZVECT FROM ZZDSOL VIA
C        VARIABLE ARGS.  IF THE VECTOR NUMBER IS CHANGED, IT MAY BE
C        NECESSARY TO CHANGE THE VALUE OF ARGS IN ZZDSOL.
 1800    CONTINUE
         V(1) = ARGS(1)
         V(2) = ARGS(1)
         V(3) = ZERO
         GOTO 10000

C        ---VECTOR 19 : V = (0,0,T,T).  USED IN FUNCTION TOIN4.
C        THE SOLUTION COMPONENT T IS PASSED INTO ZZVECT FROM ZZDSOL VIA
C        VARIABLE ARGS.  IF THE VECTOR NUMBER IS CHANGED, IT MAY BE
C        NECESSARY TO CHANGE THE VALUE OF ARGS IN ZZDSOL.
 1900    CONTINUE
         V(1) = ZERO
         V(2) = ZERO
         V(3) = ARGS(1)
         V(4) = ARGS(1)
         GOTO 10000

 2000    CONTINUE
         GOTO 10000

      ELSE IF ( ACTION .GT. 0 ) THEN
         K = 0
         DO 5000 I = 1, N
            K = MOD ( K, ACTION ) + 1
            V(I) = V(K)
5000     CONTINUE
      ENDIF
10000 GOTO 90000

C## E X I T
90000      RETURN

C## F O R M A T S:  NONE ARE DEFINED.
C##                 E N D         OF ZZVECT.
                    END
