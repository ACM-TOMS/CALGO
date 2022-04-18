      SUBROUTINE ZZSOLN ( IX, N, SOLN, X )

C## A R G U M E N T S:
                      INTEGER            IX,    N
                      DOUBLE PRECISION   X(N), SOLN(*)
C!!!!                 REAL               X(N), SOLN(*)

C## S T A T U S:
C               IGNORE LINES BEGINNING WITH  "C!!!!" .
C
C               SYSTEM  DEPENDENCE:                      NONE.
C               SINGLE/DOUBLE CONVERSION: NEEDED (SEE CONVRT).
C
C               THIS VERSION IS IN   D O U B L E   PRECISION.
C!!!!           THIS VERSION IS IN   S I N G L E   PRECISION.
C
C>RCS $HEADER: SOLN.F,V 2.1 91/11/20 10:53:07 BUCKLEY EXP $
C>RCS $LOG:     SOLN.F,V $
C>RCS REVISION 2.1  91/11/20  10:53:07  BUCKLEY
C>RCS FINAL SUBMISSION TO TOMS
C>RCS
C>RCS REVISION 2.0  90/07/31  11:34:22  BUCKLEY
C>RCS MINOR MOD TO USE BLAS
C>RCS
C>RCS REVISION 1.9  89/06/30  13:39:51  BUCKLEY
C>RCS PREPARING SUBMITTED VERSION OF MT
C>RCS
C>RCS REVISION 1.3.1.1  89/05/20  16:42:58  BUCKLEY
C>RCS TEMP. TEST OF MT BEFORE SUBMITTING
C>RCS
C>RCS REVISION 1.3  89/05/18  12:20:51  BUCKLEY
C>RCS FINAL TEST OF MT BEFORE SUBMITTING
C>RCS
C>RCS REVISION 1.2  89/05/15  21:00:03  BUCKLEY
C>RCS RECREATING INITIAL DEPOSIT FOR MT RCS.
C>RCS
C>RCS

C## D E S C R I P T I O N:
C
C    COMPUTE STANDARD SOLUTIONS FOR TEST FUNCTIONS AND PUT INTO X.
C    THE SOLUTION IS DEFINED ACCORDING TO THE CODE IX.
C
C      IF IX = 0, X(I) ARE DEFINED ON ENTRY FOR I=1,...,N.
C                ( THIS IS EQUIVALENT TO SPECIFYING IX = N. )
C      IF IX > 0, A CYCLE IS DEFINED AND MUST BE EXPANDED.
C      IF IX < 0  ON ENTRY, THEN A SELECTED FORMULA IS USED TO
C                  COMPUTE VALUES FOR THE COMPONENTS OF X.
C
C## E N T R Y   P O I N T S: THE NATURAL ENTRY ZZSOLN.

C## S U B R O U T I N E S:   NONE ARE CALLED.
C     REAL(DBLE)  ...  INTRINSIC
C     ABS         ...  INTRINSIC
C     MOD         ...  INTRINSIC
C     RD          ...  A STATEMENT FUNCTION FOR CONVERSION FROM INTEGER

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

C## L O C A L   D E C L:
                          INTEGER           I,     K
                          DOUBLE PRECISION  RD
C!!!!                     REAL              RD

C## S A V E:                 NONE SELECTED.
C## E Q U I V A L E N C E S: NONE ARE DEFINED.
C## C O M M O N:             NONE IS DEFINED.
C## D A T A:                 NO DATA VALUES SET.

C##                                           E X E C U T I O N
C##                                           E X E C U T I O N
C------            FUNCTION STATEMENT.
      RD(I) = DBLE(I)
C!!!! RD(I) = REAL(I)
C-----               NOW EXECUTE.

      IF ( IX .EQ. 0 ) THEN
         IX = N
      ENDIF

      IF ( IX .LT. 0 ) THEN
         GOTO ( 100 ) , ABS( IX )

C        ---SOLUTION 1 : X = (0, 0, . . ., N+1)
  100       DO 150 I = 1, N-1
               X(I) = ZERO
  150       CONTINUE
            X(N) = RD(N+1)
         GOTO 10000

      ELSE IF ( IX .GT. 0 ) THEN
         K = IX
         CALL ZZCOPY ( IX, SOLN, 1, X, 1 )

         DO 5000 I = IX+1, N
            K = MOD ( K, IX ) + 1
            X(I) = SOLN(K)
5000     CONTINUE
      ENDIF

10000 CONTINUE
      GOTO 90000

C## E X I T
90000      RETURN

C## F O R M A T S:  NONE ARE DEFINED.
C##                 E N D         OF ZZSOLN.
                    END
