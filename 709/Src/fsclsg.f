      SUBROUTINE ZZFSCL ( FT, FV, SCALE, FSCALE, FONLY, GONLY )

C## A R G U M E N T S:
                       INTEGER           FSCALE
                       LOGICAL           FONLY, GONLY

                       REAL              FT, FV, SCALE
C!!!!                  DOUBLE PRECISION  FT, FV, SCALE

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
C>RCS $HEADER: FSCL.F,V 2.1 91/11/20 10:46:58 BUCKLEY EXP $
C>RCS $LOG:     FSCL.F,V $
C>RCS REVISION 2.1  91/11/20  10:46:58  BUCKLEY
C>RCS FINAL SUBMISSION TO TOMS
C>RCS
C>RCS REVISION 2.0  90/07/31  11:12:20  BUCKLEY
C>RCS ADDED REVISED BLAS.
C>RCS
C>RCS REVISION 1.9  89/06/30  13:11:50  BUCKLEY
C>RCS PREPARING SUBMITTED VERSION OF MT
C>RCS
C>RCS REVISION 1.3  89/05/18  12:40:17  BUCKLEY
C>RCS FINAL TEST OF MT BEFORE SUBMITTING
C>RCS
C>RCS REVISION 1.2  89/05/15  14:49:30  BUCKLEY
C>RCS INITIAL INSTALLATION OF MT INTO RCS FORM..
C>RCS
C>RCS REVISION 1.1  89/01/07  14:36:08  BUCKLEY
C>RCS INITIAL REVISION
C>RCS
C
C## D E S C R I P T I O N:
C
C     THIS SUBROUTINE APPLIES ONE OF SEVERAL SCALINGS (LINEAR OR
C     NONLINEAR) TO A FUNCTION VALUE.
C
C-----ON ENTRY:
C
C              FT - THE PRESENT FUNCTION VALUE
C
C              FSCALE - THE CODE FOR THE TYPE OF SCALE DESIRED. WHERE
C                    THE SCALE FUNCTION IS ONE OF THE FOLLOWING:
C
C                     1:   F(Z) = 1 + Z
C                     2:   F(Z) = Z*Z
C                     3:   F(Z) = -1 / (1 + Z*Z)
C                     4:   F(Z) =  SQRT(1 + Z*Z)
C                     5:   F(Z) = Z*Z*Z
C
C              FONLY - IF TRUE ONLY THE FUNCTION IS EVALUATED.
C              GONLY - IF TRUE ONLY THE GRADIENT IS EVALUATED.
C
C-----ON EXIT:
C              FV - THE SCALED FUNCTION VALUE.
C              SCALE - GRADIENT SCALING FACTOR.
C
C## E N T R Y   P O I N T S: THE NATURAL ENTRY ZZFSCL.
C## S U B R O U T I N E S:   SQRT...  INTRINSIC
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

C## L O C A L   D E C L:     NONE ARE DEFINED.
C## S A V E:                 NONE SELECTED.
C## E Q U I V A L E N C E S: NONE ARE DEFINED.
C## C O M M O N:             NONE IS DEFINED.
C## D A T A:                 NONE ARE SET.
C##                                                E X E C U T I O N
C##                                                E X E C U T I O N

      GOTO (2100,2200,2300,2400,2500), FSCALE

C     -----FF(Z) = 1 + F(Z) -------FSCALE = 1.

 2100    IF ( .NOT. GONLY ) FV    = FT + ONE
         IF ( .NOT. FONLY ) SCALE = ONE
         GOTO 90000

C     -----FF(Z) = Z*Z ------------FSCALE = 2.

 2200    IF ( .NOT. GONLY ) FV    =  FT * FT
         IF ( .NOT. FONLY ) SCALE = TWO * FT
         GOTO 90000

C     -----FF(Z) = -1/(1+Z**2) --- FSCALE = 3.

 2300                       FV    = -ONE / ( ONE + FT**2 )
         IF ( .NOT. FONLY ) SCALE =  TWO * FT * FV**2
         GOTO 90000

C     -----FF(Z) = SQRT(1+Z**2) -- FSCALE = 4.

 2400                       FV    =  SQRT(ONE + FT**2)
         IF ( .NOT. FONLY ) SCALE =  FT/FV
         GOTO 90000

C     -----FF(Z) = Z*Z*Z --------- FSCALE = 5.

 2500    IF ( .NOT. GONLY ) FV    =    FT*FT*FT
         IF ( .NOT. FONLY ) SCALE = THREE*FT*FT
         GOTO 90000

C## E X I T
90000       RETURN

C## F O R M A T S:  NONE ARE DEFINED.
C##                 E N D         OF ZZFSCL.
                    END
