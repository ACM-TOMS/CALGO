       SUBROUTINE ZZDSOL ( N, F, X, W, SOLNS, DF, DX, SOLNF, SOLNX )

C## A R G U M E N T S:
                      INTEGER           N, SOLNF, SOLNX

                      DOUBLE PRECISION  F, X(N), W(*), SOLNS(*), DF, DX
C!!!!                 REAL              F, X(N), W(*), SOLNS(*), DF, DX

C## S T A T U S:
C               SINGLE/DOUBLE CONVERSION: NEEDED (SEE CONVRT).
C
C               IGNORE LINES BEGINNING WITH  "C!!!!" .
C
C               THIS VERSION IS IN   D O U B L E   PRECISION.
C!!!!           THIS VERSION IS IN   S I N G L E   PRECISION.
C
C               SYSTEM  DEPENDENCE:                      NONE.

C>RCS $HEADER: DSOL.F,V 2.1 91/11/20 10:52:43 BUCKLEY EXP $
C>RCS $LOG:     DSOL.F,V $
C>RCS REVISION 2.1  91/11/20  10:52:43  BUCKLEY
C>RCS FINAL SUBMISSION TO TOMS
C>RCS
C>RCS REVISION 2.0  90/07/31  11:24:09  BUCKLEY
C>RCS ADDED BLAS CALL.
C>RCS
C>RCS REVISION 1.9  89/06/30  13:40:52  BUCKLEY
C>RCS PREPARING SUBMITTED VERSION OF MT
C>RCS
C>RCS REVISION 1.1  89/01/17  16:50:08  BUCKLEY
C>RCS INITIAL REVISION
C>RCS

C## D E S C R I P T I O N:

C## E N T R Y   P O I N T S: THE NATURAL ENTRY ZZDSOL.

C## S U B R O U T I N E S: ABS, MAX, NINT   INTRINSIC
C                          ZZVECT           GENERATE VECTORS

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

C## L O C A L   D E C L:
                        INTEGER           I, K, NSOLNS, PT, NCYC, DIM
                        DOUBLE PRECISION  Z,    ZZNRM2, ARGS, TF
C!!!!                   REAL              Z,    ZZNRM2, ARGS, TF

C## S A V E:                 NONE SELECTED.
C## E Q U I V A L E N C E S: NONE ARE DEFINED.
C## C O M M O N:             NONE IS DEFINED.
C## D A T A:                 NO DATA VALUES SET.

C## E X E C U T I O N
C## E X E C U T I O N

      NSOLNS = NINT(SOLNS(1))

      PT = 2
      DF = -ONE
      DX = -ONE
      SOLNX = 0
      SOLNF = 0

      DO 500 I = 1,NSOLNS

         DIM = NINT (SOLNS(PT+1))
         IF ( DIM .EQ. N .OR. DIM .EQ. 0 ) THEN
            IF ( DF .EQ. -ONE ) THEN
               DF = ABS ( F-SOLNS(PT+2) )
               TF = SOLNS(PT+2)
               SOLNF = I
            ELSE IF ( ABS (F-SOLNS(PT+2)) .LT. DF ) THEN
               DF = ABS ( F-SOLNS(PT+2) )
               TF = SOLNS(PT+2)
               SOLNF = I
            ENDIF

            NCYC = NINT (SOLNS(PT))

            IF ( NCYC .NE. 0 ) THEN
               DO 190 K = 1,NCYC
                  W(K) = SOLNS(PT+2+K)
  190          CONTINUE
C
C  IF THE VALUE OF NCYC IS 1, THEN THAT ONE NUMBER IN THE CYCLE SPECIFIES
C  THE SPECIAL VECTOR IN ZZVECT TO BE USED AS A SOLUTION VECTOR.
C
               IF ( NCYC .EQ. 1 ) NCYC = -NINT( SOLNS(PT+3) )
C
C  SET ARGS FOR SOLUTIONS IN WHICH MORE THAN ONE COMPONENT HAS THE
C  VALUE X(1).  I.E. FUNCTIONS BRWNAL AND BOX66.
C
               IF ( NCYC .EQ. -16  .OR.
     -              NCYC .EQ. -17  .OR.
     -              NCYC .EQ. -18       ) ARGS = X(1)
C
C  SET ARGS FOR SOLUTIONS IN WHICH MORE THAN ONE COMPONENT HAS THE
C  VALUE X(3).  I.E. FUNCTION TOIN4.
C
               IF ( NCYC .EQ. -19 ) ARGS = X(3)
C
               CALL ZZVECT ( NCYC, N, SOLNS(PT+3), W, ARGS )
               CALL ZZAXPY ( N, -ONE, X, 1, W, 1 )
               Z = ZZNRM2 ( N, W )
               IF ( DX .EQ. -ONE ) THEN
                  DX = Z
                  SOLNX = I
               ELSE IF ( Z .LT. DX ) THEN
                  DX = Z
                  SOLNX = I
               ENDIF
            ENDIF

         ENDIF

         PT = PT + 3 + MAX(NCYC,0)

  500 CONTINUE

      IF ( DF .NE. -ONE ) F = TF
      GOTO 90000

C## E X I T
90000      RETURN

C## F O R M A T S:  NONE ARE DEFINED.
C##                 E N D         OF ZZDSOL.
                    END
