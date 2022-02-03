      DOUBLE PRECISION FUNCTION ZZINNR ( N, U, V, NRMFLG, IW, RW, DW )

C## A R G U M E N T S:
                       INTEGER  N, IW(*)
                       LOGICAL  NRMFLG
                       REAL              U(N), V(N)
C!!!!                  DOUBLE PRECISION  U(N), V(N)
                       DOUBLE PRECISION  DW(*)
                       REAL              RW(*)

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
C>RCS $HEADER: INNR.F,V 2.1 91/11/20 10:46:59 BUCKLEY EXP $
C>RCS $LOG:     INNR.F,V $
C>RCS REVISION 2.1  91/11/20  10:46:59  BUCKLEY
C>RCS FINAL SUBMISSION TO TOMS
C>RCS
C>RCS REVISION 2.0  90/07/31  11:12:32  BUCKLEY
C>RCS ADDED REVISED BLAS. .
C>RCS
C>RCS REVISION 1.9  89/06/30  13:11:47  BUCKLEY
C>RCS PREPARING SUBMITTED VERSION OF MT
C>RCS
C>RCS REVISION 1.3  89/05/18  12:40:14  BUCKLEY
C>RCS FINAL TEST OF MT BEFORE SUBMITTING
C>RCS
C>RCS REVISION 1.2  89/05/15  14:49:22  BUCKLEY
C>RCS INITIAL INSTALLATION OF MT INTO RCS FORM..
C>RCS
C>RCS REVISION 1.1  89/01/07  14:35:52  BUCKLEY
C>RCS INITIAL REVISION
C>RCS
C
C## D E S C R I P T I O N:
C
C     THIS ROUTINE COMPUTES THE NORMAL EUCLIDEAN INNER PRODUCT
C     OF THE VECTORS U AND V.
C
C *** NOTE THAT THE RESULT PASSED BACK IS
C
C     *ALWAYS* DOUBLE PRECISION.
C
C     IF NRMFLG IS SET ON ENTRY, THEN THE 2-NORM OF U IS COMPUTED
C     BY CALLING ZZNRM2 TO DO THE COMPUTATION WITHOUT OVERFLOW. IN
C     THIS CASE, V IS IGNORED AND THE NORM IS COMPUTED IN SINGLE OR
C     DOUBLE PRECISION AS APPROPRIATE.
C
C## E N T R Y   P O I N T S: THE NATURAL ENTRY ZZINNR
C## S U B R O U T I N E S:   ZZNRM2   FOR NO OVERFLOW 2-NORMS
C                            DBLE     ...INTRINSIC
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

C## L O C A L   D E C L:
                        INTEGER   I
                        REAL              ZZNRM2
C!!!!                   DOUBLE PRECISION  ZZNRM2

C## S A V E:                 NONE SELECTED.
C## E Q U I V A L E N C E S: NONE ARE DEFINED.
C## C O M M O N:             NONE IS DEFINED.
C## D A T A:                 NONE ARE SET.
C##                                                E X E C U T I O N
C##                                                E X E C U T I O N
      IF ( NRMFLG ) THEN
         ZZINNR = DBLE( ZZNRM2( N, U ) )
      ELSE
         ZZINNR = ZERO

         DO 500 I = 1,N
            ZZINNR = ZZINNR + DBLE(U(I)) * DBLE(V(I))
  500    CONTINUE
      ENDIF

C## E X I T
90000       RETURN

C## F O R M A T S:  NONE ARE DEFINED.
C##                 E N D         OF ZZINNR.
                    END
