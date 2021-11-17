      SUBROUTINE BBCUBC ( T, F, FP, TA, FA, FPA, LEFT, RIGHT, X, INTER )

C## A R G U M E N T S:
                     LOGICAL  INTER
                     DOUBLE PRECISION  T,F,FP,TA,FA,FPA,LEFT,RIGHT,X
C!!!!                REAL              T,F,FP,TA,FA,FPA,LEFT,RIGHT,X

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
C>RCS $HEADER: CUBC.F,V 1.11 91/11/22 11:27:36 BUCKLEY EXP $
C>RCS $LOG:     CUBC.F,V $
C>RCS REVISION 1.11  91/11/22  11:27:36  BUCKLEY
C>RCS FINAL SUBMISSION TO TOMS
C>RCS
C>RCS REVISION 1.10  90/07/31  10:49:13  BUCKLEY
C>RCS ADDED REVISED BLAS.
C>RCS
C>RCS REVISION 1.9  89/06/30  13:25:08  BUCKLEY
C>RCS PREPARING SUBMITTED VERSION OF MT
C>RCS
C>RCS REVISION 1.3.1.1  89/05/20  17:15:26  BUCKLEY
C>RCS TEMP. TEST OF MT BEFORE SUBMITTING
C>RCS
C>RCS REVISION 1.3  89/05/18  12:39:13  BUCKLEY
C>RCS FINAL TEST OF MT BEFORE SUBMITTING
C>RCS
C>RCS REVISION 1.2  89/05/15  14:55:18  BUCKLEY
C>RCS INITIAL INSTALLATION OF MT INTO RCS FORM.
C>RCS
C>RCS REVISION 1.1  89/01/17  16:54:26  BUCKLEY
C>RCS INITIAL REVISION
C>RCS
C
C## D E S C R I P T I O N:
C
C     GIVEN THE POINTS T AND TA, ALONG WITH THE FUNCTION VALUES
C     F AND FA AND SLOPES FP AND FPA AT EACH POINT,  THIS ROUTINE
C     FINDS THE POINT  X  AT WHICH THE CUBIC FITTED TO THE DATA
C     HAS ITS MINIMUM.  THE VALUES LEFT AND RIGHT DEFINE AN
C     INTERVAL. IF THERE IS NO MINIMUM OR IF IT LIES OUTSIDE THE
C     INTERVAL, X IS RETURNED AS ONE OF THE END POINTS, AS APPROPRIATE.
C     INTER IS RETURNED AS TRUE IF THE VALUE X RETURNED IS EQUAL TO
C     THAT OBTAINED FROM THE FORMULA INTERPOLATION. THE INTERPOLATION
C     IS COMPUTED FOLLOWING DETAILS GIVEN BY LEMARECHAL.
C
C## E N T R Y   P O I N T S: BBCUBC  THE NATURAL ENTRY
C                            BBSCUB  TO SET THE TRACE.
C
C## S U B R O U T I N E S:
C
C     ABS, DBLE(REAL), MAX, MIN, SQRT... INTRINSIC
C     RD... A STATEMENT FUNCTION
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

      INTEGER     XEPS,     XSMALL,     XBIG
      PARAMETER ( XEPS = 1, XSMALL = 2, XBIG = 3 )

C## L O C A L   D E C L:

      INTEGER   TRU,    STRU
      LOGICAL   EXTREM, TRACE, STRACE, ORDER, ABIGGR, FIRST, PBIGGR

      DOUBLE PRECISION  P, DISC, DUMMY
C!!!! REAL              P, DISC, DUMMY

      DOUBLE PRECISION SGN, APR, BPR, NUM, XC, RD, EPS, BIGGST, ZZMPAR
C!!!! REAL             SGN, APR, BPR, NUM, XC, RD, EPS, BIGGST, ZZMPAR

      DOUBLE PRECISION ALEFT, ARIGHT
C!!!! REAL             ALEFT, ARIGHT

C## S A V E:
           SAVE TRU, TRACE, EPS, FIRST

C## E Q U I V A L E N C E S: NONE ARE DEFINED.
C## C O M M O N:             NONE IS DEFINED.
C## D A T A:
            DATA  TRU/6/, TRACE/.FALSE./, FIRST/.TRUE./
C##                                                E X E C U T I O N
C##                                                E X E C U T I O N
C----                         A STATEMENT FUNCTION.
      RD(DUMMY) = DBLE(DUMMY)
C!!!! RD(DUMMY) = REAL(DUMMY)
C----

      IF ( FIRST ) THEN
         EPS = SQRT(ZZMPAR(XBIG))
         FIRST = .FALSE.
      ENDIF

      ALEFT  = MIN(LEFT, RIGHT)
      ARIGHT = MAX(LEFT, RIGHT)

      IF ( TRACE ) THEN
         WRITE (TRU,*) ' [CUBC] T,F,FP, TA,FA,FPA->',  T,F,FP, TA,FA,FPA
         WRITE (TRU,*) ' [CUBC] INTERVAL [',ALEFT,',',ARIGHT,']'
      ENDIF

      EXTREM = .FALSE.
      ORDER  = LEFT .LE. RIGHT .EQV.  T .LE. TA
      SGN    = SIGN(ONE,TA-T)
      IF (TRACE) WRITE(TRU,*) ' [CUBC] ORDER->',ORDER,'  SGN->',SGN

      IF ( T .EQ. TA ) THEN
         IF (TRACE) WRITE(TRU,*) ' [CUBC] POINTS EQUAL.'
         X = T
         INTER = .FALSE.
      ELSE
         P  = DBLE(FP) + DBLE(FPA) - DBLE(THREE)*DBLE(FA-F)/DBLE(TA-T)

         IF ( SIGN(ONE,FPA) .NE. SIGN(ONE,FP) ) THEN
            DISC = DBLE(ONE) - (DBLE(FP)/P)*(DBLE(FPA/P))
            DISC = ABS(P)*SQRT(DISC)
         ELSE
            IF (TRACE) WRITE(TRU,*) ' [CUBC] SIGN(FP)=SIGN(FPA).'
            BIGGST = MAX(ABS(FP),ABS(FPA),ABS(P))
            ABIGGR = BIGGST .EQ. ABS(FPA)
            PBIGGR = BIGGST .EQ. ABS( P )
            IF(TRACE)WRITE(TRU,*) ' [CUBC] P,BIGGST,EPS->',P,BIGGST,EPS
            IF (BIGGST .LE. EPS) THEN
               DISC = DBLE(P**2) - DBLE(FP)*DBLE(FPA)
               IF (TRACE) WRITE(TRU,*) ' [CUBC] P,DISC->', P, DISC
            ELSE IF ( PBIGGR ) THEN
               DISC = DBLE(P) - (DBLE(FPA)/DBLE(P))*DBLE(FP)
            ELSE IF ( ABIGGR ) THEN
               DISC = (DBLE(P)/DBLE(FPA))*DBLE(P) - DBLE(FP)
            ELSE
               DISC = (DBLE(P)/DBLE(FP))*DBLE(P) - DBLE(FPA)
            ENDIF
            IF (TRACE) WRITE(TRU,*) ' [CUBC] DISC->', DISC
            IF ( DISC .GE. 0 ) THEN
               IF (BIGGST .LE. EPS) THEN
                  DISC = SQRT(DISC)
               ELSE
                  DISC = SQRT(DISC)*SQRT(BIGGST)
               ENDIF
               IF (TRACE) WRITE(TRU,*) ' [CUBC] DISC->', DISC
            ELSE
               INTER = .FALSE.
               IF ( FP .LT. ZERO ) THEN
                     X = ARIGHT
               ELSE
                     X = ALEFT
               ENDIF
               IF (TRACE) WRITE(TRU,*) ' [CUBC] NO MINIMUM!'
               GOTO 90000
            ENDIF

         ENDIF

            DISC = SGN*DISC
            IF (TRACE) WRITE(TRU,*) ' [CUBC] DISC->',DISC

            APR = DBLE(FP) + DBLE(FPA) + DBLE(TWO*P)
            BPR = DBLE(FP) + DBLE(P)
            IF (TRACE) WRITE(TRU,*) ' [CUBC] APR,BPR->',APR,BPR

            IF ( SGN*BPR .LT. ZERO ) THEN
               IF (TRACE) WRITE(TRU,*) ' [CUBC] USING REGULAR FORM.'
               X = T + FP*(TA-T)/RD(BPR-DISC)
               IF (TRACE) WRITE(TRU,*) ' [CUBC]  PREDICT X->',X
            ELSE
               NUM = DISC + BPR
               IF (TRACE) WRITE(TRU,*) ' [CUBC] USING ALTERNATE FORM.'
               IF (TRACE) WRITE(TRU,*) ' [CUBC]  NUM->',NUM
               IF ( ABS((T-TA)*NUM) .GE. (ARIGHT-ALEFT)*ABS(APR) ) THEN
                  X = ARIGHT
                  EXTREM = .TRUE.
                  IF (TRACE) WRITE(TRU,*) ' [CUBC] CUT OFF TO X->',X
               ELSE
                  X = T + NUM*(TA-T)/APR
                  IF (TRACE) WRITE(TRU,*) ' [CUBC] PREDICT X->',X
               ENDIF
            ENDIF

         XC = X
         X  = MAX(X,ALEFT )
         X  = MIN(X,ARIGHT)

         INTER = .NOT. EXTREM .AND. XC .EQ. X

         IF (TRACE) WRITE(TRU,*) ' [CUBC] X,XC,INTER,EXTREM->',
     -                                 X,XC,INTER,EXTREM
      ENDIF

      GOTO 90000

C## E N T R Y  BBSCUB:
                      ENTRY  BBSCUB (STRACE,STRU)
         TRACE = STRACE
         TRU   = STRU
         RETURN

C## E X I T
90000       RETURN

C## F O R M A T S:  NONE ARE DEFINED.
C##                 E N D         OF BBCUBC.
                    END
