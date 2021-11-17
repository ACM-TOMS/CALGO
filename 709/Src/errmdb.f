      SUBROUTINE ZZERRM ( VALUE, *, MESSAG )

C## A R G U M E N T S:
                      CHARACTER *(*)    MESSAG

                      DOUBLE PRECISION  VALUE
C!!!!                 REAL              VALUE

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
C>RCS $HEADER: ERRM.F,V 1.10 91/11/19 16:15:44 BUCKLEY EXP $
C>RCS $LOG:     ERRM.F,V $
C>RCS REVISION 1.10  91/11/19  16:15:44  BUCKLEY
C>RCS FINAL SUBMISSION TO TOMS
C>RCS
C>RCS REVISION 1.9  89/06/30  13:28:55  BUCKLEY
C>RCS PREPARING SUBMITTED VERSION OF MT
C>RCS
C>RCS REVISION 1.3.1.1  89/05/20  13:39:31  BUCKLEY
C>RCS TEMP. TEST OF MT BEFORE SUBMITTING
C>RCS
C>RCS REVISION 1.3  89/05/18  14:26:43  BUCKLEY
C>RCS FINAL TEST OF MT BEFORE SUBMITTING
C>RCS
C>RCS REVISION 1.2  89/05/15  14:44:23  BUCKLEY
C>RCS INITIAL INSTALLATION OF MT INTO RCS FORM.
C>RCS
C>RCS REVISION 1.1  89/01/17  16:30:08  BUCKLEY
C>RCS INITIAL REVISION
C>RCS
C
C## D E S C R I P T I O N:
C
C     THIS ROUTINE PRINTS AN ERROR MESSAGE, AND THEN RETURNS.
C
C     CODE   IS MESSAGE ( 1:1 )
C     LEVEL  IS MESSAGE ( 2:2 )
C
C     IF  CODE = NONE ('N')    VALUE IS IGNORED.
C              = INTG ('I')    THE INTEGER NINT(VALUE) IS PRINTED.
C              = REEL ('R')    VALUE IS PRINTED.
C
C     LEVEL  GIVES THE SEVERITY CODE, AS FOLLOWS:
C
C           TRIVAL ('T')  NORMAL RETURN IN ANY CASE; A TRIVIAL ERROR.
C           SEVERE ('S')  NORMAL RETURN IF INTACT IS TRUE;
C                         ALTERNATE RETURN OTHERWISE.
C                         THIS ALLOWS FOR AN ABORT IN BATCH MODE.
C           FATAL  ('F')  ALTERNATE RETURN IN ANY CASE.  THIS FORCES AN
C                         ABORT EVEN IN INTERACTIVE MODE.
C
C## E N T R Y   P O I N T S: ZZERRM   THE NATURAL ENTRY.
C                            ZZETRM   DEFINE UNIT.
C
C## S U B R O U T I N E S:   NINT   ...INTRINSIC
C                            ZZLENG ...NON-BLANK LENGTH OF A STRING.
C
C## P A R A M E T E R S:     NONE ARE DEFINED.

      CHARACTER*(*) BLANK,        QUOTE,        HASH
      PARAMETER (   BLANK  = ' ', QUOTE  = '"', HASH   = '#' )

      CHARACTER*(*) PERIOD,       COMMA,        SEMICN
      PARAMETER (   PERIOD = '.', COMMA  = ',', SEMICN = ';' )

      CHARACTER*(*) COLON,        DASH,         EQUALS
      PARAMETER (   COLON  = ':', DASH   = '-', EQUALS = '=' )

      CHARACTER*(*) OBRACE,       CBRACE,       UNDERS
      PARAMETER (   OBRACE = '{', CBRACE = '}', UNDERS = '_' )

      CHARACTER*(*) PLUS,         MINUS,        EXCLAM
      PARAMETER (   PLUS   = '+', MINUS  = '-', EXCLAM = '!' )

      CHARACTER*(*) GTHAN,        LTHAN,        QUESMK
      PARAMETER (   GTHAN  = '>', LTHAN  = '<', QUESMK = '?' )

      CHARACTER*(*) SLASH,        BSLASH,       PERCNT
      PARAMETER (   SLASH  = '/', BSLASH = '\\',PERCNT = '%' )

      CHARACTER*(*) CARAT,        ATSIGN,       TILDE
      PARAMETER (   CARAT  = '^', ATSIGN = '@', TILDE = '~' )
             INTEGER         MESTRT
             PARAMETER (     MESTRT = 3 )

C## L O C A L   D E C L:

                        EXTERNAL        ZZLENG
                        INTEGER   K,    ZZLENG, OUT, DFUNIT
                        LOGICAL   MODE, INTACT

                        CHARACTER * 1   TRIVAL, SEVERE, FATAL,  CODE
                        CHARACTER * 1   INTG,   REEL,   NONE,   LEVEL

C## S A V E:
                        SAVE TRIVAL, SEVERE, FATAL, INTACT
                        SAVE    OUT,   NONE,  INTG,   REEL

C## E Q U I V A L E N C E S: NONE ARE DEFINED.
C## C O M M O N:             NONE IS DEFINED.
C## D A T A:
             DATA  OUT    /6/,   INTACT /.FALSE./
             DATA  TRIVAL /'T'/, SEVERE /'S'/, FATAL  /'F'/
             DATA  NONE   /'N'/, INTG   /'I'/, REEL   /'R'/
C##                                                E X E C U T I O N
C##                                                E X E C U T I O N
C-----                    FIRST WRITE THE ERROR MESSAGE.

      K = ZZLENG ( MESSAG )
      IF ( MESSAG(K:K) .EQ. PERIOD ) THEN
         K = K - 1
      ENDIF

      CODE  = MESSAG ( 1:1 )
      LEVEL = MESSAG ( 2:2 )

      IF      ( CODE .EQ. NONE ) THEN
         WRITE ( OUT, 99999 ) BLANK, MESSAG(MESTRT:K), PERIOD
      ELSE IF ( CODE .EQ. INTG ) THEN
         WRITE ( OUT, 99998 ) BLANK, MESSAG(MESTRT:K),
     -                        BLANK, NINT(VALUE), PERIOD
      ELSE IF ( CODE .EQ. REEL ) THEN
         WRITE ( OUT, 99997 ) BLANK, MESSAG(MESTRT:K),
     -                        BLANK, VALUE, PERIOD
      ENDIF

C-----THEN INDICATE WHETHER JOB CONTINUING OR ABORTING.

      IF ( INTACT ) THEN
         IF ( LEVEL .EQ. FATAL  ) THEN
            WRITE ( OUT, 99999 ) ' JOB BEING ABORTED.'
         ENDIF
      ELSE
         IF ( LEVEL .EQ. TRIVAL ) THEN
            WRITE ( OUT, 99999 ) ' JOB WILL CONTINUE.'
         ELSE
            WRITE ( OUT, 99999 ) ' JOB BEING ABORTED.'
         ENDIF
      ENDIF
      GOTO 90000

C## E N T R Y  ZZETRM:
                         ENTRY  ZZETRM ( MODE, DFUNIT )
      OUT    = DFUNIT
      INTACT = MODE
      RETURN

C## E X I T
C               USE ALTERNATE RETURN IF ABORT REQUIRED.

90000 IF ( (       INTACT .AND. LEVEL .EQ. FATAL  ) .OR.
     -       .NOT. INTACT .AND. LEVEL .NE. TRIVAL )      THEN
         RETURN 1
      ELSE
         RETURN
      ENDIF

C## F O R M A T S:

99997 FORMAT ( A, A, A, G14.7, A )
99998 FORMAT ( A, A, A, I8, A )
99999 FORMAT ( A, A, A )

C##                 E N D         OF ZZERRM.
                    END
