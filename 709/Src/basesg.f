      SUBROUTINE ZZBASE ( N, BASE, STRING,  * )

C## A R G U M E N T S:
                     CHARACTER *(*)  STRING
                     INTEGER         N,        BASE

C## S T A T U S:
C               SINGLE/DOUBLE CONVERSION:        NOT REQUIRED.
C
C>RCS $HEADER: BASE.F,V 1.10 91/11/19 16:04:05 BUCKLEY EXP $
C>RCS $LOG:     BASE.F,V $
C>RCS REVISION 1.10  91/11/19  16:04:05  BUCKLEY
C>RCS FINAL SUBMISSION TO TOMS
C>RCS
C>RCS REVISION 1.9  89/06/30  13:28:51  BUCKLEY
C>RCS PREPARING SUBMITTED VERSION OF MT
C>RCS
C>RCS REVISION 1.3.1.1  89/05/20  13:39:27  BUCKLEY
C>RCS TEMP. TEST OF MT BEFORE SUBMITTING
C>RCS
C>RCS REVISION 1.3  89/05/18  14:26:38  BUCKLEY
C>RCS FINAL TEST OF MT BEFORE SUBMITTING
C>RCS
C>RCS REVISION 1.2  89/05/15  14:44:03  BUCKLEY
C>RCS INITIAL INSTALLATION OF MT INTO RCS FORM.
C>RCS
C>RCS REVISION 1.1  89/01/17  16:30:07  BUCKLEY
C>RCS INITIAL REVISION
C>RCS
C
C## D E S C R I P T I O N:
C
C     THIS ROUTINE TAKES THE STRING OF CHARACTERS IN STRING AND
C     CONVERTS THEM TO AN INTEGER VALUE IN N.  THE NUMBERS ARE ASSUMED
C     TO BE IN BASE "BASE", WHERE 2 <= BASE <= 16.
C
C     THE STRING IS PROCESSED LEFT TO RIGHT.  BLANKS ARE IGNORED,
C     EVEN WITHIN THE NUMBER. A LEADING + OR - MAY BE GIVEN.  ALL
C     OTHER CHARACTERS MUST BE BETWEEN  0  AND BASE-1.  IF BASE > 10,
C     THEN THE ADDITIONAL DIGITS ALLOWED ARE AS MANY OF A,B,C,D,E,F
C     AS ARE NECESSARY.  ONLY UPPERCASE LETTERS A TO Z ARE RECOGNIZED.
C     ANY ERROR CAUSES THE ALTERNATE EXIT TO BE TAKEN.
C
C## E N T R Y   P O I N T S: THE NATURAL ENTRY ZZBASE
C## S U B R O U T I N E S:   INDEX, LEN   ... INTRINSIC
C## P A R A M E T E R S:

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

C## L O C A L   D E C L:
                             CHARACTER * 1  CH
                             CHARACTER *16  DIGITS
                             INTEGER        I, K, SIGN

C## S A V E:
                             SAVE  DIGITS

C## E Q U I V A L E N C E S: NONE ARE DEFINED.
C## C O M M O N:             NONE IS DEFINED.
C# D A T A:
          DATA   DIGITS  / '0123456789ABCDEF' /
C##                                                E X E C U T I O N
C##                                                E X E C U T I O N
      N    = 0
      SIGN = 0

      DO 1000 I = 1,LEN(STRING)

         CH = STRING(I:I)
         IF ( CH .NE. BLANK ) THEN
            IF      ( CH .EQ. PLUS ) THEN
               IF ( SIGN .EQ. 0 ) THEN
                  SIGN = 1
               ELSE
                  GOTO 91000
               ENDIF
            ELSE IF ( CH .EQ. MINUS ) THEN
               IF ( SIGN .EQ. 0 ) THEN
                  SIGN = -1
               ELSE
                  GOTO 91000
               ENDIF
            ELSE
               K = INDEX ( DIGITS(1:BASE), CH )
               IF ( K .EQ. 0 ) THEN
                  GOTO 91000
               ELSE
                  N    = N*BASE + K - 1
                  IF ( SIGN .EQ. 0 ) THEN
                     SIGN = 1
                  ENDIF
               ENDIF
            ENDIF
C                 FOR CASE-SOLUTION
         ENDIF
C               FOR NON-BLANK CHARACTER
 1000 CONTINUE

      IF ( SIGN .EQ. 0 ) THEN
         SIGN = 1
      ENDIF
      N = N * SIGN
      GOTO 90000

C## E X I T
90000      RETURN
91000      RETURN 1

C## F O R M A T S:  NONE ARE DEFINED.
C##                 E N D         OF ZZBASE.
                    END
