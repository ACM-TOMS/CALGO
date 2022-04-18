      INTEGER FUNCTION ZZLENG (LINE)

C## A R G U M E N T S:
                       CHARACTER*(*) LINE

C## S T A T U S:
C               SINGLE/DOUBLE CONVERSION:        NOT REQUIRED.
C               SYSTEM  DEPENDENCE:                      NONE.
C
C>RCS $HEADER: LENG.F,V 1.10 91/11/19 16:13:33 BUCKLEY EXP $
C>RCS $LOG:     LENG.F,V $
C>RCS REVISION 1.10  91/11/19  16:13:33  BUCKLEY
C>RCS FINAL SUBMISSION TO TOMS
C>RCS
C>RCS REVISION 1.9  89/06/30  13:28:59  BUCKLEY
C>RCS PREPARING SUBMITTED VERSION OF MT
C>RCS
C>RCS REVISION 1.3.1.1  89/05/20  13:39:34  BUCKLEY
C>RCS TEMP. TEST OF MT BEFORE SUBMITTING
C>RCS
C>RCS REVISION 1.3  89/05/18  14:26:46  BUCKLEY
C>RCS FINAL TEST OF MT BEFORE SUBMITTING
C>RCS
C>RCS REVISION 1.2  89/05/15  14:44:31  BUCKLEY
C>RCS INITIAL INSTALLATION OF MT INTO RCS FORM.
C>RCS
C>RCS REVISION 1.1  89/01/17  16:30:09  BUCKLEY
C>RCS INITIAL REVISION
C>RCS
C
C## D E S C R I P T I O N:
C
C     THIS ROUTINE DETERMINES THE POSITION OF THE LAST NONBLANK
C     CHARACTER IN THE STRING LINE. IF THE LINE IS ENTIRELY
C     BLANK, THEN ZZLENG IS SET TO 0.
C
C## E N T R Y   P O I N T S: THE NATURAL ENTRY ZZLENG
C## S U B R O U T I N E S:   LEN    ...INTRINSIC
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
                             INTEGER   I

C## S A V E:                 NONE SELECTED.
C## E Q U I V A L E N C E S: NONE ARE DEFINED.
C## C O M M O N:             NONE IS DEFINED.
C## D A T A:                 NONE ARE SET.
C##                                                E X E C U T I O N
C##                                                E X E C U T I O N
      ZZLENG = 0

      DO 1000  I = LEN(LINE), 1, -1
         IF ( LINE(I:I) .NE. BLANK ) THEN
            ZZLENG = I
            GOTO 90000
         ENDIF
 1000 CONTINUE

C## E X I T
90000      RETURN

C## F O R M A T S:  NONE ARE DEFINED.
C##                 E N D         OF ZZLENG.
                    END
