      SUBROUTINE ZZCNTR ( STRING )

C## A R G U M E N T S:
                      CHARACTER *(*) STRING

C## S T A T U S:
C               SINGLE/DOUBLE CONVERSION:        NOT REQUIRED.
C               SYSTEM  DEPENDENCE:                      NONE.
C
C>RCS $HEADER: CNTR.F,V 1.10 91/11/19 16:14:31 BUCKLEY EXP $
C>RCS $LOG:     CNTR.F,V $
C>RCS REVISION 1.10  91/11/19  16:14:31  BUCKLEY
C>RCS FINAL SUBMISSION TO TOMS
C>RCS
C>RCS REVISION 1.9  89/06/30  13:28:52  BUCKLEY
C>RCS PREPARING SUBMITTED VERSION OF MT
C>RCS
C>RCS REVISION 1.3.1.1  89/05/20  13:39:28  BUCKLEY
C>RCS TEMP. TEST OF MT BEFORE SUBMITTING
C>RCS
C>RCS REVISION 1.3  89/05/18  14:26:40  BUCKLEY
C>RCS FINAL TEST OF MT BEFORE SUBMITTING
C>RCS
C>RCS REVISION 1.2  89/05/15  14:44:13  BUCKLEY
C>RCS INITIAL INSTALLATION OF MT INTO RCS FORM.
C>RCS
C>RCS REVISION 1.1  89/01/17  16:30:07  BUCKLEY
C>RCS INITIAL REVISION
C>RCS
C
C## D E S C R I P T I O N:
C
C     THIS ROUTINE SHIFTS THE NONBLANK CHARACTERS OF STRING SO THAT
C     THERE IS A BALANCE OF BLANKS ON LEFT AND RIGHT.
C
C## E N T R Y   P O I N T S: THE NATURAL ENTRY ZZCNTR
C## S U B R O U T I N E S:   LEN    ...INTRINSIC
C                            ZZSHFT ...TO SHIFT A STRING.
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
                            INTEGER START,  ENDCH, SLEN, CLEN, LEFT, I

C## S A V E:                 NONE SELECTED.
C## E Q U I V A L E N C E S: NONE ARE DEFINED.
C## C O M M O N:             NONE IS DEFINED.
C## D A T A:                 NONE ARE SET.
C##                                                E X E C U T I O N
C##                                                E X E C U T I O N
      SLEN = LEN(STRING)

C                          FIND FIRST NONBLANK CHARACTER.
      DO 1000 I=1,SLEN
         IF ( STRING(I:I) .NE. BLANK ) THEN
            START = I
            GOTO 2000
         ENDIF
 1000 CONTINUE

C     THE STRING IS ALL BLANK IF WE FALL THROUGH THE DO LOOP.
      GOTO 90000

C                           FIND LAST NON-BLANK CHARACTER.
 2000 DO 3000 I=SLEN,1,-1
         IF ( STRING(I:I) .NE. BLANK ) THEN
            ENDCH = I
            GOTO 4000
         ENDIF
 3000 CONTINUE

C                          COMPUTE SHIFT AND DO IT.
 4000 CLEN = ENDCH - START + 1
      LEFT = 1 + (SLEN - CLEN) / 2

      CALL ZZSHFT ( STRING, START, LEFT, CLEN )
      GOTO 90000

C## E X I T
90000      RETURN

C## F O R M A T S:  NONE ARE DEFINED.
C##                 E N D         OF ZZCNTR.
                    END
