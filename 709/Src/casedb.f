      SUBROUTINE ZZCASE (STRING, TYPE )

C## A R G U M E N T S:
                      CHARACTER *(*) STRING
                      INTEGER TYPE

C## S T A T U S:
C               SINGLE/DOUBLE CONVERSION:        NOT REQUIRED.
C
C               SYSTEM  DEPENDENCE:   SYSTEM ROUTINE FOR CASE
C                                     CONVERSION OF LETTERS.
C
C             THIS IS A VERSION FOR  SUN4

C>RCS $HEADER: CASE.GL,V 2.0 90/07/05 12:44:20 BUCKLEY EXP $
C>RCS $LOG:     CASE.GL,V $
C>RCS REVISION 2.0  90/07/05  12:44:20  BUCKLEY
C>RCS COMMON VERSION FOR TOMS AND GL
C>RCS
C>RCS REVISION 1.9.1.1  89/06/30  14:59:19  BUCKLEY
C>RCS PREPARING SUBMITTED VERSION OF MT
C>RCS
C>RCS REVISION 1.9  89/06/30  13:30:11  BUCKLEY
C>RCS PREPARING SUBMITTED VERSION OF MT
C>RCS
C>RCS REVISION 1.3.1.1  89/05/20  16:07:53  BUCKLEY
C>RCS TEMP. TEST OF MT BEFORE SUBMITTING
C>RCS
C>RCS REVISION 1.3  89/05/20  13:48:39  BUCKLEY
C>RCS FINAL TEST OF MT BEFORE SUBMITTING
C>RCS

C## D E S C R I P T I O N:

C       THIS CONVERTS EACH LOWER CASE ALPHABETIC LETTER TO
C       UPPER CASE, OR VICE VERSA.
C       IF TYPE = CTOUPP, CONVERSION IS LOWER TO UPPER
C       IF TYPE = CTOLOW, CONVERSION IS UPPER TO LOWER
C       IF TYPE = CTOCAP, USE UPPER FOR FIRST LETTER; LOWER FOR REST
C       ALL OTHER CHARACTERS ARE LEFT UNCHANGED.

C## E N T R Y   P O I N T S: THE NATURAL ENTRY TTOUPPR.
C## S U B R O U T I N E S:   LEN (INTRINSIC).
C## P A R A M E T E R S:

      INTEGER     CTOUPP,     CTOLOW,     CTOCAP
      PARAMETER ( CTOUPP = 1, CTOLOW = 2, CTOCAP = 3 )

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
                        INTEGER      I, SHIFT
                        LOGICAL      FIRST
                        CHARACTER *1 CH
C## S A V E:
                        SAVE FIRST, SHIFT

C## E Q U I V A L E N C E S: NONE ARE DEFINED.
C## C O M M O N:             NONE IS DEFINED.

C## D A T A:
                        DATA FIRST/.TRUE./

C##                                                  E X E C U T I O N
C##                                                  E X E C U T I O N
        IF (FIRST) THEN
                FIRST = .FALSE.
                SHIFT = ICHAR('A') - ICHAR('a')
        ENDIF

        I = 0
  100   I = I + 1
        IF ( I .LE. LEN(STRING) ) THEN
            IF ( TYPE .EQ. CTOUPP ) THEN
               IF ('a' .LE. STRING(I:I) .AND. STRING(I:I) .LE. 'z') THEN
                  CH = CHAR( ICHAR(STRING(I:I)) + SHIFT  )
               ELSE
                  CH = STRING(I:I)
               ENDIF
            ELSE IF ( TYPE .EQ. CTOLOW .OR. TYPE .EQ. CTOCAP ) THEN
               IF ('A' .LE. STRING(I:I) .AND. STRING(I:I) .LE. 'Z') THEN
                  CH = CHAR( ICHAR(STRING(I:I)) - SHIFT  )
               ELSE
                  CH = STRING(I:I)
               ENDIF
            ENDIF
            STRING(I:I) = CH
            GOTO 100
        ENDIF
        IF ( TYPE .EQ. CTOCAP .and.
     -      'a' .LE. STRING(1:1) .AND. STRING(1:1) .LE. 'z') THEN
                CH = CHAR( ICHAR(STRING(1:1)) + SHIFT  )
                STRING(1:1) = CH
        ENDIF

C## E X I T
90000      RETURN

C## F O R M A T S:  NONE ARE DEFINED.
C##                 E N D         OF TOUPPR.
                    END
