      SUBROUTINE ZZSHFT (STRING, FROM, TO, NUMBER )

C## A R G U M E N T S:
                      INTEGER        FROM,    TO,    NUMBER
                      CHARACTER *(*) STRING
C## S T A T U S:
C               SYSTEM  DEPENDENCE:                      NONE.
C               SINGLE/DOUBLE CONVERSION:        NOT REQUIRED.
C
C>RCS $HEADER: SHFT.F,V 1.10 91/11/19 16:18:04 BUCKLEY EXP $
C>RCS $LOG:     SHFT.F,V $
C>RCS REVISION 1.10  91/11/19  16:18:04  BUCKLEY
C>RCS FINAL SUBMISSION TO TOMS
C>RCS
C>RCS REVISION 1.9  89/06/30  13:29:07  BUCKLEY
C>RCS PREPARING SUBMITTED VERSION OF MT
C>RCS
C>RCS REVISION 1.3.1.1  89/05/20  13:39:42  BUCKLEY
C>RCS TEMP. TEST OF MT BEFORE SUBMITTING
C>RCS
C>RCS REVISION 1.3  89/05/18  14:26:57  BUCKLEY
C>RCS FINAL TEST OF MT BEFORE SUBMITTING
C>RCS
C>RCS REVISION 1.2  89/05/15  14:44:49  BUCKLEY
C>RCS INITIAL INSTALLATION OF MT INTO RCS FORM.
C>RCS
C>RCS REVISION 1.1  89/01/17  16:30:10  BUCKLEY
C>RCS INITIAL REVISION
C>RCS
C## D E S C R I P T I O N:

C     THIS ROUTINE PERFORMS A SHIFT OF CHARACTERS WITHIN STRING. THE
C     NUMBER OF CHARACTERS SHIFTED IS NUMBER AND THEY ARE SHIFTED SO
C     THAT THE CHARACTER IN POSITION FROM IS MOVED TO POSITION TO.
C     CHARACTERS IN THE TO POSITION ARE OVERWRITTEN. BLANKS REPLACE
C     CHARACTERS IN THE FROM POSITION. SHIFTING MAY BE LEFT OR RIGHT,
C     AND THE FROM AND TO POSITIONS MAY OVERLAP.  CARE IS TAKEN NOT
C     TO ALTER OR USE ANY CHARACTERS BEYOND THE DEFINED LIMITS
C     OF THE STRING.

C## E N T R Y   P O I N T S: THE NATURAL ENTRY ZZSHFT.
C## S U B R O U T I N E S:   LEN  MIN  MAX (INTRINSIC)
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
                             INTEGER  N, SHIFT, INCR, I, IS, IE
                             INTEGER IBS, ETO, EFROM, K, SLEN
                             CHARACTER *1   CH

C## S A V E:                 NONE SELECTED.
C## E Q U I V A L E N C E S: NONE ARE DEFINED.
C## C O M M O N:             NONE IS DEFINED.
C## D A T A:                 NO DATA VALUES SET.

C##                                                  E X E C U T I O N
C##                                                  E X E C U T I O N

      SLEN  = LEN (STRING)
      N     = NUMBER - 1
      SHIFT = FROM - TO

      IF ( FROM .NE. TO ) THEN
         IF ( TO .LE. FROM ) THEN
            INCR = 1
            IS   = MIN( FROM+MAX(0,1-TO), SLEN+1 )
            IE   = MIN( FROM+N, SLEN )
            IBS  = MAX( IE-SHIFT+1, MAX(FROM,0) )
         ELSE
            INCR  = -1
            ETO   = TO + N
            EFROM = FROM + N
            IS    = MAX( EFROM - MAX(0,ETO-SLEN) , 0 )
            IE    = MAX(FROM , 0)
            IBS   = MIN( TO-1 , MIN(EFROM,SLEN) )
         ENDIF

         DO 1000 I=IS,IE,INCR
            K  = I - SHIFT
            CH = STRING(I:I)
            STRING(K:K) = CH
 1000    CONTINUE

         DO 2000 I=IBS,IE,INCR
            STRING(I:I) = BLANK
 2000    CONTINUE
      ENDIF
      GOTO 90000

C## E X I T
90000      RETURN

C## F O R M A T S:  NONE ARE DEFINED.
C##                 E N D         OF ZZSHFT.
                    END
