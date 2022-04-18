      SUBROUTINE ZZGETG ( UNIT,   NUM,    NAM,  SIZE, MEMBS,
     -                    NAMES,  GROUPS, ERRFLG,  *        )

C## A R G U M E N T S:
                       INTEGER           UNIT, NUM, SIZE, ERRFLG
                       INTEGER           MEMBS( * ), GROUPS( * )
                       CHARACTER * ( * ) NAM, NAMES

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
C>RCS $HEADER: GETG.F,V 1.10 91/11/20 10:52:55 BUCKLEY EXP $
C>RCS $LOG:     GETG.F,V $
C>RCS REVISION 1.10  91/11/20  10:52:55  BUCKLEY
C>RCS FINAL SUBMISSION TO TOMS
C>RCS
C>RCS REVISION 1.9  89/06/30  13:38:16  BUCKLEY
C>RCS PREPARING SUBMITTED VERSION OF MT
C>RCS
C>RCS REVISION 1.3.1.1  89/05/20  16:42:49  BUCKLEY
C>RCS TEMP. TEST OF MT BEFORE SUBMITTING
C>RCS
C>RCS REVISION 1.3  89/05/18  12:20:38  BUCKLEY
C>RCS FINAL TEST OF MT BEFORE SUBMITTING
C>RCS
C>RCS REVISION 1.2  89/05/15  14:47:52  BUCKLEY
C>RCS INITIAL INSTALLATION OF MT INTO RCS FORM.
C>RCS
C>RCS REVISION 1.1  89/01/17  16:50:10  BUCKLEY
C>RCS INITIAL REVISION
C>RCS
C
C## D E S C R I P T I O N:
C
C         THIS SUBROUTINE READS GROUP INFORMATION FROM THE
C         GIVEN FILE USING THE ARRAYS 'NAMES' AND 'GROUPS'.
C
C     PARAMETER DESCRIPTION:
C
C     ON ENTRY:
C         UNIT   - INPUT UNIT NUMBER OF DAUF
C         NUM    - GROUP NUMBER
C         NAM    - GROUP NAME
C
C     ON EXIT:
C         SIZE   - THE SIZE OF THE GROUP
C         MEMBS  - THE MEMBERS OF THE GROUP
C         NAMES  - AN STRING OF CHARACTERS THAT CONTAINS THE
C                    LIST OF GROUP NAMES
C         GROUPS - AN ARRAY CONTAINING THE POSITION OF THE
C                    GROUP INFORMATION IN THE DAUF FILE
C         ERRFLG - AN ERROR RETURN. IT IS SET TO
C                       0  IF NO ERROR OCCURRED
C                       1  IF THE GROUP NUMBER IS OUT OF RANGE
C                       2  IF THE GROUP IS UNDEFINED
C
C## E N T R Y   P O I N T S:
C
C          ZZGETG - THE NATURAL ENTRY POINT
C          ZZDEFG - DEFINES THE SIZE OF RECORDS IN THE DAUF FILE
C                     AND SOME OTHER CONSTANTS
C
C## S U B R O U T I N E S:
C
C          RD     ... STATEMENT FUNCTION TO CONVERT REALS TO INTEGERS
C          DBLE(REAL)... INTRINSIC
C          ZZRDIN ... READS INTEGERS FROM DAUF FILE
C          ZZSRCH ... SEARCHES A DICTIONARY
C          ZZERRM ... OUTPUTS ERROR MESSAGES
C
C## P A R A M E T E R S:


      LOGICAL     T,          F
      PARAMETER ( T = .TRUE., F = .FALSE. )

      CHARACTER*(*) TRUE,          QT,       FALSE,           QF
      PARAMETER (   TRUE = 'TRUE', QT = 'T', FALSE = 'FALSE', QF = 'F' )

      INTEGER     ITRUE,     IFALSE
      PARAMETER ( ITRUE = 1, IFALSE = 0 )

      DOUBLE PRECISION  RTRUE,        RFALSE
C!!!! REAL              RTRUE,        RFALSE
      PARAMETER      (  RTRUE = 1.D0, RFALSE = 0.D0 )

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

C----                      STATEMENT FUNCTION
      DOUBLE PRECISION  RD
C!!!! REAL              RD

C----                      NUMBER OF ELEMENTS PER LINE OF UNIT
      INTEGER   IPERLN, SIPRLN

C----                      GROUP NAME LENGTH

      INTEGER     PNAMLN,     FNAMLN,          GNAMLN
      PARAMETER ( PNAMLN = 8, FNAMLN = PNAMLN, GNAMLN = PNAMLN )

      INTEGER     TITLEN,      PDESCL
      PARAMETER ( TITLEN = 72, PDESCL = 72 )

C----                      MAXIMUM NUMBER OF GROUPS

      INTEGER     NOFNS,      DFPRBS,       MXGRPS,      MXGSZ
      PARAMETER ( NOFNS = 80, DFPRBS = 450, MXGRPS = 50, MXGSZ = 200 )

C----                      MISCELLANEOUS VARIABLES
      INTEGER   RECNO, I

C## S A V E:
            SAVE   IPERLN

C## E Q U I V A L E N C E S: NONE ARE DEFINED.
C## C O M M O N:             NONE IS DEFINED.
C## D A T A:                 NONE ARE SET.
C##                                                E X E C U T I O N
C##                                                E X E C U T I O N

C----                DEFINE STATEMENT FUNCTION
      RD(I) = DBLE(I)
C!!!! RD(I) = REAL(I)

C---- SET THE DEFAULT SIZE AND SEARCH START LOCATION.

      SIZE   = 0
      ERRFLG = 0

C---- CHECK NAME AND NUMBER

      IF ( NUM .LT. 1  .OR.  NUM .GT. MXGRPS ) THEN
C                                                   FIND GROUP NAME
         NUM = 0
         IF ( NAM .NE. BLANK ) THEN
            CALL ZZSRCH( NAM, GNAMLN, NAMES, MXGRPS,GNAMLN,NUM,F,T)
            IF ( NUM .EQ. 0 ) THEN
                ERRFLG = 2
                CALL ZZERRM ( RD(NUM), *91000, 'NS NO GROUP NAMED'
     -                      //NAM(1:GNAMLN) )
                GOTO 90000
            ENDIF
         ELSE
            ERRFLG = 1
            CALL ZZERRM (RD(NUM),*91000,'NS GROUP NUMBER OUT OF RANGE')
            GOTO 90000
         ENDIF

      ELSE IF ( GROUPS(NUM) .EQ. -1 ) THEN
         ERRFLG = 2
         CALL ZZERRM( RD(NUM), *91000, 'IS NO GROUP #' )
         GOTO 90000
      ENDIF

C----                       GET THE STARTING RECORD NUMBER.
      RECNO = GROUPS( NUM )

C----                       READ THE SIZE.

         CALL ZZRDIN( UNIT, SIZE, 1, IPERLN, RECNO )

C----                       READ MEMBS ARRAYS.

         CALL ZZRDIN( UNIT, MEMBS(1), SIZE, IPERLN, RECNO )
      GOTO 90000

C## E N T R Y  ZZFPAR:
                       ENTRY ZZDEFG ( SIPRLN )
      IPERLN = SIPRLN
      RETURN

C## E X I T
90000       RETURN
91000       RETURN1

C## F O R M A T S:  NONE ARE DEFINED.
C##                 E N D         OF ZZGETG.
                    END
