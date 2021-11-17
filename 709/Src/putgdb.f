      SUBROUTINE ZZPUTG ( UNIT,  NUM,    NAM,  SIZE, MEMBS,
     -                    NAMES, GROUPS, RECNO, ERRFLG,  *  )

C## A R G U M E N T S:
                      INTEGER           UNIT, NUM, SIZE, RECNO, ERRFLG
                      INTEGER           MEMBS( * ), GROUPS( * )
                      CHARACTER * ( * ) NAM, NAMES

C## S T A T U S:
C               SINGLE/DOUBLE CONVERSION:        NOT REQUIRED.
C               SYSTEM  DEPENDENCE:                      NONE.
C
C>RCS $HEADER: PUTG.F,V 1.10 91/11/20 10:53:03 BUCKLEY EXP $
C>RCS $LOG:     PUTG.F,V $
C>RCS REVISION 1.10  91/11/20  10:53:03  BUCKLEY
C>RCS FINAL SUBMISSION TO TOMS
C>RCS
C>RCS REVISION 1.9  89/06/30  13:39:47  BUCKLEY
C>RCS PREPARING SUBMITTED VERSION OF MT
C>RCS
C>RCS REVISION 1.3.1.1  89/05/20  16:42:54  BUCKLEY
C>RCS TEMP. TEST OF MT BEFORE SUBMITTING
C>RCS
C>RCS REVISION 1.3  89/05/18  12:20:47  BUCKLEY
C>RCS FINAL TEST OF MT BEFORE SUBMITTING
C>RCS
C>RCS REVISION 1.2  89/05/15  14:48:08  BUCKLEY
C>RCS INITIAL INSTALLATION OF MT INTO RCS FORM.
C>RCS
C>RCS REVISION 1.1  89/01/17  16:50:12  BUCKLEY
C>RCS INITIAL REVISION
C>RCS
C
C## D E S C R I P T I O N:
C
C     THIS SUBROUTINE WRITES GROUP INFORMATION TO THE GIVEN FILE
C     AND TO THE ARRAYS 'NAMES' AND 'GROUPS'.
C
C     PARAMETER DESCRIPTION:
C
C     ON ENTRY:
C         UNIT   - OUTPUT UNIT NUMBER OF DAUF
C         NUM    - GROUP NUMBER
C         NAM    - GROUP NAME
C         SIZE   - THE SIZE OF THE GROUP
C         MEMBS  - THE MEMBERS OF THE GROUP
C         NAMES  - AN ARRAY OF CHARACTERS THAT CONTAINS THE
C                    LIST OF GROUP NAMES
C         GROUPS - AN ARRAY CONTAINING THE POSITION OF THE
C                    GROUP INFORMATION IN THE DAUF FILE
C         RECNO  - THE NEXT FREE RECORD IN DAUF
C
C     ON EXIT:
C         ERRFLG - AN ERROR RETURN. IT IS SET TO
C                       0  IF NO ERRORS OCCURRED
C                       1  IF ILLEGAL NAME OR NUMBER IS GIVEN
C                       2  IF GROUP SET IS FULL
C                       3  IF THERE IS A BAD PROBLEM IN THE GROUP
C                       4  IF THE GROUP IS EMPTY
C
C## E N T R Y   P O I N T S:
C
C          ZZPUTG - THE NATURAL ENTRY POINT
C          ZZDEFP - DEFINES THE SIZE OF RECORDS IN THE DAUF FILE
C                     AND SOME OTHER CONSTANTS
C
C## S U B R O U T I N E S:
C
C          DBLE(REAL)... INTRINSIC
C
C          ZZWRIN ... WRITES INTEGER TO DAUF FILE
C          ZZERRM ... PRINTS ERROR MESSAGES
C          ZZSRCH ... SEARCHES A DICTIONARY
C
C          RD  ... STATEMENT FUNCTION
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

C---- ERROR MESSAGES

      CHARACTER *(*) ERR1, ERR2, ERR3, ERR4, ERR5, ERR6
      PARAMETER ( ERR1 = 'GROUP NAME USED IN GROUP #',
     -            ERR2 = 'GROUP SET FULL',
     -            ERR3 = 'GROUP NUMBER TOO LARGE',
     -            ERR4 = 'ALREADY USED GROUP #',
     -            ERR5 = 'BAD PROBLEM IN GROUP POSITION #',
     -            ERR6 = 'EMPTY GROUP' )

C## L O C A L   D E C L:

C----                       INTEGER TO REAL FUNCTION
      DOUBLE PRECISION  RD
C!!!! REAL              RD
C----                        NUMBER OF ELEMENTS PER LINE OF UNIT
      INTEGER   IPERLN
      INTEGER   SIPRLN
C----                        GROUP NAME LENGTH

      INTEGER     PNAMLN,     FNAMLN,          GNAMLN
      PARAMETER ( PNAMLN = 8, FNAMLN = PNAMLN, GNAMLN = PNAMLN )

      INTEGER     TITLEN,      PDESCL
      PARAMETER ( TITLEN = 72, PDESCL = 72 )
C----                        MAXIMUM NUMBER OF GROUPS AND PROBLEMS

      INTEGER     NOFNS,      DFPRBS,       MXGRPS,      MXGSZ
      PARAMETER ( NOFNS = 80, DFPRBS = 450, MXGRPS = 50, MXGSZ = 200 )

C----                         MISCELLANEOUS VARIABLES
      INTEGER   I, J

C## S A V E:
             SAVE   IPERLN

C## E Q U I V A L E N C E S: NONE ARE DEFINED.
C## C O M M O N:             NONE IS DEFINED.
C## D A T A:                 NONE ARE SET.
C##                                                E X E C U T I O N
C##                                                E X E C U T I O N
C----                 DEFINE FUNCTION STATEMENT
      RD(I) = DBLE(I)
C!!!! RD(I) = REAL(I)
C----

      ERRFLG = 0
C----                            CHECK THE NAME AND NUMBER
      IF ( NAM .NE. BLANK ) THEN
C                                SEE IF NAME HAS BEEN USED
         CALL ZZSRCH(NAM,GNAMLN,NAMES,MXGRPS,GNAMLN,I,F,T)

         IF ( I .NE. 0 ) THEN
             ERRFLG = 1
             CALL ZZERRM( RD(I), *91000, 'IS'//ERR1 )
             GOTO 90000
         ENDIF
      ENDIF

      IF ( NUM .LE. 0 ) THEN
C                             FIND FIRST UNUSED GROUP NUMBER
         DO 100 I = 1, MXGRPS
            IF ( GROUPS(I) .EQ. -1 ) THEN
                NUM = I
                GOTO 200
            ENDIF
 100     CONTINUE

         ERRFLG = 2
         CALL ZZERRM ( RD(I), *91000, 'NS'//ERR2 )
         GOTO 90000

      ELSE IF ( NUM .GE. MXGRPS ) THEN
         ERRFLG = 1
         CALL ZZERRM ( RD(NUM), *91000, 'NS'//ERR3 )
         GOTO 90000
      ELSE IF ( GROUPS(NUM) .NE. -1 ) THEN
         ERRFLG = 1
         CALL ZZERRM ( RD(NUM), *91000, 'IS'//ERR4 )
         GOTO 90000
      ENDIF
 200  CONTINUE

C----        CHECK THE PROBLEM NUMBERS
      J = 0
      DO 300  I = 1, SIZE
         IF ( MEMBS(I) .LE. 0 .OR. MEMBS(I) .GT. DFPRBS ) THEN
            ERRFLG = 3
            CALL ZZERRM ( RD(I), *91000, 'NS'//ERR5 )
            MEMBS(I) = -1
            J = J+1
         ENDIF
 300  CONTINUE

      IF ( SIZE-J .EQ. 0 ) THEN
         ERRFLG = 4
         CALL ZZERRM( RD(J), *91000, 'NS'//ERR6 )
         GOTO 90000
      ENDIF

C----                       UPDATE THE GROUPS ARRAY
      GROUPS( NUM ) = RECNO

C----                       UPDATE NAME ARRAY
      I = ( NUM-1 ) * GNAMLN + 1
      J = I + GNAMLN - 1
      NAMES( I : J ) = NAM

C---- WRITE THE SIZE

      CALL ZZWRIN( UNIT, SIZE, 1, IPERLN, RECNO )

C---- WRITE MEMBS ARRAY

      CALL ZZWRIN( UNIT, MEMBS(1), SIZE, IPERLN, RECNO )
      GOTO 90000

C## E N T R Y  ZZGDEF:
                      ENTRY ZZDEFP ( SIPRLN )
      IPERLN = SIPRLN
      RETURN

C## E X I T
90000        RETURN
91000        RETURN1

C## F O R M A T S:  NONE ARE DEFINED.
C##                 E N D         OF ZZPUTG.
                    END
