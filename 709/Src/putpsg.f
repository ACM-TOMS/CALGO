      SUBROUTINE ZZPUTP ( UNIT,   PNUM,  FNUM,  PNAME, FNAME, DESC,
     -                    SOLNS,  INTS,  ORDER, LOOPX, LOOPC, X0,
     -                    NAMES,  PRECNO, RECNO, ERRFLG,  *        )

C## A R G U M E N T S:

      INTEGER           UNIT, PNUM, FNUM, RECNO, ERRFLG
      INTEGER           INTS( * ), PRECNO( 3,* ), ORDER( * ), LOOPX( * )

      REAL              LOOPC( * ), X0( * ), SOLNS(*)
C!!!! DOUBLE PRECISION  LOOPC( * ), X0( * ), SOLNS(*)

      CHARACTER * ( * ) PNAME, FNAME, DESC, NAMES

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
C>RCS $HEADER: PUTP.F,V 1.10 91/11/20 10:53:04 BUCKLEY EXP $
C>RCS $LOG:     PUTP.F,V $
C>RCS REVISION 1.10  91/11/20  10:53:04  BUCKLEY
C>RCS FINAL SUBMISSION TO TOMS
C>RCS
C>RCS REVISION 1.9  89/06/30  13:39:48  BUCKLEY
C>RCS PREPARING SUBMITTED VERSION OF MT
C>RCS
C>RCS REVISION 1.3.1.1  89/05/20  16:42:55  BUCKLEY
C>RCS TEMP. TEST OF MT BEFORE SUBMITTING
C>RCS
C>RCS REVISION 1.3  89/05/18  17:08:26  BUCKLEY
C>RCS FINAL TEST OF MT BEFORE SUBMITTING
C>RCS
C>RCS
C>RCS REVISION 1.2  89/05/15  14:48:10  BUCKLEY
C>RCS INITIAL INSTALLATION OF MT INTO RCS FORM.
C>RCS
C>RCS REVISION 1.1  89/01/17  16:50:12  BUCKLEY
C>RCS INITIAL REVISION
C>RCS
C
C## D E S C R I P T I O N:
C
C     THIS SUBROUTINE WRITES PROBLEM INFORMATION TO THE GIVEN FILE
C     AND TO THE ARRAYS 'NAMES' AND 'PRECNO'.
C
C     PARAMETER DESCRIPTION:
C
C     ON ENTRY:
C         UNIT   - OUTPUT UNIT NUMBER FOR DAUF
C         PNAME  - PROBLEM NAME
C         FNAME  - FUNCTION NAME
C         PCHRS  - STRING CONTAINING CHARACTERS TO BE WRITTEN
C                    TO UNIT
C         INTS   - ARRAY CONTAINING THE INTEGERS TO BE WRITTEN
C                    TO UNIT
C         RELS   - ARRAY CONTAINING THE DOUBLE PRECISION OR REAL
C                    ( AS THE CASE MAY BE ) NUMBERS TO BE WRITTEN
C                    TO UNIT
C         NAMES  - AN ARRAY OF CHARACTERS THAT CONTAINS THE
C                    PROBLEM NAME AND THE FUNCTION NAME FOR EACH
C                    PROBLEM
C         PRECNO - AN INTEGER ARRAY THAT CONTAINS THE DIMENSION,
C                    THE FUNCTION NUMBER, AND THE RECORD NUMBER
C                    ( IN THE DAUF FILE ) FOR EACH ROBLEM
C         RECNO  - THE NEXT FREE RECORD IN UNIT DAUF
C
C     ON EXIT:
C         ERRFLG - AN ERROR RETURN. IT IS SET TO
C                       0  IF NO ERRORS OCCURRED
C                       1  IF ILLEGAL NAME OR NUMBER IS GIVEN
C                       2  IF PROBLEM SET IS FULL
C                       3  IF NO FUNCTION IS GIVEN
C                       4  IF AN ILLEGAL DIMENSION IS GIVEN
C
C## E N T R Y   P O I N T S:
C
C          ZZPUTP - THE NATURAL ENTRY POINT
C          ZZPDEF - DEFINES THE SIZE OF RECORDS IN THE DAUF FILE
C                     AND SOME OTHER CONSTANTS
C
C## S U B R O U T I N E S:
C
C          NINT, DBLE(REAL) ... INTRINSICS
C          ZZERRM   PRINTS ERROR MESSAGES
C          ZZLPCK   OBTAINS LOOP DATA
C          ZZSRCH   SEARCHES A DICTIONARY
C          ZZWRCH   WRITES CHARACTERS TO THE DAUF FILE
C          ZZWRIN   WRITES INTEGERS TO THE DAUF FILE
C          ZZWRRL   WRITES REALS TO THE DAUF FILE
C          RD ...   A STATEMENT FUNCTION
C
C## P A R A M E T E R S:

      LOGICAL     T,          F
      PARAMETER ( T = .TRUE., F = .FALSE. )

      CHARACTER*(*) TRUE,          QT,       FALSE,           QF
      PARAMETER (   TRUE = 'TRUE', QT = 'T', FALSE = 'FALSE', QF = 'F' )

      INTEGER     ITRUE,     IFALSE
      PARAMETER ( ITRUE = 1, IFALSE = 0 )

      REAL              RTRUE,        RFALSE
C!!!! DOUBLE PRECISION  RTRUE,        RFALSE
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

      CHARACTER *(*) UNDEFN
      PARAMETER (    UNDEFN = 'UNDEFN' )

C---- POSITION FLAGS FOR ENTRIES IN THE PRECNO ARRAY

C   DEFINITIONS OF THE ROWS IN THE ARRAY PRECNO. ( PRECNO HOLDS THE
C   RECORD NUMBER IN THE DAUF FILE, THE MINIMUM DIMENSION, AND THE
C   FUNCTION NUMBER OF EACH PROBLEM. )

      INTEGER     RECN,     DIMN,     FNO1
      PARAMETER ( RECN = 1, DIMN = 2, FNO1 = 3 )

C---- POSITION FLAGS FOR ENTRIES IN THE INTS ARRAY


      INTEGER     PPTMAX,     PPTIX0,     PPTDES
      PARAMETER ( PPTMAX = 1, PPTIX0 = 2, PPTDES = 3 )

      INTEGER     PPTORD,     PPTLPX,     PPTLPC,     PPTSOL
      PARAMETER ( PPTORD = 4, PPTLPX = 5, PPTLPC = 6, PPTSOL = 7 )

      INTEGER     NPNTS
      PARAMETER ( NPNTS  = 7 )

C---- INFORMATIVE MESSAGES.

      CHARACTER *(*)  ERR1, ERR2, ERR3, ERR4, ERR5, ERR6, MESS
      PARAMETER ( ERR1 = 'PROBLEM NUMBER TOO LARGE, #',
     -            ERR2 = 'PROBLEM SET FULL',
     -            ERR3 = 'PROBLEM NAME USED IN PROBLEM #',
     -            ERR4 = 'ALREADY DEFINED PROBLEM #' ,
     -            ERR5 = 'NO FUNCTION SPECIFIED IN PROBLEM DEFINITION',
     -            ERR6 = 'ILLEGAL DIMENSION',
     -            MESS = 'PROBLEM ASSIGNED PROBLEM #'  )

C## L O C A L   D E C L:

C----                    INTEGER TO REAL FUNCTION
      REAL              RD
C!!!! DOUBLE PRECISION  RD

C---- DECLARE THE NUMBER OF ELEMENTS PER LINE OF UNIT

      INTEGER   CPERLN,  IPERLN,  RPERLN
      INTEGER   SCPRLN,  SIPRLN,  SRPRLN

C---- DECLARE THE PROBLEM AND FUNCTION NAME LENGTHS


      INTEGER     PNAMLN,     FNAMLN,          GNAMLN
      PARAMETER ( PNAMLN = 8, FNAMLN = PNAMLN, GNAMLN = PNAMLN )

      INTEGER     TITLEN,      PDESCL
      PARAMETER ( TITLEN = 72, PDESCL = 72 )

C---- DECLARE MAXIMUM NUMBER OF PROBLEMS


      INTEGER     NOFNS,      DFPRBS,       MXGRPS,      MXGSZ
      PARAMETER ( NOFNS = 80, DFPRBS = 450, MXGRPS = 50, MXGSZ = 200 )

C---- POINTER TO DIMENSION DATA

      INTEGER  PTNDIM, SPTDIM

C---- MISC VARIABLES

      INTEGER  DIM, I, J, K

      LOGICAL  INFORM

      REAL              VAL
C!!!! DOUBLE PRECISION  VAL

C## S A V E:
            SAVE   CPERLN, IPERLN, RPERLN, PTNDIM

C## E Q U I V A L E N C E S: NONE ARE DEFINED.
C## C O M M O N:             NONE IS DEFINED.
C## D A T A:
             DATA   INFORM /F/
C##                                                E X E C U T I O N
C##                                                E X E C U T I O N
C-----                DEFINE FUNCTION STATEMENT.
      RD(I) = REAL(I)
C!!!! RD(I) = DBLE(I)

      ERRFLG = 0
C---- CHECK PROBLEM NAME AND NUMBER

      IF ( PNUM .GT. DFPRBS ) THEN
         ERRFLG = 1
         CALL ZZERRM ( RD(PNUM), *91000, 'IS'//ERR1 )
         GOTO 90000
      ELSE IF ( PNUM .EQ. -1 ) THEN
C                                 FIND NEXT FREE PROBLEM NUMBER
         DO 100 I = 1, DFPRBS
           IF ( PRECNO( RECN, I ) .EQ. -1 ) THEN
              PNUM   = I
              INFORM = T
              GOTO 200
           ENDIF
 100     CONTINUE

         ERRFLG = 2
         CALL ZZERRM ( RD(I), *91000, 'NS'//ERR2 )
         GOTO 90000
      ENDIF
 200  CONTINUE

      IF ( PNAME .NE. UNDEFN ) THEN
         K = 1
         CALL ZZSRCH ( PNAME, PNAMLN, NAMES, DFPRBS,
     -                 PNAMLN+FNAMLN, K, F ,T  )
         IF ( K .NE. 0 ) THEN
            ERRFLG = 1
            CALL ZZERRM( RD(K), *91000, 'IS'//ERR3 )

         ENDIF
      ENDIF

      IF ( PRECNO(RECN,PNUM) .NE. -1 ) THEN
         ERRFLG = 1
         CALL ZZERRM ( RD(PNUM), *91000, 'IS'//ERR4 )
         GOTO 90000
      ENDIF

C---- CHECK THE FUNCTION NAME AND NUMBER

      IF ( FNUM .EQ. -1 .AND. FNAME .EQ. BLANK ) THEN
         ERRFLG = 3
         CALL ZZERRM ( RD(I), *91000, 'NS'//ERR5 )
         GOTO 90000
      ENDIF

C---- FIND THE FIRST DIMENSION

         CALL ZZLPCK( LOOPC, LOOPX, ORDER, PTNDIM, VAL )
         DIM = NINT(VAL)
         IF ( DIM .LT. 2 ) THEN
            ERRFLG = 4
            CALL ZZERRM ( RD(DIM), *91000, 'NT'//ERR6 )
            GOTO 90000
         ENDIF

C---- UPDATE THE PRECNO ARRAY

      PRECNO( FNO1, PNUM ) = FNUM
      PRECNO( RECN, PNUM ) = RECNO
      PRECNO( DIMN, PNUM ) = DIM

C---- UPDATE NAME ARRAY

         I = ( PNUM-1 ) * ( 2*PNAMLN ) + 1
         J = I + PNAMLN - 1
         K = J + PNAMLN
         NAMES( I : J ) = PNAME
         NAMES( J+1 : K ) = FNAME

C---- WRITE INTS ARRAY

         CALL ZZWRIN( UNIT, INTS, NPNTS, IPERLN, RECNO )

C---- WRITE CHARACTERS

         CALL ZZWRCH( UNIT, DESC,  INTS(PPTDES), CPERLN, RECNO )

C---- WRITE INTEGER ARRAYS

         CALL ZZWRIN( UNIT, ORDER, INTS(PPTORD), IPERLN, RECNO )
         CALL ZZWRIN( UNIT, LOOPX, INTS(PPTLPX), IPERLN, RECNO )

C---- WRITE REAL ARRAY(S)

         CALL ZZWRRL( UNIT, LOOPC, INTS(PPTLPC), RPERLN, RECNO )
         IF ( INTS(PPTIX0) .GT. 0 ) THEN
            CALL ZZWRRL( UNIT, X0, INTS(PPTIX0), RPERLN, RECNO )
         ENDIF

         IF ( INTS(PPTSOL) .GT. 0 ) THEN
            CALL ZZWRRL( UNIT, SOLNS, INTS(PPTSOL), RPERLN, RECNO )
         ENDIF

C---- INFORM USER THAT DEFAULT PROBLEM NUMBER WAS USED

        IF ( INFORM ) THEN
           CALL ZZERRM( RD(PNUM), *91000, 'IT'//MESS )
        ENDIF

C---- GOTO THE EXIT POINT
        GOTO 90000

C## E N T R Y  ZZPDEF:
                       ENTRY ZZPDEF ( SCPRLN, SIPRLN, SRPRLN, SPTDIM )
      CPERLN = SCPRLN
      IPERLN = SIPRLN
      RPERLN = SRPRLN
      PTNDIM = SPTDIM
      RETURN

C## E X I T
90000       RETURN
91000       RETURN1

C## F O R M A T S:  NONE ARE DEFINED.
C
C                       THERE ARE NO FORMATS USED.
C
C##                 E N D         OF STATEMENT..
                    END
