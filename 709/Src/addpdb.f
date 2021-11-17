      SUBROUTINE ZZADDP (    ADD,   NAMES, NUMBRS, NUMBER,   LIST,
     -                    LISTNO,  LISTMX, PRECNO, PFNAMS, PROBNO,
     -                    GROUPS,  GNAMES,  NGRPS,  MEMBS,    *   )

C## A R G U M E N T S:
      INTEGER          NUMBER, LISTMX, LISTNO, PROBNO, NGRPS
      INTEGER          LIST(*), PRECNO(3,PROBNO), GROUPS(*), MEMBS(*)
      LOGICAL          ADD

      CHARACTER*(*)    NAMES(*),  PFNAMS, GNAMES

      DOUBLE PRECISION NUMBRS(*)
C!!!! REAL             NUMBRS(*)

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
C    SYSTEM  DEPENDENCE:               NONE
C
C>RCS $HEADER: ADDP.F,V 1.10 91/11/20 10:52:35 BUCKLEY EXP $
C>RCS $LOG:     ADDP.F,V $
C>RCS REVISION 1.10  91/11/20  10:52:35  BUCKLEY
C>RCS FINAL SUBMISSION TO TOMS
C>RCS
C>RCS REVISION 1.9  89/06/30  13:38:10  BUCKLEY
C>RCS PREPARING SUBMITTED VERSION OF MT
C>RCS
C>RCS REVISION 1.3.1.1  89/05/20  16:42:28  BUCKLEY
C>RCS TEMP. TEST OF MT BEFORE SUBMITTING
C>RCS
C>RCS REVISION 1.3  89/05/18  12:20:12  BUCKLEY
C>RCS FINAL TEST OF MT BEFORE SUBMITTING
C>RCS
C>RCS REVISION 1.2  89/05/15  14:46:53  BUCKLEY
C>RCS INITIAL INSTALLATION OF MT INTO RCS FORM.
C>RCS
C>RCS REVISION 1.1  89/01/17  16:50:05  BUCKLEY
C>RCS INITIAL REVISION
C>RCS
C
C## D E S C R I P T I O N:
C
C     THIS ROUTINE IS RESPONSIBLE FOR ADDING A PROBLEM, OR A COLLECTION
C     OF PROBLEMS, TO THE CURRENT LIST OF PROBLEMS TO BE EXECUTED. THE
C     PROBLEMS MAY BE DROPPED FROM THE LIST AS WELL.  IT IS POSSIBLE
C     TO SPECIFY A SINGLE PROBLEM, OR ALL PROBLEMS USING A PARTICULAR
C     FUNCTION, OR ALL PROBLEMS IN A PREDEFINED GROUP.  THE ARGUMENTS
C     HAVE THE FOLLOWING MEANINGS.
C
C        ADD     A FLAG: IF TRUE ADD THE PROBLEMS; OTHERWISE DROP THEM.
C
C        NAMES   THE ARRAY OF NAMES TO ADD OR DROP.
C        NUMBRS  THE ARRAY OF NUMBERS TO ADD OR DROP.
C        NUMBER  THE NUMBER OF ENTRIES IN NAMES AND NUMBRS.
C
C          THE FORMAT OF THESE TWO ARRAYS IS AS FOLLOWS :
C
C               THE I-TH ENTRY OF NAMES CORRESPONDS TO THE I-TH ENTRY
C               OF NUMBRS. IF THE NAME IS BLANK, THEN IT IS ASSUMED
C               THAT NUMBR CONTAINS THE NUMBER OF A PROBLEM TO ADD OR
C               DROP. IF IT IS NONBLANK AND CONTAINS NO '/' (I.E. THE
C               'RADIX' CHARACTER) IT IS ASSUMED THAT IT CONTAINS THE
C               NAME OF AN ENTITY TO ADD OR DROP. (THIS IS ALSO TRUE IF
C               THE '/' IS THE LAST NONBLANK CHARACTER.) IF IT CONTAINS
C               A '/' IN THE FIRST POSITION, THEN THE ARRAY NUMBR MUST
C               CONTAIN THE NUMBER OF THE ENTITY TO BE DROPPED OR
C               ADDED, AND THE CHARACTER FOLLOWING THE '/' INDICATES
C               THE TYPE OF THE ENTITY, I.E. 'P' FOR A PROBLEM, 'G'
C               FOR A GROUP, OR 'F' FOR A FUNCTION. IN THE REMAINING
C               CASE IT IS TAKEN THAT A NAME IS GIVEN FOR THE ENTITY
C               AND THAT THE ENTITY IDENTIFIER FOLLOWS THE SLASH.
C
C        LIST    THE CURRENT LIST OF PROBLEMS.
C        LISTNO  THE NUMBER OF ENTRIES CURRENTLY IN LIST.
C        LISTMX  THE NUMBER OF ELEMENTS ALLOWED IN LIST.
C
C        PRECNO  AN ARRAY CONTAINING INFOMATION ABOUT EACH PROBLEM IN
C                PROLOG. FOR EACH PROBLEM, THERE ARE 3 ENTRIES. THE
C                FIRST GIVES THE RECORD NUMBER IN DAUF WHERE THE PROBLEM
C                MAY BE FOUND. THE SECOND GIVES THE DIMENSION OF THE
C                PROBLEM AND THE THIRD GIVES THE NUMBER OF THE FUNCTION
C                USED BY THIS PROBLEM.
C        PFNAMS  AN ARRAY GIVING THE NAME FOR EACH PROBLEM, AS WELL AS
C                THE NAME OF THE FUNCTION USED BY THE PROBLEM.
C        PROBNO  THE TOTAL NUMBER OF PROBLEMS DEFINED.
C
C        GROUPS  AN ARRAY WHICH GIVES THE RECORD NUMBER IN DAUF WHERE
C                EACH GROUP MAY BE FOUND.
C        GNAMES  A CHARACTER STRING GIVING THE NAME OF EACH GROUP.
C        NGRPS   THE NUMBER OF GROUPS POSSIBLE.
C        MEMBS   A DUMMY ARRAY WHOSE SIZE IS THE MAXIMUM NUMBER OF
C                PROBLEMS PER GROUP
C
C          *     AN ERROR RETURN.
C
C## E N T R Y   P O I N T S: ZZADDP   THE NATURAL ENTRY POINT.
C                            ZZASET   DEFINE SPECIAL RADIX CHARS.
C
C## S U B R O U T I N E S:
C
C     INDEX, REAL(DBLE), NINT    ...GENERIC
C
C     ZZLENG  NON-BLANK LENGTH OF A STRING
C     ZZERRM  ERROR MESSAGES
C     ZZGETG  GETS GROUP OF PROBLEMS FROM DAUF
C     ZZSRCH  SEARCHES A DICTIONARY
C
C     RD...   STATEMENT FUNCTION
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

C-----I/O UNIT NUMBERS

      INTEGER     PREPRC,     DAUF,       INPTUN
      PARAMETER ( PREPRC = 1, DAUF =   2, INPTUN = 3 )

      INTEGER     TEMPUN,     STDIN,     TRMOUT
      PARAMETER ( TEMPUN = 4, STDIN = 5, TRMOUT = 6 )

      INTEGER     WRITUN,     TRACUN,     SUMMUN
      PARAMETER ( WRITUN = 7, TRACUN = 8, SUMMUN = 9 )

      INTEGER     COPYUN,     TEXUN
      PARAMETER ( COPYUN =10, TEXUN =11 )

C-----STRING LENGTHS

      INTEGER     PNAMLN,     FNAMLN,          GNAMLN
      PARAMETER ( PNAMLN = 8, FNAMLN = PNAMLN, GNAMLN = PNAMLN )

      INTEGER     TITLEN,      PDESCL
      PARAMETER ( TITLEN = 72, PDESCL = 72 )

C-----CLASSES AND MODES OF PROBLEM REQUESTS.

      INTEGER       PROBS,      FUNCS,      GRPS,      ANY
      PARAMETER  (  PROBS = 1,  FUNCS = 2,  GRPS = 3,  ANY = 4  )

      INTEGER     BYNAME,      BYNUMB
      PARAMETER ( BYNAME = 1,  BYNUMB = 2  )

C-----ROW NUMBER FOR RECORD NUMBER, DIMENSION AND FUNCTION NUMBER.

C   DEFINITIONS OF THE ROWS IN THE ARRAY PRECNO. ( PRECNO HOLDS THE
C   RECORD NUMBER IN THE DAUF FILE, THE MINIMUM DIMENSION, AND THE
C   FUNCTION NUMBER OF EACH PROBLEM. )

      INTEGER     RECN,     DIMN,     FNO1
      PARAMETER ( RECN = 1, DIMN = 2, FNO1 = 3 )

C-----DEFAULT VALUES FOR SPECIAL RADIX CHARACTERS.

      CHARACTER *(*) DRAD,      DRADP,      DRADF,      DRADG
      PARAMETER (    DRAD= '/', DRADP= 'P', DRADF= 'F', DRADG= 'G' )

C-----ERROR MESSAGES.

      CHARACTER *(PNAMLN)  NAME

      CHARACTER *(*)  ERR1, ERR2, ERR3, ERR4, ERR5, ERR6

      PARAMETER (  ERR1 = ' CAN''T FIND PROBLEM ' ,
     -             ERR2 = ' CAN''T FIND FUNCTION ' ,
     -             ERR3 = ' CAN''T FIND GROUP ' ,
     -             ERR4 = ' CAN''T FIND PROBLEM, FUNCTION OR GROUP ',
     -             ERR5 = ' PROBLEM CAPACITY EXCEEDED',
     -             ERR6 = ' LIST ALREADY CONTAINS PROBLEM')

      CHARACTER *(*)   NMD,             NUM
      PARAMETER (      NMD = ' NAMED ', NUM = ' #' )

C## L O C A L   D E C L:

C-----CHARACTERS WHICH SIGNIFY PROBLEMS AND GROUPS.

      CHARACTER*1  RADIX, RADPRB, RADFUN, RADGRP
      CHARACTER*1 SRADIX, SRADP,  SRADF,  SRADG

C-----FUNCTION DECLARATION.

      INTEGER  ZZLENG

C-----LOCAL VARIABLES.

      INTEGER  CLASS,  MODE,  PNUMB, IMARK
      INTEGER  COUNT,  LEFT,  LOC,   SIZE
      INTEGER  I,      J,     K
      INTEGER  ERRFLG

      LOGICAL        DROP,   FOUND,   TRY

      CHARACTER * 1  CH

      DOUBLE PRECISION  RD
C!!!! REAL              RD

C## S A V E:
       SAVE  RADIX, RADPRB, RADFUN, RADGRP

C## E Q U I V A L E N C E S: NONE ARE DEFINED.
C## C O M M O N:             NONE IS DEFINED.
C## D A T A:

      DATA  RADIX/DRAD/,   RADPRB/DRADP/, RADGRP/DRADG/, RADFUN/DRADF/

C##                                                E X E C U T I O N
C##                                                E X E C U T I O N
C-----DEFINE FUNCTION STATEMENT.
      RD(I) = DBLE(I)
C!!!! RD(I) = REAL(I)
C-----NOW EXECUTE.

      DROP   = .NOT. ADD
      LEFT   =  LISTNO
      COUNT  = 0

C-----LOOP THOUGH EACH PROBLEM.

 1000 FOUND = F

      IF ( COUNT .LT. NUMBER ) THEN
C                                   CONTINUE THE LOOP.
         COUNT = COUNT + 1
      ELSE
C       109     FINISHED, SO EXIT, BUT FIRST...
         IF ( DROP .AND.  LEFT .NE. LISTNO ) THEN

C           COMPRESS OUT THE -1 ENTRIES.
            I = 0
            J = 0

 1100       J = J+1
            IF ( J .LE. LISTNO ) THEN
               IF ( LIST(J) .EQ. -1 ) THEN
                  GOTO 1100
               ENDIF
               I       = I + 1
               LIST(I) = LIST(J)
               GOTO 1100
            ENDIF
            LISTNO = LEFT
         ENDIF
         GOTO 90000

      ENDIF

C-----DETERMINE CLASS AND MODE OF PROBLEM STATEMENT.

      IF ( NAMES (COUNT) .EQ. BLANK ) THEN
         CLASS = PROBS
         MODE  = BYNUMB
         PNUMB = NINT ( NUMBRS(COUNT) )
      ELSE
         IMARK = INDEX ( NAMES(COUNT), RADIX )
         IF ( IMARK .EQ. 0  .OR. IMARK .EQ. ZZLENG(NAMES(COUNT)) ) THEN
            CLASS = ANY
            MODE  = BYNAME
         ELSE
            IF ( IMARK .EQ. 1 ) THEN
               MODE = BYNUMB
            ELSE
               MODE = BYNAME
            ENDIF

            CH   = NAMES(COUNT) ( IMARK+1:IMARK+1 )
            CALL ZZCASE( CH, CTOUPP )

            IF      ( CH .EQ. RADPRB ) THEN
               CLASS = PROBS
            ELSE IF ( CH .EQ. RADFUN ) THEN
               CLASS = FUNCS
            ELSE IF ( CH .EQ. RADGRP ) THEN
               CLASS = GRPS
            ELSE
               CLASS = ANY
            ENDIF
         ENDIF

         IF ( MODE .EQ. BYNAME ) THEN

            IF ( IMARK .EQ. 0 ) THEN
               NAME = NAMES(COUNT)
            ELSE
               NAME = NAMES(COUNT) (1:IMARK-1)
            ENDIF

            CALL ZZCASE(NAME,CTOUPP)
            PNUMB = 0
         ELSE
            PNUMB = NINT( NUMBRS(COUNT) )
            NAME = BLANK
         ENDIF
      ENDIF

C-----FIRST TRY TO FIND THE NAME OR NUMBER GIVEN AS A PROBLEM.

      IF ( CLASS .EQ. PROBS .OR. CLASS .EQ. ANY ) THEN
         IF ( MODE .EQ. BYNUMB ) THEN
            IF ( PNUMB .GE. 1  .AND. PNUMB  .LE. PROBNO ) THEN
               IF ( PRECNO(RECN,PNUMB) .NE. -1 ) THEN

                  IF ( DROP ) THEN
                     DO 1140 I = 1,LISTNO
                        IF ( LIST(I) .EQ. PNUMB ) THEN
                           FOUND   = T
                           LIST(I) = -1
                           LEFT    = LEFT -1
                           GOTO 1145
                        ENDIF
 1140                CONTINUE
 1145                CONTINUE
                  ELSEIF ( LISTNO .LT. LISTMX ) THEN
                     FOUND = T
                     DO 1147 K = 1, LISTNO
                        IF ( PNUMB .EQ. LIST(K) ) THEN
                           CALL ZZERRM( RD(PNUMB), *91000,
     -                                 'IT'//ERR6//NUM    )
                           GOTO 1148
                        ENDIF
 1147                CONTINUE
                     LISTNO = LISTNO + 1
                     LIST( LISTNO ) = PNUMB
 1148                CONTINUE
                  ELSE
                     GOTO 80000
                  ENDIF
               ENDIF
            ENDIF
         ELSE
C           ---NAMED PROBLEM...

            DO 1160 PNUMB = 1, PROBNO

               IF (       PRECNO(RECN,PNUMB) .GT. 0 ) THEN
                LOC  = PNUMB * PNAMLN * 2 - PNAMLN
                IF ( PFNAMS(LOC-PNAMLN+1:LOC) .EQ. NAME ) THEN

                  IF ( DROP ) THEN
                     DO 1150 I = 1,LISTNO
                        IF ( LIST(I) .EQ. PNUMB ) THEN
                           LIST(I) = -1
                           FOUND   = T
                           LEFT    = LEFT -1
                           GOTO 1180
                        ENDIF
1150                 CONTINUE

C                 ...NOT DROP SO ADD...
                  ELSEIF ( LISTNO .LT. LISTMX ) THEN
                     FOUND = T
                     DO 1155 K = 1, LISTNO
                        IF ( PNUMB .EQ. LIST(K) ) THEN
                           CALL ZZERRM ( RD(I), *91000,
     -                                   'NT'//ERR6//NMD//NAME )
                           GOTO 1156
                        ENDIF
 1155                CONTINUE

                     LISTNO = LISTNO + 1
                     LIST( LISTNO ) = PNUMB
 1156                CONTINUE

                  ELSE
                     GOTO 80000
                  ENDIF
                ENDIF
               ENDIF

 1160       CONTINUE
 1180       CONTINUE
         ENDIF
C             ...FOR NAME OR NUMBER.
      ENDIF
C          ...FOR CLASS = PROBS.

      IF ( .NOT. FOUND .AND. CLASS .EQ. PROBS ) THEN

         IF ( MODE .EQ. BYNUMB ) THEN
            CALL ZZERRM (RD(PNUMB),*91000,'IT'//ERR1//NUM)
         ELSE
            CALL ZZERRM ( RD(I), *91000,'NT'//ERR1//NMD//NAME)
         ENDIF

      ENDIF

C-----TRY TO FIND THE NAME OR NUMBER SPECIFIED AS A FUNCTION.

      IF ( CLASS .EQ. FUNCS  .OR. (   CLASS .EQ.  ANY
     -                               .AND.  .NOT. FOUND    )  ) THEN

C        CHECK EACH PROBLEM.

         DO 1210 I = 1, PROBNO

            IF ( PRECNO(RECN,I) .NE. -1 ) THEN

               IF ( MODE .EQ. BYNAME ) THEN
                  LOC = I * PNAMLN * 2
                  TRY = PFNAMS(LOC-PNAMLN+1:LOC) .EQ. NAME
               ELSE
                  TRY = PRECNO(FNO1,I) .EQ. PNUMB
               ENDIF

               IF ( TRY ) THEN

                  IF ( DROP ) THEN

                     DO 1200 J = 1,LISTNO
                        IF ( LIST(J) .EQ. I ) THEN
                           FOUND   = T
                           LIST(J) = -1
                           LEFT    = LEFT - 1
                           GOTO 1210
                        ENDIF
 1200                CONTINUE

                  ELSEIF ( LISTNO .GE. LISTMX ) THEN
                     GOTO 80000
                  ELSE

                     FOUND = T

                     DO 1205 K = 1, LISTNO

                        IF ( I .EQ. LIST(K) ) THEN
                           CALL ZZERRM( RD(I), *91000, 'IT'//ERR6//NUM )
                           GOTO 1206
                        ENDIF

 1205                CONTINUE

                     LISTNO = LISTNO + 1
                     LIST( LISTNO ) = I

 1206                CONTINUE

                  ENDIF
C                       FOR THE "IF ADD OR DROP..."
               ENDIF
C                    FOR THE "IF TRY..."
            ENDIF

 1210    CONTINUE

         IF ( .NOT. FOUND .AND. CLASS .EQ. FUNCS ) THEN

            IF ( MODE .EQ. BYNAME ) THEN
               CALL ZZERRM ( RD(I), *91000,'NT'//ERR2//NMD//NAME)
            ELSE
               CALL ZZERRM ( RD(PNUMB), * 91000, 'IT'// ERR2 //  NUM )
            ENDIF

         ENDIF

      ENDIF

C-----TRY TO FIND THE NAME SPECIFIED AS A GROUP.

      IF (  CLASS .EQ. GRPS .OR.
     -    ( CLASS .EQ. ANY   .AND. .NOT. FOUND ) ) THEN

         IF ( MODE .EQ. BYNUMB ) THEN

            IF ( PNUMB .GT. 0  .AND.  PNUMB .LE. NGRPS ) THEN

               IF ( GROUPS(PNUMB) .NE. -1 ) THEN

                  FOUND = T
                  LOC   = GNAMLN * PNUMB
                  NAME  = GNAMES(LOC-GNAMLN+1:LOC)

               ENDIF

            ENDIF

         ELSE
C                       NAMED GROUP SO SEARCH FOR IT
            PNUMB = 1

            CALL ZZSRCH(NAME,GNAMLN,GNAMES,NGRPS,GNAMLN,PNUMB, F, T)

            FOUND = PNUMB .NE. 0

         ENDIF

         IF ( .NOT. FOUND  .AND.  CLASS .EQ. GRPS ) THEN

            IF ( MODE .EQ. BYNUMB ) THEN
               CALL ZZERRM ( RD(PNUMB), * 91000, 'IT'// ERR3 // NUM )
            ELSE
               CALL ZZERRM ( RD(PNUMB), * 91000, 'NT'//ERR3//NMD//NAME )
            ENDIF

         ENDIF

 2100    IF ( FOUND ) THEN
C                          ADD A GROUP OF PROBLEMS.
            CALL ZZGETG ( DAUF, PNUMB, NAME, SIZE, MEMBS, GNAMES,
     -                    GROUPS, ERRFLG, *91000                 )

            IF ( ERRFLG .EQ. 0 ) THEN

               IF ( ADD ) THEN

                  DO 2200 I = 1, SIZE

                     DO 2150 K = 1, LISTNO

                        IF ( MEMBS(I) .EQ. LIST(K) ) THEN

                           CALL ZZERRM ( RD(MEMBS(I)), *91000,
     -                                  'IT'//ERR6//NUM       )
                           GOTO 2160

                        ENDIF

 2150                CONTINUE

                     IF ( LISTNO .GE. LISTMX ) THEN
                        GOTO 80000
                     ELSE
                        LISTNO = LISTNO + 1
                        LIST( LISTNO ) = MEMBS(I)
                     ENDIF

 2160                CONTINUE

 2200             CONTINUE

               ELSE

                  DO 4000 I = 1, SIZE

                     DO 3500 J = 1, LISTNO

                        IF ( LIST(J) .EQ. MEMBS(I) ) THEN
                           LIST(J) = -1
                           LEFT    = LEFT -1
                           GOTO 3600
                        ENDIF

 3500                CONTINUE
 3600                CONTINUE

 4000             CONTINUE

               ENDIF

            ENDIF

         ENDIF

      ENDIF

C-----FLAG ERROR IF NOT FOUND SOMEWHERE.

      IF ( .NOT. FOUND .AND. CLASS .EQ. ANY ) THEN

         IF ( MODE .EQ. BYNUMB ) THEN
            CALL ZZERRM ( RD(PNUMB), * 91000, 'IT'// ERR4 // NUM )
         ELSE
            CALL ZZERRM ( RD(PNUMB), * 91000, 'NT'//ERR4//NMD//NAME)
         ENDIF

      ENDIF

      GOTO 1000

C## E N T R Y  ZZASET:
                     ENTRY ZZASET ( SRADIX, SRADP, SRADF, SRADG )
      RADIX  = SRADIX
      RADPRB = SRADP
      RADFUN = SRADF
      RADGRP = SRADG

      RETURN

C## E R R O R S:

80000 IF ( MODE .EQ. BYNUMB ) THEN
         CALL ZZERRM ( RD(PNUMB), *91000, 'IT'// ERR5//NUM )
      ELSE
         CALL ZZERRM ( RD(PNUMB), * 91000, 'NT'//ERR5//NMD//NAME)
      ENDIF
      GOTO 90000

C## E X I T
90000      RETURN

C-----     ERROR RETURN.
91000      RETURN 1

C## F O R M A T S:  NONE ARE DEFINED.
C##                 E N D         OF ZZADDP.
                    END
