      PROGRAM CONVRT

C============================= S T A T U S =============================
C
C        SYSTEM  DEPENDENCE:  CHECK THE USE OF FILE NAMES,
C                          AS DESCRIBED BELOW.
C              SPECIAL FORMATS ARE USED IN FORMATS 99994 ,5, 7 AND 8.
C              A "$" IS USED TO PREVENT A CARRIAGE RETURN AFTER A LINE
C              IS SENT TO A TERMINAL. ALSO, CHECK THE DESCRIPTION BELOW
C              REGARDING THE LENGTH OF FILE NAMES AND USING THE SAME
C              NAME FOR INPUT AND OUTPUT FILES. FINALLY, YOU MAY NEED
C              TO REMOVE THE "CARRIAGE CONTROL = 'LIST'" IN THE OPEN
C              STATEMENT AT LABEL 25 IF YOU ARE NOT USING VAX/VMS.
C
C        CONVERSION FOR SINGLE     NOT REQUIRED
C        OR DOUBLE PRECISION
C
C        REVISION STATUS:      DATE         AUTHOR           PURPOSE
C
C                         DEC. 8, 1981   S. VINCENT   READY TO PUBLISH
C
C======================== D E S C R I P T I O N ========================
C
C     A MORE DETAILED DESCRIPTION IS GIVEN BELOW; HERE WE GIVE FIRST A
C     BRIEF DESCRIPTION OF HOW TO USE  ZZCNVT.
C
C     THE PROGRAM IS INTENDED TO CONVERT A PROGRAM WRITTEN IN SINGLE
C     PRECISION TO DOUBLE PRECISION (OR VICE VERSA), PROVIDING CERTAIN
C     RULES HAVE BEEN FOLLOWED IN PREPARING THE PROGRAM.  SPECIFICALLY,
C     FOR ANY PROGRAM STATEMENT WHICH MUST CHANGE WHEN THE PRECISION
C     CHANGES, THERE SHOULD BE A PAIR OF STATEMENTS, ONE FOR EACH PRE-
C     CISION.  ONE OR BOTH OF THEM MAY BE AS A COMMENT.  ANY OF THE
C     FORMATS  "A" TO "F"  DESCRIBED BELOW MAY BE USED.  THE PROGRAM
C     CONVERTS FROM ONE FORMAT TO ANOTHER. THE CHOICE OF FORMAT IS
C     LEFT SOMEWHAT TO INDIVIDUAL TASTE.
C
C     THE PROGRAM WRITES ITS REQUESTS FOR DIRECTIONS ON UNIT  6  AND
C     READS YOUR DIRECTIONS FROM  UNIT  5.   IF YOU ARE RUNNING
C     CONVERT INTERACTIVELY, UNIT 5 SHOULD THEREFORE BE ASSIGNED TO
C     YOUR TERMINAL AS INPUT AND UNIT 6 SHOULD BE ASSIGNED TO THE
C     TERMINAL AS OUTPUT.  IF YOU ARE DOING THIS FROM BATCH, THE
C     PROMPTS FOR INPUT WILL APPEAR ON UNIT 5, AND CAN BE IGNORED,
C     AND YOU SHOULD PREDEFINE UNIT 6 AND INSURE THAT IT CONTAINS
C     THE CORRECT DATA FOR EACH FILE TO BE CONVERTED, AS BELOW.
C
C     YOUR RESPONSE FOR EACH FILE SHOULD BE:
C
C      1. THE NAME OF THE FILE TO CONVERT.
C      2. THE NAME OF THE FILE TO RECEIVE THE CONVERTED PROGRAM.
C      3. THE FORMAT OF THE INPUT FILE:   A LETTER FROM "A" TO "F".
C      4. THE FORMAT FOR THE OUTPUT FILE: A LETTER FROM "A" TO "F".
C            (SEE THE STYLE TABLE BELOW.)
C
C     THE FOLLOWING POINTS SHOULD BE NOTED:
C
C      1. THE FILE NAMES SHOULD BE ACCEPTABLE FOR YOUR SYSTEM.
C
C      2. THE FILES DO NOT NEED TO BE DECLARED BEFORE RUNNING THE
C         PROGRAM OR IN THE PROGRAM STATEMENT.  THEY ARE EXPLICITLY
C         OPENED AND DEFINED IN CONVRT. UNITS 1 AND 2 ARE USED.
C
C      3. THE NAME OF THE OUTPUT FILE CAN BE OMITTED AND WILL BE
C         TAKEN TO BE THE SAME AS THE INPUT FILE. DO NOT USE THIS
C         FEATURE UNLESS IT IS ACCEPTABLE TO YOUR SYSTEM.
C
C      4. THE PROGRAM LOOPS SO THAT IT CAN BE USED TO CONVERT SEVERAL
C         PROGRAMS AT A TIME.  JUST  REPEAT THE FOUR INPUTS ABOVE.
C
C      5. IF THE LETTER "R" IS SPECIFIED FOR ONE OF THE FORMATS, IT
C         IS ASSUMED TO BE THE SAME AS THE PREVIOUS ONE, AND THAT
C         ALL SUBSEQUENT FORMATS WILL BE THE SAME, SO THAT NO FURTHER
C         REQUESTS ARE MADE FOR A FORMAT TYPE.  IF "R" IS SPECIFIED
C         ON THE FIRST REQUEST, THE DEFAULTS ARE  "A" FOR INPUT AND
C         "B" FOR OUTPUT. YOU DO NOT NEED TO SEPARATELY SPECIFY "R"
C         FOR INPUT AND OUTPUT.  IN BATCH, IF "R" IS SPECIFIED FOR
C         FORMAT, REMEMBER THAT FOR SUBSEQUENT FILE CONVERSIONS ONLY
C         TWO LINES OF INPUT SHOULD BE SUPPLIED, NOT FOUR.
C
C     A MORE DETAILED DESCRIPTION NOW FOLLOWS. IT CAN SAFELY BE IGNORED
C     UNLESS YOU WISH TO MODIFY CONVRT.
C
C     THE PURPOSE OF CONVRT IS TO CONVERT FORTRAN PROGRAMS WITH REGARD
C     TO THE USAGE OF SINGLE OR DOUBLE PRECISION REAL ARITHMETIC WHERE
C     ONE OR BOTH OF THE SINGLE AND DOUBLE PRECISION USAGE STATEMENTS
C     IN THE PROGRAM HAVE BEEN BLOCKED FROM COMPILATION BY A SPECIAL
C     BLOCKING CODE.  ( WHEN WE SPEAK OF A BLOCKING CODE WE INTEND A
C     STRING OF CHARACTERS WHICH STARTS IN COLUMN 1 AND EXTENDS TO
C     COLUMN CODELN - A PARAMETER TYPICALLY SET TO 5 ).  CONVERSIONS
C     MAY BE MADE IN TERMS OF MODE ( WHETHER SINGLE PRECISION, DOUBLE
C     PRECISION, OR BOTH OF THESE ARE BLOCKED FROM COMPILATION BY A
C     BLOCKING CODE ), FORMAT ( THE FORMAT OF THE BLOCKING CODES USED),
C     OR BOTH SIMULTANEOUSLY.  AS THIS PROGRAM STANDS, THE
C     5 DECEMBER 81 VERSION, TWO FORMATS ARE INCLUDED, WHICH PRODUCES
C     SIX COMBINATIONS OF MODES AND FORMATS WHICH MAY BE IN THE CURRENT
C     PROGRAM AND WHICH MAY BE SELECTED FOR THE REVISED PROGRAM :
C
C        +-------------------------+-------------------------+
C        | (A) :       DOUBLE ...  | (D) :       DOUBLE ...  |
C        |       C!!!! REAL   ...  |       C!SNG REAL   ...  |
C        +-------------------------+-------------------------+
C        | (B) :       REAL   ...  | (E) :       REAL   ...  |
C        |       C!!!! DOUBLE ...  |       C!DBL DOUBLE ...  |
C        +-------------------------+-------------------------+
C        | (C) : C!!!! REAL   ...  | (F) : C!SNG REAL   ...  |
C        |       C!!!! DOUBLE ...  |       C!DBL DOUBLE ...  |
C        +-------------------------+-------------------------+
C
C     THE USER MUST ENTER THE NAMES OF THE SOURCE FILE AND THE
C     FILE TO HOLD THE REVISED PROGRAM.  THE NUMBER OF CHARACTERS
C     USED IN THE FILE NAMES IS SET BY THE PARAMETER LNFNAM AND
C     SHOULD BE MODIFIED AS APPROPRIATE FOR THE NAMES ENCOUNTERED
C     WITH YOUR FACILITY.  THE 5 DECEMBER 81 VERSION SETS THE
C     OUTPUT NAME TO BE THE INPUT FILE NAME IF NO RESPONSE IS
C     MADE TO THE PROMPT FOR IT.  IF THIS PRESENTS A PROBLEM FOR
C     YOUR OPERATING SYSTEM THEN THE CODE STARTING AT THE LABEL
C     13000 SHOULD BE MODIFIED TO RESEMBLE THAT AT LABEL 12000.
C
C     ONCE THE FILES HAVE BEEN SUCCESSFULLY OPENED, CONVRT
C     PROMPTS THE USER FOR THE CURRENT AND DESIRED STYLES OF
C     THE PROGRAM.  THE TASKS WHICH ARE PERFORMED SUBSEQUENT
C     TO SUCCESSFUL SPECIFICATION OF THE TWO SYLES MAY BE
C     SUMMARIZED :
C
C                      |     D E S I R E D   S T Y L E       |
C                      |                                     |
C                      |   SINGLE   |   DOUBLE   |    BOTH   |
C                      |  BLOCKED   |  BLOCKED   |  BLOCKED  |
C        --------------+------------+------------+-----------+
C        C     SINGLE  |   CHANGE   |   SWITCH   |    ADD    |
C        U     BLOCKED |   FORMAT   |  BLOCKERS  |  BLOCKER  |
C        R     --------+------------+------------+-----------+
C        R S   DOUBLE  |   SWITCH   |   CHANGE   |    ADD    |
C        E T   BLOCKED |  BLOCKERS  |   FORMAT   |  BLOCKER  |
C        N Y   --------+------------+------------+-----------+
C        T L   BOTH    |   DELETE   |   DELETE   |   CHANGE  |
C          E   BLOCKED |  BLOCKER   |  BLOCKER   |   FORMAT  |
C        --------------+------------+------------+-----------+
C
C     WE NOTE THAT FORMAT CHANGES MAY OCCUR WITH ANY OF THE NINE
C     POSSIBLE CURRENT AND DESIRED STYLE COMBINATIONS, ALTHOUGH NOT
C     EXPLICITLY STATED IN THE ABOVE TABLE. AFTER A SUCCESSFUL
C     CONVERSION CONVRT WILL RETURN TO THE BEGINNING OF THE PROGRAM
C     AND ASK FOR THE NAME OF THE INPUT FILE. IF YOU DO NOT SPECIFY A
C     FILE THEN THE PROGRAM STOPS.  ONCE THE FORMATS FOR THE CURRENT
C     AND DESIRED STYLES HAVE BEEN SET THE USER MAY SPECIFY THE
C     OPTION AS 'R' FOR REPEAT AT THE EITHER PROMPT FOR STYLE AND
C     HAVE THE PREVIOUS COMPLETELY SPECIFIED FORMATS USED FOR ALL OF
C     THE REMAINING FILES IN THE RUN.
C
C
C     O T H E R   I M P L E M E N T A T I O N   N O T E S   :
C
C     THE USER SHOULD NOTE THAT THE FOLLOWING LOGICAL UNIT NUMBERS
C     ARE USED BY CONVRT :
C
C           UNIT TYPE       INPUT NUMBER     OUTPUT NUMBER
C        ---------------    ------------     -------------
C        USER'S TERMINAL         5                 6
C           PROGRAMS             1                 2
C
C     THE UNIT NUMBERS MAY BE CHANGED BY CHANGING THE PARAMETER
C     STATEMENT BELOW.
C
C     THE NUMBER OF CHARACTERS WHICH MAY BE READ FROM THE SOURCE
C     FILE AND WRITTEN TO THE REVISED FILE IS SET BY THE PARAMETER
C     LINELN ( TYPICALLY 80 ).
C
C     SHOULD THE FORMATS OF THE BLOCKERS BE CHANGED OR THE NUMBER
C     OF FORMATS BE CHANGED THE PERSON MODIFYING THE CODE SHOULD :
C
C        (1)  CHANGE THE PARAMETER VALUE OF NUMFRM IF THE NUMBER OF
C             FORMS BE CHANGED.
C        (2)  ALTER OR ADD THE CODES FOR THE SINGLE AND DOUBLE
C             PRECISION BLOCKERS FOR THE APPROPRIATE FORMS.  THIS
C             ENTAILS DEFINING PARAMETERS LIKE COD1, COD21, COD22, ETC.
C             AND MODIFYING THE DATA STATEMENT WHICH MAKES USE OF THESE
C             VALUES.  FOR EACH FORM THE VALUE OF MONO MUST BE SET -
C             MONO ( K ) IS .TRUE. IF BOTH THE SINGLE AND DOUBLE
C             PRECISION BLOCKER CODES FOR FORM K ARE THE SAME.
C        (3)  MODIFY THE MENU FORMAT STATEMENT (99996) AND THE LIST
C             OF PERMISSIBLE OPTIONS (OPTLIS) AS APPROPRIATE.
C        (4)  ALTER THIS DESCRIPTION OF CONVRT.
C
C======================== S U B R O U T I N E S ========================
C
C     NO SUBROUTINES ARE USED.
C
C========================= P A R A M E T E R S =========================

C-----SET LOGICAL UNIT NUMBERS.

      INTEGER       INPUT    ,  OUTPUT    ,  TRMIN    ,  TRMOUT
      PARAMETER   ( INPUT = 1,  OUTPUT = 2,  TRMIN = 5,  TRMOUT = 6 )

C-----DECLARE THE LENGTH ALLOWED FOR THE INPUT AND OUTPUT FILE NAMES.

      INTEGER       LNFNAM
      PARAMETER   ( LNFNAM = 80 )

C-----DECLARE THE LENGTH OF LINES EXPECTED IN THE SOURCE FILE.

      INTEGER       LINELN
      PARAMETER   ( LINELN = 80 )

C-----SET PARAMETERS FOR NUMBER OF FORMATS, NUMBER OF MODES, AND THE
C-----NUMBER OF OPTIONS. SET THE LENGTH OF THE BLOCKING CODES, THE
C-----BLOCKING CODES AND THE MENU OPTION LIST.

      INTEGER       CODELN     , NUMFRM     ,  NUMMOD
      PARAMETER   ( CODELN = 5 , NUMFRM = 2 ,  NUMMOD = 3 )

      INTEGER       NUMOPT
      PARAMETER   ( NUMOPT = NUMMOD * NUMFRM )

      CHARACTER*(*) COD1,           COD21,           COD22
      PARAMETER   ( COD1 = 'C!!!!', COD21 = 'C!SNG', COD22 = 'C!DBL' )

      CHARACTER *(NUMOPT)  OPTLIS
      PARAMETER          ( OPTLIS = 'ABCDEF' )

C-----SET CONSTANTS FOR THE THREE MODES OF BLOCKING.

      INTEGER      SINGLE    ,  DOUBLE    ,  BOTH
      PARAMETER  ( SINGLE = 1,  DOUBLE = 2,  BOTH = 3 )

C-----SET SIMPLE CHARACTER VARIABLES AND CONSTANTS.

      CHARACTER*(*) BLANK,       REPEAT
      PARAMETER   ( BLANK = ' ', REPEAT = 'R' )

C======================= D E C L A R A T I O N S =======================

C----ARRAYS AND STRINGS.

      LOGICAL     MONO ( NUMFRM )

      CHARACTER * ( LNFNAM )  INNAM,  OUTNAM

      CHARACTER * ( CODELN )  CODDBL ( NUMFRM ),  CODSNG ( NUMFRM )

      CHARACTER * ( CODELN )  MODEL,  TARGET,  NEWTAR

      CHARACTER * ( LINELN )  LINE1,  LINE2

C.......................................................................

C----MISCELLANEOUS VARIABLE DECLARATIONS.

      INTEGER  ERRCNT,  FORM  ,  GOBACK,  K     ,  MODE  ,  I, ITMP
      INTEGER  NEWFRM,  NEWMOD,  OLDFRM,  OLDMOD,  LEN1  ,  LEN2

      INTEGER  ZZLENG, ZZLEFT

      LOGICAL  ASK   ,  FIRST ,  HAVNEW,  HAVOLD,  MONOFM,  REFORM
      LOGICAL  REVERS,  START ,  TASK1 ,  TASK3

      CHARACTER * ( 1 )       OPTION, OLDINP, OLDOUT

C=============================== S A V E ===============================
C
C                     THERE ARE NO  SAVE  VARIABLES.
C
C=============================== D A T A ===============================

C-----DATA STATEMENTS FOR BLOCKING CODES AND STATUS OF FORMS (MONO).

      DATA     CODSNG ( 1 ),  CODSNG ( 2)
     -       / COD1        ,  COD21       /

      DATA     CODDBL ( 1 ),  CODDBL ( 2 )
     -       / COD1        ,  COD22       /

      DATA     MONO ( 1 )  ,  MONO ( 2 )
     -       / .TRUE.      ,  .FALSE.     /

C========================== E X E C U T I O N ==========================

C-----GET THE INPUT AND OUTPUT FILE NAMES.

      OLDINP = 'A'
      OLDOUT = 'B'
      ASK    = .TRUE.
    1 ASSIGN 10 TO GOBACK
   10 WRITE ( TRMOUT, 99998 )
      READ  ( TRMIN , 99999, END = 90000, ERR = 11000 ) INNAM
      DO 12 I = 1, LNFNAM
         IF ( INNAM ( I:I ) .NE. BLANK ) THEN
            GO TO 15
         ENDIF
   12 CONTINUE
      GOTO 90000

   15 I = ZZLEFT(INNAM)
      OPEN  ( UNIT = INPUT, FILE = INNAM(1:I),STATUS='OLD',
     -        ERR  = 12000 )
      WRITE ( TRMOUT, '(A)') 'OPENED '//INNAM(1:I)

      ASSIGN 20 TO GOBACK
   20 WRITE ( TRMOUT, 99997 )
      READ  ( TRMIN , 99999, END = 13000, ERR = 11000 ) OUTNAM
      DO 22 I = 1, LNFNAM
         IF ( OUTNAM ( I:I ) .NE. BLANK ) THEN
            GO TO 25
         ENDIF
   22 CONTINUE
      GO TO 13000

   25 I = ZZLEFT(OUTNAM)
      OPEN  ( UNIT = OUTPUT,FILE=OUTNAM(1:I),STATUS='NEW',
     -         ERR  = 14000 )
      WRITE ( TRMOUT, '(A)') 'OPENED '//OUTNAM(1:I)

C-----INITIALIZE PROGRAM.

      REWIND INPUT
      REWIND OUTPUT

      IF ( ASK ) THEN

         HAVOLD = .FALSE.
         HAVNEW = .FALSE.

C        ----HAVE THE USER IDENTIFY THE CURRENT (OLD) AND DESIRED
C            (NEW) FORMS AND MODES OF BLOCKERS.

   30    ERRCNT =  0
         ASSIGN 40 TO GOBACK

         WRITE ( TRMOUT, 99996 )

   40    IF ( .NOT. HAVOLD ) THEN
            WRITE ( TRMOUT, 99995 )
         ELSE
            WRITE ( TRMOUT, 99994 )
         ENDIF

         READ ( TRMIN, 99999, END = 10000, ERR = 11000 )  LINE1
         ITMP = ZZLEFT(LINE1)
         OPTION = LINE1(1:1)

         IF ( OPTION .EQ. REPEAT  ) THEN
            ASK = .FALSE.
            IF ( HAVOLD ) THEN
               OPTION = OLDOUT
            ELSE
               OPTION = OLDINP
            ENDIF
         ENDIF

   45    MODE = INDEX ( OPTLIS, OPTION )

         IF ( MODE .GT. 0 ) THEN
            FORM = 1
   50       IF ( MODE .GT. NUMMOD ) THEN
               FORM = FORM + 1
               MODE = MODE - NUMMOD
               GO TO 50
            ENDIF

            IF ( .NOT. HAVOLD ) THEN
               HAVOLD = .TRUE.
               OLDFRM =  FORM
               OLDMOD =  MODE
               IF ( ASK ) THEN
                  GOTO 40
               ELSE
                  OPTION = OLDOUT
                  GOTO 45
               ENDIF
            ELSE
               HAVNEW = .TRUE.
               NEWFRM =  FORM
               NEWMOD =  MODE
               GO TO 60
            ENDIF
         ELSE
            ERRCNT = ERRCNT + 1
            WRITE ( TRMOUT, 99993 )
            IF ( ERRCNT .GT. 4 ) THEN
               GO TO 30
            ELSE
               GO TO 40
            ENDIF
         ENDIF

      ENDIF
C             FOR THE "IF ASK..."

C-----IDENTIFY THE TASK TO BE PERFORMED.

   60 MONOFM = MONO ( OLDFRM )
      TASK3  = .FALSE.
      IF ( OLDMOD .EQ. NEWMOD ) THEN

         IF ( OLDFRM .EQ. NEWFRM ) THEN

C        ---NOTHING TO DO.

            WRITE ( TRMOUT, 99992 )
            GO TO 90000

         ELSE

C        ---CHANGE FORMATS ONLY.

            IF ( OLDMOD .EQ. BOTH ) THEN
               IF ( MONOFM ) THEN
                  TARGET = CODSNG ( OLDFRM )
                  NEWTAR = CODSNG ( NEWFRM )
                  MODEL  = CODDBL ( NEWFRM )
                  GO TO 200
               ELSE
                  TASK1  = .TRUE.
                  REVERS = .FALSE.
                  REFORM = .TRUE.
                  TARGET = CODDBL ( OLDFRM )
                  NEWTAR = CODDBL ( NEWFRM )
                  MODEL  = CODSNG ( NEWFRM )
                  GO TO 100
               ENDIF

            ELSEIF ( OLDMOD .EQ. DOUBLE ) THEN
               TARGET = CODDBL ( OLDFRM )
               NEWTAR = CODDBL ( NEWFRM )
               MODEL  = NEWTAR
               GO TO 200

            ELSE
               TARGET = CODSNG ( OLDFRM )
               NEWTAR = CODSNG ( NEWFRM )
               MODEL  = NEWTAR
               GO TO 200
            ENDIF

         ENDIF

      ELSEIF ( OLDMOD .EQ. BOTH ) THEN

C     ---DELETE ONE OF THE BLOCKERS.

         IF ( MONOFM ) THEN
            TARGET = CODSNG ( OLDFRM )
            IF ( NEWMOD .EQ. DOUBLE ) THEN
               NEWTAR = BLANK
               MODEL  = CODDBL ( NEWFRM )
               GO TO 200
            ELSE
               NEWTAR = CODSNG ( NEWFRM )
               FIRST  = .TRUE.
               TASK3  = .TRUE.
               GO TO 100
            ENDIF

         ELSE
            TASK1  = .FALSE.
            TARGET = CODDBL ( OLDFRM )
            REVERS = NEWMOD .EQ. SINGLE
            IF ( REVERS ) THEN
               MODEL = CODSNG ( NEWFRM )
            ELSE
               MODEL = CODDBL ( NEWFRM )
            ENDIF
            GO TO 100
         ENDIF

      ELSEIF ( NEWMOD .EQ. BOTH ) THEN

C     ---ADD ONE OF THE BLOCKERS.

         TASK1  = .TRUE.
         REFORM = NEWFRM .NE. OLDFRM
         IF ( OLDMOD .EQ. SINGLE ) THEN
            TARGET = CODSNG ( OLDFRM )
            MODEL  = CODDBL ( NEWFRM )
            REVERS = .TRUE.
            IF ( REFORM ) THEN
               NEWTAR = CODSNG ( NEWFRM )
            ENDIF

         ELSE
            TARGET = CODDBL ( OLDFRM )
            MODEL  = CODSNG ( NEWFRM )
            REVERS = .FALSE.
            IF ( REFORM ) THEN
               NEWTAR = CODDBL ( NEWFRM )
            ENDIF
         ENDIF
         GO TO 100

      ELSE

C     ---SWITCH BLOCKERS.

         TASK1  = .FALSE.
         REVERS = .TRUE.
         IF ( OLDMOD .EQ. SINGLE ) THEN
            TARGET = CODSNG ( OLDFRM )
            MODEL  = CODDBL ( NEWFRM )
         ELSE
            TARGET = CODDBL ( OLDFRM )
            MODEL  = CODSNG ( NEWFRM )
         ENDIF
         GO TO 100
      ENDIF

C-----TWO LINE PROCESSING IS REQUIRED FOR THE TASK.

  100 START = .TRUE.

  110 IF ( START ) THEN

C     ---GET A NEW LEADING LINE.

         READ ( INPUT, 99999, END = 190, ERR = 15000 ) LINE1
         DO 120 K = LINELN, 1, -1
            IF ( LINE1 ( K:K ) .NE. BLANK ) THEN
               LEN1 = K
               GO TO 130
            ENDIF
  120    CONTINUE
         LEN1 = 1
      ENDIF

C-----(ALWAYS) GET A NEW TRAILING LINE.

  130 READ ( INPUT, 99999, END = 180, ERR = 15000 ) LINE2
      DO 140 K = LINELN, 1, -1
         IF ( LINE2 ( K:K ) .NE. BLANK ) THEN
            LEN2 = K
            GO TO 150
         ENDIF
  140 CONTINUE
      LEN2 = 1

  150 IF ( LINE2 ( 1:CODELN ) .EQ. TARGET ) THEN

C     ---PERFORM THE NECESSARY CONVERSION.

         IF ( .NOT. TASK3 ) THEN
            IF ( TASK1 ) THEN

C           ---WE ARE ADDING A BLOCKER OR CHANGING THE
C           ---FORMAT OF BOTH BLOCKERS.

               IF ( REFORM ) THEN

C              ---REVISE THE FORMAT OF THE TARGET.

                  LINE2 ( 1:CODELN ) = NEWTAR
               ENDIF

C           ---ADD OR REVISE THE FORMAT OF THE LEADING
C           ---BLOCKER.

               LINE1 ( 1:CODELN ) = MODEL

            ELSE

C           ---FOR DELETING OR SWITCHING BLOCKERS.

               IF ( REVERS ) THEN

C              ---SWITCH BLOCKERS.

                  LINE2 ( 1:CODELN ) = LINE1 ( 1:CODELN )
                  LINE1 ( 1:CODELN ) = MODEL

               ELSE

C              ---DELETE THE LEADING AND POSSIBLE REVISE
C              ---THE FORMAT OF THE TRAILING BLOCKER.

                  LINE1 ( 1:CODELN ) = BLANK
                  LINE2 ( 1:CODELN ) = MODEL

               ENDIF

            ENDIF

            IF ( REVERS ) THEN

C           ---SWITCH THE LINES.

               WRITE ( OUTPUT, 99999 )  LINE2 ( 1:LEN2 )
               WRITE ( OUTPUT, 99999 )  LINE1 ( 1:LEN1 )

            ELSE

C           ---DO NOT SWITCH THE LINES.

               WRITE ( OUTPUT, 99999 )  LINE1 ( 1:LEN1 )
               WRITE ( OUTPUT, 99999 )  LINE2 ( 1:LEN2 )

            ENDIF
            START = .TRUE.

         ELSE

C        ---TASK 3 : GOING FROM BOTH BLOCKERS WITH MONO FORM
C        ---TO SINGLE BLOCKER.

            IF ( FIRST ) THEN

C           ---LINE2 CONTAINS THE SINGLE PRECISION BLOCKER,
C           ---WRITE OUT LINE1 AND ALTER, IF NECESSARY,
C           ---THE BLOCKER CODE.  UNLIKE THE OTHER TASKS
C           ---USING TWO LINE PROCESSING, WE DO NOT NOW
C           ---READ IN TWO NEW LINES.

               LINE2 ( 1:CODELN ) = NEWTAR
               WRITE ( OUTPUT, 99999 )  LINE1 ( 1:LEN1 )
               LINE1 = LINE2
               LEN 1 = LEN 2
               START = .FALSE.

            ELSE

C           ---LINE2 CONTAINS THE DOUBLE PRECISION BLOCKER,
C           ---DELETE THIS BLOCKER AND THEN WRITE OUT LINES
C           ---1 AND 2 IN REVERSE ORDER.  WE NOW NEED TO
C           ---READ IN TWO NEW LINES AS IS USUAL FOR THE
C           ---OTHER TASKS.

               LINE2 ( 1:CODELN ) = BLANK
               WRITE ( OUTPUT, 99999 )  LINE2 ( 1:LEN2 )
               WRITE ( OUTPUT, 99999 )  LINE1 ( 1:LEN1 )
               START = .TRUE.
            ENDIF
            FIRST = .NOT. FIRST
         ENDIF

      ELSE

C     ---WRITE, BUT DO NOT ALTER THE LEADING LINE AND
C     ---SET THE LEADING LINE TO BE THE TRAILING LINE.

         WRITE ( OUTPUT, 99999 )  LINE1 ( 1:LEN1 )
         LINE1 = LINE2
         LEN 1 = LEN 2
         START = .FALSE.

      ENDIF
      GO TO 110

  180 WRITE ( OUTPUT, 99999 )  LINE1 ( 1:LEN1 )
  190 CONTINUE
      GO TO 1000

C-----ONE LINE PROCESSING FOR REVISING THE FORMAT OF SINGLE
C-----BLOCKERS OR FOR REVISING THE FORMAT OF DOUBLE BLOCKERS WITH
C-----THE SAME BLOCKING CODE.

  200 START = .TRUE.
  210 READ ( INPUT, 99999, END = 290, ERR = 15000 ) LINE1
      DO 220 K = LINELN, 1, -1
         IF ( LINE1 ( K:K ) .NE. BLANK ) THEN
            LEN1 = K
            GO TO 230
         ENDIF
  220 CONTINUE
      LEN1 = 1

  230 IF ( LINE1 ( 1:CODELN ) .EQ. TARGET ) THEN

C     ---REPLACE TARGET.

         IF ( START ) THEN
            LINE1 ( 1:CODELN ) = NEWTAR
         ELSE
            LINE1 ( 1:CODELN ) = MODEL
         ENDIF
         START = .NOT. START

      ENDIF

      WRITE ( OUTPUT, 99999 )  LINE1 ( 1:LEN1 )

      GO TO 210

  290 WRITE ( OUTPUT, 99999 )  LINE1 ( 1:LEN1 )
      GO TO 1000

C-----CLOSE FILES.

 1000 CLOSE ( UNIT = INPUT  )
      CLOSE ( UNIT = OUTPUT )
      GO TO 80000

C.......................................................................

C-----THIS BIT OF CODE IS PROVIDED FOR SYSTEMS ON WHICH NO INPUT
C-----FROM A TERMINAL IS INTERPRETED AS AN END-OF-FILE FOR THE
C-----ASSOCIATED LOGICAL UNIT.

10000 WRITE ( TRMOUT, 99991 )
10100 CLOSE ( UNIT = TRMIN )
      OPEN  ( UNIT = TRMIN, STATUS = 'UNKNOWN' )
      GO TO GOBACK

C.......................................................................

C-----AN ERROR WAS ENCOUNTERED IN READING FROM THE TERMINAL, ATTEMPT
C-----TO CONTINUE.

11000 WRITE ( TRMOUT, 99990 )
      GO TO GOBACK

C.......................................................................

C-----AN ERROR WAS ENCOUNTERED IN OPENING THE INPUT FILE.

12000 WRITE ( TRMOUT, 99989 ) INNAM
      GO TO GOBACK

C.......................................................................

C-----AN END-OF-FILE WAS ENCOUNTERED IN READING THE OUTPUT FILE NAME.
C-----THIS CODE MUST BE MODIFIED IF WE CAN NOT ASSIGN THE INPUT FILE
C-----NAME TO THE OUTPUT FILE (I.E. VERSIONS ARE POSSIBLE).

13000 OUTNAM = INNAM
      ASSIGN 25 TO GOBACK
      GO TO 25

C.......................................................................

C-----AN ERROR WAS ENCOUNTERED IN OPENING THE OUTPUT FILE.
C-----THIS SHOULD NOT OCCUR WHEN VERSIONS ARE POSSIBLE AND THE OUTPUT
C-----FILE NAME HAS BEEN SET TO THE INPUT.

14000 WRITE ( TRMOUT, 99988 ) OUTNAM
      GO TO GOBACK

C.......................................................................

C-----AN ERROR WAS ENCOUNTERED IN READING FROM THE INPUT FILE.  ABORT
C-----THE RUN

15000 WRITE ( TRMOUT, 99987 ) INNAM
      GO TO 1000

C.......................................................................

C-----CONTINUE POINT.

80000 GO TO 1

C=============================== E X I T ===============================

90000 STOP

C============================ F O R M A T S ============================

99987 FORMAT ( 1H , 'AN ERROR WAS ENCOUNTERED IN READING FROM THE',
     -              ' FILE ', A, /,
     -              ' RUN IS ABORTED.' )

99988 FORMAT ( 1H , '<', A , '>', /,
     -              ' IS NOT A PROPER NAME FOR AN OUTPUT FILE.', /,
     -              ' PLEASE TRY AGAIN.' )

99989 FORMAT ( 1H , '<', A , '>', /,
     -              ' IS NOT A PROPER NAME FOR AN INPUT FILE.', /,
     -              ' PLEASE TRY AGAIN.' )

99990 FORMAT ( 1H , 'THE FORTRAN SYSTEM DETECTED AN ERROR IN READING.',
     -         /,   ' PLEASE TRY AGAIN.' )
99991 FORMAT ( 1H , 'PLEASE RESPOND.' )

99992 FORMAT ( 1H , 'NOTHING TO DO.' )

99993 FORMAT ( 1H , 'ILLEGAL STYLE, TRY AGAIN.' )

99994 FORMAT ( 1H , /, ' ENTER THE DESIRED STYLE : ' )

99995 FORMAT ( 1H  / ' ENTER THE CURRENT STYLE : ' )

99996 FORMAT ( 1H0,
     -  'STYLES -'                                              /
     - ' +-------------------------+-------------------------+' /
     - ' | (A) :       DOUBLE ...  | (D) :       DOUBLE ...  |' /
     - ' |       C!!!! REAL   ...  |       C!SNG REAL   ...  |' /
     - ' +-------------------------+-------------------------+' /
     - ' | (B) :       REAL   ...  | (E) :       REAL   ...  |' /
     - ' |       C!!!! DOUBLE ...  |       C!DBL DOUBLE ...  |' /
     - ' +-------------------------+-------------------------+' /
     - ' | (C) : C!!!! REAL   ...  | (F) : C!SNG REAL   ...  |' /
     - ' |       C!!!! DOUBLE ...  |       C!DBL DOUBLE ...  |' /
     - ' +-------------------------+-------------------------+'   )

99997 FORMAT ( 1H , /, ' ENTER THE NAME OF THE OUTPUT FILE : ' )

99998 FORMAT ( 1H , /, ' ENTER THE NAME OF THE INPUT FILE : ' )

99999 FORMAT (  A  )

C================================ E N D ================================

      END
      INTEGER FUNCTION ZZLENG (LINE)

C============== A R G U M E N T   D E C L A R A T I O N S ==============

      CHARACTER*(*) LINE

C============================= S T A T U S =============================
C
C    SYSTEM  DEPENDENCE:               NONE
C
C    CONVERSION FOR SINGLE OR DOUBLE PRECISION:   NOT REQUIRED
C
C    REVISION STATUS:  DATE        AUTHOR           VERSION
C
C                 JUN.  1, 1985  A. BUCKLEY           1.0
C
C======================== D E S C R I P T I O N ========================
C
C     THIS ROUTINE DETERMINES THE POSITION OF THE LAST NONBLANK
C     CHARACTER IN THE STRING LINE. IF THE LINE IS ENTIRELY
C     BLANK, THEN ZZLENG IS SET TO 0.
C
C======================= E N T R Y   P O I N T S =======================
C
C                  ONLY THE NATURAL ENTRY POINT ZZLENG
C
C======================== S U B R O U T I N E S ========================
C
C     LEN    ...INTRINSIC
C
C========================= P A R A M E T E R S =========================

      CHARACTER*(*)  BLANK
      PARAMETER    ( BLANK = ' ' )

C================= L O C A L   D E C L A R A T I O N S =================

      INTEGER   I

C=============================== S A V E ===============================
C
C                     THERE ARE NO  SAVE  VARIABLES.
C
C============================= C O M M O N =============================
C
C                       THERE ARE NO COMMON BLOCKS.
C
C=============================== D A T A ===============================
C
C                        THERE ARE NO DATA VALUES.
C
C========================== E X E C U T I O N ==========================

      ZZLENG = 0

      DO 1000  I = LEN(LINE), 1, -1

         IF ( LINE(I:I) .NE. BLANK ) THEN
            ZZLENG = I
            GOTO 90000
         ENDIF

 1000 CONTINUE

C=============================== E X I T ===============================

90000 RETURN

C============================ F O R M A T S ============================
C
C                       THERE ARE NO FORMATS USED.
C
C================================ E N D ================================

      END
      INTEGER FUNCTION ZZLEFT (LINE)

C============== A R G U M E N T   D E C L A R A T I O N S ==============

      CHARACTER*(*) LINE

C============================= S T A T U S =============================
C
C    SYSTEM  DEPENDENCE:               NONE
C
C    CONVERSION FOR SINGLE OR DOUBLE PRECISION:   NOT REQUIRED
C
C    REVISION STATUS:  DATE        AUTHOR           VERSION
C
C                 JULY  14, 1987  A. BUCKLEY           1.0
C
C======================== D E S C R I P T I O N ========================
C
C     THIS ROUTINE DETERMINES THE POSITION OF THE LAST NONBLANK
C     CHARACTER IN THE STRING LINE. IF THE LINE IS ENTIRELY
C     BLANK, THEN ZZLEFT IS SET TO 0. NOTE THAT IT FIRST CLEARS
C     ANY LEADING BLANKS BY FIRST LEFT SHIFTING THE LINE IF NEEDED.
C
C======================= E N T R Y   P O I N T S =======================
C
C                  ONLY THE NATURAL ENTRY POINT ZZLEFT
C
C======================== S U B R O U T I N E S ========================
C
C     LEN    ...INTRINSIC
C
C========================= P A R A M E T E R S =========================

      CHARACTER*(*)  BLANK
      PARAMETER    ( BLANK = ' ' )

C================= L O C A L   D E C L A R A T I O N S =================

      INTEGER   I, START

C=============================== S A V E ===============================
C
C                     THERE ARE NO  SAVE  VARIABLES.
C
C======================= E Q U I V A L E N C E S =======================
C
C                       THERE ARE NO EQUIVALENCES.
C
C============================= C O M M O N =============================
C
C                       THERE ARE NO COMMON BLOCKS.
C
C=============================== D A T A ===============================
C
C                        THERE ARE NO DATA VALUES.
C
C========================== E X E C U T I O N ==========================

      ZZLEFT = 0

      DO 500 I = 1, LEN(LINE)
         IF ( LINE(I:I) .NE. BLANK ) THEN
            START = I
            GOTO 600
         ENDIF
  500 CONTINUE

      GOTO 90000

  600 DO 700  I = LEN(LINE), 1, -1

         IF ( LINE(I:I) .NE. BLANK ) THEN
            ZZLEFT = I
            GOTO 800
         ENDIF

  700 CONTINUE

  800 IF ( START .NE. 1 ) THEN
         CALL ZZSHFT ( LINE, START, 1, ZZLEFT-START+1 )
         ZZLEFT = ZZLEFT - START + 1
      ENDIF

 1000 CONTINUE


C=============================== E X I T ===============================

90000 RETURN

C============================ F O R M A T S ============================
C
C                       THERE ARE NO FORMATS USED.
C
C================================ E N D ================================

      END
      SUBROUTINE ZZSHFT (STRING, FROM, TO, NUMBER )

C============== A R G U M E N T   D E C L A R A T I O N S ==============

      INTEGER        FROM,    TO,    NUMBER

      CHARACTER *(*) STRING

C============================= S T A T U S =============================
C
C    SYSTEM  DEPENDENCE:               NONE
C
C    CONVERSION FOR SINGLE OR DOUBLE PRECISION:   NOT REQUIRED
C
C    REVISION STATUS:  DATE        AUTHOR           VERSION
C
C                 JUN.  1, 1985  A. BUCKLEY           1.0
C
C======================== D E S C R I P T I O N ========================
C
C     THIS ROUTINE PERFORMS A SHIFT OF CHARACTERS WITHIN STRING. THE
C     NUMBER OF CHARACTERS SHIFTED IS NUMBER AND THEY ARE SHIFTED SO
C     THAT THE CHARACTER IN POSITION FROM IS MOVED TO POSITION TO.
C     CHARACTERS IN THE TO POSITION ARE OVERWRITTEN. BLANKS REPLACE
C     CHARACTERS IN THE FROM POSITION. SHIFTING MAY BE LEFT OR RIGHT,
C     AND THE FROM AND TO POSITIONS MAY OVERLAP.  CARE IS TAKEN NOT
C     TO ALTER OR USE ANY CHARACTERS BEYOND THE DEFINED LIMITS
C     OF THE STRING.
C
C======================= E N T R Y   P O I N T S =======================
C
C                  ONLY THE NATURAL ENTRY POINT ZZSHFT
C
C======================== S U B R O U T I N E S ========================
C
C     LEN  MIN  MAX      ...INTRINSIC
C
C========================= P A R A M E T E R S =========================

      CHARACTER *(*)  BLANK
      PARAMETER     ( BLANK = ' ' )

C================= L O C A L   D E C L A R A T I O N S =================

      INTEGER        N, SHIFT, INCR, I, IS, IE, IBS, ETO, EFROM, K, SLEN

C=============================== S A V E ===============================
C
C                     THERE ARE NO  SAVE  VARIABLES.
C
C======================= E Q U I V A L E N C E S =======================
C
C                       THERE ARE NO EQUIVALENCES.
C
C============================= C O M M O N =============================
C
C                       THERE ARE NO COMMON BLOCKS.
C
C=============================== D A T A ===============================
C
C                        THERE ARE NO DATA VALUES.
C
C========================== E X E C U T I O N ==========================

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
            K = I - SHIFT
            STRING(K:K) = STRING(I:I)
 1000    CONTINUE

         DO 2000 I=IBS,IE,INCR
            STRING(I:I) = BLANK
 2000    CONTINUE

      ENDIF

      GOTO 90000

C=============================== E X I T ===============================

90000 RETURN

C============================ F O R M A T S ============================
C
C                       THERE ARE NO FORMATS USED.
C
C================================ E N D ================================

      END
