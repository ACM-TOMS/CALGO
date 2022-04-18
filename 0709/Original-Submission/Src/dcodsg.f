      SUBROUTINE ZZDCOD( MODE,  LL,    LINLEN, LNO,    EOPHAS, VERIFY,
     -                  OUTPT, KEYWRD, KWDICT, KWINFO, KWLENG, PWDICT,
     -                 PWINFO, PWLENG, PARSTR, PARVAL, NPARS, LOOPDL,
     -                 TRMODE, TRACE,  TRACUN, MAP,    UPCASE, *     )

C## A R G U M E N T S:

C-----ENTRY AND EXIT KEYWORD NUMBER AND MODES OF OPERATION.

      INTEGER           KEYWRD,   MODE, OUTPT, TRACUN
      LOGICAL           TRMODE, VERIFY, TRACE, UPCASE

C-----FORMAL EXTERNAL DICTIONARY PARAMETERS.

      INTEGER                   KWLENG,               PWLENG
      INTEGER              MAP( KWLENG )
      INTEGER           KWINFO( KWLENG, 2 ),  PWINFO( PWLENG, 2 )
      CHARACTER * (*)   KWDICT,  PWDICT

C-----THE LINE TO BE DECODED AND THE POSITION IN LINE OF THE LAST
C-----CHARACTER IN LINE TO BE CONSIDERED.


      INTEGER     NULLIN,   MXCRIT,    RGL
      PARAMETER ( NULLIN=0, MXCRIT=20, RGL=MXCRIT+1 )
      INTEGER     ERL,          CPL,          SVL
      PARAMETER ( ERL=MXCRIT+2, CPL=MXCRIT+3, SVL=MXCRIT+4 )
      INTEGER     NL
      PARAMETER ( NL = MXCRIT+4 )

      INTEGER           LINLEN, LNO
      CHARACTER * (*)   LL(*)
C-----                       EXIT VARIABLES.
      INTEGER           NPARS
      LOGICAL           EOPHAS
      CHARACTER * (*)   PARSTR ( * )
      CHARACTER*1       LOOPDL

      REAL              PARVAL ( * )
C!!!! DOUBLE PRECISION  PARVAL ( * )

C## S T A T U S:
C               IGNORE LINES BEGINNING WITH  "C!!!!" .
C
C               SYSTEM  DEPENDENCE:                      NONE.
C               SINGLE/DOUBLE CONVERSION: NEEDED (SEE CONVRT).
C
C               THIS VERSION IS IN   S I N G L E   PRECISION.
C!!!!           THIS VERSION IS IN   D O U B L E   PRECISION.
C
C               SINGLE/DOUBLE CONVERSION:        NOT REQUIRED.
C
C>RCS $HEADER: DCOD.F,V 2.2 91/11/20 10:52:39 BUCKLEY EXP $
C>RCS $LOG:     DCOD.F,V $
C>RCS REVISION 2.2  91/11/20  10:52:39  BUCKLEY
C>RCS FINAL SUBMISSION TO TOMS
C>RCS
C>RCS REVISION 2.1  90/08/08  17:23:12  BUCKLEY
C>RCS MINOR FIX TO DECLARATION.
C>RCS
C>RCS REVISION 2.0  90/07/31  11:23:17  BUCKLEY
C>RCS MINOR ERRORS CORRECTED, ESP START=NO
C>RCS
C>RCS REVISION 1.9  89/06/30  13:32:49  BUCKLEY
C>RCS PREPARING SUBMITTED VERSION OF MT
C>RCS
C>RCS REVISION 1.6  89/05/21  12:47:48  BUCKLEY
C>RCS FIXED DOUBLE LETTER
C>RCS
C>RCS REVISION 1.5  89/05/20  21:49:22  BUCKLEY
C>RCS FINAL DECODE I THINK
C>RCS
C>RCS REVISION 1.4  89/05/20  20:35:05  BUCKLEY
C>RCS TEMP
C>RCS
C>RCS REVISION 1.3  89/05/18  12:20:18  BUCKLEY
C>RCS FINAL TEST OF MT BEFORE SUBMITTING
C>RCS
C>RCS REVISION 1.2  89/05/15  14:47:09  BUCKLEY
C>RCS INITIAL INSTALLATION OF MT INTO RCS FORM.
C>RCS
C>RCS REVISION 1.1  89/01/17  16:50:07  BUCKLEY
C>RCS INITIAL REVISION
C>RCS

C## D E S C R I P T I O N:
C    THE PURPOSE OF ZZDCOD IS TO DECODE THE FIRST COMMAND FROM A
C    LINE OF CHARACTERS.  THE COMPLETE SYNTAX FOR A COMMAND IS
C    SHOWN IN EXTERNAL DOCUMENTATION, BUT MAY BE BRIEFLY
C    SUMMARIZED AS :
C      [KEYWORD]  [DELIMITER] [PARAMETER]  [SEPARATOR]  WHERE
C    THE ITEMS IN THE   ARE OPTIONAL AND MAY BE REPEATED AND
C    THE TERMINATING SEPARATOR IS REQUIRED IF ANOTHER COMMAND FOLLOWS.
C
C    INFORMATION CONCERNING THE LINE OF CHARACTERS IS PASSED INTO
C    ZZDCOD VIA THE ARGUMENTS :
C       LINE   - THE LINE OF CHARACTERS, OF UNKNOWN LENGTH
C       LINLEN - THE INDEX OF THE LAST NON-BLANK CHARACTER IN THE
C                LINE
C       LP - AN INDEX INTO THE LINE, SEE BELOW
C
C    AFTER A SUCCESSFUL DECODING OF A COMMAND, THE CHARACTERS
C    FROM THE COMMAND ARE REMOVED FROM LINE AND LINLEN IS
C    DECREASED.  IMPLEMENTATION NOTE : ALL CHARACTERS IN LINE ARE
C    PUT INTO UPPER CASE;  ON SYSTEMS WHICH USE TWO OR MORE
C    CHARACTER LOCATIONS TO REPRESENT A SINGLE LOWER CASE
C    CHARACTER THE LINLEN WILL DECREASE DUE TO THIS AS WELL.
C
C    PARAMETER IS A GENERAL CONCEPT AND INCLUDES THE FOLLOWING
C    FORMS :
C
C       PARAMETER TYPE     ACCEPTABLE VALUES
C       --------------     -----------------------------------
C           TLOGIC         NONE ARE EXPECTED
C           TYREAL         REAL (DOUBLE PRECISION) VALUES
C           TINTGR         INTEGER VALUES
C           TINTPW         INTEGER VALUES AND PARAMETER WORDS
C                          FROM THE PARAMETER WORD DICTIONARY
C                          DESCRIBED BELOW
C           TINTLT         INTEGER VALUES AND LITERAL STRINGS
C           TYCHAR         SINGLE CHARACTERS
C           TSTRNG         LITERAL STRINGS
C
C     INFORMATION CONCERNING VALID KEYWORDS AND PARAMETER WORDS
C     IS PASSED INTO ZZDCOD VIA THE FOLLOWING ARGUMENTS :
C
C        ARGUMENT     INTERPRETATION
C        --------     -------------------------------------------
C         KWDICT      AN ARRAY OF CHARACTER VARIABLES, EACH
C                     CONTAINING A KEYWORD
C         KWLENG      THE NUMBER OF SUCH KEYWORDS (THE DIMENSION)
C         KWINFO      AN INTEGER MATRIX WITH DIMENSIONS KWLENG
C                     ROWS BY 2 COLUMNS.  THE FIRST COLUMN
C                     CONTAINS AN INTEGER DICTATING THE MAXIMUM
C                     NUMBER OF PARAMETERS WHICH THE KEYWORD ON
C                     THE SAME ROW IN KWDICT MAY ACCEPT.  THE
C                     SECOND COLUMN CONTAINS ONE OF THE PARAMETER
C                     TYPES LISTED ABOVE AND DICTATES THE FORM
C                     OF ACCEPTABLE PARAMETERS.
C
C          PWDICT     AN ARRAY OF CHARACTER VARIABLES, EACH
C                     CONTAINING A PARAMETER WORD
C          PWLENG     THE NUMBER OF SUCH PARAMETER WORDS (THE
C                     DIMENSION)
C          PWINFO     AN INTEGER MATRIX WITH DIMENSIONS PWLENG
C                     ROWS BY 2 COLUMNS.  THE FIRST COLUMN
C                     CONTAINS THE INDEX IN KWDICT OF THE
C                     KEYWORD ASSOCIATED WITH THE PARAMETER WORD.
C                     THE SECOND COLUMN GIVES AN INTEGER VALUE FOR
C                     THE PARAMETER WORD.
C
C     THE ARGUMENTS MODE AND KEYWRD ARE USED FOR DISTINCT
C     PURPOSES ON ENTRY AND ON EXIT.
C
C            ! MODE = ! INTERPRETATION
C         ---+--------+-------------------------------------------
C            ! DCDNEW ! THE FIRST CALL TO ZZDCOD WITH THIS LINE
C            !        ! OF CHARACTERS [ LINE( 1:LINLEN ) ]
C            +--------+-------------------------------------------
C            ! DCDOLD ! A SUBSEQUENT CALL TO ZZDCOD WHERE NO
C          O !        ! ERRORS WITH ASSOCIATED REPLACEMENT
C          N !        ! REQUESTS, OR CONTINUATION REQUESTS WERE
C            !        ! ENCOUNTERED ON THE PREVIOUS CALL.  THE
C          E !        ! LINE OF CHARACTERS IS LINE( 1:LINLEN )
C          N +--------+-------------------------------------------
C          T ! DCDREP ! A CALL AFTER CHARACTERS CONTAINING AN
C          R !        ! ERROR IS REPLACED;  LINE( 1:LP )
C          Y !        ! CONTAINS THE NEW INFORMATION;  SEE DCDERR
C            +--------+-------------------------------------------
C            ! DCDCON ! A CALL IMMEDIATELY FOLLOWING A RETURN
C            !        ! FROM ZZDCOD WITH A CONTINUATION REQUEST.
C            !        ! ALL OF LINE( 1:LINLEN ) IS NEW.
C         ---+--------+-------------------------------------------
C            ! DCDDON ! A COMMAND HAS BEEN SUCCESSFULLY DECODED
C            !        ! AND THE LINE IS EXHAUSTED, I.E. LINLEN = 0
C            +--------+-------------------------------------------
C            ! DCDMOR ! A COMMAND HAS BEEN SUCCESSFULLY DECODED
C            !        ! AND THE LINE CONTAINS CHARACTERS TO BE
C          O !        ! DECODED, I.E. LINLEN > 0
C          N +--------+-------------------------------------------
C            ! DCDERR ! AN ERROR WAS FOUND IN THE LINE BETWEEN
C          E !        ! THE FIRST AND (LP)TH CHARACTER.
C          X !        ! THE CALLING PROGRAM SHOULD REPLACE THE
C          I !        ! INFORMATION.  NOTE THAT IT IS THE CALLING
C          T !        ! PROGRAM'S RESPONSIBILITY TO HANDLE ANY
C            !        ! CHARACTERS ON THE LINE WHICH MAY BE
C            !        ! LOST DURING THE REPLACEMENT.  ZZERRM WILL
C            !        ! BE CALLED BY ZZDCOD.  IF NOT TRMODE THEN
C            !        ! THE ALTERNATIVE RETURN WILL BE TAKEN.
C            +--------+-------------------------------------------
C            ! DCDBAD ! ON ENTRY MODE WAS INVALID.
C         ---+--------+-------------------------------------------
C
C
C             INTERPRETATION    ! KEYWRD !    INTERPRETATION
C                ON ENTRY       !        !       ON  EXIT
C         ----------------------+--------+----------------------
C          NO DEFAULT KEYWORD   !  < 0   ! ERROR WITH NUMBER
C                               !        ! EQUAL TO KEYWRD
C                               +--------+----------------------
C                               !  = 0   ! NO INPUT ON LINE
C         ----------------------+--------+----------------------
C          A DEFAULT KEYWORD    !  > 0   ! KEYWORD WITH INDEX
C          WITH INDEX INTO      !        ! INTO KWDICT EQUAL TO
C          KWDICT EQUAL TO      !        ! KEYWRD WAS FOUND.
C          KEYWRD IS ASSUMED    !        !
C         ----------------------+--------+----------------------
C
C     WHEN KEYWRD IS POSITIVE ON EXIT, THEN INFORMATION
C     CONCERNING THE DECODED PARAMETERS (IF ANY) IS PASSED OUT
C     USING THE FOLLOWING ARGUMENTS :
C
C        VARIABLE     INTERPRETATION
C        --------     -------------------------------------------
C         PARVAL      AN ARRAY OF REALS INTO WHICH REALS AND
C                     CONVERTED INTEGERS ARE PLACED.  INTEGERS
C                     MAY REPRESENT VALUES, KEYWORD OR PARAMETER
C                     WORD INDICES, OR THE CHARACTER SEQUENCE
C                     NUMBER OF A CHARACTER PARAMETER
C         PARSTR      AN ARRAY OF CHARACTER VARIABLES INTO WHICH
C                     LITERAL STRINGS ARE PLACED
C         NPARS        THE NUMBER OF PARAMETERS FOUND
C
C     NPARS SHOULD, ON ENTRY, BE THE DIMENSION OF PARVAL AND PARSTR.
C
C  PARAMETER TYPE        IN PARSTR                IN PARVAL
C  --------------        ---------                ---------
C
C  TLOGIC         THE STRING 'TRUE' OR F    NOTHING
C                  DEPENDS ON WHETHER KEYWORD
C                  PRECEDED WITH 'NO'.
C
C  TYREAL          NOTHING                         THE REAL VALUE
C
C  TYCHAR          THE CHARACTER, FOLLOWED BY     THE VALUE OF ICHAR,
C                  BLANKS.                        CONVERTED TO REAL.
C  TSTRNG          THE STRING, FOLLOWED BY BLANKS.  NOTHING
C
C  TINTGR  (I) WITH NO RADIX PRESENT
C
C                  BLANK                         CONVERTED VALUE
C
C          (II) RADIX  B, O, H OR X PRESENT
C
C                  BLANK                         INTEGER, CONVERTED
C                                                ACCORDING TO GIVEN
C                                                BASE, AS A REAL.
C          (III) OTHER RADIX
C
C                    , FOLLOWED BY RADIX          INTEGER IN BASE 10,
C                  CHARACTER SPECIFIED.          CONVERTED TO REAL.
C
C      EOPHAS IS TRUE WHEN (1) A SPECIAL END OF LINE
C     CHARACTER IS ENCOUNTERED, OR (2) WHEN THE ENTIRE LINE
C     IS EXHAUSTED (I.E. EXIT IN DCDDON);  TRMODE INDICATES
C     WHETHER THE SESSION IS INTERACTIVE OR NOT.
C
C     ==> FINISH
C
C## E N T R Y   P O I N T S:
C                          ZZDCOD     THE NATURAL ENTRY.
C                          ZZDSET     INITIALIZE CONTROL VARIABLES.
C
C## S U B R O U T I N E S:
C
C     ZZBASE     CONVERT STRING IN GIVEN BASE TO INTEGER.
C     ZZERRM     WRITE AN ERROR OR INFORMATIVE MESSAGE.
C     ZZSRCH     SEARCH A DICTIONARY.
C
C     ICHAR, INDEX, LEN,  REAL(DBLE)  ...INTRINSIC
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

      INTEGER     CTOUPP,     CTOLOW,     CTOCAP
      PARAMETER ( CTOUPP = 1, CTOLOW = 2, CTOCAP = 3 )
C-----CONTROL CONSTANTS

      INTEGER     SINGLE,     DOUBLE
      PARAMETER ( SINGLE = 1, DOUBLE = 2 )

      INTEGER   TLOGIC, TYREAL, TINTGR, TSTRNG
      INTEGER   TINTPW, TYCHAR, TINTLT, TYNONE

      PARAMETER ( TLOGIC = 0,             TINTGR = 3, TINTPW = 4 )
      PARAMETER ( TINTLT = 5, TYCHAR = 6, TSTRNG = 7, TYNONE = 8 )

      PARAMETER ( TYREAL = SINGLE )
C!!!! PARAMETER ( TYREAL = DOUBLE )

C  MODES FOR DECODING.

      INTEGER     DCDNEW,     DCDDON
      PARAMETER ( DCDNEW = 0, DCDDON = DCDNEW )

      INTEGER     DCDOLD,     DCDMOR
      PARAMETER ( DCDOLD = 1, DCDMOR = DCDOLD )

      INTEGER     DCDREP,     DCDERR
      PARAMETER ( DCDREP = 2, DCDERR = DCDREP )

      INTEGER     DCDCON,     DCDBAD
      PARAMETER ( DCDCON = 3, DCDBAD = 4 )

C-----SPECIAL CHARACTERS AND STRINGS, WHICH SHOULD BE STANDARD.

      CHARACTER * (*)   ALPHBT
      PARAMETER       ( ALPHBT = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ' )

      CHARACTER * (*)   DIGITS
      PARAMETER       ( DIGITS = '0123456789' )

      CHARACTER * (*)   DECPT                ,  PLUS
      PARAMETER       ( DECPT  = '.'         ,  PLUS = '+'    )

      CHARACTER * (*)   MINUS                ,  BLANK
      PARAMETER       ( MINUS  = '-'         ,  BLANK = ' '   )

C-----LENGTHS FOR  STRING ARRAYS.

      INTEGER     ERRLEN,      NERROR,      IMPERR
      PARAMETER ( ERRLEN = 40, NERROR = 23, IMPERR = 4 )

      INTEGER     TOTERR
      PARAMETER ( TOTERR = NERROR+IMPERR )

C-----DEFINE INTERNAL DECODING STATES.

      INTEGER      SBEGIN    ,  SBWORD    ,  SXPARM    ,  SBNUMR
      PARAMETER  ( SBEGIN = 1,  SBWORD = 2,  SXPARM = 3,  SBNUMR = 4 )

      INTEGER      SBPWRD    ,  SBCHAR    ,  SBSTR
      PARAMETER  ( SBPWRD = 5,  SBCHAR = 6,  SBSTR  = 7 )

      INTEGER      NOSTRT,      NOSIGN,      NOINTG,      NORADX
      PARAMETER (  NOSTRT = 0,  NOSIGN = 1,  NOINTG = 2,  NORADX = 3 )

      INTEGER      NORADC,      NODECM,      NOEXPO,      NOEXSG
      PARAMETER  ( NORADC = 4,  NODECM = 5,  NOEXPO = 6,  NOEXSG = 7 )

      INTEGER      NOEXNO,      NOHEXD,      NOERRR,      NOFRAC
      PARAMETER (  NOEXNO = 8,  NOHEXD = 9,  NOERRR = 10, NOFRAC = 11 )

C-----INTERNAL CODES FOR RESULT AFTER DECIDING.

      INTEGER     RESKEY,     RESNKY,     RESKTR
      PARAMETER ( RESKEY = 1, RESNKY = 2, RESKTR = 3 )

      INTEGER     RESPAR,     RESSTR,     RESBAD
      PARAMETER ( RESPAR = 4, RESSTR = 5, RESBAD = 6 )

C## L O C A L   D E C L:

C-----               FIRST TIME CALLED
      LOGICAL   FIRST

C-----                                 ERROR MESSAGES AND EXIT CODES.
      CHARACTER * (ERRLEN)  MESSAG(TOTERR)
      INTEGER               EXIT  (TOTERR)

C-----NOTE THAT THE FOLLOWING ARE DEFINED EXTERNALLY AND THE VALUES
C-----PASSED INTO ZZDCOD THROUGH THE ENTRY POINT ZZDSET. SEE THE
C-----SAVE VARIABLES.

C  ---SPECIAL CHARACTERS.

      INTEGER CONT, DEL, ESC, SEP, RADIX, ASSMT, COM1, COM2
      INTEGER RADPRB, RADFUN, RADGRP, STRNG1, STRNG2, LOOPA, LOOPG

      INTEGER  NDCCHS

      PARAMETER ( CONT   = 1, DEL    = 2, ESC    = 3, SEP    = 4,
     -            RADIX  = 5, ASSMT  = 6, COM1   = 7, COM2   = 8,
     -            RADPRB = 9, RADFUN =10, RADGRP =11, STRNG1 =12,
     -            STRNG2 =13, LOOPA  =14, LOOPG  =15            )

      PARAMETER ( NDCCHS = 15 )

      CHARACTER*(NDCCHS) DCC

      CHARACTER * 1    CONTIN, DELIMT, ESCAPC, SEPAR,   COMM1
      CHARACTER * 1     COMM2, RADX,  ALOOP,  GLOOP
      CHARACTER * 1     STR1,   STR2,  ASSGN

C-----RADIX CONVERSIONS.

      INTEGER    RADBAS(5)
      CHARACTER  VALSC*16,   RADCON*5

C-----LOCAL VARIABLES.

      REAL                     VALUE, RD
C!!!! DOUBLE PRECISION         VALUE, RD

      INTEGER  DEFKEY, DFPNUM, DFPTYP, DINDEX, I, PTSVL
      INTEGER  KTRIND, MXPARS, KEYPTR, MESSIN, OLDNIN, INLINE
      INTEGER  NEWKIN, NEWPIN, OLDKIN, OLDPIN, PARNUM, NEWNIN
      INTEGER  PARTYP, STATE,  SLEN, SMARKR,   I1, RESULT, LASTCP
      INTEGER  DIGF, DIGT, TRYF, TRYT, AFTDEC,  NOSTAT, B, LP
      INTEGER  OFFST1, OFFST2

      LOGICAL  CANINT, CANSTR, ESCAPE, GOTALF, GOTDEL, IGN
      LOGICAL  GOTINT, GOTREL, GOTSEP, HAVKTR, KNOWN,  MAYBEK, MAYBEP
      LOGICAL  SEPDEL, STILLK, STILLP,  TRAILD, COMENT, MAYBEN, STILLN
      LOGICAL  STRNGE, GOTRAD, HARD,   GOTNUM, EXCESS, GOTHEX, RADHEX
      LOGICAL  DEFNUM, SKIP,   HAVDEC, HAVRDX, HAVFRC, HAVHEX
      LOGICAL  INSTRG

      CHARACTER * ( 1)   NEXTCH, UPCH
      CHARACTER * ( 2)   NO
      CHARACTER * (20)   RFRMAT

C## S A V E:
C-----      ALL VARIABLES SAVED BECAUSE OF POSSIBLE CONTINUATION.
      SAVE

C## E Q U I V A L E N C E S: NONE ARE DEFINED.
C## C O M M O N:             NONE IS DEFINED.
C## D A T A:
C                                 1         2         3         4
C   MESSAGE LENGTH  =    1234567890123456789012345678901234567890

      DATA  MESSAG(24) / 'DCOD IMPL. ERR: BAD VALUE FOR RESULT' /
      DATA  MESSAG(25) / 'DCOD IMPL. ERR: NOT ENOUGH PARAMETERS.'/
      DATA  MESSAG(26) / 'DCOD IMPL. ERR: NUMBER STRINGS'    /
      DATA  MESSAG(27) / 'DCOD IMPL. ERR: STRING TOO SHORT'  /

      DATA  MESSAG(1) / 'INVALID INITIAL CHARACTER.'               /
      DATA  MESSAG(2) / 'DELIMITER NOT EXPECTED AFTER KEYWORD.'    /
      DATA  MESSAG(3) / 'MISSING AN EXPECTED PARAMETER.'           /
      DATA  MESSAG(4) / 'DELIMITER NOT EXPECTED AFTER PARAMETER.'  /
      DATA  MESSAG(5) / 'UNKNOWN KEYWORD BEFORE DELIMITER.'        /
      DATA  MESSAG(6) / 'UNRECOGNIZED WORD BEFORE SEPARATOR.'      /
      DATA  MESSAG(7) / 'BAD CHARACTER IN STRING.'                 /
      DATA  MESSAG(8) / 'PARAMETER MUST BE NUMERIC.'               /
      DATA  MESSAG(9) / 'PARAMETER MUST BEGIN WITH A LETTER.'      /
      DATA  MESSAG(10)/ 'THERE ARE TOO MANY PARAMETERS.'           /
      DATA  MESSAG(11)/ 'TOO MANY CHARACTERS FOR A NUMBER.'        /
      DATA  MESSAG(12)/ 'INVALID LETTER IN NUMBER.'                /
      DATA  MESSAG(13)/ 'INVALID SYNTAX FOR A REAL NUMBER.'        /
      DATA  MESSAG(14)/ 'INVALID SYNTAX FOR AN INTEGER.'           /
      DATA  MESSAG(15)/ 'INVALID CHARACTER FOR A NUMBER.'          /
      DATA  MESSAG(16)/ 'UNRECOGNIZED KEYWORD.'                    /
      DATA  MESSAG(17)/ 'INVALID SYNTAX IN PARAMETER.'             /
      DATA  MESSAG(18)/ 'PARAMETER MUST BE A SINGLE CHARACTER.'    /
      DATA  MESSAG(19)/ 'KEYWORD FOR SYNONYM IS UNKNOWN.'          /
      DATA  MESSAG(20)/ 'SYNONYM IS MISSING.'                      /
      DATA  MESSAG(21)/ 'INVALID CHARACTER IN KEYWORD.'            /
      DATA  MESSAG(22)/ 'PARAMETER NOT EXPECTED WITH SYNONYM.'     /
      DATA  MESSAG(23)/ 'INVALID POSITION FOR STRING.'             /

      DATA  FIRST /T/

C-----FORMAT STRINGS.

      DATA  RFRMAT / '( BN, F  40  .0 )'       /

      DATA  RADBAS / 2, 10, 8, 16, 16 /

      DATA   VALSC / '0123456789ABCDEF' /, RADCON / 'BDOHX' /

C## E X E C U T I O N
C## E X E C U T I O N

C-----                DEFINE FUNCTION STATEMENT.
      RD(I) = REAL(I)
C!!!! RD(I) = DBLE(I)

C-----                TAKE LENGTH OF INPUT LINE.
      INLINE = LEN(LL(LNO))

C-----                SET UP THE EXIT ARRAY
      IF ( FIRST ) THEN
         DO 50 I = 1, NERROR
            EXIT(I) = I
  50     CONTINUE

         TRYF = 10 * NERROR
         DO 60 I = 1, IMPERR
            EXIT(NERROR+I) =  TRYF + I
  60     CONTINUE

         FIRST = F
      ENDIF

C-----CHECK FOR CORRECT MODE.

      IF ( MODE .NE. DCDNEW .AND. MODE .NE. DCDOLD .AND.
     -     MODE .NE. DCDREP .AND. MODE .NE. DCDCON      ) THEN
C                                          => BAD VALUE FOR MODE.
         CALL ZZERRM ( RD(MODE), *91000,
     -                 'NT DECODE ENTRY WITH BAD MODE VALUE.'     )
         GOTO 91000

      ELSE IF ( MODE .EQ. DCDCON ) THEN
C                                    => PICK UP WHERE WE LEFT OFF.
         LP = 0
      ELSE
C                                    => INITIALIZE ALL VARIABLES.
         SLEN = LEN ( PARSTR (1) )
         MXPARS = NPARS
         DEFKEY = KEYWRD

         IF ( DEFKEY .GT. 0 ) THEN
            DFPNUM = KWINFO( DEFKEY, 1 )
            DFPTYP = KWINFO( DEFKEY, 2 )
            CANINT = DFPTYP .EQ. TINTGR  .OR.  DFPTYP .EQ. TINTPW
     -                                   .OR.  DFPTYP .EQ. TINTLT
            CANSTR = DFPTYP .EQ. TSTRNG  .OR.  DFPTYP .EQ. TINTLT
         ELSE
            DFPNUM = 0
            DFPTYP = TYNONE
            CANINT = F
            CANSTR = F
         ENDIF

         KEYWRD =  0
         NPARS  =  0
         LASTCP =  0
         STATE  =  SBEGIN
         NOSTAT =  NOSTRT
         MESSIN =  0
         ESCAPE = F
         EOPHAS = F
         COMENT = F
         INSTRG = F

         LP     =  0
         LL(ERL)=  BLANK
         LL(SVL)=  BLANK
         PTSVL  =  0
         SMARKR =  0

         KNOWN  = F
         TRAILD = F
         HAVDEC = F
         HAVFRC = F
         HAVRDX = F
         HAVHEX = F
         SKIP   = F
      ENDIF

C     THE FOLLOWING LOGICALS ARE USED THROUGHOUT :

C     GOTALF = NEXTCH IS ALPHABETIC
C     GOTDEL = NEXTCH IS DELIMITER OR ASSIGNMENT
C     GOTSEP = NEXTCH IS SEPARATOR
C     SEPDEL = GOTDEL OR GOTSEP
C     GOTINT = NEXTCH IS DIGIT  OR + OR - OR /
C     GOTREL = NEXTCH IS DIGIT  OR + OR - OR . OR D OR E.
C     DEFNUM = NEXTCH IS + OR - OR . OR /    ( UNAMBIGUOUS )
C     HARD   = NEXTCH IS DIGIT OR ALPHABETIC (   AMBIGUOUS )

C-----          TOP OF THE LOOP FOR CHARACTER BY CHARACTER PROCESSING.
  100 LP = LP + 1

      GOTALF = F
      GOTDEL = F
      GOTSEP = F
      GOTINT = F
      GOTREL = F
      STRNGE = F
      GOTRAD = F
      DEFNUM = F

      IF ( LP .GT. LINLEN ) THEN
C                         => SUPPLY AN ARTIFICIAL END-OF-LINE CHARACTER.
         IF ( COMENT ) THEN
            MODE = DCDCON
            GOTO 90000
         ELSE
            EOPHAS = T
            GOTSEP = T
            GOTO 150
         ENDIF
      ELSE
C                              => CONSIDER THE NEXT CHARACTER IN LINE.
         NEXTCH = LL(LNO) (LP:LP)
         UPCH   = NEXTCH
         IF ( .NOT. UPCASE ) THEN
            CALL ZZCASE(UPCH,CTOUPP)
         ENDIF
      ENDIF

      IF ( COMENT ) THEN
C                         => IGNORE IT.
         IF ( NEXTCH .EQ. COMM1 .OR. NEXTCH .EQ. COMM2 ) THEN
            COMENT = F
         ENDIF
         GOTO 100

      ELSEIF ( ESCAPE ) THEN
C                   => PREVIOUS CHAR WAS ESCAPE CHARACTER.  THIS
C                   => CHAR MUST NOT BE GIVEN ANY SPECIAL MEANING.
         GOTALF = T
         ESCAPE = F

      ELSEIF ( INSTRG ) THEN
C                            =>   DO NOTHING, BUT CHANGE STATE
         IF ( NEXTCH .EQ. STR1 .OR. NEXTCH .EQ. STR2 ) THEN
            INSTRG = .NOT. INSTRG
            GOTO 100
         ENDIF

      ELSEIF ( NEXTCH .EQ. BLANK ) THEN
C                       => IF NOT INSIDE QUOTES IGNORE A BLANK.
         IF ( STATE .NE. SBSTR ) THEN
             IF ( STATE .EQ. SBWORD ) THEN
                PTSVL = PTSVL + 1
                LL(SVL)(PTSVL:PTSVL) = BLANK
             ENDIF
             GOTO 100
         ENDIF

      ELSEIF (( NEXTCH .EQ. DELIMT ) .OR. ( NEXTCH .EQ. ASSGN )
     -   .OR. ( NEXTCH .EQ. ALOOP )  .OR. ( NEXTCH .EQ. GLOOP )) THEN
         GOTDEL = T
         LOOPDL = NEXTCH

      ELSEIF ( NEXTCH .EQ. SEPAR ) THEN
         GOTSEP = T

      ELSEIF ( NEXTCH .EQ. PLUS .OR. NEXTCH .EQ. MINUS ) THEN
         GOTREL = T
         GOTINT = T
         DEFNUM = T

         IF      ( NOSTAT .EQ. NOSTRT ) THEN
            NOSTAT = NOSIGN
         ELSE IF ( NOSTAT .EQ. NOEXPO ) THEN
            NOSTAT = NOEXSG
         ELSE
            NOSTAT = NOERRR
         ENDIF

      ELSEIF ( NEXTCH .EQ. DECPT ) THEN
         GOTREL = T
         DEFNUM = T
         HAVDEC = T

         IF ( NOSTAT .EQ. NOINTG .OR. NOSTAT .EQ. NOSTRT
     -                           .OR. NOSTAT .EQ. NOSIGN ) THEN
            NOSTAT = NODECM
         ELSE
            NOSTAT = NOERRR
         ENDIF

      ELSEIF ( NEXTCH .EQ. RADX ) THEN
         GOTINT = T
         DEFNUM = T
         HAVRDX = T

         IF ( NOSTAT .EQ. NOINTG  .OR. NOSTAT .EQ. NOHEXD ) THEN
            NOSTAT = NORADX
         ELSE
            NOSTAT = NOERRR
         ENDIF

      ELSEIF ( INDEX ( DIGITS, NEXTCH ) .NE. 0 ) THEN
         GOTINT = T
         GOTREL = T

         IF ( NOSTAT .EQ. NOSTRT .OR. NOSTAT .EQ. NOINTG
     -                           .OR. NOSTAT .EQ. NOSIGN ) THEN
            NOSTAT = NOINTG
         ELSE IF ( NOSTAT .EQ. NOHEXD ) THEN
            NOSTAT = NOHEXD
         ELSE IF ( NOSTAT .EQ. NODECM ) THEN
            NOSTAT = NOFRAC
            HAVFRC = T
         ELSE IF ( NOSTAT .EQ. NOEXSG .OR. NOSTAT .EQ. NOEXPO
     -                                .OR. NOSTAT .EQ. NOEXNO ) THEN
            NOSTAT = NOEXNO
         ELSE IF ( NOSTAT .EQ. NOFRAC ) THEN
C                 NOTHING CHANGES.
         ELSE
            NOSTAT = NOERRR
         ENDIF

      ELSEIF ( INDEX ( ALPHBT, UPCH ) .NE. 0 ) THEN
         GOTALF = T
         GOTHEX = INDEX(VALSC,UPCH) .NE. 0
         RADHEX = UPCH .EQ. 'X' .OR. UPCH .EQ. 'H'
         GOTREL = UPCH .EQ. 'D' .OR. UPCH .EQ. 'E'

         IF ( GOTREL .AND. (    NOSTAT .EQ. NODECM
     -                     .OR. NOSTAT .EQ. NOINTG
     -                     .OR. NOSTAT .EQ. NOFRAC )  ) THEN
            NOSTAT = NOEXPO
         ELSE IF ( GOTHEX  .AND.
     -            (NOSTAT .EQ. NOSTRT .OR. NOSTAT .EQ. NOINTG
     -        .OR. NOSTAT .EQ. NOHEXD .OR. NOSTAT .EQ. NOSIGN ) ) THEN
            NOSTAT = NOHEXD
            HAVHEX = T
            GOTNUM = T
         ELSE IF ( NOSTAT .EQ. NOINTG ) THEN
            NOSTAT = NORADC
         ELSE IF ( NOSTAT .EQ. NOHEXD .AND. RADHEX ) THEN
            NOSTAT = NORADC
         ELSE IF ( NOSTAT .EQ. NORADX ) THEN
            IF ( HAVHEX .AND. RADHEX ) THEN
               NOSTAT = NORADC
            ELSE IF ( .NOT. RADHEX ) THEN
               NOSTAT = NORADC
            ELSE
               NOSTAT = NOERRR
            ENDIF
         ELSE IF ( NOSTAT .NE. NOSTRT ) THEN
            NOSTAT = NOERRR
         ENDIF

      ELSEIF ( NEXTCH .EQ. ESCAPC ) THEN
         ESCAPE = T
         GO TO 100

      ELSEIF ( NEXTCH .EQ. COMM1 .OR. NEXTCH .EQ. COMM2 ) THEN
         COMENT = T
         GOTO 100

      ELSEIF ( NEXTCH .EQ. STR1 .OR. NEXTCH .EQ. STR2 ) THEN
         INSTRG = T
         GOTO 100

      ELSEIF ( NEXTCH .EQ. CONTIN ) THEN
C                                        => DETERMINE HOW TO EXIT.
         IF ( MESSIN .NE. 0 ) THEN
C                                  => IGNORE THE CONTINUATION.
            GOTSEP = T
         ELSE
C           => EXIT NOW.
            MODE = DCDCON
            GO TO 90000
         ENDIF
      ELSE
         STRNGE = T
      ENDIF

  150 SEPDEL = GOTSEP .OR. GOTDEL
      GOTNUM = GOTINT .OR. GOTREL
      HARD   = GOTALF .OR. ( GOTNUM .AND. .NOT. DEFNUM )

      IF ( .NOT. GOTSEP ) THEN
C                               ADD NEXTCH TO STRING BEING BUILT.
         EXCESS = F
C         CALL ZZBSTR ( STRING,SMARKR+1,1,SMARKR,NEXTCH,EXCESS,F )
C         TO FACILITATE THE INCLUSION OF BLANKS, DON'T USE ZZBSTR,
C         INSTEAD USE THE FOLLOWING 7 LINES :
         IF ( SMARKR+1 .LE. INLINE ) THEN
            LL(ERL)( SMARKR+1 : SMARKR+1 ) = NEXTCH
            SMARKR = SMARKR+1
         ELSE
            MESSIN = 27
            GOTO 82000
         ENDIF
      ENDIF

      IF (TRACE .AND. SMARKR .NE. 0) WRITE (TRACUN,*)
     -                           'STRING IS NOW==>',LL(ERL)(1:SMARKR)
      IF (TRACE) WRITE (TRACUN,*)
     -            'WITH: NEXTCH="',NEXTCH,'" ',
     -            'SEPDEL,GOTNUM,HARD,STATE,NOSTAT ',
     -             SEPDEL,GOTNUM,HARD,STATE,NOSTAT

C.......................................................................
C        => BRANCH TO CORRECT STATE.

  200 IF ( MESSIN .NE. 0 ) THEN
         GOTO 10000
      ELSE
         GOTO ( 1000, 2000, 3000, 4000, 5000, 6000, 7000 ) STATE
      ENDIF
C.......................................................................

C                         STATE #1 (SBEGIN) FIRST CHARACTER ENCOUNTERED.
 1000 IF ( .NOT. GOTSEP ) THEN
C                             => ACT ON THE CHARACTER.
         TRYF = SMARKR
         IF ( GOTALF ) THEN
C                          => GO DIRECTLY TO STATE SBWORD.
            STATE  =  SBWORD
            MAYBEK = T
            OLDKIN =  0
            HAVKTR = F
            MAYBEP = T
            OLDPIN =  0

         ELSEIF( ( GOTINT .AND. CANINT ) .OR.
     -           ( GOTREL .AND. DFPTYP .EQ. TYREAL ) )  THEN

C                       => PREPARE TO GO DIRECTLY TO STATE SBNUMR.
            STATE = SBNUMR
            DIGF  = SMARKR
            NPARS = 1
         ELSEIF ( DFPTYP .EQ. TYCHAR ) THEN
            STATE = SBCHAR
            NPARS = 1
         ELSEIF ( DFPTYP .EQ. TSTRNG ) THEN
            STATE = SBSTR
            PTSVL = 0
            LL(SVL) = BLANK
            NPARS = 1
         ELSE
C                        => BAD CHARACTER.
            MESSIN = 1
         ENDIF
      ELSE
C        => IGNORE MULTIPLE SEPARATORS.
         GOTO 10000
      ENDIF

      IF ( STATE .NE. SBWORD ) THEN

         KEYWRD = DEFKEY
         PARNUM = DFPNUM
         PARTYP = DFPTYP

         EXCESS = F
         IGN    = T

         OFFST1 = (KEYWRD-1) * SLEN +1
         OFFST2 = KEYWRD * SLEN
         CALL ZZBSTR( LL(ERL), 1, 0, SMARKR,
     -                KWDICT( OFFST1 : OFFST2 ), EXCESS, IGN )
         CALL ZZBSTR( LL(ERL), SMARKR, SMARKR-1, SMARKR,
     -                                    ASSGN, EXCESS, IGN )
         TRYF = SMARKR
         DIGF = SMARKR
      ENDIF
      GOTO 200

C.......................................................................
C     STATE #2 (SBWORD)  HERE WE ARE BUILDING KEYWORDS, BUT THE
C                        SITUATION MAY BE VERY AMBIGUOUS WHEN TRAILING
C         DIGITS ARE PRESENT, OR THE KEYWORD OR DELIMITER HAS BEEN
C         LEFT OUT. WE STAY IN THIS STATE FROM THE STARTING ALPHA-
C         BETIC UNTIL SOME CHARACTER CLARIFIES THE SITUATION.


 2000 IF ( SKIP ) THEN
         IF ( GOTNUM ) THEN
            DEFNUM = T
         ELSE IF ( GOTALF ) THEN
            GOTO 10000
         ENDIF
      ENDIF

      IF ( HARD ) THEN
C        => NEXTCH IS ALPHABETIC OR A DIGIT, BUT NOT + - . OR /
C        => THIS IS WHERE THE SITUATION IS STILL AMBIGUOUS AND
C        => NO RESOLUTION OF THE MEANING HAS YET BEEN MADE.

         PTSVL = PTSVL + 1
         LL(SVL)(PTSVL:PTSVL) = NEXTCH
         TRYT = SMARKR
         IF ( .NOT. TRAILD .AND. GOTNUM ) THEN
C           => START OF POSSIBLE TRAILING NUMBER.
            DIGF = SMARKR
            DIGT = SMARKR
            TRAILD = T
         ELSEIF ( TRAILD  ) THEN
C           => IF TRAILING DIGITS ARE ALREADY PRESENT ADD NEXTCH TO THE
C           => DIGIT STRING.  STILL NO RESOLUTION.

            IF ( NOSTAT .EQ. NOERRR ) THEN
               TRAILD = F
            ELSE
               DIGT = SMARKR
            ENDIF
         ELSE
C           => WHEN THERE ARE NO TRAILING DIGITS ADD NEXTCH TO THE
C           => STRING OF CHARACTERS TAKEN TO BE NONNUMERIC.
C           => ALREADY DONE, SO NOTHING EXTRA NEEDED.

         ENDIF
C              OF " IF TRAIL D..."

         IF ( MAYBEK ) THEN
C           => IF THE TRIAL STRING WITHOUT NEXTCH WAS IN THE
C           => KEYWORD DICTIONARY, SEE IF THE ADDITION OF NEXTCH
C           => HAS AFFECTED ITS PRESENCE.

            NEWKIN = OLDKIN

            CALL ZZSRCH(LL(ERL)(TRYF:TRYT),SLEN,KWDICT,KWLENG,SLEN,
     -                 NEWKIN,T,UPCASE)

            STILLK = NEWKIN .GT. 0

         ENDIF

         IF ( TRYT .GT. TRYF+1 ) THEN

            IF ( MAYBEN ) THEN

               NEWNIN = OLDNIN
               I      = TRYF+2

 2020          CALL ZZSRCH(LL(ERL)(I:TRYT),SLEN-2,KWDICT,KWLENG,SLEN,
     -                    NEWNIN, T,UPCASE)

               IF ( NEWNIN .NE. 0 ) THEN

                  IF ( KWINFO(NEWNIN,2) .EQ. TLOGIC ) THEN
                     STILLN = T
                  ELSE
                     NEWNIN = NEWNIN + 1
                     GOTO 2020
                  ENDIF

               ELSE
                  STILLN = F
               ENDIF

            ENDIF

         ELSE

            NO = LL(ERL)(TRYF:TRYF+1)
            CALL ZZCASE(NO,CTOUPP)

            IF ( TRYT .EQ. TRYF   .AND. NO(1:1) .EQ. 'N'  .OR.
     -           TRYT .EQ. TRYF+1 .AND. NO(1:2) .EQ. 'NO') THEN

               STILLN = T
               NEWNIN = 0

            ELSE

               STILLN = F

            ENDIF

         ENDIF

         IF ( MAYBEP ) THEN
C           => PERFORM A SIMILAR ANALYSIS IF THE TRIAL STRING WAS
C           => IN THE PARAMETER DICTIONARY.

            NEWPIN = OLDPIN

            IF ( PWLENG .NE. 0 ) THEN

               CALL ZZSRCH(LL(ERL)(TRYF:TRYT),SLEN,PWDICT,PWLENG,SLEN,
     -                    NEWPIN, T,UPCASE)

               STILLP = NEWPIN .GT. 0

            ELSE

               STILLP = F

            ENDIF

         ENDIF

         IF ( .NOT. ( STILLK .OR. STILLP .OR. STILLN ) ) THEN
C           => THE TRIAL STRING WITH THE ADDITION OF NEXTCH IS
C           => NEITHER A KEYWORD NOR A PARAMETER AND SO WE MUST BE
C           => ABLE TO MAKE A DECISION NOW.

            MAYBEK = F
            MAYBEP = F
            MAYBEN = F
C                                   THIS IS CASE #1.
            ASSIGN 2100 TO AFTDEC

            GOTO 81000
C             GOTO 100

         ELSE IF ( SMARKR .GT. SLEN ) THEN

C           => SET FLAG TO SKIP SUCCEEDING ALPHABETIC CHARACTERS.

C                                  (CASE # 2 WAS DROPPED.)
            SKIP = T

         ENDIF
C              FOR " NOT STILL..."

C        => UPDATE MAYBEK AND MAYBEP AND THE RELATED DICTIONARY
C        => INDEXES AS APPROPRIATE.  IF MAYBEK THEN CHECK THE
C        => KEYWORD DICTIONARY FOR AN ACCEPTABLE KEYWORD WHICH
C        => PERMITS INTEGER OR REAL PARAMETERS.

         MAYBEK = STILLK

         IF ( MAYBEK ) THEN

            OLDKIN = NEWKIN

            PARTYP = KWINFO( NEWKIN, 2 )

            IF ( PARTYP .EQ. TYREAL .OR. PARTYP .EQ. TINTGR
     -      .OR. PARTYP .EQ. TINTLT .OR. PARTYP .EQ. TINTPW) THEN
               HAVKTR = T
               KEYPTR = SMARKR
               KTRIND = NEWKIN
               NOSTAT = NOSTRT
               TRAILD = F
            ENDIF

         ELSE

            HAVKTR = F

         ENDIF

         MAYBEP = STILLP

         IF ( MAYBEP ) THEN
            OLDPIN = NEWPIN
         ENDIF

         MAYBEN = STILLN

         IF ( MAYBEN ) THEN
            OLDNIN = NEWNIN
         ENDIF

C.....=> IN ALL CASES FOLLOWING, THE SITUATION IS NO LONGER
C.....   AMBIGUOUS, ALTHOUGH IT MAY STILL NEED FINAL CLARIFICATION.

      ELSEIF ( DEFNUM ) THEN
C                              => DEFINITELY NUMERIC!.
         IF ( .NOT. KNOWN ) THEN

            TRAILD = T
            TRYT   = SMARKR

            MAYBEK = F
            MAYBEN = F
            MAYBEP = F
C                                  CASE #3.
            ASSIGN 2300 TO AFTDEC

            GOTO 81000

         ELSE

            STATE = SBNUMR
            DIGF  = SMARKR

            GOTO 200

         ENDIF

      ELSEIF ( SEPDEL ) THEN

C        => SEVERAL POSSIBILITIES STILL, SINCE NOT RESOLVED YET.
C        => IF KNOWN IS TRUE, KEYWORD HAS BEEN FOUND, SO LITTLE TO DO.

         IF ( KNOWN ) THEN

C           => WE HAVE ALREADY INTERPRETED THE STRING, TAKE THE
C           => APPROPRIATE ACTION. ONLY CASES: A LONG KEYWORD
C           => OR PARAMETER WORD.

            IF ( RESULT .EQ. RESKEY .OR. RESULT .EQ. RESNKY ) THEN

C              => STRING IS INTERPRETABLE AS A KEYWORD.

                  IF ( GOTSEP .AND. PARTYP .EQ. TLOGIC ) THEN

                     IF ( RESULT .EQ. RESKEY ) THEN
                        PARSTR(1) = 'TRUE'
                     ELSE
                        PARSTR(1) = 'FALSE'
                     ENDIF

                     ENDIF

            ELSEIF ( RESULT .EQ. RESPAR ) THEN
C              => LL(ERL) IS INTERPRETABLE AS A PARAMETER.
               IF ( GOTDEL ) THEN
                  IF ( PARNUM .GT. 1 ) THEN
                     STATE  = SXPARM
                  ELSE
                     MESSIN = 4
                  ENDIF
               ELSE
                  GO TO 10000
               ENDIF

            ELSE
C              => STRING IS NOT PROPER.  NOTE THAT WE WILL NOT
C              => ARRIVE HERE WHEN A LITERAL STRING IS PROPER,
C              => AS THIS POSSIBILITY IS HANDLED IMMEDIATELY
C              => UPON DISCOVERY THAT NEITHER OF STILLK OR
C              => STILLP IS TRUE.

               IF ( GOTDEL ) THEN
C                 => UNKNOWN KEYWORD.
                  MESSIN = 5
               ELSE
C                 => UNKNOWN KEYWORD WHICH DOES NOT PERMIT A PARAMETER
C                 => OR AN UNKNOWN PARAMETER WITHOUT ITS KEYWORD.
                  MESSIN = 6
               ENDIF

            ENDIF

         ELSE
C           NOW INTERPRET.
C                                 CASE #4.
            ASSIGN 2400 TO AFTDEC

            GOTO 81000

         ENDIF

      ELSE
C                        => MUST BE A SPECIAL CHARACTER.
         TRAILD = F

         MAYBEP = F
         MAYBEN = F

C                              CASE #5.
         ASSIGN 2500 TO AFTDEC

         GOTO 81000
      ENDIF

      GOTO 10000

C...............

C                 => THESE ARE THE RETURNS AFTER USING DECIDE.
C     CASE # 1.  RESULT MUST BE  RESKTR, RESSTR OR RESBAD.

 2100 IF ( RESULT .EQ. RESKTR ) THEN
         STATE = SBNUMR
      ELSE IF ( RESULT .EQ. RESSTR ) THEN
         PTSVL = PTSVL - 1
         STATE = SBSTR
      ELSE IF ( RESULT .EQ. RESBAD ) THEN
         MESSIN = 6
      ELSE
         MESSIN = -4
      ENDIF
      GOTO 200

C     CASE #3.  RESULT IS RESKTR OR RESBAD.

 2300 IF ( RESULT .EQ. RESKTR ) THEN
         STATE = SBNUMR
      ELSE IF ( RESULT .EQ. RESBAD ) THEN
         MESSIN = 16
      ELSE
         MESSIN = -4
      ENDIF
      GOTO 200

C     CASE #4.  ANY VALUE OF RESULT IS POSSIBLE.

 2400 IF      ( RESULT .EQ. RESKEY .OR. RESULT .EQ. RESNKY ) THEN
C        => NOTHING TO DO.
      ELSE IF ( RESULT .EQ. RESKTR ) THEN
         STATE = SBNUMR
         DIGT  = SMARKR - 1
         GOTO 200
      ELSE IF ( RESULT .EQ. RESPAR ) THEN
C        => NOTHING TO DO.
      ELSE IF ( RESULT .EQ. RESSTR ) THEN
         PTSVL = PTSVL - 1
         STATE = SBSTR
         GOTO 200
      ELSE IF ( RESULT .EQ. RESBAD ) THEN
         MESSIN = 5
      ENDIF
      GOTO 10000

C     CASE #5.  RESULT CAN BE RESSTR OR RESBAD.

 2500 IF ( RESULT .EQ. RESSTR ) THEN
         PTSVL = PTSVL - 1
         STATE = SBSTR
      ELSE IF ( RESULT .EQ. RESKEY ) THEN
         IF ( PARTYP .EQ. TYCHAR) THEN
            STATE = SBCHAR
         ELSE IF ( PARTYP .EQ. TSTRNG ) THEN
            STATE = SBSTR
            PTSVL = 0
            LL(SVL) = BLANK
         ELSE
            MESSIN = 7
         ENDIF
      ELSE IF ( RESULT .EQ. RESBAD ) THEN
         MESSIN = 21
      ELSE
         MESSIN = -4
      ENDIF
      GOTO 200

C.......................................................................

C-----STATE #3 (SXPARM) EXPECTING PARAMETER

C     => TEST TO SEE WHETHER ANOTHER PARAMETER IS VALID, AND IF SO, IF
C        THERE IS ROOM IN PARS FOR IT. THE TYPE AND KEYWORD ARE KNOWN.

 3000 IF ( NPARS  .GE. PARNUM ) THEN
C                                     => TOO MANY PARAMETERS.
         MESSIN = 10
         GOTO 10000
      ELSEIF ( NPARS .GE. MXPARS ) THEN
C                                      => IMPLEMENTATION ERROR.
         MESSIN = 25
         GOTO 10000
      ENDIF

      KNOWN = F
      TRYF  = SMARKR
      IF ( PARTYP .EQ. TYCHAR ) THEN
         STATE  =  SBCHAR
      ELSE IF ( PARTYP .EQ. TSTRNG ) THEN
         STATE = SBSTR
         PTSVL = 0
         LL(SVL) = BLANK
      ELSEIF ( GOTALF ) THEN
C        => IF AN (ALPHABETIC) PARAMETER WORD IS PROPER,
C        => GO TO STATE SBPWRD;  OTHERWISE, IF A LITERAL STRING IS
C        => PROPER, GO TO STATE SBSTR.

         IF ( PARTYP .EQ. TYREAL .OR. PARTYP .EQ. TINTGR ) THEN
            MESSIN = 8
         ELSE IF  ( PARTYP .EQ. TINTLT ) THEN
            STATE = SBSTR
            PTSVL = 0
            LL(SVL) = BLANK
         ELSE
            STATE  = SBPWRD
         ENDIF
      ELSEIF ( INSTRG ) THEN
         IF ( PARTYP .EQ. TSTRNG .OR. PARTYP .EQ. TINTLT ) THEN
            STATE = SBSTR
            PTSVL = 0
            LL(SVL) = BLANK
         ELSE
            MESSIN = 23
         ENDIF
      ELSEIF ( GOTNUM ) THEN
C        => IF A NUMBER IS PROPER, GO TO STATE SBNUMR.
         STATE  = SBNUMR
         DIGF   = SMARKR
      ELSEIF ( SEPDEL ) THEN
C        => NO SEPARATOR OR DELIMITER IS ACCEPTABLE NOW.
         MESSIN = 3
      ELSE
C        => AND NEITHER IS ANY OTHER CHARACTER.
         MESSIN = 7
      ENDIF

      IF ( MESSIN .EQ. 0 ) THEN
         NPARS = NPARS + 1
      ENDIF

      GOTO 200

C.......................................................................

C-----STATE #4 (SBNUMR) BUILDING NUMERIC PARAMETER. HERE THE KEYWORD AND
C              TYPE ARE BOTH KNOWN.

 4000 IF ( GOTNUM .OR. GOTALF  ) THEN

C        => CHARACTER IS PROPER FOR A NUMBER. CHECK SYNTAX.

         IF ( NOSTAT .EQ. NOERRR ) THEN
            MESSIN = 15
         ELSE
            DIGT = SMARKR
         ENDIF

      ELSEIF ( SEPDEL ) THEN
C                            => INTERPRET NUMBER AS A REAL OR INTEGER.
         PARSTR(NPARS) = BLANK

         IF ( HAVFRC .OR. NOSTAT .EQ. NOEXNO ) THEN

            I = DIGT - DIGF + 1
            WRITE ( RFRMAT(8:13), '(I6)'      )   I
            READ  ( LL(ERL)(DIGF:DIGT), RFRMAT ) VALUE
            IF ( PARTYP .NE. TYREAL .AND.
     -           RD( NINT(VALUE) ) .NE. VALUE ) THEN
               MESSIN = 15
            ENDIF
         ELSE
            IF ( NOSTAT .EQ. NOINTG ) THEN
               B = 10
            ELSE IF ( NOSTAT .EQ. NODECM ) THEN
               B    = 10
               DIGT = DIGT -1
            ELSE IF ( NOSTAT .EQ. NOHEXD ) THEN
               B = 16
            ELSE IF ( NOSTAT .EQ. NORADC ) THEN

               I = INDEX ( RADCON, LL(ERL)(DIGT:DIGT) )
               IF ( I .NE. 0 ) THEN
                  B = RADBAS(I)
               ELSE
                  PARSTR(1) = RADX // LL(ERL)(DIGT:DIGT)
                  DIGT = DIGT - 1
                  B    = 10
               ENDIF
               IF ( HAVRDX ) THEN
                  DIGT = DIGT - 1
               ENDIF
            ENDIF

            IF ( MESSIN .EQ. 0 ) THEN
              CALL ZZBASE(I1,B,LL(ERL)(DIGF:DIGT),*4200)
            ENDIF

            VALUE = RD(I1)
         ENDIF
         GOTO 4300

 4200    MESSIN = 15
 4300    CONTINUE

         IF ( MESSIN .EQ. 0 ) THEN
            PARVAL(NPARS) = VALUE
         ENDIF

      ELSE
         MESSIN = 15
      ENDIF
      GO TO 10000
C.......................................................................

C     STATE #5 (SBPWRD) BUILDING ALPHABETIC PARAMETER TO BE
C                       OBTAINED FROM PARAMETER DICTIONARY.

 5000 IF ( HARD ) THEN
C        => NEXTCH IS PROPER FOR AN ALPHABETIC PARAMETER.

         TRYT = SMARKR

      ELSEIF ( SEPDEL ) THEN

C        => SEARCH FOR THE STRING IN THE PARAMETER WORD DICTIONARY.
C        => NOTE THAT THE PARAMETER WORD MAY HAVE ALREADY BEEN
C        => DETERMINED FROM STATE SBEGIN IF VERY LONG.

         IF ( .NOT. KNOWN ) THEN

C           ===> CHECK THIS

            DINDEX = 0

 5100       CALL ZZSRCH ( LL(ERL)(TRYF:TRYT),SLEN,PWDICT,PWLENG,SLEN,
     -                 DINDEX,   T,UPCASE)

            IF ( DINDEX .NE. 0 ) THEN

               IF (  MAP(PWINFO( DINDEX, 1) ) .NE. KEYWRD ) THEN
                  DINDEX = DINDEX + 1
                  GO TO 5100
               ENDIF

            ENDIF

            IF ( DINDEX .NE. 0 ) THEN
C              => IT WAS FOUND IN A DICTIONARY, ACCEPT IT.

               PARVAL(NPARS) = RD ( PWINFO( DINDEX, 2 ) )
               PARSTR(NPARS) = BLANK

               CALL ZZSHFT ( LL(ERL), TRYT+1, TRYF, SMARKR )
               OFFST1 = (DINDEX-1) * SLEN +1
               OFFST2 = DINDEX * SLEN
               EXCESS = T
               IGN    = F
               CALL ZZBSTR( LL(ERL), TRYF, 0, SMARKR,
     -                      PWDICT( OFFST1 : OFFST2 ), EXCESS, IGN )
               IF (EXCESS) THEN
                  WRITE ( OUTPT, '(A)' ) ' INTERNAL DECODE STRING'//
     -                                    ' CAPACITY EXCEEDED'
               ENDIF

            ELSE
               MESSIN = 16
            ENDIF
         ENDIF

      ELSE
         MESSIN = 17
      ENDIF

      GO TO 10000

C.......................................................................

C-----STATE # 6 (SBCHAR)  A SINGLE CHARACTER PARAMETER.

 6000 IF ( .NOT. SEPDEL ) THEN

         IF ( NPARS .EQ. LASTCP ) THEN

            IF ( NPARS  .GE. PARNUM ) THEN
C              => THE USER HAS SPECIFIED TOO MANY PARAMETERS.
               MESSIN = 10
               GOTO 10000
            ELSEIF ( NPARS .GE. MXPARS ) THEN
C              => IMPLEMENTATION ERROR.
               MESSIN = -3
               GOTO 10000
            ENDIF

            NPARS = NPARS + 1

         ENDIF

         PARSTR(NPARS) =           NEXTCH
         PARVAL(NPARS) = RD( ICHAR(NEXTCH) )

         LASTCP = NPARS

      ENDIF

      GO TO 10000

C.......................................................................

C-----STATE # 7 (SBSTR) BUILD A LITERAL STRING.

 7000 IF ( .NOT. SEPDEL ) THEN
C                             => ADD NEXTCH TO THE STRING.
         TRYT = SMARKR
         PTSVL = PTSVL + 1
         LL(SVL)(PTSVL:PTSVL) = NEXTCH
      ELSE
         PARSTR(NPARS) = LL(ERL)(TRYF:TRYT)
      ENDIF
      GO TO 10000
C.......................................................................


C     END OF PROCESSING FOR THIS CHARACTER.

10000 IF ( GOTDEL ) THEN
         STATE  = SXPARM
         NOSTAT = NOSTRT
      ELSE IF ( GOTSEP ) THEN

         IF ( MESSIN .NE. 0 ) THEN
            IF ( LP .GT. LINLEN ) THEN
               MODE = DCDDON
            ELSE
               MODE   = DCDERR
               LINLEN = LP - 1
            ENDIF
            GOTO 82000
         ELSE
C             => EXIT AFTER SUCCESSFULLY ANALYZING A STRING.

            IF ( PARTYP .EQ. TLOGIC ) THEN
               IF ( RESULT .EQ. RESKEY ) THEN
                  PARSTR(1) = 'TRUE'
               ELSE
                  PARSTR(1) = 'FALSE'
               ENDIF
               NPARS = 1
            ENDIF

            IF ( LP .LT. LINLEN ) THEN
C                                           => DROP ANALYZED PORTION.
               CALL ZZSHFT ( LL(LNO), LP+1, 1, LINLEN )
               LINLEN = LINLEN - LP
               MODE   = DCDMOR
            ELSE
               LINLEN   = 0
               LL(LNO)  = BLANK
               MODE     = DCDDON
            ENDIF
            GOTO 90000
         ENDIF
      ENDIF
      GOTO 100

C>>>>>>>>>>>>>>>>>>>>> R E M O T E   B L O C K   1 <<<<<<<<<<<<<<<<<<<<<

C        => DECIDE ON THE INTERPRETATION OF THE STRING.

C     RESULT:  RESKEY.  A KEYWORD.
C              RESNKY.  A KEYWORD PRECEDED BY "NO".
C              RESKTR.  A KEYWORD WITH A TRAILING NUMBER.
C              RESPAR.  A PARAMETER WORD FROM A DICTIONARY.
C              RESSTR.  A LITERAL STRING.
C              RESBAD.  INVALID INTERPRETATION.

81000    KNOWN  = T
         IGN    = F

         IF (TRACE) WRITE (TRACUN,*)
     -                     'MAYBEK,MAYBEN,MAYBEP,HAVKTR ',
     -                      MAYBEK,MAYBEN,MAYBEP,HAVKTR

         IF ( MAYBEP ) THEN
            I = MAP( PWINFO( OLDPIN, 1 ))

            IF ( DEFKEY .GT. 0  .AND.  DFPTYP .EQ. TINTPW
     -                          .AND.  I .NE. DEFKEY ) THEN
C              => SEARCH THE PARAMETER DICTIONARY.

               NEWPIN = OLDPIN

81100          NEWPIN = NEWPIN + 1
               CALL ZZSRCH(LL(ERL)(TRYF:TRYT),SLEN,PWDICT,PWLENG,SLEN,
     -                    NEWPIN, T,UPCASE )

               IF ( NEWPIN .GT. 0 ) THEN
                I = MAP( PWINFO( NEWPIN, 1 ))

                  NEWKIN = MAP( PWINFO( NEWPIN, 1 ) )
                  IF ( NEWKIN .NE. DEFKEY ) THEN
C                                            => TRY AGAIN.
                     GOTO 81100

                  ELSE
C                    => WE HAVE FOUND A SUITABLE PARAMETER WORD SO
C                    => WE MUST USE THIS PARAMETER FOR THE DEFAULT.

                     OLDKIN = DEFKEY
                     OLDPIN = NEWPIN

                  ENDIF
               ELSE
C                 => NOTHING BETTER FOUND SO
C                 => RETAIN ORIGINAL INTERPRETATION.
               ENDIF

            ELSE IF ( I .EQ. DEFKEY ) THEN
               OLDKIN = DEFKEY
            ENDIF
C               FOR IF FOUND ANOTHER.
         ENDIF

         IF ( MAYBEP .AND. OLDKIN .EQ. DEFKEY ) THEN
C           => ACCEPT A PARAMETER.  GIVE PREFERENCE TO A PARAMETER
C           => WHICH FITS STRING AND HAS DEFKEY AS ITS KEYWORD.

            KEYWRD    = OLDKIN
            NPARS     = 1
            PARVAL(1) = RD ( PWINFO(OLDPIN,2) )

            RESULT = RESPAR

            OFFST1 = (KEYWRD-1) * SLEN +1
            OFFST2 = KEYWRD * SLEN
            CALL ZZBSTR( LL(ERL), 1, 0, SMARKR,
     -                   KWDICT( OFFST1 : OFFST2 ), EXCESS, IGN )
            CALL ZZBSTR( LL(ERL), SLEN+1, 0, SMARKR,
     -                                         ASSGN, EXCESS, IGN )

            TRYF = TRYF + SLEN +1
            TRYT = TRYT + SLEN +1

            CALL ZZSHFT ( LL(ERL), TRYT+1, TRYF, SMARKR )
            OFFST1 = (OLDPIN-1) * SLEN +1
            OFFST2 = OLDPIN * SLEN
            CALL ZZBSTR( LL(ERL), TRYF, 0, SMARKR,
     -                   PWDICT( OFFST1 : OFFST2 ), EXCESS, IGN )

         ELSE IF ( MAYBEK .AND. MAYBEN ) THEN
C                               => ACCEPT THE KEYWORD FOUND
C                                  WITHOUT THE "NO".
            KEYWRD = OLDKIN
            RESULT = RESKEY

            CALL ZZSHFT ( LL(ERL), TRYT+1, TRYF, SMARKR )
            OFFST1 = (KEYWRD-1) * SLEN +1
            OFFST2 = KEYWRD * SLEN
            CALL ZZBSTR( LL(ERL), TRYF, 0, SMARKR,
     -                   KWDICT( OFFST1 : OFFST2 ), EXCESS, IGN )

         ELSE IF ( MAYBEK ) THEN
C                               KEYWORD.
            KEYWRD = OLDKIN
            RESULT = RESKEY

            CALL ZZSHFT ( LL(ERL), TRYT+1, TRYF, SMARKR )
            OFFST1 = (KEYWRD-1) * SLEN +1
            OFFST2 = KEYWRD * SLEN
            CALL ZZBSTR( LL(ERL), TRYF, 0, SMARKR,
     -                   KWDICT( OFFST1 : OFFST2 ), EXCESS, IGN )

         ELSE IF ( MAYBEN .AND. OLDNIN .NE. 0 ) THEN
C                                                   "NO" KEYWORD.
            KEYWRD = OLDNIN
            RESULT = RESNKY

            CALL ZZSHFT ( LL(ERL), TRYT+1, TRYF+2, SMARKR )
            OFFST1 = (KEYWRD-1) * SLEN +1
            OFFST2 = KEYWRD * SLEN
            CALL ZZBSTR( LL(ERL), TRYF+2, 0, SMARKR,
     -                   KWDICT( OFFST1 : OFFST2 ), EXCESS, IGN )

         ELSE IF ( HAVKTR .AND. TRAILD .AND. TRYT .GT. KEYPTR ) THEN

            KEYWRD = KTRIND
            DIGF   = KEYPTR + 1
            NPARS  = 1
            RESULT = RESKTR

            CALL ZZSHFT ( LL(ERL), KEYPTR+1, TRYF, SMARKR )
            OFFST1 = (KEYWRD-1) * SLEN +1
            OFFST2 = KEYWRD * SLEN
            CALL ZZBSTR( LL(ERL), TRYF, 0, SMARKR,
     -                   KWDICT( OFFST1 : OFFST2 ), EXCESS, IGN )
            CALL ZZBSTR( LL(ERL), SMARKR, 0, SMARKR, ASSGN, EXCESS,IGN)

            DIGF = DIGF + SLEN - TRYT +2
            DIGT = DIGT + SLEN - TRYT +2

         ELSE IF ( MAYBEP ) THEN

C           => ACCEPT A PARAMETER.  GIVE PREFERENCE TO A PARAMETER
C           => WHICH FITS STRING AND HAS DEFKEY AS ITS KEYWORD.

            KEYWRD    = MAP(PWINFO(OLDPIN,1))
            NPARS     = 1
            PARVAL(1) = RD ( PWINFO(OLDPIN,2) )

            RESULT = RESPAR

            OFFST1 = (KEYWRD-1) * SLEN +1
            OFFST2 = KEYWRD * SLEN
            CALL ZZBSTR( LL(ERL), 1, 0, SMARKR,
     -                   KWDICT( OFFST1 : OFFST2 ), EXCESS, IGN )
            CALL ZZBSTR( LL(ERL), SLEN+1, 0, SMARKR,
     -                                         ASSGN, EXCESS, IGN )

            TRYF = TRYF + SLEN +1
            TRYT = TRYT + SLEN +1

            CALL ZZSHFT ( LL(ERL), TRYT+1, TRYF, SMARKR )
            OFFST1 = (OLDPIN-1) * SLEN +1
            OFFST2 = OLDPIN * SLEN
            CALL ZZBSTR( LL(ERL), TRYF, 0, SMARKR,
     -                   PWDICT( OFFST1 : OFFST2 ), EXCESS, IGN )

         ELSE
C             => TAKE IT AS A LITERAL STRING, IF POSSIBLE.

            IF ( CANSTR ) THEN
               KEYWRD = DEFKEY
               NPARS  = 1
               RESULT = RESSTR
            ELSE
               RESULT = RESBAD
            ENDIF

         ENDIF

         IF ( RESULT .NE. RESBAD ) THEN
            PARNUM = KWINFO(KEYWRD,1)
            PARTYP = KWINFO(KEYWRD,2)
         ENDIF

         IF (TRACE) WRITE (TRACUN,*) 'RESULT,KEYWRD,PARNUM,PARTYP ',
     -                                RESULT,KEYWRD,PARNUM,PARTYP

         GOTO AFTDEC

C## E R R O R S:

C                                 EXIT WITH AN ERROR CONDITION.
82000 MODE   = DCDERR
      KEYWRD = -EXIT( MESSIN )

      IF ( KEYWRD .LT. -NERROR ) THEN
C                                    => ALWAYS A FATAL ERROR.
         CALL ZZERRM ( VALUE, *91000,
     -                 'NT DECODE FATAL ERROR : ' // MESSAG( MESSIN ) )
      ELSE
C                                    => LEAVE IT TO USER TO DECIDE.
         CALL ZZERRM ( VALUE, *91000,
     -                 'NT DECODE ERROR : ' // MESSAG( MESSIN )  )
      ENDIF
      GO TO 90000

C## E N T R Y  ZZZZZZ:
                                    ENTRY  ZZDSET ( DCC )
         CONTIN = DCC(CONT  :CONT  )
         DELIMT = DCC(DEL   :DEL   )
         ESCAPC = DCC(ESC   :ESC   )
         SEPAR  = DCC(SEP   :SEP   )
         RADX   = DCC(RADIX :RADIX )
         COMM1  = DCC(COM1  :COM1  )
         COMM2  = DCC(COM2  :COM2  )
         GLOOP  = DCC(LOOPG :LOOPG )
         ALOOP  = DCC(LOOPA :LOOPA )
         STR1   = DCC(STRNG1:STRNG1)
         STR2   = DCC(STRNG2:STRNG2)
         ASSGN  = DCC(ASSMT :ASSMT )
      RETURN

C## E X I T
90000      CONTINUE
      IF ( MODE .NE. DCDCON ) THEN
         IF (TRACE) WRITE (TRACUN,*) 'MODE, KEYWRD, NPARS==> ',
     -                                MODE,KEYWRD,NPARS
         IF (TRACE) WRITE (TRACUN,*) 'STRINGS==> ',
     -                               (PARSTR(I),I=1,NPARS)
         IF (TRACE) WRITE (TRACUN,*) 'VALUES ==> ',
     -                               (PARVAL(I),I=1,NPARS)
         IF ( VERIFY .AND. SMARKR .NE. 0 ) THEN
            I = 1
90500       IF ( I+INLINE .LT. SMARKR ) THEN
                WRITE ( OUTPT, '(A,A)' ) ' ',LL(ERL)(I:I+INLINE-1)
                I = I + INLINE
                GOTO 90500
            ENDIF
            WRITE ( OUTPT, '(A,A)' ) ' ', LL(ERL)(I:SMARKR)
         ENDIF
      ENDIF
      RETURN

C-----RETURN ON ERROR.
91000                 RETURN 1

C## F O R M A T S:  NONE ARE DEFINED.
C##                 E N D         OF ZZDCOD.
                    END
