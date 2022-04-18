C     PROGRAM  GO
C
C----------------------------------------------------------------------
C
C     FAMILY
C     ------
C     SYSTEM/USER INTERFACE
C
C     PURPOSE
C     -------
C     THIS IS A SAMPLE MAIN PROGRAM TO CALL THE
C     DRIVING ROUTINE OF THE MACRO PROCESSOR.
C
C----------------------------------------------------------------------
C
C         GLOBAL CONSTANTS
C
      CHARACTER*1         CA,      CBLANK,  CC,      CI,      CLEFT,
     A                    CMINUS,  CPLUS,   CPOINT,  CQUOTE,  CRIGHT,
     B                    CZ,      C0,      C9
      COMMON  / GLCOMC /  CA,      CBLANK,  CC,      CI,      CLEFT,
     A                    CMINUS,  CPLUS,   CPOINT,  CQUOTE,  CRIGHT,
     B                    CZ,      C0,      C9
C
C         INPUT / OUTPUT CONTROL INTERFACE
C
      CHARACTER*1         CBUFFR(2000)
      LOGICAL             LBREAK,  LFORT,   LISTI,   LISTO
      INTEGER             ICBADD,  ICBEND,  ICBEOL,  ICBSUB,  ICB0,
     A                    ICB1,    ICB2,    ICB3,    ICBDIM,  ICPLI,
     B                    ICPLO,   ILCTR,   ILNMBR,  ILPP,    IPAGE,
     A                    IUNITE,  IUNITI,  IUNITL,  IUNITO
      COMMON  / IOCOMC /  CBUFFR
      COMMON  / IOCOMI /  ICBADD,  ICBEND,  ICBEOL,  ICBSUB,  ICB0,
     A                    ICB1,    ICB2,    ICB3,    ICBDIM,  ICPLI,
     B                    ICPLO,   ILCTR,   ILNMBR,  ILPP,    IPAGE,
     A                    IUNITE,  IUNITI,  IUNITL,  IUNITO
      COMMON  / IOCOML /  LBREAK,  LFORT,   LISTI,   LISTO
C
C         MEMORY MANAGER INTERFACE
C
      CHARACTER*1         CSTORE(20000)
      INTEGER             IHASH(601),       ISTORE(6000)
      INTEGER             ICSDIM,  ICSP1,   ICSP2,   IHADIM,
     A                    ISFREE,  ISTDIM,  IS2HDC,  IS2HDS
      COMMON  / MMCOMC /  CSTORE
      COMMON  / MMCOMH /  IHASH
      COMMON  / MMCOMS /  ISTORE
      COMMON  / MMCOMI /  ICSDIM,  ICSP1,   ICSP2,   IHADIM,
     A                    ISFREE,  ISTDIM,  IS2HDC,  IS2HDS
C
C         MACRO PROCESSOR INTERFACE
C
      CHARACTER*1         CDIV,    CEOL,     CEOR,    CONC,    CSUB,
     A                    CTOP
      LOGICAL             LEMPTY,  LSUB
      COMMON  / MPCOMC /  CDIV,    CEOL,    CEOR,    CONC,    CSUB,
     A                    CTOP
      COMMON  / MPCOML /  LEMPTY,  LSUB
C
C         TEMPLATE PROCESSOR INTERFACE
C
      CHARACTER*1         CDIR,    CSTAR
      INTEGER             ICBP1(4),         ICBP2(4)
      INTEGER             ITOPDO,  IARGS, INESTD, INESTF
      LOGICAL             LCOL1,   LDIRL,   LEND,    LINITM,  L1TRIP
      COMMON  / TPCOMC /  CDIR,    CSTAR
      COMMON  / TPCOMI /  ICBP1,   ITOPDO,  IARGS,   ICBP2,   INESTD,
     B                    INESTF
      COMMON  / TPCOML /  LCOL1,   LDIRL,   LEND,    LINITM,  L1TRIP
      EXTERNAL TPDRV, TPMMIN
C
C   SET DIMENSIONS FOR ARRAYS
C
      ICBDIM = 2000
      ICSDIM = 20000
      IHADIM = 601
      ISTDIM = 6000
C
C   INITIALIZE TEMPLATE PROCESSOR
C
      CALL TPMMIN
C
C   CALL DRIVER
C   USING UNIX STANDARD ERROR, INPUT, AND OUTPUT UNITS
C
      CALL  TPDRV  (0, 5, 0, 6)
C
      STOP 0
      END
      SUBROUTINE  TPDRV  (IUE0, IUI0, IUL0, IUO0)
C
C-----------------------------------------------------------------------
C
C     FAMILY
C     ------
C     TEMPLATE PROCESSOR
C
C     PURPOSE
C     -------
C     THIS IS THE DRIVING ROUTINE OF THE TEMPLATE PROCESSOR.
C     IT CALLS ROUTINES TO READ, EVALUATE, AND WRITE LINES
C     UNTIL AN END DIRECTIVE IS ENCOUNTERED
C
C     PARAMETERS
C     ----------
C     IUE0    -I-  UNIT NUMBER FOR THE ERROR FILE
C     IUI0    -I-  UNIT NUMBER FOR THE INPUT FILE
C     IUL0    -I-  UNIT NUMBER FOR THE LISTING FILE
C     IUO0    -I-  UNIT NUMBER FOR THE OUTPUT FILE
C
C     COMMON VARIABLES AND DATA STRUCTURES
C     ------------------------------------
C     THE COMMENTS BELOW GIVE A BRIEF DESCRIPTION OF THE COMMON
C     VARIABLES USED BY THE ROUTINES OF THE TEMPLATE PROCESSOR.
C     A MORE DETAILED LOOK AT THE MAIN DATA STRUCTURES IS ALSO
C     INCLUDED.
C
C         GLOBAL CONSTANTS
C
C     COMMON  / GLCOMC /
C         CA       -    'A'           CPOINT   -    '.'
C         CBLANK   -    ' '           CQUOTE   -    '''
C         CC       -    'C'           CRIGHT   -    '('
C         CI       -    'I'           CZ       -    'Z'
C         CLEFT    -    '('           C0       -    '0'
C         CMINUS   -    '-'           C9       -    '9'
C         CPLUS    -    '+'
C
C         INPUT / OUTPUT CONTROL INTERFACE
C
C     COMMON  / IOCOMC /
C         CBUFFR   -    I/O BUFFER
C     COMMON  / IOCOMI /
C         ICBADD   -    NUMBER OF SPACES TO SKIP BEFORE THE CONTINUATION
C                       OF A BROKEN LINE
C         ICBEND   -    BUFFER POSITION OF END OF CURRENT LOGICAL LINE
C                       (LOGICAL LINE MAY INCLUDE SEVERAL ACTUAL LINES)
C         ICBEOL   -    BUFFER POSITION OF CURRENT EOL.
C         ICBSUB   -    BUFFER POSITION OF CURRENT SUB. PREF. CHARACTER
C         ICB0     -    BUFFER POSITION OF START OF CURRENT LINE
C         ICB1     -    BUFFER POSITION WHERE CURRENT PROCESSING BEGINS
C         ICB2     -    BUFFER POSITION WHERE CURRENT PROCESSING ENDS
C         ICB3     -    BUFFER POSITION OF END OF CURRENT LINE
C         ICBDIM   -    DIMENSION OF CBUFFR
C         ICPLI    -    INPUT LINE LENGTH
C         ICPLO    -    OUPUT LINE LENGTH
C         ILCTR    -    LINE NUMBER ON CURRENT LISTING PAGE
C         ILNMBR   -    LINE NUMBER FOR LISTING (OVER ALL PAGES)
C         ILPP     -    MAX NUMBER OF LINES PER LISTING PAGE
C         IPAGE    -    PAGE NUMBER ON LISTING
C         IUNITE   -    ERROR OUTPUT UNIT
C         IUNITI   -    INPUT UNIT
C         IUNITL   -    LISTING OUTPUT UNIT
C         IUNITO   -    STANDARD OUTPUT UNIT
C     COMMON  / IOCOML /
C         LBREAK   -    BREAK LONG LINES AT NICE PLACE IF TRUE
C         LFORT    -    USE FORTRAN CONTINUATION CHAR. IF TRUE
C         LISTI    -    LIST INPUT IF TRUE
C         LISTO    -    LIST OUTPUT IF TRUE
C
C         MEMORY MANAGER INTERFACE
C
C     COMMON  / MMCOMC /
C         CSTORE   -    CHARACTER STORAGE
C     COMMON  / MMCOMH /
C         IHASH    -    HASH TABLE (IHASH(I) IS AN INDEX INTO ISTORE)
C     COMMON  / MMCOMS /
C         ISTORE   -    INTEGER STORAGE; HOLDS THE POINTERS WHICH
C                       IMPLEMENT THE SYMBOL TABLE AND THE STACK
C     COMMON  / MMCOMI /
C         ICSDIM   -    DIMENSION OF ICSDIM
C         ICSP1    -    PTR. TO TOP CHARACTER IN SUBSTITUTION STACK
C         ICSP2    -    PTR. TO LAST CHAR. IN FIRST STRING ON STACK
C         IHADIM   -    DIMENSION OF IHASH
C         ISFREE   -    PTR. TO HEAD OF ISTORE FREELIST
C         ISTDIM   -    DIMENSION OF ISTORE
C         IS2HDC   -    PTR. TO HEAD OF FREE CHARACTER STORAGE BLOCKS
C                       (ACTUALLY AN INDEX INTO ISTORE)
C         IS2HDS   -    PTR. TO TOP OF STACK
C                       (ACTUALLY AN INDEX INTO ISTORE)
C
C         MACRO PROCESSOR INTERFACE
C
C     COMMON  / MPCOMC /
C         CDIV     -    '/'
C         CEOL     -    '-'
C         CEOR     -    '/'
C         CONC     -    '+'
C         CSUB     -    DOLLAR SIGN
C         CTOP     -    TOP CHAR. IN STACK
C     COMMON  / MPCOML /
C         LEMPTY   -    TRUE IF SUBSTITUTION STACK EMPTY
C         LSUB     -    TRUE IF SUBSTITUTIONS ARE TO BE PERFORMED
C
C         TEMPLATE PROCESSOR INTERFACE
C
C     COMMON  / TPCOMC /
C         CDIR     -    '*'
C         CSTAR    -    '*'
C     COMMON  / TPCOMI /
C         ICBP1    -    ICBP1(I) IS BUFF. POSITION OF START OF
C                       ITH ARGUMENT
C         ITOPDO   -    PTR. TO 'TOP' (INNERMOST) DO LOOP ENTRY
C                       IN ISTORE
C         IARGS    -    NUMBER OF ARGUMENTS IN A DIRECTIVE
C         ICBP2    -    ICBP2(I) IS BUFF. POSITION OF END OF
C                       ITH ARGUMENT
C         INESTD   -    DO LOOP NESTING DEPTH
C         INESTF   -    IF-ELSE-ENDIF NESTING DEPTH
C     COMMON  / TPCOML /
C         LCOL1    -    TRUE IF DIRECTIVES MUST BEGIN IN COL 1
C         LDIRL    -    TRUE IF A DIRECTIVE HAS BEEN FOUND
C         LEND     -    TRUE IF AN END DIRECTIVE HAS BEEN FOUND
C         LINITM   -    TRUE IF MMINIT HAS BEEN CALLED
C         L1TRIP   -    TRUE IF ONE TRIP DO-LOOPS SHOULD BE ASSUMED
C
C
C     DATA STRUCTURES
C     ---------------
C
C     I/O BUFFER
C       THE ARRAY CBUFFR HOLDS THE I/O BUFFER.  INPUT LINES ARE READ
C       IN, MACRO SUBSTITUTIONS PERFORMED, AND LISTING  AND OUTPUT
C       (WHEN APPROPRIATE) ARE DONE FROM THE I/O BUFFER.
C
C     INTEGER STORAGE
C       THE ARRAY ISTORE IS USED TO HOLD THE POINTERS WHICH IMPLEMENT
C       THE SYMBOL TABLE AND  THE SUBSTITUTION STACK. IT IS USED IN
C       BLOCKS OF 3 ELEMENTS AT A TIME.  THE VARIABLE ISFREE POINTS
C       TO THE HEAD OF A LINKED LIST OF FREE ISTORE BLOCKS. INITIALLY
C       ALL BLOCKS ARE FREE (THE 3RD ELEMENT IN A BLOCK POINTS TO THE
C       NEXT FREE BLOCK).
C
C     CHARACTER STORAGE
C       THE ARRAY CSTORE PROVIDES A POOL OF CHARACTER STORAGE.  IT
C       IS USED TO RECORD MACRO NAMES AND VALUES, AS WELL AS STRINGS
C       WHICH MUST BE PUSHED ONTO THE SUBSTITUTION STACK.  THE VARIABLE
C       IS2HDC POINTS TO THE HEAD OF A FREELIST OF CHARACTER STORAGE
C       BLOCKS.  THIS FREELIST IS MADE UP OF ISTORE BLOCKS OF THE
C       FOLLOWING FORMAT:
C               ISTORE(I)  = CSTORE INDEX OF FIRST CHAR. IN BLOCK
C               ISTORE(I+1)= CSTORE INDEX OF LAST CHAR. IN BLOCK
C               ISTORE(I+2)= POINTER TO NEXT BLOCK
C
C     SYMBOL TABLE
C       THE SYMBOL TABLE KEEPS TRACK OF MACRO NAMES AND VALUES.  IT
C       IS BUILT OUT OF ISTORE BLOCKS WHICH CONTAIN POINTERS TO
C       OTHER ISTORE BLOCKS OR INDEXES INTO CSTORE.  GIVEN A MACRO
C       NAME, ROUTINE MMHASH COMPUTES ITS HASH INDEX IH.  THEN
C       IHASH(IH) IS THE ISTORE INDEX OF THE SYMBOL TABLE ENTRY FOR
C       THAT NAME.  IF IHASH(IH)=I SAY, THE ISTORE BLOCK AT I HOLDS
C       THE FOLLOWING:
C               ISTORE(I)   = PTR. TO ISTORE BLOCK FOR VARIABLE NAME
C               ISTORE(I+1) = PTR. TO HEAD OF LINKED LIST OF ISTORE
C                             BLOCKS FOR VALUE OF VARIABLE
C               ISTORE(I+2) = PTR. TO TAIL OF THE LINKED LIST FOR THE
C                             VALUE
C
C       AN ISTORE BLOCK FOR THE NAME OF A VARIABLE CONTAINS:
C               ISTORE(J)   = CSTORE INDEX OF FIRST CHAR. IN NAME
C               ISTORE(J+1) = CSTORE INDEX OF LAST CHAR. IN NAME
C               ISTORE(J+2) = 0
C
C       AN ISTORE BLOCK IN THE LINKED LIST WHICH KEEPS TRACK OF
C       THE VALUE OF A VARIABLE LOOKS LIKE:
C               ISTORE(K)   = CSTORE INDEX OF FIRST CHAR. ASSOCIATED
C                             WITH THIS BLOCK
C               ISTORE(K+1) = CSTORE INDEX OF LAST CHAR. ASSOCIATED
C                             WITH THIS BLOCK
C               ISTORE(K+2) = ISTORE INDEX OF NEXT BLOCK IN LIST
C                             (0 IF LAST ONE)
C
C     SUBSTITUTION STACK
C       WHEN A MACRO SUBSTITUTION IS FOUND, IT AND THE REST OF THE
C       CURRENT LINE ARE PUSHED ONTO THE SUBSTITUTION STACK.  THE
C       MACRO NAME IS POPPED OFF AND REPLACED BY ITS VALUE.  CHARACTERS
C       ARE THEN POPPED OFF THE STACK, INTO THE I/O BUFFER, UNTIL
C       THE STACK IS EMPTY OR ANOTHER SUBSTITUTION IS CALLED FOR.
C       IF ANOTHER MACRO SUBSTITUTION IS NEEDED THE SAME PROCESS IS
C       REPEATED--THE MACRO NAME IS REPLACED BY ITS VALUE, AND THE
C       STACK POPPING RESUMES.
C
C       THE STACK IS IMPLEMENTED AS A LINKED LIST OF ISTORE BLOCKS.
C       THE VARIABLE IS2HDS POINTS TO THE TOP BLOCK ON THE STACK.
C       A BLOCK AT INDEX I CONTAINS:
C               ISTORE(I)   = PTR. TO ISTORE BLOCK WHICH POINTS TO A
C                             STRING ON THE STACK
C               ISTORE(I+1) = CSTORE INDEX OF 1ST CHAR. OF
C                             CORRESPONDING STRING
C               ISTORE(I+2) = LINK TO NEXT ISTORE BLOCK ON STACK
C                             (0 IF THERE IS NONE)
C
C       THE FORMAT OF AN ISTORE BLOCK WHICH POINTS TO A STRING ON THE
C       STACK IS LIKE THAT OF ONE WHICH POINTS TO A VARIABLE NAME:
C               ISTORE(J)   = CSTORE INDEX OF FIRST CHAR. IN STRING
C               ISTORE(J+1) = CSTORE INDEX OF LAST CHAR. IN STRING
C               ISTORE(J+2) = 0
C
C
C
C-----------------------------------------------------------------------
C
C         GLOBAL CONSTANTS
C
      CHARACTER*1         CA,      CBLANK,  CC,      CI,      CLEFT,
     A                    CMINUS,  CPLUS,   CPOINT,  CQUOTE,  CRIGHT,
     B                    CZ,      C0,      C9
      COMMON  / GLCOMC /  CA,      CBLANK,  CC,      CI,      CLEFT,
     A                    CMINUS,  CPLUS,   CPOINT,  CQUOTE,  CRIGHT,
     B                    CZ,      C0,      C9
C
C         INPUT / OUTPUT CONTROL INTERFACE
C
      CHARACTER*1         CBUFFR(2000)
      LOGICAL             LBREAK,  LFORT,   LISTI,   LISTO
      INTEGER             ICBADD,  ICBEND,  ICBEOL,  ICBSUB,  ICB0,
     A                    ICB1,    ICB2,    ICB3,    ICBDIM,  ICPLI,
     B                    ICPLO,   ILCTR,   ILNMBR,  ILPP,    IPAGE,
     A                    IUNITE,  IUNITI,  IUNITL,  IUNITO
      COMMON  / IOCOMC /  CBUFFR
      COMMON  / IOCOMI /  ICBADD,  ICBEND,  ICBEOL,  ICBSUB,  ICB0,
     A                    ICB1,    ICB2,    ICB3,    ICBDIM,  ICPLI,
     B                    ICPLO,   ILCTR,   ILNMBR,  ILPP,    IPAGE,
     A                    IUNITE,  IUNITI,  IUNITL,  IUNITO
      COMMON  / IOCOML /  LBREAK,  LFORT,   LISTI,   LISTO
C
C         MEMORY MANAGER INTERFACE
C
      CHARACTER*1         CSTORE(20000)
      INTEGER             IHASH(601),       ISTORE(6000)
      INTEGER             ICSDIM,  ICSP1,   ICSP2,   IHADIM,
     A                    ISFREE,  ISTDIM,  IS2HDC,  IS2HDS
      COMMON  / MMCOMC /  CSTORE
      COMMON  / MMCOMH /  IHASH
      COMMON  / MMCOMS /  ISTORE
      COMMON  / MMCOMI /  ICSDIM,  ICSP1,   ICSP2,   IHADIM,
     A                    ISFREE,  ISTDIM,  IS2HDC,  IS2HDS
C
C         MACRO PROCESSOR INTERFACE
C
      CHARACTER*1         CDIV,    CEOL,     CEOR,    CONC,    CSUB,
     A                    CTOP
      LOGICAL             LEMPTY,  LSUB
      COMMON  / MPCOMC /  CDIV,    CEOL,    CEOR,    CONC,    CSUB,
     A                    CTOP
      COMMON  / MPCOML /  LEMPTY,  LSUB
C
C         TEMPLATE PROCESSOR INTERFACE
C
      CHARACTER*1         CDIR,    CSTAR
      INTEGER             ICBP1(4),         ICBP2(4)
      INTEGER             ITOPDO,  IARGS, INESTD, INESTF
      LOGICAL             LCOL1,   LDIRL,   LEND,    LINITM,  L1TRIP
      COMMON  / TPCOMC /  CDIR,    CSTAR
      COMMON  / TPCOMI /  ICBP1,   ITOPDO,  IARGS,   ICBP2,   INESTD,
     B                    INESTF
      COMMON  / TPCOML /  LCOL1,   LDIRL,   LEND,    LINITM,  L1TRIP
      INTEGER IUE0, IUI0, IUL0, IUO0
      EXTERNAL TPINIT, MPLINE, TPEVAL, IOWRIT
C
      CALL  TPINIT  (IUE0, IUI0, IUL0, IUO0)
C
   10 CONTINUE
          ICBEOL  =  0
          CALL  MPLINE  (.TRUE.)
          CALL  TPEVAL
          IF  (.NOT. LDIRL)  CALL  IOWRIT
      IF  (.NOT. LEND)  GO  TO  10
C
      RETURN
      END
      SUBROUTINE  IOERRM  (LFATAL, CFMT)
C
C----------------------------------------------------------------------
C
C     FAMILY
C     ------
C     INPUT/OUTPUT
C
C     PURPOSE
C     -------
C     TO PRINT OUT THE OFFENDING LINE AND AN ERROR MESSAGE BENEATH IT.
C     IF THE ERROR IS FATAL, PROCESSOR EXECUTION IS TERMINATED.
C
C     PARAMETERS
C     ----------
C     LFATAL  -I-  TRUE FOR FATAL ERRORS
C     CFMT    -I-  FORMAT FOR ERROR MESSAGE
C
C----------------------------------------------------------------------
C
C         INPUT / OUTPUT CONTROL INTERFACE
C
      CHARACTER*1         CBUFFR(2000)
      LOGICAL             LBREAK,  LFORT,   LISTI,   LISTO
      INTEGER             ICBADD,  ICBEND,  ICBEOL,  ICBSUB,  ICB0,
     A                    ICB1,    ICB2,    ICB3,    ICBDIM,  ICPLI,
     B                    ICPLO,   ILCTR,   ILNMBR,  ILPP,    IPAGE,
     A                    IUNITE,  IUNITI,  IUNITL,  IUNITO
      COMMON  / IOCOMC /  CBUFFR
      COMMON  / IOCOMI /  ICBADD,  ICBEND,  ICBEOL,  ICBSUB,  ICB0,
     A                    ICB1,    ICB2,    ICB3,    ICBDIM,  ICPLI,
     B                    ICPLO,   ILCTR,   ILNMBR,  ILPP,    IPAGE,
     A                    IUNITE,  IUNITI,  IUNITL,  IUNITO
      COMMON  / IOCOML /  LBREAK,  LFORT,   LISTI,   LISTO
C
C         LOCAL VARIABLES AND PARAMETERS
C
      CHARACTER*(*)      CFMT
      LOGICAL             LFATAL
      INTEGER             I
      EXTERNAL            IOPAGE
C
      IF  (IUNITE .EQ. IUNITL)  CALL  IOPAGE  (2)
      IF  (ICB0 .GT. ICB2)  GO  TO 10
      WRITE  (IUNITE, 1010)  (CBUFFR(I), I=ICB0,ICB2)
 1010 FORMAT(' ********   ', 117A1)
   10 CONTINUE
      WRITE  (IUNITE, CFMT)
      IF  (LFATAL)  STOP 1
C
      RETURN
      END
      SUBROUTINE  IOLIST  (LNUMBR)
C
C----------------------------------------------------------------------
C
C     INPUT/OUTPUT
C
C     PURPOSE
C     -------
C     TO LIST THE LINE CURRENTLY IN THE INPUT/OUTPUT BUFFER.
C
C     PARAMETER
C     ---------
C     LNUMBR  -I-  TRUE IF THE LINE SHOULD BE NUMBERED
C
C----------------------------------------------------------------------
C
C         INPUT / OUTPUT CONTROL INTERFACE
C
      CHARACTER*1         CBUFFR(2000)
      LOGICAL             LBREAK,  LFORT,   LISTI,   LISTO
      INTEGER             ICBADD,  ICBEND,  ICBEOL,  ICBSUB,  ICB0,
     A                    ICB1,    ICB2,    ICB3,    ICBDIM,  ICPLI,
     B                    ICPLO,   ILCTR,   ILNMBR,  ILPP,    IPAGE,
     A                    IUNITE,  IUNITI,  IUNITL,  IUNITO
      COMMON  / IOCOMC /  CBUFFR
      COMMON  / IOCOMI /  ICBADD,  ICBEND,  ICBEOL,  ICBSUB,  ICB0,
     A                    ICB1,    ICB2,    ICB3,    ICBDIM,  ICPLI,
     B                    ICPLO,   ILCTR,   ILNMBR,  ILPP,    IPAGE,
     A                    IUNITE,  IUNITI,  IUNITL,  IUNITO
      COMMON  / IOCOML /  LBREAK,  LFORT,   LISTI,   LISTO
C
C         LOCAL VARIABLES AND PARAMETERS
C
      LOGICAL             LNUMBR
      INTEGER             I
      EXTERNAL            IOPAGE
C
      CALL  IOPAGE  (1)
      IF  (.NOT. LNUMBR)  GO  TO  20
      ILNMBR  =  ILNMBR + 1
C
      IF  (ICB1 .LE. ICB2)  GO  TO  10
      WRITE  (IUNITL, 1010)  ILNMBR
      GO  TO  999
C
   10 CONTINUE
      WRITE  (IUNITL, 1020)  ILNMBR, (CBUFFR(I), I=ICB1,ICB2)
      GO  TO  999
C
   20 CONTINUE
      IF  (ICB1 .LE. ICB2)  GO  TO  30
      WRITE  (IUNITL, 1030)
      GO  TO  999
C
   30 CONTINUE
      WRITE  (IUNITL, 1040)  (CBUFFR(I), I=ICB1,ICB2)
C
  999 CONTINUE
      RETURN
 1010 FORMAT(' ', I8)
 1020 FORMAT(' ', I8, 3X, 117A1)
 1030 FORMAT(' ')
 1040 FORMAT(' ', 11X, 117A1)
      END
      SUBROUTINE  IOPAGE  (IL)
C
C----------------------------------------------------------------------
C
C     FAMILY
C     ------
C     INPUT/OUTPUT
C
C     PURPOSE
C     -------
C     TO DETERMINE IF THERE IS ROOM TO PRINT THE SPECIFIED NUMBER
C     OF LINES ON THE CURRENT PAGE. IF THERE IS NOT, A NEW PAGE
C     IS BEGUN AND A HEADING IS PRINTED.
C
C     PARAMETERS
C     ----------
C     IL      -I-  NUMBER OF LINES TO BE PRINTED
C
C----------------------------------------------------------------------
      INTEGER IL
C
C         INPUT / OUTPUT CONTROL INTERFACE
C
      CHARACTER*1         CBUFFR(2000)
      LOGICAL             LBREAK,  LFORT,   LISTI,   LISTO
      INTEGER             ICBADD,  ICBEND,  ICBEOL,  ICBSUB,  ICB0,
     A                    ICB1,    ICB2,    ICB3,    ICBDIM,  ICPLI,
     B                    ICPLO,   ILCTR,   ILNMBR,  ILPP,    IPAGE,
     A                    IUNITE,  IUNITI,  IUNITL,  IUNITO
      COMMON  / IOCOMC /  CBUFFR
      COMMON  / IOCOMI /  ICBADD,  ICBEND,  ICBEOL,  ICBSUB,  ICB0,
     A                    ICB1,    ICB2,    ICB3,    ICBDIM,  ICPLI,
     B                    ICPLO,   ILCTR,   ILNMBR,  ILPP,    IPAGE,
     A                    IUNITE,  IUNITI,  IUNITL,  IUNITO
      COMMON  / IOCOML /  LBREAK,  LFORT,   LISTI,   LISTO
C
      ILCTR  =  ILCTR + IL
      IF  (ILCTR .LE. ILPP)  GO  TO  999
      IPAGE  =  IPAGE + 1
      ILCTR  =  3 + IL
      WRITE (IUNITL,1010) IPAGE
C
  999 CONTINUE
      RETURN
 1010 FORMAT('1', 'PURDUE  UNIVERSITY  TEMPLATE  PROCESSOR  ',
     A            '(V2 - 07/31/83)  PAGE', I6 //)
      END
      SUBROUTINE  IORDLN (CLINE, ICL1, ICL2, IUNIT)
C
C----------------------------------------------------------------------
C
C     FAMILY
C     ------
C     INPUT/OUTPUT
C
C     PURPOSE
C     -------
C     TO READ A LINE INTO THE INPUT/OUTPUT BUFFER. THIS
C     MAY BE REPLACED BY A MORE EFFICIENT LOCAL I/O ROUTINE
C
C     PARAMETERS
C     ----------
C     CLINE   -I-  I/O BUFFER
C     ICL1    -I-  INDEX OF THE FIRST CHARACTER TO BE READ
C     ICL2    -I-  INDEX OF THE LAST CHARACTER TO BE READ
C     IUNIT   -I-  INPUT UNIT NUMBER
C
C----------------------------------------------------------------------
C
      INTEGER ICL1, ICL2, IUNIT
      CHARACTER*(*)      CLINE(ICL2)
      INTEGER             I, IBOT
C
C  ACCESS CDIR DIRECTIVE PREFIX
C
C
C         TEMPLATE PROCESSOR INTERFACE
C
      CHARACTER*1         CDIR,    CSTAR
      INTEGER             ICBP1(4),         ICBP2(4)
      INTEGER             ITOPDO,  IARGS, INESTD, INESTF
      LOGICAL             LCOL1,   LDIRL,   LEND,    LINITM,  L1TRIP
      COMMON  / TPCOMC /  CDIR,    CSTAR
      COMMON  / TPCOMI /  ICBP1,   ITOPDO,  IARGS,   ICBP2,   INESTD,
     B                    INESTF
      COMMON  / TPCOML /  LCOL1,   LDIRL,   LEND,    LINITM,  L1TRIP
C
      CHARACTER*1         STREND(5)
      SAVE STREND
      DATA STREND(1)/'*'/,STREND(2)/'E'/,STREND(3)/'N'/
      DATA STREND(4)/'D'/,STREND(5)/' '/
C
      READ (IUNIT, 1010, END=999) (CLINE(I), I=ICL1,ICL2)
      RETURN
 999  CONTINUE
      STREND(1) = CDIR
      DO 10 I=1,4
      CLINE(ICL1+I-1)=STREND(I)
 10   CONTINUE
      IBOT=ICL1+4
      DO 20 I=IBOT,ICL2
      CLINE(I)=STREND(5)
 20   CONTINUE
 1010 FORMAT(132A1)
      END
      SUBROUTINE  IOREAD
C
C----------------------------------------------------------------------
C
C     FAMILY
C     ------
C     SUBSTITUTION PROCESSING
C
C     PURPOSE
C     -------
C     TO FILL THE BUFFER WITH A LINE, REMOVE THE TRAILING BLANKS,
C     SET THE BUFFER POINTERS, AND APPEND AN END-OF-LINE MARKER.
C
C----------------------------------------------------------------------
C
C         GLOBAL CONSTANTS
C
      CHARACTER*1         CA,      CBLANK,  CC,      CI,      CLEFT,
     A                    CMINUS,  CPLUS,   CPOINT,  CQUOTE,  CRIGHT,
     B                    CZ,      C0,      C9
      COMMON  / GLCOMC /  CA,      CBLANK,  CC,      CI,      CLEFT,
     A                    CMINUS,  CPLUS,   CPOINT,  CQUOTE,  CRIGHT,
     B                    CZ,      C0,      C9
C
C         INPUT / OUTPUT CONTROL INTERFACE
C
      CHARACTER*1         CBUFFR(2000)
      LOGICAL             LBREAK,  LFORT,   LISTI,   LISTO
      INTEGER             ICBADD,  ICBEND,  ICBEOL,  ICBSUB,  ICB0,
     A                    ICB1,    ICB2,    ICB3,    ICBDIM,  ICPLI,
     B                    ICPLO,   ILCTR,   ILNMBR,  ILPP,    IPAGE,
     A                    IUNITE,  IUNITI,  IUNITL,  IUNITO
      COMMON  / IOCOMC /  CBUFFR
      COMMON  / IOCOMI /  ICBADD,  ICBEND,  ICBEOL,  ICBSUB,  ICB0,
     A                    ICB1,    ICB2,    ICB3,    ICBDIM,  ICPLI,
     B                    ICPLO,   ILCTR,   ILNMBR,  ILPP,    IPAGE,
     A                    IUNITE,  IUNITI,  IUNITL,  IUNITO
      COMMON  / IOCOML /  LBREAK,  LFORT,   LISTI,   LISTO
C
C         MEMORY MANAGER INTERFACE
C
      CHARACTER*1         CSTORE(20000)
      INTEGER             IHASH(601),       ISTORE(6000)
      INTEGER             ICSDIM,  ICSP1,   ICSP2,   IHADIM,
     A                    ISFREE,  ISTDIM,  IS2HDC,  IS2HDS
      COMMON  / MMCOMC /  CSTORE
      COMMON  / MMCOMH /  IHASH
      COMMON  / MMCOMS /  ISTORE
      COMMON  / MMCOMI /  ICSDIM,  ICSP1,   ICSP2,   IHADIM,
     A                    ISFREE,  ISTDIM,  IS2HDC,  IS2HDS
C
C         MACRO PROCESSOR INTERFACE
C
      CHARACTER*1         CDIV,    CEOL,     CEOR,    CONC,    CSUB,
     A                    CTOP
      LOGICAL             LEMPTY,  LSUB
      COMMON  / MPCOMC /  CDIV,    CEOL,    CEOR,    CONC,    CSUB,
     A                    CTOP
      COMMON  / MPCOML /  LEMPTY,  LSUB
      EXTERNAL   IORDLN, IOLIST, IOERRM
C
C         IF THERE IS ENOUGH SPACE IN THE BUFFER
C         READ A LINE FROM THE INPUT FILE
C
      ICB1  =  ICB2 + 1
      ICB2  =  ICB2 + ICPLI
      IF  (ICB2+2 .GT. ICBDIM)  GO  TO  30
      CALL  IORDLN (CBUFFR, ICB1, ICB2, IUNITI)
      IF  (LISTI)  CALL  IOLIST  (.TRUE.)
C
C         REMOVE TRAILING BLANKS
C
   10 CONTINUE
          IF  (CBUFFR(ICB2) .NE. CBLANK)  GO  TO  20
          ICB2  =  ICB2 - 1
      IF  (ICB2 .GE. ICB1)  GO  TO  10
C
C         ADD THE END-OF-LINE MARKER
C
   20 CONTINUE
      CBUFFR(ICB2+1)  =  CSUB
      CBUFFR(ICB2+2)  =  CEOL
      ICB3            =  ICB2
      ICBEOL          =  ICB2 + 2
      ICBEND          =  ICBEOL
      GO  TO  999
C
   30 CONTINUE
      CALL  IOERRM  (.TRUE.,
     A '('' ********   IOREAD - BUFFER SPACE EXCEEDED'')')
C
  999 CONTINUE
      RETURN
      END
      SUBROUTINE  IOWRIT
C
C----------------------------------------------------------------------
C
C     FAMILY
C     ------
C     SUBSTITUTION PROCESSING
C
C     PURPOSE
C     -------
C     TO WRITE THE LINE CURRENTLY IN THE BUFFER TO THE OUTPUT FILE.
C     IF THE -BREAK- OPTION IS SPECIFIED, AN ATTEMPT WILL BE MADE TO
C     BREAK LONG LINES AT A BLANK, RIGHT PARENTHESIS, COMMA, OR AN
C     ARITHMETIC OPERATOR. IF THE -FORTRAN- OPTION IS SPECIFIED,
C     CONTINUATION LINES WILL BE WRITTEN WITH CONTINUATION CHARACTERS
C     IN COLUMN SIX UNLESS THE LINE IS A COMMENT.
C
C----------------------------------------------------------------------
C
C         GLOBAL CONSTANTS
C
      CHARACTER*1         CA,      CBLANK,  CC,      CI,      CLEFT,
     A                    CMINUS,  CPLUS,   CPOINT,  CQUOTE,  CRIGHT,
     B                    CZ,      C0,      C9
      COMMON  / GLCOMC /  CA,      CBLANK,  CC,      CI,      CLEFT,
     A                    CMINUS,  CPLUS,   CPOINT,  CQUOTE,  CRIGHT,
     B                    CZ,      C0,      C9
C
C         INPUT / OUTPUT CONTROL INTERFACE
C
      CHARACTER*1         CBUFFR(2000)
      LOGICAL             LBREAK,  LFORT,   LISTI,   LISTO
      INTEGER             ICBADD,  ICBEND,  ICBEOL,  ICBSUB,  ICB0,
     A                    ICB1,    ICB2,    ICB3,    ICBDIM,  ICPLI,
     B                    ICPLO,   ILCTR,   ILNMBR,  ILPP,    IPAGE,
     A                    IUNITE,  IUNITI,  IUNITL,  IUNITO
      COMMON  / IOCOMC /  CBUFFR
      COMMON  / IOCOMI /  ICBADD,  ICBEND,  ICBEOL,  ICBSUB,  ICB0,
     A                    ICB1,    ICB2,    ICB3,    ICBDIM,  ICPLI,
     B                    ICPLO,   ILCTR,   ILNMBR,  ILPP,    IPAGE,
     A                    IUNITE,  IUNITI,  IUNITL,  IUNITO
      COMMON  / IOCOML /  LBREAK,  LFORT,   LISTI,   LISTO
C
C         MEMORY MANAGER INTERFACE
C
      CHARACTER*1         CSTORE(20000)
      INTEGER             IHASH(601),       ISTORE(6000)
      INTEGER             ICSDIM,  ICSP1,   ICSP2,   IHADIM,
     A                    ISFREE,  ISTDIM,  IS2HDC,  IS2HDS
      COMMON  / MMCOMC /  CSTORE
      COMMON  / MMCOMH /  IHASH
      COMMON  / MMCOMS /  ISTORE
      COMMON  / MMCOMI /  ICSDIM,  ICSP1,   ICSP2,   IHADIM,
     A                    ISFREE,  ISTDIM,  IS2HDC,  IS2HDS
C
C         MACRO PROCESSOR INTERFACE
C
      CHARACTER*1         CDIV,    CEOL,     CEOR,    CONC,    CSUB,
     A                    CTOP
      LOGICAL             LEMPTY,  LSUB
      COMMON  / MPCOMC /  CDIV,    CEOL,    CEOR,    CONC,    CSUB,
     A                    CTOP
      COMMON  / MPCOML /  LEMPTY,  LSUB
C
C         LOCAL VARIABLES AND PARAMETERS
C
      CHARACTER*1         C(7), CBI, COL1
      INTEGER             ICDIM, I, IC, ICB
      SAVE ICDIM, C
      EXTERNAL            IOWRLN, IOLIST
      DATA ICDIM / 7 /
      DATA
     A     C(1),  C(2),  C(3),  C(4),  C(5),  C(6),  C(7)
     B  /  ' ',   ')',   ',',   '/',   '*',   '-',   '+'  /
C
      ICB1  =  ICB0
      COL1  =  CBUFFR(ICB1)
      IF  (ICB1 .LE. ICB3)  GO  TO  10
          CBUFFR(ICB1)  =  CBLANK
          ICB2          =  ICB1
          GO  TO  60
C
   10 CONTINUE
          ICB2    =  MIN0(ICB1+ICPLO-1,ICB3)
          IF  (ICB2 .EQ. ICB3)  GO  TO  60
          IF  (.NOT. LBREAK)    GO  TO  40
C
C             FIND A PLACE TO BREAK THE LINE.
C
          DO  30  I=1,10
              CBI  =  CBUFFR(ICB2)
              DO  20  IC=1,ICDIM
                  IF  (C(IC) .EQ. CBI)  GO  TO 30
   20         CONTINUE
              ICB2  =  ICB2 - 1
   30     CONTINUE
C
C             WRITE THE LINE
C
   40     CONTINUE
          CALL  IOWRLN (CBUFFR, ICB1, ICB2, IUNITO)
          IF  (LISTO)  CALL  IOLIST  (.NOT.LISTI)
          ICB1    =  ICB2 + ICBADD
          IF  (.NOT. LFORT)  GO  TO  10
C
C             PAD THE BEGINNING OF THE THE LINE
C             WITH THE STRING BBBBBZBBBB (B=BLANK)
C
          DO  50  ICB=ICB1,ICB2
              CBUFFR(ICB)  =  CBLANK
   50     CONTINUE
          IF  (COL1 .EQ. CC)  CBUFFR(ICB1)    =  CC
          IF  (COL1 .NE. CC)  CBUFFR(ICB1+5)  =  CZ
      GO  TO  10
C
   60 CONTINUE
      CALL  IOWRLN (CBUFFR, ICB1, ICB2, IUNITO)
      IF  (LISTO)  CALL  IOLIST  (.NOT.LISTI)
C
      RETURN
      END
      SUBROUTINE  IOWRLN  (CLINE, ICL1, ICL2, IUNIT)
C
C----------------------------------------------------------------------
C
C     FAMILY
C     ------
C     INPUT/OUTPUT
C
C     PURPOSE
C     -------
C     TO WRITE A LINE FROM THE INPUT/OUTPUT BUFFER.  THIS
C     MAY BE REPLACED BY A MORE EFFICIENT LOCAL I/O ROUTINE
C
C     PARAMETERS
C     ----------
C     CLINE   -I-  I/O BUFFER
C     ICL1    -I-  INDEX OF THE FIRST CHARACTER TO BE WRITTEN
C     ICL2    -I-  INDEX OF THE LAST CHARACTER TO BE WRITTEN
C     IUNIT   -I-  OUTPUT UNIT NUMBER
C
C----------------------------------------------------------------------
      INTEGER ICL1, ICL2, IUNIT
      CHARACTER*(*)      CLINE(ICL2)
      INTEGER I
C
      WRITE (IUNIT, 1010) (CLINE(I), I=ICL1,ICL2)
C
      RETURN
 1010 FORMAT(132A1)
      END
      SUBROUTINE  MMAPPV  (CNAME, ICN1, ICN2, CVALUE, ICV1, ICV2)
C
C----------------------------------------------------------------------
C
C     FAMILY
C     ------
C     MEMORY MANAGER
C
C     PURPOSE
C     -------
C     TO APPEND A STRING TO A VARIABLE
C
C     PARAMETERS
C     ----------
C     CNAME    -I-  ARRAY CONTAINING THE NAME OF THE VARIABLE
C     ICN1     -I-  INDEX OF THE FIRST CHARACTER IN THE NAME
C     ICN2     -I-  INDEX OF THE LAST CHARACTER IN THE NAME
C     CVALUE   -I-  ARRAY CONTAINING THE STRING TO BE APPENDED
C     ICV1     -I-  INDEX OF THE FIRST CHARACTER IN THE STRING
C     ICV2     -I-  INDEX OF THE LAST CHARACTER IN THE STRING
C
C----------------------------------------------------------------------
C
C         GLOBAL CONSTANTS
C
      CHARACTER*1         CA,      CBLANK,  CC,      CI,      CLEFT,
     A                    CMINUS,  CPLUS,   CPOINT,  CQUOTE,  CRIGHT,
     B                    CZ,      C0,      C9
      COMMON  / GLCOMC /  CA,      CBLANK,  CC,      CI,      CLEFT,
     A                    CMINUS,  CPLUS,   CPOINT,  CQUOTE,  CRIGHT,
     B                    CZ,      C0,      C9
C
C         INPUT / OUTPUT CONTROL INTERFACE
C
      CHARACTER*1         CBUFFR(2000)
      LOGICAL             LBREAK,  LFORT,   LISTI,   LISTO
      INTEGER             ICBADD,  ICBEND,  ICBEOL,  ICBSUB,  ICB0,
     A                    ICB1,    ICB2,    ICB3,    ICBDIM,  ICPLI,
     B                    ICPLO,   ILCTR,   ILNMBR,  ILPP,    IPAGE,
     A                    IUNITE,  IUNITI,  IUNITL,  IUNITO
      COMMON  / IOCOMC /  CBUFFR
      COMMON  / IOCOMI /  ICBADD,  ICBEND,  ICBEOL,  ICBSUB,  ICB0,
     A                    ICB1,    ICB2,    ICB3,    ICBDIM,  ICPLI,
     B                    ICPLO,   ILCTR,   ILNMBR,  ILPP,    IPAGE,
     A                    IUNITE,  IUNITI,  IUNITL,  IUNITO
      COMMON  / IOCOML /  LBREAK,  LFORT,   LISTI,   LISTO
C
C         MEMORY MANAGER INTERFACE
C
      CHARACTER*1         CSTORE(20000)
      INTEGER             IHASH(601),       ISTORE(6000)
      INTEGER             ICSDIM,  ICSP1,   ICSP2,   IHADIM,
     A                    ISFREE,  ISTDIM,  IS2HDC,  IS2HDS
      COMMON  / MMCOMC /  CSTORE
      COMMON  / MMCOMH /  IHASH
      COMMON  / MMCOMS /  ISTORE
      COMMON  / MMCOMI /  ICSDIM,  ICSP1,   ICSP2,   IHADIM,
     A                    ISFREE,  ISTDIM,  IS2HDC,  IS2HDS
C
C         LOCAL VARIABLES AND PARAMETERS
C
      INTEGER             ICN1, ICN2, ICV1, ICV2
      CHARACTER*(*)      CNAME(ICN2), CVALUE(ICV2)
      LOGICAL             LFOUND
      INTEGER             IH, IS1, I, IS2
      EXTERNAL            MMHASH, MMNEWI, MMPUT1
C
C         HASH THE VARIABLE NAME TO SEE IF IT EXISTS.
C         IF IT DOES NOT, CREATE IT AND RETURN.
C
      CALL  MMHASH  (CNAME, ICN1, ICN2, IH, LFOUND)
      IF  (LFOUND)  GO  TO  10
          CALL  MMNEWI  (IS1)
          IHASH(IH)  =  IS1
          CALL  MMPUT1  (CNAME,  ICN1, ICN2, ISTORE(IS1), I)
          CALL  MMPUT1  (CVALUE, ICV1, ICV2, ISTORE(IS1+1),
     A                                       ISTORE(IS1+2))
          GO  TO  999
C
C         THE VARIABLE ALREADY EXISTS. APPEND THE VALUE.
C
   10 CONTINUE
      IS1  =  IHASH(IH)
      IS2  =  ISTORE(IS1+2)
      CALL  MMPUT1  (CVALUE, ICV1, ICV2, ISTORE(IS2+2), ISTORE(IS1+2))
C
  999 CONTINUE
      RETURN
      END
 
      SUBROUTINE  MMDELV  (CNAME, ICN1, ICN2, LFOUND)
C
C----------------------------------------------------------------------
C
C     FAMILY
C     ------
C     MEMORY MANAGER
C
C     PURPOSE
C     -------
C     TO DELETE A VARIABLE
C
C     PARAMETERS
C     ----------
C     CNAME   -I-  ARRAY CONTAINING THE NAME OF THE VARIABLE
C     ICN1    -I-  INDEX OF THE FIRST CHARACTER IN THE NAME
C     ICN2    -I-  INDEX OF THE LAST CHARACTER IN THE NAME
C     LFOUND  -O-  TRUE IF THE VARIABLE EXISTED
C
C----------------------------------------------------------------------
C
C         GLOBAL CONSTANTS
C
      CHARACTER*1         CA,      CBLANK,  CC,      CI,      CLEFT,
     A                    CMINUS,  CPLUS,   CPOINT,  CQUOTE,  CRIGHT,
     B                    CZ,      C0,      C9
      COMMON  / GLCOMC /  CA,      CBLANK,  CC,      CI,      CLEFT,
     A                    CMINUS,  CPLUS,   CPOINT,  CQUOTE,  CRIGHT,
     B                    CZ,      C0,      C9
C
C         INPUT / OUTPUT CONTROL INTERFACE
C
      CHARACTER*1         CBUFFR(2000)
      LOGICAL             LBREAK,  LFORT,   LISTI,   LISTO
      INTEGER             ICBADD,  ICBEND,  ICBEOL,  ICBSUB,  ICB0,
     A                    ICB1,    ICB2,    ICB3,    ICBDIM,  ICPLI,
     B                    ICPLO,   ILCTR,   ILNMBR,  ILPP,    IPAGE,
     A                    IUNITE,  IUNITI,  IUNITL,  IUNITO
      COMMON  / IOCOMC /  CBUFFR
      COMMON  / IOCOMI /  ICBADD,  ICBEND,  ICBEOL,  ICBSUB,  ICB0,
     A                    ICB1,    ICB2,    ICB3,    ICBDIM,  ICPLI,
     B                    ICPLO,   ILCTR,   ILNMBR,  ILPP,    IPAGE,
     A                    IUNITE,  IUNITI,  IUNITL,  IUNITO
      COMMON  / IOCOML /  LBREAK,  LFORT,   LISTI,   LISTO
C
C         MEMORY MANAGER INTERFACE
C
      CHARACTER*1         CSTORE(20000)
      INTEGER             IHASH(601),       ISTORE(6000)
      INTEGER             ICSDIM,  ICSP1,   ICSP2,   IHADIM,
     A                    ISFREE,  ISTDIM,  IS2HDC,  IS2HDS
      COMMON  / MMCOMC /  CSTORE
      COMMON  / MMCOMH /  IHASH
      COMMON  / MMCOMS /  ISTORE
      COMMON  / MMCOMI /  ICSDIM,  ICSP1,   ICSP2,   IHADIM,
     A                    ISFREE,  ISTDIM,  IS2HDC,  IS2HDS
C
C         LOCAL VARIABLES AND PARAMETERS
C
      INTEGER             ICN1, ICN2
      CHARACTER*(*)      CNAME(ICN2)
      LOGICAL             LFOUND
      INTEGER             IH, IS1
      EXTERNAL            MMHASH, MMDEL1, MMRETI
C
C         IF THE VARIABLE EXISTS, DELETE IT BY RETURNING THE SPACE
C         TAKEN UP BY IT-S NAME AND VALUE, RETURNING THE SPACE POINTER,
C         AND ZEROING OUT THE HASH TABLE ENTRY.
C
      CALL  MMHASH  (CNAME, ICN1, ICN2, IH, LFOUND)
      IF  (.NOT. LFOUND)  GO  TO  999
      IS1        =  IHASH(IH)
      CALL  MMDEL1  (ISTORE(IS1))
      CALL  MMDEL1  (ISTORE(IS1+1))
      CALL  MMRETI  (IS1)
      IHASH(IH)  =  0
C
  999 CONTINUE
      RETURN
      END
      SUBROUTINE  MMDEL1  (IS2)
C
C----------------------------------------------------------------------
C
C     FAMILY
C     ------
C     MEMORY MANAGER
C
C     PURPOSE
C     -------
C     TO RETURN BLOCKS OF CHARACTER STORAGE TO THE FREE SPACE POOL
C
C     PARAMETERS
C     ----------
C     IS2     -I-  POINTER TO THE FIRST LINK IN A LIST
C                  OF CHARACTER STORAGE BLOCKS
C
C----------------------------------------------------------------------
      INTEGER IS2
C
C         GLOBAL CONSTANTS
C
      CHARACTER*1         CA,      CBLANK,  CC,      CI,      CLEFT,
     A                    CMINUS,  CPLUS,   CPOINT,  CQUOTE,  CRIGHT,
     B                    CZ,      C0,      C9
      COMMON  / GLCOMC /  CA,      CBLANK,  CC,      CI,      CLEFT,
     A                    CMINUS,  CPLUS,   CPOINT,  CQUOTE,  CRIGHT,
     B                    CZ,      C0,      C9
C
C         INPUT / OUTPUT CONTROL INTERFACE
C
      CHARACTER*1         CBUFFR(2000)
      LOGICAL             LBREAK,  LFORT,   LISTI,   LISTO
      INTEGER             ICBADD,  ICBEND,  ICBEOL,  ICBSUB,  ICB0,
     A                    ICB1,    ICB2,    ICB3,    ICBDIM,  ICPLI,
     B                    ICPLO,   ILCTR,   ILNMBR,  ILPP,    IPAGE,
     A                    IUNITE,  IUNITI,  IUNITL,  IUNITO
      COMMON  / IOCOMC /  CBUFFR
      COMMON  / IOCOMI /  ICBADD,  ICBEND,  ICBEOL,  ICBSUB,  ICB0,
     A                    ICB1,    ICB2,    ICB3,    ICBDIM,  ICPLI,
     B                    ICPLO,   ILCTR,   ILNMBR,  ILPP,    IPAGE,
     A                    IUNITE,  IUNITI,  IUNITL,  IUNITO
      COMMON  / IOCOML /  LBREAK,  LFORT,   LISTI,   LISTO
C
C         MEMORY MANAGER INTERFACE
C
      CHARACTER*1         CSTORE(20000)
      INTEGER             IHASH(601),       ISTORE(6000)
      INTEGER             ICSDIM,  ICSP1,   ICSP2,   IHADIM,
     A                    ISFREE,  ISTDIM,  IS2HDC,  IS2HDS
      COMMON  / MMCOMC /  CSTORE
      COMMON  / MMCOMH /  IHASH
      COMMON  / MMCOMS /  ISTORE
      COMMON  / MMCOMI /  ICSDIM,  ICSP1,   ICSP2,   IHADIM,
     A                    ISFREE,  ISTDIM,  IS2HDC,  IS2HDS
      INTEGER IS
C
      IS  =  IS2
      IF  (IS .EQ. 0)  GO  TO  999
C
C         LOOP THROUGH EVERY LINK TO FIND THE TAIL
C
   10 CONTINUE
          IF  (ISTORE(IS+2) .EQ. 0)  GO  TO  20
          IS  =  ISTORE(IS+2)
      GO  TO  10
C
C         ATTACH THE LIST TO THE FREE SPACE POOL AND
C         RESET THE FREE SPACE HEAD POINTER
C
   20 CONTINUE
      ISTORE(IS+2)  =  IS2HDC
      IS2HDC        =  IS2
C
C
  999 CONTINUE
      RETURN
      END
      SUBROUTINE  MMGETV  (CNAME,  ICN1, ICN2,
     A                     CVALUE, ICV1, ICV2, ICVDIM, LFOUND)
C
C----------------------------------------------------------------------
C
C     FAMILY
C     ------
C     MEMORY MANAGER
C
C     PURPOSE
C     -------
C     TO GET THE VALUE OF THE NAMED VARIABLE FROM THE STORAGE
C     POOL AND COPY IT INTO THE SPECIFIED ARRAY.
C
C     PARAMETERS
C     ----------
C     CNAME   -I-  ARRAY CONTAINING THE NAME OF THE VARIABLE
C     ICN1    -I-  INDEX OF THE FIRST CHARACTER IN THE NAME
C     ICN2     -I-  INDEX OF THE LAST CHARACTER IN THE NAME
C     CVALUE   -O-  ARRAY TO CONTAIN THE VALUE OF THE VARIABLE
C     ICV1     -O-  INDEX OF THE FIRST CHARACTER IN THE VALUE
C     ICV2     -O-  INDEX OF THE LAST CHARACTER IN THE VALUE
C     ICVDIM   -I-  LENGTH OF ARRAY CVALUE
C     LFOUND   -O-  TRUE IF THE VARIABLE EXISTS
C
C----------------------------------------------------------------------
C
C         GLOBAL CONSTANTS
C
      CHARACTER*1         CA,      CBLANK,  CC,      CI,      CLEFT,
     A                    CMINUS,  CPLUS,   CPOINT,  CQUOTE,  CRIGHT,
     B                    CZ,      C0,      C9
      COMMON  / GLCOMC /  CA,      CBLANK,  CC,      CI,      CLEFT,
     A                    CMINUS,  CPLUS,   CPOINT,  CQUOTE,  CRIGHT,
     B                    CZ,      C0,      C9
C
C         INPUT / OUTPUT CONTROL INTERFACE
C
      CHARACTER*1         CBUFFR(2000)
      LOGICAL             LBREAK,  LFORT,   LISTI,   LISTO
      INTEGER             ICBADD,  ICBEND,  ICBEOL,  ICBSUB,  ICB0,
     A                    ICB1,    ICB2,    ICB3,    ICBDIM,  ICPLI,
     B                    ICPLO,   ILCTR,   ILNMBR,  ILPP,    IPAGE,
     A                    IUNITE,  IUNITI,  IUNITL,  IUNITO
      COMMON  / IOCOMC /  CBUFFR
      COMMON  / IOCOMI /  ICBADD,  ICBEND,  ICBEOL,  ICBSUB,  ICB0,
     A                    ICB1,    ICB2,    ICB3,    ICBDIM,  ICPLI,
     B                    ICPLO,   ILCTR,   ILNMBR,  ILPP,    IPAGE,
     A                    IUNITE,  IUNITI,  IUNITL,  IUNITO
      COMMON  / IOCOML /  LBREAK,  LFORT,   LISTI,   LISTO
C
C         MEMORY MANAGER INTERFACE
C
      CHARACTER*1         CSTORE(20000)
      INTEGER             IHASH(601),       ISTORE(6000)
      INTEGER             ICSDIM,  ICSP1,   ICSP2,   IHADIM,
     A                    ISFREE,  ISTDIM,  IS2HDC,  IS2HDS
      COMMON  / MMCOMC /  CSTORE
      COMMON  / MMCOMH /  IHASH
      COMMON  / MMCOMS /  ISTORE
      COMMON  / MMCOMI /  ICSDIM,  ICSP1,   ICSP2,   IHADIM,
     A                    ISFREE,  ISTDIM,  IS2HDC,  IS2HDS
C
C         LOCAL VARIABLES AND PARAMETERS
C
      INTEGER             ICN1, ICN2, ICV1, ICV2, ICVDIM
      CHARACTER*(*)      CNAME(ICN2), CVALUE(ICVDIM)
      LOGICAL             LFOUND
      INTEGER             IH, IS1, IS2H
      EXTERNAL            MMHASH, MMGET1
C
C         IF THE VARIABLE EXISTS, COPY ITS VALUE
C
      ICV2  =  0
      CALL  MMHASH  (CNAME, ICN1, ICN2, IH, LFOUND)
      IF  (.NOT. LFOUND)  GO  TO  999
      IS1   =  IHASH(IH)
      IS2H  =  ISTORE(IS1+1)
      CALL  MMGET1  (CVALUE, ICV1, ICV2, ICVDIM, IS2H)
C
  999 CONTINUE
      RETURN
      END
      SUBROUTINE  MMGET1  (CVALUE, ICV1, ICV2, ICVDIM, IS2H)
C
C----------------------------------------------------------------------
C
C     FAMILY
C     ------
C     MEMORY MANAGER
C
C     PURPOSE
C     -------
C     TO COPY THE STRING SPECIFIED BY THE POINTER IS2H
C     AND COPY IT INTO A SPECIFIED ARRAY.
C
C     PARAMETERS
C     ----------
C     CVALUE  -O-  ARRAY TO CONTAIN THE VALUE OF THE VARIABLE
C     ICV1    -O-  INDEX OF THE FIRST CHARACTER IN THE VALUE
C     ICV2    -O-  INDEX OF THE LAST CHARACTER IN THE VALUE
C     ICVDIM  -I-  LENGTH OF ARRAY CVALUE
C     IS2H    -I-  HEAD POINTER TO THE LINKED LIST OF
C                  BLOCKS CONTAINING THE STRING VALUE
C
C----------------------------------------------------------------------
C
C         GLOBAL CONSTANTS
C
      CHARACTER*1         CA,      CBLANK,  CC,      CI,      CLEFT,
     A                    CMINUS,  CPLUS,   CPOINT,  CQUOTE,  CRIGHT,
     B                    CZ,      C0,      C9
      COMMON  / GLCOMC /  CA,      CBLANK,  CC,      CI,      CLEFT,
     A                    CMINUS,  CPLUS,   CPOINT,  CQUOTE,  CRIGHT,
     B                    CZ,      C0,      C9
C
C         INPUT / OUTPUT CONTROL INTERFACE
C
      CHARACTER*1         CBUFFR(2000)
      LOGICAL             LBREAK,  LFORT,   LISTI,   LISTO
      INTEGER             ICBADD,  ICBEND,  ICBEOL,  ICBSUB,  ICB0,
     A                    ICB1,    ICB2,    ICB3,    ICBDIM,  ICPLI,
     B                    ICPLO,   ILCTR,   ILNMBR,  ILPP,    IPAGE,
     A                    IUNITE,  IUNITI,  IUNITL,  IUNITO
      COMMON  / IOCOMC /  CBUFFR
      COMMON  / IOCOMI /  ICBADD,  ICBEND,  ICBEOL,  ICBSUB,  ICB0,
     A                    ICB1,    ICB2,    ICB3,    ICBDIM,  ICPLI,
     B                    ICPLO,   ILCTR,   ILNMBR,  ILPP,    IPAGE,
     A                    IUNITE,  IUNITI,  IUNITL,  IUNITO
      COMMON  / IOCOML /  LBREAK,  LFORT,   LISTI,   LISTO
C
C         MEMORY MANAGER INTERFACE
C
      CHARACTER*1         CSTORE(20000)
      INTEGER             IHASH(601),       ISTORE(6000)
      INTEGER             ICSDIM,  ICSP1,   ICSP2,   IHADIM,
     A                    ISFREE,  ISTDIM,  IS2HDC,  IS2HDS
      COMMON  / MMCOMC /  CSTORE
      COMMON  / MMCOMH /  IHASH
      COMMON  / MMCOMS /  ISTORE
      COMMON  / MMCOMI /  ICSDIM,  ICSP1,   ICSP2,   IHADIM,
     A                    ISFREE,  ISTDIM,  IS2HDC,  IS2HDS
C
C         LOCAL VARIABLES AND PARAMETERS
C
      INTEGER             ICV1, ICV2, ICVDIM, IS2H
      CHARACTER*(*)      CVALUE(ICVDIM)
      INTEGER             ICS1, ICS2, ICS, IS2
      EXTERNAL            IOERRM
C
      IS2   =  IS2H
      ICV2  =  ICV1 - 1
C
C         LOOP THROUGH EACH BLOCK IN WHICH THE STRING IS STORED
C
   10 CONTINUE
          IF  (IS2 .EQ. 0)  GO  TO  999
          ICS1  =  ISTORE(IS2)
          ICS2  =  ISTORE(IS2+1)
          IS2   =  ISTORE(IS2+2)
          IF  (ICV2+ICS2-ICS1 .GE. ICVDIM)  GO  TO  30
C
C             LOOP OVER EACH CHARACTER IN THIS BLOCK
C
          DO  20  ICS=ICS1,ICS2
              ICV2  =  ICV2 + 1
              CVALUE(ICV2)  =  CSTORE(ICS)
   20     CONTINUE
      GO  TO  10
C
   30 CONTINUE
      CALL  IOERRM  (.TRUE.,
     A '('' ********   MMGET1 - STRING TOO LONG FOR CVALUE(*)'')')
C
  999 CONTINUE
      RETURN
      END
      SUBROUTINE  MMHASH  (CNAME, ICN1, ICN2, IH, LFOUND)
C
C----------------------------------------------------------------------
C
C     FAMILY
C     ------
C     MEMORY MANAGER
C
C     PURPOSE
C     -------
C     TO HASH A NAME AND RETURN IT-S HASH TABLE INDEX
C
C     PARAMETERS
C     ----------
C     CNAME   -I-  ARRAY CONTAINING THE NAME OF THE VARIABLE
C     ICN1    -I-  INDEX OF THE FIRST CHARACTER IN THE NAME
C     ICN2    -I-  INDEX OF THE LAST CHARACTER IN THE NAME
C     IH      -O-  HASH INDEX INTO ARRAY IHASH
C     LFOUND  -O-  TRUE IF THE VARIABLE IS ALREADY IN THE TABLE
C
C----------------------------------------------------------------------
C
C         GLOBAL CONSTANTS
C
      CHARACTER*1         CA,      CBLANK,  CC,      CI,      CLEFT,
     A                    CMINUS,  CPLUS,   CPOINT,  CQUOTE,  CRIGHT,
     B                    CZ,      C0,      C9
      COMMON  / GLCOMC /  CA,      CBLANK,  CC,      CI,      CLEFT,
     A                    CMINUS,  CPLUS,   CPOINT,  CQUOTE,  CRIGHT,
     B                    CZ,      C0,      C9
C
C         INPUT / OUTPUT CONTROL INTERFACE
C
      CHARACTER*1         CBUFFR(2000)
      LOGICAL             LBREAK,  LFORT,   LISTI,   LISTO
      INTEGER             ICBADD,  ICBEND,  ICBEOL,  ICBSUB,  ICB0,
     A                    ICB1,    ICB2,    ICB3,    ICBDIM,  ICPLI,
     B                    ICPLO,   ILCTR,   ILNMBR,  ILPP,    IPAGE,
     A                    IUNITE,  IUNITI,  IUNITL,  IUNITO
      COMMON  / IOCOMC /  CBUFFR
      COMMON  / IOCOMI /  ICBADD,  ICBEND,  ICBEOL,  ICBSUB,  ICB0,
     A                    ICB1,    ICB2,    ICB3,    ICBDIM,  ICPLI,
     B                    ICPLO,   ILCTR,   ILNMBR,  ILPP,    IPAGE,
     A                    IUNITE,  IUNITI,  IUNITL,  IUNITO
      COMMON  / IOCOML /  LBREAK,  LFORT,   LISTI,   LISTO
C
C         MEMORY MANAGER INTERFACE
C
      CHARACTER*1         CSTORE(20000)
      INTEGER             IHASH(601),       ISTORE(6000)
      INTEGER             ICSDIM,  ICSP1,   ICSP2,   IHADIM,
     A                    ISFREE,  ISTDIM,  IS2HDC,  IS2HDS
      COMMON  / MMCOMC /  CSTORE
      COMMON  / MMCOMH /  IHASH
      COMMON  / MMCOMS /  ISTORE
      COMMON  / MMCOMI /  ICSDIM,  ICSP1,   ICSP2,   IHADIM,
     A                    ISFREE,  ISTDIM,  IS2HDC,  IS2HDS
C
C         LOCAL VARIABLES AND PARAMETERS
C
      INTEGER             ICN1, ICN2, IH
      CHARACTER*(*)      CNAME(ICN2)
      LOGICAL             LERROR, LFOUND
      INTEGER             INAME, IADD, I, IS1
      EXTERNAL            UTCVNI, MMTEST, IOERRM
C
C         ENCODE THE NAME INTO AN INTEGER
C
      CALL  UTCVNI  (CNAME, ICN1, ICN2, INAME, LERROR)
      INAME   =  MOD(INAME, IHADIM)
      IADD    =  MAX0(1, INAME)
      LFOUND  =  .FALSE.
C
C         LOOP THROUGH ENTRIES IN THE TABLE UNTIL THE
C         NAME IS FOUND OR AN EMPTY BUCKET IS REACHED
C
      DO  10  I=1,IHADIM
          IH     =  INAME + 1
          IS1    =  IHASH(IH)
          IF  (IS1 .EQ. 0)  GO  TO  999
          CALL  MMTEST  (CNAME, ICN1, ICN2, ISTORE(IS1), LFOUND)
          IF  (LFOUND)  GO  TO  999
          INAME  =  MOD(INAME+IADD, IHADIM)
   10 CONTINUE
C
C         EXIT FROM THE ABOVE LOOP INDICATES THAT THE HASH
C         TABLE IS FULL. TO OBTAIN MORE SPACE THE PROCESSOR
C         MUST BE RECOMPILED WITH A LARGER DIMENSION -IHADIM-
C         FOR ARRAY IHASH. IHADIM SHOULD BE A PRIME NUMBER.
C
      CALL  IOERRM  (.TRUE.,
     A '('' ********   MMHASH - HASH TABLE ARRAY IHASH(*) IS FULL'')')
C
  999 CONTINUE
      RETURN
      END
      SUBROUTINE  MMINIT
C
C-----------------------------------------------------------------------
C
C     FAMILY
C     ------
C     MEMORY MANAGER
C
C     PURPOSE
C     -------
C     TO INITIALIZE MEMORY MANAGER VARIABLES
C
C-----------------------------------------------------------------------
C
C         GLOBAL CONSTANTS
C
      CHARACTER*1         CA,      CBLANK,  CC,      CI,      CLEFT,
     A                    CMINUS,  CPLUS,   CPOINT,  CQUOTE,  CRIGHT,
     B                    CZ,      C0,      C9
      COMMON  / GLCOMC /  CA,      CBLANK,  CC,      CI,      CLEFT,
     A                    CMINUS,  CPLUS,   CPOINT,  CQUOTE,  CRIGHT,
     B                    CZ,      C0,      C9
C
C         INPUT / OUTPUT CONTROL INTERFACE
C
      CHARACTER*1         CBUFFR(2000)
      LOGICAL             LBREAK,  LFORT,   LISTI,   LISTO
      INTEGER             ICBADD,  ICBEND,  ICBEOL,  ICBSUB,  ICB0,
     A                    ICB1,    ICB2,    ICB3,    ICBDIM,  ICPLI,
     B                    ICPLO,   ILCTR,   ILNMBR,  ILPP,    IPAGE,
     A                    IUNITE,  IUNITI,  IUNITL,  IUNITO
      COMMON  / IOCOMC /  CBUFFR
      COMMON  / IOCOMI /  ICBADD,  ICBEND,  ICBEOL,  ICBSUB,  ICB0,
     A                    ICB1,    ICB2,    ICB3,    ICBDIM,  ICPLI,
     B                    ICPLO,   ILCTR,   ILNMBR,  ILPP,    IPAGE,
     A                    IUNITE,  IUNITI,  IUNITL,  IUNITO
      COMMON  / IOCOML /  LBREAK,  LFORT,   LISTI,   LISTO
C
C         MEMORY MANAGER INTERFACE
C
      CHARACTER*1         CSTORE(20000)
      INTEGER             IHASH(601),       ISTORE(6000)
      INTEGER             ICSDIM,  ICSP1,   ICSP2,   IHADIM,
     A                    ISFREE,  ISTDIM,  IS2HDC,  IS2HDS
      COMMON  / MMCOMC /  CSTORE
      COMMON  / MMCOMH /  IHASH
      COMMON  / MMCOMS /  ISTORE
      COMMON  / MMCOMI /  ICSDIM,  ICSP1,   ICSP2,   IHADIM,
     A                    ISFREE,  ISTDIM,  IS2HDC,  IS2HDS
      INTEGER  I
      EXTERNAL MMNEWI
C
      DO  10  I=1,IHADIM
          IHASH(I)      =  0
   10 CONTINUE
C
      DO  20  I=1,ISTDIM,3
          ISTORE(I)     =  0
          ISTORE(I+1)   =  0
          ISTORE(I+2)   =  I + 3
   20 CONTINUE
C
      ISTORE(ISTDIM)    =  0
      ISFREE            =  1
C
      CALL  MMNEWI  (IS2HDC)
      ISTORE(IS2HDC)    =  1
      ISTORE(IS2HDC+1)  =  ICSDIM
      ISTORE(IS2HDC+2)  =  0
      IS2HDS            =  0
      ICSP1             =  1
      ICSP2             =  0
C
      RETURN
      END
      SUBROUTINE  MMNEWI  (IS)
C
C-----------------------------------------------------------------------
C
C     FAMILY
C     ------
C     MEMORY MANAGER
C
C     PURPOSE
C     -------
C     TO RETURN A POINTER TO AN AVAILABLE BLOCK FROM THE INTEGER
C     STORAGE POOL
C
C     PARAMETERS
C     ----------
C     IS      -O-  INDEX INTO ARRAY ISTORE OF THE FREE BLOCK
C
C-----------------------------------------------------------------------
C
C         GLOBAL CONSTANTS
C
      CHARACTER*1         CA,      CBLANK,  CC,      CI,      CLEFT,
     A                    CMINUS,  CPLUS,   CPOINT,  CQUOTE,  CRIGHT,
     B                    CZ,      C0,      C9
      COMMON  / GLCOMC /  CA,      CBLANK,  CC,      CI,      CLEFT,
     A                    CMINUS,  CPLUS,   CPOINT,  CQUOTE,  CRIGHT,
     B                    CZ,      C0,      C9
C
C         INPUT / OUTPUT CONTROL INTERFACE
C
      CHARACTER*1         CBUFFR(2000)
      LOGICAL             LBREAK,  LFORT,   LISTI,   LISTO
      INTEGER             ICBADD,  ICBEND,  ICBEOL,  ICBSUB,  ICB0,
     A                    ICB1,    ICB2,    ICB3,    ICBDIM,  ICPLI,
     B                    ICPLO,   ILCTR,   ILNMBR,  ILPP,    IPAGE,
     A                    IUNITE,  IUNITI,  IUNITL,  IUNITO
      COMMON  / IOCOMC /  CBUFFR
      COMMON  / IOCOMI /  ICBADD,  ICBEND,  ICBEOL,  ICBSUB,  ICB0,
     A                    ICB1,    ICB2,    ICB3,    ICBDIM,  ICPLI,
     B                    ICPLO,   ILCTR,   ILNMBR,  ILPP,    IPAGE,
     A                    IUNITE,  IUNITI,  IUNITL,  IUNITO
      COMMON  / IOCOML /  LBREAK,  LFORT,   LISTI,   LISTO
C
C         MEMORY MANAGER INTERFACE
C
      CHARACTER*1         CSTORE(20000)
      INTEGER             IHASH(601),       ISTORE(6000)
      INTEGER             ICSDIM,  ICSP1,   ICSP2,   IHADIM,
     A                    ISFREE,  ISTDIM,  IS2HDC,  IS2HDS
      COMMON  / MMCOMC /  CSTORE
      COMMON  / MMCOMH /  IHASH
      COMMON  / MMCOMS /  ISTORE
      COMMON  / MMCOMI /  ICSDIM,  ICSP1,   ICSP2,   IHADIM,
     A                    ISFREE,  ISTDIM,  IS2HDC,  IS2HDS
        INTEGER IS
        EXTERNAL IOERRM
C
      IF  (ISFREE .EQ. 0)  GO  TO  10
      IS      =  ISFREE
      ISFREE  =  ISTORE(ISFREE+2)
      GO  TO  999
C
   10 CONTINUE
      CALL  IOERRM  (.TRUE.,
     A '('' ********   MMNEWI - STORAGE ARRAY ISTORE(*) IS FULL'')')
C
  999 CONTINUE
      RETURN
      END
      SUBROUTINE  MMPOPC  (CTEST, IPOP, CTOP, LEMPTY)
C
C-----------------------------------------------------------------------
C
C     FAMILY
C     ------
C     MEMORY MANAGER
C
C     PURPOSE
C     -------
C     TO POP CHARACTERS OFF THE SUBSTITUTION STACK
C
C     PARAMETERS
C     ----------
C     CTEST   -I-  CHARACTER WHOSE PURPOSE DEPENDS ON IPOP
C     IPOP    -I-  INDICATES THE OPERATION TO BE PERFORMED
C                  1 - LOOK AT THE TOP CHARACTER
C                  2 - POP ONE CHARACTER OFF THE STACK
C                  3 - POP ONE VARIABLE OFF THE STACK
C                  4 - POP UNTIL TOP .NE. CTEST
C                  5 - POP UNTIL TOP .EQ. CTEST
C                  6 - POP ALL ALPHNUMERICS
C     CTOP    -O-  TOP CHARACTER ON STACK
C     LEMPTY  -I-  TRUE IF STACK IS EMPTY
C
C-----------------------------------------------------------------------
C
C         GLOBAL CONSTANTS
C
      CHARACTER*1         CA,      CBLANK,  CC,      CI,      CLEFT,
     A                    CMINUS,  CPLUS,   CPOINT,  CQUOTE,  CRIGHT,
     B                    CZ,      C0,      C9
      COMMON  / GLCOMC /  CA,      CBLANK,  CC,      CI,      CLEFT,
     A                    CMINUS,  CPLUS,   CPOINT,  CQUOTE,  CRIGHT,
     B                    CZ,      C0,      C9
C
C         INPUT / OUTPUT CONTROL INTERFACE
C
      CHARACTER*1         CBUFFR(2000)
      LOGICAL             LBREAK,  LFORT,   LISTI,   LISTO
      INTEGER             ICBADD,  ICBEND,  ICBEOL,  ICBSUB,  ICB0,
     A                    ICB1,    ICB2,    ICB3,    ICBDIM,  ICPLI,
     B                    ICPLO,   ILCTR,   ILNMBR,  ILPP,    IPAGE,
     A                    IUNITE,  IUNITI,  IUNITL,  IUNITO
      COMMON  / IOCOMC /  CBUFFR
      COMMON  / IOCOMI /  ICBADD,  ICBEND,  ICBEOL,  ICBSUB,  ICB0,
     A                    ICB1,    ICB2,    ICB3,    ICBDIM,  ICPLI,
     B                    ICPLO,   ILCTR,   ILNMBR,  ILPP,    IPAGE,
     A                    IUNITE,  IUNITI,  IUNITL,  IUNITO
      COMMON  / IOCOML /  LBREAK,  LFORT,   LISTI,   LISTO
C
C         MEMORY MANAGER INTERFACE
C
      CHARACTER*1         CSTORE(20000)
      INTEGER             IHASH(601),       ISTORE(6000)
      INTEGER             ICSDIM,  ICSP1,   ICSP2,   IHADIM,
     A                    ISFREE,  ISTDIM,  IS2HDC,  IS2HDS
      COMMON  / MMCOMC /  CSTORE
      COMMON  / MMCOMH /  IHASH
      COMMON  / MMCOMS /  ISTORE
      COMMON  / MMCOMI /  ICSDIM,  ICSP1,   ICSP2,   IHADIM,
     A                    ISFREE,  ISTDIM,  IS2HDC,  IS2HDS
C
C         LOCAL VARIABLES AND PARAMETERS
C
      INTEGER                IPOP
      CHARACTER*(*)         CTEST, CTOP
      INTEGER                ICS
      LOGICAL             LEMPTY
      EXTERNAL            MMPOP1, MMPOPV, IOERRM
C
   10 CONTINUE
      CTOP  =  CBLANK
C
C         CHECK FOR NULL ENTRIES ON STACK
C
      IF  (ICSP1 .GT. ICSP2)  CALL  MMPOP1  (LEMPTY)
      IF  (LEMPTY)  GO  TO  999
      GO  TO  (20, 30, 40, 50, 70, 90), IPOP
C
C         IPOP = 1  -  LOOK AT THE TOP OF THE STACK
C
   20 CONTINUE
      CTOP  =  CSTORE(ICSP1)
      GO  TO  999
C
C         IPOP = 2  -  POP ONE CHARACTER OFF THE STACK
C
   30 CONTINUE
      ICB2              =  ICB2 + 1
      IF  (ICB2 .GT. ICBDIM)  GO  TO  130
      CBUFFR(ICB2)      =  CSTORE(ICSP1)
      ICSP1             =  ICSP1 + 1
      ISTORE(IS2HDS+1)  =  ICSP1
      IF  (ICSP1 .GT. ICSP2)  CALL  MMPOP1  (LEMPTY)
      IF  (.NOT. LEMPTY)   CTOP  =  CSTORE(ICSP1)
      GO  TO  999
C
C         IPOP = 3  -  POP ONE VARIABLE OFF THE STACK
C
   40 CONTINUE
      ICB2              =  ICB2 + 1
      IF  (ICB2 .GT. ICBDIM)  GO  TO  130
      CBUFFR(ICB2)      =  CSTORE(ICSP1)
      ISTORE(IS2HDS+1)  =  ICSP1 + 1
      CALL  MMPOPV  (LEMPTY)
      CALL  MMPOP1  (LEMPTY)
      IF  (.NOT. LEMPTY)    CTOP  =  CSTORE(ICSP1)
      GO  TO  999
C
C         IPOP = 4  -  POP UNTIL TOP CHAR .NE. CTEST
C
   50 CONTINUE
      IF  (ICSP2-ICSP1 .GE. ICBDIM-ICB2)  GO  TO  130
      DO  60  ICS=ICSP1,ICSP2
          IF  (CSTORE(ICS) .NE. CTEST)  GO  TO  120
          ICB2          =  ICB2 + 1
          CBUFFR(ICB2)  =  CSTORE(ICS)
   60 CONTINUE
      GO  TO  110
C
C         IPOP = 5  -  POP UNTIL TOP CHAR .EQ. CTEST
C
   70 CONTINUE
      IF  (ICSP2-ICSP1 .GE. ICBDIM-ICB2)  GO  TO  130
      DO  80  ICS=ICSP1,ICSP2
          IF  (CSTORE(ICS) .EQ. CTEST)   GO  TO  120
          ICB2          =  ICB2 + 1
          CBUFFR(ICB2)  =  CSTORE(ICS)
   80 CONTINUE
      GO  TO  110
C
C         IPOP = 6  -  POP ALL ALPHANUMERICS OFF THE STACK
C
   90 CONTINUE
      IF  (ICSP2-ICSP1 .GE. ICBDIM-ICB2)  GO  TO  130
      DO  100  ICS=ICSP1,ICSP2
          IF  (.NOT. ((LLE(CA,CSTORE(ICS))
     A         .AND.   LLE(CSTORE(ICS),CZ))
     B         .OR.   (LLE(C0,CSTORE(ICS))
     C         .AND.   LLE(CSTORE(ICS),C9)))) GO TO 120
          ICB2          =  ICB2 + 1
          CBUFFR(ICB2)  =  CSTORE(ICS)
  100 CONTINUE
C
C         THE SPECIFIED CONDITION HAS NOT BEEN MET.
C         GET ANOTHER PIECE OF THE STACK AND TRY AGAIN.
C
  110 CONTINUE
      ICSP1             =  ICSP2 + 1
      ISTORE(IS2HDS+1)  =  ICSP1
      GO  TO  10
C
C         THE SPECIFIED CONDITION HAS BEEN MET.
C         SAVE THE STACK POINTER AND RETURN.
C
  120 CONTINUE
      ICSP1             =  ICS
      ISTORE(IS2HDS+1)  =  ICS
      CTOP              =  CSTORE(ICS)
      GO  TO  999
C
C         THE BUFFER SPACE HAS BEEN EXCEEDED
C
  130 CONTINUE
      CALL  IOERRM  (.TRUE.,
     A '('' ********   MMPOPC - STRING TOO LONG FOR BUFFER'')')
C
  999 CONTINUE
      RETURN
      END
      SUBROUTINE  MMPOPV  (LEMPTY)
C
C-----------------------------------------------------------------------
C
C     FAMILY
C     ------
C     MEMORY MANAGER
C
C     PURPOSE
C     -------
C     TO POP A VARIABLE OFF THE SUBSTITUTION STACK
C
C     PARAMETERS
C     ----------
C     LEMPTY  -O-  TRUE IF THE STACK IS EMPTY
C
C-----------------------------------------------------------------------
C
C         GLOBAL CONSTANTS
C
      CHARACTER*1         CA,      CBLANK,  CC,      CI,      CLEFT,
     A                    CMINUS,  CPLUS,   CPOINT,  CQUOTE,  CRIGHT,
     B                    CZ,      C0,      C9
      COMMON  / GLCOMC /  CA,      CBLANK,  CC,      CI,      CLEFT,
     A                    CMINUS,  CPLUS,   CPOINT,  CQUOTE,  CRIGHT,
     B                    CZ,      C0,      C9
C
C         INPUT / OUTPUT CONTROL INTERFACE
C
      CHARACTER*1         CBUFFR(2000)
      LOGICAL             LBREAK,  LFORT,   LISTI,   LISTO
      INTEGER             ICBADD,  ICBEND,  ICBEOL,  ICBSUB,  ICB0,
     A                    ICB1,    ICB2,    ICB3,    ICBDIM,  ICPLI,
     B                    ICPLO,   ILCTR,   ILNMBR,  ILPP,    IPAGE,
     A                    IUNITE,  IUNITI,  IUNITL,  IUNITO
      COMMON  / IOCOMC /  CBUFFR
      COMMON  / IOCOMI /  ICBADD,  ICBEND,  ICBEOL,  ICBSUB,  ICB0,
     A                    ICB1,    ICB2,    ICB3,    ICBDIM,  ICPLI,
     B                    ICPLO,   ILCTR,   ILNMBR,  ILPP,    IPAGE,
     A                    IUNITE,  IUNITI,  IUNITL,  IUNITO
      COMMON  / IOCOML /  LBREAK,  LFORT,   LISTI,   LISTO
C
C         MEMORY MANAGER INTERFACE
C
      CHARACTER*1         CSTORE(20000)
      INTEGER             IHASH(601),       ISTORE(6000)
      INTEGER             ICSDIM,  ICSP1,   ICSP2,   IHADIM,
     A                    ISFREE,  ISTDIM,  IS2HDC,  IS2HDS
      COMMON  / MMCOMC /  CSTORE
      COMMON  / MMCOMH /  IHASH
      COMMON  / MMCOMS /  ISTORE
      COMMON  / MMCOMI /  ICSDIM,  ICSP1,   ICSP2,   IHADIM,
     A                    ISFREE,  ISTDIM,  IS2HDC,  IS2HDS
C
C         LOCAL VARIABLES AND PARAMETERS
C
      LOGICAL             LEMPTY
      INTEGER             IS2
      EXTERNAL            MMRETI
C
      LEMPTY         =  IS2HDS .EQ. 0
      IF  (LEMPTY)  GO  TO  999
      IS2            =  IS2HDS
      IS2HDS         =  ISTORE(IS2+2)
      ISTORE(IS2+2)  =  0
      LEMPTY         =  IS2HDS .EQ. 0
      IF  (ISTORE(IS2) .GT. 0)  CALL  MMRETI  (IS2)
C
  999 CONTINUE
      RETURN
      END
      SUBROUTINE  MMPOP1  (LEMPTY)
C
C-----------------------------------------------------------------------
C
C     FAMILY
C     ------
C     MEMORY MANAGER
C
C     PURPOSE
C     -------
C     TO POP NULL ENTRIES OFF THE SUBSTITUTION STACK
C
C     PARAMETERS
C     ----------
C     LEMPTY  -O-  TRUE IF THE STACK IS EMPTY
C
C-----------------------------------------------------------------------
C
C         GLOBAL CONSTANTS
C
      CHARACTER*1         CA,      CBLANK,  CC,      CI,      CLEFT,
     A                    CMINUS,  CPLUS,   CPOINT,  CQUOTE,  CRIGHT,
     B                    CZ,      C0,      C9
      COMMON  / GLCOMC /  CA,      CBLANK,  CC,      CI,      CLEFT,
     A                    CMINUS,  CPLUS,   CPOINT,  CQUOTE,  CRIGHT,
     B                    CZ,      C0,      C9
C
C         INPUT / OUTPUT CONTROL INTERFACE
C
      CHARACTER*1         CBUFFR(2000)
      LOGICAL             LBREAK,  LFORT,   LISTI,   LISTO
      INTEGER             ICBADD,  ICBEND,  ICBEOL,  ICBSUB,  ICB0,
     A                    ICB1,    ICB2,    ICB3,    ICBDIM,  ICPLI,
     B                    ICPLO,   ILCTR,   ILNMBR,  ILPP,    IPAGE,
     A                    IUNITE,  IUNITI,  IUNITL,  IUNITO
      COMMON  / IOCOMC /  CBUFFR
      COMMON  / IOCOMI /  ICBADD,  ICBEND,  ICBEOL,  ICBSUB,  ICB0,
     A                    ICB1,    ICB2,    ICB3,    ICBDIM,  ICPLI,
     B                    ICPLO,   ILCTR,   ILNMBR,  ILPP,    IPAGE,
     A                    IUNITE,  IUNITI,  IUNITL,  IUNITO
      COMMON  / IOCOML /  LBREAK,  LFORT,   LISTI,   LISTO
C
C         MEMORY MANAGER INTERFACE
C
      CHARACTER*1         CSTORE(20000)
      INTEGER             IHASH(601),       ISTORE(6000)
      INTEGER             ICSDIM,  ICSP1,   ICSP2,   IHADIM,
     A                    ISFREE,  ISTDIM,  IS2HDC,  IS2HDS
      COMMON  / MMCOMC /  CSTORE
      COMMON  / MMCOMH /  IHASH
      COMMON  / MMCOMS /  ISTORE
      COMMON  / MMCOMI /  ICSDIM,  ICSP1,   ICSP2,   IHADIM,
     A                    ISFREE,  ISTDIM,  IS2HDC,  IS2HDS
C
C         LOCAL VARIABLES AND PARAMETERS
C
      LOGICAL             LEMPTY
      INTEGER             IS2
      EXTERNAL            MMPOPV
C
   10 CONTINUE
          LEMPTY            =  IS2HDS .EQ. 0
          IF  (LEMPTY)  GO  TO  999
          IS2               =  IABS(ISTORE(IS2HDS))
          IF  (IS2 .NE. 0)  GO  TO  30
   20     CONTINUE
          CALL  MMPOPV  (LEMPTY)
      GO  TO  10
C
   30 CONTINUE
          ICSP1              =  ISTORE(IS2HDS+1)
          ICSP2              =  ISTORE(IS2+1)
          IF  (ICSP1 .LE. ICSP2)  GO  TO  999
          IS2               =  ISTORE(IS2+2)
          IF  (IS2 .EQ. 0)  GO  TO  20
          ISTORE(IS2HDS)    =  ISIGN(IS2, ISTORE(IS2HDS))
          ISTORE(IS2HDS+1)  =  ISTORE(IS2)
      GO  TO  30
C
  999 CONTINUE
      RETURN
      END
      SUBROUTINE  MMPSHV  (CNAME, ICN1, ICN2, IPUSH, LEMPTY, LFOUND)
C
C-----------------------------------------------------------------------
C
C     FAMILY
C     ------
C     MEMORY MANAGER
C
C     PURPOSE
C     -------
C     TO PUSH A VARIABLE ONTO THE SUBSTITUTION STACK
C
C     PARAMETERS
C     ----------
C     CNAME   -I-  THE NAME OF THE VARIABLE TO PUSH ONTO THE STACK
C     ICN1    -I-  INDEX OF THE FIRST CHARACTER IN THE NAME
C     ICN2    -I-  INDEX OF THE LAST CHARACTER IN THE NAME
C     IPUSH   -I-  INDICATES THE OPERATION TO BE PERFORMED
C                  1 - PUSH A VARIABLE ONTO THE STACK
C                  2 - PUSH A POINTER ONTO THE STACK
C                  3 - PUSH THE ACTUAL POINTER ONTO THE STACK
C
C-----------------------------------------------------------------------
C
C         GLOBAL CONSTANTS
C
      CHARACTER*1         CA,      CBLANK,  CC,      CI,      CLEFT,
     A                    CMINUS,  CPLUS,   CPOINT,  CQUOTE,  CRIGHT,
     B                    CZ,      C0,      C9
      COMMON  / GLCOMC /  CA,      CBLANK,  CC,      CI,      CLEFT,
     A                    CMINUS,  CPLUS,   CPOINT,  CQUOTE,  CRIGHT,
     B                    CZ,      C0,      C9
C
C         INPUT / OUTPUT CONTROL INTERFACE
C
      CHARACTER*1         CBUFFR(2000)
      LOGICAL             LBREAK,  LFORT,   LISTI,   LISTO
      INTEGER             ICBADD,  ICBEND,  ICBEOL,  ICBSUB,  ICB0,
     A                    ICB1,    ICB2,    ICB3,    ICBDIM,  ICPLI,
     B                    ICPLO,   ILCTR,   ILNMBR,  ILPP,    IPAGE,
     A                    IUNITE,  IUNITI,  IUNITL,  IUNITO
      COMMON  / IOCOMC /  CBUFFR
      COMMON  / IOCOMI /  ICBADD,  ICBEND,  ICBEOL,  ICBSUB,  ICB0,
     A                    ICB1,    ICB2,    ICB3,    ICBDIM,  ICPLI,
     B                    ICPLO,   ILCTR,   ILNMBR,  ILPP,    IPAGE,
     A                    IUNITE,  IUNITI,  IUNITL,  IUNITO
      COMMON  / IOCOML /  LBREAK,  LFORT,   LISTI,   LISTO
C
C         MEMORY MANAGER INTERFACE
C
      CHARACTER*1         CSTORE(20000)
      INTEGER             IHASH(601),       ISTORE(6000)
      INTEGER             ICSDIM,  ICSP1,   ICSP2,   IHADIM,
     A                    ISFREE,  ISTDIM,  IS2HDC,  IS2HDS
      COMMON  / MMCOMC /  CSTORE
      COMMON  / MMCOMH /  IHASH
      COMMON  / MMCOMS /  ISTORE
      COMMON  / MMCOMI /  ICSDIM,  ICSP1,   ICSP2,   IHADIM,
     A                    ISFREE,  ISTDIM,  IS2HDC,  IS2HDS
C
C         LOCAL VARIABLES AND PARAMETERS
C
      INTEGER             ICN1, ICN2, IPUSH
      CHARACTER*(*)      CNAME(ICN2)
      LOGICAL             LEMPTY, LFOUND
      INTEGER             IH, IS1, IS2, ITEMP
      EXTERNAL            MMHASH, MMNEWI, MMPOP1
C
      CALL  MMHASH  (CNAME, ICN1, ICN2, IH, LFOUND)
      IF  (.NOT. LFOUND)  GO  TO  999
      IS1              =  IHASH(IH)
      IS2              =  ISTORE(IS1+1)
      IF  (IS2 .EQ. 0)    GO  TO  999
      IF  (IPUSH .EQ. 1)  GO  TO  10
      IF  (IPUSH .EQ. 2)  GO  TO  20
      GO  TO  30
C
C         PUSH A VARIABLE ONTO THE STACK; NEW ENTRY WILL POINT TO
C         VALUE OF THE VARIABLE
C
   10 CONTINUE
      CALL  MMNEWI  (ITEMP)
      ISTORE(ITEMP)    =  IS2
      ISTORE(ITEMP+1)  =  ISTORE(IS2)
      ISTORE(ITEMP+2)  =  IS2HDS
      IS2HDS           =  ITEMP
      GO  TO  40
C
C         PUSH A POINTER ONTO THE STACK
C
   20 CONTINUE
      CALL  MMNEWI  (ITEMP)
      ISTORE(ITEMP)    =  ISTORE(IS2)
      ISTORE(ITEMP+1)  =  ISTORE(IS2+1)
      ISTORE(ITEMP+2)  =  IS2HDS
      IS2HDS           =  ITEMP
      GO  TO  40
C
C         PUSH THE ACTUAL POINTER ONTO THE STACK
C
   30 CONTINUE
      ISTORE(IS2)      =  -IABS(ISTORE(IS2))
      ISTORE(IS2+2)    =  IS2HDS
      IS2HDS           =  IS2
C
C         CALL MMPOP1 TO SET THE POINTERS (ICSP1, ICSP2) INTO CSTORE
C
   40 CONTINUE
      CALL  MMPOP1  (LEMPTY)
C
  999 CONTINUE
      RETURN
      END
      SUBROUTINE  MMPUTP  (CNAME, ICN1, ICN2, CPTR, ICP1, ICP2, LFOUND)
C
C-----------------------------------------------------------------------
C
C     FAMILY
C     ------
C     MEMORY MANAGER
C
C     PURPOSE
C     -------
C     TO PUT A POINTER TO A VARIABLE IN THE SYMBOL TABLE
C
C     PARAMETERS
C     ----------
C     CNAME   -I-  NAME OF THE VARIABLE
C     ICN1    -I-  INDEX OF THE FIRST CHARACTER IN THE NAME
C     ICN2    -I-  INDEX OF THE LAST CHARACTER IN THE NAME
C     CPTR    -I-  NAME OF THE POINTER
C     ICP1    -I-  INDEX OF THE FIRST CHARACTER IN THE POINTER NAME
C     ICP2    -I-  INDEX OF THE LAST CHARACTER IN THE POINTER NAME
C     LFOUND  -O-  TRUE IF THE VARIABLE WAS FOUND IN THE TABLE
C
C-----------------------------------------------------------------------
C
C         GLOBAL CONSTANTS
C
      CHARACTER*1         CA,      CBLANK,  CC,      CI,      CLEFT,
     A                    CMINUS,  CPLUS,   CPOINT,  CQUOTE,  CRIGHT,
     B                    CZ,      C0,      C9
      COMMON  / GLCOMC /  CA,      CBLANK,  CC,      CI,      CLEFT,
     A                    CMINUS,  CPLUS,   CPOINT,  CQUOTE,  CRIGHT,
     B                    CZ,      C0,      C9
C
C         INPUT / OUTPUT CONTROL INTERFACE
C
      CHARACTER*1         CBUFFR(2000)
      LOGICAL             LBREAK,  LFORT,   LISTI,   LISTO
      INTEGER             ICBADD,  ICBEND,  ICBEOL,  ICBSUB,  ICB0,
     A                    ICB1,    ICB2,    ICB3,    ICBDIM,  ICPLI,
     B                    ICPLO,   ILCTR,   ILNMBR,  ILPP,    IPAGE,
     A                    IUNITE,  IUNITI,  IUNITL,  IUNITO
      COMMON  / IOCOMC /  CBUFFR
      COMMON  / IOCOMI /  ICBADD,  ICBEND,  ICBEOL,  ICBSUB,  ICB0,
     A                    ICB1,    ICB2,    ICB3,    ICBDIM,  ICPLI,
     B                    ICPLO,   ILCTR,   ILNMBR,  ILPP,    IPAGE,
     A                    IUNITE,  IUNITI,  IUNITL,  IUNITO
      COMMON  / IOCOML /  LBREAK,  LFORT,   LISTI,   LISTO
C
C         MEMORY MANAGER INTERFACE
C
      CHARACTER*1         CSTORE(20000)
      INTEGER             IHASH(601),       ISTORE(6000)
      INTEGER             ICSDIM,  ICSP1,   ICSP2,   IHADIM,
     A                    ISFREE,  ISTDIM,  IS2HDC,  IS2HDS
      COMMON  / MMCOMC /  CSTORE
      COMMON  / MMCOMH /  IHASH
      COMMON  / MMCOMS /  ISTORE
      COMMON  / MMCOMI /  ICSDIM,  ICSP1,   ICSP2,   IHADIM,
     A                    ISFREE,  ISTDIM,  IS2HDC,  IS2HDS
C
C         LOCAL VARIABLES AND PARAMETERS
C
      INTEGER             ICN1, ICN2, ICP1, ICP2
      CHARACTER*(*)      CNAME(ICN2), CPTR(ICP2)
      LOGICAL             L, LFOUND
      INTEGER             IH, IS1, IS2CN, I, IS2
      EXTERNAL            MMHASH, MMNEWI, MMPUT1
C
      CALL  MMHASH  (CNAME, ICN1, ICN2, IH, LFOUND)
      IF  (.NOT. LFOUND)  GO  TO  999
      IS1            =  IHASH(IH)
      IS2CN          =  ISTORE(IS1+1)
      CALL  MMHASH  (CPTR, ICP1, ICP2, IH, L)
      IF  (L)  GO  TO  10
          CALL  MMNEWI  (IS1)
          IHASH(IH)      =  IS1
          CALL  MMPUT1  (CPTR, ICP1, ICP2, ISTORE(IS1), I)
          CALL  MMNEWI  (IS2)
          ISTORE(IS1+1)  =  IS2
          ISTORE(IS1+2)  =  0
          GO  TO  20
C
   10 CONTINUE
      IS1            =  IHASH(IH)
      IS2            =  ISTORE(IS1+1)
C
   20 CONTINUE
      IF  (IS2CN .NE. 0)  GO  TO  30
      ISTORE(IS2)    =  0
      ISTORE(IS2+1)  =  0
      ISTORE(IS2+2)  =  0
      GO  TO  999
C
   30 CONTINUE
      ISTORE(IS2)    =  IS2CN
      ISTORE(IS2+1)  =  ISTORE(IS2CN)
      ISTORE(IS2+2)  =  0
C
  999 CONTINUE
      RETURN
      END
      SUBROUTINE  MMPUTV  (CNAME, ICN1, ICN2, CVALUE, ICV1, ICV2)
C
C-----------------------------------------------------------------------
C
C     FAMILY
C     ------
C     MEMORY MANAGER
C
C     PURPOSE
C     -------
C     TO PUT A VARIABLE INTO THE SYMBOL TABLE
C
C     PARAMETERS
C     ----------
C     CNAME   -I-  NAME OF THE VARIABLE
C     ICN1    -I-  INDEX OF THE FIRST CHARACTER IN THE NAME
C     ICN2    -I-  INDEX OF THE LAST CHARACTER IN THE NAME
C     CVALUE  -I-  VALUE OF THE VARIABLE
C     ICV1    -I-  INDEX OF THE FIRST CHARACTER IN THE VALUE
C     ICV2    -I-  INDEX OF THE LAST CHARACTER IN THE VALUE
C
C-----------------------------------------------------------------------
C
C         GLOBAL CONSTANTS
C
      CHARACTER*1         CA,      CBLANK,  CC,      CI,      CLEFT,
     A                    CMINUS,  CPLUS,   CPOINT,  CQUOTE,  CRIGHT,
     B                    CZ,      C0,      C9
      COMMON  / GLCOMC /  CA,      CBLANK,  CC,      CI,      CLEFT,
     A                    CMINUS,  CPLUS,   CPOINT,  CQUOTE,  CRIGHT,
     B                    CZ,      C0,      C9
C
C         INPUT / OUTPUT CONTROL INTERFACE
C
      CHARACTER*1         CBUFFR(2000)
      LOGICAL             LBREAK,  LFORT,   LISTI,   LISTO
      INTEGER             ICBADD,  ICBEND,  ICBEOL,  ICBSUB,  ICB0,
     A                    ICB1,    ICB2,    ICB3,    ICBDIM,  ICPLI,
     B                    ICPLO,   ILCTR,   ILNMBR,  ILPP,    IPAGE,
     A                    IUNITE,  IUNITI,  IUNITL,  IUNITO
      COMMON  / IOCOMC /  CBUFFR
      COMMON  / IOCOMI /  ICBADD,  ICBEND,  ICBEOL,  ICBSUB,  ICB0,
     A                    ICB1,    ICB2,    ICB3,    ICBDIM,  ICPLI,
     B                    ICPLO,   ILCTR,   ILNMBR,  ILPP,    IPAGE,
     A                    IUNITE,  IUNITI,  IUNITL,  IUNITO
      COMMON  / IOCOML /  LBREAK,  LFORT,   LISTI,   LISTO
C
C         MEMORY MANAGER INTERFACE
C
      CHARACTER*1         CSTORE(20000)
      INTEGER             IHASH(601),       ISTORE(6000)
      INTEGER             ICSDIM,  ICSP1,   ICSP2,   IHADIM,
     A                    ISFREE,  ISTDIM,  IS2HDC,  IS2HDS
      COMMON  / MMCOMC /  CSTORE
      COMMON  / MMCOMH /  IHASH
      COMMON  / MMCOMS /  ISTORE
      COMMON  / MMCOMI /  ICSDIM,  ICSP1,   ICSP2,   IHADIM,
     A                    ISFREE,  ISTDIM,  IS2HDC,  IS2HDS
C
C         LOCAL VARIABLES AND PARAMETERS
C
      INTEGER             ICN1, ICN2, ICV1, ICV2
      CHARACTER*(*)      CNAME(ICN2), CVALUE(ICV2)
      LOGICAL             LFOUND
      INTEGER             IH, IS1, I
      EXTERNAL            MMHASH, MMNEWI, MMPUT1, MMDEL1
C
C         HASH THE NAME TO SEE IF IT IS IN THE TABLE.
C         IF IT IS NOT, STORE A NEW NAME IN THE TABLE.
C
      CALL  MMHASH  (CNAME, ICN1, ICN2, IH, LFOUND)
      IF  (LFOUND)  GO  TO  10
          CALL  MMNEWI  (IS1)
          IHASH(IH)  =  IS1
          CALL  MMPUT1  (CNAME, ICN1, ICN2, ISTORE(IS1), I)
          GO  TO  20
C
C         RETURN THE SPACE ALLOCATED TO THE OLD VALUE
C
   10 CONTINUE
      IS1   =  IHASH(IH)
      CALL  MMDEL1  (ISTORE(IS1+1))
C
C         STORE THE NEW VALUE IN THE TABLE
C
   20 CONTINUE
      CALL  MMPUT1  (CVALUE, ICV1, ICV2, ISTORE(IS1+1), ISTORE(IS1+2))
C
      RETURN
      END
 
      SUBROUTINE  MMPUT1  (CVALUE, ICV1, ICV2, IS2H, IS2T)
C
C-----------------------------------------------------------------------
C
C     FAMILY
C     ------
C     MEMORY MANAGER
C
C     PURPOSE
C     -------
C     TO PUT A STRING VALUE INTO CHARACTER STORAGE
C     AND RETURN POINTERS TO ITS LOCATION
C
C     PARAMETERS
C     ----------
C     CVALUE  -I-  CONTAINS THE CHARACTER STRING
C     ICV1    -I-  INDEX OF THE FIRST CHARACTER IN THE STRING
C     ICV2    -I-  INDEX OF THE LAST CHARACTER IN THE STRING
C     IS2H    -O-  POINTER TO THE FIRST BLOCK CONTAINING THE STRING
C     IS2T    -O-  POINTER TO THE LAST BLOCK CONTAINING THE STRING
C
C     LOCAL VARIABLES
C     ---------------
C     ICV      -   INDEX OF THE CURRENT CHARACTER IN THE STRING
C     IS2      -   POINTER TO CURRENT BLOCK FOR THE STRING
C     ICS      -   INDEX OF CURRENT STORE POSITION IN CSTORE
C     ICS1     -   INDEX OF BEGINNING OF CURRENT CSTORE BLOCK
C     ICS2     -   INDEX OF END OF CURRENT CSTORE BLOCK
C     ICSTST   -   INDEX OF LAST CSTORE POSITION NEEDED
C     ICSMIN   -   INDEX OF LAST CSTORE POSITION NEEDED IN CURRENT BLOCK
C
C-----------------------------------------------------------------------
C
C         GLOBAL CONSTANTS
C
      CHARACTER*1         CA,      CBLANK,  CC,      CI,      CLEFT,
     A                    CMINUS,  CPLUS,   CPOINT,  CQUOTE,  CRIGHT,
     B                    CZ,      C0,      C9
      COMMON  / GLCOMC /  CA,      CBLANK,  CC,      CI,      CLEFT,
     A                    CMINUS,  CPLUS,   CPOINT,  CQUOTE,  CRIGHT,
     B                    CZ,      C0,      C9
C
C         INPUT / OUTPUT CONTROL INTERFACE
C
      CHARACTER*1         CBUFFR(2000)
      LOGICAL             LBREAK,  LFORT,   LISTI,   LISTO
      INTEGER             ICBADD,  ICBEND,  ICBEOL,  ICBSUB,  ICB0,
     A                    ICB1,    ICB2,    ICB3,    ICBDIM,  ICPLI,
     B                    ICPLO,   ILCTR,   ILNMBR,  ILPP,    IPAGE,
     A                    IUNITE,  IUNITI,  IUNITL,  IUNITO
      COMMON  / IOCOMC /  CBUFFR
      COMMON  / IOCOMI /  ICBADD,  ICBEND,  ICBEOL,  ICBSUB,  ICB0,
     A                    ICB1,    ICB2,    ICB3,    ICBDIM,  ICPLI,
     B                    ICPLO,   ILCTR,   ILNMBR,  ILPP,    IPAGE,
     A                    IUNITE,  IUNITI,  IUNITL,  IUNITO
      COMMON  / IOCOML /  LBREAK,  LFORT,   LISTI,   LISTO
C
C         MEMORY MANAGER INTERFACE
C
      CHARACTER*1         CSTORE(20000)
      INTEGER             IHASH(601),       ISTORE(6000)
      INTEGER             ICSDIM,  ICSP1,   ICSP2,   IHADIM,
     A                    ISFREE,  ISTDIM,  IS2HDC,  IS2HDS
      COMMON  / MMCOMC /  CSTORE
      COMMON  / MMCOMH /  IHASH
      COMMON  / MMCOMS /  ISTORE
      COMMON  / MMCOMI /  ICSDIM,  ICSP1,   ICSP2,   IHADIM,
     A                    ISFREE,  ISTDIM,  IS2HDC,  IS2HDS
C
C         LOCAL VARIABLES AND PARAMETERS
C
      INTEGER             ICV1, ICV2, IS2H, IS2T
      CHARACTER*(*)      CVALUE(ICV2)
      INTEGER             ICV, IS2, ICS1, ICS2, ICSTST, ICSMIN, ICS
      EXTERNAL            MMNEWI, IOERRM
C
      IF  ((ICV1 .GT. 0) .AND. (ICV1 .LE. ICV2))  GO  TO  10
          IS2H  =  0
          IS2T  =  0
          GO  TO  999
C
   10 CONTINUE
      ICV   =  ICV1
      IS2   =  IS2HDC
      IS2H  =  IS2HDC
C
C         LOOP THROUGH THE LINKED LIST OF AVAILABLE MEMORY BLOCKS
C
   20 CONTINUE
          IF  (IS2 .EQ. 0)  GO  TO  60
          IS2T    =  IS2
          ICS1    =  ISTORE(IS2T)
          ICS2    =  ISTORE(IS2T+1)
          IS2     =  ISTORE(IS2T+2)
          ICSTST  =  ICS1 + ICV2 - ICV
          ICSMIN  =  MIN0(ICS2, ICSTST)
C
C             STORE CHARACTERS INTO A PARTICULAR BLOCK
C
          DO  30  ICS=ICS1,ICSMIN
              CSTORE(ICS)  =  CVALUE(ICV)
              ICV  =  ICV + 1
   30     CONTINUE
      IF  (ICSTST .GT. ICS2)  GO  TO  20
C
C         IF THE LAST BLOCK USED WAS COMPLETELY FILLED, GO TO 40
C
      IF  (ICSTST .NE. ICS2)  GO  TO  40
          IS2HDC  =  IS2
          GO  TO  50
C
C         THE LAST BLOCK OF MEMORY WAS NOT COMPLETELY USED.
C         PUT A NEW BLOCK ON THE AVAILABLE MEMORY STACK
C         CORRESPONDING TO THE REMAINING CHARACTERS.
C
   40 CONTINUE
      CALL  MMNEWI  (IS2HDC)
      ISTORE(IS2HDC)    =  ICSMIN+1
      ISTORE(IS2HDC+1)  =  ICS2
      ISTORE(IS2HDC+2)  =  IS2
C
   50 CONTINUE
      ISTORE(IS2T+1)    =  ICSTST
      ISTORE(IS2T+2)    =  0
      GO  TO  999
C
C         FATAL ERROR - NO MORE CHARACTER STORAGE SPACE
C
   60 CONTINUE
      CALL  IOERRM  (.TRUE.,
     A '('' ********   MMPUT1 - STORAGE ARRAY CSTORE(*) FULL'')')
C
  999 CONTINUE
      RETURN
      END
      SUBROUTINE  MMRETI  (IS)
C
C-----------------------------------------------------------------------
C
C     FAMILY
C     ------
C     MEMORY MANAGER
C
C     PURPOSE
C     -------
C     TO RETURN AN INTEGER BLOCK TO THE FREE LIST
C
C     PARAMETERS
C     ----------
C     IS      -I-  POINTER TO THE BLOCK TO BE RETURNED
C
C-----------------------------------------------------------------------
C
C         GLOBAL CONSTANTS
C
      CHARACTER*1         CA,      CBLANK,  CC,      CI,      CLEFT,
     A                    CMINUS,  CPLUS,   CPOINT,  CQUOTE,  CRIGHT,
     B                    CZ,      C0,      C9
      COMMON  / GLCOMC /  CA,      CBLANK,  CC,      CI,      CLEFT,
     A                    CMINUS,  CPLUS,   CPOINT,  CQUOTE,  CRIGHT,
     B                    CZ,      C0,      C9
C
C         INPUT / OUTPUT CONTROL INTERFACE
C
      CHARACTER*1         CBUFFR(2000)
      LOGICAL             LBREAK,  LFORT,   LISTI,   LISTO
      INTEGER             ICBADD,  ICBEND,  ICBEOL,  ICBSUB,  ICB0,
     A                    ICB1,    ICB2,    ICB3,    ICBDIM,  ICPLI,
     B                    ICPLO,   ILCTR,   ILNMBR,  ILPP,    IPAGE,
     A                    IUNITE,  IUNITI,  IUNITL,  IUNITO
      COMMON  / IOCOMC /  CBUFFR
      COMMON  / IOCOMI /  ICBADD,  ICBEND,  ICBEOL,  ICBSUB,  ICB0,
     A                    ICB1,    ICB2,    ICB3,    ICBDIM,  ICPLI,
     B                    ICPLO,   ILCTR,   ILNMBR,  ILPP,    IPAGE,
     A                    IUNITE,  IUNITI,  IUNITL,  IUNITO
      COMMON  / IOCOML /  LBREAK,  LFORT,   LISTI,   LISTO
C
C         MEMORY MANAGER INTERFACE
C
      CHARACTER*1         CSTORE(20000)
      INTEGER             IHASH(601),       ISTORE(6000)
      INTEGER             ICSDIM,  ICSP1,   ICSP2,   IHADIM,
     A                    ISFREE,  ISTDIM,  IS2HDC,  IS2HDS
      COMMON  / MMCOMC /  CSTORE
      COMMON  / MMCOMH /  IHASH
      COMMON  / MMCOMS /  ISTORE
      COMMON  / MMCOMI /  ICSDIM,  ICSP1,   ICSP2,   IHADIM,
     A                    ISFREE,  ISTDIM,  IS2HDC,  IS2HDS
C
      INTEGER IS
      EXTERNAL IOERRM
      IF  (IS .EQ. 0)  GO  TO  999
      IF  ((IS .LT. 0) .OR. (IS .GT. ISTDIM)
     A                 .OR. (MOD(IS,3) .NE. 1))  GO  TO  10
      ISTORE(IS+2)  =  ISFREE
      ISFREE        =  IS
      GO  TO  999
C
   10 CONTINUE
      CALL  IOERRM  (.FALSE.,
     A '('' ********   MMRETI - ATTEMPT TO RETURN INVALID POINTER'')')
C
  999 CONTINUE
      RETURN
      END
      SUBROUTINE  MMSETP  (CPTR, ICP1, ICP2)
C
C-----------------------------------------------------------------------
C
C     FAMILY
C     ------
C     MEMORY MANAGER
C
C     PURPOSE
C     -------
C     TO SAVE A POINTER TO THE CURRENT TOP OF THE SUBSTITUTION STACK
C
C     PARAMETERS
C     ----------
C     CPTR    -I-  NAME OF THE POINTER
C     ICP1    -I-  INDEX OF THE FIRST CHARACTER IN THE NAME
C     ICP2    -I-  INDEX OF THE LAST CHARACTER IN THE NAME
C
C-----------------------------------------------------------------------
C
C         GLOBAL CONSTANTS
C
      CHARACTER*1         CA,      CBLANK,  CC,      CI,      CLEFT,
     A                    CMINUS,  CPLUS,   CPOINT,  CQUOTE,  CRIGHT,
     B                    CZ,      C0,      C9
      COMMON  / GLCOMC /  CA,      CBLANK,  CC,      CI,      CLEFT,
     A                    CMINUS,  CPLUS,   CPOINT,  CQUOTE,  CRIGHT,
     B                    CZ,      C0,      C9
C
C         INPUT / OUTPUT CONTROL INTERFACE
C
      CHARACTER*1         CBUFFR(2000)
      LOGICAL             LBREAK,  LFORT,   LISTI,   LISTO
      INTEGER             ICBADD,  ICBEND,  ICBEOL,  ICBSUB,  ICB0,
     A                    ICB1,    ICB2,    ICB3,    ICBDIM,  ICPLI,
     B                    ICPLO,   ILCTR,   ILNMBR,  ILPP,    IPAGE,
     A                    IUNITE,  IUNITI,  IUNITL,  IUNITO
      COMMON  / IOCOMC /  CBUFFR
      COMMON  / IOCOMI /  ICBADD,  ICBEND,  ICBEOL,  ICBSUB,  ICB0,
     A                    ICB1,    ICB2,    ICB3,    ICBDIM,  ICPLI,
     B                    ICPLO,   ILCTR,   ILNMBR,  ILPP,    IPAGE,
     A                    IUNITE,  IUNITI,  IUNITL,  IUNITO
      COMMON  / IOCOML /  LBREAK,  LFORT,   LISTI,   LISTO
C
C         MEMORY MANAGER INTERFACE
C
      CHARACTER*1         CSTORE(20000)
      INTEGER             IHASH(601),       ISTORE(6000)
      INTEGER             ICSDIM,  ICSP1,   ICSP2,   IHADIM,
     A                    ISFREE,  ISTDIM,  IS2HDC,  IS2HDS
      COMMON  / MMCOMC /  CSTORE
      COMMON  / MMCOMH /  IHASH
      COMMON  / MMCOMS /  ISTORE
      COMMON  / MMCOMI /  ICSDIM,  ICSP1,   ICSP2,   IHADIM,
     A                    ISFREE,  ISTDIM,  IS2HDC,  IS2HDS
C
C         LOCAL VARIABLES AND PARAMETERS
C
      INTEGER             ICP1, ICP2
      CHARACTER*(*)      CPTR(ICP2)
      LOGICAL             LFOUND
      INTEGER             I, IH, IS1, IS2
      EXTERNAL            MMHASH, MMNEWI, MMPUT1
C
      CALL  MMHASH  (CPTR, ICP1, ICP2, IH, LFOUND)
      IF  (.NOT. LFOUND)  GO  TO  10
      IS1            =  IHASH(IH)
      IS2            =  ISTORE(IS1+1)
      GO  TO  20
C
   10 CONTINUE
      CALL  MMNEWI  (IS1)
      IHASH(IH)      =  IS1
      CALL  MMPUT1  (CPTR, ICP1, ICP2, ISTORE(IS1), I)
      CALL  MMNEWI  (IS2)
      ISTORE(IS1+1)  =  IS2
      ISTORE(IS1+2)  =  0
C
   20 CONTINUE
      IF  (IS2HDS .NE. 0)  GO  TO  30
      ISTORE(IS2)    =  0
      ISTORE(IS2+1)  =  0
      ISTORE(IS2+2)  =  0
      GO  TO  999
C
   30 CONTINUE
      ISTORE(IS2)    =  ISTORE(IS2HDS)
      ISTORE(IS2+1)  =  ISTORE(IS2HDS+1)
      ISTORE(IS2+2)  =  0
C
  999 CONTINUE
      RETURN
      END
      SUBROUTINE  MMTEST  (CVALUE, ICV1, ICV2, IS2H, LEQUAL)
C
C-----------------------------------------------------------------------
C
C     FAMILY
C     ------
C     MEMORY MANAGER
C
C     PURPOSE
C     -------
C     TO SEE IF A GIVEN STRING IS EQUAL TO ONE IN THE SYMBOL TABLE
C
C     PARAMETERS
C     ----------
C     CVALUE  -I-  CONTAINS THE STRING TO BE TESTED
C     ICV1    -I-  INDEX OF THE FIRST CHARACTER IN THE STRING
C     ICV2    -I-  INDEX OF THE LAST CHARACTER IN THE STRING
C     IS2H    -I-  POINTER TO THE STRING IN THE SYMBOL TABLE
C     LEQUAL  -O-  TRUE IF THE STRINGS ARE EQUAL
C
C-----------------------------------------------------------------------
C
C         GLOBAL CONSTANTS
C
      CHARACTER*1         CA,      CBLANK,  CC,      CI,      CLEFT,
     A                    CMINUS,  CPLUS,   CPOINT,  CQUOTE,  CRIGHT,
     B                    CZ,      C0,      C9
      COMMON  / GLCOMC /  CA,      CBLANK,  CC,      CI,      CLEFT,
     A                    CMINUS,  CPLUS,   CPOINT,  CQUOTE,  CRIGHT,
     B                    CZ,      C0,      C9
C
C         INPUT / OUTPUT CONTROL INTERFACE
C
      CHARACTER*1         CBUFFR(2000)
      LOGICAL             LBREAK,  LFORT,   LISTI,   LISTO
      INTEGER             ICBADD,  ICBEND,  ICBEOL,  ICBSUB,  ICB0,
     A                    ICB1,    ICB2,    ICB3,    ICBDIM,  ICPLI,
     B                    ICPLO,   ILCTR,   ILNMBR,  ILPP,    IPAGE,
     A                    IUNITE,  IUNITI,  IUNITL,  IUNITO
      COMMON  / IOCOMC /  CBUFFR
      COMMON  / IOCOMI /  ICBADD,  ICBEND,  ICBEOL,  ICBSUB,  ICB0,
     A                    ICB1,    ICB2,    ICB3,    ICBDIM,  ICPLI,
     B                    ICPLO,   ILCTR,   ILNMBR,  ILPP,    IPAGE,
     A                    IUNITE,  IUNITI,  IUNITL,  IUNITO
      COMMON  / IOCOML /  LBREAK,  LFORT,   LISTI,   LISTO
C
C         MEMORY MANAGER INTERFACE
C
      CHARACTER*1         CSTORE(20000)
      INTEGER             IHASH(601),       ISTORE(6000)
      INTEGER             ICSDIM,  ICSP1,   ICSP2,   IHADIM,
     A                    ISFREE,  ISTDIM,  IS2HDC,  IS2HDS
      COMMON  / MMCOMC /  CSTORE
      COMMON  / MMCOMH /  IHASH
      COMMON  / MMCOMS /  ISTORE
      COMMON  / MMCOMI /  ICSDIM,  ICSP1,   ICSP2,   IHADIM,
     A                    ISFREE,  ISTDIM,  IS2HDC,  IS2HDS
C
C         LOCAL VARIABLES AND PARAMETERS
C
      INTEGER             ICV1, ICV2, IS2H
      CHARACTER*(*)      CVALUE(ICV2)
      LOGICAL             LEQUAL
      INTEGER             ICS, ICS1, ICS2, ICV, IS2
C
      ICV     =  ICV1
      IS2     =  IS2H
      LEQUAL  =  .FALSE.
C
   10 CONTINUE
          IF  (IS2 .EQ. 0)  GO  TO  30
          ICS1  =  ISTORE(IS2)
          ICS2  =  ISTORE(IS2+1)
          IS2   =  ISTORE(IS2+2)
          IF  (ICS2-ICS1 .GT. ICV2-ICV)  GO  TO  999
          DO  20  ICS=ICS1,ICS2
              IF  (CSTORE(ICS) .NE. CVALUE(ICV))  GO  TO  999
              ICV  =  ICV + 1
   20     CONTINUE
      GO  TO  10
C
   30 CONTINUE
      LEQUAL  =  ICV .GT. ICV2
C
  999 CONTINUE
      RETURN
      END
      SUBROUTINE  MPEOL
C
C-----------------------------------------------------------------------
C
C     FAMILY
C     ------
C     MACRO PROCESSOR
C
C     PURPOSE
C     -------
C     TO REMOVE TRAILING BLANKS AND ADD AN END-OF-LINE MARKER
C     TO THE LINE IN THE INPUT/OUTPUT BUFFER
C
C-----------------------------------------------------------------------
C
C         GLOBAL CONSTANTS
C
      CHARACTER*1         CA,      CBLANK,  CC,      CI,      CLEFT,
     A                    CMINUS,  CPLUS,   CPOINT,  CQUOTE,  CRIGHT,
     B                    CZ,      C0,      C9
      COMMON  / GLCOMC /  CA,      CBLANK,  CC,      CI,      CLEFT,
     A                    CMINUS,  CPLUS,   CPOINT,  CQUOTE,  CRIGHT,
     B                    CZ,      C0,      C9
C
C         INPUT / OUTPUT CONTROL INTERFACE
C
      CHARACTER*1         CBUFFR(2000)
      LOGICAL             LBREAK,  LFORT,   LISTI,   LISTO
      INTEGER             ICBADD,  ICBEND,  ICBEOL,  ICBSUB,  ICB0,
     A                    ICB1,    ICB2,    ICB3,    ICBDIM,  ICPLI,
     B                    ICPLO,   ILCTR,   ILNMBR,  ILPP,    IPAGE,
     A                    IUNITE,  IUNITI,  IUNITL,  IUNITO
      COMMON  / IOCOMC /  CBUFFR
      COMMON  / IOCOMI /  ICBADD,  ICBEND,  ICBEOL,  ICBSUB,  ICB0,
     A                    ICB1,    ICB2,    ICB3,    ICBDIM,  ICPLI,
     B                    ICPLO,   ILCTR,   ILNMBR,  ILPP,    IPAGE,
     A                    IUNITE,  IUNITI,  IUNITL,  IUNITO
      COMMON  / IOCOML /  LBREAK,  LFORT,   LISTI,   LISTO
C
C         MEMORY MANAGER INTERFACE
C
      CHARACTER*1         CSTORE(20000)
      INTEGER             IHASH(601),       ISTORE(6000)
      INTEGER             ICSDIM,  ICSP1,   ICSP2,   IHADIM,
     A                    ISFREE,  ISTDIM,  IS2HDC,  IS2HDS
      COMMON  / MMCOMC /  CSTORE
      COMMON  / MMCOMH /  IHASH
      COMMON  / MMCOMS /  ISTORE
      COMMON  / MMCOMI /  ICSDIM,  ICSP1,   ICSP2,   IHADIM,
     A                    ISFREE,  ISTDIM,  IS2HDC,  IS2HDS
C
C         MACRO PROCESSOR INTERFACE
C
      CHARACTER*1         CDIV,    CEOL,     CEOR,    CONC,    CSUB,
     A                    CTOP
      LOGICAL             LEMPTY,  LSUB
      COMMON  / MPCOMC /  CDIV,    CEOL,    CEOR,    CONC,    CSUB,
     A                    CTOP
      COMMON  / MPCOML /  LEMPTY,  LSUB
C
C         LOCAL VARIABLES AND PARAMETERS
C
      CHARACTER*1         C
      EXTERNAL            MMPOPC
C
      CALL  MMPOPC  (C, 2, CTOP, LEMPTY)
      ICBEOL  =  ICB2
      ICBEND  =  ICBEOL
      ICB2    =  ICB2 - 2
      ICB3    =  ICB2
      IF  (ICB1 .GT. ICB2)            GO  TO  999
      IF  (CBUFFR(ICB2) .NE. CBLANK)  GO  TO  999
C
C         REMOVE TRAILING BLANKS
C
   10 CONTINUE
          ICB2  =  ICB2 - 1
          IF  (ICB1 .GT. ICB2)        GO  TO  20
      IF  (CBUFFR(ICB2) .EQ. CBLANK)  GO  TO  10
C
C         ADD THE END-OF-LINE MARKER
C
   20 CONTINUE
      CBUFFR(ICB2+1)  =  CSUB
      CBUFFR(ICB2+2)  =  CEOL
      ICB3            =  ICB2
      ICBEOL          =  ICB2 + 2
      ICBEND          =  ICBEOL
C
  999 CONTINUE
      RETURN
      END
      SUBROUTINE  MPITEM
C
C-----------------------------------------------------------------------
C
C     FAMILY
C     ------
C     MACRO PROCESSOR
C
C     PURPOSE
C     -------
C     TO PUSH THE NEXT ITEM IN A LIST ONTO THE SUBSTITUTION STACK
C
C-----------------------------------------------------------------------
C
C         GLOBAL CONSTANTS
C
      CHARACTER*1         CA,      CBLANK,  CC,      CI,      CLEFT,
     A                    CMINUS,  CPLUS,   CPOINT,  CQUOTE,  CRIGHT,
     B                    CZ,      C0,      C9
      COMMON  / GLCOMC /  CA,      CBLANK,  CC,      CI,      CLEFT,
     A                    CMINUS,  CPLUS,   CPOINT,  CQUOTE,  CRIGHT,
     B                    CZ,      C0,      C9
C
C         INPUT / OUTPUT CONTROL INTERFACE
C
      CHARACTER*1         CBUFFR(2000)
      LOGICAL             LBREAK,  LFORT,   LISTI,   LISTO
      INTEGER             ICBADD,  ICBEND,  ICBEOL,  ICBSUB,  ICB0,
     A                    ICB1,    ICB2,    ICB3,    ICBDIM,  ICPLI,
     B                    ICPLO,   ILCTR,   ILNMBR,  ILPP,    IPAGE,
     A                    IUNITE,  IUNITI,  IUNITL,  IUNITO
      COMMON  / IOCOMC /  CBUFFR
      COMMON  / IOCOMI /  ICBADD,  ICBEND,  ICBEOL,  ICBSUB,  ICB0,
     A                    ICB1,    ICB2,    ICB3,    ICBDIM,  ICPLI,
     B                    ICPLO,   ILCTR,   ILNMBR,  ILPP,    IPAGE,
     A                    IUNITE,  IUNITI,  IUNITL,  IUNITO
      COMMON  / IOCOML /  LBREAK,  LFORT,   LISTI,   LISTO
C
C         MEMORY MANAGER INTERFACE
C
      CHARACTER*1         CSTORE(20000)
      INTEGER             IHASH(601),       ISTORE(6000)
      INTEGER             ICSDIM,  ICSP1,   ICSP2,   IHADIM,
     A                    ISFREE,  ISTDIM,  IS2HDC,  IS2HDS
      COMMON  / MMCOMC /  CSTORE
      COMMON  / MMCOMH /  IHASH
      COMMON  / MMCOMS /  ISTORE
      COMMON  / MMCOMI /  ICSDIM,  ICSP1,   ICSP2,   IHADIM,
     A                    ISFREE,  ISTDIM,  IS2HDC,  IS2HDS
C
C         MACRO PROCESSOR INTERFACE
C
      CHARACTER*1         CDIV,    CEOL,     CEOR,    CONC,    CSUB,
     A                    CTOP
      LOGICAL             LEMPTY,  LSUB
      COMMON  / MPCOMC /  CDIV,    CEOL,    CEOR,    CONC,    CSUB,
     A                    CTOP
      COMMON  / MPCOML /  LEMPTY,  LSUB
C
C         LOCAL VARIABLES AND PARAMETERS
C
      LOGICAL             LERROR, LFOUND
      INTEGER             ICN1, ICN2, ICP1, ICP2
      EXTERNAL            MPPOPN, UTBLDN, MMPSHV, MMPUTP, IOERRM
C
      CALL  MPPOPN  (ICN1, ICN2, LERROR)
      ICP1  =  ICB2 + 1
      CALL  UTBLDN  (CDIV, CBUFFR, ICN1, ICN2, 1,
     A               CBUFFR, ICP1, ICP2, ICBDIM, LERROR)
      CALL  MMPSHV  (CBUFFR, ICP1, ICP2, 3, LEMPTY, LFOUND)
      IF  (LFOUND)  GO  TO  10
      CALL  MMPUTP  (CBUFFR, ICN1, ICN2, CBUFFR, ICP1, ICP2, LFOUND)
      IF  (.NOT. LFOUND)  GO  TO  20
      CALL  MMPSHV  (CBUFFR, ICP1, ICP2, 3, LEMPTY, LFOUND)
C
   10 CONTINUE
      ICB2  =  ICBSUB - 1
      GO  TO  999
C
   20 CONTINUE
      CALL  IOERRM  (.FALSE.,
     A '('' ********   MPITEM - VARIABLE NOT DEFINED'')')
C
  999 CONTINUE
      RETURN
      END
      SUBROUTINE  MPLABL
C
C-----------------------------------------------------------------------
C
C     FAMILY
C     ------
C     MACRO PROCESSOR
C
C     PURPOSE
C     -------
C     TO COPY THE CURRENT LABEL TO THE BUFFER AND THEN
C     INCREMENT AND SAVE ITS VALUE
C
C-----------------------------------------------------------------------
C
C         GLOBAL CONSTANTS
C
      CHARACTER*1         CA,      CBLANK,  CC,      CI,      CLEFT,
     A                    CMINUS,  CPLUS,   CPOINT,  CQUOTE,  CRIGHT,
     B                    CZ,      C0,      C9
      COMMON  / GLCOMC /  CA,      CBLANK,  CC,      CI,      CLEFT,
     A                    CMINUS,  CPLUS,   CPOINT,  CQUOTE,  CRIGHT,
     B                    CZ,      C0,      C9
C
C         INPUT / OUTPUT CONTROL INTERFACE
C
      CHARACTER*1         CBUFFR(2000)
      LOGICAL             LBREAK,  LFORT,   LISTI,   LISTO
      INTEGER             ICBADD,  ICBEND,  ICBEOL,  ICBSUB,  ICB0,
     A                    ICB1,    ICB2,    ICB3,    ICBDIM,  ICPLI,
     B                    ICPLO,   ILCTR,   ILNMBR,  ILPP,    IPAGE,
     A                    IUNITE,  IUNITI,  IUNITL,  IUNITO
      COMMON  / IOCOMC /  CBUFFR
      COMMON  / IOCOMI /  ICBADD,  ICBEND,  ICBEOL,  ICBSUB,  ICB0,
     A                    ICB1,    ICB2,    ICB3,    ICBDIM,  ICPLI,
     B                    ICPLO,   ILCTR,   ILNMBR,  ILPP,    IPAGE,
     A                    IUNITE,  IUNITI,  IUNITL,  IUNITO
      COMMON  / IOCOML /  LBREAK,  LFORT,   LISTI,   LISTO
C
C         MEMORY MANAGER INTERFACE
C
      CHARACTER*1         CSTORE(20000)
      INTEGER             IHASH(601),       ISTORE(6000)
      INTEGER             ICSDIM,  ICSP1,   ICSP2,   IHADIM,
     A                    ISFREE,  ISTDIM,  IS2HDC,  IS2HDS
      COMMON  / MMCOMC /  CSTORE
      COMMON  / MMCOMH /  IHASH
      COMMON  / MMCOMS /  ISTORE
      COMMON  / MMCOMI /  ICSDIM,  ICSP1,   ICSP2,   IHADIM,
     A                    ISFREE,  ISTDIM,  IS2HDC,  IS2HDS
C
C         MACRO PROCESSOR INTERFACE
C
      CHARACTER*1         CDIV,    CEOL,     CEOR,    CONC,    CSUB,
     A                    CTOP
      LOGICAL             LEMPTY,  LSUB
      COMMON  / MPCOMC /  CDIV,    CEOL,    CEOR,    CONC,    CSUB,
     A                    CTOP
      COMMON  / MPCOML /  LEMPTY,  LSUB
C
C         LOCAL VARIABLES AND PARAMETERS
C
      INTEGER             I, ICV1, ICV2, ICV2M4, ICV1M1
      CHARACTER*1         C(5)
      LOGICAL             L
      INTEGER             ICNDIM
      EXTERNAL            MMGETV, UTCVCI, UTCVIC, MMPUTV, IOERRM
      SAVE ICNDIM, C
      DATA                ICNDIM / 5 /
      DATA                C(1), C(2), C(3), C(4), C(5)
     A                 /  'L',  'A',  'B',  'E',  'L'  /
C
C         GET CURRENT VALUE OF LABEL, CONVERT TO AN INTEGER AND
C         CHECK IT.  ADD ONE, CONVERT BACK TO CHARACTERS AND PLACE
C         IN CBUFFR
C
      ICV1  =  ICB2 + 5
      CALL  MMGETV  (C, 1, ICNDIM, CBUFFR, ICV1, ICV2, ICBDIM, L)
      IF  (.NOT. L)  GO  TO  40
      CALL  UTCVCI  (CBUFFR, ICV1, ICV2, I, L)
      IF  ((I .LT. 0) .OR. (99999 .LT. I))  GO  TO  40
      I       =  I + 1
      CALL  UTCVIC  (CBUFFR, ICV1, ICV2, ICBDIM, I, L)
      ICV2M4  =  ICV2 - 4
      ICV1M1  =  ICV1 - 1
      IF  (ICV2M4 .GT. ICV1M1)  GO  TO  20
C
      DO  10  I=ICV2M4,ICV1M1
          CBUFFR(I)  =  CBLANK
   10 CONTINUE
C
   20 CONTINUE
      CALL  MMPUTV  (C, 1, ICNDIM, CBUFFR, ICV2M4, ICV2)
      ICB2  =  ICBSUB - 1
      DO  30  I=ICV2M4,ICV2
          ICB2  =  ICB2 + 1
          CBUFFR(ICB2)  =  CBUFFR(I)
   30 CONTINUE
      GO  TO  999
C
C         WARNING - INVALID LABEL VALUE
C
   40 CONTINUE
      CALL  IOERRM  (.FALSE.,
     A '('' ********   MPLABL - ILLEGAL LABEL VALUE'')')
C
  999 CONTINUE
      RETURN
      END
 
      SUBROUTINE  MPLINE  (LSUBL)
C
C-----------------------------------------------------------------------
C
C     FAMILY
C     ------
C     MACRO PROCESSOR
C
C     PURPOSE
C     -------
C     TO BUILD THE NEXT LINE IN THE I/O BUFFER
C
C     PARAMETERS
C     ----------
C     LSUBL   -I-  LOCAL SUBSTITUTION FLAG INDICATING WHETHER
C                  OR NOT MACROS ON THIS LINE ARE TO BE EXPANDED
C
C-----------------------------------------------------------------------
C
C         GLOBAL CONSTANTS
C
      CHARACTER*1         CA,      CBLANK,  CC,      CI,      CLEFT,
     A                    CMINUS,  CPLUS,   CPOINT,  CQUOTE,  CRIGHT,
     B                    CZ,      C0,      C9
      COMMON  / GLCOMC /  CA,      CBLANK,  CC,      CI,      CLEFT,
     A                    CMINUS,  CPLUS,   CPOINT,  CQUOTE,  CRIGHT,
     B                    CZ,      C0,      C9
C
C         INPUT / OUTPUT CONTROL INTERFACE
C
      CHARACTER*1         CBUFFR(2000)
      LOGICAL             LBREAK,  LFORT,   LISTI,   LISTO
      INTEGER             ICBADD,  ICBEND,  ICBEOL,  ICBSUB,  ICB0,
     A                    ICB1,    ICB2,    ICB3,    ICBDIM,  ICPLI,
     B                    ICPLO,   ILCTR,   ILNMBR,  ILPP,    IPAGE,
     A                    IUNITE,  IUNITI,  IUNITL,  IUNITO
      COMMON  / IOCOMC /  CBUFFR
      COMMON  / IOCOMI /  ICBADD,  ICBEND,  ICBEOL,  ICBSUB,  ICB0,
     A                    ICB1,    ICB2,    ICB3,    ICBDIM,  ICPLI,
     B                    ICPLO,   ILCTR,   ILNMBR,  ILPP,    IPAGE,
     A                    IUNITE,  IUNITI,  IUNITL,  IUNITO
      COMMON  / IOCOML /  LBREAK,  LFORT,   LISTI,   LISTO
C
C         MEMORY MANAGER INTERFACE
C
      CHARACTER*1         CSTORE(20000)
      INTEGER             IHASH(601),       ISTORE(6000)
      INTEGER             ICSDIM,  ICSP1,   ICSP2,   IHADIM,
     A                    ISFREE,  ISTDIM,  IS2HDC,  IS2HDS
      COMMON  / MMCOMC /  CSTORE
      COMMON  / MMCOMH /  IHASH
      COMMON  / MMCOMS /  ISTORE
      COMMON  / MMCOMI /  ICSDIM,  ICSP1,   ICSP2,   IHADIM,
     A                    ISFREE,  ISTDIM,  IS2HDC,  IS2HDS
C
C         MACRO PROCESSOR INTERFACE
C
      CHARACTER*1         CDIV,    CEOL,     CEOR,    CONC,    CSUB,
     A                    CTOP
      LOGICAL             LEMPTY,  LSUB
      COMMON  / MPCOMC /  CDIV,    CEOL,    CEOR,    CONC,    CSUB,
     A                    CTOP
      COMMON  / MPCOML /  LEMPTY,  LSUB
C
C         LOCAL VARIABLES AND PARAMETERS
C
      INTEGER             ICB
      CHARACTER*1         C(2)
      LOGICAL             LEOL, LFOUND, LSUBL
      INTEGER             ICNDIM
      EXTERNAL            MMPOPC, MPSUBS, IOREAD, MMPUTV, MMPSHV
      SAVE  ICNDIM, C
      DATA                ICNDIM / 2 /
      DATA                C(1), C(2)  /'*', 'L'  /
C
C         SET I/O BUFFER POINTERS
C
      ICB0  =  ICBEOL + 1
      ICB1  =  ICBEOL + 1
      ICB2  =  ICBEOL
      IF  (LEMPTY)  GO  TO  20
C
C         IF THE STACK IS NONEMPTY, POP INTO THE I/O BUFFER UNTIL A
C         SUB. CHAR. IS FOUND. THEN CALL MPSUBS.
C
   10 CONTINUE
          CALL  MMPOPC  (CSUB, 5, CTOP, LEMPTY)
          IF  (LEMPTY)  GO  TO  20
          CALL  MPSUBS  (LEOL, LSUBL)
          IF  (LEOL)  GO  TO  999
      GO  TO  10
C
C         IF THE STACK IS EMPTY, GET MORE INPUT
C
   20 CONTINUE
      CALL  IOREAD
      IF  ((.NOT. LSUBL) .OR. (ICB1 .GT. ICB2))  GO  TO  999
C
C         LOOK FOR SUBSTITUTION CHARACTER
C
      DO  30  ICB=ICB1,ICB2
          IF  (CBUFFR(ICB) .EQ. CSUB)  GO  TO  40
   30 CONTINUE
      GO  TO  999
C
C         WHEN A SUB. CHAR IS FOUND, PUT THE VARIABLE '*L' IN THE
C         SYMBOL TABLE.  THE VALUE OF THIS SPECIAL VARIABLE IS THE
C         REST OF THE LINE.  ALSO PUSH *L ONTO THE SUBST. STACK.
C
   40 CONTINUE
      CALL  MMPUTV  (C, 1, ICNDIM, CBUFFR, ICB, ICBEOL)
      CALL  MMPSHV  (C, 1, ICNDIM, 1, LEMPTY, LFOUND)
      ICB2  =  ICB - 1
      CALL  MPSUBS  (LEOL, LSUBL)
      IF  (.NOT. LEOL)  GO  TO  10
C
  999 CONTINUE
      RETURN
      END
      SUBROUTINE  MPMAC
C
C-----------------------------------------------------------------------
C
C     FAMILY
C     ------
C     MACRO PROCESSOR
C
C     PURPOSE
C     -------
C     TO DETERMINE THE TYPE OF MACRO EXPANSION INDICATED
C     AND CALL THE APPROPRIATE ROUTINES
C
C-----------------------------------------------------------------------
C
C         GLOBAL CONSTANTS
C
      CHARACTER*1         CA,      CBLANK,  CC,      CI,      CLEFT,
     A                    CMINUS,  CPLUS,   CPOINT,  CQUOTE,  CRIGHT,
     B                    CZ,      C0,      C9
      COMMON  / GLCOMC /  CA,      CBLANK,  CC,      CI,      CLEFT,
     A                    CMINUS,  CPLUS,   CPOINT,  CQUOTE,  CRIGHT,
     B                    CZ,      C0,      C9
C
C         INPUT / OUTPUT CONTROL INTERFACE
C
      CHARACTER*1         CBUFFR(2000)
      LOGICAL             LBREAK,  LFORT,   LISTI,   LISTO
      INTEGER             ICBADD,  ICBEND,  ICBEOL,  ICBSUB,  ICB0,
     A                    ICB1,    ICB2,    ICB3,    ICBDIM,  ICPLI,
     B                    ICPLO,   ILCTR,   ILNMBR,  ILPP,    IPAGE,
     A                    IUNITE,  IUNITI,  IUNITL,  IUNITO
      COMMON  / IOCOMC /  CBUFFR
      COMMON  / IOCOMI /  ICBADD,  ICBEND,  ICBEOL,  ICBSUB,  ICB0,
     A                    ICB1,    ICB2,    ICB3,    ICBDIM,  ICPLI,
     B                    ICPLO,   ILCTR,   ILNMBR,  ILPP,    IPAGE,
     A                    IUNITE,  IUNITI,  IUNITL,  IUNITO
      COMMON  / IOCOML /  LBREAK,  LFORT,   LISTI,   LISTO
C
C         MEMORY MANAGER INTERFACE
C
      CHARACTER*1         CSTORE(20000)
      INTEGER             IHASH(601),       ISTORE(6000)
      INTEGER             ICSDIM,  ICSP1,   ICSP2,   IHADIM,
     A                    ISFREE,  ISTDIM,  IS2HDC,  IS2HDS
      COMMON  / MMCOMC /  CSTORE
      COMMON  / MMCOMH /  IHASH
      COMMON  / MMCOMS /  ISTORE
      COMMON  / MMCOMI /  ICSDIM,  ICSP1,   ICSP2,   IHADIM,
     A                    ISFREE,  ISTDIM,  IS2HDC,  IS2HDS
C
C         MACRO PROCESSOR INTERFACE
C
      CHARACTER*1         CDIV,    CEOL,     CEOR,    CONC,    CSUB,
     A                    CTOP
      LOGICAL             LEMPTY,  LSUB
      COMMON  / MPCOMC /  CDIV,    CEOL,    CEOR,    CONC,    CSUB,
     A                    CTOP
      COMMON  / MPCOML /  LEMPTY,  LSUB
C
C         LOCAL VARIABLES AND PARAMETERS
C
      CHARACTER*1         C(12)
      INTEGER             IK(3)
      INTEGER             IKYDIM, ICKDIM, ICN1, ICN2, ICN1SV, I,
     A                    ICL1NA, ICL2NA, IH
      LOGICAL             LERROR, LFOUND
      EXTERNAL            MMHASH, MPPOPN, UTRDKY, UTRDNA, MMPSHV,
     A                    UTCVLC, MPLABL, MPITEM, IOERRM
      SAVE IKYDIM, IK, ICKDIM, C
      DATA
     A     IKYDIM, ICKDIM,
     B     IK(1),
     C     IK(2),
     D     IK(3)
     E  /  3,      12,
     F     3,
     G     5,
     H     4  /
      DATA
     B     C(1),  C(2),  C(3),
     C     C(4),  C(5),  C(6),  C(7),  C(8),
     D     C(9),  C(10), C(11), C(12)
     E  /
     F     'D',   'E',   'F',
     G     'L',   'A',   'B',   'E',   'L',
     H     'L',   'I',   'S',   'T'  /
C
C         MARK THE CURRENT LINE POSITION, POP THE NAME OFF THE STACK,
C         AND THEN DETERMINE THE TYPE OF SUBSTITUTION
C
      ICBSUB  =  ICB2
      CALL  MPPOPN  (ICN1, ICN2, LERROR)
      IF  (LERROR)  GO  TO  999
      ICN1SV = ICN1
      CALL  UTRDKY  (CBUFFR, ICN1, ICN2, IK, IKYDIM, C, ICKDIM, I)
      IF(I.GT.IKYDIM) GO TO 10
C
C CHECK IF WE ONLY HAPPENED TO MATCH A PREFIX
C
      CALL UTRDNA (CBUFFR, ICN1, ICN2, ICL1NA, ICL2NA, LERROR)
      IF(LERROR)
     1       GO  TO  (20, 30,  40),  I
C
C IF SO, RESTORE ORIGINAL STATE AND PROCESS VARIABLE SUBSTITUTION
C
      ICN1 = ICN1SV
      GO TO 10
C
C         A SIMPLE MACRO SUBSTITUTION HAS BEEN FOUND.
C         PUSH THE NEW NAME ONTO THE STACK
C         AND RESET THE CURRENT LINE POINTER.
C         SUBSEQUENT POPPING OF THE STACK WILL PUT THE VALUE
C         OF THE MACRO INTO THE I/O BUFFER.
C
   10 CONTINUE
      CALL  MMPSHV  (CBUFFR, ICN1, ICN2, 1, LEMPTY, LFOUND)
      IF  (.NOT. LFOUND)  GO  TO  50
      ICB2    =  ICBSUB - 1
      GO  TO  999
C
C         A DEF SUBSTITUTION HAS BEEN ENCOUNTERED
C
   20 CONTINUE
      CALL  MPPOPN  (ICN1, ICN2, LERROR)
      CALL  MMHASH  (CBUFFR, ICN1, ICN2, IH, LFOUND)
      CALL  UTCVLC  (CBUFFR, ICBSUB, ICB2, ICBDIM, LFOUND, LERROR)
      GO  TO  999
C
C         A LABEL SUBSTITUTION HAS BEEN ENCOUNTERED
C
   30 CONTINUE
      CALL  MPLABL
      GO  TO  999
C
C         A LIST SUBSTITUTION HAS POSSIBLY BEEN ENCOUNTERED
C
   40 CONTINUE
      CALL  MPITEM
      GO  TO  999
C
C         WARNING - NAME NOT FOUND IN SYMBOL TABLE
C
   50 CONTINUE
      CALL  IOERRM  (.FALSE.,
     A '('' ********   MPMAC  - VARIABLE NOT DEFINED'')')
C
  999 CONTINUE
      RETURN
      END
      SUBROUTINE  MPPOPN  (ICN1, ICN2, LERROR)
C
C-----------------------------------------------------------------------
C
C     FAMILY
C     ------
C     MACRO PROCESSOR
C
C     PURPOSE
C     -------
C     TO POP A NAME OFF THE SUBSTITUTION STACK INTO THE I/O BUFFER
C
C     PARAMETERS
C     ----------
C     ICN1    -O-  INDEX IN THE BUFFER OF THE FIRST CHARACTER
C                  IN THE NAME
C     ICN2    -O-  INDEX OF THE LAST CHARACTER
C     LERROR  -O-  TRUE IF THE NAME WAS INVALID
C
C-----------------------------------------------------------------------
C
C         GLOBAL CONSTANTS
C
      CHARACTER*1         CA,      CBLANK,  CC,      CI,      CLEFT,
     A                    CMINUS,  CPLUS,   CPOINT,  CQUOTE,  CRIGHT,
     B                    CZ,      C0,      C9
      COMMON  / GLCOMC /  CA,      CBLANK,  CC,      CI,      CLEFT,
     A                    CMINUS,  CPLUS,   CPOINT,  CQUOTE,  CRIGHT,
     B                    CZ,      C0,      C9
C
C         INPUT / OUTPUT CONTROL INTERFACE
C
      CHARACTER*1         CBUFFR(2000)
      LOGICAL             LBREAK,  LFORT,   LISTI,   LISTO
      INTEGER             ICBADD,  ICBEND,  ICBEOL,  ICBSUB,  ICB0,
     A                    ICB1,    ICB2,    ICB3,    ICBDIM,  ICPLI,
     B                    ICPLO,   ILCTR,   ILNMBR,  ILPP,    IPAGE,
     A                    IUNITE,  IUNITI,  IUNITL,  IUNITO
      COMMON  / IOCOMC /  CBUFFR
      COMMON  / IOCOMI /  ICBADD,  ICBEND,  ICBEOL,  ICBSUB,  ICB0,
     A                    ICB1,    ICB2,    ICB3,    ICBDIM,  ICPLI,
     B                    ICPLO,   ILCTR,   ILNMBR,  ILPP,    IPAGE,
     A                    IUNITE,  IUNITI,  IUNITL,  IUNITO
      COMMON  / IOCOML /  LBREAK,  LFORT,   LISTI,   LISTO
C
C         MEMORY MANAGER INTERFACE
C
      CHARACTER*1         CSTORE(20000)
      INTEGER             IHASH(601),       ISTORE(6000)
      INTEGER             ICSDIM,  ICSP1,   ICSP2,   IHADIM,
     A                    ISFREE,  ISTDIM,  IS2HDC,  IS2HDS
      COMMON  / MMCOMC /  CSTORE
      COMMON  / MMCOMH /  IHASH
      COMMON  / MMCOMS /  ISTORE
      COMMON  / MMCOMI /  ICSDIM,  ICSP1,   ICSP2,   IHADIM,
     A                    ISFREE,  ISTDIM,  IS2HDC,  IS2HDS
C
C         MACRO PROCESSOR INTERFACE
C
      CHARACTER*1         CDIV,    CEOL,     CEOR,    CONC,    CSUB,
     A                    CTOP
      LOGICAL             LEMPTY,  LSUB
      COMMON  / MPCOMC /  CDIV,    CEOL,    CEOR,    CONC,    CSUB,
     A                    CTOP
      COMMON  / MPCOML /  LEMPTY,  LSUB
C
C         LOCAL VARIABLES AND PARAMETERS
C
      INTEGER             ICN1, ICN2
      CHARACTER*1         C
      LOGICAL             LEFT, LERROR
      EXTERNAL            MMPOPC, IOERRM
C
C         POP BLANKS; LOOK FOR LEFT PAREN.
C
      CALL  MMPOPC  (CBLANK, 4, CTOP, LEMPTY)
      LEFT  =  CLEFT .EQ. CTOP
      IF  (.NOT. LEFT)  GO  TO  10
          CALL  MMPOPC  (C, 2, CTOP, LEMPTY)
          CALL  MMPOPC  (CBLANK, 4, CTOP, LEMPTY)
C
   10 CONTINUE
      ICN1    =  ICB2 + 1
C
C         CHECK FOR A LEGAL NAME
C
      LERROR  =  .NOT. (LLE(CA,CTOP)
     A                  .AND. LLE(CTOP,CZ))
      IF  (LERROR)  GO  TO  20
C
C         POP THE CHAR'S OF THE NAME OFF
C
      CALL  MMPOPC  (C, 6, CTOP, LEMPTY)
      ICN2  =  ICB2
      IF  (.NOT. LEFT)  GO  TO  999
          CALL  MMPOPC  (CBLANK, 4, CTOP, LEMPTY)
          LERROR  =  CRIGHT .NE. CTOP
          IF  (LERROR)  GO  TO  30
          CALL  MMPOPC  (C, 2, CTOP, LEMPTY)
          GO  TO  999
C
C         WARNING - ILLEGAL NAME
C
   20 CONTINUE
      CALL  IOERRM  (.FALSE.,
     A '('' ********   MPPOPN - ILLEGAL VARIABLE NAME'')')
      GO  TO  999
C
C         WARNING - NO CLOSING RIGHT PARENTHESIS
C
   30 CONTINUE
      CALL  IOERRM  (.FALSE.,
     A '('' ********   MPPOPN - MISSING RIGHT PARENTHESIS'')')
C
  999 CONTINUE
      RETURN
      END
      SUBROUTINE  MPSUBS  (LEOL, LSUBL)
C
C-----------------------------------------------------------------------
C
C     FAMILY
C     ------
C     MACRO PROCESSOR
C
C     PURPOSE
C     -------
C     TO EVALUATE OF THE SUBSTITUTION ESCAPE CHARACTER
C     AND DECIDE WHAT ACTION IS TO BE TAKEN
C
C     PARAMETERS
C     ----------
C     LEOL    -O-  TRUE IF AN END-OF-LINE MARKER WAS FOUND
C     LSUBL   -I-  TRUE IF NO SUBSTITUTION IS TO BE PERFORMED
C                  ON THIS LINE
C
C-----------------------------------------------------------------------
C
C         GLOBAL CONSTANTS
C
      CHARACTER*1         CA,      CBLANK,  CC,      CI,      CLEFT,
     A                    CMINUS,  CPLUS,   CPOINT,  CQUOTE,  CRIGHT,
     B                    CZ,      C0,      C9
      COMMON  / GLCOMC /  CA,      CBLANK,  CC,      CI,      CLEFT,
     A                    CMINUS,  CPLUS,   CPOINT,  CQUOTE,  CRIGHT,
     B                    CZ,      C0,      C9
C
C         INPUT / OUTPUT CONTROL INTERFACE
C
      CHARACTER*1         CBUFFR(2000)
      LOGICAL             LBREAK,  LFORT,   LISTI,   LISTO
      INTEGER             ICBADD,  ICBEND,  ICBEOL,  ICBSUB,  ICB0,
     A                    ICB1,    ICB2,    ICB3,    ICBDIM,  ICPLI,
     B                    ICPLO,   ILCTR,   ILNMBR,  ILPP,    IPAGE,
     A                    IUNITE,  IUNITI,  IUNITL,  IUNITO
      COMMON  / IOCOMC /  CBUFFR
      COMMON  / IOCOMI /  ICBADD,  ICBEND,  ICBEOL,  ICBSUB,  ICB0,
     A                    ICB1,    ICB2,    ICB3,    ICBDIM,  ICPLI,
     B                    ICPLO,   ILCTR,   ILNMBR,  ILPP,    IPAGE,
     A                    IUNITE,  IUNITI,  IUNITL,  IUNITO
      COMMON  / IOCOML /  LBREAK,  LFORT,   LISTI,   LISTO
C
C         MEMORY MANAGER INTERFACE
C
      CHARACTER*1         CSTORE(20000)
      INTEGER             IHASH(601),       ISTORE(6000)
      INTEGER             ICSDIM,  ICSP1,   ICSP2,   IHADIM,
     A                    ISFREE,  ISTDIM,  IS2HDC,  IS2HDS
      COMMON  / MMCOMC /  CSTORE
      COMMON  / MMCOMH /  IHASH
      COMMON  / MMCOMS /  ISTORE
      COMMON  / MMCOMI /  ICSDIM,  ICSP1,   ICSP2,   IHADIM,
     A                    ISFREE,  ISTDIM,  IS2HDC,  IS2HDS
C
C         MACRO PROCESSOR INTERFACE
C
      CHARACTER*1         CDIV,    CEOL,     CEOR,    CONC,    CSUB,
     A                    CTOP
      LOGICAL             LEMPTY,  LSUB
      COMMON  / MPCOMC /  CDIV,    CEOL,    CEOR,    CONC,    CSUB,
     A                    CTOP
      COMMON  / MPCOML /  LEMPTY,  LSUB
C
C         LOCAL VARIABLES AND PARAMETERS
C
      CHARACTER*1         C
      LOGICAL             LEOL, LSUBL
      EXTERNAL            MPEOL, MPMAC, MMPOPC
C
C         DETERMINE WHAT FOLLOWS THE SUBSTITUTION PREFIX CHARACTER
C
      LEOL  =  .FALSE.
      CALL  MMPOPC  (C, 2, CTOP, LEMPTY)
      IF  (CEOL .EQ. CTOP)              GO  TO  10
      IF  (CEOR .EQ. CTOP)              GO  TO  20
      IF  (CONC .EQ. CTOP)              GO  TO  30
      IF  (.NOT. (LSUB .AND. LSUBL))    GO  TO  999
      IF  (CSUB .EQ. CTOP)              GO  TO  40
      GO  TO  50
C
C         PROCESS AN END-OF-LINE MARKER
C
   10 CONTINUE
      CALL  MPEOL
      LEOL    =  .TRUE.
      GO  TO  999
C
C         PROCESS AN END-OF-RECORD MARKER
C
   20 CONTINUE
      CALL  MMPOPC  (C, 3, CTOP, LEMPTY)
      ICB2  =  ICB2 - 2
      GO  TO  999
C
C         PROCESS A CONTINUATION CHARACTER
C
   30 CONTINUE
      CALL  MMPOPC  (C, 2, CTOP, LEMPTY)
      CALL  MMPOPC  (C, 2, CTOP, LEMPTY)
      CALL  MMPOPC  (C, 2, CTOP, LEMPTY)
      ICB2    =  ICB2 - 4
      GO  TO  999
C
C         PROCESS AN EMBEDDED SUBSTITUTION PREFIX CHARACTER
C
   40 CONTINUE
      CALL  MMPOPC  (C, 2, CTOP, LEMPTY)
      ICB2  =  ICB2 - 1
      GO  TO  999
C
C         A LIST OR MACRO SUBSTITUTION HAS BEEN ENCOUNTERED
C
   50 CONTINUE
      CALL  MPMAC
C
  999 CONTINUE
      RETURN
      END
      SUBROUTINE  TPAPPE
C
C-----------------------------------------------------------------------
C
C     FAMILY
C     ------
C     TEMPLATE PROCESSOR
C
C     PURPOSE
C     -------
C     TO PROCESS APPEND DIRECTIVES
C
C-----------------------------------------------------------------------
C
C         GLOBAL CONSTANTS
C
      CHARACTER*1         CA,      CBLANK,  CC,      CI,      CLEFT,
     A                    CMINUS,  CPLUS,   CPOINT,  CQUOTE,  CRIGHT,
     B                    CZ,      C0,      C9
      COMMON  / GLCOMC /  CA,      CBLANK,  CC,      CI,      CLEFT,
     A                    CMINUS,  CPLUS,   CPOINT,  CQUOTE,  CRIGHT,
     B                    CZ,      C0,      C9
C
C         INPUT / OUTPUT CONTROL INTERFACE
C
      CHARACTER*1         CBUFFR(2000)
      LOGICAL             LBREAK,  LFORT,   LISTI,   LISTO
      INTEGER             ICBADD,  ICBEND,  ICBEOL,  ICBSUB,  ICB0,
     A                    ICB1,    ICB2,    ICB3,    ICBDIM,  ICPLI,
     B                    ICPLO,   ILCTR,   ILNMBR,  ILPP,    IPAGE,
     A                    IUNITE,  IUNITI,  IUNITL,  IUNITO
      COMMON  / IOCOMC /  CBUFFR
      COMMON  / IOCOMI /  ICBADD,  ICBEND,  ICBEOL,  ICBSUB,  ICB0,
     A                    ICB1,    ICB2,    ICB3,    ICBDIM,  ICPLI,
     B                    ICPLO,   ILCTR,   ILNMBR,  ILPP,    IPAGE,
     A                    IUNITE,  IUNITI,  IUNITL,  IUNITO
      COMMON  / IOCOML /  LBREAK,  LFORT,   LISTI,   LISTO
C
C         MEMORY MANAGER INTERFACE
C
      CHARACTER*1         CSTORE(20000)
      INTEGER             IHASH(601),       ISTORE(6000)
      INTEGER             ICSDIM,  ICSP1,   ICSP2,   IHADIM,
     A                    ISFREE,  ISTDIM,  IS2HDC,  IS2HDS
      COMMON  / MMCOMC /  CSTORE
      COMMON  / MMCOMH /  IHASH
      COMMON  / MMCOMS /  ISTORE
      COMMON  / MMCOMI /  ICSDIM,  ICSP1,   ICSP2,   IHADIM,
     A                    ISFREE,  ISTDIM,  IS2HDC,  IS2HDS
C
C         MACRO PROCESSOR INTERFACE
C
      CHARACTER*1         CDIV,    CEOL,     CEOR,    CONC,    CSUB,
     A                    CTOP
      LOGICAL             LEMPTY,  LSUB
      COMMON  / MPCOMC /  CDIV,    CEOL,    CEOR,    CONC,    CSUB,
     A                    CTOP
      COMMON  / MPCOML /  LEMPTY,  LSUB
C
C         TEMPLATE PROCESSOR INTERFACE
C
      CHARACTER*1         CDIR,    CSTAR
      INTEGER             ICBP1(4),         ICBP2(4)
      INTEGER             ITOPDO,  IARGS, INESTD, INESTF
      LOGICAL             LCOL1,   LDIRL,   LEND,    LINITM,  L1TRIP
      COMMON  / TPCOMC /  CDIR,    CSTAR
      COMMON  / TPCOMI /  ICBP1,   ITOPDO,  IARGS,   ICBP2,   INESTD,
     B                    INESTF
      COMMON  / TPCOML /  LCOL1,   LDIRL,   LEND,    LINITM,  L1TRIP
C
C         LOCAL VARIABLES AND PARAMETERS
C
      CHARACTER*1         C(15)
      INTEGER             ID(7), IK(3)
      INTEGER             IKYDIM, ICKDIM, IDSDIM
      LOGICAL             LERROR
      EXTERNAL            TPSYNT, TPRDBL, MMAPPV, IOERRM
      SAVE IKYDIM, IK, ICKDIM, C, IDSDIM, ID
      DATA
     A     IKYDIM, ICKDIM,
     B     IK(1),
     C     IK(2),
     D     IK(3)
     E  /  3,      15,
     F     6,
     G     6,
     H     3  /
      DATA
     B     C(1),  C(2),  C(3),  C(4),  C(5),  C(6),
     C     C(7),  C(8),  C(9),  C(10), C(11), C(12),
     D     C(13), C(14), C(15)
     E  /
     F     'A',   'P',   'P',   'E',   'N',   'D',
     G     'E',   'N',   'D',   'A',   'P',   'P',
     H     'E',   'N',   'D'  /
      DATA
     A     IDSDIM,
     B     ID(1), ID(2), ID(3), ID(4), ID(5), ID(6), ID(7)
     C  /  7,
     D     1,      5,      -3,     6,      6,      2,      7  /
C
C         CHECK SYNTAX
C
      CALL  TPSYNT  (ID, IDSDIM, .FALSE., LERROR)
      IF  (LERROR)        GO  TO  999
      IF  (IARGS .EQ. 2)  GO  TO  10
C
C         PROCESS A MULTI-LINE APPEND STATEMENT
C
      ICBP1(2)  =  ICBEOL + 1
C
C         READ A BLOCK; UNTIL *ENDAPP
C
      CALL  TPRDBL  (IK, IKYDIM, C, ICKDIM, ID, IDSDIM,
     A                           .TRUE., .FALSE., .TRUE., LERROR)
      IF  (LEND)    GO  TO  20
      IF  (LERROR)  GO  TO  999
      ICBP2(2)  =  ICB0 - 1
C
C         APPEND THE VALUE
C
   10 CONTINUE
      CALL  MMAPPV  (CBUFFR, ICBP1(1), ICBP2(1),
     A               CBUFFR, ICBP1(2), ICBP2(2))
      GO  TO  999
C
   20 CONTINUE
      CALL  IOERRM  (.FALSE.,
     A '('' ********   TPAPPE - APPEND HAS NO MATCHING ENDAPP'')')
C
  999 CONTINUE
      RETURN
      END
      SUBROUTINE  TPCHKD
C
C-----------------------------------------------------------------------
C
C     FAMILY
C     ------
C     TEMPLATE PROCESSOR
C
C     PURPOSE
C     -------
C     TO DETERMINE IF A LINE CONTAINS A DIRECTIVE
C
C-----------------------------------------------------------------------
C
C         GLOBAL CONSTANTS
C
      CHARACTER*1         CA,      CBLANK,  CC,      CI,      CLEFT,
     A                    CMINUS,  CPLUS,   CPOINT,  CQUOTE,  CRIGHT,
     B                    CZ,      C0,      C9
      COMMON  / GLCOMC /  CA,      CBLANK,  CC,      CI,      CLEFT,
     A                    CMINUS,  CPLUS,   CPOINT,  CQUOTE,  CRIGHT,
     B                    CZ,      C0,      C9
C
C         INPUT / OUTPUT CONTROL INTERFACE
C
      CHARACTER*1         CBUFFR(2000)
      LOGICAL             LBREAK,  LFORT,   LISTI,   LISTO
      INTEGER             ICBADD,  ICBEND,  ICBEOL,  ICBSUB,  ICB0,
     A                    ICB1,    ICB2,    ICB3,    ICBDIM,  ICPLI,
     B                    ICPLO,   ILCTR,   ILNMBR,  ILPP,    IPAGE,
     A                    IUNITE,  IUNITI,  IUNITL,  IUNITO
      COMMON  / IOCOMC /  CBUFFR
      COMMON  / IOCOMI /  ICBADD,  ICBEND,  ICBEOL,  ICBSUB,  ICB0,
     A                    ICB1,    ICB2,    ICB3,    ICBDIM,  ICPLI,
     B                    ICPLO,   ILCTR,   ILNMBR,  ILPP,    IPAGE,
     A                    IUNITE,  IUNITI,  IUNITL,  IUNITO
      COMMON  / IOCOML /  LBREAK,  LFORT,   LISTI,   LISTO
C
C         MEMORY MANAGER INTERFACE
C
      CHARACTER*1         CSTORE(20000)
      INTEGER             IHASH(601),       ISTORE(6000)
      INTEGER             ICSDIM,  ICSP1,   ICSP2,   IHADIM,
     A                    ISFREE,  ISTDIM,  IS2HDC,  IS2HDS
      COMMON  / MMCOMC /  CSTORE
      COMMON  / MMCOMH /  IHASH
      COMMON  / MMCOMS /  ISTORE
      COMMON  / MMCOMI /  ICSDIM,  ICSP1,   ICSP2,   IHADIM,
     A                    ISFREE,  ISTDIM,  IS2HDC,  IS2HDS
C
C         MACRO PROCESSOR INTERFACE
C
      CHARACTER*1         CDIV,    CEOL,     CEOR,    CONC,    CSUB,
     A                    CTOP
      LOGICAL             LEMPTY,  LSUB
      COMMON  / MPCOMC /  CDIV,    CEOL,    CEOR,    CONC,    CSUB,
     A                    CTOP
      COMMON  / MPCOML /  LEMPTY,  LSUB
C
C         TEMPLATE PROCESSOR INTERFACE
C
      CHARACTER*1         CDIR,    CSTAR
      INTEGER             ICBP1(4),         ICBP2(4)
      INTEGER             ITOPDO,  IARGS, INESTD, INESTF
      LOGICAL             LCOL1,   LDIRL,   LEND,    LINITM,  L1TRIP
      COMMON  / TPCOMC /  CDIR,    CSTAR
      COMMON  / TPCOMI /  ICBP1,   ITOPDO,  IARGS,   ICBP2,   INESTD,
     B                    INESTF
      COMMON  / TPCOML /  LCOL1,   LDIRL,   LEND,    LINITM,  L1TRIP
C
C         LOCAL VARIABLES AND PARAMETERS
C
      LOGICAL             LEOL
      EXTERNAL            UTRDBL
C
      LDIRL  =  .FALSE.
      IF  (ICB1 .GT. ICB2)          GO  TO  999
      IF  (CBUFFR(ICB1) .EQ. CDIR)  GO  TO  10
      IF  (LCOL1)                   GO  TO  999
      CALL  UTRDBL  (CBUFFR, ICB1, ICB2, LEOL)
      IF  (LEOL)                    GO  TO  999
      IF  (CBUFFR(ICB1) .NE. CDIR)  GO  TO  999
C
   10 CONTINUE
      ICB1   =  ICB1 + 1
      LDIRL  =  .TRUE.
C
  999 CONTINUE
      RETURN
      END
      SUBROUTINE  TPCOMM
C
C-----------------------------------------------------------------------
C
C     FAMILY
C     ------
C     TEMPLATE PROCESSOR
C
C     PURPOSE
C     -------
C     TO PROCESS COMMENT DIRECTIVES
C
C-----------------------------------------------------------------------
C
C         GLOBAL CONSTANTS
C
      CHARACTER*1         CA,      CBLANK,  CC,      CI,      CLEFT,
     A                    CMINUS,  CPLUS,   CPOINT,  CQUOTE,  CRIGHT,
     B                    CZ,      C0,      C9
      COMMON  / GLCOMC /  CA,      CBLANK,  CC,      CI,      CLEFT,
     A                    CMINUS,  CPLUS,   CPOINT,  CQUOTE,  CRIGHT,
     B                    CZ,      C0,      C9
C
C         INPUT / OUTPUT CONTROL INTERFACE
C
      CHARACTER*1         CBUFFR(2000)
      LOGICAL             LBREAK,  LFORT,   LISTI,   LISTO
      INTEGER             ICBADD,  ICBEND,  ICBEOL,  ICBSUB,  ICB0,
     A                    ICB1,    ICB2,    ICB3,    ICBDIM,  ICPLI,
     B                    ICPLO,   ILCTR,   ILNMBR,  ILPP,    IPAGE,
     A                    IUNITE,  IUNITI,  IUNITL,  IUNITO
      COMMON  / IOCOMC /  CBUFFR
      COMMON  / IOCOMI /  ICBADD,  ICBEND,  ICBEOL,  ICBSUB,  ICB0,
     A                    ICB1,    ICB2,    ICB3,    ICBDIM,  ICPLI,
     B                    ICPLO,   ILCTR,   ILNMBR,  ILPP,    IPAGE,
     A                    IUNITE,  IUNITI,  IUNITL,  IUNITO
      COMMON  / IOCOML /  LBREAK,  LFORT,   LISTI,   LISTO
C
C         MEMORY MANAGER INTERFACE
C
      CHARACTER*1         CSTORE(20000)
      INTEGER             IHASH(601),       ISTORE(6000)
      INTEGER             ICSDIM,  ICSP1,   ICSP2,   IHADIM,
     A                    ISFREE,  ISTDIM,  IS2HDC,  IS2HDS
      COMMON  / MMCOMC /  CSTORE
      COMMON  / MMCOMH /  IHASH
      COMMON  / MMCOMS /  ISTORE
      COMMON  / MMCOMI /  ICSDIM,  ICSP1,   ICSP2,   IHADIM,
     A                    ISFREE,  ISTDIM,  IS2HDC,  IS2HDS
C
C         MACRO PROCESSOR INTERFACE
C
      CHARACTER*1         CDIV,    CEOL,     CEOR,    CONC,    CSUB,
     A                    CTOP
      LOGICAL             LEMPTY,  LSUB
      COMMON  / MPCOMC /  CDIV,    CEOL,    CEOR,    CONC,    CSUB,
     A                    CTOP
      COMMON  / MPCOML /  LEMPTY,  LSUB
C
C         TEMPLATE PROCESSOR INTERFACE
C
      CHARACTER*1         CDIR,    CSTAR
      INTEGER             ICBP1(4),         ICBP2(4)
      INTEGER             ITOPDO,  IARGS, INESTD, INESTF
      LOGICAL             LCOL1,   LDIRL,   LEND,    LINITM,  L1TRIP
      COMMON  / TPCOMC /  CDIR,    CSTAR
      COMMON  / TPCOMI /  ICBP1,   ITOPDO,  IARGS,   ICBP2,   INESTD,
     B                    INESTF
      COMMON  / TPCOML /  LCOL1,   LDIRL,   LEND,    LINITM,  L1TRIP
C
C         LOCAL VARIABLES AND PARAMETERS
C
      CHARACTER*1         C(16)
      INTEGER             ID(1), IK(3)
      INTEGER             IKYDIM, ICKDIM, IDSDIM
      LOGICAL             LERROR
      EXTERNAL            TPSYNT, TPRDBL, IOERRM
      SAVE IKYDIM, IK, ICKDIM, C, IDSDIM, ID
      DATA                IDSDIM, ID(1)  /  1, 7  /
      DATA
     A     IKYDIM, ICKDIM,
     B     IK(1),
     C     IK(2),
     D     IK(3)
     E  /  3,      16,
     F     7,
     G     6,
     H     3  /
      DATA
     B     C(1),  C(2),  C(3),  C(4),  C(5),  C(6),  C(7),
     C     C(8),  C(9),  C(10), C(11), C(12), C(13),
     D     C(14), C(15), C(16)
     E  /
     F     'C',   'O',   'M',   'M',   'E',   'N',   'T',
     G     'E',   'N',   'D',   'C',   'O',   'M',
     H     'E',   'N',   'D'  /
C
C         CHECK SYNTAX
C
      CALL  TPSYNT  (ID, IDSDIM, .FALSE., LERROR)
      IF  (LERROR)  GO  TO  999
C
C         READ A BLOCK; UNTIL *ENDCOM
C
      CALL  TPRDBL  (IK, IKYDIM, C, ICKDIM, ID, IDSDIM,
     A                           .FALSE., .TRUE., .FALSE., LERROR)
      IF  (.NOT. LEND)  GO  TO  999
C
C         AN -END- HAS POSSIBLY BEEN ENCOUNTERED
C
      CALL  IOERRM  (.FALSE.,
     A '('' ********   TPCOMM - COMMENT HAS NO MATCHING ENDCOM'')')
C
  999 CONTINUE
      RETURN
      END
      SUBROUTINE  TPDELE
C
C-----------------------------------------------------------------------
C
C     FAMILY
C     ------
C     TEMPLATE PROCESSOR
C
C     PURPOSE
C     -------
C     TO PROCESS DELETE DIRECTIVES
C
C-----------------------------------------------------------------------
C
C         GLOBAL CONSTANTS
C
      CHARACTER*1         CA,      CBLANK,  CC,      CI,      CLEFT,
     A                    CMINUS,  CPLUS,   CPOINT,  CQUOTE,  CRIGHT,
     B                    CZ,      C0,      C9
      COMMON  / GLCOMC /  CA,      CBLANK,  CC,      CI,      CLEFT,
     A                    CMINUS,  CPLUS,   CPOINT,  CQUOTE,  CRIGHT,
     B                    CZ,      C0,      C9
C
C         INPUT / OUTPUT CONTROL INTERFACE
C
      CHARACTER*1         CBUFFR(2000)
      LOGICAL             LBREAK,  LFORT,   LISTI,   LISTO
      INTEGER             ICBADD,  ICBEND,  ICBEOL,  ICBSUB,  ICB0,
     A                    ICB1,    ICB2,    ICB3,    ICBDIM,  ICPLI,
     B                    ICPLO,   ILCTR,   ILNMBR,  ILPP,    IPAGE,
     A                    IUNITE,  IUNITI,  IUNITL,  IUNITO
      COMMON  / IOCOMC /  CBUFFR
      COMMON  / IOCOMI /  ICBADD,  ICBEND,  ICBEOL,  ICBSUB,  ICB0,
     A                    ICB1,    ICB2,    ICB3,    ICBDIM,  ICPLI,
     B                    ICPLO,   ILCTR,   ILNMBR,  ILPP,    IPAGE,
     A                    IUNITE,  IUNITI,  IUNITL,  IUNITO
      COMMON  / IOCOML /  LBREAK,  LFORT,   LISTI,   LISTO
C
C         MEMORY MANAGER INTERFACE
C
      CHARACTER*1         CSTORE(20000)
      INTEGER             IHASH(601),       ISTORE(6000)
      INTEGER             ICSDIM,  ICSP1,   ICSP2,   IHADIM,
     A                    ISFREE,  ISTDIM,  IS2HDC,  IS2HDS
      COMMON  / MMCOMC /  CSTORE
      COMMON  / MMCOMH /  IHASH
      COMMON  / MMCOMS /  ISTORE
      COMMON  / MMCOMI /  ICSDIM,  ICSP1,   ICSP2,   IHADIM,
     A                    ISFREE,  ISTDIM,  IS2HDC,  IS2HDS
C
C         MACRO PROCESSOR INTERFACE
C
      CHARACTER*1         CDIV,    CEOL,     CEOR,    CONC,    CSUB,
     A                    CTOP
      LOGICAL             LEMPTY,  LSUB
      COMMON  / MPCOMC /  CDIV,    CEOL,    CEOR,    CONC,    CSUB,
     A                    CTOP
      COMMON  / MPCOML /  LEMPTY,  LSUB
C
C         TEMPLATE PROCESSOR INTERFACE
C
      CHARACTER*1         CDIR,    CSTAR
      INTEGER             ICBP1(4),         ICBP2(4)
      INTEGER             ITOPDO,  IARGS, INESTD, INESTF
      LOGICAL             LCOL1,   LDIRL,   LEND,    LINITM,  L1TRIP
      COMMON  / TPCOMC /  CDIR,    CSTAR
      COMMON  / TPCOMI /  ICBP1,   ITOPDO,  IARGS,   ICBP2,   INESTD,
     B                    INESTF
      COMMON  / TPCOML /  LCOL1,   LDIRL,   LEND,    LINITM,  L1TRIP
C
C         LOCAL VARIABLES AND PARAMETERS
C
      INTEGER             ID(4)
      INTEGER             IDSDIM
      LOGICAL             LERROR, LFOUND
      EXTERNAL            TPSYNT, MMDELV
      SAVE IDSDIM, ID
      DATA                IDSDIM, ID(1), ID(2), ID(3), ID(4)
     A                 /  4,      1,     5,     2,     7  /
C
C         CHECK SYNTAX
C
      CALL  TPSYNT  (ID, IDSDIM, .FALSE., LERROR)
      IF  (LERROR)  GO  TO  999
C
C         DELETE THE VARIABLE
C
      CALL  MMDELV  (CBUFFR, ICBP1(1), ICBP2(1), LFOUND)
C
  999 CONTINUE
      RETURN
      END
 
      SUBROUTINE  TPDO
C
C-----------------------------------------------------------------------
C
C     FAMILY
C     ------
C     TEMPLATE PROCESSOR
C
C     PURPOSE
C     -------
C     TO PROCESS DO DIRECTIVES
C
C-----------------------------------------------------------------------
C
C         GLOBAL CONSTANTS
C
      CHARACTER*1         CA,      CBLANK,  CC,      CI,      CLEFT,
     A                    CMINUS,  CPLUS,   CPOINT,  CQUOTE,  CRIGHT,
     B                    CZ,      C0,      C9
      COMMON  / GLCOMC /  CA,      CBLANK,  CC,      CI,      CLEFT,
     A                    CMINUS,  CPLUS,   CPOINT,  CQUOTE,  CRIGHT,
     B                    CZ,      C0,      C9
C
C         INPUT / OUTPUT CONTROL INTERFACE
C
      CHARACTER*1         CBUFFR(2000)
      LOGICAL             LBREAK,  LFORT,   LISTI,   LISTO
      INTEGER             ICBADD,  ICBEND,  ICBEOL,  ICBSUB,  ICB0,
     A                    ICB1,    ICB2,    ICB3,    ICBDIM,  ICPLI,
     B                    ICPLO,   ILCTR,   ILNMBR,  ILPP,    IPAGE,
     A                    IUNITE,  IUNITI,  IUNITL,  IUNITO
      COMMON  / IOCOMC /  CBUFFR
      COMMON  / IOCOMI /  ICBADD,  ICBEND,  ICBEOL,  ICBSUB,  ICB0,
     A                    ICB1,    ICB2,    ICB3,    ICBDIM,  ICPLI,
     B                    ICPLO,   ILCTR,   ILNMBR,  ILPP,    IPAGE,
     A                    IUNITE,  IUNITI,  IUNITL,  IUNITO
      COMMON  / IOCOML /  LBREAK,  LFORT,   LISTI,   LISTO
C
C         MEMORY MANAGER INTERFACE
C
      CHARACTER*1         CSTORE(20000)
      INTEGER             IHASH(601),       ISTORE(6000)
      INTEGER             ICSDIM,  ICSP1,   ICSP2,   IHADIM,
     A                    ISFREE,  ISTDIM,  IS2HDC,  IS2HDS
      COMMON  / MMCOMC /  CSTORE
      COMMON  / MMCOMH /  IHASH
      COMMON  / MMCOMS /  ISTORE
      COMMON  / MMCOMI /  ICSDIM,  ICSP1,   ICSP2,   IHADIM,
     A                    ISFREE,  ISTDIM,  IS2HDC,  IS2HDS
C
C         MACRO PROCESSOR INTERFACE
C
      CHARACTER*1         CDIV,    CEOL,     CEOR,    CONC,    CSUB,
     A                    CTOP
      LOGICAL             LEMPTY,  LSUB
      COMMON  / MPCOMC /  CDIV,    CEOL,    CEOR,    CONC,    CSUB,
     A                    CTOP
      COMMON  / MPCOML /  LEMPTY,  LSUB
C
C         TEMPLATE PROCESSOR INTERFACE
C
      CHARACTER*1         CDIR,    CSTAR
      INTEGER             ICBP1(4),         ICBP2(4)
      INTEGER             ITOPDO,  IARGS, INESTD, INESTF
      LOGICAL             LCOL1,   LDIRL,   LEND,    LINITM,  L1TRIP
      COMMON  / TPCOMC /  CDIR,    CSTAR
      COMMON  / TPCOMI /  ICBP1,   ITOPDO,  IARGS,   ICBP2,   INESTD,
     B                    INESTF
      COMMON  / TPCOML /  LCOL1,   LDIRL,   LEND,    LINITM,  L1TRIP
C
C         LOCAL VARIABLES AND PARAMETERS
C
      CHARACTER*1         C(10), CD(1), CIN(1)
      INTEGER             ID(11), IK(3)
      LOGICAL             LERROR, LFOUND
      INTEGER             IKYDIM, ICKDIM, ICN1, ICN2, I1, I2, I3, ITEMP,
     A                    ICV1
      INTEGER             ICDDIM, ICIDIM, IDSDIM
      EXTERNAL            TPSYNT, UTBLDN, MMPUTV, UTCVCI, TPRDBL,
     A                    MMNEWI, MMPSHV, MMSETP, IOERRM
      SAVE ICDDIM, CD, ICIDIM, CIN
      SAVE IKYDIM, IK, ICKDIM, C, IDSDIM, ID
      DATA                ICDDIM / 1 /
      DATA                CD(1)  / 'D'  /
      DATA                ICIDIM / 1 /
      DATA                CIN(1)  / 'I'  /
      DATA
     A    IKYDIM,ICKDIM,IK(1),IK(2),IK(3)
     B  / 3, 10, 2, 5, 3 /
      DATA
     A    C(1), C(2),
     B    C(3), C(4), C(5), C(6), C(7),
     C    C(8), C(9), C(10)
     D /  'D', 'O',
     E    'E', 'N', 'D', 'D', 'O',
     F    'E', 'N', 'D' /
      DATA
     A     IDSDIM,
     B     ID(1),  ID(2),  ID(3),  ID(4),  ID(5),
     C     ID(6),  ID(7),  ID(8),  ID(9),  ID(10), ID(11)
     D  /  11,
     E     1,       5,       4,       6,       3,
     F     6,       -3,      10,      6,       2,       7  /
C
C         CHECK SYNTAX
C
      CALL  TPSYNT  (ID, IDSDIM, .FALSE., LERROR)
      IF  (LERROR)  GO  TO  999
      INESTD  =  INESTD + 1
      ICN1    =  ICBEND + 1
C
C         GET LOOP INDEX AND LOOP PARAMETERS
C
      CALL  UTBLDN  (CSTAR, CIN, 1, ICIDIM, INESTD,
     A               CBUFFR, ICN1, ICN2, ICBDIM, LERROR)
      ICBEND  =  ICN2
      CALL  MMPUTV  (CBUFFR, ICN1,     ICN2,
     A               CBUFFR, ICBP1(1), ICBP2(1))
      CALL  MMPUTV  (CBUFFR, ICBP1(1), ICBP2(1),
     A               CBUFFR, ICBP1(2), ICBP2(2))
      CALL  UTCVCI  (CBUFFR, ICBP1(2), ICBP2(2), I1, LERROR)
      CALL  UTCVCI  (CBUFFR, ICBP1(3), ICBP2(3), I2, LERROR)
      I3      =  1
      IF  (IARGS .EQ. 4)
     A    CALL  UTCVCI  (CBUFFR, ICBP1(4), ICBP2(4), I3, LERROR)
      IF  (L1TRIP .OR. ((I2-I1)*ISIGN(1,I3) .GE. 0))  GO  TO  10
          CALL  TPRDBL  (IK, IKYDIM, C, ICKDIM, ID, IDSDIM,
     A                               .FALSE., .TRUE., .FALSE., LERROR)
          IF  (LEND)  GO  TO  30
          INESTD  =  INESTD - 1
          GO  TO  999
C
   10 CONTINUE
      CALL  MMNEWI  (ITEMP)
      ISTORE(ITEMP)    =  I2
      ISTORE(ITEMP+1)  =  I3
      ISTORE(ITEMP+2)  =  ITOPDO
      ITOPDO           =  ITEMP
      IF  (INESTD .GT. 1)  GO  TO  20
          CALL  UTBLDN  (CSTAR, CD, 1, ICDDIM, -1,
     A                   CBUFFR, 1, ICN2, ICBDIM, LERROR)
          ICBEOL  =  ICN2
          ICV1    =  ICN2 + 1
C
C         READ A BLOCK UNTIL *ENDDO.  PUSH CONTENTS OF DO RANGE
C         ONTO STACK.
C
          CALL  TPRDBL  (IK, IKYDIM, C, ICKDIM, ID, IDSDIM,
     A                               .FALSE., .FALSE., .FALSE., LERROR)
          IF  (LEND)  GO  TO  30
          CALL  MMPUTV  (CBUFFR, 1, ICN2, CBUFFR, ICV1, ICBEOL)
          CALL  MMPSHV  (CBUFFR, 1, ICN2, 1, LEMPTY, LFOUND)
C
   20 CONTINUE
      CALL  UTBLDN  (CSTAR, CD, 1, ICDDIM, INESTD,
     A               CBUFFR, 1, ICN2, ICBDIM, LERROR)
      CALL  MMSETP  (CBUFFR, 1, ICN2)
      GO  TO  999
C
C         WARNING - MATCHING ENDDO NOT FOUND
C
   30 CONTINUE
      CALL  IOERRM  (.FALSE.,
     A '('' ********   TPDO   - DO HAS NO MATCHING ENDDO'')')
C
  999 CONTINUE
      RETURN
      END
      SUBROUTINE  TPELSE
C
C-----------------------------------------------------------------------
C
C     FAMILY
C     ------
C     TEMPLATE PROCESSOR
C
C     PURPOSE
C     -------
C     TO PROCESS ELSE DIRECTIVES
C
C-----------------------------------------------------------------------
C
C         GLOBAL CONSTANTS
C
      CHARACTER*1         CA,      CBLANK,  CC,      CI,      CLEFT,
     A                    CMINUS,  CPLUS,   CPOINT,  CQUOTE,  CRIGHT,
     B                    CZ,      C0,      C9
      COMMON  / GLCOMC /  CA,      CBLANK,  CC,      CI,      CLEFT,
     A                    CMINUS,  CPLUS,   CPOINT,  CQUOTE,  CRIGHT,
     B                    CZ,      C0,      C9
C
C         INPUT / OUTPUT CONTROL INTERFACE
C
      CHARACTER*1         CBUFFR(2000)
      LOGICAL             LBREAK,  LFORT,   LISTI,   LISTO
      INTEGER             ICBADD,  ICBEND,  ICBEOL,  ICBSUB,  ICB0,
     A                    ICB1,    ICB2,    ICB3,    ICBDIM,  ICPLI,
     B                    ICPLO,   ILCTR,   ILNMBR,  ILPP,    IPAGE,
     A                    IUNITE,  IUNITI,  IUNITL,  IUNITO
      COMMON  / IOCOMC /  CBUFFR
      COMMON  / IOCOMI /  ICBADD,  ICBEND,  ICBEOL,  ICBSUB,  ICB0,
     A                    ICB1,    ICB2,    ICB3,    ICBDIM,  ICPLI,
     B                    ICPLO,   ILCTR,   ILNMBR,  ILPP,    IPAGE,
     A                    IUNITE,  IUNITI,  IUNITL,  IUNITO
      COMMON  / IOCOML /  LBREAK,  LFORT,   LISTI,   LISTO
C
C         MEMORY MANAGER INTERFACE
C
      CHARACTER*1         CSTORE(20000)
      INTEGER             IHASH(601),       ISTORE(6000)
      INTEGER             ICSDIM,  ICSP1,   ICSP2,   IHADIM,
     A                    ISFREE,  ISTDIM,  IS2HDC,  IS2HDS
      COMMON  / MMCOMC /  CSTORE
      COMMON  / MMCOMH /  IHASH
      COMMON  / MMCOMS /  ISTORE
      COMMON  / MMCOMI /  ICSDIM,  ICSP1,   ICSP2,   IHADIM,
     A                    ISFREE,  ISTDIM,  IS2HDC,  IS2HDS
C
C         MACRO PROCESSOR INTERFACE
C
      CHARACTER*1         CDIV,    CEOL,     CEOR,    CONC,    CSUB,
     A                    CTOP
      LOGICAL             LEMPTY,  LSUB
      COMMON  / MPCOMC /  CDIV,    CEOL,    CEOR,    CONC,    CSUB,
     A                    CTOP
      COMMON  / MPCOML /  LEMPTY,  LSUB
C
C         TEMPLATE PROCESSOR INTERFACE
C
      CHARACTER*1         CDIR,    CSTAR
      INTEGER             ICBP1(4),         ICBP2(4)
      INTEGER             ITOPDO,  IARGS, INESTD, INESTF
      LOGICAL             LCOL1,   LDIRL,   LEND,    LINITM,  L1TRIP
      COMMON  / TPCOMC /  CDIR,    CSTAR
      COMMON  / TPCOMI /  ICBP1,   ITOPDO,  IARGS,   ICBP2,   INESTD,
     B                    INESTF
      COMMON  / TPCOML /  LCOL1,   LDIRL,   LEND,    LINITM,  L1TRIP
C
C         LOCAL VARIABLES AND PARAMETERS
C
      CHARACTER*1         C(10)
      INTEGER             ID(1), IDIF(4), IK(3)
      LOGICAL             LERROR
      INTEGER             IDSDIM, IDSDIF, IKYDIM, ICKDIM
      EXTERNAL            TPSYNT, TPRDBL, IOERRM
      SAVE IKYDIM, IK, ICKDIM, C, IDSDIM, ID, IDSDIF, IDIF
      DATA                IDSDIM, ID(1)  /  1, 7  /
      DATA
     A     IDSDIF,
     B     IDIF(1), IDIF(2), IDIF(3), IDIF(4)
     C  /  4,
     D     1,        6,        2,        8  /
      DATA
     A     IKYDIM, ICKDIM, IK(1), IK(2), IK(3)
     B  /  3, 10, 2, 5, 3 /
      DATA
     A     C(1),  C(2),
     B     C(3),  C(4),  C(5),  C(6),  C(7),
     C     C(8),  C(9),  C(10) /
     D     'I',   'F',
     E     'E',   'N',   'D',   'I',   'F',
     F     'E',   'N',   'D'  /
C
C         CHECK SYNTAX
C
      CALL  TPSYNT  (ID, IDSDIM, .FALSE., LERROR)
      IF  (LERROR)  GO  TO  999
      INESTF  =  INESTF - 1
      CALL  TPRDBL  (IK, IKYDIM, C, ICKDIM, IDIF, IDSDIF,
     A                           .TRUE., .TRUE., .FALSE., LERROR)
      IF  (.NOT. LEND)  GO  TO  999
C
C         AN -END- HAS BEEN ENCOUNTERED
C
      CALL  IOERRM  (.FALSE.,
     A '('' ********   TPELSE - IF HAS NO MATCHING ENDIF'')')
C
  999 CONTINUE
      RETURN
      END
      SUBROUTINE  TPENDO
C
C-----------------------------------------------------------------------
C
C     FAMILY
C     ------
C     TEMPLATE PROCESSOR
C
C     PURPOSE
C     -------
C     TO PROCESS ENDDO DIRECTIVES
C
C-----------------------------------------------------------------------
C
C         GLOBAL CONSTANTS
C
      CHARACTER*1         CA,      CBLANK,  CC,      CI,      CLEFT,
     A                    CMINUS,  CPLUS,   CPOINT,  CQUOTE,  CRIGHT,
     B                    CZ,      C0,      C9
      COMMON  / GLCOMC /  CA,      CBLANK,  CC,      CI,      CLEFT,
     A                    CMINUS,  CPLUS,   CPOINT,  CQUOTE,  CRIGHT,
     B                    CZ,      C0,      C9
C
C         INPUT / OUTPUT CONTROL INTERFACE
C
      CHARACTER*1         CBUFFR(2000)
      LOGICAL             LBREAK,  LFORT,   LISTI,   LISTO
      INTEGER             ICBADD,  ICBEND,  ICBEOL,  ICBSUB,  ICB0,
     A                    ICB1,    ICB2,    ICB3,    ICBDIM,  ICPLI,
     B                    ICPLO,   ILCTR,   ILNMBR,  ILPP,    IPAGE,
     A                    IUNITE,  IUNITI,  IUNITL,  IUNITO
      COMMON  / IOCOMC /  CBUFFR
      COMMON  / IOCOMI /  ICBADD,  ICBEND,  ICBEOL,  ICBSUB,  ICB0,
     A                    ICB1,    ICB2,    ICB3,    ICBDIM,  ICPLI,
     B                    ICPLO,   ILCTR,   ILNMBR,  ILPP,    IPAGE,
     A                    IUNITE,  IUNITI,  IUNITL,  IUNITO
      COMMON  / IOCOML /  LBREAK,  LFORT,   LISTI,   LISTO
C
C         MEMORY MANAGER INTERFACE
C
      CHARACTER*1         CSTORE(20000)
      INTEGER             IHASH(601),       ISTORE(6000)
      INTEGER             ICSDIM,  ICSP1,   ICSP2,   IHADIM,
     A                    ISFREE,  ISTDIM,  IS2HDC,  IS2HDS
      COMMON  / MMCOMC /  CSTORE
      COMMON  / MMCOMH /  IHASH
      COMMON  / MMCOMS /  ISTORE
      COMMON  / MMCOMI /  ICSDIM,  ICSP1,   ICSP2,   IHADIM,
     A                    ISFREE,  ISTDIM,  IS2HDC,  IS2HDS
C
C         MACRO PROCESSOR INTERFACE
C
      CHARACTER*1         CDIV,    CEOL,     CEOR,    CONC,    CSUB,
     A                    CTOP
      LOGICAL             LEMPTY,  LSUB
      COMMON  / MPCOMC /  CDIV,    CEOL,    CEOR,    CONC,    CSUB,
     A                    CTOP
      COMMON  / MPCOML /  LEMPTY,  LSUB
C
C         TEMPLATE PROCESSOR INTERFACE
C
      CHARACTER*1         CDIR,    CSTAR
      INTEGER             ICBP1(4),         ICBP2(4)
      INTEGER             ITOPDO,  IARGS, INESTD, INESTF
      LOGICAL             LCOL1,   LDIRL,   LEND,    LINITM,  L1TRIP
      COMMON  / TPCOMC /  CDIR,    CSTAR
      COMMON  / TPCOMI /  ICBP1,   ITOPDO,  IARGS,   ICBP2,   INESTD,
     B                    INESTF
      COMMON  / TPCOML /  LCOL1,   LDIRL,   LEND,    LINITM,  L1TRIP
C
C         LOCAL VARIABLES
C
      CHARACTER*1         CD(1), CIN(1)
      INTEGER             ID(1)
      INTEGER             ICDDIM, ICIDIM, IDSDIM, ICN2, ICV2, ICV1,
     A                    I, I2, I3, ITEMP
      LOGICAL             LERROR, LFOUND
      EXTERNAL            TPSYNT, UTBLDN, MMGETV, UTCVCI, UTCVIC,
     A                    MMPUTV, MMPOPV, MMPSHV, MMRETI, IOERRM
      SAVE  ICDDIM, CD, ICIDIM, CIN, IDSDIM, ID
      DATA                ICDDIM /  1 /
      DATA                CD(1) / 'D' /
      DATA                ICIDIM /  1 /
      DATA                CIN(1) / 'I' /
      DATA                IDSDIM, ID(1)  /  1, 7  /
C
C         CHECK SYNTAX
C
      CALL  TPSYNT  (ID, IDSDIM, .FALSE., LERROR)
      IF  (LERROR)  GO  TO  999
      IF  (INESTD .LE. 0)  GO  TO  20
      CALL  UTBLDN  (CSTAR, CIN, 1, ICIDIM, INESTD,
     A               CBUFFR, 1, ICN2, ICBDIM, LERROR)
      CALL  MMGETV  (CBUFFR, 1, ICN2, CBUFFR, 1, ICV2, ICBDIM, LFOUND)
      ICN2  =  ICV2
      ICV1  =  ICV2 + 1
      CALL  MMGETV  (CBUFFR, 1, ICN2,
     A               CBUFFR, ICV1, ICV2, ICBDIM, LFOUND)
      CALL  UTCVCI  (CBUFFR, ICV1, ICV2, I, LERROR)
      I2    =  ISTORE(ITOPDO)
      I3    =  ISTORE(ITOPDO+1)
      I     =  I + I3
      IF  ((I2-I)*ISIGN(1,I3) .LT. 0)  GO  TO  10
          CALL  UTCVIC  (CBUFFR, ICV1, ICV2, ICBDIM, I, LERROR)
          CALL  MMPUTV  (CBUFFR, 1, ICN2, CBUFFR, ICV1, ICV2)
          CALL  UTBLDN  (CSTAR, CD, 1, ICDDIM, INESTD,
     A                   CBUFFR, 1, ICN2, ICBDIM, LERROR)
          IF  (INESTD .GT. 1)  CALL  MMPOPV  (LEMPTY)
          CALL  MMPSHV  (CBUFFR, 1, ICN2, 2, LEMPTY, LFOUND)
          GO  TO  999
C
   10 CONTINUE
      INESTD  =  INESTD - 1
      ITEMP   =  ITOPDO
      ITOPDO  =  ISTORE(ITOPDO+2)
      CALL  MMRETI  (ITEMP)
      GO  TO  999
C
   20 CONTINUE
      CALL  IOERRM  (.FALSE.,
     A '('' ********   TPENDO - MISPLACED ENDDO'')')
C
C
  999 CONTINUE
      RETURN
      END
      SUBROUTINE  TPENDF
C
C-----------------------------------------------------------------------
C
C     FAMILY
C     ------
C     TEMPLATE PROCESSOR
C
C     PURPOSE
C     -------
C     TO PROCESS ENDIF DIRECTIVES
C
C-----------------------------------------------------------------------
C
C         GLOBAL CONSTANTS
C
      CHARACTER*1         CA,      CBLANK,  CC,      CI,      CLEFT,
     A                    CMINUS,  CPLUS,   CPOINT,  CQUOTE,  CRIGHT,
     B                    CZ,      C0,      C9
      COMMON  / GLCOMC /  CA,      CBLANK,  CC,      CI,      CLEFT,
     A                    CMINUS,  CPLUS,   CPOINT,  CQUOTE,  CRIGHT,
     B                    CZ,      C0,      C9
C
C         INPUT / OUTPUT CONTROL INTERFACE
C
      CHARACTER*1         CBUFFR(2000)
      LOGICAL             LBREAK,  LFORT,   LISTI,   LISTO
      INTEGER             ICBADD,  ICBEND,  ICBEOL,  ICBSUB,  ICB0,
     A                    ICB1,    ICB2,    ICB3,    ICBDIM,  ICPLI,
     B                    ICPLO,   ILCTR,   ILNMBR,  ILPP,    IPAGE,
     A                    IUNITE,  IUNITI,  IUNITL,  IUNITO
      COMMON  / IOCOMC /  CBUFFR
      COMMON  / IOCOMI /  ICBADD,  ICBEND,  ICBEOL,  ICBSUB,  ICB0,
     A                    ICB1,    ICB2,    ICB3,    ICBDIM,  ICPLI,
     B                    ICPLO,   ILCTR,   ILNMBR,  ILPP,    IPAGE,
     A                    IUNITE,  IUNITI,  IUNITL,  IUNITO
      COMMON  / IOCOML /  LBREAK,  LFORT,   LISTI,   LISTO
C
C         MEMORY MANAGER INTERFACE
C
      CHARACTER*1         CSTORE(20000)
      INTEGER             IHASH(601),       ISTORE(6000)
      INTEGER             ICSDIM,  ICSP1,   ICSP2,   IHADIM,
     A                    ISFREE,  ISTDIM,  IS2HDC,  IS2HDS
      COMMON  / MMCOMC /  CSTORE
      COMMON  / MMCOMH /  IHASH
      COMMON  / MMCOMS /  ISTORE
      COMMON  / MMCOMI /  ICSDIM,  ICSP1,   ICSP2,   IHADIM,
     A                    ISFREE,  ISTDIM,  IS2HDC,  IS2HDS
C
C         MACRO PROCESSOR INTERFACE
C
      CHARACTER*1         CDIV,    CEOL,     CEOR,    CONC,    CSUB,
     A                    CTOP
      LOGICAL             LEMPTY,  LSUB
      COMMON  / MPCOMC /  CDIV,    CEOL,    CEOR,    CONC,    CSUB,
     A                    CTOP
      COMMON  / MPCOML /  LEMPTY,  LSUB
C
C         TEMPLATE PROCESSOR INTERFACE
C
      CHARACTER*1         CDIR,    CSTAR
      INTEGER             ICBP1(4),         ICBP2(4)
      INTEGER             ITOPDO,  IARGS, INESTD, INESTF
      LOGICAL             LCOL1,   LDIRL,   LEND,    LINITM,  L1TRIP
      COMMON  / TPCOMC /  CDIR,    CSTAR
      COMMON  / TPCOMI /  ICBP1,   ITOPDO,  IARGS,   ICBP2,   INESTD,
     B                    INESTF
      COMMON  / TPCOML /  LCOL1,   LDIRL,   LEND,    LINITM,  L1TRIP
C
C         LOCAL VARIABLES
C
      INTEGER             ID(1)
      LOGICAL             LERROR
      INTEGER             IDSDIM
      EXTERNAL            TPSYNT, IOERRM
      SAVE IDSDIM, ID
      DATA                IDSDIM, ID(1)  /  1, 7  /
C
C         CHECK SYNTAX
C
      CALL  TPSYNT  (ID, IDSDIM, .FALSE., LERROR)
      IF  (LERROR)  GO  TO  999
      IF  (INESTF .LE. 0)  GO  TO  10
      INESTF  =  INESTF - 1
      GO  TO  999
C
   10 CONTINUE
      CALL  IOERRM  (.FALSE.,
     A '('' ********   TPENDF - MISPLACED ENDIF'')')
C
  999 CONTINUE
      RETURN
      END
      SUBROUTINE  TPEVAL
C
C-----------------------------------------------------------------------
C
C     FAMILY
C     ------
C     TEMPLATE PROCESSOR
C
C     PURPOSE
C     -------
C     TO CALL ROUTINES TO PROCESS DIRECTIVES
C
C-----------------------------------------------------------------------
C
C         GLOBAL CONSTANTS
C
      CHARACTER*1         CA,      CBLANK,  CC,      CI,      CLEFT,
     A                    CMINUS,  CPLUS,   CPOINT,  CQUOTE,  CRIGHT,
     B                    CZ,      C0,      C9
      COMMON  / GLCOMC /  CA,      CBLANK,  CC,      CI,      CLEFT,
     A                    CMINUS,  CPLUS,   CPOINT,  CQUOTE,  CRIGHT,
     B                    CZ,      C0,      C9
C
C         INPUT / OUTPUT CONTROL INTERFACE
C
      CHARACTER*1         CBUFFR(2000)
      LOGICAL             LBREAK,  LFORT,   LISTI,   LISTO
      INTEGER             ICBADD,  ICBEND,  ICBEOL,  ICBSUB,  ICB0,
     A                    ICB1,    ICB2,    ICB3,    ICBDIM,  ICPLI,
     B                    ICPLO,   ILCTR,   ILNMBR,  ILPP,    IPAGE,
     A                    IUNITE,  IUNITI,  IUNITL,  IUNITO
      COMMON  / IOCOMC /  CBUFFR
      COMMON  / IOCOMI /  ICBADD,  ICBEND,  ICBEOL,  ICBSUB,  ICB0,
     A                    ICB1,    ICB2,    ICB3,    ICBDIM,  ICPLI,
     B                    ICPLO,   ILCTR,   ILNMBR,  ILPP,    IPAGE,
     A                    IUNITE,  IUNITI,  IUNITL,  IUNITO
      COMMON  / IOCOML /  LBREAK,  LFORT,   LISTI,   LISTO
C
C         MEMORY MANAGER INTERFACE
C
      CHARACTER*1         CSTORE(20000)
      INTEGER             IHASH(601),       ISTORE(6000)
      INTEGER             ICSDIM,  ICSP1,   ICSP2,   IHADIM,
     A                    ISFREE,  ISTDIM,  IS2HDC,  IS2HDS
      COMMON  / MMCOMC /  CSTORE
      COMMON  / MMCOMH /  IHASH
      COMMON  / MMCOMS /  ISTORE
      COMMON  / MMCOMI /  ICSDIM,  ICSP1,   ICSP2,   IHADIM,
     A                    ISFREE,  ISTDIM,  IS2HDC,  IS2HDS
C
C         MACRO PROCESSOR INTERFACE
C
      CHARACTER*1         CDIV,    CEOL,     CEOR,    CONC,    CSUB,
     A                    CTOP
      LOGICAL             LEMPTY,  LSUB
      COMMON  / MPCOMC /  CDIV,    CEOL,    CEOR,    CONC,    CSUB,
     A                    CTOP
      COMMON  / MPCOML /  LEMPTY,  LSUB
C
C         TEMPLATE PROCESSOR INTERFACE
C
      CHARACTER*1         CDIR,    CSTAR
      INTEGER             ICBP1(4),         ICBP2(4)
      INTEGER             ITOPDO,  IARGS, INESTD, INESTF
      LOGICAL             LCOL1,   LDIRL,   LEND,    LINITM,  L1TRIP
      COMMON  / TPCOMC /  CDIR,    CSTAR
      COMMON  / TPCOMI /  ICBP1,   ITOPDO,  IARGS,   ICBP2,   INESTD,
     B                    INESTF
      COMMON  / TPCOML /  LCOL1,   LDIRL,   LEND,    LINITM,  L1TRIP
C
C         LOCAL VARIABLES
C
      CHARACTER*1         C(79)
      INTEGER             IK(16)
      INTEGER             IKYDIM, ICKDIM, I
      EXTERNAL            TPCHKD, UTRDKY, TPAPPE, TPCOMM, TPDELE, TPDO,
     A                    TPELSE, TPENDO, TPENDF, TPIF, TPINCL, TPOPT,
     B                    TPRSET, TPSET, IOERRM
      SAVE IKYDIM, IK, ICKDIM, C
      DATA
     A     IKYDIM, ICKDIM,
     B     IK(1),
     C     IK(2),
     D     IK(3),
     E     IK(4),
     F     IK(5),
     G     IK(6),
     H     IK(7),
     I     IK(8)
     J  /  16,     79,
     K     6,
     L     7,
     M     6,
     N     2,
     O     4,
     P     5,
     Q     5,
     R     2 /
      DATA
     B     C(1),  C(2),  C(3),  C(4),  C(5),  C(6),
     C     C(7),  C(8),  C(9),  C(10), C(11), C(12), C(13),
     D     C(14), C(15), C(16), C(17), C(18), C(19),
     E     C(20), C(21),
     F     C(22), C(23), C(24), C(25),
     G     C(26), C(27), C(28), C(29), C(30),
     H     C(31), C(32), C(33), C(34), C(35),
     I     C(36), C(37)
     K  /  'A',   'P',   'P',   'E',   'N',   'D',
     L     'C',   'O',   'M',   'M',   'E',   'N',   'T',
     M     'D',   'E',   'L',   'E',   'T',   'E',
     N     'D',   'O',
     O     'E',   'L',   'S',   'E',
     P     'E',   'N',   'D',   'D',   'O',
     Q     'E',   'N',   'D',   'I',   'F',
     R     'I',   'F'  /
      DATA
     A     IK(9),
     B     IK(10),
     C     IK(11),
     D     IK(12),
     E     IK(13),
     F     IK(14),
     G     IK(15),
     H     IK(16)
     I  /  7,
     J     6,
     K     5,
     L     3,
     M     6,
     N     6,
     O     6,
     P     3  /
      DATA
     A     C(38), C(39), C(40), C(41), C(42), C(43), C(44),
     B     C(45), C(46), C(47), C(48), C(49), C(50),
     C     C(51), C(52), C(53), C(54), C(55),
     D     C(56), C(57), C(58),
     E     C(59), C(60), C(61), C(62), C(63), C(64),
     F     C(65), C(66), C(67), C(68), C(69), C(70),
     G     C(71), C(72), C(73), C(74), C(75), C(76),
     H     C(77), C(78), C(79)
     I  /  'I',   'N',   'C',   'L',   'U',   'D',   'E',
     J     'O',   'P',   'T',   'I',   'O',   'N',
     K     'R',   'E',   'S',   'E',   'T',
     L     'S',   'E',   'T',
     M     'E',   'N',   'D',   'A',   'P',   'P',
     N     'E',   'N',   'D',   'C',   'O',   'M',
     O     'E',   'N',   'D',   'S',   'E',   'T',
     P     'E',   'N',   'D'  /
C
      ICB1    =  ICB0
      ICB2    =  ICB3
      CALL  TPCHKD
      IF  (.NOT. LDIRL)  GO  TO  999
C
C         A DIRECTIVE LINE HAS BEEN FOUND.  CHECK WHICH ONE IT IS
C
      CALL  UTRDKY  (CBUFFR, ICB1, ICB2, IK, IKYDIM, C, ICKDIM, I)
      GO  TO  (10,  20,  30,  40,  50,  60,  70,  80,  90,
     A              100, 110, 120, 130, 140, 150, 160, 170),  I
C
C         PROCESS -APPEND-
C
   10 CONTINUE
      CALL  TPAPPE
      GO  TO  999
C
C         PROCESS -COMMENT-
C
   20 CONTINUE
      CALL  TPCOMM
      GO  TO  999
C
C         PROCESS -DELETE-
C
   30 CONTINUE
      CALL  TPDELE
      GO  TO  999
C
C         PROCESS -DO-
C
   40 CONTINUE
      CALL  TPDO
      GO  TO  999
C
C         PROCESS -ELSE-
C
   50 CONTINUE
      CALL  TPELSE
      GO  TO  999
C
C         PROCESS -ENDDO-
C
   60 CONTINUE
      CALL  TPENDO
      GO  TO  999
C
C         PROCESS -ENDIF-
C
   70 CONTINUE
      CALL  TPENDF
      GO  TO  999
C
C         PROCESS -IF-
C
   80 CONTINUE
      CALL  TPIF
      GO  TO  999
C
C         PROCESS -INCLUDE-
C
   90 CONTINUE
      CALL  TPINCL
      GO  TO  999
C
C         PROCESS -OPTION-
C
  100 CONTINUE
      CALL  TPOPT
      GO  TO  999
C
C         PROCESS -RESET-
C
  110 CONTINUE
      CALL  TPRSET
      GO  TO  999
C
C         PROCESS -SET-
C
  120 CONTINUE
      CALL  TPSET
      GO  TO  999
C
C         PROCESS -ENDAPP-
C
  130 CONTINUE
      CALL  IOERRM  (.FALSE.,
     A '('' ********   TPEVAL - MISPLACED ENDAPP'')')
      GO  TO  999
C
C         PROCESS -ENDCOM-
C
  140 CONTINUE
      CALL  IOERRM  (.FALSE.,
     A '('' ********   TPEVAL - MISPLACED ENDCOM'')')
      GO  TO  999
C
C         PROCESS -ENDSET-
C
  150 CONTINUE
      CALL  IOERRM  (.FALSE.,
     A '('' ********   TPEVAL - MISPLACED ENDSET'')')
      GO  TO  999
C
C         PROCESS -END-
C
  160 CONTINUE
      LEND  =  .TRUE.
      GO  TO  999
C
C         PROCESS UNRECOGNIZED DIRECTIVES
C
  170 CONTINUE
      CALL  IOERRM  (.FALSE.,
     A '('' ********   TPEVAL - ILLEGAL OR MISSPELLED DIRECTIVE'')')
C
  999 CONTINUE
      RETURN
      END
      SUBROUTINE  TPEXPR  (ICV1, ICV2, LSCAN, LERROR)
C
C-----------------------------------------------------------------------
C
C     FAMILY
C     ------
C     TEMPLATE PROCESSOR
C
C     PURPOSE
C     -------
C     TO DETERMINE IF AN EXPRESSION IS VALID AND RETURN ITS VALUE.
C     CURRENTLY, EXPRESSIONS MAY CONSIST ONLY OF VARIABLES OR CONSTANTS.
C
C     PARAMETERS
C     ----------
C     ICV1    -I-  INDEX INTO CBUFFR OF THE FIRST
C                  CHARACTER IN THE EXPRESSION
C     ICV2    -I-  INDEX OF THE LAST CHARACTER
C     LSCAN   -I-  IF TRUE, THEN VALIDATE (SCAN) BUT DO NOT EVALUATE
C     LERROR  -O-  TRUE IF THE EXPRESSION WAS INVALID
C
C-----------------------------------------------------------------------
C
C         GLOBAL CONSTANTS
C
      CHARACTER*1         CA,      CBLANK,  CC,      CI,      CLEFT,
     A                    CMINUS,  CPLUS,   CPOINT,  CQUOTE,  CRIGHT,
     B                    CZ,      C0,      C9
      COMMON  / GLCOMC /  CA,      CBLANK,  CC,      CI,      CLEFT,
     A                    CMINUS,  CPLUS,   CPOINT,  CQUOTE,  CRIGHT,
     B                    CZ,      C0,      C9
C
C         INPUT / OUTPUT CONTROL INTERFACE
C
      CHARACTER*1         CBUFFR(2000)
      LOGICAL             LBREAK,  LFORT,   LISTI,   LISTO
      INTEGER             ICBADD,  ICBEND,  ICBEOL,  ICBSUB,  ICB0,
     A                    ICB1,    ICB2,    ICB3,    ICBDIM,  ICPLI,
     B                    ICPLO,   ILCTR,   ILNMBR,  ILPP,    IPAGE,
     A                    IUNITE,  IUNITI,  IUNITL,  IUNITO
      COMMON  / IOCOMC /  CBUFFR
      COMMON  / IOCOMI /  ICBADD,  ICBEND,  ICBEOL,  ICBSUB,  ICB0,
     A                    ICB1,    ICB2,    ICB3,    ICBDIM,  ICPLI,
     B                    ICPLO,   ILCTR,   ILNMBR,  ILPP,    IPAGE,
     A                    IUNITE,  IUNITI,  IUNITL,  IUNITO
      COMMON  / IOCOML /  LBREAK,  LFORT,   LISTI,   LISTO
C
C         MEMORY MANAGER INTERFACE
C
      CHARACTER*1         CSTORE(20000)
      INTEGER             IHASH(601),       ISTORE(6000)
      INTEGER             ICSDIM,  ICSP1,   ICSP2,   IHADIM,
     A                    ISFREE,  ISTDIM,  IS2HDC,  IS2HDS
      COMMON  / MMCOMC /  CSTORE
      COMMON  / MMCOMH /  IHASH
      COMMON  / MMCOMS /  ISTORE
      COMMON  / MMCOMI /  ICSDIM,  ICSP1,   ICSP2,   IHADIM,
     A                    ISFREE,  ISTDIM,  IS2HDC,  IS2HDS
C
C         MACRO PROCESSOR INTERFACE
C
      CHARACTER*1         CDIV,    CEOL,     CEOR,    CONC,    CSUB,
     A                    CTOP
      LOGICAL             LEMPTY,  LSUB
      COMMON  / MPCOMC /  CDIV,    CEOL,    CEOR,    CONC,    CSUB,
     A                    CTOP
      COMMON  / MPCOML /  LEMPTY,  LSUB
C
C         TEMPLATE PROCESSOR INTERFACE
C
      CHARACTER*1         CDIR,    CSTAR
      INTEGER             ICBP1(4),         ICBP2(4)
      INTEGER             ITOPDO,  IARGS, INESTD, INESTF
      LOGICAL             LCOL1,   LDIRL,   LEND,    LINITM,  L1TRIP
      COMMON  / TPCOMC /  CDIR,    CSTAR
      COMMON  / TPCOMI /  ICBP1,   ITOPDO,  IARGS,   ICBP2,   INESTD,
     B                    INESTF
      COMMON  / TPCOML /  LCOL1,   LDIRL,   LEND,    LINITM,  L1TRIP
C
C         LOCAL VARIABLES AND PARAMETERS
C
      INTEGER             ICV1, ICV2
      LOGICAL             L, LERROR, LSCAN
      INTEGER             ICN1, ICN2
      EXTERNAL            UTRDNA, MMGETV, UTRDNU, UTRDQS
C
      IF  (LLE(CA,CBUFFR(ICB1)) .AND.
     A     LLE(CBUFFR(ICB1),CZ))    GO  TO  10
      IF  (LLE(C0,CBUFFR(ICB1)) .AND.
     A     LLE(CBUFFR(ICB1),C9))    GO  TO  20
      IF  (CBUFFR(ICB1) .EQ. CMINUS)  GO  TO  20
      IF  (CBUFFR(ICB1) .EQ. CPLUS)   GO  TO  20
      IF  (CBUFFR(ICB1) .EQ. CQUOTE)  GO  TO  30
      IF  (CBUFFR(ICB1) .EQ. CPOINT)  GO  TO  40
      LERROR  =  .TRUE.
      GO  TO  999
C
C         PROCESS A NAME
C
   10 CONTINUE
      CALL  UTRDNA  (CBUFFR, ICB1, ICB2, ICN1, ICN2, LERROR)
      IF  (LERROR)  GO  TO  999
      IF  (LSCAN)   GO  TO  999
      ICV1    =  ICBEND + 1
      CALL  MMGETV  (CBUFFR, ICN1, ICN2, CBUFFR, ICV1, ICV2, ICBDIM, L)
      ICBEND  =  ICV2
      LERROR  =  .NOT. L
      GO  TO  999
C
C         PROCESS A NUMBER
C
   20 CONTINUE
      CALL  UTRDNU  (CBUFFR, ICB1, ICB2, ICV1, ICV2, LERROR)
      GO  TO  999
C
C         PROCESS A QUOTED STRING
C
   30 CONTINUE
      CALL  UTRDQS  (CBUFFR, ICB1, ICB2, ICV1, ICV2, LERROR)
      GO  TO  999
C
C         PROCESS A LOGICAL CONSTANT
C
   40 CONTINUE
      CALL  UTRDQS  (CBUFFR, ICB1, ICB2, ICV1, ICV2, LERROR)
      IF  (LERROR)  GO  TO  999
      ICV1  =  ICV1 - 1
      ICV2  =  ICV2 + 1
C
  999 CONTINUE
      RETURN
      END
      SUBROUTINE  TPIF
C
C-----------------------------------------------------------------------
C
C     FAMILY
C     ------
C     TEMPLATE PROCESSOR
C
C     PURPOSE
C     -------
C     TO PROCESS IF DIRECTIVES
C
C-----------------------------------------------------------------------
C
C         GLOBAL CONSTANTS
C
      CHARACTER*1         CA,      CBLANK,  CC,      CI,      CLEFT,
     A                    CMINUS,  CPLUS,   CPOINT,  CQUOTE,  CRIGHT,
     B                    CZ,      C0,      C9
      COMMON  / GLCOMC /  CA,      CBLANK,  CC,      CI,      CLEFT,
     A                    CMINUS,  CPLUS,   CPOINT,  CQUOTE,  CRIGHT,
     B                    CZ,      C0,      C9
C
C         INPUT / OUTPUT CONTROL INTERFACE
C
      CHARACTER*1         CBUFFR(2000)
      LOGICAL             LBREAK,  LFORT,   LISTI,   LISTO
      INTEGER             ICBADD,  ICBEND,  ICBEOL,  ICBSUB,  ICB0,
     A                    ICB1,    ICB2,    ICB3,    ICBDIM,  ICPLI,
     B                    ICPLO,   ILCTR,   ILNMBR,  ILPP,    IPAGE,
     A                    IUNITE,  IUNITI,  IUNITL,  IUNITO
      COMMON  / IOCOMC /  CBUFFR
      COMMON  / IOCOMI /  ICBADD,  ICBEND,  ICBEOL,  ICBSUB,  ICB0,
     A                    ICB1,    ICB2,    ICB3,    ICBDIM,  ICPLI,
     B                    ICPLO,   ILCTR,   ILNMBR,  ILPP,    IPAGE,
     A                    IUNITE,  IUNITI,  IUNITL,  IUNITO
      COMMON  / IOCOML /  LBREAK,  LFORT,   LISTI,   LISTO
C
C         MEMORY MANAGER INTERFACE
C
      CHARACTER*1         CSTORE(20000)
      INTEGER             IHASH(601),       ISTORE(6000)
      INTEGER             ICSDIM,  ICSP1,   ICSP2,   IHADIM,
     A                    ISFREE,  ISTDIM,  IS2HDC,  IS2HDS
      COMMON  / MMCOMC /  CSTORE
      COMMON  / MMCOMH /  IHASH
      COMMON  / MMCOMS /  ISTORE
      COMMON  / MMCOMI /  ICSDIM,  ICSP1,   ICSP2,   IHADIM,
     A                    ISFREE,  ISTDIM,  IS2HDC,  IS2HDS
C
C         MACRO PROCESSOR INTERFACE
C
      CHARACTER*1         CDIV,    CEOL,     CEOR,    CONC,    CSUB,
     A                    CTOP
      LOGICAL             LEMPTY,  LSUB
      COMMON  / MPCOMC /  CDIV,    CEOL,    CEOR,    CONC,    CSUB,
     A                    CTOP
      COMMON  / MPCOML /  LEMPTY,  LSUB
C
C         TEMPLATE PROCESSOR INTERFACE
C
      CHARACTER*1         CDIR,    CSTAR
      INTEGER             ICBP1(4),         ICBP2(4)
      INTEGER             ITOPDO,  IARGS, INESTD, INESTF
      LOGICAL             LCOL1,   LDIRL,   LEND,    LINITM,  L1TRIP
      COMMON  / TPCOMC /  CDIR,    CSTAR
      COMMON  / TPCOMI /  ICBP1,   ITOPDO,  IARGS,   ICBP2,   INESTD,
     B                    INESTF
      COMMON  / TPCOML /  LCOL1,   LDIRL,   LEND,    LINITM,  L1TRIP
C
C         LOCAL VARIABLES AND PARAMETERS
C
      CHARACTER*1         C(14), CN(2)
      INTEGER             ID(7), IK(4)
      INTEGER             ICNDIM, IKYDIM, ICKDIM, IDSDIM, LEN1, LEN2,
     A                    I1, I2, I, INEST
      LOGICAL             LERROR, LFOUND, LVALUE
      EXTERNAL            TPSYNT, UTCVCL, MMPUTV, MMPSHV, MPLINE,
     A                    TPCHKD, UTRDKY, IOERRM
      SAVE IKYDIM, IK, ICKDIM, C, ICNDIM, CN, IDSDIM, ID
      DATA                ICNDIM / 2 /
      DATA                CN(1), CN(2)  /  '*', 'F'  /
      DATA
     A     IKYDIM, ICKDIM,
     B     IK(1),
     C     IK(2),
     D     IK(3),
     E     IK(4)
     F  /  4,      14,
     G     2,
     H     4,
     I     5,
     J     3 /
      DATA
     B     C(1),  C(2),
     C     C(3),  C(4),  C(5),  C(6),
     D     C(7),  C(8),  C(9),  C(10), C(11),
     E     C(12), C(13), C(14)
     G  /  'I',   'F',
     H     'E',   'L',   'S',   'E',
     I     'E',   'N',   'D',   'I',   'F',
     J     'E',   'N',   'D'  /
      DATA    IDSDIM, ID(1), ID(2), ID(3), ID(4), ID(5), ID(6), ID(7)
     A          /  7,    1,     6,    -4,     6,     6,     2,     8 /
C
C         CHECK SYNTAX
C
      CALL  TPSYNT  (ID, IDSDIM, .FALSE., LERROR)
      IF  (LERROR)  GO  TO  999
      IF ((IARGS .EQ. 3) .OR.
     A    ((IARGS .EQ. 2) .AND. (ICBEOL.NE.ICBP2(2))))  GO TO 3
C
C       HAVE FOUND FORM: '*IF(L)'
C
      CALL  UTCVCL  (CBUFFR, ICBP1(1), ICBP2(1), LVALUE, LERROR)
      IF (IARGS .EQ. 1)  GO TO 10
      GO TO 9
    3 CONTINUE
C
C       HAVE FOUND FORM: '*IF(EXP1=EXP2)'
C
      LVALUE = .FALSE.
      LEN1 = ICBP2(1) - ICBP1(1) + 1
      LEN2 = ICBP2(2) - ICBP1(2) + 1
      IF (LEN1 .NE. LEN2)  GO TO 8
      I1 = ICBP1(1)
      I2 = ICBP1(2)
      DO 5 I = 1, LEN1
        IF (CBUFFR(I1) .NE. CBUFFR(I2)) GO TO 8
        I1 = I1 + 1
        I2 = I2 + 1
    5 CONTINUE
      LVALUE = .TRUE.
    8 CONTINUE
      IF (IARGS .EQ. 2) GO TO 10
C
C
C         PROCESS A ONE-LINE IF STATEMENT
C
    9 CONTINUE
      IF  (.NOT. LVALUE)       GO  TO  999
      CALL  MMPUTV  (CN, 1, ICNDIM, CBUFFR, ICBP1(IARGS), ICBP2(IARGS))
      CALL  MMPSHV  (CN, 1, ICNDIM, 1, LEMPTY, LFOUND)
      GO  TO  999
C
C         PROCESS A MULTI-LINE IF STATEMENT
C
   10 CONTINUE
      INESTF  =  INESTF + 1
      IF  (LVALUE)  GO  TO  999
      INEST   =  INESTF
C
   20 CONTINUE
          CALL  MPLINE   (.FALSE.)
          CALL  TPCHKD
          IF  (.NOT. LDIRL)  GO  TO  20
          CALL  UTRDKY  (CBUFFR, ICB1, ICB2, IK, IKYDIM, C, ICKDIM, I)
      GO  TO  (30,  40,  50,  60,  20),  I
C
C         AN -IF- HAS BEEN ENCOUNTERED
C
   30 CONTINUE
      CALL  TPSYNT  (ID, IDSDIM, .TRUE., LERROR)
      IF  (IARGS .EQ. 1)  INESTF  =  INESTF + 1
C
C             IF IARGS=2 AND ICB1 > ICB2, ASSUME
C             DIRECTIVE IS OF FORM   '*IF(ARG1 = ARG2)'
C
      IF  ((IARGS .EQ. 2) .AND. (ICB1 .GT. ICB2))  INESTF = INESTF + 1
      GO  TO  20
C
C         AN -ELSE- HAS BEEN ENCOUNTERED
C
   40 CONTINUE
      IF  (INESTF .LE. INEST)  GO  TO  999
      GO  TO  20
C
C         AN -ENDIF- HAS BEEN ENCOUNTERED
C
   50 CONTINUE
      INESTF  =  INESTF - 1
      IF  (INESTF .LT. INEST)  GO  TO  999
      GO  TO  20
C
C         AN -END- HAS POSSIBLY BEEN ENCOUNTERED
C
   60 CONTINUE
      LEND  =  ICB1 .GT. ICB2
      IF  (.NOT. LEND)  GO  TO  20
      CALL  IOERRM  (.FALSE.,
     A '('' ********   TPIF   - IF HAS NO MATCHING ENDIF'')')
C
  999 CONTINUE
      RETURN
      END
      SUBROUTINE  TPINCL
C
C-----------------------------------------------------------------------
C
C     FAMILY
C     ------
C     TEMPLATE PROCESSOR
C
C     PURPOSE
C     -------
C     TO PROCESS INCLUDE DIRECTIVES
C
C-----------------------------------------------------------------------
C
C         GLOBAL CONSTANTS
C
      CHARACTER*1         CA,      CBLANK,  CC,      CI,      CLEFT,
     A                    CMINUS,  CPLUS,   CPOINT,  CQUOTE,  CRIGHT,
     B                    CZ,      C0,      C9
      COMMON  / GLCOMC /  CA,      CBLANK,  CC,      CI,      CLEFT,
     A                    CMINUS,  CPLUS,   CPOINT,  CQUOTE,  CRIGHT,
     B                    CZ,      C0,      C9
C
C         INPUT / OUTPUT CONTROL INTERFACE
C
      CHARACTER*1         CBUFFR(2000)
      LOGICAL             LBREAK,  LFORT,   LISTI,   LISTO
      INTEGER             ICBADD,  ICBEND,  ICBEOL,  ICBSUB,  ICB0,
     A                    ICB1,    ICB2,    ICB3,    ICBDIM,  ICPLI,
     B                    ICPLO,   ILCTR,   ILNMBR,  ILPP,    IPAGE,
     A                    IUNITE,  IUNITI,  IUNITL,  IUNITO
      COMMON  / IOCOMC /  CBUFFR
      COMMON  / IOCOMI /  ICBADD,  ICBEND,  ICBEOL,  ICBSUB,  ICB0,
     A                    ICB1,    ICB2,    ICB3,    ICBDIM,  ICPLI,
     B                    ICPLO,   ILCTR,   ILNMBR,  ILPP,    IPAGE,
     A                    IUNITE,  IUNITI,  IUNITL,  IUNITO
      COMMON  / IOCOML /  LBREAK,  LFORT,   LISTI,   LISTO
C
C         MEMORY MANAGER INTERFACE
C
      CHARACTER*1         CSTORE(20000)
      INTEGER             IHASH(601),       ISTORE(6000)
      INTEGER             ICSDIM,  ICSP1,   ICSP2,   IHADIM,
     A                    ISFREE,  ISTDIM,  IS2HDC,  IS2HDS
      COMMON  / MMCOMC /  CSTORE
      COMMON  / MMCOMH /  IHASH
      COMMON  / MMCOMS /  ISTORE
      COMMON  / MMCOMI /  ICSDIM,  ICSP1,   ICSP2,   IHADIM,
     A                    ISFREE,  ISTDIM,  IS2HDC,  IS2HDS
C
C         MACRO PROCESSOR INTERFACE
C
      CHARACTER*1         CDIV,    CEOL,     CEOR,    CONC,    CSUB,
     A                    CTOP
      LOGICAL             LEMPTY,  LSUB
      COMMON  / MPCOMC /  CDIV,    CEOL,    CEOR,    CONC,    CSUB,
     A                    CTOP
      COMMON  / MPCOML /  LEMPTY,  LSUB
C
C         TEMPLATE PROCESSOR INTERFACE
C
      CHARACTER*1         CDIR,    CSTAR
      INTEGER             ICBP1(4),         ICBP2(4)
      INTEGER             ITOPDO,  IARGS, INESTD, INESTF
      LOGICAL             LCOL1,   LDIRL,   LEND,    LINITM,  L1TRIP
      COMMON  / TPCOMC /  CDIR,    CSTAR
      COMMON  / TPCOMI /  ICBP1,   ITOPDO,  IARGS,   ICBP2,   INESTD,
     B                    INESTF
      COMMON  / TPCOML /  LCOL1,   LDIRL,   LEND,    LINITM,  L1TRIP
C
C         LOCAL VARIABLES AND PARAMETERS
C
      INTEGER             ID(4), IDSDIM
      LOGICAL             LERROR, LFOUND
      EXTERNAL            TPSYNT, MMPSHV, IOERRM
      SAVE IDSDIM, ID
      DATA                IDSDIM, ID(1), ID(2), ID(3), ID(4)
     A                 /  4,      1,     5,     2,     8  /
C
C         CHECK SYNTAX, AND PUSH THE VARIABLE ON THE STACK
C
      CALL  TPSYNT  (ID, IDSDIM, .FALSE., LERROR)
      IF  (LERROR)  GO  TO  999
      CALL  MMPSHV  (CBUFFR, ICBP1(1), ICBP2(1), 1, LEMPTY, LFOUND)
      IF  (LFOUND)  GO  TO  999
      CALL  IOERRM  (.FALSE.,
     A '('' ********   TPINCL - VARIABLE NOT DEFINED'')')
C
  999 CONTINUE
      RETURN
      END
      SUBROUTINE  TPINIT  (IUE0, IUI0, IUL0, IUO0)
C
C-----------------------------------------------------------------------
C
C     FAMILY
C     ------
C     TEMPLATE PROCESSOR
C
C     PURPOSE
C     -------
C     TO INITIALIZE TEMPLATE PROCESSOR VATIABLES
C
C     PARAMETERS
C     ----------
C     IUE0    -I-  UNIT NUMBER OF THE ERROR FILE
C     IUI0    -I-  UNIT NUMBER OF THE INPUT FILE
C     IUL0    -I-  UNIT NUMBER OF THE LISTING FILE
C     IUO0    -I-  UNIT NUMBER OF THE OUTPUT FILE
C
C-----------------------------------------------------------------------
C
C         GLOBAL CONSTANTS
C
      CHARACTER*1         CA,      CBLANK,  CC,      CI,      CLEFT,
     A                    CMINUS,  CPLUS,   CPOINT,  CQUOTE,  CRIGHT,
     B                    CZ,      C0,      C9
      COMMON  / GLCOMC /  CA,      CBLANK,  CC,      CI,      CLEFT,
     A                    CMINUS,  CPLUS,   CPOINT,  CQUOTE,  CRIGHT,
     B                    CZ,      C0,      C9
C
C         INPUT / OUTPUT CONTROL INTERFACE
C
      CHARACTER*1         CBUFFR(2000)
      LOGICAL             LBREAK,  LFORT,   LISTI,   LISTO
      INTEGER             ICBADD,  ICBEND,  ICBEOL,  ICBSUB,  ICB0,
     A                    ICB1,    ICB2,    ICB3,    ICBDIM,  ICPLI,
     B                    ICPLO,   ILCTR,   ILNMBR,  ILPP,    IPAGE,
     A                    IUNITE,  IUNITI,  IUNITL,  IUNITO
      COMMON  / IOCOMC /  CBUFFR
      COMMON  / IOCOMI /  ICBADD,  ICBEND,  ICBEOL,  ICBSUB,  ICB0,
     A                    ICB1,    ICB2,    ICB3,    ICBDIM,  ICPLI,
     B                    ICPLO,   ILCTR,   ILNMBR,  ILPP,    IPAGE,
     A                    IUNITE,  IUNITI,  IUNITL,  IUNITO
      COMMON  / IOCOML /  LBREAK,  LFORT,   LISTI,   LISTO
C
C         MEMORY MANAGER INTERFACE
C
      CHARACTER*1         CSTORE(20000)
      INTEGER             IHASH(601),       ISTORE(6000)
      INTEGER             ICSDIM,  ICSP1,   ICSP2,   IHADIM,
     A                    ISFREE,  ISTDIM,  IS2HDC,  IS2HDS
      COMMON  / MMCOMC /  CSTORE
      COMMON  / MMCOMH /  IHASH
      COMMON  / MMCOMS /  ISTORE
      COMMON  / MMCOMI /  ICSDIM,  ICSP1,   ICSP2,   IHADIM,
     A                    ISFREE,  ISTDIM,  IS2HDC,  IS2HDS
C
C         MACRO PROCESSOR INTERFACE
C
      CHARACTER*1         CDIV,    CEOL,     CEOR,    CONC,    CSUB,
     A                    CTOP
      LOGICAL             LEMPTY,  LSUB
      COMMON  / MPCOMC /  CDIV,    CEOL,    CEOR,    CONC,    CSUB,
     A                    CTOP
      COMMON  / MPCOML /  LEMPTY,  LSUB
C
C         TEMPLATE PROCESSOR INTERFACE
C
      CHARACTER*1         CDIR,    CSTAR
      INTEGER             ICBP1(4),         ICBP2(4)
      INTEGER             ITOPDO,  IARGS, INESTD, INESTF
      LOGICAL             LCOL1,   LDIRL,   LEND,    LINITM,  L1TRIP
      COMMON  / TPCOMC /  CDIR,    CSTAR
      COMMON  / TPCOMI /  ICBP1,   ITOPDO,  IARGS,   ICBP2,   INESTD,
     B                    INESTF
      COMMON  / TPCOML /  LCOL1,   LDIRL,   LEND,    LINITM,  L1TRIP
      INTEGER    IUE0, IUI0, IUL0, IUO0
      EXTERNAL            MMINIT
C
      IUNITE  =  IUE0
      IUNITI  =  IUI0
      IUNITL  =  IUL0
      IUNITO  =  IUO0
C
      ILNMBR  =  0
      ILCTR   =  ILPP
      INESTD  =  0
      INESTF  =  0
      IPAGE   =  0
      ITOPDO  =  0
      LEMPTY  =  .TRUE.
      LEND    =  .FALSE.
      IF  (.NOT. LINITM)  CALL  MMINIT
      LINITM  =  .TRUE.
C
      RETURN
      END
      SUBROUTINE  TPMMIN
C
C----------------------------------------------------------------------
C
C     FAMILY
C     ------
C     TEMPLATE PROCESSOR
C
C     PURPOSE
C     -------
C     THIS ROUTINE INITIALIZES TEMPLATE PROCESSOR CONSTANTS.
C
C----------------------------------------------------------------------
C
C         GLOBAL CONSTANTS
C
      CHARACTER*1         CA,      CBLANK,  CC,      CI,      CLEFT,
     A                    CMINUS,  CPLUS,   CPOINT,  CQUOTE,  CRIGHT,
     B                    CZ,      C0,      C9
      COMMON  / GLCOMC /  CA,      CBLANK,  CC,      CI,      CLEFT,
     A                    CMINUS,  CPLUS,   CPOINT,  CQUOTE,  CRIGHT,
     B                    CZ,      C0,      C9
C
C         INPUT / OUTPUT CONTROL INTERFACE
C
      CHARACTER*1         CBUFFR(2000)
      LOGICAL             LBREAK,  LFORT,   LISTI,   LISTO
      INTEGER             ICBADD,  ICBEND,  ICBEOL,  ICBSUB,  ICB0,
     A                    ICB1,    ICB2,    ICB3,    ICBDIM,  ICPLI,
     B                    ICPLO,   ILCTR,   ILNMBR,  ILPP,    IPAGE,
     A                    IUNITE,  IUNITI,  IUNITL,  IUNITO
      COMMON  / IOCOMC /  CBUFFR
      COMMON  / IOCOMI /  ICBADD,  ICBEND,  ICBEOL,  ICBSUB,  ICB0,
     A                    ICB1,    ICB2,    ICB3,    ICBDIM,  ICPLI,
     B                    ICPLO,   ILCTR,   ILNMBR,  ILPP,    IPAGE,
     A                    IUNITE,  IUNITI,  IUNITL,  IUNITO
      COMMON  / IOCOML /  LBREAK,  LFORT,   LISTI,   LISTO
C
C         MEMORY MANAGER INTERFACE
C
      CHARACTER*1         CSTORE(20000)
      INTEGER             IHASH(601),       ISTORE(6000)
      INTEGER             ICSDIM,  ICSP1,   ICSP2,   IHADIM,
     A                    ISFREE,  ISTDIM,  IS2HDC,  IS2HDS
      COMMON  / MMCOMC /  CSTORE
      COMMON  / MMCOMH /  IHASH
      COMMON  / MMCOMS /  ISTORE
      COMMON  / MMCOMI /  ICSDIM,  ICSP1,   ICSP2,   IHADIM,
     A                    ISFREE,  ISTDIM,  IS2HDC,  IS2HDS
C
C         MACRO PROCESSOR INTERFACE
C
      CHARACTER*1         CDIV,    CEOL,     CEOR,    CONC,    CSUB,
     A                    CTOP
      LOGICAL             LEMPTY,  LSUB
      COMMON  / MPCOMC /  CDIV,    CEOL,    CEOR,    CONC,    CSUB,
     A                    CTOP
      COMMON  / MPCOML /  LEMPTY,  LSUB
C
C         TEMPLATE PROCESSOR INTERFACE
C
      CHARACTER*1         CDIR,    CSTAR
      INTEGER             ICBP1(4),         ICBP2(4)
      INTEGER             ITOPDO,  IARGS, INESTD, INESTF
      LOGICAL             LCOL1,   LDIRL,   LEND,    LINITM,  L1TRIP
      COMMON  / TPCOMC /  CDIR,    CSTAR
      COMMON  / TPCOMI /  ICBP1,   ITOPDO,  IARGS,   ICBP2,   INESTD,
     B                    INESTF
      COMMON  / TPCOML /  LCOL1,   LDIRL,   LEND,    LINITM,  L1TRIP
C
      CHARACTER*1       CA0,      CBLAN0,   CC0,      CI0,      CLEFT0,
     A                  CMINU0,   CPLUS0,   CPOIN0,   CQUOT0,   CRIGH0,
     B                  CZ0,      C00,      C90,
     C                  CDIR0,    CDIV0,    CEOL0,    CEOR0,   CONC0,
     D                  CSTAR0,   CSUB0
      DATA              CA0,      CBLAN0,   CC0,      CI0,      CLEFT0,
     A                  CMINU0,   CPLUS0,   CPOIN0,   CQUOT0,   CRIGH0,
     B                  CZ0,      C00,      C90
     C               /  'A',      ' ',      'C',      'I',      '(',
     D                  '-',      '+',      '.',     '''',      ')',
     E                  'Z',      '0',      '9'  /
      DATA      CDIR0       /  '*'     /
      DATA      CDIV0       /  '/'     /
      DATA      CEOL0       /  '-'     /
      DATA      CEOR0       /  '/'     /
      DATA      CONC0       /  '+'     /
      DATA      CSTAR0      /  '*'     /
      DATA      CSUB0       /  '$'     /
      CA       =  CA0
      CBLANK   =  CBLAN0
      CC       =  CC0
      CI       =  CI0
      CLEFT    =  CLEFT0
      CMINUS   =  CMINU0
      CPLUS    =  CPLUS0
      CPOINT   =  CPOIN0
      CQUOTE   =  CQUOT0
      CRIGHT   =  CRIGH0
      CZ       =  CZ0
      C0       =  C00
      C9       =  C90
C
      CDIR     =  CDIR0
      CDIV     =  CDIV0
      CEOL     =  CEOL0
      CEOR     =  CEOR0
      CONC     =  CONC0
      CSTAR    =  CSTAR0
      CSUB     =  CSUB0
C
      ICBADD   =  1
      ICPLI    =  72
      ICPLO    =  72
      ILPP     =  58
      LBREAK   =  .FALSE.
      LCOL1    =  .TRUE.
      LFORT    =  .FALSE.
      LINITM   =  .FALSE.
      LISTI    =  .FALSE.
      LISTO    =  .FALSE.
      LSUB     =  .TRUE.
      L1TRIP   =  .FALSE.
C
      RETURN
      END
      SUBROUTINE  TPOPT
C
C-----------------------------------------------------------------------
C
C     FAMILY
C     ------
C     TEMPLATE PROCESSOR
C
C     PURPOSE
C     -------
C     TO PROCESS OPTION DIRECTIVES
C
C-----------------------------------------------------------------------
C
C         GLOBAL CONSTANTS
C
      CHARACTER*1         CA,      CBLANK,  CC,      CI,      CLEFT,
     A                    CMINUS,  CPLUS,   CPOINT,  CQUOTE,  CRIGHT,
     B                    CZ,      C0,      C9
      COMMON  / GLCOMC /  CA,      CBLANK,  CC,      CI,      CLEFT,
     A                    CMINUS,  CPLUS,   CPOINT,  CQUOTE,  CRIGHT,
     B                    CZ,      C0,      C9
C
C         INPUT / OUTPUT CONTROL INTERFACE
C
      CHARACTER*1         CBUFFR(2000)
      LOGICAL             LBREAK,  LFORT,   LISTI,   LISTO
      INTEGER             ICBADD,  ICBEND,  ICBEOL,  ICBSUB,  ICB0,
     A                    ICB1,    ICB2,    ICB3,    ICBDIM,  ICPLI,
     B                    ICPLO,   ILCTR,   ILNMBR,  ILPP,    IPAGE,
     A                    IUNITE,  IUNITI,  IUNITL,  IUNITO
      COMMON  / IOCOMC /  CBUFFR
      COMMON  / IOCOMI /  ICBADD,  ICBEND,  ICBEOL,  ICBSUB,  ICB0,
     A                    ICB1,    ICB2,    ICB3,    ICBDIM,  ICPLI,
     B                    ICPLO,   ILCTR,   ILNMBR,  ILPP,    IPAGE,
     A                    IUNITE,  IUNITI,  IUNITL,  IUNITO
      COMMON  / IOCOML /  LBREAK,  LFORT,   LISTI,   LISTO
C
C         MEMORY MANAGER INTERFACE
C
      CHARACTER*1         CSTORE(20000)
      INTEGER             IHASH(601),       ISTORE(6000)
      INTEGER             ICSDIM,  ICSP1,   ICSP2,   IHADIM,
     A                    ISFREE,  ISTDIM,  IS2HDC,  IS2HDS
      COMMON  / MMCOMC /  CSTORE
      COMMON  / MMCOMH /  IHASH
      COMMON  / MMCOMS /  ISTORE
      COMMON  / MMCOMI /  ICSDIM,  ICSP1,   ICSP2,   IHADIM,
     A                    ISFREE,  ISTDIM,  IS2HDC,  IS2HDS
C
C         MACRO PROCESSOR INTERFACE
C
      CHARACTER*1         CDIV,    CEOL,     CEOR,    CONC,    CSUB,
     A                    CTOP
      LOGICAL             LEMPTY,  LSUB
      COMMON  / MPCOMC /  CDIV,    CEOL,    CEOR,    CONC,    CSUB,
     A                    CTOP
      COMMON  / MPCOML /  LEMPTY,  LSUB
C
C         TEMPLATE PROCESSOR INTERFACE
C
      CHARACTER*1         CDIR,    CSTAR
      INTEGER             ICBP1(4),         ICBP2(4)
      INTEGER             ITOPDO,  IARGS, INESTD, INESTF
      LOGICAL             LCOL1,   LDIRL,   LEND,    LINITM,  L1TRIP
      COMMON  / TPCOMC /  CDIR,    CSTAR
      COMMON  / TPCOMI /  ICBP1,   ITOPDO,  IARGS,   ICBP2,   INESTD,
     B                    INESTF
      COMMON  / TPCOML /  LCOL1,   LDIRL,   LEND,    LINITM,  L1TRIP
C
C         LOCAL VARIABLES AND PARAMETERS
C
      CHARACTER*1         C(84), CVALUE
      INTEGER             ID(6), IK(17)
      INTEGER             IKYDIM, ICKDIM, IDSDIM, ICB, I, IVALUE
      LOGICAL             LERROR, LVALUE
      EXTERNAL            TPSYNT, UTRDKY, UTCVCI, UTCVCL, IOERRM
      SAVE IKYDIM, IK, ICKDIM, C, IDSDIM, ID
      DATA
     A     IKYDIM, ICKDIM,
     B     IK(1),
     C     IK(2),
     D     IK(3),
     E     IK(4),
     F     IK(5)
     G  /  17,     84,
     H     4,
     I     4,
     J     4,
     K     4,
     L     4  /
      DATA
     B     C(1),  C(2),  C(3),  C(4),
     C     C(5),  C(6),  C(7),  C(8),
     D     C(9),  C(10), C(11), C(12),
     E     C(13), C(14), C(15), C(16),
     F     C(17), C(18), C(19), C(20)
     G  /
     H     'C',   'D',   'I',   'R',
     I     'C',   'E',   'O',   'L',
     J     'C',   'E',   'O',   'R',
     K     'C',   'O',   'N',   'C',
     L     'C',   'S',   'U',   'B'  /
      DATA
     A     IK(6),
     B     IK(7),
     C     IK(8),
     D     IK(9),
     E     IK(10)
     F  /  5,
     G     5,
     H     6,
     I     6,
     J     6  /
      DATA
     A     C(21), C(22), C(23), C(24), C(25),
     B     C(26), C(27), C(28), C(29), C(30),
     C     C(31), C(32), C(33), C(34), C(35), C(36),
     D     C(37), C(38), C(39), C(40), C(41), C(42),
     E     C(43), C(44), C(45), C(46), C(47), C(48)
     F  /  'I',   'C',   'P',   'L',   'I',
     G     'I',   'C',   'P',   'L',   'O',
     H     'I',   'U',   'N',   'I',   'T',   'I',
     I     'I',   'U',   'N',   'I',   'T',   'L',
     J     'I',   'U',   'N',   'I',   'T',   'O'  /
      DATA
     A     IK(11),
     B     IK(12),
     C     IK(13),
     D     IK(14),
     E     IK(15),
     F     IK(16),
     G     IK(17)
     H  /  6,
     I     5,
     K     5,
     L     5,
     M     5,
     N     4,
     O     6  /
      DATA
     A     C(49), C(50), C(51), C(52), C(53), C(54),
     B     C(55), C(56), C(57), C(58), C(59),
     C     C(60), C(61), C(62), C(63), C(64),
     D     C(65), C(66), C(67), C(68), C(69),
     E     C(70), C(71), C(72), C(73), C(74),
     F     C(75), C(76), C(77), C(78),
     G     C(79), C(80), C(81), C(82), C(83), C(84)
     H  /  'L',   'B',   'R',   'E',   'A',   'K',
     I     'L',   'C',   'O',   'L',   '1',
     K     'L',   'F',   'O',   'R',   'T',
     L     'L',   'I',   'S',   'T',   'I',
     M     'L',   'I',   'S',   'T',   'O',
     N     'L',   'S',   'U',   'B',
     O     'L',   '1',   'T',   'R',   'I',   'P'  /
      DATA
     A     IDSDIM, ID(1), ID(2), ID(3), ID(4), ID(5), ID(6)
     B  /  6,      1,     5,     4,     6,     2,     7  /
C
C         CHECK SYNTAX
C
      CALL  TPSYNT  (ID, IDSDIM, .FALSE., LERROR)
      IF  (LERROR)  GO  TO  999
      ICB  =  ICBP1(1)
      CALL  UTRDKY  (CBUFFR, ICBP1(1), ICBP2(1), IK, IKYDIM,
     A                                           C,  ICKDIM, I)
      IF  (I .GT. IKYDIM)        GO  TO  220
      IF  (CBUFFR(ICB) .EQ. CC)  GO  TO  10
      IF  (CBUFFR(ICB) .EQ. CI)  GO  TO  20
      GO  TO  30
C
   10 CONTINUE
      IF  (ICBP1(2) .NE. ICBP2(2))  GO  TO  230
      ICB     =  ICBP1(2)
      CVALUE  =  CBUFFR(ICB)
      GO  TO  40
C
   20 CONTINUE
      CALL  UTCVCI  (CBUFFR, ICBP1(2), ICBP2(2), IVALUE, LERROR)
      IF  (LERROR)  GO  TO  240
      GO  TO  40
C
   30 CONTINUE
      CALL  UTCVCL  (CBUFFR, ICBP1(2), ICBP2(2), LVALUE, LERROR)
      IF  (LERROR)  GO  TO  250
C
   40 CONTINUE
      GO  TO  (50,  60,  70,  80,  90,  100, 110, 120,
     A         130, 140, 150, 160, 170, 180, 190, 200, 210),  I
C
C         PROCESS -CDIR-
C
   50 CONTINUE
      CDIR    =  CVALUE
      GO  TO  999
C
C         PROCESS -CEOL-
C
   60 CONTINUE
      CEOL    =  CVALUE
      GO  TO  999
C
C         PROCESS -CEOR-
C
   70 CONTINUE
      CEOR    =  CVALUE
      GO  TO  999
C
C         PROCESS -CONC-
C
   80 CONTINUE
      CONC    =  CVALUE
      GO  TO  999
C
C         PROCESS -CSUB-
C
   90 CONTINUE
      CSUB    =  CVALUE
      GO  TO  999
C
C         PROCESS -ICPLI-
C
  100 CONTINUE
      ICPLI   =  IVALUE
      GO  TO  999
C
C         PROCESS -ICPLO-
C
  110 CONTINUE
      ICPLO   =  IVALUE
      GO  TO  999
C
C         PROCESS -IUNITI-
C
  120 CONTINUE
      IUNITI  =  IVALUE
      GO  TO  999
C
C         PROCESS -IUNITL-
C
  130 CONTINUE
      IUNITL  =  IVALUE
      GO  TO  999
C
C         PROCESS -IUNITO-
C
  140 CONTINUE
      IUNITO  =  IVALUE
      GO  TO  999
C
C         PROCESS -LBREAK-
C
  150 CONTINUE
      LBREAK  =  LVALUE
      ICBADD  =  1
      IF  (LFORT)               ICBADD  =  -5
      IF  (LFORT .AND. LBREAK)  ICBADD  =  -9
      GO  TO  999
C
C         PROCESS -LCOL1-
C
  160 CONTINUE
      LCOL1   =  LVALUE
      GO  TO  999
C
C         PROCESS -LFORT-
C
  170 CONTINUE
      LFORT   =  LVALUE
      ICBADD  =  1
      IF  (LFORT)               ICBADD  =  -5
      IF  (LFORT .AND. LBREAK)  ICBADD  =  -9
      GO  TO  999
C
C         PROCESS -LISTI-
C
  180 CONTINUE
      LISTI   =  LVALUE
      GO  TO  999
C
C         PROCESS -LISTO-
C
  190 CONTINUE
      LISTO   =  LVALUE
      GO  TO  999
C
C         PROCESS -LSUB-
C
  200 CONTINUE
      LSUB    =  LVALUE
      GO  TO  999
C
C         PROCESS -L1TRIP-
C
  210 CONTINUE
      L1TRIP  =  LVALUE
      GO  TO  999
C
C         ERROR   - UNKNOWN OPTION NAME
C
  220 CONTINUE
      CALL  IOERRM  (.FALSE.,
     A '('' ********   TPOPT  - ILLEGAL OR MISSPELLED OPTION'')')
      GO  TO  999
C
C         ERROR   - SINGLE CHARACTER EXPECTED
C
  230 CONTINUE
      CALL  IOERRM  (.FALSE.,
     A '('' ********   TPOPT  - OPTION REQUIRES SINGLE CHARACTER'')')
      GO  TO  999
C
C         ERROR   - INTEGER EXPECTED
C
  240 CONTINUE
      CALL  IOERRM  (.FALSE.,
     A '('' ********   TPOPT  - OPTION REQUIRES AN INTEGER'')')
      GO  TO  999
C
C         ERROR   - LOGICAL VALUE EXPECTED
C
  250 CONTINUE
      CALL  IOERRM  (.FALSE.,
     A '('' ********   TPOPT  - OPTION REQUIRES A LOGICAL VALUE'')')
C
  999 CONTINUE
      RETURN
      END
      SUBROUTINE  TPRDBL  (IK, IKYDIM, C, ICKDIM, ID, IDSDIM,
     A                                 LSCAN, LSKIP, LSUBL, LERROR)
C
C-----------------------------------------------------------------------
C
C     FAMILY
C     ------
C     TEMPLATE PROCESSOR
C
C     PURPOSE
C     -------
C     TO READ A BLOCK OF STATEMENTS DELIMITED BY
C     DIRECTIVES OF THE FORM -NAME- AND -ENDNAME-.
C     THESE DIRECTIVES MAY BE NESTED.
C
C     PARAMETERS
C     ----------
C     IK      -I-  INDEXES OF DIRECTIVES IN ARRAY C
C     IKYDIM  -I-  DIMENSION OF IK (SHOULD BE 3)
C     C       -I-  CONTAINS DIRECTIVE NAMES. DIRECTIVE 1 IS -NAME-,
C                  2 IS -ENDNAME, AND 3 IS -END-.
C     ICKDIM  -I-  DIMENSION OF C (TOTAL NUMBER OF CHARACTERS)
C     ID      -I-  CONTAINS THE SYNTAX PATTERN FOR DIRECTIVE -NAME-
C     IDSDIM  -I-  DIMENSION OF ID
C     LSCAN   -I-  IF TRUE, EXPRESSIONS WILL BE SCANNED FOR ERRORS
C                  BUT NOT EVALUATED
C     LSKIP   -I-  IF TRUE, INPUT LINES ARE SKIPPED, NOT SAVED
C     LSUBL   -I-  IF TRUE, MACRO SUBSTITUTIONS WILL BE PERFORMED
C                  WHEN ENCOUNTERED WITHIN THE BLOCK
C     LERROR  -O-  TRUE IF AN ERROR WAS ENCOUNTERED
C
C-----------------------------------------------------------------------
C
C         GLOBAL CONSTANTS
C
      CHARACTER*1         CA,      CBLANK,  CC,      CI,      CLEFT,
     A                    CMINUS,  CPLUS,   CPOINT,  CQUOTE,  CRIGHT,
     B                    CZ,      C0,      C9
      COMMON  / GLCOMC /  CA,      CBLANK,  CC,      CI,      CLEFT,
     A                    CMINUS,  CPLUS,   CPOINT,  CQUOTE,  CRIGHT,
     B                    CZ,      C0,      C9
C
C         INPUT / OUTPUT CONTROL INTERFACE
C
      CHARACTER*1         CBUFFR(2000)
      LOGICAL             LBREAK,  LFORT,   LISTI,   LISTO
      INTEGER             ICBADD,  ICBEND,  ICBEOL,  ICBSUB,  ICB0,
     A                    ICB1,    ICB2,    ICB3,    ICBDIM,  ICPLI,
     B                    ICPLO,   ILCTR,   ILNMBR,  ILPP,    IPAGE,
     A                    IUNITE,  IUNITI,  IUNITL,  IUNITO
      COMMON  / IOCOMC /  CBUFFR
      COMMON  / IOCOMI /  ICBADD,  ICBEND,  ICBEOL,  ICBSUB,  ICB0,
     A                    ICB1,    ICB2,    ICB3,    ICBDIM,  ICPLI,
     B                    ICPLO,   ILCTR,   ILNMBR,  ILPP,    IPAGE,
     A                    IUNITE,  IUNITI,  IUNITL,  IUNITO
      COMMON  / IOCOML /  LBREAK,  LFORT,   LISTI,   LISTO
C
C         MEMORY MANAGER INTERFACE
C
      CHARACTER*1         CSTORE(20000)
      INTEGER             IHASH(601),       ISTORE(6000)
      INTEGER             ICSDIM,  ICSP1,   ICSP2,   IHADIM,
     A                    ISFREE,  ISTDIM,  IS2HDC,  IS2HDS
      COMMON  / MMCOMC /  CSTORE
      COMMON  / MMCOMH /  IHASH
      COMMON  / MMCOMS /  ISTORE
      COMMON  / MMCOMI /  ICSDIM,  ICSP1,   ICSP2,   IHADIM,
     A                    ISFREE,  ISTDIM,  IS2HDC,  IS2HDS
C
C         MACRO PROCESSOR INTERFACE
C
      CHARACTER*1         CDIV,    CEOL,     CEOR,    CONC,    CSUB,
     A                    CTOP
      LOGICAL             LEMPTY,  LSUB
      COMMON  / MPCOMC /  CDIV,    CEOL,    CEOR,    CONC,    CSUB,
     A                    CTOP
      COMMON  / MPCOML /  LEMPTY,  LSUB
C
C         TEMPLATE PROCESSOR INTERFACE
C
      CHARACTER*1         CDIR,    CSTAR
      INTEGER             ICBP1(4),         ICBP2(4)
      INTEGER             ITOPDO,  IARGS, INESTD, INESTF
      LOGICAL             LCOL1,   LDIRL,   LEND,    LINITM,  L1TRIP
      COMMON  / TPCOMC /  CDIR,    CSTAR
      COMMON  / TPCOMI /  ICBP1,   ITOPDO,  IARGS,   ICBP2,   INESTD,
     B                    INESTF
      COMMON  / TPCOML /  LCOL1,   LDIRL,   LEND,    LINITM,  L1TRIP
C
C         LOCAL VARIABLES AND PARAMETERS
C
      INTEGER             IKYDIM, ICKDIM, IDSDIM
      CHARACTER*(*)      C(ICKDIM)
      INTEGER             ID(IDSDIM), IK(IKYDIM)
      LOGICAL             LERROR, LSCAN, LSKIP, LSUBL
      INTEGER             I, INEST
      EXTERNAL            MPLINE, TPCHKD, UTRDKY, TPSYNT
C
      INEST  =  1
C
   10 CONTINUE
          IF  (LSKIP)  ICBEOL  =  0
          CALL  MPLINE   (LSUBL)
          CALL  TPCHKD
          IF  (.NOT. LDIRL)  GO  TO  10
          CALL  UTRDKY  (CBUFFR, ICB1, ICB2, IK, IKYDIM, C, ICKDIM, I)
      GO  TO  (20,  30,  40,  10),  I
C
C         A -NAME- DIRECTIVE HAS BEEN ENCOUNTERED
C
   20 CONTINUE
      IF  (LSCAN)  CALL  TPSYNT  (ID, IDSDIM, LSCAN, LERROR)
      IF  (LSCAN .AND. (IARGS .GE. 2))  GO  TO  10
      INEST  =  INEST + 1
      GO  TO  10
C
C         AN -ENDNAME- DIRECTIVE HAS BEEN ENCOUNTERED
C
   30 CONTINUE
      INEST  =  INEST - 1
      IF  (INEST .GT. 0)  GO  TO  10
      GO  TO  999
C
C         AN -END- DIRECTIVE HAS POSSIBLY BEEN ENCOUNTERED
C
   40 CONTINUE
      LEND  =  ICB1 .GT. ICB2
      IF  (.NOT. LEND)  GO  TO  10
C
  999 CONTINUE
      RETURN
      END
      SUBROUTINE  TPRSET
C
C-----------------------------------------------------------------------
C
C     FAMILY
C     ------
C     TEMPLATE PROCESSOR
C
C     PURPOSE
C     -------
C     TO PROCESS RESET DIRECTIVES
C
C-----------------------------------------------------------------------
C
C         GLOBAL CONSTANTS
C
      CHARACTER*1         CA,      CBLANK,  CC,      CI,      CLEFT,
     A                    CMINUS,  CPLUS,   CPOINT,  CQUOTE,  CRIGHT,
     B                    CZ,      C0,      C9
      COMMON  / GLCOMC /  CA,      CBLANK,  CC,      CI,      CLEFT,
     A                    CMINUS,  CPLUS,   CPOINT,  CQUOTE,  CRIGHT,
     B                    CZ,      C0,      C9
C
C         INPUT / OUTPUT CONTROL INTERFACE
C
      CHARACTER*1         CBUFFR(2000)
      LOGICAL             LBREAK,  LFORT,   LISTI,   LISTO
      INTEGER             ICBADD,  ICBEND,  ICBEOL,  ICBSUB,  ICB0,
     A                    ICB1,    ICB2,    ICB3,    ICBDIM,  ICPLI,
     B                    ICPLO,   ILCTR,   ILNMBR,  ILPP,    IPAGE,
     A                    IUNITE,  IUNITI,  IUNITL,  IUNITO
      COMMON  / IOCOMC /  CBUFFR
      COMMON  / IOCOMI /  ICBADD,  ICBEND,  ICBEOL,  ICBSUB,  ICB0,
     A                    ICB1,    ICB2,    ICB3,    ICBDIM,  ICPLI,
     B                    ICPLO,   ILCTR,   ILNMBR,  ILPP,    IPAGE,
     A                    IUNITE,  IUNITI,  IUNITL,  IUNITO
      COMMON  / IOCOML /  LBREAK,  LFORT,   LISTI,   LISTO
C
C         MEMORY MANAGER INTERFACE
C
      CHARACTER*1         CSTORE(20000)
      INTEGER             IHASH(601),       ISTORE(6000)
      INTEGER             ICSDIM,  ICSP1,   ICSP2,   IHADIM,
     A                    ISFREE,  ISTDIM,  IS2HDC,  IS2HDS
      COMMON  / MMCOMC /  CSTORE
      COMMON  / MMCOMH /  IHASH
      COMMON  / MMCOMS /  ISTORE
      COMMON  / MMCOMI /  ICSDIM,  ICSP1,   ICSP2,   IHADIM,
     A                    ISFREE,  ISTDIM,  IS2HDC,  IS2HDS
C
C         MACRO PROCESSOR INTERFACE
C
      CHARACTER*1         CDIV,    CEOL,     CEOR,    CONC,    CSUB,
     A                    CTOP
      LOGICAL             LEMPTY,  LSUB
      COMMON  / MPCOMC /  CDIV,    CEOL,    CEOR,    CONC,    CSUB,
     A                    CTOP
      COMMON  / MPCOML /  LEMPTY,  LSUB
C
C         TEMPLATE PROCESSOR INTERFACE
C
      CHARACTER*1         CDIR,    CSTAR
      INTEGER             ICBP1(4),         ICBP2(4)
      INTEGER             ITOPDO,  IARGS, INESTD, INESTF
      LOGICAL             LCOL1,   LDIRL,   LEND,    LINITM,  L1TRIP
      COMMON  / TPCOMC /  CDIR,    CSTAR
      COMMON  / TPCOMI /  ICBP1,   ITOPDO,  IARGS,   ICBP2,   INESTD,
     B                    INESTF
      COMMON  / TPCOML /  LCOL1,   LDIRL,   LEND,    LINITM,  L1TRIP
C
C         LOCAL VARIABLES AND PARAMETERS
C
      INTEGER             IDSDIM, ICP1,ICP2
      INTEGER             ID(4)
      LOGICAL             LERROR, LFOUND
      EXTERNAL            TPSYNT, UTBLDN, MMPUTP
      SAVE IDSDIM, ID
      DATA                IDSDIM, ID(1), ID(2), ID(3), ID(4)
     A                 /  4,      1,     5,     2,     7  /
C
C         CHECK SYNTAX
C
      CALL  TPSYNT  (ID, IDSDIM, .FALSE., LERROR)
      IF  (LERROR)  GO  TO  999
      ICP1    =  ICBEND + 1
      CALL  UTBLDN  (CDIV, CBUFFR, ICBP1(1), ICBP2(1), 1,
     A               CBUFFR, ICP1, ICP2, ICBDIM, LERROR)
      ICBEND  =  ICP2
      CALL  MMPUTP  (CBUFFR, ICBP1(1), ICBP2(1),
     A               CBUFFR, ICP1, ICP2, LFOUND)
C
  999 CONTINUE
      RETURN
      END
      SUBROUTINE  TPSET
C
C-----------------------------------------------------------------------
C
C     FAMILY
C     ------
C     TEMPLATE PROCESSOR
C
C     PURPOSE
C     -------
C     TO PROCESS SET DIRECTIVES
C
C-----------------------------------------------------------------------
C
C         GLOBAL CONSTANTS
C
      CHARACTER*1         CA,      CBLANK,  CC,      CI,      CLEFT,
     A                    CMINUS,  CPLUS,   CPOINT,  CQUOTE,  CRIGHT,
     B                    CZ,      C0,      C9
      COMMON  / GLCOMC /  CA,      CBLANK,  CC,      CI,      CLEFT,
     A                    CMINUS,  CPLUS,   CPOINT,  CQUOTE,  CRIGHT,
     B                    CZ,      C0,      C9
C
C         INPUT / OUTPUT CONTROL INTERFACE
C
      CHARACTER*1         CBUFFR(2000)
      LOGICAL             LBREAK,  LFORT,   LISTI,   LISTO
      INTEGER             ICBADD,  ICBEND,  ICBEOL,  ICBSUB,  ICB0,
     A                    ICB1,    ICB2,    ICB3,    ICBDIM,  ICPLI,
     B                    ICPLO,   ILCTR,   ILNMBR,  ILPP,    IPAGE,
     A                    IUNITE,  IUNITI,  IUNITL,  IUNITO
      COMMON  / IOCOMC /  CBUFFR
      COMMON  / IOCOMI /  ICBADD,  ICBEND,  ICBEOL,  ICBSUB,  ICB0,
     A                    ICB1,    ICB2,    ICB3,    ICBDIM,  ICPLI,
     B                    ICPLO,   ILCTR,   ILNMBR,  ILPP,    IPAGE,
     A                    IUNITE,  IUNITI,  IUNITL,  IUNITO
      COMMON  / IOCOML /  LBREAK,  LFORT,   LISTI,   LISTO
C
C         MEMORY MANAGER INTERFACE
C
      CHARACTER*1         CSTORE(20000)
      INTEGER             IHASH(601),       ISTORE(6000)
      INTEGER             ICSDIM,  ICSP1,   ICSP2,   IHADIM,
     A                    ISFREE,  ISTDIM,  IS2HDC,  IS2HDS
      COMMON  / MMCOMC /  CSTORE
      COMMON  / MMCOMH /  IHASH
      COMMON  / MMCOMS /  ISTORE
      COMMON  / MMCOMI /  ICSDIM,  ICSP1,   ICSP2,   IHADIM,
     A                    ISFREE,  ISTDIM,  IS2HDC,  IS2HDS
C
C         MACRO PROCESSOR INTERFACE
C
      CHARACTER*1         CDIV,    CEOL,     CEOR,    CONC,    CSUB,
     A                    CTOP
      LOGICAL             LEMPTY,  LSUB
      COMMON  / MPCOMC /  CDIV,    CEOL,    CEOR,    CONC,    CSUB,
     A                    CTOP
      COMMON  / MPCOML /  LEMPTY,  LSUB
C
C         TEMPLATE PROCESSOR INTERFACE
C
      CHARACTER*1         CDIR,    CSTAR
      INTEGER             ICBP1(4),         ICBP2(4)
      INTEGER             ITOPDO,  IARGS, INESTD, INESTF
      LOGICAL             LCOL1,   LDIRL,   LEND,    LINITM,  L1TRIP
      COMMON  / TPCOMC /  CDIR,    CSTAR
      COMMON  / TPCOMI /  ICBP1,   ITOPDO,  IARGS,   ICBP2,   INESTD,
     B                    INESTF
      COMMON  / TPCOML /  LCOL1,   LDIRL,   LEND,    LINITM,  L1TRIP
C
C         LOCAL VARIABLES AND PARAMETERS
C
      INTEGER             IKYDIM, ICKDIM, IDSDIM
      CHARACTER*1         C(12)
      INTEGER             ID(8), IK(3)
      LOGICAL             LERROR
      EXTERNAL            TPSYNT, TPSETM, TPRDBL, MMPUTV, IOERRM
      SAVE IKYDIM, IK, ICKDIM, C, IDSDIM, ID
      DATA
     A     IKYDIM, ICKDIM,
     B     IK(1),
     C     IK(2),
     D     IK(3)
     E  /  3,      12,
     F     3,
     G     6,
     H     3  /
      DATA
     B     C(1),  C(2),  C(3),
     C     C(4),  C(5),  C(6),  C(7),  C(8),  C(9),
     D     C(10), C(11), C(12)
     E  /
     F     'S',   'E',   'T',
     G     'E',   'N',   'D',   'S',   'E',   'T',
     H     'E',   'N',   'D'  /
      DATA
     A     IDSDIM,
     B     ID(1), ID(2), ID(3), ID(4), ID(5), ID(6), ID(7), ID(8)
     C  /  8,
     D     -1,    8,     5,     -4,    7,     6,     2,     7  /
C
C         CHECK SYNTAX
C
      CALL  TPSYNT  (ID, IDSDIM, .FALSE., LERROR)
      IF  (LERROR)  GO  TO  999
      IF  (IARGS .EQ. 2)  GO  TO  20
      IF  (IARGS .EQ. 1)  GO  TO  10
      CALL  TPSETM
      IF  (LEND)  GO  TO  30
      GO  TO  999
C
C         PROCESS A MULTI-LINE SET STATEMENT
C
   10 CONTINUE
      ICBP1(2)  =  ICBEOL + 1
      CALL  TPRDBL  (IK, IKYDIM, C, ICKDIM, ID, IDSDIM,
     A                           .TRUE., .FALSE., .TRUE., LERROR)
      IF  (LEND)    GO  TO  30
      IF  (LERROR)  GO  TO  999
      ICBP2(2)  =  ICB0 - 1
C
C         SET THE VALUE
C
   20 CONTINUE
      CALL  MMPUTV  (CBUFFR, ICBP1(1), ICBP2(1),
     A               CBUFFR, ICBP1(2), ICBP2(2))
      GO  TO  999
C
   30 CONTINUE
      CALL  IOERRM  (.FALSE.,
     A '('' ********   TPSET  - SET HAS NO MATCHING ENDSET'')')
C
  999 CONTINUE
      RETURN
      END
      SUBROUTINE  TPSETM
C
C-----------------------------------------------------------------------
C
C     FAMILY
C     ------
C     TEMPLATE PROCESSOR
C
C     PURPOSE
C     -------
C     TO PROCESS MULTILINE SET DIRECTIVES
C
C-----------------------------------------------------------------------
C
C         GLOBAL CONSTANTS
C
      CHARACTER*1         CA,      CBLANK,  CC,      CI,      CLEFT,
     A                    CMINUS,  CPLUS,   CPOINT,  CQUOTE,  CRIGHT,
     B                    CZ,      C0,      C9
      COMMON  / GLCOMC /  CA,      CBLANK,  CC,      CI,      CLEFT,
     A                    CMINUS,  CPLUS,   CPOINT,  CQUOTE,  CRIGHT,
     B                    CZ,      C0,      C9
C
C         INPUT / OUTPUT CONTROL INTERFACE
C
      CHARACTER*1         CBUFFR(2000)
      LOGICAL             LBREAK,  LFORT,   LISTI,   LISTO
      INTEGER             ICBADD,  ICBEND,  ICBEOL,  ICBSUB,  ICB0,
     A                    ICB1,    ICB2,    ICB3,    ICBDIM,  ICPLI,
     B                    ICPLO,   ILCTR,   ILNMBR,  ILPP,    IPAGE,
     A                    IUNITE,  IUNITI,  IUNITL,  IUNITO
      COMMON  / IOCOMC /  CBUFFR
      COMMON  / IOCOMI /  ICBADD,  ICBEND,  ICBEOL,  ICBSUB,  ICB0,
     A                    ICB1,    ICB2,    ICB3,    ICBDIM,  ICPLI,
     B                    ICPLO,   ILCTR,   ILNMBR,  ILPP,    IPAGE,
     A                    IUNITE,  IUNITI,  IUNITL,  IUNITO
      COMMON  / IOCOML /  LBREAK,  LFORT,   LISTI,   LISTO
C
C         MEMORY MANAGER INTERFACE
C
      CHARACTER*1         CSTORE(20000)
      INTEGER             IHASH(601),       ISTORE(6000)
      INTEGER             ICSDIM,  ICSP1,   ICSP2,   IHADIM,
     A                    ISFREE,  ISTDIM,  IS2HDC,  IS2HDS
      COMMON  / MMCOMC /  CSTORE
      COMMON  / MMCOMH /  IHASH
      COMMON  / MMCOMS /  ISTORE
      COMMON  / MMCOMI /  ICSDIM,  ICSP1,   ICSP2,   IHADIM,
     A                    ISFREE,  ISTDIM,  IS2HDC,  IS2HDS
C
C         MACRO PROCESSOR INTERFACE
C
      CHARACTER*1         CDIV,    CEOL,     CEOR,    CONC,    CSUB,
     A                    CTOP
      LOGICAL             LEMPTY,  LSUB
      COMMON  / MPCOMC /  CDIV,    CEOL,    CEOR,    CONC,    CSUB,
     A                    CTOP
      COMMON  / MPCOML /  LEMPTY,  LSUB
C
C         TEMPLATE PROCESSOR INTERFACE
C
      CHARACTER*1         CDIR,    CSTAR
      INTEGER             ICBP1(4),         ICBP2(4)
      INTEGER             ITOPDO,  IARGS, INESTD, INESTF
      LOGICAL             LCOL1,   LDIRL,   LEND,    LINITM,  L1TRIP
      COMMON  / TPCOMC /  CDIR,    CSTAR
      COMMON  / TPCOMI /  ICBP1,   ITOPDO,  IARGS,   ICBP2,   INESTD,
     B                    INESTF
      COMMON  / TPCOML /  LCOL1,   LDIRL,   LEND,    LINITM,  L1TRIP
C
C         LOCAL VARIABLES AND PARAMETERS
C
      INTEGER             IKYDIM,ICKDIM, IDSDIM, I
      CHARACTER*1         C(9)
      INTEGER             ID(5), IK(2)
      LOGICAL             LERROR, LSKIP
      EXTERNAL            MPLINE, TPCHKD, UTRDKY, TPSYNT, MMPUTV
      SAVE IKYDIM, IK, ICKDIM, C, IDSDIM, ID
      DATA
     A     IKYDIM, ICKDIM,
     B     IK(1),
     C     IK(2)
     D  /  2,      9,
     E     6,
     F     3  /
      DATA
     B     C(1),  C(2),  C(3),  C(4),  C(5),  C(6),
     C     C(7),  C(8),  C(9)
     D  /
     E     'E',   'N',   'D',   'S',   'E',   'T',
     F     'E',   'N',   'D'  /
      DATA                IDSDIM, ID(1), ID(2), ID(3), ID(4), ID(5)
     A                 /  5,      5,     4,     -6,    5,     7  /
C
      LSKIP  =  .TRUE.
C
   10 CONTINUE
          IF  (LSKIP)  ICBEOL  =  0
          CALL  MPLINE   (.TRUE.)
          CALL  TPCHKD
          IF  (.NOT. LDIRL)     GO  TO  20
          IF  (ICB1 .GT. ICB2)  GO  TO  30
          CALL  UTRDKY  (CBUFFR, ICB1, ICB2, IK, IKYDIM, C, ICKDIM, I)
          IF  (I .EQ. 1)  GO  TO  50
          IF  (I .EQ. 2)  GO  TO  60
      IF  (I .EQ. 3)  GO  TO  10
C
C         A TEXT LINE HAS BEEN ENCOUNTERED
C
   20 CONTINUE
      IF  (.NOT. LSKIP)  GO  TO  10
      CALL  TPSYNT  (ID, IDSDIM, .FALSE., LERROR)
      IF  (LERROR)  GO  TO  10
      IF  (IARGS .EQ. 2)  GO  TO  40
      ICBEOL  =  ICBP2(1)
      LSKIP   =  .FALSE.
      GO  TO  10
C
C         A DIRECTIVE PREFIX CHARACTER HAS BEEN
C         ENCOUNTERED ON A LINE BY ITSELF
C
   30 CONTINUE
      IF  (LSKIP)  GO  TO  10
      ICBP1(2)  =  ICBP2(1) + 1
      ICBP2(2)  =  ICB0 - 1
      LSKIP     =  .TRUE.
C
C         SAVE THE VALUE
C
   40 CONTINUE
      CALL  MMPUTV  (CBUFFR, ICBP1(1), ICBP2(1),
     A               CBUFFR, ICBP1(2), ICBP2(2))
      GO  TO  10
C
C         AN ENDSET DIRECTIVE HAS BEEN ENCOUNTERED
C
   50 CONTINUE
      IF  (.NOT. LSKIP)  GO  TO  10
      GO  TO  999
C
C         AN -END- DIRECTIVE HAS POSSIBLY BEEN ENCOUNTERED
C
   60 CONTINUE
      LEND  =  ICB1 .GT. ICB2
      IF  (.NOT. LEND)  GO  TO  10
C
  999 CONTINUE
      RETURN
      END
      SUBROUTINE  TPSYNT  (IDSYNT, IDSDIM, LSCAN, LERROR)
C
C-----------------------------------------------------------------------
C
C     FAMILY
C     ------
C     TEMPLATE PROCESSOR
C
C     PURPOSE
C     -------
C     TO CHECK A DIRECTIVE LINE FOR CORRECT SYNTAX
C
C     PARAMETERS
C     ----------
C     IDSYNT  -I-  CONTAINS THE DIRECTIVE SYNTAX PATTERN.
C                  THE VECTOR IDSYNT DESCRIBES THE TOKENS THAT
C                  ARE ALLOWED.  POSSIBLE VALUES OF IDSYNT(I):
C                          ABS(IDSYNT(I))     TOKEN
C                          --------------     -----
C                                1              (
C                                2              )
C                                3              ,
C                                4              =
C                                5              ID
C                                6              EXP
C                                7              EOL
C                                8              EOL
C                  WHEN IDSYNT(I) < 0, TWO THINGS CAN HAPPEN:
C                    - IF ABS(IDSYNT(I)) 'MATCHES' CURRENT TOKEN,
C                       SKIP TO IDSYNT(I+2) FOR NEXT MATCH.
C                    - IF NOT, SKIP TO IDSYNT(IDSYNT(I+1))
C                      FOR NEXT MATCH.
C
C     IDSDIM  -I-  DIMENSION OF IDSYNT
C     LSCAN   -I-  IF TRUE, DIRECTIVES ARE TO BE SCANNED
C                  BUT NOT EXECUTED
C     LERROR  -O-  TRUE IF THE DIRECTIVE HAS A SYNTAX ERROR
C
C-----------------------------------------------------------------------
C
C         GLOBAL CONSTANTS
C
      CHARACTER*1         CA,      CBLANK,  CC,      CI,      CLEFT,
     A                    CMINUS,  CPLUS,   CPOINT,  CQUOTE,  CRIGHT,
     B                    CZ,      C0,      C9
      COMMON  / GLCOMC /  CA,      CBLANK,  CC,      CI,      CLEFT,
     A                    CMINUS,  CPLUS,   CPOINT,  CQUOTE,  CRIGHT,
     B                    CZ,      C0,      C9
C
C         INPUT / OUTPUT CONTROL INTERFACE
C
      CHARACTER*1         CBUFFR(2000)
      LOGICAL             LBREAK,  LFORT,   LISTI,   LISTO
      INTEGER             ICBADD,  ICBEND,  ICBEOL,  ICBSUB,  ICB0,
     A                    ICB1,    ICB2,    ICB3,    ICBDIM,  ICPLI,
     B                    ICPLO,   ILCTR,   ILNMBR,  ILPP,    IPAGE,
     A                    IUNITE,  IUNITI,  IUNITL,  IUNITO
      COMMON  / IOCOMC /  CBUFFR
      COMMON  / IOCOMI /  ICBADD,  ICBEND,  ICBEOL,  ICBSUB,  ICB0,
     A                    ICB1,    ICB2,    ICB3,    ICBDIM,  ICPLI,
     B                    ICPLO,   ILCTR,   ILNMBR,  ILPP,    IPAGE,
     A                    IUNITE,  IUNITI,  IUNITL,  IUNITO
      COMMON  / IOCOML /  LBREAK,  LFORT,   LISTI,   LISTO
C
C         MEMORY MANAGER INTERFACE
C
      CHARACTER*1         CSTORE(20000)
      INTEGER             IHASH(601),       ISTORE(6000)
      INTEGER             ICSDIM,  ICSP1,   ICSP2,   IHADIM,
     A                    ISFREE,  ISTDIM,  IS2HDC,  IS2HDS
      COMMON  / MMCOMC /  CSTORE
      COMMON  / MMCOMH /  IHASH
      COMMON  / MMCOMS /  ISTORE
      COMMON  / MMCOMI /  ICSDIM,  ICSP1,   ICSP2,   IHADIM,
     A                    ISFREE,  ISTDIM,  IS2HDC,  IS2HDS
C
C         MACRO PROCESSOR INTERFACE
C
      CHARACTER*1         CDIV,    CEOL,     CEOR,    CONC,    CSUB,
     A                    CTOP
      LOGICAL             LEMPTY,  LSUB
      COMMON  / MPCOMC /  CDIV,    CEOL,    CEOR,    CONC,    CSUB,
     A                    CTOP
      COMMON  / MPCOML /  LEMPTY,  LSUB
C
C         TEMPLATE PROCESSOR INTERFACE
C
      CHARACTER*1         CDIR,    CSTAR
      INTEGER             ICBP1(4),         ICBP2(4)
      INTEGER             ITOPDO,  IARGS, INESTD, INESTF
      LOGICAL             LCOL1,   LDIRL,   LEND,    LINITM,  L1TRIP
      COMMON  / TPCOMC /  CDIR,    CSTAR
      COMMON  / TPCOMI /  ICBP1,   ITOPDO,  IARGS,   ICBP2,   INESTD,
     B                    INESTF
      COMMON  / TPCOML /  LCOL1,   LDIRL,   LEND,    LINITM,  L1TRIP
C
C         LOCAL VARIABLES AND PARAMETERS
C
      INTEGER             IDSDIM
      CHARACTER*1         C(4)
      INTEGER             IDSYNT(IDSDIM)
      LOGICAL             L, LEOL, LERROR, LSCAN
      INTEGER             I, ICV1, IJUMP, ICV2
      EXTERNAL            UTRDBL, UTRDNA, TPEXPR, IOERRM
      SAVE C
      DATA                C(1), C(2), C(3), C(4)
     A                 /  '(',  ')',  ',',  '='   /
C
      I       =  1
      IARGS   =  0
      LERROR  =  .FALSE.
C
C         DETERMINE WHICH TOKEN TO CHECK FOR
C
   10 CONTINUE
      ICV1    =  ICB1
      CALL  UTRDBL  (CBUFFR, ICB1, ICB2, LEOL)
      IJUMP   =  IABS(IDSYNT(I))
      GO  TO  (20, 20, 20, 20, 30, 40, 50, 50),  IJUMP
C
C         CHECK FOR DELIMITERS AND SEPARATERS
C
   20 CONTINUE
      IF  (LEOL)  GO  TO  80
      IF  (C(IJUMP) .NE. CBUFFR(ICB1))  GO  TO  80
      ICB1    =  ICB1 + 1
      GO  TO  70
C
C         CHECK FOR A NAME
C
   30 CONTINUE
      IF  (LEOL)  GO  TO  80
      CALL  UTRDNA  (CBUFFR, ICB1, ICB2, ICV1, ICV2, L)
      IF  (L)     GO  TO  80
      GO  TO  60
C
C         CHECK FOR AN EXPRESSION
C
   40 CONTINUE
      IF  (LEOL)  GO  TO  80
      CALL  TPEXPR  (ICV1, ICV2, LSCAN, L)
      IF  (L)     GO  TO  80
      GO  TO  60
C
C         CHECK FOR END OF LINE
C
   50 CONTINUE
      IF  (LEOL)  GO  TO  999
      IF  (IJUMP .NE. 8)  GO  TO  80
      ICV2  =  ICBEOL
C
   60 CONTINUE
      IARGS  =  IARGS + 1
      IF  (LSCAN)  GO  TO  70
      ICBP1(IARGS)  =  ICV1
      ICBP2(IARGS)  =  ICV2
C
   70 CONTINUE
      IF  (IDSYNT(I) .LT. 0)  I  =  I + 1
      I   =  I + 1
      IF  (I .LE. IDSDIM)  GO  TO  10
      GO  TO  999
C
C         IF THERE IS AN ALTERNATE SYNTAX FOR THIS STATEMENT
C         THEN TRY IT, OTHERWISE PRINT AN ERROR MESSAGE
C
   80 CONTINUE
      IF  (IDSYNT(I).GT. 0)  GO  TO  90
      I  =  IDSYNT(I+1)
      IF  (I .LE. IDSDIM)  GO  TO  10
C
C         ERROR EXITS
C
   90 CONTINUE
      LERROR  =  .TRUE.
      GO  TO  (100, 110, 120, 130, 140, 150, 160, 160), IJUMP
C
  100 CONTINUE
      CALL  IOERRM  (.FALSE.,
     A '('' ********   TPSYNT - LEFT PARENTHESIS EXPECTED'')')
      GO  TO  999
C
  110 CONTINUE
      CALL  IOERRM  (.FALSE.,
     A '('' ********   TPSYNT - RIGHT PARENTHESIS EXPECTED'')')
      GO  TO  999
C
  120 CONTINUE
      CALL  IOERRM  (.FALSE.,
     A '('' ********   TPSYNT - COMMA EXPECTED'')')
      GO  TO  999
C
  130 CONTINUE
      CALL  IOERRM  (.FALSE.,
     A '('' ********   TPSYNT - EQUALS SIGN EXPECTED'')')
      GO  TO  999
C
  140 CONTINUE
      CALL  IOERRM  (.FALSE.,
     A '('' ********   TPSYNT - VARIABLE EXPECTED'')')
      GO  TO  999
C
  150 CONTINUE
      CALL  IOERRM  (.FALSE.,
     A '('' ********   TPSYNT - MISSING OR UNRECOGNIZED EXPRESSION'')')
      GO  TO  999
C
  160 CONTINUE
      CALL  IOERRM  (.FALSE.,
     A '('' ********   TPSYNT - ILLEGAL CHARACTERS AT END OF LINE'')')
C
  999 CONTINUE
      RETURN
      END
      SUBROUTINE  UTBLDN  (CPREFX, CROOT, ICR1, ICR2, ISUFFX,
     A                     CNAME, ICN1, ICN2, ICNDIM, LERROR)
C
C-----------------------------------------------------------------------
C
C     FAMILY
C     ------
C     UTILITY
C
C     PURPOSE
C     -------
C     TO BUILD A NAME GIVEN A PREFIX, ROOT, AND SUFFIX
C
C     PARAMETERS
C     ----------
C     CPREFX  -I-  A ONE CHARACTER PREFIX
C     CROOT   -I-  ROOT OF THE NAME
C     ICR1    -I-  INDEX OF THE FIRST CHARACTER IN THE ROOT
C     ICR2    -I-  INDEX OF THE LAST CHARACTER IN THE ROOT
C     ISUFFX  -I-  INTEGER SUFFIX
C     CNAME   -O-  THE NAME
C     ICN1    -I-  INDEX OF THE FIRST CHARACTER IN THE NAME
C     ICN2    -O-  INDEX OF THE LAST CHARACTER IN THE NAME
C     ICNDIM  -I-  DIMENSION OF CNAME
C     LERROR  -O-  TRUE IF AN ERROR OCCURED
C
C-----------------------------------------------------------------------
C
C         GLOBAL CONSTANTS
C
      CHARACTER*1         CA,      CBLANK,  CC,      CI,      CLEFT,
     A                    CMINUS,  CPLUS,   CPOINT,  CQUOTE,  CRIGHT,
     B                    CZ,      C0,      C9
      COMMON  / GLCOMC /  CA,      CBLANK,  CC,      CI,      CLEFT,
     A                    CMINUS,  CPLUS,   CPOINT,  CQUOTE,  CRIGHT,
     B                    CZ,      C0,      C9
C
C         LOCAL VARIABLES AND PARAMETERS
C
      INTEGER             ICR1, ICR2, ISUFFX, ICN1, ICN2, ICNDIM
      CHARACTER*(*)      CNAME(ICNDIM), CPREFX, CROOT(ICR2)
      LOGICAL             LERROR
      INTEGER             I
      EXTERNAL            UTCVIC
C
      LERROR       =  ICN1 + ICR2 - ICR1 + 1 .GT. ICNDIM
      IF  (LERROR)  GO  TO  999
      ICN2         =  ICN1
      CNAME(ICN2)  =  CPREFX
      IF  (ICR1 .GT. ICR2)  GO  TO  20
      DO  10  I=ICR1,ICR2
          ICN2         =  ICN2 + 1
          CNAME(ICN2)  =  CROOT(I)
   10 CONTINUE
   20 CONTINUE
      IF  (ISUFFX .LT. 0)  GO  TO  999
      I           =  ICN2 + 1
      CALL  UTCVIC  (CNAME, I, ICN2, ICNDIM, ISUFFX, LERROR)
C
  999 CONTINUE
      RETURN
      END
      SUBROUTINE  UTCVCI  (CLINE, ICL1, ICL2, IVALUE, LERROR)
C
C-----------------------------------------------------------------------
C
C     FAMILY
C     ------
C     UTILITY
C
C     PURPOSE
C     -------
C     TO CONVERT A CHARACTER STRING INTO AN INTEGER
C
C     PARAMETERS
C     ----------
C     CLINE   -I-  STRING TO BE CONVERTED
C     ICL1    -I-  INDEX OF THE FIRST CHARACTER IN THE STRING
C     ICL2    -I-  INDEX OF THE LAST CHARACTER IN THE STRING
C     IVALUE  -O-  INTEGER RESULT
C     LERROR  -O-  TRUE IF AN ERROR OCCURED
C
C-----------------------------------------------------------------------
C
C         GLOBAL CONSTANTS
C
      CHARACTER*1         CA,      CBLANK,  CC,      CI,      CLEFT,
     A                    CMINUS,  CPLUS,   CPOINT,  CQUOTE,  CRIGHT,
     B                    CZ,      C0,      C9
      COMMON  / GLCOMC /  CA,      CBLANK,  CC,      CI,      CLEFT,
     A                    CMINUS,  CPLUS,   CPOINT,  CQUOTE,  CRIGHT,
     B                    CZ,      C0,      C9
C
C         LOCAL VARIABLES AND PARAMETERS
C
      INTEGER             ICL1, ICL2, IVALUE
      CHARACTER*(*)      CLINE(ICL2)
      CHARACTER*1         C(10), CLINEI
      LOGICAL             LERROR, LMINUS
      INTEGER             I, I1, IC
      EXTERNAL            UTRDBL
      SAVE C
      DATA                C(1),  C(2),  C(3),  C(4),  C(5),
     A                    C(6),  C(7),  C(8),  C(9),  C(10)
     B                 /  '0',   '1',   '2',   '3',   '4',
     C                    '5',   '6',   '7',   '8',   '9'  /
C
      IVALUE  =  0
      I       =  ICL1
      CALL  UTRDBL  (CLINE, I, ICL2, LERROR)
      IF  (LERROR)  GO  TO  999
      LMINUS  =  CLINE(I) .EQ. CMINUS
      IF  ((.NOT. LMINUS) .AND. (CLINE(I) .NE. CPLUS))  GO  TO  10
          I  =  I + 1
          CALL  UTRDBL  (CLINE, I, ICL2, LERROR)
          IF  (LERROR)  GO  TO  999
C
   10 CONTINUE
      I1  =  I
      DO  40  I=I1,ICL2
          CLINEI  =  CLINE(I)
          DO  20  IC=1,10
              IF  (CLINEI .EQ. C(IC))  GO  TO  30
   20     CONTINUE
          IF  (I .GT. I1)  GO  TO  50
          GO  TO  999
   30     CONTINUE
          IVALUE  =  IVALUE*10 + IC - 1
   40 CONTINUE
C
   50 CONTINUE
      LERROR  =  .FALSE.
      IF  (LMINUS)  IVALUE  =  -IVALUE
C
  999 CONTINUE
      RETURN
      END
      SUBROUTINE  UTCVCL  (CLINE, ICL1, ICL2, LVALUE, LERROR)
C
C-----------------------------------------------------------------------
C
C     FAMILY
C     ------
C     UTILITY
C
C     PURPOSE
C     -------
C     TO CONVERT A CHARACTER STRING TO A LOGICAL VALUE
C
C     PARAMETERS
C     ----------
C     CLINE   -I-  STRING TO BE CONVERTED
C     ICL1    -I-  INDEX OF THE FIRST CHARACTER IN THE STRING
C     ICL2    -I-  INDEX OF THE LAST CHARACTER IN THE STRING
C     LVALUE  -O-  THE LOGICAL RESULT
C     LERROR  -I-  TRUE IF AN ERROR OCCURED
C
C-----------------------------------------------------------------------
C
C         GLOBAL CONSTANTS
C
      CHARACTER*1         CA,      CBLANK,  CC,      CI,      CLEFT,
     A                    CMINUS,  CPLUS,   CPOINT,  CQUOTE,  CRIGHT,
     B                    CZ,      C0,      C9
      COMMON  / GLCOMC /  CA,      CBLANK,  CC,      CI,      CLEFT,
     A                    CMINUS,  CPLUS,   CPOINT,  CQUOTE,  CRIGHT,
     B                    CZ,      C0,      C9
C
C         LOCAL VARIABLES AND PARAMETERS
C
      INTEGER             ICL1, ICL2
      CHARACTER*(*)      CLINE(ICL2)
      CHARACTER*1         C(13)
      INTEGER             IK(2)
      LOGICAL             LERROR, LV(3), LVALUE
      INTEGER             IKYDIM, ICKDIM, I
      EXTERNAL            UTRDKY
      SAVE IKYDIM, IK, ICKDIM, C, LV
      DATA
     A     IKYDIM, ICKDIM,
     B     IK(1),
     C     IK(2)
     D  /  6,      13,
     E     6,
     F     7  /
      DATA
     B     C(1),  C(2),  C(3),  C(4),  C(5),  C(6),
     C     C(7),  C(8),  C(9),  C(10), C(11), C(12), C(13)
     D  /
     E     '.',   'T',   'R',   'U',   'E',   '.',
     F     '.',   'F',   'A',   'L',   'S',   'E',   '.'  /
      DATA
     A     LV(1),   LV(2),   LV(3)
     B  /  .TRUE.,  .FALSE., .TRUE.  /
C
      LERROR  =  .TRUE.
      IF  (ICL1 .GT. ICL2)  GO  TO  999
      DO  10  I=ICL1,ICL2
          IF  (CLINE(I) .NE. CBLANK)  GO  TO  20
   10 CONTINUE
      GO  TO  999
C
   20 CONTINUE
      CALL  UTRDKY  (CLINE, ICL1, ICL2, IK, IKYDIM, C, ICKDIM, I)
      LERROR  =  I .GT. IKYDIM
      LVALUE  =  LV(I)
C
  999 CONTINUE
      RETURN
      END
      SUBROUTINE  UTCVIC  (CLINE, ICL1, ICL2, ICLDIM, IVALUE, LERROR)
C
C-----------------------------------------------------------------------
C
C     FAMILY
C     ------
C     UTILITY
C
C     PURPOSE
C     -------
C     TO CONVERT AN INTEGER INTO A CHARACTER STRING
C
C     PARAMETERS
C     ----------
C     CLINE   -O-  STRING RESULT
C     ICL1    -I-  INDEX OF THE FIRST CHARACTER IN THE STRING
C     ICL2    -O-  INDEX OF THE LAST CHARACTER IN THE STRING
C     ICLDIM  -I-  DIMENSION OF CLINE
C     IVALUE  -I-  INTEGER TO BE CONVERTED
C     LERROR  -O-  TRUE IF AN ERROR OCCURED
C
C-----------------------------------------------------------------------
C
C         GLOBAL CONSTANTS
C
      CHARACTER*1         CA,      CBLANK,  CC,      CI,      CLEFT,
     A                    CMINUS,  CPLUS,   CPOINT,  CQUOTE,  CRIGHT,
     B                    CZ,      C0,      C9
      COMMON  / GLCOMC /  CA,      CBLANK,  CC,      CI,      CLEFT,
     A                    CMINUS,  CPLUS,   CPOINT,  CQUOTE,  CRIGHT,
     B                    CZ,      C0,      C9
C
C         LOCAL VARIABLES AND PARAMETERS
C
      INTEGER             ICL1, ICL2, ICLDIM, IVALUE
      CHARACTER*(*)      CLINE(ICLDIM)
      CHARACTER*1         C(10), CTEMP
      LOGICAL             LERROR
      INTEGER             I1, I2, ICL2MD
      SAVE C
      DATA                C(1),  C(2),  C(3),  C(4),  C(5),
     A                    C(6),  C(7),  C(8),  C(9),  C(10)
     B                 /  '0',   '1',   '2',   '3',   '4',
     C                    '5',   '6',   '7',   '8',   '9'  /
C
      I1      =  IABS(IVALUE)
      LERROR  =  .TRUE.
      ICL2     =  ICL1 - 1
C
C         CONVERT AND THEN REMOVE THE LEAST SIGNIFICANT DIGITS FIRST
C
   10 CONTINUE
          I2          =  I1
          I1          =  I1 / 10
          I2          =  I2 - I1*10
          ICL2         =  ICL2 + 1
          IF  (ICL2 .GT. ICLDIM)  GO  TO  999
          CLINE(ICL2)  =  C(I2+1)
      IF  (I1 .GT. 0)  GO  TO  10
C
C         IF NECESSARY, ADD THE MINUS SIGN
C
      IF  (IVALUE .GE. 0)  GO  TO  20
          ICL2         =  ICL2 + 1
          IF  (ICL2 .GT. ICLDIM)  GO  TO  999
          CLINE(ICL2)  =  CMINUS
C
C         REVERSE THE STRING TO PUT THE DIGITS IN THE PROPER ORDER
C
   20 CONTINUE
      LERROR  =  .FALSE.
      IF  (ICL1 .GE. ICL2)  GO  TO  999
      ICL2MD  =  (ICL1 + ICL2 - 1) / 2
      I2      =  ICL2
      DO  30  I1=ICL1,ICL2MD
          CTEMP      =  CLINE(I1)
          CLINE(I1)  =  CLINE(I2)
          CLINE(I2)  =  CTEMP
          I2         =  I2 - 1
   30 CONTINUE
C
  999 CONTINUE
      RETURN
      END
      SUBROUTINE  UTCVLC  (CLINE, ICL1, ICL2, ICLDIM, LVALUE, LERROR)
C
C-----------------------------------------------------------------------
C
C     FAMILY
C     ------
C     UTILITY
C
C     PURPOSE
C     -------
C     TO CONVERT A LOGICAL VALUE TO A CHARACTER
C
C     PARAMETERS
C     ----------
C     CLINE   -O-  STRING RESULT
C     ICL1    -I-  INDEX OF THE FIRST CHARACTER IN THE STRING
C     ICL2    -O-  INDEX OF THE LAST CHARACTER IN THE STRING
C     ICLDIM  -I-  DIMENSION OF CLINE
C     LVALUE  -I-  LOGICAL VALUE TO BE CONVERTED
C     LERROR  -O-  TRUE IF AN ERROR OCCURED
C
C-----------------------------------------------------------------------
C
C         GLOBAL CONSTANTS
C
      CHARACTER*1         CA,      CBLANK,  CC,      CI,      CLEFT,
     A                    CMINUS,  CPLUS,   CPOINT,  CQUOTE,  CRIGHT,
     B                    CZ,      C0,      C9
      COMMON  / GLCOMC /  CA,      CBLANK,  CC,      CI,      CLEFT,
     A                    CMINUS,  CPLUS,   CPOINT,  CQUOTE,  CRIGHT,
     B                    CZ,      C0,      C9
C
C         LOCAL VARIABLES AND PARAMETERS
C
      INTEGER             ICL1, ICL2, ICLDIM
      CHARACTER*(*)      CLINE(ICLDIM)
      CHARACTER*1         CF(7), CT(6)
      LOGICAL             LERROR, LVALUE
      INTEGER             ICFDIM, ICTDIM, I
      SAVE ICFDIM, CF, ICTDIM, CT
      DATA
     A     ICFDIM,
     B     ICTDIM
     C  /  7,
     D     6  /
      DATA
     A     CF(1),  CF(2),  CF(3),  CF(4),  CF(5),  CF(6),  CF(7),
     B     CT(1),  CT(2),  CT(3),  CT(4),  CT(5),  CT(6)
     C  /  '.',    'F',    'A',    'L',    'S',    'E',    '.',
     D     '.',    'T',    'R',    'U',    'E',    '.'  /
C
      ICL2    =  ICL1 - 1
      IF  (LVALUE)  GO  TO  20
      LERROR  =  (ICL2 + ICFDIM) .GT. ICLDIM
      IF  (LERROR)  GO  TO  999
      DO  10  I=1,ICFDIM
          ICL2         =  ICL2 + 1
          CLINE(ICL2)  =  CF(I)
   10 CONTINUE
      GO  TO  999
C
   20 CONTINUE
      LERROR  =  (ICL2 + ICTDIM) .GT. ICLDIM
      IF  (LERROR)  GO  TO  999
      DO  30  I=1,ICTDIM
          ICL2         =  ICL2 + 1
          CLINE(ICL2)  =  CT(I)
   30 CONTINUE
C
  999 CONTINUE
      RETURN
      END
      SUBROUTINE  UTCVNI  (CNAME, ICN1, ICN2, INAME, LERROR)
C-----------------------------------------------------------------------
C
C     FAMILY
C     ------
C     UTILITY
C
C     PURPOSE
C     -------
C     TO CONVERT (HASH) A NAME INTO AN INTEGER
C
C     PARAMETERS
C     ----------
C     CNAME   -I-  THE NAME TO BE HASHED
C     ICN1    -I-  INDEX OF THE FIRST CHARACTER IN THE NAME
C     ICN2    -I-  INDEX OF THE LAST CHARACTER IN THE NAME
C     INAME   -O-  INTEGER RESULT
C     LERROR  -O-  TRUE IF AN ERROR OCCURED
C
C-----------------------------------------------------------------------
C
C         GLOBAL CONSTANTS
C
      CHARACTER*1         CA,      CBLANK,  CC,      CI,      CLEFT,
     A                    CMINUS,  CPLUS,   CPOINT,  CQUOTE,  CRIGHT,
     B                    CZ,      C0,      C9
      COMMON  / GLCOMC /  CA,      CBLANK,  CC,      CI,      CLEFT,
     A                    CMINUS,  CPLUS,   CPOINT,  CQUOTE,  CRIGHT,
     B                    CZ,      C0,      C9
C
C         LOCAL VARIABLES AND PARAMETERS
C
      INTEGER                ICN1, ICN2, INAME
      CHARACTER*(*)         CNAME(ICN2)
      INTEGER ICHDIM, ICIDIM, I, ICNMIN, ICN, ICH
      INTEGER             IC(6)
      LOGICAL             LERROR
      CHARACTER*1         C(48), CNAMEI
      SAVE ICIDIM, IC, ICHDIM, C
      DATA
     A              ICHDIM,
     B              C(1),   C(2),   C(3),   C(4),   C(5),   C(6),
     C              C(7),   C(8),   C(9),   C(10),  C(11),  C(12),
     D              C(13),  C(14),  C(15),  C(16),  C(17),  C(18),
     E              C(19),  C(20),  C(21),  C(22),  C(23),  C(24),
     F              C(25),  C(26),  C(27),  C(28),  C(29),  C(30),
     G              C(31),  C(32),  C(33),  C(34),  C(35),  C(36),
     H              C(37),  C(38),  C(39),  C(40),  C(41),  C(42),
     I              C(43),  C(44),  C(45),  C(46),  C(47),  C(48)
     J           /  48,
     K              'A',    'B',    'C',    'D',    'E',    'F',
     L              'G',    'H',    'I',    'J',    'K',    'L',
     M              'M',    'N',    'O',    'P',    'Q',    'R',
     N              'S',    'T',    'U',    'V',    'W',    'X',
     O              'Y',    'Z',    '0',    '1',    '2',    '3',
     P              '4',    '5',    '6',    '7',    '8',    '9',
     Q              '+',    '-',    '*',    ',',    '=',    '(',
     R              ')',    '.',    ',',    '''',   '$',    ' ' /
      DATA
     A              ICIDIM,
     B              IC(1),  IC(2),  IC(3),  IC(4),  IC(5),  IC(6)
     C           /  6,
     D              61,     1,      47,     61,     1,      47  /
C
      LERROR  =  ICN1 .GT. ICN2
      IF  (LERROR)  GO  TO  999
      I       =  0
      INAME   =  0
      ICNMIN  =  MIN0(ICN2, ICN1+ICIDIM-1)
C
C
      DO  30  ICN=ICN1,ICNMIN
          CNAMEI  =  CNAME(ICN)
          DO  10  ICH=1,ICHDIM
              IF  (CNAMEI .EQ. C(ICH))  GO  TO  20
   10     CONTINUE
          ICH     =  ICHDIM + 1
   20     CONTINUE
          I       =  I + 1
          INAME   =  INAME + IC(I)*ICH
   30 CONTINUE
C
  999 CONTINUE
      RETURN
      END
      SUBROUTINE  UTRDBL  (CLINE, ICL1, ICL2, LEOL)
C
C-----------------------------------------------------------------------
C
C     FAMILY
C     ------
C     UTILITY
C
C     PURPOSE
C     -------
C     TO READ (SKIP) BLANKS IN A LINE
C
C     PARAMETERS
C     ----------
C     CLINE   -I-  LINE OF CHARACTERS
C     ICL1    -I-  INDEX OF THE FIRST CHARACTER IN THE LINE
C     ICL2    -I-  INDEX OF THE LAST CHARACTER IN THE LINE
C     LEOL    -O-  TRUE IF THE END OF THE LINE WAS REACHED
C
C-----------------------------------------------------------------------
C
C         GLOBAL CONSTANTS
C
      CHARACTER*1         CA,      CBLANK,  CC,      CI,      CLEFT,
     A                    CMINUS,  CPLUS,   CPOINT,  CQUOTE,  CRIGHT,
     B                    CZ,      C0,      C9
      COMMON  / GLCOMC /  CA,      CBLANK,  CC,      CI,      CLEFT,
     A                    CMINUS,  CPLUS,   CPOINT,  CQUOTE,  CRIGHT,
     B                    CZ,      C0,      C9
C
C         LOCAL VARIABLES AND PARAMETERS
C
      INTEGER             ICL1, ICL2
      CHARACTER*(*)         CLINE(ICL2)
      LOGICAL             LEOL
      INTEGER             I
C
      IF  (ICL1 .GT. ICL2)  GO  TO  20
C
      DO  10  I=ICL1,ICL2
          IF  (CBLANK .NE. CLINE(I))  GO  TO  30
   10 CONTINUE
C
   20 CONTINUE
      ICL1  =  ICL2 + 1
      LEOL  =  .TRUE.
      GO  TO  999
C
   30 CONTINUE
      ICL1  =  I
      LEOL  =  .FALSE.
C
  999 CONTINUE
      RETURN
      END
      SUBROUTINE  UTRDKY  (CLINE, ICL1, ICL2, IKEY, IKYDIM,
     A                                        CKEY, ICKDIM, IK)
C
C-----------------------------------------------------------------------
C
C     FAMILY
C     ------
C     UTILITY
C
C     PURPOSE
C     -------
C     TO MATCH CHARACTERS WITH ONE OF A GIVEN SET OF KEYS
C
C     PARAMETERS
C     ----------
C     CLINE   -I-  LINE OF CHARACTERS
C     ICL1    -I-  INDEX OF THE FIRST CHARACTER IN THE LINE
C     ICL2    -I-  INDEX OF THE LAST CHARACTER IN THE LINE
C     IKEY    -I-  CONTAINS THE LENGTH OF EACH KEY
C     IKYDIM  -I-  NUMBER OF KEYS
C     CKEY    -I-  CONTAINS THE KEYS
C     CKYDIM  -I-  DIMENSION OF CKEY
C     IK      -O-  NUMBER OF THE MATCHED KEY
C
C-----------------------------------------------------------------------
C
C         GLOBAL CONSTANTS
C
      CHARACTER*1         CA,      CBLANK,  CC,      CI,      CLEFT,
     A                    CMINUS,  CPLUS,   CPOINT,  CQUOTE,  CRIGHT,
     B                    CZ,      C0,      C9
      COMMON  / GLCOMC /  CA,      CBLANK,  CC,      CI,      CLEFT,
     A                    CMINUS,  CPLUS,   CPOINT,  CQUOTE,  CRIGHT,
     B                    CZ,      C0,      C9
C
C         LOCAL VARIABLES AND PARAMETERS
C
      INTEGER             ICL1, ICL2, IKYDIM, ICKDIM, IK
      CHARACTER*(*)      CLINE(ICL2), CKEY(ICKDIM)
      INTEGER             IKEY(IKYDIM)
      INTEGER             ICK2, ICLDIF, ICK1, I, ICK
C
      IF  (ICL1 .GT. ICL2)  GO  TO  30
      ICK2    =  0
      ICLDIF  =  ICL2 - ICL1 + 1
C
      DO  20  IK=1,IKYDIM
          ICK1   =  ICK2 + 1
          ICK2   =  ICK2 + IKEY(IK)
          IF  (ICLDIF .LT. IKEY(IK))  GO  TO  20
          I      =  ICL1
          DO  10  ICK=ICK1,ICK2
              IF  (CLINE(I) .NE. CKEY(ICK))   GO  TO  20
              I  =  I + 1
   10     CONTINUE
          GO  TO  40
   20 CONTINUE
C
   30 CONTINUE
      IK    =  IKYDIM + 1
      GO  TO  999
C
   40 CONTINUE
      ICL1  =  ICL1 + IKEY(IK)
C
  999 CONTINUE
      RETURN
      END
      SUBROUTINE  UTRDNA  (CLINE, ICL1, ICL2, ICL1NA, ICL2NA, LERROR)
C
C-----------------------------------------------------------------------
C
C     FAMILY
C     ------
C     UTILITY
C
C     PURPOSE
C     -------
C     TO READ A NAME ON A LINE
C
C     PARAMETERS
C     ----------
C     CLINE   -I-  LINE OF CHARACTERS
C     ICL1    -I-  INDEX OF THE FIRST CHARACTER IN THE LINE
C     ICL2    -I-  INDEX OF THE LAST CHARACTER IN THE LINE
C     ICL1NA  -O-  INDEX OF THE FIRST CHARACTER IN THE NAME
C     ICL2NA  -O-  INDEX OF THE LAST CHARACTER IN THE NAME
C     LERROR  -O-  TRUE IF AN ERROR OCCURED
C
C-----------------------------------------------------------------------
C
C         GLOBAL CONSTANTS
C
      CHARACTER*1         CA,      CBLANK,  CC,      CI,      CLEFT,
     A                    CMINUS,  CPLUS,   CPOINT,  CQUOTE,  CRIGHT,
     B                    CZ,      C0,      C9
      COMMON  / GLCOMC /  CA,      CBLANK,  CC,      CI,      CLEFT,
     A                    CMINUS,  CPLUS,   CPOINT,  CQUOTE,  CRIGHT,
     B                    CZ,      C0,      C9
C
C         LOCAL VARIABLES AND PARAMETERS
C
      INTEGER             ICL1, ICL2, ICL1NA, ICL2NA
      CHARACTER*(*)      CLINE(ICL2)
      LOGICAL             LERROR
      INTEGER             I
C
      LERROR  =  .TRUE.
      IF  (ICL1 .GT. ICL2)  GO  TO  999
C
      DO  10  I=ICL1,ICL2
          IF  (.NOT. ((LLE(CA,CLINE(I))
     A         .AND.   LLE(CLINE(I),CZ))
     B         .OR.   (LLE(C0,CLINE(I))
     C         .AND.   LLE(CLINE(I),C9))))  GO  TO  20
   10 CONTINUE
C
      I       =  ICL2 + 1
C
   20 CONTINUE
      IF  (.NOT. (LLE(CA,CLINE(ICL1))
     A     .AND.  LLE(CLINE(ICL1),CZ)))  GO  TO  999
      ICL1NA   =  ICL1
      ICL1     =  I
      ICL2NA  =  I - 1
      LERROR  =  .FALSE.
C
  999 CONTINUE
      RETURN
      END
      SUBROUTINE  UTRDNU  (CLINE, ICL1, ICL2, ICL1NU, ICL2NU, LERROR)
C
C-----------------------------------------------------------------------
C
C     FAMILY
C     ------
C     UTILITY
C
C     PURPOSE
C     -------
C     TO READ A NUMBER ON A LINE
C
C     PARAMETERS
C     ----------
C     CLINE   -I-  LINE OF CHARACTERS
C     ICL1    -I-  INDEX OF THE FIRST CHARACTER IN THE LINE
C     ICL2    -I-  INDEX OF THE LAST CHARACTER IN THE LINE
C     ICL1NU  -O-  INDEX OF THE FIRST CHARACTER IN THE NUMBER
C     ICL2NU  -O-  INDEX OF THE LAST CHARACTER IN THE NUMBER
C     LERROR  -O-  TRUE IF AN ERROR OCCURED
C
C-----------------------------------------------------------------------
C
C         GLOBAL CONSTANTS
C
      CHARACTER*1         CA,      CBLANK,  CC,      CI,      CLEFT,
     A                    CMINUS,  CPLUS,   CPOINT,  CQUOTE,  CRIGHT,
     B                    CZ,      C0,      C9
      COMMON  / GLCOMC /  CA,      CBLANK,  CC,      CI,      CLEFT,
     A                    CMINUS,  CPLUS,   CPOINT,  CQUOTE,  CRIGHT,
     B                    CZ,      C0,      C9
C
C         LOCAL VARIABLES AND PARAMETERS
C
      INTEGER             ICL1, ICL2, ICL1NU, ICL2NU
      CHARACTER*(*)      CLINE(ICL2)
      LOGICAL             LERROR
      INTEGER             I, ICL
      EXTERNAL            UTRDBL
C
      LERROR  =  ICL1 .GT. ICL2
      IF  (LERROR)  GO  TO  999
      I       =  ICL1
      IF  ((CLINE(I) .NE. CMINUS)
     A    .AND. (CLINE(I) .NE. CPLUS))  GO  TO  10
          I   =  I + 1
          CALL  UTRDBL  (CLINE, I, ICL2, LERROR)
          IF  (LERROR)  GO  TO  999
C
   10 CONTINUE
      ICL     =  I
      DO  20  I=ICL,ICL2
          IF  (.NOT. (LLE(C0,CLINE(I))
     A         .AND.  LLE(CLINE(I),C9)))  GO  TO  30
   20 CONTINUE
C
      I       =  ICL2 + 1
C
   30 CONTINUE
      ICL1NU  =  ICL1
      ICL1    =  I
      ICL2NU  =  I - 1
      LERROR  =  ICL1NU .GT. ICL2NU
C
  999 CONTINUE
      RETURN
      END
      SUBROUTINE  UTRDQS  (CLINE, ICL1, ICL2, ICL1QS, ICL2QS, LERROR)
C
C-----------------------------------------------------------------------
C
C     FAMILY
C     ------
C     UTILITY
C
C     PURPOSE
C     -------
C     TO READ A QUOTED STRING ON A LINE
C
C     PARAMETERS
C     ----------
C     CLINE   -I-  LINE OF CHARACTERS
C     ICL1    -I-  INDEX OF THE FIRST CHARACTER IN THE LINE
C     ICL2    -I-  INDEX OF THE LAST CHARACTER IN THE LINE
C     ICL1QS  -O-  INDEX OF THE FIRST CHARACTER IN THE STRING
C     ICL2QS  -O-  INDEX OF THE LAST CHARACTER IN THE STRING
C     LERROR  -O-  TRUE IF AN ERROR OCCURED
C
C-----------------------------------------------------------------------
C
C         GLOBAL CONSTANTS
C
      CHARACTER*1         CA,      CBLANK,  CC,      CI,      CLEFT,
     A                    CMINUS,  CPLUS,   CPOINT,  CQUOTE,  CRIGHT,
     B                    CZ,      C0,      C9
      COMMON  / GLCOMC /  CA,      CBLANK,  CC,      CI,      CLEFT,
     A                    CMINUS,  CPLUS,   CPOINT,  CQUOTE,  CRIGHT,
     B                    CZ,      C0,      C9
C
C         LOCAL VARIABLES AND PARAMETERS
C
      INTEGER             ICL1, ICL2, ICL1QS, ICL2QS
      CHARACTER*(*)      CLINE(ICL2)
      CHARACTER*1         CQTEMP
      LOGICAL             LERROR
      INTEGER             I
C
      LERROR  =  .TRUE.
      ICL1QS  =  ICL1 + 1
      IF  (ICL1QS .GT. ICL2)  GO  TO  999
      CQTEMP  =  CLINE(ICL1)
C
      DO  10  I=ICL1QS,ICL2
          IF  (CLINE(I) .EQ. CQTEMP)  GO  TO  20
   10 CONTINUE
C
      GO  TO  999
C
   20 CONTINUE
      ICL1    =  I + 1
      ICL2QS  =  I - 1
      LERROR  =  .FALSE.
C
  999 CONTINUE
      RETURN
      END
