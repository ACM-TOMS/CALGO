      SUBROUTINE ZZLINK ( N, X,  F, G,  ACC, RELF, RELG, STATUS, SUBR,
     -          TRACUN,  TRACES, NTR,   PRINTL, MAX,    DERVMD, USER,
     -                   NU,     RWORK, LR, IWORK,  LI, IW, RW, DW   )

C## A R G U M E N T S:

      INTEGER   N, LR, LI, NU, SUBR, STATUS
      INTEGER   TRACUN, NTR,   MAX   , PRINTL, DERVMD
      INTEGER   IWORK(LI), IW(*)
      LOGICAL   TRACES(NTR), RELF, RELG

      DOUBLE PRECISION F, ACC, X(N), G(N), RWORK(LR), USER(NU)
C!!!! REAL             F, ACC, X(N), G(N), RWORK(LR), USER(NU)
      REAL             RW(*)
      DOUBLE PRECISION DW(*)

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
C>RCS $HEADER: LNK.F,V 1.11 91/11/20 10:53:00 BUCKLEY EXP $
C>RCS $LOG:     LNK.F,V $
C>RCS REVISION 1.11  91/11/20  10:53:00  BUCKLEY
C>RCS FINAL SUBMISSION TO TOMS
C>RCS
C>RCS REVISION 1.10  90/08/06  16:25:31  BUCKLEY
C>RCS MODIFIED WITH FIX TO BBLNIR TO CHECK ALL ARGUMENTS.
C>RCS
C>RCS REVISION 1.9  89/06/30  13:39:45  BUCKLEY
C>RCS PREPARING SUBMITTED VERSION OF MT
C>RCS
C>RCS REVISION 1.3.1.1  89/05/20  16:42:53  BUCKLEY
C>RCS TEMP. TEST OF MT BEFORE SUBMITTING
C>RCS
C>RCS REVISION 1.3  89/05/18  12:20:45  BUCKLEY
C>RCS FINAL TEST OF MT BEFORE SUBMITTING
C>RCS
C>RCS REVISION 1.2  89/05/15  14:48:04  BUCKLEY
C>RCS INITIAL INSTALLATION OF MT INTO RCS FORM.
C>RCS
C>RCS REVISION 1.1  89/01/17  16:50:11  BUCKLEY
C>RCS INITIAL REVISION
C>RCS
C
C## D E S C R I P T I O N:
C
C     THIS ROUTINE ACTS AS A LINKAGE BETWEEN THE MAIN DRIVER ROUTINE
C     ZZTP AND THE MINIMIZATION ROUTINE TO BE TESTED.
C
C     THE ACTUAL CODE REQUIRED IS QUITE SIMPLE, BUT IT MUST BE PROVIDED
C     BY THE PERSON WHO HAS A ROUTINE TO BE TESTED.  THE STANDARD
C     DECLARATIONS MARKED BELOW SHOULD NOT BE CHANGED, BUT OTHERS WHICH
C     ARE PARTICULAR TO THE ALGORITHM BEING TESTED MAY BE APPROPRIATE
C     AND MAY BE ADDED.
C
C     THIS PARTICULAR SAMPLE IS SET UP TO TEST ONE OF TWO DIFFERENT
C     ALGORITHMS; THESE ARE DOCUMENTED BELOW.  EITHER OF THESE ILLUS-
C     TRATES WHAT IS REQUIRED TO IMPLEMENT AN ALGORITHM FOR TESTING.
C
C     ESSENTIALLY, THIS ROUTINE IS USED TO DO TWO THINGS:
C
C        1. IT DECLARES AND INITIALIZES ANY VALUES WHICH
C           WILL BE NEEDED BY THE TESTER'S ALGORITHM.
C
C        2. IT BREAKS UP THE WORKING STORAGE PROVIDED IN  RWORK AND
C           IWORK, AND SETS UP THE APPROPRIATE CALL TO THE SUBROUTINE.
C
C     IN ADDITION, THIS ROUTINE TAKES CARE OF A COUPLE OF MINOR TASKS.
C     IT MAY CHECK THAT SUFFICIENT STORAGE HAS BEEN PROVIDED.
C     THIS IS ILLUSTRATED IN THE SAMPLE PROVIDED BELOW.
C
C     MOST ROUTINES WILL REQUIRE VERY LITTLE CODE HERE AND THE CODE
C     WILL BE SHORT AND SIMPLE.
C
C-----THE VARIABLES IN THE CALLING SEQUENCE HAVE THE FOLLOWING MEANINGS.
C
C     N       THE DIMENSION OF THE PROBLEM.
C     X       NORMALLY, THE INITIAL GUESS AT THE SOLUTION (ON ENTRY) AND
C               THE SOLUTION FOUND (ON EXIT). ALSO SEE THE COMMENTS
C               BELOW ON REVERSE COMMUNICATION.
C     F       THE FUNCTION VALUE AT THE SOLUTION (ON EXIT, ALSO SEE RC).
C     G       THE GRADIENT VALUE AT THE SOLUTION (ON EXIT, ALSO SEE RC).
C     ACC     THE DESIRED ACCURACY.
C
C     STATUS  IT IS NORMALLY 0 ON ENTRY AND IT CONTAINS AN ERROR CODE ON
C             EXIT. ALSO SEE THE NOTE BELOW ON REVERSE COMMUNICATION.
C             THE VALUES FOR STATUS REFER TO THE VALUES WHICH WILL BE
C             SET BY ZZTP AND PASSED INTO ZZLINK, AND TO THE VALUES
C             WHICH ZZTP EXPECTS TO BE RETURNED BY ZZLINK. IT MAY
C             BE NECESSARY TO USE ALTERED VALUES TO COMMUNICATE WITH
C             THE TEST ALGORITHM. SEE FOR EXAMPLE THE VALUE OF STATUS
C             RETURNED FROM BBLNIR. THE STATUS CODES ARE DEFINED BELOW.
C
C     SUBR    THE NUMBER IDENTIFYING THE SUBROUTINE TO BE USED.
C
C     TRACES, NTR
C             THESE ARE FOR PROVIDING TRACES. NTR IS THE NUMBER OF TRACE
C             FLAGS CONTAINED IN THE ARRAY TRACES. SEE ZZTP FOR MORE
C             INFORMATION.
C
C     PRINTL  SEE ZZPRNT - SET BY "PRINT = _"
C     MAX     SEE ZZEVAL - MAX NUMBER OF FUNCTION EVALUATIONS, SET
C                          BY "MAX = _" IN PROBLEM DEFINITION.
C     DERVMD  SEE ZZEVAL - SET BY "DERIVATIVE = _"
C
C     USER(NU) SEE ZZTP AND THE DESCRIPTION OF THE USER NAME FEATURE.
C             THIS IS THE ARRAY OF GENERIC USER PARAMETERS.
C
C     RWORK)  A REAL (OR DOUBLE PRECISION) WORK ARRAY. LR IS
C     LR   )  THE DIMENSION OF RWORK. THIS MEANS THAT LWORK ELEMENTS OF
C     EXTRA)  RWORK MAY BE TAKEN FOR THE USER'S TEST ROUTINE. THE REAL
C             LENGTH OF RWORK IS  LWORK + EXTRA, AND IT IS ASSUMED IN
C             ZZEVAL WHERE THE TEST FUNCTIONS ARE EVALUATED THAT THE
C             THE REMAINING EXTRA LOCATIONS CAN BE USED IF NECESSARY IN
C             EVALUATING THE FUNCTIONS.
C
C     IWORK   AN INTEGER WORK ARRAY.
C     LI      THE DIMENSION OF IWORK.
C
C-----STATUS CODES.
C
C            STATUS IS 'NORMAL' ON ENTRY UNLESS REVERSE COMMUNICATION
C     IS BEING USED. STATUS SHOULD BE 'NORMAL' ON RETURN UNLESS THERE
C     IS AN ERROR OR REVERSE COMMUNICATION IS USED.
C
C     ON ENTRY, THE FOLLOWING VALUES ARE USED:
C
C     NORMAL        START NORMAL RUN
C     NORMFG        JUST LIKE NORMAL, BUT F AND G COMPUTED AT X ALREADY.
C     RCSTRT        START REVERSE COMMUNICATION
C     RCRPT         REVERSE COMMUNICATION REENTRY
C
C     ON RETURN TO ZZTP, THE FOLLOWING VALUES ARE EXPECTED BY ZZTP.
C
C     NORMAL        COMPLETION OF A TEST WITH NO ERRORS.
C     RCRPT         REPEAT REVERSE COMMUNICATION
C     RCF,RCFG,RCG  RESERVED FOR REVERSE COMMUNICATION.
C     XSFUNC        TERMINATION OCCURRED BECAUSE TOO MANY FUNCTION
C                         EVALUATIONS WERE DONE.
C     NOSTOR        EXECUTION OF THE TEST NEVER BEGAN BECAUSE TOO LITTLE
C                         STORAGE WAS AVAILABLE.
C     IPMIN         TERMINATION OCCURRED BECAUSE THE INITIAL POINT WAS
C                         FOUND TO BE A MINIMUM.
C     NOF0          THE FUNCTION OR GRADIENT WAS UNDEFINED AT X0.
C     RABORT        AN ABORT WAS REQUESTED BY THE EVALUATION ROUTINE.
C
C-----REVERSE COMMUNICATION.
C
C             THE MEANING AND USE OF THIS FEATURE IS DESCRIBED IN THE
C     EXTERNAL DOCUMENTATION. HERE WE JUST DESCRIBE ITS USE WITH RESPECT
C     TO THE VALUE OF STATUS.
C
C             NORMALLY, STATUS IS SET TO 0 BEFORE CALLING ZZLINK. IF
C     FUNCTION VALUES ARE TO BE OBTAINED THROUGH REVERSE COMMUNICATION,
C     THEN STATUS SHOULD BE SET TO 1 BEFORE CALLING ZZLINK. IN THIS
C     CASE, BOTH THE FUNCTION VALUE F AND (IF REQUIRED) THE GRADIENT
C     VALUE G AT THE INITIAL POINT X MUST BE DEFINED ON ENTRY. THESE
C     VALUES CAN THEN BE PASSED BY ZZLINK TO THE TEST ALGORITHM.
C
C             WHEN THE TEST ALGORITHM REQUIRES A FUNCTION AND/OR
C     GRADIENT VALUE, IT SHOULD RETURN TO ZZLINK. THE POINT AT WHICH
C     THE FUNCTION AND/OR GRADIENT VALUE IS TO BE EVALUATED MUST BE IN
C     X. AS ILLUSTRATED WITH BBLNIR, A RETURN TO ZZLINK SHOULD BE MADE
C     WITH X DEFINED AS JUST STATED. BEFORE RETURNING FROM ZZLINK TO
C     ZZTP, STATUS SHOULD BE SET AS FOLLOWS :
C
C             STATUS = RCF     IF ONLY F(X) IS REQUIRED.
C                    = RCFG    IF BOTH F(X) AND G(X) ARE REQUIRED.
C                    = RCG     IF ONLY G(X) IS REQUIRED.
C
C     IT IS UP TO ZZTP TO DETERMINE F(X) AND G(X). ZZTP WILL THEN CALL
C     ZZLINK WITH STATUS = RCRPT. WHEN STATUS = RCRPT ON ENTRY TO
C     ZZLINK, CONTROL SHOULD BE PASSED DIRECTLY TO THE TEST ALGORITHM,
C     AS IS DONE FOR BBLNIR. NOTE THAT I1,I2,... ARE DECLARED AS SAVE
C     SO THAT THEIR VALUES ARE RETAINED.
C
C             NOTE THAT WE HAVE MADE NO STATEMENTS ABOUT ARGUMENTS
C     PASSED TO OR FROM THE TEST ALGORITHM. HERE WE ARE ONLY CONCERNED
C     WITH THE VALUES PASSED BETWEEN ZZTP AND ZZLINK.
C
C## E N T R Y   P O I N T S: THE NATURAL ENTRY ZZLINK
C## S U B R O U T I N E S:
C
C     THE ONLY CALLS ARE TO THE MINIMIZATION ROUTINES
C
C        BBLNIR   CONMIN
C
C     (WHICH INCLUDES CALLS TO THE ENTRY BBLSET IN BBLNIR)
C
C     NINT...  INTRINSIC TO GET NEAREST INTEGER
C
C## P A R A M E T E R S:

      DOUBLE PRECISION  ZERO,       ONE,       TWO,       THREE
C!!!! REAL              ZERO,       ONE,       TWO,       THREE
      PARAMETER (       ZERO = 0D0, ONE = 1D0, TWO = 2D0, THREE = 3D0)

      DOUBLE PRECISION  FOUR,       FIVE,      SIX,       SEVEN
C!!!! REAL              FOUR,       FIVE,      SIX,       SEVEN
      PARAMETER (       FOUR = 4D0, FIVE = 5D0,SIX = 6D0, SEVEN = 7D0)

      DOUBLE PRECISION  EIGHT,         NINE,          TEN
C!!!! REAL              EIGHT,         NINE,          TEN
      PARAMETER (       EIGHT = 8D0,   NINE = 9D0,    TEN = 10D0     )
C-----                DEFINITION OF STATUS CODES.

C--ON ENTRY:
      INTEGER     CNORML,     CRCSTR,     CRCRPT,     CRCNFG
      PARAMETER ( CNORML = 0, CRCSTR = 1, CRCRPT = 2, CRCNFG = 3 )

      INTEGER     CNRMFG,      CPSTHR
      PARAMETER ( CNRMFG = -1, CPSTHR = -2 )

C--ON EXIT:
      INTEGER     CDONE,     CRCF,     CRCFG,     CRCG
      PARAMETER ( CDONE = 0, CRCF = 1, CRCFG = 2, CRCG = 3 )

      INTEGER     CNSTOR,      CIPMIN,      CIPUNF,      CBDMTH
      PARAMETER ( CNSTOR = -1, CIPMIN = -2, CIPUNF = -3, CBDMTH = -4 )

      INTEGER     CLSFAL,      CNODSC,      CXSFNC,      CPSBCK
      PARAMETER ( CLSFAL = -5, CNODSC = -6, CXSFNC = -7, CPSBCK = -8  )

      INTEGER     CRABRT,      CUSERV
      PARAMETER ( CRABRT = -9, CUSERV = -10 )

C## L O C A L   D E C L:

C-----STANDARD DECLARATIONS.

      EXTERNAL  ZZFNS, ZZINNR

      INTEGER      I1,  I2,  I3,  I4,  I5

      DOUBLE PRECISION  ZZINNR

C-----STATUS CODES

      INTEGER SNRMFG, SNORML, SRCSTR, SRCRPT, SRCNFG, SPSTHR
      INTEGER NORMFG, NORMAL, RCSTRT,  RCRPT, RCNOFG, PSTHRU

      INTEGER SDONE,  SRCF,   SRCFG,    SRCG, SNSTOR, SIPMIN, SPSBCK
      INTEGER DONE,   RCF,    RCFG,     RCG,  NOSTOR, IPMIN,  PSBACK
      INTEGER SIPUNF, SBDMTH, SLSFAL, SNODSC, SXSFNC, SRABRT, SUSERV
      INTEGER IPUNDF, BDMETH, LSFAIL, NODESC, XSFUNC, RABORT, USERV

C-----DECLARATIONS FOR BBLNIR.

      DOUBLE PRECISION  DECRF
C!!!! REAL              DECRF

C## S A V E:

      SAVE  I1, I2, I3, I4, I5
      SAVE  NORMFG, NORMAL, RCSTRT, RCRPT, RCNOFG, PSTHRU
      SAVE  DONE,   RCF,    RCFG,   RCG,  NOSTOR,  IPMIN, IPUNDF, BDMETH
      SAVE  LSFAIL, NODESC, RABORT, XSFUNC, USERV, PSBACK

C## E Q U I V A L E N C E S: NONE ARE DEFINED.
C## C O M M O N:             NONE IS DEFINED.
C## D A T A:                 NONE ARE SET.

      DATA  NORMFG/CNRMFG/, NORMAL/CNORML/, RCSTRT/CRCSTR/,
     -      RCRPT /CRCRPT/, RCNOFG/CRCNFG/, PSTHRU/CPSTHR/

      DATA  DONE  /CDONE/,  RCF   /CRCF/,   RCFG  /CRCFG/,  RCG/CRCG/
     -      NOSTOR/CNSTOR/, IPMIN /CIPMIN/, IPUNDF/CIPUNF/,
     -      BDMETH/CBDMTH/, LSFAIL/CLSFAL/, NODESC/CNODSC/,
     -      RABORT/CRABRT/, XSFUNC/CXSFNC/, USERV /CUSERV/,
     -      PSBACK/CPSBCK/
C##                                                E X E C U T I O N
C##                                                E X E C U T I O N
      GOTO ( 1000, 2000 )  SUBR

C-----         CALL BUCKLEY-LENIR ROUTINE BBLNIR.
 1000 CONTINUE

      IF ( STATUS .NE. RCRPT .AND. STATUS .NE. RCNOFG ) THEN

C        INITIALIZE, EXCEPT FOR A REVERSE COMMUNICATION RE-ENTRANCE.
         I1 = 1
         I2 = I1 + N
         I3 = I2 + N
         I4 = I3 + N
         I5 = LR - 3*N
         IF ( I5 .LT. 0 ) THEN
            STATUS = NOSTOR
            GOTO 90000
         ENDIF

         CALL BBLSET ( NINT(USER(1)),  NINT(USER(2)), NINT(USER(3)),
     -                 NINT(USER(13)), NINT(USER(4)), NINT(USER(5)),
     -                 NINT(USER(14)),
     -                 USER(6), USER(7),
     -                 USER(8)  .NE. ZERO , USER(9)  .NE. ZERO ,
     -                 USER(10) .NE. ZERO , USER(11) .NE. ZERO ,
     -                 USER(12) .NE. ZERO ,
     -                 RELF, RELG,
     -                 TRACUN,  TRACES                            )
      ENDIF

      CALL BBLNIR ( ZZFNS, N, X, F, USER(15), G, ACC, STATUS, ZZINNR,
     -    RWORK(I1), RWORK(I2), RWORK(I3), RWORK(I4),  I5, IW, RW, DW )

      IF ( STATUS .EQ. USERV ) THEN
         STATUS = USERV - NINT(RWORK(I2))
      ENDIF
      GOTO 90000

C-----         CALL SHANNO'S CG ALGORITHM CONMIN.
 2000 CONTINUE

      CALL CONMIN(ZZFNS,N,X,F,G,ACC,STATUS,RWORK,LR,NINT(USER(19)),
     -                  TRACES, TRACUN, NTR )

C       REVISE STATUS CODES TO MATCH THOSE FOR MINTEST
      IF      ( STATUS .EQ.  0 ) THEN
         STATUS = DONE
      ELSE IF ( STATUS .EQ.  1 ) THEN
         STATUS = XSFUNC
      ELSE IF ( STATUS .EQ.  2 ) THEN
         STATUS = LSFAIL
      ELSE IF ( STATUS .EQ.  3 ) THEN
         STATUS = NODESC
      ENDIF
      GOTO 90000

C## E X I T
90000      RETURN

C## F O R M A T S:  NONE ARE DEFINED.
C##                 E N D         OF ZZLINK.
                    END
