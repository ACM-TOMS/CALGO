      SUBROUTINE ZZEVAL ( ZZUFNC, N, X, F, G, INDIC, IW, RW, DW )

C## A R G U M E N T S:
                       INTEGER          INDIC,   N,   IW(*)

                       EXTERNAL         ZZUFNC
                       REAL             F,      X(N), G(N)
C!!!!                  DOUBLE PRECISION F,      X(N), G(N)

                       DOUBLE PRECISION DW(*)
                       REAL             RW(*)

C## S T A T U S:
C               IGNORE LINES BEGINNING WITH  "C!!!!" .
C
C               SYSTEM  DEPENDENCE:                      NONE.
C               SINGLE/DOUBLE CONVERSION: NEEDED (SEE CONVRT).
C
C               THIS VERSION IS IN   S I N G L E   PRECISION.
C!!!!           THIS VERSION IS IN   D O U B L E   PRECISION.
C
C>RCS $HEADER: EVAL.F,V 2.1 91/11/20 10:46:57 BUCKLEY EXP $
C>RCS $LOG:     EVAL.F,V $
C>RCS REVISION 2.1  91/11/20  10:46:57  BUCKLEY
C>RCS FINAL SUBMISSION TO TOMS
C>RCS
C>RCS REVISION 2.0  90/07/31  11:15:06  BUCKLEY
C>RCS ADDED REVISED BLAS.
C>RCS
C>RCS REVISION 1.9  89/06/30  13:11:45  BUCKLEY
C>RCS PREPARING SUBMITTED VERSION OF MT
C>RCS
C>RCS REVISION 1.3  89/05/18  12:47:25  BUCKLEY
C>RCS FINAL TEST OF MT BEFORE SUBMITTING
C>RCS
C>RCS REVISION 1.2  89/05/15  14:53:33  BUCKLEY
C>RCS INITIAL INSTALLATION OF MT INTO RCS FORM.
C>RCS
C>RCS REVISION 1.1  89/01/07  14:35:47  BUCKLEY
C>RCS INITIAL REVISION
C>RCS
C
C## D E S C R I P T I O N:
C
C     THIS SUBROUTINE EVALUATES A TEST FUNCTION AT THE GIVEN
C     POINT "X".  IT RETURNS THE VALUE OF THE FUNCTION AND / OR
C     THE VALUE OF THE GRADIENT AT X.  IT ALLOWS THE APPLICATION OF A
C     NONLINEAR SCALING TO THE FUNCTION IF DESIRED (SEE FSCALE BELOW).
C     IT ALSO ALLOWS THE USE OF FINITE DIFFERENCES (SEE SDRVMD BELOW).
C     IT CAN ALSO ACT AS A NOOP, I.E. AS A DO NOTHING ROUTINE; (SEE
C     INDIC BELOW).
C
C-----ON ENTRY:
C
C        ZZUFNC THE NAME OF THE FUNCTION TO EVALUATE.  THERE
C               MUST BE A SUBROUTINE PROVIDED OF THE FORM
C
C                    SUBROUTINE ZZUFNC(INDIC,N,X,F,G,IW,RW,DW)
C
C                (WHERE N, X, F, G, INDIC, IW, RW AND DW HAVE THE
C                 SAME MEANING AS IN THIS SUBROUTINE ZZEVAL.)
C
C        N     THE DIMENSION OF THE PROBLEM, I.E. THE
C              NUMBER OF VARIABLES IN THE FUNCTION ZZUFNC.
C
C        X     CONTAINS THE VALUE OF THE N-COORDINATES X[1],...,X[N]
C              AT WHICH TO EVALUATE THE FUNCTION.
C
C        INDIC = DOF   ONLY EVALUATE THE FUNCTION.
C              = DOFG  EVALUATE BOTH.
C              = DOG   ONLY EVALUATE THE GRADIENT.
C              = NONE  ACTUALLY, IF INDIC HAS ANY VALUE OTHER THAN
C                      ONE OF THE FIRST THREE, THEN  JUST CALL ZZUFNC
C              WITH THIS SAME CODE FOR INDIC; I.E. ZZEVAL SHOULD DO
C              NOTHING. THIS IS INTENDED FOR THE CONVENIENCE OF THE
C              WRITER OF ZZUFNC.
C
C              NOTE THAT THE VALUES OF THESE CODES CAN BE REDEFINED
C              THROUGH THE ENTRY POINT ZZFDEF BELOW. DEFAULT VALUES
C              ARE GIVEN IN THE PARAMETER SECTION BELOW.
C
C        IW    THESE ARE 3 WORK ARRAYS WHICH ARE NOT USED AT ALL BY
C        RW    ZZEVAL, BUT WHICH ARE JUST PASSED TO THE USER'S
C        DW    ROUTINE ZZUFNC TO BE USED AS DESIRED. WITH THESE ARRAYS
C              AVAILABLE, IT IS OFTEN NOT NECESSARY TO USE REVERSE
C              COMMUNICATION. NOTE THAT THERE IS ONE AVAILABLE OF
C              EACH BASIC NUMERIC TYPE.
C
C-----ON EXIT:
C
C        F           CONTAINS THE FUNCTION VALUE (WITH THE SCALING
C                    APPLIED IF REQUIRED).
C
C        G           CONTAINS THE GRADIENT VALUE (WITH THE SCALING
C                    APPLIED IF REQUIRED).
C
C     NEITHER F NOR G IS REFERENCED UNLESS ITS VALUE IS REQUESTED.
C
C    INDIC = OK    THE REQUEST MADE ON THE CALL WAS COMPLETED SATIS-
C                  FACTORILY. F AND/OR G ARE AVAILABLE AS REQUESTED.
C            ABORT THE MINIMIZATION ROUTINE WHICH CALLED ZZEVAL IS
C                  HEREBY REQUESTED TO EXIT IMMEDIATELY TO THE ROUTINE
C                  WHICH CALLED IT. THIS CAN BE USED BY THE ROUTINE
C                  ZZUFNC TO TRIGGER PREMATURE TERMINATION DUE TO
C                  CIRCUMSTANCES OF WHICH THE MINIMIZATION ROUTINE MAY
C                  NOT BE AWARE.
C            LIMIT TERMINATE THE MINIMIZATION; THE PRESET LIMIT ON THE
C                  NUMBER OF FUNCTION EVALUATIONS ALLOWED HAS BEEN
C                  EXCEEDED. SEE MAXFN BELOW.
C            NOF   THE FUNCTION VALUE COULD NOT BE DETERMINED.
C            NOG   THE GRADIENT VALUE COULD NOT BE DETERMINED.
C            NOFG  NEITHER F NOR G COULD BE EVALUATED.
C
C            THESE CODES CAN BE REDEFINED THROUGH AN ENTRY POINT BELOW,
C            AND HAVE DEFAULT VALUES SPECIFIED IN THE PARAMETER SECTION.
C
C-----SET THROUGH ENTRY POINT CALLS.
C
C  ZZESRT ( FSCALE, SDRVMD, MAXFN )      M A N D A T O R Y
C
C           THIS IS CALLED BEFORE MINIMIZING EACH FUNCTION. THIS CALL
C                 IS   M A N D A T O R Y.
C
C        FSCALE      CONTROLS THE NONLINEAR SCALING OF ZZUFNC.
C
C              = 0   NO EFFECT.
C
C              = K>0 THIS ROUTINE COMPUTES AND RETURNS  FF( ZZUFNC(X) ),
C                    WHERE FF IS THE  K-TH  OF THE NONLINEAR FUNCTIONS
C                    OF ONE VARIABLE DEFINED IN THE ROUTINE ZZFSCL.
C
C                NOTE THAT FOR CERTAIN SCALINGS, IF YOU CALL ZZEVAL
C                JUST FOR A GRADIENT VALUE, IT MAY BE NECESSARY TO
C                REQUEST A FUNCTION VALUE AS WELL IN ORDER TO DO THE
C                SCALING.  THAT FUNCTION VALUE WILL NOT BE PASSED BACK.
C                THOSE WHICH DO NOT REQUIRE F FOR THE SCALING ARE THOSE
C                WITH K = 1,2,..,REQF - 1. FOR K = REQF,..., THE VALUE
C                OF F IS NECESSARY.
C
C        SDRVMD  THIS SPECIFIES THE METHOD BY WHICH DERIVATIVES ARE
C                TO BE COMPUTED, WHEN REQUESTED. THE CHOICE IS BETWEEN
C
C            CANAL   USE ANALYTIC FORMULAE WHICH MUST BE CODED AND
C                 AVAILABLE IN THE USER ROUTINE  ZZUFNC.
C
C            CDIFF   USE FINITE DIFFERENCE APPROXIMATIONS. IN THIS CASE,
C                 THE USER ROUTINE MAY IGNORE CALLS WITH INDIC <> JUSTF,
C                 AND NEED ONLY BE ABLE TO COMPUTE FUNCTION VALUES.
C                 FURTHER COMMENTS APPEAR IN THE DISCUSSION OF FINITE
C                 DIFFERENCE COMPUTATIONS (BELOW).
C
C            CTEST   IN THIS CASE BOTH ANALYTIC AND FINITE DIFFERENCES
C                 ARE COMPUTED.  THEY ARE THEN COMPARED AND A RECORD
C                 IS KEPT TO SEE TO WHAT EXTENT THEY DISAGREE. A
C                 RECORD OF THE LEVEL OF AGREEMENT IS AVAILABLE
C                 THROUGH THE ENTRY POINT ZZECHK GIVEN BELOW. A MORE
C                 COMPLETE DESCRIPTION IS ALSO GIVEN WHERE ZZECHK IS
C                 DISCUSSED BELOW.
C
C           CFIRST   THIS CASE IS PRECISELY THE SAME AS FOR CTEST, WITH
C                 THE SOLE EXCEPTION THAT THE TESTING ONLY TAKES PLACE
C                 ON THE FIRST CALL TO ZZEVAL.
C
C                 THE INTEGER VALUES OF THE CODES FOR CANAL, ETC ARE
C                 SET IN THE PARAMETER SECTION BELOW. THEY MAY BE
C                 RESET VIA THE ENTRY POINT ZZEDEF DESCRIBED BELOW.
C
C        MAXFN       THE MAXIMUM VALUE ALLOWED FOR THE COUNT IFNCT.
C
C             <= 0   ON ENTRY SPECIFIES NO MAXIMUM, I.E. MAXFN IS
C                    IGNORED.
C
C           = K>0    SPECIFIES THE MAXIMUM NUMBER OF TIMES THAT ZZUFNC
C                    MAY BE CALLED.  IF THE FUNCTION EVALUATION COUNT
C                    IN  IFNCT  IS GREATER THAN OR EQUAL TO  MAXFN ON
C                    ENTRY TO ZZEVAL, THEN THE FUNCTION IS NOT
C                    EVALUATED AND THE RETURN CODE INDIC IS SET AS
C                    ABOVE. NOTE THAT THE COUNT IN IFNCT DOES   N O T
C                    INCLUDE FUNCTION EVALUATIONS USED FOR COMPUTING
C                    FINITE DIFFERENCE GRADIENTS.
C
C
C      THE NEXT FOUR PARAMETERS ARE NOT IN THE CALLING SEQUENCE OF
C        ZZESRT, BUT THEY ARE INITIALIZED WHEN ZZESRT IS CALLED.
C
C        IFNCT       COUNTS THE NUMBER OF TIMES THE ROUTINE IS CALLED
C                    TO EVALUATE THE FUNCTION.  IT IS INITIALIZED TO 0
C                    DURING THE CALL TO  ZZESRT.
C
C        IGRCT       COUNTS THE NUMBER OF TIMES THE ROUTINE IS CALLED
C                    TO EVALUATE THE GRADIENT.  IT IS INITIALIZED TO 0
C                    DURING THE CALL TO ZZESRT.
C
C        FTIME       RECORDS THE TIME ACCUMULATED IN EVALUATING THE
C                    FUNCTION AND/OR THE GRADIENT.  IT IS PRESET TO ZERO
C                    WHEN ZZESRT IS CALLED.  THE TIME USED IN THE FINAL
C                    SCALING IS INCLUDED IN THE TIMING WHEN THE VALUE OF
C                    FSCALE IS NON-ZERO. TIMING COMMENCES ON ENTRY TO
C                    ZZEVAL, AND ENDS JUST BEFORE RETURN FROM ZZEVAL.
C
C        ERR         THE ESTIMATE OF THE ERROR BETWEEN THE ANALYTIC
C                    AND DIFFERENCE VALUES FOR THE GRADIENT IS RECORDED
C                    IN A SET OF VARIABLES ERR, SERR, DCNT, INDEX AND
C                    GCNT, SO THESE ARE INITIALIZED TO 0.
C
C  ZZESET ( TRF, TRG, TRTEST, ITRUN )
C
C           THIS IS CALLED BEFORE USING ZZEVAL (THESE VALUES ALSO HAVE
C           INTERNALLY SET DEFAULT VALUES GIVEN IN [..], SO THE CALL TO
C           ZZESET IS NOT MANDATORY.)
C
C        TRF  = TRUE IF THE FUNCTION VALUE IS TO BE PRINTED [FALSE]
C        TRG  = TRUE IF THE GRADIENT VALUE IS TO BE PRINTED [FALSE]
C        TRTEST  = TRUE IF THE FUNCTION VALUES ARE TO BE PRINTED
C                  DURING TESTING MODE [FALSE]
C
C        ITRUN  THE UNIT NUMBER FOR OUTPUT OF COMPUTED VALUES [6]
C
C     NOTE THAT AN ERROR MESSAGE IS PRINTED WHEN THE MAXIMUM NUMBER
C     OF FUNCTION EVALUATIONS IS EXCEEDED, PROVIDED EITHER TRF OR
C     TRG IS TRUE.
C
C  ZZEDDF ( SANAL, SDIFF, STEST, SFIRST )
C
C           THIS MAY BE CALLED BEFORE USING ZZEVAL, AS FOR ZZESET. THIS
C           ALLOWS THE CODES FOR ANAL, ETC., TO BE REDEFINED. ALL HAVE
C           DEFAULTS, SO THIS CALL IS NOT MANDATORY.
C
C        SANAL   THE INTEGER VALUE FOR THE CODE FOR USING ANALYTIC
C                   DERIVATIVES [CANAL].
C        SDIFF   THE INTEGER VALUE FOR THE CODE FOR USING FINITE
C                   DIFFERENCES TO APPROXIMATE DERIVATIVES [CDIFF].
C        STEST   THE INTEGER VALUE FOR THE CODE FOR USING BOTH CANAL
C                   AND CDIFF AND DOING A TEST FOR AGREEMENT [CTEST].
C        SFIRST  THE INTEGER VALUE FOR THE CODE FOR USING BOTH CANAL
C                   AND CDIFF ON THE FIRST ITERATION ONLY.
C
C  ZZEFDF ( SDOF, SDOG, SDOFG, SNONE )
C
C     THIS MAY BE CALLED BEFORE USING ZZEVAL, JUST AS FOR ZZEDEF.
C
C        DOF   THE CODE INDICATING THAT JUST THE FUNCTION VALUE IS
C              DESIRED. [JUSTF]
C        DOG   THE CODE INDICATING THAT JUST THE GRADIENT VALUE IS
C              DESIRED. [JUSTG]
C        DOFG  THE CODE INDICATING THAT BOTH THE FUNCTION AND GRADIENT
C              VALUES ARE DESIRED. [BOTH]
C        SNONE THE CODE INDICATING THAT NO ACTION IS TO BE TAKEN AND
C              THAT ZZUFNC SHOULD BE CALLED WITH NO OTHER PROCESSING.
C
C  ZZERDF ( OK, LIMIT, ABORT, NOF, NOG, NOFG )
C
C     THIS MAY BE CALLED BEFORE USING ZZEVAL, JUST AS FOR ZZEDEF.
C
C     OK    THIS CODE INDICATES THAT THE REQUEST WAS SUCCESSFULLY DONE.
C     ABORT THIS MEANS THAT THE CALLING ROUTINE SHOULD IMMEDIATELY
C           TERMINATE THE MINIMIZATION AND RETURN TO THE ROUTINE WHICH
C           CALLED IT.
C     LIMIT THIS MEANS THAT THE ALLOWED NUMBER OF FUNCTION EVALUATIONS
C           HAS BEEN EXCEEDED.
C     NOF   THIS MEANS THAT ZZEVAL WAS UNABLE TO SUCCESSFULLY EVALUATE
C           THE FUNCTION.
C     NOG   THIS MEANS THAT ZZEVAL WAS UNABLE TO SUCCESSFULLY EVALUATE
C           THE GRADIENT.
C     NOFG  THIS MEANS THAT ZZEVAL WAS UNABLE TO OBTAIN EITHER A
C           FUNCTION OR GRADIENT VALUE.
C
C-----AVAILABLE THROUGH ENTRY POINT CALLS AFTER A FUNCTION HAS BEEN
C     MINIMIZED:
C
C  ZZEGET ( FNCT, GRCT, TIME )
C
C      THIS MAY BE CALLED AFTER MINIMIZING A FUNCTION TO OBTAIN SOME
C      SIMPLE STATISTICS WHICH HAVE BEEN ACCUMULATED SINCE THE LAST
C      CALL TO  ZZESRT. THESE ARE
C
C         FNCT  THE NUMBER OF CALLS TO EVALUATE THE FUNCTION, I.E.
C                 CALLS WITH INDIC = JUSTF OR BOTH.
C
C         GRCT  THE NUMBER OF CALLS TO EVALUATE THE GRADIENT, I.E.
C                 CALLS WITH INDIC = JUSTG OR BOTH.
C
C         TIME  THE AMOUNT OF CPU TIME SPENT IN  ZZEVAL.
C
C  ZZECHK ( ERR, AVERR, INDX, ITERAT )
C
C
C         THIS MAY ALSO BE CALLED AFTER A SEQUENCE OF CALLS TO ZZEVAL.
C     IT GIVES AN ESTIMATE OF THE AGREEMENT BETWEEN ANALYTIC AND
C     FINITE DIFFERENCE DERIVATIVES. OF COURSE THESE VALUES ARE ONLY
C     DEFINED IF SDRVMD = CTEST.
C
C         ERR WILL BE RETURNED AS AN ESTIMATE OF THE LARGEST ERROR
C     WHICH OCCURRED AND AVERR IS AN ESTIMATE OF THE AVERAGE NUMBER OF
C     DECIMAL DIGITS OF AGREEMENT BETWEEN THE COMPONENTS OF THE ANALYTIC
C     AND DIFFERENCE DERIVATIVES.
C
C         TO BE SPECIFIC, WHEN IN TEST MODE, EACH COMPONENT OF THE
C     ANALYTIC DERIVATIVE IS COMPUTED, AND THAT IS RETURNED IN G AS THE
C     GRADIENT.  AS WELL, FOR EACH COMPONENT, A FINITE DIFFERENCE
C     APPROXIMATION IS COMPUTED (AS DESCRIBED BELOW) AND THE RELATIVE
C     DIFFERENCE BETWEEN THAT AND THE ANALYTIC COMPONENT IS DETERMINED.
C     THIS QUANTITY IS MONITORED, AND THE LARGEST SUCH VALUE IS
C     RECORDED. IN ADDITION, INDX INDICATES IN WHICH COMPONENT OF THAT
C     GRADIENT THE ERROR OCCURRED AND ITERAT TELLS WHICH GRADIENT
C     EVALUATION WAS IN PROGRESS WHEN THE ERROR OCCURRED; I.E. ITERAT
C     JUST RECORDS THE CURRENT VALUE OF IGRCT. NOTE THAT INDX
C     AND ITERAT ONLY REFER TO THE POINT AT WHICH THE LARGEST
C     ERROR OCCURRED.
C
C         IF THE FUNCTION AND GRADIENT EVALUATIONS ARE CORRECT, ONE
C     WOULD NORMALLY EXPECT THE RELATIVE ERROR TO BE OF THE ORDER OF
C     10**-(T/2), WHERE  T  IS THE NUMBER OF FIGURES OF RELATIVE
C     ACCURACY OF THE MACHINE IN USE.  HOWEVER, AS THE MINIMUM IS
C     APPROACHED AND THE GRADIENT COMPONENTS GENERALLY BECOME VERY
C     SMALL, THIS RELATIVE ACCURACY MAY BE MUCH WORSE THAN EXPECTED.
C     THEREFORE WE ALSO MAINTAIN AN ESTIMATE OF THE AVERAGE AGREEMENT.
C     HERE, FOR EACH COMPONENT OF EACH GRADIENT COMPUTATION, WE COMPUTE
C     THE BASE 10 LOG OF THE RELATIVE ACCURACY; THIS IS ROUGHLY THE
C     NUMBER OF SIGNIFICANT FIGURES OF AGREEMENT BETWEEN THE TWO VALUES.
C     THIS QUANTITY IS MONITORED AND  AVERR IS RETURNED AS THE AVERAGE
C     VALUE OF THE NUMBER OF SIGNIFICANT FIGURES OF AGREEMENT.
C
C          WHEN FUNCTION AND GRADIENT COMPUTATIONS ARE CORRECT, ERR WILL
C     GENERALLY BE AT LEAST AS SMALL AS  10**(-T/2), ALTHOUGH IT CAN BE
C     MORE LIKE 10**(-T/4).  GROSS BLUNDERS WILL USUALLY GIVE ERR A
C     VALUE VERY NEAR TO 1, BUT NOT ALWAYS.  IF ALL IS WELL, AVERR WILL
C     USUALLY BE ABOUT T/2; BLUNDERS WILL OFTEN RESULT IN AVERR BEING
C     NEAR 0 OR 1.
C
C-----FINITE DIFFERENCE COMPUTATIONS
C
C         FOR FIRST DERIVATIVES, SIMPLE FORWARD DIFFERENCES ARE USED.
C     TO ESTIMATE THE I-TH COMPONENT OF THE GRADIENT OF F, WE COMPUTE
C
C                 ( F(X + H*E[I]) - F(X) ) / H,
C
C     WHERE H = EPS * ABS(X[I]).  WHEN X[I] = 0, WE JUST CHOOSE H = EPS.
C     HERE EPS IS THE ROOT OF ETA, WHERE ETA DEFINES THE RELATIVE
C     MACHINE ACCURACY. THIS IS USED WHEN SDRVMD = CDIFF OR CTEST.
C
C         WHEN SDRVMD = CTEST, MORE INFORMATION IS REQUIRED; THUS WE
C     ALSO COMPUTE F(X + SQRT(H)*E[I]). THIS MEANS THAT WHEN IN TEST
C     MODE, TWICE AS MANY FUNCTION EVALUATIONS ARE NEEDED.  THIS IS
C     REQUIRED TO ELIMINATE SCALING EFFECTS IN THE ESTIMATE OF FIGURES
C     OF AGREEMENT.
C
C## E N T R Y   P O I N T S: THE NATURAL ENTRY ZZEVAL.
C
C     ZZESRT    TO INITIALIZE FOR TESTING EACH FUNCTION.
C     ZZESET    TO INITIALIZE CONTROL PARAMETERS.

C     ZZEDDF    TO REDEFINE CODES FOR DERIVATIVE CALCULATIONS.
C     ZZEFDF    TO REDEFINE CODES FOR FUNCTION EVALUATIONS.
C     ZZERDF    TO REDEFINE CODES FOR RETURN CODES.
C
C     ZZEGET    TO RETURN FUNCTION/GRADIENT COUNTS AND TIME.
C     ZZECHK    RETURNS THE ERROR VALUE IF IN TESTMODE.
C
C## S U B R O U T I N E S:
C     ZZSECS  ...FOR FUNCTION TIMING.
C     ZZMPAR  ...FOR MACHINE PARAMETERS.
C     ZZUFNC  ...THE USER ROUTINE.
C     ZZFSCL  ...PERFORMS SCALING.
C
C     SQRT,  MAX, ABS    ...INTRINSIC FUNCTIONS.
C     LOG10, MIN, SIGN   ...INTRINSIC FUNCTIONS.
C
C## P A R A M E T E R S:
      REAL              ZERO,       ONE,       TWO,       THREE
C!!!! DOUBLE PRECISION  ZERO,       ONE,       TWO,       THREE
      PARAMETER (       ZERO = 0D0, ONE = 1D0, TWO = 2D0, THREE = 3D0)

      REAL              FOUR,       FIVE,      SIX,       SEVEN
C!!!! DOUBLE PRECISION  FOUR,       FIVE,      SIX,       SEVEN
      PARAMETER (       FOUR = 4D0, FIVE = 5D0,SIX = 6D0, SEVEN = 7D0)

      REAL              EIGHT,         NINE,          TEN
C!!!! DOUBLE PRECISION  EIGHT,         NINE,          TEN
      PARAMETER (       EIGHT = 8D0,   NINE = 9D0,    TEN = 10D0     )
      REAL              R11,        R12,        R13,       R14
C!!!! DOUBLE PRECISION  R11,        R12,        R13,       R14
      PARAMETER (       R11 = 11D0, R12 = 12D0, R13 = 13D0,R14 = 14D0)

      REAL              R15,        R16,        R17,       R18
C!!!! DOUBLE PRECISION  R15,        R16,        R17,       R18
      PARAMETER (       R15 = 15D0, R16 = 16D0, R17 = 17D0,R18 = 18D0)

      REAL              R19,        R20,        R25,       R29
C!!!! DOUBLE PRECISION  R19,        R20,        R25,       R29
      PARAMETER (       R19 = 19D0, R20 = 20D0, R25 = 25D0,R29 = 29D0)

      REAL              R32,        R36,        R40,       R42
C!!!! DOUBLE PRECISION  R32,        R36,        R40,       R42
      PARAMETER (       R32 = 32D0, R36 = 36D0, R40 = 40D0,R42 = 42D0)

      REAL              R45,        R49
C!!!! DOUBLE PRECISION  R45,        R49
      PARAMETER (       R45 = 45D0, R49 = 49D0 )

      REAL              R50,        R56,        R84,       R90
C!!!! DOUBLE PRECISION  R50,        R56,        R84,       R90
      PARAMETER (       R50 = 50D0, R56 = 56D0, R84 = 84D0,R90 = 90D0)

      REAL              R100,            R180,           R200
C!!!! DOUBLE PRECISION  R100,            R180,           R200
      PARAMETER (       R100 = 100D0,    R180 = 180D0,   R200 = 200D0 )

      REAL              R256,            R360,           R400
C!!!! DOUBLE PRECISION  R256,            R360,           R400
      PARAMETER (       R256 = 256D0,    R360 = 360D0,   R400 = 400D0 )

      REAL              R600,            R681,           R991
C!!!! DOUBLE PRECISION  R600,            R681,           R991
      PARAMETER (       R600 = 600D0,    R681 = 681D0,   R991 = 991D0 )

      REAL              R1162,                 R2324
C!!!! DOUBLE PRECISION  R1162,                 R2324
      PARAMETER (       R1162 = 1162D0,        R2324 = 2324D0         )

      REAL              R10000,                R40000
C!!!! DOUBLE PRECISION  R10000,                R40000
      PARAMETER (       R10000 = 10000D0,      R40000 = 40000D0       )

      INTEGER     XEPS,     XSMALL,     XBIG
      PARAMETER ( XEPS = 1, XSMALL = 2, XBIG = 3 )

      INTEGER          REQF
      PARAMETER (      REQF = 2 )

C                DEFINE THE DERIVATIVE CODES

      INTEGER     CANAL,     CDIFF,     CTEST,     CFIRST
      PARAMETER ( CANAL = 1, CDIFF = 2, CTEST = 3, CFIRST = 4 )
C                DEFINE THE FUNCTION CODES

      INTEGER     JUSTF,     BOTH,     JUSTG,      NOOP
      PARAMETER ( JUSTF = 1, BOTH = 0, JUSTG = -1, NOOP = 2 )
C                DEFINE THE RETURN CODES
C   THE RETURN CODES TO BE USED BY THE FUNCTION EVALUATION ROUTINE
C   TO INDICATE TO THE MINIMIZATION ROUTINE WHETHER OR NOT THE CALL
C   WAS SUCCESSFUL.

      INTEGER     COK,     CABORT,      CLIMIT
      PARAMETER ( COK = 0, CABORT = -1, CLIMIT = -2 )

      INTEGER     CNOF,      CNOG,      CNOFG
      PARAMETER ( CNOF = -3, CNOG = -4, CNOFG = -5 )

C## L O C A L   D E C L:

      LOGICAL  TRF, TRG, FIRST, FONLY, GONLY, VALID, BAD, TRTEST, FCALL

      INTEGER         IFNCT,  FSCALE, IGRCT,  SDRVMD, REMBAD
      INTEGER         ITRUN,  MAXFN,  DERVMD, EXPENS, KK
      INTEGER   CASE, CALLS,  COUNT,  INDEX,  GCNT,  DCNT
      INTEGER  ZDCNT, ZINDEX, ZGCNT

      REAL              FT, FV, FTIME,  TT,  SCALE,  SERR, RH,ZERR,ZSERR
C!!!! DOUBLE PRECISION  FT, FV, FTIME,  TT,  SCALE,  SERR, RH,ZERR,ZSERR

      REAL              FVAL, FVAL2, ERR, ETA, EPS, H, ZZMPAR, TERR
C!!!! DOUBLE PRECISION  FVAL, FVAL2, ERR, ETA, EPS, H, ZZMPAR, TERR

C-----DECLARATIONS FOR ENTRY POINT DUMMY ARGUMENTS.

      INTEGER   DITRUN,  FSCAL,  MAXM,  FNCT, GRCT, SEXPEN
      INTEGER   DDERV, INDX, ITERAT

      INTEGER   ANAL,  DIFF,  TEST,  DOF,  DOG,  DOFG,  NONE, TFIRST
      INTEGER  SANAL, SDIFF, STEST, SDOF, SDOG, SDOFG, SNONE, SFIRST

      INTEGER     OK,  ABORT,  LIMIT,  NOF,  NOG,  NOFG
      INTEGER    SOK, SABORT, SLIMIT, SNOF, SNOG, SNOFG

      LOGICAL     DTRF,   DTRG, DTRTST

      REAL              TIME, ERROR, AVERR
C!!!! DOUBLE PRECISION  TIME, ERROR, AVERR

C## S A V E:

      SAVE   ITRUN, FSCALE, IFNCT, IGRCT, SERR, DCNT, EXPENS
      SAVE   TRF, TRG, FTIME, SDRVMD,  MAXFN, TRTEST
      SAVE   FIRST, ERR, INDEX, GCNT, EPS, ETA, FCALL
      SAVE   ANAL, DIFF, TEST, TFIRST, DOF, DOG, DOFG, NONE
      SAVE   OK, ABORT, LIMIT, NOF, NOG, NOFG

C## E Q U I V A L E N C E S: NONE ARE DEFINED.
C## C O M M O N:             NONE IS DEFINED.
C## D A T A:

      DATA  FIRST /.TRUE./, FCALL/.TRUE./

      DATA  TRF, TRG, TRTEST /  3 * .FALSE. /, ITRUN / 6 /

      DATA  SDRVMD/CANAL/, FSCALE/0/, MAXFN/0/

      DATA  ANAL/CANAL/,  DIFF/ CDIFF/,  TEST/ CTEST/, TFIRST/CFIRST/
      DATA   DOF/JUSTF/,   DOG/ JUSTG/,  DOFG/  BOTH/,   NONE/  NOOP/
      DATA    OK/  COK/, ABORT/CABORT/, LIMIT/CLIMIT/
      DATA   NOF/ CNOF/,   NOG/  CNOG/,  NOFG/ CNOFG/

C##                                                 E X E C U T I O N
C##                                                 E X E C U T I O N

C-----STATEMENT FUNCTION.

      BAD() =  CASE .EQ. ABORT .OR. CASE .EQ. LIMIT .OR. CASE .EQ. NOF
     -    .OR. CASE .EQ.  NOG  .OR. CASE .EQ.  NOFG

C-----FIRST TEST FOR NOOP CALL.

      IF ( TRF .OR. TRG ) WRITE ( ITRUN,99992 )INDIC,DOF,DOG,DOFG

      VALID = INDIC .EQ. DOF .OR. INDIC .EQ. DOG .OR. INDIC .EQ. DOFG

      IF ( .NOT. VALID ) THEN
         CALL ZZUFNC ( INDIC, N, X, F, G, IW, RW, DW )
         GOTO 90500
      ENDIF

      IF ( MAXFN .GT. 0 .AND. IFNCT .GE. MAXFN ) THEN
         GOTO 91000
      ENDIF

      DERVMD = SDRVMD

      IF ( FIRST ) THEN
         FIRST = .FALSE.
         ETA  = ZZMPAR(XEPS)
         EPS  = SQRT  (ETA)
      ENDIF

      IF ( FCALL ) THEN
         IF ( DERVMD .EQ. TFIRST ) DERVMD = TEST
         FCALL = .FALSE.
      ELSE
         IF ( DERVMD .EQ. TFIRST ) DERVMD = ANAL
      ENDIF

      FONLY = INDIC .EQ. DOF
      GONLY = INDIC .EQ. DOG

      CALL    ZZSECS (TT)
      FTIME = FTIME - TT

      IF ( .NOT. GONLY ) IFNCT  = IFNCT + 1
      IF ( .NOT. FONLY ) IGRCT  = IGRCT + 1

C-----FIRST COMPUTE REQUIRED FUNCTION AND/OR GRADIENT VALUES.

C     REPEAT IF REQUIRED TO SIMULATE EXPENSIVE CALL.

      ZERR   = ERR
      ZSERR  = SERR
      ZDCNT  = DCNT
      ZINDEX = INDEX
      ZGCNT  = GCNT
      REMBAD = OK

      DO 9900 KK = 1, EXPENS
          CASE = INDIC

C         NO OF EXTRA CALLS TO USER ROUTINE WHICH WILL BE NEEDED.

          IF ( DERVMD .EQ. ANAL  .OR.  FONLY ) THEN
             CALLS = 0
          ELSE
             CALLS = N
          ENDIF

C         FORCE FUNCTION EVALUATION IF REQUIRED FOR SCALING.

          IF ( FSCALE .GE. REQF .AND.  GONLY ) THEN
             CASE = DOFG
          ENDIF

C         FIRST COMPUTE  F(X) --- AND  G(X) IF NEEDED.

          CALL ZZUFNC ( CASE, N, X, FVAL, G, IW, RW, DW )

          IF ( BAD() ) THEN
             REMBAD = CASE
             GOTO 9900
          ENDIF

          IF ( INDIC .NE. DOG ) THEN
             FT = FVAL
          ENDIF

C         -----DO EXTRA CALLS, IF REQUIRED.
C              AFTER FIRST CALL, FUNCTION VALUES ONLY.

          DO 1500 COUNT = 1, CALLS

             TT = X(COUNT)

             IF ( TT .EQ. ZERO ) THEN
                H = EPS
             ELSE
                H = EPS * ABS( TT )
             ENDIF

             X(COUNT) = TT + H
C            COMPUTE  F( X + H * E[COUNT] )

             CASE = DOF
             CALL  ZZUFNC ( CASE, N, X, FVAL, G, IW, RW, DW )

             IF ( BAD() ) THEN
                REMBAD = CASE
                GOTO 9900
             ENDIF

             X(COUNT) = TT
             IF ( DERVMD .EQ. TEST ) THEN
C               ---IF TRACE REQUESTED, PRINT ESTIMATED AND ANALYTIC VALUES.

                IF ( TRTEST ) WRITE(ITRUN,99995)
     -                          G(COUNT),COUNT,(FVAL-FT)/H
C               ---ESTIMATE ERROR, AND LEAVE COMPUTED
C                  ANALYTIC GRADIENTS IN G.  USE F AT
C                  X + A * E[COUNT], FOR A = H AND SQRT(H).
                RH       = SQRT(H)
                X(COUNT) = TT + RH

                CASE = DOF
                CALL  ZZUFNC ( CASE, N, X, FVAL2, G, IW, RW, DW )
                IF ( BAD() ) THEN
                   REMBAD = CASE
                   GOTO 9900
                ENDIF

                X(COUNT) = TT
                IF ( ABS(FVAL2-FT) .GT. R100*ETA*ABS(FT) ) THEN
                   TERR = (FVAL-FT - H*G(COUNT))/
     -                    (FVAL2-FT - RH*G(COUNT))
                   IF (TT .GT. ONE) TERR = TERR / TT

C                  TRUNCATE TO INTERVAL [ETA,1].
                   TERR = MAX( MIN(ONE,ABS(TERR)), ETA )

C                  ESTIMATE NUMBER OF FIGURES OF AGREEMENT.
                   ZSERR = ZSERR - LOG10 (TERR)
                   ZDCNT = ZDCNT + 1

                   IF (TRTEST) WRITE(ITRUN,99994) TERR,-LOG10(TERR)
                   IF ( TERR .GT. ABS(ZERR) ) THEN
                      ZINDEX = COUNT
                      ZGCNT  = IGRCT
                      ZERR   = SIGN (TERR, ZERR)
                   ENDIF
                ELSE
C                  FLAG CASE WHERE THERE IS EXCESSIVE CANCELLATION.

                   ZERR = - ABS(ZERR)
                   IF (TRTEST) WRITE(ITRUN,99993)
                ENDIF
             ELSE
C               ---ESTIMATE GRADIENTS USING FORWARD FINITE DIFFERENCE
C                  FORMULAE AND STORE IN G.
                G(COUNT) = ( FVAL - FT ) / H
             ENDIF
 1500     CONTINUE

C         -----DO SCALING: DEFINE FV AND SCALE.  NOTE THAT IN SOME
C              INSTANCES THIS MAY REQUIRE AN EXTRA CALL TO GET THE
C              FUNCTION VALUE WHEN INDIC = DOG;  THIS WAS DONE IN
C              THE CALLS ABOVE.

          IF ( FSCALE .NE. 0 ) THEN
              CALL ZZFSCL( FT, FV, SCALE, FSCALE, FONLY, GONLY )
          ELSE
              FV    = FT
              SCALE = ONE
          ENDIF

C         -----NOW REVISE THE FUNCTION AND GRADIENT AS NECESSARY.

          IF ( .NOT. GONLY ) THEN
             F  = FV
          ENDIF

          IF ( .NOT. FONLY .AND. SCALE .NE. ONE ) THEN
             CALL ZZSCAL ( N, SCALE, G, 1 )
          ENDIF

 9900 CONTINUE
      ERR   = ZERR
      INDEX = ZINDEX
      GCNT  = ZGCNT
      DCNT  = ZDCNT
      SERR  = ZSERR
      INDIC = REMBAD
      GOTO 90000

C## E N T R Y  ZZESRT:
                      ENTRY  ZZESRT ( FSCAL, DDERV, MAXM, SEXPEN )
      FCALL  = .TRUE.
      FSCALE = FSCAL
      SDRVMD = DDERV
      MAXFN  = MAXM
      EXPENS = SEXPEN
      IFNCT  = 0
      IGRCT  = 0
      FTIME  = ZERO

      IF ( SDRVMD .EQ. TEST .OR. SDRVMD .EQ. TFIRST ) THEN
         ERR   = ZERO
         SERR  = ZERO
         DCNT  = 0
         INDEX = 0
         GCNT  = 0
      ENDIF
      RETURN

C## E N T R Y  ZZESET:
                     ENTRY  ZZESET ( DTRF, DTRG, DTRTST, DITRUN )
      TRF    = DTRF
      TRG    = DTRG
      TRTEST = DTRTST
      ITRUN  = DITRUN
      RETURN

C## E N T R Y  ZZEDDF:
                      ENTRY  ZZEDDF (  SANAL, SDIFF, STEST, SFIRST )
      ANAL   = SANAL
      DIFF   = SDIFF
      TEST   = STEST
      TFIRST = SFIRST
      RETURN

C## E N T R Y  ZZEFDF:
                      ENTRY  ZZEFDF ( SDOF, SDOG, SDOFG, SNONE )
      DOF  = SDOF
      DOG  = SDOG
      DOFG = SDOFG
      NONE = SNONE
      RETURN

C## E N T R Y  ZZERDF:
                      ENTRY ZZERDF(SOK,SABORT,SLIMIT,SNOF,SNOG,SNOFG)
      OK    = SOK
      ABORT = SABORT
      LIMIT = SLIMIT
      NOF   = SNOF
      NOG   = SNOG
      NOFG  = SNOFG
      RETURN

C## E N T R Y  ZZEGET:
                      ENTRY  ZZEGET ( FNCT, GRCT, TIME )
      FNCT = IFNCT
      GRCT = IGRCT
      TIME = FTIME
      RETURN

C## E N T R Y  ZZECHK:
                      ENTRY  ZZECHK ( ERROR, AVERR, INDX, ITERAT )
                      ERROR  = ERR
                      INDX   = INDEX
                      ITERAT = GCNT
                      AVERR = SERR / DCNT
                      RETURN
C## E X I T
90000      CALL ZZSECS (TT)

      FTIME = FTIME + TT
      IF ( TRF .AND. .NOT. BAD() ) WRITE (ITRUN,99998) F
      IF ( TRG .AND. .NOT. BAD() ) THEN
         WRITE (ITRUN,99997)
         WRITE (ITRUN,99996) G
      ENDIF

90500 RETURN

C     ALTERNATE RETURN IF MAXIMUM NUMBER OF FUNCTION EVALUATIONS
C     EXCEEDED.

91000 IF ( TRF .OR. TRG ) WRITE ( ITRUN,99999 )
      INDIC = LIMIT
      RETURN

C## F O R M A T S:

99992 FORMAT( ' [EVAL] INDIC (F,G,FG)=',4I3)
99993 FORMAT( ' EXCESSIVE ERROR IN GRADIENT ESTIMATION.')
99994 FORMAT( ' ERROR ESTIMATE IN GRADIENT ESTIMATION: ', G15.7/
     -        ' ESTIMATED FIGURES OF AGREEMENT:        ', G9.2 )
99995 FORMAT( ' ANALYTIC GRADIENT   ', G22.15, ' (COMPONENT ',I3,')'/
     -        ' ESTIMATED DERIVATIVE', G22.15                       )
99996 FORMAT( ' ', 5 G15.8 )
99997 FORMAT( ' (ZZEVAL) GRADIENT = ' )
99998 FORMAT( ' (ZZEVAL) FUNCTION = ', G26.16 )
99999 FORMAT(/' THE NUMBER OF FUNCTION EVALUATIONS ALLOWED HAS ',
     -         'BEEN EXCEEDED.')

C##                 E N D OF ZZEVAL.
                    END
