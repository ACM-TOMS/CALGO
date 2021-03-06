      SUBROUTINE BBDFLT ( PFREQ,  MAXF )

C## A R G U M E N T S:
                     INTEGER    PFREQ,  MAXF

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
C>RCS $HEADER: DFLT.F,V 1.12 91/12/31 14:52:48 BUCKLEY EXP $
C>RCS $LOG:     DFLT.F,V $
C>RCS REVISION 1.12  91/12/31  14:52:48  BUCKLEY
C>RCS FINAL SUBMISSION TO TOMS
C>RCS
C>RCS REVISION 1.11  91/11/22  11:28:33  BUCKLEY
C>RCS FINAL SUBMISSION TO TOMS
C>RCS
C>RCS REVISION 1.10  90/07/31  10:49:33  BUCKLEY
C>RCS ADDED REVISED BLAS.
C>RCS
C>RCS REVISION 1.9  89/06/30  13:12:45  BUCKLEY
C>RCS PREPARING SUBMITTED VERSION OF MT
C>RCS
C>RCS REVISION 1.3  89/05/18  12:39:15  BUCKLEY
C>RCS FINAL TEST OF MT BEFORE SUBMITTING
C>RCS
C>RCS REVISION 1.2  89/05/15  14:55:21  BUCKLEY
C>RCS INITIAL INSTALLATION OF MT INTO RCS FORM.
C>RCS
C>RCS REVISION 1.1  89/01/17  16:54:27  BUCKLEY
C>RCS INITIAL REVISION
C>RCS
C
C## D E S C R I P T I O N:
C
C     THIS ROUTINE FIRST OBTAINS THE DEFAULT VALUES FOR INITIALIZING
C     THE ROUTINES  ZZPRNT,  ZZEVAL,  ZZTERM  AND  BBLNIR.  IT THEN
C     CALLS ENTRY POINTS IN EACH OF THESE ROUTINES TO SET THE INITIAL
C     VALUES NEEDED IN THOSE ROUTINES TO THOSE DEFAULT VALUES.
C
C## E N T R Y   P O I N T S: THE NATURAL ENTRY BBDFLT
C
C## S U B R O U T I N E S:
C
C     BBVALS                TO OBTAIN THE DEFAULT VALUES
C     ZZP1ST ZZP2ST         ENTRY POINTS TO ZZPRNT
C     ZZTSET ZZESET ZZESRT: ENTRY POINTS TO ZZTERM, ZZEVAL
C     BBLSET                ENTRY POINT  TO BBLNIR.
C
C## P A R A M E T E R S:

      INTEGER     NL1,     NL2,     NLINF
      PARAMETER ( NL1 = 1, NL2 = 2, NLINF = 3 )

      INTEGER     NQUITS
      PARAMETER ( NQUITS = 4 )

      INTEGER     PGRAD,     PSTEP,     PSHXG,     PFUNC
      PARAMETER ( PGRAD = 1, PSTEP = 2, PSHXG = 3, PFUNC = 4 )

      INTEGER     NINTS,      NLOGS,      NREALS,     NTRACF
      PARAMETER ( NINTS = 14, NLOGS = 32, NREALS = 2, NTRACF = 15 )

      INTEGER     XDRVMD,     XNORM,      XSCALE,     XLTRCU
      PARAMETER ( XDRVMD = 1, XNORM = 2,  XSCALE = 3, XLTRCU = 4 )

      INTEGER     XETRCU,     XPTRCU,     XTTRCU
      PARAMETER ( XETRCU = 5, XPTRCU = 6, XTTRCU = 7 )

      INTEGER     XMETH,      XQUADN,     XALPS1,      XSCGMM
      PARAMETER ( XMETH = 8,  XQUADN = 9, XALPS1 = 10, XSCGMM = 11 )

      INTEGER     XHTEST,     XUPDTT,      XSTSTP
      PARAMETER ( XHTEST = 12,XUPDTT = 13, XSTSTP = 14 )

      INTEGER     XTRACE
      PARAMETER ( XTRACE = 1 )

      INTEGER     XTRF,       XTRG,        XTTRCE,      XTRTST
      PARAMETER ( XTRF =  16, XTRG =   17, XTTRCE = 18, XTRTST = 19 )

      INTEGER     XGRAD,      XPOINT,      XTGRAD
      PARAMETER ( XGRAD = 20, XPOINT = 21, XTGRAD = 22 )

      INTEGER     XTSTEP,     XTSHXG,      XTFUNC,      XRELF
      PARAMETER ( XTSTEP = 23,XTSHXG = 24, XTFUNC = 25, XRELF = 26 )

      INTEGER     XRELG,      XFQUAD,      XDIAGL
      PARAMETER ( XRELG = 27, XFQUAD = 28, XDIAGL = 29 )

      INTEGER     XSHNNO,     XFRMRS,      XFRCEF
      PARAMETER ( XSHNNO = 30,XFRMRS = 31, XFRCEF = 32 )

      INTEGER     XRO,        XBETA
      PARAMETER ( XRO = 1,    XBETA = 2 )

C## L O C A L   D E C L:
                        CHARACTER*(NQUITS) TESTS
                        INTEGER           INTS(NINTS)
                        LOGICAL           LOGS(NLOGS)
                        DOUBLE PRECISION  REALS(NREALS)
C!!!!                   REAL              REALS(NREALS)

C## S A V E:                 NONE SELECTED.
C## E Q U I V A L E N C E S: NONE ARE DEFINED.
C## C O M M O N:             NONE IS DEFINED.
C## D A T A:                 NONE ARE SET.
C##                                                E X E C U T I O N
C##                                                E X E C U T I O N

C-----OBTAIN DEFAULTS.
                     CALL BBVALS ( INTS, LOGS, REALS )
C-ZZEVAL.
        CALL ZZESET ( LOGS(XTRF),LOGS(XTRG),LOGS(XTRTST),INTS(XETRCU))
        CALL ZZESRT ( INTS(XSCALE), INTS(XDRVMD), MAXF, 1 )

C-ZZPRNT.
        CALL ZZP1ST (INTS(XPTRCU),LOGS(XGRAD),LOGS(XPOINT), PFREQ )
        CALL ZZP2ST (INTS(XPTRCU),LOGS(XGRAD),LOGS(XPOINT),0,.FALSE.,0)
C-ZZTERM.
        TESTS = 'FFFF'
        IF (LOGS(XTGRAD)) TESTS(PGRAD:PGRAD) = 'T'
        IF (LOGS(XTSTEP)) TESTS(PSTEP:PSTEP) = 'T'
        IF (LOGS(XTSHXG)) TESTS(PSHXG:PSHXG) = 'T'
        IF (LOGS(XTFUNC)) TESTS(PFUNC:PFUNC) = 'T'
        CALL ZZTSET ( INTS(XNORM), TESTS, LOGS(XTTRCE), INTS(XTTRCU) )

C-BBLNIR.
        CALL BBLSET ( INTS(XMETH),  INTS(XQUADN), INTS(XALPS1),
     -                INTS(XSTSTP), INTS(XSCGMM), INTS(XHTEST),
     -                INTS(XUPDTT),
     -                REALS(XRO),   REALS(XBETA),
     -                LOGS(XFQUAD), LOGS(XDIAGL), LOGS(XSHNNO),
     -                LOGS(XFRMRS), LOGS(XFRCEF), LOGS(XRELF),
     -                LOGS(XRELG),
     -                INTS(XLTRCU), LOGS(XTRACE)              )
      GOTO 90000

C## E X I T
90000      RETURN

C## F O R M A T S:  NONE ARE DEFINED.
C##                 E N D         OF BBDFLT.
                    END
