C     BEGIN OF PROGRAM *** DRIVER ***
      PROGRAM DRIVER
C
C     #################################################################
C
C     FIRST VERSION:    04.01.1998
C     PREVIOUS VERSION: 08.02.1999
C     LATEST VERSION:   11.02.1999
C
C     PROGRAMMED BY:
C     THOMAS WIEDER
C     E-MAIL: WIEDER@HRZ.UNI-KASSEL.DE
C     HTTP://WWW.UNI-KASSEL.DE/~WIEDER
C
C     ACKNOWLEDGMENTS:
C
C     MANY THANKS TO DR. T.R. HOPKINS, THE UNIVERSITY OF KENT, FOR
C     DETECTING A BUG AND FOR SEVERAL SUGGESTIONS TO IMPROVE THE FORTRAN
C     CODE (1998).
C     MANY THANKS TO DR. RICHARD HANSON, ACM-TOMS, FOR SEVERAL
C     SUGGESTIONS TO IMPROVE THE PROGRAM DESCRIPTION "NUMERICAL
C     HANKEL TRANSFORM BY THE FORTRAN PROGRAM HANKEL" (1998).
C
C     PURPOSE:
C
C     DRIVER ROUTINE.
C     EXAMPLE FOR USAGE OF SUBROUTINE *** HANKEL ***.
C
C     MOST IMPORTANT VARIABLES:
C
C     XI = POSITION AT WHICH THE HANKEL TRANSFORM *** HT ***
C          OF FUNCTION *** FUN *** HAS TO BE EVALUATED.
C     NU = ORDER OF HANKEL TRANSFORM
C        = ORDER OF BESSEL FUNCTION *** DHJNUX ***.
C
C                     ( infty
C     HT(FUN,NU,XI) =  |        (XI * X)**(1/2) JNUX(XI * X) FUN(X) dX
C                   0 )
C
C     HT = HANKEL TRANSFORM OF FUNCTION *** FUN ***.
C
C     THE FUNCTION *** FUN *** IS OF THE FORM: Y = FUN(X)
C
C     IMODE = LOGICAL FLAG
C     IMODE EQUAL TO 1 --> FUNCTION *** FUN *** IS GIVEN ANALYTICALLY
C                          WITHIN FUNCTION *** DHFUNC ***.
C     IMODE EQUAL TO 2 --> FUNCTION *** FUN *** IS GIVEN AS TABULATED
C                          DATA:
C                          X_1,Y_1, ..., X_N,Y_N, ..., X_NEND,Y_NEND.
C                          THESE DATA ARE PASSED TO *** HANKEL ***
C                          VIA AN INPUT FILE AND A CALL TO SUBROUTINE
C                          *** FUNK2 ***
C                          PRIOR TO THE CALL OF *** HANKEL ***.
C     IMODE EQUAL TO 3 --> FUNCTION *** FUN *** IS GIVEN AS TABULATED
C                          DATA
C                          X_1,Y_1, ..., X_N,Y_N, ..., X_NEND,Y_NEND.
C                          THESE DATA ARE PASSED TO *** HANKEL ***
C                          VIA A CALL TO SUBROUTINE *** FUNK3 ***
C                          PRIOR TO THE CALL OF *** HANKEL ***.
C
C     IVERBO = 0 --> SUBROUTINE *** HANKEL *** RUNS IN QUIET MODE.
C     IVERBO = 1 --> *** SUBROUTINE HANKEL *** SENDS SOME REPORT
C                    TO UNIT * (STANDARD OUTPUT).
C     IVERBO > 1 --> *** SUBROUTINE HANKEL *** SENDS SOME REPORT
C                    TO UNIT * AND ADDITIONAL REPORT TO
C                    FILE *** FILEER ***.
C     USUALLY YOU WILL CALL *** HANKEL *** WITH *** IVERBO *** = 0.
C     *** IVERBO *** 1 IS FOR DEBUGGING PURPOSES ONLY!
C
C     #################################################################
C
C      IMPLICIT NONE
C
C     -----------------------------------------------------------------
C
C     .. Scalars in Common ..
      CHARACTER*80 FILEER,FILEIN,FILEOU
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION DUMMY,HT,JNUXV,NU,SLTN,XI
      INTEGER ICON,IERR,IERRJ,IMODE,IVERBO,LUERR,LUIN,LUOUT,N,NEND
C     ..
C     .. Local Arrays ..
      DOUBLE PRECISION XDAT(1000),YDAT(1000)
C     ..
C     .. External Subroutines ..
      EXTERNAL DHINIT,DHJNUX,ERRMSS,FUNK1,FUNK2,FUNK3,HANKEL
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC DBLE,FLOAT
C     ..
C     .. Common blocks ..
      COMMON /NAMEN/FILEIN,FILEOU,FILEER
      integer i1mach, nout, nin
      nin =i1mach(1)
      nout = i1mach(2)
C     ..
   10 WRITE(NOUT,FMT=9000)
      WRITE(NOUT,FMT=9010)
C
C    -------------------------------------------------------------------
C
C     SET SOME PARAMETERS:
C
C     *** IVERBO *** CONTROLS THE VERBOSITY OF OUTPUT
      IVERBO = 1
C
C     *** LUIN *** IS THE LOGICAL NUMBER OF THE INPUT FILE
      LUIN = 66
C
C     *** LUOUT *** IS THE LOGICAL NUBER OF THE OUTPUT FILE
      LUOUT = 67
C
C     *** LUERR *** IS THE LOGICAL NUMBER OF THE ERROR LOG FILE
      LUERR = 68
C
C     SET *** DUMMY *** TO AN ARBITRARY REAL NUMBER
      DUMMY = 0.0D0
C
C     ------------------------------------------------------------------
C
C     OPEN FILES
C
      WRITE(NOUT,FMT=9020)
      READ(NIN,FMT=*,ERR=20,END=20) FILEER
      GO TO 30

   20 FILEER = 'hankel.err'
   30 WRITE(NOUT,FMT=9030) FILEER
      WRITE(NOUT,FMT=9040)
      READ(NIN,FMT=9050,ERR=40,END=40) FILEOU
      GO TO 50

   40 FILEOU = 'hankel.out'
   50 WRITE(NOUT,FMT=9060) FILEOU
C     FILE *** FILEER *** STORES MESSAGES FROM *** PROGRAM HANKEL ***.
C     SEARCH FOR ERROR MESSAGES IN FILE *** FILEER *** !
      OPEN (UNIT=LUERR,FILE=FILEER,STATUS='UNKNOWN',ACCESS='SEQUENTIAL',
     +     ERR=90)
C     FILE *** FILEOU *** STORES RESULTS FROM *** PROGRAM HANKEL ***.
      OPEN (UNIT=LUOUT,FILE=FILEOU,STATUS='UNKNOWN',ACCESS='SEQUENTIAL',
     +     ERR=90)
C
C    -------------------------------------------------------------------
C
C     INITIALIZE:
C
C     INITIALIZE THE LOGICAL NUMBERS FOR INPUT AND OUTPUT FILES
C
      CALL DHINIT(LUIN,LUOUT,LUERR)
C
C    -------------------------------------------------------------------
C
C     READ SOME INPUT DATA:
C
C     TELL WHERE THE FUNCTION *** FUN *** IS:
   60 WRITE(NOUT,FMT=9070)
      READ(NIN,FMT=*) IMODE
      IF (IMODE.EQ.0) STOP
C
C     *** NU *** IS THE ORDER OF THE HANKEL TRANSFORM:
      WRITE(NOUT,FMT=9080)
      READ(NIN,FMT=*) NU
C
C     SET PARAMETER *** XI ***:
      WRITE(NOUT,FMT=9090)
      READ(NIN,FMT=*) XI
C
C     ------------------------------------------------------------------
C
C     DO THE HANKEL TRANSFORM:
C
C     ------------------------------------------------------------------
C
      IF (IMODE.EQ.1) THEN
C
C     FUNCTION *** FUN *** IS CODED WITHIN SUBROUTINE *** DHFUNC ***,
C     THEREFORE, FIRST CALL SUBROUTINE *** FUNK1 *** !!!!!!!!!!!!!!!!!
C
          CALL FUNK1(DUMMY,DUMMY,IERR)
C
C     DO NOT CHANGE *** IMODE *** AFTER CALL OF *** FUNK1 *** !
C     THEN CALL SUBROUTINE *** HANKEL ***,
C     TO DO THE HANKEL TRANSFORM !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C
          CALL HANKEL(XI,NU,IVERBO,IERR,HT)
C
          IF (IERR.NE.0) THEN
              WRITE(NOUT,FMT=9100)
          END IF
C
C     ------------------------------------------------------------------
C
      ELSE IF (IMODE.EQ.2) THEN
C
C     DATA ARE PROVIDED IN AN INPUT FILE.
C
   70     WRITE(NOUT,FMT=9110)
          READ(NIN,FMT=9120,ERR=70) FILEIN
C
C     FIRST CALL SUBROUTINE *** FUNK2 ***
C     TO READ THE DATA FROM FILE!!!!!!!!!!!!!!!!!
C
          CALL FUNK2(DUMMY,DUMMY,IERR)
C
C     DO NOT CHANGE *** IMODE *** AFTER CALL OF *** FUNK2 *** !
C     THEN CALL SUBROUTINE *** HANKEL ***,
C     TO DO THE HANKEL TRANSFORM !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C
          CALL HANKEL(XI,NU,IVERBO,IERR,HT)
C
          IF (IERR.NE.0) THEN
              WRITE(NOUT,FMT=9100)
          END IF
C
C     ------------------------------------------------------------------
C
      ELSE IF (IMODE.EQ.3) THEN
C
C     DATA ARE TRANSFERED VIA CALL.
C
C     GENERATE DATA.
C
          NEND = 200
C
C     AS AN EXAMPLE:
C     FUN(X)=X**(NU+0.5D0) FOR X < 1
C     AND
C     FUN(X)=0 FOR X=> 1
          DO 80 N = 1,NEND
              XDAT(N) = DBLE(FLOAT(N)/FLOAT(NEND))
C     JUST AN EXAMPLE:
              YDAT(N) = XDAT(N)** (NU+0.5D0)
   80     CONTINUE
C
C     DATA ARE PASSED VIA CALL.
C     THEREFORE, FIRST CALL SUBROUTINE *** FUNK3 *** !!!!!!!!!!!!!!!!!
C
          CALL FUNK3(DUMMY,DUMMY,XDAT,YDAT,NEND,IERR)
C
C     DO NOT CHANGE *** IMODE *** AFTER CALL OF *** FUNK3 *** !
C     THEN CALL SUBROUTINE *** HANKEL ***,
C     TO DO THE HANKEL TRANSFORM !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C
          CALL HANKEL(XI,NU,IVERBO,IERR,HT)
C
          IF (IERR.NE.0) THEN
              WRITE(NOUT,FMT=9100)
          END IF
C
      ELSE
C     WRONG INPUT FOR *** IMODE ***, READ AGAIN:
          GO TO 60

      END IF
C
C     COMPARE ANALYTICAL AND NUMERICAL SOLUTION:
C
      IF (IERR.EQ.0 .AND. IVERBO.NE.0) THEN
C     CALCULATE THE ANALYTICAL SOLUTION FOR COMPARISION:
          IERRJ = 0
C
          CALL DHJNUX(XI,NU+1.0D0,IERRJ,JNUXV)
          IF (IERRJ.NE.0) THEN
              WRITE(NOUT,FMT=9130)
          END IF
C
          IF (IVERBO.GT.1) THEN
              SLTN = XI** (-0.5D0)*JNUXV
              WRITE(NOUT,FMT=9140) NU,XI,SLTN,HT
          END IF
C
      END IF
C
C    -------------------------------------------------------------------
C
      WRITE(NOUT,FMT=9150)
      READ(NIN,FMT=*) ICON
      IF (ICON.EQ.1) GO TO 10
C
      WRITE(NOUT,FMT=9160)
      WRITE(NOUT,FMT=9170)
C
C    -------------------------------------------------------------------
C
C     ERROR HANDLING:
C
      GO TO 100

   90 WRITE(NOUT,FMT=9180) FILEER
      IERR = -10
      CALL ERRMSS(IERR)
C
  100 STOP

 9000 FORMAT (/,' WELCOME TO   H A N K E L !',/,
     +       ' YOU ARE COMMUNICATING WITH   D R I V E R,',/,
     +       ' THE DRIVER ROUTINE TO   H A N K E L.')
 9010 FORMAT (/,' AS AN EXAMPLE, WE WILL NUMERICALLY HANKEL TRANSFORM',
     +       /,' FUN(X) = X**(NU+1/2) FOR X < 1 AND = 0 ELSEWHERE',/,
     +       ' AT POSITION XI FOR TRANSFORMATION ORDER NU.',/,
     +       ' THE ANALYTICAL SOLUTION IS',/,
     +       ' HT(FUN(X),NU,XI) = XI**(-1/2)*J(NU+1,XI)',/,
     +       ' (J(NU,X) IS THE BESSEL FUNCTION OF THE FIRST KIND.)')
 9020 FORMAT (/,' MESSAGE DRIVER: GIVE NAME FOR PROTOCOL FILE',/,
     +       ' (E.G. hankel.err):')
 9030 FORMAT (' MESSAGE DRIVER: ERROR LOG FILE =',A80)
 9040 FORMAT (/,' MESSAGE DRIVER: GIVE NAME FOR OUTPUT FILE',/,
     +       ' (E.G. hankel.out):')
 9050 FORMAT (A80)
 9060 FORMAT (' MESSAGE DRIVER: RESULT FILE =',A80)
 9070 FORMAT (/,' FUNCTION *** FUN(X) *** IS EITHER CODED WITHIN',/,
     +       '   FUNCTION *** DHFUNC *** = 1',/,' OR',/,
     +       ' TABULATED IN A FILE = 2 (VERY SLOW!)',/,' OR',/,
     +       ' PASSED AS TABLE VIA CALL TO HANKEL = 3 (SLOW!)',/,' OR ',
     +       /,' STOP = 0',/,' GIVE CHOICE (1 OR 2 OR 3 OR 0):',/)
 9080 FORMAT (/,' GIVE THE ORDER *** NU ***:',/,' (E.G. 1.0):',/)
 9090 FORMAT (/,' GIVE THE ARGUMENT *** XI ***:',/,' (E.G. 5.0):',/)
 9100 FORMAT (1X,' SERIOUS ERROR! HANKEL TRANSFORMATION FAILED!')
 9110 FORMAT (/,' MESSAGE DRIVER: GIVE NAME OF INPUT FILE:',/)
 9120 FORMAT (A80)
 9130 FORMAT (/,
     +' AN ERROR DURING FUNCTION EVALUATION                    OCCURED!'
     +       ,/,' THE HANKEL TRANSFORM HAS FAILED!')
 9140 FORMAT (/,
     +    ' RESULTS FOR TEST FUNCTION FUN(X)=X**(NU+0.5D0) FOR X = 1.0:'
     +       ,/,
     +       ' (IT IS ASSUMED THAT *** IFUNC=0 *** IN *** DHFUNC *** )',
     +       /,
     +      ' NU     XI     ANALYTICAL SOLUTION     NUMERICAL SOLUTION:'
     +       ,/,D17.9,3X,D17.9,3X,D17.9,3X,D17.9)
 9150 FORMAT (/,' DO YOU WANT ANOTHER CALCULATION = 1',/,
     +       '                            STOP = 2',/,' GIVE 1 OR 2:',/)
 9160 FORMAT (/,' A PROTOCOL HAS BEEN WRITTEN INTO A FILE,',/,
     +       ' RESULTS HAVE BEEN WRITTEN INTO A FILE.')
 9170 FORMAT (/,' GOODBYE!')
 9180 FORMAT (/,' MESSAGE DHINIT: ERROR ON OPENING FILE:',/,A80)
      END
C     END OF PROGRAM *** DRIVER ***
C
C     BEGIN OF SUBROUTINE *** DHFUNC ***
      SUBROUTINE DHFUNC(NU,X,Y,XLIMIT,IERR)
C
C     FIRST VERSION:    04.01.1998
C     PREVIOUS VERSION: 05.20.1999
C     LATEST VERSION:   06.02.1999
C
C     #################################################################
C
C     PROGRAMMED BY:
C     THOMAS WIEDER
C     E-MAIL: WIEDER@HRZ.UNI-KASSEL.DE
C     HTTP://WWW.UNI-KASSEL.DE/~WIEDER
C
C     PURPOSE:
C     SUBROUTINE *** DHFUNC *** PROVIDES THE FUNCTION *** FUN ***
C     WHICH SHOULD BE HANKEL-INVERTED BY SUBROUTINE *** HANKEL ***.
C
C     THE FORM OF SUBROUTINE *** DHFUNC *** IS:
C
C        SUBROUTINE DHFUNC(NU,X,Y,XLIMIT,LUERR,IERR)
C        IMPLICIT NONE
C        INTEGER IERR,LUIN,LUOUT,LUERR
C        DOUBLE PRECISION X,Y,NU,XASYMP,XLIMIT
C        COMMON /INOUTE/ LUIN,LUOUT,LUERR
C
C        ##############################################################
C
C        INPUT VARIABLES: NU, X,LUERR
C        OUPUT VARIABLES: Y,XLIMIT
C
C        ##############################################################
C
C        INITIALIZE THE ERROR FLAG.
C        *** IERR=0 *** MEANS THAT NO ERRORS OCCURRED SO FAR.
C        (OF COURSE YOU OMIT THE 'C'):
C
C        IERR=0
C
C        THE INFINITE INTEGRAL OVER THE FUNCTION *** FUN ***
C        WILL BE SPLIT INTO A FINITE INTEGRAL FROM 0 TO
C        *** XLIMIT *** AND AN INFINITE INTEGRAL FROM
C        *** XLIMIT *** TO INFINITY.
C        *** XLIMIT *** IS THE UPPER LIMIT OF THE FINITE INTEGRAL.
C        YOU NEED NOT TO SPECIFY *** XLIMIT ***. IN THAT CASE YOU
C        JUST GIVE A NEGATIVE VALUE, E.G. *** XLIMIT=-1.0D0 ***
C        (OF COURSE YOU OMIT THE 'C'):
C
C        XLIMIT=-1.0D0
C
C        HERE YOU GIVE YOUR FUNCTION *** FUN ***
C        Y=FUN(X,NU)
C        (OF COURSE YOU OMIT THE 'C'):
C
C        Y=...
C
C        IF AN ERROR HAS OCCURRED, THEN PUT *** IERR=1 *** !
C
C        RETURN
C        END
C
C     CALLING SEQUENCE:
C     CALLED BY FUNCTION *** FUNINT ***.
C
C     MOST IMPORTANT VARIABLES:
C
C     X = INPUT VALUE = INDEPENDENT VARIABLE OF *** FUN ***.
C     Y = OUTPUT VALUE = DEPENDENT VARIABLE OF *** FUN ***.
C     Y = FUN(X)
C     XLIMIT = UPPER LIMIT OF THE FINITE INTEGRAL.
C     LUERR = LOGICAL NUMBER OF ERROR LOG FILE
C     IERR = ERROR FLAG
C
C     IERR = 0 --> NO ERROR OCCURED DURING EVALUATION OF
C                  SUBROUTINE *** DHFUNC ***.
C
C     #################################################################
C
C      IMPLICIT NONE
C
C    -------------------------------------------------------------------
C
C     THE FOLLOWING VARIBALES ARE NECCESSARY
C     FOR DEBUGGING PURPOSES ONLY.
C     THE COMMON BLOCK *** SOLUTI *** IS NECEESSARY
C     FOR DEBUGGING PURPOSES ONLY.
C
C     ------------------------------------------------------------------
C
C     SET SOME PARAMETERS:
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION NU,X,XLIMIT,Y
      INTEGER IERR
C     ..
C     .. Scalars in Common ..
      DOUBLE PRECISION SLTN,XISLTN
      INTEGER IFUNC,LUERR,LUIN,LUOUT
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION JNUXV,PI,XASYMP
      INTEGER IERRJ
C     ..
C     .. External Subroutines ..
      EXTERNAL DHJNUX
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ASIN,DCOS,DSIN,EXP
      integer i1mach, nout
C     ..
C     .. Common blocks ..
      COMMON /INOUTE/LUIN,LUOUT,LUERR
      COMMON /SOLUTI/XISLTN,SLTN,IFUNC
C     ..
      nout = i1mach(2)
      IERR = 0
      XLIMIT = -1.0D0
C
      Y = 0.0D0
C
C     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C     *** 24.02.1998 *** A GOOD VALUE FOR *** XASYMP *** = 100.0
C     USUALLY, YOU DO NOT CHANGE THE FOLLOWING LINE:
      XASYMP = 100.0D0
C     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C
C     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C     !!! IFUNC MUST BE EQUAL TO ZERO IN ORDER TO HAVE YOUR FUNCTION
C     FUN HANKEL TRANSFORMED !!!
C     FOR DEBUGGING *** HANKEL ** ONE MAY CHOOSE ONE OF THE TEST
C     FUNCTIONS GIVEN BELOW BY SETTING *** IFUNC *** EQUAL TO VALUES
C     FROM 1 TO 4.
C     USUALLY, YOU DO NOT CHANGE THE FOLLOWING LINE:
      IFUNC = 0
C     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C
C     -----------------------------------------------------------------
C
C     GIVE THE FUNCTION *** FUN ***:
C
      IF (IFUNC.EQ.0) THEN
C     FOR *** IFUNC=0 *** YOUR FUNCTION WILL BE CALCULATED!
C
C     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C     HERE YOU HAVE TO TYPE YOUR FUNCTION:
C     HERE YOU HAVE TO TYPE YOUR FUNCTION:
C     HERE YOU HAVE TO TYPE YOUR FUNCTION:
C
C     WE PUT X**(NU+0.5D0) FOR X < 1 AS AN EXAMPLE:
C
C     WE SET *** XLIMIT *** = 1 BECAUSE OF THE DEFINITION RANGE.
          XLIMIT = 1.0D0
C
C     THIS IS THE FUNCTION *** FUN(X) ***:
          IF (X.LE.1.0D0) THEN
C     BE  V E R Y  CAUTIOUS ABOUT 0**0!
              IF (X.NE.0.0D0) Y = X** (NU+0.5D0)

          ELSE IF (X.GT.1.0D0) THEN
              Y = 0.0D0
          END IF
C
C     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C
C     ------------------------------------------------------------------
C
C     HERE SOME TEST FUNCTIONS:
C
      ELSE IF (IFUNC.EQ.1) THEN
C     AN ALGEBRAIC TEST FUNCTION.
C     SEE: ERDELYI, W. MAGNUS, F. OBERHETTINGER, F.G. TRICOMI: 8.5. (6)
C     THIS TEST FUNCTION IS EQUAL TO ZERO ABOVE *** X = 1.0 ***.
C     THEREFORE WE SET *** XLIMIT = 1.0D0 ***:
          XLIMIT = 1.0D0
          IF (NU.LE.-1.0D0) THEN
              WRITE (LUERR,FMT=9000) IFUNC,NU
              WRITE(NOUT,FMT=9000) IFUNC,NU
              STOP

          END IF

          IF (NU.LT.0.0D0 .AND. X.EQ.0.0D0) THEN
              Y = 0.0D0
              RETURN

          END IF

          IF (X.LE.1.0D0) THEN
C     BE  V E R Y  CAUTIOUS ABOUT 0**0!
              IF (X.GT.0.0D0) Y = X** (NU+0.5D0)

          ELSE IF (X.GT.1.0D0) THEN
              Y = 0.0D0
          END IF
C
C     SOLUTION:
C     XI**(-0.5) * JNUX(XI,NU+1)
          IERRJ = 0
C
          CALL DHJNUX(XISLTN,NU+1.0D0,IERRJ,JNUXV)
C
          IF (IERRJ.NE.0) THEN
              IERR = 1
              WRITE(NOUT,FMT=9010)
              WRITE (LUERR,FMT=9010)
          END IF

          SLTN = XISLTN** (-0.5D0)*JNUXV
C
C     ------------------------------------------------------------------
C
      ELSE IF (IFUNC.EQ.2) THEN
C     AN EXPONENTIAL TEST FUNCTION.
C     SEE: ERDELYI, W. MAGNUS, F. OBERHETTINGER, F.G. TRICOMI: 8.6. (1)
          XLIMIT = 100.0D0
          IF (NU.LE.-1.0D0) THEN
              WRITE (LUERR,FMT=9000) IFUNC,NU
              WRITE(NOUT,FMT=9000) IFUNC,NU
              STOP

          END IF

          IF (X.EQ.0.0D0) THEN
              Y = 0.0D0
              RETURN

          END IF
C     BE  V E R Y  CAUTIOUS ABOUT 0**0!
          IF (X.GT.0.0D0) Y = X** (-0.5D0)*EXP(-X)
C
C     SOLUTION:
C     XI**(0.5-NU) * (1.0+XI**2)**(-0.5) * ((1.0+XI**2)**(-0.5) - 1)**N
          SLTN = XISLTN** (0.5D0-NU)* (1.0D0+XISLTN**2)** (-0.5D0)*
     +           ((1.0D0+XISLTN**2)** (0.5D0)-1.0D0)**NU
C
C     ------------------------------------------------------------------
C
      ELSE IF (IFUNC.EQ.3) THEN
C     AN TRIGONOMETRIC TEST FUNCTION.
C     SEE: ERDELYI, W. MAGNUS, F. OBERHETTINGER, F.G. TRICOMI: 8.7. (1)
          PI = 3.141592653D0
          IF (XISLTN.EQ.1.0D0) THEN
              WRITE(NOUT,FMT=9020) IFUNC,XISLTN
              WRITE (LUERR,FMT=9020) IFUNC,XISLTN
              STOP

          END IF

          IF (NU.LE.-2.0D0) THEN
              WRITE (LUERR,FMT=9000) IFUNC,NU
              WRITE(NOUT,FMT=9000) IFUNC,NU
              STOP

          END IF

          IF (X.EQ.0.0D0) THEN
              Y = 0.0D0
              RETURN

          END IF

          IF (X.LE.XASYMP) THEN
C     BE  V E R Y  CAUTIOUS ABOUT 0**0!
              IF (X.GT.0.0D0) Y = X** (-0.5D0)*DSIN(X)

          ELSE IF (X.GT.XASYMP) THEN
              Y = 0.0D0
          END IF
C
C     SOLUTION:
C     FOR 0 < XI < 1:
C     COS(0.5 * PI * NU) * XI**(NU+0.5) * (1 - XI**2)**(-0.5) *
C     (1 + (1 - XI**2)**0.5)**(-NU)
C     FOR 1 < XI < INFINITY:
C     (1/NU) * XI**0.5 * SIN(NU * (ASIN(1/XI)))
          IF (XISLTN.LT.1.0D0) THEN
              SLTN = DCOS(0.5D0*PI*NU)*XISLTN** (NU+0.5D0)*
     +               (1.0D0-XISLTN**2)** (-0.5D0)*
     +               (1.0D0+ (1.0D0-XISLTN**2)**0.5D0)** (-NU)

          ELSE IF (XISLTN.GT.1.0D0) THEN
              SLTN = XISLTN**0.5D0* (XISLTN**2-1.0D0)** (-0.5D0)*
     +               DSIN(NU*ASIN(1.0D0/XISLTN))
          END IF
C
C     ------------------------------------------------------------------
C
      ELSE IF (IFUNC.EQ.4) THEN
C     A TEST FUNCTION CONTAINING A BESSEL FUNCTION.
C     SEE: ERDELYI, W. MAGNUS, F. OBERHETTINGER, F.G. TRICOMI:
C     8.11. (1)
          IF (NU.LE.-1.0D0) THEN
              WRITE (LUERR,FMT=9000) IFUNC,NU
              WRITE(NOUT,FMT=9000) IFUNC,NU
              STOP

          END IF

          IF (XISLTN.LE.1.0D0) THEN
              Y = 0.0D0
              SLTN = 0.0D0
              RETURN

          END IF

          IF (X.EQ.0.0D0) THEN
              Y = 0.0D0
              RETURN

          END IF
C
C      *** AMENDMENT SUGGESTED BY DR. T.R. HOPKINS, 15.04.1998 ***
          IERRJ = 0
C
          IF (X.LE.XASYMP) THEN
C
              CALL DHJNUX(X,NU-1.0D0,IERRJ,JNUXV)
C
C     BE  V E R Y  CAUTIOUS ABOUT 0**0!
              IF (X.GT.0.0D0) Y = X** (-0.5D0)*JNUXV
C
          ELSE IF (X.GT.XASYMP) THEN
              Y = 0.0D0
          END IF
C
          IF (IERRJ.NE.0) THEN
              WRITE(NOUT,FMT=9030) IERRJ
              WRITE(NOUT,FMT=9040) X
          END IF
C
C     SOLUTION:
C     ASSUME A=1 IN FORMULA 8.11.(1)!
C     FOR 0 < XI < 1:
C     0.0
C     FOR 1 < XI < INFTY:
C     Y**(-NU+0.5)
          IF (XISLTN.LT.1.0D0) THEN
              SLTN = 0.0D0

          ELSE IF (XISLTN.GT.1.0D0) THEN
              SLTN = XISLTN** (-NU+0.5D0)
          END IF
C
C     END OF IF-CLAUSE ***  IFUNC ***
      END IF
C
      RETURN

 9000 FORMAT (/,' MESSAGE FUNC:',/,' TEST FUNCTION',I6,
     +       ' NOT DEFINED FOR *** NU *** =',D17.9)
 9010 FORMAT (/,' MESSAGE FUNC:',/,
     +       ' AN ERROR DURING FUNCTION EVALUATIO  OCCURED!',/,
     +       ' THE HANKEL TRANSFORM HAS FAILED!')
 9020 FORMAT (/,' MESSAGE FUNC:',/,' HANKEL TRANSFORM OF TEST FUNCTION '
     +       ,I6,/,' NOT DEFINED FOR  *** XI *** =',D17.9)
 9030 FORMAT (/,' AN ERROR OCCURED DURING CALL OF *** DHJNUX ***',/,
     +       ' *** IERRJ *** =',I6)
 9040 FORMAT (' ARGUMENT *** X *** =',D17.9)
      END
C     END OF SUBROUTINE *** DHFUNC ***
C
C     BEGIN OF SUBROUTINE *** FUNK1 ***
      SUBROUTINE FUNK1(X,Y,IERR4)
C
C     FIRST VERSION:    04.01.1998
C     PREVIOUS VERSION: 06.10.1998
C     LATEST VERSION:   07.10.1998
C
C     ##################################################################
C
C     PROGRAMMED BY:
C     THOMAS WIEDER
C     E-MAIL: WIEDER@HRZ.UNI-KASSEL.DE
C
C     PURPOSE:
C     SUBROUTINE *** FUNK *** PROVIDES THE FUNCTION *** FUN ***
C     WHICH SHOULD BE HANKEL-INVERTED BY SUBROUTINE *** HANKEL ***.
C
C     *** FUN *** HAS TO BE PROVIDED ANALYTICALLY BY
C     SUBROUTINE *** DHFUNC *** !!!!!!!!!!!!!!!!!!!!
C
C     CALLING SEQUENCE:
C     CALLED BY SUBROUTINE *** HANKEL ***.
C     CALLS SUBROUTINE *** DHFUNC ***.
C
C     MOST IMPORTANT VARIABLES:
C
C     X = INDEPENDENT REAL VARIABLE OF *** FUN *** (INPUT).
C     Y = DEPENDENT REAL VARIABLE OF *** FUN *** (OUTPUT).
C     Y = FUN(X)
C
C     IDAT (OUTPUT)
C     IDAT EQUAL TO 1 --> FUNCTION *** FUN *** IS GIVEN ANALYTICALLY
C     IDAT EQUAL TO 2 --> FUNCTION *** FUN *** IS GIVEN AS TABULATED
C                         DATA IN A FILE
C     IDAT EQUAL TO 3 --> FUNCTION *** FUN *** IS GIVEN AS TABULATED
C                         DATA PASSED BY CALL TO *** HANKEL ***
C
C     IERR4 = ERROR FLAG (OUTPUT)
C     IERR4 EQUAL TO 0 ON RETURN --> NO ERROR OCCURED.
C     IERR4 EQUAL TO 1 ON RETURN --> AN ERROR OCCURED DURING
C                                   EVALUATION OF FUNCTION *** FUN ***
C                                   WITHIN SUBROUTINE *** DHFUNC ***.
C
C     ##################################################################
C
C      IMPLICIT NONE
C
C     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C     FUNCTION *** FUN *** IS GIVEN ANALYTICALLY WITHIN
C     SUBROUTINE *** DHFUNC ***.
C     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C
C     .. Parameters ..
      INTEGER ITEST
      PARAMETER (ITEST=0)
C     ..
C     .. Scalar Arguments ..
      DOUBLE PRECISION X,Y
      INTEGER IERR4
C     ..
C     .. Scalars in Common ..
      DOUBLE PRECISION NU,XI,XLIMIT
      INTEGER IDAT,IREAD,LUERR,LUIN,LUOUT
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION RESULT
      INTEGER IERR5,N,NEND
C     ..
C     .. External Subroutines ..
      EXTERNAL DHFUNC
C     ..
C     .. Common blocks ..
      COMMON /ARGUME/XI,NU,XLIMIT
      COMMON /INOUTE/LUIN,LUOUT,LUERR
      COMMON /LOGIK/IDAT,IREAD
      integer i1mach, nout
      nout = i1mach(2)
C     ..
      IF (IREAD.NE.1) THEN
C
C     SAVE INPUT
C
          IREAD = 1
      END IF
C
      IERR4 = 0
      IERR5 = 0
C
      CALL DHFUNC(NU,X,RESULT,XLIMIT,IERR5)
C
      IF (IERR5.NE.0) THEN
          WRITE (LUERR,FMT=9000) IERR5
          WRITE(NOUT,FMT=9000) IERR5
          WRITE (LUERR,FMT=9010) X
          IERR4 = 1
          RETURN

      END IF
C
      IF (ITEST.EQ.1) THEN
          NEND = 100
          WRITE(NOUT,FMT=9020)
          X = 0.1D0
          DO 10 N = 1,NEND
              WRITE(NOUT,FMT=*) N,X,RESULT
              X = X + 0.1D0
   10     CONTINUE
      END IF
C
C     SUCCESSFUL RETURN FROM *** FUNK1 ***
      IDAT = 1
      Y = RESULT
C
      RETURN

 9000 FORMAT (/,' MESSAGE FUNK1:',/,
     +       ' ERROR ON CALCULATING FUNCTION VA  LUE!',/,
     +      ' ERROR IN SUBROUTINE *** DHFUNC *** OR SIMILAR SUBROUTINE!'
     +       ,/,' IERR5 = ',I6)
 9010 FORMAT (1X,' FUNCTION ARGUMENT X =',D17.9)
 9020 FORMAT (/,
     +       ' MESSAGE FUNK1: TEST RUN FOR SUBROUTINE *** DHFUNC **  *:'
     +       )
      END
C     END OF SUBROUTINE *** FUNK1 ***
C
C     BEGIN OF SUBROUTINE *** FUNK2 ***
      SUBROUTINE FUNK2(X,Y,IERR4)
C
C     FIRST VERSION:    04.01.1998
C     PREVIOUS VERSION: 10.10.1998
C     LATEST VERSION:   13.10.1998
C
C     ##################################################################
C
C     PROGRAMMED BY:
C     THOMAS WIEDER
C     E-MAIL: WIEDER@HRZ.UNI-KASSEL.DE
C
C     PURPOSE:
C     SUBROUTINE *** FUNK *** PROVIDES THE FUNCTION *** FUN ***
C     WHICH SHOULD BE HANKEL-INVERTED BY SUBROUTINE *** HANKEL ***.
C
C     *** FUN *** HAS TO BE PROVIDED AS TABULATED DATA X,Y
C     IN A FILE. !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C
C     CALLING SEQUENCE:
C     CALLED BY SUBROUTINE *** HANKEL ***.
C     CALLS SUBROUTINE *** DIVDIF ***.
C
C     MOST IMPORTANT VARIABLES:
C
C     X = INDEPENDENT REAL VARIABLE OF *** FUN *** (INPUT).
C     Y = DEPENDENT REAL VARIABLE OF *** FUN *** (OUTPUT).
C     Y = FUN(X)
C
C     NEND = NUMBER OF DATA PAIRS (INPUT).
C     NMAX = MAXIMUM NUMBER OF ALLOWED FUNCTION VALUES (OUTPUT).
C
C     IDAT (OUTPUT)
C     IDAT EQUAL TO 1 --> FUNCTION *** FUN *** IS GIVEN ANALYTICALLY
C     IDAT EQUAL TO 2 --> FUNCTION *** FUN *** IS GIVEN AS TABULATED
C                         DATA IN A FILE
C     IDAT EQUAL TO 3 --> FUNCTION *** FUN *** IS GIVEN AS TABULATED
C                         DATA PASSED BY CALL TO *** HANKEL ***
C
C     IERR4 = ERROR FLAG (OUTPUT)
C     IERR4 EQUAL TO 0 ON RETURN --> NO ERROR OCCURED.
C     IERR4 EQUAL TO 1 ON RETURN --> AN ERROR OCCURED DURING
C                                   EVALUATION OF FUNCTION *** FUN ***
C
C     ##################################################################
C
C      IMPLICIT NONE
C
C     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C     FUNCTION *** FUN *** IS GIVEN AS DATA TABLE IN AN INPUT FILE.
C     THE FORMAT OF THE DATA IS:
C     X_1      Y_1
C     ...      ...
C     X_N      Y_N
C     ...      ...
C     X_NEND   Y_NEND
C     I.E. EACH LINE CONTAINS ONE DATA POINT (X,Y).
C     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C
C     .. Parameters ..
      INTEGER NMAX
      PARAMETER (NMAX=1000)
      INTEGER ITEST
      PARAMETER (ITEST=0)
C     ..
C     .. Scalar Arguments ..
      DOUBLE PRECISION X,Y
      INTEGER IERR4
C     ..
C     .. Scalars in Common ..
      DOUBLE PRECISION NU,XDATMA,XDATMI,XI,XLIMIT,YDATMA,YDATMI
      INTEGER IDAT,IREAD,LUERR,LUIN,LUOUT,NEND
C     ..
C     .. Arrays in Common ..
      DOUBLE PRECISION XDAT(NMAX),YDAT(NMAX)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION RESULT
      INTEGER IDEG,IERR6,IERRRW,N
C     ..
C     .. Local Arrays ..
      DOUBLE PRECISION TABLE(NMAX,NMAX)
C     ..
C     .. External Subroutines ..
      EXTERNAL DHRW,DIVDIF
C     ..
C     .. Common blocks ..
      COMMON /ARGUME/XI,NU,XLIMIT
      COMMON /DATEN/XDAT,YDAT,XDATMI,XDATMA,YDATMI,YDATMA,NEND
      COMMON /INOUTE/LUIN,LUOUT,LUERR
      COMMON /LOGIK/IDAT,IREAD
      integer i1mach, nout
      nout = i1mach(2)
C     ..
      IERR4 = 0
C
C     READ DATA FROM FILE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C
      IF (IREAD.NE.1) THEN
C
C     *** LUIN *** = LOGICAL NUMBER OF INPUT FILE
C
          CALL DHRW(IERRRW)
C     READ DATA ONLY ONCE:
          IREAD = 1
C
          IF (IERRRW.NE.0) THEN
              WRITE(NOUT,FMT=9000) IERRRW
              WRITE (LUERR,FMT=9000) IERRRW
C     AN ERROR OCCURED, DO NOT CONTINUE, GO TO STOP!
              WRITE (LUERR,FMT=9010)
              DO 10 N = 1,NEND
                  WRITE (LUERR,FMT=*) N,XDAT(N),YDAT(N)
   10         CONTINUE
              STOP

          END IF
C
          IF (ITEST.EQ.1) THEN
              WRITE(NOUT,FMT=9020)
              DO 20 N = 1,NEND
                  WRITE(NOUT,FMT=*) N,XDAT(N),YDAT(N)
   20         CONTINUE
          END IF
C
C     SUCCESSFUL RETURN FROM *** FUNK2 ***, DATA HAVE BEEN READ!
          IDAT = 2
          RETURN
C
C     END OF IF-CLAUSE *** IREAD ***
      END IF
C
C     INTERPOLATE FUNCTION VALUE FOR ARGUMENT *** X *** !!!!!!!!!!!!!!!!
C
C     DO NOT PERFORM THE DEFINITE INTEGRAL ABOVE THE TABULATED
C     DATA RANGE:
      XLIMIT = XDATMA
C
      IF (X.GE.XDATMI .AND. X.LE.XDATMA) THEN
C
C     *** X*** LIES WITHIN TABULATED DATA RANGE. TRY THE INTERPOLATION!
C
C     CHOOSE *** IDEG *** NOT TOO HIGH!!!
          IDEG = 5
          IERR6 = 0
C
          CALL DIVDIF(XDAT,YDAT,NEND,IDEG,IDEG,TABLE,X,RESULT,IERR6)
C      CALL DIVDIF (XDAT,YDAT,NEND,IDEG,IDEG,X,RESULT,IERR6)
C
          IF (IERR6.NE.0) THEN
              WRITE (LUERR,FMT=9030) IERR6
              WRITE(NOUT,FMT=9030) IERR6
              WRITE (LUERR,FMT=9040) X
              WRITE (LUERR,FMT=9050) NEND
              IERR4 = 2
              WRITE (LUERR,FMT=9060) IERR4
              RETURN

          END IF
C
C     SUCCESSFUL RETURN FROM *** DIVDIF ***
          Y = RESULT
C
      ELSE
C
C     *** X *** LIES OUTSIDE THE TABULATED DATA RANGE!
C     NO INTERPOLATION POSSIBLE!
C     THIS CASE WILL LEGALLY HAPPEN, SINCE THE HANKEL TRANSFORM
C     INCLUDES AN INTEGRATION TO INFINITY.
C
          Y = 0.0D0
C
C     END IF-CLAUSE *** X.LE.XDATMA ***.
      END IF
C
      RETURN

 9000 FORMAT (/,' MESSAGE FUNK2:',/,' ERROR ON READING DATA! IERRRW =',
     +       I6)
 9010 FORMAT (/,' MESSAGE FUNK2: DATA READ FROM INPUT FILE:')
 9020 FORMAT (' MESSAGE FUNK2: DATA READ FROM FILE:')
 9030 FORMAT (/,' MESSAGE FUNK2:',/,
     +       ' ERROR ON INTERPOLATING FUNCTION   VALUE!',/,
     +       ' ERROR IN SUBROUTINE *** DIVDIF ***!',/,' IERR6 = ',I6)
 9040 FORMAT (1X,' ARGUMENT *** X *** =',D17.9)
 9050 FORMAT (1X,' NUMBER OF DATA PAIRS *** NEND *** =',I6)
 9060 FORMAT (/,' MESSAGE FUNK2: IERR4 =',I6)
      END
C     END OF SUBROUTINE FUNK2 ***
C
C     BEGIN OF SUBROUTINE *** FUNK3 ***
      SUBROUTINE FUNK3(X,Y,XDATIN,YDATIN,NENDIN,IERR4)
C
C     FIRST VERSION:    04.01.1998
C     PREVIOUS VERSION: 10.10.1998
C     LATEST VERSION:   13.10.1998
C
C     ##################################################################
C
C     PROGRAMMED BY:
C     THOMAS WIEDER
C     E-MAIL: WIEDER@HRZ.UNI-KASSEL.DE
C
C     PURPOSE:
C     SUBROUTINE *** FUNK *** PROVIDES THE FUNCTION *** FUN ***
C     WHICH SHOULD BE HANKEL-INVERTED BY SUBROUTINE *** HANKEL ***.
C
C     *** FUN *** HAS TO BE PROVIDED AS FIELDS X(N),Y(N) PASSED TO
C     *** FUNK3 *** BY THE CALLING PROGRAM. !!!!!!!!!!!!!!!!!!!!
C
C     CALLING SEQUENCE:
C     CALLED BY SUBROUTINE *** HANKEL ***.
C     CALLS SUBROUTINE *** DIVDIF ***.
C
C     MOST IMPORTANT VARIABLES:
C
C     X = INDEPENDENT REAL VARIABLE OF *** FUN *** (INPUT).
C     Y = DEPENDENT REAL VARIABLE OF *** FUN *** (OUTPUT).
C     Y = FUN(X)
C
C     NEND = NUMBER OF DATA PAIRS (INPUT).
C     NMAX = MAXIMUM NUMBER OF ALLOWED FUNCTION VALUES (OUTPUT).
C
C     IDAT (OUTPUT)
C     IDAT EQUAL TO 1 --> FUNCTION *** FUN *** IS GIVEN ANALYTICALLY
C     IDAT EQUAL TO 2 --> FUNCTION *** FUN *** IS GIVEN AS TABULATED
C                         DATA IN A FILE
C     IDAT EQUAL TO 3 --> FUNCTION *** FUN *** IS GIVEN AS TABULATED
C                         DATA PASSED BY CALL TO *** HANKEL ***
C
C     IERR4 = ERROR FLAG (OUTPUT)
C     IERR4 EQUAL TO 0 ON RETURN --> NO ERROR OCCURED.
C     IERR4 EQUAL TO 1 ON RETURN --> AN ERROR OCCURED DURING
C                                   EVALUATION OF FUNCTION *** FUN ***
C
C     ##################################################################
C
C      IMPLICIT NONE
C
C     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C     FUNCTION *** FUN *** IS GIVEN AS DATA TABLE IN AN INPUT FILE.
C     THE FORMAT OF THE DATA IS:
C     X_1      Y_1
C     ...      ...
C     X_N      Y_N
C     ...      ...
C     X_NEND   Y_NEND
C     I.E. EACH LINE CONTAINS ONE DATA POINT (X,Y).
C     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C
C     .. Parameters ..
      INTEGER NMAX
      PARAMETER (NMAX=1000)
      INTEGER ITEST
      PARAMETER (ITEST=0)
C     ..
C     .. Scalar Arguments ..
      DOUBLE PRECISION X,Y
      INTEGER IERR4,NENDIN
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION XDATIN(NMAX),YDATIN(NMAX)
C     ..
C     .. Scalars in Common ..
      DOUBLE PRECISION NU,XDATMA,XDATMI,XI,XLIMIT,YDATMA,YDATMI
      INTEGER IDAT,IREAD,LUERR,LUIN,LUOUT,NEND
C     ..
C     .. Arrays in Common ..
      DOUBLE PRECISION XDAT(NMAX),YDAT(NMAX)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION RESULT
      INTEGER IDEG,IERR6,N
C     ..
C     .. Local Arrays ..
      DOUBLE PRECISION TABLE(NMAX,NMAX)
C     ..
C     .. External Functions ..
      DOUBLE PRECISION D1MACH
      EXTERNAL D1MACH
C     ..
C     .. External Subroutines ..
      EXTERNAL DIVDIF
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC MAX,MIN
C     ..
C     .. Common blocks ..
      COMMON /ARGUME/XI,NU,XLIMIT
      COMMON /DATEN/XDAT,YDAT,XDATMI,XDATMA,YDATMI,YDATMA,NEND
      COMMON /INOUTE/LUIN,LUOUT,LUERR
      COMMON /LOGIK/IDAT,IREAD
      integer i1mach, nout
      nout = i1mach(2)
C     ..
      IERR4 = 0
C
C     ACCEPT DATA FROM CALLING PROGRAM !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C
      IF (IREAD.NE.1) THEN
C
C     SAVE INPUT DATA
C
          NEND = NENDIN
          DO 10 N = 1,NEND
              XDAT(N) = XDATIN(N)
              YDAT(N) = YDATIN(N)
   10     CONTINUE
C
          WRITE (LUERR,FMT=9000)
C     READ DATA ONLY ONCE:
          IREAD = 1
C
C     CHECK DATA
C
          IF (NEND.GT.NMAX) THEN
              WRITE(NOUT,FMT=9010) NEND
              STOP

          END IF
C
          XDATMI = D1MACH(2)
          XDATMA = -D1MACH(2)
          YDATMI = D1MACH(2)
          YDATMA = -D1MACH(2)
          DO 20 N = 1,NEND
              XDATMI = MIN(XDATMI,XDAT(N))
              XDATMA = MAX(XDATMA,XDAT(N))
              YDATMI = MIN(YDATMI,YDAT(N))
              YDATMA = MAX(YDATMA,YDAT(N))
   20     CONTINUE
C
          IF (XDATMI.LT.-D1MACH(2)) THEN
              WRITE(NOUT,FMT=9020)
              WRITE (LUERR,FMT=9020)
              STOP

          END IF

          IF (XDATMI.GT.D1MACH(2)) THEN
              WRITE(NOUT,FMT=9030)
              WRITE (LUERR,FMT=9030)
              STOP

          END IF

          IF (YDATMI.LT.-D1MACH(2)) THEN
              WRITE(NOUT,FMT=9040)
              WRITE (LUERR,FMT=9040)
              STOP

          END IF

          IF (YDATMI.GT.D1MACH(2)) THEN
              WRITE(NOUT,FMT=9050)
              WRITE (LUERR,FMT=9050)
              STOP

          END IF
C
C     SUCCESSFUL RETURN FROM *** DHFUNC3 ***, DATA HAVE BEEN READ!
          IDAT = 3
          RETURN
C
C     END OF IF-CLAUSE *** IREAD ***
      END IF
C
      IF (ITEST.EQ.1) THEN
          WRITE(NOUT,FMT=9060)
          DO 30 N = 1,NEND
              WRITE(NOUT,FMT=*) N,XDAT(N),YDAT(N)
   30     CONTINUE
      END IF
C
C     INTERPOLATE FUNCTION VALUE FOR ARGUMENT *** X *** !!!!!!!!!!!!!!!!
C
C     DO NOT PERFORM THE DEFINITE INTEGRAL ABOVE THE TABULATED
C     DATA RANGE:
      XLIMIT = XDATMA
C
      IF (X.GE.XDATMI .AND. X.LE.XDATMA) THEN
C
C     *** X*** LIES WITHIN TABULATED DATA RANGE. TRY THE INTERPOLATION!
C
C     CHOOSE *** IDEG *** NOT TOO HIGH!!!
          IDEG = 5
          IERR6 = 0
C
          CALL DIVDIF(XDAT,YDAT,NEND,IDEG,IDEG,TABLE,X,RESULT,IERR6)
C      CALL DIVDIF (XDAT,YDAT,NEND,IDEG,IDEG,X,RESULT,IERR6)
C
          IF (IERR6.NE.0) THEN
              WRITE (LUERR,FMT=9070) IERR6
              WRITE(NOUT,FMT=9070) IERR6
              WRITE (LUERR,FMT=9080) X
              WRITE (LUERR,FMT=9090) NEND
              IERR4 = 1
              WRITE (LUERR,FMT=9100) IERR4
              RETURN

          END IF
C
C     SUCCESSFUL RETURN FROM *** DIVDIF ***
          Y = RESULT
C
      ELSE
C
C     *** X *** LIES OUTSIDE THE TABULATED DATA RANGE!
C     NO INTERPOLATION POSSIBLE!
C     THIS CASE WILL LEGALLY HAPPEN, SINCE THE HANKEL TRANSFORM
C     INCLUDES AN INTEGRATION TO INFINITY.
C
          Y = 0.0D0
C
C     END IF-CLAUSE *** X.LE.XDATMA ***.
      END IF
C
      RETURN

 9000 FORMAT (/,' MESSAGE FUNK3: ',/,
     +  ' DATA ** X(N), Y(N) *** HAVE BEE  N READ FROM CALLING PROGRAM!'
     +       )
 9010 FORMAT (/,' MELDUNG RW:',/,' ERROR! MAXIMAL NUMBER =',I6,
     +       ' OF DATA PAIRS EXCEEDED!')
 9020 FORMAT (/,' MESSAGE RW: TOO SMALL VALUE FOR ABCISSA *** X ***!',/,
     +       ' PLEASE CHECK YOUR DATA!')
 9030 FORMAT (/,' MESSAGE FUNK3: TOO BIG VALUE FOR ABCISSA *** X ***!',
     +       /,'      PLEASE CHECK YOUR DATA!')
 9040 FORMAT (/,
     +       ' MESSAGE FUNK3: TOO SMALL VALUE FOR ORDINATE *** Y **  *!'
     +       ,/,' PLEASE CHECK YOUR DATA!')
 9050 FORMAT (/,
     +       ' MESSAGE FUNK3: TOO BIG VALUE FOR ORDINATE *** X ***!  ',
     +       /,' PLEASE CHECK YOUR DATA!')
 9060 FORMAT (' MESSAGE FUNK3: DATA RECEIVED FROM CALLING PROGRAM:')
 9070 FORMAT (/,' MESSAGE FUNK3:',/,
     +       ' ERROR ON INTERPOLATING FUNCTION   VALUE!',/,
     +       ' ERROR IN SUBROUTINE *** DIVDIF ***!',/,' IERR6 = ',I6)
 9080 FORMAT (1X,' ARGUMENT *** X *** =',D17.9)
 9090 FORMAT (1X,' NUMBER OF DATA PAIRS *** NEND *** =',I6)
 9100 FORMAT (/,' MESSAGE FUNK3: IERR4 =',I6)
      END
