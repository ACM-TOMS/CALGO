C     BEGIN OF SUBROUTINE *** DHINIT ***
      SUBROUTINE DHINIT(LUININ,LUOUIN,LUERIN)
C
C     #################################################################
C
C     FIRST VERSION:    10.10.1998
C     PREVIOUS VERSION: 10.12.1998
C     LATEST VERSION:   03.02.1999
C
C     PROGRAMMED BY:
C     THOMAS WIEDER
C     E-MAIL: WIEDER@HRZ.UNI-KASSEL.DE
C     HTTP://WWW.UNI-KASSEL.DE/~WIEDER
C
C     PURPOSE:
C
C     STORE LOGICAL FILE NUMBERS IN COMMON BLOCK /INOUTE/
C     TO TRANSFER THESE NUMBERS TO ALL OTHER SUBROUTINES
C     VIA COMMON BLOCK /INOUTE/.
C
C     #################################################################
C
C      IMPLICIT NONE
C
C     -----------------------------------------------------------------
C
C     INITIALIZE SOME VARAIBLES:
C
C     .. Scalar Arguments ..
      INTEGER LUERIN,LUININ,LUOUIN
C     ..
C     .. Scalars in Common ..
      INTEGER IDAT,IREAD,LUERR,LUIN,LUOUT
      CHARACTER*80 FILEER,FILEIN,FILEOU
C     ..
C     .. Common blocks ..
      COMMON /INOUTE/LUIN,LUOUT,LUERR
      COMMON /LOGIK/IDAT,IREAD
      COMMON /NAMEN/FILEIN,FILEOU,FILEER
C     ..
      IREAD = 1
C
C     SAVE INPUT:
C
      LUIN = LUININ
      LUERR = LUERIN
      LUOUT = LUOUIN
C
      RETURN

      END
C     END OF SUBROUTINE *** DHINIT ***
C
C     BEGIN OF SUBROUTINE *** HANKEL ***
      SUBROUTINE HANKEL(XI,NU,IVERBO,IERR1,HT)
C
C     ##################################################################
C
C     FIRST VERSION:    04.01.1998
C     PREVIOUS VERSION: 14.10.1998
C     LATEST VERSION:   12.12.1998
C
C     PROGRAMMED BY:
C     THOMAS WIEDER
C     TECHNICAL UNIVERSITY DARMSTADT
C     GERMANY
C     E-MAIL: WIEDER@HRZ.UNI-KASSEL.DE
C     HTTP://WWW.UNI-KASSEL.DE/~WIEDER
C
C     PURPOSE:
C
C     SUBROUTINE *** HANKEL *** CALCULATES THE HANKEL TRANSFROM
C     OF A REAL-VALUED FUNCTION *** FUN ***.
C
C                     ( infty
C     HT(FUN,NU,XI) =  |        (XI * X)**(1/2) JNUX(XI * X) FUN(X) dX
C                   0 )
C
C     HT = HANKEL TRANSFORM OF FUNCTION *** FUN ***.
C
C     THE FORM OF ** FUN *** IS: Y = FUN(X).
C
C     ------------------------------------------------------------------
C
C     INPUT VARIABLES:
C
C     ------------------------------------------------------------------
C
C     *** IDATIN ***:
C
C     FOR *** IDAT=1 *** THE USER HAS TO PROVIDE A
C     SUBROUTINE *** DHFUNC ***
C     WHICH GIVES THE FUNCTION *** FUN ***: Y = FUN(X).
C
C     THE FORM OF SUBROUTINE *** DHFUNC *** IS:
C
C        SUBROUTINE DHFUNC(NU,X,Y,XLIMIT,IERR)
C        IMPLICIT NONE
C        INTEGER IERR,LUIN,LUOUT,LUERR
C        DOUBLE PRECISION X,Y,NU,XASYMP,XLIMIT
C        COMMON /INOUTE/ LUIN,LUOUT,LUERR
C
C        ###############################################################
C
C        INPUT VARIABLES: NU, X
C        OUPUT VARIABLES: Y,XLIMIT,IERR
C
C        ###############################################################
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
C     ------------------------------------------------------------------
C
C     FOR *** IDAT=2 *** THE USER HAS TO PROVIDE AN INPUT FILE
C     WHICH CONTAINS *** FUN *** AS A TABLE OF (X,Y)-PAIRS.
C
C     THE FORM OF THE INPUT FILE IS:
C
C     X_1 Y_1
C     ...
C     X_N Y_N
C     ...
C     X_NEND Y_NEND
C
C     THESE DATA WILL BE READ BY A CALL TO SUBROUTINE *** FUNK2 ***
C     PRIOR TO THE CALL TO SUBROUTINE *** HANKEL *** !!! !!! !!!
C
C     ------------------------------------------------------------------
C
C     FOR *** IDAT=3 *** THE USER HAS TO PASS FUNCTION *** FUN ***
C     IN FORM OF TWO ARRAYS *** XDAT(N) ***, **** YDAT(N) ***
C     VIA A CALL TO SUBROUTINE *** FUNK3 ***
C     PRIOR TO THE CALL TO SUBROUTINE *** HANKEL *** !!! !!! !!!
C
C     THE ARRAY *** XDAT(N) ***, *** YDAT(N) *** CONTAIN THE VALUES
C     DATA X_1,Y_1, ..., X_N,Y_N, ..., X_NEND,Y_NEND.
C
C     *** XI ***:
C
C     POSITION AT WHICH THE HANKEL TRANSFORM HAS TO BE EVALUATED.
C     *** XI *** IS A POSITIVE REAL VARIABLE (IN DOUBLE PRECISION).
C
C     *** NU ***:
C
C     *** NU *** IS THE ORDER OF THE HANKEL TRANSFORM.
C     *** NU *** MAY BE A POSITIVE OR NEGATIVE REAL VALUE
C     IN DOUBLE PRECISION.
C     HOWEVER, A MAXIMUM MAGNITUDE HAS BEEN SET FOR *** NU ***
C     WITHIN SUBROUTINE *** HANKEL ***. YOU MAY ALTER THIS SETTING.
C
C     ------------------------------------------------------------------
C
C     *** IVERBO ***:
C
C     A SWITCH BY WHICH THE USER DECIDES HOW MUCH INTERMEDIATE
C     OUTPUT WILL BE GENERATED BY SUBROUTINE *** HANKEL ***.
C     PLEASE SET *** IVERBO *** EQUAL TO 0, 1 OR 2.
C
C     ------------------------------------------------------------------
C
C     *** LUOUT *** = LOGICAL NUMBER OF OUTPUT FILE
C     *** LUERR *** = LOGICAL NUMBER OF ERROR LOG FILE
C
C     ------------------------------------------------------------------
C
C     OUTPUT VARIBLES:
C
C     ------------------------------------------------------------------
C
C     *** IERR1 ***:
C
C     *** IERR1 *** IS AN ERROR FLAG.
C     *** IRERR1 *** = 0 --> NO ERROR HAS BEEN DETECTED,
C                            SUCCESSFULL CALCULATION.
C     *** IERR1 *** NOT EQUAL TO 0 --> AN ARROR OCCURED,
C                                      PLEASE CONSULT THE OUTPUT FILES!
C
C     ------------------------------------------------------------------
C
C     *** HT ***:
C
C     *** HT *** IS THE VALUE OF THE HANKEL TRANSFORM.
C     *** HT *** IS A POSITIVE OR NEGATIVE REAL VARIABLE
C     (IN DOUBLE PRECISION).
C
C     ------------------------------------------------------------------
C
C     LITERATURE:
C
C     HANDBUCH DER PHYSIK, S. FLUEGGE (HRSG.),
C     SPRINGER VERLAG, BERLIN, 1955, BAND II
C     I.N. SNEDDON: FUNCTIONAL ANALYSIS, PP. 198-348,
C     SECTION IV: THE HANKEL TRANSFORM, PP. 298-307.
C
C     A. ERDELYI, W. MAGNUS, F. OBERHETTINGER, F.G. TRICOMI:
C     TABLES OF INTEGRAL TRANSFORMS, VOL. II,
C     McGRAW-HILL BOOK COMPANY, NEW YORK, 1954.
C     CHAPTER VIII: HANKEL TRANSFORMS.
C
C     ------------------------------------------------------------------
C
C     SUCCESSFULLY COMPILED BY:
C
C     1.) AIX FORTRAN XLF (UNIX)
C     2.) NAG FTN90 (DOS)
C     3.) F90 NAG (UNIX)
C     4.) DIGITAL VISUAL FORTRAN (WINDOWS NT)
C
C     ------------------------------------------------------------------
C
C     CALLING SEQUENCE:
C
C     CALLS SUBROUTINE *** FUNCHE ***.
C     CALLS SUBROUTINE *** HNKLNM ***.
C     CALLS SUBROUTINE *** DHINFO ***.
C
C     ------------------------------------------------------------------
C
C     MOST IMPORTANT VARIABLES:
C
C     XI = PARAMETER OF THE HANKEL TRANSFORM
C        = POSITION AT WHICH THE HANKEL TRANSFORM *** HT ***
C          OF FUNCTION *** FUN *** HAS TO BE EVALUATED.
C
C     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C     XI IS AN INPUT VARIABLE TO SUBROUTINE *** HANKEL ***.
C     XI IS A POSITIVE, REAL-VALUED VARIABLE!
C     XI > 0
C     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C
C     NU = ORDER OF HANKEL TRANSFORM
C        = ORDER OF BESSEL FUNCTION *** DHJNUX ***.
C
C     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C     NU IS AN INPUT VARIABLE TO SUBROUTINE *** HANKEL ***.
C     NU IS A REAL-VALUED VARIABLE!
C     NU MAY BE POSITIVE, ZERO, OR NEGATIVE.
C     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C
C     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C     HT IS AN OUTPUT VARIABLE FROM SUBROUTINE *** HANKEL ***.
C     HT IS A REAL-VALUED VARIBALE!
C     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C
C     IDAT EQUAL TO 1 --> FUNCTION *** FUN *** IS GIVEN ANALYTICALLY.
C     IDAT EQUAL TO 2 --> FUNCTION *** FUN *** IS GIVEN AS TABULATED
C                         DATA WITHIN A FILE.
C     IDAT EQUAL TO 3 --> FUNCTION *** FUN *** IS PASSED AS TABULATED
C                         DATA VIA CALL.
C
C     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C     IDAT IS AN INPUT VARIABLE TO SUBROUTINE *** HANKEL ***.
C     IDAT IS A INTEGER-VALUED VARIABLE!
C     IDAT MAY TAKE ON THE VALUES 1, 2 OR 3.
C     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C
C     IERR1 = ERROR FLAG
C     IERR1 = 0 --> NO ERROR OCCURED DURING CALCULATION.
C
C     NMAX = MAXIMUM NUMBER OF ALLOWED FUNCTION VALUES IN CASE OF
C            A TABULATED FUNCTION (*** IDAT *** = 2).
C
C     ##################################################################
C
C      IMPLICIT NONE
C
C     .. Parameters ..
      DOUBLE PRECISION ZERO
      PARAMETER (ZERO=0.0D0)
C     ..
C     .. Scalar Arguments ..
      DOUBLE PRECISION HT,NU,XI
      INTEGER IERR1,IVERBO
C     ..
C     .. Scalars in Common ..
      DOUBLE PRECISION SLTN,XISLTN
      INTEGER IDAT,IFUNC,IREAD,LUERR,LUIN,LUOUT
      CHARACTER*80 FILEER,FILEIN,FILEOU
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION DPRCSN,NUMAX,NUMIN,XIMAX,XIMIN
      INTEGER IERR7,IERRH1,IHPKNS,IMOD
      CHARACTER*40 DATE
C     ..
C     .. External Functions ..
      DOUBLE PRECISION D1MACH
      EXTERNAL D1MACH
C     ..
C     .. External Subroutines ..
      EXTERNAL FUNCHE,HNKLNM
C     ..
C     .. Common blocks ..
      COMMON /INOUTE/LUIN,LUOUT,LUERR
      COMMON /LOGIK/IDAT,IREAD
      COMMON /NAMEN/FILEIN,FILEOU,FILEER
      COMMON /SOLUTI/XISLTN,SLTN,IFUNC
      integer i1mach, nout
C     ..
C     .. Save statement ..
      SAVE IHPKNS
C     ..
      nout = i1mach(2)
      DATE = ' FEBRUARY 1999'
C
C    -------------------------------------------------------------------
C
C     SET SOME PARAMETERS:
C
      IERR1 = 0
C     CODE MODIFIED BY DR. T.R. HOPKINS, 28.04.1998:
      IHPKNS = 0
C
C     XISLTN IS NECCESSARY FOR DEBUGGING PURPOSES ONLY.
      XISLTN = XI
C
      XIMIN = ZERO
      XIMAX = D1MACH(2)
C
C     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C     THE CHOICE OF *** NUMIN *** AND *** NUMAX ***
C     IS SOMEWHAT ARBITRARY. YOU MAY MEET ANOTHER CHOICE.
      NUMIN = -100.0D0
      NUMAX = 100.0D0
C     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C
      IF (IVERBO.NE.0) WRITE(NOUT,FMT=9000)
      IF (IVERBO.GT.1) WRITE (LUERR,FMT=9000)
      IF (IVERBO.NE.0) WRITE(NOUT,FMT=9010) DATE
      IF (IVERBO.GT.1) WRITE (LUERR,FMT=9010) DATE
C
C    -------------------------------------------------------------------
C
C     CHECK THE INPUT VALUES:
C
C     IS *** XI *** POSITIVE?
      IF (XI.LE.ZERO) THEN
          WRITE(NOUT,FMT=9020)
          WRITE (LUERR,FMT=9020)
          WRITE(NOUT,FMT=9030) XI
          WRITE (LUERR,FMT=9030) XI
          IERR1 = 1
          WRITE(NOUT,FMT=9130)
C     AN ERROR OCCURED, DO NOT CONTINUE, GO TO STOP!
          GO TO 20

      END IF
C
C     CHECK ABSOLUT VALUE OF *** XI ***
      IF (XI.LT.XIMIN .OR. XI.GT.XIMAX) THEN
          WRITE(NOUT,FMT=9020)
          WRITE (LUERR,FMT=9020)
          WRITE(NOUT,FMT=9040) XI
          WRITE (LUERR,FMT=9040) XI
          WRITE (LUERR,FMT=9050) XIMIN,XIMAX
          IERR1 = 1
          WRITE(NOUT,FMT=9130)
C     AN ERROR OCCURED, DO NOT CONTINUE, GO TO STOP!
          GO TO 20

      END IF
C
C     CHECK ABSOLUT VALUE OF *** NU ***
      IF (NU.LT.NUMIN .OR. NU.GT.NUMAX) THEN
          WRITE(NOUT,FMT=9020)
          WRITE (LUERR,FMT=9020)
          WRITE(NOUT,FMT=9060) NU
          WRITE (LUERR,FMT=9060) NU
          WRITE (LUERR,FMT=9070) NUMIN,NUMAX
          IERR1 = 1
          WRITE(NOUT,FMT=9130)
C     AN ERROR OCCURED, DO NOT CONTINUE, GO TO STOP!
          GO TO 20

      END IF
C
      IF (IDAT.EQ.1) THEN
          IF (IVERBO.NE.0) WRITE(NOUT,FMT=9080)
          IF (IVERBO.GT.1) WRITE (LUERR,FMT=9080)
C
      ELSE IF (IDAT.EQ.2) THEN
          IF (IVERBO.NE.0) WRITE(NOUT,FMT=9090)
          IF (IVERBO.GT.1) WRITE (LUERR,FMT=9090)
C
C     THIS IS FOR INTERNAL LOGIC, NEVER CHANGE THE FOLLOWING LINE:
          IFUNC = 0
C
      ELSE IF (IDAT.EQ.3) THEN
          IF (IVERBO.NE.0) WRITE(NOUT,FMT=9100)
          IF (IVERBO.GT.1) WRITE (LUERR,FMT=9100)
C
C     THIS IS FOR INTERNAL LOGIC, NEVER CHANGE THE FOLLOWING LINE:
          IFUNC = 0
C
      ELSE
          WRITE(NOUT,FMT=9110) IDAT
          WRITE (LUERR,FMT=9110) IDAT
          WRITE(NOUT,FMT=9130)
C     AN ERROR OCCURED, DO NOT CONTINUE, GO TO STOP!
          GO TO 20

      END IF
C
C    -------------------------------------------------------------------
C
C     CHECK THE FUNCTION *** DHFUNC ***:
C
      CALL FUNCHE(IERR7)
C
      IF (IERR7.NE.0) THEN
          WRITE(NOUT,FMT=9120) FILEER
          WRITE (LUERR,FMT=9120) FILEER
          WRITE(NOUT,FMT=9130)
          WRITE(NOUT,FMT=9140) FILEER
C
C         CALL DHINFO
C
          STOP

      END IF
C
C    -------------------------------------------------------------------
C
C     BEGIN OF CALCULATION!
C
      IF (XI.EQ.ZERO) THEN
          HT = ZERO
          GO TO 10

      END IF
C
      IMOD = 1
C
      IF (IMOD.EQ.1) THEN
C
C     !!! BEGIN OF MODE 1: !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C
C     MODE 1: HANKEL TRANSFORM BY EXPLICIT CALCULATION OF THE
C             TRANSFORMATION FORMULA.
C             I.E. HANKEL TRANSFORM BY INDFINITE INTEGRATION.
C
          IF (IVERBO.NE.0) WRITE(NOUT,FMT=9150)
          IF (IVERBO.GT.1) WRITE (LUERR,FMT=9160)
C
          CALL HNKLNM(XI,NU,IVERBO,IERRH1,HT)
C
          IF (IVERBO.GT.1) WRITE (LUERR,FMT=9170) XI,NU,HT,IERRH1
C
          IF (IERRH1.NE.0) THEN
              WRITE(NOUT,FMT=9180)
              WRITE (LUERR,FMT=9180)
              IERR1 = 1
              GO TO 20

          END IF
C
C     !!! END OF MODE 1: !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C
      ELSE
          WRITE (LUERR,FMT=9190) IMOD
          IERR1 = -1
          GO TO 20
C
C     END IF IF-CLAUSE IMOD
      END IF
C
C     END OF CALCULATION!
C
C     ------------------------------------------------------------------
C
C     WRITE OUT RESULTS:
C
   10 IF (IHPKNS.NE.1) THEN
C
C      *** AMENDMENT SUGGESTED BY DR. T.R. HOPKINS, 15.04.1998 ***
C     WRITE OUT HEAD LINE ONLY ONCE:
          IHPKNS = 1
C
          IF (IFUNC.EQ.0 .OR. IDAT.EQ.2 .OR. IDAT.EQ.3) THEN
              IF (IVERBO.GT.0) WRITE (LUOUT,FMT=9200)

          ELSE IF (IFUNC.NE.0) THEN
              IF (IVERBO.GT.0) WRITE (LUOUT,FMT=9210)
          END IF
C
C     END IF IF_CLAUSE *** IHPKNS ***
      END IF
C
      IF (IFUNC.EQ.0 .OR. IDAT.EQ.2 .OR. IDAT.EQ.3) THEN
C
C     IFUNC = 0 --> A USER SUPPLIED FUNCTION HAS BEEN TRANSFORMED:
          IF (IVERBO.GT.0) WRITE(NOUT,FMT=9220) XI,NU,HT
          IF (IVERBO.GT.0) WRITE (LUERR,FMT=9220) XI,NU,HT
          IF (IVERBO.GT.0) WRITE (LUOUT,FMT=9230) XI,HT
C
      ELSE IF (IFUNC.NE.0) THEN
C
C     IFUNC NOT 0 --> A TEST FUNCTION FROM WITHIN FUNCTION *** DHFUNC **
C     HAS BEEN TRANSFORMED:
          WRITE(NOUT,FMT=9240) XISLTN,SLTN,HT
          WRITE (LUOUT,FMT=9250) XISLTN,SLTN,HT
          IF (IVERBO.GT.0) THEN
C     CALCULATE THE MACHINE PRECISION *** DPRSCN ***
	      DPRCSN = D1mach(4)
              IF (SLTN.NE.0.0D0) THEN
                  WRITE (LUERR,FMT=9260) XISLTN,SLTN,HT,HT - SLTN,
     +              (HT-SLTN)/SLTN, (HT-SLTN)/DPRCSN,
     +              (HT-SLTN)/SLTN/DPRCSN

              ELSE
                  WRITE (LUERR,FMT=9270) XISLTN,SLTN,HT
              END IF

          END IF
C
C     END IF-CLAUSE *** IFUNC ***
      END IF
C
   20 IF (IVERBO.NE.0) WRITE(NOUT,FMT=9280) FILEER
      IF (IVERBO.NE.0) WRITE(NOUT,FMT=9290) FILEOU
C
      IF (IVERBO.NE.0) WRITE(NOUT,FMT=9300) IERR1
      IF (IVERBO.GT.1) WRITE (LUERR,FMT=9300) IERR1
      GO TO 30
C
C     ------------------------------------------------------------------
C
C     ERROR HANDLING:
C
   30 IF (IERR1.EQ.0) THEN
          IF (IVERBO.NE.0) WRITE(NOUT,FMT=9310)
          IF (IVERBO.GT.1) WRITE (LUERR,FMT=9310)
C
      ELSE IF (IERR1.NE.0 .AND. IERR1.NE.-1) THEN
          WRITE(NOUT,FMT=9320)
          WRITE (LUERR,FMT=9320)
C
      ELSE IF (IERR1.EQ.-1) THEN
C
C         CALL DHINFO
C
          WRITE(NOUT,FMT=9330)
          CLOSE (LUOUT)
          CLOSE (LUERR)
          STOP
C
C     ENDE DER ABFRAGE AUF *** IERR1 ***.
      END IF
C
      IF (IVERBO.GT.1) WRITE(NOUT,FMT=9340)
      IF (IVERBO.GT.1) WRITE (LUERR,FMT=9340)
C
C    -------------------------------------------------------------------
C
C     *** SUBROUTINE DHINFO *** PROVIDES FURTHER INFORMATIONS:
C
C     IF (IVERBO.GT.0) CALL DHINFO
C
C    -------------------------------------------------------------------
C
      RETURN

 9000 FORMAT (/,/,' ####################################',/,
     +       ' MESSAGE HANKEL: BEGIN OF CALCULATION',/,
     +       ' ####################################')
 9010 FORMAT (/,' *** PROGRAM HANKEL *** ',/,
     +       ' WRITTEN BY THOMAS WIEDER',/,
     +       ' E-MAIL: WIEDER@HRZ.UNI-KASSEL.DE',/,
     +       ' DATE OF FIRST VERSION: A.D. 1998',/,
     +       ' DATE OF LAST CHANGE:',A40)
 9020 FORMAT (/,' MESSAGE HANKEL: WRONG INPUT VALUE!')
 9030 FORMAT (/,' MESSAGE HANKEL:',/,
     +       ' PARAMETER *** XI *** IS EQUAL TO ZERO OR NEGATIVE!',/,
     +       ' *** XI *** =',D17.9)
 9040 FORMAT (/,' MESSAGE HANKEL:',/,
     +       ' PARAMETER *** XI *** OUTSIDE ALLOWED RANGE!',/,
     +       ' *** XI *** =',D17.9)
 9050 FORMAT (/,' ALLOWED INTERVAL [XI_MIN,XI_MAX]:',/,'[',D17.9,',',
     +       D17.9,']')
 9060 FORMAT (/,' MESSAGE HANKEL:',/,
     +       ' PARAMETER *** NU *** OUTSIDE ALLOWED RANGE!',/,
     +       ' *** NU *** =',D17.9)
 9070 FORMAT (/,' ALLOWED INTERVAL [NU_MIN,NU_MAX]:',/,' [',D17.9,',',
     +       D17.9,']')
 9080 FORMAT (/,' MESSAGE HANKEL:',/,
     +       ' WILL TRY TO TRANSFORM *** FUNCTION FUNC ***!')
 9090 FORMAT (/,' MESSAGE HANKEL:',/,
     +       ' WILL TRY TO TRANSFORM TABULATED DATA FROM INPUT FILE!')
 9100 FORMAT (/,' MESSAGE HANKEL:',/,
     +' WILL TRY TO TRANSFORM TABULATED DATA TRANSMITTED BY CALL TO HANK
     +EL!')
 9110 FORMAT (/,' MESSAGE HANKEL',/,
     +       ' WRONG VALUE FOR CALL PARAMETER *** IDAT *** !',/,
     +       ' IDAT = ',I6)
 9120 FORMAT (/,' MESSAGE HANKEL:',/,
     +' A TEST OF YOUR FUNCTION *** DHFUNC       *** GAVE AN ERRONEOUS R
     +ESULT!',/,' PLEASE HAVE A LOOK INTO FILE',/,A80)
 9130 FORMAT (/,' MESSAGE HANKEL: NO TRANSFORMATION CALCULATED! STOP!')
 9140 FORMAT (/,' MESSAGE HANKEL: HAVE A LOOK INTO FILE:',/,A80)
 9150 FORMAT (/,
     +' MESSAGE HANKEL: WILL TRY NUMERICAL HANKEL TRANSFORM...PLEASE WAI
     +T!')
 9160 FORMAT (/,
     +      ' MESSAGE HANKEL: WILL CALL SUBROUTINE *** HNKLNM **    * !'
     +       )
 9170 FORMAT (/,
     +    ' MESSAGE HANKEL: XI          NU          HT      ERROR FLAG:'
     +       ,/,D17.9,3X,D17.9,3X,D17.9,3X,I6)
 9180 FORMAT (/,' MESSAGE HANKEL: ERROR DURING TRANSFORMATION!',/,
     +       ' ERROR IN *** SUBROUTINE HNKLNM *** !')
 9190 FORMAT (/,' MESSAGE HANKEL:',/,' WRONG CALCULATION MODUS CHOOSEN!'
     +       ,/,' IMOD = ',I6)
 9200 FORMAT ('# XI                    HANKEL_TRANSFORM(XI)')
 9210 FORMAT (
     +'#  XI                    EXACT SOLUTION             NUMERICAL SOL
     +UTION')
 9220 FORMAT (/,' MESSAGE HANKEL:',/,
     +'     ARGUMENT XI              ORDER NU            NUMERICAL HANKE
     +L TRANSFORM HT:',/,1X,D17.9,8X,D17.9,8X,D17.9)
 9230 FORMAT (E17.9,5X,E17.9)
 9240 FORMAT (/,
     +  ' ARGUMENT XI,        ANALYTICAL TRANSFORM, NUMERICAL SOLUTION:'
     +       ,/,D17.9,3X,D17.9,3X,D17.9)
 9250 FORMAT (E17.9,5X,E17.9,5X,E19.7)
 9260 FORMAT (1X,D17.9,1X,D17.9,1X,D17.9,1X,D17.9,1X,D17.9,1X,D17.9,1X,
     +       D17.9)
 9270 FORMAT (3X,D17.9,3X,D17.9,3X,D17.9)
 9280 FORMAT (/,' MESSAGE HANKEL: FOR ERROR MESSAGES LOOK INTO FILE:',/,
     +       A80)
 9290 FORMAT (' MESSAGE HANKEL: FOR RESULTS LOOK INTO FILE:',/,A80)
 9300 FORMAT (' IERR1 =',I6)
 9310 FORMAT (/,' MESSAGE HANKEL: SUCCESSFUL END OF CALCULATION!',/,
     +     ' NO ERRORS HAVE BEEN DETECTED DURING CALCULATION! GOOD BYE!'
     +       )
 9320 FORMAT (/,' MESSAGE HANKEL:',/,' UNSUCCESSFUL END OF CALCULATION!'
     +       ,/,' ERRORS HAVE BEEN DETECTED DURING CALCULATION!',/,
     +       ' GOOD BYE!')
 9330 FORMAT (/,/,' MESSAGE HANKEL:',/,
     +       ' END OF CALCULATION BECAUSE OF UNRCOVERABLE ERROR!',/,
     +       ' GOOD BYE!')
 9340 FORMAT (/,' ####################################',/,
     +       ' MESSAGE HANKEL: END OF CALCULATION',/,
     +       ' ####################################',/,/)
      END
C     END OF SUBROUTINE *** HANKEL ***
C
C     BEGIN OF SUBROUTINE *** HNKLNM ***
      SUBROUTINE HNKLNM(XIIN,NUIN,IVERBO,IERR2,HT)
C
C     ##################################################################
C
C     FIRST VERSION:    21.01.1998
C     PREVIOUS VERSION: 10.12.1998
C     LATEST VERSION:   05.02.1999
C
C     PROGRAMMED BY:
C     THOMAS WIEDER
C     E-MAIL: WIEDER@HRZ.UNI-KASSEL.DE
C     HTTP://WWW.UNI-KASSEL.DE/~WIEDER
C
C     PURPOSE:
C
C     SUBROUTINE *** HNKLNM *** CALCULATES THE HANKEL TRANSFROM
C     OF A REAL-VALUED FUNCTION *** DHFUNC ***
C     BY EXPLICIT NUMERICAL EVALUATION OF THE TRANSFORMATION FORMULA.
C
C     LITERATURE:
C
C     A. ERDELYI, W. MAGNUS, F. OBERHETTINGER, F.G. TRICOMI:
C     TABLES OF INTEGRAL TRANSFORMS, VOL. II,
C     McGRAW-HILL BOOK COMPANY, NEW YORK, 1954.
C     CHAPTER VIII: HANKEL TRANSFORMS.
C
C     MOST IMPORTANT VARIABLES:
C
C     XI = PARAMETER OF THE HANKEL TRANSFORM.
C     THE HANKEL TRANSFORM OF *** FUN *** WILL BE CALCULATED
C     FOR ARGUMENT *** XI ***.
C     *** FUN *** DEPENDS ON VARIABLE *** X ***.
C     HT = HANKEL TRANSFORM OF FUNCTION *** FUN ***.
C
C     CALLING SEQUENCE:
C
C     CALLED BY SUBROUTINE *** HANKEL ***.
C     CALLS FUNCTION *** FUNINT ***.
C     CALLS SUBROUTINE *** INTHP ***.
C     CALLS SUBROUTINE *** OSCINT ***.
C
C     ##################################################################
C
C      IMPLICIT NONE
C
C     INTEGER VARIABLES FOR SUBROUTINE *** INTHP ***
C
C     INTEGER VARIABLES FOR SUBROUTINE *** OSCINT ***
C     PARAMETERS FOR SUBROUTINE *** OSCINT ***
C     *** 19.02.1998 *** SUCCESSFUL VALUE FOR *** NDIM1 *** = 1000
C     *** 19.02.1998 *** SUCCESSFUL VALUE FOR *** NDIM2 *** = 28
C
C
C
C
C     REAL VARIABLES FOR SUBROUTINE *** INTHP ***
C
C     REAL VARAIBLES FOR SUBROUTINE *** OSCINT ***
C
C     VARIABLES FOR TESTING SUBROUTINE *** INTHP ***
C     1 ,TFUNC
C
C
C     COMMON BLOCK FOR SUBROUTINE *** INTHP ***
C
C
C     EXTERNAL FUNCTION FOR SUBROUTINE *** INTHP ***
C      EXTERNAL TFUNC
C
C     EXTERNAL FUNCTION FOR SUBROUTINE *** OSCINT ***
C
C    -------------------------------------------------------------------
C
C     SAVE INPUT
C
C     .. Parameters ..
      INTEGER NDIM1,NDIM2
      PARAMETER (NDIM1=1000,NDIM2=28)
      INTEGER NMAX
      PARAMETER (NMAX=1000)
C     ..
C     .. Scalar Arguments ..
      DOUBLE PRECISION HT,NUIN,XIIN
      INTEGER IERR2,IVERBO
C     ..
C     .. Scalars in Common ..
      DOUBLE PRECISION EXACT,NU,XDATMA,XDATMI,XI,XLIMIT,XSPLIT,YDATMA,
     +                 YDATMI
      INTEGER IDAT,IERR3,IINT,IREAD,LUERR,LUIN,LUOUT,NEND
C     ..
C     .. Arrays in Common ..
      DOUBLE PRECISION XDAT(NMAX),YDAT(NMAX)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION ALIMIT,BLIMIT,DCLASS,DUMMY,EPS,PCLASS,PERIOD,PI,
     +                 RESULT,RFIRST,XASYMP
      INTEGER IERR4,IERRME,INF,INTPAR,INTSUB,MCALL,N1,N2,NQUAD
C     ..
C     .. Local Arrays ..
      DOUBLE PRECISION ABSCIS(36),QLIST(0:NDIM1),SAVPER(0:256),
     +                 WEIGHT(36),WORK(0:NDIM1,0:NDIM2)
      INTEGER ISTATE(6)
C     ..
C     .. External Functions ..
      DOUBLE PRECISION FUNINT,GPER,HFUN
      EXTERNAL FUNINT,GPER,HFUN
C     ..
C     .. External Subroutines ..
      EXTERNAL ERRMSS,FUNK1,FUNK2,FUNK3,G5AND9,INTHP,OSCINT
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC DACOS
C     ..
C     .. Common blocks ..
      COMMON /ARGUME/XI,NU,XLIMIT
      COMMON /DATEN/XDAT,YDAT,XDATMI,XDATMA,YDATMI,YDATMA,NEND
      COMMON /FHNKL1/IERR3
      COMMON /INOUTE/LUIN,LUOUT,LUERR
      COMMON /LOGIK/IDAT,IREAD
      COMMON /LOSUNG/XSPLIT,EXACT,IINT
      integer i1mach, nout
      nout = i1mach(2)
C     ..
      XI = XIIN
      NU = NUIN
C
C    -------------------------------------------------------------------
C
C     SET SOME PARAMETERS:
C
      IERR2 = 0
      IERR3 = 0
      IERR4 = 0
      HT = 0.0D0
      RESULT = 0.0D0
      PI = 2.0D0*DACOS(0.0D0)
      DUMMY = 0.0D0
C
C     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C     *** 24.02.1998 *** A GOOD VALUE FOR *** XASYMP *** = 100.0
      XASYMP = 100.0D0
C     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C
C     CHOOSE A TEST FUNCTION WITHIN SUBROUTINE *** TINTHP ***
C     BY SETTING *** IINT ***:
      IINT = 9
C
C    -------------------------------------------------------------------
C
C     BEGIN THE CALCULATION:
C
C     WE SPLIT THE INFINITE INTEGRAL INTO TWO PARTS.
C     THE FIRST PART IS A FINITE INTEGRAL FROM 0 TO *** XSPLIT ***.
C     THE SECOND PART IS AN INFINITE INTEGRAL FROM *** XSPLIT ***
C     TO INFINITY.
C
C     WE CHOOSE *** XSPLIT *** IN A WAY THAT THE OSCILLATING PART
C     IN THE INFINITE INTEGRAL HAS A CONSTANT PERIOD OF OSCILLATION!!!
C
      XSPLIT = XASYMP
C
C     IF THE FUNCTION *** FUN *** IS EQUAL TO ZERO FOR *** X ***
C     GREATER THAN *** XLIMIT ***
C     OR THE FUNCTION *** FUN *** MAY BETTER BE INTEGRATED FROM
C     *** XLIMIT *** TO INFINITY THAN FROM *** XASYMP ***
C     TO INFINITY,
C     THEN WE SET *** XSPLIT = XLIMIT ***.
C
C     FIRST, CALL SUBROUTINE *** FUNK ** TO GET THE VALUE OF
C     *** XLIMIT ***.
C     *** XLIMIT *** APPEARS IN COMMON BLOCK *** /ARGUMENTE/ ***.
C
      XLIMIT = -1.0D0
C
      IF (IDAT.EQ.1) THEN
          CALL FUNK1(DUMMY,DUMMY,IERR4)

      ELSE IF (IDAT.EQ.2) THEN
          CALL FUNK2(DUMMY,DUMMY,IERR4)

      ELSE IF (IDAT.EQ.3) THEN
          CALL FUNK3(DUMMY,DUMMY,XDAT,YDAT,NEND,IERR4)
      END IF
C
      IF (IERR4.NE.0) THEN
          WRITE (LUERR,FMT=9000)
      END IF
C
      IF (XLIMIT.LT.0.0D0) THEN
          IF (IVERBO.GT.1) THEN
              WRITE(NOUT,FMT=9010) XSPLIT
              WRITE (LUERR,FMT=9010) XSPLIT
          END IF

      END IF
C
C     NOW WE KNOW *** XLIMIT ***.
C
      IF (XLIMIT.GT.0.0D0) THEN
          IF (XLIMIT.GT.XASYMP) THEN
              IF (IVERBO.GT.1) WRITE (LUERR,FMT=9030) XLIMIT
              IF (IVERBO.GT.1) WRITE (LUERR,FMT=9020) XLIMIT,XASYMP,
     +            XASYMP

          ELSE
              XSPLIT = XLIMIT
              IF (IVERBO.GT.1) WRITE (LUERR,FMT=9030) XLIMIT
C     END OF IF-CLAUSE *** XASYMP ***
          END IF
C     END OF IF-CLAUSE *** XLIMIT ***
      END IF
C
      IF (IVERBO.GT.1) WRITE (LUERR,FMT=9040) XSPLIT
C
C     *** INPART *** = 1 --> DO THE FINITE INTEGRAL!
      INTPAR = 1
C
   10 IF (INTPAR.LE.2) THEN
          IF (INTPAR.EQ.1) THEN
C     DO THE FINITE INTEGRAL:
C     CHOOSE INTEGRATION SUBROUTINE BY SETTING *** INTSUB ***:
              INTSUB = 1

          ELSE IF (INTPAR.EQ.2) THEN
C     DO THE INFINITE INTEGRAL:
              INTSUB = 2

          ELSE
              WRITE (LUERR,FMT=9050) INTPAR
              STOP
C     END OF IF-CLAUSE *** INTPAR ***.
          END IF

      ELSE
          WRITE (LUERR,FMT=9050) INTPAR
C     END OF IF-CLAUSE *** INTPAR *** .LT. 2
      END IF
C
      IF (INTSUB.EQ.1) THEN
C     *** SBROUTINE INTHP *** CALCULATES A FINITE OR
C     AN INFINITE INTEGRAL.
          ALIMIT = 0.0D0
          BLIMIT = XSPLIT
          DCLASS = 10.0D0
          MCALL = 200000
          PCLASS = 0.0D0
C     01.10.1998 EPS=1.0D-06 --> 1.0D-09
          EPS = 1.0D-09
C     *** INF *** = 3 --> INFINITE OSZILLATING TAIL
C     *** INF *** = 4 --> FINITE INTEGRAL
          INF = 4
C
          CALL INTHP(ALIMIT,BLIMIT,DCLASS,FUNINT,MCALL,PCLASS,EPS,INF,
     +               RESULT)
C
          IF (IVERBO.GT.1) WRITE (LUERR,FMT=9060) RESULT
C
          IF (INF.GT.3) THEN
              WRITE (LUERR,FMT=9070) INF
              IERR2 = 1
          END IF
C
      ELSE IF (INTSUB.EQ.2) THEN
C     SUBROUTINE *** OSCINT *** CALCULATES AN INFINITE, OSCILLATING
C     INTEGRAL.
C
C     THE FOLLOWING FEW LINES ARE TAKEN FROM SOURCE CODE OF ACM 639:
C
C     SET THE INPUT PARAMETERS FOR OSCINT.
C
C     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C     THE PERIOD OF THE ASYMPTOTIC EXPANSION OF THE BESSEL FUNCTION:
C     THIS EXPANSION IS THE COSINUS:
          PERIOD = 2.0D0*PI/XI
C     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C
C     ALIMIT (AZERO) IS THE LOWER LIMIT OF INTEGRATION
          ALIMIT = XSPLIT
C
C     RFIRST IS THE RIGHT-HAND ENDPOINT OF THE FIRST INTERVAL
          RFIRST = 0.0D0
C
C     NQUAD REPRESENTS THE TYPE OF QUADRATURE RULE USED
C           NQUAD = -5    5 PT GAUSS RULE
C           NQUAD = -9    9 PT GAUSS RULE
C           NQUAD = N,N>0 (N ODD) N POINT TRAPAZOIDAL RULE(N-1 PANELS)
          NQUAD = -9
C
C     SET EPSILON, THE REQUESTED TOLERANCE
C     *** T. WIEDER, 06.03.98 *** A GOOD VALUE FOR *** EPS *** = 1.0D-03
          EPS = 1.0D-03
C
C     NDIM1 AND NDIM2 ARE THE DIMENSIONS OF THE WORK TABLE. THEY
C     MUST BE LESS THAN OR EQUAL TO THE DIMENSIONS GIVEN TO WORK
C     IN THE DIMENSION STATEMENT
C
          DO 30 N1 = 0,NDIM1
              QLIST(N1) = 0.0D0
              DO 20 N2 = 0,NDIM2
                  WORK(N1,N2) = 0.0D0
   20         CONTINUE
   30     CONTINUE
          DO 40 N1 = 1,36
              SAVPER(N1) = 0.0D0
              WEIGHT(N1) = 0.0D0
              ABSCIS(N1) = 0.0D0
   40     CONTINUE
C
          CALL OSCINT(ALIMIT,PERIOD,RFIRST,EPS,NQUAD,NDIM1,NDIM2,G5AND9,
     +                HFUN,GPER,WORK,SAVPER,WEIGHT,ABSCIS,QLIST,RESULT,
     +                ISTATE)
C
          IF (IVERBO.GT.1) WRITE (LUERR,FMT=9080) RESULT
C
          IF (ISTATE(1).LT.0) THEN
C     AN ERROR HAS OCCURRED IN SUBROUTINE *** OSCINT ***.
              WRITE (LUERR,FMT=9090)
              WRITE (LUERR,FMT=9100) ISTATE(1),ISTATE(2),ISTATE(3),
     +          ISTATE(4),ISTATE(5),ISTATE(6)
              IERR2 = 1
C     RESET THE UNRELIABLE RESULT:
              RESULT = 0.0D0
          END IF
C
      ELSE IF (INTSUB.EQ.3) THEN
C     SUBROUTINE *** DQAINF *** CALCULATES AN INFINITE INTEGRAL.
          WRITE(NOUT,FMT=9110)
          WRITE (LUERR,FMT=*) 55
C
C      CALL DTESTI (IERR3,RESULTS)
C
          IERR3 = 1
          IF (IERR3.GT.0) THEN
              WRITE (LUERR,FMT=9120) IERR3
              IERR2 = 1
          END IF
C
C     END OF IF-CLAUSE *** INTSUB ***
      END IF
C
      HT = HT + RESULT
      INTPAR = INTPAR + 1
C     CHECK WHETHER THE INFINITE INTEGRAL STILL HAS TO BE EVALUATED:
      IF (INTPAR.LE.2) GO TO 10
C
      IF (IERR3.NE.0) THEN
C     AN ERROR IN FUNCTION *** FUNINT *** HAS OCCURED.
          WRITE(NOUT,FMT=9130)
          WRITE (LUERR,FMT=9130)
          WRITE (LUERR,FMT=9140) XI
          WRITE (LUERR,FMT=9150) NU
          IERRME = 10
C
          CALL ERRMSS(IERRME)
C
      END IF
C
      IF (IERR2.NE.0) THEN
C     AN ERROR DURING INTEGRATION HAS OCCURED!
          WRITE(NOUT,FMT=9160)
          WRITE (LUERR,FMT=9160)
          WRITE (LUERR,FMT=9140) XI
          WRITE (LUERR,FMT=9150) NU
          IERRME = 11
C
          CALL ERRMSS(IERRME)
C
      END IF
C
 9000 FORMAT (/,' MESSAGE HNKLNM:',/,
     +       ' ERROR MESSAGE FROM SUBROUTINE *  ** FUNK *** ',/,
     +       ' WHILE TRYING TO GET THE VALUE FOR *** XLIMIT *  ** !',/,
     +       ' HOWEVER, THE CALCULATION MAY STILL TURN OUT WELL...')
 9010 FORMAT (/,' MESSAGE HNKLNM:',/,
     +       ' *** XLIMIT *** LESS OR EQUAL TO   ZERO!',/,
     +' *** XLIMIT HAS BEEN REJECTED AS UPPER LIMIT *** XSPLIT *** FOR
     +THE FINITE INTEGRAL',/,' *** XSPLIT *** =',D17.9)
 9020 FORMAT (/,' MESSAGE HNKLNM: *** XLIMIT.GT.XASYMP *** !',/,
     +       ' *** XLIMIT *** =',D17.9,3X,' XASYMP =',D17.9,/,
     +       ' WILL NOT ACCEPT REQUESTED VALUE FOR *** XLIMIT ***!',/,
     + ' WILL USE AS UPPER LIMIT *** XSPLIT *** FOR THE FINITE INTEGRAL'
     +       ,/,' THE VALUE OF *** XASYMP *** =',D17.9)
 9030 FORMAT (/,' MESSAGE HNKLNM:',/,
     +' SUBROUTINE *** FUNK (OR DHFUNC)   *** REQUESTS TO TAKE THE FINIT
     +E',/,' INTEGRAL FROM 0 TO *** XLI  MIT *** =',D17.9)
 9040 FORMAT (/,' MESSAGE HNKLNM: ',/,
     +' WILL CALCULTATE THE FINITE       INTEGRAL FROM 0 TO *** XSPLIT *
     +** =',D17.9)
 9050 FORMAT (/,' MESSAGE HNKLNM: ERROR IN  INTERNAL LOGIC!',/,
     +       ' ***     INTPAR *** =',I6)
 9060 FORMAT (/,
     +     ' MESSAGE HNKLNM: *** INTHP *** RETURNS *** RESULT *    ** ='
     +       ,D17.9)
 9070 FORMAT (/,' MESSAGE HNKLNM:',/,
     + ' ERROR ON RETURN FROM INTEGRAT    ION SUBROUTINE *** INTHP ***!'
     +       ,/,' *** INF *** =',I6)
 9080 FORMAT (/,
     +    ' MESSAGE HNKLNM: *** OSCINT *** RETURNS *** RESULT     *** ='
     +       ,D17.9)
 9090 FORMAT (/,' MESSAGE HNKLNM:',/,
     +' ERROR ON RETURN FROM INTEGRAT    ION SUBROUTINE *** OSCINT ***!'
     +       ,/,' *** ISTATE *** =',I6)
 9100 FORMAT (/,' ISTATE(1)=',I6,2X,' ISTATE(2)=',I6,2X,' ISTATE(3)=',
     +       I6,2X,' ISTATE(4)=',I6,2X,' ISTATE(5)=',I6,2X,
     +       ' ISTATE(6)=',I6)
 9110 FORMAT (/,' MESSAGE HNKLNM:',/,
     +       ' SUBROUTINE *** DQAINF *** IS     NOT IMPLEMENTED!')
 9120 FORMAT (/,' MESSAGE HNKLNM:',/,
     +' ERROR ON RETURN FROM INTEGRAT    ION SUBROUTINE *** DQAINF ***!'
     +       ,/,' *** IFAIL *** =',I6)
 9130 FORMAT (/,' MESSAGE HNKLNM:',/,
     +       ' AN ERROR OCCURED DURING EVALU    ATION OF THE INTEGRAND!'
     +       )
 9140 FORMAT (/,' MESSAGE HNKLNM: PARAMETER XI = ',D19.7)
 9150 FORMAT (/,' MESSAGE HNKLNM: PARAMETER NU = ',D19.7)
 9160 FORMAT (/,
     +       ' MESSAGE HNKLNM: AN ERROR OCCURED DURING INTEGRATIO    N!'
     +       )
      END
C     END OF SUBROUTINE *** HNKLNM ***
C
C     BEGIN OF FUNCTION *** FUNINT ***
      DOUBLE PRECISION FUNCTION FUNINT(XIN)
C
C     FIRST VERSION:    21.01.1998
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
C     FUNCTION *** FUNINT *** CALCULATES THE INTEGRAND
C     (XI * X)**(1/2) JNUX(XI * X) FUN(X)
C     OF THE HANKEL TRANSFORMATION INTEGRAL.
C
C     CALLING SEQUENCE:
C     CALLED BY SUBROUTINE *** HNKLNM ***.
C     CALLS SUBROUTINE *** FUNK ***
C
C     MOST IMPORTANT VARIABLES:
C
C     X = INDEPENDENT REAL VARIABLE OF *** FUN ***.
C     Y = DEPENDENT REAL VARIABLE OF *** FUN ***.
C     Y = FUN(X)
C
C     XI = PARAMETER OF THE HANKEL TRANSFORM.
C     THE HANKEL TRANSFORM OF *** FUN *** WILL BE CALCULATED
C     FOR PARAMETER *** XI ***.
C
C     NU = ORDER OF BESSEL FUNCTION *** DHJNUX ***.
C
C     FUNINT = OUTPUT VALUE = INTEGRAND.
C
C     IERR3 = ERROR FLAG
C     IERR3 = 0 --> NO ERROR OCCURED DURING EVALUATION OF
C                  FUNCTION *** FUNINT ***.
C
C     ##################################################################
C
C      IMPLICIT NONE
C
C     SAVE INPUT
C
C     .. Parameters ..
      INTEGER NMAX
      PARAMETER (NMAX=1000)
C     ..
C     .. Scalar Arguments ..
      DOUBLE PRECISION XIN
C     ..
C     .. Scalars in Common ..
      DOUBLE PRECISION NU,XDATMA,XDATMI,XI,XLIMIT,YDATMA,YDATMI
      INTEGER IDAT,IERR3,IREAD,LUERR,LUIN,LUOUT,NEND
C     ..
C     .. Arrays in Common ..
      DOUBLE PRECISION XDAT(NMAX),YDAT(NMAX)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION JNUXV,NUOUT,X,XIX,Y
      INTEGER IERR4,IERRJ
C     ..
C     .. External Subroutines ..
      EXTERNAL DHJNUX,FUNK1,FUNK2,FUNK3
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC DSQRT
C     ..
C     .. Common blocks ..
      COMMON /ARGUME/XI,NU,XLIMIT
      COMMON /DATEN/XDAT,YDAT,XDATMI,XDATMA,YDATMI,YDATMA,NEND
      COMMON /FHNKL1/IERR3
      COMMON /INOUTE/LUIN,LUOUT,LUERR
      COMMON /LOGIK/IDAT,IREAD
      integer i1mach, nout
      nout = i1mach(2)
C     ..
      X = XIN
      NUOUT = NU
C
      XIX = XI*X
C
C     CALCULATE THE BESSEL FUNCTION *** JNUX ***
C     OF ORDER *** NU *** FOR ARGUMENT *** XIX ***.
      IERRJ = 0
C
      CALL DHJNUX(XIX,NUOUT,IERRJ,JNUXV)
C
      IF (IERRJ.NE.0) THEN
          WRITE (LUERR,FMT=9000)
          WRITE (LUERR,FMT=9010) XIX,XI,X,NU
          IERR3 = 1
          WRITE(NOUT,FMT=9020)
          WRITE (LUERR,FMT=9020)
          STOP

      END IF
C
C     CALCULATE *** FUN *** FOR ARGUMENT *** X ***.
C
      IF (IDAT.EQ.1) THEN
          CALL FUNK1(X,Y,IERR4)

      ELSE IF (IDAT.EQ.2) THEN
          CALL FUNK2(X,Y,IERR4)

      ELSE IF (IDAT.EQ.3) THEN
          CALL FUNK3(X,Y,XDAT,YDAT,NEND,IERR4)
      END IF
C
      IF (IERR4.NE.0) THEN
          WRITE (LUERR,FMT=9030)
          WRITE (LUERR,FMT=9040) X
          WRITE (LUERR,FMT=9050) IERR4
          IERR3 = 1
          FUNINT = 0.0D0
          WRITE(NOUT,FMT=9020)
          WRITE (LUERR,FMT=9020)
          STOP

      END IF
C
C     FINALLY CALCULATE THE INTEGRAND:
      FUNINT = DSQRT(XIX)*JNUXV*Y
C
      RETURN

 9000 FORMAT (/,' MESSAGE FUNINT:',/,
     +       ' ERROR DURING BESSEL FUNCTION CALC ULATION!')
 9010 FORMAT (' MESSAGE FUNINT:',/,' ARGUMENTS TO BESSEL FUNCTION:',/,
     +       ' XI * X =',D17.9,/,' XI = ',D17.9,'     X =',D17.9,/,
     +       ' NU =',D17.9)
 9020 FORMAT (/,' MESSAGE FUNINT:',/,
     +       ' AN ERROR DURING FUNCTION EVALUATION HAS OCCURED!',/,
     +       ' THE HANKEL TRANSFORM HAS FAILED! SORRY! STOP!')
 9030 FORMAT (/,' MESSAGE FUNINT:',/,
     +       ' ERROR DURING CALCULATION OF FUNCTION *** FUN *** !')
 9040 FORMAT (/,' MESSAGE FUNINT:',/,
     +       ' ARGUMENT TO FUNCTION *** FUN **  * :',/,' X =',D17.9)
 9050 FORMAT (1X,' *** IERR4 *** =',I6)
      END
C     END OF FUNCTION *** FUNINT ***
C
C     BEGIN OF SUBROUTINE *** DHRW ***
      SUBROUTINE DHRW(IERR)
C
C     ##################################################################
C
C     BEGONNEN AM:         09.12.96
C     VORLETZTE AENDERUNG: 01.02.99
C     LETZE AENDERUNG AM:  08.02.99
C
C     PROGRAMMED BY:
C     THOMAS WIEDER
C     E-MAIL: WIEDER@HRZ.UNI-KASSEL.DE
C
C     PURPOSE:
C     LESEN VON X-Y-DATEN AUS EINER EINGABEDATEI.
C     SORTIEREN VON X-Y-DATEN.
C     SCHREIBEN VON X-Y-DATEN IN EINE AUSGABEDATEI.
C
C     EINGABEDATEN IM FORMAT:
C     X_1 Y_1
C     X_2 Y_2
C       ...
C     X_N Y_N
C
C     LUIN = LOGICAL NUMBER OF INPUT FILE
C
C     IERR = 0 --> NO ERROR OCCURED.
C     IERR = 1 --> ERROR WHILE TRYING TO OPEN THE INPUT FILE.
C     IERR = 2 --> ERROR WHILE READING FROM INPUT FILE.
C
C     ##################################################################
C
C      IMPLICIT NONE
C
C     .. Parameters ..
      INTEGER NMAX
      PARAMETER (NMAX=1000)
C     ..
C     .. Scalar Arguments ..
      INTEGER IERR
C     ..
C     .. Scalars in Common ..
      DOUBLE PRECISION XDATMA,XDATMI,YDATMA,YDATMI
      INTEGER LUERR,LUIN,LUOUT,NEND
      CHARACTER*80 FILEER,FILEIN,FILEOU
C     ..
C     .. Arrays in Common ..
      DOUBLE PRECISION XDAT(NMAX),YDAT(NMAX)
C     ..
C     .. Local Scalars ..
      INTEGER N
C     ..
C     .. External Functions ..
      DOUBLE PRECISION D1MACH
      EXTERNAL D1MACH
C     ..
C     .. External Subroutines ..
      EXTERNAL DHSORT
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC MAX,MIN
C     ..
C     .. Common blocks ..
      COMMON /DATEN/XDAT,YDAT,XDATMI,XDATMA,YDATMI,YDATMA,NEND
      COMMON /INOUTE/LUIN,LUOUT,LUERR
      COMMON /NAMEN/FILEIN,FILEOU,FILEER
      integer i1mach, nout
      nout = i1mach(2)
C     ..
      IERR = 0
      XDATMI = D1MACH(2)
      XDATMA = -D1MACH(2)
      YDATMI = D1MACH(2)
      YDATMA = -D1MACH(2)
C
      WRITE(NOUT,FMT=9000)
C
C     EINLESEN:
C
C20    WRITE (*,30)
C30    FORMAT (/,' MESSAGE RW: GIVE NAME OF INPUT FILE:',/)
C      READ (*,40,ERR=20) NAME
C40    FORMAT (A40)
C      FILEIN=NAME
C     DATEI OEFFNEN:
      OPEN (LUIN,FILE=FILEIN,STATUS='OLD',ACCESS='SEQUENTIAL',ERR=40)
C
      NEND = 0
      DO 10 N = 1,NMAX
          READ (LUIN,FMT=*,ERR=50,END=20) XDAT(N),YDAT(N)
          NEND = NEND + 1
          IF (NEND.EQ.NMAX) THEN
              WRITE(NOUT,FMT=9010) NEND
              GO TO 20

          END IF

          XDATMI = MIN(XDATMI,XDAT(N))
          XDATMA = MAX(XDATMA,XDAT(N))
          YDATMI = MIN(YDATMI,YDAT(N))
          YDATMA = MAX(YDATMA,YDAT(N))
   10 CONTINUE
C
   20 WRITE(NOUT,FMT=9020) NEND
C
      IF (XDATMI.LT.-D1MACH(2)) THEN
          WRITE(NOUT,FMT=9030)
          WRITE (LUERR,FMT=9030)
          IERR = 1
      END IF

      IF (XDATMI.GT.D1MACH(2)) THEN
          WRITE(NOUT,FMT=9040)
          WRITE (LUERR,FMT=9040)
          IERR = 1
      END IF

      IF (YDATMI.LT.-D1MACH(2)) THEN
          WRITE(NOUT,FMT=9050)
          WRITE (LUERR,FMT=9050)
          IERR = 1
      END IF

      IF (YDATMI.GT.D1MACH(2)) THEN
          WRITE(NOUT,FMT=9060)
          WRITE (LUERR,FMT=9060)
          IERR = 1
      END IF
C
C     DATEI SCHLIESSEN:
      CLOSE (LUIN)
C
C     NACH AUFSTEIGENDEN X-WERTEN SORTIEREN:
C
      CALL DHSORT(XDAT,YDAT,NEND)
C
      DO 30 N = 1,NEND
          WRITE (LUERR,FMT=9070) N,XDAT(N),YDAT(N)
   30 CONTINUE
C
      WRITE(NOUT,FMT=9080)
C
      GO TO 60

   40 WRITE(NOUT,FMT=*) 'ERROR ON OPENING FILE!'
      WRITE (LUERR,FMT=*) 'ERROR ON OPENING FILE! FILE DOES NOT EXIST?'
      IERR = 1
      RETURN

   50 WRITE(NOUT,FMT=*) 'ERROR ON READING DATA!'
      IERR = 2
      RETURN
C
   60 RETURN

 9000 FORMAT (/,' MESSAGE RW: BEGIN OF DATA INPUT FROM FILE!')
 9010 FORMAT (/,' MELDUNG RW:',/,' WARNING! MAXIMAL NUMBER =',I6,
     +       ' OF DATA PAIRS HAVE BEEN READ!')
 9020 FORMAT (/,' MESSAGE RW:',I6,' DATA PAIRS HAVE BEEN READ!')
 9030 FORMAT (/,' MESSAGE RW: TOO SMALL VALUE FOR ABCISSA *** X ***!',/,
     +       ' PLEASE CHECK YOUR DATA!')
 9040 FORMAT (/,' MESSAGE RW: TOO BIG VALUE FOR ABCISSA *** X ***!',/,
     +       ' PLEASE CHECK YOUR DATA!')
 9050 FORMAT (/,' MESSAGE RW: TOO SMALL VALUE FOR ORDINATE *** Y ***!',
     +       /,' PLEASE CHECK YOUR DATA!')
 9060 FORMAT (/,' MESSAGE RW: TOO BIG VALUE FOR ORDINATE *** X ***!',/,
     +       ' PLEASE CHECK YOUR DATA!')
 9070 FORMAT (' MESSAGE RW: N, X(N), Y(N):',I6,3X,D17.9,3X,D17.9)
 9080 FORMAT (/,' MESSAGE RW: END OF READING DATA')
      END
C     END OF SUBROUTINE *** DHRW ***
C
C     BEGIN OF SUBROUTINE *** DHSORT ***
      SUBROUTINE DHSORT(X,Y,N)
C
C     ##################################################################
C
C     BEGONNEN AM:      09.12.96
C     LETZTE AENDERUNG: 04.01.98
C
C     PROGRAMMED BY:
C     THOMAS WIEDER
C     E-MAIL: WIEDER@HRZ.UNI-KASSEL.DE
C
C     PURPOSE:
C     SORTIEREN VON X-Y-DATEN NACH AUFSTEIGENDEN X-WERTEN.
C
C     ##################################################################
C
C      IMPLICIT NONE
C
C     .. Scalar Arguments ..
      INTEGER N
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION X(N),Y(N)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION XHILF,YHILF
      INTEGER I,J
C     ..
      DO 20 I = 1,N
          DO 10 J = I + 1,N
              IF (X(J).LT.X(I)) THEN
C     DIE BEIDEN DATENPAARE AUSTAUSCHEN:
                  XHILF = X(I)
                  YHILF = Y(I)
                  X(I) = X(J)
                  Y(I) = Y(J)
                  X(J) = XHILF
                  Y(J) = YHILF
              END IF

   10     CONTINUE
   20 CONTINUE
C
      RETURN

      END
C     END OF SUBROUTINE *** DHSORT ***
C
C     BEGIN OF SUBROUTINE *** DIVDIF ***
      SUBROUTINE DIVDIF(X,Y,N,M,IDEG,TABLE,XARG,VALUE,TRUBLE)
C      SUBROUTINE DIVDIF (X,Y,N,M,IDEG,XARG,VALUE,TRUBLE)
C
C     FIRST VERSION:    24.03.1988
C     PREVIOUS VERSION: 13.02.1996
C     LATEST VERSION:   15.04.1998
C
C     SOURCE: ???
C     I AM TRYING TO FIND THE SOURCE OF SUBROUTINE **** DIVDIF ***.
C     IN THE MEANTIME, I ASK THE AUTHOR(S) FOR EXCUSE.
C     I THINK, **** SUBROUTINE DIVDIF *** WAS GIVEN IN A BOOK
C     ON FORTRAN PROGRAMMING WHICH I READ DURING MY TIME AS PHD STUDENT.
C
C     PURPOSE:
C     INTERPOLATION
C
C     METHOD:
C     NEWTON'S DIVIDED DIFFERENCES
C
C      IMPLICIT NONE
C     .. Scalar Arguments ..
      DOUBLE PRECISION VALUE,XARG
      INTEGER IDEG,M,N,TRUBLE
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION TABLE(N,N),X(N),Y(N)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION YEST
      INTEGER I,IDEGM1,ISUB,ISUB1,ISUB2,J,MAX,NM1
C     ..
      TRUBLE = 0
C     CHECK FOR ARGUMENT CONSISTENCY
      IF (M.LT.N) GO TO 10
      TRUBLE = 1
      RETURN
C     CALCULATE FIRST-ORDER DIFFERENCES
   10 NM1 = N - 1
      DO 20 I = 1,NM1
          TABLE(I,1) = (Y(I+1)-Y(I))/ (X(I+1)-X(I))
   20 CONTINUE
      IF (M.LE.1) GO TO 50
C     CALCULATE HIGHER ORDER DIFFERENCES
      DO 40 J = 2,M
          DO 30 I = J,NM1
              ISUB = I + 1 - J
              TABLE(I,J) = (TABLE(I,J-1)-TABLE(I-1,J-1))/
     +                     (X(I+1)-X(ISUB))
   30     CONTINUE
   40 CONTINUE
   50 TRUBLE = 0
C     CHECK FOR ARGUMENT CONSISTENCY
      IF (IDEG.LE.M) GO TO 60
      TRUBLE = 2
      VALUE = 0.0D0
      RETURN
C     SEARCH X VECTOR FOR ELEMENT .GE. XARG
   60 DO 70 I = 1,N
          IF (I.EQ.N .OR. XARG.LE.X(I)) GO TO 80
   70 CONTINUE
   80 MAX = I + IDEG/2
C     INSURE THAT ALL REQUIRED DIFFERENCES ARE IN TABLE
      IF (MAX.LE.IDEG) MAX = IDEG + 1
      IF (MAX.GT.N) MAX = N
C     COMPUTE INTERPOLANT VALUE
      YEST = TABLE(MAX-1,IDEG)
      IF (IDEG.LE.1) GO TO 100
      IDEGM1 = IDEG - 1
      DO 90 I = 1,IDEGM1
          ISUB1 = MAX - I
          ISUB2 = IDEG - I
          YEST = YEST* (XARG-X(ISUB1)) + TABLE(ISUB1-1,ISUB2)
   90 CONTINUE
  100 ISUB1 = MAX - IDEG
      TRUBLE = 0
      VALUE = YEST* (XARG-X(ISUB1)) + Y(ISUB1)
      RETURN

      END
C     END OF SUBROUTINE *** DIVDIF ***
C
C     BEGIN OF SUBROUTINE *** FUNCHE ***
      SUBROUTINE FUNCHE(IERR7)
C
C     ##################################################################
C
C     BEGONNEN AM:         28.01.98
C     VORLETZTE AENDERUNG: 13.10.98
C     LATEST VERSION:      03.02.99
C
C     PROGRAMMED BY:
C     THOMAS WIEDER
C     E-MAIL: WIEDER@HRZ.UNI-KASSEL.DE
C
C     PURPOSE:
C     SUBROUTINE *** FUNCHE *** CALCULATES FUNCTION *** FUN ***
C     FOR VARIUS VALUES OF ARGUMENT *** X ***
C     IN ORDER TO LEARN ABOUT *** FUN ***.
C
C     CALLING SEQUENCE:
C     CALLED BY SUBROUTINE *** HANKEL ***.
C     CALLS SUBROUTINE *** FUNK ***.
C
C     MOST IMPORTANT VARIABLES:
C
C     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C     X IS THE INPUT VARIABLE TO SUBROUTINE *** FUNK ***.
C     X IS A REAL-VALUED VARIBALE.
C     X IS THE INDEPENDENT VARIBALE OF FUNCTION *** FUN ***.
C     X >= 0.0D0
C     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C
C     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C     Y IS THE OUTPUT VARIABLE FROM SUBROUTINE *** FUNK ***.
C     Y IS THE DEPENDENT VARIBALE OF FUNCTION *** FUN ***.
C     Y IS A REAL-VALUED VARIBALE.
C     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C
C     XMAX = LARGEST POSITVE MAGNITUDE OF DOUBLE PRECISION NUMBERS
C            PROVDIED BY THE COMPUTER
C     XTINY = SMALLEST POSITIVE MAGNITUDE OF DOUBLE PRECISION NUMBERS
C             PROVDIED BY THE COMPUTER
C
C     ##################################################################
C
C      IMPLICIT NONE
C     .. Parameters ..
      DOUBLE PRECISION ZERO
      PARAMETER (ZERO=0.0D0)
C     ..
C     .. Scalar Arguments ..
      INTEGER IERR7
C     ..
C     .. Scalars in Common ..
      INTEGER LUERR,LUIN,LUOUT
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION NU,X,XLIMIT,XMAX,Y
      INTEGER IERR5,N
C     ..
C     .. Local Arrays ..
      DOUBLE PRECISION NUDAT(4)
C     ..
C     .. External Functions ..
      DOUBLE PRECISION D1MACH
      EXTERNAL D1MACH
C     ..
C     .. External Subroutines ..
      EXTERNAL DHFUNC
C     ..
C     .. Common blocks ..
      COMMON /INOUTE/LUIN,LUOUT,LUERR
C     ..
C     .. Data statements ..
      DATA NUDAT/-0.5D0,0.0D0,3.5D0,10.0D0/
C     ..
C
C     SET SOME PARAMETERS:
C
      IERR7 = 0
C      XTINY=D1MACH(1)
      XMAX = D1MACH(2)
C
C     DO THE CALCLULATION:
C
C     WE SELECT SOME VALUES FOR *** NU ***
      DO 10 N = 1,4
          NU = NUDAT(N)
C
C     CHECK *** FUN *** FOR *** X *** = 0.
C
          X = ZERO
          Y = ZERO
C
          IERR5 = 0
          CALL DHFUNC(NU,X,Y,XLIMIT,IERR5)
C
          IF (IERR5.NE.0) THEN
              WRITE (LUERR,FMT=9000) IERR5
              IERR7 = 1
          END IF

          IF (Y.LT.-XMAX) THEN
              WRITE (LUERR,FMT=9010) X,Y
              WRITE (LUERR,FMT=9020)
          END IF

          IF (Y.GT.XMAX) THEN
              WRITE (LUERR,FMT=9030) X,Y
              WRITE (LUERR,FMT=9020)
          END IF
C
C     CHECK *** FUN *** FOR *** X *** = *** XMAX ***.
C
          X = XMAX
          Y = ZERO
C
          CALL DHFUNC(NU,X,Y,XLIMIT,IERR5)
C
          IF (IERR5.NE.0) THEN
              WRITE (LUERR,FMT=9000) IERR5
              IERR7 = 1
          END IF

          IF (Y.LT.-XMAX) THEN
              WRITE (LUERR,FMT=9040) X,Y
              WRITE (LUERR,FMT=9020)
          END IF

          IF (Y.GT.XMAX) THEN
              WRITE (LUERR,FMT=9030) X,Y
              WRITE (LUERR,FMT=9020)
          END IF
C
          IF (IERR7.NE.0) THEN
C     AN ERROR DURING CALL OF FUNCTION *** DHFUNC *** HAS OCCURRED.
              WRITE (LUERR,FMT=9050) NU,X,Y
          END IF
C
   10 CONTINUE
C
      RETURN

 9000 FORMAT (/,' MESSAGE FUNCHE:',/,
     +       ' AN ERROR OCCURED DURING CALL TO  *** DHFUNC ***!',/,
     +       ' ERROR FLAG FROM *** DHFUNC *** : ***  IERR *** =,',I6)
 9010 FORMAT (/,' MESSAGE FUNCHE:',/,' Y = FUN(X=',D19.7,') = ',D19.7,/,
     +       ' Y TOO SMALL, OVERFLOW!')
 9020 FORMAT (/,' MESSAGE FUNCHE:',/,
     +       ' THE HANKEL TRANSFORM POSSIBLY WI LL FAIL!',/,
     +       ' PLEASE CHECK YOUR FUNCTION *** FUN *** !')
 9030 FORMAT (/,' MESSAGE FUNCHE:',/,' Y = FUN(X=',D19.7,') = ',D19.7,/,
     +       ' Y TOO BIG, OVERFLOW!')
 9040 FORMAT (/,' MESSAGE FUNCHE:',/,' Y = FUN(X=',D19.7,') = ',D19.7,/,
     +       ' Y TOO SMALL, OVERFLOW!')
 9050 FORMAT (/,' MESSAGE FUNCHE:',/,
     +       ' AN ERROR DURING CALL OF SUBROUTINE *** FUNC ***',/,
     +       ' HAS BEEN REPORTED:',/,' NU, X, Y =',I6,5X,D17.9,5X,D17.9)
      END
C     END OF SUBROUTINE *** FUNCHE ***
C
C     BEGIN OF SUBROUTINE *** ERRMSS ***
      SUBROUTINE ERRMSS(IERRME)
C
C     ##################################################################
C
C     FIRST VERSION:    14.01.1998
C     PREVIOUS VERSION: 14.02.1998
C     LATEST VERSION:   13.10.1998
C
C     PROGRAMMED BY:
C     THOMAS WIEDER
C     E-MAIL: WIEDER@HRZ.UNI-KASSEL.DE
C
C     PURPOSE:
C     OUTPUTS ERROR MESSAGES.
C
C     MOST IMPORTANT VARIABLES:
C
C     IERRME = FLAG FOR ERROR MESSAGES
C
C     ##################################################################
C
C      IMPLICIT NONE
C
C     .. Scalar Arguments ..
      INTEGER IERRME
C     ..
C     .. Scalars in Common ..
      INTEGER LUERR,LUIN,LUOUT
      CHARACTER*80 FILEER,FILEIN,FILEOU
C     ..
C     .. Common blocks ..
      COMMON /INOUTE/LUIN,LUOUT,LUERR
      COMMON /NAMEN/FILEIN,FILEOU,FILEER
      integer i1mach, nout
      nout = i1mach(2)
C     ..
      WRITE(NOUT,FMT=9000) FILEER
      WRITE (LUERR,FMT=9000) FILEER
      WRITE (LUERR,FMT=9010) IERRME
C
      IF (IERRME.EQ.-10) THEN
          WRITE (LUERR,FMT=9020)

      ELSE IF (IERRME.EQ.10) THEN
          WRITE (LUERR,FMT=9030)
          WRITE (LUERR,FMT=9040)

      ELSE IF (IERRME.EQ.11) THEN
          WRITE (LUERR,FMT=9050)
          WRITE (LUERR,FMT=9060)

      ELSE
          WRITE(NOUT,FMT=9070) IERRME
          WRITE (LUERR,FMT=9070) IERRME
      END IF
C
      RETURN

 9000 FORMAT (/,' MESSAGE ERRMSS:',/,
     +       ' AN ERROR OCCURED DURING HANKEL TRANSFORMATION!!!',/,
     +       ' RESULTS WILL BE UNRELIABLE OR WRONG!!!',/,
     +       ' PLEASE LOOK INTO FILE:',/,A80)
 9010 FORMAT (/,' MESSAGE ERRMS: ERROR FLAG *** IERRME *** =',I6)
 9020 FORMAT (/,' MESSAGE ERRMS: PLEASE CHECK FILE NAMES, PATH ETC.!')
 9030 FORMAT (/,
     +       ' MESSAGE ERRMS: THE INTEGRAND COULD NOT BE CALCULATED!!!')
 9040 FORMAT (' PLEASE CHECK YOUR FUNCTION *** FUN ***!',/,
     +       'ARE THERE    ANY SINGULARITIES?',/,
     +       'DOES THE FUNCTION GO TO INFINITY (FOR LARGE ARGUMENTS)?')
 9050 FORMAT (/,' MESSAGE ERRMS: THE INTEGRATION COULD NOT BE DONE!!!')
 9060 FORMAT (' THE FUNCTION *** FUN *** MUST GO TO ZERO ',/,
     +       ' FOR LA   RGE ARGUMENTS!')
 9070 FORMAT (/,' MESSAGE ERRMSS:',/,
     +       ' UNKNOWN ERROR IDENTIFIER *** IERRMESS *** =',I6)
      END
C     END OF SUBROUTINE *** ERRMSS ***
C
C     BEGIN OF SUBROUTINE *** DHJNUX ***
      SUBROUTINE DHJNUX(XIN,NUIN,IERR,JNUXV)
C     ##################################################################
C
C     FIRST VERSION:    22.01.98
C     PREVIOUS VERSION: 04.02.99
C     LATEST VERSION:   04.02.99
C
C     PROGRAMMED BY:
C     THOMAS WIEDER
C     E-MAIL: WIEDER@HRZ.UNI-KASSEL.DE
C
C     PURPOSE:
C
C     SUBROUTINE *** DHJNUX *** CALCULATES THE BESSEL FUNCTION
C     OF ORDER *** NU *** FOR ARGUMENT *** X ***.
C
C     CALLING SEQUENCE:
C
C     CALLED BY SUBROUTINE *** FUNINT ***.
C     CALLS SUBROUTINE *** RJBESL ***.
C     CALLS SUBROUTINE *** BSSJYA ***.
C
C     MOST IMPORTANT VARIABLES:
C
C     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C     X IS AN INPUT VARIABLE TO FUNCTION *** DHJNUX ***.
C     X IS A REAL-VALUED VARIBALE!
C     X >= 0.0D0
C     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C
C     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C     NU IS AN INPUT VARIABLE TO FUNCTION *** JNUX ***.
C     NU IS A REAL-VALUED VARIBALE!
C     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C
C     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C     JNUXV IS AN OUTPUT VARIABLE OF FUNCTION *** JNUX ***.
C     JNUXV IS A REAL-VALUED VARIBALE!
C     JNUXV IS THE RESULT: JNUXV=JNUX(X)
C     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C
C     ##################################################################
C
C      IMPLICIT NONE
C
C     SET SOME PARAMETERS:
C
C     .. Parameters ..
      DOUBLE PRECISION PI,WINZIG,ZERO
      PARAMETER (PI=3.141592653D0,WINZIG=1.0D-09,ZERO=0.0D0)
      INTEGER NUMIMA
      PARAMETER (NUMIMA=101)
C     ..
C     .. Scalar Arguments ..
      DOUBLE PRECISION JNUXV,NUIN,XIN
      INTEGER IERR
C     ..
C     .. Scalars in Common ..
      INTEGER LUERR,LUIN,LUOUT
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION ALPHA,ARG,NU,NUMAX,NUMIN,X,XASYMP,XMAX,XMIN,YNUXV
      INTEGER NB,NCALC
C     ..
C     .. Local Arrays ..
      DOUBLE PRECISION RESULT(NUMIMA)
C     ..
C     .. External Functions ..
      DOUBLE PRECISION D1MACH
      EXTERNAL D1MACH
C     ..
C     .. External Subroutines ..
      EXTERNAL BSSJYA,RJBESL,RYBESL
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC DBLE,DCOS,DSIN,IDINT
C     ..
C     .. Common blocks ..
      COMMON /INOUTE/LUIN,LUOUT,LUERR
      integer i1mach, nout
      nout = i1mach(2)
C     ..
      IERR = 0
C
C     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C     *** 24.02.1998 *** A GOOD VALUE FOR *** XASYMP *** = 100.0
      XASYMP = 100.0D0
C     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C
C     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      NUMIN = -DBLE(NUMIMA)
      NUMAX = DBLE(NUMIMA)
C     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C
C     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      XMIN = ZERO
      XMAX = D1MACH(2)
C     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C
C     CHECK INPUT:
C
      X = XIN
      NU = NUIN
C
C     IS *** X *** POSITIVE?
      IF (X.LT.ZERO) THEN
          WRITE(NOUT,FMT=9000)
          WRITE (LUERR,FMT=9000)
          WRITE(NOUT,FMT=9010) X
          WRITE (LUERR,FMT=9010) X
          IERR = 1
          GO TO 10

      END IF
C
C     CHECK ABSOLUT VALUE OF *** X ***
      IF (X.LT.XMIN .OR. X.GT.XMAX) THEN
          WRITE(NOUT,FMT=9000)
          WRITE (LUERR,FMT=9000)
          WRITE(NOUT,FMT=9020) X
          WRITE (LUERR,FMT=9020) X
          WRITE (LUERR,FMT=9030) XMIN,XMAX
          IERR = 1
          GO TO 10

      END IF
C
C     CHECK ABSOLUT VALUE OF *** NU ***
      IF (NU.LT.NUMIN .OR. NU.GT.NUMAX) THEN
          WRITE(NOUT,FMT=9000)
          WRITE (LUERR,FMT=9000)
          WRITE(NOUT,FMT=9040) NU
          WRITE (LUERR,FMT=9040) NU
          WRITE (LUERR,FMT=9050) NUMIN,NUMAX
          IERR = 1
          GO TO 10

      END IF
C
C     DO THE CALCULATION:
C
C     ###############################################################
C     X = 0.0
C     ###############################################################
C
C     THE BESSEL FUNCTION IS DEFINED FOR X=0.0 ALSO:
      IF (X.EQ.ZERO) THEN
          IF (NU.EQ.ZERO) JNUXV = 1.0D0
          IF (NU.GT.ZERO) JNUXV = 0.0D0
          IF (NU.LT.ZERO) JNUXV = 0.0D0
          RETURN

      END IF
C
C     ###############################################################
C     X > 0.0
C     ###############################################################
C
      IF (NU.GE.0) THEN
C
          IF (X.LE.XASYMP) THEN
C
C     PREVENT NUMERICAL INSTABILITIES FOR VERY SMALL ARGUMENTS:
              IF (X.LT.WINZIG) THEN
                  IF (NU.EQ.ZERO) JNUXV = 1.0D0
                  IF (NU.GT.ZERO) JNUXV = 0.0D0
                  RETURN

              END IF
C
              ALPHA = NU - IDINT(NU)
              NB = IDINT(NU) + 1
C
              CALL RJBESL(X,ALPHA,NB,RESULT,NCALC)
C
              IF (NCALC.NE.NB) THEN
                  WRITE (LUERR,FMT=9060)
                  WRITE (LUERR,FMT=9070) NCALC
                  WRITE (LUERR,FMT=9080) X,NB,RESULT(NB)
C      IERR=1
C      WRITE (LUERR,*) 'IERR =,IERR
              END IF

              JNUXV = RESULT(NB)
C
          ELSE IF (X.GT.XASYMP) THEN
C     USE AN ASYMPTOTIC FORMULA FOR THE BESSEL FUNCTION:
C
              CALL BSSJYA(X,NU,JNUXV)
C
C     END OF IF-CLAUSE *** X ***
          END IF
C
C     #################################################################
C     NU < 0
C     #################################################################
C
      ELSE IF (NU.LT.0.0D0) THEN
C
          IF (X.LE.XASYMP) THEN
C
C     PREVENT NUMERICAL INSTABILITIES FOR VERY SMALL ARGUMENTS:
              IF (X.LT.WINZIG) THEN
                  JNUXV = 0.0D0
                  RETURN

              END IF
C
              NU = -NU
              ALPHA = NU - IDINT(NU)
              NB = IDINT(NU) + 1
C
              CALL RJBESL(X,ALPHA,NB,RESULT,NCALC)
C
              IF (NCALC.NE.NB) THEN
                  WRITE (LUERR,FMT=9060)
                  WRITE (LUERR,FMT=9070) NCALC
                  WRITE (LUERR,FMT=9080) X,NB,RESULT(NB)
                  IERR = 2
                  WRITE (LUERR,FMT=*) 'IERR =',IERR
              END IF

              NU = -NU
              JNUXV = RESULT(NB)
C
          ELSE IF (X.GT.XASYMP) THEN
C     USE AN ASYMPTOTIC FORMULA FOR THE BESSEL FUNCTION:
C
              CALL BSSJYA(X,-NU,JNUXV)
C
C     END OF IF-CLAUSE *** X ***
          END IF
C
C     CALCULATE THE BESSEL FUNCTION OF THE SECOND KIND Y(NU,X):
C
          NU = -NU
          ALPHA = NU - IDINT(NU)
          NB = IDINT(NU) + 1
C
C     TO CIRCUMVENT PROBLEMS WITH SMALL ARGUMENTS *** X ***
C     IN CASE OF NEGATIVE *** NU ***:
C
          IF (X.LT.1.0D-03) THEN
C      YNUXV=-1.0D0+X
              YNUXV = ZERO
C
          ELSE
C
              CALL RYBESL(X,ALPHA,NB,RESULT,NCALC)
C
              IF (NCALC.LT.0) THEN
                  WRITE (LUERR,FMT=9090)
                  WRITE (LUERR,FMT=9100) NCALC
                  WRITE (LUERR,FMT=9080) X,NB,RESULT(NB)
                  IERR = 3
                  WRITE (LUERR,FMT=*) 'IERR =',IERR
              END IF

              NU = -NU
              YNUXV = RESULT(NB)
C
C     END OF IF-CLAUSE *** X ***
          END IF
C
C     FOR NEGATIVE *** NU *** SEE:
C     G. KORN AND T. KORN, MC GRAW HILL, 1968, CHAPTER 21.8-4.
C     W.H. PRESS, S.A. TEUKOLSKY, W.T. VETTERLING, B.P. FLANNERY:
C     NUMERICAL RECIPES, SECOND EDITION 1994, CAMBRDIGE UNIVERSITY PRESS
C     EQUATION (6.7.19):
C
C      ARG=DMOD(-NU*PI,2.0D0*PI)
          ARG = -NU*PI
          JNUXV = DCOS(ARG)*JNUXV - DSIN(ARG)*YNUXV
C
C     END OF IF-CLAUSE *** NU.LT.0.0D0 ***.
      END IF
C
   10 RETURN

 9000 FORMAT (/,' MESSAGE JNUX: WRONG INPUT VALUE!')
 9010 FORMAT (/,' MESSAGE JNUX:',/,' PARAMETER *** X *** IS NEGATIVE!',
     +       /,' *** X *** =',D17.9)
 9020 FORMAT (/,' MESSAGE JNUX:',/,
     +       ' PARAMETER *** X *** OUTSIDE ALLOWED RANGE!',/,
     +       ' *** X *** =',D17.9)
 9030 FORMAT (/,' ALLOWED INTERVAL [X_MIN,X_MAX]:',/,'[',D17.9,',',
     +       D17.9,']')
 9040 FORMAT (/,' MESSAGE HANKEL:',/,
     +       ' PARAMETER *** NU *** OUTSIDE ALLOWED RANGE!',/,
     +       ' *** NU *** =',D17.9)
 9050 FORMAT (/,' ALLOWED INTERVAL [NU_MIN,NU_MAX]:',/,'[',D17.9,',',
     +       D17.9,']')
 9060 FORMAT (/,'MESSAGE JNUX: NCALC NOT EQUAL NB! WARNING!',/,
     +       ' THE ACCURACY OF THE RESULTS MAY BE LOW!')
 9070 FORMAT (/,' MESSAGE JNUX: ',/,
     +       ' NCALC ON RETURN FORM SUBROUTINE    *** RJBESL *** =',I6)
 9080 FORMAT (1X,'X, NB, RESULT(NB):',D17.9,3X,I6,3X,D17.9)
 9090 FORMAT (/,'MESSAGE JNUX: NCALC NOT EQUAL NB! ERROR!')
 9100 FORMAT (/,' MESSAGE JNUX: ',/,
     +       ' NCALC ON RETURN FORM SUBROUTINE   *** RYBESL *** =',I6)
      END
C     END OF SUBROUTINE *** DHJNUX ***
C
C     BEGIN OF SUBROUTINE *** BSSJYA ***
      SUBROUTINE BSSJYA(X,NU,RESULT)
C
C     ##################################################################
C
C     FIRST VERSION: 12.11.1997
C     PREVIOUS VERSION: 29.01.1998
C     LATEST VERSION: 03.02.1998
C
C     PROGRAMMED BY:
C     THOMAS WIEDER
C     E-MAIL: WIEDER@HRZ.UNI-KASSEL.DE
C
C     THIS PROGRAM CALCULATES THE ASYMPTOTIC EXPANSION OF
C     THE FIRST BESSEL FUNCTION J(NU,X) FOR LARGE ARGUMENT X.
C
C     LITERATURE:
C
C     1.) E.T. WHITTAKER AND G.N. WATSON:
C         ACOURSE OF MODERN ANALYSIS,
C         CAMBRIDGE UNIVERSITY PRESS, LONDON,
C         4TH EDITION, REPRINT 1978,
C         P. 368.
C     2.) R. COURANT AND D. HILBERT:
C         METHODEN DER MATHEMATISCHEN PHYSIK,
C         SPRINGER VERLAG, BERLIN, 3TE AUFLAGE, 1968,
C         KAPITEL VII.
C
C     ##################################################################
C
C      IMPLICIT NONE
C
C     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C     *** 24.02.1998 *** A GOOD VALUE FOR *** XASYMP *** = 100.0
C     .. Parameters ..
      DOUBLE PRECISION PI,ZERO
      PARAMETER (PI=3.141592653D0,ZERO=0.0D0)
C     ..
C     .. Scalar Arguments ..
      DOUBLE PRECISION NU,RESULT,X
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION ARG,XASYMP
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC DCOS,DSQRT
      integer i1mach, nout
      nout = i1mach(2)
C     ..
      XASYMP = 100.0D0
C     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C
      IF (X.LT.XASYMP) THEN
          WRITE (68,FMT=9000) X
          WRITE(NOUT,FMT=9000) X
          STOP

      END IF
C
      IF (NU.LT.ZERO) THEN
          WRITE (68,FMT=9010) NU
          WRITE(NOUT,FMT=9010) NU
          STOP

      END IF
C
      ARG = X - (PI*NU/2.0D0) - PI/4.0D0
      RESULT = DSQRT(2.0D0/ (PI*X))*DCOS(ARG)
C
      RETURN

 9000 FORMAT (/,'MESSAGE BSSJYA: ARGUMENT *** X *** TOO LOW:',D17.9)
 9010 FORMAT (/,'MESSAGE BSSJYA: ARGUMENT *** NU *** NEGATIVE:',D17.9)
      END
C     END OF SUBROUTINE *** BSSJYA ***
C
C     BEGIN OF FUNCTION *** TFUNC ***
      DOUBLE PRECISION FUNCTION TFUNC(XIN)
C
C     #################################################################
C
C     FIRST VERSION: 05.02.1998
C     PREVIOUS VERSION: 05.02.1998
C     LATEST VERSION: 06.02.1998
C
C     PROGRAMMED BY:
C     THOMAS WIEDER
C     E-MAIL: WIEDER@HRZ.UNI-KASSEL.DE
C
C     PURPOSE:
C     TEST SUBROUTINE *** INTHP ***.
C
C     #################################################################
C
C      IMPLICIT NONE
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION XIN
C     ..
C     .. Scalars in Common ..
      DOUBLE PRECISION EXACT,XSPLIT
      INTEGER IINT
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION JNUXV,NU,X
      INTEGER IERRJ
C     ..
C     .. External Subroutines ..
      EXTERNAL DHJNUX
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC DSIN,EXP
C     ..
C     .. Common blocks ..
      COMMON /LOSUNG/XSPLIT,EXACT,IINT
      integer i1mach, nout
      nout = i1mach(2)
C     ..
      IERRJ = 0
      X = XIN
C
      IF (IINT.EQ.1) THEN
          TFUNC = EXP(-X)
C     SOLUTION: 1.0
          EXACT = 1.0D0
C
      ELSE IF (IINT.EQ.2) THEN
          TFUNC = 1.0D0/ (X**2+1.0D0)
C     SOLUTION: PI/2
          EXACT = 1.570796327D0
C
      ELSE IF (IINT.EQ.3) THEN
          TFUNC = 1.0D0/ (X**1.1D0+1.0D0)
C     SOLUTION: 10.13724986
          EXACT = 10.13724986D0
C
      ELSE IF (IINT.EQ.4) THEN
          TFUNC = DSIN(X)/ (X+1.0D0)
C     SOLUTION:0.6214496243
          EXACT = 0.6214496243D0
C
      ELSE IF (IINT.EQ.5) THEN
          NU = 1.0D0
C
          CALL DHJNUX(X,NU,IERRJ,JNUXV)
C
          TFUNC = JNUXV*EXP(-X)
          IF (IERRJ.NE.0) THEN
              WRITE(NOUT,FMT=9000) IERRJ
          END IF
C     SOLUTION: 1/(2+sqrt(2)) = 0.292893218
          EXACT = 0.292893218D0
C
      ELSE IF (IINT.EQ.6) THEN
          NU = 1.0D0
C
          CALL DHJNUX(X,NU,IERRJ,JNUXV)
C
          TFUNC = JNUXV
          IF (IERRJ.NE.0) THEN
              WRITE(NOUT,FMT=9000) IERRJ
          END IF
C     SOLUTION: 1.0
          EXACT = 1.0D0
C
      ELSE IF (IINT.EQ.7) THEN
          NU = 1.0D0
C
          CALL DHJNUX(X,NU,IERRJ,JNUXV)
C
          TFUNC = (1.0D0/ (X+1.0D0))*JNUXV
          IF (IERRJ.NE.0) THEN
              WRITE(NOUT,FMT=9000) IERRJ
          END IF
C     SOLUTION: 0.3980927698
          EXACT = 0.3980927698D0
C
      ELSE IF (IINT.EQ.8) THEN
          NU = 3.0D0
C
          CALL DHJNUX(X,NU,IERRJ,JNUXV)
C
          TFUNC = (1.0D0/ (X+1.0D0))*JNUXV
          IF (IERRJ.NE.0) THEN
              WRITE(NOUT,FMT=9000) IERRJ
          END IF
C     SOLUTION: 0.2464041090
          EXACT = 0.2464041090D0
C
      ELSE IF (IINT.EQ.9) THEN
          NU = 5.0D0
C
          CALL DHJNUX(X,NU,IERRJ,JNUXV)
C
          TFUNC = (1.0D0/ (X+1.0D0))*JNUXV
          IF (IERRJ.NE.0) THEN
              WRITE(NOUT,FMT=9000) IERRJ
          END IF
C     SOLUTION: 0.1659094729
          EXACT = 0.1659094729D0
C
      ELSE
          WRITE(NOUT,FMT=*) 'UNKNOWN INTEGRAND IDENTIFICATION NUMBER!'
          STOP

      END IF
C
      RETURN

 9000 FORMAT (/,' AN ERROR OCCURED DURING CALL OF *** DHJNUX ***',I6)
      END
C     END OF FUNCTION *** TFUNC ***
C
C     BEGIN OF FUNCTION *** HFUN ***
      DOUBLE PRECISION FUNCTION HFUN(X)
C
C     #################################################################
C
C     FIRST VERSION:    05.02.1998
C     PREVIOUS VERSION: 15.04.1998
C     LATEST VERSION:   08.10.1998
C
C     PROGRAMMED BY:
C     THOMAS WIEDER
C     E-MAIL: WIEDER@HRZ.UNI-KASSEL.DE
C
C     PURPOSE:
C     CALCULATE SUBROUTINE *** DHFUNC ***.
C
C     CALLING SEQUENCE:
C     CALLED BY SUBROUTINE *** OSCINT ***.
C
C     #################################################################
C
C      IMPLICIT NONE
C
C     .. Parameters ..
      INTEGER NMAX
      PARAMETER (NMAX=1000)
C     ..
C     .. Scalar Arguments ..
      DOUBLE PRECISION X
C     ..
C     .. Scalars in Common ..
      INTEGER IDAT,IREAD
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION Y
      INTEGER IDUMMY,IERR4
C     ..
C     .. Local Arrays ..
      DOUBLE PRECISION DUMMYA(NMAX)
C     ..
C     .. External Subroutines ..
      EXTERNAL FUNK1,FUNK2,FUNK3
C     ..
C     .. Common blocks ..
      COMMON /LOGIK/IDAT,IREAD
C     ..
      IDUMMY = 12
      DUMMYA(1) = 0.0D0
C
      IF (IDAT.EQ.1) THEN
          CALL FUNK1(X,Y,IERR4)

      ELSE IF (IDAT.EQ.2) THEN
          CALL FUNK2(X,Y,IERR4)

      ELSE IF (IDAT.EQ.3) THEN
          CALL FUNK3(X,Y,DUMMYA,DUMMYA,IDUMMY,IERR4)
      END IF
C
      IF (IERR4.NE.0) THEN
          WRITE (68,FMT=9000)
      END IF
C
      HFUN = Y
C
      RETURN

 9000 FORMAT (/,' MESSAGE HFUN:',/,
     +       ' ERROR MESSAGE FROM SUBROUTINE *** DHFUNK *** ',/,
     +       ' *** IERR4 *** = ',I6)
      END
C     END OF FUNCTION *** HFUN ***
C
C     BEGIN OF FUNCTION *** GPER ***
      DOUBLE PRECISION FUNCTION GPER(X)
C
C     #################################################################
C
C     FIRST VERSION:    05.02.1998
C     PREVIOUS VERSION: 19.02.1998
C     LATEST VERSION:   24.02.1998
C
C     PROGRAMMED BY:
C     THOMAS WIEDER
C     E-MAIL: WIEDER@HRZ.UNI-KASSEL.DE
C
C     PURPOSE:
C     CALCULATE THE OSCILLATING PART OF THE INTEGRAND.
C
C     CALLING SEQUENCE:
C     CALLED BY SUBROUTINE *** OSCINT ***.
C
C     #################################################################
C
C      IMPLICIT NONE
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION X
C     ..
C     .. Scalars in Common ..
      DOUBLE PRECISION NU,XI,XLIMIT
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION JNUXV,XIX
      INTEGER IERRJ
C     ..
C     .. External Subroutines ..
      EXTERNAL DHJNUX
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC DSQRT
C     ..
C     .. Common blocks ..
      COMMON /ARGUME/XI,NU,XLIMIT
      integer i1mach, nout
      nout = i1mach(2)
C     ..
      IERRJ = 0
C
      XIX = XI*X
C
      CALL DHJNUX(XIX,NU,IERRJ,JNUXV)
C
      GPER = DSQRT(XIX)*JNUXV
C
      IF (IERRJ.NE.0) THEN
          WRITE(NOUT,FMT=9000) IERRJ
          WRITE(NOUT,FMT=9010) X
      END IF
C
      RETURN

 9000 FORMAT (/,' MESSAGE GPER:',/,
     +       ' AN ERROR OCCURED DURING CALL OF *** DHJNUX ***',/,
     +       ' *** IERRJ *** =',I6)
 9010 FORMAT (' ARGUMENT *** X *** =',D17.9)
      END
C     END OF FUNCTION *** GPER ***
C
C     BEGIN OF SUBROUTINE *** DGAMMA ***
      DOUBLE PRECISION FUNCTION DGAMMA(XIN)
C
C     #################################################################
C
C     FIRST VERSION: 05.02.1998
C     PREVIOUS VERSION: 05.02.1998
C     LATEST VERSION: 06.02.1998
C
C     PROGRAMMED BY:
C     THOMAS WIEDER
C     E-MAIL: WIEDER@HRZ.UNI-KASSEL.DE
C
C     PURPOSE:
C
C     SUPPORT TO FUNCTION *** DLGAMA ***.
C
C     #################################################################
C
C      IMPLICIT NONE
C     .. Parameters ..
      DOUBLE PRECISION ZERO
      PARAMETER (ZERO=0.0D0)
C     ..
C     .. Scalar Arguments ..
      DOUBLE PRECISION XIN
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION X
C     ..
C     .. External Functions ..
      DOUBLE PRECISION DLGAMA
      EXTERNAL DLGAMA
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC EXP
      integer i1mach, nout
      nout = i1mach(2)
C     ..
      X = XIN
      IF (X.GT.ZERO) THEN
          DGAMMA = DLGAMA(X)
          DGAMMA = EXP(DGAMMA)

      ELSE IF (X.EQ.ZERO) THEN
          DGAMMA = ZERO

      ELSE
          WRITE(NOUT,FMT=9000) X
          STOP

      END IF

      RETURN

 9000 FORMAT (/,' MESSAGE DGAMMA: WRONG INPUT VALUE X =',D17.9)
      END
