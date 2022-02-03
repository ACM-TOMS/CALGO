C  PROGRAM DRIVE

C     **********
C     MARCH 1, 2001; P.B. BAILEY, W.N. EVERITT AND A. ZETTL
C     VERSION 1.2
C     **********

      PROGRAM DRIVE
C
C     THIS PROGRAM IS AN "INTERACTIVE" DRIVER FOR THE USE OF
C     SLEIGN2. IT CAN BE COMPILED WITH SLEIGN2.F AND XAMPLES.F,
C     FOR EXAMPLE, WHICH MAKES IT EASY TO RUN PROBLEMS WITH ANY
C     ONE OF 32 DIFFERENT DIFFERENTIAL EQUATIONS.  OR, IN PLACE
C     OF XAMPLES.F, ANY OTHER STURM-LIOUVILLE DIFFERENTIAL
C     EQUATION CAN BE EMPLOYED.  THE EASIEST WAY TO CREATE A FILE
C     CONTAINING AN ARBITRARY DIFFERENTIAL EQUATION, IN PLACE OF
C     THE FILE XAMPLES.F, IS BY USING THE INTERACTIVE PROGRAM
C     MAKEPQW.F.
C          IN EITHER CASE, AS SOON AS THE "EXECUTABLE" IS ACTIVATED,
C     THIS PROGRAM DISPLAYS PROMPTS ON THE SCREEN WHICH INVITE THE
C     USER TO SUPPLY, VIA THE KEYBOARD, THE DATA WHICH DEFINE THE
C     PARTICULAR EIGENVALUE PROBLEM WANTED.  DATA SUCH AS
C        THE INTERVAL (a,b),
C        WHETHER THE ENDPOINTS ARE REGULAR, LIMIT CIRCLE, LIMIT
C          POINT, ETC,
C        WHAT KIND OF BOUNDARY CONDITIONS ARE WANTED AT EACH END,
C        WHICH EIGENVALUES ARE WANTED,
C        AND WHETHER OR NOT A PLOT OF AN EIGENFUNCTION IS WANTED.
C        ETC.
C
C     THERE IS ALSO ANOTHER WAY OF USING THIS DRIVER, WHICH AVOIDS
C     THE "QUESTION & ANSWER" FORMAT, IF A USER WOULD PREFER.  IT
C     REQUIRES PUTTING THE PROBLEM DATA (SUCH AS a, b, Regular,
C     Limit Circle, Limit Point, Boundary Conditions, etc.) IN A
C     VERY BRIEF TEXT FILE, CALLED auto.in .  FOR MORE DETAILS
C     ABOUT THIS METHOD, SEE THE COMMENTS AT THE BEGINNING OF THE
C     SUBROUTINE AUTO() WHICH IS IN THE SAME FILE AS THE ONE
C     CONTAINING THIS DRIVER.  THIS "AUTOMATIC" MODE IS, OF COURSE,
C     MUCH FASTER, BUT PROBABLY SHOULD BE USED ONLY BY SOMEONE WHO
C     IS ALREADY FAMILIAR WITH USING THE MORE USUAL "QUESTION &
C     ANSWER" MODE.
C
C     .. Scalars in Common ..
      REAL A,A1,A1S,A2,A2S,AA,ALF,ALFA1,ALFA2,ALPH0,
     +     ASAV,B,B1,B1S,B2,
     +     B2S,BB,BETA0,BETA1,BETA2,BSAV,DTHDAA,DTHDBB,EIG,EIGSAV,
     +     EPSMIN,FA,FB,GAMM0,GQA,GQB,GWA,GWB,H0,HPI,K0,K11,K12,K21,K22,
     +     L0,LPQA,LPQB,LPWA,LPWB,NU0,P0ATA,P0ATAS,P0ATB,P0ATBS,P10,P20,
     +     P30,P40,P50,P60,PI,QFATA,QFATAS,QFATB,QFATBS,SLF9,TMID,TOL,
     +     TOLL,TSAVEL,TSAVER,TWOPI,Z
      INTEGER ILAST,IND,INTAB,INTSAV,ISAVE,ISLFN,JFLAG,MDTHZ,MFS,MLS,
     +        MMWD,NCA,NCB,NEIG1,NEIG2,NEND,NF,NIVP,NLAST,NUMB0,NUMEIG,
     +        NV,T21,T22,T23,T24,T25
      LOGICAL ADDD,LCA,LCB,PEIGF,PR,REGA,REGB,RITE
      CHARACTER*9 INFM,INFP
      CHARACTER*19 FILLA,FILLB,FILLC
      CHARACTER*32 CH1,CH2,CH3,CH4,CH5,CH6
C     ..
C     .. Arrays in Common ..
      REAL EES(50),SLFUN(9),TEE(100),TT(7,2),TTS(50),
     +     YS(200),YY(7,3,2),
     +     ZEE(100)
      INTEGER IIS(50),JAY(100),MMW(100),NT(2)
      CHARACTER*39 BLNK(2),STAR(2),STR(2)
      CHARACTER*55 XC(8)
C     ..
C     .. Local Scalars ..
      REAL ONE
      INTEGER I,IFLAG,NEIGS,NTMP,NUMB,RESP
      LOGICAL EIGV,PERIOD,SKIPB
      CHARACTER*32 TAPE23
C     ..
C     .. Local Arrays ..
      REAL PLOTF(1000,6),XT(1000,2)
C     ..
C     .. External Subroutines ..
      EXTERNAL AUTO,DRAW,DSPLAY,EXAMP,PERIO,QPLOT,SLEIGN,STAGE
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ATAN,MIN
C     ..
C     .. Common blocks ..
      COMMON /ALBE/LPWA,LPQA,FA,GWA,GQA,LPWB,LPQB,FB,GWB,GQB
      COMMON /BCDATA/A1S,A2S,P0ATAS,QFATAS,B1S,B2S,P0ATBS,QFATBS
      COMMON /BCONDS/A1,A2,B1,B2,ALFA1,ALFA2,BETA1,BETA2,NUMEIG,EIG,TOL,
     +       TOLL,NEIG1,NEIG2,ALF,K11,K12,K21,K22
      COMMON /DATADT/ASAV,BSAV,INTSAV
      COMMON /DATAF/EIGSAV,IND
      COMMON /EIGS/EES,TTS,IIS,ILAST,NLAST,JFLAG,SLF9,NIVP,NEND,PEIGF,
     +       SLFUN,RITE,ISLFN,NF,NV
      COMMON /FLAG/NUMB0
      COMMON /LP/MFS,MLS
      COMMON /PAR/NU0,H0,K0,L0,ALPH0,BETA0,GAMM0,P10,P20,P30,P40,P50,P60
      COMMON /PARAM/INTAB,A,NCA,P0ATA,QFATA,B,NCB,P0ATB,QFATB,REGA,LCA,
     +       REGB,LCB
      COMMON /PASS/YS,MMW,MMWD
      COMMON /PIE/PI,TWOPI,HPI
      COMMON /PRIN/PR,T21
      COMMON /RNDOFF/EPSMIN
      COMMON /SHAR/INFM,INFP,CH1,CH2,CH3,CH4,CH5,CH6,STAR,BLNK,STR,
     +       FILLA,FILLB,FILLC,XC
      COMMON /TAPES/T22,T23,T24,T25
      COMMON /TDATA/AA,TMID,BB,DTHDAA,DTHDBB,MDTHZ,ADDD
      COMMON /TEEZ/TEE
      COMMON /TEMP/TT,YY,NT
      COMMON /TSAVE/TSAVEL,TSAVER,ISAVE
      COMMON /Z1/Z
      COMMON /ZEEZ/JAY,ZEE
C     ..
      ONE = 1.0D0
      HPI = 2.0D0*ATAN(ONE)
      PI = 2.0D0*HPI
      TWOPI = 2.0D0*PI
C
      T21 = 21
      T22 = 22
      T23 = 23
      T24 = 24
      T25 = 25
C
C         DEFINITIONS OF SOME STRINGS.
C      (TO BE USED IN SUBROUTINE DSPLAY.)
C
      INFP = '+INFINITY'
      INFM = '-INFINITY'
      CH1 = 'REGULAR                       * '
      CH2 = 'WEAKLY REGULAR                * '
      CH3 = 'LIMIT CIRCLE, NON-OSCILLATORY * '
      CH4 = 'LIMIT CIRCLE, OSCILLATORY     * '
      CH5 = 'LIMIT POINT                   * '
      CH6 = 'UNSPEC.(NOT LCO), DEFAULT B.C.* '
      STAR(1) = ' **************************************'
      STAR(2) = '***************************************'
      BLNK(1) = ' *                                     '
      BLNK(2) = '                                      *'
      FILLA = '*******************'
      FILLB = '                  *'
      FILLC = '------------------ '
C
      XC(1) = ' * (1) THE SOLUTION Y                               * '
      XC(2) = ' * (2) THE QUASI-DERIVATIVE p*Y''                    *'
      XC(3) = ' * (3) THE BOUNDARY CONDITION FUNCTION Y OR [Y,U]   * '
      XC(4) = ' * (4) THE BOUNDARY CONDITION FUNCTION p*Y'' OR [Y,V]*'
      XC(5) = ' * (5) THE PRUFER ANGLE, THETA                      * '
      XC(6) = ' * (6) THE PRUFER MODULUS, RHO                      * '
      XC(7) = ' * (1) x IN THE INTERVAL (a,b)                      * '
      XC(8) = ' * (2) t IN THE INTERVAL (-1,1)                     * '
C
C
      OPEN (T21,FILE='test.out')
C
C     INTRODUCTION :-
      CALL DSPLAY(1,RESP)
      IF (RESP.EQ.0) CALL AUTO

      CALL EXAMP

C     IS MORE INFORMATION REQUIRED? :-
      CALL DSPLAY(2,RESP)

C     SHOULD THE RESULTS BE RECORDED IN A FILE? :=
      CALL DSPLAY(3,RESP)
C     IF SO, GET THE NAME OF THE FILE :-
      IF (RESP.EQ.1) CALL DSPLAY(4,RESP)

  100 CONTINUE
      SKIPB = .FALSE.
      P0ATA = -1.0D0
      QFATA = 1.0D0
      P0ATB = -1.0D0
      QFATB = 1.0D0
C     WHAT KIND OF INTERVAL IS THE PROBLEM ON? :-
      CALL DSPLAY(5,RESP)

  120 CONTINUE
C     GET ENDPOINT A, IF FINITE :-
      IF (INTAB.EQ.1 .OR. INTAB.EQ.2) CALL DSPLAY(6,RESP)

C     CLASSIFICATION OF A :-
      CALL DSPLAY(7,RESP)

C       (IS P0 AT A, OR QF AT A? :- )
      IF (.NOT.REGA .AND. INTAB.LE.2) CALL DSPLAY(8,RESP)

      IF (SKIPB) GO TO 150
  140 CONTINUE
C     (GET ENDPOINT B, IF FINITE :- )
      IF (INTAB.EQ.1 .OR. INTAB.EQ.3) CALL DSPLAY(9,RESP)

C     (CLASSIFICATION OF B :- )
      CALL DSPLAY(10,RESP)

C       (IS P0 AT B, OR QF AT B :- )
      IF (.NOT.REGB .AND. (INTAB.EQ.1.OR.INTAB.EQ.3)) CALL DSPLAY(11,
     +    RESP)

  150 CONTINUE
C     BRIEF SUMMARY OF PROBLEM PARAMETERS :-
C     (IS THIS CORRECT SO FAR?  :- )
      CALL DSPLAY(12,RESP)

C     (IF THIS IS NOT THE RIGHT PROBLEM, DO IT OVER :- )
      IF (RESP.NE.1) THEN
          CALL DSPLAY(13,RESP)
          IF (RESP.EQ.1) THEN
              SKIPB = .TRUE.
              GO TO 120
          ELSE IF (RESP.EQ.2) THEN
              GO TO 140
          ELSE IF (RESP.EQ.3) THEN
              GO TO 100
          END IF
      END IF

C ----------------------------------------------------C
C
C     AT THIS POINT THE DIFFERENTIAL EQUATION AND THE INTERVAL OF
C       INTEREST HAVE BEEN DEFINED AND CHARACTERIZED.
C

C     SHOULD WE COMPUTE AN EIVENVALUE? OR THE SOLUTION
C       TO AN INITIAL VALUE PROBLEM? :-
      CALL DSPLAY(14,RESP)
      EIGV = (RESP.EQ.1)
      IF (EIGV) THEN
          NIVP = 0
          NEND = 0
          PERIOD = .FALSE.
      END IF

  200 CONTINUE
      IF (EIGV .AND. NCA.LE.4 .AND. NCB.LE.4) THEN
C        (ARE THE BC'S SEPARATED?, OR COUPLED?)
          CALL DSPLAY(15,RESP)
          IF (RESP.EQ.1) THEN
C           (THIS MEANS THE BC'S ARE COUPLED)
              PERIOD = .TRUE.
          ELSE
              PERIOD = .FALSE.
          END IF
      END IF
C

  300 CONTINUE
C       SET BOUNDARY CONDITIONS, OR INITIAL CONDITIONS
      IF (EIGV) THEN
          IF (PERIOD) THEN
C           (SET COUPLED BC'S :- )
              CALL DSPLAY(18,RESP)

          ELSE
C           (SET SEPARATED BC'S AT A :- )
              CALL DSPLAY(16,RESP)

C           (SET SEPARATED BC'S AT B :- )
              CALL DSPLAY(17,RESP)
          END IF
      ELSE
C        (SET INITIAL CONDITIONS FOR IVP :- )
          PERIOD = .FALSE.
          CALL DSPLAY(19,RESP)
      END IF
C ----------------------------------------------------C
C
C     PRESUMABLY WE NOW HAVE ASSEMBLED ALL THE DATA
C     NEEDED FOR DEFINING
C          THE EIGENVALUE PROBLEM,
C                     OR
C          THE INITIAL VALUE PROBLEM,
C     WHICHEVER IS WANTED.

  400 CONTINUE
C          COMPUTE A SINGLE EIGENVALUE?, OR SERIES? :-
      NEIGS = 0
      IF (EIGV .AND. .NOT.PERIOD) THEN
          CALL DSPLAY(20,RESP)
          IF (RESP.EQ.1) THEN
              CALL DSPLAY(21,RESP)
              NEIGS = 1
          ELSE
              CALL DSPLAY(23,RESP)
              NEIGS = NEIG2 - NEIG1 + 1
          END IF
      ELSE IF (EIGV) THEN
          CALL DSPLAY(24,RESP)
          IF (RESP.EQ.1) THEN
              CALL DSPLAY(25,RESP)
              NEIGS = 1
          ELSE
              CALL DSPLAY(27,RESP)
              NEIGS = NEIG2 - NEIG1 + 1
              EIG = 0.0D0
          END IF
      END IF
C
      IF (EIGV) THEN
C        (COMPUTE REQUESTED EIGENVALUES)

          IF (NEIGS.EQ.1) THEN
              NEIG1 = NUMEIG
              NEIG2 = NUMEIG
              NLAST = NUMEIG
          END IF
          ILAST = NEIG2 - NEIG1 + 1

          IF (.NOT.PERIOD) THEN
              DO 410 I = 1,ILAST
                  NUMB = NEIG1 + I - 1
                  TOL = TOLL
                  NTMP = NUMB
                  IF (NEIGS.GT.1) EIG = 0.0D0
                  PEIGF = .FALSE.

                  CALL SLEIGN(A,B,INTAB,P0ATA,QFATA,P0ATB,QFATB,A1,A2,
     +                        B1,B2,NTMP,EIG,TOL,IFLAG,0,SLFUN,NCA,NCB)

                  IF (IFLAG.EQ.0 .OR. IFLAG.GE.15) THEN
                      IF (IFLAG.EQ.0) THEN
                          WRITE (*,FMT=*) ' IMPROPER INPUT PARAMETERS '
                          IF (RITE) THEN
                              WRITE (T22,FMT=*)
     +                          ' IMPROPER INPUT PARAMETERS '
                              WRITE (T22,FMT=*)
     +                        '----------------------------------------'
                          END IF
                      ELSE IF (IFLAG.EQ.15) THEN
                          WRITE (*,FMT=*)
     +                      ' WE CANNOT HANDLE THIS KIND OF ENDPOINT '
                          IF (RITE) THEN
                              WRITE (T22,FMT=*)
     +                        ' WE CANNOT HANDLE THIS KIND OF ENDPOINT '
                              WRITE (T22,FMT=*)
     +                        '----------------------------------------'
                          END IF
                      ELSE IF (IFLAG.EQ.16) THEN
                          WRITE (*,FMT=*) ' COULD NOT GET STARTED '
                          IF (RITE) THEN
                              WRITE (T22,FMT=*)
     +                          ' COULD NOT GET STARTED '
                              WRITE (T22,FMT=*)
     +                        '----------------------------------------'
                          END IF
                      ELSE IF (IFLAG.EQ.17) THEN
                          WRITE (*,FMT=*) ' FAILED TO GET A BRACKET '
                          IF (RITE) THEN
                              WRITE (T22,FMT=*)
     +                          ' FAILED TO GET A BRACKET '
                              WRITE (T22,FMT=*)
     +                        '----------------------------------------'
                          END IF
                      ELSE IF (IFLAG.EQ.18) THEN
                          WRITE (*,FMT=*) ' ESTIMATOR FAILED '
                          IF (RITE) THEN
                              WRITE (T22,FMT=*) ' ESTIMATOR FAILED '
                              WRITE (T22,FMT=*)
     +                        '----------------------------------------'
                          END IF
                      END IF
                      GO TO 700
                  END IF

                  IFLAG = MIN(IFLAG,4)
                  JFLAG = 0

C   JFLAG = 1 OR 2 MEANS ONE OF THE END-POINTS IS LP .
C   JFLAG = 1 MEANS THERE IS NO CONTINUOUS SPECTRUM .
C         = 2 MEANS THERE IS A CONTINUOUS SPECTRUM .

                  SLF9 = SLFUN(9)
                  IF (SLF9.GT.-10000.0D0) JFLAG = 1
                  IF (SLF9.LT.10000.0D0 .AND. JFLAG.EQ.1) JFLAG = 2

                  EES(I) = EIG
                  TTS(I) = TOL
                  IIS(I) = IFLAG

                  IF (IFLAG.EQ.3) THEN
C  (THIS MEANS THAT THERE IS NO EIGENVALUE )
C  (  WITH THIS INDEX                      )
                      ILAST = I
                      NLAST = NTMP
                      GO TO 420
                  END IF
  410         CONTINUE
  420         CONTINUE
              CALL DSPLAY(28,RESP)

          ELSE

              DO 430 I = 1,ILAST
                  NUMEIG = NEIG1 + I - 1
                  TOL = TOLL
                  EIG = 0.0D0
                  CALL PERIO(A,B,INTAB,P0ATA,QFATA,P0ATB,QFATB,A1,A2,B1,
     +                       B2,NUMEIG,EIG,TOL,IFLAG,SLFUN,NCA,NCB,ALF,
     +                       K11,K12,K21,K22)
                  EES(I) = EIG
                  TTS(I) = TOL
                  IIS(I) = IFLAG
C      (IFLAG = 2 MEANS THAT THE EIGENVALUE IS A DOUBLE)
  430         CONTINUE
              CALL DSPLAY(29,RESP)

          END IF
      END IF
C                                                              C
C -------------------------------------------------------------C
  500 CONTINUE
C
C  WE SHOULD NOW HAVE THE EIGENVALUES COMPUTED, OR HAVE        C
C  THE INITIAL VALUE PROBLEM DEFINED, AS THE CASE MAY BE.      C
C  IF ONLY ONE EIGENVALUE HAS BEEN COMPUTED, THE CORRESPONDING C
C  EIGENFUNCTION CAN NOW BE PLOTTED.  OR IF AN INITIAL VALUE   C
C  PROBLEM HAS BEEN DEFINED, IT MAY BE PLOTTED NOW.            C
C
      IF (.NOT.EIGV) THEN
          CALL DSPLAY(32,RESP)
          CALL STAGE(A,B,INTAB,P0ATA,QFATA,P0ATB,QFATB,A1,A2,B1,B2,EIG,
     +               IFLAG,SLFUN,NCA,NCB,NIVP,NEND)
      ELSE IF (NEIGS.EQ.1 .AND. (IFLAG.EQ.1.OR.IFLAG.EQ.2)) THEN
          IF (.NOT.PERIOD) THEN
              CALL DSPLAY(22,RESP)
          ELSE
              CALL DSPLAY(26,RESP)
          END IF
      ELSE IF (IFLAG.EQ.0) THEN
          PEIGF = .FALSE.
      END IF

      IF (PEIGF) THEN
          CALL DRAW(A1,A2,B1,B2,NUMEIG,EIG,SLFUN,NIVP,NEND,EIGV,NCA,NCB,
     +              ISLFN,XT,PLOTF,K11,K12,K21,K22,PERIOD)
  600     CONTINUE
C         PLOT WHICH FUNCTION ?
          CALL DSPLAY(33,RESP)
          IF (RESP.EQ.1) CALL QPLOT(ISLFN,XT,NV,PLOTF,NF)
C         SAVE THE PLOT FILE ?
          CALL DSPLAY(34,RESP)
          IF (RESP.EQ.1) THEN
              WRITE (*,FMT=*) ' SPECIFY NAME OF FILE FOR PLOTTING '
              READ (*,FMT='(A)') TAPE23
              OPEN (T23,FILE=TAPE23,STATUS='NEW')
              DO 610 I = 1,ISLFN
                  WRITE (T23,FMT=*) XT(9+I,NV),PLOTF(9+I,NF)
  610         CONTINUE
              CLOSE (T23)
              WRITE (*,FMT=*) ' THE PLOT FILE HAS BEEN WRITTEN TO ',
     +          TAPE23
              IF (RITE) WRITE (T22,FMT=*)
     +            ' THE PLOT FILE HAS BEEN WRITTEN TO ',TAPE23
          END IF
C         PLOT ANOTHER FUNCTION ?
          CALL DSPLAY(35,RESP)
          IF (RESP.EQ.1) GO TO 600
      END IF

  700 CONTINUE
C  ARE WE FINISHED?, OR DO WE HAVE MORE PROBLEMS TO DO?         C
C
      IF (EIGV) THEN
          CALL DSPLAY(30,RESP)
          IF (RESP.EQ.1) GO TO 400
          IF (RESP.EQ.2) GO TO 200
          IF (RESP.EQ.3) GO TO 100
          EIGV = .FALSE.
          IF (RESP.EQ.4) GO TO 300
      ELSE
          CALL DSPLAY(31,RESP)
          IF (RESP.EQ.1) GO TO 500
          IF (RESP.EQ.2) GO TO 300
          IF (RESP.EQ.3) GO TO 100
          EIGV = .TRUE.
          IF (RESP.EQ.4) GO TO 200
      END IF
C
      CLOSE (T22)
      CLOSE (T21)
      STOP
      END
C
      SUBROUTINE DSPLAY(DIS,RESP)
C     .. Scalars in Common ..
      REAL A,A1,A2,ALF,ALFA1,ALFA2,B,B1,B2,BETA1,BETA2,
     +     EIG,HPI,K11,K12,
     +     K21,K22,P0ATA,P0ATB,PI,QFATA,QFATB,SLF9,TOL,TOLL,TWOPI
      INTEGER ILAST,INTAB,ISLFN,JFLAG,NCA,NCB,NEIG1,NEIG2,NEND,NF,NIVP,
     +        NLAST,NUMEIG,NV,T22,T23,T24,T25
      LOGICAL LCA,LCB,PEIGF,REGA,REGB,RITE
      CHARACTER*9 INFM,INFP
      CHARACTER*19 FILLA,FILLB,FILLC
      CHARACTER*32 CH1,CH2,CH3,CH4,CH5,CH6
C     ..
C     .. Local Scalars ..
      REAL CC,DETK,R1,R2,RHO,THA,THB,TMP
      INTEGER I,I1,I2,IFLAG,NANS
      LOGICAL DUBBLE,YEH
      CHARACTER ANSCH,HQ,YN
      CHARACTER*32 CHA,CHANS,CHB,FMAT,TAPE22
      CHARACTER*70 CHTXT
C     ..
C     .. Local Arrays ..
      INTEGER ICOL(2)
      CHARACTER*2 COL(32)
C     ..
C     .. External Subroutines ..
      EXTERNAL HELP,LST,LSTDIR
C     ..
C     .. External Functions ..
      CHARACTER*32 FMT2
      EXTERNAL FMT2
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS,ATAN2,COS,SIN,SQRT
C     ..

C
C     .. Scalar Arguments ..
      INTEGER DIS,RESP
C     ..
C     .. Arrays in Common ..
      REAL EES(50),SLFUN(9),TTS(50)
      INTEGER IIS(50)
      CHARACTER*39 BLNK(2),STAR(2),STR(2)
      CHARACTER*55 XC(8)
C     ..
C     .. Common blocks ..
      COMMON /BCONDS/A1,A2,B1,B2,ALFA1,ALFA2,BETA1,BETA2,NUMEIG,EIG,TOL,
     +       TOLL,NEIG1,NEIG2,ALF,K11,K12,K21,K22
      COMMON /EIGS/EES,TTS,IIS,ILAST,NLAST,JFLAG,SLF9,NIVP,NEND,PEIGF,
     +       SLFUN,RITE,ISLFN,NF,NV
      COMMON /PARAM/INTAB,A,NCA,P0ATA,QFATA,B,NCB,P0ATB,QFATB,REGA,LCA,
     +       REGB,LCB
      COMMON /PIE/PI,TWOPI,HPI
      COMMON /SHAR/INFM,INFP,CH1,CH2,CH3,CH4,CH5,CH6,STAR,BLNK,STR,
     +       FILLA,FILLB,FILLC,XC
      COMMON /TAPES/T22,T23,T24,T25
C     ..
C     .. Data statements ..
      DATA COL/'01','02','03','04','05','06','07','08','09','10','11',
     +     '12','13','14','15','16','17','18','19','20','21','22','23',
     +     '24','25','26','27','28','29','30','31','32'/
C     ..
      RESP = 1
      GO TO (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,
     +       23,24,25,26,27,28,29,30,31,32,33,34,35) DIS
C
    1 CONTINUE
      WRITE (*,FMT=*)
      WRITE (*,FMT=*) ' This program solves the boundary value problem '
      WRITE (*,FMT=*) '   defined by the differential equation '
      WRITE (*,FMT=*)
      WRITE (*,FMT=*) '   -(py'')'' + q*y = lambda*w*y  '
      WRITE (*,FMT=*)
      WRITE (*,FMT=*) ' together with appropriate boundary conditions. '
      WRITE (*,FMT=*)
      WRITE (*,FMT=*)
     +  ' HELP may be called at any point where the program   '
      WRITE (*,FMT=*)
     +  '    halts and displays (h?) by pressing "h <ENTER>". '
      WRITE (*,FMT=*)
     +  '    To RETURN from HELP, press "r <ENTER>".          '
      WRITE (*,FMT=*)
     +  '    To QUIT at any program halt, press "q <ENTER>".  '
      WRITE (*,FMT=*)
     +  ' WOULD YOU LIKE AN OVERVIEW OF HELP ? (Y/N) (h?) '
      READ (*,FMT=9010) CHANS
      READ (CHANS,FMT=9020) HQ
      IF (HQ.EQ.'q' .OR. HQ.EQ.'Q') STOP
      IF (HQ.EQ.'h' .OR. HQ.EQ.'H' .OR. HQ.EQ.'y' .OR.
     +    HQ.EQ.'Y') CALL HELP(1)
      WRITE (*,FMT=*)
      YEH = HQ .EQ. 'y' .OR. HQ .EQ. 'Y' .OR. HQ .EQ. 'n' .OR.
     +      HQ .EQ. 'N' .OR. HQ .EQ. 'a' .OR. HQ .EQ. 'A'
      IF (.NOT.YEH) GO TO 1
      IF (HQ.EQ.'a' .OR. HQ.EQ.'A') RESP = 0
      RETURN
C
    2 CONTINUE
      WRITE (*,FMT=*) ' DO YOU REQUIRE INFORMATION ON THE RANGE OF '
      WRITE (*,FMT=*) '    BOUNDARY CONDITIONS AVAILABLE ? (Y/N) (h?) '
      READ (*,FMT=9010) CHANS
      READ (CHANS,FMT=9020) HQ
      IF (HQ.EQ.'q' .OR. HQ.EQ.'Q') STOP
      IF (HQ.EQ.'h' .OR. HQ.EQ.'H' .OR. HQ.EQ.'y' .OR.
     +    HQ.EQ.'Y') CALL HELP(7)
      YEH = HQ .EQ. 'y' .OR. HQ .EQ. 'Y' .OR. HQ .EQ. 'n' .OR.
     +      HQ .EQ. 'N' .OR. HQ .EQ. 'h' .OR. HQ .EQ. 'H'
      IF (.NOT.YEH) GO TO 2
      IF (HQ.EQ.'h' .OR. HQ.EQ.'H') CALL HELP(7)
      RETURN
C
    3 CONTINUE
  105 CONTINUE
      WRITE (*,FMT=*) ' DO YOU WANT A RECORD KEPT OF THE PROBLEMS '
      WRITE (*,FMT=*) '    AND RESULTS ? (Y/N) (h?) '
      READ (*,FMT=9010) CHANS
      READ (CHANS,FMT=9020) HQ
      IF (HQ.EQ.'q' .OR. HQ.EQ.'Q') STOP
      YEH = HQ .EQ. 'y' .OR. HQ .EQ. 'Y' .OR. HQ .EQ. 'n' .OR.
     +      HQ .EQ. 'N' .OR. HQ .EQ. 'h' .OR. HQ .EQ. 'H'
      IF (.NOT.YEH) GO TO 105
      IF (HQ.EQ.'h' .OR. HQ.EQ.'H') THEN
          CALL HELP(8)
          GO TO 105
      END IF
      RITE = HQ .EQ. 'y' .OR. HQ .EQ. 'Y' .OR. HQ .EQ. 'yes' .OR.
     +       HQ .EQ. 'YES'
      IF (.NOT.RITE) RESP = 0
      RETURN
C
    4 CONTINUE
  205 CONTINUE
      WRITE (*,FMT=*) '  SPECIFY NAME OF THE OUTPUT RECORD FILE: (h?) '
      READ (*,FMT=9010) CHANS
      IF (CHANS.EQ.'q' .OR. CHANS.EQ.'Q') THEN
          STOP
      ELSE IF (CHANS.EQ.'h' .OR. CHANS.EQ.'H') THEN
          CALL HELP(8)
          GO TO 205
      ELSE
          TAPE22 = CHANS
          WRITE (*,FMT=*) '  YOU MAY ENTER SOME HEADER LINE FOR THE '
          WRITE (*,FMT=*) '  OUTPUT RECORD FILE (<=70 CHARACTERS)   '
          WRITE (*,FMT=*) '  IF YOU WISH; OTHERWISE JUST HIT "RETURN".'
          WRITE (*,FMT=*)
          READ (*,FMT=9030) CHTXT
      END IF

C
      OPEN (T22,FILE=TAPE22,STATUS='NEW')
C
      WRITE (T22,FMT=*) '                         ',TAPE22
      WRITE (T22,FMT=*)
      WRITE (T22,FMT=*) CHTXT
      WRITE (T22,FMT=*)
C
      RETURN
C
    5 CONTINUE
  405 CONTINUE
      WRITE (*,FMT=*)
     +  ' ************************************************** '
      WRITE (*,FMT=*)
     +  ' * INDICATE THE KIND OF PROBLEM INTERVAL (a,b):   * '
      WRITE (*,FMT=*)
     +  ' *                                                * '
      WRITE (*,FMT=*)
     +  ' *  (CHECK THAT THE COEFFICIENTS p,q,w ARE WELL   * '
      WRITE (*,FMT=*)
     +  ' *  DEFINED THROUGHOUT THE INTERVAL open(a,b).)   * '
      WRITE (*,FMT=*)
     +  ' *                                                * '
      WRITE (*,FMT=*)
     +  ' *    (1) FINITE, (a,b)                           * '
      WRITE (*,FMT=*)
     +  ' *                                                * '
      WRITE (*,FMT=*)
     +  ' *    (2) SEMI-INFINITE, (a,+INFINITY)            * '
      WRITE (*,FMT=*)
     +  ' *                                                * '
      WRITE (*,FMT=*)
     +  ' *    (3) SEMI-INFINITE, (-INFINITY,b)            * '
      WRITE (*,FMT=*)
     +  ' *                                                * '
      WRITE (*,FMT=*)
     +  ' *    (4) DOUBLY INFINITE, (-INFINITY,+INFINITY)  * '
      WRITE (*,FMT=*)
     +  ' *                                                * '
      WRITE (*,FMT=*)
     +  ' * ENTER THE NUMBER OF YOUR CHOICE: (h?)          * '
      WRITE (*,FMT=*)
     +  ' ************************************************** '
      WRITE (*,FMT=*)
      READ (*,FMT=9010) CHANS
      READ (CHANS,FMT=9020) HQ
      IF (HQ.EQ.'q' .OR. HQ.EQ.'Q') STOP
      IF (HQ.EQ.'h' .OR. HQ.EQ.'H') THEN
          CALL HELP(9)
          GO TO 405
      END IF
      YEH = HQ .EQ. '1' .OR. HQ .EQ. '2' .OR. HQ .EQ. '3' .OR.
     +      HQ .EQ. '4'
      IF (.NOT.YEH) GO TO 405
      READ (CHANS,FMT='(I32)') INTAB
      IF (RITE) THEN
          IF (INTAB.EQ.1) WRITE (T22,FMT=*) ' The interval is (a,b) .'
          IF (INTAB.EQ.2) WRITE (T22,FMT=*) ' The interval is (a,+inf).'
          IF (INTAB.EQ.3) WRITE (T22,FMT=*
     +        ) ' The interval is (-inf,b). '
          IF (INTAB.EQ.4) WRITE (T22,FMT=*
     +        ) ' The interval is (-inf,+inf).'
      END IF
      RESP = INTAB
      RETURN
C
    6 CONTINUE
   60 CONTINUE
      WRITE (*,FMT=*)
      WRITE (*,FMT=*) ' ********************************************* '
      WRITE (*,FMT=*) ' * INPUT a: (h?)                             * '
      WRITE (*,FMT=*) ' ********************************************* '
      WRITE (*,FMT=*) ' a = '
      READ (*,FMT=9010) CHANS
      READ (CHANS,FMT=9020) HQ
      IF (HQ.EQ.'q' .OR. HQ.EQ.'Q') STOP
      IF (HQ.EQ.'h' .OR. HQ.EQ.'H') THEN
          CALL HELP(10)
          GO TO 60
      END IF
      CALL LST(CHANS,A)
      IF (RITE) WRITE (T22,FMT=*) ' a = ',A
      IF (INTAB.EQ.2) B = A + 1.0D0

      WRITE (*,FMT=*)
      WRITE (*,FMT=*)
      WRITE (*,FMT=*)
      WRITE (*,FMT=*)
      WRITE (*,FMT=*)
      RETURN
C
    7 CONTINUE
   70 CONTINUE
      STR(1) = ' * IS THIS PROBLEM:                    '
      STR(2) = '                                      *'
      WRITE (*,FMT=9110) STAR
      WRITE (*,FMT=9110) STR
      WRITE (*,FMT=9110) BLNK
      STR(1) = ' *   (1) REGULAR AT a ?                '
      STR(2) = '                                      *'
      WRITE (*,FMT=9110) STR
      STR(1) = ' *       (I.E., THE FUNCTIONS p, q, & w'
      STR(2) = ' ARE BOUNDED CONTINUOUS NEAR a;       *'
      WRITE (*,FMT=9110) STR
      STR(1) = ' *        p & w ARE POSITIVE AT a.)    '
      STR(2) = '                                      *'
      WRITE (*,FMT=9110) STR
      WRITE (*,FMT=9110) BLNK
      STR(1) = ' *   (2) WEAKLY REGULAR AT a ?         '
      STR(2) = '                                      *'
      WRITE (*,FMT=9110) STR
      STR(1) = ' *       (I.E., THE FUNCTIONS 1/p, q, &'
      STR(2) = ' w ALL ARE FINITELY INTEGRABLE ON     *'
      WRITE (*,FMT=9110) STR
      STR(1) = ' *        SOME INTERVAL [a,a+e] FOR e >'
      STR(2) = ' 0; p & w ARE POSITIVE NEAR a.)       *'
      WRITE (*,FMT=9110) STR
      WRITE (*,FMT=9110) BLNK
      STR(1) = ' *   (3) LIMIT CIRCLE, NON-OSCILLATORY '
      STR(2) = 'AT a ?                                *'
      WRITE (*,FMT=9110) STR
      WRITE (*,FMT=9110) BLNK
      STR(1) = ' *   (4) LIMIT CIRCLE, OSCILLATORY AT a'
      STR(2) = ' ?                                    *'
      WRITE (*,FMT=9110) STR
      WRITE (*,FMT=9110) BLNK
      STR(1) = ' *   (5) LIMIT POINT AT a ?            '
      STR(2) = '                                      *'
      WRITE (*,FMT=9110) STR
      WRITE (*,FMT=9110) BLNK
      STR(1) = ' *   (6) NOT SPECIFIED (BUT NOT LIMIT C'
      STR(2) = 'IRCLE OSCILLATORY) WITH DEFAULT       *'
      WRITE (*,FMT=9110) STR
      STR(1) = ' *       BOUNDARY CONDITION AT a ?     '
      STR(2) = '                                      *'
      WRITE (*,FMT=9110) STR
      WRITE (*,FMT=9110) BLNK
      STR(1) = ' * ENTER THE NUMBER OF YOUR CHOICE: (h)'
      STR(2) = '?                                     *'
      WRITE (*,FMT=9110) STR
      WRITE (*,FMT=9110) STAR
      WRITE (*,FMT=*)
C
C     SPECIFY TYPE OF BOUNDARY CONDITION AT a.
C
      READ (*,FMT=9010) CHANS
      READ (CHANS,FMT=9020) HQ
      IF (HQ.EQ.'q' .OR. HQ.EQ.'Q') STOP
      IF (HQ.EQ.'h' .OR. HQ.EQ.'H') THEN
          CALL HELP(4)
          GO TO 70
      END IF
      YEH = HQ .EQ. '1' .OR. HQ .EQ. '2' .OR. HQ .EQ. '3' .OR.
     +      HQ .EQ. '4' .OR. HQ .EQ. '5' .OR. HQ .EQ. '6'
      IF (.NOT.YEH) GO TO 70
      READ (CHANS,FMT='(I32)') NANS
C
C     SET CHARACTER STRING CHA ACCORDING TO BOUNDARY CONDITION AT a.
C
      IF (NANS.EQ.1) THEN
          REGA = .TRUE.
          NCA = 1
          P0ATA = -1.0D0
          QFATA = 1.0D0
          CHA = CH1
          IF (RITE) WRITE (T22,FMT=*) ' Endpoint a is Regular. '
      ELSE IF (NANS.EQ.2) THEN
          NCA = 2
          CHA = CH2
          IF (RITE) WRITE (T22,FMT=*) ' Endpoint a is Weakly Regular. '
      ELSE IF (NANS.EQ.3) THEN
          LCA = .TRUE.
          CHA = CH3
          NCA = 3
          IF (RITE) WRITE (T22,FMT=*) ' Endpoint a is Limit Circle,',
     +        ' Non-Oscillatory. '
      ELSE IF (NANS.EQ.4) THEN
          LCA = .TRUE.
          CHA = CH4
          NCA = 4
          IF (RITE) WRITE (T22,FMT=*) ' Endpoint a is Limit Circle,',
     +        ' Oscillatory. '
      ELSE IF (NANS.EQ.5) THEN
          CHA = CH5
          NCA = 5
          IF (RITE) WRITE (T22,FMT=*) ' Endpoint a is Limit Point. '
      ELSE
          CHA = CH6
          NCA = 6
          IF (RITE) WRITE (T22,FMT=*)
     +        ' Endpoint a is Singular, unspecified. '
      END IF
      RESP = NANS
      WRITE (*,FMT=*)
      WRITE (*,FMT=*)
      WRITE (*,FMT=*)
      WRITE (*,FMT=*)
      RETURN
C
    8 CONTINUE
   80 CONTINUE
      WRITE (*,FMT=*)
      WRITE (*,FMT=*) ' ********************************************* '
      WRITE (*,FMT=*) ' * IS p = 0. AT a ? (Y/N) (h?)               * '
      WRITE (*,FMT=*) ' ********************************************* '
      WRITE (*,FMT=*)
      READ (*,FMT=9010) CHANS
      READ (CHANS,FMT=9020) HQ
      IF (HQ.EQ.'q' .OR. HQ.EQ.'Q') STOP
      IF (HQ.EQ.'h' .OR. HQ.EQ.'H') THEN
          CALL HELP(11)
          GO TO 80
      END IF
      YEH = HQ .EQ. 'y' .OR. HQ .EQ. 'Y' .OR. HQ .EQ. 'n' .OR.
     +      HQ .EQ. 'N'
      IF (.NOT.YEH) GO TO 80
      YEH = HQ .EQ. 'y' .OR. HQ .EQ. 'Y'
      P0ATA = -1.0D0
      IF (YEH) P0ATA = 1.0D0
      IF (RITE) THEN
          IF (YEH) THEN
              WRITE (T22,FMT=*) ' p is zero at a. '
          ELSE
              WRITE (T22,FMT=*) ' p is not zero at a. '
          END IF
      END IF

   90 CONTINUE
      WRITE (*,FMT=*)
      WRITE (*,FMT=*) ' ********************************************* '
      WRITE (*,FMT=*) ' * IS EITHER OF THE COEFFICIENT FUNCTIONS    * '
      WRITE (*,FMT=*) ' *                                           * '
      WRITE (*,FMT=*) ' *  q OR w UNBOUNDED NEAR a ? (Y/N) (h?)     * '
      WRITE (*,FMT=*) ' ********************************************* '
      WRITE (*,FMT=*)
      READ (*,FMT=9010) CHANS
      READ (CHANS,FMT=9020) HQ
      IF (HQ.EQ.'q' .OR. HQ.EQ.'Q') STOP
      IF (HQ.EQ.'h' .OR. HQ.EQ.'H') THEN
          CALL HELP(11)
          GO TO 90
      END IF
      YEH = HQ .EQ. 'y' .OR. HQ .EQ. 'Y' .OR. HQ .EQ. 'n' .OR.
     +      HQ .EQ. 'N'
      IF (.NOT.YEH) GO TO 90
      YEH = HQ .EQ. 'y' .OR. HQ .EQ. 'Y'
      QFATA = 1.0D0
      IF (YEH) QFATA = -1.0D0
      IF (RITE) THEN
          IF (YEH) THEN
              WRITE (T22,FMT=*) ' either q or w is unbounded near a. '
          ELSE
              WRITE (T22,FMT=*) ' both q and w are bounded near a.  '
          END IF
      END IF
      IF (.NOT.YEH) RESP = 0
      RETURN
C
    9 CONTINUE
  110 CONTINUE
      WRITE (*,FMT=*)
      WRITE (*,FMT=*) ' ********************************************* '
      WRITE (*,FMT=*) ' * INPUT b: (h?)                             * '
      WRITE (*,FMT=*) ' ********************************************* '
      WRITE (*,FMT=*) ' b = '
      READ (*,FMT=9010) CHANS
      READ (CHANS,FMT=9020) HQ
      IF (HQ.EQ.'q' .OR. HQ.EQ.'Q') STOP
      IF (HQ.EQ.'h' .OR. HQ.EQ.'H') THEN
          CALL HELP(10)
          GO TO 110
      END IF
      CALL LST(CHANS,B)
      IF (RITE) WRITE (T22,FMT=*) ' b = ',B
      IF (INTAB.EQ.3) A = B - 1.0D0
      WRITE (*,FMT=*)
      WRITE (*,FMT=*)
      WRITE (*,FMT=*)
      WRITE (*,FMT=*)
      WRITE (*,FMT=*)
      RETURN
C
   10 CONTINUE
  120 CONTINUE
      STR(1) = ' * IS THIS PROBLEM:                    '
      STR(2) = '                                      *'
      WRITE (*,FMT=9110) STAR
      WRITE (*,FMT=9110) STR
      WRITE (*,FMT=9110) BLNK
      STR(1) = ' *   (1) REGULAR AT b ?                '
      STR(2) = '                                      *'
      WRITE (*,FMT=9110) STR
      STR(1) = ' *       (I.E., THE FUNCTIONS p, q, & w'
      STR(2) = ' ARE BOUNDED CONTINUOUS NEAR b;       *'
      WRITE (*,FMT=9110) STR
      STR(1) = ' *        p & w ARE POSITIVE AT b.)    '
      STR(2) = '                                      *'
      WRITE (*,FMT=9110) STR
      WRITE (*,FMT=9110) BLNK
      STR(1) = ' *   (2) WEAKLY REGULAR AT b ?         '
      STR(2) = '                                      *'
      WRITE (*,FMT=9110) STR
      STR(1) = ' *       (I.E., THE FUNCTIONS 1/p, q, &'
      STR(2) = ' w ALL ARE FINITELY INTEGRABLE ON     *'
      WRITE (*,FMT=9110) STR
      STR(1) = ' *        SOME INTERVAL [b-e,b] FOR e >'
      STR(2) = ' 0; p & w ARE POSITIVE NEAR b.)       *'
      WRITE (*,FMT=9110) STR
      WRITE (*,FMT=9110) BLNK
      STR(1) = ' *   (3) LIMIT CIRCLE, NON-OSCILLATORY '
      STR(2) = 'AT b ?                                *'
      WRITE (*,FMT=9110) STR
      WRITE (*,FMT=9110) BLNK
      STR(1) = ' *   (4) LIMIT CIRCLE, OSCILLATORY AT b'
      STR(2) = ' ?                                    *'
      WRITE (*,FMT=9110) STR
      WRITE (*,FMT=9110) BLNK
      STR(1) = ' *   (5) LIMIT POINT AT b ?            '
      STR(2) = '                                      *'
      WRITE (*,FMT=9110) STR
      WRITE (*,FMT=9110) BLNK
      STR(1) = ' *   (6) NOT SPECIFIED (BUT NOT LIMIT C'
      STR(2) = 'IRCLE OSCILLATORY) WITH DEFAULT       *'
      WRITE (*,FMT=9110) STR
      STR(1) = ' *       BOUNDARY CONDITION AT b ?     '
      STR(2) = '                                      *'
      WRITE (*,FMT=9110) STR
      WRITE (*,FMT=9110) BLNK
      STR(1) = ' * ENTER THE NUMBER OF YOUR CHOICE: (h?'
      STR(2) = ')                                     *'
      WRITE (*,FMT=9110) STR
      WRITE (*,FMT=9110) STAR
      WRITE (*,FMT=*)
C
C     SPECIFY TYPE OF BOUNDARY CONDITION AT b.
C
      READ (*,FMT=9010) CHANS
      READ (CHANS,FMT=9020) HQ
      IF (HQ.EQ.'q' .OR. HQ.EQ.'Q') STOP
      IF (HQ.EQ.'h' .OR. HQ.EQ.'H') THEN
          CALL HELP(4)
          GO TO 120
      END IF
      YEH = HQ .EQ. '1' .OR. HQ .EQ. '2' .OR. HQ .EQ. '3' .OR.
     +      HQ .EQ. '4' .OR. HQ .EQ. '5' .OR. HQ .EQ. '6'
      IF (.NOT.YEH) GO TO 120
      READ (CHANS,FMT='(I32)') NANS
C
C     SET CHARACTER STRING CHB ACCORDING TO BOUNDARY CONDITION AT b.
C
      WRITE (*,FMT=*)
      IF (NANS.EQ.1) THEN
          REGB = .TRUE.
          CHB = CH1
          NCB = 1
          P0ATB = -1.0D0
          QFATB = 1.0D0
          IF (RITE) WRITE (T22,FMT=*) ' Endpoint b is Regular. '
      ELSE IF (NANS.EQ.2) THEN
          CHB = CH2
          NCB = 2
          IF (RITE) WRITE (T22,FMT=*) ' Endpoint b is Weakly Regular. '
      ELSE IF (NANS.EQ.3) THEN
          LCB = .true.
          CHB = CH3
          NCB = 3
          IF (RITE) WRITE (T22,FMT=*) ' Endpoint b is Limit Circle,',
     +        ' Non-Oscillatory. '
      ELSE IF (NANS.EQ.4) THEN
          LCB = .true.
          CHB = CH4
          NCB = 4
          IF (RITE) WRITE (T22,FMT=*) ' Endpoint b is Limit Circle,',
     +        ' Oscillatory. '
      ELSE IF (NANS.EQ.5) THEN
          CHB = CH5
          NCB = 5
          IF (RITE) WRITE (T22,FMT=*) ' Endpoint b is Limit Point. '
      ELSE
          CHB = CH6
          NCB = 6
          IF (RITE) WRITE (T22,FMT=*)
     +        ' Endpoint b is Singular, unspecified. '
      END IF
      RESP = NANS
      WRITE (*,FMT=*)
      WRITE (*,FMT=*)
      WRITE (*,FMT=*)
      WRITE (*,FMT=*)
      RETURN
C
   11 CONTINUE
  130 CONTINUE
      WRITE (*,FMT=*)
      WRITE (*,FMT=*) ' ********************************************* '
      WRITE (*,FMT=*) ' * IS p = 0. AT b ? (Y/N) (h?)               * '
      WRITE (*,FMT=*) ' ********************************************* '
      WRITE (*,FMT=*)
      READ (*,FMT=9010) CHANS
      READ (CHANS,FMT=9020) HQ
      IF (HQ.EQ.'q' .OR. HQ.EQ.'Q') STOP
      IF (HQ.EQ.'h' .OR. HQ.EQ.'H') THEN
          CALL HELP(11)
          GO TO 130
      END IF
      READ (CHANS,FMT=9020) YN
      YEH = YN .EQ. 'y' .OR. YN .EQ. 'Y'
      IF (.NOT. (YEH.OR.YN.EQ.'n'.OR.YN.EQ.'N')) GO TO 130
      P0ATB = -1.0D0
      IF (YEH) P0ATB = 1.0D0
      IF (RITE) THEN
          IF (YEH) THEN
              WRITE (T22,FMT=*) ' p is zero at b. '
          ELSE
              WRITE (T22,FMT=*) ' p is not zero at b. '
          END IF
      END IF

  140 CONTINUE
      WRITE (*,FMT=*)
      WRITE (*,FMT=*) ' ********************************************* '
      WRITE (*,FMT=*) ' * IS EITHER OF THE COEFFICIENT FUNCTIONS    * '
      WRITE (*,FMT=*) ' *                                           * '
      WRITE (*,FMT=*) ' *  q OR w UNBOUNDED NEAR b ? (Y/N) (h?)     * '
      WRITE (*,FMT=*) ' ********************************************* '
      WRITE (*,FMT=*)
      READ (*,FMT=9010) CHANS
      READ (CHANS,FMT=9020) HQ
      IF (HQ.EQ.'q' .OR. HQ.EQ.'Q') STOP
      IF (HQ.EQ.'h' .OR. HQ.EQ.'H') THEN
          CALL HELP(11)
          GO TO 140
      END IF
      YEH = HQ .EQ. 'y' .OR. HQ .EQ. 'Y' .OR. HQ .EQ. 'n' .OR.
     +      HQ .EQ. 'N'
      IF (.NOT.YEH) GO TO 140
      YEH = HQ .EQ. 'y' .OR. HQ .EQ. 'Y'
      QFATB = 1.0D0
      IF (YEH) QFATB = -1.0D0
      IF (RITE) THEN
          IF (YEH) THEN
              WRITE (T22,FMT=*) ' either q or w is unbounded near b. '
          ELSE
              WRITE (T22,FMT=*) ' both q and w are bounded near b. '
          END IF
      END IF
      IF (.NOT.YEH) RESP = 0
      RETURN
C
   12 CONTINUE
  150 CONTINUE
      WRITE (*,FMT=*)
      WRITE (*,FMT=*)
      WRITE (*,FMT=*)
      WRITE (*,FMT=*)
      WRITE (*,FMT=*)
      WRITE (*,FMT=*)
     +  ' ************************************************ '
      WRITE (*,FMT=*)
     +  ' * THIS PROBLEM IS ON THE INTERVAL              * '
      WRITE (*,FMT=*)
     +  ' *                                              * '
      IF (INTAB.EQ.1) THEN
          WRITE (*,FMT=9070) A,B
          WRITE (*,FMT=*) ' *                            ',FILLB
          WRITE (*,FMT=*) ' * ENDPOINT a IS  ',CHA
          WRITE (*,FMT=*) ' *                            ',FILLB
          IF (P0ATA.GT.0.0D0) THEN
              WRITE (*,FMT=*) ' * p IS ZERO AT a             ',FILLB
              WRITE (*,FMT=*) ' *                            ',FILLB
          END IF
          IF (QFATA.LT.0.0D0) THEN
              WRITE (*,FMT=*) ' * q OR w IS NOT BOUNDED AT a ',FILLB
              WRITE (*,FMT=*) ' *                            ',FILLB
          END IF
          WRITE (*,FMT=*) ' * ENDPOINT b IS  ',CHB
          WRITE (*,FMT=*) ' *                            ',FILLB
          IF (P0ATB.GT.0.0D0) THEN
              WRITE (*,FMT=*) ' * p IS ZERO AT b             ',FILLB
              WRITE (*,FMT=*) ' *                            ',FILLB
          END IF
          IF (QFATB.LT.0.0D0) THEN
              WRITE (*,FMT=*) ' * q OR w IS NOT BOUNDED AT b ',FILLB
              WRITE (*,FMT=*) ' *                            ',FILLB
          END IF
      ELSE IF (INTAB.EQ.2) THEN
          WRITE (*,FMT=9080) A,INFP
          WRITE (*,FMT=*) ' *                            ',FILLB
          WRITE (*,FMT=*) ' * ENDPOINT a IS  ',CHA
          WRITE (*,FMT=*) ' *                            ',FILLB
          IF (P0ATA.GT.0.0D0) THEN
              WRITE (*,FMT=*) ' * p IS ZERO AT a             ',FILLB
              WRITE (*,FMT=*) ' *                            ',FILLB
          END IF
          IF (QFATA.LT.0.0D0) THEN
              WRITE (*,FMT=*) ' * q OR w IS NOT BOUNDED AT a ',FILLB
              WRITE (*,FMT=*) ' *                            ',FILLB
          END IF
          WRITE (*,FMT=*) ' * ENDPT +INF IS  ',CHB
          WRITE (*,FMT=*) ' *                            ',FILLB
      ELSE IF (INTAB.EQ.3) THEN
          WRITE (*,FMT=9090) INFM,B
          WRITE (*,FMT=*) ' *                            ',FILLB
          WRITE (*,FMT=*) ' * ENDPT -INF IS  ',CHA
          WRITE (*,FMT=*) ' *                            ',FILLB
          WRITE (*,FMT=*) ' * ENDPOINT b IS  ',CHB
          WRITE (*,FMT=*) ' *                            ',FILLB
          IF (P0ATB.GT.0.0D0) THEN
              WRITE (*,FMT=*) ' * p IS ZERO AT b             ',FILLB
              WRITE (*,FMT=*) ' *                            ',FILLB
          END IF
          IF (QFATB.LT.0.0D0) THEN
              WRITE (*,FMT=*) ' * q OR w IS NOT BOUNDED AT b ',FILLB
              WRITE (*,FMT=*) ' *                            ',FILLB
          END IF
      ELSE
          WRITE (*,FMT=9100) INFM,INFP
          WRITE (*,FMT=*) ' *                            ',FILLB
          WRITE (*,FMT=*) ' * ENDPT -INF IS  ',CHA
          WRITE (*,FMT=*) ' *                            ',FILLB
          WRITE (*,FMT=*) ' * ENDPT +INF IS  ',CHB
          WRITE (*,FMT=*) ' *                            ',FILLB
      END IF
      WRITE (*,FMT=*)
     +  ' * IS THIS THE PROBLEM YOU WANT ? (Y/N) (h?)    * '
      WRITE (*,FMT=*)
     +  ' ************************************************ '
      WRITE (*,FMT=*)
      READ (*,FMT=9010) CHANS
      READ (CHANS,FMT=9020) HQ
      IF (HQ.EQ.'q' .OR. HQ.EQ.'Q') STOP
      IF (HQ.EQ.'h' .OR. HQ.EQ.'H') THEN
          CALL HELP(3)
          GO TO 150
      END IF
      YEH = HQ .EQ. 'y' .OR. HQ .EQ. 'Y' .OR. HQ .EQ. 'n' .OR.
     +      HQ .EQ. 'N'
      IF (.NOT.YEH) GO TO 150
      IF (HQ.NE.'y' .AND. HQ.NE.'Y') RESP = 0
      RETURN
C
   13 CONTINUE
  160 CONTINUE
      WRITE (*,FMT=*)
      WRITE (*,FMT=*) ' ********************************************* '
      WRITE (*,FMT=*) ' * DO YOU WANT TO RE-DO                      * '
      WRITE (*,FMT=*) ' *                                           * '
      WRITE (*,FMT=*) ' *    (1)  ENDPOINT a                        * '
      WRITE (*,FMT=*) ' *                                           * '
      WRITE (*,FMT=*) ' *    (2)  ENDPOINT b                        * '
      WRITE (*,FMT=*) ' *                                           * '
      WRITE (*,FMT=*) ' *    (3)  BOTH ENDPOINTS a AND b ?          * '
      WRITE (*,FMT=*) ' *                                           * '
      WRITE (*,FMT=*) ' * ENTER THE NUMBER OF YOUR CHOICE: (h?)     * '
      WRITE (*,FMT=*) ' ********************************************* '
      WRITE (*,FMT=*)
      READ (*,FMT=9010) CHANS
      READ (CHANS,FMT=9020) HQ
      IF (HQ.EQ.'q' .OR. HQ.EQ.'Q') STOP
      IF (HQ.EQ.'h' .OR. HQ.EQ.'H') THEN
          CALL HELP(3)
          GO TO 160
      END IF
      YEH = HQ .EQ. '1' .OR. HQ .EQ. '2' .OR. HQ .EQ. '3'
      IF (.NOT.YEH) GO TO 160
      READ (CHANS,FMT='(I32)') NANS
      RESP = NANS
      WRITE (*,FMT=*)
      WRITE (*,FMT=*)
      WRITE (*,FMT=*)
      WRITE (*,FMT=*)
      WRITE (*,FMT=*)
      RETURN
C
   14 CONTINUE
  170 CONTINUE
      IF (RITE) WRITE (T22,FMT=*)
     +    '----------------------------------------'
      IF (RITE) WRITE (T22,FMT=*)
      WRITE (*,FMT=*)
      WRITE (*,FMT=*)
     +  ' *************************************************'
      WRITE (*,FMT=*)
     +  ' * DO YOU WANT TO COMPUTE                        *'
      WRITE (*,FMT=*)
     +  ' *                                               *'
      WRITE (*,FMT=*)
     +  ' *   (1) AN EIGENVALUE, OR SERIES OF EIGENVALUES *'
      WRITE (*,FMT=*)
     +  ' *                                               *'
      WRITE (*,FMT=*)
     +  ' *   (2) SOLUTION TO AN INITIAL VALUE PROBLEM ?  *'
      WRITE (*,FMT=*)
     +  ' *                                               *'
      WRITE (*,FMT=*)
     +  ' * ENTER THE NUMBER OF YOUR CHOICE: (h?)         *'
      WRITE (*,FMT=*)
     +  ' *************************************************'
      WRITE (*,FMT=*)
      READ (*,FMT=9010) CHANS
      READ (CHANS,FMT=9020) HQ
      IF (HQ.EQ.'q' .OR. HQ.EQ.'Q') STOP
      IF (HQ.EQ.'h' .OR. HQ.EQ.'H') THEN
          CALL HELP(12)
          GO TO 170
      END IF
      YEH = HQ .EQ. '1' .OR. HQ .EQ. '2'
      IF (.NOT.YEH) GO TO 170
      READ (CHANS,FMT='(I32)') NANS
      RESP = NANS
      RETURN
C
   15 CONTINUE
  200 CONTINUE
      WRITE (*,FMT=*) ' ******************************************** '
      WRITE (*,FMT=*) ' * IS THE BOUNDARY CONDITION PERIODIC ?     * '
      WRITE (*,FMT=*) ' *   OR COUPLED ?                           * '
      WRITE (*,FMT=*) ' *                                          * '
      WRITE (*,FMT=*) ' *   (I.E., y(b) = c*y(a)                   * '
      WRITE (*,FMT=*) ' *    &   p(b)*y''(b) = (1/c)*p(a)*y''(a)     *'
      WRITE (*,FMT=*) ' *                                          * '
      WRITE (*,FMT=*) ' *   OR SOME OTHER COUPLED CONDITION)       * '
      WRITE (*,FMT=*) ' *                                          * '
      WRITE (*,FMT=*) ' * ANSWER (Y/N): (h?)                       * '
      WRITE (*,FMT=*) ' ******************************************** '
      WRITE (*,FMT=*)
      READ (*,FMT=9010) CHANS
      READ (CHANS,FMT=9020) HQ
      IF (HQ.EQ.'q' .OR. HQ.EQ.'Q') STOP
      IF (HQ.EQ.'h' .OR. HQ.EQ.'H') THEN
          CALL HELP(7)
          GO TO 200
      END IF
      YEH = HQ .EQ. 'y' .OR. HQ .EQ. 'Y' .OR. HQ .EQ. 'n' .OR.
     +      HQ .EQ. 'N'
      IF (.NOT.YEH) GO TO 200
      YEH = HQ .EQ. 'y' .OR. HQ .EQ. 'Y'
      IF (.NOT.YEH) RESP = 0
      WRITE (*,FMT=*)
      WRITE (*,FMT=*)
      WRITE (*,FMT=*)
      WRITE (*,FMT=*)
      WRITE (*,FMT=*)
      RETURN
C
   16 CONTINUE
      IF (NCA.LE.2 .OR. CHA.EQ.CH6) THEN
  190     CONTINUE
          WRITE (*,FMT=*) ' ****************************************** '
          WRITE (*,FMT=*) ' * IS THE BOUNDARY CONDITION AT a         * '
          WRITE (*,FMT=*) ' *                                        * '
          WRITE (*,FMT=*) ' *    (1) THE DIRICHLET CONDITION         * '
          WRITE (*,FMT=*) ' *        (I.E., y(a) = 0.0)              * '
          WRITE (*,FMT=*) ' *                                        * '
          WRITE (*,FMT=*) ' *    (2) THE NEUMANN CONDITION           * '
          WRITE (*,FMT=*) ' *        (I.E., y''(a) = 0.0)             *'
          WRITE (*,FMT=*) ' *                                        * '
          WRITE (*,FMT=*) ' *    (3) A MORE GENERAL LINEAR           * '
          WRITE (*,FMT=*) ' *        BOUNDARY CONDITION              * '
          WRITE (*,FMT=*) ' *        A1*y(a) + A2*(py'')(a) = 0 ?   *'
          WRITE (*,FMT=*) ' *                                        * '
          WRITE (*,FMT=*) ' * ENTER THE NUMBER OF YOUR CHOICE: (h?)  * '
          WRITE (*,FMT=*) ' ****************************************** '
          WRITE (*,FMT=*)
          READ (*,FMT=9010) CHANS
          READ (CHANS,FMT=9020) HQ
          IF (HQ.EQ.'q' .OR. HQ.EQ.'Q') STOP
          IF (HQ.EQ.'h' .OR. HQ.EQ.'H') THEN
              CALL HELP(7)
              GO TO 190
          END IF
          YEH = HQ .EQ. '1' .OR. HQ .EQ. '2' .OR. HQ .EQ. '3'
          IF (.NOT.YEH) GO TO 190
          READ (CHANS,FMT='(I32)') NANS
          IF (NANS.EQ.1) THEN
              A1 = 1.0D0
              A2 = 0.0D0
              IF (RITE) WRITE (T22,FMT=*) ' Dirichlet B.C. at a. '
          ELSE IF (NANS.EQ.2) THEN
              A1 = 0.0D0
              A2 = 1.0D0
              IF (RITE) WRITE (T22,FMT=*) ' Neumann B.C. at a. '
          ELSE
  207         CONTINUE
              WRITE (*,FMT=*)
     +          ' *************************************** '
              WRITE (*,FMT=*)
     +          ' * CHOOSE A1,A2: (h?)                  * '
              WRITE (*,FMT=*)
     +          ' *************************************** '
              WRITE (*,FMT=*)
              WRITE (*,FMT=*) ' A1,A2 = '
              READ (*,FMT=9010) CHANS
              READ (CHANS,FMT=9020) HQ
              IF (HQ.EQ.'q' .OR. HQ.EQ.'Q') STOP
              IF (HQ.EQ.'h' .OR. HQ.EQ.'H') THEN
                  CALL HELP(7)
                  GO TO 207
              END IF
              CALL LSTDIR(CHANS,I,ICOL)
              I1 = ICOL(1) - 1
              FMAT = FMT2(I1)
              READ (CHANS,FMT=FMAT) A1,A2
              IF (RITE) WRITE (T22,FMT=*) ' A1,A2 = ',A1,A2
          END IF
      ELSE IF (LCA) THEN
  210     CONTINUE
          IF (RITE) THEN
              WRITE (T22,FMT=*) ' The B.C. at a is '
              WRITE (T22,FMT=*) ' A1*[y,u](a) + A2*[y,v](a) = 0. '
          END IF
          WRITE (*,FMT=*) ' ****************************************** '
          WRITE (*,FMT=*) ' * THE BOUNDARY CONDITION AT a IS         * '
          WRITE (*,FMT=*) ' *                                        * '
          WRITE (*,FMT=*) ' *    A1*[y,u](a) + A2*[y,v](a) = 0,      * '
          WRITE (*,FMT=*) ' *                                        * '
          WRITE (*,FMT=*) ' *  WHERE THE CONSTANTS A1 AND A2         * '
          WRITE (*,FMT=*) ' *    MAY BE CHOSEN ARBITRARILY.          * '
          WRITE (*,FMT=*) ' *                                        * '
          WRITE (*,FMT=*) ' * CHOOSE A1,A2: (h?)                     * '
          WRITE (*,FMT=*) ' ****************************************** '
          WRITE (*,FMT=*)
          WRITE (*,FMT=*) ' A1,A2 = '
          READ (*,FMT=9010) CHANS
          READ (CHANS,FMT=9020) HQ
          IF (HQ.EQ.'q' .OR. HQ.EQ.'Q') STOP
          IF (HQ.EQ.'h' .OR. HQ.EQ.'H') THEN
              CALL HELP(7)
              GO TO 210
          END IF
          CALL LSTDIR(CHANS,I,ICOL)
          I1 = ICOL(1) - 1
          FMAT = FMT2(I1)
          READ (CHANS,FMT=FMAT) A1,A2
          IF (RITE) WRITE (T22,FMT=*) ' A1,A2 = ',A1,A2
      END IF
      RETURN
C
   17 CONTINUE
      IF (NCB.LE.2 .OR. CHB.EQ.CH6) THEN
  220     CONTINUE
          WRITE (*,FMT=*) ' ****************************************** '
          WRITE (*,FMT=*) ' * IS THE BOUNDARY CONDITION AT b         * '
          WRITE (*,FMT=*) ' *                                        * '
          WRITE (*,FMT=*) ' *    (1) THE DIRICHLET CONDITION         * '
          WRITE (*,FMT=*) ' *        (I.E., y(b) = 0.0)              * '
          WRITE (*,FMT=*) ' *                                        * '
          WRITE (*,FMT=*) ' *    (2) THE NEUMANN CONDITION           * '
          WRITE (*,FMT=*) ' *        (I.E., y''(b) = 0.0)             *'
          WRITE (*,FMT=*) ' *                                        * '
          WRITE (*,FMT=*) ' *    (3) A MORE GENERAL LINEAR           * '
          WRITE (*,FMT=*) ' *        BOUNDARY CONDITION              * '
          WRITE (*,FMT=*) ' *        B1*y(b) + B2*(py'')(b) = 0 ?   *'
          WRITE (*,FMT=*) ' *                                        * '
          WRITE (*,FMT=*) ' * ENTER THE NUMBER OF YOUR CHOICE: (h?)  * '
          WRITE (*,FMT=*) ' ****************************************** '
          WRITE (*,FMT=*)
          READ (*,FMT=9010) CHANS
          READ (CHANS,FMT=9020) HQ
          IF (HQ.EQ.'q' .OR. HQ.EQ.'Q') STOP
          IF (HQ.EQ.'h' .OR. HQ.EQ.'H') THEN
              CALL HELP(7)
              GO TO 220
          END IF
          YEH = HQ .EQ. '1' .OR. HQ .EQ. '2' .OR. HQ .EQ. '3'
          IF (.NOT.YEH) GO TO 220
          READ (CHANS,FMT='(I32)') NANS
          IF (NANS.EQ.1) THEN
              B1 = 1.0D0
              B2 = 0.0D0
              IF (RITE) WRITE (T22,FMT=*) ' Dirichlet B.C. at b. '
          ELSE IF (NANS.EQ.2) THEN
              B1 = 0.0D0
              B2 = 1.0D0
              IF (RITE) WRITE (T22,FMT=*) ' Neumann B.C. at b. '
          ELSE
  230         CONTINUE
              WRITE (*,FMT=*)
     +          ' *************************************** '
              WRITE (*,FMT=*)
     +          ' * CHOOSE B1,B2: (h?)                  * '
              WRITE (*,FMT=*)
     +          ' *************************************** '
              WRITE (*,FMT=*)
              WRITE (*,FMT=*) ' B1,B2 = '
              READ (*,FMT=9010) CHANS
              READ (CHANS,FMT=9020) HQ
              IF (HQ.EQ.'q' .OR. HQ.EQ.'Q') STOP
              IF (HQ.EQ.'h' .OR. HQ.EQ.'H') THEN
                  CALL HELP(7)
                  GO TO 230
              END IF
              CALL LSTDIR(CHANS,I,ICOL)
              I1 = ICOL(1) - 1
              FMAT = FMT2(I1)
              READ (CHANS,FMT=FMAT) B1,B2
              IF (RITE) WRITE (T22,FMT=*) ' B1,B2 = ',B1,B2
          END IF
      ELSE IF (LCB) THEN
  240     CONTINUE
          IF (RITE) THEN
              WRITE (T22,FMT=*) ' The B.C. at b is '
              WRITE (T22,FMT=*) ' B1*[y,u](b) + B2*[y,v](b) = 0. '
          END IF
          WRITE (*,FMT=*) ' ****************************************** '
          WRITE (*,FMT=*) ' * THE BOUNDARY CONDITION AT b IS         * '
          WRITE (*,FMT=*) ' *                                        * '
          WRITE (*,FMT=*) ' *    B1*[y,u](b) + B2*[y,v](b) = 0,      * '
          WRITE (*,FMT=*) ' *                                        * '
          WRITE (*,FMT=*) ' *  WHERE THE CONSTANTS B1 AND B2         * '
          WRITE (*,FMT=*) ' *    MAY BE CHOSEN ARBITRARILY.          * '
          WRITE (*,FMT=*) ' *                                        * '
          WRITE (*,FMT=*) ' * CHOOSE B1,B2: (h?)                     * '
          WRITE (*,FMT=*) ' ****************************************** '
          WRITE (*,FMT=*)
          WRITE (*,FMT=*) ' B1,B2 = '
          READ (*,FMT=9010) CHANS
          READ (CHANS,FMT=9020) HQ
          IF (HQ.EQ.'q' .OR. HQ.EQ.'Q') STOP
          IF (HQ.EQ.'h' .OR. HQ.EQ.'H') THEN
              CALL HELP(7)
              GO TO 240
          END IF
          CALL LSTDIR(CHANS,I,ICOL)
          I1 = ICOL(1) - 1
          FMAT = FMT2(I1)
          READ (CHANS,FMT=FMAT) B1,B2
          IF (RITE) WRITE (T22,FMT=*) ' B1,B2 = ',B1,B2
      END IF
      RETURN
C
   18 CONTINUE
  310 CONTINUE
      WRITE (*,FMT=*) ' ******************************************** '
      WRITE (*,FMT=*) ' * IS THIS PROBLEM:                         * '
      WRITE (*,FMT=*) ' *                                          * '
      WRITE (*,FMT=*) ' *   (1) PERIODIC ?                         * '
      WRITE (*,FMT=*) ' *         (I.E., y(b) = y(a)               * '
      WRITE (*,FMT=*) ' *          &   p(b)*y''(b) = p(a)*y''(a) )   *'
      WRITE (*,FMT=*) ' *                                          * '
      WRITE (*,FMT=*) ' *   (2) SEMI-PERIODIC ?                    * '
      WRITE (*,FMT=*) ' *         (I.E., y(b) = -y(a)              * '
      WRITE (*,FMT=*) ' *          &   p(b)*y''(b) = -p(a)*y''(a) )  *'
      WRITE (*,FMT=*) ' *                                          * '
      WRITE (*,FMT=*) ' *   (3) GENERAL PERIODIC TYPE ?            * '
      WRITE (*,FMT=*) ' *         (I.E., y(b) = c*y(a)             * '
      WRITE (*,FMT=*) ' *          &   p(b)*y''(b) = p(a)*y''(a)/c   *'
      WRITE (*,FMT=*) ' *            for some number c .NE. 0. )   * '
      WRITE (*,FMT=*) ' *                                          * '
      WRITE (*,FMT=*) ' *   (4) MORE GENERAL COUPLED TYPE ?        * '
      WRITE (*,FMT=*) ' *                                          * '
      WRITE (*,FMT=*) ' * ENTER THE NUMBER OF YOUR CHOICE: (h)?    * '
      WRITE (*,FMT=*) ' ******************************************** '
      WRITE (*,FMT=*)
      READ (*,FMT=9010) CHANS
      READ (CHANS,FMT=9020) HQ
      IF (HQ.EQ.'q' .OR. HQ.EQ.'Q') STOP
      IF (HQ.EQ.'h' .OR. HQ.EQ.'H') THEN
          CALL HELP(7)
          GO TO 310
      END IF
      YEH = HQ .EQ. '1' .OR. HQ .EQ. '2' .OR. HQ .EQ. '3' .OR.
     +      HQ .EQ. '4'
      IF (.NOT.YEH) GO TO 310
      READ (CHANS,FMT='(I32)') NANS
      IF (NANS.EQ.1) THEN
          CC = 1.0D0
          ALF = 0.0D0
          K11 = 1.0D0
          K12 = 0.0D0
          K21 = 0.0D0
          K22 = 1.0D0
          IF (RITE) WRITE (T22,FMT=*) ' The B.C. is Periodic. '
      ELSE IF (NANS.EQ.2) THEN
          CC = -1.0D0
          ALF = 0.0D0
          K11 = -1.0D0
          K12 = 0.0D0
          K21 = 0.0D0
          K22 = -1.0D0
          IF (RITE) WRITE (T22,FMT=*) ' The B.C. is Semi-Periodic. '
      ELSE IF (NANS.EQ.3) THEN
  320     CONTINUE
          WRITE (*,FMT=*) ' ****************************************** '
          WRITE (*,FMT=*) ' * INPUT c: (h?)                          * '
          WRITE (*,FMT=*) ' ****************************************** '
          WRITE (*,FMT=*) ' c = '
          READ (*,FMT=9010) CHANS
          READ (CHANS,FMT=9020) HQ
          IF (HQ.EQ.'q' .OR. HQ.EQ.'Q') STOP
          IF (HQ.EQ.'h' .OR. HQ.EQ.'H') THEN
              CALL HELP(7)
              GO TO 320
          END IF
          READ (CHANS,FMT='(F32.0)') CC
          ALF = 0.0D0
          K11 = CC
          K12 = 0.0D0
          K21 = 0.0D0
          K22 = 1.0D0/CC
          IF (RITE) WRITE (T22,FMT=*)
     +        ' The B.C. is General Periodic type. '
          IF (RITE) WRITE (T22,FMT=*) ' Parameter c = ',CC
      ELSE
  322     CONTINUE
          WRITE (*,FMT=*)
     +      ' ************************************************** '
          WRITE (*,FMT=*)
     +      ' * FOR THIS PROBLEM, THE GENERAL COUPLED          * '
          WRITE (*,FMT=*)
     +      ' *   BOUNDARY CONDITIONS ARE:                     * '
          WRITE (*,FMT=*)
     +      ' *                                                * '
          IF (NCA.LE.2 .AND. NCB.LE.2) THEN
              WRITE (*,FMT=*)
     +          ' * y( b)  =  EX*{k11*y(a) + k12*py''(a)}           *'
              WRITE (*,FMT=*)
     +          ' * py''(b) = EX*{k21*y(a) + k22*py''(a)}            *'
          ELSE IF (NCA.GE.3 .AND. NCB.LE.2) THEN
              WRITE (*,FMT=*)
     +          ' * y(b) = EX*{k11*n*[y,u](a)+k12[y,v](a)}/d       *'
              WRITE (*,FMT=*)
     +          ' * py''(b) = EX*{k21*n*[y,u](a) + k22*[y,v](a)}/d  *'
              WRITE (*,FMT=*)
     +          ' *                                                * '
              WRITE (*,FMT=*)
     +          ' * WHERE d = sqrt(abs([u,v](a)))                  * '
              WRITE (*,FMT=*)
     +          ' *       n = +1 or -1                             * '
          ELSE IF (NCA.LE.2 .AND. NCB.GE.3) THEN
              WRITE (*,FMT=*)
     +          ' * [y,U](b)*N/D = EX*{k11*y(a) + k12*py''(a)}      *'
              WRITE (*,FMT=*)
     +          ' * [y,V](b)/D = EX*{k21*y(a) + k22*py''(a)}        *'
              WRITE (*,FMT=*)
     +          ' *                                                * '
              WRITE (*,FMT=*)
     +          ' * WHERE D = sqrt(abs([U,V](b)))                  * '
              WRITE (*,FMT=*)
     +          ' *       N = +1 OR -1                             * '
          ELSE
              WRITE (*,FMT=*)
     +          ' * [y,U](b)*N/D=EX*{k11*n*[y,u](a)+k12*[y,v](a)}/d*'
              WRITE (*,FMT=*)
     +          ' * [y,V](b)/D=EX*{k21*n*[y,u](a)+k22*[y,v](a)}/d  *'
              WRITE (*,FMT=*)
     +          ' *                                                * '
              WRITE (*,FMT=*)
     +          ' *       D = sqrt(abs([U,V](b)))                  * '
              WRITE (*,FMT=*)
     +          ' *       N = +1 OR -1                             * '
              WRITE (*,FMT=*)
     +          ' *       D = SQRT(ABS([U,V](A)))                  * '
              WRITE (*,FMT=*)
     +          ' *       N = +1 OR -1                             * '
          END IF
          WRITE (*,FMT=*)
     +      ' *                                                * '
          WRITE (*,FMT=*)
     +      ' * WHERE EX = EXP(I*ALFA)                         * '
          WRITE (*,FMT=*)
     +      ' * AND WHERE     ALFA                             * '
          WRITE (*,FMT=*)
     +      ' * AND THE       K11,K12                          * '
          WRITE (*,FMT=*)
     +      ' *               K21,K22                          * '
          WRITE (*,FMT=*)
     +      ' *                                                * '
          WRITE (*,FMT=*)
     +      ' * NEED TO BE CHOSEN. IT IS NECESSARY THAT        * '
          WRITE (*,FMT=*)
     +      ' *                                                * '
          WRITE (*,FMT=*)
     +      ' *      0.0 .LE. ALFA .LT. PI                     * '
          WRITE (*,FMT=*)
     +      ' *                                                * '
          WRITE (*,FMT=*)
     +      ' * INPUT ALFA : (H?)                              * '
          WRITE (*,FMT=*)
     +      ' ************************************************** '
          WRITE (*,FMT=*) ' ALFA = '
          READ (*,FMT=9010) CHANS
          READ (CHANS,FMT=9020) HQ
          IF (HQ.EQ.'q' .OR. HQ.EQ.'Q') STOP
          IF (HQ.EQ.'h' .OR. HQ.EQ.'H') THEN
              CALL HELP(7)
              GO TO 322
          END IF
          READ (CHANS,FMT='(1F32.0)') ALF
          WRITE (*,FMT=*)
     +      ' ********************************************* '
          WRITE (*,FMT=*)
     +      ' * FOR SELF ADJOINTNESS IT IS NECESSARY      * '
          WRITE (*,FMT=*)
     +      ' * THAT                                      * '
          WRITE (*,FMT=*)
     +      ' *                                           * '
          WRITE (*,FMT=*)
     +      ' *      K11*K22 - K12*K21 =  1               * '
          WRITE (*,FMT=*)
     +      ' *                                           * '
          WRITE (*,FMT=*)
     +      ' * INPUT K11,K12 :                           * '
          WRITE (*,FMT=*)
     +      ' ********************************************* '
          WRITE (*,FMT=*) ' K11,K12 = '
          READ (*,FMT=9010) CHANS
          READ (CHANS,FMT=9020) HQ
          IF (HQ.EQ.'q' .OR. HQ.EQ.'Q') STOP
          IF (HQ.EQ.'h' .OR. HQ.EQ.'H') THEN
              CALL HELP(7)
              GO TO 322
          END IF
          CALL LSTDIR(CHANS,I,ICOL)
          I1 = ICOL(1) - 1
          FMAT = FMT2(I1)
          READ (CHANS,FMT=FMAT) K11,K12
          WRITE (*,FMT=*)
     +      ' ********************************************* '
          WRITE (*,FMT=*)
     +      ' * INPUT K21,K22 :                           * '
          WRITE (*,FMT=*)
     +      ' ********************************************* '
          WRITE (*,FMT=*) ' K21,K22 = '
          READ (*,FMT=9010) CHANS
          READ (CHANS,FMT=9020) HQ
          IF (HQ.EQ.'q' .OR. HQ.EQ.'Q') STOP
          IF (HQ.EQ.'h' .OR. HQ.EQ.'H') THEN
              CALL HELP(7)
              GO TO 322
          END IF
          CALL LSTDIR(CHANS,I,ICOL)
          I1 = ICOL(1) - 1
          FMAT = FMT2(I1)
          READ (CHANS,FMT=FMAT) K21,K22
          WRITE (*,FMT=*)
          WRITE (*,FMT=*)
          WRITE (*,FMT=*)
          WRITE (*,FMT=*) ' ALFA = ',ALF
          WRITE (*,FMT=*) ' K11,K12 = ',K11,K12
          WRITE (*,FMT=*) ' K21,K22 = ',K21,K22
      END IF
      DETK = K11*K22 - K21*K12
      TMP = ABS(DETK-1.0D0)
      IF (TMP.GT.0.01D0) THEN
          WRITE (*,FMT=*) ' WARNING: K11*K22-K12*K21 IS NOT = 1 '
      END IF
      WRITE (*,FMT=*)
      WRITE (*,FMT=*)
      IF (RITE) THEN
          WRITE (T22,FMT=*)
          WRITE (T22,FMT=*) '   ALFA = ',ALF
          WRITE (T22,FMT=*)
          WRITE (T22,FMT=*) ' K11,K12 = ',K11,K12
          WRITE (T22,FMT=*) ' K21,K22 = ',K21,K22
          WRITE (T22,FMT=*) '   DET(K) = ',DETK
          IF (TMP.GT.0.01D0) THEN
              WRITE (T22,FMT=*) ' WARNING: K11*K22-K12*K21 IS NOT = 1 '
          END IF
      END IF
      IF (RITE) WRITE (T22,FMT=*) ' -------------------------------'//
     +    FILLC
      RETURN
C
   19 CONTINUE
  350 CONTINUE
      WRITE (*,FMT=*)
      WRITE (*,FMT=*)
      WRITE (*,FMT=*)
      WRITE (*,FMT=*)
      WRITE (*,FMT=*)
      WRITE (*,FMT=*)
     +  ' *********************************************** '
      WRITE (*,FMT=*)
     +  ' * DO YOU WANT TO COMPUTE THE SOLUTION TO:     * '
      WRITE (*,FMT=*)
     +  ' *                                             * '
      WRITE (*,FMT=*)
     +  ' *    (1) AN INITIAL VALUE PROBLEM FROM ONE    * '
      WRITE (*,FMT=*)
     +  ' *        END OF THE INTERVAL TO THE OTHER     * '
      WRITE (*,FMT=*)
     +  ' *                                             * '
      WRITE (*,FMT=*)
     +  ' *    (2) INITIAL VALUE PROBLEMS FROM BOTH     * '
      WRITE (*,FMT=*)
     +  ' *        ENDS TO A MIDPOINT ?                 * '
      WRITE (*,FMT=*)
     +  ' *                                             * '
      WRITE (*,FMT=*)
     +  ' * ENTER THE NUMBER OF YOUR CHOICE: (h?)       * '
      WRITE (*,FMT=*)
     +  ' *********************************************** '
      WRITE (*,FMT=*)
      READ (*,FMT=9010) CHANS
      READ (CHANS,FMT=9020) HQ
      IF (HQ.EQ.'q' .OR. HQ.EQ.'Q') STOP
      IF (HQ.EQ.'h' .OR. HQ.EQ.'H') THEN
          CALL HELP(12)
          GO TO 350
      END IF
      YEH = HQ .EQ. '1' .OR. HQ .EQ. '2'
      IF (.NOT.YEH) GO TO 350
      READ (CHANS,FMT='(I32)') NIVP
      IF (NIVP.EQ.1) THEN
          WRITE (*,FMT=*)
          WRITE (*,FMT=*)
          WRITE (*,FMT=*)
          WRITE (*,FMT=*)
  360     CONTINUE
          WRITE (*,FMT=*)
          WRITE (*,FMT=*)
     +      ' ********************************************* '
          WRITE (*,FMT=*)
     +      ' * WHICH IS THE INITIAL POINT: a OR b ? (h?) * '
          WRITE (*,FMT=*)
     +      ' ********************************************* '
          WRITE (*,FMT=*)
          WRITE (*,FMT=*) ' INITIAL POINT IS '
          READ (*,FMT=9010) CHANS
          READ (CHANS,FMT=9020) HQ
          IF (HQ.EQ.'q' .OR. HQ.EQ.'Q') STOP
          IF (HQ.EQ.'h' .OR. HQ.EQ.'H') THEN
              CALL HELP(12)
              GO TO 360
          END IF
          READ (CHANS,FMT=9020) ANSCH
          IF (ANSCH.EQ.'a' .OR. ANSCH.EQ.'A') THEN
              NEND = 1
              IF (RITE) WRITE (T22,FMT=*) ' The Initial Point for this',
     +            ' Initial Value Problem is a. '
              IF (NCA.LE.4 .OR. NCA.EQ.6) THEN
                  WRITE (*,FMT=*)
                  WRITE (*,FMT=*)
                  WRITE (*,FMT=*)
                  WRITE (*,FMT=*)
  370             CONTINUE
                  WRITE (*,FMT=*)
                  WRITE (*,FMT=*)
     +              ' ************************************ '
                  WRITE (*,FMT=*)
     +              ' * THE INITIAL CONDITIONS AT a ARE  * '
                  WRITE (*,FMT=*)
     +              ' *                                  * '
                  IF (NCA.LE.2 .OR. NCA.EQ.6) THEN
                      WRITE (*,FMT=*)
     +                  ' *    y(a)=alfa1, py''(a)=alfa2      *'
                  ELSE
                      WRITE (*,FMT=*)
     +                  ' *  [y,u](a)=alfa1, [y,v](a)=alfa2  * '
                  END IF
                  WRITE (*,FMT=*)
     +              ' *                                  * '
                  WRITE (*,FMT=*)
     +              ' * WHERE THE CONSTANTS alfa1, alfa2 * '
                  WRITE (*,FMT=*)
     +              ' *    MAY BE CHOSEN ARBITRARILY.    * '
                  WRITE (*,FMT=*)
     +              ' *                                  * '
                  WRITE (*,FMT=*)
     +              ' * CHOOSE alfa1,alfa2: (h?)         * '
                  WRITE (*,FMT=*)
     +              ' ************************************ '
                  WRITE (*,FMT=*)
                  WRITE (*,FMT=*) ' alfa1,alfa2 = '
                  READ (*,FMT=9010) CHANS
                  READ (CHANS,FMT=9020) HQ
                  IF (HQ.EQ.'q' .OR. HQ.EQ.'Q') STOP
                  IF (HQ.EQ.'h' .OR. HQ.EQ.'H') THEN
                      CALL HELP(12)
                      GO TO 370
                  END IF
                  CALL LSTDIR(CHANS,I,ICOL)
                  I1 = ICOL(1) - 1
                  FMAT = FMT2(I1)
                  READ (CHANS,FMT=FMAT) ALFA1,ALFA2
                  IF (RITE) WRITE (T22,FMT=*) ' alfa1,alfa2 = ',ALFA1,
     +                ALFA2
                  A1 = ALFA2
                  A2 = -ALFA1
              END IF
          ELSE IF (ANSCH.EQ.'b' .OR. ANSCH.EQ.'B') THEN
              NEND = 2
              IF (RITE) WRITE (T22,FMT=*) ' The Initial Point for this',
     +            ' Initial Value Problem is b. '
              IF (NCB.LE.4 .OR. NCB.EQ.6) THEN
                  WRITE (*,FMT=*)
                  WRITE (*,FMT=*)
                  WRITE (*,FMT=*)
                  WRITE (*,FMT=*)
  380             CONTINUE
                  WRITE (*,FMT=*)
                  WRITE (*,FMT=*)
     +              ' ************************************ '
                  WRITE (*,FMT=*)
     +              ' * THE INITIAL CONDITIONS AT b ARE  * '
                  WRITE (*,FMT=*)
     +              ' *                                  * '
                  IF (NCB.LE.2 .OR. NCB.EQ.6) THEN
                      WRITE (*,FMT=*)
     +                  ' *    y(b)=beta1, py''(b)=beta2      *'
                  ELSE
                      WRITE (*,FMT=*)
     +                  ' *  [y,u](b)=beta1, [y,v](b)=beta2  * '
                  END IF
                  WRITE (*,FMT=*)
     +              ' *                                  * '
                  WRITE (*,FMT=*)
     +              ' * WHERE THE CONSTANTS beta1, beta2 * '
                  WRITE (*,FMT=*)
     +              ' *    MAY BE CHOSEN ARBITRARILY.    * '
                  WRITE (*,FMT=*)
     +              ' *                                  * '
                  WRITE (*,FMT=*)
     +              ' * CHOOSE beta1,beta2: (h?)         * '
                  WRITE (*,FMT=*)
     +              ' ************************************ '
                  WRITE (*,FMT=*)
                  WRITE (*,FMT=*) ' beta1,beta2 = '
                  READ (*,FMT=9010) CHANS
                  READ (CHANS,FMT=9020) HQ
                  IF (HQ.EQ.'q' .OR. HQ.EQ.'Q') STOP
                  IF (HQ.EQ.'h' .OR. HQ.EQ.'H') THEN
                      CALL HELP(12)
                      GO TO 380
                  END IF
                  CALL LSTDIR(CHANS,I,ICOL)
                  I1 = ICOL(1) - 1
                  FMAT = FMT2(I1)
                  READ (CHANS,FMT=FMAT) BETA1,BETA2
                  IF (RITE) WRITE (T22,FMT=*) ' beta1,beta2 = ',BETA1,
     +                BETA2
                  B1 = BETA2
                  B2 = -BETA1
              END IF
          END IF
      ELSE IF (NIVP.EQ.2) THEN
          NEND = 3
          IF (NCA.LE.4 .OR. NCA.EQ.6) THEN
              WRITE (*,FMT=*)
              WRITE (*,FMT=*)
              WRITE (*,FMT=*)
              WRITE (*,FMT=*)
  390         CONTINUE
              WRITE (*,FMT=*)
              WRITE (*,FMT=*) ' ************************************ '
              WRITE (*,FMT=*) ' * THE INITIAL CONDITIONS AT a ARE  * '
              WRITE (*,FMT=*) ' *                                  * '
              IF (NCA.LE.2 .OR. NCA.EQ.6) THEN
                  WRITE (*,FMT=*)
     +              ' *    y(a)=alfa1, py''(a)=alfa2      *'
              ELSE
                  WRITE (*,FMT=*)
     +              ' *  [y,u](a)=alfa1, [y,v](a)=alfa2  * '
              END IF
              WRITE (*,FMT=*) ' *                                  * '
              WRITE (*,FMT=*) ' * WHERE THE CONSTANTS alfa1, alfa2 * '
              WRITE (*,FMT=*) ' *    MAY BE CHOSEN ARBITRARILY.    * '
              WRITE (*,FMT=*) ' *                                  * '
              WRITE (*,FMT=*) ' * CHOOSE alfa1,alfa2: (h?)         * '
              WRITE (*,FMT=*) ' ************************************ '
              WRITE (*,FMT=*)
              WRITE (*,FMT=*) ' alfa1,alfa2 = '
              READ (*,FMT=9010) CHANS
              READ (CHANS,FMT=9020) HQ
              IF (HQ.EQ.'q' .OR. HQ.EQ.'Q') STOP
              IF (HQ.EQ.'h' .OR. HQ.EQ.'H') THEN
                  CALL HELP(12)
                  GO TO 390
              END IF
              CALL LSTDIR(CHANS,I,ICOL)
              I1 = ICOL(1) - 1
              FMAT = FMT2(I1)
              READ (CHANS,FMT=FMAT) ALFA1,ALFA2
              A1 = ALFA2
              A2 = -ALFA1
          END IF
          IF (NCB.LE.4 .OR. NCB.EQ.6) THEN
              WRITE (*,FMT=*)
              WRITE (*,FMT=*)
              WRITE (*,FMT=*)
              WRITE (*,FMT=*)
  400         CONTINUE
              WRITE (*,FMT=*)
              WRITE (*,FMT=*) ' ************************************ '
              WRITE (*,FMT=*) ' * THE INITIAL CONDITIONS AT b ARE  * '
              WRITE (*,FMT=*) ' *                                  * '
              IF (NCB.LE.2 .OR. NCB.EQ.6) THEN
                  WRITE (*,FMT=*)
     +              ' *    y(b)=beta1, py''(b)=beta2      *'
              ELSE
                  WRITE (*,FMT=*)
     +              ' *  [y,u](b)=beta1, [y,v](b)=beta2  * '
              END IF
              WRITE (*,FMT=*) ' *                                  * '
              WRITE (*,FMT=*) ' * WHERE THE CONSTANTS beta1, beta2 * '
              WRITE (*,FMT=*) ' *    MAY BE CHOSEN ARBITRARILY.    * '
              WRITE (*,FMT=*) ' *                                  * '
              WRITE (*,FMT=*) ' * CHOOSE beta1,beta2: (h?)         * '
              WRITE (*,FMT=*) ' ************************************ '
              WRITE (*,FMT=*)
              WRITE (*,FMT=*) ' beta1,beta2 = '
              READ (*,FMT=9010) CHANS
              READ (CHANS,FMT=9020) HQ
              IF (HQ.EQ.'q' .OR. HQ.EQ.'Q') STOP
              IF (HQ.EQ.'h' .OR. HQ.EQ.'H') THEN
                  CALL HELP(12)
                  GO TO 400
              END IF
              CALL LSTDIR(CHANS,I,ICOL)
              I1 = ICOL(1) - 1
              FMAT = FMT2(I1)
              READ (CHANS,FMT=FMAT) BETA1,BETA2
              B1 = BETA2
              B2 = -BETA1
          END IF
      END IF
      WRITE (*,FMT=*)
      WRITE (*,FMT=*)
      WRITE (*,FMT=*)
      WRITE (*,FMT=*)
      RETURN
C
   20 CONTINUE
  250 CONTINUE
      WRITE (*,FMT=*)
      WRITE (*,FMT=*) ' ********************************************* '
      WRITE (*,FMT=*) ' * DO YOU WANT TO COMPUTE                    * '
      WRITE (*,FMT=*) ' *                                           * '
      WRITE (*,FMT=*) ' *    (1) A SINGLE EIGENVALUE                * '
      WRITE (*,FMT=*) ' *                                           * '
      WRITE (*,FMT=*) ' *    (2) A SERIES OF EIGENVALUES ?          * '
      WRITE (*,FMT=*) ' *                                           * '
      WRITE (*,FMT=*) ' * ENTER THE NUMBER OF YOUR CHOICE: (h?)     * '
      WRITE (*,FMT=*) ' ********************************************* '
      WRITE (*,FMT=*)
      READ (*,FMT=9010) CHANS
      READ (CHANS,FMT=9020) HQ
      IF (HQ.EQ.'q' .OR. HQ.EQ.'Q') STOP
      IF (HQ.EQ.'h' .OR. HQ.EQ.'H') THEN
          CALL HELP(13)
          GO TO 250
      END IF
      YEH = HQ .EQ. '1' .OR. HQ .EQ. '2'
      IF (.NOT.YEH) GO TO 250
      READ (CHANS,FMT='(I32)') NANS
      RESP = NANS
      RETURN
   21 CONTINUE
  260 CONTINUE
      WRITE (*,FMT=*) ' ****************************************** '
      WRITE (*,FMT=*) ' * INPUT NUMEIG, EIG, TOL: (h?)           * '
      WRITE (*,FMT=*) ' ****************************************** '
      WRITE (*,FMT=*)
      WRITE (*,FMT=*) ' NUMEIG,EIG,TOL = '
      READ (*,FMT=9010) CHANS
      READ (CHANS,FMT=9020) HQ
      IF (HQ.EQ.'q' .OR. HQ.EQ.'Q') STOP
      IF (HQ.EQ.'h' .OR. HQ.EQ.'H') THEN
          CALL HELP(14)
          GO TO 260
      END IF
      EIG = 0.0D0
      TOL = 0.0D0
      CALL LSTDIR(CHANS,I,ICOL)
      I1 = ICOL(1) - 1
      IF (I.EQ.3) THEN
          I2 = ICOL(2) - ICOL(1) - 1
          FMAT = '(I'//COL(I1)//',1X,F'//COL(I2)//'.0,1X,F'//
     +           COL(30-I1-I2)//'.0)'
          READ (CHANS,FMT=FMAT) NUMEIG,EIG,TOL
      ELSE IF (I.EQ.2) THEN
          FMAT = '(I'//COL(I1)//',1X,F'//COL(31-I1)//'.0)'
          READ (CHANS,FMT=FMAT) NUMEIG,EIG
      ELSE
          FMAT = '(I'//COL(I1)//')'
          READ (CHANS,FMT=FMAT) NUMEIG
      END IF
      IF (RITE) WRITE (T22,FMT=*) ' NUMEIG,EIG,TOL = ',NUMEIG,EIG,TOL
      WRITE (*,FMT=*)
      WRITE (*,FMT=*)
      WRITE (*,FMT=*)
      WRITE (*,FMT=*)
      WRITE (*,FMT=*)
      RETURN
C
   22 CONTINUE
  270 CONTINUE
      WRITE (*,FMT=*) ' *******************************'//FILLA
      WRITE (*,FMT=*) ' *                              '//FILLB
      WRITE (*,FMT=*) ' * PLOT EIGENFUNCTION ? (Y/N) (h?)'//
     +  '                *'
      WRITE (*,FMT=*) ' *******************************'//FILLA
      WRITE (*,FMT=*)
      READ (*,FMT=9010) CHANS
      READ (CHANS,FMT=9020) HQ
      IF (HQ.EQ.'q' .OR. HQ.EQ.'Q') STOP
      IF (HQ.EQ.'h' .OR. HQ.EQ.'H') THEN
          CALL HELP(16)
          GO TO 270
      END IF
      YEH = HQ .EQ. 'y' .OR. HQ .EQ. 'Y' .OR. HQ .EQ. 'n' .OR.
     +      HQ .EQ. 'N'
      IF (.NOT.YEH) GO TO 270
      PEIGF = (HQ.EQ.'y' .OR. HQ.EQ.'Y')
      IF (.NOT.PEIGF) RESP = 0
      RETURN
C
   23 CONTINUE
  280 CONTINUE
      WRITE (*,FMT=*) ' ****************************************** '
      WRITE (*,FMT=*) ' * INPUT numeig1, numeig2, TOL (h?)       * '
      WRITE (*,FMT=*) ' ****************************************** '
      WRITE (*,FMT=*)
      WRITE (*,FMT=*) ' numeig1,numeig2,TOL = '
      READ (*,FMT=9010) CHANS
      READ (CHANS,FMT=9020) HQ
      IF (HQ.EQ.'q' .OR. HQ.EQ.'Q') STOP
      IF (HQ.EQ.'h' .OR. HQ.EQ.'H') THEN
          CALL HELP(14)
          GO TO 280
      END IF
      TOLL = 0.0D0
c           TOLL = 1.E-5
      CALL LSTDIR(CHANS,I,ICOL)
      I1 = ICOL(1) - 1
      IF (I.EQ.3) THEN
          I2 = ICOL(2) - ICOL(1) - 1
          FMAT = '(I'//COL(I1)//',1X,I'//COL(I2)//',1X,F'//
     +           COL(30-I1-I2)//'.0)'
          READ (CHANS,FMT=FMAT) NEIG1,NEIG2,TOLL
      ELSE
          FMAT = '(I'//COL(I1)//',1X,I'//COL(31-I1)//')'
          READ (CHANS,FMT=FMAT) NEIG1,NEIG2
      END IF
      IF (NEIG1.NE.NEIG2) PEIGF = .FALSE.
      IF (RITE) WRITE (T22,FMT=*) ' numeig1,numeig2,TOL = ',NEIG1,NEIG2,
     +    TOLL
      RETURN
C
   24 CONTINUE
  253 CONTINUE
      WRITE (*,FMT=*)
      WRITE (*,FMT=*) ' ********************************************* '
      WRITE (*,FMT=*) ' * DO YOU WANT TO COMPUTE                    * '
      WRITE (*,FMT=*) ' *                                           * '
      WRITE (*,FMT=*) ' *    (1) A SINGLE EIGENVALUE                * '
      WRITE (*,FMT=*) ' *                                           * '
      WRITE (*,FMT=*) ' *    (2) A SERIES OF EIGENVALUES ?          * '
      WRITE (*,FMT=*) ' *                                           * '
      WRITE (*,FMT=*) ' * ENTER THE NUMBER OF YOUR CHOICE: (h?)     * '
      WRITE (*,FMT=*) ' ********************************************* '
      WRITE (*,FMT=*)
      READ (*,FMT=9010) CHANS
      READ (CHANS,FMT=9020) HQ
      IF (HQ.EQ.'q' .OR. HQ.EQ.'Q') STOP
      IF (HQ.EQ.'h' .OR. HQ.EQ.'H') THEN
          CALL HELP(13)
          GO TO 253
      END IF
      YEH = HQ .EQ. '1' .OR. HQ .EQ. '2'
      IF (.NOT.YEH) GO TO 253
      READ (CHANS,FMT='(I32)') NANS
      RESP = NANS
      RETURN
C
   25 CONTINUE
  330 CONTINUE
      WRITE (*,FMT=*)
      WRITE (*,FMT=*) ' ********************************************* '
      WRITE (*,FMT=*) ' * INPUT NUMEIG,TOL: (h?)                    * '
      WRITE (*,FMT=*) ' ********************************************* '
      WRITE (*,FMT=*)
      WRITE (*,FMT=*) ' NUMEIG,TOL = '
      READ (*,FMT=9010) CHANS
      READ (CHANS,FMT=9020) HQ
      IF (HQ.EQ.'q' .OR. HQ.EQ.'Q') STOP
      IF (HQ.EQ.'h' .OR. HQ.EQ.'H') THEN
          CALL HELP(17)
          GO TO 330
      END IF
      TOL = 0.0D0
      CALL LSTDIR(CHANS,I,ICOL)
      I1 = ICOL(1) - 1
      IF (I.EQ.2) THEN
          FMAT = '(I'//COL(I1)//',1X,F'//COL(31-I1)//'.0)'
          READ (CHANS,FMT=FMAT) NUMEIG,TOL
      ELSE
          FMAT = '(I'//COL(I1)//')'
          READ (CHANS,FMT=FMAT) NUMEIG
      END IF
      IF (RITE) WRITE (T22,FMT=*) ' NUMEIG,TOL = ',NUMEIG,TOL
      IF (NUMEIG.LT.0 .AND. (NCA.NE.4.AND.NCB.NE.4)) THEN
          WRITE (*,FMT=*) ' NUMEIG MUST BE .GE. 0 '
          GO TO 330
      END IF
      WRITE (*,FMT=*)
      WRITE (*,FMT=*)
      WRITE (*,FMT=*)
      WRITE (*,FMT=*)
      WRITE (*,FMT=*)
      WRITE (*,FMT=*)
      RETURN
C
   26 CONTINUE
  340 CONTINUE
      WRITE (*,FMT=*)
      WRITE (*,FMT=*)
     +  ' *************************************************'
      WRITE (*,FMT=*)
     +  ' * PLOT EIGENFUNCTION ? (Y/N) (h?)               *'
      WRITE (*,FMT=*)
     +  ' *************************************************'
      WRITE (*,FMT=*)
      READ (*,FMT=9010) CHANS
      READ (CHANS,FMT=9020) HQ
      IF (HQ.EQ.'q' .OR. HQ.EQ.'Q') STOP
      IF (HQ.EQ.'h' .OR. HQ.EQ.'H') THEN
          CALL HELP(16)
          GO TO 340
      END IF
      YEH = HQ .EQ. 'y' .OR. HQ .EQ. 'Y' .OR. HQ .EQ. 'n' .OR.
     +      HQ .EQ. 'N'
      IF (.NOT.YEH) GO TO 340
      PEIGF = HQ .EQ. 'y' .OR. HQ .EQ. 'Y'
      IF (.NOT.PEIGF) RESP = 0
      IF (PEIGF) THEN
          THA = ATAN2(A2,-A1)
          IF (THA.LT.0.0D0) THA = THA + PI
          THB = ATAN2(B2,-B1)
          IF (THB.LE.0.0D0) THB = THB + PI
          A1 = COS(THA)
          A2 = -SIN(THA)
          B1 = COS(THB)
          B2 = -SIN(THB)
          R1 = K11*SIN(THA) + K12*COS(THA)
          R2 = K21*SIN(THA) + K22*COS(THA)
          RHO = SQRT(R1**2+R2**2)
          B1 = RHO*B1
          B2 = RHO*B2
          NIVP = 2
          NEND = 3
          SLFUN(1) = 0.0D0
          SLFUN(2) = -1.0D0
          SLFUN(3) = THA
          SLFUN(5) = 1.0D0
          SLFUN(6) = THB + NUMEIG*PI
          SLFUN(8) = 0.00001D0
          SLFUN(9) = 1.0D0
      END IF
      RETURN
C
   27 CONTINUE
  283 CONTINUE
      WRITE (*,FMT=*) ' ****************************************** '
      WRITE (*,FMT=*) ' * INPUT numeig1, numeig2, TOL (h?)       * '
      WRITE (*,FMT=*) ' ****************************************** '
      WRITE (*,FMT=*)
      WRITE (*,FMT=*) ' numeig1,numeig2,TOL = '
      READ (*,FMT=9010) CHANS
      READ (CHANS,FMT=9020) HQ
      IF (HQ.EQ.'q' .OR. HQ.EQ.'Q') STOP
      IF (HQ.EQ.'h' .OR. HQ.EQ.'H') THEN
          CALL HELP(14)
          GO TO 283
      END IF
      CALL LSTDIR(CHANS,I,ICOL)
      I1 = ICOL(1) - 1
      IF (I.EQ.3) THEN
          I2 = ICOL(2) - ICOL(1) - 1
          FMAT = '(I'//COL(I1)//',1X,I'//COL(I2)//',1X,F'//
     +           COL(30-I1-I2)//'.0)'
          READ (CHANS,FMT=FMAT) NEIG1,NEIG2,TOLL
      ELSE
          FMAT = '(I'//COL(I1)//',1X,I'//COL(31-I1)//')'
          READ (CHANS,FMT=FMAT) NEIG1,NEIG2
      END IF
      IF (NEIG1.NE.NEIG2) PEIGF = .FALSE.
      IF (RITE) WRITE (T22,FMT=*) ' numeig1,numeig2,TOL = ',NEIG1,NEIG2,
     +    TOLL
      RETURN
C
   28 CONTINUE
      WRITE (*,FMT=*)
     +  ' **************************************************'
      DO 285 I = 1,ILAST
          NUMEIG = NEIG1 + (I-1)
          EIG = EES(I)
          TOL = TTS(I)
          IFLAG = IIS(I)
          IF (IFLAG.LE.2) THEN
              WRITE (*,FMT=9055) NUMEIG,EIG
              WRITE (*,FMT=9050) TOL,IFLAG

              IF (RITE) THEN
                  WRITE (T22,FMT=9055) NUMEIG,EIG
                  WRITE (T22,FMT=9050) TOL,IFLAG
              END IF

          ELSE IF (IFLAG.EQ.4) THEN
              WRITE (*,FMT=9045) NUMEIG
              IF (RITE) WRITE (T22,FMT=9045) NUMEIG
          ELSE
              WRITE (*,FMT=9040) NUMEIG,IFLAG
              IF (RITE) WRITE (T22,FMT=9040) NUMEIG,IFLAG
              IF (IFLAG.EQ.3) THEN
                IF(NLAST.GE.0) THEN
                  WRITE (*,FMT=9230) NLAST
                  IF (RITE) WRITE (T22,FMT=9230) NLAST
                ELSE IF(.NOT.(NCA.EQ.4 .OR. NCB.EQ.4)) THEN
                  WRITE (*,FMT=9235) 
                  IF (RITE) WRITE (T22,FMT=9235)
                ELSE
                  WRITE (*,FMT=9240)
                  IF (RITE) WRITE (T22,FMT=9240)
                ENDIF
              END IF
          END IF

  285 CONTINUE
      WRITE (*,FMT=*)
     +  ' **************************************************'

      WRITE (*,FMT=*)
      IF (RITE) WRITE (T22,FMT=*)
      IF (JFLAG.EQ.1) WRITE (*,FMT=9260)
      IF (JFLAG.EQ.1 .AND. RITE) WRITE (T22,FMT=9260)
      IF (JFLAG.EQ.2) WRITE (*,FMT=9250) SLF9
      IF (JFLAG.EQ.2 .AND. RITE) WRITE (T22,FMT=9250) SLF9
      IF (RITE) WRITE (T22,FMT=*) ' *******************************'//
     +    FILLA
      WRITE (*,FMT=*)
      RETURN
C
   29 CONTINUE
      WRITE (*,FMT=*)
     +  ' **************************************************'
      DO 303 I = 1,ILAST
          NUMEIG = NEIG1 + I - 1
          EIG = EES(I)
          TOL = TTS(I)
          IFLAG = IIS(I)
          DUBBLE = .FALSE.
          IF (IFLAG.EQ.2) THEN
              DUBBLE = .TRUE.
              IFLAG = 1
          END IF
          WRITE (*,FMT=9055) NUMEIG,EIG
          WRITE (*,FMT=9050) TOL,IFLAG
          IF (DUBBLE) THEN
              WRITE (*,FMT=*)
     +          ' * This eigenvalue appears to be double.          *'
          END IF
          IF (RITE) THEN
              WRITE (T22,FMT=9055) NUMEIG,EIG
              WRITE (T22,FMT=9050) TOL,IFLAG
              IF (DUBBLE) WRITE (T22,FMT=*
     +            ) ' * This eigenvalue appears to be double. * '
          END IF
          WRITE (*,FMT=*)
  303 CONTINUE
      WRITE (*,FMT=*)
     +  ' **************************************************'
      WRITE (*,FMT=*)
      WRITE (*,FMT=*)
C
      IF (RITE) WRITE (T22,FMT=*) ' *******************************'//
     +    FILLA
      RETURN
C
   30 CONTINUE
      WRITE (*,FMT=*) ' Press any key to continue. '
      READ (*,FMT=9010) CHANS
      WRITE (*,FMT=*)
     +  ' *********************************************** '
      WRITE (*,FMT=*)
     +  ' * WHAT WOULD YOU LIKE TO DO NOW ?             * '
      WRITE (*,FMT=*)
     +  ' *                                             * '
      WRITE (*,FMT=*)
     +  ' *    (1)  SAME EIGENVALUE PROBLEM, DIFFERENT  * '
      WRITE (*,FMT=*)
     +  ' *           NUMEIG, EIG, OR TOL               * '
      WRITE (*,FMT=*)
     +  ' *                                             * '
      WRITE (*,FMT=*)
     +  ' *    (2)  SAME EIGENVALUE PROBLEM, SAME (a,b) * '
      WRITE (*,FMT=*)
     +  ' *           AND p,q,w,u,v BUT DIFFERENT       * '
      WRITE (*,FMT=*)
     +  ' *           BOUNDARY CONDITIONS A1,A2,B1,B2   * '
      WRITE (*,FMT=*)
     +  ' *                                             * '
      WRITE (*,FMT=*)
     +  ' *    (3)  INTERVAL CHANGE, PROBLEM RESTART    * '
      WRITE (*,FMT=*)
     +  ' *                                             * '
      WRITE (*,FMT=*)
     +  ' *    (4)  AN INITIAL VALUE PROBLEM            * '
      WRITE (*,FMT=*)
     +  ' *                                             * '
      WRITE (*,FMT=*)
     +  ' *    (5)  QUIT                                * '
      WRITE (*,FMT=*)
     +  ' *                                             * '
      WRITE (*,FMT=*)
     +  ' * ENTER THE NUMBER OF YOUR CHOICE: (h?)       * '
      WRITE (*,FMT=*)
     +  ' *********************************************** '
      WRITE (*,FMT=*)
      READ (*,FMT=9010) CHANS
      READ (CHANS,FMT=9020) HQ
      IF (HQ.EQ.'q' .OR. HQ.EQ.'Q') STOP
      IF (HQ.EQ.'h' .OR. HQ.EQ.'H') THEN
          CALL HELP(3)
          GO TO 30
      END IF
      YEH = HQ .EQ. '1' .OR. HQ .EQ. '2' .OR. HQ .EQ. '3' .OR.
     +      HQ .EQ. '4' .OR. HQ .EQ. '5'
      IF (.NOT.YEH) GO TO 30
      READ (CHANS,FMT='(I32)') NANS
      RESP = NANS
      RETURN
C
   31 CONTINUE
      WRITE (*,FMT=*) ' Press any key to continue. '
      READ (*,FMT=9010) CHANS
      WRITE (*,FMT=*)
     +  ' *********************************************** '
      WRITE (*,FMT=*)
     +  ' *    (1)  SAME INITIAL VALUE PROBLEM,         * '
      WRITE (*,FMT=*)
     +  ' *           DIFFERENT LAMBDA                  * '
      WRITE (*,FMT=*)
     +  ' *                                             * '
      WRITE (*,FMT=*)
     +  ' *    (2)  NEW INITIAL VALUE PROBLEM           * '
      WRITE (*,FMT=*)
     +  ' *                                             * '
      WRITE (*,FMT=*)
     +  ' *    (3)  INTERVAL CHANGE, PROBLEM RESTART    * '
      WRITE (*,FMT=*)
     +  ' *                                             * '
      WRITE (*,FMT=*)
     +  ' *    (4)  AN EIGENVALUE PROBLEM               * '
      WRITE (*,FMT=*)
     +  ' *                                             * '
      WRITE (*,FMT=*)
     +  ' *    (5)  QUIT                                * '
      WRITE (*,FMT=*)
     +  ' *                                             * '
      WRITE (*,FMT=*)
     +  ' * ENTER THE NUMBER OF YOUR CHOICE: (h?)       * '
      WRITE (*,FMT=*)
     +  ' *********************************************** '
      WRITE (*,FMT=*)
      READ (*,FMT=9010) CHANS
      READ (CHANS,FMT=9020) HQ
      IF (HQ.EQ.'q' .OR. HQ.EQ.'Q') STOP
      IF (HQ.EQ.'h' .OR. HQ.EQ.'H') THEN
          CALL HELP(3)
          GO TO 31
      END IF
      YEH = HQ .EQ. '1' .OR. HQ .EQ. '2' .OR. HQ .EQ. '3' .OR.
     +      HQ .EQ. '4' .OR. HQ .EQ. '5'
      IF (.NOT.YEH) GO TO 31
      READ (CHANS,FMT='(I32)') NANS
      RESP = NANS
      RETURN
C
   32 CONTINUE
  410 CONTINUE
      WRITE (*,FMT=*)
      WRITE (*,FMT=*)
     +  ' ************************************************ '
      WRITE (*,FMT=*)
     +  ' * WHAT VALUE SHOULD BE USED FOR THE            * '
      WRITE (*,FMT=*)
     +  ' *   EIGENPARAMETER, EIG ?                      * '
      WRITE (*,FMT=*)
     +  ' *                                              * '
      WRITE (*,FMT=*)
     +  ' * INPUT EIG = (h?)                             * '
      WRITE (*,FMT=*)
     +  ' ************************************************ '
      WRITE (*,FMT=*)
      WRITE (*,FMT=*) ' EIG = '
      READ (*,FMT=9010) CHANS
      READ (CHANS,FMT=9020) HQ
      IF (HQ.EQ.'q' .OR. HQ.EQ.'Q') STOP
      IF (HQ.EQ.'h' .OR. HQ.EQ.'H') THEN
          CALL HELP(12)
          GO TO 410
      END IF
      READ (CHANS,FMT='(F32.0)') EIG
      PEIGF = .TRUE.
      RETURN
C
   33 CONTINUE
  450 CONTINUE
      WRITE (*,FMT=*)
     +  ' ****************************************************'
      WRITE (*,FMT=*)
     +  ' * WHICH FUNCTION DO YOU WANT TO PLOT ?             *'
      WRITE (*,FMT=*)
     +  ' *                                                  *'
      WRITE (*,FMT=*) XC(1)
      WRITE (*,FMT=*) XC(2)
      WRITE (*,FMT=*) XC(3)
      WRITE (*,FMT=*) XC(4)
      WRITE (*,FMT=*) XC(5)
      WRITE (*,FMT=*) XC(6)
      WRITE (*,FMT=*)
     +  ' *                                                  *'
      WRITE (*,FMT=*)
     +  ' * ENTER THE NUMBER OF YOUR CHOICE: (h?)            *'
      WRITE (*,FMT=*)
     +  ' ****************************************************'
      WRITE (*,FMT=*)
      READ (*,FMT=9010) CHANS
      READ (CHANS,FMT=9020) HQ
      IF (HQ.EQ.'q' .OR. HQ.EQ.'Q') RESP = 0
      IF (HQ.EQ.'h' .OR. HQ.EQ.'H') THEN
          CALL HELP(16)
          GO TO 450
      END IF
      YEH = HQ .EQ. '1' .OR. HQ .EQ. '2' .OR. HQ .EQ. '3' .OR.
     +      HQ .EQ. '4' .OR. HQ .EQ. '5' .OR. HQ .EQ. '6'
      IF (.NOT.YEH) GO TO 450
      READ (CHANS,FMT='(I32)') NF
  470 CONTINUE
      WRITE (*,FMT=*)
     +  ' ****************************************************'
      WRITE (*,FMT=*)
     +  ' * WHICH DO YOU WANT AS THE INDEPENDENT VARIABLE ?  *'
      WRITE (*,FMT=*)
     +  ' *                                                  *'
      WRITE (*,FMT=*) XC(7)
      WRITE (*,FMT=*) XC(8)
      WRITE (*,FMT=*)
     +  ' *                                                  *'
      WRITE (*,FMT=*)
     +  ' * ENTER THE NUMBER OF YOUR CHOICE: (h?)            *'
      WRITE (*,FMT=*)
     +  ' ****************************************************'
      WRITE (*,FMT=*)
C
      READ (*,FMT=9010) CHANS
      READ (CHANS,FMT=9020) HQ
      IF (HQ.EQ.'q' .OR. HQ.EQ.'Q') RESP = 0
      IF (HQ.EQ.'h' .OR. HQ.EQ.'H') THEN
          CALL HELP(16)
          GO TO 470
      END IF
      YEH = HQ .EQ. '1' .OR. HQ .EQ. '2'
      IF (.NOT.YEH) GO TO 470
      READ (CHANS,FMT='(I32)') NV
      RETURN
C
   34 CONTINUE
  480 CONTINUE
      WRITE (*,FMT=*)
      WRITE (*,FMT=*) ' DO YOU WANT TO SAVE THE PLOT FILE ? (Y/N) (h?)'
      READ (*,FMT=9010) CHANS
      READ (CHANS,FMT=9020) HQ
      IF (HQ.EQ.'q' .OR. HQ.EQ.'Q') STOP
      IF (HQ.EQ.'h' .OR. HQ.EQ.'H') THEN
          CALL HELP(16)
          GO TO 480
      END IF
      YEH = HQ .EQ. 'y' .OR. HQ .EQ. 'Y' .OR. HQ .EQ. 'n' .OR.
     +      HQ .EQ. 'N'
      IF (.NOT.YEH) GO TO 480
      RESP = 0
      IF (HQ.EQ.'y' .OR. HQ.EQ.'Y') RESP = 1
      RETURN
C
   35 CONTINUE
  490 CONTINUE
      WRITE (*,FMT=*) ' PLOT ANOTHER FUNCTION ? (Y/N) (h?)'
      READ (*,FMT=9010) CHANS
      READ (CHANS,FMT=9020) HQ
      IF (HQ.EQ.'q' .OR. HQ.EQ.'Q') STOP
      IF (HQ.EQ.'h' .OR. HQ.EQ.'H') THEN
          CALL HELP(16)
          GO TO 490
      END IF
      YEH = HQ .EQ. 'y' .OR. HQ .EQ. 'Y' .OR. HQ .EQ. 'n' .OR.
     +      HQ .EQ. 'N'
      IF (.NOT.YEH) GO TO 490
      RESP = 0
      IF (HQ.EQ.'y' .OR. HQ.EQ.'Y') RESP = 1
      RETURN
C
 9010 FORMAT (A32)
 9020 FORMAT (A1)
 9030 FORMAT (A70)
 9040 FORMAT (1X,'   NUMEIG = ',I5,'   IFLAG = ',I3)
 9045 FORMAT (1X,'   NUMEIG = ',I5,' : COMPUTATION FAILED ')
 9050 FORMAT (1X,' * TOL = ',D14.5,2X,' IFLAG = ',I3,'             *')
 9055 FORMAT (1X,' * NUMEIG = ',I5,' EIG = ',D18.9,'        *')
 9070 FORMAT (1X,1X,'*   (',F12.7,',',F12.7,')','                *')
 9080 FORMAT (1X,1X,'*   (',F12.7,',',A12,')','                *')
 9090 FORMAT (1X,1X,'*   (',A12,',',F12.7,')','                *')
 9100 FORMAT (1X,1X,'*   (',A12,',',A12,')','                *')
 9110 FORMAT (1X,2A39)
 9230 FORMAT (1X,' * THERE SEEMS TO BE NO EIGENVALUE OF INDEX       *',
     +       /,1X,' * GREATER THAN',I5,'                              *'
     +       )
 9235 FORMAT (1X,' * THERE SEEM TO BE NO EIGENVALUES               *')
 9240 FORMAT (1X,' * THERE SEEMS TO BE NO EIGENVALUE OF THIS INDEX *')
 9250 FORMAT (1X,'* THERE SEEMS TO BE CONTINUOUS SPECTRUM BEGINNING*',/,
     +       1X,'* AT ABOUT',1P,D8.1,'                               *')
 9260 FORMAT (1X,'* THERE SEEMS TO BE NO CONTINUOUS SPECTRUM       *')
      END
C
      SUBROUTINE STAGE(A,B,INTAB,P0ATA,QFATA,P0ATB,QFATB,A1,A2,B1,B2,
     +                 EIG,IFLAG,SLFUN,NCA,NCB,NIVP,NEND)

C
C     THE FOLLOWING CALL SETS THE STAGE IN sleign -- I.E., SAMPLES
C     THE COEFFICIENTS AND SETS THE INITIAL INTERVAL.
C
C     .. Scalar Arguments ..
      REAL A,A1,A2,B,B1,B2,EIG,P0ATA,P0ATB,QFATA,QFATB
      INTEGER IFLAG,INTAB,NCA,NCB,NEND,NIVP
C     ..
C     .. Array Arguments ..
      REAL SLFUN(9)
C     ..
C     .. Scalars in Common ..
      REAL AA,BB,DTHDAA,DTHDBB,HPI,PI,TMID,TWOPI
      INTEGER MDTHZ
      LOGICAL ADDD
C     ..
C     .. Local Scalars ..
      REAL ALFA,BETA,TOL
      INTEGER NUMEIG
C     ..
C     .. External Subroutines ..
      EXTERNAL SLEIGN
C     ..
C     .. Common blocks ..
      COMMON /PIE/PI,TWOPI,HPI
      COMMON /TDATA/AA,TMID,BB,DTHDAA,DTHDBB,MDTHZ,ADDD
C     ..
      NUMEIG = 0
      TOL = .001D0
      CALL SLEIGN(A,B,INTAB,P0ATA,QFATA,P0ATB,QFATB,A1,A2,B1,B2,NUMEIG,
     +            EIG,TOL,IFLAG,-1,SLFUN,NCA,NCB)
      IF (NIVP.EQ.1 .AND. NEND.EQ.1) TMID = BB
      IF (NIVP.EQ.1 .AND. NEND.EQ.2) TMID = AA
C
C     ACTUALLY, WE MAY HAVE TO AVOID AA OR BB IN SOME CASES,
C     WHICH WILL BE TAKEN CARE OF LATER.
C
      ALFA = 0.0D0
      IF ((NIVP.EQ.1.AND.NEND.EQ.1.AND.NCA.LE.2) .OR.
     +    NIVP.EQ.2) ALFA = SLFUN(3)
      BETA = PI
      IF ((NIVP.EQ.1.AND.NEND.EQ.2.AND.NCB.LE.2) .OR.
     +    NIVP.EQ.2) BETA = SLFUN(6)
      SLFUN(3) = ALFA
      SLFUN(6) = BETA
      SLFUN(4) = 0.0D0
      SLFUN(7) = 0.0D0
      SLFUN(8) = .01D0
      RETURN
      END
C
      SUBROUTINE LSTDIR(CHANS,I,ICOL)
C     .. Local Scalars ..
      INTEGER J
C     ..
C     .. Scalar Arguments ..
      INTEGER I
      CHARACTER*32 CHANS
C     ..
C     .. Array Arguments ..
      INTEGER ICOL(2)
C     ..
      I = 1
      DO 10 J = 1,32
          IF (CHANS(J:J).EQ.',' .OR. CHANS(J:J).EQ.'/') THEN
              ICOL(I) = J
              IF (CHANS(J:J).EQ.'/') THEN
                  CHANS(J:J) = ' '
                  RETURN
              END IF
              I = I + 1
          END IF
   10 CONTINUE
      RETURN
      END
      CHARACTER*32 FUNCTION FMT2(I1)
C     .. Local Arrays ..
      CHARACTER*2 COL(32)
C     ..
C
C     .. Scalar Arguments ..
      INTEGER I1
C     ..
C     .. Data statements ..
      DATA COL/'01','02','03','04','05','06','07','08','09','10','11',
     +     '12','13','14','15','16','17','18','19','20','21','22','23',
     +     '24','25','26','27','28','29','30','31','32'/
C     ..
      FMT2 = '(F'//COL(I1)//'.0,1X,F'//COL(31-I1)//'.0)'
      RETURN
      END
C
      SUBROUTINE CONVT(CHIN,CHINN)
C   THIS PROGRAM CONVERTS AN 8 CHARACTER STRING WITHOUT A ':'
C     TO A BLANK STRING.  BUT IF IT HAS A ':', THEN IT KEEPS THE
C     STRING UP TO AND INCLUDING THE ':' .
C   THE INPUT IS STRING CHIN, AND THE OUTPUT IS CHINN.

C     .. Scalar Arguments ..
      CHARACTER*8 CHIN,CHINN
C     ..
C     .. Local Scalars ..
      INTEGER I,J,K
C     ..
      CHINN = ' '
      J = 0
      K = 0
      DO 10 I = 1,8
          IF (CHIN(I:I).EQ.' ' .AND. J.EQ.0) K = I
          IF (CHIN(I:I).EQ.':') J = I
   10 CONTINUE
      IF (J.GE.1 .AND. J.LE.8) THEN
          DO 20 I = 1,J - K
              CHINN(I:I) = CHIN(I+K:I+K)
   20     CONTINUE
      ELSE
          CHINN = ' '
      END IF
      RETURN
      END
C
      SUBROUTINE AUTO
C  THIS SUBROUTINE IS, IN EFFECT, A VERY SPECIAL KIND OF DRIVER
C    FOR SLEIGN2. IT ENABLES ONE TO BYPASS THE QUESTION
C    AND ANSWER FORMAT IN DRIVE (FOR EXPERIENCED USERS ONLY!).
C    INSTEAD OF ENTERING THE NEEDED INPUT FOR DRIVE FROM THE
C    KEYBOARD, A USER CAN INSTEAD CREATE A VERY BRIEF FILE,
C    CALLED auto.in, IN THE SAME DIRECTORY, WHICH CONTAINS ALL
C    THE DATA THAT WOULD OTHERWISE BE ENTERED FROM THE
C    KEYBOARD.
C
C    ONCE SUCH A FILE HAS BEEN CREATED, THE USER BEGINS AS USUAL,
C    TYPING THE NAME OF THE EXECUTABLE, AS IN

C       XAMPLES.X  <ENTER>

C    FOR EXAMPLE.  THIS WOULD CAUSE THE PROMPT

C        WOULD YOU LIKE AN OVERVIEW OF HELP?(Y/N)(h?)

C    TO APPEAR, AS USUAL.  BUT INSTEAD OF REPLYING TO THE
C    QUESTION WITH y, OR n, OR h, ONE SIMPLY TYPES IN THE
C    RESPONSE
C       a  <ENTER>
C    (The "a" here stands for "automatic" operation of SLEIGN2.);
C    and at this point the computation of the requested
C    eigenvalues(s) proceeds without further action by the user,
C    taking the needed data from the auto.in file instead.
C         The construction of the file "auto.in" consists of simply
C    defining a number of "KEYWORDS", each on a separate line, which
C    together constitute a complete set of input parameters defining
C    the eigenvalue problem to be solved. (The Differential Equation
C    coefficients have, presumably, been already defined in either
C    XAMPLES.F, or BLOGGS.F, or.. .)  The order in which the
C    needed keywords are defined is of no importance.
C         -------------------------------------------
C    The KEYWORDS, all of which end in a colon and must be
C    followed by at least one space, are:

C    a:  -- The left endpoint of the interval (a,b);
C         Value is any real number;
C         Default value is -infinity.
C    b:  -- The right endpoint of the interval (a,b);
C         Value is any real number;
C         Default value is +infinity.
C    classa:  -- The endpoint classification of a;
C              Value is one of { r, wr, lcno, lco, lp, d }.
C    classb:  -- The endpoint classification of b:
C              Value is one of { r, wr, lcno, lco, lp, d }.
C    bca: -- Boundary Condition for the endpoint a;
C          Value is either d (for Dirichlet),
C                          n (for Neuman),
C                       or two real numbers A1, A2 .
C    bcb: -- Boundary Condition for the endpoint b;
C          Value is either d (for Dirichlet),
C                          n (for Neuman),
C                       or two real numbers B1, B2 .
C    bcc:  -- Coupled Boundary Condition;
C           Value is either p (for Periodic),
C                           s (for Semi-Periodic),
C                       or five real numbers
C                          alpha, k11, k12, k21, k22 .
C    numeig:  -- Index (or range of Indices) of desired Eigenvalue;
C              Value is an integer N1, or pair of integers N1,N2 .
C    param:  -- Parameter(s) appropriate for the problem;
C             Value is one or two real numbers, param1, param2 .
C    np:  --  Problem Number;
C           (Appropriate only if one of the Differential Equations
C             in XAMPLES.F is being used.);
C           Value is an integer from 1 to 32 .
C    output:  -- Name of output file;
C              Value is a character string, the name of the output
C                file where the results of the computation are to
C                be written;
C              Default value is "auto.out" .
C    end:  -- Last line of file "auto.in" ; no value set.

C    again:  --  Terminates the input for one eigenvalue problem and
C                begins the input for another ; no value set.

C Although the KEYWORDS used to define any one problem can be defined
C in any order, there are a few rules to be observed:  Namely,
C Only those KEYWORDS whose values are necessary need be mentioned.
C Any KEYWORD definition remains in effect until redefined;
C To erase a previous definition of a KEYWORD, redefine it to
C have the value "null". (For instance, if one problem has
C the endpoint a defined as in "a: 0.0" , and a following
C problem needs to have a undefined (so that a = -Infinity),
C then it would be necessary to set "a: null" .)

C  An example of such a file is exhibited in the box below.

C   ________________________
C  |                        |
C  |output: Bessel.rep      |
C  |np: 2                   |
C  |param: 0.75             |
C  |a: 0.0                  |
C  |b: 1.0                  |
C  |classa: lcno            |
C  |classb: r               |
C  |bca: 1.0,0.0            |
C  |bcb: d                  |
C  |numeig: 2,5             |
C  |end:                    |
C  |________________________|

C     This file (which must be called "auto.in", of course), would
C be suitable for running Bessel's equation in XAMPLES.F .
C Evidently the problem selected is #2 of XAMPLES.F (Bessel's
C equation), with the parameter nu = 0.75, on the interval (0.0,1.0).
C The endpoint a is asserted to be LCNO; endpoint b is R; the
C Boundary Condition at a is defined by A1 = 1.0 & A2 = 0.0;
C Boundary Condition at b is Dirichlet; and the eigenvalues
C lambda(2), lambda(3), lambda(4), lambda(5) are to be computed.

C  --------------------------------------------------------------- C
C     .. Scalars in Common ..
      REAL A,ALPHA,B,BETA,GAMMA,H,KK,L,NU,P1,P2,P3,P4,P5,P6
      INTEGER INTAB,NUMBER,T21,T22,T23,T24,T25
      LOGICAL PR
C     ..
C     .. Local Scalars ..
      REAL A1,A2,ALF,AS,B1,B2,BS,DET,EIG,
     +     K11,K12,K21,K22,ONE,P0ATA,
     +     P0ATB,PARAM1,PARAM2,QFATA,QFATB,TOL,ZER
      INTEGER I,I1,IFLAG,II,INTABS,J,JFLAG,K,LAST,LASTI,M,NCA,NCB,NEIG1,
     +        NEIG2,NP,NUMEIG
      CHARACTER*8 CH8,CHH8
      CHARACTER*32 FMAT
      CHARACTER*62 BCC,BLANK,CH1,CHANS,CHEND,CHIN,CHT,CLASSA,CLASSB,
     +             TAPE25
C     ..
C     .. Local Arrays ..
      REAL PP(5),SLFUN(9)
      INTEGER ICOL(2),VAL(60)
      CHARACTER*2 COL(62)
      CHARACTER*62 CH(60)
C     ..
C     .. External Functions ..
      CHARACTER*32 FMT2
      EXTERNAL FMT2
C     ..
C     .. External Subroutines ..
      EXTERNAL CHAR,CONVT,LST,LSTDIR,PERIO,PQ,SLEIGN,STR2R
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC MIN
C     ..
C     .. Common blocks ..
      COMMON /DATADT/A,B,INTAB
      COMMON /FLAG/NUMBER
      COMMON /PAR/NU,H,KK,L,ALPHA,BETA,GAMMA,P1,P2,P3,P4,P5,P6
      COMMON /PRIN/PR,T21
      COMMON /TAPES/T22,T23,T24,T25
C     ..
C     .. Data statements ..
c

      DATA COL/'01','02','03','04','05','06','07','08','09','10','11',
     +     '12','13','14','15','16','17','18','19','20','21','22','23',
     +     '24','25','26','27','28','29','30','31','32','33','34','35',
     +     '36','37','38','39','40','41','42','43','44','45','46','47',
     +     '48','49','50','51','52','53','54','55','56','57','58','59',
     +     '60','61','62'/
C     ..
c
      OPEN (T21,FILE='test.out')
      OPEN (T24,FILE='auto.in')
c  Set the name of the Output File to 'auto.out' as default.
      TAPE25 = 'auto.out'
c
      ONE = 1.0D0
      ZER = 0.0D0
      ALF = 0.0D0
      CHEND = 'end:'
      BLANK = '                                '
C
C  READ THE auto.in FILE (ASSUMED LESS THAN 50 LINES
C       FOR ANY ONE PROBLEM):

      READ (T24,FMT=9010,END=75) (CH(I),I=1,50)
   75 CONTINUE
C
C  THE VALUES OF THE FUNCTION VAL(*) BELOW INDICATE WHETHER
C    OR NOT CERTAIN KEYWORDS ARE PRESENT IN THE FILE auto.in.
C  WHEN I = 1, VAL(I) REFERS TO THE PROBLEM # IN XAMPLES.F;
C           2, REFERS TO THE POSSIBLE PARAMETERS ;
C              (THERE MAY BE UP TO 2 INPUT PARAMETERS IN
C               THE PARTICULAR DIFFERENTIAL EQUATION.)
C           3, REFERS TO THE ENDPOINT a;
C           4, REFERS TO THE ENDPOINT b;
C           5, REFERS TO THE CLASSIFICATION OF ENDPOINT a;
C           6, REFERS TO THE CLASSIFICATION OF ENDPOINT b;
C           7, REFERS TO A SEPARATED BOUNDARY CONDITION AT a;
C           8, REFERS TO A SEPARATED BOUNDARY CONDITION AT b;
C           9, REFERS TO COUPLED BOUNDARY CONDITIONS;
C          10, REFERS TO THE INDEX(ES) OF THE EIGENVALUES WANTED;
C          11, REFERS TO THE OUTPUT FILE FOR THE RESULTS.
C          12, REFERS TO WHETHER OR NOT THERE IS ANOTHER
C                PROBLEM TO BE RUN.

      DO 15 I = 1,50
          VAL(I) = 0
   15 CONTINUE
      LAST = 0
C
  300 CONTINUE
      A1 = ONE
      A2 = ZER
      B1 = ONE
      B2 = ZER
      DO 60 I = 1,50
          IF (CH(I).NE.BLANK) THEN
              READ (CH(I),FMT=9010) CH1
              READ (CH1,FMT='(A8)') CH8
          ELSE
              CH1 = ' '
              CH8 = ' '
          END IF
          CALL CONVT(CH8,CHH8)
          IF (CHH8.EQ.'np:') THEN
              CALL CHAR(CH(I),K,M,CHANS)
              READ (CHANS,FMT='(I7)') NP
              VAL(1) = 1
          ELSE IF (CHH8.EQ.'param:') THEN
              CALL CHAR(CH(I),K,M,CHANS)
              IF (CHANS.EQ.'null') THEN
                  VAL(2) = 0
              ELSE IF (M.EQ.0) THEN
                  READ (CHANS,FMT='(F32.0)') PARAM1
                  VAL(2) = 1
              ELSE IF (M.EQ.1) THEN
                  CALL LSTDIR(CHANS,II,ICOL)
                  I1 = ICOL(1) - 1
                  FMAT = FMT2(I1)
                  READ (CHANS,FMT=FMAT) PARAM1,PARAM2
                  VAL(2) = 2
              ELSE IF (M.EQ.4) THEN
                  CHIN = CHANS
                  IF (CHIN.NE.' ') THEN
                      DO 400 J = 1,5
                          IF (CHIN.NE.' ') CALL STR2R(CHIN,PP(J))
  400                 CONTINUE
                  END IF
                  P1 = PP(1)
                  P2 = PP(2)
                  P3 = PP(3)
                  P4 = PP(4)
                  P5 = PP(5)
                  P6 = P2 + P3 + 1.0D0 - P4 - P5
C  THESE NUMBERS ARE THE PARAMETERS s,a,b,c,d,e in the Heun eqn.
C    WITH e = a+b+1-c-d
                  VAL(2) = 5
              END IF
          ELSE IF (CHH8.EQ.'a:') THEN
              CALL CHAR(CH(I),K,M,CHANS)
              IF (CHANS.EQ.'null') THEN
                  VAL(3) = 0
              ELSE
                  CALL LST(CHANS,A)
                  VAL(3) = 1
              END IF
          ELSE IF (CHH8.EQ.'b:') THEN
              CALL CHAR(CH(I),K,M,CHANS)
              IF (CHANS.EQ.'null') THEN
                  VAL(4) = 0
              ELSE
                  CALL LST(CHANS,B)
                  VAL(4) = 1
              END IF
          ELSE IF (CHH8.EQ.'classa:') THEN
              CALL CHAR(CH(I),K,M,CHANS)
              READ (CHANS,FMT='(A32)') CLASSA
              VAL(5) = 1
          ELSE IF (CHH8.EQ.'classb:') THEN
              CALL CHAR(CH(I),K,M,CHANS)
              READ (CHANS,FMT='(A32)') CLASSB
              VAL(6) = 1
          ELSE IF (CHH8.EQ.'bca:') THEN
              CALL CHAR(CH(I),K,M,CHANS)
              VAL(7) = 2
              IF (CHANS.EQ.'null') THEN
                  VAL(7) = 0
              ELSE IF (M.EQ.0) THEN
                  IF (CHANS.EQ.'d') THEN
                      A1 = ONE
                      A2 = ZER
                  ELSE IF (CHANS.EQ.'n') THEN
                      A1 = ZER
                      A2 = ONE
                  END IF
              ELSE
                  CALL LSTDIR(CHANS,II,ICOL)
                  I1 = ICOL(1) - 1
                  FMAT = FMT2(I1)
                  READ (CHANS,FMT=FMAT) A1,A2
              END IF
          ELSE IF (CHH8.EQ.'bcb:') THEN
              CALL CHAR(CH(I),K,M,CHANS)
              VAL(8) = 2
              IF (CHANS.EQ.'null') THEN
                  VAL(8) = 0
              ELSE IF (M.EQ.0) THEN
                  IF (CHANS.EQ.'d') THEN
                      B1 = ONE
                      B2 = ZER
                  ELSE IF (CHANS.EQ.'n') THEN
                      B1 = ZER
                      B2 = ONE
                  END IF
              ELSE
                  CALL LSTDIR(CHANS,II,ICOL)
                  I1 = ICOL(1) - 1
                  FMAT = FMT2(I1)
                  READ (CHANS,FMT=FMAT) B1,B2
              END IF
          ELSE IF (CHH8.EQ.'bcc:') THEN
              CALL CHAR(CH(I),K,M,CHANS)
              IF (CHANS.EQ.'null') THEN
                  VAL(9) = 0
              ELSE IF (M.EQ.0) THEN
                  VAL(9) = 1
                  READ (CHANS,FMT='(A32)') BCC
                  IF (BCC.EQ.'p') THEN
                      K11 = ONE
                      K12 = ZER
                      K21 = ZER
                      K22 = ONE
                  ELSE IF (BCC.EQ.'s') THEN
                      K11 = -ONE
                      K12 = ZER
                      K21 = ZER
                      K22 = -ONE
                  END IF
              ELSE
                  CHIN = CHANS
                  DO 405 J = 1,5
                     CALL STR2R(CHIN,PP(J))
 405              CONTINUE
                  ALF = PP(1)
                  K11 = PP(2)
                  K12 = PP(3)
                  K21 = PP(4)
                  K22 = PP(5)
                  VAL(9) = 5
              END IF
          ELSE IF (CHH8.EQ.'numeig:') THEN
              CALL CHAR(CH(I),K,M,CHANS)
              IF (M.EQ.0) THEN
                  READ (CHANS,FMT='(I32)') NUMEIG
                  NEIG1 = NUMEIG
                  NEIG2 = NEIG1
                  VAL(10) = 1
              ELSE
                  CALL LSTDIR(CHANS,II,ICOL)
                  I1 = ICOL(1) - 1
                  FMAT = '(I'//COL(I1)//',1X,I'//COL(61-I1)//')'
                  READ (CHANS,FMT=FMAT) NEIG1,NEIG2
                  VAL(10) = 2
              END IF
          ELSE IF (CHH8.EQ.'output:') THEN
              CLOSE (T25)
              CALL CHAR(CH(I),K,M,CHANS)
              READ (CHANS,FMT='(A32)') TAPE25
              OPEN (T25,FILE=TAPE25)
              VAL(11) = 1
          ELSE IF (CHH8.EQ.'again:') THEN
              LASTI = I
              VAL(12) = VAL(12) + 1
              GO TO 70
          ELSE IF (CHH8.EQ.CHEND) THEN
              GO TO 70
          END IF
   60 CONTINUE
   70 CONTINUE
c
C  WHEN AN ENDPOINT IS LP, NO BOUNDARY CONDITION CAN BE GIVEN;
C    AND COUPLED BOUNDARY CONDITIONS ARE NOT PERMISSIBLE,
C    SO VAL(9) = 0.
      IF (CLASSA.EQ.'lp') THEN
          VAL(7) = 0
          VAL(9) = 0
          NCA = 5
      END IF
      IF (CLASSB.EQ.'lp') THEN
          VAL(8) = 0
          VAL(9) = 0
          NCB = 5
      END IF
      IF ((CLASSA.EQ.'n') .OR. (CLASSA.EQ.'d')) THEN
          VAL(9) = 0
      END IF
      IF ((CLASSB.EQ.'n') .OR. (CLASSB.EQ.'d')) THEN
          VAL(9) = 0
      END IF
c
C  IF ENDPOINT a = -INF, THEN a WAS OMITTED, SO VAL(3) = 0)
C  IF ENDPOINT b = +INF, THEN b WAS OMITTED, SO VAL(4) = 0)
      IF (VAL(3).NE.0 .AND. VAL(4).NE.0) THEN
          INTAB = 1
      ELSE IF (VAL(3).EQ.0 .AND. VAL(4).NE.0) THEN
          INTAB = 3
      ELSE IF (VAL(3).NE.0 .AND. VAL(4).EQ.0) THEN
          INTAB = 2
      ELSE
          INTAB = 4
      END IF
      P0ATA = -1.0D0
      QFATA = 1.0D0
      P0ATB = -1.0D0
      QFATB = 1.0D0
      NCA = 1
      NCB = 1
      NUMBER = NP

C  N.B. FOR classa, or classb, THE VALUE 'd' MEANS DEFAULT.
C    (NOT DIRICHLET)

      IF (VAL(5).NE.0) THEN
          IF (CLASSA.EQ.'r') THEN
              NCA = 1
          ELSE IF (CLASSA.EQ.'wr') THEN
              NCA = 2
          ELSE IF (CLASSA.EQ.'lcno') THEN
              NCA = 3
          ELSE IF (CLASSA.EQ.'lco') THEN
              NCA = 4
          ELSE IF (CLASSA.EQ.'lp') THEN
              NCA = 5
          ELSE IF (CLASSA.EQ.'d') THEN
              NCA = 6
          END IF
      END IF
      IF (VAL(6).NE.0) THEN
          IF (CLASSB.EQ.'r') THEN
              NCB = 1
          ELSE IF (CLASSB.EQ.'wr') THEN
              NCB = 2
          ELSE IF (CLASSB.EQ.'lcno') THEN
              NCB = 3
          ELSE IF (CLASSB.EQ.'lco') THEN
              NCB = 4
          ELSE IF (CLASSB.EQ.'lp') THEN
              NCB = 5
          ELSE IF (CLASSB.EQ.'d') THEN
              NCB = 6
          END IF
      END IF
C
C   WRITE THE RESULTS TO THE OUTPUT FILE:
C
      IF (VAL(11).EQ.0) OPEN (T25,FILE=TAPE25)
C
      IF (VAL(1).NE.0) THEN
          WRITE (T25,FMT=*) ' np = ',NP
      END IF
      IF (VAL(2).EQ.1) THEN
          WRITE (T25,FMT=*) ' param = ',PARAM1
          NU = PARAM1
          KK = PARAM1
          ALPHA = PARAM1
          GAMMA = PARAM1
      ELSE IF (VAL(2).EQ.2) THEN
          WRITE (T25,FMT=*) ' param = ',PARAM1,PARAM2
          NU = PARAM1
          KK = PARAM1
          ALPHA = PARAM1
          GAMMA = PARAM1
          H = PARAM2
          BETA = PARAM2
      ELSE IF (VAL(2).EQ.5) THEN
          WRITE (T25,FMT=*) ' param = ',P1,P2,P3
          WRITE (T25,FMT=*) '         ',P4,P5,P6
      END IF
      IF (VAL(3).NE.0) THEN
          WRITE (T25,FMT=*) ' a = ',A
          IF (CLASSA.NE.'r') THEN
              CALL PQ(-ONE,P0ATA,QFATA)
              WRITE (T25,FMT=*) '   P0ATA,QFATA = ',P0ATA,QFATA
          END IF
      ELSE
          WRITE (T25,FMT=*) ' A = ',' -INF'
      END IF
      IF (VAL(5).NE.0) WRITE (T25,FMT=*) '   CLASSA = ',CLASSA
      IF (VAL(7).EQ.2) WRITE (T25,FMT=*) '   A1,A2 = ',A1,A2
      IF (VAL(4).NE.0) THEN
          WRITE (T25,FMT=*) ' b = ',B
          IF (CLASSB.NE.'r') THEN
              CALL PQ(ONE,P0ATB,QFATB)
              WRITE (T25,FMT=*) '   P0ATB,QFATB = ',P0ATB,QFATB
          END IF
      ELSE
          WRITE (T25,FMT=*) ' B = ','+INF'
      END IF
      IF (VAL(6).NE.0) WRITE (T25,FMT=*) '   CLASSB = ',CLASSB
      IF (VAL(8).EQ.2) WRITE (T25,FMT=*) '   B1,B2 = ',B1,B2
      IF (VAL(9).NE.0) THEN
          IF (VAL(9).EQ.5) WRITE (T25,FMT=*) ' BCC = G '
          WRITE (T25,FMT=*) '    ALFA = ',ALF
          WRITE (T25,FMT=*) ' K11,K12 = ',K11,K12
          WRITE (T25,FMT=*) ' K21,K22 = ',K21,K22
          DET = K11*K22 - K12*K21
          WRITE (T25,FMT=*) ' DET = ',DET
      END IF
      IF (VAL(10).EQ.1) THEN
          WRITE (T25,FMT=*) ' NUMEIG = ',NEIG1
      ELSE IF (VAL(10).EQ.2) THEN
          WRITE (T25,FMT=*) ' NUMEIG1,NUMEIG2 = ',NEIG1,NEIG2
      END IF
C
      WRITE (T25,FMT=*)

C  HERE, VAL(9) = 0 MEANS THE BOUNDARY CONDITIONS ARE SEPARATED.

      IF (VAL(9).EQ.0) THEN
          DO 10 I = NEIG1,NEIG2
              NUMEIG = I
              EIG = 0.0D0
              TOL = ZER
              IFLAG = 1
              AS = A
              BS = B
              INTABS = INTAB
c             CALL SLEIGN(A,B,INTAB,P0ATA,QFATA,P0ATB,QFATB,A1,A2,B1,B2,
              CALL SLEIGN(AS,BS,INTABS,P0ATA,QFATA,P0ATB,QFATB,A1,A2,B1,
     +                    B2,NUMEIG,EIG,TOL,IFLAG,0,SLFUN,NCA,NCB)
c    +                    NUMEIG,EIG,TOL,IFLAG,0,SLFUN,NCA,NCB)
              IF (IFLAG.EQ.15) THEN
                  WRITE (T25,FMT=*)
     +              ' WE CANNOT HANDLE THIS KIND OF ENDPOINT. '
                  WRITE (*,FMT=*)
                  GO TO 12
              END IF
              IFLAG = MIN(IFLAG,4)
              JFLAG = 0
              IF (SLFUN(9).GT.-10000.0D0) JFLAG = 1
              IF (SLFUN(9).LT.10000.0D0 .AND. JFLAG.EQ.1) JFLAG = 2
              IF (IFLAG.EQ.4) WRITE (*,FMT=9045) NUMEIG
              IF (IFLAG.EQ.4) WRITE (T25,FMT=9045) NUMEIG
              IF (IFLAG.LE.3) WRITE (*,FMT=*) ' IFLAG = ',IFLAG
              IF (IFLAG.LE.3) WRITE (T25,FMT=*) ' IFLAG = ',IFLAG
              IF (IFLAG.EQ.3) THEN
                IF (NUMEIG.GE.0) THEN
                  WRITE (T25,FMT=9230) NUMEIG
                ELSE IF( .NOT.(NCA.EQ.4 .OR. NCB.EQ.4)) THEN
                  WRITE (T25,FMT=9235)
                ELSE
                  WRITE (T25,FMT=9240)
                ENDIF
              ELSE IF (IFLAG.LE.2) THEN
                  WRITE (*,FMT=*) ' NUMEIG,EIG,TOL = ',NUMEIG,EIG,TOL
                  WRITE (T25,FMT=*) ' NUMEIG,EIG,TOL = ',NUMEIG,EIG,TOL
              END IF
   10     CONTINUE
   12     CONTINUE
          WRITE (*,FMT=*)
          WRITE (T25,FMT=*)
          IF (JFLAG.EQ.1) WRITE (*,FMT=9260)
          IF (JFLAG.EQ.2) WRITE (*,FMT=9250) SLFUN(9)
          IF (JFLAG.EQ.1) WRITE (T25,FMT=9260)
          IF (JFLAG.EQ.2) WRITE (T25,FMT=9250) SLFUN(9)
      ELSE
          DO 20 I = NEIG1,NEIG2
              NUMEIG = I
              EIG = 0.0D0
              TOL = ZER
              AS = A
              BS = B
              INTABS = INTAB
              IFLAG = 1
              CALL PERIO(AS,BS,INTABS,P0ATA,QFATA,P0ATB,QFATB,A1,A2,B1,
     +                   B2,NUMEIG,EIG,TOL,IFLAG,SLFUN,NCA,NCB,ALF,K11,
     +                   K12,K21,K22)
              IFLAG = MIN(IFLAG,4)
              IF (IFLAG.EQ.0) IFLAG = 4
              IF (IFLAG.EQ.4) WRITE (*,FMT=9045) NUMEIG
              IF (IFLAG.EQ.4) WRITE (T25,FMT=9045) NUMEIG
              IF (IFLAG.LE.2) THEN
                  WRITE (*,FMT=*) ' IFLAG = ',IFLAG
                  WRITE (T25,FMT=*) ' IFLAG = ',IFLAG
                  WRITE (*,FMT=*) ' NUMEIG,EIG,TOL = ',NUMEIG,EIG,TOL
                  WRITE (T25,FMT=*) ' NUMEIG,EIG,TOL = ',NUMEIG,EIG,TOL
              END IF
              IF (IFLAG.EQ.2) WRITE (T25,FMT=*
     +            ) ' THIS EIGENVALUE APPEARS TO BE DOUBLE. '
   20     CONTINUE
      END IF
      WRITE (T25,FMT=*) '______________________________________________'
C
C  HERE, VAL(12) .NE. 0 MEANS THAT ANOTHER PROBLEM IS TO BE RUN.
C  SO ANY CHANGES IN THE ENDPOINTS, OR BOUNDARY CONDITIONS, OR
C  ETC., MUST BE READ IN.  ANY PARAMETERS THAT ARE NOT CHANGED
C  ARE KEPT AS BEFORE.

      IF (VAL(12).NE.0) THEN
          REWIND T24
          DO 80 I = 1,50
              CH(I) = BLANK
   80     CONTINUE
          LAST = LAST + LASTI
          READ (T24,FMT=9010) (CHT,I=1,LAST)
          READ (T24,FMT=9010,END=85) (CH(I),I=1,50)
   85     CONTINUE
          VAL(12) = 0
          GO TO 300
      END IF
C
      CLOSE (T21)
      CLOSE (T24)
      WRITE (T25,FMT=*) 'end ',TAPE25
      CLOSE (T25)
C
      STOP
 9010 FORMAT (A62)
 9045 FORMAT (1X,'   NUMEIG = ',I5,' : COMPUTATION FAILED ')
 9230 FORMAT (1X,' * THERE SEEMS TO BE NO EIGENVALUE OF INDEX       *',
     +       /,1X,' * GREATER THAN',I5,'                              *'
     +       )
 9235 FORMAT (1X,' * THERE SEEM TO BE NO EIGENVALUES                *')
 9240 FORMAT (1X,' * THERE SEEMS TO BE NO EIGENVALUE OF THIS INDEX  *')
 9250 FORMAT (1X,' * THERE SEEMS TO BE CONTINUOUS SPECTRUM BEGINNING*',
     +       /,1X,' * AT ABOUT',1P,D8.1,
     +       '                               *')
 9260 FORMAT (1X,' * THERE SEEMS TO BE NO CONTINUOUS SPECTRUM       *')
      END
c
      SUBROUTINE STR2R(CHIN,P)
C  THIS PROGRAM EXPECTS A CHARACTER STRING, CHIN, CONTAINING
C    DIGITS AND, POSSIBLY, COMMAS.  IT CONVERTS THE FIRST SUBSTRING
C    TO BE THE CORRESPONDING REAL NUMBER P.
C    THUS THE STRING
C        4.0,2.0,1.5,2.5
C    WOULD CAUSE P TO BE SET EQUAL TO THE REAL NUMBER 4.0,
C    AND CHIN ON OUTPUT WOULD BE THE STRING
C        2.0,1.5,2.5
C     .. Scalar Arguments ..
      REAL P
      CHARACTER*62 CHIN
C     ..
C     .. Local Scalars ..
      INTEGER I,J,K
      CHARACTER CH1
      CHARACTER*62 CHINN,CHOUT
C     ..
      CH1 = ','
      CHINN = '                               '
      CHOUT = '                               '
      DO 10 I = 1,60
          IF ((CHIN(I:I).EQ.CH1) .OR. (I.GT.1.AND.CHIN(I:I).EQ.' '))
     +        THEN
              J = I - 1
              GO TO 20
          END IF
   10 CONTINUE
   20 CONTINUE
      DO 30 I = 1,J
          CHINN(I:I) = CHIN(I:I)
   30 CONTINUE
      READ (CHINN,FMT=9020) P
      DO 40 I = J + 2,62
          K = I - J - 1
          CHOUT(K:K) = CHIN(I:I)
   40 CONTINUE
      CHIN = CHOUT
 9020 FORMAT (5F12.7)
      RETURN
      END
c
      SUBROUTINE CHAR(CHIN,K,M,CHOUT)
C
C  THIS PROGRAM IS INTENDED TO READ A LINE LIKE
C     '  abcd:  0.532,12.34,-57.000,0.7693 '
C  AND RETURN WITH
C    K = 9     (the number of characters up to the first non-blank
C                one after the : )
C    M = 3     (the number of commas after the : )
C  AND WITH
C    CHOUT = '0.532,12.34,-57.000,0.7693'
C  THUS READING THE LINE
C    'param:  1.23,-0.43,22.7,0.0037,-11.21'
C  THE RESULT WOULD BE 
C     K = 8
C     M = 4
C     CHOUT = 1.23,-0.43,22.7,0.0037,-11.21
C
C     .. Scalar Arguments ..
      INTEGER K,M
      CHARACTER*62 CHIN,CHOUT
C     ..
C     .. Local Scalars ..
      INTEGER I,J,L
      CHARACTER*16 CH2
C     ..
C     .. Local Arrays ..
      CHARACTER CH(62)
C     ..
      M = 0
      DO 10 I = 1,8
          READ (CHIN,FMT=100) (CH(J),J=1,I)
          IF (CH(I).EQ.':') GO TO 20
          M = I
   10 CONTINUE
   20 CONTINUE
      M = 0
      J = 0
      K = 0
      DO 30 I = 1,12
          READ (CHIN,FMT=100) (CH(L),L=1,I)
          IF (J.NE.0 .AND. K.EQ.0 .AND. CH(I).NE.' ') K = I
          IF (M.NE.0 .AND. K.EQ.0 .AND. CH(I).EQ.' ') J = I
          IF (CH(I).EQ.':') M = I
   30 CONTINUE
      M = 0
      DO 45 I = 1,62
          READ (CHIN,FMT=100) (CH(L),L=1,I)
          IF (CH(I).EQ.',') M = M + 1
   45 CONTINUE
      K = J
      IF (K.EQ.1) THEN
          READ (CHIN,FMT=101) CH2,CHOUT
      ELSE IF (K.EQ.2) THEN
          READ (CHIN,FMT=102) CH2,CHOUT
      ELSE IF (K.EQ.3) THEN
          READ (CHIN,FMT=103) CH2,CHOUT
      ELSE IF (K.EQ.4) THEN
          READ (CHIN,FMT=104) CH2,CHOUT
      ELSE IF (K.EQ.5) THEN
          READ (CHIN,FMT=105) CH2,CHOUT
      ELSE IF (K.EQ.6) THEN
          READ (CHIN,FMT=106) CH2,CHOUT
      ELSE IF (K.EQ.7) THEN
          READ (CHIN,FMT=107) CH2,CHOUT
      ELSE IF (K.EQ.8) THEN
          READ (CHIN,FMT=108) CH2,CHOUT
      ELSE IF (K.EQ.9) THEN
          READ (CHIN,FMT=109) CH2,CHOUT
      ELSE IF (K.EQ.10) THEN
          READ (CHIN,FMT=110) CH2,CHOUT
      ELSE IF (K.EQ.11) THEN
          READ (CHIN,FMT=111) CH2,CHOUT
      ELSE IF (K.EQ.12) THEN
          READ (CHIN,FMT=112) CH2,CHOUT
      END IF
      RETURN
  100 FORMAT (62A1)
  101 FORMAT (A1,A61)
  102 FORMAT (A2,A60)
  103 FORMAT (A3,A59)
  104 FORMAT (A4,A58)
  105 FORMAT (A5,A57)
  106 FORMAT (A6,A56)
  107 FORMAT (A7,A55)
  108 FORMAT (A8,A54)
  109 FORMAT (A9,A53)
  110 FORMAT (A10,A52)
  111 FORMAT (A11,A51)
  112 FORMAT (A12,A50)
      END
C
      SUBROUTINE LST(CHANS,A)
C  THIS PROGRAM CONVERTS A CHARACTER STRING CHANS TO A REAL
C    NUMBER A.  IT JUST ALLOWS THE NUMBERS PI, PI/2, PI/4,
C    2PI TO BE READ IN AS CHARACTERS INSTEAD OF HAVING TO
C    ENTER THEM AS DECIMAL DIGITS.
C
C     .. Scalar Arguments ..
      REAL A
      CHARACTER*32 CHANS
C     ..
C     .. Local Scalars ..
      REAL ONE,PI,PI4
      INTEGER I
C     ..
C     .. Local Arrays ..
      REAL Y(48)
      CHARACTER*8 X(48)
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ATAN
C     ..
      ONE = 1.0D0
      PI4 = ATAN(ONE)
      PI = 4.0D0*PI4
c
      X(1) = '-PI'
      Y(1) = -PI
      X(2) = '-2.0*PI'
      Y(2) = -2.0D0*PI
      X(3) = '-2.*PI'
      Y(3) = -2.0D0*PI
      X(4) = '-.5*PI'
      Y(4) = -0.5D0*PI
      X(5) = '-0.5*PI'
      Y(5) = -0.5D0*PI
      X(6) = '-.25*PI'
      Y(6) = -0.25D0*PI
      X(7) = '-0.25*PI'
      Y(7) = -0.25D0*PI
      X(8) = '-PI/2'
      Y(8) = -0.5D0*PI
      X(9) = '-PI/4'
      Y(9) = -0.25D0*PI
c
      X(10) = '-pi'
      Y(10) = -PI
      X(11) = '-2.0*pi'
      Y(11) = -2.0D0*PI
      X(12) = '-2.*pi'
      Y(12) = -2.0D0*PI
      X(13) = '-.5*pi'
      Y(13) = -0.5D0*PI
      X(14) = '-0.5*pi'
      Y(14) = -0.5D0*PI
      X(15) = '-.25*pi'
      Y(15) = -0.25D0*PI
      X(16) = '-0.25*pi'
      Y(16) = -0.25D0*PI
      X(17) = '-pi/2'
      Y(17) = -0.5D0*PI
      X(18) = '-pi/4'
      Y(18) = -0.25D0*PI
c
      X(19) = 'PI'
      Y(19) = PI
      X(20) = '2.0*PI'
      Y(20) = 2.0D0*PI
      X(21) = '2.*PI'
      Y(21) = 2.0D0*PI
      X(22) = '.5*PI'
      Y(22) = 0.5D0*PI
      X(23) = '.25*PI'
      Y(23) = 0.25D0*PI
      X(24) = 'PI/2'
      Y(24) = 0.5D0*PI
      X(25) = 'PI/4'
      Y(25) = 0.25D0*PI
c
      X(26) = 'pi'
      Y(26) = PI
      X(27) = '2.0*pi'
      Y(27) = 2.0D0*PI
      X(28) = '2.*pi'
      Y(28) = 2.0D0*PI
      X(29) = '.5*pi'
      Y(29) = 0.5D0*PI
      X(30) = '.25*pi'
      Y(30) = 0.25D0*PI
      X(31) = 'pi/2'
      Y(31) = 0.5D0*PI
      X(32) = 'pi/4'
      Y(32) = 0.25D0*PI
c
      X(33) = '0.5*PI'
      Y(33) = 0.5D0*PI
      X(34) = '0.5*pi'
      Y(34) = 0.5D0*PI
      X(35) = '-0.5*PI'
      Y(35) = -0.5D0*PI
      X(36) = '-0.5*pi'
      Y(36) = -0.5D0*PI
      X(37) = 'PI/2.0'
      Y(37) = PI/2.0D0
      X(38) = 'pi/2.0'
      Y(38) = PI/2.0D0
      X(39) = '-PI/2.0'
      Y(39) = -PI/2.0D0
      X(40) = '-pi/2.0'
      Y(40) = -PI/2.0D0
      X(41) = 'PI/4.0'
      Y(41) = PI/4.0D0
      X(42) = 'pi/4.0'
      Y(42) = PI/4.0D0
      X(43) = '-PI/4.0'
      Y(43) = -PI/4.0D0
      X(44) = '-pi/4.0'
      Y(44) = -PI/4.0D0
      X(45) = '-0.25*PI'
      Y(45) = -0.25D0*PI
      X(46) = '-0.25*pi'
      Y(46) = -0.25D0*PI
      X(47) = '0.25*PI'
      Y(47) = 0.25D0*PI
      X(48) = '0.25*pi'
      Y(48) = 0.25D0*PI
C
      DO 10 I = 1,48
          IF (CHANS.EQ.X(I)) THEN
              A = Y(I)
              RETURN
          END IF
   10 CONTINUE
      READ (CHANS,FMT='(F32.0)') A
C
      RETURN
      END
C
      SUBROUTINE HELP(NH)
C
C     .. Scalar Arguments ..
      INTEGER NH
C     ..
C     .. Local Scalars ..
      INTEGER I,N
      CHARACTER ANS
C     ..
C     .. Local Arrays ..
      CHARACTER*36 X(23),Y(23)
C     ..
      GO TO (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17) NH
c
    1 CONTINUE
      WRITE (*,FMT=*)
      WRITE (*,FMT=*) 'H1:  Overview of HELP.'
      X(1) = '  This ASCII text file is supplied a'
      Y(1) = 's a separate file with the sleign2  '
      X(2) = 'package; it can be accessed on-line '
      Y(2) = 'in both MAKEPQW (if used) and DRIVE.'
      X(3) = '  HELP contains information to aid t'
      Y(3) = 'he user in entering data on the co- '
      X(4) = 'efficient functions p,q,w; on the se'
      Y(4) = 'lf-adjoint, separated and coupled,  '
      X(5) = 'regular and singular, boundary condi'
      Y(5) = 'tions; on the limit circle boundary '
      X(6) = 'condition functions u,v at a and U,V'
      Y(6) = ' at b; on the end-point classifica- '
      X(7) = 'tions of the differential equation; '
      Y(7) = 'on DEFAULT entry; on eigenvalue in- '
      X(8) = 'dexes; on IFLAG information; and on '
      Y(8) = 'the general use of the program      '
      X(9) = 'sleign2.                            '
      Y(9) = ' '
      X(10) = '  The 17 sections of HELP are:      '
      Y(10) = ' '
      X(11) = ' '
      Y(11) = ' '
      X(12) = '    H1: Overview of HELP.           '
      Y(12) = ' '
      X(13) = '    H2: File name entry.            '
      Y(13) = ' '
      X(14) = '    H3: The differential equation.  '
      Y(14) = ' '
      X(15) = '    H4: End-point classification.   '
      Y(15) = ' '
      X(16) = '    H5: DEFAULT entry.              '
      Y(16) = ' '
      X(17) = '    H6: Self-adjoint limit-circle bo'
      Y(17) = 'undary conditions.                  '
      DO 101 I = 1,17
          WRITE (*,FMT=*) X(I),Y(I)
  101 CONTINUE
      WRITE (*,FMT=*)
      WRITE (*,FMT=*)
      WRITE (*,FMT=*)
      READ (*,FMT=999) ANS,N
      IF (ANS.EQ.'r' .OR. ANS.EQ.'R') RETURN
      GO TO (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17) N
      X(1) = '    H7: General self-adjoint boundar'
      Y(1) = 'y conditions. '
      X(2) = '    H8: Recording the results.      '
      Y(2) = ' '
      X(3) = '    H9: Type and choice of interval.'
      Y(3) = ' '
      X(4) = '   H10: Entry of end-points.        '
      Y(4) = ' '
      X(5) = '   H11: End-point values of p,q,w.  '
      Y(5) = ' '
      X(6) = '   H12: Initial value problems.     '
      Y(6) = ' '
      X(7) = '   H13: Indexing of eigenvalues for '
      Y(7) = 'separated boundary conditions.      '
      X(8) = '   H14: Entry of eigenvalue index, i'
      Y(8) = 'nitial guess, and tolerance.     '
      X(9) = '   H15: IFLAG information.          '
      Y(9) = ' '
      X(10) = '   H16: Plotting.                  '
      Y(10) = ' '
      X(11) = '   H17: Indexing of eigenvalues for'
      Y(11) = 'coupled boundary conditions.       '
      X(12) = ' '
      Y(12) = ' '
      X(13) = '  HELP can be accessed at each point'
      Y(13) = ' in MAKEPQW and DRIVE where the user'
      X(14) = 'is asked for input, by pressing "h <'
      Y(14) = 'ENTER>"; this places the user at the'
      X(15) = 'appropriate HELP section.  Once in H'
      Y(15) = 'ELP, the user can scroll the further'
      X(16) = 'HELP sections by repeatedly pressing'
      Y(16) = ' "h <ENTER>", or jump to a specific '
      X(17) = 'HELP section Hn (n=1,2,...17) by typ'
      Y(17) = 'ing "Hn <ENTER>"; to return to the  '
      X(18) = 'place in the program from which HELP'
      Y(18) = ' is called, press "r <ENTER>".      '
      DO 102 I = 1,18
          WRITE (*,FMT=*) X(I),Y(I)
  102 CONTINUE
      WRITE (*,FMT=*)
      WRITE (*,FMT=*)
      WRITE (*,FMT=*)
      WRITE (*,FMT=*)
      READ (*,FMT=999) ANS,N
      IF (ANS.EQ.'r' .OR. ANS.EQ.'R') RETURN
      GO TO (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17) N
    2 CONTINUE
      WRITE (*,FMT=*)
      WRITE (*,FMT=*) 'H2:  File name entry.'
      X(1) = '  MAKEPQW is used to create a FORTRA'
      Y(1) = 'N file containing the coefficients  '
      X(2) = 'p(x),q(x),w(x), defining the differe'
      Y(2) = 'ntial equation, and the boundary    '
      X(3) = 'condition functions u(x),v(x) and U('
      Y(3) = 'x),V(x) if required.  The file must '
      X(4) = 'be given a NEW filename which is acc'
      Y(4) = 'eptable to your FORTRAN compiler.   '
      X(5) = 'For example, it might be called bess'
      Y(5) = 'el.f or bessel.for, depending upon  '
      X(6) = 'your compiler.                      '
      Y(6) = ' '
      X(7) = '  The same naming considerations app'
      Y(7) = 'ly if the FORTRAN file is prepared  '
      X(8) = 'other than with the use of MAKEPQW. '
      Y(8) = ' '
      DO 201 I = 1,8
          WRITE (*,FMT=*) X(I),Y(I)
  201 CONTINUE
      WRITE (*,FMT=*)
      WRITE (*,FMT=*)
      WRITE (*,FMT=*)
      WRITE (*,FMT=*)
      WRITE (*,FMT=*)
      WRITE (*,FMT=*)
      WRITE (*,FMT=*)
      WRITE (*,FMT=*)
      WRITE (*,FMT=*)
      READ (*,FMT=999) ANS,N
      IF (ANS.EQ.'r' .OR. ANS.EQ.'R') RETURN
      GO TO (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17) N
    3 CONTINUE
      WRITE (*,FMT=*) 'H3:  The differential equation.'
      X(1) = '  The prompt "Input p (or q or w) ="'
      Y(1) = ' requests you to type in a FORTRAN  '
      X(2) = 'expression defining the function p(x'
      Y(2) = '), which is one of the three coeffi-'
      X(3) = 'cient functions defining the Sturm-L'
      Y(3) = 'iouville differential equation      '
      X(4) = ' '
      Y(4) = ' '
      X(5) = '                   -(p*y'')'' + q*y = '
      Y(5) = ' lambda*w*y              (*)        '
      X(6) = ' '
      Y(6) = ' '
      X(7) = 'to be considered on some interval (a'
      Y(7) = ',b) of the real line.  The actual   '
      X(8) = 'interval used in a particular proble'
      Y(8) = 'm can be chosen later, and may be   '
      X(9) = 'either the whole interval (a,b) wher'
      Y(9) = 'e the coefficient functions p,q,w,  '
      X(10) = 'etc. are defined or any sub-interval'
      Y(10) = ' (a'',b'') of (a,b); a = -infinity  '
      X(11) = 'and/or b = +infinity are allowable c'
      Y(11) = 'hoices for the end-points.          '
      X(12) = '  The coefficient functions p,q,w of'
      Y(12) = ' the differential equation may be   '
      X(13) = 'chosen arbitrarily but must satisfy '
      Y(13) = 'the following conditions:           '
      X(14) = '  1. p,q,w are real-valued throughou'
      Y(14) = 't (a,b).                            '
      X(15) = '  2. p,q,w are piece-wise continuous'
      Y(15) = ' and defined throughout the         '
      X(16) = '     interior of the interval (a,b).'
      Y(16) = '                                    '
      X(17) = '  3. p and w are strictly positive i'
      Y(17) = 'n (a,b).                            '
      X(18) = '     For better error analysis in th'
      Y(18) = 'e numerical procedures, condition 2.'
      X(19) = '     above is often replaced with   '
      Y(19) = ' '
      X(20) = '  2''. p,q,w are four times continuou'
      Y(20) = 'sly differentiable on (a,b).        '
      X(21) = '  The behavior of p,q,w near the end'
      Y(21) = '-points a and b is critical to the  '
      X(22) = 'classification of the differential e'
      Y(22) = 'quation (see H4 and H11).           '
      DO 301 I = 1,22
          WRITE (*,FMT=*) X(I),Y(I)
  301 CONTINUE
      WRITE (*,FMT=*)
      READ (*,FMT=999) ANS,N
      IF (ANS.EQ.'r' .OR. ANS.EQ.'R') RETURN
      GO TO (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17) N
    4 CONTINUE
      WRITE (*,FMT=*) 'H4:  End-point classification.'
      X(1) = '  The correct classification of the '
      Y(1) = 'end-points a and b is essential to  '
      X(2) = 'the working of the sleign2 program. '
      Y(2) = ' To classify the end-points, it is  '
      X(3) = 'convenient to choose a point c in (a'
      Y(3) = ',b); i.e.,  a < c < b.  Subject to  '
      X(4) = 'the general conditions on the coeffi'
      Y(4) = 'cient functions p,q,w (see H3):     '
      X(5) = '  1. a is REGULAR (say R) if -infini'
      Y(5) = 'ty < a, p,q,w are piece-wise        '
      X(6) = '     continuous on [a,c], and p(a) >'
      Y(6) = ' 0 and w(a) > 0.                    '
      X(7) = '  2. a is WEAKLY REGULAR (say WR) if'
      Y(7) = ' -infinity < a, a is not R, and     '
      X(8) = '                     |c             '
      Y(8) = ' '
      X(9) = '            integral | {1/p+|q|+w} <'
      Y(9) = ' +infinity.                         '
      X(10) = '                     |a             '
      Y(10) = ' '
      X(11) = '     If end-point a is neither R nor'
      Y(11) = ' WR, then a is SINGULAR; that is,   '
      X(12) = '     either -infinity = a, or -infin'
      Y(12) = 'ity < a and                         '
      X(13) = '                     |c             '
      Y(13) = ' '
      X(14) = '            integral | {1/p+|q|+w} ='
      Y(14) = ' +infinity.                         '
      X(15) = '                     |a             '
      Y(15) = ' '
      X(16) = '  3. The SINGULAR end-point a is LIM'
      Y(16) = 'IT-CIRCLE NON-OSCILLATORY (say      '
      X(17) = '     LCNO) if for some real lambda A'
      Y(17) = 'LL real-valued solutions y of the   '
      X(18) = '     differential equation          '
      Y(18) = ' '
      X(19) = ' '
      Y(19) = ' '
      X(20) = '                   -(p*y'')'' + q*y = '
      Y(20) = 'lambda*w*y on (a,c]     (*)         '
      X(21) = ' '
      Y(21) = ' '
      X(22) = '     satisfy the conditions:        '
      Y(22) = ' '
      DO 401 I = 1,22
          WRITE (*,FMT=*) X(I),Y(I)
  401 CONTINUE
      WRITE (*,FMT=*)
      READ (*,FMT=999) ANS,N
      IF (ANS.EQ.'r' .OR. ANS.EQ.'R') RETURN
      GO TO (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17) N
      X(1) = '                     |c             '
      Y(1) = ' '
      X(2) = '            integral | { w*y*y } < +'
      Y(2) = 'infinity, and                       '
      X(3) = '                     |a             '
      Y(3) = ' '
      X(4) = '     y has at most a finite number o'
      Y(4) = 'f zeros in (a,c].                   '
      X(5) = '  4. The SINGULAR end-point a is LIM'
      Y(5) = 'IT-CIRCLE OSCILLATORY (say LCO) if  '
      X(6) = '     for some real lambda ALL real-v'
      Y(6) = ' alued solutions of the differential'
      X(7) = '     equation (*) satisfy the condit'
      Y(7) = 'ions:                               '
      X(8) = '                     |c             '
      Y(8) = ' '
      X(9) = '            integral | { w*y*y } < +'
      Y(9) = 'infinity, and                       '
      X(10) = '                     |a             '
      Y(10) = ' '
      X(11) = '     y has an infinite number of zer'
      Y(11) = 'os in (a,c].                        '
      X(12) = '  5. The SINGULAR end-point a is LIM'
      Y(12) = 'IT POINT (say LP) if for some real  '
      X(13) = '     lambda at least one solution of'
      Y(13) = ' the differential equation (*) sat- '
      X(14) = '     isfies the condition:          '
      Y(14) = ' '
      X(15) = '                     |c             '
      Y(15) = ' '
      X(16) = '            integral | {w*y*y} = +in'
      Y(16) = 'finity.                             '
      X(17) = '                     |a             '
      Y(17) = ' '
      X(18) = '  There is a similar classification '
      Y(18) = 'of the end-point b into one of the  '
      X(19) = 'five distinct cases R, WR, LCNO, LCO'
      Y(19) = ', LP.                               '
      DO 402 I = 1,19
          WRITE (*,FMT=*) X(I),Y(I)
  402 CONTINUE
      WRITE (*,FMT=*)
      WRITE (*,FMT=*)
      WRITE (*,FMT=*)
      WRITE (*,FMT=*)
      READ (*,FMT=999) ANS,N
      IF (ANS.EQ.'r' .OR. ANS.EQ.'R') RETURN
      GO TO (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17) N
      X(1) = '  Although the classification of sin'
      Y(1) = 'gular end-points invokes a real     '
      X(2) = 'value of the parameter lambda, this '
      Y(2) = 'classification is invariant in      '
      X(3) = 'lambda; all real choices of lambda l'
      Y(3) = 'ead to the same classification.     '
      X(4) = '  In determining the classification '
      Y(4) = 'of singular end-points for the      '
      X(5) = 'differential equation (*), it is oft'
      Y(5) = 'en convenient to start with the     '
      X(6) = 'choice lambda = 0 in attempting to f'
      Y(6) = 'ind solutions (particularly when    '
      X(7) = 'q = 0 on (a,b)); however, see exampl'
      Y(7) = 'e 7 below.                          '
      X(8) = '  See H6 on the use of maximal domai'
      Y(8) = 'n functions to determine the        '
      X(9) = 'classification at singular end-point'
      Y(9) = 's.                                  '
      DO 403 I = 1,9
          WRITE (*,FMT=*) X(I),Y(I)
  403 CONTINUE
      WRITE (*,FMT=*)
      WRITE (*,FMT=*)
      WRITE (*,FMT=*)
      WRITE (*,FMT=*)
      WRITE (*,FMT=*)
      WRITE (*,FMT=*)
      WRITE (*,FMT=*)
      WRITE (*,FMT=*)
      WRITE (*,FMT=*)
      WRITE (*,FMT=*)
      READ (*,FMT=999) ANS,N
      IF (ANS.EQ.'r' .OR. ANS.EQ.'R') RETURN
      GO TO (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17) N
      WRITE (*,FMT=*) ' EXAMPLES: '
      X(1) = '  1. -y'''' = lambda*y is R at both en'
      Y(1) = 'd-points of (a,b) when a and b are  '
      X(2) = '     finite.                        '
      Y(2) = ' '
      X(3) = '  2. -y'''' = lambda*y on (-infinity,i'
      Y(3) = 'nfinity) is LP at both end-points.  '
      X(4) = '  3. -(sqrt(x)*y''(x))'' = lambda*(1./'
      Y(4) = 'sqrt(x))*y(x) on (0,infinity) is    '
      X(5) = '     WR at 0 and LP at +infinity (ta'
      Y(5) = 'ke lambda = 0 in (*)).  See         '
      X(6) = '     examples.f, #10 (Weakly Regular'
      Y(6) = ').                                  '
      X(7) = '  4. -((1-x*x)*y''(x))'' = lambda*y(x)'
      Y(7) = ' on (-1,1) is LCNO at both ends     '
      X(8) = '     (take lambda = 0 in (*)).  See '
      Y(8) = 'xamples.f, #1 (Legendre).           '
      X(9) = '  5. -y''''(x) + C*(1/(x*x))*y(x) = la'
      Y(9) = 'mbda*y(x) on (0,infinity) is LP at  '
      X(10) = '     infinity and at 0 is (take lamb'
      Y(10) = 'da = 0 in (*)):                     '
      X(11) = '       LP for C .ge. 3/4 ;          '
      Y(11) = ' '
      X(12) = '       LCNO for -1/4 .le. C .lt. 3/4'
      Y(12) = ' (but C .ne. 0);                    '
      X(13) = '       LCO for C .lt. -1/4.         '
      Y(13) = ' '
      X(14) = '  6. -(x*y''(x))'' - (1/x)*y(x) = lamb'
      Y(14) = 'da*y(x) on (0,infinity) is LCO at 0 '
      X(15) = '     and LP at +infinity (take lambd'
      Y(15) = 'a = 0 in (*) with solutions         '
      X(16) = '     cos(ln(x)) and sin(ln(x))).  Se'
      Y(16) = 'e xamples.f, #7 (BEZ).              '
      X(17) = '  7. -(x*y''(x))'' - x*y(x) = lambda*('
      Y(17) = '1/x)*y(x) on (0,infinity) is LP at 0'
      X(18) = '     and LCO at infinity (take lambd'
      Y(18) = 'a = -1/4 in (*) with solutions      '
      X(19) = '     cos(x)/sqrt(x) and sin(x)/sqrt('
      Y(19) = 'x)).  See xamples.f,                '
      X(20) = '     #6 (Sears-Titchmarsh).         '
      Y(20) = ' '
      X(21) = '  8. -y''''(x) + x*sin(x)*y(x) = lambd'
      Y(21) = 'a*y(x) on (0,infinity) is R at 0 and'
      X(22) = '     LP at infinity.  See examples.f'
      Y(22) = ' #30 (Littlewood-McLeod).           '
      DO 404 I = 1,22
          WRITE (*,FMT=*) X(I),Y(I)
  404 CONTINUE
      WRITE (*,FMT=*) ' '
      READ (*,FMT=999) ANS,N
      IF (ANS.EQ.'r' .OR. ANS.EQ.'R') RETURN
      GO TO (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17) N
    5 CONTINUE
      WRITE (*,FMT=*) 'H5:  DEFAULT entry.'
      X(1) = '  The complete range of problems for'
      Y(1) = ' which sleign2 is applicable can    '
      X(2) = 'only be reached by appropriate entri'
      Y(2) = 'es under end-point classification   '
      X(3) = 'and boundary conditions.  However, t'
      Y(3) = 'here is a DEFAULT application which '
      X(4) = 'requires no detailed entry of end-po'
      Y(4) = 'int classification or boundary      '
      X(5) = 'conditions, subject to:             '
      Y(5) = ' '
      X(6) = '  1. The DEFAULT application CANNOT '
      Y(6) = 'be used at a LCO end-point.         '
      X(7) = '  2. If an end-point a is R, then th'
      Y(7) = 'e Dirichlet boundary condition      '
      X(8) = '     y(a) = 0 is automatically used.'
      Y(8) = ' '
      X(9) = '  3. If an end-point a is WR, then t'
      Y(9) = 'he following boundary condition     '
      X(10) = '     is automatically applied:      '
      Y(10) = ' '
      X(11) = '       if p(a) = 0, and both q(a),w('
      Y(11) = 'a) are bounded, then the Neumann    '
      X(12) = '       boundary condition (py'')(a) '
      Y(12) = '= 0 is used, or                     '
      X(13) = '       if p(a) > 0, and q(a) and/or '
      Y(13) = 'w(a)) are not bounded, then the     '
      X(14) = '       Dirichlet boundary condition '
      Y(14) = 'y(a) = 0 is used.                   '
      X(15) = '       If p(a) = 0, and q(a) and/or '
      Y(15) = 'w(a) are not bounded, then no simple'
      X(16) = '       information, in general, can '
      Y(16) = ' be given on the DEFAULT boundary   '
      X(17) = '       condition.                   '
      Y(17) = ' '
      X(18) = '  4. If an end-point is LCNO, then i'
      Y(18) = 'n most cases the principal or       '
      X(19) = '     Friedrichs boundary condition i'
      Y(19) = 's applied (see H6).                 '
      X(20) = '  5. If an end-point is LP, then the'
      Y(20) = ' normal LP procedure is applied     '
      X(21) = '     (see H7(1.)).                  '
      Y(21) = ' '
      X(22) = 'If the DEFAULT condition is chosen, '
      Y(22) = 'then no entry is required for the   '
      X(23) = 'u,v and U,V boundary condition funct'
      Y(23) = 'ions. '
      DO 501 I = 1,23
          WRITE (*,FMT=*) X(I),Y(I)
  501 CONTINUE
      READ (*,FMT=999) ANS,N
      IF (ANS.EQ.'r' .OR. ANS.EQ.'R') RETURN
      GO TO (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17) N
    6 CONTINUE
      WRITE (*,FMT=*) 'H6:  Limit-circle boundary conditions.'
      X(1) = '  At an end-point a, the limit-circl'
      Y(1) = 'e type separated boundary condition '
      X(2) = 'is of the form (similar remarks thro'
      Y(2) = 'ughout apply to the end-point b with'
      X(3) = ' U,V being boundary condition functi'
      Y(3) = 'ons at b)                           '
      X(4) = ' '
      Y(4) = ' '
      X(5) = '       A1*[y,u](a) + A2*[y,v](a) = 0'
      Y(5) = ',                          (**)     '
      X(6) = ' '
      Y(6) = ' '
      X(7) = 'where y is a solution of the differe'
      Y(7) = 'ntial equation                      '
      X(8) = ' '
      Y(8) = ' '
      X(9) = '     -(p*y'')'' + q*y = lambda*w*y  on'
      Y(9) = ' (a,b).                     (*)     '
      X(10) = ' '
      Y(10) = ' '
      X(11) = 'Here A1, A2 are real numbers; u and '
      Y(11) = 'v are boundary condition functions  '
      X(12) = 'at a; and for real-valued y and u th'
      Y(12) = 'e form [y,u] is defined by          '
      X(13) = '                                    '
      Y(13) = '                                    '
      X(14) = '   [y,u](x) = y(x)*(pu'')(x) - u(x)*'
      Y(14) = '(py'')(x)   for x in (a,b).         '
      X(15) = '                                    '
      Y(15) = '                                    '
      X(16) = '  If neither end-point is LP then th'
      Y(16) = 'ere are also self-adjoint coupled   '
      X(17) = 'boundary conditions.  These have a c'
      Y(17) = 'anonical form given by              '
      X(18) = '                                    '
      Y(18) = '                                    '
      X(19) = '    Y(b) = exp(i*alpha)*K*Y(a),     '
      Y(19) = ' '
      X(20) = '                                    '
      Y(20) = '                                    '
      X(21) = 'where K is a real 2 by 2 matrix with'
      Y(21) = ' determinant 1, alpha is restricted '
      X(22) = 'to the interval (-pi,pi], and Y is  '
      Y(22) = '  '
      DO 551 I = 1,22
          WRITE (*,FMT=*) X(I),Y(I)
  551 CONTINUE
      READ (*,FMT=999) ANS,N
      IF (ANS.EQ.'r' .OR. ANS.EQ.'R') RETURN
      GO TO (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17) N
      X(1) = ' (i) the solution vector Y = transpo'
      Y(1) = 'se [y(a), (py'')(a)] at a regular   '
      X(2) = '     end-point a, or                '
      Y(2) = ' '
      X(3) = ' (ii) the "singular solution vector"'
      Y(3) = ' Y(a) = transpose [[y,u](a),        '
      X(4) = '      [y,v](a)] at a singular LC end'
      Y(4) = ' -point a.  Similarly at the right  '
      X(5) = '      end-point b with U,V.         '
      Y(5) = ' '
      X(6) = '  The object of this section is to p'
      Y(6) = 'rovide help in choosing appropriate '
      X(7) = 'functions u and v in (**) (or U,V), '
      Y(7) = 'given the differential equation (*).'
      X(8) = 'Full details of the boundary conditi'
      Y(8) = 'ons for (*) are discussed in H7;    '
      X(9) = 'here it is sufficient to say that th'
      Y(9) = 'e limit-circle type boundary condi- '
      X(10) = 'tion (**) can be applied at any end-'
      Y(10) = 'point in the LCNO, LCO classifica-  '
      X(11) = 'tion, but also in the R, WR classifi'
      Y(11) = 'cation subject to the appropriate   '
      X(12) = 'choice of u,v or U,V.               '
      Y(12) = ' '
      X(13) = '  Let (*) be R, WR, LCNO, or LCO at '
      Y(13) = 'end-point a and choose c in (a,b).  '
      X(14) = 'Then either                         '
      Y(14) = ' '
      X(15) = '  u and v are a pair of linearly ind'
      Y(15) = 'ependent solutions of (*) on (a,c]  '
      X(16) = '  for any chosen real values of lamb'
      Y(16) = 'da, or                              '
      X(17) = '  u and v are a pair of real-valued '
      Y(17) = ' maximal domain functions defined on'
      X(18) = '  (a,c] satisfying [u,v](a) .ne. 0. '
      Y(18) = ' The maximal domain D(a,c] is de-   '
      X(19) = '  defined by                        '
      Y(19) = ' '
      X(20) = ' '
      Y(20) = ' '
      X(21) = '       D(a,c] = {f:(a,c]->R:: f,pf'' '
      Y(21) = 'in AC(a,c];                         '
      X(22) = '                   f, ((-pf'')''+qf)/w'
      Y(22) = ' in L2((a,c;w)}  .                  '
      X(23) = ' '
      Y(23) = ' '
      DO 601 I = 1,23
          WRITE (*,FMT=*) X(I),Y(I)
  601 CONTINUE
      READ (*,FMT=999) ANS,N
      WRITE (*,FMT=*)
      IF (ANS.EQ.'r' .OR. ANS.EQ.'R') RETURN
      GO TO (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17) N
      X(1) = '  It is known that for all f,g in D('
      Y(1) = 'a,c] the limit                      '
      X(2) = ' '
      Y(2) = ' '
      X(3) = '       [f,g](a) = lim[f,g](x) as x->'
      Y(3) = 'a                                   '
      X(4) = ' '
      Y(4) = ' '
      X(5) = 'exists and is finite.  If (*) is LCN'
      Y(5) = 'O or LCO at a, then all solutions of'
      X(6) = 'of (*) belong to D(a,c] for all valu'
      Y(6) = 'es of lambda.                       '
      X(7) = '  The boundary condition (**) is ess'
      Y(7) = 'ential in the LCNO and LCO cases but'
      X(8) = 'can also be used with advantage in s'
      Y(8) = 'ome R and WR cases.  In the R, WR, '
      X(9) = 'and LCNO cases, but not in the LCO c'
      Y(9) = 'ase, the boundary condition         '
      X(10) = 'functions can always be chosen so th'
      Y(10) = 'at                                  '
      X(11) = '        lim u(x)/v(x) = 0 as x->a,  '
      Y(11) = ' '
      X(12) = 'and it is recommended that this norm'
      Y(12) = 'alisation be effected, but this is  '
      X(13) = 'not essential; this normalisation ha'
      Y(13) = 's been entered in the examples given'
      X(14) = 'below.  In this case, the boundary c'
      Y(14) = 'ondition [y,u](a) = 0 (i.e., A1 = 1,'
      X(15) = 'A2 = 0 in (**) is called the princip'
      Y(15) = 'al or Friedrichs boundary condition '
      X(16) = 'at a.                              '
      Y(16) = ' '
      X(17) = '  In the case when end-points a and '
      Y(17) = 'b are, independently, in R, WR,     '
      X(18) = 'LCNO, or LCO classification, it may '
      Y(18) = 'be that symmetry or other reasons   '
      X(19) = 'permit one set of boundary condition'
      Y(19) = ' functions to be used at both end-  '
      X(20) = 'points (see xamples.f, #1 (Legendre)'
      Y(20) = ').  In other cases, different pairs '
      X(21) = 'must be chosen for each end-point (s'
      Y(21) = 'ee xamples.f: #16 (Jacobi),         '
      X(22) = '#18 (Dunsch), and #19 (Donsch)).    '
      Y(22) = ' '
      DO 602 I = 1,22
          WRITE (*,FMT=*) X(I),Y(I)
  602 CONTINUE
      WRITE (*,FMT=*)
      READ (*,FMT=999) ANS,N
      IF (ANS.EQ.'r' .OR. ANS.EQ.'R') RETURN
      GO TO (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17) N
      X(1) = ' '
      Y(1) = ' '
      X(2) = '  Note that a solution pair u,v is a'
      Y(2) = 'lways a maximal domain pair, but not'
      X(3) = 'necessarily vice versa.             '
      Y(3) = ' '
      X(4) = ' '
      Y(4) = ' '
      X(5) = 'EXAMPLES:                           '
      Y(5) = ' '
      X(6) = '1. -y''''(x) = lambda*y(x) on [0,pi] i'
      Y(6) = 's R at 0 and R at pi.               '
      X(7) = '   At 0, with lambda = 0, a solution'
      Y(7) = ' pair is u(x) = x, v(x) = 1.        '
      X(8) = '   At pi, with lambda = 1, a solutio'
      Y(8) = 'n pair is                           '
      X(9) = '     u(x) = sin(x), v(x) = cos(x).  '
      Y(9) = ' '
      X(10) = '2. -(sqrt(x)*y''(x))'' = lambda*y(x)/s'
      Y(10) = 'qrt(x) on (0,1] is                  '
      X(11) = '     WR at 0 and R at 1.            '
      Y(11) = ' '
      X(12) = '   (The general solutions of this eq'
      Y(12) = 'uation are                          '
      X(13) = '     u(x) = cos(2*sqrt(x*lambda)), v'
      Y(13) = '(x) = sin(2*sqrt(x*lambda)).)       '
      X(14) = '   At 0, with lambda = 0, a solution'
      Y(14) = ' pair is                            '
      X(15) = '     u(x) = 2*sqrt(x), v(x) = 1.    '
      Y(15) = ' '
      X(16) = '   At 1, with lambda = pi*pi/4, a so'
      Y(16) = 'lution pair is                      '
      X(17) = '     u(x) = sin(pi*sqrt(x)), v(x) = '
      Y(17) = 'cos(pi*sqrt(x)).                    '
      X(18) = '   At 1, with lambda = 0, a solution'
      Y(18) = ' pair is                            '
      X(19) = '     u(x) = 2*(1-sqrt(x)), v(x) = 1.'
      Y(19) = ' '
      X(20) = '   See also xamples.f, #10 (Weakly R'
      Y(20) = 'egular).                            '
      X(21) = '3. -((1-x*x)*y''(x))'' = lambda*y(x) o'
      Y(21) = 'n (-1,1) is LCNO at both ends.      '
      X(22) = '   At +-1, with lambda = 0, a soluti'
      Y(22) = 'on pair is                          '
      X(23) = '     u(x) = 1, v(x) = 0.5*log((1+x)/'
      Y(23) = '(1-x)).                             '
      DO 603 I = 1,23
          WRITE (*,FMT=*) X(I),Y(I)
  603 CONTINUE
      READ (*,FMT=999) ANS,N
      IF (ANS.EQ.'r' .OR. ANS.EQ.'R') RETURN
      GO TO (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17) N
      X(1) = '   At 1, a maximal domain pair is u('
      Y(1) = 'x) = 1, v(x) = log(1-x)             '
      X(2) = '   At -1, a maximal domain pair is u'
      Y(2) = '(x) = 1, v(x) = log(1+x).           '
      X(3) = '   See also xamples.f, #1 (Legendre)'
      Y(3) = '.                                   '
      X(4) = '4. -y''''(x) - (1/(4x*x))*y(x) = lambd'
      Y(4) = 'a*y(x) on (0,infinity) is           '
      X(5) = '     LCNO at 0 and LP at +infinity. '
      Y(5) = ' '
      X(6) = '   At 0, a maximal domain pair is   '
      Y(6) = ' '
      X(7) = '     u(x) = sqrt(x), v(x) = sqrt(x)*'
      Y(7) = 'log(x).                             '
      X(8) = '   See also xamples.f, #2 (Bessel). '
      Y(8) = ' '
      X(9) = '5. -y''''(x) - 5*(1/(4*x*x))*y(x) = la'
      Y(9) = 'mbda*y(x) on (0,infinity) is        '
      X(10) = '     LCO at 0 and LP at +infinity.  '
      Y(10) = ' '
      X(11) = '   At 0, with lambda = 0, a solution'
      Y(11) = ' pair is                            '
      X(12) = '     u(x) = sqrt(x)*cos(log(x)), v(x'
      Y(12) = ') = sqrt(x)*sin(log(x))             '
      X(13) = '   See also xamples.f, #20 (Krall). '
      Y(13) = ' '
      X(14) = '6. -y''''(x) - (1/x)*y(x) = lambda*y(x'
      Y(14) = ') on (0,infinity) is               '
      X(15) = '     LCNO at 0 and LP at +infinity.'
      Y(15) = ' '
      X(16) = '   At 0, a maximal domain pair is   '
      Y(16) = ' '
      X(17) = '     u(x) = x, v(x) = 1 -x*log(x).  '
      Y(17) = ' '
      X(18) = '   See also xamples.f, #4(Boyd).    '
      Y(18) = ' '
      X(19) = '7. -((1/x)*y''(x))'' + (k/(x*x) + k*k/'
      Y(19) = 'x)*y(x) = lambda*y(x) on (0,1],     '
      X(20) = '     with k real and .ne. 0, is LCNO'
      Y(20) = ' at 0 and R at 1.                   '
      X(21) = '   At 0, a maximal domain pair is   '
      Y(21) = ' '
      X(22) = '     u(x) = x*x, v(x) = x - 1/k.    '
      Y(22) = ' '
      X(23) = '   See also xamples.f, #8 (Laplace T'
      Y(23) = 'idal Wave).                         '
      DO 604 I = 1,23
          WRITE (*,FMT=*) X(I),Y(I)
  604 CONTINUE
      READ (*,FMT=999) ANS,N
      IF (ANS.EQ.'r' .OR. ANS.EQ.'R') RETURN
      GO TO (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17) N
    7 CONTINUE
      WRITE (*,FMT=*)
      WRITE (*,FMT=*)
      WRITE (*,FMT=*)
      WRITE (*,FMT=*)
      WRITE (*,FMT=*) 'H7:  General self-adjoint boundary conditions.'
      X(1) = '  Boundary conditions for Sturm-Liou'
      Y(1) = 'ville boundary value problems       '
      X(2) = ' '
      Y(2) = ' '
      X(3) = '                   -(p*y'')'' + q*y = '
      Y(3) = 'lambda*w*y              (*)         '
      X(4) = ' '
      Y(4) = ' '
      X(5) = 'on an interval (a,b) are either     '
      Y(5) = ' '
      X(6) = '  SEPARATED, with at most one condit'
      Y(6) = 'ion at end-point a and at most one  '
      X(7) = '    condition at end-point b, or    '
      Y(7) = ' '
      X(8) = '  COUPLED, when both a and b are, in'
      Y(8) = 'dependently, in one of the end-point'
      X(9) = '    classifications R, WR, LCNO, LCO'
      Y(9) = ', in which case two independent     '
      X(10) = '    boundary conditions are required'
      Y(10) = ' which link the solution values near'
      X(11) = '    a to those near b.              '
      Y(11) = ' '
      X(12) = '  The sleign2 program allows for all'
      Y(12) = ' self-adjoint boundary conditions:  '
      X(13) = 'separated self-adjoint conditions an'
      Y(13) = 'd all cases of coupled self-adjoint '
      X(14) = 'conditions.                         '
      Y(14) = ' '
      DO 701 I = 1,14
          WRITE (*,FMT=*) X(I),Y(I)
  701 CONTINUE
      WRITE (*,FMT=*)
      WRITE (*,FMT=*)
      WRITE (*,FMT=*)
      WRITE (*,FMT=*)
      WRITE (*,FMT=*)
      READ (*,FMT=999) ANS,N
      IF (ANS.EQ.'r' .OR. ANS.EQ.'R') RETURN
      GO TO (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17) N
      WRITE (*,FMT=*) '    Separated Conditions:      '
      WRITE (*,FMT=*) '    ---------------------      '
      X(1) = '  The boundary conditions to be sele'
      Y(1) = 'cted depend upon the classification '
      X(2) = 'of the differential equation at the '
      Y(2) = 'end-point, say, a:                  '
      X(3) = '  1. If the end-point a is LP, then '
      Y(3) = 'no boundary condition is required or'
      X(4) = '     allowed.                       '
      Y(4) = ' '
      X(5) = '  2. If the end-point a is R or WR, '
      Y(5) = 'then a separated boundary condition '
      X(6) = '     is of the form                 '
      Y(6) = ' '
      X(7) = '         A1*y(a) + A2*(py'')(a) = 0,'
      Y(7) = ' '
      X(8) = '     where A1, A2 are real constants'
      Y(8) = ' the user must choose, not both zero.  '
      X(9) = '  3. If the end-point a is LCNO or'
      Y(9) = ' LCO, then a separated boundary     '
      X(10) = '  condition is of the form          '
      Y(10) = ' '
      X(11) = '         A1*[y,u](a) + A2*[y,v](a) ='
      Y(11) = ' 0,                                 '
      X(12) = '  where A1, A2 are real constants th'
      Y(12) = 'e user must choose, not both zero;  '
      X(13) = '  here, u,v are the pair of boundary'
      Y(13) = ' condition functions you have       '
      X(14) = '  previously selected when the input'
      Y(14) = ' FORTRAN file was being prepared    '
      X(15) = '  with makepqw.f.                   '
      Y(15) = ' '
      X(16) = '  4. If the end-point a is LCNO an'
      Y(16) = 'd the boundary condition pair       '
      X(17) = '  u,v has been chosen so that       '
      Y(17) = ' '
      X(18) = '         lim u(x)/v(x) = 0  as x->a '
      Y(18) = ' '
      X(19) = '  (which is always possible), then A'
      Y(19) = '1 = 1, A2 = 0 (i.e., [y,u](a) = 0)  '
      X(20) = '  gives the principal (Friedrichs) b'
      Y(20) = 'oundary condition at a.             '
      DO 702 I = 1,20
          WRITE (*,FMT=*) X(I),Y(I)
  702 CONTINUE
      WRITE (*,FMT=*)
      WRITE (*,FMT=*)
      READ (*,FMT=999) ANS,N
      IF (ANS.EQ.'r' .OR. ANS.EQ.'R') RETURN
      GO TO (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17) N
      X(1) = '  5. If a is R or WR and boundary '
      Y(1) = 'condition functions u,v have been   '
      X(2) = '  entered in the FORTRAN input file,'
      Y(2) = ' then (3.,4.) above apply to        '
      X(3) = '  entering separated boundary condit'
      Y(3) = 'ions at such an end-point; the      '
      X(4) = '  boundary conditions in this form a'
      Y(4) = 're equivalent to the point-wise     '
      X(5) = '  conditions in (2.) (subject to car'
      Y(5) = 'e in choosing A1, A2).  This        '
      X(6) = '  singular form of a regular boundar'
      Y(6) = 'y condition may be particularly     '
      X(7) = '  effective in the WR case if the bo'
      Y(7) = 'undary condition form in (2.) leads '
      X(8) = '  to numerical difficulties.        '
      Y(8) = ' '
      X(9) = '  Conditions (2.,3.,4.,5.) apply sim'
      Y(9) = 'ilarly at end-point b (with U,V as  '
      X(10) = 'the boundary condition functions at '
      Y(10) = 'b.                                  '
      X(11) = '    6. If a is R, WR, LCNO, or LCO a'
      Y(11) = 'nd b is LP, then only a separated   '
      X(12) = '  condition at a is required and all'
      Y(12) = 'owed (or instead at b if a and b    '
      X(13) = '  are interchanged).                '
      Y(13) = ' '
      X(14) = '    7. If both end-points a and b ar'
      Y(14) = 'e LP, then no boundary conditions   '
      X(15) = '  are required or allowed.          '
      Y(15) = ' '
      X(16) = '  The indexing of eigenvalues for bo'
      Y(16) = 'undary value problems with separated'
      X(17) = '  conditions is discussed in H13.   '
      Y(17) = ' '
      X(18) = '    Coupled Conditions:             '
      Y(18) = ' '
      X(19) = '    -------------------             '
      Y(19) = ' '
      X(20) = '    8. Coupled regular self-adjoint '
      Y(20) = 'boundary conditions on (a,b) apply  '
      X(21) = 'only when both end-points a and b ar'
      Y(21) = 'e R or WR.                          '
      X(22) = ' '
      Y(22) = ' '
      DO 704 I = 1,22
          WRITE (*,FMT=*) X(I),Y(I)
  704 CONTINUE
      READ (*,FMT=999) ANS,N
      IF (ANS.EQ.'r' .OR. ANS.EQ.'R') RETURN
      GO TO (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17) N
    8 CONTINUE
      WRITE (*,FMT=*) 'H8:  Recording the results.'
      X(1) = '  If you choose to have a record kep'
      Y(1) = 't of the results, then the following'
      X(2) = 'information is stored in a file with'
      Y(2) = ' the name you select:               '
      X(3) = ' '
      Y(3) = ' '
      X(4) = '  1. The file name.                 '
      Y(4) = ' '
      X(5) = '  2. The header line prompted for (u'
      Y(5) = 'p to 32 characters of your choice). '
      X(6) = '  3. The interval (a,b) which was us'
      Y(6) = 'ed.                                 '
      X(7) = ' '
      Y(7) = ' '
      X(8) = '  For SEPARATED boundary conditions:'
      Y(8) = ' '
      X(9) = '  4. The end-point classification.  '
      Y(9) = ' '
      X(10) = '  5. A summary of coefficient inform'
      Y(10) = 'ation at WR, LCNO, LCO end-points.  '
      X(11) = '  6. The boundary condition constant'
      Y(11) = 's (A1,A2), (B1,B2) if entered.      '
      X(12) = '  7. (NUMEIG,EIG,TOL) or (NUMEIG1,NU'
      Y(12) = 'MEIG2,TOL), as entered.             '
      X(13) = ' '
      Y(13) = ' '
      X(14) = '  For COUPLED boundary conditions:  '
      Y(14) = ' '
      X(15) = '  8. The boundary condition paramete'
      Y(15) = 'r alpha and the coupling matrix K;  '
      X(16) = '     see H6.                        '
      Y(16) = ' '
      X(17) = '  For ALL self-adjoint boundary cond'
      Y(17) = 'itions: '
      X(18) = '  9. The computed eigenvalue, EIG, a'
      Y(18) = 'nd its estimated accuracy, TOL.     '
      X(19) = ' 10. IFLAG reported (see H15).      '
      Y(19) = ' '
      DO 801 I = 1,19
          WRITE (*,FMT=*) X(I),Y(I)
  801 CONTINUE
      WRITE (*,FMT=*)
      WRITE (*,FMT=*)
      WRITE (*,FMT=*)
      WRITE (*,FMT=*) ' '
      READ (*,FMT=999) ANS,N
      IF (ANS.EQ.'r' .OR. ANS.EQ.'R') RETURN
      GO TO (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17) N
    9 CONTINUE
      WRITE (*,FMT=*)
      WRITE (*,FMT=*)
      WRITE (*,FMT=*)
      WRITE (*,FMT=*)
      WRITE (*,FMT=*)
      WRITE (*,FMT=*) 'H9:  Type and choice of interval.'
      X(1) = '  You may enter any interval (a,b) f'
      Y(1) = 'or which the coefficients p,q,w are '
      X(2) = 'well-defined by your FORTRAN stateme'
      Y(2) = 'nts in the input file, provided that'
      X(3) = '(a,b) contains no interior singulari'
      Y(3) = 'ties.                               '
      DO 901 I = 1,3
          WRITE (*,FMT=*) X(I),Y(I)
  901 CONTINUE
      WRITE (*,FMT=*)
      WRITE (*,FMT=*)
      WRITE (*,FMT=*)
      WRITE (*,FMT=*)
      WRITE (*,FMT=*)
      WRITE (*,FMT=*)
      WRITE (*,FMT=*)
      WRITE (*,FMT=*)
      WRITE (*,FMT=*)
      WRITE (*,FMT=*)
      WRITE (*,FMT=*)
      WRITE (*,FMT=*)
      WRITE (*,FMT=*)
      WRITE (*,FMT=*)
      WRITE (*,FMT=*) ' '
      READ (*,FMT=999) ANS,N
      IF (ANS.EQ.'r' .OR. ANS.EQ.'R') RETURN
      GO TO (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17) N
   10 CONTINUE
      WRITE (*,FMT=*)
      WRITE (*,FMT=*)
      WRITE (*,FMT=*)
      WRITE (*,FMT=*)
      WRITE (*,FMT=*)
      WRITE (*,FMT=*)
      WRITE (*,FMT=*)
      WRITE (*,FMT=*)
      WRITE (*,FMT=*) 'H10:  Entry of end-points.'
      X(1) = '  End-points a and b should generall'
      Y(1) = 'y be entered as real numbers to an  '
      X(2) = 'appropriate number of decimal places'
      Y(2) = '. '
      DO 1001 I = 1,2
          WRITE (*,FMT=*) X(I),Y(I)
 1001 CONTINUE
      WRITE (*,FMT=*)
      WRITE (*,FMT=*)
      WRITE (*,FMT=*)
      WRITE (*,FMT=*)
      WRITE (*,FMT=*)
      WRITE (*,FMT=*)
      WRITE (*,FMT=*)
      WRITE (*,FMT=*)
      WRITE (*,FMT=*)
      WRITE (*,FMT=*)
      WRITE (*,FMT=*)
      WRITE (*,FMT=*)
      WRITE (*,FMT=*) ' '
      READ (*,FMT=999) ANS,N
      IF (ANS.EQ.'r' .OR. ANS.EQ.'R') RETURN
      GO TO (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17) N
   11 CONTINUE
      WRITE (*,FMT=*)
      WRITE (*,FMT=*)
      WRITE (*,FMT=*)
      WRITE (*,FMT=*)
      WRITE (*,FMT=*)
      WRITE (*,FMT=*)
      WRITE (*,FMT=*)
      WRITE (*,FMT=*)
      WRITE (*,FMT=*) 'H11:  End-point values of p,q,w.'
      X(1) = '  The program sleign2 needs to know '
      Y(1) = 'whether the coefficient functions   '
      X(2) = 'p(x),q(x),w(x) defined by the FORTRA'
      Y(2) = 'N expressions entered in the input  '
      X(3) = 'file can be evaluated numerically wi'
      Y(3) = 'thout running into difficulty.  If, '
      X(4) = 'for example, either q or w is unboun'
      Y(4) = 'ded at a, or p(a) is 0, then sleign2'
      X(5) = 'needs to know this so that a is not '
      Y(5) = 'chosen for functional evaluation.   '
      DO 1101 I = 1,5
          WRITE (*,FMT=*) X(I),Y(I)
 1101 CONTINUE
      WRITE (*,FMT=*)
      WRITE (*,FMT=*)
      WRITE (*,FMT=*)
      WRITE (*,FMT=*)
      WRITE (*,FMT=*)
      WRITE (*,FMT=*)
      WRITE (*,FMT=*)
      WRITE (*,FMT=*)
      WRITE (*,FMT=*)
      WRITE (*,FMT=*) ' '
      READ (*,FMT=999) ANS,N
      IF (ANS.EQ.'r' .OR. ANS.EQ.'R') RETURN
      GO TO (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17) N
   12 CONTINUE
      WRITE (*,FMT=*) 'H12:  Initial value problems.'
      X(1) = '  The initial value problem facility'
      Y(1) = ' for Sturm-Liouville problems       '
      X(2) = ' '
      Y(2) = ' '
      X(3) = '                   -(p*y'')'' + q*y = '
      Y(3) = ' lambda*w*y              (*)        '
      X(4) = ' '
      Y(4) = ' '
      X(5) = 'allows for the computation of a solu'
      Y(5) = 'tion of (*) with a user-chosen      '
      X(6) = 'value lambda and any one of the foll'
      Y(6) = 'owing initial conditions:           '
      X(7) = '  1. From end-point a of any classif'
      Y(7) = 'ication except LP towards           '
      X(8) = 'end-point b of any classification,  '
      Y(8) = ' '
      X(9) = '  2. From end-point b of any classif'
      Y(9) = 'ication except LP back towards      '
      X(10) = 'end-point a of any classification,  '
      Y(10) = ' '
      X(11) = '  3. From end-points a and b of any '
      Y(11) = 'classifications except LP towards an'
      X(12) = 'interior point of (a,b) selected by '
      Y(12) = 'the program.                        '
      X(13) = ' '
      Y(13) = ' '
      X(14) = '  Initial values at a are of the for'
      Y(14) = 'm y(a) = alpha1, (p*y'')(a) =alpha2,'
      X(15) = 'when a is R or WR; and [y,u](a) = al'
      Y(15) = 'pha1, [y,v](a) = alpha2, when a is  '
      X(16) = 'LCNO or LCO.                        '
      Y(16) = ' '
      X(17) = '  Initial values at b are of the for'
      Y(17) = 'm y(b) = beta1, (p*y'')(b) = beta2, '
      X(18) = 'when b is R or WR; and [y,u](b) = be'
      Y(18) = 'ta1, [y,v](b) = beta2, when b is    '
      X(19) = 'LCNO or LCO.                        '
      Y(19) = ' '
      X(20) = '  In (*), lambda is a user-chosen re'
      Y(20) = 'al number; while in the above ini-  '
      X(21) = 'tial values, (alpha1,alpha2) and (be'
      Y(21) = 'ta1,beta2) are user-chosen pairs of '
      X(22) = 'real numbers not both zero.         '
      Y(22) = ' '
      DO 1201 I = 1,22
          WRITE (*,FMT=*) X(I),Y(I)
 1201 CONTINUE
      WRITE (*,FMT=*)
      READ (*,FMT=999) ANS,N
      IF (ANS.EQ.'r' .OR. ANS.EQ.'R') RETURN
      GO TO (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17) N
      WRITE (*,FMT=*)
      WRITE (*,FMT=*)
      WRITE (*,FMT=*)
      WRITE (*,FMT=*)
      WRITE (*,FMT=*)
      X(1) = '  In the initial value case (3.) abo'
      Y(1) = 've when the interval (a,b) is       '
      X(2) = 'finite, the interior point selected '
      Y(2) = 'by the program is generally near the'
      X(3) = 'midpoint of (a,b); when (a,b) is inf'
      Y(3) = 'inite, no general rule can be given.'
      X(4) = 'Also if, given (alpha1,alpha2) and ('
      Y(4) = 'beta1,beta2), the lambda chosen is  '
      X(5) = 'an eigenvalue of the associated boun'
      Y(5) = 'dary value problem, the computed    '
      X(6) = 'solution may not be the correspondin'
      Y(6) = 'g eigenfunction -- the signs of the '
      X(7) = 'computed solutions on either side of'
      Y(7) = ' the interior point may be opposite.'
      X(8) = '  The output for a solution of an in'
      Y(8) = 'itial value problem is in the form  '
      X(9) = 'of stored numerical data which can b'
      Y(9) = 'e plotted on the screen (see H16),  '
      X(10) = 'or printed out in graphical form if '
      Y(10) = 'graphics software is available.     '
      DO 1202 I = 1,10
          WRITE (*,FMT=*) X(I),Y(I)
 1202 CONTINUE
      WRITE (*,FMT=*)
      WRITE (*,FMT=*)
      WRITE (*,FMT=*)
      WRITE (*,FMT=*)
      WRITE (*,FMT=*)
      WRITE (*,FMT=*)
      WRITE (*,FMT=*)
      WRITE (*,FMT=*) ' '
      READ (*,FMT=999) ANS,N
      IF (ANS.EQ.'r' .OR. ANS.EQ.'R') RETURN
      GO TO (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17) N
   13 CONTINUE
      WRITE (*,FMT=*) 'H13:  Indexing of eigenvalues.'
      X(1) = '  The indexing of eigenvalues is an '
      Y(1) = 'automatic facility in sleign2.  The '
      X(2) = 'following general results hold for t'
      Y(2) = 'he separated boundary condition     '
      X(3) = 'problem (see H7):                   '
      Y(3) = ' '
      X(4) = '  1. If neither end-point a or b is '
      Y(4) = 'LP or LCO, then the spectrum of the '
      X(5) = 'eigenvalue problem is discrete (eige'
      Y(5) = 'nvalues only), simple (eigenvalues  '
      X(6) = 'all of multiplicity 1), and bounded '
      Y(6) = 'below with a single cluster point at'
      X(7) = '+infinity.  The eigenvalues are inde'
      Y(7) = 'xed as {lambda(n): n=0,1,2,...},    '
      X(8) = 'where lambda(n) < lambda(n+1) (n=0,1'
      Y(8) = ',2,...), lim lambda(n) -> +infinity;'
      X(9) = 'and if {psi(n): n=0,1,2,...} are the'
      Y(9) = ' corresponding eigenfunctions, then '
      X(10) = 'psi(n) has exactly n zeros in the op'
      Y(10) = 'en interval (a,b).                  '
      X(11) = '  2. If neither end-point a or b is '
      Y(11) = 'LP but at least one end-point is    '
      X(12) = 'LCO, then the spectrum is discrete a'
      Y(12) = 'nd simple as for (1.), but with     '
      X(13) = 'cluster points at both +infinity and'
      Y(13) = ' -infinity.  The eigenvalues are in-'
      X(14) = 'dexed as {lambda(n): n=0,1,-1,2,-2,.'
      Y(14) = '..}, where                          '
      X(15) = 'lambda(n) < lambda(n+1) (n=...-2,-1,'
      Y(15) = '0,1,2,...) with lambda(0) the small-'
      X(16) = 'est non-negative eigenvalue and lim '
      Y(16) = 'lambda(n) -> +infinity or -> -infi- '
      X(17) = 'nity with n; and if {psi(n): n=0,1,-'
      Y(17) = '1,2,-2,...} are the corresponding   '
      X(18) = 'eigenfunctions, then every psi(n) ha'
      Y(18) = 's infinitely many zeros in (a,b).   '
      DO 1301 I = 1,18
          WRITE (*,FMT=*) X(I),Y(I)
 1301 CONTINUE
      WRITE (*,FMT=*)
      WRITE (*,FMT=*)
      WRITE (*,FMT=*)
      WRITE (*,FMT=*)
      READ (*,FMT=999) ANS,N
      IF (ANS.EQ.'r' .OR. ANS.EQ.'R') RETURN
      GO TO (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17) N
      X(1) = '  3. If one or both end-points is LP'
      Y(1) = ', then there can be one or more in- '
      X(2) = 'tervals of continuous spectrum for t'
      Y(2) = 'he boundary value problem in addi-  '
      X(3) = 'tion to some (necessarily simple) ei'
      Y(3) = 'genvalues.  For these essentially   '
      X(4) = 'more difficult problems, sleign2 can'
      Y(4) = ' be used as an investigative tool to'
      X(5) = 'give qualitative and possibly quanti'
      Y(5) = 'tative information on the spectrum. '
      X(6) = '     For example, if a problem has a'
      Y(6) = ' continuous spectrum starting at L, '
      X(7) = 'then there may be no eigenvalue belo'
      Y(7) = 'w L, any finite number of eigenval- '
      X(8) = 'ues below L, or an infinite (but cou'
      Y(8) = 'ntable) number of eigenvalues below '
      X(9) = 'L. sleign2 can be used to compute L '
      Y(9) = '(see the paper bewz on the sleign2  '
      X(10) = 'home page for an algorithm to comput'
      Y(10) = 'e L), and determine the number of   '
      X(11) = 'these eigenvalues and compute them. '
      Y(11) = ' In this respect, see xamples.f: #13'
      X(12) = '(Hydrogen Atom), #17 (Morse Oscillat'
      Y(12) = 'or), #21 (Fourier), #27 (Joergens)  '
      X(13) = 'as examples of success; and #12 (Mat'
      Y(13) = 'hieu), #14 (Marletta), and #28      '
      X(14) = '(Behnke-Goerisch) as examples of fai'
      Y(14) = 'lure.                               '
      X(15) = '     The problem need not have a con'
      Y(15) = 'tinuous spectrum, in which case if  '
      X(16) = 'its discrete spectrum is bounded bel'
      Y(16) = 'ow, then the eigenvalues are indexed'
      X(17) = 'and the eigenfunctions have zero cou'
      Y(17) = 'nts as in (1.).  If, on the other   '
      X(18) = 'hand, the discrete spectrum is unbou'
      Y(18) = 'nded below, then all the eigenfunc- '
      X(19) = 'tions have infinitely many zeros in '
      Y(19) = 'the open interval (a,b).  sleign2   '
      X(20) = 'can, in principle, compute these eig'
      Y(20) = 'envalues if neither end-point is LP,'
      X(21) = 'although this is a computationally d'
      Y(21) = 'ifficult problem.  Note, however,   '
      X(22) = 'that sleign2 has no algorithm when t'
      Y(22) = 'he spectrum is discrete, unbounded  '
      X(23) = 'above and below, and one end-point i'
      Y(23) = 's LP, as in xamples.f #30.          '
      DO 1302 I = 1,23
          WRITE (*,FMT=*) X(I),Y(I)
 1302 CONTINUE
      WRITE (*,FMT=*)
      READ (*,FMT=999) ANS,N
      IF (ANS.EQ.'r' .OR. ANS.EQ.'R') RETURN
      GO TO (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17) N
      X(1) = '  In respect to the three different '
      Y(1) = 'types of indexing discussed above,  '
      X(2) = 'the following identified examples fr'
      Y(2) = 'om xamples.f illustrate the spectral'
      X(3) = 'property of these boundary problems:'
      Y(3) = '                                    '
      X(4) = '  1. Neither end-point is LP or LCO.'
      Y(4) = ' '
      X(5) = '       #1 (Legendre)                '
      Y(5) = ' '
      X(6) = '       #2 (Bessel) with -1/4 < c < 3'
      Y(6) = '/4                                  '
      X(7) = '       #4 (Boyd)                    '
      Y(7) = ' '
      X(8) = '       #5 (Latzko)                  '
      Y(8) = ' '
      X(9) = '  2. Neither end-point is LP, but at'
      Y(9) = ' least one is LCO.                  '
      X(10) = '       #6 (Sears-Titchmarsh)        '
      Y(10) = ' '
      X(11) = '       #7 (BEZ)                     '
      Y(11) = ' '
      X(12) = '      #19 (Donsch)                  '
      Y(12) = ' '
      X(13) = '  3. At least one end-point is LP.  '
      Y(13) = ' '
      X(14) = '      #13 (Hydrogen Atom)           '
      Y(14) = ' '
      X(15) = '      #14 (Marletta)                '
      Y(15) = ' '
      X(16) = '      #20 (Krall)                   '
      Y(16) = ' '
      X(17) = '      #21 (Fourier) on [0,infinity) '
      Y(17) = ' '
      DO 1303 I = 1,17
          WRITE (*,FMT=*) X(I),Y(I)
 1303 CONTINUE
      WRITE (*,FMT=*)
      WRITE (*,FMT=*)
      WRITE (*,FMT=*)
      WRITE (*,FMT=*)
      WRITE (*,FMT=*)
      WRITE (*,FMT=*) ' '
      READ (*,FMT=999) ANS,N
      IF (ANS.EQ.'r' .OR. ANS.EQ.'R') RETURN
      GO TO (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17) N
   14 CONTINUE
      WRITE (*,FMT=*) 'H14:  Entry of eigenvalue index, initial guess,'
     +  //' and tolerance.'
      X(1) = '  For all self-adjoint boundary cond'
      Y(1) = 'ition problems (see H7), sleign2    '
      X(2) = 'calls for input information options '
      Y(2) = 'to compute either                   '
      X(3) = '  1. a single eigenvalue, or        '
      Y(3) = ' '
      X(4) = '  2. a series of eigenvalues.       '
      Y(4) = ' '
      X(5) = 'In each case indexing of eigenvalues'
      Y(5) = ' is called for (see H13).           '
      X(6) = '  (1.) above asks for data triples N'
      Y(6) = 'UMEIG, EIG, TOL separated by commas.'
      X(7) = 'Here NUMEIG is the integer index of '
      Y(7) = 'the desired eigenvalue; NUMEIG can  '
      X(8) = 'be negative only when the problem is'
      Y(8) = ' LCO at one or both end-points.     '
      X(9) = 'EIG allows for the entry of an initi'
      Y(9) = 'al guess for the requested          '
      X(10) = 'eigenvalue (if an especially good on'
      Y(10) = 'e is available), or can be set to 0 '
      X(11) = 'in which case an initial guess is ge'
      Y(11) = 'nerated by sleign2 itself.          '
      X(12) = 'TOL is the desired accuracy of the c'
      Y(12) = 'omputed eigenvalue.  It is an       '
      X(13) = 'absolute accuracy if the magnitude o'
      Y(13) = 'f the eigenvalue is 1 or less, and  '
      X(14) = 'is a relative accuracy otherwise.  T'
      Y(14) = 'ypical values might be .001 for     '
      X(15) = 'moderate accuracy and .0000001 for h'
      Y(15) = 'igh accuracy in single precision.   '
      X(16) = 'If TOL is set to 0, the maximum achi'
      Y(16) = 'evable accuracy is requested.       '
      X(17) = '  If the input data list is truncate'
      Y(17) = 'd with a "/" after NUMEIG or EIG,   '
      X(18) = 'then the remaining elements default '
      Y(18) = 'to 0.                               '
      DO 1401 I = 1,18
          WRITE (*,FMT=*) X(I),Y(I)
 1401 CONTINUE
      WRITE (*,FMT=*)
      WRITE (*,FMT=*)
      WRITE (*,FMT=*)
      WRITE (*,FMT=*)
      WRITE (*,FMT=*)
      READ (*,FMT=999) ANS,N
      IF (ANS.EQ.'r' .OR. ANS.EQ.'R') RETURN
      GO TO (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17) N
      X(1) = '  (2.) above asks for data triples N'
      Y(1) = 'UMEIG1, neig2, TOL separated by   '
      X(2) = 'commas.  Here numeig1 and numeig2 ar'
      Y(2) = 'e the first and last integer indices'
      X(3) = 'of the sequence of desired eigenvalu'
      Y(3) = 'es, numeig1 < numeig2; they can be  '
      X(4) = 'negative only when the problem is LC'
      Y(4) = 'O at one or both end-points.        '
      X(5) = 'TOL is the desired accuracy of the c'
      Y(5) = 'omputed eigenvalues.  It is an      '
      X(6) = 'absolute accuracy if the magnitude o'
      Y(6) = 'f an eigenvalue is 1 or less, and   '
      X(7) = 'is a relative accuracy otherwise.  T'
      Y(7) = 'ypical values might be .001 for     '
      X(8) = 'moderate accuracy and .0000001 for h'
      Y(8) = 'igh accuracy in single precision.   '
      X(9) = 'If TOL is set to 0, the maximum achi'
      Y(9) = 'evable accuracy is requested.       '
      X(10) = '  If the input data list is truncate'
      Y(10) = 'd with a "/" after neig2, then TOL'
      X(11) = 'defaults to 0.                      '
      Y(11) = ' '
      X(12) = '  For COUPLED self-adjoint boundary '
      Y(12) = 'condition problems (see H7 and H17),'
      X(13) = 'sleign2 also reports which eigenvalu'
      Y(13) = 'es are double.  Double eigenvalues  '
      X(14) = 'can occur only for coupled boundary '
      Y(14) = 'conditions with alpha = 0 or pi.    '
      DO 1402 I = 1,14
          WRITE (*,FMT=*) X(I),Y(I)
 1402 CONTINUE
      WRITE (*,FMT=*)
      WRITE (*,FMT=*)
      WRITE (*,FMT=*)
      WRITE (*,FMT=*)
      WRITE (*,FMT=*)
      WRITE (*,FMT=*)
      WRITE (*,FMT=*)
      WRITE (*,FMT=*) ' '
      READ (*,FMT=999) ANS,N
      IF (ANS.EQ.'r' .OR. ANS.EQ.'R') RETURN
      GO TO (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17) N
   15 CONTINUE
      WRITE (*,FMT=*)
      WRITE (*,FMT=*)
      WRITE (*,FMT=*)
      WRITE (*,FMT=*)
      WRITE (*,FMT=*)
      WRITE (*,FMT=*) 'H15:  IFLAG information.'
      X(1) = '  All results are reported by sleign2'
      Y(1) = '2 with a flag identification.  There'
      X(2) = 'are four values of IFLAG:           '
      Y(2) = ' '
      X(3) = ' '
      Y(3) = ' '
      X(4) = '  1 - The computed eigenvalue has an'
      Y(4) = ' estimated accuracy within the      '
      X(5) = '      tolerance requested.          '
      Y(5) = ' '
      X(6) = ' '
      Y(6) = ' '
      X(7) = '  2 - The computed eigenvalue does n'
      Y(7) = 'ot have an estimated accuracy within'
      X(8) = '      the tolerance requested, but i'
      Y(8) = 's the best the program could obtain.'
      X(9) = ' '
      Y(9) = ' '
      X(10) = '  3 - There seems to be no eigenvalu'
      Y(10) = 'e of index equal to NUMEIG.         '
      X(11) = ' '
      Y(11) = ' '
      X(12) = '  4 - The program has been unable to'
      Y(12) = ' compute the requested eigenvalue.  '
      DO 1501 I = 1,12
          WRITE (*,FMT=*) X(I),Y(I)
 1501 CONTINUE
      WRITE (*,FMT=*)
      WRITE (*,FMT=*)
      WRITE (*,FMT=*)
      WRITE (*,FMT=*)
      WRITE (*,FMT=*) ' '
      READ (*,FMT=999) ANS,N
      IF (ANS.EQ.'r' .OR. ANS.EQ.'R') RETURN
      GO TO (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17) N
   16 CONTINUE
      WRITE (*,FMT=*) 'H16:  Plotting.'
      X(1) = '  After computing a single eigenvalu'
      Y(1) = 'e (see H14(1.)), but not a sequence '
      X(2) = 'of eigenvalues (see H14(2.)), the ei'
      Y(2) = 'genfunction can be plotted for sepa-'
      X(3) = 'rated conditions and for coupled one'
      Y(3) = 's with alpha = 0 or pi.  When alpha '
      X(4) = 'is not zero or pi the eigenfunctions'
      Y(4) = ' are not real-valued and so no plot '
      X(5) = 'is provided.  If this is desired, re'
      Y(5) = 'spond "y" when asked so that sleign2'
      X(6) = 'will compute some eigenfunction data'
      Y(6) = ' and store them.                    '
      X(7) = '  One can ask that the eigenfunction'
      Y(7) = ' data be in the form of either      '
      X(8) = 'points (x,y) for x in (a,b), or poin'
      Y(8) = 'ts (t,y) for t in the standardized  '
      X(9) = 'interval (-1,1) mapped onto from (a,'
      Y(9) = 'b); the t- choice can be especially '
      X(10) = 'helpful when the original interval i'
      Y(10) = 's infinite.  Additionally, one can  '
      X(11) = 'ask for a plot of the so-called Prue'
      Y(11) = 'fer angle, in x- or t- variables.   '
      X(12) = '  In both forms, once the choice has'
      Y(12) = ' been made of the function to be    '
      X(13) = 'plotted, a crude plot is displayed o'
      Y(13) = 'n the monitor screen and you are    '
      X(14) = 'asked whether you wish to save the c'
      Y(14) = 'omputed plot points in a file.      '
      DO 1601 I = 1,14
          WRITE (*,FMT=*) X(I),Y(I)
 1601 CONTINUE
      WRITE (*,FMT=*)
      WRITE (*,FMT=*)
      WRITE (*,FMT=*)
      WRITE (*,FMT=*)
      WRITE (*,FMT=*)
      WRITE (*,FMT=*)
      WRITE (*,FMT=*)
      WRITE (*,FMT=*)
      WRITE (*,FMT=*) ' '
      READ (*,FMT=999) ANS,N
      IF (ANS.EQ.'r' .OR. ANS.EQ.'R') RETURN
      GO TO (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17) N
   17 CONTINUE
      WRITE (*,FMT=*) 'H17:  Indexing of eigenvalues for'//
     +  ' coupled self-adjoint problems.'
      X(1) = '  The indexing of eigenvalues is an '
      Y(1) = 'automatic facility in sleign2.  The '
      X(2) = 'following general result holds for c'
      Y(2) = 'oupled boundary condition problems  '
      X(3) = '(see H7):                           '
      Y(3) = ' '
      X(4) = '  The spectrum of the eigenvalue pro'
      Y(4) = 'blem is discrete (eigenvalues only).'
      X(5) = 'In general the spectrum is not simpl'
      Y(5) = 'e, but no eigenvalue exceeds multi- '
      X(6) = 'plicity 2.The eigenvalues are indexe'
      Y(6) = 'd as {lambda(n): n=0,1,2,...}, where'
      X(7) = 'lambda(n) .le. lambda(n+1) (n=0,1,2,'
      Y(7) = '...), lim lambda(n) -> +infinity if '
      X(8) = 'neither end-point is LCO.  If one or'
      Y(8) = ' both end-points are LCO the eigen- '
      X(9) = 'values cluster at both +infinity and'
      Y(9) = ' at -infinity, and all eigenfunc-   '
      X(10) = 'tions have infinitely many zeros.   '
      Y(10) = '                                    '
      X(11) = '  If neither end-point is LCO and al'
      Y(11) = 'pha = 0 or pi, then the n-th eigen- '
      X(12) = 'function has n-1, n, or n+1 zeros in'
      Y(12) = ' the half-open interval [a,b) (also '
      X(13) = 'in (a,b]).  All three possibilities '
      Y(13) = 'occur.  Recall that in the case of  '
      X(14) = 'double eigenvalues, although the n-t'
      Y(14) = 'h eigenvalue is well-defined, there '
      X(15) = 'is an ambiguity about which solution'
      Y(15) = ' is declared the n-th eigenfunction.'
      X(16) = '  If alpha is not 0 or pi, then the '
      Y(16) = 'eigenfunction is non-real and has no'
      X(17) = 'zero in (a,b); but each of the real '
      Y(17) = 'and imaginary parts of the n-th ei- '
      X(18) = 'genfunction have the same zero prope'
      Y(18) = 'rties as mentioned above when alpha '
      X(19) = '= 0 or pi.                          '
      Y(19) = ' '
      DO 1651 I = 1,19
          WRITE (*,FMT=*) X(I),Y(I)
 1651 CONTINUE
      READ (*,FMT=999) ANS,N
      IF (ANS.EQ.'r' .OR. ANS.EQ.'R') RETURN
      GO TO (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17) N
      WRITE (*,FMT=*)
      WRITE (*,FMT=*)
      WRITE (*,FMT=*)
      WRITE (*,FMT=*)
      WRITE (*,FMT=*)
      WRITE (*,FMT=*)
      X(1) = '  The following identified examples '
      Y(1) = 'from xamples.f are of special inter-'
      X(2) = 'est:                                '
      Y(2) = ' '
      X(3) = '    #11 (Plum) on [0,pi]           '
      Y(3) = ' '
      X(4) = '    #21 (Fourier) on [0,pi]        '
      Y(4) = ' '
      X(5) = '    #25 (Meissner) on [-0.5,0.5]   '
      Y(5) = ' '
      DO 1701 I = 1,5
          WRITE (*,FMT=*) X(I),Y(I)
 1701 CONTINUE
      WRITE (*,FMT=*) ' '
      WRITE (*,FMT=*)
      WRITE (*,FMT=*)
      WRITE (*,FMT=*)
      WRITE (*,FMT=*)
      WRITE (*,FMT=*)
      WRITE (*,FMT=*)
      WRITE (*,FMT=*)
      WRITE (*,FMT=*)
      WRITE (*,FMT=*)
      WRITE (*,FMT=*)
      WRITE (*,FMT=*)
      READ (*,FMT=999) ANS,N
      IF (ANS.EQ.'r' .OR. ANS.EQ.'R') RETURN
      GO TO (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17) N
      GO TO 1
  999 FORMAT (A1,I2)
      END
C
      SUBROUTINE PQ(END,P0,QF)
C  THIS PROGRAM IS INTENDED TO DETERMINE NUMERICALLY WHETHER OR
C    NOT THE FUNCTION P(X) IS ZERO AT A NON-REGULAR ENDPOINT,
C    AND WHETHER OR NOT THE FUNCTION Q(X)/W(X) IS INFINITE THERE.
C    IF P IS ZERO, P0 IS SET TO +1.0; ELSE TO -1.0.
C    IF Q AND/OR W IS FINITE, QF IS SET TO +1.0; ELSE TO -1.0.
C     .. Scalar Arguments ..
      REAL END,P0,QF
C     ..
C     .. Local Scalars ..
      REAL DT,HALF,ONE,PV,QV,QVP,QVSAV,T,TMP,TSAV,X,ZER
      INTEGER I
C     ..
C     .. External Functions ..
      REAL P,Q,W
      EXTERNAL P,Q,W
C     ..
C     .. External Subroutines ..
      EXTERNAL DXDT
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS
C     ..
      ONE = 1.0D0
      HALF = 0.5D0
      ZER = 0.0D0
c
      P0 = -1.0D0
      QF = 1.0D0
      DT = 0.1D0
      PV = 1.0D0
      IF (END.LT.0.0D0) DT = -DT
      DO 10 I = 1,25
          DT = HALF*DT
          T = END - DT
          IF (T.EQ.-1.0D0) THEN
              IF (PV.LT.0.0001D0) THEN
                  P0 = 1.0D0
                  GO TO 15
              ELSE
                  GO TO 15
              END IF
          END IF
          CALL DXDT(T,TMP,X)
          PV = ABS(P(X))
          IF (PV.GT. (1.D+5)) THEN
              GO TO 15
          ELSE IF (PV.LT. (.0001D0) .AND. I.GT.15) THEN
              P0 = 1.0D0
              GO TO 15
          END IF
   10 CONTINUE
   15 CONTINUE
c
      DT = 0.1D0
      IF (END.LT.0.0D0) DT = -DT
      QVP = ONE
      QVSAV = ZER
      TSAV = END
      QV = 1.0D0
      DO 20 I = 1,17
          DT = HALF*DT
          T = END - DT
          IF (T.EQ.1.0D0 .AND. (QV.GT. (1.D+5))) THEN
              QF = -1.0D0
              GO TO 25
          END IF
          CALL DXDT(T,TMP,X)
          QV = ABS(Q(X)) + ABS(W(X))
          QVP = ABS((QV-QVSAV)/ (T-TSAV))
          QVSAV = QV
          IF (QVP.GT. (1.D+5)) THEN
              QF = -1.0D0
              GO TO 25
          END IF
   20 CONTINUE
   25 CONTINUE
c
      RETURN
      END
