      PROGRAM DRIVE
C
C     **********
C     OCTOBER 15, 1995; P.B. BAILEY, W.N. EVERITT, B. GARBOW AND A. ZETTL
C     **********
C     This program solves the boundary value problem defined by
C     the differential equation
C
C           -(py')' + q*y = lambda*w*y
C
C     and appropriate boundary conditions.  If the problem is regular,
C     the boundary conditions can be of periodic type.
C
C     Program usage is facilitated by a large HELP subroutine.
C
C     **********
C     .. Scalars in Common ..
      INTEGER IND,MDTHZ
      LOGICAL ADDD
      DOUBLE PRECISION AA,BB,CC,DTHDAA,DTHDBB,EPSMIN,HPI,PI,THETU,THETV,
     1     TMID,TWOPI,UB,UL,UR,VB,VL,VR,Z
C     ..
C     .. Local Scalars ..
      INTEGER I,I1,I2,IFLAG,INTAB,NANS,NEND,NIVP,NUMEIG,NUMEIG1,NUMEIG2
      LOGICAL REGA,REGB,EIGV,PEIGF,PERIOD,RITE,SINGA,SINGB,
     1        SKIPB,WREGA,WREGB,YEH
      CHARACTER*1 ANSCH,HQ,YN
      CHARACTER*9 INFM,INFP
      CHARACTER*19 FILLA,FILLB
      CHARACTER*32 CHA,CHB,CHANS,CH1,CH2,CH3,CH4,CH5,CH6,FMT,TAPE2
      CHARACTER*70 CHTXT
      DOUBLE PRECISION A,ALFA,ALFA1,ALFA2,A1,A2,B,BETA,BETA1,BETA2,
     1     B1,B2,CIRCLA,CIRCLB,EIG,OSCILA,OSCILB,P0ATA,P0ATB,
     2     QFATA,QFATB,SINGATA,SINGATB,THA,THB,TOL,TOLL
C     ..
C     .. Local Arrays ..
      INTEGER ICOL(2),IIS(50)
      CHARACTER*2 COL(32)
      CHARACTER*39 BLNK(2),STAR(2),STR(2)
      DOUBLE PRECISION EES(50),SLFN(9),TTS(50)
C     ..
C     .. External Subroutines ..
      EXTERNAL DRAW,EXAMP,HELP,LSTDIR,SLEIGN2,PERIO
C     ..
C     .. External Functions ..
      CHARACTER*32 FMT2
      EXTERNAL FMT2
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS,ATAN,ATAN2,MIN,MOD
C     ..
C     .. Common blocks ..
      COMMON /EPP2/CC,UL,UR,VL,VR,UB,VB,IND
      COMMON /RNDOFF/EPSMIN
      COMMON /ZEE/Z
      COMMON /PIE/PI,TWOPI,HPI
      COMMON /TDATA/AA,TMID,BB,DTHDAA,DTHDBB,ADDD,MDTHZ
      COMMON /THET/THETU,THETV
C     ..
      DATA COL/'01','02','03','04','05','06','07','08','09','10','11',
     1         '12','13','14','15','16','17','18','19','20','21','22',
     2         '23','24','25','26','27','28','29','30','31','32'/
C
      CALL EXAMP()
C
      WRITE(*,*)
      WRITE(*,*) ' This program solves the boundary value problem '
      WRITE(*,*) '   defined by the differential equation '
      WRITE(*,*)
      WRITE(*,*) '   -(py'')'' + q*y = lambda*w*y  '
      WRITE(*,*)
      WRITE(*,*) ' together with appropriate boundary conditions. '
      WRITE(*,*)
      WRITE(*,*) ' HELP may be called at any point where the program   '
      WRITE(*,*) '    halts and displays (h?) by pressing "h <ENTER>". '
      WRITE(*,*) '    To RETURN from HELP, press "r <ENTER>".          '
      WRITE(*,*) '    To QUIT at any program halt, press "q <ENTER>".  '
      WRITE(*,*) ' WOULD YOU LIKE AN OVERVIEW OF HELP ? (Y/N) (h?) '
      READ(*,9010) CHANS
      READ(CHANS,9020) HQ
      IF (HQ.EQ.'q' .OR. HQ.EQ.'Q') GO TO 460
      IF (HQ.EQ.'h' .OR. HQ.EQ.'H' .OR. HQ.EQ.'y' .OR. HQ.EQ.'Y')
     1   CALL HELP(1)
      WRITE(*,*)
      WRITE(*,*) ' DO YOU REQUIRE INFORMATION ON THE RANGE OF '
      WRITE(*,*) '    BOUNDARY CONDITIONS AVAILABLE ? (Y/N) (h?) '
      READ(*,9010) CHANS
      READ(CHANS,9020) HQ
      IF (HQ.EQ.'q' .OR. HQ.EQ.'Q') GO TO 460
      IF (HQ.EQ.'h' .OR. HQ.EQ.'H' .OR. HQ.EQ.'y' .OR. HQ.EQ.'Y')
     1   CALL HELP(7)
   10 CONTINUE
         WRITE(*,*) ' DO YOU WANT A RECORD KEPT OF THE PROBLEMS '
         WRITE(*,*) '    AND RESULTS ? (Y/N) (h?) '
         READ(*,9010) CHANS
         READ(CHANS,9020) HQ
         IF (HQ.EQ.'q' .OR. HQ.EQ.'Q') GO TO 460
         IF (HQ.EQ.'h' .OR. HQ.EQ.'H') THEN
            CALL HELP(8)
            GO TO 10
            END IF
      READ(CHANS,9020) YN
      RITE = YN.EQ.'y' .OR. YN.EQ.'Y' .OR. YN.EQ.'yes' .OR. YN.EQ.'YES'
      IF (RITE) THEN
   20    CONTINUE
            WRITE(*,*) '  SPECIFY NAME OF THE OUTPUT RECORD FILE: (h?) '
            READ(*,9010) CHANS
            IF (CHANS.EQ.'q' .OR. CHANS.EQ.'Q') THEN
               GO TO 460
            ELSE IF (CHANS.EQ.'h' .OR. CHANS.EQ.'H') THEN
               CALL HELP(8)
               GO TO 20
            ELSE
               TAPE2 = CHANS
               WRITE(*,*) '  ENTER SOME HEADER LINE (<=70 CHARACTERS): '
               WRITE(*,*)
               READ(*,9030) CHTXT
               END IF
C
         OPEN(2,FILE=TAPE2,STATUS='NEW')
C
         WRITE(2,*) '                         ',TAPE2
         WRITE(2,*)
         WRITE(2,*) CHTXT
         WRITE(2,*)
         END IF
C
      OPEN(21,FILE='test.out')
C
C     DEFINITIONS OF SOME STRINGS.
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
C
   30 CONTINUE
      SKIPB = .FALSE.
C
   40 CONTINUE
      WRITE(*,*) ' ************************************************** '
      WRITE(*,*) ' * INDICATE THE KIND OF PROBLEM INTERVAL (a,b):   * '
      WRITE(*,*) ' *                                                * '
      WRITE(*,*) ' *  (CHECK THAT THE COEFFICIENTS p,q,w ARE WELL   * '
      WRITE(*,*) ' *  DEFINED THROUGHOUT THE INTERVAL open(a,b).)   * '
      WRITE(*,*) ' *                                                * '
      WRITE(*,*) ' *    (1) FINITE, (a,b)                           * '
      WRITE(*,*) ' *                                                * '
      WRITE(*,*) ' *    (2) SEMI-INFINITE, (a,+INFINITY)            * '
      WRITE(*,*) ' *                                                * '
      WRITE(*,*) ' *    (3) SEMI-INFINITE, (-INFINITY,b)            * '
      WRITE(*,*) ' *                                                * '
      WRITE(*,*) ' *    (4) DOUBLY INFINITE, (-INFINITY,+INFINITY)  * '
      WRITE(*,*) ' *                                                * '
      WRITE(*,*) ' * ENTER THE NUMBER OF YOUR CHOICE: (h?)          * '
      WRITE(*,*) ' ************************************************** '
      WRITE(*,*)
      READ(*,9010) CHANS
      READ(CHANS,9020) HQ
      IF (HQ.EQ.'q' .OR. HQ.EQ.'Q') GO TO 460
      IF (HQ.EQ.'h' .OR. HQ.EQ.'H') THEN
         CALL HELP(9)
         GO TO 40
         END IF
      READ(CHANS,'(I32)') INTAB
      IF (RITE) THEN
         IF (INTAB.EQ.1) WRITE(2,*) ' The interval is (a,b) .'
         IF (INTAB.EQ.2) WRITE(2,*) ' The interval is (a,+inf).'
         IF (INTAB.EQ.3) WRITE(2,*) ' The interval is (-inf,b). '
         IF (INTAB.EQ.4) WRITE(2,*) ' The interval is (-inf,+inf).'
         END IF
      IF (INTAB.LT.1 .OR. INTAB.GT.4) GO TO 40
C
   50 CONTINUE
      WRITE(*,*)
      P0ATA = -1.0D0
      QFATA =  1.0D0
      SINGATA = -1.0D0
      CIRCLA = -1.0D0
      OSCILA = -1.0D0
      REGA = .FALSE.
      WREGA = .FALSE.
      SINGA = .FALSE.
C
      IF (INTAB.EQ.1 .OR. INTAB.EQ.2) THEN
         WRITE(*,*)
         WRITE(*,*)
         WRITE(*,*)
         WRITE(*,*)
C
   60    CONTINUE
            WRITE(*,*)
            WRITE(*,*) ' ********************************************* '
            WRITE(*,*) ' * INPUT a: (h?)                             * '
            WRITE(*,*) ' ********************************************* '
            WRITE(*,*) ' a = '
            READ(*,9010) CHANS
            READ(CHANS,9020) HQ
            IF (HQ.EQ.'q' .OR. HQ.EQ.'Q') GO TO 460
            IF (HQ.EQ.'h' .OR. HQ.EQ.'H') THEN
               CALL HELP(10)
               GO TO 60
               END IF
         READ(CHANS,'(F32.0)') A
         IF (RITE) WRITE(2,*) ' a = ',A
         IF (INTAB.EQ.2) B = A + 1.0D0
         END IF
      WRITE(*,*)
      WRITE(*,*)
      WRITE(*,*)
      WRITE(*,*)
      WRITE(*,*)
C
   70 CONTINUE
         STR(1) = ' * IS THIS PROBLEM:                    '
         STR(2) = '                                      *'
         WRITE(*,9110) STAR
         WRITE(*,9110) STR
         WRITE(*,9110) BLNK
         STR(1) = ' *   (1) REGULAR AT a ?                '
         STR(2) = '                                      *'
         WRITE(*,9110) STR
         STR(1) = ' *       (I.E., THE FUNCTIONS p, q, & w'
         STR(2) = ' ARE BOUNDED CONTINUOUS NEAR a;       *'
         WRITE(*,9110) STR
         STR(1) = ' *        p & w ARE POSITIVE AT a.)    '
         STR(2) = '                                      *'
         WRITE(*,9110) STR
         WRITE(*,9110) BLNK
         STR(1) = ' *   (2) WEAKLY REGULAR AT a ?         '
         STR(2) = '                                      *'
         WRITE(*,9110) STR
         STR(1) = ' *       (I.E., THE FUNCTIONS 1/p, q, &'
         STR(2) = ' w ALL ARE FINITELY INTEGRABLE ON     *'
         WRITE(*,9110) STR
         STR(1) = ' *        SOME INTERVAL [a,a+e] FOR e >'
         STR(2) = ' 0; p & w ARE POSITIVE NEAR a.)       *'
         WRITE(*,9110) STR
         WRITE(*,9110) BLNK
         STR(1) = ' *   (3) LIMIT CIRCLE, NON-OSCILLATORY '
         STR(2) = 'AT a ?                                *'
         WRITE(*,9110) STR
         WRITE(*,9110) BLNK
         STR(1) = ' *   (4) LIMIT CIRCLE, OSCILLATORY AT a'
         STR(2) = ' ?                                    *'
         WRITE(*,9110) STR
         WRITE(*,9110) BLNK
         STR(1) = ' *   (5) LIMIT POINT AT a ?            '
         STR(2) = '                                      *'
         WRITE(*,9110) STR
         WRITE(*,9110) BLNK
         STR(1) = ' *   (6) NOT SPECIFIED (BUT NOT LIMIT C'
         STR(2) = 'IRCLE OSCILLATORY) WITH DEFAULT       *'
         WRITE(*,9110) STR
         STR(1) = ' *       BOUNDARY CONDITION AT a ?     '
         STR(2) = '                                      *'
         WRITE(*,9110) STR
         WRITE(*,9110) BLNK
         STR(1) = ' * ENTER THE NUMBER OF YOUR CHOICE: (h)'
         STR(2) = '?                                     *'
         WRITE(*,9110) STR
         WRITE(*,9110) STAR
         WRITE(*,*)
C
C     SPECIFY TYPE OF BOUNDARY CONDITION AT a.
C
         READ(*,9010) CHANS
         READ(CHANS,9020) HQ
         IF (HQ.EQ.'q' .OR. HQ.EQ.'Q') GO TO 460
         IF (HQ.EQ.'h' .OR. HQ.EQ.'H') THEN
            CALL HELP(4)
            GO TO 70
            END IF
         READ(CHANS,'(I32)') NANS
         IF (NANS.LT.1 .OR. NANS.GT.6) GO TO 70
C
C     SET CHARACTER STRING CHA ACCORDING TO BOUNDARY CONDITION AT a.
C
      IF (NANS.EQ.1) THEN
         REGA = .TRUE.
         CHA = CH1
         IF (RITE) WRITE(2,*) ' Endpoint a is Regular. '
      ELSE IF (NANS.EQ.2) THEN
         WREGA = .TRUE.
         CHA = CH2
         IF (RITE) WRITE(2,*) ' Endpoint a is Weakly Regular. '
      ELSE IF (NANS.EQ.3) THEN
         CIRCLA = 1.0D0
         CHA = CH3
         IF (RITE) WRITE(2,*) ' Endpoint a is Limit Circle,',
     1                        ' Non-Oscillatory. '
      ELSE IF (NANS.EQ.4) THEN
         CIRCLA = 1.0D0
         OSCILA = 1.0D0
         CHA = CH4
         IF (RITE) WRITE(2,*) ' Endpoint a is Limit Circle,',
     1                        ' Oscillatory. '
      ELSE IF (NANS.EQ.5) THEN
         CHA = CH5
         IF (RITE) WRITE(2,*) ' Endpoint a is Limit Point. '
      ELSE
         CHA = CH6
         IF (RITE) WRITE(2,*) ' Endpoint a is Singular, unspecified. '
         END IF
      WRITE(*,*)
      WRITE(*,*)
      WRITE(*,*)
      WRITE(*,*)
C
      IF (.NOT.REGA .AND. INTAB.LE.2) THEN
   80    CONTINUE
            WRITE(*,*)
            WRITE(*,*) ' ********************************************* '
            WRITE(*,*) ' * IS p = 0. AT a ? (Y/N) (h?)               * '
            WRITE(*,*) ' ********************************************* '
            WRITE(*,*)
            READ(*,9010) CHANS
            READ(CHANS,9020) HQ
            IF (HQ.EQ.'q' .OR. HQ.EQ.'Q') GO TO 460
            IF (HQ.EQ.'h' .OR. HQ.EQ.'H') THEN
               CALL HELP(11)
               GO TO 80
               END IF
            READ(CHANS,9020) YN
            YEH = YN.EQ.'y' .OR. YN.EQ.'Y'
            IF (.NOT.(YEH .OR. YN.EQ.'n' .OR. YN.EQ.'N')) GO TO 80
         IF (YEH) P0ATA = 1.0D0
         IF (RITE) THEN
            IF (YEH) THEN
               WRITE(2,*) ' p is zero at a. '
            ELSE
               WRITE(2,*) ' p is not zero at a. '
               END IF
            END IF
   90    CONTINUE
            WRITE(*,*)
            WRITE(*,*) ' ********************************************* '
            WRITE(*,*) ' * ARE THE COEFFICIENT FUNCTIONS q & w       * '
            WRITE(*,*) ' *                                           * '
            WRITE(*,*) ' *  FINITE AT THE ENDPOINT a ? (Y/N) (h?)    * '
            WRITE(*,*) ' ********************************************* '
            WRITE(*,*)
            READ(*,9010) CHANS
            READ(CHANS,9020) HQ
            IF (HQ.EQ.'q' .OR. HQ.EQ.'Q') GO TO 460
            IF (HQ.EQ.'h' .OR. HQ.EQ.'H') THEN
               CALL HELP(11)
               GO TO 90
               END IF
            READ(CHANS,9020) YN
            YEH = YN.EQ.'y' .OR. YN.EQ.'Y'
            IF (.NOT.(YEH .OR. YN.EQ.'n' .OR. YN.EQ.'N')) GO TO 90
         IF (.NOT.YEH) QFATA = -1.0D0
         IF (RITE) THEN
            IF (YEH) THEN
               WRITE(2,*) ' q and w are both finite at a. '
            ELSE
               WRITE(2,*) ' q and w are not both finite at a. '
               END IF
            END IF
         IF (P0ATA.LT.0.0D0 .AND. QFATA.GT.0.0D0) THEN
            IF (RITE)
     1         WRITE(2,*) 'This problem appears to be Regular at a.'
            WRITE(*,*)
            WRITE(*,*) ' ********************************************* '
            WRITE(*,*) ' * THIS PROBLEM APPEARS TO BE REGULAR AT a.  * '
            WRITE(*,*) ' ********************************************* '
            WRITE(*,*)
            CHA = CH1
            END IF
         END IF
      SINGA = .NOT.(REGA .OR. WREGA)
      IF (SINGA) SINGATA = 1.0D0
      IF (SKIPB) GO TO 150
C
  100 CONTINUE
      WRITE(*,*)
      P0ATB = -1.0D0
      QFATB =  1.0D0
      SINGATB = -1.0D0
      CIRCLB = -1.0D0
      OSCILB = -1.0D0
      REGB = .FALSE.
      WREGB = .FALSE.
      SINGB = .FALSE.
C
      IF (INTAB.EQ.1 .OR. INTAB.EQ.3) THEN
  110    CONTINUE
            WRITE(*,*)
            WRITE(*,*) ' ********************************************* '
            WRITE(*,*) ' * INPUT b: (h?)                             * '
            WRITE(*,*) ' ********************************************* '
            WRITE(*,*) ' b = '
            READ(*,9010) CHANS
            READ(CHANS,9020) HQ
            IF (HQ.EQ.'q' .OR. HQ.EQ.'Q') GO TO 460
            IF (HQ.EQ.'h' .OR. HQ.EQ.'H') THEN
               CALL HELP(10)
               GO TO 110
               END IF
         READ(CHANS,'(F32.0)') B
         IF (RITE) WRITE(2,*) ' b = ',B
         IF (INTAB.EQ.3) A = B - 1.0D0
         END IF
      WRITE(*,*)
      WRITE(*,*)
      WRITE(*,*)
      WRITE(*,*)
      WRITE(*,*)
C
  120 CONTINUE
         STR(1) = ' * IS THIS PROBLEM:                    '
         STR(2) = '                                      *'
         WRITE(*,9110) STAR
         WRITE(*,9110) STR
         WRITE(*,9110) BLNK
         STR(1) = ' *   (1) REGULAR AT b ?                '
         STR(2) = '                                      *'
         WRITE(*,9110) STR
         STR(1) = ' *       (I.E., THE FUNCTIONS p, q, & w'
         STR(2) = ' ARE BOUNDED CONTINUOUS NEAR b;       *'
         WRITE(*,9110) STR
         STR(1) = ' *        p & w ARE POSITIVE AT b.)    '
         STR(2) = '                                      *'
         WRITE(*,9110) STR
         WRITE(*,9110) BLNK
         STR(1) = ' *   (2) WEAKLY REGULAR AT b ?         '
         STR(2) = '                                      *'
         WRITE(*,9110) STR
         STR(1) = ' *       (I.E., THE FUNCTIONS 1/p, q, &'
         STR(2) = ' w ALL ARE FINITELY INTEGRABLE ON     *'
         WRITE(*,9110) STR
         STR(1) = ' *        SOME INTERVAL [b-e,b] FOR e >'
         STR(2) = ' 0; p & w ARE POSITIVE NEAR b.)       *'
         WRITE(*,9110) STR
         WRITE(*,9110) BLNK
         STR(1) = ' *   (3) LIMIT CIRCLE, NON-OSCILLATORY '
         STR(2) = 'AT b ?                                *'
         WRITE(*,9110) STR
         WRITE(*,9110) BLNK
         STR(1) = ' *   (4) LIMIT CIRCLE, OSCILLATORY AT b'
         STR(2) = ' ?                                    *'
         WRITE(*,9110) STR
         WRITE(*,9110) BLNK
         STR(1) = ' *   (5) LIMIT POINT AT b ?            '
         STR(2) = '                                      *'
         WRITE(*,9110) STR
         WRITE(*,9110) BLNK
         STR(1) = ' *   (6) NOT SPECIFIED (BUT NOT LIMIT C'
         STR(2) = 'IRCLE OSCILLATORY) WITH DEFAULT       *'
         WRITE(*,9110) STR
         STR(1) = ' *       BOUNDARY CONDITION AT b ?     '
         STR(2) = '                                      *'
         WRITE(*,9110) STR
         WRITE(*,9110) BLNK
         STR(1) = ' * ENTER THE NUMBER OF YOUR CHOICE: (h?'
         STR(2) = ')                                     *'
         WRITE(*,9110) STR
         WRITE(*,9110) STAR
         WRITE(*,*)
C
C     SPECIFY TYPE OF BOUNDARY CONDITION AT b.
C
         READ(*,9010) CHANS
         READ(CHANS,9020) HQ
         IF (HQ.EQ.'q' .OR. HQ.EQ.'Q') GO TO 460
         IF (HQ.EQ.'h' .OR. HQ.EQ.'H') THEN
            CALL HELP(4)
            GO TO 120
            END IF
         READ(CHANS,'(I32)') NANS
         IF (NANS.LT.1 .OR. NANS.GT.6) GO TO 120
C
C     SET CHARACTER STRING CHB ACCORDING TO BOUNDARY CONDITION AT b.
C
      WRITE(*,*)
      IF (NANS.EQ.1) THEN
         REGB = .TRUE.
         CHB = CH1
         IF (RITE) WRITE(2,*) ' Endpoint b is Regular. '
      ELSE IF (NANS.EQ.2) THEN
         WREGB = .TRUE.
         CHB = CH2
         IF (RITE) WRITE(2,*) ' Endpoint b is Weakly Regular. '
      ELSE IF (NANS.EQ.3) THEN
         CIRCLB = 1.0D0
         CHB = CH3
         IF (RITE) WRITE(2,*) ' Endpoint b is Limit Circle,',
     1                        ' Non-Oscillatory. '
      ELSE IF (NANS.EQ.4) THEN
         CIRCLB = 1.0D0
         OSCILB = 1.0D0
         CHB = CH4
         IF (RITE) WRITE(2,*) ' Endpoint b is Limit Circle,',
     1                        ' Oscillatory. '
      ELSE IF (NANS.EQ.5) THEN
         CHB = CH5
         IF (RITE) WRITE(2,*) ' Endpoint b is Limit Point. '
      ELSE
         CHB = CH6
         IF (RITE) WRITE(2,*) ' Endpoint b is Singular, unspecified. '
         END IF
      WRITE(*,*)
      WRITE(*,*)
      WRITE(*,*)
      WRITE(*,*)
C
      IF (.NOT.REGB .AND. (INTAB.EQ.1 .OR. INTAB.EQ.3)) THEN
  130    CONTINUE
            WRITE(*,*)
            WRITE(*,*) ' ********************************************* '
            WRITE(*,*) ' * IS p = 0. AT b ? (Y/N) (h?)               * '
            WRITE(*,*) ' ********************************************* '
            WRITE(*,*)
            READ(*,9010) CHANS
            READ(CHANS,9020) HQ
            IF (HQ.EQ.'q' .OR. HQ.EQ.'Q') GO TO 460
            IF (HQ.EQ.'h' .OR. HQ.EQ.'H') THEN
               CALL HELP(11)
               GO TO 130
               END IF
            READ(CHANS,9020) YN
            YEH = YN.EQ.'y' .OR. YN.EQ.'Y'
            IF (.NOT.(YEH .OR. YN.EQ.'n' .OR. YN.EQ.'N')) GO TO 130
         IF (YEH) P0ATB = 1.0D0
         IF (RITE) THEN
            IF (YEH) THEN
               WRITE(2,*) ' p is zero at b. '
            ELSE
               WRITE(2,*) ' p is not zero at b. '
               END IF
            END IF
  140    CONTINUE
            WRITE(*,*)
            WRITE(*,*) ' ********************************************* '
            WRITE(*,*) ' * ARE THE COEFFICIENT FUNCTIONS q & w       * '
            WRITE(*,*) ' *                                           * '
            WRITE(*,*) ' *   FINITE AT THE ENDPOINT b ? (Y/N) (h?)   * '
            WRITE(*,*) ' ********************************************* '
            WRITE(*,*)
            READ(*,9010) CHANS
            READ(CHANS,9020) HQ
            IF (HQ.EQ.'q' .OR. HQ.EQ.'Q') GO TO 460
            IF (HQ.EQ.'h' .OR. HQ.EQ.'H') THEN
               CALL HELP(11)
               GO TO 140
               END IF
            READ(CHANS,9020) YN
            YEH = YN.EQ.'y' .OR. YN.EQ.'Y'
            IF (.NOT.(YEH .OR. YN.EQ.'n' .OR. YN.EQ.'N')) GO TO 140
         IF (.NOT.YEH) QFATB = -1.0D0
         IF (RITE) THEN
            IF (YEH) THEN
               WRITE(2,*) ' q and w are both finite at b. '
            ELSE
               WRITE(2,*) ' q and w are not both finite at b. '
               END IF
            END IF
         IF (P0ATB.LT.0.0D0 .AND. QFATB.GT.0.0D0) THEN
            IF (RITE)
     1         WRITE(2,*) 'This problem appears to be Regular at b.'
            WRITE(*,*)
            WRITE(*,*) ' ********************************************* '
            WRITE(*,*) ' * THIS PROBLEM APPEARS TO BE REGULAR AT b.  * '
            WRITE(*,*) ' ********************************************* '
            CHB = CH1
            END IF
         END IF
      SINGB = .NOT.(REGB .OR. WREGB)
      IF (SINGB) SINGATB = 1.0D0
  150 CONTINUE
         WRITE(*,*)
         WRITE(*,*)
         WRITE(*,*)
         WRITE(*,*)
         WRITE(*,*)
         WRITE(*,*) ' ************************************************ '
         WRITE(*,*) ' * THIS PROBLEM IS ON THE INTERVAL              * '
         WRITE(*,*) ' *                                              * '
         IF (INTAB.EQ.1) THEN
            WRITE(*,9070) A,B
            WRITE(*,*) ' *                            ',FILLB
            WRITE(*,*) ' * ENDPOINT a IS  ',CHA
            WRITE(*,*) ' *                            ',FILLB
            IF (P0ATA.GT.0.0D0) THEN
               WRITE(*,*) ' * p IS ZERO AT a             ',FILLB
               WRITE(*,*) ' *                            ',FILLB
               END IF
            IF (QFATA.LT.0.0D0) THEN
               WRITE(*,*) ' * q & w ARE NOT BOUNDED AT a ',FILLB
               WRITE(*,*) ' *                            ',FILLB
               END IF
            WRITE(*,*) ' * ENDPOINT b IS  ',CHB
            WRITE(*,*) ' *                            ',FILLB
            IF (P0ATB.GT.0.0D0) THEN
               WRITE(*,*) ' * p IS ZERO AT b             ',FILLB
               WRITE(*,*) ' *                            ',FILLB
               END IF
            IF (QFATB.LT.0.0D0) THEN
               WRITE(*,*) ' * q & w ARE NOT BOUNDED AT b ',FILLB
               WRITE(*,*) ' *                            ',FILLB
               END IF
         ELSE IF (INTAB.EQ.2) THEN
            WRITE(*,9080) A,INFP
            WRITE(*,*) ' *                            ',FILLB
            WRITE(*,*) ' * ENDPOINT a IS  ',CHA
            WRITE(*,*) ' *                            ',FILLB
            IF (P0ATA.GT.0.0D0) THEN
               WRITE(*,*) ' * p IS ZERO AT a             ',FILLB
               WRITE(*,*) ' *                            ',FILLB
               END IF
            IF (QFATA.LT.0.0D0) THEN
               WRITE(*,*) ' * q & w ARE NOT BOUNDED AT a ',FILLB
               WRITE(*,*) ' *                            ',FILLB
               END IF
            WRITE(*,*) ' * ENDPT +INF IS  ',CHB
            WRITE(*,*) ' *                            ',FILLB
         ELSE IF (INTAB.EQ.3) THEN
            WRITE(*,9090) INFM,B
            WRITE(*,*) ' *                            ',FILLB
            WRITE(*,*) ' * ENDPT -INF IS  ',CHA
            WRITE(*,*) ' *                            ',FILLB
            WRITE(*,*) ' * ENDPOINT b IS  ',CHB
            WRITE(*,*) ' *                            ',FILLB
            IF (P0ATB.GT.0.0D0) THEN
               WRITE(*,*) ' * p IS ZERO AT b             ',FILLB
               WRITE(*,*) ' *                            ',FILLB
               END IF
            IF (QFATB.LT.0.0D0) THEN
               WRITE(*,*) ' * q & w ARE NOT BOUNDED AT b ',FILLB
               WRITE(*,*) ' *                            ',FILLB
               END IF
         ELSE
            WRITE(*,9100) INFM,INFP
            WRITE(*,*) ' *                            ',FILLB
            WRITE(*,*) ' * ENDPT -INF IS  ',CHA
            WRITE(*,*) ' *                            ',FILLB
            WRITE(*,*) ' * ENDPT +INF IS  ',CHB
            WRITE(*,*) ' *                            ',FILLB
            END IF
         WRITE(*,*) ' * IS THIS THE PROBLEM YOU WANT ? (Y/N) (h?)    * '
         WRITE(*,*) ' ************************************************ '
         WRITE(*,*)
         READ(*,9010) CHANS
         READ(CHANS,9020) HQ
         IF (HQ.EQ.'q' .OR. HQ.EQ.'Q') GO TO 460
         IF (HQ.EQ.'h' .OR. HQ.EQ.'H') THEN
            CALL HELP(3)
            GO TO 150
            END IF
      READ(CHANS,9020) YN
      IF (YN.NE.'y' .AND. YN.NE.'Y') THEN
  160    CONTINUE
            WRITE(*,*)
            WRITE(*,*) ' ********************************************* '
            WRITE(*,*) ' * DO YOU WANT TO RE-DO                      * '
            WRITE(*,*) ' *                                           * '
            WRITE(*,*) ' *    (1)  ENDPOINT a                        * '
            WRITE(*,*) ' *                                           * '
            WRITE(*,*) ' *    (2)  ENDPOINT b                        * '
            WRITE(*,*) ' *                                           * '
            WRITE(*,*) ' *    (3)  BOTH ENDPOINTS a AND b ?          * '
            WRITE(*,*) ' *                                           * '
            WRITE(*,*) ' * ENTER THE NUMBER OF YOUR CHOICE: (h?)     * '
            WRITE(*,*) ' ********************************************* '
            WRITE(*,*)
            READ(*,9010) CHANS
            READ(CHANS,9020) HQ
            IF (HQ.EQ.'q' .OR. HQ.EQ.'Q') GO TO 460
            IF (HQ.EQ.'h' .OR. HQ.EQ.'H') THEN
               CALL HELP(3)
               GO TO 160
               END IF
            READ(CHANS,'(I32)') NANS
            IF (NANS.LT.1 .OR. NANS.GT.3) GO TO 160
         IF (NANS.EQ.1) THEN
            SKIPB = .TRUE.
            IF (RITE) WRITE(2,*) ' Redo endpoint a. '
            GO TO 50
         ELSE IF (NANS.EQ.2) THEN
            IF (RITE) WRITE(2,*) ' Redo endpoint b. '
            GO TO 100
         ELSE
            IF (RITE) WRITE(2,*) ' Redo both endpoints a and b. '
            GO TO 30
            END IF
         END IF
      WRITE(*,*)
      WRITE(*,*)
      WRITE(*,*)
      WRITE(*,*)
      WRITE(*,*)
C
C     AT THIS POINT THE DIFFERENTIAL EQUATION AND THE INTERVAL OF
C       INTEREST HAVE BEEN DEFINED AND CHARACTERIZED.
C
  170 CONTINUE
         IF (RITE) WRITE(2,*) '----------------------------------------'
         IF (RITE) WRITE(2,*)
         WRITE(*,*)
         WRITE(*,*) ' *************************************************'
         WRITE(*,*) ' * DO YOU WANT TO COMPUTE                        *'
         WRITE(*,*) ' *                                               *'
         WRITE(*,*) ' *   (1) AN EIGENVALUE, OR SERIES OF EIGENVALUES *'
         WRITE(*,*) ' *                                               *'
         WRITE(*,*) ' *   (2) SOLUTION TO AN INITIAL VALUE PROBLEM ?  *'
         WRITE(*,*) ' *                                               *'
         WRITE(*,*) ' * ENTER THE NUMBER OF YOUR CHOICE: (h?)         *'
         WRITE(*,*) ' *************************************************'
         WRITE(*,*)
         READ(*,9010) CHANS
         READ(CHANS,9020) HQ
         IF (HQ.EQ.'q' .OR. HQ.EQ.'Q') GO TO 460
         IF (HQ.EQ.'h' .OR. HQ.EQ.'H') THEN
            CALL HELP(12)
            GO TO 170
            END IF
         READ(CHANS,'(I32)') NANS
         IF (NANS.LT.1 .OR. NANS.GT.2) GO TO 170
      EIGV = .TRUE.
      IF (NANS.EQ.2) THEN
         EIGV = .FALSE.
         GO TO 350
         END IF
C
C     THIS IS THE ENTRY POINT FOR COMPUTING A NEW EIGENVALUE
C     OR A SERIES OF EIGENVALUES.
C
  180 CONTINUE
         IF ((REGA .OR. WREGA) .AND. (REGB .OR. WREGB)) THEN
            WRITE(*,*) ' ******************************************** '
            WRITE(*,*) ' * IS THE BOUNDARY CONDITION PERIODIC ?     * '
            WRITE(*,*) ' *                                          * '
            WRITE(*,*) ' *   (I.E., y(b) = c*y(a)                   * '
            WRITE(*,*) ' *    &   p(b)*y''(b) = (1/c)*p(a)*y''(a) )   *'
            WRITE(*,*) ' *                                          * '
            WRITE(*,*) ' * ANSWER (Y/N): (h?)                       * '
            WRITE(*,*) ' ******************************************** '
            WRITE(*,*)
            READ(*,9010) CHANS
            READ(CHANS,9020) HQ
            IF (HQ.EQ.'q' .OR. HQ.EQ.'Q') GO TO 460
            IF (HQ.EQ.'h' .OR. HQ.EQ.'H') THEN
               CALL HELP(7)
               GO TO 180
               END IF
            READ(CHANS,9020) YN
            PERIOD = YN.EQ.'y' .OR. YN.EQ.'Y'
            END IF
      WRITE(*,*)
      WRITE(*,*)
      WRITE(*,*)
      WRITE(*,*)
      WRITE(*,*)
      IF (PERIOD) GO TO 310
         IF (SINGATA.LT.0.0D0) THEN
  190       CONTINUE
               WRITE(*,*) ' ****************************************** '
               WRITE(*,*) ' * IS THE BOUNDARY CONDITION AT a         * '
               WRITE(*,*) ' *                                        * '
               WRITE(*,*) ' *    (1) THE DIRICHLET CONDITION         * '
               WRITE(*,*) ' *        (I.E., y(a) = 0.0)              * '
               WRITE(*,*) ' *                                        * '
               WRITE(*,*) ' *    (2) THE NEUMANN CONDITION           * '
               WRITE(*,*) ' *        (I.E., y''(a) = 0.0)             *'
               WRITE(*,*) ' *                                        * '
               WRITE(*,*) ' *    (3) A MORE GENERAL LINEAR           * '
               WRITE(*,*) ' *        BOUNDARY CONDITION              * '
               WRITE(*,*) ' *        A1*[y(a)] + A2*[py''](a) = 0 ?   *'
               WRITE(*,*) ' *                                        * '
               WRITE(*,*) ' * ENTER THE NUMBER OF YOUR CHOICE: (h?)  * '
               WRITE(*,*) ' ****************************************** '
               WRITE(*,*)
               READ(*,9010) CHANS
               READ(CHANS,9020) HQ
               IF (HQ.EQ.'q' .OR. HQ.EQ.'Q') GO TO 460
               IF (HQ.EQ.'h' .OR. HQ.EQ.'H') THEN
                  CALL HELP(7)
                  GO TO 190
                  END IF
               READ(CHANS,'(I32)') NANS
               IF (NANS.LT.1 .OR. NANS.GT.3) GO TO 190
            IF (NANS.EQ.1) THEN
               A1 = 1.0D0
               A2 = 0.0D0
               IF (RITE) WRITE(2,*) ' Dirichlet B.C. at a. '
            ELSE IF (NANS.EQ.2) THEN
               A1 = 0.0D0
               A2 = 1.0D0
               IF (RITE) WRITE(2,*) ' Neumann B.C. at a. '
            ELSE
  200          CONTINUE
                  WRITE(*,*) ' *************************************** '
                  WRITE(*,*) ' * CHOOSE A1,A2: (h?)                  * '
                  WRITE(*,*) ' *************************************** '
                  WRITE(*,*)
                  WRITE(*,*) ' A1,A2 = '
                  READ(*,9010) CHANS
                  READ(CHANS,9020) HQ
                  IF (HQ.EQ.'q' .OR. HQ.EQ.'Q') GO TO 460
                  IF (HQ.EQ.'h' .OR. HQ.EQ.'H') THEN
                     CALL HELP(7)
                     GO TO 200
                     END IF
               CALL LSTDIR(CHANS,I,ICOL)
               I1 = ICOL(1) - 1
               FMT = FMT2(I1)
               READ(CHANS,FMT) A1,A2
               IF (RITE) WRITE(2,*) ' A1,A2 = ',a1,a2
               END IF
         ELSE IF (CIRCLA.GT.0.0D0) THEN
  210       CONTINUE
               IF (RITE) THEN
                  WRITE(2,*) ' The B.C. at a is '
                  WRITE(2,*) ' A1*[y,u](a) + A2*[y,v](a) = 0. '
                  END IF
               WRITE(*,*) ' ****************************************** '
               WRITE(*,*) ' * THE BOUNDARY CONDITION AT a IS         * '
               WRITE(*,*) ' *                                        * '
               WRITE(*,*) ' *    A1*[y,u](a) + A2*[y,v](a) = 0,      * '
               WRITE(*,*) ' *                                        * '
               WRITE(*,*) ' *  WHERE THE CONSTANTS A1 AND A2         * '
               WRITE(*,*) ' *    MAY BE CHOSEN ARBITRARILY.          * '
               WRITE(*,*) ' *                                        * '
               WRITE(*,*) ' * CHOOSE A1,A2: (h?)                     * '
               WRITE(*,*) ' ****************************************** '
               WRITE(*,*)
               WRITE(*,*) ' A1,A2 = '
               READ(*,9010) CHANS
               READ(CHANS,9020) HQ
               IF (HQ.EQ.'q' .OR. HQ.EQ.'Q') GO TO 460
               IF (HQ.EQ.'h' .OR. HQ.EQ.'H') THEN
                  CALL HELP(7)
                  GO TO 210
                  END IF
            CALL LSTDIR(CHANS,I,ICOL)
            I1 = ICOL(1) - 1
            FMT = FMT2(I1)
            READ(CHANS,FMT) A1,A2
            IF (RITE) WRITE(2,*) ' A1,A2 = ',A1,A2
            END IF
         IF (SINGATB.LT.0.0D0) THEN
  220       CONTINUE
               WRITE(*,*) ' ****************************************** '
               WRITE(*,*) ' * IS THE BOUNDARY CONDITION AT b         * '
               WRITE(*,*) ' *                                        * '
               WRITE(*,*) ' *    (1) THE DIRICHLET CONDITION         * '
               WRITE(*,*) ' *        (I.E., y(b) = 0.0)              * '
               WRITE(*,*) ' *                                        * '
               WRITE(*,*) ' *    (2) THE NEUMANN CONDITION           * '
               WRITE(*,*) ' *        (I.E., y''(b) = 0.0)             *'
               WRITE(*,*) ' *                                        * '
               WRITE(*,*) ' *    (3) A MORE GENERAL LINEAR           * '
               WRITE(*,*) ' *        BOUNDARY CONDITION              * '
               WRITE(*,*) ' *        B1*[y(b)] + B2*[py''](b) = 0 ?   *'
               WRITE(*,*) ' *                                        * '
               WRITE(*,*) ' * ENTER THE NUMBER OF YOUR CHOICE: (h?)  * '
               WRITE(*,*) ' ****************************************** '
               WRITE(*,*)
               READ(*,9010) CHANS
               READ(CHANS,9020) HQ
               IF (HQ.EQ.'q' .OR. HQ.EQ.'Q') GO TO 460
               IF (HQ.EQ.'h' .OR. HQ.EQ.'H') THEN
                  CALL HELP(7)
                  GO TO 220
                  END IF
               READ(CHANS,'(I32)') NANS
               IF (NANS.LT.1 .OR. NANS.GT.3) GO TO 220
            IF (NANS.EQ.1) THEN
               B1 = 1.0D0
               B2 = 0.0D0
               IF (RITE) WRITE(2,*) ' Dirichlet B.C. at b. '
            ELSE IF (NANS.EQ.2) THEN
               B1 = 0.0D0
               B2 = 1.0D0
               IF (RITE) WRITE(2,*) ' Neumann B.C. at b. '
            ELSE
  230          CONTINUE
                  WRITE(*,*) ' *************************************** '
                  WRITE(*,*) ' * CHOOSE B1,B2: (h?)                  * '
                  WRITE(*,*) ' *************************************** '
                  WRITE(*,*)
                  WRITE(*,*) ' B1,B2 = '
                  READ(*,9010) CHANS
                  READ(CHANS,9020) HQ
                  IF (HQ.EQ.'q' .OR. HQ.EQ.'Q') GO TO 460
                  IF (HQ.EQ.'h' .OR. HQ.EQ.'H') THEN
                     CALL HELP(7)
                     GO TO 230
                     END IF
               CALL LSTDIR(CHANS,I,ICOL)
               I1 = ICOL(1) - 1
               FMT = FMT2(I1)
               READ(CHANS,FMT) B1,B2
               IF (RITE) WRITE(2,*) ' B1,B2 = ',b1,b2
               END IF
         ELSE IF (CIRCLB.GT.0.0D0) THEN
  240       CONTINUE
               IF (RITE) THEN
                  WRITE(2,*) ' The B.C. at b is '
                  WRITE(2,*) ' B1*[y,u](b) + B2*[y,v](b) = 0. '
                  END IF
               WRITE(*,*) ' ****************************************** '
               WRITE(*,*) ' * THE BOUNDARY CONDITION AT b IS         * '
               WRITE(*,*) ' *                                        * '
               WRITE(*,*) ' *    B1*[y,u](b) + B2*[y,v](b) = 0,      * '
               WRITE(*,*) ' *                                        * '
               WRITE(*,*) ' *  WHERE THE CONSTANTS B1 AND B2         * '
               WRITE(*,*) ' *    MAY BE CHOSEN ARBITRARILY.          * '
               WRITE(*,*) ' *                                        * '
               WRITE(*,*) ' * CHOOSE B1,B2: (h?)                     * '
               WRITE(*,*) ' ****************************************** '
               WRITE(*,*)
               WRITE(*,*) ' B1,B2 = '
               READ(*,9010) CHANS
               READ(CHANS,9020) HQ
               IF (HQ.EQ.'q' .OR. HQ.EQ.'Q') GO TO 460
               IF (HQ.EQ.'h' .OR. HQ.EQ.'H') THEN
                  CALL HELP(7)
                  GO TO 240
                  END IF
            CALL LSTDIR(CHANS,I,ICOL)
            I1 = ICOL(1) - 1
            FMT = FMT2(I1)
            READ(CHANS,FMT) B1,B2
            IF (RITE) WRITE(2,*) ' B1,B2 = ',B1,B2
            END IF
         WRITE(*,*)
         WRITE(*,*)
         WRITE(*,*)
         WRITE(*,*)
         WRITE(*,*)
C
C     THIS IS THE ENTRY POINT FOR COMPUTING AN EIGENVALUE
C     OR A SERIES OF EIGENVALUES IN THE NON-PERIODIC CASE.
C
  250    CONTINUE
            IF (RITE)
     1         WRITE(2,*) ' *******************************'//FILLA
            WRITE(*,*)
            WRITE(*,*) ' ********************************************* '
            WRITE(*,*) ' * DO YOU WANT TO COMPUTE                    * '
            WRITE(*,*) ' *                                           * '
            WRITE(*,*) ' *    (1) A SINGLE EIGENVALUE                * '
            WRITE(*,*) ' *                                           * '
            WRITE(*,*) ' *    (2) A SERIES OF EIGENVALUES ?          * '
            WRITE(*,*) ' *                                           * '
            WRITE(*,*) ' * ENTER THE NUMBER OF YOUR CHOICE: (h?)     * '
            WRITE(*,*) ' ********************************************* '
            WRITE(*,*)
            READ(*,9010) CHANS
            READ(CHANS,9020) HQ
            IF (HQ.EQ.'q' .OR. HQ.EQ.'Q') GO TO 460
            IF (HQ.EQ.'h' .OR. HQ.EQ.'H') THEN
               CALL HELP(13)
               GO TO 250
               END IF
            READ(CHANS,'(I32)') NANS
            IF (NANS.LT.1 .OR. NANS.GT.2) GO TO 250
         IF (NANS.EQ.1) THEN
  260       CONTINUE
               WRITE(*,*) ' ****************************************** '
               WRITE(*,*) ' * INPUT NUMEIG, EIG, TOL: (h?)           * '
               WRITE(*,*) ' ****************************************** '
               WRITE(*,*)
               WRITE(*,*) ' NUMEIG,EIG,TOL = '
               READ(*,9010) CHANS
               READ(CHANS,9020) HQ
               IF (HQ.EQ.'q' .OR. HQ.EQ.'Q') GO TO 460
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
               FMT = '(I'//COL(I1)//',1X,F'//COL(I2)
     1               //'.0,1X,F'//COL(30-I1-I2)//'.0)'
               READ(CHANS,FMT) NUMEIG,EIG,TOL
            ELSE IF (I.EQ.2) THEN
               FMT = '(I'//COL(I1)//',1X,F'//COL(31-I1)//'.0)'
               READ(CHANS,FMT) NUMEIG,EIG
            ELSE
               FMT = '(I'//COL(I1)//')'
               READ(CHANS,FMT) NUMEIG
               END IF
            IF (RITE) WRITE(2,*) ' NUMEIG,EIG,TOL = ',NUMEIG,EIG,TOL
            WRITE(*,*)
            WRITE(*,*)
            WRITE(*,*)
            WRITE(*,*)
            WRITE(*,*)
C
            I2 = NUMEIG
            CALL SLEIGN2(A,B,INTAB,P0ATA,QFATA,P0ATB,QFATB,A1,A2,
     1                   B1,B2,I2,EIG,TOL,IFLAG,0,SLFN,
     2                   SINGATA,SINGATB,CIRCLA,CIRCLB,OSCILA,OSCILB)
            IFLAG = MIN(IFLAG,4)
            WRITE(*,*) ' *******************************'//FILLA
            IF (IFLAG.LE.2) THEN
               WRITE(*,*) '  * NUMEIG = ',NUMEIG,'   EIG =',EIG
               WRITE(*,9050) TOL,IFLAG
               IF (RITE) THEN
                  WRITE(2,*) '  * NUMEIG = ',NUMEIG,'   EIG =',EIG
                  WRITE(2,9050) TOL,IFLAG
                  WRITE(2,*) ' *******************************'//FILLA
                  END IF
            ELSE
               WRITE(*,9040) NUMEIG,IFLAG
               IF (RITE) WRITE(2,9040) NUMEIG,IFLAG
               IF (IFLAG.EQ.3) THEN
                  WRITE(*,9230) I2
                  WRITE(*,9250) EIG
                  IF (RITE) THEN
                     WRITE(2,9230) I2
                     WRITE(2,9250) EIG
                     END IF
                  END IF
               WRITE(*,*) ' *******************************'//FILLA
               END IF
  270       CONTINUE
               IF (IFLAG.LE.2) THEN
                  WRITE(*,*) ' *                              '//FILLB
                  WRITE(*,*) ' * PLOT EIGENFUNCTION ? (Y/N) (h?)'//
     1                       '                *'
                  WRITE(*,*) ' *******************************'//FILLA
                  WRITE(*,*)
                  READ(*,9010) CHANS
                  READ(CHANS,9020) HQ
                  IF (HQ.EQ.'q' .OR. HQ.EQ.'Q') GO TO 460
                  IF (HQ.EQ.'h' .OR. HQ.EQ.'H') THEN
                     CALL HELP(16)
                     GO TO 270
                     END IF
                  READ(CHANS,9020) YN
                  PEIGF = YN.EQ.'y' .OR. YN.EQ.'Y'
                  IF (PEIGF) GO TO 420
                  END IF
         ELSE
  280       CONTINUE
               WRITE(*,*) ' ****************************************** '
               WRITE(*,*) ' * INPUT NUMEIG1, NUMEIG2, TOL (h?)       * '
               WRITE(*,*) ' ****************************************** '
               WRITE(*,*)
               WRITE(*,*) ' NUMEIG1,NUMEIG2,TOL = '
               READ(*,9010) CHANS
               READ(CHANS,9020) HQ
               IF (HQ.EQ.'q' .OR. HQ.EQ.'Q') GO TO 460
               IF (HQ.EQ.'h' .OR. HQ.EQ.'H') THEN
                  CALL HELP(14)
                  GO TO 280
                  END IF
            TOLL = 0.0D0
            CALL LSTDIR(CHANS,I,ICOL)
            I1 = ICOL(1) - 1
            IF (I.EQ.3) THEN
               I2 = ICOL(2) - ICOL(1) - 1
               FMT = '(I'//COL(I1)//',1X,I'//COL(I2)
     1               //',1X,F'//COL(30-I1-I2)//'.0)'
               READ(CHANS,FMT) NUMEIG1,NUMEIG2,TOLL
            ELSE
               FMT = '(I'//COL(I1)//',1X,I'//COL(31-I1)//')'
               READ(CHANS,FMT) NUMEIG1,NUMEIG2
            END IF
            IF (RITE) WRITE(2,*) ' NUMEIG1,NUMEIG2,TOL = ',
     1                             NUMEIG1,NUMEIG2,TOLL
C
            I = 0
            DO 290 NUMEIG = NUMEIG1,NUMEIG2
               TOL = TOLL
               EIG = 0.0D0
               I2 = NUMEIG
               CALL SLEIGN2(A,B,INTAB,P0ATA,QFATA,P0ATB,QFATB,A1,A2,
     1                      B1,B2,I2,EIG,TOL,IFLAG,0,SLFN,
     2                      SINGATA,SINGATB,CIRCLA,CIRCLB,OSCILA,OSCILB)
               IFLAG = MIN(IFLAG,4)
               IF (IFLAG.LE.2) THEN
                  WRITE(*,*) '  * NUMEIG = ',NUMEIG,'   EIG =',EIG
                  WRITE(*,9050) TOL,IFLAG
                  IF (RITE) THEN
                     WRITE(2,*) '  * NUMEIG = ',NUMEIG,'   EIG =',EIG
                     WRITE(2,9050) TOL,IFLAG
                     END IF
               ELSE
                  WRITE(*,9040) NUMEIG,IFLAG
                  IF (RITE) WRITE(2,9040) NUMEIG,IFLAG
                  IF (IFLAG.EQ.3 .AND. NUMEIG.EQ.NUMEIG2) THEN
                     WRITE(*,9230) I2
                     WRITE(*,9250) EIG
                     IF (RITE) THEN
                        WRITE(2,9230) I2
                        WRITE(2,9250) EIG
                        END IF
                     END IF
                  END IF
               I = I + 1
               EES(I) = EIG
               TTS(I) = TOL
               IIS(I) = IFLAG
  290          CONTINUE
            IF (RITE)
     1         WRITE(2,*) ' *******************************'//FILLA
            WRITE(*,*)
            I = 0
            DO 300 NUMEIG = NUMEIG1,NUMEIG2
               I = I + 1
               WRITE(*,9060) I+NUMEIG1-1,EES(I),TTS(I),IIS(I)
  300          CONTINUE
            END IF
         WRITE(*,*)
         GO TO 430
C
C     THIS IS THE ENTRY POINT FOR COMPUTING AN EIGENVALUE
C     IN A PERIODIC TYPE PROBLEM.
C
  310    CONTINUE
            WRITE(*,*) ' ******************************************** '
            WRITE(*,*) ' * IS THIS PROBLEM:                         * '
            WRITE(*,*) ' *                                          * '
            WRITE(*,*) ' *   (1) PERIODIC ?                         * '
            WRITE(*,*) ' *         (I.E., y(b) = y(a)               * '
            WRITE(*,*) ' *          &   p(b)*y''(b) = p(a)*y''(a) )   *'
            WRITE(*,*) ' *                                          * '
            WRITE(*,*) ' *   (2) SEMI-PERIODIC ?                    * '
            WRITE(*,*) ' *         (I.E., y(b) = -y(a)              * '
            WRITE(*,*) ' *          &   p(b)*y''(b) = -p(a)*y''(a) )  *'
            WRITE(*,*) ' *                                          * '
            WRITE(*,*) ' *   (3) GENERAL PERIODIC TYPE ?            * '
            WRITE(*,*) ' *         (I.E., y(b) = c*y(a)             * '
            WRITE(*,*) ' *          &   p(b)*y''(b) = p(a)*y''(a)/c   *'
            WRITE(*,*) ' *            for some number c .NE. 0. )   * '
            WRITE(*,*) ' *                                          * '
            WRITE(*,*) ' * ENTER THE NUMBER OF YOUR CHOICE: (h)?    * '
            WRITE(*,*) ' ******************************************** '
            WRITE(*,*)
            READ(*,9010) CHANS
            READ(CHANS,9020) HQ
            IF (HQ.EQ.'q' .OR. HQ.EQ.'Q') GO TO 460
            IF (HQ.EQ.'h' .OR. HQ.EQ.'H') THEN
               CALL HELP(7)
               GO TO 310
               END IF
            READ(CHANS,'(I32)') NANS
            IF (NANS.LT.1 .OR. NANS.GT.3) GO TO 310
         IF (NANS.EQ.1) THEN
            CC = 1.0D0
            IF (RITE) WRITE(2,*) ' The B.C. is Periodic. '
         ELSE IF (NANS.EQ.2) THEN
            CC = -1.0D0
            IF (RITE) WRITE(2,*) ' The B.C. is Semi-Periodic. '
         ELSE
  320       CONTINUE
               WRITE(*,*) ' ****************************************** '
               WRITE(*,*) ' * INPUT c: (h?)                          * '
               WRITE(*,*) ' ****************************************** '
               WRITE(*,*) ' c = '
               READ(*,9010) CHANS
               READ(CHANS,9020) HQ
               IF (HQ.EQ.'q' .OR. HQ.EQ.'Q') GO TO 460
               IF (HQ.EQ.'h' .OR. HQ.EQ.'H') THEN
                  CALL HELP(7)
                  GO TO 320
                  END IF
            READ(CHANS,'(F32.0)') CC
            IF (RITE) WRITE(2,*) ' The B.C. is General Periodic type. '
            IF (RITE) WRITE(2,*) ' Parameter c = ',CC
            END IF
         WRITE(*,*)
         WRITE(*,*)
         WRITE(*,*)
         WRITE(*,*)
  330    CONTINUE
            WRITE(*,*)
            WRITE(*,*) ' ********************************************* '
            WRITE(*,*) ' * INPUT NUMEIG,TOL: (h?)                    * '
            WRITE(*,*) ' ********************************************* '
            WRITE(*,*)
            WRITE(*,*) ' NUMEIG,TOL = '
            READ(*,9010) CHANS
            READ(CHANS,9020) HQ
            IF (HQ.EQ.'q' .OR. HQ.EQ.'Q') GO TO 460
            IF (HQ.EQ.'h' .OR. HQ.EQ.'H') THEN
               CALL HELP(17)
               GO TO 330
               END IF
         TOL = 0.0D0
         CALL LSTDIR(CHANS,I,ICOL)
         I1 = ICOL(1) - 1
         IF (I.EQ.2) THEN
            FMT = '(I'//COL(I1)//',1X,F'//COL(31-I1)//'.0)'
            READ(CHANS,FMT) NUMEIG,TOL
         ELSE
            FMT = '(I'//COL(I1)//')'
            READ(CHANS,FMT) NUMEIG
            END IF
         IF (RITE) WRITE(2,*) ' NUMEIG,TOL = ',NUMEIG,TOL
         IF (NUMEIG.LT.0) THEN
            WRITE(*,*) ' NUMEIG MUST BE .GE. 0 '
            GO TO 330
            END IF
         CALL PERIO(A,B,INTAB,P0ATA,QFATA,P0ATB,QFATB,
     1              A1,A2,B1,B2,NUMEIG,EIG,TOL,IFLAG,SLFN,
     2              SINGATA,SINGATB,CIRCLA,CIRCLB,OSCILA,OSCILB)
         IF (IFLAG.NE.1) IFLAG = 2
         WRITE(*,*)
         WRITE(*,*)
         WRITE(*,*)
         WRITE(*,*)
         WRITE(*,*)
         WRITE(*,*)' **************************************************'
         WRITE(*,*)'  * NUMEIG = ',NUMEIG,'   EIG =',EIG
         WRITE(*,9050) TOL,IFLAG
         WRITE(*,*)' **************************************************'
         IF (RITE) WRITE(2,*) '  * NUMEIG = ',NUMEIG,'   EIG =',EIG
         IF (RITE) WRITE(2,9050) TOL,IFLAG
         WRITE(*,*)
         WRITE(*,*)
         WRITE(*,*)
         WRITE(*,*)
  340    CONTINUE
         WRITE(*,*)
         WRITE(*,*) ' *************************************************'
         WRITE(*,*) ' * PLOT EIGENFUNCTION ? (Y/N) (h?)               *'
         WRITE(*,*) ' *************************************************'
         WRITE(*,*)
         READ(*,9010) CHANS
         READ(CHANS,9020) HQ
         IF (HQ.EQ.'q' .OR. HQ.EQ.'Q') GO TO 460
         IF (HQ.EQ.'h' .OR. HQ.EQ.'H') THEN
            CALL HELP(16)
            GO TO 340
            END IF
         READ(CHANS,9020) YN
         PEIGF = YN.EQ.'y' .OR. YN.EQ.'Y'
         IF (PEIGF) THEN
            THA = ATAN2(VB,CC-UB)
            IF (THA.LT.0.0D0) THA = THA + PI
            THB = ATAN(CC*CC*TAN(THA))
            IF (THB.LE.0.0D0) THB = THB + PI
            I = MOD(NUMEIG,2)
            IF ((CC.GT.0.0D0 .AND. I.EQ.1) .OR.
     1          (CC.LT.0.0D0 .AND. I.EQ.0)) NUMEIG = NUMEIG + 1
            A1 = COS(THA)
            A2 = -SIN(THA)
            B1 = COS(THB)
            B2 = -SIN(THB)
            SLFN(1) = 0.0D0
            SLFN(2) = -1.0D0
            SLFN(3) = THA
            SLFN(4) = 0.0D0
            SLFN(5) = 1.0D0
            SLFN(6) = THB
            SLFN(7) = 0.0D0
            SLFN(8) = .00001
            SLFN(9) = Z
            GO TO 420
            END IF
C
C     END OF COMPUTING AN EIGENVALUE.
C
         IF (RITE) WRITE(2,*) ' *******************************'//FILLA
         GO TO 430
C
C     THIS IS THE ENTRY POINT FOR PLOTTING A NEW INITIAL VALUE PROBLEM.
C
  350 CONTINUE
         WRITE(*,*)
         WRITE(*,*)
         WRITE(*,*)
         WRITE(*,*)
         WRITE(*,*)
         WRITE(*,*) ' *********************************************** '
         WRITE(*,*) ' * DO YOU WANT TO COMPUTE THE SOLUTION TO:     * '
         WRITE(*,*) ' *                                             * '
         WRITE(*,*) ' *    (1) AN INITIAL VALUE PROBLEM FROM ONE    * '
         WRITE(*,*) ' *        END OF THE INTERVAL TO THE OTHER     * '
         WRITE(*,*) ' *                                             * '
         WRITE(*,*) ' *    (2) INITIAL VALUE PROBLEMS FROM BOTH     * '
         WRITE(*,*) ' *        ENDS TO A MIDPOINT ?                 * '
         WRITE(*,*) ' *                                             * '
         WRITE(*,*) ' * ENTER THE NUMBER OF YOUR CHOICE: (h?)       * '
         WRITE(*,*) ' *********************************************** '
         WRITE(*,*)
         READ(*,9010) CHANS
         READ(CHANS,9020) HQ
         IF (HQ.EQ.'q' .OR. HQ.EQ.'Q') GO TO 460
         IF (HQ.EQ.'h' .OR. HQ.EQ.'H') THEN
            CALL HELP(12)
            GO TO 350
            END IF
         READ(CHANS,'(I32)') NIVP
         IF (NIVP.LT.1 .OR. NIVP.GT.2) GO TO 350
      IF (NIVP.EQ.1) THEN
         WRITE(*,*)
         WRITE(*,*)
         WRITE(*,*)
         WRITE(*,*)
  360    CONTINUE
            WRITE(*,*)
            WRITE(*,*) ' ********************************************* '
            WRITE(*,*) ' * WHICH IS THE INITIAL POINT: a OR b ? (h?) * '
            WRITE(*,*) ' ********************************************* '
            WRITE(*,*)
            WRITE(*,*) ' INITIAL POINT IS '
            READ(*,9010) CHANS
            READ(CHANS,9020) HQ
            IF (HQ.EQ.'q' .OR. HQ.EQ.'Q') GO TO 460
            IF (HQ.EQ.'h' .OR. HQ.EQ.'H') THEN
               CALL HELP(12)
               GO TO 360
               END IF
         READ(CHANS,9020) ANSCH
         IF (ANSCH.EQ.'a' .OR. ANSCH.EQ.'A') THEN
            NEND = 1
            IF (RITE) WRITE(2,*) ' The Initial Point for this',
     1                           ' Initial Value Problem is a. '
            IF (SINGATA.LT.0.0D0 .OR. CIRCLA.GT.0.0D0) THEN
               WRITE(*,*)
               WRITE(*,*)
               WRITE(*,*)
               WRITE(*,*)
  370          CONTINUE
                  WRITE(*,*)
                  WRITE(*,*) ' ************************************ '
                  WRITE(*,*) ' * THE INITIAL CONDITIONS AT a ARE  * '
                  WRITE(*,*) ' *                                  * '
                  IF (SINGATA.LT.0.0D0) THEN
                     WRITE(*,*) ' *    y(a)=alfa1, py''(a)=alfa2      *'
                  ELSE
                     WRITE(*,*) ' *  [y,u](a)=alfa1, [y,v](a)=alfa2  * '
                     END IF
                  WRITE(*,*) ' *                                  * '
                  WRITE(*,*) ' * WHERE THE CONSTANTS alfa1, alfa2 * '
                  WRITE(*,*) ' *    MAY BE CHOSEN ARBITRARILY.    * '
                  WRITE(*,*) ' *                                  * '
                  WRITE(*,*) ' * CHOOSE alfa1,alfa2: (h?)         * '
                  WRITE(*,*) ' ************************************ '
                  WRITE(*,*)
                  WRITE(*,*) ' alfa1,alfa2 = '
                  READ(*,9010) CHANS
                  READ(CHANS,9020) HQ
                  IF (HQ.EQ.'q' .OR. HQ.EQ.'Q') GO TO 460
                  IF (HQ.EQ.'h' .OR. HQ.EQ.'H') THEN
                    CALL HELP(12)
                    GO TO 370
                    END IF
               CALL LSTDIR(CHANS,I,ICOL)
               I1 = ICOL(1) - 1
               FMT = FMT2(I1)
               READ(CHANS,FMT) ALFA1,ALFA2
               IF (RITE) WRITE(2,*) ' alfa1,alfa2 = ',ALFA1,ALFA2
               A1 = ALFA2
               A2 = -ALFA1
               END IF
         ELSE IF (ANSCH.EQ.'b' .OR. ANSCH.EQ.'B') THEN
            NEND = 2
            IF (RITE) WRITE(2,*) ' The Initial Point for this',
     1                           ' Initial Value Problem is b. '
            IF (SINGATB.LT.0.0D0 .OR. CIRCLB.GT.0.0D0) THEN
               WRITE(*,*)
               WRITE(*,*)
               WRITE(*,*)
               WRITE(*,*)
  380          CONTINUE
                  WRITE(*,*)
                  WRITE(*,*) ' ************************************ '
                  WRITE(*,*) ' * THE INITIAL CONDITIONS AT b ARE  * '
                  WRITE(*,*) ' *                                  * '
                  IF (SINGATB.LT.0.0D0) THEN
                     WRITE(*,*) ' *    y(b)=beta1, py''(b)=beta2      *'
                  ELSE
                     WRITE(*,*) ' *  [y,u](b)=beta1, [y,v](b)=beta2  * '
                     END IF
                  WRITE(*,*) ' *                                  * '
                  WRITE(*,*) ' * WHERE THE CONSTANTS beta1, beta2 * '
                  WRITE(*,*) ' *    MAY BE CHOSEN ARBITRARILY.    * '
                  WRITE(*,*) ' *                                  * '
                  WRITE(*,*) ' * CHOOSE beta1,beta2: (h?)         * '
                  WRITE(*,*) ' ************************************ '
                  WRITE(*,*)
                  WRITE(*,*) ' beta1,beta2 = '
                  READ(*,9010) CHANS
                  READ(CHANS,9020) HQ
                  IF (HQ.EQ.'q' .OR. HQ.EQ.'Q') GO TO 460
                  IF (HQ.EQ.'h' .OR. HQ.EQ.'H') THEN
                     CALL HELP(12)
                     GO TO 380
                     END IF
               CALL LSTDIR(CHANS,I,ICOL)
               I1 = ICOL(1) - 1
               FMT = FMT2(I1)
               READ(CHANS,FMT) BETA1,BETA2
               IF (RITE) WRITE(2,*) ' beta1,beta2 = ',BETA1,BETA2
               B1 = BETA2
               B2 = -BETA1
               END IF
            END IF
      ELSE IF (NIVP.EQ.2) THEN
         NEND = 3
         IF (SINGATA.LT.0.0D0 .OR. CIRCLA.GT.0.0D0) THEN
            WRITE(*,*)
            WRITE(*,*)
            WRITE(*,*)
            WRITE(*,*)
  390       CONTINUE
               WRITE(*,*)
               WRITE(*,*) ' ************************************ '
               WRITE(*,*) ' * THE INITIAL CONDITIONS AT a ARE  * '
               WRITE(*,*) ' *                                  * '
               IF (SINGATA.LT.0.0D0) THEN
                  WRITE(*,*) ' *    y(a)=alfa1, py''(a)=alfa2      *'
               ELSE
                  WRITE(*,*) ' *  [y,u](a)=alfa1, [y,v](a)=alfa2  * '
                  END IF
               WRITE(*,*) ' *                                  * '
               WRITE(*,*) ' * WHERE THE CONSTANTS alfa1, alfa2 * '
               WRITE(*,*) ' *    MAY BE CHOSEN ARBITRARILY.    * '
               WRITE(*,*) ' *                                  * '
               WRITE(*,*) ' * CHOOSE alfa1,alfa2: (h?)         * '
               WRITE(*,*) ' ************************************ '
               WRITE(*,*)
               WRITE(*,*) ' alfa1,alfa2 = '
               READ(*,9010) CHANS
               READ(CHANS,9020) HQ
               IF (HQ.EQ.'q' .OR. HQ.EQ.'Q') GO TO 460
               IF (HQ.EQ.'h' .OR. HQ.EQ.'H') THEN
                 CALL HELP(12)
                 GO TO 390
                 END IF
            CALL LSTDIR(CHANS,I,ICOL)
            I1 = ICOL(1) - 1
            FMT = FMT2(I1)
            READ(CHANS,FMT) ALFA1,ALFA2
            A1 = ALFA2
            A2 = -ALFA1
            END IF
         IF (SINGATB.LT.0.0D0 .OR. CIRCLB.GT.0.0D0) THEN
            WRITE(*,*)
            WRITE(*,*)
            WRITE(*,*)
            WRITE(*,*)
  400       CONTINUE
               WRITE(*,*)
               WRITE(*,*) ' ************************************ '
               WRITE(*,*) ' * THE INITIAL CONDITIONS AT b ARE  * '
               WRITE(*,*) ' *                                  * '
               IF (SINGATB.LT.0.0D0) THEN
                  WRITE(*,*) ' *    y(b)=beta1, py''(b)=beta2      *'
               ELSE
                  WRITE(*,*) ' *  [y,u](b)=beta1, [y,v](b)=beta2  * '
                  END IF
               WRITE(*,*) ' *                                  * '
               WRITE(*,*) ' * WHERE THE CONSTANTS beta1, beta2 * '
               WRITE(*,*) ' *    MAY BE CHOSEN ARBITRARILY.    * '
               WRITE(*,*) ' *                                  * '
               WRITE(*,*) ' * CHOOSE beta1,beta2: (h?)         * '
               WRITE(*,*) ' ************************************ '
               WRITE(*,*)
               WRITE(*,*) ' beta1,beta2 = '
               READ(*,9010) CHANS
               READ(CHANS,9020) HQ
               IF (HQ.EQ.'q' .OR. HQ.EQ.'Q') GO TO 460
               IF (HQ.EQ.'h' .OR. HQ.EQ.'H') THEN
                  CALL HELP(12)
                  GO TO 400
                  END IF
            CALL LSTDIR(CHANS,I,ICOL)
            I1 = ICOL(1) - 1
            FMT = FMT2(I1)
            READ(CHANS,FMT) BETA1,BETA2
            B1 = BETA2
            B2 = -BETA1
            END IF
         END IF
      WRITE(*,*)
      WRITE(*,*)
      WRITE(*,*)
      WRITE(*,*)
C
C     THIS IS THE ENTRY POINT FOR PLOTTING A SOLUTION OF THE INITIAL
C     VALUE PROBLEM AFTER THE EIGENPARAMETER HAS BEEN CHOSEN.
C
  410 CONTINUE
         WRITE(*,*)
         WRITE(*,*) ' ************************************************ '
         WRITE(*,*) ' * WHAT VALUE SHOULD BE USED FOR THE            * '
         WRITE(*,*) ' *   EIGENPARAMETER, EIG ?                      * '
         WRITE(*,*) ' *                                              * '
         WRITE(*,*) ' * INPUT EIG = (h?)                             * '
         WRITE(*,*) ' ************************************************ '
         WRITE(*,*)
         WRITE(*,*) ' EIG = '
         READ(*,9010) CHANS
         READ(CHANS,9020) HQ
         IF (HQ.EQ.'q' .OR. HQ.EQ.'Q') GO TO 460
         IF (HQ.EQ.'h' .OR. HQ.EQ.'H') THEN
            CALL HELP(12)
            GO TO 410
            END IF
      READ(CHANS,'(F32.0)') EIG
C
C     THE FOLLOWING CALL SETS THE STAGE IN SLEIGN2 -- I.E., SAMPLES
C     THE COEFFICIENTS AND SETS THE INITIAL INTERVAL.
C
      NUMEIG = 0
      TOL = .001
      CALL SLEIGN2(A,B,INTAB,P0ATA,QFATA,P0ATB,QFATB,A1,A2,
     1             B1,B2,NUMEIG,EIG,TOL,IFLAG,-1,SLFN,
     2             SINGATA,SINGATB,CIRCLA,CIRCLB,OSCILA,OSCILB)
      IF (NIVP.EQ.1 .AND. NEND.EQ.1) TMID = BB
      IF (NIVP.EQ.1 .AND. NEND.EQ.2) TMID = AA
C
C     ACTUALLY, WE MAY HAVE TO AVOID AA OR BB IN SOME CASES,
C     WHICH WILL BE TAKEN CARE OF LATER.
C
      ALFA = 0.0D0
      IF ((NIVP.EQ.1 .AND. NEND.EQ.1 .AND. .NOT.SINGA) .OR. NIVP.EQ.2)
     1   ALFA = SLFN(3)
      BETA = PI
      IF ((NIVP.EQ.1 .AND. NEND.EQ.2 .AND. .NOT.SINGB) .OR. NIVP.EQ.2)
     1   BETA = SLFN(6)
      SLFN(3) = ALFA
      SLFN(6) = BETA
      SLFN(4) = 0.0D0
      SLFN(7) = 0.0D0
      SLFN(8) = .01
C
C     THE FOLLOWING CALL COMPUTES VALUES OF THE SOLUTION FOR PLOTTING.
C
  420 CONTINUE
      WRITE(*,*)
      CALL DRAW(A1,A2,B1,B2,NUMEIG,EIG,SLFN,SINGATA,SINGATB,CIRCLA,
     1          CIRCLB,OSCILA,OSCILB,REGA,REGB,NIVP,NEND,EIGV,RITE)
      WRITE(*,*)
C
  430 CONTINUE
      IF (EIGV) THEN
         WRITE(*,*) ' Press any key to continue. '
         READ(*,9010) CHANS
  440    CONTINUE
         CALL CHOICE(1)
         READ(*,9010) CHANS
         READ(CHANS,9020) HQ
         IF (HQ.EQ.'q' .OR. HQ.EQ.'Q') GO TO 460
         IF (HQ.EQ.'h' .OR. HQ.EQ.'H') THEN
            CALL HELP(3)
            GO TO 440
            END IF
         READ(CHANS,'(I32)') NANS
         IF (NANS.LT.1 .OR. NANS.GT.5) GO TO 440
         IF (NANS.EQ.1 .AND. PERIOD) GO TO 330
         IF (NANS.EQ.1 .AND. .NOT.PERIOD) GO TO 250
         IF (NANS.EQ.2) GO TO 180
         IF (NANS.EQ.3) GO TO 30
         EIGV = .FALSE.
         IF (NANS.EQ.4) GO TO 350
         GO TO 460
      ELSE
  450    CONTINUE
         WRITE(*,*) ' Press any key to continue. '
         READ(*,9010) CHANS
         CALL CHOICE(2)
         READ(*,9010) CHANS
         READ(CHANS,9020) HQ
         IF (HQ.EQ.'q' .OR. HQ.EQ.'Q') GO TO 460
         IF (HQ.EQ.'h' .OR. HQ.EQ.'H') THEN
            CALL HELP(3)
            GO TO 450
            END IF
         READ(CHANS,'(I32)') NANS
         IF (NANS.LT.1 .OR. NANS.GT.5) GO TO 450
         IF (NANS.EQ.1) GO TO 410
         IF (NANS.EQ.2) GO TO 350
         IF (NANS.EQ.3) GO TO 30
         EIGV = .TRUE.
         IF (NANS.EQ.4) GO TO 180
         END IF
  460 CONTINUE
      CLOSE(21)
      CLOSE(2)
      STOP
 9010 FORMAT(A32)
 9020 FORMAT(A1)
 9030 FORMAT(A70)
 9040 FORMAT(1X,'   NUMEIG = ',I5,'   IFLAG = ',I3)
 9050 FORMAT(1X,' * TOL = ',E14.5,2X,' IFLAG = ',I3,'             *')
 9060 FORMAT(1X,' NUMEIG = ',I5,2X,' EIG = ',E18.9,2X,' TOL = ',
     1  E14.5,' IFLAG = ',I3)
 9070 FORMAT(1X,1X,6H*    (,F12.7,1H,,F12.7,1H),'               *')
 9080 FORMAT(1X,1X,6H*    (,F12.7,1H,,A12  ,1H),'               *')
 9090 FORMAT(1X,1X,6H*    (,A12  ,1H,,F12.7,1H),'               *')
 9100 FORMAT(1X,1X,6H*    (,A12  ,1H,,A12  ,1H),'               *')
 9110 FORMAT(1X,2A39)
 9230 FORMAT(1X,' * THERE SEEMS TO BE NO EIGENVALUE OF INDEX       *'/
     1       1X,' * GREATER THAN',I5,'                              *')
 9250 FORMAT(1X,' * THERE MAY BE A CONTINUOUS SPECTRUM BEGINNING   *'/
     1       1X,' * AT ABOUT',1PE8.1,'                               *')
      END
      SUBROUTINE CHOICE(I)
      INTEGER I
C     **********
C     THIS PROGRAM DISPLAYS THE CHOICES FOR PROBLEM CONTINUATION.
C     **********
      WRITE(*,*) ' *********************************************** '
      WRITE(*,*) ' * WHAT WOULD YOU LIKE TO DO NOW ?             * '
      WRITE(*,*) ' *                                             * '
      IF (I.EQ.1) THEN
         WRITE(*,*) ' *    (1)  SAME EIGENVALUE PROBLEM, DIFFERENT  * '
         WRITE(*,*) ' *           NUMEIG, EIG, OR TOL               * '
         WRITE(*,*) ' *                                             * '
         WRITE(*,*) ' *    (2)  SAME EIGENVALUE PROBLEM, SAME (a,b) * '
         WRITE(*,*) ' *           AND p,q,w,u,v BUT DIFFERENT       * '
         WRITE(*,*) ' *           BOUNDARY CONDITIONS A1,A2,B1,B2   * '
         WRITE(*,*) ' *                                             * '
         WRITE(*,*) ' *    (3)  INTERVAL CHANGE, PROBLEM RESTART    * '
         WRITE(*,*) ' *                                             * '
         WRITE(*,*) ' *    (4)  AN INITIAL VALUE PROBLEM            * '
      ELSE
         WRITE(*,*) ' *    (1)  SAME INITIAL VALUE PROBLEM,         * '
         WRITE(*,*) ' *           DIFFERENT LAMBDA                  * '
         WRITE(*,*) ' *                                             * '
         WRITE(*,*) ' *    (2)  NEW INITIAL VALUE PROBLEM           * '
         WRITE(*,*) ' *                                             * '
         WRITE(*,*) ' *    (3)  INTERVAL CHANGE, PROBLEM RESTART    * '
         WRITE(*,*) ' *                                             * '
         WRITE(*,*) ' *    (4)  AN EIGENVALUE PROBLEM               * '
         END IF
      WRITE(*,*) ' *                                             * '
      WRITE(*,*) ' *    (5)  QUIT                                * '
      WRITE(*,*) ' *                                             * '
      WRITE(*,*) ' * ENTER THE NUMBER OF YOUR CHOICE: (h?)       * '
      WRITE(*,*) ' *********************************************** '
      WRITE(*,*)
      RETURN
      END
      SUBROUTINE DRAW(A1,A2,B1,B2,
     1                NUMEIG,EIG,SLFN,SINGATA,SINGATB,CIRCLA,CIRCLB,
     2                OSCILA,OSCILB,REGA,REGB,NIVP,NEND,EIGV,RITE)
      INTEGER NUMEIG,NIVP,NEND
      LOGICAL REGA,REGB,EIGV,RITE
      DOUBLE PRECISION A1,A2,B1,B2,EIG,
     1     SINGATA,SINGATB,CIRCLA,CIRCLB,OSCILA,OSCILB
      DOUBLE PRECISION SLFN(10)
C     **********
C     **********
C     .. Scalars in Common ..
      INTEGER IND,MDTHZ
      LOGICAL ADDD
      DOUBLE PRECISION AA,BB,DTHDAA,DTHDBB,EIGSAV,HPI,PI,TMID,TWOPI
C     ..
C     .. Arrays in Common ..
      INTEGER NT(2)
      DOUBLE PRECISION TT(7,2),YY(7,3,2)
C     ..
C     .. Local Scalars ..
      INTEGER I,ISLFUN,JJ,K,KFLAG,MM,NF,NPTS,NV
      LOGICAL ENDA,ENDB,WREGA,WREGB,LCOA,LCOB,LCIRCA,LCIRCB,
     1        OSCA,OSCB,SINGA,SINGB
      CHARACTER*1 HQ,YN
      CHARACTER*32 CHANS,TAPE
      DOUBLE PRECISION BRYU,BRYV,CA,CB,DD,EIGPI,FA,FB,FAC,FZ,GEE,
     1     HU,HUI,HV,HVI,PHI,PUP,PUPI,PVP,PVPI,PYP,RHO,RHOZ,
     2     SIG,SQ,T,THETAZ,TH,TI,TMP,U,UI,V,VI,X,XI,XJMP,XJMPS,Y,Z
C     ..
C     .. Local Arrays ..
      CHARACTER*55 XC(8)
      DOUBLE PRECISION SLFUN(1000,2),PLOTF(1000,6),XT(1000,2),
     1     UTHZ(3),UTH(3)
C     ..
C     .. External Subroutines ..
      EXTERNAL DXDT,EIGENF,HELP,MESH,QPLOT,THZTOTH,UV
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS,INT
C     ..
C     .. Common blocks ..
      COMMON /TEMP/TT,YY,NT
      COMMON /DATAF/EIGSAV,IND
      COMMON /PIE/PI,TWOPI,HPI
      COMMON /TDATA/AA,TMID,BB,DTHDAA,DTHDBB,ADDD,MDTHZ
C     ..
C     Definition of some logicals.
C
      SINGA = SINGATA.GT.0.0D0
      WREGA = .NOT.REGA .AND. .NOT.SINGA
      SINGB = SINGATB.GT.0.0D0
      WREGB = .NOT.REGB .AND. .NOT.SINGB
      LCIRCA = CIRCLA.GT.0.0D0
      LCIRCB = CIRCLB.GT.0.0D0
      OSCA = OSCILA.GT.0.0D0
      LCOA = LCIRCA .AND. OSCA
      OSCB = OSCILB.GT.0.0D0
      LCOB = LCIRCB .AND. OSCB
      EIGPI = NUMEIG*PI
C
C     Definition of some character strings.
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
C     NV = 1 MEANS THE INDEPENDENT VARIABLE IS X.
C     NV = 2 MEANS THE INDEPENDENT VARIABLE IS T.
C
C     NF = 1 MEANS THE EIGENFUNCTION Y IS WANTED.
C     NF = 2 MEANS THE QUASI-DERIVATIVE P*Y' IS WANTED.
C     NF = 3 MEANS BOUNDARY CONDITION FUNCTION Y OR [Y,U] IS WANTED.
C     NF = 4 MEANS BOUNDARY CONDITION FUNCTION P*Y' OR [Y,V] IS WANTED.
C     NF = 5 MEANS THE PRUFER ANGLE THETA IS WANTED.
C     NF = 6 MEANS THE PRUFER MODULUS RHO IS WANTED.
C
C     IF AN EIGENVALUE HAS BEEN COMPUTED (EIGV = .TRUE.),
C     THE RELEVANT VALUES HAVE BEEN STORED IN ARRAY SLFN, AND
C     NEED TO BE COPIED INTO THE FIRST COLUMN OF ARRAY SLFUN.
C
      DO 10 I = 1,9
         SLFUN(I,1) = SLFN(I)
   10    CONTINUE
C
      CALL MESH(EIG,SLFUN(10,1),ISLFUN)
C
C     THE POINTS GENERATED BY MESH MAY NOT BE WITHIN THE INTERVAL
C     (AA,BB).  WE WANT TO ENSURE THAT WE USE ONLY SUCH POINTS
C     AS ARE WITHIN THE INTERVAL.
C
      K = 0
      DO 20 I = 1,ISLFUN
         IF (SLFUN(9+I,1).GT.AA .AND. SLFUN(9+I,1).LT.BB) THEN
            K = K + 1
            SLFUN(9+K,1) = SLFUN(9+I,1)
            END IF
   20    CONTINUE
      ISLFUN = K
C
C     WE ALSO WANT TO BE SURE AN INITIAL ENDPOINT IS AA OR BB UNLESS THE
C     POINT IS LIMIT CIRCLE, AND THAT THE LAST POINT IS NOT AA OR BB.
C
      IF ((EIGV .OR. NEND.EQ.1 .OR. NIVP.EQ.2) .AND. .NOT.LCIRCA) THEN
         ISLFUN = ISLFUN + 1
         DO 30 I = ISLFUN,2,-1
            SLFUN(9+I,1) = SLFUN(8+I,1)
   30       CONTINUE
         SLFUN(9+1,1) = AA
         END IF
      IF ((EIGV .OR. NEND.EQ.2 .OR. NIVP.EQ.2) .AND. .NOT.LCIRCB) THEN
         SLFUN(ISLFUN+10,1) = BB
         ISLFUN = ISLFUN + 1
         END IF
C
C     NEAR AN OSCILLATORY ENDPOINT THE POINTS MAY BE SO CLOSE THAT WE
C     COULDN'T SEE THE ACTUAL CURVE EVEN IF WE PLOTTED IT.  SO WE WANT
C     TO REMOVE POINTS FROM THAT END UP TO WHERE THEY ARE NOT SO CLOSE.
C
      IF (OSCA) THEN
         JJ = 0
         DO 40 I = 1,ISLFUN
            IF (ABS(SLFUN(9+I,1)-SLFUN(10+I,1)).LT.0.0D001) JJ = JJ + 1
   40       CONTINUE
         IF (JJ.GT.0) THEN
            ISLFUN = ISLFUN - JJ
            DO 50 I = 1,ISLFUN
            SLFUN(9+I,1) = SLFUN(9+I+JJ,1)
   50       CONTINUE
            END IF
         END IF
      IF (OSCB) THEN
         JJ = 0
         DO 60 I = ISLFUN,1,-1
            IF (ABS(SLFUN(9+I,1)-SLFUN(8+I,1)).LT.0.0D001) JJ = JJ + 1
   60       CONTINUE
         IF (JJ.GT.0) ISLFUN = ISLFUN - JJ
         END IF
C
C     FINALLY, IN THE CASE OF AN INITIAL VALUE PROBLEM,
C     WE CANNOT AFFORD TO INTEGRATE TO THE OTHER END B (OR A)
C     UNLESS BOTH ENDS ARE REGULAR.
C
      IF (NIVP.EQ.1 .AND. .NOT.(REGA .AND. REGB)) THEN
         ISLFUN = ISLFUN - 1
         IF (NEND.EQ.2) THEN
            DO 70 I = 1,ISLFUN
            SLFUN(9+I,1) = SLFUN(10+I,1)
   70       CONTINUE
            END IF
         END IF
C
      DO 80 I = 1,ISLFUN
         XT(9+I,2) = SLFUN(9+I,1)
   80    CONTINUE
C
C     WARNING: THE VALUES RETURNED IN SLFUN BY EIGENF DEPEND
C              ON THE VALUE OF KFLAG IN THE CALL TO EIGENF:
C
C              IF KFLAG = 1, THE VALUES IN SLFUN(9+I,1) ARE THE
C                            EIGENFUNCTION Y ITSELF.
C
C              IF KFLAG = 2, THE VALUES IN THE TWO COLUMNS OF SLFUN ARE
C                               SLFUN(9+I,1) = THETA(9+I)
C                               SLFUN(9+I,2) = EFF(9+I)
C                     WHERE
C                               RHO  = EXP(EFF)
C                                Y   = RHO*SIN(THETA)
C                               P*Y' = RHO*COS(THETA)*Z
C
C     HERE, WE SET KFLAG = 2 SO THAT EIGENF WILL RETURN EFF, ENABLING
C     US TO DEAL WITH IT BEFORE FORMING THE FUNCTION Y = RHO*SIN .
C
      KFLAG = 2
      EIGSAV = EIG
      CALL EIGENF(EIGPI,A1,A2,B1,B2,REGA,SINGA,LCIRCA,OSCA,
     1            REGB,SINGB,LCIRCB,OSCB,SLFUN,ISLFUN,KFLAG)
C
C     IT MAY HAPPEN THAT THETA(I) HAS A JUMP OF MM*PI AT TMID.
C     IN THIS CASE, WE NEED TO SUBTRACT MM*PI FROM THETA(I).
C
      IF (EIGV) THEN
         JJ = 0
         XJMPS = HPI
         DO 90 I = 1,ISLFUN-1
            XJMP = SLFUN(10+I,1) - SLFUN(9+I,1)
            IF (ABS(XJMP).GT.XJMPS) THEN
               XJMPS = ABS(XJMP)
               JJ = I
               END IF
   90       CONTINUE
         IF (JJ.NE.0) THEN
            XJMP = SLFUN(10+JJ,1) - SLFUN(9+JJ,1)
            MM = INT(XJMP/PI)
            IF ((SLFUN(10+JJ,1)-MM*PI).LT.SLFUN(9+JJ,1)) MM = MM - 1
            DO 100 I = JJ,ISLFUN-1
               SLFUN(10+I,1) = SLFUN(10+I,1) - MM*PI
  100          CONTINUE
            END IF
         END IF
C***********************************************************************
C                                                                      *
C  SLEIGN2 NORMALIZES THE WAVEFUNCTION TO HAVE L2-NORM 1.0D0, BUT FOR AN
C  INITIAL VALUE PROBLEM WE WANT TO HAVE THE VALUE (Y) AND SLOPE (P*Y')
C  (OR [Y,U] AND [Y,V]) AT THE END TO BE THOSE SPECIFIED FOR THE
C  INITIAL CONDITIONS.  SO HERE WE MUST RE-NORMALIZE FOR THIS PURPOSE.
C
C     THE VALUES IN THE ARRAYS TT AND YY COME FROM THE INTEGRATIONS
C     IN SUBROUTINE WR.
C
      ENDA = (NEND.EQ.1 .OR. NIVP.EQ.2) .AND. (LCIRCA .OR. WREGA)
      IF (ENDA) THEN
         NPTS = 7
         IF (LCOA) NPTS = NT(1)
         DO 110 I = 2,NPTS
            T = TT(I,1)
            PHI = YY(I,1,1)
            GEE = YY(I,3,1)
            SIG = EXP(GEE)
            IF (LCIRCA) THEN
               CALL DXDT(T,TMP,X)
               CALL UV(X,U,PUP,V,PVP,HU,HV)
               DD = U*PVP - V*PUP
               SQ = SQRT(A1**2+A2**2)
               BRYU = -DD*SIN(PHI)*SIG*SQ
               BRYV = DD*COS(PHI)*SIG*SQ
               END IF
  110       CONTINUE
         END IF
      ENDB = (NEND.EQ.2 .OR. NIVP.EQ.2) .AND. (LCIRCB .OR. WREGB)
      IF (ENDB) THEN
         NPTS = 7
         IF (LCOB) NPTS = NT(2)
         DO 120 I = 2,NPTS
            T = TT(I,2)
            PHI = YY(I,1,2)
            GEE = YY(I,3,2)
            SIG = EXP(GEE)
            IF (LCIRCB) THEN
               CALL DXDT(T,TMP,X)
               CALL UV(X,U,PUP,V,PVP,HU,HV)
               DD = U*PVP - V*PUP
               SQ = SQRT(B1**2+B2**2)
               BRYU = -DD*SIN(PHI)*SIG*SQ
               BRYV = DD*COS(PHI)*SIG*SQ
               END IF
  120       CONTINUE
         END IF
C
      Z = SLFUN(9,1)
      ENDA = (NEND.EQ.1 .OR. NIVP.EQ.2) .AND. LCIRCA
      CA = 0.0D0
      IF (.NOT.(EIGV .OR. ENDA) .AND. NEND.NE.2) THEN
         FA = SLFUN(4,1)
         CA = 0.5*LOG(A2**2+(A1/Z)**2) - FA
         END IF
      ENDB = (NEND.EQ.2 .OR. NIVP.EQ.2) .AND. LCIRCB
      CB = 0.0D0
      IF (.NOT.(EIGV .OR. ENDB) .AND. NEND.NE.1) THEN
         FB = SLFUN(7,1)
         CB = 0.5*LOG(B2**2+(B1/Z)**2) - FB
         END IF
C
      DO 130 I = 1,ISLFUN
         IF (XT(9+I,2).LE.TMID .AND. .NOT.LCIRCA) THEN
            SLFUN(9+I,2) = SLFUN(9+I,2) + CA
         ELSE IF (XT(9+I,2).GT.TMID .AND. .NOT.LCIRCB) THEN
            SLFUN(9+I,2) = SLFUN(9+I,2) + CB
            END IF
  130    CONTINUE
C
         UTHZ(2) = 0.0D0
         UTHZ(3) = 0.0D0
         DO 140 I = 1,ISLFUN
            TI = XT(9+I,2)
            CALL DXDT(TI,TMP,XI)
            XT(9+I,1) = XI
            THETAZ = SLFUN(9+I,1)
            FZ = SLFUN(9+I,2)
            RHOZ = EXP(FZ)
            Y = RHOZ*SIN(THETAZ)
            PYP = Z*RHOZ*COS(THETAZ)
            RHO = SQRT(Y**2+PYP**2)
            PLOTF(9+I,1) = Y
            PLOTF(9+I,2) = PYP
            PLOTF(9+I,3) = Y
            PLOTF(9+I,4) = PYP
            UTHZ(1) = THETAZ
            CALL THZTOTH(UTHZ,Z,UTH)
            TH = UTH(1)
            PLOTF(9+I,5) = TH
            PLOTF(9+I,6) = RHO
C
            IF ((ENDA .AND. TI.LE.TMID) .OR.
     1          (ENDB .AND .TI.GT.TMID)) THEN
               CALL UV(XI,UI,PUPI,VI,PVPI,HUI,HVI)
               DD = UI*PVPI - VI*PUPI
               BRYU = PUPI*Y - PYP*UI
               BRYV = PVPI*Y - PYP*VI
C
C     RENORMALIZE.
C
               IF (ENDA .AND. TI.LE.TMID) THEN
                  FAC = SQRT(A1**2+A2**2)
               ELSE
                  FAC = SQRT(B1**2+B2**2)
                  END IF
               BRYU = FAC*BRYU
               BRYV = FAC*BRYV
               PLOTF(9+I,3) = BRYU
               PLOTF(9+I,4) = BRYV
               END IF
  140       CONTINUE
C
C  IN THE CASE OF NIVP = 2  WE HAVE COMPUTED THE SOLUTIONS
C  TO TWO INITIAL VALUE PROBLEMS FROM THE TWO ENDS.
C  IF (.NOT.ENDA .AND. .NOT.ENDB) WE SHOULD NOW HAVE
C
C      Y(A) = -A2 (ALFA1)  ;  Y(B) = -B2 (BETA1)
C     PY'(A) = A1 (ALFA2)  ;  PY'(B) = B1 (BETA2)
C
C  IT SOMETIMES HAPPENS THAT THE VALUES AT THE FINAL END HAVE GROWN
C  SO LARGE THAT THE REST OF THE CURVE IS TOTALLY SUBMERGED.
C  IN SUCH CASES WE WANT TO REMOVE THE POINTS WITH THE LARGE VALUES.
C                                                                    *
C*********************************************************************
  150 CONTINUE
      WRITE(*,*) ' ****************************************************'
      WRITE(*,*) ' * WHICH FUNCTION DO YOU WANT TO PLOT ?             *'
      WRITE(*,*) ' *                                                  *'
      WRITE(*,*) XC(1)
      WRITE(*,*) XC(2)
      WRITE(*,*) XC(3)
      WRITE(*,*) XC(4)
      WRITE(*,*) XC(5)
      WRITE(*,*) XC(6)
      WRITE(*,*) ' *                                                  *'
      WRITE(*,*) ' * ENTER THE NUMBER OF YOUR CHOICE: (h?)            *'
      WRITE(*,*) ' ****************************************************'
      WRITE(*,*)
      READ(*,1) CHANS
      READ(CHANS,2) HQ
      IF (HQ.EQ.'q' .OR. HQ.EQ.'Q') RETURN
      IF (HQ.EQ.'h' .OR. HQ.EQ.'H') THEN
         CALL HELP(16)
         GO TO 150
         END IF
      READ(CHANS,'(I32)') NF
  160 CONTINUE
      WRITE(*,*) ' ****************************************************'
      WRITE(*,*) ' * WHICH DO YOU WANT AS THE INDEPENDENT VARIABLE ?  *'
      WRITE(*,*) ' *                                                  *'
      WRITE(*,*) XC(7)
      WRITE(*,*) XC(8)
      WRITE(*,*) ' *                                                  *'
      WRITE(*,*) ' * ENTER THE NUMBER OF YOUR CHOICE: (h?)            *'
      WRITE(*,*) ' ****************************************************'
      WRITE(*,*)
C
      READ(*,1) CHANS
      READ(CHANS,2) HQ
      IF (HQ.EQ.'q' .OR. HQ.EQ.'Q') RETURN
      IF (HQ.EQ.'h' .OR. HQ.EQ.'H') THEN
         CALL HELP(16)
         GO TO 160
         END IF
      READ(CHANS,'(I32)') NV
      CALL QPLOT(ISLFUN,XT,NV,PLOTF,NF)
C
  170 CONTINUE
         WRITE(*,*)
         WRITE(*,*) ' DO YOU WANT TO SAVE THE PLOT FILE ? (Y/N) (h?)'
         READ(*,1) CHANS
         READ(CHANS,2) HQ
         IF (HQ.EQ.'q' .OR. HQ.EQ.'Q') STOP
         IF (HQ.EQ.'h' .OR. HQ.EQ.'H') THEN
            CALL HELP(16)
            GO TO 170
            END IF
      READ(CHANS,2) YN
      IF (YN.EQ.'y' .OR. YN.EQ.'Y') THEN
         WRITE(*,*) ' SPECIFY NAME OF FILE FOR PLOTTING '
         READ(*,'(A)') TAPE
         OPEN(1,FILE=TAPE,STATUS='NEW')
         DO 180 I = 1,ISLFUN
            WRITE(1,*) XT(9+I,NV),PLOTF(9+I,NF)
  180       CONTINUE
         CLOSE(1)
         WRITE(*,*) ' THE PLOT FILE HAS BEEN WRITTEN TO ',TAPE
         IF (RITE) WRITE(2,*) ' The plot file has been written to ',TAPE
         END IF
      WRITE(*,*)
  190 CONTINUE
         WRITE(*,*) ' PLOT ANOTHER FUNCTION ? (Y/N) (h?)'
         READ(*,1) CHANS
         READ(CHANS,2) HQ
         IF (HQ.EQ.'q' .OR. HQ.EQ.'Q') STOP
         IF (HQ.EQ.'h' .OR. HQ.EQ.'H') THEN
           CALL HELP(16)
           GO TO 190
           END IF
      READ(CHANS,2) YN
      IF (YN.EQ.'y' .OR. YN.EQ.'Y') GO TO 150
      RETURN
   1  FORMAT(A32)
   2  FORMAT(A1)
      END
      SUBROUTINE EIGENF(EIGPI,A1,A2,B1,B2,AOK,SINGA,LCIRCA,OSCA,
     1                  BOK,SINGB,LCIRCB,OSCB,SLFUN,ISLFUN,KFLAG)
      INTEGER ISLFUN,KFLAG
      LOGICAL AOK,SINGA,LCIRCA,OSCA,BOK,SINGB,LCIRCB,OSCB
      DOUBLE PRECISION EIGPI,A1,A2,B1,B2
      DOUBLE PRECISION SLFUN(1000,2)
C     **********
C     WARNING: DANGER! The array SLFUN here is two-dimensional, whereas
C                      it is one-dimensional in subroutine SLEIGN2.
C     **********
C     .. Scalars in Common ..
      INTEGER MDTHZ
      LOGICAL ADDD
      DOUBLE PRECISION AA,BB,DTHDAA,DTHDBB,HPI,PI,TMID,TWOPI
C     ..
C     .. Local Scalars ..
      INTEGER I,IFLAG,J,NMID
      LOGICAL LCIRC,OK,SING
      DOUBLE PRECISION DTHDAT,DTHDBT,DTHDET,EFF,T,THT,TM
C     ..
C     .. Local Arrays ..
      DOUBLE PRECISION ERL(3),ERR(3),YL(3),YR(3)
C     ..
C     .. External Subroutines ..
      EXTERNAL INTEG
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC EXP,SIN
C     ..
C     .. Common blocks ..
      COMMON /PIE/PI,TWOPI,HPI
      COMMON /TDATA/AA,TMID,BB,DTHDAA,DTHDBB,ADDD,MDTHZ
C     ..
C
C     WARNING: In this program it is assumed that the points T
C              in SLFUN all lie within the interval (AA,BB).
C
C     Calculate selected eigenfunction values by integration (over T).
C
      NMID = 0
      DO 10 I=1,ISLFUN
         IF (SLFUN(9+I,1).LE.TMID) NMID = I
   10    CONTINUE
      IF (NMID.GT.0) THEN
         T = AA
         YL(1) = SLFUN(3,1)
         YL(2) = 0.0D0
         YL(3) = 0.0D0
         LCIRC = LCIRCA
         OK = AOK
         SING = SINGA
         EFF = 0.0D0
         DO 20 J=1,NMID
            TM = SLFUN(J+9,1)
            IF (TM.LT.AA .OR. TM.GT.BB) THEN
               WRITE(*,*) ' t.lt.aa .or. t.gt.bb '
               STOP
               END IF
            THT = YL(1)
            DTHDAT = DTHDAA*EXP(-2.0*EFF)
            DTHDET = YL(2)
            IF (TM.GT.AA) THEN
               CALL INTEG(T,THT,DTHDAT,DTHDET,TM,A1,A2,SLFUN(8,1),
     1                    YL,ERL,LCIRC,OK,SING,OSCA,IFLAG)
               IF (OSCA) THEN
                  EFF = YL(3)
               ELSE
                  LCIRC = .FALSE.
                  SING = .FALSE.
                  EFF = EFF + YL(3)
                  END IF
               END IF
            IF (KFLAG.EQ.1) THEN
               SLFUN(J+9,1) = SIN(YL(1))*EXP(EFF+SLFUN(4,1))
            ELSE
               SLFUN(J+9,1) = YL(1)
               SLFUN(J+9,2) = EFF + SLFUN(4,1)
               END IF
            T = TM
            IF (T.GT.-1.0D0) OK = .TRUE.
            IF (T.LT.-0.9 .AND .OSCA) THEN
               OK = .FALSE.
               T = AA
               YL(1) = SLFUN(3,1)
               YL(2) = 0.0D0
               YL(3) = 0.0D0
               END IF
   20       CONTINUE
         END IF
      IF (NMID.LT.ISLFUN) THEN
         T = BB
         YR(1) = SLFUN(6,1) - EIGPI
         YR(2) = 0.0D0
         YR(3) = 0.0D0
         LCIRC = LCIRCB
         OK = BOK
         SING = SINGB
         EFF = 0.0D0
         DO 30 J=ISLFUN,NMID+1,-1
            TM = SLFUN(J+9,1)
            IF (TM.LT.AA .OR. TM.GT.BB) THEN
               WRITE(*,*) ' t.lt.aa .or. t.gt.bb '
               STOP
               END IF
            THT = YR(1)
            DTHDBT = DTHDBB*EXP(-2.0*EFF)
            DTHDET = YR(2)
            IF (TM.LT.BB) THEN
               CALL INTEG(T,THT,DTHDBT,DTHDET,TM,B1,B2,SLFUN(8,1),
     1                    YR,ERR,LCIRC,OK,SING,OSCB,IFLAG)
               IF (OSCB) THEN
                  EFF = YR(3)
               ELSE
                  LCIRC = .FALSE.
                  SING = .FALSE.
                  EFF = EFF + YR(3)
                  END IF
               END IF
            IF (KFLAG.EQ.1) THEN
               SLFUN(J+9,1) = SIN(YR(1)+EIGPI)*EXP(EFF+SLFUN(7,1))
               IF (ADDD) SLFUN(J+9,1) = -SLFUN(J+9,1)
            ELSE
               SLFUN(J+9,1) = YR(1) + EIGPI
               IF (ADDD) SLFUN(J+9,1) = SLFUN(J+9,1) + PI
               IF (OSCA .OR. OSCB) SLFUN(J+9,1) =
     1            SLFUN(J+9,1) - MDTHZ*PI
               SLFUN(J+9,2) = EFF + SLFUN(7,1)
               END IF
            T = TM
            IF (T.LT.1.0D0) OK = .TRUE.
            IF (T.GT.0.9 .AND. OSCB) THEN
               OK = .FALSE.
               T = BB
               YR(1) = SLFUN(6,1) - EIGPI
               YR(2) = 0.0D0
               YR(3) = 0.0D0
               END IF
   30       CONTINUE
         END IF
      RETURN
      END
      SUBROUTINE FZERO(F,B,C,R,RE,AE,IFLAG)
      INTEGER IFLAG
      DOUBLE PRECISION F,B,C,R,RE,AE
      EXTERNAL F
C     **********
C     **********
C     .. Local Scalars ..
      INTEGER IC,KOUNT
      DOUBLE PRECISION A,ACBS,ACMB,AW,CMB,DIF,DIFS,FA,FB,FC,FX,FZ,
     1     P,Q,RW,TOL,Z
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS,MAX,MIN,SIGN
C     ..
      DIF = 1000.0D0
      Z = R
      IF (R.LE.MIN(B,C).OR.R.GE.MAX(B,C)) Z = C
      RW = MAX(RE,0.0D0)
      AW = MAX(AE,0.0D0)
      IC = 0
      FZ = F(Z)
      FB = F(B)
      KOUNT = 2
      IF (FZ*FB.LT.0.0D0) THEN
         C = Z
         FC = FZ
      ELSE IF (Z.NE.C) THEN
         FC = F(C)
         KOUNT = 3
         IF (FZ*FC.LT.0.0D0) THEN
            B = Z
            FB = FZ
            END IF
         END IF
      A = C
      FA = FC
      ACBS = ABS(B-C)
      FX = MAX(ABS(FB),ABS(FC))
C
   10 CONTINUE
         IF (ABS(FC).LT.ABS(FB)) THEN
            A = B
            FA = FB
            B = C
            FB = FC
            C = A
            FC = FA
            END IF
         CMB = 0.5*(C-B)
         ACMB = ABS(CMB)
         TOL = RW*ABS(B) + AW
         IFLAG = 1
         IF (ACMB.LE.TOL) THEN
            IF (FB*FC.GE.0.0D0) IFLAG = 4
            IF (ABS(FB).GT.FX) IFLAG = 3
            RETURN
            END IF
         IFLAG = 2
         IF (FB.EQ.0.0D0) RETURN
         IFLAG = 5
         IF (KOUNT.GE.500) RETURN
C
         P = (B-A)*FB
         Q = FA - FB
         IF (P.LT.0.0D0) THEN
            P = -P
            Q = -Q
            END IF
         A = B
         FA = FB
         IC = IC + 1
         IF (IC.GE.4) THEN
            IF (8.0*ACMB.GE.ACBS) B = 0.5*(C+B)
         ELSE
            IC = 0
            ACBS = ACMB
            IF (P.LE.ABS(Q)*TOL) THEN
               B = B + SIGN(TOL,CMB)
            ELSE IF (P.GE.CMB*Q) THEN
               B = 0.5*(C+B)
            ELSE
               B = B + P/Q
               END IF
            END IF
         FB = F(B)
         DIFS = DIF
         DIF = FB - FC
         IF (DIF.EQ.DIFS) THEN
            IFLAG = 6
            RETURN
            END IF
         KOUNT = KOUNT + 1
         IF (FB*FC.GE.0.0D0) THEN
            C = A
            FC = FA
            END IF
         GO TO 10
      END
      subroutine help(nh)
      integer i,n,nh
      character*36 x(23),y(23)
      character*1 ans
c
      GO TO (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17),NH
c
    1 CONTINUE
      write(*,*) 'H1:  Overview of HELP.'
      x(1)='  This ASCII text file is supplied a'
      y(1)='s a separate file with the SLEIGN2  '
      x(2)='package; it can be accessed on-line '
      y(2)='in both MAKEPQW (if used) and DRIVE.'
      x(3)='  HELP contains information to aid t'
      y(3)='he user in entering data on the     '
      x(4)='coefficient functions p,q,w; on the '
      y(4)='limit circle boundary condition     '
      x(5)='functions u,v; on the end-point clas'
      y(5)='sifications of the differential     '
      x(6)='equation; on DEFAULT entry; on eigen'
      y(6)='value indexes; on IFLAG information;'
      x(7)='and on the general use of the progra'
      y(7)='m SLEIGN2.                          '
      x(8)='  The 17 sections of HELP are:      '
      y(8)=' '
      x(9)=' '
      y(9)=' '
      x(10)='    H1: Overview of HELP.           '
      y(10)=' '
      x(11)='    H2: File name entry.            '
      y(11)=' '
      x(12)='    H3: The differential equation.  '
      y(12)=' '
      x(13)='    H4: End-point classification.   '
      y(13)=' '
      x(14)='    H5: DEFAULT entry.              '
      y(14)=' '
      x(15)='    H6: Limit-circle boundary condit'
      y(15)='ions.                               '
      do 101 i = 1,15
        write(*,*) x(i),y(i)
  101   continue
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*)
      read(*,999) ans,n
      if (ans.eq.'r' .or. ans.eq.'R') return
      GO TO (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17),N
      x(1)='    H7: General boundary conditions.'
      y(1)=' '
      x(2)='    H8: Recording the results.      '
      y(2)=' '
      x(3)='    H9: Type and choice of interval.'
      y(3)=' '
      x(4)='   H10: Entry of end-points.        '
      y(4)=' '
      x(5)='   H11: End-point values of p,q,w.  '
      y(5)=' '
      x(6)='   H12: Initial value problems.     '
      y(6)=' '
      x(7)='   H13: Indexing of eigenvalues.    '
      y(7)=' '
      x(8)='   H14: Entry of eigenvalue index, i'
      y(8)='nitial guess, and tolerance.     '
      x(9)='   H15: IFLAG information.          '
      y(9)=' '
      x(10)='   H16: Plotting.                   '
      y(10)=' '
      x(11)='   H17: Indexing of eigenvalues for '
      y(11)='periodic-type problems.             '
      x(12)=' '
      y(12)=' '
      x(13)='  HELP can be accessed at each point'
      y(13)=' in MAKEPQW and DRIVE where the user'
      x(14)='is asked for input, by pressing "h <'
      y(14)='ENTER>"; this places the user at the'
      x(15)='appropriate HELP section.  Once in H'
      y(15)='ELP, the user can scroll the further'
      x(16)='HELP sections by repeatedly pressing'
      y(16)=' "h <ENTER>", or jump to a specific '
      x(17)='HELP section Hn (n=1,2,...17) by typ'
      y(17)='ing "Hn <ENTER>"; to return to the  '
      x(18)='place in the program from which HELP'
      y(18)=' is called, press "r <ENTER>".      '
      do 102 i = 1,18
        write(*,*) x(i),y(i)
  102   continue
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*) '-----------------------------------------------'
      read(*,999) ans,n
      if (ans.eq.'r' .or. ans.eq.'R') return
      GO TO (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17),N
    2 CONTINUE
      write(*,*) 'H2:  File name entry.'
      x(1)='  MAKEPQW is used to create a FORTRA'
      y(1)='N file containing the coefficients  '
      x(2)='p(x),q(x),w(x), defining the differe'
      y(2)='ntial equation, and the boundary    '
      x(3)='condition functions u(x),v(x) if req'
      y(3)='uired.  The file must be given a NEW'
      x(4)='filename which is acceptable to your'
      y(4)=' FORTRAN compiler.  For example, it '
      x(5)='might be called bessel.f or bessel.f'
      y(5)='or depending upon your compiler.    '
      x(6)='  The same naming considerations app'
      y(6)='ly if the FORTRAN file is prepared  '
      x(7)='other than with the use of MAKEPQW. '
      y(7)=' '
      do 201 i = 1,7
        write(*,*) x(i),y(i)
  201   continue
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*) '-----------------------------------------------'
      read(*,999) ans,n
      if (ans.eq.'r' .or. ans.eq.'R') return
      GO TO (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17),N
    3 CONTINUE
      write(*,*) 'H3:  The differential equation.'
      x(1)='  The prompt "Input p (or q or w) ="'
      y(1)=' requests you to type in a FORTRAN  '
      x(2)='expression defining the function p(x'
      y(2)='), which is one of the three coeffi-'
      x(3)='cient functions defining the Sturm-L'
      y(3)='iouville differential equation      '
      x(4)=' '
      y(4)=' '
      x(5)='                   -(p*y'')'' + q*y = '
      y(5)=' lambda*w*y              (*)        '
      x(6)=' '
      y(6)=' '
      x(7)='to be considered on some interval (a'
      y(7)=',b) of the real line.  The actual   '
      x(8)='interval used in a particular proble'
      y(8)='m can be chosen later, and may be   '
      x(9)='either the whole interval (a,b) wher'
      y(9)='e the coefficient functions p,q,w,  '
      x(10)='etc. are defined or any sub-interval'
      y(10)=' (a'',b'') of (a,b); a = -infinity  '
      x(11)='and/or b = +infinity are allowable c'
      y(11)='hoices for the end-points.          '
      x(12)='  The coefficient functions p,q,w of'
      y(12)=' the differential equation may be   '
      x(13)='chosen arbitrarily but must satisfy '
      y(13)='the following conditions:           '
      x(14)='  (1) p,q,w are real-valued througho'
      y(14)='ut (a,b).                           '
      x(15)='  (2) p,q,w are piece-wise continuou'
      y(15)='s and defined throughout the        '
      x(16)='      interior of the interval (a,b)'
      y(16)='.                                   '
      x(17)='  (3) p and w are strictly positive '
      y(17)='in (a,b).                           '
      x(18)='  For better error analysis in the n'
      y(18)='umerical procedures, condition      '
      x(19)='(2) above is often replaced with    '
      y(19)=' '
      x(20)='  (2'') p,q,w are four times continuo'
      y(20)='usly differentiable on (a,b).       '
      x(21)='  The behavior of p,q,w near the end'
      y(21)='-points a and b is critical to the  '
      x(22)='classification of the differential e'
      y(22)='quation (see H4 and H11).           '
      do 301 i = 1,22
        write(*,*) x(i),y(i)
  301   continue
      write(*,*) '-----------------------------------------------'
      read(*,999) ans,n
      if (ans.eq.'r' .or. ans.eq.'R') return
      GO TO (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17),N
    4 CONTINUE
      write(*,*) 'H4:  End-point classification.'
      x(1)='  The correct classification of the '
      y(1)='end-points a and b is essential to  '
      x(2)='the working of the SLEIGN2 program. '
      y(2)=' To classify the end-points, it is  '
      x(3)='convenient to choose a point c in (a'
      y(3)=',b); i.e.,  a < c < b.  Subject to  '
      x(4)='the general conditions on the coeffi'
      y(4)='cient functions p,q,w (see H3):     '
      x(5)='  (1) a is REGULAR (say R) if -infin'
      y(5)='ity < a, p,q,w are piece-wise       '
      x(6)='      continuous on [a,c], and p(a) '
      y(6)='> 0 and w(a) > 0.                   '
      x(7)='  (2) a is WEAKLY REGULAR (say WR) i'
      y(7)='f -infinity < a, a is not R, and    '
      x(8)='               |c                   '
      y(8)=' '
      x(9)='      integral | {1/p+|q|+w} < +infi'
      y(9)='nity.                               '
      x(10)='               |a                   '
      y(10)=' '
      x(11)=' '
      y(11)=' '
      x(12)='    If end-point a is neither R nor '
      y(12)='WR, then a is SINGULAR; that is,    '
      x(13)='  either -infinity = a, or -infinity'
      y(13)=' < a and                            '
      x(14)='               |c                   '
      y(14)=' '
      x(15)='      integral | {1/p+|q|+w} = +infi'
      y(15)='nity.                               '
      x(16)='               |a                   '
      y(16)=' '
      x(17)='  (3) The SINGULAR end-point a is LI'
      y(17)='MIT-CIRCLE NON-OSCILLATORY (say     '
      x(18)='      LCNO) if for some real lambda '
      y(18)='ALL real-valued solutions y of the  '
      x(19)='      differential equation         '
      y(19)=' '
      x(20)=' '
      y(20)=' '
      x(21)='                   -(p*y'')'' + q*y = '
      y(21)=' lambda*w*y on (a,c]     (*)        '
      x(22)=' '
      y(22)=' '
      x(23)='      satisfy the conditions:       '
      y(23)=' '
      do 401 i = 1,23
        write(*,*) x(i),y(i)
  401   continue
      read(*,999) ans,n
      if (ans.eq.'r' .or. ans.eq.'R') return
      GO TO (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17),N
      x(1)='               |c                   '
      y(1)=' '
      x(2)='      integral | { w*y*y } < +infini'
      y(2)='ty, and                             '
      x(3)='               |a                   '
      y(3)=' '
      x(4)='      y has at most a finite number '
      y(4)='of zeros in (a,c].                  '
      x(5)='  (4) The SINGULAR end-point a is LI'
      y(5)='MIT-CIRCLE OSCILLATORY (say LCO) if '
      x(6)='for some real lambda ALL real-valued'
      y(6)=' solutions of the differential      '
      x(7)='equation (*) satisfy the conditions:'
      y(7)=' '
      x(8)='               |c                   '
      y(8)=' '
      x(9)='      integral | { w*y*y } < +infini'
      y(9)='ty, and                             '
      x(10)='               |a                   '
      y(10)=' '
      x(11)='      y has an infinite number of ze'
      y(11)='ros in (a,c].                       '
      x(12)='  (5) The SINGULAR end-point a is LI'
      y(12)='MIT POINT (say LP) if for some real '
      x(13)='lambda at least one solution of the '
      y(13)='differential equation (*) satisfies '
      x(14)='the condition:                      '
      y(14)=' '
      x(15)='               |c                   '
      y(15)=' '
      x(16)='      integral | {w*y*y} = +infinity'
      y(16)='.                                   '
      x(17)='               |a                   '
      y(17)=' '
      x(18)='  There is a similar classification '
      y(18)='of the end-point b into one of the  '
      x(19)='five distinct cases R, WR, LCNO, LCO'
      y(19)=', LP.                               '
      do 402 i = 1,19
        write(*,*) x(i),y(i)
  402   continue
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*)
      read(*,999) ans,n
      if (ans.eq.'r' .or. ans.eq.'R') return
      GO TO (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17),N
      x(1)='  Although the classification of sin'
      y(1)='gular end-points invokes a real     '
      x(2)='value of the parameter lambda, this '
      y(2)='classification is invariant in      '
      x(3)='lambda; all real choices of lambda l'
      y(3)='ead to the same classification.     '
      x(4)='  In determining the classification '
      y(4)='of singular end-points for the      '
      x(5)='differential equation (*), it is oft'
      y(5)='en convenient to start with the     '
      x(6)='choice lambda = 0 in attempting to f'
      y(6)='ind solutions (particularly when    '
      x(7)='q = 0 on (a,b)); however, see exampl'
      y(7)='e 7 below.                          '
      x(8)='  See H6 on the use of maximal domai'
      y(8)='n functions to determine the        '
      x(9)='classification at singular end-point'
      y(9)='s.                                  '
      do 403 i = 1,9
        write(*,*) x(i),y(i)
  403   continue
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*)
      read(*,999) ans,n
      if (ans.eq.'r' .or. ans.eq.'R') return
      GO TO (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17),N
      write(*,*) ' EXAMPLES: '
      x(1)='  1. -y'''' = lambda*y is R at both en'
      y(1)='d-points of (a,b) when a and b are  '
      x(2)='     finite.                        '
      y(2)=' '
      x(3)='  2. -y'''' = lambda*y on (-infinity,i'
      y(3)='nfinity) is LP at both end-points.  '
      x(4)='  3. -(sqrt(x)*y''(x))'' = lambda*(1./'
      y(4)='sqrt(x))*y(x) on (0,infinity) is    '
      x(5)='     WR at 0 and LP at +infinity (ta'
      y(5)='ke lambda = 0 in (*)).  See         '
      x(6)='     examples.f, #10 (Weakly Regular'
      y(6)=').                                  '
      x(7)='  4. -((1-x*x)*y''(x))'' = lambda*y(x)'
      y(7)=' on (-1,1) is LCNO at both ends     '
      x(8)='     (take lambda = 0 in (*)).  See '
      y(8)='xamples.f, #1 (Legendre).           '
      x(9)='  5. -y''''(x) + C*(1/(x*x))*y(x) = la'
      y(9)='mbda*y(x) on (0,infinity) is LP at  '
      x(10)='     infinity and at 0 is (take lamb'
      y(10)='da = 0 in (*)):                     '
      x(11)='       LP for C .ge. 3/4 ;          '
      y(11)=' '
      x(12)='       LCNO for -1/4 .le. C .lt. 3/4'
      y(12)=' (but C .ne. 0);                    '
      x(13)='       LCO for C .lt. -1/4.         '
      y(13)=' '
      x(14)='  6. -(x*y''(x))'' - (1/x)*y(x) = lamb'
      y(14)='da*y(x) on (0,infinity) is LCO at 0 '
      x(15)='     and LP at +infinity (take lambd'
      y(15)='a = 0 in (*) with solutions         '
      x(16)='     cos(ln(x)) and sin(ln(x))).  Se'
      y(16)='e xamples.f, #7 (BEZ).              '
      x(17)='  7. -(x*y''(x))'' - x*y(x) = lambda*('
      y(17)='1/x)*y(x) on (0,infinity) is LP at 0'
      x(18)='     and LCO at infinity (take lambd'
      y(18)='a = -1/4 in (*) with solutions      '
      x(19)='     cos(x)/sqrt(x) and sin(x)/sqrt('
      y(19)='x)).  See xamples.f,                '
      x(20)='     #6 (Sears-Titchmarsh).         '
      y(20)=' '
      do 404 i = 1,20
        write(*,*) x(i),y(i)
  404   continue
      write(*,*)
      write(*,*)
      write(*,*) '-----------------------------------------------'
      read(*,999) ans,n
      if (ans.eq.'r' .or. ans.eq.'R') return
      GO TO (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17),N
    5 CONTINUE
      write(*,*) 'H5:  DEFAULT entry.'
      x(1)='  The complete range of problems for'
      y(1)=' which SLEIGN2 is applicable can    '
      x(2)='only be reached by appropriate entri'
      y(2)='es under end-point classification   '
      x(3)='and boundary conditions.  However, t'
      y(3)='here is a DEFAULT application which '
      x(4)='requires no detailed entry of end-po'
      y(4)='int classification or boundary      '
      x(5)='conditions, subject to:             '
      y(5)=' '
      x(6)='  1) The DEFAULT application CANNOT '
      y(6)='be used at a LCO end-point.         '
      x(7)='  2) If an end-point a is R, then th'
      y(7)='e Dirichlet boundary condition      '
      x(8)='     y(a) = 0 is automatically used.'
      y(8)=' '
      x(9)='  3) If an end-point a is WR, then t'
      y(9)='he following boundary condition     '
      x(10)='     is automatically applied:      '
      y(10)=' '
      x(11)='       if p(a) = 0, and both q(a),w('
      y(11)='a) are bounded, then the Dirichlet  '
      x(12)='       boundary condition y(a) = 0 i'
      y(12)='s used, or                          '
      x(13)='       if p(a) > 0, and q(a) and/or '
      y(13)='w(a)) are not bounded, then the     '
      x(14)='       Neumann boundary condition (p'
      y(14)='y'')(a) = 0 is used.                 '
      x(15)='     If p(a) = 0, and q(a) and/or w('
      y(15)='a) are not bounded, then no reliable'
      x(16)='     information can be given on the'
      y(16)=' DEFAULT boundary condition.        '
      x(17)='  4) If an end-point is LCNO, then i'
      y(17)='n most cases the principal or       '
      x(18)='     Friedrichs boundary condition i'
      y(18)='s applied (see H6).                 '
      x(19)='  5) If an end-point is LP, then the'
      y(19)=' normal LP procedure is applied     '
      x(20)='     (see H7(1.)).                  '
      y(20)=' '
      x(21)='If you choose the DEFAULT condition,'
      y(21)=' then no entry is required for the  '
      x(22)='u,v boundary condition functions.   '
      y(22)=' '
      do 501 i = 1,22
        write(*,*) x(i),y(i)
  501   continue
      write(*,*) '-----------------------------------------------'
      read(*,999) ans,n
      if (ans.eq.'r' .or. ans.eq.'R') return
      GO TO (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17),N
    6 CONTINUE
      write(*,*) 'H6:  Limit-circle boundary conditions.'
      x(1)='  At an end-point a, the limit-circl'
      y(1)='e type separated boundary condition '
      x(2)='is of the form (similar remarks thro'
      y(2)='ughout apply to the end-point b)    '
      x(3)=' '
      y(3)=' '
      x(4)='       A1*[y,u](a) + A2*[y,v](a) = 0'
      y(4)=',                          (**)     '
      x(5)=' '
      y(5)=' '
      x(6)='where y is a solution of the differe'
      y(6)='ntial equation                      '
      x(7)=' '
      y(7)=' '
      x(8)='     -(p*y'')'' + q*y = lambda*w*y  on'
      y(8)=' (a,b).                     (*)     '
      x(9)=' '
      y(9)=' '
      x(10)='Here A1, A2 are real numbers; u and '
      y(10)='v are boundary condition functions; '
      x(11)='and for real-valued y and u the form'
      y(11)=' [y,u] is defined by                '
      x(12)=' '
      y(12)=' '
      x(13)='      [y,u](x) = y(x)*(pu'')(x) - u(x'
      y(13)=')*(py'')(x)   for x in (a,b).        '
      x(14)=' '
      y(14)=' '
      x(15)='  The object of this section is to p'
      y(15)='rovide help in choosing appropriate '
      x(16)='functions u and v in (**), given the'
      y(16)=' differential equation (*).  Full   '
      x(17)='details of the boundary conditions f'
      y(17)='or (*) are discussed in H7; here it '
      x(18)='is sufficient to say that the limit-'
      y(18)='circle type boundary condition (**) '
      x(19)='can be applied at any end-point in t'
      y(19)='he LCNO, LCO classification, but    '
      x(20)='also in the R, WR classification sub'
      y(20)='ject to the appropriate choice of   '
      x(21)='u and v.                            '
      y(21)=' '
      do 601 i = 1,21
        write(*,*) x(i),y(i)
  601   continue
      write(*,*)
      write(*,*)
      read(*,999) ans,n
      if (ans.eq.'r' .or. ans.eq.'R') return
      GO TO (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17),N
      x(1)='  Let (*) be R, WR, LCNO, or LCO at '
      y(1)='end-point a and choose c in (a,b).  '
      x(2)='Then either                         '
      y(2)=' '
      x(3)='    u and v are a pair of linearly i'
      y(3)='ndependent solutions of (*) on (a,c]'
      x(4)='  for any chosen real values of lamb'
      y(4)='da, or                              '
      x(5)='    u and v are a pair of real-value'
      y(5)='d maximal domain functions defined  '
      x(6)='  on (a,c] satisfying [u,v](a) .ne. '
      y(6)='0.  The maximal domain D(a,c] is    '
      x(7)='  defined by                        '
      y(7)=' '
      x(8)=' '
      y(8)=' '
      x(9)='       D(a,c] = {f:(a,c]->R:: f,pf'' '
      y(9)='in AC(a,c];                         '
      x(10)='                   f, ((-pf'')''+qf)/w'
      y(10)=' in L2((a,c;w)}                     '
      x(11)=' '
      y(11)=' '
      x(12)='  It is known that for all f,g in D('
      y(12)='a,c] the limit                      '
      x(13)=' '
      y(13)=' '
      x(14)='       [f,g](a) = lim[f,g](x) as x->'
      y(14)='a                                   '
      x(15)=' '
      y(15)=' '
      x(16)='  exists and is finite.  If (*) is L'
      y(16)='CNO or LCO at a, then all solutions '
      x(17)='  of (*) belong to D(a,c] for all va'
      y(17)='lues of lambda.                     '
      do 602 i = 1,17
        write(*,*) x(i),y(i)
  602   continue
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*)
      read(*,999) ans,n
      if (ans.eq.'r' .or. ans.eq.'R') return
      GO TO (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17),N
      x(1)='  The boundary condition (**) is ess'
      y(1)='ential in the LCNO and LCO cases but'
      x(2)='can also be used with advantage in s'
      y(2)='ome R and WR cases.  In the R, WR, '
      x(3)='and LCNO cases, but not in the LCO c'
      y(3)='ase, the boundary condition         '
      x(4)='functions can always be chosen so th'
      y(4)='at                                  '
      x(5)='        lim u(x)/v(x) = 0 as x->a,  '
      y(5)=' '
      x(6)='and it is recommended that this norm'
      y(6)='alisation be effected; this has been'
      x(7)='done in the examples given below. In'
      y(7)=' this case, the boundary condition  '
      x(8)='[y,u](a) = 0 (i.e., A1 = 1, A2 = 0 i'
      y(8)='n (**)) is called the principal or  '
      x(9)='Friedrichs boundary condition.      '
      y(9)=' '
      x(10)=' '
      y(10)=' '
      x(11)='  In the case when end-points a and '
      y(11)='b are, independently, in R, WR,     '
      x(12)='LCNO, or LCO classification, it may '
      y(12)='be that symmetry or other reasons   '
      x(13)='permit one set of boundary condition'
      y(13)=' functions to be used at both end-  '
      x(14)='points (see xamples.f, #1 (Legendre)'
      y(14)=').  In other cases, different pairs '
      x(15)='must be chosen for each end-point (s'
      y(15)='ee xamples.f: #16 (Jacobi),         '
      x(16)='#18 (Dunsch), and #19 (Donsch)).    '
      y(16)=' '
      x(17)=' '
      y(17)=' '
      x(18)='  Note that a solution pair u,v is a'
      y(18)='lways a maximal domain pair, but not'
      x(19)='necessarily vice versa.             '
      y(19)=' '
      do 603 i = 1,19
        write(*,*) x(i),y(i)
  603   continue
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*)
      read(*,999) ans,n
      if (ans.eq.'r' .or. ans.eq.'R') return
      GO TO (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17),N
      x(1)=' EXAMPLES:                          '
      y(1)=' '
      x(2)='1. -y''''(x) = lambda*y(x) on [0,pi] i'
      y(2)='s R at 0 and R at pi.               '
      x(3)='   At 0, with lambda = 0, a solution'
      y(3)=' pair is u(x) = x, v(x) = 1.        '
      x(4)='   At pi, with lambda = 1, a solutio'
      y(4)='n pair is                           '
      x(5)='     u(x) = sin(x), v(x) = cos(x).  '
      y(5) =' '
      x(6)='2. -(sqrt(x)*y''(x))'' = lambda*y(x)/s'
      y(6)='qrt(x) on (0,1] is                  '
      x(7)='     WR at 0 and R at 1.            '
      y(7)=' '
      x(8)='   (The general solutions of this eq'
      y(8)='uation are                          '
      x(9)='     u(x) = cos(2*sqrt(x*lambda)), v'
      y(9)='(x) = sin(2*sqrt(x*lambda)).)       '
      x(10)='   At 0, with lambda = 0, a solution'
      y(10)=' pair is                            '
      x(11)='     u(x) = 2*sqrt(x), v(x) = 1.    '
      y(11)=' '
      x(12)='   At 1, with lambda = pi*pi/4, a so'
      y(12)='lution pair is                      '
      x(13)='     u(x) = sin(pi*sqrt(x)), v(x) = '
      y(13)='cos(pi*sqrt(x)).                    '
      x(14)='   At 1, with lambda = 0, a solution'
      y(14)=' pair is                            '
      x(15)='     u(x) = 2*(1-sqrt(x)), v(x) = 1.'
      y(15)=' '
      x(16)='   See also xamples.f, #10 (Weakly R'
      y(16)='egular).                            '
      x(17)='3. -((1-x*x)*y''(x))'' = lambda*y(x) o'
      y(17)='n (-1,1) is LCNO at both ends.      '
      x(18)='   At +-1, with lambda = 0, a soluti'
      y(18)='on pair is                          '
      x(19)='     u(x) = 1, v(x) = 0.5*log((1+x)/'
      y(19)='(1-x)).                             '
      x(20)='   At 1, a maximal domain pair is u('
      y(20)='x) = 1, v(x) = log(1-x)             '
      x(21)='   At -1, a maximal domain pair is u'
      y(21)='(x) = 1, v(x) = log(1+x).           '
      x(22)='   See also xamples.f, #1 (Legendre)'
      y(22)='.                                   '
      do 604 i = 1,22
        write(*,*) x(i),y(i)
  604   continue
      write(*,*)
      write(*,*)
      read(*,999) ans,n
      if (ans.eq.'r' .or. ans.eq.'R') return
      GO TO (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17),N
      x(1)='4. -y''''(x) - (1/(4x*x))*y(x) = lambd'
      y(1)='a*y(x) on (0,infinity) is           '
      x(2)='     LCNO at 0 and LP at +infinity. '
      y(2)=' '
      x(3)='   At 0, a maximal domain pair is   '
      y(3)=' '
      x(4)='     u(x) = sqrt(x), v(x) = sqrt(x)*'
      y(4)='log(x).                             '
      x(5)='   See also xamples.f, #2 (Bessel). '
      y(5)=' '
      x(6)='5. -y''''(x) - 5*(1/(4*x*x))*y(x) = la'
      y(6)='mbda*y(x) on (0,infinity) is        '
      x(7)='     LCO at 0 and LP at +infinity.  '
      y(7)=' '
      x(8)='   At 0, with lambda = 0, a solution'
      y(8)=' pair is                            '
      x(9)='     u(x) = sqrt(x)*cos(log(x)), v(x'
      y(9)=') = sqrt(x)*sin(log(x))             '
      x(10)='   See also xamples.f, #20 (Krall). '
      y(10)=' '
      x(11)='6. -y''''(x) - (1/x)*y(x) = lambda*y(x'
      y(11)=') on (0,infinity) is               '
      x(12)='     LCNO at 0 and LP at +infinity.'
      y(12)=' '
      x(13)='   At 0, a maximal domain pair is   '
      y(13)=' '
      x(14)='     u(x) = x, v(x) = 1 -x*log(x).  '
      y(14)=' '
      x(15)='   See also xamples.f, #4(Boyd).    '
      y(15)=' '
      x(16)='7. -((1/x)*y''(x))'' + (k/(x*x) + k*k/'
      y(16)='x)*y(x) = lambda*y(x) on (0,1],     '
      x(17)='     with k real and .ne. 0, is LCNO'
      y(17)=' at 0 and R at 1.                   '
      x(18)='   At 0, a maximal domain pair is   '
      y(18)=' '
      x(19)='     u(x) = x*x, v(x) = x - 1/k.    '
      y(19)=' '
      x(20)='   See also xamples.f, #8 (Laplace T'
      y(20)='idal Wave).                         '
      do 605 i = 1,20
        write(*,*) x(i),y(i)
  605   continue
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*) '-----------------------------------------------'
      read(*,999) ans,n
      if (ans.eq.'r' .or. ans.eq.'R') return
      GO TO (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17),N
    7 CONTINUE
      write(*,*) 'H7:  General boundary conditions.'
      x(1)='  Boundary conditions for Sturm-Liou'
      y(1)='ville boundary value problems       '
      x(2)=' '
      y(2)=' '
      x(3)='                   -(p*y'')'' + q*y = '
      y(3)=' lambda*w*y              (*)        '
      x(4)=' '
      y(4)=' '
      x(5)='on an interval (a,b) are either     '
      y(5)=' '
      x(6)='    SEPARATED, with at most one cond'
      y(6)='ition at end-point a and at most    '
      x(7)='  one condition at end-point b, or  '
      y(7)=' '
      x(8)='    COUPLED, when both a and b are, '
      y(8)='independently, in one of the end-   '
      x(9)='  point classifications R, WR, LCNO,'
      y(9)=' LCO, in which case two independent '
      x(10)='  ent boundary conditions are requir'
      y(10)='ed which link the solution values   '
      x(11)='  near a to those near b.           '
      y(11)=' '
      x(12)='The SLEIGN2 program allows for all s'
      y(12)='eparated conditions; and special    '
      x(13)='cases of the coupled conditions -- t'
      y(13)='he so-called periodic boundary      '
      x(14)='conditions applicable only when the '
      y(14)='interval (a,b) is finite and both   '
      x(15)='a and b are R.                      '
      y(15)=' '
      do 701 i = 1,15
        write(*,*) x(i),y(i)
  701   continue
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*)
      read(*,999) ans,n
      if (ans.eq.'r' .or. ans.eq.'R') return
      GO TO (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17),N
      write(*,*) '    Separated Conditions:      '
      write(*,*) '    ---------------------      '
      x(1)='  The boundary conditions to be sele'
      y(1)='cted depend upon the classification '
      x(2)='of the differential equation at the '
      y(2)='end-point, say, a:                  '
      x(3)='    1. If the end-point a is LP, the'
      y(3)='n no boundary condition is required '
      x(4)='  or allowed.                       '
      y(4)=' '
      x(5)='    2. If the end-point a is R or WR'
      y(5)=', then a separated boundary         '
      x(6)='  condition is of the form          '
      y(6)=' '
      x(7)='         A1*y(a) + A2*(py'')(a) = 0,'
      y(7)=' '
      x(8)='  where A1, A2 are real constants yo'
      y(8)='u must choose, not both zero.       '
      x(9)='    3. If the end-point a is LCNO or'
      y(9)=' LCO, then a separated boundary     '
      x(10)='  condition is of the form          '
      y(10)=' '
      x(11)='         A1*[y,u](a) + A2*[y,v](a) ='
      y(11)=' 0,                                 '
      x(12)='  where A1, A2 are real constants yo'
      y(12)='u must choose, not both zero;       '
      x(13)='  here, u,v are the pair of boundary'
      y(13)=' condition functions you have       '
      x(14)='  previously selected when the input'
      y(14)=' FORTRAN file was being prepared.   '
      x(15)='    4. If the end-point a is LCNO an'
      y(15)='d the boundary condition pair       '
      x(16)='  u,v has been chosen so that       '
      y(16)=' '
      x(17)='         lim u(x)/v(x) = 0  as x->a '
      y(17)=' '
      x(18)='  (which is always possible), then A'
      y(18)='1 = 1, A2 = 0 (i.e., [y,u](a) = 0)  '
      x(19)='  gives the principal (Friedrichs) b'
      y(19)='oundary condition at a.             '
      do 702 i = 1,19
        write(*,*) x(i),y(i)
  702   continue
      write(*,*)
      write(*,*)
      write(*,*)
      read(*,999) ans,n
      if (ans.eq.'r' .or. ans.eq.'R') return
      GO TO (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17),N
      x(1)='    5. If a is R or WR and boundary '
      y(1)='condition functions u,v have been   '
      x(2)='  entered in the FORTRAN input file,'
      y(2)=' then (3.,4.) above apply to        '
      x(3)='  entering separated boundary condit'
      y(3)='ions at such an end-point; the      '
      x(4)='  boundary conditions in this form a'
      y(4)='re equivalent to the point-wise     '
      x(5)='  conditions in (2.) (subject to car'
      y(5)='e in choosing A1, A2).  This        '
      x(6)='  singular form of a regular boundar'
      y(6)='y condition may be particularly     '
      x(7)='  effective in the WR case if the bo'
      y(7)='undary condition form in (2.) leads '
      x(8)='  to numerical difficulties.        '
      y(8)=' '
      x(9)=' '
      y(9)=' '
      x(10)='  Conditions (2.,3.,4.,5.) apply sim'
      y(10)='ilarly at end-point b.              '
      x(11)=' '
      y(11)=' '
      x(12)='    6. If a is R, WR, LCNO, or LCO a'
      y(12)='nd b is LP, then only a separated   '
      x(13)='  condition at a is required and all'
      y(13)='owed (or instead at b if a and b    '
      x(14)='  are interchanged).                '
      y(14)=' '
      x(15)='    7. If both end-points a and b ar'
      y(15)='e LP, then no boundary conditions   '
      x(16)='  are required or allowed.          '
      y(16)=' '
      x(17)=' '
      y(17)=' '
      x(18)='  The indexing of eigenvalues for bo'
      y(18)='undary value problems with separated'
      x(19)='  conditions is discussed in H13.   '
      y(19)=' '
      do 703 i = 1,19
        write(*,*) x(i),y(i)
  703   continue
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*)
      read(*,999) ans,n
      if (ans.eq.'r' .or. ans.eq.'R') return
      GO TO (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17),n
      write(*,*) '    Coupled Conditions:        '
      write(*,*) '    -------------------        '
      x(1)='    8. Periodic-type boundary condit'
      y(1)='ions on (a,b) apply only when       '
      x(2)='  both end-points a and b are R; the'
      y(2)='se conditions are of the form       '
      x(3)=' '
      y(3)=' '
      x(4)='        y(b) = c*y(a),  (py'')(b) = ('
      y(4)='py'')(a)/c,                          '
      x(5)=' '
      y(5)=' '
      x(6)='  where c may be chosen to be any re'
      y(6)='al number not equal to 0.  The case '
      x(7)='  c = 1 is called periodic, the case'
      y(7)=' c = -1 is called semi-periodic.    '
      x(8)=' '
      y(8)=' '
      x(9)='  The indexing of eigenvalues for pe'
      y(9)='riodic-type boundary conditions is  '
      x(10)='  discussed in H17.                 '
      y(10)=' '
      do 704 i = 1,10
        write(*,*) x(i),y(i)
  704   continue
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*) '-----------------------------------------------'
      read(*,999) ans,n
      if (ans.eq.'r' .or. ans.eq.'R') return
      GO TO (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17),N
    8 CONTINUE
      write(*,*) 'H8:  Recording the results.'
      x(1)='  If you choose to have a record kep'
      y(1)='t of the results, then the following'
      x(2)='information is stored in a file with'
      y(2)=' the name you select:               '
      x(3)=' '
      y(3)=' '
      x(4)='  1. The file name.                 '
      y(4)=' '
      x(5)='  2. The header line prompted for (u'
      y(5)='p to 32 characters of your choice). '
      x(6)='  3. The interval (a,b) which was us'
      y(6)='ed.                                 '
      x(7)=' '
      y(7)=' '
      x(8)='  For SEPARATED boundary conditions:'
      y(8)=' '
      x(9)='  4. The end-point classification.  '
      y(9)=' '
      x(10)='  5. A summary of coefficient inform'
      y(10)='ation at WR, LCNO, LCO end-points.  '
      x(11)='  6. The boundary condition constant'
      y(11)='s (A1,A2), (B1,B2) if entered.      '
      x(12)='  7. (NUMEIG,EIG,TOL) or (NUMEIG1,NU'
      y(12)='MEIG2,TOL), as entered.             '
      x(13)=' '
      y(13)=' '
      x(14)='  For COUPLED boundary conditions:  '
      y(14)=' '
      x(15)='  8. The boundary condition paramete'
      y(15)='r, c.                               '
      x(16)=' '
      y(16)=' '
      x(17)='  For ALL boundary conditions:      '
      y(17)=' '
      x(18)='  9. The computed eigenvalue, EIG, a'
      y(18)='nd its estimated accuracy, TOL.     '
      x(19)=' 10. IFLAG reported (see H15).      '
      y(19)=' '
      do 801 i = 1,19
        write(*,*) x(i),y(i)
  801   continue
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*) '-----------------------------------------------'
      read(*,999) ans,n
      if (ans.eq.'r' .or. ans.eq.'R') return
      GO TO (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17),N
    9 CONTINUE
      write(*,*) 'H9:  Type and choice of interval.'
      x(1)='  You may enter any interval (a,b) f'
      y(1)='or which the coefficients p,q,w are '
      x(2)='well defined by your FORTRAN stateme'
      y(2)='nts in the input file, provided that'
      x(3)='(a,b) contains no interior singulari'
      y(3)='ities.                              '
      do 901 i = 1,3
        write(*,*) x(i),y(i)
  901   continue
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*) '-----------------------------------------------'
      read(*,999) ans,n
      if (ans.eq.'r' .or. ans.eq.'R') return
      GO TO (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17),N
   10 CONTINUE
      write(*,*) 'H10:  Entry of end-points.'
      x(1)='  End-points a and b must be entered'
      y(1)=' as real numbers.  There is no      '
      x(2)='symbolic entry; e.g., pi must be ent'
      y(2)='ered as 3.14159... to an appropriate'
      x(3)='number of decimal places.           '
      y(3)=' '
      do 1001 i = 1,3
        write(*,*) x(i),y(i)
 1001   continue
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*) '-----------------------------------------------'
      read(*,999) ans,n
      if (ans.eq.'r' .or. ans.eq.'R') return
      GO TO (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17),N
   11 CONTINUE
      write(*,*) 'H11:  End-point values of p,q,w.'
      x(1)='  The program SLEIGN2 needs to know '
      y(1)='whether the coefficient functions   '
      x(2)='p(x),q(x),w(x) defined by the FORTRA'
      y(2)='N expressions entered in the input  '
      x(3)='file can be evaluated numerically wi'
      y(3)='thout running into difficulty.  If, '
      x(4)='for example, either q or w is unboun'
      y(4)='ded at a, or p(a) is 0, then SLEIGN2'
      x(5)='needs to know this so that a is not '
      y(5)='chosen for functional evaluation.   '
      do 1101 i = 1,5
        write(*,*) x(i),y(i)
 1101   continue
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*) '-----------------------------------------------'
      read(*,999) ans,n
      if (ans.eq.'r' .or. ans.eq.'R') return
      GO TO (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17),N
   12 CONTINUE
      write(*,*) 'H12:  Initial value problems.'
      x(1)='  The initial value problem facility'
      y(1)=' for Sturm-Liouville problems       '
      x(2)=' '
      y(2)=' '
      x(3)='                   -(p*y'')'' + q*y = '
      y(3)=' lambda*w*y              (*)        '
      x(4)=' '
      y(4)=' '
      x(5)='allows for the computation of a solu'
      y(5)='tion of (*) with a user-chosen      '
      x(6)='value lambda and any one of the foll'
      y(6)='owing initial conditions:           '
      x(7)='  1. From end-point a of any classif'
      y(7)='ication except LP towards           '
      x(8)='end-point b of any classification,  '
      y(8)=' '
      x(9)='  2. From end-point b of any classif'
      y(9)='ication except LP back towards      '
      x(10)='end-point a of any classification,  '
      y(10)=' '
      x(11)='  3. From end-points a and b of any '
      y(11)='classifications except LP towards an'
      x(12)='interior point of (a,b) selected by '
      y(12)='the program.                        '
      x(13)=' '
      y(13)=' '
      x(14)='  Initial values at a are of the for'
      y(14)='m y(a) = alpha1, (p*y'')a = alpha2, '
      x(15)='when a is R or WR; and [y,u](a) = al'
      y(15)='pha1, [y,v](a) = alpha2, when a is  '
      x(16)='LCNO or LCO.                        '
      y(16)=' '
      x(17)='  Initial values at b are of the for'
      y(17)='m y(b) = beta1, (p*y'')b = beta2,   '
      x(18)='when b is R or WR; and [y,u](b) = be'
      y(18)='ta1, [y,v](b) = beta2, when b is  '
      x(19)='LCNO or LCO.                        '
      y(19)=' '
      x(20)='  In (*), lambda is a user-chosen re'
      y(20)='al number; while in the above       '
      x(21)='initial values, (alpha1,alpha2) and '
      y(21)='(beta1,beta2) are user-chosen pairs '
      x(22)='of real numbers not both zero.      '
      y(22)=' '
      do 1201 i = 1,22
        write(*,*) x(i),y(i)
 1201   continue
      write(*,*)
      read(*,999) ans,n
      if (ans.eq.'r' .or. ans.eq.'R') return
      GO TO (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17),N
      x(1)='  In the initial value case (3.) abo'
      y(1)='ve when the interval (a,b) is       '
      x(2)='finite, the interior point selected '
      y(2)='by the program is generally near the'
      x(3)='midpoint of (a,b); when (a,b) is inf'
      y(3)='inite, no general rule can be given.'
      x(4)='Also if, given (alpha1,alpha2) and ('
      y(4)='beta1,beta2), the lambda chosen is  '
      x(5)='an eigenvalue of the associated boun'
      y(5)='dary value problem, the computed    '
      x(6)='solution may not be the correspondin'
      y(6)='g eigenfunction -- the signs of the '
      x(7)='computed solutions on either side of'
      y(7)=' the interior point may be opposite.'
      x(8)='  The output for a solution of an in'
      y(8)='itial value problem is in the form  '
      x(9)='of stored numerical data which can b'
      y(9)='e plotted on the screen (see H16),  '
      x(10)='or printed out in graphical form if '
      y(10)='graphics software is available.     '
      do 1202 i = 1,10
        write(*,*) x(i),y(i)
 1202   continue
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*) '-----------------------------------------------'
      read(*,999) ans,n
      if (ans.eq.'r' .or. ans.eq.'R') return
      GO TO (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17),N
   13 continue
      write(*,*) 'H13:  Indexing of eigenvalues.'
      x(1)='  The indexing of eigenvalues is an '
      y(1)='automatic facility in SLEIGN2.  The '
      x(2)='following general results hold for t'
      y(2)='he separated boundary condition     '
      x(3)='problem (see H7):                   '
      y(3)=' '
      x(4)='  1. If neither end-point a or b is '
      y(4)='LP or LCO, then the spectrum of the '
      x(5)='eigenvalue problem is discrete (eige'
      y(5)='nvalues only), simple (eigenvalues  '
      x(6)='all of multiplicity 1), and bounded '
      y(6)='below with a single cluster point at'
      x(7)='+infinity.  The eigenvalues are inde'
      y(7)='xed as {lambda(n): n=0,1,2,...},    '
      x(8)='where lambda(n) < lambda(n+1) (n=0,1'
      y(8)=',2,...), lim lambda(n) -> +infinity;'
      x(9)='and if {psi(n): n=0,1,2,...} are the'
      y(9)=' corresponding eigenfunctions, then '
      x(10)='psi(n) has exactly n zeros in the op'
      y(10)='en interval (a,b).                 '
      x(11)='  2. If neither end-point a or b is '
      y(11)='LP but at least one end-point is    '
      x(12)='LCO, then the spectrum is discrete a'
      y(12)='nd simple as for (1.), but with     '
      x(13)='cluster points at both +infinity and'
      y(13)=' -infinity.  The eigenvalues are    '
      x(14)='indexed as {lambda(n): n=0,1,-1,2,-2'
      y(14)=',...}, where                        '
      x(15)='lambda(n) < lambda(n+1) (n=...-2,-1,'
      y(15)=',0,1,2,...) with lambda(0) the      '
      x(16)='smallest non-negative eigenvalue and'
      y(16)=' lim lambda(n) -> +infinity or      '
      x(17)='-> -infinity with n; and if {psi(n):'
      y(17)=' n=0,1,-1,2,-2,...} are the         '
      x(18)='corresponding eigenfunctions, then e'
      y(18)='very psi(n) has infinitely many     '
      x(19)='zeros in (a,b).                     '
      y(19)=' '
      do 1301 i = 1,19
        write(*,*) x(i),y(i)
 1301   continue
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*)
      read(*,999) ans,n
      if (ans.eq.'r' .or. ans.eq.'R') return
      GO TO (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17),N
      x(1)='  3. If one or both end-points is LP'
      y(1)=', then there can be one or more     '
      x(2)='intervals of continuous spectrum for'
      y(2)=' the boundary value problem in      '
      x(3)='addition to some (necessarily simple'
      y(3)=') eigenvalues.  For these           '
      x(4)='essentially more difficult problems,'
      y(4)=' SLEIGN2 can be used as an          '
      x(5)='investigative tool to give qualitati'
      y(5)='ve and possibly quantitative        '
      x(6)='information on the spectrum.        '
      y(6)=' '
      x(7)='     For example, if a problem has a'
      y(7)=' single interval of continuous      '
      x(8)='spectrum bounded below by K, then th'
      y(8)='ere may be any number of eigenvalues'
      x(9)='below K.  In some cases, SLEIGN2 can'
      y(9)=' compute K, and determine the number'
      x(10)='of these eigenvalues and compute the'
      y(10)='m.  In this respect, see xamples.f: '
      x(11)='#13 (Hydrogen Atom), #17 (Morse Osci'
      y(11)='llator), #21 (Fourier), and         '
      x(12)='#27 (Joergens) as examples of succes'
      y(12)='s; and #2 (Mathieu), #14 (Marletta),'
      x(13)='and #28 (Behnke-Goerisch) as example'
      y(13)='s of failure.                       '
      x(14)='     The problem need not have a con'
      y(14)='tinuous spectrum, in which case if  '
      x(15)='its discrete spectrum is bounded bel'
      y(15)='ow, then the eigenvalues are indexed'
      x(16)='and the eigenfunctions have zero cou'
      y(16)='nts as in (1.).  If, on the other   '
      x(17)='hand, the discrete spectrum is unbou'
      y(17)='nded below, then all the            '
      x(18)='eigenfunctions have infinitely many '
      y(18)='zeros in the interval.              '
      do 1302 i = 1,18
        write(*,*) x(i),y(i)
 1302   continue
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*)
      read(*,999) ans,n
      if (ans.eq.'r' .or. ans.eq.'R') return
      GO TO (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17),N
      x(1)='  In respect to the three classes of'
      y(1)=' end-points, the following          '
      x(2)='identified examples from xamples.f i'
      y(2)='llustrate the spectral property     '
      x(3)='of these boundary value problems:   '
      y(3)=' '
      x(4)='  1. Neither end-point is LP or LCO.'
      y(4)=' '
      x(5)='       #1 (Legendre)                '
      y(5)=' '
      x(6)='       #2 (Bessel) with -1/4 < c < 3'
      y(6)='/4                                  '
      x(7)='       #4 (Boyd)                    '
      y(7)=' '
      x(8)='       #5 (Latzko)                  '
      y(8)=' '
      x(9)='  2. Neither end-point is LP, but at'
      y(9)=' least one is LCO.                  '
      x(10)='       #6 (Sears-Titchmarsh)        '
      y(10)=' '
      x(11)='       #7 (BEZ)                     '
      y(11)=' '
      x(12)='      #19 (Donsch)                  '
      y(12)=' '
      x(13)='  3. At least one end-point is LP.  '
      y(13)=' '
      x(14)='      #13 (Hydrogen Atom)           '
      y(14)=' '
      x(15)='      #14 (Marletta)                '
      y(15)=' '
      x(16)='      #20 (Krall)                   '
      y(16)=' '
      x(17)='      #21 (Fourier) on [0,infinity) '
      y(17)=' '
      do 1303 i = 1,17
        write(*,*) x(i),y(i)
 1303   continue
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*) '-----------------------------------------------'
      read(*,999) ans,n
      if (ans.eq.'r' .or. ans.eq.'R') return
      GO TO (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17),N
   14 CONTINUE
      write(*,*) 'H14:  Entry of eigenvalue index, initial guess,'//
     1           ' and tolerance.'
      x(1)='  For SEPARATED boundary condition p'
      y(1)='roblems (see H7), SLEIGN2 calls for '
      x(2)='input information options to compute'
      y(2)=' either                             '
      x(3)='  1. a single eigenvalue, or        '
      y(3)=' '
      x(4)='  2. a series of eigenvalues.       '
      y(4)=' '
      x(5)='In each case indexing of eigenvalues'
      y(5)=' is called for (see H13).           '
      x(6)='  (1.) above asks for data triples N'
      y(6)='UMEIG, EIG, TOL separated by commas.'
      x(7)='Here NUMEIG is the integer index of '
      y(7)='the desired eigenvalue; NUMEIG can  '
      x(8)='be negative only when the problem is'
      y(8)=' LCO at one or both end-points.     '
      x(9)='EIG allows for the entry of an initi'
      y(9)='al guess for the requested          '
      x(10)='eigenvalue (if an especially good on'
      y(10)='e is available), or can be set to 0 '
      x(11)='in which case an initial guess is ge'
      y(11)='nerated by SLEIGN2 itself.          '
      x(12)='TOL is the desired accuracy of the c'
      y(12)='omputed eigenvalue.  It is an       '
      x(13)='absolute accuracy if the magnitude o'
      y(13)='f the eigenvalue is 1 or less, and  '
      x(14)='is a relative accuracy otherwise.  T'
      y(14)='ypical values might be .001 for     '
      x(15)='moderate accuracy and .0000001 for h'
      y(15)='igh accuracy in single precision.   '
      x(16)='If TOL is set to 0, the maximum achi'
      y(16)='evable accuracy is requested.       '
      x(17)='  If the input data list is truncate'
      y(17)='d with a "/" after NUMEIG or EIG,   '
      x(18)='then the remaining elements default '
      y(18)='to 0.                               '
      do 1401 i = 1,18
        write(*,*) x(i),y(i)
 1401   continue
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*)
      read(*,999) ans,n
      if (ans.eq.'r' .or. ans.eq.'R') return
      GO TO (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17),N
      x(1)='  (2.) above asks for data triples N'
      y(1)='UMEIG1, NUMEIG2, TOL separated by   '
      x(2)='commas.  Here NUMEIG1 and NUMEIG2 ar'
      y(2)='e the first and last integer indices'
      x(3)='of the sequence of desired eigenvalu'
      y(3)='es, NUMEIG1 < NUMEIG2; they can be  '
      x(4)='negative only when the problem is LC'
      y(4)='O at one or both end-points.        '
      x(5)='TOL is the desired accuracy of the c'
      y(5)='omputed eigenvalues.  It is an      '
      x(6)='absolute accuracy if the magnitude o'
      y(6)='f an eigenvalue is 1 or less, and   '
      x(7)='is a relative accuracy otherwise.  T'
      y(7)='ypical values might be .001 for     '
      x(8)='moderate accuracy and .0000001 for h'
      y(8)='igh accuracy in single precision.   '
      x(9)='If TOL is set to 0, the maximum achi'
      y(9)='evable accuracy is requested.       '
      x(10)='  If the input data list is truncate'
      y(10)='d with a "/" after NUMEIG2, then TOL'
      x(11)='defaults to 0.                      '
      y(11)=' '
      x(12)=' '
      y(12)=' '
      x(13)='  For COUPLED periodic-type boundary'
      y(13)=' condition problems (see H7 and     '
      x(14)='H17), SLEIGN2 asks only for NUMEIG, '
      y(14)='the non-negative integer index of   '
      x(15)='the desired eigenvalue; TOL is set i'
      y(15)='nternally.                          '
      do 1402 i = 1,15
        write(*,*) x(i),y(i)
 1402   continue
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*) '-----------------------------------------------'
      read(*,999) ans,n
      if (ans.eq.'r' .or. ans.eq.'R') return
      GO TO (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17),N
   15 CONTINUE
      write(*,*) 'H15:  IFLAG information.'
      x(1)='  All results are reported by SLEIGN'
      y(1)='2 with a flag identification.  There'
      x(2)='are four values of IFLAG:           '
      y(2)=' '
      x(3)=' '
      y(3)=' '
      x(4)='  1 - The computed eigenvalue has an'
      y(4)=' estimated accuracy within the      '
      x(5)='      tolerance requested.          '
      y(5)=' '
      x(6)=' '
      y(6)=' '
      x(7)='  2 - The computed eigenvalue does n'
      y(7)='ot have an estimated accuracy within'
      x(8)='      the tolerance requested, but i'
      y(8)='s the best the program could obtain.'
      x(9)=' '
      y(9)=' '
      x(10)='  3 - There seems to be no eigenvalu'
      y(10)='e of index equal to NUMEIG.         '
      x(11)=' '
      y(11)=' '
      x(12)='  4 - The program has been unable to'
      y(12)=' compute the requested eigenvalue.  '
      do 1501 i = 1,12
        write(*,*) x(i),y(i)
 1501   continue
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*) '-----------------------------------------------'
      read(*,999) ans,n
      if (ans.eq.'r' .or. ans.eq.'R') return
      GO TO (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17),N
   16 CONTINUE
      write(*,*) 'H16:  Plotting.'
      x(1)='  After computing a single eigenvalu'
      y(1)='e (see H14(1.)), but not a sequence '
      x(2)='of eigenvalues (see H14(2.)), the ei'
      y(2)='genfunction can be plotted.  If this'
      x(3)='is desired, respond "y" when asked s'
      y(3)='o that SLEIGN2 will compute some    '
      x(4)='eigenfunction data and store them.  '
      y(4)=' '
      x(5)='  One can ask that the eigenfunction'
      y(5)=' data be in the form of either      '
      x(6)='points (x,y) for x in (a,b), or poin'
      y(6)='ts (t,y) for t in the standardized  '
      x(7)='interval (-1,1) mapped onto from (a,'
      y(7)='b); the t- choice can be especially '
      x(8)='helpful when the original interval i'
      y(8)='s infinite.  Additionally, one can  '
      x(9)='ask for a plot of the so-called Pruf'
      y(9)='er angle, in x- or t- variables.    '
      x(10)='  In both forms, once the choice has'
      y(10)=' been made of the function to be    '
      x(11)='plotted, a crude plot is displayed o'
      y(11)='n the monitor screen and you are    '
      x(12)='asked whether you wish to save the c'
      y(12)='omputed plot points in a file. '
      do 1601 i = 1,12
        write(*,*) x(i),y(i)
 1601   continue
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*) '-----------------------------------------------'
      read(*,999) ans,n
      if (ans.eq.'r' .or. ans.eq.'R') return
      GO TO (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17),N
   17 CONTINUE
      write(*,*) 'H17:  Indexing of eigenvalues for'//
     1           ' periodic-type problems.'
      x(1)='  The indexing of eigenvalues is an '
      y(1)='automatic facility in SLEIGN2.  The '
      x(2)='following general result holds for t'
      y(2)='he periodic-type boundary condition '
      x(3)='problem (see H7):                   '
      y(3)=' '
      x(4)='  The spectrum of the eigenvalue pro'
      y(4)='blem is discrete (eigenvalues only),'
      x(5)='and bounded below with a single clus'
      y(5)='ter point at +infinity.  In general,'
      x(6)='the spectrum is not simple, but no e'
      y(6)='igenvalue exceeds multiplicity 2.   '
      x(7)='The eigenvalues are indexed as {lamb'
      y(7)='da(n): n=0,1,2,...}, where          '
      x(8)='lambda(n) .le. lambda(n+1) (n=0,1,2,'
      y(8)='...), lim lambda(n) -> +infinity.   '
      x(9)='  The connection between the index n'
      y(9)=' and the number of zeros of the     '
      x(10)='corresponding eigenfunction psi(x,n)'
      y(10)=' is not as simple as for separated  '
      x(11)='conditions; that is, psi(x,n) need n'
      y(11)='ot have exactly n zeros in (a,b).   '
      x(12)='  The following identified examples '
      y(12)='from xamples.f are of special       '
      x(13)='interest:                           '
      y(13)=' '
      x(14)='    #11 (Plum) on [0,pi]            '
      y(14)=' '
      x(15)='    #21 (Fourier) on [0,pi]         '
      y(15)=' '
      x(16)='    #25 (Meissner) on [-0.5,0.5]    '
      y(16)=' '
      do 1701 i = 1,16
        write(*,*) x(i),y(i)
 1701   continue
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*) '-----------------------------------------------'
      read(*,999) ans,n
      if (ans.eq.'r' .or. ans.eq.'R') return
      GO TO (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17),N
      go to 1
  999 FORMAT(A1,I2)
      end
      SUBROUTINE MESH(EIG,TS,NN)
      INTEGER NN
      DOUBLE PRECISION EIG
      DOUBLE PRECISION TS(991)
C     **********
C     **********
C     .. Scalars in Common ..
      DOUBLE PRECISION PI,TWOPI,HPI
C     ..
C     .. Local Scalars ..
      INTEGER I,IJ,IMAX,IMAXP,J,JJ,JM,JP,M
      PARAMETER (IMAX = 100)
      PARAMETER (IMAXP = IMAX + 1)
      DOUBLE PRECISION DELL,DELL0,DEN,DT,FACTOR,PX,QQ,QX,
     1     SAV,SUM,S1,T1,WX,X
C     ..
C     .. Local Arrays ..
      INTEGER JL(5)
      DOUBLE PRECISION TT(IMAXP),DELT(IMAXP),T(1000),S(1000)
C     ..
C     .. External Subroutines ..
      EXTERNAL DXDT
C     ..
C     .. External Functions ..
      DOUBLE PRECISION P,Q,W
      EXTERNAL P,Q,W
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC MIN
C     ..
C     .. Common blocks ..
      COMMON /PIE/PI,TWOPI,HPI
C     ..
      TT(1) = -1.0D0
      SAV = -1.0D0
      M = 0
      DO 10 I = 2,IMAX
         TT(I) = -1.0D0 + (2.0*(I-1))/IMAX
         CALL DXDT(TT(I),DT,X)
         PX = P(X)
         QX = Q(X)
         WX = W(X)
         QQ = DT**2*(EIG*WX-QX)/PX
         DEN = 3.0*HPI
         IF (QQ.GT.HPI**2) DEN = 3.0*SQRT(QQ)
         DELT(I) = 1.0D0/DEN
         IF (QQ*SAV.LT.0.0D0) THEN
            M = M + 1
            JL(M) = I
            IF (QQ.LT.0.0D0) JL(M) = -I
            SAV = -SAV
            END IF
   10    CONTINUE
      TT(IMAXP) = 1.0D0
      DELT(1) = DELT(2)
      DELT(IMAXP) = DELT(IMAX)
      DO 30 I = 1,M
         DO 20 J = 1,8
            IJ = J + JL(I) - 6
            IF (JL(I).LT.0) IJ = J - JL(I) - 4
            IF (IJ.GE.1.AND.IJ.LE.IMAXP) DELT(IJ) = MIN(.033D0,DELT(IJ))
   20       CONTINUE
   30    CONTINUE
      SUM = 0.0D0
      DO 40 I = 1,IMAXP
         SUM = SUM + DELT(I)
   40    CONTINUE
      FACTOR = 5.0/SUM
C     ------------------------------------------------------
C     GENERATE MESH POINTS IN (0,1):
C
      J = 1
      T(J) = .001
      DO 50 I = 1,IMAX
         IF (TT(I).LE.T(J) .AND. TT(I+1).GE.T(J)) THEN
            JJ = I
            GO TO 60
            END IF
   50    CONTINUE
   60 CONTINUE
      T1 = T(1)
   70 CONTINUE
         DELL0 = MIN(DELT(JJ),DELT(JJ+1))
         DELL = FACTOR*DELL0
   80    CONTINUE
            IF (T1.LE.TT(JJ+1)) T1 = T1 + DELL
            IF (T1.LE.TT(JJ+1)) THEN
               J = J + 1
               T(J) = T1
               IF (J.LE.495) GO TO 80
            ELSE IF (T1.GT.TT(JJ+1) .AND. T1.LE.TT(JJ+2)) THEN
               J = J + 1
               T(J) = T1
               END IF
         JJ = JJ + 1
         IF (JJ.LT.IMAXP .AND. J.LE.495) GO TO 70
      JP = J
C     ------------------------------------------------------
C     GENERATE MESH POINTS IN (-1,0):
C
      J = 1
      S(J) = -.001
      DO 90 I = 1,IMAX
         IF (TT(I).LE.S(J) .AND. TT(I+1).GE.S(J)) THEN
            JJ = I + 1
            GO TO 100
            END IF
   90    CONTINUE
  100 CONTINUE
      S1 = S(J)
  110 CONTINUE
         DELL0 = MIN(DELT(JJ-1),DELT(JJ))
         DELL = FACTOR*DELL0
  120    CONTINUE
            IF (S1.GE.TT(JJ-1)) S1 = S1 - DELL
            IF (S1.GE.TT(JJ-1)) THEN
               J = J + 1
               S(J) = S1
               IF (J.LE.495) GO TO 120
            ELSE IF (S1.LT.TT(JJ-1) .AND. S1.GE.TT(JJ-2)) THEN
               J = J + 1
               S(J) = S1
               END IF
         JJ = JJ - 1
         IF (JJ.GT.1 .AND. J.LE.495) GO TO 110
      JM = J
C     ------------------------------------------------------
      DO 130 I = 1,JM
         TS(I) = S(JM+1-I)
 130     CONTINUE
      DO 140 I = 1,JP
         TS(JM+I) = T(I)
 140     CONTINUE
C
      NN = JM + JP
      RETURN
      END
      SUBROUTINE QPLOT(ISLFUN,XT,NV,PLOTF,NF)
      INTEGER ISLFUN,NV,NF
      DOUBLE PRECISION PLOTF(1000,6),XT(1000,2)
C     **********
C     THIS PROGRAM TRIES TO DRAW GRAPHS.
C     **********
C     .. Local Scalars ..
      INTEGER I,II,IZ,J,K,L,MMAX,NMAX
      PARAMETER (NMAX = 75, MMAX = 22)
      DOUBLE PRECISION DZ,REM,X,XK,XKP,XMAX,XMIN,Y,YK,YKP,YMAX,YMIN
C     ..
C     .. Local Arrays ..
      DOUBLE PRECISION A(1000,2)
      CHARACTER*1 AX(NMAX,MMAX)
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS,INT,MAX,MIN
C     ..
      XMAX = -1000000.
      XMIN = 1000000.
      YMAX = -1000000.
      YMIN = 1000000.
      DZ = YMIN
      DO 10 I = 1,ISLFUN
         X = XT(9+I,NV)
         Y = PLOTF(9+I,NF)
         XMAX = MAX(XMAX,X)
         XMIN = MIN(XMIN,X)
         YMAX = MAX(YMAX,Y)
         YMIN = MIN(YMIN,Y)
         A(I,1) = X
         A(I,2) = Y
         IF (ABS(Y).LT.DZ) THEN
            DZ = ABS(Y)
            IZ = I
            END IF
   10    CONTINUE
C
      IF (YMIN*YMAX.LE.0.0D0) THEN
         Y = MAX(1.0D0,YMAX-YMIN)
      ELSE
         Y = MAX(1.0D0,ABS(YMIN),ABS(YMAX))
         END IF
      DO 20 I = 1,ISLFUN
         A(I,2) = A(I,2)/Y
   20    CONTINUE
      YMAX = YMAX/Y
      YMIN = YMIN/Y
C
      DO 30 I = 1,ISLFUN
         A(I,1) = A(I,1) - XMIN
         IF (YMIN*YMAX.LE.0.0D0) A(I,2) = A(I,2) - YMIN
   30    CONTINUE
C
C     NOW MIN(X) = 0. AND MIN(Y) = 0.
C
      X = XMAX - XMIN
      DO 40 I = 1,ISLFUN
         A(I,1) = NMAX*A(I,1)/X
         A(I,2) = MMAX*A(I,2)
   40    CONTINUE
C
      DO 60 J = 1,NMAX
         DO 50 K = 1,MMAX
            AX(J,K) = ' '
   50       CONTINUE
   60    CONTINUE
C
      DO 80 J = 2,NMAX
         II = 0
         X = J - 0.5
         DO 70 I = 1,ISLFUN
            IF (A(I,1).LE.X) II = I
   70       CONTINUE
C
C     LINE PK,PKP IS: Y-YK = (X-XK)*(YKP-YK)/(XKP-XK)
C     THIS LINE MEETS THE LINE X = J - 0.5 WHERE:
C
         XK = A(II,1)
         XKP = A(II+1,1)
         YK = A(II,2)
         YKP = A(II+1,2)
         Y = YK + (X-XK)*(YKP-YK)/(XKP-XK)
C
         K = MMAX - INT(Y)
         REM = Y + (K-MMAX)
         IF (REM.LE.0.25) THEN
            AX(J,K) = '_'
         ELSE IF (REM.LE.0.50) THEN
            AX(J,K) = '.'
         ELSE IF (REM.LE.0.75) THEN
            AX(J,K) = '+'
         ELSE
            AX(J,K) = '"'
            END IF
   80    CONTINUE
C
      IF (YMIN*YMAX.LT.0.0D0) THEN
         L = INT(A(IZ,2))
      ELSE IF (YMAX.GT.0.0D0) THEN
         L = 0
      ELSE
         L = 22
         END IF
      K = MMAX - L
      DO 90 J = 1,NMAX
         IF (AX(J,K).EQ.' ') AX(J,K) = '.'
   90    CONTINUE
      WRITE(*,*)
      DO 100 K = 1,MMAX
         WRITE(*,'(1X,80A1)') (AX(J,K),J=1,NMAX)
  100    CONTINUE
      RETURN
      END
      SUBROUTINE LSTDIR(CHANS,I,ICOL)
      INTEGER I,ICOL(2)
      CHARACTER*32 CHANS
C     .. Local Scalars ..
      INTEGER J
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
   10    CONTINUE
      RETURN
      END
      CHARACTER*32 FUNCTION FMT2(I1)
      INTEGER I1
C     .. Local Arrays ..
      CHARACTER*2 COL(32)
C     ..
      DATA COL/'01','02','03','04','05','06','07','08','09','10','11',
     1         '12','13','14','15','16','17','18','19','20','21','22',
     2         '23','24','25','26','27','28','29','30','31','32'/
C
      FMT2 = '(F'//COL(I1)//'.0,1X,F'//COL(31-I1)//'.0)'
      RETURN
      END
