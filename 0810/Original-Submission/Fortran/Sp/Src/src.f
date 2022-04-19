c-------------------------------------------------------------------
C  SUBROUTINE SLEIGN

C     **********
C     MARCH 1, 2001; P.B. BAILEY, W.N. EVERITT AND A. ZETTL
C     VERSION 1.2
C     **********

      SUBROUTINE SLEIGN(A,B,INTAB,P0ATA,QFATA,P0ATB,QFATB,A1,A2,B1,B2,
     +                  NUMEIG,EIG,TOL,IFLAG,ISLFUN,SLFUN,NCA,NCB)

C     **********
C
C     This subroutine is designed for the calculation of a specified
C     eigenvalue, EIG, of a Sturm-Liouville problem for the equation
C
C       -(p(x)*y'(x))' + q(x)*y(x) = eig*w(x)*y(x)  on (a,b)
C
C     with user-supplied coefficient functions p, q, and w,
C     and with separated boundary conditions.  (For coupled
C     boundary conditions, see the companion subroutine SLCOUP.)
C          The problem may be either nonsingular or singular.  In
C     the nonsingular case, boundary conditions are of the form
C
C        A1*y(a) + A2*p(a)*y'(a) = 0
C        B1*y(b) + B2*p(b)*y'(b) = 0,
C
C     but are of the form
C
C        A1*[y,u](a) + A2*[y,v](a) = 0
C        B1*[y,U](b) + B2*[y,V](a) = 0,
C
C     when the endpoints are singular, of type Limit Circle.
C     In either case the boundary conditions are prescribed by
C     specifying the numbers A1, A2, B1, B2.  In the singular
C     case the user must also supply the "boundary condition
C     functions" u, v near a and/or U, V near b, whichever is
C     singular, of type Limit Circle.
C     The index of the desired eigenvalue is specified in NUMEIG
C     and its requested accuracy in TOL.  Initial data for the
C     associated eigenfunction are also computed along with values
C     at selected points, if desired, in array SLFUN.
C
C     In addition to the coefficient functions p, q, and w, the user
C     must supply subroutine UV to describe the boundary condition
C     when the problem is limit circle.  UV can be a dummy subroutine
C     if the problem is not limit circle.
C
C     The SUBROUTINE statement is
C
C     SUBROUTINE SLEIGN(A,B,INTAB,P0ATA,QFATA,P0ATB,QFATB,A1,A2,B1,B2,
C    1    NUMEIG,EIG,TOL,IFLAG,ISLFUN,SLFUN,NCA,NCB)
C
C     where
C
C       A and B are input variables defining the interval.  If the
C         interval is finite, A must be less than B.  (See INTAB below.)
C
C       INTAB is an integer input variable specifying the nature of the
C         interval.  It can have four values.
C
C         INTAB = 1 - A and B are finite.
C         INTAB = 2 - A is finite and B is infinite (+).
C         INTAB = 3 - A is infinite (-) and B is finite.
C         INTAB = 4 - A is infinite (-) and B is infinite (+).
C
C         If either A or B is infinite, it is classified singular and
C         its value is ignored.
C
C       P0ATA, QFATA, P0ATB, and QFATB are input variables set to
C         1.0 or -1.0 as the following properties of p, q, and w at
C         the interval endpoints are true or false, respectively.
C
C         P0ATA -  p(a) is zero.              (If true, A is singular.)
C         QFATA -  q(a) and w(a) are finite.  (If false, A is singular.)
C         P0ATB -  p(b) is zero.              (If true, B is singular.)
C         QFATB -  q(b) and w(b) are finite.  (If false, B is singular.)
C
C       A1 and A2 are input variables set to prescribe the boundary
C         condition at A.
C
C       B1 and B2 are input variables set to prescribe the boundary
C         condition at B.
C
C       NUMEIG is an integer variable.  On input, it should be set to
C         the index of the desired eigenvalue (increasing sequence where
C         index 0 corresponds to the lowest eigenvalue -- if the
C         eigenvalues are bounded below -- or to the smallest nonegative
C         eigenvalue otherwise).  On output, it is unchanged unless the
C         problem (apparently) lacks eigenvalue NUMEIG, in which case it
C         is reset to the index of the largest eigenvalue that seems to
C         exist.
C
C       EIG is a variable set on input to 0.0 or to an initial guess of
C         the eigenvalue.  If EIG is set to 0.0, SLEIGN2 will generate
C         the initial guess.  On output, EIG holds the calculated
C         eigenvalue if IFLAG (see below) signals success.
C
C       TOL is a variable set on input to the desired accuracy of the
C         eigenvalue.  On output, TOL is reset to the accuracy estimated
C         to have been achieved if IFLAG (see below) signals success.
C         This accuracy estimate is absolute if EIG is less than one
C         in magnitude, and relative otherwise.  In addition, prefixing
C         TOL with a negative sign, removed after interrogation, serves
C         as a flag to request trace output from the calculation.
C
C       IFLAG is an integer output variable set as follows:
C
C         IFLAG = 0 -  improper input parameters.
C         IFLAG = 1  - successful problem solution, within tolerance.
C         IFLAG = 2  - best problem result, not within tolerance.
C         IFLAG = 3  - NUMEIG exceeds actual highest eigenvalue index.
C         IFLAG = 4  - RAY and EIG fail to agree after 5 tries.
C         IFLAG = 6  - in SECANT-METHOD, ABS(DE) .LT. EPSMIN .
C         IFLAG = 7  - iterations are stuck in a loop.
C         IFLAG = 8  - number of iterations has reached the set limit.
C         IFLAG = 9  - residual truncation error dominates.
C         IFLAG = 10 - integrator tolerance cannot be reduced.
C         IFLAG = 11 - no more improvement.
C         IFLAG = 13 - AA cannot be moved in any further.
C         IFLAG = 14 - BB cannot be moved in any further.
c         iflag = 15 - Bad behavior of coefficients at an endpoint.
C         IFLAG = 16 - Could not get started.
C         IFLAG = 17 - Failed to get a bracket.
C         IFLAG = 18 - Estimator failed.
C         IFLAG = 51 - integration failure after 1st call to INTEG.
C         IFLAG = 52 - integration failure after 2nd call to INTEG.
C         IFLAG = 53 - integration failure after 3rd call to INTEG.
C         IFLAG = 54 - integration failure after 4th call to INTEG.
C
C       ISLFUN is an integer input variable set to the number of
C         selected eigenfunction values desired.  If no values are
C         desired, set ISLFUN to zero.
c         (If ISLFUN is set to -1, the result will be that SLEIGN2
c           will return to the calling program directly after sampling
c           the coefficients p,q,w .  This device is used only once,
c           in SUBROUTINE PERIO.)
C
C       SLFUN is an array of length at least 9.  On output, the first 9
C         locations contain the integration interval and initial data
C         that completely determine the eigenfunction.
C
C         SLFUN(1) - point where two pieces of eigenfunction Y match.
C         SLFUN(2) - left endpoint XAA of the (truncated) interval.
C         SLFUN(3) - value of THETA at XAA.  (Y = RHO*sin(THETA))
C         SLFUN(4) - value of F at XAA.  (RHO = exp(F))
C         SLFUN(5) - right endpoint XBB of the (truncated) interval.
C         SLFUN(6) - value of THETA at XBB.
C         SLFUN(7) - value of F at XBB.
C         SLFUN(8) - final value of integration accuracy parameter EPS.
C         SLFUN(9) - the constant Z in the polar form transformation.
C
C         F(XAA) and F(XBB) are chosen so that the eigenfunction is
C         continuous in the interval (XAA,XBB) and has weighted (by W)
C         L2-norm of 1.0 on the interval.  If ISLFUN is positive, then
C         on input the further ISLFUN locations of SLFUN specify the
C         points, in ascending order, where the eigenfunction values
C         are desired and on output contain the values themselves.
C
c       nca & ncb are integers which indicate the nature of the
c         endpoints a & b, respectively.  Namely:
c       nca = 1 : Endpoint a is REGULAR .
C             2 :     "         WEAKLY REGULAR .
C             3 :     "         LIMIT CIRCLE, NON-OSC .
C             4 :     "         LIMIT CIRCLE, OSC .
C             5 :     "         LIMIT POINT, REGULAR, AT FINITE PT.
C             6 :     "         LIMIT POINT, DEFAULT .
C             7 :     "         LIMIT POINT, AT INFINITY OR IRREG.
C             8 :     "         LIMIT POINT, BAD BEHAVIOR AT ENDPT.
C        ncb = 1 : Endpoint b is REGULAR .
C                  etc. as for nca .
C
C     **********
C  INPUT QUANTITIES: A,B,INTAB,P0ATA,QFATA,P0ATB,QFATB,A1,A2,B1,B2,
C                   NUMEIG,EIG,TOL,IFLAG,ISLFUN,SLFUN,NCA,NCB,PR
C  OUTPUT QUANTITIES: NUMEIG,EIG,TOL,IFLAG,SLFUN,
C                   MFS,MLS,PI,TWOPI,HPI,EPSMIN,Z,JAY,ZEE,
C                   AA,TMID,BB,DTHDAA,DTHDBB,MDTHZ,ADDD
C
C     .. Scalar Arguments ..
      REAL A,A1,A2,B,B1,B2,EIG,P0ATA,P0ATB,QFATA,QFATB,TOL
      INTEGER IFLAG,INTAB,ISLFUN,NCA,NCB,NUMEIG
C     ..
C     .. Array Arguments ..
      REAL SLFUN(9)
C     ..
C     .. Scalars in Common ..
      REAL A1S,A2S,AA,ASAV,B1S,B2S,BB,BSAV,DTHDAA,DTHDBB,
     +                 EIGSAV,EPSMIN,FA,FB,GQA,GQB,GWA,GWB,HPI,LPQA,
     +                 LPQB,LPWA,LPWB,P0ATAS,P0ATBS,PI,QFATAS,QFATBS,
     +                 TMID,TSAVEL,TSAVER,TWOPI,Z
      INTEGER IND,INTSAV,ISAVE,MDTHZ,MFS,MLS,MMWD,T21
      LOGICAL ADDD,PR
C     ..
C     .. Arrays in Common ..
      REAL TEE(100),TT(7,2),YS(200),YY(7,3,2),ZEE(100)
      INTEGER JAY(100),MMW(100),NT(2)
C     ..
C     .. Local Scalars ..
      REAL AA1,AAA,AAF,AAL,AAS,ALFA,ATHETA,BALLPK,BB1,BBB,
     +                 BBF,BBL,BBS,BESTAA,BESTBB,BETA,BSTEIG,BSTEPS,
     +                 BSTEST,BSTMID,CHLIM,CHNG,CONVC,DE,DEDW,DIST,
     +                 DTHDA,DTHDB,DTHDE,DTHDEA,DTHDEB,DTHETA,DTHOLD,
     +                 DTHOLY,EEE,EIGLO,EIGLT,EIGPI,EIGRT,EIGUP,EL,
     +                 ELIMUP,EMAX,EMIN,EOLD,EPS,EPSL,EPSM,ER1,ER1M,
     +                 ESTERR,FLO,FLOUP,FMAX,FUP,GUESS,OLDEST,OLDRAY,
     +                 OLRAYS,ONE,PIN,PT2,PT3,RAY,RLX,SAVAA,SAVBB,
     +                 SAVERR,SUM,T1,T2,T3,TAU,TAU0,TAUM,TMID1,TMP,U,V,
     +                 WL,ZZ
      INTEGER I,IMAX,IMID,IMIN,J,JFLAG,JJL,JJR,K,LOOP2,LOOP3,MF,ML,NEIG,
     +        NEIGST,NITER,NRAY,NTMP
      LOGICAL AOK,BOK,BRACKT,BRS,BRSS,CHEPS,CONVRG,ENDA,ENDB,EXIT,
     +        FIRSTT,GESS0,IOSC,LCOA,LCOB,LIMUP,LOGIC,NEWTF,NEWTON,
     +        OLNEWT,ONEDIG,THEGT0,THELT0,TRUNKA,TRUNKB
C     ..
C     .. Local Arrays ..
      REAL DELT(100),DS(100),EIGEST(200,4),PS(100),PSS(100),
     +                 QS(100),WS(100),XS(100)
C     ..
C     .. External Functions ..
      REAL EPSLON
      EXTERNAL EPSLON
C     ..
C     .. External Subroutines ..
      EXTERNAL AABB,EFDATA,EIGFCN,ESTIM,OBTAIN,PRELIM,RESET,SAMPLE,SAVE
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS,ATAN,INT,MAX,MIN,SIGN
C     ..
C     .. Common blocks ..
c     COMMON /ALBE/LPWA,LPQA,LPWB,LPQB
      COMMON /ALBE/LPWA,LPQA,FA,GWA,GQA,LPWB,LPQB,FB,GWB,GQB
      COMMON /BCDATA/A1S,A2S,P0ATAS,QFATAS,B1S,B2S,P0ATBS,QFATBS
      COMMON /DATADT/ASAV,BSAV,INTSAV
      COMMON /DATAF/EIGSAV,IND
      COMMON /LP/MFS,MLS
      COMMON /PASS/YS,MMW,MMWD
      COMMON /PIE/PI,TWOPI,HPI
      COMMON /PRIN/PR,T21
      COMMON /RNDOFF/EPSMIN
      COMMON /TDATA/AA,TMID,BB,DTHDAA,DTHDBB,MDTHZ,ADDD
      COMMON /TEEZ/TEE
      COMMON /TEMP/TT,YY,NT
      COMMON /TSAVE/TSAVEL,TSAVER,ISAVE
      COMMON /Z1/Z
      COMMON /ZEEZ/JAY,ZEE
C     ..
C     To produce printout, set PR = .TRUE.
      PR = .true.
C
      IFLAG = 1
      ONE = 1.0D0
      EPSMIN = EPSLON(ONE)
        if(epsmin .le. 1.0d-12) epsmin = 1.0d-12
      PI = 4.0D0*ATAN(ONE)
      TWOPI = 2.0D0*PI
      HPI = 0.5D0*PI
      Z = 1.0D0
      NEIG = NUMEIG - 1
C
      LOGIC = 1 .LE. INTAB .AND. INTAB .LE. 4 .AND.
     +        P0ATA*QFATA*P0ATB*QFATB .NE. 0.0D0
      IF (INTAB.EQ.1) LOGIC = LOGIC .AND. A .LT. B
      IF (.NOT.LOGIC) THEN
          IFLAG = 0
          GO TO 150
      END IF
C
C     (JFLAG = 0 INDICATES FAILURE OF THE ESTIMATOR)
      JFLAG = 1
C
      CALL SAVE(A,B,INTAB,P0ATA,QFATA,P0ATB,QFATB,A1,A2,B1,B2)
C
      DO 5 I = 1,100
          JAY(I) = 0
          ZEE(I) = 1.0D0
    5 CONTINUE
C
      CALL SAMPLE(NCA,NCB,MF,ML,XS,PS,QS,WS,DS,DELT,EMIN,EMAX,IMIN,IMAX,
     +            AAA,BBB)

C       (MAKE MF,ML AVALABLE THROUGH COMMON/LP/ )
      MFS = MF
      MLS = ML
C
      EIGPI = NUMEIG*PI
      PIN = EIGPI + PI
      TAU = ABS(TOL)
      TAU0 = 0.001D0
      TAUM = MAX(TAU,EPSMIN)
      LIMUP = .FALSE.
      ELIMUP = EMAX
      GUESS = EIG
      GESS0 = ABS(EIG) .LE. (1D-7)
      LCOA = NCA .EQ. 4
      LCOB = NCB .EQ. 4
      IOSC = LCOA .OR. LCOB

      IF (GESS0) THEN
          CALL ESTIM(IOSC,PIN,MF,ML,PS,QS,WS,PSS,DS,DELT,TAU0,JJL,JJR,
     +               SUM,LIMUP,ELIMUP,EMAX,IMAX,IMIN,EL,WL,DEDW,BALLPK,
     +               EEE,JFLAG)
      END IF

      IF (JFLAG.EQ.0) THEN
          IF (PR) WRITE (T21,FMT=*) ' ESTIMATOR FAILED '
          IFLAG = 18
          RETURN
      END IF

      CALL PRELIM(MF,ML,EEE,BALLPK,XS,QS,WS,DS,DELT,PS,PSS,TAU0,JJL,JJR,
     +            LIMUP,ELIMUP,EL,WL,DEDW,GUESS,EIG,AAA,AA,BBB,BB,NCA,
     +            NCB,EIGPI,ALFA,BETA,DTHDEA,DTHDEB,IMID,TMID)

      IF (JAY(3).EQ.0 .AND. ZEE(2).NE.1.0D0) THEN
          IF (PR) WRITE (T21,FMT=*) ' COULD NOT GET STARTED. '
          IFLAG = 16
          RETURN
      END IF
C
      IF (NCA.EQ.8 .OR. NCB.EQ.8) THEN
C       (WE CAN'T HANDLE THIS KIND OF LIMIT POINT)
          IFLAG = 15
          RETURN
      END IF
C
      IF (ISLFUN.EQ.-1) THEN
          SLFUN(1) = TMID
          SLFUN(2) = AA
          SLFUN(3) = ALFA
          SLFUN(5) = BB
          SLFUN(6) = BETA + EIGPI
          SLFUN(9) = Z
          ADDD = .FALSE.
          MDTHZ = 0
          IFLAG = 1
          RETURN
      END IF
C
C     SET LOGICAL VARIABLES:
      AOK = INTAB .LT. 3.D0 .AND. P0ATA .LT. 0.0D0 .AND.
     +      QFATA .GT. 0.0D0
      BOK = (INTAB.EQ.1 .OR. INTAB.EQ.3) .AND. P0ATB .LT. 0.0D0 .AND.
     +      QFATB .GT. 0.0D0
      LCOA = NCA .EQ. 4
      LCOB = NCB .EQ. 4
      TRUNKA = LCOA .OR. (NCA.EQ.6)
      TRUNKB = LCOB .OR. (NCB.EQ.6)
C
C     END PRELIMINARY WORK, BEGIN MAIN TASK OF COMPUTING EIG.
C
C     LOGICAL VARIABLES HAVE THE FOLLOWING MEANINGS IF TRUE.
C        AOK    - ENDPOINT A IS NOT SINGULAR.
C        BOK    - ENDPOINT B IS NOT SINGULAR.
C        BRACKT - EIG HAS BEEN BRACKETED.
C        CONVRG - CONVERGENCE TEST FOR EIG HAS BEEN SUCCESSFULLY PASSED.
C        NEWTON - NEWTON ITERATION MAY BE EMPLOYED.
C        THELT0 - LOWER BOUND FOR EIG HAS BEEN FOUND.
C        THEGT0 - UPPER BOUND FOR EIG HAS BEEN FOUND.
C        LIMIT  - UPPER BOUND EXISTS WITH BOUNDARY CONDITIONS SATISFIED.
C        ONEDIG - MOST SIGNIFICANT DIGIT CAN BE EXPECTED TO BE CORRECT.
C
C     INITIALIZE SOME OF THE CONTROL VARIABLES
      EXIT = .FALSE.
      FIRSTT = .TRUE.
      LOOP2 = 0
      LOOP3 = 0
      PT2 = 0.0D0
      PT3 = 0.0D0
      ENDA = .FALSE.
      ENDB = .FALSE.
      EPSM = EPSMIN
      CHEPS = .FALSE.
      NEWTF = .FALSE.
      BRS = .FALSE.
      BRSS = .FALSE.
      SAVERR = 1.0D+9
      ATHETA = 1.0D+9
      BSTEST = 1.0D+9
      OLDEST = 1.0D+9
      OLDRAY = 1.0D+9
      NRAY = 1
      AAL = AA
      BBL = BB
      AAS = AA
      BBS = BB
      AAF = AAA
      BBF = BBB
      NEIGST = 0
      EIG = EEE
      BSTEIG = EIG
      EPS = 0.0001D0
      EPSL = EPS
      BESTAA = AA
      BESTBB = BB
      BSTMID = TMID
      BSTEPS = EPS
C
  110 CONTINUE
C       (INITIAL-IZE)
      BRACKT = .FALSE.
      CONVRG = .FALSE.
      THELT0 = .FALSE.
      THEGT0 = .FALSE.
      NEWTON = .FALSE.
      EIGLO = EMIN - 1.0D0
      FLO = -5.0D0
      FUP = 5.0D0
      EIGLT = 0.0D0
      EIGRT = 0.0D0
      EIGUP = EMAX + 1.0D0
      IF (LIMUP) EIGUP = MIN(EMAX,ELIMUP)
      DTHOLD = 1.0D0
C
      IF (PR) WRITE (T21,FMT=*)
      IF (PR) WRITE (T21,FMT=*)
     +    '---------------------------------------------'
      IF (PR) WRITE (T21,FMT=*) ' INITIAL GUESS FOR EIG = ',EIG
      IF (PR) WRITE (T21,FMT=*) ' AA,BB = ',AA,BB

C           (UNTIL(CONVRG .OR. EXIT)
      DO 120 NITER = 1,40
          IF (PR) WRITE (*,FMT=*)
          IF (PR) WRITE (*,FMT=*) ' ******************** '
          IF (PR) WRITE (*,FMT=*) ' EIGENVALUE ',NUMEIG
C
          CALL RESET(AA,BB,ALFA,BETA,EIG,GUESS,MF,ML,QS,WS,IMID,TMID,
     +               LIMUP,ELIMUP,DTHDEA,DTHDEB,THELT0,EIGLO,EIGUP,
     +               BRACKT,NCA,NCB,IFLAG)
          IF (IFLAG.EQ.11) THEN
c             EXIT = .TRUE.
c             GO TO 130
              GO TO 333
          END IF
C
          AA1 = AA
          BB1 = BB
          TMID1 = TMID
          IFLAG = 1
          CALL OBTAIN(NCA,ALFA,DTHDEA,NCB,BETA,DTHDEB,AA1,BB1,EIG,EIGPI,
     +                TMID1,EPS,DTHETA,DTHDE,ONEDIG,DTHDA,DTHDB,ER1,
     +                ER1M,IFLAG)
          IF (IFLAG.EQ.2) THEN
              AAS = AA
              BBS = BB
          END IF
          ATHETA = ABS(DTHETA)
          IF (51.LE.IFLAG .AND. IFLAG.LE.54) THEN
              EXIT = .TRUE.
              GO TO 130
          ELSE
              FIRSTT = .FALSE.
          END IF
C-----------------------------------------------------------
          CHEPS = .FALSE.
          CONVRG = .FALSE.
          OLNEWT = NEWTON
          NEWTON = ABS(DTHETA) .LT. 0.06D0 .AND. BRACKT
          IF (NEWTON) ONEDIG = ONEDIG .OR.
     +                         ABS(DTHETA+ER1) .LT. 0.5D0*DTHOLD
          IF (.NOT.ONEDIG .AND. BRS) THEN
              EXIT = .TRUE.
              IF (PR) WRITE (T21,FMT=*) ' NOT ONEDIG '
              GO TO 130
          END IF
C
          IF (PR) WRITE (*,FMT=*) ' SET BRACKET '
          IF (DTHETA.GT.0.0D0) THEN
              IF (.NOT.THEGT0 .OR. EIG.LE.EIGUP) THEN
                  THEGT0 = .TRUE.
                  EIGUP = EIG
                  FUP = DTHETA
                  EIGRT = EIG - DTHETA/DTHDE
              END IF
          ELSE
              IF (.NOT.THELT0 .OR. EIG.GE.EIGLO) THEN
                  THELT0 = .TRUE.
                  EIGLO = EIG
                  FLO = DTHETA
                  EIGLT = EIG - DTHETA/DTHDE
              END IF
          END IF
C
C     EIG IS BRACKETED WHEN BOTH THEGT0=.TRUE. AND THELT0=.TRUE.
C
          BRACKT = THELT0 .AND. THEGT0
          IF (BRACKT) LOOP2 = 0
C              (TEST-FOR-CONVERGENCE)
C
C     MEASURE CONVERGENCE AFTER ADDING SEPARATE CONTRIBUTIONS TO ERROR.
C
          FLOUP = MIN(ABS(FLO),ABS(FUP))
          T1 = (ABS(DTHETA)+ER1M)/ABS(DTHDE)
          IF (TRUNKA) THEN
              T2 = (1.0D0+AA)*ABS(DTHDA)/ABS(DTHDE)
              PT2 = (AAF-AA)*DTHDA/DTHDE
          ELSE
              T2 = 0.0D0
          END IF
          IF (TRUNKB) THEN
              T3 = (1.0D0-BB)*ABS(DTHDB)/ABS(DTHDE)
              PT3 = (BBF-BB)*DTHDB/DTHDE
          ELSE
              T3 = 0.0D0
          END IF
          IF (PR) WRITE (T21,FMT=*) ' FLO,FUP,FLOUP = ',FLO,FUP,FLOUP
          IF (PR) WRITE (T21,FMT=*) ' DTHDE,DTHDA,DTHDB = ',DTHDE,DTHDA,
     +        DTHDB
          ESTERR = T1 + T2 + T3
          CONVC = T1 + T2 + T3

          IF (BRACKT) THEN
              TMP = EIGUP - EIGLO
c             ESTERR = MIN(ESTERR,TMP)
              CONVC = MIN(ESTERR,TMP)
          END IF

 333      CONTINUE
          ESTERR = ESTERR/MAX(ONE,ABS(EIG))
          CONVC = CONVC/MAX(ONE,ABS(EIG))
          NEIGST = NEIGST + 1
          EIGEST(NEIGST,1) = EIG
          EIGEST(NEIGST,2) = DTHETA
          EIGEST(NEIGST,3) = ESTERR
          EIGEST(NEIGST,4) = 0.0D0
c         CONVRG = ESTERR .LE. TAUM .AND. NEWTON
          CONVRG = CONVC .LE. TAUM .AND. NEWTON
          IF (PR) WRITE (T21,FMT=*) ' T1,T2,T3 = ',T1,T2,T3
          IF (PR) WRITE (T21,FMT=*) ' PT2,PT3 = ',PT2,PT3
          IF (PR) WRITE (T21,FMT=*) ' TMID,EPS = ',TMID,EPS
          IF (PR) WRITE (T21,FMT=*) ' ONEDIG,BRACKT,NEWTON,CONVRG = ',
     +        ONEDIG,BRACKT,NEWTON,CONVRG
          IF (PR) WRITE (T21,FMT=*) ' EIG,DTHETA,ESTERR = ',EIG,DTHETA,
     +        ESTERR

          IF (BRACKT .AND. (ESTERR.LT.BSTEST.OR..NOT.BRS) .AND.
     +        T1.LT.0.1D0) THEN
              BESTAA = AA
              BESTBB = BB
              BSTMID = TMID
              BSTEPS = EPS
              BSTEIG = EIG
              BSTEST = ESTERR
              BRS = BRACKT
              IF (BRS) BRSS = BRS
              IF (PR) WRITE (T21,FMT=*) ' BSTEIG,BSTEST = ',BSTEIG,
     +            BSTEST
              IF (PR) WRITE (T21,FMT=*) ' BRS = ',BRS
          END IF

          IF (THEGT0 .AND. PR) WRITE (T21,FMT=*) '         EIGUP = ',
     +        EIGUP
          IF (THELT0 .AND. PR) WRITE (T21,FMT=*) '         EIGLO = ',
     +        EIGLO
          IF (PR) WRITE (T21,FMT=*)
     +        '-------------------------------------------'
          IF (CONVRG) THEN
              IF (PR) WRITE (*,FMT=*) ' NUMBER OF ITERATIONS WAS ',NITER
              IF (PR) WRITE (T21,FMT=*) ' NUMBER OF ITERATIONS WAS ',
     +            NITER
              IF (PR) WRITE (*,FMT=*)
     +            '-----------------------------------------'
              GO TO 130
          ELSE
              IF (NEWTON) THEN
                  IF (OLNEWT .AND. ATHETA.GT.0.8D0*ABS(DTHOLD)) THEN
                      IF (PR) WRITE (T21,FMT=*) ' ATHETA,DTHOLD = ',
     +                    ATHETA,DTHOLD
                      IF (PR) WRITE (T21,FMT=*
     +                    ) ' NEWTON DID NOT IMPROVE EIG '
                      NEWTF = .TRUE.
                      LOOP3 = LOOP3 + 1
                  ELSE IF (TRUNKA .OR. TRUNKB) THEN
                      ENDA = ATHETA .LT. 1.0D0 .AND. TRUNKA .AND.
     +                       ABS(PT2) .GT. MAX(TAUM,T1)
                      ENDB = ATHETA .LT. 1.0D0 .AND. TRUNKB .AND.
     +                       ABS(PT3) .GT. MAX(TAUM,T1)
                      IF (ENDA .OR. ENDB) THEN
                          NEWTON = .FALSE.
                      ELSE IF (((T2+T3).GT.T1) .AND.
     +                         (ATHETA.LT.1.0D0) .AND.
     +                         (AA.LE.AAF.AND.BB.GE.BBF)) THEN
                          IF (PR) WRITE (*,FMT=*
     +                        ) ' RESIDUAL TRUNCATION ERROR DOMINATES '
                          EXIT = .TRUE.
                          IFLAG = 9
                          IF (PR) WRITE (T21,FMT=*) ' IFLAG = 9 '
                          GO TO 130
                      END IF
                  END IF
                  IF (NEWTF .OR. ENDA .OR. ENDB) THEN
                      IF (PR) WRITE (T21,FMT=*) ' NEWTF,ENDA,ENDB = ',
     +                    NEWTF,ENDA,ENDB
                      EXIT = .TRUE.
                      GO TO 130
                  END IF
C                  (NEWTON'S-METHOD)
                  IF (PR) WRITE (*,FMT=*) ' NEWTON''S METHOD '
                  RLX = 1.2D0
                  IF (BRACKT) RLX = 1.0D0
                  EIG = EIG - RLX*DTHETA/DTHDE
                  IF (EIG.LE.EIGLO .OR. EIG.GE.EIGUP) EIG = 0.5D0*
     +                (EIGLO+EIGUP)
                  IF (PR) WRITE (T21,FMT=*) ' NEWTON: EIG = ',EIG
              ELSE IF (BRACKT) THEN
                  IF (PR) WRITE (*,FMT=*) ' BRACKET '
C                  (SECANT-METHOD)
                  IF (PR) WRITE (*,FMT=*) ' DO SECANT METHOD '
                  FMAX = MAX(-FLO,FUP)
                  EOLD = EIG
                  EIG = 0.5D0* (EIGLO+EIGUP)
                  IF (FMAX.LE.1.5D0) THEN
                      U = -FLO/ (FUP-FLO)
                      DIST = EIGUP - EIGLO
                      EIG = EIGLO + U*DIST
                      V = MIN(EIGLT,EIGRT)
                      IF (EIG.LE.V) EIG = 0.5D0* (EIG+V)
                      V = MAX(EIGLT,EIGRT)
                      IF (EIG.GE.V) EIG = 0.5D0* (EIG+V)
                      DE = EIG - EOLD
                      IF (ABS(DE).LT.EPSMIN) THEN
                          TOL = ABS(DE)/MAX(ONE,ABS(EIG))
                          IFLAG = 6
                          EXIT = .TRUE.
                          GO TO 130
                      END IF
                  END IF
                  IF (PR) WRITE (T21,FMT=*) ' SECANT: EIG = ',EIG
              ELSE
C                  (TRY-FOR-BRACKET)
                  LOOP2 = LOOP2 + 1
                  IF (LOOP2.GT.9 .AND. .NOT.LIMUP) THEN
                      IFLAG = 12
                      IF (PR) WRITE (T21,FMT=*) ' IFLAG = 12 '
                      EXIT = .TRUE.
                      GO TO 130
                  END IF
                  IF (EIG.EQ.EEE) THEN
                      IF (GUESS.NE.0.0D0) DEDW = 1.0D0/DTHDE
                      CHNG = -0.6D0* (DEDW+1.0D0/DTHDE)*DTHETA
                      IF (EIG.NE.0.0D0 .AND. ABS(CHNG).GT.
     +                    0.1D0*ABS(EIG)) CHNG = -0.1D0*SIGN(EIG,DTHETA)
                  ELSE
                      CHNG = -1.2D0*DTHETA/DTHDE
                      IF (CHNG.EQ.0.D0) CHNG = 0.1D0*MAX(ONE,ABS(EIG))
                      IF (PR) WRITE (T21,FMT=*
     +                    ) ' IN BRACKET, 1,CHNG = ',CHNG
C
C     LIMIT CHANGE IN EIG TO A FACTOR OF 2.
C
                      IF (ABS(CHNG).GT. (1.0D0+2.0D0*ABS(EIG))) THEN
                          TMP = 1.0D0+2.0D0*ABS(EIG)
                          CHNG = SIGN(TMP,CHNG)
                          IF (PR) WRITE (T21,FMT=*
     +                        ) ' IN BRACKET, 2,CHNG = ',CHNG
                      ELSE IF (ABS(EIG).GE.1.0D0 .AND.
     +                         ABS(CHNG).LT.0.1D0*ABS(EIG)) THEN
                          CHNG = 0.1D0*SIGN(EIG,CHNG)
                          IF (PR) WRITE (T21,FMT=*
     +                        ) ' IN BRACKET, 3,CHNG = ',CHNG
                      END IF
                      IF (DTHOLD.LT.0.0D0 .AND. LIMUP .AND.
     +                    CHNG.GT. (ELIMUP-EIG)) THEN
                          CHNG = 0.95D0* (ELIMUP-EIG)
                          IF (PR) WRITE (T21,FMT=*
     +                        ) ' ELIMUP,EIG,CHNG = ',ELIMUP,EIG,CHNG
                          IF (CHNG.LT.EPSMIN) THEN
                              IF (PR) WRITE (*,FMT=*) ' ELIMUP,EIG = ',
     +                            ELIMUP,EIG
                              IF (PR) WRITE (T21,FMT=*
     +                            ) ' IN BRACKET, CHNG.LT.EPSMIN '
                              ENDA = TRUNKA .AND. AA .GT. AAF .AND.
     +                               (0.5D0*DTHETA+ (AAF-AA)*DTHDA) .GT.
     +                               0.0D0
                              ENDB = TRUNKB .AND. BB .LT. BBF .AND.
     +                               (0.5D0*DTHETA- (BBF-BB)*DTHDB) .GT.
     +                               0.0D0
                              IF (.NOT. (ENDA.OR.ENDB)) THEN
                                  NUMEIG = NEIG - INT(-DTHETA/PI)
                                  IF (PR) WRITE (*,FMT=*
     +                                ) ' NEW NUMEIG = ',NUMEIG
                                  IF (PR) WRITE (T21,
     +                                FMT=*) ' NEW NUMEIG = ',NUMEIG
                              END IF
                              IFLAG = 3
                              GO TO 150
                          END IF
                      END IF
                  END IF
                  EOLD = EIG
                  CHLIM = 2.0D0*ESTERR*MAX(ONE,ABS(EIG))
                  IF (ATHETA.LT.0.06D0 .AND. ABS(CHNG).GT.CHLIM .AND.
     +                CHLIM.NE.0.0D0) CHNG = SIGN(CHLIM,CHNG)
                  IF ((THELT0.AND.CHNG.LT.0.0D0) .OR.
     +                (THEGT0.AND.CHNG.GT.0.0D0)) CHNG = -CHNG
                  EIG = EIG + CHNG
                  IF (PR) WRITE (T21,FMT=*) ' BRACKET: EIG = ',EIG
              END IF
          END IF
          IF (IFLAG.EQ.3) GO TO 130
          IF (NITER.GE.3 .AND. DTHOLY.EQ.DTHETA) THEN
              IFLAG = 7
              IF (PR) WRITE (T21,FMT=*) ' IFLAG = 7 '
              BSTEIG = EIG
              BSTEST = ESTERR
              EXIT = .TRUE.
              GO TO 130
          END IF
          DTHOLY = DTHOLD
          DTHOLD = DTHETA
          IF (PR) WRITE (*,FMT=*) ' NUMBER OF ITERATIONS WAS ',NITER
          IF (PR) WRITE (*,FMT=*)
     +        '-----------------------------------------------'
  120 CONTINUE
      IFLAG = 8
      IF (PR) WRITE (T21,FMT=*) ' IFLAG = 8 '
      EXIT = .TRUE.
  130 CONTINUE
      IF (AA.EQ.AAL .AND. BB.EQ.BBL .AND. EPS.LT.EPSL .AND.
     +    ESTERR.GE.0.5D0*SAVERR .AND. .NOT.EXIT) GO TO 140
      EPSL = EPS
      AAL = AA
      BBL = BB
      TOL = BSTEST
      EIG = BSTEIG
      IF (EXIT) THEN
          IF (PR) WRITE (T21,FMT=*) ' EXIT '
          IF (FIRSTT) THEN
              IF (IFLAG.EQ.51 .OR. IFLAG.EQ.53) THEN
                  IF (AA.LT.-0.71D0) THEN
                      IF (PR) WRITE (T21,FMT=*
     +                    ) ' FIRST COMPLETE INTEGRATION FAILED. '
                      IF (AA.EQ.-1.0D0) GO TO 150
                      AAF = AA
                      CALL AABB(AA,-ONE)
                      IF (PR) WRITE (T21,FMT=*) ' AA MOVED FROM ',AAF,
     +                    ' IN TO ',AA
                      EXIT = .FALSE.
                      GO TO 110
                  ELSE
                      IF (PR) WRITE (T21,FMT=*) ' AA.GE.-0.71 '
                      IFLAG = 13
                      GO TO 150
                  END IF
              ELSE IF (IFLAG.EQ.52 .OR. IFLAG.EQ.54) THEN
                  IF (BB.GT.0.71D0) THEN
                      IF (PR) WRITE (T21,FMT=*
     +                    ) ' FIRST COMPLETE INTEGRATION FAILED. '
                      IF (BB.EQ.1.0D0) GO TO 150
                      BBF = BB
                      CALL AABB(BB,-ONE)
                      IF (PR) WRITE (T21,FMT=*) ' BB MOVED FROM ',BBF,
     +                    ' IN TO ',BB
                      EXIT = .FALSE.
                      GO TO 110
                  ELSE
                      IF (PR) WRITE (T21,FMT=*) ' BB.LE.0.71 '
                      IFLAG = 14
                      GO TO 150
                  END IF
              END IF
          ELSE IF (IFLAG.EQ.51 .OR. IFLAG.EQ.53) THEN
              IF (PR) WRITE (*,FMT=*) ' A COMPLETE INTEGRATION FAILED. '
              IF (PR) WRITE (T21,FMT=*)
     +            ' A COMPLETE INTEGRATION FAILED. '
              IF (CHEPS) THEN
                  EPS = 5.0D0*EPS
                  EPSM = EPS
                  IF (PR) WRITE (T21,FMT=*) ' EPS INCREASED TO ',EPS
              ELSE
                  AAF = AA
                  CALL AABB(AA,-ONE)
                  IF (PR) WRITE (T21,FMT=*) ' AA MOVED FROM ',AAF,
     +                ' IN TO ',AA
              END IF
              EXIT = .FALSE.
              GO TO 110
          ELSE IF (IFLAG.EQ.52 .OR. IFLAG.EQ.54) THEN
              IF (PR) WRITE (*,FMT=*) ' A COMPLETE INTEGRATION FAILED. '
              IF (PR) WRITE (T21,FMT=*)
     +            ' A COMPLETE INTEGRATION FAILED. '
              IF (CHEPS) THEN
                  EPS = 5.0D0*EPS
                  EPSM = EPS
                  IF (PR) WRITE (T21,FMT=*) ' EPS INCREASED TO ',EPS
              ELSE
                  BBF = BB
                  CALL AABB(BB,-ONE)
                  IF (PR) WRITE (T21,FMT=*) ' BB MOVED FROM ',BBF,
     +                ' IN TO ',BB
              END IF
              EXIT = .FALSE.
              GO TO 110
          ELSE IF (IFLAG.EQ.6) THEN
              IF (PR) WRITE (*,FMT=*) ' IN SECANT, CHNG.LT.EPSMIN '
              IF (PR) WRITE (T21,FMT=*) ' IN SECANT, CHNG.LT.EPSMIN '
              GO TO 140
          ELSE IF (IFLAG.EQ.7) THEN
              IF (PR) WRITE (*,FMT=*) ' DTHETA IS REPEATING '
              IF (PR) WRITE (T21,FMT=*) ' DTHETA IS REPEATING '
              GO TO 140
          ELSE IF (IFLAG.EQ.8) THEN
              IF (PR) WRITE (*,FMT=*)
     +            ' NUMBER OF ITERATIONS REACHED SET LIMIT '
              IF (PR) WRITE (T21,FMT=*)
     +            ' NUMBER OF ITERATIONS REACHED SET LIMIT '
              GO TO 140
          ELSE IF (IFLAG.EQ.9) THEN
              IF (PR) WRITE (T21,FMT=*)
     +            ' RESIDUAL TRUNCATION ERROR DOMINATES '
              GO TO 140
          ELSE IF (IFLAG.EQ.11) THEN
              IF (PR) WRITE (*,FMT=*)
     +            ' IN TRY FOR BRACKET, CHNG.LT.EPSMIN '
              IF (PR) WRITE (T21,FMT=*)
     +            ' IN TRY FOR BRACKET, CHNG.LT.EPSMIN '
              GO TO 140
          ELSE IF (IFLAG.EQ.12) THEN
              IF (PR) WRITE (*,FMT=*) ' FAILED TO GET A BRACKET. '
              IF (PR) WRITE (T21,FMT=*) ' FAILED TO GET A BRACKET. '
              GO TO 140
          ELSE IF (NEWTF .OR. .NOT.ONEDIG) THEN
              IF (LOOP3.GE.3) THEN
                  IF (PR) WRITE (T21,FMT=*)
     +                ' NEWTON IS NOT GETTING ANYWHERE '
                  NEWTF = .FALSE.
                  GO TO 140
              END IF
              IF (PR) WRITE (T21,FMT=*) ' BSTEST,OLDEST = ',BSTEST,
     +            OLDEST
              IF (EPS.GT.EPSM .AND. BSTEST.LT.OLDEST) THEN
                  CHEPS = .TRUE.
                  SAVERR = ESTERR
                  EPS = 0.2D0*EPS
                  IF (BSTEST.LT.0.001D0) EPS = EPSM
                  IF (PR) WRITE (T21,FMT=*) ' EPS REDUCED TO ',EPS
                  EXIT = .FALSE.
                  NEWTON = .FALSE.
                  OLDEST = BSTEST
                  GO TO 110
              ELSE
                  IF (EPS.LE.EPSM) THEN
                      IF (PR) WRITE (T21,FMT=*
     +                    ) ' EPS CANNOT BE REDUCED FURTHER. '
                      IFLAG = 2
                      GO TO 140
                  ELSE
                      IF (PR) WRITE (*,FMT=*) ' NO MORE IMPROVEMENT '
                      IF (PR) WRITE (T21,FMT=*) ' NO MORE IMPROVEMENT '
                      GO TO 140
                  END IF
              END IF
          ELSE IF (ENDA) THEN
              CALL AABB(AA,ONE)
              AA = MAX(AA,AAA)
              IF (AA.LE.AAF) THEN
                  IF (PR) WRITE (T21,FMT=*) ' NO MORE IMPROVEMENT '
                  CALL AABB(AA,-ONE)
                  GO TO 140
              END IF
              IF (PR) WRITE (*,FMT=*) ' AA MOVED OUT TO ',AA
              IF (PR) WRITE (T21,FMT=*) ' AA MOVED OUT TO ',AA
              EXIT = .FALSE.
              GO TO 110
          ELSE IF (ENDB) THEN
              CALL AABB(BB,ONE)
              BB = MIN(BB,BBB)
              IF (BB.GE.BBF) THEN
                  IF (PR) WRITE (T21,FMT=*) ' NO MORE IMPROVEMENT '
                  CALL AABB(BB,-ONE)
                  GO TO 140
              END IF
              IF (PR) WRITE (*,FMT=*) ' BB MOVED OUT TO ',BB
              IF (PR) WRITE (T21,FMT=*) ' BB MOVED OUT TO ',BB
              EXIT = .FALSE.
              GO TO 110
          END IF
      END IF
  140 CONTINUE
C
C     IF CONVRG IS FALSE, CHECK THAT ANY TRUNCATION ERROR MIGHT POSSIBLY
C     BE REDUCED OR THAT THE INTEGRATIONS MIGHT BE DONE MORE ACCURATELY.
C
      IF (.NOT.CONVRG .AND. IFLAG.LT.50 .AND. IFLAG.NE.11) THEN
          SAVAA = AA
          SAVBB = BB
          IF (EPS.GT.EPSM .AND. ESTERR.LT.0.5D0*SAVERR) THEN
              IF (PR) WRITE (T21,FMT=*) ' SAVERR,ESTERR = ',SAVERR,
     +            ESTERR
              SAVERR = ESTERR
              EPS = 0.2D0*EPS
              IF (ESTERR.LT.0.001D0) EPS = EPSM
              IF (PR) WRITE (T21,FMT=*) ' 2,EPS REDUCED TO ',EPS
              EXIT = .FALSE.
              NEWTON = .FALSE.
              IF (DTHOLY.EQ.DTHETA) THEN
                  BSTEST = ESTERR
                  BSTEIG = EIG
              END IF
              OLDEST = BSTEST
              GO TO 110
          ELSE IF (ABS(PT2).GT.TAUM .OR. ABS(PT3).GT.TAUM) THEN
              IF ((AAS-AAF).GT.2.0D0*EPSMIN .AND. ABS(PT2).GT.TAUM) THEN
                  CALL AABB(AA,ONE)
                  AA = MAX(AA,AAA)
                  IF (AA.GT.AAF .AND. AA.LT.SAVAA) THEN
                      IF (PR) WRITE (*,FMT=*) ' AA MOVED OUT TO ',AA
                      IF (PR) WRITE (T21,FMT=*) ' 3,AA MOVED OUT TO ',AA
                      EXIT = .FALSE.
                  END IF
              END IF
              IF ((BBF-BBS).GT.2.0D0*EPSMIN .AND. ABS(PT3).GT.TAUM) THEN
                  CALL AABB(BB,ONE)
                  BB = MIN(BB,BBB)
                  IF (BB.GT.SAVBB .AND. BB.LT.BBF) THEN
                      IF (PR) WRITE (*,FMT=*) ' BB MOVED OUT TO ',BB
                      IF (PR) WRITE (T21,FMT=*) ' 3,BB MOVED OUT TO ',BB
                      EXIT = .FALSE.
                  END IF
              END IF
              IF (.NOT.EXIT .AND. (AA.NE.SAVAA.OR.BB.NE.SAVBB))
     +            GO TO 110
          END IF
      END IF
      IF (PR) WRITE (T21,FMT=*) 'NUMEIG = ',NUMEIG,' EIG = ',EIG,
     +    ' TOL = ',TOL
      IF (BRSS .AND. BSTEST.LT.0.05D0) THEN

C          DO (COMPUTE-EIGENFUNCTION-DATA)
          EPS = BSTEPS
          CALL EFDATA(ALFA,DTHDEA,A1,A2,BETA,DTHDEB,B1,B2,EIG,EIGPI,EPS,
     +                AOK,NCA,BOK,NCB,RAY,SLFUN)
C
C     IF NEXT CONDITION IS .TRUE., THEN SOMETHING IS APPARENTLY WRONG
C     WITH THE ACCURACY OF EIG. BISECT AND GO THROUGH THE LOOP AGAIN.
C
          IF (PR) WRITE (T21,FMT=*) ' EIG,RAY = ',EIG,RAY
          IF (ABS(RAY-EIG).GT.2.0D0*TAUM*MAX(ONE,ABS(EIG))) THEN
              NRAY = NRAY + 1
              IF (PR) WRITE (*,FMT=*) ' NRAY = ',NRAY
              IF (PR) WRITE (T21,FMT=*) ' NRAY,RAY,OLDRAY = ',NRAY,RAY,
     +            OLDRAY
              IF (ESTERR.GE.0.5D0*SAVERR) THEN
                  IFLAG = 2
                  GO TO 150
              END IF
              TMP = EIG
              EIG = 0.5D0* (EIG+RAY)
              OLRAYS = OLDRAY
              OLDRAY = RAY
              IF (OLDRAY.NE.OLRAYS .AND. NRAY.LT.2) THEN
                  GO TO 110
              ELSE
                  EIG = TMP
              END IF
          END IF
C     DO (GENERATE-EIGENFUNCTION-VALUES)
          CALL EIGFCN(EIGPI,A1,A2,B1,B2,AOK,BOK,NCA,NCB,SLFUN,ISLFUN)
          IFLAG = 1
      END IF
C
C     IF THE ESTIMATED ACCURACY IMPLIES THAT THE COMPUTED VALUE
C     OF THE EIGENVALUE IS UNCERTAIN, SIGNAL BY IFLAG = 2.
C
      IF ((ABS(EIG).LE.ONE.AND.TOL.GE.ABS(EIG)) .OR.
     +    (ABS(EIG).GT.ONE.AND.TOL.GE.ONE)) IFLAG = 2
  150 CONTINUE
      IF ((IFLAG.EQ.2) .OR. (6.LE.IFLAG.AND.IFLAG.LE.11)) THEN
          NTMP = 1
          IF (NEIGST.GE.10) NTMP = NEIGST - 9
          K = NTMP
          TMP = EIGEST(NTMP,4)
          DO 160 I = NTMP,NEIGST
              IF (EIGEST(I,4).LT.TMP) THEN
                  K = I
                  TMP = EIGEST(I,4)
              END IF
  160     CONTINUE
          J = K
          TMP = EIGEST(K,3)
          DO 155 I = K,NEIGST
              IF (EIGEST(I,3).LT.TMP) THEN
                  J = I
                  TMP = EIGEST(I,3)
              END IF
  155     CONTINUE
          EIG = EIGEST(J,1)
          TOL = EIGEST(J,3)
          IFLAG = 2
      END IF
C
C     TO BE SAFE, RESET AA,BB,EPS:-
      AA = BESTAA
      BB = BESTBB
      EPS = BSTEPS
C
      ZZ = -100000.0D0
      IF (NCA.GE.5 .OR. NCB.GE.5) THEN
          ZZ = 1.0D+19
          IF (LIMUP) ZZ = ELIMUP
      END IF
      SLFUN(9) = ZZ

      IF (PR) WRITE (T21,FMT=*) ' BEST AA,BB,TMID,EPS = ',BESTAA,BESTBB,
     +    BSTMID,BSTEPS
      IF (PR) WRITE (T21,FMT=*) ' BRS,BSTEIG,BSTEST = ',BRS,BSTEIG,
     +    BSTEST
      IF (PR) WRITE (T21,FMT=*) ' IFLAG = ',IFLAG
      IF (PR) WRITE (T21,FMT=*)
     +    '********************************************'
      IF (PR) WRITE (T21,FMT=*)
      IF (PR) WRITE (T21,FMT=*) ' EIGEST = '
      DO 1211 I = 1,NEIGST
          IF (PR) WRITE (T21,FMT=*) EIGEST(I,1),EIGEST(I,2),EIGEST(I,3)
 1211 CONTINUE
C-----------------------------------------------------------
C
      RETURN
      END
C
      SUBROUTINE SAVE(A,B,INTAB,P0ATA,QFATA,P0ATB,QFATB,A1,A2,B1,B2)

C  THIS PROGRAM SIMPLY SAVES THE CALLING ARGUMENTS IN
C    TWO COMMON BLOCKS FOR USE WHEREVER NEEDED.
C    IT IS CALLED BY SLEIGN.
C
C     .. Scalar Arguments ..
      REAL A,A1,A2,B,B1,B2,P0ATA,P0ATB,QFATA,QFATB
      INTEGER INTAB
C     ..
C     .. Scalars in Common ..
      REAL A1S,A2S,ASAV,B1S,B2S,BSAV,P0ATAS,P0ATBS,QFATAS,
     +                 QFATBS
      INTEGER INTSAV
C     ..
C     .. Common blocks ..
      COMMON /BCDATA/A1S,A2S,P0ATAS,QFATAS,B1S,B2S,P0ATBS,QFATBS
      COMMON /DATADT/ASAV,BSAV,INTSAV
C     ..
      ASAV = A
      BSAV = B
      INTSAV = INTAB
      P0ATAS = P0ATA
      QFATAS = QFATA
      P0ATBS = P0ATB
      QFATBS = QFATB
      A1S = A1
      A2S = A2
      B1S = B1
      B2S = B2
      IF (A1S.LT.0.0D0) THEN
          A1S = -A1S
          A2S = -A2S
      END IF
      IF (B1S.LT.0.0D0) THEN
          B1S = -B1S
          B2S = -B2S
      END IF
      RETURN
      END
C
      SUBROUTINE THUM(MF,ML,XS)
C     **********
C
C    THIS PROGRAM DETERMINES THE NUMBER OF CHANGES OF SIGN OF
C    THE BOUNDARY CONDITION FUNCTION U (OF THE PAIR (U,V))
C    WHICH THE USER SUPPLIES.  IT IS NEEDED ONLY IF ONE OF THE
C    ENDPOINTS OF THE INTERVAL (A,B) IS LCO.
C    IT IS CALLED BY SLEIGN.
C
C  YS IS LIKE XS, BUT HAS TWICE AS MANY POINTS.
C  MMW(N) IS THE VALUE OF THE INDEX I OF U(I), MF .LE. I .LE. 2*ML-2,
C    WHERE U FOR THE NTH TIME CHANGES SIGN FROM - TO +
C    AND WHERE P*U' IS POSITIVE.
C  MMWD IS THE NUMBER OF SUCH POINTS OF U.
C  THE QUANTITIES YS, MMW, MMWD ARE NEEDED BY SUBROUTINE SETTHU.
C
C  INPUT QUANTITIES: MF,ML,XS
C  OUTPUT QUANTITIES: YS,MMW,MMWD
C     **********
C     .. Scalars in Common ..
      INTEGER MMWD
C     ..
C     .. Arrays in Common ..
      REAL YS(200)
      INTEGER MMW(100)
C     ..
C     .. Local Scalars ..
C
      REAL PUP,PUP1,TMP1,TMP2,TMP3,TMP4,U,U1
      INTEGER I,N
C     ..
C     .. Common blocks ..
      COMMON /PASS/YS,MMW,MMWD
C     ..
C     .. Scalar Arguments ..
      INTEGER MF,ML
C     ..
C     .. Array Arguments ..
      REAL XS(*)
C     ..
C     .. External Subroutines ..
      EXTERNAL UV
C     ..
      DO 10 I = 1,99
          YS(2*I-1) = XS(I)
          YS(2*I) = 0.5D0* (XS(I)+XS(I+1))
   10 CONTINUE
      YS(199) = XS(100)
      N = 0
      U1 = 0.0D0
      PUP1 = 0.0D0
      DO 20 I = 2*MF - 1,2*ML - 1
          CALL UV(YS(I),U,PUP,TMP1,TMP2,TMP3,TMP4)
          IF (U1.LT.0.0D0 .AND. U.GT.0.0D0 .AND. PUP1.GT.0.0D0) THEN
              N = N + 1
              MMW(N) = I - 1
          END IF
          U1 = U
          PUP1 = PUP
   20 CONTINUE
      MMWD = N
      RETURN
      END
C
      SUBROUTINE AABB(TEND,OUT)

C  THIS PROGRAM IS FOR THE PURPOSE OF MOVING THE TRUNCATED ENDPOINTS,
C    AA OR BB, EITHER FURTHER OUT TOWARDS -1 OR +1, OR CLOSER IN,
C    DEPENDING ON THE SIGN OF OUT.
C    IT IS CALLED BY SLEIGN AND BY EFDATA.

C  INPUT QUANTITIES: TEND,OUT
C  OUTPUT QUANTITIES: TEND
C
C     .. Scalar Arguments ..
      REAL OUT,TEND
C     ..
C     .. Local Scalars ..
      REAL S,SEND
      INTEGER I,J
C     ..
C     .. Local Arrays ..
      REAL U(15)
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS
C     ..
      U(1) = 0.7D0
      U(2) = 0.8D0
      U(3) = 0.9D0
      U(4) = 0.95D0
      U(5) = 0.99D0
      U(6) = 0.999D0
      U(7) = 0.9999D0
      U(8) = 0.99999D0
      U(9) = 0.999999D0
      U(10) = 0.9999999D0
      U(11) = 1.0D0
c
      S = ABS(TEND)
      J = 9
      DO 10 I = 1,9
          IF (U(I).LT.S .AND. S.LE.U(I+1)) J = I
   10 CONTINUE
      IF (OUT.GT.0.0D0) THEN
          SEND = U(J+1)
          IF (S.EQ.SEND .AND. J.LT.10) SEND = U(J+2)
      ELSE
          SEND = U(J)
      END IF
      IF (SEND*TEND.LT.0.0D0) SEND = -SEND
      TEND = SEND
      RETURN
      END
C
      SUBROUTINE OBTAIN(NCA,ALFA,DTHDEA,NCB,BETA,DTHDEB,AA,BB,EIG,EIGPI,
     +                  TMID,EPS,DTHETA,DTHDE,ONEDIG,DTHDA,DTHDB,ER1,
     +                  ER1M,JFLAG)
C
C  THIS PROGRAM OBTAINS THE DIFFERENCE, DTHETA, BETWEEN THE VALUES
C    OF THETA OBTAINED BY INTEGRATING THE INITIAL VALUE PROBLEMS FOR
C    THETA FROM (OR NEAR) THE TWO ENDS OF THE INTERVAL (A,B) TO XMID.
C    THIS DIFFERENCE VANISHES WHEN THE VALUE BEING USED FOR THE
C     EIGENPARAMETER, EIG, IS EQUAL TO AN EIGENVALUE FOR THE PROBLEM.
C    IT IS CALLED BY SLEIGN.
C
C  INPUT QUANTITIES: NCA,ALFA,DTHDEA,NCB,BETA,DTHDEB,
C                    AA,BB,EIG,EIGPI,TMID,EPS,PI,TWOPI,HPI,
C                    A1,A2,P0ATA,QFATA,B1,B2,P0ATB,QFATB,
C                    A,B,INTAB,PR
C  OUTPUT QUANTITIES: EIGSAV,IND,ISAVE,TSAVEL,TSAVER,DTHDAA,
C                     DTHDBB,DTHETA,DTHDE,ONEDIG,DTHZ,
C                     MDTHZ,TSAVEL,TSAVER,ISAVE,JFLAG,
C                     DTHDA,DTHDB,ER1,ER1M
C     .. Scalar Arguments ..
      REAL AA,ALFA,BB,BETA,DTHDA,DTHDB,DTHDE,DTHDEA,DTHDEB,
     +                 DTHETA,EIG,EIGPI,EPS,ER1,ER1M,TMID
      INTEGER JFLAG,NCA,NCB
      LOGICAL ONEDIG
C     ..
C     .. Scalars in Common ..
      REAL A,A1,A2,AAS,B,B1,B2,BBS,DTHDAA,DTHDBB,EIGSAV,HPI,
     +                 P0ATA,P0ATB,PI,QFATA,QFATB,TMIDS,TSAVEL,TSAVER,
     +                 TWOPI
      INTEGER IND,INTAB,ISAVE,MDTHZ,T21
      LOGICAL ADDD,PR
C     ..
C     .. Local Scalars ..
      REAL ATHETA,C,DT,DTHZ,ER2,PX,QX,REMZ,THA,THB,WX,X
      INTEGER IFLAG
      LOGICAL AOK,BOK,LCA,LCB,LCOA,LCOB,SINGA,SINGB
C     ..
C     .. Local Arrays ..
      REAL ERL(3),ERR(3),YL(3),YR(3),YZL(3),YZR(3)
C     ..
C     .. External Functions ..
      REAL P,Q,W
      EXTERNAL P,Q,W
C     ..
C     .. External Subroutines ..
      EXTERNAL DXDT,INTEG
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS,COS,EXP,MAX,SIN
C     ..
C     .. Common blocks ..
      COMMON /BCDATA/A1,A2,P0ATA,QFATA,B1,B2,P0ATB,QFATB
      COMMON /DATADT/A,B,INTAB
      COMMON /DATAF/EIGSAV,IND
      COMMON /PIE/PI,TWOPI,HPI
      COMMON /PRIN/PR,T21
      COMMON /TDATA/AAS,TMIDS,BBS,DTHDAA,DTHDBB,MDTHZ,ADDD
      COMMON /TSAVE/TSAVEL,TSAVER,ISAVE
C     ..
      AOK = INTAB .LT. 3.D0 .AND. P0ATA .LT. 0.0D0 .AND.
     +      QFATA .GT. 0.0D0
      BOK = (INTAB.EQ.1 .OR. INTAB.EQ.3) .AND. P0ATB .LT. 0.0D0 .AND.
     +      QFATB .GT. 0.0D0
      LCA = NCA .EQ. 3 .OR. NCA .EQ. 4
      LCB = NCB .EQ. 3 .OR. NCB .EQ. 4
      LCOA = NCA .EQ. 4
      LCOB = NCB .EQ. 4
      SINGA = NCA .GE. 3
      SINGB = NCB .GE. 3
c     JFLAG = 0
      jflag = 0
      IND = 1
C  INITIALIZE TSAVEL, TSAVER, DTHDAA
      TSAVEL = -1.0D0
      TSAVER = 1.0D0
      IF (PR) WRITE (*,FMT=*) ' OBTAIN DTHETA '
      THA = ALFA
      DTHDAA = 0.0D0
      IF (SINGA .AND. .NOT.LCA) THEN
          CALL DXDT(AA,DT,X)
          PX = P(X)
          QX = Q(X)
          WX = W(X)
          C = EIG*WX - QX
          DTHDAA = - (COS(ALFA)**2/PX+C*SIN(ALFA)**2)*DT
C
C        TWO SPECIAL CASES FOR DTHDAA .
C
          IF (C.GE.0.0D0 .AND. P0ATA.LT.0.0D0 .AND.
     +        QFATA.LT.0.0D0) DTHDAA = DTHDAA + ALFA*DT/ (X-A)
          IF (C.GE.0.0D0 .AND. P0ATA.GT.0.0D0 .AND.
     +        QFATA.GT.0.0D0) DTHDAA = DTHDAA +
     +                                 (ALFA-0.5D0*PI)*DT/ (X-A)
      END IF

      THB = BETA
      DTHDBB = 0.0D0
      IF (SINGB .AND. .NOT.LCB) THEN
          CALL DXDT(BB,DT,X)
          PX = P(X)
          QX = Q(X)
          WX = W(X)
          C = EIG*WX - QX
          DTHDBB = - (COS(BETA)**2/PX+C*SIN(BETA)**2)*DT
C
C        TWO SPECIAL CASES FOR DTHDBB .
C
          IF (C.GE.0.0D0 .AND. P0ATB.LT.0.0D0 .AND.
     +        QFATB.LT.0.0D0) DTHDBB = DTHDBB + (PI-BETA)*DT/ (B-X)
          IF (C.GE.0.0D0 .AND. P0ATB.GT.0.0D0 .AND.
     +        QFATB.GT.0.0D0) DTHDBB = DTHDBB +
     +                                 (0.5D0*PI-BETA)*DT/ (B-X)
      END IF
C
C  PASS EIG TO SUBROUTINE F VIA THE COMMON/DATAF/ :
      EIGSAV = EIG

C  INTEGRATE FOR YL:
C
C     YL = (THETA,D(THETA)/D(EIG),D(THETA)/DA)
C
      YL(1) = 0.0D0
      YL(2) = 0.0D0
      YL(3) = 0.0D0
C
      ISAVE = 0
      CALL INTEG(AA,THA,DTHDAA,DTHDEA,TMID,A1,A2,EPS,YL,ERL,AOK,NCA,
     +           IFLAG)
      IF (IFLAG.EQ.5) THEN
          JFLAG = 51
          IF (PR) WRITE (T21,FMT=*) ' JFLAG = 51 '
          GO TO 130
      END IF

C  SET DTHDA:
      DTHDA = DTHDAA*EXP(-2.0D0*YL(3))

C  INTEGRATE FOR YR:
C
C     YR = (THETA,D(THETA)/D(EIG),D(THETA)/DB)
C
      YR(1) = 0.0D0
      YR(2) = 0.0D0
      YR(3) = 0.0D0
C
      ISAVE = 0
      CALL INTEG(BB,THB,DTHDBB,DTHDEB,TMID,B1,B2,EPS,YR,ERR,BOK,NCB,
     +           IFLAG)
      IF (IFLAG.EQ.5) THEN
          JFLAG = 52
          IF (PR) WRITE (T21,FMT=*) ' JFLAG = 52 '
          GO TO 130
      END IF

C  SET DTHDB:
      DTHDB = DTHDBB*EXP(-2.0D0*YR(3))
C
C  SET ER1, ER2, ER1M:
      ER1 = ERL(1) - ERR(1)
      ER2 = ERL(2) - ERR(2)
      ER1M = MAX(ABS(ERL(1)),ABS(ERR(1)))
C
C  IF EITHER ENDPOINT IS LCO, THEN INTEGRATIONS WITH EIG = 0
C    MUST ALSO BE CARRIED OUT.
C  THE VALUE OF THETA WHEN EIG=0 IS USED TO "CALIBRATE" THE
C  EIGENVALUE INDEXING -- IN THE LCO CASE, ONLY.
C  IN THIS CASE, SET ISAVE = 1 FOR USE IN SUBROUTINE LCO,
C    AND INTEGRATE FOR YZL:
C
      IF (LCOA .OR. LCOB) THEN
          EIGSAV = 0.0D0
          ISAVE = 1
          CALL INTEG(AA,THA,DTHDAA,DTHDEA,TMID,A1,A2,EPS,YZL,ERL,AOK,
     +               NCA,IFLAG)
          IF (IFLAG.EQ.5) THEN
              JFLAG = 53
              IF (PR) WRITE (T21,FMT=*) ' JFLAG = 53 '
              GO TO 130
          END IF
          ISAVE = 1
          CALL INTEG(BB,THB,DTHDBB,DTHDEB,TMID,B1,B2,EPS,YZR,ERR,BOK,
     +               NCB,IFLAG)
          IF (IFLAG.EQ.5) THEN
              JFLAG = 54
              IF (PR) WRITE (T21,FMT=*) ' JFLAG = 54 '
              GO TO 130
          END IF

          EIGSAV = EIG
C  SET DTHZ, MDTHZ:
          DTHZ = YZR(1) - YZL(1)
          MDTHZ = DTHZ/PI
          REMZ = DTHZ - MDTHZ*PI
          IF (DTHZ.LT.0.0D0 .AND. REMZ.LT.0.0D0) THEN
              MDTHZ = MDTHZ - 1
              REMZ = REMZ + PI
          END IF
          IF (REMZ.GT.3.14D0) MDTHZ = MDTHZ + 1
C  RESET ISAVE TO 0:
          ISAVE = 0
      END IF
C
C     DTHETA MEASURES THETA DIFFERENCE FROM LEFT AND RIGHT INTEGRATIONS.
C
C  SET DTHETA, DTHDE:
      DTHETA = YL(1) - YR(1) - EIGPI
      IF (LCOA .OR. LCOB) DTHETA = DTHETA + MDTHZ*PI
      DTHDE = YL(2) - YR(2)
      IF (PR) WRITE (T21,FMT=*) ' EIG = ',EIG
      IF (PR) WRITE (T21,FMT=*) ' YL(1),YR(1) = ',YL(1),YR(1)
      IF (PR) WRITE (T21,FMT=*) ' DTHETA,DTHDE = ',DTHETA,DTHDE
      IF (PR) WRITE (T21,FMT=*) ' MDTHZ = ',MDTHZ
C
      ATHETA = ABS(DTHETA)
      ONEDIG = (ABS(ER1).LE.0.5D0*ABS(DTHETA) .AND.
     +         ABS(ER2).LE.0.5D0*ABS(DTHDE)) .OR.
     +         MAX(ATHETA,ABS(ER1)) .LT. 1.0D-6

C     RIGHT AFTER LEAVING THIS SUBROUTINE:
C     WE NEED TO SET AAS = AA, AND BBS = BB ,
C     AND SET ATHETA = ABS(DTHETA).
C
C     WE ALSO NEED TO SET FIRSTT = .FALSE., UNLESS
C     JFLAG.EQ. 51, 52, 53, OR 54, WHEN WE NEED TO "GOTO 130"
C     AND SET          EXIT = .TRUE.
C     THIS SUBROUTINE NEEDS TO BE CALLED WITH THE INPUT FLAG
C     BEING CALLED "IFLAG" -- SO THAT EXTERNALLY THE "IFLAG"
C     WILL RETURN AS "IFLAG = 52, ETC".

  130 CONTINUE
      IFLAG = JFLAG
      RETURN
      END
C
      SUBROUTINE RESET(AA,BB,ALFA,BETA,EIG,GUESS,MF,ML,QS,WS,IMID,TMID,
     +                 LIMUP,ELIMUP,DTHDEA,DTHDEB,THELT0,EIGLO,EIGUP,
     +                 BRACKT,NCA,NCB,IFLAG)
C
C  THIS PROGRAM IS USED TO RESET THE MATCHING POINT, TMID, AND THE
C    BOUNDARY VALUES, ALFA & BETA, WHEN NECESSARY.  THESE
C    QUANTITIES CAN DEPEND ON THE CURRENT VALUE BEING USED FOR EIG.
C  IT IS CALLED BY SLEIGN.
C
C  INPUT QUANTITIES: AA,BB,ALFA,BETA,EIG,GUESS,MF,ML,
C                    QS,WS,IMID,TMID,LIMUP,ELIMUP,DTHDEA,DTHDEB,
C                    THELT0,EIGLO,EIGUP,BRACKT,NCA,NCB,
C                    A,B,INTAB,A1,A2,P0ATA,QFATA,B1,B2,P0ATB,QFATB,
C                    PI,TWOPI,HPI,EPSMIN,PR
C  OUTPUT QUANTITIES: IFLAG,ALFA,BETA,DTHDEA,DTHDEB,EIG,IMID,TMID,
C                    EIGUP
C
C     .. Scalar Arguments ..
      REAL AA,ALFA,BB,BETA,DTHDEA,DTHDEB,EIG,EIGLO,EIGUP,
     +                 ELIMUP,GUESS,TMID
      INTEGER IFLAG,IMID,MF,ML,NCA,NCB
      LOGICAL BRACKT,LIMUP,THELT0
C     ..
C     .. Array Arguments ..
      REAL QS(100),WS(100)
C     ..
C     .. Scalars in Common ..
      REAL A,A1,A2,B,B1,B2,EPSMIN,HPI,P0ATA,P0ATB,PI,QFATA,
     +                 QFATB,TWOPI
      INTEGER INTAB,T21
      LOGICAL PR
C     ..
C     .. Local Scalars ..
      REAL DERIVL,DERIVR,TMP1,TMP2,V
      INTEGER JFLAG,KFLAG
      LOGICAL LCA,LCB,SINGA,SINGB
C     ..
C     .. External Subroutines ..
      EXTERNAL ALFBET,SETMID
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC MIN
C     ..
C     .. Common blocks ..
      COMMON /BCDATA/A1,A2,P0ATA,QFATA,B1,B2,P0ATB,QFATB
      COMMON /DATADT/A,B,INTAB
      COMMON /PIE/PI,TWOPI,HPI
      COMMON /PRIN/PR,T21
      COMMON /RNDOFF/EPSMIN
C     ..
      LCA = NCA .EQ. 3 .OR. NCA .EQ. 4
      LCB = NCB .EQ. 3 .OR. NCB .EQ. 4
      SINGA = NCA .GE. 3
      SINGB = NCB .GE. 3

C        (SET-TMID-AND-BOUNDARY-CONDITIONS)
      IF (PR) WRITE (*,FMT=*) ' SET TMID AND BOUNDARY CONDITIONS '
      V = EIG*WS(IMID) - QS(IMID)
      IF (V.LE.0.0D0) CALL SETMID(MF,ML,EIG,QS,WS,IMID,TMID)

C        (RESET-BOUNDARY-CONDITIONS)
      DERIVL = 0.0D0
      IF (SINGA) CALL ALFBET(A,INTAB,AA,A1,A2,EIG,P0ATA,QFATA,.TRUE.,
     +                       ALFA,KFLAG,DERIVL)
      DTHDEA = DERIVL
      DERIVR = 0.0D0

      IF (SINGB) THEN
          CALL ALFBET(B,INTAB,BB,B1,B2,EIG,P0ATB,QFATB,.TRUE.,BETA,
     +                JFLAG,DERIVR)
          BETA = PI - BETA
      END IF
      DTHDEB = -DERIVR
      IF (PR) WRITE (T21,FMT='(A,E22.14,A,E22.14)') '  ALFA=',ALFA,
     +    '   BETA=',BETA
C
C     CHECK THAT BOUNDARY CONDITIONS CAN BE SATISFIED AT SINGULAR
C     ENDPOINTS.  IF NOT, TRY FOR SLIGHTLY ALTERED EIG CONSISTENT
C     WITH BOUNDARY CONDITIONS.
C     LIMUP = .TRUE. MEANS THAT THE EIGENVALUES ARE BOUNDED ABOVE,
C       BY APPROX. ELIMUP.
C
      IF (LIMUP .AND. EIG.NE.GUESS .AND. .NOT.BRACKT) THEN
          KFLAG = 1
          IF (SINGA .AND. .NOT.LCA) CALL ALFBET(A,INTAB,AA,A1,A2,EIG,
     +        P0ATA,QFATA,.TRUE.,TMP1,KFLAG,TMP2)
          JFLAG = 1
          IF (SINGB .AND. .NOT.LCB) CALL ALFBET(B,INTAB,BB,B1,B2,EIG,
     +        P0ATB,QFATB,.TRUE.,TMP1,JFLAG,TMP2)
          IFLAG = 1
          IF ((KFLAG.NE.1.OR.JFLAG.NE.1) .AND.
     +        (THELT0.AND.EIGLO.LT.ELIMUP)) THEN
              TMP1 = ELIMUP+2.0D0*EPSMIN
              EIGUP = MIN(TMP1,EIGUP)
c             EIGUP = MIN(ELIMUP+2.0D0*EPSMIN,EIGUP)
              IF (EIG.NE.EIGLO .AND. EIG.NE.EIGUP) THEN
                  EIG = 0.05D0*EIGLO + 0.95D0*EIGUP
              ELSE
                  IFLAG = 11
                  IF (PR) WRITE (T21,FMT=*) ' IFLAG = 11 '
                  GO TO 130
              END IF
          END IF
      END IF

C        (IMMEDIATELY UPON LEAVING THIS SUBROUTINE WE NEED
C         TO CHECK FOR IFLAG = 11.  IF SO, SET EXIT = .TRUE.
C         AND GOTO 130.)

  130 CONTINUE
      RETURN
      END
C
      SUBROUTINE DXDT(T,DT,X)
C     **********
C
C     THIS SUBROUTINE TRANSFORMS COORDINATES FROM T ON (-1,1) TO
C     X ON (A,B) .
C
C  INPUT QUANTITIES: T,A,B INTAB
C  OUTPUT QUANTITIES: DT,X
C     **********
C     .. Scalars in Common ..
      REAL A,B
      INTEGER INTAB
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS
C     ..
C     .. Common blocks ..
      COMMON /DATADT/A,B,INTAB
C     ..
C     .. Scalar Arguments ..
      REAL DT,T,X
C     ..
      GO TO (10,20,30,40) INTAB
   10 CONTINUE
      DT = 0.5D0* (B-A)
      X = 0.5D0* ((B+A)+ (B-A)*T)
      RETURN
   20 CONTINUE
      DT = 2.0D0/ (1.0D0-T)**2
      X = A + (1.0D0+T)/ (1.0D0-T)
      RETURN
   30 CONTINUE
      DT = 2.0D0/ (1.0D0+T)**2
      X = B - (1.0D0-T)/ (1.0D0+T)
      RETURN
   40 CONTINUE
      DT = 1.0D0/ (1.0D0-ABS(T))**2
      X = T/ (1.0D0-ABS(T))
      RETURN
      END
C
      REAL FUNCTION TFROMX(X)

C  THIS FUNCTION DETERMINES THE VALUE OF T IN (-1,1) WHICH
C    CORRESPONDS TO ANY VALUE OF X IN (A,B).

C  INPUT QUANTITIES: X,A,B,INTAB
C  OUTPUT QUANTITIES: TFROMX

C     .. Scalar Arguments ..
      REAL X
C     ..
C     .. Scalars in Common ..
      REAL A,B
      INTEGER INTAB
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS
C     ..
C     .. Common blocks ..
      COMMON /DATADT/A,B,INTAB
C     ..
      GO TO (10,20,30,40) INTAB
   10 CONTINUE
      TFROMX = (2.0D0*X- (B+A))/ (B-A)
      RETURN
   20 CONTINUE
      TFROMX = (X-A-1.0D0)/ (X-A+1.0D0)
      RETURN
   30 CONTINUE
      TFROMX = (1.0D0+X-B)/ (1.0D0-X+B)
      RETURN
   40 CONTINUE
      TFROMX = X/ (1.0D0+ABS(X))
      RETURN
      END
C
      REAL FUNCTION TFROMI(I)
C     **********
C
C  THIS FUNCTION DETERMINES THE VALUE OF THE VARIABLE T IN (-1,1)
C    WHICH CORRESPONDS TO THE SAMPLE POINT INDEX I.
C
C  INPUT QUANTITIES: I
C  OUTPUT QUANTITIES: TFROMI
C     **********
C     .. Scalar Arguments ..
      INTEGER I
C     ..
      IF (I.LT.8) THEN
          TFROMI = -1.0D0 + 0.1D0/4.0D0** (8-I)
      ELSE IF (I.GT.92) THEN
          TFROMI = 1.0D0 - 0.1D0/4.0D0** (I-92)
      ELSE
          TFROMI = 0.0227D0* (I-50)
      END IF
      RETURN
      END
C
      SUBROUTINE EXTRAP(T,TT,EIG,VALUE,DERIV,IFLAG)
C     **********
C
C     THIS SUBROUTINE IS CALLED FROM ALFBET IN DETERMINING BOUNDARY
C     VALUES AT A SINGULAR ENDPOINT OF THE INTERVAL FOR A
C     STURM-LIOUVILLE PROBLEM IN THE FORM
C
C       -(P(X)*Y'(X))' + Q(X)*Y(X) = EIG*W(X)*Y(X)  ON (A,B)
C
C     FOR USER-SUPPLIED COEFFICIENT FUNCTIONS P, Q, AND W.
C
C     EXTRAP, WHICH IN TURN CALLS INTPOL, EXTRAPOLATES THE FUNCTION
C
C        ARCTAN(1.0/SQRT(-P*(EIG*W-Q)))
C
C     FROM ITS VALUES FOR T WITHIN (-1,1) TO AN ENDPOINT.
C
C  INPUT QUANTITIES: T,TT,EIG,PR
C  OUTPUT QUANTITIES: VALUE,DERIV,IFLAG

C     SUBPROGRAMS CALLED
C
C       USER-SUPPLIED ..... P,Q,W
C
C       SLEIGN-SUPPLIED .. DXDT,INTPOL
C
C     **********
C     .. Local Scalars ..
      REAL ANS,CTN,ERROR,PROD,PX,QX,T1,TEMP,WX,X
      INTEGER KGOOD
C     ..
C     .. Local Arrays ..
      REAL FN1(5),XN(5)
C     ..
C     .. External Functions ..
      REAL P,Q,W
      EXTERNAL P,Q,W
C     ..
C     .. External Subroutines ..
      EXTERNAL DXDT,INTPOL
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS,ATAN,SQRT,TAN
C     ..
C     .. Scalar Arguments ..
      REAL DERIV,EIG,T,TT,VALUE
      INTEGER IFLAG
C     ..
C     .. Scalars in Common ..
      INTEGER T21
      LOGICAL PR
C     ..
C     .. Common blocks ..
      COMMON /PRIN/PR,T21
C     ..
C     (JUST TO MAKE SURE VALUE IS DEFINED, EVEN IF NOT
C     (NEEDED, SET IT TO 0.
      VALUE = 0.0D0
C
      IFLAG = 1
      KGOOD = 0
C  HERE, COMING FROM ALFBET, TT IS (PROBABLY) AA OR BB
      T1 = TT
      XN(1) = TT
   10 CONTINUE
      CALL DXDT(T1,TEMP,X)
      PX = P(X)
      QX = Q(X)
      WX = W(X)
      PROD = -PX* (EIG*WX-QX)
      IF (PROD.GT. (1.0D+10)) THEN
          ANS = 1.D-5
          GO TO 20
      END IF
      IF (PROD.LE.0.0D0) THEN
          T1 = 0.5D0* (T1+T)
          IF ((1.0D0+ (T1-T)**2).GT.1.0D0) GO TO 10
          IF (PROD.GT.-1.0D-6) VALUE = 2.0D0*ATAN(1.0D0)
          IFLAG = 5
          RETURN
      ELSE
          KGOOD = KGOOD + 1
          XN(KGOOD) = T1
          FN1(KGOOD) = ATAN(1.0D0/SQRT(PROD))
          T1 = 0.5D0* (T+T1)
          IF (KGOOD.LT.5) GO TO 10
      END IF

C  AT THIS POINT, THE XN(I) ARE VALUES OF T BETWEEN TT
C   (AA OR BB) AND 1.0 OR -1.0, OBTAINED BY BISECTION,
C   AND THE VALUES OF FN1 ARE THE CORRESPONDING VALUES
C   OF ATAN(1.0/SQRT(PROD)).
C  IN THIS CALL TO INTPOL, T IS 1.0 OR -1.0 AND
C    THE T1 SERVES AS ABSERR IN INTPOL.  THE RETURNED
C    VALUE ANS IS THE EXTRAPOLATED VALUE OF
C      ATAN(1.0/SQRT(PROD)) AT 1.0 OR -1.0   .

      IF (KGOOD.EQ.5) THEN
          IFLAG = 1
      ELSE
          IF (PR) WRITE (T21,FMT=*) ' IN EXTRAP, FAILED TO '
          IF (PR) WRITE (T21,FMT=*) ' GET 5 VALUES OF XN,FN '
      END IF
      T1 = 0.00001D0
      CALL INTPOL(5,XN,FN1,T,T1,3,ANS,ERROR)
   20 CONTINUE
      VALUE = ABS(ANS)
      CTN = 1.0D0/TAN(VALUE)
      DERIV = 0.5D0*PX*WX/CTN/ (1.0D0+CTN**2)
C  RESTORE TT TO ITS ORIGINAL VALUE.
      TT = XN(1)
      RETURN
      END
C
      SUBROUTINE INTPOL(N,XN,FN,X,ABSERR,MAXDEG,ANS,ERROR)
C     **********
C
C     THIS SUBROUTINE FORMS AN INTERPOLATING POLYNOMIAL FOR DATA PAIRS.
C     IT IS CALLED FROM EXTRAP AND FROM EXTR.
C     IT IS BEING USED TO EXTRAPOLATE THE VALUES FN AT POINTS XN
C     TO THE VALUE OF F AT X.
C     THE PARAMETER ABSERR IS THE REQUESTED ACCURACY IN ANS .
C
C   INPUT QUANTITIES: N,XN,FN,X,ABSERR,MAXDEG
C   OUTPUT QUANTITIES: ANS, ERROR
C     **********
C     .. Local Scalars ..
      REAL PROD
      INTEGER I,I1,II,IJ,IK,IKM1,J,K,L,LIMIT
C     ..
C     .. Local Arrays ..
      REAL V(10,10)
      INTEGER INDEX(10)
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS,MIN
C     ..
C     .. Scalar Arguments ..
      REAL ABSERR,ANS,ERROR,X
      INTEGER MAXDEG,N
C     ..
C     .. Array Arguments ..
      REAL FN(N),XN(N)
C     ..
      L = MIN(MAXDEG,N-2) + 2
      LIMIT = MIN(L,N-1)
      DO 10 I = 1,N
          V(I,1) = ABS(XN(I)-X)
          INDEX(I) = I
   10 CONTINUE
      DO 30 I = 1,LIMIT
          DO 20 J = I + 1,N
              II = INDEX(I)
              IJ = INDEX(J)
              IF (V(II,1).GT.V(IJ,1)) THEN
                  INDEX(I) = IJ
                  INDEX(J) = II
              END IF
   20     CONTINUE
   30 CONTINUE
      PROD = 1.0D0
      I1 = INDEX(1)
      ANS = FN(I1)
      V(1,1) = FN(I1)
      DO 50 K = 2,L
          IK = INDEX(K)
          V(K,1) = FN(IK)
          DO 40 I = 1,K - 1
              II = INDEX(I)
              V(K,I+1) = (V(I,I)-V(K,I))/ (XN(II)-XN(IK))
   40     CONTINUE
          IKM1 = INDEX(K-1)
          PROD = (X-XN(IKM1))*PROD
          ERROR = PROD*V(K,K)
          IF (ABS(ERROR).LE.ABSERR) RETURN
          ANS = ANS + ERROR
   50 CONTINUE
      ANS = ANS - ERROR
      RETURN
      END
C
      REAL FUNCTION EPSLON(X)
C
C     ESTIMATE UNIT ROUNDOFF IN QUANTITIES OF SIZE X.
C
C     THIS PROGRAM SHOULD FUNCTION PROPERLY ON ALL SYSTEMS
C     SATISFYING THE FOLLOWING TWO ASSUMPTIONS,
C        1.  THE BASE USED IN REPRESENTING FLOATING POINT
C            NUMBERS IS NOT A POWER OF THREE.
C        2.  THE QUANTITY  A  IN STATEMENT 10 IS REPRESENTED TO
C            THE ACCURACY USED IN FLOATING POINT VARIABLES
C            THAT ARE STORED IN MEMORY.
C     THE STATEMENT NUMBER 10 AND THE GO TO 10 ARE INTENDED TO
C     FORCE OPTIMIZING COMPILERS TO GENERATE CODE SATISFYING
C     ASSUMPTION 2.
C     UNDER THESE ASSUMPTIONS, IT SHOULD BE TRUE THAT,
C            A  IS NOT EXACTLY EQUAL TO FOUR-THIRDS,
C            B  HAS A ZERO FOR ITS LAST BIT OR DIGIT,
C            C  IS NOT EXACTLY EQUAL TO ONE,
C            EPS  MEASURES THE SEPARATION OF 1.0 FROM
C                 THE NEXT LARGER FLOATING POINT NUMBER.
C
C     .. Scalar Arguments ..
      REAL X
C     ..
C     .. Local Scalars ..
      REAL A,B,C,EPS,FOUR,THREE
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS
C     ..
      FOUR = 4.0D0
      THREE = 3.0D0
      A = FOUR/THREE
   10 B = A - 1.0D0
      C = B + B + B
      EPS = ABS(C-1.0D0)
      IF (EPS.EQ.0.0D0) GO TO 10
      EPSLON = EPS*ABS(X)
      RETURN
      END
C
      SUBROUTINE ESTPAC(IOSC,MF,ML,EEE,SUM0,QS,WS,DS,DELT,PS,PSS,TAU,
     +                  JJL,JJR,SUM,U)
C     **********
C
C     THIS SUBROUTINE ESTIMATES THE CHANGE IN 'PHASE ANGLE' IN THE
C     EIGENVALUE DETERMINATION OF A STURM-LIOUVILLE PROBLEM IN THE FORM
C
C       -(P(X)*Y'(X))' + Q(X)*Y(X) = EIG*W(X)*Y(X)  ON (A,B)
C
C     FOR USER-SUPPLIED COEFFICIENT FUNCTIONS P, Q, AND W.
C
C     THE SUBROUTINE APPROXIMATES (BY TRAPEZOIDAL RULE) THE INTEGRAL OF
C
C        SQRT((EIG*W-Q)/P)
C
C     WHERE THE INTEGRAL IS TAKEN OVER THOSE X IN (A,B) FOR WHICH
C
C        (EIG*W-Q)/P .GT. 0
C
C     ESTPAC IS CALLED BY SUBROUTINES ESTIM AND PRELIM.
C
C   IN THE MAIN IF..THEN..ELSE, THE FIRST PART DOES THE ESTIMATING OF
C     THE PHASE ANGLE CHANGE, AND THE SECOND PART DETERMINES THE ZEE'S
C     TO BE USED IN THE PRUFER TRANSFORMATION BY SUBROUTINE GERKZ.
C
C  INPUT QUANTITIES: IOSC,MF,ML,EEE,SUM0,QS,WS,DS,DELT,PS,TAU
C  OUTPUT QUANTITIES: PSS,JJL,JJR,SUM,U,JAY,ZEE
C     **********
C     .. Local Scalars ..
      REAL C1,C2,C3,C4,DPSUM,DPSUMT,ONE,PSUM,RT,RTSAV,UT,V,
     +                 WSAV,WW,ZAV,ZAVJ
      INTEGER J,JJ,MF1
C     ..
C     .. Arrays in Common ..
      REAL ZEE(100)
      INTEGER JAY(100)
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS,MAX,MIN,MOD,SIGN,SQRT
C     ..
C     .. Common blocks ..
      COMMON /ZEEZ/JAY,ZEE
C     ..
C     .. Scalar Arguments ..
      REAL EEE,SUM,SUM0,TAU,U
      INTEGER JJL,JJR,MF,ML
      LOGICAL IOSC
C     ..
C     .. Array Arguments ..
      REAL DELT(100),DS(100),PS(100),PSS(100),QS(100),
     +                 WS(100)
C     ..
      ONE = 1.0D0
      C1 = 0.1D0
      C2 = 0.25D0
      C3 = 0.5D0
      C4 = 0.001D0
C
C     SUM ACCUMULATES THE INTEGRAL APPROXIMATION.  U MEASURES THE TOTAL
C     LENGTH OF SUBINTERVALS WHERE (EIG*W-Q)/P .GT. 0.0.  ZAV IS THE
C     AVERAGE VALUE OF SQRT((EIG*W-Q)*P) OVER THOSE SUBINTERVALS.
C
      IF (.NOT.IOSC) THEN

          DO 5 J = 1,100
              PSS(J) = 0.0D0
    5     CONTINUE

          JJL = 99
          JJR = 1
          SUM = 0.0D0
          U = 0.0D0
          UT = 0.0D0
          ZAV = 0.0D0
          WSAV = EEE*WS(MF) - QS(MF)
          IF (WSAV.GT.0.0D0) THEN
              RTSAV = SIGN(SQRT(WSAV),PS(MF))
          ELSE
              RTSAV = 0.0D0
          END IF
          DO 10 J = MF + 1,ML
              WW = EEE*WS(J) - QS(J)
              IF (WW.GT.0.0D0) THEN
                  U = U + DS(J-1)
                  UT = UT + DELT(J-1)
                  RT = SIGN(SQRT(WW),PS(J))
              ELSE
                  RT = 0.0D0
              END IF
              IF (WW.EQ.0.0D0 .OR. WSAV.EQ.0.0D0 .OR.
     +            WW.EQ.SIGN(WW,WSAV)) THEN
                  V = RT + RTSAV
              ELSE
                  V = (WW*RT+WSAV*RTSAV)/ABS(WW-WSAV)
              END IF
              WSAV = WW
              RTSAV = RT
              PSUM = DS(J-1)*V
              IF (EEE.EQ.0.0D0) THEN
                  PSS(J) = PSUM
              ELSE
                  DPSUM = PSUM - PSS(J)
                  DPSUMT = DPSUM*DELT(J-1)/DS(J-1)
                  IF (DPSUMT.GT.C4*TAU) THEN
                      JJL = MIN(JJL,J)
                      JJR = MAX(JJR,J)
                  END IF
              END IF
              SUM = SUM + PSUM
              IF (U.GT.0.0D0) ZAV = ZAV + DELT(J-1)*V*ABS(PS(J)+PS(J-1))
   10     CONTINUE
          SUM = C3*SUM - SUM0
          ZAV = C2*ZAV
      ELSE
          JJ = 1
          JAY(1) = MF
   20     CONTINUE
          SUM = 0.0D0
          U = 0.0D0
          UT = 0.0D0
          ZAV = 0.0D0
          ZAVJ = 0.0D0
          MF1 = JAY(JJ)
          WSAV = EEE*WS(MF1) - QS(MF1)
          IF (WSAV.GT.0.0D0) THEN
              RTSAV = SIGN(SQRT(WSAV),PS(MF1))
          ELSE
              RTSAV = 0.0D0
          END IF
          DO 30 J = MF1 + 1,ML
              WW = EEE*WS(J) - QS(J)
              IF (WW.GT.0.0D0) THEN
                  U = U + DS(J-1)
                  UT = UT + DELT(J-1)
                  RT = SIGN(SQRT(WW),PS(J))
              ELSE
                  RT = 0.0D0
              END IF
              IF (WW.EQ.0.0D0 .OR. WSAV.EQ.0.0D0 .OR.
     +            WW.EQ.SIGN(WW,WSAV)) THEN
                  V = RT + RTSAV
              ELSE
                  V = (WW*RT+WSAV*RTSAV)/ABS(WW-WSAV)
              END IF
              WSAV = WW
              RTSAV = RT
              PSUM = DS(J-1)*V
              SUM = SUM + PSUM
              IF (U.GT.0.0D0) ZAV = ZAV + DELT(J-1)*V*ABS(PS(J)+PS(J-1))
              IF (U.NE.0.0D0) THEN
                  ZAVJ = C2*ZAV/UT
                  IF ((MOD(J-MF1,7).EQ.0) .OR. J.EQ.ML) THEN
                      JJ = JJ + 1
                      JAY(JJ) = J
                      ZEE(JJ) = ZAVJ
                      IF (ZEE(JJ).NE.0.0D0) ZEE(JJ) = MAX(ZAVJ,C1)
                      IF (ZEE(JJ).EQ.0.0D0) ZEE(JJ) = ONE
                      GO TO 40
                  END IF
              END IF
   30     CONTINUE
   40     CONTINUE
          IF (J.GT.ML) THEN
              JJ = JJ + 1
              JAY(JJ) = ML
              ZEE(JJ) = ZAVJ
              IF (ZEE(JJ).NE.0.0D0) ZEE(JJ) = MAX(ZAVJ,C1)
              IF (ZEE(JJ).EQ.0.D0) ZEE(JJ) = ONE
          END IF
          IF (J.LT.ML) GO TO 20
          SUM = C3*SUM
          ZAV = C2*ZAV
      END IF
      RETURN
      END
C
      SUBROUTINE ALFBET(XEND,INTAB,TT,COEF1,COEF2,EIG,P0,QF,SING,VALUE,
     +                  IFLAG,DERIV)
C     **********
C
C  THIS SUBROUTINE COMPUTES A BOUNDARY VALUE OF THE PRUEFER ANGLE,
C    THETA, FOR A SPECIFIED ENDPOINT OF THE INTERVAL FOR A
C    STURM-LIOUVILLE PROBLEM IN THE FORM
C
C       -(P(X)*Y'(X))' + Q(X)*Y(X) = EIG*W(X)*Y(X)  ON (A,B)
C
C     FOR USER-SUPPLIED COEFFICIENT FUNCTIONS P, Q, AND W.  IT IS CALLED
C     FROM RESET AND PRELIM.
C     BOTH REGULAR AND SINGULAR ENDPOINTS ARE TREATED.
C
C     SUBPROGRAMS CALLED
C
C       USER-SUPPLIED ..... P,Q,W
C
C       SLEIGN-SUPPLIED .. DXDT,EXTRAP
C  INPUT QUANTITIES: XEND,INTAB,TT,COEF1,COEF2,EIG,P0,QF,
C                    SING,PI,TWOPI,HPI
C  OUTPUT QUANTITIES: VALUE,DERIV,IFLAG
C
C     **********
C     .. Local Scalars ..
      REAL C,CD,D,HH,ONE,PX,QX,T,TEMP,TTS,WX,X
      LOGICAL LOGIC
C     ..
C     .. External Functions ..
      REAL P,Q,W
      EXTERNAL P,Q,W
C     ..
C     .. External Subroutines ..
      EXTERNAL DXDT,EXTRAP
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS,ATAN,SIGN,SQRT
C     ..
C     .. Common blocks ..
      COMMON /PIE/PI,TWOPI,HPI
C     ..
C     SET MACHINE DEPENDENT CONSTANT.
C
C     .. Scalar Arguments ..
      REAL COEF1,COEF2,DERIV,EIG,P0,QF,TT,VALUE,XEND
      INTEGER IFLAG,INTAB
      LOGICAL SING
C     ..
C     .. Scalars in Common ..
      REAL HPI,PI,TWOPI
C     ..
      ONE = 1.0D0
      IFLAG = 1
      DERIV = 0.0D0
      IF (.NOT.SING) THEN
          VALUE = 0.5D0*PI
          IF (COEF1.NE.0.0D0) VALUE = ATAN(-COEF2/COEF1)
          LOGIC = (TT.LT.0.0D0 .AND. VALUE.LT.0.0D0) .OR.
     +            (TT.GT.0.0D0 .AND. VALUE.LE.0.0D0)
          IF (LOGIC) VALUE = VALUE + PI
      ELSE
          LOGIC = (INTAB.EQ.2 .AND. TT.GT.0.0D0) .OR.
     +            (INTAB.EQ.3 .AND. TT.LT.0.0D0) .OR. INTAB .EQ. 4 .OR.
     +            (P0.GT.0.0D0 .AND. QF.LT.0.0D0)
          IF (LOGIC) THEN
              T = SIGN(ONE,TT)
              TTS = TT

C  HERE, T IS 1.0 OR -1.0, AND
C    TTS IS AA OR BB (PROBABLY)
C  THIS CALL TO EXTRAP EXTRAPOLATES THE VALUES OF
C    ATAN(1.0/SQRT(-P(EIG*W - Q)) AT POINTS T BETWEEN
C    AA OR BB AND -1.0 OR 1.0 TO THE ENDPOINT, -1.0 OR 1.0
C    IT ALSO RETURNS THE ASSOCIATED VALUE OF DERIV AT THAT
C      ENDPOINT.  DERIV IS D(VALUE)/D(EIG) .

              CALL EXTRAP(T,TTS,EIG,VALUE,DERIV,IFLAG)
          ELSE
              CALL DXDT(TT,TEMP,X)
              PX = P(X)
              QX = Q(X)
              WX = W(X)
              C = 2.0D0* (EIG*WX-QX)
              IF (C.LT.0.0D0) THEN
                  VALUE = 0.0D0
                  IF (P0.GT.0.0D0) VALUE = 0.5D0*PI
              ELSE
                  HH = ABS(XEND-X)
                  D = 2.0D0*HH/PX
                  CD = C*D*HH
                  IF (P0.GT.0.0D0) THEN
                      VALUE = C*HH
                      IF (CD.LT.1.0D0) VALUE = VALUE/
     +                    (1.0D0+SQRT(1.0D0-CD))
                      VALUE = VALUE + 0.5D0*PI
                  ELSE
                      VALUE = D
                      IF (CD.LT.1.0D0) VALUE = VALUE/
     +                    (1.0D0+SQRT(1.0D0-CD))
                  END IF
              END IF
          END IF
      END IF
      RETURN
      END
C
      SUBROUTINE SAMPLE(NCA,NCB,MF,ML,XS,PS,QS,WS,DS,DELT,EMIN,EMAX,
     +                  IMIN,IMAX,AAA,BBB)

C  THIS PROGRAM IS FOR THE PURPOSE OF SAMPLING THE COEFFICIENT
C    FUNCTIONS P, Q, W, PRIMARILY IN ORDER TO BE ABLE TO MAKE A
C    FIRST ESTIMATE OF THE DESIRED EIGENVALUE.
C    IT IS CALLED FROM SLEIGN.
C  THE ARRAY XS CONTAINS THE VALUES OF X IN (A,B) WHICH ARE
C    USED IN THE SAMPLING PROCESS.  THE ARRAY PS CONTAINS THE
C    CORRESPONDING VALUES OF P.  BUT THE ARRAYS QS AND WS CONTAIN
C    THE CORRESPONDING VALUES NOT OF Q AND W BUT OF
C    Q/P AND W/P.
C  ARRAYS DS AND DELT CONTAIN THE CORRESPONDING SAMPLING POINT
C    INTERVALS OF X AND T.
C  ALSO COMPUTED ARE MINIMUM AND MAXIMUM VALUES OF Q/W, CALLED
C    EMIN AND EMAX, WHICH OCCUR AT THE INDEX I EQUAL TO IMIN AND IMAX.
C  MF AND ML ARE THE FIRST AND LAST VALUES OF THE INDEX I,
C    1 .LE. MF .LT. ML .LE. 99.
C
C  INPUT QUANTITIES: NCA,NCB,EPSMIN
C  OUTPUT QUANTITIES: MF,ML,XS,PS,QS,WS,DS,DELT,EMIN,EMAX,
C                 IMIN,IMAX,AAA,BBB

C     .. Scalar Arguments ..
      REAL AAA,BBB,EMAX,EMIN
      INTEGER IMAX,IMIN,MF,ML,NCA,NCB
C     ..
C     .. Array Arguments ..
      REAL DELT(100),DS(100),PS(100),QS(100),WS(100),XS(100)
C     ..
C     .. Scalars in Common ..
      REAL EPSMIN
C     ..
C     .. Local Scalars ..
      REAL ONE,PX,QX,T,T50,THRESH,TMP,TS,WX,X,X50,XSAV
      INTEGER I
      LOGICAL FYNYT,LCOA,LCOB,SINGA,SINGB
C     ..
C     .. External Functions ..
      REAL P,Q,TFROMI,W
      EXTERNAL P,Q,TFROMI,W
C     ..
C     .. External Subroutines ..
      EXTERNAL DXDT,THUM
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS
C     ..
C     .. Common blocks ..
      COMMON /RNDOFF/EPSMIN
C     ..
      SINGA = NCA .GE. 3
      SINGB = NCB .GE. 3
      LCOA = NCA .EQ. 4
      LCOB = NCB .EQ. 4
      ONE = 1.0D0
C     DO (SAMPLE-COEFFICIENTS)
      THRESH = 1.D+35
      IF (NCA.GE.5) THRESH = 1.D+17
      T50 = 0.0D0
      CALL DXDT(T50,TMP,X50)
      XS(50) = X50
      TS = T50
      PX = P(X50)
      QX = Q(X50)
      WX = W(X50)
      PS(50) = PX
      QS(50) = QX/PX
      WS(50) = WX/PX
C
C     EMIN = MIN(Q/W), ACHIEVED AT X FOR INDEX VALUE IMIN.
C     EMAX = MAX(Q/W), ACHIEVED AT X FOR INDEX VALUE IMAX.
C     MF AND ML ARE THE LEAST AND GREATEST INDEX VALUES, RESPECTIVELY.
C
      XSAV = X50
      EMIN = 0.0D0
      IF (QX.NE.0.0D0) EMIN = QX/WX
      EMAX = EMIN
      IMIN = 50
      IMAX = 50
      DO 20 I = 49,1,-1
          T = TFROMI(I)
          CALL DXDT(T,TMP,X)
          XS(I) = X
          PX = P(X)
          QX = Q(X)
          WX = W(X)
          PS(I) = PX
          QS(I) = QX/PX
          WS(I) = WX/PX
          DS(I) = XSAV - X
          DELT(I) = 0.5D0* (TS-T)
          XSAV = X
          TS = T
C
C     TRY TO AVOID OVERFLOW BY STOPPING WHEN FUNCTIONS ARE LARGE NEAR A
C     OR WHEN W IS SMALL NEAR A.
C
          FYNYT = (ABS(WX)+ABS(QX)+1.0D0/ABS(PX)) .LE. THRESH .AND.
     +            WX .GT. EPSMIN
          IF (QX.NE.0.0D0 .AND. QX/WX.LT.EMIN) THEN
              EMIN = QX/WX
              IMIN = I
          END IF
          IF (QX.NE.0.0D0 .AND. QX/WX.GT.EMAX) THEN
              EMAX = QX/WX
              IMAX = I
          END IF
          MF = I
          IF (.NOT.FYNYT) GO TO 30
   20 CONTINUE
      DO 25 I = 0,-14,-1
          T = TFROMI(I)
          IF (T.LE.-ONE+EPSMIN) GO TO 30
          CALL DXDT(T,TMP,X)
          PX = P(X)
          QX = Q(X)
          WX = W(X)
          FYNYT = (ABS(WX)+ABS(QX)+1.0D0/ABS(PX)) .LE. THRESH .AND.
     +            WX .GT. EPSMIN
          IF (.NOT.FYNYT) GO TO 30
   25 CONTINUE
   30 CONTINUE
      THRESH = 1.D+35
      IF (NCB.GE.5) THRESH = 1.D+17
      AAA = T
      IF (.NOT.SINGA) AAA = -ONE
      TS = T50
      XSAV = X50
      DO 40 I = 51,99
          T = TFROMI(I)
          CALL DXDT(T,TMP,X)
          XS(I) = X
          PX = P(X)
          QX = Q(X)
          WX = W(X)
          PS(I) = PX
          QS(I) = QX/PX
          WS(I) = WX/PX
          DS(I-1) = X - XSAV
          DELT(I-1) = 0.5D0* (T-TS)
          XSAV = X
          TS = T
C
C     TRY TO AVOID OVERFLOW BY STOPPING WHEN FUNCTIONS ARE LARGE NEAR B
C     OR WHEN W IS SMALL NEAR B.
C
          FYNYT = (ABS(QX)+ABS(WX)+1.0D0/ABS(PX)) .LE. THRESH .AND.
     +            WX .GT. EPSMIN
          IF (QX.NE.0.0D0 .AND. QX/WX.LT.EMIN) THEN
              EMIN = QX/WX
              IMIN = I
          END IF
          IF (QX.NE.0.0D0 .AND. QX/WX.GT.EMAX) THEN
              EMAX = QX/WX
              IMAX = I
          END IF
          ML = I - 1
          IF (.NOT.FYNYT) GO TO 50
   40 CONTINUE
      XS(100) = 0.0D0
      DO 45 I = 100,114
          T = TFROMI(I)
          IF (T.GE.ONE-EPSMIN) GO TO 50
          CALL DXDT(T,TMP,X)
          PX = P(X)
          QX = Q(X)
          WX = W(X)
          FYNYT = (ABS(WX)+ABS(QX)+1.0D0/ABS(PX)) .LE. THRESH .AND.
     +            WX .GT. EPSMIN
          IF (.NOT.FYNYT) GO TO 50
   45 CONTINUE
   50 CONTINUE
      BBB = T
      IF (.NOT.SINGB) BBB = ONE
C
      IF (LCOA .OR. LCOB) CALL THUM(MF,ML,XS)
C
      RETURN
      END
C
      SUBROUTINE ESTIM(IOSC,PIN,MF,ML,PS,QS,WS,PSS,DS,DELT,TAU,JJL,JJR,
     +                 SUM,LIMUP,ELIMUP,EMAX,IMAX,IMIN,EL,WL,DEDW,
     +                 BALLPK,EEE,JFLAG)

C  THIS PROGRAM MAKES A FIRST ESTIMATE OF THE DESIRED
C    EIGENVALUE, USING THE ARRAYS PS,QS,WS OF SAMPLED VALUES OF
C    THE COEFFICIENT FUNCTIONS P,Q,W OBTAINED BY SUBROUTINE SAMPLE.
C    IT IS CALLED BY SLEIGN.
C    THIS ESTIMATE IS RETURNED IN THE VARIABLE EEE.
C  JJL AND JJR ARE THE MIN AND MAX OF THE INDEX I FOR WHICH
C    EIG*W - Q IS NONNEGATIVE (OF THE SAMPLED VALUES).
C  IMIN AND IMAX ARE THE SAMPLE INDICES WHERE THE FUNCTION QS/WS
C  ATTAINS ITS MINIMUM AND MAXIMUM OF EMIN AND EMAX.
C  WHEN THE ESTIMATING PROCESS BREAKS DOWN, A "BALLPARK" ESTIMATE,
C    BALLPK, IS DETERMINED.
C
C  BASICALLY, THE ESTIMATE, EEE, IS THE SOLUTION OF THE
C  EQUATION
C                  |b
C          INTEGRAL|SQRT((EEE*W-Q)/P)*DX = (NUMEIG+1)*PI
C                  |a
C
C  WHERE THE INTEGRATION IS OVER THOSE X FOR WHICH THE INTEGRAND
C  IS REAL.
C  INPUT QUANTITIES: IOSC,PIN,MF,ML,PS,QS,WS,DS,DELT,TAU,
C                    JJL,JJR,LIMUP,ELIMUP,EMIN,EMAX,IMIN,IMAX
C  OUTPUT QUANTITIES: PSS,JJL,JJR,SUM,IMIN,IMAX,EL,WL,DEDW,
C                     BALLPK,EEE,JFLAG
C
C     .. Scalar Arguments ..
      REAL BALLPK,DEDW,EEE,EL,ELIMUP,EMAX,PIN,SUM,TAU,WL
      INTEGER IMAX,IMIN,JFLAG,JJL,JJR,MF,ML
      LOGICAL IOSC,LIMUP
C     ..
C     .. Array Arguments ..
      REAL DELT(100),DS(100),PS(100),PSS(100),QS(100),
     +                 WS(100)
C     ..
C     .. Scalars in Common ..
      REAL A,B
      INTEGER INTAB,T21
      LOGICAL PR
C     ..
C     .. Local Scalars ..
      REAL EU,FNEW,FOLD,ONE,SUM0,U,ULO,UUP,WU
      INTEGER IE,JLOOP
      LOGICAL LOGIC
C     ..
C     .. External Subroutines ..
      EXTERNAL ESTPAC
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS,MAX,MIN
C     ..
C     .. Common blocks ..
      COMMON /DATADT/A,B,INTAB
      COMMON /PRIN/PR,T21
C     ..
      JFLAG = 1
      ONE = 1.0D0
      SUM0 = 0.0D0
C           ------------------------
C  WHEN ESTPAC IS CALLED, THE RETURNED VALUE OF SUM (OR SUM0) IS
C  THE APPROXIMATION TO THE ABOVE INTEGRAL,
C  THE VALUE OF U IS THE TOTAL LENGTH OF THOSE SUBINTERVALS FOR
C  WHICH THE INTEGRAND IS REAL.

      IF (IOSC) THEN
          EEE = 0.0D0
          CALL ESTPAC(.FALSE.,MF,ML,EEE,SUM0,QS,WS,DS,DELT,PS,PSS,TAU,
     +                JJL,JJR,SUM,U)
          SUM0 = SUM
      END IF
C           ------------------------
      EEE = MIN(ELIMUP,EMAX)
      CALL ESTPAC(.FALSE.,MF,ML,EEE,SUM0,QS,WS,DS,DELT,PS,PSS,TAU,JJL,
     +            JJR,SUM,U)
C           ------------------------
      EU = EEE
      WU = SUM
      UUP = U
C          --------------------------
C   PIN IS NORMALLY SET = (NUMEIG+1)*PI, EXCEPT WHEN IOSC = .TRUE.

   55 CONTINUE
      IF (.NOT.LIMUP .AND. ABS(SUM).GE.10.0D0*MAX(ONE,ABS(PIN))) THEN
          IF (SUM.GE.10.0D0*PIN) THEN
              IF (EEE.GE.ONE) THEN
                  EEE = EEE/10.0D0
              ELSE IF (EEE.LT.-ONE) THEN
                  EEE = 10.0D0*EEE
              ELSE
                  EEE = EEE - ONE
              END IF
          ELSE
              IF (EEE.LE.-ONE) THEN
                  EEE = EEE/10.0D0
              ELSE IF (EEE.GT.ONE) THEN
                  EEE = 10.0D0*EEE
              ELSE
                  EEE = EEE + ONE
              END IF
          END IF
          CALL ESTPAC(.FALSE.,MF,ML,EEE,SUM0,QS,WS,DS,DELT,PS,PSS,TAU,
     +                JJL,JJR,SUM,U)
          GO TO 55
      END IF
C           -------------------------------
C   THE LOCAL VARIABLES EL, WL ARE USED TO INDICATE VALUE OF EIG
C   AND VALUE OF SUM WHEN SUM .LT. PIN;
C   EU, WU ARE USED TO INDICATE VALUE OF EIG AND VALUE OF SUM
C   WHEN SUM .GT. PIN.
C   ULO, UUP ARE CORRESPONDING VALUES FOR U, THE TOTAL LENGTH
C   OF INTERVAL FOR WHICH THE INTEGRAND IS POSITIVE.
      EU = EEE
      WU = SUM
      UUP = U
      IF (SUM.GE.PIN) THEN
          EL = EU
          WL = WU
          ULO = UUP
C           ----------------------
          JLOOP = 0
   60     CONTINUE
          IF (WL.GE.PIN) THEN
C           REDUCE EEE:
              EU = EL
              WU = WL
              UUP = ULO
              EEE = EL - ((WL-PIN+3.0D0)/U)**2 - ONE
              CALL ESTPAC(.FALSE.,MF,ML,EEE,SUM0,QS,WS,DS,DELT,PS,PSS,
     +                    TAU,JJL,JJR,SUM,U)
              EL = EEE
              WL = SUM
              ULO = U
              IF (SUM.EQ.0.0D0) JLOOP = JLOOP + 1
              IF (JLOOP.GE.5) THEN
C                   (ESTIMATOR FAILURE)
                  JFLAG = 0
                  RETURN
              END IF
              GO TO 60
          END IF
C           -----------------------

      ELSE
C       INCREASE EEE:
          EL = EEE
          WL = SUM
          ULO = U
      END IF
      IF (LIMUP .AND. WU.LT.PIN) THEN
          EEE = ELIMUP
      ELSE
          IF (UUP.EQ.0.0D0) THEN
              EEE = EMAX + ONE
              CALL ESTPAC(.FALSE.,MF,ML,EEE,SUM0,QS,WS,DS,DELT,PS,PSS,
     +                    TAU,JJL,JJR,SUM,U)
              EU = EEE
              WU = SUM
              UUP = U
          END IF
   70     CONTINUE
          IF (WU.LE.PIN) THEN
C           INCREASE EEE:
              EL = EU
              WL = WU
              ULO = UUP
              EEE = EU + ((PIN-WU+3.0D0)/U)**2 + ONE
              CALL ESTPAC(.FALSE.,MF,ML,EEE,SUM0,QS,WS,DS,DELT,PS,PSS,
     +                    TAU,JJL,JJR,SUM,U)
              EU = EEE
              WU = SUM
              UUP = U
              GO TO 70
          END IF
C            --------------------------------
   80     CONTINUE
C    DETERMINE THE INDICES IMIN, IMAX, WHERE THE FUNCTION
C    QS/WS ATTAINS ITS MINIMUM AND MAXIMUM VALUES OF EMIN, EMAX:
          IF (ABS(IMAX-IMIN).GE.2 .AND. EU.LE.EMAX) THEN
              IE = (IMAX+IMIN)/2
              EEE = QS(IE)/WS(IE)
              CALL ESTPAC(.FALSE.,MF,ML,EEE,SUM0,QS,WS,DS,DELT,PS,PSS,
     +                    TAU,JJL,JJR,SUM,U)
              IF (SUM.GT.PIN) THEN
                  IMAX = IE
                  WU = SUM
                  EU = EEE
                  UUP = U
              ELSE
                  IMIN = IE
                  WL = SUM
                  EL = EEE
                  ULO = U
              END IF
              GO TO 80
          END IF
C            ---------------------------------
C
C     IMPROVE APPROXIMATION FOR EIG USING BISECTION OR SECANT METHOD.
C     SUBSTITUTE 'BALLPARK' ESTIMATE IF APPROXIMATION GROWS TOO LARGE.
C
          DEDW = (EU-EL)/ (WU-WL)
          FOLD = 0.0D0
          IF (INTAB.EQ.1) BALLPK = (PIN/ (A-B))**2
          IF (INTAB.EQ.1 .AND. PR) WRITE (T21,FMT=*) ' BALLPK = ',BALLPK
C             ---------------------------------
          LOGIC = .TRUE.
   90     CONTINUE
C      NOW TRY TO REFINE THE ESTIMATE EEE:
C       USE A SECANT METHOD.
          IF (LOGIC) THEN
              LOGIC = (WL.LT.PIN-ONE .OR. WU.GT.PIN+ONE)
              EEE = EL + DEDW* (PIN-WL)
              FNEW = MIN(PIN-WL,WU-PIN)
              IF (FNEW.GT.0.4D0*FOLD .OR. FNEW.LE.ONE) EEE = 0.5D0*
     +            (EL+EU)
              IF (INTAB.EQ.1 .AND. ABS(EEE).GT.1.0D3*BALLPK) THEN
                  EEE = BALLPK
                  GO TO 100
              ELSE IF (INTAB.NE.1 .AND. ABS(EEE).GT.1.0D6) THEN
                  EEE = ONE
                  GO TO 100
              ELSE
                  FOLD = FNEW
                  CALL ESTPAC(.FALSE.,MF,ML,EEE,SUM0,QS,WS,DS,DELT,PS,
     +                        PSS,TAU,JJL,JJR,SUM,U)
                  IF (SUM.LT.PIN) THEN
                      EL = EEE
                      WL = SUM
                      ULO = U
                  ELSE
                      EU = EEE
                      WU = SUM
                      UUP = U
                  END IF
                  DEDW = (EU-EL)/ (WU-WL)
                  GO TO 90
              END IF
          END IF
C             -----------------------------
      END IF
  100 CONTINUE
      RETURN
      END
C
      SUBROUTINE PRELIM(MF,ML,EEE,BALLPK,XS,QS,WS,DS,DELT,PS,PSS,TAU,
     +                  JJL,JJR,LIMUP,ELIMUP,EL,WL,DEDW,GUESS,EIG,AAA,
     +                  AA,BBB,BB,NCA,NCB,EIGPI,ALFA,BETA,DTHDEA,DTHDEB,
     +                  IMID,TMID)

C  THIS PROGRAM, USING THE FIRST ESTIMATE FOR THE EIGENVALUE
C    OBTAINED BY ESTIM,
C    SETS THE INITIAL INTERVAL (POSSIBLY A SUBINTERVAL OF (A,B))
C    TO BE USED FOR THE PROBLEM.  IT IS CALLED BY SLEIGN.
C    IT DETERMINES ALFA, BETA (THE
C    BOUNDARY VALUES OF THE PRUEFER ANGLE THETA), TRIES TO IMPROVE
C    THE ESTIMATE OF THE EIGENVALUE (TAKING ALFA, BETA INTO ACCOUNT)
C    AND SETS THE MIDPOINT, TMID, IN (-1,1) TO BE USED IN THE
C    COMPUTATIONS.
C
C    IT ALSO CHECKS TO SEE IF AN ENDPOINT IS LP, AND IF SO, WHETHER
C    THE POINT SPECTRUM MIGHT BE BOUNDED ABOVE.

C  INPUT QUANTITIES: MF,ML,EEE,BALLPK,XS,QS,WS,DS,DELT,PS,PSS,
C                    TAU,JJL,JJR,LIMUP,ELIMUP,EL,WL,DEDW,GUESS,
C                    EIG,NCA,NCB,EIGPI
C  OUTPUT QUANTITIES: EEE,PSS,JJL,JJR,TEE,JAY,ZEE,LIMUP,ELIMUP,
C                     AAA,AA,BBB,BB,ALFA,BETA,DTHDEA,DTHDEB,
C                     IMID,TMID
C     .. Scalar Arguments ..
      REAL AA,AAA,ALFA,BALLPK,BB,BBB,BETA,DEDW,DTHDEA,
     +                 DTHDEB,EEE,EIG,EIGPI,EL,ELIMUP,GUESS,TAU,TMID,WL
      INTEGER IMID,JJL,JJR,MF,ML,NCA,NCB
      LOGICAL LIMUP
C     ..
C     .. Array Arguments ..
      REAL DELT(100),DS(100),PS(100),PSS(100),QS(100),
     +                 WS(100),XS(100)
C     ..
C     .. Scalars in Common ..
      REAL A,A1,A2,B,B1,B2,HPI,P0ATA,P0ATB,PI,QFATA,QFATB,
     +                 TWOPI
      INTEGER INTAB,MFS,MLS,T21
      LOGICAL PR
C     ..
C     .. Arrays in Common ..
      REAL TEE(100),ZEE(100)
      INTEGER JAY(100)
C     ..
C     .. Local Scalars ..
      REAL BOUNDA,BOUNDB,DERIVL,DERIVR,ELIMA,ELIMB,ONE,PIN,
     +                 SUM,SUM0,THOUS,U,TMP
      INTEGER I,IPA,IPB,JFLAG,KFLAG
      LOGICAL GESS0,LCOA,LCOB,LIMA,LIMB,SINGA,SINGB
C     ..
C     .. External Functions ..
      REAL TFROMI
      EXTERNAL TFROMI
C     ..
C     .. External Subroutines ..
      EXTERNAL ALFBET,ENDPT,ESTPAC,SETMID
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS,MAX,MIN,SIGN
C     ..
C     .. Common blocks ..
      COMMON /BCDATA/A1,A2,P0ATA,QFATA,B1,B2,P0ATB,QFATB
      COMMON /DATADT/A,B,INTAB
      COMMON /LP/MFS,MLS
      COMMON /PIE/PI,TWOPI,HPI
      COMMON /PRIN/PR,T21
      COMMON /TEEZ/TEE
      COMMON /ZEEZ/JAY,ZEE
C     ..
      THOUS = 1000.0D0
      PIN = EIGPI + PI
      SINGA = NCA .GE. 3
      SINGB = NCB .GE. 3
      LCOA = NCA .EQ. 4
      LCOB = NCB .EQ. 4
      GESS0 = ABS(GUESS) .LE. (1D-7)
      ONE = 1.0D0
      SUM0 = 0.0D0
C  LIMUP = .TRUE. MEANS IT APPEARS THE EIGENVALUES ARE BOUNDED ABOVE.
C  IN THIS CASE, ELIMUP IS THE ESTIMATED UPPER BOUND -- OR, RATHER,
C  THE LOWER BOUND OF THE CONTINUOUS SPECTRUM.
      IF (LIMUP .AND. EEE.GE.ELIMUP) EEE = ELIMUP - 1.0D0

C        SET-INITIAL-INTERVAL-AND-MATCHPOINT
      IF (.NOT.GESS0) THEN
          EEE = EIG
          CALL ESTPAC(.FALSE.,MF,ML,EEE,SUM0,QS,WS,DS,DELT,PS,PSS,TAU,
     +                JJL,JJR,SUM,U)
      END IF
C
C     CHOOSE INITIAL INTERVAL AS LARGE AS POSSIBLE THAT AVOIDS OVERFLOW.
C     JJL AND JJR ARE BOUNDARY INDICES FOR NONNEGATIVITY OF EIG*W-Q.
C
      IF (.NOT.SINGA) THEN
          AA = -ONE
      ELSE IF (LCOA) THEN
          AA = -0.9999D0
          AAA = AA
      ELSE
          AA = TFROMI(JJL)
          AAA = MIN(AAA,AA)
          AA = MAX(AA,AAA)
          TMP = -0.99D0
          AA = MIN(AA,TMP)
          AA = MAX(AA,AAA)
      END IF
      IF (.NOT.SINGB) THEN
          BB = ONE
      ELSE IF (LCOB) THEN
          BB = 0.9999D0
          BBB = BB
      ELSE
          BB = TFROMI(JJR)
          BBB = MAX(BBB,BB)
          BB = MIN(BB,BBB)
          TMP = 0.99D0
          BB = MAX(BB,TMP)
          BB = MIN(BB,BBB)
      END IF
C
C     DETERMINE BOUNDARY VALUES ALFA AND BETA FOR THETA AT A AND B.
C
      CALL ALFBET(A,INTAB,AA,A1,A2,EEE,P0ATA,QFATA,SINGA,ALFA,KFLAG,
     +            DERIVL)
      CALL ALFBET(B,INTAB,BB,B1,B2,EEE,P0ATB,QFATB,SINGB,BETA,JFLAG,
     +            DERIVR)
      IF (SINGB) BETA = PI - BETA
C
C     TAKE BOUNDARY CONDITIONS INTO ACCOUNT IN ESTIMATION OF EIG.
C
      PIN = EIGPI + BETA - ALFA
      IF (LCOA) PIN = PIN + ALFA
      IF (LCOB) PIN = PIN + PI - BETA

      IF (GESS0) THEN
          EEE = EL + DEDW* (PIN-WL)
          IF (.NOT. (LCOA.OR.LCOB) .AND.
     +        ABS(EEE).GT.THOUS) EEE = SIGN(THOUS,EEE)
          IF (INTAB.EQ.1 .AND. ABS(EEE).GT.1.0D3*BALLPK) EEE = BALLPK
      END IF

      CALL ESTPAC(.TRUE.,MF,ML,EEE,SUM0,QS,WS,DS,DELT,PS,PSS,TAU,JJL,
     +            JJR,SUM,U)
C
C     RESET BOUNDARY VALUES ALFA AND BETA, WHICH MAY DEPEND UPON
C     THE CURRENT VALUE FOR THE EIGENPARAMETER.
C
      CALL ALFBET(A,INTAB,AA,A1,A2,EEE,P0ATA,QFATA,SINGA,ALFA,KFLAG,
     +            DERIVL)
      CALL ALFBET(B,INTAB,BB,B1,B2,EEE,P0ATB,QFATB,SINGB,BETA,JFLAG,
     +            DERIVR)
      IF (SINGB) BETA = PI - BETA
      DTHDEA = DERIVL
      DTHDEB = -DERIVR
      IF (PR) WRITE (T21,FMT='(A,E22.14,A,E22.14)') '  ALFA=',ALFA,
     +    '   BETA=',BETA
C
C     CHOOSE INITIAL MATCHING POINT TMID .
C
      IMID = 50
      TMID = 0.5D0* (AA+BB) + 0.031D0*PI
      IF (PR) WRITE (T21,FMT='(A,E15.7,A,F11.8,A,E15.7)') ' ESTIM=',EEE,
     +    '  TMID=',TMID
      IF (PR) WRITE (T21,FMT='(A,F11.8,A,F11.8,A,F11.8,A,F11.8)')
     +    ' AAA=',AAA,'  AA=',AA,'  BB=',BB,'  BBB=',BBB

C        END (SET-INITIAL-INTERVAL-AND-MATCHPOINT)
C
C  IF EITHER ENDPOINT IS LP, USE SUBROUTINE ENDPT TO SEE IF THAT
C  ENDPOINT REQUIRES THE EIGENVALUES TO BE BOUNDED ABOVE, AND IF SO,
C  WHAT IS AN ESTIMATED UPPER BOUND.
      LIMA = .FALSE.
      IF (NCA.EQ.5 .OR. NCA.EQ.6) THEN
          CALL ENDPT(MF,ML,XS,PS,QS,WS,-ONE,IPA,NCA,BOUNDA)
          IF (PR) WRITE (T21,FMT=*) ' IPA,NCA,BOUNDA = ',IPA,NCA,BOUNDA
          IF (BOUNDA.LT. (10000.0D0)) THEN
              LIMA = .TRUE.
              ELIMA = BOUNDA
          END IF
      END IF
      LIMB = .FALSE.
      IF (NCB.EQ.5 .OR. NCB.EQ.6) THEN
          CALL ENDPT(MF,ML,XS,PS,QS,WS,ONE,IPB,NCB,BOUNDB)
          IF (PR) WRITE (T21,FMT=*) ' IPB,NCB,BOUNDB = ',IPB,NCB,BOUNDB
          IF (BOUNDB.LT. (10000.0D0)) THEN
              LIMB = .TRUE.
              ELIMB = BOUNDB
          END IF
      END IF
C
C  IF BOTH ENDPOINTS CAUSE THE EIGENVALUES TO BE BOUNDED ABOVE,
C  SET ELIMUP TO BE THE LOWEST UPPER BOUND.
C
      IF (LIMA .OR. LIMB) THEN
          LIMUP = .TRUE.

          IF (.NOT.LIMB) THEN
              ELIMUP = ELIMA
          ELSE IF (.NOT.LIMA) THEN
              ELIMUP = ELIMB
          ELSE
              ELIMUP = MIN(ELIMA,ELIMB)
          END IF

          IF (PR) WRITE (T21,FMT=*)
     +        ' THE CONTINUOUS SPECTRUM HAS A LOWER '
          IF (PR) WRITE (T21,FMT=*) '   BOUND, SIGMA0 = ',ELIMUP
      END IF
C
      IF (EIG.EQ.0.0D0 .AND. LIMUP .AND. EEE.GE.ELIMUP) EEE = ELIMUP -
     +    0.01D0

      CALL SETMID(MF,ML,EEE,QS,WS,IMID,TMID)
C  SET THE VALUES OF T, TEE(*), CORRESPONDING TO THE PLACES
C  WHERE THE Z's, ZEE(*), CHANGE WHEN USING SUBROUTINE GERKZ
C  FOR THE INTEGRATIONS FOR THETA.
C  THE VALUES OF JAY(*) WERE SET IN SUBROUTINE ESTPAC.
      DO 85 I = 1,100
          TEE(I) = ONE
          IF (JAY(I).NE.0) TEE(I) = TFROMI(JAY(I))
          IF (ZEE(I).GT.100.0D0) ZEE(I) = 100.0D0
          IF (JAY(I).NE.0 .AND. PR) WRITE (T21,FMT=*) ' I,T,Z = ',I,
     +        TEE(I),ZEE(I)
   85 CONTINUE
      IF (JAY(3).EQ.0 .AND. ZEE(2).NE.1.0D0) THEN
          TEE(3) = TEE(2)
          ZEE(3) = ZEE(2)
          TEE(2) = 0.0D0
          ZEE(2) = 1.0D0
      ENDIF

      RETURN
      END
C
      SUBROUTINE EIGFCN(EIGPI,A1,A2,B1,B2,AOK,BOK,NCA,NCB,SLFUN,ISLFUN)
C     **********
C  THIS PROGRAM CALCULATES SELECTED EIGENFUNCTION VALUES BY
C    INTEGRATION (OVER T).  THE SELECTED VALUES OF T ARE IN
C    THE ARRAY SLFUN, WHICH ARE REPLACED BY THE CALCULATED
C    VALUES OF THE EIGENFUNCTION.
C    IT IS CALLED FROM SLEIGN.
C  INPUT QUANTITIES: EIGPI,A1,A2,B1,B2,AOK,BOK,NCA,NCB,
C                    SLFUN,ISLFUN
C  OUTPUT QUANTITIES: SLFUN
C
C     N.B.: IN THIS PROGRAM IT IS ASSUMED THAT THE POINTS T
C              IN SLFUN ALL LIE WITHIN THE INTERVAL (AA,BB).
C
C     .. Scalars in Common ..
      REAL AA,BB,DTHDAA,DTHDBB,TMID
      INTEGER MDTHZ,T21
      LOGICAL ADDD,PR
C     ..
C     .. Local Scalars ..
      REAL DTHDAT,DTHDBT,DTHDET,EFF,T,THT,TM
      INTEGER I,IFLAG,J,NC,NMID
      LOGICAL LCOA,LCOB,OK
C     ..
C     .. Local Arrays ..
      REAL ERL(3),ERR(3),YL(3),YR(3)
C     ..
C     .. External Subroutines ..
      EXTERNAL INTEG
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC EXP,SIN
C     ..
C     .. Common blocks ..
      COMMON /PRIN/PR,T21
      COMMON /TDATA/AA,TMID,BB,DTHDAA,DTHDBB,MDTHZ,ADDD
C     ..
C     .. Scalar Arguments ..
      REAL A1,A2,B1,B2,EIGPI
      INTEGER ISLFUN,NCA,NCB
      LOGICAL AOK,BOK
C     ..
C     .. Array Arguments ..
      REAL SLFUN(ISLFUN+9)
C     ..
      NMID = 0
      DO 10 I = 1,ISLFUN
          IF (SLFUN(9+I).LE.TMID) NMID = I
   10 CONTINUE
      IF (NMID.GT.0) THEN
          T = AA
          YL(1) = SLFUN(3)
          YL(2) = 0.0D0
          YL(3) = 0.0D0
          LCOA = NCA .EQ. 4
          OK = AOK
          NC = NCA
          EFF = 0.0D0
          DO 20 J = 1,NMID
              TM = SLFUN(J+9)
              IF (TM.LT.AA .OR. TM.GT.BB) THEN
                  IF (PR) WRITE (*,FMT=*) ' T.LT.AA .OR. T.GT.BB '
                  STOP
              END IF
              THT = YL(1)
              DTHDAT = DTHDAA*EXP(-2.0D0*EFF)
              DTHDET = YL(2)
              IF (TM.GT.AA) THEN
                  CALL INTEG(T,THT,DTHDAT,DTHDET,TM,A1,A2,SLFUN(8),YL,
     +                       ERL,OK,NC,IFLAG)
                  IF (LCOA) THEN
                      EFF = YL(3)
                  ELSE
                      NC = 1
                      OK = .TRUE.
                      EFF = EFF + YL(3)
                  END IF
              END IF
              SLFUN(J+9) = SIN(YL(1))*EXP(EFF+SLFUN(4))
              T = TM
              IF (T.GT.-1.0D0) NC = 1
              IF (T.LT.-0.9D0 .AND. LCOA) THEN
                  NC = 4
                  T = AA
                  YL(1) = SLFUN(3)
                  YL(2) = 0.0D0
                  YL(3) = 0.0D0
              END IF
   20     CONTINUE
      END IF
      IF (NMID.LT.ISLFUN) THEN
          T = BB
          YR(1) = SLFUN(6)
          YR(2) = 0.0D0
          YR(3) = 0.0D0
          LCOB = NCB .EQ. 4
          OK = BOK
          NC = NCB
          EFF = 0.0D0
          DO 30 J = ISLFUN,NMID + 1,-1
              TM = SLFUN(J+9)
              IF (TM.LT.AA .OR. TM.GT.BB) THEN
                  IF (PR) WRITE (*,FMT=*) ' T.LT.AA .OR. T.GT.BB '
                  STOP
              END IF
              THT = YR(1)
              DTHDBT = DTHDBB*EXP(-2.0D0*EFF)
              DTHDET = YR(2)
              IF (TM.LT.BB) THEN
                  CALL INTEG(T,THT,DTHDBT,DTHDET,TM,B1,B2,SLFUN(8),YR,
     +                       ERR,OK,NC,IFLAG)
                  IF (LCOB) THEN
                      EFF = YR(3)
                  ELSE
                      OK = .TRUE.
                      NC = 1
                      EFF = EFF + YR(3)
                  END IF
              END IF
              SLFUN(J+9) = SIN(YR(1)+EIGPI)*EXP(EFF+SLFUN(7))
              IF (ADDD) SLFUN(J+9) = -SLFUN(J+9)
              T = TM
              IF (T.LT.1.0D0) NC = 1
              IF (T.GT.0.9D0 .AND. LCOB) THEN
                  NC = 4
                  T = BB
                  YR(1) = SLFUN(6)
                  YR(2) = 0.0D0
                  YR(3) = 0.0D0
              END IF
   30     CONTINUE
      END IF
      RETURN
      END
C
      SUBROUTINE EFDATA(ALFA,DTHDEA,A1,A2,BETA,DTHDEB,B1,B2,EIG,EIGPI,
     +                  EPS,AOK,NCA,BOK,NCB,RAY,SLFUN)
C
C  THIS PROGRAM IS USED ONLY AFTER THE WANTED EIGENVALUE HAS BEEN
C    SUCCESSFULLY COMPUTED. IT IS CALLED FROM SLEIGN.
C    IT CONVERTS THE T-DATA AA,TMID,BB TO CORRESPONDING X-DATA,
C    FILLS 7 OF THE FIRST 9 LOCATIONS OF SLFUN,
C    COMPUTES VALUES OF LOG(RHO(A)), LOG(RHO(B)) SUCH THAT THE
C    CORRESPONDING SOLUTION Y(X) = RHO(X)*SIN(THETA(X)) IS
C    CONTINUOUS AT XMID, AND HAS L2-NORM EQUAL TO 1.0 .
C    THESE VALUES ARE PLACED IN SLFUN(4) AND SLFUN(7).
C  THE COMPUTED VALUE OF EIG IS FINALLY CHECKED ONE LAST TIME
C    USING THE JUMPS IN Y AND PY' AT XMID IN A FORM OF
C    RALEIGH QUOTIENT CORRECTION.  (HERE CALLED "RAY".)
C
C  INPUT QUANTITIES: ALFA,DTHDEA,A1,A2,BETA,DTHDEB,B1,B2,EIG,
C                    EIGPI,EPS,AOK,BOK,NCA,NCB,SLFUN
C  OUTPUT QUANTITIES: XAA,XMID,XBB,EIGSAV,ISAVE,AA,BB,RAY,ADDD
C                     SLFUN
C
C     .. Scalar Arguments ..
      REAL A1,A2,ALFA,B1,B2,BETA,DTHDEA,DTHDEB,EIG,EIGPI,
     +                 EPS,RAY
      INTEGER NCA,NCB
      LOGICAL AOK,BOK
C     ..
C     .. Array Arguments ..
      REAL SLFUN(9)
C     ..
C     .. Scalars in Common ..
      REAL AA,BB,DTHDAA,DTHDBB,EIGSAV,TMID,TSAVEL,TSAVER,Z
      INTEGER IND,ISAVE,MDTHZ,T21
      LOGICAL ADDD,PR
C     ..
C     .. Local Scalars ..
      REAL AAF,BBF,CL,CR,DEN,DUM,E,ONE,PSIL,PSIPL,PSIPR,
     +                 PSIR,SL,SQL,SQR,SR,THDAAX,THDBBX,TMP,UL,UR,XAA,
     +                 XBB,XMID
      INTEGER JFLAG
C     ..
C     .. Local Arrays ..
      REAL ERL(3),ERR(3),YL(3),YR(3)
C     ..
C     .. External Subroutines ..
      EXTERNAL AABB,DXDT,INTEG
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS,COS,EXP,LOG,SIN,SQRT
C     ..
C     .. Common blocks ..
      COMMON /DATAF/EIGSAV,IND
      COMMON /PRIN/PR,T21
      COMMON /TDATA/AA,TMID,BB,DTHDAA,DTHDBB,MDTHZ,ADDD
      COMMON /TSAVE/TSAVEL,TSAVER,ISAVE
      COMMON /Z1/Z
C     ..
      ONE = 1.0D0
      Z = 1.0D0
   50 CONTINUE
      CALL DXDT(TMID,TMP,XMID)
      CALL DXDT(AA,TMP,XAA)
      SLFUN(1) = XMID
      SLFUN(2) = XAA
      SLFUN(3) = ALFA
      SLFUN(6) = BETA + EIGPI
      SLFUN(8) = EPS
      SLFUN(9) = Z
C
C     INTEGRATE FROM BOTH ENDPOINTS TO THE MIDPOINT:
C
      EIGSAV = EIG
      THDAAX = 0.0D0
      YL(1) = 0.0D0
      YL(2) = 0.0D0
      YL(3) = 0.0D0
      ISAVE = 0
      CALL INTEG(AA,ALFA,THDAAX,DTHDEA,TMID,A1,A2,EPS,YL,ERL,AOK,NCA,
     +           JFLAG)
      IF (JFLAG.EQ.5) THEN
          AAF = AA
          CALL AABB(AA,-ONE)
          IF (PR) WRITE (T21,FMT=*) ' AA MOVED FROM ',AAF,'IN TO',AA
          GO TO 50
      END IF
  100 CONTINUE
      CALL DXDT(BB,TMP,XBB)
      SLFUN(5) = XBB
      THDBBX = 0.0D0
      YR(1) = 0.0D0
      YR(2) = 0.0D0
      YR(3) = 0.0D0
      ISAVE = 0
      CALL INTEG(BB,BETA,THDBBX,DTHDEB,TMID,B1,B2,EPS,YR,ERR,BOK,NCB,
     +           JFLAG)
      IF (JFLAG.EQ.5) THEN
          BBF = BB
          CALL AABB(BB,-ONE)
          IF (PR) WRITE (T21,FMT=*) ' BB MOVED FROM ',BBF,'IN TO',BB
          GO TO 100
      END IF
      YR(1) = YR(1) + EIGPI
      SL = SIN(YL(1))
      SR = SIN(YR(1))
      CL = COS(YL(1))
      CR = COS(YR(1))
C  UL AND UR ARE THE VALUES OF THE INTEGRAL OF Y**2*W FROM THE
C  ENDPOINTS A AND B, RESPECTIVELY.
      UL = (YL(2)-DTHDEA*EXP(-2.0D0*YL(3)))
      UR = (YR(2)-DTHDEB*EXP(-2.0D0*YR(3)))
      DEN = 0.5D0*LOG(UL*SR*SR-UR*SL*SL)
      DUM = 0.5D0*LOG(UL-UR)
C     COMPUTE SLFUN(4), SLFUN(7) TOWARDS NORMALIZING THE EIGENFUNCTION.
      SLFUN(4) = -YL(3) - DUM
      SLFUN(7) = -YR(3) - DUM

C     DO (CHECK-MATCHING-VALUES-OF-EIGENFUNCTION)
C     PERFORM FINAL CHECK ON EIG.
C
      DEN = UL*SR*SR - UR*SL*SL
      E = ABS(SR)/SQRT(DEN)
      PSIL = E*SL
      PSIPL = E*CL
      SQL = E*E*UL
      E = ABS(SL)/SQRT(DEN)
      PSIR = E*SR
      PSIPR = E*CR
      SQR = E*E*UR
      ADDD = PSIL*PSIR .LT. 0.0D0 .AND. PSIPL*PSIPR .LT. 0.0D0
      RAY = EIG + (PSIL*PSIPL-PSIR*PSIPR)/ (SQL-SQR)
      IF (PR) WRITE (T21,FMT=*) ' RAY,ADDD = ',RAY,ADDD
C           END (CHECK-MATCHING-VALUES-OF-EIGENFUNCTION)
      RETURN
      END
C
      SUBROUTINE SETMID(MF,ML,EIG,QS,WS,IMID,TMID)
C     **********
C
C     THIS PROGRAM SELECTS A POINT IN (-1,1) AS TMID.
C     IT IS CALLED BY RESET AND PRELIM.
C     IT TESTS THE INTERVAL SAMPLE POINTS IN THE ORDER
C     50,51,49,52,48,...,ETC. FOR THE FIRST ONE WHERE THE EXPRESSION
C     (LAMBDA*W-Q) IS POSITIVE.  THIS T-POINT IS DESIGNATED TMID.
C  THEORETICALLY IT DOESN'T MATTER WHICH POINT T IN (-1,1) IS USED
C    FOR TMID, BUT IN PRACTICE IT CAN MAKE A GREAT DEAL OF
C    DIFFERENCE. THIS SCHEME IS JUST A RULE OF THUMB WHICH SEEMS
C    TO WORK FAIRLY WELL IN AVOIDING BAD CHOICES FOR TMID.
C
C   INPUT QUANTITIES: MF, ML, EIG, QS, WS
C   OUTPUT QUANTITIES: IMID, TMID
C     **********
C     .. Local Scalars ..
      REAL S
      INTEGER I,J
C     ..
C     .. External Functions ..
      REAL TFROMI
      EXTERNAL TFROMI
C     ..
C     .. Scalar Arguments ..
      REAL EIG,TMID
      INTEGER IMID,MF,ML
C     ..
C     .. Array Arguments ..
      REAL QS(*),WS(*)
C     ..
C     .. Scalars in Common ..
      INTEGER T21
      LOGICAL PR
C     ..
C     .. Common blocks ..
      COMMON /PRIN/PR,T21
C     ..
      S = -1.0D0
      DO 10 J = 1,100
          I = 50 + S* (J/2)
          S = -S
          IF (I.LT.MF .OR. I.GT.ML) GO TO 20
          IF (EIG*WS(I)-QS(I).GT.0.0D0) THEN
              IMID = I
              TMID = TFROMI(IMID)
              GO TO 20
          END IF
   10 CONTINUE
   20 CONTINUE
      IF (PR) WRITE (T21,FMT=*) ' NEW TMID = ',TMID
      RETURN
      END
C
      SUBROUTINE EXTR(MM,F,GQ,GW,XS,PS,QS,WS)

C  THIS SUBROUTINE IS CALLED FROM SUBROUTINE ENDPT ONLY,
C    FOR THE PURPOSE OF DETERMINING THE ENDPOINT VALUES
C    OF F (OR Q, OR W) IN THE INDICIAL EQUATION AT A FINITE
C    LP ENDPOINT.
C    IT PUTS VALUES OF F (OR Q OR W), CORRESPONDING TO
C    EQUALLY SPACED VALUES OF X, IN
C    AN ARRAY DF, AND EXTRAPOLATES THE VALUES OF F (OR Q, OR W)
C    TO THE ENDPOINT, XEND, OF THE INTERVAL.
C
C  INPUT QUANTITIES: MM,XS,PS,QS,WS,A,B,EPSMIN
C  OUTPUT QUANTITIES: F,GQ,GW

C     .. Scalar Arguments ..
      REAL F,GQ,GW
      INTEGER MM
C     ..
C     .. Array Arguments ..
      REAL PS(100),QS(100),WS(100),XS(100)
C     ..
C     .. Scalars in Common ..
      REAL A,B,EPSMIN
      INTEGER INTAB
C     ..
C     .. Local Scalars ..
      REAL TI,TIE,TMP,XEND,XI,XIE,EPST,ERR
      INTEGER EP,I,J
C     ..
C     .. Local Arrays ..
      REAL DF(6),DQ(6),DW(6),XX(6)
C     ..
C     .. External Functions ..
      REAL P,TFROMI
      EXTERNAL P,TFROMI
C     ..
C     .. External Subroutines ..
      EXTERNAL DXDT,INTPOL
C     ..
C     .. Common blocks ..
      COMMON /DATADT/A,B,INTAB
      COMMON /RNDOFF/EPSMIN
C     ..
      IF (MM.LT.50) THEN
          EP = 1
          XEND = A
      ELSE
          EP = -1
          XEND = B
      END IF

C  HERE WE WANT TO FIND THE LIMIT OF
C   (X-A)*P'(X)/P(X) AT X=A
C  AND OF
C    (X-A)**2*Q(X)/P(X)
C  AND OF
C    (X-A)**2*W(X)/P(X)
C  OR THE SAME SORT OF THINGS AT X=B.

      DO 10 J = 1,5
          I = MM + (J-1)*EP
          TI = TFROMI(I)
          TIE = TI + 2.0D0*EPSMIN*EP
          CALL DXDT(TI,TMP,XI)
          CALL DXDT(TIE,TMP,XIE)
          DF(J) = (XS(I)-XEND)* (P(XIE)-PS(I))/ ((XIE-XS(I))*PS(I))
          DQ(J) = (XS(I)-XEND)**2*QS(I)
          DW(J) = (XS(I)-XEND)**2*WS(I)
          XX(J) = XS(I)
   10 CONTINUE
           EPST = 0.00001D0
           CALL INTPOL(5,XX,DF,XEND,EPST,3,F,ERR)
           CALL INTPOL(5,XX,DQ,XEND,EPST,3,GQ,ERR)
           CALL INTPOL(5,XX,DW,XEND,EPST,3,GW,ERR)
      IF(ABS(F).LE.ERR) F = 0.0D0
      IF(ABS(GQ).LE.ERR) GQ = 0.0D0
      IF(GW.LE.ERR) GW = 0.0D0 
      RETURN
      END
C
      SUBROUTINE ENDPT(MF,ML,XS,PS,QS,WS,TEND,IPM,NC,BOUND)
C  THIS PROGRAM IS CALLED FROM PRELIM.
C  MAINLY IT COMPUTES THE QUANTITIES IN /COMMON/ALBE/ FOR
C    USE WHEN A FINITE ENDPOINT IS LP.
C    IF THE ENDPOINT IS SUCH THAT THE QUANTITIES INVOLVED
C    APPEAR TO NOT HAVE LIMITING VALUES, (AS HAPPENS IN THE
C    PROBLEM
C         p(x) = 1, q(x) = cos(x)**2, w(x) = 1  on (0,+Inf)
C    FOR EXAMPLE),  THEN SLEIGN2 CANNOT CONTINUE.
C
C    AT THE SAME TIME, IT TRIES TO DETERMINE WHETHER OR NOT
C    THE EIGENVALUES HAVE A FINITE UPPER BOUND.
C
C  INPUT QUANTITIES: MF,ML,XS,PS,QS,WS,TEND,INTAB,EPSMIN
C  OUTPUT QUANTITIES: LPWA,LPQA,LPWB,LPQB,NC,BOUND
C
C  RECALL THAT THE STORED ARRAYS QS AND WS ARE REALLY THE
C    VALUES OF Q/P AND W/P.
C  N.B.: INDICIAL EQUATION IS (AT A REGULAR SINGULAR POINT):
C            S*S + (F-1)*S + G = 0
C      WHERE  F IS FA OR FB
C      AND    G IS LAMBDA*GWA - GQA
C               OR LAMBDA*GWB - GQB
C  (HERE, THE INDEPENDENT VARIABLE IS X IN (A,B)

C     .. Scalar Arguments ..
      REAL BOUND,TEND
      INTEGER IPM,MF,ML,NC
C     ..
C     .. Array Arguments ..
      REAL PS(100),QS(100),WS(100),XS(100)
C     ..
C     .. Scalars in Common ..
      REAL A,B,EPSMIN,LPQA,LPQB,LPWA,LPWB
      INTEGER INTAB,T21
      LOGICAL PR
C     ..
C     .. Local Scalars ..
      REAL APQ1,APQ2,APW1,APW2,PQ1,PQ2,PW1,PW2,TMP,
     +     FA,FB,GQA,GQB,GWA,GWB
      INTEGER I,IQ,IW,MM,MM6
      LOGICAL LOGIC
C     ..
C     .. External Subroutines ..
      EXTERNAL EXTR
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS,MIN
C     ..
C     .. Common blocks ..
c     COMMON /ALBE/LPWA,LPQA,LPWB,LPQB
      COMMON /ALBE/LPWA,LPQA,FA,GWA,GQA,LPWB,LPQB,FB,GWB,GQB
      COMMON /DATADT/A,B,INTAB
      COMMON /PRIN/PR,T21
      COMMON /RNDOFF/EPSMIN
C     ..
      IF (PR) WRITE (T21,FMT=*) ' IN ENDPT, NC = ',NC
      IF (TEND.LT.0.0D0) THEN
          MM = MF
      ELSE
          MM = ML
      END IF
      BOUND = 1.D+19
      IQ = 0
      IW = 0
      IF (MM.LT.50) THEN
          MM6 = MM + 6
          DO 10 I = MM,MM6
              PW1 = PS(I+1)**2*WS(I+1) - PS(I)**2*WS(I)
              PW2 = PS(I+2)**2*WS(I+2) - PS(I+1)**2*WS(I+1)
              PQ1 = PS(I+1)**2*QS(I+1) - PS(I)**2*QS(I)
              PQ2 = PS(I+2)**2*QS(I+2) - PS(I+1)**2*QS(I+1)
              APW1 = ABS(PW1)
              APW2 = ABS(PW2)
              APQ1 = ABS(PQ1)
              APQ2 = ABS(PQ2)
              IF ((APW1.LE.EPSMIN.AND.APW2.LE.EPSMIN) .OR.
     +            PW1*PW2.GT.0.D0) IW = IW + 1
              IF ((APQ1.LE.EPSMIN.AND.APQ2.LE.EPSMIN) .OR.
     +            PQ1*PQ2.GT.0.D0) IQ = IQ + 1
   10     CONTINUE
          IPM = MIN(IW,IQ)
          IF (IPM.GE.5) THEN
              LPWA = PS(MM)**2*WS(MM)
              LPQA = PS(MM)**2*QS(MM)
              IF (INTAB.EQ.1 .OR. INTAB.EQ.2) THEN
                  CALL EXTR(MM,FA,GQA,GWA,XS,PS,QS,WS)
                  IF (PR) WRITE (T21,FMT=*) ' FA,GWA,GQA = ',FA,GWA,GQA
                  LOGIC = (ABS(FA)+ABS(GWA)+ABS(GQA)) .LT. 1000.0D0
                  IF (LOGIC) THEN

C  INDICIAL EQUATION IS:  S*S + (FA-1)*S + (LAMBDA*GWA-GQA)=0

                      IF (GWA.GT.EPSMIN) THEN
C               FOR REAL ROOTS, REQUIRE LAMBDA .LE. BOUND, WHERE:
                          TMP = (0.25D0*(FA-1.0D0)**2+GQA)/GWA
                          BOUND = MIN(BOUND,TMP)
                      END IF
                  ELSE
C  THIS CASE MEANS THAT THE ENDPOINT IS PROBABLY AN "IRREGULAR"
C    POINT.
                      NC = 6
                  END IF
              ELSE
C  THIS CASE MEANS THAT THE ENDPOINT A IS AT -INF.  SO
                  IF (LPWA.GT.EPSMIN) BOUND = MIN(BOUND,LPQA/LPWA)
C  CHOOSE NC=7 ONLY IF THE USER HAS NOT ALREADY CHOSEN NC=6.
                  IF (NC.NE.6) NC = 7
C    TEMPORARILY SET (NC = 7 IS NOT BEING USED YET)
                  NC = 6
              END IF
          ELSE
C           ARRIVING HERE MEANS THAT IPM.LT.7, WHICH MEANS THAT
C           P*W AND P*Q DO NOT APPEAR TO BE MONOTONELY TENDING
C           TO A LIMIT (OR +INF OR -INF).  I SUPPOSE THIS MAY
C           MEAN THAT THERE ARE BANDS & GAPS.  I DON'T KNOW.
              IF (PR) WRITE (T21,FMT=*)
     +            ' P*W & P*Q BEHAVE BADLY NEAR THE END A. '
              NC = 8
          END IF
      ELSE
          MM6 = MM - 6
          DO 20 I = MM,MM6,-1
              PW1 = PS(I-1)**2*WS(I-1) - PS(I)**2*WS(I)
              PW2 = PS(I-2)**2*WS(I-2) - PS(I-1)**2*WS(I-1)
              PQ1 = PS(I-1)**2*QS(I-1) - PS(I)**2*QS(I)
              PQ2 = PS(I-2)**2*QS(I-2) - PS(I-1)**2*QS(I-1)
              APW1 = ABS(PW1)
              APW2 = ABS(PW2)
              APQ1 = ABS(PQ1)
              APQ2 = ABS(PQ2)
              IF ((APW1.LE.EPSMIN.AND.APW2.LE.EPSMIN) .OR.
     +            PW1*PW2.GT.0.D0) IW = IW + 1
              IF ((APQ1.LE.EPSMIN.AND.APQ2.LE.EPSMIN) .OR.
     +            PQ1*PQ2.GT.0.D0) IQ = IQ + 1
   20     CONTINUE
          IPM = MIN(IQ,IW)
          IF (IPM.GE.5) THEN
              LPWB = PS(MM)**2*WS(MM)
              LPQB = PS(MM)**2*QS(MM)
              IF (INTAB.EQ.1 .OR. INTAB.EQ.3) THEN
                  CALL EXTR(MM,FB,GQB,GWB,XS,PS,QS,WS)
                  IF (PR) WRITE (T21,FMT=*) ' FB,GWB,GQB = ',FB,GWB,GQB
                  LOGIC = (ABS(FB)+ABS(GWB)+ABS(GQB)) .LT. 1000.0D0
                  IF (LOGIC) THEN
C  INDICIAL EQUATION IS:  S*S + (FB-1)*S + (LAMBDA*GWB-GQB)=0
                      IF (GWB.GT.EPSMIN) THEN
C    FOR REAL ROOTS, REQUIRE LAMBDA .LE. BOUND, WHERE:
                          TMP = (0.25D0*(FB-1.0D0)**2+GQB)/GWB
                          BOUND = MIN(BOUND, TMP)
                      END IF
                  ELSE
C  THIS CASE MEANS THAT THE ENDPOINT IS PROBABLY AN "IRREGULAR"
C    POINT.
                      NC = 6
                  END IF
              ELSE
C  THIS CASE MEANS THAT THE ENDPOINT B IS AT +INF.  SO
                  IF (LPWB.GT.EPSMIN) BOUND = MIN(BOUND,LPQB/LPWB)
C  CHOOSE NC=7 ONLY IF THE USER HAS NOT ALREADY CHOSEN NC=6.
                  IF (NC.NE.6) NC = 7
C  TEMPORARILY SET (NC = 7 IS NOT BEING USED YET)
                  NC = 6
              END IF
          ELSE
C         ARRIVING HERE MEANS THAT IPM.LT.7, WHICH MEANS THAT
C         P*W AND P*Q DO NOT APPEAR TO BE MONOTONELY TENDING
C         TO A LIMIT (OR +INF OR -INF).  I SUPPOSE THIS MAY
C         MEAN THAT THERE ARE BANDS & GAPS.  I DON'T KNOW.
              IF (PR) WRITE (T21,FMT=*)
     +            ' P*W & P*Q BEHAVE BADLY NEAR THE END B. '
              NC = 8
          END IF
      END IF
      RETURN
      END
C
      SUBROUTINE THZ2TH(U,ERU,Z,Y,ERY)
C     **********

C  THIS PROGRAM IS CALLED FROM GERKZ.  THERE THE PRUEFER
C    ANGLE HAS A CONSTANT Z IN ITS DEFINITION, AND IS HERE
C    DENOTED BY THZ.  THE USUAL THETA (EQUIVALENT TO Z = 1)
C    IS DENOTED BY TH.
C  THIS PROGRAM CONVERTS FROM THZ TO TH WHERE
C    TAN(TH)=Y/(P*Y')
C  AND
C    TAN(THZ)=Z*Y/(P*Y'),
C  OR
C    TAN(THZ)=Z*TAN(TH) .
C  SO WE HAVE
C    DTHZ=Z*(COS(THZ)/COS(TH))**2 * DTH  ,
C  OR
C    DTH=(1/Z)*(COS(TH)/COS(THZ))**2 * DTHZ .
C
C  INPUT QUANTITIES: U,ERU,Z
C  OUTPUT QUANTITIES: Y,ERY
C     **********
C     .. Scalars in Common ..
      REAL HPI,PI,TWOPI
C     ..
C     .. Local Scalars ..
      REAL DTH,DTHZ,DUM,FAC,PIK,REMTHZ,TH,THZ
      INTEGER K
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS,ATAN,COS,LOG,SIN,TAN
C     ..
C     .. Common blocks ..
      COMMON /PIE/PI,TWOPI,HPI
C     ..
C     .. Scalar Arguments ..
      REAL Z
C     ..
C     .. Array Arguments ..
      REAL ERU(3),ERY(3),U(3),Y(3)
C     ..
      THZ = U(1)
      DTHZ = U(2)
      K = THZ/PI
      IF (THZ.LT.0.0D0) K = K - 1
      PIK = K*PI
      REMTHZ = THZ - PIK
      IF (4.0D0*REMTHZ.LE.PI) THEN
          TH = ATAN(TAN(REMTHZ)/Z) + PIK
      ELSE IF (4.0D0*REMTHZ.GE.3.0D0*PI) THEN
          TH = ATAN(TAN(REMTHZ)/Z) + PIK + PI
      ELSE
          TH = ATAN(Z*TAN(REMTHZ-HPI)) + PIK + HPI
      END IF
C  THE TWO DIFFERENT APPEARING FORMULAS BELOW FOR FAC ARE
C    IN FACT EQUIVALENT.  WE USE WHICHEVER ONE AVOIDS
C    DIVIDING BY A SMALL NUMBER.
      DUM = ABS(COS(TH))
      IF (DUM.GE.0.5D0) THEN
          FAC = (DUM/COS(THZ))**2/Z
      ELSE
          FAC = Z* (SIN(TH)/SIN(THZ))**2
      END IF
      DTH = FAC*DTHZ
      Y(1) = TH
      Y(2) = DTH
      Y(3) = U(3) + 0.5D0*LOG(Z/FAC)

C  ALSO CONVERT THE ESTIMATED ERRORS IN THE THZ COMPUTATIONS
C  TO CORRESPONDING ERRORS IN TH.
      ERY(1) = FAC*ERU(1)
      ERY(2) = FAC*ERU(2)
      ERY(3) = FAC*ERU(3)
      RETURN
      END
C
      SUBROUTINE TH2THZ(Y,Z,U)
C     **********

C  THIS PROGRAM IS CALLED FROM GERKZ.  THERE THE PRUEFER
C    ANGLE HAS A CONSTANT Z IN ITS DEFINITION, AND IS HERE
C    DENOTED BY THZ.  THE USUAL THETA (EQUIVALENT TO Z = 1)
C    IS DENOTED BY TH.
C  THIS PROGRAM CONVERTS FROM TH TO THZ WHERE
C    TAN(TH)=Y/(P*Y')
C  AND
C    TAN(THZ)=Z*Y/(P*Y'),
C  OR
C    TAN(THZ)=Z*TAN(TH) .
C  SO WE HAVE
C    DTHZ=Z*(COS(THZ)/COS(TH))**2 * DTH  ,
C  OR
C    DTH=(1/Z)*(COS(TH)/COS(THZ))**2 * DTHZ .
C
C  INPUT QUANTITIES: Y,Z
C  OUTPUT QUANTITIES: U
C     **********

C     .. Scalars in Common ..
      REAL HPI,PI,TWOPI
C     ..
C     .. Local Scalars ..
      REAL DTH,DTHZ,DUM,FAC,PIK,REMTH,TH,THZ
      INTEGER K
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS,ATAN,COS,LOG,SIN,TAN
C     ..
C     .. Common blocks ..
      COMMON /PIE/PI,TWOPI,HPI
C     ..
C     .. Scalar Arguments ..
      REAL Z
C     ..
C     .. Array Arguments ..
      REAL U(3),Y(3)
C     ..
      TH = Y(1)
      DTH = Y(2)
      K = TH/PI
      IF (TH.LT.0.0D0) K = K - 1
      PIK = K*PI
      REMTH = TH - PIK
      IF (4.0D0*REMTH.LE.PI) THEN
          THZ = ATAN(Z*TAN(REMTH)) + PIK
      ELSE IF (4.0D0*REMTH.GE.3.0D0*PI) THEN
          THZ = ATAN(Z*TAN(REMTH)) + PIK + PI
      ELSE
          THZ = ATAN(TAN(REMTH-HPI)/Z) + PIK + HPI
      END IF
C  THE TWO DIFFERENT APPEARING FORMULAS BELOW FOR FAC ARE
C    IN FACT EQUIVALENT.  WE USE WHICHEVER ONE AVOIDS
C    DIVIDING BY A SMALL NUMBER.
      DUM = ABS(COS(THZ))
      IF (DUM.GE.0.5D0) THEN
          FAC = Z* (DUM/COS(TH))**2
      ELSE
          FAC = (SIN(THZ)/SIN(TH))**2/Z
      END IF
      DTHZ = FAC*DTH
      U(1) = THZ
      U(2) = DTHZ
      U(3) = Y(3) - 0.5D0*LOG(Z*FAC)
      RETURN
      END
C
      SUBROUTINE PQEXT(F)

C  THIS SUBROUTINE IS CALLED ONLY BY SUBROUTINE WR.	
C    IT IS USED TO EXTRAPOLATE THE VALUES F(I),
C      I=2,5, TO F(1) . IT ASSUMES THE INDEPENDENT
C      VARIABLE VALUES ARE EQUALLY SPACED.
C
C  INPUT QUANTITIES: F(2),F(3),F(4),F(5)
C  OUTPUT QUANTITIES: F(1)

C     .. Array Arguments ..
      REAL F(6)
C     ..
C     .. Local Scalars ..
      REAL D2F1,D2F2,D2F3,D3F2
      INTEGER I
C     ..
C     .. Local Arrays ..
      REAL DF(6)
C     ..
      DO 10 I = 2,4
          DF(I) = F(I+1) - F(I)
   10 CONTINUE
      D2F3 = DF(4) - DF(3)
      D2F2 = DF(3) - DF(2)
      D3F2 = D2F3 - D2F2
      D2F1 = D2F2 + D3F2
      DF(1) = DF(2) + D2F1
      F(1) = F(2) + DF(1)
      RETURN
      END
C
      SUBROUTINE FIT(TH1,TH,TH2)
C     **********
C
C  THIS PROGRAM IS CALLED ONLY FROM SUBROUTINE UVPHI, WHICH
C    IS CALLED BY SUBROUTINE LCO.
C     IT CONVERTS TH INTO AN 'EQUIVALENT' ANGLE BETWEEN
C     TH1 AND TH2.  WE ASSUME TH1.LT.TH2 AND PI.LE.(TH2-TH1).
C
C  INPUT QUANTITIES: TH1,TH2
C  OUTPUT QUANTITIES: TH
C
C     **********
C     .. Scalars in Common ..
      REAL HPI,PI,TWOPI
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC AINT
C     ..
C     .. Common blocks ..
      COMMON /PIE/PI,TWOPI,HPI
C     ..
C     .. Scalar Arguments ..
      REAL TH,TH1,TH2
C     ..
      IF (TH.LT.TH1) TH = TH + AINT((TH1-TH+PI)/PI)*PI
      IF (TH.GT.TH2) TH = TH - AINT((TH-TH2+PI)/PI)*PI
      RETURN
      END
C
      SUBROUTINE F(U,Y,YP)
C     **********
C
C     THIS SUBROUTINE EVALUATES THE DERIVATIVE FUNCTIONS FOR USE WITH
C     INTEGRATOR GERK IN SOLVING THE SYSTEM OF DIFFERENTIAL
C     EQUATIONS FOR THETA, RHO (OR, RATHER, LOG(RHO)), AND
C     DTHDE = D(THETA)/D(LAMBDA), WHERE
C        Y = RHO*SIN(THETA)
C     AND
C       PY' = Z*RHO*COS(THETA) .
C     THE DIFFERENTIAL EQUATIONS ARE:
C       THETA' = (Z/P)*COS(THETA)**2 + (EIG*W - Q)*SIN(THETA)**2/Z
C       LOG(RHO)' = (Z/P - (EIG*W - Q)/Z)*SIN(THETA)*COS(THETA)
C       DTHDE' = -2*(Z/P - (EIG*W - Q)/Z)*DTHDE + W*SIN(THETA)**2/Z
C
C     EXCEPT WHEN CALLED FROM GERKZ, Z IS ALWAYS = 1.0 .
C
C  INPUT QUANTITIES: U,Y,EIG,IND,Z
C  OUTPUT QUANTITIES: YP
C     **********
C     .. Scalars in Common ..
      REAL EIG,Z
      INTEGER IND
C     ..
C     .. Local Scalars ..
      REAL C,C2,DT,QX,S,S2,T,TH,V,WW,WX,X,XP
C     ..
C     .. External Functions ..
      REAL P,Q,W
      EXTERNAL P,Q,W
C     ..
C     .. External Subroutines ..
      EXTERNAL DXDT
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC COS,MOD,SIN
C     ..
C     .. Common blocks ..
      COMMON /DATAF/EIG,IND
      COMMON /Z1/Z
C     ..
C     .. Scalar Arguments ..
      REAL U
C     ..
C     .. Array Arguments ..
      REAL Y(2),YP(3)
C     ..

      IF (MOD(IND,2).EQ.1) THEN
          T = U
          TH = Y(1)
      ELSE
          T = Y(1)
          TH = U
      END IF
      CALL DXDT(T,DT,X)
      XP = Z/P(X)
      QX = Q(X)/Z
      WX = W(X)/Z
      V = EIG*WX - QX
      S = SIN(TH)
      C = COS(TH)
      S2 = S*S
      C2 = C*C
      YP(1) = DT* (XP*C2+V*S2)
      IF (IND.EQ.1) THEN
          WW = (XP-V)*S*C
          YP(2) = DT* (-2.0D0*WW*Y(2)+WX*S2)
          YP(3) = DT*WW
      ELSE IF (IND.EQ.2) THEN
C  IN THIS CASE THE INDEPENDENT AND DEPENDENT VARIABLES
C    (t and dtheta, etc.) HAVE BEEN INTERCHANGED.
          YP(2) = YP(2)/YP(1)
          YP(3) = YP(3)/YP(1)
          YP(1) = 1.0D0/YP(1)
      ELSE IF (IND.EQ.3) THEN
      ELSE
          YP(1) = 1.0D0/YP(1)
      END IF
      RETURN
      END
C
      SUBROUTINE UVPHI(U,PUP,V,PVP,THU,THV,PHI,TH)
C     **********
C
C  THIS PROGRAM IS CALLED BY SUBROUTINE LCO.
C     IT FINDS TH (= THETA) APPROPRIATE TO THU, THV, AND PHI, WHERE
C     THU IS THE PHASE ANGLE FOR U, AND THV IS THE PHASE ANGLE FOR V.
C
C  INPUT QUANTITIES: U,PUP,V,PVP,THU,THV
C  OUTPUT QUANTITIES: TH
C
C     **********
C     .. Scalars in Common ..
      REAL HPI,PI,TWOPI
C     ..
C     .. Local Scalars ..
      REAL C,D,PYP,S,Y
C     ..
C     .. External Subroutines ..
      EXTERNAL FIT
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ATAN2,COS,SIN
C     ..
C     .. Common blocks ..
      COMMON /PIE/PI,TWOPI,HPI
C     ..
C     .. Scalar Arguments ..
      REAL PHI,PUP,PVP,TH,THU,THV,U,V
C     ..
      TH = THU
      IF (PHI.EQ.0.0D0) RETURN
      IF (THV-THU.LT.PI) THEN
          TH = THV
          IF (PHI.EQ.-HPI) RETURN
          TH = THV - PI
          IF (PHI.EQ.HPI) RETURN
      ELSE
          TH = THV - PI
          IF (PHI.EQ.-HPI) RETURN
          TH = THV - TWOPI
          IF (PHI.EQ.HPI) RETURN
      END IF
      C = COS(PHI)
      S = SIN(PHI)
      Y = U*C + V*S
      PYP = PUP*C + PVP*S
      TH = ATAN2(Y,PYP)
      IF (Y.LT.0.0D0) TH = TH + TWOPI
      D = U*PVP - V*PUP
      IF (D*PHI.GT.0.0D0) THEN
          CALL FIT(THU-PI,TH,THU)
      ELSE
          CALL FIT(THU,TH,THU+PI)
      END IF
      RETURN
      END

      SUBROUTINE WR(FG,NC,TSTAR,YSTAR,THEND,DTHDE,TOUT,Y,TT,YY,IFLAG,
     +              ERR,WORK,IWORK)
C     **********
C
C     THIS PROGRAM IS CALLED BY WREG AND LCNO.
C     IT INTEGRATES Y' = F(T,Y) FROM T = +/-1.0 TO THEND
C     WHEN F CANNOT BE EVALUATED AT T.  (T*,Y*) IS CHOSEN AS A
C     NEARBY POINT, AND THE EQUATION IS INTEGRATED FROM THERE AND
C     CHECKED FOR CONSISTENCY WITH HAVING INTEGRATED FROM T.  IF NOT,
C     A DIFFERENT (T*,Y*) IS CHOSEN UNTIL CONSISTENCY IS ACHIEVED.
C
C  INPUT QUANTITIES: NC,TSTAR,YSTAR,THEND,DTHDE
C  OUTPUT QUANTITIES: IFLAG,TT,YY,IND,Y,TOUT,ERR,WORK,IWORK,TSTAR
C
C   The criterion to be used for "consistency" here is a so-called
C   "correction formula", obtained by integrating Newton's backward
C   difference formula - which can be written as
C       Dx0 = h(q4 - (7/2)Dq3 + (53/12)DDq2 - (55/24)DDDq1 +
C                     + (251/720)DDDDq0 ....)
C   (Here Dx0 = x1-x0, Dq1 = q2-q1, etc., and the qi terms are the
C    derivatives of the function x(*) at the points yi.)  When the
C   numerical integration for x has reached y4 and there is some
C   uncertainty about the earlier value of x1 at y1, this formula
C   can be used to effect a correction.
C
C  THE BASIC IDEA IS THIS:
C    SUPPOSE THE PROBLEM IS TO INTEGRATE
C            y'(x) = f(x, y(x)),  y(a) = A
C    FROM a TO c, BUT f(a,A) is +Infinity.
C    INTERCHANGE INDEPENDENT AND DEPENDENT VARIABLES SO THAT THE
C    PROBLEM BECOMES
C           x'(y) = 1/f(x,y), x(A) = a.
C    Let yi, i=1,2,3,4, be equally spaced points, distance h apart.
C    Let x1 be an initial "guess" for the value of the solution at y1,
C    and integrate the d.e. for the variable x from y1 to y4, obtaining
C    values xi for i=1,2,3,4.  Let qi be the values for the derivative
C    function for x at the points yi.  Then, according to the above
C    correction formula, a "corrected" value of the initial guess x1 is
C      x1 = x0 + h(q4 - (7/2)Dq3 + (53/12)DDq2 - (55/24)DDDq1 +
C                   (251/720)DDDDq0 + ... )
C    This process can be repeated until, hopefully, x1 no longer
C    changes.
C    IN THE EVENT THAT THE EQUATION y'(x) = f(x, y(x)) IS NOT
C    +Infinity, BUT SIMPLY CANNOT BE EVALUATED AT x=a, THE ONLY
C    DIFFERENCE IS THAT THERE IS NO NEED TO INTERCHANGE INDEPENDENT
C    AND DEPENDENT VARIABLES.
C
C     **********
C     .. Scalars in Common ..
      REAL EIG,EPSMIN
      INTEGER IND,T21
      LOGICAL PR
C     ..
C     .. Local Scalars ..
      REAL CHM,CHNG,D2F,D2G,D3F,D3G,D4F,D4G,EPST,HT,HU,
     +                 OLDSS2,OLDYY2,ONE,RR,SLO,SOUT,SUMM,SUP,T,TEN5,
     +                 TIN,U,UOUT,USTAR,YLO,YOUT,YUP
      INTEGER I,K,KFLAG,NN3
C     ..
C     .. Local Arrays ..
      REAL DF(4),DG(4),FF(6),GG(5),S(3),SS(6,3),UU(6)
C     ..
C     .. External Subroutines ..
      EXTERNAL GERK,PQEXT
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS,MAX,SIGN
C     ..
C     .. Scalar Arguments ..
      REAL DTHDE,THEND,TOUT,TSTAR,YSTAR
      INTEGER IFLAG,NC
C     ..
C     .. Array Arguments ..
      REAL ERR(3),TT(7),WORK(27),Y(3),YY(7,3)
      INTEGER IWORK(5)
C     ..
C     .. Subroutine Arguments ..
      EXTERNAL FG
C     ..
C     .. Common blocks ..
      COMMON /DATAF/EIG,IND
      COMMON /PRIN/PR,T21
      COMMON /RNDOFF/EPSMIN
C     ..
      IFLAG = 1
      ONE = 1.0D0
      TEN5 = 100000.0D0
      EPST = 1.D-7
C
C     INTEGRATE Y' = F(T,Y,YP), STARTING AT T = TSTAR.
C
      TIN = SIGN(ONE,TSTAR)
      TT(1) = TIN
      YY(1,1) = THEND
      YY(1,2) = DTHDE
      YY(1,3) = 0.0D0

C  THE NEXT THREE VALUES FOR YY(2,*) ARE JUST INITIAL GUESSES.
      YY(2,1) = YSTAR
      YY(2,2) = DTHDE
      YY(2,3) = 0.0D0
C
      YLO = -TEN5
      YUP = TEN5
      CHM = MAX(0.01D0*ONE,0.1D0*ABS(THEND))
C
C     NORMALLY IND = 1 OR 3;  IND IS SET TO 2 OR 4 WHEN Y IS TO BE USED
C     AS THE INDEPENDENT VARIABLE IN FG(T,Y,YP), INSTEAD OF THE USUAL T.
C     BEFORE LEAVING THIS SUBROUTINE, IND IS RESET TO 1.
C
   10 CONTINUE
      T = TSTAR
      HT = TSTAR - TIN
      TT(2) = TSTAR
      Y(1) = YY(2,1)
      Y(2) = YY(2,2)
      Y(3) = YY(2,3)
      KFLAG = 1
      IND = 1
      DO 30 K = 3,6
          TOUT = T + HT
          NN3 = 0
   20     CONTINUE
          CALL GERK(FG,3,Y,T,TOUT,EPST,EPST,KFLAG,ERR,WORK,IWORK)
          IF (KFLAG.GT.3 .AND. ABS(TSTAR).GE.0.89D0) THEN
              TSTAR = TIN + 5.0D0* (TSTAR-TIN)
              IF (PR) WRITE (T21,FMT=*) ' KFLAG,NEW TSTAR = ',KFLAG,
     +            TSTAR
              GO TO 10
          END IF
          IF (KFLAG.EQ.3) THEN
              NN3 = NN3 + 1
              IF (NN3.GE.6) THEN
                  IFLAG = 0
                  RETURN
              END IF
              GO TO 20
          END IF
          YOUT = Y(1)
          TT(K) = T
          YY(K,1) = YOUT
          YY(K,2) = Y(2)
          YY(K,3) = Y(3)
   30 CONTINUE
      IND = 3
      DO 40 I = 2,6
          CALL FG(TT(I),YY(I,1),FF(I))
   40 CONTINUE
C SET FF(1) = 0., JUST TO HAVE IT DEFINED BEFORE GOING INTO PQEXT() :
C SUBROUTINE PQEXT, IF USED, TRIES TO DETERMINE WHETHER THE VALUES
C   OF FF(I) EXTRAPOLATE TO A FINITE OR INFINITE VALUE AT FF(1).
      FF(1) = 0.0D0
      IF (NC.LE.4 .OR. NC.EQ.7) THEN
          CALL PQEXT(FF(1))
      ELSE
          CALL FG(TT(1),YY(1,1),FF(1))
      END IF

      WRITE (T21,FMT=*) ' in wr, ff1 = ',FF(1)
      IF (ABS(FF(1)).LE.50.0D0) THEN
C
C     NOW WE WANT TO APPLY THE ABOVE CONSISTENCY CRITERION,
C     TO SEE IF THESE RESULTS ARE CONSISTENT WITH HAVING
C     INTEGRATED FROM (TT(1),YY(1,1).
C
          WRITE (T21,FMT=*) ' in wr, finite '
          DO 50 I = 1,4
              DF(I) = FF(I+1) - FF(I)
   50     CONTINUE
          D2F = DF(4) - DF(3)
          D3F = DF(4) - 2.0D0*DF(3) + DF(2)
          D4F = DF(4) - 3.0D0*DF(3) + 3.0D0*DF(2) - DF(1)
          SUMM = HT* (FF(5)-3.5D0*DF(4)+53.0D0*D2F/12.0D0-
     +           55.0D0*D3F/24.0D0+251.0D0*D4F/720.0D0)
C
C     PRESUMABLY, YY(2,1) SHOULD BE YY(1,1) + SUMM.
C
          OLDYY2 = YY(2,1)
C     TO COUNTER THE TENDENCY TO OVERCORRECT, USE A FACTOR
C     OF RR < 1.0 :
          RR = 0.95D0
          YY(2,1) = YY(1,1) + RR*SUMM
C
C     ALSO IMPROVE THE VALUE OF Y(2) AT TSTAR.
C
          YY(2,2) = 0.5D0* (YY(1,2)+YY(3,2))
          YY(2,3) = 0.5D0* (YY(1,3)+YY(3,3))
          CHNG = YY(2,1) - OLDYY2
          IF (CHNG.GE.0.0D0 .AND. OLDYY2.GT.YLO) YLO = OLDYY2
          IF (CHNG.LE.0.0D0 .AND. OLDYY2.LT.YUP) YUP = OLDYY2
          IF ((YY(2,1).GE.YUP.AND.YLO.GT.-TEN5) .OR.
     +        (YY(2,1).LE.YLO.AND.YUP.LT.TEN5)) YY(2,1) = 0.5D0*
     +        (YLO+YUP)
          CHNG = YY(2,1) - OLDYY2
          IF (ABS(CHNG).GT.CHM) CHNG = SIGN(CHM,CHNG)
          YY(2,1) = OLDYY2 + CHNG
c            IF(PR)WRITE(T21,*) ' YY2, CHNG = ',YY(2,1),CHNG
          IF (ABS(YY(2,1)-OLDYY2).GT.EPST) GO TO 10
          TOUT = TT(6)
      ELSE
C
C     HERE, Y' IS ASSUMED INFINITE AT T = TIN.  IN THIS CASE,
C     IT CANNOT BE EXPECTED TO APPROXIMATE Y WITH A POLYNOMIAL,
C     SO THE INDEPENDENT AND DEPENDENT VARIABLES ARE INTERCHANGED.
C     THE POINTS ARE ASSUMED EQUALLY SPACED.
C
          WRITE (T21,FMT=*) ' in wr, infinite '
          HU = (YY(6,1)-YY(1,1))/5.0D0
          UU(1) = YY(1,1)
          SS(1,1) = TT(1)
          SS(1,2) = DTHDE
          SS(1,3) = 0.0D0
          UU(2) = UU(1) + HU
          SS(2,1) = TSTAR
          SS(2,2) = DTHDE
          SS(2,3) = 0.0D0
          USTAR = UU(2)
          SLO = -TEN5
          SUP = TEN5
   60     CONTINUE
          U = USTAR
          S(1) = SS(2,1)
          S(2) = SS(2,2)
          S(3) = SS(2,3)
          KFLAG = 1
          IND = 2
          DO 80 K = 3,6
              UOUT = U + HU
              NN3 = 0
   70         CONTINUE
              CALL GERK(FG,3,S,U,UOUT,EPST,EPST,KFLAG,ERR,WORK,IWORK)
              IF (KFLAG.GT.3 .AND. ABS(TSTAR).GE.0.89D0) THEN
                  TSTAR = TIN + 5.0D0* (TSTAR-TIN)
                  GO TO 10
              END IF
              IF (KFLAG.EQ.3) THEN
                  NN3 = NN3 + 1
                  IF (NN3.GE.6) THEN
                      IFLAG = 0
                      RETURN
                  END IF
                  GO TO 70
              END IF
              SOUT = S(1)
              UU(K) = U
              SS(K,1) = SOUT
              SS(K,2) = S(2)
              SS(K,3) = S(3)
   80     CONTINUE
          IND = 4
          DO 90 I = 2,5
              CALL FG(UU(I),SS(I,1),GG(I))
   90     CONTINUE
          GG(1) = 0.0D0
          DO 100 I = 1,4
              DG(I) = GG(I+1) - GG(I)
  100     CONTINUE
          D2G = DG(4) - DG(3)
          D3G = DG(4) - 2.0D0*DG(3) + DG(2)
          D4G = DG(4) - 3.0D0*DG(3) + 3.0D0*DG(2) - DG(1)
          SUMM = HU* (GG(5)-3.5D0*DG(4)+53.0D0*D2G/12.0D0-
     +           55.0D0*D3G/24.0D0+251.0D0*D4G/720.0D0)
C
C     PRESUMABLY, SS(2,1) SHOULD BE SS(1,1) + SUMM.
C
          OLDSS2 = SS(2,1)
          SS(2,1) = SS(1,1) + SUMM
          IF (SS(2,1).LE.-1.0D0) SS(2,1) = -1.0D0 + EPSMIN
          IF (SS(2,1).GE.1.0D0) SS(2,1) = 1.0D0 - EPSMIN
C
C     ALSO IMPROVE THE VALUE OF Y(2) AT TSTAR.
C
          SS(2,2) = 0.5D0* (SS(1,2)+SS(3,2))
          SS(2,3) = 0.5D0* (SS(1,3)+SS(3,3))
          CHNG = SS(2,1) - OLDSS2
C                IF(PR)WRITE(T21,*) ' CHNG = ',CHNG
          IF (CHNG.GE.0.0D0 .AND. OLDSS2.GT.SLO) SLO = OLDSS2
          IF (CHNG.LE.0.0D0 .AND. OLDSS2.LT.SUP) SUP = OLDSS2
          IF ((SS(2,1).GE.SUP.AND.SLO.GT.-TEN5) .OR.
     +        (SS(2,1).LE.SLO.AND.SUP.LT.TEN5)) SS(2,1) = 0.5D0*
     +        (SLO+SUP)
          IF (ABS(SS(2,1)-OLDSS2).GT.EPST) GO TO 60
      END IF

      IF (IND.EQ.4) THEN

C  NOW INTEGRATE AGAIN BUT WITH T AS THE INDEPENDENT VARIABLE:
C  THIS IS USEFUL FOR OBTAINING THE USUAL GLOBAL ERROR
C  ESTIMATES FOR THE QUANTITIES Y(1),Y(2),Y(3).
  120     CONTINUE
          TT(6) = SS(6,1)
          YY(6,1) = UU(6)
          YY(6,2) = SS(6,2)
          YY(6,3) = SS(6,3)
          T = TT(6)
          HT = T - TIN
          Y(1) = YY(6,1)
          Y(2) = YY(6,2)
          Y(3) = YY(6,3)
          KFLAG = 1
          IND = 1

          DO 130 K = 7,12
              TOUT = T + HT
              CALL GERK(FG,3,Y,T,TOUT,EPST,EPST,KFLAG,ERR,WORK,IWORK)
              IF (KFLAG.EQ.5) THEN
                  EPST = 5.0D0*EPST
                  WRITE (T21,FMT=*) ' in wr, epst = ',EPST
                  GO TO 120
              END IF
  130     CONTINUE

          TOUT = T
      END IF

C     TOUT = TT(6)
      IND = 1
      RETURN
      END
C
      SUBROUTINE SETTHU(X,THU)
C     **********
C
C  THIS SUBROUTINE IS CALLED ONLY FROM SUBROUTINE LCO.
C    IT ESTABLISHES A DEFINITE VALUE FOR THU,
C    THE PHASE ANGLE FOR THE FUNCTION U, INCLUDING AN
C    APPROPRIATE INTEGER MULTIPLE OF PI
C    IT NEEDS THE NUMBERS MMW(*) FOUND IN THUM
C
C  INPUT QUANTITIES: X,THU,YS,MMW,MMWD
C  OUTPUT QUANTITIES: THU
C
C     **********
C     .. Scalars in Common ..
      REAL HPI,PI,TWOPI
      INTEGER MMWD
C     ..
C     .. Arrays in Common ..
      REAL YS(200)
      INTEGER MMW(100)
C     ..
C     .. Local Scalars ..
      INTEGER I
C     ..
C     .. Common blocks ..
      COMMON /PASS/YS,MMW,MMWD
      COMMON /PIE/PI,TWOPI,HPI
C     ..
C     .. Scalar Arguments ..
      REAL THU,X
C     ..
      DO 10 I = 1,MMWD
          IF (X.GE.YS(MMW(I)) .AND. X.LE.YS(MMW(I)+1)) THEN
              IF (THU.GT.PI) THEN
                  THU = THU + (I-1)*TWOPI
                  RETURN
              ELSE
                  THU = THU + I*TWOPI
                  RETURN
              END IF
          END IF
   10 CONTINUE
      DO 20 I = 1,MMWD
          IF (X.GE.YS(MMW(I))) THU = THU + TWOPI
   20 CONTINUE
      RETURN
      END
C
      SUBROUTINE FZ(UU,Y,YP)
C     **********
C
C  THIS SUBROUTINE EVALUATES THE DERIVATIVE OF THE FUNCTION PHI,
C    PLUS THOSE OF ITS COMPANION FUNCTIONS SIGMA AND D(PHI)/D(LAMBDA),
C    FOR INTEGRATION BY GERK.  THESE INTEGRATIONS ARE CALLED FOR
C    BY THE SUBROUTINES LCNO AND LCO.
C
C  INPUT QUANTITIES: EIG,IND,UU,Y
C  OUTPUT QUANTITIES: YP
C
C  THE FUNCTIONS PHI AND SIGMA RESULT FROM REGULARIZING THE GIVEN
C    STURM-LIOUVILLE DIFFERENTIAL EQUATION NEAR A LIMIT CIRCLE
C    ENDPOINT.  NAMELY, WRITING
C      Y = Z1*U + Z2*V
C      PY' = Z1*PU' + Z2*PV',
C    WHERE U,V ARE THE BOUNDARY CONDITION FUNCTIONS NEAR THE LC
C    ENDPOINT, FOLLOWED BY
C      Z1 = SIGMA*COS(PHI)
C      Z2 = SIGMA*SIN(PHI),
C    GIVES THE FOLLOWING DIFFERENTIAL EQUATIONS FOR PHI, SIGMA,
C    AND D(PHI)/D(LAMBDA):
C    [u,v]*phi' = -a1122*sin(phi)*cos(phi) +
C                   a21*cos(phi)**2 -a12*sin(phi)**2
C    [u,v]*sigma'/sigma = (a12+a21)*sin(phi)*cos(phi) +
C                   a11*cos(phi)**2 + a22*sin(phi)**2
C    [u,v]*dphide' = -{a1122*(cos(phi)**2 - sin(phi)**2) +
C                   2*(a12+a21)*sin(phi)*cos(phi)}*dphide -
C                   {2*w*u*v*sin(phi)*cos(phi) +
C                    w*v**2*sin(phi)**2 + w*u**2*cos(phi)**2}
C    WHERE
C       a11 = eig*w*u*v - v*Hu,  a12 = eig*w*v**2 - v*Hv
C       a21 = -eig*w*u**2 + u*Hu,  a22 = -eig*w*u*v + uHv .
C       a1122 = u*(eig*w*v - Hv) + v*(eig*w*u - Hu)
C
C    HERE, Hu MEANS -(pu')' + qu,
C    AND [u,v] means u*pv' - v*pu' .
C
C     **********
C     .. Scalars in Common ..
      REAL EIG
      INTEGER IND
C     ..
C     .. Local Scalars ..
      REAL A1122,A12,A21,AU,AV,B1122,B12,B21,C,C2,D,DT,HU,
     +                 HV,PHI,PUP,PVP,S,S2,SC,T,U,V,WW,WX,X
C     ..
C     .. External Functions ..
      REAL W
      EXTERNAL W
C     ..
C     .. External Subroutines ..
      EXTERNAL DXDT,UV
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC COS,MOD,SIN
C     ..
C     .. Common blocks ..
      COMMON /DATAF/EIG,IND
C     ..
C     .. Scalar Arguments ..
      REAL UU
C     ..
C     .. Array Arguments ..
      REAL Y(2),YP(3)
C     ..
      IF (MOD(IND,2).EQ.1) THEN
          T = UU
          PHI = Y(1)
      ELSE
          T = Y(1)
          PHI = UU
      END IF
      CALL DXDT(T,DT,X)
      CALL UV(X,U,PUP,V,PVP,HU,HV)
      D = U*PVP - V*PUP
      WX = W(X)
      S = SIN(PHI)
      C = COS(PHI)
      S2 = S*S
      C2 = C*C
      SC = S*C
      B1122 = WX*U*V
      B12 = WX*V*V
      B21 = -WX*U*U
      AU = EIG*WX*U - HU
      AV = EIG*WX*V - HV
      A1122 = U*AV + V*AU
      A12 = V*AV
      A21 = -U*AU
      YP(1) = -DT* (A1122*SC+A12*S2-A21*C2)/D
      IF (IND.EQ.1) THEN
          WW = 2.0D0* (A12+A21)*SC + A1122* (C2-S2)
          YP(2) = -DT* (WW*Y(2)+2.0D0*B1122*SC+B12*S2-B21*C2)/D
          YP(3) = 0.5D0*DT*WW/D
      ELSE IF (IND.EQ.2) THEN
C  IN THIS CASE THE INDEPENDENT AND DEPENDENT VARIABLES
C    (t and phi, etc.) HAVE BEEN INTERCHANGED.
          YP(2) = YP(2)/YP(1)
          YP(3) = YP(3)/YP(1)
          YP(1) = 1.0D0/YP(1)
      ELSE IF (IND.EQ.3) THEN
      ELSE
          YP(1) = 1.0D0/YP(1)
      END IF
      RETURN
      END
C
      SUBROUTINE GERKZ(F,NEQ,Y,TIN,TOUT,REPS,AEPS,LFLAG,ERY,WORK,IWORK)
C     **********

C  THIS PROGRAM CONTROLS THE INTEGRATION OF THE DIFFERENTIAL EQUATIONS
C    FOR THETA, RHO (OR F = LOG(RHO)), AND DTHDE = D(THETA)/D(LAMBDA)
C    WHERE Y = RHO*SIN(THETA) AND PY' = Z*RHO*COS(THETA), WITH
C    SUITABLY CHOSEN CONSTANT Z.
C    IT IS CALLED FROM INTEGZ.
C
C  INPUT QUANTITIES: TIN,TOUT,Y,REPS,AEPS,ERY,EPSMIN,TEE,ZEE,
C                    WORK,IWORK
C  OUTPUT QUANTITIES: Y,LFLAG,ERY
C
C  ACTUALLY, DIFFERENT CONSTANTS Z ARE USED IN DIFFERENT SUBINTERVALS
C    OF THE INTEGRATION.  THE Z's HAVE BEEN CHOSEN ELSEWHERE AND HAVE
C    BEEN STORED IN THE ARRAY ZEE.  THE CONSTANT ZEE(J) IS USED IN
C    THE T-SUBINTERVAL (TEE(J),TEE(J+1)).
C  N.B. THE EXACT T-SUBINTERVALS USED ARE NOT AT ALL CRITICAL.
C  THE PRUEFER ANGLE, CALLED THETA WHEN Z = 1 (THE USUAL PRUEFER
C    ANGLE), IS HERE CALLED THZ WHEN Z MAY NOT BE 1.
C    THE SUBROUTINES TH2THZ AND THZ2TH ARE USED TO CONVERT FROM THE
C    ONE ANGLE TO THE OTHER.
C  THE VALUE OF Z IN COMMON/Z1/Z IS ALWAYS EQUAL TO 1 EXCEPT
C    DURING THE INTEGRATIONS IN THIS SUBROUTINE.  SO WE ALWAYS
C    SET Z = 1.0 BEFORE LEAVING HERE.
C  IT IS ALWAYS ASSUMED THAT WHEN THIS PROGRAM IS CALLED, THE
C    ARRAY ERY HAS MEANINGFUL VALUES IN IT.

C     **********
C     .. Scalars in Common ..
      REAL EPSMIN,Z
      INTEGER T21
      LOGICAL PR
C     ..
C     .. Local Scalars ..
      REAL T,TOUTS,Y3,Y3S
      INTEGER I,J,K,L,LLFLAG,NK
C     ..
C     .. Arrays in Common ..
      REAL TEE(100),ZEE(100)
      INTEGER JAY(100)
C     ..
C     .. Local Arrays ..
      REAL ERT(3),ERU(3),U(3)
C     ..
C     .. External Subroutines ..
      EXTERNAL ERRZ,GERK,TH2THZ,THZ2TH
C     ..
C     .. Common blocks ..
      COMMON /PRIN/PR,T21
      COMMON /RNDOFF/EPSMIN
      COMMON /TEEZ/TEE
      COMMON /Z1/Z
      COMMON /ZEEZ/JAY,ZEE
C     ..
C     .. Scalar Arguments ..
      REAL AEPS,REPS,TIN,TOUT
      INTEGER LFLAG,NEQ
C     ..
C     .. Array Arguments ..
      REAL ERY(3),WORK(27),Y(3)
      INTEGER IWORK(5)
C     ..
C     .. Subroutine Arguments ..
      EXTERNAL F
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC MAX,MIN
C     ..
C  SET ARRAY ERT TO THE INCOMING ARRAY ERY:
      ERT(1) = ERY(1)
      ERT(2) = ERY(2)
      ERT(3) = ERY(3)

      T = TIN
      J = 1
      L = 1
      IF (TIN.LT.TOUT) THEN
          DO 10 I = 1,19
              IF (TEE(I)-EPSMIN.LE.TIN .AND.
     +            TIN.LT.TEE(I+1)+EPSMIN) J = I
              IF (TEE(I)-EPSMIN.LT.TOUT .AND.
     +            TOUT.LE.TEE(I+1)+EPSMIN) L = I
   10     CONTINUE
          DO 30 K = J,L
              TOUTS = MIN(TOUT,TEE(K+1))
              Z = ZEE(K+1)
              IF (Z.EQ.0.0D0) Z = 1.0D0
              CALL TH2THZ(Y,Z,U)
              LLFLAG = 1
              NK = 0
   20         CONTINUE
              CALL GERK(F,NEQ,U,T,TOUTS,REPS,AEPS,LLFLAG,ERU,WORK,IWORK)
              IF (LLFLAG.GT.3) THEN
                  IF (PR) WRITE (T21,FMT=*) ' LLFLAG = ',LLFLAG
                  LFLAG = 5
                  RETURN
              END IF
              IF (LLFLAG.EQ.3 .OR. LLFLAG.EQ.-2) THEN
                  NK = NK + 1
                  IF (NK.GE.10) THEN
                      LFLAG = 5
                      RETURN
                  END IF
                  GO TO 20
              END IF
              Y3S = Y(3)
              CALL THZ2TH(U,ERU,Z,Y,ERY)
              Y3 = Y(3)
              CALL ERRZ(ERT,Y3S,Y3,ERY)
   30     CONTINUE
      ELSE
          DO 40 I = 20,2,-1
              IF (TEE(I-1)-EPSMIN.LT.TIN .AND.
     +            TIN.LE.TEE(I)+EPSMIN) J = I
              IF (TEE(I-1)-EPSMIN.LE.TOUT .AND.
     +            TOUT.LT.TEE(I)+EPSMIN) L = I
   40     CONTINUE
          DO 60 K = J,L,-1
              TOUTS = MAX(TOUT,TEE(K-1))
              Z = ZEE(K)
              IF (Z.EQ.0.0D0) Z = 1.0D0
              CALL TH2THZ(Y,Z,U)
              LLFLAG = 1
              NK = 0
   50         CONTINUE
              CALL GERK(F,NEQ,U,T,TOUTS,REPS,AEPS,LLFLAG,ERU,WORK,IWORK)
              IF (LLFLAG.GT.3) THEN
                  IF (PR) WRITE (T21,FMT=*) ' LLFLAG = ',LLFLAG
                  LFLAG = 5
                  RETURN
              END IF
              IF (LLFLAG.EQ.3 .OR. LLFLAG.EQ.-2) THEN
                  NK = NK + 1
                  IF (NK.GE.10) THEN
                      LFLAG = 5
                      RETURN
                  END IF
                  GO TO 50
              END IF
              Y3S = Y(3)
              CALL THZ2TH(U,ERU,Z,Y,ERY)
              Y3 = Y(3)
              CALL ERRZ(ERT,Y3S,Y3,ERY)
   60     CONTINUE
      END IF
      TIN = T

C  BEFORE LEAVING THIS ROUTINE WE WANT TO RESET THE PARAMETER Z
C    (STORED IN COMMON/Z1/Z) TO BE 1.0 AGAIN.

      Z = 1.0D0
      LFLAG = LLFLAG
      RETURN
      END
C
      SUBROUTINE ERRZ(ERT,Y3S,Y3,ERY)
C  THIS PROGRAM COMPUTES AN ESTIMATE OF THE GLOBAL ERROR,
C  ERY, IN Y IN INTEGRATING FROM X1 TO X2, WHEN THERE
C  IS AN ERROR, ERT, IN THE INITIAL VALUE OF Y AT X1.
C  IT IS CALLED FROM INTEGZ.
C
C  INPUT QUANTITIES: ERT,Y3S,Y3,ERY
C  OUTPUT QUANTITIES: ERY
C
C  USING THE FACT THAT F (= LOG(RHO)) SATISFIES
C      F' = f
C  WHILE IF G REPRESENTS D(THETA)/D(INITIAL VALUE), THEN
C      G' = -2*f*G,
C  SO
C      AN ERROR OF ERT AT X1 BECOMES AN ERROR OF
C            ERT*EXP(-2*(F(X2)-F(X1)))
C  WHICH MUST BE ADDED TO THE GLOBAL ERROR, ERY.
C  WHEN THIS SUBROUTINE IS CALLED, Y3 IS THE CURRENT
C  VALUE OF Y(3) (OR F), AND Y3S IS THE VALUE IT
C  HAD AT THE BEGINNING OF THE LAST INTEGRATION, AT X1.
C  ON EXIT, ERY IS THE CURRENT ESTIMATE OF THE GLOBAL
C  ERROR.
C
C     .. Scalar Arguments ..
      REAL Y3,Y3S
C     ..
C     .. Array Arguments ..
      REAL ERT(3),ERY(3)
C     ..
C     .. Local Scalars ..
      INTEGER I
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC EXP
C     ..
      DO 10 I = 1,3
          ERY(I) = ERT(I)*EXP(-2.0D0* (Y3-Y3S)) + ERY(I)
   10 CONTINUE
      RETURN
      END
C
      SUBROUTINE LCNO(TEND,THEND,DTHDAA,TM,COF1,COF2,EPS,Y,ER,IFLAG)
C  THIS PROGRAM IS CALLED FROM INTEG.
C
C  INPUT QUANTITIES: TEND,THEND,TM,COF1,COF2,EPS
C  OUTPUT QUANTITIES: Y,THEND,DTHDAA,TT,YY,ER,IFLAG
C
C  IT CONTROLS THE INTEGRATIONS WHICH BEGIN AT AN
C    LCNO ENDPOINT.  NEAR SUCH AN ENDPOINT, THE EQUATION
C         -(py')' + q*y = lambda*w*y
C    IS TRANSFORMED BY SETTING
C        y = sigma*(u*cos(phi) + v*sin(phi))
C       py' = sigma*(pu'*cos(phi) + pv'*sin(phi))
C    WHERE u,v ARE THE BOUNDARY CONDITION FUNCTIONS.  THE RESULTING
C    DIFFERENTIAL EQUATIONS FOR sigma, phi, d(phi)/d(lambda) ARE
C    FIRST INTEGRATED FROM THE ENDPOINT TO A POINT FAR ENOUGH AWAY,
C    AND ARE THEN TRANSLATED TO THE USUAL PRUEFER VARIABLES
C    THETA, RHO, ETC., WHICH ARE THEN INTEGRATED TO TM.


C     .. Scalar Arguments ..
      REAL COF1,COF2,DTHDAA,EPS,TEND,THEND,TM
      INTEGER IFLAG
C     ..
C     .. Array Arguments ..
      REAL ER(3),Y(3)
C     ..
C     .. Scalars in Common ..
      REAL HPI,PI,TWOPI
      INTEGER T21
      LOGICAL PR
C     ..
C     .. Arrays in Common ..
      REAL TT(7,2),YY(7,3,2)
      INTEGER NT(2)
C     ..
C     .. Local Scalars ..
      REAL C,D,DDD,DPHIDE,DTHIN,DUM,EFF,FAC2,FRA,HU,HV,ONE,
     +                 PHI,PHI0,PHIDAA,PUP,PVP,PYPZ,PYPZ0,RHOSQ,S,T,TH,
     +                 TH0,THBAR,THIN,TIN,TMP,TOUT,TSTAR,U,V,XT,YZ,YZ0
      INTEGER I,J,K2PI,KFLAG,M,NC
      CHARACTER*16 PREC
C     ..
C     .. Local Arrays ..
      REAL WORK(27),YP(3)
      INTEGER IWORK(5)
C     ..
C     .. External Subroutines ..
      EXTERNAL DXDT,FZ,GERK,INTEGZ,UV,WR
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS,ATAN2,COS,EXP,LOG,MIN,SIN
C     ..
C     .. Common blocks ..
      COMMON /PIE/PI,TWOPI,HPI
      COMMON /PRIN/PR,T21
      COMMON /TEMP/TT,YY,NT
C     ..
      ONE = 1.0D0
      PREC = '  REAL'
      NC = 3

      T = TEND

C  THE LCNO BOUNDARY CONDITION A1*[u,y] + A2*[v,y] = 0
C  IS EQUIVALENT TO COF1*sin(phi) - COF2*cos(phi) = 0.  SO
      PHI0 = ATAN2(COF2,COF1)
C
C     WE WANT -PI/2 .LT. PHI0 .LE. PI/2.
C
      Y(1) = PHI0
      Y(2) = 0.0D0
      Y(3) = 0.0D0
      CALL DXDT(T,TMP,XT)
      CALL FZ(T,Y,YP)
      PHIDAA = -YP(1)
      CALL UV(XT,U,PUP,V,PVP,HU,HV)
      D = U*PVP - V*PUP
      YZ0 = U*COF1 + V*COF2
      PYPZ0 = (PUP*COF1+PVP*COF2)

C   USING THE RELATION BETWEEN theta AND phi, NAMELY
C  tan(theta) = (u*cos(phi) + v*sin(phi))/(pup*cos(phi)+pv'*sin(phi))
C  AND, DIFFERENTIATING W.R.T. lambda or a or b,
C  d(theta)*sec(theta)**2 = -[u,v]*d(phi)/(pu'*cos(phi)+pv'*sin(phi))**2
C
C     SET TH0 AND COPY INTO THEND, OVERWRITING ITS INPUT VALUE.
C     ALSO, REDEFINE DTHDAA, OVERWRITING ITS INPUT VALUE.
C
      TH0 = ATAN2(YZ0,PYPZ0)
      IF (YZ0.LT.0.0D0) TH0 = TH0 + TWOPI
      THEND = TH0
      IF (PR) WRITE (T21,FMT=*) ' PHI0,TH0 = ',PHI0,TH0
      DUM = ABS(COS(TH0))

      IF (DUM.GE.0.5D0) THEN
          FAC2 = (-D)* (DUM/PYPZ0)**2
      ELSE
          FAC2 = (-D)* (SIN(TH0)/YZ0)**2
      END IF

C  PHIDAA REPRESENTS d(phi)/da  EVALUATED AT THE ENDPOINT (A or B).
      DTHDAA = FAC2*PHIDAA
C
C     IN THE NEXT PIECE, WE ASSUME TH0 .GE. 0.
C
      M = 0
      IF (TH0.EQ.0.0D0) M = -1
      IF (TH0.GT.PI .OR. (TH0.EQ.PI.AND.T.LT.TM)) M = 1
      PHI = PHI0
      K2PI = 0
      YZ0 = U*COS(PHI0) + V*SIN(PHI0)

C  BEGIN THE INTEGRATION FOR PHI,SIGMA, ETC. AT THE LCNO ENDPOINT
C    USING SUBROUTINE WR.
C  THE FIRST 6 VALUES OF PHI,ETC. AT T = TT(K,J) ARE YY(1,K,J),
C    K = 1,..,6.

      IF (TM.GT.TEND) THEN
          J = 1
          TSTAR = -0.99999D0
          IF (PREC(3:3).NE.'R') TSTAR = -0.9999999999D0
      ELSE
          J = 2
          TSTAR = 0.99999D0
          IF (PREC(3:3).NE.'R') TSTAR = 0.9999999999D0
      END IF
      DPHIDE = 0.0D0

        DO 23 I = 1,27
           WORK(I) = 0.0D0
           IF(I.LE.5) IWORK(I) = 0
  23    CONTINUE

      CALL WR(FZ,NC,TSTAR,PHI0,PHI0,DPHIDE,TOUT,Y,TT(1,J),YY(1,1,J),
     +        KFLAG,ER,WORK,IWORK)
      IF (KFLAG.NE.1) THEN
          IF (PR) WRITE (T21,FMT=*) ' KFLAG = 0 '
          IFLAG = 5
          RETURN
      END IF

C  THE INTEGRATION FOR PHI, ETC. HAS BEEN STARTED, BUT NEEDS TO
C    CONTINUE A LITTLE FARTHER TO GET WELL AWAY FROM THE ENDPOINT.
C    FOR THIS, WE CAN JUST CONTINUE USING KFLAG = 2.
      TIN = TOUT
      TMP = 0.01D0
      DDD = MIN(TMP,ABS(TOUT))
      TOUT = TOUT + DDD* (-TOUT)/ABS(TOUT)
      T = TIN
      KFLAG = 2
   60 CONTINUE
      CALL GERK(FZ,3,Y,T,TOUT,EPS,EPS,KFLAG,ER,WORK,IWORK)
      IF (KFLAG.GT.3) THEN
          IF (PR) WRITE (T21,FMT=*) ' KFLAG2 = ',KFLAG
          IFLAG = 5
          RETURN
      END IF
      IF (KFLAG.EQ.3) THEN
          FRA = ABS(T-TIN)/ABS(T-TOUT)
          IF (FRA.LT.0.001D0) THEN
              IF (PR) WRITE (T21,FMT=*) ' KFLAG2 = ',KFLAG
              IFLAG = 5
              RETURN
          END IF
          GO TO 60
      END IF
      IF (T.EQ.TOUT) THEN
C  STORE THE 7th VALUES FOR T, PHI,ETC.
          TT(7,J) = T
          YY(7,1,J) = Y(1)
          YY(7,2,J) = Y(2)
          YY(7,3,J) = Y(3)
      END IF
      T = TOUT
      CALL DXDT(T,TMP,XT)
      CALL UV(XT,U,PUP,V,PVP,HU,HV)
      PHI = Y(1)
      S = SIN(PHI)
      C = COS(PHI)
      YZ = U*C + V*S
      IF (YZ*YZ0.LT.0.0D0) K2PI = K2PI + 1
      YZ0 = YZ
C
C     CONVERT FROM PHI TO THETA.
C
      DPHIDE = Y(2)
      D = U*PVP - V*PUP
      PYPZ = (PUP*C+PVP*S)
      THBAR = ATAN2(YZ,PYPZ)
      IF (TM.GT.TEND .AND. THBAR.LT.TH0 .AND.
     +    PHI.LT.PHI0) THBAR = THBAR + TWOPI
      IF (TM.LT.TEND .AND. THBAR.GT.TH0 .AND.
     +    PHI.GT.PHI0) THBAR = THBAR - TWOPI
      TH = THBAR - M*PI
      IF (TH.LT.-PI) TH = TH + TWOPI
      IF (TH.GT.TWOPI) TH = TH - TWOPI
      IF (TM.LT.TEND .AND. K2PI.GT.1) TH = TH - (K2PI-1)*TWOPI
      IF (TM.GT.TEND .AND. K2PI.GT.1) TH = TH + (K2PI-1)*TWOPI
      IF (TM.GT.TEND .AND. TH*TH0.LT.0.0D0) TH = TH + TWOPI
      IF (PR) WRITE (T21,FMT=*) ' PHI,TH = ',PHI,TH
C
C     WE NOW HAVE YZ, PYPZ, PHI AND TH, SO WE CAN GET DTHDE
C       FROM DPHIDE.
C
      DUM = ABS(COS(TH))
      IF (DUM.GE.0.5D0) THEN
          FAC2 = - (D)* (DUM/PYPZ)**2
      ELSE
          FAC2 = - (D)* (SIN(TH)/YZ)**2
      END IF
      DTHIN = FAC2*DPHIDE
C  ALSO CONVERT THE ESTIMATED INTEGRATION ERRORS FOR PHI,ETC.,
C    INTO CORRESPONDING ERRORS FOR THETA,ETC.
      ER(1) = FAC2*ER(1)
      ER(2) = FAC2*ER(2)
      ER(3) = FAC2*ER(2)
C
      RHOSQ = EXP(2.0D0*Y(3))* (YZ**2+PYPZ**2)
      EFF = 0.5D0*LOG(RHOSQ)
C           END (INTEGRATE-FOR-PHI-NONOSC)

C  NOW COMPLETE THE INTEGRATION FOR THETA, RHO, ETC. TO TM.
      TIN = TOUT
      TOUT = TM
      THIN = TH
C        DO (INTEGRATE-FOR-THETA)
      CALL INTEGZ(TIN,TOUT,THIN,DTHIN,EPS,Y,IFLAG,ER,WORK,IWORK)
C
C           END (INTEGRATE-FOR-THETA)
      Y(3) = Y(3) + EFF
      RETURN
      END
C
      SUBROUTINE REG(TEND,THEND,TM,DTHDE,EPS,Y,ER,IFLAG)
C
C     THIS PROGRAM IS CALLED FROM INTEG.
C     THIS IS THE REGULAR (NOT WEAKLY REGULAR) CASE,
C     SO THERE IS NO ENDPOINT PROBLEM AT ALL, AND INTEGZ
C     CAN BE CALLED DIRECTLY.
C
C  INPUT QUANTITIES: TEND,THEND,TM,DTHDE,EPS,Y
C  OUTPUT QUANTITIES: Y,ER,IFLAG
C
C     .. Scalar Arguments ..
      REAL DTHDE,EPS,TEND,THEND,TM
      INTEGER IFLAG
C     ..
C     .. Array Arguments ..
      REAL ER(3),Y(3)
C     ..
C     .. Local Scalars ..
      REAL DTHIN,THIN,TIN,TOUT
C     ..
C     .. Local Arrays ..
      REAL WORK(27)
      INTEGER IWORK(5)
C     ..
C     .. External Subroutines ..
      EXTERNAL INTEGZ
C     ..
      TIN = TEND
      TOUT = TM
      THIN = THEND
      DTHIN = DTHDE
C        DO (INTEGRATE-FOR-THETA)
C        INITIALIZE ER BEFORE CALLING INTEGZ:
      ER(1) = 0.0D0
      ER(2) = 0.0D0
      ER(3) = 0.0D0
      CALL INTEGZ(TIN,TOUT,THIN,DTHIN,EPS,Y,IFLAG,ER,WORK,IWORK)
C           END (INTEGRATE-FOR-THETA)
      RETURN
      END
C
      SUBROUTINE WREG(TEND,THEND,DTHDE,TM,EPS,ER,Y,IFLAG)
C    THIS PROGRAM IS CALLED FROM INTEG.
C
C  INPUT QUANTITIES: TEND,THEND,DTHDE,TM,EPS,ER
C  OUTPUT QUANTITIES: Y,ER,IFLAG,TT,YY,NT
C
C    THIS IS THE 'WEAKLY REGULAR' CASE, SO SUBROUTINE WR MUST
C    BE USED FOR THE INTEGRATIONS STARTING AT THE WR ENDPOINT
C    UNTIL THE INTEGRATIONS HAVE REACHED A POINT FAR ENOUGH
C    AWAY THAT THE WEAKLY REGULAR ENDPOINT NO LONGER PRESENTS
C    A DIFFICULTY, AND THE SUBROUTINE INTEGZ CAN BE USED FROM
C    THERE ON.
C     .. Scalar Arguments ..
      REAL DTHDE,EPS,TEND,THEND,TM
      INTEGER IFLAG
C     ..
C     .. Array Arguments ..
      REAL ER(3),Y(3)
C     ..
C     .. Scalars in Common ..
      REAL EIG,HPI,PI,TWOPI
      INTEGER IND,T21
      LOGICAL PR
C     ..
C     .. Arrays in Common ..
      REAL TT(7,2),YY(7,3,2)
      INTEGER NT(2)
C     ..
C     .. Local Scalars ..
      REAL DDD,DTHIN,EFF,ONE,P1,Q1,T,THIN,TIN,TMP,TOUT,
     +                 TSTAR,W1,XSTAR,YSTAR
      INTEGER I,J,KFLAG,NC
      CHARACTER*16 PREC
C     ..
C     .. Local Arrays ..
      REAL WORK(27),YU(3)
      INTEGER IWORK(5)
C     ..
C     .. External Functions ..
      REAL P,Q,W
      EXTERNAL P,Q,W
C     ..
C     .. External Subroutines ..
      EXTERNAL DXDT,F,GERK,INTEGZ,WR
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS,COS,MIN,SIN
C     ..
C     .. Common blocks ..
      COMMON /DATAF/EIG,IND
      COMMON /PIE/PI,TWOPI,HPI
      COMMON /PRIN/PR,T21
      COMMON /TEMP/TT,YY,NT
C     ..
      PREC = '  REAL'
      ONE = 1.0D0
C
      NC = 2
      IF (TM.GT.TEND) THEN
          J = 1
          TSTAR = -0.99999D0
          IF (PREC(3:3).NE.'R') TSTAR = -0.9999999999D0
      ELSE
          J = 2
          TSTAR = 0.99999D0
          IF (PREC(3:3).NE.'R') TSTAR = 0.9999999999D0
      END IF
C  FORM AN ESTIMATE FOR YSTAR TO BE USED IN THE CALL TO WR.
      CALL DXDT(TSTAR,TMP,XSTAR)
      P1 = 1.0D0/P(XSTAR)
      Q1 = Q(XSTAR)
      W1 = W(XSTAR)
      YSTAR = THEND + 0.5D0* (TSTAR-TEND)*
     +        (P1*COS(THEND)**2+ (EIG*W1-Q1)*SIN(THEND)**2)

        DO 23 I = 1,27
           WORK(I) = 0.0D0
           IF(I.LE.5) IWORK(I) = 0
  23    CONTINUE

      CALL WR(F,NC,TSTAR,YSTAR,THEND,DTHDE,TOUT,YU,TT(1,J),YY(1,1,J),
     +        KFLAG,ER,WORK,IWORK)
C  USE THE SAME ARRAY, ER, AS IN THE NEXT CALL TO GERK, SO THAT THE
C  ERROR MEASUREMENT IS CONTINUED.
      IF (KFLAG.NE.1) THEN
          IF (PR) WRITE (T21,FMT=*) ' KFLAG = 0 '
          IFLAG = 5
          RETURN
      END IF
C  CONTINUE THE INTEGRATIONS TO A POINT FARTHER AWAY FROM THE
C  WR ENDPOINT.
      T = TOUT
      TMP = 0.01D0
      DDD = MIN(TMP,ABS(TOUT))
      TOUT = TOUT + DDD* (-TOUT)/ABS(TOUT)
C  IN ORDER TO JUST CONTINUE THE INTEGRATION THAT WAS BEING
C  DONE IN SUBROUTINE WR(), SET KFLAG = 2.
      KFLAG = 2
      CALL GERK(F,3,YU,T,TOUT,EPS,EPS,KFLAG,ER,WORK,IWORK)
      IF (KFLAG.EQ.5) THEN
          IFLAG = 5
          GO TO 50
      END IF
C  THE FIRST 6 VALUES OF PHI,ETC. AT T = TT(K,J) ARE YY(1,K,J),
C    K = 1,..,6.
C  STORE THE 7th VALUE FOR T, THETA, ETC.
      TT(7,J) = T
      YY(7,1,J) = YU(1)
      YY(7,2,J) = YU(2)
      YY(7,3,J) = YU(3)
C  NOW CONTINUE THE INTEGRATIONS FOR THETA, ETC., TO TM.
      TIN = TOUT
      TOUT = TM
      THIN = YU(1)
      DTHIN = YU(2)
      EFF = YU(3)
      CALL INTEGZ(TIN,TOUT,THIN,DTHIN,EPS,Y,IFLAG,ER,WORK,IWORK)
      Y(3) = EFF + Y(3)
   50 CONTINUE
      RETURN
      END
C
      SUBROUTINE LCO(TEND,THEND,DTHDAA,DTHDE,TM,COF1,COF2,EPS,Y,ER,OK,
     +               IFLAG)

C  THIS PROGRAM CONTROLS THE INTEGRATIONS WHICH BEGIN NEAR AN
C    LCO ENDPOINT.  IT IS CALLED FROM INTEG.
C
C  INPUT QUANTITIES: OK,TEND,THEND,COF1,COF2,DTHDAA,TM,EPS,
C                    TSAVEL,TSAVER,EIG,ISAVE,Y
C  OUTPUT QUANTITIES: THEND,DTHDAA,TT,YY,NT,TSAVEL,TSAVER,
C                     Y,ER,IFLAG
C
C    NEAR SUCH AN ENDPOINT, THE EQUATION
C         -(py')' + q*y = lambda*w*y
C    IS TRANSFORMED BY SETTING
C        y = sigma*(u*cos(phi) + v*sin(phi))
C       py' = sigma*(pu'*cos(phi) + pv'*sin(phi))
C    WHERE u,v ARE THE BOUNDARY CONDITION FUNCTIONS.  THE RESULTING
C    DIFFERENTIAL EQUATIONS FOR sigma, phi, d(phi)/d(lambda) ARE
C    FIRST INTEGRATED FROM NEAR THE ENDPOINT TO A POINT FAR ENOUGH AWAY,
C    AND ARE THEN TRANSLATED TO THE USUAL PRUEFER VARIABLES
C    THETA, RHO, ETC., WHICH ARE THEN INTEGRATED TO TM.
C  BECAUSE THE INTEGRATIONS MUST BE CARRIED OUT NOT ONLY WHEN LAMBDA
C    IS SET EQUAL TO EIG, BUT ALSO WHEN LAMBDA IS SET EQUAL TO 0.0,
C    IT IS NECESSARY TO MAKE SURE THE INTEGRATIONS IN BOTH CASES
C    ARE CARRIED OUT OVER EXACTLY THE SAME SUBINTERVALS.
C    THE T-VALUES WHERE THE CHANGE FROM PHI TO THETA TAKES PLACE
C    ARE DENOTED AND SAVED AS THE VARIABLES TSAVEL AND TSAVER.


C     .. Scalar Arguments ..
      REAL COF1,COF2,DTHDAA,DTHDE,EPS,TEND,THEND,TM
      INTEGER IFLAG
      LOGICAL OK
C     ..
C     .. Array Arguments ..
      REAL ER(3),Y(3)
C     ..
C     .. Scalars in Common ..
      REAL EIG,HPI,PI,TSAVEL,TSAVER,TWOPI
      INTEGER IND,ISAVE,T21
      LOGICAL PR
C     ..
C     .. Arrays in Common ..
      REAL TT(7,2),YY(7,3,2)
      INTEGER NT(2)
C     ..
C     .. Local Scalars ..
      REAL C,D,DPHIDE,DTHIN,DUM,EFF,FAC2,HU,HV,PHI,PHI0,
     +                 PHIDAA,PUP,PVP,PYPZ,PYPZ0,RHOSQ,S,T,TH,TH0,THIN,
     +                 THU,THU0,THV,THV0,TIN,TINTHZ,TMP,TOUT,U,V,XT,XT0,
     +                 YZ,YZ0
      INTEGER I,J,KFLAG,LFLAG,NK,NNN
      LOGICAL LOGIC
C     ..
C     .. Local Arrays ..
      REAL WORK(27),YP(3)
      INTEGER IWORK(5)
C     ..
C     .. External Subroutines ..
      EXTERNAL DXDT,FZ,GERK,INTEGZ,SETTHU,UV,UVPHI
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS,ATAN2,COS,EXP,LOG,SIGN,SIN
C     ..
C     .. Common blocks ..
      COMMON /DATAF/EIG,IND
      COMMON /PIE/PI,TWOPI,HPI
      COMMON /PRIN/PR,T21
      COMMON /TEMP/TT,YY,NT
      COMMON /TSAVE/TSAVEL,TSAVER,ISAVE
C     ..
      IF (.NOT.OK) THEN
          LOGIC = .FALSE.
      ELSE IF (TEND.LE.TM) THEN
          LOGIC = TEND .GE. TSAVEL
      ELSE
          LOGIC = TEND .LE. TSAVER
      END IF
C
      IF (LOGIC) THEN
          TINTHZ = TEND
          TH = THEND
          DTHIN = DTHDE
          Y(1) = TH
          Y(2) = DTHIN
          EFF = Y(3)
      ELSE
C           DO (INTEGRATE-FOR-PHI-OSC)
          T = TEND
C  THE LCO BOUNDARY CONDITION A1*[u,y] + A2*[v,y] = 0
C  IS EQUIVALENT TO COF1*sin(phi) - COF2*cos(phi) = 0.  SO
          PHI0 = ATAN2(COF2,COF1)

          IF (TM.GT.TEND) THEN
              J = 1
          ELSE
              J = 2
          END IF
C
C     WE WANT -PI/2 .LT. PHI0 .LE. PI/2.
C
          IF (COF1.LT.0.0D0) PHI0 = PHI0 - SIGN(PI,COF2)
          Y(1) = PHI0
          Y(2) = 0.0D0
          Y(3) = 0.0D0
C  PHIDAA REPRESENTS d(phi)/da  EVALUATED AT THE ENDPOINT (A or B).
          CALL DXDT(T,TMP,XT0)
          CALL FZ(T,Y,YP)
          PHIDAA = -YP(1)
          CALL UV(XT0,U,PUP,V,PVP,HU,HV)
          D = U*PVP - V*PUP
C
C  DETERMINE A DEFINITE VALUE FOR ATAN(u/pu') AND ATAN(v/pv') AT TEND:
C     SET THU0 AND THV0.
C
          THU0 = ATAN2(U,PUP)
          IF (U.LT.0.0D0) THU0 = THU0 + TWOPI
          CALL SETTHU(XT0,THU0)
          THV0 = ATAN2(V,PVP)
          IF (V.LT.0.0D0) THV0 = THV0 + TWOPI
   10     CONTINUE
          IF (THV0.LT.THU0) THEN
              THV0 = THV0 + TWOPI
              GO TO 10
          END IF
C
C   USING THE RELATION BETWEEN theta AND phi, NAMELY
C  tan(theta) = (u*cos(phi) + v*sin(phi))/(pup*cos(phi)+pv'*sin(phi))
C  AND, DIFFERENTIATING W.R.T. lambda or a or b,
C  d(theta)*sec(theta)**2 = -[u,v]*d(phi)/(pu'*cos(phi)+pv'*sin(phi))**2

C     SET TH0 AND COPY INTO THEND, OVERWRITING ITS INPUT VALUE.
C     ALSO, REDEFINE DTHDAA, OVERWRITING ITS INPUT VALUE.
C
          CALL UVPHI(U,PUP,V,PVP,THU0,THV0,PHI0,TH0)
          THEND = TH0
          C = COS(PHI0)
          S = SIN(PHI0)
          YZ0 = U*C + V*S
          PYPZ0 = PUP*C + PVP*S
          DUM = ABS(COS(TH0))

          IF (DUM.GE.0.5D0) THEN
              FAC2 = -D* (DUM/PYPZ0)**2
          ELSE
              FAC2 = -D* (SIN(TH0)/YZ0)**2
          END IF

C  DTHDAA REPRESENTS d(theta)/da  EVALUATED AT THE ENDPOINT (A or B).
          DTHDAA = PHIDAA*FAC2
          TOUT = TM
          I = 2
          TT(I,J) = T
          YY(I,1,J) = Y(1)
          YY(I,2,J) = Y(2)
          YY(I,3,J) = Y(3)

          DO 23 I = 1,27
             WORK(I) = 0.0D0
             IF(I.LE.5) IWORK(I) = 0
  23      CONTINUE

   18     CONTINUE
          NNN = 0
          KFLAG = -1
C
C     TSAVEL, TSAVER PRESUMED SET BY AN EARLIER CALL WITH EIG .NE. 0.
C
          IF (EIG.EQ.0.0D0 .AND. ISAVE.EQ.1) THEN
              IF (T.LE.TSAVEL) TOUT = TSAVEL
              IF (T.GT.TSAVER) TOUT = TSAVER
              KFLAG = 1
          END IF

          NK = 0
   20     CONTINUE
          CALL GERK(FZ,3,Y,T,TOUT,EPS,EPS,KFLAG,ER,WORK,IWORK)
          IF (KFLAG.GT.3) THEN
              IF (PR) WRITE (T21,FMT=*) ' KFLAG1 = 5 '
              IFLAG = 5
              RETURN
          END IF
          IF (KFLAG.EQ.3) THEN
              NK = NK + 1
              IF (NK.GE.10) THEN
                  KFLAG = 5
                  RETURN
              END IF
              GO TO 20
          END IF
          NNN = NNN + 1
C  PREVENT EXP(Y(3)) FROM BECOMING UNNECESSARILY SMALL.
          IF (Y(3).LT.-15.0D0) Y(3) = -15.0D0
          PHI = Y(1)
C
C     STORE UP TO SEVEN VALUES OF (T,PHI) FOR LATER REFERENCE.
C
          I = I + 1
          IF (I.LE.7) THEN
              TT(I,J) = T
              YY(I,1,J) = PHI
              YY(I,2,J) = Y(2)
              YY(I,3,J) = Y(3)
          END IF
C  STOP THE INTEGRATION FROM GOING TOO FAR AWAY FROM THE ENDPOINT.
          IF (10.0D0*ABS(PHI-PHI0).GE.PI) GO TO 30
          IF (KFLAG.EQ.-2) THEN
              IF (NNN.GT.500) THEN
                  EPS = 10.D0*EPS
                  GO TO 18
              END IF
              GO TO 20
          END IF
   30     CONTINUE
          IF (T.LE.TOUT) THEN
              TSAVEL = T
          ELSE
              TSAVER = T
          END IF
          NT(J) = I
          DPHIDE = Y(2)
          TINTHZ = T
          CALL DXDT(T,TMP,XT)
          CALL UV(XT,U,PUP,V,PVP,HU,HV)
          D = U*PVP - V*PUP
C
C     SET THU AND THV.
C
          THU = ATAN2(U,PUP)
          IF (U.LT.0.0D0) THU = THU + TWOPI
          CALL SETTHU(XT,THU)
          THV = ATAN2(V,PVP)
          IF (V.LT.0.0D0) THV = THV + TWOPI
   40     CONTINUE
          IF (THV.LT.THU) THEN
              THV = THV + TWOPI
              GO TO 40
          END IF
C
C  TRANSLATE FROM THE VARIABLES PHI, SIGMA, ETC. TO THETA, RHO,ETC.
C    DEFINE TH IN TERMS OF PHI, THU, AND THV.
C
          CALL UVPHI(U,PUP,V,PVP,THU,THV,PHI,TH)
          YZ = U*COS(PHI) + V*SIN(PHI)
          PYPZ = PUP*COS(PHI) + PVP*SIN(PHI)
          Y(1) = TH
          S = SIN(TH)
          C = COS(TH)
          DUM = ABS(C)

          IF (DUM.GE.0.5D0) THEN
              FAC2 = -D* (DUM/PYPZ)**2
          ELSE
              FAC2 = -D* (S/YZ)**2
          END IF

C  d(theta)/d(lambda) = FAC2*d(phi)/d(lambda)
          DTHIN = FAC2*DPHIDE
C  ALSO CONVERT THE GLOBAL ERROR ESTIMATES FROM THOSE FOR
C  PHI,ETC., TO THOSE CORRESPONDING TO THETA, ETC.
          ER(1) = FAC2*ER(1)
          ER(2) = FAC2*ER(2)
          ER(3) = FAC2*ER(3)
          Y(2) = DTHIN
          RHOSQ = EXP(2.0D0*Y(3))* (YZ**2+PYPZ**2)
          EFF = 0.5D0*LOG(RHOSQ* (S**2+C**2))
          Y(3) = EFF
C              END (INTEGRATE-FOR-PHI-OSC)
      END IF
      IF (TINTHZ.NE.TM) THEN
C  NOW INTEGRATE THE REST OF THE WAY TO TM FOR theta, etc.
          TIN = TINTHZ
          TOUT = TM
          THIN = TH
C           DO (INTEGRATE-FOR-THETA)
          CALL INTEGZ(TIN,TOUT,THIN,DTHIN,EPS,Y,LFLAG,ER,WORK,IWORK)
          Y(3) = EFF + Y(3)
C
      END IF
      RETURN
      END
C
      SUBROUTINE LP5(TEND,THEND,DTHDE,TM,EPS,ER,Y,IFLAG)
C  THIS PROGRAM IS CALLED FROM INTEG.
C
C  INPUT QUANTITIES:
C  OUTPUT QUANTITIES:
C
C    IT CONTROLS THE INTEGRATIONS FOR THE USUAL PRUEFER
C    VARIABLES THETA, RHO, ETC. FROM AN ENDPOINT THAT IS LP
C    BUT IS NOT AT +INF OR  -INF.  IN THIS CASE THE ENDPOINT IS
C    A REGULAR SINGULAR POINT AND THE INDICIAL EQUATION IS USED
C    TO OBTAIN THE INITIAL VALUE OF THETA AND IT'S SLOPE THERE
C    (WHICH IS COMPUTED BY THE SUBROUTINE FINELP.)


C     .. Scalar Arguments ..
      REAL DTHDE,EPS,TEND,THEND,TM
      INTEGER IFLAG
C     ..
C     .. Array Arguments ..
      REAL ER(3),Y(3)
C     ..
C     .. Scalars in Common ..
      REAL EIG
      INTEGER IND,T21
      LOGICAL PR
C     ..
C     .. Local Scalars ..
      REAL DDD,DTHIN,EFF,ONE,T,THIN,TIN,TMP,TOUT
      INTEGER I,KFLAG
C     ..
C     .. Local Arrays ..
      REAL WORK(27),YU(3)
      INTEGER IWORK(5)
C     ..
C     .. External Subroutines ..
      EXTERNAL F,GERK,INTEGZ
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS,MIN
C     ..
C     .. Common blocks ..
      COMMON /DATAF/EIG,IND
      COMMON /PRIN/PR,T21
C     ..
      ONE = 1.0D0
C
      IF (PR) WRITE (T21,FMT=*) ' THIS IS THE LP CASE AT '
      IF (PR) WRITE (T21,FMT=*) ' A FINITE ENDPOINT. '

C
      T = TEND
      YU(1) = THEND
      YU(2) = DTHDE
      YU(3) = 0.0D0

      DO 23 I = 1,27
         WORK(I) = 0.0D0
         IF(I.LE.5) IWORK(I) = 0
  23  CONTINUE

      TMP = 0.01D0
      DDD = MIN(TMP,ABS(TEND))
      TOUT = TEND + DDD* (-TEND)/ABS(TEND)
      KFLAG = 1
      CALL GERK(F,3,YU,T,TOUT,EPS,EPS,KFLAG,ER,WORK,IWORK)
C  NOW CONTINUE THE INTEGRATIONS TO TM.
      TIN = TOUT
      TOUT = TM
      THIN = YU(1)
      DTHIN = YU(2)
      EFF = YU(3)
C          DO (INTEGRATE-FOR-THETA)
      CALL INTEGZ(TIN,TOUT,THIN,DTHIN,EPS,Y,IFLAG,ER,WORK,IWORK)
      Y(3) = EFF + Y(3)
      RETURN
      END
C
      SUBROUTINE LP6(TEND,THEND,TM,DTHDE,EPS,Y,ER,IFLAG)
C  THIS PROGRAM IS CALLED FROM INTEG.
C
C  INPUT QUANTITIES: TEND,TM,THEND,DTHDE,EPS,Y
C  OUTPUT QUANTITIES: Y,ER,IFLAG
C
C    THIS CASE IS EITHER LP AT AN INFINITE ENDPOINT,
C    OR IS IRREGULAR LP AT A FINITE ENDPOINT,
C    OR IS DEFAULT BY THE USER .

C     .. Scalar Arguments ..
      REAL DTHDE,EPS,TEND,THEND,TM
      INTEGER IFLAG
C     ..
C     .. Array Arguments ..
      REAL ER(3),Y(3)
C     ..
C     .. Local Scalars ..
      REAL DTHIN,THIN,TIN,TOUT
C     ..
C     .. Local Arrays ..
      REAL WORK(27)
      INTEGER IWORK(5)
C     ..
C     .. External Subroutines ..
      EXTERNAL INTEGZ
C     ..
      TIN = TEND
      TOUT = TM
      THIN = THEND
      DTHIN = DTHDE
C         DO (INTEGRATE-FOR-THETA)
C        INITIALIZE ER BEFORE CALLING INTEGZ:
      ER(1) = 0.0D0
      ER(2) = 0.0D0
      ER(3) = 0.0D0
      CALL INTEGZ(TIN,TOUT,THIN,DTHIN,EPS,Y,IFLAG,ER,WORK,IWORK)
C            END (INTEGRATE-FOR-THETA)
C
      RETURN
      END
C
      SUBROUTINE INTEG(TEND,THEND,DTHDAA,DTHDE,TM,COEF1,COEF2,EPS,Y,ER,
     +                 OK,NC,IFLAG)
C     **********

C  THIS PROGRAM ACTUALLY DOES NOTHING BUT CALL THE
C    APPROPRIATE SUBROUTINE FROM THE LIST
C    REG,WREG,LCNO,LCO,LP5,LP6,
C    DEPENDING UPON THE VALUE OF NC.
C
C  INPUT QUANTITIES: TEND,THEND,DTHDAA,DTHDE,TM,COEF1,
C                    COEF2,EPS,Y,ER,OK,NC
C  OUTPUT QUANTITIES: DTHDAA,DTHDE,Y,ER,IFLAG
C
C     **********
C     .. Local Scalars ..
      REAL COF1,COF2
C     ..
C     .. External Subroutines ..
      EXTERNAL LCNO,LCO,LP5,LP6,REG,WREG
C     ..
C     NOTE: THE INPUT VALUES OF THEND AND DTHDAA ARE OVERWRITTEN
C           WHEN INTEGRATING FROM A LIMIT CIRCLE ENDPOINT.
C
C     .. Scalar Arguments ..
      REAL COEF1,COEF2,DTHDAA,DTHDE,EPS,TEND,THEND,TM
      INTEGER IFLAG,NC
      LOGICAL OK
C     ..
C     .. Array Arguments ..
      REAL ER(3),Y(3)
C     ..
      IFLAG = 1
      COF1 = COEF1
      COF2 = COEF2

      IF (NC.EQ.1) THEN

C        THIS IS THE REGULAR (NOT WEAKLY REGULAR) CASE.
          CALL REG(TEND,THEND,TM,DTHDE,EPS,Y,ER,IFLAG)
C
      ELSE IF (NC.EQ.2) THEN
C
C        THIS IS THE 'WEAKLY REGULAR' CASE.
          CALL WREG(TEND,THEND,DTHDE,TM,EPS,ER,Y,IFLAG)
C
      ELSE IF (NC.EQ.3) THEN

          CALL LCNO(TEND,THEND,DTHDAA,TM,COF1,COF2,EPS,Y,ER,IFLAG)

      ELSE IF (NC.EQ.4) THEN

          CALL LCO(TEND,THEND,DTHDAA,DTHDE,TM,COF1,COF2,EPS,Y,ER,OK,
     +             IFLAG)

      ELSE IF (NC.EQ.5) THEN
C
C        THIS IS THE CASE OF A REGULAR LP AT A FINITE END.
          CALL LP5(TEND,THEND,DTHDE,TM,EPS,ER,Y,IFLAG)
C
      ELSE IF (NC.EQ.6) THEN

C        THIS CASE IS EITHER LP AT AN INFINITE ENDPOINT,
C        OR ELSE IS IRREGULAR LP AT A FINITE ENDPOINT,
C        OR ELSE IS DEFAULT BY THE USER .
          CALL LP6(TEND,THEND,TM,DTHDE,EPS,Y,ER,IFLAG)

      END IF

      RETURN
      END
C
C
      SUBROUTINE INTEGZ(TIN,TOUT,THIN,DTHIN,EPS,Y,IFLAG,ER,WORK,IWORK)
C     **********
C
C  THIS PROGRAM IS USED FOR CALLING SUBROUTINE GERKZ WHEN THERE
C    IS NO SINGULARITY OF ANY KIND EITHER AT TIN OR AT TOUT.
C  IT IS CALLED FROM LCNO, REG, WREG, LCO, LP5, LP6.
C
C  INPUT QUANTITIES: TIN,TOUT,THIN,DTHIN,EPS,Y,ER
C  OUTPUT QUANTITIES: Y,IFLAG,ER
C
C     **********
C     .. Local Scalars ..
      INTEGER I,LFLAG,NK
C     ..
C     .. External Subroutines ..
      EXTERNAL F,GERKZ
C     ..
C     DO (INTEGRATE-FOR-TH)
C     .. Scalar Arguments ..
      REAL DTHIN,EPS,THIN,TIN,TOUT
      INTEGER IFLAG
C     ..
C     .. Array Arguments ..
      REAL ER(3),WORK(27),Y(3)
      INTEGER IWORK(5)
C     ..
C     .. Scalars in Common ..
      INTEGER T21
      LOGICAL PR
C     ..
C     .. Common blocks ..
      COMMON /PRIN/PR,T21
C     ..

C     IT IS ALWAYS ASSUMED THAT THE ARRAY ER HAS BEEN
C     INITIALIZED BEFORE THIS PROGRAM IS CALLED.
C     EITHER TO ZEROS, OR ELSE IT HAS THE VALUES IT ACQUIRED
C     THE LAST TIME GERK WAS CALLED.  THE PROGRAM INTEGZ
C     WILL BE DEPENDING UPON RECEIVING MEANINGFUL VALUES IN ER.

        DO 23 I = 1,27
           WORK(I) = 0.0D0
           IF(I.LE.5) IWORK(I) = 0
  23    CONTINUE

      Y(1) = THIN
      Y(2) = DTHIN
      Y(3) = 0.0D0
      LFLAG = 1
      NK = 0
   10 CONTINUE
      CALL GERKZ(F,3,Y,TIN,TOUT,EPS,EPS,LFLAG,ER,WORK,IWORK)
      IF (LFLAG.EQ.3) THEN
          NK = NK + 1
          IF (NK.GE.10) THEN
              LFLAG = 5
              RETURN
          END IF
          GO TO 10
      END IF
      IF (LFLAG.GT.3) THEN
          IF (PR) WRITE (T21,FMT=*) ' LFLAG = ',LFLAG
          IFLAG = 5
          RETURN
      END IF
C        END (INTEGRATE-FOR-TH)
      RETURN
      END

C
      SUBROUTINE FZERO(F,B,C,R,RE,AE,IFLAG)
C     **********
C  THIS PROGRAM IS CALLED FROM PERIO.
C  IT SEARCHES FOR A ZERO OF A FUNCTION F(X)
C  BETWEEN THE GIVEN VALUES B AND C UNTIL THE WIDTH OF
C  THE INTERVAL (B,C) HAS COLLAPSED TO WITHIN A TOLERANCE
C  SPECIFIED BY THE STOPPING CRITERION,
C     ABS(B-C) .LE. 2.*(RW*ABS(B)+AE).  THE METHOD USED IS
C  A COMBINATION OF BISECTION AND THE SECANT RULE.
C
C  THE MEANING OF IFLAG IS:-
C  IFLAG = 1, B IS WITHIN THE REQUESTED TOLERANCE OF A ZERO.
C             THE INTERVAL (B,C) COLLAPSED TO THE REQUESTED
C             TOLERANCE, THE FUNCTION CHANGES SIGN IN (B,C),
C             AND F(X) DECREASED IN MAGNITUDE AS (B,C) COLLAPSED.
C          2, F(B) = 0. HOWEVER, THE INTERVAL (B,C) MAY NOT HAVE
C             HAVE COLLAPSED TO THE REQUESTED TOLERANCE.
C          3, B MAY BE NEAR A SINGULAR POINT OF F(X).  THE
C             INTERVAL (B,C) COLLAPSED TO THE REQUESTED TOLERANCE
C             AND THE FUNCTION CHANGES SIGN IN (B,C), BUT F(X)
C             INCREASED IN MAGNITUDE AS (B,C) COLLAPSED.
C          4, NO CHANGE IN SIGN OF F(X) WAS FOUND ALTHOUGH THE
C             INTERVAL (B,C) COLLAPSED TO THE REQUESTED TOLERANCE.
C          5, TOO MANY (.GT. 500) FUNCTION EVALUATIONS USED.
C          6, NO MORE PROGRESS IS BEING MADE.
C
C     .. Local Scalars ..
      REAL A,ACBS,ACMB,AW,CMB,DIF,DIFS,FA,FB,FC,FX,FZ,P,Q,
     +                 RW,TOL,Z,ZER
      INTEGER IC,KOUNT
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS,MAX,MIN,SIGN
C     ..
C     .. Scalar Arguments ..
      REAL AE,B,C,R,RE
      INTEGER IFLAG
C     ..
C     .. Function Arguments ..
      REAL F
      EXTERNAL F
C     ..
      IFLAG = 1
      ZER = 0.0D0
      DIF = 1.D+9
      Z = R
      IF (R.LE.MIN(B,C) .OR. R.GE.MAX(B,C)) Z = C
      RW = MAX(RE,ZER)
      AW = MAX(AE,ZER)
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
      CMB = 0.5D0* (C-B)
      ACMB = ABS(CMB)
      TOL = RW*ABS(B) + AW
      IF (ACMB.LE.TOL) THEN
          IFLAG = 1
          IF (FB*FC.GE.0.0D0) IFLAG = 4
          IF (ABS(FB).GT.FX) IFLAG = 3
          RETURN
      END IF
      IF (FB.EQ.0.0D0) THEN
          IFLAG = 2
          RETURN
      END IF
      IF (KOUNT.GE.500) THEN
          IFLAG = 5
          RETURN
      END IF
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
          IF (8.0D0*ACMB.GE.ACBS) B = 0.5D0* (C+B)
      ELSE
          IC = 0
          ACBS = ACMB
          IF (P.LE.ABS(Q)*TOL) THEN
              B = B + SIGN(TOL,CMB)
          ELSE IF (P.GE.CMB*Q) THEN
              B = 0.5D0* (C+B)
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
C
      SUBROUTINE PERFUN(AA,BB,TMID,NCA,NCB,K11,K12,K21,K22,LAMBDA,ISLFN,
     +                  XT,PLOTF)
C     **********
C  THIS PROGRAM IS USED ONLY TO COMPUTE EIGENFUNCTIONS FOR PROBLEMS
C  WITH COUPLED BOUNDARY CONDITIONS.  IT IS CALLED BY SUBROUTINE
C  SLCOUP AND BY DRAW.
C
C  INPUT QUANTITIES: AA,BB,TMID,NCA,NCB,K11,K12,K21,K22,LAMBDA,
C                    ISLFN,XT
C  OUTPUT QUANTITIES: XT,PLOTF
C
C     .. Scalars in Common ..
C
      REAL EIGSAV,EPSMIN,HPI,PI,TSAVEL,TSAVER,TWOPI,Z
      INTEGER IND,ISAVE,T21
      LOGICAL PR
C     ..
C     .. Local Scalars ..
      REAL A1,A2,ADTH,AEE,B1,B2,BEE,BRA,BRB,C,DA,DB,DELTMP,
     +                 DTH,DTHDAA,DTHDBB,DTHDEA,DTHDEB,DTHM,DTHS,E,EFF,
     +                 EPS,ETAA,ETAB,F,HU,HV,PUP,PVP,PY1PL,PY1PR,PY2PL,
     +                 PY2PR,RHOL,RHOR,T,THA,THB,THL,THR,THS,THT,TM,TMP,
     +                 TMP0,TMP1,TMP2,TMPS,TT,U,V,XX,Y1L,Y1R,Y2L,Y2R,YM,
     +                 YM1,YM2,ZM
      INTEGER I,J,K,LFLAG,NC,NI,NM,NMID,NMIDM,NMIDP,NPI
      LOGICAL AOK,BOK,OK
C     ..
C     .. Local Arrays ..
      REAL ERR(3),PYLI(500),PYRI(500),THLI(500),THRI(500),
     +                 TI(1000),Y(3),YL(500,2,2),YLI(500),YR(500,2,2),
     +                 YRI(500)
C     ..
C     .. External Subroutines ..
      EXTERNAL DXDT,INTEG,UV
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS,ATAN2,COS,EXP,MAX,MIN,SIN,SQRT
C     ..
C     .. Common blocks ..
      COMMON /DATAF/EIGSAV,IND
      COMMON /PIE/PI,TWOPI,HPI
      COMMON /PRIN/PR,T21
      COMMON /RNDOFF/EPSMIN
      COMMON /TSAVE/TSAVEL,TSAVER,ISAVE
      COMMON /Z1/Z
C     ..
C
C     .. Scalar Arguments ..
      REAL AA,BB,K11,K12,K21,K22,LAMBDA,TMID
      INTEGER ISLFN,NCA,NCB
C     ..
C     .. Array Arguments ..
      REAL PLOTF(1000,6),XT(1000,2)
C     ..
      IND = 1
      EIGSAV = LAMBDA
      EPS = EPSMIN
      AOK = NCA .EQ. 1
      BOK = NCB .EQ. 1
C
C  JUST IN CASE THAT THE USER SUPPLIED BOUNDARY CONDITION
C  FUNCTIONS u,v AND/OR U,V HAVE NOT BEEN "NORMALIZED" TO HAVE
C   [u,v](a) = 1.0 AND/OR [U,V](b) = 1.0,
C  THE VALUES OF THE WRONSKIANS BRA = [u,v](a) & BRB = [U,V](b),
C  WILL BE NEEDED, AS ARE THE QUANTITIES ETAA, DA, ETAB,DB
C  WHERE BRA = DA*ETAA AND BRB = DB*ETAB, AND DA,DB ARE POSITIVE.
      BRA = 1.0D0
      BRB = 1.0D0
      ETAA = 1.0D0
      ETAB = 1.0D0
      DA = 1.0D0
      DB = 1.0D0
      IF (NCA.EQ.3 .OR. NCA.EQ.4) THEN
          TT = AA
          CALL DXDT(TT,TMP,XX)
          CALL UV(XX,U,PUP,V,PVP,HU,HV)
          BRA = U*PVP - V*PUP
          DA = SQRT(ABS(BRA))
          ETAA = BRA/ (DA*DA)
      END IF
      IF (NCB.EQ.3 .OR. NCB.EQ.4) THEN
          TT = BB
          CALL DXDT(TT,TMP,XX)
          CALL UV(XX,U,PUP,V,PVP,HU,HV)
          BRB = U*PVP - V*PUP
          DB = SQRT(ABS(BRB))
          ETAB = BRB/ (DB*DB)
      END IF
      IF (PR) WRITE (T21,FMT=*) ' BRA,BRB = ',BRA,BRB
C
C  GET THE ARRAY OF POINTS T IN (-1,1) WHERE THE EIGENFUNCTION
C    IS WANTED.  THESE POINTS CORRESPOND TO POINTS X IN (a,b).

      DO 3 I = 1,ISLFN
          TI(I) = XT(9+I,2)
          CALL DXDT(TI(I),TMP,XX)
          XT(9+I,1) = XX
    3 CONTINUE
      NI = ISLFN

C  FIND THE INDEX, NMID, IN (1,NI) CORRESPONDING NEARLY TO
C    TMID, OR XMID.

      NMID = 0
      DO 5 I = 1,NI
          IF (TI(I).LE.TMID) NMID = I
    5 CONTINUE
      NMIDP = NMID + 5
      NMIDM = NMID - 5
C
C  NOW COMPUTE VALUES OF THE FUNDAMENTAL SYSTEM OF FUNCTIONS
C    Y1L, Y1R, Y2L, Y2R
C  AT THE POINTS TI.
C
      Z = 1.0D0
C
C  N.B.: IN LIMIT CIRCLE CASE, IF Y IS A SOLUTION,
C
C         [Y,V] =  SIGMA*[U,V]*COS(PHI)
C         [Y,U] = -SIGMA*[U,V]*SIN(PHI)
C
C     FOR Y2: NEUMANN PROBLEM :
C
C  N.B. IF ENDPOINT A IS REGULAR, INTEG USES THA;
C       BUT IF A IS LC, INTEG USES A1,A2.
C
      A1 = 0.0D0
      A2 = -1.0D0
      THA = HPI
      Y(1) = THA
      Y(2) = 1.0D0
      Y(3) = 0.0D0
      DTHDAA = 0.0D0
      DTHDEA = 1.0D0
      NC = NCA
      OK = AOK
      T = AA
      EFF = 0.0D0
      LFLAG = 1
      ISAVE = 0
      DO 12 I = 1,NMIDP
          THT = Y(1)
          TM = TI(I)
          IF (TM.GT.-1.0D0) CALL INTEG(T,THT,DTHDAA,DTHDEA,TM,A1,A2,EPS,
     +                                 Y,ERR,OK,NC,LFLAG)
          IF (NC.EQ.4) THEN
              EFF = Y(3)
          ELSE IF (TM.GT.-1.0D0) THEN
              OK = .TRUE.
              NC = 1
              EFF = EFF + Y(3)
          END IF
          RHOL = EXP(EFF)
          THL = Y(1)
          Y2L = RHOL*SIN(THL)
          PY2PL = RHOL*COS(THL)
C
          Y2L = Y2L*DA
          PY2PL = PY2PL*DA
          YL(I,2,1) = Y2L
          YL(I,2,2) = PY2PL
          T = TM
          IF (TM.GT.-1.0D0) THEN
              OK = .TRUE.
              NC = 1
          END IF
          IF (T.LT.-0.9D0 .AND. NCA.EQ.4) THEN
C  IN THIS CASE, INTEGRATE FROM AA AGAIN, AS AN OSC PROBLEM.
              NC = 4
              OK = .FALSE.
              T = AA
              Y(1) = THA
              Y(2) = 0.0D0
              Y(3) = 0.0D0
          END IF
   12 CONTINUE
C
C  N.B. IF ENDPOINT B IS REGULAR, INTEG USES THB;
C       BUT IF B IS LC, INTEG USES B1,B2.
      B1 = 0.0D0
      B2 = -1.0D0
      THB = HPI
      Y(1) = THB
      Y(2) = 1.0D0
      Y(3) = 0.0D0
      DTHDBB = 0.0D0
      DTHDEB = 1.0D0
      NC = NCB
      OK = BOK
      T = BB
      EFF = 0.0D0
      LFLAG = 1
      ISAVE = 0
      DO 22 I = NI,NMIDM,-1
          THT = Y(1)
          TM = TI(I)
          IF (TM.LT.1.0D0) CALL INTEG(T,THT,DTHDBB,DTHDEB,TM,B1,B2,EPS,
     +                                Y,ERR,OK,NC,LFLAG)
          IF (NC.EQ.4) THEN
              EFF = Y(3)
          ELSE IF (TM.LT.1.0D0) THEN
              OK = .TRUE.
              NC = 1
              EFF = EFF + Y(3)
          END IF
          RHOR = EXP(EFF)
          THR = Y(1)
          Y2R = RHOR*SIN(THR)
          PY2PR = RHOR*COS(THR)
C
          Y2R = Y2R*DB
          PY2PR = PY2PR*DB
          YR(I,2,1) = Y2R
          YR(I,2,2) = PY2PR
          T = TM
          IF (T.LT.1.0D0) THEN
              OK = .TRUE.
              NC = 1
          END IF
          IF (T.GT.0.9D0 .AND. NCB.EQ.4) THEN
C  IN THIS CASE, INTEGRATE FROM BB AGAIN, AS AN OSC PROBLEM.
              NC = 4
              OK = .FALSE.
              T = BB
              Y(1) = THB
              Y(2) = 0.0D0
              Y(3) = 0.0D0
          END IF
   22 CONTINUE
C
C     FOR Y1: DIRICHLET PROBLEM :
C
C  N.B. IF ENDPOINT A IS REGULAR, INTEG USES THA;
C       BUT IF A IS LC, INTEG USES A1,A2.
      A1 = 1.0D0
      A2 = 0.0D0
      THA = 0.0D0
      Y(1) = THA
      Y(2) = 1.0D0
      Y(3) = 0.0D0
      DTHDAA = 0.0D0
      DTHDEA = 1.0D0
      NC = NCA
      OK = AOK
      T = AA
      EFF = 0.0D0
      LFLAG = 1
      ISAVE = 0
      DO 32 I = 1,NMIDP
          THT = Y(1)
          TM = TI(I)
          IF (TM.GT.-1.0D0) CALL INTEG(T,THT,DTHDAA,DTHDEA,TM,A1,A2,EPS,
     +                                 Y,ERR,OK,NC,LFLAG)
          IF (NC.EQ.4) THEN
              EFF = Y(3)
          ELSE IF (TM.GT.-1.0D0) THEN
              NC = 1
              OK = .TRUE.
              EFF = EFF + Y(3)
          END IF
          RHOL = EXP(EFF)
          THL = Y(1)
          Y1L = RHOL*SIN(THL)
          PY1PL = RHOL*COS(THL)
C
          Y1L = Y1L*DA*ETAA
          PY1PL = PY1PL*DA*ETAA
          YL(I,1,1) = Y1L
          YL(I,1,2) = PY1PL
          T = TM
          IF (T.GT.-1.0D0) THEN
              OK = .TRUE.
              NC = 1
          END IF
          IF (T.LT.-0.9D0 .AND. NCA.EQ.4) THEN
C  IN THIS CASE, INTEGRATE FROM AA AGAIN, AS AN OSC PROBLEM.
              NC = 4
              OK = .FALSE.
              T = AA
              Y(1) = THA
              Y(2) = 0.0D0
              Y(3) = 0.0D0
          END IF
   32 CONTINUE
C
C  N.B. IF ENDPOINT B IS REGULAR, INTEG USES THB;
C       BUT IF B IS LC, INTEG USES B1,B2.
      B1 = 1.0D0
      B2 = 0.0D0
      THB = 0.0D0
      Y(1) = THB
      Y(2) = 1.0D0
      Y(3) = 0.0D0
      DTHDBB = 0.0D0
      DTHDEB = 1.0D0
      NC = NCB
      OK = BOK
      T = BB
      EFF = 0.0D0
      LFLAG = 1
      ISAVE = 0
      DO 42 I = NI,NMIDM,-1
          THT = Y(1)
          TM = TI(I)
          IF (TM.LT.1.0D0) CALL INTEG(T,THT,DTHDBB,DTHDEB,TM,B1,B2,EPS,
     +                                Y,ERR,OK,NC,LFLAG)
          IF (NC.EQ.4) THEN
              EFF = Y(3)
          ELSE IF (TM.LT.1.0D0) THEN
              NC = 1
              OK = .TRUE.
              EFF = EFF + Y(3)
          END IF
          RHOR = EXP(EFF)
          THR = Y(1)
          Y1R = RHOR*SIN(THR)
          PY1PR = RHOR*COS(THR)
          IF (NCB.EQ.3) THEN
              Y1R = -Y1R
              PY1PR = -PY1PR
          END IF
C
          Y1R = Y1R*DB*ETAB
          PY1PR = PY1PR*DB*ETAB
          YR(I,1,1) = Y1R
          YR(I,1,2) = PY1PR
          T = TM
          IF (T.LT.1.0D0) THEN
              OK = .TRUE.
              NC = 1
          END IF
          IF (T.GT.0.9D0 .AND. NCB.EQ.4) THEN
C  IN THIS CASE, INTEGRATE FROM BB AGAIN, AS AN OSC PROBLEM.
              NC = 4
              OK = .FALSE.
              T = BB
              Y(1) = THB
              Y(2) = 0.0D0
              Y(3) = 0.0D0
          END IF
   42 CONTINUE
C
C  NOW WE WANT TO FIND CONSTANTS, AEE, BEE, AND
C  E, F SO THAT THE EIGENFUNCTION , Y, CAN BE OBTAINED AS
C      Y = AEE*Y1L + BEE*Y2L = E*Y1R + F*Y2R.
C  THE CONSTANTS E, F ARE OBTAINED IN TERMS OF AEE, BEE
C  BY MEANS OF THE COUPLED BOUNDARY CONDITIONS, AND ARE
C      E = K21*BEE + K22*AEE
C      F = K11*BEE + K12*AEE,
C
C   WHILE THE AEE, BEE ARE GOING TO BE CHOSEN BELOW BY THE
C   REQUIREMENT THAT Y AND PY' ARE CONTINUOUS.
C
C
      NM = NMID + 5

      DTHM = 5.0D0
      DELTMP = PI/100.0D0
      TMP0 = 0.0D0
      K = 0
  150 CONTINUE
      DO 110 J = 1,100
          TMP0 = TMP0 + DELTMP
          AEE = SIN(TMP0)
          BEE = COS(TMP0)
          E = K21*BEE + K22*AEE
          F = K11*BEE + K12*AEE
          TMP1 = AEE*YL(NM,1,1) + BEE*YL(NM,2,1)
          TMP2 = AEE*YL(NM,1,2) + BEE*YL(NM,2,2)
          TMP = ATAN2(TMP1,TMP2)
          IF (TMP.LT.0.0D0) TMP = TMP + PI
          THL = TMP
          TMP1 = E*YR(NM,1,1) + F*YR(NM,2,1)
          TMP2 = E*YR(NM,1,2) + F*YR(NM,2,2)
          TMP = ATAN2(TMP1,TMP2)
          IF (TMP.LE.0.0D0) TMP = TMP + PI
          THR = TMP
          DTH = THR - THL
          ADTH = ABS(DTH)
          IF (ADTH.LT.DTHM) THEN
              DTHM = ADTH
              TMPS = TMP0
              DTHS = DTH
          END IF
  110 CONTINUE
      K = K + 1
      IF (K.LT.2 .OR. (K.LE.4.AND.ABS(DTHS).GT. (1.0D-6))) THEN
          TMP0 = TMPS - DELTMP
          DELTMP = DELTMP/50.0D0
          GO TO 150
      END IF
C  WE NOW HAVE THE "BEST" CHOICE OF THE RATIO OF AEE:BEE .
C
C  NOW RESET NM WHERE YM IS THE LARGEST VALUE OF Y.
      NM = NMIDM
      YM = 0.0D0
      DO 56 I = NMIDM,NMIDP
          TMP1 = AEE*YL(I,1,1) + BEE*YL(I,2,1)
          TMP = ABS(TMP1)
          TMP2 = AEE*YL(I,1,2) + BEE*YL(I,2,2)
          TMP2 = ATAN2(TMP1,TMP2)
          IF (TMP2.LT.0.0D0) TMP2 = TMP2 + PI
          IF (TMP.GT.YM) THEN
              YM = TMP
              NM = I
          END IF
   56 CONTINUE
C
C  DEFINE THLI(*) AND SET YM1 WHERE YL IS THE GREATEST.
      YM1 = 0.0D0
      DO 62 I = 1,NMIDP
          YLI(I) = (AEE*YL(I,1,1)+BEE*YL(I,2,1))
          PYLI(I) = (AEE*YL(I,1,2)+BEE*YL(I,2,2))
          THLI(I) = ATAN2(YLI(I),PYLI(I))
          IF (THLI(I).LT.0.0D0) THLI(I) = THLI(I) + PI
          TMP = ABS(YLI(I))
          IF (TMP.GT.YM1) YM1 = TMP
   62 CONTINUE

C  DEFINE THRI(*) AND SET YM2 WHERE YR IS THE GREATEST.
      YM2 = 0.0D0
      DO 72 I = NMIDM,NI
          YRI(I) = (E*YR(I,1,1)+F*YR(I,2,1))
          PYRI(I) = (E*YR(I,1,2)+F*YR(I,2,2))
          THRI(I) = ATAN2(YRI(I),PYRI(I))
          IF (THRI(I).LT.0.0D0) THRI(I) = THRI(I) + PI
          TMP = ABS(YRI(I))
          IF (TMP.GT.YM2) YM2 = TMP
   72 CONTINUE
C
C  ADJUST THE THLI(*) TO BE CONTINUOUS; I.E. HAVE NO "JUMPS".
      NPI = 0
      THS = THLI(1)
      DO 82 I = 1,NMIDP
          TMP = THLI(I) - THS
          THS = THLI(I)
          IF (TMP.LT.-2.0D0) THEN
              NPI = NPI + 1
          END IF
          THLI(I) = THLI(I) + NPI*PI
   82 CONTINUE
C
C  ADJUST THE THRI(*) TO BE CONTINUOUS; I.E. HAVE NO "JUMPS".
      NPI = 0
      THS = THRI(NMIDM)
      DO 87 I = NMIDM,NI
          TMP = THRI(I) - THS
          THS = THRI(I)
          IF (TMP.LT.-2.0D0) THEN
              NPI = NPI + 1
          END IF
          THRI(I) = THRI(I) + NPI*PI
   87 CONTINUE
C
C  FIND J IN (NMIDM,NMIDP) WHERE TMP, BELOW, IS GREATEST.
      J = NMIDM
      ZM = 0.0D0
      DO 177 I = NMIDM,NMIDP
          TMP = MIN(ABS(YLI(I)),ABS(YRI(I)))
          IF (TMP.GT.ZM) THEN
              ZM = TMP
              J = I
          END IF
  177 CONTINUE
C
      C = YLI(J)/YRI(J)
      YM = MAX(YM1,C*YM2)
C  WITH THIS C AND YM RESCALE THE YLI, PYLI, YRI, PYRI
C  SO AS TO MAKE YL AND YR "CONTINUOUS".
C  ALSO FILL THE ARRAY PLOTF.

      DO 182 I = 1,NM
          YLI(I) = YLI(I)/YM
          PYLI(I) = PYLI(I)/YM
          PLOTF(9+I,1) = YLI(I)
          PLOTF(9+I,2) = PYLI(I)
          PLOTF(9+I,3) = YLI(I)
          PLOTF(9+I,4) = PYLI(I)
          PLOTF(9+I,5) = THLI(I)
          PLOTF(9+I,6) = SQRT(YLI(I)**2+PYLI(I)**2)
  182 CONTINUE
      DO 192 I = NM,NI
          YRI(I) = C*YRI(I)/YM
          PYRI(I) = C*PYRI(I)/YM
          PLOTF(9+I,1) = YRI(I)
          PLOTF(9+I,2) = PYRI(I)
          PLOTF(9+I,3) = YRI(I)
          PLOTF(9+I,4) = PYRI(I)
          PLOTF(9+I,5) = THRI(I)
          PLOTF(9+I,6) = SQRT(YRI(I)**2+PYRI(I)**2)
  192 CONTINUE
C
      RETURN
      END
C
      REAL FUNCTION FF(ALFLAM)
C     **********
C  THIS PROGRAM COMPUTES THE FUNCTION WHOSE ZEROS ARE THE
C  EIGENVALUES OF PROBLEMS WITH COUPLED BOUNDARY CONDITIONS.
C
C  INPUT QUANTITIES: ALFLAM,AA,BB,TMID,EIGSAV,IND,AOK,BOK,
C                    NCA,NCB,ALFA,K11,K12,K21,K22,ETAA,ETAB,
C                    DA,DB,EPSMIN
C  OUTPUT QUANTITIES: FF
C
C       THE STURM-LIOUVILLE DIFFERENTIAL EQUATION IS
C           -(py')' + q*y = lambda*w*y ,
C  AND THE COUPLED BOUNDARY CONDITIONS ARE OF THE FORM
C            y(b) = exp(-i*alfa)*(k11*y(a) + k12*py'(a))
C           py'(b) = exp(-i*alfa)*(k21*y(a) + k22*py'(a))
C  (IN THE REGULAR CASE).
C
C  WHEN y1,y2 IS A SYSTEM OF SOLUTIONS SUITABLY DEFINED AT A,
C  AND Y1,Y2 IS A SYSTEM OF SOLUTIONS SUITABLY DEFINED AT B,
C  THEN
C  D(LAMBDA) IS DEFINED BY THE FORMULA:
C    D(LAMBDA) := K11[Y2,y1]+K12[y2,Y2)+K21[Y1,y1]+K22[y2,Y1]
C  AND THE EIGENVALUES OF THE COUPLED BOUNDARY CONDITION
C  PROBLEM ARE THE ZEROS OF
C        D(LAMBDA) - 2.0*COS(ALFA)
C
C  THIS PROGRAM RECEIVES DATA FROM SUBROUTINE PERIO VIA
C  COMMON/PASS2, COMMON/ENCEE, AND COMMON/ABEE,
C  AND RETURNS DATA TO PERIO VIA COMMON/TERM
C
C     **********
C     .. Scalars in Common ..
C
      REAL AA,ALFA,BB,BB1,BB2,DA,DB,DMU,DNU,EIGSAV,EPSMIN,
     +                 ETAA,ETAB,FFMAX,FFMIN,FFV,HPI,K11,K12,K21,K22,PI,
     +                 TMID,TRM11,TRM12,TRM21,TRM22,TSAVEL,TSAVER,TWOPI,
     +                 Z
      INTEGER IND,ISAVE,NCA,NCB,T21
      LOGICAL AOK,BOK,PR
C     ..
C     .. Local Scalars ..
      REAL A1,A2,B1,B2,DTHDAA,DTHDBB,DTHDEA,DTHDEB,EPS,
     +                 LAMBDA,PY1PL,PY1PR,PY2PL,PY2PR,RHOL,RHOR,THA,THB,
     +                 THL,THR,Y1L,Y1R,Y2L,Y2R
      INTEGER LFLAG
C     ..
C     .. Local Arrays ..
      REAL ERR(3),Y(3)
C     ..
C     .. External Subroutines ..
      EXTERNAL INTEG
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC COS,EXP,MAX,MIN,SIN
C     ..
C     .. Common blocks ..
      COMMON /ABEE/AA,BB,TMID
      COMMON /DATAF/EIGSAV,IND
      COMMON /ENCEE/AOK,BOK,NCA,NCB
      COMMON /PASS2/ALFA,K11,K12,K21,K22,ETAA,ETAB,DA,DB
      COMMON /PIE/PI,TWOPI,HPI
      COMMON /PRIN/PR,T21
      COMMON /RNDOFF/EPSMIN
      COMMON /TERM/TRM11,TRM12,TRM21,TRM22,FFMIN,FFMAX,DMU,DNU,BB1,BB2,
     +       FFV
      COMMON /TSAVE/TSAVEL,TSAVER,ISAVE
      COMMON /Z1/Z
C     ..
C
C     .. Scalar Arguments ..
      REAL ALFLAM
C     ..
      IND = 1
      LAMBDA = ALFLAM
      EIGSAV = LAMBDA
      EPS = EPSMIN
C
   50 CONTINUE
      Z = 1.0D0

C  LET y1,y2 DENOTE SOLUTIONS OF
C      -(py')' + q*y = lambda*w*y
C  SATISFYING THE BOUNDARY CONDITIONS AT a
C       y1(a) = 0,  py1'(a) = 1
C       y2(a) = 1,  py2'(a) = 0
C  IF a IS REGULAR,
C  OR
C       [y1,u](a) = 0, [y1,v](a) = 1
C       [y2,u](a) = 1, [y2,v](a) = 0
C  IF a IS LC.
C  DENOTE THEIR VALUES AT THE MIDPOINT BY
C    Y1L, PY1PL AND Y2L, PY2PL, RESPECTIVELY.

C  SIMILARLY, LET Y1,Y2 DENOTE SOLUTIONS OF
C      -(py')' + q*y = lambda*w*y
C  SATISFYING THE BOUNDARY CONDITIONS AT b
C       Y1(b) = 0, pY1'(b) = 1
C       Y2(b) = 1, pY2'(b) = 0
C  IF b IS REGULAR,
C  OR
C       [Y1,U](b) = 0, [Y1,V](b) = 1
C       [Y2,U](b) = 1, [Y2,V](b) = 0
C  IF b IS LC.
C  DENOTE THEIR VALUES AT THE MIDPOINT BY
C    Y1R, PY1PR AND Y2R, PY2PR, RESPECTIVELY.
C
C  THUS:
C     FOR y2 & Y2: NEUMANN PROBLEM :
C
C  N.B. IF ENDPOINT A IS REGULAR, INTEG USES THA;
C       BUT IF A IS LC, INTEG USES A1,A2.
C
      A1 = 0.0D0
      A2 = -1.0D0
      THA = HPI
      Y(1) = THA
      Y(2) = 1.0D0
      Y(3) = 0.0D0
      DTHDAA = 0.0D0
      DTHDEA = 1.0D0
      LFLAG = 1
      ISAVE = 0
      CALL INTEG(AA,THA,DTHDAA,DTHDEA,TMID,A1,A2,EPS,Y,ERR,AOK,NCA,
     +           LFLAG)
      IF (LFLAG.EQ.5) THEN
          EPS = 10.0D0*EPS
          GO TO 50
      END IF
      RHOL = EXP(Y(3))
      THL = Y(1)
      Y2L = RHOL*SIN(THL)
      PY2PL = RHOL*COS(THL)
C
C  IF [u,v](a) .NE. 1, LET [u,v](a) = DA*ETAA (ABS(ETAA) = 1).
C  THEN NORMALIZE:
      Y2L = Y2L*DA
      PY2PL = PY2PL*DA
C
C  N.B. IF ENDPOINT B IS REGULAR, INTEG USES THB;
C       BUT IF B IS LC, INTEG USES B1,B2.
      B1 = 0.0D0
      B2 = -1.0D0
      THB = HPI
      Y(1) = THB
      Y(2) = 1.0D0
      Y(3) = 0.0D0
      DTHDBB = 0.0D0
      DTHDEB = 1.0D0
      LFLAG = 1
      ISAVE = 0
      CALL INTEG(BB,THB,DTHDBB,DTHDEB,TMID,B1,B2,EPS,Y,ERR,BOK,NCB,
     +           LFLAG)
      IF (LFLAG.EQ.5) THEN
          EPS = 10.0D0*EPS
          GO TO 50
      END IF
      THR = Y(1)
      RHOR = EXP(Y(3))
      Y2R = RHOR*SIN(THR)
      PY2PR = RHOR*COS(THR)
C
C  IF [U,V](b) .NE. 1, LET [U,V](b) = DB*ETAB (ABS(ETAB) = 1).
C  THEN NORMALIZE:
      Y2R = Y2R*DB
      PY2PR = PY2PR*DB

C
C     FOR y1 & Y1: DIRICHLET PROBLEM :
C
C  N.B. IF ENDPOINT A IS REGULAR, INTEG USES THA;
C       BUT IF A IS LC, INTEG USES A1,A2.
      A1 = 1.0D0
      A2 = 0.0D0
      THA = 0.0D0
      Y(1) = THA
      Y(2) = 1.0D0
      Y(3) = 0.0D0
      DTHDAA = 0.0D0
      DTHDEA = 1.0D0
      LFLAG = 1
      ISAVE = 0
      CALL INTEG(AA,THA,DTHDAA,DTHDEA,TMID,A1,A2,EPS,Y,ERR,AOK,NCA,
     +           LFLAG)
      IF (LFLAG.EQ.5) THEN
          EPS = 10.0D0*EPS
          GO TO 50
      END IF
      RHOL = EXP(Y(3))
      THL = Y(1)
      Y1L = RHOL*SIN(THL)
      PY1PL = RHOL*COS(THL)
C
C  IF [u,v](a) .NE. 1, LET [u,v](a) = DA*ETAA (ABS(ETAA) = 1).
C  THEN NORMALIZE:
      Y1L = Y1L*DA*ETAA
      PY1PL = PY1PL*DA*ETAA
C
C  N.B. IF ENDPOINT B IS REGULAR, INTEG USES THB;
C       BUT IF B IS LC, INTEG USES B1,B2.
      B1 = 1.0D0
      B2 = 0.0D0
      THB = 0.0D0
      Y(1) = THB
      Y(2) = 1.0D0
      Y(3) = 0.0D0
      DTHDBB = 0.0D0
      DTHDEB = 1.0D0
      LFLAG = 1
      ISAVE = 0
      CALL INTEG(BB,THB,DTHDBB,DTHDEB,TMID,B1,B2,EPS,Y,ERR,BOK,NCB,
     +           LFLAG)
      IF (LFLAG.EQ.5) THEN
          EPS = 10.0D0*EPS
          GO TO 50
      END IF
      RHOR = EXP(Y(3))
      THR = Y(1)
      Y1R = RHOR*SIN(THR)
      PY1PR = RHOR*COS(THR)
      IF (NCB.EQ.2 .OR. NCB.EQ.3) THEN
          Y1R = -Y1R
          PY1PR = -PY1PR
      END IF
C
C  IF [U,V](b) .NE. 1, LET [U,V](b) = DB*ETAB (ABS(ETAB) = 1).
C  THEN NORMALIZE:
      Y1R = Y1R*DB*ETAB
      PY1PR = PY1PR*DB*ETAB

C  N.B. D(LAMBDA) IS DEFINED BY THE FORMULA:
C    D(LAMBDA) := K11[Y2,y1]+K12[y2,Y2)+K21[Y1,y1]+K22[y2,Y1]
C  AND THE EIGENVALUES OF THE COUPLED BOUNDARY CONDITION
C  PROBLEM ARE THE ZEROS OF
C        D(LAMBDA) - 2.0*COS(ALFA)
C
      TRM11 = Y2R*PY1PL - Y1L*PY2PR
      TRM22 = Y2L*PY1PR - Y1R*PY2PL
      TRM21 = Y1R*PY1PL - Y1L*PY1PR
      TRM12 = Y2L*PY2PR - Y2R*PY2PL
C
      BB1 = K12*TRM12 + K22*TRM22
      BB2 = K21*TRM21 + K11*TRM11
C
      FF = BB1 + BB2 - 2.0D0*COS(ALFA)
c     IF(PR)WRITE(T21,*) ' TRMS =',TRM11,TRM12,TRM21,TRM22
c     IF(PR)WRITE(T21,*) ' Y1L,Y1R = ',Y1L,Y1R
c     IF(PR)WRITE(T21,*) ' PY1PL,PY1PR = ',PY1PL,PY1PR
c     IF(PR)WRITE(T21,*) ' Y2L,Y2R = ',Y2L,Y2R
c     IF(PR)WRITE(T21,*) ' PY2PL,PY2PR = ',PY2PL,PY2PR
C     IF(PR)WRITE(T21,*) ' BB1,BB2 = ',BB1,BB2
c      dmu = k22*trm21 + k12*trm11
c      dnu = k11*trm12 + k21*trm22
c      dmu = dmu/sqrt(1.0+dmu**2)
c      dnu = dnu/sqrt(1.0+dnu**2)

      FFV = FF
      IF (PR) WRITE (T21,FMT=*) ' EIG,FF = ',EIGSAV,FF
      FFMIN = MIN(FFMIN,FF)
      FFMAX = MAX(FFMAX,FF)

      RETURN
      END
C
      SUBROUTINE SLCOUP(A,B,INTAB,P0ATA,QFATA,P0ATB,QFATB,A1,A2,B1,B2,
     +                  NUMEIG,EIG,TOL,IFLAG,ICPFUN,CPFUN,NCA,NCB,ALFA,
     +                  K11,K12,K21,K22)
C
C     This routine is for the sole purpose of making it as simple as
C     possible to obtain the eigenvalues and eigenfunctions of a
C     Sturm-Liouville problem with coupled boundary conditions.
C
C  INPUT QUANTITIES: A,B,INTAB,P0ATA,QFATA,P0ATB,QFATB,A1,A2,B1,B2,
C                    NUMEIG,EIG,TOL,ICPFUN,CPFUN,NCA,NCB,ALFA,
C                    K11,K12,K21,K22
C  OUTPUT QUANTITIES: EIG,TOL,CPFUN
C
C     The differential equation is of the form
C
C       -(p(x)*y'(x))' + q(x)*y(x) = eig*w(x)*y(x)  on (a,b)
C
C     with user-supplied coefficient functions p, q, and w,
C     and the coupled boundary conditions are of the form
C
C        Yb = exp(-i*alfa)*Ya,
C
C     where
C              (  y(b)
C        Yb =  (             if b is Regular,
C              (py'(b)
C
C              (  [y,U](b)
C        Yb =  (             if b is Limit Circle.
C              (  [y,V](b)
C
C              (  k11*y(a) + k12*py'(a)
C        Ya =  (                         if a is Regular,
C              (  k21*y(a) + k22*py'(a)
C
C              (  k11*[y,u](a) + k12*[y,v](a)
C        Ya =  (                          if a is Limit Circle,
C              (  k21*[y,u](a) + k22*[y,v](a)
C
C     where alfa is any real number in [0,pi), and the
C     "boundary condition" functions u, v and/or U, V (if any)
C     are "maximal domain functions".
C     The real numbers k11, k12, k21, k22 are required to
C     satisfy the condition
C
C            k11*k22 - k12*k21 = 1
C
C     in order that the resulting eigenvalue problem be
C     Self-Adjoint.
C     Here the expression
C
C        [y,u](x)
C     represents the Wronskian
C           [y,u](x) = y(x)*pu'(x) - u(x)*py'(x).
C
C     THE INPUT ARGUMENTS A, B, INTAB, P0ATA, QFATA, P0ATB,
C     QFATB, NUMEIG, EIG, TOL, NCA, NCB
C     HAVE THE SAME MEANINGS HERE AS IN SUBROUTINE SLEIGN.
C     THE INPUT ARGUMENTS A1, A2, B1, B2 ARE THE BOUNDARY
C     CONDITION CONSTANTS, AS ARE ALFA, K11, K12, K21, K22.
C     If ICPFUN on entry to this routine is .le. 0, then no
C     eigenvalues are wanted.
C     If ICPFUN is positive, then on input the ICPFUN values in
C     CPFUN specify the x-coordinates, in ascending order, where
C     the eigenfunction values are desired, and on output
C     those x-values will have been replaced by the eigenfunction
C     values.  The x-values must be interior to the interval (A,B).
C
C     .. Scalar Arguments ..
      REAL A,A1,A2,ALFA,B,B1,B2,EIG,K11,K12,K21,K22,P0ATA,
     +                 P0ATB,QFATA,QFATB,TOL
      INTEGER ICPFUN,IFLAG,INTAB,NCA,NCB,NUMEIG
C     ..
C     .. Array Arguments ..
      REAL CPFUN(100)
C     ..
C     .. Local Scalars ..
      REAL AA,BB,TMID
      INTEGER I
C     ..
C     .. Local Arrays ..
      REAL PLOTF(1000,6),SLFUN(9),XT(1000,2)
C     ..
C     .. External Functions ..
      REAL TFROMX
      EXTERNAL TFROMX
C     ..
C     .. External Subroutines ..
      EXTERNAL PERFUN,PERIO
C     ..
C  ----------------------------------------------------------------C
C     OBTAIN THE WANTED EIGENVALUE:
C  ----------------------------------------------------------------C
      CALL PERIO(A,B,INTAB,P0ATA,QFATA,P0ATB,QFATB,A1,A2,B1,B2,NUMEIG,
     +           EIG,TOL,IFLAG,SLFUN,NCA,NCB,ALFA,K11,K12,K21,K22)

C  ----------------------------------------------------------------C
      IF (ICPFUN.GT.0) THEN

C        IN THIS CASE, THE EIGENFUNCTION IS ALSO WANTED.
C        IT MUST BE COMPUTED SEPARATELY, USING PERFUN.
C
C        CONVERT VALUES OF X IN (A,B) INTO VALUES FOR T IN (-1,1):

          DO 10 I = 1,ICPFUN
              XT(9+I,2) = TFROMX(CPFUN(I))
   10     CONTINUE
C
C        AFTER THE CALL TO PERIO, THE VALUES IN SLFUN ARE T-VALUES.
          TMID = SLFUN(1)
          AA = SLFUN(2)
          BB = SLFUN(5)
          CALL PERFUN(AA,BB,TMID,NCA,NCB,K11,K12,K21,K22,EIG,ICPFUN,XT,
     +                PLOTF)

C      REPLACE THE INPUT VALUES OF X BY THE OUTPUT VALUES OF Y(X),
C      WHICH HAVE BEEN PLACED IN PLOTF BY SUBROUTINE PERFUN.
          DO 20 I = 1,ICPFUN
              CPFUN(I) = PLOTF(9+I,1)
   20     CONTINUE

      END IF
C  ----------------------------------------------------------------C
      RETURN
      END
C
      SUBROUTINE PERIO(A,B,INTAB,P0ATA,QFATA,P0ATB,QFATB,A1,A2,B1,B2,
     +                 NUMEIG,EIG,TOL,IFLAG,SLFUN,NCA,NCB,ALFA,K11,K12,
     +                 K21,K22)
C     **********
C  THIS PROGRAM IS USED TO FIND THE EIGENVALUES OF PROBLEMS WITH
C  COUPLED BOUNDARY CONDITIONS, BY FINDING THE APPROPRIATE ZEROS OF
C  THE FUNCTION FF.
C  IT IS CALLED BY THE PROGRAM DRIVE, OR BY SUBROUTINE SLCOUP,
C  OR CAN BE CALLED BY ANY OTHER "DRIVER".
C
C  INPUT QUANTITIES: A,B,INTAB,P0ATA,QFATA,P0ATB,QFATB,A1,A2,B1,B2,
C                    NUMEIG,EIG,TOL,NCA,NCB,ALFA,K11,K12,K21,K22
C
C  OUTPUT QUANTITIES: EIG,TOL,IFLAG
C
C  FUNCTION FF COMPUTES THE VALUES OF THE FUNCTION
C
C             D(LAMBDA) - 2.0*COS(ALFA).
C
C  NECESSARY DATA IS SENT TO FF VIA ENCEE, ABEE,
C    AND PASS2.  THE PROGRAM RECEIVES DATA FROM FF VIA TERM.
C
C     **********
C     .. Scalars in Common ..
      REAL AA,BB,BB1,BB2,DA,DB,DMU,DNU,EPSMIN,ETAA,ETAB,
     +                 FFMAX,FFMIN,FFV,GAMMA,H11,H12,H21,H22,HPI,PI,
     +                 TMID,TRM11,TRM12,TRM21,TRM22,TWOPI
      INTEGER MCA,MCB,T21
      LOGICAL AOK,BOK,PR
C     ..
C     .. Local Scalars ..
      REAL A1D,A1N,A2D,A2N,AE,B1D,B1N,B2D,B2N,BRA,BRB,EIGLO,
     +                 EIGLOD,EIGMU,EIGNU,EIGUP,EIGUPD,HU,HV,LAMBDA,
     +                 LAMUP,ONE,PUP,PVP,RE,STEP,TMP,TOLL,TOLS,TT,U,V,
     +                 VFFMU,VFFNU,XX
      INTEGER I,J,JFLAG,KFLAG,LCASE,NTRY,NUM
      LOGICAL LL1,LL2,LL2S,LL3
C     ..
C     .. External Subroutines ..
      EXTERNAL CERRZ,DXDT,FZERO,SLEIGN,UV
C     ..
C     .. External Functions ..
      REAL FF
      EXTERNAL FF
C     ..
C     .. Common blocks ..
      COMMON /ABEE/AA,BB,TMID
      COMMON /ENCEE/AOK,BOK,MCA,MCB
      COMMON /PASS2/GAMMA,H11,H12,H21,H22,ETAA,ETAB,DA,DB
      COMMON /PIE/PI,TWOPI,HPI
      COMMON /PRIN/PR,T21
      COMMON /RNDOFF/EPSMIN
      COMMON /TERM/TRM11,TRM12,TRM21,TRM22,FFMIN,FFMAX,DMU,DNU,BB1,BB2,
     +       FFV
C     ..
C  THE FUNCTION FF(LAMBDA) HAS BEEN DEFINED SO THAT ITS
C    ZEROS ARE THE EIGENVALUES OF THE COUPLED BOUNDARY
C    CONDITION PROBLEM. BUT IN ORDER TO COMPUTE A PARTICULAR
C    EIGENVALUE, IT IS NECESSARY TO ISOLATE IT FROM ALL THE
C    OTHERS.  FOR THIS PURPOSE WE CALCULATE THE EIGENVALUES
C    OF TWO RELATED PROBLEMS WITH SEPARATED BOUNDARY
C    CONDITIONS, WHICH EIGENVALUES DEFINE AN INTERVAL
C    CONTAINING ONLY THE ONE ZERO OF FF(LAMBDA) WHICH IS
C    WANTED.
C    THE PARTICULAR RELATED PROBLEMS WHICH ARE MOST SUITABLE FOR
C    THIS PURPOSE CAN BE SELECTED BY EXAMINING THE FOUR
C    NUMBERS K11,K12,K21,K22 IN THE GIVEN COUPLED BOUNDARY
C    CONDITIONS.
C
C     .. Scalar Arguments ..
      REAL A,A1,A2,ALFA,B,B1,B2,EIG,K11,K12,K21,K22,P0ATA,
     +                 P0ATB,QFATA,QFATB,TOL
      INTEGER IFLAG,INTAB,NCA,NCB,NUMEIG
C     ..
C     .. Array Arguments ..
      REAL SLFUN(9)
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS,MAX,SQRT
C     ..
C     .. Local Arrays ..
      REAL FTT(5)
C     ..
      ONE = 1.0D0
      IFLAG = 1
C
      MCA = NCA
      MCB = NCB

      AOK = (INTAB.EQ.1 .OR. INTAB.EQ.2) .AND. P0ATA .LT. 0.0D0 .AND.
     +      QFATA .GT. 0.0D0
      BOK = (INTAB.EQ.1 .OR. INTAB.EQ.3) .AND. P0ATB .LT. 0.0D0 .AND.
     +      QFATB .GT. 0.0D0
C
      A1 = 1.0D0
      A2 = 0.0D0
      B1 = 1.0D0
      B2 = 0.0D0
C     THE FOLLOWING CALL TO SLEIGN SETS THE STAGE FOR INTEG.
C     THIS IS NECESSARY BECAUSE, OTHERWISE, THE USUAL
C     SAMPLING IN SLEIGN WOULD NOT TAKE PLACE.
C       (INTEG IS USED IN FUNCTION FF.)
C     BY SETTING THE CALLING ARGUMENT ISLFUN EQUAL TO -1,
C     SLEIGN RETURNS RIGHT AFTER OBTAINING THE QUANTITIES
C       AA, TMID, BB.
C
      LAMBDA = 10.D0
      TOLS = .001D0
      CALL SLEIGN(A,B,INTAB,P0ATA,QFATA,P0ATB,QFATB,A1,A2,B1,B2,NUMEIG,
     +            LAMBDA,TOLS,KFLAG,-1,SLFUN,NCA,NCB)
C     (RETURNING FROM SLEIGN WITH ISLFUN = -1, SLFUN HAS THE
C        T-VARIABLES TMID, AA, BB.)
      TMID = SLFUN(1)
      AA = SLFUN(2)
      BB = SLFUN(5)
C
C  THE FUNCTION FF(LAMBDA) NEEDS TO HAVE THE QUANTITIES
C    BRA,DA,ETAA AND BRB,DB,ETAB FOR THIS PROBLEM, WHICH
C    ARE PASSED TO IT VIA COMMON/PASS2/ .
C  THESE WILL ALL BE = 1.0 IF THE BOUNDARY CONDITION FUNCTIONS
C  u,v AND U,V HAVE BEEN "NORMALIZED" SO THAT
C    [u,v](a) = 1 AND/OR [U,V](b) = 1.
C
      BRA = 1.0D0
      BRB = 1.0D0
      ETAA = 1.0D0
      ETAB = 1.0D0
      DA = 1.0D0
      DB = 1.0D0
      IF (NCA.EQ.3 .OR. NCA.EQ.4) THEN
          TT = AA
          CALL DXDT(TT,TMP,XX)
          CALL UV(XX,U,PUP,V,PVP,HU,HV)
          BRA = U*PVP - V*PUP
          DA = SQRT(ABS(BRA))
          ETAA = BRA/ (DA*DA)
      END IF
      IF (NCB.EQ.3 .OR. NCB.EQ.4) THEN
          TT = BB
          CALL DXDT(TT,TMP,XX)
          CALL UV(XX,U,PUP,V,PVP,HU,HV)
          BRB = U*PVP - V*PUP
          DB = SQRT(ABS(BRB))
          ETAB = BRB/ (DB*DB)
      END IF
      IF (PR) WRITE (T21,FMT=*) ' BRA,BRB = ',BRA,BRB
C
C  WE DISTINGUISH 4 CASES:
C  CASE(1): K11 .GT. 0 & K12 .LE. 0
C      (2): K11 .LT. 0 & K12 .LT. 0
C      (3): K11 .LT. 0 & K12 .GE. 0
C      (4): K11 .GT. 0 & K12 .GT. 0
C  IN CASE(1) WE HAVE
C    NU(0).LE.EIG(0).LE.NU(1),MU(1).LE.EIG(1).LE.NU(2),MU(2)...
C  IN CASE(2) WE HAVE
C    EIG(0).LE.NU(0),MU(0).LE.EIG(1).LE.NU(1),MU(1).LE.EIG(2)...
C  IN CASE(3) WE SET
C    HIJ = -KIJ & GAMMA = PI-ALFA
C    SO THAT THE HIJ ARE LIKE THE KIJ IN CASE(1).
C  IN CASE(4) WE SET
C    HIJ = -KIJ & GAMMA = PI-ALFA
C    SO THAT THE HIJ ARE LIKE THE KIJ IN CASE(2).
C  NOTICE THAT IT IS GAMMA AND THE HIJ THAT ARE
C    PASSED TO FF(LAMBDA) VIA PASS2, RATHER THAN
C    THE GIVEN ALFA AND THE KIJ.
      GAMMA = ALFA
      H11 = K11
      H12 = K12
      H21 = K21
      H22 = K22

      IF (K11.GT.0.0D0 .AND. K12.LE.0.0D0) THEN
          LCASE = 1
      ELSE IF (K11.LT.0.0D0 .AND. K12.LT.0.0D0) THEN
          LCASE = 2
      ELSE IF (K11.LT.0.0D0 .AND. K12.GE.0.0D0) THEN
          LCASE = 3
      ELSE IF (K11.GT.0.0D0 .AND. K12.GT.0.0D0) THEN
          LCASE = 4
      END IF

      IF (LCASE.EQ.3 .OR. LCASE.EQ.4) THEN
          H11 = -K11
          H12 = -K12
          H21 = -K21
          H22 = -K22
          GAMMA = PI
          IF (ALFA.NE.0.0D0) GAMMA = PI-ALFA
      END IF
      NTRY = 0
      A1N = 0.0D0
      A2N = 1.0D0
      A1D = 1.0D0
      A2D = 0.0D0
      B1N = -H21
      B2N = H11
      B1D = H22
      B2D = -H12
C
C   SETTING TOLL = 0.0 CAUSES SLEIGN TO GET THE MOST ACCURACY
C   IT IS CAPABLE OF.
      TOLL = 0.0D0
C
C  PROBLEM FOR NU:
      NUM = NUMEIG
      IF ((LCASE.EQ.2.OR.LCASE.EQ.4) .AND. NUMEIG.GE.1) NUM = NUMEIG - 1
      TOLS = TOLL
      EIGNU = 0.0D0
      CALL SLEIGN(A,B,INTAB,P0ATA,QFATA,P0ATB,QFATB,A1N,A2N,B1N,B2N,NUM,
     +            EIGNU,TOLS,KFLAG,0,SLFUN,NCA,NCB)
      IF (PR) WRITE (T21,FMT=*) ' EIGNU,KFLAG = ',EIGNU,KFLAG
      EIGLO = EIGNU
   35 CONTINUE
      IF ((LCASE.EQ.2.OR.LCASE.EQ.4) .AND. NUMEIG.EQ.0) THEN
          NTRY = NTRY + 1
          EIGLO = EIGLO - 10.0D0
      END IF

C
C  PROBLEM FOR MU:
      NUM = NUMEIG
      IF (NCA.EQ.4 .OR. NCB.EQ.4) NUM = NUMEIG + 1
      TOLS = TOLL
      EIGMU = 0.0D0
      CALL SLEIGN(A,B,INTAB,P0ATA,QFATA,P0ATB,QFATB,A1D,A2D,B1D,B2D,NUM,
     +            EIGMU,TOLS,KFLAG,0,SLFUN,NCA,NCB)
      IF (PR) WRITE (T21,FMT=*) ' EIGMU,KFLAG = ',EIGMU,KFLAG
C
      VFFNU = FF(EIGNU)
      IF (PR) WRITE (T21,FMT=*) ' EIGNU,BB1,BB2 = ',EIGNU,BB1,BB2
      VFFMU = FF(EIGMU)
      IF (PR) WRITE (T21,FMT=*) ' EIGMU,BB1,BB2 = ',EIGMU,BB1,BB2
      IF (PR) WRITE (T21,FMT=*) ' EIGNU,VFFNU,EIGMU,VFFMU = ',EIGNU,
     +    VFFNU,EIGMU,VFFMU
C
      EIGUP = EIGMU
C
      IF (PR) WRITE (T21,FMT=*) ' EIGLO,EIGUP = ',EIGLO,EIGUP
C
      EIG = EIGLO
      LAMUP = EIGUP
      RE = EPSMIN
      AE = RE
C
C  THE WANTED EIGENVALUE IS PRESUMABLY BRACKETED BY (EIGLO,EIGUP),
C  BUT BECAUSE OF COMPUTING IMPERFECTIONS WE WILL USE A SLIGHTLY
C  LARGER BRACKET.
      EIGLOD = 0.01D0*MAX(ONE,ABS(EIGLO))
      EIGUPD = 0.01D0*MAX(ONE,ABS(EIGUP))

      EIGLO = EIGLO - EIGLOD
      EIGUP = EIGUP + EIGUPD
      IF (PR) WRITE (T21,FMT=*) ' 2,EIGLO,EIGUP = ',EIGLO,EIGUP

      EIG = EIGLO
      LAMUP = EIGUP
      LAMBDA = 0.5D0* (EIG+LAMUP)

      CALL FZERO(FF,EIG,LAMUP,LAMBDA,RE,AE,JFLAG)

C  ESTIMATE ERROR IN EIG, USING SUBROUTINE CERRZ:
      LL2S = .FALSE.
      TOL = 1.0D-1
      STEP = TOL
      DO 78 J = 1,5
          STEP = STEP* (1.0D-1)
          DO 73 I = 1,5
              TMP = EIG + (I-3)*STEP*MAX(ONE,ABS(EIG))
              FTT(I) = FF(TMP)
              WRITE (T21,FMT=*) TMP,FTT(I)
   73     CONTINUE

          CALL CERRZ(FTT,LL1,LL2,LL3)
C  LL1 = .TRUE. MEANS EIG IS SIMPLE.
C  LL2 = .TRUE. MEANS EIG IS PROBABLY DOUBLE.
          IF (LL1 .AND. LL3) THEN
              TOL = STEP
          ELSE IF (LL2 .AND. LL3) THEN
              LL2S = .TRUE.
              TOL = STEP
          ELSE
              GO TO 80
          END IF
   78 CONTINUE
   80 CONTINUE
      TOL = TOL/10.0D0
      WRITE (*,FMT=*) ' tol = ',TOL

      KFLAG = 0
      IF (EIG.GT.EIGLO .AND. EIG.LT.EIGUP) KFLAG = 1
      IF (KFLAG.EQ.1 .AND. LL2S) KFLAG = 2
      IF (PR) WRITE (T21,FMT=*) ' EIG,JFLAG,KFLAG = ',EIG,JFLAG,KFLAG
      TMP = BB1*BB2
      IF (PR) WRITE (T21,FMT=*) ' BB1*BB2 = ',TMP
C
C  KFLAG=0 MEANS WE HAVE NOT FOUND A ROOT.
C
      IF (KFLAG.EQ.0 .AND. NTRY.EQ.1) THEN
          IF (PR) WRITE (T21,FMT=*) ' FAILED ONCE WITH EIGLO =',EIGLO
          GO TO 35
      ELSE IF (KFLAG.EQ.0 .AND. NTRY.GT.1) THEN
          IFLAG = 4
          IF (PR) WRITE (T21,FMT=*) ' REQUESTED EIGENVALUE NOT FOUND. '
          IF (PR) WRITE (*,FMT=*) ' REQUESTED EIGENVALUE NOT FOUND. '
          GO TO 100
      ELSE
          IFLAG = KFLAG
      END IF
C
C  KFLAG=1 MEANS WE HAVE FOUND A ROOT IN (EIGLO,EIGUP).
C  KFLAG=2 MEANS THE ROOT MAY BE A DOUBLE.
C
C  N.B. JFLAG = 1 MEANS EIG IS WITHIN THE REQUESTED
C       ERROR TOLERANCE AND ALL IS WELL;
C       JFLAG = 2 MEANS FF(EIG) = 0 BUT INTERVAL
C       (EIG,LAMUP) MAY NOT HAVE COLLAPSED TO REQUESTED
C       SIZE AS SPECIFIED BY RE,AE.
C  SO :
      IF (KFLAG.EQ.0) TOL = 0.0D0
      IF (PR) WRITE (T21,FMT=*) ' EIGLO,EIG,EIGUP = ',EIGLO,EIG,EIGUP
C
C  NOTE THAT IF ALFA = 0 (IN THE BOUNDARY CONDITION),
C    THEN IT IS POSSIBLE THAT AN EIGENVALUE BE A DOUBLE.
C    ON THE OTHER HAND, IF ALFA .NE. 0, THEN ANY EIGEN-
C    VALUE MUST BE SIMPLE;  MOREOVER, IN THIS CASE THE
C    EIGENFUNCTION IS COMPLEX !, AND WE HAVE NO PROVISION
C    (AT THE MOMENT) FOR PLOTTING IT.  SO, IN THIS CASE,
C    WE DO NOT NEED TO PREPARE FOR PLOTTING.
C
C  NOTE THAT FACTORS DA,DB,ETAA,ETAB HAVE BEEN USED IN
C  FUNCTION FF
C    WHEN THE GIVEN MAX.DOMAIN FUNCTIONS U,V GIVEN
C         BY THE USER HAVE [U,V] = 1 -- THEN
C
C     YL = AEE*Y1L + BEE*Y2L
C    AND
C     YR = E*Y1R + F*Y2R,
C    (WHERE Y1L,Y1R,Y2L,Y2R ARE THE FUNCTIONS COMPUTED IN FF)
C    WITH
C    AEE = Y2L - K21*Y1R - K11*Y2R
C    BEE = -(Y1L - K22*Y1R - K12*Y2R).

  100 CONTINUE
C
      RETURN
      END
C
      SUBROUTINE CERRZ(FTT,LOGIC1,LOGIC2,LOGIC3)
C
C  THIS PROGRAM ATTEMPTS TO DETERMINE WHETHER OR NOT THE
C  NUMBERS IN THE ARRAY FTT ARE TOO CONTAMINATED WITH ERRORS
C  TO BE USEFUL.  IT IS CALLED BY PERIO.
C
C  INPUT QUANTITIES: FTT
C  OUTPUT QUANTITIES: LOGIC1,LOGIC2,LOGIC3
C
C  IT IS ASSUMED THAT THE NUMBERS FTT(I), I=1,5, CORRESPOND TO
C  EQUALLY SPACED ARGUMENTS IN THE INDEPENDENT VARIABLE, AND THAT
C  FTT(3) IS NEAR A SIMPLE ZERO OF THE FUNCTION, OR NEAR A
C  DOUBLE ZERO, WHERE THE FUNCTION HAS A MAXIMUM.
C  THE THREE LOGICAL VARIABLES INDICATE
C
C  LOGIC1 = .TRUE. MEANS  THERE IS A SIMPLE ZERO NEAR FTT(3),
C  LOGIC2 = .TRUE. MEANS  THERE IS A MAXIMUM NEAR FTT(3),
C  LOGIC3 = .TRUE. MEANS  THE NUMBERS IN ARRAY FTT ARE
C                         PROBABLY USEFUL.
C
C     .. Scalars in Common ..
      INTEGER T21
      LOGICAL PR
C     ..
C     .. Scalar Arguments ..
      LOGICAL LOGIC1,LOGIC2,LOGIC3
C     ..
C     .. Array Arguments ..
      REAL FTT(5)
C     ..
C     .. Local Scalars ..
      REAL TMP1,TMP2
      INTEGER I
C     ..
C     .. Local Arrays ..
      REAL D2F(3),DF(4)
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS,MAX,MIN
C     ..
C     Common blocks
C     .. Common blocks ..
      COMMON /PRIN/PR,T21
C     ..
      DO 10 I = 1,4
          DF(I) = FTT(I+1) - FTT(I)
   10 CONTINUE
      DO 20 I = 1,3
          D2F(I) = DF(I+1) - DF(I)
   20 CONTINUE
      LOGIC3 = (FTT(1)*FTT(2).GT.0.0D0) .AND. (FTT(4)*FTT(5).GT.0.0D0)
      LOGIC1 = LOGIC3 .AND. (FTT(1)*FTT(4).LT.0.0D0)
      LOGIC2 = LOGIC3 .AND. (FTT(1).LT.0.0D0) .AND. (FTT(4).LT.0.0D0)

C  LOGIC1 MEANS FUNCTION FTT HAS A SIMPLE ZERO NEAR FTT(3)
C  LOGIC2 MEANS FUNCTION FTT HAS A MAXIMUM NEAR FTT(3)

      LOGIC3 = (D2F(1).LT.0.0D0) .AND. (D2F(2).LT.0.0D0) .AND.
     +         (D2F(3).LT.0.0D0)
      TMP1 = MIN(ABS(D2F(1)),ABS(D2F(2)),ABS(D2F(3)))
      TMP2 = MAX(ABS(D2F(1)),ABS(D2F(2)),ABS(D2F(3)))
      LOGIC3 = LOGIC3 .AND. (TMP2/TMP1.LE.1.2D0)
      IF (PR) WRITE (T21,FMT=*) 'LOGIC1, LOGIC2, LOGIC3 = ',LOGIC1,
     +    LOGIC2,LOGIC2
      IF (PR .AND. LOGIC1) WRITE (T21,FMT=*) ' FTT HAS A SIMPLE ZERO '
      IF (PR .AND. LOGIC2) WRITE (T21,FMT=*) ' FTT HAS A MAXIMUM '

C  LOGIC3 MEANS ALL THE SECOND DIFFERENCES ARE NEGATIVE
C    AND ARE NEARLY CONSTANT.  THIS IS INTERPRETED TO MEAN
C    THAT THE FUNCTION VALUES FTT(I), I=1,5, ARE PROBABLY
C    USEFUL NUMBERS.  I.E. NOT OVERLY CONTAMINATED BY ERRORS.

      RETURN
      END
C
      SUBROUTINE MESH(EIG,TS,NN)
C  THIS PROGRAM IS CALLED FROM DRAW.
C
C  INPUT QUANTITIES: EIG
C  OUTPUT QUANTITIES: TS,NN
C
C  IT GENERATES A MESH OF POINTS IN (-1,1), TO BE USED
C    FOR PLOTTING SOLUTIONS OF THE D.EQUATION
C           -(py')' + q*y = lambda*w*y
C    ON THE TRANSFORMED INTERVAL (-1,1).  THE POINTS ARE
C    GENERALLY NOT EVENLY SPACED, BUT ARE INTENDED TO BE
C    SPACED SO THAT THEY ARE CLOSER TOGETHER WHERE THE
C    D.E. IS "TRIGONOMETRIC", AND FARTHER APART WHERE IT
C    IS "EXPONENTIAL".  THE SPACING IS BASED ON THE SIGN
C    AND SIZE OF THE QUANTITY
C            (EIG*WX - QX)/PX
C
C
C     .. Scalar Arguments ..
      REAL EIG
      INTEGER NN
C     ..
C     .. Array Arguments ..
      REAL TS(1000)
C     ..
C     .. Local Scalars ..
      REAL DELT,H,T,TM,TP
      INTEGER I,NMAX,NN1,NN2
C     ..
C     .. Local Arrays ..
      REAL TS1(1000),TS2(1000)
C     ..
C     .. External Subroutines ..
      EXTERNAL SIZEH
C     ..
      NMAX = 300
      DELT = 0.01D0

      TP = 1.0D0 - (1.0D-7)
      TM = -TP
      T = 0.0D0
      H = DELT
      DO 10 I = 1,NMAX
          IF (T.GE.TP) GO TO 15
          CALL SIZEH(EIG,T,H)
          TS1(I) = T
          IF (T.GE.TP) TS1(I) = TP
          NN1 = I
   10 CONTINUE
   15 CONTINUE
      T = 0.0D0
      H = -DELT
      DO 20 I = 1,NMAX
          IF (T.LE.TM) GO TO 25
          CALL SIZEH(EIG,T,H)
          TS2(I) = T
          IF (T.LE.TM) TS2(I) = TM
          NN2 = I
   20 CONTINUE
   25 CONTINUE
      TS1(1) = 1.0D-7
      TS2(1) = -1.0D-7
      NN = NN1 + NN1
      DO 30 I = 1,NN2
          TS(NN2-I+1) = TS2(I)
   30 CONTINUE
      DO 40 I = 1,NN1
          TS(NN2+I) = TS1(I)
   40 CONTINUE
      NN = NN1 + NN2
      RETURN
      END
C
      SUBROUTINE SIZEH(EIG,T,H)
C  INPUT QUANTITIES: EIG,T,H,HPI
C  OUTPUT QUANTITIES: T,H
C
C  N.B.  FK = 3.0 GIVES ABOUT THE "MINIMAL" NUMBER OF POINTS
C        FOR A DECENT GRAPH.  DOUBLING FK MORE OR LESS DOUBLES
C        THE NUMBER OF POINTS.
C     FK = 3.0
C     .. Scalar Arguments ..
      REAL EIG,H,T
C     ..
C     .. Scalars in Common ..
      REAL HPI,PI,TWOPI
C     ..
C     .. Local Scalars ..
      REAL DT,FK,H1,H2,PX,Q0,Q0SQ,QQ,QX,SQ,T1,TMP,WX,X,X1
C     ..
C     .. External Functions ..
      REAL P,Q,W
      EXTERNAL P,Q,W
C     ..
C     .. External Subroutines ..
      EXTERNAL DXDT
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS,MAX,MIN,SQRT
C     ..
C     .. Common blocks ..
      COMMON /PIE/PI,TWOPI,HPI
C     ..
      FK = 5.0D0
      Q0 = HPI
      Q0SQ = Q0*Q0
      CALL DXDT(T,DT,X)
      PX = P(X)
      QX = Q(X)
      WX = W(X)
      QQ = DT*DT* (EIG*WX-QX)/PX
      IF (QQ.GE.Q0SQ) THEN
          SQ = FK*SQRT(QQ)
      ELSE IF (QQ.LT.Q0SQ .AND. QQ.GE.-Q0SQ) THEN
          SQ = FK*SQRT(0.9D0*Q0SQ+0.1D0*QQ)
      ELSE
          SQ = FK*SQRT(0.8D0)*Q0
      END IF
      H1 = 1.0D0/SQ
C  DON'T LET H CHANGE TOO MUCH IN SIZE FROM THE LAST TIME.
      TMP = 2.0D0*ABS(H)
      H1 = MIN(H1,TMP)
      TMP = 0.5D0*ABS(H)
      H1 = MAX(H1,TMP)
      IF (H.LT.0.0D0) H1 = -H1
      T1 = T + H1
      CALL DXDT(T1,DT,X1)
      PX = P(X1)
      QX = Q(X1)
      WX = W(X1)
      QQ = DT*DT* (EIG*WX-QX)/PX
      IF (QQ.GE.Q0SQ) THEN
          SQ = FK*SQRT(QQ)
      ELSE IF (QQ.LT.Q0SQ .AND. QQ.GE.-Q0SQ) THEN
          SQ = FK*SQRT(0.9D0*Q0SQ+0.1D0*QQ)
      ELSE
          SQ = FK*SQRT(0.8D0)*Q0
      END IF
      H2 = 1.0D0/SQ
C  DON'T LET H CHANGE TOO MUCH IN SIZE FROM THE LAST TIME.
      TMP = 2.0D0*ABS(H)
      H2 = MIN(H2,TMP)
      TMP = 0.5D0*ABS(H)
      H2 = MAX(H2,TMP)
      IF (H.LT.0.0D0) H2 = -H2
      H = 0.5D0* (H1+H2)
      T = T + H
      RETURN
      END
C
      SUBROUTINE QPLOT(ISLFN,XT,NV,PLOTF,NF)
C     **********
C  THIS PROGRAM IS CALLED FROM "DRIVE".
C
C  INPUT QUANTITIES: ISLFN,XT,NV,PLOTF,NF
C  OUTPUT QUANTITIES: A "GRAPH"
C
C    IT IS INTENDED TO PROVIDE A ROUGH PLOT
C    OF EIGENFUNCTIONS.  IN ORDER TO BE INDEPENDENT OF
C    GRAPHICS PACKAGES, WHICH ARE EXTREMELY PLATFORM DEPENDENT,
C    IT USES ONLY THE SET OF CHARACTERS AVAILABLE ON THE USUAL
C    KEYBOARD.  RATHER THAN USING JUST ONE SYMBOL, LIKE A "."
C    OR A "x", IT USES THE SET OF SYMBOLS ".", "+", "_" (underscore),
C    AND """ (double quote).  THIS MAKES THE RESULTING "graph"
C    A LITTLE BIT SMOOTHER.
C
C    NOTE THAT THE NUMBER OF POINTS RESULTING HERE ARE NOT
C    NECESSARILY THE SAME AS THE NUMBER OF POINTS, ISLFUN,
C    IN THE INPUT, XT.  HERE, THERE WILL BE 75, ONE FOR EACH
C    COLUMN OF A TYPEWRITTEN PAGE (NOT COUNTING THE MARGIN).
C    IF ISLFUN IS LESS THAN 75, OR IF THE GIVEN X-COORDS. ARE
C    NOT UNIFORMLY SPACED, THE SCHEME HERE WILL "INTERPOLATE"
C    TO GET ONE POINT FOR EACH COLUMN.
C     **********
C     .. Parameters ..
      INTEGER NMAX,MMAX
      PARAMETER (NMAX=75,MMAX=22)
C     ..
C     .. Local Scalars ..
      REAL DZ,ONE,REM,X,XK,XKP,XMAX,XMIN,Y,YK,YKP,YMAX,YMIN
      INTEGER I,II,IZ,J,K,L
C     ..
C     .. Local Arrays ..
      REAL A(1000,2)
      CHARACTER AX(NMAX,MMAX)
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS,INT,MAX,MIN
C     ..
C     .. Scalar Arguments ..
      INTEGER ISLFN,NF,NV
C     ..
C     .. Array Arguments ..
      REAL PLOTF(1000,6),XT(1000,2)
C     ..
      ONE = 1.0D0
C  FIRST, DETERMINE SCALING FACTORS FOR THE X-COORDS.
C    AND FOR THE Y-CORRDS.
      XMAX = -1000000.D0
      XMIN = 1000000.D0
      YMAX = -1000000.D0
      YMIN = 1000000.D0
      DZ = YMIN
      DO 10 I = 1,ISLFN
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
   10 CONTINUE
C  IZ IS THE INDEX FOR WHICH DZ = ABS(Y) IS LEAST.
C
      IF (YMIN*YMAX.LE.0.0D0) THEN
          Y = MAX(ONE,YMAX-YMIN)
      ELSE
          Y = MAX(ONE,ABS(YMIN),ABS(YMAX))
      END IF
C  RESCALE THE Y-COORDS.
      DO 20 I = 1,ISLFN
          A(I,2) = A(I,2)/Y
   20 CONTINUE
      YMAX = YMAX/Y
      YMIN = YMIN/Y
C
C  SHIFT BOTH COORDS. SO THAT THEIR MINIMA ARE EQUAL TO 0.0
      DO 30 I = 1,ISLFN
          A(I,1) = A(I,1) - XMIN
          A(I,2) = A(I,2) - YMIN
   30 CONTINUE
C
C     NOW MIN(X) = 0. AND MIN(Y) = 0.
C
      X = XMAX - XMIN
      DO 40 I = 1,ISLFN
          A(I,1) = NMAX*A(I,1)/X
          A(I,2) = MMAX*A(I,2)
   40 CONTINUE
C
C  INITIALIZE THE ARRAY ENTRIES TO BE BLANKS.
      DO 60 J = 1,NMAX
          DO 50 K = 1,MMAX
              AX(J,K) = ' '
   50     CONTINUE
   60 CONTINUE
C
C  FOR EACH J, FIND INDEX II SUCH THAT A(II,1).LE.J-1/2
C    AND A(II+1,1).GT.J-1/2 .
      DO 80 J = 2,NMAX
          II = 0
          X = J - 0.5D0
          DO 70 I = 1,ISLFN
              IF (A(I,1).LE.X) II = I
   70     CONTINUE
C
C   POINT PK IS (XK,YK); POINT PKP IS (XKP,YKP).
C     LINE PK,PKP IS: Y-YK = (X-XK)*(YKP-YK)/(XKP-XK)
C     THIS LINE MEETS THE LINE X = J - 0.5 WHERE:
C
          XK = A(II,1)
          XKP = A(II+1,1)
          YK = A(II,2)
          YKP = A(II+1,2)
          Y = YK + (X-XK)* (YKP-YK)/ (XKP-XK)
C
C  THE CHARACTERS _ . + " , IN ORDER, ARE ASCENDING.
C  USE WHICHEVER MARK IS CLOSEST TO Y FOR THE ENTRY AX(J, )
          K = MMAX - INT(Y)
          REM = Y + (K-MMAX)
          IF (REM.LE.0.25D0) THEN
              AX(J,K) = '_'
          ELSE IF (REM.LE.0.50D0) THEN
              AX(J,K) = '.'
          ELSE IF (REM.LE.0.75D0) THEN
              AX(J,K) = '+'
          ELSE
              AX(J,K) = '"'
          END IF
   80 CONTINUE
C
C  DRAW A HORIZONTAL LINE WHICH IS EITHER THE ZERO LINE,
C    OR IS BELOW THE LOWEST POINT ON THE CURVE,
C    OR IS ABOVE THE HIGHEST POINT ON THE CURVE.

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
   90 CONTINUE
      WRITE (*,FMT=*)
      DO 100 K = 1,MMAX
          WRITE (*,FMT='(1X,80A1)') (AX(J,K),J=1,NMAX)
  100 CONTINUE
      RETURN
      END
C
      SUBROUTINE SCREEN(AA,BB,NCA,NCB,NIVP,NEND,ISLFN,SLFN)
C  THIS PROGRAM IS CALLED FROM DRAW.
C
C  INPUT QUANTITIES: AA,BB,NCA,NCB,NIVP,NEND,ISLFN,SLFN
C  OUTPUT QUANTITIES: ISLFN,SLFN
C
C    IT TAKES THE POINTS GENERATED BY SUBROUTINE MESH
C    AND REMOVES ANY WHICH ARE LIKELY TO CAUSE TROUBLE IN
C    EITHER COMPUTING FUNCTION VALUES THERE OR ARE TOO CLOSE
C    TOGETHER FOR PLOTTING.
C
C     THE POINTS GENERATED BY MESH MAY NOT BE WITHIN THE INTERVAL
C     (AA,BB).  WE WANT TO ENSURE THAT WE USE ONLY SUCH POINTS
C     AS ARE WITHIN THE INTERVAL.
C
C     IN ADDITION, IF THE ENDPOINT IS NOT REG, WE WILL NOT TRY TO
C      PLOT POINTS WITH T.LT.-.95 OR T.GT..95 (IT COULD CAUSE TROUBLE
C      IN TRYING TO INTEGRATE FROM TOO NEAR A SINGULAR POINT). SO:

C     .. Scalar Arguments ..
      REAL AA,BB
      INTEGER ISLFN,NCA,NCB,NEND,NIVP
C     ..
C     .. Array Arguments ..
      REAL SLFN(1000,2)
C     ..
C     .. Local Scalars ..
      REAL AAM,BBM,ONE,TMP
      INTEGER I,JJ,K
      LOGICAL REGA,REGB
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS,MAX,MIN
C     ..
      ONE = 1.0D0
      REGA = NCA .EQ. 1
      REGB = NCB .EQ. 1
      AAM = AA
      BBM = BB
      TMP = -0.995D0
      IF (.NOT.REGA) AAM = MAX(AA,TMP)
      TMP = 0.995D0
      IF (.NOT.REGB) BBM = MIN(BB,TMP)
      K = 0
      DO 20 I = 1,ISLFN
          IF (SLFN(9+I,1).GT.AAM .AND. SLFN(9+I,1).LT.BBM) THEN
              K = K + 1
              SLFN(9+K,1) = SLFN(9+I,1)
          END IF
   20 CONTINUE
      ISLFN = K
C
C     WE ALSO WANT TO BE SURE AN INITIAL ENDPOINT IS AA OR BB UNLESS THE
C     POINT IS LIMIT CIRCLE, AND THAT THE LAST POINT IS NOT AA OR BB.
C
C     NEAR AN OSCILLATORY ENDPOINT THE POINTS MAY BE SO CLOSE THAT WE
C     COULDN'T SEE THE ACTUAL CURVE EVEN IF WE PLOTTED IT.  SO WE WANT
C     TO REMOVE POINTS FROM THAT END UP TO WHERE THEY ARE NOT SO CLOSE.
C
      IF (NCA.EQ.4) THEN
          JJ = 0
          DO 40 I = 1,ISLFN
c             IF (ABS(SLFN(9+I,1)-SLFN(10+I,1)).LT.0.001) JJ = JJ + 1
              IF (ABS(SLFN(9+I,1)-SLFN(10+I,1)).LT.0.0001D0) JJ = JJ + 1

   40     CONTINUE
          IF (JJ.GT.0) THEN
              ISLFN = ISLFN - JJ
              DO 50 I = 1,ISLFN
                  SLFN(9+I,1) = SLFN(9+I+JJ,1)
   50         CONTINUE
          END IF
      END IF
      IF (NCB.EQ.4) THEN
          JJ = 0
          DO 60 I = ISLFN,1,-1
              IF (ABS(SLFN(9+I,1)-SLFN(8+I,1)).LT.0.001D0) JJ = JJ + 1
   60     CONTINUE
          IF (JJ.GT.0) ISLFN = ISLFN - JJ
      END IF
C
C     FINALLY, IN THE CASE OF AN INITIAL VALUE PROBLEM,
C     WE CANNOT AFFORD TO INTEGRATE TO THE OTHER END B (OR A)
C     UNLESS BOTH ENDS ARE REGULAR.

C  RECALL: NIVP=1 MEANS IVP FROM ONE END.
C          NIVP=2 MEANS IVP FROM BOTH ENDS.
C          NEND=1 MEANS INITIAL POINT A.
C          NEND=2 MEANS INITIAL POINT B.
C
      IF (NIVP.EQ.1 .AND. .NOT. (REGA.AND.REGB)) THEN
          ISLFN = ISLFN - 1
          IF (NEND.EQ.2) THEN
              DO 70 I = 1,ISLFN
                  SLFN(9+I,1) = SLFN(10+I,1)
   70         CONTINUE
          END IF
      END IF
C
      RETURN
      END
C
      SUBROUTINE RENORM(EIGV,NEND,NIVP,A1,A2,B1,B2,NCA,NCB,ISLFN,SLFN,
     +                  CCYA,CCYB,ENDA,ENDB,PERIOD)
C  THIS PROGRAM IS CALLED FROM DRAW.
C
C  INPUT QUANTITIES: EIGV,NEND,NIVP,A1,A2,B1,B2,NCA,NCB,ISLFN,SLFN,
C                    ENDA,ENDB,PERIOD
C  OUTPUT QUANTITIES: CCYA,CCYB
C
C  SLEIGN NORMALIZES THE WAVEFUNCTION TO HAVE L2-NORM 1.0, BUT FOR AN
C  INITIAL VALUE PROBLEM WE WANT TO HAVE THE VALUE (Y) AND SLOPE (P*Y')
C  (OR [Y,U] AND [Y,V]) AT THE END TO BE THOSE SPECIFIED FOR THE
C  INITIAL CONDITIONS.  SO HERE WE MUST RE-NORMALIZE FOR THIS PURPOSE.
C  N.B. A1 = ALFA2 & A2 = -ALFA1
C
C    THE RESULT OF THIS RE-NORMALIZATION IS THE PAIR OF
C    QUANTITIES CCYA AND CCYB.
C
C     THE VALUES IN THE ARRAYS TT AND YY COME FROM THE INTEGRATIONS
C     IN SUBROUTINE WR, EXCEPT WHEN THE ENDPOINT IS OSCILLATORY.

C     .. Scalar Arguments ..
      REAL A1,A2,B1,B2,CCYA,CCYB
      INTEGER ISLFN,NCA,NCB,NEND,NIVP
      LOGICAL EIGV,ENDA,ENDB,PERIOD
C     ..
C     .. Array Arguments ..
      REAL SLFN(1000,2)
C     ..
C     .. Arrays in Common ..
      REAL TT(7,2),YY(7,3,2)
      INTEGER NT(2)
C     ..
C     .. Local Scalars ..
      REAL BRYU,BRYUA,BRYUB,BRYV,BRYVA,BRYVB,DD,FZ,GEE,HU,
     +                 HV,PHI,PUP,PVP,PYPA,PYPB,RHOZ,SIG,T,THETAZ,TMP,U,
     +                 V,X,YA,YB
      INTEGER I,NPTS
      LOGICAL LCA,LCB,REGA,REGB,WREGA,WREGB
C     ..
C     .. External Subroutines ..
      EXTERNAL DXDT,UV
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC COS,EXP,SIN
C     ..
C     .. Common blocks ..
      COMMON /TEMP/TT,YY,NT
C     ..
C
C  RECALL: EIGV MEANS AN EIGENFUNCTION WAS COMPUTED.
C          PERIOD MEANS COUPLED BOUNDARY CONDITIONS.
C          NIVP=1 MEANS IVP FROM ONE END.
C          NIVP=2 MEANS IVP FROM BOTH ENDS.
C          NEND=1 MEANS INITIAL POINT A.
C          NEND=2 MEANS INITIAL POINT B.
C
      REGA = NCA .EQ. 1
      REGB = NCB .EQ. 1
      WREGA = NCA .EQ. 2
      WREGB = NCB .EQ. 2
      LCA = NCA .EQ. 3 .OR. NCA .EQ. 4
      LCB = NCB .EQ. 3 .OR. NCB .EQ. 4
C
C  THE ENTRIES IN THE ARRAY TT ARE
C     TT(I,1) VALUES OF T NEAR A,
C     TT(I,2) VALUES OF T NEAR B.
C  THE ENTRIES IN THE ARRAY SLFN ARE (FOR I = 1,ISLFN):
C         SLFN(9+I,1) IS THETA(I)
C         SLFN(9+I,2) IS F(I), F = LN(RHO)
C         Y = RHO*SIN(THETA)
C       PY' = RHO*COS(THETA)
C
C  AT OR NEAR THE ENDPOINT A:
      RHOZ = 1.0D0
      CCYA = 1.0D0
      IF (.NOT.EIGV .AND. (NEND.EQ.1.OR.NIVP.EQ.2) .AND. REGA) THEN
          THETAZ = SLFN(9+1,1)
          FZ = SLFN(9+1,2)
          RHOZ = EXP(FZ)
          YA = RHOZ*SIN(THETAZ)
          PYPA = RHOZ*COS(THETAZ)
          IF (A1.EQ.0.0D0) THEN
              CCYA = -A2/YA
          ELSE IF (A2.EQ.0.0D0) THEN
              CCYA = A1/PYPA
          ELSE
              CCYA = -A2/YA
          END IF
      END IF
      ENDA = (NEND.EQ.1 .OR. NIVP.EQ.2) .AND. (LCA .OR. WREGA)
      IF ((.NOT.EIGV.OR.PERIOD) .AND. ENDA) THEN
          NPTS = 6
          IF (NCA.EQ.4) NPTS = NT(1)
          SIG = 1.0D0
          DO 110 I = 2,NPTS
              T = TT(I,1)
              PHI = YY(I,1,1)
              GEE = YY(I,3,1)
              SIG = EXP(GEE)
              IF (WREGA .AND. I.EQ.2) THEN
                  YA = SIG*SIN(PHI)
                  PYPA = SIG*COS(PHI)
                  IF (A1.EQ.0.0D0) THEN
                      CCYA = -A2/YA
                  ELSE IF (A2.EQ.0.0D0) THEN
                      CCYA = A1/PYPA
                  ELSE
                      CCYA = -A2/YA
                  END IF
              END IF
              IF (LCA) THEN
                  CALL DXDT(T,TMP,X)
                  CALL UV(X,U,PUP,V,PVP,HU,HV)
                  DD = U*PVP - V*PUP
                  BRYU = -DD*SIN(PHI)*SIG
                  BRYV = DD*COS(PHI)*SIG
                  IF (I.EQ.2) THEN
                      BRYUA = BRYU
                      BRYVA = BRYV
                      IF (A1.EQ.0.0D0) THEN
                          CCYA = -A2/BRYUA
                      ELSE IF (A2.EQ.0.0D0) THEN
                          CCYA = A1/BRYVA
                      ELSE
                          CCYA = -A2/BRYUA
                      END IF
                  END IF
                  BRYU = BRYU*CCYA
                  BRYV = BRYV*CCYA
              END IF
  110     CONTINUE
      END IF

C  AT OR NEAR THE ENDPOINT B:
      RHOZ = 1.0D0
      CCYB = 1.0D0
      IF (.NOT.EIGV .AND. (NEND.EQ.2.OR.NIVP.EQ.2) .AND. REGB) THEN
          THETAZ = SLFN(9+ISLFN,1)
          FZ = SLFN(9+ISLFN,2)
          RHOZ = EXP(FZ)
          YB = RHOZ*SIN(THETAZ)
          PYPB = RHOZ*COS(THETAZ)
          IF (B1.EQ.0.0D0) THEN
              CCYB = -B2/YB
          ELSE IF (B2.EQ.0.0D0) THEN
              CCYB = B1/PYPB
          ELSE
              CCYB = -B2/YB
          END IF
      END IF
      ENDB = (NEND.EQ.2 .OR. NIVP.EQ.2) .AND. (LCB .OR. WREGB)
      IF ((.NOT.EIGV.OR.PERIOD) .AND. ENDB) THEN
          NPTS = 6
          IF (NCB.EQ.4) NPTS = NT(2)
          SIG = 1.0D0
          DO 120 I = 2,NPTS
              T = TT(I,2)
              PHI = YY(I,1,2)
              GEE = YY(I,3,2)
              SIG = EXP(GEE)
              IF (WREGB .AND. I.EQ.2) THEN
                  YB = SIG*SIN(PHI)
                  PYPB = SIG*COS(PHI)
                  IF (B1.EQ.0.0D0) THEN
                      CCYB = -B2/YB
                  ELSE IF (B2.EQ.0.0D0) THEN
                      CCYB = B1/PYPB
                  ELSE
                      CCYB = -B2/YB
                  END IF
              END IF
              IF (LCB) THEN
                  CALL DXDT(T,TMP,X)
                  CALL UV(X,U,PUP,V,PVP,HU,HV)
                  DD = U*PVP - V*PUP
                  BRYU = -DD*SIN(PHI)*SIG
                  BRYV = DD*COS(PHI)*SIG
                  IF (I.EQ.2) THEN
                      BRYUB = BRYU
                      BRYVB = BRYV
                      IF (B1.EQ.0.0D0) THEN
                          CCYB = -B2/BRYUB
                      ELSE IF (B2.EQ.0.0D0) THEN
                          CCYB = B1/BRYVB
                      ELSE
                          CCYB = -B2/BRYUB
                      END IF
                  END IF
                  BRYU = BRYU*CCYB
                  BRYV = BRYV*CCYB
              END IF
  120     CONTINUE
      END IF
C
      ENDA = (NEND.EQ.1 .OR. NIVP.EQ.2) .AND. NCA .GE. 3
      ENDB = (NEND.EQ.2 .OR. NIVP.EQ.2) .AND. NCB .GE. 3
C
C     RENORMALIZE.
C
      IF (NCA.EQ.3) THEN
          IF (A1.EQ.0.0D0) THEN
              IF (A2.GT.0.0D0) CCYA = -CCYA
          ELSE IF (A2.EQ.0.0D0) THEN
              IF (A1.LT.0.0D0) CCYA = -CCYA
          ELSE
              IF (A2.GT.0.0D0) CCYA = -CCYA
          END IF
      END IF
      IF (NCB.EQ.3) THEN
          IF (B1.EQ.0.0D0) THEN
              IF (B2.LT.0.0D0) CCYB = -CCYB
          ELSE IF (B2.EQ.0.0D0) THEN
              IF (B1.LT.0.0D0) CCYB = -CCYB
          ELSE
              IF (B2.LT.0.0D0) CCYB = -CCYB
          END IF
      END IF
      RETURN
      END
C
      SUBROUTINE DRAW(A1,A2,B1,B2,NUMEIG,EIG,SLFUN,NIVP,NEND,EIGV,NCA,
     +                NCB,ISLFN,XT,PLOTF,K11,K12,K21,K22,PERIOD)
C     **********
C  THIS PROGRAM PROVIDES POINTS OF EITHER AN EIGENFUNCTION
C  OR THE SOLUTION OF AN INITIAL VALUE PROBLEM FOR PLOTTING.
C  IT IS CALLED FROM "DRIVE".
C
C  INPUT QUANTITIES: A1,A2,B1,B2,NUMEIG,EIG,SLFUN,NIVP,NEND,
C                    EIGV,NCA,NCB,K11,K12,K21,K22,PERIOD
C  OUTPUT QUANTITIES: ISLFN,XT,PLOTF
C
C  IT FIRST OBTAINS A SET OF MESH POINTS T IN (-1,1)
C  WHICH ARE DEEMED SUITABLE FOR PLOTTING THE FUNCTION
C  WANTED, AND THEN OBTAINS THE VALUES OF THE FUNCTION AT
C  THOSE POINTS.  THE RESULTS FOR PLOTTING ARE STORED IN
C  ARRAYS XT AND PLOTF.
C     .. Scalars in Common ..
      REAL AA,BB,DTHDAA,DTHDBB,EIGSAV,HPI,PI,TMID,TWOPI
      INTEGER IND,MDTHZ,T21
      LOGICAL ADDD,PR
C     ..
C     .. Local Scalars ..
      REAL AXJMP,BRYU,BRYV,CCY,CCYA,CCYB,FZ,HUI,HVI,PUPI,
     +                 PVPI,PYP,RHO,RHOZ,TH,THETAZ,TI,TMP,UI,VI,XI,XJMP,
     +                 XJMPS,Y
      INTEGER I,JJ,KFLAG,MM
      LOGICAL ENDA,ENDB,LCOA,LCOB,REGA,REGB
C     ..
C     .. Local Arrays ..
      REAL SLFN(1000,2)
C     ..
C     .. External Subroutines ..
      EXTERNAL DXDT,EIGENF,MESH,PERFUN,RENORM,SCREEN,UV
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS,COS,EXP,INT,SIN,SQRT
C     ..
C     .. Common blocks ..
C
      COMMON /DATAF/EIGSAV,IND
      COMMON /PIE/PI,TWOPI,HPI
      COMMON /PRIN/PR,T21
      COMMON /TDATA/AA,TMID,BB,DTHDAA,DTHDBB,MDTHZ,ADDD
C     ..
C     DEFINITION OF SOME LOGICALS.
C
C     .. Scalar Arguments ..
      REAL A1,A2,B1,B2,EIG,K11,K12,K21,K22
      INTEGER ISLFN,NCA,NCB,NEND,NIVP,NUMEIG
      LOGICAL EIGV,PERIOD
C     ..
C     .. Array Arguments ..
      REAL PLOTF(1000,6),SLFUN(9),XT(1000,2)
C     ..
      REGA = NCA .EQ. 1
      REGB = NCB .EQ. 1
      LCOA = NCA .EQ. 4
      LCOB = NCB .EQ. 4
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
C     THE RELEVANT VALUES HAVE BEEN STORED IN ARRAY SLFUN, AND
C     NEED TO BE COPIED INTO THE FIRST COLUMN OF ARRAY SLFN.
C
C  THE FIRST 9 POSITIONS IN SLFUN() ARE RESERVED FOR
C  INITIAL DATA.  THEY ARE NOT USED FOR EIGENFUNCTION VALUES.
C
      IF (PR) WRITE (T21,FMT=*) ' SLFUN = '
      DO 10 I = 1,9
          SLFN(I,1) = SLFUN(I)
          IF (PR) WRITE (T21,FMT=*) SLFUN(I)
   10 CONTINUE

C ---------------------------------------------------------C
C
      CALL MESH(EIG,SLFN(10,1),ISLFN)
C
C     THE POINTS GENERATED BY MESH MAY NOT BE WITHIN THE INTERVAL
C     (AA,BB).  WE WANT TO ENSURE THAT WE USE ONLY SUCH POINTS
C     AS ARE WITHIN THE INTERVAL.
C
C     IN ADDITION, IF THE ENDPOINT IS NOT REG, WE WILL NOT TRY TO
C     PLOT POINTS WITH T.LT.-.95 OR T.GT..95 (IT COULD CAUSE TROUBLE
C     IN TRYING TO INTEGRATE FROM TOO NEAR A SINGULAR POINT).
C     SUBROUTINE SCREEN REMOVES POINTS GENERATED BY SUBROUTINE
C     MESH WHICH PROBABLY WOULD BE A PROBLEM.
C
C
      CALL SCREEN(AA,BB,NCA,NCB,NIVP,NEND,ISLFN,SLFN)
C
C  THE ENTRIES IN SLFN(9+I,1) INITIALLY ARE THE POINTS T IN (-1,1)
C     WHERE THE FUNCTION IS TO BE EVALUATED.
C
      DO 80 I = 1,ISLFN
          XT(9+I,2) = SLFN(9+I,1)
          IF (PR) WRITE (T21,FMT=*) ' I,XT = ',I,XT(9+I,2)
   80 CONTINUE
C
C     WARNING: THE VALUES RETURNED IN SLFN BY EIGENF DEPEND
C              ON THE VALUE OF KFLAG IN THE CALL TO EIGENF:
C
C              IF KFLAG = 1, THE VALUES IN SLFN(9+I,1) ARE THE
C                            EIGENFUNCTION Y ITSELF.
C
C              IF KFLAG = 2, THE VALUES IN THE TWO COLUMNS OF SLFN ARE
C                               SLFN(9+I,1) = THETA(9+I)
C                               SLFN(9+I,2) = EFF(9+I)
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

c     IF(PR)WRITE(T21,*) ' A1,A2,B1,B2 = ',A1,A2,B1,B2
c     IF(PR)WRITE(T21,*) ' REGA,NCA,REGB,NCB = ',REGA,NCA,REGB,NCB
c     IF(PR)WRITE(T21,*) ' ISLFN,KFLAG,EIGSAV = ',ISLFN,KFLAG,EIGSAV
c     IF(PR)WRITE(T21,*) ' THA,THB = ',SLFN(3,1),SLFN(6,1)


      IF (PERIOD) THEN
          CALL PERFUN(AA,BB,TMID,NCA,NCB,K11,K12,K21,K22,EIG,ISLFN,XT,
     +                PLOTF)

      ELSE
          CALL EIGENF(NUMEIG,A1,A2,B1,B2,REGA,NCA,REGB,NCB,SLFN,ISLFN,
     +                KFLAG)

          DO 113 I = 1,ISLFN
              IF (PR) WRITE (T21,FMT=*) I,SLFN(9+I,1),SLFN(9+I,2)
  113     CONTINUE
C
C     IT MAY HAPPEN THAT THETA(I) HAS A JUMP OF MM*PI AT TMID.
C     IN THIS CASE, WE NEED TO SUBTRACT MM*PI FROM THETA(I).
C

          IF (EIGV .AND. (LCOA.OR.LCOB)) THEN
              JJ = 0
              XJMPS = 2.5D0
              DO 90 I = 1,ISLFN - 1
                  XJMP = SLFN(10+I,1) - SLFN(9+I,1)
                  IF (ABS(XJMP).GT.XJMPS) THEN
                      XJMPS = ABS(XJMP)
                      JJ = I
                  END IF
   90         CONTINUE
              IF (JJ.NE.0) THEN
                  XJMP = SLFN(10+JJ,1) - SLFN(9+JJ,1)
                  AXJMP = ABS(XJMP)
                  MM = INT(AXJMP/PI)
                  IF (AXJMP-MM*PI.GT.HPI) MM = MM + 1
                  IF (XJMP.LT.0.0D0) MM = -MM
                  DO 100 I = JJ,ISLFN - 1
                      SLFN(10+I,1) = SLFN(10+I,1) - MM*PI
  100             CONTINUE
              END IF
          END IF
C
C  SLEIGN NORMALIZES THE WAVEFUNCTION TO HAVE L2-NORM 1.0, BUT FOR AN
C  INITIAL VALUE PROBLEM WE WANT TO HAVE THE VALUE (Y) AND SLOPE (P*Y')
C  (OR [Y,U] AND [Y,V]) AT THE END TO BE THOSE SPECIFIED FOR THE
C  INITIAL CONDITIONS.  SO HERE WE MUST RE-NORMALIZE FOR THIS PURPOSE.
C  N.B. A1 = ALFA2 & A2 = -ALFA1

          CALL RENORM(EIGV,NEND,NIVP,A1,A2,B1,B2,NCA,NCB,ISLFN,SLFN,
     +                CCYA,CCYB,ENDA,ENDB,PERIOD)

C  THE MAIN OUTPUT OF SUBROUTINE RENORM IS THE QUANTITIES
C    CCYA AND CCYB.

          CCY = 1.0D0
          DO 140 I = 1,ISLFN
              TI = XT(9+I,2)
              IF (.NOT.EIGV) THEN
                  IF (TI.LE.TMID) THEN
                      CCY = CCYA
                  ELSE
                      CCY = CCYB
                  END IF
              END IF

              CALL DXDT(TI,TMP,XI)
              XT(9+I,1) = XI
              THETAZ = SLFN(9+I,1)
              FZ = SLFN(9+I,2)
              RHOZ = EXP(FZ)
              Y = RHOZ*SIN(THETAZ)
              PYP = RHOZ*COS(THETAZ)
              Y = Y*CCY
              PYP = PYP*CCY
              RHO = SQRT(Y**2+PYP**2)
              IF ((.NOT.ENDA.AND.TI.LE.TMID) .OR.
     +            (.NOT.ENDB.AND.TI.GT.TMID)) THEN
              END IF
              PLOTF(9+I,1) = Y
              PLOTF(9+I,2) = PYP
              PLOTF(9+I,3) = Y
              PLOTF(9+I,4) = PYP
              TH = THETAZ
              PLOTF(9+I,5) = TH
              PLOTF(9+I,6) = RHO
C
              IF ((ENDA.AND.TI.LE.TMID) .OR. (ENDB.AND.TI.GT.TMID)) THEN
                  CALL UV(XI,UI,PUPI,VI,PVPI,HUI,HVI)
                  BRYU = PUPI*Y - PYP*UI
                  BRYV = PVPI*Y - PYP*VI
C
                  PLOTF(9+I,3) = BRYU
                  PLOTF(9+I,4) = BRYV
              END IF
  140     CONTINUE

C
C        IN THE CASE OF NIVP = 2  WE HAVE COMPUTED
C        THE SOLUTIONS TO TWO INITIAL VALUE PROBLEMS.
C        FROM THE TWO ENDS.  IF (.NOT.ENDA .AND. .NOT.ENDB)
C        WE SHOULD NOW HAVE
C
C        Y(A) = -A2 (ALFA1)  ;  Y(B) = -B2 (BETA1)
C        PY'(A) = A1 (ALFA2)  ;  PY'(B) = B1 (BETA2)
C
C  ------------------------------------------------------C
      END IF

      RETURN
      END
C
      SUBROUTINE EIGENF(NUMEIG,A1,A2,B1,B2,AOK,NCA,BOK,NCB,SLFN,ISLFN,
     +                  KFLAG)
C     **********
C  THIS PROGRAM CALCULATES SELECTED EIGENFUNCTION VALUES BY
C  INTEGRATION (OVER T IN (-1,1)).  IT IS CALLED FROM DRAW.
C
C  INPUT QUANTITIES: NUMEIG,A1,A2,B1,B2,AOK,NCA,BOK,NCB,SLFN,ISLFN,
C                    KFLAG,AA,TMID,BB,DTHDAA,DTHDBB
C  OUTPUT QUANTITIES: SLFN
C
C  N.B.: IT IS ASSUMED THAT THE POINTS T (=SLFN(9+I,1)) IN THE ARRAY
C     SLFN ALL LIE WITHIN THE INTERVAL (AA,BB).
C
C     WARNING:  THE ARRAY SLFN HERE IS TWO-DIMENSIONAL,
C            WHEREAS IT IS ONE-DIMENSIONAL IN SUBROUTINE SLEIGN.
C     **********
C     .. Scalars in Common ..
      REAL AA,BB,DTHDAA,DTHDBB,HPI,PI,TMID,TSAVEL,TSAVER,
     +                 TWOPI,Z
      INTEGER ISAVE,MDTHZ,T21
      LOGICAL ADDD,PR
C     ..
C     .. Local Scalars ..
      REAL DTHDAT,DTHDBT,DTHDET,EFF,EIGPI,T,THT,TM
      INTEGER I,IFLAG,J,NC,NMID
      LOGICAL LCOA,LCOB,OK
C     ..
C     .. Local Arrays ..
      REAL ERL(3),ERR(3),YL(3),YR(3)
C     ..
C     .. External Subroutines ..
      EXTERNAL INTEG
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC EXP,SIN
C     ..
C     .. Common blocks ..
      COMMON /PIE/PI,TWOPI,HPI
      COMMON /PRIN/PR,T21
      COMMON /TDATA/AA,TMID,BB,DTHDAA,DTHDBB,MDTHZ,ADDD
      COMMON /TSAVE/TSAVEL,TSAVER,ISAVE
      COMMON /Z1/Z
C     ..
C
C     .. Scalar Arguments ..
      REAL A1,A2,B1,B2
      INTEGER ISLFN,KFLAG,NCA,NCB,NUMEIG
      LOGICAL AOK,BOK
C     ..
C     .. Array Arguments ..
      REAL SLFN(1000,2)
C     ..
      Z = 1.0D0
      EIGPI = NUMEIG*PI
      NMID = 0
      DO 10 I = 1,ISLFN
          IF (SLFN(9+I,1).LE.TMID) NMID = I
   10 CONTINUE
      IF (NMID.GT.0) THEN
          LCOA = NCA .EQ. 4
          T = AA
          YL(1) = SLFN(3,1)
          YL(2) = 0.0D0
          YL(3) = 0.0D0
          OK = AOK
          NC = NCA
          EFF = 0.0D0
          DO 20 J = 1,NMID
              TM = SLFN(J+9,1)
              IF (TM.LT.AA .OR. TM.GT.BB) THEN
                  IF (PR) WRITE (*,FMT=*) ' T.LT.AA .OR. T.GT.BB '
                  STOP
              END IF
              THT = YL(1)
              DTHDAT = DTHDAA*EXP(-2.0D0*EFF)
              DTHDET = YL(2)
              IF (TM.GT.AA) THEN
                  ISAVE = 0
                  CALL INTEG(T,THT,DTHDAT,DTHDET,TM,A1,A2,SLFN(8,1),YL,
     +                       ERL,OK,NC,IFLAG)
                  IF (NC.EQ.4) THEN
                      EFF = YL(3)
                  ELSE
                      NC = 1
                      OK = .TRUE.
                      EFF = EFF + YL(3)
                  END IF
              END IF
              IF (KFLAG.EQ.1) THEN
                  SLFN(J+9,1) = SIN(YL(1))*EXP(EFF+SLFN(4,1))
              ELSE
                  SLFN(J+9,1) = YL(1)
                  SLFN(J+9,2) = EFF + SLFN(4,1)
              END IF
              T = TM
              IF (T.GT.-1.0D0) THEN
                  OK = .TRUE.
                  NC = 1
              END IF
              IF (T.LT.-0.9D0 .AND. LCOA) THEN
C  IN THIS CASE, INTEGRATE FROM AA AGAIN, AS AN OSC PROBLEM.
                  NC = 4
                  OK = .FALSE.
                  T = AA
                  YL(1) = SLFN(3,1)
                  YL(2) = 0.0D0
                  YL(3) = 0.0D0
              END IF
   20     CONTINUE
      END IF
      EFF = 0.0D0
      IF (NMID.LT.ISLFN) THEN
          LCOB = NCB .EQ. 4
          T = BB
          YR(1) = SLFN(6,1) - EIGPI
          YR(2) = 0.0D0
          YR(3) = 0.0D0
          OK = BOK
          NC = NCB
          EFF = 0.0D0
          DO 30 J = ISLFN,NMID + 1,-1
              TM = SLFN(J+9,1)
              IF (TM.LT.AA .OR. TM.GT.BB) THEN
                  IF (PR) WRITE (*,FMT=*) ' T.LT.AA .OR. T.GT.BB '
                  STOP
              END IF
              THT = YR(1)
              DTHDBT = DTHDBB*EXP(-2.0D0*EFF)
              DTHDET = YR(2)
              IF (TM.LT.BB) THEN
                  ISAVE = 0
                  CALL INTEG(T,THT,DTHDBT,DTHDET,TM,B1,B2,SLFN(8,1),YR,
     +                       ERR,OK,NC,IFLAG)
                  IF (NC.EQ.4) THEN
                      EFF = YR(3)
                  ELSE
                      NC = 1
                      OK = .TRUE.
                      EFF = EFF + YR(3)
                  END IF
              END IF
              IF (KFLAG.EQ.1) THEN
                  SLFN(J+9,1) = SIN(YR(1)+EIGPI)*EXP(EFF+SLFN(7,1))
              ELSE
                  SLFN(J+9,1) = YR(1) + EIGPI
                  SLFN(J+9,2) = EFF + SLFN(7,1)
              END IF
              T = TM
              IF (T.LT.1.0D0) THEN
                  OK = .TRUE.
                  NC = 1
              END IF
              IF (T.GT.0.9D0 .AND. LCOB) THEN
C  IN THIS CASE, INTEGRATE FROM BB AGAIN, AS AN OSC PROBLEM.
                  NC = 4
                  OK = .FALSE.
                  T = BB
                  YR(1) = SLFN(6,1) - EIGPI
                  YR(2) = 0.0D0
                  YR(3) = 0.0D0
              END IF
   30     CONTINUE
      END IF
      RETURN
      END
C
      SUBROUTINE GERK(F,NEQN,Y,T,TOUT,RELERR,ABSERR,IFLAG,GERROR,WORK,
     +                IWORK)
C
C     FEHLBERG FOURTH(FIFTH) ORDER RUNGE-KUTTA METHOD WITH
C     GLOBAL ERROR ASSESSMENT
C
C     WRITTEN BY H.A.WATTS AND L.F.SHAMPINE
C                   SANDIA LABORATORIES
C
C    GERK IS DESIGNED TO SOLVE SYSTEMS OF DIFFERENTIAL EQUATIONS
C    WHEN IT IS IMPORTANT TO HAVE A READILY AVAILABLE GLOBAL ERROR
C    ESTIMATE.  PARALLEL INTEGRATION IS PERFORMED TO YIELD TWO
C    SOLUTIONS ON DIFFERENT MESH SPACINGS AND GLOBAL EXTRAPOLATION
C    IS APPLIED TO PROVIDE AN ESTIMATE OF THE GLOBAL ERROR IN THE
C    MORE ACCURATE SOLUTION.
C
C    FOR IBM SYSTEM 360 AND 370 AND OTHER MACHINES OF SIMILAR
C    ARITHMETIC CHARACTERISTICS, THIS CODE SHOULD BE CONVERTED TO
C    REAL.
C
C*******************************************************************
C ABSTRACT
C*******************************************************************
C
C    SUBROUTINE  GERK  INTEGRATES A SYSTEM OF NEQN FIRST ORDER
C    ORDINARY DIFFERENTIAL EQUATIONS OF THE FORM
C             DY(I)/DT = F(T,Y(1),Y(2),...,Y(NEQN))
C              WHERE THE Y(I) ARE GIVEN AT T .
C    TYPICALLY THE SUBROUTINE IS USED TO INTEGRATE FROM T TO TOUT
C    BUT IT CAN BE USED AS A ONE-STEP INTEGRATOR TO ADVANCE THE
C    SOLUTION A SINGLE STEP IN THE DIRECTION OF TOUT. ON RETURN, AN
C    ESTIMATE OF THE GLOBAL ERROR IN THE SOLUTION AT T IS PROVIDED
C    AND THE PARAMETERS IN THE CALL LIST ARE SET FOR CONTINUING THE
C    INTEGRATION. THE USER HAS ONLY TO CALL GERK AGAIN (AND PERHAPS
C    DEFINE A NEW VALUE FOR TOUT). ACTUALLY, GERK  IS MERELY AN
C    INTERFACING ROUTINE WHICH ALLOCATES VIRTUAL STORAGE IN THE
C    ARRAYS WORK, IWORK AND CALLS SUBROUTINE GERKS FOR THE SOLUTION.
C    GERKS IN TURN CALLS SUBROUTINE FEHL WHICH COMPUTES AN APPROX-
C    IMATE SOLUTION OVER ONE STEP.
C
C    GERK  USES THE RUNGE-KUTTA-FEHLBERG (4,5)  METHOD DESCRIBED
C    IN THE REFERENCE
C    E.FEHLBERG , LOW-ORDER CLASSICAL RUNGE-KUTTA FORMULAS WITH
C                 STEPSIZE CONTROL , NASA TR R-315
C
C
C    THE PARAMETERS REPRESENT-
C      F -- SUBROUTINE F(T,Y,YP) TO EVALUATE DERIVATIVES
C           YP(I)=DY(I)/DT
C      NEQN -- NUMBER OF EQUATIONS TO BE INTEGRATED
C      Y(*) -- SOLUTION VECTOR AT T
C      T -- INDEPENDENT VARIABLE
C      TOUT -- OUTPUT POINT AT WHICH SOLUTION IS DESIRED
C      RELERR,ABSERR -- RELATIVE AND ABSOLUTE ERROR TOLERANCES FOR
C            LOCAL ERROR TEST. AT EACH STEP THE CODE REQUIRES THAT
C                 ABS(LOCAL ERROR) .LE. RELERR*ABS(Y) + ABSERR
C            FOR EACH COMPONENT OF THE LOCAL ERROR AND SOLUTION
C            VECTORS.
C      IFLAG -- INDICATOR FOR STATUS OF INTEGRATION
C      GERROR(*) -- VECTOR WHICH ESTIMATES THE GLOBAL ERROR AT T.
C                   THAT IS, GERROR(I) APPROXIMATES  Y(I)-TRUE
C                   SOLUTION(I).
C      WORK(*) -- ARRAY TO HOLD INFORMATION INTERNAL TO GERK WHICH
C            IS NECESSARY FOR SUBSEQUENT CALLS. MUST BE DIMENSIONED
C            AT LEAST  3+8*NEQN
C      IWORK(*) -- INTEGER ARRAY USED TO HOLD INFORMATION INTERNAL
C            TO GERK WHICH IS NECESSARY FOR SUBSEQUENT CALLS. MUST
C            BE DIMENSIONED AT LEAST  5
C
C
C*******************************************************************
C  FIRST CALL TO GERK
C*******************************************************************
C
C    THE USER MUST PROVIDE STORAGE IN HIS CALLING PROGRAM FOR THE
C    ARRAYS IN THE CALL LIST  -   Y(NEQN), WORK(3+8*NEQN), IWORK(5),
C    DECLARE F IN AN EXTERNAL STATEMENT, SUPPLY SUBROUTINE F(T,Y,YP)
C    AND INITIALIZE THE FOLLOWING PARAMETERS-
C
C      NEQN -- NUMBER OF EQUATIONS TO BE INTEGRATED.  (NEQN .GE. 1)
C      Y(*) -- VECTOR OF INITIAL CONDITIONS
C      T -- STARTING POINT OF INTEGRATION , MUST BE A VARIABLE
C      TOUT -- OUTPUT POINT AT WHICH SOLUTION IS DESIRED.
C            T=TOUT IS ALLOWED ON THE FIRST CALL ONLY,IN WHICH CASE
C            GERK RETURNS WITH IFLAG=2 IF CONTINUATION IS POSSIBLE.
C      RELERR,ABSERR -- RELATIVE AND ABSOLUTE LOCAL ERROR TOLERANCES
C            WHICH MUST BE NON-NEGATIVE BUT MAY BE CONSTANTS. WE CAN
C            USUALLY EXPECT THE GLOBAL ERRORS TO BE SOMEWHAT SMALLER
C            THAN THE REQUESTED LOCAL ERROR TOLERANCES. TO AVOID
C            LIMITING PRECISION DIFFICULTIES THE CODE ALWAYS USES
C            THE LARGER OF RELERR AND AN INTERNAL RELATIVE ERROR
C            PARAMETER WHICH IS MACHINE DEPENDENT.
C      IFLAG -- +1,-1  INDICATOR TO INITIALIZE THE CODE FOR EACH NEW
C            PROBLEM. NORMAL INPUT IS +1. THE USER SHOULD SET IFLAG=
C            -1 ONLY WHEN ONE-STEP INTEGRATOR CONTROL IS ESSENTIAL.
C            IN THIS CASE, GERK ATTEMPTS TO ADVANCE THE SOLUTION A
C            SINGLE STEP IN THE DIRECTION OF TOUT EACH TIME IT IS
C            CALLED. SINCE THIS MODE OF OPERATION RESULTS IN EXTRA
C            COMPUTING OVERHEAD, IT SHOULD BE AVOIDED UNLESS NEEDED.
C
C
C*******************************************************************
C  OUTPUT FROM GERK
C*******************************************************************
C
C      Y(*) -- SOLUTION AT T
C      T -- LAST POINT REACHED IN INTEGRATION.
C      IFLAG = 2 -- INTEGRATION REACHED TOUT.  INDICATES SUCCESSFUL
C                   RETURN AND IS THE NORMAL MODE FOR CONTINUING
C                   INTEGRATION.
C            =-2 -- A SINGLE SUCCESSFUL STEP IN THE DIRECTION OF
C                   TOUT HAS BEEN TAKEN. NORMAL MODE FOR CONTINUING
C                   INTEGRATION ONE STEP AT A TIME.
C            = 3 -- INTEGRATION WAS NOT COMPLETED BECAUSE MORE THAN
C                   9000 DERIVATIVE EVALUATIONS WERE NEEDED. THIS
C                   IS APPROXIMATELY 500 STEPS.
C            = 4 -- INTEGRATION WAS NOT COMPLETED BECAUSE SOLUTION
C                   VANISHED MAKING A PURE RELATIVE ERROR TEST
C                   IMPOSSIBLE. MUST USE NON-ZERO ABSERR TO CONTINUE.
C                   USING THE ONE-STEP INTEGRATION MODE FOR ONE STEP
C                   IS A GOOD WAY TO PROCEED.
C            = 5 -- INTEGRATION WAS NOT COMPLETED BECAUSE REQUESTED
C                   ACCURACY COULD NOT BE ACHIEVED USING SMALLEST
C                   ALLOWABLE STEPSIZE. USER MUST INCREASE THE ERROR
C                   TOLERANCE BEFORE CONTINUED INTEGRATION CAN BE
C                   ATTEMPTED.
C            = 6 -- GERK IS BEING USED INEFFICIENTLY IN SOLVING
C                   THIS PROBLEM. TOO MUCH OUTPUT IS RESTRICTING THE
C                   NATURAL STEPSIZE CHOICE. USE THE ONE-STEP
C                   INTEGRATOR MODE.
C            = 7 -- INVALID INPUT PARAMETERS
C                   THIS INDICATOR OCCURS IF ANY OF THE FOLLOWING IS
C                   SATISFIED - NEQN .LE. 0
C                               T=TOUT  AND  IFLAG .NE. +1 OR -1
C                               RELERR OR ABSERR .LT. 0.
C                               IFLAG .EQ. 0 OR .LT. -2 OR .GT. 7
C      GERROR(*) -- ESTIMATE OF THE GLOBAL ERROR IN THE SOLUTION AT T
C      WORK(*),IWORK(*) -- INFORMATION WHICH IS USUALLY OF NO
C                   INTEREST TO THE USER BUT NECESSARY FOR SUBSEQUENT
C                   CALLS.  WORK(1),...,WORK(NEQN) CONTAIN THE FIRST
C                   DERIVATIVES OF THE SOLUTION VECTOR Y AT T.
C                   WORK(NEQN+1) CONTAINS THE STEPSIZE H TO BE
C                   ATTEMPTED ON THE NEXT STEP.  IWORK(1) CONTAINS
C                   THE DERIVATIVE EVALUATION COUNTER.
C
C
C*******************************************************************
C  SUBSEQUENT CALLS TO GERK
C*******************************************************************
C
C    SUBROUTINE GERK RETURNS WITH ALL INFORMATION NEEDED TO CONTINUE
C    THE INTEGRATION. IF THE INTEGRATION REACHED TOUT, THE USER NEED
C    ONLY DEFINE A NEW TOUT AND CALL GERK AGAIN. IN THE ONE-STEP
C    INTEGRATOR MODE (IFLAG=-2) THE USER MUST KEEP IN MIND THAT EACH
C    STEP TAKEN IS IN THE DIRECTION OF THE CURRENT TOUT.  UPON
C    REACHING TOUT (INDICATED BY CHANGING IFLAG TO 2), THE USER MUST
C    THEN DEFINE A NEW TOUT AND RESET IFLAG TO -2 TO CONTINUE IN THE
C    ONE-STEP INTEGRATOR MODE.
C
C    IF THE INTEGRATION WAS NOT COMPLETED BUT THE USER STILL WANTS
C    TO CONTINUE (IFLAG=3 CASE), HE JUST CALLS GERK AGAIN.  THE
C    FUNCTION COUNTER IS THEN RESET TO 0 AND ANOTHER 9000 FUNCTION
C    EVALUATIONS ARE ALLOWED.
C
C    HOWEVER, IN THE CASE IFLAG=4, THE USER MUST FIRST ALTER THE
C    ERROR CRITERION TO USE A POSITIVE VALUE OF ABSERR BEFORE
C    INTEGRATION CAN PROCEED. IF HE DOES NOT,EXECUTION IS TERMINATED.
C
C    ALSO, IN THE CASE IFLAG=5, IT IS NECESSARY FOR THE USER TO
C    RESET IFLAG TO 2 (OR -2 WHEN THE ONE-STEP INTEGRATION MODE IS
C    BEING USED) AS WELL AS INCREASING EITHER ABSERR,RELERR OR BOTH
C    BEFORE THE INTEGRATION CAN BE CONTINUED. IF THIS IS NOT DONE,
C    EXECUTION WILL BE TERMINATED.  THE OCCURRENCE OF IFLAG=5
C    INDICATES A TROUBLE SPOT (SOLUTION IS CHANGING RAPIDLY,
C    SINGULARITY MAY BE PRESENT) AND IT OFTEN IS INADVISABLE TO
C    CONTINUE.
C
C    IF IFLAG=6 IS ENCOUNTERED, THE USER SHOULD USE THE ONE-STEP
C    INTEGRATION MODE WITH THE STEPSIZE DETERMINED BY THE CODE. IF
C    THE USER INSISTS UPON CONTINUING THE INTEGRATION WITH GERK IN
C    THE INTERVAL MODE, HE MUST RESET IFLAG TO 2 BEFORE CALLING GERK
C    AGAIN.  OTHERWISE,EXECUTION WILL BE TERMINATED.
C
C    IF IFLAG=7 IS OBTAINED, INTEGRATION CAN NOT BE CONTINUED UNLESS
C    THE INVALID INPUT PARAMETERS ARE CORRECTED.
C
C    IT SHOULD BE NOTED THAT THE ARRAYS WORK,IWORK CONTAIN
C    INFORMATION REQUIRED FOR SUBSEQUENT INTEGRATION. ACCORDINGLY,
C    WORK AND IWORK SHOULD NOT BE ALTERED.
C
C*******************************************************************
C
C     .. Scalar Arguments ..
      REAL ABSERR,RELERR,T,TOUT
      INTEGER IFLAG,NEQN
C     ..
C     .. Array Arguments ..
      REAL GERROR(NEQN),WORK(3+8*NEQN),Y(NEQN)
      INTEGER IWORK(5)
C     ..
C     .. Subroutine Arguments ..
      EXTERNAL F
C     ..
C     .. Local Scalars ..
      REAL H,SAVAE,SAVRE
      INTEGER I,IM,INIT,JFLAG,K1,K1M,K2,K3,K4,K5,K6,K7,K8,KFLAG,KOP,NFE
C     ..
C     .. External Subroutines ..
      EXTERNAL GERKS
C     ..
C COMPUTE INDICES FOR THE SPLITTING OF THE WORK ARRAY
C     .. Local Arrays ..
      REAL F1(3),F2(3),F3(3),F4(3),F5(3),YG(3),YGP(3),YP(3)
C     ..
      K1M = NEQN + 1
      K1 = K1M + 1
      K2 = K1 + NEQN
      K3 = K2 + NEQN
      K4 = K3 + NEQN
      K5 = K4 + NEQN
      K6 = K5 + NEQN
      K7 = K6 + NEQN
      K8 = K7 + NEQN
C  THE FOLLOWING SECTION DEFINES LOCAL VARIABLES, F1,F2,...,SO
C  THAT GERKS CAN BE CALLED IN A WAY THAT IS MORE "PORTABLE"
C  THAN THE ORIGINAL ARRANGEMENT.  NOTE THE NEW ARGUMENT LIST
C  FOR GERKS.
      DO 13 I = 1,3
          IM = I - 1
          YP(I) = WORK(I)
          F1(I) = WORK(K1+IM)
          F2(I) = WORK(K2+IM)
          F3(I) = WORK(K3+IM)
          F4(I) = WORK(K4+IM)
          F5(I) = WORK(K5+IM)
          YG(I) = WORK(K6+IM)
          YGP(I) = WORK(K7+IM)
   13 CONTINUE
      H = WORK(K1M)
      SAVRE = WORK(K8)
      SAVAE = WORK(K8+1)
      NFE = IWORK(1)
      KOP = IWORK(2)
      INIT = IWORK(3)
      JFLAG = IWORK(4)
      KFLAG = IWORK(5)
C *******************************************************************
C      THIS INTERFACING ROUTINE MERELY RELIEVES THE USER OF A LONG
C      CALLING LIST VIA THE SPLITTING APART OF TWO WORKING STORAGE
C      ARRAYS. IF THIS IS NOT COMPATIBLE WITH THE USERS COMPILER,
C      HE MUST USE GERKS DIRECTLY.
C *******************************************************************
c     CALL GERKS(F,NEQN,Y,T,TOUT,RELERR,ABSERR,IFLAG,GERROR,WORK(1),
c    +           WORK(K1M),WORK(K1),WORK(K2),WORK(K3),WORK(K4),WORK(K5),
c    +           WORK(K6),WORK(K7),WORK(K8),WORK(K8+1),IWORK(1),
c    +           IWORK(2),IWORK(3),IWORK(4),IWORK(5))

      CALL GERKS(F,NEQN,Y,T,TOUT,RELERR,ABSERR,IFLAG,GERROR,YP,H,F1,F2,
     +           F3,F4,F5,YG,YGP,SAVRE,SAVAE,NFE,KOP,INIT,JFLAG,KFLAG)

C  NOW WE HAVE TO REPLACE THE LOCAL ARRAYS, F1,F2,..., WITH THEIR
C  EQUIVALENT LOCATIONS IN THE ARRAY WORK.
      DO 15 I = 1,3
          IM = I - 1
          WORK(I) = YP(I)
          WORK(K1+IM) = F1(I)
          WORK(K2+IM) = F2(I)
          WORK(K3+IM) = F3(I)
          WORK(K4+IM) = F4(I)
          WORK(K5+IM) = F5(I)
          WORK(K6+IM) = YG(I)
          WORK(K7+IM) = YGP(I)
   15 CONTINUE
      WORK(K1M) = H
      WORK(K8) = SAVRE
      WORK(K8+1) = SAVAE
      IWORK(1) = NFE
      IWORK(2) = KOP
      IWORK(3) = INIT
      IWORK(4) = JFLAG
      IWORK(5) = KFLAG
      RETURN
      END
      SUBROUTINE GERKS(F,NEQN,Y,T,TOUT,RELERR,ABSERR,IFLAG,GERROR,YP,H,
     +                 F1,F2,F3,F4,F5,YG,YGP,SAVRE,SAVAE,NFE,KOP,INIT,
     +                 JFLAG,KFLAG)
C      FEHLBERG FOURTH(FIFTH) ORDER RUNGE-KUTTA METHOD WITH
C      GLOBAL ERROR ASSESSMENT
C *******************************************************************
C      GERKS INTEGRATES A SYSTEM OF FIRST ORDER ORDINARY DIFFERENTIAL
C      EQUATIONS AS DESCRIBED IN THE COMMENTS FOR GERK.  THE ARRAYS
C      YP,F1,F2,F3,F4,F5,YG AND YGP (OF DIMENSION AT LEAST NEQN) AND
C      THE VARIABLES H,SAVRE,SAVAE,NFE,KOP,INIT,JFLAG,AND KFLAG ARE
C      USED INTERNALLY BY THE CODE AND APPEAR IN THE CALL LIST TO
C      ELIMINATE LOCAL RETENTION OF VARIABLES BETWEEN CALLS.
C      ACCORDINGLY, THEY SHOULD NOT BE ALTERED. ITEMS OF POSSIBLE
C      INTEREST ARE
C          YP - DERIVATIVE OF SOLUTION VECTOR AT T
C          H  - AN APPROPRIATE STEPSIZE TO BE USED FOR THE NEXT STEP
C          NFE- COUNTER ON THE NUMBER OF DERIVATIVE FUNCTION
C               EVALUATIONS.
C *******************************************************************
C     .. Scalar Arguments ..
      REAL ABSERR,H,RELERR,SAVAE,SAVRE,T,TOUT
      INTEGER IFLAG,INIT,JFLAG,KFLAG,KOP,NEQN,NFE
C     ..
C     .. Array Arguments ..
      REAL F1(NEQN),F2(NEQN),F3(NEQN),F4(NEQN),F5(NEQN),
     +                 GERROR(NEQN),Y(NEQN),YG(NEQN),YGP(NEQN),YP(NEQN)
C     ..
C     .. Subroutine Arguments ..
      EXTERNAL F
C     ..
C     .. Local Scalars ..
      REAL A,AE,DT,EE,EEOET,ESTTOL,ET,HH,HMIN,ONE,REMIN,RER,
     +                 S,SCALE,TMP,TOL,TOLN,TS,U,U26,YPK
      INTEGER K,MAXNFE,MFLAG
      LOGICAL HFAILD,OUTPUT
C     ..
C     .. External Functions ..
      REAL EPSLON
      EXTERNAL EPSLON
C     ..
C     .. External Subroutines ..
      EXTERNAL FEHL
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS,MAX,MIN,SIGN
C     ..
C *******************************************************************
C   REMIN IS A TOLERANCE THRESHOLD WHICH IS ALSO DETERMINED BY THE
C   INTEGRATION METHOD. IN PARTICULAR, A FIFTH ORDER METHOD WILL
C   GENERALLY NOT BE CAPABLE OF DELIVERING ACCURACIES NEAR LIMITING
C   PRECISION ON COMPUTERS WITH LONG WORDLENGTHS.
C *******************************************************************
C      THE EXPENSE IS CONTROLLED BY RESTRICTING THE NUMBER
C      OF FUNCTION EVALUATIONS TO BE APPROXIMATELY MAXNFE.
C      AS SET,THIS CORRESPONDS TO ABOUT 500 STEPS.
C *******************************************************************
C     U - THE COMPUTER UNIT ROUNDOFF ERROR U IS THE SMALLEST POSITIVE
C         VALUE REPRESENTABLE IN THE MACHINE SUCH THAT  1.+ U .GT. 1.
C     (VARIABLE ONE SET TO 1.0 EASES PRECISION CONVERSION.)
C
C     .. Data statements ..
      DATA REMIN/3.0D-11/
      DATA MAXNFE/9000/
C     ..
      ONE = 1.0D0
      U = EPSLON(ONE)
C *******************************************************************
C      CHECK INPUT PARAMETERS
      IF (NEQN.LT.1) GO TO 10
      IF ((RELERR.LT.0.D0) .OR. (ABSERR.LT.0.D0)) GO TO 10
      MFLAG = ABS(IFLAG)
      IF ((MFLAG.GE.1) .AND. (MFLAG.LE.7)) GO TO 20
C INVALID INPUT
   10 IFLAG = 7
      RETURN
C IS THIS THE FIRST CALL
   20 IF (MFLAG.EQ.1) GO TO 70
C CHECK CONTINUATION POSSIBILITIES
      IF (T.EQ.TOUT) GO TO 10
      IF (MFLAG.NE.2) GO TO 30
C IFLAG = +2 OR -2
      IF (INIT.EQ.0) GO TO 60
      IF (KFLAG.EQ.3) GO TO 50
      IF ((KFLAG.EQ.4) .AND. (ABSERR.EQ.0.D0)) GO TO 40
      IF ((KFLAG.EQ.5) .AND. (RELERR.LE.SAVRE) .AND.
     +    (ABSERR.LE.SAVAE)) GO TO 40
      GO TO 70
C IFLAG = 3,4,5,6, OR 7
   30 IF (IFLAG.EQ.3) GO TO 50
      IF ((IFLAG.EQ.4) .AND. (ABSERR.GT.0.D0)) GO TO 60
C INTEGRATION CANNOT BE CONTINUED SINCE USER DID NOT RESPOND TO
C THE INSTRUCTIONS PERTAINING TO IFLAG=4,5,6 OR 7
   40 STOP
C *******************************************************************
C      RESET FUNCTION EVALUATION COUNTER
   50 NFE = 0
      IF (MFLAG.EQ.2) GO TO 70
C RESET FLAG VALUE FROM PREVIOUS CALL
   60 IFLAG = JFLAG
C SAVE INPUT IFLAG AND SET CONTINUATION FLAG VALUE FOR SUBSEQUENT
C INPUT CHECKING
   70 JFLAG = IFLAG
      KFLAG = 0
C SAVE RELERR AND ABSERR FOR CHECKING INPUT ON SUBSEQUENT CALLS
      SAVRE = RELERR
      SAVAE = ABSERR
C RESTRICT RELATIVE ERROR TOLERANCE TO BE AT LEAST AS LARGE AS
C 32U+REMIN TO AVOID LIMITING PRECISION DIFFICULTIES ARISING
C FROM IMPOSSIBLE ACCURACY REQUESTS
      TMP = 32.0D0*U+REMIN
      RER = MAX(RELERR,TMP)
      U26 = 26.D0*U
      DT = TOUT - T
      IF (MFLAG.EQ.1) GO TO 80
      IF (INIT.EQ.0) GO TO 90
      GO TO 110
C *******************************************************************
C      INITIALIZATION --
C                        SET INITIALIZATION COMPLETION INDICATOR,INIT
C                        SET INDICATOR FOR TOO MANY OUTPUT POINTS,KOP
C                        EVALUATE INITIAL DERIVATIVES
C                        COPY INITIAL VALUES AND DERIVATIVES FOR THE
C                              PARALLEL SOLUTION
C                        SET COUNTER FOR FUNCTION EVALUATIONS,NFE
C                        ESTIMATE STARTING STEPSIZE
   80 INIT = 0
      KOP = 0
      A = T
      CALL F(A,Y,YP)
      NFE = 1
      IF (T.NE.TOUT) GO TO 90
      IFLAG = 2
      RETURN
   90 INIT = 1
      H = ABS(DT)
      TOLN = 0.D0
      DO 100 K = 1,NEQN
          YG(K) = Y(K)
          YGP(K) = YP(K)
          TOL = RER*ABS(Y(K)) + ABSERR
          IF (TOL.LE.0.D0) GO TO 100
          TOLN = TOL
          YPK = ABS(YP(K))
          IF (YPK*H**5.GT.TOL) H = (TOL/YPK)**0.2D0
  100 CONTINUE
      IF (TOLN.LE.0.D0) H = 0.D0
      H = MAX(H,U26*MAX(ABS(T),ABS(DT)))
C *******************************************************************
C      SET STEPSIZE FOR INTEGRATION IN THE DIRECTION FROM T TO TOUT
  110 H = SIGN(H,DT)
C TEST TO SEE IF GERK IS BEING SEVERELY IMPACTED BY TOO MANY
C OUTPUT POINTS
      IF (ABS(H).GT.2.D0*ABS(DT)) KOP = KOP + 1
      IF (KOP.NE.100) GO TO 120
      KOP = 0
      IFLAG = 6
      RETURN
  120 IF (ABS(DT).GT.U26*ABS(T)) GO TO 140
C IF TOO CLOSE TO OUTPUT POINT,EXTRAPOLATE AND RETURN
      DO 130 K = 1,NEQN
          YG(K) = YG(K) + DT*YGP(K)
          Y(K) = Y(K) + DT*YP(K)
  130 CONTINUE
      A = TOUT
      CALL F(A,YG,YGP)
      CALL F(A,Y,YP)
      NFE = NFE + 2
      GO TO 230
C INITIALIZE OUTPUT POINT INDICATOR
  140 OUTPUT = .FALSE.
C TO AVOID PREMATURE UNDERFLOW IN THE ERROR TOLERANCE FUNCTION,
C SCALE THE ERROR TOLERANCES
      SCALE = 2.D0/RER
      AE = SCALE*ABSERR
C *******************************************************************
C *******************************************************************
C      STEP BY STEP INTEGRATION
  150 HFAILD = .FALSE.
C SET SMALLEST ALLOWABLE STEPSIZE
      HMIN = U26*ABS(T)
C ADJUST STEPSIZE IF NECESSARY TO HIT THE OUTPUT POINT.
C LOOK AHEAD TWO STEPS TO AVOID DRASTIC CHANGES IN THE STEPSIZE
C AND THUS LESSEN THE IMPACT OF OUTPUT POINTS ON THE CODE.
      DT = TOUT - T
      IF (ABS(DT).GE.2.D0*ABS(H)) GO TO 170
      IF (ABS(DT).GT.ABS(H)) GO TO 160
C THE NEXT SUCCESSFUL STEP WILL COMPLETE THE INTEGRATION TO THE
C OUTPUT POINT
      OUTPUT = .TRUE.
      H = DT
      GO TO 170
  160 H = 0.5D0*DT
C *******************************************************************
C      CORE INTEGRATOR FOR TAKING A SINGLE STEP
C *******************************************************************
C      THE TOLERANCES HAVE BEEN SCALED TO AVOID PREMATURE UNDERFLOW
C      IN COMPUTING THE ERROR TOLERANCE FUNCTION ET.
C      TO AVOID PROBLEMS WITH ZERO CROSSINGS, RELATIVE ERROR IS
C      MEASURED USING THE AVERAGE OF THE MAGNITUDES OF THE SOLUTION
C      AT THE BEGINNING AND END OF A STEP.
C      THE ERROR ESTIMATE FORMULA HAS BEEN GROUPED TO CONTROL LOSS OF
C      SIGNIFICANCE.
C      TO DISTINGUISH THE VARIOUS ARGUMENTS, H IS NOT PERMITTED
C      TO BECOME SMALLER THAN 26 UNITS OF ROUNDOFF IN T.
C      PRACTICAL LIMITS ON THE CHANGE IN THE STEPSIZE ARE ENFORCED TO
C      SMOOTH THE STEPSIZE SELECTION PROCESS AND TO AVOID EXCESSIVE
C      CHATTERING ON PROBLEMS HAVING DISCONTINUITIES.
C      TO PREVENT UNNECESSARY FAILURES, THE CODE USES 9/10 THE
C      STEPSIZE IT ESTIMATES WILL SUCCEED.
C      AFTER A STEP FAILURE, THE STEPSIZE IS NOT ALLOWED TO INCREASE
C      FOR THE NEXT ATTEMPTED STEP.  THIS MAKES THE CODE MORE
C      EFFICIENT ON PROBLEMS HAVING DISCONTINUITIES AND MORE
C      EFFECTIVE IN GENERAL SINCE LOCAL EXTRAPOLATION IS BEING USED
C      AND THE ERROR ESTIMATE MAY BE UNRELIABLE OR UNACCEPTABLE WHEN
C      A STEP FAILS.
C *******************************************************************
C      TEST NUMBER OF DERIVATIVE FUNCTION EVALUATIONS.
C      IF OKAY,TRY TO ADVANCE THE INTEGRATION FROM T TO T+H
  170 IF (NFE.LE.MAXNFE) GO TO 180
C TOO MUCH WORK
      IFLAG = 3
      KFLAG = 3
      RETURN
C ADVANCE AN APPROXIMATE SOLUTION OVER ONE STEP OF LENGTH H
  180 CALL FEHL(F,NEQN,YG,T,H,YGP,F1,F2,F3,F4,F5,F1)
      NFE = NFE + 5
C COMPUTE AND TEST ALLOWABLE TOLERANCES VERSUS LOCAL ERROR
C ESTIMATES AND REMOVE SCALING OF TOLERANCES. NOTE THAT RELATIVE
C ERROR IS MEASURED WITH RESPECT TO THE AVERAGE MAGNITUDES OF THE
C OF THE SOLUTION AT THE BEGINNING AND END OF THE STEP.
      EEOET = 0.D0
      DO 200 K = 1,NEQN
          ET = ABS(YG(K)) + ABS(F1(K)) + AE
          IF (ET.GT.0.D0) GO TO 190
C INAPPROPRIATE ERROR TOLERANCE
          IFLAG = 4
          KFLAG = 4
          RETURN
  190     EE = ABS((-2090.D0*YGP(K)+ (21970.D0*F3(K)-15048.D0*F4(K)))+
     +         (22528.D0*F2(K)-27360.D0*F5(K)))
          EEOET = MAX(EEOET,EE/ET)
  200 CONTINUE
      ESTTOL = ABS(H)*EEOET*SCALE/752400.D0
      IF (ESTTOL.LE.1.D0) GO TO 210
C UNSUCCESSFUL STEP
C                   REDUCE THE STEPSIZE , TRY AGAIN
C                   THE DECREASE IS LIMITED TO A FACTOR OF 1/10
      HFAILD = .TRUE.
      OUTPUT = .FALSE.
      S = 0.1D0
      IF (ESTTOL.LT.59049.D0) S = 0.9D0/ESTTOL**0.2D0
      H = S*H
      IF (ABS(H).GT.HMIN) GO TO 170
C REQUESTED ERROR UNATTAINABLE AT SMALLEST ALLOWABLE STEPSIZE
      IFLAG = 5
      KFLAG = 5
      RETURN
C SUCCESSFUL STEP
C                    STORE ONE-STEP SOLUTION YG AT T+H
C                    AND EVALUATE DERIVATIVES THERE
  210 TS = T
      T = T + H
      DO 220 K = 1,NEQN
          YG(K) = F1(K)
  220 CONTINUE
      A = T
      CALL F(A,YG,YGP)
      NFE = NFE + 1
C NOW ADVANCE THE Y SOLUTION OVER TWO STEPS OF
C LENGTH H/2 AND EVALUATE DERIVATIVES THERE
      HH = 0.5D0*H
      CALL FEHL(F,NEQN,Y,TS,HH,YP,F1,F2,F3,F4,F5,Y)
      TS = TS + HH
      A = TS
      CALL F(A,Y,YP)
      CALL FEHL(F,NEQN,Y,TS,HH,YP,F1,F2,F3,F4,F5,Y)
      A = T
      CALL F(A,Y,YP)
      NFE = NFE + 12
C CHOOSE NEXT STEPSIZE
C THE INCREASE IS LIMITED TO A FACTOR OF 5
C IF STEP FAILURE HAS JUST OCCURRED, NEXT
C    STEPSIZE IS NOT ALLOWED TO INCREASE
      S = 5.D0
      IF (ESTTOL.GT.1.889568D-4) S = 0.9D0/ESTTOL**0.2D0
      IF (HFAILD) S = MIN(S,ONE)
      H = SIGN(MAX(S*ABS(H),HMIN),H)
C *******************************************************************
C      END OF CORE INTEGRATOR
C *******************************************************************
C      SHOULD WE TAKE ANOTHER STEP
      IF (OUTPUT) GO TO 230
      IF (IFLAG.GT.0) GO TO 150
C *******************************************************************
C *******************************************************************
C      INTEGRATION SUCCESSFULLY COMPLETED
C      ONE-STEP MODE
      IFLAG = -2
      GO TO 240
C INTERVAL MODE
  230 T = TOUT
      IFLAG = 2
  240 DO 250 K = 1,NEQN
          GERROR(K) = (YG(K)-Y(K))/31.D0
  250 CONTINUE
      RETURN
      END
      SUBROUTINE FEHL(F,NEQN,Y,T,H,YP,F1,F2,F3,F4,F5,S)
C      FEHLBERG FOURTH-FIFTH ORDER RUNGE-KUTTA METHOD
C *******************************************************************
C     FEHL INTEGRATES A SYSTEM OF NEQN FIRST ORDER
C     ORDINARY DIFFERENTIAL EQUATIONS OF THE FORM
C              DY(I)/DT=F(T,Y(1),---,Y(NEQN))
C     WHERE THE INITIAL VALUES Y(I) AND THE INITIAL DERIVATIVES
C     YP(I) ARE SPECIFIED AT THE STARTING POINT T. FEHL ADVANCES
C     THE SOLUTION OVER THE FIXED STEP H AND RETURNS
C     THE FIFTH ORDER (SIXTH ORDER ACCURATE LOCALLY) SOLUTION
C     APPROXIMATION AT T+H IN ARRAY S(I).
C     F1,---,F5 ARE ARRAYS OF DIMENSION NEQN WHICH ARE NEEDED
C     FOR INTERNAL STORAGE.
C     THE FORMULAS HAVE BEEN GROUPED TO CONTROL LOSS OF SIGNIFICANCE.
C     FEHL SHOULD BE CALLED WITH AN H NOT SMALLER THAN 13 UNITS OF
C     ROUNDOFF IN T SO THAT THE VARIOUS INDEPENDENT ARGUMENTS CAN BE
C     DISTINGUISHED.
C *******************************************************************
C     .. Scalar Arguments ..
      REAL H,T
      INTEGER NEQN
C     ..
C     .. Array Arguments ..
      REAL F1(NEQN),F2(NEQN),F3(NEQN),F4(NEQN),F5(NEQN),
     +                 S(NEQN),Y(NEQN),YP(NEQN)
C     ..
C     .. Subroutine Arguments ..
      EXTERNAL F
C     ..
C     .. Local Scalars ..
      REAL CH,T1
      INTEGER K
C     ..
      CH = 0.25D0*H
      DO 10 K = 1,NEQN
          F5(K) = Y(K) + CH*YP(K)
   10 CONTINUE
      T1 = T + 0.25D0*H
c     CALL F(T+0.25*H,F5,F1)
      CALL F(T1,F5,F1)
      CH = 0.09375D0*H
      DO 20 K = 1,NEQN
          F5(K) = Y(K) + CH* (YP(K)+3.D0*F1(K))
   20 CONTINUE
      T1 = T + 0.375D0*H
c     CALL F(T+0.375*H,F5,F2)
      CALL F(T1,F5,F2)
      CH = H/2197.D0
      DO 30 K = 1,NEQN
          F5(K) = Y(K) + CH* (1932.D0*YP(K)+
     +            (7296.D0*F2(K)-7200.D0*F1(K)))
   30 CONTINUE
      T1 = T + 12.D0/13.D0*H
c     CALL F(T+12./13.*H,F5,F3)
      CALL F(T1,F5,F3)
      CH = H/4104.D0
      DO 40 K = 1,NEQN
          F5(K) = Y(K) + CH* ((8341.D0*YP(K)-845.D0*F3(K))+
     +            (29440.D0*F2(K)-32832.D0*F1(K)))
   40 CONTINUE
      T1 = T + H
c     CALL F(T+H,F5,F4)
      CALL F(T1,F5,F4)
      CH = H/20520.D0
      DO 50 K = 1,NEQN
          F1(K) = Y(K) + CH* ((-6080.D0*YP(K)+ (9295.D0*F3(K)-
     +            5643.D0*F4(K)))+ (41040.D0*F1(K)-28352.D0*F2(K)))
   50 CONTINUE
      T1 = T + 0.5D0*H
c     CALL F(T+0.5*H,F1,F5)
      CALL F(T1,F1,F5)
C COMPUTE APPROXIMATE SOLUTION AT T+H
      CH = H/7618050.D0
      DO 60 K = 1,NEQN
          S(K) = Y(K) + CH* ((902880.D0*YP(K)+ (3855735.D0*F3(K)-
     +           1371249.D0*F4(K)))+ (3953664.D0*F2(K)+277020.D0*F5(K)))
   60 CONTINUE
      RETURN
      END
