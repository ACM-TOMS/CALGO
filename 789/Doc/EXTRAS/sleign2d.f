C     OCTOBER 15, 1995; P.B. BAILEY, W.N. EVERITT, B. GARBOW AND A. ZETTL
C
C     This program is for an equation of the form
C
C               -(p(x)*y'(x))' + q(x)*y(x) = eig*w(x)*y(x)
C
      SUBROUTINE SLEIGN2(A,B,INTAB,P0ATA,QFATA,P0ATB,QFATB,A1,A2,
     1                   B1,B2,NUMEIG,EIG,TOL,IFLAG,ISLFUN,SLFUN,
     2                   SINGATA,SINGATB,CIRCLA,CIRCLB,OSCILA,OSCILB)
      INTEGER INTAB,NUMEIG,IFLAG,ISLFUN
      DOUBLE PRECISION A,B,P0ATA,QFATA,P0ATB,QFATB,A1,A2,B1,B2,EIG,TOL,
     1     SINGATA,SINGATB,CIRCLA,CIRCLB,OSCILA,OSCILB
      DOUBLE PRECISION SLFUN(9)
C     **********
C
C     This subroutine is designed for the calculation of a specified
C     eigenvalue, EIG, of a Sturm-Liouville problem for the equation
C
C       -(p(x)*y'(x))' + q(x)*y(x) = eig*w(x)*y(x)  on (a,b)
C
C     with user-supplied coefficient functions p, q, and w.
C     The problem may be either nonsingular or singular.  In the
C     nonsingular case, boundary conditions of the form
C
C        A1*y(a) + A2*p(a)*y'(a) = 0
C        B1*y(b) + B2*p(b)*y'(b) = 0
C
C     are prescribed by specifying the numbers A1, A2, B1, and B2.
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
C       SUBROUTINE sleign2(a,b,intab,p0ata,qfata,p0atb,qfatb,a1,a2,
C                          b1,b2,numeig,eig,tol,iflag,islfun,slfun,
C                          singata,singatb,circla,circlb,oscila,oscilb)
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
C         IFLAG = 12 - failed to get a bracket.
C         IFLAG = 13 - AA cannot be moved in any further.
C         IFLAG = 14 - BB cannot be moved in any further.
C         IFLAG = 51 - integration failure after 1st call to INTEG.
C         IFLAG = 52 - integration failure after 2nd call to INTEG.
C         IFLAG = 53 - integration failure after 3rd call to INTEG.
C         IFLAG = 54 - integration failure after 4th call to INTEG.
C
C       ISLFUN is an integer input variable set to the number of
C         selected eigenfunction values desired.  If no values are
C         desired, set ISLFUN to zero.
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
C       SINGATA is an input variable set positive if endpoint A
C         is singular, and negative or zero if A is nonsingular.
C
C       SINGATB is an input variable set positive if endpoint B
C         is singular, and negative or zero if B is nonsingular.
C
C       CIRCLA is an input variable set positive if endpoint A has
C         a limit-circle singularity, and negative or zero if not.
C
C       CIRCLB is an input variable set positive if endpoint B has
C         a limit-circle singularity, and negative or zero if not.
C
C       OSCILA is an input variable set positive if the limit-circle
C         singularity at A is oscillatory, and negative or zero if not.
C
C       OSCILB is an input variable set positive if the limit-circle
C         singularity at B is oscillatory, and negative or zero if not.
C
C     Subprograms called
C
C       user-supplied ..... p,q,w,uv
C
C       sleign2-supplied .. aabb,alfbet,dxdt,eigenf,epslon,estpac,integ,
C                           setmid,tfromi,thum
C
C     This version dated 7/23/95.
C     Paul B. Bailey, Albuquerque, New Mexico
C     Burton S. Garbow, Park Forest, Illinois
C
C     **********
C     .. Scalars in Common ..
      INTEGER INTSAV,IND
      DOUBLE PRECISION ASAV,BSAV,C1,C2,EIGSAV,EPSMIN,Z,
     1     PI,TWOPI,HPI,TSAVEL,TSAVER
C     ..
C     .. Arrays in Common ..
      INTEGER JAY(100)
      DOUBLE PRECISION TEE(100),ZEE(100)
C     ..
C     .. Local Scalars ..
      INTEGER I,IA,IB,IMAX,IMIN,IOUT,JFLAG,KFLAG,MF,ML,
     1        NEIG,K,JJL,JJR,IE,IMID,NITER,NRAY,LOOP2,LOOP3,
     2        MDTHZ
      LOGICAL AOK,BOK,BRACKT,CONVRG,FYNYT,FYNYT1,LOGIC,
     1        NEWTON,ONEDIG,PRIN,THEGT0,THELT0,OLDNEWT,SINGA,
     2        SINGB,LCIRCA,LCIRCB,OSCA,OSCB,ENDA,ENDB,LIMUP,
     3        ADDD,EXIT,FIRSTT,NEWTONF,LIMA,LIMB,BRS,
     4        BRSS,CHNGEPS,TRUNKA,TRUNKB
      DOUBLE PRECISION AA,AAA,AAF,ALFA,BB,BBB,BBF,BETA,
     1     C,CHNG,CL,CR,DE,DEDW,DEN,DERIVL,DERIVR,DIST,
     2     DT,DTHDA,DTHDAA,DTHDB,
     3     DTHDBB,DTHDE,DTHDEA,DTHDEB,DTHETA,DTHOLD,E,EEE,
     4     EIGLO,EIGLT,EIGRT,EIGUP,EL,EMAX,EMIN,EOLD,EPS,
     5     ER1,ER2,ESTERR,FLO,FMAX,FUP,GUESS,ONE,PIN,
     6     PSIL,PSIPL,PSIPR,PSIR,PX,QX,RAY,WX,
     7     SL,SQL,SQR,SR,T,T1,T2,T3,TAU,THRESH,TMID,TMP,
     8     U,UL,UR,V,WL,X,X50,XAA,XBB,XMID,XSAV,ZAV,
     9     TS,ELIMA,ELIMB,ELIMUP,SUM,SUM0,UT,EU,WU,FOLD,BALLPK,
     A     BESTEIG,BESTEST,OLDEST,EPSM,THA,THB,XT,PUP,PVP,HU,HV,
     B     DTHZ,REMZ,OLDRAY,SAVRAY,RLX,EIGPI,FNEW,AAL,BBL,EPSL,
     C     CHNGLIM,DTHOLDY,AAS,BBS,DTHDAAX,DTHDBBX,DUM,
     D     RATL1,RATL2,RATL3,RATR1,RATR2,RATR3,SL1,SL2,SL3,
     E     TAUM,ADTHETA,FLOUP,PT2,PT3,SAVAA,SAVBB,SAVERR,
     F     BESTAA,BESTBB,BESTMID,BESTEPS,SR1,SR2,SR3
C     ..
C     .. Local Arrays ..
      DOUBLE PRECISION DS(99),PS(99),QS(99),WS(99),DELT(99),PSS(99),
     1     XS(99),YL(3),YR(3),ERL(3),ERR(3),YZL(3),YZR(3)
C     ..
C     .. External Functions ..
      DOUBLE PRECISION EPSLON,P,Q,W,TFROMI
      EXTERNAL EPSLON,P,Q,W,TFROMI
C     ..
C     .. External Subroutines ..
      EXTERNAL AABB,ALFBET,DXDT,EIGFCN,ESTPAC,INTEG,SETMID,THUM,UV
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS,ATAN,COS,EXP,INT,LOG,MAX,MIN,SIGN,SIN,SQRT,TAN
C     ..
C     .. Common blocks ..
      COMMON /DATADT/ASAV,BSAV,C1,C2,INTSAV
      COMMON /DATAF/EIGSAV,IND
      COMMON /RNDOFF/EPSMIN
      COMMON /ZEE/Z
      COMMON /TEEZ/TEE
      COMMON /ZEEZ/JAY,ZEE
      COMMON /PIE/PI,TWOPI,HPI
      COMMON /TDATA/AA,TMID,BB,DTHDAA,DTHDBB,ADDD,MDTHZ
      COMMON /TSAVE/TSAVEL,TSAVER
C     ..
C     Set constants EPSMIN, the computer unit roundoff error, and PI.
C     (Variable ONE set to 1.0 eases precision conversion.)
C
      ONE = 1.0
      EPSMIN = EPSLON(ONE)
      PI = 4.0*ATAN(ONE)
      TWOPI = 2.0*PI
      HPI = 0.5*PI
      do 3 I=1,99
        PSS(I) = 0.0D0
  3   continue
C
C     Set output device number.
C
      IOUT = 6
C
C     Check input parameters for errors.  If errors, return IFLAG = 0.
C
      LOGIC = 1.LE.INTAB .AND. INTAB.LE.4 .AND.
     1        P0ATA*QFATA*P0ATB*QFATB.NE.0.0
      IF (INTAB.EQ.1) LOGIC = LOGIC .AND. A.LT.B
      IF (.NOT.LOGIC) THEN
         IFLAG = 0
         GO TO 150
         END IF
C
C     Set PRIN = .true. to trigger trace printout of successive steps.
C
      PRIN = .FALSE.
      IF (TOL.LT.0.0) PRIN = .TRUE.
C
C     Set EPS to the (initial) integration accuracy.
C
      EPS = 0.0001
C
C     Set logical variables.
C
      AOK = INTAB.LT.3 .AND. P0ATA.LT.0.0 .AND. QFATA.GT.0.0
      BOK = (INTAB.EQ.1 .OR. INTAB.EQ.3) .AND.
     1       P0ATB.LT.0.0 .AND. QFATB.GT.0.0
      SINGA = SINGATA.GT.0.0
      SINGB = SINGATB.GT.0.0
      LCIRCA = CIRCLA.GT.0.0
      LCIRCB = CIRCLB.GT.0.0
      OSCA = OSCILA.GT.0.0
      OSCB = OSCILB.GT.0.0
      TRUNKA = (SINGA .AND. .NOT.LCIRCA) .OR. OSCA
      TRUNKB = (SINGB .AND. .NOT.LCIRCB) .OR. OSCB
      EIGPI = NUMEIG*PI
      NEIG = NUMEIG - 1
      WRITE(21,*) ' NUMEIG = ',NUMEIG
C
C     Initial C1 and C2, used in the mapping between X and T intervals.
C
      C1 = 1.0
      C2 = 0.0
C     DO (SAVE-INPUT-DATA)
         ASAV = A
         BSAV = B
         INTSAV = INTAB
         TAU = ABS(TOL)
         TAUM = MAX(TAU,EPSMIN)
C        END (SAVE-INPUT-DATA)
C
C     Initialize the arrays JAY and ZEE if either end is oscillatory.
C
      IF (OSCA .OR. OSCB) THEN
         DO 5 K=1,100
            JAY(K) = 0
            ZEE(K) = 1.0
    5       CONTINUE
         END IF
C
C     Evaluate P, Q, W to obtain preliminary information about the
C     differential equation.
C
C     DO (SAMPLE-COEFFICIENTS)
         THRESH = 1.0E+17
   10    CONTINUE
            CALL DXDT(EPSMIN,TMP,X50)
            XS(50) = X50
            TS = EPSMIN
            PX = P(X50)
            QX = Q(X50)
            WX = W(X50)
            PS(50) = PX
            QS(50) = QX/PX
            WS(50) = WX/PX
C
C     EMIN = min(Q/W), achieved at X for index value IMIN.
C     EMAX = max(Q/W), achieved at X for index value IMAX.
C     MF and ML are the least and greatest index values, respectively.
C
            XSAV = X50
            EMIN = 0.0
            IF (QX.NE.0.0) EMIN = QX/WX
            EMAX = EMIN
            IMIN = 50
            IMAX = 50
            DO 20 I=49,1,-1
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
               DELT(I) = 0.5*(TS-T)
               XSAV = X
               TS = T
C
C     Try to avoid overflow by stopping when functions are large near A
C     or when w is small near A.
C
               FYNYT = (ABS(WX)+ABS(QX)+1.0/ABS(PX)).LE.THRESH
     1                  .AND. WX.GT.EPSMIN
               IF (QX.NE.0.0 .AND. QX/WX.LT.EMIN) THEN
                  EMIN = QX/WX
                  IMIN = I
                  END IF
               IF (QX.NE.0.0 .AND. QX/WX.GT.EMAX) THEN
                  EMAX = QX/WX
                  IMAX = I
                  END IF
               MF = I
               IF (.NOT.FYNYT) GO TO 30
   20          CONTINUE
   30       CONTINUE
            AAA=T
            IF (.NOT.SINGA) AAA = -1.0
            XSAV = X50
            DO 40 I=51,99
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
               DELT(I-1) = 0.5*(T-TS)
               XSAV = X
               TS = T
C
C     Try to avoid overflow by stopping when functions are large near B
C     or when w is small near A.
C
               FYNYT1 = (ABS(QX)+ABS(WX)+1.0/ABS(PX)).LE.THRESH
     1                   .AND. WX.GT.EPSMIN
               IF (QX.NE.0.0 .AND. QX/WX.LT.EMIN) THEN
                  EMIN = QX/WX
                  IMIN = I
                  END IF
               IF (QX.NE.0.0 .AND. QX/WX.GT.EMAX) THEN
                  EMAX = QX/WX
                  IMAX = I
                  END IF
               ML = I - 1
               IF (.NOT.FYNYT1) GO TO 50
   40          CONTINUE
   50       CONTINUE
            BBB = T
            IF (.NOT.SINGB) BBB = 1.0
            LOGIC = C1.EQ.1.0 .AND. (.NOT.FYNYT .OR. .NOT.FYNYT1)
C
C     Modify (T,X) transformation corresponding to truncated interval.
C
            IF (LOGIC) THEN
               C1 = 0.5*(BBB-AAA)
               C2 = 0.5*(AAA+BBB)
               GO TO 10
               END IF
         IF (OSCA .OR. OSCB) CALL THUM(MF,ML,XS)
C
C     Here we try to determine 'sigma0'.  Initially, we will be
C     satisfied to determine eliml and elimr, the limiting values
C     of q/w, if they exist.
C
         LIMA = .FALSE.
         IF (SINGA .AND. .NOT.LCIRCA) THEN
            RATL1 = QS(MF)/WS(MF)
            RATL2 = QS(MF+1)/WS(MF+1)
            RATL3 = QS(MF+2)/WS(MF+2)
            SL1 = RATL1/(XS(MF+1)-XS(MF))
            SL2 = RATL2/(XS(MF+2)-XS(MF+1))
            SL3 = RATL3/(XS(MF+3)-XS(MF+2))
            IF (ABS(SL2).GE.ABS(SL1) .AND. ABS(SL2).LE.ABS(SL3)) THEN
               ELIMA = RATL1
               LIMA = PS(MF).EQ.PS(MF+1)
               WRITE(*,*) ' There is a limit at a, LIMIT = ',ELIMA
               WRITE(21,*) ' THERE IS A LIMIT AT a, LIMIT = ',ELIMA
               END IF
            END IF
         LIMB = .FALSE.
         IF (SINGB .AND. .NOT.LCIRCB) THEN
            RATR1 = QS(ML)/WS(ML)
            RATR2 = QS(ML-1)/WS(ML-1)
            RATR3 = QS(ML-2)/WS(ML-2)
            SR1 = RATR1/(XS(ML)-XS(ML-1))
            SR2 = RATR2/(XS(ML-1)-XS(ML-2))
            SR3 = RATR3/(XS(ML-2)-XS(ML-3))
            IF (ABS(SR2).GE.ABS(SR1) .AND. ABS(SR2).LE.ABS(SR3)) THEN
               ELIMB = RATR1
               LIMB = PS(ML).EQ.PS(ML-1)
               WRITE(*,*) ' There is a limit at b, LIMIT = ',ELIMB
               WRITE(21,*) ' THERE IS A LIMIT AT b, LIMIT = ',ELIMB
               END IF
            END IF
         LIMUP = .FALSE.
         ELIMUP = EMAX
         IF (LIMA .OR. LIMB) THEN
            LIMUP = .TRUE.
            IF (.NOT.LIMB) THEN
               ELIMUP = ELIMA
            ELSE IF (.NOT.LIMA) THEN
               ELIMUP = ELIMB
            ELSE
               ELIMUP = MIN(ELIMA,ELIMB)
               END IF
            WRITE(21,*) ' THE CONTINUOUS SPECTRUM HAS A LOWER '
            WRITE(21,*) '   BOUND, SIGMA0 = ',ELIMUP
            END IF
C        END (SAMPLE-COEFFICIENTS)
      PIN = EIGPI + PI
      IF (EIG.EQ.0.0) THEN
C        DO (ESTIMATE-EIG)
            SUM0 = 0.0
            IF (OSCA .OR. OSCB) THEN
               EEE = 0.0
C              DO (ESTIMATE-PHASE-ANGLE-CHANGE)
                  CALL ESTPAC(.FALSE.,MF,ML,EEE,SUM0,QS,WS,DS,DELT,
     1                        PS,PSS,TAU,IA,IB,JJL,JJR,SUM,U,UT,ZAV)
C                 END (ESTIMATE-PHASE-ANGLE-CHANGE)
               SUM0 = SUM
               END IF
            EEE = MIN(ELIMUP,EMAX)
C           DO (ESTIMATE-PHASE-ANGLE-CHANGE)
               CALL ESTPAC(.FALSE.,MF,ML,EEE,SUM0,QS,WS,DS,DELT,
     1                     PS,PSS,TAU,IA,IB,JJL,JJR,SUM,U,UT,ZAV)
C              END (ESTIMATE-PHASE-ANGLE-CHANGE)
   55       CONTINUE
               IF (.NOT.LIMUP .AND.
     1             ABS(SUM).GE.10.0*MAX(1.0,ABS(PIN))) THEN
                  IF (SUM.GE.10.0*PIN) THEN
                     IF (EEE.GE.1.0) THEN
                        EEE = EEE/10.0
                     ELSE IF (EEE.LT.-1.0) THEN
                        EEE = 10.0*EEE
                     ELSE
                        EEE = EEE - 1.0
                        END IF
                  ELSE
                     IF (EEE.LE.-1.0) THEN
                        EEE = EEE/10.0
                     ELSE IF (EEE.GT.1.0) THEN
                        EEE = 10.0*EEE
                     ELSE
                        EEE = EEE + 1.0
                        END IF
                     END IF
C                 DO (ESTIMATE-PHASE-ANGLE-CHANGE)
                     CALL ESTPAC(.FALSE.,MF,ML,EEE,SUM0,QS,WS,DS,DELT,
     1                           PS,PSS,TAU,IA,IB,JJL,JJR,SUM,U,UT,ZAV)
C                    END (ESTIMATE-PHASE-ANGLE-CHANGE)
                  GO TO 55
                  END IF
            EU = EEE
            WU = SUM
            IF (SUM.GE.PIN) THEN
               EL = EU
               WL = WU
   60          CONTINUE
                  IF (WL.GE.PIN) THEN
                     EU = EL
                     WU = WL
                     EEE = EL - ((WL-PIN+3.0)/U)**2 - 1.0
C                    DO (ESTIMATE-PHASE-ANGLE-CHANGE)
                        CALL ESTPAC(.FALSE.,MF,ML,EEE,SUM0,QS,WS,DS,
     1                    DELT,PS,PSS,TAU,IA,IB,JJL,JJR,SUM,U,UT,ZAV)
C                       END (ESTIMATE-PHASE-ANGLE-CHANGE)
                     EL = EEE
                     WL = SUM
                     GO TO 60
                     END IF
            ELSE
               EL = EEE
               WL = SUM
               END IF
            IF (LIMUP .AND. WU.LT.PIN) THEN
               EEE = ELIMUP
            ELSE
               IF (U.EQ.0.0) THEN
                  EEE = EMAX + 1.0
C                 DO (ESTIMATE-PHASE-ANGLE-CHANGE)
                     CALL ESTPAC(.FALSE.,MF,ML,EEE,SUM0,QS,WS,DS,DELT,
     1                           PS,PSS,TAU,IA,IB,JJL,JJR,SUM,U,UT,ZAV)
C                    END (ESTIMATE-PHASE-ANGLE-CHANGE)
                  EU = EEE
                  WU = SUM
                  END IF
   70          CONTINUE
                  IF (WU.LE.PIN) THEN
                     EL = EU
                     WL = WU
                     EEE = EU + ((PIN-WU+3.0)/U)**2 + 1.0
C                    DO (ESTIMATE-PHASE-ANGLE-CHANGE)
                        CALL ESTPAC(.FALSE.,MF,ML,EEE,SUM0,QS,WS,DS,
     1                    DELT,PS,PSS,TAU,IA,IB,JJL,JJR,SUM,U,UT,ZAV)
C                       END (ESTIMATE-PHASE-ANGLE-CHANGE)
                     EU = EEE
                     WU = SUM
                     GO TO 70
                     END IF
   80          CONTINUE
                  IF (ABS(IMAX-IMIN).GE.2 .AND. EU.LE.EMAX) THEN
                     IE = (IMAX+IMIN)/2
                     EEE =  QS(IE)/WS(IE)
C                    DO (ESTIMATE-PHASE-ANGLE-CHANGE)
                        CALL ESTPAC(.FALSE.,MF,ML,EEE,SUM0,QS,WS,DS,
     1                    DELT,PS,PSS,TAU,IA,IB,JJL,JJR,SUM,U,UT,ZAV)
C                       END (ESTIMATE-PHASE-ANGLE-CHANGE)
                     IF (SUM.GT.PIN) THEN
                        IMAX = IE
                        WU = SUM
                        EU = EEE
                     ELSE
                        IMIN = IE
                        WL = SUM
                        EL = EEE
                        END IF
                     GO TO 80
                     END IF
C
C     Improve approximation for EIG using bisection or secant method.
C     Substitute 'ballpark' estimate if approximation grows too large.
C
               DEDW = (EU-EL)/(WU-WL)
               FOLD = 0.0
               IF (INTAB.EQ.1) BALLPK = (PIN/(A-B))**2
               IF (INTAB.EQ.1) WRITE(21,*) ' BALLPK = ',BALLPK
               LOGIC = .TRUE.
   90          CONTINUE
                  IF (LOGIC) THEN
                     LOGIC = (WL.LT.PIN-1.0 .OR. WU.GT.PIN+1.0)
                     EEE = EL + DEDW*(PIN-WL)
                     FNEW = MIN(PIN-WL,WU-PIN)
                     IF (FNEW.GT.0.4*FOLD .OR. FNEW.LE.1.0)
     1                  EEE = 0.5*(EL+EU)
                     IF (INTAB.EQ.1 .AND. ABS(EEE).GT.1.0E3*BALLPK) THEN
                        EEE = BALLPK
                        GO TO 100
                     ELSE IF (INTAB.NE.1 .AND. ABS(EEE).GT.1.0E6) THEN
                        EEE = 1.0
                        GO TO 100
                     ELSE
                        FOLD = FNEW
C                       DO (ESTIMATE-PHASE-ANGLE-CHANGE)
                           CALL ESTPAC(.FALSE.,MF,ML,EEE,SUM0,QS,WS,DS,
     1                       DELT,PS,PSS,TAU,IA,IB,JJL,JJR,SUM,U,UT,ZAV)
C                          END (ESTIMATE-PHASE-ANGLE-CHANGE)
                        IF (SUM.LT.PIN) THEN
                           EL = EEE
                           WL = SUM
                        ELSE
                           EU = EEE
                           WU = SUM
                           END IF
                        DEDW = (EU-EL)/(WU-WL)
                        GO TO 90
                        END IF
                     END IF
               END IF
C           END (ESTIMATE-EIG)
         END IF
  100 CONTINUE
      GUESS = EIG
      IF (LIMUP .AND. EEE.GE.ELIMUP) EEE = ELIMUP - 0.01
C     DO (SET-INITIAL-INTERVAL-AND-MATCHPOINT)
         IF (GUESS.NE.0.0) THEN
            EEE=EIG
C           DO (ESTIMATE-PHASE-ANGLE-CHANGE)
               CALL ESTPAC(.FALSE.,MF,ML,EEE,SUM0,QS,WS,DS,DELT,
     1                     PS,PSS,TAU,IA,IB,JJL,JJR,SUM,U,UT,ZAV)
C              END (ESTIMATE-PHASE-ANGLE-CHANGE)
            END IF
C
C     Choose initial interval as large as possible that avoids overflow.
C     JJL and JJR are boundary indices for nonnegativity of EIG*W-Q.
C
         AA = -1.0
         IF (SINGA) AA = TFROMI(JJL)
         BB = 1.0
         IF (SINGB) BB = TFROMI(JJR)
         AA = MIN(-0.01,AA)
         BB = MAX(0.01,BB)
         AA = MIN(AA,-0.95)
         BB = MAX(BB,0.95)
         IF (OSCA) AA = -0.9999
         IF (OSCB) BB = 0.9999
         AAF = AAA
         BBF = BBB
C
C     Determine boundary values ALFA and BETA for theta at A and B.
C
         Z = 1.0
         CALL ALFBET(A,INTAB,AA,A1,A2,EEE,P0ATA,QFATA,SINGA,LCIRCA,
     1               ALFA,KFLAG,DERIVL)
         CALL ALFBET(B,INTAB,BB,B1,B2,EEE,P0ATB,QFATB,SINGB,LCIRCB,
     1               BETA,JFLAG,DERIVR)
         IF (SINGB) BETA = PI - BETA
C
C     Take boundary conditions into account in estimation of EIG.
C
         PIN = EIGPI + BETA - ALFA
         IF (OSCA) PIN = PIN + ALFA
         IF (OSCB) PIN = PIN + PI - BETA
         IF (GUESS.EQ.0.0) THEN
            EEE = EL + DEDW*(PIN-WL)
            IF (.NOT.(OSCA .OR. OSCB) .AND. ABS(EEE).GT.1000.0)
     1         EEE = SIGN(1000.0,EEE)
            IF (INTAB.EQ.1 .AND. ABS(EEE).GT.1.0E3*BALLPK) EEE = BALLPK
            END IF
C        DO (ESTIMATE-PHASE-ANGLE-CHANGE)
            CALL ESTPAC(OSCA.OR.OSCB,MF,ML,EEE,SUM0,QS,WS,DS,DELT,
     1                  PS,PSS,TAU,IA,IB,JJL,JJR,SUM,U,UT,ZAV)
C           END (ESTIMATE-PHASE-ANGLE-CHANGE)
C
C     Choose the constant Z .
C
         IF (U.GT.0.0) Z = ZAV/UT
C
C     Reset boundary values ALFA and BETA .
C
         CALL ALFBET(A,INTAB,AA,A1,A2,EEE,P0ATA,QFATA,SINGA,LCIRCA,
     1               ALFA,KFLAG,DERIVL)
         CALL ALFBET(B,INTAB,BB,B1,B2,EEE,P0ATB,QFATB,SINGB,LCIRCB,
     1               BETA,JFLAG,DERIVR)
         IF (SINGB) BETA = PI - BETA
         IF (PRIN) WRITE(IOUT,'(A,E22.14,A,E22.14)')
     1             ' alfa=',ALFA,'   beta=',BETA
         WRITE(21,'(A,E22.14,A,E22.14)')
     1             ' ALFA=',ALFA,'   BETA=',BETA
C
C     Choose initial matching point TMID .
C
         IMID = 50
         TMID = 0.5*(AA+BB)
         IF (PRIN) WRITE(IOUT,'(A,E15.7,A,F11.8,A,E15.7)')
     1             ' estim=',EEE,'  tmid=',TMID,'  z=',Z
         IF (PRIN) WRITE(IOUT,'(A,F11.8,A,F11.8,A,F11.8,A,F11.8)')
     1             ' aaa=',AAA,'  aa=',AA,'  bb=',BB,'  bbb=',BBB
         WRITE(21,'(A,E15.7,A,F11.8,A,E15.7)')
     1             ' estim=',EEE,'  tmid=',TMID,'  z=',Z
         WRITE(21,'(A,F11.8,A,F11.8,A,F11.8,A,F11.8)')
     1             ' aaa=',AAA,'  aa=',AA,'  bb=',BB,'  bbb=',BBB
C        END (SET-INITIAL-INTERVAL-AND-MATCHPOINT)
      IF (EIG.EQ.0.0 .AND. LIMUP .AND. EEE.GE.ELIMUP)
     1   EEE = ELIMUP - 0.01
C     DO (RESET-TMID)
         CALL SETMID(MF,ML,EEE,QS,WS,IMID,TMID)
C        END (RESET-TMID)
      IF (OSCA .OR. OSCB) THEN
         Z = 1.0
C        DO (PREP-ZEEZ)
            DO 85 I=1,100
               TEE(I) = 1.0
               IF (JAY(I).NE.0) TEE(I) = TFROMI(JAY(I))
   85          CONTINUE
C           END (PREP-ZEEZ)
         END IF
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
C     End preliminary work, begin main task of computing EIG.
C
C     Logical variables have the following meanings if true.
C        AOK    - endpoint A is not singular.
C        BOK    - endpoint B is not singular.
C        BRACKT - EIG has been bracketed.
C        CONVRG - convergence test for EIG has been successfully passed.
C        NEWTON - Newton iteration may be employed.
C        THELT0 - lower bound for EIG has been found.
C        THEGT0 - upper bound for EIG has been found.
C        LIMIT  - upper bound exists with boundary conditions satisfied.
C        ONEDIG - most significant digit can be expected to be correct.
C
      EIG = EEE
      NRAY = 1
      OLDRAY = 1.0E+9
      EXIT = .FALSE.
      FIRSTT = .TRUE.
      LOOP2 = 0
      LOOP3 = 0
      BESTEIG = EIG
      BESTEST = 1.0E+9
      OLDEST = BESTEST
      ADTHETA = 1.0E+9
      NEWTONF = .FALSE.
      CHNGEPS = .FALSE.
      EPSM = EPSMIN
      ENDA = .FALSE.
      ENDB = .FALSE.
      TSAVEL = -1.0
      TSAVER = 1.0
      BRS = .FALSE.
      BRSS = .FALSE.
      SAVERR = 1.0E+9
      AAL = AA
      BBL = BB
      EPSL = EPS
  110 CONTINUE
C     DO (INITIAL-IZE)
         BRACKT = .FALSE.
         CONVRG = .FALSE.
         THELT0 = .FALSE.
         THEGT0 = .FALSE.
         EIGLO = EMIN - 1.0
         FLO = -5.0
         FUP = 5.0
         EIGLT = 0.0
         EIGRT = 0.0
         EIGUP = EMAX + 1.0
         IF (LIMUP) EIGUP = MIN(EMAX,ELIMUP)
         DTHOLD = 1.0
C        END (INITIAL-IZE)
      WRITE(21,*)
      WRITE(21,*) '---------------------------------------------'
      WRITE(21,*) ' INITIAL GUESS FOR EIG = ',EIG
      WRITE(21,*) ' aa,bb = ',AA,BB
C     DO UNTIL(CONVRG .OR. EXIT)
      DO 120 NITER = 1,40
         WRITE(*,*)
         WRITE(*,*) ' ******************** '
C        DO (SET-TMID-AND-BOUNDARY-CONDITIONS)
            WRITE(*,*) ' set tmid and boundary conditions '
            V = EIG*WS(IMID) - QS(IMID)
C           IF (V.LE.0.0) DO (RESET-TMID)
               IF (V.LE.0.0) CALL SETMID(MF,ML,EIG,QS,WS,IMID,TMID)
C              END (RESET-TMID)
C           DO (RESET-BOUNDARY-CONDITIONS)
               DERIVL = 0.0
               IF (SINGA) CALL ALFBET(A,INTAB,AA,A1,A2,EIG,
     1             P0ATA,QFATA,.TRUE.,LCIRCA,ALFA,KFLAG,DERIVL)
               DERIVR = 0.0
               IF (SINGB) THEN
                  CALL ALFBET(B,INTAB,BB,B1,B2,EIG,P0ATB,QFATB,
     1                        .TRUE.,LCIRCB,BETA,JFLAG,DERIVR)
                  BETA = PI - BETA
                  END IF
               IF (PRIN) WRITE(IOUT,'(A,E22.14,A,E22.14)')
     1                   ' alfa=',ALFA,'   beta=',BETA
               WRITE(21,'(A,E22.14,A,E22.14)')
     1                   ' ALFA=',ALFA,'   BETA=',BETA
C              END (RESET-BOUNDARY-CONDITIONS)
C
C     Check that boundary conditions can be satisfied at singular
C     endpoints.  If not, try for slightly altered EIG consistent
C     with boundary conditions.
C
            IF (LIMUP .AND. EIG.NE.GUESS .AND. .NOT.BRACKT) THEN
               KFLAG = 1
               IF (SINGA .AND. .NOT.LCIRCA) CALL ALFBET(A,INTAB,AA,
     1            A1,A2,EIG,P0ATA,QFATA,.TRUE.,.FALSE.,TMP,KFLAG,TMP)
               JFLAG = 1
               IF (SINGB .AND. .NOT.LCIRCB) CALL ALFBET(B,INTAB,BB,
     1            B1,B2,EIG,P0ATB,QFATB,.TRUE.,.FALSE.,TMP,JFLAG,TMP)
               IF ((KFLAG.NE.1 .OR. JFLAG.NE.1) .AND.
     1             (THELT0 .AND. EIGLO.LT.ELIMUP)) THEN
                  EIGUP = MIN(ELIMUP+2.0*EPSMIN,EIGUP)
                  IF (EIG.NE.EIGLO .AND. EIG.NE.EIGUP) THEN
                     EIG = 0.05*EIGLO + 0.95*EIGUP
                  ELSE
                     IFLAG = 11
                     WRITE(21,*) ' IFLAG = 11 '
                     EXIT = .TRUE.
                     GO TO 130
                     END IF
                  END IF
               END IF
C           END (SET-TMID-AND-BOUNDARY-CONDITIONS)
C        DO (OBTAIN-DTHETA-WITH-ONE-CORRECT-DIGIT)
            IF (PRIN) WRITE(IOUT,'(/A,E22.14,A,E10.3,A,E10.3)')
     1                ' guess=',EIG,'  eps=',EPS,'  tmid=',TMID
C           DO (INTEGRATE-FOR-DTHETA)
C              DO (SET-INITIAL-CONDITIONS)
                  THA = ALFA
                  DTHDEA = DERIVL
                  DTHDAA = 0.0
                  IF (SINGA .AND. .NOT.LCIRCA) THEN
                     CALL DXDT(AA,DT,X)
                     PX = P(X)/Z
                     QX = Q(X)/Z
                     WX = W(X)/Z
                     C = EIG*WX - QX
                     DTHDAA = -(COS(ALFA)**2/PX + C*SIN(ALFA)**2)*DT
C
C     Two special cases for DTHDAA .
C
                     IF (C.GE.0.0 .AND. P0ATA.LT.0.0 .AND. QFATA.LT.0.0)
     1                   DTHDAA = DTHDAA + ALFA*DT/(X-A)
                     IF (C.GE.0.0 .AND. P0ATA.GT.0.0 .AND. QFATA.GT.0.0)
     1                   DTHDAA = DTHDAA + (ALFA-0.5*PI)*DT/(X-A)
                     END IF
                  THB = BETA
                  DTHDEB = -DERIVR
                  DTHDBB = 0.0
                  IF (SINGB .AND. .NOT.LCIRCB) THEN
                     CALL DXDT(BB,DT,X)
                     PX = P(X)/Z
                     QX = Q(X)/Z
                     WX = W(X)/Z
                     C = EIG*WX - QX
                     DTHDBB = -(COS(BETA)**2/PX + C*SIN(BETA)**2)*DT
C
C     Two special cases for DTHDBB .
C
                     IF (C.GE.0.0 .AND. P0ATB.LT.0.0 .AND. QFATB.LT.0.0)
     1                   DTHDBB = DTHDBB + (PI-BETA)*DT/(B-X)
                     IF (C.GE.0.0 .AND. P0ATB.GT.0.0 .AND. QFATB.GT.0.0)
     1                   DTHDBB = DTHDBB + (0.5*PI-BETA)*DT/(B-X)
                     END IF
C                 END (SET-INITIAL-CONDITIONS)
               EIGSAV = EIG
C
C     YL = (theta,d(theta)/d(eig),d(theta)/da)
C
               YL(1) = ALFA
               YL(2) = DTHDEA
               YL(3) = 0.0
C
               CALL INTEG(AA,THA,DTHDAA,DTHDEA,TMID,A1,A2,EPS,YL,ERL,
     1                    LCIRCA,AOK,SINGA,OSCA,IFLAG)
               IF (IFLAG.EQ.5) THEN
                  IFLAG = 51
                  WRITE(21,*) ' IFLAG = 51 '
                  EXIT = .TRUE.
                  GO TO 130
                  END IF
               DTHDA = DTHDAA*EXP(-2.0*YL(3))
C
C     YR = (theta,d(theta)/d(eig),d(theta)/db)
C
               YR(1) = BETA + EIGPI - PI
               YR(2) = DTHDEB
               YR(3) = 0.0
C
               CALL INTEG(BB,THB,DTHDBB,DTHDEB,TMID,B1,B2,EPS,YR,ERR,
     1                    LCIRCB,BOK,SINGB,OSCB,IFLAG)
               IF (IFLAG.EQ.5) THEN
                  IFLAG = 52
                  WRITE(21,*) ' IFLAG = 52 '
                  EXIT = .TRUE.
                  GO TO 130
                  END IF
               DTHDB = DTHDBB*EXP(-2.0*YR(3))
C
               ER1 = ERL(1) - ERR(1)
               ER2 = ERL(2) - ERR(2)
C
               IF (OSCA .OR. OSCB) THEN
                  Z = 1.0
                  CALL DXDT(TMID,TMP,XT)
                  CALL UV(XT,U,PUP,V,PVP,HU,HV)
                  EIGSAV = 0.0
                  CALL INTEG(AA,THA,DTHDAA,DTHDEA,TMID,A1,A2,EPS,
     1                       YZL,ERL,LCIRCA,AOK,SINGA,OSCA,IFLAG)
                  IF (IFLAG.EQ.5) THEN
                     IFLAG = 53
                     WRITE(21,*) ' IFLAG = 53 '
                     EXIT = .TRUE.
                     GO TO 130
                     END IF
                  CALL INTEG(BB,THB,DTHDBB,DTHDEB,TMID,B1,B2,EPS,
     1                       YZR,ERR,LCIRCB,BOK,SINGB,OSCB,IFLAG)
                  IF (IFLAG.EQ.5) THEN
                     IFLAG = 54
                     WRITE(21,*) ' IFLAG = 54 '
                     EXIT = .TRUE.
                     GO TO 130
                     END IF
                  EIGSAV = EIG
                  DTHZ = YZR(1) - YZL(1)
                  MDTHZ = DTHZ/PI
                  REMZ = DTHZ - MDTHZ*PI
                  IF (DTHZ.LT.0.0 .AND. REMZ.LT.0.0) THEN
                     MDTHZ = MDTHZ - 1
                     REMZ = REMZ + PI
                     END IF
                  IF (REMZ.GT.3.14) MDTHZ = MDTHZ + 1
                  END IF
C
C     Record the environment parameters of this most recent
C     successful integration.
C
               AAS = AA
               BBS = BB
C
C     DTHETA measures theta difference from left and right integrations.
C
C              DO (FORM-DTHETA)
                  DTHETA = YL(1) - YR(1) - EIGPI
                  IF (OSCA .OR. OSCB) DTHETA = DTHETA + MDTHZ*PI
                  DTHDE = YL(2) - YR(2)
C                 END (FORM-DTHETA)
               ADTHETA = ABS(DTHETA)
               ONEDIG = (ABS(ER1).LE.0.5*ABS(DTHETA) .AND.
     1                  ABS(ER2).LE.0.5*ABS(DTHDE)) .OR.
     2                  MAX(ADTHETA,ABS(ER1)).LT.1.0E-6
               FIRSTT = .FALSE.
C              END (INTEGRATE-FOR-DTHETA)
            CHNGEPS = .FALSE.
            CONVRG = .FALSE.
            OLDNEWT = NEWTON
            NEWTON = ABS(DTHETA).LT.0.06 .AND. BRACKT
            IF (NEWTON)
     1          ONEDIG = ONEDIG .OR. ABS(DTHETA+ER1).LT.0.5*DTHOLD
            IF (PRIN) WRITE(IOUT,'(A,E15.7,A,E15.7)')
     1                ' dtheta=',DTHETA,'   dthde=',DTHDE
            IF (PRIN) WRITE(IOUT,'(/A,E15.7,A,E15.7)')
     1                ' thetal=',YL(1),'   thetar=',YR(1)
C           END (OBTAIN-DTHETA-WITH-ONE-CORRECT-DIGIT)
         IF (.NOT.ONEDIG .AND. BRS) THEN
            EXIT = .TRUE.
            WRITE(21,*) ' NOT ONEDIG '
            GO TO 130
            END IF
C        DO (SET-BRACKET-DATA)
            WRITE(*,*) ' set-bracket '
            IF (DTHETA.GT.0.0) THEN
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
C     EIG is bracketed when both THEGT0=.true. and THELT0=.true.
C
            BRACKT = THELT0 .AND. THEGT0
            IF (PRIN) WRITE(IOUT,'(A,E22.14,A,E22.14)')
     1                ' eigrt=',EIGRT,'  eigup=',EIGUP
            IF (PRIN) WRITE(IOUT,'(A,E22.14,A,E22.14)')
     1                ' eiglt=',EIGLT,'  eiglo=',EIGLO
C           END (SET-BRACKET-DATA)
         IF (BRACKT) LOOP2 = 0
C        DO (TEST-FOR-CONVERGENCE)
C
C     Measure convergence after adding separate contributions to error.
C
            FLOUP = MIN(ABS(FLO),ABS(FUP))
            T1 = (ABS(DTHETA)+MAX(ABS(ERL(1)),ABS(ERR(1))))/ABS(DTHDE)
            T2 = (1.0+AA)*ABS(DTHDA)/ABS(DTHDE)
            T3 = (1.0-BB)*ABS(DTHDB)/ABS(DTHDE)
            WRITE(21,*) ' FLO,FUP,FLOUP = ',FLO,FUP,FLOUP
            WRITE(21,*) ' DTHDE,DTHDA,DTHDB = ',DTHDE,DTHDA,DTHDB
            PT2 = (AAF-AA)*DTHDA/DTHDE
            PT3 = (BBF-BB)*DTHDB/DTHDE
            IF (.NOT.TRUNKA) THEN
               T2 = 0.0
               PT2 = 0.0
               END IF
            IF (.NOT.TRUNKB) THEN
               T3 = 0.0
               PT3 = 0.0
               END IF
            ESTERR = T1 + T2 + T3
            ESTERR = ESTERR/MAX(ONE,ABS(EIG))
            CONVRG = ESTERR.LE.TAUM .AND. NEWTON
            WRITE(21,*) ' T1,T2,T3 = ',T1,T2,T3
            WRITE(21,*) ' PT2,PT3 = ',PT2,PT3
            WRITE(21,*) ' TMID,EPS = ',TMID,EPS
            WRITE(21,*) ' ONEDIG,BRACKT,NEWTON,CONVRG = ',
     1                    ONEDIG,BRACKT,NEWTON,CONVRG
            WRITE(21,*) ' EIG,DTHETA,ESTERR = ',EIG,DTHETA,ESTERR
            IF (BRACKT .AND. (ESTERR.LT.BESTEST .OR. .NOT.BRS) .AND.
     1          ADTHETA.LT.0.1) THEN
               BESTAA = AA
               BESTBB = BB
               BESTMID = TMID
               BESTEPS = EPS
               BESTEIG = EIG
               BESTEST = ESTERR
               BRS = BRACKT
               IF (BRS) BRSS = BRS
               WRITE(21,*) ' BESTEIG,BESTEST = ',BESTEIG,BESTEST
               WRITE(21,*) ' BRS = ',BRS
               END IF
            IF (THEGT0) WRITE(21,*) '         EIGUP = ',EIGUP
            IF (THELT0) WRITE(21,*) '         EIGLO = ',EIGLO
            WRITE(21,*) '---------------------------------------------'
            IF (PRIN) WRITE(IOUT,'(A,L2)') ' converge=',CONVRG
            IF (PRIN .AND. .NOT.CONVRG) WRITE(IOUT,'(A,E15.7)')
     1                                  ' estim. acc.=',ESTERR
C           END (TEST-FOR-CONVERGENCE)
         IF (CONVRG) THEN
            WRITE(*,*) ' number of iterations was ',NITER
            WRITE(21,*) ' NUMBER OF ITERATIONS WAS ',NITER
            WRITE(*,*) '-----------------------------------------------'
            GO TO 130
         ELSE
            IF (NEWTON) THEN
               IF (OLDNEWT .AND. ADTHETA.GT.0.8*ABS(DTHOLD)) THEN
                  WRITE(21,*) ' ADTHETA,DTHOLD = ',ADTHETA,DTHOLD
                  WRITE(21,*) ' NEWTON DID NOT IMPROVE EIG '
                  NEWTONF = .TRUE.
                  LOOP3 = LOOP3 + 1
               ELSE IF (TRUNKA .OR. TRUNKB) THEN
                  ENDA = ADTHETA.LT.1.0 .AND. ABS(PT2).GT.MAX(TAUM,T1)
     1                   .AND. TRUNKA
                  ENDB = ADTHETA.LT.1.0 .AND. ABS(PT3).GT.MAX(TAUM,T1)
     1                   .AND. TRUNKB
                  IF (ENDA .OR. ENDB) THEN
                     NEWTON = .FALSE.
                  ELSE IF ((T2+T3).GT.T1 .AND. ADTHETA.LT.1.0 .AND.
     1               (AA.LE.AAF .AND. BB.GE.BBF)) THEN
                     WRITE(*,*) ' RESIDUAL TRUNCATION ERROR DOMINATES '
                     EXIT = .TRUE.
                     IFLAG = 9
                     WRITE(21,*) ' IFLAG = 9 '
                     GO TO 130
                     END IF
                  END IF
               IF (NEWTONF .OR. ENDA .OR. ENDB) THEN
                  WRITE(21,*) ' NEWTONF,ENDA,ENDB = ',NEWTONF,ENDA,ENDB
                  EXIT = .TRUE.
                  GO TO 130
                  END IF
C              DO (NEWTON'S-METHOD)
                  WRITE(*,*) ' Newton''s method '
                  RLX = 1.2
                  IF (BRACKT) RLX = 1.0
                  EIG = EIG - RLX*DTHETA/DTHDE
                  IF (EIG.LE.EIGLO .OR. EIG.GE.EIGUP)
     1               EIG = 0.5*(EIGLO+EIGUP)
                  WRITE(21,*) ' NEWTON: EIG = ',EIG
C                 END (NEWTON'S-METHOD)
            ELSE IF (BRACKT) THEN
               WRITE(*,*) ' bracket '
C              DO (SECANT-METHOD)
                  WRITE(*,*) ' do secant method '
                  FMAX = MAX(-FLO,FUP)
                  EOLD = EIG
                  EIG = 0.5*(EIGLO+EIGUP)
                  IF (FMAX.LE.1.5) THEN
                     U = -FLO/(FUP-FLO)
                     DIST = EIGUP - EIGLO
                     EIG = EIGLO + U*DIST
                     V = MIN(EIGLT,EIGRT)
                     IF (EIG.LE.V) EIG = 0.5*(EIG+V)
                     V = MAX(EIGLT,EIGRT)
                     IF (EIG.GE.V) EIG = 0.5*(EIG+V)
                     DE = EIG - EOLD
                     IF (ABS(DE).LT.EPSMIN) THEN
                        TOL = ABS(DE)/MAX(ONE,ABS(EIG))
                        IFLAG = 6
                        EXIT = .TRUE.
                        GO TO 130
                        END IF
                     END IF
                     WRITE(21,*) ' SECANT: EIG = ',EIG
C                 END (SECANT-METHOD)
            ELSE
C              DO (TRY-FOR-BRACKET)
                  LOOP2 = LOOP2 + 1
                  IF (LOOP2.GT.9 .AND. .NOT.LIMUP) THEN
                     IFLAG = 12
                     WRITE(21,*) ' IFLAG = 12 '
                     EXIT = .TRUE.
                     GO TO 130
                     END IF
                  IF (EIG.EQ.EEE) THEN
                     IF (GUESS.NE.0.0) DEDW = 1.0/DTHDE
                     CHNG = -0.6*(DEDW+1.0/DTHDE)*DTHETA
                     IF (EIG.NE.0.0 .AND. ABS(CHNG).GT.0.1*ABS(EIG))
     1                  CHNG = -0.1*SIGN(EIG,DTHETA)
                  ELSE
                     CHNG = -1.2*DTHETA/DTHDE
                     IF (CHNG.EQ.0.) CHNG = 0.1*MAX(1.0,ABS(EIG))
                     WRITE(21,*) ' IN BRACKET, 1,CHNG = ',CHNG
C
C     Limit change in EIG to a factor of 10.
C
                     IF (ABS(CHNG).GT.(1.0+10.0*ABS(EIG))) THEN
                        CHNG = SIGN(1.0+10.0*ABS(EIG),CHNG)
                     WRITE(21,*) ' IN BRACKET, 2,CHNG = ',CHNG
                     ELSE IF (ABS(EIG).GE.1.0 .AND.
     1                  ABS(CHNG).LT.0.1*ABS(EIG)) THEN
                        CHNG = 0.1*SIGN(EIG,CHNG)
                     WRITE(21,*) ' IN BRACKET, 3,CHNG = ',CHNG
                        END IF
                     IF (DTHOLD.LT.0.0 .AND. LIMUP .AND.
     1                  CHNG.GT.(ELIMUP-EIG)) THEN
                        CHNG = 0.95*(ELIMUP-EIG)
                        WRITE(21,*) ' ELIMUP,EIG,CHNG = ',
     1                                ELIMUP,EIG,CHNG
                        IF (CHNG.LT.EPSMIN) THEN
                           WRITE(*,*) ' ELIMUP,EIG = ',ELIMUP,EIG
                           WRITE(21,*) ' IN BRACKET, CHNG.LT.EPSMIN '
                           ENDA = TRUNKA .AND. AA.GT.AAF .AND.
     1                            (0.5*DTHETA+(AAF-AA)*DTHDA).GT.0.0
                           ENDB = TRUNKB .AND. BB.LT.BBF .AND.
     1                            (0.5*DTHETA-(BBF-BB)*DTHDB).GT.0.0
                           IF (.NOT.(ENDA .OR. ENDB)) THEN
                              NUMEIG = NEIG - INT(-DTHETA/PI)
                              WRITE(*,*) ' new numeig = ',NUMEIG
                              WRITE(21,*) ' NEW NUMEIG = ',NUMEIG
                              END IF
                              IFLAG = 3
                              GO TO 150
                           END IF
                        END IF
                     END IF
                  EOLD = EIG
                  CHNGLIM = 2.0*ESTERR*MAX(ONE,ABS(EIG))
                  WRITE(21,*) ' CHNGLIM = ',CHNGLIM
                  IF (ADTHETA.LT.0.06 .AND. ABS(CHNG).GT.CHNGLIM .AND.
     1                CHNGLIM.NE.0.0) CHNG = SIGN(CHNGLIM,CHNG)
                  IF ((THELT0 .AND. CHNG.LT.0.0) .OR.
     1                (THEGT0 .AND. CHNG.GT.0.0)) CHNG = -CHNG
                  EIG = EIG + CHNG
                  WRITE(21,*) ' BRACKET: EIG = ',EIG
C                 END (TRY-FOR-BRACKET)
               END IF
            END IF
         IF (IFLAG.EQ.3) GO TO 130
         IF (NITER.GE.3 .AND. DTHOLDY.EQ.DTHETA) THEN
            IFLAG = 7
            WRITE(21,*) ' IFLAG = 7 '
            EXIT = .TRUE.
            GO TO 130
            END IF
         DTHOLDY = DTHOLD
         DTHOLD = DTHETA
         WRITE(*,*) ' number of iterations was ',NITER
         WRITE(*,*) '-----------------------------------------------'
  120    CONTINUE
      IFLAG = 8
      WRITE(21,*) ' IFLAG = 8 '
      EXIT = .TRUE.
  130 CONTINUE
      IF (AA.EQ.AAL .AND. BB.EQ.BBL .AND. EPS.LT.EPSL .AND.
     1    ESTERR.GE.0.5*SAVERR) GO TO 140
      EPSL = EPS
      AAL = AA
      BBL = BB
      TOL = BESTEST
      EIG = BESTEIG
      IF (EXIT) THEN
         WRITE(21,*) ' EXIT '
         IF (FIRSTT) THEN
            IF (IFLAG.EQ.51 .OR. IFLAG.EQ.53) THEN
               IF (AA.LT.-0.1) THEN
                  WRITE(*,*) ' FIRST COMPLETE INTEGRATION FAILED. '
                  WRITE(21,*) ' FIRST COMPLETE INTEGRATION FAILED. '
                  IF (AA.EQ.-1.0) GO TO 150
                  AAF = AA
                  CALL AABB(AA,-ONE)
                  WRITE(21,*) ' aa MOVED FROM ',AAf,' IN TO ',AA
                  EXIT = .FALSE.
                  GO TO 110
               ELSE
                  WRITE(21,*) ' aa.GE.-0.1 '
                  IFLAG = 13
                  GO TO 150
                  END IF
            ELSE IF (IFLAG.EQ.52 .OR. IFLAG.EQ.54) THEN
               IF (BB.GT.0.1) THEN
                  WRITE(*,*) ' FIRST COMPLETE INTEGRATION FAILED. '
                  WRITE(21,*) ' FIRST COMPLETE INTEGRATION FAILED. '
                  IF (BB.EQ.1.0) GO TO 150
                  BBF = BB
                  CALL AABB(BB,-ONE)
                  WRITE(21,*) ' bb MOVED FROM ',BBf,' IN TO ',BB
                  EXIT = .FALSE.
                  GO TO 110
               ELSE
                  WRITE(21,*) ' bb.LE.0.1 '
                  IFLAG = 14
                  GO TO 150
                  END IF
               END IF
         ELSE IF (IFLAG.EQ.51 .OR. IFLAG.EQ.53) THEN
            WRITE(*,*) ' A COMPLETE INTEGRATION FAILED. '
            WRITE(21,*) ' A COMPLETE INTEGRATION FAILED. '
            IF (CHNGEPS) THEN
               EPS = 5.0*EPS
               EPSM = EPS
               WRITE(21,*) ' EPS INCREASED TO ',EPS
            ELSE
               AAF = AA
               CALL AABB(AA,-ONE)
               WRITE(21,*) ' aa MOVED FROM ',AAf,' IN TO ',AA
               END IF
            EXIT = .FALSE.
            GO TO 110
         ELSE IF (IFLAG.EQ.52 .OR. IFLAG.EQ.54) THEN
            WRITE(*,*) ' A COMPLETE INTEGRATION FAILED. '
            WRITE(21,*) ' A COMPLETE INTEGRATION FAILED. '
            IF (CHNGEPS) THEN
               EPS = 5.0*EPS
               EPSM = EPS
               WRITE(21,*) ' EPS INCREASED TO ',EPS
            ELSE
               BBF = BB
               CALL AABB(BB,-ONE)
               WRITE(21,*) ' bb MOVED FROM ',BBf,' IN TO ',BB
               END IF
            EXIT = .FALSE.
            GO TO 110
         ELSE IF (IFLAG.EQ.6) THEN
            WRITE(*,*) ' IN SECANT, CHNG.LT.EPSMIN '
            WRITE(21,*) ' IN SECANT, CHNG.LT.EPSMIN '
            GO TO 140
         ELSE IF (IFLAG.EQ.7) THEN
            WRITE(*,*) ' DTHETA IS REPEATING '
            WRITE(21,*) ' DTHETA IS REPEATING '
            GO TO 140
         ELSE IF (IFLAG.EQ.8) THEN
            WRITE(*,*) ' NUMBER OF ITERATIONS REACHED SET LIMIT '
            WRITE(21,*) ' NUMBER OF ITERATIONS REACHED SET LIMIT '
            GO TO 140
         ELSE IF (IFLAG.EQ.9) THEN
            WRITE(21,*) ' RESIDUAL TRUNCATION ERROR DOMINATES '
            GO TO 140
         ELSE IF (IFLAG.EQ.11) THEN
            WRITE(*,*) ' IN TRY FOR BRACKET, CHNG.LT.EPSMIN '
            WRITE(21,*) ' IN TRY FOR BRACKET, CHNG.LT.EPSMIN '
            GO TO 140
         ELSE IF (IFLAG.EQ.12) THEN
            WRITE(*,*) ' FAILED TO GET A BRACKET. '
            WRITE(21,*) ' FAILED TO GET A BRACKET. '
            GO TO 140
         ELSE IF (NEWTONF .OR. .NOT.ONEDIG) THEN
            IF (LOOP3.GE.3) THEN
               WRITE(21,*) ' NEWTON IS NOT GETTING ANYWHERE '
               NEWTONF = .FALSE.
               GO TO 140
               END IF
            WRITE(21,*) ' BESTEST,OLDEST = ',BESTEST,OLDEST
            IF (EPS.GT.EPSM .AND. BESTEST.LT.OLDEST) THEN
               CHNGEPS = .TRUE.
               SAVERR = ESTERR
               EPS = 0.2*EPS
               WRITE(*,*) ' EPS REDUCED TO ',EPS
               WRITE(21,*) ' EPS REDUCED TO ',EPS
               EXIT = .FALSE.
               NEWTON = .FALSE.
               OLDEST = BESTEST
               GO TO 110
            ELSE
               IF (EPS.LE.EPSM) THEN
                  WRITE(*,*) ' EPS CANNOT BE REDUCED FURTHER. '
                  WRITE(21,*) ' EPS CANNOT BE REDUCED FURTHER. '
                  IFLAG = 2
                  GO TO 140
               ELSE
                  WRITE(*,*) ' no more improvement '
                  WRITE(21,*) ' NO MORE IMPROVEMENT '
                  GO TO 140
                  END IF
               END IF
         ELSE IF (ENDA) THEN
            CALL AABB(AA,ONE)
            AA = MAX(AA,AAA)
            IF (AA.LE.AAF) THEN
               WRITE(21,*) ' NO MORE IMPROVEMENT '
               GO TO 140
               END IF
            WRITE(*,*) ' aa MOVED OUT TO ',AA
            WRITE(21,*) ' aa MOVED OUT TO ',AA
            EXIT = .FALSE.
            GO TO 110
         ELSE IF (ENDB) THEN
            CALL AABB(BB,ONE)
            BB = MIN(BB,BBB)
            IF (BB.GE.BBF) THEN
               WRITE(21,*) ' NO MORE IMPROVEMENT '
               GO TO 140
               END IF
            WRITE(*,*) ' bb MOVED OUT TO ',BB
            WRITE(21,*) ' bb MOVED OUT TO ',BB
            EXIT = .FALSE.
            GO TO 110
            END IF
         END IF
  140 CONTINUE
C
C     If CONVRG is false, check that any truncation error might possibly
C     be reduced or that the integrations might be done more accurately.
C
      IF (.NOT.CONVRG .AND. IFLAG.LT.50 .AND. IFLAG.NE.11) THEN
         SAVAA = AA
         SAVBB = BB
         IF (EPS.GT.EPSM .AND. ESTERR.LT.0.5*SAVERR) THEN
            WRITE(21,*) ' SAVERR,ESTERR = ',SAVERR,ESTERR
            SAVERR = ESTERR
            EPS = 0.2*EPS
            WRITE(*,*) ' EPS REDUCED TO ',EPS
            WRITE(21,*) ' 2,EPS REDUCED TO ',EPS
            EXIT = .FALSE.
            NEWTON = .FALSE.
            OLDEST = BESTEST
            GO TO 110
         ELSE IF (ABS(PT2).GT.TAUM .OR. ABS(PT3).GT.TAUM) THEN
            IF ((AAS-AAF).GT.2.0*EPSMIN .AND. ABS(PT2).GT.TAUM) THEN
               CALL AABB(AA,ONE)
               AA = MAX(AA,AAA)
               IF (AA.GT.AAF .AND. AA.LT.SAVAA) THEN
                  WRITE(*,*) ' aa MOVED OUT TO ',AA
                  WRITE(21,*) ' 3,aa MOVED OUT TO ',AA
                  EXIT = .FALSE.
                  END IF
               END IF
            IF ((BBF-BBS).GT.2.0*EPSMIN .AND. ABS(PT3).GT.TAUM) THEN
               CALL AABB(BB,ONE)
               BB = MIN(BB,BBB)
               IF (BB.GT.SAVBB .AND. BB.LT.BBF) THEN
                  WRITE(*,*) ' bb MOVED OUT TO ',BB
                  WRITE(21,*) ' 3,bb MOVED OUT TO ',BB
                  EXIT = .FALSE.
                  END IF
               END IF
            IF (.NOT.EXIT .AND. (AA.NE.SAVAA .OR. BB.NE.SAVBB))
     1         GO TO 110
            END IF
         END IF
      IF (PRIN) WRITE(IOUT,'(A,I7,A,E22.14,A,E10.3)')
     1          ' numeig=',NUMEIG,'  eig=',EIG,'  tol=',TOL
      WRITE(21,*) 'NUMEIG = ',NUMEIG,' EIG = ',EIG,' TOL = ',TOL
      IF (BRSS .AND. BESTEST.LT.0.05) THEN
C        DO (COMPUTE-EIGENFUNCTION-DATA)
C
C     Convert from T to X values, fill 7 of first 9 locations of SLFUN.
C
            CALL DXDT(TMID,TMP,XMID)
            CALL DXDT(AA,TMP,XAA)
            CALL DXDT(BB,TMP,XBB)
            SLFUN(1) = XMID
            SLFUN(2) = XAA
            SLFUN(3) = ALFA
            SLFUN(5) = XBB
            SLFUN(6) = BETA + EIGPI
            SLFUN(8) = EPS
            SLFUN(9) = Z
C
C     Compute SLFUN(4), SLFUN(7) towards normalizing the eigenfunction.
C
            EIGSAV = EIG
            THA = ALFA
            DTHDAAX = 0.0
            YL(1) = 0.0
            YL(2) = 0.0
            YL(3) = 0.0
            CALL INTEG(AA,THA,DTHDAAX,DTHDEA,TMID,A1,A2,EPS,YL,ERL,
     1                 LCIRCA,AOK,SINGA,OSCA,JFLAG)
            THB = BETA
            DTHDBBX = 0.0
            CALL INTEG(BB,THB,DTHDBBX,DTHDEB,TMID,B1,B2,EPS,YR,ERR,
     1                 LCIRCB,BOK,SINGB,OSCB,JFLAG)
            YR(1) = YR(1) + EIGPI
            SL = SIN(YL(1))
            SR = SIN(YR(1))
            CL = COS(YL(1))
            CR = COS(YR(1))
            UL = (YL(2)-DTHDEA*EXP(-2.0*YL(3)))*Z
            UR = (YR(2)-DTHDEB*EXP(-2.0*YR(3)))*Z
            DUM = 0.5*LOG(UL-UR)
            SLFUN(4) = -YL(3) - DUM
            SLFUN(7) = -YR(3) - DUM
C           END (COMPUTE-EIGENFUNCTION-DATA)
C        DO (CHECK-MATCHING-VALUES-OF-EIGENFUNCTION)
C
C     Perform final check on EIG. Return IFLAG = 2
C     if not accurate enough.
C
            DEN = UL*SR*SR - UR*SL*SL
            E = ABS(SR)/SQRT(DEN)
            PSIL = E*SL
            PSIPL = E*CL*Z
            SQL = E*E*UL
            E = ABS(SL)/SQRT(DEN)
            PSIR = E*SR
            PSIPR = E*CR*Z
            SQR = E*E*UR
            ADDD = PSIL*PSIR.LT.0.0 .AND. PSIPL*PSIPR.LT.0.0
            RAY = EIG + (PSIL*PSIPL-PSIR*PSIPR)/(SQL-SQR)
            IF (PRIN) THEN
               WRITE(IOUT,'(A,E22.14)') ' ray=',RAY
               WRITE(IOUT,'(A,E22.14,A,E22.14)')
     1         ' psil=',PSIL,'  psir=',PSIR
               WRITE(IOUT,'(A,E22.14,A,E22.14)')
     1         ' psipl=',PSIPL,'  psipr=',PSIPR
               WRITE(IOUT,'(A,E22.14,A,E22.14)')
     1         ' sql=',SQL,'  sqr=',SQR
               END IF
C           END (CHECK-MATCHING-VALUES-OF-EIGENFUNCTION)
C
C     If next condition is .true., then something is apparently wrong
C     with the accuracy of EIG. Bisect and go through the loop again.
C
         WRITE(21,*) ' EIG,RAY = ',EIG,RAY
         IF (ABS(RAY-EIG).GT.2.0*TAUM*MAX(ONE,ABS(EIG))) THEN
            NRAY = NRAY + 1
            WRITE(*,*) ' nray = ',nray
            WRITE(21,*) ' NRAY,RAY,OLDRAY = ',NRAY,RAY,OLDRAY
            IF (ESTERR.GE.0.5*SAVERR) THEN
               IFLAG = 2
               GO TO 150
               END IF
            EIG = 0.5*(EIG+RAY)
            SAVRAY = OLDRAY
            OLDRAY = RAY
            IF (OLDRAY.NE.SAVRAY .AND. NRAY.LT.3) GO TO 110
            END IF
C        DO (GENERATE-EIGENFUNCTION-VALUES)
            CALL EIGFCN(EIGPI,A1,A2,B1,B2,AOK,SINGA,LCIRCA,OSCA,
     1                  BOK,SINGB,LCIRCB,OSCB,SLFUN,ISLFUN)
C           END (GENERATE-EIGENFUNCTION-VALUES)
         IFLAG = 1
         END IF
C
C     IF THE ESTIMATED ACCURACY IMPLIES THAT THE COMPUTED VALUE
C     OF THE EIGENVALUE IS UNCERTAIN, SIGNAL BY IFLAG = 2.
C
      IF ((ABS(EIG).LE.ONE .AND. TOL.GE.ABS(EIG)) .OR.
     1    (ABS(EIG).GT.ONE .AND. TOL.GE.ONE)) IFLAG = 2
C
  150 CONTINUE
      WRITE(21,*) ' BEST aa,bb,TMID,EPS = ',
     1              BESTAA,BESTBB,BESTMID,BESTEPS
      WRITE(21,*) ' BRS,BESTEIG,BESTEST = ',BRS,BESTEIG,BESTEST
      WRITE(21,*) ' IFLAG = ',IFLAG
      WRITE(21,*) '********************************************'
      WRITE(21,*)
      RETURN
      END
      SUBROUTINE AABB(TEND,OUT)
      DOUBLE PRECISION TEND,OUT
C     **********
C
C     This subroutine moves aa or bb further out or closer in.
C     TEND is either aa or bb.  If OUT is positive, TEND is moved
C     further out, and if OUT is negative, TEND is moved closer in.
C
C     **********
C     .. Local Scalars ..
      INTEGER J
      DOUBLE PRECISION DIF,REM,TTEND
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS,LOG10,SIGN
C     ..
      REM = 1.0 - ABS(TEND)
      J = LOG10(5.0/REM)
      IF (OUT.GT.0.0) J = J + 1
      DIF = 10.0**(-J)
      TTEND = SIGN(1.0-DIF,TEND)
      IF (TTEND.EQ.TEND) THEN
         IF (OUT.LT.0.0) TTEND = SIGN(1.0-10.0*DIF,TEND)
         IF (OUT.GT.0.0) TTEND = SIGN(1.0-0.1*DIF,TEND)
         END IF
      TEND = TTEND
      RETURN
      END
      SUBROUTINE ALFBET(XEND,INTAB,TT,COEF1,COEF2,EIG,P0,QF,SING,
     1                  LCIRC,VALUE,IFLAG,DERIV)
      INTEGER INTAB,IFLAG
      LOGICAL SING,LCIRC
      DOUBLE PRECISION XEND,TT,COEF1,COEF2,EIG,P0,QF,VALUE,DERIV
C     **********
C
C     This subroutine computes a boundary value for a specified endpoint
C     of the interval for a Sturm-Liouville problem in the form
C
C       -(p(x)*y'(x))' + q(x)*y(x) = eig*w(x)*y(x)  on (a,b)
C
C     for user-supplied coefficient functions P, Q, and W.  It is called
C     from SLEIGN2.  Both regular and singular endpoints are treated.
C
C     Subprograms called
C
C       user-supplied ..... p,q,w
C
C       sleign2-supplied .. dxdt,extrap
C
C     **********
C     .. Scalars in Common ..
      DOUBLE PRECISION Z
C     ..
C     .. Local Scalars ..
      LOGICAL LOGIC
      DOUBLE PRECISION C,CD,D,HH,ONE,PI,PUP,PVP,PX,QX,WX,T,TEMP,TTS,
     1     U,V,X,XDENOM,XNUM
C     ..
C     .. External Functions ..
      DOUBLE PRECISION P,Q,W
      EXTERNAL P,Q,W
C     ..
C     .. External Subroutines ..
      EXTERNAL DXDT,EXTRAP,UV
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS,ATAN,SIGN,SQRT
C     ..
C     .. Common blocks ..
      COMMON /ZEE/Z
C     ..
C     Set machine dependent constant.
C
C     PI (variable ONE set to 1.0 eases precision conversion).
      ONE = 1.0
      PI = 4.0*ATAN(ONE)
C
      IFLAG = 1
      DERIV = 0.0
      IF (.NOT.SING) THEN
         VALUE = 0.5*PI
         IF (COEF1.NE.0.0) VALUE = ATAN(-Z*COEF2/COEF1)
         LOGIC = (TT.LT.0.0 .AND. VALUE.LT.0.0) .OR.
     1           (TT.GT.0.0 .AND. VALUE.LE.0.0)
         IF (LOGIC) VALUE = VALUE + PI
      ELSE IF (LCIRC) THEN
         CALL DXDT(TT,TEMP,X)
         CALL UV(X,U,PUP,V,PVP,TEMP,TEMP)
         XNUM = COEF1*U+COEF2*V
         XDENOM = (COEF1*PUP+COEF2*PVP)/Z
         VALUE = ATAN2(XNUM,XDENOM)
         IF (XNUM.LT.0.0) VALUE = VALUE + 2.0*PI
      ELSE
         LOGIC = (INTAB.EQ.2 .AND. TT.GT.0.0) .OR.
     1           (INTAB.EQ.3 .AND. TT.LT.0.0) .OR.
     2           INTAB.EQ.4 .OR. (P0.GT.0.0 .AND. QF.LT.0.0)
         IF (LOGIC) THEN
            T = SIGN(ONE,TT)
            TTS = TT
            CALL EXTRAP(T,TTS,EIG,VALUE,DERIV,IFLAG)
         ELSE
            CALL DXDT(TT,TEMP,X)
            PX = P(X)/Z
            QX = Q(X)/Z
            WX = W(X)/Z
            C = 2.0*(EIG*WX-QX)
            IF (C.LT.0.0) THEN
               VALUE = 0.0
               IF (P0.GT.0.0) VALUE = 0.5*PI
            ELSE
               HH = ABS(XEND-X)
               D = 2.0*HH/PX
               CD = C*D*HH
               IF (P0.GT.0.0) THEN
                  VALUE = C*HH
                  IF (CD.LT.1.0) VALUE = VALUE/(1.0+SQRT(1.0-CD))
                  VALUE = VALUE + 0.5*PI
               ELSE
                  VALUE = D
                  IF (CD.LT.1.0) VALUE = VALUE/(1.0+SQRT(1.0-CD))
                  END IF
               END IF
            END IF
         END IF
      RETURN
      END
      SUBROUTINE DXDT(T,DT,X)
      DOUBLE PRECISION T,DT,X
C     **********
C
C     This subroutine transforms coordinates from T on (-1,1) to
C     X on (A,B) in the solution of a Sturm-Liouville problem.
C     It is called from subroutines SLEIGN2, ALFBET, F, and EXTRAP.
C
C     **********
C     .. Scalars in Common ..
      INTEGER INTAB
      DOUBLE PRECISION A,B,C1,C2
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION U
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS
C     ..
C     .. Common blocks ..
      COMMON /DATADT/A,B,C1,C2,INTAB
C     ..
      U = C1*T + C2
      GO TO (10,20,30,40), INTAB
   10 CONTINUE
         DT = C1*0.5*(B-A)
         X = 0.5*((B+A)+(B-A)*U)
         RETURN
   20 CONTINUE
         DT = C1*2.0/(1.0-U)**2
         X = A + (1.0+U)/(1.0-U)
         RETURN
   30 CONTINUE
         DT = C1*2.0/(1.0+U)**2
         X = B - (1.0-U)/(1.0+U)
         RETURN
   40 CONTINUE
         DT = C1/(1.0-ABS(U))**2
         X = U/(1.0-ABS(U))
         RETURN
      END
      SUBROUTINE EIGFCN(EIGPI,A1,A2,B1,B2,AOK,SINGA,LCIRCA,OSCA,
     1                  BOK,SINGB,LCIRCB,OSCB,SLFUN,ISLFUN)
      INTEGER ISLFUN
      LOGICAL AOK,SINGA,LCIRCA,OSCA,BOK,SINGB,LCIRCB,OSCB
      DOUBLE PRECISION EIGPI,A1,A2,B1,B2
      DOUBLE PRECISION SLFUN(ISLFUN+9)
C     **********
C     **********
C     .. Scalars in Common ..
      INTEGER MDTHZ
      LOGICAL ADDD
      DOUBLE PRECISION AA,BB,DTHDAA,DTHDBB,TMID
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
         IF (SLFUN(9+I).LE.TMID) NMID = I
   10    CONTINUE
      IF (NMID.GT.0) THEN
         T = AA
         YL(1) = SLFUN(3)
         YL(2) = 0.0
         YL(3) = 0.0
         LCIRC = LCIRCA
         OK = AOK
         SING = SINGA
         EFF = 0.0
         DO 20 J=1,NMID
            TM = SLFUN(J+9)
            IF (TM.LT.AA .OR. TM.GT.BB) THEN
               WRITE(*,*) ' t.lt.aa .or. t.gt.bb '
               STOP
               END IF
            THT = YL(1)
            DTHDAT = DTHDAA*EXP(-2.0*EFF)
            DTHDET = YL(2)
            IF (TM.GT.AA) THEN
               CALL INTEG(T,THT,DTHDAT,DTHDET,TM,A1,A2,SLFUN(8),
     1                    YL,ERL,LCIRC,OK,SING,OSCA,IFLAG)
               IF (OSCA) THEN
                  EFF = YL(3)
               ELSE
                  LCIRC = .FALSE.
                  SING = .FALSE.
                  EFF = EFF + YL(3)
                  END IF
               END IF
            SLFUN(J+9) = SIN(YL(1))*EXP(EFF+SLFUN(4))
            T = TM
            IF (T.GT.-1.0) OK = .TRUE.
            IF (T.LT.-0.9 .AND .OSCA) THEN
               OK = .FALSE.
               T = AA
               YL(1) = SLFUN(3)
               YL(2) = 0.0
               YL(3) = 0.0
               END IF
   20       CONTINUE
         END IF
      IF (NMID.LT.ISLFUN) THEN
         T = BB
         YR(1) = SLFUN(6) - EIGPI
         YR(2) = 0.0
         YR(3) = 0.0
         LCIRC = LCIRCB
         OK = BOK
         SING = SINGB
         EFF = 0.0
         DO 30 J=ISLFUN,NMID+1,-1
            TM = SLFUN(J+9)
            IF (TM.LT.AA .OR. TM.GT.BB) THEN
               WRITE(*,*) ' t.lt.aa .or. t.gt.bb '
               STOP
               END IF
            THT = YR(1)
            DTHDBT = DTHDBB*EXP(-2.0*EFF)
            DTHDET = YR(2)
            IF (TM.LT.BB) THEN
               CALL INTEG(T,THT,DTHDBT,DTHDET,TM,B1,B2,SLFUN(8),
     1                    YR,ERR,LCIRC,OK,SING,OSCB,IFLAG)
               IF (OSCB) THEN
                  EFF = YR(3)
               ELSE
                  LCIRC = .FALSE.
                  SING = .FALSE.
                  EFF = EFF + YR(3)
                  END IF
               END IF
            SLFUN(J+9) = SIN(YR(1)+EIGPI)*EXP(EFF+SLFUN(7))
            IF (ADDD) SLFUN(J+9) = -SLFUN(J+9)
            T = TM
            IF (T.LT.1.0) OK = .TRUE.
            IF (T.GT.0.9 .AND. OSCB) THEN
               OK = .FALSE.
               T = BB
               YR(1) = SLFUN(6) - EIGPI
               YR(2) = 0.0
               YR(3) = 0.0
               END IF
   30       CONTINUE
         END IF
      RETURN
      END
      SUBROUTINE ESTPAC(IOSC,MF,ML,EEE,SUM0,QS,WS,DS,DELT,PS,PSS,TAU,
     1                  IA,IB,JJL,JJR,SUM,U,UT,ZAV)
      INTEGER MF,ML,IA,IB,JJL,JJR
      LOGICAL IOSC
      DOUBLE PRECISION EEE,SUM0,TAU,SUM,U,UT,ZAV
      DOUBLE PRECISION QS(ML),WS(ML),DS(ML),DELT(ML),PS(ML),PSS(ML)
C     **********
C
C     This subroutine estimates the change in 'phase angle' in the
C     eigenvalue determination of a Sturm-Liouville problem in the form
C
C       -(p(x)*y'(x))' + q(x)*y(x) = eig*w(x)*y(x)  on (a,b)
C
C     for user-supplied coefficient functions P, Q, and W.
C
C     The subroutine approximates the (trapezoidal rule) integral of
C
C        sqrt((eig*w-q)/p)
C
C     where the integral is taken over those X in (A,B) for which
C
C        (eig*w-q)/p .gt. 0
C
C     **********
C     .. Local Scalars ..
      INTEGER J,JJ,JSAV,MF1
      DOUBLE PRECISION PSUM,RT,RTSAV,V,WW,WSAV,DPSUM,DPSUMT,ZAVJ,ZAVSAV
C     ..
C     .. Arrays in Common ..
      INTEGER JAY(100)
      DOUBLE PRECISION ZEE(100)
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS,MAX,MIN,SIGN,SQRT
C     ..
C     .. Common blocks ..
      COMMON /ZEEZ/JAY,ZEE
C     ..
      IA = MF
      IB = 80
C
C     SUM accumulates the integral approximation.  U measures the total
C     length of subintervals where (EIG*W-Q)/P .gt. 0.0.  ZAV is the
C     average value of sqrt((EIG*W-Q)*P) over those subintervals.
C
      IF (.NOT.IOSC) THEN
         JJL = 99
         JJR = 1
         SUM = 0.0
         U = 0.0
         UT = 0.0
         ZAV = 0.0
         WSAV = EEE*WS(MF) - QS(MF)
         IF (WSAV.GT.0.0) THEN
            RTSAV = SIGN(SQRT(WSAV),PS(MF))
         ELSE
            RTSAV = 0.0
            END IF
         DO 10 J=MF+1,ML
            WW = EEE*WS(J) - QS(J)
            IF (WW.GT.0.0) THEN
               IF (J.GT.80) IB = J
               U = U + DS(J-1)
               UT = UT + DELT(J-1)
               RT = SIGN(SQRT(WW),PS(J))
            ELSE
               RT = 0.0
               IF (U.EQ.0.0 .AND. RTSAV.EQ.0.0 .AND. IA.LE.19)
     1             IA = IA + 1
               END IF
            IF (WW.EQ.0.0 .OR.WSAV.EQ.0.0 .OR. WW.EQ.SIGN(WW,WSAV)) THEN
               V = RT + RTSAV
            ELSE
               V = (WW*RT+WSAV*RTSAV)/ABS(WW-WSAV)
               END IF
            WSAV = WW
            RTSAV = RT
            PSUM = DS(J-1)*V
            IF (EEE.EQ.0.0) THEN
               PSS(J) = PSUM
            ELSE
               DPSUM = PSUM - PSS(J)
               DPSUMT = DPSUM*DELT(J-1)/DS(J-1)
               IF (DPSUMT.GT.0.001*TAU) THEN
                  JJL = MIN(JJL,J)
                  JJR = MAX(JJR,J)
                  END IF
               END IF
            SUM = SUM + PSUM
            IF (U.GT.0.0) ZAV = ZAV + DELT(J-1)*V*ABS(PS(J)+PS(J-1))
   10       CONTINUE
         SUM = 0.5*SUM - SUM0
         ZAV = 0.25*ZAV
      ELSE
         JJ = 1
         JAY(1) = MF
   20    CONTINUE
            SUM = 0.0
            U = 0.0
            UT  =  0.0
            ZAV = 0.0
            ZAVJ = 0.0
            MF1 = JAY(JJ)
            WSAV = EEE*WS(MF1) - QS(MF1)
            IF (WSAV.GT.0.0) THEN
               RTSAV = SIGN(SQRT(WSAV),PS(MF1))
            ELSE
               RTSAV = 0.0
               END IF
            DO 30 J=MF1+1,ML
               WW = EEE*WS(J) - QS(J)
               IF (WW.GT.0.0) THEN
                  IF (J.GT.80) IB = J
                  U = U + DS(J-1)
                  UT = UT + DELT(J-1)
                  RT = SIGN(SQRT(WW),PS(J))
               ELSE
                  RT = 0.0
                  IF (U.EQ.0.0 .AND. RTSAV.EQ.0.0 .AND. IA.LE.19)
     1               IA = IA + 1
                  END IF
               IF (WW.EQ.0.0 .OR. WSAV.EQ.0.0 .OR.
     1             WW.EQ.SIGN(WW,WSAV)) THEN
                  V = RT + RTSAV
               ELSE
                  V = (WW*RT+WSAV*RTSAV)/ABS(WW-WSAV)
                  END IF
               WSAV = WW
               RTSAV = RT
               PSUM = DS(J-1)*V
               SUM = SUM + PSUM
               IF (U.GT.0.0) ZAV = ZAV + DELT(J-1)*V*ABS(PS(J)+PS(J-1))
               IF (U.NE.0.0) THEN
                  IF (ZAVJ.EQ.0.0) JSAV = J
                  ZAVJ = 0.25*ZAV/UT
                  IF (J.EQ.JSAV) ZAVSAV = ZAVJ
                  IF (2.0*ZAVJ.LT.ZAVSAV .OR. ZAVJ.GT.2.0*ZAVSAV) THEN
                     JJ = JJ + 1
                     JAY(JJ) = J
                     ZEE(JJ) = 0.5*(ZAVJ+ZAVSAV)
                     GO TO 40
                     END IF
                  END IF
   30          CONTINUE
   40       CONTINUE
            IF (J.GT.ML) THEN
               JJ = JJ + 1
               JAY(JJ) = ML
               ZEE(JJ) = 0.5*(ZAVJ+ZAVSAV)
               END IF
            IF (J.LT.ML) GO TO 20
         SUM = 0.5*SUM
         ZAV = 0.25*ZAV
         END IF
      IB = IB + 1
      RETURN
      END
      SUBROUTINE EXTRAP(T,TT,EIG,VALUE,DERIV,IFLAG)
      INTEGER IFLAG
      DOUBLE PRECISION T,TT,EIG,VALUE,DERIV
C     **********
C
C     This subroutine is called from ALFBET in determining boundary
C     values at a singular endpoint of the interval for a
C     Sturm-Liouville problem in the form
C
C       -(p(x)*y'(x))' + q(x)*y(x) = eig*w(x)*y(x)  on (a,b)
C
C     for user-supplied coefficient functions P, Q, and W.
C
C     EXTRAP, which in turn calls INTPOL, extrapolates the function
C
C        arctan(1.0/sqrt(-p*(eig*w-q)))
C
C     from its values for T within (-1,1) to an endpoint.
C
C     Subprograms called
C
C       user-supplied ..... p,q,w
C
C       sleign2-supplied .. dxdt,intpol
C
C     **********
C     .. Scalars in Common ..
      DOUBLE PRECISION Z
C     ..
C     .. Local Scalars ..
      INTEGER KGOOD
      DOUBLE PRECISION ANS,CTN,ERROR,PROD,PX,QX,WX,T1,TEMP,X
C     ..
C     .. Local Arrays ..
      DOUBLE PRECISION FN1(5),XN(5)
C     ..
C     .. External Functions ..
      DOUBLE PRECISION P,Q,W
      EXTERNAL P,Q,W
C     ..
C     .. External Subroutines ..
      EXTERNAL DXDT,INTPOL
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS,ATAN,SQRT,TAN
C     ..
C     .. Common blocks ..
      COMMON /ZEE/Z
C     ..
      IFLAG = 1
      KGOOD = 0
      T1 = TT
   10 CONTINUE
         CALL DXDT(T1,TEMP,X)
         PX = P(X)/Z
         QX = Q(X)/Z
         WX = W(X)/Z
         PROD = -PX*(EIG*WX-QX)
         IF (PROD.LE.0.0) THEN
            T1 = 0.5*(T1+T)
            IF ((1.0+(T1-T)**2).GT.1.0) GO TO 10
            IF (PROD.GT.-1.0E-6) VALUE = 2.0*ATAN(1.0)
            IFLAG = 5
            WRITE(21,*) ' In EXTRAP, iflag = 5'
            RETURN
         ELSE
            KGOOD = KGOOD + 1
            XN(KGOOD) = T1
            FN1(KGOOD) = ATAN(1.0/SQRT(PROD))
            T1 = 0.5*(T+T1)
            IF (KGOOD.LT.5) GO TO 10
            END IF
      T1 = 0.01
      CALL INTPOL(5,XN,FN1,T,T1,3,ANS,ERROR)
      VALUE = ABS(ANS)
      CTN = 1.0/TAN(VALUE)
      DERIV = 0.5*PX*WX/CTN/(1.0+CTN**2)
      TT = XN(1)
      RETURN
      END
      SUBROUTINE F(U,Y,YP)
      DOUBLE PRECISION U
      DOUBLE PRECISION Y(2),YP(3)
C     **********
C
C     This subroutine evaluates the derivative functions for use with
C     integrator GERK in solving a Sturm-Liouville problem in the form
C
C       -(p(x)*y'(x))' + q(x)*y(x) = eig*w(x)*y(x)  on (a,b)
C
C     for user-supplied coefficient functions P, Q, and W.
C
C     Subprograms called
C
C       user-supplied ..... p,q,w
C
C       sleign2-supplied .. dxdt
C
C     **********
C     .. Scalars in Common ..
      INTEGER IND
      DOUBLE PRECISION EIG,Z
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION C,C2,DT,QX,WX,S,S2,T,TH,V,WW,X,XP
C     ..
C     .. External Functions ..
      DOUBLE PRECISION P,Q,W
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
      COMMON /ZEE/Z
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
      YP(1) = DT*(XP*C2+V*S2)
      IF (IND.EQ.1) THEN
         WW = (XP-V)*S*C
         YP(2) = DT*(-2.0*WW*Y(2)+WX*S2)
         YP(3) = DT*WW
      ELSE IF (IND.EQ.2) THEN
         YP(2) = YP(2)/YP(1)
         YP(3) = YP(3)/YP(1)
         YP(1) = 1.0/YP(1)
      ELSE IF (IND.EQ.3) THEN
      ELSE
         YP(1) = 1.0/YP(1)
         END IF
      RETURN
      END
      DOUBLE PRECISION FUNCTION FF(ALFLAM)
      DOUBLE PRECISION ALFLAM
C     **********
C     **********
C     .. Scalars in Common ..
      INTEGER IND,INDD
      DOUBLE PRECISION CC,EIGSAV,HPI,PI,THETU,THETV,TWOPI,
     1     UL,UR,VL,VR,UB,VB,Z
C     ..
C     .. Local Scalars ..
      INTEGER LFLAG
      LOGICAL AOK,BOK,LCIRCA,LCIRCB,OSCA,OSCB,SINGA,SINGB
      DOUBLE PRECISION AA,BB,DTHDAA,DTHDBB,DTHDEA,DTHDEB,DUM,EPS,LAMBDA,
     1     PVPB,PVPL,PVPR,RHOB,RHOL,RHOR,THA,THB,THL,THR,TMID,
     1     EPSMIN
C     ..
C     .. Local Arrays ..
      DOUBLE PRECISION ERR(3),Y(3)
C     ..
C     .. External Subroutines ..
      EXTERNAL INTEG
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC COS,EXP,SIN
C     ..
C     .. Common blocks ..
      COMMON /EPP2/CC,UL,UR,VL,VR,UB,VB,IND
      COMMON /DATAF/EIGSAV,INDD
      COMMON /ZEE/Z
      COMMON /PIE/PI,TWOPI,HPI
      COMMON /THET/THETU,THETV
      COMMON /RNDOFF/EPSMIN
C     ..
C
C     (THIS ROUTINE IS BEING MODIFIED FOR PERIODIC PROBLEMS WHICH
C     ARE NOT NECESSARILY REGULAR, BUT IS NOT YET COMPLETE.)
C
      AOK = .TRUE.
      LCIRCA = .FALSE.
      SINGA = .FALSE.
      OSCA = .FALSE.
      BOK = .TRUE.
      LCIRCB = .FALSE.
      SINGB = .FALSE.
      OSCB = .FALSE.
C
      AA = -1.0
      BB = 1.0
C
C     SET TMID SO THAT IT IS NOT IN THE EXACT MIDDLE.
C
      TMID = 0.1*PI
C
      INDD = 1
      LAMBDA = ALFLAM
      EIGSAV = LAMBDA
      EPS = EPSMIN
C
      IF (IND.GE.2) GO TO 10
C
  50  CONTINUE
C
C     FOR U:
C
         THA = HPI
         Y(1) = THA
         Y(2) = 1.0
         Y(3) = 0.0
         DTHDAA = 0.0
         DTHDEA = 1.0
         LFLAG = 1
         CALL INTEG(AA,THA,DTHDAA,DTHDEA,TMID,DUM,DUM,EPS,Y,ERR,
     1              LCIRCA,AOK,SINGA,OSCA,LFLAG)
         IF (LFLAG.EQ.5) THEN
            EPS = 10.0*EPS
            GO TO 50
            END IF
         RHOL = EXP(Y(3))
         THL = Y(1)
         UL = RHOL*SIN(THL)
C
         THB = HPI
         Y(1) = THB
         Y(2) = 1.0
         Y(3) = 0.0
         DTHDBB = 0.0
         DTHDEB = 1.0
         LFLAG = 1
         CALL INTEG(BB,THB,DTHDBB,DTHDEB,TMID,DUM,DUM,EPS,Y,ERR,
     1              LCIRCB,BOK,SINGB,OSCB,LFLAG)
         IF (LFLAG.EQ.5) THEN
            EPS = 10.0*EPS
            GO TO 50
            END IF
         RHOR = EXP(Y(3))
         THR = Y(1)
         UR = RHOR*SIN(THR)
C
C     FOR V:
C
         THA = 0.0
         Y(1) = THA
         Y(2) = 1.0
         Y(3) = 0.0
         DTHDAA = 0.0
         DTHDEA = 1.0
         LFLAG = 1
         CALL INTEG(AA,THA,DTHDAA,DTHDEA,TMID,DUM,DUM,EPS,Y,ERR,
     1              LCIRCA,AOK,SINGA,OSCA,LFLAG)
         IF (LFLAG.EQ.5) THEN
            EPS = 10.0*EPS
            GO TO 50
            END IF
         RHOL = EXP(Y(3))
         THL = Y(1)
         VL = RHOL*SIN(THL)
         PVPL = Z*RHOL*COS(THL)
C
         THB = 0.0
         Y(1) = THB
         Y(2) = 1.0
         Y(3) = 0.0
         DTHDBB = 0.0
         DTHDEB = 1.0
         LFLAG = 1
         CALL INTEG(BB,THB,DTHDBB,DTHDEB,TMID,DUM,DUM,EPS,Y,ERR,
     1              LCIRCB,BOK,SINGB,OSCB,LFLAG)
         IF (LFLAG.EQ.5) THEN
            EPS = 10.0*EPS
            GO TO 50
            END IF
         RHOR = EXP(Y(3))
         THR = Y(1)
         VR = RHOR*SIN(THR)
         PVPR = Z*RHOR*COS(THR)
      FF = (VR*(UL*PVPL-1.0)-CC*VL*(UR*PVPR-1.0))*(VL-VR/CC) -
     1      VL*VR*(UL-CC*UR)*(PVPL-PVPR/CC)
      RETURN
C
  10  CONTINUE
C
C     FOR U:
C
         THA = HPI
         Y(1) = THA
         Y(2) = 1.0
         Y(3) = 0.0
         DTHDAA = 0.0
         DTHDEA = 1.0
         LFLAG = 1
         CALL INTEG(AA,THA,DTHDAA,DTHDEA,BB,DUM,DUM,EPS,Y,ERR,
     1              LCIRCA,AOK,SINGA,OSCA,LFLAG)
         IF (LFLAG.EQ.5) THEN
            EPS = 10.0*EPS
            GO TO 10
            END IF
         RHOB = EXP(Y(3))
         THB = Y(1)
         THETU = THB
         UB = RHOB*SIN(THB)
C
C     FOR V:
C
         THA = 0.0
         Y(1) = THA
         Y(2) = 1.0
         Y(3) = 0.0
         DTHDAA = 0.0
         DTHDEA = 1.0
         LFLAG = 1
         CALL INTEG(AA,THA,DTHDAA,DTHDEA,BB,DUM,DUM,EPS,Y,ERR,
     1              LCIRCA,AOK,SINGA,OSCA,LFLAG)
         IF (LFLAG.EQ.5) THEN
            EPS = 10.0*EPS
            GO TO 10
            END IF
         RHOB = EXP(Y(3))
         THB = Y(1)
         THETV = THB
         VB = RHOB*SIN(THB)
         PVPB = Z*RHOB*COS(THB)
      FF = CC*PVPB + UB/CC - 2.0
      RETURN
      END
      SUBROUTINE FIT(TH1,TH,TH2)
      DOUBLE PRECISION TH1,TH,TH2
C     **********
C
C     This program converts TH into an 'equivalent' angle between
C     TH1 and TH2.  We assume TH1.LT.TH2 and PI.LE.(TH2-TH1).
C
C     **********
C     .. Scalars in Common ..
      DOUBLE PRECISION PI,TWOPI,HPI
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC AINT
C     ..
C     .. Common blocks ..
      COMMON /PIE/PI,TWOPI,HPI
C     ..
      IF (TH.LT.TH1) TH = TH + AINT((TH1-TH+PI)/PI)*PI
      IF (TH.GT.TH2) TH = TH - AINT((TH-TH2+PI)/PI)*PI
      RETURN
      END
      SUBROUTINE FZ(UU,Y,YP)
      DOUBLE PRECISION UU
      DOUBLE PRECISION Y(2),YP(3)
C     **********
C
C     This subroutine evaluates the derivative of the function PHI
C     in the regularization at a singular endpoint.
C
C     Here, HU means -(PUP)' + QU.
C
C     **********
C     .. Scalars in Common ..
      INTEGER IND
      DOUBLE PRECISION EIG
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION AU,AV,A1122,A12,A21,B1122,B12,B21,C,C2,D,DT,
     1     HU,HV,PHI,PUP,PVP,S,SC,S2,T,U,V,WW,WX,X
C     ..
C     .. External Functions ..
      DOUBLE PRECISION W
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
      YP(1) = -DT*(A1122*SC+A12*S2-A21*C2)/D
      IF (IND.EQ.1) THEN
         WW = 2.0*(A12+A21)*SC + A1122*(C2-S2)
         YP(2) = -DT*(WW*Y(2)+2.0*B1122*SC+B12*S2-B21*C2)/D
         YP(3) = 0.5*DT*WW/D
      ELSE IF (IND.EQ.2) THEN
         YP(2) = YP(2)/YP(1)
         YP(3) = YP(3)/YP(1)
         YP(1) = 1.0/YP(1)
      ELSE IF (IND.EQ.3) THEN
      ELSE
         YP(1) = 1.0/YP(1)
         END IF
      RETURN
      END
      SUBROUTINE GERKZ(F,NEQ,Y,TIN,TOUT,REPS,AEPS,LFLAG,ER,WORK,IWORK)
      INTEGER NEQ,LFLAG
      INTEGER IWORK(5)
      DOUBLE PRECISION TIN,TOUT,REPS,AEPS
      DOUBLE PRECISION Y(3),ER(3),WORK(27)
      EXTERNAL F
C     **********
C     **********
C     .. Scalars in Common ..
      DOUBLE PRECISION EPSMIN,Z
C     ..
C     .. Local Scalars ..
      INTEGER I,J,K,L,LLFLAG
      DOUBLE PRECISION T,TOUTS
C     ..
C     .. Arrays in Common ..
      INTEGER JAY(100)
      DOUBLE PRECISION TEE(100),ZEE(100)
C     ..
C     .. Local Arrays ..
      DOUBLE PRECISION U(3)
C     ..
C     .. External Subroutines ..
      EXTERNAL GERK,THTOTHZ,THZTOTH
C     ..
C     .. Intrinsic Functions
      INTRINSIC MAX,MIN
C     ..
C     .. Common blocks ..
      COMMON /RNDOFF/EPSMIN
      COMMON /ZEE/Z
      COMMON /TEEZ/TEE
      COMMON /ZEEZ/JAY,ZEE
C     ..
      T = TIN
      IF (TIN.LT.TOUT) THEN
         DO 10 I=1,19
            IF (TEE(I)-EPSMIN.LE.TIN .AND. TIN.LT.TEE(I+1)+EPSMIN) J = I
            IF (TEE(I)-EPSMIN.LT.TOUT.AND.TOUT.LE.TEE(I+1)+EPSMIN) L = I
   10       CONTINUE
         DO 30 K=J,L
            TOUTS = MIN(TOUT,TEE(K+1))
            Z = ZEE(K+1)
            IF (Z.EQ.0.0) Z = 1.0
            CALL THTOTHZ(Y,Z,U)
            LLFLAG = 1
   20       CONTINUE
               CALL GERK(F,NEQ,U,T,TOUTS,REPS,AEPS,LLFLAG,ER,WORK,IWORK)
               IF (LLFLAG.GT.3) THEN
                  WRITE(21,*) ' llflag = ',LLFLAG
                  LFLAG = 5
                  RETURN
                  END IF
               IF (LLFLAG.EQ.3 .OR. LLFLAG.EQ.-2) GO TO 20
            CALL THZTOTH(U,Z,Y)
   30       CONTINUE
      ELSE
         DO 40 I=20,2,-1
            IF (TEE(I-1)-EPSMIN.LT.TIN .AND. TIN.LE.TEE(I)+EPSMIN) J = I
            IF (TEE(I-1)-EPSMIN.LE.TOUT.AND.TOUT.LT.TEE(I)+EPSMIN) L = I
   40       CONTINUE
         DO 60 K=J,L,-1
            TOUTS = MAX(TOUT,TEE(K-1))
            Z = ZEE(K)
            IF (Z.EQ.0.0) Z = 1.0
            CALL THTOTHZ(Y,Z,U)
            LLFLAG = 1
   50       CONTINUE
               CALL GERK(F,NEQ,U,T,TOUTS,REPS,AEPS,LLFLAG,ER,WORK,IWORK)
               IF (LLFLAG.GT.3) THEN
                  WRITE(21,*) ' llflag = ',LLFLAG
                  LFLAG = 5
                  RETURN
                  END IF
               IF (LLFLAG.EQ.3 .OR. LLFLAG.EQ.-2) GO TO 50
            CALL THZTOTH(U,Z,Y)
   60       CONTINUE
         END IF
      TIN = T
      LFLAG = LLFLAG
      RETURN
      END
      SUBROUTINE INTEG(TEND,THEND,DTHDAA,DTHDE,TMID,COEF1,COEF2,
     1                 EPS,Y,ER,LCIRC,OK,SING,OSC,IFLAG)
      INTEGER IFLAG
      LOGICAL LCIRC,OK,SING,OSC
      DOUBLE PRECISION TEND,THEND,DTHDAA,DTHDE,TMID,COEF1,COEF2,EPS
      DOUBLE PRECISION Y(3),ER(3)
      CHARACTER*16 PREC
C     **********
C     **********
C     .. Scalars in Common ..
      INTEGER IND
      DOUBLE PRECISION EIG,HPI,PI,TSAVEL,TSAVER,TWOPI,Z
C     ..
C     .. Arrays in Common ..
      INTEGER NT(2)
      DOUBLE PRECISION TT(7,2),YY(7,3,2)
C     ..
C     .. Local Scalars ..
      INTEGER I,J,KFLAG,K2PI,LFLAG,M
      LOGICAL LOGIC
      DOUBLE PRECISION D,DDD,EFF,HU,HV,PHI,PHI0,PUP,PVP,T,TMP,U,V,XT0,
     1     ZSAV,C,DTHIN,DUM,DPHIDE,DPHIDAA,FAC2,P1,PYPZ,PYPZ0,
     2     Q1,RHOSQ,S,TSTAR,TINTHZ,TH0,TH,THBAR,THETA,
     3     THIN,THU,THV,THU0,THV0,TIN,TOUT,W1,XSTAR,XT,YSTAR,YZ,YZ0
C     ..
C     .. Local Arrays ..
      INTEGER IWORK(5)
      DOUBLE PRECISION ERZ(3),WORK(27),YP(3),YU(3)
C     ..
C     .. External Subroutines ..
      EXTERNAL DXDT,F,FZ,INTEGT,GERK,GERKZ,SETTHU,UV,UVPHI,WR
C     ..
C     .. External Functions ..
      DOUBLE PRECISION P,Q,W
      EXTERNAL P,Q,W
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ATAN2,EXP,LOG,MIN,SIGN
C     ..
C     .. Common blocks ..
      COMMON /DATAF/EIG,IND
      COMMON /ZEE/Z
      COMMON /PIE/PI,TWOPI,HPI
      COMMON /TEMP/TT,YY,NT
      COMMON /TSAVE/TSAVEL,TSAVER
C     ..
C
C     Note: The input values of THEND and DTHDAA are overwritten
C           when integrating from a limit circle endpoint.
C
      PREC = 'DOUBLE PRECISION'
      IFLAG = 1
      IND = 1
      IF (OSC) THEN
         IF (.NOT.OK) THEN
            LOGIC = .FALSE.
         ELSE IF (TEND.LE.TMID) THEN
            LOGIC = TEND.GE.TSAVEL
         ELSE
            LOGIC = TEND.LE.TSAVER
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
               ZSAV = Z
               Z = 1.0
               T = TEND
               PHI0 = ATAN2(COEF2,COEF1)
               IF (TMID.GT.TEND) THEN
                  J = 1
               ELSE
                  J = 2
                  END IF
C
C     We want -PI/2 .lt. PHI0 .le. PI/2.
C
               IF (COEF1.LT.0.0) PHI0 = PHI0 - SIGN(PI,COEF2)
               Y(1) = PHI0
               Y(2) = 0.0
               Y(3) = 0.0
               CALL DXDT(T,TMP,XT0)
               CALL FZ(T,Y,YP)
               DPHIDAA = -YP(1)
               CALL UV(XT0,U,PUP,V,PVP,HU,HV)
               D = U*PVP - V*PUP
C
C     Set THU0 and THV0.
C
               THU0 = ATAN2(U,PUP)
               IF (U.LT.0.0) THU0 = THU0 + TWOPI
               CALL SETTHU(XT0,THU0)
               THV0 = ATAN2(V,PVP)
               IF (V.LT.0.0) THV0 = THV0 + TWOPI
   10          CONTINUE
                  IF (THV0.LT.THU0) THEN
                     THV0 = THV0 + TWOPI
                     GO TO 10
                     END IF
C
C     Set TH0 and copy into THEND, overwriting its input value.
C     Also, redefine DTHDAA, overwriting its input value.
C
               CALL UVPHI(U,PUP,V,PVP,THU0,THV0,PHI0,TH0)
               THEND = TH0
               C = COS(PHI0)
               S = SIN(PHI0)
               YZ0 = U*C + V*S
               PYPZ0 = PUP*C + PVP*S
               DUM = ABS(COS(TH0))
               IF (DUM.GE.0.5) THEN
                  FAC2 = -D*(DUM/PYPZ0)**2
               ELSE
                  FAC2 = -D*(SIN(TH0)/YZ0)**2
                  END IF
               DTHDAA = DPHIDAA*FAC2
               TOUT = TMID
               I = 2
               TT(I,J) = T
               YY(I,1,J) = Y(1)
               YY(I,2,J) = Y(2)
               YY(I,3,J) = Y(3)
               KFLAG = -1
C
C     TSAVEL, TSAVER presumed set by an earlier call with EIG .ne. 0.
C
               IF (EIG.EQ.0.0) THEN
                  IF (T.LE.TSAVEL) TOUT = TSAVEL
                  IF (T.GT.TSAVER) TOUT = TSAVER
                  KFLAG = 1
                  END IF
   20          CONTINUE
                  CALL GERK(FZ,3,Y,T,TOUT,EPS,EPS,KFLAG,ER,WORK,IWORK)
                  IF (KFLAG.GT.3) THEN
                     WRITE(*,*) ' KFLAG1 = 5 '
                     IFLAG = 5
                     RETURN
                     END IF
                  IF (KFLAG.EQ.3) GO TO 20
                  IF (Y(3).LT.-15.0)  Y(3) = -15.0
                  PHI = Y(1)
C
C     Store up to seven values of (T,PHI) for later reference.
C
                  I = I + 1
                  IF (I.LE.7) THEN
                     TT(I,J) = T
                     YY(I,1,J) = PHI
                     YY(I,2,J) = Y(2)
                     YY(I,3,J) = Y(3)
                     END IF
                  IF (10.0*ABS(PHI-PHI0).GE.PI) GO TO 30
                  IF (KFLAG.EQ.-2) GO TO 20
   30          CONTINUE
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
C     Set THU and THV.
C
               THU = ATAN2(U,PUP)
               IF (U.LT.0.0) THU = THU + TWOPI
               CALL SETTHU(XT,THU)
               THV = ATAN2(V,PVP)
               IF (V.LT.0.0) THV = THV + TWOPI
   40          CONTINUE
                  IF (THV.LT.THU) THEN
                    THV = THV + TWOPI
                    GO TO 40
                    END IF
C
C     Now define TH in terms of PHI, THU, and THV.
C
               CALL UVPHI(U,PUP,V,PVP,THU,THV,PHI,TH)
               YZ = U*COS(PHI) + V*SIN(PHI)
               PYPZ = PUP*COS(PHI) + PVP*SIN(PHI)
               Y(1) = TH
               S = SIN(TH)
               C = COS(TH)
               DUM = ABS(C)
               IF (DUM.GE.0.5) THEN
                  FAC2 = -D*(DUM/PYPZ)**2
               ELSE
                  FAC2 = -D*(S/YZ)**2
                  END IF
               DTHIN = FAC2*DPHIDE
               Y(2) = DTHIN
               Z = ZSAV
               RHOSQ = EXP(2.0*Y(3))*(YZ**2+PYPZ**2)
               EFF = 0.5*LOG(RHOSQ*(S**2+(C/Z)**2))
               Y(3) = EFF
C              END (INTEGRATE-FOR-PHI-OSC)
            END IF
         IF (TINTHZ.NE.TMID) THEN
            TIN = TINTHZ
            TOUT = TMID
            THIN = TH
C           DO (INTEGRATE-FOR-THETAZ)
C
C     The following block of code replaces what was used for
C     non-limit-circle problems, with a few changes.
C
               T = TIN
               Y(1) = THIN
               Y(2) = DTHIN
               Y(3) = 0.0
   50          CONTINUE
                  CALL GERKZ(F,3,Y,T,TOUT,EPS,EPS,LFLAG,ERZ,WORK,IWORK)
                  IF (LFLAG.GT.3) THEN
                     WRITE(21,*) ' LFLAGZ = ',LFLAG
                     IFLAG = 5
                     RETURN
                     END IF
                  IF (LFLAG.EQ.3 .OR. LFLAG.EQ.-2) GO TO 50
C              END (INTEGRATE-FOR-THETAZ)
            Y(3) = Y(3) + EFF
            Z = 1.0
            END IF
      ELSE IF (LCIRC) THEN
C        DO (INTEGRATE-FOR-PHI-NONOSC)
            T = TEND
            IF (COEF1.LT.0.0) THEN
               COEF1 = -COEF1
               COEF2 = -COEF2
               END IF
            PHI0 = ATAN2(COEF2,COEF1)
C
C     We want -PI/2 .lt. PHI0 .le. PI/2.
C
            Y(1) = PHI0
            Y(2) = 0.0
            Y(3) = 0.0
            CALL DXDT(T,TMP,XT)
            CALL FZ(T,Y,YP)
            DPHIDAA = -YP(1)
            CALL UV(XT,U,PUP,V,PVP,HU,HV)
            D = U*PVP - V*PUP
            YZ0 = U*COEF1 + V*COEF2
            PYPZ0 = (PUP*COEF1+PVP*COEF2)/Z
C
C     Set TH0 and copy into THEND, overwriting its input value.
C     Also, redefine DTHDAA, overwriting its input value.
C
            TH0 = ATAN2(YZ0,PYPZ0)
            IF (YZ0.LT.0.0) TH0 = TH0 + TWOPI
            THEND = TH0
            DUM = ABS(COS(TH0))
            IF (DUM.GE.0.5) THEN
               FAC2 = (-D/Z)*(DUM/PYPZ0)**2
            ELSE
               FAC2 = (-D/Z)*(SIN(TH0)/YZ0)**2
               END IF
            DTHDAA = FAC2*DPHIDAA
C
C     In the next piece, we assume TH0 .ge. 0.
C
            M = 0
            IF (TH0.EQ.0.0) M = -1
            IF (TH0.GT.PI .OR. (TH0.EQ.PI. AND. T.LT.TMID)) M = 1
            PHI = PHI0
            K2PI = 0
            YZ0 = U*COS(PHI0) + V*SIN(PHI0)
            IF (TMID.GT.TEND) THEN
               J = 1
               TSTAR = -0.99999
               IF (PREC(1:1).NE.'R') TSTAR = -0.9999999999D0
            ELSE
               J = 2
               TSTAR = 0.99999
               IF (PREC(1:1).NE.'R') TSTAR = 0.9999999999D0
               END IF
            DPHIDE = 0.0
            CALL WR(FZ,EPS,TSTAR,PHI0,PHI0,DPHIDE,TOUT,Y,
     1              TT(1,J),YY(1,1,J),ERZ,WORK,IWORK)
            T = TOUT
            DDD = MIN(0.01,ABS(TMID-TOUT))
            TOUT = TOUT + DDD*(TMID-TOUT)/ABS(TMID-TOUT)
            KFLAG = 1
            CALL GERK(FZ,3,Y,T,TOUT,EPS,EPS,KFLAG,ER,WORK,IWORK)
            TT(7,J) = T
            YY(7,1,J) = Y(1)
            YY(7,2,J) = Y(2)
            YY(7,3,J) = Y(3)
   60       CONTINUE
               T = TOUT
               CALL DXDT(T,TMP,XT)
               CALL UV(XT,U,PUP,V,PVP,HU,HV)
               PHI = Y(1)
               S = SIN(PHI)
               C = COS(PHI)
               YZ = U*C + V*S
               IF (YZ*YZ0 .LT. 0.0) K2PI = K2PI + 1
               YZ0 = YZ
               IF (KFLAG.GT.3) THEN
                  WRITE(*,*) ' KFLAG2 = ',KFLAG
                  IFLAG = 5
                  RETURN
                  END IF
               IF (KFLAG.EQ.3 .OR. KFLAG.EQ.-2) GO TO 60
C
C     Convert from PHI to THETA.
C
            DPHIDE = Y(2)
            D = U*PVP-V*PUP
            PYPZ = (PUP*C+PVP*S)/Z
            THBAR = ATAN2(YZ,PYPZ)
            IF (TMID.GT.TEND .AND. THBAR.LT.TH0 .AND. PHI.LT.PHI0)
     1         THBAR = THBAR + TWOPI
            IF (TMID.LT.TEND .AND. THBAR.GT.TH0 .AND. PHI.GT.PHI0)
     1         THBAR = THBAR - TWOPI
            TH = THBAR - M*PI
            IF (TH.LT.-PI) TH = TH + TWOPI
            IF (TH.GT.TWOPI) TH = TH - TWOPI
            IF (TMID.LT.TEND .AND. K2PI.GT.1) TH = TH - (K2PI-1)*TWOPI
            IF (TMID.GT.TEND .AND. K2PI.GT.1) TH = TH + (K2PI-1)*TWOPI
            IF (TMID.GT.TEND .AND. TH*TH0.LT.0.0) TH = TH + TWOPI
C
C     We now have YZ, PYPZ, PHI and TH.
C
            DUM = ABS(COS(TH))
            IF (DUM.GE.0.5) THEN
               FAC2 = -(D/Z)*(DUM/PYPZ)**2
            ELSE
               FAC2 = -(D/Z)*(SIN(TH)/YZ)**2
               END IF
            DTHIN = FAC2*DPHIDE
            THETA = ATAN2(YZ,Z*PYPZ)
            S = SIN(THETA)
            C = COS(THETA)
            RHOSQ = EXP(2.0*Y(3))*(YZ**2+(Z*PYPZ)**2)
            EFF = 0.5*LOG(RHOSQ*(S**2+(C/Z)**2))
C           END (INTEGRATE-FOR-PHI-NONOSC)
         TIN = TOUT
         TOUT = TMID
         THIN = TH
C        DO (INTEGRATE-FOR-THETA)
            CALL INTEGT(TIN,TOUT,THIN,DTHIN,EPS,Y,IFLAG,ER,WORK,IWORK)
C           END (INTEGRATE-FOR-THETA)
         Y(3) = Y(3) + EFF
      ELSE IF (.NOT.OK .AND. .NOT.SING) THEN
C
C     This is the 'weakly regular' case.
C
         IF (TMID.GT.TEND) THEN
            J = 1
            TSTAR = -0.99999
            IF (PREC(1:1).NE.'R') TSTAR = -0.9999999999D0
         ELSE
            J = 2
            TSTAR = 0.99999
            IF (PREC(1:1).NE.'R') TSTAR = 0.9999999999D0
            END IF
         CALL DXDT(TSTAR,TMP,XSTAR)
         P1 = 1.0/P(XSTAR)
         Q1 = Q(XSTAR)
         W1 = W(XSTAR)
         YSTAR = THEND + 0.5*(TSTAR-TEND)*
     1           (P1*COS(THEND)**2+(EIG*W1-Q1)*SIN(THEND)**2)
         CALL WR(F,EPS,TSTAR,YSTAR,THEND,DTHDE,TOUT,YU,
     1           TT(1,J),YY(1,1,J),ERZ,WORK,IWORK)
         T = TOUT
         DDD = MIN(0.01,ABS(TMID-TOUT))
         TOUT = TOUT + DDD*(TMID-TOUT)/ABS(TMID-TOUT)
         KFLAG = 1
         CALL GERK(F,3,YU,T,TOUT,EPS,EPS,KFLAG,ER,WORK,IWORK)
         TT(7,J) = T
         YY(7,1,J) = YU(1)
         YY(7,2,J) = YU(2)
         YY(7,3,J) = YU(3)
         TIN = TOUT
         TOUT = TMID
         THIN = YU(1)
         DTHIN = YU(2)
C        DO (INTEGRATE-FOR-THETA)
            CALL INTEGT(TIN,TOUT,THIN,DTHIN,EPS,Y,IFLAG,ER,WORK,IWORK)
C           END (INTEGRATE-FOR-THETA)
      ELSE
C
C     This is the regular (not weakly regular) or limit point case.
C
         TIN = TEND
         TOUT = TMID
         THIN = THEND
         DTHIN = DTHDE
C        DO (INTEGRATE-FOR-THETA)
            CALL INTEGT(TIN,TOUT,THIN,DTHIN,EPS,Y,IFLAG,ER,WORK,IWORK)
C           END (INTEGRATE-FOR-THETA)
         END IF
      RETURN
      END
      SUBROUTINE INTEGT(TIN,TOUT,THIN,DTHIN,EPS,Y,IFLAG,ER,WORK,IWORK)
      INTEGER IFLAG,IWORK(5)
      DOUBLE PRECISION TIN,TOUT,THIN,DTHIN,EPS
      DOUBLE PRECISION Y(3),ER(3),WORK(27)
C     **********
C
C     Integrate for theta.
C
C     **********
C     .. Local Scalars ..
      INTEGER LFLAG
      DOUBLE PRECISION T
C     ..
C     .. External Subroutines ..
      EXTERNAL F
C     ..
C     DO (INTEGRATE-FOR-TH)
         T = TIN
         Y(1) = THIN
         Y(2) = DTHIN
         Y(3) = 0.0
         LFLAG = 1
   10    CONTINUE
            CALL GERK(F,3,Y,T,TOUT,EPS,EPS,LFLAG,ER,WORK,IWORK)
            IF (LFLAG.EQ.3) GO TO 10
         IF (LFLAG.GT.3) THEN
            WRITE(21,*) ' LFLAG = ',LFLAG
            IFLAG = 5
            RETURN
            END IF
         THIN = Y(1)
         DTHIN = Y(2)
C        END (INTEGRATE-FOR-TH)
      RETURN
      END
      SUBROUTINE INTPOL(N,XN,FN,X,ABSERR,MAXDEG,ANS,ERROR)
      INTEGER N,MAXDEG
      DOUBLE PRECISION X,ABSERR,ANS,ERROR
      DOUBLE PRECISION XN(N),FN(N)
C     **********
C
C     This subroutine forms an interpolating polynomial for data pairs.
C     It is called from EXTRAP in solving a Sturm-Liouville problem.
C
C     **********
C     .. Local Scalars ..
      INTEGER I,I1,II,IJ,IK,IKM1,J,K,L,LIMIT
      DOUBLE PRECISION PROD
C     ..
C     .. Local Arrays ..
      INTEGER INDEX(10)
      DOUBLE PRECISION V(10,10)
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS,MIN
C     ..
      L = MIN(MAXDEG,N-2) + 2
      LIMIT = MIN(L,N-1)
      DO 10 I = 1,N
         V(I,1) = ABS(XN(I)-X)
         INDEX(I) = I
   10    CONTINUE
      DO 30 I=1,LIMIT
         DO 20 J=I+1,N
            II = INDEX(I)
            IJ = INDEX(J)
            IF (V(II,1).GT.V(IJ,1)) THEN
               INDEX(I) = IJ
               INDEX(J) = II
               END IF
   20       CONTINUE
   30    CONTINUE
      PROD = 1.0
      I1 = INDEX(1)
      ANS = FN(I1)
      V(1,1) = FN(I1)
      DO 50 K=2,L
         IK = INDEX(K)
         V(K,1) = FN(IK)
         DO 40 I=1,K-1
            II = INDEX(I)
            V(K,I+1) = (V(I,I)-V(K,I))/(XN(II)-XN(IK))
   40       CONTINUE
         IKM1 = INDEX(K-1)
         PROD = (X-XN(IKM1))*PROD
         ERROR = PROD*V(K,K)
         IF(ABS(ERROR).LE.ABSERR) RETURN
         ANS = ANS + ERROR
   50    CONTINUE
      ANS = ANS - ERROR
      RETURN
      END
      SUBROUTINE LPOL(KEIGS,XN,FN,X,F)
      INTEGER KEIGS
      DOUBLE PRECISION X,F
      DOUBLE PRECISION XN(KEIGS),FN(KEIGS)
C     **********
C     **********
C     ..
C     .. Local variables ..
      DOUBLE PRECISION X12,X13,X14,X15,X21,X23,X24,X25,X31,X32,X34,X35,
     1     X41,X42,X43,X45,X51,X52,X53,X54,XX1,XX2,XX3,XX4,XX5
C     ..
      X12 = XN(1) - XN(2)
      X13 = XN(1) - XN(3)
      X23 = XN(2) - XN(3)
      X21 = -X12
      X31 = -X13
      X32 = -X23
      XX1 = X - XN(1)
      XX2 = X - XN(2)
      XX3 = X - XN(3)
      IF (KEIGS.GE.4) THEN
         X14 = XN(1) - XN(4)
         X24 = XN(2) - XN(4)
         X34 = XN(3) - XN(4)
         X41 = -X14
         X42 = -X24
         X43 = -X34
         XX4 = X - XN(4)
         END IF
      IF (KEIGS.EQ.5) THEN
         X15 = XN(1) - XN(5)
         X25 = XN(2) - XN(5)
         X35 = XN(3) - XN(5)
         X45 = XN(4) - XN(5)
         X51 = -X15
         X52 = -X25
         X53 = -X35
         X54 = -X45
         XX5 = X - XN(5)
         END IF
      IF (KEIGS.EQ.3) THEN
         F = XX2*XX3*FN(1)/(X12*X13) +
     1       XX1*XX3*FN(2)/(X21*X23) +
     2       XX1*XX2*FN(3)/(X31*X32)
      ELSE IF (KEIGS.EQ.4) THEN
         F = XX2*XX3*XX4*FN(1)/(X12*X13*X14) +
     1       XX1*XX3*XX4*FN(2)/(X21*X23*X24) +
     2       XX1*XX2*XX4*FN(3)/(X31*X32*X34) +
     3       XX1*XX2*XX3*FN(4)/(X41*X42*X43)
      ELSE
         F = XX2*XX3*XX4*XX5*FN(1)/(X12*X13*X14*X15) +
     1       XX1*XX3*XX4*XX5*FN(2)/(X21*X23*X24*X25) +
     2       XX1*XX2*XX4*XX5*FN(3)/(X31*X32*X34*X35) +
     3       XX1*XX2*XX3*XX5*FN(4)/(X41*X42*X43*X45) +
     4       XX1*XX2*XX3*XX4*FN(5)/(X51*X52*X53*X54)
         END IF
      RETURN
      END
C
      SUBROUTINE PERIO(A,B,INTAB,P0ATA,QFATA,P0ATB,QFATB,
     1                 A1,A2,B1,B2,NUMEIG,EIG,TOL,IFLAG,SLFN,
     2                 SINGATA,SINGATB,CIRCLA,CIRCLB,OSCILA,OSCILB)
      INTEGER INTAB,NUMEIG,IFLAG
      DOUBLE PRECISION A,B,P0ATA,QFATA,P0ATB,QFATB,A1,A2,B1,B2,EIG,TOL,
     1     SINGATA,SINGATB,CIRCLA,CIRCLB,OSCILA,OSCILB
      DOUBLE PRECISION SLFN(9)
C     **********
C     **********
C     .. Scalars in Common ..
      INTEGER IND
      DOUBLE PRECISION EPSMIN,CC,UB,UL,UR,VB,VL,VR,Z
C     ..
C     .. Local Scalars ..
      INTEGER JFLAG
      DOUBLE PRECISION A1D,A1N,A2D,A2N,B1N,B1D,B2N,B2D
      DOUBLE PRECISION AE,EIGLO,EIGUP,LAMBDA,LAMUP,RE,TOLL,TOLS
C     ..
C     .. External Subroutines ..
      EXTERNAL FZERO,SLEIGN2
C     ..
C     .. External Functions ..
      EXTERNAL FF
C     ..
C     .. Intrinsic Functions ..
C     ..
C     .. Common blocks ..
      COMMON /EPP2/CC,UL,UR,VL,VR,UB,VB,IND
      COMMON /RNDOFF/EPSMIN
      COMMON /ZEE/Z
C     ..
C
C     Get upper and lower bounds as accurately as possible.
C
      TOLL = 0.0
      TOLS = TOL
C
      A1N = 0.0
      A2N = 1.0
      B1N = 0.0
      B2N = 1.0
      TOL = TOLL
      EIG = 0.0
      CALL SLEIGN2(A,B,INTAB,P0ATA,QFATA,P0ATB,QFATB,A1N,A2N,
     1             B1N,B2N,NUMEIG,EIG,TOL,IFLAG,0,SLFN,
     2             SINGATA,SINGATB,CIRCLA,CIRCLB,OSCILA,OSCILB)
      WRITE(*,*) ' eiglo,iflag = ',EIG,IFLAG
      EIGLO = EIG - 0.1
C
      A1D = 1.0
      A2D = 0.0
      B1D = 1.0
      B2D = 0.0
      TOL = TOLL
      EIG = 0.0
      CALL SLEIGN2(A,B,INTAB,P0ATA,QFATA,P0ATB,QFATB,A1D,A2D,
     1             B1D,B2D,NUMEIG,EIG,TOL,IFLAG,0,SLFN,
     2             SINGATA,SINGATB,CIRCLA,CIRCLB,OSCILA,OSCILB)
      WRITE(*,*) ' eigup,iflag = ',EIG,IFLAG
      EIGUP = EIG + 0.1
      LAMBDA = 0.5*(EIGLO+EIGUP)
C
C     The following call to SLEIGN2 sets the stage for INTEG.
C
      TOL = .001
      CALL SLEIGN2(A,B,INTAB,P0ATA,QFATA,P0ATB,QFATB,A1,A2,
     1             B1,B2,NUMEIG,LAMBDA,TOL,IFLAG,-1,SLFN,
     2             SINGATA,SINGATB,CIRCLA,CIRCLB,OSCILA,OSCILB)
      Z = 1.0
C
      EIG = EIGLO
      LAMUP = EIGUP
      RE = TOLS
      AE = RE
C
C     IND = 2 selects the classical function d(LAMBDA) in function FF.
C
      IND = 2
      CALL FZERO(FF,EIG,LAMUP,LAMBDA,RE,AE,JFLAG)
      IF (JFLAG.NE.1) WRITE(*,*) ' jflag = ',JFLAG
      IFLAG = JFLAG
      TOL = ABS(EIG-LAMUP)/MAX(1.0,ABS(EIG))
      WRITE(*,*) ' EIGLO,EIGUP = ',EIGLO,EIGUP
C
      RETURN
      END
      SUBROUTINE SETMID(MF,ML,EIG,QS,WS,IMID,TMID)
      INTEGER MF,ML,IMID
      DOUBLE PRECISION EIG,QS(*),WS(*),TMID
C     **********
C
C     This procedures tests the interval sample points in the order
C     50,51,49,52,48,...,etc. for the first one where the expression
C     (lambda*w-q) is positive.  This point is designated TMID.
C
C     **********
C     .. Local Scalars ..
      INTEGER I,J
      DOUBLE PRECISION S
C     ..
C     .. External Functions ..
      DOUBLE PRECISION TFROMI
      EXTERNAL TFROMI
C     ..
      S = -1.0
      DO 10 J=1,100
         I = 50 + S*(J/2)
         S = -S
         IF (I.LT.MF .OR. I.GT.ML) GO TO 20
         IF (EIG*WS(I)-QS(I).GT.0.0) THEN
            IMID = I
            TMID = TFROMI(IMID)
            GO TO 20
            END IF
   10    CONTINUE
   20 CONTINUE
      WRITE(*,*) ' new tmid = ',TMID
      RETURN
      END
      SUBROUTINE SETTHU(X,THU)
      DOUBLE PRECISION X,THU
C     **********
c
c  This subroutine establishes a definite value for THU,
c    the phase angle for the function U, including an
c    appropriate integer multiple of pi
c  It needs the numbers MMW(*) found in THUM
c
C     **********
C     .. Scalars in Common ..
      INTEGER MMWD
      DOUBLE PRECISION PI,TWOPI,HPI
C     ..
C     .. Arrays in Common ..
      INTEGER MMW(98)
      DOUBLE PRECISION YS(197)
C     ..
C     .. Local Scalars ..
      INTEGER I
C     ..
C     .. Common blocks ..
      COMMON /PASS/YS,MMW,MMWD
      COMMON /PIE/PI,TWOPI,HPI
C     ..
      DO 10 I=1,MMWD
         IF (X.GE.YS(MMW(I)) .AND. X.LE.YS(MMW(I)+1)) THEN
            IF (THU.GT.PI) THEN
               THU = THU + (I-1)*TWOPI
               RETURN
            ELSE
               THU = THU + I*TWOPI
               RETURN
               END IF
            END IF
   10    CONTINUE
      DO 20 I=1,MMWD
         IF (X.GE.YS(MMW(I))) THU = THU + TWOPI
   20    CONTINUE
      RETURN
      END
      DOUBLE PRECISION FUNCTION TFROMI(I)
      INTEGER I
C     **********
C
C     This function associates the value of an interval sample point
C     with its index.
C
C     **********
      IF (I.LT.8) THEN
         TFROMI = -1.0 + 0.1/4.0**(8-I)
      ELSE IF (I.GT.92) THEN
         TFROMI = 1.0 - 0.1/4.0**(I-92)
      ELSE
         TFROMI = 0.0227*(I-50)
         END IF
      RETURN
      END
      SUBROUTINE THTOTHZ(Y,Z,U)
      DOUBLE PRECISION Z
      DOUBLE PRECISION Y(3),U(3)
C     **********
C  THIS PROGRAM CONVERTS FROM TH TO THZ WHERE
C    TAN(TH)=PSI/(P*PSI')
C  AND
C    TAN(THZ)=Z*PSI/(P*PSI'),
C  OR
C    TAN(THZ)=Z*TAN(TH) .
C  SO WE HAVE
C    DTHZ=Z*(COS(THZ)/COS(TH))**2 * DTH  ,
C  OR
C    DTH=(1/Z)*(COS(TH)/COS(THZ))**2 * DTHZ .
C     **********
C     .. Scalars in Common ..
      DOUBLE PRECISION PI,TWOPI,HPI
C     ..
C     .. Local Scalars ..
      INTEGER K
      DOUBLE PRECISION DTH,DTHZ,DUM,FAC,PIK,REMTH,TH,THZ
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS,ATAN,COS,LOG,SIN,TAN
C     ..
C     .. Common blocks ..
      COMMON /PIE/PI,TWOPI,HPI
C     ..
      TH = Y(1)
      DTH = Y(2)
      K = TH/PI
      IF (TH.LT.0.0) K = K - 1
      PIK = K*PI
      REMTH = TH - PIK
      IF (4.0*REMTH.LE.PI) THEN
         THZ = ATAN(Z*TAN(REMTH)) + PIK
      ELSE IF (4.0*REMTH.GE.3.0*PI) THEN
         THZ = ATAN(Z*TAN(REMTH)) + PIK + PI
      ELSE
         THZ = ATAN(TAN(REMTH-HPI)/Z) + PIK + HPI
         END IF
      DUM = ABS(COS(THZ))
      IF (DUM.GE.0.5) THEN
         FAC = Z*(DUM/COS(TH))**2
      ELSE
         FAC = (SIN(THZ)/SIN(TH))**2/Z
         END IF
      DTHZ = FAC*DTH
      U(1) = THZ
      U(2) = DTHZ
      U(3) = Y(3) - 0.5*LOG(Z*FAC)
      RETURN
      END
      SUBROUTINE THUM(MF,ML,XS)
      INTEGER MF,ML
      DOUBLE PRECISION XS(*)
C     **********
C
C  YS IS LIKE XS, BUT HAS TWICE AS MANY POINTS.
C  MMW(N) IS THE VALUE OF THE INDEX I OF U(I), MF .LE. I .LE. 2*ML-2,
C    WHERE U FOR THE Nth TIME CHANGES SIGN FROM - TO +
C    AND WHERE P*U' IS POSITIVE.
C  MMWD IS THE NUMBER OF SUCH POINTS OF U.
C
C     **********
C     .. Scalars in Common ..
      INTEGER MMWD
C     ..
C     .. Arrays in Common ..
      INTEGER MMW(98)
      DOUBLE PRECISION YS(197)
C     ..
C     .. Local Scalars ..
      INTEGER I,N
      DOUBLE PRECISION PUP,PUP1,U,U1,tmp
C     ..
C     .. Common blocks ..
      COMMON /PASS/YS,MMW,MMWD
C     ..
      DO 10 I=1,98
         YS(2*I-1) = XS(I)
         YS(2*I) = 0.5*(XS(I)+XS(I+1))
  10     CONTINUE
      YS(197) = XS(99)
      N = 0
      U1 = 0.0
      PUP1 = 0.0
      DO 20 I=2*MF-1,2*ML-1
         CALL UV(YS(I),U,PUP,tmp,tmp,tmp,tmp)
         IF (U1.LT.0.0 .AND. U.GT.0.0 .AND. PUP1.GT.0.0) THEN
            N = N + 1
            MMW(N) = I - 1
            END IF
         U1 = U
         PUP1 = PUP
  20     CONTINUE
      MMWD = N
      RETURN
      END
      SUBROUTINE THZTOTH(U,Z,Y)
      DOUBLE PRECISION Z
      DOUBLE PRECISION U(3),Y(3)
C     **********
C  THIS PROGRAM CONVERTS FROM THZ TO TH WHERE
C    TAN(TH)=PSI/(P*PSI')
C  AND
C    TAN(THZ)=Z*PSI/(P*PSI'),
C  OR
C    TAN(THZ)=Z*TAN(TH) .
C  SO WE HAVE
C    DTHZ=Z*(COS(THZ)/COS(TH))**2 * DTH  ,
C  OR
C    DTH=(1/Z)*(COS(TH)/COS(THZ))**2 * DTHZ .
C     **********
C     .. Scalars in Common ..
      DOUBLE PRECISION PI,TWOPI,HPI
C     ..
C     .. Local Scalars ..
      INTEGER K
      DOUBLE PRECISION DTH,DTHZ,DUM,FAC,PIK,REMTHZ,TH,THZ
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS,ATAN,COS,LOG,SIN,TAN
C     ..
C     .. Common blocks ..
      COMMON /PIE/PI,TWOPI,HPI
C     ..
      THZ = U(1)
      DTHZ = U(2)
      K = THZ/PI
      IF (THZ.LT.0.0) K = K - 1
      PIK = K*PI
      REMTHZ = THZ - PIK
      IF (4.0*REMTHZ.LE.PI) THEN
         TH = ATAN(TAN(REMTHZ)/Z) + PIK
      ELSE IF (4.0*REMTHZ.GE.3.0*PI) THEN
         TH = ATAN(TAN(REMTHZ)/Z) + PIK + PI
      ELSE
         TH = ATAN(Z*TAN(REMTHZ-HPI)) + PIK + HPI
         END IF
      DUM = ABS(COS(TH))
      IF (DUM.GE.0.5) THEN
         FAC = (DUM/COS(THZ))**2/Z
      ELSE
         FAC = Z*(SIN(TH)/SIN(THZ))**2
         END IF
      DTH = FAC*DTHZ
      Y(1) = TH
      Y(2) = DTH
      Y(3) = U(3) + 0.5*LOG(Z/FAC)
      RETURN
      END
      SUBROUTINE UVPHI(U,PUP,V,PVP,THU,THV,PHI,TH)
      DOUBLE PRECISION U,PUP,V,PVP,THU,THV,PHI,TH
C     **********
C
C     This program finds TH appropriate to THU, THV, and PHI, where
C     THU is the phase angle for U, and THV is the phase angle for V.
C
C     **********
C     .. Scalars in Common ..
      DOUBLE PRECISION PI,TWOPI,HPI
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION C,D,PYP,S,Y
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
      TH = THU
      IF (PHI.EQ.0.0) RETURN
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
      IF (Y.LT.0.0) TH = TH + TWOPI
      D = U*PVP - V*PUP
      IF (D*PHI.GT.0.0) THEN
         CALL FIT(THU-PI,TH,THU)
      ELSE
         CALL FIT(THU,TH,THU+PI)
         END IF
      RETURN
      END
      SUBROUTINE WR(FG,EPS,TSTAR,YSTAR,THEND,DTHDE,TOUT,Y,
     1              TT,YY,ERR,WORK,IWORK)
      INTEGER IWORK(5)
      DOUBLE PRECISION EPS,TSTAR,YSTAR,THEND,DTHDE,TOUT
      DOUBLE PRECISION Y(3),TT(7),YY(7,3),ERR(3),WORK(27)
      EXTERNAL FG
C     **********
C
C     This subroutine integrates Y' = F(T,Y) from t = +/-1.0 to THEND
C     even when F cannot be evaluated at t.  (T*,Y*) is chosen as a
C     nearby point, and the equation is integrated from there and
C     checked for consistency with having integrated from t.  If not,
C     a different (T*,Y*) is chosen until consistency is achieved.
C
C     **********
C     .. Scalars in Common ..
      INTEGER IND
      DOUBLE PRECISION EIG,EPSMIN
C     ..
C     .. Local Scalars ..
      INTEGER I,K,KFLAG
      DOUBLE PRECISION CHNG,D2F,D2G,D3F,D3G,D4F,D4G,HT,HU,OLDSS2,OLDYY2,
     1     ONE,SLO,SOUT,SUMM,SUP,T,TEN5,TIN,U,UOUT,USTAR,YLO,YOUT,YUP
C     ..
C     .. Local Arrays ..
      DOUBLE PRECISION DF(4),DG(4),FF(6),GG(5),S(3),SS(6,3),UU(6)
C     ..
C     .. External Subroutines ..
      EXTERNAL GERK,LPOL
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS,SIGN
C     ..
C     .. Common blocks
      COMMON /DATAF/EIG,IND
      COMMON /RNDOFF/EPSMIN
C     ..
      ONE = 1.0
      TEN5 = 100000.0
C
C     Integrate Y' = F(T,Y,YP).
C
      TIN = SIGN(ONE,TSTAR)
      HT = TSTAR - TIN
      TT(1) = TIN
      YY(1,1) = THEND
      YY(1,2) = DTHDE
      YY(1,3) = 0.0
      TT(2) = TSTAR
      YY(2,1) = YSTAR
      YY(2,2) = DTHDE
      YY(2,3) = 0.0
      YLO = -TEN5
      YUP = TEN5
C
C     Normally IND = 1 or 3;  IND is set to 2 or 4 when Y is to be used
C     as the independent variable in FG(T,Y,YP), instead of the usual T.
C     Before leaving this subroutine, IND is reset to 1.
C
   10 CONTINUE
         T = TSTAR
         Y(1) = YY(2,1)
         Y(2) = YY(2,2)
         Y(3) = YY(2,3)
         KFLAG = 1
         IND = 1
         DO 30 K = 3,6
            TOUT = T + HT
   20       CONTINUE
               CALL GERK(FG,3,Y,T,TOUT,EPS,EPS,KFLAG,ERR,WORK,IWORK)
               IF (KFLAG.GT.3) THEN
                  WRITE(*,*) ' KFLAG3 = 5 '
                  END IF
               IF (KFLAG.EQ.3) GO TO 20
            YOUT = Y(1)
            TT(K) = T
            YY(K,1) = YOUT
            YY(K,2) = Y(2)
            YY(K,3) = Y(3)
   30       CONTINUE
         IND = 3
         DO 40 I = 2,6
            CALL FG(TT(I),YY(I,1),FF(I))
   40       CONTINUE
         CALL LPOL(5,TT(2),FF(2),TT(1),FF(1))
         IF (ABS(FF(1)).LE.500.0) THEN
C
C     Now we want to apply some criterion to see if these results are
C     consistent with having integrated from (TT(1),YY(1,1).
C
            DO 50 I = 1,4
               DF(I) = FF(I+1) - FF(I)
   50          CONTINUE
            D2F = DF(4) - DF(3)
            D3F = DF(4) - 2.0*DF(3) + DF(2)
            D4F = DF(4) - 3.0*DF(3) + 3.0*DF(2) - DF(1)
            SUMM = HT*(FF(5)-3.5*DF(4)+53.0*D2F/12.0
     1                 -55.0*D3F/24.0+251.0*D4F/720.0)
C
C     Presumably, YY(2,1) should be YY(1,1) + SUMM.
C
            OLDYY2 = YY(2,1)
            YY(2,1) = YY(1,1) + SUMM
C
C     Also improve the value of Y(2) at TSTAR.
C
            YY(2,2) = 0.5*(YY(1,2)+YY(3,2))
            YY(2,3) = 0.5*(YY(1,3)+YY(3,3))
            CHNG = YY(2,1) - OLDYY2
            IF (CHNG.GE.0.0 .AND. OLDYY2.GT.YLO) YLO = OLDYY2
            IF (CHNG.LE.0.0 .AND. OLDYY2.LT.YUP) YUP = OLDYY2
            IF ((YY(2,1).GE.YUP .AND. YLO.GT.-TEN5) .OR.
     1          (YY(2,1).LE.YLO .AND. YUP.LT.TEN5))
     2         YY(2,1) = 0.5*(YLO+YUP)
            IF (ABS(YY(2,1)-OLDYY2).GT.EPSMIN) GO TO 10
         ELSE
C
C     Here, Y' is assumed infinite at T = TIN.  In this case,
C     it cannot be expected to approximate Y with a polynomial,
C     so the independent and dependent variables are interchanged.
C     The points are assumed equally spaced.
C
            HU = (YY(6,1)-YY(1,1))/5.0
            UU(1) = YY(1,1)
            SS(1,1) = TT(1)
            SS(1,2) = DTHDE
            SS(1,3) = 0.0
            UU(2) = UU(1) + HU
            SS(2,1) = TSTAR
            SS(2,2) = DTHDE
            SS(2,3) = 0.0
            USTAR = UU(2)
            SLO = -TEN5
            SUP = TEN5
   60       CONTINUE
               U = USTAR
               S(1) = SS(2,1)
               S(2) = SS(2,2)
               S(3) = SS(2,3)
               KFLAG = 1
               IND = 2
               DO 80 K = 3,6
                  UOUT = U + HU
   70             CONTINUE
                     CALL GERK(FG,3,S,U,UOUT,EPS,EPS,KFLAG,ERR,
     1                         WORK,IWORK)
                     IF (KFLAG.GT.3) THEN
                        WRITE(*,*) ' KFLAG4 = 5 '
                        END IF
                     IF (KFLAG.EQ.3) GO TO 70
                  SOUT = S(1)
                  UU(K) = U
                  SS(K,1) = SOUT
                  SS(K,2) = S(2)
                  SS(K,3) = S(3)
   80             CONTINUE
               IND = 4
               DO 90 I = 2,5
                  CALL FG(UU(I),SS(I,1),GG(I))
   90             CONTINUE
               GG(1) = 0.0
               DO 100 I = 1,4
                  DG(I) = GG(I+1) - GG(I)
  100             CONTINUE
               D2G = DG(4) - DG(3)
               D3G = DG(4) - 2.0*DG(3) + DG(2)
               D4G = DG(4) - 3.0*DG(3) + 3.0*DG(2) - DG(1)
               SUMM = HU*(GG(5)-3.5*DG(4)+53.0*D2G/12.0
     1                 -55.0*D3G/24.0+251.0*D4G/720.0)
C
C     Presumably, SS(2,1) should be SS(1,1) + SUMM.
C
               OLDSS2 = SS(2,1)
               SS(2,1) = SS(1,1) + SUMM
               IF (SS(2,1).LE.-1.0) SS(2,1) = -1.0 + EPSMIN
               IF (SS(2,1).GE.1.0) SS(2,1) = 1.0 - EPSMIN
C
C     Also improve the value of Y(2) at TSTAR.
C
               SS(2,2) = 0.5*(SS(1,2)+SS(3,2))
               SS(2,3) = 0.5*(SS(1,3)+SS(3,3))
               CHNG = SS(2,1) - OLDSS2
               IF (CHNG.GE.0.0 .AND. OLDSS2.GT.SLO) SLO = OLDSS2
               IF (CHNG.LE.0.0 .AND. OLDSS2.LT.SUP) SUP = OLDSS2
               IF ((SS(2,1).GE.SUP .AND. SLO.GT.-TEN5) .OR.
     1             (SS(2,1).LE.SLO .AND. SUP.LT.TEN5))
     2            SS(2,1) = 0.5*(SLO+SUP)
               IF (ABS(SS(2,1)-OLDSS2).GT.EPSMIN) GO TO 60
            END IF
      IF (IND.EQ.4) THEN
         DO 110 I = 1,6
            TT(I) = SS(I,1)
            YY(I,1) = UU(I)
            YY(I,2) = SS(I,2)
            YY(I,3) = SS(I,3)
  110       CONTINUE
         END IF
      TOUT = TT(6)
      IND = 1
      RETURN
      END
      SUBROUTINE GERK(F, NEQN, Y, T, TOUT, RELERR, ABSERR, IFLAG,
     * GERROR, WORK, IWORK)
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
C    DOUBLE PRECISION.
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
C     .. SCALAR ARGUMENTS ..
      INTEGER IFLAG,NEQN
      DOUBLE PRECISION ABSERR,RELERR,T,TOUT
C     ..
C     .. ARRAY ARGUMENTS ..
      INTEGER IWORK(5)
      DOUBLE PRECISION GERROR(NEQN),WORK(3+8*NEQN),Y(NEQN)
C     ..
C     .. SUBROUTINE ARGUMENTS ..
      EXTERNAL F
C     ..
C     .. LOCAL SCALARS ..
      INTEGER K1,K1M,K2,K3,K4,K5,K6,K7,K8
C     ..
C     .. EXTERNAL SUBROUTINES ..
      EXTERNAL GERKS
C     ..
C COMPUTE INDICES FOR THE SPLITTING OF THE WORK ARRAY
      K1M = NEQN + 1
      K1 = K1M + 1
      K2 = K1 + NEQN
      K3 = K2 + NEQN
      K4 = K3 + NEQN
      K5 = K4 + NEQN
      K6 = K5 + NEQN
      K7 = K6 + NEQN
      K8 = K7 + NEQN
C *******************************************************************
C      THIS INTERFACING ROUTINE MERELY RELIEVES THE USER OF A LONG
C      CALLING LIST VIA THE SPLITTING APART OF TWO WORKING STORAGE
C      ARRAYS. IF THIS IS NOT COMPATIBLE WITH THE USERS COMPILER,
C      HE MUST USE GERKS DIRECTLY.
C *******************************************************************
      CALL GERKS(F, NEQN, Y, T, TOUT, RELERR, ABSERR, IFLAG,
     * GERROR, WORK(1), WORK(K1M), WORK(K1), WORK(K2), WORK(K3),
     * WORK(K4), WORK(K5), WORK(K6), WORK(K7), WORK(K8),
     * WORK(K8+1), IWORK(1), IWORK(2), IWORK(3), IWORK(4), IWORK(5))
      RETURN
      END
      SUBROUTINE GERKS(F, NEQN, Y, T, TOUT, RELERR, ABSERR, IFLAG,
     * GERROR, YP, H, F1, F2, F3, F4, F5, YG, YGP, SAVRE, SAVAE,
     * NFE, KOP, INIT, JFLAG, KFLAG)
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
C     .. SCALAR ARGUMENTS ..
      INTEGER IFLAG,INIT,JFLAG,KFLAG,KOP,NEQN,NFE
      DOUBLE PRECISION ABSERR,H,RELERR,SAVAE,SAVRE,T,TOUT
C     ..
C     .. ARRAY ARGUMENTS ..
      DOUBLE PRECISION F1(NEQN),F2(NEQN),F3(NEQN),F4(NEQN),F5(NEQN),
     1     GERROR(NEQN),Y(NEQN),YG(NEQN),YGP(NEQN),YP(NEQN)
C     ..
C     .. SUBROUTINE ARGUMENTS ..
      EXTERNAL F
C     ..
C     .. LOCAL SCALARS ..
      INTEGER K,MAXNFE,MFLAG
      LOGICAL HFAILD,OUTPUT
      DOUBLE PRECISION A,AE,DT,EE,EEOET,ESTTOL,ET,HH,HMIN,ONE,REMIN,RER,
     1     S,SCALE,TOL,TOLN,TS,U,U26,YPK
C     ..
C     .. EXTERNAL FUNCTIONS ..
      DOUBLE PRECISION EPSLON
      EXTERNAL EPSLON
C     ..
C     .. EXTERNAL SUBROUTINES ..
      EXTERNAL FEHL
C     ..
C     .. INTRINSIC FUNCTIONS ..
      INTRINSIC ABS,MAX,MIN,SIGN
C     ..
C *******************************************************************
C   REMIN IS A TOLERANCE THRESHOLD WHICH IS ALSO DETERMINED BY THE
C   INTEGRATION METHOD. IN PARTICULAR, A FIFTH ORDER METHOD WILL
C   GENERALLY NOT BE CAPABLE OF DELIVERING ACCURACIES NEAR LIMITING
C   PRECISION ON COMPUTERS WITH LONG WORDLENGTHS.
      DATA REMIN /3.E-11/
C *******************************************************************
C      THE EXPENSE IS CONTROLLED BY RESTRICTING THE NUMBER
C      OF FUNCTION EVALUATIONS TO BE APPROXIMATELY MAXNFE.
C      AS SET,THIS CORRESPONDS TO ABOUT 500 STEPS.
      DATA MAXNFE /9000/
C *******************************************************************
C     U - THE COMPUTER UNIT ROUNDOFF ERROR U IS THE SMALLEST POSITIVE
C         VALUE REPRESENTABLE IN THE MACHINE SUCH THAT  1.+ U .GT. 1.
C     (VARIABLE ONE SET TO 1.0 EASES PRECISION CONVERSION.)
C
      ONE = 1.0
      U = EPSLON(ONE)
C *******************************************************************
C      CHECK INPUT PARAMETERS
      IF (NEQN.LT.1) GO TO 10
      IF ((RELERR.LT.0.) .OR. (ABSERR.LT.0.)) GO TO 10
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
      IF ((KFLAG.EQ.4) .AND. (ABSERR.EQ.0.)) GO TO 40
      IF ((KFLAG.EQ.5) .AND. (RELERR.LE.SAVRE) .AND.
     * (ABSERR.LE.SAVAE)) GO TO 40
      GO TO 70
C IFLAG = 3,4,5,6, OR 7
   30 IF (IFLAG.EQ.3) GO TO 50
      IF ((IFLAG.EQ.4) .AND. (ABSERR.GT.0.)) GO TO 60
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
      RER = MAX(RELERR,32.*U+REMIN)
      U26 = 26.*U
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
      CALL F(A, Y, YP)
      NFE = 1
      IF (T.NE.TOUT) GO TO 90
      IFLAG = 2
      RETURN
   90 INIT = 1
      H = ABS(DT)
      TOLN = 0.
      DO 100 K=1,NEQN
        YG(K) = Y(K)
        YGP(K) = YP(K)
        TOL = RER*ABS(Y(K)) + ABSERR
        IF (TOL.LE.0.) GO TO 100
        TOLN = TOL
        YPK = ABS(YP(K))
        IF (YPK*H**5.GT.TOL) H = (TOL/YPK)**0.2
  100 CONTINUE
      IF (TOLN.LE.0.) H = 0.
      H = MAX(H,U26*MAX(ABS(T),ABS(DT)))
C *******************************************************************
C      SET STEPSIZE FOR INTEGRATION IN THE DIRECTION FROM T TO TOUT
  110 H = SIGN(H,DT)
C TEST TO SEE IF GERK IS BEING SEVERELY IMPACTED BY TOO MANY
C OUTPUT POINTS
      IF (ABS(H).GT.2.*ABS(DT)) KOP = KOP + 1
      IF (KOP.NE.100) GO TO 120
      KOP = 0
      IFLAG = 6
      RETURN
  120 IF (ABS(DT).GT.U26*ABS(T)) GO TO 140
C IF TOO CLOSE TO OUTPUT POINT,EXTRAPOLATE AND RETURN
      DO 130 K=1,NEQN
        YG(K) = YG(K) + DT*YGP(K)
        Y(K) = Y(K) + DT*YP(K)
  130 CONTINUE
      A = TOUT
      CALL F(A, YG, YGP)
      CALL F(A, Y, YP)
      NFE = NFE + 2
      GO TO 230
C INITIALIZE OUTPUT POINT INDICATOR
  140 OUTPUT = .FALSE.
C TO AVOID PREMATURE UNDERFLOW IN THE ERROR TOLERANCE FUNCTION,
C SCALE THE ERROR TOLERANCES
      SCALE = 2./RER
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
      IF (ABS(DT).GE.2.*ABS(H)) GO TO 170
      IF (ABS(DT).GT.ABS(H)) GO TO 160
C THE NEXT SUCCESSFUL STEP WILL COMPLETE THE INTEGRATION TO THE
C OUTPUT POINT
      OUTPUT = .TRUE.
      H = DT
      GO TO 170
  160 H = 0.5*DT
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
  180 CALL FEHL(F, NEQN, YG, T, H, YGP, F1, F2, F3, F4, F5, F1)
      NFE = NFE + 5
C COMPUTE AND TEST ALLOWABLE TOLERANCES VERSUS LOCAL ERROR
C ESTIMATES AND REMOVE SCALING OF TOLERANCES. NOTE THAT RELATIVE
C ERROR IS MEASURED WITH RESPECT TO THE AVERAGE MAGNITUDES OF THE
C OF THE SOLUTION AT THE BEGINNING AND END OF THE STEP.
      EEOET = 0.
      DO 200 K=1,NEQN
        ET = ABS(YG(K)) + ABS(F1(K)) + AE
        IF (ET.GT.0.) GO TO 190
C INAPPROPRIATE ERROR TOLERANCE
        IFLAG = 4
        KFLAG = 4
        RETURN
  190   EE = ABS((-2090.*YGP(K)+(21970.*F3(K)-15048.*F4(K)))
     *   +(22528.*F2(K)-27360.*F5(K)))
        EEOET = MAX(EEOET,EE/ET)
  200 CONTINUE
      ESTTOL = ABS(H)*EEOET*SCALE/752400.
      IF (ESTTOL.LE.1.) GO TO 210
C UNSUCCESSFUL STEP
C                   REDUCE THE STEPSIZE , TRY AGAIN
C                   THE DECREASE IS LIMITED TO A FACTOR OF 1/10
      HFAILD = .TRUE.
      OUTPUT = .FALSE.
      S = 0.1
      IF (ESTTOL.LT.59049.) S = 0.9/ESTTOL**0.2
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
      DO 220 K=1,NEQN
        YG(K) = F1(K)
  220 CONTINUE
      A = T
      CALL F(A, YG, YGP)
      NFE = NFE + 1
C NOW ADVANCE THE Y SOLUTION OVER TWO STEPS OF
C LENGTH H/2 AND EVALUATE DERIVATIVES THERE
      HH = 0.5*H
      CALL FEHL(F, NEQN, Y, TS, HH, YP, F1, F2, F3, F4, F5, Y)
      TS = TS + HH
      A = TS
      CALL F(A, Y, YP)
      CALL FEHL(F, NEQN, Y, TS, HH, YP, F1, F2, F3, F4, F5, Y)
      A = T
      CALL F(A, Y, YP)
      NFE = NFE + 12
C CHOOSE NEXT STEPSIZE
C THE INCREASE IS LIMITED TO A FACTOR OF 5
C IF STEP FAILURE HAS JUST OCCURRED, NEXT
C    STEPSIZE IS NOT ALLOWED TO INCREASE
      S = 5.
      IF (ESTTOL.GT.1.889568E-4) S = 0.9/ESTTOL**0.2
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
  240 DO 250 K=1,NEQN
        GERROR(K) = (YG(K)-Y(K))/31.
  250 CONTINUE
      RETURN
      END
      SUBROUTINE FEHL(F, NEQN, Y, T, H, YP, F1, F2, F3, F4, F5, S)
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
C     .. SCALAR ARGUMENTS ..
      INTEGER NEQN
      DOUBLE PRECISION H,T
C     ..
C     .. ARRAY ARGUMENTS ..
      DOUBLE PRECISION F1(NEQN),F2(NEQN),F3(NEQN),F4(NEQN),F5(NEQN),
     1     S(NEQN),Y(NEQN),YP(NEQN)
C     ..
C     .. SUBROUTINE ARGUMENTS ..
      EXTERNAL F
C     ..
C     .. LOCAL SCALARS ..
      INTEGER K
      DOUBLE PRECISION CH
C     ..
      CH = 0.25*H
      DO 10 K=1,NEQN
        F5(K) = Y(K) + CH*YP(K)
   10 CONTINUE
      CALL F(T+0.25*H, F5, F1)
      CH = 0.09375*H
      DO 20 K=1,NEQN
        F5(K) = Y(K) + CH*(YP(K)+3.*F1(K))
   20 CONTINUE
      CALL F(T+0.375*H, F5, F2)
      CH = H/2197.
      DO 30 K=1,NEQN
        F5(K) = Y(K) + CH*(1932.*YP(K)+(7296.*F2(K)-7200.*F1(K)))
   30 CONTINUE
      CALL F(T+12./13.*H, F5, F3)
      CH = H/4104.
      DO 40 K=1,NEQN
        F5(K) = Y(K) + CH*((8341.*YP(K)-845.*F3(K))+(29440.*F2(K)
     *   -32832.*F1(K)))
   40 CONTINUE
      CALL F(T+H, F5, F4)
      CH = H/20520.
      DO 50 K=1,NEQN
        F1(K) = Y(K) + CH*((-6080.*YP(K)+(9295.*F3(K)-5643.*F4(K)))
     *   +(41040.*F1(K)-28352.*F2(K)))
   50 CONTINUE
      CALL F(T+0.5*H, F1, F5)
C COMPUTE APPROXIMATE SOLUTION AT T+H
      CH = H/7618050.
      DO 60 K=1,NEQN
        S(K) = Y(K) + CH*((902880.*YP(K)+(3855735.*F3(K)-1371249.*
     *   F4(K)))+(3953664.*F2(K)+277020.*F5(K)))
   60 CONTINUE
      RETURN
      END
      DOUBLE PRECISION FUNCTION EPSLON (X)
      DOUBLE PRECISION X
C
C     ESTIMATE UNIT ROUNDOFF IN QUANTITIES OF SIZE X.
C
      DOUBLE PRECISION A,B,C,EPS,FOUR,THREE
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
      FOUR = 4.0
      THREE = 3.0
      A = FOUR/THREE
   10 B = A - 1.0
      C = B + B + B
      EPS = ABS(C-1.0)
      IF (EPS .EQ. 0.0) GO TO 10
      EPSLON = EPS*ABS(X)
      RETURN
      END
