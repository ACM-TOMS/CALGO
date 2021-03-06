      module sleignmd
      private
      public:: SLEIGN

      contains

      SUBROUTINE SLEIGN(A,B,INTAB,P0ATA,QFATA,P0ATB,QFATB,A1,A2,B1,B2,
     1                  NUMEIG,EIG,TOL,IFLAG,ISLFUN,SLFUN)
      INTEGER INTAB,NUMEIG,IFLAG,ISLFUN
      DOUBLE PRECISION A,B,P0ATA,QFATA,P0ATB,QFATB,A1,A2,B1,B2,EIG,TOL
      DOUBLE PRECISION SLFUN(*)
C     **********
C
C     This subroutine is designed for the calculation of a specified
C     eigenvalue, EIG, of a Sturm-Liouville problem in the form
C
C        (p(x)*y'(x))' + (q(x) + eig*r(x))*y(x) = 0  on  (a,b)
C
C     for user-supplied coefficient functions P, Q, and R.
C     The problem may be either regular or singular.  In the
C     regular case, boundary conditions of the form
C
C        a1*y(a) + a2*p(a)*y'(a) = 0
C        b1*y(b) + b2*p(b)*y'(b) = 0
C
C     are prescribed by specifying the numbers A1, A2, B1, and B2.
C     The index of the desired eigenvalue is specified in NUMEIG
C     and its requested accuracy in TOL.  Initial data for the
C     associated eigenfunction are also computed along with values
C     at selected points, if desired, in array SLFUN.
C
C     The SUBROUTINE statement is
C
C       SUBROUTINE sleign(a,b,intab,p0ata,qfata,p0atb,qfatb,a1,a2,b1,b2,
C                         numeig,eig,tol,iflag,islfun,slfun)
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
C         1.0 or -1.0 as the following properties of P, Q, and R at
C         the interval endpoints are true or false, respectively.
C
C         P0ATA -  p(a) is zero.              (If true, A is singular.)
C         QFATA -  q(a) and r(a) are finite.  (If false, A is singular.)
C         P0ATB -  p(b) is zero.              (If true, B is singular.)
C         QFATB -  q(b) and r(b) are finite.  (If false, B is singular.)
C
C       A1 and A2 are input variables set to prescribe the boundary
C         condition at A.  Their values are ignored if A is singular.
C
C       B1 and B2 are input variables set to prescribe the boundary
C         condition at B.  Their values are ignored if B is singular.
C
C       NUMEIG is an integer variable.  On input, it should be set to
C         the index of the desired eigenvalue (increasing sequence).
C         On output, it is unchanged unless the problem (apparently)
C         lacks NUMEIG eigenvalues, in which case it is reset to the
C         index of the largest eigenvalue.
C
C       EIG is a variable set on input to 0.0 or to an initial guess of
C         the eigenvalue.  If EIG is set to 0.0, SLEIGN will generate
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
C       IFLAG is an integer output variable set as follows.
C
C         IFLAG = 1 - successful problem solution.
C         IFLAG = 2 - improper input parameters.
C         IFLAG = 3 - NUMEIG exceeds actual number of eigenvalues.
C         IFLAG = 4 - some uncertainty about accuracy estimate TOL.
C         IFLAG = 5 - convergence too slow, best results returned.
C         IFLAG = 6 - failure, integrator could not complete.
C
C       ISLFUN is an integer input variable set to the number of
C         selected eigenfunction values desired.  If no values are
C         desired, set ISLFUN zero or negative.
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
C         continuous in the interval (XAA,XBB) and has weighted (by R)
C         L2-norm of 1.0 on the interval.  If ISLFUN is positive, then
C         on input the further ISLFUN locations of SLFUN specify the
C         points, in ascending order, where the eigenfunction values
C         are desired and on output contain the values themselves.
C
C     Subprograms called
C
C       user-supplied ..... p,q,r
C
C       sleign-supplied ... alfbet,dxdt,epslon,esteig,estpac,f,gerk
C
C     This version dated 5/18/89.
C     Paul B. Bailey, Sandia Laboratories, Albuquerque
C     Burton S. Garbow, Argonne National Laboratory
C
C     **********
C     .. Scalars in Common ..
      INTEGER INTSAV
      LOGICAL EIGF
      DOUBLE PRECISION ASAV,BSAV,C1,C2,EIGSAV,Z
C     ..
C     .. Local Scalars ..
      INTEGER I,IA,IB,IMAX,IMIN,IOUT,IP,J,JFLAG,KFLAG,LFLAG,MF,ML,
     1        NEIG,NMID
      LOGICAL AOK,BOK,BRACKT,CHNGAB,CHNGM,CONVRG,FYNYT,FYNYT1,
     1        LIMIT,LOGIC,NEWTON,ONEDIG,PRIN,THEGT0,THELT0
      DOUBLE PRECISION AA,AAA,ALFA,ASL,ASR,BB,BBB,BETA,
     1     C,CHNG,CL,CR,DAV,DE,DEDW,DEN,DERIVL,DERIVR,DIST,
     2     DPSIL,DPSIPL,DPSIPR,DPSIR,DT,DTHDA,DTHDAA,DTHDB,
     3     DTHDBB,DTHDE,DTHDEA,DTHDEB,DTHETA,DTHOLD,E,EEE,
     4     EIGLO,EIGLT,EIGRT,EIGUP,EL,ELIM,EMAX,EMIN,EOLD,EPS,EPSMIN,
     5     ER1,ER2,ESTERR,FLO,FMAX,FUP,GMAX,GUESS,H,ONE,PI,PIN,
     6     PSIL,PSIPL,PSIPR,PSIR,PX,QAV,QX,RATIO,RAV,RAY,RX,
     7     SL,SQL,SQR,SR,T,T1,T2,T3,TAU,TEMP,THRESH,TMAX,TMID,TMP,
     8     U,UL,UR,V,WL,X,X50,XAA,XBB,XBC,XMID,XSAV,ZAV
C     ..
C     .. Local Arrays ..
      INTEGER IWORK(5)
      DOUBLE PRECISION DS(98),ERL(3),ERR(3),PS(99),QS(99),RS(99),
     1     WORK(27),YL(3),YR(3)
C     ..
C     .. External Functions ..
      DOUBLE PRECISION P,Q,R
      EXTERNAL P,Q,R
C     ..
C     .. External Subroutines ..
cc      EXTERNAL ALFBET,DXDT,ESTEIG,ESTPAC,F,GERK
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS,ATAN,COS,EXP,INT,LOG,MAX,MIN,SIGN,SIN,TAN
C     ..
C     .. Common blocks ..
      COMMON /DATADT/ASAV,BSAV,C1,C2,INTSAV
      COMMON /DATAF/EIGSAV,EIGF
      COMMON /ZEE/Z
C     ..
C     Set constants EPSMIN, the computer unit roundoff error, and PI.
C     (Variable ONE set to 1.0 eases precision conversion.)
C
      ONE = 1.0
      EPSMIN = EPSLON(ONE)
      PI = 4.0*ATAN(ONE)
C
C     Set output device number.
C
      IOUT = 8
C
C     Check input parameters for errors.  If errors, return IFLAG=2.
C
      LOGIC = TOL.NE.0.0 .AND. 1.LE.INTAB .AND. INTAB.LE.4 .AND.
     1        P0ATA*QFATA*P0ATB*QFATB*NUMEIG.NE.0.0
      IF (INTAB.EQ.1) LOGIC = LOGIC .AND. A.LT.B
      IF (.NOT.LOGIC) THEN
         IFLAG = 2
         RETURN
         END IF
C
C     Set PRIN = .true. to trigger trace printout of successive steps.
C
      PRIN = .FALSE.
      IF (TOL.LT.0.0) PRIN = .TRUE.
C
C     Set EPS to the (initial) integration accuracy.
C
      EPS = 0.001
C
C     AOK (BOK) signals, if true, that endpoint A (B) is not singular.
C
      AOK = INTAB.LT.3 .AND. P0ATA.LT.0.0 .AND. QFATA.GT.0.0
      BOK = (INTAB.EQ.1 .OR. INTAB.EQ.3) .AND.
     1       P0ATB.LT.0.0 .AND. QFATB.GT.0.0
      NEIG = ABS(NUMEIG) - 1
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
C        END (SAVE-INPUT-DATA)
C
C     Evaluate P, Q, R to obtain preliminary information about the
C     differential equation.
C
C     DO (SAMPLE-COEFFICIENTS)
         THRESH = 1.0E+17
   10    CONTINUE
            CALL DXDT(EPSMIN,TEMP,X50)
            PX = P(X50)
            QX = Q(X50)
            RX = R(X50)
            PS(50) = PX
            QS(50) = QX/PX
            RS(50) = RX/PX
C
C     DAV,QAV,RAV are used in special case estimation when NUMEIG = 1,2.
C     EMIN = min(-Q/R), achieved at X for index value IMIN.
C     EMAX = max(-Q/R), achieved at X for index value IMAX.
C     MF and ML are the least and greatest index values, respectively.
C
            DAV = 0.0
            QAV = 0.0
            RAV = 0.0
            XSAV = X50
            EMIN = THRESH
            EMAX = -THRESH
            IF (RX.NE.0.0) THEN
               EMIN = -QX/RX
               EMAX = EMIN
               IMIN = 50
               IMAX = 50
               END IF
            H = 0.9/40.0
            DO 20 I=49,1,-1
               IF (I.GE.10) T = H*(I-50)
               IF (I.LT.10) T = T - 0.75*(1.0+T)
               CALL DXDT(T,TEMP,X)
               PX = P(X)
               QX = Q(X)
               RX = R(X)
               PS(I) = PX
               QS(I) = QX/PX
               RS(I) = RX/PX
               DS(I) = XSAV - X
               DAV = DAV + DS(I)
               QAV = QAV + DS(I)*(0.5*(QS(I)+QS(I+1))-QAV)/DAV
               RAV = RAV + DS(I)*(0.5*(RS(I)+RS(I+1))-RAV)/DAV
               XSAV = X
C
C     Try to avoid overflow by stopping when functions are large near A.
C
               FYNYT = (ABS(RX)+ABS(QX)+1.0/PX).LE.THRESH
               IF (RX.NE.0.0) THEN
                  IF (-QX/RX.LT.EMIN) THEN
                     EMIN = -QX/RX
                     IMIN = I
                     END IF
                  IF (-QX/RX.GT.EMAX) THEN
                     EMAX = -QX/RX
                     IMAX = I
                     END IF
                  END IF
               MF = I
               IF (.NOT.FYNYT) GO TO 30
   20          CONTINUE
   30       CONTINUE
            AAA = T
            IF (AOK) AAA = -1.0
            XSAV = X50
            DO 40 I=51,99
               IF (I.LE.90) T = H*(I-50)
               IF (I.GT.90) T = T + 0.75*(1.0-T)
               CALL DXDT(T,TEMP,X)
               PX = P(X)
               QX = Q(X)
               RX = R(X)
               PS(I) = PX
               QS(I) = QX/PX
               RS(I) = RX/PX
               DS(I-1) = X - XSAV
               DAV = DAV + DS(I-1)
               QAV = QAV + DS(I-1)*(0.5*(QS(I-1)+QS(I))-QAV)/DAV
               RAV = RAV + DS(I-1)*(0.5*(RS(I-1)+RS(I))-RAV)/DAV
               XSAV = X
C
C     Try to avoid overflow by stopping when functions are large near B.
C
               FYNYT1 = (ABS(QX)+ABS(RX)+1.0/PX).LE.THRESH
               IF (RX.NE.0.0) THEN
                  IF (-QX/RX.LT.EMIN) THEN
                     EMIN = -QX/RX
                     IMIN = I
                     END IF
                  IF (-QX/RX.GT.EMAX) THEN
                     EMAX = -QX/RX
                     IMAX = I
                     END IF
                  END IF
               ML = I - 1
               IF (.NOT.FYNYT1) GO TO 50
   40          CONTINUE
   50       CONTINUE
            BBB = T
            IF (BOK) BBB = 1.0
            LOGIC = C1.EQ.1.0 .AND. (.NOT.FYNYT .OR. .NOT.FYNYT1)
C
C     Modify (T,X) transformation corresponding to truncated interval.
C
            IF (LOGIC) THEN
               C1 = 0.5*(BBB-AAA)
               C2 = 0.5*(AAA+BBB)
               GO TO 10
               END IF
C
C     Estimate upper bound ELIM for EIG such that boundary conditions
C     can be satisfied.
C
         ELIM = EMAX + 1.0
         IF (INTAB.EQ.3 .OR. (P0ATA.GT.0.0 .AND. QFATA.LT.0.0)) THEN
            IF (-QS(MF)/RS(MF).LE.ELIM) THEN
               ELIM = -QS(MF)/RS(MF)
               IMAX = MF
               END IF
            END IF
         IF (INTAB.EQ.2 .OR. (P0ATB.GT.0.0 .AND. QFATB.LT.0.0)) THEN
            IF (-QS(ML)/RS(ML).LE.ELIM) THEN
               ELIM = -QS(ML)/RS(ML)
               IMAX = ML
               END IF
            END IF
         IF (INTAB.EQ.4) THEN
            ELIM = MIN(ELIM,-QS(MF)/RS(MF),-QS(ML)/RS(ML))
            IF (-QS(MF)/RS(MF).EQ.ELIM) IMAX = MF
            IF (-QS(ML)/RS(ML).EQ.ELIM) IMAX = ML
            END IF
         ELIM = ELIM - EPSMIN
         IF (ELIM.EQ.0.0) ELIM = -EPSMIN
         LIMIT = ELIM.LE.EMAX
C        END (SAMPLE-COEFFICIENTS)
      PIN = (NEIG+1)*PI
      IF (EIG.EQ.0.0) THEN
C        DO (ESTIMATE-EIG)
            CALL ESTEIG(MF,ML,LIMIT,ELIM,EMAX,EMIN,PIN,QS,RS,DS,PS,
     1                  IMAX,IMIN,EEE,EIG,IA,IB,EL,WL,DEDW)
C           END (ESTIMATE-EIG)
         END IF
      GUESS = EIG
C     DO (SET-INITIAL-INTERVAL-AND-MATCHPOINT)
         IF (GUESS.NE.0.0) THEN
C
C     Reduce overly large guess for EIG to upper bound if necessary.
C
            IF (LIMIT .AND. EIG.GT.ELIM) EIG = ELIM
            EEE = EIG
C           DO (ESTIMATE-PHASE-ANGLE-CHANGE)
               CALL ESTPAC(MF,ML,EEE,PIN,QS,RS,DS,PS,
     1                     IA,IB,IP,TEMP,U,TMP)
C              END (ESTIMATE-PHASE-ANGLE-CHANGE)
            END IF
C
C     Choose initial interval as large as possible that avoids overflow.
C     IA and IB are boundary indices for nonnegativity of EIG*R+Q.
C
         AA = -1.0
         IF (.NOT.AOK) THEN
            AA = H*(IA-50)
            IF (IA.LT.10) AA = -1.0 + 0.1/4.0**(10-IA)
            END IF
         BB = 1.0
         IF (.NOT.BOK) THEN
            BB = H*(IB-50)
            IF (IB.GT.90) BB = 1.0 - 0.1/4.0**(IB-90)
            END IF
         AA = AA + 0.6*(AAA-AA)
         BB = BB + 0.6*(BBB-BB)
C
C     Determine boundary values ALFA and BETA for theta at A and B.
C
         Z = 1.0
         CALL ALFBET(A,INTAB,AA,A1,A2,EEE,P0ATA,QFATA,AOK,
     1               ALFA,KFLAG,DERIVL)
         CALL ALFBET(B,INTAB,BB,B1,B2,EEE,P0ATB,QFATB,BOK,
     1               BETA,JFLAG,DERIVR)
         IF (.NOT.BOK) BETA = PI - BETA
C
C     Take boundary conditions into account in estimation of EIG.
C
         PIN = PIN + BETA - ALFA - PI
         IF (GUESS.EQ.0.0) EEE = EL + DEDW*(PIN-WL)
C
C     Subroutine ESTPAC must be called again because PIN has changed.
C
C        DO (ESTIMATE-PHASE-ANGLE-CHANGE)
            CALL ESTPAC(MF,ML,EEE,PIN,QS,RS,DS,PS,IA,IB,IP,TEMP,U,ZAV)
C           END (ESTIMATE-PHASE-ANGLE-CHANGE)
C
C     Choose the constant Z.
C
         IF (U.GT.0.0) Z = ZAV/U
C
C     Reset boundary values ALFA and BETA.
C
         CALL ALFBET(A,INTAB,AA,A1,A2,EEE,P0ATA,QFATA,AOK,
     1               ALFA,KFLAG,DERIVL)
         CALL ALFBET(B,INTAB,BB,B1,B2,EEE,P0ATB,QFATB,BOK,
     1               BETA,JFLAG,DERIVR)
         IF (.NOT.BOK) BETA = PI - BETA
         IF (PRIN) WRITE(IOUT,'(A,E22.14,A,E22.14)')
     1             ' alfa=',ALFA,'   beta=',BETA
C
C     Special case formula for estimation of EIG when NUMEIG = 1,2.
C
         IF (U.EQ.0.0 .AND. NEIG.LE.1 .AND. (BETA+NEIG*PI).LT.ALFA) THEN
            XBC = MAX(-1.0/TAN(ALFA),1.0/TAN(BETA))
            EEE = -(XBC*XBC-QAV)/RAV
            DEDW = XBC*(1.0+XBC*XBC)/RAV
            END IF
C
C     Choose initial matching point TMID.
C
         TMID = H*(IP-50)
         IF (TMID.LT.-0.8) TMID = -0.4
         IF (TMID.GT.0.8) TMID = 0.4
         IF (PRIN) WRITE(IOUT,'(A,E15.7,A,F11.8,A,E15.7)')
     1             ' estim=',EEE,'  tmid=',TMID,'  z=',Z
         IF (PRIN) WRITE(IOUT,'(A,F11.8,A,F11.8,A,F11.8,A,F11.8)')
     1             ' aaa=',AAA,'  aa=',AA,'  bb=',BB,'  bbb=',BBB
C        END (SET-INITIAL-INTERVAL-AND-MATCHPOINT)
C
C     End preliminary work, begin main task of computing EIG.
C
C     Logical variables have the following meanings if true.
C        AOK    - endpoint A is not singular.
C        BOK    - endpoint B is not singular.
C        CHNGM  - matching point TMID should be changed.
C        CHNGAB - one or both endpoints should be moved farther out.
C        BRACKT - EIG has been bracketed.
C        CONVRG - convergence test for EIG has been successfully passed.
C        NEWTON - Newton iteration may be employed.
C        THELT0 - lower bound for EIG has been found.
C        THEGT0 - upper bound for EIG has been found.
C        LIMIT  - upper bound exists with boundary conditions satisfied.
C        ONEDIG - most significant digit can be expected to be correct.
C        EIGF   - derivative argument is in original coordinate system.
C
      EIG = EEE
      EIGF = .FALSE.
      CHNGM = .FALSE.
      CHNGAB = .TRUE.
   60 CONTINUE
         IF (CHNGAB) THEN
C           DO (INITIAL-IZE)
               CHNGAB = .FALSE.
               BRACKT = .FALSE.
               CONVRG = .FALSE.
               THELT0 = .FALSE.
               THEGT0 = .FALSE.
               EIGLO = 0.0
               EIGLT = 0.0
               EIGRT = 0.0
               EIGUP = 0.0
               DTHOLD = 0.0
C              END (INITIAL-IZE)
            END IF
C
C     Recompute boundary conditions at singular endpoint(s).
C
C        DO (RESET-BOUNDARY-CONDITIONS)
            DERIVL = 0.0
            IF (.NOT.AOK) CALL ALFBET(A,INTAB,AA,A1,A2,EIG,
     1          P0ATA,QFATA,.FALSE.,ALFA,KFLAG,DERIVL)
            DERIVR = 0.0
            IF (.NOT.BOK) THEN
               CALL ALFBET(B,INTAB,BB,B1,B2,EIG,P0ATB,QFATB,.FALSE.,
     1                     BETA,JFLAG,DERIVR)
               BETA = PI - BETA
               END IF
            IF (PRIN) WRITE(IOUT,'(A,E22.14,A,E22.14)')
     1                ' alfa=',ALFA,'   beta=',BETA
C           END (RESET-BOUNDARY-CONDITIONS)
   70    CONTINUE
            IF (EIG.NE.GUESS .AND. .NOT.BRACKT) THEN
C
C     If initial guess was supplied, check that boundary conditions
C     can be satisfied at singular endpoints.  If not, try for
C     slightly lower EIG consistent with boundary conditions.
C
   80          CONTINUE
                  IF (.NOT.AOK) CALL ALFBET(A,INTAB,AA,A1,A2,EIG,
     1                P0ATA,QFATA,.FALSE.,TMP,KFLAG,TEMP)
                  IF (.NOT.BOK) CALL ALFBET(B,INTAB,BB,B1,B2,EIG,
     1                P0ATB,QFATB,.FALSE.,TMP,JFLAG,TEMP)
                  IF (KFLAG*JFLAG.NE.1) THEN
                     IF (THEGT0) EIG = 0.6*EIG + 0.4*EIGUP
                     IF (THELT0) EIG = 0.6*EIG + 0.4*EIGLO
                     IF (THELT0 .AND. EIGLO.LT.ELIM) EIGUP = ELIM
                     GO TO 80
                     END IF
               END IF
C           DO (OBTAIN-DTHETA-WITH-ONE-CORRECT-DIGIT)
   90          CONTINUE
                  IF (PRIN) WRITE(IOUT,'(/A,E22.14,A,E10.3,A,E10.3)')
     1                      ' guess=',EIG,'  eps=',EPS,'  tmid=',TMID
C                 DO (INTEGRATE-FOR-DTHETA)
C                    DO (SET-INITIAL-CONDITIONS)
                        DTHDEA = DERIVL
                        DTHDAA = 0.0
                        IF (.NOT.AOK) THEN
                           CALL DXDT(AA,DT,X)
                           PX = P(X)/Z
                           QX = Q(X)/Z
                           RX = R(X)/Z
                           C = EIG*RX + QX
                           DTHDAA = -(COS(ALFA)**2/PX +
     1                                C*SIN(ALFA)**2)*DT
C
C     Two special cases for DTHDAA.
C
                           IF (C.GE.0.0 .AND. P0ATA.LT.0.0 .AND.
     1                         QFATA.LT.0.0) DTHDAA = DTHDAA +
     1                         ALFA*DT/(X-A)
                           IF (C.GE.0.0 .AND. P0ATA.GT.0.0 .AND.
     1                         QFATA.GT.0.0) DTHDAA = DTHDAA +
     2                         (ALFA-0.5*PI)*DT/(X-A)
                           END IF
                        DTHDEB = -DERIVR
                        DTHDBB = 0.0
                        IF (.NOT.BOK) THEN
                           CALL DXDT(BB,DT,X)
                           PX = P(X)/Z
                           QX = Q(X)/Z
                           RX = R(X)/Z
                           C = EIG*RX + QX
                           DTHDBB = -(COS(BETA)**2/PX +
     1                                C*SIN(BETA)**2)*DT
C
C     Two special cases for DTHDBB.
C
                           IF (C.GE.0.0 .AND. P0ATB.LT.0.0 .AND.
     1                         QFATB.LT.0.0) DTHDBB = DTHDBB +
     2                         (PI-BETA)*DT/(B-X)
                           IF (C.GE.0.0 .AND. P0ATB.GT.0.0 .AND.
     1                         QFATB.GT.0.0) DTHDBB = DTHDBB +
     2                         (0.5*PI-BETA)*DT/(B-X)
                           END IF
                        TMAX = TMID
                        GMAX = ABS(DTHDEA)
                        EIGSAV = EIG
C                       END (SET-INITIAL-CONDITIONS)
C                                             T
C     YL = (theta,d(theta)/d(eig),d(theta)/da)
C
                     T = AA
                     YL(1) = ALFA
                     YL(2) = DTHDEA
                     YL(3) = 0.0
C
C     Use integrator in one-step mode towards change to different TMID.
C
                     LFLAG = 1
                     IF (CHNGM) LFLAG = -1
  100                CONTINUE
                        CALL GERK(F,3,YL,T,TMID,EPS,EPS,LFLAG,ERL,
     1                            WORK,IWORK)
                        IF (LFLAG.EQ.-2 .AND. T.GT.-0.8 .AND.
     1                      ABS(YL(2)).GT.GMAX) THEN
                           TMAX = T
                           GMAX = ABS(YL(2))
                           END IF
        IF (LFLAG.EQ.3) PRINT*,
     1  'After 9000 function evaluations GERK reached T=',T
                        IF (LFLAG.EQ.3 .OR. LFLAG.EQ.-2) GO TO 100
                     IF (LFLAG.GT.3) THEN
                        IFLAG = 6
                        RETURN
                        END IF
                     DTHDA = DTHDAA*EXP(-2.0*YL(3))
C                                             T
C     YR = (theta,d(theta)/d(eig),d(theta)/db)
C
                     T = BB
                     YR(1) = BETA + PI*NEIG
                     YR(2) = DTHDEB
                     YR(3) = 0.0
C
C     Use integrator in one-step mode towards change to different TMID.
C
                     LFLAG = 1
                     IF (CHNGM) LFLAG = -1
  110                CONTINUE
                        CALL GERK(F,3,YR,T,TMID,EPS,EPS,LFLAG,ERR,
     1                            WORK,IWORK)
                        IF (LFLAG.EQ.-2 .AND. T.LT.0.8 .AND.
     1                      ABS(YR(2)).GT.GMAX) THEN
                           TMAX = T
                           GMAX = ABS(YR(2))
                           END IF
        IF (LFLAG.EQ.3) PRINT*,
     1  'After 9000 function evaluations GERK reached T=',T
                        IF (LFLAG.EQ.3 .OR. LFLAG.EQ.-2) GO TO 110
                     IF (LFLAG.GT.3) THEN
                        IFLAG = 6
                        RETURN
                        END IF
                     DTHDB = DTHDBB*EXP(-2.0*YR(3))
C
C     DTHETA measures theta difference from left and right integrations.
C
                     DTHETA = YL(1) - YR(1)
                     DTHDE = YL(2) - YR(2)
                     ER1 = ERL(1) - ERR(1)
                     ER2 = ERL(2) - ERR(2)
                     TMID = TMAX
C                    END (INTEGRATE-FOR-DTHETA)
C
C     Define ONEDIG to try to be sure of one correct digit in DTHETA.
C     Redo integrations with tighter tolerance until ONEDIG is true.
C
                  ONEDIG = ABS(ER1).LE.0.1*ABS(DTHETA) .AND.
     1                     ABS(ER2).LE.0.1*ABS(DTHDE)
                  NEWTON = ABS(DTHETA).LT.0.06
                  IF (NEWTON) THEN
C                    DO (COMPUTE-CONVRG)
C
C     Measure convergence after adding separate contributions to error.
C
                        T1 = ABS(DTHETA)+50.0*ABS(ER1)
                        T2 = (1.0+AA)*ABS(DTHDA)
                        T3 = (1.0-BB)*ABS(DTHDB)
                        ESTERR = (T1+T2+T3)/ABS(DTHDE)/MAX(ONE,ABS(EIG))
                        CONVRG = ESTERR.LE.TAU
                        IF (PRIN) WRITE(IOUT,'(A,L2)')
     1                            ' converge=',CONVRG
                        IF (PRIN .AND. .NOT.CONVRG)
     1                            WRITE(IOUT,'(A,E15.7)')
     2                            ' estim. acc.=',ESTERR
C                       END (COMPUTE-CONVRG)
                     END IF
                  IF (.NOT.ONEDIG .OR.
     1                ABS(ER1).GT.0.01*ABS(DTHETA)) THEN
C
C     Reduce local error criterion, but return IFLAG=5 if too small.
C
                     EPS = 0.05*EPS
                     IF (EPS.LE.EPSMIN) THEN
                        IFLAG = 5
                        RETURN
                        END IF
                     END IF
                  IF (.NOT.(ONEDIG .OR. CONVRG)) GO TO 90
               IF (ABS(DTHETA).LT.0.1) CHNGM = .FALSE.
               IF (PRIN) WRITE(IOUT,'(A,E15.7,A,E15.7)')
     1                   ' dtheta=',DTHETA,'   dthde=',DTHDE
               IF (PRIN) WRITE(IOUT,'(/A,E15.7,A,E15.7)')
     1                   ' thetal=',YL(1),'   thetar=',YR(1)
C              END (OBTAIN-DTHETA-WITH-ONE-CORRECT-DIGIT)
C           DO (SET-BRACKET-AND-FORM-NEWTON-ITERATES)
C
C     EIG is bracketed when both THEGT0=.true. and THELT0=.true.
C
               IF (DTHETA*DTHDE.GT.0.0) THEN
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
               BRACKT = THELT0 .AND. THEGT0
               IF (PRIN) WRITE(IOUT,'(A,E22.14,A,E22.14)')
     1                   ' eigrt=',EIGRT,'  eigup=',EIGUP
               IF (PRIN) WRITE(IOUT,'(A,E22.14,A,E22.14)')
     1                   ' eiglt=',EIGLT,'  eiglo=',EIGLO
C              END (SET-BRACKET-AND-FORM-NEWTON-ITERATES)
            IF (NEWTON) THEN
               CHNGM = .TRUE.
               IF (CONVRG) THEN
                  CHNG = DTHDA*(-1.0-AA) - DTHDB*(1.0-BB)
                  TEMP = (DTHETA+CHNG)/DTHDE
                  EIG = EIG - TEMP
                  TOL = ABS(TEMP)/MAX(ONE,ABS(EIG))
               ELSE
                  CHNGAB = T1.LT.0.5*(T2+T3)
C
C     Move endpoint(s) out or take Newton step, according to CHNGAB.
C
                  IF (CHNGAB) THEN
C                    DO (MOVE-ENDPOINTS)
                        IF (T2.GT.T1 .AND. AA.GT.AAA)
     1                      AA = AA + 0.8*(-1.0-AA)
                        IF (T3.GE.T1 .AND. BB.LT.BBB)
     1                      BB = BB + 0.8*(1.0-BB)
                        AA = MAX(AA,AAA)
                        BB = MIN(BB,BBB)
                        IF ((AAA-AA).EQ.(BBB-BB)) THEN
C
C     Cannot move endpoint(s) again. Store estimates and return IFLAG=5.
C
                           CHNG = (DTHDA-DTHDB)*(AAA-AA)
                           TEMP = (DTHETA+CHNG)/DTHDE
                           EIG = EIG - TEMP
                           TOL = ABS(TEMP)/MAX(ONE,ABS(EIG))
                           IFLAG = 5
                           RETURN
                           END IF
                        EEE = EIG
                        IF (PRIN) WRITE(IOUT,'(A,2E15.8)')
     1                            ' new endpoints    ',AA,BB
C                       END (MOVE-ENDPOINTS)
                  ELSE
                     EIG = EIG - DTHETA/DTHDE
                     END IF
                  END IF
            ELSE
               IF (BRACKT) THEN
C
C     Obtain next estimate of EIG by bisection or linear interpolation.
C
                  FMAX = MAX(-FLO,FUP)
                  EOLD = EIG
                  RATIO = DTHETA/DTHOLD
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
                     CHNGAB = RATIO.GE.0.4 .AND. .NOT.(AOK.AND.BOK)
                     IF (ABS(DE).LT.EPSMIN) THEN
                        TOL = ABS(DE)/MAX(ONE,ABS(EIG))
                        IFLAG = 5
                        RETURN
                        END IF
                     END IF
                  CHNGM = .NOT.CHNGM .AND. RATIO.GE.0.4
               ELSE
C                 DO (TRY-FOR-BRACKET)
C
C     Take twice the Newton step in trying for a bracket.
C
                     IF (EIG.EQ.EEE) THEN
                        IF (GUESS.NE.0.0) DEDW = 1.0/DTHDE
                        CHNG = -(DEDW+1.0/DTHDE)*DTHETA
                        IF (ABS(CHNG).GT.0.1*ABS(EIG))
     1                      CHNG = -0.1*SIGN(EIG,DTHETA)
                     ELSE
                        CHNG = -2.0*DTHETA/DTHDE
                        END IF
                     LOGIC = EIG.NE.EEE .AND. DTHOLD.LT.0.0 .AND.
     1                       LIMIT .AND. CHNG.GT.(ELIM-EIG)
                     IF (LOGIC) THEN
                        CHNG = 0.99*(ELIM-EIG+EPSMIN)
                        IF (CHNG.LT.EPSMIN) THEN
C
C     If change is too small, EIG is presumed not to exist (IFLAG=3).
C
                           NUMEIG = NEIG - INT(-DTHETA/PI)
                           IFLAG = 3
                           RETURN
                           END IF
C
C     Limit change in EIG to a factor of 10.
C
                     ELSE IF (ABS(CHNG).GT.(1.0+10.0*ABS(EIG))) THEN
                        CHNG = SIGN(1.0+10.0*ABS(EIG),CHNG)
                     ELSE IF (ABS(EIG).GE.1.0 .AND.
     1                        ABS(CHNG).LT.0.1*ABS(EIG)) THEN
                        CHNG = 0.1*SIGN(EIG,CHNG)
                        END IF
                     EOLD = EIG
                     EIG = EIG + CHNG
C                    END (TRY-FOR-BRACKET)
                  END IF
               END IF
            DTHOLD = DTHETA
            IF (.NOT.(CONVRG .OR. CHNGAB .OR. NEWTON)) GO TO 70
         IF (.NOT.CONVRG) GO TO 60
      IFLAG = 1
      IF (PRIN) WRITE(IOUT,'(A,I7,A,E22.14,A,E10.3)')
     1          ' numeig=',NUMEIG,'  eig=',EIG,'  tol=',TOL
C     DO (COMPUTE-EIGENFUNCTION-DATA)
C
C     Convert from T to X values, fill 7 of first 9 locations of SLFUN.
C
         CALL DXDT(TMID,TEMP,XMID)
         CALL DXDT(AA,TEMP,XAA)
         CALL DXDT(BB,TEMP,XBB)
         SLFUN(1) = XMID
         SLFUN(2) = XAA
         SLFUN(3) = ALFA
         SLFUN(5) = XBB
         SLFUN(6) = BETA + PI*NEIG
         SLFUN(8) = EPS
         SLFUN(9) = Z
C
C     Compute SLFUN(4), SLFUN(7) towards normalizing the eigenfunction.
C
         EIGSAV = EIG
         Z = -Z
         T = AA
         YL(1) = ALFA
         YL(2) = DTHDEA
         YL(3) = 0.0
         LFLAG = 1
  120    CONTINUE
            CALL GERK(F,3,YL,T,TMID,EPS,EPS,LFLAG,ERL,WORK,IWORK)
        IF (LFLAG.EQ.3) PRINT*,
     1  'After 9000 function evaluations GERK reached T=',T
            IF (LFLAG.EQ.3) GO TO 120
         T = BB
         YR(1) = BETA + PI*NEIG
         YR(2) = DTHDEB
         YR(3) = 0.0
         LFLAG = 1
  130    CONTINUE
            CALL GERK(F,3,YR,T,TMID,EPS,EPS,LFLAG,ERR,WORK,IWORK)
        IF (LFLAG.EQ.3) PRINT*,
     1  'After 9000 function evaluations GERK reached T=',T
            IF (LFLAG.EQ.3) GO TO 130
         Z = -Z
         SL = SIN(YL(1))
         SR = SIN(YR(1))
         CL = COS(YL(1))
         CR = COS(YR(1))
         UL = (YL(2)-DTHDEA*EXP(-2.0*YL(3)))*Z
         UR = (YR(2)-DTHDEB*EXP(-2.0*YR(3)))*Z
         ASL = ABS(SL)
         ASR = ABS(SR)
         DEN = 0.5*LOG(UL*ASR*ASR-UR*ASL*ASL)
         SLFUN(4) = LOG(ASR) - YL(3) - DEN
         SLFUN(7) = LOG(ASL) - YR(3) - DEN
C        END (COMPUTE-EIGENFUNCTION-DATA)
C     DO (CHECK-MATCHING-VALUES-OF-EIGENFUNCTION)
C
C     Perform final check on EIG. Return IFLAG=4 if not accurate enough.
C
         E = ASR*EXP(-DEN)
         PSIL = E*SL
         PSIPL = E*CL
         SQL = E*E*UL
         DPSIL = PSIL*ERL(3) + PSIPL*ERL(1)
         DPSIPL = PSIL*ERL(1) + PSIPL*ERL(3)
         PSIPL = PSIPL*Z
         E = ASL*EXP(-DEN)
         PSIR = E*SR
         PSIPR = E*CR
         SQR = E*E*UR
         DPSIR = PSIR*ERR(3) + PSIPR*ERR(1)
         DPSIPR = PSIR*ERR(1) + PSIPR*ERR(3)
         PSIPR = PSIPR*Z
         RAY = EIG + (PSIL*PSIPL-PSIR*PSIPR)/(SQL-SQR)
         IF (PRIN) THEN
            WRITE(IOUT,'(A,E22.14)') ' ray=',RAY
            WRITE(IOUT,'(A,E22.14,A,E22.14)')
     1      ' psil=',PSIL,'  psir=',PSIR
            WRITE(IOUT,'(A,E22.14,A,E22.14)')
     1      ' psipl=',PSIPL,'  psipr=',PSIPR
            WRITE(IOUT,'(A,E22.14,A,E22.14)')
     1      ' sql=',SQL,'  sqr=',SQR
            WRITE(IOUT,'(A,E22.14,A,E22.14)')
     1      ' dpsil=',DPSIL,'  dpsir=',DPSIR
            WRITE(IOUT,'(A,E22.14,A,E22.14)')
     1      ' dpsipl=',DPSIPL,'  dpsipr=',DPSIPR
            END IF
C        END (CHECK-MATCHING-VALUES-OF-EIGENFUNCTION)
      IF (ABS(RAY-EIG).GT.TAU*MAX(ONE,ABS(EIG))) IFLAG = 4
      IF (ISLFUN.GT.0) THEN
C
C     Calculate selected eigenfunction values by integration.
C
C        DO (GENERATE-EIGENFUNCTION-VALUES)
            EIGF = .TRUE.
            NMID = 0
            DO 140 I=1,ISLFUN
               IF (SLFUN(9+I).LE.SLFUN(1)) NMID = I
  140          CONTINUE
            IF (NMID.GT.0) THEN
               X = SLFUN(2)
               YL(1) = SLFUN(3)
               YL(2) = 0.0
               YL(3) = SLFUN(4)
               LFLAG = 1
               DO 160 J=1,NMID
C Move any eigenfunction output points to lie within the range XAA..XBB
CJDP (an attempt to stop errors at singular endpoints)
                  SLFUN(J+9)=min(max(SLFUN(J+9),SLFUN(2)),SLFUN(5))
  150             CONTINUE
                     IF (X.NE.SLFUN(J+9))
     1                CALL GERK(F,3,YL,X,SLFUN(J+9),SLFUN(8),SLFUN(8),
     1                         LFLAG,ERL,WORK,IWORK)
        IF (LFLAG.EQ.3) PRINT*,
     1  'After 9000 function evaluations GERK reached X=',X
                     IF (LFLAG.EQ.3) GO TO 150
                  IF (LFLAG.EQ.6) LFLAG = 2
                  SLFUN(J+9) = EXP(YL(3))*SIN(YL(1))
  160             CONTINUE
               END IF
            IF (NMID.LT.ISLFUN) THEN
               X = SLFUN(5)
               YR(1) = SLFUN(6)
               YR(2) = 0.0
               YR(3) = SLFUN(7)
               LFLAG = 1
               DO 180 J=ISLFUN,NMID+1,-1
C Move any eigenfunction output points to lie within the range XAA..XBB
CJDP (see above)
                  SLFUN(J+9)=min(max(SLFUN(J+9),SLFUN(2)),SLFUN(5))
  170             CONTINUE
                     IF (X.NE.SLFUN(J+9))
     1                CALL GERK(F,3,YR,X,SLFUN(J+9),SLFUN(8),SLFUN(8),
     1                         LFLAG,ERR,WORK,IWORK)
        IF (LFLAG.EQ.3) PRINT*,
     1  'After 9000 function evaluations GERK reached X=',X
                     IF (LFLAG.EQ.3) GO TO 170
                  IF (LFLAG.EQ.6) LFLAG = 2
                  SLFUN(J+9) = EXP(YR(3))*SIN(YR(1))
  180             CONTINUE
               END IF
C           END (GENERATE-EIGENFUNCTION-VALUES)
         END IF
      RETURN
      END SUBROUTINE SLEIGN
C ----------------------------------------------------------------------
      SUBROUTINE ALFBET(XEND,INTAB,TT,COEF1,COEF2,EIG,P0,QF,OK,
     1                  VALUE,IFLAG,DERIV)
      INTEGER INTAB,IFLAG
      LOGICAL OK
      DOUBLE PRECISION XEND,TT,COEF1,COEF2,EIG,P0,QF,VALUE,DERIV
C     **********
C
C     This subroutine computes a boundary value for a specified endpoint
C     of the interval for a Sturm-Liouville problem in the form
C
C        (p(x)*y'(x))' + (q(x) + eig*r(x))*y(x) = 0  on  (a,b)
C
C     for user-supplied coefficient functions P, Q, and R.  It is called
C     from SLEIGN.  Both regular and singular endpoints are treated.
C
C     Subprograms called
C
C       user-supplied ..... p,q,r
C
C       sleign-supplied ... dxdt,extrap
C
C     **********
C     .. Scalars in Common ..
      DOUBLE PRECISION Z
C     ..
C     .. Local Scalars ..
      LOGICAL LOGIC
      DOUBLE PRECISION C,CD,D,HH,ONE,PI,PX,QX,RX,T,TEMP,X
C     ..
C     .. External Functions ..
      DOUBLE PRECISION P,Q,R
      EXTERNAL P,Q,R
C     ..
C     .. External Subroutines ..
cc      EXTERNAL DXDT,EXTRAP
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
      IF (OK) THEN
         VALUE = 0.5*PI
         IF (COEF1.NE.0.0) VALUE = ATAN(-Z*COEF2/COEF1)
         LOGIC = (TT.LT.0.0 .AND. VALUE.LT.0.0) .OR.
     1           (TT.GT.0.0 .AND. VALUE.LE.0.0)
         IF (LOGIC) VALUE = VALUE + PI
      ELSE
         LOGIC = (INTAB.EQ.2 .AND. TT.GT.0.0) .OR.
     1           (INTAB.EQ.3 .AND. TT.LT.0.0) .OR.
     2           INTAB.EQ.4 .OR. (P0.GT.0.0 .AND. QF.LT.0.0)
         IF (LOGIC) THEN
            T = SIGN(ONE,TT)
            CALL EXTRAP(T,TT,EIG,VALUE,DERIV,IFLAG)
         ELSE
            CALL DXDT(TT,TEMP,X)
            PX = P(X)/Z
            QX = Q(X)/Z
            RX = R(X)/Z
            C = 2.0*(EIG*RX+QX)
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
      END SUBROUTINE ALFBET
C ----------------------------------------------------------------------
      SUBROUTINE DXDT(T,DT,X)
      DOUBLE PRECISION T,DT,X
C     **********
C
C     This subroutine transforms coordinates from T on (-1,1) to
C     X on (A,B) in the solution of a Sturm-Liouville problem.
C     It is called from subroutines SLEIGN, ALFBET, F, and EXTRAP.
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
      END SUBROUTINE DXDT
C ----------------------------------------------------------------------
      SUBROUTINE ESTEIG(MF,ML,LIMIT,ELIM,EMAX,EMIN,PIN,QS,RS,DS,PS,
     1                  IMAX,IMIN,EEE,EIG,IA,IB,EL,WL,DEDW)
      INTEGER MF,ML,IMAX,IMIN,IA,IB
      LOGICAL LIMIT
      DOUBLE PRECISION ELIM,EMAX,EMIN,PIN,EEE,EIG,EL,WL,DEDW
      DOUBLE PRECISION QS(ML),RS(ML),DS(ML),PS(ML)
C     **********
C
C     This subroutine generates an initial guess for a specified
C     eigenvalue of a Sturm-Liouville problem in the form
C
C        (p(x)*y'(x))' + (q(x) + eig*r(x))*y(x) = 0  on  (a,b)
C
C     for user-supplied coefficient functions P, Q, and R.  It is
C     called from SLEIGN when no initial guess is provided by the user.
C
C     The method used is to approximately solve the equation for EIG
C
C        Integral (sqrt((eig*r+q)/p)) = numeig*pi
C
C     where the integral is taken over those X in (A,B) for which
C
C        (eig*r+q)/p .gt. 0
C
C     and NUMEIG is the index of the desired eigenvalue (PIN=NUMEIG*pi).
C
C     Subprograms called
C
C       sleign-supplied ... estpac
C
C     **********
C     .. Scalars in Common ..
      INTEGER INTAB
      DOUBLE PRECISION A,B,C1,C2
C     ..
C     .. Local Scalars ..
      INTEGER IE,IP
      LOGICAL LOGIC
      DOUBLE PRECISION BALLPK,EU,FNEW,FOLD,SUM,TEMP,U,WU
C     ..
C     .. External Subroutines ..
cc      EXTERNAL ESTPAC
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS,MIN
C     ..
C     .. Common blocks ..
      COMMON /DATADT/A,B,C1,C2,INTAB
C     ..
      EEE = MIN(ELIM,EMAX)
C     DO (ESTIMATE-PHASE-ANGLE-CHANGE)
         CALL ESTPAC(MF,ML,EEE,PIN,QS,RS,DS,PS,IA,IB,IP,SUM,U,TEMP)
C        END (ESTIMATE-PHASE-ANGLE-CHANGE)
C
C     Choose bounds for EIG and associate function (integral) values.
C
      EL = EMIN
      WL = 0.0
      EU = EEE
      WU = SUM
      IF (LIMIT .AND. WU.LT.PIN) THEN
         EIG = ELIM
      ELSE
         IF (U.EQ.0.0) THEN
            EL = EMAX
            EEE = EMAX + 1.0
C           DO (ESTIMATE-PHASE-ANGLE-CHANGE)
               CALL ESTPAC(MF,ML,EEE,PIN,QS,RS,DS,PS,
     1                     IA,IB,IP,SUM,U,TEMP)
C              END (ESTIMATE-PHASE-ANGLE-CHANGE)
            EU = EEE
            WU = SUM
            END IF
   10    CONTINUE
            IF (WU.LE.PIN) THEN
C
C     Increase trial value if integral is still too small.
C
               EL = EU
               WL = WU
               EEE = EU + ((PIN-WU+3.0)/U)**2
C              DO (ESTIMATE-PHASE-ANGLE-CHANGE)
                  CALL ESTPAC(MF,ML,EEE,PIN,QS,RS,DS,PS,
     1                        IA,IB,IP,SUM,U,TEMP)
C                 END (ESTIMATE-PHASE-ANGLE-CHANGE)
               EU = EEE
               WU = SUM
               GO TO 10
               END IF
C
C     EIG is bracketed.  Now try to reduce the size of the bracket
C     by searching among the saved values of -QS()/RS().
C
   20    CONTINUE
            IF (ABS(IMAX-IMIN).GE.2 .AND. EU.LE.EMAX) THEN
               IE = (IMAX+IMIN)/2
               IF (RS(IE).NE.0.0) THEN
                  EEE = -QS(IE)/RS(IE)
C                 DO (ESTIMATE-PHASE-ANGLE-CHANGE)
                     CALL ESTPAC(MF,ML,EEE,PIN,QS,RS,DS,PS,
     1                           IA,IB,IP,SUM,U,TEMP)
C                    END (ESTIMATE-PHASE-ANGLE-CHANGE)
                  IF (SUM.GT.PIN) THEN
                     IMAX = IE
                     WU = SUM
                     EU = EEE
                  ELSE
                     IMIN = IE
                     WL = SUM
                     EL = EEE
                     END IF
                  GO TO 20
                  END IF
               END IF
C
C     Improve approximation for EIG using bisection or secant method.
C     Substitute 'ballpark' estimate if approximation grows too large.
C
         DEDW = (EU-EL)/(WU-WL)
         FOLD = 0.0
         IF (INTAB.EQ.1) BALLPK = (PIN/(A-B))**2
         LOGIC = .TRUE.
   30    CONTINUE
            IF (LOGIC) THEN
               LOGIC = (WL.LT.(PIN-1.0) .OR. WU.GT.(PIN+1.0))
               EEE = EL + DEDW*(PIN-WL)
               FNEW = MIN(PIN-WL,WU-PIN)
               IF (FNEW.GT.0.4*FOLD .OR. FNEW.LE.1.0)
     1             EEE = 0.5*(EL+EU)
               IF (INTAB.EQ.1 .AND. ABS(EEE).GT.1.0E+3*BALLPK) THEN
                  EIG = BALLPK
                  RETURN
               ELSE IF (INTAB.NE.1 .AND. ABS(EEE).GT.1.0E+6) THEN
                  EIG = 1.0
                  RETURN
               ELSE
                  FOLD = FNEW
C                 DO (ESTIMATE-PHASE-ANGLE-CHANGE)
                     CALL ESTPAC(MF,ML,EEE,PIN,QS,RS,DS,PS,
     1                           IA,IB,IP,SUM,U,TEMP)
C                    END (ESTIMATE-PHASE-ANGLE-CHANGE)
                  IF (SUM.LT.PIN) THEN
                     EL = EEE
                     WL = SUM
                  ELSE
                     EU = EEE
                     WU = SUM
                     END IF
                  DEDW = (EU-EL)/(WU-WL)
                  GO TO 30
                  END IF
               END IF
         END IF
      RETURN
      END SUBROUTINE ESTEIG
C ----------------------------------------------------------------------
      SUBROUTINE ESTPAC(MF,ML,EEE,PIN,QS,RS,DS,PS,IA,IB,IP,SUM,U,ZAV)
      INTEGER MF,ML,IA,IB,IP
      DOUBLE PRECISION EEE,PIN,SUM,U,ZAV
      DOUBLE PRECISION QS(ML),RS(ML),DS(ML),PS(ML)
C     **********
C
C     This subroutine estimates the change in 'phase angle' in the
C     eigenvalue determination of a Sturm-Liouville problem in the form
C
C        (p(x)*y'(x))' + (q(x) + eig*r(x))*y(x) = 0  on  (a,b)
C
C     for user-supplied coefficient functions P, Q, and R.  It is
C     called from SLEIGN, and also from ESTEIG when no initial guess
C     is provided by the user.
C
C     The subroutine approximates the (trapezoidal rule) integral of
C
C        sqrt((eig*r+q)/p)
C
C     where the integral is taken over those X in (A,B) for which
C
C        (eig*r+q)/p .gt. 0
C
C     **********
C     .. Local Scalars ..
      INTEGER J
      DOUBLE PRECISION PSUM,RT,RTSAV,V,W,WSAV
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS,SIGN,SQRT
C     ..
      IA = MF
      IB = 80
      IP = MF
C
C     SUM accumulates the integral approximation.  U measures the total
C     length of subintervals where (EIG*R+Q)/P .gt. 0.0.  ZAV is the
C     average value of sqrt((EIG*R+Q)*P) over those subintervals.
C
      SUM = 0.0
      U = 0.0
      ZAV = 0.0
      WSAV = EEE*RS(MF) + QS(MF)
      IF (WSAV.GT.0.0) THEN
         RTSAV = SQRT(WSAV)
      ELSE
         RTSAV = 0.0
         END IF
      DO 10 J=MF+1,ML
         W = EEE*RS(J) + QS(J)
         IF (W.GT.0.0) THEN
            IF (J.GE.80) IB = J
            U = U + DS(J-1)
            RT = SQRT(W)
         ELSE
            RT = 0.0
            IF (U.EQ.0.0 .AND. RTSAV.EQ.0.0 .AND. IA.LE.19) IA = IA + 1
            END IF
         IF (W.EQ.0.0 .OR. WSAV.EQ.0.0 .OR. W.EQ.SIGN(W,WSAV)) THEN
            V = RT + RTSAV
         ELSE
            V = (W*RT+WSAV*RTSAV)/ABS(W-WSAV)
            END IF
         WSAV = W
         RTSAV = RT
         PSUM = DS(J-1)*V
         IF (PSUM.LT.(PIN-SUM)) IP = J
         SUM = SUM + PSUM
         IF (U.GT.0.0) ZAV = ZAV + PSUM*(PS(J)+PS(J-1))
   10    CONTINUE
      SUM = 0.5*SUM
      ZAV = 0.25*ZAV
      RETURN
      END SUBROUTINE ESTPAC
C ----------------------------------------------------------------------
      SUBROUTINE EXTRAP(T,TT,EIG,VALUE,DERIV,IFLAG)
      INTEGER IFLAG
      DOUBLE PRECISION T,TT,EIG,VALUE,DERIV
C     **********
C
C     This subroutine is called from ALFBET in determining boundary
C     values at a singular endpoint of the interval for a
C     Sturm-Liouville problem in the form
C
C        (p(x)*y'(x))' + (q(x) + eig*r(x))*y(x) = 0  on  (a,b)
C
C     for user-supplied coefficient functions P, Q, and R.
C
C     EXTRAP, which in turn calls INTPOL, extrapolates the function
C
C        arctan(1.0/sqrt(-p*(eig*r+q)))
C
C     from its values for T within (-1,1) to an endpoint.
C
C     Subprograms called
C
C       user-supplied ..... p,q,r
C
C       sleign-supplied ... dxdt,intpol
C
C     **********
C     .. Scalars in Common ..
      DOUBLE PRECISION Z
C     ..
C     .. Local Scalars ..
      INTEGER KGOOD
      DOUBLE PRECISION ANS,CTN,ERROR,PROD,PX,QX,RX,T1,TEMP,X
C     ..
C     .. Local Arrays ..
      DOUBLE PRECISION FN1(5),XN(5)
C     ..
C     .. External Functions ..
      DOUBLE PRECISION P,Q,R
      EXTERNAL P,Q,R
C     ..
C     .. External Subroutines ..
cc      EXTERNAL DXDT,INTPOL
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
         RX = R(X)/Z
         PROD = -PX*(EIG*RX+QX)
         IF (PROD.LE.0.0) THEN
            T1 = 0.5*(T1+T)
            IF ((1.0+(T1-T)**2).GT.1.0) GO TO 10
            IFLAG = 5
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
      DERIV = 0.5*PX*RX/CTN/(1.0+CTN**2)
      TT = XN(1)
      RETURN
      END SUBROUTINE EXTRAP
C ----------------------------------------------------------------------
      SUBROUTINE F(T,Y,YP)
      DOUBLE PRECISION T
      DOUBLE PRECISION Y(2),YP(3)
C     **********
C
C     This subroutine evaluates the derivative functions for use with
C     integrator GERK in solving a Sturm-Liouville problem in the form
C
C        (p(x)*y'(x))' + (q(x) + eig*r(x))*y(x) = 0  on  (a,b)
C
C     for user-supplied coefficient functions P, Q, and R.
C
C     Subprograms called
C
C       user-supplied ..... p,q,r
C
C       sleign-supplied ... dxdt
C
C     **********
C     .. Scalars in Common ..
      LOGICAL EIGF
      DOUBLE PRECISION EIG,Z
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION C,C2,DT,QX,RX,S,S2,V,W,X,XP,ZP
C     ..
C     .. External Functions ..
      DOUBLE PRECISION P,Q,R
      EXTERNAL P,Q,R
C     ..
C     .. External Subroutines ..
cc      EXTERNAL DXDT
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS,COS,SIN
C     ..
C     .. Common blocks ..
      COMMON /DATAF/EIG,EIGF
      COMMON /ZEE/Z
C     ..
      DT = 1.0
      X = T
      IF (.NOT.EIGF) CALL DXDT(T,DT,X)
      ZP = ABS(Z)
      XP = ZP/P(X)
      QX = Q(X)/ZP
      RX = R(X)/ZP
      V = EIG*RX + QX
      S = SIN(Y(1))
      C = COS(Y(1))
      S2 = S*S
      C2 = C*C
      YP(1) = DT*(XP*C2+V*S2)
      W = (XP-V)*S*C
      IF (Z.LT.0.0) RX = ABS(RX)
      YP(2) = DT*(-2.0*W*Y(2)+RX*S2)
      YP(3) = DT*W
      RETURN
      END SUBROUTINE F
C ----------------------------------------------------------------------
      SUBROUTINE GS(X,Y,YP)
      DOUBLE PRECISION X
      DOUBLE PRECISION Y(1),YP(1)
C     **********
C
C     This subroutine evaluates the derivative function for use with
C     integrator GERK in solving a differential equation in the form
C
C        (p(x)*y'(x))' + q(x)*y(x) = 0  on  (a,b)
C
C     for user-supplied coefficient functions P and Q.
C
C     Subprograms called
C
C       user-supplied ..... p,q
C
C     **********
C     .. Scalars in Common ..
      DOUBLE PRECISION Z
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION C,C2,QX,S,S2,XP
C     ..
C     .. External Functions ..
      DOUBLE PRECISION P,Q
      EXTERNAL P,Q
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC COS,SIN
C     ..
C     .. Common blocks ..
      COMMON /ZEE/Z
C     ..
      XP = Z/P(X)
      QX = Q(X)/Z
      S = SIN(Y(1))
      C = COS(Y(1))
      S2 = S*S
      C2 = C*C
      YP(1) = XP*C2+QX*S2
      RETURN
      END SUBROUTINE GS
C ----------------------------------------------------------------------
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
      END SUBROUTINE INTPOL
C ----------------------------------------------------------------------
      SUBROUTINE ZCOUNT(A,B,A1,A2,B1,B2,JPAIRS,PAIRS,JSUM)
      INTEGER JPAIRS,JSUM
      DOUBLE PRECISION A,B,A1,A2,B1,B2
      DOUBLE PRECISION PAIRS(2*JPAIRS)
C     **********
C
C     This subroutine counts the zeros, over specified subintervals, of
C     the solutions of a second order differential equation in the form
C
C        (p(x)*y'(x))' + q(x)*y(x) = 0  on  (a,b)
C
C     for user-supplied coefficient functions P and Q.  This count in
C     turn corresponds to the number of zeros, in the interior of (A,B),
C     of the first eigenfunction of the related Sturm-Liouville problem
C     whose (semidefinite) weight function vanishes identically in the
C     subintervals.  The problem is restricted to be regular.
C
C     The applicable initial condition depends upon three cases.
C
C        Case 1 -- On a subinterval with left endpoint A,
C                     A1*Y(A) + A2*Y'(A)*P(A) = 0.
C
C        Case 2 -- On a subinterval with right endpoint B,
C                     B1*Y(B) + B2*Y'(B)*P(B) = 0.
C
C        Case 3 -- On a subinterval with neither A nor B as endpoint,
C                     Y(XAA) = 0, where XAA is the left endpoint.
C
C     The SUBROUTINE statement is
C
C       SUBROUTINE zcount(a,b,a1,a2,b1,b2,jpairs,pairs,jsum)
C
C     where
C
C       A and B are input variables defining the full interval.
C         A must be less than B.
C
C       A1 and A2 are input variables set to prescribe the initial
C         condition at A (Case 1).
C
C       B1 and B2 are input variables set to prescribe the initial
C         condition at B (Case 2).
C
C       JPAIRS is an integer input variable set to the number of
C         specified subintervals of (A,B).
C
C       PAIRS is an input array of length 2*JPAIRS whose successive
C         ordered element pairs specify the subintervals.
C
C       JSUM is an integer output variable set to the total zero count.
C
C     Subprograms called
C
C       sleign-supplied ... epslon,g,gerk
C
C       user-supplied ..... p,q
C
C     This version dated 5/18/89.
C     Burton S. Garbow, Argonne National Laboratory
C
C     **********
C     .. Scalars in Common ..
      DOUBLE PRECISION Z
C     ..
C     .. Local Scalars ..
      INTEGER I,J,LFLAG,MF,ML
      DOUBLE PRECISION EPS,EPSMIN,H,ONE,PI,PSUM,PX,QX,RT,RTSAV,
     1     T,U,V,W,WSAV,X,X50,XAA,XBB,XSAV,ZAV
C     ..
C     .. Local Arrays ..
      INTEGER IWORK(5)
      DOUBLE PRECISION DS(98),GERROR(1),PS(99),QS(99),WORK(11),Y(1)
C     ..
C     .. External Functions ..
      DOUBLE PRECISION P,Q
      EXTERNAL P,Q
C     ..
C     .. External Subroutines ..
cc      EXTERNAL GS,GERK
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS,ATAN,INT,SIGN,SQRT
C     ..
C     .. Common blocks ..
      COMMON /ZEE/Z
C     ..
C     Set constants EPSMIN, the computer unit roundoff error, and PI.
C     (Variable ONE set to 1.0 eases precision conversion.)
C
      ONE = 1.0
      EPSMIN = EPSLON(ONE)
      PI = 4.0*ATAN(ONE)
C
C     Set relative and absolute error tolerances for GERK.
C
      EPS = SQRT(EPSMIN)
C
      JSUM = 0
      DO 70 J = 1,JPAIRS
         XAA = PAIRS(2*J-1)
         XBB = PAIRS(2*J)
C        DO (CALCULATE-MODIFIED-PRUFER-TRANSFORMATION-CONSTANT)
C
C     Evaluate P, Q to obtain preliminary information about the
C     differential equation.
C
C           DO (SAMPLE-COEFFICIENTS)
               X50 = 0.5*((XBB+XAA)+(XBB-XAA)*EPSMIN)
               PX = P(X50)
               QX = Q(X50)
               PS(50) = PX
               QS(50) = QX/PX
C
C     MF and ML are the least and greatest index values, respectively.
C
               XSAV = X50
               H = 0.9/40.0
               DO 10 I=49,1,-1
                  IF (I.GE.10) T = H*(I-50)
                  IF (I.LT.10) T = T - 0.75*(1.0+T)
                  X = 0.5*((XBB+XAA)+(XBB-XAA)*T)
                  PX = P(X)
                  QX = Q(X)
                  PS(I) = PX
                  QS(I) = QX/PX
                  DS(I) = XSAV - X
                  XSAV = X
                  MF = I
   10             CONTINUE
               XSAV = X50
               DO 20 I=51,99
                  IF (I.LE.90) T = H*(I-50)
                  IF (I.GT.90) T = T + 0.75*(1.0-T)
                  X = 0.5*((XBB+XAA)+(XBB-XAA)*T)
                  PX = P(X)
                  QX = Q(X)
                  PS(I) = PX
                  QS(I) = QX/PX
                  DS(I-1) = X - XSAV
                  XSAV = X
                  ML = I - 1
   20             CONTINUE
C              END (SAMPLE-COEFFICIENTS)
C           DO (ESTIMATE-PHASE-ANGLE-CHANGE)
C
C     U measures the total length of subintervals where Q/P .gt. 0.0.
C     ZAV is the average value of sqrt(Q*P) over those subintervals.
C
               U = 0.0
               ZAV = 0.0
               WSAV = QS(MF)
               IF (WSAV.GT.0.0) THEN
                  RTSAV = SQRT(WSAV)
               ELSE
                  RTSAV = 0.0
                  END IF
               DO 30 I=MF+1,ML
                  W = QS(I)
                  IF (W.GT.0.0) THEN
                     U = U + DS(I-1)
                     RT = SQRT(W)
                  ELSE
                     RT = 0.0
                     END IF
                  IF (W.EQ.0.0 .OR. WSAV.EQ.0.0 .OR.
     1                W.EQ.SIGN(W,WSAV)) THEN
                     V = RT + RTSAV
                  ELSE
                     V = (W*RT+WSAV*RTSAV)/ABS(W-WSAV)
                     END IF
                  WSAV = W
                  RTSAV = RT
                  PSUM = DS(I-1)*V
                  IF (U.GT.0.0) ZAV = ZAV + PSUM*(PS(I)+PS(I-1))
   30             CONTINUE
               ZAV = 0.25*ZAV
C              END (ESTIMATE-PHASE-ANGLE-CHANGE)
            Z = 1.0
            IF (U.GT.0.0) Z = ZAV/U
C           END (CALCULATE-MODIFIED-PRUFER-TRANSFORMATION-CONSTANT)
         LFLAG = 1
         IF (XAA.EQ.A) THEN
C
C           Case 1 ----------
C
            Y(1) = PI/2.0
            IF (A1.NE.0.0) Y(1) = ATAN(-Z*A2/A1)
            IF (Y(1).LT.0.0) Y(1) = Y(1) + PI
   40       CONTINUE
               CALL GERK(GS,1,Y,XAA,XBB,EPS,EPS,LFLAG,GERROR,WORK,IWORK)
               IF (LFLAG.EQ.3) GO TO 40
            JSUM = JSUM + INT((Y(1)+ABS(EPS))/PI)
         ELSE IF (XBB.EQ.B) THEN
C
C           Case 2 ----------
C
            Y(1) = PI/2.0
            IF (B1.NE.0.0) Y(1) = ATAN(-Z*B2/B1)
            IF (Y(1).GT.0.0) Y(1) = Y(1) - PI
   50       CONTINUE
               CALL GERK(GS,1,Y,XBB,XAA,EPS,EPS,LFLAG,GERROR,WORK,IWORK)
               IF (LFLAG.EQ.3) GO TO 50
            JSUM = JSUM - INT((Y(1)-ABS(EPS))/PI)
         ELSE
C
C           Case 3 ----------
C
            Y(1) = 0.0
   60       CONTINUE
               CALL GERK(GS,1,Y,XAA,XBB,EPS,EPS,LFLAG,GERROR,WORK,IWORK)
               IF (LFLAG.EQ.3) GO TO 60
            JSUM = JSUM + INT((Y(1)+ABS(EPS))/PI)
            END IF
   70    CONTINUE
      RETURN
      END SUBROUTINE ZCOUNT

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
cc      EXTERNAL GERKS
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
      END SUBROUTINE GERK
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
cc      DOUBLE PRECISION EPSLON
cc      EXTERNAL EPSLON
C     ..
C     .. EXTERNAL SUBROUTINES ..
cc      EXTERNAL FEHL
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
      END SUBROUTINE GERKS
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
      END SUBROUTINE FEHL
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
      END FUNCTION EPSLON

      end module sleignmd
