C Revision of sleig2md.f1 4/9/98
C Removes some EXTERNAL statements that refer to module procedures
C Marked with cc? at start of line
      module sleig2md
      private
      public::sleign2

      contains

C     THIS VERSION OF SLEIGN2 IS DATED SEPT. 8, 1994.
C DOUBLE PRECISION version (created using NAG F77 Tools)
C
C     This program is for an equation of the form
C
C               -(p(x)*y'(x))' + q(x)*y(x) = eig*w(x)*y(x)
C
      subroutine SLEIGN2(A,B,INTAB,P0ATA,QFATA,P0ATB,QFATB,A1,A2,B1,B2,
     +                   NUMEIG,EIG,TOL,IFLAG,ISLFUN,SLFUN,SINGATA,
     +                   SINGATB,CIRCLA,CIRCLB,OSCILA,OSCILB)
C     **********
C
C     This subroutine is designed for the calculation of a specified
C     eigenvalue, EIG, of a Sturm-Liouville problem in the form
C
C       -(p(x)*y'(x))' + q(x)*y(x) = eig*w(x)*y(x)  on (a,b)
C
C     for user-supplied coefficient functions p, q, and w.
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
C       SUBROUTINE sleign2(a,b,intab,p0ata,qfata,p0atb,qfatb,a1,a2,b1,b2,
C                         numeig,eig,tol,iflag,islfun,slfun,
C                         singata,singatb,circla,circlb,oscila,oscilb)
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
C         condition at A.  Their values are ignored if A is singular.
C
C       B1 and B2 are input variables set to prescribe the boundary
C         condition at B.  Their values are ignored if B is singular.
C
C       NUMEIG is an integer variable.  On input, it should be set to
C         the index of the desired eigenvalue (increasing sequence where
C         index 0 corresponds to the smallest nonnegative eigenvalue).
C         On output, it is unchanged unless the problem (apparently)
C         lacks eigenvalue NUMEIG, in which case it is reset to the
C         index of the largest eigenvalue.
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
C       IFLAG is an integer output variable set as follows.
C
C         IFLAG = 1  - successful problem solution.
C         IFLAG = 2  - integrator tolerance cannot be reduced.
C         IFLAG = 3  - no more improvement.
C         IFLAG = 4  - RAY and EIG fail to agree after 5 tries.
C         IFLAG = 6  - in SECANT-METHOD, ABS(DE) .LT. EPSMIN .
C         IFLAG = 7  - iterations are stuck in a loop.
C         IFLAG = 8  - number of iterations has reached the set limit.
C         IFLAG = 9  - residual truncation error dominates.
C         IFLAG = 10 - improper input parameters.
C         IFLAG = 11 - NUMEIG exceeds actual highest eigenvalue index.
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
C     This version dated 6/21/94.
C     Paul B. Bailey, Albuquerque, New Mexico
C     Burton S. Garbow, Park Forest, Illinois
C
C     **********
C     .. Scalars in Common ..
      double precision AA,ASAV,BB,BSAV,C1,C2,DTHDAA,DTHDBB,EIGSAV,
     +                 EPSMIN,HPI,PI,TMID,TSAVEL,TSAVER,TWOPI,Z
      integer IND,INTSAV,MDTHZ
      logical ADDD
C     ..
C     .. Arrays in Common ..
      double precision TEE(100),ZEE(100)
      integer JAY(100)
C     ..
C     .. Local Scalars ..
      double precision AAA,AAAA,AAS,ALFA,BALLPK,BBB,BBBB,BBS,BESTEIG,
     +                 BESTEST,BETA,C,CHNG,CHNGLIM,CL,CR,DAV,DE,DEDW,
     +                 DEN,DERIVL,DERIVR,DIST,DT,DTHDA,DTHDAAX,DTHDB,
     +                 DTHDBBX,DTHDE,DTHDEA,DTHDEB,DTHETA,DTHOLD,
     +                 DTHOLDY,DTHZ,DUM,E,EEE,EIGLO,EIGLT,EIGPI,EIGRT,
     +                 EIGUP,EL,ELIMA,ELIMB,ELIMUP,EMAX,EMIN,EOLD,EPS,
     +                 EPSM,ER1,ER2,ESTERR,EU,FLO,FMAX,FNEW,FOLD,FUP,
     +                 GUESS,HU,HV,OLDEST,OLDRAY,ONE,PIN,PSIL,PSIPL,
     +                 PSIPR,PSIR,PUP,PVP,PX,QAV,QX,RATL1,RATL2,RATL3,
     +                 RATR1,RATR2,RATR3,RAY,REMZ,RLX,SL,SL1,SL2,SL3,
     +                 SQL,SQR,SR,SR1,SR2,SR3,SUM,SUM0,T,T1,T2,T3,TAU,
     +                 THA,THB,THRESH,TMP,TS,U,UL,UR,UT,V,WAV,WL,WU,WX,
     +                 X,X50,XAA,XBB,XMID,XSAV,XT,ZAV
      integer I,IA,IB,IE,IMAX,IMID,IMIN,IOUT,JFLAG,JJL,JJR,K,KFLAG,
     +        LOOP2,LOOP3,MF,ML,MP,NEIG,NITER,NRAY
      logical AOK,BOK,BRACKT,CHNGEPS,CONVRG,ENDA,ENDB,EXIT,FIRSTT,FYNYT,
     +        FYNYT1,LCIRCA,LCIRCB,LIMA,LIMB,LIMUP,LOGIC,NEWTON,NEWTONF,
     +        OLDNEWT,ONEDIG,OSCA,OSCB,PRIN,SINGA,SINGB,THEGT0,THELT0
C     ..
C     .. Local Arrays ..
      double precision DELT(99),DS(99),ERL(3),ERR(3),PS(99),PSS(99),
     +                 QS(99),WS(99),XS(99),YL(3),YR(3),YZL(3),YZR(3)
C     ..
C     .. External Functions ..
      double precision P,Q,W
      external P,Q,W
C     ..
C     .. External Subroutines ..
cc      external AABB,ALFBET,DXDT,EIGFCN,ESTPAC,INTEG,SETMID,THUM,UV
C     ..
C     .. Intrinsic Functions ..
      intrinsic ABS,ATAN,COS,EXP,INT,LOG,MAX,MIN,SIGN,SIN,SQRT
C     ..
C     .. Common blocks ..
      common /DATADT/ASAV,BSAV,C1,C2,INTSAV
      common /DATAF/EIGSAV,IND
      common /PIE/PI,TWOPI,HPI
      common /RNDOFF/EPSMIN
      common /TDATA/AA,TMID,BB,DTHDAA,DTHDBB,ADDD,MDTHZ
      common /TEEZ/TEE
      common /TSAVE/TSAVEL,TSAVER
      common /ZEE/Z
      common /ZEEZ/JAY,ZEE
C     ..
C     .. Scalar Arguments ..
      double precision A,A1,A2,B,B1,B2,CIRCLA,CIRCLB,EIG,OSCILA,OSCILB,
     +                 P0ATA,P0ATB,QFATA,QFATB,SINGATA,SINGATB,TOL
      integer IFLAG,INTAB,ISLFUN,NUMEIG
C     ..
C     .. Array Arguments ..
      double precision SLFUN(9)
C     ..
C     Set constants EPSMIN, the computer unit roundoff error, and PI.
C     (Variable ONE set to 1.0 eases precision conversion.)
C Initialize PSS (John Pryce)
C
      ONE = 1.0D0
      EPSMIN = EPSLON(ONE)
      PI = 4.0D0*ATAN(ONE)
      TWOPI = 2.0D0*PI
      HPI = 0.5D0*PI
      do 3 I=1,99
        PSS(I) = 0.0D0
  3   continue
C
C     Set output device number.
C
      IOUT = 6
C
C     Check input parameters for errors.  If errors, return IFLAG=2.
C
      LOGIC = 1 .le. INTAB .and. INTAB .le. 4 .and.
     +        P0ATA*QFATA*P0ATB*QFATB .ne. 0.0D0
      if (INTAB.eq.1) LOGIC = LOGIC .and. A .lt. B
      if (.not.LOGIC) then
         IFLAG = 10
         go to 150

      end if
C
C     Set PRIN = .true. to trigger trace printout of successive steps.
C
      PRIN = .FALSE.
      if (TOL.lt.0.0D0) PRIN = .TRUE.
C
C     Set EPS to the (initial) integration accuracy.
C
      EPS = 0.0001D0
C
C     Set logical variables.
C
      AOK = INTAB .lt. 3 .and. P0ATA .lt. 0.0D0 .and. QFATA .gt. 0.0D0
      BOK = (INTAB.eq.1 .or. INTAB.eq.3) .and. P0ATB .lt. 0.0D0 .and.
     +      QFATB .gt. 0.0D0
      SINGA = SINGATA .gt. 0.0D0
      SINGB = SINGATB .gt. 0.0D0
      LCIRCA = CIRCLA .gt. 0.0D0
      LCIRCB = CIRCLB .gt. 0.0D0
      OSCA = OSCILA .gt. 0.0D0
      OSCB = OSCILB .gt. 0.0D0
      EIGPI = NUMEIG*PI
      NEIG = NUMEIG - 1
C
C     Initial C1 and C2, used in the mapping between X and T intervals.
C
      C1 = 1.0D0
      C2 = 0.0D0
C     DO (SAVE-INPUT-DATA)
      ASAV = A
      BSAV = B
      INTSAV = INTAB
      TAU = ABS(TOL)
C        END (SAVE-INPUT-DATA)
C
C     Initialize the arrays JAY and ZEE if either end is oscillatory.
C
      if (OSCA .or. OSCB) then
         do 5 K = 1,100
            JAY(K) = 0
            ZEE(K) = 1.0D0
    5    continue
      end if
C
C     Evaluate P, Q, W to obtain preliminary information about the
C     differential equation.
C
C     DO (SAMPLE-COEFFICIENTS)
      THRESH = 1.0D+17
   10 continue
      call DXDT(EPSMIN,TMP,X50)
      XS(50) = X50
      TS = EPSMIN
      PX = P(X50)
      QX = Q(X50)
      WX = W(X50)
      PS(50) = PX
      QS(50) = QX/PX
      WS(50) = WX/PX
C
C     DAV,QAV,WAV are used in special case estimation when NUMEIG = 0,1.
C     EMIN = min(Q/W), achieved at X for index value IMIN.
C     EMAX = max(Q/W), achieved at X for index value IMAX.
C     MF and ML are the least and greatest index values, respectively.
C
      DAV = 0.0D0
      QAV = 0.0D0
      WAV = 0.0D0
      XSAV = X50
      EMIN = 0.0D0
      if (QX.ne.0.0D0) EMIN = QX/WX
      EMAX = EMIN
      IMIN = 50
      IMAX = 50
      do 20 I = 49,1,-1
         T = TFROMI(I)
         call DXDT(T,TMP,X)
         XS(I) = X
         PX = P(X)
         QX = Q(X)
         WX = W(X)
         PS(I) = PX
         QS(I) = QX/PX
         WS(I) = WX/PX
         DS(I) = XSAV - X
         DELT(I) = 0.5D0* (TS-T)
         DAV = DAV + DS(I)
         QAV = QAV + DS(I)* (0.5D0* (QS(I)+QS(I+1))-QAV)/DAV
         WAV = WAV + DS(I)* (0.5D0* (WS(I)+WS(I+1))-WAV)/DAV
         XSAV = X
         TS = T
C
C     Try to avoid overflow by stopping when functions are large near A.
C
         FYNYT = (ABS(WX)+ABS(QX)+1.0D0/ABS(PX)) .le. THRESH
         if (QX.ne.0.0D0 .and. QX/WX.lt.EMIN) then
            EMIN = QX/WX
            IMIN = I
         end if

         if (QX.ne.0.0D0 .and. QX/WX.gt.EMAX) then
            EMAX = QX/WX
            IMAX = I
         end if

         MF = I
         if (.not.FYNYT) go to 30
   20 continue
   30 continue
      AAA = T
      if (.not.SINGA) AAA = -1.0D0
      XSAV = X50
      do 40 I = 51,99
         T = TFROMI(I)
         call DXDT(T,TMP,X)
         XS(I) = X
         PX = P(X)
         QX = Q(X)
         WX = W(X)
         PS(I) = PX
         QS(I) = QX/PX
         WS(I) = WX/PX
         DS(I-1) = X - XSAV
         DELT(I-1) = 0.5D0* (T-TS)
         DAV = DAV + DS(I-1)
         QAV = QAV + DS(I-1)* (0.5D0* (QS(I-1)+QS(I))-QAV)/DAV
         WAV = WAV + DS(I-1)* (0.5D0* (WS(I-1)+WS(I))-WAV)/DAV
         XSAV = X
         TS = T
C
C     Try to avoid overflow by stopping when functions are large near B.
C
         FYNYT1 = (ABS(QX)+ABS(WX)+1.0D0/ABS(PX)) .le. THRESH
         if (QX.ne.0.0D0 .and. QX/WX.lt.EMIN) then
            EMIN = QX/WX
            IMIN = I
         end if

         if (QX.ne.0.0D0 .and. QX/WX.gt.EMAX) then
            EMAX = QX/WX
            IMAX = I
         end if

         ML = I - 1
         if (.not.FYNYT1) go to 50
   40 continue
   50 continue
      BBB = T
      if (.not.SINGB) BBB = 1.0D0
      LOGIC = C1 .eq. 1.0D0 .and. (.not.FYNYT .or. .not.FYNYT1)
C
C     Modify (T,X) transformation corresponding to truncated interval.
C
      if (LOGIC) then
         C1 = 0.5D0* (BBB-AAA)
         C2 = 0.5D0* (AAA+BBB)
         go to 10

      end if

      if (OSCA .or. OSCB) call THUM(MF,ML,XS)
C
C     Here we try to determine 'sigma0'.  Initially, we will be
C     satisfied to determine eliml and elimr, the limiting values
C     of q/w, if they exist.
C
      LIMA = .FALSE.
      if (SINGA .and. .not.LCIRCA) then
         RATL1 = QS(MF)/WS(MF)
         RATL2 = QS(MF+1)/WS(MF+1)
         RATL3 = QS(MF+2)/WS(MF+2)
         SL1 = RATL1/ (XS(MF+1)-XS(MF))
         SL2 = RATL2/ (XS(MF+2)-XS(MF+1))
         SL3 = RATL3/ (XS(MF+3)-XS(MF+2))
         if (ABS(SL2).ge.ABS(SL1) .and. ABS(SL2).le.ABS(SL3)) then
            ELIMA = RATL1
            LIMA = PS(MF) .eq. PS(MF+1)
            WRITE (021,FMT=*) ' There is a limit at a = ',ELIMA
            write (21,FMT=*) ' There is a limit at a = ',ELIMA
         end if

      end if

      LIMB = .FALSE.
      if (SINGB .and. .not.LCIRCB) then
         RATR1 = QS(ML)/WS(ML)
         RATR2 = QS(ML-1)/WS(ML-1)
         RATR3 = QS(ML-2)/WS(ML-2)
         SR1 = RATR1/ (XS(ML)-XS(ML-1))
         SR2 = RATR2/ (XS(ML-1)-XS(ML-2))
         SR3 = RATR3/ (XS(ML-2)-XS(ML-3))
         if (ABS(SR2).ge.ABS(SR1) .and. ABS(SR2).le.ABS(SR3)) then
            ELIMB = RATR1
            LIMB = PS(ML) .eq. PS(ML-1)
            WRITE (021,FMT=*) ' There is a limit at b = ',ELIMB
            write (21,FMT=*) ' There is a limit at b = ',ELIMB
         end if

      end if

      LIMUP = .FALSE.
      ELIMUP = EMAX
      if (LIMA .or. LIMB) then
         LIMUP = .TRUE.
         if (.not.LIMB) then
            ELIMUP = ELIMA
         else if (.not.LIMA) then
            ELIMUP = ELIMB
         else
            ELIMUP = MIN(ELIMA,ELIMB)
         end if

         write (21,FMT=*) ' The continuous spectrum has a lower '
         write (21,FMT=*) '   bound, sigma0 = ',ELIMUP
      end if
C        END (SAMPLE-COEFFICIENTS)
      PIN = EIGPI + PI
      if (EIG.eq.0.0D0) then
C        DO (ESTIMATE-EIG)
         SUM0 = 0.0D0
         if (OSCA .or. OSCB) then
            EEE = 0.0D0
C              DO (ESTIMATE-PHASE-ANGLE-CHANGE)
            call ESTPAC(.FALSE.,MF,ML,EEE,SUM0,QS,WS,DS,DELT,PS,PSS,TAU,
     +                  IA,IB,JJL,JJR,SUM,U,UT,ZAV)
C                 END (ESTIMATE-PHASE-ANGLE-CHANGE)
            SUM0 = SUM
         end if

         EEE = MIN(ELIMUP,EMAX)
C           DO (ESTIMATE-PHASE-ANGLE-CHANGE)
         call ESTPAC(.FALSE.,MF,ML,EEE,SUM0,QS,WS,DS,DELT,PS,PSS,TAU,IA,
     +               IB,JJL,JJR,SUM,U,UT,ZAV)
C              END (ESTIMATE-PHASE-ANGLE-CHANGE)
   17    continue
         if (.not.LIMUP .and. ABS(SUM).ge.
     +       10.D0*MAX(1.0D0,ABS(PIN))) then
            if (SUM.ge.10.D0*PIN) then
               if (EEE.ge.1.0D0) then
                  EEE = EEE/10.D0
               else if (EEE.lt.-1.0D0) then
                  EEE = 10.D0*EEE
               else
                  EEE = EEE - 1.0D0
               end if

            else if (SUM.le.10.D0*PIN) then
               if (EEE.le.-1.0D0) then
                  EEE = EEE/10.D0
               else if (EEE.gt.1.0D0) then
                  EEE = 10.D0*EEE
               else
                  EEE = EEE + 1.0D0
               end if

            else
               go to 27

            end if
C             DO (ESTIMATE-PHASE-ANGLE-CHANGE)
            call ESTPAC(.FALSE.,MF,ML,EEE,SUM0,QS,WS,DS,DELT,PS,PSS,TAU,
     +                  IA,IB,JJL,JJR,SUM,U,UT,ZAV)
C                END (ESTIMATE-PHASE-ANGLE-CHANGE)
            go to 17

         end if

   27    continue
C
         EU = EEE
         WU = SUM
         if (SUM.ge.PIN) then
            EL = EU
            WL = WU
   60       continue
            if (WL.ge.PIN) then
               EU = EL
               WU = WL
               EEE = EL - ((WL-PIN+3.0D0)/U)**2 - 1.0D0
C                    DO (ESTIMATE-PHASE-ANGLE-CHANGE)
               call ESTPAC(.FALSE.,MF,ML,EEE,SUM0,QS,WS,DS,DELT,PS,PSS,
     +                     TAU,IA,IB,JJL,JJR,SUM,U,UT,ZAV)
C                       END (ESTIMATE-PHASE-ANGLE-CHANGE)
               EL = EEE
               WL = SUM
               go to 60

            end if

         else
            EL = EEE
            WL = SUM
         end if

         if (LIMUP .and. WU.lt.PIN) then
            EEE = ELIMUP
         else
            if (U.eq.0.0D0) then
               EEE = EMAX + 1.0D0
C                 DO (ESTIMATE-PHASE-ANGLE-CHANGE)
               call ESTPAC(.FALSE.,MF,ML,EEE,SUM0,QS,WS,DS,DELT,PS,PSS,
     +                     TAU,IA,IB,JJL,JJR,SUM,U,UT,ZAV)
C                    END (ESTIMATE-PHASE-ANGLE-CHANGE)
               EU = EEE
               WU = SUM
            end if

   70       continue
            if (WU.le.PIN) then
               EL = EU
               WL = WU
               EEE = EU + ((PIN-WU+3.0D0)/U)**2 + 1.0D0
C                    DO (ESTIMATE-PHASE-ANGLE-CHANGE)
               call ESTPAC(.FALSE.,MF,ML,EEE,SUM0,QS,WS,DS,DELT,PS,PSS,
     +                     TAU,IA,IB,JJL,JJR,SUM,U,UT,ZAV)
C                       END (ESTIMATE-PHASE-ANGLE-CHANGE)
               EU = EEE
               WU = SUM
               go to 70

            end if

   80       continue
            if (ABS(IMAX-IMIN).ge.2 .and. EU.le.EMAX) then
               IE = (IMAX+IMIN)/2
               EEE = QS(IE)/WS(IE)
C                    DO (ESTIMATE-PHASE-ANGLE-CHANGE)
               call ESTPAC(.FALSE.,MF,ML,EEE,SUM0,QS,WS,DS,DELT,PS,PSS,
     +                     TAU,IA,IB,JJL,JJR,SUM,U,UT,ZAV)
C                       END (ESTIMATE-PHASE-ANGLE-CHANGE)
               if (SUM.gt.PIN) then
                  IMAX = IE
                  WU = SUM
                  EU = EEE
               else
                  IMIN = IE
                  WL = SUM
                  EL = EEE
               end if

               go to 80

            end if
C
C     Improve approximation for EIG using bisection or secant method.
C     Substitute 'ballpark' estimate if approximation grows too large.
C
            DEDW = (EU-EL)/ (WU-WL)
            FOLD = 0.0D0
            if (INTAB.eq.1) BALLPK = (PIN/ (A-B))**2
            if (INTAB.eq.1) write (21,FMT=*) ' BALLPK = ',BALLPK
            LOGIC = .TRUE.
   90       continue
            if (LOGIC) then
               LOGIC = (WL.lt.PIN-1.0D0 .or. WU.gt.PIN+1.0D0)
               EEE = EL + DEDW* (PIN-WL)
               FNEW = MIN(PIN-WL,WU-PIN)
               if (FNEW.gt.0.4D0*FOLD .or. FNEW.le.1.0D0) EEE = 0.5D0*
     +             (EL+EU)
               if (INTAB.eq.1 .and. ABS(EEE).gt.1.0D3*BALLPK) then
                  EEE = BALLPK
                  go to 100

               else if (INTAB.ne.1 .and. ABS(EEE).gt.1.0D6) then
                  EEE = 1.0D0
                  go to 100

               else
                  FOLD = FNEW
C                       DO (ESTIMATE-PHASE-ANGLE-CHANGE)
                  call ESTPAC(.FALSE.,MF,ML,EEE,SUM0,QS,WS,DS,DELT,PS,
     +                        PSS,TAU,IA,IB,JJL,JJR,SUM,U,UT,ZAV)
C                          END (ESTIMATE-PHASE-ANGLE-CHANGE)
                  if (SUM.lt.PIN) then
                     EL = EEE
                     WL = SUM
                  else
                     EU = EEE
                     WU = SUM
                  end if

                  DEDW = (EU-EL)/ (WU-WL)
                  go to 90

               end if

            end if

         end if
C           END (ESTIMATE-EIG)
      end if

  100 continue
      GUESS = EIG
      if (LIMUP .and. EEE.ge.ELIMUP) EEE = ELIMUP - 0.01D0
C     DO (SET-INITIAL-INTERVAL-AND-MATCHPOINT)
      if (GUESS.ne.0.0D0) then
         EEE = EIG
C           DO (ESTIMATE-PHASE-ANGLE-CHANGE)
         call ESTPAC(.FALSE.,MF,ML,EEE,SUM0,QS,WS,DS,DELT,PS,PSS,TAU,IA,
     +               IB,JJL,JJR,SUM,U,UT,ZAV)
C              END (ESTIMATE-PHASE-ANGLE-CHANGE)
      end if
C
C     Choose initial interval as large as possible that avoids overflow.
C     JJL and JJR are boundary indices for nonnegativity of EIG*W-Q.
C
      AA = -1.0D0
      if (SINGA) AA = TFROMI(JJL)
      BB = 1.0D0
      if (SINGB) BB = TFROMI(JJR)
      AA = MIN(-0.01D0,AA)
      BB = MAX(0.01D0,BB)
      AA = MIN(AA,-0.9D0)
      BB = MAX(BB,0.9D0)
      if (OSCA) AA = -0.9999D0
      if (OSCB) BB = 0.9999D0
C
C     Determine boundary values ALFA and BETA for theta at A and B.
C
      Z = 1.0D0
      call ALFBET(A,INTAB,AA,A1,A2,EEE,P0ATA,QFATA,SINGA,LCIRCA,ALFA,
     +            KFLAG,DERIVL)
      call ALFBET(B,INTAB,BB,B1,B2,EEE,P0ATB,QFATB,SINGB,LCIRCB,BETA,
     +            JFLAG,DERIVR)
      if (SINGB) BETA = PI - BETA
C
C     Take boundary conditions into account in estimation of EIG.
C
      PIN = EIGPI + BETA - ALFA
      if (OSCA) PIN = PIN + ALFA
      if (OSCB) PIN = PIN + PI - BETA
      if (GUESS.eq.0.0D0) then
         EEE = EL + DEDW* (PIN-WL)
         if (.not. (OSCA.or.OSCB) .and.
     +       ABS(EEE).gt.1000.0D0) EEE = SIGN(1000.0D0,EEE)
         if (INTAB.eq.1 .and. ABS(EEE).gt.1.0D3*BALLPK) EEE = BALLPK
      end if
C        DO (ESTIMATE-PHASE-ANGLE-CHANGE)
      call ESTPAC(OSCA .or. OSCB,MF,ML,EEE,SUM0,QS,WS,DS,DELT,PS,PSS,
     +            TAU,IA,IB,JJL,JJR,SUM,U,UT,ZAV)
C           END (ESTIMATE-PHASE-ANGLE-CHANGE)
C
C     Choose the constant Z .
C
      if (U.gt.0.0D0) Z = ZAV/UT
C
C     Reset boundary values ALFA and BETA .
C
      call ALFBET(A,INTAB,AA,A1,A2,EEE,P0ATA,QFATA,SINGA,LCIRCA,ALFA,
     +            KFLAG,DERIVL)
      call ALFBET(B,INTAB,BB,B1,B2,EEE,P0ATB,QFATB,SINGB,LCIRCB,BETA,
     +            JFLAG,DERIVR)
      if (SINGB) BETA = PI - BETA
      if (PRIN) write (IOUT,FMT='(A,E22.14,A,E22.14)') ' alfa=',ALFA,
     +    '   beta=',BETA
      write (21,FMT='(A,E22.14,A,E22.14)') ' alfa=',ALFA,'   beta=',BETA
C
C     Special case formula for estimation of EIG when NUMEIG = 0,1.
C
      if (U.eq.0.0D0 .and. NUMEIG.le.0 .and. (BETA+EIGPI).lt.ALFA) then
C           XBC = MAX(-1.0/TAN(ALFA),1.0/TAN(BETA))
C           EEE = -(XBC*XBC-QAV)/WAV
C           DEDW = XBC*(1.0+XBC*XBC)/WAV
      end if
C
C     Choose initial matching point TMID .
C
      IMID = 50
      TMID = 0.5D0* (AA+BB)
      if (PRIN) write (IOUT,FMT='(A,E15.7,A,F11.8,A,E15.7)') ' estim=',
     +    EEE,'  tmid=',TMID,'  z=',Z
      if (PRIN) write (IOUT,FMT='(A,F11.8,A,F11.8,A,F11.8,A,F11.8)')
     +    ' aaa=',AAA,'  aa=',AA,'  bb=',BB,'  bbb=',BBB
      write (21,FMT='(A,E15.7,A,F11.8,A,E15.7)') ' estim=',EEE,
     +  '  tmid=',TMID,'  z=',Z
      write (21,FMT='(A,F11.8,A,F11.8,A,F11.8,A,F11.8)') ' aaa=',AAA,
     +  '  aa=',AA,'  bb=',BB,'  bbb=',BBB
C        END (SET-INITIAL-INTERVAL-AND-MATCHPOINT)
      if (EIG.eq.0.0D0 .and. LIMUP .and. EEE.ge.ELIMUP) EEE = ELIMUP -
     +    0.01D0
C     DO (RESET-TMID)
      call SETMID(MF,ML,EEE,QS,WS,IMID,TMID)
C        END (RESET-TMID)
      if (OSCA .or. OSCB) then
         Z = 1.0D0
C        DO (PREP-ZEEZ)
         do 85 I = 1,100
            TEE(I) = 1.0D0
            if (JAY(I).ne.0) TEE(I) = TFROMI(JAY(I))
   85    continue
C           END (PREP-ZEEZ)
      end if

      if (IFLAG.eq.-1) then
         SLFUN(1) = TMID
         SLFUN(2) = AA
         SLFUN(3) = ALFA
         SLFUN(5) = BB
         SLFUN(6) = BETA + EIGPI
         SLFUN(9) = Z
         ADDD = .FALSE.
         MDTHZ = 0
         IFLAG = 1
         return

      end if
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
      EXIT = .FALSE.
      FIRSTT = .TRUE.
      LOOP2 = 0
      LOOP3 = 0
      MP = 0
      BESTEIG = EIG
      BESTEST = 1.D+9
      OLDEST = BESTEST
      NEWTONF = .FALSE.
      CHNGEPS = .FALSE.
      EPSM = EPSMIN
      AAAA = AAA
      BBBB = BBB
      ENDA = .FALSE.
      ENDB = .FALSE.
      TSAVEL = -1.0D0
      TSAVER = 1.0D0
  110 continue
C     DO (INITIAL-IZE)
      BRACKT = .FALSE.
      CONVRG = .FALSE.
      THELT0 = .FALSE.
      THEGT0 = .FALSE.
      EIGLO = EMIN - 1.0D0
      EIGLT = 0.0D0
      EIGRT = 0.0D0
      EIGUP = EMAX + 1.0D0
      if (LIMUP) EIGUP = MIN(EMAX,ELIMUP)
      DTHOLD = 1.0D0
      IFLAG = 1
C        END (INITIAL-IZE)
      write (21,FMT=*)
      write (21,FMT=*) '---------------------------------------------'
      write (21,FMT=*) ' INITIAL GUESS FOR EIG = ',EIG
      write (21,FMT=*) ' aa,bb = ',AA,BB
C     DO UNTIL(CONVRG .OR. EXIT)
      do 120 NITER = 1,40
         WRITE (021,FMT=*)
         WRITE (021,FMT=*) ' ******************** '
C        DO (SET-TMID-AND-BOUNDARY-CONDITIONS)
         WRITE (021,FMT=*) ' set tmid and boundary conditions '
         V = EIG*WS(IMID) - QS(IMID)
C           IF (V.LE.0.0) DO (RESET-TMID)
         if (V.le.0.0D0) call SETMID(MF,ML,EIG,QS,WS,IMID,TMID)
C              END (RESET-TMID)
C           DO (RESET-BOUNDARY-CONDITIONS)
         DERIVL = 0.0D0
         if (SINGA) call ALFBET(A,INTAB,AA,A1,A2,EIG,P0ATA,QFATA,.TRUE.,
     +                          LCIRCA,ALFA,KFLAG,DERIVL)
         DERIVR = 0.0D0
         if (SINGB) then
            call ALFBET(B,INTAB,BB,B1,B2,EIG,P0ATB,QFATB,.TRUE.,LCIRCB,
     +                  BETA,JFLAG,DERIVR)
            BETA = PI - BETA
         end if

         if (PRIN) write (IOUT,FMT='(A,E22.14,A,E22.14)') ' alfa=',ALFA,
     +       '   beta=',BETA
C              END (RESET-BOUNDARY-CONDITIONS)
C
C     Check that boundary conditions can be satisfied at singular
C     endpoints.  If not, try for slightly altered EIG consistent
C     with boundary conditions.
C
         if (LIMUP .and. EIG.ne.GUESS .and. .not.BRACKT) then
  115       continue
            KFLAG = 1
            if (SINGA .and. .not.LCIRCA) call ALFBET(A,INTAB,AA,A1,A2,
     +          EIG,P0ATA,QFATA,.TRUE.,.FALSE.,TMP,KFLAG,TMP)
            JFLAG = 1
            if (SINGB .and. .not.LCIRCB) call ALFBET(B,INTAB,BB,B1,B2,
     +          EIG,P0ATB,QFATB,.TRUE.,.FALSE.,TMP,JFLAG,TMP)
            if ((KFLAG.ne.1.or.JFLAG.ne.1) .and.
     +          (THELT0.and.EIGLO.lt.ELIMUP)) then
               EIGUP = ELIMUP
               EIG = ELIMUP - EPSMIN
               go to 115

            end if

         end if
C           END (SET-TMID-AND-BOUNDARY-CONDITIONS)
C        DO (OBTAIN-DTHETA-WITH-ONE-CORRECT-DIGIT)
         if (PRIN) write (IOUT,FMT='(/A,E22.14,A,E10.3,A,E10.3)')
     +       ' guess=',EIG,'  eps=',EPS,'  tmid=',TMID
C           DO (INTEGRATE-FOR-DTHETA)
C              DO (SET-INITIAL-CONDITIONS)
         THA = ALFA
         DTHDEA = DERIVL
         DTHDAA = 0.0D0
         if (SINGA .and. .not.LCIRCA) then
            call DXDT(AA,DT,X)
            PX = P(X)/Z
            QX = Q(X)/Z
            WX = W(X)/Z
            C = EIG*WX - QX
            DTHDAA = - (COS(ALFA)**2/PX+C*SIN(ALFA)**2)*DT
C
C     Two special cases for DTHDAA .
C
            if (C.ge.0.0D0 .and. P0ATA.lt.0.0D0 .and.
     +          QFATA.lt.0.0D0) DTHDAA = DTHDAA + ALFA*DT/ (X-A)
            if (C.ge.0.0D0 .and. P0ATA.gt.0.0D0 .and.
     +          QFATA.gt.0.0D0) DTHDAA = DTHDAA + (ALFA-0.5D0*PI)*DT/
     +                                   (X-A)
         end if

         THB = BETA
         DTHDEB = -DERIVR
         DTHDBB = 0.0D0
         if (SINGB .and. .not.LCIRCB) then
            call DXDT(BB,DT,X)
            PX = P(X)/Z
            QX = Q(X)/Z
            WX = W(X)/Z
            C = EIG*WX - QX
            DTHDBB = - (COS(BETA)**2/PX+C*SIN(BETA)**2)*DT
C
C     Two special cases for DTHDBB .
C
            if (C.ge.0.0D0 .and. P0ATB.lt.0.0D0 .and.
     +          QFATB.lt.0.0D0) DTHDBB = DTHDBB + (PI-BETA)*DT/ (B-X)
            if (C.ge.0.0D0 .and. P0ATB.gt.0.0D0 .and.
     +          QFATB.gt.0.0D0) DTHDBB = DTHDBB + (0.5D0*PI-BETA)*DT/
     +                                   (B-X)
         end if
C                 END (SET-INITIAL-CONDITIONS)
         EIGSAV = EIG
C                                             T
C     YL = (theta,d(theta)/d(eig),d(theta)/da)
C
         YL(1) = ALFA
         YL(2) = DTHDEA
         YL(3) = 0.0D0
C
         call INTEG(AA,THA,DTHDAA,DTHDEA,TMID,A1,A2,EPS,YL,ERL,LCIRCA,
     +              AOK,SINGA,OSCA,IFLAG)
         if (IFLAG.eq.5) then
            IFLAG = 51
            WRITE (021,FMT=*) ' IFLAG = 51 '
            EXIT = .TRUE.
            go to 130

         end if

         DTHDA = DTHDAA*EXP(-2.0D0*YL(3))
C                                             T
C     YR = (theta,d(theta)/d(eig),d(theta)/db)
C
         YR(1) = BETA + EIGPI - PI
         YR(2) = DTHDEB
         YR(3) = 0.0D0
C
         call INTEG(BB,THB,DTHDBB,DTHDEB,TMID,B1,B2,EPS,YR,ERR,LCIRCB,
     +              BOK,SINGB,OSCB,IFLAG)
         if (IFLAG.eq.5) then
            IFLAG = 52
            WRITE (021,FMT=*) ' IFLAG = 52 '
            EXIT = .TRUE.
            go to 130

         end if

         DTHDB = DTHDBB*EXP(-2.0D0*YR(3))
C
         ER1 = ERL(1) - ERR(1)
         ER2 = ERL(2) - ERR(2)
C
         if (OSCA .or. OSCB) then
            Z = 1.0D0
            call DXDT(TMID,TMP,XT)
            call UV(XT,U,PUP,V,PVP,HU,HV)
            EIGSAV = 0.0D0
            call INTEG(AA,THA,DTHDAA,DTHDEA,TMID,A1,A2,EPS,YZL,ERL,
     +                 LCIRCA,AOK,SINGA,OSCA,IFLAG)
            if (IFLAG.eq.5) then
               IFLAG = 53
               WRITE (021,FMT=*) ' IFLAG = 53 '
               EXIT = .TRUE.
               go to 130

            end if

            call INTEG(BB,THB,DTHDBB,DTHDEB,TMID,B1,B2,EPS,YZR,ERR,
     +                 LCIRCB,BOK,SINGB,OSCB,IFLAG)
            if (IFLAG.eq.5) then
               IFLAG = 54
               WRITE (021,FMT=*) ' IFLAG = 54 '
               EXIT = .TRUE.
               go to 130

            end if

            EIGSAV = EIG
            DTHZ = YZR(1) - YZL(1)
            MDTHZ = DTHZ/PI
            REMZ = DTHZ - MDTHZ*PI
            if (DTHZ.lt.0.0D0 .and. REMZ.lt.0.0D0) then
               MDTHZ = MDTHZ - 1
               REMZ = REMZ + PI
            end if

            if (REMZ.gt.3.14D0) MDTHZ = MDTHZ + 1
         end if
C
C     DTHETA measures theta difference from left and right integrations.
C
C              DO (FORM-DTHETA)
         DTHETA = YL(1) - YR(1) - EIGPI
         if (OSCA .or. OSCB) DTHETA = DTHETA + MDTHZ*PI
         DTHDE = YL(2) - YR(2)
C                 END (FORM-DTHETA)
         ONEDIG = ABS(ER1) .le. 0.5D0*ABS(DTHETA) .and.
     +            ABS(ER2) .le. 0.5D0*ABS(DTHDE)
         write (21,FMT=*)
         write (21,FMT=*) ' EIG,DTHETA,ONEDIG = ',EIG,DTHETA,ONEDIG
         FIRSTT = .FALSE.
C              END (INTEGRATE-FOR-DTHETA)
         CHNGEPS = .FALSE.
         CONVRG = .FALSE.
         OLDNEWT = NEWTON
         NEWTON = ABS(DTHETA) .lt. 0.06D0 .and. BRACKT
         if (NEWTON) then
            ONEDIG = ONEDIG .or. ABS(DTHETA+ER1) .lt. 0.5D0*DTHOLD
            if (.not.ONEDIG .and. EPS.gt.EPSM) then
               EPS = MAX(0.01D0*EPS,EPSM)
               CHNGEPS = .TRUE.
            end if

         end if

         if (PRIN) write (IOUT,FMT='(A,E15.7,A,E15.7)') ' dtheta=',
     +       DTHETA,'   dthde=',DTHDE
         if (PRIN) write (IOUT,FMT='(/A,E15.7,A,E15.7)') ' thetal=',
     +       YL(1),'   thetar=',YR(1)
C           END (OBTAIN-DTHETA-WITH-ONE-CORRECT-DIGIT)
         if (.not.ONEDIG) then
            EXIT = .TRUE.
            write (21,FMT=*) ' NOT ONEDIG '
            go to 130

         end if
C        DO (SET-BRACKET-DATA)
         WRITE (021,FMT=*) ' set-bracket '
         if (DTHETA*DTHDE.gt.0.0D0) then
            if (.not.THEGT0 .or. EIG.le.EIGUP) then
               THEGT0 = .TRUE.
               EIGUP = EIG
               FUP = DTHETA
               EIGRT = EIG - DTHETA/DTHDE
            end if

         else
            if (.not.THELT0 .or. EIG.ge.EIGLO) then
               THELT0 = .TRUE.
               EIGLO = EIG
               FLO = DTHETA
               EIGLT = EIG - DTHETA/DTHDE
            end if

         end if
C
C     EIG is bracketed when both THEGT0=.true. and THELT0=.true.
C
         BRACKT = THELT0 .and. THEGT0
         if (PRIN) write (IOUT,FMT='(A,E22.14,A,E22.14)') ' eigrt=',
     +       EIGRT,'  eigup=',EIGUP
         if (PRIN) write (IOUT,FMT='(A,E22.14,A,E22.14)') ' eiglt=',
     +       EIGLT,'  eiglo=',EIGLO
C           END (SET-BRACKET-DATA)
         if (BRACKT) LOOP2 = 0
C        DO (TEST-FOR-CONVERGENCE)
C
C     Measure convergence after adding separate contributions to error.
C
         T1 = ABS(DTHETA)/ABS(DTHDE)
         T2 = (1.0D0+AA)*ABS(DTHDA)/ABS(DTHDE)
         T3 = (1.0D0-BB)*ABS(DTHDB)/ABS(DTHDE)
         if (.not. (AOK.or.SINGA) .or. (LCIRCA.and.
     +       .not.OSCA)) T2 = 0.0D0
         if (.not. (BOK.or.SINGB) .or. (LCIRCB.and.
     +       .not.OSCB)) T3 = 0.0D0
         ESTERR = T1 + T2 + T3
         ESTERR = ESTERR/MAX(ONE,ABS(EIG))
         CONVRG = ESTERR .le. TAU .and. NEWTON
         write (21,FMT=*) ' T1,T2,T3 = ',T1,T2,T3
         write (21,FMT=*) ' EPS,ESTERR = ',EPS,ESTERR
         write (21,FMT=*) ' ONEDIG,BRACKT,NEWTON,CONVRG = ',ONEDIG,
     +     BRACKT,NEWTON,CONVRG
         if (ESTERR.lt.BESTEST) then
            BESTEIG = EIG
            BESTEST = ESTERR
            write (21,FMT=*) ' BESTEIG,BESTEST = ',BESTEIG,BESTEST
         end if

         if (THEGT0) write (21,FMT=*) '         EIGUP = ',EIGUP
         if (THELT0) write (21,FMT=*) '         EIGLO = ',EIGLO
         if (PRIN) write (IOUT,FMT='(A,L2)') ' converge=',CONVRG
         if (PRIN .and. .not.CONVRG) write (IOUT,
     +       FMT='(A,E15.7)') ' estim. acc.=',ESTERR
C           END (TEST-FOR-CONVERGENCE)
         if (CONVRG) then
            WRITE (021,FMT=*) ' number of iterations was ',NITER
            WRITE (021,FMT=*)
     +        '-----------------------------------------------'
            go to 130

         else
            if (NEWTON) then
               if (OLDNEWT .and. ABS(DTHETA).ge.0.5D0*ABS(DTHOLD)) then
                  write (21,FMT=*) ' NEWTON DID NOT IMPROVE EIG '
                  NEWTONF = .TRUE.
                  LOOP3 = LOOP3 + 1
               else
                  ENDA = T2 .gt. T1 .and. AA .gt. AAAA
                  ENDB = T3 .gt. T1 .and. BB .lt. BBBB
                  if (ENDA .or. ENDB) then
                     NEWTON = .FALSE.
                  else if ((T2+T3).gt.T1 .and.
     +                     (AA.le.AAAA.and.BB.ge.BBBB)) then
                     WRITE (021,FMT=*)
     +                 ' RESIDUAL TRUNCATION ERROR DOMINATES '
                     EXIT = .TRUE.
                     IFLAG = 9
                     go to 130

                  end if

               end if

               if (NEWTONF .or. ENDA .or. ENDB) then
                  EXIT = .TRUE.
                  go to 130

               end if
C              DO (NEWTON'S-METHOD)
               WRITE (021,FMT=*) ' Newton''s method '
               RLX = 1.2D0
               if (BRACKT) RLX = 1.0D0
               EIG = EIG - RLX*DTHETA/DTHDE
C                 END (NEWTON'S-METHOD)
            else if (BRACKT) then
               WRITE (021,FMT=*) ' bracket '
C              DO (SECANT-METHOD)
               WRITE (021,FMT=*) ' do secant method '
               FMAX = MAX(-FLO,FUP)
               EOLD = EIG
               EIG = 0.5D0* (EIGLO+EIGUP)
               if (FMAX.le.1.5D0) then
                  U = -FLO/ (FUP-FLO)
                  DIST = EIGUP - EIGLO
                  EIG = EIGLO + U*DIST
                  V = MIN(EIGLT,EIGRT)
                  if (EIG.le.V) EIG = 0.5D0* (EIG+V)
                  V = MAX(EIGLT,EIGRT)
                  if (EIG.ge.V) EIG = 0.5D0* (EIG+V)
                  DE = EIG - EOLD
                  if (ABS(DE).lt.EPSMIN) then
                     TOL = ABS(DE)/MAX(ONE,ABS(EIG))
                     IFLAG = 6
                     EXIT = .TRUE.
                     go to 130

                  end if

               end if
C                 END (SECANT-METHOD)
            else
C              DO (TRY-FOR-BRACKET)
               LOOP2 = LOOP2 + 1
               if (LOOP2.le.9 .and. DTHETA.lt.0.0D0) then
                  MP = MIN(MP,INT(-DTHETA/PI))
                  write (21,FMT=*) ' MP = ',MP
               end if

               if (LOOP2.gt.9) then
                  IFLAG = 12
                  EXIT = .TRUE.
                  go to 130

               end if

               if (EIG.eq.EEE) then
                  if (GUESS.ne.0.0D0) DEDW = 1.0D0/DTHDE
                  CHNG = -0.6D0* (DEDW+1.0D0/DTHDE)*DTHETA
                  if (EIG.ne.0.0D0 .and. ABS(CHNG).gt.
     +                0.1D0*ABS(EIG)) CHNG = -0.1D0*SIGN(EIG,DTHETA)
               else
                  CHNG = -1.2D0*DTHETA/DTHDE
C
C     Limit change in EIG to a factor of 10.
C
                  if (ABS(CHNG).gt. (1.0D0+10.0D0*ABS(EIG))) then
                     CHNG = SIGN(1.0D0+10.0D0*ABS(EIG),CHNG)
                  else if (ABS(EIG).ge.1.0D0 .and.
     +                     ABS(CHNG).lt.0.1D0*ABS(EIG)) then
                     CHNG = 0.1D0*SIGN(EIG,CHNG)
                  end if

                  if (DTHOLD.lt.0.0D0 .and. LIMUP .and.
     +                CHNG.gt. (ELIMUP-EIG)) then
                     CHNG = 0.95D0* (ELIMUP-EIG)
                     if (CHNG.lt.EPSMIN) then
                        WRITE (021,FMT=*) ' elimup,eig = ',ELIMUP,EIG
                        write (21,FMT=*) ' IN BRACKET, CHNG.LT.EPSMIN '
                        NUMEIG = NEIG - INT(-DTHETA/PI)
                        WRITE (021,FMT=*) ' new numeig = ',NUMEIG
                        write (21,FMT=*) ' new numeig = ',NUMEIG
                        IFLAG = 11
                        EXIT = .TRUE.
                     end if

                  end if

               end if

               EOLD = EIG
               CHNGLIM = 2.0D0*ESTERR*MAX(ONE,ABS(EIG))
               if (ABS(DTHETA).lt.0.06D0 .and.
     +             ABS(CHNG).gt.CHNGLIM) CHNG = SIGN(CHNGLIM,CHNG)
               EIG = EIG + CHNG
C                 END (TRY-FOR-BRACKET)
            end if

         end if

         if (IFLAG.eq.11) go to 130
         if (NITER.ge.3 .and. DTHOLDY.eq.DTHETA) then
            IFLAG = 7
            EXIT = .TRUE.
            go to 130

         end if

         DTHOLDY = DTHOLD
         DTHOLD = DTHETA
         WRITE (021,FMT=*) ' number of iterations was ',NITER
         WRITE (021,FMT=*)
     +     '-----------------------------------------------'
  120 continue
      IFLAG = 8
      EXIT = .TRUE.
  130 continue
      TOL = BESTEST
      EIG = BESTEIG
      if (EXIT) then
         write (21,FMT=*) ' EXIT '
         if (FIRSTT) then
            if (IFLAG.eq.51 .or. IFLAG.eq.53) then
               if (AA.lt.-0.1D0) then
                  WRITE (021,FMT=*)' FIRST COMPLETE INTEGRATION FAILED.'
                  write (21,FMT=*)
     +              ' FIRST COMPLETE INTEGRATION FAILED. '
                  if (AA.eq.-1.0D0) go to 150
                  AAS = AA
                  call AABB(AA,-ONE)
                  write (21,FMT=*) ' aa MOVED FROM ',AAS,' IN TO ',AA
                  EXIT = .FALSE.
                  go to 110

               else
                  write (21,FMT=*) ' aa.ge.-0.1 '
                  IFLAG = 13
                  go to 150

               end if

            else if (IFLAG.eq.52 .or. IFLAG.eq.54) then
               if (BB.gt.0.1D0) then
                  WRITE (021,FMT=*)' FIRST COMPLETE INTEGRATION FAILED.'
                  write (21,FMT=*)
     +              ' FIRST COMPLETE INTEGRATION FAILED. '
                  if (BB.eq.1.0D0) go to 150
                  BBS = BB
                  call AABB(BB,-ONE)
                  write (21,FMT=*) ' bb MOVED FROM ',BBS,' IN TO ',BB
                  EXIT = .FALSE.
                  go to 110

               else
                  write (21,FMT=*) ' bb.le.0.1 '
                  IFLAG = 14
                  go to 150

               end if

            end if

         else if (IFLAG.eq.51 .or. IFLAG.eq.53) then
            WRITE (021,FMT=*) ' A COMPLETE INTEGRATION FAILED. '
            write (21,FMT=*) ' A COMPLETE INTEGRATION FAILED. '
            if (CHNGEPS .and. EPS.lt.0.002D0) then
               EPS = 5.0D0*EPS
               EPSM = EPS
               write (21,FMT=*) ' EPS INCREASED TO ',EPS
            else
               AAS = AA
               call AABB(AA,-ONE)
               AAAA = AA
               write (21,FMT=*) ' aa MOVED FROM ',AAS,' IN TO ',AA
            end if

            EXIT = .FALSE.
            go to 110

         else if (IFLAG.eq.52 .or. IFLAG.eq.54) then
            WRITE (021,FMT=*) ' A COMPLETE INTEGRATION FAILED. '
            write (21,FMT=*) ' A COMPLETE INTEGRATION FAILED. '
            if (CHNGEPS .and. EPS.lt.0.002D0) then
               EPS = 5.0D0*EPS
               EPSM = EPS
               write (21,FMT=*) ' EPS INCREASED TO ',EPS
            else
               BBS = BB
               call AABB(BB,-ONE)
               BBBB = BB
               write (21,FMT=*) ' bb MOVED FROM ',BBS,' IN TO ',BB
            end if

            EXIT = .FALSE.
            go to 110

         else if (IFLAG.eq.6) then
            WRITE (021,FMT=*) ' IN SECANT, CHNG.LT.EPSMIN '
            write (21,FMT=*) ' IN SECANT, CHNG.LT.EPSMIN '
            go to 140

         else if (IFLAG.eq.7) then
            WRITE (021,FMT=*) ' DTHETA IS REPEATING '
            write (21,FMT=*) ' DTHETA IS REPEATING '
            go to 140

         else if (IFLAG.eq.8) then
            WRITE (021,FMT=*) ' NUMBER OF ITERATIONS REACHED SET LIMIT '
            write (21,FMT=*) ' NUMBER OF ITERATIONS REACHED SET LIMIT '
            go to 140

         else if (IFLAG.eq.9) then
            write (21,FMT=*) ' RESIDUAL TRUNCATION ERROR DOMINATES '
            go to 140

         else if (IFLAG.eq.11) then
            WRITE (021,FMT=*) ' IN TRY FOR BRACKET, CHNG.LT.EPSMIN '
            write (21,FMT=*) ' IN TRY FOR BRACKET, CHNG.LT.EPSMIN '
            go to 150

         else if (IFLAG.eq.12) then
            WRITE (021,FMT=*) ' FAILED TO GET A BRACKET. '
            write (21,FMT=*) ' FAILED TO GET A BRACKET. '
            go to 140

         else if (NEWTONF .or. .not.ONEDIG) then
            if (LOOP3.ge.3) then
               write (21,FMT=*) ' NEWTON IS NOT GETTING ANYWHERE '
               NEWTONF = .FALSE.
               IFLAG = 3
               go to 140

            end if

            if (EPS.gt.EPSM .and. BESTEST.lt.OLDEST) then
               CHNGEPS = .TRUE.
               EPS = 0.2D0*EPS
               WRITE (021,FMT=*) ' EPS REDUCED TO ',EPS
               write (21,FMT=*) ' EPS REDUCED TO ',EPS
               EXIT = .FALSE.
               NEWTON = .FALSE.
               OLDEST = BESTEST
               go to 110

            else
               if (EPS.le.EPSM) then
                  WRITE (021,FMT=*) ' EPS CANNOT BE REDUCED FURTHER. '
                  write (21,FMT=*) ' EPS CANNOT BE REDUCED FURTHER. '
                  IFLAG = 2
                  go to 140

               else
                  WRITE (021,FMT=*) ' no more improvement '
                  write (21,FMT=*) ' NO MORE IMPROVEMENT '
                  IFLAG = 3
                  go to 140

               end if

            end if

         else if (ENDA) then
            AAS = AA
            call AABB(AA,ONE)
            AA = MAX(AA,AAA)
            WRITE (021,FMT=*) ' aa MOVED OUT TO ',AA
            write (21,FMT=*) ' aa MOVED FROM ',AAS,' OUT TO ',AA
            ENDA = .FALSE.
            EXIT = .FALSE.
            go to 110

         else if (ENDB) then
            BBS = BB
            call AABB(BB,ONE)
            BB = MIN(BB,BBB)
            WRITE (021,FMT=*) ' bb MOVED OUT TO ',BB
            write (21,FMT=*) ' bb MOVED FROM ',BBS,' OUT TO ',BB
            ENDB = .FALSE.
            EXIT = .FALSE.
            go to 110

         end if

      end if

  140 continue
C
C     If CONVRG is false, check that any truncation error might possibly
C     be reduced or that the integrations might be done more accurately.
C
      if (.not.CONVRG) then
         if (EPS.gt.EPSM) then
            CHNGEPS = .TRUE.
            EPS = 0.2D0*EPS
            WRITE (021,FMT=*) ' EPS REDUCED TO ',EPS
            write (21,FMT=*) ' EPS REDUCED TO ',EPS
            EXIT = .FALSE.
            NEWTON = .FALSE.
            OLDEST = BESTEST
            go to 110

         else if (AA.gt.AAAA .and. T2.gt.0.0D0 .and. T2.ge.T3) then
            AAS = AA
            call AABB(AA,ONE)
            AA = MAX(AA,AAA)
            WRITE (021,FMT=*) ' aa MOVED OUT TO ',AA
            write (21,FMT=*) ' aa MOVED FROM ',AAS,' OUT TO ',AA
            EXIT = .FALSE.
            go to 110

         else if (BB.lt.BBBB .and. T3.gt.0.0D0 .and. T3.ge.T2) then
            BBS = BB
            call AABB(BB,ONE)
            BB = MIN(BB,BBB)
            WRITE (021,FMT=*) ' bb MOVED OUT TO ',BB
            write (21,FMT=*) ' bb MOVED FROM ',BBS,' OUT TO ',BB
            EXIT = .FALSE.
            go to 110

         end if

      end if

      if (PRIN) write (IOUT,FMT='(A,I7,A,E22.14,A,E10.3)') ' numeig=',
     +    NUMEIG,'  eig=',EIG,'  tol=',TOL
C     DO (COMPUTE-EIGENFUNCTION-DATA)
C
C     Convert from T to X values, fill 7 of first 9 locations of SLFUN.
C
      call DXDT(TMID,TMP,XMID)
      call DXDT(AA,TMP,XAA)
      call DXDT(BB,TMP,XBB)
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
      DTHDAAX = 0.0D0
      YL(1) = 0.0D0
      YL(2) = 0.0D0
      YL(3) = 0.0D0
      call INTEG(AA,THA,DTHDAAX,DTHDEA,TMID,A1,A2,EPS,YL,ERL,LCIRCA,AOK,
     +           SINGA,OSCA,JFLAG)
      THB = BETA
      DTHDBBX = 0.0D0
      call INTEG(BB,THB,DTHDBBX,DTHDEB,TMID,B1,B2,EPS,YR,ERR,LCIRCB,BOK,
     +           SINGB,OSCB,JFLAG)
      YR(1) = YR(1) + EIGPI
      SL = SIN(YL(1))
      SR = SIN(YR(1))
      CL = COS(YL(1))
      CR = COS(YR(1))
      UL = (YL(2)-DTHDEA*EXP(-2.0D0*YL(3)))*Z
      UR = (YR(2)-DTHDEB*EXP(-2.0D0*YR(3)))*Z
      DUM = 0.5D0*LOG(UL-UR)
      SLFUN(4) = -YL(3) - DUM
      SLFUN(7) = -YR(3) - DUM
C        END (COMPUTE-EIGENFUNCTION-DATA)
C     DO (CHECK-MATCHING-VALUES-OF-EIGENFUNCTION)
C
C     Perform final check on EIG. Return IFLAG = 400+
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
      ADDD = PSIL*PSIR .lt. 0.0D0 .and. PSIPL*PSIPR .lt. 0.0D0
      RAY = EIG + (PSIL*PSIPL-PSIR*PSIPR)/ (SQL-SQR)
      if (PRIN) then
         write (IOUT,FMT='(A,E22.14)') ' ray=',RAY
         write (IOUT,FMT='(A,E22.14,A,E22.14)') ' psil=',PSIL,'  psir=',
     +     PSIR
         write (IOUT,FMT='(A,E22.14,A,E22.14)') ' psipl=',PSIPL,
     +     '  psipr=',PSIPR
         write (IOUT,FMT='(A,E22.14,A,E22.14)') ' sql=',SQL,'  sqr=',SQR
      end if
C        END (CHECK-MATCHING-VALUES-OF-EIGENFUNCTION)
C
C     If next condition is .true., then something is apparently wrong
C     with the accuracy of EIG. Bisect and go through the loop again.
C
      if (ABS(RAY-EIG).gt.2.0D0*TOL*MAX(ONE,ABS(EIG))) then
         NRAY = NRAY + 1
         write (21,FMT=*) ' NRAY, RAY = ',NRAY,RAY
         if (RAY.eq.OLDRAY .or. NRAY.gt.5) then
            IFLAG = 400 + IFLAG
            go to 150

         end if

         EIG = 0.5D0* (EIG+RAY)
         OLDRAY = RAY
         go to 110

      end if

      if (NRAY.le.5) IFLAG = IFLAG + 100
C     DO (GENERATE-EIGENFUNCTION-VALUES)
      call EIGFCN(EIGPI,A1,A2,B1,B2,AOK,SINGA,LCIRCA,OSCA,BOK,SINGB,
     +            LCIRCB,OSCB,SLFUN,ISLFUN)
C        END (GENERATE-EIGENFUNCTION-VALUES)
  150 continue
      write (21,FMT=*) ' IFLAG = ',IFLAG
      return

      end subroutine SLEIGN2
      subroutine AABB(TEND,OUT)
C     **********
C
C     This subroutine moves aa or bb further out or closer in.
C     TEND is either aa or bb.  If OUT is positive, TEND is moved
C     further out, and if OUT is negative, TEND is moved closer in.
C
C     **********
C     .. Local Scalars ..
      double precision D,TST,TSTS
      integer I
C     ..
C     .. Intrinsic Functions ..
      intrinsic ABS,SIGN
C     ..
C     .. Scalar Arguments ..
      double precision OUT,TEND
C     ..
      if (OUT.lt.0.0D0 .and. ABS(TEND).le.0.9D0) then
         TST = 0.5D0*ABS(TEND)
      else
         D = 0.1D0
         TSTS = 0.9D0
         do 10 I = 1,20
            TST = 1.0D0 - D
            if (TST.gt.ABS(TEND) .or. (TST.eq.ABS(TEND).and.
     +          OUT.lt.0.0D0)) go to 20
            D = 0.1D0*D
            TSTS = TST
   10    continue
   20    continue
         if (OUT.lt.0.0D0) TST = TSTS
      end if

      TEND = SIGN(TST,TEND)
      return

      end subroutine AABB
      subroutine ALFBET(XEND,INTAB,TT,COEF1,COEF2,EIG,P0,QF,SING,LCIRC,
     +                  VALUE,IFLAG,DERIV)
C     **********
C
C     This subroutine computes a boundary value for a specified endpoint
C     of the interval for a Sturm-Liouville problem in the form
C
C       -(p(x)*y'(x))' + q(x)*y(x) = eig*w(x)*y(x)  on (a,b)
C
C     for user-supplied coefficient functions P, Q, and W.  It is called
C     from SLEIGN.  Both regular and singular endpoints are treated.
C
C     Subprograms called
C
C       user-supplied ..... p,q,w
C
C       sleign-supplied ... dxdt,extrap
C
C     **********
C     .. Scalars in Common ..
      double precision Z
C     ..
C     .. Local Scalars ..
      double precision C,CD,D,HH,ONE,PI,PUP,PVP,PX,QX,T,TEMP,TTS,U,V,WX,
     +                 X,XDENOM,XNUM
      logical LOGIC
C     ..
C     .. External Functions ..
      double precision P,Q,W
      external P,Q,W
C     ..
C     .. External Subroutines ..
cc      external DXDT,EXTRAP,UV
C     ..
C     .. Intrinsic Functions ..
      intrinsic ABS,ATAN,ATAN2,SIGN,SQRT
C     ..
C     .. Common blocks ..
      common /ZEE/Z
C     ..
C     Set machine dependent constant.
C
C     PI (variable ONE set to 1.0 eases precision conversion).
C     .. Scalar Arguments ..
      double precision COEF1,COEF2,DERIV,EIG,P0,QF,TT,VALUE,XEND
      integer IFLAG,INTAB
      logical LCIRC,SING
C     ..
      ONE = 1.0D0
      PI = 4.0D0*ATAN(ONE)
C
      IFLAG = 1
      DERIV = 0.0D0
      if (.not.SING) then
         VALUE = 0.5D0*PI
         if (COEF1.ne.0.0D0) VALUE = ATAN(-Z*COEF2/COEF1)
         LOGIC = (TT.lt.0.0D0 .and. VALUE.lt.0.0D0) .or.
     +           (TT.gt.0.0D0 .and. VALUE.le.0.0D0)
         if (LOGIC) VALUE = VALUE + PI
      else if (LCIRC) then
         call DXDT(TT,TEMP,X)
         call UV(X,U,PUP,V,PVP,TEMP,TEMP)
         XNUM = COEF1*U + COEF2*V
         XDENOM = (COEF1*PUP+COEF2*PVP)/Z
         VALUE = ATAN2(XNUM,XDENOM)
         if (XNUM.lt.0.0D0) VALUE = VALUE + 2.0D0*PI
      else
         LOGIC = (INTAB.eq.2 .and. TT.gt.0.0D0) .or.
     +           (INTAB.eq.3 .and. TT.lt.0.0D0) .or. INTAB .eq. 4 .or.
     +           (P0.gt.0.0D0 .and. QF.lt.0.0D0)
         if (LOGIC) then
            T = SIGN(ONE,TT)
            TTS = TT
            call EXTRAP(T,TTS,EIG,VALUE,DERIV,IFLAG)
         else
            call DXDT(TT,TEMP,X)
            PX = P(X)/Z
            QX = Q(X)/Z
            WX = W(X)/Z
            C = 2.0D0* (EIG*WX-QX)
            if (C.lt.0.0D0) then
               VALUE = 0.0D0
               if (P0.gt.0.0D0) VALUE = 0.5D0*PI
            else
               HH = ABS(XEND-X)
               D = 2.0D0*HH/PX
               CD = C*D*HH
               if (P0.gt.0.0D0) then
                  VALUE = C*HH
                  if (CD.lt.1.0D0) VALUE = VALUE/ (1.0D0+SQRT(1.0D0-CD))
                  VALUE = VALUE + 0.5D0*PI
               else
                  VALUE = D
                  if (CD.lt.1.0D0) VALUE = VALUE/ (1.0D0+SQRT(1.0D0-CD))
               end if

            end if

         end if

      end if

      return

      end subroutine ALFBET
      subroutine DXDT(T,DT,X)
C     **********
C
C     This subroutine transforms coordinates from T on (-1,1) to
C     X on (A,B) in the solution of a Sturm-Liouville problem.
C     It is called from subroutines SLEIGN, ALFBET, F, and EXTRAP.
C
C     **********
C     .. Scalars in Common ..
      double precision A,B,C1,C2
      integer INTAB
C     ..
C     .. Local Scalars ..
      double precision U
C     ..
C     .. Intrinsic Functions ..
      intrinsic ABS
C     ..
C     .. Common blocks ..
      common /DATADT/A,B,C1,C2,INTAB
C     ..
C     .. Scalar Arguments ..
      double precision DT,T,X
C     ..
      U = C1*T + C2
      go to (10,20,30,40) INTAB

   10 continue
      DT = C1*0.5D0* (B-A)
      X = 0.5D0* ((B+A)+ (B-A)*U)
      return

   20 continue
      DT = C1*2.0D0/ (1.0D0-U)**2
      X = A + (1.0D0+U)/ (1.0D0-U)
      return

   30 continue
      DT = C1*2.0D0/ (1.0D0+U)**2
      X = B - (1.0D0-U)/ (1.0D0+U)
      return

   40 continue
      DT = C1/ (1.0D0-ABS(U))**2
      X = U/ (1.0D0-ABS(U))
      return

      end subroutine DXDT
      subroutine EIGFCN(EIGPI,A1,A2,B1,B2,AOK,SINGA,LCIRCA,OSCA,BOK,
     +                  SINGB,LCIRCB,OSCB,SLFUN,ISLFUN)
C     **********
C     **********
C     .. Scalars in Common ..
      double precision AA,BB,DTHDAA,DTHDBB,TMID
      integer MDTHZ
      logical ADDD
C     ..
C     .. Local Scalars ..
      double precision DTHDAT,DTHDBT,DTHDET,EFF,T,THT,TM
      integer I,IFLAG,J,NMID
      logical LCIRC,OK,SING
C     ..
C     .. Local Arrays ..
      double precision ERL(3),ERR(3),YL(3),YR(3)
C     ..
C     .. External Subroutines ..
cc      external INTEG
C     ..
C     .. Intrinsic Functions ..
      intrinsic EXP,SIN
C     ..
C     .. Common blocks ..
      common /TDATA/AA,TMID,BB,DTHDAA,DTHDBB,ADDD,MDTHZ
C     ..
C
C     WARNING: In this program it is assumed that the points T
C              in SLFUN all lie within the interval (AA,BB).
C
C     Calculate selected eigenfunction values by integration (over T).
C
C     .. Scalar Arguments ..
      double precision A1,A2,B1,B2,EIGPI
      integer ISLFUN
      logical AOK,BOK,LCIRCA,LCIRCB,OSCA,OSCB,SINGA,SINGB
C     ..
C     .. Array Arguments ..
      double precision SLFUN(ISLFUN+9)
C     ..
      NMID = 0
      do 10 I = 1,ISLFUN
         if (SLFUN(9+I).le.TMID) NMID = I
   10 continue
      if (NMID.gt.0) then
         T = AA
         YL(1) = SLFUN(3)
         YL(2) = 0.0D0
         YL(3) = 0.0D0
         LCIRC = LCIRCA
         OK = AOK
         SING = SINGA
         EFF = 0.0D0
         do 20 J = 1,NMID
            TM = SLFUN(J+9)
            if (TM.lt.AA .or. TM.gt.BB) then
C               WRITE (021,FMT=*) ' t.lt.aa .or. t.gt.bb '
C               stop
                SLFUN(J+9) = 0D0
            end if

            THT = YL(1)
            DTHDAT = DTHDAA*EXP(-2.0D0*EFF)
            DTHDET = YL(2)
            if (TM.gt.AA) then
               call INTEG(T,THT,DTHDAT,DTHDET,TM,A1,A2,SLFUN(8),YL,ERL,
     +                    LCIRC,OK,SING,OSCA,IFLAG)
               if (OSCA) then
                  EFF = YL(3)
               else
                  LCIRC = .FALSE.
                  SING = .FALSE.
                  EFF = EFF + YL(3)
               end if

            end if

            SLFUN(J+9) = SIN(YL(1))*EXP(EFF+SLFUN(4))
            T = TM
            if (T.gt.-1.0D0) OK = .TRUE.
            if (T.lt.-0.9D0 .and. OSCA) then
               OK = .FALSE.
               T = AA
               YL(1) = SLFUN(3)
               YL(2) = 0.0D0
               YL(3) = 0.0D0
            end if

   20    continue
      end if

      if (NMID.lt.ISLFUN) then
         T = BB
         YR(1) = SLFUN(6) - EIGPI
         YR(2) = 0.0D0
         YR(3) = 0.0D0
         LCIRC = LCIRCB
         OK = BOK
         SING = SINGB
         EFF = 0.0D0
         do 30 J = ISLFUN,NMID + 1,-1
            TM = SLFUN(J+9)
            if (TM.lt.AA .or. TM.gt.BB) then
C               WRITE (021,FMT=*) ' t.lt.aa .or. t.gt.bb '
C               stop
               SLFUN(J+9) = 0D0
            end if

            THT = YR(1)
            DTHDBT = DTHDBB*EXP(-2.0D0*EFF)
            DTHDET = YR(2)
            if (TM.lt.BB) then
               call INTEG(T,THT,DTHDBT,DTHDET,TM,B1,B2,SLFUN(8),YR,ERR,
     +                    LCIRC,OK,SING,OSCB,IFLAG)
               if (OSCB) then
                  EFF = YR(3)
               else
                  LCIRC = .FALSE.
                  SING = .FALSE.
                  EFF = EFF + YR(3)
               end if

            end if

            SLFUN(J+9) = SIN(YR(1)+EIGPI)*EXP(EFF+SLFUN(7))
            if (ADDD) SLFUN(J+9) = -SLFUN(J+9)
            T = TM
            if (T.lt.1.0D0) OK = .TRUE.
            if (T.gt.0.9D0 .and. OSCB) then
               OK = .FALSE.
               T = BB
               YR(1) = SLFUN(6) - EIGPI
               YR(2) = 0.0D0
               YR(3) = 0.0D0
            end if

   30    continue
      end if

      return

      end subroutine EIGFCN
      subroutine ESTPAC(IOSC,MF,ML,EEE,SUM0,QS,WS,DS,DELT,PS,PSS,TAU,IA,
     +                  IB,JJL,JJR,SUM,U,UT,ZAV)
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
      double precision DPSUM,DPSUMT,PSUM,RT,RTSAV,V,WSAV,WW,ZAVJ,ZAVSAV
      integer J,JJ,JSAV,MF1
C     ..
C     .. Arrays in Common ..
      double precision ZEE(100)
      integer JAY(100)
C     ..
C     .. Intrinsic Functions ..
      intrinsic ABS,MAX,MIN,SIGN,SQRT
C     ..
C     .. Common blocks ..
      common /ZEEZ/JAY,ZEE
C     ..
C     .. Scalar Arguments ..
      double precision EEE,SUM,SUM0,TAU,U,UT,ZAV
      integer IA,IB,JJL,JJR,MF,ML
      logical IOSC
C     ..
C     .. Array Arguments ..
      double precision DELT(ML),DS(ML),PS(ML),PSS(ML),QS(ML),WS(ML)
C     ..
      IA = MF
      IB = 80
C
C     SUM accumulates the integral approximation.  U measures the total
C     length of subintervals where (EIG*W-Q)/P .gt. 0.0.  ZAV is the
C     average value of sqrt((EIG*W-Q)*P) over those subintervals.
C
      if (.not.IOSC) then
         JJL = 99
         JJR = 1
         SUM = 0.0D0
         U = 0.0D0
         UT = 0.0D0
         ZAV = 0.0D0
         WSAV = EEE*WS(MF) - QS(MF)
         if (WSAV.gt.0.0D0) then
            RTSAV = SIGN(SQRT(WSAV),PS(MF))
         else
            RTSAV = 0.0D0
         end if

         do 10 J = MF + 1,ML
            WW = EEE*WS(J) - QS(J)
            if (WW.gt.0.0D0) then
               if (J.gt.80) IB = J
               U = U + DS(J-1)
               UT = UT + DELT(J-1)
               RT = SIGN(SQRT(WW),PS(J))
            else
               RT = 0.0D0
               if (U.eq.0.0D0 .and. RTSAV.eq.0.0D0 .and.
     +             IA.le.19) IA = IA + 1
            end if

            if (WW.eq.0.0D0 .or. WSAV.eq.0.0D0 .or.
     +          WW.eq.SIGN(WW,WSAV)) then
               V = RT + RTSAV
            else
               V = (WW*RT+WSAV*RTSAV)/ABS(WW-WSAV)
            end if

            WSAV = WW
            RTSAV = RT
            PSUM = DS(J-1)*V
            if (EEE.eq.0.0D0) then
               PSS(J) = PSUM
            else
               DPSUM = PSUM - PSS(J)
               DPSUMT = DPSUM*DELT(J-1)/DS(J-1)
               if (DPSUMT.gt.0.001D0*TAU) then
                  JJL = MIN(JJL,J)
                  JJR = MAX(JJR,J)
               end if

            end if

            SUM = SUM + PSUM
            if (U.gt.0.0D0) ZAV = ZAV + DELT(J-1)*V*ABS(PS(J)+PS(J-1))
   10    continue
         SUM = 0.5D0*SUM - SUM0
         ZAV = 0.25D0*ZAV
      else
         JJ = 1
         JAY(1) = MF
   20    continue
         SUM = 0.0D0
         U = 0.0D0
         UT = 0.0D0
         ZAV = 0.0D0
         ZAVJ = 0.0D0
         MF1 = JAY(JJ)
         WSAV = EEE*WS(MF1) - QS(MF1)
         if (WSAV.gt.0.0D0) then
            RTSAV = SIGN(SQRT(WSAV),PS(MF1))
         else
            RTSAV = 0.0D0
         end if

         do 30 J = MF1 + 1,ML
            WW = EEE*WS(J) - QS(J)
            if (WW.gt.0.0D0) then
               if (J.gt.80) IB = J
               U = U + DS(J-1)
               UT = UT + DELT(J-1)
               RT = SIGN(SQRT(WW),PS(J))
            else
               RT = 0.0D0
               if (U.eq.0.0D0 .and. RTSAV.eq.0.0D0 .and.
     +             IA.le.19) IA = IA + 1
            end if

            if (WW.eq.0.0D0 .or. WSAV.eq.0.0D0 .or.
     +          WW.eq.SIGN(WW,WSAV)) then
               V = RT + RTSAV
            else
               V = (WW*RT+WSAV*RTSAV)/ABS(WW-WSAV)
            end if

            WSAV = WW
            RTSAV = RT
            PSUM = DS(J-1)*V
            SUM = SUM + PSUM
            if (U.gt.0.0D0) ZAV = ZAV + DELT(J-1)*V*ABS(PS(J)+PS(J-1))
            if (U.ne.0.0D0) then
               if (ZAVJ.eq.0.0D0) JSAV = J
               ZAVJ = 0.25D0*ZAV/UT
               if (J.eq.JSAV) ZAVSAV = ZAVJ
               if (2.0D0*ZAVJ.lt.ZAVSAV .or. ZAVJ.gt.2.0D0*ZAVSAV) then
                  JJ = JJ + 1
                  JAY(JJ) = J
                  ZEE(JJ) = 0.5D0* (ZAVJ+ZAVSAV)
                  go to 40

               end if

            end if

   30    continue
   40    continue
         if (J.gt.ML) then
            JJ = JJ + 1
            JAY(JJ) = ML
            ZEE(JJ) = 0.5D0* (ZAVJ+ZAVSAV)
         end if

         if (J.lt.ML) go to 20
         SUM = 0.5D0*SUM
         ZAV = 0.25D0*ZAV
      end if

      IB = IB + 1
      return

      end subroutine ESTPAC
      subroutine EXTRAP(T,TT,EIG,VALUE,DERIV,IFLAG)
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
C       sleign-supplied ... dxdt,intpol
C
C     **********
C     .. Scalars in Common ..
      double precision Z
C     ..
C     .. Local Scalars ..
      double precision ANS,CTN,ERROR,PROD,PX,QX,T1,TEMP,WX,X
      integer KGOOD
C     ..
C     .. Local Arrays ..
      double precision FN1(5),XN(5)
C     ..
C     .. External Functions ..
      double precision P,Q,W
      external P,Q,W
C     ..
C     .. External Subroutines ..
cc      external DXDT,INTPOL
C     ..
C     .. Intrinsic Functions ..
      intrinsic ABS,ATAN,SQRT,TAN
C     ..
C     .. Common blocks ..
      common /ZEE/Z
C     ..
C     .. Scalar Arguments ..
      double precision DERIV,EIG,T,TT,VALUE
      integer IFLAG
C     ..
      IFLAG = 1
      KGOOD = 0
      T1 = TT
   10 continue
      call DXDT(T1,TEMP,X)
      PX = P(X)/Z
      QX = Q(X)/Z
      WX = W(X)/Z
      PROD = -PX* (EIG*WX-QX)
      if (PROD.le.0.0D0) then
         T1 = 0.5D0* (T1+T)
         if ((1.0D0+ (T1-T)**2).gt.1.0D0) go to 10
         IFLAG = 5
         return

      else
         KGOOD = KGOOD + 1
         XN(KGOOD) = T1
         FN1(KGOOD) = ATAN(1.0D0/SQRT(PROD))
         T1 = 0.5D0* (T+T1)
         if (KGOOD.lt.5) go to 10
      end if

      T1 = 0.01D0
      call INTPOL(5,XN,FN1,T,T1,3,ANS,ERROR)
      VALUE = ABS(ANS)
      CTN = 1.0D0/TAN(VALUE)
      DERIV = 0.5D0*PX*WX/CTN/ (1.0D0+CTN**2)
      TT = XN(1)
      return

      end subroutine EXTRAP
      subroutine F(U,Y,YP)
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
C       sleign-supplied ... dxdt
C
C     **********
C     .. Scalars in Common ..
      double precision EIG,Z
      integer IND
C     ..
C     .. Local Scalars ..
      double precision C,C2,DT,QX,S,S2,T,TH,V,WW,WX,X,XP
C     ..
C     .. External Functions ..
      double precision P,Q,W
      external P,Q,W
C     ..
C     .. External Subroutines ..
cc      external DXDT
C     ..
C     .. Intrinsic Functions ..
      intrinsic COS,MOD,SIN
C     ..
C     .. Common blocks ..
      common /DATAF/EIG,IND
      common /ZEE/Z
C     ..
C     .. Scalar Arguments ..
      double precision U
C     ..
C     .. Array Arguments ..
      double precision Y(2),YP(3)
C     ..
      if (MOD(IND,2).eq.1) then
         T = U
         TH = Y(1)
      else
         T = Y(1)
         TH = U
      end if

      call DXDT(T,DT,X)
      XP = Z/P(X)
      QX = Q(X)/Z
      WX = W(X)/Z
      V = EIG*WX - QX
      S = SIN(TH)
      C = COS(TH)
      S2 = S*S
      C2 = C*C
      YP(1) = DT* (XP*C2+V*S2)
      if (IND.eq.1) then
         WW = (XP-V)*S*C
         YP(2) = DT* (-2.0D0*WW*Y(2)+WX*S2)
         YP(3) = DT*WW
      else if (IND.eq.2) then
         YP(2) = YP(2)/YP(1)
         YP(3) = YP(3)/YP(1)
         YP(1) = 1.0D0/YP(1)
      else if (IND.eq.3) then
      else
         YP(1) = 1.0D0/YP(1)
      end if

      return

      end subroutine F
      double precision function FF(ALFLAM)
C     **********
C     **********
C     .. Scalars in Common ..
      double precision CC,EIGSAV,HPI,PI,THETU,THETV,TWOPI,UB,UL,UR,VB,
     +                 VL,VR,Z
      integer IND,INDD
C     ..
C     .. Local Scalars ..
      double precision AA,BB,DTHDAA,DTHDBB,DTHDEA,DTHDEB,DUM,EPS,LAMBDA,
     +                 PVPB,PVPL,PVPR,RHOB,RHOL,RHOR,THA,THB,THL,THR,
     +                 TMID
      integer LFLAG
      logical AOK,BOK,LCIRCA,LCIRCB,OSCA,OSCB,SINGA,SINGB
C     ..
C     .. Local Arrays ..
      double precision ERR(3),Y(3)
C     ..
C     .. External Subroutines ..
cc      external INTEG
C     ..
C     .. Intrinsic Functions ..
      intrinsic COS,EXP,SIN
C     ..
C     .. Common blocks ..
      common /DATAF/EIGSAV,INDD
      common /EPP2/CC,UL,UR,VL,VR,UB,VB,IND
      common /PIE/PI,TWOPI,HPI
      common /THET/THETU,THETV
      common /ZEE/Z
C     ..
C
C     (THIS ROUTINE IS BEING MODIFIED FOR PERIODIC PROBLEMS WHICH
C     ARE NOT NECESSARILY REGULAR, BUT IS NOT YET COMPLETE.)
C
C     .. Scalar Arguments ..
      double precision ALFLAM
C     ..
      AOK = .TRUE.
      LCIRCA = .FALSE.
      SINGA = .FALSE.
      OSCA = .FALSE.
      BOK = .TRUE.
      LCIRCB = .FALSE.
      SINGB = .FALSE.
      OSCB = .FALSE.
C
      AA = -1.0D0
      BB = 1.0D0
C
C     SET TMID SO THAT IT IS NOT IN THE EXACT MIDDLE.
C
      TMID = 0.1D0*PI
C
      INDD = 1
      LAMBDA = ALFLAM
      EIGSAV = LAMBDA
      EPS = 1.D-7
C
      if (IND.ge.2) go to 10
C
   50 continue
C
C     FOR U:
C
      THA = HPI
      Y(1) = THA
      Y(2) = 1.0D0
      Y(3) = 0.0D0
      DTHDAA = 0.0D0
      DTHDEA = 1.0D0
      LFLAG = 1
      call INTEG(AA,THA,DTHDAA,DTHDEA,TMID,DUM,DUM,EPS,Y,ERR,LCIRCA,AOK,
     +           SINGA,OSCA,LFLAG)
      if (LFLAG.eq.5) then
         EPS = 10.0D0*EPS
         go to 50

      end if

      RHOL = EXP(Y(3))
      THL = Y(1)
      UL = RHOL*SIN(THL)
C
      THB = HPI
      Y(1) = THB
      Y(2) = 1.0D0
      Y(3) = 0.0D0
      DTHDBB = 0.0D0
      DTHDEB = 1.0D0
      LFLAG = 1
      call INTEG(BB,THB,DTHDBB,DTHDEB,TMID,DUM,DUM,EPS,Y,ERR,LCIRCB,BOK,
     +           SINGB,OSCB,LFLAG)
      if (LFLAG.eq.5) then
         EPS = 10.0D0*EPS
         go to 50

      end if

      RHOR = EXP(Y(3))
      THR = Y(1)
      UR = RHOR*SIN(THR)
C
C     FOR V:
C
      THA = 0.0D0
      Y(1) = THA
      Y(2) = 1.0D0
      Y(3) = 0.0D0
      DTHDAA = 0.0D0
      DTHDEA = 1.0D0
      LFLAG = 1
      call INTEG(AA,THA,DTHDAA,DTHDEA,TMID,DUM,DUM,EPS,Y,ERR,LCIRCA,AOK,
     +           SINGA,OSCA,LFLAG)
      if (LFLAG.eq.5) then
         EPS = 10.0D0*EPS
         go to 50

      end if

      RHOL = EXP(Y(3))
      THL = Y(1)
      VL = RHOL*SIN(THL)
      PVPL = Z*RHOL*COS(THL)
C
      THB = 0.0D0
      Y(1) = THB
      Y(2) = 1.0D0
      Y(3) = 0.0D0
      DTHDBB = 0.0D0
      DTHDEB = 1.0D0
      LFLAG = 1
      call INTEG(BB,THB,DTHDBB,DTHDEB,TMID,DUM,DUM,EPS,Y,ERR,LCIRCB,BOK,
     +           SINGB,OSCB,LFLAG)
      if (LFLAG.eq.5) then
         EPS = 10.0D0*EPS
         go to 50

      end if

      RHOR = EXP(Y(3))
      THR = Y(1)
      VR = RHOR*SIN(THR)
      PVPR = Z*RHOR*COS(THR)
      FF = (VR* (UL*PVPL-1.0D0)-CC*VL* (UR*PVPR-1.0D0))* (VL-VR/
     +     CC) - VL*VR* (UL-CC*UR)* (PVPL-PVPR/CC)
      return
C
   10 continue
C
C     FOR U:
C
      THA = HPI
      Y(1) = THA
      Y(2) = 1.0D0
      Y(3) = 0.0D0
      DTHDAA = 0.0D0
      DTHDEA = 1.0D0
      LFLAG = 1
      call INTEG(AA,THA,DTHDAA,DTHDEA,BB,DUM,DUM,EPS,Y,ERR,LCIRCA,AOK,
     +           SINGA,OSCA,LFLAG)
      if (LFLAG.eq.5) then
         EPS = 10.0D0*EPS
         go to 10

      end if

      RHOB = EXP(Y(3))
      THB = Y(1)
      THETU = THB
      UB = RHOB*SIN(THB)
C
C     FOR V:
C
      THA = 0.0D0
      Y(1) = THA
      Y(2) = 1.0D0
      Y(3) = 0.0D0
      DTHDAA = 0.0D0
      DTHDEA = 1.0D0
      LFLAG = 1
      call INTEG(AA,THA,DTHDAA,DTHDEA,BB,DUM,DUM,EPS,Y,ERR,LCIRCA,AOK,
     +           SINGA,OSCA,LFLAG)
      if (LFLAG.eq.5) then
         EPS = 10.0D0*EPS
         go to 10

      end if

      RHOB = EXP(Y(3))
      THB = Y(1)
      THETV = THB
      VB = RHOB*SIN(THB)
      PVPB = Z*RHOB*COS(THB)
      FF = CC*PVPB + UB/CC - 2.0D0
      return

      end function FF
      subroutine FIT(TH1,TH,TH2)
C     **********
C
C     This program converts TH into an 'equivalent' angle between
C     TH1 and TH2.  We assume TH1.LT.TH2 and PI.LE.(TH2-TH1).
C
C     **********
C     .. Scalars in Common ..
      double precision HPI,PI,TWOPI
C     ..
C     .. Intrinsic Functions ..
      intrinsic AINT
C     ..
C     .. Common blocks ..
      common /PIE/PI,TWOPI,HPI
C     ..
C     .. Scalar Arguments ..
      double precision TH,TH1,TH2
C     ..
      if (TH.lt.TH1) TH = TH + AINT((TH1-TH+PI)/PI)*PI
      if (TH.gt.TH2) TH = TH - AINT((TH-TH2+PI)/PI)*PI
      return

      end subroutine FIT
      subroutine FZ(UU,Y,YP)
C     **********
C
C     This subroutine evaluates the derivative of the function PHI
C     in the regularization at a singular endpoint.
C
C     Here, HU means -(PUP)' + QU.
C
C     **********
C     .. Scalars in Common ..
      double precision EIG
      integer IND
C     ..
C     .. Local Scalars ..
      double precision A1122,A12,A21,AU,AV,B1122,B12,B21,C,C2,D,DT,HU,
     +                 HV,PHI,PUP,PVP,S,S2,SC,T,U,V,WW,WX,X
C     ..
C     .. External Functions ..
      double precision W
      external W
C     ..
C     .. External Subroutines ..
cc      external DXDT,UV
C     ..
C     .. Intrinsic Functions ..
      intrinsic COS,MOD,SIN
C     ..
C     .. Common blocks ..
      common /DATAF/EIG,IND
C     ..
C     .. Scalar Arguments ..
      double precision UU
C     ..
C     .. Array Arguments ..
      double precision Y(2),YP(3)
C     ..
      if (MOD(IND,2).eq.1) then
         T = UU
         PHI = Y(1)
      else
         T = Y(1)
         PHI = UU
      end if

      call DXDT(T,DT,X)
      call UV(X,U,PUP,V,PVP,HU,HV)
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
      if (IND.eq.1) then
         WW = 2.0D0* (A12+A21)*SC + A1122* (C2-S2)
         YP(2) = -DT* (WW*Y(2)+2.0D0*B1122*SC+B12*S2-B21*C2)/D
         YP(3) = 0.5D0*DT*WW/D
      else if (IND.eq.2) then
         YP(2) = YP(2)/YP(1)
         YP(3) = YP(3)/YP(1)
         YP(1) = 1.0D0/YP(1)
      else if (IND.eq.3) then
      else
         YP(1) = 1.0D0/YP(1)
      end if

      return

      end subroutine FZ
      subroutine GERKZ(F,NEQ,Y,TIN,TOUT,REPS,AEPS,LFLAG,ER,WORK,IWORK)
C     **********
C     **********
C     .. Scalars in Common ..
      double precision EPSMIN,Z
C     ..
C     .. Local Scalars ..
      double precision T,TOUTS
      integer I,J,K,L,LLFLAG
C     ..
C     .. Arrays in Common ..
      double precision TEE(100),ZEE(100)
      integer JAY(100)
C     ..
C     .. Local Arrays ..
      double precision U(3)
C     ..
C     .. External Subroutines ..
cc      external GERK,THTOTHZ,THZTOTH
C     ..
C     .. Intrinsic Functions
C     .. Common blocks ..
      common /RNDOFF/EPSMIN
      common /TEEZ/TEE
      common /ZEE/Z
      common /ZEEZ/JAY,ZEE
C     ..
C     .. Scalar Arguments ..
      double precision AEPS,REPS,TIN,TOUT
      integer LFLAG,NEQ
C     ..
C     .. Array Arguments ..
      double precision ER(3),WORK(27),Y(3)
      integer IWORK(5)
C     ..
C     .. Subroutine Arguments ..
      external F
C     ..
C     .. Intrinsic Functions ..
      intrinsic MAX,MIN
C     ..
      T = TIN
      if (TIN.lt.TOUT) then
         do 10 I = 1,19
            if (TEE(I)-EPSMIN.le.TIN .and. TIN.lt.TEE(I+1)+EPSMIN) J = I
            if (TEE(I)-EPSMIN.lt.TOUT .and.
     +          TOUT.le.TEE(I+1)+EPSMIN) L = I
   10    continue
         do 30 K = J,L
            TOUTS = MIN(TOUT,TEE(K+1))
            Z = ZEE(K+1)
            if (Z.eq.0.0D0) Z = 1.0D0
            call THTOTHZ(Y,Z,U)
            LLFLAG = 1
   20       continue
            call GERK(F,NEQ,U,T,TOUTS,REPS,AEPS,LLFLAG,ER,WORK,IWORK)
            if (LLFLAG.gt.3) then
               WRITE (021,FMT=*) ' llflag = ',LLFLAG
               LFLAG = 5
               return

            end if

            if (LLFLAG.eq.3 .or. LLFLAG.eq.-2) go to 20
            call THZTOTH(U,Z,Y)
   30    continue
      else
         do 40 I = 20,2,-1
            if (TEE(I-1)-EPSMIN.lt.TIN .and. TIN.le.TEE(I)+EPSMIN) J = I
            if (TEE(I-1)-EPSMIN.le.TOUT .and.
     +          TOUT.lt.TEE(I)+EPSMIN) L = I
   40    continue
         do 60 K = J,L,-1
            TOUTS = MAX(TOUT,TEE(K-1))
            Z = ZEE(K)
            if (Z.eq.0.0D0) Z = 1.0D0
            call THTOTHZ(Y,Z,U)
            LLFLAG = 1
   50       continue
            call GERK(F,NEQ,U,T,TOUTS,REPS,AEPS,LLFLAG,ER,WORK,IWORK)
            if (LLFLAG.gt.3) then
               WRITE (021,FMT=*) ' llflag = ',LLFLAG
               LFLAG = 5
               return

            end if

            if (LLFLAG.eq.3 .or. LLFLAG.eq.-2) go to 50
            call THZTOTH(U,Z,Y)
   60    continue
      end if

      TIN = T
      LFLAG = LLFLAG
      return

      end subroutine GERKZ
      subroutine INTEG(TEND,THEND,DTHDAA,DTHDE,TMID,COEF1,COEF2,EPS,Y,
     +                 ER,LCIRC,OK,SING,OSC,IFLAG)
C     **********
C     **********
C     .. Scalars in Common ..
      double precision EIG,HPI,PI,TSAVEL,TSAVER,TWOPI,Z
      integer IND
C     ..
C     .. Arrays in Common ..
      double precision TT(7,2),YY(7,3,2)
      integer NT(2)
C     ..
C     .. Local Scalars ..
      double precision C,D,DDD,DPHIDAA,DPHIDE,DTHIN,DUM,EFF,FAC2,HU,HV,
     +                 P1,PHI,PHI0,PUP,PVP,PYPZ,PYPZ0,Q1,RHOSQ,S,T,TH,
     +                 TH0,THBAR,THETA,THIN,THU,THU0,THV,THV0,TIN,
     +                 TINTHZ,TMP,TOUT,TSTAR,U,V,W1,XSTAR,XT,XT0,YSTAR,
     +                 YZ,YZ0,ZSAV
      integer I,J,K2PI,KFLAG,LFLAG,M
      logical LOGIC
C     ..
C     .. Local Arrays ..
      double precision ERZ(3),WORK(27),YP(3),YU(3)
      integer IWORK(5)
C     ..
C     .. External Subroutines ..
cc      external DXDT,F,FZ,GERK,GERKZ,INTEGT,SETTHU,UV,UVPHI,WR
cc?      external F,FZ
C     ..
C     .. External Functions ..
      double precision P,Q,W
      external P,Q,W
C     ..
C     .. Intrinsic Functions ..
      intrinsic ABS,ATAN2,COS,EXP,LOG,MIN,SIGN,SIN
C     ..
C     .. Common blocks ..
      common /DATAF/EIG,IND
      common /PIE/PI,TWOPI,HPI
      common /TEMP/TT,YY,NT
      common /TSAVE/TSAVEL,TSAVER
      common /ZEE/Z
C     ..
C
C     Note: The input values of THEND and DTHDAA are overwritten
C           when integrating from a limit circle endpoint.
C
C     .. Scalar Arguments ..
      double precision COEF1,COEF2,DTHDAA,DTHDE,EPS,TEND,THEND,TMID
      integer IFLAG
      logical LCIRC,OK,OSC,SING
C     ..
C     .. Array Arguments ..
      double precision ER(3),Y(3)
C     ..
      IFLAG = 1
      IND = 1
      if (OSC) then
         if (.not.OK) then
            LOGIC = .FALSE.
         else if (TEND.le.TMID) then
            LOGIC = TEND .ge. TSAVEL
         else
            LOGIC = TEND .le. TSAVER
         end if
C
         if (LOGIC) then
            TINTHZ = TEND
            TH = THEND
            DTHIN = DTHDE
            Y(1) = TH
            Y(2) = DTHIN
            EFF = Y(3)
         else
C           DO (INTEGRATE-FOR-PHI-OSC)
            ZSAV = Z
            Z = 1.0D0
            T = TEND
            PHI0 = ATAN2(COEF2,COEF1)
            if (TMID.gt.TEND) then
               J = 1
            else
               J = 2
            end if
C
C     We want -PI/2 .lt. PHI0 .le. PI/2.
C
            if (COEF1.lt.0.0D0) PHI0 = PHI0 - SIGN(PI,COEF2)
            Y(1) = PHI0
            Y(2) = 0.0D0
            Y(3) = 0.0D0
            call DXDT(T,TMP,XT0)
            call FZ(T,Y,YP)
            DPHIDAA = -YP(1)
            call UV(XT0,U,PUP,V,PVP,HU,HV)
            D = U*PVP - V*PUP
C
C     Set THU0 and THV0.
C
            THU0 = ATAN2(U,PUP)
            if (U.lt.0.0D0) THU0 = THU0 + TWOPI
            call SETTHU(XT0,THU0)
            THV0 = ATAN2(V,PVP)
            if (V.lt.0.0D0) THV0 = THV0 + TWOPI
   10       continue
            if (THV0.lt.THU0) then
               THV0 = THV0 + TWOPI
               go to 10

            end if
C
C     Set TH0 and copy into THEND, overwriting its input value.
C     Also, redefine DTHDAA, overwriting its input value.
C
            call UVPHI(U,PUP,V,PVP,THU0,THV0,PHI0,TH0)
            THEND = TH0
            C = COS(PHI0)
            S = SIN(PHI0)
            YZ0 = U*C + V*S
            PYPZ0 = PUP*C + PVP*S
            DUM = ABS(COS(TH0))
            if (DUM.ge.0.5D0) then
               FAC2 = -D* (DUM/PYPZ0)**2
            else
               FAC2 = -D* (SIN(TH0)/YZ0)**2
            end if

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
            if (EIG.eq.0.0D0) then
               if (T.le.TSAVEL) TOUT = TSAVEL
               if (T.gt.TSAVER) TOUT = TSAVER
               KFLAG = 1
            end if

   20       continue
            call GERK(FZ,3,Y,T,TOUT,EPS,EPS,KFLAG,ER,WORK,IWORK)
            if (KFLAG.gt.3) then
               WRITE (021,FMT=*) ' KFLAG1 = 5 '
               IFLAG = 5
               return

            end if

            if (KFLAG.eq.3) go to 20
            if (Y(3).lt.-15.0D0) Y(3) = -15.0D0
            PHI = Y(1)
C
C     Store up to seven values of (T,PHI) for later reference.
C
            I = I + 1
            if (I.le.7) then
               TT(I,J) = T
               YY(I,1,J) = PHI
               YY(I,2,J) = Y(2)
               YY(I,3,J) = Y(3)
            end if

            if (10.0D0*ABS(PHI-PHI0).ge.PI) go to 30
            if (KFLAG.eq.-2) go to 20
   30       continue
            if (T.le.TOUT) then
               TSAVEL = T
            else
               TSAVER = T
            end if

            NT(J) = I
            DPHIDE = Y(2)
            TINTHZ = T
            call DXDT(T,TMP,XT)
            call UV(XT,U,PUP,V,PVP,HU,HV)
            D = U*PVP - V*PUP
C
C     Set THU and THV.
C
            THU = ATAN2(U,PUP)
            if (U.lt.0.0D0) THU = THU + TWOPI
            call SETTHU(XT,THU)
            THV = ATAN2(V,PVP)
            if (V.lt.0.0D0) THV = THV + TWOPI
   40       continue
            if (THV.lt.THU) then
               THV = THV + TWOPI
               go to 40

            end if
C
C     Now define TH in terms of PHI, THU, and THV.
C
            call UVPHI(U,PUP,V,PVP,THU,THV,PHI,TH)
            YZ = U*COS(PHI) + V*SIN(PHI)
            PYPZ = PUP*COS(PHI) + PVP*SIN(PHI)
            Y(1) = TH
            S = SIN(TH)
            C = COS(TH)
            DUM = ABS(C)
            if (DUM.ge.0.5D0) then
               FAC2 = -D* (DUM/PYPZ)**2
            else
               FAC2 = -D* (S/YZ)**2
            end if

            DTHIN = FAC2*DPHIDE
            Y(2) = DTHIN
            Z = ZSAV
            RHOSQ = EXP(2.0D0*Y(3))* (YZ**2+PYPZ**2)
            EFF = 0.5D0*LOG(RHOSQ* (S**2+ (C/Z)**2))
            Y(3) = EFF
C              END (INTEGRATE-FOR-PHI-OSC)
         end if

         if (TINTHZ.ne.TMID) then
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
            Y(3) = 0.0D0
   50       continue
            call GERKZ(F,3,Y,T,TOUT,EPS,EPS,LFLAG,ERZ,WORK,IWORK)
            if (LFLAG.gt.3) then
               WRITE (021,FMT=*) ' LFLAGZ = ',LFLAG
               IFLAG = 5
               return

            end if

            if (LFLAG.eq.3 .or. LFLAG.eq.-2) go to 50
C              END (INTEGRATE-FOR-THETAZ)
            Y(3) = Y(3) + EFF
            Z = 1.0D0
         end if

      else if (LCIRC) then
C        DO (INTEGRATE-FOR-PHI-NONOSC)
         T = TEND
         PHI0 = ATAN2(COEF2,COEF1)
C
C     We want -PI/2 .lt. PHI0 .le. PI/2.
C
         if (COEF1.lt.0.0D0) PHI0 = PHI0 - SIGN(PI,COEF2)
         Y(1) = PHI0
         Y(2) = 0.0D0
         Y(3) = 0.0D0
         call DXDT(T,TMP,XT)
         call FZ(T,Y,YP)
         DPHIDAA = -YP(1)
         call UV(XT,U,PUP,V,PVP,HU,HV)
         D = U*PVP - V*PUP
         YZ0 = U*COEF1 + V*COEF2
         PYPZ0 = (PUP*COEF1+PVP*COEF2)/Z
C
C     Set TH0 and copy into THEND, overwriting its input value.
C     Also, redefine DTHDAA, overwriting its input value.
C
         TH0 = ATAN2(YZ0,PYPZ0)
         if (YZ0.lt.0.0D0) TH0 = TH0 + TWOPI
         THEND = TH0
         DUM = ABS(COS(TH0))
         if (DUM.ge.0.5D0) then
            FAC2 = (-D/Z)* (DUM/PYPZ0)**2
         else
            FAC2 = (-D/Z)* (SIN(TH0)/YZ0)**2
         end if

         DTHDAA = FAC2*DPHIDAA
C
C     In the next piece, we assume TH0 .ge. 0.
C
         M = 0
         if (TH0.eq.0.0D0) M = -1
         if (TH0.gt.PI .or. (TH0.eq.PI.and.T.lt.TMID)) M = 1
         PHI = PHI0
         K2PI = 0
         YZ0 = U*COS(PHI0) + V*SIN(PHI0)
         if (TMID.gt.TEND) then
            J = 1
            TSTAR = -0.9999D0
         else
            J = 2
            TSTAR = 0.9999D0
         end if

         DPHIDE = 0.0D0
         call WR(FZ,EPS,TSTAR,PHI0,PHI0,DPHIDE,TOUT,Y,TT(1,J),YY(1,1,J),
     +           ERZ,WORK,IWORK)
         T = TOUT
         DDD = MIN(0.01D0,ABS(TMID-TOUT))
         TOUT = TOUT + DDD* (TMID-TOUT)/ABS(TMID-TOUT)
         KFLAG = 1
         call GERK(FZ,3,Y,T,TOUT,EPS,EPS,KFLAG,ER,WORK,IWORK)
         TT(7,J) = T
         YY(7,1,J) = Y(1)
         YY(7,2,J) = Y(2)
         YY(7,3,J) = Y(3)
   60    continue
         T = TOUT
         call DXDT(T,TMP,XT)
         call UV(XT,U,PUP,V,PVP,HU,HV)
         PHI = Y(1)
         S = SIN(PHI)
         C = COS(PHI)
         YZ = U*C + V*S
         if (YZ*YZ0.lt.0.0D0) K2PI = K2PI + 1
         YZ0 = YZ
         if (KFLAG.gt.3) then
            WRITE (021,FMT=*) ' KFLAG2 = ',KFLAG
            IFLAG = 5
            return

         end if

         if (KFLAG.eq.3 .or. KFLAG.eq.-2) go to 60
C
C     Convert from PHI to THETA.
C
         DPHIDE = Y(2)
         D = U*PVP - V*PUP
         PYPZ = (PUP*C+PVP*S)/Z
         THBAR = ATAN2(YZ,PYPZ)
         if (TMID.gt.TEND .and. THBAR.lt.TH0 .and.
     +       PHI.lt.PHI0) THBAR = THBAR + TWOPI
         if (TMID.lt.TEND .and. THBAR.gt.TH0 .and.
     +       PHI.gt.PHI0) THBAR = THBAR - TWOPI
         TH = THBAR - M*PI
         if (TH.lt.-PI) TH = TH + TWOPI
         if (TH.gt.TWOPI) TH = TH - TWOPI
         if (TMID.lt.TEND .and. K2PI.gt.1) TH = TH - (K2PI-1)*TWOPI
         if (TMID.gt.TEND .and. K2PI.gt.1) TH = TH + (K2PI-1)*TWOPI
         if (TMID.gt.TEND .and. TH*TH0.lt.0.0D0) TH = TH + TWOPI
C
C     We now have YZ, PYPZ, PHI and TH.
C
         DUM = ABS(COS(TH))
         if (DUM.ge.0.5D0) then
            FAC2 = - (D/Z)* (DUM/PYPZ)**2
         else
            FAC2 = - (D/Z)* (SIN(TH)/YZ)**2
         end if

         DTHIN = FAC2*DPHIDE
         THETA = ATAN2(YZ,Z*PYPZ)
         S = SIN(THETA)
         C = COS(THETA)
         RHOSQ = EXP(2.0D0*Y(3))* (YZ**2+ (Z*PYPZ)**2)
         EFF = 0.5D0*LOG(RHOSQ* (S**2+ (C/Z)**2))
C           END (INTEGRATE-FOR-PHI-NONOSC)
         TIN = TOUT
         TOUT = TMID
         THIN = TH
C        DO (INTEGRATE-FOR-THETA)
         call INTEGT(TIN,TOUT,THIN,DTHIN,EPS,Y,IFLAG,ER,WORK,IWORK)
C           END (INTEGRATE-FOR-THETA)
         Y(3) = Y(3) + EFF
      else if (.not.OK .and. .not.SING) then
C
C     This is the 'weakly regular' case.
C
         if (TMID.gt.TEND) then
            J = 1
            TSTAR = -0.9999D0
         else
            J = 2
            TSTAR = 0.9999D0
         end if

         call DXDT(TSTAR,TMP,XSTAR)
         P1 = 1.0D0/P(XSTAR)
         Q1 = Q(XSTAR)
         W1 = W(XSTAR)
         YSTAR = THEND + 0.5D0* (TSTAR-TEND)*
     +           (P1*COS(THEND)**2+ (EIG*W1-Q1)*SIN(THEND)**2)
         call WR(F,EPS,TSTAR,YSTAR,THEND,DTHDE,TOUT,YU,TT(1,J),
     +           YY(1,1,J),ERZ,WORK,IWORK)
         T = TOUT
         DDD = MIN(0.01D0,ABS(TMID-TOUT))
         TOUT = TOUT + DDD* (TMID-TOUT)/ABS(TMID-TOUT)
         KFLAG = 1
         call GERK(F,3,YU,T,TOUT,EPS,EPS,KFLAG,ER,WORK,IWORK)
         TT(7,J) = T
         YY(7,1,J) = YU(1)
         YY(7,2,J) = YU(2)
         YY(7,3,J) = YU(3)
         TIN = TOUT
         TOUT = TMID
         THIN = YU(1)
         DTHIN = YU(2)
C        DO (INTEGRATE-FOR-THETA)
         call INTEGT(TIN,TOUT,THIN,DTHIN,EPS,Y,IFLAG,ER,WORK,IWORK)
C           END (INTEGRATE-FOR-THETA)
      else
C
C     This is the regular (not weakly regular) or limit point case.
C
         TIN = TEND
         TOUT = TMID
         THIN = THEND
         DTHIN = DTHDE
C        DO (INTEGRATE-FOR-THETA)
         call INTEGT(TIN,TOUT,THIN,DTHIN,EPS,Y,IFLAG,ER,WORK,IWORK)
C           END (INTEGRATE-FOR-THETA)
      end if

      return

      end subroutine INTEG
      subroutine INTEGT(TIN,TOUT,THIN,DTHIN,EPS,Y,IFLAG,ER,WORK,IWORK)
C     **********
C
C     Integrate for theta.
C
C     **********
C     .. Local Scalars ..
      double precision T
      integer LFLAG
C     ..
C     .. External Subroutines ..
cc      external F,GERK
C     ..
C     DO (INTEGRATE-FOR-TH)
C     .. Scalar Arguments ..
      double precision DTHIN,EPS,THIN,TIN,TOUT
      integer IFLAG
C     ..
C     .. Array Arguments ..
      double precision ER(3),WORK(27),Y(3)
      integer IWORK(5)
C     ..
      T = TIN
      Y(1) = THIN
      Y(2) = DTHIN
      Y(3) = 0.0D0
      LFLAG = 1
   10 continue
      call GERK(F,3,Y,T,TOUT,EPS,EPS,LFLAG,ER,WORK,IWORK)
      if (LFLAG.eq.3) go to 10
      if (LFLAG.gt.3) then
         WRITE (021,FMT=*) ' LFLAG = ',LFLAG
         IFLAG = 5
         return

      end if

      THIN = Y(1)
      DTHIN = Y(2)
C        END (INTEGRATE-FOR-TH)
      return

      end subroutine INTEGT
      subroutine INTPOL(N,XN,FN,X,ABSERR,MAXDEG,ANS,ERROR)
C     **********
C
C     This subroutine forms an interpolating polynomial for data pairs.
C     It is called from EXTRAP in solving a Sturm-Liouville problem.
C
C     **********
C     .. Local Scalars ..
      double precision PROD
      integer I,I1,II,IJ,IK,IKM1,J,K,L,LIMIT
C     ..
C     .. Local Arrays ..
      double precision V(10,10)
      integer INDEX(10)
C     ..
C     .. Intrinsic Functions ..
      intrinsic ABS,MIN
C     ..
C     .. Scalar Arguments ..
      double precision ABSERR,ANS,ERROR,X
      integer MAXDEG,N
C     ..
C     .. Array Arguments ..
      double precision FN(N),XN(N)
C     ..
      L = MIN(MAXDEG,N-2) + 2
      LIMIT = MIN(L,N-1)
      do 10 I = 1,N
         V(I,1) = ABS(XN(I)-X)
         INDEX(I) = I
   10 continue
      do 30 I = 1,LIMIT
         do 20 J = I + 1,N
            II = INDEX(I)
            IJ = INDEX(J)
            if (V(II,1).gt.V(IJ,1)) then
               INDEX(I) = IJ
               INDEX(J) = II
            end if

   20    continue
   30 continue
      PROD = 1.0D0
      I1 = INDEX(1)
      ANS = FN(I1)
      V(1,1) = FN(I1)
      do 50 K = 2,L
         IK = INDEX(K)
         V(K,1) = FN(IK)
         do 40 I = 1,K - 1
            II = INDEX(I)
            V(K,I+1) = (V(I,I)-V(K,I))/ (XN(II)-XN(IK))
   40    continue
         IKM1 = INDEX(K-1)
         PROD = (X-XN(IKM1))*PROD
         ERROR = PROD*V(K,K)
         if (ABS(ERROR).le.ABSERR) return
         ANS = ANS + ERROR
   50 continue
      ANS = ANS - ERROR
      return

      end subroutine INTPOL
      subroutine LPOL(KEIGS,XN,FN,X,F)
C     **********
C     **********
C     .. Local variables ..
C     .. Scalar Arguments ..
      double precision F,X
      integer KEIGS
C     ..
C     .. Array Arguments ..
      double precision FN(KEIGS),XN(KEIGS)
C     ..
C     .. Local Scalars ..
      double precision X12,X13,X14,X15,X21,X23,X24,X25,X31,X32,X34,X35,
     +                 X41,X42,X43,X45,X51,X52,X53,X54,XX1,XX2,XX3,XX4,
     +                 XX5
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
      if (KEIGS.ge.4) then
         X14 = XN(1) - XN(4)
         X24 = XN(2) - XN(4)
         X34 = XN(3) - XN(4)
         X41 = -X14
         X42 = -X24
         X43 = -X34
         XX4 = X - XN(4)
      end if

      if (KEIGS.eq.5) then
         X15 = XN(1) - XN(5)
         X25 = XN(2) - XN(5)
         X35 = XN(3) - XN(5)
         X45 = XN(4) - XN(5)
         X51 = -X15
         X52 = -X25
         X53 = -X35
         X54 = -X45
         XX5 = X - XN(5)
      end if

      if (KEIGS.eq.3) then
         F = XX2*XX3*FN(1)/ (X12*X13) + XX1*XX3*FN(2)/
     +       (X21*X23) + XX1*XX2*FN(3)/ (X31*X32)
      else if (KEIGS.eq.4) then
         F = XX2*XX3*XX4*FN(1)/ (X12*X13*X14) + XX1*XX3*XX4*FN(2)/
     +       (X21*X23*X24) + XX1*XX2*XX4*FN(3)/
     +       (X31*X32*X34) + XX1*XX2*XX3*FN(4)/ (X41*X42*X43)
      else
         F = XX2*XX3*XX4*XX5*FN(1)/ (X12*X13*X14*X15) +
     +       XX1*XX3*XX4*XX5*FN(2)/ (X21*X23*X24*X25) +
     +       XX1*XX2*XX4*XX5*FN(3)/ (X31*X32*X34*X35) +
     +       XX1*XX2*XX3*XX5*FN(4)/ (X41*X42*X43*X45) +
     +       XX1*XX2*XX3*XX4*FN(5)/ (X51*X52*X53*X54)
      end if

      return

      end subroutine LPOL
C
C
      subroutine PERIO(A,B,INTAB,P0ATA,QFATA,P0ATB,QFATB,A1,A2,B1,B2,
     +                 NUMEIG,EIG,TOL,IFLAG,SLFN,SINGATA,SINGATB,CIRCLA,
     +                 CIRCLB,OSCILA,OSCILB)
C     **********
C     **********
C     .. Scalars in Common ..
      double precision CC,UB,UL,UR,VB,VL,VR,Z
C      double precision EPSMIN
      integer IND
C     ..
C     .. Local Scalars ..
      double precision A1D,A1N,A2D,A2N,AE,B1D,B1N,B2D,B2N,EIGLO,EIGUP,
     +                 LAMBDA,LAMUP,RE,TOLL,TOLS
      integer JFLAG
C     ..
C     .. External Subroutines ..
cc      external FZERO,SLEIGN2
C     ..
C     .. External Functions ..
cc?      external FF
cc?      double precision FF
C     ..
C     .. Intrinsic Functions ..
      intrinsic ABS,MAX
C     ..
C     .. Common blocks ..
      common /EPP2/CC,UL,UR,VL,VR,UB,VB,IND
C      common /RNDOFF/EPSMIN
      common /ZEE/Z
C     ..
C
C     Get upper and lower bounds as accurately as possible.
C
C     .. Scalar Arguments ..
      double precision A,A1,A2,B,B1,B2,CIRCLA,CIRCLB,EIG,OSCILA,OSCILB,
     +                 P0ATA,P0ATB,QFATA,QFATB,SINGATA,SINGATB,TOL
      integer IFLAG,INTAB,NUMEIG
C     ..
C     .. Array Arguments ..
      double precision SLFN(9)
C     ..
      TOLL = 0.0D0
      TOLS = TOL
C
      A1N = 0.0D0
      A2N = 1.0D0
      B1N = 0.0D0
      B2N = 1.0D0
      TOL = TOLL
      EIG = 0.0D0
      IFLAG = 0
      call SLEIGN2(A,B,INTAB,P0ATA,QFATA,P0ATB,QFATB,A1N,A2N,B1N,B2N,
     +             NUMEIG,EIG,TOL,IFLAG,0,SLFN,SINGATA,SINGATB,CIRCLA,
     +             CIRCLB,OSCILA,OSCILB)
      WRITE (021,FMT=*) ' eiglo,iflag = ',EIG,IFLAG
      EIGLO = EIG - 0.1D0
C
      A1D = 1.0D0
      A2D = 0.0D0
      B1D = 1.0D0
      B2D = 0.0D0
      TOL = TOLL
      EIG = 0.0D0
      IFLAG = 0
      call SLEIGN2(A,B,INTAB,P0ATA,QFATA,P0ATB,QFATB,A1D,A2D,B1D,B2D,
     +             NUMEIG,EIG,TOL,IFLAG,0,SLFN,SINGATA,SINGATB,CIRCLA,
     +             CIRCLB,OSCILA,OSCILB)
      WRITE (021,FMT=*) ' eigup,iflag = ',EIG,IFLAG
      EIGUP = EIG + 0.1D0
      LAMBDA = 0.5D0* (EIGLO+EIGUP)
C
C     The following call to SLEIGN2 sets the stage for INTEG.
C
      TOL = .001D0
      IFLAG = -1
      call SLEIGN2(A,B,INTAB,P0ATA,QFATA,P0ATB,QFATB,A1,A2,B1,B2,NUMEIG,
     +             LAMBDA,TOL,IFLAG,0,SLFN,SINGATA,SINGATB,CIRCLA,
     +             CIRCLB,OSCILA,OSCILB)
      Z = 1.0D0
C
      EIG = EIGLO
      LAMUP = EIGUP
      RE = TOLS
      AE = RE
C
C     IND = 2 selects the classical function d(LAMBDA) in function FF.
C
      IND = 2
      call FZERO(FF,EIG,LAMUP,LAMBDA,RE,AE,JFLAG)
      if (JFLAG.ne.1) WRITE (021,FMT=*) ' jflag = ',JFLAG
      IFLAG = JFLAG
      TOL = ABS(EIG-LAMUP)/MAX(1.0D0,ABS(EIG))
      WRITE (021,FMT=*) ' EIGLO,EIGUP = ',EIGLO,EIGUP
C
      return

      end subroutine PERIO
      subroutine SETMID(MF,ML,EIG,QS,WS,IMID,TMID)
C     **********
C
C     This procedures tests the interval sample points in the order
C     50,51,49,52,48,...,etc. for the first one where the expression
C     (lambda*w-q) is positive.  This point is designated TMID.
C
C     **********
C     .. Local Scalars ..
      double precision S
      integer I,J
C     ..
C     .. External Functions ..
cc      double precision TFROMI
cc      external TFROMI
C     ..
C     .. Scalar Arguments ..
      double precision EIG,TMID
      integer IMID,MF,ML
C     ..
C     .. Array Arguments ..
      double precision QS(*),WS(*)
C     ..
      S = -1.0D0
      do 10 J = 1,100
         I = 50 + S* (J/2)
         S = -S
         if (I.lt.MF .or. I.gt.ML) go to 20
         if (EIG*WS(I)-QS(I).gt.0.0D0) then
            IMID = I
            TMID = TFROMI(IMID)
            go to 20

         end if

   10 continue
   20 continue
      WRITE (021,FMT=*) ' new tmid = ',TMID
      return

      end subroutine SETMID
      subroutine SETTHU(X,THU)
C     **********
c
c  This subroutine establishes a definite value for THU,
c    the phase angle for the function U, including an
c    appropriate integer multiple of pi
c  It needs the numbers MMW(*) found in THUM
c
C     **********
C     .. Scalars in Common ..
      double precision HPI,PI,TWOPI
      integer MMWD
C     ..
C     .. Arrays in Common ..
      double precision YS(197)
      integer MMW(98)
C     ..
C     .. Local Scalars ..
      integer I
C     ..
C     .. Common blocks ..
      common /PASS/YS,MMW,MMWD
      common /PIE/PI,TWOPI,HPI
C     ..
C     .. Scalar Arguments ..
      double precision THU,X
C     ..
      do 10 I = 1,MMWD
         if (X.ge.YS(MMW(I)) .and. X.le.YS(MMW(I)+1)) then
            if (THU.gt.PI) then
               THU = THU + (I-1)*TWOPI
               return

            else
               THU = THU + I*TWOPI
               return

            end if

         end if

   10 continue
      do 20 I = 1,MMWD
         if (X.ge.YS(MMW(I))) THU = THU + TWOPI
   20 continue
      return

      end subroutine SETTHU
      double precision function TFROMI(I)
C     **********
C
C     This function associates the value of an interval sample point
C     with its index.
C
C     **********
C     .. Scalar Arguments ..
      integer I
C     ..
      if (I.lt.8) then
         TFROMI = -1.0D0 + 0.1D0/4.0D0** (8-I)
      else if (I.gt.92) then
         TFROMI = 1.0D0 - 0.1D0/4.0D0** (I-92)
      else
         TFROMI = 0.0227D0* (I-50)
      end if

      return

      end function TFROMI
      subroutine THTOTHZ(Y,Z,U)
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
      double precision HPI,PI,TWOPI
C     ..
C     .. Local Scalars ..
      double precision DTH,DTHZ,DUM,FAC,PIK,REMTH,TH,THZ
      integer K
C     ..
C     .. Intrinsic Functions ..
      intrinsic ABS,ATAN,COS,LOG,SIN,TAN
C     ..
C     .. Common blocks ..
      common /PIE/PI,TWOPI,HPI
C     ..
C     .. Scalar Arguments ..
      double precision Z
C     ..
C     .. Array Arguments ..
      double precision U(3),Y(3)
C     ..
      TH = Y(1)
      DTH = Y(2)
      K = TH/PI
      if (TH.lt.0.0D0) K = K - 1
      PIK = K*PI
      REMTH = TH - PIK
      if (4.0D0*REMTH.le.PI) then
         THZ = ATAN(Z*TAN(REMTH)) + PIK
      else if (4.0D0*REMTH.ge.3.0D0*PI) then
         THZ = ATAN(Z*TAN(REMTH)) + PIK + PI
      else
         THZ = ATAN(TAN(REMTH-HPI)/Z) + PIK + HPI
      end if

      DUM = ABS(COS(THZ))
      if (DUM.ge.0.5D0) then
         FAC = Z* (DUM/COS(TH))**2
      else
         FAC = (SIN(THZ)/SIN(TH))**2/Z
      end if

      DTHZ = FAC*DTH
      U(1) = THZ
      U(2) = DTHZ
      U(3) = Y(3) - 0.5D0*LOG(Z*FAC)
      return

      end subroutine THTOTHZ
      subroutine THUM(MF,ML,XS)
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
      integer MMWD
C     ..
C     .. Arrays in Common ..
      double precision YS(197)
      integer MMW(98)
C     ..
C     .. Local Scalars ..
      double precision PUP,PUP1,TMP,U,U1
      integer I,N
C     ..
C     .. Common blocks ..
      common /PASS/YS,MMW,MMWD
C     ..
C     .. Scalar Arguments ..
      integer MF,ML
C     ..
C     .. Array Arguments ..
      double precision XS(*)
C     ..
C     .. External Subroutines ..
cc      external UV
C     ..
      do 10 I = 1,98
         YS(2*I-1) = XS(I)
         YS(2*I) = 0.5D0* (XS(I)+XS(I+1))
   10 continue
      YS(197) = XS(99)
      N = 0
      U1 = 0.0D0
      PUP1 = 0.0D0
      do 20 I = 2*MF - 1,2*ML - 1
         call UV(YS(I),U,PUP,TMP,TMP,TMP,TMP)
         if (U1.lt.0.0D0 .and. U.gt.0.0D0 .and. PUP1.gt.0.0D0) then
            N = N + 1
            MMW(N) = I - 1
         end if

         U1 = U
         PUP1 = PUP
   20 continue
      MMWD = N
      return

      end subroutine THUM
      subroutine THZTOTH(U,Z,Y)
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
      double precision HPI,PI,TWOPI
C     ..
C     .. Local Scalars ..
      double precision DTH,DTHZ,DUM,FAC,PIK,REMTHZ,TH,THZ
      integer K
C     ..
C     .. Intrinsic Functions ..
      intrinsic ABS,ATAN,COS,LOG,SIN,TAN
C     ..
C     .. Common blocks ..
      common /PIE/PI,TWOPI,HPI
C     ..
C     .. Scalar Arguments ..
      double precision Z
C     ..
C     .. Array Arguments ..
      double precision U(3),Y(3)
C     ..
      THZ = U(1)
      DTHZ = U(2)
      K = THZ/PI
      if (THZ.lt.0.0D0) K = K - 1
      PIK = K*PI
      REMTHZ = THZ - PIK
      if (4.0D0*REMTHZ.le.PI) then
         TH = ATAN(TAN(REMTHZ)/Z) + PIK
      else if (4.0D0*REMTHZ.ge.3.0D0*PI) then
         TH = ATAN(TAN(REMTHZ)/Z) + PIK + PI
      else
         TH = ATAN(Z*TAN(REMTHZ-HPI)) + PIK + HPI
      end if

      DUM = ABS(COS(TH))
      if (DUM.ge.0.5D0) then
         FAC = (DUM/COS(THZ))**2/Z
      else
         FAC = Z* (SIN(TH)/SIN(THZ))**2
      end if

      DTH = FAC*DTHZ
      Y(1) = TH
      Y(2) = DTH
      Y(3) = U(3) + 0.5D0*LOG(Z/FAC)
      return

      end subroutine THZTOTH
      subroutine UVPHI(U,PUP,V,PVP,THU,THV,PHI,TH)
C     **********
C
C     This program finds TH appropriate to THU, THV, and PHI, where
C     THU is the phase angle for U, and THV is the phase angle for V.
C
C     **********
C     .. Scalars in Common ..
      double precision HPI,PI,TWOPI
C     ..
C     .. Local Scalars ..
      double precision C,D,PYP,S,Y
C     ..
C     .. External Subroutines ..
cc      external FIT
C     ..
C     .. Intrinsic Functions ..
      intrinsic ATAN2,COS,SIN
C     ..
C     .. Common blocks ..
      common /PIE/PI,TWOPI,HPI
C     ..
C     .. Scalar Arguments ..
      double precision PHI,PUP,PVP,TH,THU,THV,U,V
C     ..
      TH = THU
      if (PHI.eq.0.0D0) return
      if (THV-THU.lt.PI) then
         TH = THV
         if (PHI.eq.-HPI) return
         TH = THV - PI
         if (PHI.eq.HPI) return
      else
         TH = THV - PI
         if (PHI.eq.-HPI) return
         TH = THV - TWOPI
         if (PHI.eq.HPI) return
      end if

      C = COS(PHI)
      S = SIN(PHI)
      Y = U*C + V*S
      PYP = PUP*C + PVP*S
      TH = ATAN2(Y,PYP)
      if (Y.lt.0.0D0) TH = TH + TWOPI
      D = U*PVP - V*PUP
      if (D*PHI.gt.0.0D0) then
         call FIT(THU-PI,TH,THU)
      else
         call FIT(THU,TH,THU+PI)
      end if

      return

      end subroutine UVPHI
      subroutine WR(FG,EPS,TSTAR,YSTAR,THEND,DTHDE,TOUT,Y,TT,YY,ERR,
     +              WORK,IWORK)
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
      double precision EIG,EPSMIN
      integer IND
C     ..
C     .. Local Scalars ..
      double precision CHNG,D2F,D2G,D3F,D3G,D4F,D4G,HT,HU,OLDSS2,OLDYY2,
     +                 ONE,SLO,SOUT,SUMM,SUP,T,TEN5,TIN,U,UOUT,USTAR,
     +                 YLO,YOUT,YUP
      integer I,K,KFLAG
C     ..
C     .. Local Arrays ..
      double precision DF(4),DG(4),FF(6),GG(5),S(3),SS(6,3),UU(6)
C     ..
C     .. External Subroutines ..
cc      external GERK,LPOL
C     ..
C     .. Intrinsic Functions ..
      intrinsic ABS,SIGN
C     ..
C     .. Common blocks
C     .. Scalar Arguments ..
      double precision DTHDE,EPS,THEND,TOUT,TSTAR,YSTAR
C     ..
C     .. Array Arguments ..
      double precision ERR(3),TT(7),WORK(27),Y(3),YY(7,3)
      integer IWORK(5)
C     ..
C     .. Subroutine Arguments ..
      external FG
C     ..
C     .. Common blocks ..
      common /DATAF/EIG,IND
      common /RNDOFF/EPSMIN
C     ..
      ONE = 1.0D0
      TEN5 = 100000.0D0
C
C     Integrate Y' = F(T,Y,YP).
C
      TIN = SIGN(ONE,TSTAR)
      HT = TSTAR - TIN
      TT(1) = TIN
      YY(1,1) = THEND
      YY(1,2) = DTHDE
      YY(1,3) = 0.0D0
      TT(2) = TSTAR
      YY(2,1) = YSTAR
      YY(2,2) = DTHDE
      YY(2,3) = 0.0D0
      YLO = -TEN5
      YUP = TEN5
C
C     Normally IND = 1 or 3;  IND is set to 2 or 4 when Y is to be used
C     as the independent variable in FG(T,Y,YP), instead of the usual T.
C     Before leaving this subroutine, IND is reset to 1.
C
   10 continue
      T = TSTAR
      Y(1) = YY(2,1)
      Y(2) = YY(2,2)
      Y(3) = YY(2,3)
      KFLAG = 1
      IND = 1
      do 30 K = 3,6
         TOUT = T + HT
   20    continue
         call GERK(FG,3,Y,T,TOUT,EPS,EPS,KFLAG,ERR,WORK,IWORK)
         if (KFLAG.gt.3) then
            WRITE (021,FMT=*) ' KFLAG3 = 5 '
         end if

         if (KFLAG.eq.3) go to 20
         YOUT = Y(1)
         TT(K) = T
         YY(K,1) = YOUT
         YY(K,2) = Y(2)
         YY(K,3) = Y(3)
   30 continue
      IND = 3
      do 40 I = 2,6
         call FG(TT(I),YY(I,1),FF(I))
   40 continue
      call LPOL(5,TT(2),FF(2),TT(1),FF(1))
      if (ABS(FF(1)).le.500.0D0) then
C
C     Now we want to apply some criterion to see if these results are
C     consistent with having integrated from (TT(1),YY(1,1).
C
         do 50 I = 1,4
            DF(I) = FF(I+1) - FF(I)
   50    continue
         D2F = DF(4) - DF(3)
         D3F = DF(4) - 2.0D0*DF(3) + DF(2)
         D4F = DF(4) - 3.0D0*DF(3) + 3.0D0*DF(2) - DF(1)
         SUMM = HT* (FF(5)-3.5D0*DF(4)+53.0D0*D2F/12.0D0-55.0D0*D3F/
     +          24.0D0+251.0D0*D4F/720.0D0)
C
C     Presumably, YY(2,1) should be YY(1,1) + SUMM.
C
         OLDYY2 = YY(2,1)
         YY(2,1) = YY(1,1) + SUMM
C
C     Also improve the value of Y(2) at TSTAR.
C
         YY(2,2) = 0.5D0* (YY(1,2)+YY(3,2))
         YY(2,3) = 0.5D0* (YY(1,3)+YY(3,3))
         CHNG = YY(2,1) - OLDYY2
         if (CHNG.ge.0.0D0 .and. OLDYY2.gt.YLO) YLO = OLDYY2
         if (CHNG.le.0.0D0 .and. OLDYY2.lt.YUP) YUP = OLDYY2
         if ((YY(2,1).ge.YUP.and.YLO.gt.-TEN5) .or.
     +       (YY(2,1).le.YLO.and.YUP.lt.TEN5)) YY(2,1) = 0.5D0*
     +       (YLO+YUP)
         if (ABS(YY(2,1)-OLDYY2).gt.EPSMIN) go to 10
      else
C
C     Here, Y' is assumed infinite at T = TIN.  In this case,
C     it cannot be expected to approximate Y with a polynomial,
C     so the independent and dependent variables are interchanged.
C     The points are assumed equally spaced.
C
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
   60    continue
         U = USTAR
         S(1) = SS(2,1)
         S(2) = SS(2,2)
         S(3) = SS(2,3)
         KFLAG = 1
         IND = 2
         do 80 K = 3,6
            UOUT = U + HU
   70       continue
            call GERK(FG,3,S,U,UOUT,EPS,EPS,KFLAG,ERR,WORK,IWORK)
            if (KFLAG.gt.3) then
               WRITE (021,FMT=*) ' KFLAG4 = 5 '
            end if

            if (KFLAG.eq.3) go to 70
            SOUT = S(1)
            UU(K) = U
            SS(K,1) = SOUT
            SS(K,2) = S(2)
            SS(K,3) = S(3)
   80    continue
         IND = 4
         do 90 I = 2,5
            call FG(UU(I),SS(I,1),GG(I))
   90    continue
         GG(1) = 0.0D0
         do 100 I = 1,4
            DG(I) = GG(I+1) - GG(I)
  100    continue
         D2G = DG(4) - DG(3)
         D3G = DG(4) - 2.0D0*DG(3) + DG(2)
         D4G = DG(4) - 3.0D0*DG(3) + 3.0D0*DG(2) - DG(1)
         SUMM = HU* (GG(5)-3.5D0*DG(4)+53.0D0*D2G/12.0D0-55.0D0*D3G/
     +          24.0D0+251.0D0*D4G/720.0D0)
C
C     Presumably, SS(2,1) should be SS(1,1) + SUMM.
C
         OLDSS2 = SS(2,1)
         SS(2,1) = SS(1,1) + SUMM
         if (SS(2,1).le.-1.0D0) SS(2,1) = -1.0D0 + EPSMIN
         if (SS(2,1).ge.1.0D0) SS(2,1) = 1.0D0 - EPSMIN
C
C     Also improve the value of Y(2) at TSTAR.
C
         SS(2,2) = 0.5D0* (SS(1,2)+SS(3,2))
         SS(2,3) = 0.5D0* (SS(1,3)+SS(3,3))
         CHNG = SS(2,1) - OLDSS2
         if (CHNG.ge.0.0D0 .and. OLDSS2.gt.SLO) SLO = OLDSS2
         if (CHNG.le.0.0D0 .and. OLDSS2.lt.SUP) SUP = OLDSS2
         if ((SS(2,1).ge.SUP.and.SLO.gt.-TEN5) .or.
     +       (SS(2,1).le.SLO.and.SUP.lt.TEN5)) SS(2,1) = 0.5D0*
     +       (SLO+SUP)
         if (ABS(SS(2,1)-OLDSS2).gt.EPSMIN) go to 60
      end if

      if (IND.eq.4) then
         do 110 I = 1,6
            TT(I) = SS(I,1)
            YY(I,1) = UU(I)
            YY(I,2) = SS(I,2)
            YY(I,3) = SS(I,3)
  110    continue
      end if

      TOUT = TT(6)
      IND = 1
      return

      end subroutine WR
      subroutine GERK(F,NEQN,Y,T,TOUT,RELERR,ABSERR,IFLAG,GERROR,WORK,
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
      double precision ABSERR,RELERR,T,TOUT
      integer IFLAG,NEQN
C     ..
C     .. Array Arguments ..
      double precision GERROR(NEQN),WORK(3+8*NEQN),Y(NEQN)
      integer IWORK(5)
C     ..
C     .. Subroutine Arguments ..
      external F
C     ..
C     .. Local Scalars ..
      integer K1,K1M,K2,K3,K4,K5,K6,K7,K8
C     ..
C     .. External Subroutines ..
cc      external GERKS
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
      call GERKS(F,NEQN,Y,T,TOUT,RELERR,ABSERR,IFLAG,GERROR,WORK(1),
     +           WORK(K1M),WORK(K1),WORK(K2),WORK(K3),WORK(K4),WORK(K5),
     +           WORK(K6),WORK(K7),WORK(K8),WORK(K8+1),IWORK(1),
     +           IWORK(2),IWORK(3),IWORK(4),IWORK(5))
      return

      end subroutine GERK
      subroutine GERKS(F,NEQN,Y,T,TOUT,RELERR,ABSERR,IFLAG,GERROR,YP,H,
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
      double precision ABSERR,H,RELERR,SAVAE,SAVRE,T,TOUT
      integer IFLAG,INIT,JFLAG,KFLAG,KOP,NEQN,NFE
C     ..
C     .. Array Arguments ..
      double precision F1(NEQN),F2(NEQN),F3(NEQN),F4(NEQN),F5(NEQN),
     +                 GERROR(NEQN),Y(NEQN),YG(NEQN),YGP(NEQN),YP(NEQN)
C     ..
C     .. Subroutine Arguments ..
      external F
C     ..
C     .. Local Scalars ..
      double precision A,AE,DT,EE,EEOET,ESTTOL,ET,HH,HMIN,ONE,REMIN,RER,
     +                 S,SCALE,TOL,TOLN,TS,U,U26,YPK
      integer K,MAXNFE,MFLAG
      logical HFAILD,OUTPUT
C     ..
C     .. External Functions ..
cc      double precision EPSLON
cc      external EPSLON
C     ..
C     .. External Subroutines ..
cc      external FEHL
C     ..
C     .. Intrinsic Functions ..
      intrinsic ABS,MAX,MIN,SIGN
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
      data REMIN/3.D-11/
      data MAXNFE/9000/
C     ..
      ONE = 1.0D0
      U = EPSLON(ONE)
C *******************************************************************
C      CHECK INPUT PARAMETERS
      if (NEQN.lt.1) go to 10
      if ((RELERR.lt.0.D0) .or. (ABSERR.lt.0.D0)) go to 10
      MFLAG = ABS(IFLAG)
      if ((MFLAG.ge.1) .and. (MFLAG.le.7)) go to 20
C INVALID INPUT
   10 IFLAG = 7
      return
C IS THIS THE FIRST CALL
   20 if (MFLAG.eq.1) go to 70
C CHECK CONTINUATION POSSIBILITIES
      if (T.eq.TOUT) go to 10
      if (MFLAG.ne.2) go to 30
C IFLAG = +2 OR -2
      if (INIT.eq.0) go to 60
      if (KFLAG.eq.3) go to 50
      if ((KFLAG.eq.4) .and. (ABSERR.eq.0.D0)) go to 40
      if ((KFLAG.eq.5) .and. (RELERR.le.SAVRE) .and.
     +    (ABSERR.le.SAVAE)) go to 40
      go to 70
C IFLAG = 3,4,5,6, OR 7
   30 if (IFLAG.eq.3) go to 50
      if ((IFLAG.eq.4) .and. (ABSERR.gt.0.D0)) go to 60
C INTEGRATION CANNOT BE CONTINUED SINCE USER DID NOT RESPOND TO
C THE INSTRUCTIONS PERTAINING TO IFLAG=4,5,6 OR 7
   40 stop
C *******************************************************************
C      RESET FUNCTION EVALUATION COUNTER
   50 NFE = 0
      if (MFLAG.eq.2) go to 70
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
      RER = MAX(RELERR,32.D0*U+REMIN)
      U26 = 26.D0*U
      DT = TOUT - T
      if (MFLAG.eq.1) go to 80
      if (INIT.eq.0) go to 90
      go to 110
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
      call F(A,Y,YP)
      NFE = 1
      if (T.ne.TOUT) go to 90
      IFLAG = 2
      return

   90 INIT = 1
      H = ABS(DT)
      TOLN = 0.D0
      do 100 K = 1,NEQN
         YG(K) = Y(K)
         YGP(K) = YP(K)
         TOL = RER*ABS(Y(K)) + ABSERR
         if (TOL.le.0.D0) go to 100
         TOLN = TOL
         YPK = ABS(YP(K))
         if (YPK*H**5.gt.TOL) H = (TOL/YPK)**0.2D0
  100 continue
      if (TOLN.le.0.D0) H = 0.D0
      H = MAX(H,U26*MAX(ABS(T),ABS(DT)))
C *******************************************************************
C      SET STEPSIZE FOR INTEGRATION IN THE DIRECTION FROM T TO TOUT
  110 H = SIGN(H,DT)
C TEST TO SEE IF GERK IS BEING SEVERELY IMPACTED BY TOO MANY
C OUTPUT POINTS
      if (ABS(H).gt.2.D0*ABS(DT)) KOP = KOP + 1
      if (KOP.ne.100) go to 120
      KOP = 0
      IFLAG = 6
      return

  120 if (ABS(DT).gt.U26*ABS(T)) go to 140
C IF TOO CLOSE TO OUTPUT POINT,EXTRAPOLATE AND RETURN
      do 130 K = 1,NEQN
         YG(K) = YG(K) + DT*YGP(K)
         Y(K) = Y(K) + DT*YP(K)
  130 continue
      A = TOUT
      call F(A,YG,YGP)
      call F(A,Y,YP)
      NFE = NFE + 2
      go to 230
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
      if (ABS(DT).ge.2.D0*ABS(H)) go to 170
      if (ABS(DT).gt.ABS(H)) go to 160
C THE NEXT SUCCESSFUL STEP WILL COMPLETE THE INTEGRATION TO THE
C OUTPUT POINT
      OUTPUT = .TRUE.
      H = DT
      go to 170

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
  170 if (NFE.le.MAXNFE) go to 180
C TOO MUCH WORK
      IFLAG = 3
      KFLAG = 3
      return
C ADVANCE AN APPROXIMATE SOLUTION OVER ONE STEP OF LENGTH H
  180 call FEHL(F,NEQN,YG,T,H,YGP,F1,F2,F3,F4,F5,F1)
      NFE = NFE + 5
C COMPUTE AND TEST ALLOWABLE TOLERANCES VERSUS LOCAL ERROR
C ESTIMATES AND REMOVE SCALING OF TOLERANCES. NOTE THAT RELATIVE
C ERROR IS MEASURED WITH RESPECT TO THE AVERAGE MAGNITUDES OF THE
C OF THE SOLUTION AT THE BEGINNING AND END OF THE STEP.
      EEOET = 0.D0
      do 200 K = 1,NEQN
         ET = ABS(YG(K)) + ABS(F1(K)) + AE
         if (ET.gt.0.D0) go to 190
C INAPPROPRIATE ERROR TOLERANCE
         IFLAG = 4
         KFLAG = 4
         return

  190    EE = ABS((-2090.D0*YGP(K)+ (21970.D0*F3(K)-15048.D0*F4(K)))+
     +        (22528.D0*F2(K)-27360.D0*F5(K)))
         EEOET = MAX(EEOET,EE/ET)
  200 continue
      ESTTOL = ABS(H)*EEOET*SCALE/752400.D0
      if (ESTTOL.le.1.D0) go to 210
C UNSUCCESSFUL STEP
C                   REDUCE THE STEPSIZE , TRY AGAIN
C                   THE DECREASE IS LIMITED TO A FACTOR OF 1/10
      HFAILD = .TRUE.
      OUTPUT = .FALSE.
      S = 0.1D0
      if (ESTTOL.lt.59049.D0) S = 0.9D0/ESTTOL**0.2D0
      H = S*H
      if (ABS(H).gt.HMIN) go to 170
C REQUESTED ERROR UNATTAINABLE AT SMALLEST ALLOWABLE STEPSIZE
      IFLAG = 5
      KFLAG = 5
      return
C SUCCESSFUL STEP
C                    STORE ONE-STEP SOLUTION YG AT T+H
C                    AND EVALUATE DERIVATIVES THERE
  210 TS = T
      T = T + H
      do 220 K = 1,NEQN
         YG(K) = F1(K)
  220 continue
      A = T
      call F(A,YG,YGP)
      NFE = NFE + 1
C NOW ADVANCE THE Y SOLUTION OVER TWO STEPS OF
C LENGTH H/2 AND EVALUATE DERIVATIVES THERE
      HH = 0.5D0*H
      call FEHL(F,NEQN,Y,TS,HH,YP,F1,F2,F3,F4,F5,Y)
      TS = TS + HH
      A = TS
      call F(A,Y,YP)
      call FEHL(F,NEQN,Y,TS,HH,YP,F1,F2,F3,F4,F5,Y)
      A = T
      call F(A,Y,YP)
      NFE = NFE + 12
C CHOOSE NEXT STEPSIZE
C THE INCREASE IS LIMITED TO A FACTOR OF 5
C IF STEP FAILURE HAS JUST OCCURRED, NEXT
C    STEPSIZE IS NOT ALLOWED TO INCREASE
      S = 5.D0
      if (ESTTOL.gt.1.889568D-4) S = 0.9D0/ESTTOL**0.2D0
      if (HFAILD) S = MIN(S,ONE)
      H = SIGN(MAX(S*ABS(H),HMIN),H)
C *******************************************************************
C      END OF CORE INTEGRATOR
C *******************************************************************
C      SHOULD WE TAKE ANOTHER STEP
      if (OUTPUT) go to 230
      if (IFLAG.gt.0) go to 150
C *******************************************************************
C *******************************************************************
C      INTEGRATION SUCCESSFULLY COMPLETED
C      ONE-STEP MODE
      IFLAG = -2
      go to 240
C INTERVAL MODE
  230 T = TOUT
      IFLAG = 2
  240 do 250 K = 1,NEQN
         GERROR(K) = (YG(K)-Y(K))/31.D0
  250 continue
      return

      end subroutine GERKS
      subroutine FEHL(F,NEQN,Y,T,H,YP,F1,F2,F3,F4,F5,S)
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
      double precision H,T
      integer NEQN
C     ..
C     .. Array Arguments ..
      double precision F1(NEQN),F2(NEQN),F3(NEQN),F4(NEQN),F5(NEQN),
     +                 S(NEQN),Y(NEQN),YP(NEQN)
C     ..
C     .. Subroutine Arguments ..
      external F
C     ..
C     .. Local Scalars ..
      double precision CH
      integer K
C     ..
      CH = 0.25D0*H
      do 10 K = 1,NEQN
         F5(K) = Y(K) + CH*YP(K)
   10 continue
      call F(T+0.25D0*H,F5,F1)
      CH = 0.09375D0*H
      do 20 K = 1,NEQN
         F5(K) = Y(K) + CH* (YP(K)+3.D0*F1(K))
   20 continue
      call F(T+0.375D0*H,F5,F2)
      CH = H/2197.D0
      do 30 K = 1,NEQN
         F5(K) = Y(K) + CH* (1932.D0*YP(K)+
     +           (7296.D0*F2(K)-7200.D0*F1(K)))
   30 continue
      call F(T+12.D0/13.D0*H,F5,F3)
      CH = H/4104.D0
      do 40 K = 1,NEQN
         F5(K) = Y(K) + CH* ((8341.D0*YP(K)-845.D0*F3(K))+
     +           (29440.D0*F2(K)-32832.D0*F1(K)))
   40 continue
      call F(T+H,F5,F4)
      CH = H/20520.D0
      do 50 K = 1,NEQN
         F1(K) = Y(K) + CH* ((-6080.D0*YP(K)+ (9295.D0*F3(K)-
     +           5643.D0*F4(K)))+ (41040.D0*F1(K)-28352.D0*F2(K)))
   50 continue
      call F(T+0.5D0*H,F1,F5)
C COMPUTE APPROXIMATE SOLUTION AT T+H
      CH = H/7618050.D0
      do 60 K = 1,NEQN
         S(K) = Y(K) + CH* ((902880.D0*YP(K)+ (3855735.D0*F3(K)-
     +          1371249.D0*F4(K)))+ (3953664.D0*F2(K)+277020.D0*F5(K)))
   60 continue
      return

      end subroutine FEHL
      double precision function EPSLON(X)
C
C     ESTIMATE UNIT ROUNDOFF IN QUANTITIES OF SIZE X.
C
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
      double precision X
C     ..
C     .. Local Scalars ..
      double precision A,B,C,EPS,FOUR,THREE
C     ..
C     .. Intrinsic Functions ..
      intrinsic ABS
C     ..
      FOUR = 4.0D0
      THREE = 3.0D0
      A = FOUR/THREE
   10 B = A - 1.0D0
      C = B + B + B
      EPS = ABS(C-1.0D0)
      if (EPS.eq.0.0D0) go to 10
      EPSLON = EPS*ABS(X)
      return

      end function EPSLON
      subroutine FZERO(F,B,C,R,RE,AE,IFLAG)
C     **********
C     **********
C     .. Local Scalars ..
      double precision A,ACBS,ACMB,AW,CMB,DIF,DIFS,FA,FB,FC,FX,FZ,P,Q,
     +                 RW,TOL,Z
      integer IC,KOUNT
C     ..
C     .. Intrinsic Functions ..
      intrinsic ABS,MAX,MIN,SIGN
C     ..
C     .. Scalar Arguments ..
      double precision AE,B,C,R,RE
      integer IFLAG
C     ..
C     .. Function Arguments ..
      double precision F
      external F
C     ..
      DIF = 1000.0D0
      Z = R
      if (R.le.MIN(B,C) .or. R.ge.MAX(B,C)) Z = C
      RW = MAX(RE,0.0D0)
      AW = MAX(AE,0.0D0)
      IC = 0
      FZ = F(Z)
      FB = F(B)
      KOUNT = 2
      if (FZ*FB.lt.0.0D0) then
         C = Z
         FC = FZ
      else if (Z.ne.C) then
         FC = F(C)
         KOUNT = 3
         if (FZ*FC.lt.0.0D0) then
            B = Z
            FB = FZ
         end if

      end if

      A = C
      FA = FC
      ACBS = ABS(B-C)
      FX = MAX(ABS(FB),ABS(FC))
C
   10 continue
      if (ABS(FC).lt.ABS(FB)) then
         A = B
         FA = FB
         B = C
         FB = FC
         C = A
         FC = FA
      end if

      CMB = 0.5D0* (C-B)
      ACMB = ABS(CMB)
      TOL = RW*ABS(B) + AW
      IFLAG = 1
      if (ACMB.le.TOL) then
         if (FB*FC.ge.0.0D0) IFLAG = 4
         if (ABS(FB).gt.FX) IFLAG = 3
         return

      end if

      IFLAG = 2
      if (FB.eq.0.0D0) return
      IFLAG = 5
      if (KOUNT.ge.500) return
C
      P = (B-A)*FB
      Q = FA - FB
      if (P.lt.0.0D0) then
         P = -P
         Q = -Q
      end if

      A = B
      FA = FB
      IC = IC + 1
      if (IC.ge.4) then
         if (8.0D0*ACMB.ge.ACBS) B = 0.5D0* (C+B)
      else
         IC = 0
         ACBS = ACMB
         if (P.le.ABS(Q)*TOL) then
            B = B + SIGN(TOL,CMB)
         else if (P.ge.CMB*Q) then
            B = 0.5D0* (C+B)
         else
            B = B + P/Q
         end if

      end if

      FB = F(B)
      DIFS = DIF
      DIF = FB - FC
      if (DIF.eq.DIFS) then
         IFLAG = 6
         return

      end if

      KOUNT = KOUNT + 1
      if (FB*FC.ge.0.0D0) then
         C = A
         FC = FA
      end if

      go to 10

      end subroutine FZERO

      end module sleig2md
