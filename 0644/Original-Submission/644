C      REMARK ON ALGORITHM 644, COLLECTED ALGORITHMS FROM ACM.
C      THIS WORK PUBLISHED IN TRANSACTIONS ON MATHEMATICAL SOFTWARE,
C      VOL. 21, NO. 4, December, 1995, P.  388--393.
C
C This file contains 6 files separated by lines of the form
C         C*** filename
C
C The filenames in this file are:
C
C Readme               cbsubs.f             cqcbes.f            
C machcon.f            zbsubs.f             zqcbes.f            
C                                                               
C
C*** Readme
             INSTALLATION HINTS FOR ALGORITHM 644

                          D. E. AMOS
                 Sandia National Laboratories
-----------------------------------------------------------------

Algorithm 644 [2],[3] computes all major Bessel functions of a
complex argument and nonnegative order.  The single precision
callable routines are denoted by CBESH, CBESI, CBESJ, CBESK,
CBESY, CAIRY AND CBIRY for the H, I, J, K, Y, and Airy functions
Ai and Bi. The corresponding double precision routines start with
Z in place of C. CBESH and ZBESH compute H functions of kinds 1
and 2 and an option exists for the derivatives of the Airy
functions in CAIRY, ZAIRY, CBIRY, and ZBIRY. Scaling to remove
exponential underflow or overflow and generation of sequences in
increasing orders are auxilliary options where appropriate.

It is convenient to consider the package in four parts:

   (1) Quick check PROGRAMS (CQC?? or ZQC??) which exercise the
       main callable routines in single and double precision
       arithmetic:

          CQCBH and ZQCBH  to test CBESH and ZBESH
          CQCBI and ZQCBI  to test CBESI and ZBESI
          CQCBJ and ZQCBJ  to test CBESJ and ZBESJ
          CQCBK and ZQCBK  to test CBESK and ZBESK
          CQCBY and ZQCBY  to test CBESY and ZBESY
          CQCAI and ZQCAI  to test CAIRY, ZAIRY, CBIRY and ZBIRY
          and auxilliary subroutines CBESYH and ZBESYH

       The documentation for the latest versions of these quick
       checks is found in reference [4]. These routines need not
       be a permanent part of the installation once they have
       been compiled and successfully executed.

   (2) The main callable subroutines which comprise the algorithm
       in single and double precision arithmetic:

          CBESH and ZBESH  for H functions of kinds 1 and 2
          CBESI and ZBESI  for I functions
          CBESJ and ZBESJ  for J functions
          CBESK and ZBESK  for K functions
          CBESY and ZBESY  for Y functions
          CAIRY and ZAIRY  for Airy function Ai and Ai'
          CBIRY and ZBIRY  for Airy function Bi and Bi'
          GAMLN and DGAMLN for lngamma(x), x>0

          The Z routines manipulate double precision ordered
          pairs as double precision complex numbers and return
          double precison ordered pairs as answers.  The lower
          level routines

          ZMLT,ZDIV,ZSQRT,ZEXP,ZLOG,ZSHCH,ZABS

          provide a means of doing double precision complex
          arithmetic on double precision ordered pairs. These can
          also be called by user-coded routines to accomplish
          double precision complex arithmetic. There are no ZADD
          or ZSUB routines -- double precision complex additions
          and subtractions are simple enough to do in-line.

          Each subroutine or function has a prologue which
          defines the calling sequence. These routines, along
          with those of (3) and (4), make up the permanent
          installation.

   (3) Included in the algorithm are three (callable) function
       subroutines which define the machine to the Bessel
       function package.  These subroutines consist mainly of
       comment lines which contain (inactive) FORTRAN statements
       for machine constants such as real and double precision
       unit roundoff, real and double precision overflow and
       underflow limits, largest integer constant, base of
       integer, real and double precision arithmetic, etc. for a
       variety of machines. The prologue of each of these
       functions

       INTEGER          I1MACH(I)   I=1,16  for integer constants
       REAL             R1MACH(I)   I=1,5   for real    constants
       DOUBLE PRECISION D1MACH(I)   I=1,5   for double  precision
                                            constants

       defines the return value for each of the calling indices.

       TO USE THE PACKAGE, YOU MUST DEFINE THE MACHINE ON WHICH
       THE PACKAGE IS TO BE USED BY EDITING EACH OF THESE
       FUNCTIONS.

       With an editor, search out your machine from among those
       listed in each ?1MACH function. A keyword like IBM, SUN,
       VAX, CDC or PC is often helpful in locating the right
       machine.  When your machine has been located, replace the
       C in column 1 with a space to make each of the FORTRAN
       statements active for the selected machine.

       If your machine is not found, the correct constants can be
       computed from the definitions in the prologue. It often
       happens that many PC constants are close enough to the IBM
       or SUN constants to work with the package. Reference [5]
       should help in generating a set of machine constants. One
       does not need constants accurate to the last bit and
       leaving a little slack on arithmetic limits can sometimes
       avoid hardware or software problems.

   (4) Lower level non-callable routines

The latest version containing all improvements is dated 930101.

REFERENCES

1. Abramowitz, M. and Stegun, I. A., Handbook of Mathematical
   Functions, NBS Applied Math Series 55, U.S. Dept. of Commerce,
   Washington, D.C., 1955

2. Amos, D. E., Algorithm 644, A Portable Package For Bessel
   Functions of a Complex Argument and Nonnegative Order, ACM
   Transactions on Mathematical Software, Vol. 12, No. 3,
   September 1986, Pages 265-273

3. Amos, D. E., Remark on Algorithm 644, ACM Transactions on
   Mathematical Software, Vol. 16, No. 4, December 1990, Page 404

4. Amos, D. E., Remark on Algorithm 644 (Improvements in
   Algorithm 644), ACM Transactions on Mathematical Software,
   (Estimated date 1995)

5. Cody, W. J., Algorithm 665, MACHAR: A Subroutine to
   Dynamically Determine Machine Parameters, ACM Transactions on
   Mathematical Software, Vol. 14, No. 4, December 1988, Pages
   303-311

C*** cbsubs.f
-----------------------------------------------------------------
C>>>  CBSUBS.FOR:  Single precision subroutines
-----------------------------------------------------------------

      SUBROUTINE CBESH(Z, FNU, KODE, M, N, CY, NZ, IERR)
C***BEGIN PROLOGUE  CBESH
C***DATE WRITTEN   830501   (YYMMDD)
C***REVISION DATE  890801, 930101   (YYMMDD)
C***CATEGORY NO.  B5K
C***KEYWORDS  H-BESSEL FUNCTIONS,BESSEL FUNCTIONS OF COMPLEX ARGUMENT,
C             BESSEL FUNCTIONS OF THIRD KIND,HANKEL FUNCTIONS
C***AUTHOR  AMOS, DONALD E., SANDIA NATIONAL LABORATORIES
C***PURPOSE  TO COMPUTE THE H-BESSEL FUNCTIONS OF A COMPLEX ARGUMENT
C***DESCRIPTION
C
C         ON KODE=1, CBESH COMPUTES AN N MEMBER SEQUENCE OF COMPLEX
C         HANKEL (BESSEL) FUNCTIONS CY(J)=H(M,FNU+J-1,Z) FOR KINDS M=1
C         OR 2, REAL, NONNEGATIVE ORDERS FNU+J-1, J=1,...,N, AND COMPLEX
C         Z.NE.CMPLX(0.0E0,0.0E0) IN THE CUT PLANE -PI.LT.ARG(Z).LE.PI.
C         ON KODE=2, CBESH COMPUTES THE SCALED HANKEL FUNCTIONS
C
C         CY(I)=H(M,FNU+J-1,Z)*EXP(-MM*Z*I)       MM=3-2M,      I**2=-1.
C
C         WHICH REMOVES THE EXPONENTIAL BEHAVIOR IN BOTH THE UPPER
C         AND LOWER HALF PLANES. DEFINITIONS AND NOTATION ARE FOUND IN
C         THE NBS HANDBOOK OF MATHEMATICAL FUNCTIONS (REF. 1).
C
C         INPUT
C           Z      - Z=CMPLX(X,Y), Z.NE.CMPLX(0.,0.),-PI.LT.ARG(Z).LE.PI
C           FNU    - ORDER OF INITIAL H FUNCTION, FNU.GE.0.0E0
C           KODE   - A PARAMETER TO INDICATE THE SCALING OPTION
C                    KODE= 1  RETURNS
C                             CY(J)=H(M,FNU+J-1,Z),      J=1,...,N
C                        = 2  RETURNS
C                             CY(J)=H(M,FNU+J-1,Z)*EXP(-I*Z*(3-2M))
C                                  J=1,...,N  ,  I**2=-1
C           M      - KIND OF HANKEL FUNCTION, M=1 OR 2
C           N      - NUMBER OF MEMBERS OF THE SEQUENCE, N.GE.1
C
C         OUTPUT
C           CY     - A COMPLEX VECTOR WHOSE FIRST N COMPONENTS CONTAIN
C                    VALUES FOR THE SEQUENCE
C                    CY(J)=H(M,FNU+J-1,Z)  OR
C                    CY(J)=H(M,FNU+J-1,Z)*EXP(-I*Z*(3-2M))  J=1,...,N
C                    DEPENDING ON KODE, I**2=-1.
C           NZ     - NUMBER OF COMPONENTS SET TO ZERO DUE TO UNDERFLOW,
C                    NZ= 0   , NORMAL RETURN
C                    NZ.GT.0 , FIRST NZ COMPONENTS OF CY SET TO ZERO
C                              DUE TO UNDERFLOW, CY(J)=CMPLX(0.0,0.0)
C                              J=1,...,NZ WHEN Y.GT.0.0 AND M=1 OR
C                              Y.LT.0.0 AND M=2. FOR THE COMPLMENTARY
C                              HALF PLANES, NZ STATES ONLY THE NUMBER
C                              OF UNDERFLOWS.
C           IERR    -ERROR FLAG
C                    IERR=0, NORMAL RETURN - COMPUTATION COMPLETED
C                    IERR=1, INPUT ERROR   - NO COMPUTATION
C                    IERR=2, OVERFLOW      - NO COMPUTATION, FNU+N-1 TOO
C                            LARGE OR CABS(Z) TOO SMALL OR BOTH
C                    IERR=3, CABS(Z) OR FNU+N-1 LARGE - COMPUTATION DONE
C                            BUT LOSSES OF SIGNIFCANCE BY ARGUMENT
C                            REDUCTION PRODUCE LESS THAN HALF OF MACHINE
C                            ACCURACY
C                    IERR=4, CABS(Z) OR FNU+N-1 TOO LARGE - NO COMPUTA-
C                            TION BECAUSE OF COMPLETE LOSSES OF SIGNIFI-
C                            CANCE BY ARGUMENT REDUCTION
C                    IERR=5, ERROR              - NO COMPUTATION,
C                            ALGORITHM TERMINATION CONDITION NOT MET
C
C***LONG DESCRIPTION
C
C         THE COMPUTATION IS CARRIED OUT BY THE RELATION
C
C         H(M,FNU,Z)=(1/MP)*EXP(-MP*FNU)*K(FNU,Z*EXP(-MP))
C             MP=MM*HPI*I,  MM=3-2*M,  HPI=PI/2,  I**2=-1
C
C         FOR M=1 OR 2 WHERE THE K BESSEL FUNCTION IS COMPUTED FOR THE
C         RIGHT HALF PLANE RE(Z).GE.0.0. THE K FUNCTION IS CONTINUED
C         TO THE LEFT HALF PLANE BY THE RELATION
C
C         K(FNU,Z*EXP(MP)) = EXP(-MP*FNU)*K(FNU,Z)-MP*I(FNU,Z)
C         MP=MR*PI*I, MR=+1 OR -1, RE(Z).GT.0, I**2=-1
C
C         WHERE I(FNU,Z) IS THE I BESSEL FUNCTION.
C
C         EXPONENTIAL DECAY OF H(M,FNU,Z) OCCURS IN THE UPPER HALF Z
C         PLANE FOR M=1 AND THE LOWER HALF Z PLANE FOR M=2.  EXPONENTIAL
C         GROWTH OCCURS IN THE COMPLEMENTARY HALF PLANES.  SCALING
C         BY EXP(-MM*Z*I) REMOVES THE EXPONENTIAL BEHAVIOR IN THE
C         WHOLE Z PLANE FOR Z TO INFINITY.
C
C         FOR NEGATIVE ORDERS,THE FORMULAE
C
C               H(1,-FNU,Z) = H(1,FNU,Z)*CEXP( PI*FNU*I)
C               H(2,-FNU,Z) = H(2,FNU,Z)*CEXP(-PI*FNU*I)
C                         I**2=-1
C
C         CAN BE USED.
C
C         IN MOST COMPLEX VARIABLE COMPUTATION, ONE MUST EVALUATE ELE-
C         MENTARY FUNCTIONS. WHEN THE MAGNITUDE OF Z OR FNU+N-1 IS
C         LARGE, LOSSES OF SIGNIFICANCE BY ARGUMENT REDUCTION OCCUR.
C         CONSEQUENTLY, IF EITHER ONE EXCEEDS U1=SQRT(0.5/UR), THEN
C         LOSSES EXCEEDING HALF PRECISION ARE LIKELY AND AN ERROR FLAG
C         IERR=3 IS TRIGGERED WHERE UR=R1MACH(4)=UNIT ROUNDOFF. ALSO
C         IF EITHER IS LARGER THAN U2=0.5/UR, THEN ALL SIGNIFICANCE IS
C         LOST AND IERR=4. IN ORDER TO USE THE INT FUNCTION, ARGUMENTS
C         MUST BE FURTHER RESTRICTED NOT TO EXCEED THE LARGEST MACHINE
C         INTEGER, U3=I1MACH(9). THUS, THE MAGNITUDE OF Z AND FNU+N-1 IS
C         RESTRICTED BY MIN(U2,U3). ON 32 BIT MACHINES, U1,U2, AND U3
C         ARE APPROXIMATELY 2.0E+3, 4.2E+6, 2.1E+9 IN SINGLE PRECISION
C         ARITHMETIC AND 1.3E+8, 1.8E+16, 2.1E+9 IN DOUBLE PRECISION
C         ARITHMETIC RESPECTIVELY. THIS MAKES U2 AND U3 LIMITING IN
C         THEIR RESPECTIVE ARITHMETICS. THIS MEANS THAT ONE CAN EXPECT
C         TO RETAIN, IN THE WORST CASES ON 32 BIT MACHINES, NO DIGITS
C         IN SINGLE AND ONLY 7 DIGITS IN DOUBLE PRECISION ARITHMETIC.
C         SIMILAR CONSIDERATIONS HOLD FOR OTHER MACHINES.
C
C         THE APPROXIMATE RELATIVE ERROR IN THE MAGNITUDE OF A COMPLEX
C         BESSEL FUNCTION CAN BE EXPRESSED BY P*10**S WHERE P=MAX(UNIT
C         ROUNDOFF,1.0E-18) IS THE NOMINAL PRECISION AND 10**S REPRE-
C         SENTS THE INCREASE IN ERROR DUE TO ARGUMENT REDUCTION IN THE
C         ELEMENTARY FUNCTIONS. HERE, S=MAX(1,ABS(LOG10(CABS(Z))),
C         ABS(LOG10(FNU))) APPROXIMATELY (I.E. S=MAX(1,ABS(EXPONENT OF
C         CABS(Z),ABS(EXPONENT OF FNU)) ). HOWEVER, THE PHASE ANGLE MAY
C         HAVE ONLY ABSOLUTE ACCURACY. THIS IS MOST LIKELY TO OCCUR WHEN
C         ONE COMPONENT (IN ABSOLUTE VALUE) IS LARGER THAN THE OTHER BY
C         SEVERAL ORDERS OF MAGNITUDE. IF ONE COMPONENT IS 10**K LARGER
C         THAN THE OTHER, THEN ONE CAN EXPECT ONLY MAX(ABS(LOG10(P))-K,
C         0) SIGNIFICANT DIGITS; OR, STATED ANOTHER WAY, WHEN K EXCEEDS
C         THE EXPONENT OF P, NO SIGNIFICANT DIGITS REMAIN IN THE SMALLER
C         COMPONENT. HOWEVER, THE PHASE ANGLE RETAINS ABSOLUTE ACCURACY
C         BECAUSE, IN COMPLEX ARITHMETIC WITH PRECISION P, THE SMALLER
C         COMPONENT WILL NOT (AS A RULE) DECREASE BELOW P TIMES THE
C         MAGNITUDE OF THE LARGER COMPONENT. IN THESE EXTREME CASES,
C         THE PRINCIPAL PHASE ANGLE IS ON THE ORDER OF +P, -P, PI/2-P,
C         OR -PI/2+P.
C
C***REFERENCES  HANDBOOK OF MATHEMATICAL FUNCTIONS BY M. ABRAMOWITZ
C                 AND I. A. STEGUN, NBS AMS SERIES 55, U.S. DEPT. OF
C                 COMMERCE, 1955.
C
C               COMPUTATION OF BESSEL FUNCTIONS OF COMPLEX ARGUMENT
C                 BY D. E. AMOS, SAND83-0083, MAY, 1983.
C
C               COMPUTATION OF BESSEL FUNCTIONS OF COMPLEX ARGUMENT
C                 AND LARGE ORDER BY D. E. AMOS, SAND83-0643, MAY, 1983
C
C               A SUBROUTINE PACKAGE FOR BESSEL FUNCTIONS OF A COMPLEX
C                 ARGUMENT AND NONNEGATIVE ORDER BY D. E. AMOS, SAND85-
C                 1018, MAY, 1985
C
C               A PORTABLE PACKAGE FOR BESSEL FUNCTIONS OF A COMPLEX
C                 ARGUMENT AND NONNEGATIVE ORDER BY D. E. AMOS, ACM
C                 TRANS. MATH. SOFTWARE, VOL. 12, NO. 3, SEPTEMBER 1986,
C                 PP 265-273.
C
C***ROUTINES CALLED  CACON,CBKNU,CBUNK,CUOIK,I1MACH,R1MACH
C***END PROLOGUE  CBESH
C
      COMPLEX CY, Z, ZN, ZT, CSGN
      REAL AA, ALIM, ALN, ARG, AZ, CPN, DIG, ELIM, FMM, FN, FNU, FNUL,
     * HPI, RHPI, RL, R1M5, SGN, SPN, TOL, UFL, XN, XX, YN, YY, R1MACH,
     * BB, ASCLE, RTOL, ATOL
      INTEGER I, IERR, INU, INUH, IR, K, KODE, K1, K2, M,
     * MM, MR, N, NN, NUF, NW, NZ, I1MACH
      DIMENSION CY(N)
C
      DATA HPI /1.57079632679489662E0/
C
C***FIRST EXECUTABLE STATEMENT  CBESH
      NZ=0
      XX = REAL(Z)
      YY = AIMAG(Z)
      IERR = 0
      IF (XX.EQ.0.0E0 .AND. YY.EQ.0.0E0) IERR=1
      IF (FNU.LT.0.0E0) IERR=1
      IF (M.LT.1 .OR. M.GT.2) IERR=1
      IF (KODE.LT.1 .OR. KODE.GT.2) IERR=1
      IF (N.LT.1) IERR=1
      IF (IERR.NE.0) RETURN
      NN = N
C-----------------------------------------------------------------------
C     SET PARAMETERS RELATED TO MACHINE CONSTANTS.
C     TOL IS THE APPROXIMATE UNIT ROUNDOFF LIMITED TO 1.0E-18.
C     ELIM IS THE APPROXIMATE EXPONENTIAL OVER- AND UNDERFLOW LIMIT.
C     EXP(-ELIM).LT.EXP(-ALIM)=EXP(-ELIM)/TOL    AND
C     EXP(ELIM).GT.EXP(ALIM)=EXP(ELIM)*TOL       ARE INTERVALS NEAR
C     UNDERFLOW AND OVERFLOW LIMITS WHERE SCALED ARITHMETIC IS DONE.
C     RL IS THE LOWER BOUNDARY OF THE ASYMPTOTIC EXPANSION FOR LARGE Z.
C     DIG = NUMBER OF BASE 10 DIGITS IN TOL = 10**(-DIG).
C     FNUL IS THE LOWER BOUNDARY OF THE ASYMPTOTIC SERIES FOR LARGE FNU
C-----------------------------------------------------------------------
      TOL = AMAX1(R1MACH(4),1.0E-18)
      K1 = I1MACH(12)
      K2 = I1MACH(13)
      R1M5 = R1MACH(5)
      K = MIN0(IABS(K1),IABS(K2))
      ELIM = 2.303E0*(FLOAT(K)*R1M5-3.0E0)
      K1 = I1MACH(11) - 1
      AA = R1M5*FLOAT(K1)
      DIG = AMIN1(AA,18.0E0)
      AA = AA*2.303E0
      ALIM = ELIM + AMAX1(-AA,-41.45E0)
      FNUL = 10.0E0 + 6.0E0*(DIG-3.0E0)
      RL = 1.2E0*DIG + 3.0E0
      FN = FNU + FLOAT(NN-1)
      MM = 3 - M - M
      FMM = FLOAT(MM)
      ZN = Z*CMPLX(0.0E0,-FMM)
      XN = REAL(ZN)
      YN = AIMAG(ZN)
      AZ = CABS(Z)
C-----------------------------------------------------------------------
C     TEST FOR RANGE
C-----------------------------------------------------------------------
      AA = 0.5E0/TOL
      BB=FLOAT(I1MACH(9))*0.5E0
      AA=AMIN1(AA,BB)
      IF(AZ.GT.AA) GO TO 240
      IF(FN.GT.AA) GO TO 240
      AA=SQRT(AA)
      IF(AZ.GT.AA) IERR=3
      IF(FN.GT.AA) IERR=3
C-----------------------------------------------------------------------
C     OVERFLOW TEST ON THE LAST MEMBER OF THE SEQUENCE
C-----------------------------------------------------------------------
      UFL = R1MACH(1)*1.0E+3
      IF (AZ.LT.UFL) GO TO 220
      IF (FNU.GT.FNUL) GO TO 90
      IF (FN.LE.1.0E0) GO TO 70
      IF (FN.GT.2.0E0) GO TO 60
      IF (AZ.GT.TOL) GO TO 70
      ARG = 0.5E0*AZ
      ALN = -FN*ALOG(ARG)
      IF (ALN.GT.ELIM) GO TO 220
      GO TO 70
   60 CONTINUE
      CALL CUOIK(ZN, FNU, KODE, 2, NN, CY, NUF, TOL, ELIM, ALIM)
      IF (NUF.LT.0) GO TO 220
      NZ = NZ + NUF
      NN = NN - NUF
C-----------------------------------------------------------------------
C     HERE NN=N OR NN=0 SINCE NUF=0,NN, OR -1 ON RETURN FROM CUOIK
C     IF NUF=NN, THEN CY(I)=CZERO FOR ALL I
C-----------------------------------------------------------------------
      IF (NN.EQ.0) GO TO 130
   70 CONTINUE
      IF ((XN.LT.0.0E0) .OR. (XN.EQ.0.0E0 .AND. YN.LT.0.0E0 .AND.
     * M.EQ.2)) GO TO 80
C-----------------------------------------------------------------------
C     RIGHT HALF PLANE COMPUTATION, XN.GE.0. .AND. (XN.NE.0. .OR.
C     YN.GE.0. .OR. M=1)
C-----------------------------------------------------------------------
      CALL CBKNU(ZN, FNU, KODE, NN, CY, NZ, TOL, ELIM, ALIM)
      GO TO 110
C-----------------------------------------------------------------------
C     LEFT HALF PLANE COMPUTATION
C-----------------------------------------------------------------------
   80 CONTINUE
      MR = -MM
      CALL CACON(ZN, FNU, KODE, MR, NN, CY, NW, RL, FNUL, TOL, ELIM,
     * ALIM)
      IF (NW.LT.0) GO TO 230
      NZ=NW
      GO TO 110
   90 CONTINUE
C-----------------------------------------------------------------------
C     UNIFORM ASYMPTOTIC EXPANSIONS FOR FNU.GT.FNUL
C-----------------------------------------------------------------------
      MR = 0
      IF ((XN.GE.0.0E0) .AND. (XN.NE.0.0E0 .OR. YN.GE.0.0E0 .OR.
     * M.NE.2)) GO TO 100
      MR = -MM
      IF (XN.EQ.0.0E0 .AND. YN.LT.0.0E0) ZN = -ZN
  100 CONTINUE
      CALL CBUNK(ZN, FNU, KODE, MR, NN, CY, NW, TOL, ELIM, ALIM)
      IF (NW.LT.0) GO TO 230
      NZ = NZ + NW
  110 CONTINUE
C-----------------------------------------------------------------------
C     H(M,FNU,Z) = -FMM*(I/HPI)*(ZT**FNU)*K(FNU,-Z*ZT)
C
C     ZT=EXP(-FMM*HPI*I) = CMPLX(0.0,-FMM), FMM=3-2*M, M=1,2
C-----------------------------------------------------------------------
      SGN = SIGN(HPI,-FMM)
C-----------------------------------------------------------------------
C     CALCULATE EXP(FNU*HPI*I) TO MINIMIZE LOSSES OF SIGNIFICANCE
C     WHEN FNU IS LARGE
C-----------------------------------------------------------------------
      INU = INT(FNU)
      INUH = INU/2
      IR = INU - 2*INUH
      ARG = (FNU-FLOAT(INU-IR))*SGN
      RHPI = 1.0E0/SGN
      CPN = RHPI*COS(ARG)
      SPN = RHPI*SIN(ARG)
C     ZN = CMPLX(-SPN,CPN)
      CSGN = CMPLX(-SPN,CPN)
C     IF (MOD(INUH,2).EQ.1) ZN = -ZN
      IF (MOD(INUH,2).EQ.1) CSGN = -CSGN
      ZT = CMPLX(0.0E0,-FMM)
      RTOL = 1.0E0/TOL
      ASCLE = UFL*RTOL
      DO 120 I=1,NN
C       CY(I) = CY(I)*ZN
C       ZN = ZN*ZT
        ZN=CY(I)
        AA=REAL(ZN)
        BB=AIMAG(ZN)
        ATOL=1.0E0
        IF (AMAX1(ABS(AA),ABS(BB)).GT.ASCLE) GO TO 125
          ZN = ZN*CMPLX(RTOL,0.0E0)
          ATOL = TOL
  125   CONTINUE
        ZN = ZN*CSGN
        CY(I) = ZN*CMPLX(ATOL,0.0E0)
        CSGN = CSGN*ZT
  120 CONTINUE
      RETURN
  130 CONTINUE
      IF (XN.LT.0.0E0) GO TO 220
      RETURN
  220 CONTINUE
      IERR=2
      NZ=0
      RETURN
  230 CONTINUE
      IF(NW.EQ.(-1)) GO TO 220
      NZ=0
      IERR=5
      RETURN
  240 CONTINUE
      NZ=0
      IERR=4
      RETURN
      END
      SUBROUTINE CBESI(Z, FNU, KODE, N, CY, NZ, IERR)
C***BEGIN PROLOGUE  CBESI
C***DATE WRITTEN   830501   (YYMMDD)
C***REVISION DATE  890801, 930101   (YYMMDD)
C***CATEGORY NO.  B5K
C***KEYWORDS  I-BESSEL FUNCTION,COMPLEX BESSEL FUNCTION,
C             MODIFIED BESSEL FUNCTION OF THE FIRST KIND
C***AUTHOR  AMOS, DONALD E., SANDIA NATIONAL LABORATORIES
C***PURPOSE  TO COMPUTE I-BESSEL FUNCTIONS OF COMPLEX ARGUMENT
C***DESCRIPTION
C
C         ON KODE=1, CBESI COMPUTES AN N MEMBER SEQUENCE OF COMPLEX
C         BESSEL FUNCTIONS CY(J)=I(FNU+J-1,Z) FOR REAL, NONNEGATIVE
C         ORDERS FNU+J-1, J=1,...,N AND COMPLEX Z IN THE CUT PLANE
C         -PI.LT.ARG(Z).LE.PI. ON KODE=2, CBESI RETURNS THE SCALED
C         FUNCTIONS
C
C         CY(J)=EXP(-ABS(X))*I(FNU+J-1,Z)   J = 1,...,N , X=REAL(Z)
C
C         WITH THE EXPONENTIAL GROWTH REMOVED IN BOTH THE LEFT AND
C         RIGHT HALF PLANES FOR Z TO INFINITY. DEFINITIONS AND
C         NOTATION ARE FOUND IN THE NBS HANDBOOK OF MATHEMATICAL
C         FUNCTIONS (REF.1)
C
C         INPUT
C           Z      - Z=CMPLX(X,Y),  -PI.LT.ARG(Z).LE.PI
C           FNU    - ORDER OF INITIAL I FUNCTION, FNU.GE.0.0E0
C           KODE   - A PARAMETER TO INDICATE THE SCALING OPTION
C                    KODE= 1  RETURNS
C                             CY(J)=I(FNU+J-1,Z), J=1,...,N
C                        = 2  RETURNS
C                             CY(J)=I(FNU+J-1,Z)*EXP(-ABS(X)), J=1,...,N
C           N      - NUMBER OF MEMBERS OF THE SEQUENCE, N.GE.1
C
C         OUTPUT
C           CY     - A COMPLEX VECTOR WHOSE FIRST N COMPONENTS CONTAIN
C                    VALUES FOR THE SEQUENCE
C                    CY(J)=I(FNU+J-1,Z)  OR
C                    CY(J)=I(FNU+J-1,Z)*EXP(-ABS(X))  J=1,...,N
C                    DEPENDING ON KODE, X=REAL(Z)
C           NZ     - NUMBER OF COMPONENTS SET TO ZERO DUE TO UNDERFLOW,
C                    NZ= 0   , NORMAL RETURN
C                    NZ.GT.0 , LAST NZ COMPONENTS OF CY SET TO ZERO
C                              DUE TO UNDERFLOW, CY(J)=CMPLX(0.0,0.0),
C                              J = N-NZ+1,...,N
C           IERR   - ERROR FLAG
C                    IERR=0, NORMAL RETURN - COMPUTATION COMPLETED
C                    IERR=1, INPUT ERROR   - NO COMPUTATION
C                    IERR=2, OVERFLOW      - NO COMPUTATION, REAL(Z) TOO
C                            LARGE ON KODE=1
C                    IERR=3, CABS(Z) OR FNU+N-1 LARGE - COMPUTATION DONE
C                            BUT LOSSES OF SIGNIFCANCE BY ARGUMENT
C                            REDUCTION PRODUCE LESS THAN HALF OF MACHINE
C                            ACCURACY
C                    IERR=4, CABS(Z) OR FNU+N-1 TOO LARGE - NO COMPUTA-
C                            TION BECAUSE OF COMPLETE LOSSES OF SIGNIFI-
C                            CANCE BY ARGUMENT REDUCTION
C                    IERR=5, ERROR              - NO COMPUTATION,
C                            ALGORITHM TERMINATION CONDITION NOT MET
C
C***LONG DESCRIPTION
C
C         THE COMPUTATION IS CARRIED OUT BY THE POWER SERIES FOR
C         SMALL CABS(Z), THE ASYMPTOTIC EXPANSION FOR LARGE CABS(Z),
C         THE MILLER ALGORITHM NORMALIZED BY THE WRONSKIAN AND A
C         NEUMANN SERIES FOR IMTERMEDIATE MAGNITUDES, AND THE
C         UNIFORM ASYMPTOTIC EXPANSIONS FOR I(FNU,Z) AND J(FNU,Z)
C         FOR LARGE ORDERS. BACKWARD RECURRENCE IS USED TO GENERATE
C         SEQUENCES OR REDUCE ORDERS WHEN NECESSARY.
C
C         THE CALCULATIONS ABOVE ARE DONE IN THE RIGHT HALF PLANE AND
C         CONTINUED INTO THE LEFT HALF PLANE BY THE FORMULA
C
C         I(FNU,Z*EXP(M*PI)) = EXP(M*PI*FNU)*I(FNU,Z)  REAL(Z).GT.0.0
C                       M = +I OR -I,  I**2=-1
C
C         FOR NEGATIVE ORDERS,THE FORMULA
C
C              I(-FNU,Z) = I(FNU,Z) + (2/PI)*SIN(PI*FNU)*K(FNU,Z)
C
C         CAN BE USED. HOWEVER,FOR LARGE ORDERS CLOSE TO INTEGERS, THE
C         THE FUNCTION CHANGES RADICALLY. WHEN FNU IS A LARGE POSITIVE
C         INTEGER,THE MAGNITUDE OF I(-FNU,Z)=I(FNU,Z) IS A LARGE
C         NEGATIVE POWER OF TEN. BUT WHEN FNU IS NOT AN INTEGER,
C         K(FNU,Z) DOMINATES IN MAGNITUDE WITH A LARGE POSITIVE POWER OF
C         TEN AND THE MOST THAT THE SECOND TERM CAN BE REDUCED IS BY
C         UNIT ROUNDOFF FROM THE COEFFICIENT. THUS, WIDE CHANGES CAN
C         OCCUR WITHIN UNIT ROUNDOFF OF A LARGE INTEGER FOR FNU. HERE,
C         LARGE MEANS FNU.GT.CABS(Z).
C
C         IN MOST COMPLEX VARIABLE COMPUTATION, ONE MUST EVALUATE ELE-
C         MENTARY FUNCTIONS. WHEN THE MAGNITUDE OF Z OR FNU+N-1 IS
C         LARGE, LOSSES OF SIGNIFICANCE BY ARGUMENT REDUCTION OCCUR.
C         CONSEQUENTLY, IF EITHER ONE EXCEEDS U1=SQRT(0.5/UR), THEN
C         LOSSES EXCEEDING HALF PRECISION ARE LIKELY AND AN ERROR FLAG
C         IERR=3 IS TRIGGERED WHERE UR=R1MACH(4)=UNIT ROUNDOFF. ALSO
C         IF EITHER IS LARGER THAN U2=0.5/UR, THEN ALL SIGNIFICANCE IS
C         LOST AND IERR=4. IN ORDER TO USE THE INT FUNCTION, ARGUMENTS
C         MUST BE FURTHER RESTRICTED NOT TO EXCEED THE LARGEST MACHINE
C         INTEGER, U3=I1MACH(9). THUS, THE MAGNITUDE OF Z AND FNU+N-1 IS
C         RESTRICTED BY MIN(U2,U3). ON 32 BIT MACHINES, U1,U2, AND U3
C         ARE APPROXIMATELY 2.0E+3, 4.2E+6, 2.1E+9 IN SINGLE PRECISION
C         ARITHMETIC AND 1.3E+8, 1.8E+16, 2.1E+9 IN DOUBLE PRECISION
C         ARITHMETIC RESPECTIVELY. THIS MAKES U2 AND U3 LIMITING IN
C         THEIR RESPECTIVE ARITHMETICS. THIS MEANS THAT ONE CAN EXPECT
C         TO RETAIN, IN THE WORST CASES ON 32 BIT MACHINES, NO DIGITS
C         IN SINGLE AND ONLY 7 DIGITS IN DOUBLE PRECISION ARITHMETIC.
C         SIMILAR CONSIDERATIONS HOLD FOR OTHER MACHINES.
C
C         THE APPROXIMATE RELATIVE ERROR IN THE MAGNITUDE OF A COMPLEX
C         BESSEL FUNCTION CAN BE EXPRESSED BY P*10**S WHERE P=MAX(UNIT
C         ROUNDOFF,1.0E-18) IS THE NOMINAL PRECISION AND 10**S REPRE-
C         SENTS THE INCREASE IN ERROR DUE TO ARGUMENT REDUCTION IN THE
C         ELEMENTARY FUNCTIONS. HERE, S=MAX(1,ABS(LOG10(CABS(Z))),
C         ABS(LOG10(FNU))) APPROXIMATELY (I.E. S=MAX(1,ABS(EXPONENT OF
C         CABS(Z),ABS(EXPONENT OF FNU)) ). HOWEVER, THE PHASE ANGLE MAY
C         HAVE ONLY ABSOLUTE ACCURACY. THIS IS MOST LIKELY TO OCCUR WHEN
C         ONE COMPONENT (IN ABSOLUTE VALUE) IS LARGER THAN THE OTHER BY
C         SEVERAL ORDERS OF MAGNITUDE. IF ONE COMPONENT IS 10**K LARGER
C         THAN THE OTHER, THEN ONE CAN EXPECT ONLY MAX(ABS(LOG10(P))-K,
C         0) SIGNIFICANT DIGITS; OR, STATED ANOTHER WAY, WHEN K EXCEEDS
C         THE EXPONENT OF P, NO SIGNIFICANT DIGITS REMAIN IN THE SMALLER
C         COMPONENT. HOWEVER, THE PHASE ANGLE RETAINS ABSOLUTE ACCURACY
C         BECAUSE, IN COMPLEX ARITHMETIC WITH PRECISION P, THE SMALLER
C         COMPONENT WILL NOT (AS A RULE) DECREASE BELOW P TIMES THE
C         MAGNITUDE OF THE LARGER COMPONENT. IN THESE EXTREME CASES,
C         THE PRINCIPAL PHASE ANGLE IS ON THE ORDER OF +P, -P, PI/2-P,
C         OR -PI/2+P.
C
C***REFERENCES  HANDBOOK OF MATHEMATICAL FUNCTIONS BY M. ABRAMOWITZ
C                 AND I. A. STEGUN, NBS AMS SERIES 55, U.S. DEPT. OF
C                 COMMERCE, 1955.
C
C               COMPUTATION OF BESSEL FUNCTIONS OF COMPLEX ARGUMENT
C                 BY D. E. AMOS, SAND83-0083, MAY, 1983.
C
C               COMPUTATION OF BESSEL FUNCTIONS OF COMPLEX ARGUMENT
C                 AND LARGE ORDER BY D. E. AMOS, SAND83-0643, MAY, 1983
C
C               A SUBROUTINE PACKAGE FOR BESSEL FUNCTIONS OF A COMPLEX
C                 ARGUMENT AND NONNEGATIVE ORDER BY D. E. AMOS, SAND85-
C                 1018, MAY, 1985
C
C               A PORTABLE PACKAGE FOR BESSEL FUNCTIONS OF A COMPLEX
C                 ARGUMENT AND NONNEGATIVE ORDER BY D. E. AMOS, ACM
C                 TRANS. MATH. SOFTWARE, VOL. 12, NO. 3, SEPTEMBER 1986,
C                 PP 265-273.
C
C***ROUTINES CALLED  CBINU,I1MACH,R1MACH
C***END PROLOGUE  CBESI
      COMPLEX CONE, CSGN, CY, Z, ZN
      REAL AA, ALIM, ARG, DIG, ELIM, FNU, FNUL, PI, RL, R1M5, S1, S2,
     * TOL, XX, YY, R1MACH, AZ, FN, BB, ASCLE, RTOL, ATOL
      INTEGER I, IERR, INU, K, KODE, K1, K2, N, NN, NZ, I1MACH
      DIMENSION CY(N)
      DATA PI /3.14159265358979324E0/
      DATA CONE / (1.0E0,0.0E0) /
C
C***FIRST EXECUTABLE STATEMENT  CBESI
      IERR = 0
      NZ=0
      IF (FNU.LT.0.0E0) IERR=1
      IF (KODE.LT.1 .OR. KODE.GT.2) IERR=1
      IF (N.LT.1) IERR=1
      IF (IERR.NE.0) RETURN
      XX = REAL(Z)
      YY = AIMAG(Z)
C-----------------------------------------------------------------------
C     SET PARAMETERS RELATED TO MACHINE CONSTANTS.
C     TOL IS THE APPROXIMATE UNIT ROUNDOFF LIMITED TO 1.0E-18.
C     ELIM IS THE APPROXIMATE EXPONENTIAL OVER- AND UNDERFLOW LIMIT.
C     EXP(-ELIM).LT.EXP(-ALIM)=EXP(-ELIM)/TOL    AND
C     EXP(ELIM).GT.EXP(ALIM)=EXP(ELIM)*TOL       ARE INTERVALS NEAR
C     UNDERFLOW AND OVERFLOW LIMITS WHERE SCALED ARITHMETIC IS DONE.
C     RL IS THE LOWER BOUNDARY OF THE ASYMPTOTIC EXPANSION FOR LARGE Z.
C     DIG = NUMBER OF BASE 10 DIGITS IN TOL = 10**(-DIG).
C     FNUL IS THE LOWER BOUNDARY OF THE ASYMPTOTIC SERIES FOR LARGE FNU.
C-----------------------------------------------------------------------
      TOL = AMAX1(R1MACH(4),1.0E-18)
      K1 = I1MACH(12)
      K2 = I1MACH(13)
      R1M5 = R1MACH(5)
      K = MIN0(IABS(K1),IABS(K2))
      ELIM = 2.303E0*(FLOAT(K)*R1M5-3.0E0)
      K1 = I1MACH(11) - 1
      AA = R1M5*FLOAT(K1)
      DIG = AMIN1(AA,18.0E0)
      AA = AA*2.303E0
      ALIM = ELIM + AMAX1(-AA,-41.45E0)
      RL = 1.2E0*DIG + 3.0E0
      FNUL = 10.0E0 + 6.0E0*(DIG-3.0E0)
      AZ = CABS(Z)
C-----------------------------------------------------------------------
C     TEST FOR RANGE
C-----------------------------------------------------------------------
      AA = 0.5E0/TOL
      BB=FLOAT(I1MACH(9))*0.5E0
      AA=AMIN1(AA,BB)
      IF(AZ.GT.AA) GO TO 140
      FN=FNU+FLOAT(N-1)
      IF(FN.GT.AA) GO TO 140
      AA=SQRT(AA)
      IF(AZ.GT.AA) IERR=3
      IF(FN.GT.AA) IERR=3
      ZN = Z
      CSGN = CONE
      IF (XX.GE.0.0E0) GO TO 40
      ZN = -Z
C-----------------------------------------------------------------------
C     CALCULATE CSGN=EXP(FNU*PI*I) TO MINIMIZE LOSSES OF SIGNIFICANCE
C     WHEN FNU IS LARGE
C-----------------------------------------------------------------------
      INU = INT(FNU)
      ARG = (FNU-FLOAT(INU))*PI
      IF (YY.LT.0.0E0) ARG = -ARG
      S1 = COS(ARG)
      S2 = SIN(ARG)
      CSGN = CMPLX(S1,S2)
      IF (MOD(INU,2).EQ.1) CSGN = -CSGN
   40 CONTINUE
      CALL CBINU(ZN, FNU, KODE, N, CY, NZ, RL, FNUL, TOL, ELIM, ALIM)
      IF (NZ.LT.0) GO TO 120
      IF (XX.GE.0.0E0) RETURN
C-----------------------------------------------------------------------
C     ANALYTIC CONTINUATION TO THE LEFT HALF PLANE
C-----------------------------------------------------------------------
      NN = N - NZ
      IF (NN.EQ.0) RETURN
      RTOL = 1.0E0/TOL
      ASCLE = R1MACH(1)*RTOL*1.0E+3
      DO 50 I=1,NN
C       CY(I) = CY(I)*CSGN
        ZN=CY(I)
        AA=REAL(ZN)
        BB=AIMAG(ZN)
        ATOL=1.0E0
        IF (AMAX1(ABS(AA),ABS(BB)).GT.ASCLE) GO TO 55
          ZN = ZN*CMPLX(RTOL,0.0E0)
          ATOL = TOL
   55   CONTINUE
        ZN = ZN*CSGN
        CY(I) = ZN*CMPLX(ATOL,0.0E0)
        CSGN = -CSGN
   50 CONTINUE
      RETURN
  120 CONTINUE
      IF(NZ.EQ.(-2)) GO TO 130
      NZ = 0
      IERR=2
      RETURN
  130 CONTINUE
      NZ=0
      IERR=5
      RETURN
  140 CONTINUE
      NZ=0
      IERR=4
      RETURN
      END
      SUBROUTINE CBESJ(Z, FNU, KODE, N, CY, NZ, IERR)
C***BEGIN PROLOGUE  CBESJ
C***DATE WRITTEN   830501   (YYMMDD)
C***REVISION DATE  890801, 930101   (YYMMDD)
C***CATEGORY NO.  B5K
C***KEYWORDS  J-BESSEL FUNCTION,BESSEL FUNCTION OF COMPLEX ARGUMENT,
C             BESSEL FUNCTION OF FIRST KIND
C***AUTHOR  AMOS, DONALD E., SANDIA NATIONAL LABORATORIES
C***PURPOSE  TO COMPUTE THE J-BESSEL FUNCTION OF A COMPLEX ARGUMENT
C***DESCRIPTION
C
C         ON KODE=1, CBESJ COMPUTES AN N MEMBER  SEQUENCE OF COMPLEX
C         BESSEL FUNCTIONS CY(I)=J(FNU+I-1,Z) FOR REAL, NONNEGATIVE
C         ORDERS FNU+I-1, I=1,...,N AND COMPLEX Z IN THE CUT PLANE
C         -PI.LT.ARG(Z).LE.PI. ON KODE=2, CBESJ RETURNS THE SCALED
C         FUNCTIONS
C
C         CY(I)=EXP(-ABS(Y))*J(FNU+I-1,Z)   I = 1,...,N , Y=AIMAG(Z)
C
C         WHICH REMOVE THE EXPONENTIAL GROWTH IN BOTH THE UPPER AND
C         LOWER HALF PLANES FOR Z TO INFINITY. DEFINITIONS AND NOTATION
C         ARE FOUND IN THE NBS HANDBOOK OF MATHEMATICAL FUNCTIONS
C         (REF. 1).
C
C         INPUT
C           Z      - Z=CMPLX(X,Y),  -PI.LT.ARG(Z).LE.PI
C           FNU    - ORDER OF INITIAL J FUNCTION, FNU.GE.0.0E0
C           KODE   - A PARAMETER TO INDICATE THE SCALING OPTION
C                    KODE= 1  RETURNS
C                             CY(I)=J(FNU+I-1,Z), I=1,...,N
C                        = 2  RETURNS
C                             CY(I)=J(FNU+I-1,Z)*EXP(-ABS(Y)), I=1,...
C           N      - NUMBER OF MEMBERS OF THE SEQUENCE, N.GE.1
C
C         OUTPUT
C           CY     - A COMPLEX VECTOR WHOSE FIRST N COMPONENTS CONTAIN
C                    VALUES FOR THE SEQUENCE
C                    CY(I)=J(FNU+I-1,Z)  OR
C                    CY(I)=J(FNU+I-1,Z)*EXP(-ABS(Y))  I=1,...,N
C                    DEPENDING ON KODE, Y=AIMAG(Z).
C           NZ     - NUMBER OF COMPONENTS SET TO ZERO DUE TO UNDERFLOW,
C                    NZ= 0   , NORMAL RETURN
C                    NZ.GT.0 , LAST NZ COMPONENTS OF CY SET TO ZERO
C                              DUE TO UNDERFLOW, CY(I)=CMPLX(0.0,0.0),
C                              I = N-NZ+1,...,N
C           IERR   - ERROR FLAG
C                    IERR=0, NORMAL RETURN - COMPUTATION COMPLETED
C                    IERR=1, INPUT ERROR   - NO COMPUTATION
C                    IERR=2, OVERFLOW      - NO COMPUTATION, AIMAG(Z)
C                            TOO LARGE ON KODE=1
C                    IERR=3, CABS(Z) OR FNU+N-1 LARGE - COMPUTATION DONE
C                            BUT LOSSES OF SIGNIFCANCE BY ARGUMENT
C                            REDUCTION PRODUCE LESS THAN HALF OF MACHINE
C                            ACCURACY
C                    IERR=4, CABS(Z) OR FNU+N-1 TOO LARGE - NO COMPUTA-
C                            TION BECAUSE OF COMPLETE LOSSES OF SIGNIFI-
C                            CANCE BY ARGUMENT REDUCTION
C                    IERR=5, ERROR              - NO COMPUTATION,
C                            ALGORITHM TERMINATION CONDITION NOT MET
C
C***LONG DESCRIPTION
C
C         THE COMPUTATION IS CARRIED OUT BY THE FORMULA
C
C         J(FNU,Z)=EXP( FNU*PI*I/2)*I(FNU,-I*Z)    AIMAG(Z).GE.0.0
C
C         J(FNU,Z)=EXP(-FNU*PI*I/2)*I(FNU, I*Z)    AIMAG(Z).LT.0.0
C
C         WHERE I**2 = -1 AND I(FNU,Z) IS THE I BESSEL FUNCTION.
C
C         FOR NEGATIVE ORDERS,THE FORMULA
C
C              J(-FNU,Z) = J(FNU,Z)*COS(PI*FNU) - Y(FNU,Z)*SIN(PI*FNU)
C
C         CAN BE USED. HOWEVER,FOR LARGE ORDERS CLOSE TO INTEGERS, THE
C         THE FUNCTION CHANGES RADICALLY. WHEN FNU IS A LARGE POSITIVE
C         INTEGER,THE MAGNITUDE OF J(-FNU,Z)=J(FNU,Z)*COS(PI*FNU) IS A
C         LARGE NEGATIVE POWER OF TEN. BUT WHEN FNU IS NOT AN INTEGER,
C         Y(FNU,Z) DOMINATES IN MAGNITUDE WITH A LARGE POSITIVE POWER OF
C         TEN AND THE MOST THAT THE SECOND TERM CAN BE REDUCED IS BY
C         UNIT ROUNDOFF FROM THE COEFFICIENT. THUS, WIDE CHANGES CAN
C         OCCUR WITHIN UNIT ROUNDOFF OF A LARGE INTEGER FOR FNU. HERE,
C         LARGE MEANS FNU.GT.CABS(Z).
C
C         IN MOST COMPLEX VARIABLE COMPUTATION, ONE MUST EVALUATE ELE-
C         MENTARY FUNCTIONS. WHEN THE MAGNITUDE OF Z OR FNU+N-1 IS
C         LARGE, LOSSES OF SIGNIFICANCE BY ARGUMENT REDUCTION OCCUR.
C         CONSEQUENTLY, IF EITHER ONE EXCEEDS U1=SQRT(0.5/UR), THEN
C         LOSSES EXCEEDING HALF PRECISION ARE LIKELY AND AN ERROR FLAG
C         IERR=3 IS TRIGGERED WHERE UR=R1MACH(4)=UNIT ROUNDOFF. ALSO
C         IF EITHER IS LARGER THAN U2=0.5/UR, THEN ALL SIGNIFICANCE IS
C         LOST AND IERR=4. IN ORDER TO USE THE INT FUNCTION, ARGUMENTS
C         MUST BE FURTHER RESTRICTED NOT TO EXCEED THE LARGEST MACHINE
C         INTEGER, U3=I1MACH(9). THUS, THE MAGNITUDE OF Z AND FNU+N-1 IS
C         RESTRICTED BY MIN(U2,U3). ON 32 BIT MACHINES, U1,U2, AND U3
C         ARE APPROXIMATELY 2.0E+3, 4.2E+6, 2.1E+9 IN SINGLE PRECISION
C         ARITHMETIC AND 1.3E+8, 1.8E+16, 2.1E+9 IN DOUBLE PRECISION
C         ARITHMETIC RESPECTIVELY. THIS MAKES U2 AND U3 LIMITING IN
C         THEIR RESPECTIVE ARITHMETICS. THIS MEANS THAT ONE CAN EXPECT
C         TO RETAIN, IN THE WORST CASES ON 32 BIT MACHINES, NO DIGITS
C         IN SINGLE AND ONLY 7 DIGITS IN DOUBLE PRECISION ARITHMETIC.
C         SIMILAR CONSIDERATIONS HOLD FOR OTHER MACHINES.
C
C         THE APPROXIMATE RELATIVE ERROR IN THE MAGNITUDE OF A COMPLEX
C         BESSEL FUNCTION CAN BE EXPRESSED BY P*10**S WHERE P=MAX(UNIT
C         ROUNDOFF,1.0E-18) IS THE NOMINAL PRECISION AND 10**S REPRE-
C         SENTS THE INCREASE IN ERROR DUE TO ARGUMENT REDUCTION IN THE
C         ELEMENTARY FUNCTIONS. HERE, S=MAX(1,ABS(LOG10(CABS(Z))),
C         ABS(LOG10(FNU))) APPROXIMATELY (I.E. S=MAX(1,ABS(EXPONENT OF
C         CABS(Z),ABS(EXPONENT OF FNU)) ). HOWEVER, THE PHASE ANGLE MAY
C         HAVE ONLY ABSOLUTE ACCURACY. THIS IS MOST LIKELY TO OCCUR WHEN
C         ONE COMPONENT (IN ABSOLUTE VALUE) IS LARGER THAN THE OTHER BY
C         SEVERAL ORDERS OF MAGNITUDE. IF ONE COMPONENT IS 10**K LARGER
C         THAN THE OTHER, THEN ONE CAN EXPECT ONLY MAX(ABS(LOG10(P))-K,
C         0) SIGNIFICANT DIGITS; OR, STATED ANOTHER WAY, WHEN K EXCEEDS
C         THE EXPONENT OF P, NO SIGNIFICANT DIGITS REMAIN IN THE SMALLER
C         COMPONENT. HOWEVER, THE PHASE ANGLE RETAINS ABSOLUTE ACCURACY
C         BECAUSE, IN COMPLEX ARITHMETIC WITH PRECISION P, THE SMALLER
C         COMPONENT WILL NOT (AS A RULE) DECREASE BELOW P TIMES THE
C         MAGNITUDE OF THE LARGER COMPONENT. IN THESE EXTREME CASES,
C         THE PRINCIPAL PHASE ANGLE IS ON THE ORDER OF +P, -P, PI/2-P,
C         OR -PI/2+P.
C
C***REFERENCES  HANDBOOK OF MATHEMATICAL FUNCTIONS BY M. ABRAMOWITZ
C                 AND I. A. STEGUN, NBS AMS SERIES 55, U.S. DEPT. OF
C                 COMMERCE, 1955.
C
C               COMPUTATION OF BESSEL FUNCTIONS OF COMPLEX ARGUMENT
C                 BY D. E. AMOS, SAND83-0083, MAY, 1983.
C
C               COMPUTATION OF BESSEL FUNCTIONS OF COMPLEX ARGUMENT
C                 AND LARGE ORDER BY D. E. AMOS, SAND83-0643, MAY, 1983
C
C               A SUBROUTINE PACKAGE FOR BESSEL FUNCTIONS OF A COMPLEX
C                 ARGUMENT AND NONNEGATIVE ORDER BY D. E. AMOS, SAND85-
C                 1018, MAY, 1985
C
C               A PORTABLE PACKAGE FOR BESSEL FUNCTIONS OF A COMPLEX
C                 ARGUMENT AND NONNEGATIVE ORDER BY D. E. AMOS, ACM
C                 TRANS. MATH. SOFTWARE, VOL. 12, NO. 3, SEPTEMBER 1986,
C                 PP 265-273.
C
C***ROUTINES CALLED  CBINU,I1MACH,R1MACH
C***END PROLOGUE  CBESJ
C
      COMPLEX CI, CSGN, CY, Z, ZN
      REAL AA, ALIM, ARG, DIG, ELIM, FNU, FNUL, HPI, RL, R1, R1M5, R2,
     * TOL, YY, R1MACH, AZ, FN, BB, ASCLE, RTOL, ATOL
      INTEGER I, IERR, INU, INUH, IR, KODE, K1, K2, N, NL, NZ, I1MACH, K
      DIMENSION CY(N)
      DATA HPI /1.57079632679489662E0/
C
C***FIRST EXECUTABLE STATEMENT  CBESJ
      IERR = 0
      NZ=0
      IF (FNU.LT.0.0E0) IERR=1
      IF (KODE.LT.1 .OR. KODE.GT.2) IERR=1
      IF (N.LT.1) IERR=1
      IF (IERR.NE.0) RETURN
C-----------------------------------------------------------------------
C     SET PARAMETERS RELATED TO MACHINE CONSTANTS.
C     TOL IS THE APPROXIMATE UNIT ROUNDOFF LIMITED TO 1.0E-18.
C     ELIM IS THE APPROXIMATE EXPONENTIAL OVER- AND UNDERFLOW LIMIT.
C     EXP(-ELIM).LT.EXP(-ALIM)=EXP(-ELIM)/TOL    AND
C     EXP(ELIM).GT.EXP(ALIM)=EXP(ELIM)*TOL       ARE INTERVALS NEAR
C     UNDERFLOW AND OVERFLOW LIMITS WHERE SCALED ARITHMETIC IS DONE.
C     RL IS THE LOWER BOUNDARY OF THE ASYMPTOTIC EXPANSION FOR LARGE Z.
C     DIG = NUMBER OF BASE 10 DIGITS IN TOL = 10**(-DIG).
C     FNUL IS THE LOWER BOUNDARY OF THE ASYMPTOTIC SERIES FOR LARGE FNU.
C-----------------------------------------------------------------------
      TOL = AMAX1(R1MACH(4),1.0E-18)
      K1 = I1MACH(12)
      K2 = I1MACH(13)
      R1M5 = R1MACH(5)
      K = MIN0(IABS(K1),IABS(K2))
      ELIM = 2.303E0*(FLOAT(K)*R1M5-3.0E0)
      K1 = I1MACH(11) - 1
      AA = R1M5*FLOAT(K1)
      DIG = AMIN1(AA,18.0E0)
      AA = AA*2.303E0
      ALIM = ELIM + AMAX1(-AA,-41.45E0)
      RL = 1.2E0*DIG + 3.0E0
      FNUL = 10.0E0 + 6.0E0*(DIG-3.0E0)
      CI = CMPLX(0.0E0,1.0E0)
      YY = AIMAG(Z)
      AZ = CABS(Z)
C-----------------------------------------------------------------------
C     TEST FOR RANGE
C-----------------------------------------------------------------------
      AA = 0.5E0/TOL
      BB=FLOAT(I1MACH(9))*0.5E0
      AA=AMIN1(AA,BB)
      FN=FNU+FLOAT(N-1)
      IF(AZ.GT.AA) GO TO 140
      IF(FN.GT.AA) GO TO 140
      AA=SQRT(AA)
      IF(AZ.GT.AA) IERR=3
      IF(FN.GT.AA) IERR=3
C-----------------------------------------------------------------------
C     CALCULATE CSGN=EXP(FNU*HPI*I) TO MINIMIZE LOSSES OF SIGNIFICANCE
C     WHEN FNU IS LARGE
C-----------------------------------------------------------------------
      INU = INT(FNU)
      INUH = INU/2
      IR = INU - 2*INUH
      ARG = (FNU-FLOAT(INU-IR))*HPI
      R1 = COS(ARG)
      R2 = SIN(ARG)
      CSGN = CMPLX(R1,R2)
      IF (MOD(INUH,2).EQ.1) CSGN = -CSGN
C-----------------------------------------------------------------------
C     ZN IS IN THE RIGHT HALF PLANE
C-----------------------------------------------------------------------
      ZN = -Z*CI
      IF (YY.GE.0.0E0) GO TO 40
      ZN = -ZN
      CSGN = CONJG(CSGN)
      CI = CONJG(CI)
   40 CONTINUE
      CALL CBINU(ZN, FNU, KODE, N, CY, NZ, RL, FNUL, TOL, ELIM, ALIM)
      IF (NZ.LT.0) GO TO 120
      NL = N - NZ
      IF (NL.EQ.0) RETURN
      RTOL = 1.0E0/TOL
      ASCLE = R1MACH(1)*RTOL*1.0E+3
      DO 50 I=1,NL
C       CY(I)=CY(I)*CSGN
        ZN=CY(I)
        AA=REAL(ZN)
        BB=AIMAG(ZN)
        ATOL=1.0E0
        IF (AMAX1(ABS(AA),ABS(BB)).GT.ASCLE) GO TO 55
          ZN = ZN*CMPLX(RTOL,0.0E0)
          ATOL = TOL
   55   CONTINUE
        ZN = ZN*CSGN
        CY(I) = ZN*CMPLX(ATOL,0.0E0)
        CSGN = CSGN*CI
   50 CONTINUE
      RETURN
  120 CONTINUE
      IF(NZ.EQ.(-2)) GO TO 130
      NZ = 0
      IERR = 2
      RETURN
  130 CONTINUE
      NZ=0
      IERR=5
      RETURN
  140 CONTINUE
      NZ=0
      IERR=4
      RETURN
      END
      SUBROUTINE CBESK(Z, FNU, KODE, N, CY, NZ, IERR)
C***BEGIN PROLOGUE  CBESK
C***DATE WRITTEN   830501   (YYMMDD)
C***REVISION DATE  890801, 930101   (YYMMDD)
C***CATEGORY NO.  B5K
C***KEYWORDS  K-BESSEL FUNCTION,COMPLEX BESSEL FUNCTION,
C             MODIFIED BESSEL FUNCTION OF THE SECOND KIND,
C             BESSEL FUNCTION OF THE THIRD KIND
C***AUTHOR  AMOS, DONALD E., SANDIA NATIONAL LABORATORIES
C***PURPOSE  TO COMPUTE K-BESSEL FUNCTIONS OF COMPLEX ARGUMENT
C***DESCRIPTION
C
C         ON KODE=1, CBESK COMPUTES AN N MEMBER SEQUENCE OF COMPLEX
C         BESSEL FUNCTIONS CY(J)=K(FNU+J-1,Z) FOR REAL, NONNEGATIVE
C         ORDERS FNU+J-1, J=1,...,N AND COMPLEX Z.NE.CMPLX(0.0,0.0)
C         IN THE CUT PLANE -PI.LT.ARG(Z).LE.PI. ON KODE=2, CBESK
C         RETURNS THE SCALED K FUNCTIONS,
C
C         CY(J)=EXP(Z)*K(FNU+J-1,Z) , J=1,...,N,
C
C         WHICH REMOVE THE EXPONENTIAL BEHAVIOR IN BOTH THE LEFT AND
C         RIGHT HALF PLANES FOR Z TO INFINITY. DEFINITIONS AND
C         NOTATION ARE FOUND IN THE NBS HANDBOOK OF MATHEMATICAL
C         FUNCTIONS (REF. 1).
C
C         INPUT
C           Z      - Z=CMPLX(X,Y),Z.NE.CMPLX(0.,0.),-PI.LT.ARG(Z).LE.PI
C           FNU    - ORDER OF INITIAL K FUNCTION, FNU.GE.0.0E0
C           N      - NUMBER OF MEMBERS OF THE SEQUENCE, N.GE.1
C           KODE   - A PARAMETER TO INDICATE THE SCALING OPTION
C                    KODE= 1  RETURNS
C                             CY(I)=K(FNU+I-1,Z), I=1,...,N
C                        = 2  RETURNS
C                             CY(I)=K(FNU+I-1,Z)*EXP(Z), I=1,...,N
C
C         OUTPUT
C           CY     - A COMPLEX VECTOR WHOSE FIRST N COMPONENTS CONTAIN
C                    VALUES FOR THE SEQUENCE
C                    CY(I)=K(FNU+I-1,Z), I=1,...,N OR
C                    CY(I)=K(FNU+I-1,Z)*EXP(Z), I=1,...,N
C                    DEPENDING ON KODE
C           NZ     - NUMBER OF COMPONENTS SET TO ZERO DUE TO UNDERFLOW.
C                    NZ= 0   , NORMAL RETURN
C                    NZ.GT.0 , FIRST NZ COMPONENTS OF CY SET TO ZERO
C                              DUE TO UNDERFLOW, CY(I)=CMPLX(0.0,0.0),
C                              I=1,...,N WHEN X.GE.0.0. WHEN X.LT.0.0
C                              NZ STATES ONLY THE NUMBER OF UNDERFLOWS
C                              IN THE SEQUENCE.
C           IERR   - ERROR FLAG
C                    IERR=0, NORMAL RETURN - COMPUTATION COMPLETED
C                    IERR=1, INPUT ERROR   - NO COMPUTATION
C                    IERR=2, OVERFLOW      - NO COMPUTATION, FNU+N-1 IS
C                            TOO LARGE OR CABS(Z) IS TOO SMALL OR BOTH
C                    IERR=3, CABS(Z) OR FNU+N-1 LARGE - COMPUTATION DONE
C                            BUT LOSSES OF SIGNIFCANCE BY ARGUMENT
C                            REDUCTION PRODUCE LESS THAN HALF OF MACHINE
C                            ACCURACY
C                    IERR=4, CABS(Z) OR FNU+N-1 TOO LARGE - NO COMPUTA-
C                            TION BECAUSE OF COMPLETE LOSSES OF SIGNIFI-
C                            CANCE BY ARGUMENT REDUCTION
C                    IERR=5, ERROR              - NO COMPUTATION,
C                            ALGORITHM TERMINATION CONDITION NOT MET
C
C***LONG DESCRIPTION
C
C         EQUATIONS OF THE REFERENCE ARE IMPLEMENTED FOR SMALL ORDERS
C         DNU AND DNU+1.0 IN THE RIGHT HALF PLANE X.GE.0.0. FORWARD
C         RECURRENCE GENERATES HIGHER ORDERS. K IS CONTINUED TO THE LEFT
C         HALF PLANE BY THE RELATION
C
C         K(FNU,Z*EXP(MP)) = EXP(-MP*FNU)*K(FNU,Z)-MP*I(FNU,Z)
C         MP=MR*PI*I, MR=+1 OR -1, RE(Z).GT.0, I**2=-1
C
C         WHERE I(FNU,Z) IS THE I BESSEL FUNCTION.
C
C         FOR LARGE ORDERS, FNU.GT.FNUL, THE K FUNCTION IS COMPUTED
C         BY MEANS OF ITS UNIFORM ASYMPTOTIC EXPANSIONS.
C
C         FOR NEGATIVE ORDERS, THE FORMULA
C
C                       K(-FNU,Z) = K(FNU,Z)
C
C         CAN BE USED.
C
C         CBESK ASSUMES THAT A SIGNIFICANT DIGIT SINH(X) FUNCTION IS
C         AVAILABLE.
C
C         IN MOST COMPLEX VARIABLE COMPUTATION, ONE MUST EVALUATE ELE-
C         MENTARY FUNCTIONS. WHEN THE MAGNITUDE OF Z OR FNU+N-1 IS
C         LARGE, LOSSES OF SIGNIFICANCE BY ARGUMENT REDUCTION OCCUR.
C         CONSEQUENTLY, IF EITHER ONE EXCEEDS U1=SQRT(0.5/UR), THEN
C         LOSSES EXCEEDING HALF PRECISION ARE LIKELY AND AN ERROR FLAG
C         IERR=3 IS TRIGGERED WHERE UR=R1MACH(4)=UNIT ROUNDOFF. ALSO
C         IF EITHER IS LARGER THAN U2=0.5/UR, THEN ALL SIGNIFICANCE IS
C         LOST AND IERR=4. IN ORDER TO USE THE INT FUNCTION, ARGUMENTS
C         MUST BE FURTHER RESTRICTED NOT TO EXCEED THE LARGEST MACHINE
C         INTEGER, U3=I1MACH(9). THUS, THE MAGNITUDE OF Z AND FNU+N-1 IS
C         RESTRICTED BY MIN(U2,U3). ON 32 BIT MACHINES, U1,U2, AND U3
C         ARE APPROXIMATELY 2.0E+3, 4.2E+6, 2.1E+9 IN SINGLE PRECISION
C         ARITHMETIC AND 1.3E+8, 1.8E+16, 2.1E+9 IN DOUBLE PRECISION
C         ARITHMETIC RESPECTIVELY. THIS MAKES U2 AND U3 LIMITING IN
C         THEIR RESPECTIVE ARITHMETICS. THIS MEANS THAT ONE CAN EXPECT
C         TO RETAIN, IN THE WORST CASES ON 32 BIT MACHINES, NO DIGITS
C         IN SINGLE AND ONLY 7 DIGITS IN DOUBLE PRECISION ARITHMETIC.
C         SIMILAR CONSIDERATIONS HOLD FOR OTHER MACHINES.
C
C         THE APPROXIMATE RELATIVE ERROR IN THE MAGNITUDE OF A COMPLEX
C         BESSEL FUNCTION CAN BE EXPRESSED BY P*10**S WHERE P=MAX(UNIT
C         ROUNDOFF,1.0E-18) IS THE NOMINAL PRECISION AND 10**S REPRE-
C         SENTS THE INCREASE IN ERROR DUE TO ARGUMENT REDUCTION IN THE
C         ELEMENTARY FUNCTIONS. HERE, S=MAX(1,ABS(LOG10(CABS(Z))),
C         ABS(LOG10(FNU))) APPROXIMATELY (I.E. S=MAX(1,ABS(EXPONENT OF
C         CABS(Z),ABS(EXPONENT OF FNU)) ). HOWEVER, THE PHASE ANGLE MAY
C         HAVE ONLY ABSOLUTE ACCURACY. THIS IS MOST LIKELY TO OCCUR WHEN
C         ONE COMPONENT (IN ABSOLUTE VALUE) IS LARGER THAN THE OTHER BY
C         SEVERAL ORDERS OF MAGNITUDE. IF ONE COMPONENT IS 10**K LARGER
C         THAN THE OTHER, THEN ONE CAN EXPECT ONLY MAX(ABS(LOG10(P))-K,
C         0) SIGNIFICANT DIGITS; OR, STATED ANOTHER WAY, WHEN K EXCEEDS
C         THE EXPONENT OF P, NO SIGNIFICANT DIGITS REMAIN IN THE SMALLER
C         COMPONENT. HOWEVER, THE PHASE ANGLE RETAINS ABSOLUTE ACCURACY
C         BECAUSE, IN COMPLEX ARITHMETIC WITH PRECISION P, THE SMALLER
C         COMPONENT WILL NOT (AS A RULE) DECREASE BELOW P TIMES THE
C         MAGNITUDE OF THE LARGER COMPONENT. IN THESE EXTREME CASES,
C         THE PRINCIPAL PHASE ANGLE IS ON THE ORDER OF +P, -P, PI/2-P,
C         OR -PI/2+P.
C
C***REFERENCES  HANDBOOK OF MATHEMATICAL FUNCTIONS BY M. ABRAMOWITZ
C                 AND I. A. STEGUN, NBS AMS SERIES 55, U.S. DEPT. OF
C                 COMMERCE, 1955.
C
C               COMPUTATION OF BESSEL FUNCTIONS OF COMPLEX ARGUMENT
C                 BY D. E. AMOS, SAND83-0083, MAY, 1983.
C
C               COMPUTATION OF BESSEL FUNCTIONS OF COMPLEX ARGUMENT
C                 AND LARGE ORDER BY D. E. AMOS, SAND83-0643, MAY, 1983.
C
C               A SUBROUTINE PACKAGE FOR BESSEL FUNCTIONS OF A COMPLEX
C                 ARGUMENT AND NONNEGATIVE ORDER BY D. E. AMOS, SAND85-
C                 1018, MAY, 1985
C
C               A PORTABLE PACKAGE FOR BESSEL FUNCTIONS OF A COMPLEX
C                 ARGUMENT AND NONNEGATIVE ORDER BY D. E. AMOS, ACM
C                 TRANS. MATH. SOFTWARE, VOL. 12, NO. 3, SEPTEMBER 1986,
C                 PP 265-273.
C
C***ROUTINES CALLED  CACON,CBKNU,CBUNK,CUOIK,I1MACH,R1MACH
C***END PROLOGUE  CBESK
C
      COMPLEX CY, Z
      REAL AA, ALIM, ALN, ARG, AZ, DIG, ELIM, FN, FNU, FNUL, RL, R1M5,
     * TOL, UFL, XX, YY, R1MACH, BB
      INTEGER IERR, K, KODE, K1, K2, MR, N, NN, NUF, NW, NZ, I1MACH
      DIMENSION CY(N)
C***FIRST EXECUTABLE STATEMENT  CBESK
      IERR = 0
      NZ=0
      XX = REAL(Z)
      YY = AIMAG(Z)
      IF (YY.EQ.0.0E0 .AND. XX.EQ.0.0E0) IERR=1
      IF (FNU.LT.0.0E0) IERR=1
      IF (KODE.LT.1 .OR. KODE.GT.2) IERR=1
      IF (N.LT.1) IERR=1
      IF (IERR.NE.0) RETURN
      NN = N
C-----------------------------------------------------------------------
C     SET PARAMETERS RELATED TO MACHINE CONSTANTS.
C     TOL IS THE APPROXIMATE UNIT ROUNDOFF LIMITED TO 1.0E-18.
C     ELIM IS THE APPROXIMATE EXPONENTIAL OVER- AND UNDERFLOW LIMIT.
C     EXP(-ELIM).LT.EXP(-ALIM)=EXP(-ELIM)/TOL    AND
C     EXP(ELIM).GT.EXP(ALIM)=EXP(ELIM)*TOL       ARE INTERVALS NEAR
C     UNDERFLOW AND OVERFLOW LIMITS WHERE SCALED ARITHMETIC IS DONE.
C     RL IS THE LOWER BOUNDARY OF THE ASYMPTOTIC EXPANSION FOR LARGE Z.
C     DIG = NUMBER OF BASE 10 DIGITS IN TOL = 10**(-DIG).
C     FNUL IS THE LOWER BOUNDARY OF THE ASYMPTOTIC SERIES FOR LARGE FNU
C-----------------------------------------------------------------------
      TOL = AMAX1(R1MACH(4),1.0E-18)
      K1 = I1MACH(12)
      K2 = I1MACH(13)
      R1M5 = R1MACH(5)
      K = MIN0(IABS(K1),IABS(K2))
      ELIM = 2.303E0*(FLOAT(K)*R1M5-3.0E0)
      K1 = I1MACH(11) - 1
      AA = R1M5*FLOAT(K1)
      DIG = AMIN1(AA,18.0E0)
      AA = AA*2.303E0
      ALIM = ELIM + AMAX1(-AA,-41.45E0)
      FNUL = 10.0E0 + 6.0E0*(DIG-3.0E0)
      RL = 1.2E0*DIG + 3.0E0
      AZ = CABS(Z)
      FN = FNU + FLOAT(NN-1)
C-----------------------------------------------------------------------
C     TEST FOR RANGE
C-----------------------------------------------------------------------
      AA = 0.5E0/TOL
      BB=FLOAT(I1MACH(9))*0.5E0
      AA=AMIN1(AA,BB)
      IF(AZ.GT.AA) GO TO 210
      IF(FN.GT.AA) GO TO 210
      AA=SQRT(AA)
      IF(AZ.GT.AA) IERR=3
      IF(FN.GT.AA) IERR=3
C-----------------------------------------------------------------------
C     OVERFLOW TEST ON THE LAST MEMBER OF THE SEQUENCE
C-----------------------------------------------------------------------
C     UFL = EXP(-ELIM)
      UFL = R1MACH(1)*1.0E+3
      IF (AZ.LT.UFL) GO TO 180
      IF (FNU.GT.FNUL) GO TO 80
      IF (FN.LE.1.0E0) GO TO 60
      IF (FN.GT.2.0E0) GO TO 50
      IF (AZ.GT.TOL) GO TO 60
      ARG = 0.5E0*AZ
      ALN = -FN*ALOG(ARG)
      IF (ALN.GT.ELIM) GO TO 180
      GO TO 60
   50 CONTINUE
      CALL CUOIK(Z, FNU, KODE, 2, NN, CY, NUF, TOL, ELIM, ALIM)
      IF (NUF.LT.0) GO TO 180
      NZ = NZ + NUF
      NN = NN - NUF
C-----------------------------------------------------------------------
C     HERE NN=N OR NN=0 SINCE NUF=0,NN, OR -1 ON RETURN FROM CUOIK
C     IF NUF=NN, THEN CY(I)=CZERO FOR ALL I
C-----------------------------------------------------------------------
      IF (NN.EQ.0) GO TO 100
   60 CONTINUE
      IF (XX.LT.0.0E0) GO TO 70
C-----------------------------------------------------------------------
C     RIGHT HALF PLANE COMPUTATION, REAL(Z).GE.0.
C-----------------------------------------------------------------------
      CALL CBKNU(Z, FNU, KODE, NN, CY, NW, TOL, ELIM, ALIM)
      IF (NW.LT.0) GO TO 200
      NZ=NW
      RETURN
C-----------------------------------------------------------------------
C     LEFT HALF PLANE COMPUTATION
C     PI/2.LT.ARG(Z).LE.PI AND -PI.LT.ARG(Z).LT.-PI/2.
C-----------------------------------------------------------------------
   70 CONTINUE
      IF (NZ.NE.0) GO TO 180
      MR = 1
      IF (YY.LT.0.0E0) MR = -1
      CALL CACON(Z, FNU, KODE, MR, NN, CY, NW, RL, FNUL, TOL, ELIM,
     * ALIM)
      IF (NW.LT.0) GO TO 200
      NZ=NW
      RETURN
C-----------------------------------------------------------------------
C     UNIFORM ASYMPTOTIC EXPANSIONS FOR FNU.GT.FNUL
C-----------------------------------------------------------------------
   80 CONTINUE
      MR = 0
      IF (XX.GE.0.0E0) GO TO 90
      MR = 1
      IF (YY.LT.0.0E0) MR = -1
   90 CONTINUE
      CALL CBUNK(Z, FNU, KODE, MR, NN, CY, NW, TOL, ELIM, ALIM)
      IF (NW.LT.0) GO TO 200
      NZ = NZ + NW
      RETURN
  100 CONTINUE
      IF (XX.LT.0.0E0) GO TO 180
      RETURN
  180 CONTINUE
      NZ = 0
      IERR=2
      RETURN
  200 CONTINUE
      IF(NW.EQ.(-1)) GO TO 180
      NZ=0
      IERR=5
      RETURN
  210 CONTINUE
      NZ=0
      IERR=4
      RETURN
      END
      SUBROUTINE CBESY(Z, FNU, KODE, N, CY, NZ, CWRK, IERR)
C***BEGIN PROLOGUE  CBESY
C***DATE WRITTEN   830501   (YYMMDD)
C***REVISION DATE  890801, 930101  (YYMMDD)
C***CATEGORY NO.  B5K
C***KEYWORDS  Y-BESSEL FUNCTION,BESSEL FUNCTION OF COMPLEX ARGUMENT,
C             BESSEL FUNCTION OF SECOND KIND
C***AUTHOR  AMOS, DONALD E., SANDIA NATIONAL LABORATORIES
C***PURPOSE  TO COMPUTE THE Y-BESSEL FUNCTION OF A COMPLEX ARGUMENT
C***DESCRIPTION
C
C         ON KODE=1, CBESY COMPUTES AN N MEMBER SEQUENCE OF COMPLEX
C         BESSEL FUNCTIONS CY(I)=Y(FNU+I-1,Z) FOR REAL, NONNEGATIVE
C         ORDERS FNU+I-1, I=1,...,N AND COMPLEX Z IN THE CUT PLANE
C         -PI.LT.ARG(Z).LE.PI. ON KODE=2, CBESY RETURNS THE SCALED
C         FUNCTIONS
C
C         CY(I)=EXP(-ABS(Y))*Y(FNU+I-1,Z)   I = 1,...,N , Y=AIMAG(Z)
C
C         WHICH REMOVE THE EXPONENTIAL GROWTH IN BOTH THE UPPER AND
C         LOWER HALF PLANES FOR Z TO INFINITY. DEFINITIONS AND NOTATION
C         ARE FOUND IN THE NBS HANDBOOK OF MATHEMATICAL FUNCTIONS
C         (REF. 1).
C
C         INPUT
C           Z      - Z=CMPLX(X,Y), Z.NE.CMPLX(0.,0.),-PI.LT.ARG(Z).LE.PI
C           FNU    - ORDER OF INITIAL Y FUNCTION, FNU.GE.0.0E0
C           KODE   - A PARAMETER TO INDICATE THE SCALING OPTION
C                    KODE= 1  RETURNS
C                             CY(I)=Y(FNU+I-1,Z), I=1,...,N
C                        = 2  RETURNS
C                             CY(I)=Y(FNU+I-1,Z)*EXP(-ABS(Y)), I=1,...,N
C                             WHERE Y=AIMAG(Z)
C           N      - NUMBER OF MEMBERS OF THE SEQUENCE, N.GE.1
C           CWRK   - A COMPLEX WORK VECTOR OF DIMENSION AT LEAST N
C
C         OUTPUT
C           CY     - A COMPLEX VECTOR WHOSE FIRST N COMPONENTS CONTAIN
C                    VALUES FOR THE SEQUENCE
C                    CY(I)=Y(FNU+I-1,Z)  OR
C                    CY(I)=Y(FNU+I-1,Z)*EXP(-ABS(Y))  I=1,...,N
C                    DEPENDING ON KODE.
C           NZ     - NZ=0 , A NORMAL RETURN
C                    NZ.GT.0 , NZ COMPONENTS OF CY SET TO ZERO DUE TO
C                    UNDERFLOW (GENERALLY ON KODE=2)
C           IERR   - ERROR FLAG
C                    IERR=0, NORMAL RETURN - COMPUTATION COMPLETED
C                    IERR=1, INPUT ERROR   - NO COMPUTATION
C                    IERR=2, OVERFLOW      - NO COMPUTATION, FNU+N-1 IS
C                            TOO LARGE OR CABS(Z) IS TOO SMALL OR BOTH
C                    IERR=3, CABS(Z) OR FNU+N-1 LARGE - COMPUTATION DONE
C                            BUT LOSSES OF SIGNIFCANCE BY ARGUMENT
C                            REDUCTION PRODUCE LESS THAN HALF OF MACHINE
C                            ACCURACY
C                    IERR=4, CABS(Z) OR FNU+N-1 TOO LARGE - NO COMPUTA-
C                            TION BECAUSE OF COMPLETE LOSSES OF SIGNIFI-
C                            CANCE BY ARGUMENT REDUCTION
C                    IERR=5, ERROR              - NO COMPUTATION,
C                            ALGORITHM TERMINATION CONDITION NOT MET
C
C***LONG DESCRIPTION
C
C         THE COMPUTATION IS CARRIED OUT IN TERMS OF THE I(FNU,Z) AND
C         K(FNU,Z) BESSEL FUNCTIONS IN THE RIGHT HALF PLANE BY
C
C             Y(FNU,Z) = I*CC*I(FNU,ARG) - (2/PI)*CONJG(CC)*K(FNU,ARG)
C
C             Y(FNU,Z) = CONJG(Y(FNU,CONJG(Z)))
C
C         FOR AIMAG(Z).GE.0 AND AIMAG(Z).LT.0 RESPECTIVELY, WHERE
C         CC=EXP(I*PI*FNU/2), ARG=Z*EXP(-I*PI/2) AND I**2=-1.
C
C         FOR NEGATIVE ORDERS,THE FORMULA
C
C             Y(-FNU,Z) = Y(FNU,Z)*COS(PI*FNU) + J(FNU,Z)*SIN(PI*FNU)
C
C         CAN BE USED. HOWEVER,FOR LARGE ORDERS CLOSE TO HALF ODD
C         INTEGERS THE FUNCTION CHANGES RADICALLY. WHEN FNU IS A LARGE
C         POSITIVE HALF ODD INTEGER,THE MAGNITUDE OF Y(-FNU,Z)=J(FNU,Z)*
C         SIN(PI*FNU) IS A LARGE NEGATIVE POWER OF TEN. BUT WHEN FNU IS
C         NOT A HALF ODD INTEGER, Y(FNU,Z) DOMINATES IN MAGNITUDE WITH A
C         LARGE POSITIVE POWER OF TEN AND THE MOST THAT THE SECOND TERM
C         CAN BE REDUCED IS BY UNIT ROUNDOFF FROM THE COEFFICIENT. THUS,
C         WIDE CHANGES CAN OCCUR WITHIN UNIT ROUNDOFF OF A LARGE HALF
C         ODD INTEGER. HERE, LARGE MEANS FNU.GT.CABS(Z).
C
C         IN MOST COMPLEX VARIABLE COMPUTATION, ONE MUST EVALUATE ELE-
C         MENTARY FUNCTIONS. WHEN THE MAGNITUDE OF Z OR FNU+N-1 IS
C         LARGE, LOSSES OF SIGNIFICANCE BY ARGUMENT REDUCTION OCCUR.
C         CONSEQUENTLY, IF EITHER ONE EXCEEDS U1=SQRT(0.5/UR), THEN
C         LOSSES EXCEEDING HALF PRECISION ARE LIKELY AND AN ERROR FLAG
C         IERR=3 IS TRIGGERED WHERE UR=R1MACH(4)=UNIT ROUNDOFF. ALSO
C         IF EITHER IS LARGER THAN U2=0.5/UR, THEN ALL SIGNIFICANCE IS
C         LOST AND IERR=4. IN ORDER TO USE THE INT FUNCTION, ARGUMENTS
C         MUST BE FURTHER RESTRICTED NOT TO EXCEED THE LARGEST MACHINE
C         INTEGER, U3=I1MACH(9). THUS, THE MAGNITUDE OF Z AND FNU+N-1 IS
C         RESTRICTED BY MIN(U2,U3). ON 32 BIT MACHINES, U1,U2, AND U3
C         ARE APPROXIMATELY 2.0E+3, 4.2E+6, 2.1E+9 IN SINGLE PRECISION
C         ARITHMETIC AND 1.3E+8, 1.8E+16, 2.1E+9 IN DOUBLE PRECISION
C         ARITHMETIC RESPECTIVELY. THIS MAKES U2 AND U3 LIMITING IN
C         THEIR RESPECTIVE ARITHMETICS. THIS MEANS THAT ONE CAN EXPECT
C         TO RETAIN, IN THE WORST CASES ON 32 BIT MACHINES, NO DIGITS
C         IN SINGLE AND ONLY 7 DIGITS IN DOUBLE PRECISION ARITHMETIC.
C         SIMILAR CONSIDERATIONS HOLD FOR OTHER MACHINES.
C
C         THE APPROXIMATE RELATIVE ERROR IN THE MAGNITUDE OF A COMPLEX
C         BESSEL FUNCTION CAN BE EXPRESSED BY P*10**S WHERE P=MAX(UNIT
C         ROUNDOFF,1.0E-18) IS THE NOMINAL PRECISION AND 10**S REPRE-
C         SENTS THE INCREASE IN ERROR DUE TO ARGUMENT REDUCTION IN THE
C         ELEMENTARY FUNCTIONS. HERE, S=MAX(1,ABS(LOG10(CABS(Z))),
C         ABS(LOG10(FNU))) APPROXIMATELY (I.E. S=MAX(1,ABS(EXPONENT OF
C         CABS(Z),ABS(EXPONENT OF FNU)) ). HOWEVER, THE PHASE ANGLE MAY
C         HAVE ONLY ABSOLUTE ACCURACY. THIS IS MOST LIKELY TO OCCUR WHEN
C         ONE COMPONENT (IN ABSOLUTE VALUE) IS LARGER THAN THE OTHER BY
C         SEVERAL ORDERS OF MAGNITUDE. IF ONE COMPONENT IS 10**K LARGER
C         THAN THE OTHER, THEN ONE CAN EXPECT ONLY MAX(ABS(LOG10(P))-K,
C         0) SIGNIFICANT DIGITS; OR, STATED ANOTHER WAY, WHEN K EXCEEDS
C         THE EXPONENT OF P, NO SIGNIFICANT DIGITS REMAIN IN THE SMALLER
C         COMPONENT. HOWEVER, THE PHASE ANGLE RETAINS ABSOLUTE ACCURACY
C         BECAUSE, IN COMPLEX ARITHMETIC WITH PRECISION P, THE SMALLER
C         COMPONENT WILL NOT (AS A RULE) DECREASE BELOW P TIMES THE
C         MAGNITUDE OF THE LARGER COMPONENT. IN THESE EXTREME CASES,
C         THE PRINCIPAL PHASE ANGLE IS ON THE ORDER OF +P, -P, PI/2-P,
C         OR -PI/2+P.
C
C***REFERENCES  HANDBOOK OF MATHEMATICAL FUNCTIONS BY M. ABRAMOWITZ
C                 AND I. A. STEGUN, NBS AMS SERIES 55, U.S. DEPT. OF
C                 COMMERCE, 1955.
C
C               COMPUTATION OF BESSEL FUNCTIONS OF COMPLEX ARGUMENT
C                 BY D. E. AMOS, SAND83-0083, MAY, 1983.
C
C               COMPUTATION OF BESSEL FUNCTIONS OF COMPLEX ARGUMENT
C                 AND LARGE ORDER BY D. E. AMOS, SAND83-0643, MAY, 1983
C
C               A SUBROUTINE PACKAGE FOR BESSEL FUNCTIONS OF A COMPLEX
C                 ARGUMENT AND NONNEGATIVE ORDER BY D. E. AMOS, SAND85-
C                 1018, MAY, 1985
C
C               A PORTABLE PACKAGE FOR BESSEL FUNCTIONS OF A COMPLEX
C                 ARGUMENT AND NONNEGATIVE ORDER BY D. E. AMOS, ACM
C                 TRANS. MATH. SOFTWARE, VOL. 12, NO. 3, SEPTEMBER 1986,
C                 PP 265-273.
C
C***ROUTINES CALLED  CBESI,CBESK,I1MACH,R1MACH
C***END PROLOGUE  CBESY
C
      COMPLEX CWRK, CY, CI, CIP, CSGN, CSPN, EX, Z, ZU, ZV, ZZ, ZN
      REAL ARG, ELIM, EY, FNU, R1, R2, TAY, XX, YY, R1MACH, ASCLE, RTOL,
     * ATOL, TOL, AA, BB, FFNU, HPI, RHPI, R1M5
      INTEGER I, IERR, IFNU, K, KODE, K1, K2, N, NZ, NZ1, NZ2, I1MACH,
     *I4
      DIMENSION CY(N), CWRK(N), CIP(4)
      DATA CIP(1),CIP(2),CIP(3),CIP(4)/
     * (1.0E0,0.0E0) , (0.0E0,1.0E0) , (-1.0E0,0.0E0) , (0.0E0,-1.0E0) /
      DATA HPI / 1.57079632679489662E0 /
C***FIRST EXECUTABLE STATEMENT  CBESY
      XX = REAL(Z)
      YY = AIMAG(Z)
      IERR = 0
      NZ=0
      IF (XX.EQ.0.0E0 .AND. YY.EQ.0.0E0) IERR=1
      IF (FNU.LT.0.0E0) IERR=1
      IF (KODE.LT.1 .OR. KODE.GT.2) IERR=1
      IF (N.LT.1) IERR=1
      IF (IERR.NE.0) RETURN
      CI = CMPLX(0.0E0,1.0E0)
      ZZ=Z
      IF (YY.LT.0.0E0) ZZ=CONJG(Z)
      ZN = -CI*ZZ
      CALL CBESI(ZN, FNU, KODE, N, CY, NZ1, IERR)
      IF (IERR.NE.0.AND.IERR.NE.3) GO TO 90
      CALL CBESK(ZN, FNU, KODE, N, CWRK, NZ2, IERR)
      IF (IERR.NE.0.AND.IERR.NE.3) GO TO 90
      NZ = MIN0(NZ1,NZ2)
      IFNU = INT(FNU)
      FFNU = FNU - FLOAT(IFNU)
      ARG = HPI*FFNU
      CSGN = CMPLX(COS(ARG),SIN(ARG))
      I4 = MOD(IFNU,4) + 1
      CSGN = CSGN*CIP(I4)
      RHPI = 1.0E0/HPI
      CSPN = CONJG(CSGN)*CMPLX(RHPI,0.0E0)
      CSGN = CSGN*CI
      IF (KODE.EQ.2) GO TO 60
      DO 50 I=1,N
        CY(I) = CSGN*CY(I)-CSPN*CWRK(I)
        CSGN =  CI*CSGN
        CSPN = -CI*CSPN
   50 CONTINUE
      IF (YY.LT.0.0E0) THEN
        DO 55 I=1,N
          CY(I)=CONJG(CY(I))
   55   CONTINUE
      ENDIF
      RETURN
   60 CONTINUE
      R1 = COS(XX)
      R2 = SIN(XX)
      EX = CMPLX(R1,R2)
      TOL = AMAX1(R1MACH(4),1.0E-18)
      K1 = I1MACH(12)
      K2 = I1MACH(13)
      K = MIN0(IABS(K1),IABS(K2))
      R1M5 = R1MACH(5)
C-----------------------------------------------------------------------
C     ELIM IS THE APPROXIMATE EXPONENTIAL UNDER- AND OVERFLOW LIMIT
C-----------------------------------------------------------------------
      ELIM = 2.303E0*(FLOAT(K)*R1M5-3.0E0)
      EY = 0.0E0
      TAY = ABS(YY+YY)
      IF (TAY.LT.ELIM) EY = EXP(-TAY)
      CSPN = EX*CMPLX(EY,0.0E0)*CSPN
      NZ = 0
      RTOL = 1.0E0/TOL
      ASCLE = R1MACH(1)*RTOL*1.0E+3
      DO 80 I=1,N
C----------------------------------------------------------------------
C       CY(I) = CSGN*CY(I)-CSPN*CWRK(I): PRODUCTS ARE COMPUTED IN
C       SCALED MODE IF CY(I) OR CWRK(I) ARE CLOSE TO UNDERFLOW TO 
C       PREVENT UNDERFLOW IN AN INTERMEDIATE COMPUTATION.
C----------------------------------------------------------------------
        ZV = CWRK(I)
        AA=REAL(ZV)
        BB=AIMAG(ZV)
        ATOL=1.0E0
        IF (AMAX1(ABS(AA),ABS(BB)).GT.ASCLE) GO TO 75
          ZV = ZV*CMPLX(RTOL,0.0E0)
          ATOL = TOL
   75   CONTINUE
        ZV = ZV*CSPN
        ZV = ZV*CMPLX(ATOL,0.0E0)
        ZU = CY(I)
        AA=REAL(ZU)
        BB=AIMAG(ZU)
        ATOL=1.0E0
        IF (AMAX1(ABS(AA),ABS(BB)).GT.ASCLE) GO TO 85
          ZU = ZU*CMPLX(RTOL,0.0E0)
          ATOL = TOL
   85   CONTINUE
        ZU = ZU*CSGN
        ZU = ZU*CMPLX(ATOL,0.0E0)
        CY(I) = ZU - ZV
        IF (YY.LT.0.0E0) CY(I)=CONJG(CY(I))
        IF (CY(I).EQ.CMPLX(0.0E0,0.0E0) .AND. EY.EQ.0.0E0) NZ = NZ + 1
        CSGN =  CI*CSGN
        CSPN = -CI*CSPN
   80 CONTINUE
      RETURN
   90 CONTINUE
      NZ = 0
      RETURN
      END
      SUBROUTINE CAIRY(Z, ID, KODE, AI, NZ, IERR)
C***BEGIN PROLOGUE  CAIRY
C***DATE WRITTEN   830501   (YYMMDD)
C***REVISION DATE  890801, 930101   (YYMMDD)
C***CATEGORY NO.  B5K
C***KEYWORDS  AIRY FUNCTION,BESSEL FUNCTIONS OF ORDER ONE THIRD
C***AUTHOR  AMOS, DONALD E., SANDIA NATIONAL LABORATORIES
C***PURPOSE  TO COMPUTE AIRY FUNCTIONS AI(Z) AND DAI(Z) FOR COMPLEX Z
C***DESCRIPTION
C
C         ON KODE=1, CAIRY COMPUTES THE COMPLEX AIRY FUNCTION AI(Z) OR
C         ITS DERIVATIVE DAI(Z)/DZ ON ID=0 OR ID=1 RESPECTIVELY. ON
C         KODE=2, A SCALING OPTION CEXP(ZTA)*AI(Z) OR CEXP(ZTA)*
C         DAI(Z)/DZ IS PROVIDED TO REMOVE THE EXPONENTIAL DECAY IN
C         -PI/3.LT.ARG(Z).LT.PI/3 AND THE EXPONENTIAL GROWTH IN
C         PI/3.LT.ABS(ARG(Z)).LT.PI WHERE ZTA=(2/3)*Z*CSQRT(Z)
C
C         WHILE THE AIRY FUNCTIONS AI(Z) AND DAI(Z)/DZ ARE ANALYTIC IN
C         THE WHOLE Z PLANE, THE CORRESPONDING SCALED FUNCTIONS DEFINED
C         FOR KODE=2 HAVE A CUT ALONG THE NEGATIVE REAL AXIS.
C         DEFINITIONS AND NOTATION ARE FOUND IN THE NBS HANDBOOK OF
C         MATHEMATICAL FUNCTIONS (REF. 1).
C
C         INPUT
C           Z      - Z=CMPLX(X,Y)
C           ID     - ORDER OF DERIVATIVE, ID=0 OR ID=1
C           KODE   - A PARAMETER TO INDICATE THE SCALING OPTION
C                    KODE= 1  RETURNS
C                             AI=AI(Z)                ON ID=0 OR
C                             AI=DAI(Z)/DZ            ON ID=1
C                        = 2  RETURNS
C                             AI=CEXP(ZTA)*AI(Z)       ON ID=0 OR
C                             AI=CEXP(ZTA)*DAI(Z)/DZ   ON ID=1 WHERE
C                             ZTA=(2/3)*Z*CSQRT(Z)
C
C         OUTPUT
C           AI     - COMPLEX ANSWER DEPENDING ON THE CHOICES FOR ID AND
C                    KODE
C           NZ     - UNDERFLOW INDICATOR
C                    NZ= 0   , NORMAL RETURN
C                    NZ= 1   , AI=CMPLX(0.0,0.0) DUE TO UNDERFLOW IN
C                              -PI/3.LT.ARG(Z).LT.PI/3 ON KODE=1
C           IERR   - ERROR FLAG
C                    IERR=0, NORMAL RETURN - COMPUTATION COMPLETED
C                    IERR=1, INPUT ERROR   - NO COMPUTATION
C                    IERR=2, OVERFLOW      - NO COMPUTATION, REAL(ZTA)
C                            TOO LARGE WITH KODE=1.
C                    IERR=3, CABS(Z) LARGE      - COMPUTATION COMPLETED
C                            LOSSES OF SIGNIFCANCE BY ARGUMENT REDUCTION
C                            PRODUCE LESS THAN HALF OF MACHINE ACCURACY
C                    IERR=4, CABS(Z) TOO LARGE  - NO COMPUTATION
C                            COMPLETE LOSS OF ACCURACY BY ARGUMENT
C                            REDUCTION
C                    IERR=5, ERROR              - NO COMPUTATION,
C                            ALGORITHM TERMINATION CONDITION NOT MET
C
C
C***LONG DESCRIPTION
C
C         AI AND DAI ARE COMPUTED FOR CABS(Z).GT.1.0 FROM THE K BESSEL
C         FUNCTIONS BY
C
C            AI(Z)=C*SQRT(Z)*K(1/3,ZTA) , DAI(Z)=-C*Z*K(2/3,ZTA)
C                           C=1.0/(PI*SQRT(3.0))
C                           ZTA=(2/3)*Z**(3/2)
C
C         WITH THE POWER SERIES FOR CABS(Z).LE.1.0.
C
C         IN MOST COMPLEX VARIABLE COMPUTATION, ONE MUST EVALUATE ELE-
C         MENTARY FUNCTIONS. WHEN THE MAGNITUDE OF Z IS LARGE, LOSSES
C         OF SIGNIFICANCE BY ARGUMENT REDUCTION OCCUR. CONSEQUENTLY, IF
C         THE MAGNITUDE OF ZETA=(2/3)*Z**1.5 EXCEEDS U1=SQRT(0.5/UR),
C         THEN LOSSES EXCEEDING HALF PRECISION ARE LIKELY AND AN ERROR
C         FLAG IERR=3 IS TRIGGERED WHERE UR=R1MACH(4)=UNIT ROUNDOFF.
C         ALSO, IF THE MAGNITUDE OF ZETA IS LARGER THAN U2=0.5/UR, THEN
C         ALL SIGNIFICANCE IS LOST AND IERR=4. IN ORDER TO USE THE INT
C         FUNCTION, ZETA MUST BE FURTHER RESTRICTED NOT TO EXCEED THE
C         LARGEST INTEGER, U3=I1MACH(9). THUS, THE MAGNITUDE OF ZETA
C         MUST BE RESTRICTED BY MIN(U2,U3). ON 32 BIT MACHINES, U1,U2,
C         AND U3 ARE APPROXIMATELY 2.0E+3, 4.2E+6, 2.1E+9 IN SINGLE
C         PRECISION ARITHMETIC AND 1.3E+8, 1.8E+16, 2.1E+9 IN DOUBLE
C         PRECISION ARITHMETIC RESPECTIVELY. THIS MAKES U2 AND U3 LIMIT-
C         ING IN THEIR RESPECTIVE ARITHMETICS. THIS MEANS THAT THE MAG-
C         NITUDE OF Z CANNOT EXCEED 3.1E+4 IN SINGLE AND 2.1E+6 IN
C         DOUBLE PRECISION ARITHMETIC. THIS ALSO MEANS THAT ONE CAN
C         EXPECT TO RETAIN, IN THE WORST CASES ON 32 BIT MACHINES,
C         NO DIGITS IN SINGLE PRECISION AND ONLY 7 DIGITS IN DOUBLE
C         PRECISION ARITHMETIC. SIMILAR CONSIDERATIONS HOLD FOR OTHER
C         MACHINES.
C
C         THE APPROXIMATE RELATIVE ERROR IN THE MAGNITUDE OF A COMPLEX
C         BESSEL FUNCTION CAN BE EXPRESSED BY P*10**S WHERE P=MAX(UNIT
C         ROUNDOFF,1.0E-18) IS THE NOMINAL PRECISION AND 10**S REPRE-
C         SENTS THE INCREASE IN ERROR DUE TO ARGUMENT REDUCTION IN THE
C         ELEMENTARY FUNCTIONS. HERE, S=MAX(1,ABS(LOG10(CABS(Z))),
C         ABS(LOG10(FNU))) APPROXIMATELY (I.E. S=MAX(1,ABS(EXPONENT OF
C         CABS(Z),ABS(EXPONENT OF FNU)) ). HOWEVER, THE PHASE ANGLE MAY
C         HAVE ONLY ABSOLUTE ACCURACY. THIS IS MOST LIKELY TO OCCUR WHEN
C         ONE COMPONENT (IN ABSOLUTE VALUE) IS LARGER THAN THE OTHER BY
C         SEVERAL ORDERS OF MAGNITUDE. IF ONE COMPONENT IS 10**K LARGER
C         THAN THE OTHER, THEN ONE CAN EXPECT ONLY MAX(ABS(LOG10(P))-K,
C         0) SIGNIFICANT DIGITS; OR, STATED ANOTHER WAY, WHEN K EXCEEDS
C         THE EXPONENT OF P, NO SIGNIFICANT DIGITS REMAIN IN THE SMALLER
C         COMPONENT. HOWEVER, THE PHASE ANGLE RETAINS ABSOLUTE ACCURACY
C         BECAUSE, IN COMPLEX ARITHMETIC WITH PRECISION P, THE SMALLER
C         COMPONENT WILL NOT (AS A RULE) DECREASE BELOW P TIMES THE
C         MAGNITUDE OF THE LARGER COMPONENT. IN THESE EXTREME CASES,
C         THE PRINCIPAL PHASE ANGLE IS ON THE ORDER OF +P, -P, PI/2-P,
C         OR -PI/2+P.
C
C***REFERENCES  HANDBOOK OF MATHEMATICAL FUNCTIONS BY M. ABRAMOWITZ
C                 AND I. A. STEGUN, NBS AMS SERIES 55, U.S. DEPT. OF
C                 COMMERCE, 1955.
C
C               COMPUTATION OF BESSEL FUNCTIONS OF COMPLEX ARGUMENT
C                 AND LARGE ORDER BY D. E. AMOS, SAND83-0643, MAY, 1983
C
C               A SUBROUTINE PACKAGE FOR BESSEL FUNCTIONS OF A COMPLEX
C                 ARGUMENT AND NONNEGATIVE ORDER BY D. E. AMOS, SAND85-
C                 1018, MAY, 1985
C
C               A PORTABLE PACKAGE FOR BESSEL FUNCTIONS OF A COMPLEX
C                 ARGUMENT AND NONNEGATIVE ORDER BY D. E. AMOS, ACM
C                 TRANS. MATH. SOFTWARE, VOL. 12, NO. 3, SEPTEMBER 1986,
C                 PP 265-273.
C
C***ROUTINES CALLED  CACAI,CBKNU,I1MACH,R1MACH
C***END PROLOGUE  CAIRY
      COMPLEX AI, CONE, CSQ, CY, S1, S2, TRM1, TRM2, Z, ZTA, Z3
      REAL AA, AD, AK, ALIM, ATRM, AZ, AZ3, BK, CK, COEF, C1, C2, DIG,
     * DK, D1, D2, ELIM, FID, FNU, RL, R1M5, SFAC, TOL, TTH, ZI, ZR,
     * Z3I, Z3R, R1MACH, BB, ALAZ
      INTEGER ID, IERR, IFLAG, K, KODE, K1, K2, MR, NN, NZ, I1MACH
      DIMENSION CY(1)
      DATA TTH, C1, C2, COEF /6.66666666666666667E-01,
     * 3.55028053887817240E-01,2.58819403792806799E-01,
     * 1.83776298473930683E-01/
      DATA  CONE / (1.0E0,0.0E0) /
C***FIRST EXECUTABLE STATEMENT  CAIRY
      IERR = 0
      NZ=0
      IF (ID.LT.0 .OR. ID.GT.1) IERR=1
      IF (KODE.LT.1 .OR. KODE.GT.2) IERR=1
      IF (IERR.NE.0) RETURN
      AZ = CABS(Z)
      TOL = AMAX1(R1MACH(4),1.0E-18)
      FID = FLOAT(ID)
      IF (AZ.GT.1.0E0) GO TO 60
C-----------------------------------------------------------------------
C     POWER SERIES FOR CABS(Z).LE.1.
C-----------------------------------------------------------------------
      S1 = CONE
      S2 = CONE
      IF (AZ.LT.TOL) GO TO 160
      AA = AZ*AZ
      IF (AA.LT.TOL/AZ) GO TO 40
      TRM1 = CONE
      TRM2 = CONE
      ATRM = 1.0E0
      Z3 = Z*Z*Z
      AZ3 = AZ*AA
      AK = 2.0E0 + FID
      BK = 3.0E0 - FID - FID
      CK = 4.0E0 - FID
      DK = 3.0E0 + FID + FID
      D1 = AK*DK
      D2 = BK*CK
      AD = AMIN1(D1,D2)
      AK = 24.0E0 + 9.0E0*FID
      BK = 30.0E0 - 9.0E0*FID
      Z3R = REAL(Z3)
      Z3I = AIMAG(Z3)
      DO 30 K=1,25
        TRM1 = TRM1*CMPLX(Z3R/D1,Z3I/D1)
        S1 = S1 + TRM1
        TRM2 = TRM2*CMPLX(Z3R/D2,Z3I/D2)
        S2 = S2 + TRM2
        ATRM = ATRM*AZ3/AD
        D1 = D1 + AK
        D2 = D2 + BK
        AD = AMIN1(D1,D2)
        IF (ATRM.LT.TOL*AD) GO TO 40
        AK = AK + 18.0E0
        BK = BK + 18.0E0
   30 CONTINUE
   40 CONTINUE
      IF (ID.EQ.1) GO TO 50
      AI = S1*CMPLX(C1,0.0E0) - Z*S2*CMPLX(C2,0.0E0)
      IF (KODE.EQ.1) RETURN
      ZTA = Z*CSQRT(Z)*CMPLX(TTH,0.0E0)
      AI = AI*CEXP(ZTA)
      RETURN
   50 CONTINUE
      AI = -S2*CMPLX(C2,0.0E0)
      IF (AZ.GT.TOL) AI = AI + Z*Z*S1*CMPLX(C1/(1.0E0+FID),0.0E0)
      IF (KODE.EQ.1) RETURN
      ZTA = Z*CSQRT(Z)*CMPLX(TTH,0.0E0)
      AI = AI*CEXP(ZTA)
      RETURN
C-----------------------------------------------------------------------
C     CASE FOR CABS(Z).GT.1.0
C-----------------------------------------------------------------------
   60 CONTINUE
      FNU = (1.0E0+FID)/3.0E0
C-----------------------------------------------------------------------
C     SET PARAMETERS RELATED TO MACHINE CONSTANTS.
C     TOL IS THE APPROXIMATE UNIT ROUNDOFF LIMITED TO 1.0E-18.
C     ELIM IS THE APPROXIMATE EXPONENTIAL OVER- AND UNDERFLOW LIMIT.
C     EXP(-ELIM).LT.EXP(-ALIM)=EXP(-ELIM)/TOL    AND
C     EXP(ELIM).GT.EXP(ALIM)=EXP(ELIM)*TOL       ARE INTERVALS NEAR
C     UNDERFLOW AND OVERFLOW LIMITS WHERE SCALED ARITHMETIC IS DONE.
C     RL IS THE LOWER BOUNDARY OF THE ASYMPTOTIC EXPANSION FOR LARGE Z.
C     DIG = NUMBER OF BASE 10 DIGITS IN TOL = 10**(-DIG).
C-----------------------------------------------------------------------
      K1 = I1MACH(12)
      K2 = I1MACH(13)
      R1M5 = R1MACH(5)
      K = MIN0(IABS(K1),IABS(K2))
      ELIM = 2.303E0*(FLOAT(K)*R1M5-3.0E0)
      K1 = I1MACH(11) - 1
      AA = R1M5*FLOAT(K1)
      DIG = AMIN1(AA,18.0E0)
      AA = AA*2.303E0
      ALIM = ELIM + AMAX1(-AA,-41.45E0)
      RL = 1.2E0*DIG + 3.0E0
      ALAZ=ALOG(AZ)
C-----------------------------------------------------------------------
C     TEST FOR RANGE
C-----------------------------------------------------------------------
      AA=0.5E0/TOL
      BB=FLOAT(I1MACH(9))*0.5E0
      AA=AMIN1(AA,BB)
      AA=AA**TTH
      IF (AZ.GT.AA) GO TO 260
      AA=SQRT(AA)
      IF (AZ.GT.AA) IERR=3
      CSQ=CSQRT(Z)
      ZTA=Z*CSQ*CMPLX(TTH,0.0E0)
C-----------------------------------------------------------------------
C     RE(ZTA).LE.0 WHEN RE(Z).LT.0, ESPECIALLY WHEN IM(Z) IS SMALL
C-----------------------------------------------------------------------
      IFLAG = 0
      SFAC = 1.0E0
      ZI = AIMAG(Z)
      ZR = REAL(Z)
      AK = AIMAG(ZTA)
      IF (ZR.GE.0.0E0) GO TO 70
      BK = REAL(ZTA)
      CK = -ABS(BK)
      ZTA = CMPLX(CK,AK)
   70 CONTINUE
      IF (ZI.NE.0.0E0) GO TO 80
      IF (ZR.GT.0.0E0) GO TO 80
      ZTA = CMPLX(0.0E0,AK)
   80 CONTINUE
      AA = REAL(ZTA)
      IF (AA.GE.0.0E0 .AND. ZR.GT.0.0E0) GO TO 100
      IF (KODE.EQ.2) GO TO 90
C-----------------------------------------------------------------------
C     OVERFLOW TEST
C-----------------------------------------------------------------------
      IF (AA.GT.(-ALIM)) GO TO 90
      AA = -AA + 0.25E0*ALAZ
      IFLAG = 1
      SFAC = TOL
      IF (AA.GT.ELIM) GO TO 240
   90 CONTINUE
C-----------------------------------------------------------------------
C     CBKNU AND CACAI RETURN EXP(ZTA)*K(FNU,ZTA) ON KODE=2
C-----------------------------------------------------------------------
      MR = 1
      IF (ZI.LT.0.0E0) MR = -1
      CALL CACAI(ZTA, FNU, KODE, MR, 1, CY, NN, RL, TOL, ELIM, ALIM)
      IF (NN.LT.0) GO TO 250
      NZ = NZ + NN
      GO TO 120
  100 CONTINUE
      IF (KODE.EQ.2) GO TO 110
C-----------------------------------------------------------------------
C     UNDERFLOW TEST
C-----------------------------------------------------------------------
      IF (AA.LT.ALIM) GO TO 110
      AA = -AA - 0.25E0*ALAZ
      IFLAG = 2
      SFAC = 1.0E0/TOL
      IF (AA.LT.(-ELIM)) GO TO 180
  110 CONTINUE
      CALL CBKNU(ZTA, FNU, KODE, 1, CY, NZ, TOL, ELIM, ALIM)
  120 CONTINUE
      S1 = CY(1)*CMPLX(COEF,0.0E0)
      IF (IFLAG.NE.0) GO TO 140
      IF (ID.EQ.1) GO TO 130
      AI = CSQ*S1
      RETURN
  130 AI = -Z*S1
      RETURN
  140 CONTINUE
      S1 = S1*CMPLX(SFAC,0.0E0)
      IF (ID.EQ.1) GO TO 150
      S1 = S1*CSQ
      AI = S1*CMPLX(1.0E0/SFAC,0.0E0)
      RETURN
  150 CONTINUE
      S1 = -S1*Z
      AI = S1*CMPLX(1.0E0/SFAC,0.0E0)
      RETURN
  160 CONTINUE
      AA = 1.0E+3*R1MACH(1)
      S1 = CMPLX(0.0E0,0.0E0)
      IF (ID.EQ.1) GO TO 170
      IF (AZ.GT.AA) S1 = CMPLX(C2,0.0E0)*Z
      AI = CMPLX(C1,0.0E0) - S1
      RETURN
  170 CONTINUE
      AI = -CMPLX(C2,0.0E0)
      AA = SQRT(AA)
      IF (AZ.GT.AA) S1 = Z*Z*CMPLX(0.5E0,0.0E0)
      AI = AI + S1*CMPLX(C1,0.0E0)
      RETURN
  180 CONTINUE
      NZ = 1
      AI = CMPLX(0.0E0,0.0E0)
      RETURN
  240 CONTINUE
      NZ = 0
      IERR=2
      RETURN
  250 CONTINUE
      IF(NN.EQ.(-1)) GO TO 240
      NZ=0
      IERR=5
      RETURN
  260 CONTINUE
      IERR=4
      NZ=0
      RETURN
      END
      SUBROUTINE CBIRY(Z, ID, KODE, BI, IERR)
C***BEGIN PROLOGUE  CBIRY
C***DATE WRITTEN   830501   (YYMMDD)
C***REVISION DATE  890801, 930101   (YYMMDD)
C***CATEGORY NO.  B5K
C***KEYWORDS  AIRY FUNCTION,BESSEL FUNCTIONS OF ORDER ONE THIRD
C***AUTHOR  AMOS, DONALD E., SANDIA NATIONAL LABORATORIES
C***PURPOSE  TO COMPUTE AIRY FUNCTIONS BI(Z) AND DBI(Z) FOR COMPLEX Z
C***DESCRIPTION
C
C         ON KODE=1, CBIRY COMPUTES THE COMPLEX AIRY FUNCTION BI(Z) OR
C         ITS DERIVATIVE DBI(Z)/DZ ON ID=0 OR ID=1 RESPECTIVELY. ON
C         KODE=2, A SCALING OPTION CEXP(-AXZTA)*BI(Z) OR CEXP(-AXZTA)*
C         DBI(Z)/DZ IS PROVIDED TO REMOVE THE EXPONENTIAL BEHAVIOR IN
C         BOTH THE LEFT AND RIGHT HALF PLANES WHERE
C         ZTA=(2/3)*Z*CSQRT(Z)=CMPLX(XZTA,YZTA) AND AXZTA=ABS(XZTA).
C         DEFINITIONS AND NOTATION ARE FOUND IN THE NBS HANDBOOK OF
C         MATHEMATICAL FUNCTIONS (REF. 1).
C
C         INPUT
C           Z      - Z=CMPLX(X,Y)
C           ID     - ORDER OF DERIVATIVE, ID=0 OR ID=1
C           KODE   - A PARAMETER TO INDICATE THE SCALING OPTION
C                    KODE= 1  RETURNS
C                             BI=BI(Z)                 ON ID=0 OR
C                             BI=DBI(Z)/DZ             ON ID=1
C                        = 2  RETURNS
C                             BI=CEXP(-AXZTA)*BI(Z)     ON ID=0 OR
C                             BI=CEXP(-AXZTA)*DBI(Z)/DZ ON ID=1 WHERE
C                             ZTA=(2/3)*Z*CSQRT(Z)=CMPLX(XZTA,YZTA)
C                             AND AXZTA=ABS(XZTA)
C
C         OUTPUT
C           BI     - COMPLEX ANSWER DEPENDING ON THE CHOICES FOR ID AND
C                    KODE
C           IERR   - ERROR FLAG
C                    IERR=0, NORMAL RETURN - COMPUTATION COMPLETED
C                    IERR=1, INPUT ERROR   - NO COMPUTATION
C                    IERR=2, OVERFLOW      - NO COMPUTATION, REAL(Z)
C                            TOO LARGE WITH KODE=1
C                    IERR=3, CABS(Z) LARGE      - COMPUTATION COMPLETED
C                            LOSSES OF SIGNIFCANCE BY ARGUMENT REDUCTION
C                            PRODUCE LESS THAN HALF OF MACHINE ACCURACY
C                    IERR=4, CABS(Z) TOO LARGE  - NO COMPUTATION
C                            COMPLETE LOSS OF ACCURACY BY ARGUMENT
C                            REDUCTION
C                    IERR=5, ERROR              - NO COMPUTATION,
C                            ALGORITHM TERMINATION CONDITION NOT MET
C
C***LONG DESCRIPTION
C
C         BI AND DBI ARE COMPUTED FOR CABS(Z).GT.1.0 FROM THE I BESSEL
C         FUNCTIONS BY
C
C                BI(Z)=C*SQRT(Z)*( I(-1/3,ZTA) + I(1/3,ZTA) )
C               DBI(Z)=C *  Z  * ( I(-2/3,ZTA) + I(2/3,ZTA) )
C                               C=1.0/SQRT(3.0)
C                               ZTA=(2/3)*Z**(3/2)
C
C         WITH THE POWER SERIES FOR CABS(Z).LE.1.0.
C
C         IN MOST COMPLEX VARIABLE COMPUTATION, ONE MUST EVALUATE ELE-
C         MENTARY FUNCTIONS. WHEN THE MAGNITUDE OF Z IS LARGE, LOSSES
C         OF SIGNIFICANCE BY ARGUMENT REDUCTION OCCUR. CONSEQUENTLY, IF
C         THE MAGNITUDE OF ZETA=(2/3)*Z**1.5 EXCEEDS U1=SQRT(0.5/UR),
C         THEN LOSSES EXCEEDING HALF PRECISION ARE LIKELY AND AN ERROR
C         FLAG IERR=3 IS TRIGGERED WHERE UR=R1MACH(4)=UNIT ROUNDOFF.
C         ALSO, IF THE MAGNITUDE OF ZETA IS LARGER THAN U2=0.5/UR, THEN
C         ALL SIGNIFICANCE IS LOST AND IERR=4. IN ORDER TO USE THE INT
C         FUNCTION, ZETA MUST BE FURTHER RESTRICTED NOT TO EXCEED THE
C         LARGEST INTEGER, U3=I1MACH(9). THUS, THE MAGNITUDE OF ZETA
C         MUST BE RESTRICTED BY MIN(U2,U3). ON 32 BIT MACHINES, U1,U2,
C         AND U3 ARE APPROXIMATELY 2.0E+3, 4.2E+6, 2.1E+9 IN SINGLE
C         PRECISION ARITHMETIC AND 1.3E+8, 1.8E+16, 2.1E+9 IN DOUBLE
C         PRECISION ARITHMETIC RESPECTIVELY. THIS MAKES U2 AND U3 LIMIT-
C         ING IN THEIR RESPECTIVE ARITHMETICS. THIS MEANS THAT THE MAG-
C         NITUDE OF Z CANNOT EXCEED 3.1E+4 IN SINGLE AND 2.1E+6 IN
C         DOUBLE PRECISION ARITHMETIC. THIS ALSO MEANS THAT ONE CAN
C         EXPECT TO RETAIN, IN THE WORST CASES ON 32 BIT MACHINES,
C         NO DIGITS IN SINGLE PRECISION AND ONLY 7 DIGITS IN DOUBLE
C         PRECISION ARITHMETIC.
C
C         THE APPROXIMATE RELATIVE ERROR IN THE MAGNITUDE OF A COMPLEX
C         BESSEL FUNCTION CAN BE EXPRESSED BY P*10**S WHERE P=MAX(UNIT
C         ROUNDOFF,1.0E-18) IS THE NOMINAL PRECISION AND 10**S REPRE-
C         SENTS THE INCREASE IN ERROR DUE TO ARGUMENT REDUCTION IN THE
C         ELEMENTARY FUNCTIONS. HERE, S=MAX(1,ABS(LOG10(CABS(Z))),
C         ABS(LOG10(FNU))) APPROXIMATELY (I.E. S=MAX(1,ABS(EXPONENT OF
C         CABS(Z),ABS(EXPONENT OF FNU)) ). HOWEVER, THE PHASE ANGLE MAY
C         HAVE ONLY ABSOLUTE ACCURACY. THIS IS MOST LIKELY TO OCCUR WHEN
C         ONE COMPONENT (IN ABSOLUTE VALUE) IS LARGER THAN THE OTHER BY
C         SEVERAL ORDERS OF MAGNITUDE. IF ONE COMPONENT IS 10**K LARGER
C         THAN THE OTHER, THEN ONE CAN EXPECT ONLY MAX(ABS(LOG10(P))-K,
C         0) SIGNIFICANT DIGITS; OR, STATED ANOTHER WAY, WHEN K EXCEEDS
C         THE EXPONENT OF P, NO SIGNIFICANT DIGITS REMAIN IN THE SMALLER
C         COMPONENT. HOWEVER, THE PHASE ANGLE RETAINS ABSOLUTE ACCURACY
C         BECAUSE, IN COMPLEX ARITHMETIC WITH PRECISION P, THE SMALLER
C         COMPONENT WILL NOT (AS A RULE) DECREASE BELOW P TIMES THE
C         MAGNITUDE OF THE LARGER COMPONENT. IN THESE EXTREME CASES,
C         THE PRINCIPAL PHASE ANGLE IS ON THE ORDER OF +P, -P, PI/2-P,
C         OR -PI/2+P.
C
C***REFERENCES  HANDBOOK OF MATHEMATICAL FUNCTIONS BY M. ABRAMOWITZ
C                 AND I. A. STEGUN, NBS AMS SERIES 55, U.S. DEPT. OF
C                 COMMERCE, 1955.
C
C               COMPUTATION OF BESSEL FUNCTIONS OF COMPLEX ARGUMENT
C                 AND LARGE ORDER BY D. E. AMOS, SAND83-0643, MAY, 1983
C
C               A SUBROUTINE PACKAGE FOR BESSEL FUNCTIONS OF A COMPLEX
C                 ARGUMENT AND NONNEGATIVE ORDER BY D. E. AMOS, SAND85-
C                 1018, MAY, 1985
C
C               A PORTABLE PACKAGE FOR BESSEL FUNCTIONS OF A COMPLEX
C                 ARGUMENT AND NONNEGATIVE ORDER BY D. E. AMOS, ACM
C                 TRANS. MATH. SOFTWARE, VOL. 12, NO. 3, SEPTEMBER 1986,
C                 PP 265-273.
C
C***ROUTINES CALLED  CBINU,I1MACH,R1MACH
C***END PROLOGUE  CBIRY
      COMPLEX BI, CONE, CSQ, CY, S1, S2, TRM1, TRM2, Z, ZTA, Z3
      REAL AA, AD, AK, ALIM, ATRM, AZ, AZ3, BB, BK, CK, COEF, C1, C2,
     * DIG, DK, D1, D2, ELIM, FID, FMR, FNU, FNUL, PI, RL, R1M5, SFAC,
     * TOL, TTH, ZI, ZR, Z3I, Z3R, R1MACH
      INTEGER ID, IERR, K, KODE, K1, K2, NZ, I1MACH
      DIMENSION CY(2)
      DATA TTH, C1, C2, COEF, PI /6.66666666666666667E-01,
     * 6.14926627446000736E-01,4.48288357353826359E-01,
     * 5.77350269189625765E-01,3.14159265358979324E+00/
      DATA  CONE / (1.0E0,0.0E0) /
C***FIRST EXECUTABLE STATEMENT  CBIRY
      IERR = 0
      NZ=0
      IF (ID.LT.0 .OR. ID.GT.1) IERR=1
      IF (KODE.LT.1 .OR. KODE.GT.2) IERR=1
      IF (IERR.NE.0) RETURN
      AZ = CABS(Z)
      TOL = AMAX1(R1MACH(4),1.0E-18)
      FID = FLOAT(ID)
      IF (AZ.GT.1.0E0) GO TO 60
C-----------------------------------------------------------------------
C     POWER SERIES FOR CABS(Z).LE.1.
C-----------------------------------------------------------------------
      S1 = CONE
      S2 = CONE
      IF (AZ.LT.TOL) GO TO 110
      AA = AZ*AZ
      IF (AA.LT.TOL/AZ) GO TO 40
      TRM1 = CONE
      TRM2 = CONE
      ATRM = 1.0E0
      Z3 = Z*Z*Z
      AZ3 = AZ*AA
      AK = 2.0E0 + FID
      BK = 3.0E0 - FID - FID
      CK = 4.0E0 - FID
      DK = 3.0E0 + FID + FID
      D1 = AK*DK
      D2 = BK*CK
      AD = AMIN1(D1,D2)
      AK = 24.0E0 + 9.0E0*FID
      BK = 30.0E0 - 9.0E0*FID
      Z3R = REAL(Z3)
      Z3I = AIMAG(Z3)
      DO 30 K=1,25
        TRM1 = TRM1*CMPLX(Z3R/D1,Z3I/D1)
        S1 = S1 + TRM1
        TRM2 = TRM2*CMPLX(Z3R/D2,Z3I/D2)
        S2 = S2 + TRM2
        ATRM = ATRM*AZ3/AD
        D1 = D1 + AK
        D2 = D2 + BK
        AD = AMIN1(D1,D2)
        IF (ATRM.LT.TOL*AD) GO TO 40
        AK = AK + 18.0E0
        BK = BK + 18.0E0
   30 CONTINUE
   40 CONTINUE
      IF (ID.EQ.1) GO TO 50
      BI = S1*CMPLX(C1,0.0E0) + Z*S2*CMPLX(C2,0.0E0)
      IF (KODE.EQ.1) RETURN
      ZTA = Z*CSQRT(Z)*CMPLX(TTH,0.0E0)
      AA = REAL(ZTA)
      AA = -ABS(AA)
      BI = BI*CMPLX(EXP(AA),0.0E0)
      RETURN
   50 CONTINUE
      BI = S2*CMPLX(C2,0.0E0)
      IF (AZ.GT.TOL) BI = BI + Z*Z*S1*CMPLX(C1/(1.0E0+FID),0.0E0)
      IF (KODE.EQ.1) RETURN
      ZTA = Z*CSQRT(Z)*CMPLX(TTH,0.0E0)
      AA = REAL(ZTA)
      AA = -ABS(AA)
      BI = BI*CMPLX(EXP(AA),0.0E0)
      RETURN
C-----------------------------------------------------------------------
C     CASE FOR CABS(Z).GT.1.0
C-----------------------------------------------------------------------
   60 CONTINUE
      FNU = (1.0E0+FID)/3.0E0
C-----------------------------------------------------------------------
C     SET PARAMETERS RELATED TO MACHINE CONSTANTS.
C     TOL IS THE APPROXIMATE UNIT ROUNDOFF LIMITED TO 1.0E-18.
C     ELIM IS THE APPROXIMATE EXPONENTIAL OVER- AND UNDERFLOW LIMIT.
C     EXP(-ELIM).LT.EXP(-ALIM)=EXP(-ELIM)/TOL    AND
C     EXP(ELIM).GT.EXP(ALIM)=EXP(ELIM)*TOL       ARE INTERVALS NEAR
C     UNDERFLOW AND OVERFLOW LIMITS WHERE SCALED ARITHMETIC IS DONE.
C     RL IS THE LOWER BOUNDARY OF THE ASYMPTOTIC EXPANSION FOR LARGE Z.
C     DIG = NUMBER OF BASE 10 DIGITS IN TOL = 10**(-DIG).
C     FNUL IS THE LOWER BOUNDARY OF THE ASYMPTOTIC SERIES FOR LARGE FNU.
C-----------------------------------------------------------------------
      K1 = I1MACH(12)
      K2 = I1MACH(13)
      R1M5 = R1MACH(5)
      K = MIN0(IABS(K1),IABS(K2))
      ELIM = 2.303E0*(FLOAT(K)*R1M5-3.0E0)
      K1 = I1MACH(11) - 1
      AA = R1M5*FLOAT(K1)
      DIG = AMIN1(AA,18.0E0)
      AA = AA*2.303E0
      ALIM = ELIM + AMAX1(-AA,-41.45E0)
      RL = 1.2E0*DIG + 3.0E0
      FNUL = 10.0E0 + 6.0E0*(DIG-3.0E0)
C-----------------------------------------------------------------------
C     TEST FOR RANGE
C-----------------------------------------------------------------------
      AA=0.5E0/TOL
      BB=FLOAT(I1MACH(9))*0.5E0
      AA=AMIN1(AA,BB)
      AA=AA**TTH
      IF (AZ.GT.AA) GO TO 190
      AA=SQRT(AA)
      IF (AZ.GT.AA) IERR=3
      CSQ=CSQRT(Z)
      ZTA=Z*CSQ*CMPLX(TTH,0.0E0)
C-----------------------------------------------------------------------
C     RE(ZTA).LE.0 WHEN RE(Z).LT.0, ESPECIALLY WHEN IM(Z) IS SMALL
C-----------------------------------------------------------------------
      SFAC = 1.0E0
      ZI = AIMAG(Z)
      ZR = REAL(Z)
      AK = AIMAG(ZTA)
      IF (ZR.GE.0.0E0) GO TO 70
      BK = REAL(ZTA)
      CK = -ABS(BK)
      ZTA = CMPLX(CK,AK)
   70 CONTINUE
      IF (ZI.EQ.0.0E0 .AND. ZR.LE.0.0E0) ZTA = CMPLX(0.0E0,AK)
      AA = REAL(ZTA)
      IF (KODE.EQ.2) GO TO 80
C-----------------------------------------------------------------------
C     OVERFLOW TEST
C-----------------------------------------------------------------------
      BB = ABS(AA)
      IF (BB.LT.ALIM) GO TO 80
      BB = BB + 0.25E0*ALOG(AZ)
      SFAC = TOL
      IF (BB.GT.ELIM) GO TO 170
   80 CONTINUE
      FMR = 0.0E0
      IF (AA.GE.0.0E0 .AND. ZR.GT.0.0E0) GO TO 90
      FMR = PI
      IF (ZI.LT.0.0E0) FMR = -PI
      ZTA = -ZTA
   90 CONTINUE
C-----------------------------------------------------------------------
C     AA=FACTOR FOR ANALYTIC CONTINUATION OF I(FNU,ZTA)
C     KODE=2 RETURNS EXP(-ABS(XZTA))*I(FNU,ZTA) FROM CBINU
C-----------------------------------------------------------------------
      CALL CBINU(ZTA, FNU, KODE, 1, CY, NZ, RL, FNUL, TOL, ELIM, ALIM)
      IF (NZ.LT.0) GO TO 180
      AA = FMR*FNU
      Z3 = CMPLX(SFAC,0.0E0)
      S1 = CY(1)*CMPLX(COS(AA),SIN(AA))*Z3
      FNU = (2.0E0-FID)/3.0E0
      CALL CBINU(ZTA, FNU, KODE, 2, CY, NZ, RL, FNUL, TOL, ELIM, ALIM)
      CY(1) = CY(1)*Z3
      CY(2) = CY(2)*Z3
C-----------------------------------------------------------------------
C     BACKWARD RECUR ONE STEP FOR ORDERS -1/3 OR -2/3
C-----------------------------------------------------------------------
      S2 = CY(1)*CMPLX(FNU+FNU,0.0E0)/ZTA + CY(2)
      AA = FMR*(FNU-1.0E0)
      S1 = (S1+S2*CMPLX(COS(AA),SIN(AA)))*CMPLX(COEF,0.0E0)
      IF (ID.EQ.1) GO TO 100
      S1 = CSQ*S1
      BI = S1*CMPLX(1.0E0/SFAC,0.0E0)
      RETURN
  100 CONTINUE
      S1 = Z*S1
      BI = S1*CMPLX(1.0E0/SFAC,0.0E0)
      RETURN
  110 CONTINUE
      AA = C1*(1.0E0-FID) + FID*C2
      BI = CMPLX(AA,0.0E0)
      RETURN
  170 CONTINUE
      NZ=0
      IERR=2
      RETURN
  180 CONTINUE
      IF(NZ.EQ.(-1)) GO TO 170
      NZ=0
      IERR=5
      RETURN
  190 CONTINUE
      IERR=4
      NZ=0
      RETURN
      END
      SUBROUTINE CUNIK(ZR, FNU, IKFLG, IPMTR, TOL, INIT, PHI, ZETA1,
     * ZETA2, SUM, CWRK)
C***BEGIN PROLOGUE  CUNIK
C***REFER TO  CBESI,CBESK
C
C        CUNIK COMPUTES PARAMETERS FOR THE UNIFORM ASYMPTOTIC
C        EXPANSIONS OF THE I AND K FUNCTIONS ON IKFLG= 1 OR 2
C        RESPECTIVELY BY
C
C        W(FNU,ZR) = PHI*EXP(ZETA)*SUM
C
C        WHERE       ZETA=-ZETA1 + ZETA2       OR
C                          ZETA1 - ZETA2
C
C        THE FIRST CALL MUST HAVE INIT=0. SUBSEQUENT CALLS WITH THE
C        SAME ZR AND FNU WILL RETURN THE I OR K FUNCTION ON IKFLG=
C        1 OR 2 WITH NO CHANGE IN INIT. CWRK IS A COMPLEX WORK
C        ARRAY. IPMTR=0 COMPUTES ALL PARAMETERS. IPMTR=1 COMPUTES PHI,
C        ZETA1,ZETA2.
C
C***ROUTINES CALLED  (NONE)
C***END PROLOGUE  CUNIK
      COMPLEX CFN, CON, CONE, CRFN, CWRK, CZERO, PHI, S, SR, SUM, T,
     * T2, ZETA1, ZETA2, ZN, ZR
      REAL AC, C, FNU, RFN, TEST, TOL, TSTR, TSTI
      INTEGER I, IKFLG, INIT, IPMTR, J, K, L
      DIMENSION C(120), CWRK(16), CON(2)
      DATA CZERO, CONE / (0.0E0,0.0E0), (1.0E0,0.0E0) /
      DATA CON(1), CON(2)  /
     1(3.98942280401432678E-01,0.0E0),(1.25331413731550025E+00,0.0E0)/
      DATA C(1), C(2), C(3), C(4), C(5), C(6), C(7), C(8), C(9), C(10),
     1     C(11), C(12), C(13), C(14), C(15), C(16), C(17), C(18),
     2     C(19), C(20), C(21), C(22), C(23), C(24)/
     3     1.00000000000000000E+00,    -2.08333333333333333E-01,
     4     1.25000000000000000E-01,     3.34201388888888889E-01,
     5    -4.01041666666666667E-01,     7.03125000000000000E-02,
     6    -1.02581259645061728E+00,     1.84646267361111111E+00,
     7    -8.91210937500000000E-01,     7.32421875000000000E-02,
     8     4.66958442342624743E+00,    -1.12070026162229938E+01,
     9     8.78912353515625000E+00,    -2.36408691406250000E+00,
     A     1.12152099609375000E-01,    -2.82120725582002449E+01,
     B     8.46362176746007346E+01,    -9.18182415432400174E+01,
     C     4.25349987453884549E+01,    -7.36879435947963170E+00,
     D     2.27108001708984375E-01,     2.12570130039217123E+02,
     E    -7.65252468141181642E+02,     1.05999045252799988E+03/
      DATA C(25), C(26), C(27), C(28), C(29), C(30), C(31), C(32),
     1     C(33), C(34), C(35), C(36), C(37), C(38), C(39), C(40),
     2     C(41), C(42), C(43), C(44), C(45), C(46), C(47), C(48)/
     3    -6.99579627376132541E+02,     2.18190511744211590E+02,
     4    -2.64914304869515555E+01,     5.72501420974731445E-01,
     5    -1.91945766231840700E+03,     8.06172218173730938E+03,
     6    -1.35865500064341374E+04,     1.16553933368645332E+04,
     7    -5.30564697861340311E+03,     1.20090291321635246E+03,
     8    -1.08090919788394656E+02,     1.72772750258445740E+00,
     9     2.02042913309661486E+04,    -9.69805983886375135E+04,
     A     1.92547001232531532E+05,    -2.03400177280415534E+05,
     B     1.22200464983017460E+05,    -4.11926549688975513E+04,
     C     7.10951430248936372E+03,    -4.93915304773088012E+02,
     D     6.07404200127348304E+00,    -2.42919187900551333E+05,
     E     1.31176361466297720E+06,    -2.99801591853810675E+06/
      DATA C(49), C(50), C(51), C(52), C(53), C(54), C(55), C(56),
     1     C(57), C(58), C(59), C(60), C(61), C(62), C(63), C(64),
     2     C(65), C(66), C(67), C(68), C(69), C(70), C(71), C(72)/
     3     3.76327129765640400E+06,    -2.81356322658653411E+06,
     4     1.26836527332162478E+06,    -3.31645172484563578E+05,
     5     4.52187689813627263E+04,    -2.49983048181120962E+03,
     6     2.43805296995560639E+01,     3.28446985307203782E+06,
     7    -1.97068191184322269E+07,     5.09526024926646422E+07,
     8    -7.41051482115326577E+07,     6.63445122747290267E+07,
     9    -3.75671766607633513E+07,     1.32887671664218183E+07,
     A    -2.78561812808645469E+06,     3.08186404612662398E+05,
     B    -1.38860897537170405E+04,     1.10017140269246738E+02,
     C    -4.93292536645099620E+07,     3.25573074185765749E+08,
     D    -9.39462359681578403E+08,     1.55359689957058006E+09,
     E    -1.62108055210833708E+09,     1.10684281682301447E+09/
      DATA C(73), C(74), C(75), C(76), C(77), C(78), C(79), C(80),
     1     C(81), C(82), C(83), C(84), C(85), C(86), C(87), C(88),
     2     C(89), C(90), C(91), C(92), C(93), C(94), C(95), C(96)/
     3    -4.95889784275030309E+08,     1.42062907797533095E+08,
     4    -2.44740627257387285E+07,     2.24376817792244943E+06,
     5    -8.40054336030240853E+04,     5.51335896122020586E+02,
     6     8.14789096118312115E+08,    -5.86648149205184723E+09,
     7     1.86882075092958249E+10,    -3.46320433881587779E+10,
     8     4.12801855797539740E+10,    -3.30265997498007231E+10,
     9     1.79542137311556001E+10,    -6.56329379261928433E+09,
     A     1.55927986487925751E+09,    -2.25105661889415278E+08,
     B     1.73951075539781645E+07,    -5.49842327572288687E+05,
     C     3.03809051092238427E+03,    -1.46792612476956167E+10,
     D     1.14498237732025810E+11,    -3.99096175224466498E+11,
     E     8.19218669548577329E+11,    -1.09837515608122331E+12/
      DATA C(97), C(98), C(99), C(100), C(101), C(102), C(103), C(104),
     1     C(105), C(106), C(107), C(108), C(109), C(110), C(111),
     2     C(112), C(113), C(114), C(115), C(116), C(117), C(118)/
     3     1.00815810686538209E+12,    -6.45364869245376503E+11,
     4     2.87900649906150589E+11,    -8.78670721780232657E+10,
     5     1.76347306068349694E+10,    -2.16716498322379509E+09,
     6     1.43157876718888981E+08,    -3.87183344257261262E+06,
     7     1.82577554742931747E+04,     2.86464035717679043E+11,
     8    -2.40629790002850396E+12,     9.10934118523989896E+12,
     9    -2.05168994109344374E+13,     3.05651255199353206E+13,
     A    -3.16670885847851584E+13,     2.33483640445818409E+13,
     B    -1.23204913055982872E+13,     4.61272578084913197E+12,
     C    -1.19655288019618160E+12,     2.05914503232410016E+11,
     D    -2.18229277575292237E+10,     1.24700929351271032E+09/
      DATA C(119), C(120)/
     1    -2.91883881222208134E+07,     1.18838426256783253E+05/
C
      IF (INIT.NE.0) GO TO 40
C-----------------------------------------------------------------------
C     INITIALIZE ALL VARIABLES
C-----------------------------------------------------------------------
      RFN = 1.0E0/FNU
      CRFN = CMPLX(RFN,0.0E0)
C     T = ZR*CRFN
C-----------------------------------------------------------------------
C     OVERFLOW TEST (ZR/FNU TOO SMALL)
C-----------------------------------------------------------------------
      TSTR = REAL(ZR)
      TSTI = AIMAG(ZR)
      TEST = R1MACH(1)*1.0E+3
      AC = FNU*TEST
      IF (ABS(TSTR).GT.AC .OR. ABS(TSTI).GT.AC) GO TO 15
      AC = 2.0E0*ABS(ALOG(TEST))+FNU
      ZETA1 = CMPLX(AC,0.0E0)
      ZETA2 = CMPLX(FNU,0.0E0)
      PHI=CONE
      RETURN
   15 CONTINUE
      T=ZR*CRFN
      S = CONE + T*T
      SR = CSQRT(S)
      CFN = CMPLX(FNU,0.0E0)
      ZN = (CONE+SR)/T
      ZETA1 = CFN*CLOG(ZN)
      ZETA2 = CFN*SR
      T = CONE/SR
      SR = T*CRFN
      CWRK(16) = CSQRT(SR)
      PHI = CWRK(16)*CON(IKFLG)
      IF (IPMTR.NE.0) RETURN
      T2 = CONE/S
      CWRK(1) = CONE
      CRFN = CONE
      AC = 1.0E0
      L = 1
      DO 20 K=2,15
        S = CZERO
        DO 10 J=1,K
          L = L + 1
          S = S*T2 + CMPLX(C(L),0.0E0)
   10   CONTINUE
        CRFN = CRFN*SR
        CWRK(K) = CRFN*S
        AC = AC*RFN
        TSTR = REAL(CWRK(K))
        TSTI = AIMAG(CWRK(K))
        TEST = ABS(TSTR) + ABS(TSTI)
        IF (AC.LT.TOL .AND. TEST.LT.TOL) GO TO 30
   20 CONTINUE
      K = 15
   30 CONTINUE
      INIT = K
   40 CONTINUE
      IF (IKFLG.EQ.2) GO TO 60
C-----------------------------------------------------------------------
C     COMPUTE SUM FOR THE I FUNCTION
C-----------------------------------------------------------------------
      S = CZERO
      DO 50 I=1,INIT
        S = S + CWRK(I)
   50 CONTINUE
      SUM = S
      PHI = CWRK(16)*CON(1)
      RETURN
   60 CONTINUE
C-----------------------------------------------------------------------
C     COMPUTE SUM FOR THE K FUNCTION
C-----------------------------------------------------------------------
      S = CZERO
      T = CONE
      DO 70 I=1,INIT
        S = S + T*CWRK(I)
        T = -T
   70 CONTINUE
      SUM = S
      PHI = CWRK(16)*CON(2)
      RETURN
      END
      SUBROUTINE CUOIK(Z, FNU, KODE, IKFLG, N, Y, NUF, TOL, ELIM, ALIM)
C***BEGIN PROLOGUE  CUOIK
C***REFER TO  CBESI,CBESK,CBESH
C
C     CUOIK COMPUTES THE LEADING TERMS OF THE UNIFORM ASYMPTOTIC
C     EXPANSIONS FOR THE I AND K FUNCTIONS AND COMPARES THEM
C     (IN LOGARITHMIC FORM) TO ALIM AND ELIM FOR OVER AND UNDERFLOW
C     WHERE ALIM.LT.ELIM. IF THE MAGNITUDE, BASED ON THE LEADING
C     EXPONENTIAL, IS LESS THAN ALIM OR GREATER THAN -ALIM, THEN
C     THE RESULT IS ON SCALE. IF NOT, THEN A REFINED TEST USING OTHER
C     MULTIPLIERS (IN LOGARITHMIC FORM) IS MADE BASED ON ELIM. HERE
C     EXP(-ELIM)=SMALLEST MACHINE NUMBER*1.0E+3 AND EXP(-ALIM)=
C     EXP(-ELIM)/TOL
C
C     IKFLG=1 MEANS THE I SEQUENCE IS TESTED
C          =2 MEANS THE K SEQUENCE IS TESTED
C     NUF = 0 MEANS THE LAST MEMBER OF THE SEQUENCE IS ON SCALE
C         =-1 MEANS AN OVERFLOW WOULD OCCUR
C     IKFLG=1 AND NUF.GT.0 MEANS THE LAST NUF Y VALUES WERE SET TO ZERO
C             THE FIRST N-NUF VALUES MUST BE SET BY ANOTHER ROUTINE
C     IKFLG=2 AND NUF.EQ.N MEANS ALL Y VALUES WERE SET TO ZERO
C     IKFLG=2 AND 0.LT.NUF.LT.N NOT CONSIDERED. Y MUST BE SET BY
C             ANOTHER ROUTINE
C
C***ROUTINES CALLED  CUCHK,CUNHJ,CUNIK,R1MACH
C***END PROLOGUE  CUOIK
      COMPLEX ARG, ASUM, BSUM, CWRK, CZ, CZERO, PHI, SUM, Y, Z, ZB,
     * ZETA1, ZETA2, ZN, ZR
      REAL AARG, AIC, ALIM, APHI, ASCLE, AX, AY, ELIM, FNN, FNU, GNN,
     * GNU, RCZ, TOL, X, YY
      INTEGER I, IFORM, IKFLG, INIT, KODE, N, NN, NUF, NW
      DIMENSION Y(N), CWRK(16)
      DATA CZERO / (0.0E0,0.0E0) /
      DATA AIC / 1.265512123484645396E+00 /
      NUF = 0
      NN = N
      X = REAL(Z)
      ZR = Z
      IF (X.LT.0.0E0) ZR = -Z
      ZB = ZR
      YY = AIMAG(ZR)
      AX = ABS(X)*1.7321E0
      AY = ABS(YY)
      IFORM = 1
      IF (AY.GT.AX) IFORM = 2
      GNU = AMAX1(FNU,1.0E0)
      IF (IKFLG.EQ.1) GO TO 10
      FNN = FLOAT(NN)
      GNN = FNU + FNN - 1.0E0
      GNU = AMAX1(GNN,FNN)
   10 CONTINUE
C-----------------------------------------------------------------------
C     ONLY THE MAGNITUDE OF ARG AND PHI ARE NEEDED ALONG WITH THE
C     REAL PARTS OF ZETA1, ZETA2 AND ZB. NO ATTEMPT IS MADE TO GET
C     THE SIGN OF THE IMAGINARY PART CORRECT.
C-----------------------------------------------------------------------
      IF (IFORM.EQ.2) GO TO 20
      INIT = 0
      CALL CUNIK(ZR, GNU, IKFLG, 1, TOL, INIT, PHI, ZETA1, ZETA2, SUM,
     * CWRK)
      CZ = -ZETA1 + ZETA2
      GO TO 40
   20 CONTINUE
      ZN = -ZR*CMPLX(0.0E0,1.0E0)
      IF (YY.GT.0.0E0) GO TO 30
      ZN = CONJG(-ZN)
   30 CONTINUE
      CALL CUNHJ(ZN, GNU, 1, TOL, PHI, ARG, ZETA1, ZETA2, ASUM, BSUM)
      CZ = -ZETA1 + ZETA2
      AARG = CABS(ARG)
   40 CONTINUE
      IF (KODE.EQ.2) CZ = CZ - ZB
      IF (IKFLG.EQ.2) CZ = -CZ
      APHI = CABS(PHI)
      RCZ = REAL(CZ)
C-----------------------------------------------------------------------
C     OVERFLOW TEST
C-----------------------------------------------------------------------
      IF (RCZ.GT.ELIM) GO TO 170
      IF (RCZ.LT.ALIM) GO TO 50
      RCZ = RCZ + ALOG(APHI)
      IF (IFORM.EQ.2) RCZ = RCZ - 0.25E0*ALOG(AARG) - AIC
      IF (RCZ.GT.ELIM) GO TO 170
      GO TO 100
   50 CONTINUE
C-----------------------------------------------------------------------
C     UNDERFLOW TEST
C-----------------------------------------------------------------------
      IF (RCZ.LT.(-ELIM)) GO TO 60
      IF (RCZ.GT.(-ALIM)) GO TO 100
      RCZ = RCZ + ALOG(APHI)
      IF (IFORM.EQ.2) RCZ = RCZ - 0.25E0*ALOG(AARG) - AIC
      IF (RCZ.GT.(-ELIM)) GO TO 80
   60 CONTINUE
      DO 70 I=1,NN
        Y(I) = CZERO
   70 CONTINUE
      NUF = NN
      RETURN
   80 CONTINUE
      ASCLE = 1.0E+3*R1MACH(1)/TOL
      CZ = CZ + CLOG(PHI)
      IF (IFORM.EQ.1) GO TO 90
      CZ = CZ - CMPLX(0.25E0,0.0E0)*CLOG(ARG) - CMPLX(AIC,0.0E0)
   90 CONTINUE
      AX = EXP(RCZ)/TOL
      AY = AIMAG(CZ)
      CZ = CMPLX(AX,0.0E0)*CMPLX(COS(AY),SIN(AY))
      CALL CUCHK(CZ, NW, ASCLE, TOL)
      IF (NW.EQ.1) GO TO 60
  100 CONTINUE
      IF (IKFLG.EQ.2) RETURN
      IF (N.EQ.1) RETURN
C-----------------------------------------------------------------------
C     SET UNDERFLOWS ON I SEQUENCE
C-----------------------------------------------------------------------
  110 CONTINUE
      GNU = FNU + FLOAT(NN-1)
      IF (IFORM.EQ.2) GO TO 120
      INIT = 0
      CALL CUNIK(ZR, GNU, IKFLG, 1, TOL, INIT, PHI, ZETA1, ZETA2, SUM,
     * CWRK)
      CZ = -ZETA1 + ZETA2
      GO TO 130
  120 CONTINUE
      CALL CUNHJ(ZN, GNU, 1, TOL, PHI, ARG, ZETA1, ZETA2, ASUM, BSUM)
      CZ = -ZETA1 + ZETA2
      AARG = CABS(ARG)
  130 CONTINUE
      IF (KODE.EQ.2) CZ = CZ - ZB
      APHI = CABS(PHI)
      RCZ = REAL(CZ)
      IF (RCZ.LT.(-ELIM)) GO TO 140
      IF (RCZ.GT.(-ALIM)) RETURN
      RCZ = RCZ + ALOG(APHI)
      IF (IFORM.EQ.2) RCZ = RCZ - 0.25E0*ALOG(AARG) - AIC
      IF (RCZ.GT.(-ELIM)) GO TO 150
  140 CONTINUE
      Y(NN) = CZERO
      NN = NN - 1
      NUF = NUF + 1
      IF (NN.EQ.0) RETURN
      GO TO 110
  150 CONTINUE
      ASCLE = 1.0E+3*R1MACH(1)/TOL
      CZ = CZ + CLOG(PHI)
      IF (IFORM.EQ.1) GO TO 160
      CZ = CZ - CMPLX(0.25E0,0.0E0)*CLOG(ARG) - CMPLX(AIC,0.0E0)
  160 CONTINUE
      AX = EXP(RCZ)/TOL
      AY = AIMAG(CZ)
      CZ = CMPLX(AX,0.0E0)*CMPLX(COS(AY),SIN(AY))
      CALL CUCHK(CZ, NW, ASCLE, TOL)
      IF (NW.EQ.1) GO TO 140
      RETURN
  170 CONTINUE
      NUF = -1
      RETURN
      END
      SUBROUTINE CWRSK(ZR, FNU, KODE, N, Y, NZ, CW, TOL, ELIM, ALIM)
C***BEGIN PROLOGUE  CWRSK
C***REFER TO  CBESI,CBESK
C
C     CWRSK COMPUTES THE I BESSEL FUNCTION FOR RE(Z).GE.0.0 BY
C     NORMALIZING THE I FUNCTION RATIOS FROM CRATI BY THE WRONSKIAN
C
C***ROUTINES CALLED  CBKNU,CRATI,R1MACH
C***END PROLOGUE  CWRSK
      COMPLEX CINU, CSCL, CT, CW, C1, C2, RCT, ST, Y, ZR
      REAL ACT, ACW, ALIM, ASCLE, ELIM, FNU, S1, S2, TOL, YY
      INTEGER I, KODE, N, NW, NZ
      DIMENSION Y(N), CW(2)
C-----------------------------------------------------------------------
C     I(FNU+I-1,Z) BY BACKWARD RECURRENCE FOR RATIOS
C     Y(I)=I(FNU+I,Z)/I(FNU+I-1,Z) FROM CRATI NORMALIZED BY THE
C     WRONSKIAN WITH K(FNU,Z) AND K(FNU+1,Z) FROM CBKNU.
C-----------------------------------------------------------------------
      NZ = 0
      CALL CBKNU(ZR, FNU, KODE, 2, CW, NW, TOL, ELIM, ALIM)
      IF (NW.NE.0) GO TO 50
      CALL CRATI(ZR, FNU, N, Y, TOL)
C-----------------------------------------------------------------------
C     RECUR FORWARD ON I(FNU+1,Z) = R(FNU,Z)*I(FNU,Z),
C     R(FNU+J-1,Z)=Y(J),  J=1,...,N
C-----------------------------------------------------------------------
      CINU = CMPLX(1.0E0,0.0E0)
      IF (KODE.EQ.1) GO TO 10
      YY = AIMAG(ZR)
      S1 = COS(YY)
      S2 = SIN(YY)
      CINU = CMPLX(S1,S2)
   10 CONTINUE
C-----------------------------------------------------------------------
C     ON LOW EXPONENT MACHINES THE K FUNCTIONS CAN BE CLOSE TO BOTH
C     THE UNDER AND OVERFLOW LIMITS AND THE NORMALIZATION MUST BE
C     SCALED TO PREVENT OVER OR UNDERFLOW. CUOIK HAS DETERMINED THAT
C     THE RESULT IS ON SCALE.
C-----------------------------------------------------------------------
      ACW = CABS(CW(2))
      ASCLE = 1.0E+3*R1MACH(1)/TOL
      CSCL = CMPLX(1.0E0,0.0E0)
      IF (ACW.GT.ASCLE) GO TO 20
      CSCL = CMPLX(1.0E0/TOL,0.0E0)
      GO TO 30
   20 CONTINUE
      ASCLE = 1.0E0/ASCLE
      IF (ACW.LT.ASCLE) GO TO 30
      CSCL = CMPLX(TOL,0.0E0)
   30 CONTINUE
      C1 = CW(1)*CSCL
      C2 = CW(2)*CSCL
      ST = Y(1)
C-----------------------------------------------------------------------
C     CINU=CINU*(CONJG(CT)/CABS(CT))*(1.0E0/CABS(CT) PREVENTS
C     UNDER- OR OVERFLOW PREMATURELY BY SQUARING CABS(CT)
C-----------------------------------------------------------------------
      CT = ZR*(C2+ST*C1)
      ACT = CABS(CT)
      RCT = CMPLX(1.0E0/ACT,0.0E0)
      CT = CONJG(CT)*RCT
      CINU = CINU*RCT*CT
      Y(1) = CINU*CSCL
      IF (N.EQ.1) RETURN
      DO 40 I=2,N
        CINU = ST*CINU
        ST = Y(I)
        Y(I) = CINU*CSCL
   40 CONTINUE
      RETURN
   50 CONTINUE
      NZ = -1
      IF(NW.EQ.(-2)) NZ=-2
      RETURN
      END
      SUBROUTINE CMLRI(Z, FNU, KODE, N, Y, NZ, TOL)
C***BEGIN PROLOGUE  CMLRI
C***REFER TO  CBESI,CBESK
C
C     CMLRI COMPUTES THE I BESSEL FUNCTION FOR RE(Z).GE.0.0 BY THE
C     MILLER ALGORITHM NORMALIZED BY A NEUMANN SERIES.
C
C***ROUTINES CALLED  GAMLN,R1MACH
C***END PROLOGUE  CMLRI
      COMPLEX CK, CNORM, CONE, CTWO, CZERO, PT, P1, P2, RZ, SUM, Y, Z
      REAL ACK, AK, AP, AT, AZ, BK, FKAP, FKK, FLAM, FNF, FNU, RHO,
     * RHO2, SCLE, TFNF, TOL, TST, X, GAMLN, R1MACH
      INTEGER I, IAZ, IDUM, IFNU, INU, ITIME, K, KK, KM, KODE, M, N
      DIMENSION Y(N)
      DATA CZERO,CONE,CTWO /(0.0E0,0.0E0),(1.0E0,0.0E0),(2.0E0,0.0E0)/
      SCLE = 1.0E+3*R1MACH(1)/TOL
      NZ=0
      AZ = CABS(Z)
      X = REAL(Z)
      IAZ = INT(AZ)
      IFNU = INT(FNU)
      INU = IFNU + N - 1
      AT = FLOAT(IAZ) + 1.0E0
      CK = CMPLX(AT,0.0E0)/Z
      RZ = CTWO/Z
      P1 = CZERO
      P2 = CONE
      ACK = (AT+1.0E0)/AZ
      RHO = ACK + SQRT(ACK*ACK-1.0E0)
      RHO2 = RHO*RHO
      TST = (RHO2+RHO2)/((RHO2-1.0E0)*(RHO-1.0E0))
      TST = TST/TOL
C-----------------------------------------------------------------------
C     COMPUTE RELATIVE TRUNCATION ERROR INDEX FOR SERIES
C-----------------------------------------------------------------------
      AK = AT
      DO 10 I=1,80
        PT = P2
        P2 = P1 - CK*P2
        P1 = PT
        CK = CK + RZ
        AP = CABS(P2)
        IF (AP.GT.TST*AK*AK) GO TO 20
        AK = AK + 1.0E0
   10 CONTINUE
      GO TO 110
   20 CONTINUE
      I = I + 1
      K = 0
      IF (INU.LT.IAZ) GO TO 40
C-----------------------------------------------------------------------
C     COMPUTE RELATIVE TRUNCATION ERROR FOR RATIOS
C-----------------------------------------------------------------------
      P1 = CZERO
      P2 = CONE
      AT = FLOAT(INU) + 1.0E0
      CK = CMPLX(AT,0.0E0)/Z
      ACK = AT/AZ
      TST = SQRT(ACK/TOL)
      ITIME = 1
      DO 30 K=1,80
        PT = P2
        P2 = P1 - CK*P2
        P1 = PT
        CK = CK + RZ
        AP = CABS(P2)
        IF (AP.LT.TST) GO TO 30
        IF (ITIME.EQ.2) GO TO 40
        ACK = CABS(CK)
        FLAM = ACK + SQRT(ACK*ACK-1.0E0)
        FKAP = AP/CABS(P1)
        RHO = AMIN1(FLAM,FKAP)
        TST = TST*SQRT(RHO/(RHO*RHO-1.0E0))
        ITIME = 2
   30 CONTINUE
      GO TO 110
   40 CONTINUE
C-----------------------------------------------------------------------
C     BACKWARD RECURRENCE AND SUM NORMALIZING RELATION
C-----------------------------------------------------------------------
      K = K + 1
      KK = MAX0(I+IAZ,K+INU)
      FKK = FLOAT(KK)
      P1 = CZERO
C-----------------------------------------------------------------------
C     SCALE P2 AND SUM BY SCLE
C-----------------------------------------------------------------------
      P2 = CMPLX(SCLE,0.0E0)
      FNF = FNU - FLOAT(IFNU)
      TFNF = FNF + FNF
      BK = GAMLN(FKK+TFNF+1.0E0,IDUM) - GAMLN(FKK+1.0E0,IDUM)
     *     -GAMLN(TFNF+1.0E0,IDUM)
      BK = EXP(BK)
      SUM = CZERO
      KM = KK - INU
      DO 50 I=1,KM
        PT = P2
        P2 = P1 + CMPLX(FKK+FNF,0.0E0)*RZ*P2
        P1 = PT
        AK = 1.0E0 - TFNF/(FKK+TFNF)
        ACK = BK*AK
        SUM = SUM + CMPLX(ACK+BK,0.0E0)*P1
        BK = ACK
        FKK = FKK - 1.0E0
   50 CONTINUE
      Y(N) = P2
      IF (N.EQ.1) GO TO 70
      DO 60 I=2,N
        PT = P2
        P2 = P1 + CMPLX(FKK+FNF,0.0E0)*RZ*P2
        P1 = PT
        AK = 1.0E0 - TFNF/(FKK+TFNF)
        ACK = BK*AK
        SUM = SUM + CMPLX(ACK+BK,0.0E0)*P1
        BK = ACK
        FKK = FKK - 1.0E0
        M = N - I + 1
        Y(M) = P2
   60 CONTINUE
   70 CONTINUE
      IF (IFNU.LE.0) GO TO 90
      DO 80 I=1,IFNU
        PT = P2
        P2 = P1 + CMPLX(FKK+FNF,0.0E0)*RZ*P2
        P1 = PT
        AK = 1.0E0 - TFNF/(FKK+TFNF)
        ACK = BK*AK
        SUM = SUM + CMPLX(ACK+BK,0.0E0)*P1
        BK = ACK
        FKK = FKK - 1.0E0
   80 CONTINUE
   90 CONTINUE
      PT = Z
      IF (KODE.EQ.2) PT = PT - CMPLX(X,0.0E0)
      P1 = -CMPLX(FNF,0.0E0)*CLOG(RZ) + PT
      AP = GAMLN(1.0E0+FNF,IDUM)
      PT = P1 - CMPLX(AP,0.0E0)
C-----------------------------------------------------------------------
C     THE DIVISION CEXP(PT)/(SUM+P2) IS ALTERED TO AVOID OVERFLOW
C     IN THE DENOMINATOR BY SQUARING LARGE QUANTITIES
C-----------------------------------------------------------------------
      P2 = P2 + SUM
      AP = CABS(P2)
      P1 = CMPLX(1.0E0/AP,0.0E0)
      CK = CEXP(PT)*P1
      PT = CONJG(P2)*P1
      CNORM = CK*PT
      DO 100 I=1,N
        Y(I) = Y(I)*CNORM
  100 CONTINUE
      RETURN
  110 CONTINUE
      NZ=-2
      RETURN
      END
      SUBROUTINE CUNHJ(Z, FNU, IPMTR, TOL, PHI, ARG, ZETA1, ZETA2,
     * ASUM, BSUM)
C***BEGIN PROLOGUE  CUNHJ
C***REFER TO  CBESI,CBESK
C
C     REFERENCES
C         HANDBOOK OF MATHEMATICAL FUNCTIONS BY M. ABRAMOWITZ AND I.A.
C         STEGUN, AMS55, NATIONAL BUREAU OF STANDARDS, 1965, CHAPTER 9.
C
C         ASYMPTOTICS AND SPECIAL FUNCTIONS BY F.W.J. OLVER, ACADEMIC
C         PRESS, N.Y., 1974, PAGE 420
C
C     ABSTRACT
C         CUNHJ COMPUTES PARAMETERS FOR BESSEL FUNCTIONS C(FNU,Z) =
C         J(FNU,Z), Y(FNU,Z) OR H(I,FNU,Z) I=1,2 FOR LARGE ORDERS FNU
C         BY MEANS OF THE UNIFORM ASYMPTOTIC EXPANSION
C
C         C(FNU,Z)=C1*PHI*( ASUM*AIRY(ARG) + C2*BSUM*DAIRY(ARG) )
C
C         FOR PROPER CHOICES OF C1, C2, AIRY AND DAIRY WHERE AIRY IS
C         AN AIRY FUNCTION AND DAIRY IS ITS DERIVATIVE.
C
C               (2/3)*FNU*ZETA**1.5 = ZETA1-ZETA2,
C
C         ZETA1=0.5*FNU*CLOG((1+W)/(1-W)), ZETA2=FNU*W FOR SCALING
C         PURPOSES IN AIRY FUNCTIONS FROM CAIRY OR CBIRY.
C
C         MCONJ=SIGN OF AIMAG(Z), BUT IS AMBIGUOUS WHEN Z IS REAL AND
C         MUST BE SPECIFIED. IPMTR=0 RETURNS ALL PARAMETERS. IPMTR=
C         1 COMPUTES ALL EXCEPT ASUM AND BSUM.
C
C***ROUTINES CALLED  (NONE)
C***END PROLOGUE  CUNHJ
      COMPLEX ARG, ASUM, BSUM, CFNU, CONE, CR, CZERO, DR, P, PHI,
     * PRZTH, PTFN, RFN13, RTZTA, RZTH, SUMA, SUMB, TFN, T2, UP, W, W2,
     * Z, ZA, ZB, ZC, ZETA, ZETA1, ZETA2, ZTH
      REAL ALFA, ANG, AP, AR, ATOL, AW2, AZTH, BETA, BR, BTOL, C, EX1,
     * EX2, FNU, FN13, FN23, GAMA, HPI, PI, PP, RFNU, RFNU2, THPI, TOL,
     * WI, WR, ZCI, ZCR, ZETAI, ZETAR, ZTHI, ZTHR, ASUMR, ASUMI, BSUMR,
     * BSUMI, TEST, TSTR, TSTI, AC
      INTEGER IAS, IBS, IPMTR, IS, J, JR, JU, K, KMAX, KP1, KS, L, LR,
     * LRP1, L1, L2, M
      DIMENSION AR(14), BR(14), C(105), ALFA(180), BETA(210), GAMA(30),
     * AP(30), P(30), UP(14), CR(14), DR(14)
      DATA AR(1), AR(2), AR(3), AR(4), AR(5), AR(6), AR(7), AR(8),
     1     AR(9), AR(10), AR(11), AR(12), AR(13), AR(14)/
     2     1.00000000000000000E+00,     1.04166666666666667E-01,
     3     8.35503472222222222E-02,     1.28226574556327160E-01,
     4     2.91849026464140464E-01,     8.81627267443757652E-01,
     5     3.32140828186276754E+00,     1.49957629868625547E+01,
     6     7.89230130115865181E+01,     4.74451538868264323E+02,
     7     3.20749009089066193E+03,     2.40865496408740049E+04,
     8     1.98923119169509794E+05,     1.79190200777534383E+06/
      DATA BR(1), BR(2), BR(3), BR(4), BR(5), BR(6), BR(7), BR(8),
     1     BR(9), BR(10), BR(11), BR(12), BR(13), BR(14)/
     2     1.00000000000000000E+00,    -1.45833333333333333E-01,
     3    -9.87413194444444444E-02,    -1.43312053915895062E-01,
     4    -3.17227202678413548E-01,    -9.42429147957120249E-01,
     5    -3.51120304082635426E+00,    -1.57272636203680451E+01,
     6    -8.22814390971859444E+01,    -4.92355370523670524E+02,
     7    -3.31621856854797251E+03,    -2.48276742452085896E+04,
     8    -2.04526587315129788E+05,    -1.83844491706820990E+06/
      DATA C(1), C(2), C(3), C(4), C(5), C(6), C(7), C(8), C(9), C(10),
     1     C(11), C(12), C(13), C(14), C(15), C(16), C(17), C(18),
     2     C(19), C(20), C(21), C(22), C(23), C(24)/
     3     1.00000000000000000E+00,    -2.08333333333333333E-01,
     4     1.25000000000000000E-01,     3.34201388888888889E-01,
     5    -4.01041666666666667E-01,     7.03125000000000000E-02,
     6    -1.02581259645061728E+00,     1.84646267361111111E+00,
     7    -8.91210937500000000E-01,     7.32421875000000000E-02,
     8     4.66958442342624743E+00,    -1.12070026162229938E+01,
     9     8.78912353515625000E+00,    -2.36408691406250000E+00,
     A     1.12152099609375000E-01,    -2.82120725582002449E+01,
     B     8.46362176746007346E+01,    -9.18182415432400174E+01,
     C     4.25349987453884549E+01,    -7.36879435947963170E+00,
     D     2.27108001708984375E-01,     2.12570130039217123E+02,
     E    -7.65252468141181642E+02,     1.05999045252799988E+03/
      DATA C(25), C(26), C(27), C(28), C(29), C(30), C(31), C(32),
     1     C(33), C(34), C(35), C(36), C(37), C(38), C(39), C(40),
     2     C(41), C(42), C(43), C(44), C(45), C(46), C(47), C(48)/
     3    -6.99579627376132541E+02,     2.18190511744211590E+02,
     4    -2.64914304869515555E+01,     5.72501420974731445E-01,
     5    -1.91945766231840700E+03,     8.06172218173730938E+03,
     6    -1.35865500064341374E+04,     1.16553933368645332E+04,
     7    -5.30564697861340311E+03,     1.20090291321635246E+03,
     8    -1.08090919788394656E+02,     1.72772750258445740E+00,
     9     2.02042913309661486E+04,    -9.69805983886375135E+04,
     A     1.92547001232531532E+05,    -2.03400177280415534E+05,
     B     1.22200464983017460E+05,    -4.11926549688975513E+04,
     C     7.10951430248936372E+03,    -4.93915304773088012E+02,
     D     6.07404200127348304E+00,    -2.42919187900551333E+05,
     E     1.31176361466297720E+06,    -2.99801591853810675E+06/
      DATA C(49), C(50), C(51), C(52), C(53), C(54), C(55), C(56),
     1     C(57), C(58), C(59), C(60), C(61), C(62), C(63), C(64),
     2     C(65), C(66), C(67), C(68), C(69), C(70), C(71), C(72)/
     3     3.76327129765640400E+06,    -2.81356322658653411E+06,
     4     1.26836527332162478E+06,    -3.31645172484563578E+05,
     5     4.52187689813627263E+04,    -2.49983048181120962E+03,
     6     2.43805296995560639E+01,     3.28446985307203782E+06,
     7    -1.97068191184322269E+07,     5.09526024926646422E+07,
     8    -7.41051482115326577E+07,     6.63445122747290267E+07,
     9    -3.75671766607633513E+07,     1.32887671664218183E+07,
     A    -2.78561812808645469E+06,     3.08186404612662398E+05,
     B    -1.38860897537170405E+04,     1.10017140269246738E+02,
     C    -4.93292536645099620E+07,     3.25573074185765749E+08,
     D    -9.39462359681578403E+08,     1.55359689957058006E+09,
     E    -1.62108055210833708E+09,     1.10684281682301447E+09/
      DATA C(73), C(74), C(75), C(76), C(77), C(78), C(79), C(80),
     1     C(81), C(82), C(83), C(84), C(85), C(86), C(87), C(88),
     2     C(89), C(90), C(91), C(92), C(93), C(94), C(95), C(96)/
     3    -4.95889784275030309E+08,     1.42062907797533095E+08,
     4    -2.44740627257387285E+07,     2.24376817792244943E+06,
     5    -8.40054336030240853E+04,     5.51335896122020586E+02,
     6     8.14789096118312115E+08,    -5.86648149205184723E+09,
     7     1.86882075092958249E+10,    -3.46320433881587779E+10,
     8     4.12801855797539740E+10,    -3.30265997498007231E+10,
     9     1.79542137311556001E+10,    -6.56329379261928433E+09,
     A     1.55927986487925751E+09,    -2.25105661889415278E+08,
     B     1.73951075539781645E+07,    -5.49842327572288687E+05,
     C     3.03809051092238427E+03,    -1.46792612476956167E+10,
     D     1.14498237732025810E+11,    -3.99096175224466498E+11,
     E     8.19218669548577329E+11,    -1.09837515608122331E+12/
      DATA C(97), C(98), C(99), C(100), C(101), C(102), C(103), C(104),
     1     C(105)/
     2     1.00815810686538209E+12,    -6.45364869245376503E+11,
     3     2.87900649906150589E+11,    -8.78670721780232657E+10,
     4     1.76347306068349694E+10,    -2.16716498322379509E+09,
     5     1.43157876718888981E+08,    -3.87183344257261262E+06,
     6     1.82577554742931747E+04/
      DATA ALFA(1), ALFA(2), ALFA(3), ALFA(4), ALFA(5), ALFA(6),
     1     ALFA(7), ALFA(8), ALFA(9), ALFA(10), ALFA(11), ALFA(12),
     2     ALFA(13), ALFA(14), ALFA(15), ALFA(16), ALFA(17), ALFA(18),
     3     ALFA(19), ALFA(20), ALFA(21), ALFA(22)/
     4    -4.44444444444444444E-03,    -9.22077922077922078E-04,
     5    -8.84892884892884893E-05,     1.65927687832449737E-04,
     6     2.46691372741792910E-04,     2.65995589346254780E-04,
     7     2.61824297061500945E-04,     2.48730437344655609E-04,
     8     2.32721040083232098E-04,     2.16362485712365082E-04,
     9     2.00738858762752355E-04,     1.86267636637545172E-04,
     A     1.73060775917876493E-04,     1.61091705929015752E-04,
     B     1.50274774160908134E-04,     1.40503497391269794E-04,
     C     1.31668816545922806E-04,     1.23667445598253261E-04,
     D     1.16405271474737902E-04,     1.09798298372713369E-04,
     E     1.03772410422992823E-04,     9.82626078369363448E-05/
      DATA ALFA(23), ALFA(24), ALFA(25), ALFA(26), ALFA(27), ALFA(28),
     1     ALFA(29), ALFA(30), ALFA(31), ALFA(32), ALFA(33), ALFA(34),
     2     ALFA(35), ALFA(36), ALFA(37), ALFA(38), ALFA(39), ALFA(40),
     3     ALFA(41), ALFA(42), ALFA(43), ALFA(44)/
     4     9.32120517249503256E-05,     8.85710852478711718E-05,
     5     8.42963105715700223E-05,     8.03497548407791151E-05,
     6     7.66981345359207388E-05,     7.33122157481777809E-05,
     7     7.01662625163141333E-05,     6.72375633790160292E-05,
     8     6.93735541354588974E-04,     2.32241745182921654E-04,
     9    -1.41986273556691197E-05,    -1.16444931672048640E-04,
     A    -1.50803558053048762E-04,    -1.55121924918096223E-04,
     B    -1.46809756646465549E-04,    -1.33815503867491367E-04,
     C    -1.19744975684254051E-04,    -1.06184319207974020E-04,
     D    -9.37699549891194492E-05,    -8.26923045588193274E-05,
     E    -7.29374348155221211E-05,    -6.44042357721016283E-05/
      DATA ALFA(45), ALFA(46), ALFA(47), ALFA(48), ALFA(49), ALFA(50),
     1     ALFA(51), ALFA(52), ALFA(53), ALFA(54), ALFA(55), ALFA(56),
     2     ALFA(57), ALFA(58), ALFA(59), ALFA(60), ALFA(61), ALFA(62),
     3     ALFA(63), ALFA(64), ALFA(65), ALFA(66)/
     4    -5.69611566009369048E-05,    -5.04731044303561628E-05,
     5    -4.48134868008882786E-05,    -3.98688727717598864E-05,
     6    -3.55400532972042498E-05,    -3.17414256609022480E-05,
     7    -2.83996793904174811E-05,    -2.54522720634870566E-05,
     8    -2.28459297164724555E-05,    -2.05352753106480604E-05,
     9    -1.84816217627666085E-05,    -1.66519330021393806E-05,
     A    -1.50179412980119482E-05,    -1.35554031379040526E-05,
     B    -1.22434746473858131E-05,    -1.10641884811308169E-05,
     C    -3.54211971457743841E-04,    -1.56161263945159416E-04,
     D     3.04465503594936410E-05,     1.30198655773242693E-04,
     E     1.67471106699712269E-04,     1.70222587683592569E-04/
      DATA ALFA(67), ALFA(68), ALFA(69), ALFA(70), ALFA(71), ALFA(72),
     1     ALFA(73), ALFA(74), ALFA(75), ALFA(76), ALFA(77), ALFA(78),
     2     ALFA(79), ALFA(80), ALFA(81), ALFA(82), ALFA(83), ALFA(84),
     3     ALFA(85), ALFA(86), ALFA(87), ALFA(88)/
     4     1.56501427608594704E-04,     1.36339170977445120E-04,
     5     1.14886692029825128E-04,     9.45869093034688111E-05,
     6     7.64498419250898258E-05,     6.07570334965197354E-05,
     7     4.74394299290508799E-05,     3.62757512005344297E-05,
     8     2.69939714979224901E-05,     1.93210938247939253E-05,
     9     1.30056674793963203E-05,     7.82620866744496661E-06,
     A     3.59257485819351583E-06,     1.44040049814251817E-07,
     B    -2.65396769697939116E-06,    -4.91346867098485910E-06,
     C    -6.72739296091248287E-06,    -8.17269379678657923E-06,
     D    -9.31304715093561232E-06,    -1.02011418798016441E-05,
     E    -1.08805962510592880E-05,    -1.13875481509603555E-05/
      DATA ALFA(89), ALFA(90), ALFA(91), ALFA(92), ALFA(93), ALFA(94),
     1     ALFA(95), ALFA(96), ALFA(97), ALFA(98), ALFA(99), ALFA(100),
     2     ALFA(101), ALFA(102), ALFA(103), ALFA(104), ALFA(105),
     3     ALFA(106), ALFA(107), ALFA(108), ALFA(109), ALFA(110)/
     4    -1.17519675674556414E-05,    -1.19987364870944141E-05,
     5     3.78194199201772914E-04,     2.02471952761816167E-04,
     6    -6.37938506318862408E-05,    -2.38598230603005903E-04,
     7    -3.10916256027361568E-04,    -3.13680115247576316E-04,
     8    -2.78950273791323387E-04,    -2.28564082619141374E-04,
     9    -1.75245280340846749E-04,    -1.25544063060690348E-04,
     A    -8.22982872820208365E-05,    -4.62860730588116458E-05,
     B    -1.72334302366962267E-05,     5.60690482304602267E-06,
     C     2.31395443148286800E-05,     3.62642745856793957E-05,
     D     4.58006124490188752E-05,     5.24595294959114050E-05,
     E     5.68396208545815266E-05,     5.94349820393104052E-05/
      DATA ALFA(111), ALFA(112), ALFA(113), ALFA(114), ALFA(115),
     1     ALFA(116), ALFA(117), ALFA(118), ALFA(119), ALFA(120),
     2     ALFA(121), ALFA(122), ALFA(123), ALFA(124), ALFA(125),
     3     ALFA(126), ALFA(127), ALFA(128), ALFA(129), ALFA(130)/
     4     6.06478527578421742E-05,     6.08023907788436497E-05,
     5     6.01577894539460388E-05,     5.89199657344698500E-05,
     6     5.72515823777593053E-05,     5.52804375585852577E-05,
     7     5.31063773802880170E-05,     5.08069302012325706E-05,
     8     4.84418647620094842E-05,     4.60568581607475370E-05,
     9    -6.91141397288294174E-04,    -4.29976633058871912E-04,
     A     1.83067735980039018E-04,     6.60088147542014144E-04,
     B     8.75964969951185931E-04,     8.77335235958235514E-04,
     C     7.49369585378990637E-04,     5.63832329756980918E-04,
     D     3.68059319971443156E-04,     1.88464535514455599E-04/
      DATA ALFA(131), ALFA(132), ALFA(133), ALFA(134), ALFA(135),
     1     ALFA(136), ALFA(137), ALFA(138), ALFA(139), ALFA(140),
     2     ALFA(141), ALFA(142), ALFA(143), ALFA(144), ALFA(145),
     3     ALFA(146), ALFA(147), ALFA(148), ALFA(149), ALFA(150)/
     4     3.70663057664904149E-05,    -8.28520220232137023E-05,
     5    -1.72751952869172998E-04,    -2.36314873605872983E-04,
     6    -2.77966150694906658E-04,    -3.02079514155456919E-04,
     7    -3.12594712643820127E-04,    -3.12872558758067163E-04,
     8    -3.05678038466324377E-04,    -2.93226470614557331E-04,
     9    -2.77255655582934777E-04,    -2.59103928467031709E-04,
     A    -2.39784014396480342E-04,    -2.20048260045422848E-04,
     B    -2.00443911094971498E-04,    -1.81358692210970687E-04,
     C    -1.63057674478657464E-04,    -1.45712672175205844E-04,
     D    -1.29425421983924587E-04,    -1.14245691942445952E-04/
      DATA ALFA(151), ALFA(152), ALFA(153), ALFA(154), ALFA(155),
     1     ALFA(156), ALFA(157), ALFA(158), ALFA(159), ALFA(160),
     2     ALFA(161), ALFA(162), ALFA(163), ALFA(164), ALFA(165),
     3     ALFA(166), ALFA(167), ALFA(168), ALFA(169), ALFA(170)/
     4     1.92821964248775885E-03,     1.35592576302022234E-03,
     5    -7.17858090421302995E-04,    -2.58084802575270346E-03,
     6    -3.49271130826168475E-03,    -3.46986299340960628E-03,
     7    -2.82285233351310182E-03,    -1.88103076404891354E-03,
     8    -8.89531718383947600E-04,     3.87912102631035228E-06,
     9     7.28688540119691412E-04,     1.26566373053457758E-03,
     A     1.62518158372674427E-03,     1.83203153216373172E-03,
     B     1.91588388990527909E-03,     1.90588846755546138E-03,
     C     1.82798982421825727E-03,     1.70389506421121530E-03,
     D     1.55097127171097686E-03,     1.38261421852276159E-03/
      DATA ALFA(171), ALFA(172), ALFA(173), ALFA(174), ALFA(175),
     1     ALFA(176), ALFA(177), ALFA(178), ALFA(179), ALFA(180)/
     2     1.20881424230064774E-03,     1.03676532638344962E-03,
     3     8.71437918068619115E-04,     7.16080155297701002E-04,
     4     5.72637002558129372E-04,     4.42089819465802277E-04,
     5     3.24724948503090564E-04,     2.20342042730246599E-04,
     6     1.28412898401353882E-04,     4.82005924552095464E-05/
      DATA BETA(1), BETA(2), BETA(3), BETA(4), BETA(5), BETA(6),
     1     BETA(7), BETA(8), BETA(9), BETA(10), BETA(11), BETA(12),
     2     BETA(13), BETA(14), BETA(15), BETA(16), BETA(17), BETA(18),
     3     BETA(19), BETA(20), BETA(21), BETA(22)/
     4     1.79988721413553309E-02,     5.59964911064388073E-03,
     5     2.88501402231132779E-03,     1.80096606761053941E-03,
     6     1.24753110589199202E-03,     9.22878876572938311E-04,
     7     7.14430421727287357E-04,     5.71787281789704872E-04,
     8     4.69431007606481533E-04,     3.93232835462916638E-04,
     9     3.34818889318297664E-04,     2.88952148495751517E-04,
     A     2.52211615549573284E-04,     2.22280580798883327E-04,
     B     1.97541838033062524E-04,     1.76836855019718004E-04,
     C     1.59316899661821081E-04,     1.44347930197333986E-04,
     D     1.31448068119965379E-04,     1.20245444949302884E-04,
     E     1.10449144504599392E-04,     1.01828770740567258E-04/
      DATA BETA(23), BETA(24), BETA(25), BETA(26), BETA(27), BETA(28),
     1     BETA(29), BETA(30), BETA(31), BETA(32), BETA(33), BETA(34),
     2     BETA(35), BETA(36), BETA(37), BETA(38), BETA(39), BETA(40),
     3     BETA(41), BETA(42), BETA(43), BETA(44)/
     4     9.41998224204237509E-05,     8.74130545753834437E-05,
     5     8.13466262162801467E-05,     7.59002269646219339E-05,
     6     7.09906300634153481E-05,     6.65482874842468183E-05,
     7     6.25146958969275078E-05,     5.88403394426251749E-05,
     8    -1.49282953213429172E-03,    -8.78204709546389328E-04,
     9    -5.02916549572034614E-04,    -2.94822138512746025E-04,
     A    -1.75463996970782828E-04,    -1.04008550460816434E-04,
     B    -5.96141953046457895E-05,    -3.12038929076098340E-05,
     C    -1.26089735980230047E-05,    -2.42892608575730389E-07,
     D     8.05996165414273571E-06,     1.36507009262147391E-05,
     E     1.73964125472926261E-05,     1.98672978842133780E-05/
      DATA BETA(45), BETA(46), BETA(47), BETA(48), BETA(49), BETA(50),
     1     BETA(51), BETA(52), BETA(53), BETA(54), BETA(55), BETA(56),
     2     BETA(57), BETA(58), BETA(59), BETA(60), BETA(61), BETA(62),
     3     BETA(63), BETA(64), BETA(65), BETA(66)/
     4     2.14463263790822639E-05,     2.23954659232456514E-05,
     5     2.28967783814712629E-05,     2.30785389811177817E-05,
     6     2.30321976080909144E-05,     2.28236073720348722E-05,
     7     2.25005881105292418E-05,     2.20981015361991429E-05,
     8     2.16418427448103905E-05,     2.11507649256220843E-05,
     9     2.06388749782170737E-05,     2.01165241997081666E-05,
     A     1.95913450141179244E-05,     1.90689367910436740E-05,
     B     1.85533719641636667E-05,     1.80475722259674218E-05,
     C     5.52213076721292790E-04,     4.47932581552384646E-04,
     D     2.79520653992020589E-04,     1.52468156198446602E-04,
     E     6.93271105657043598E-05,     1.76258683069991397E-05/
      DATA BETA(67), BETA(68), BETA(69), BETA(70), BETA(71), BETA(72),
     1     BETA(73), BETA(74), BETA(75), BETA(76), BETA(77), BETA(78),
     2     BETA(79), BETA(80), BETA(81), BETA(82), BETA(83), BETA(84),
     3     BETA(85), BETA(86), BETA(87), BETA(88)/
     4    -1.35744996343269136E-05,    -3.17972413350427135E-05,
     5    -4.18861861696693365E-05,    -4.69004889379141029E-05,
     6    -4.87665447413787352E-05,    -4.87010031186735069E-05,
     7    -4.74755620890086638E-05,    -4.55813058138628452E-05,
     8    -4.33309644511266036E-05,    -4.09230193157750364E-05,
     9    -3.84822638603221274E-05,    -3.60857167535410501E-05,
     A    -3.37793306123367417E-05,    -3.15888560772109621E-05,
     B    -2.95269561750807315E-05,    -2.75978914828335759E-05,
     C    -2.58006174666883713E-05,    -2.41308356761280200E-05,
     D    -2.25823509518346033E-05,    -2.11479656768912971E-05,
     E    -1.98200638885294927E-05,    -1.85909870801065077E-05/
      DATA BETA(89), BETA(90), BETA(91), BETA(92), BETA(93), BETA(94),
     1     BETA(95), BETA(96), BETA(97), BETA(98), BETA(99), BETA(100),
     2     BETA(101), BETA(102), BETA(103), BETA(104), BETA(105),
     3     BETA(106), BETA(107), BETA(108), BETA(109), BETA(110)/
     4    -1.74532699844210224E-05,    -1.63997823854497997E-05,
     5    -4.74617796559959808E-04,    -4.77864567147321487E-04,
     6    -3.20390228067037603E-04,    -1.61105016119962282E-04,
     7    -4.25778101285435204E-05,     3.44571294294967503E-05,
     8     7.97092684075674924E-05,     1.03138236708272200E-04,
     9     1.12466775262204158E-04,     1.13103642108481389E-04,
     A     1.08651634848774268E-04,     1.01437951597661973E-04,
     B     9.29298396593363896E-05,     8.40293133016089978E-05,
     C     7.52727991349134062E-05,     6.69632521975730872E-05,
     D     5.92564547323194704E-05,     5.22169308826975567E-05,
     E     4.58539485165360646E-05,     4.01445513891486808E-05/
      DATA BETA(111), BETA(112), BETA(113), BETA(114), BETA(115),
     1     BETA(116), BETA(117), BETA(118), BETA(119), BETA(120),
     2     BETA(121), BETA(122), BETA(123), BETA(124), BETA(125),
     3     BETA(126), BETA(127), BETA(128), BETA(129), BETA(130)/
     4     3.50481730031328081E-05,     3.05157995034346659E-05,
     5     2.64956119950516039E-05,     2.29363633690998152E-05,
     6     1.97893056664021636E-05,     1.70091984636412623E-05,
     7     1.45547428261524004E-05,     1.23886640995878413E-05,
     8     1.04775876076583236E-05,     8.79179954978479373E-06,
     9     7.36465810572578444E-04,     8.72790805146193976E-04,
     A     6.22614862573135066E-04,     2.85998154194304147E-04,
     B     3.84737672879366102E-06,    -1.87906003636971558E-04,
     C    -2.97603646594554535E-04,    -3.45998126832656348E-04,
     D    -3.53382470916037712E-04,    -3.35715635775048757E-04/
      DATA BETA(131), BETA(132), BETA(133), BETA(134), BETA(135),
     1     BETA(136), BETA(137), BETA(138), BETA(139), BETA(140),
     2     BETA(141), BETA(142), BETA(143), BETA(144), BETA(145),
     3     BETA(146), BETA(147), BETA(148), BETA(149), BETA(150)/
     4    -3.04321124789039809E-04,    -2.66722723047612821E-04,
     5    -2.27654214122819527E-04,    -1.89922611854562356E-04,
     6    -1.55058918599093870E-04,    -1.23778240761873630E-04,
     7    -9.62926147717644187E-05,    -7.25178327714425337E-05,
     8    -5.22070028895633801E-05,    -3.50347750511900522E-05,
     9    -2.06489761035551757E-05,    -8.70106096849767054E-06,
     A     1.13698686675100290E-06,     9.16426474122778849E-06,
     B     1.56477785428872620E-05,     2.08223629482466847E-05,
     C     2.48923381004595156E-05,     2.80340509574146325E-05,
     D     3.03987774629861915E-05,     3.21156731406700616E-05/
      DATA BETA(151), BETA(152), BETA(153), BETA(154), BETA(155),
     1     BETA(156), BETA(157), BETA(158), BETA(159), BETA(160),
     2     BETA(161), BETA(162), BETA(163), BETA(164), BETA(165),
     3     BETA(166), BETA(167), BETA(168), BETA(169), BETA(170)/
     4    -1.80182191963885708E-03,    -2.43402962938042533E-03,
     5    -1.83422663549856802E-03,    -7.62204596354009765E-04,
     6     2.39079475256927218E-04,     9.49266117176881141E-04,
     7     1.34467449701540359E-03,     1.48457495259449178E-03,
     8     1.44732339830617591E-03,     1.30268261285657186E-03,
     9     1.10351597375642682E-03,     8.86047440419791759E-04,
     A     6.73073208165665473E-04,     4.77603872856582378E-04,
     B     3.05991926358789362E-04,     1.60315694594721630E-04,
     C     4.00749555270613286E-05,    -5.66607461635251611E-05,
     D    -1.32506186772982638E-04,    -1.90296187989614057E-04/
      DATA BETA(171), BETA(172), BETA(173), BETA(174), BETA(175),
     1     BETA(176), BETA(177), BETA(178), BETA(179), BETA(180),
     2     BETA(181), BETA(182), BETA(183), BETA(184), BETA(185),
     3     BETA(186), BETA(187), BETA(188), BETA(189), BETA(190)/
     4    -2.32811450376937408E-04,    -2.62628811464668841E-04,
     5    -2.82050469867598672E-04,    -2.93081563192861167E-04,
     6    -2.97435962176316616E-04,    -2.96557334239348078E-04,
     7    -2.91647363312090861E-04,    -2.83696203837734166E-04,
     8    -2.73512317095673346E-04,    -2.61750155806768580E-04,
     9     6.38585891212050914E-03,     9.62374215806377941E-03,
     A     7.61878061207001043E-03,     2.83219055545628054E-03,
     B    -2.09841352012720090E-03,    -5.73826764216626498E-03,
     C    -7.70804244495414620E-03,    -8.21011692264844401E-03,
     D    -7.65824520346905413E-03,    -6.47209729391045177E-03/
      DATA BETA(191), BETA(192), BETA(193), BETA(194), BETA(195),
     1     BETA(196), BETA(197), BETA(198), BETA(199), BETA(200),
     2     BETA(201), BETA(202), BETA(203), BETA(204), BETA(205),
     3     BETA(206), BETA(207), BETA(208), BETA(209), BETA(210)/
     4    -4.99132412004966473E-03,    -3.45612289713133280E-03,
     5    -2.01785580014170775E-03,    -7.59430686781961401E-04,
     6     2.84173631523859138E-04,     1.10891667586337403E-03,
     7     1.72901493872728771E-03,     2.16812590802684701E-03,
     8     2.45357710494539735E-03,     2.61281821058334862E-03,
     9     2.67141039656276912E-03,     2.65203073395980430E-03,
     A     2.57411652877287315E-03,     2.45389126236094427E-03,
     B     2.30460058071795494E-03,     2.13684837686712662E-03,
     C     1.95896528478870911E-03,     1.77737008679454412E-03,
     D     1.59690280765839059E-03,     1.42111975664438546E-03/
      DATA GAMA(1), GAMA(2), GAMA(3), GAMA(4), GAMA(5), GAMA(6),
     1     GAMA(7), GAMA(8), GAMA(9), GAMA(10), GAMA(11), GAMA(12),
     2     GAMA(13), GAMA(14), GAMA(15), GAMA(16), GAMA(17), GAMA(18),
     3     GAMA(19), GAMA(20), GAMA(21), GAMA(22)/
     4     6.29960524947436582E-01,     2.51984209978974633E-01,
     5     1.54790300415655846E-01,     1.10713062416159013E-01,
     6     8.57309395527394825E-02,     6.97161316958684292E-02,
     7     5.86085671893713576E-02,     5.04698873536310685E-02,
     8     4.42600580689154809E-02,     3.93720661543509966E-02,
     9     3.54283195924455368E-02,     3.21818857502098231E-02,
     A     2.94646240791157679E-02,     2.71581677112934479E-02,
     B     2.51768272973861779E-02,     2.34570755306078891E-02,
     C     2.19508390134907203E-02,     2.06210828235646240E-02,
     D     1.94388240897880846E-02,     1.83810633800683158E-02,
     E     1.74293213231963172E-02,     1.65685837786612353E-02/
      DATA GAMA(23), GAMA(24), GAMA(25), GAMA(26), GAMA(27), GAMA(28),
     1     GAMA(29), GAMA(30)/
     2     1.57865285987918445E-02,     1.50729501494095594E-02,
     3     1.44193250839954639E-02,     1.38184805735341786E-02,
     4     1.32643378994276568E-02,     1.27517121970498651E-02,
     5     1.22761545318762767E-02,     1.18338262398482403E-02/
      DATA EX1, EX2, HPI, PI, THPI /
     1     3.33333333333333333E-01,     6.66666666666666667E-01,
     2     1.57079632679489662E+00,     3.14159265358979324E+00,
     3     4.71238898038468986E+00/
      DATA CZERO, CONE / (0.0E0,0.0E0), (1.0E0,0.0E0) /
C
      RFNU = 1.0E0/FNU
C     ZB = Z*CMPLX(RFNU,0.0E0)
C-----------------------------------------------------------------------
C     OVERFLOW TEST (Z/FNU TOO SMALL)
C-----------------------------------------------------------------------
      TSTR = REAL(Z)
      TSTI = AIMAG(Z)
      TEST = R1MACH(1)*1.0E+3
      AC = FNU*TEST
      IF (ABS(TSTR).GT.AC .OR. ABS(TSTI).GT.AC) GO TO 15
      AC = 2.0E0*ABS(ALOG(TEST))+FNU
      ZETA1 = CMPLX(AC,0.0E0)
      ZETA2 = CMPLX(FNU,0.0E0)
      PHI=CONE
      ARG=CONE
      RETURN
   15 CONTINUE
      ZB = Z*CMPLX(RFNU,0.0E0)
      RFNU2 = RFNU*RFNU
C-----------------------------------------------------------------------
C     COMPUTE IN THE FOURTH QUADRANT
C-----------------------------------------------------------------------
      FN13 = FNU**EX1
      FN23 = FN13*FN13
      RFN13 = CMPLX(1.0E0/FN13,0.0E0)
      W2 = CONE - ZB*ZB
      AW2 = CABS(W2)
      IF (AW2.GT.0.25E0) GO TO 130
C-----------------------------------------------------------------------
C     POWER SERIES FOR CABS(W2).LE.0.25E0
C-----------------------------------------------------------------------
      K = 1
      P(1) = CONE
      SUMA = CMPLX(GAMA(1),0.0E0)
      AP(1) = 1.0E0
      IF (AW2.LT.TOL) GO TO 20
      DO 10 K=2,30
        P(K) = P(K-1)*W2
        SUMA = SUMA + P(K)*CMPLX(GAMA(K),0.0E0)
        AP(K) = AP(K-1)*AW2
        IF (AP(K).LT.TOL) GO TO 20
   10 CONTINUE
      K = 30
   20 CONTINUE
      KMAX = K
      ZETA = W2*SUMA
      ARG = ZETA*CMPLX(FN23,0.0E0)
      ZA = CSQRT(SUMA)
      ZETA2 = CSQRT(W2)*CMPLX(FNU,0.0E0)
      ZETA1 = ZETA2*(CONE+ZETA*ZA*CMPLX(EX2,0.0E0))
      ZA = ZA + ZA
      PHI = CSQRT(ZA)*RFN13
      IF (IPMTR.EQ.1) GO TO 120
C-----------------------------------------------------------------------
C     SUM SERIES FOR ASUM AND BSUM
C-----------------------------------------------------------------------
      SUMB = CZERO
      DO 30 K=1,KMAX
        SUMB = SUMB + P(K)*CMPLX(BETA(K),0.0E0)
   30 CONTINUE
      ASUM = CZERO
      BSUM = SUMB
      L1 = 0
      L2 = 30
      BTOL = TOL*CABS(BSUM)
      ATOL = TOL
      PP = 1.0E0
      IAS = 0
      IBS = 0
      IF (RFNU2.LT.TOL) GO TO 110
      DO 100 IS=2,7
        ATOL = ATOL/RFNU2
        PP = PP*RFNU2
        IF (IAS.EQ.1) GO TO 60
        SUMA = CZERO
        DO 40 K=1,KMAX
          M = L1 + K
          SUMA = SUMA + P(K)*CMPLX(ALFA(M),0.0E0)
          IF (AP(K).LT.ATOL) GO TO 50
   40   CONTINUE
   50   CONTINUE
        ASUM = ASUM + SUMA*CMPLX(PP,0.0E0)
        IF (PP.LT.TOL) IAS = 1
   60   CONTINUE
        IF (IBS.EQ.1) GO TO 90
        SUMB = CZERO
        DO 70 K=1,KMAX
          M = L2 + K
          SUMB = SUMB + P(K)*CMPLX(BETA(M),0.0E0)
          IF (AP(K).LT.ATOL) GO TO 80
   70   CONTINUE
   80   CONTINUE
        BSUM = BSUM + SUMB*CMPLX(PP,0.0E0)
        IF (PP.LT.BTOL) IBS = 1
   90   CONTINUE
        IF (IAS.EQ.1 .AND. IBS.EQ.1) GO TO 110
        L1 = L1 + 30
        L2 = L2 + 30
  100 CONTINUE
  110 CONTINUE
      ASUM = ASUM + CONE
      PP = RFNU*REAL(RFN13)
      BSUM = BSUM*CMPLX(PP,0.0E0)
  120 CONTINUE
      RETURN
C-----------------------------------------------------------------------
C     CABS(W2).GT.0.25E0
C-----------------------------------------------------------------------
  130 CONTINUE
      W = CSQRT(W2)
      WR = REAL(W)
      WI = AIMAG(W)
      IF (WR.LT.0.0E0) WR = 0.0E0
      IF (WI.LT.0.0E0) WI = 0.0E0
      W = CMPLX(WR,WI)
      ZA = (CONE+W)/ZB
      ZC = CLOG(ZA)
      ZCR = REAL(ZC)
      ZCI = AIMAG(ZC)
      IF (ZCI.LT.0.0E0) ZCI = 0.0E0
      IF (ZCI.GT.HPI) ZCI = HPI
      IF (ZCR.LT.0.0E0) ZCR = 0.0E0
      ZC = CMPLX(ZCR,ZCI)
      ZTH = (ZC-W)*CMPLX(1.5E0,0.0E0)
      CFNU = CMPLX(FNU,0.0E0)
      ZETA1 = ZC*CFNU
      ZETA2 = W*CFNU
      AZTH = CABS(ZTH)
      ZTHR = REAL(ZTH)
      ZTHI = AIMAG(ZTH)
      ANG = THPI
      IF (ZTHR.GE.0.0E0 .AND. ZTHI.LT.0.0E0) GO TO 140
      ANG = HPI
      IF (ZTHR.EQ.0.0E0) GO TO 140
      ANG = ATAN(ZTHI/ZTHR)
      IF (ZTHR.LT.0.0E0) ANG = ANG + PI
  140 CONTINUE
      PP = AZTH**EX2
      ANG = ANG*EX2
      ZETAR = PP*COS(ANG)
      ZETAI = PP*SIN(ANG)
      IF (ZETAI.LT.0.0E0) ZETAI = 0.0E0
      ZETA = CMPLX(ZETAR,ZETAI)
      ARG = ZETA*CMPLX(FN23,0.0E0)
      RTZTA = ZTH/ZETA
      ZA = RTZTA/W
      PHI = CSQRT(ZA+ZA)*RFN13
      IF (IPMTR.EQ.1) GO TO 120
      TFN = CMPLX(RFNU,0.0E0)/W
      RZTH = CMPLX(RFNU,0.0E0)/ZTH
      ZC = RZTH*CMPLX(AR(2),0.0E0)
      T2 = CONE/W2
      UP(2) = (T2*CMPLX(C(2),0.0E0)+CMPLX(C(3),0.0E0))*TFN
      BSUM = UP(2) + ZC
      ASUM = CZERO
      IF (RFNU.LT.TOL) GO TO 220
      PRZTH = RZTH
      PTFN = TFN
      UP(1) = CONE
      PP = 1.0E0
      BSUMR = REAL(BSUM)
      BSUMI = AIMAG(BSUM)
      BTOL = TOL*(ABS(BSUMR)+ABS(BSUMI))
      KS = 0
      KP1 = 2
      L = 3
      IAS = 0
      IBS = 0
      DO 210 LR=2,12,2
        LRP1 = LR + 1
C-----------------------------------------------------------------------
C     COMPUTE TWO ADDITIONAL CR, DR, AND UP FOR TWO MORE TERMS IN
C     NEXT SUMA AND SUMB
C-----------------------------------------------------------------------
        DO 160 K=LR,LRP1
          KS = KS + 1
          KP1 = KP1 + 1
          L = L + 1
          ZA = CMPLX(C(L),0.0E0)
          DO 150 J=2,KP1
            L = L + 1
            ZA = ZA*T2 + CMPLX(C(L),0.0E0)
  150     CONTINUE
          PTFN = PTFN*TFN
          UP(KP1) = PTFN*ZA
          CR(KS) = PRZTH*CMPLX(BR(KS+1),0.0E0)
          PRZTH = PRZTH*RZTH
          DR(KS) = PRZTH*CMPLX(AR(KS+2),0.0E0)
  160   CONTINUE
        PP = PP*RFNU2
        IF (IAS.EQ.1) GO TO 180
        SUMA = UP(LRP1)
        JU = LRP1
        DO 170 JR=1,LR
          JU = JU - 1
          SUMA = SUMA + CR(JR)*UP(JU)
  170   CONTINUE
        ASUM = ASUM + SUMA
        ASUMR = REAL(ASUM)
        ASUMI = AIMAG(ASUM)
        TEST = ABS(ASUMR) + ABS(ASUMI)
        IF (PP.LT.TOL .AND. TEST.LT.TOL) IAS = 1
  180   CONTINUE
        IF (IBS.EQ.1) GO TO 200
        SUMB = UP(LR+2) + UP(LRP1)*ZC
        JU = LRP1
        DO 190 JR=1,LR
          JU = JU - 1
          SUMB = SUMB + DR(JR)*UP(JU)
  190   CONTINUE
        BSUM = BSUM + SUMB
        BSUMR = REAL(BSUM)
        BSUMI = AIMAG(BSUM)
        TEST = ABS(BSUMR) + ABS(BSUMI)
        IF (PP.LT.BTOL .AND. TEST.LT.TOL) IBS = 1
  200   CONTINUE
        IF (IAS.EQ.1 .AND. IBS.EQ.1) GO TO 220
  210 CONTINUE
  220 CONTINUE
      ASUM = ASUM + CONE
      BSUM = -BSUM*RFN13/RTZTA
      GO TO 120
      END
      SUBROUTINE CSERI(Z, FNU, KODE, N, Y, NZ, TOL, ELIM, ALIM)
C***BEGIN PROLOGUE  CSERI
C***REFER TO  CBESI,CBESK
C
C     CSERI COMPUTES THE I BESSEL FUNCTION FOR REAL(Z).GE.0.0 BY
C     MEANS OF THE POWER SERIES FOR LARGE CABS(Z) IN THE
C     REGION CABS(Z).LE.2*SQRT(FNU+1). NZ=0 IS A NORMAL RETURN.
C     NZ.GT.0 MEANS THAT THE LAST NZ COMPONENTS WERE SET TO ZERO
C     DUE TO UNDERFLOW. NZ.LT.0 MEANS UNDERFLOW OCCURRED, BUT THE
C     CONDITION CABS(Z).LE.2*SQRT(FNU+1) WAS VIOLATED AND THE
C     COMPUTATION MUST BE COMPLETED IN ANOTHER ROUTINE WITH N=N-ABS(NZ).
C
C***ROUTINES CALLED  CUCHK,GAMLN,R1MACH
C***END PROLOGUE  CSERI
      COMPLEX AK1, CK, COEF, CONE, CRSC, CZ, CZERO, HZ, RZ, S1, S2, W,
     * Y, Z
      REAL AA, ACZ, AK, ALIM, ARM, ASCLE, ATOL, AZ, DFNU, ELIM, FNU,
     * FNUP, RAK1, RS, RTR1, S, SS, TOL, X, GAMLN, R1MACH
      INTEGER I, IB, IDUM, IFLAG, IL, K, KODE, L, M, N, NN, NW, NZ
      DIMENSION Y(N), W(2)
      DATA CZERO, CONE / (0.0E0,0.0E0), (1.0E0,0.0E0) /
C
      NZ = 0
      AZ = CABS(Z)
      IF (AZ.EQ.0.0E0) GO TO 150
      X = REAL(Z)
      ARM = 1.0E+3*R1MACH(1)
      RTR1 = SQRT(ARM)
      CRSC = CMPLX(1.0E0,0.0E0)
      IFLAG = 0
      IF (AZ.LT.ARM) GO TO 140
      HZ = Z*CMPLX(0.5E0,0.0E0)
      CZ = CZERO
      IF (AZ.GT.RTR1) CZ = HZ*HZ
      ACZ = CABS(CZ)
      NN = N
      CK = CLOG(HZ)
   10 CONTINUE
      DFNU = FNU + FLOAT(NN-1)
      FNUP = DFNU + 1.0E0
C-----------------------------------------------------------------------
C     UNDERFLOW TEST
C-----------------------------------------------------------------------
      AK1 = CK*CMPLX(DFNU,0.0E0)
      AK = GAMLN(FNUP,IDUM)
      AK1 = AK1 - CMPLX(AK,0.0E0)
      IF (KODE.EQ.2) AK1 = AK1 - CMPLX(X,0.0E0)
      RAK1 = REAL(AK1)
      IF (RAK1.GT.(-ELIM)) GO TO 30
   20 CONTINUE
      NZ = NZ + 1
      Y(NN) = CZERO
      IF (ACZ.GT.DFNU) GO TO 170
      NN = NN - 1
      IF (NN.EQ.0) RETURN
      GO TO 10
   30 CONTINUE
      IF (RAK1.GT.(-ALIM)) GO TO 40
      IFLAG = 1
      SS = 1.0E0/TOL
      CRSC = CMPLX(TOL,0.0E0)
      ASCLE = ARM*SS
   40 CONTINUE
      AK = AIMAG(AK1)
      AA = EXP(RAK1)
      IF (IFLAG.EQ.1) AA = AA*SS
      COEF = CMPLX(AA,0.0E0)*CMPLX(COS(AK),SIN(AK))
      ATOL = TOL*ACZ/FNUP
      IL = MIN0(2,NN)
      DO 80 I=1,IL
        DFNU = FNU + FLOAT(NN-I)
        FNUP = DFNU + 1.0E0
        S1 = CONE
        IF (ACZ.LT.TOL*FNUP) GO TO 60
        AK1 = CONE
        AK = FNUP + 2.0E0
        S = FNUP
        AA = 2.0E0
   50   CONTINUE
        RS = 1.0E0/S
        AK1 = AK1*CZ*CMPLX(RS,0.0E0)
        S1 = S1 + AK1
        S = S + AK
        AK = AK + 2.0E0
        AA = AA*ACZ*RS
        IF (AA.GT.ATOL) GO TO 50
   60   CONTINUE
        M = NN - I + 1
        S2 = S1*COEF
        W(I) = S2
        IF (IFLAG.EQ.0) GO TO 70
        CALL CUCHK(S2, NW, ASCLE, TOL)
        IF (NW.NE.0) GO TO 20
   70   CONTINUE
        Y(M) = S2*CRSC
        IF (I.NE.IL) COEF = COEF*CMPLX(DFNU,0.0E0)/HZ
   80 CONTINUE
      IF (NN.LE.2) RETURN
      K = NN - 2
      AK = FLOAT(K)
      RZ = (CONE+CONE)/Z
      IF (IFLAG.EQ.1) GO TO 110
      IB = 3
   90 CONTINUE
      DO 100 I=IB,NN
        Y(K) = CMPLX(AK+FNU,0.0E0)*RZ*Y(K+1) + Y(K+2)
        AK = AK - 1.0E0
        K = K - 1
  100 CONTINUE
      RETURN
C-----------------------------------------------------------------------
C     RECUR BACKWARD WITH SCALED VALUES
C-----------------------------------------------------------------------
  110 CONTINUE
C-----------------------------------------------------------------------
C     EXP(-ALIM)=EXP(-ELIM)/TOL=APPROX. ONE PRECISION ABOVE THE
C     UNDERFLOW LIMIT = ASCLE = R1MACH(1)*CSCL*1.0E+3
C-----------------------------------------------------------------------
      S1 = W(1)
      S2 = W(2)
      DO 120 L=3,NN
        CK = S2
        S2 = S1 + CMPLX(AK+FNU,0.0E0)*RZ*S2
        S1 = CK
        CK = S2*CRSC
        Y(K) = CK
        AK = AK - 1.0E0
        K = K - 1
        IF (CABS(CK).GT.ASCLE) GO TO 130
  120 CONTINUE
      RETURN
  130 CONTINUE
      IB = L + 1
      IF (IB.GT.NN) RETURN
      GO TO 90
  140 CONTINUE
      NZ = N
      IF (FNU.EQ.0.0E0) NZ = NZ - 1
  150 CONTINUE
      Y(1) = CZERO
      IF (FNU.EQ.0.0E0) Y(1) = CONE
      IF (N.EQ.1) RETURN
      DO 160 I=2,N
        Y(I) = CZERO
  160 CONTINUE
      RETURN
C-----------------------------------------------------------------------
C     RETURN WITH NZ.LT.0 IF CABS(Z*Z/4).GT.FNU+N-NZ-1 COMPLETE
C     THE CALCULATION IN CBINU WITH N=N-IABS(NZ)
C-----------------------------------------------------------------------
  170 CONTINUE
      NZ = -NZ
      RETURN
      END
      SUBROUTINE CASYI(Z, FNU, KODE, N, Y, NZ, RL, TOL, ELIM, ALIM)
C***BEGIN PROLOGUE  CASYI
C***REFER TO  CBESI,CBESK
C
C     CASYI COMPUTES THE I BESSEL FUNCTION FOR REAL(Z).GE.0.0 BY
C     MEANS OF THE ASYMPTOTIC EXPANSION FOR LARGE CABS(Z) IN THE
C     REGION CABS(Z).GT.MAX(RL,FNU*FNU/2). NZ=0 IS A NORMAL RETURN.
C     NZ.LT.0 INDICATES AN OVERFLOW ON KODE=1.
C
C***ROUTINES CALLED  R1MACH
C***END PROLOGUE  CASYI
      COMPLEX AK1, CK, CONE, CS1, CS2, CZ, CZERO, DK, EZ, P1, RZ, S2,
     * Y, Z
      REAL AA, ACZ, AEZ, AK, ALIM, ARG, ARM, ATOL, AZ, BB, BK, DFNU,
     * DNU2, ELIM, FDN, FNU, PI, RL, RTPI, RTR1, S, SGN, SQK, TOL, X,
     * YY, R1MACH
      INTEGER I, IB, IL, INU, J, JL, K, KODE, KODED, M, N, NN, NZ
      DIMENSION Y(N)
      DATA PI, RTPI  /3.14159265358979324E0 , 0.159154943091895336E0 /
      DATA CZERO, CONE / (0.0E0,0.0E0), (1.0E0,0.0E0) /
C
      NZ = 0
      AZ = CABS(Z)
      X = REAL(Z)
      ARM = 1.0E+3*R1MACH(1)
      RTR1 = SQRT(ARM)
      IL = MIN0(2,N)
      DFNU = FNU + FLOAT(N-IL)
C-----------------------------------------------------------------------
C     OVERFLOW TEST
C-----------------------------------------------------------------------
      AK1 = CMPLX(RTPI,0.0E0)/Z
      AK1 = CSQRT(AK1)
      CZ = Z
      IF (KODE.EQ.2) CZ = Z - CMPLX(X,0.0E0)
      ACZ = REAL(CZ)
      IF (ABS(ACZ).GT.ELIM) GO TO 80
      DNU2 = DFNU + DFNU
      KODED = 1
      IF ((ABS(ACZ).GT.ALIM) .AND. (N.GT.2)) GO TO 10
      KODED = 0
      AK1 = AK1*CEXP(CZ)
   10 CONTINUE
      FDN = 0.0E0
      IF (DNU2.GT.RTR1) FDN = DNU2*DNU2
      EZ = Z*CMPLX(8.0E0,0.0E0)
C-----------------------------------------------------------------------
C     WHEN Z IS IMAGINARY, THE ERROR TEST MUST BE MADE RELATIVE TO THE
C     FIRST RECIPROCAL POWER SINCE THIS IS THE LEADING TERM OF THE
C     EXPANSION FOR THE IMAGINARY PART.
C-----------------------------------------------------------------------
      AEZ = 8.0E0*AZ
      S = TOL/AEZ
      JL = INT(RL+RL) + 2
      YY = AIMAG(Z)
      P1 = CZERO
      IF (YY.EQ.0.0E0) GO TO 20
C-----------------------------------------------------------------------
C     CALCULATE EXP(PI*(0.5+FNU+N-IL)*I) TO MINIMIZE LOSSES OF
C     SIGNIFICANCE WHEN FNU OR N IS LARGE
C-----------------------------------------------------------------------
      INU = INT(FNU)
      ARG = (FNU-FLOAT(INU))*PI
      INU = INU + N - IL
      AK = -SIN(ARG)
      BK = COS(ARG)
      IF (YY.LT.0.0E0) BK = -BK
      P1 = CMPLX(AK,BK)
      IF (MOD(INU,2).EQ.1) P1 = -P1
   20 CONTINUE
      DO 50 K=1,IL
        SQK = FDN - 1.0E0
        ATOL = S*ABS(SQK)
        SGN = 1.0E0
        CS1 = CONE
        CS2 = CONE
        CK = CONE
        AK = 0.0E0
        AA = 1.0E0
        BB = AEZ
        DK = EZ
        DO 30 J=1,JL
          CK = CK*CMPLX(SQK,0.0E0)/DK
          CS2 = CS2 + CK
          SGN = -SGN
          CS1 = CS1 + CK*CMPLX(SGN,0.0E0)
          DK = DK + EZ
          AA = AA*ABS(SQK)/BB
          BB = BB + AEZ
          AK = AK + 8.0E0
          SQK = SQK - AK
          IF (AA.LE.ATOL) GO TO 40
   30   CONTINUE
        GO TO 90
   40   CONTINUE
        S2 = CS1
        IF (X+X.LT.ELIM) S2 = S2 + P1*CS2*CEXP(-Z-Z)
        FDN = FDN + 8.0E0*DFNU + 4.0E0
        P1 = -P1
        M = N - IL + K
        Y(M) = S2*AK1
   50 CONTINUE
      IF (N.LE.2) RETURN
      NN = N
      K = NN - 2
      AK = FLOAT(K)
      RZ = (CONE+CONE)/Z
      IB = 3
      DO 60 I=IB,NN
        Y(K) = CMPLX(AK+FNU,0.0E0)*RZ*Y(K+1) + Y(K+2)
        AK = AK - 1.0E0
        K = K - 1
   60 CONTINUE
      IF (KODED.EQ.0) RETURN
      CK = CEXP(CZ)
      DO 70 I=1,NN
        Y(I) = Y(I)*CK
   70 CONTINUE
      RETURN
   80 CONTINUE
      NZ = -1
      RETURN
   90 CONTINUE
      NZ=-2
      RETURN
      END
      SUBROUTINE CBUNK(Z, FNU, KODE, MR, N, Y, NZ, TOL, ELIM, ALIM)
C***BEGIN PROLOGUE  CBUNK
C***REFER TO  CBESK,CBESH
C
C     CBUNK COMPUTES THE K BESSEL FUNCTION FOR FNU.GT.FNUL.
C     ACCORDING TO THE UNIFORM ASYMPTOTIC EXPANSION FOR K(FNU,Z)
C     IN CUNK1 AND THE EXPANSION FOR H(2,FNU,Z) IN CUNK2
C
C***ROUTINES CALLED  CUNK1,CUNK2
C***END PROLOGUE  CBUNK
      COMPLEX Y, Z
      REAL ALIM, AX, AY, ELIM, FNU, TOL, XX, YY
      INTEGER KODE, MR, N, NZ
      DIMENSION Y(N)
      NZ = 0
      XX = REAL(Z)
      YY = AIMAG(Z)
      AX = ABS(XX)*1.7321E0
      AY = ABS(YY)
      IF (AY.GT.AX) GO TO 10
C-----------------------------------------------------------------------
C     ASYMPTOTIC EXPANSION FOR K(FNU,Z) FOR LARGE FNU APPLIED IN
C     -PI/3.LE.ARG(Z).LE.PI/3
C-----------------------------------------------------------------------
      CALL CUNK1(Z, FNU, KODE, MR, N, Y, NZ, TOL, ELIM, ALIM)
      GO TO 20
   10 CONTINUE
C-----------------------------------------------------------------------
C     ASYMPTOTIC EXPANSION FOR H(2,FNU,Z*EXP(M*HPI)) FOR LARGE FNU
C     APPLIED IN PI/3.LT.ABS(ARG(Z)).LE.PI/2 WHERE M=+I OR -I
C     AND HPI=PI/2
C-----------------------------------------------------------------------
      CALL CUNK2(Z, FNU, KODE, MR, N, Y, NZ, TOL, ELIM, ALIM)
   20 CONTINUE
      RETURN
      END
      SUBROUTINE CUNK1(Z, FNU, KODE, MR, N, Y, NZ, TOL, ELIM, ALIM)
C***BEGIN PROLOGUE  CUNK1
C***REFER TO  CBESK
C
C     CUNK1 COMPUTES K(FNU,Z) AND ITS ANALYTIC CONTINUATION FROM THE
C     RIGHT HALF PLANE TO THE LEFT HALF PLANE BY MEANS OF THE
C     UNIFORM ASYMPTOTIC EXPANSION.
C     MR INDICATES THE DIRECTION OF ROTATION FOR ANALYTIC CONTINUATION.
C     NZ=-1 MEANS AN OVERFLOW WILL OCCUR
C
C***ROUTINES CALLED  CS1S2,CUCHK,CUNIK,R1MACH
C***END PROLOGUE  CUNK1
      COMPLEX CFN, CK, CONE, CRSC, CS, CSCL, CSGN, CSPN, CSR, CSS,
     * CWRK, CY, CZERO, C1, C2, PHI,  RZ, SUM,  S1, S2, Y, Z,
     * ZETA1,  ZETA2,  ZR, PHID, ZETA1D, ZETA2D, SUMD
      REAL ALIM, ANG, APHI, ASC, ASCLE, BRY, CPN, C2I, C2M, C2R, ELIM,
     * FMR, FN, FNF, FNU, PI, RS1, SGN, SPN, TOL, X, R1MACH
      INTEGER I, IB, IFLAG, IFN, IL, INIT, INU, IUF, K, KDFLG, KFLAG,
     * KK, KODE, MR, N, NW, NZ, J, IPARD, INITD, IC
      DIMENSION BRY(3), INIT(2), Y(N), SUM(2), PHI(2), ZETA1(2),
     * ZETA2(2), CY(2), CWRK(16,3), CSS(3), CSR(3)
      DATA CZERO, CONE / (0.0E0,0.0E0) , (1.0E0,0.0E0) /
      DATA PI / 3.14159265358979324E0 /
C
      KDFLG = 1
      NZ = 0
C-----------------------------------------------------------------------
C     EXP(-ALIM)=EXP(-ELIM)/TOL=APPROX. ONE PRECISION GREATER THAN
C     THE UNDERFLOW LIMIT
C-----------------------------------------------------------------------
      CSCL = CMPLX(1.0E0/TOL,0.0E0)
      CRSC = CMPLX(TOL,0.0E0)
      CSS(1) = CSCL
      CSS(2) = CONE
      CSS(3) = CRSC
      CSR(1) = CRSC
      CSR(2) = CONE
      CSR(3) = CSCL
      BRY(1) = 1.0E+3*R1MACH(1)/TOL
      BRY(2) = 1.0E0/BRY(1)
      BRY(3) = R1MACH(2)
      X = REAL(Z)
      ZR = Z
      IF (X.LT.0.0E0) ZR = -Z
      J=2
      DO 70 I=1,N
C-----------------------------------------------------------------------
C     J FLIP FLOPS BETWEEN 1 AND 2 IN J = 3 - J
C-----------------------------------------------------------------------
        J = 3 - J
        FN = FNU + FLOAT(I-1)
        INIT(J) = 0
        CALL CUNIK(ZR, FN, 2, 0, TOL, INIT(J), PHI(J), ZETA1(J),
     *   ZETA2(J), SUM(J), CWRK(1,J))
        IF (KODE.EQ.1) GO TO 20
        CFN = CMPLX(FN,0.0E0)
        S1 = ZETA1(J) - CFN*(CFN/(ZR+ZETA2(J)))
        GO TO 30
   20   CONTINUE
        S1 = ZETA1(J) - ZETA2(J)
   30   CONTINUE
C-----------------------------------------------------------------------
C     TEST FOR UNDERFLOW AND OVERFLOW
C-----------------------------------------------------------------------
        RS1 = REAL(S1)
        IF (ABS(RS1).GT.ELIM) GO TO 60
        IF (KDFLG.EQ.1) KFLAG = 2
        IF (ABS(RS1).LT.ALIM) GO TO 40
C-----------------------------------------------------------------------
C     REFINE  TEST AND SCALE
C-----------------------------------------------------------------------
        APHI = CABS(PHI(J))
        RS1 = RS1 + ALOG(APHI)
        IF (ABS(RS1).GT.ELIM) GO TO 60
        IF (KDFLG.EQ.1) KFLAG = 1
        IF (RS1.LT.0.0E0) GO TO 40
        IF (KDFLG.EQ.1) KFLAG = 3
   40   CONTINUE
C-----------------------------------------------------------------------
C     SCALE S1 TO KEEP INTERMEDIATE ARITHMETIC ON SCALE NEAR
C     EXPONENT EXTREMES
C-----------------------------------------------------------------------
        S2 = PHI(J)*SUM(J)
        C2R = REAL(S1)
        C2I = AIMAG(S1)
        C2M = EXP(C2R)*REAL(CSS(KFLAG))
        S1 = CMPLX(C2M,0.0E0)*CMPLX(COS(C2I),SIN(C2I))
        S2 = S2*S1
        IF (KFLAG.NE.1) GO TO 50
        CALL CUCHK(S2, NW, BRY(1), TOL)
        IF (NW.NE.0) GO TO 60
   50   CONTINUE
        CY(KDFLG) = S2
        Y(I) = S2*CSR(KFLAG)
        IF (KDFLG.EQ.2) GO TO 75
        KDFLG = 2
        GO TO 70
   60   CONTINUE
        IF (RS1.GT.0.0E0) GO TO 290
C-----------------------------------------------------------------------
C     FOR X.LT.0.0, THE I FUNCTION TO BE ADDED WILL OVERFLOW
C-----------------------------------------------------------------------
        IF (X.LT.0.0E0) GO TO 290
        KDFLG = 1
        Y(I) = CZERO
        NZ=NZ+1
        IF (I.EQ.1) GO TO 70
        IF (Y(I-1).EQ.CZERO) GO TO 70
        Y(I-1) = CZERO
        NZ=NZ+1
   70 CONTINUE
      I=N
   75 CONTINUE
      RZ = CMPLX(2.0E0,0.0E0)/ZR
      CK = CMPLX(FN,0.0E0)*RZ
      IB = I+1
      IF (N.LT.IB) GO TO 160
C-----------------------------------------------------------------------
C     TEST LAST MEMBER FOR UNDERFLOW AND OVERFLOW, SET SEQUENCE TO ZERO
C     ON UNDERFLOW
C-----------------------------------------------------------------------
      FN = FNU+FLOAT(N-1)
      IPARD = 1
      IF (MR.NE.0) IPARD = 0
      INITD = 0
      CALL CUNIK(ZR,FN,2,IPARD,TOL,INITD,PHID,ZETA1D,ZETA2D,SUMD,
     *CWRK(1,3))
      IF (KODE.EQ.1) GO TO 80
      CFN=CMPLX(FN,0.0E0)
      S1=ZETA1D-CFN*(CFN/(ZR+ZETA2D))
      GO TO 90
   80 CONTINUE
      S1=ZETA1D-ZETA2D
   90 CONTINUE
      RS1=REAL(S1)
      IF (ABS(RS1).GT.ELIM) GO TO 95
      IF (ABS(RS1).LT.ALIM) GO TO 100
C-----------------------------------------------------------------------
C     REFINE ESTIMATE AND TEST
C-----------------------------------------------------------------------
      APHI=CABS(PHID)
      RS1=RS1+ALOG(APHI)
      IF (ABS(RS1).LT.ELIM) GO TO 100
   95 CONTINUE
      IF (RS1.GT.0.0E0) GO TO 290
C-----------------------------------------------------------------------
C     FOR X.LT.0.0, THE I FUNCTION TO BE ADDED WILL OVERFLOW
C-----------------------------------------------------------------------
      IF (X.LT.0.0E0) GO TO 290
      NZ=N
      DO 96 I=1,N
        Y(I) = CZERO
   96 CONTINUE
      RETURN
  100 CONTINUE
C-----------------------------------------------------------------------
C     RECUR FORWARD FOR REMAINDER OF THE SEQUENCE
C-----------------------------------------------------------------------
      S1 = CY(1)
      S2 = CY(2)
      C1 = CSR(KFLAG)
      ASCLE = BRY(KFLAG)
      DO 120 I=IB,N
        C2 = S2
        S2 = CK*S2 + S1
        S1 = C2
        CK = CK + RZ
        C2 = S2*C1
        Y(I) = C2
        IF (KFLAG.GE.3) GO TO 120
        C2R = REAL(C2)
        C2I = AIMAG(C2)
        C2R = ABS(C2R)
        C2I = ABS(C2I)
        C2M = AMAX1(C2R,C2I)
        IF (C2M.LE.ASCLE) GO TO 120
        KFLAG = KFLAG + 1
        ASCLE = BRY(KFLAG)
        S1 = S1*C1
        S2 = C2
        S1 = S1*CSS(KFLAG)
        S2 = S2*CSS(KFLAG)
        C1 = CSR(KFLAG)
  120 CONTINUE
  160 CONTINUE
      IF (MR.EQ.0) RETURN
C-----------------------------------------------------------------------
C     ANALYTIC CONTINUATION FOR RE(Z).LT.0.0E0
C-----------------------------------------------------------------------
      NZ = 0
      FMR = FLOAT(MR)
      SGN = -SIGN(PI,FMR)
C-----------------------------------------------------------------------
C     CSPN AND CSGN ARE COEFF OF K AND I FUNCIONS RESP.
C-----------------------------------------------------------------------
      CSGN = CMPLX(0.0E0,SGN)
      INU = INT(FNU)
      FNF = FNU - FLOAT(INU)
      IFN = INU + N - 1
      ANG = FNF*SGN
      CPN = COS(ANG)
      SPN = SIN(ANG)
      CSPN = CMPLX(CPN,SPN)
      IF (MOD(IFN,2).EQ.1) CSPN = -CSPN
      ASC = BRY(1)
      KK = N
      IUF = 0
      KDFLG = 1
      IB = IB-1
      IC = IB-1
      DO 260 K=1,N
        FN = FNU + FLOAT(KK-1)
C-----------------------------------------------------------------------
C     LOGIC TO SORT OUT CASES WHOSE PARAMETERS WERE SET FOR THE K
C     FUNCTION ABOVE
C-----------------------------------------------------------------------
        M=3
        IF (N.GT.2) GO TO 175
  170   CONTINUE
        INITD = INIT(J)
        PHID = PHI(J)
        ZETA1D = ZETA1(J)
        ZETA2D = ZETA2(J)
        SUMD = SUM(J)
        M = J
        J = 3 - J
        GO TO 180
  175   CONTINUE
        IF ((KK.EQ.N).AND.(IB.LT.N)) GO TO 180
        IF ((KK.EQ.IB).OR.(KK.EQ.IC)) GO TO 170
        INITD = 0
  180   CONTINUE
        CALL CUNIK(ZR, FN, 1, 0, TOL, INITD, PHID, ZETA1D,
     *   ZETA2D, SUMD, CWRK(1,M))
        IF (KODE.EQ.1) GO TO 190
        CFN = CMPLX(FN,0.0E0)
        S1 = -ZETA1D + CFN*(CFN/(ZR+ZETA2D))
        GO TO 200
  190   CONTINUE
        S1 = -ZETA1D + ZETA2D
  200   CONTINUE
C-----------------------------------------------------------------------
C     TEST FOR UNDERFLOW AND OVERFLOW
C-----------------------------------------------------------------------
        RS1 = REAL(S1)
        IF (ABS(RS1).GT.ELIM) GO TO 250
        IF (KDFLG.EQ.1) IFLAG = 2
        IF (ABS(RS1).LT.ALIM) GO TO 210
C-----------------------------------------------------------------------
C     REFINE  TEST AND SCALE
C-----------------------------------------------------------------------
        APHI = CABS(PHID)
        RS1 = RS1 + ALOG(APHI)
        IF (ABS(RS1).GT.ELIM) GO TO 250
        IF (KDFLG.EQ.1) IFLAG = 1
        IF (RS1.LT.0.0E0) GO TO 210
        IF (KDFLG.EQ.1) IFLAG = 3
  210   CONTINUE
        S2 = CSGN*PHID*SUMD
        C2R = REAL(S1)
        C2I = AIMAG(S1)
        C2M = EXP(C2R)*REAL(CSS(IFLAG))
        S1 = CMPLX(C2M,0.0E0)*CMPLX(COS(C2I),SIN(C2I))
        S2 = S2*S1
        IF (IFLAG.NE.1) GO TO 220
        CALL CUCHK(S2, NW, BRY(1), TOL)
        IF (NW.NE.0) S2 = CMPLX(0.0E0,0.0E0)
  220   CONTINUE
        CY(KDFLG) = S2
        C2 = S2
        S2 = S2*CSR(IFLAG)
C-----------------------------------------------------------------------
C     ADD I AND K FUNCTIONS, K SEQUENCE IN Y(I), I=1,N
C-----------------------------------------------------------------------
        S1 = Y(KK)
        IF (KODE.EQ.1) GO TO 240
        CALL CS1S2(ZR, S1, S2, NW, ASC, ALIM, IUF)
        NZ = NZ + NW
  240   CONTINUE
        Y(KK) = S1*CSPN + S2
        KK = KK - 1
        CSPN = -CSPN
        IF (C2.NE.CZERO) GO TO 245
        KDFLG = 1
        GO TO 260
  245   CONTINUE
        IF (KDFLG.EQ.2) GO TO 265
        KDFLG = 2
        GO TO 260
  250   CONTINUE
        IF (RS1.GT.0.0E0) GO TO 290
        S2 = CZERO
        GO TO 220
  260 CONTINUE
      K = N
  265 CONTINUE
      IL = N - K
      IF (IL.EQ.0) RETURN
C-----------------------------------------------------------------------
C     RECUR BACKWARD FOR REMAINDER OF I SEQUENCE AND ADD IN THE
C     K FUNCTIONS, SCALING THE I SEQUENCE DURING RECURRENCE TO KEEP
C     INTERMEDIATE ARITHMETIC ON SCALE NEAR EXPONENT EXTREMES.
C-----------------------------------------------------------------------
      S1 = CY(1)
      S2 = CY(2)
      CS = CSR(IFLAG)
      ASCLE = BRY(IFLAG)
      FN = FLOAT(INU+IL)
      DO 280 I=1,IL
        C2 = S2
        S2 = S1 + CMPLX(FN+FNF,0.0E0)*RZ*S2
        S1 = C2
        FN = FN - 1.0E0
        C2 = S2*CS
        CK = C2
        C1 = Y(KK)
        IF (KODE.EQ.1) GO TO 270
        CALL CS1S2(ZR, C1, C2, NW, ASC, ALIM, IUF)
        NZ = NZ + NW
  270   CONTINUE
        Y(KK) = C1*CSPN + C2
        KK = KK - 1
        CSPN = -CSPN
        IF (IFLAG.GE.3) GO TO 280
        C2R = REAL(CK)
        C2I = AIMAG(CK)
        C2R = ABS(C2R)
        C2I = ABS(C2I)
        C2M = AMAX1(C2R,C2I)
        IF (C2M.LE.ASCLE) GO TO 280
        IFLAG = IFLAG + 1
        ASCLE = BRY(IFLAG)
        S1 = S1*CS
        S2 = CK
        S1 = S1*CSS(IFLAG)
        S2 = S2*CSS(IFLAG)
        CS = CSR(IFLAG)
  280 CONTINUE
      RETURN
  290 CONTINUE
      NZ = -1
      RETURN
      END
      SUBROUTINE CUNK2(Z, FNU, KODE, MR, N, Y, NZ, TOL, ELIM, ALIM)
C***BEGIN PROLOGUE  CUNK2
C***REFER TO  CBESK
C
C     CUNK2 COMPUTES K(FNU,Z) AND ITS ANALYTIC CONTINUATION FROM THE
C     RIGHT HALF PLANE TO THE LEFT HALF PLANE BY MEANS OF THE
C     UNIFORM ASYMPTOTIC EXPANSIONS FOR H(KIND,FNU,ZN) AND J(FNU,ZN)
C     WHERE ZN IS IN THE RIGHT HALF PLANE, KIND=(3-MR)/2, MR=+1 OR
C     -1. HERE ZN=ZR*I OR -ZR*I WHERE ZR=Z IF Z IS IN THE RIGHT
C     HALF PLANE OR ZR=-Z IF Z IS IN THE LEFT HALF PLANE. MR INDIC-
C     ATES THE DIRECTION OF ROTATION FOR ANALYTIC CONTINUATION.
C     NZ=-1 MEANS AN OVERFLOW WILL OCCUR
C
C***ROUTINES CALLED  CAIRY,CS1S2,CUCHK,CUNHJ,R1MACH
C***END PROLOGUE  CUNK2
      COMPLEX AI, ARG, ASUM, BSUM, CFN, CI, CIP,
     * CK, CONE, CRSC, CR1, CR2, CS, CSCL, CSGN, CSPN, CSR, CSS, CY,
     * CZERO, C1, C2, DAI, PHI,  RZ, S1, S2, Y, Z, ZB, ZETA1,
     * ZETA2, ZN, ZR, PHID, ARGD, ZETA1D, ZETA2D, ASUMD, BSUMD
      REAL AARG, AIC, ALIM, ANG, APHI, ASC, ASCLE, BRY, CAR, CPN, C2I,
     * C2M, C2R, ELIM, FMR, FN, FNF, FNU, HPI, PI, RS1, SAR, SGN, SPN,
     * TOL, X, YY, R1MACH
      INTEGER I, IB, IFLAG, IFN, IL, IN, INU, IUF, K, KDFLG, KFLAG, KK,
     * KODE, MR, N, NAI, NDAI, NW, NZ, IDUM, J, IPARD, IC
      DIMENSION BRY(3), Y(N), ASUM(2), BSUM(2), PHI(2), ARG(2),
     * ZETA1(2), ZETA2(2), CY(2), CIP(4), CSS(3), CSR(3)
      DATA CZERO, CONE, CI, CR1, CR2 /
     1         (0.0E0,0.0E0),(1.0E0,0.0E0),(0.0E0,1.0E0),
     1(1.0E0,1.73205080756887729E0),(-0.5E0,-8.66025403784438647E-01)/
      DATA HPI, PI, AIC /
     1     1.57079632679489662E+00,     3.14159265358979324E+00,
     1     1.26551212348464539E+00/
      DATA CIP(1),CIP(2),CIP(3),CIP(4)/
     1 (1.0E0,0.0E0), (0.0E0,-1.0E0), (-1.0E0,0.0E0), (0.0E0,1.0E0)/
C
      KDFLG = 1
      NZ = 0
C-----------------------------------------------------------------------
C     EXP(-ALIM)=EXP(-ELIM)/TOL=APPROX. ONE PRECISION GREATER THAN
C     THE UNDERFLOW LIMIT
C-----------------------------------------------------------------------
      CSCL = CMPLX(1.0E0/TOL,0.0E0)
      CRSC = CMPLX(TOL,0.0E0)
      CSS(1) = CSCL
      CSS(2) = CONE
      CSS(3) = CRSC
      CSR(1) = CRSC
      CSR(2) = CONE
      CSR(3) = CSCL
      BRY(1) = 1.0E+3*R1MACH(1)/TOL
      BRY(2) = 1.0E0/BRY(1)
      BRY(3) = R1MACH(2)
      X = REAL(Z)
      ZR = Z
      IF (X.LT.0.0E0) ZR = -Z
      YY = AIMAG(ZR)
      ZN = -ZR*CI
      ZB = ZR
      INU = INT(FNU)
      FNF = FNU - FLOAT(INU)
      ANG = -HPI*FNF
      CAR = COS(ANG)
      SAR = SIN(ANG)
      CPN = -HPI*CAR
      SPN = -HPI*SAR
      C2 = CMPLX(-SPN,CPN)
      KK = MOD(INU,4) + 1
      CS = CR1*C2*CIP(KK)
      IF (YY.GT.0.0E0) GO TO 10
      ZN = CONJG(-ZN)
      ZB = CONJG(ZB)
   10 CONTINUE
C-----------------------------------------------------------------------
C     K(FNU,Z) IS COMPUTED FROM H(2,FNU,-I*Z) WHERE Z IS IN THE FIRST
C     QUADRANT. FOURTH QUADRANT VALUES (YY.LE.0.0E0) ARE COMPUTED BY
C     CONJUGATION SINCE THE K FUNCTION IS REAL ON THE POSITIVE REAL AXIS
C-----------------------------------------------------------------------
      J = 2
      DO 70 I=1,N
C-----------------------------------------------------------------------
C     J FLIP FLOPS BETWEEN 1 AND 2 IN J = 3 - J
C-----------------------------------------------------------------------
        J = 3 - J
        FN = FNU + FLOAT(I-1)
        CALL CUNHJ(ZN, FN, 0, TOL, PHI(J), ARG(J), ZETA1(J), ZETA2(J),
     *   ASUM(J), BSUM(J))
        IF (KODE.EQ.1) GO TO 20
        CFN = CMPLX(FN,0.0E0)
        S1 = ZETA1(J) - CFN*(CFN/(ZB+ZETA2(J)))
        GO TO 30
   20   CONTINUE
        S1 = ZETA1(J) - ZETA2(J)
   30   CONTINUE
C-----------------------------------------------------------------------
C     TEST FOR UNDERFLOW AND OVERFLOW
C-----------------------------------------------------------------------
        RS1 = REAL(S1)
        IF (ABS(RS1).GT.ELIM) GO TO 60
        IF (KDFLG.EQ.1) KFLAG = 2
        IF (ABS(RS1).LT.ALIM) GO TO 40
C-----------------------------------------------------------------------
C     REFINE  TEST AND SCALE
C-----------------------------------------------------------------------
        APHI = CABS(PHI(J))
        AARG = CABS(ARG(J))
        RS1 = RS1 + ALOG(APHI) - 0.25E0*ALOG(AARG) - AIC
        IF (ABS(RS1).GT.ELIM) GO TO 60
        IF (KDFLG.EQ.1) KFLAG = 1
        IF (RS1.LT.0.0E0) GO TO 40
        IF (KDFLG.EQ.1) KFLAG = 3
   40   CONTINUE
C-----------------------------------------------------------------------
C     SCALE S1 TO KEEP INTERMEDIATE ARITHMETIC ON SCALE NEAR
C     EXPONENT EXTREMES
C-----------------------------------------------------------------------
        C2 = ARG(J)*CR2
        CALL CAIRY(C2, 0, 2, AI, NAI, IDUM)
        CALL CAIRY(C2, 1, 2, DAI, NDAI, IDUM)
        S2 = CS*PHI(J)*(AI*ASUM(J)+CR2*DAI*BSUM(J))
        C2R = REAL(S1)
        C2I = AIMAG(S1)
        C2M = EXP(C2R)*REAL(CSS(KFLAG))
        S1 = CMPLX(C2M,0.0E0)*CMPLX(COS(C2I),SIN(C2I))
        S2 = S2*S1
        IF (KFLAG.NE.1) GO TO 50
        CALL CUCHK(S2, NW, BRY(1), TOL)
        IF (NW.NE.0) GO TO 60
   50   CONTINUE
        IF (YY.LE.0.0E0) S2 = CONJG(S2)
        CY(KDFLG) = S2
        Y(I) = S2*CSR(KFLAG)
        CS = -CI*CS
        IF (KDFLG.EQ.2) GO TO 75
        KDFLG = 2
        GO TO 70
   60   CONTINUE
        IF (RS1.GT.0.0E0) GO TO 300
C-----------------------------------------------------------------------
C     FOR X.LT.0.0, THE I FUNCTION TO BE ADDED WILL OVERFLOW
C-----------------------------------------------------------------------
        IF (X.LT.0.0E0) GO TO 300
        KDFLG = 1
        Y(I) = CZERO
        CS = -CI*CS
        NZ=NZ+1
        IF (I.EQ.1) GO TO 70
        IF (Y(I-1).EQ.CZERO) GO TO 70
        Y(I-1) = CZERO
        NZ=NZ+1
   70 CONTINUE
      I=N
   75 CONTINUE
      RZ = CMPLX(2.0E0,0.0E0)/ZR
      CK = CMPLX(FN,0.0E0)*RZ
      IB = I + 1
      IF (N.LT.IB) GO TO 170
C-----------------------------------------------------------------------
C     TEST LAST MEMBER FOR UNDERFLOW AND OVERFLOW, SET SEQUENCE TO ZERO
C     ON UNDERFLOW
C-----------------------------------------------------------------------
      FN = FNU+FLOAT(N-1)
      IPARD = 1
      IF (MR.NE.0) IPARD = 0
      CALL CUNHJ(ZN,FN,IPARD,TOL,PHID,ARGD,ZETA1D,ZETA2D,ASUMD,BSUMD)
      IF (KODE.EQ.1) GO TO 80
      CFN=CMPLX(FN,0.0E0)
      S1=ZETA1D-CFN*(CFN/(ZB+ZETA2D))
      GO TO 90
   80 CONTINUE
      S1=ZETA1D-ZETA2D
   90 CONTINUE
      RS1=REAL(S1)
      IF (ABS(RS1).GT.ELIM) GO TO 95
      IF (ABS(RS1).LT.ALIM) GO TO 100
C-----------------------------------------------------------------------
C     REFINE ESTIMATE AND TEST
C-----------------------------------------------------------------------
      APHI=CABS(PHID)
      AARG = CABS(ARGD)
      RS1=RS1+ALOG(APHI)-0.25E0*ALOG(AARG)-AIC
      IF (ABS(RS1).LT.ELIM) GO TO 100
   95 CONTINUE
      IF (RS1.GT.0.0E0) GO TO 300
C-----------------------------------------------------------------------
C     FOR X.LT.0.0, THE I FUNCTION TO BE ADDED WILL OVERFLOW
C-----------------------------------------------------------------------
      IF (X.LT.0.0E0) GO TO 300
      NZ=N
      DO 96 I=1,N
        Y(I) = CZERO
   96 CONTINUE
      RETURN
  100 CONTINUE
C-----------------------------------------------------------------------
C     SCALED FORWARD RECURRENCE FOR REMAINDER OF THE SEQUENCE
C-----------------------------------------------------------------------
      S1 = CY(1)
      S2 = CY(2)
      C1 = CSR(KFLAG)
      ASCLE = BRY(KFLAG)
      DO 120 I=IB,N
        C2 = S2
        S2 = CK*S2 + S1
        S1 = C2
        CK = CK + RZ
        C2 = S2*C1
        Y(I) = C2
        IF (KFLAG.GE.3) GO TO 120
        C2R = REAL(C2)
        C2I = AIMAG(C2)
        C2R = ABS(C2R)
        C2I = ABS(C2I)
        C2M = AMAX1(C2R,C2I)
        IF (C2M.LE.ASCLE) GO TO 120
        KFLAG = KFLAG + 1
        ASCLE = BRY(KFLAG)
        S1 = S1*C1
        S2 = C2
        S1 = S1*CSS(KFLAG)
        S2 = S2*CSS(KFLAG)
        C1 = CSR(KFLAG)
  120 CONTINUE
  170 CONTINUE
      IF (MR.EQ.0) RETURN
C-----------------------------------------------------------------------
C     ANALYTIC CONTINUATION FOR RE(Z).LT.0.0E0
C-----------------------------------------------------------------------
      NZ = 0
      FMR = FLOAT(MR)
      SGN = -SIGN(PI,FMR)
C-----------------------------------------------------------------------
C     CSPN AND CSGN ARE COEFF OF K AND I FUNCTIONS RESP.
C-----------------------------------------------------------------------
      CSGN = CMPLX(0.0E0,SGN)
      IF (YY.LE.0.0E0) CSGN = CONJG(CSGN)
      IFN = INU + N - 1
      ANG = FNF*SGN
      CPN = COS(ANG)
      SPN = SIN(ANG)
      CSPN = CMPLX(CPN,SPN)
      IF (MOD(IFN,2).EQ.1) CSPN = -CSPN
C-----------------------------------------------------------------------
C     CS=COEFF OF THE J FUNCTION TO GET THE I FUNCTION. I(FNU,Z) IS
C     COMPUTED FROM EXP(I*FNU*HPI)*J(FNU,-I*Z) WHERE Z IS IN THE FIRST
C     QUADRANT. FOURTH QUADRANT VALUES (YY.LE.0.0E0) ARE COMPUTED BY
C     CONJUGATION SINCE THE I FUNCTION IS REAL ON THE POSITIVE REAL AXIS
C-----------------------------------------------------------------------
      CS = CMPLX(CAR,-SAR)*CSGN
      IN = MOD(IFN,4) + 1
      C2 = CIP(IN)
      CS = CS*CONJG(C2)
      ASC = BRY(1)
      KK = N
      KDFLG = 1
      IB = IB-1
      IC = IB-1
      IUF = 0
      DO 270 K=1,N
C-----------------------------------------------------------------------
C     LOGIC TO SORT OUT CASES WHOSE PARAMETERS WERE SET FOR THE K
C     FUNCTION ABOVE
C-----------------------------------------------------------------------
        FN = FNU+FLOAT(KK-1)
        IF (N.GT.2) GO TO 180
  175   CONTINUE
        PHID = PHI(J)
        ARGD = ARG(J)
        ZETA1D = ZETA1(J)
        ZETA2D = ZETA2(J)
        ASUMD = ASUM(J)
        BSUMD = BSUM(J)
        J = 3 - J
        GO TO 190
  180   CONTINUE
        IF ((KK.EQ.N).AND.(IB.LT.N)) GO TO 190
        IF ((KK.EQ.IB).OR.(KK.EQ.IC)) GO TO 175
        CALL CUNHJ(ZN, FN, 0, TOL, PHID, ARGD, ZETA1D, ZETA2D,
     *   ASUMD, BSUMD)
  190   CONTINUE
        IF (KODE.EQ.1) GO TO 200
        CFN = CMPLX(FN,0.0E0)
        S1 = -ZETA1D + CFN*(CFN/(ZB+ZETA2D))
        GO TO 210
  200   CONTINUE
        S1 = -ZETA1D + ZETA2D
  210   CONTINUE
C-----------------------------------------------------------------------
C     TEST FOR UNDERFLOW AND OVERFLOW
C-----------------------------------------------------------------------
        RS1 = REAL(S1)
        IF (ABS(RS1).GT.ELIM) GO TO 260
        IF (KDFLG.EQ.1) IFLAG = 2
        IF (ABS(RS1).LT.ALIM) GO TO 220
C-----------------------------------------------------------------------
C     REFINE  TEST AND SCALE
C-----------------------------------------------------------------------
        APHI = CABS(PHID)
        AARG = CABS(ARGD)
        RS1 = RS1 + ALOG(APHI) - 0.25E0*ALOG(AARG) - AIC
        IF (ABS(RS1).GT.ELIM) GO TO 260
        IF (KDFLG.EQ.1) IFLAG = 1
        IF (RS1.LT.0.0E0) GO TO 220
        IF (KDFLG.EQ.1) IFLAG = 3
  220   CONTINUE
        CALL CAIRY(ARGD, 0, 2, AI, NAI, IDUM)
        CALL CAIRY(ARGD, 1, 2, DAI, NDAI, IDUM)
        S2 = CS*PHID*(AI*ASUMD+DAI*BSUMD)
        C2R = REAL(S1)
        C2I = AIMAG(S1)
        C2M = EXP(C2R)*REAL(CSS(IFLAG))
        S1 = CMPLX(C2M,0.0E0)*CMPLX(COS(C2I),SIN(C2I))
        S2 = S2*S1
        IF (IFLAG.NE.1) GO TO 230
        CALL CUCHK(S2, NW, BRY(1), TOL)
        IF (NW.NE.0) S2 = CMPLX(0.0E0,0.0E0)
  230   CONTINUE
        IF (YY.LE.0.0E0) S2 = CONJG(S2)
        CY(KDFLG) = S2
        C2 = S2
        S2 = S2*CSR(IFLAG)
C-----------------------------------------------------------------------
C     ADD I AND K FUNCTIONS, K SEQUENCE IN Y(I), I=1,N
C-----------------------------------------------------------------------
        S1 = Y(KK)
        IF (KODE.EQ.1) GO TO 250
        CALL CS1S2(ZR, S1, S2, NW, ASC, ALIM, IUF)
        NZ = NZ + NW
  250   CONTINUE
        Y(KK) = S1*CSPN + S2
        KK = KK - 1
        CSPN = -CSPN
        CS = -CS*CI
        IF (C2.NE.CZERO) GO TO 255
        KDFLG = 1
        GO TO 270
  255   CONTINUE
        IF (KDFLG.EQ.2) GO TO 275
        KDFLG = 2
        GO TO 270
  260   CONTINUE
        IF (RS1.GT.0.0E0) GO TO 300
        S2 = CZERO
        GO TO 230
  270 CONTINUE
      K = N
  275 CONTINUE
      IL = N-K
      IF (IL.EQ.0) RETURN
C-----------------------------------------------------------------------
C     RECUR BACKWARD FOR REMAINDER OF I SEQUENCE AND ADD IN THE
C     K FUNCTIONS, SCALING THE I SEQUENCE DURING RECURRENCE TO KEEP
C     INTERMEDIATE ARITHMETIC ON SCALE NEAR EXPONENT EXTREMES.
C-----------------------------------------------------------------------
      S1 = CY(1)
      S2 = CY(2)
      CS = CSR(IFLAG)
      ASCLE = BRY(IFLAG)
      FN = FLOAT(INU+IL)
      DO 290 I=1,IL
        C2 = S2
        S2 = S1 + CMPLX(FN+FNF,0.0E0)*RZ*S2
        S1 = C2
        FN = FN - 1.0E0
        C2 = S2*CS
        CK = C2
        C1 = Y(KK)
        IF (KODE.EQ.1) GO TO 280
        CALL CS1S2(ZR, C1, C2, NW, ASC, ALIM, IUF)
        NZ = NZ + NW
  280   CONTINUE
        Y(KK) = C1*CSPN + C2
        KK = KK - 1
        CSPN = -CSPN
        IF (IFLAG.GE.3) GO TO 290
        C2R = REAL(CK)
        C2I = AIMAG(CK)
        C2R = ABS(C2R)
        C2I = ABS(C2I)
        C2M = AMAX1(C2R,C2I)
        IF (C2M.LE.ASCLE) GO TO 290
        IFLAG = IFLAG + 1
        ASCLE = BRY(IFLAG)
        S1 = S1*CS
        S2 = CK
        S1 = S1*CSS(IFLAG)
        S2 = S2*CSS(IFLAG)
        CS = CSR(IFLAG)
  290 CONTINUE
      RETURN
  300 CONTINUE
      NZ = -1
      RETURN
      END
      SUBROUTINE CBUNI(Z, FNU, KODE, N, Y, NZ, NUI, NLAST, FNUL, TOL,
     * ELIM, ALIM)
C***BEGIN PROLOGUE  CBUNI
C***REFER TO  CBESI,CBESK
C
C     CBUNI COMPUTES THE I BESSEL FUNCTION FOR LARGE CABS(Z).GT.
C     FNUL AND FNU+N-1.LT.FNUL. THE ORDER IS INCREASED FROM
C     FNU+N-1 GREATER THAN FNUL BY ADDING NUI AND COMPUTING
C     ACCORDING TO THE UNIFORM ASYMPTOTIC EXPANSION FOR I(FNU,Z)
C     ON IFORM=1 AND THE EXPANSION FOR J(FNU,Z) ON IFORM=2
C
C***ROUTINES CALLED  CUNI1,CUNI2,R1MACH
C***END PROLOGUE  CBUNI
      COMPLEX CSCL, CSCR, CY, RZ, ST, S1, S2, Y, Z
      REAL ALIM, AX, AY, DFNU, ELIM, FNU, FNUI, FNUL, GNU, TOL, XX, YY,
     * ASCLE, BRY, STR, STI, STM, R1MACH
      INTEGER I, IFLAG, IFORM, K, KODE, N, NL, NLAST, NUI, NW, NZ
      DIMENSION Y(N), CY(2), BRY(3)
      NZ = 0
      XX = REAL(Z)
      YY = AIMAG(Z)
      AX = ABS(XX)*1.7321E0
      AY = ABS(YY)
      IFORM = 1
      IF (AY.GT.AX) IFORM = 2
      IF (NUI.EQ.0) GO TO 60
      FNUI = FLOAT(NUI)
      DFNU = FNU + FLOAT(N-1)
      GNU = DFNU + FNUI
      IF (IFORM.EQ.2) GO TO 10
C-----------------------------------------------------------------------
C     ASYMPTOTIC EXPANSION FOR I(FNU,Z) FOR LARGE FNU APPLIED IN
C     -PI/3.LE.ARG(Z).LE.PI/3
C-----------------------------------------------------------------------
      CALL CUNI1(Z, GNU, KODE, 2, CY, NW, NLAST, FNUL, TOL, ELIM, ALIM)
      GO TO 20
   10 CONTINUE
C-----------------------------------------------------------------------
C     ASYMPTOTIC EXPANSION FOR J(FNU,Z*EXP(M*HPI)) FOR LARGE FNU
C     APPLIED IN PI/3.LT.ABS(ARG(Z)).LE.PI/2 WHERE M=+I OR -I
C     AND HPI=PI/2
C-----------------------------------------------------------------------
      CALL CUNI2(Z, GNU, KODE, 2, CY, NW, NLAST, FNUL, TOL, ELIM, ALIM)
   20 CONTINUE
      IF (NW.LT.0) GO TO 50
      IF (NW.NE.0) GO TO 90
      AY = CABS(CY(1))
C----------------------------------------------------------------------
C     SCALE BACKWARD RECURRENCE, BRY(3) IS DEFINED BUT NEVER USED
C----------------------------------------------------------------------
      BRY(1) = 1.0E+3*R1MACH(1)/TOL
      BRY(2) = 1.0E0/BRY(1)
      BRY(3) = BRY(2)
      IFLAG = 2
      ASCLE = BRY(2)
      AX = 1.0E0
      CSCL = CMPLX(AX,0.0E0)
      IF (AY.GT.BRY(1)) GO TO 21
      IFLAG = 1
      ASCLE = BRY(1)
      AX = 1.0E0/TOL
      CSCL = CMPLX(AX,0.0E0)
      GO TO 25
   21 CONTINUE
      IF (AY.LT.BRY(2)) GO TO 25
      IFLAG = 3
      ASCLE = BRY(3)
      AX = TOL
      CSCL = CMPLX(AX,0.0E0)
   25 CONTINUE
      AY = 1.0E0/AX
      CSCR = CMPLX(AY,0.0E0)
      S1 = CY(2)*CSCL
      S2 = CY(1)*CSCL
      RZ = CMPLX(2.0E0,0.0E0)/Z
      DO 30 I=1,NUI
        ST = S2
        S2 = CMPLX(DFNU+FNUI,0.0E0)*RZ*S2 + S1
        S1 = ST
        FNUI = FNUI - 1.0E0
        IF (IFLAG.GE.3) GO TO 30
        ST = S2*CSCR
        STR = REAL(ST)
        STI = AIMAG(ST)
        STR = ABS(STR)
        STI = ABS(STI)
        STM = AMAX1(STR,STI)
        IF (STM.LE.ASCLE) GO TO 30
        IFLAG = IFLAG+1
        ASCLE = BRY(IFLAG)
        S1 = S1*CSCR
        S2 = ST
        AX = AX*TOL
        AY = 1.0E0/AX
        CSCL = CMPLX(AX,0.0E0)
        CSCR = CMPLX(AY,0.0E0)
        S1 = S1*CSCL
        S2 = S2*CSCL
   30 CONTINUE
      Y(N) = S2*CSCR
      IF (N.EQ.1) RETURN
      NL = N - 1
      FNUI = FLOAT(NL)
      K = NL
      DO 40 I=1,NL
        ST = S2
        S2 = CMPLX(FNU+FNUI,0.0E0)*RZ*S2 + S1
        S1 = ST
        ST = S2*CSCR
        Y(K) = ST
        FNUI = FNUI - 1.0E0
        K = K - 1
        IF (IFLAG.GE.3) GO TO 40
        STR = REAL(ST)
        STI = AIMAG(ST)
        STR = ABS(STR)
        STI = ABS(STI)
        STM = AMAX1(STR,STI)
        IF (STM.LE.ASCLE) GO TO 40
        IFLAG = IFLAG+1
        ASCLE = BRY(IFLAG)
        S1 = S1*CSCR
        S2 = ST
        AX = AX*TOL
        AY = 1.0E0/AX
        CSCL = CMPLX(AX,0.0E0)
        CSCR = CMPLX(AY,0.0E0)
        S1 = S1*CSCL
        S2 = S2*CSCL
   40 CONTINUE
      RETURN
   50 CONTINUE
      NZ = -1
      IF(NW.EQ.(-2)) NZ=-2
      RETURN
   60 CONTINUE
      IF (IFORM.EQ.2) GO TO 70
C-----------------------------------------------------------------------
C     ASYMPTOTIC EXPANSION FOR I(FNU,Z) FOR LARGE FNU APPLIED IN
C     -PI/3.LE.ARG(Z).LE.PI/3
C-----------------------------------------------------------------------
      CALL CUNI1(Z, FNU, KODE, N, Y, NW, NLAST, FNUL, TOL, ELIM, ALIM)
      GO TO 80
   70 CONTINUE
C-----------------------------------------------------------------------
C     ASYMPTOTIC EXPANSION FOR J(FNU,Z*EXP(M*HPI)) FOR LARGE FNU
C     APPLIED IN PI/3.LT.ABS(ARG(Z)).LE.PI/2 WHERE M=+I OR -I
C     AND HPI=PI/2
C-----------------------------------------------------------------------
      CALL CUNI2(Z, FNU, KODE, N, Y, NW, NLAST, FNUL, TOL, ELIM, ALIM)
   80 CONTINUE
      IF (NW.LT.0) GO TO 50
      NZ = NW
      RETURN
   90 CONTINUE
      NLAST = N
      RETURN
      END
      SUBROUTINE CUNI1(Z, FNU, KODE, N, Y, NZ, NLAST, FNUL, TOL, ELIM,
     * ALIM)
C***BEGIN PROLOGUE  CUNI1
C***REFER TO  CBESI,CBESK
C
C     CUNI1 COMPUTES I(FNU,Z)  BY MEANS OF THE UNIFORM ASYMPTOTIC
C     EXPANSION FOR I(FNU,Z) IN -PI/3.LE.ARG Z.LE.PI/3.
C
C     FNUL IS THE SMALLEST ORDER PERMITTED FOR THE ASYMPTOTIC
C     EXPANSION. NLAST=0 MEANS ALL OF THE Y VALUES WERE SET.
C     NLAST.NE.0 IS THE NUMBER LEFT TO BE COMPUTED BY ANOTHER
C     FORMULA FOR ORDERS FNU TO FNU+NLAST-1 BECAUSE FNU+NLAST-1.LT.FNUL.
C     Y(I)=CZERO FOR I=NLAST+1,N
C
C***ROUTINES CALLED  CUCHK,CUNIK,CUOIK,R1MACH
C***END PROLOGUE  CUNI1
      COMPLEX CFN, CONE, CRSC, CSCL, CSR, CSS, CWRK, CZERO, C1, C2,
     * PHI, RZ, SUM, S1, S2, Y, Z, ZETA1, ZETA2, CY
      REAL ALIM, APHI, ASCLE, BRY, C2I, C2M, C2R, ELIM, FN, FNU, FNUL,
     * RS1, TOL, YY, R1MACH
      INTEGER I, IFLAG, INIT, K, KODE, M, N, ND, NLAST, NN, NUF, NW, NZ
      DIMENSION BRY(3), Y(N), CWRK(16), CSS(3), CSR(3), CY(2)
      DATA CZERO, CONE / (0.0E0,0.0E0), (1.0E0,0.0E0) /
C
      NZ = 0
      ND = N
      NLAST = 0
C-----------------------------------------------------------------------
C     COMPUTED VALUES WITH EXPONENTS BETWEEN ALIM AND ELIM IN MAG-
C     NITUDE ARE SCALED TO KEEP INTERMEDIATE ARITHMETIC ON SCALE,
C     EXP(ALIM)=EXP(ELIM)*TOL
C-----------------------------------------------------------------------
      CSCL = CMPLX(1.0E0/TOL,0.0E0)
      CRSC = CMPLX(TOL,0.0E0)
      CSS(1) = CSCL
      CSS(2) = CONE
      CSS(3) = CRSC
      CSR(1) = CRSC
      CSR(2) = CONE
      CSR(3) = CSCL
      BRY(1) = 1.0E+3*R1MACH(1)/TOL
C-----------------------------------------------------------------------
C     CHECK FOR UNDERFLOW AND OVERFLOW ON FIRST MEMBER
C-----------------------------------------------------------------------
      FN = AMAX1(FNU,1.0E0)
      INIT = 0
      CALL CUNIK(Z, FN, 1, 1, TOL, INIT, PHI, ZETA1, ZETA2, SUM, CWRK)
      IF (KODE.EQ.1) GO TO 10
      CFN = CMPLX(FN,0.0E0)
      S1 = -ZETA1 + CFN*(CFN/(Z+ZETA2))
      GO TO 20
   10 CONTINUE
      S1 = -ZETA1 + ZETA2
   20 CONTINUE
      RS1 = REAL(S1)
      IF (ABS(RS1).GT.ELIM) GO TO 130
   30 CONTINUE
      NN = MIN0(2,ND)
      DO 80 I=1,NN
        FN = FNU + FLOAT(ND-I)
        INIT = 0
        CALL CUNIK(Z, FN, 1, 0, TOL, INIT, PHI, ZETA1, ZETA2, SUM, CWRK)
        IF (KODE.EQ.1) GO TO 40
        CFN = CMPLX(FN,0.0E0)
        YY = AIMAG(Z)
        S1 = -ZETA1 + CFN*(CFN/(Z+ZETA2)) + CMPLX(0.0E0,YY)
        GO TO 50
   40   CONTINUE
        S1 = -ZETA1 + ZETA2
   50   CONTINUE
C-----------------------------------------------------------------------
C     TEST FOR UNDERFLOW AND OVERFLOW
C-----------------------------------------------------------------------
        RS1 = REAL(S1)
        IF (ABS(RS1).GT.ELIM) GO TO 110
        IF (I.EQ.1) IFLAG = 2
        IF (ABS(RS1).LT.ALIM) GO TO 60
C-----------------------------------------------------------------------
C     REFINE  TEST AND SCALE
C-----------------------------------------------------------------------
        APHI = CABS(PHI)
        RS1 = RS1 + ALOG(APHI)
        IF (ABS(RS1).GT.ELIM) GO TO 110
        IF (I.EQ.1) IFLAG = 1
        IF (RS1.LT.0.0E0) GO TO 60
        IF (I.EQ.1) IFLAG = 3
   60   CONTINUE
C-----------------------------------------------------------------------
C     SCALE S1 IF CABS(S1).LT.ASCLE
C-----------------------------------------------------------------------
        S2 = PHI*SUM
        C2R = REAL(S1)
        C2I = AIMAG(S1)
        C2M = EXP(C2R)*REAL(CSS(IFLAG))
        S1 = CMPLX(C2M,0.0E0)*CMPLX(COS(C2I),SIN(C2I))
        S2 = S2*S1
        IF (IFLAG.NE.1) GO TO 70
        CALL CUCHK(S2, NW, BRY(1), TOL)
        IF (NW.NE.0) GO TO 110
   70   CONTINUE
        M = ND - I + 1
        CY(I) = S2
        Y(M) = S2*CSR(IFLAG)
   80 CONTINUE
      IF (ND.LE.2) GO TO 100
      RZ = CMPLX(2.0E0,0.0E0)/Z
      BRY(2) = 1.0E0/BRY(1)
      BRY(3) = R1MACH(2)
      S1 = CY(1)
      S2 = CY(2)
      C1 = CSR(IFLAG)
      ASCLE = BRY(IFLAG)
      K = ND - 2
      FN = FLOAT(K)
      DO 90 I=3,ND
        C2 = S2
        S2 = S1 + CMPLX(FNU+FN,0.0E0)*RZ*S2
        S1 = C2
        C2 = S2*C1
        Y(K) = C2
        K = K - 1
        FN = FN - 1.0E0
        IF (IFLAG.GE.3) GO TO 90
        C2R = REAL(C2)
        C2I = AIMAG(C2)
        C2R = ABS(C2R)
        C2I = ABS(C2I)
        C2M = AMAX1(C2R,C2I)
        IF (C2M.LE.ASCLE) GO TO 90
        IFLAG = IFLAG + 1
        ASCLE = BRY(IFLAG)
        S1 = S1*C1
        S2 = C2
        S1 = S1*CSS(IFLAG)
        S2 = S2*CSS(IFLAG)
        C1 = CSR(IFLAG)
   90 CONTINUE
  100 CONTINUE
      RETURN
C-----------------------------------------------------------------------
C     SET UNDERFLOW AND UPDATE PARAMETERS
C-----------------------------------------------------------------------
  110 CONTINUE
      IF (RS1.GT.0.0E0) GO TO 120
      Y(ND) = CZERO
      NZ = NZ + 1
      ND = ND - 1
      IF (ND.EQ.0) GO TO 100
      CALL CUOIK(Z, FNU, KODE, 1, ND, Y, NUF, TOL, ELIM, ALIM)
      IF (NUF.LT.0) GO TO 120
      ND = ND - NUF
      NZ = NZ + NUF
      IF (ND.EQ.0) GO TO 100
      FN = FNU + FLOAT(ND-1)
      IF (FN.GE.FNUL) GO TO 30
      NLAST = ND
      RETURN
  120 CONTINUE
      NZ = -1
      RETURN
  130 CONTINUE
      IF (RS1.GT.0.0E0) GO TO 120
      NZ = N
      DO 140 I=1,N
        Y(I) = CZERO
  140 CONTINUE
      RETURN
      END
      SUBROUTINE CUNI2(Z, FNU, KODE, N, Y, NZ, NLAST, FNUL, TOL, ELIM,
     * ALIM)
C***BEGIN PROLOGUE  CUNI2
C***REFER TO  CBESI,CBESK
C
C     CUNI2 COMPUTES I(FNU,Z) IN THE RIGHT HALF PLANE BY MEANS OF
C     UNIFORM ASYMPTOTIC EXPANSION FOR J(FNU,ZN) WHERE ZN IS Z*I
C     OR -Z*I AND ZN IS IN THE RIGHT HALF PLANE ALSO.
C
C     FNUL IS THE SMALLEST ORDER PERMITTED FOR THE ASYMPTOTIC
C     EXPANSION. NLAST=0 MEANS ALL OF THE Y VALUES WERE SET.
C     NLAST.NE.0 IS THE NUMBER LEFT TO BE COMPUTED BY ANOTHER
C     FORMULA FOR ORDERS FNU TO FNU+NLAST-1 BECAUSE FNU+NLAST-1.LT.FNUL.
C     Y(I)=CZERO FOR I=NLAST+1,N
C
C***ROUTINES CALLED  CAIRY,CUCHK,CUNHJ,CUOIK,R1MACH
C***END PROLOGUE  CUNI2
      COMPLEX AI, ARG, ASUM, BSUM, CFN, CI, CID, CIP, CONE, CRSC, CSCL,
     * CSR, CSS, CY, CZERO, C1, C2, DAI, PHI, RZ, S1, S2, Y, Z, ZB,
     * ZETA1, ZETA2, ZN, ZAR
      REAL AARG, AIC, ALIM, ANG, APHI, ASCLE, AY, BRY, CAR, C2I, C2M,
     * C2R, ELIM, FN, FNU, FNUL, HPI, RS1, SAR, TOL, YY, R1MACH
      INTEGER I, IFLAG, IN, INU, J, K, KODE, N, NAI, ND, NDAI, NLAST,
     * NN, NUF, NW, NZ, IDUM
      DIMENSION BRY(3), Y(N), CIP(4), CSS(3), CSR(3), CY(2)
      DATA CZERO,CONE,CI/(0.0E0,0.0E0),(1.0E0,0.0E0),(0.0E0,1.0E0)/
      DATA CIP(1),CIP(2),CIP(3),CIP(4)/
     1 (1.0E0,0.0E0), (0.0E0,1.0E0), (-1.0E0,0.0E0), (0.0E0,-1.0E0)/
      DATA HPI, AIC  /
     1      1.57079632679489662E+00,     1.265512123484645396E+00/
C
      NZ = 0
      ND = N
      NLAST = 0
C-----------------------------------------------------------------------
C     COMPUTED VALUES WITH EXPONENTS BETWEEN ALIM AND ELIM IN MAG-
C     NITUDE ARE SCALED TO KEEP INTERMEDIATE ARITHMETIC ON SCALE,
C     EXP(ALIM)=EXP(ELIM)*TOL
C-----------------------------------------------------------------------
      CSCL = CMPLX(1.0E0/TOL,0.0E0)
      CRSC = CMPLX(TOL,0.0E0)
      CSS(1) = CSCL
      CSS(2) = CONE
      CSS(3) = CRSC
      CSR(1) = CRSC
      CSR(2) = CONE
      CSR(3) = CSCL
      BRY(1) = 1.0E+3*R1MACH(1)/TOL
      YY = AIMAG(Z)
C-----------------------------------------------------------------------
C     ZN IS IN THE RIGHT HALF PLANE AFTER ROTATION BY CI OR -CI
C-----------------------------------------------------------------------
      ZN = -Z*CI
      ZB = Z
      CID = -CI
      INU = INT(FNU)
      ANG = HPI*(FNU-FLOAT(INU))
      CAR = COS(ANG)
      SAR = SIN(ANG)
      C2 = CMPLX(CAR,SAR)
      ZAR = C2
      IN = INU + N - 1
      IN = MOD(IN,4)
      C2 = C2*CIP(IN+1)
      IF (YY.GT.0.0E0) GO TO 10
      ZN = CONJG(-ZN)
      ZB = CONJG(ZB)
      CID = -CID
      C2 = CONJG(C2)
   10 CONTINUE
C-----------------------------------------------------------------------
C     CHECK FOR UNDERFLOW AND OVERFLOW ON FIRST MEMBER
C-----------------------------------------------------------------------
      FN = AMAX1(FNU,1.0E0)
      CALL CUNHJ(ZN, FN, 1, TOL, PHI, ARG, ZETA1, ZETA2, ASUM, BSUM)
      IF (KODE.EQ.1) GO TO 20
      CFN = CMPLX(FNU,0.0E0)
      S1 = -ZETA1 + CFN*(CFN/(ZB+ZETA2))
      GO TO 30
   20 CONTINUE
      S1 = -ZETA1 + ZETA2
   30 CONTINUE
      RS1 = REAL(S1)
      IF (ABS(RS1).GT.ELIM) GO TO 150
   40 CONTINUE
      NN = MIN0(2,ND)
      DO 90 I=1,NN
        FN = FNU + FLOAT(ND-I)
        CALL CUNHJ(ZN, FN, 0, TOL, PHI, ARG, ZETA1, ZETA2, ASUM, BSUM)
        IF (KODE.EQ.1) GO TO 50
        CFN = CMPLX(FN,0.0E0)
        AY = ABS(YY)
        S1 = -ZETA1 + CFN*(CFN/(ZB+ZETA2)) + CMPLX(0.0E0,AY)
        GO TO 60
   50   CONTINUE
        S1 = -ZETA1 + ZETA2
   60   CONTINUE
C-----------------------------------------------------------------------
C     TEST FOR UNDERFLOW AND OVERFLOW
C-----------------------------------------------------------------------
        RS1 = REAL(S1)
        IF (ABS(RS1).GT.ELIM) GO TO 120
        IF (I.EQ.1) IFLAG = 2
        IF (ABS(RS1).LT.ALIM) GO TO 70
C-----------------------------------------------------------------------
C     REFINE  TEST AND SCALE
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
        APHI = CABS(PHI)
        AARG = CABS(ARG)
        RS1 = RS1 + ALOG(APHI) - 0.25E0*ALOG(AARG) - AIC
        IF (ABS(RS1).GT.ELIM) GO TO 120
        IF (I.EQ.1) IFLAG = 1
        IF (RS1.LT.0.0E0) GO TO 70
        IF (I.EQ.1) IFLAG = 3
   70   CONTINUE
C-----------------------------------------------------------------------
C     SCALE S1 TO KEEP INTERMEDIATE ARITHMETIC ON SCALE NEAR
C     EXPONENT EXTREMES
C-----------------------------------------------------------------------
        CALL CAIRY(ARG, 0, 2, AI, NAI, IDUM)
        CALL CAIRY(ARG, 1, 2, DAI, NDAI, IDUM)
        S2 = PHI*(AI*ASUM+DAI*BSUM)
        C2R = REAL(S1)
        C2I = AIMAG(S1)
        C2M = EXP(C2R)*REAL(CSS(IFLAG))
        S1 = CMPLX(C2M,0.0E0)*CMPLX(COS(C2I),SIN(C2I))
        S2 = S2*S1
        IF (IFLAG.NE.1) GO TO 80
        CALL CUCHK(S2, NW, BRY(1), TOL)
        IF (NW.NE.0) GO TO 120
   80   CONTINUE
        IF (YY.LE.0.0E0) S2 = CONJG(S2)
        J = ND - I + 1
        S2 = S2*C2
        CY(I) = S2
        Y(J) = S2*CSR(IFLAG)
        C2 = C2*CID
   90 CONTINUE
      IF (ND.LE.2) GO TO 110
      RZ = CMPLX(2.0E0,0.0E0)/Z
      BRY(2) = 1.0E0/BRY(1)
      BRY(3) = R1MACH(2)
      S1 = CY(1)
      S2 = CY(2)
      C1 = CSR(IFLAG)
      ASCLE = BRY(IFLAG)
      K = ND - 2
      FN = FLOAT(K)
      DO 100 I=3,ND
        C2 = S2
        S2 = S1 + CMPLX(FNU+FN,0.0E0)*RZ*S2
        S1 = C2
        C2 = S2*C1
        Y(K) = C2
        K = K - 1
        FN = FN - 1.0E0
        IF (IFLAG.GE.3) GO TO 100
        C2R = REAL(C2)
        C2I = AIMAG(C2)
        C2R = ABS(C2R)
        C2I = ABS(C2I)
        C2M = AMAX1(C2R,C2I)
        IF (C2M.LE.ASCLE) GO TO 100
        IFLAG = IFLAG + 1
        ASCLE = BRY(IFLAG)
        S1 = S1*C1
        S2 = C2
        S1 = S1*CSS(IFLAG)
        S2 = S2*CSS(IFLAG)
        C1 = CSR(IFLAG)
  100 CONTINUE
  110 CONTINUE
      RETURN
  120 CONTINUE
      IF (RS1.GT.0.0E0) GO TO 140
C-----------------------------------------------------------------------
C     SET UNDERFLOW AND UPDATE PARAMETERS
C-----------------------------------------------------------------------
      Y(ND) = CZERO
      NZ = NZ + 1
      ND = ND - 1
      IF (ND.EQ.0) GO TO 110
      CALL CUOIK(Z, FNU, KODE, 1, ND, Y, NUF, TOL, ELIM, ALIM)
      IF (NUF.LT.0) GO TO 140
      ND = ND - NUF
      NZ = NZ + NUF
      IF (ND.EQ.0) GO TO 110
      FN = FNU + FLOAT(ND-1)
      IF (FN.LT.FNUL) GO TO 130
C      FN = AIMAG(CID)
C      J = NUF + 1
C      K = MOD(J,4) + 1
C      S1 = CIP(K)
C      IF (FN.LT.0.0E0) S1 = CONJG(S1)
C      C2 = C2*S1
      IN = INU + ND - 1
      IN = MOD(IN,4) + 1
      C2 = ZAR*CIP(IN)
      IF (YY.LE.0.0E0)C2=CONJG(C2)
      GO TO 40
  130 CONTINUE
      NLAST = ND
      RETURN
  140 CONTINUE
      NZ = -1
      RETURN
  150 CONTINUE
      IF (RS1.GT.0.0E0) GO TO 140
      NZ = N
      DO 160 I=1,N
        Y(I) = CZERO
  160 CONTINUE
      RETURN
      END
      SUBROUTINE CS1S2(ZR, S1, S2, NZ, ASCLE, ALIM, IUF)
C***BEGIN PROLOGUE  CS1S2
C***REFER TO  CBESK,CAIRY
C
C     CS1S2 TESTS FOR A POSSIBLE UNDERFLOW RESULTING FROM THE
C     ADDITION OF THE I AND K FUNCTIONS IN THE ANALYTIC CON-
C     TINUATION FORMULA WHERE S1=K FUNCTION AND S2=I FUNCTION.
C     ON KODE=1 THE I AND K FUNCTIONS ARE DIFFERENT ORDERS OF
C     MAGNITUDE, BUT FOR KODE=2 THEY CAN BE OF THE SAME ORDER
C     OF MAGNITUDE AND THE MAXIMUM MUST BE AT LEAST ONE
C     PRECISION ABOVE THE UNDERFLOW LIMIT.
C
C***ROUTINES CALLED  (NONE)
C***END PROLOGUE  CS1S2
      COMPLEX CZERO, C1, S1, S1D, S2, ZR
      REAL AA, ALIM, ALN, ASCLE, AS1, AS2, XX
      INTEGER IUF, NZ
      DATA CZERO / (0.0E0,0.0E0) /
      NZ = 0
      AS1 = CABS(S1)
      AS2 = CABS(S2)
      AA = REAL(S1)
      ALN = AIMAG(S1)
      IF (AA.EQ.0.0E0 .AND. ALN.EQ.0.0E0) GO TO 10
      IF (AS1.EQ.0.0E0) GO TO 10
      XX = REAL(ZR)
      ALN = -XX - XX + ALOG(AS1)
      S1D = S1
      S1 = CZERO
      AS1 = 0.0E0
      IF (ALN.LT.(-ALIM)) GO TO 10
      C1 = CLOG(S1D) - ZR - ZR
      S1 = CEXP(C1)
      AS1 = CABS(S1)
      IUF = IUF + 1
   10 CONTINUE
      AA = AMAX1(AS1,AS2)
      IF (AA.GT.ASCLE) RETURN
      S1 = CZERO
      S2 = CZERO
      NZ = 1
      IUF = 0
      RETURN
      END
      SUBROUTINE CSHCH(Z, CSH, CCH)
C***BEGIN PROLOGUE  CSHCH
C***REFER TO  CBESK,CBESH
C
C     CSHCH COMPUTES THE COMPLEX HYPERBOLIC FUNCTIONS CSH=SINH(X+I*Y)
C     AND CCH=COSH(X+I*Y), WHERE I**2=-1.
C
C***ROUTINES CALLED  (NONE)
C***END PROLOGUE  CSHCH
      COMPLEX CCH, CSH, Z
      REAL CCHI, CCHR, CH, CN, CSHI, CSHR, SH, SN, X, Y, COSH, SINH
      X = REAL(Z)
      Y = AIMAG(Z)
      SH = SINH(X)
      CH = COSH(X)
      SN = SIN(Y)
      CN = COS(Y)
      CSHR = SH*CN
      CSHI = CH*SN
      CSH = CMPLX(CSHR,CSHI)
      CCHR = CH*CN
      CCHI = SH*SN
      CCH = CMPLX(CCHR,CCHI)
      RETURN
      END
      SUBROUTINE CRATI(Z, FNU, N, CY, TOL)
C***BEGIN PROLOGUE  CRATI
C***REFER TO  CBESI,CBESK,CBESH
C
C     CRATI COMPUTES RATIOS OF I BESSEL FUNCTIONS BY BACKWARD
C     RECURRENCE.  THE STARTING INDEX IS DETERMINED BY FORWARD
C     RECURRENCE AS DESCRIBED IN J. RES. OF NAT. BUR. OF STANDARDS-B,
C     MATHEMATICAL SCIENCES, VOL 77B, P111-114, SEPTEMBER, 1973,
C     BESSEL FUNCTIONS I AND J OF COMPLEX ARGUMENT AND INTEGER ORDER,
C     BY D. J. SOOKNE.
C
C***ROUTINES CALLED  (NONE)
C***END PROLOGUE  CRATI
      COMPLEX CDFNU, CONE, CY, CZERO, PT, P1, P2, RZ, T1, Z
      REAL AK, AMAGZ, AP1, AP2, ARG, AZ, DFNU, FDNU, FLAM, FNU, FNUP,
     * RAP1, RHO, TEST, TEST1, TOL
      INTEGER I, ID, IDNU, INU, ITIME, K, KK, MAGZ, N
      DIMENSION CY(N)
      DATA CZERO, CONE / (0.0E0,0.0E0), (1.0E0,0.0E0) /
      AZ = CABS(Z)
      INU = INT(FNU)
      IDNU = INU + N - 1
      FDNU = FLOAT(IDNU)
      MAGZ = INT(AZ)
      AMAGZ = FLOAT(MAGZ+1)
      FNUP = AMAX1(AMAGZ,FDNU)
      ID = IDNU - MAGZ - 1
      ITIME = 1
      K = 1
      RZ = (CONE+CONE)/Z
      T1 = CMPLX(FNUP,0.0E0)*RZ
      P2 = -T1
      P1 = CONE
      T1 = T1 + RZ
      IF (ID.GT.0) ID = 0
      AP2 = CABS(P2)
      AP1 = CABS(P1)
C-----------------------------------------------------------------------
C     THE OVERFLOW TEST ON K(FNU+I-1,Z) BEFORE THE CALL TO CBKNX
C     GUARANTEES THAT P2 IS ON SCALE. SCALE TEST1 AND ALL SUBSEQUENT
C     P2 VALUES BY AP1 TO ENSURE THAT AN OVERFLOW DOES NOT OCCUR
C     PREMATURELY.
C-----------------------------------------------------------------------
      ARG = (AP2+AP2)/(AP1*TOL)
      TEST1 = SQRT(ARG)
      TEST = TEST1
      RAP1 = 1.0E0/AP1
      P1 = P1*CMPLX(RAP1,0.0E0)
      P2 = P2*CMPLX(RAP1,0.0E0)
      AP2 = AP2*RAP1
   10 CONTINUE
      K = K + 1
      AP1 = AP2
      PT = P2
      P2 = P1 - T1*P2
      P1 = PT
      T1 = T1 + RZ
      AP2 = CABS(P2)
      IF (AP1.LE.TEST) GO TO 10
      IF (ITIME.EQ.2) GO TO 20
      AK = CABS(T1)*0.5E0
      FLAM = AK + SQRT(AK*AK-1.0E0)
      RHO = AMIN1(AP2/AP1,FLAM)
      TEST = TEST1*SQRT(RHO/(RHO*RHO-1.0E0))
      ITIME = 2
      GO TO 10
   20 CONTINUE
      KK = K + 1 - ID
      AK = FLOAT(KK)
      DFNU = FNU + FLOAT(N-1)
      CDFNU = CMPLX(DFNU,0.0E0)
      T1 = CMPLX(AK,0.0E0)
      P1 = CMPLX(1.0E0/AP2,0.0E0)
      P2 = CZERO
      DO 30 I=1,KK
        PT = P1
        P1 = RZ*(CDFNU+T1)*P1 + P2
        P2 = PT
        T1 = T1 - CONE
   30 CONTINUE
      IF (REAL(P1).NE.0.0E0 .OR. AIMAG(P1).NE.0.0E0) GO TO 40
      P1 = CMPLX(TOL,TOL)
   40 CONTINUE
      CY(N) = P2/P1
      IF (N.EQ.1) RETURN
      K = N - 1
      AK = FLOAT(K)
      T1 = CMPLX(AK,0.0E0)
      CDFNU = CMPLX(FNU,0.0E0)*RZ
      DO 60 I=2,N
        PT = CDFNU + T1*RZ + CY(K+1)
        IF (REAL(PT).NE.0.0E0 .OR. AIMAG(PT).NE.0.0E0) GO TO 50
        PT = CMPLX(TOL,TOL)
   50   CONTINUE
        CY(K) = CONE/PT
        T1 = T1 - CONE
        K = K - 1
   60 CONTINUE
      RETURN
      END
      SUBROUTINE CBKNU(Z, FNU, KODE, N, Y, NZ, TOL, ELIM, ALIM)
C***BEGIN PROLOGUE  CBKNU
C***REFER TO  CBESI,CBESK,CAIRY,CBESH
C
C     CBKNU COMPUTES THE K BESSEL FUNCTION IN THE RIGHT HALF Z PLANE
C
C***ROUTINES CALLED  CKSCL,CSHCH,GAMLN,I1MACH,R1MACH,CUCHK
C***END PROLOGUE  CBKNU
C
      COMPLEX CCH, CK, COEF, CONE, CRSC, CS, CSCL, CSH, CSR, CSS, CTWO,
     * CZ, CZERO, F, FMU, P, PT, P1, P2, Q, RZ, SMU, ST, S1, S2, Y, Z,
     * ZD, CELM, CY
      REAL AA, AK, ALIM, ASCLE, A1, A2, BB, BK, BRY, CAZ, CC, DNU,
     * DNU2, ELIM, ETEST, FC, FHS, FK, FKS, FNU, FPI, G1, G2, HPI, PI,
     * P2I, P2M, P2R, RK, RTHPI, R1, S, SPI, TM, TOL, TTH, T1, T2, XX,
     * YY, GAMLN, R1MACH, HELIM, ELM, XD, YD, ALAS, AS
      INTEGER I, IDUM, IFLAG, INU, K, KFLAG, KK, KMAX, KODE, KODED, N,
     * NZ, I1MACH, NW, J, IC, INUB
      DIMENSION BRY(3), CC(8), CSS(3), CSR(3), Y(N), CY(2)
C
      DATA KMAX / 30 /
      DATA R1 / 2.0E0 /
      DATA CZERO,CONE,CTWO /(0.0E0,0.0E0),(1.0E0,0.0E0),(2.0E0,0.0E0)/
C
      DATA PI, RTHPI, SPI ,HPI, FPI, TTH /
     1     3.14159265358979324E0,       1.25331413731550025E0,
     2     1.90985931710274403E0,       1.57079632679489662E0,
     3     1.89769999331517738E0,       6.66666666666666666E-01/
C
      DATA CC(1), CC(2), CC(3), CC(4), CC(5), CC(6), CC(7), CC(8)/
     1     5.77215664901532861E-01,    -4.20026350340952355E-02,
     2    -4.21977345555443367E-02,     7.21894324666309954E-03,
     3    -2.15241674114950973E-04,    -2.01348547807882387E-05,
     4     1.13302723198169588E-06,     6.11609510448141582E-09/
C
      XX = REAL(Z)
      YY = AIMAG(Z)
      CAZ = CABS(Z)
      CSCL = CMPLX(1.0E0/TOL,0.0E0)
      CRSC = CMPLX(TOL,0.0E0)
      CSS(1) = CSCL
      CSS(2) = CONE
      CSS(3) = CRSC
      CSR(1) = CRSC
      CSR(2) = CONE
      CSR(3) = CSCL
      BRY(1) = 1.0E+3*R1MACH(1)/TOL
      BRY(2) = 1.0E0/BRY(1)
      BRY(3) = R1MACH(2)
      NZ = 0
      IFLAG = 0
      KODED = KODE
      RZ = CTWO/Z
      INU = INT(FNU+0.5E0)
      DNU = FNU - FLOAT(INU)
      IF (ABS(DNU).EQ.0.5E0) GO TO 110
      DNU2 = 0.0E0
      IF (ABS(DNU).GT.TOL) DNU2 = DNU*DNU
      IF (CAZ.GT.R1) GO TO 110
C-----------------------------------------------------------------------
C     SERIES FOR CABS(Z).LE.R1
C-----------------------------------------------------------------------
      FC = 1.0E0
      SMU = CLOG(RZ)
      FMU = SMU*CMPLX(DNU,0.0E0)
      CALL CSHCH(FMU, CSH, CCH)
      IF (DNU.EQ.0.0E0) GO TO 10
      FC = DNU*PI
      FC = FC/SIN(FC)
      SMU = CSH*CMPLX(1.0E0/DNU,0.0E0)
   10 CONTINUE
      A2 = 1.0E0 + DNU
C-----------------------------------------------------------------------
C     GAM(1-Z)*GAM(1+Z)=PI*Z/SIN(PI*Z), T1=1/GAM(1-DNU), T2=1/GAM(1+DNU)
C-----------------------------------------------------------------------
      T2 = EXP(-GAMLN(A2,IDUM))
      T1 = 1.0E0/(T2*FC)
      IF (ABS(DNU).GT.0.1E0) GO TO 40
C-----------------------------------------------------------------------
C     SERIES FOR F0 TO RESOLVE INDETERMINACY FOR SMALL ABS(DNU)
C-----------------------------------------------------------------------
      AK = 1.0E0
      S = CC(1)
      DO 20 K=2,8
        AK = AK*DNU2
        TM = CC(K)*AK
        S = S + TM
        IF (ABS(TM).LT.TOL) GO TO 30
   20 CONTINUE
   30 G1 = -S
      GO TO 50
   40 CONTINUE
      G1 = (T1-T2)/(DNU+DNU)
   50 CONTINUE
      G2 = 0.5E0*(T1+T2)*FC
      G1 = G1*FC
      F = CMPLX(G1,0.0E0)*CCH + SMU*CMPLX(G2,0.0E0)
      PT = CEXP(FMU)
      P = CMPLX(0.5E0/T2,0.0E0)*PT
      Q = CMPLX(0.5E0/T1,0.0E0)/PT
      S1 = F
      S2 = P
      AK = 1.0E0
      A1 = 1.0E0
      CK = CONE
      BK = 1.0E0 - DNU2
      IF (INU.GT.0 .OR. N.GT.1) GO TO 80
C-----------------------------------------------------------------------
C     GENERATE K(FNU,Z), 0.0D0 .LE. FNU .LT. 0.5D0 AND N=1
C-----------------------------------------------------------------------
      IF (CAZ.LT.TOL) GO TO 70
      CZ = Z*Z*CMPLX(0.25E0,0.0E0)
      T1 = 0.25E0*CAZ*CAZ
   60 CONTINUE
      F = (F*CMPLX(AK,0.0E0)+P+Q)*CMPLX(1.0E0/BK,0.0E0)
      P = P*CMPLX(1.0E0/(AK-DNU),0.0E0)
      Q = Q*CMPLX(1.0E0/(AK+DNU),0.0E0)
      RK = 1.0E0/AK
      CK = CK*CZ*CMPLX(RK,0.0)
      S1 = S1 + CK*F
      A1 = A1*T1*RK
      BK = BK + AK + AK + 1.0E0
      AK = AK + 1.0E0
      IF (A1.GT.TOL) GO TO 60
   70 CONTINUE
      Y(1) = S1
      IF (KODED.EQ.1) RETURN
      Y(1) = S1*CEXP(Z)
      RETURN
C-----------------------------------------------------------------------
C     GENERATE K(DNU,Z) AND K(DNU+1,Z) FOR FORWARD RECURRENCE
C-----------------------------------------------------------------------
   80 CONTINUE
      IF (CAZ.LT.TOL) GO TO 100
      CZ = Z*Z*CMPLX(0.25E0,0.0E0)
      T1 = 0.25E0*CAZ*CAZ
   90 CONTINUE
      F = (F*CMPLX(AK,0.0E0)+P+Q)*CMPLX(1.0E0/BK,0.0E0)
      P = P*CMPLX(1.0E0/(AK-DNU),0.0E0)
      Q = Q*CMPLX(1.0E0/(AK+DNU),0.0E0)
      RK = 1.0E0/AK
      CK = CK*CZ*CMPLX(RK,0.0E0)
      S1 = S1 + CK*F
      S2 = S2 + CK*(P-F*CMPLX(AK,0.0E0))
      A1 = A1*T1*RK
      BK = BK + AK + AK + 1.0E0
      AK = AK + 1.0E0
      IF (A1.GT.TOL) GO TO 90
  100 CONTINUE
      KFLAG = 2
      BK = REAL(SMU)
      A1 = FNU + 1.0E0
      AK = A1*ABS(BK)
      IF (AK.GT.ALIM) KFLAG = 3
      P2 = S2*CSS(KFLAG)
      S2 = P2*RZ
      S1 = S1*CSS(KFLAG)
      IF (KODED.EQ.1) GO TO 210
      F = CEXP(Z)
      S1 = S1*F
      S2 = S2*F
      GO TO 210
C-----------------------------------------------------------------------
C     IFLAG=0 MEANS NO UNDERFLOW OCCURRED
C     IFLAG=1 MEANS AN UNDERFLOW OCCURRED- COMPUTATION PROCEEDS WITH
C     KODED=2 AND A TEST FOR ON SCALE VALUES IS MADE DURING FORWARD
C     RECURSION
C-----------------------------------------------------------------------
  110 CONTINUE
      COEF = CMPLX(RTHPI,0.0E0)/CSQRT(Z)
      KFLAG = 2
      IF (KODED.EQ.2) GO TO 120
      IF (XX.GT.ALIM) GO TO 290
C     BLANK LINE
      A1 = EXP(-XX)*REAL(CSS(KFLAG))
      PT = CMPLX(A1,0.0E0)*CMPLX(COS(YY),-SIN(YY))
      COEF = COEF*PT
  120 CONTINUE
      IF (ABS(DNU).EQ.0.5E0) GO TO 300
C-----------------------------------------------------------------------
C     MILLER ALGORITHM FOR CABS(Z).GT.R1
C-----------------------------------------------------------------------
      AK = COS(PI*DNU)
      AK = ABS(AK)
      IF (AK.EQ.0.0E0) GO TO 300
      FHS = ABS(0.25E0-DNU2)
      IF (FHS.EQ.0.0E0) GO TO 300
C-----------------------------------------------------------------------
C     COMPUTE R2=F(E). IF CABS(Z).GE.R2, USE FORWARD RECURRENCE TO
C     DETERMINE THE BACKWARD INDEX K. R2=F(E) IS A STRAIGHT LINE ON
C     12.LE.E.LE.60. E IS COMPUTED FROM 2**(-E)=B**(1-I1MACH(11))=
C     TOL WHERE B IS THE BASE OF THE ARITHMETIC.
C-----------------------------------------------------------------------
      T1 = FLOAT(I1MACH(11)-1)*R1MACH(5)*3.321928094E0
      T1 = AMAX1(T1,12.0E0)
      T1 = AMIN1(T1,60.0E0)
      T2 = TTH*T1 - 6.0E0
      IF (XX.NE.0.0E0) GO TO 130
      T1 = HPI
      GO TO 140
  130 CONTINUE
      T1 = ATAN(YY/XX)
      T1 = ABS(T1)
  140 CONTINUE
      IF (T2.GT.CAZ) GO TO 170
C-----------------------------------------------------------------------
C     FORWARD RECURRENCE LOOP WHEN CABS(Z).GE.R2
C-----------------------------------------------------------------------
      ETEST = AK/(PI*CAZ*TOL)
      FK = 1.0E0
      IF (ETEST.LT.1.0E0) GO TO 180
      FKS = 2.0E0
      RK = CAZ + CAZ + 2.0E0
      A1 = 0.0E0
      A2 = 1.0E0
      DO 150 I=1,KMAX
        AK = FHS/FKS
        BK = RK/(FK+1.0E0)
        TM = A2
        A2 = BK*A2 - AK*A1
        A1 = TM
        RK = RK + 2.0E0
        FKS = FKS + FK + FK + 2.0E0
        FHS = FHS + FK + FK
        FK = FK + 1.0E0
        TM = ABS(A2)*FK
        IF (ETEST.LT.TM) GO TO 160
  150 CONTINUE
      GO TO 310
  160 CONTINUE
      FK = FK + SPI*T1*SQRT(T2/CAZ)
      FHS = ABS(0.25E0-DNU2)
      GO TO 180
  170 CONTINUE
C-----------------------------------------------------------------------
C     COMPUTE BACKWARD INDEX K FOR CABS(Z).LT.R2
C-----------------------------------------------------------------------
      A2 = SQRT(CAZ)
      AK = FPI*AK/(TOL*SQRT(A2))
      AA = 3.0E0*T1/(1.0E0+CAZ)
      BB = 14.7E0*T1/(28.0E0+CAZ)
      AK = (ALOG(AK)+CAZ*COS(AA)/(1.0E0+0.008E0*CAZ))/COS(BB)
      FK = 0.12125E0*AK*AK/CAZ + 1.5E0
  180 CONTINUE
      K = INT(FK)
C-----------------------------------------------------------------------
C     BACKWARD RECURRENCE LOOP FOR MILLER ALGORITHM
C-----------------------------------------------------------------------
      FK = FLOAT(K)
      FKS = FK*FK
      P1 = CZERO
      P2 = CMPLX(TOL,0.0E0)
      CS = P2
      DO 190 I=1,K
        A1 = FKS - FK
        A2 = (FKS+FK)/(A1+FHS)
        RK = 2.0E0/(FK+1.0E0)
        T1 = (FK+XX)*RK
        T2 = YY*RK
        PT = P2
        P2 = (P2*CMPLX(T1,T2)-P1)*CMPLX(A2,0.0E0)
        P1 = PT
        CS = CS + P2
        FKS = A1 - FK + 1.0E0
        FK = FK - 1.0E0
  190 CONTINUE
C-----------------------------------------------------------------------
C     COMPUTE (P2/CS)=(P2/CABS(CS))*(CONJG(CS)/CABS(CS)) FOR BETTER
C     SCALING
C-----------------------------------------------------------------------
      TM = CABS(CS)
      PT = CMPLX(1.0E0/TM,0.0E0)
      S1 = PT*P2
      CS = CONJG(CS)*PT
      S1 = COEF*S1*CS
      IF (INU.GT.0 .OR. N.GT.1) GO TO 200
      ZD = Z
      IF(IFLAG.EQ.1) GO TO 270
      GO TO 240
  200 CONTINUE
C-----------------------------------------------------------------------
C     COMPUTE P1/P2=(P1/CABS(P2)*CONJG(P2)/CABS(P2) FOR SCALING
C-----------------------------------------------------------------------
      TM = CABS(P2)
      PT = CMPLX(1.0E0/TM,0.0E0)
      P1 = PT*P1
      P2 = CONJG(P2)*PT
      PT = P1*P2
      S2 = S1*(CONE+(CMPLX(DNU+0.5E0,0.0E0)-PT)/Z)
C-----------------------------------------------------------------------
C     FORWARD RECURSION ON THE THREE TERM RECURSION RELATION WITH
C     SCALING NEAR EXPONENT EXTREMES ON KFLAG=1 OR KFLAG=3
C-----------------------------------------------------------------------
  210 CONTINUE
      CK = CMPLX(DNU+1.0E0,0.0E0)*RZ
      IF (N.EQ.1) INU = INU - 1
      IF (INU.GT.0) GO TO 220
      IF (N.EQ.1) S1=S2
      ZD = Z
      IF(IFLAG.EQ.1) GO TO 270
      GO TO 240
  220 CONTINUE
      INUB = 1
      IF (IFLAG.EQ.1) GO TO 261
  225 CONTINUE
      P1 = CSR(KFLAG)
      ASCLE = BRY(KFLAG)
      DO 230 I=INUB,INU
        ST = S2
        S2 = CK*S2 + S1
        S1 = ST
        CK = CK + RZ
        IF (KFLAG.GE.3) GO TO 230
        P2 = S2*P1
        P2R = REAL(P2)
        P2I = AIMAG(P2)
        P2R = ABS(P2R)
        P2I = ABS(P2I)
        P2M = AMAX1(P2R,P2I)
        IF (P2M.LE.ASCLE) GO TO 230
        KFLAG = KFLAG + 1
        ASCLE = BRY(KFLAG)
        S1 = S1*P1
        S2 = P2
        S1 = S1*CSS(KFLAG)
        S2 = S2*CSS(KFLAG)
        P1 = CSR(KFLAG)
  230 CONTINUE
      IF (N.EQ.1) S1 = S2
  240 CONTINUE
      Y(1) = S1*CSR(KFLAG)
      IF (N.EQ.1) RETURN
      Y(2) = S2*CSR(KFLAG)
      IF (N.EQ.2) RETURN
      KK = 2
  250 CONTINUE
      KK = KK + 1
      IF (KK.GT.N) RETURN
      P1 = CSR(KFLAG)
      ASCLE = BRY(KFLAG)
      DO 260 I=KK,N
        P2 = S2
        S2 = CK*S2 + S1
        S1 = P2
        CK = CK + RZ
        P2 = S2*P1
        Y(I) = P2
        IF (KFLAG.GE.3) GO TO 260
        P2R = REAL(P2)
        P2I = AIMAG(P2)
        P2R = ABS(P2R)
        P2I = ABS(P2I)
        P2M = AMAX1(P2R,P2I)
        IF (P2M.LE.ASCLE) GO TO 260
        KFLAG = KFLAG + 1
        ASCLE = BRY(KFLAG)
        S1 = S1*P1
        S2 = P2
        S1 = S1*CSS(KFLAG)
        S2 = S2*CSS(KFLAG)
        P1 = CSR(KFLAG)
  260 CONTINUE
      RETURN
C-----------------------------------------------------------------------
C     IFLAG=1 CASES, FORWARD RECURRENCE ON SCALED VALUES ON UNDERFLOW
C-----------------------------------------------------------------------
  261 CONTINUE
      HELIM = 0.5E0*ELIM
      ELM = EXP(-ELIM)
      CELM = CMPLX(ELM,0.0)
      ASCLE = BRY(1)
      ZD = Z
      XD = XX
      YD = YY
      IC = -1
      J = 2
      DO 262 I=1,INU
        ST = S2
        S2 = CK*S2+S1
        S1 = ST
        CK = CK+RZ
        AS = CABS(S2)
        ALAS = ALOG(AS)
        P2R = -XD+ALAS
        IF(P2R.LT.(-ELIM)) GO TO 263
        P2 = -ZD+CLOG(S2)
        P2R = REAL(P2)
        P2I = AIMAG(P2)
        P2M = EXP(P2R)/TOL
        P1 = CMPLX(P2M,0.0E0)*CMPLX(COS(P2I),SIN(P2I))
        CALL CUCHK(P1,NW,ASCLE,TOL)
        IF(NW.NE.0) GO TO 263
        J=3-J
        CY(J) = P1
        IF(IC.EQ.(I-1)) GO TO 264
        IC = I
        GO TO 262
  263   CONTINUE
        IF(ALAS.LT.HELIM) GO TO 262
        XD = XD-ELIM
        S1 = S1*CELM
        S2 = S2*CELM
        ZD = CMPLX(XD,YD)
  262 CONTINUE
      IF(N.EQ.1) S1 = S2
      GO TO 270
  264 CONTINUE
      KFLAG = 1
      INUB = I+1
      S2 = CY(J)
      J = 3 - J
      S1 = CY(J)
      IF(INUB.LE.INU) GO TO 225
      IF(N.EQ.1) S1 = S2
      GO TO 240
  270 CONTINUE
      Y(1) = S1
      IF (N.EQ.1) GO TO 280
      Y(2) = S2
  280 CONTINUE
      ASCLE = BRY(1)
      CALL CKSCL(ZD, FNU, N, Y, NZ, RZ, ASCLE, TOL, ELIM)
      INU = N - NZ
      IF (INU.LE.0) RETURN
      KK = NZ + 1
      S1 = Y(KK)
      Y(KK) = S1*CSR(1)
      IF (INU.EQ.1) RETURN
      KK = NZ + 2
      S2 = Y(KK)
      Y(KK) = S2*CSR(1)
      IF (INU.EQ.2) RETURN
      T2 = FNU + FLOAT(KK-1)
      CK = CMPLX(T2,0.0E0)*RZ
      KFLAG = 1
      GO TO 250
  290 CONTINUE
C-----------------------------------------------------------------------
C     SCALE BY EXP(Z), IFLAG = 1 CASES
C-----------------------------------------------------------------------
      KODED = 2
      IFLAG = 1
      KFLAG = 2
      GO TO 120
C-----------------------------------------------------------------------
C     FNU=HALF ODD INTEGER CASE, DNU=-0.5
C-----------------------------------------------------------------------
  300 CONTINUE
      S1 = COEF
      S2 = COEF
      GO TO 210
  310 CONTINUE
      NZ=-2
      RETURN
      END
      SUBROUTINE CKSCL(ZR, FNU, N, Y, NZ, RZ, ASCLE, TOL, ELIM)
C***BEGIN PROLOGUE  CKSCL
C***REFER TO  CBKNU,CUNK1,CUNK2
C
C     SET K FUNCTIONS TO ZERO ON UNDERFLOW, CONTINUE RECURRENCE
C     ON SCALED FUNCTIONS UNTIL TWO MEMBERS COME ON SCALE, THEN
C     RETURN WITH MIN(NZ+2,N) VALUES SCALED BY 1/TOL.
C
C***ROUTINES CALLED  CUCHK
C***END PROLOGUE  CKSCL
      COMPLEX CK, CS, CY, CZERO, RZ, S1, S2, Y, ZR, ZD, CELM
      REAL AA, ASCLE, ACS, AS, CSI, CSR, ELIM, FN, FNU, TOL, XX, ZRI,
     * ELM, ALAS, HELIM
      INTEGER I, IC, K, KK, N, NN, NW, NZ
      DIMENSION Y(N), CY(2)
      DATA CZERO / (0.0E0,0.0E0) /
C
      NZ = 0
      IC = 0
      XX = REAL(ZR)
      NN = MIN0(2,N)
      DO 10 I=1,NN
        S1 = Y(I)
        CY(I) = S1
        AS = CABS(S1)
        ACS = -XX + ALOG(AS)
        NZ = NZ + 1
        Y(I) = CZERO
        IF (ACS.LT.(-ELIM)) GO TO 10
        CS = -ZR + CLOG(S1)
        CSR = REAL(CS)
        CSI = AIMAG(CS)
        AA = EXP(CSR)/TOL
        CS = CMPLX(AA,0.0E0)*CMPLX(COS(CSI),SIN(CSI))
        CALL CUCHK(CS, NW, ASCLE, TOL)
        IF (NW.NE.0) GO TO 10
        Y(I) = CS
        NZ = NZ - 1
        IC = I
   10 CONTINUE
      IF (N.EQ.1) RETURN
      IF (IC.GT.1) GO TO 20
      Y(1) = CZERO
      NZ = 2
   20 CONTINUE
      IF (N.EQ.2) RETURN
      IF (NZ.EQ.0) RETURN
      FN = FNU + 1.0E0
      CK = CMPLX(FN,0.0E0)*RZ
      S1 = CY(1)
      S2 = CY(2)
      HELIM = 0.5E0*ELIM
      ELM = EXP(-ELIM)
      CELM = CMPLX(ELM,0.0E0)
      ZRI =AIMAG(ZR)
      ZD = ZR
C
C     FIND TWO CONSECUTIVE Y VALUES ON SCALE. SCALE RECURRENCE IF
C     S2 GETS LARGER THAN EXP(ELIM/2)
C
      DO 30 I=3,N
        KK = I
        CS = S2
        S2 = CK*S2 + S1
        S1 = CS
        CK = CK + RZ
        AS = CABS(S2)
        ALAS = ALOG(AS)
        ACS = -XX + ALAS
        NZ = NZ + 1
        Y(I) = CZERO
        IF (ACS.LT.(-ELIM)) GO TO 25
        CS = -ZD + CLOG(S2)
        CSR = REAL(CS)
        CSI = AIMAG(CS)
        AA = EXP(CSR)/TOL
        CS = CMPLX(AA,0.0E0)*CMPLX(COS(CSI),SIN(CSI))
        CALL CUCHK(CS, NW, ASCLE, TOL)
        IF (NW.NE.0) GO TO 25
        Y(I) = CS
        NZ = NZ - 1
        IF (IC.EQ.(KK-1)) GO TO 40
        IC = KK
        GO TO 30
   25   CONTINUE
        IF(ALAS.LT.HELIM) GO TO 30
        XX = XX-ELIM
        S1 = S1*CELM
        S2 = S2*CELM
        ZD = CMPLX(XX,ZRI)
   30 CONTINUE
      NZ = N
      IF(IC.EQ.N) NZ=N-1
      GO TO 45
   40 CONTINUE
      NZ = KK - 2
   45 CONTINUE
      DO 50 K=1,NZ
        Y(K) = CZERO
   50 CONTINUE
      RETURN
      END
      SUBROUTINE CACON(Z, FNU, KODE, MR, N, Y, NZ, RL, FNUL, TOL, ELIM,
     * ALIM)
C***BEGIN PROLOGUE  CACON
C***REFER TO  CBESK,CBESH
C
C     CACON APPLIES THE ANALYTIC CONTINUATION FORMULA
C
C         K(FNU,ZN*EXP(MP))=K(FNU,ZN)*EXP(-MP*FNU) - MP*I(FNU,ZN)
C                 MP=PI*MR*CMPLX(0.0,1.0)
C
C     TO CONTINUE THE K FUNCTION FROM THE RIGHT HALF TO THE LEFT
C     HALF Z PLANE
C
C***ROUTINES CALLED  CBINU,CBKNU,CS1S2,R1MACH
C***END PROLOGUE  CACON
      COMPLEX CK, CONE, CS, CSCL, CSCR, CSGN, CSPN, CSS, CSR, C1, C2,
     * RZ, SC1, SC2, ST, S1, S2, Y, Z, ZN, CY
      REAL ALIM, ARG, ASCLE, AS2, BSCLE, BRY, CPN, C1I, C1M, C1R, ELIM,
     * FMR, FNU, FNUL, PI, RL, SGN, SPN, TOL, YY, R1MACH
      INTEGER I, INU, IUF, KFLAG, KODE, MR, N, NN, NW, NZ
      DIMENSION Y(N), CY(2), CSS(3), CSR(3), BRY(3)
      DATA PI / 3.14159265358979324E0 /
      DATA CONE / (1.0E0,0.0E0) /
      NZ = 0
      ZN = -Z
      NN = N
      CALL CBINU(ZN, FNU, KODE, NN, Y, NW, RL, FNUL, TOL, ELIM, ALIM)
      IF (NW.LT.0) GO TO 80
C-----------------------------------------------------------------------
C     ANALYTIC CONTINUATION TO THE LEFT HALF PLANE FOR THE K FUNCTION
C-----------------------------------------------------------------------
      NN = MIN0(2,N)
      CALL CBKNU(ZN, FNU, KODE, NN, CY, NW, TOL, ELIM, ALIM)
      IF (NW.NE.0) GO TO 80
      S1 = CY(1)
      FMR = FLOAT(MR)
      SGN = -SIGN(PI,FMR)
      CSGN = CMPLX(0.0E0,SGN)
      IF (KODE.EQ.1) GO TO 10
      YY = -AIMAG(ZN)
      CPN = COS(YY)
      SPN = SIN(YY)
      CSGN = CSGN*CMPLX(CPN,SPN)
   10 CONTINUE
C-----------------------------------------------------------------------
C     CALCULATE CSPN=EXP(FNU*PI*I) TO MINIMIZE LOSSES OF SIGNIFICANCE
C     WHEN FNU IS LARGE
C-----------------------------------------------------------------------
      INU = INT(FNU)
      ARG = (FNU-FLOAT(INU))*SGN
      CPN = COS(ARG)
      SPN = SIN(ARG)
      CSPN = CMPLX(CPN,SPN)
      IF (MOD(INU,2).EQ.1) CSPN = -CSPN
      IUF = 0
      C1 = S1
      C2 = Y(1)
      ASCLE = 1.0E+3*R1MACH(1)/TOL
      IF (KODE.EQ.1) GO TO 20
      CALL CS1S2(ZN, C1, C2, NW, ASCLE, ALIM, IUF)
      NZ = NZ + NW
      SC1 = C1
   20 CONTINUE
      Y(1) = CSPN*C1 + CSGN*C2
      IF (N.EQ.1) RETURN
      CSPN = -CSPN
      S2 = CY(2)
      C1 = S2
      C2 = Y(2)
      IF (KODE.EQ.1) GO TO 30
      CALL CS1S2(ZN, C1, C2, NW, ASCLE, ALIM, IUF)
      NZ = NZ + NW
      SC2 = C1
   30 CONTINUE
      Y(2) = CSPN*C1 + CSGN*C2
      IF (N.EQ.2) RETURN
      CSPN = -CSPN
      RZ = CMPLX(2.0E0,0.0E0)/ZN
      CK = CMPLX(FNU+1.0E0,0.0E0)*RZ
C-----------------------------------------------------------------------
C     SCALE NEAR EXPONENT EXTREMES DURING RECURRENCE ON K FUNCTIONS
C-----------------------------------------------------------------------
      CSCL = CMPLX(1.0E0/TOL,0.0E0)
      CSCR = CMPLX(TOL,0.0E0)
      CSS(1) = CSCL
      CSS(2) = CONE
      CSS(3) = CSCR
      CSR(1) = CSCR
      CSR(2) = CONE
      CSR(3) = CSCL
      BRY(1) = ASCLE
      BRY(2) = 1.0E0/ASCLE
      BRY(3) = R1MACH(2)
      AS2 = CABS(S2)
      KFLAG = 2
      IF (AS2.GT.BRY(1)) GO TO 40
      KFLAG = 1
      GO TO 50
   40 CONTINUE
      IF (AS2.LT.BRY(2)) GO TO 50
      KFLAG = 3
   50 CONTINUE
      BSCLE = BRY(KFLAG)
      S1 = S1*CSS(KFLAG)
      S2 = S2*CSS(KFLAG)
      CS = CSR(KFLAG)
      DO 70 I=3,N
        ST = S2
        S2 = CK*S2 + S1
        S1 = ST
        C1 = S2*CS
        ST = C1
        C2 = Y(I)
        IF (KODE.EQ.1) GO TO 60
        IF (IUF.LT.0) GO TO 60
        CALL CS1S2(ZN, C1, C2, NW, ASCLE, ALIM, IUF)
        NZ = NZ + NW
        SC1 = SC2
        SC2 = C1
        IF (IUF.NE.3) GO TO 60
        IUF = -4
        S1 = SC1*CSS(KFLAG)
        S2 = SC2*CSS(KFLAG)
        ST = SC2
   60   CONTINUE
        Y(I) = CSPN*C1 + CSGN*C2
        CK = CK + RZ
        CSPN = -CSPN
        IF (KFLAG.GE.3) GO TO 70
        C1R = REAL(C1)
        C1I = AIMAG(C1)
        C1R = ABS(C1R)
        C1I = ABS(C1I)
        C1M = AMAX1(C1R,C1I)
        IF (C1M.LE.BSCLE) GO TO 70
        KFLAG = KFLAG + 1
        BSCLE = BRY(KFLAG)
        S1 = S1*CS
        S2 = ST
        S1 = S1*CSS(KFLAG)
        S2 = S2*CSS(KFLAG)
        CS = CSR(KFLAG)
   70 CONTINUE
      RETURN
   80 CONTINUE
      NZ = -1
      IF(NW.EQ.(-2)) NZ=-2
      RETURN
      END
      SUBROUTINE CBINU(Z, FNU, KODE, N, CY, NZ, RL, FNUL, TOL, ELIM,
     * ALIM)
C***BEGIN PROLOGUE  CBINU
C***REFER TO  CBESH,CBESI,CBESJ,CBESK,CAIRY,CBIRY
C
C     CBINU COMPUTES THE I FUNCTION IN THE RIGHT HALF Z PLANE
C
C***ROUTINES CALLED  CASYI,CBUNI,CMLRI,CSERI,CUOIK,CWRSK
C***END PROLOGUE  CBINU
      COMPLEX CW, CY, CZERO, Z
      REAL ALIM, AZ, DFNU, ELIM, FNU, FNUL, RL, TOL
      INTEGER I, INW, KODE, N, NLAST, NN, NUI, NW, NZ
      DIMENSION CY(N), CW(2)
      DATA CZERO / (0.0E0,0.0E0) /
C
      NZ = 0
      AZ = CABS(Z)
      NN = N
      DFNU = FNU + FLOAT(N-1)
      IF (AZ.LE.2.0E0) GO TO 10
      IF (AZ*AZ*0.25E0.GT.DFNU+1.0E0) GO TO 20
   10 CONTINUE
C-----------------------------------------------------------------------
C     POWER SERIES
C-----------------------------------------------------------------------
      CALL CSERI(Z, FNU, KODE, NN, CY, NW, TOL, ELIM, ALIM)
      INW = IABS(NW)
      NZ = NZ + INW
      NN = NN - INW
      IF (NN.EQ.0) RETURN
      IF (NW.GE.0) GO TO 120
      DFNU = FNU + FLOAT(NN-1)
   20 CONTINUE
      IF (AZ.LT.RL) GO TO 40
      IF (DFNU.LE.1.0E0) GO TO 30
      IF (AZ+AZ.LT.DFNU*DFNU) GO TO 50
C-----------------------------------------------------------------------
C     ASYMPTOTIC EXPANSION FOR LARGE Z
C-----------------------------------------------------------------------
   30 CONTINUE
      CALL CASYI(Z, FNU, KODE, NN, CY, NW, RL, TOL, ELIM, ALIM)
      IF (NW.LT.0) GO TO 130
      GO TO 120
   40 CONTINUE
      IF (DFNU.LE.1.0E0) GO TO 70
   50 CONTINUE
C-----------------------------------------------------------------------
C     OVERFLOW AND UNDERFLOW TEST ON I SEQUENCE FOR MILLER ALGORITHM
C-----------------------------------------------------------------------
      CALL CUOIK(Z, FNU, KODE, 1, NN, CY, NW, TOL, ELIM, ALIM)
      IF (NW.LT.0) GO TO 130
      NZ = NZ + NW
      NN = NN - NW
      IF (NN.EQ.0) RETURN
      DFNU = FNU+FLOAT(NN-1)
      IF (DFNU.GT.FNUL) GO TO 110
      IF (AZ.GT.FNUL) GO TO 110
   60 CONTINUE
      IF (AZ.GT.RL) GO TO 80
   70 CONTINUE
C-----------------------------------------------------------------------
C     MILLER ALGORITHM NORMALIZED BY THE SERIES
C-----------------------------------------------------------------------
      CALL CMLRI(Z, FNU, KODE, NN, CY, NW, TOL)
      IF(NW.LT.0) GO TO 130
      GO TO 120
   80 CONTINUE
C-----------------------------------------------------------------------
C     MILLER ALGORITHM NORMALIZED BY THE WRONSKIAN
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C     OVERFLOW TEST ON K FUNCTIONS USED IN WRONSKIAN
C-----------------------------------------------------------------------
      CALL CUOIK(Z, FNU, KODE, 2, 2, CW, NW, TOL, ELIM, ALIM)
      IF (NW.GE.0) GO TO 100
      NZ = NN
      DO 90 I=1,NN
        CY(I) = CZERO
   90 CONTINUE
      RETURN
  100 CONTINUE
      IF (NW.GT.0) GO TO 130
      CALL CWRSK(Z, FNU, KODE, NN, CY, NW, CW, TOL, ELIM, ALIM)
      IF (NW.LT.0) GO TO 130
      GO TO 120
  110 CONTINUE
C-----------------------------------------------------------------------
C     INCREMENT FNU+NN-1 UP TO FNUL, COMPUTE AND RECUR BACKWARD
C-----------------------------------------------------------------------
      NUI = INT(FNUL-DFNU) + 1
      NUI = MAX0(NUI,0)
      CALL CBUNI(Z, FNU, KODE, NN, CY, NW, NUI, NLAST, FNUL, TOL, ELIM,
     * ALIM)
      IF (NW.LT.0) GO TO 130
      NZ = NZ + NW
      IF (NLAST.EQ.0) GO TO 120
      NN = NLAST
      GO TO 60
  120 CONTINUE
      RETURN
  130 CONTINUE
      NZ = -1
      IF(NW.EQ.(-2)) NZ=-2
      RETURN
      END
      FUNCTION GAMLN(Z,IERR)
C***BEGIN PROLOGUE  GAMLN
C***DATE WRITTEN   830501   (YYMMDD)
C***REVISION DATE  830501   (YYMMDD)
C***CATEGORY NO.  B5F
C***KEYWORDS  GAMMA FUNCTION,LOGARITHM OF GAMMA FUNCTION
C***AUTHOR  AMOS, DONALD E., SANDIA NATIONAL LABORATORIES
C***PURPOSE  TO COMPUTE THE LOGARITHM OF THE GAMMA FUNCTION
C***DESCRIPTION
C
C         GAMLN COMPUTES THE NATURAL LOG OF THE GAMMA FUNCTION FOR
C         Z.GT.0.  THE ASYMPTOTIC EXPANSION IS USED TO GENERATE VALUES
C         GREATER THAN ZMIN WHICH ARE ADJUSTED BY THE RECURSION
C         G(Z+1)=Z*G(Z) FOR Z.LE.ZMIN.  THE FUNCTION WAS MADE AS
C         PORTABLE AS POSSIBLE BY COMPUTIMG ZMIN FROM THE NUMBER OF BASE
C         10 DIGITS IN A WORD, RLN=AMAX1(-ALOG10(R1MACH(4)),0.5E-18)
C         LIMITED TO 18 DIGITS OF (RELATIVE) ACCURACY.
C
C         SINCE INTEGER ARGUMENTS ARE COMMON, A TABLE LOOK UP ON 100
C         VALUES IS USED FOR SPEED OF EXECUTION.
C
C     DESCRIPTION OF ARGUMENTS
C
C         INPUT
C           Z      - REAL ARGUMENT, Z.GT.0.0E0
C
C         OUTPUT
C           GAMLN  - NATURAL LOG OF THE GAMMA FUNCTION AT Z
C           IERR   - ERROR FLAG
C                    IERR=0, NORMAL RETURN, COMPUTATION COMPLETED
C                    IERR=1, Z.LE.0.0E0,    NO COMPUTATION
C
C***REFERENCES  COMPUTATION OF BESSEL FUNCTIONS OF COMPLEX ARGUMENT
C                 BY D. E. AMOS, SAND83-0083, MAY, 1983.
C***ROUTINES CALLED  I1MACH,R1MACH
C***END PROLOGUE  GAMLN
C
      INTEGER I, I1M, K, MZ, NZ, IERR, I1MACH
      REAL CF, CON, FLN, FZ, GLN, RLN, S, TLG, TRM, TST, T1, WDTOL, Z,
     * ZDMY, ZINC, ZM, ZMIN, ZP, ZSQ
      REAL R1MACH
      DIMENSION CF(22), GLN(100)
C           LNGAMMA(N), N=1,100
      DATA GLN(1), GLN(2), GLN(3), GLN(4), GLN(5), GLN(6), GLN(7),
     1     GLN(8), GLN(9), GLN(10), GLN(11), GLN(12), GLN(13), GLN(14),
     2     GLN(15), GLN(16), GLN(17), GLN(18), GLN(19), GLN(20),
     3     GLN(21), GLN(22)/
     4     0.00000000000000000E+00,     0.00000000000000000E+00,
     5     6.93147180559945309E-01,     1.79175946922805500E+00,
     6     3.17805383034794562E+00,     4.78749174278204599E+00,
     7     6.57925121201010100E+00,     8.52516136106541430E+00,
     8     1.06046029027452502E+01,     1.28018274800814696E+01,
     9     1.51044125730755153E+01,     1.75023078458738858E+01,
     A     1.99872144956618861E+01,     2.25521638531234229E+01,
     B     2.51912211827386815E+01,     2.78992713838408916E+01,
     C     3.06718601060806728E+01,     3.35050734501368889E+01,
     D     3.63954452080330536E+01,     3.93398841871994940E+01,
     E     4.23356164607534850E+01,     4.53801388984769080E+01/
      DATA GLN(23), GLN(24), GLN(25), GLN(26), GLN(27), GLN(28),
     1     GLN(29), GLN(30), GLN(31), GLN(32), GLN(33), GLN(34),
     2     GLN(35), GLN(36), GLN(37), GLN(38), GLN(39), GLN(40),
     3     GLN(41), GLN(42), GLN(43), GLN(44)/
     4     4.84711813518352239E+01,     5.16066755677643736E+01,
     5     5.47847293981123192E+01,     5.80036052229805199E+01,
     6     6.12617017610020020E+01,     6.45575386270063311E+01,
     7     6.78897431371815350E+01,     7.12570389671680090E+01,
     8     7.46582363488301644E+01,     7.80922235533153106E+01,
     9     8.15579594561150372E+01,     8.50544670175815174E+01,
     A     8.85808275421976788E+01,     9.21361756036870925E+01,
     B     9.57196945421432025E+01,     9.93306124547874269E+01,
     C     1.02968198614513813E+02,     1.06631760260643459E+02,
     D     1.10320639714757395E+02,     1.14034211781461703E+02,
     E     1.17771881399745072E+02,     1.21533081515438634E+02/
      DATA GLN(45), GLN(46), GLN(47), GLN(48), GLN(49), GLN(50),
     1     GLN(51), GLN(52), GLN(53), GLN(54), GLN(55), GLN(56),
     2     GLN(57), GLN(58), GLN(59), GLN(60), GLN(61), GLN(62),
     3     GLN(63), GLN(64), GLN(65), GLN(66)/
     4     1.25317271149356895E+02,     1.29123933639127215E+02,
     5     1.32952575035616310E+02,     1.36802722637326368E+02,
     6     1.40673923648234259E+02,     1.44565743946344886E+02,
     7     1.48477766951773032E+02,     1.52409592584497358E+02,
     8     1.56360836303078785E+02,     1.60331128216630907E+02,
     9     1.64320112263195181E+02,     1.68327445448427652E+02,
     A     1.72352797139162802E+02,     1.76395848406997352E+02,
     B     1.80456291417543771E+02,     1.84533828861449491E+02,
     C     1.88628173423671591E+02,     1.92739047287844902E+02,
     D     1.96866181672889994E+02,     2.01009316399281527E+02,
     E     2.05168199482641199E+02,     2.09342586752536836E+02/
      DATA GLN(67), GLN(68), GLN(69), GLN(70), GLN(71), GLN(72),
     1     GLN(73), GLN(74), GLN(75), GLN(76), GLN(77), GLN(78),
     2     GLN(79), GLN(80), GLN(81), GLN(82), GLN(83), GLN(84),
     3     GLN(85), GLN(86), GLN(87), GLN(88)/
     4     2.13532241494563261E+02,     2.17736934113954227E+02,
     5     2.21956441819130334E+02,     2.26190548323727593E+02,
     6     2.30439043565776952E+02,     2.34701723442818268E+02,
     7     2.38978389561834323E+02,     2.43268849002982714E+02,
     8     2.47572914096186884E+02,     2.51890402209723194E+02,
     9     2.56221135550009525E+02,     2.60564940971863209E+02,
     A     2.64921649798552801E+02,     2.69291097651019823E+02,
     B     2.73673124285693704E+02,     2.78067573440366143E+02,
     C     2.82474292687630396E+02,     2.86893133295426994E+02,
     D     2.91323950094270308E+02,     2.95766601350760624E+02,
     E     3.00220948647014132E+02,     3.04686856765668715E+02/
      DATA GLN(89), GLN(90), GLN(91), GLN(92), GLN(93), GLN(94),
     1     GLN(95), GLN(96), GLN(97), GLN(98), GLN(99), GLN(100)/
     2     3.09164193580146922E+02,     3.13652829949879062E+02,
     3     3.18152639620209327E+02,     3.22663499126726177E+02,
     4     3.27185287703775217E+02,     3.31717887196928473E+02,
     5     3.36261181979198477E+02,     3.40815058870799018E+02,
     6     3.45379407062266854E+02,     3.49954118040770237E+02,
     7     3.54539085519440809E+02,     3.59134205369575399E+02/
C             COEFFICIENTS OF ASYMPTOTIC EXPANSION
      DATA CF(1), CF(2), CF(3), CF(4), CF(5), CF(6), CF(7), CF(8),
     1     CF(9), CF(10), CF(11), CF(12), CF(13), CF(14), CF(15),
     2     CF(16), CF(17), CF(18), CF(19), CF(20), CF(21), CF(22)/
     3     8.33333333333333333E-02,    -2.77777777777777778E-03,
     4     7.93650793650793651E-04,    -5.95238095238095238E-04,
     5     8.41750841750841751E-04,    -1.91752691752691753E-03,
     6     6.41025641025641026E-03,    -2.95506535947712418E-02,
     7     1.79644372368830573E-01,    -1.39243221690590112E+00,
     8     1.34028640441683920E+01,    -1.56848284626002017E+02,
     9     2.19310333333333333E+03,    -3.61087712537249894E+04,
     A     6.91472268851313067E+05,    -1.52382215394074162E+07,
     B     3.82900751391414141E+08,    -1.08822660357843911E+10,
     C     3.47320283765002252E+11,    -1.23696021422692745E+13,
     D     4.88788064793079335E+14,    -2.13203339609193739E+16/
C
C             LN(2*PI)
      DATA CON                    /     1.83787706640934548E+00/
C
C***FIRST EXECUTABLE STATEMENT  GAMLN
      IERR=0
      IF (Z.LE.0.0E0) GO TO 70
      IF (Z.GT.101.0E0) GO TO 10
      NZ = INT(Z)
      FZ = Z - FLOAT(NZ)
      IF (FZ.GT.0.0E0) GO TO 10
      IF (NZ.GT.100) GO TO 10
      GAMLN = GLN(NZ)
      RETURN
   10 CONTINUE
      WDTOL = R1MACH(4)
      WDTOL = AMAX1(WDTOL,0.5E-18)
      I1M = I1MACH(11)
      RLN = R1MACH(5)*FLOAT(I1M)
      FLN = AMIN1(RLN,20.0E0)
      FLN = AMAX1(FLN,3.0E0)
      FLN = FLN - 3.0E0
      ZM = 1.8000E0 + 0.3875E0*FLN
      MZ = INT(ZM) + 1
      ZMIN = FLOAT(MZ)
      ZDMY = Z
      ZINC = 0.0E0
      IF (Z.GE.ZMIN) GO TO 20
      ZINC = ZMIN - FLOAT(NZ)
      ZDMY = Z + ZINC
   20 CONTINUE
      ZP = 1.0E0/ZDMY
      T1 = CF(1)*ZP
      S = T1
      IF (ZP.LT.WDTOL) GO TO 40
      ZSQ = ZP*ZP
      TST = T1*WDTOL
      DO 30 K=2,22
        ZP = ZP*ZSQ
        TRM = CF(K)*ZP
        IF (ABS(TRM).LT.TST) GO TO 40
        S = S + TRM
   30 CONTINUE
   40 CONTINUE
      IF (ZINC.NE.0.0E0) GO TO 50
      TLG = ALOG(Z)
      GAMLN = Z*(TLG-1.0E0) + 0.5E0*(CON-TLG) + S
      RETURN
   50 CONTINUE
      ZP = 1.0E0
      NZ = INT(ZINC)
      DO 60 I=1,NZ
        ZP = ZP*(Z+FLOAT(I-1))
   60 CONTINUE
      TLG = ALOG(ZDMY)
      GAMLN = ZDMY*(TLG-1.0E0) - ALOG(ZP) + 0.5E0*(CON-TLG) + S
      RETURN
C
C
   70 CONTINUE
      IERR=1
      RETURN
      END
      SUBROUTINE CUCHK(Y, NZ, ASCLE, TOL)
C***BEGIN PROLOGUE  CUCHK
C***REFER TO CSERI,CUOIK,CUNK1,CUNK2,CUNI1,CUNI2,CKSCL
C
C      Y ENTERS AS A SCALED QUANTITY WHOSE MAGNITUDE IS GREATER THAN
C      EXP(-ALIM)=ASCLE=1.0E+3*R1MACH(1)/TOL. THE TEST IS MADE TO SEE
C      IF THE MAGNITUDE OF THE REAL OR IMAGINARY PART WOULD UNDER FLOW
C      WHEN Y IS SCALED (BY TOL) TO ITS PROPER VALUE. Y IS ACCEPTED
C      IF THE UNDERFLOW IS AT LEAST ONE PRECISION BELOW THE MAGNITUDE
C      OF THE LARGEST COMPONENT; OTHERWISE THE PHASE ANGLE DOES NOT HAVE
C      ABSOLUTE ACCURACY AND AN UNDERFLOW IS ASSUMED.
C
C***ROUTINES CALLED  (NONE)
C***END PROLOGUE  CUCHK
C
      COMPLEX Y
      REAL ASCLE, SS, ST, TOL, YR, YI
      INTEGER NZ
      NZ = 0
      YR = REAL(Y)
      YI = AIMAG(Y)
      YR = ABS(YR)
      YI = ABS(YI)
      ST = AMIN1(YR,YI)
      IF (ST.GT.ASCLE) RETURN
      SS = AMAX1(YR,YI)
      ST=ST/TOL
      IF (SS.LT.ST) NZ = 1
      RETURN
      END
      SUBROUTINE CACAI(Z, FNU, KODE, MR, N, Y, NZ, RL, TOL, ELIM, ALIM)
C***BEGIN PROLOGUE  CACAI
C***REFER TO  CAIRY
C
C     CACAI APPLIES THE ANALYTIC CONTINUATION FORMULA
C
C         K(FNU,ZN*EXP(MP))=K(FNU,ZN)*EXP(-MP*FNU) - MP*I(FNU,ZN)
C                 MP=PI*MR*CMPLX(0.0,1.0)
C
C     TO CONTINUE THE K FUNCTION FROM THE RIGHT HALF TO THE LEFT
C     HALF Z PLANE FOR USE WITH CAIRY WHERE FNU=1/3 OR 2/3 AND N=1.
C     CACAI IS THE SAME AS CACON WITH THE PARTS FOR LARGER ORDERS AND
C     RECURRENCE REMOVED. A RECURSIVE CALL TO CACON CAN RESULT IF CACON
C     IS CALLED FROM CAIRY.
C
C***ROUTINES CALLED  CASYI,CBKNU,CMLRI,CSERI,CS1S2,R1MACH
C***END PROLOGUE  CACAI
      COMPLEX CSGN, CSPN, C1, C2, Y, Z, ZN, CY
      REAL ALIM, ARG, ASCLE, AZ, CPN, DFNU, ELIM, FMR, FNU, PI, RL,
     * SGN, SPN, TOL, YY, R1MACH
      INTEGER INU, IUF, KODE, MR, N, NN, NW, NZ
      DIMENSION Y(N), CY(2)
      DATA PI / 3.14159265358979324E0 /
      NZ = 0
      ZN = -Z
      AZ = CABS(Z)
      NN = N
      DFNU = FNU + FLOAT(N-1)
      IF (AZ.LE.2.0E0) GO TO 10
      IF (AZ*AZ*0.25E0.GT.DFNU+1.0E0) GO TO 20
   10 CONTINUE
C-----------------------------------------------------------------------
C     POWER SERIES FOR THE I FUNCTION
C-----------------------------------------------------------------------
      CALL CSERI(ZN, FNU, KODE, NN, Y, NW, TOL, ELIM, ALIM)
      GO TO 40
   20 CONTINUE
      IF (AZ.LT.RL) GO TO 30
C-----------------------------------------------------------------------
C     ASYMPTOTIC EXPANSION FOR LARGE Z FOR THE I FUNCTION
C-----------------------------------------------------------------------
      CALL CASYI(ZN, FNU, KODE, NN, Y, NW, RL, TOL, ELIM, ALIM)
      IF (NW.LT.0) GO TO 70
      GO TO 40
   30 CONTINUE
C-----------------------------------------------------------------------
C     MILLER ALGORITHM NORMALIZED BY THE SERIES FOR THE I FUNCTION
C-----------------------------------------------------------------------
      CALL CMLRI(ZN, FNU, KODE, NN, Y, NW, TOL)
      IF(NW.LT.0) GO TO 70
   40 CONTINUE
C-----------------------------------------------------------------------
C     ANALYTIC CONTINUATION TO THE LEFT HALF PLANE FOR THE K FUNCTION
C-----------------------------------------------------------------------
      CALL CBKNU(ZN, FNU, KODE, 1, CY, NW, TOL, ELIM, ALIM)
      IF (NW.NE.0) GO TO 70
      FMR = FLOAT(MR)
      SGN = -SIGN(PI,FMR)
      CSGN = CMPLX(0.0E0,SGN)
      IF (KODE.EQ.1) GO TO 50
      YY = -AIMAG(ZN)
      CPN = COS(YY)
      SPN = SIN(YY)
      CSGN = CSGN*CMPLX(CPN,SPN)
   50 CONTINUE
C-----------------------------------------------------------------------
C     CALCULATE CSPN=EXP(FNU*PI*I) TO MINIMIZE LOSSES OF SIGNIFICANCE
C     WHEN FNU IS LARGE
C-----------------------------------------------------------------------
      INU = INT(FNU)
      ARG = (FNU-FLOAT(INU))*SGN
      CPN = COS(ARG)
      SPN = SIN(ARG)
      CSPN = CMPLX(CPN,SPN)
      IF (MOD(INU,2).EQ.1) CSPN = -CSPN
      C1 = CY(1)
      C2 = Y(1)
      IF (KODE.EQ.1) GO TO 60
      IUF = 0
      ASCLE = 1.0E+3*R1MACH(1)/TOL
      CALL CS1S2(ZN, C1, C2, NW, ASCLE, ALIM, IUF)
      NZ = NZ + NW
   60 CONTINUE
      Y(1) = CSPN*C1 + CSGN*C2
      RETURN
   70 CONTINUE
      NZ = -1
      IF(NW.EQ.(-2)) NZ=-2
      RETURN
      END
      SUBROUTINE XERROR(MESS,NMESS,L1,L2)
C
C     THIS IS A DUMMY XERROR ROUTINE TO PRINT ERROR MESSAGES WITH NMESS
C     CHARACTERS. L1 AND L2 ARE DUMMY PARAMETERS TO MAKE THIS CALL
C     COMPATIBLE WITH THE SLATEC XERROR ROUTINE. THIS IS A FORTRAN 77
C     ROUTINE.
C
      CHARACTER*(*) MESS
      NN=NMESS/70
      NR=NMESS-70*NN
      IF(NR.NE.0) NN=NN+1
      K=1
      PRINT 900
  900 FORMAT(/)
      DO 10 I=1,NN
        KMIN=MIN0(K+69,NMESS)
        PRINT *, MESS(K:KMIN)
        K=K+70
   10 CONTINUE
      PRINT 900
      RETURN
      END
C*** cqcbes.f
-----------------------------------------------------------------
C>>>  CQCBES.FOR:  Single precision quick check programs
-----------------------------------------------------------------

      PROGRAM CQCBH
C
C***DATE WRITTEN   830501   (YYMMDD)
C***REVISION DATE  890801, 930101   (YYMMDD)
C
C     CQCBH IS A QUICK CHECK ROUTINE FOR THE COMPLEX H BESSEL FUNCTIONS
C     GENERATED BY SUBROUTINE CBESH.
C
C     CQCBH GENERATES SEQUENCES OF H BESSEL FUNCTIONS FOR KIND 2 FROM
C     CBESH AND CHECKS THEM AGAINST ANALYTIC CONTINUATION FORMULAS
C     IN THE (Z,FNU) SPACE:
C
C     KODE = 1 TESTS (ANALYTIC CONTINUATION FORMULAE, I**2 = -1):
C
C     H(FNU,2,Z)=-EXP(I*PI*FNU)*H(FNU,1,-Z),       -PI.LT.ARG(Z).LE.0
C
C               = 2*COS(PI*FNU)*H(FNU,2,-Z) + EXP(I*PI*FNU)*H(FNU,1,-Z),
C
C                                                   0.LT.ARG(Z).LE.PI
C
C     KODE = 2 TESTS FOR KINDS 1 AND 2:
C
C            EXP(-I*Z)*H(FNU,1,Z) = [EXP(-I*Z)*H(FNU,1,Z)]
C
C            EXP( I*Z)*H(FNU,2,Z) = [EXP( I*Z)*H(FNU,2,Z)]
C
C     WHERE THE LEFT SIDE IS COMPUTED WITH KODE = 1 AND THE RIGHT SIDE
C     WITH KODE = 2.
C
C     THE PARAMETER MQC CAN HAVE VALUES 1 (THE DEFAULT) FOR A FASTER,
C     LESS DEFINITIVE TEST OR 2 FOR A SLOWER, MORE DEFINITIVE TEST.
C
C     MACHINE CONSTANTS ARE DEFINED IN FUNCTIONS I1MACH, R1MACH, AND
C     D1MACH. THESE MUST BE SELECTED BY THE USER OR SET ACCORDING TO
C     PROLOGUE INSTRUCTIONS.
C
      COMPLEX CW, CI, U, V, W, Y, Z, ZN, CSGN
      REAL AA, AB, AER, ALIM, ATOL, AV, CT, DIG, ERR, XX, YY,
     * ELIM, EPS, ER, ERTOL, FNU, FNUL, PI, R, RL,
     * RM, R1M4, R1M5, R2, ST, T, TOL, TS, XNU, R1MACH, SLAK, FILM
      INTEGER I, ICASE, IERR, IHP, IL, IR, IRB, IT, ITL, K, KODE, KK,
     *K1, K2, LFLG, LUN, MFLG, M, N, NU, NZ1, NZ2, NZ3, I1MACH, KEPS,
     *MQC, NL, NUL, KDO
      DIMENSION T(20), AER(20), XNU(20), U(20), V(20), W(20), Y(20),
     *KEPS(20), KDO(20)
      DATA LUN /7/
      PARAMETER (MQC=1)
      OPEN(LUN,FILE='CQCBH.OUT')
C-----------------------------------------------------------------------
C     SET PARAMETERS RELATED TO MACHINE CONSTANTS.
C     TOL IS THE APPROXIMATE UNIT ROUNDOFF LIMITED TO 1.0E-18.
C     ELIM IS THE APPROXIMATE EXPONENTIAL OVER- AND UNDERFLOW LIMIT.
C     EXP(-ELIM).LT.EXP(-ALIM)=EXP(-ELIM)/TOL    AND
C     EXP(ELIM).GT.EXP(ALIM)=EXP(ELIM)*TOL       ARE INTERVALS NEAR
C     UNDERFLOW AND OVERFLOW LIMITS WHERE SCALED ARITHMETIC IS DONE.
C     RL IS THE LOWER BOUNDARY OF THE ASYMPTOTIC EXPANSION FOR LARGE Z.
C     DIG = NUMBER OF BASE 10 DIGITS IN TOL = 10**(-DIG).
C     FNUL IS THE LOWER BOUNDARY OF THE ASYMPTOTIC SERIES FOR LARGE FNU
C-----------------------------------------------------------------------
      R1M4 = R1MACH(4)
      TOL = AMAX1(R1M4,1.0E-18)
      AA = -ALOG10(R1M4)
      K1 = I1MACH(12)
      K2 = I1MACH(13)
      R1M5 = R1MACH(5)
      K = MIN0(IABS(K1),IABS(K2))
      ELIM = 2.303E0*(FLOAT(K)*R1M5-3.0E0)
      AB = AA*2.303E0
      ALIM = ELIM + AMAX1(-AB,-41.45E0)
      DIG = AMIN1(AA,18.0E0)
      FNUL = 10.0E0 + 6.0E0*(DIG-3.0E0)
      RL = 1.2E0*DIG + 3.0E0
      SLAK = 3.0E0+4.0E0*(-ALOG10(TOL)-7.0E0)/11.0E0
      SLAK = AMAX1(SLAK,3.0E0)
      ERTOL = TOL*10.0E0**SLAK
      RM = 0.5E0*(ALIM + ELIM)
      RM = AMIN1(RM,200.0E0)
      RM = AMAX1(RM,RL+10.0E0)
      R2 = AMIN1(FNUL,RM)
C-----------------------------------------------------------------------
      WRITE (LUN,99999)
99999 FORMAT (' QUICK CHECK ROUTINE FOR THE H BESSEL FUNCTIONS FROM CBES
     *H'/)
      WRITE (LUN,99998)
99998 FORMAT (' PARAMETERS TOL,ELIM,ALIM,RL,FNUL,DIG')
      WRITE (LUN,99997) TOL, ELIM, ALIM, RL, FNUL, DIG
99997 FORMAT (6E12.4/)
      ATOL = 100.0E0*TOL
      PI = 4.0E0*ATAN(1.0E0)
      CI = CMPLX(0.0E0,1.0E0)
      WRITE (LUN,99996) MQC
99996 FORMAT (/' CHECKS IN THE (Z,FNU) SPACE WITH MQC = ',I2/)
C-----------------------------------------------------------------------
C     TEST VALUES OF Z IN -PI.LT.ARG(Z).LE.PI NEAR FORMULA BOUNDARIES
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C     KDO(K), K=1,IL  DETERMINES WHICH OF THE IL ANGLES IN -PI TO PI
C     ARE USE TO COMPUTE VALUES OF Z
C       KDO(K) = 0  MEANS THAT THE INDEX K WILL BE USED FOR ONE OR TWO
C                   VALUES OF Z, DEPENDING ON THE CHOICE OF KEPS(K)
C              = 1  MEANS THAT THE INDEX K AND THE CORRESPONDING ANGLE
C                   WILL BE SKIPPED
C     KEPS(K), K=1,IL DETERMINES WHICH OF THE ANGLES GET INCREMENTED
C     UP AND DOWN TO PUT VALUES OF Z IN REGIONS WHERE DIFFERENT
C     FORMULAE ARE USED.
C       KEPS(K) =0  MEANS THAT THE ANGLE WILL BE USED WITHOUT CHANGE
C               =1  MEANS THAT THE ANGLE WILL BE INCREMENTED UP AND
C                   DOWN BY EPS
C     THE ANGLES TO BE USED ARE STORED IN THE T(I) ARRAY, I=1,ITL
C-----------------------------------------------------------------------
      IF (MQC.NE.2) THEN
        NL=2
        IL=5
        DO 5 I=1,IL
          KEPS(I)=0
          KDO(I)=0
    5   CONTINUE
        NUL=5
        XNU(1) = 0.0E0
        XNU(2) = 1.0E0
        XNU(3) = 2.0E0
        XNU(4) = 0.5E0*FNUL
        XNU(5) = FNUL + 1.1E0
      ELSE
        NL=4
        IL=13
        DO 6 I=1,IL
          KDO(I)=0
          KEPS(I)=0
    6   CONTINUE
        KDO(2)=1
        KDO(6)=1
        KDO(8)=1
        KDO(12)=1
        KEPS(3)=1
        KEPS(4)=1
        KEPS(5)=1
        KEPS(9)=1
        KEPS(10)=1
        KEPS(11)=1
        NUL=6
        XNU(1) = 0.0E0
        XNU(2) = 0.6E0
        XNU(3) = 1.3E0
        XNU(4) = 2.0E0
        XNU(5) = 0.5E0*FNUL
        XNU(6) = FNUL + 1.1E0
      ENDIF
      I = 2
      EPS = 0.01E0
      FILM=FLOAT(IL-1)
      T(1) = -PI + EPS
      DO 30 K=2,IL
        IF (KDO(K).EQ.0) THEN
          T(I) = PI*FLOAT(-IL+2*K-1)/FILM
          IF (KEPS(K).EQ.0) GO TO 20
          TS=T(I)
          T(I) = TS - EPS
          I = I + 1
          T(I) = TS + EPS
   20     CONTINUE
          I = I + 1
        ENDIF
   30 CONTINUE
      ITL = I - 1
      LFLG = 0
      DO 170 KODE=1,2
        DO 160 N=1,NL
          DO 150 NU=1,NUL
            FNU = XNU(NU)
            DO 140 ICASE=1,3
              IRB = MIN0(ICASE,2)
              DO 130 IR=IRB,3
                GO TO (50, 60, 70), ICASE
   50           CONTINUE
                R = (EPS*FLOAT(3-IR)+2.0E0*FLOAT(IR-1))/2.0E0
                GO TO 80
   60           CONTINUE
                R = (2.0E0*FLOAT(3-IR)+R2*FLOAT(IR-1))/2.0E0
                GO TO 80
   70           CONTINUE
                IF (R2.GE.RM) GO TO 140
                R = (R2*FLOAT(3-IR)+RM*FLOAT(IR-1))/2.0E0
   80           CONTINUE
                DO 120 IT=1,ITL
                  CT = COS(T(IT))
                  ST = SIN(T(IT))
                  IF (ABS(CT).LT.ATOL) CT = 0.0E0
                  IF (ABS(ST).LT.ATOL) ST = 0.0E0
                  XX = R*CT
                  YY = R*ST
                  Z = CMPLX(XX,YY)
                  IF (KODE.EQ.1) THEN
                    M=2
                    CALL CBESH(Z,FNU,KODE,M,N,Y,NZ1,IERR)
                    IF (IERR.NE.0.OR.NZ1.NE.0) GO TO 120
                    IF (ST.LT.0.0E0 .OR. (ST.EQ.0.0E0.AND.CT.GT.0.0E0))
     *              THEN
                      IHP = 1
                      ZN=-Z
                      M=1
                      CALL CBESH(ZN,FNU,KODE,M,N,W,NZ2,IERR)
                      IF (IERR.NE.0.OR.NZ2.NE.0) GO TO 120
                    ELSE
                      IHP = 2
                      ZN=-Z
                      M=2
                      CALL CBESH(ZN,FNU,KODE,M,N,W,NZ3,IERR)
                      IF (IERR.NE.0.OR.NZ3.NE.0) GO TO 120
                      M=1
                      CALL CBESH(ZN,FNU,KODE,M,N,V,NZ2,IERR)
                      IF (IERR.NE.0.OR.NZ2.NE.0) GO TO 120
                    ENDIF
                    AB=AMOD(FNU,2.0E0)*PI
                    CSGN = CMPLX(COS(AB),SIN(AB))
                    MFLG = 0
                    DO 100 I=1,N
                      AB = FNU+FLOAT(I-1)
                      AA = MAX(0.5E0,AB)
                      IF (IHP.EQ.1) THEN
                        V(I) = -CSGN*W(I)
                        CW = Y(I) - V(I)
                      ELSE
                        V(I) = 2.0E0*REAL(CSGN)*W(I) + CSGN*V(I)
                        CW = Y(I) - V(I)
                      ENDIF
                      AV = CABS(Y(I))
                      ER = CABS(CW)
                      IF (YY.EQ.0.0E0) THEN
                        IF (ABS(XX).LT.AA) ER = ER/AV
                      ELSE
                        ER = ER/AV
                      ENDIF
                      AER(I) = ER
                      IF (ER.GT.ERTOL) MFLG = 1
                      CSGN = -CSGN
  100               CONTINUE
                  ELSE
                    M=1
                    KK=1
                    CALL CBESH(Z,FNU,KK,M,N,U,NZ1,IERR)
                    IF (IERR.NE.0.OR.NZ1.NE.0) GO TO 120
                    CALL CBESH(Z,FNU,KODE,M,N,V,NZ2,IERR)
                    IF (IERR.NE.0.OR.NZ2.NE.0) GO TO 120
                    M=2
                    KK=1
                    CALL CBESH(Z,FNU,KK,M,N,W,NZ1,IERR)
                    IF (IERR.NE.0.OR.NZ1.NE.0) GO TO 120
                    CALL CBESH(Z,FNU,KODE,M,N,Y,NZ2,IERR)
                    IF (IERR.NE.0.OR.NZ2.NE.0) GO TO 120
                    ZN=CI*Z
                    ZN=EXP(ZN)
                    MFLG = 0
                    DO 105 I=1,N
                      AB = FNU+FLOAT(I-1)
                      AA = MAX(0.5E0,AB)
                      CW = U(I)/ZN-V(I)
                      ER = CABS(CW)
                      AV = CABS(V(I))
                      IF (YY.EQ.0.0E0) THEN
                        IF (ABS(XX).LT.AA) ER = ER/AV
                      ELSE
                        ER = ER/AV
                      ENDIF
                      ERR = ER
                      IF (ER.GT.ERTOL) MFLG = 1
                      CW=ZN*W(I)-Y(I)
                      ER = CABS(CW)
                      AV = CABS(Y(I))
                      IF (YY.EQ.0.0E0) THEN
                        IF (ABS(XX).LT.AA) ER = ER/AV
                      ELSE
                        ER = ER/AV
                      ENDIF
                      IF (ER.GT.ERTOL) MFLG = 1
                      AER(I) = ER+ERR
  105               CONTINUE
                  ENDIF
                  IF (MFLG.EQ.0) GO TO 120
                  IF (LFLG.EQ.1) GO TO 110
                  WRITE (LUN,99995) ERTOL
99995             FORMAT (/' CASES WHICH VIOLATE THE RELATIVE ERROR TEST
     * WITH ERTOL =', E12.4/)
                  WRITE (LUN,99994)
99994             FORMAT (/' OUTPUT FORMAT'/' KODE,N,IR,IT,ICASE')
                  WRITE (LUN,99993)
99993             FORMAT (' ER(K),K=1,N'/' Z,FNU,V(1),Y(1)')
                  LFLG = 1
  110             CONTINUE
                  WRITE (LUN,99992) KODE, N, IR, IT, ICASE
99992             FORMAT (5I5)
                  WRITE (LUN,99991) (AER(K),K=1,N)
                  WRITE (LUN,99991) Z, FNU, V(1), Y(1)
99991             FORMAT (7E12.4)
  120           CONTINUE
  130         CONTINUE
  140       CONTINUE
  150     CONTINUE
  160   CONTINUE
  170 CONTINUE
      IF (LFLG.EQ.0) WRITE (LUN,99990)
99990 FORMAT (/' QUICK CHECKS OK'/)
      STOP
      END
      PROGRAM CQCBI
C
C***DATE WRITTEN   830501   (YYMMDD)
C***REVISION DATE  890801, 930101   (YYMMDD)
C
C     CQCBI IS A QUICK CHECK ROUTINE FOR THE COMPLEX I BESSEL FUNCTION
C     GENERATED BY SUBROUTINE CBESI.
C
C     CQCBK GENERATES SEQUENCES OF I AND K BESSEL FUNCTIONS FROM
C     CBESI AND CBESK AND CHECKS THE WRONSKIAN EVALUATION
C
C           I(FNU,Z)*K(FNU+1,Z) + I(FNU+1,Z)*K(FNU,Z) = 1/Z
C
C     IN THE RIGHT HALF PLANE AND A MODIFIED FORM
C
C          I(FNU+1,Z)*K(FNU,ZR) - I(FNU,Z)*K(FNU+1,ZR) = C/Z
C
C     IN THE LEFT HALF PLANE WHERE ZR=-Z AND C=EXP(I*FNU*SGN) WITH
C     SGN=+1 FOR IM(Z).GE.0 AND SGN=-1 FOR IM(Z).LT.0.
C
C     THE PARAMETER MQC CAN HAVE VALUES 1 (THE DEFAULT) FOR A FASTER,
C     LESS DEFINITIVE TEST OR 2 FOR A SLOWER, MORE DEFINITIVE TEST.
C
C     MACHINE CONSTANTS ARE DEFINED IN FUNCTIONS I1MACH, R1MACH, AND
C     D1MACH. THESE MUST BE SELECTED BY THE USER OR SET ACCORDING TO
C     PROLOGUE INSTRUCTIONS.
C
      COMPLEX  CONE, CSGN, CC, CV, CW, CY, W, Y, Z, ZR
      REAL AA, AB, AER, ALIM, ARG, ATOL, DIG, ELIM, EPS,
     * ER, ERTOL, FFNU, FNU, FNUL, GNU, HPI, PI, R, RL, R1M4, R1M5,
     * R2, RM, T, TOL, XX, YY, R1MACH, SLAK, TS, ST, CT, FILM, XNU
      INTEGER I, ICASE, IFNU, IL, IPRNT, IR, IT, ITL, IRB, K, KK, KODE,
     * K1, K2, LFLG, LUN, MFLG, N, NU, NZ, N1, NUL, IERR,
     * MQC, NL, KEPS, KDO, I1MACH
      DIMENSION T(20), AER(20), Y(20), W(20), XNU(20), KEPS(20), KDO(20)
      DATA LUN /7/
      PARAMETER (MQC=1)
      OPEN(LUN,FILE='CQCBI.OUT')
C-----------------------------------------------------------------------
C     SET PARAMETERS RELATED TO MACHINE CONSTANTS.
C     TOL IS THE APPROXIMATE UNIT ROUNDOFF LIMITED TO 1.0E-18.
C     ELIM IS THE APPROXIMATE EXPONENTIAL OVER- AND UNDERFLOW LIMIT.
C     EXP(-ELIM).LT.EXP(-ALIM)=EXP(-ELIM)/TOL    AND
C     EXP(ELIM).GT.EXP(ALIM)=EXP(ELIM)*TOL       ARE INTERVALS NEAR
C     UNDERFLOW AND OVERFLOW LIMITS WHERE SCALED ARITHMETIC IS DONE.
C     RL IS THE LOWER BOUNDARY OF THE ASYMPTOTIC EXPANSION FOR LARGE Z.
C     DIG = NUMBER OF BASE 10 DIGITS IN TOL = 10**(-DIG).
C     FNUL IS THE LOWER BOUNDARY OF THE ASYMPTOTIC SERIES FOR LARGE FNU
C-----------------------------------------------------------------------
      R1M4 = R1MACH(4)
      TOL = MAX(R1M4,1.0E-18)
      AA = -ALOG10(R1M4)
      K1 = I1MACH(12)
      K2 = I1MACH(13)
      R1M5 = R1MACH(5)
      K = MIN(ABS(K1),ABS(K2))
      ELIM = 2.303E0*(FLOAT(K)*R1M5-3.0E0)
      AB = AA*2.303E0
      ALIM = ELIM + MAX(-AB,-41.45E0)
      DIG = MIN(AA,18.0E0)
      FNUL = 10.0E0 + 6.0E0*(DIG-3.0E0)
      RL = 1.2E0*DIG + 3.0E0
      SLAK = 3.0E0+4.0E0*(-ALOG10(TOL)-7.0E0)/11.0E0
      SLAK = MAX(SLAK,3.0E0)
      ERTOL = TOL*10.0E0**SLAK
      RM = 0.5E0*(ALIM + ELIM)
      RM = MIN(RM,200.0E0)
      RM = MAX(RM,RL+10.0E0)
      R2 = MIN(FNUL,RM)
C-----------------------------------------------------------------------
      WRITE (LUN,99999)
99999 FORMAT (' QUICK CHECK ROUTINE FOR THE I BESSEL FUNCTION FROM CBESI
     *'/)
      WRITE (LUN,99998)
99998 FORMAT (' PARAMETERS TOL,ELIM,ALIM,RL,FNUL,DIG')
      WRITE (LUN,99997) TOL, ELIM, ALIM, RL, FNUL, DIG
99997 FORMAT (6E12.4/)
      CONE = CMPLX(1.0E0,0.0E0)
      ATOL = 100.0E0*TOL
      HPI = 2.0E0*ATAN(1.0E0)
      PI = HPI + HPI
      WRITE (LUN,99996) MQC
99996 FORMAT (/' CHECKS IN THE (Z,FNU) SPACE WITH MQC = ',I2/)
C-----------------------------------------------------------------------
C     TEST VALUES OF Z IN -PI.LT.ARG(Z).LE.PI NEAR FORMULA BOUNDARIES
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C     KDO(K), K=1,IL  DETERMINES WHICH OF THE IL ANGLES IN -PI TO PI
C     ARE USE TO COMPUTE VALUES OF Z
C       KDO(K) = 0  MEANS THAT THE INDEX K WILL BE USED FOR ONE OR TWO
C                   VALUES OF Z, DEPENDING ON THE CHOICE OF KEPS(K)
C              = 1  MEANS THAT THE INDEX K AND THE CORRESPONDING ANGLE
C                   WILL BE SKIPPED
C     KEPS(K), K=1,IL DETERMINES WHICH OF THE ANGLES GET INCREMENTED
C     UP AND DOWN TO PUT VALUES OF Z IN REGIONS WHERE DIFFERENT
C     FORMULAE ARE USED.
C       KEPS(K) =0  MEANS THAT THE ANGLE WILL BE USED WITHOUT CHANGE
C               =1  MEANS THAT THE ANGLE WILL BE INCREMENTED UP AND
C                   DOWN BY EPS
C     THE ANGLES TO BE USED ARE STORED IN THE T(I) ARRAY, I=1,ITL
C-----------------------------------------------------------------------
      IF (MQC.NE.2) THEN
        NL=2
        IL=5
        DO 5 I=1,IL
          KEPS(I)=0
          KDO(I)=0
    5   CONTINUE
        NUL=5
        XNU(1) = 0.0E0
        XNU(2) = 1.0E0
        XNU(3) = 2.0E0
        XNU(4) = 0.5E0*FNUL
        XNU(5) = FNUL + 1.1E0
      ELSE
        NL=4
        IL=13
        DO 6 I=1,IL
          KDO(I)=0
          KEPS(I)=0
    6   CONTINUE
        KDO(2)=1
        KDO(6)=1
        KDO(8)=1
        KDO(12)=1
        KEPS(3)=1
        KEPS(4)=1
        KEPS(5)=1
        KEPS(9)=1
        KEPS(10)=1
        KEPS(11)=1
        NUL=6
        XNU(1) = 0.0E0
        XNU(2) = 0.6E0
        XNU(3) = 1.3E0
        XNU(4) = 2.0E0
        XNU(5) = 0.5E0*FNUL
        XNU(6) = FNUL + 1.1E0
      ENDIF
      I = 2
      EPS = 0.01E0
      FILM=FLOAT(IL-1)
      T(1) = -PI + EPS
      DO 30 K=2,IL
        IF (KDO(K).EQ.0) THEN
          T(I) = PI*FLOAT(-IL+2*K-1)/FILM
          IF (KEPS(K).EQ.0) GO TO 20
          TS=T(I)
          T(I) = TS - EPS
          I = I + 1
          T(I) = TS + EPS
   20     CONTINUE
          I = I + 1
        ENDIF
   30 CONTINUE
      ITL = I - 1
      LFLG = 0
      DO 200 KODE=1,2
        DO 190 N=1,NL
          N1 = N + 1
          DO 180 NU=1,NUL
            FNU = XNU(NU)
            IFNU = INT(FNU)
            FFNU = FNU - FLOAT(IFNU)
            ARG = PI*FFNU
            CSGN = CMPLX(COS(ARG),SIN(ARG))
            IF (MOD(IFNU,2).EQ.1) CSGN = -CSGN
            DO 170 ICASE=1,3
              IRB = MIN(2,ICASE)
              DO 160 IR=IRB,4
                GO TO (50, 60, 70), ICASE
   50           CONTINUE
                R = (0.2E0*FLOAT(4-IR)+2.0E0*FLOAT(IR-1))/3.0E0
                GO TO 80
   60           CONTINUE
                R = (2.0E0*FLOAT(4-IR)+R2*FLOAT(IR-1))/3.0E0
                GO TO 80
   70           CONTINUE
                IF (R2.GE.RM) GO TO 170
                R = (R2*FLOAT(4-IR)+RM*FLOAT(IR-1))/3.0E0
   80           CONTINUE
                DO 150 IT=1,ITL
                  CT = COS(T(IT))
                  ST = SIN(T(IT))
                  IF (ABS(CT).LT.ATOL) CT = 0.0E0
                  IF (ABS(ST).LT.ATOL) ST = 0.0E0
                  XX = R*CT
                  YY = R*ST
                  Z = CMPLX(XX,YY)
                  IF (CT.GE.0.0E0) THEN
C-----------------------------------------------------------------------
C     WRONSKIAN CHECKS IN THE RIGHT HALF PLANE
C-----------------------------------------------------------------------
                    CALL CBESI(Z, FNU, KODE, N1, Y, NZ, IERR)
                    IF (NZ.NE.0.OR.IERR.NE.0) GO TO 150
                    CALL CBESK(Z, FNU, KODE, N1, W, NZ, IERR)
                    IF (NZ.NE.0.OR.IERR.NE.0) GO TO 150
C-----------------------------------------------------------------------
C     ADJUSTMENTS TO WRONSKIAN DUE TO SCALING OF I AND K FUNCTIONS
C     ON KODE=2
C-----------------------------------------------------------------------
                    CV = CONE/Z
                    IF (KODE.EQ.2) THEN
                      CW = CMPLX(COS(YY),SIN(YY))
                      CV = CW*CV
                    ENDIF
                    CC = CONE
                  ELSE
C-----------------------------------------------------------------------
C     WRONSKIAN CHECKS IN THE LEFT HALF PLANE
C-----------------------------------------------------------------------
                    ZR = -Z
                    CALL CBESI(Z, FNU, KODE, N1, Y, NZ, IERR)
                    IF (NZ.NE.0.OR.IERR.NE.0) GO TO 150
                    CALL CBESK(ZR, FNU, KODE, N1, W, NZ, IERR)
                    IF (NZ.NE.0.OR.IERR.NE.0) GO TO 150
                    CV = CSGN
                    IF (YY.LT.0.0E0) THEN
                      CV = CONJG(CV)
                    ENDIF
                    CV = CV/Z
                    IF (KODE.EQ.2) THEN
C-----------------------------------------------------------------------
C     ADJUSTMENTS TO WRONSKIAN DUE TO SCALING OF I AND K FUNCTIONS
C     ON KODE=2. SCALE FACTOR = EXP(-I*YY) FOR RE(Z).LT.0
C-----------------------------------------------------------------------
                      CW = CMPLX(COS(YY),-SIN(YY))
                      CV = CV*CW
                    ENDIF
                    CC = -CONE
                  ENDIF
                  MFLG = 0
                  KK=0
                  DO 130 I=1,N
                    CW = W(I)*Y(I+1)
                    CY = CC*W(I+1)*Y(I)
                    CY = CY + CW - CV
                    ER = CABS(CY)/CABS(CV)
                    AER(I) = ER
                    IF (ER.GT.ERTOL) THEN
                      IF(KK.EQ.0) THEN
                        MFLG = 1
                        KK=I
                      ENDIF
                    ENDIF
                    IF (CT.LT.0.0E0) CV = -CV
  130             CONTINUE
                  IF (MFLG.EQ.0) GO TO 150
                  IF (LFLG.EQ.1) GO TO 140
                  WRITE (LUN,99995) ERTOL
99995             FORMAT (/' CASES WHICH VIOLATE THE RELATIVE ERROR TEST
     * WITH ERTOL =', E12.4/)
                  WRITE (LUN,99994)
99994             FORMAT (/' OUTPUT FORMAT'/' KODE,N,IR,IT,ICASE,KK')
                  WRITE (LUN,99993)
99993             FORMAT (' ER(K),K=1,N'/' Z,FNU,Y(KK)        KK=INDEX O
     *F FIRST NON-ZERO PAIR'/)
                  LFLG = 1
  140             CONTINUE
                  WRITE (LUN,99992) KODE, N, IR, IT, ICASE, KK
99992             FORMAT (6I5)
                  WRITE (LUN,99991) (AER(K),K=1,N)
                  WRITE (LUN,99991) Z, FNU, Y(KK)
99991             FORMAT (6E12.4)
  150           CONTINUE
  160         CONTINUE
  170       CONTINUE
  180     CONTINUE
  190   CONTINUE
  200 CONTINUE
      IF (LFLG.EQ.0) WRITE (LUN,99990)
99990 FORMAT (/' QUICK CHECKS OK'/)
      IF (MQC.EQ.1) STOP
C-----------------------------------------------------------------------
C     CHECKS NEAR UNDERFLOW LIMITS ON SERIES(I=1) AND UNIFORM
C     ASYMPTOTIC EXPANSION(I=2)
C-----------------------------------------------------------------------
      WRITE (LUN,99989)
99989 FORMAT (/' CHECKS NEAR UNDERFLOW AND OVERFLOW LIMITS'/)
      Z = CMPLX(1.4E0,1.4E0)
      IPRNT = 0
      DO 280 I=1,2
        FNU = 10.2E0
        KODE = 1
        N = 20
  230   CONTINUE
        CALL CBESI(Z, FNU, KODE, N, Y, NZ, IERR)
        IF (NZ.NE.0) GO TO 240
        FNU = FNU + 5.0E0
        GO TO 230
  240   CONTINUE
        IF (NZ.LT.10) GO TO 250
        FNU = FNU - 1.0E0
        GO TO 230
  250   CONTINUE
        CALL CBESK(Z, FNU, KODE, 2, W, NZ, IERR)
        CV = CONE/Z
        CY = W(1)*Y(2)
        CW = W(2)*Y(1)
        CW = CW + CY - CV
        ER = CABS(CW)/CABS(CV)
        IF (ER.LT.ERTOL) GO TO 270
        IF (IPRNT.EQ.1) GO TO 260
        WRITE (LUN,99988)
99988   FORMAT (/' OUTPUT FORMAT'/' ERROR,Z,FNU,KODE,N'/)
        IPRNT = 1
  260   CONTINUE
        WRITE (LUN,99987) ER, Z, FNU, KODE, N
99987   FORMAT (4E12.4, 2I5)
  270   CONTINUE
        XX = RL + RL
        Z = CMPLX(XX,0.0E0)
  280 CONTINUE
C-----------------------------------------------------------------------
C     CHECK NEAR OVERFLOW LIMITS
C-----------------------------------------------------------------------
      Z = CMPLX(ELIM,0.0E0)
      FNU = 0.0
  290 CONTINUE
      CALL CBESK(Z, FNU, KODE, N, Y, NZ, IERR)
      IF (NZ.LT.10) GO TO 300
      IF (NZ.EQ.N) FNU = FNU + 3.0E0
      FNU = FNU + 2.0E0
      GO TO 290
  300 CONTINUE
      GNU = FNU + FLOAT(N-2)
      CALL CBESI(Z, GNU, KODE, 2, W, NZ, IERR)
      CV = CONE/Z
      CY = Y(N-1)*W(2)
      CW = Y(N)*W(1)
      CW = CW + CY - CV
      ER = CABS(CW)/CABS(CV)
      IF (ER.LT.ERTOL) GO TO 320
      IF (IPRNT.EQ.1) GO TO 310
      WRITE (LUN,99988)
      IPRNT = 1
  310 CONTINUE
      WRITE (LUN,99987) ER, Z, FNU, KODE, N
  320 CONTINUE
      IF (IPRNT.EQ.0) WRITE (LUN,99990)
      STOP
      END
      PROGRAM CQCBJ
C
C***DATE WRITTEN   830501   (YYMMDD)
C***REVISION DATE  890801, 930101   (YYMMDD)
C
C     CQCBJ IS A QUICK CHECK ROUTINE FOR THE COMPLEX J BESSEL FUNCTION
C     GENERATED BY SUBROUTINE CBESJ.
C
C     CQCBJ GENERATES SEQUENCES OF J AND H BESSEL FUNCTIONS FROM CBESJ
C     AND CBESH AND CHECKS THE WRONSKIANS
C
C     J(FNU,Z)*H(FNU+1,1,Z)-J(FNU+1,Z)*H(FNU,1,Z)=2/(PI*I*Z)   Y.GE.0
C
C     J(FNU,Z)*H(FNU+1,2,Z)-J(FNU+1,Z)*H(FNU,2,Z)=-2/(PI*I*Z)  Y.LT.0
C
C     IN THEIR RESPECTIVE HALF PLANES.
C
C     THE PARAMETER MQC CAN HAVE VALUES 1 (THE DEFAULT) FOR A FASTER,
C     LESS DEFINITIVE TEST OR 2 FOR A SLOWER, MORE DEFINITIVE TEST.
C
C     MACHINE CONSTANTS ARE DEFINED IN FUNCTIONS I1MACH, R1MACH, AND
C     D1MACH. THESE MUST BE SELECTED BY THE USER OR SET ACCORDING TO
C     PROLOGUE INSTRUCTIONS.
C
      COMPLEX Z, WR, CJ, CH, CON, T1, T2, CER
      REAL AA, AB, AER, ALIM, DIG, ELIM, EPS, ER, ERTOL, FNU, FNUL,
     * GNU, HPI, R, RL, RM, R1M4, R1M5, R2, T, TOL, XNU, XX, YY,
     * R1MACH, SLAK, FILM, ST, TS, CT, SGN
      INTEGER I, ICASE, IL, IR, IRB, IT, ITL, K, KK, KODE, K1, K2,
     * LFLG, LUN, M, MFLG, N, NU, NZJ, NZH, IERRJ, IERRH, I1MACH, KEPS,
     * KDO, NL, NUL, MQC
      DIMENSION T(20), AER(20), XNU(20), CJ(20), CH(20), KEPS(20),
     * KDO(20)
      DATA LUN /7/
      PARAMETER (MQC=1)
      OPEN(LUN,FILE='CQCBJ.OUT')
C-----------------------------------------------------------------------
C     SET PARAMETERS RELATED TO MACHINE CONSTANTS.
C     TOL IS THE APPROXIMATE UNIT ROUNDOFF LIMITED TO 1.0E-18.
C     ELIM IS THE APPROXIMATE EXPONENTIAL OVER- AND UNDERFLOW LIMIT.
C     EXP(-ELIM).LT.EXP(-ALIM)=EXP(-ELIM)/TOL    AND
C     EXP(ELIM).GT.EXP(ALIM)=EXP(ELIM)*TOL       ARE INTERVALS NEAR
C     UNDERFLOW AND OVERFLOW LIMITS WHERE SCALED ARITHMETIC IS DONE.
C     RL IS THE LOWER BOUNDARY OF THE ASYMPTOTIC EXPANSION FOR LARGE Z.
C     DIG = NUMBER OF BASE 10 DIGITS IN TOL = 10**(-DIG).
C     FNUL IS THE LOWER BOUNDARY OF THE ASYMPTOTIC SERIES FOR LARGE FNU
C-----------------------------------------------------------------------
      R1M4 = R1MACH(4)
      TOL = MAX(R1M4,1.0E-18)
      AA = -LOG10(R1M4)
      K1 = I1MACH(12)
      K2 = I1MACH(13)
      R1M5 = R1MACH(5)
      K = MIN(ABS(K1),ABS(K2))
      ELIM = 2.303E0*(FLOAT(K)*R1M5-3.0E0)
      AB = AA*2.303E0
      ALIM = ELIM + MAX(-AB,-41.45E0)
      DIG = MIN(AA,18.0E0)
      FNUL = 10.0E0 + 6.0E0*(DIG-3.0E0)
      RL = 1.2E0*DIG + 3.0E0
      SLAK = 3.0E0+4.0E0*(-LOG10(TOL)-7.0E0)/11.0E0
      SLAK = MAX(SLAK,3.0E0)
      ERTOL = TOL*10.0E0**SLAK
      RM = 0.5E0*(ALIM + ELIM)
      RM = MIN(RM,200.0E0)
      RM = MAX(RM,RL+10.0E0)
      R2 = MIN(FNUL,RM)
C-----------------------------------------------------------------------
      WRITE (LUN,99999)
99999 FORMAT (' QUICK CHECK ROUTINE FOR THE J BESSEL FUNCTION FROM CBESJ
     *'/)
      WRITE (LUN,99998)
99998 FORMAT (' PARAMETERS TOL,ELIM,ALIM,RL,FNUL,DIG')
      WRITE (LUN,99997) TOL, ELIM, ALIM, RL, FNUL, DIG
99997 FORMAT (6E12.4/)
      ATOL=100.0E0*TOL
      HPI = 2.0E0*ATAN(1.0E0)
      PI = HPI + HPI
      CON=CMPLX(0.0,-1.0/HPI)
      WRITE (LUN,99996) MQC
99996 FORMAT (/' CHECKS IN THE (Z,FNU) SPACE WITH MQC = ',I2/)
C-----------------------------------------------------------------------
C     TEST VALUES OF Z IN -PI.LT.ARG(Z).LE.PI
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C     KDO(K), K=1,IL  DETERMINES WHICH OF THE IL ANGLES IN -PI TO PI
C     ARE USE TO COMPUTE VALUES OF Z
C       KDO(K) = 0  MEANS THAT THE INDEX K WILL BE USED FOR ONE OR TWO
C                   VALUES OF Z, DEPENDING ON THE CHOICE OF KEPS(K)
C              = 1  MEANS THAT THE INDEX K AND THE CORRESPONDING ANGLE
C                   WILL BE SKIPPED
C     KEPS(K), K=1,IL DETERMINES WHICH OF THE ANGLES GET INCREMENTED
C     UP AND DOWN TO PUT VALUES OF Z IN REGIONS WHERE DIFFERENT
C     FORMULAE ARE USED.
C       KEPS(K) =0  MEANS THAT THE ANGLE WILL BE USED WITHOUT CHANGE
C               =1  MEANS THAT THE ANGLE WILL BE INCREMENTED UP AND
C                   DOWN BY EPS
C     THE ANGLES TO BE USED ARE STORED IN THE T(I) ARRAY, I=1,ITL
C-----------------------------------------------------------------------
      IF (MQC.NE.2) THEN
        NL=2
        IL=5
        DO 5 I=1,IL
          KEPS(I)=0
          KDO(I)=0
    5   CONTINUE
        NUL=5
        XNU(1) = 0.0E0
        XNU(2) = 1.0E0
        XNU(3) = 2.0E0
        XNU(4) = 0.5E0*FNUL
        XNU(5) = FNUL + 1.1E0
      ELSE
        NL=4
        IL=13
        DO 6 I=1,IL
          KDO(I)=0
          KEPS(I)=0
    6   CONTINUE
        KDO(2)=1
        KDO(6)=1
        KDO(8)=1
        KDO(12)=1
        KEPS(3)=1
        KEPS(4)=1
        KEPS(5)=1
        KEPS(9)=1
        KEPS(10)=1
        KEPS(11)=1
        NUL=6
        XNU(1) = 0.0E0
        XNU(2) = 0.6E0
        XNU(3) = 1.3E0
        XNU(4) = 2.0E0
        XNU(5) = 0.5E0*FNUL
        XNU(6) = FNUL + 1.1E0
      ENDIF
      I = 2
      EPS = 0.01E0
      FILM=FLOAT(IL-1)
      T(1) = -PI + EPS
      DO 30 K=2,IL
        IF (KDO(K).EQ.0) THEN
          T(I) = PI*FLOAT(-IL+2*K-1)/FILM
          IF (KEPS(K).EQ.0) GO TO 20
          TS=T(I)
          T(I) = TS - EPS
          I = I + 1
          T(I) = TS + EPS
   20     CONTINUE
          I = I + 1
        ENDIF
   30 CONTINUE
      ITL = I - 1
      LFLG = 0
      DO 260 KODE=1,2
        DO 250 N=1,NL
          NP=N+1
          DO 240 NU=1,NUL
            FNU = XNU(NU)
            GNU = FNU + FLOAT(N-1) + 1.0E0
            GNU=SQRT(GNU)
            GNU=MIN(GNU,0.5*RL)
            DO 230 ICASE=1,3
              IRB = MIN(2,ICASE)
              DO 220 IR=IRB,4
                GO TO (50, 60, 70), ICASE
   50           CONTINUE
                R = (GNU*FLOAT(4-IR)+2.0E0*FLOAT(IR-1))/3.0E0
                GO TO 80
   60           CONTINUE
                R = (2.0E0*FLOAT(4-IR)+R2*FLOAT(IR-1))/3.0E0
                GO TO 80
   70           CONTINUE
                IF (R2.GE.RM) GO TO 230
                R = (R2*FLOAT(4-IR)+RM*FLOAT(IR-1))/3.0E0
   80           CONTINUE
                DO 210 IT=1,ITL
                  CT = COS(T(IT))
                  ST = SIN(T(IT))
                  IF (ABS(CT).LT.ATOL) CT = 0.0E0
                  IF (ABS(ST).LT.ATOL) ST = 0.0E0
                  Z = CMPLX(R*CT,R*ST)
                  XX = REAL(Z)
                  YY = AIMAG(Z)
                  IF(XX.EQ.0.0.AND.YY.EQ.0.0) GO TO 210
                  WR=CON/Z
                  M=1
                  IF(YY.LT.0.0) THEN
                    M=2
                    WR=-WR
                  ENDIF
                  CALL CBESJ(Z,FNU,KODE,NP,CJ,NZJ,IERRJ)
                  CALL CBESH(Z,FNU,KODE,M,NP,CH,NZH,IERRH)
                  IF(NZJ.NE.0.OR.NZH.NE.0) GO TO 210
                  IF(IERRJ.NE.0.OR.IERRH.NE.0) GO TO 210
                  IF(KODE.EQ.2) THEN
                    SGN=3.0-2.0*FLOAT(M)
                    WR=WR*CMPLX(COS(XX),-SGN*SIN(XX))
                  ENDIF
                  KK = 0
                  MFLG = 0
                  DO 190 I=1,N
                    T1=CJ(I)*CH(I+1)
                    T2=CJ(I+1)*CH(I)
                    CER=T1-T2-WR
                    ER=CABS(CER)/CABS(WR)
                    IF (ER.GT.ERTOL) THEN
                      IF(MFLG.EQ.0) THEN
                        MFLG = 1
                        KK=I
                      ENDIF
                    ENDIF
                    AER(I)=ER
  190             CONTINUE
                  IF (MFLG.EQ.0) GO TO 210
                  IF (LFLG.EQ.1) GO TO 200
                  WRITE (LUN,99995) ERTOL
99995             FORMAT (/' CASES WHICH VIOLATE THE RELATIVE ERROR TEST
     * WITH ERTOL =', E12.4/)
                  WRITE (LUN,99994)
99994             FORMAT (/' OUTPUT FORMAT'/' KODE,N,IR,IT,NZJ,NZH,ICASE
     *')
                  WRITE (LUN,99993)
99993             FORMAT (' ER(K),K=1,N'/' Z,FNU,CJ(KK),CH(KK), KK=INDEX
     * OF FIRST NON-ZERO CJ,CH PAIR'/)
                  LFLG = 1
  200             CONTINUE
                  WRITE (LUN,99992) KODE, N, IR, IT, NZJ, NZH, ICASE
99992             FORMAT (8I5)
                  WRITE (LUN,99991) (AER(K),K=1,N)
                  WRITE (LUN,99991) Z, FNU, CJ(KK), CH(KK)
99991             FORMAT (9E12.4)
  210           CONTINUE
  220         CONTINUE
  230       CONTINUE
  240     CONTINUE
  250   CONTINUE
  260 CONTINUE
      IF (LFLG.EQ.0) WRITE (LUN,99990)
99990 FORMAT (/' QUICK CHECKS OK'/)
      STOP
      END
      PROGRAM CQCBK
C
C***DATE WRITTEN   830501   (YYMMDD)
C***REVISION DATE  890801, 930101   (YYMMDD)
C
C     CQCBK IS A QUICK CHECK ROUTINE FOR THE COMPLEX K BESSEL FUNCTION
C     GENERATED BY SUBROUTINE CBESK.
C
C     CQCBK GENERATES SEQUENCES OF I AND K BESSEL FUNCTIONS FROM
C     CBESI AND CBESK AND CHECKS THEM AGAINST THE WRONSKIAN EVALUATION
C
C           I(FNU,Z)*K(FNU+1,Z) + I(FNU+1,Z)*K(FNU,Z) = 1/Z
C
C     IN THE RIGHT HALF PLANE AND THE ANALYTIC CONTINUATION FORMULA
C     FOR H(FNU,2,Z) IN TERMS OF THE K FUNCTION
C
C           K(FNU,Z) = C3*H(FNU,2,ZR) + C4*H(FNU,1,ZR)    IM(Z).GE.0
C
C                    = CONJG(K(FNU,CONJG(Z)))             IM(Z).LT.0
C
C     IN THE LEFT HALF PLANE WHERE C3=C1*CONJG(C2)*C5, C4 = C2*C5
C     C1=2*COS(PI*FNU),   C2=EXP(PI*FNU*I/2),   C5 =-PI*I/2   AND
C     ZR = Z*EXP(-3*PI*I/2) = Z*I
C
C     THE PARAMETER MQC CAN HAVE VALUES 1 (THE DEFAULT) FOR A FASTER,
C     LESS DEFINITIVE TEST OR 2 FOR A SLOWER, MORE DEFINITIVE TEST.
C
C     MACHINE CONSTANTS ARE DEFINED IN FUNCTIONS I1MACH, R1MACH, AND
C     D1MACH. THESE MUST BE SELECTED BY THE USER OR SET ACCORDING TO
C     PROLOGUE INSTRUCTIONS.
C
      COMPLEX CONE, CSGN, CV, CW, CY, C1, C2, C3, C4, V, W, Y, Z, ZR,
     * ZZ, CIP, COE
      REAL AA, AB, AER, ALIM, ARG, ATOL, DIG, ELIM, EPS, ER,
     * ERTOL, FFNU, FNU, FNUL, HPI, PI, R, RL, RM, R1M4, R1M5, R2,
     * T, TOL, XNU, XX, R1MACH, FILM, ST, CT, TS, SLAK, YY
      INTEGER I, ICASE, IFNU, IL, IR, IRB, IT, ITL, I4, K, KK, KODE,
     * K1, K2, LFLG, LUN, M, MFLG, N, NU, NZ, N1, IERR, I1MACH,
     * KEPS, KDO, MQC, NL, NUL
      DIMENSION T(20), AER(20), XNU(20), V(20), Y(20), W(20), KEPS(20),
     * KDO(20), CIP(4)
      DATA LUN /7/
      DATA CIP(1),CIP(2),CIP(3),CIP(4) / (1.0E0,0.0E0),(0.0E0,1.0E0),
     * (-1.0E0,0.0E0),(0.0E0,-1.0E0) /
      PARAMETER (MQC=1)
      OPEN(LUN,FILE='CQCBK.OUT')
C-----------------------------------------------------------------------
C     SET PARAMETERS RELATED TO MACHINE CONSTANTS.
C     TOL IS THE APPROXIMATE UNIT ROUNDOFF LIMITED TO 1.0E-18.
C     ELIM IS THE APPROXIMATE EXPONENTIAL OVER- AND UNDERFLOW LIMIT.
C     EXP(-ELIM).LT.EXP(-ALIM)=EXP(-ELIM)/TOL    AND
C     EXP(ELIM).GT.EXP(ALIM)=EXP(ELIM)*TOL       ARE INTERVALS NEAR
C     UNDERFLOW AND OVERFLOW LIMITS WHERE SCALED ARITHMETIC IS DONE.
C     RL IS THE LOWER BOUNDARY OF THE ASYMPTOTIC EXPANSION FOR LARGE Z.
C     DIG = NUMBER OF BASE 10 DIGITS IN TOL = 10**(-DIG).
C     FNUL IS THE LOWER BOUNDARY OF THE ASYMPTOTIC SERIES FOR LARGE FNU
C-----------------------------------------------------------------------
      R1M4 = R1MACH(4)
      TOL = MAX(R1M4,1.0E-18)
      AA = -ALOG10(R1M4)
      K1 = I1MACH(12)
      K2 = I1MACH(13)
      R1M5 = R1MACH(5)
      K = MIN(ABS(K1),ABS(K2))
      ELIM = 2.303E0*(FLOAT(K)*R1M5-3.0E0)
      AB = AA*2.303E0
      ALIM = ELIM + MAX(-AB,-41.45E0)
      DIG = MIN(AA,18.0E0)
      FNUL = 10.0E0 + 6.0E0*(DIG-3.0E0)
      RL = 1.2E0*DIG + 3.0E0
      SLAK = 3.0E0+4.0E0*(-ALOG10(TOL)-7.0E0)/11.0E0
      SLAK = MAX(SLAK,3.0E0)
      ERTOL = TOL*10.0E0**SLAK
      RM = 0.5E0*(ALIM + ELIM)
      RM = MIN(RM,200.0E0)
      RM = MAX(RM,RL+10.0E0)
      R2 = MIN(FNUL,RM)
C-----------------------------------------------------------------------
      WRITE (LUN,99999)
99999 FORMAT (' QUICK CHECK ROUTINE FOR THE K BESSEL FUNCTION FROM CBESK
     *'/)
      WRITE (LUN,99998)
99998 FORMAT (' PARAMETERS TOL,ELIM,ALIM,RL,FNUL,DIG')
      WRITE (LUN,99997) TOL, ELIM, ALIM, RL, FNUL, DIG
99997 FORMAT (6E12.4/)
      CONE = CMPLX(1.0E0,0.0E0)
      ATOL = 100.0E0*TOL
      HPI = 2.0E0*ATAN(1.0E0)
      PI = HPI + HPI
      WRITE (LUN,99996) MQC
99996 FORMAT (/' CHECKS IN THE (Z,FNU) SPACE WITH MQC = ',I2/)
C-----------------------------------------------------------------------
C     TEST VALUES OF Z IN -PI.LT.ARG(Z).LE.PI NEAR FORMULA BOUNDARIES
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C     KDO(K), K=1,IL  DETERMINES WHICH OF THE IL ANGLES IN -PI TO PI
C     ARE USE TO COMPUTE VALUES OF Z
C       KDO(K) = 0  MEANS THAT THE INDEX K WILL BE USED FOR ONE OR TWO
C                   VALUES OF Z, DEPENDING ON THE CHOICE OF KEPS(K)
C              = 1  MEANS THAT THE INDEX K AND THE CORRESPONDING ANGLE
C                   WILL BE SKIPPED
C     KEPS(K), K=1,IL DETERMINES WHICH OF THE ANGLES GET INCREMENTED
C     UP AND DOWN TO PUT VALUES OF Z IN REGIONS WHERE DIFFERENT
C     FORMULAE ARE USED.
C       KEPS(K) =0  MEANS THAT THE ANGLE WILL BE USED WITHOUT CHANGE
C               =1  MEANS THAT THE ANGLE WILL BE INCREMENTED UP AND
C                   DOWN BY EPS
C     THE ANGLES TO BE USED ARE STORED IN THE T(I) ARRAY, I=1,ITL
C-----------------------------------------------------------------------
      IF (MQC.NE.2) THEN
        NL=2
        IL=5
        DO 5 I=1,IL
          KEPS(I)=0
          KDO(I)=0
    5   CONTINUE
        NUL=5
        XNU(1) = 0.0E0
        XNU(2) = 1.0E0
        XNU(3) = 2.0E0
        XNU(4) = 0.5E0*FNUL
        XNU(5) = FNUL + 1.1E0
      ELSE
        NL=4
        IL=13
        DO 6 I=1,IL
          KDO(I)=0
          KEPS(I)=0
    6   CONTINUE
        KDO(2)=1
        KDO(6)=1
        KDO(8)=1
        KDO(12)=1
        KEPS(3)=1
        KEPS(4)=1
        KEPS(5)=1
        KEPS(9)=1
        KEPS(10)=1
        KEPS(11)=1
        NUL=6
        XNU(1) = 0.0E0
        XNU(2) = 0.6E0
        XNU(3) = 1.3E0
        XNU(4) = 2.0E0
        XNU(5) = 0.5E0*FNUL
        XNU(6) = FNUL + 1.1E0
      ENDIF
      I = 2
      EPS = 0.01E0
      FILM=FLOAT(IL-1)
      T(1) = -PI + EPS
      DO 30 K=2,IL
        IF (KDO(K).EQ.0) THEN
          T(I) = PI*FLOAT(-IL+2*K-1)/FILM
          IF (KEPS(K).EQ.0) GO TO 20
          TS=T(I)
          T(I) = TS - EPS
          I = I + 1
          T(I) = TS + EPS
   20     CONTINUE
          I = I + 1
        ENDIF
   30 CONTINUE
      ITL = I - 1
      LFLG = 0
      DO 200 KODE=1,2
        DO 190 N=1,NL
          N1 = N + 1
          DO 180 NU=1,NUL
            FNU = XNU(NU)
            IFNU = INT(FNU)
            FFNU = FNU - FLOAT(IFNU)
            ARG = HPI*FFNU
            CSGN = CMPLX(COS(ARG),SIN(ARG))
            I4 = MOD(IFNU,4)+1
            CSGN = CSGN*CIP(I4)
            DO 170 ICASE=1,3
              IRB = MIN(2,ICASE)
              DO 160 IR=IRB,4
                GO TO (50, 60, 70), ICASE
   50           CONTINUE
                R = (0.2E0*FLOAT(4-IR)+2.0E0*FLOAT(IR-1))/3.0E0
                GO TO 80
   60           CONTINUE
                R = (2.0E0*FLOAT(4-IR)+R2*FLOAT(IR-1))/3.0E0
                GO TO 80
   70           CONTINUE
                IF (R2.GE.RM) GO TO 170
                R = (R2*FLOAT(4-IR)+RM*FLOAT(IR-1))/3.0E0
   80           CONTINUE
                DO 150 IT=1,ITL
                  CT = COS(T(IT))
                  ST = SIN(T(IT))
                  IF (ABS(CT).LT.ATOL) CT = 0.0E0
                  IF (ABS(ST).LT.ATOL) ST = 0.0E0
                  XX = R*CT
                  YY = R*ST
                  Z = CMPLX(XX,YY)
                  IF (XX.GE.0.0E0) THEN
C-----------------------------------------------------------------------
C     WRONSKIAN CHECKS IN THE RIGHT HALF PLANE
C-----------------------------------------------------------------------
                    CALL CBESI(Z, FNU, KODE, N1, W, NZ, IERR)
                    IF (NZ.NE.0.OR.IERR.NE.0) GO TO 150
                    CALL CBESK(Z, FNU, KODE, N1, Y, NZ, IERR)
                    IF (NZ.NE.0.OR.IERR.NE.0) GO TO 150
C-----------------------------------------------------------------------
C     ADJUSTMENTS TO WRONSKIAN DUE TO SCALING OF I AND K FUNCTIONS
C     ON KODE=2
C-----------------------------------------------------------------------
                    CV = CONE/Z
                    IF (KODE.EQ.2) THEN
                      CV = CV*CMPLX(COS(YY),SIN(YY))
                    ENDIF
                    MFLG = 0
                    KK=0
                    DO 130 I=1,N
                      CW = W(I)*Y(I+1)
                      CY = W(I+1)*Y(I)
                      CY = CY + CW - CV
                      ER = CABS(CY)/CABS(CV)
                      AER(I) = ER
                      IF (ER.GT.ERTOL) THEN
                        IF(KK.EQ.0) THEN
                          MFLG = 1
                          KK=I
                        ENDIF
                      ENDIF
  130               CONTINUE
                  ELSE
C-----------------------------------------------------------------------
C     ANALYTIC CONTINUATION FORMULA CHECKS FOR LEFT HALF PLANE IN TERMS
C     OF H(FNU,1,Z) AND H(FNU,2,Z)
C-----------------------------------------------------------------------
                    ZZ = Z
                    IF (YY.LT.0.0E0) THEN
                      ZZ=CONJG(Z)
                    ENDIF
                    ZR = ZZ*CMPLX(0.0E0,1.0E0)
                    M=1
                    CALL CBESH(ZR, FNU, KODE, M, N, W, NZ, IERR)
                    IF (IERR.NE.0) GO TO 150
                    M=2
                    CALL CBESH(ZR, FNU, KODE, M, N, V, NZ, IERR)
                    IF (IERR.NE.0) GO TO 150
                    CALL CBESK(Z, FNU, KODE, N, Y, NZ, IERR)
                    IF (NZ.NE.0.OR.IERR.NE.0) GO TO 150
                    COE=CMPLX(0.0E0,-HPI)
                    MFLG = 0
                    KK = 0
                    AA = 2.0E0*COS(PI*FFNU)
                    IF(MOD(IFNU,2).NE.0) AA = -AA
                    C1 = CMPLX(AA,0.0E0)
                    C2 = CSGN
                    DO 135 I=1,N
                      C3 = C1
                      C4 = C2
                      IF (KODE.EQ.2) THEN
C-----------------------------------------------------------------------
C     ADJUSTMENTS TO COEFICIENTS DUE TO SCALING OF H(FNU,1,Z) AND
C     H(FNU,2,Z) FUNCTIONS ON KODE = 2.
C-----------------------------------------------------------------------
                        AB=CABS(V(I))
                        AA=ALOG(AB)+XX+XX
                        IF (AA.GT.ELIM) GO TO 150
                        IF (AA.LT.-ELIM) THEN
                           C3 = CMPLX(0.0E0,0.0E0)
                        ELSE
                          CW = ZZ+ZZ
                          C3 = C3*CEXP(CW)
                        ENDIF
                      ENDIF
                      CY = (C3*CONJG(C2)*V(I)+C4*W(I))*COE
                      IF (YY.LT.0.0E0) THEN
                        CY=CONJG(CY)
                      ENDIF
                      ER = CABS(CY-Y(I))/CABS(Y(I))
                      AER(I) = ER
                      IF (ER.GT.ERTOL) THEN
                        IF(KK.EQ.0) THEN
                          MFLG = 1
                          KK=I
                        ENDIF
                      ENDIF
                      C2 = C2*CMPLX(0.0E0,1.0E0)
                      C1 = -C1
  135               CONTINUE
                  ENDIF
                  IF (MFLG.EQ.0) GO TO 150
                  IF (LFLG.EQ.1) GO TO 140
                  WRITE (LUN,99995) ERTOL
99995             FORMAT (/' CASES WHICH VIOLATE THE RELATIVE ERROR TEST
     * WITH ERTOL =', E12.4/)
                  WRITE (LUN,99994)
99994             FORMAT (/' OUTPUT FORMAT'/' KODE,N,IR,IT,ICASE,KK')
                  WRITE (LUN,99993)
99993             FORMAT (' ER(K),K=1,N'/' Z,FNU,Y(KK)        KK=INDEX O
     *F FIRST NON-ZERO PAIR'/)
                  LFLG = 1
  140             CONTINUE
                  WRITE (LUN,99992) KODE, N, IR, IT, ICASE, KK
99992             FORMAT (6I5)
                  WRITE (LUN,99991) (AER(K),K=1,N)
                  WRITE (LUN,99991) Z, FNU, Y(KK)
99991             FORMAT (6E12.4)
  150           CONTINUE
  160         CONTINUE
  170       CONTINUE
  180     CONTINUE
  190   CONTINUE
  200 CONTINUE
      IF (LFLG.EQ.0) WRITE (LUN,99990)
99990 FORMAT (/' QUICK CHECKS OK'/)
      STOP
      END
      PROGRAM CQCBY
C
C***DATE WRITTEN   830501   (YYMMDD)
C***REVISION DATE  890801, 930101   (YYMMDD)
C
C     CQCBY IS A QUICK CHECK ROUTINE FOR THE COMPLEX Y BESSEL FUNCTION
C     GENERATED BY SUBROUTINE CBESY.
C
C     CQCBY GENERATES SEQUENCES OF Y BESSEL FUNCTIONS FROM CBESY AND
C     CBESYH AND COMPARES THEM FOR A VARIETY OF VALUES IN THE (Z,FNU)
C     SPACE. CBESYH IS AN OLD VERSION OF CBESY WHICH COMPUTES THE Y
C     FUNCTION FROM THE H FUNCTIONS OF KINDS 1 AND 2.
C
C     THE PARAMETER MQC CAN HAVE VALUES 1 (THE DEFAULT) FOR A FASTER,
C     LESS DEFINITIVE TEST OR 2 FOR A SLOWER, MORE DEFINITIVE TEST.
C
C     MACHINE CONSTANTS ARE DEFINED IN FUNCTIONS I1MACH, R1MACH, AND
C     D1MACH. THESE MUST BE SELECTED BY THE USER OR SET ACCORDING TO
C     PROLOGUE INSTRUCTIONS.
C
      COMPLEX  CW, CWRK, V, W, Z
      REAL AA, AB, AER, ALIM, ATOL, AV, CABS, DIG, ELIM, EPS, ER,
     * ERTOL, FFNU, FNU, FNUL, PI, R, RL, RM, R1M4, R1M5, R2,
     * T, TOL, XNU, XX, YY, R1MACH, CT, ST, TS, SLAK, FILM
      INTEGER I, ICASE, IFNU, IL, IR, IRB, IT, ITL, K, KK,
     * KODE, K1, K2, LFLG, LUN, MFLG, N, NU, NZ1, NZ2, IERR, I1MACH,
     * KEPS, KDO, NL, NUL, MQC
      DIMENSION T(20), AER(20), XNU(20), W(20), V(20),
     * CWRK(20), KEPS(20), KDO(20)
      DATA LUN /6/
      PARAMETER (MQC=1)
      OPEN(LUN,FILE='CQCBY.OUT')
C-----------------------------------------------------------------------
C     SET PARAMETERS RELATED TO MACHINE CONSTANTS.
C     TOL IS THE APPROXIMATE UNIT ROUNDOFF LIMITED TO 1.0E-18.
C     ELIM IS THE APPROXIMATE EXPONENTIAL OVER- AND UNDERFLOW LIMIT.
C     EXP(-ELIM).LT.EXP(-ALIM)=EXP(-ELIM)/TOL    AND
C     EXP(ELIM).GT.EXP(ALIM)=EXP(ELIM)*TOL       ARE INTERVALS NEAR
C     UNDERFLOW AND OVERFLOW LIMITS WHERE SCALED ARITHMETIC IS DONE.
C     RL IS THE LOWER BOUNDARY OF THE ASYMPTOTIC EXPANSION FOR LARGE Z.
C     DIG = NUMBER OF BASE 10 DIGITS IN TOL = 10**(-DIG).
C     FNUL IS THE LOWER BOUNDARY OF THE ASYMPTOTIC SERIES FOR LARGE FNU
C-----------------------------------------------------------------------
      R1M4 = R1MACH(4)
      TOL = AMAX1(R1M4,1.0E-18)
      AA = -ALOG10(R1M4)
      K1 = I1MACH(12)
      K2 = I1MACH(13)
      R1M5 = R1MACH(5)
      K = MIN0(IABS(K1),IABS(K2))
      ELIM = 2.303E0*(FLOAT(K)*R1M5-3.0E0)
      AB = AA*2.303E0
      ALIM = ELIM + AMAX1(-AB,-41.45E0)
      DIG = AMIN1(AA,18.0E0)
      FNUL = 10.0E0 + 6.0E0*(DIG-3.0E0)
      RL = 1.2E0*DIG + 3.0E0
      SLAK = 3.0E0+4.0E0*(-ALOG10(TOL)-7.0E0)/11.0E0
      SLAK = AMAX1(SLAK,3.0E0)
      ERTOL = TOL*10.0E0**SLAK
      RM = 0.5E0*(ALIM + ELIM)
      RM = AMIN1(RM,200.0E0)
      RM = AMAX1(RM,RL+10.0E0)
      R2 = AMIN1(FNUL,RM)
C-----------------------------------------------------------------------
      WRITE (LUN,99999)
99999 FORMAT (' QUICK CHECK ROUTINE FOR THE Y BESSEL FUNCTION FROM CBESY
     * '/)
      WRITE (LUN,99998)
99998 FORMAT (' PARAMETERS TOL,ELIM,ALIM,RL,FNUL,DIG')
      WRITE (LUN,99997) TOL, ELIM, ALIM, RL, FNUL, DIG
99997 FORMAT (6E12.4/)
      ATOL = 100.0E0*TOL
      PI = 4.0E0*ATAN(1.0E0)
      WRITE (LUN,99996) MQC
99996 FORMAT (/' CHECKS IN THE (Z,FNU) SPACE WITH MQC = ',I2/)
C-----------------------------------------------------------------------
C     TEST VALUES OF Z IN -PI/2.LT.ARG(Z).LE.PI
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C     KDO(K), K=1,IL  DETERMINES WHICH OF THE IL ANGLES IN -PI TO PI
C     ARE USE TO COMPUTE VALUES OF Z
C       KDO(K) = 0  MEANS THAT THE INDEX K WILL BE USED FOR ONE OR TWO
C                   VALUES OF Z, DEPENDING ON THE CHOICE OF KEPS(K)
C              = 1  MEANS THAT THE INDEX K AND THE CORRESPONDING ANGLE
C                   WILL BE SKIPPED
C     KEPS(K), K=1,IL DETERMINES WHICH OF THE ANGLES GET INCREMENTED
C     UP AND DOWN TO PUT VALUES OF Z IN REGIONS WHERE DIFFERENT
C     FORMULAE ARE USED.
C       KEPS(K) =0  MEANS THAT THE ANGLE WILL BE USED WITHOUT CHANGE
C               =1  MEANS THAT THE ANGLE WILL BE INCREMENTED UP AND
C                   DOWN BY EPS
C     THE ANGLES TO BE USED ARE STORED IN THE T(I) ARRAY, I=1,ITL
C-----------------------------------------------------------------------
      IF (MQC.NE.2) THEN
        NL=2
        IL=5
        DO 5 I=1,IL
          KEPS(I)=0
          KDO(I)=0
    5   CONTINUE
        NUL=5
        XNU(1) = 0.0E0
        XNU(2) = 1.0E0
        XNU(3) = 2.0E0
        XNU(4) = 0.5E0*FNUL
        XNU(5) = FNUL + 1.1E0
      ELSE
        NL=4
        IL=13
        DO 6 I=1,IL
          KDO(I)=0
          KEPS(I)=0
    6   CONTINUE
        KDO(2)=1
        KDO(6)=1
        KDO(8)=1
        KDO(12)=1
        KEPS(3)=1
        KEPS(4)=1
        KEPS(5)=1
        KEPS(9)=1
        KEPS(10)=1
        KEPS(11)=1
        NUL=6
        XNU(1) = 0.0E0
        XNU(2) = 0.6E0
        XNU(3) = 1.3E0
        XNU(4) = 2.0E0
        XNU(5) = 0.5E0*FNUL
        XNU(6) = FNUL + 1.1E0
      ENDIF
      I = 2
      EPS = 0.01E0
      FILM=FLOAT(IL-1)
      T(1) = -PI + EPS
      DO 30 K=2,IL
        IF (KDO(K).EQ.0) THEN
          T(I) = PI*FLOAT(-IL+2*K-1)/FILM
          IF (KEPS(K).EQ.0) GO TO 20
          TS=T(I)
          T(I) = TS - EPS
          I = I + 1
          T(I) = TS + EPS
   20     CONTINUE
          I = I + 1
        ENDIF
   30 CONTINUE
      ITL = I - 1
      LFLG = 0
      DO 190 KODE=1,2
        DO 180 N=1,NL
          DO 170 NU=1,NUL
            FNU = XNU(NU)
            IFNU = INT(FNU)
            FFNU = FNU-FLOAT(IFNU)
            DO 160 ICASE=1,3
              IRB = MIN0(2,ICASE)
              DO 150 IR=IRB,4
                GO TO (50, 60, 70), ICASE
   50           CONTINUE
                R = (EPS*FLOAT(4-IR)+2.0E0*FLOAT(IR-1))/3.0E0
                GO TO 80
   60           CONTINUE
                R = (2.0E0*FLOAT(4-IR)+R2*FLOAT(IR-1))/3.0E0
                GO TO 80
   70           CONTINUE
                IF (R2.EQ.RM) GO TO 160
                R = (R2*FLOAT(4-IR)+RM*FLOAT(IR-1))/3.0E0
   80           CONTINUE
                DO 140 IT=1,ITL
                  CT = COS(T(IT))
                  ST = SIN(T(IT))
                  IF (ABS(CT).LT.ATOL) CT = 0.0E0
                  IF (ABS(ST).LT.ATOL) ST = 0.0E0
                  Z = CMPLX(R*CT,R*ST)
                  XX = REAL(Z)
                  YY = AIMAG(Z)
                  CALL CBESY(Z, FNU, KODE, N, V, NZ2, CWRK, IERR)
                  IF (NZ2.NE.0.OR.IERR.NE.0) GO TO 140
                  CALL CBESYH(Z, FNU, KODE, N, W, NZ1, CWRK, IERR)
                  IF (NZ1.NE.0.OR.IERR.NE.0) GO TO 140
                  MFLG = 0
                  DO 120 I=1,N
                    AB = FNU+FLOAT(I-1)
                    AA = MAX(0.5E0,AB)
                    CW = W(I) - V(I)
                    AV = CABS(V(I))
                    ER = CABS(CW)
                    IF (AV.NE.0.0E0) THEN
                      IF (YY.EQ.0.0E0) THEN
                        IF (XX.GT.0.0E0) THEN
                          IF (ABS(XX).LT.AA) ER = ER/AV
                        ELSE
                          IF (ABS(FFNU-0.5E0).LT.0.125E0) THEN
                            IF (ABS(XX).LT.AA) ER = ER/AV
                          ELSE
                            ER = ER/AV
                          ENDIF
                        ENDIF
                      ELSE
                        ER = ER/AV
                      ENDIF
                    ENDIF
                    AER(I) = ER
                    IF (ER.GT.ERTOL) MFLG = 1
  120             CONTINUE
                  IF (MFLG.EQ.0) GO TO 140
                  IF (LFLG.EQ.1) GO TO 130
                  WRITE (LUN,99995) ERTOL
99995             FORMAT (/' CASES WHICH VIOLATE THE RELATIVE ERROR TEST
     * WITH ERTOL = ', E12.4/)
                  WRITE (LUN,99994)
99994             FORMAT (/' OUTPUT FORMAT'/' KODE,N,IR,IT,NZ1,NZ2,ICASE
     *')
                  WRITE (LUN,99993)
99993             FORMAT (' ER(K),K=1,N'/' Z,FNU,W(KK),V(KK), KK=INDEX O
     *F FIRST NON-ZERO W,V PAIR'/)
                  LFLG = 1
  130             CONTINUE
                  KK = MAX0(NZ1,NZ2) + 1
                  KK = MIN0(N,KK)
                  WRITE (LUN,99992) KODE, N, IR, IT, NZ1, NZ2, ICASE
99992             FORMAT (8I5)
                  WRITE (LUN,99991) (AER(K),K=1,N)
                  WRITE (LUN,99991) Z, FNU, W(KK), V(KK)
99991             FORMAT (7E12.4)
  140           CONTINUE
  150         CONTINUE
  160       CONTINUE
  170     CONTINUE
  180   CONTINUE
  190 CONTINUE
      IF (LFLG.EQ.0) WRITE (LUN,99990)
99990 FORMAT (/' QUICK CHECKS OK'/)
      STOP
      END
      SUBROUTINE CBESYH(Z, FNU, KODE, N, CY, NZ, CWRK, IERR)
C***BEGIN PROLOGUE  CBESYH
C***DATE WRITTEN   830501   (YYMMDD)
C***REVISION DATE  890801, 930101   (YYMMDD)
C***CATEGORY NO.  B5K
C***KEYWORDS  Y-BESSEL FUNCTION,BESSEL FUNCTION OF COMPLEX ARGUMENT,
C             BESSEL FUNCTION OF SECOND KIND
C***AUTHOR  AMOS, DONALD E., SANDIA NATIONAL LABORATORIES
C***PURPOSE  TO COMPUTE THE Y-BESSEL FUNCTION OF A COMPLEX ARGUMENT
C***DESCRIPTION
C
C         ON KODE=1, CBESYH COMPUTES AN N MEMBER SEQUENCE OF COMPLEX
C         BESSEL FUNCTIONS CY(I)=Y(FNU+I-1,Z) FOR REAL, NONNEGATIVE
C         ORDERS FNU+I-1, I=1,...,N AND COMPLEX Z IN THE CUT PLANE
C         -PI.LT.ARG(Z).LE.PI. ON KODE=2, CBESYH RETURNS THE SCALED
C         FUNCTIONS
C
C         CY(I)=EXP(-ABS(Y))*Y(FNU+I-1,Z)   I = 1,...,N , Y=AIMAG(Z)
C
C         WHICH REMOVE THE EXPONENTIAL GROWTH IN BOTH THE UPPER AND
C         LOWER HALF PLANES FOR Z TO INFINITY. DEFINITIONS AND NOTATION
C         ARE FOUND IN THE NBS HANDBOOK OF MATHEMATICAL FUNCTIONS
C         (REF. 1).
C
C         INPUT
C           Z      - Z=CMPLX(X,Y), Z.NE.CMPLX(0.,0.),-PI.LT.ARG(Z).LE.PI
C           FNU    - ORDER OF INITIAL Y FUNCTION, FNU.GE.0.0E0
C           KODE   - A PARAMETER TO INDICATE THE SCALING OPTION
C                    KODE= 1  RETURNS
C                             CY(I)=Y(FNU+I-1,Z), I=1,...,N
C                        = 2  RETURNS
C                             CY(I)=Y(FNU+I-1,Z)*EXP(-ABS(Y)), I=1,...,N
C                             WHERE Y=AIMAG(Z)
C           N      - NUMBER OF MEMBERS OF THE SEQUENCE, N.GE.1
C           CWRK   - A COMPLEX WORK VECTOR OF DIMENSION AT LEAST N
C
C         OUTPUT
C           CY     - A COMPLEX VECTOR WHOSE FIRST N COMPONENTS CONTAIN
C                    VALUES FOR THE SEQUENCE
C                    CY(I)=Y(FNU+I-1,Z)  OR
C                    CY(I)=Y(FNU+I-1,Z)*EXP(-ABS(Y))  I=1,...,N
C                    DEPENDING ON KODE.
C           NZ     - NZ=0 , A NORMAL RETURN
C                    NZ.GT.0 , NZ COMPONENTS OF CY SET TO ZERO DUE TO
C                    UNDERFLOW (GENERALLY ON KODE=2)
C           IERR   - ERROR FLAG
C                    IERR=0, NORMAL RETURN - COMPUTATION COMPLETED
C                    IERR=1, INPUT ERROR   - NO COMPUTATION
C                    IERR=2, OVERFLOW      - NO COMPUTATION, FNU+N-1 IS
C                            TOO LARGE OR CABS(Z) IS TOO SMALL OR BOTH
C                    IERR=3, CABS(Z) OR FNU+N-1 LARGE - COMPUTATION DONE
C                            BUT LOSSES OF SIGNIFCANCE BY ARGUMENT
C                            REDUCTION PRODUCE LESS THAN HALF OF MACHINE
C                            ACCURACY
C                    IERR=4, CABS(Z) OR FNU+N-1 TOO LARGE - NO COMPUTA-
C                            TION BECAUSE OF COMPLETE LOSSES OF SIGNIFI-
C                            CANCE BY ARGUMENT REDUCTION
C                    IERR=5, ERROR              - NO COMPUTATION,
C                            ALGORITHM TERMINATION CONDITION NOT MET
C
C***LONG DESCRIPTION
C
C         THE COMPUTATION IS CARRIED OUT BY THE FORMULA
C
C              Y(FNU,Z) = 0.5*(H(1,FNU,Z) - H(2,FNU,Z))/I
C
C         WHERE I**2 = -1 AND THE HANKEL BESSEL FUNCTIONS H(1,FNU,Z)
C         AND H(2,FNU,Z) ARE CALCULATED IN CBESH.
C
C         FOR NEGATIVE ORDERS,THE FORMULA
C
C              Y(-FNU,Z) = Y(FNU,Z)*COS(PI*FNU) + J(FNU,Z)*SIN(PI*FNU)
C
C         CAN BE USED. HOWEVER,FOR LARGE ORDERS CLOSE TO HALF ODD
C         INTEGERS THE FUNCTION CHANGES RADICALLY. WHEN FNU IS A LARGE
C         POSITIVE HALF ODD INTEGER,THE MAGNITUDE OF Y(-FNU,Z)=J(FNU,Z)*
C         SIN(PI*FNU) IS A LARGE NEGATIVE POWER OF TEN. BUT WHEN FNU IS
C         NOT A HALF ODD INTEGER, Y(FNU,Z) DOMINATES IN MAGNITUDE WITH A
C         LARGE POSITIVE POWER OF TEN AND THE MOST THAT THE SECOND TERM
C         CAN BE REDUCED IS BY UNIT ROUNDOFF FROM THE COEFFICIENT. THUS,
C         WIDE CHANGES CAN OCCUR WITHIN UNIT ROUNDOFF OF A LARGE HALF
C         ODD INTEGER. HERE, LARGE MEANS FNU.GT.CABS(Z).
C
C         IN MOST COMPLEX VARIABLE COMPUTATION, ONE MUST EVALUATE ELE-
C         MENTARY FUNCTIONS. WHEN THE MAGNITUDE OF Z OR FNU+N-1 IS
C         LARGE, LOSSES OF SIGNIFICANCE BY ARGUMENT REDUCTION OCCUR.
C         CONSEQUENTLY, IF EITHER ONE EXCEEDS U1=SQRT(0.5/UR), THEN
C         LOSSES EXCEEDING HALF PRECISION ARE LIKELY AND AN ERROR FLAG
C         IERR=3 IS TRIGGERED WHERE UR=R1MACH(4)=UNIT ROUNDOFF. ALSO
C         IF EITHER IS LARGER THAN U2=0.5/UR, THEN ALL SIGNIFICANCE IS
C         LOST AND IERR=4. IN ORDER TO USE THE INT FUNCTION, ARGUMENTS
C         MUST BE FURTHER RESTRICTED NOT TO EXCEED THE LARGEST MACHINE
C         INTEGER, U3=I1MACH(9). THUS, THE MAGNITUDE OF Z AND FNU+N-1 IS
C         RESTRICTED BY MIN(U2,U3). ON 32 BIT MACHINES, U1,U2, AND U3
C         ARE APPROXIMATELY 2.0E+3, 4.2E+6, 2.1E+9 IN SINGLE PRECISION
C         ARITHMETIC AND 1.3E+8, 1.8E+16, 2.1E+9 IN DOUBLE PRECISION
C         ARITHMETIC RESPECTIVELY. THIS MAKES U2 AND U3 LIMITING IN
C         THEIR RESPECTIVE ARITHMETICS. THIS MEANS THAT ONE CAN EXPECT
C         TO RETAIN, IN THE WORST CASES ON 32 BIT MACHINES, NO DIGITS
C         IN SINGLE AND ONLY 7 DIGITS IN DOUBLE PRECISION ARITHMETIC.
C         SIMILAR CONSIDERATIONS HOLD FOR OTHER MACHINES.
C
C         THE APPROXIMATE RELATIVE ERROR IN THE MAGNITUDE OF A COMPLEX
C         BESSEL FUNCTION CAN BE EXPRESSED BY P*10**S WHERE P=MAX(UNIT
C         ROUNDOFF,1.0E-18) IS THE NOMINAL PRECISION AND 10**S REPRE-
C         SENTS THE INCREASE IN ERROR DUE TO ARGUMENT REDUCTION IN THE
C         ELEMENTARY FUNCTIONS. HERE, S=MAX(1,ABS(LOG10(CABS(Z))),
C         ABS(LOG10(FNU))) APPROXIMATELY (I.E. S=MAX(1,ABS(EXPONENT OF
C         CABS(Z),ABS(EXPONENT OF FNU)) ). HOWEVER, THE PHASE ANGLE MAY
C         HAVE ONLY ABSOLUTE ACCURACY. THIS IS MOST LIKELY TO OCCUR WHEN
C         ONE COMPONENT (IN ABSOLUTE VALUE) IS LARGER THAN THE OTHER BY
C         SEVERAL ORDERS OF MAGNITUDE. IF ONE COMPONENT IS 10**K LARGER
C         THAN THE OTHER, THEN ONE CAN EXPECT ONLY MAX(ABS(LOG10(P))-K,
C         0) SIGNIFICANT DIGITS; OR, STATED ANOTHER WAY, WHEN K EXCEEDS
C         THE EXPONENT OF P, NO SIGNIFICANT DIGITS REMAIN IN THE SMALLER
C         COMPONENT. HOWEVER, THE PHASE ANGLE RETAINS ABSOLUTE ACCURACY
C         BECAUSE, IN COMPLEX ARITHMETIC WITH PRECISION P, THE SMALLER
C         COMPONENT WILL NOT (AS A RULE) DECREASE BELOW P TIMES THE
C         MAGNITUDE OF THE LARGER COMPONENT. IN THESE EXTREME CASES,
C         THE PRINCIPAL PHASE ANGLE IS ON THE ORDER OF +P, -P, PI/2-P,
C         OR -PI/2+P.
C
C***REFERENCES  HANDBOOK OF MATHEMATICAL FUNCTIONS BY M. ABRAMOWITZ
C                 AND I. A. STEGUN, NBS AMS SERIES 55, U.S. DEPT. OF
C                 COMMERCE, 1955.
C
C               COMPUTATION OF BESSEL FUNCTIONS OF COMPLEX ARGUMENT
C                 BY D. E. AMOS, SAND83-0083, MAY, 1983.
C
C               COMPUTATION OF BESSEL FUNCTIONS OF COMPLEX ARGUMENT
C                 AND LARGE ORDER BY D. E. AMOS, SAND83-0643, MAY, 1983
C
C               A SUBROUTINE PACKAGE FOR BESSEL FUNCTIONS OF A COMPLEX
C                 ARGUMENT AND NONNEGATIVE ORDER BY D. E. AMOS, SAND85-
C                 1018, MAY, 1985
C
C               A PORTABLE PACKAGE FOR BESSEL FUNCTIONS OF A COMPLEX
C                 ARGUMENT AND NONNEGATIVE ORDER BY D. E. AMOS, TRANS.
C                 MATH. SOFTWARE, 12, NO. 3, SEPTEMBER 1986, PP 265-273.
C
C***ROUTINES CALLED  CBESH,I1MACH,R1MACH
C***END PROLOGUE  CBESYH
C
      COMPLEX CWRK, CY, C1, C2, EX, HCI, Z, ZU, ZV
      REAL ELIM, EY, FNU, R1, R2, TAY, XX, YY, R1MACH, ASCLE, RTOL,
     * ATOL, AA, BB, R1M5, TOL
      INTEGER I, IERR, K, KODE, K1, K2, N, NZ, NZ1, NZ2, I1MACH
      DIMENSION CY(N), CWRK(N)
C***FIRST EXECUTABLE STATEMENT  CBESYH
      XX = REAL(Z)
      YY = AIMAG(Z)
      IERR = 0
      NZ=0
      IF (XX.EQ.0.0E0 .AND. YY.EQ.0.0E0) IERR=1
      IF (FNU.LT.0.0E0) IERR=1
      IF (KODE.LT.1 .OR. KODE.GT.2) IERR=1
      IF (N.LT.1) IERR=1
      IF (IERR.NE.0) RETURN
      HCI = CMPLX(0.0E0,0.5E0)
      CALL CBESH(Z, FNU, KODE, 1, N, CY, NZ1, IERR)
      IF (IERR.NE.0.AND.IERR.NE.3) GO TO 170
      CALL CBESH(Z, FNU, KODE, 2, N, CWRK, NZ2, IERR)
      IF (IERR.NE.0.AND.IERR.NE.3) GO TO 170
      NZ = MIN0(NZ1,NZ2)
      IF (KODE.EQ.2) GO TO 60
      DO 50 I=1,N
        CY(I) = HCI*(CWRK(I)-CY(I))
   50 CONTINUE
      RETURN
   60 CONTINUE
      TOL = AMAX1(R1MACH(4),1.0E-18)
      K1 = I1MACH(12)
      K2 = I1MACH(13)
      K = MIN0(IABS(K1),IABS(K2))
      R1M5 = R1MACH(5)
C-----------------------------------------------------------------------
C     ELIM IS THE APPROXIMATE EXPONENTIAL UNDER- AND OVERFLOW LIMIT
C-----------------------------------------------------------------------
      ELIM = 2.303E0*(FLOAT(K)*R1M5-3.0E0)
      R1 = COS(XX)
      R2 = SIN(XX)
      EX = CMPLX(R1,R2)
      EY = 0.0E0
      TAY = ABS(YY+YY)
      IF (TAY.LT.ELIM) EY = EXP(-TAY)
      IF (YY.LT.0.0E0) GO TO 90
      C1 = EX*CMPLX(EY,0.0E0)
      C2 = CONJG(EX)
   70 CONTINUE
      NZ = 0
      RTOL = 1.0E0/TOL
      ASCLE = R1MACH(1)*RTOL*1.0E+3
      DO 80 I=1,N
C       CY(I) = HCI*(C2*CWRK(I)-C1*CY(I))
        ZV = CWRK(I)
        AA=REAL(ZV)
        BB=AIMAG(ZV)
        ATOL=1.0E0
        IF (AMAX1(ABS(AA),ABS(BB)).GT.ASCLE) GO TO 75
          ZV = ZV*CMPLX(RTOL,0.0E0)
          ATOL = TOL
   75   CONTINUE
        ZV = ZV*C2*HCI
        ZV = ZV*CMPLX(ATOL,0.0E0)
        ZU=CY(I)
        AA=REAL(ZU)
        BB=AIMAG(ZU)
        ATOL=1.0E0
        IF (AMAX1(ABS(AA),ABS(BB)).GT.ASCLE) GO TO 85
          ZU = ZU*CMPLX(RTOL,0.0E0)
          ATOL = TOL
   85   CONTINUE
        ZU = ZU*C1*HCI
        ZU = ZU*CMPLX(ATOL,0.0E0)
        CY(I) = ZV - ZU
        IF (CY(I).EQ.CMPLX(0.0E0,0.0E0) .AND. EY.EQ.0.0E0) NZ = NZ + 1
   80 CONTINUE
      RETURN
   90 CONTINUE
      C1 = EX
      C2 = CONJG(EX)*CMPLX(EY,0.0E0)
      GO TO 70
  170 CONTINUE
      NZ = 0
      RETURN
      END
      PROGRAM CQCAI
C
C***DATE WRITTEN   830501   (YYMMDD)
C***REVISION DATE  890801, 930101   (YYMMDD)
C
C     CQCAI IS A QUICK CHECK ROUTINE FOR THE COMPLEX AIRY FUNCTIONS
C     GENERATED BY SUBROUTINES CAIRY AND CBIRY.
C
C     CQCAI GENERATES AIRY FUNCTIONS AND THEIR DERIVATIVES FROM CAIRY
C     AND CBIRY AND CHECKS THEM AGAINST THE WRONSKIAN EVALUATION IN THE
C     REGION -PI/3 .LE. ARG(Z) .LE. PI/3:
C
C                 AI(Z)*BI'(Z)-AI'(Z)*BI(Z)=1/PI.
C
C     IN THE REMAINDER OF THE CUT PLANE, THE IDENTITIES
C
C              AI(Z)  = SQRT(-Z)*( J(-1/3,ZR) + J(1/3,ZR) )/3
C
C              AI'(Z) =        Z*( J(-2/3,ZR) - J(2/3,ZR) )/3
C
C       BI(Z)  =   I*SQRT(-Z/3)*( C1*H(1/3,1,ZR) - C2*H(1/3,2,ZR) )/2
C
C       BI'(Z) = I*(-Z)/SQRT(3)*( C2*H(2/3,1,ZR) - C1*H(2/3,2,ZR) )/2
C
C     ARE CHECKED WHERE ZR = (2/3)(-Z)**(3/2) WITH C1 = EXP(PI*I/6),
C     C2 = CONJG(C1) AND I**2 = -1.
C
C     THE PARAMETER MQC CAN HAVE VALUES 1 (THE DEFAULT) FOR A FASTER,
C     LESS DEFINITIVE TEST OR 2 FOR A SLOWER, MORE DEFINITIVE TEST.
C
C     MACHINE CONSTANTS ARE DEFINED IN FUNCTIONS I1MACH, R1MACH, AND
C     D1MACH. THESE MUST BE SELECTED BY THE USER OR SET ACCORDING TO
C     PROLOGUE INSTRUCTIONS.
C
      COMPLEX CA, CAV, CHI, CI, CONA, CONB, CONC, COND, CON1, CON2, CV,
     *  CW, CY, W, Y, YY, Z, ZR, ZW, SC
      REAL AA, AB, ALIM, ATOL, AV, CT, C13, C23, C43, ZX, ZY,
     * DIG, ELIM, EPS, ER, ERTOL, FNUL, FPI, HPI, PI, R, RL, RPI,
     * R1M4, R1M5, SPI, ST, SLAK, T, TOL, RM, FILM, PI3, R1MACH, TS
      INTEGER I, ICASE, ICL, IL, IR, IRSET, IT, ITL, J, JB, JL, K,
     * KEPS, KODE, K1, K2, LFLG, LUN, NZ, IERR, I1MACH, MQC, IRB, KDO
      DIMENSION  ER(5), T(20), Y(20), YY(20), W(20), KDO(20), KEPS(20)
      DATA LUN /7/
      PARAMETER (MQC=1)
      OPEN(LUN,FILE='CQCAI.OUT')
C-----------------------------------------------------------------------
C     SET PARAMETERS RELATED TO MACHINE CONSTANTS.
C     TOL IS THE APPROXIMATE UNIT ROUNDOFF LIMITED TO 1.0E-18.
C     ELIM IS THE APPROXIMATE EXPONENTIAL OVER- AND UNDERFLOW LIMIT.
C     EXP(-ELIM).LT.EXP(-ALIM)=EXP(-ELIM)/TOL    AND
C     EXP(ELIM).GT.EXP(ALIM)=EXP(ELIM)*TOL       ARE INTERVALS NEAR
C     UNDERFLOW AND OVERFLOW LIMITS WHERE SCALED ARITHMETIC IS DONE.
C     RL IS THE LOWER BOUNDARY OF THE ASYMPTOTIC EXPANSION FOR LARGE Z.
C     DIG = NUMBER OF BASE 10 DIGITS IN TOL = 10**(-DIG).
C     FNUL IS THE LOWER BOUNDARY OF THE ASYMPTOTIC SERIES FOR LARGE FNU
C-----------------------------------------------------------------------
      R1M4 = R1MACH(4)
      TOL = MAX(R1M4,1.0E-18)
      AA = -ALOG10(R1M4)
      K1 = I1MACH(12)
      K2 = I1MACH(13)
      R1M5 = R1MACH(5)
      K = MIN(ABS(K1),ABS(K2))
      ELIM = 2.303E0*(FLOAT(K)*R1M5-3.0E0)
      AB = AA*2.303E0
      ALIM = ELIM + MAX(-AB,-41.45E0)
      DIG = MIN(AA,18.0E0)
      FNUL = 10.0E0 + 6.0E0*(DIG-3.0E0)
      RL = 1.2E0*DIG + 3.0E0
      SLAK = 3.0E0+4.0E0*(-ALOG10(TOL)-7.0E0)/11.0E0
      SLAK = MAX(SLAK,3.0E0)
      ERTOL = TOL*10.0E0**SLAK
      RM = 0.5E0*(ALIM + ELIM)
      RM = MIN(RM,200.0E0)
      RM = MAX(RM,RL+10.0E0)
C-----------------------------------------------------------------------
      WRITE (LUN,99999)
99999 FORMAT (' QUICK CHECK ROUTINE FOR THE AIRY FUNCTIONS FROM CAIRY AN
     *D CBIRY'/)
      WRITE (LUN,99998)
99998 FORMAT (' PARAMETERS TOL,ELIM,ALIM,RL,FNUL,DIG')
      WRITE (LUN,99997) TOL, ELIM, ALIM, RL, FNUL, DIG
99997 FORMAT (6E12.4/)
      ATOL = 100.0E0*TOL
      FPI = ATAN(1.0E0)
      HPI = FPI + FPI
      PI = HPI + HPI
      RPI = 1.0E0/PI
      SPI = PI/6.0E0
      CON1 = CMPLX(COS(SPI),SIN(SPI))
      CON2 = CONJG(CON1)
      PI3 = SPI+SPI
      C13 = 1.0E0/3.0E0
      C23 = C13+C13
      C43 = C23+C23
      AV = SQRT(C13)
      CAV = CMPLX(AV,0.0E0)
      CHI = CMPLX(0.0E0,0.5E0)
      CI = CMPLX(0.0E0,1.0E0)
      WRITE (LUN,99996) MQC
99996 FORMAT (/' CHECKS IN THE (Z,FNU) SPACE WITH MQC = ',I2/)
C-----------------------------------------------------------------------
C     TEST VALUES OF Z IN -PI.LT.ARG(Z).LE.PI
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C     KDO(K), K=1,IL  DETERMINES WHICH OF THE IL ANGLES IN -PI TO PI
C     ARE USE TO COMPUTE VALUES OF Z
C       KDO(K) = 0  MEANS THAT THE INDEX K WILL BE USED FOR ONE OR TWO
C                   VALUES OF Z, DEPENDING ON THE CHOICE OF KEPS(K)
C              = 1  MEANS THAT THE INDEX K AND THE CORRESPONDING ANGLE
C                   WILL BE SKIPPED
C     KEPS(K), K=1,IL DETERMINES WHICH OF THE ANGLES GET INCREMENTED
C     UP AND DOWN TO PUT VALUES OF Z IN REGIONS WHERE DIFFERENT
C     FORMULAE ARE USED.
C       KEPS(K) =0  MEANS THAT THE ANGLE WILL BE USED WITHOUT CHANGE
C               =1  MEANS THAT THE ANGLE WILL BE INCREMENTED UP AND
C                   DOWN BY EPS
C     THE ANGLES TO BE USED ARE STORED IN THE T(I) ARRAY, I=1,ITL
C-----------------------------------------------------------------------
      IF (MQC.NE.2) THEN
        ICL=1
        IL=5
        DO 5 I=1,IL
          KDO(I)=0
          KEPS(I)=0
    5   CONTINUE
      ELSE
        ICL=2
        IL=7
        DO 6 I=1,IL
          KDO(I)=0
          KEPS(I)=0
    6   CONTINUE
        KEPS(2)=1
        KEPS(3)=1
        KEPS(5)=1
        KEPS(6)=1
      ENDIF
      I = 2
      EPS = 0.01E0
      FILM=FLOAT(IL-1)
      T(1) = -PI + EPS
      DO 30 K=2,IL
        IF(KDO(K).EQ.0) THEN
          T(I) = PI*FLOAT(-IL+2*K-1)/FILM
          IF (KEPS(K).EQ.0) GO TO 20
          TS=T(I)
          T(I) = TS - EPS
          I = I + 1
          T(I) = TS + EPS
   20     CONTINUE
          I = I + 1
        ENDIF
   30 CONTINUE
      ITL = I - 1
      LFLG = 0
      DO 180 ICASE=1,ICL
        DO 170 KODE=1,2
          DO 160 IRSET=1,3
            IRB = MIN(IRSET,2)
            DO 150 IR=IRB,4
              GO TO (40, 50, 60), IRSET
   40         CONTINUE
              R = (0.2E0*FLOAT(4-IR)+2.0E0*FLOAT(IR-1))/3.0E0
              GO TO 70
   50         CONTINUE
              R = (2.0E0*FLOAT(4-IR)+RL*FLOAT(IR-1))/3.0E0
              GO TO 70
   60         CONTINUE
              R = (RL*FLOAT(4-IR)+RM*FLOAT(IR-1))/3.0E0
   70         CONTINUE
              DO 140 IT=1,ITL
                CT = COS(T(IT))
                ST = SIN(T(IT))
                IF (ABS(CT).LT.ATOL) CT = 0.0E0
                IF (ABS(ST).LT.ATOL) ST = 0.0E0
                ZX = R*CT
                ZY = R*ST
                Z = CMPLX(ZX,ZY)
                IF(ABS(T(IT)).LE.PI3) THEN
C-----------------------------------------------------------------------
C     WRONSKIAN CHECK IN -PI/3.LT.ARG(Z).LT.PI/3, TEST #1
C-----------------------------------------------------------------------
                  CALL CAIRY(Z, 0, KODE, Y(1), NZ, IERR)
                  IF (NZ.NE.0.OR.IERR.NE.0) GO TO 140
                  CALL CAIRY(Z, 1, KODE, Y(2), NZ, IERR)
                  IF (NZ.NE.0.OR.IERR.NE.0) GO TO 140
                  CALL CBIRY(Z, 0, KODE, W(1), IERR)
                  CALL CBIRY(Z, 1, KODE, W(2), IERR)
                  CW = Y(1)*W(2)
                  CY = Y(2)*W(1)
                  CV = CMPLX(RPI,0.0E0)
                  IF (KODE.EQ.2) THEN
                    ZR=Z
                    CA=CSQRT(ZR)
                    ZR=ZR*CA*CMPLX(C23,0.0E0)
                    AA=REAL(ZR)
                    AA=ABS(AA)
                    CA = ZR - CMPLX(AA,0.0E0)
                    CV = CEXP(CA)*CV
                  ENDIF
                  CY = CW - CY - CV
                  ER(1) = CABS(CY)/CABS(CV)
                  JB = 1
                  JL = 1
                ELSE
C-----------------------------------------------------------------------
C     CHECKS IN -PI.LT.ARG(Z).LT.-PI/3 AND PI/3.LT.ARG(Z).LE.PI
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C     CHECK AI    TEST #2
C-----------------------------------------------------------------------
                  CALL CAIRY(Z, 0, KODE, Y(2), NZ, IERR)
                  IF (NZ.NE.0.OR.IERR.NE.0) GO TO 140
                  ZR=-Z
                  CV=CSQRT(ZR)
                  ZR=ZR*CV*CMPLX(C23,0.0E0)
                  CALL CBESJ(ZR,C23,KODE,2,YY,NZ,IERR)
                  IF (NZ.NE.0.OR.IERR.NE.0) GO TO 140
                  CY=CMPLX(C43,0.0E0)*YY(1)/ZR - YY(2)
                  CA = YY(1)
                  CALL CBESJ(ZR,C13,KODE,2,YY,NZ,IERR)
                  IF (NZ.NE.0.OR.IERR.NE.0) GO TO 140
                  IF (KODE.EQ.2) THEN
                    AB = AIMAG(ZR)
                    AB = ABS(AB)
                    CW = CSQRT(Z)
                    ZW = Z*CW*CMPLX(C23,0.0E0)
                    CW = ZW+CMPLX(AB,0.0E0)
                    CW = CEXP(CW)
                    YY(1) = YY(1)*CW
                    YY(2) = YY(2)*CW
                    CY = CY*CW
                    CA = CA*CW
                    SC = CW
                  ENDIF
                  CW = CV*CMPLX(C13,0.0E0)
                  W(2) = CW*(YY(1)+CY)
                  ER(2) = CABS(Y(2)-W(2))
                  IF (ZY.NE.0.0D0.OR.ZX.GE.0.0D0) THEN
                    ER(2)=ER(2)/CABS(Y(2))
                  ELSE
                    IF (KODE.EQ.2) THEN
                      ER(2) = ER(2)/CABS(SC)
                    ENDIF
                  ENDIF
C-----------------------------------------------------------------------
C     CHECK AI'   TEST #3
C-----------------------------------------------------------------------
                  CY=CMPLX(C23,0.0E0)*YY(1)/ZR - YY(2)
                  W(3) = Z*CMPLX(C13,0.0E0)*(CY-CA)
                  CALL CAIRY(Z,1,KODE,Y(3),NZ,IERR)
                  ER(3) = CABS(Y(3)-W(3))
                  IF (ZY.NE.0.0D0.OR.ZX.GE.0.0D0) THEN
                    ER(3)=ER(3)/CABS(Y(3))
                  ELSE
                    IF (KODE.EQ.2) THEN
                      ER(3) = ER(3)/CABS(SC)
                    ENDIF
                  ENDIF
C-----------------------------------------------------------------------
C     CHECK BI    TEST #4
C-----------------------------------------------------------------------
                  CALL CBESH(ZR,C13,KODE,1,1,Y,NZ,IERR)
                  IF (NZ.NE.0.OR.IERR.NE.0) GO TO 140
                  CALL CBESH(ZR,C13,KODE,2,1,YY,NZ,IERR)
                  IF (NZ.NE.0.OR.IERR.NE.0) GO TO 140
                  CONA = CON1
                  CONB = CON2
                  CONC = CON2
                  COND = CON1
                  IF (KODE.EQ.2) THEN
                    AA = REAL(ZW)
                    AA = ABS(AA)
                    ZW =  CI*ZR - CMPLX(AA,0.0E0)
                    CW = CEXP(ZW)
                    CONA = CONA*CW
                    CONC = CONC*CW
                    ZW = -CI*ZR - CMPLX(AA,0.0E0)
                    CW = CEXP(ZW)
                    CONB = CONB*CW
                    COND = COND*CW
                    SC = CW
                  ENDIF
                  CW = CONA*Y(1)-CONB*YY(1)
                  CW = CV*CAV*CW
                  W(4) = CW*CHI
                  CALL CBIRY(Z,0,KODE,Y(4),IERR)
                  ER(4) = CABS(Y(4)-W(4))
                  IF (ZY.NE.0.0D0.OR.ZX.GE.0.0D0) THEN
                    ER(4)=ER(4)/CABS(Y(4))
                  ELSE
                    IF (KODE.EQ.2) THEN
                      ER(4) = ER(4)/CABS(SC)
                    ENDIF
                  ENDIF
C-----------------------------------------------------------------------
C     CHECK BI'   TEST #5
C-----------------------------------------------------------------------
                  CALL CBESH(ZR,C23,KODE,1,1,Y,NZ,IERR)
                  IF (NZ.NE.0.OR.IERR.NE.0) GO TO 140
                  CALL CBESH(ZR,C23,KODE,2,1,YY,NZ,IERR)
                  IF (NZ.NE.0.OR.IERR.NE.0) GO TO 140
                  CW = CONC*Y(1)-COND*YY(1)
                  CW = -Z*CAV*CW
                  W(5) = CW*CHI
                  CALL CBIRY(Z,1,KODE,Y(5),IERR)
                  ER(5) = CABS(Y(5)-W(5))
                  IF (ZY.NE.0.0D0.OR.ZX.GE.0.0D0) THEN
                    ER(5)=ER(5)/CABS(Y(5))
                  ELSE
                    IF (KODE.EQ.2) THEN
                      ER(5) = ER(5)/CABS(SC)
                    ENDIF
                  ENDIF
                  JB = 2
                  JL = 5
                ENDIF
                DO 190 J=JB,JL
                IF (ER(J).LT.ERTOL) GO TO 190
                IF (LFLG.EQ.1) GO TO 130
                WRITE (LUN,99995) ERTOL
99995           FORMAT (/' CASES WHICH VIOLATE THE RELATIVE ERROR TEST W
     *ITH ERTOL =', E12.4/)
                WRITE (LUN,99994)
99994           FORMAT (/' OUTPUT FORMAT'/' KODE,IR,IT,IRSET,ICASE')
                WRITE (LUN,99993)
99993           FORMAT (' ER'/' I, Z, Y(I), W(I), ON THE I-TH TEST, I=1,
     *5'/)
                LFLG = 1
  130           CONTINUE
                WRITE (LUN,99992) KODE, IR, IT, IRSET, ICASE
99992           FORMAT (5I5)
                WRITE (LUN,99991) ER(J)
                WRITE (LUN,99989) J, Z, Y(J), W(J)
99991           FORMAT (E12.4)
99989           FORMAT (I5,6E12.4)
  190           CONTINUE
  140         CONTINUE
  150       CONTINUE
  160     CONTINUE
  170   CONTINUE
  180 CONTINUE
      IF (LFLG.EQ.0) WRITE (LUN,99990)
99990 FORMAT (/' QUICK CHECKS OK'/)
      STOP
      END

C*** machcon.f
-----------------------------------------------------------------
C>>>  MACHCON.FOR:  Machine constants
-----------------------------------------------------------------

*DECK I1MACH
      INTEGER FUNCTION I1MACH(I)
C***BEGIN PROLOGUE  I1MACH
C***DATE WRITTEN   750101   (YYMMDD)
C***REVISION DATE  890213   (YYMMDD)
C***CATEGORY NO.  R1
C***KEYWORDS  LIBRARY=SLATEC,TYPE=INTEGER(I1MACH-I),MACHINE CONSTANTS
C***AUTHOR  FOX, P. A., (BELL LABS)
C           HALL, A. D., (BELL LABS)
C           SCHRYER, N. L., (BELL LABS)
C***PURPOSE  Returns integer machine dependent constants
C***DESCRIPTION
C
C     I1MACH can be used to obtain machine-dependent parameters
C     for the local machine environment.  It is a function
C     subroutine with one (input) argument, and can be called
C     as follows, for example
C
C          K = I1MACH(I)
C
C     where I=1,...,16.  The (output) value of K above is
C     determined by the (input) value of I.  The results for
C     various values of I are discussed below.
C
C  I/O unit numbers.
C    I1MACH( 1) = the standard input unit.
C    I1MACH( 2) = the standard output unit.
C    I1MACH( 3) = the standard punch unit.
C    I1MACH( 4) = the standard error message unit.
C
C  Words.
C    I1MACH( 5) = the number of bits per integer storage unit.
C    I1MACH( 6) = the number of characters per integer storage unit.
C
C  Integers.
C    assume integers are represented in the S-digit, base-A form
C
C               sign ( X(S-1)*A**(S-1) + ... + X(1)*A + X(0) )
C
C               where 0 .LE. X(I) .LT. A for I=0,...,S-1.
C    I1MACH( 7) = A, the base.
C    I1MACH( 8) = S, the number of base-A digits.
C    I1MACH( 9) = A**S - 1, the largest magnitude.
C
C  Floating-Point Numbers.
C    Assume floating-point numbers are represented in the T-digit,
C    base-B form
C               sign (B**E)*( (X(1)/B) + ... + (X(T)/B**T) )
C
C               where 0 .LE. X(I) .LT. B for I=1,...,T,
C               0 .LT. X(1), and EMIN .LE. E .LE. EMAX.
C    I1MACH(10) = B, the base.
C
C  Single-Precision
C    I1MACH(11) = T, the number of base-B digits.
C    I1MACH(12) = EMIN, the smallest exponent E.
C    I1MACH(13) = EMAX, the largest exponent E.
C
C  Double-Precision
C    I1MACH(14) = T, the number of base-B digits.
C    I1MACH(15) = EMIN, the smallest exponent E.
C    I1MACH(16) = EMAX, the largest exponent E.
C
C  To alter this function for a particular environment,
C  the desired set of DATA statements should be activated by
C  removing the C from column 1.  Also, the values of
C  I1MACH(1) - I1MACH(4) should be checked for consistency
C  with the local operating system.
C
C***REFERENCES  FOX P.A., HALL A.D., SCHRYER N.L.,*FRAMEWORK FOR A
C                 PORTABLE LIBRARY*, ACM TRANSACTIONS ON MATHEMATICAL
C                 SOFTWARE, VOL. 4, NO. 2, JUNE 1978, PP. 177-188.
C***ROUTINES CALLED  (NONE)
C***END PROLOGUE  I1MACH
C
      INTEGER IMACH(16),OUTPUT
      SAVE IMACH
      EQUIVALENCE (IMACH(4),OUTPUT)
C
C     MACHINE CONSTANTS FOR THE AMIGA
C     ABSOFT COMPILER
C
C     DATA IMACH(1) /    5 /
C     DATA IMACH(2) /    6 /
C     DATA IMACH(3) /    5 /
C     DATA IMACH(4) /    6 /
C     DATA IMACH(5) /   32 /
C     DATA IMACH(6) /    4 /
C     DATA IMACH(7) /    2 /
C     DATA IMACH(8) /   31 /
C     DATA IMACH(9) / 2147483647 /
C     DATA IMACH(10)/    2 /
C     DATA IMACH(11)/   24 /
C     DATA IMACH(12)/ -126 /
C     DATA IMACH(13)/  127 /
C     DATA IMACH(14)/   53 /
C     DATA IMACH(15)/ -1022 /
C     DATA IMACH(16)/  1023 /
C
C     MACHINE CONSTANTS FOR THE APOLLO
C
C     DATA IMACH(1) /    5 /
C     DATA IMACH(2) /    6 /
C     DATA IMACH(3) /    6 /
C     DATA IMACH(4) /    6 /
C     DATA IMACH(5) /   32 /
C     DATA IMACH(6) /    4 /
C     DATA IMACH(7) /    2 /
C     DATA IMACH(8) /   31 /
C     DATA IMACH(9) / 2147483647 /
C     DATA IMACH(10)/    2 /
C     DATA IMACH(11)/   24 /
C     DATA IMACH(12)/ -125 /
C     DATA IMACH(13)/  129 /
C     DATA IMACH(14)/   53 /
C     DATA IMACH(15)/ -1021 /
C     DATA IMACH(16)/  1025 /
C
C     MACHINE CONSTANTS FOR THE BURROUGHS 1700 SYSTEM
C
C     DATA IMACH( 1) /    7 /
C     DATA IMACH( 2) /    2 /
C     DATA IMACH( 3) /    2 /
C     DATA IMACH( 4) /    2 /
C     DATA IMACH( 5) /   36 /
C     DATA IMACH( 6) /    4 /
C     DATA IMACH( 7) /    2 /
C     DATA IMACH( 8) /   33 /
C     DATA IMACH( 9) / Z1FFFFFFFF /
C     DATA IMACH(10) /    2 /
C     DATA IMACH(11) /   24 /
C     DATA IMACH(12) / -256 /
C     DATA IMACH(13) /  255 /
C     DATA IMACH(14) /   60 /
C     DATA IMACH(15) / -256 /
C     DATA IMACH(16) /  255 /
C
C     MACHINE CONSTANTS FOR THE BURROUGHS 5700 SYSTEM
C
C     DATA IMACH( 1) /   5 /
C     DATA IMACH( 2) /   6 /
C     DATA IMACH( 3) /   7 /
C     DATA IMACH( 4) /   6 /
C     DATA IMACH( 5) /  48 /
C     DATA IMACH( 6) /   6 /
C     DATA IMACH( 7) /   2 /
C     DATA IMACH( 8) /  39 /
C     DATA IMACH( 9) / O0007777777777777 /
C     DATA IMACH(10) /   8 /
C     DATA IMACH(11) /  13 /
C     DATA IMACH(12) / -50 /
C     DATA IMACH(13) /  76 /
C     DATA IMACH(14) /  26 /
C     DATA IMACH(15) / -50 /
C     DATA IMACH(16) /  76 /
C
C     MACHINE CONSTANTS FOR THE BURROUGHS 6700/7700 SYSTEMS
C
C     DATA IMACH( 1) /   5 /
C     DATA IMACH( 2) /   6 /
C     DATA IMACH( 3) /   7 /
C     DATA IMACH( 4) /   6 /
C     DATA IMACH( 5) /  48 /
C     DATA IMACH( 6) /   6 /
C     DATA IMACH( 7) /   2 /
C     DATA IMACH( 8) /  39 /
C     DATA IMACH( 9) / O0007777777777777 /
C     DATA IMACH(10) /   8 /
C     DATA IMACH(11) /  13 /
C     DATA IMACH(12) / -50 /
C     DATA IMACH(13) /  76 /
C     DATA IMACH(14) /  26 /
C     DATA IMACH(15) / -32754 /
C     DATA IMACH(16) /  32780 /
C
C     MACHINE CONSTANTS FOR THE CDC 170/180 SERIES USING NOS/VE
C
C     DATA IMACH( 1) /     5 /
C     DATA IMACH( 2) /     6 /
C     DATA IMACH( 3) /     7 /
C     DATA IMACH( 4) /     6 /
C     DATA IMACH( 5) /    64 /
C     DATA IMACH( 6) /     8 /
C     DATA IMACH( 7) /     2 /
C     DATA IMACH( 8) /    63 /
C     DATA IMACH( 9) / 9223372036854775807 /
C     DATA IMACH(10) /     2 /
C     DATA IMACH(11) /    47 /
C     DATA IMACH(12) / -4095 /
C     DATA IMACH(13) /  4094 /
C     DATA IMACH(14) /    94 /
C     DATA IMACH(15) / -4095 /
C     DATA IMACH(16) /  4094 /
C
C     MACHINE CONSTANTS FOR THE CDC 6000/7000 SERIES
C
C     DATA IMACH( 1) /    5 /
C     DATA IMACH( 2) /    6 /
C     DATA IMACH( 3) /    7 /
C     DATA IMACH( 4) /6LOUTPUT/
C     DATA IMACH( 5) /   60 /
C     DATA IMACH( 6) /   10 /
C     DATA IMACH( 7) /    2 /
C     DATA IMACH( 8) /   48 /
C     DATA IMACH( 9) / 00007777777777777777B /
C     DATA IMACH(10) /    2 /
C     DATA IMACH(11) /   47 /
C     DATA IMACH(12) / -929 /
C     DATA IMACH(13) / 1070 /
C     DATA IMACH(14) /   94 /
C     DATA IMACH(15) / -929 /
C     DATA IMACH(16) / 1069 /
C
C     MACHINE CONSTANTS FOR THE CELERITY C1260
C
C     DATA IMACH(1) /    5 /
C     DATA IMACH(2) /    6 /
C     DATA IMACH(3) /    6 /
C     DATA IMACH(4) /    0 /
C     DATA IMACH(5) /   32 /
C     DATA IMACH(6) /    4 /
C     DATA IMACH(7) /    2 /
C     DATA IMACH(8) /   31 /
C     DATA IMACH(9) / Z'7FFFFFFF' /
C     DATA IMACH(10)/    2 /
C     DATA IMACH(11)/   24 /
C     DATA IMACH(12)/ -126 /
C     DATA IMACH(13)/  127 /
C     DATA IMACH(14)/   53 /
C     DATA IMACH(15)/ -1022 /
C     DATA IMACH(16)/  1023 /
C
C     MACHINE CONSTANTS FOR THE CONVEX C-1
C
C     DATA IMACH( 1) /     5/
C     DATA IMACH( 2) /     6/
C     DATA IMACH( 3) /     7/
C     DATA IMACH( 4) /     6/
C     DATA IMACH( 5) /    32/
C     DATA IMACH( 6) /     4/
C     DATA IMACH( 7) /     2/
C     DATA IMACH( 8) /    31/
C     DATA IMACH( 9) /2147483647/
C     DATA IMACH(10) /     2/
C     DATA IMACH(11) /    24/
C     DATA IMACH(12) /  -128/
C     DATA IMACH(13) /   127/
C     DATA IMACH(14) /    53/
C     DATA IMACH(15) / -1024/
C     DATA IMACH(16) /  1023/
C
C     MACHINE CONSTANTS FOR THE CRAY-1
C
C     DATA IMACH( 1) /   100 /
C     DATA IMACH( 2) /   101 /
C     DATA IMACH( 3) /   102 /
C     DATA IMACH( 4) /   101 /
C     DATA IMACH( 5) /    64 /
C     DATA IMACH( 6) /     8 /
C     DATA IMACH( 7) /     2 /
C     DATA IMACH( 8) /    63 /
C     DATA IMACH( 9) /  777777777777777777777B /
C     DATA IMACH(10) /     2 /
C     DATA IMACH(11) /    47 /
C     DATA IMACH(12) / -8189 /
C     DATA IMACH(13) /  8190 /
C     DATA IMACH(14) /    94 /
C     DATA IMACH(15) / -8099 /
C     DATA IMACH(16) /  8190 /
C
C     MACHINE CONSTANTS FOR THE DATA GENERAL ECLIPSE S/200
C
C     DATA IMACH( 1) /   11 /
C     DATA IMACH( 2) /   12 /
C     DATA IMACH( 3) /    8 /
C     DATA IMACH( 4) /   10 /
C     DATA IMACH( 5) /   16 /
C     DATA IMACH( 6) /    2 /
C     DATA IMACH( 7) /    2 /
C     DATA IMACH( 8) /   15 /
C     DATA IMACH( 9) /32767 /
C     DATA IMACH(10) /   16 /
C     DATA IMACH(11) /    6 /
C     DATA IMACH(12) /  -64 /
C     DATA IMACH(13) /   63 /
C     DATA IMACH(14) /   14 /
C     DATA IMACH(15) /  -64 /
C     DATA IMACH(16) /   63 /
C
C     MACHINE CONSTANTS FOR THE ELXSI 6400
C
C     DATA IMACH( 1) /     5/
C     DATA IMACH( 2) /     6/
C     DATA IMACH( 3) /     6/
C     DATA IMACH( 4) /     6/
C     DATA IMACH( 5) /    32/
C     DATA IMACH( 6) /     4/
C     DATA IMACH( 7) /     2/
C     DATA IMACH( 8) /    32/
C     DATA IMACH( 9) /2147483647/
C     DATA IMACH(10) /     2/
C     DATA IMACH(11) /    24/
C     DATA IMACH(12) /  -126/
C     DATA IMACH(13) /   127/
C     DATA IMACH(14) /    53/
C     DATA IMACH(15) / -1022/
C     DATA IMACH(16) /  1023/
C
C     MACHINE CONSTANTS FOR THE HARRIS 220
C
C     DATA IMACH( 1) /       5 /
C     DATA IMACH( 2) /       6 /
C     DATA IMACH( 3) /       0 /
C     DATA IMACH( 4) /       6 /
C     DATA IMACH( 5) /      24 /
C     DATA IMACH( 6) /       3 /
C     DATA IMACH( 7) /       2 /
C     DATA IMACH( 8) /      23 /
C     DATA IMACH( 9) / 8388607 /
C     DATA IMACH(10) /       2 /
C     DATA IMACH(11) /      23 /
C     DATA IMACH(12) /    -127 /
C     DATA IMACH(13) /     127 /
C     DATA IMACH(14) /      38 /
C     DATA IMACH(15) /    -127 /
C     DATA IMACH(16) /     127 /
C
C     MACHINE CONSTANTS FOR THE HONEYWELL 600/6000 SERIES
C
C     DATA IMACH( 1) /    5 /
C     DATA IMACH( 2) /    6 /
C     DATA IMACH( 3) /   43 /
C     DATA IMACH( 4) /    6 /
C     DATA IMACH( 5) /   36 /
C     DATA IMACH( 6) /    6 /
C     DATA IMACH( 7) /    2 /
C     DATA IMACH( 8) /   35 /
C     DATA IMACH( 9) / O377777777777 /
C     DATA IMACH(10) /    2 /
C     DATA IMACH(11) /   27 /
C     DATA IMACH(12) / -127 /
C     DATA IMACH(13) /  127 /
C     DATA IMACH(14) /   63 /
C     DATA IMACH(15) / -127 /
C     DATA IMACH(16) /  127 /
C
C     MACHINE CONSTANTS FOR THE HP 2100
C     3 WORD DOUBLE PRECISION OPTION WITH FTN4
C
C     DATA IMACH(1) /      5/
C     DATA IMACH(2) /      6 /
C     DATA IMACH(3) /      4 /
C     DATA IMACH(4) /      1 /
C     DATA IMACH(5) /     16 /
C     DATA IMACH(6) /      2 /
C     DATA IMACH(7) /      2 /
C     DATA IMACH(8) /     15 /
C     DATA IMACH(9) /  32767 /
C     DATA IMACH(10)/      2 /
C     DATA IMACH(11)/     23 /
C     DATA IMACH(12)/   -128 /
C     DATA IMACH(13)/    127 /
C     DATA IMACH(14)/     39 /
C     DATA IMACH(15)/   -128 /
C     DATA IMACH(16)/    127 /
C
C     MACHINE CONSTANTS FOR THE HP 2100
C     4 WORD DOUBLE PRECISION OPTION WITH FTN4
C
C     DATA IMACH(1) /      5 /
C     DATA IMACH(2) /      6 /
C     DATA IMACH(3) /      4 /
C     DATA IMACH(4) /      1 /
C     DATA IMACH(5) /     16 /
C     DATA IMACH(6) /      2 /
C     DATA IMACH(7) /      2 /
C     DATA IMACH(8) /     15 /
C     DATA IMACH(9) /  32767 /
C     DATA IMACH(10)/      2 /
C     DATA IMACH(11)/     23 /
C     DATA IMACH(12)/   -128 /
C     DATA IMACH(13)/    127 /
C     DATA IMACH(14)/     55 /
C     DATA IMACH(15)/   -128 /
C     DATA IMACH(16)/    127 /
C
C     MACHINE CONSTANTS FOR THE HP 9000
C
C     DATA IMACH(1)  /    5 /
C     DATA IMACH(2)  /    6 /
C     DATA IMACH(3)  /    6 /
C     DATA IMACH(3)  /    7 /
C     DATA IMACH(5)  /   32 /
C     DATA IMACH(6)  /    4 /
C     DATA IMACH(7)  /    2 /
C     DATA IMACH(8)  /   32 /
C     DATA IMACH(9)  /2147483647 /
C     DATA IMACH(10) /    2 /
C     DATA IMACH(11) /   24 /
C     DATA IMACH(12) / -126 /
C     DATA IMACH(13) /  127 /
C     DATA IMACH(14) /   53 /
C     DATA IMACH(15) /-1015 /
C     DATA IMACH(16) / 1017 /
C
C     MACHINE CONSTANTS FOR THE IBM 360/370 SERIES,
C     THE XEROX SIGMA 5/7/9, THE SEL SYSTEMS 85/86, AND
C     THE PERKIN ELMER (INTERDATA) 7/32.
C
C     DATA IMACH( 1) /   5 /
C     DATA IMACH( 2) /   6 /
C     DATA IMACH( 3) /   7 /
C     DATA IMACH( 4) /   6 /
C     DATA IMACH( 5) /  32 /
C     DATA IMACH( 6) /   4 /
C     DATA IMACH( 7) /  16 /
C     DATA IMACH( 8) /  31 /
C     DATA IMACH( 9) / Z7FFFFFFF /
C     DATA IMACH(10) /  16 /
C     DATA IMACH(11) /   6 /
C     DATA IMACH(12) / -64 /
C     DATA IMACH(13) /  63 /
C     DATA IMACH(14) /  14 /
C     DATA IMACH(15) / -64 /
C     DATA IMACH(16) /  63 /
C
C     MACHINE CONSTANTS FOR THE IBM PC
C
C     DATA IMACH( 1) /     5 /
C     DATA IMACH( 2) /     6 /
C     DATA IMACH( 3) /     0 /
C     DATA IMACH( 4) /     0 /
C     DATA IMACH( 5) /    32 /
C     DATA IMACH( 6) /     4 /
C     DATA IMACH( 7) /     2 /
C     DATA IMACH( 8) /    31 /
C     DATA IMACH( 9) / 2147483647 /
C     DATA IMACH(10) /     2 /
C     DATA IMACH(11) /    24 /
C     DATA IMACH(12) /  -125 /
C     DATA IMACH(13) /   127 /
C     DATA IMACH(14) /    53 /
C     DATA IMACH(15) / -1021 /
C     DATA IMACH(16) /  1023 /
C
C     MACHINE CONSTANTS FOR THE IBM RS 6000
C
C     DATA IMACH( 1) /          5 /
C     DATA IMACH( 2) /          6 /
C     DATA IMACH( 3) /          6 /
C     DATA IMACH( 4) /          0 /
C     DATA IMACH( 5) /         32 /
C     DATA IMACH( 6) /          4 /
C     DATA IMACH( 7) /          2 /
C     DATA IMACH( 8) /         31 /
C     DATA IMACH( 9) / 2147483647 /
C     DATA IMACH(10) /          2 /
C     DATA IMACH(11) /         24 /
C     DATA IMACH(12) /       -125 /
C     DATA IMACH(13) /        128 /
C     DATA IMACH(14) /         53 /
C     DATA IMACH(15) /      -1021 /
C     DATA IMACH(16) /       1024 /
C
C     MACHINE CONSTANTS FOR THE PDP-10 (KA PROCESSOR)
C
C     DATA IMACH( 1) /    5 /
C     DATA IMACH( 2) /    6 /
C     DATA IMACH( 3) /    5 /
C     DATA IMACH( 4) /    6 /
C     DATA IMACH( 5) /   36 /
C     DATA IMACH( 6) /    5 /
C     DATA IMACH( 7) /    2 /
C     DATA IMACH( 8) /   35 /
C     DATA IMACH( 9) / "377777777777 /
C     DATA IMACH(10) /    2 /
C     DATA IMACH(11) /   27 /
C     DATA IMACH(12) / -128 /
C     DATA IMACH(13) /  127 /
C     DATA IMACH(14) /   54 /
C     DATA IMACH(15) / -101 /
C     DATA IMACH(16) /  127 /
C
C     MACHINE CONSTANTS FOR THE PDP-10 (KI PROCESSOR)
C
C     DATA IMACH( 1) /    5 /
C     DATA IMACH( 2) /    6 /
C     DATA IMACH( 3) /    5 /
C     DATA IMACH( 4) /    6 /
C     DATA IMACH( 5) /   36 /
C     DATA IMACH( 6) /    5 /
C     DATA IMACH( 7) /    2 /
C     DATA IMACH( 8) /   35 /
C     DATA IMACH( 9) / "377777777777 /
C     DATA IMACH(10) /    2 /
C     DATA IMACH(11) /   27 /
C     DATA IMACH(12) / -128 /
C     DATA IMACH(13) /  127 /
C     DATA IMACH(14) /   62 /
C     DATA IMACH(15) / -128 /
C     DATA IMACH(16) /  127 /
C
C     MACHINE CONSTANTS FOR PDP-11 FORTRAN SUPPORTING
C     32-BIT INTEGER ARITHMETIC.
C
C     DATA IMACH( 1) /    5 /
C     DATA IMACH( 2) /    6 /
C     DATA IMACH( 3) /    5 /
C     DATA IMACH( 4) /    6 /
C     DATA IMACH( 5) /   32 /
C     DATA IMACH( 6) /    4 /
C     DATA IMACH( 7) /    2 /
C     DATA IMACH( 8) /   31 /
C     DATA IMACH( 9) / 2147483647 /
C     DATA IMACH(10) /    2 /
C     DATA IMACH(11) /   24 /
C     DATA IMACH(12) / -127 /
C     DATA IMACH(13) /  127 /
C     DATA IMACH(14) /   56 /
C     DATA IMACH(15) / -127 /
C     DATA IMACH(16) /  127 /
C
C     MACHINE CONSTANTS FOR PDP-11 FORTRAN SUPPORTING
C     16-BIT INTEGER ARITHMETIC.
C
C     DATA IMACH( 1) /    5 /
C     DATA IMACH( 2) /    6 /
C     DATA IMACH( 3) /    5 /
C     DATA IMACH( 4) /    6 /
C     DATA IMACH( 5) /   16 /
C     DATA IMACH( 6) /    2 /
C     DATA IMACH( 7) /    2 /
C     DATA IMACH( 8) /   15 /
C     DATA IMACH( 9) / 32767 /
C     DATA IMACH(10) /    2 /
C     DATA IMACH(11) /   24 /
C     DATA IMACH(12) / -127 /
C     DATA IMACH(13) /  127 /
C     DATA IMACH(14) /   56 /
C     DATA IMACH(15) / -127 /
C     DATA IMACH(16) /  127 /
C
C     MACHINE CONSTANTS FOR THE SILICON GRAPHICS IRIS
C
C     DATA IMACH( 1) /     5 /
C     DATA IMACH( 2) /     6 /
C     DATA IMACH( 3) /     6 /
C     DATA IMACH( 4) /     0 /
C     DATA IMACH( 5) /    32 /
C     DATA IMACH( 6) /     4 /
C     DATA IMACH( 7) /     2 /
C     DATA IMACH( 8) /    31 /
C     DATA IMACH( 9) / 2147483647 /
C     DATA IMACH(10) /     2 /
C     DATA IMACH(11) /    23 /
C     DATA IMACH(12) /  -126 /
C     DATA IMACH(13) /   127 /
C     DATA IMACH(14) /    52 /
C     DATA IMACH(15) / -1022 /
C     DATA IMACH(16) /  1023 /
C
C     MACHINE CONSTANTS FOR THE SUN
C
C     DATA IMACH(1) /    5 /
C     DATA IMACH(2) /    6 /
C     DATA IMACH(3) /    6 /
C     DATA IMACH(4) /    6 /
C     DATA IMACH(5) /   32 /
C     DATA IMACH(6) /    4 /
C     DATA IMACH(7) /    2 /
C     DATA IMACH(8) /   31 /
C     DATA IMACH(9) /2147483647 /
C     DATA IMACH(10)/    2 /
C     DATA IMACH(11)/   24 /
C     DATA IMACH(12)/ -125 /
C     DATA IMACH(13)/  128 /
C     DATA IMACH(14)/   53 /
C     DATA IMACH(15)/ -1021 /
C     DATA IMACH(16)/  1024 /
C
C     MACHINE CONSTANTS FOR THE UNIVAC 1100 SERIES FTN COMPILER
C
C
C     DATA IMACH( 1) /    5 /
C     DATA IMACH( 2) /    6 /
C     DATA IMACH( 3) /    1 /
C     DATA IMACH( 4) /    6 /
C     DATA IMACH( 5) /   36 /
C     DATA IMACH( 6) /    4 /
C     DATA IMACH( 7) /    2 /
C     DATA IMACH( 8) /   35 /
C     DATA IMACH( 9) / O377777777777 /
C     DATA IMACH(10) /    2 /
C     DATA IMACH(11) /   27 /
C     DATA IMACH(12) / -128 /
C     DATA IMACH(13) /  127 /
C     DATA IMACH(14) /   60 /
C     DATA IMACH(15) /-1024 /
C     DATA IMACH(16) / 1023 /
C
C     MACHINE CONSTANTS FOR THE VAX 11/780
C
C     DATA IMACH(1) /    5 /
C     DATA IMACH(2) /    6 /
C     DATA IMACH(3) /    5 /
C     DATA IMACH(4) /    6 /
C     DATA IMACH(5) /   32 /
C     DATA IMACH(6) /    4 /
C     DATA IMACH(7) /    2 /
C     DATA IMACH(8) /   31 /
C     DATA IMACH(9) /2147483647 /
C     DATA IMACH(10)/    2 /
C     DATA IMACH(11)/   24 /
C     DATA IMACH(12)/ -127 /
C     DATA IMACH(13)/  127 /
C     DATA IMACH(14)/   56 /
C     DATA IMACH(15)/ -127 /
C     DATA IMACH(16)/  127 /
C
C     MACHINE CONSTANTS FOR THE VAX 11/780, G-FLOAT OPTION
C
C     DATA IMACH(1) /    5 /
C     DATA IMACH(2) /    6 /
C     DATA IMACH(3) /    5 /
C     DATA IMACH(4) /    6 /
C     DATA IMACH(5) /   32 /
C     DATA IMACH(6) /    4 /
C     DATA IMACH(7) /    2 /
C     DATA IMACH(8) /   31 /
C     DATA IMACH(9) /2147483647 /
C     DATA IMACH(10)/    2 /
C     DATA IMACH(11)/   24 /
C     DATA IMACH(12)/ -127 /
C     DATA IMACH(13)/  127 /
C     DATA IMACH(14)/   53 /
C     DATA IMACH(15)/ -1022 /
C     DATA IMACH(16)/  1023 /
C
C     MACHINE CONSTANTS FOR THE Z80 MICROPROCESSOR
C
C     DATA IMACH( 1) /     1/
C     DATA IMACH( 2) /     1/
C     DATA IMACH( 3) /     0/
C     DATA IMACH( 4) /     1/
C     DATA IMACH( 5) /    16/
C     DATA IMACH( 6) /     2/
C     DATA IMACH( 7) /     2/
C     DATA IMACH( 8) /    15/
C     DATA IMACH( 9) / 32767/
C     DATA IMACH(10) /     2/
C     DATA IMACH(11) /    24/
C     DATA IMACH(12) /  -127/
C     DATA IMACH(13) /   127/
C     DATA IMACH(14) /    56/
C     DATA IMACH(15) /  -127/
C     DATA IMACH(16) /   127/
C
C
C***FIRST EXECUTABLE STATEMENT  I1MACH
      IF (I .LT. 1  .OR.  I .GT. 16) GO TO 10
C
      I1MACH = IMACH(I)
      RETURN
C
   10 CONTINUE
      WRITE (UNIT = OUTPUT, FMT = 9000)
 9000 FORMAT ('1ERROR    1 IN I1MACH - I OUT OF BOUNDS')
C
C     CALL FDUMP
C
C
      STOP
      END
*DECK R1MACH
      REAL FUNCTION R1MACH(I)
C***BEGIN PROLOGUE  R1MACH
C***DATE WRITTEN   790101   (YYMMDD)
C***REVISION DATE  890213   (YYMMDD)
C***CATEGORY NO.  R1
C***KEYWORDS  LIBRARY=SLATEC,TYPE=SINGLE PRECISION(R1MACH-S D1MACH-D),
C             MACHINE CONSTANTS
C***AUTHOR  FOX, P. A., (BELL LABS)
C           HALL, A. D., (BELL LABS)
C           SCHRYER, N. L., (BELL LABS)
C***PURPOSE  Returns single precision machine dependent constants
C***DESCRIPTION
C
C   R1MACH can be used to obtain machine-dependent parameters
C   for the local machine environment.  It is a function
C   subroutine with one (input) argument, and can be called
C   as follows, for example
C
C        A = R1MACH(I)
C
C   where I=1,...,5.  The (output) value of A above is
C   determined by the (input) value of I.  The results for
C   various values of I are discussed below.
C
C   Single-Precision Machine Constants
C   R1MACH(1) = B**(EMIN-1), the smallest positive magnitude.
C   R1MACH(2) = B**EMAX*(1 - B**(-T)), the largest magnitude.
C   R1MACH(3) = B**(-T), the smallest relative spacing.
C   R1MACH(4) = B**(1-T), the largest relative spacing.
C   R1MACH(5) = LOG10(B)
C
C   Assume single precision numbers are represented in the T-digit,
C   base-B form
C
C              sign (B**E)*( (X(1)/B) + ... + (X(T)/B**T) )
C
C   where 0 .LE. X(I) .LT. B for I=1,...,T, 0 .LT. X(1), and
C   EMIN .LE. E .LE. EMAX.
C
C   The values of B, T, EMIN and EMAX are provided in I1MACH as
C   follows:
C   I1MACH(10) = B, the base.
C   I1MACH(11) = T, the number of base-B digits.
C   I1MACH(12) = EMIN, the smallest exponent E.
C   I1MACH(13) = EMAX, the largest exponent E.
C
C   To alter this function for a particular environment,
C   the desired set of DATA statements should be activated by
C   removing the C from column 1.  Also, the values of
C   R1MACH(1) - R1MACH(4) should be checked for consistency
C   with the local operating system.
C
C***REFERENCES  FOX, P.A., HALL, A.D., SCHRYER, N.L, *FRAMEWORK FOR
C                 A PORTABLE LIBRARY*, ACM TRANSACTIONS ON MATHE-
C                 MATICAL SOFTWARE, VOL. 4, NO. 2, JUNE 1978,
C                 PP. 177-188.
C***ROUTINES CALLED  XERROR
C***END PROLOGUE  R1MACH
C
      INTEGER SMALL(2)
      INTEGER LARGE(2)
      INTEGER RIGHT(2)
      INTEGER DIVER(2)
      INTEGER LOG10(2)
C
      REAL RMACH(5)
      SAVE RMACH
C
      EQUIVALENCE (RMACH(1),SMALL(1))
      EQUIVALENCE (RMACH(2),LARGE(1))
      EQUIVALENCE (RMACH(3),RIGHT(1))
      EQUIVALENCE (RMACH(4),DIVER(1))
      EQUIVALENCE (RMACH(5),LOG10(1))
C
C     MACHINE CONSTANTS FOR THE AMIGA
C     ABSOFT FORTRAN COMPILER USING THE 68020/68881 COMPILER OPTION
C
C     DATA SMALL(1) / Z'00800000' /
C     DATA LARGE(1) / Z'7F7FFFFF' /
C     DATA RIGHT(1) / Z'33800000' /
C     DATA DIVER(1) / Z'34000000' /
C     DATA LOG10(1) / Z'3E9A209B' /
C
C     MACHINE CONSTANTS FOR THE AMIGA
C     ABSOFT FORTRAN COMPILER USING SOFTWARE FLOATING POINT
C
C     DATA SMALL(1) / Z'00800000' /
C     DATA LARGE(1) / Z'7EFFFFFF' /
C     DATA RIGHT(1) / Z'33800000' /
C     DATA DIVER(1) / Z'34000000' /
C     DATA LOG10(1) / Z'3E9A209B' /
C
C     MACHINE CONSTANTS FOR THE APOLLO
C
C     DATA SMALL(1) / 16#00800000 /
C     DATA LARGE(1) / 16#7FFFFFFF /
C     DATA RIGHT(1) / 16#33800000 /
C     DATA DIVER(1) / 16#34000000 /
C     DATA LOG10(1) / 16#3E9A209B /
C
C     MACHINE CONSTANTS FOR THE BURROUGHS 1700 SYSTEM
C
C     DATA RMACH(1) / Z400800000 /
C     DATA RMACH(2) / Z5FFFFFFFF /
C     DATA RMACH(3) / Z4E9800000 /
C     DATA RMACH(4) / Z4EA800000 /
C     DATA RMACH(5) / Z500E730E8 /
C
C     MACHINE CONSTANTS FOR THE BURROUGHS 5700/6700/7700 SYSTEMS
C
C     DATA RMACH(1) / O1771000000000000 /
C     DATA RMACH(2) / O0777777777777777 /
C     DATA RMACH(3) / O1311000000000000 /
C     DATA RMACH(4) / O1301000000000000 /
C     DATA RMACH(5) / O1157163034761675 /
C
C     MACHINE CONSTANTS FOR THE CDC 170/180 SERIES USING NOS/VE
C
C     DATA RMACH(1) / Z"3001800000000000" /
C     DATA RMACH(2) / Z"4FFEFFFFFFFFFFFE" /
C     DATA RMACH(3) / Z"3FD2800000000000" /
C     DATA RMACH(4) / Z"3FD3800000000000" /
C     DATA RMACH(5) / Z"3FFF9A209A84FBCF" /
C
C     MACHINE CONSTANTS FOR THE CDC 6000/7000 SERIES
C
C     DATA RMACH(1) / 00564000000000000000B /
C     DATA RMACH(2) / 37767777777777777776B /
C     DATA RMACH(3) / 16414000000000000000B /
C     DATA RMACH(4) / 16424000000000000000B /
C     DATA RMACH(5) / 17164642023241175720B /
C
C     MACHINE CONSTANTS FOR THE CELERITY C1260
C
C     DATA SMALL(1) / Z'00800000' /
C     DATA LARGE(1) / Z'7F7FFFFF' /
C     DATA RIGHT(1) / Z'33800000' /
C     DATA DIVER(1) / Z'34000000' /
C     DATA LOG10(1) / Z'3E9A209B' /
C
C     MACHINE CONSTANTS FOR THE CONVEX C-1
C
C     DATA SMALL(1) / '00800000'X /
C     DATA LARGE(1) / '7FFFFFFF'X /
C     DATA RIGHT(1) / '34800000'X /
C     DATA DIVER(1) / '35000000'X /
C     DATA LOG10(1) / '3F9A209B'X /
C
C     MACHINE CONSTANTS FOR THE CRAY-1
C
C     DATA RMACH(1) / 200034000000000000000B /
C     DATA RMACH(2) / 577767777777777777776B /
C     DATA RMACH(3) / 377224000000000000000B /
C     DATA RMACH(4) / 377234000000000000000B /
C     DATA RMACH(5) / 377774642023241175720B /
C
C     MACHINE CONSTANTS FOR THE DATA GENERAL ECLIPSE S/200
C
C     NOTE - IT MAY BE APPROPRIATE TO INCLUDE THE FOLLOWING CARD -
C     STATIC RMACH(5)
C
C     DATA SMALL /    20K,       0 /
C     DATA LARGE / 77777K, 177777K /
C     DATA RIGHT / 35420K,       0 /
C     DATA DIVER / 36020K,       0 /
C     DATA LOG10 / 40423K,  42023K /
C
C     MACHINE CONSTANTS FOR THE HARRIS 220
C
C     DATA SMALL(1), SMALL(2) / '20000000, '00000201 /
C     DATA LARGE(1), LARGE(2) / '37777777, '00000177 /
C     DATA RIGHT(1), RIGHT(2) / '20000000, '00000352 /
C     DATA DIVER(1), DIVER(2) / '20000000, '00000353 /
C     DATA LOG10(1), LOG10(2) / '23210115, '00000377 /
C
C     MACHINE CONSTANTS FOR THE HONEYWELL 600/6000 SERIES
C
C     DATA RMACH(1) / O402400000000 /
C     DATA RMACH(2) / O376777777777 /
C     DATA RMACH(3) / O714400000000 /
C     DATA RMACH(4) / O716400000000 /
C     DATA RMACH(5) / O776464202324 /
C
C     MACHINE CONSTANTS FOR THE HP 2100
C     3 WORD DOUBLE PRECISION WITH FTN4
C
C     DATA SMALL(1), SMALL(2) / 40000B,       1 /
C     DATA LARGE(1), LARGE(2) / 77777B, 177776B /
C     DATA RIGHT(1), RIGHT(2) / 40000B,    325B /
C     DATA DIVER(1), DIVER(2) / 40000B,    327B /
C     DATA LOG10(1), LOG10(2) / 46420B,  46777B /
C
C     MACHINE CONSTANTS FOR THE HP 2100
C     4 WORD DOUBLE PRECISION WITH FTN4
C
C     DATA SMALL(1), SMALL(2) / 40000B,       1 /
C     DATA LARGE91), LARGE(2) / 77777B, 177776B /
C     DATA RIGHT(1), RIGHT(2) / 40000B,    325B /
C     DATA DIVER(1), DIVER(2) / 40000B,    327B /
C     DATA LOG10(1), LOG10(2) / 46420B,  46777B /
C
C     MACHINE CONSTANTS FOR THE HP 9000
C
C     DATA SMALL(1) / 00004000000B /
C     DATA LARGE(1) / 17677777777B /
C     DATA RIGHT(1) / 06340000000B /
C     DATA DIVER(1) / 06400000000B /
C     DATA LOG10(1) / 07646420233B /
C
C     MACHINE CONSTANTS FOR THE ELXSI 6400
C       ASSUMING REAL*4 IS THE DEFAULT REAL
C
C     DATA SMALL(1) / '00800000'X /
C     DATA LARGE(1) / '7F7FFFFF'X /
C     DATA RIGHT(1) / '33800000'X /
C     DATA DIVER(1) / '34000000'X /
C     DATA LOG10(1) / '3E9A209B'X /
C
C     MACHINE CONSTANTS FOR THE IBM 360/370 SERIES,
C     THE XEROX SIGMA 5/7/9, THE SEL SYSTEMS 85/86  AND
C     THE PERKIN ELMER (INTERDATA) 7/32.
C
C     DATA RMACH(1) / Z00100000 /
C     DATA RMACH(2) / Z7FFFFFFF /
C     DATA RMACH(3) / Z3B100000 /
C     DATA RMACH(4) / Z3C100000 /
C     DATA RMACH(5) / Z41134413 /
C
C     MACHINE CONSTANTS FOR THE IBM PC
C
C     DATA SMALL(1) /     8420761 /
C     DATA LARGE(1) /  2139081118 /
C     DATA RIGHT(1) /   863997169 /
C     DATA DIVER(1) /   872385777 /
C     DATA LOG10(1) /  1050288283 /
C
C     MACHINE CONSTANTS FOR THE IBM RS 6000
C
C     DATA SMALL(1) / Z'00800000' /
C     DATA LARGE(1) / Z'7F7FFFFF' /
C     DATA RIGHT(1) / Z'33800000' /
C     DATA DIVER(1) / Z'34000000' /
C     DATA LOG10(1) / Z'3E9A209B' /
C
C     MACHINE CONSTANTS FOR THE PDP-10 (KA OR KI PROCESSOR)
C
C     DATA RMACH(1) / "000400000000 /
C     DATA RMACH(2) / "377777777777 /
C     DATA RMACH(3) / "146400000000 /
C     DATA RMACH(4) / "147400000000 /
C     DATA RMACH(5) / "177464202324 /
C
C     MACHINE CONSTANTS FOR PDP-11 FORTRAN SUPPORTING
C     32-BIT INTEGERS (EXPRESSED IN INTEGER AND OCTAL).
C
C     DATA SMALL(1) /    8388608 /
C     DATA LARGE(1) / 2147483647 /
C     DATA RIGHT(1) /  880803840 /
C     DATA DIVER(1) /  889192448 /
C     DATA LOG10(1) / 1067065499 /
C
C     DATA RMACH(1) / O00040000000 /
C     DATA RMACH(2) / O17777777777 /
C     DATA RMACH(3) / O06440000000 /
C     DATA RMACH(4) / O06500000000 /
C     DATA RMACH(5) / O07746420233 /
C
C     MACHINE CONSTANTS FOR PDP-11 FORTRAN SUPPORTING
C     16-BIT INTEGERS  (EXPRESSED IN INTEGER AND OCTAL).
C
C     DATA SMALL(1), SMALL(2) /   128,     0 /
C     DATA LARGE(1), LARGE(2) / 32767,    -1 /
C     DATA RIGHT(1), RIGHT(2) / 13440,     0 /
C     DATA DIVER(1), DIVER(2) / 13568,     0 /
C     DATA LOG10(1), LOG10(2) / 16282,  8347 /
C
C     DATA SMALL(1), SMALL(2) / O000200, O000000 /
C     DATA LARGE(1), LARGE(2) / O077777, O177777 /
C     DATA RIGHT(1), RIGHT(2) / O032200, O000000 /
C     DATA DIVER(1), DIVER(2) / O032400, O000000 /
C     DATA LOG10(1), LOG10(2) / O037632, O020233 /
C
C     MACHINE CONSTANTS FOR THE SILICON GRAPHICS IRIS
C
c     data rmach(1) / 1.17549 424 e-38 /
c     data rmach(2) / 3.40282 356 e+38 /
c     data rmach(3) / 1.19209 290 e-07 /
c     data rmach(4) / 2.38418 579 e-07 /
c     data rmach(5) / 0.30103 001 /
C
C     DATA SMALL(1) / Z'00800000' /
C     DATA LARGE(1) / Z'7F7FFFFF' /
C     DATA RIGHT(1) / Z'34000000' /
C     DATA DIVER(1) / Z'34800000' /
C     DATA LOG10(1) / Z'3E9A209B' /
C
C     MACHINE CONSTANTS FOR THE SUN
C
C     DATA SMALL(1) / Z'00800000' /
C     DATA LARGE(1) / Z'7F7FFFFF' /
C     DATA RIGHT(1) / Z'33800000' /
C     DATA DIVER(1) / Z'34000000' /
C     DATA LOG10(1) / Z'3E9A209B' /
C
C     MACHINE CONSTANTS FOR THE UNIVAC 1100 SERIES
C
C     DATA RMACH(1) / O000400000000 /
C     DATA RMACH(2) / O377777777777 /
C     DATA RMACH(3) / O146400000000 /
C     DATA RMACH(4) / O147400000000 /
C     DATA RMACH(5) / O177464202324 /
C
C     MACHINE CONSTANTS FOR THE VAX 11/780
C     (EXPRESSED IN INTEGER AND HEXADECIMAL)
C     THE HEX FORMAT BELOW MAY NOT BE SUITABLE FOR UNIX SYSTEMS
C     THE INTEGER FORMAT SHOULD BE OK FOR UNIX SYSTEMS
C
C     DATA SMALL(1) /       128 /
C     DATA LARGE(1) /    -32769 /
C     DATA RIGHT(1) /     13440 /
C     DATA DIVER(1) /     13568 /
C     DATA LOG10(1) / 547045274 /
C
C     DATA SMALL(1) / Z00000080 /
C     DATA LARGE(1) / ZFFFF7FFF /
C     DATA RIGHT(1) / Z00003480 /
C     DATA DIVER(1) / Z00003500 /
C     DATA LOG10(1) / Z209B3F9A /
C
C     MACHINE CONSTANTS FOR THE Z80 MICROPROCESSOR
C
C     DATA SMALL(1), SMALL(2) /     0,    256/
C     DATA LARGE(1), LARGE(2) /    -1,   -129/
C     DATA RIGHT(1), RIGHT(2) /     0,  26880/
C     DATA DIVER(1), DIVER(2) /     0,  27136/
C     DATA LOG10(1), LOG10(2) /  8347,  32538/
C
C
C***FIRST EXECUTABLE STATEMENT  R1MACH
      IF (I .LT. 1  .OR.  I .GT. 5)
     1   CALL XERROR ('R1MACH -- I OUT OF BOUNDS', 25, 1, 2)
C
      R1MACH = RMACH(I)
      RETURN
C
      END
*DECK D1MACH
      DOUBLE PRECISION FUNCTION D1MACH(I)
C***BEGIN PROLOGUE  D1MACH
C***DATE WRITTEN   750101   (YYMMDD)
C***REVISION DATE  890213   (YYMMDD)
C***CATEGORY NO.  R1
C***KEYWORDS  LIBRARY=SLATEC,TYPE=DOUBLE PRECISION(R1MACH-S D1MACH-D),
C             MACHINE CONSTANTS
C***AUTHOR  FOX, P. A., (BELL LABS)
C           HALL, A. D., (BELL LABS)
C           SCHRYER, N. L., (BELL LABS)
C***PURPOSE  Returns double precision machine dependent constants
C***DESCRIPTION
C
C   D1MACH can be used to obtain machine-dependent parameters
C   for the local machine environment.  It is a function
C   subprogram with one (input) argument, and can be called
C   as follows, for example
C
C        D = D1MACH(I)
C
C   where I=1,...,5.  The (output) value of D above is
C   determined by the (input) value of I.  The results for
C   various values of I are discussed below.
C
C   D1MACH( 1) = B**(EMIN-1), the smallest positive magnitude.
C   D1MACH( 2) = B**EMAX*(1 - B**(-T)), the largest magnitude.
C   D1MACH( 3) = B**(-T), the smallest relative spacing.
C   D1MACH( 4) = B**(1-T), the largest relative spacing.
C   D1MACH( 5) = LOG10(B)
C
C   Assume double precision numbers are represented in the T-digit,
C   base-B form
C
C              sign (B**E)*( (X(1)/B) + ... + (X(T)/B**T) )
C
C   where 0 .LE. X(I) .LT. B for I=1,...,T, 0 .LT. X(1), and
C   EMIN .LE. E .LE. EMAX.
C
C   The values of B, T, EMIN and EMAX are provided in I1MACH as
C   follows:
C   I1MACH(10) = B, the base.
C   I1MACH(14) = T, the number of base-B digits.
C   I1MACH(15) = EMIN, the smallest exponent E.
C   I1MACH(16) = EMAX, the largest exponent E.
C
C   To alter this function for a particular environment,
C   the desired set of DATA statements should be activated by
C   removing the C from column 1.  Also, the values of
C   D1MACH(1) - D1MACH(4) should be checked for consistency
C   with the local operating system.
C
C***REFERENCES  FOX P.A., HALL A.D., SCHRYER N.L.,*FRAMEWORK FOR A
C                 PORTABLE LIBRARY*, ACM TRANSACTIONS ON MATHEMATICAL
C                 SOFTWARE, VOL. 4, NO. 2, JUNE 1978, PP. 177-188.
C***ROUTINES CALLED  XERROR
C***END PROLOGUE  D1MACH
C
      INTEGER SMALL(4)
      INTEGER LARGE(4)
      INTEGER RIGHT(4)
      INTEGER DIVER(4)
      INTEGER LOG10(4)
C
      DOUBLE PRECISION DMACH(5)
      SAVE DMACH
C
      EQUIVALENCE (DMACH(1),SMALL(1))
      EQUIVALENCE (DMACH(2),LARGE(1))
      EQUIVALENCE (DMACH(3),RIGHT(1))
      EQUIVALENCE (DMACH(4),DIVER(1))
      EQUIVALENCE (DMACH(5),LOG10(1))
C
C     MACHINE CONSTANTS FOR THE AMIGA
C     ABSOFT FORTRAN COMPILER USING THE 68020/68881 COMPILER OPTION
C
C     DATA SMALL(1), SMALL(2) / Z'00100000', Z'00000000' /
C     DATA LARGE(1), LARGE(2) / Z'7FEFFFFF', Z'FFFFFFFF' /
C     DATA RIGHT(1), RIGHT(2) / Z'3CA00000', Z'00000000' /
C     DATA DIVER(1), DIVER(2) / Z'3CB00000', Z'00000000' /
C     DATA LOG10(1), LOG10(2) / Z'3FD34413', Z'509F79FF' /
C
C     MACHINE CONSTANTS FOR THE AMIGA
C     ABSOFT FORTRAN COMPILER USING SOFTWARE FLOATING POINT
C
C     DATA SMALL(1), SMALL(2) / Z'00100000', Z'00000000' /
C     DATA LARGE(1), LARGE(2) / Z'7FDFFFFF', Z'FFFFFFFF' /
C     DATA RIGHT(1), RIGHT(2) / Z'3CA00000', Z'00000000' /
C     DATA DIVER(1), DIVER(2) / Z'3CB00000', Z'00000000' /
C     DATA LOG10(1), LOG10(2) / Z'3FD34413', Z'509F79FF' /
C
C     MACHINE CONSTANTS FOR THE APOLLO
C
C     DATA SMALL(1), SMALL(2) / 16#00100000, 16#00000000 /
C     DATA LARGE(1), LARGE(2) / 16#7FFFFFFF, 16#FFFFFFFF /
C     DATA RIGHT(1), RIGHT(2) / 16#3CA00000, 16#00000000 /
C     DATA DIVER(1), DIVER(2) / 16#3CB00000, 16#00000000 /
C     DATA LOG10(1), LOG10(2) / 16#3FD34413, 16#509F79FF /
C
C     MACHINE CONSTANTS FOR THE BURROUGHS 1700 SYSTEM
C
C     DATA SMALL(1) / ZC00800000 /
C     DATA SMALL(2) / Z000000000 /
C     DATA LARGE(1) / ZDFFFFFFFF /
C     DATA LARGE(2) / ZFFFFFFFFF /
C     DATA RIGHT(1) / ZCC5800000 /
C     DATA RIGHT(2) / Z000000000 /
C     DATA DIVER(1) / ZCC6800000 /
C     DATA DIVER(2) / Z000000000 /
C     DATA LOG10(1) / ZD00E730E7 /
C     DATA LOG10(2) / ZC77800DC0 /
C
C     MACHINE CONSTANTS FOR THE BURROUGHS 5700 SYSTEM
C
C     DATA SMALL(1) / O1771000000000000 /
C     DATA SMALL(2) / O0000000000000000 /
C     DATA LARGE(1) / O0777777777777777 /
C     DATA LARGE(2) / O0007777777777777 /
C     DATA RIGHT(1) / O1461000000000000 /
C     DATA RIGHT(2) / O0000000000000000 /
C     DATA DIVER(1) / O1451000000000000 /
C     DATA DIVER(2) / O0000000000000000 /
C     DATA LOG10(1) / O1157163034761674 /
C     DATA LOG10(2) / O0006677466732724 /
C
C     MACHINE CONSTANTS FOR THE BURROUGHS 6700/7700 SYSTEMS
C
C     DATA SMALL(1) / O1771000000000000 /
C     DATA SMALL(2) / O7770000000000000 /
C     DATA LARGE(1) / O0777777777777777 /
C     DATA LARGE(2) / O7777777777777777 /
C     DATA RIGHT(1) / O1461000000000000 /
C     DATA RIGHT(2) / O0000000000000000 /
C     DATA DIVER(1) / O1451000000000000 /
C     DATA DIVER(2) / O0000000000000000 /
C     DATA LOG10(1) / O1157163034761674 /
C     DATA LOG10(2) / O0006677466732724 /
C
C     MACHINE CONSTANTS FOR THE CDC 170/180 SERIES USING NOS/VE
C
C     DATA SMALL(1) / Z"3001800000000000" /
C     DATA SMALL(2) / Z"3001000000000000" /
C     DATA LARGE(1) / Z"4FFEFFFFFFFFFFFE" /
C     DATA LARGE(2) / Z"4FFE000000000000" /
C     DATA RIGHT(1) / Z"3FD2800000000000" /
C     DATA RIGHT(2) / Z"3FD2000000000000" /
C     DATA DIVER(1) / Z"3FD3800000000000" /
C     DATA DIVER(2) / Z"3FD3000000000000" /
C     DATA LOG10(1) / Z"3FFF9A209A84FBCF" /
C     DATA LOG10(2) / Z"3FFFF7988F8959AC" /
C
C     MACHINE CONSTANTS FOR THE CDC 6000/7000 SERIES
C
C     DATA SMALL(1) / 00564000000000000000B /
C     DATA SMALL(2) / 00000000000000000000B /
C     DATA LARGE(1) / 37757777777777777777B /
C     DATA LARGE(2) / 37157777777777777777B /
C     DATA RIGHT(1) / 15624000000000000000B /
C     DATA RIGHT(2) / 00000000000000000000B /
C     DATA DIVER(1) / 15634000000000000000B /
C     DATA DIVER(2) / 00000000000000000000B /
C     DATA LOG10(1) / 17164642023241175717B /
C     DATA LOG10(2) / 16367571421742254654B /
C
C     MACHINE CONSTANTS FOR THE CELERITY C1260
C
C     DATA SMALL(1), SMALL(2) / Z'00100000', Z'00000000' /
C     DATA LARGE(1), LARGE(2) / Z'7FEFFFFF', Z'FFFFFFFF' /
C     DATA RIGHT(1), RIGHT(2) / Z'3CA00000', Z'00000000' /
C     DATA DIVER(1), DIVER(2) / Z'3CB00000', Z'00000000' /
C     DATA LOG10(1), LOG10(2) / Z'3FD34413', Z'509F79FF' /
C
C     MACHINE CONSTANTS FOR THE CONVEX C-1
C
C     DATA SMALL(1), SMALL(2) / '00100000'X,'00000000'X /
C     DATA LARGE(1), LARGE(2) / '7FFFFFFF'X,'FFFFFFFF'X /
C     DATA RIGHT(1), RIGHT(2) / '3CC00000'X,'00000000'X /
C     DATA DIVER(1), DIVER(2) / '3CD00000'X,'00000000'X /
C     DATA LOG10(1), LOG10(2) / '3FF34413'X,'509F79FF'X /
C
C     MACHINE CONSTANTS FOR THE CRAY-1
C
C     DATA SMALL(1) / 201354000000000000000B /
C     DATA SMALL(2) / 000000000000000000000B /
C     DATA LARGE(1) / 577767777777777777777B /
C     DATA LARGE(2) / 000007777777777777774B /
C     DATA RIGHT(1) / 376434000000000000000B /
C     DATA RIGHT(2) / 000000000000000000000B /
C     DATA DIVER(1) / 376444000000000000000B /
C     DATA DIVER(2) / 000000000000000000000B /
C     DATA LOG10(1) / 377774642023241175717B /
C     DATA LOG10(2) / 000007571421742254654B /
C
C     MACHINE CONSTANTS FOR THE DATA GENERAL ECLIPSE S/200
C
C     NOTE - IT MAY BE APPROPRIATE TO INCLUDE THE FOLLOWING CARD -
C     STATIC DMACH(5)
C
C     DATA SMALL /    20K, 3*0 /
C     DATA LARGE / 77777K, 3*177777K /
C     DATA RIGHT / 31420K, 3*0 /
C     DATA DIVER / 32020K, 3*0 /
C     DATA LOG10 / 40423K, 42023K, 50237K, 74776K /
C
C     MACHINE CONSTANTS FOR THE ELXSI 6400
C     (ASSUMING REAL*8 IS THE DEFAULT DOUBLE PRECISION)
C
C     DATA SMALL(1), SMALL(2) / '00100000'X,'00000000'X /
C     DATA LARGE(1), LARGE(2) / '7FEFFFFF'X,'FFFFFFFF'X /
C     DATA RIGHT(1), RIGHT(2) / '3CB00000'X,'00000000'X /
C     DATA DIVER(1), DIVER(2) / '3CC00000'X,'00000000'X /
C     DATA LOG10(1), LOG10(2) / '3FD34413'X,'509F79FF'X /
C
C     MACHINE CONSTANTS FOR THE HARRIS 220
C
C     DATA SMALL(1), SMALL(2) / '20000000, '00000201 /
C     DATA LARGE(1), LARGE(2) / '37777777, '37777577 /
C     DATA RIGHT(1), RIGHT(2) / '20000000, '00000333 /
C     DATA DIVER(1), DIVER(2) / '20000000, '00000334 /
C     DATA LOG10(1), LOG10(2) / '23210115, '10237777 /
C
C     MACHINE CONSTANTS FOR THE HONEYWELL 600/6000 SERIES
C
C     DATA SMALL(1), SMALL(2) / O402400000000, O000000000000 /
C     DATA LARGE(1), LARGE(2) / O376777777777, O777777777777 /
C     DATA RIGHT(1), RIGHT(2) / O604400000000, O000000000000 /
C     DATA DIVER(1), DIVER(2) / O606400000000, O000000000000 /
C     DATA LOG10(1), LOG10(2) / O776464202324, O117571775714 /
C
C     MACHINE CONSTANTS FOR THE HP 2100
C     THREE WORD DOUBLE PRECISION OPTION WITH FTN4
C
C     DATA SMALL(1), SMALL(2), SMALL(3) / 40000B,       0,       1 /
C     DATA LARGE(1), LARGE(2), LARGE(3) / 77777B, 177777B, 177776B /
C     DATA RIGHT(1), RIGHT(2), RIGHT(3) / 40000B,       0,    265B /
C     DATA DIVER(1), DIVER(2), DIVER(3) / 40000B,       0,    276B /
C     DATA LOG10(1), LOG10(2), LOG10(3) / 46420B,  46502B,  77777B /
C
C     MACHINE CONSTANTS FOR THE HP 2100
C     FOUR WORD DOUBLE PRECISION OPTION WITH FTN4
C
C     DATA SMALL(1), SMALL(2) /  40000B,       0 /
C     DATA SMALL(3), SMALL(4) /       0,       1 /
C     DATA LARGE(1), LARGE(2) /  77777B, 177777B /
C     DATA LARGE(3), LARGE(4) / 177777B, 177776B /
C     DATA RIGHT(1), RIGHT(2) /  40000B,       0 /
C     DATA RIGHT(3), RIGHT(4) /       0,    225B /
C     DATA DIVER(1), DIVER(2) /  40000B,       0 /
C     DATA DIVER(3), DIVER(4) /       0,    227B /
C     DATA LOG10(1), LOG10(2) /  46420B,  46502B /
C     DATA LOG10(3), LOG10(4) /  76747B, 176377B /
C
C     MACHINE CONSTANTS FOR THE HP 9000
C
C     DATA SMALL(1), SMALL(2) / 00040000000B, 00000000000B /
C     DATA LARGE(1), LARGE(2) / 17737777777B, 37777777777B /
C     DATA RIGHT(1), RIGHT(2) / 07454000000B, 00000000000B /
C     DATA DIVER(1), DIVER(2) / 07460000000B, 00000000000B /
C     DATA LOG10(1), LOG10(2) / 07764642023B, 12047674777B /
C
C     MACHINE CONSTANTS FOR THE IBM 360/370 SERIES,
C     THE XEROX SIGMA 5/7/9, THE SEL SYSTEMS 85/86, AND
C     THE PERKIN ELMER (INTERDATA) 7/32.
C
C     DATA SMALL(1), SMALL(2) / Z00100000, Z00000000 /
C     DATA LARGE(1), LARGE(2) / Z7FFFFFFF, ZFFFFFFFF /
C     DATA RIGHT(1), RIGHT(2) / Z33100000, Z00000000 /
C     DATA DIVER(1), DIVER(2) / Z34100000, Z00000000 /
C     DATA LOG10(1), LOG10(2) / Z41134413, Z509F79FF /
C
C     MACHINE CONSTANTS FOR THE IBM PC
C     ASSUMES THAT ALL ARITHMETIC IS DONE IN DOUBLE PRECISION
C     ON 8088, I.E., NOT IN 80 BIT FORM FOR THE 8087.
C
C     DATA SMALL(1),SMALL(2) /  2002288515,    1050897 /
C     DATA LARGE(1),LARGE(2) /  1487780761, 2146426097 /
C     DATA RIGHT(1),RIGHT(2) / -1209488034, 1017118298 /
C     DATA DIVER(1),DIVER(2) / -1209488034, 1018166874 /
C     DATA LOG10(1),LOG10(2) /  1352628735, 1070810131 /
C
C     MACHINE CONSTANTS FOR THE IBM RS 6000
C
C     DATA SMALL(1), SMALL(2) / Z'00100000', Z'00000000' /
C     DATA LARGE(1), LARGE(2) / Z'7FEFFFFF', Z'FFFFFFFF' /
C     DATA RIGHT(1), RIGHT(2) / Z'3CA00000', Z'00000000' /
C     DATA DIVER(1), DIVER(2) / Z'3CB00000', Z'00000000' /
C     DATA LOG10(1), LOG10(2) / Z'3FD34413', Z'509F79FF' /
C
C     MACHINE CONSTANTS FOR THE PDP-10 (KA PROCESSOR)
C
C     DATA SMALL(1), SMALL(2) / "033400000000, "000000000000 /
C     DATA LARGE(1), LARGE(2) / "377777777777, "344777777777 /
C     DATA RIGHT(1), RIGHT(2) / "113400000000, "000000000000 /
C     DATA DIVER(1), DIVER(2) / "114400000000, "000000000000 /
C     DATA LOG10(1), LOG10(2) / "177464202324, "144117571776 /
C
C     MACHINE CONSTANTS FOR THE PDP-10 (KI PROCESSOR)
C
C     DATA SMALL(1), SMALL(2) / "000400000000, "000000000000 /
C     DATA LARGE(1), LARGE(2) / "377777777777, "377777777777 /
C     DATA RIGHT(1), RIGHT(2) / "103400000000, "000000000000 /
C     DATA DIVER(1), DIVER(2) / "104400000000, "000000000000 /
C     DATA LOG10(1), LOG10(2) / "177464202324, "476747767461 /
C
C     MACHINE CONSTANTS FOR PDP-11 FORTRAN SUPPORTING
C     32-BIT INTEGERS (EXPRESSED IN INTEGER AND OCTAL).
C
C     DATA SMALL(1), SMALL(2) /    8388608,           0 /
C     DATA LARGE(1), LARGE(2) / 2147483647,          -1 /
C     DATA RIGHT(1), RIGHT(2) /  612368384,           0 /
C     DATA DIVER(1), DIVER(2) /  620756992,           0 /
C     DATA LOG10(1), LOG10(2) / 1067065498, -2063872008 /
C
C     DATA SMALL(1), SMALL(2) / O00040000000, O00000000000 /
C     DATA LARGE(1), LARGE(2) / O17777777777, O37777777777 /
C     DATA RIGHT(1), RIGHT(2) / O04440000000, O00000000000 /
C     DATA DIVER(1), DIVER(2) / O04500000000, O00000000000 /
C     DATA LOG10(1), LOG10(2) / O07746420232, O20476747770 /
C
C     MACHINE CONSTANTS FOR PDP-11 FORTRAN SUPPORTING
C     16-BIT INTEGERS (EXPRESSED IN INTEGER AND OCTAL).
C
C     DATA SMALL(1), SMALL(2) /    128,      0 /
C     DATA SMALL(3), SMALL(4) /      0,      0 /
C     DATA LARGE(1), LARGE(2) /  32767,     -1 /
C     DATA LARGE(3), LARGE(4) /     -1,     -1 /
C     DATA RIGHT(1), RIGHT(2) /   9344,      0 /
C     DATA RIGHT(3), RIGHT(4) /      0,      0 /
C     DATA DIVER(1), DIVER(2) /   9472,      0 /
C     DATA DIVER(3), DIVER(4) /      0,      0 /
C     DATA LOG10(1), LOG10(2) /  16282,   8346 /
C     DATA LOG10(3), LOG10(4) / -31493, -12296 /
C
C     DATA SMALL(1), SMALL(2) / O000200, O000000 /
C     DATA SMALL(3), SMALL(4) / O000000, O000000 /
C     DATA LARGE(1), LARGE(2) / O077777, O177777 /
C     DATA LARGE(3), LARGE(4) / O177777, O177777 /
C     DATA RIGHT(1), RIGHT(2) / O022200, O000000 /
C     DATA RIGHT(3), RIGHT(4) / O000000, O000000 /
C     DATA DIVER(1), DIVER(2) / O022400, O000000 /
C     DATA DIVER(3), DIVER(4) / O000000, O000000 /
C     DATA LOG10(1), LOG10(2) / O037632, O020232 /
C     DATA LOG10(3), LOG10(4) / O102373, O147770 /
C
C     MACHINE CONSTANTS FOR THE SILICON GRAPHICS IRIS
C
c     data dmach(1) / 2.22507 38585 072012 d-308 /
c     data dmach(2) / 1.79769 31348 623158 d+308 /
c     data dmach(3) / 2.22044 60492 503131 d-16  /
c     data dmach(4) / 4.44089 20985 006262 d-16  /
c     data dmach(5) / 0.30102 99956 639812       /
C
C     DATA SMALL(1), SMALL(2) / Z'00100000',Z'00000000' /
C     DATA LARGE(1), LARGE(2) / Z'7FEFFFFF',Z'FFFFFFFF' /
C     DATA RIGHT(1), RIGHT(2) / Z'3CB00000',Z'00000000' /
C     DATA DIVER(1), DIVER(2) / Z'3CC00000',Z'00000000' /
C     DATA LOG10(1), LOG10(2) / Z'3FD34413',Z'509F79FF' /
C
C     MACHINE CONSTANTS FOR THE SUN
C
C     from SLATEC CML committee - work for Sun3, Sun4, and Sparc
C
C     DATA SMALL(1), SMALL(2) / Z'00100000', Z'00000000' /
C     DATA LARGE(1), LARGE(2) / Z'7FEFFFFF', Z'FFFFFFFF' /
C     DATA RIGHT(1), RIGHT(2) / Z'3CA00000', Z'00000000' /
C     DATA DIVER(1), DIVER(2) / Z'3CB00000', Z'00000000' /
C     DATA LOG10(1), LOG10(2) / Z'3FD34413', Z'509F79FF' /
C
C     from Sun Microsystems - work for Sun 386i
C
C     DATA SMALL(1), SMALL(2) / Z'00000000', Z'00100000' /
C     DATA LARGE(1), LARGE(2) / Z'FFFFFFFF', Z'7FEFFFFF' /
C     DATA RIGHT(1), RIGHT(2) / Z'00000000', Z'3CA00000' /
C     DATA DIVER(1), DIVER(2) / Z'00000000', Z'3CB00000' /
C     DATA LOG10(1), LOG10(2) / Z'509F79FF', Z'3FD34413' /
C
C     MACHINE CONSTANTS FOR THE UNIVAC 1100 SERIES FTN COMPILER
C
C     DATA SMALL(1), SMALL(2) / O000040000000, O000000000000 /
C     DATA LARGE(1), LARGE(2) / O377777777777, O777777777777 /
C     DATA RIGHT(1), RIGHT(2) / O170540000000, O000000000000 /
C     DATA DIVER(1), DIVER(2) / O170640000000, O000000000000 /
C     DATA LOG10(1), LOG10(2) / O177746420232, O411757177572 /
C
C     MACHINE CONSTANTS FOR VAX 11/780
C     (EXPRESSED IN INTEGER AND HEXADECIMAL)
C     THE HEX FORMAT BELOW MAY NOT BE SUITABLE FOR UNIX SYSYEMS
C     THE INTEGER FORMAT SHOULD BE OK FOR UNIX SYSTEMS
C
C     DATA SMALL(1), SMALL(2) /        128,           0 /
C     DATA LARGE(1), LARGE(2) /     -32769,          -1 /
C     DATA RIGHT(1), RIGHT(2) /       9344,           0 /
C     DATA DIVER(1), DIVER(2) /       9472,           0 /
C     DATA LOG10(1), LOG10(2) /  546979738,  -805796613 /
C
C     DATA SMALL(1), SMALL(2) / Z00000080, Z00000000 /
C     DATA LARGE(1), LARGE(2) / ZFFFF7FFF, ZFFFFFFFF /
C     DATA RIGHT(1), RIGHT(2) / Z00002480, Z00000000 /
C     DATA DIVER(1), DIVER(2) / Z00002500, Z00000000 /
C     DATA LOG10(1), LOG10(2) / Z209A3F9A, ZCFF884FB /
C
C     MACHINE CONSTANTS FOR VAX 11/780 (G-FLOATING)
C     (EXPRESSED IN INTEGER AND HEXADECIMAL)
C     THE HEX FORMAT BELOW MAY NOT BE SUITABLE FOR UNIX SYSYEMS
C     THE INTEGER FORMAT SHOULD BE OK FOR UNIX SYSTEMS
C
C     DATA SMALL(1), SMALL(2) /         16,           0 /
C     DATA LARGE(1), LARGE(2) /     -32769,          -1 /
C     DATA RIGHT(1), RIGHT(2) /      15552,           0 /
C     DATA DIVER(1), DIVER(2) /      15568,           0 /
C     DATA LOG10(1), LOG10(2) /  1142112243, 2046775455 /
C
C     DATA SMALL(1), SMALL(2) / Z00000010, Z00000000 /
C     DATA LARGE(1), LARGE(2) / ZFFFF7FFF, ZFFFFFFFF /
C     DATA RIGHT(1), RIGHT(2) / Z00003CC0, Z00000000 /
C     DATA DIVER(1), DIVER(2) / Z00003CD0, Z00000000 /
C     DATA LOG10(1), LOG10(2) / Z44133FF3, Z79FF509F /
C
C
C***FIRST EXECUTABLE STATEMENT  D1MACH
      IF (I .LT. 1  .OR.  I .GT. 5)
     1   CALL XERROR ('D1MACH -- I OUT OF BOUNDS', 25, 1, 2)
C
      D1MACH = DMACH(I)
      RETURN
C
      END
      SUBROUTINE FDUMP
C***BEGIN PROLOGUE  FDUMP
C***DATE WRITTEN   790801   (YYMMDD)
C***REVISION DATE  861211   (YYMMDD)
C***CATEGORY NO.  R3
C***KEYWORDS  LIBRARY=SLATEC(XERROR),TYPE=ALL(FDUMP-A),ERROR
C***AUTHOR  JONES, R. E., (SNLA)
C***PURPOSE  Symbolic dump (should be locally written).
C***DESCRIPTION
C
C        ***Note*** Machine Dependent Routine
C        FDUMP is intended to be replaced by a locally written
C        version which produces a symbolic dump.  Failing this,
C        it should be replaced by a version which prints the
C        subprogram nesting list.  Note that this dump must be
C        printed on each of up to five files, as indicated by the
C        XGETUA routine.  See XSETUA and XGETUA for details.
C
C     Written by Ron Jones, with SLATEC Common Math Library Subcommittee
C***REFERENCES  (NONE)
C***ROUTINES CALLED  (NONE)
C***END PROLOGUE  FDUMP
C***FIRST EXECUTABLE STATEMENT  FDUMP
      RETURN
      END
C*** zbsubs.f
-----------------------------------------------------------------
C>>>  ZBSUBS.FOR:  Double precision subroutines
-----------------------------------------------------------------

      SUBROUTINE ZBESH(ZR, ZI, FNU, KODE, M, N, CYR, CYI, NZ, IERR)
C***BEGIN PROLOGUE  ZBESH
C***DATE WRITTEN   830501   (YYMMDD)
C***REVISION DATE  890801, 930101   (YYMMDD)
C***CATEGORY NO.  B5K
C***KEYWORDS  H-BESSEL FUNCTIONS,BESSEL FUNCTIONS OF COMPLEX ARGUMENT,
C             BESSEL FUNCTIONS OF THIRD KIND,HANKEL FUNCTIONS
C***AUTHOR  AMOS, DONALD E., SANDIA NATIONAL LABORATORIES
C***PURPOSE  TO COMPUTE THE H-BESSEL FUNCTIONS OF A COMPLEX ARGUMENT
C***DESCRIPTION
C
C                      ***A DOUBLE PRECISION ROUTINE***
C         ON KODE=1, ZBESH COMPUTES AN N MEMBER SEQUENCE OF COMPLEX
C         HANKEL (BESSEL) FUNCTIONS CY(J)=H(M,FNU+J-1,Z) FOR KINDS M=1
C         OR 2, REAL, NONNEGATIVE ORDERS FNU+J-1, J=1,...,N, AND COMPLEX
C         Z.NE.CMPLX(0.0,0.0) IN THE CUT PLANE -PI.LT.ARG(Z).LE.PI.
C         ON KODE=2, ZBESH RETURNS THE SCALED HANKEL FUNCTIONS
C
C         CY(I)=EXP(-MM*Z*I)*H(M,FNU+J-1,Z)       MM=3-2*M,   I**2=-1.
C
C         WHICH REMOVES THE EXPONENTIAL BEHAVIOR IN BOTH THE UPPER AND
C         LOWER HALF PLANES. DEFINITIONS AND NOTATION ARE FOUND IN THE
C         NBS HANDBOOK OF MATHEMATICAL FUNCTIONS (REF. 1).
C
C         INPUT      ZR,ZI,FNU ARE DOUBLE PRECISION
C           ZR,ZI  - Z=CMPLX(ZR,ZI), Z.NE.CMPLX(0.0D0,0.0D0),
C                    -PT.LT.ARG(Z).LE.PI
C           FNU    - ORDER OF INITIAL H FUNCTION, FNU.GE.0.0D0
C           KODE   - A PARAMETER TO INDICATE THE SCALING OPTION
C                    KODE= 1  RETURNS
C                             CY(J)=H(M,FNU+J-1,Z),   J=1,...,N
C                        = 2  RETURNS
C                             CY(J)=H(M,FNU+J-1,Z)*EXP(-I*Z*(3-2M))
C                                  J=1,...,N  ,  I**2=-1
C           M      - KIND OF HANKEL FUNCTION, M=1 OR 2
C           N      - NUMBER OF MEMBERS IN THE SEQUENCE, N.GE.1
C
C         OUTPUT     CYR,CYI ARE DOUBLE PRECISION
C           CYR,CYI- DOUBLE PRECISION VECTORS WHOSE FIRST N COMPONENTS
C                    CONTAIN REAL AND IMAGINARY PARTS FOR THE SEQUENCE
C                    CY(J)=H(M,FNU+J-1,Z)  OR
C                    CY(J)=H(M,FNU+J-1,Z)*EXP(-I*Z*(3-2M))  J=1,...,N
C                    DEPENDING ON KODE, I**2=-1.
C           NZ     - NUMBER OF COMPONENTS SET TO ZERO DUE TO UNDERFLOW,
C                    NZ= 0   , NORMAL RETURN
C                    NZ.GT.0 , FIRST NZ COMPONENTS OF CY SET TO ZERO DUE
C                              TO UNDERFLOW, CY(J)=CMPLX(0.0D0,0.0D0)
C                              J=1,...,NZ WHEN Y.GT.0.0 AND M=1 OR
C                              Y.LT.0.0 AND M=2. FOR THE COMPLMENTARY
C                              HALF PLANES, NZ STATES ONLY THE NUMBER
C                              OF UNDERFLOWS.
C           IERR   - ERROR FLAG
C                    IERR=0, NORMAL RETURN - COMPUTATION COMPLETED
C                    IERR=1, INPUT ERROR   - NO COMPUTATION
C                    IERR=2, OVERFLOW      - NO COMPUTATION, FNU TOO
C                            LARGE OR CABS(Z) TOO SMALL OR BOTH
C                    IERR=3, CABS(Z) OR FNU+N-1 LARGE - COMPUTATION DONE
C                            BUT LOSSES OF SIGNIFCANCE BY ARGUMENT
C                            REDUCTION PRODUCE LESS THAN HALF OF MACHINE
C                            ACCURACY
C                    IERR=4, CABS(Z) OR FNU+N-1 TOO LARGE - NO COMPUTA-
C                            TION BECAUSE OF COMPLETE LOSSES OF SIGNIFI-
C                            CANCE BY ARGUMENT REDUCTION
C                    IERR=5, ERROR              - NO COMPUTATION,
C                            ALGORITHM TERMINATION CONDITION NOT MET
C
C***LONG DESCRIPTION
C
C         THE COMPUTATION IS CARRIED OUT BY THE RELATION
C
C         H(M,FNU,Z)=(1/MP)*EXP(-MP*FNU)*K(FNU,Z*EXP(-MP))
C             MP=MM*HPI*I,  MM=3-2*M,  HPI=PI/2,  I**2=-1
C
C         FOR M=1 OR 2 WHERE THE K BESSEL FUNCTION IS COMPUTED FOR THE
C         RIGHT HALF PLANE RE(Z).GE.0.0. THE K FUNCTION IS CONTINUED
C         TO THE LEFT HALF PLANE BY THE RELATION
C
C         K(FNU,Z*EXP(MP)) = EXP(-MP*FNU)*K(FNU,Z)-MP*I(FNU,Z)
C         MP=MR*PI*I, MR=+1 OR -1, RE(Z).GT.0, I**2=-1
C
C         WHERE I(FNU,Z) IS THE I BESSEL FUNCTION.
C
C         EXPONENTIAL DECAY OF H(M,FNU,Z) OCCURS IN THE UPPER HALF Z
C         PLANE FOR M=1 AND THE LOWER HALF Z PLANE FOR M=2.  EXPONENTIAL
C         GROWTH OCCURS IN THE COMPLEMENTARY HALF PLANES.  SCALING
C         BY EXP(-MM*Z*I) REMOVES THE EXPONENTIAL BEHAVIOR IN THE
C         WHOLE Z PLANE FOR Z TO INFINITY.
C
C         FOR NEGATIVE ORDERS,THE FORMULAE
C
C               H(1,-FNU,Z) = H(1,FNU,Z)*CEXP( PI*FNU*I)
C               H(2,-FNU,Z) = H(2,FNU,Z)*CEXP(-PI*FNU*I)
C                         I**2=-1
C
C         CAN BE USED.
C
C         IN MOST COMPLEX VARIABLE COMPUTATION, ONE MUST EVALUATE ELE-
C         MENTARY FUNCTIONS. WHEN THE MAGNITUDE OF Z OR FNU+N-1 IS
C         LARGE, LOSSES OF SIGNIFICANCE BY ARGUMENT REDUCTION OCCUR.
C         CONSEQUENTLY, IF EITHER ONE EXCEEDS U1=SQRT(0.5/UR), THEN
C         LOSSES EXCEEDING HALF PRECISION ARE LIKELY AND AN ERROR FLAG
C         IERR=3 IS TRIGGERED WHERE UR=DMAX1(D1MACH(4),1.0D-18) IS
C         DOUBLE PRECISION UNIT ROUNDOFF LIMITED TO 18 DIGITS PRECISION.
C         IF EITHER IS LARGER THAN U2=0.5/UR, THEN ALL SIGNIFICANCE IS
C         LOST AND IERR=4. IN ORDER TO USE THE INT FUNCTION, ARGUMENTS
C         MUST BE FURTHER RESTRICTED NOT TO EXCEED THE LARGEST MACHINE
C         INTEGER, U3=I1MACH(9). THUS, THE MAGNITUDE OF Z AND FNU+N-1 IS
C         RESTRICTED BY MIN(U2,U3). ON 32 BIT MACHINES, U1,U2, AND U3
C         ARE APPROXIMATELY 2.0E+3, 4.2E+6, 2.1E+9 IN SINGLE PRECISION
C         ARITHMETIC AND 1.3E+8, 1.8E+16, 2.1E+9 IN DOUBLE PRECISION
C         ARITHMETIC RESPECTIVELY. THIS MAKES U2 AND U3 LIMITING IN
C         THEIR RESPECTIVE ARITHMETICS. THIS MEANS THAT ONE CAN EXPECT
C         TO RETAIN, IN THE WORST CASES ON 32 BIT MACHINES, NO DIGITS
C         IN SINGLE AND ONLY 7 DIGITS IN DOUBLE PRECISION ARITHMETIC.
C         SIMILAR CONSIDERATIONS HOLD FOR OTHER MACHINES.
C
C         THE APPROXIMATE RELATIVE ERROR IN THE MAGNITUDE OF A COMPLEX
C         BESSEL FUNCTION CAN BE EXPRESSED BY P*10**S WHERE P=MAX(UNIT
C         ROUNDOFF,1.0D-18) IS THE NOMINAL PRECISION AND 10**S REPRE-
C         SENTS THE INCREASE IN ERROR DUE TO ARGUMENT REDUCTION IN THE
C         ELEMENTARY FUNCTIONS. HERE, S=MAX(1,ABS(LOG10(CABS(Z))),
C         ABS(LOG10(FNU))) APPROXIMATELY (I.E. S=MAX(1,ABS(EXPONENT OF
C         CABS(Z),ABS(EXPONENT OF FNU)) ). HOWEVER, THE PHASE ANGLE MAY
C         HAVE ONLY ABSOLUTE ACCURACY. THIS IS MOST LIKELY TO OCCUR WHEN
C         ONE COMPONENT (IN ABSOLUTE VALUE) IS LARGER THAN THE OTHER BY
C         SEVERAL ORDERS OF MAGNITUDE. IF ONE COMPONENT IS 10**K LARGER
C         THAN THE OTHER, THEN ONE CAN EXPECT ONLY MAX(ABS(LOG10(P))-K,
C         0) SIGNIFICANT DIGITS; OR, STATED ANOTHER WAY, WHEN K EXCEEDS
C         THE EXPONENT OF P, NO SIGNIFICANT DIGITS REMAIN IN THE SMALLER
C         COMPONENT. HOWEVER, THE PHASE ANGLE RETAINS ABSOLUTE ACCURACY
C         BECAUSE, IN COMPLEX ARITHMETIC WITH PRECISION P, THE SMALLER
C         COMPONENT WILL NOT (AS A RULE) DECREASE BELOW P TIMES THE
C         MAGNITUDE OF THE LARGER COMPONENT. IN THESE EXTREME CASES,
C         THE PRINCIPAL PHASE ANGLE IS ON THE ORDER OF +P, -P, PI/2-P,
C         OR -PI/2+P.
C
C***REFERENCES  HANDBOOK OF MATHEMATICAL FUNCTIONS BY M. ABRAMOWITZ
C                 AND I. A. STEGUN, NBS AMS SERIES 55, U.S. DEPT. OF
C                 COMMERCE, 1955.
C
C               COMPUTATION OF BESSEL FUNCTIONS OF COMPLEX ARGUMENT
C                 BY D. E. AMOS, SAND83-0083, MAY, 1983.
C
C               COMPUTATION OF BESSEL FUNCTIONS OF COMPLEX ARGUMENT
C                 AND LARGE ORDER BY D. E. AMOS, SAND83-0643, MAY, 1983
C
C               A SUBROUTINE PACKAGE FOR BESSEL FUNCTIONS OF A COMPLEX
C                 ARGUMENT AND NONNEGATIVE ORDER BY D. E. AMOS, SAND85-
C                 1018, MAY, 1985
C
C               A PORTABLE PACKAGE FOR BESSEL FUNCTIONS OF A COMPLEX
C                 ARGUMENT AND NONNEGATIVE ORDER BY D. E. AMOS, ACM
C                 TRANS. MATH. SOFTWARE, VOL. 12, NO. 3, SEPTEMBER 1986,
C                 PP 265-273.
C
C***ROUTINES CALLED  ZACON,ZBKNU,ZBUNK,ZUOIK,ZABS,I1MACH,D1MACH
C***END PROLOGUE  ZBESH
C
C     COMPLEX CY,Z,ZN,ZT,CSGN
      EXTERNAL ZABS
      DOUBLE PRECISION AA, ALIM, ALN, ARG, AZ, CYI, CYR, DIG, ELIM,
     * FMM, FN, FNU, FNUL, HPI, RHPI, RL, R1M5, SGN, STR, TOL, UFL, ZI,
     * ZNI, ZNR, ZR, ZTI, D1MACH, ZABS, BB, ASCLE, RTOL, ATOL, STI,
     * CSGNR, CSGNI
      INTEGER I, IERR, INU, INUH, IR, K, KODE, K1, K2, M,
     * MM, MR, N, NN, NUF, NW, NZ, I1MACH
      DIMENSION CYR(N), CYI(N)
C
      DATA HPI /1.57079632679489662D0/
C
C***FIRST EXECUTABLE STATEMENT  ZBESH
      IERR = 0
      NZ=0
      IF (ZR.EQ.0.0D0 .AND. ZI.EQ.0.0D0) IERR=1
      IF (FNU.LT.0.0D0) IERR=1
      IF (M.LT.1 .OR. M.GT.2) IERR=1
      IF (KODE.LT.1 .OR. KODE.GT.2) IERR=1
      IF (N.LT.1) IERR=1
      IF (IERR.NE.0) RETURN
      NN = N
C-----------------------------------------------------------------------
C     SET PARAMETERS RELATED TO MACHINE CONSTANTS.
C     TOL IS THE APPROXIMATE UNIT ROUNDOFF LIMITED TO 1.0E-18.
C     ELIM IS THE APPROXIMATE EXPONENTIAL OVER- AND UNDERFLOW LIMIT.
C     EXP(-ELIM).LT.EXP(-ALIM)=EXP(-ELIM)/TOL    AND
C     EXP(ELIM).GT.EXP(ALIM)=EXP(ELIM)*TOL       ARE INTERVALS NEAR
C     UNDERFLOW AND OVERFLOW LIMITS WHERE SCALED ARITHMETIC IS DONE.
C     RL IS THE LOWER BOUNDARY OF THE ASYMPTOTIC EXPANSION FOR LARGE Z.
C     DIG = NUMBER OF BASE 10 DIGITS IN TOL = 10**(-DIG).
C     FNUL IS THE LOWER BOUNDARY OF THE ASYMPTOTIC SERIES FOR LARGE FNU
C-----------------------------------------------------------------------
      TOL = DMAX1(D1MACH(4),1.0D-18)
      K1 = I1MACH(15)
      K2 = I1MACH(16)
      R1M5 = D1MACH(5)
      K = MIN0(IABS(K1),IABS(K2))
      ELIM = 2.303D0*(DBLE(FLOAT(K))*R1M5-3.0D0)
      K1 = I1MACH(14) - 1
      AA = R1M5*DBLE(FLOAT(K1))
      DIG = DMIN1(AA,18.0D0)
      AA = AA*2.303D0
      ALIM = ELIM + DMAX1(-AA,-41.45D0)
      FNUL = 10.0D0 + 6.0D0*(DIG-3.0D0)
      RL = 1.2D0*DIG + 3.0D0
      FN = FNU + DBLE(FLOAT(NN-1))
      MM = 3 - M - M
      FMM = DBLE(FLOAT(MM))
      ZNR = FMM*ZI
      ZNI = -FMM*ZR
C-----------------------------------------------------------------------
C     TEST FOR PROPER RANGE
C-----------------------------------------------------------------------
      AZ = ZABS(ZR,ZI)
      AA = 0.5D0/TOL
      BB=DBLE(FLOAT(I1MACH(9)))*0.5D0
      AA = DMIN1(AA,BB)
      IF (AZ.GT.AA) GO TO 260
      IF (FN.GT.AA) GO TO 260
      AA = DSQRT(AA)
      IF (AZ.GT.AA) IERR=3
      IF (FN.GT.AA) IERR=3
C-----------------------------------------------------------------------
C     OVERFLOW TEST ON THE LAST MEMBER OF THE SEQUENCE
C-----------------------------------------------------------------------
      UFL = D1MACH(1)*1.0D+3
      IF (AZ.LT.UFL) GO TO 230
      IF (FNU.GT.FNUL) GO TO 90
      IF (FN.LE.1.0D0) GO TO 70
      IF (FN.GT.2.0D0) GO TO 60
      IF (AZ.GT.TOL) GO TO 70
      ARG = 0.5D0*AZ
      ALN = -FN*DLOG(ARG)
      IF (ALN.GT.ELIM) GO TO 230
      GO TO 70
   60 CONTINUE
      CALL ZUOIK(ZNR, ZNI, FNU, KODE, 2, NN, CYR, CYI, NUF, TOL, ELIM,
     * ALIM)
      IF (NUF.LT.0) GO TO 230
      NZ = NZ + NUF
      NN = NN - NUF
C-----------------------------------------------------------------------
C     HERE NN=N OR NN=0 SINCE NUF=0,NN, OR -1 ON RETURN FROM CUOIK
C     IF NUF=NN, THEN CY(I)=CZERO FOR ALL I
C-----------------------------------------------------------------------
      IF (NN.EQ.0) GO TO 140
   70 CONTINUE
      IF ((ZNR.LT.0.0D0) .OR. (ZNR.EQ.0.0D0 .AND. ZNI.LT.0.0D0 .AND.
     * M.EQ.2)) GO TO 80
C-----------------------------------------------------------------------
C     RIGHT HALF PLANE COMPUTATION, XN.GE.0. .AND. (XN.NE.0. .OR.
C     YN.GE.0. .OR. M=1)
C-----------------------------------------------------------------------
      CALL ZBKNU(ZNR, ZNI, FNU, KODE, NN, CYR, CYI, NZ, TOL, ELIM, ALIM)
      GO TO 110
C-----------------------------------------------------------------------
C     LEFT HALF PLANE COMPUTATION
C-----------------------------------------------------------------------
   80 CONTINUE
      MR = -MM
      CALL ZACON(ZNR, ZNI, FNU, KODE, MR, NN, CYR, CYI, NW, RL, FNUL,
     * TOL, ELIM, ALIM)
      IF (NW.LT.0) GO TO 240
      NZ=NW
      GO TO 110
   90 CONTINUE
C-----------------------------------------------------------------------
C     UNIFORM ASYMPTOTIC EXPANSIONS FOR FNU.GT.FNUL
C-----------------------------------------------------------------------
      MR = 0
      IF ((ZNR.GE.0.0D0) .AND. (ZNR.NE.0.0D0 .OR. ZNI.GE.0.0D0 .OR.
     * M.NE.2)) GO TO 100
      MR = -MM
      IF (ZNR.NE.0.0D0 .OR. ZNI.GE.0.0D0) GO TO 100
      ZNR = -ZNR
      ZNI = -ZNI
  100 CONTINUE
      CALL ZBUNK(ZNR, ZNI, FNU, KODE, MR, NN, CYR, CYI, NW, TOL, ELIM,
     * ALIM)
      IF (NW.LT.0) GO TO 240
      NZ = NZ + NW
  110 CONTINUE
C-----------------------------------------------------------------------
C     H(M,FNU,Z) = -FMM*(I/HPI)*(ZT**FNU)*K(FNU,-Z*ZT)
C
C     ZT=EXP(-FMM*HPI*I) = CMPLX(0.0,-FMM), FMM=3-2*M, M=1,2
C-----------------------------------------------------------------------
      SGN = DSIGN(HPI,-FMM)
C-----------------------------------------------------------------------
C     CALCULATE EXP(FNU*HPI*I) TO MINIMIZE LOSSES OF SIGNIFICANCE
C     WHEN FNU IS LARGE
C-----------------------------------------------------------------------
      INU = INT(SNGL(FNU))
      INUH = INU/2
      IR = INU - 2*INUH
      ARG = (FNU-DBLE(FLOAT(INU-IR)))*SGN
      RHPI = 1.0D0/SGN
C     ZNI = RHPI*DCOS(ARG)
C     ZNR = -RHPI*DSIN(ARG)
      CSGNI = RHPI*DCOS(ARG)
      CSGNR = -RHPI*DSIN(ARG)
      IF (MOD(INUH,2).EQ.0) GO TO 120
C     ZNR = -ZNR
C     ZNI = -ZNI
      CSGNR = -CSGNR
      CSGNI = -CSGNI
  120 CONTINUE
      ZTI = -FMM
      RTOL = 1.0D0/TOL
      ASCLE = UFL*RTOL
      DO 130 I=1,NN
C       STR = CYR(I)*ZNR - CYI(I)*ZNI
C       CYI(I) = CYR(I)*ZNI + CYI(I)*ZNR
C       CYR(I) = STR
C       STR = -ZNI*ZTI
C       ZNI = ZNR*ZTI
C       ZNR = STR
        AA = CYR(I)
        BB = CYI(I)
        ATOL = 1.0D0
        IF (DMAX1(DABS(AA),DABS(BB)).GT.ASCLE) GO TO 135
          AA = AA*RTOL
          BB = BB*RTOL
          ATOL = TOL
  135 CONTINUE
      STR = AA*CSGNR - BB*CSGNI
      STI = AA*CSGNI + BB*CSGNR
      CYR(I) = STR*ATOL
      CYI(I) = STI*ATOL
      STR = -CSGNI*ZTI
      CSGNI = CSGNR*ZTI
      CSGNR = STR
  130 CONTINUE
      RETURN
  140 CONTINUE
      IF (ZNR.LT.0.0D0) GO TO 230
      RETURN
  230 CONTINUE
      NZ=0
      IERR=2
      RETURN
  240 CONTINUE
      IF(NW.EQ.(-1)) GO TO 230
      NZ=0
      IERR=5
      RETURN
  260 CONTINUE
      NZ=0
      IERR=4
      RETURN
      END
      SUBROUTINE ZBESI(ZR, ZI, FNU, KODE, N, CYR, CYI, NZ, IERR)
C***BEGIN PROLOGUE  ZBESI
C***DATE WRITTEN   830501   (YYMMDD)
C***REVISION DATE  890801, 930101   (YYMMDD)
C***CATEGORY NO.  B5K
C***KEYWORDS  I-BESSEL FUNCTION,COMPLEX BESSEL FUNCTION,
C             MODIFIED BESSEL FUNCTION OF THE FIRST KIND
C***AUTHOR  AMOS, DONALD E., SANDIA NATIONAL LABORATORIES
C***PURPOSE  TO COMPUTE I-BESSEL FUNCTIONS OF COMPLEX ARGUMENT
C***DESCRIPTION
C
C                    ***A DOUBLE PRECISION ROUTINE***
C         ON KODE=1, ZBESI COMPUTES AN N MEMBER SEQUENCE OF COMPLEX
C         BESSEL FUNCTIONS CY(J)=I(FNU+J-1,Z) FOR REAL, NONNEGATIVE
C         ORDERS FNU+J-1, J=1,...,N AND COMPLEX Z IN THE CUT PLANE
C         -PI.LT.ARG(Z).LE.PI. ON KODE=2, ZBESI RETURNS THE SCALED
C         FUNCTIONS
C
C         CY(J)=EXP(-ABS(X))*I(FNU+J-1,Z)   J = 1,...,N , X=REAL(Z)
C
C         WITH THE EXPONENTIAL GROWTH REMOVED IN BOTH THE LEFT AND
C         RIGHT HALF PLANES FOR Z TO INFINITY. DEFINITIONS AND NOTATION
C         ARE FOUND IN THE NBS HANDBOOK OF MATHEMATICAL FUNCTIONS
C         (REF. 1).
C
C         INPUT      ZR,ZI,FNU ARE DOUBLE PRECISION
C           ZR,ZI  - Z=CMPLX(ZR,ZI),  -PI.LT.ARG(Z).LE.PI
C           FNU    - ORDER OF INITIAL I FUNCTION, FNU.GE.0.0D0
C           KODE   - A PARAMETER TO INDICATE THE SCALING OPTION
C                    KODE= 1  RETURNS
C                             CY(J)=I(FNU+J-1,Z), J=1,...,N
C                        = 2  RETURNS
C                             CY(J)=I(FNU+J-1,Z)*EXP(-ABS(X)), J=1,...,N
C           N      - NUMBER OF MEMBERS OF THE SEQUENCE, N.GE.1
C
C         OUTPUT     CYR,CYI ARE DOUBLE PRECISION
C           CYR,CYI- DOUBLE PRECISION VECTORS WHOSE FIRST N COMPONENTS
C                    CONTAIN REAL AND IMAGINARY PARTS FOR THE SEQUENCE
C                    CY(J)=I(FNU+J-1,Z)  OR
C                    CY(J)=I(FNU+J-1,Z)*EXP(-ABS(X))  J=1,...,N
C                    DEPENDING ON KODE, X=REAL(Z)
C           NZ     - NUMBER OF COMPONENTS SET TO ZERO DUE TO UNDERFLOW,
C                    NZ= 0   , NORMAL RETURN
C                    NZ.GT.0 , LAST NZ COMPONENTS OF CY SET TO ZERO
C                              TO UNDERFLOW, CY(J)=CMPLX(0.0D0,0.0D0)
C                              J = N-NZ+1,...,N
C           IERR   - ERROR FLAG
C                    IERR=0, NORMAL RETURN - COMPUTATION COMPLETED
C                    IERR=1, INPUT ERROR   - NO COMPUTATION
C                    IERR=2, OVERFLOW      - NO COMPUTATION, REAL(Z) TOO
C                            LARGE ON KODE=1
C                    IERR=3, CABS(Z) OR FNU+N-1 LARGE - COMPUTATION DONE
C                            BUT LOSSES OF SIGNIFCANCE BY ARGUMENT
C                            REDUCTION PRODUCE LESS THAN HALF OF MACHINE
C                            ACCURACY
C                    IERR=4, CABS(Z) OR FNU+N-1 TOO LARGE - NO COMPUTA-
C                            TION BECAUSE OF COMPLETE LOSSES OF SIGNIFI-
C                            CANCE BY ARGUMENT REDUCTION
C                    IERR=5, ERROR              - NO COMPUTATION,
C                            ALGORITHM TERMINATION CONDITION NOT MET
C
C***LONG DESCRIPTION
C
C         THE COMPUTATION IS CARRIED OUT BY THE POWER SERIES FOR
C         SMALL CABS(Z), THE ASYMPTOTIC EXPANSION FOR LARGE CABS(Z),
C         THE MILLER ALGORITHM NORMALIZED BY THE WRONSKIAN AND A
C         NEUMANN SERIES FOR IMTERMEDIATE MAGNITUDES, AND THE
C         UNIFORM ASYMPTOTIC EXPANSIONS FOR I(FNU,Z) AND J(FNU,Z)
C         FOR LARGE ORDERS. BACKWARD RECURRENCE IS USED TO GENERATE
C         SEQUENCES OR REDUCE ORDERS WHEN NECESSARY.
C
C         THE CALCULATIONS ABOVE ARE DONE IN THE RIGHT HALF PLANE AND
C         CONTINUED INTO THE LEFT HALF PLANE BY THE FORMULA
C
C         I(FNU,Z*EXP(M*PI)) = EXP(M*PI*FNU)*I(FNU,Z)  REAL(Z).GT.0.0
C                       M = +I OR -I,  I**2=-1
C
C         FOR NEGATIVE ORDERS,THE FORMULA
C
C              I(-FNU,Z) = I(FNU,Z) + (2/PI)*SIN(PI*FNU)*K(FNU,Z)
C
C         CAN BE USED. HOWEVER,FOR LARGE ORDERS CLOSE TO INTEGERS, THE
C         THE FUNCTION CHANGES RADICALLY. WHEN FNU IS A LARGE POSITIVE
C         INTEGER,THE MAGNITUDE OF I(-FNU,Z)=I(FNU,Z) IS A LARGE
C         NEGATIVE POWER OF TEN. BUT WHEN FNU IS NOT AN INTEGER,
C         K(FNU,Z) DOMINATES IN MAGNITUDE WITH A LARGE POSITIVE POWER OF
C         TEN AND THE MOST THAT THE SECOND TERM CAN BE REDUCED IS BY
C         UNIT ROUNDOFF FROM THE COEFFICIENT. THUS, WIDE CHANGES CAN
C         OCCUR WITHIN UNIT ROUNDOFF OF A LARGE INTEGER FOR FNU. HERE,
C         LARGE MEANS FNU.GT.CABS(Z).
C
C         IN MOST COMPLEX VARIABLE COMPUTATION, ONE MUST EVALUATE ELE-
C         MENTARY FUNCTIONS. WHEN THE MAGNITUDE OF Z OR FNU+N-1 IS
C         LARGE, LOSSES OF SIGNIFICANCE BY ARGUMENT REDUCTION OCCUR.
C         CONSEQUENTLY, IF EITHER ONE EXCEEDS U1=SQRT(0.5/UR), THEN
C         LOSSES EXCEEDING HALF PRECISION ARE LIKELY AND AN ERROR FLAG
C         IERR=3 IS TRIGGERED WHERE UR=DMAX1(D1MACH(4),1.0D-18) IS
C         DOUBLE PRECISION UNIT ROUNDOFF LIMITED TO 18 DIGITS PRECISION.
C         IF EITHER IS LARGER THAN U2=0.5/UR, THEN ALL SIGNIFICANCE IS
C         LOST AND IERR=4. IN ORDER TO USE THE INT FUNCTION, ARGUMENTS
C         MUST BE FURTHER RESTRICTED NOT TO EXCEED THE LARGEST MACHINE
C         INTEGER, U3=I1MACH(9). THUS, THE MAGNITUDE OF Z AND FNU+N-1 IS
C         RESTRICTED BY MIN(U2,U3). ON 32 BIT MACHINES, U1,U2, AND U3
C         ARE APPROXIMATELY 2.0E+3, 4.2E+6, 2.1E+9 IN SINGLE PRECISION
C         ARITHMETIC AND 1.3E+8, 1.8E+16, 2.1E+9 IN DOUBLE PRECISION
C         ARITHMETIC RESPECTIVELY. THIS MAKES U2 AND U3 LIMITING IN
C         THEIR RESPECTIVE ARITHMETICS. THIS MEANS THAT ONE CAN EXPECT
C         TO RETAIN, IN THE WORST CASES ON 32 BIT MACHINES, NO DIGITS
C         IN SINGLE AND ONLY 7 DIGITS IN DOUBLE PRECISION ARITHMETIC.
C         SIMILAR CONSIDERATIONS HOLD FOR OTHER MACHINES.
C
C         THE APPROXIMATE RELATIVE ERROR IN THE MAGNITUDE OF A COMPLEX
C         BESSEL FUNCTION CAN BE EXPRESSED BY P*10**S WHERE P=MAX(UNIT
C         ROUNDOFF,1.0E-18) IS THE NOMINAL PRECISION AND 10**S REPRE-
C         SENTS THE INCREASE IN ERROR DUE TO ARGUMENT REDUCTION IN THE
C         ELEMENTARY FUNCTIONS. HERE, S=MAX(1,ABS(LOG10(CABS(Z))),
C         ABS(LOG10(FNU))) APPROXIMATELY (I.E. S=MAX(1,ABS(EXPONENT OF
C         CABS(Z),ABS(EXPONENT OF FNU)) ). HOWEVER, THE PHASE ANGLE MAY
C         HAVE ONLY ABSOLUTE ACCURACY. THIS IS MOST LIKELY TO OCCUR WHEN
C         ONE COMPONENT (IN ABSOLUTE VALUE) IS LARGER THAN THE OTHER BY
C         SEVERAL ORDERS OF MAGNITUDE. IF ONE COMPONENT IS 10**K LARGER
C         THAN THE OTHER, THEN ONE CAN EXPECT ONLY MAX(ABS(LOG10(P))-K,
C         0) SIGNIFICANT DIGITS; OR, STATED ANOTHER WAY, WHEN K EXCEEDS
C         THE EXPONENT OF P, NO SIGNIFICANT DIGITS REMAIN IN THE SMALLER
C         COMPONENT. HOWEVER, THE PHASE ANGLE RETAINS ABSOLUTE ACCURACY
C         BECAUSE, IN COMPLEX ARITHMETIC WITH PRECISION P, THE SMALLER
C         COMPONENT WILL NOT (AS A RULE) DECREASE BELOW P TIMES THE
C         MAGNITUDE OF THE LARGER COMPONENT. IN THESE EXTREME CASES,
C         THE PRINCIPAL PHASE ANGLE IS ON THE ORDER OF +P, -P, PI/2-P,
C         OR -PI/2+P.
C
C***REFERENCES  HANDBOOK OF MATHEMATICAL FUNCTIONS BY M. ABRAMOWITZ
C                 AND I. A. STEGUN, NBS AMS SERIES 55, U.S. DEPT. OF
C                 COMMERCE, 1955.
C
C               COMPUTATION OF BESSEL FUNCTIONS OF COMPLEX ARGUMENT
C                 BY D. E. AMOS, SAND83-0083, MAY, 1983.
C
C               COMPUTATION OF BESSEL FUNCTIONS OF COMPLEX ARGUMENT
C                 AND LARGE ORDER BY D. E. AMOS, SAND83-0643, MAY, 1983
C
C               A SUBROUTINE PACKAGE FOR BESSEL FUNCTIONS OF A COMPLEX
C                 ARGUMENT AND NONNEGATIVE ORDER BY D. E. AMOS, SAND85-
C                 1018, MAY, 1985
C
C               A PORTABLE PACKAGE FOR BESSEL FUNCTIONS OF A COMPLEX
C                 ARGUMENT AND NONNEGATIVE ORDER BY D. E. AMOS, ACM
C                 TRANS. MATH. SOFTWARE, VOL. 12, NO. 3, SEPTEMBER 1986,
C                 PP 265-273.
C
C***ROUTINES CALLED  ZBINU,ZABS,I1MACH,D1MACH
C***END PROLOGUE  ZBESI
C     COMPLEX CONE,CSGN,CW,CY,CZERO,Z,ZN
      EXTERNAL ZABS
      DOUBLE PRECISION AA, ALIM, ARG, CONEI, CONER, CSGNI, CSGNR, CYI,
     * CYR, DIG, ELIM, FNU, FNUL, PI, RL, R1M5, STR, TOL, ZI, ZNI, ZNR,
     * ZR, D1MACH, AZ, BB, FN, ZABS, ASCLE, RTOL, ATOL, STI
      INTEGER I, IERR, INU, K, KODE, K1,K2,N,NZ,NN, I1MACH
      DIMENSION CYR(N), CYI(N)
      DATA PI /3.14159265358979324D0/
      DATA CONER, CONEI /1.0D0,0.0D0/
C
C***FIRST EXECUTABLE STATEMENT  ZBESI
      IERR = 0
      NZ=0
      IF (FNU.LT.0.0D0) IERR=1
      IF (KODE.LT.1 .OR. KODE.GT.2) IERR=1
      IF (N.LT.1) IERR=1
      IF (IERR.NE.0) RETURN
C-----------------------------------------------------------------------
C     SET PARAMETERS RELATED TO MACHINE CONSTANTS.
C     TOL IS THE APPROXIMATE UNIT ROUNDOFF LIMITED TO 1.0E-18.
C     ELIM IS THE APPROXIMATE EXPONENTIAL OVER- AND UNDERFLOW LIMIT.
C     EXP(-ELIM).LT.EXP(-ALIM)=EXP(-ELIM)/TOL    AND
C     EXP(ELIM).GT.EXP(ALIM)=EXP(ELIM)*TOL       ARE INTERVALS NEAR
C     UNDERFLOW AND OVERFLOW LIMITS WHERE SCALED ARITHMETIC IS DONE.
C     RL IS THE LOWER BOUNDARY OF THE ASYMPTOTIC EXPANSION FOR LARGE Z.
C     DIG = NUMBER OF BASE 10 DIGITS IN TOL = 10**(-DIG).
C     FNUL IS THE LOWER BOUNDARY OF THE ASYMPTOTIC SERIES FOR LARGE FNU.
C-----------------------------------------------------------------------
      TOL = DMAX1(D1MACH(4),1.0D-18)
      K1 = I1MACH(15)
      K2 = I1MACH(16)
      R1M5 = D1MACH(5)
      K = MIN0(IABS(K1),IABS(K2))
      ELIM = 2.303D0*(DBLE(FLOAT(K))*R1M5-3.0D0)
      K1 = I1MACH(14) - 1
      AA = R1M5*DBLE(FLOAT(K1))
      DIG = DMIN1(AA,18.0D0)
      AA = AA*2.303D0
      ALIM = ELIM + DMAX1(-AA,-41.45D0)
      RL = 1.2D0*DIG + 3.0D0
      FNUL = 10.0D0 + 6.0D0*(DIG-3.0D0)
C-----------------------------------------------------------------------------
C     TEST FOR PROPER RANGE
C-----------------------------------------------------------------------
      AZ = ZABS(ZR,ZI)
      FN = FNU+DBLE(FLOAT(N-1))
      AA = 0.5D0/TOL
      BB=DBLE(FLOAT(I1MACH(9)))*0.5D0
      AA = DMIN1(AA,BB)
      IF (AZ.GT.AA) GO TO 260
      IF (FN.GT.AA) GO TO 260
      AA = DSQRT(AA)
      IF (AZ.GT.AA) IERR=3
      IF (FN.GT.AA) IERR=3
      ZNR = ZR
      ZNI = ZI
      CSGNR = CONER
      CSGNI = CONEI
      IF (ZR.GE.0.0D0) GO TO 40
      ZNR = -ZR
      ZNI = -ZI
C-----------------------------------------------------------------------
C     CALCULATE CSGN=EXP(FNU*PI*I) TO MINIMIZE LOSSES OF SIGNIFICANCE
C     WHEN FNU IS LARGE
C-----------------------------------------------------------------------
      INU = INT(SNGL(FNU))
      ARG = (FNU-DBLE(FLOAT(INU)))*PI
      IF (ZI.LT.0.0D0) ARG = -ARG
      CSGNR = DCOS(ARG)
      CSGNI = DSIN(ARG)
      IF (MOD(INU,2).EQ.0) GO TO 40
      CSGNR = -CSGNR
      CSGNI = -CSGNI
   40 CONTINUE
      CALL ZBINU(ZNR, ZNI, FNU, KODE, N, CYR, CYI, NZ, RL, FNUL, TOL,
     * ELIM, ALIM)
      IF (NZ.LT.0) GO TO 120
      IF (ZR.GE.0.0D0) RETURN
C-----------------------------------------------------------------------
C     ANALYTIC CONTINUATION TO THE LEFT HALF PLANE
C-----------------------------------------------------------------------
      NN = N - NZ
      IF (NN.EQ.0) RETURN
      RTOL = 1.0D0/TOL
      ASCLE = D1MACH(1)*RTOL*1.0D+3
      DO 50 I=1,NN
C       STR = CYR(I)*CSGNR - CYI(I)*CSGNI
C       CYI(I) = CYR(I)*CSGNI + CYI(I)*CSGNR
C       CYR(I) = STR
        AA = CYR(I)
        BB = CYI(I)
        ATOL = 1.0D0
        IF (DMAX1(DABS(AA),DABS(BB)).GT.ASCLE) GO TO 55
          AA = AA*RTOL
          BB = BB*RTOL
          ATOL = TOL
   55   CONTINUE
        STR = AA*CSGNR - BB*CSGNI
        STI = AA*CSGNI + BB*CSGNR
        CYR(I) = STR*ATOL
        CYI(I) = STI*ATOL
        CSGNR = -CSGNR
        CSGNI = -CSGNI
   50 CONTINUE
      RETURN
  120 CONTINUE
      IF(NZ.EQ.(-2)) GO TO 130
      NZ = 0
      IERR=2
      RETURN
  130 CONTINUE
      NZ=0
      IERR=5
      RETURN
  260 CONTINUE
      NZ=0
      IERR=4
      RETURN
      END
      SUBROUTINE ZBESJ(ZR, ZI, FNU, KODE, N, CYR, CYI, NZ, IERR)
C***BEGIN PROLOGUE  ZBESJ
C***DATE WRITTEN   830501   (YYMMDD)
C***REVISION DATE  890801, 930101   (YYMMDD)
C***CATEGORY NO.  B5K
C***KEYWORDS  J-BESSEL FUNCTION,BESSEL FUNCTION OF COMPLEX ARGUMENT,
C             BESSEL FUNCTION OF FIRST KIND
C***AUTHOR  AMOS, DONALD E., SANDIA NATIONAL LABORATORIES
C***PURPOSE  TO COMPUTE THE J-BESSEL FUNCTION OF A COMPLEX ARGUMENT
C***DESCRIPTION
C
C                      ***A DOUBLE PRECISION ROUTINE***
C         ON KODE=1, ZBESJ COMPUTES AN N MEMBER  SEQUENCE OF COMPLEX
C         BESSEL FUNCTIONS CY(I)=J(FNU+I-1,Z) FOR REAL, NONNEGATIVE
C         ORDERS FNU+I-1, I=1,...,N AND COMPLEX Z IN THE CUT PLANE
C         -PI.LT.ARG(Z).LE.PI. ON KODE=2, ZBESJ RETURNS THE SCALED
C         FUNCTIONS
C
C         CY(I)=EXP(-ABS(Y))*J(FNU+I-1,Z)   I = 1,...,N , Y=AIMAG(Z)
C
C         WHICH REMOVE THE EXPONENTIAL GROWTH IN BOTH THE UPPER AND
C         LOWER HALF PLANES FOR Z TO INFINITY. DEFINITIONS AND NOTATION
C         ARE FOUND IN THE NBS HANDBOOK OF MATHEMATICAL FUNCTIONS
C         (REF. 1).
C
C         INPUT      ZR,ZI,FNU ARE DOUBLE PRECISION
C           ZR,ZI  - Z=CMPLX(ZR,ZI),  -PI.LT.ARG(Z).LE.PI
C           FNU    - ORDER OF INITIAL J FUNCTION, FNU.GE.0.0D0
C           KODE   - A PARAMETER TO INDICATE THE SCALING OPTION
C                    KODE= 1  RETURNS
C                             CY(I)=J(FNU+I-1,Z), I=1,...,N
C                        = 2  RETURNS
C                             CY(I)=J(FNU+I-1,Z)EXP(-ABS(Y)), I=1,...,N
C           N      - NUMBER OF MEMBERS OF THE SEQUENCE, N.GE.1
C
C         OUTPUT     CYR,CYI ARE DOUBLE PRECISION
C           CYR,CYI- DOUBLE PRECISION VECTORS WHOSE FIRST N COMPONENTS
C                    CONTAIN REAL AND IMAGINARY PARTS FOR THE SEQUENCE
C                    CY(I)=J(FNU+I-1,Z)  OR
C                    CY(I)=J(FNU+I-1,Z)EXP(-ABS(Y))  I=1,...,N
C                    DEPENDING ON KODE, Y=AIMAG(Z).
C           NZ     - NUMBER OF COMPONENTS SET TO ZERO DUE TO UNDERFLOW,
C                    NZ= 0   , NORMAL RETURN
C                    NZ.GT.0 , LAST NZ COMPONENTS OF CY SET  ZERO DUE
C                              TO UNDERFLOW, CY(I)=CMPLX(0.0D0,0.0D0),
C                              I = N-NZ+1,...,N
C           IERR   - ERROR FLAG
C                    IERR=0, NORMAL RETURN - COMPUTATION COMPLETED
C                    IERR=1, INPUT ERROR   - NO COMPUTATION
C                    IERR=2, OVERFLOW      - NO COMPUTATION, AIMAG(Z)
C                            TOO LARGE ON KODE=1
C                    IERR=3, CABS(Z) OR FNU+N-1 LARGE - COMPUTATION DONE
C                            BUT LOSSES OF SIGNIFCANCE BY ARGUMENT
C                            REDUCTION PRODUCE LESS THAN HALF OF MACHINE
C                            ACCURACY
C                    IERR=4, CABS(Z) OR FNU+N-1 TOO LARGE - NO COMPUTA-
C                            TION BECAUSE OF COMPLETE LOSSES OF SIGNIFI-
C                            CANCE BY ARGUMENT REDUCTION
C                    IERR=5, ERROR              - NO COMPUTATION,
C                            ALGORITHM TERMINATION CONDITION NOT MET
C
C***LONG DESCRIPTION
C
C         THE COMPUTATION IS CARRIED OUT BY THE FORMULA
C
C         J(FNU,Z)=EXP( FNU*PI*I/2)*I(FNU,-I*Z)    AIMAG(Z).GE.0.0
C
C         J(FNU,Z)=EXP(-FNU*PI*I/2)*I(FNU, I*Z)    AIMAG(Z).LT.0.0
C
C         WHERE I**2 = -1 AND I(FNU,Z) IS THE I BESSEL FUNCTION.
C
C         FOR NEGATIVE ORDERS,THE FORMULA
C
C              J(-FNU,Z) = J(FNU,Z)*COS(PI*FNU) - Y(FNU,Z)*SIN(PI*FNU)
C
C         CAN BE USED. HOWEVER,FOR LARGE ORDERS CLOSE TO INTEGERS, THE
C         THE FUNCTION CHANGES RADICALLY. WHEN FNU IS A LARGE POSITIVE
C         INTEGER,THE MAGNITUDE OF J(-FNU,Z)=J(FNU,Z)*COS(PI*FNU) IS A
C         LARGE NEGATIVE POWER OF TEN. BUT WHEN FNU IS NOT AN INTEGER,
C         Y(FNU,Z) DOMINATES IN MAGNITUDE WITH A LARGE POSITIVE POWER OF
C         TEN AND THE MOST THAT THE SECOND TERM CAN BE REDUCED IS BY
C         UNIT ROUNDOFF FROM THE COEFFICIENT. THUS, WIDE CHANGES CAN
C         OCCUR WITHIN UNIT ROUNDOFF OF A LARGE INTEGER FOR FNU. HERE,
C         LARGE MEANS FNU.GT.CABS(Z).
C
C         IN MOST COMPLEX VARIABLE COMPUTATION, ONE MUST EVALUATE ELE-
C         MENTARY FUNCTIONS. WHEN THE MAGNITUDE OF Z OR FNU+N-1 IS
C         LARGE, LOSSES OF SIGNIFICANCE BY ARGUMENT REDUCTION OCCUR.
C         CONSEQUENTLY, IF EITHER ONE EXCEEDS U1=SQRT(0.5/UR), THEN
C         LOSSES EXCEEDING HALF PRECISION ARE LIKELY AND AN ERROR FLAG
C         IERR=3 IS TRIGGERED WHERE UR=DMAX1(D1MACH(4),1.0D-18) IS
C         DOUBLE PRECISION UNIT ROUNDOFF LIMITED TO 18 DIGITS PRECISION.
C         IF EITHER IS LARGER THAN U2=0.5/UR, THEN ALL SIGNIFICANCE IS
C         LOST AND IERR=4. IN ORDER TO USE THE INT FUNCTION, ARGUMENTS
C         MUST BE FURTHER RESTRICTED NOT TO EXCEED THE LARGEST MACHINE
C         INTEGER, U3=I1MACH(9). THUS, THE MAGNITUDE OF Z AND FNU+N-1 IS
C         RESTRICTED BY MIN(U2,U3). ON 32 BIT MACHINES, U1,U2, AND U3
C         ARE APPROXIMATELY 2.0E+3, 4.2E+6, 2.1E+9 IN SINGLE PRECISION
C         ARITHMETIC AND 1.3E+8, 1.8E+16, 2.1E+9 IN DOUBLE PRECISION
C         ARITHMETIC RESPECTIVELY. THIS MAKES U2 AND U3 LIMITING IN
C         THEIR RESPECTIVE ARITHMETICS. THIS MEANS THAT ONE CAN EXPECT
C         TO RETAIN, IN THE WORST CASES ON 32 BIT MACHINES, NO DIGITS
C         IN SINGLE AND ONLY 7 DIGITS IN DOUBLE PRECISION ARITHMETIC.
C         SIMILAR CONSIDERATIONS HOLD FOR OTHER MACHINES.
C
C         THE APPROXIMATE RELATIVE ERROR IN THE MAGNITUDE OF A COMPLEX
C         BESSEL FUNCTION CAN BE EXPRESSED BY P*10**S WHERE P=MAX(UNIT
C         ROUNDOFF,1.0E-18) IS THE NOMINAL PRECISION AND 10**S REPRE-
C         SENTS THE INCREASE IN ERROR DUE TO ARGUMENT REDUCTION IN THE
C         ELEMENTARY FUNCTIONS. HERE, S=MAX(1,ABS(LOG10(CABS(Z))),
C         ABS(LOG10(FNU))) APPROXIMATELY (I.E. S=MAX(1,ABS(EXPONENT OF
C         CABS(Z),ABS(EXPONENT OF FNU)) ). HOWEVER, THE PHASE ANGLE MAY
C         HAVE ONLY ABSOLUTE ACCURACY. THIS IS MOST LIKELY TO OCCUR WHEN
C         ONE COMPONENT (IN ABSOLUTE VALUE) IS LARGER THAN THE OTHER BY
C         SEVERAL ORDERS OF MAGNITUDE. IF ONE COMPONENT IS 10**K LARGER
C         THAN THE OTHER, THEN ONE CAN EXPECT ONLY MAX(ABS(LOG10(P))-K,
C         0) SIGNIFICANT DIGITS; OR, STATED ANOTHER WAY, WHEN K EXCEEDS
C         THE EXPONENT OF P, NO SIGNIFICANT DIGITS REMAIN IN THE SMALLER
C         COMPONENT. HOWEVER, THE PHASE ANGLE RETAINS ABSOLUTE ACCURACY
C         BECAUSE, IN COMPLEX ARITHMETIC WITH PRECISION P, THE SMALLER
C         COMPONENT WILL NOT (AS A RULE) DECREASE BELOW P TIMES THE
C         MAGNITUDE OF THE LARGER COMPONENT. IN THESE EXTREME CASES,
C         THE PRINCIPAL PHASE ANGLE IS ON THE ORDER OF +P, -P, PI/2-P,
C         OR -PI/2+P.
C
C***REFERENCES  HANDBOOK OF MATHEMATICAL FUNCTIONS BY M. ABRAMOWITZ
C                 AND I. A. STEGUN, NBS AMS SERIES 55, U.S. DEPT. OF
C                 COMMERCE, 1955.
C
C               COMPUTATION OF BESSEL FUNCTIONS OF COMPLEX ARGUMENT
C                 BY D. E. AMOS, SAND83-0083, MAY, 1983.
C
C               COMPUTATION OF BESSEL FUNCTIONS OF COMPLEX ARGUMENT
C                 AND LARGE ORDER BY D. E. AMOS, SAND83-0643, MAY, 1983
C
C               A SUBROUTINE PACKAGE FOR BESSEL FUNCTIONS OF A COMPLEX
C                 ARGUMENT AND NONNEGATIVE ORDER BY D. E. AMOS, SAND85-
C                 1018, MAY, 1985
C
C               A PORTABLE PACKAGE FOR BESSEL FUNCTIONS OF A COMPLEX
C                 ARGUMENT AND NONNEGATIVE ORDER BY D. E. AMOS, ACM
C                 TRANS. MATH. SOFTWARE, VOL. 12, NO. 3, SEPTEMBER 1986,
C                 PP 265-273.
C
C***ROUTINES CALLED  ZBINU,ZABS,I1MACH,D1MACH
C***END PROLOGUE  ZBESJ
C
C     COMPLEX CI,CSGN,CY,Z,ZN
      EXTERNAL ZABS
      DOUBLE PRECISION AA, ALIM, ARG, CII, CSGNI, CSGNR, CYI, CYR, DIG,
     * ELIM, FNU, FNUL, HPI, RL, R1M5, STR, TOL, ZI, ZNI, ZNR, ZR,
     * D1MACH, BB, FN, AZ, ZABS, ASCLE, RTOL, ATOL, STI
      INTEGER I, IERR, INU, INUH, IR, K, KODE, K1, K2, N, NL, NZ, I1MACH
      DIMENSION CYR(N), CYI(N)
      DATA HPI /1.57079632679489662D0/
C
C***FIRST EXECUTABLE STATEMENT  ZBESJ
      IERR = 0
      NZ=0
      IF (FNU.LT.0.0D0) IERR=1
      IF (KODE.LT.1 .OR. KODE.GT.2) IERR=1
      IF (N.LT.1) IERR=1
      IF (IERR.NE.0) RETURN
C-----------------------------------------------------------------------
C     SET PARAMETERS RELATED TO MACHINE CONSTANTS.
C     TOL IS THE APPROXIMATE UNIT ROUNDOFF LIMITED TO 1.0E-18.
C     ELIM IS THE APPROXIMATE EXPONENTIAL OVER- AND UNDERFLOW LIMIT.
C     EXP(-ELIM).LT.EXP(-ALIM)=EXP(-ELIM)/TOL    AND
C     EXP(ELIM).GT.EXP(ALIM)=EXP(ELIM)*TOL       ARE INTERVALS NEAR
C     UNDERFLOW AND OVERFLOW LIMITS WHERE SCALED ARITHMETIC IS DONE.
C     RL IS THE LOWER BOUNDARY OF THE ASYMPTOTIC EXPANSION FOR LARGE Z.
C     DIG = NUMBER OF BASE 10 DIGITS IN TOL = 10**(-DIG).
C     FNUL IS THE LOWER BOUNDARY OF THE ASYMPTOTIC SERIES FOR LARGE FNU.
C-----------------------------------------------------------------------
      TOL = DMAX1(D1MACH(4),1.0D-18)
      K1 = I1MACH(15)
      K2 = I1MACH(16)
      R1M5 = D1MACH(5)
      K = MIN0(IABS(K1),IABS(K2))
      ELIM = 2.303D0*(DBLE(FLOAT(K))*R1M5-3.0D0)
      K1 = I1MACH(14) - 1
      AA = R1M5*DBLE(FLOAT(K1))
      DIG = DMIN1(AA,18.0D0)
      AA = AA*2.303D0
      ALIM = ELIM + DMAX1(-AA,-41.45D0)
      RL = 1.2D0*DIG + 3.0D0
      FNUL = 10.0D0 + 6.0D0*(DIG-3.0D0)
C-----------------------------------------------------------------------
C     TEST FOR PROPER RANGE
C-----------------------------------------------------------------------
      AZ = ZABS(ZR,ZI)
      FN = FNU+DBLE(FLOAT(N-1))
      AA = 0.5D0/TOL
      BB=DBLE(FLOAT(I1MACH(9)))*0.5D0
      AA = DMIN1(AA,BB)
      IF (AZ.GT.AA) GO TO 260
      IF (FN.GT.AA) GO TO 260
      AA = DSQRT(AA)
      IF (AZ.GT.AA) IERR=3
      IF (FN.GT.AA) IERR=3
C-----------------------------------------------------------------------
C     CALCULATE CSGN=EXP(FNU*HPI*I) TO MINIMIZE LOSSES OF SIGNIFICANCE
C     WHEN FNU IS LARGE
C-----------------------------------------------------------------------
      CII = 1.0D0
      INU = INT(SNGL(FNU))
      INUH = INU/2
      IR = INU - 2*INUH
      ARG = (FNU-DBLE(FLOAT(INU-IR)))*HPI
      CSGNR = DCOS(ARG)
      CSGNI = DSIN(ARG)
      IF (MOD(INUH,2).EQ.0) GO TO 40
      CSGNR = -CSGNR
      CSGNI = -CSGNI
   40 CONTINUE
C-----------------------------------------------------------------------
C     ZN IS IN THE RIGHT HALF PLANE
C-----------------------------------------------------------------------
      ZNR = ZI
      ZNI = -ZR
      IF (ZI.GE.0.0D0) GO TO 50
      ZNR = -ZNR
      ZNI = -ZNI
      CSGNI = -CSGNI
      CII = -CII
   50 CONTINUE
      CALL ZBINU(ZNR, ZNI, FNU, KODE, N, CYR, CYI, NZ, RL, FNUL, TOL,
     * ELIM, ALIM)
      IF (NZ.LT.0) GO TO 130
      NL = N - NZ
      IF (NL.EQ.0) RETURN
      RTOL = 1.0D0/TOL
      ASCLE = D1MACH(1)*RTOL*1.0D+3
      DO 60 I=1,NL
C       STR = CYR(I)*CSGNR - CYI(I)*CSGNI
C       CYI(I) = CYR(I)*CSGNI + CYI(I)*CSGNR
C       CYR(I) = STR
        AA = CYR(I)
        BB = CYI(I)
        ATOL = 1.0D0
        IF (DMAX1(DABS(AA),DABS(BB)).GT.ASCLE) GO TO 55
          AA = AA*RTOL
          BB = BB*RTOL
          ATOL = TOL
   55   CONTINUE
        STR = AA*CSGNR - BB*CSGNI
        STI = AA*CSGNI + BB*CSGNR
        CYR(I) = STR*ATOL
        CYI(I) = STI*ATOL
        STR = -CSGNI*CII
        CSGNI = CSGNR*CII
        CSGNR = STR
   60 CONTINUE
      RETURN
  130 CONTINUE
      IF(NZ.EQ.(-2)) GO TO 140
      NZ = 0
      IERR = 2
      RETURN
  140 CONTINUE
      NZ=0
      IERR=5
      RETURN
  260 CONTINUE
      NZ=0
      IERR=4
      RETURN
      END
      SUBROUTINE ZBESK(ZR, ZI, FNU, KODE, N, CYR, CYI, NZ, IERR)
C***BEGIN PROLOGUE  ZBESK
C***DATE WRITTEN   830501   (YYMMDD)
C***REVISION DATE  890801, 930101   (YYMMDD)
C***CATEGORY NO.  B5K
C***KEYWORDS  K-BESSEL FUNCTION,COMPLEX BESSEL FUNCTION,
C             MODIFIED BESSEL FUNCTION OF THE SECOND KIND,
C             BESSEL FUNCTION OF THE THIRD KIND
C***AUTHOR  AMOS, DONALD E., SANDIA NATIONAL LABORATORIES
C***PURPOSE  TO COMPUTE K-BESSEL FUNCTIONS OF COMPLEX ARGUMENT
C***DESCRIPTION
C
C                      ***A DOUBLE PRECISION ROUTINE***
C
C         ON KODE=1, ZBESK COMPUTES AN N MEMBER SEQUENCE OF COMPLEX
C         BESSEL FUNCTIONS CY(J)=K(FNU+J-1,Z) FOR REAL, NONNEGATIVE
C         ORDERS FNU+J-1, J=1,...,N AND COMPLEX Z.NE.CMPLX(0.0,0.0)
C         IN THE CUT PLANE -PI.LT.ARG(Z).LE.PI. ON KODE=2, ZBESK
C         RETURNS THE SCALED K FUNCTIONS,
C
C         CY(J)=EXP(Z)*K(FNU+J-1,Z) , J=1,...,N,
C
C         WHICH REMOVE THE EXPONENTIAL BEHAVIOR IN BOTH THE LEFT AND
C         RIGHT HALF PLANES FOR Z TO INFINITY. DEFINITIONS AND
C         NOTATION ARE FOUND IN THE NBS HANDBOOK OF MATHEMATICAL
C         FUNCTIONS (REF. 1).
C
C         INPUT      ZR,ZI,FNU ARE DOUBLE PRECISION
C           ZR,ZI  - Z=CMPLX(ZR,ZI), Z.NE.CMPLX(0.0D0,0.0D0),
C                    -PI.LT.ARG(Z).LE.PI
C           FNU    - ORDER OF INITIAL K FUNCTION, FNU.GE.0.0D0
C           N      - NUMBER OF MEMBERS OF THE SEQUENCE, N.GE.1
C           KODE   - A PARAMETER TO INDICATE THE SCALING OPTION
C                    KODE= 1  RETURNS
C                             CY(I)=K(FNU+I-1,Z), I=1,...,N
C                        = 2  RETURNS
C                             CY(I)=K(FNU+I-1,Z)*EXP(Z), I=1,...,N
C
C         OUTPUT     CYR,CYI ARE DOUBLE PRECISION
C           CYR,CYI- DOUBLE PRECISION VECTORS WHOSE FIRST N COMPONENTS
C                    CONTAIN REAL AND IMAGINARY PARTS FOR THE SEQUENCE
C                    CY(I)=K(FNU+I-1,Z), I=1,...,N OR
C                    CY(I)=K(FNU+I-1,Z)*EXP(Z), I=1,...,N
C                    DEPENDING ON KODE
C           NZ     - NUMBER OF COMPONENTS SET TO ZERO DUE TO UNDERFLOW.
C                    NZ= 0   , NORMAL RETURN
C                    NZ.GT.0 , FIRST NZ COMPONENTS OF CY SET TO ZERO DUE
C                              TO UNDERFLOW, CY(I)=CMPLX(0.0D0,0.0D0),
C                              I=1,...,N WHEN X.GE.0.0. WHEN X.LT.0.0
C                              NZ STATES ONLY THE NUMBER OF UNDERFLOWS
C                              IN THE SEQUENCE.
C
C           IERR   - ERROR FLAG
C                    IERR=0, NORMAL RETURN - COMPUTATION COMPLETED
C                    IERR=1, INPUT ERROR   - NO COMPUTATION
C                    IERR=2, OVERFLOW      - NO COMPUTATION, FNU IS
C                            TOO LARGE OR CABS(Z) IS TOO SMALL OR BOTH
C                    IERR=3, CABS(Z) OR FNU+N-1 LARGE - COMPUTATION DONE
C                            BUT LOSSES OF SIGNIFCANCE BY ARGUMENT
C                            REDUCTION PRODUCE LESS THAN HALF OF MACHINE
C                            ACCURACY
C                    IERR=4, CABS(Z) OR FNU+N-1 TOO LARGE - NO COMPUTA-
C                            TION BECAUSE OF COMPLETE LOSSES OF SIGNIFI-
C                            CANCE BY ARGUMENT REDUCTION
C                    IERR=5, ERROR              - NO COMPUTATION,
C                            ALGORITHM TERMINATION CONDITION NOT MET
C
C***LONG DESCRIPTION
C
C         EQUATIONS OF THE REFERENCE ARE IMPLEMENTED FOR SMALL ORDERS
C         DNU AND DNU+1.0 IN THE RIGHT HALF PLANE X.GE.0.0. FORWARD
C         RECURRENCE GENERATES HIGHER ORDERS. K IS CONTINUED TO THE LEFT
C         HALF PLANE BY THE RELATION
C
C         K(FNU,Z*EXP(MP)) = EXP(-MP*FNU)*K(FNU,Z)-MP*I(FNU,Z)
C         MP=MR*PI*I, MR=+1 OR -1, RE(Z).GT.0, I**2=-1
C
C         WHERE I(FNU,Z) IS THE I BESSEL FUNCTION.
C
C         FOR LARGE ORDERS, FNU.GT.FNUL, THE K FUNCTION IS COMPUTED
C         BY MEANS OF ITS UNIFORM ASYMPTOTIC EXPANSIONS.
C
C         FOR NEGATIVE ORDERS, THE FORMULA
C
C                       K(-FNU,Z) = K(FNU,Z)
C
C         CAN BE USED.
C
C         ZBESK ASSUMES THAT A SIGNIFICANT DIGIT SINH(X) FUNCTION IS
C         AVAILABLE.
C
C         IN MOST COMPLEX VARIABLE COMPUTATION, ONE MUST EVALUATE ELE-
C         MENTARY FUNCTIONS. WHEN THE MAGNITUDE OF Z OR FNU+N-1 IS
C         LARGE, LOSSES OF SIGNIFICANCE BY ARGUMENT REDUCTION OCCUR.
C         CONSEQUENTLY, IF EITHER ONE EXCEEDS U1=SQRT(0.5/UR), THEN
C         LOSSES EXCEEDING HALF PRECISION ARE LIKELY AND AN ERROR FLAG
C         IERR=3 IS TRIGGERED WHERE UR=DMAX1(D1MACH(4),1.0D-18) IS
C         DOUBLE PRECISION UNIT ROUNDOFF LIMITED TO 18 DIGITS PRECISION.
C         IF EITHER IS LARGER THAN U2=0.5/UR, THEN ALL SIGNIFICANCE IS
C         LOST AND IERR=4. IN ORDER TO USE THE INT FUNCTION, ARGUMENTS
C         MUST BE FURTHER RESTRICTED NOT TO EXCEED THE LARGEST MACHINE
C         INTEGER, U3=I1MACH(9). THUS, THE MAGNITUDE OF Z AND FNU+N-1 IS
C         RESTRICTED BY MIN(U2,U3). ON 32 BIT MACHINES, U1,U2, AND U3
C         ARE APPROXIMATELY 2.0E+3, 4.2E+6, 2.1E+9 IN SINGLE PRECISION
C         ARITHMETIC AND 1.3E+8, 1.8E+16, 2.1E+9 IN DOUBLE PRECISION
C         ARITHMETIC RESPECTIVELY. THIS MAKES U2 AND U3 LIMITING IN
C         THEIR RESPECTIVE ARITHMETICS. THIS MEANS THAT ONE CAN EXPECT
C         TO RETAIN, IN THE WORST CASES ON 32 BIT MACHINES, NO DIGITS
C         IN SINGLE AND ONLY 7 DIGITS IN DOUBLE PRECISION ARITHMETIC.
C         SIMILAR CONSIDERATIONS HOLD FOR OTHER MACHINES.
C
C         THE APPROXIMATE RELATIVE ERROR IN THE MAGNITUDE OF A COMPLEX
C         BESSEL FUNCTION CAN BE EXPRESSED BY P*10**S WHERE P=MAX(UNIT
C         ROUNDOFF,1.0E-18) IS THE NOMINAL PRECISION AND 10**S REPRE-
C         SENTS THE INCREASE IN ERROR DUE TO ARGUMENT REDUCTION IN THE
C         ELEMENTARY FUNCTIONS. HERE, S=MAX(1,ABS(LOG10(CABS(Z))),
C         ABS(LOG10(FNU))) APPROXIMATELY (I.E. S=MAX(1,ABS(EXPONENT OF
C         CABS(Z),ABS(EXPONENT OF FNU)) ). HOWEVER, THE PHASE ANGLE MAY
C         HAVE ONLY ABSOLUTE ACCURACY. THIS IS MOST LIKELY TO OCCUR WHEN
C         ONE COMPONENT (IN ABSOLUTE VALUE) IS LARGER THAN THE OTHER BY
C         SEVERAL ORDERS OF MAGNITUDE. IF ONE COMPONENT IS 10**K LARGER
C         THAN THE OTHER, THEN ONE CAN EXPECT ONLY MAX(ABS(LOG10(P))-K,
C         0) SIGNIFICANT DIGITS; OR, STATED ANOTHER WAY, WHEN K EXCEEDS
C         THE EXPONENT OF P, NO SIGNIFICANT DIGITS REMAIN IN THE SMALLER
C         COMPONENT. HOWEVER, THE PHASE ANGLE RETAINS ABSOLUTE ACCURACY
C         BECAUSE, IN COMPLEX ARITHMETIC WITH PRECISION P, THE SMALLER
C         COMPONENT WILL NOT (AS A RULE) DECREASE BELOW P TIMES THE
C         MAGNITUDE OF THE LARGER COMPONENT. IN THESE EXTREME CASES,
C         THE PRINCIPAL PHASE ANGLE IS ON THE ORDER OF +P, -P, PI/2-P,
C         OR -PI/2+P.
C
C***REFERENCES  HANDBOOK OF MATHEMATICAL FUNCTIONS BY M. ABRAMOWITZ
C                 AND I. A. STEGUN, NBS AMS SERIES 55, U.S. DEPT. OF
C                 COMMERCE, 1955.
C
C               COMPUTATION OF BESSEL FUNCTIONS OF COMPLEX ARGUMENT
C                 BY D. E. AMOS, SAND83-0083, MAY, 1983.
C
C               COMPUTATION OF BESSEL FUNCTIONS OF COMPLEX ARGUMENT
C                 AND LARGE ORDER BY D. E. AMOS, SAND83-0643, MAY, 1983.
C
C               A SUBROUTINE PACKAGE FOR BESSEL FUNCTIONS OF A COMPLEX
C                 ARGUMENT AND NONNEGATIVE ORDER BY D. E. AMOS, SAND85-
C                 1018, MAY, 1985
C
C               A PORTABLE PACKAGE FOR BESSEL FUNCTIONS OF A COMPLEX
C                 ARGUMENT AND NONNEGATIVE ORDER BY D. E. AMOS, ACM
C                 TRANS. MATH. SOFTWARE, VOL. 12, NO. 3, SEPTEMBER 1986,
C                 PP 265-273.
C
C***ROUTINES CALLED  ZACON,ZBKNU,ZBUNK,ZUOIK,ZABS,I1MACH,D1MACH
C***END PROLOGUE  ZBESK
C
C     COMPLEX CY,Z
      EXTERNAL ZABS
      DOUBLE PRECISION AA, ALIM, ALN, ARG, AZ, CYI, CYR, DIG, ELIM, FN,
     * FNU, FNUL, RL, R1M5, TOL, UFL, ZI, ZR, D1MACH, ZABS, BB
      INTEGER IERR, K, KODE, K1, K2, MR, N, NN, NUF, NW, NZ, I1MACH
      DIMENSION CYR(N), CYI(N)
C***FIRST EXECUTABLE STATEMENT  ZBESK
      IERR = 0
      NZ=0
      IF (ZI.EQ.0.0E0 .AND. ZR.EQ.0.0E0) IERR=1
      IF (FNU.LT.0.0D0) IERR=1
      IF (KODE.LT.1 .OR. KODE.GT.2) IERR=1
      IF (N.LT.1) IERR=1
      IF (IERR.NE.0) RETURN
      NN = N
C-----------------------------------------------------------------------
C     SET PARAMETERS RELATED TO MACHINE CONSTANTS.
C     TOL IS THE APPROXIMATE UNIT ROUNDOFF LIMITED TO 1.0E-18.
C     ELIM IS THE APPROXIMATE EXPONENTIAL OVER- AND UNDERFLOW LIMIT.
C     EXP(-ELIM).LT.EXP(-ALIM)=EXP(-ELIM)/TOL    AND
C     EXP(ELIM).GT.EXP(ALIM)=EXP(ELIM)*TOL       ARE INTERVALS NEAR
C     UNDERFLOW AND OVERFLOW LIMITS WHERE SCALED ARITHMETIC IS DONE.
C     RL IS THE LOWER BOUNDARY OF THE ASYMPTOTIC EXPANSION FOR LARGE Z.
C     DIG = NUMBER OF BASE 10 DIGITS IN TOL = 10**(-DIG).
C     FNUL IS THE LOWER BOUNDARY OF THE ASYMPTOTIC SERIES FOR LARGE FNU
C-----------------------------------------------------------------------
      TOL = DMAX1(D1MACH(4),1.0D-18)
      K1 = I1MACH(15)
      K2 = I1MACH(16)
      R1M5 = D1MACH(5)
      K = MIN0(IABS(K1),IABS(K2))
      ELIM = 2.303D0*(DBLE(FLOAT(K))*R1M5-3.0D0)
      K1 = I1MACH(14) - 1
      AA = R1M5*DBLE(FLOAT(K1))
      DIG = DMIN1(AA,18.0D0)
      AA = AA*2.303D0
      ALIM = ELIM + DMAX1(-AA,-41.45D0)
      FNUL = 10.0D0 + 6.0D0*(DIG-3.0D0)
      RL = 1.2D0*DIG + 3.0D0
C-----------------------------------------------------------------------------
C     TEST FOR PROPER RANGE
C-----------------------------------------------------------------------
      AZ = ZABS(ZR,ZI)
      FN = FNU + DBLE(FLOAT(NN-1))
      AA = 0.5D0/TOL
      BB=DBLE(FLOAT(I1MACH(9)))*0.5D0
      AA = DMIN1(AA,BB)
      IF (AZ.GT.AA) GO TO 260
      IF (FN.GT.AA) GO TO 260
      AA = DSQRT(AA)
      IF (AZ.GT.AA) IERR=3
      IF (FN.GT.AA) IERR=3
C-----------------------------------------------------------------------
C     OVERFLOW TEST ON THE LAST MEMBER OF THE SEQUENCE
C-----------------------------------------------------------------------
C     UFL = DEXP(-ELIM)
      UFL = D1MACH(1)*1.0D+3
      IF (AZ.LT.UFL) GO TO 180
      IF (FNU.GT.FNUL) GO TO 80
      IF (FN.LE.1.0D0) GO TO 60
      IF (FN.GT.2.0D0) GO TO 50
      IF (AZ.GT.TOL) GO TO 60
      ARG = 0.5D0*AZ
      ALN = -FN*DLOG(ARG)
      IF (ALN.GT.ELIM) GO TO 180
      GO TO 60
   50 CONTINUE
      CALL ZUOIK(ZR, ZI, FNU, KODE, 2, NN, CYR, CYI, NUF, TOL, ELIM,
     * ALIM)
      IF (NUF.LT.0) GO TO 180
      NZ = NZ + NUF
      NN = NN - NUF
C-----------------------------------------------------------------------
C     HERE NN=N OR NN=0 SINCE NUF=0,NN, OR -1 ON RETURN FROM CUOIK
C     IF NUF=NN, THEN CY(I)=CZERO FOR ALL I
C-----------------------------------------------------------------------
      IF (NN.EQ.0) GO TO 100
   60 CONTINUE
      IF (ZR.LT.0.0D0) GO TO 70
C-----------------------------------------------------------------------
C     RIGHT HALF PLANE COMPUTATION, REAL(Z).GE.0.
C-----------------------------------------------------------------------
      CALL ZBKNU(ZR, ZI, FNU, KODE, NN, CYR, CYI, NW, TOL, ELIM, ALIM)
      IF (NW.LT.0) GO TO 200
      NZ=NW
      RETURN
C-----------------------------------------------------------------------
C     LEFT HALF PLANE COMPUTATION
C     PI/2.LT.ARG(Z).LE.PI AND -PI.LT.ARG(Z).LT.-PI/2.
C-----------------------------------------------------------------------
   70 CONTINUE
      IF (NZ.NE.0) GO TO 180
      MR = 1
      IF (ZI.LT.0.0D0) MR = -1
      CALL ZACON(ZR, ZI, FNU, KODE, MR, NN, CYR, CYI, NW, RL, FNUL,
     * TOL, ELIM, ALIM)
      IF (NW.LT.0) GO TO 200
      NZ=NW
      RETURN
C-----------------------------------------------------------------------
C     UNIFORM ASYMPTOTIC EXPANSIONS FOR FNU.GT.FNUL
C-----------------------------------------------------------------------
   80 CONTINUE
      MR = 0
      IF (ZR.GE.0.0D0) GO TO 90
      MR = 1
      IF (ZI.LT.0.0D0) MR = -1
   90 CONTINUE
      CALL ZBUNK(ZR, ZI, FNU, KODE, MR, NN, CYR, CYI, NW, TOL, ELIM,
     * ALIM)
      IF (NW.LT.0) GO TO 200
      NZ = NZ + NW
      RETURN
  100 CONTINUE
      IF (ZR.LT.0.0D0) GO TO 180
      RETURN
  180 CONTINUE
      NZ = 0
      IERR=2
      RETURN
  200 CONTINUE
      IF(NW.EQ.(-1)) GO TO 180
      NZ=0
      IERR=5
      RETURN
  260 CONTINUE
      NZ=0
      IERR=4
      RETURN
      END
      SUBROUTINE ZBESY(ZR, ZI, FNU, KODE, N, CYR, CYI, NZ, CWRKR,
     *           CWRKI, IERR)
C***BEGIN PROLOGUE  ZBESY
C***DATE WRITTEN   830501   (YYMMDD)
C***REVISION DATE  890801, 930101   (YYMMDD)
C***CATEGORY NO.  B5K
C***KEYWORDS  Y-BESSEL FUNCTION,BESSEL FUNCTION OF COMPLEX ARGUMENT,
C             BESSEL FUNCTION OF SECOND KIND
C***AUTHOR  AMOS, DONALD E., SANDIA NATIONAL LABORATORIES
C***PURPOSE  TO COMPUTE THE Y-BESSEL FUNCTION OF A COMPLEX ARGUMENT
C***DESCRIPTION
C
C                      ***A DOUBLE PRECISION ROUTINE***
C
C         ON KODE=1, ZBESY COMPUTES AN N MEMBER SEQUENCE OF COMPLEX
C         BESSEL FUNCTIONS CY(I)=Y(FNU+I-1,Z) FOR REAL, NONNEGATIVE
C         ORDERS FNU+I-1, I=1,...,N AND COMPLEX Z IN THE CUT PLANE
C         -PI.LT.ARG(Z).LE.PI. ON KODE=2, ZBESY RETURNS THE SCALED
C         FUNCTIONS
C
C         CY(I)=EXP(-ABS(Y))*Y(FNU+I-1,Z)   I = 1,...,N , Y=AIMAG(Z)
C
C         WHICH REMOVE THE EXPONENTIAL GROWTH IN BOTH THE UPPER AND
C         LOWER HALF PLANES FOR Z TO INFINITY. DEFINITIONS AND NOTATION
C         ARE FOUND IN THE NBS HANDBOOK OF MATHEMATICAL FUNCTIONS
C         (REF. 1).
C
C         INPUT      ZR,ZI,FNU ARE DOUBLE PRECISION
C           ZR,ZI  - Z=CMPLX(ZR,ZI), Z.NE.CMPLX(0.0D0,0.0D0),
C                    -PI.LT.ARG(Z).LE.PI
C           FNU    - ORDER OF INITIAL Y FUNCTION, FNU.GE.0.0D0
C           KODE   - A PARAMETER TO INDICATE THE SCALING OPTION
C                    KODE= 1  RETURNS
C                             CY(I)=Y(FNU+I-1,Z), I=1,...,N
C                        = 2  RETURNS
C                             CY(I)=Y(FNU+I-1,Z)*EXP(-ABS(Y)), I=1,...,N
C                             WHERE Y=AIMAG(Z)
C           N      - NUMBER OF MEMBERS OF THE SEQUENCE, N.GE.1
C           CWRKR, - DOUBLE PRECISION WORK VECTORS OF DIMENSION AT
C           CWRKI    AT LEAST N
C
C         OUTPUT     CYR,CYI ARE DOUBLE PRECISION
C           CYR,CYI- DOUBLE PRECISION VECTORS WHOSE FIRST N COMPONENTS
C                    CONTAIN REAL AND IMAGINARY PARTS FOR THE SEQUENCE
C                    CY(I)=Y(FNU+I-1,Z)  OR
C                    CY(I)=Y(FNU+I-1,Z)*EXP(-ABS(Y))  I=1,...,N
C                    DEPENDING ON KODE.
C           NZ     - NZ=0 , A NORMAL RETURN
C                    NZ.GT.0 , NZ COMPONENTS OF CY SET TO ZERO DUE TO
C                    UNDERFLOW (GENERALLY ON KODE=2)
C           IERR   - ERROR FLAG
C                    IERR=0, NORMAL RETURN - COMPUTATION COMPLETED
C                    IERR=1, INPUT ERROR   - NO COMPUTATION
C                    IERR=2, OVERFLOW      - NO COMPUTATION, FNU IS
C                            TOO LARGE OR CABS(Z) IS TOO SMALL OR BOTH
C                    IERR=3, CABS(Z) OR FNU+N-1 LARGE - COMPUTATION DONE
C                            BUT LOSSES OF SIGNIFCANCE BY ARGUMENT
C                            REDUCTION PRODUCE LESS THAN HALF OF MACHINE
C                            ACCURACY
C                    IERR=4, CABS(Z) OR FNU+N-1 TOO LARGE - NO COMPUTA-
C                            TION BECAUSE OF COMPLETE LOSSES OF SIGNIFI-
C                            CANCE BY ARGUMENT REDUCTION
C                    IERR=5, ERROR              - NO COMPUTATION,
C                            ALGORITHM TERMINATION CONDITION NOT MET
C
C***LONG DESCRIPTION
C
C         THE COMPUTATION IS CARRIED OUT IN TERMS OF THE I(FNU,Z) AND
C         K(FNU,Z) BESSEL FUNCTIONS IN THE RIGHT HALF PLANE BY
C
C             Y(FNU,Z) = I*CC*I(FNU,ARG) - (2/PI)*CONJG(CC)*K(FNU,ARG)
C
C             Y(FNU,Z) = CONJG(Y(FNU,CONJG(Z)))
C
C         FOR AIMAG(Z).GE.0 AND AIMAG(Z).LT.0 RESPECTIVELY, WHERE
C         CC=EXP(I*PI*FNU/2), ARG=Z*EXP(-I*PI/2) AND I**2=-1.
C
C         FOR NEGATIVE ORDERS,THE FORMULA
C
C              Y(-FNU,Z) = Y(FNU,Z)*COS(PI*FNU) + J(FNU,Z)*SIN(PI*FNU)
C
C         CAN BE USED. HOWEVER,FOR LARGE ORDERS CLOSE TO HALF ODD
C         INTEGERS THE FUNCTION CHANGES RADICALLY. WHEN FNU IS A LARGE
C         POSITIVE HALF ODD INTEGER,THE MAGNITUDE OF Y(-FNU,Z)=J(FNU,Z)*
C         SIN(PI*FNU) IS A LARGE NEGATIVE POWER OF TEN. BUT WHEN FNU IS
C         NOT A HALF ODD INTEGER, Y(FNU,Z) DOMINATES IN MAGNITUDE WITH A
C         LARGE POSITIVE POWER OF TEN AND THE MOST THAT THE SECOND TERM
C         CAN BE REDUCED IS BY UNIT ROUNDOFF FROM THE COEFFICIENT. THUS,
C         WIDE CHANGES CAN OCCUR WITHIN UNIT ROUNDOFF OF A LARGE HALF
C         ODD INTEGER. HERE, LARGE MEANS FNU.GT.CABS(Z).
C
C         IN MOST COMPLEX VARIABLE COMPUTATION, ONE MUST EVALUATE ELE-
C         MENTARY FUNCTIONS. WHEN THE MAGNITUDE OF Z OR FNU+N-1 IS
C         LARGE, LOSSES OF SIGNIFICANCE BY ARGUMENT REDUCTION OCCUR.
C         CONSEQUENTLY, IF EITHER ONE EXCEEDS U1=SQRT(0.5/UR), THEN
C         LOSSES EXCEEDING HALF PRECISION ARE LIKELY AND AN ERROR FLAG
C         IERR=3 IS TRIGGERED WHERE UR=DMAX1(D1MACH(4),1.0D-18) IS
C         DOUBLE PRECISION UNIT ROUNDOFF LIMITED TO 18 DIGITS PRECISION.
C         IF EITHER IS LARGER THAN U2=0.5/UR, THEN ALL SIGNIFICANCE IS
C         LOST AND IERR=4. IN ORDER TO USE THE INT FUNCTION, ARGUMENTS
C         MUST BE FURTHER RESTRICTED NOT TO EXCEED THE LARGEST MACHINE
C         INTEGER, U3=I1MACH(9). THUS, THE MAGNITUDE OF Z AND FNU+N-1 IS
C         RESTRICTED BY MIN(U2,U3). ON 32 BIT MACHINES, U1,U2, AND U3
C         ARE APPROXIMATELY 2.0E+3, 4.2E+6, 2.1E+9 IN SINGLE PRECISION
C         ARITHMETIC AND 1.3E+8, 1.8E+16, 2.1E+9 IN DOUBLE PRECISION
C         ARITHMETIC RESPECTIVELY. THIS MAKES U2 AND U3 LIMITING IN
C         THEIR RESPECTIVE ARITHMETICS. THIS MEANS THAT ONE CAN EXPECT
C         TO RETAIN, IN THE WORST CASES ON 32 BIT MACHINES, NO DIGITS
C         IN SINGLE AND ONLY 7 DIGITS IN DOUBLE PRECISION ARITHMETIC.
C         SIMILAR CONSIDERATIONS HOLD FOR OTHER MACHINES.
C
C         THE APPROXIMATE RELATIVE ERROR IN THE MAGNITUDE OF A COMPLEX
C         BESSEL FUNCTION CAN BE EXPRESSED BY P*10**S WHERE P=MAX(UNIT
C         ROUNDOFF,1.0E-18) IS THE NOMINAL PRECISION AND 10**S REPRE-
C         SENTS THE INCREASE IN ERROR DUE TO ARGUMENT REDUCTION IN THE
C         ELEMENTARY FUNCTIONS. HERE, S=MAX(1,ABS(LOG10(CABS(Z))),
C         ABS(LOG10(FNU))) APPROXIMATELY (I.E. S=MAX(1,ABS(EXPONENT OF
C         CABS(Z),ABS(EXPONENT OF FNU)) ). HOWEVER, THE PHASE ANGLE MAY
C         HAVE ONLY ABSOLUTE ACCURACY. THIS IS MOST LIKELY TO OCCUR WHEN
C         ONE COMPONENT (IN ABSOLUTE VALUE) IS LARGER THAN THE OTHER BY
C         SEVERAL ORDERS OF MAGNITUDE. IF ONE COMPONENT IS 10**K LARGER
C         THAN THE OTHER, THEN ONE CAN EXPECT ONLY MAX(ABS(LOG10(P))-K,
C         0) SIGNIFICANT DIGITS; OR, STATED ANOTHER WAY, WHEN K EXCEEDS
C         THE EXPONENT OF P, NO SIGNIFICANT DIGITS REMAIN IN THE SMALLER
C         COMPONENT. HOWEVER, THE PHASE ANGLE RETAINS ABSOLUTE ACCURACY
C         BECAUSE, IN COMPLEX ARITHMETIC WITH PRECISION P, THE SMALLER
C         COMPONENT WILL NOT (AS A RULE) DECREASE BELOW P TIMES THE
C         MAGNITUDE OF THE LARGER COMPONENT. IN THESE EXTREME CASES,
C         THE PRINCIPAL PHASE ANGLE IS ON THE ORDER OF +P, -P, PI/2-P,
C         OR -PI/2+P.
C
C***REFERENCES  HANDBOOK OF MATHEMATICAL FUNCTIONS BY M. ABRAMOWITZ
C                 AND I. A. STEGUN, NBS AMS SERIES 55, U.S. DEPT. OF
C                 COMMERCE, 1955.
C
C               COMPUTATION OF BESSEL FUNCTIONS OF COMPLEX ARGUMENT
C                 BY D. E. AMOS, SAND83-0083, MAY, 1983.
C
C               COMPUTATION OF BESSEL FUNCTIONS OF COMPLEX ARGUMENT
C                 AND LARGE ORDER BY D. E. AMOS, SAND83-0643, MAY, 1983
C
C               A SUBROUTINE PACKAGE FOR BESSEL FUNCTIONS OF A COMPLEX
C                 ARGUMENT AND NONNEGATIVE ORDER BY D. E. AMOS, SAND85-
C                 1018, MAY, 1985
C
C               A PORTABLE PACKAGE FOR BESSEL FUNCTIONS OF A COMPLEX
C                 ARGUMENT AND NONNEGATIVE ORDER BY D. E. AMOS, ACM
C                 TRANS. MATH. SOFTWARE, VOL. 12, NO. 3, SEPTEMBER 1986,
C                 PP 265-273.
C
C***ROUTINES CALLED  ZBESI,ZBESK,I1MACH,D1MACH
C***END PROLOGUE  ZBESY
C
C     COMPLEX CWRK,CY,C1,C2,EX,HCI,Z,ZU,ZV
      DOUBLE PRECISION ARG, ASCLE, CIPI, CIPR, CSGNI, CSGNR, CSPNI,
     * CSPNR, CWRKI, CWRKR, CYI, CYR, D1M5, D1MACH, ELIM, EXI, EXR, EY,
     * FNU, FFNU, HPI, RHPI, STR, STI, TAY, TOL, ATOL, RTOL, ZI, ZR,
     * ZNI, ZNR, ZUI, ZUR, ZVI, ZVR, ZZI, ZZR
      INTEGER I, IERR, IFNU, I4, K, KODE, K1, K2, N, NZ, NZ1, NZ2,
     * I1MACH
      DIMENSION CYR(N), CYI(N), CWRKR(N), CWRKI(N), CIPR(4), CIPI(4)
      DATA CIPR(1),CIPR(2),CIPR(3),CIPR(4)/1.0D0, 0.0D0, -1.0D0, 0.0D0/
      DATA CIPI(1),CIPI(2),CIPI(3),CIPI(4)/0.0D0, 1.0D0, 0.0D0, -1.0D0/
      DATA HPI / 1.57079632679489662D0 /
C***FIRST EXECUTABLE STATEMENT  ZBESY
      IERR = 0
      NZ=0
      IF (ZR.EQ.0.0D0 .AND. ZI.EQ.0.0D0) IERR=1
      IF (FNU.LT.0.0D0) IERR=1
      IF (KODE.LT.1 .OR. KODE.GT.2) IERR=1
      IF (N.LT.1) IERR=1
      IF (IERR.NE.0) RETURN
      ZZR = ZR
      ZZI = ZI
      IF (ZI.LT.0.0D0) ZZI = -ZZI
      ZNR = ZZI
      ZNI = -ZZR
      CALL ZBESI(ZNR, ZNI, FNU, KODE, N, CYR, CYI, NZ1, IERR)
      IF (IERR.NE.0.AND.IERR.NE.3) GO TO 90
      CALL ZBESK(ZNR, ZNI, FNU, KODE, N, CWRKR, CWRKI, NZ2, IERR)
      IF (IERR.NE.0.AND.IERR.NE.3) GO TO 90
      NZ = MIN(NZ1,NZ2)
      IFNU = INT(SNGL(FNU))
      FFNU = FNU - DBLE(FLOAT(IFNU))
      ARG = HPI*FFNU
      CSGNR = COS(ARG)
      CSGNI = SIN(ARG)
      I4 = MOD(IFNU,4) + 1
      STR = CSGNR*CIPR(I4) - CSGNI*CIPI(I4)
      CSGNI = CSGNR*CIPI(I4) + CSGNI*CIPR(I4)
      CSGNR = STR
      RHPI = 1.0D0/HPI
      CSPNR = CSGNR*RHPI
      CSPNI = -CSGNI*RHPI
      STR = -CSGNI
      CSGNI = CSGNR
      CSGNR = STR
      IF (KODE.EQ.2) GO TO 60
      DO 50 I=1,N
C       CY(I) = CSGN*CY(I)-CSPN*CWRK(I)
        STR = CSGNR*CYR(I) - CSGNI*CYI(I)
        STR = STR - (CSPNR*CWRKR(I) - CSPNI*CWRKI(I))
        STI = CSGNR*CYI(I) + CSGNI*CYR(I)
        STI = STI - (CSPNR*CWRKI(I) + CSPNI*CWRKR(I))
        CYR(I) = STR
        CYI(I) = STI
        STR = - CSGNI
        CSGNI = CSGNR
        CSGNR = STR
        STR = CSPNI
        CSPNI = -CSPNR
        CSPNR = STR
   50 CONTINUE
      IF (ZI.LT.0.0D0) THEN
        DO 55 I=1,N
          CYI(I) = -CYI(I)
   55   CONTINUE
      ENDIF
      RETURN
   60 CONTINUE
      EXR = COS(ZR)
      EXI = SIN(ZR)
      TOL = MAX(D1MACH(4),1.0D-18)
      K1 = I1MACH(15)
      K2 = I1MACH(16)
      K = MIN(IABS(K1),IABS(K2))
      D1M5 = D1MACH(5)
C-----------------------------------------------------------------------
C     ELIM IS THE APPROXIMATE EXPONENTIAL UNDER- AND OVERFLOW LIMIT
C-----------------------------------------------------------------------
      ELIM = 2.303D0*(DBLE(FLOAT(K))*D1M5-3.0D0)
      EY = 0.0D0
      TAY = ABS(ZI+ZI)
      IF (TAY.LT.ELIM) EY = EXP(-TAY)
      STR = (EXR*CSPNR - EXI*CSPNI)*EY
      CSPNI = (EXR*CSPNI + EXI*CSPNR)*EY
      CSPNR = STR
      NZ = 0
      RTOL = 1.0D0/TOL
      ASCLE = D1MACH(1)*RTOL*1.0D+3
      DO 80 I=1,N
C----------------------------------------------------------------------
C       CY(I) = CSGN*CY(I)-CSPN*CWRK(I): PRODUCTS ARE COMPUTED IN
C       SCALED MODE IF CY(I) OR CWRK(I) ARE CLOSE TO UNDERFLOW TO 
C       PREVENT UNDERFLOW IN AN INTERMEDIATE COMPUTATION.
C----------------------------------------------------------------------
        ZVR = CWRKR(I)
        ZVI = CWRKI(I)
        ATOL=1.0D0
        IF (MAX(ABS(ZVR),ABS(ZVI)).GT.ASCLE) GO TO 75
          ZVR = ZVR*RTOL
          ZVI = ZVI*RTOL
          ATOL = TOL
   75   CONTINUE
        STR = (ZVR*CSPNR - ZVI*CSPNI)*ATOL
        ZVI = (ZVR*CSPNI + ZVI*CSPNR)*ATOL
        ZVR = STR
        ZUR = CYR(I)
        ZUI = CYI(I)
        ATOL=1.0D0
        IF (MAX(ABS(ZUR),ABS(ZUI)).GT.ASCLE) GO TO 85
          ZUR = ZUR*RTOL
          ZUI = ZUI*RTOL
          ATOL = TOL
   85   CONTINUE
        STR = (ZUR*CSGNR - ZUI*CSGNI)*ATOL
        ZUI = (ZUR*CSGNI + ZUI*CSGNR)*ATOL
        ZUR = STR
        CYR(I) = ZUR - ZVR
        CYI(I) = ZUI - ZVI
        IF (ZI.LT.0.0D0) CYI(I) = -CYI(I)
        IF (CYR(I).EQ.0.0D0 .AND. CYI(I).EQ.0.0D0 .AND. EY.EQ.0.0D0)
     &      NZ = NZ + 1
        STR = -CSGNI
        CSGNI = CSGNR
        CSGNR = STR
        STR = CSPNI
        CSPNI = -CSPNR
        CSPNR = STR
   80 CONTINUE
      RETURN
   90 CONTINUE
      NZ = 0
      RETURN
      END
      SUBROUTINE ZAIRY(ZR, ZI, ID, KODE, AIR, AII, NZ, IERR)
C***BEGIN PROLOGUE  ZAIRY
C***DATE WRITTEN   830501   (YYMMDD)
C***REVISION DATE  890801, 930101   (YYMMDD)
C***CATEGORY NO.  B5K
C***KEYWORDS  AIRY FUNCTION,BESSEL FUNCTIONS OF ORDER ONE THIRD
C***AUTHOR  AMOS, DONALD E., SANDIA NATIONAL LABORATORIES
C***PURPOSE  TO COMPUTE AIRY FUNCTIONS AI(Z) AND DAI(Z) FOR COMPLEX Z
C***DESCRIPTION
C
C                      ***A DOUBLE PRECISION ROUTINE***
C         ON KODE=1, ZAIRY COMPUTES THE COMPLEX AIRY FUNCTION AI(Z) OR
C         ITS DERIVATIVE DAI(Z)/DZ ON ID=0 OR ID=1 RESPECTIVELY. ON
C         KODE=2, A SCALING OPTION CEXP(ZTA)*AI(Z) OR CEXP(ZTA)*
C         DAI(Z)/DZ IS PROVIDED TO REMOVE THE EXPONENTIAL DECAY IN
C         -PI/3.LT.ARG(Z).LT.PI/3 AND THE EXPONENTIAL GROWTH IN
C         PI/3.LT.ABS(ARG(Z)).LT.PI WHERE ZTA=(2/3)*Z*CSQRT(Z).
C
C         WHILE THE AIRY FUNCTIONS AI(Z) AND DAI(Z)/DZ ARE ANALYTIC IN
C         THE WHOLE Z PLANE, THE CORRESPONDING SCALED FUNCTIONS DEFINED
C         FOR KODE=2 HAVE A CUT ALONG THE NEGATIVE REAL AXIS.
C         DEFINTIONS AND NOTATION ARE FOUND IN THE NBS HANDBOOK OF
C         MATHEMATICAL FUNCTIONS (REF. 1).
C
C         INPUT      ZR,ZI ARE DOUBLE PRECISION
C           ZR,ZI  - Z=CMPLX(ZR,ZI)
C           ID     - ORDER OF DERIVATIVE, ID=0 OR ID=1
C           KODE   - A PARAMETER TO INDICATE THE SCALING OPTION
C                    KODE= 1  RETURNS
C                             AI=AI(Z)                ON ID=0 OR
C                             AI=DAI(Z)/DZ            ON ID=1
C                        = 2  RETURNS
C                             AI=CEXP(ZTA)*AI(Z)       ON ID=0 OR
C                             AI=CEXP(ZTA)*DAI(Z)/DZ   ON ID=1 WHERE
C                             ZTA=(2/3)*Z*CSQRT(Z)
C
C         OUTPUT     AIR,AII ARE DOUBLE PRECISION
C           AIR,AII- COMPLEX ANSWER DEPENDING ON THE CHOICES FOR ID AND
C                    KODE
C           NZ     - UNDERFLOW INDICATOR
C                    NZ= 0   , NORMAL RETURN
C                    NZ= 1   , AI=CMPLX(0.0D0,0.0D0) DUE TO UNDERFLOW IN
C                              -PI/3.LT.ARG(Z).LT.PI/3 ON KODE=1
C           IERR   - ERROR FLAG
C                    IERR=0, NORMAL RETURN - COMPUTATION COMPLETED
C                    IERR=1, INPUT ERROR   - NO COMPUTATION
C                    IERR=2, OVERFLOW      - NO COMPUTATION, REAL(ZTA)
C                            TOO LARGE ON KODE=1
C                    IERR=3, CABS(Z) LARGE      - COMPUTATION COMPLETED
C                            LOSSES OF SIGNIFCANCE BY ARGUMENT REDUCTION
C                            PRODUCE LESS THAN HALF OF MACHINE ACCURACY
C                    IERR=4, CABS(Z) TOO LARGE  - NO COMPUTATION
C                            COMPLETE LOSS OF ACCURACY BY ARGUMENT
C                            REDUCTION
C                    IERR=5, ERROR              - NO COMPUTATION,
C                            ALGORITHM TERMINATION CONDITION NOT MET
C
C***LONG DESCRIPTION
C
C         AI AND DAI ARE COMPUTED FOR CABS(Z).GT.1.0 FROM THE K BESSEL
C         FUNCTIONS BY
C
C            AI(Z)=C*SQRT(Z)*K(1/3,ZTA) , DAI(Z)=-C*Z*K(2/3,ZTA)
C                           C=1.0/(PI*SQRT(3.0))
C                            ZTA=(2/3)*Z**(3/2)
C
C         WITH THE POWER SERIES FOR CABS(Z).LE.1.0.
C
C         IN MOST COMPLEX VARIABLE COMPUTATION, ONE MUST EVALUATE ELE-
C         MENTARY FUNCTIONS. WHEN THE MAGNITUDE OF Z IS LARGE, LOSSES
C         OF SIGNIFICANCE BY ARGUMENT REDUCTION OCCUR. CONSEQUENTLY, IF
C         THE MAGNITUDE OF ZETA=(2/3)*Z**1.5 EXCEEDS U1=SQRT(0.5/UR),
C         THEN LOSSES EXCEEDING HALF PRECISION ARE LIKELY AND AN ERROR
C         FLAG IERR=3 IS TRIGGERED WHERE UR=DMAX1(D1MACH(4),1.0D-18) IS
C         DOUBLE PRECISION UNIT ROUNDOFF LIMITED TO 18 DIGITS PRECISION.
C         ALSO, IF THE MAGNITUDE OF ZETA IS LARGER THAN U2=0.5/UR, THEN
C         ALL SIGNIFICANCE IS LOST AND IERR=4. IN ORDER TO USE THE INT
C         FUNCTION, ZETA MUST BE FURTHER RESTRICTED NOT TO EXCEED THE
C         LARGEST INTEGER, U3=I1MACH(9). THUS, THE MAGNITUDE OF ZETA
C         MUST BE RESTRICTED BY MIN(U2,U3). ON 32 BIT MACHINES, U1,U2,
C         AND U3 ARE APPROXIMATELY 2.0E+3, 4.2E+6, 2.1E+9 IN SINGLE
C         PRECISION ARITHMETIC AND 1.3E+8, 1.8E+16, 2.1E+9 IN DOUBLE
C         PRECISION ARITHMETIC RESPECTIVELY. THIS MAKES U2 AND U3 LIMIT-
C         ING IN THEIR RESPECTIVE ARITHMETICS. THIS MEANS THAT THE MAG-
C         NITUDE OF Z CANNOT EXCEED 3.1E+4 IN SINGLE AND 2.1E+6 IN
C         DOUBLE PRECISION ARITHMETIC. THIS ALSO MEANS THAT ONE CAN
C         EXPECT TO RETAIN, IN THE WORST CASES ON 32 BIT MACHINES,
C         NO DIGITS IN SINGLE PRECISION AND ONLY 7 DIGITS IN DOUBLE
C         PRECISION ARITHMETIC. SIMILAR CONSIDERATIONS HOLD FOR OTHER
C         MACHINES.
C
C         THE APPROXIMATE RELATIVE ERROR IN THE MAGNITUDE OF A COMPLEX
C         BESSEL FUNCTION CAN BE EXPRESSED BY P*10**S WHERE P=MAX(UNIT
C         ROUNDOFF,1.0E-18) IS THE NOMINAL PRECISION AND 10**S REPRE-
C         SENTS THE INCREASE IN ERROR DUE TO ARGUMENT REDUCTION IN THE
C         ELEMENTARY FUNCTIONS. HERE, S=MAX(1,ABS(LOG10(CABS(Z))),
C         ABS(LOG10(FNU))) APPROXIMATELY (I.E. S=MAX(1,ABS(EXPONENT OF
C         CABS(Z),ABS(EXPONENT OF FNU)) ). HOWEVER, THE PHASE ANGLE MAY
C         HAVE ONLY ABSOLUTE ACCURACY. THIS IS MOST LIKELY TO OCCUR WHEN
C         ONE COMPONENT (IN ABSOLUTE VALUE) IS LARGER THAN THE OTHER BY
C         SEVERAL ORDERS OF MAGNITUDE. IF ONE COMPONENT IS 10**K LARGER
C         THAN THE OTHER, THEN ONE CAN EXPECT ONLY MAX(ABS(LOG10(P))-K,
C         0) SIGNIFICANT DIGITS; OR, STATED ANOTHER WAY, WHEN K EXCEEDS
C         THE EXPONENT OF P, NO SIGNIFICANT DIGITS REMAIN IN THE SMALLER
C         COMPONENT. HOWEVER, THE PHASE ANGLE RETAINS ABSOLUTE ACCURACY
C         BECAUSE, IN COMPLEX ARITHMETIC WITH PRECISION P, THE SMALLER
C         COMPONENT WILL NOT (AS A RULE) DECREASE BELOW P TIMES THE
C         MAGNITUDE OF THE LARGER COMPONENT. IN THESE EXTREME CASES,
C         THE PRINCIPAL PHASE ANGLE IS ON THE ORDER OF +P, -P, PI/2-P,
C         OR -PI/2+P.
C
C***REFERENCES  HANDBOOK OF MATHEMATICAL FUNCTIONS BY M. ABRAMOWITZ
C                 AND I. A. STEGUN, NBS AMS SERIES 55, U.S. DEPT. OF
C                 COMMERCE, 1955.
C
C               COMPUTATION OF BESSEL FUNCTIONS OF COMPLEX ARGUMENT
C                 AND LARGE ORDER BY D. E. AMOS, SAND83-0643, MAY, 1983
C
C               A SUBROUTINE PACKAGE FOR BESSEL FUNCTIONS OF A COMPLEX
C                 ARGUMENT AND NONNEGATIVE ORDER BY D. E. AMOS, SAND85-
C                 1018, MAY, 1985
C
C               A PORTABLE PACKAGE FOR BESSEL FUNCTIONS OF A COMPLEX
C                 ARGUMENT AND NONNEGATIVE ORDER BY D. E. AMOS, ACM
C                 TRANS. MATH. SOFTWARE, VOL. 12, NO. 3, SEPTEMBER 1986,
C                 PP 265-273.
C
C***ROUTINES CALLED  ZACAI,ZBKNU,ZEXP,ZSQRT,ZABS,I1MACH,D1MACH
C***END PROLOGUE  ZAIRY
C     COMPLEX AI,CONE,CSQ,CY,S1,S2,TRM1,TRM2,Z,ZTA,Z3
      EXTERNAL ZABS
      DOUBLE PRECISION AA, AD, AII, AIR, AK, ALIM, ATRM, AZ, AZ3, BK,
     * CC, CK, COEF, CONEI, CONER, CSQI, CSQR, CYI, CYR, C1, C2, DIG,
     * DK, D1, D2, ELIM, FID, FNU, PTR, RL, R1M5, SFAC, STI, STR,
     * S1I, S1R, S2I, S2R, TOL, TRM1I, TRM1R, TRM2I, TRM2R, TTH, ZEROI,
     * ZEROR, ZI, ZR, ZTAI, ZTAR, Z3I, Z3R, D1MACH, ZABS, ALAZ, BB
      INTEGER ID, IERR, IFLAG, K, KODE, K1, K2, MR, NN, NZ, I1MACH
      DIMENSION CYR(1), CYI(1)
      DATA TTH, C1, C2, COEF /6.66666666666666667D-01,
     * 3.55028053887817240D-01,2.58819403792806799D-01,
     * 1.83776298473930683D-01/
      DATA ZEROR, ZEROI, CONER, CONEI /0.0D0,0.0D0,1.0D0,0.0D0/
C***FIRST EXECUTABLE STATEMENT  ZAIRY
      IERR = 0
      NZ=0
      IF (ID.LT.0 .OR. ID.GT.1) IERR=1
      IF (KODE.LT.1 .OR. KODE.GT.2) IERR=1
      IF (IERR.NE.0) RETURN
      AZ = ZABS(ZR,ZI)
      TOL = DMAX1(D1MACH(4),1.0D-18)
      FID = DBLE(FLOAT(ID))
      IF (AZ.GT.1.0D0) GO TO 70
C-----------------------------------------------------------------------
C     POWER SERIES FOR CABS(Z).LE.1.
C-----------------------------------------------------------------------
      S1R = CONER
      S1I = CONEI
      S2R = CONER
      S2I = CONEI
      IF (AZ.LT.TOL) GO TO 170
      AA = AZ*AZ
      IF (AA.LT.TOL/AZ) GO TO 40
      TRM1R = CONER
      TRM1I = CONEI
      TRM2R = CONER
      TRM2I = CONEI
      ATRM = 1.0D0
      STR = ZR*ZR - ZI*ZI
      STI = ZR*ZI + ZI*ZR
      Z3R = STR*ZR - STI*ZI
      Z3I = STR*ZI + STI*ZR
      AZ3 = AZ*AA
      AK = 2.0D0 + FID
      BK = 3.0D0 - FID - FID
      CK = 4.0D0 - FID
      DK = 3.0D0 + FID + FID
      D1 = AK*DK
      D2 = BK*CK
      AD = DMIN1(D1,D2)
      AK = 24.0D0 + 9.0D0*FID
      BK = 30.0D0 - 9.0D0*FID
      DO 30 K=1,25
        STR = (TRM1R*Z3R-TRM1I*Z3I)/D1
        TRM1I = (TRM1R*Z3I+TRM1I*Z3R)/D1
        TRM1R = STR
        S1R = S1R + TRM1R
        S1I = S1I + TRM1I
        STR = (TRM2R*Z3R-TRM2I*Z3I)/D2
        TRM2I = (TRM2R*Z3I+TRM2I*Z3R)/D2
        TRM2R = STR
        S2R = S2R + TRM2R
        S2I = S2I + TRM2I
        ATRM = ATRM*AZ3/AD
        D1 = D1 + AK
        D2 = D2 + BK
        AD = DMIN1(D1,D2)
        IF (ATRM.LT.TOL*AD) GO TO 40
        AK = AK + 18.0D0
        BK = BK + 18.0D0
   30 CONTINUE
   40 CONTINUE
      IF (ID.EQ.1) GO TO 50
      AIR = S1R*C1 - C2*(ZR*S2R-ZI*S2I)
      AII = S1I*C1 - C2*(ZR*S2I+ZI*S2R)
      IF (KODE.EQ.1) RETURN
      CALL ZSQRT(ZR, ZI, STR, STI)
      ZTAR = TTH*(ZR*STR-ZI*STI)
      ZTAI = TTH*(ZR*STI+ZI*STR)
      CALL ZEXP(ZTAR, ZTAI, STR, STI)
      PTR = AIR*STR - AII*STI
      AII = AIR*STI + AII*STR
      AIR = PTR
      RETURN
   50 CONTINUE
      AIR = -S2R*C2
      AII = -S2I*C2
      IF (AZ.LE.TOL) GO TO 60
      STR = ZR*S1R - ZI*S1I
      STI = ZR*S1I + ZI*S1R
      CC = C1/(1.0D0+FID)
      AIR = AIR + CC*(STR*ZR-STI*ZI)
      AII = AII + CC*(STR*ZI+STI*ZR)
   60 CONTINUE
      IF (KODE.EQ.1) RETURN
      CALL ZSQRT(ZR, ZI, STR, STI)
      ZTAR = TTH*(ZR*STR-ZI*STI)
      ZTAI = TTH*(ZR*STI+ZI*STR)
      CALL ZEXP(ZTAR, ZTAI, STR, STI)
      PTR = STR*AIR - STI*AII
      AII = STR*AII + STI*AIR
      AIR = PTR
      RETURN
C-----------------------------------------------------------------------
C     CASE FOR CABS(Z).GT.1.0
C-----------------------------------------------------------------------
   70 CONTINUE
      FNU = (1.0D0+FID)/3.0D0
C-----------------------------------------------------------------------
C     SET PARAMETERS RELATED TO MACHINE CONSTANTS.
C     TOL IS THE APPROXIMATE UNIT ROUNDOFF LIMITED TO 1.0D-18.
C     ELIM IS THE APPROXIMATE EXPONENTIAL OVER- AND UNDERFLOW LIMIT.
C     EXP(-ELIM).LT.EXP(-ALIM)=EXP(-ELIM)/TOL    AND
C     EXP(ELIM).GT.EXP(ALIM)=EXP(ELIM)*TOL       ARE INTERVALS NEAR
C     UNDERFLOW AND OVERFLOW LIMITS WHERE SCALED ARITHMETIC IS DONE.
C     RL IS THE LOWER BOUNDARY OF THE ASYMPTOTIC EXPANSION FOR LARGE Z.
C     DIG = NUMBER OF BASE 10 DIGITS IN TOL = 10**(-DIG).
C-----------------------------------------------------------------------
      K1 = I1MACH(15)
      K2 = I1MACH(16)
      R1M5 = D1MACH(5)
      K = MIN0(IABS(K1),IABS(K2))
      ELIM = 2.303D0*(DBLE(FLOAT(K))*R1M5-3.0D0)
      K1 = I1MACH(14) - 1
      AA = R1M5*DBLE(FLOAT(K1))
      DIG = DMIN1(AA,18.0D0)
      AA = AA*2.303D0
      ALIM = ELIM + DMAX1(-AA,-41.45D0)
      RL = 1.2D0*DIG + 3.0D0
      ALAZ = DLOG(AZ)
C--------------------------------------------------------------------------
C     TEST FOR PROPER RANGE
C-----------------------------------------------------------------------
      AA=0.5D0/TOL
      BB=DBLE(FLOAT(I1MACH(9)))*0.5D0
      AA=DMIN1(AA,BB)
      AA=AA**TTH
      IF (AZ.GT.AA) GO TO 260
      AA=DSQRT(AA)
      IF (AZ.GT.AA) IERR=3
      CALL ZSQRT(ZR, ZI, CSQR, CSQI)
      ZTAR = TTH*(ZR*CSQR-ZI*CSQI)
      ZTAI = TTH*(ZR*CSQI+ZI*CSQR)
C-----------------------------------------------------------------------
C     RE(ZTA).LE.0 WHEN RE(Z).LT.0, ESPECIALLY WHEN IM(Z) IS SMALL
C-----------------------------------------------------------------------
      IFLAG = 0
      SFAC = 1.0D0
      AK = ZTAI
      IF (ZR.GE.0.0D0) GO TO 80
      BK = ZTAR
      CK = -DABS(BK)
      ZTAR = CK
      ZTAI = AK
   80 CONTINUE
      IF (ZI.NE.0.0D0) GO TO 90
      IF (ZR.GT.0.0D0) GO TO 90
      ZTAR = 0.0D0
      ZTAI = AK
   90 CONTINUE
      AA = ZTAR
      IF (AA.GE.0.0D0 .AND. ZR.GT.0.0D0) GO TO 110
      IF (KODE.EQ.2) GO TO 100
C-----------------------------------------------------------------------
C     OVERFLOW TEST
C-----------------------------------------------------------------------
      IF (AA.GT.(-ALIM)) GO TO 100
      AA = -AA + 0.25D0*ALAZ
      IFLAG = 1
      SFAC = TOL
      IF (AA.GT.ELIM) GO TO 270
  100 CONTINUE
C-----------------------------------------------------------------------
C     CBKNU AND CACON RETURN EXP(ZTA)*K(FNU,ZTA) ON KODE=2
C-----------------------------------------------------------------------
      MR = 1
      IF (ZI.LT.0.0D0) MR = -1
      CALL ZACAI(ZTAR, ZTAI, FNU, KODE, MR, 1, CYR, CYI, NN, RL, TOL,
     * ELIM, ALIM)
      IF (NN.LT.0) GO TO 280
      NZ = NZ + NN
      GO TO 130
  110 CONTINUE
      IF (KODE.EQ.2) GO TO 120
C-----------------------------------------------------------------------
C     UNDERFLOW TEST
C-----------------------------------------------------------------------
      IF (AA.LT.ALIM) GO TO 120
      AA = -AA - 0.25D0*ALAZ
      IFLAG = 2
      SFAC = 1.0D0/TOL
      IF (AA.LT.(-ELIM)) GO TO 210
  120 CONTINUE
      CALL ZBKNU(ZTAR, ZTAI, FNU, KODE, 1, CYR, CYI, NZ, TOL, ELIM,
     * ALIM)
  130 CONTINUE
      S1R = CYR(1)*COEF
      S1I = CYI(1)*COEF
      IF (IFLAG.NE.0) GO TO 150
      IF (ID.EQ.1) GO TO 140
      AIR = CSQR*S1R - CSQI*S1I
      AII = CSQR*S1I + CSQI*S1R
      RETURN
  140 CONTINUE
      AIR = -(ZR*S1R-ZI*S1I)
      AII = -(ZR*S1I+ZI*S1R)
      RETURN
  150 CONTINUE
      S1R = S1R*SFAC
      S1I = S1I*SFAC
      IF (ID.EQ.1) GO TO 160
      STR = S1R*CSQR - S1I*CSQI
      S1I = S1R*CSQI + S1I*CSQR
      S1R = STR
      AIR = S1R/SFAC
      AII = S1I/SFAC
      RETURN
  160 CONTINUE
      STR = -(S1R*ZR-S1I*ZI)
      S1I = -(S1R*ZI+S1I*ZR)
      S1R = STR
      AIR = S1R/SFAC
      AII = S1I/SFAC
      RETURN
  170 CONTINUE
      AA = 1.0D+3*D1MACH(1)
      S1R = ZEROR
      S1I = ZEROI
      IF (ID.EQ.1) GO TO 190
      IF (AZ.LE.AA) GO TO 180
      S1R = C2*ZR
      S1I = C2*ZI
  180 CONTINUE
      AIR = C1 - S1R
      AII = -S1I
      RETURN
  190 CONTINUE
      AIR = -C2
      AII = 0.0D0
      AA = DSQRT(AA)
      IF (AZ.LE.AA) GO TO 200
      S1R = 0.5D0*(ZR*ZR-ZI*ZI)
      S1I = ZR*ZI
  200 CONTINUE
      AIR = AIR + C1*S1R
      AII = AII + C1*S1I
      RETURN
  210 CONTINUE
      NZ = 1
      AIR = ZEROR
      AII = ZEROI
      RETURN
  270 CONTINUE
      NZ = 0
      IERR=2
      RETURN
  280 CONTINUE
      IF(NN.EQ.(-1)) GO TO 270
      NZ=0
      IERR=5
      RETURN
  260 CONTINUE
      IERR=4
      NZ=0
      RETURN
      END
      SUBROUTINE ZBIRY(ZR, ZI, ID, KODE, BIR, BII, IERR)
C***BEGIN PROLOGUE  ZBIRY
C***DATE WRITTEN   830501   (YYMMDD)
C***REVISION DATE  890801, 930101   (YYMMDD)
C***CATEGORY NO.  B5K
C***KEYWORDS  AIRY FUNCTION,BESSEL FUNCTIONS OF ORDER ONE THIRD
C***AUTHOR  AMOS, DONALD E., SANDIA NATIONAL LABORATORIES
C***PURPOSE  TO COMPUTE AIRY FUNCTIONS BI(Z) AND DBI(Z) FOR COMPLEX Z
C***DESCRIPTION
C
C                      ***A DOUBLE PRECISION ROUTINE***
C         ON KODE=1, CBIRY COMPUTES THE COMPLEX AIRY FUNCTION BI(Z) OR
C         ITS DERIVATIVE DBI(Z)/DZ ON ID=0 OR ID=1 RESPECTIVELY. ON
C         KODE=2, A SCALING OPTION CEXP(-AXZTA)*BI(Z) OR CEXP(-AXZTA)*
C         DBI(Z)/DZ IS PROVIDED TO REMOVE THE EXPONENTIAL BEHAVIOR IN
C         BOTH THE LEFT AND RIGHT HALF PLANES WHERE
C         ZTA=(2/3)*Z*CSQRT(Z)=CMPLX(XZTA,YZTA) AND AXZTA=ABS(XZTA).
C         DEFINTIONS AND NOTATION ARE FOUND IN THE NBS HANDBOOK OF
C         MATHEMATICAL FUNCTIONS (REF. 1).
C
C         INPUT      ZR,ZI ARE DOUBLE PRECISION
C           ZR,ZI  - Z=CMPLX(ZR,ZI)
C           ID     - ORDER OF DERIVATIVE, ID=0 OR ID=1
C           KODE   - A PARAMETER TO INDICATE THE SCALING OPTION
C                    KODE= 1  RETURNS
C                             BI=BI(Z)                 ON ID=0 OR
C                             BI=DBI(Z)/DZ             ON ID=1
C                        = 2  RETURNS
C                             BI=CEXP(-AXZTA)*BI(Z)     ON ID=0 OR
C                             BI=CEXP(-AXZTA)*DBI(Z)/DZ ON ID=1 WHERE
C                             ZTA=(2/3)*Z*CSQRT(Z)=CMPLX(XZTA,YZTA)
C                             AND AXZTA=ABS(XZTA)
C
C         OUTPUT     BIR,BII ARE DOUBLE PRECISION
C           BIR,BII- COMPLEX ANSWER DEPENDING ON THE CHOICES FOR ID AND
C                    KODE
C           IERR   - ERROR FLAG
C                    IERR=0, NORMAL RETURN - COMPUTATION COMPLETED
C                    IERR=1, INPUT ERROR   - NO COMPUTATION
C                    IERR=2, OVERFLOW      - NO COMPUTATION, REAL(Z)
C                            TOO LARGE ON KODE=1
C                    IERR=3, CABS(Z) LARGE      - COMPUTATION COMPLETED
C                            LOSSES OF SIGNIFCANCE BY ARGUMENT REDUCTION
C                            PRODUCE LESS THAN HALF OF MACHINE ACCURACY
C                    IERR=4, CABS(Z) TOO LARGE  - NO COMPUTATION
C                            COMPLETE LOSS OF ACCURACY BY ARGUMENT
C                            REDUCTION
C                    IERR=5, ERROR              - NO COMPUTATION,
C                            ALGORITHM TERMINATION CONDITION NOT MET
C
C***LONG DESCRIPTION
C
C         BI AND DBI ARE COMPUTED FOR CABS(Z).GT.1.0 FROM THE I BESSEL
C         FUNCTIONS BY
C
C                BI(Z)=C*SQRT(Z)*( I(-1/3,ZTA) + I(1/3,ZTA) )
C               DBI(Z)=C *  Z  * ( I(-2/3,ZTA) + I(2/3,ZTA) )
C                               C=1.0/SQRT(3.0)
C                             ZTA=(2/3)*Z**(3/2)
C
C         WITH THE POWER SERIES FOR CABS(Z).LE.1.0.
C
C         IN MOST COMPLEX VARIABLE COMPUTATION, ONE MUST EVALUATE ELE-
C         MENTARY FUNCTIONS. WHEN THE MAGNITUDE OF Z IS LARGE, LOSSES
C         OF SIGNIFICANCE BY ARGUMENT REDUCTION OCCUR. CONSEQUENTLY, IF
C         THE MAGNITUDE OF ZETA=(2/3)*Z**1.5 EXCEEDS U1=SQRT(0.5/UR),
C         THEN LOSSES EXCEEDING HALF PRECISION ARE LIKELY AND AN ERROR
C         FLAG IERR=3 IS TRIGGERED WHERE UR=DMAX1(D1MACH(4),1.0D-18) IS
C         DOUBLE PRECISION UNIT ROUNDOFF LIMITED TO 18 DIGITS PRECISION.
C         ALSO, IF THE MAGNITUDE OF ZETA IS LARGER THAN U2=0.5/UR, THEN
C         ALL SIGNIFICANCE IS LOST AND IERR=4. IN ORDER TO USE THE INT
C         FUNCTION, ZETA MUST BE FURTHER RESTRICTED NOT TO EXCEED THE
C         LARGEST INTEGER, U3=I1MACH(9). THUS, THE MAGNITUDE OF ZETA
C         MUST BE RESTRICTED BY MIN(U2,U3). ON 32 BIT MACHINES, U1,U2,
C         AND U3 ARE APPROXIMATELY 2.0E+3, 4.2E+6, 2.1E+9 IN SINGLE
C         PRECISION ARITHMETIC AND 1.3E+8, 1.8E+16, 2.1E+9 IN DOUBLE
C         PRECISION ARITHMETIC RESPECTIVELY. THIS MAKES U2 AND U3 LIMIT-
C         ING IN THEIR RESPECTIVE ARITHMETICS. THIS MEANS THAT THE MAG-
C         NITUDE OF Z CANNOT EXCEED 3.1E+4 IN SINGLE AND 2.1E+6 IN
C         DOUBLE PRECISION ARITHMETIC. THIS ALSO MEANS THAT ONE CAN
C         EXPECT TO RETAIN, IN THE WORST CASES ON 32 BIT MACHINES,
C         NO DIGITS IN SINGLE PRECISION AND ONLY 7 DIGITS IN DOUBLE
C         PRECISION ARITHMETIC. SIMILAR CONSIDERATIONS HOLD FOR OTHER
C         MACHINES.
C
C         THE APPROXIMATE RELATIVE ERROR IN THE MAGNITUDE OF A COMPLEX
C         BESSEL FUNCTION CAN BE EXPRESSED BY P*10**S WHERE P=MAX(UNIT
C         ROUNDOFF,1.0E-18) IS THE NOMINAL PRECISION AND 10**S REPRE-
C         SENTS THE INCREASE IN ERROR DUE TO ARGUMENT REDUCTION IN THE
C         ELEMENTARY FUNCTIONS. HERE, S=MAX(1,ABS(LOG10(CABS(Z))),
C         ABS(LOG10(FNU))) APPROXIMATELY (I.E. S=MAX(1,ABS(EXPONENT OF
C         CABS(Z),ABS(EXPONENT OF FNU)) ). HOWEVER, THE PHASE ANGLE MAY
C         HAVE ONLY ABSOLUTE ACCURACY. THIS IS MOST LIKELY TO OCCUR WHEN
C         ONE COMPONENT (IN ABSOLUTE VALUE) IS LARGER THAN THE OTHER BY
C         SEVERAL ORDERS OF MAGNITUDE. IF ONE COMPONENT IS 10**K LARGER
C         THAN THE OTHER, THEN ONE CAN EXPECT ONLY MAX(ABS(LOG10(P))-K,
C         0) SIGNIFICANT DIGITS; OR, STATED ANOTHER WAY, WHEN K EXCEEDS
C         THE EXPONENT OF P, NO SIGNIFICANT DIGITS REMAIN IN THE SMALLER
C         COMPONENT. HOWEVER, THE PHASE ANGLE RETAINS ABSOLUTE ACCURACY
C         BECAUSE, IN COMPLEX ARITHMETIC WITH PRECISION P, THE SMALLER
C         COMPONENT WILL NOT (AS A RULE) DECREASE BELOW P TIMES THE
C         MAGNITUDE OF THE LARGER COMPONENT. IN THESE EXTREME CASES,
C         THE PRINCIPAL PHASE ANGLE IS ON THE ORDER OF +P, -P, PI/2-P,
C         OR -PI/2+P.
C
C***REFERENCES  HANDBOOK OF MATHEMATICAL FUNCTIONS BY M. ABRAMOWITZ
C                 AND I. A. STEGUN, NBS AMS SERIES 55, U.S. DEPT. OF
C                 COMMERCE, 1955.
C
C               COMPUTATION OF BESSEL FUNCTIONS OF COMPLEX ARGUMENT
C                 AND LARGE ORDER BY D. E. AMOS, SAND83-0643, MAY, 1983
C
C               A SUBROUTINE PACKAGE FOR BESSEL FUNCTIONS OF A COMPLEX
C                 ARGUMENT AND NONNEGATIVE ORDER BY D. E. AMOS, SAND85-
C                 1018, MAY, 1985
C
C               A PORTABLE PACKAGE FOR BESSEL FUNCTIONS OF A COMPLEX
C                 ARGUMENT AND NONNEGATIVE ORDER BY D. E. AMOS, ACM
C                 TRANS. MATH. SOFTWARE, VOL. 12, NO. 3, SEPTEMBER 1986,
C                 PP 265-273.
C
C***ROUTINES CALLED  ZBINU,ZABS,ZDIV,ZSQRT,D1MACH,I1MACH
C***END PROLOGUE  ZBIRY
C     COMPLEX BI,CONE,CSQ,CY,S1,S2,TRM1,TRM2,Z,ZTA,Z3
      EXTERNAL ZABS
      DOUBLE PRECISION AA, AD, AK, ALIM, ATRM, AZ, AZ3, BB, BII, BIR,
     * BK, CC, CK, COEF, CONEI, CONER, CSQI, CSQR, CYI, CYR, C1, C2,
     * DIG, DK, D1, D2, EAA, ELIM, FID, FMR, FNU, FNUL, PI, RL, R1M5,
     * SFAC, STI, STR, S1I, S1R, S2I, S2R, TOL, TRM1I, TRM1R, TRM2I,
     * TRM2R, TTH, ZI, ZR, ZTAI, ZTAR, Z3I, Z3R, D1MACH, ZABS
      INTEGER ID, IERR, K, KODE, K1, K2, NZ, I1MACH
      DIMENSION CYR(2), CYI(2)
      DATA TTH, C1, C2, COEF, PI /6.66666666666666667D-01,
     * 6.14926627446000736D-01,4.48288357353826359D-01,
     * 5.77350269189625765D-01,3.14159265358979324D+00/
      DATA CONER, CONEI /1.0D0,0.0D0/
C***FIRST EXECUTABLE STATEMENT  ZBIRY
      IERR = 0
      NZ=0
      IF (ID.LT.0 .OR. ID.GT.1) IERR=1
      IF (KODE.LT.1 .OR. KODE.GT.2) IERR=1
      IF (IERR.NE.0) RETURN
      AZ = ZABS(ZR,ZI)
      TOL = DMAX1(D1MACH(4),1.0D-18)
      FID = DBLE(FLOAT(ID))
      IF (AZ.GT.1.0E0) GO TO 70
C-----------------------------------------------------------------------
C     POWER SERIES FOR CABS(Z).LE.1.
C-----------------------------------------------------------------------
      S1R = CONER
      S1I = CONEI
      S2R = CONER
      S2I = CONEI
      IF (AZ.LT.TOL) GO TO 130
      AA = AZ*AZ
      IF (AA.LT.TOL/AZ) GO TO 40
      TRM1R = CONER
      TRM1I = CONEI
      TRM2R = CONER
      TRM2I = CONEI
      ATRM = 1.0D0
      STR = ZR*ZR - ZI*ZI
      STI = ZR*ZI + ZI*ZR
      Z3R = STR*ZR - STI*ZI
      Z3I = STR*ZI + STI*ZR
      AZ3 = AZ*AA
      AK = 2.0D0 + FID
      BK = 3.0D0 - FID - FID
      CK = 4.0D0 - FID
      DK = 3.0D0 + FID + FID
      D1 = AK*DK
      D2 = BK*CK
      AD = DMIN1(D1,D2)
      AK = 24.0D0 + 9.0D0*FID
      BK = 30.0D0 - 9.0D0*FID
      DO 30 K=1,25
        STR = (TRM1R*Z3R-TRM1I*Z3I)/D1
        TRM1I = (TRM1R*Z3I+TRM1I*Z3R)/D1
        TRM1R = STR
        S1R = S1R + TRM1R
        S1I = S1I + TRM1I
        STR = (TRM2R*Z3R-TRM2I*Z3I)/D2
        TRM2I = (TRM2R*Z3I+TRM2I*Z3R)/D2
        TRM2R = STR
        S2R = S2R + TRM2R
        S2I = S2I + TRM2I
        ATRM = ATRM*AZ3/AD
        D1 = D1 + AK
        D2 = D2 + BK
        AD = DMIN1(D1,D2)
        IF (ATRM.LT.TOL*AD) GO TO 40
        AK = AK + 18.0D0
        BK = BK + 18.0D0
   30 CONTINUE
   40 CONTINUE
      IF (ID.EQ.1) GO TO 50
      BIR = C1*S1R + C2*(ZR*S2R-ZI*S2I)
      BII = C1*S1I + C2*(ZR*S2I+ZI*S2R)
      IF (KODE.EQ.1) RETURN
      CALL ZSQRT(ZR, ZI, STR, STI)
      ZTAR = TTH*(ZR*STR-ZI*STI)
      ZTAI = TTH*(ZR*STI+ZI*STR)
      AA = ZTAR
      AA = -DABS(AA)
      EAA = DEXP(AA)
      BIR = BIR*EAA
      BII = BII*EAA
      RETURN
   50 CONTINUE
      BIR = S2R*C2
      BII = S2I*C2
      IF (AZ.LE.TOL) GO TO 60
      CC = C1/(1.0D0+FID)
      STR = S1R*ZR - S1I*ZI
      STI = S1R*ZI + S1I*ZR
      BIR = BIR + CC*(STR*ZR-STI*ZI)
      BII = BII + CC*(STR*ZI+STI*ZR)
   60 CONTINUE
      IF (KODE.EQ.1) RETURN
      CALL ZSQRT(ZR, ZI, STR, STI)
      ZTAR = TTH*(ZR*STR-ZI*STI)
      ZTAI = TTH*(ZR*STI+ZI*STR)
      AA = ZTAR
      AA = -DABS(AA)
      EAA = DEXP(AA)
      BIR = BIR*EAA
      BII = BII*EAA
      RETURN
C-----------------------------------------------------------------------
C     CASE FOR CABS(Z).GT.1.0
C-----------------------------------------------------------------------
   70 CONTINUE
      FNU = (1.0D0+FID)/3.0D0
C-----------------------------------------------------------------------
C     SET PARAMETERS RELATED TO MACHINE CONSTANTS.
C     TOL IS THE APPROXIMATE UNIT ROUNDOFF LIMITED TO 1.0E-18.
C     ELIM IS THE APPROXIMATE EXPONENTIAL OVER- AND UNDERFLOW LIMIT.
C     EXP(-ELIM).LT.EXP(-ALIM)=EXP(-ELIM)/TOL    AND
C     EXP(ELIM).GT.EXP(ALIM)=EXP(ELIM)*TOL       ARE INTERVALS NEAR
C     UNDERFLOW AND OVERFLOW LIMITS WHERE SCALED ARITHMETIC IS DONE.
C     RL IS THE LOWER BOUNDARY OF THE ASYMPTOTIC EXPANSION FOR LARGE Z.
C     DIG = NUMBER OF BASE 10 DIGITS IN TOL = 10**(-DIG).
C     FNUL IS THE LOWER BOUNDARY OF THE ASYMPTOTIC SERIES FOR LARGE FNU.
C-----------------------------------------------------------------------
      K1 = I1MACH(15)
      K2 = I1MACH(16)
      R1M5 = D1MACH(5)
      K = MIN0(IABS(K1),IABS(K2))
      ELIM = 2.303D0*(DBLE(FLOAT(K))*R1M5-3.0D0)
      K1 = I1MACH(14) - 1
      AA = R1M5*DBLE(FLOAT(K1))
      DIG = DMIN1(AA,18.0D0)
      AA = AA*2.303D0
      ALIM = ELIM + DMAX1(-AA,-41.45D0)
      RL = 1.2D0*DIG + 3.0D0
      FNUL = 10.0D0 + 6.0D0*(DIG-3.0D0)
C-----------------------------------------------------------------------
C     TEST FOR RANGE
C-----------------------------------------------------------------------
      AA=0.5D0/TOL
      BB=DBLE(FLOAT(I1MACH(9)))*0.5D0
      AA=DMIN1(AA,BB)
      AA=AA**TTH
      IF (AZ.GT.AA) GO TO 260
      AA=DSQRT(AA)
      IF (AZ.GT.AA) IERR=3
      CALL ZSQRT(ZR, ZI, CSQR, CSQI)
      ZTAR = TTH*(ZR*CSQR-ZI*CSQI)
      ZTAI = TTH*(ZR*CSQI+ZI*CSQR)
C-----------------------------------------------------------------------
C     RE(ZTA).LE.0 WHEN RE(Z).LT.0, ESPECIALLY WHEN IM(Z) IS SMALL
C-----------------------------------------------------------------------
      SFAC = 1.0D0
      AK = ZTAI
      IF (ZR.GE.0.0D0) GO TO 80
      BK = ZTAR
      CK = -DABS(BK)
      ZTAR = CK
      ZTAI = AK
   80 CONTINUE
      IF (ZI.NE.0.0D0 .OR. ZR.GT.0.0D0) GO TO 90
      ZTAR = 0.0D0
      ZTAI = AK
   90 CONTINUE
      AA = ZTAR
      IF (KODE.EQ.2) GO TO 100
C-----------------------------------------------------------------------
C     OVERFLOW TEST
C-----------------------------------------------------------------------
      BB = DABS(AA)
      IF (BB.LT.ALIM) GO TO 100
      BB = BB + 0.25D0*DLOG(AZ)
      SFAC = TOL
      IF (BB.GT.ELIM) GO TO 190
  100 CONTINUE
      FMR = 0.0D0
      IF (AA.GE.0.0D0 .AND. ZR.GT.0.0D0) GO TO 110
      FMR = PI
      IF (ZI.LT.0.0D0) FMR = -PI
      ZTAR = -ZTAR
      ZTAI = -ZTAI
  110 CONTINUE
C-----------------------------------------------------------------------
C     AA=FACTOR FOR ANALYTIC CONTINUATION OF I(FNU,ZTA)
C     KODE=2 RETURNS EXP(-ABS(XZTA))*I(FNU,ZTA) FROM ZBESI
C-----------------------------------------------------------------------
      CALL ZBINU(ZTAR, ZTAI, FNU, KODE, 1, CYR, CYI, NZ, RL, FNUL, TOL,
     * ELIM, ALIM)
      IF (NZ.LT.0) GO TO 200
      AA = FMR*FNU
      Z3R = SFAC
      STR = DCOS(AA)
      STI = DSIN(AA)
      S1R = (STR*CYR(1)-STI*CYI(1))*Z3R
      S1I = (STR*CYI(1)+STI*CYR(1))*Z3R
      FNU = (2.0D0-FID)/3.0D0
      CALL ZBINU(ZTAR, ZTAI, FNU, KODE, 2, CYR, CYI, NZ, RL, FNUL, TOL,
     * ELIM, ALIM)
      CYR(1) = CYR(1)*Z3R
      CYI(1) = CYI(1)*Z3R
      CYR(2) = CYR(2)*Z3R
      CYI(2) = CYI(2)*Z3R
C-----------------------------------------------------------------------
C     BACKWARD RECUR ONE STEP FOR ORDERS -1/3 OR -2/3
C-----------------------------------------------------------------------
      CALL ZDIV(CYR(1), CYI(1), ZTAR, ZTAI, STR, STI)
      S2R = (FNU+FNU)*STR + CYR(2)
      S2I = (FNU+FNU)*STI + CYI(2)
      AA = FMR*(FNU-1.0D0)
      STR = DCOS(AA)
      STI = DSIN(AA)
      S1R = COEF*(S1R+S2R*STR-S2I*STI)
      S1I = COEF*(S1I+S2R*STI+S2I*STR)
      IF (ID.EQ.1) GO TO 120
      STR = CSQR*S1R - CSQI*S1I
      S1I = CSQR*S1I + CSQI*S1R
      S1R = STR
      BIR = S1R/SFAC
      BII = S1I/SFAC
      RETURN
  120 CONTINUE
      STR = ZR*S1R - ZI*S1I
      S1I = ZR*S1I + ZI*S1R
      S1R = STR
      BIR = S1R/SFAC
      BII = S1I/SFAC
      RETURN
  130 CONTINUE
      AA = C1*(1.0D0-FID) + FID*C2
      BIR = AA
      BII = 0.0D0
      RETURN
  190 CONTINUE
      IERR=2
      NZ=0
      RETURN
  200 CONTINUE
      IF(NZ.EQ.(-1)) GO TO 190
      NZ=0
      IERR=5
      RETURN
  260 CONTINUE
      IERR=4
      NZ=0
      RETURN
      END
      SUBROUTINE ZMLT(AR, AI, BR, BI, CR, CI)
C***BEGIN PROLOGUE  ZMLT
C***REFER TO  ZBESH,ZBESI,ZBESJ,ZBESK,ZBESY,ZAIRY,ZBIRY
C
C     DOUBLE PRECISION COMPLEX MULTIPLY, C=A*B.
C
C***ROUTINES CALLED  (NONE)
C***END PROLOGUE  ZMLT
      DOUBLE PRECISION AR, AI, BR, BI, CR, CI, CA, CB
      CA = AR*BR - AI*BI
      CB = AR*BI + AI*BR
      CR = CA
      CI = CB
      RETURN
      END
      SUBROUTINE ZDIV(AR, AI, BR, BI, CR, CI)
C***BEGIN PROLOGUE  ZDIV
C***REFER TO  ZBESH,ZBESI,ZBESJ,ZBESK,ZBESY,ZAIRY,ZBIRY
C
C     DOUBLE PRECISION COMPLEX DIVIDE C=A/B.
C
C***ROUTINES CALLED  ZABS
C***END PROLOGUE  ZDIV
      EXTERNAL ZABS
      DOUBLE PRECISION AR, AI, BR, BI, CR, CI, BM, CA, CB, CC, CD
      DOUBLE PRECISION ZABS
      BM = 1.0D0/ZABS(BR,BI)
      CC = BR*BM
      CD = BI*BM
      CA = (AR*CC+AI*CD)*BM
      CB = (AI*CC-AR*CD)*BM
      CR = CA
      CI = CB
      RETURN
      END
      SUBROUTINE ZSQRT(AR, AI, BR, BI)
C***BEGIN PROLOGUE  ZSQRT
C***REFER TO  ZBESH,ZBESI,ZBESJ,ZBESK,ZBESY,ZAIRY,ZBIRY
C
C     DOUBLE PRECISION COMPLEX SQUARE ROOT, B=CSQRT(A)
C
C***ROUTINES CALLED  ZABS
C***END PROLOGUE  ZSQRT
      EXTERNAL ZABS
      DOUBLE PRECISION AR, AI, BR, BI, ZM, DTHETA, DPI, DRT
      DOUBLE PRECISION ZABS
      DATA DRT , DPI / 7.071067811865475244008443621D-1,
     1                 3.141592653589793238462643383D+0/
      ZM = ZABS(AR,AI)
      ZM = DSQRT(ZM)
      IF (AR.EQ.0.0D+0) GO TO 10
      IF (AI.EQ.0.0D+0) GO TO 20
      DTHETA = DATAN(AI/AR)
      IF (DTHETA.LE.0.0D+0) GO TO 40
      IF (AR.LT.0.0D+0) DTHETA = DTHETA - DPI
      GO TO 50
   10 IF (AI.GT.0.0D+0) GO TO 60
      IF (AI.LT.0.0D+0) GO TO 70
      BR = 0.0D+0
      BI = 0.0D+0
      RETURN
   20 IF (AR.GT.0.0D+0) GO TO 30
      BR = 0.0D+0
      BI = DSQRT(DABS(AR))
      RETURN
   30 BR = DSQRT(AR)
      BI = 0.0D+0
      RETURN
   40 IF (AR.LT.0.0D+0) DTHETA = DTHETA + DPI
   50 DTHETA = DTHETA*0.5D+0
      BR = ZM*DCOS(DTHETA)
      BI = ZM*DSIN(DTHETA)
      RETURN
   60 BR = ZM*DRT
      BI = ZM*DRT
      RETURN
   70 BR = ZM*DRT
      BI = -ZM*DRT
      RETURN
      END
      SUBROUTINE ZEXP(AR, AI, BR, BI)
C***BEGIN PROLOGUE  ZEXP
C***REFER TO  ZBESH,ZBESI,ZBESJ,ZBESK,ZBESY,ZAIRY,ZBIRY
C
C     DOUBLE PRECISION COMPLEX EXPONENTIAL FUNCTION B=EXP(A)
C
C***ROUTINES CALLED  (NONE)
C***END PROLOGUE  ZEXP
      DOUBLE PRECISION AR, AI, BR, BI, ZM, CA, CB
      ZM = DEXP(AR)
      CA = ZM*DCOS(AI)
      CB = ZM*DSIN(AI)
      BR = CA
      BI = CB
      RETURN
      END
      SUBROUTINE ZLOG(AR, AI, BR, BI, IERR)
C***BEGIN PROLOGUE  ZLOG
C***REFER TO  ZBESH,ZBESI,ZBESJ,ZBESK,ZBESY,ZAIRY,ZBIRY
C
C     DOUBLE PRECISION COMPLEX LOGARITHM B=CLOG(A)
C     IERR=0,NORMAL RETURN      IERR=1, Z=CMPLX(0.0,0.0)
C***ROUTINES CALLED  ZABS
C***END PROLOGUE  ZLOG
      EXTERNAL ZABS
      DOUBLE PRECISION AR, AI, BR, BI, ZM, DTHETA, DPI, DHPI
      DOUBLE PRECISION ZABS
      INTEGER IERR
      DATA DPI , DHPI  / 3.141592653589793238462643383D+0,
     1                   1.570796326794896619231321696D+0/
C
      IERR=0
      IF (AR.EQ.0.0D+0) GO TO 10
      IF (AI.EQ.0.0D+0) GO TO 20
      DTHETA = DATAN(AI/AR)
      IF (DTHETA.LE.0.0D+0) GO TO 40
      IF (AR.LT.0.0D+0) DTHETA = DTHETA - DPI
      GO TO 50
   10 IF (AI.EQ.0.0D+0) GO TO 60
      BI = DHPI
      BR = DLOG(DABS(AI))
      IF (AI.LT.0.0D+0) BI = -BI
      RETURN
   20 IF (AR.GT.0.0D+0) GO TO 30
      BR = DLOG(DABS(AR))
      BI = DPI
      RETURN
   30 BR = DLOG(AR)
      BI = 0.0D+0
      RETURN
   40 IF (AR.LT.0.0D+0) DTHETA = DTHETA + DPI
   50 ZM = ZABS(AR,AI)
      BR = DLOG(ZM)
      BI = DTHETA
      RETURN
   60 CONTINUE
      IERR=1
      RETURN
      END
      DOUBLE PRECISION FUNCTION ZABS(ZR, ZI)
C***BEGIN PROLOGUE  ZABS
C***REFER TO  ZBESH,ZBESI,ZBESJ,ZBESK,ZBESY,ZAIRY,ZBIRY
C
C     ZABS COMPUTES THE ABSOLUTE VALUE OR MAGNITUDE OF A DOUBLE
C     PRECISION COMPLEX VARIABLE CMPLX(ZR,ZI)
C
C***ROUTINES CALLED  (NONE)
C***END PROLOGUE  ZABS
      DOUBLE PRECISION ZR, ZI, U, V, Q, S
      U = DABS(ZR)
      V = DABS(ZI)
      S = U + V
C-----------------------------------------------------------------------
C     S*1.0D0 MAKES AN UNNORMALIZED UNDERFLOW ON CDC MACHINES INTO A
C     TRUE FLOATING ZERO
C-----------------------------------------------------------------------
      S = S*1.0D+0
      IF (S.EQ.0.0D+0) GO TO 20
      IF (U.GT.V) GO TO 10
      Q = U/V
      ZABS = V*DSQRT(1.D+0+Q*Q)
      RETURN
   10 Q = V/U
      ZABS = U*DSQRT(1.D+0+Q*Q)
      RETURN
   20 ZABS = 0.0D+0
      RETURN
      END
      SUBROUTINE ZBKNU(ZR, ZI, FNU, KODE, N, YR, YI, NZ, TOL, ELIM,
     * ALIM)
C***BEGIN PROLOGUE  ZBKNU
C***REFER TO  ZBESI,ZBESK,ZAIRY,ZBESH
C
C     ZBKNU COMPUTES THE K BESSEL FUNCTION IN THE RIGHT HALF Z PLANE.
C
C***ROUTINES CALLED  DGAMLN,I1MACH,D1MACH,ZKSCL,ZSHCH,ZUCHK,ZABS,ZDIV,
C                    ZEXP,ZLOG,ZMLT,ZSQRT
C***END PROLOGUE  ZBKNU
C
      EXTERNAL ZABS
      DOUBLE PRECISION AA, AK, ALIM, ASCLE, A1, A2, BB, BK, BRY, CAZ,
     * CBI, CBR, CC, CCHI, CCHR, CKI, CKR, COEFI, COEFR, CONEI, CONER,
     * CRSCR, CSCLR, CSHI, CSHR, CSI, CSR, CSRR, CSSR, CTWOR,
     * CZEROI, CZEROR, CZI, CZR, DNU, DNU2, DPI, ELIM, ETEST, FC, FHS,
     * FI, FK, FKS, FMUI, FMUR, FNU, FPI, FR, G1, G2, HPI, PI, PR, PTI,
     * PTR, P1I, P1R, P2I, P2M, P2R, QI, QR, RAK, RCAZ, RTHPI, RZI,
     * RZR, R1, S, SMUI, SMUR, SPI, STI, STR, S1I, S1R, S2I, S2R, TM,
     * TOL, TTH, T1, T2, YI, YR, ZI, ZR, DGAMLN, D1MACH, ZABS, ELM,
     * CELMR, ZDR, ZDI, AS, ALAS, HELIM, CYR, CYI
      INTEGER I, IFLAG, INU, K, KFLAG, KK, KMAX, KODE, KODED, N, NZ,
     * IDUM, I1MACH, J, IC, INUB, NW
      DIMENSION YR(N), YI(N), CC(8), CSSR(3), CSRR(3), BRY(3), CYR(2),
     * CYI(2)
C     COMPLEX Z,Y,A,B,RZ,SMU,FU,FMU,F,FLRZ,CZ,S1,S2,CSH,CCH
C     COMPLEX CK,P,Q,COEF,P1,P2,CBK,PT,CZERO,CONE,CTWO,ST,EZ,CS,DK
C
      DATA KMAX / 30 /
      DATA CZEROR,CZEROI,CONER,CONEI,CTWOR,R1/
     1  0.0D0 , 0.0D0 , 1.0D0 , 0.0D0 , 2.0D0 , 2.0D0 /
      DATA DPI, RTHPI, SPI ,HPI, FPI, TTH /
     1     3.14159265358979324D0,       1.25331413731550025D0,
     2     1.90985931710274403D0,       1.57079632679489662D0,
     3     1.89769999331517738D0,       6.66666666666666666D-01/
      DATA CC(1), CC(2), CC(3), CC(4), CC(5), CC(6), CC(7), CC(8)/
     1     5.77215664901532861D-01,    -4.20026350340952355D-02,
     2    -4.21977345555443367D-02,     7.21894324666309954D-03,
     3    -2.15241674114950973D-04,    -2.01348547807882387D-05,
     4     1.13302723198169588D-06,     6.11609510448141582D-09/
C
      CAZ = ZABS(ZR,ZI)
      CSCLR = 1.0D0/TOL
      CRSCR = TOL
      CSSR(1) = CSCLR
      CSSR(2) = 1.0D0
      CSSR(3) = CRSCR
      CSRR(1) = CRSCR
      CSRR(2) = 1.0D0
      CSRR(3) = CSCLR
      BRY(1) = 1.0D+3*D1MACH(1)/TOL
      BRY(2) = 1.0D0/BRY(1)
      BRY(3) = D1MACH(2)
      NZ = 0
      IFLAG = 0
      KODED = KODE
      RCAZ = 1.0D0/CAZ
      STR = ZR*RCAZ
      STI = -ZI*RCAZ
      RZR = (STR+STR)*RCAZ
      RZI = (STI+STI)*RCAZ
      INU = INT(FNU+0.5D0)
      DNU = FNU - DBLE(FLOAT(INU))
      IF (DABS(DNU).EQ.0.5D0) GO TO 110
      DNU2 = 0.0D0
      IF (DABS(DNU).GT.TOL) DNU2 = DNU*DNU
      IF (CAZ.GT.R1) GO TO 110
C-----------------------------------------------------------------------
C     SERIES FOR CABS(Z).LE.R1
C-----------------------------------------------------------------------
      FC = 1.0D0
      CALL ZLOG(RZR, RZI, SMUR, SMUI, IDUM)
      FMUR = SMUR*DNU
      FMUI = SMUI*DNU
      CALL ZSHCH(FMUR, FMUI, CSHR, CSHI, CCHR, CCHI)
      IF (DNU.EQ.0.0D0) GO TO 10
      FC = DNU*DPI
      FC = FC/DSIN(FC)
      SMUR = CSHR/DNU
      SMUI = CSHI/DNU
   10 CONTINUE
      A2 = 1.0D0 + DNU
C-----------------------------------------------------------------------
C     GAM(1-Z)*GAM(1+Z)=PI*Z/SIN(PI*Z), T1=1/GAM(1-DNU), T2=1/GAM(1+DNU)
C-----------------------------------------------------------------------
      T2 = DEXP(-DGAMLN(A2,IDUM))
      T1 = 1.0D0/(T2*FC)
      IF (DABS(DNU).GT.0.1D0) GO TO 40
C-----------------------------------------------------------------------
C     SERIES FOR F0 TO RESOLVE INDETERMINACY FOR SMALL ABS(DNU)
C-----------------------------------------------------------------------
      AK = 1.0D0
      S = CC(1)
      DO 20 K=2,8
        AK = AK*DNU2
        TM = CC(K)*AK
        S = S + TM
        IF (DABS(TM).LT.TOL) GO TO 30
   20 CONTINUE
   30 G1 = -S
      GO TO 50
   40 CONTINUE
      G1 = (T1-T2)/(DNU+DNU)
   50 CONTINUE
      G2 = (T1+T2)*0.5D0
      FR = FC*(CCHR*G1+SMUR*G2)
      FI = FC*(CCHI*G1+SMUI*G2)
      CALL ZEXP(FMUR, FMUI, STR, STI)
      PR = 0.5D0*STR/T2
      PI = 0.5D0*STI/T2
      CALL ZDIV(0.5D0, 0.0D0, STR, STI, PTR, PTI)
      QR = PTR/T1
      QI = PTI/T1
      S1R = FR
      S1I = FI
      S2R = PR
      S2I = PI
      AK = 1.0D0
      A1 = 1.0D0
      CKR = CONER
      CKI = CONEI
      BK = 1.0D0 - DNU2
      IF (INU.GT.0 .OR. N.GT.1) GO TO 80
C-----------------------------------------------------------------------
C     GENERATE K(FNU,Z), 0.0D0 .LE. FNU .LT. 0.5D0 AND N=1
C-----------------------------------------------------------------------
      IF (CAZ.LT.TOL) GO TO 70
      CALL ZMLT(ZR, ZI, ZR, ZI, CZR, CZI)
      CZR = 0.25D0*CZR
      CZI = 0.25D0*CZI
      T1 = 0.25D0*CAZ*CAZ
   60 CONTINUE
      FR = (FR*AK+PR+QR)/BK
      FI = (FI*AK+PI+QI)/BK
      STR = 1.0D0/(AK-DNU)
      PR = PR*STR
      PI = PI*STR
      STR = 1.0D0/(AK+DNU)
      QR = QR*STR
      QI = QI*STR
      STR = CKR*CZR - CKI*CZI
      RAK = 1.0D0/AK
      CKI = (CKR*CZI+CKI*CZR)*RAK
      CKR = STR*RAK
      S1R = CKR*FR - CKI*FI + S1R
      S1I = CKR*FI + CKI*FR + S1I
      A1 = A1*T1*RAK
      BK = BK + AK + AK + 1.0D0
      AK = AK + 1.0D0
      IF (A1.GT.TOL) GO TO 60
   70 CONTINUE
      YR(1) = S1R
      YI(1) = S1I
      IF (KODED.EQ.1) RETURN
      CALL ZEXP(ZR, ZI, STR, STI)
      CALL ZMLT(S1R, S1I, STR, STI, YR(1), YI(1))
      RETURN
C-----------------------------------------------------------------------
C     GENERATE K(DNU,Z) AND K(DNU+1,Z) FOR FORWARD RECURRENCE
C-----------------------------------------------------------------------
   80 CONTINUE
      IF (CAZ.LT.TOL) GO TO 100
      CALL ZMLT(ZR, ZI, ZR, ZI, CZR, CZI)
      CZR = 0.25D0*CZR
      CZI = 0.25D0*CZI
      T1 = 0.25D0*CAZ*CAZ
   90 CONTINUE
      FR = (FR*AK+PR+QR)/BK
      FI = (FI*AK+PI+QI)/BK
      STR = 1.0D0/(AK-DNU)
      PR = PR*STR
      PI = PI*STR
      STR = 1.0D0/(AK+DNU)
      QR = QR*STR
      QI = QI*STR
      STR = CKR*CZR - CKI*CZI
      RAK = 1.0D0/AK
      CKI = (CKR*CZI+CKI*CZR)*RAK
      CKR = STR*RAK
      S1R = CKR*FR - CKI*FI + S1R
      S1I = CKR*FI + CKI*FR + S1I
      STR = PR - FR*AK
      STI = PI - FI*AK
      S2R = CKR*STR - CKI*STI + S2R
      S2I = CKR*STI + CKI*STR + S2I
      A1 = A1*T1*RAK
      BK = BK + AK + AK + 1.0D0
      AK = AK + 1.0D0
      IF (A1.GT.TOL) GO TO 90
  100 CONTINUE
      KFLAG = 2
      A1 = FNU + 1.0D0
      AK = A1*DABS(SMUR)
      IF (AK.GT.ALIM) KFLAG = 3
      STR = CSSR(KFLAG)
      P2R = S2R*STR
      P2I = S2I*STR
      CALL ZMLT(P2R, P2I, RZR, RZI, S2R, S2I)
      S1R = S1R*STR
      S1I = S1I*STR
      IF (KODED.EQ.1) GO TO 210
      CALL ZEXP(ZR, ZI, FR, FI)
      CALL ZMLT(S1R, S1I, FR, FI, S1R, S1I)
      CALL ZMLT(S2R, S2I, FR, FI, S2R, S2I)
      GO TO 210
C-----------------------------------------------------------------------
C     IFLAG=0 MEANS NO UNDERFLOW OCCURRED
C     IFLAG=1 MEANS AN UNDERFLOW OCCURRED- COMPUTATION PROCEEDS WITH
C     KODED=2 AND A TEST FOR ON SCALE VALUES IS MADE DURING FORWARD
C     RECURSION
C-----------------------------------------------------------------------
  110 CONTINUE
      CALL ZSQRT(ZR, ZI, STR, STI)
      CALL ZDIV(RTHPI, CZEROI, STR, STI, COEFR, COEFI)
      KFLAG = 2
      IF (KODED.EQ.2) GO TO 120
      IF (ZR.GT.ALIM) GO TO 290
C     BLANK LINE
      STR = DEXP(-ZR)*CSSR(KFLAG)
      STI = -STR*DSIN(ZI)
      STR = STR*DCOS(ZI)
      CALL ZMLT(COEFR, COEFI, STR, STI, COEFR, COEFI)
  120 CONTINUE
      IF (DABS(DNU).EQ.0.5D0) GO TO 300
C-----------------------------------------------------------------------
C     MILLER ALGORITHM FOR CABS(Z).GT.R1
C-----------------------------------------------------------------------
      AK = DCOS(DPI*DNU)
      AK = DABS(AK)
      IF (AK.EQ.CZEROR) GO TO 300
      FHS = DABS(0.25D0-DNU2)
      IF (FHS.EQ.CZEROR) GO TO 300
C-----------------------------------------------------------------------
C     COMPUTE R2=F(E). IF CABS(Z).GE.R2, USE FORWARD RECURRENCE TO
C     DETERMINE THE BACKWARD INDEX K. R2=F(E) IS A STRAIGHT LINE ON
C     12.LE.E.LE.60. E IS COMPUTED FROM 2**(-E)=B**(1-I1MACH(14))=
C     TOL WHERE B IS THE BASE OF THE ARITHMETIC.
C-----------------------------------------------------------------------
      T1 = DBLE(FLOAT(I1MACH(14)-1))
      T1 = T1*D1MACH(5)*3.321928094D0
      T1 = DMAX1(T1,12.0D0)
      T1 = DMIN1(T1,60.0D0)
      T2 = TTH*T1 - 6.0D0
      IF (ZR.NE.0.0D0) GO TO 130
      T1 = HPI
      GO TO 140
  130 CONTINUE
      T1 = DATAN(ZI/ZR)
      T1 = DABS(T1)
  140 CONTINUE
      IF (T2.GT.CAZ) GO TO 170
C-----------------------------------------------------------------------
C     FORWARD RECURRENCE LOOP WHEN CABS(Z).GE.R2
C-----------------------------------------------------------------------
      ETEST = AK/(DPI*CAZ*TOL)
      FK = CONER
      IF (ETEST.LT.CONER) GO TO 180
      FKS = CTWOR
      CKR = CAZ + CAZ + CTWOR
      P1R = CZEROR
      P2R = CONER
      DO 150 I=1,KMAX
        AK = FHS/FKS
        CBR = CKR/(FK+CONER)
        PTR = P2R
        P2R = CBR*P2R - P1R*AK
        P1R = PTR
        CKR = CKR + CTWOR
        FKS = FKS + FK + FK + CTWOR
        FHS = FHS + FK + FK
        FK = FK + CONER
        STR = DABS(P2R)*FK
        IF (ETEST.LT.STR) GO TO 160
  150 CONTINUE
      GO TO 310
  160 CONTINUE
      FK = FK + SPI*T1*DSQRT(T2/CAZ)
      FHS = DABS(0.25D0-DNU2)
      GO TO 180
  170 CONTINUE
C-----------------------------------------------------------------------
C     COMPUTE BACKWARD INDEX K FOR CABS(Z).LT.R2
C-----------------------------------------------------------------------
      A2 = DSQRT(CAZ)
      AK = FPI*AK/(TOL*DSQRT(A2))
      AA = 3.0D0*T1/(1.0D0+CAZ)
      BB = 14.7D0*T1/(28.0D0+CAZ)
      AK = (DLOG(AK)+CAZ*DCOS(AA)/(1.0D0+0.008D0*CAZ))/DCOS(BB)
      FK = 0.12125D0*AK*AK/CAZ + 1.5D0
  180 CONTINUE
C-----------------------------------------------------------------------
C     BACKWARD RECURRENCE LOOP FOR MILLER ALGORITHM
C-----------------------------------------------------------------------
      K = INT(SNGL(FK))
      FK = DBLE(FLOAT(K))
      FKS = FK*FK
      P1R = CZEROR
      P1I = CZEROI
      P2R = TOL
      P2I = CZEROI
      CSR = P2R
      CSI = P2I
      DO 190 I=1,K
        A1 = FKS - FK
        AK = (FKS+FK)/(A1+FHS)
        RAK = 2.0D0/(FK+CONER)
        CBR = (FK+ZR)*RAK
        CBI = ZI*RAK
        PTR = P2R
        PTI = P2I
        P2R = (PTR*CBR-PTI*CBI-P1R)*AK
        P2I = (PTI*CBR+PTR*CBI-P1I)*AK
        P1R = PTR
        P1I = PTI
        CSR = CSR + P2R
        CSI = CSI + P2I
        FKS = A1 - FK + CONER
        FK = FK - CONER
  190 CONTINUE
C-----------------------------------------------------------------------
C     COMPUTE (P2/CS)=(P2/CABS(CS))*(CONJG(CS)/CABS(CS)) FOR BETTER
C     SCALING
C-----------------------------------------------------------------------
      TM = ZABS(CSR,CSI)
      PTR = 1.0D0/TM
      S1R = P2R*PTR
      S1I = P2I*PTR
      CSR = CSR*PTR
      CSI = -CSI*PTR
      CALL ZMLT(COEFR, COEFI, S1R, S1I, STR, STI)
      CALL ZMLT(STR, STI, CSR, CSI, S1R, S1I)
      IF (INU.GT.0 .OR. N.GT.1) GO TO 200
      ZDR = ZR
      ZDI = ZI
      IF(IFLAG.EQ.1) GO TO 270
      GO TO 240
  200 CONTINUE
C-----------------------------------------------------------------------
C     COMPUTE P1/P2=(P1/CABS(P2)*CONJG(P2)/CABS(P2) FOR SCALING
C-----------------------------------------------------------------------
      TM = ZABS(P2R,P2I)
      PTR = 1.0D0/TM
      P1R = P1R*PTR
      P1I = P1I*PTR
      P2R = P2R*PTR
      P2I = -P2I*PTR
      CALL ZMLT(P1R, P1I, P2R, P2I, PTR, PTI)
      STR = DNU + 0.5D0 - PTR
      STI = -PTI
      CALL ZDIV(STR, STI, ZR, ZI, STR, STI)
      STR = STR + 1.0D0
      CALL ZMLT(STR, STI, S1R, S1I, S2R, S2I)
C-----------------------------------------------------------------------
C     FORWARD RECURSION ON THE THREE TERM RECURSION WITH RELATION WITH
C     SCALING NEAR EXPONENT EXTREMES ON KFLAG=1 OR KFLAG=3
C-----------------------------------------------------------------------
  210 CONTINUE
      STR = DNU + 1.0D0
      CKR = STR*RZR
      CKI = STR*RZI
      IF (N.EQ.1) INU = INU - 1
      IF (INU.GT.0) GO TO 220
      IF (N.GT.1) GO TO 215
      S1R = S2R
      S1I = S2I
  215 CONTINUE
      ZDR = ZR
      ZDI = ZI
      IF(IFLAG.EQ.1) GO TO 270
      GO TO 240
  220 CONTINUE
      INUB = 1
      IF(IFLAG.EQ.1) GO TO 261
  225 CONTINUE
      P1R = CSRR(KFLAG)
      ASCLE = BRY(KFLAG)
      DO 230 I=INUB,INU
        STR = S2R
        STI = S2I
        S2R = CKR*STR - CKI*STI + S1R
        S2I = CKR*STI + CKI*STR + S1I
        S1R = STR
        S1I = STI
        CKR = CKR + RZR
        CKI = CKI + RZI
        IF (KFLAG.GE.3) GO TO 230
        P2R = S2R*P1R
        P2I = S2I*P1R
        STR = DABS(P2R)
        STI = DABS(P2I)
        P2M = DMAX1(STR,STI)
        IF (P2M.LE.ASCLE) GO TO 230
        KFLAG = KFLAG + 1
        ASCLE = BRY(KFLAG)
        S1R = S1R*P1R
        S1I = S1I*P1R
        S2R = P2R
        S2I = P2I
        STR = CSSR(KFLAG)
        S1R = S1R*STR
        S1I = S1I*STR
        S2R = S2R*STR
        S2I = S2I*STR
        P1R = CSRR(KFLAG)
  230 CONTINUE
      IF (N.NE.1) GO TO 240
      S1R = S2R
      S1I = S2I
  240 CONTINUE
      STR = CSRR(KFLAG)
      YR(1) = S1R*STR
      YI(1) = S1I*STR
      IF (N.EQ.1) RETURN
      YR(2) = S2R*STR
      YI(2) = S2I*STR
      IF (N.EQ.2) RETURN
      KK = 2
  250 CONTINUE
      KK = KK + 1
      IF (KK.GT.N) RETURN
      P1R = CSRR(KFLAG)
      ASCLE = BRY(KFLAG)
      DO 260 I=KK,N
        P2R = S2R
        P2I = S2I
        S2R = CKR*P2R - CKI*P2I + S1R
        S2I = CKI*P2R + CKR*P2I + S1I
        S1R = P2R
        S1I = P2I
        CKR = CKR + RZR
        CKI = CKI + RZI
        P2R = S2R*P1R
        P2I = S2I*P1R
        YR(I) = P2R
        YI(I) = P2I
        IF (KFLAG.GE.3) GO TO 260
        STR = DABS(P2R)
        STI = DABS(P2I)
        P2M = DMAX1(STR,STI)
        IF (P2M.LE.ASCLE) GO TO 260
        KFLAG = KFLAG + 1
        ASCLE = BRY(KFLAG)
        S1R = S1R*P1R
        S1I = S1I*P1R
        S2R = P2R
        S2I = P2I
        STR = CSSR(KFLAG)
        S1R = S1R*STR
        S1I = S1I*STR
        S2R = S2R*STR
        S2I = S2I*STR
        P1R = CSRR(KFLAG)
  260 CONTINUE
      RETURN
C-----------------------------------------------------------------------
C     IFLAG=1 CASES, FORWARD RECURRENCE ON SCALED VALUES ON UNDERFLOW
C-----------------------------------------------------------------------
  261 CONTINUE
      HELIM = 0.5D0*ELIM
      ELM = DEXP(-ELIM)
      CELMR = ELM
      ASCLE = BRY(1)
      ZDR = ZR
      ZDI = ZI
      IC = -1
      J = 2
      DO 262 I=1,INU
        STR = S2R
        STI = S2I
        S2R = STR*CKR-STI*CKI+S1R
        S2I = STI*CKR+STR*CKI+S1I
        S1R = STR
        S1I = STI
        CKR = CKR+RZR
        CKI = CKI+RZI
        AS = ZABS(S2R,S2I)
        ALAS = DLOG(AS)
        P2R = -ZDR+ALAS
        IF(P2R.LT.(-ELIM)) GO TO 263
        CALL ZLOG(S2R,S2I,STR,STI,IDUM)
        P2R = -ZDR+STR
        P2I = -ZDI+STI
        P2M = DEXP(P2R)/TOL
        P1R = P2M*DCOS(P2I)
        P1I = P2M*DSIN(P2I)
        CALL ZUCHK(P1R,P1I,NW,ASCLE,TOL)
        IF(NW.NE.0) GO TO 263
        J = 3 - J
        CYR(J) = P1R
        CYI(J) = P1I
        IF(IC.EQ.(I-1)) GO TO 264
        IC = I
        GO TO 262
  263   CONTINUE
        IF(ALAS.LT.HELIM) GO TO 262
        ZDR = ZDR-ELIM
        S1R = S1R*CELMR
        S1I = S1I*CELMR
        S2R = S2R*CELMR
        S2I = S2I*CELMR
  262 CONTINUE
      IF(N.NE.1) GO TO 270
      S1R = S2R
      S1I = S2I
      GO TO 270
  264 CONTINUE
      KFLAG = 1
      INUB = I+1
      S2R = CYR(J)
      S2I = CYI(J)
      J = 3 - J
      S1R = CYR(J)
      S1I = CYI(J)
      IF(INUB.LE.INU) GO TO 225
      IF(N.NE.1) GO TO 240
      S1R = S2R
      S1I = S2I
      GO TO 240
  270 CONTINUE
      YR(1) = S1R
      YI(1) = S1I
      IF(N.EQ.1) GO TO 280
      YR(2) = S2R
      YI(2) = S2I
  280 CONTINUE
      ASCLE = BRY(1)
      CALL ZKSCL(ZDR,ZDI,FNU,N,YR,YI,NZ,RZR,RZI,ASCLE,TOL,ELIM)
      INU = N - NZ
      IF (INU.LE.0) RETURN
      KK = NZ + 1
      S1R = YR(KK)
      S1I = YI(KK)
      YR(KK) = S1R*CSRR(1)
      YI(KK) = S1I*CSRR(1)
      IF (INU.EQ.1) RETURN
      KK = NZ + 2
      S2R = YR(KK)
      S2I = YI(KK)
      YR(KK) = S2R*CSRR(1)
      YI(KK) = S2I*CSRR(1)
      IF (INU.EQ.2) RETURN
      T2 = FNU + DBLE(FLOAT(KK-1))
      CKR = T2*RZR
      CKI = T2*RZI
      KFLAG = 1
      GO TO 250
  290 CONTINUE
C-----------------------------------------------------------------------
C     SCALE BY DEXP(Z), IFLAG = 1 CASES
C-----------------------------------------------------------------------
      KODED = 2
      IFLAG = 1
      KFLAG = 2
      GO TO 120
C-----------------------------------------------------------------------
C     FNU=HALF ODD INTEGER CASE, DNU=-0.5
C-----------------------------------------------------------------------
  300 CONTINUE
      S1R = COEFR
      S1I = COEFI
      S2R = COEFR
      S2I = COEFI
      GO TO 210
C
C
  310 CONTINUE
      NZ=-2
      RETURN
      END
      SUBROUTINE ZKSCL(ZRR,ZRI,FNU,N,YR,YI,NZ,RZR,RZI,ASCLE,TOL,ELIM)
C***BEGIN PROLOGUE  ZKSCL
C***REFER TO  ZBESK
C
C     SET K FUNCTIONS TO ZERO ON UNDERFLOW, CONTINUE RECURRENCE
C     ON SCALED FUNCTIONS UNTIL TWO MEMBERS COME ON SCALE, THEN
C     RETURN WITH MIN(NZ+2,N) VALUES SCALED BY 1/TOL.
C
C***ROUTINES CALLED  ZUCHK,ZABS,ZLOG
C***END PROLOGUE  ZKSCL
C     COMPLEX CK,CS,CY,CZERO,RZ,S1,S2,Y,ZR,ZD,CELM
      EXTERNAL ZABS
      DOUBLE PRECISION ACS, AS, ASCLE, CKI, CKR, CSI, CSR, CYI,
     * CYR, ELIM, FN, FNU, RZI, RZR, STR, S1I, S1R, S2I,
     * S2R, TOL, YI, YR, ZEROI, ZEROR, ZRI, ZRR, ZABS,
     * ZDR, ZDI, CELMR, ELM, HELIM, ALAS
      INTEGER I, IC, IDUM, KK, N, NN, NW, NZ
      DIMENSION YR(N), YI(N), CYR(2), CYI(2)
      DATA ZEROR,ZEROI / 0.0D0 , 0.0D0 /
C
      NZ = 0
      IC = 0
      NN = MIN0(2,N)
      DO 10 I=1,NN
        S1R = YR(I)
        S1I = YI(I)
        CYR(I) = S1R
        CYI(I) = S1I
        AS = ZABS(S1R,S1I)
        ACS = -ZRR + DLOG(AS)
        NZ = NZ + 1
        YR(I) = ZEROR
        YI(I) = ZEROI
        IF (ACS.LT.(-ELIM)) GO TO 10
        CALL ZLOG(S1R, S1I, CSR, CSI, IDUM)
        CSR = CSR - ZRR
        CSI = CSI - ZRI
        STR = DEXP(CSR)/TOL
        CSR = STR*DCOS(CSI)
        CSI = STR*DSIN(CSI)
        CALL ZUCHK(CSR, CSI, NW, ASCLE, TOL)
        IF (NW.NE.0) GO TO 10
        YR(I) = CSR
        YI(I) = CSI
        IC = I
        NZ = NZ - 1
   10 CONTINUE
      IF (N.EQ.1) RETURN
      IF (IC.GT.1) GO TO 20
      YR(1) = ZEROR
      YI(1) = ZEROI
      NZ = 2
   20 CONTINUE
      IF (N.EQ.2) RETURN
      IF (NZ.EQ.0) RETURN
      FN = FNU + 1.0D0
      CKR = FN*RZR
      CKI = FN*RZI
      S1R = CYR(1)
      S1I = CYI(1)
      S2R = CYR(2)
      S2I = CYI(2)
      HELIM = 0.5D0*ELIM
      ELM = DEXP(-ELIM)
      CELMR = ELM
      ZDR = ZRR
      ZDI = ZRI
C
C     FIND TWO CONSECUTIVE Y VALUES ON SCALE. SCALE RECURRENCE IF
C     S2 GETS LARGER THAN EXP(ELIM/2)
C
      DO 30 I=3,N
        KK = I
        CSR = S2R
        CSI = S2I
        S2R = CKR*CSR - CKI*CSI + S1R
        S2I = CKI*CSR + CKR*CSI + S1I
        S1R = CSR
        S1I = CSI
        CKR = CKR + RZR
        CKI = CKI + RZI
        AS = ZABS(S2R,S2I)
        ALAS = DLOG(AS)
        ACS = -ZDR + ALAS
        NZ = NZ + 1
        YR(I) = ZEROR
        YI(I) = ZEROI
        IF (ACS.LT.(-ELIM)) GO TO 25
        CALL ZLOG(S2R, S2I, CSR, CSI, IDUM)
        CSR = CSR - ZDR
        CSI = CSI - ZDI
        STR = DEXP(CSR)/TOL
        CSR = STR*DCOS(CSI)
        CSI = STR*DSIN(CSI)
        CALL ZUCHK(CSR, CSI, NW, ASCLE, TOL)
        IF (NW.NE.0) GO TO 25
        YR(I) = CSR
        YI(I) = CSI
        NZ = NZ - 1
        IF (IC.EQ.KK-1) GO TO 40
        IC = KK
        GO TO 30
   25   CONTINUE
        IF(ALAS.LT.HELIM) GO TO 30
        ZDR = ZDR - ELIM
        S1R = S1R*CELMR
        S1I = S1I*CELMR
        S2R = S2R*CELMR
        S2I = S2I*CELMR
   30 CONTINUE
      NZ = N
      IF(IC.EQ.N) NZ=N-1
      GO TO 45
   40 CONTINUE
      NZ = KK - 2
   45 CONTINUE
      DO 50 I=1,NZ
        YR(I) = ZEROR
        YI(I) = ZEROI
   50 CONTINUE
      RETURN
      END
      SUBROUTINE ZSHCH(ZR, ZI, CSHR, CSHI, CCHR, CCHI)
C***BEGIN PROLOGUE  ZSHCH
C***REFER TO  ZBESK,ZBESH
C
C     ZSHCH COMPUTES THE COMPLEX HYPERBOLIC FUNCTIONS CSH=SINH(X+I*Y)
C     AND CCH=COSH(X+I*Y), WHERE I**2=-1.
C
C***ROUTINES CALLED  (NONE)
C***END PROLOGUE  ZSHCH
C
      DOUBLE PRECISION CCHI, CCHR, CH, CN, CSHI, CSHR, SH, SN, ZI, ZR,
     * DCOSH, DSINH
      SH = DSINH(ZR)
      CH = DCOSH(ZR)
      SN = DSIN(ZI)
      CN = DCOS(ZI)
      CSHR = SH*CN
      CSHI = CH*SN
      CCHR = CH*CN
      CCHI = SH*SN
      RETURN
      END
      SUBROUTINE ZRATI(ZR, ZI, FNU, N, CYR, CYI, TOL)
C***BEGIN PROLOGUE  ZRATI
C***REFER TO  ZBESI,ZBESK,ZBESH
C
C     ZRATI COMPUTES RATIOS OF I BESSEL FUNCTIONS BY BACKWARD
C     RECURRENCE.  THE STARTING INDEX IS DETERMINED BY FORWARD
C     RECURRENCE AS DESCRIBED IN J. RES. OF NAT. BUR. OF STANDARDS-B,
C     MATHEMATICAL SCIENCES, VOL 77B, P111-114, SEPTEMBER, 1973,
C     BESSEL FUNCTIONS I AND J OF COMPLEX ARGUMENT AND INTEGER ORDER,
C     BY D. J. SOOKNE.
C
C***ROUTINES CALLED  ZABS,ZDIV
C***END PROLOGUE  ZRATI
C     COMPLEX Z,CY(1),CONE,CZERO,P1,P2,T1,RZ,PT,CDFNU
      EXTERNAL ZABS
      DOUBLE PRECISION AK, AMAGZ, AP1, AP2, ARG, AZ, CDFNUI, CDFNUR,
     * CONEI, CONER, CYI, CYR, CZEROI, CZEROR, DFNU, FDNU, FLAM, FNU,
     * FNUP, PTI, PTR, P1I, P1R, P2I, P2R, RAK, RAP1, RHO, RT2, RZI,
     * RZR, TEST, TEST1, TOL, TTI, TTR, T1I, T1R, ZI, ZR, ZABS
      INTEGER I, ID, IDNU, INU, ITIME, K, KK, MAGZ, N
      DIMENSION CYR(N), CYI(N)
      DATA CZEROR,CZEROI,CONER,CONEI,RT2/
     1 0.0D0, 0.0D0, 1.0D0, 0.0D0, 1.41421356237309505D0 /
      AZ = ZABS(ZR,ZI)
      INU = INT(SNGL(FNU))
      IDNU = INU + N - 1
      MAGZ = INT(SNGL(AZ))
      AMAGZ = DBLE(FLOAT(MAGZ+1))
      FDNU = DBLE(FLOAT(IDNU))
      FNUP = DMAX1(AMAGZ,FDNU)
      ID = IDNU - MAGZ - 1
      ITIME = 1
      K = 1
      PTR = 1.0D0/AZ
      RZR = PTR*(ZR+ZR)*PTR
      RZI = -PTR*(ZI+ZI)*PTR
      T1R = RZR*FNUP
      T1I = RZI*FNUP
      P2R = -T1R
      P2I = -T1I
      P1R = CONER
      P1I = CONEI
      T1R = T1R + RZR
      T1I = T1I + RZI
      IF (ID.GT.0) ID = 0
      AP2 = ZABS(P2R,P2I)
      AP1 = ZABS(P1R,P1I)
C-----------------------------------------------------------------------
C     THE OVERFLOW TEST ON K(FNU+I-1,Z) BEFORE THE CALL TO CBKNU
C     GUARANTEES THAT P2 IS ON SCALE. SCALE TEST1 AND ALL SUBSEQUENT
C     P2 VALUES BY AP1 TO ENSURE THAT AN OVERFLOW DOES NOT OCCUR
C     PREMATURELY.
C-----------------------------------------------------------------------
      ARG = (AP2+AP2)/(AP1*TOL)
      TEST1 = DSQRT(ARG)
      TEST = TEST1
      RAP1 = 1.0D0/AP1
      P1R = P1R*RAP1
      P1I = P1I*RAP1
      P2R = P2R*RAP1
      P2I = P2I*RAP1
      AP2 = AP2*RAP1
   10 CONTINUE
      K = K + 1
      AP1 = AP2
      PTR = P2R
      PTI = P2I
      P2R = P1R - (T1R*PTR-T1I*PTI)
      P2I = P1I - (T1R*PTI+T1I*PTR)
      P1R = PTR
      P1I = PTI
      T1R = T1R + RZR
      T1I = T1I + RZI
      AP2 = ZABS(P2R,P2I)
      IF (AP1.LE.TEST) GO TO 10
      IF (ITIME.EQ.2) GO TO 20
      AK = ZABS(T1R,T1I)*0.5D0
      FLAM = AK + DSQRT(AK*AK-1.0D0)
      RHO = DMIN1(AP2/AP1,FLAM)
      TEST = TEST1*DSQRT(RHO/(RHO*RHO-1.0D0))
      ITIME = 2
      GO TO 10
   20 CONTINUE
      KK = K + 1 - ID
      AK = DBLE(FLOAT(KK))
      T1R = AK
      T1I = CZEROI
      DFNU = FNU + DBLE(FLOAT(N-1))
      P1R = 1.0D0/AP2
      P1I = CZEROI
      P2R = CZEROR
      P2I = CZEROI
      DO 30 I=1,KK
        PTR = P1R
        PTI = P1I
        RAP1 = DFNU + T1R
        TTR = RZR*RAP1
        TTI = RZI*RAP1
        P1R = (PTR*TTR-PTI*TTI) + P2R
        P1I = (PTR*TTI+PTI*TTR) + P2I
        P2R = PTR
        P2I = PTI
        T1R = T1R - CONER
   30 CONTINUE
      IF (P1R.NE.CZEROR .OR. P1I.NE.CZEROI) GO TO 40
      P1R = TOL
      P1I = TOL
   40 CONTINUE
      CALL ZDIV(P2R, P2I, P1R, P1I, CYR(N), CYI(N))
      IF (N.EQ.1) RETURN
      K = N - 1
      AK = DBLE(FLOAT(K))
      T1R = AK
      T1I = CZEROI
      CDFNUR = FNU*RZR
      CDFNUI = FNU*RZI
      DO 60 I=2,N
        PTR = CDFNUR + (T1R*RZR-T1I*RZI) + CYR(K+1)
        PTI = CDFNUI + (T1R*RZI+T1I*RZR) + CYI(K+1)
        AK = ZABS(PTR,PTI)
        IF (AK.NE.CZEROR) GO TO 50
        PTR = TOL
        PTI = TOL
        AK = TOL*RT2
   50   CONTINUE
        RAK = CONER/AK
        CYR(K) = RAK*PTR*RAK
        CYI(K) = -RAK*PTI*RAK
        T1R = T1R - CONER
        K = K - 1
   60 CONTINUE
      RETURN
      END
      SUBROUTINE ZS1S2(ZRR, ZRI, S1R, S1I, S2R, S2I, NZ, ASCLE, ALIM,
     * IUF)
C***BEGIN PROLOGUE  ZS1S2
C***REFER TO  ZBESK,ZAIRY
C
C     ZS1S2 TESTS FOR A POSSIBLE UNDERFLOW RESULTING FROM THE
C     ADDITION OF THE I AND K FUNCTIONS IN THE ANALYTIC CON-
C     TINUATION FORMULA WHERE S1=K FUNCTION AND S2=I FUNCTION.
C     ON KODE=1 THE I AND K FUNCTIONS ARE DIFFERENT ORDERS OF
C     MAGNITUDE, BUT FOR KODE=2 THEY CAN BE OF THE SAME ORDER
C     OF MAGNITUDE AND THE MAXIMUM MUST BE AT LEAST ONE
C     PRECISION ABOVE THE UNDERFLOW LIMIT.
C
C***ROUTINES CALLED  ZABS,ZEXP,ZLOG
C***END PROLOGUE  ZS1S2
C     COMPLEX CZERO,C1,S1,S1D,S2,ZR
      EXTERNAL ZABS
      DOUBLE PRECISION AA, ALIM, ALN, ASCLE, AS1, AS2, C1I, C1R, S1DI,
     * S1DR, S1I, S1R, S2I, S2R, ZEROI, ZEROR, ZRI, ZRR, ZABS
      INTEGER IUF, IDUM, NZ
      DATA ZEROR,ZEROI  / 0.0D0 , 0.0D0 /
      NZ = 0
      AS1 = ZABS(S1R,S1I)
      AS2 = ZABS(S2R,S2I)
      IF (S1R.EQ.0.0D0 .AND. S1I.EQ.0.0D0) GO TO 10
      IF (AS1.EQ.0.0D0) GO TO 10
      ALN = -ZRR - ZRR + DLOG(AS1)
      S1DR = S1R
      S1DI = S1I
      S1R = ZEROR
      S1I = ZEROI
      AS1 = ZEROR
      IF (ALN.LT.(-ALIM)) GO TO 10
      CALL ZLOG(S1DR, S1DI, C1R, C1I, IDUM)
      C1R = C1R - ZRR - ZRR
      C1I = C1I - ZRI - ZRI
      CALL ZEXP(C1R, C1I, S1R, S1I)
      AS1 = ZABS(S1R,S1I)
      IUF = IUF + 1
   10 CONTINUE
      AA = DMAX1(AS1,AS2)
      IF (AA.GT.ASCLE) RETURN
      S1R = ZEROR
      S1I = ZEROI
      S2R = ZEROR
      S2I = ZEROI
      NZ = 1
      IUF = 0
      RETURN
      END
      SUBROUTINE ZBUNK(ZR, ZI, FNU, KODE, MR, N, YR, YI, NZ, TOL, ELIM,
     * ALIM)
C***BEGIN PROLOGUE  ZBUNK
C***REFER TO  ZBESK,ZBESH
C
C     ZBUNK COMPUTES THE K BESSEL FUNCTION FOR FNU.GT.FNUL.
C     ACCORDING TO THE UNIFORM ASYMPTOTIC EXPANSION FOR K(FNU,Z)
C     IN ZUNK1 AND THE EXPANSION FOR H(2,FNU,Z) IN ZUNK2
C
C***ROUTINES CALLED  ZUNK1,ZUNK2
C***END PROLOGUE  ZBUNK
C     COMPLEX Y,Z
      DOUBLE PRECISION ALIM, AX, AY, ELIM, FNU, TOL, YI, YR, ZI, ZR
      INTEGER KODE, MR, N, NZ
      DIMENSION YR(N), YI(N)
      NZ = 0
      AX = DABS(ZR)*1.7321D0
      AY = DABS(ZI)
      IF (AY.GT.AX) GO TO 10
C-----------------------------------------------------------------------
C     ASYMPTOTIC EXPANSION FOR K(FNU,Z) FOR LARGE FNU APPLIED IN
C     -PI/3.LE.ARG(Z).LE.PI/3
C-----------------------------------------------------------------------
      CALL ZUNK1(ZR, ZI, FNU, KODE, MR, N, YR, YI, NZ, TOL, ELIM, ALIM)
      GO TO 20
   10 CONTINUE
C-----------------------------------------------------------------------
C     ASYMPTOTIC EXPANSION FOR H(2,FNU,Z*EXP(M*HPI)) FOR LARGE FNU
C     APPLIED IN PI/3.LT.ABS(ARG(Z)).LE.PI/2 WHERE M=+I OR -I
C     AND HPI=PI/2
C-----------------------------------------------------------------------
      CALL ZUNK2(ZR, ZI, FNU, KODE, MR, N, YR, YI, NZ, TOL, ELIM, ALIM)
   20 CONTINUE
      RETURN
      END
      SUBROUTINE ZMLRI(ZR, ZI, FNU, KODE, N, YR, YI, NZ, TOL)
C***BEGIN PROLOGUE  ZMLRI
C***REFER TO  ZBESI,ZBESK
C
C     ZMLRI COMPUTES THE I BESSEL FUNCTION FOR RE(Z).GE.0.0 BY THE
C     MILLER ALGORITHM NORMALIZED BY A NEUMANN SERIES.
C
C***ROUTINES CALLED  DGAMLN,D1MACH,ZABS,ZEXP,ZLOG,ZMLT
C***END PROLOGUE  ZMLRI
C     COMPLEX CK,CNORM,CONE,CTWO,CZERO,PT,P1,P2,RZ,SUM,Y,Z
      EXTERNAL ZABS
      DOUBLE PRECISION ACK, AK, AP, AT, AZ, BK, CKI, CKR, CNORMI,
     * CNORMR, CONEI, CONER, FKAP, FKK, FLAM, FNF, FNU, PTI, PTR, P1I,
     * P1R, P2I, P2R, RAZ, RHO, RHO2, RZI, RZR, SCLE, STI, STR, SUMI,
     * SUMR, TFNF, TOL, TST, YI, YR, ZEROI, ZEROR, ZI, ZR, DGAMLN,
     * D1MACH, ZABS
      INTEGER I, IAZ, IDUM, IFNU, INU, ITIME, K, KK, KM, KODE, M, N, NZ
      DIMENSION YR(N), YI(N)
      DATA ZEROR,ZEROI,CONER,CONEI / 0.0D0, 0.0D0, 1.0D0, 0.0D0 /
      SCLE = D1MACH(1)/TOL
      NZ=0
      AZ = ZABS(ZR,ZI)
      IAZ = INT(SNGL(AZ))
      IFNU = INT(SNGL(FNU))
      INU = IFNU + N - 1
      AT = DBLE(FLOAT(IAZ)) + 1.0D0
      RAZ = 1.0D0/AZ
      STR = ZR*RAZ
      STI = -ZI*RAZ
      CKR = STR*AT*RAZ
      CKI = STI*AT*RAZ
      RZR = (STR+STR)*RAZ
      RZI = (STI+STI)*RAZ
      P1R = ZEROR
      P1I = ZEROI
      P2R = CONER
      P2I = CONEI
      ACK = (AT+1.0D0)*RAZ
      RHO = ACK + DSQRT(ACK*ACK-1.0D0)
      RHO2 = RHO*RHO
      TST = (RHO2+RHO2)/((RHO2-1.0D0)*(RHO-1.0D0))
      TST = TST/TOL
C-----------------------------------------------------------------------
C     COMPUTE RELATIVE TRUNCATION ERROR INDEX FOR SERIES
C-----------------------------------------------------------------------
      AK = AT
      DO 10 I=1,80
        PTR = P2R
        PTI = P2I
        P2R = P1R - (CKR*PTR-CKI*PTI)
        P2I = P1I - (CKI*PTR+CKR*PTI)
        P1R = PTR
        P1I = PTI
        CKR = CKR + RZR
        CKI = CKI + RZI
        AP = ZABS(P2R,P2I)
        IF (AP.GT.TST*AK*AK) GO TO 20
        AK = AK + 1.0D0
   10 CONTINUE
      GO TO 110
   20 CONTINUE
      I = I + 1
      K = 0
      IF (INU.LT.IAZ) GO TO 40
C-----------------------------------------------------------------------
C     COMPUTE RELATIVE TRUNCATION ERROR FOR RATIOS
C-----------------------------------------------------------------------
      P1R = ZEROR
      P1I = ZEROI
      P2R = CONER
      P2I = CONEI
      AT = DBLE(FLOAT(INU)) + 1.0D0
      STR = ZR*RAZ
      STI = -ZI*RAZ
      CKR = STR*AT*RAZ
      CKI = STI*AT*RAZ
      ACK = AT*RAZ
      TST = DSQRT(ACK/TOL)
      ITIME = 1
      DO 30 K=1,80
        PTR = P2R
        PTI = P2I
        P2R = P1R - (CKR*PTR-CKI*PTI)
        P2I = P1I - (CKR*PTI+CKI*PTR)
        P1R = PTR
        P1I = PTI
        CKR = CKR + RZR
        CKI = CKI + RZI
        AP = ZABS(P2R,P2I)
        IF (AP.LT.TST) GO TO 30
        IF (ITIME.EQ.2) GO TO 40
        ACK = ZABS(CKR,CKI)
        FLAM = ACK + DSQRT(ACK*ACK-1.0D0)
        FKAP = AP/ZABS(P1R,P1I)
        RHO = DMIN1(FLAM,FKAP)
        TST = TST*DSQRT(RHO/(RHO*RHO-1.0D0))
        ITIME = 2
   30 CONTINUE
      GO TO 110
   40 CONTINUE
C-----------------------------------------------------------------------
C     BACKWARD RECURRENCE AND SUM NORMALIZING RELATION
C-----------------------------------------------------------------------
      K = K + 1
      KK = MAX0(I+IAZ,K+INU)
      FKK = DBLE(FLOAT(KK))
      P1R = ZEROR
      P1I = ZEROI
C-----------------------------------------------------------------------
C     SCALE P2 AND SUM BY SCLE
C-----------------------------------------------------------------------
      P2R = SCLE
      P2I = ZEROI
      FNF = FNU - DBLE(FLOAT(IFNU))
      TFNF = FNF + FNF
      BK = DGAMLN(FKK+TFNF+1.0D0,IDUM) - DGAMLN(FKK+1.0D0,IDUM) -
     * DGAMLN(TFNF+1.0D0,IDUM)
      BK = DEXP(BK)
      SUMR = ZEROR
      SUMI = ZEROI
      KM = KK - INU
      DO 50 I=1,KM
        PTR = P2R
        PTI = P2I
        P2R = P1R + (FKK+FNF)*(RZR*PTR-RZI*PTI)
        P2I = P1I + (FKK+FNF)*(RZI*PTR+RZR*PTI)
        P1R = PTR
        P1I = PTI
        AK = 1.0D0 - TFNF/(FKK+TFNF)
        ACK = BK*AK
        SUMR = SUMR + (ACK+BK)*P1R
        SUMI = SUMI + (ACK+BK)*P1I
        BK = ACK
        FKK = FKK - 1.0D0
   50 CONTINUE
      YR(N) = P2R
      YI(N) = P2I
      IF (N.EQ.1) GO TO 70
      DO 60 I=2,N
        PTR = P2R
        PTI = P2I
        P2R = P1R + (FKK+FNF)*(RZR*PTR-RZI*PTI)
        P2I = P1I + (FKK+FNF)*(RZI*PTR+RZR*PTI)
        P1R = PTR
        P1I = PTI
        AK = 1.0D0 - TFNF/(FKK+TFNF)
        ACK = BK*AK
        SUMR = SUMR + (ACK+BK)*P1R
        SUMI = SUMI + (ACK+BK)*P1I
        BK = ACK
        FKK = FKK - 1.0D0
        M = N - I + 1
        YR(M) = P2R
        YI(M) = P2I
   60 CONTINUE
   70 CONTINUE
      IF (IFNU.LE.0) GO TO 90
      DO 80 I=1,IFNU
        PTR = P2R
        PTI = P2I
        P2R = P1R + (FKK+FNF)*(RZR*PTR-RZI*PTI)
        P2I = P1I + (FKK+FNF)*(RZR*PTI+RZI*PTR)
        P1R = PTR
        P1I = PTI
        AK = 1.0D0 - TFNF/(FKK+TFNF)
        ACK = BK*AK
        SUMR = SUMR + (ACK+BK)*P1R
        SUMI = SUMI + (ACK+BK)*P1I
        BK = ACK
        FKK = FKK - 1.0D0
   80 CONTINUE
   90 CONTINUE
      PTR = ZR
      PTI = ZI
      IF (KODE.EQ.2) PTR = ZEROR
      CALL ZLOG(RZR, RZI, STR, STI, IDUM)
      P1R = -FNF*STR + PTR
      P1I = -FNF*STI + PTI
      AP = DGAMLN(1.0D0+FNF,IDUM)
      PTR = P1R - AP
      PTI = P1I
C-----------------------------------------------------------------------
C     THE DIVISION CEXP(PT)/(SUM+P2) IS ALTERED TO AVOID OVERFLOW
C     IN THE DENOMINATOR BY SQUARING LARGE QUANTITIES
C-----------------------------------------------------------------------
      P2R = P2R + SUMR
      P2I = P2I + SUMI
      AP = ZABS(P2R,P2I)
      P1R = 1.0D0/AP
      CALL ZEXP(PTR, PTI, STR, STI)
      CKR = STR*P1R
      CKI = STI*P1R
      PTR = P2R*P1R
      PTI = -P2I*P1R
      CALL ZMLT(CKR, CKI, PTR, PTI, CNORMR, CNORMI)
      DO 100 I=1,N
        STR = YR(I)*CNORMR - YI(I)*CNORMI
        YI(I) = YR(I)*CNORMI + YI(I)*CNORMR
        YR(I) = STR
  100 CONTINUE
      RETURN
  110 CONTINUE
      NZ=-2
      RETURN
      END
      SUBROUTINE ZWRSK(ZRR, ZRI, FNU, KODE, N, YR, YI, NZ, CWR, CWI,
     * TOL, ELIM, ALIM)
C***BEGIN PROLOGUE  ZWRSK
C***REFER TO  ZBESI,ZBESK
C
C     ZWRSK COMPUTES THE I BESSEL FUNCTION FOR RE(Z).GE.0.0 BY
C     NORMALIZING THE I FUNCTION RATIOS FROM ZRATI BY THE WRONSKIAN
C
C***ROUTINES CALLED  D1MACH,ZBKNU,ZRATI,ZABS
C***END PROLOGUE  ZWRSK
C     COMPLEX CINU,CSCL,CT,CW,C1,C2,RCT,ST,Y,ZR
      EXTERNAL ZABS
      DOUBLE PRECISION ACT, ACW, ALIM, ASCLE, CINUI, CINUR, CSCLR, CTI,
     * CTR, CWI, CWR, C1I, C1R, C2I, C2R, ELIM, FNU, PTI, PTR, RACT,
     * STI, STR, TOL, YI, YR, ZRI, ZRR, ZABS, D1MACH
      INTEGER I, KODE, N, NW, NZ
      DIMENSION YR(N), YI(N), CWR(2), CWI(2)
C-----------------------------------------------------------------------
C     I(FNU+I-1,Z) BY BACKWARD RECURRENCE FOR RATIOS
C     Y(I)=I(FNU+I,Z)/I(FNU+I-1,Z) FROM CRATI NORMALIZED BY THE
C     WRONSKIAN WITH K(FNU,Z) AND K(FNU+1,Z) FROM CBKNU.
C-----------------------------------------------------------------------
      NZ = 0
      CALL ZBKNU(ZRR, ZRI, FNU, KODE, 2, CWR, CWI, NW, TOL, ELIM, ALIM)
      IF (NW.NE.0) GO TO 50
      CALL ZRATI(ZRR, ZRI, FNU, N, YR, YI, TOL)
C-----------------------------------------------------------------------
C     RECUR FORWARD ON I(FNU+1,Z) = R(FNU,Z)*I(FNU,Z),
C     R(FNU+J-1,Z)=Y(J),  J=1,...,N
C-----------------------------------------------------------------------
      CINUR = 1.0D0
      CINUI = 0.0D0
      IF (KODE.EQ.1) GO TO 10
      CINUR = DCOS(ZRI)
      CINUI = DSIN(ZRI)
   10 CONTINUE
C-----------------------------------------------------------------------
C     ON LOW EXPONENT MACHINES THE K FUNCTIONS CAN BE CLOSE TO BOTH
C     THE UNDER AND OVERFLOW LIMITS AND THE NORMALIZATION MUST BE
C     SCALED TO PREVENT OVER OR UNDERFLOW. CUOIK HAS DETERMINED THAT
C     THE RESULT IS ON SCALE.
C-----------------------------------------------------------------------
      ACW = ZABS(CWR(2),CWI(2))
      ASCLE = 1.0D+3*D1MACH(1)/TOL
      CSCLR = 1.0D0
      IF (ACW.GT.ASCLE) GO TO 20
      CSCLR = 1.0D0/TOL
      GO TO 30
   20 CONTINUE
      ASCLE = 1.0D0/ASCLE
      IF (ACW.LT.ASCLE) GO TO 30
      CSCLR = TOL
   30 CONTINUE
      C1R = CWR(1)*CSCLR
      C1I = CWI(1)*CSCLR
      C2R = CWR(2)*CSCLR
      C2I = CWI(2)*CSCLR
      STR = YR(1)
      STI = YI(1)
C-----------------------------------------------------------------------
C     CINU=CINU*(CONJG(CT)/CABS(CT))*(1.0D0/CABS(CT) PREVENTS
C     UNDER- OR OVERFLOW PREMATURELY BY SQUARING CABS(CT)
C-----------------------------------------------------------------------
      PTR = STR*C1R - STI*C1I
      PTI = STR*C1I + STI*C1R
      PTR = PTR + C2R
      PTI = PTI + C2I
      CTR = ZRR*PTR - ZRI*PTI
      CTI = ZRR*PTI + ZRI*PTR
      ACT = ZABS(CTR,CTI)
      RACT = 1.0D0/ACT
      CTR = CTR*RACT
      CTI = -CTI*RACT
      PTR = CINUR*RACT
      PTI = CINUI*RACT
      CINUR = PTR*CTR - PTI*CTI
      CINUI = PTR*CTI + PTI*CTR
      YR(1) = CINUR*CSCLR
      YI(1) = CINUI*CSCLR
      IF (N.EQ.1) RETURN
      DO 40 I=2,N
        PTR = STR*CINUR - STI*CINUI
        CINUI = STR*CINUI + STI*CINUR
        CINUR = PTR
        STR = YR(I)
        STI = YI(I)
        YR(I) = CINUR*CSCLR
        YI(I) = CINUI*CSCLR
   40 CONTINUE
      RETURN
   50 CONTINUE
      NZ = -1
      IF(NW.EQ.(-2)) NZ=-2
      RETURN
      END
      SUBROUTINE ZSERI(ZR, ZI, FNU, KODE, N, YR, YI, NZ, TOL, ELIM,
     * ALIM)
C***BEGIN PROLOGUE  ZSERI
C***REFER TO  ZBESI,ZBESK
C
C     ZSERI COMPUTES THE I BESSEL FUNCTION FOR REAL(Z).GE.0.0 BY
C     MEANS OF THE POWER SERIES FOR LARGE CABS(Z) IN THE
C     REGION CABS(Z).LE.2*SQRT(FNU+1). NZ=0 IS A NORMAL RETURN.
C     NZ.GT.0 MEANS THAT THE LAST NZ COMPONENTS WERE SET TO ZERO
C     DUE TO UNDERFLOW. NZ.LT.0 MEANS UNDERFLOW OCCURRED, BUT THE
C     CONDITION CABS(Z).LE.2*SQRT(FNU+1) WAS VIOLATED AND THE
C     COMPUTATION MUST BE COMPLETED IN ANOTHER ROUTINE WITH N=N-ABS(NZ).
C
C***ROUTINES CALLED  DGAMLN,D1MACH,ZUCHK,ZABS,ZDIV,ZLOG,ZMLT
C***END PROLOGUE  ZSERI
C     COMPLEX AK1,CK,COEF,CONE,CRSC,CSCL,CZ,CZERO,HZ,RZ,S1,S2,Y,Z
      EXTERNAL ZABS
      DOUBLE PRECISION AA, ACZ, AK, AK1I, AK1R, ALIM, ARM, ASCLE, ATOL,
     * AZ, CKI, CKR, COEFI, COEFR, CONEI, CONER, CRSCR, CZI, CZR, DFNU,
     * ELIM, FNU, FNUP, HZI, HZR, RAZ, RS, RTR1, RZI, RZR, S, SS, STI,
     * STR, S1I, S1R, S2I, S2R, TOL, YI, YR, WI, WR, ZEROI, ZEROR, ZI,
     * ZR, DGAMLN, D1MACH, ZABS
      INTEGER I, IB, IDUM, IFLAG, IL, K, KODE, L, M, N, NN, NZ, NW
      DIMENSION YR(N), YI(N), WR(2), WI(2)
      DATA ZEROR,ZEROI,CONER,CONEI / 0.0D0, 0.0D0, 1.0D0, 0.0D0 /
C
      NZ = 0
      AZ = ZABS(ZR,ZI)
      IF (AZ.EQ.0.0D0) GO TO 160
      ARM = 1.0D+3*D1MACH(1)
      RTR1 = DSQRT(ARM)
      CRSCR = 1.0D0
      IFLAG = 0
      IF (AZ.LT.ARM) GO TO 150
      HZR = 0.5D0*ZR
      HZI = 0.5D0*ZI
      CZR = ZEROR
      CZI = ZEROI
      IF (AZ.LE.RTR1) GO TO 10
      CALL ZMLT(HZR, HZI, HZR, HZI, CZR, CZI)
   10 CONTINUE
      ACZ = ZABS(CZR,CZI)
      NN = N
      CALL ZLOG(HZR, HZI, CKR, CKI, IDUM)
   20 CONTINUE
      DFNU = FNU + DBLE(FLOAT(NN-1))
      FNUP = DFNU + 1.0D0
C-----------------------------------------------------------------------
C     UNDERFLOW TEST
C-----------------------------------------------------------------------
      AK1R = CKR*DFNU
      AK1I = CKI*DFNU
      AK = DGAMLN(FNUP,IDUM)
      AK1R = AK1R - AK
      IF (KODE.EQ.2) AK1R = AK1R - ZR
      IF (AK1R.GT.(-ELIM)) GO TO 40
   30 CONTINUE
      NZ = NZ + 1
      YR(NN) = ZEROR
      YI(NN) = ZEROI
      IF (ACZ.GT.DFNU) GO TO 190
      NN = NN - 1
      IF (NN.EQ.0) RETURN
      GO TO 20
   40 CONTINUE
      IF (AK1R.GT.(-ALIM)) GO TO 50
      IFLAG = 1
      SS = 1.0D0/TOL
      CRSCR = TOL
      ASCLE = ARM*SS
   50 CONTINUE
      AA = DEXP(AK1R)
      IF (IFLAG.EQ.1) AA = AA*SS
      COEFR = AA*DCOS(AK1I)
      COEFI = AA*DSIN(AK1I)
      ATOL = TOL*ACZ/FNUP
      IL = MIN0(2,NN)
      DO 90 I=1,IL
        DFNU = FNU + DBLE(FLOAT(NN-I))
        FNUP = DFNU + 1.0D0
        S1R = CONER
        S1I = CONEI
        IF (ACZ.LT.TOL*FNUP) GO TO 70
        AK1R = CONER
        AK1I = CONEI
        AK = FNUP + 2.0D0
        S = FNUP
        AA = 2.0D0
   60   CONTINUE
        RS = 1.0D0/S
        STR = AK1R*CZR - AK1I*CZI
        STI = AK1R*CZI + AK1I*CZR
        AK1R = STR*RS
        AK1I = STI*RS
        S1R = S1R + AK1R
        S1I = S1I + AK1I
        S = S + AK
        AK = AK + 2.0D0
        AA = AA*ACZ*RS
        IF (AA.GT.ATOL) GO TO 60
   70   CONTINUE
        S2R = S1R*COEFR - S1I*COEFI
        S2I = S1R*COEFI + S1I*COEFR
        WR(I) = S2R
        WI(I) = S2I
        IF (IFLAG.EQ.0) GO TO 80
        CALL ZUCHK(S2R, S2I, NW, ASCLE, TOL)
        IF (NW.NE.0) GO TO 30
   80   CONTINUE
        M = NN - I + 1
        YR(M) = S2R*CRSCR
        YI(M) = S2I*CRSCR
        IF (I.EQ.IL) GO TO 90
        CALL ZDIV(COEFR, COEFI, HZR, HZI, STR, STI)
        COEFR = STR*DFNU
        COEFI = STI*DFNU
   90 CONTINUE
      IF (NN.LE.2) RETURN
      K = NN - 2
      AK = DBLE(FLOAT(K))
      RAZ = 1.0D0/AZ
      STR = ZR*RAZ
      STI = -ZI*RAZ
      RZR = (STR+STR)*RAZ
      RZI = (STI+STI)*RAZ
      IF (IFLAG.EQ.1) GO TO 120
      IB = 3
  100 CONTINUE
      DO 110 I=IB,NN
        YR(K) = (AK+FNU)*(RZR*YR(K+1)-RZI*YI(K+1)) + YR(K+2)
        YI(K) = (AK+FNU)*(RZR*YI(K+1)+RZI*YR(K+1)) + YI(K+2)
        AK = AK - 1.0D0
        K = K - 1
  110 CONTINUE
      RETURN
C-----------------------------------------------------------------------
C     RECUR BACKWARD WITH SCALED VALUES
C-----------------------------------------------------------------------
  120 CONTINUE
C-----------------------------------------------------------------------
C     EXP(-ALIM)=EXP(-ELIM)/TOL=APPROX. ONE PRECISION ABOVE THE
C     UNDERFLOW LIMIT = ASCLE = D1MACH(1)*SS*1.0D+3
C-----------------------------------------------------------------------
      S1R = WR(1)
      S1I = WI(1)
      S2R = WR(2)
      S2I = WI(2)
      DO 130 L=3,NN
        CKR = S2R
        CKI = S2I
        S2R = S1R + (AK+FNU)*(RZR*CKR-RZI*CKI)
        S2I = S1I + (AK+FNU)*(RZR*CKI+RZI*CKR)
        S1R = CKR
        S1I = CKI
        CKR = S2R*CRSCR
        CKI = S2I*CRSCR
        YR(K) = CKR
        YI(K) = CKI
        AK = AK - 1.0D0
        K = K - 1
        IF (ZABS(CKR,CKI).GT.ASCLE) GO TO 140
  130 CONTINUE
      RETURN
  140 CONTINUE
      IB = L + 1
      IF (IB.GT.NN) RETURN
      GO TO 100
  150 CONTINUE
      NZ = N
      IF (FNU.EQ.0.0D0) NZ = NZ - 1
  160 CONTINUE
      YR(1) = ZEROR
      YI(1) = ZEROI
      IF (FNU.NE.0.0D0) GO TO 170
      YR(1) = CONER
      YI(1) = CONEI
  170 CONTINUE
      IF (N.EQ.1) RETURN
      DO 180 I=2,N
        YR(I) = ZEROR
        YI(I) = ZEROI
  180 CONTINUE
      RETURN
C-----------------------------------------------------------------------
C     RETURN WITH NZ.LT.0 IF CABS(Z*Z/4).GT.FNU+N-NZ-1 COMPLETE
C     THE CALCULATION IN CBINU WITH N=N-IABS(NZ)
C-----------------------------------------------------------------------
  190 CONTINUE
      NZ = -NZ
      RETURN
      END
      SUBROUTINE ZASYI(ZR, ZI, FNU, KODE, N, YR, YI, NZ, RL, TOL, ELIM,
     * ALIM)
C***BEGIN PROLOGUE  ZASYI
C***REFER TO  ZBESI,ZBESK
C
C     ZASYI COMPUTES THE I BESSEL FUNCTION FOR REAL(Z).GE.0.0 BY
C     MEANS OF THE ASYMPTOTIC EXPANSION FOR LARGE CABS(Z) IN THE
C     REGION CABS(Z).GT.MAX(RL,FNU*FNU/2). NZ=0 IS A NORMAL RETURN.
C     NZ.LT.0 INDICATES AN OVERFLOW ON KODE=1.
C
C***ROUTINES CALLED  D1MACH,ZABS,ZDIV,ZEXP,ZMLT,ZSQRT
C***END PROLOGUE  ZASYI
C     COMPLEX AK1,CK,CONE,CS1,CS2,CZ,CZERO,DK,EZ,P1,RZ,S2,Y,Z
      EXTERNAL ZABS
      DOUBLE PRECISION AA, AEZ, AK, AK1I, AK1R, ALIM, ARG, ARM, ATOL,
     * AZ, BB, BK, CKI, CKR, CONEI, CONER, CS1I, CS1R, CS2I, CS2R, CZI,
     * CZR, DFNU, DKI, DKR, DNU2, ELIM, EZI, EZR, FDN, FNU, PI, P1I,
     * P1R, RAZ, RL, RTPI, RTR1, RZI, RZR, S, SGN, SQK, STI, STR, S2I,
     * S2R, TOL, TZI, TZR, YI, YR, ZEROI, ZEROR, ZI, ZR, D1MACH, ZABS
      INTEGER I, IB, IL, INU, J, JL, K, KODE, KODED, M, N, NN, NZ
      DIMENSION YR(N), YI(N)
      DATA PI, RTPI  /3.14159265358979324D0 , 0.159154943091895336D0 /
      DATA ZEROR,ZEROI,CONER,CONEI / 0.0D0, 0.0D0, 1.0D0, 0.0D0 /
C
      NZ = 0
      AZ = ZABS(ZR,ZI)
      ARM = 1.0D+3*D1MACH(1)
      RTR1 = DSQRT(ARM)
      IL = MIN0(2,N)
      DFNU = FNU + DBLE(FLOAT(N-IL))
C-----------------------------------------------------------------------
C     OVERFLOW TEST
C-----------------------------------------------------------------------
      RAZ = 1.0D0/AZ
      STR = ZR*RAZ
      STI = -ZI*RAZ
      AK1R = RTPI*STR*RAZ
      AK1I = RTPI*STI*RAZ
      CALL ZSQRT(AK1R, AK1I, AK1R, AK1I)
      CZR = ZR
      CZI = ZI
      IF (KODE.NE.2) GO TO 10
      CZR = ZEROR
      CZI = ZI
   10 CONTINUE
      IF (DABS(CZR).GT.ELIM) GO TO 100
      DNU2 = DFNU + DFNU
      KODED = 1
      IF ((DABS(CZR).GT.ALIM) .AND. (N.GT.2)) GO TO 20
      KODED = 0
      CALL ZEXP(CZR, CZI, STR, STI)
      CALL ZMLT(AK1R, AK1I, STR, STI, AK1R, AK1I)
   20 CONTINUE
      FDN = 0.0D0
      IF (DNU2.GT.RTR1) FDN = DNU2*DNU2
      EZR = ZR*8.0D0
      EZI = ZI*8.0D0
C-----------------------------------------------------------------------
C     WHEN Z IS IMAGINARY, THE ERROR TEST MUST BE MADE RELATIVE TO THE
C     FIRST RECIPROCAL POWER SINCE THIS IS THE LEADING TERM OF THE
C     EXPANSION FOR THE IMAGINARY PART.
C-----------------------------------------------------------------------
      AEZ = 8.0D0*AZ
      S = TOL/AEZ
      JL = INT(SNGL(RL+RL)) + 2
      P1R = ZEROR
      P1I = ZEROI
      IF (ZI.EQ.0.0D0) GO TO 30
C-----------------------------------------------------------------------
C     CALCULATE EXP(PI*(0.5+FNU+N-IL)*I) TO MINIMIZE LOSSES OF
C     SIGNIFICANCE WHEN FNU OR N IS LARGE
C-----------------------------------------------------------------------
      INU = INT(SNGL(FNU))
      ARG = (FNU-DBLE(FLOAT(INU)))*PI
      INU = INU + N - IL
      AK = -DSIN(ARG)
      BK = DCOS(ARG)
      IF (ZI.LT.0.0D0) BK = -BK
      P1R = AK
      P1I = BK
      IF (MOD(INU,2).EQ.0) GO TO 30
      P1R = -P1R
      P1I = -P1I
   30 CONTINUE
      DO 70 K=1,IL
        SQK = FDN - 1.0D0
        ATOL = S*DABS(SQK)
        SGN = 1.0D0
        CS1R = CONER
        CS1I = CONEI
        CS2R = CONER
        CS2I = CONEI
        CKR = CONER
        CKI = CONEI
        AK = 0.0D0
        AA = 1.0D0
        BB = AEZ
        DKR = EZR
        DKI = EZI
        DO 40 J=1,JL
          CALL ZDIV(CKR, CKI, DKR, DKI, STR, STI)
          CKR = STR*SQK
          CKI = STI*SQK
          CS2R = CS2R + CKR
          CS2I = CS2I + CKI
          SGN = -SGN
          CS1R = CS1R + CKR*SGN
          CS1I = CS1I + CKI*SGN
          DKR = DKR + EZR
          DKI = DKI + EZI
          AA = AA*DABS(SQK)/BB
          BB = BB + AEZ
          AK = AK + 8.0D0
          SQK = SQK - AK
          IF (AA.LE.ATOL) GO TO 50
   40   CONTINUE
        GO TO 110
   50   CONTINUE
        S2R = CS1R
        S2I = CS1I
        IF (ZR+ZR.GE.ELIM) GO TO 60
        TZR = ZR + ZR
        TZI = ZI + ZI
        CALL ZEXP(-TZR, -TZI, STR, STI)
        CALL ZMLT(STR, STI, P1R, P1I, STR, STI)
        CALL ZMLT(STR, STI, CS2R, CS2I, STR, STI)
        S2R = S2R + STR
        S2I = S2I + STI
   60   CONTINUE
        FDN = FDN + 8.0D0*DFNU + 4.0D0
        P1R = -P1R
        P1I = -P1I
        M = N - IL + K
        YR(M) = S2R*AK1R - S2I*AK1I
        YI(M) = S2R*AK1I + S2I*AK1R
   70 CONTINUE
      IF (N.LE.2) RETURN
      NN = N
      K = NN - 2
      AK = DBLE(FLOAT(K))
      STR = ZR*RAZ
      STI = -ZI*RAZ
      RZR = (STR+STR)*RAZ
      RZI = (STI+STI)*RAZ
      IB = 3
      DO 80 I=IB,NN
        YR(K) = (AK+FNU)*(RZR*YR(K+1)-RZI*YI(K+1)) + YR(K+2)
        YI(K) = (AK+FNU)*(RZR*YI(K+1)+RZI*YR(K+1)) + YI(K+2)
        AK = AK - 1.0D0
        K = K - 1
   80 CONTINUE
      IF (KODED.EQ.0) RETURN
      CALL ZEXP(CZR, CZI, CKR, CKI)
      DO 90 I=1,NN
        STR = YR(I)*CKR - YI(I)*CKI
        YI(I) = YR(I)*CKI + YI(I)*CKR
        YR(I) = STR
   90 CONTINUE
      RETURN
  100 CONTINUE
      NZ = -1
      RETURN
  110 CONTINUE
      NZ=-2
      RETURN
      END
      SUBROUTINE ZUOIK(ZR, ZI, FNU, KODE, IKFLG, N, YR, YI, NUF, TOL,
     * ELIM, ALIM)
C***BEGIN PROLOGUE  ZUOIK
C***REFER TO  ZBESI,ZBESK,ZBESH
C
C     ZUOIK COMPUTES THE LEADING TERMS OF THE UNIFORM ASYMPTOTIC
C     EXPANSIONS FOR THE I AND K FUNCTIONS AND COMPARES THEM
C     (IN LOGARITHMIC FORM) TO ALIM AND ELIM FOR OVER AND UNDERFLOW
C     WHERE ALIM.LT.ELIM. IF THE MAGNITUDE, BASED ON THE LEADING
C     EXPONENTIAL, IS LESS THAN ALIM OR GREATER THAN -ALIM, THEN
C     THE RESULT IS ON SCALE. IF NOT, THEN A REFINED TEST USING OTHER
C     MULTIPLIERS (IN LOGARITHMIC FORM) IS MADE BASED ON ELIM. HERE
C     EXP(-ELIM)=SMALLEST MACHINE NUMBER*1.0E+3 AND EXP(-ALIM)=
C     EXP(-ELIM)/TOL
C
C     IKFLG=1 MEANS THE I SEQUENCE IS TESTED
C          =2 MEANS THE K SEQUENCE IS TESTED
C     NUF = 0 MEANS THE LAST MEMBER OF THE SEQUENCE IS ON SCALE
C         =-1 MEANS AN OVERFLOW WOULD OCCUR
C     IKFLG=1 AND NUF.GT.0 MEANS THE LAST NUF Y VALUES WERE SET TO ZERO
C             THE FIRST N-NUF VALUES MUST BE SET BY ANOTHER ROUTINE
C     IKFLG=2 AND NUF.EQ.N MEANS ALL Y VALUES WERE SET TO ZERO
C     IKFLG=2 AND 0.LT.NUF.LT.N NOT CONSIDERED. Y MUST BE SET BY
C             ANOTHER ROUTINE
C
C***ROUTINES CALLED  ZUCHK,ZUNHJ,ZUNIK,D1MACH,ZABS,ZLOG
C***END PROLOGUE  ZUOIK
C     COMPLEX ARG,ASUM,BSUM,CWRK,CZ,CZERO,PHI,SUM,Y,Z,ZB,ZETA1,ZETA2,ZN,
C    *ZR
      EXTERNAL ZABS
      DOUBLE PRECISION AARG, AIC, ALIM, APHI, ARGI, ARGR, ASUMI, ASUMR,
     * ASCLE, AX, AY, BSUMI, BSUMR, CWRKI, CWRKR, CZI, CZR, ELIM, FNN,
     * FNU, GNN, GNU, PHII, PHIR, RCZ, STR, STI, SUMI, SUMR, TOL, YI,
     * YR, ZBI, ZBR, ZEROI, ZEROR, ZETA1I, ZETA1R, ZETA2I, ZETA2R, ZI,
     * ZNI, ZNR, ZR, ZRI, ZRR, D1MACH, ZABS
      INTEGER I, IDUM, IFORM, IKFLG, INIT, KODE, N, NN, NUF, NW
      DIMENSION YR(N), YI(N), CWRKR(16), CWRKI(16)
      DATA ZEROR,ZEROI / 0.0D0, 0.0D0 /
      DATA AIC / 1.265512123484645396D+00 /
      NUF = 0
      NN = N
      ZRR = ZR
      ZRI = ZI
      IF (ZR.GE.0.0D0) GO TO 10
      ZRR = -ZR
      ZRI = -ZI
   10 CONTINUE
      ZBR = ZRR
      ZBI = ZRI
      AX = DABS(ZR)*1.7321D0
      AY = DABS(ZI)
      IFORM = 1
      IF (AY.GT.AX) IFORM = 2
      GNU = DMAX1(FNU,1.0D0)
      IF (IKFLG.EQ.1) GO TO 20
      FNN = DBLE(FLOAT(NN))
      GNN = FNU + FNN - 1.0D0
      GNU = DMAX1(GNN,FNN)
   20 CONTINUE
C-----------------------------------------------------------------------
C     ONLY THE MAGNITUDE OF ARG AND PHI ARE NEEDED ALONG WITH THE
C     REAL PARTS OF ZETA1, ZETA2 AND ZB. NO ATTEMPT IS MADE TO GET
C     THE SIGN OF THE IMAGINARY PART CORRECT.
C-----------------------------------------------------------------------
      IF (IFORM.EQ.2) GO TO 30
      INIT = 0
      CALL ZUNIK(ZRR, ZRI, GNU, IKFLG, 1, TOL, INIT, PHIR, PHII,
     * ZETA1R, ZETA1I, ZETA2R, ZETA2I, SUMR, SUMI, CWRKR, CWRKI)
      CZR = -ZETA1R + ZETA2R
      CZI = -ZETA1I + ZETA2I
      GO TO 50
   30 CONTINUE
      ZNR = ZRI
      ZNI = -ZRR
      IF (ZI.GT.0.0D0) GO TO 40
      ZNR = -ZNR
   40 CONTINUE
      CALL ZUNHJ(ZNR, ZNI, GNU, 1, TOL, PHIR, PHII, ARGR, ARGI, ZETA1R,
     * ZETA1I, ZETA2R, ZETA2I, ASUMR, ASUMI, BSUMR, BSUMI)
      CZR = -ZETA1R + ZETA2R
      CZI = -ZETA1I + ZETA2I
      AARG = ZABS(ARGR,ARGI)
   50 CONTINUE
      IF (KODE.EQ.1) GO TO 60
      CZR = CZR - ZBR
      CZI = CZI - ZBI
   60 CONTINUE
      IF (IKFLG.EQ.1) GO TO 70
      CZR = -CZR
      CZI = -CZI
   70 CONTINUE
      APHI = ZABS(PHIR,PHII)
      RCZ = CZR
C-----------------------------------------------------------------------
C     OVERFLOW TEST
C-----------------------------------------------------------------------
      IF (RCZ.GT.ELIM) GO TO 210
      IF (RCZ.LT.ALIM) GO TO 80
      RCZ = RCZ + DLOG(APHI)
      IF (IFORM.EQ.2) RCZ = RCZ - 0.25D0*DLOG(AARG) - AIC
      IF (RCZ.GT.ELIM) GO TO 210
      GO TO 130
   80 CONTINUE
C-----------------------------------------------------------------------
C     UNDERFLOW TEST
C-----------------------------------------------------------------------
      IF (RCZ.LT.(-ELIM)) GO TO 90
      IF (RCZ.GT.(-ALIM)) GO TO 130
      RCZ = RCZ + DLOG(APHI)
      IF (IFORM.EQ.2) RCZ = RCZ - 0.25D0*DLOG(AARG) - AIC
      IF (RCZ.GT.(-ELIM)) GO TO 110
   90 CONTINUE
      DO 100 I=1,NN
        YR(I) = ZEROR
        YI(I) = ZEROI
  100 CONTINUE
      NUF = NN
      RETURN
  110 CONTINUE
      ASCLE = 1.0D+3*D1MACH(1)/TOL
      CALL ZLOG(PHIR, PHII, STR, STI, IDUM)
      CZR = CZR + STR
      CZI = CZI + STI
      IF (IFORM.EQ.1) GO TO 120
      CALL ZLOG(ARGR, ARGI, STR, STI, IDUM)
      CZR = CZR - 0.25D0*STR - AIC
      CZI = CZI - 0.25D0*STI
  120 CONTINUE
      AX = DEXP(RCZ)/TOL
      AY = CZI
      CZR = AX*DCOS(AY)
      CZI = AX*DSIN(AY)
      CALL ZUCHK(CZR, CZI, NW, ASCLE, TOL)
      IF (NW.NE.0) GO TO 90
  130 CONTINUE
      IF (IKFLG.EQ.2) RETURN
      IF (N.EQ.1) RETURN
C-----------------------------------------------------------------------
C     SET UNDERFLOWS ON I SEQUENCE
C-----------------------------------------------------------------------
  140 CONTINUE
      GNU = FNU + DBLE(FLOAT(NN-1))
      IF (IFORM.EQ.2) GO TO 150
      INIT = 0
      CALL ZUNIK(ZRR, ZRI, GNU, IKFLG, 1, TOL, INIT, PHIR, PHII,
     * ZETA1R, ZETA1I, ZETA2R, ZETA2I, SUMR, SUMI, CWRKR, CWRKI)
      CZR = -ZETA1R + ZETA2R
      CZI = -ZETA1I + ZETA2I
      GO TO 160
  150 CONTINUE
      CALL ZUNHJ(ZNR, ZNI, GNU, 1, TOL, PHIR, PHII, ARGR, ARGI, ZETA1R,
     * ZETA1I, ZETA2R, ZETA2I, ASUMR, ASUMI, BSUMR, BSUMI)
      CZR = -ZETA1R + ZETA2R
      CZI = -ZETA1I + ZETA2I
      AARG = ZABS(ARGR,ARGI)
  160 CONTINUE
      IF (KODE.EQ.1) GO TO 170
      CZR = CZR - ZBR
      CZI = CZI - ZBI
  170 CONTINUE
      APHI = ZABS(PHIR,PHII)
      RCZ = CZR
      IF (RCZ.LT.(-ELIM)) GO TO 180
      IF (RCZ.GT.(-ALIM)) RETURN
      RCZ = RCZ + DLOG(APHI)
      IF (IFORM.EQ.2) RCZ = RCZ - 0.25D0*DLOG(AARG) - AIC
      IF (RCZ.GT.(-ELIM)) GO TO 190
  180 CONTINUE
      YR(NN) = ZEROR
      YI(NN) = ZEROI
      NN = NN - 1
      NUF = NUF + 1
      IF (NN.EQ.0) RETURN
      GO TO 140
  190 CONTINUE
      ASCLE = 1.0D+3*D1MACH(1)/TOL
      CALL ZLOG(PHIR, PHII, STR, STI, IDUM)
      CZR = CZR + STR
      CZI = CZI + STI
      IF (IFORM.EQ.1) GO TO 200
      CALL ZLOG(ARGR, ARGI, STR, STI, IDUM)
      CZR = CZR - 0.25D0*STR - AIC
      CZI = CZI - 0.25D0*STI
  200 CONTINUE
      AX = DEXP(RCZ)/TOL
      AY = CZI
      CZR = AX*DCOS(AY)
      CZI = AX*DSIN(AY)
      CALL ZUCHK(CZR, CZI, NW, ASCLE, TOL)
      IF (NW.NE.0) GO TO 180
      RETURN
  210 CONTINUE
      NUF = -1
      RETURN
      END
      SUBROUTINE ZACON(ZR, ZI, FNU, KODE, MR, N, YR, YI, NZ, RL, FNUL,
     * TOL, ELIM, ALIM)
C***BEGIN PROLOGUE  ZACON
C***REFER TO  ZBESK,ZBESH
C
C     ZACON APPLIES THE ANALYTIC CONTINUATION FORMULA
C
C         K(FNU,ZN*EXP(MP))=K(FNU,ZN)*EXP(-MP*FNU) - MP*I(FNU,ZN)
C                 MP=PI*MR*CMPLX(0.0,1.0)
C
C     TO CONTINUE THE K FUNCTION FROM THE RIGHT HALF TO THE LEFT
C     HALF Z PLANE
C
C***ROUTINES CALLED  ZBINU,ZBKNU,ZS1S2,D1MACH,ZABS,ZMLT
C***END PROLOGUE  ZACON
C     COMPLEX CK,CONE,CSCL,CSCR,CSGN,CSPN,CY,CZERO,C1,C2,RZ,SC1,SC2,ST,
C    *S1,S2,Y,Z,ZN
      EXTERNAL ZABS
      DOUBLE PRECISION ALIM, ARG, ASCLE, AS2, AZN, BRY, BSCLE, CKI,
     * CKR, CONER, CPN, CSCL, CSCR, CSGNI, CSGNR, CSPNI, CSPNR,
     * CSR, CSRR, CSSR, CYI, CYR, C1I, C1M, C1R, C2I, C2R, ELIM, FMR,
     * FN, FNU, FNUL, PI, PTI, PTR, RAZN, RL, RZI, RZR, SC1I, SC1R,
     * SC2I, SC2R, SGN, SPN, STI, STR, S1I, S1R, S2I, S2R, TOL, YI, YR,
     * YY, ZEROR, ZI, ZNI, ZNR, ZR, D1MACH, ZABS
      INTEGER I, INU, IUF, KFLAG, KODE, MR, N, NN, NW, NZ
      DIMENSION YR(N), YI(N), CYR(2), CYI(2), CSSR(3), CSRR(3), BRY(3)
      DATA PI / 3.14159265358979324D0 /
      DATA ZEROR,CONER / 0.0D0,1.0D0 /
      NZ = 0
      ZNR = -ZR
      ZNI = -ZI
      NN = N
      CALL ZBINU(ZNR, ZNI, FNU, KODE, NN, YR, YI, NW, RL, FNUL, TOL,
     * ELIM, ALIM)
      IF (NW.LT.0) GO TO 90
C-----------------------------------------------------------------------
C     ANALYTIC CONTINUATION TO THE LEFT HALF PLANE FOR THE K FUNCTION
C-----------------------------------------------------------------------
      NN = MIN0(2,N)
      CALL ZBKNU(ZNR, ZNI, FNU, KODE, NN, CYR, CYI, NW, TOL, ELIM, ALIM)
      IF (NW.NE.0) GO TO 90
      S1R = CYR(1)
      S1I = CYI(1)
      FMR = DBLE(FLOAT(MR))
      SGN = -DSIGN(PI,FMR)
      CSGNR = ZEROR
      CSGNI = SGN
      IF (KODE.EQ.1) GO TO 10
      YY = -ZNI
      CPN = DCOS(YY)
      SPN = DSIN(YY)
      CALL ZMLT(CSGNR, CSGNI, CPN, SPN, CSGNR, CSGNI)
   10 CONTINUE
C-----------------------------------------------------------------------
C     CALCULATE CSPN=EXP(FNU*PI*I) TO MINIMIZE LOSSES OF SIGNIFICANCE
C     WHEN FNU IS LARGE
C-----------------------------------------------------------------------
      INU = INT(SNGL(FNU))
      ARG = (FNU-DBLE(FLOAT(INU)))*SGN
      CPN = DCOS(ARG)
      SPN = DSIN(ARG)
      CSPNR = CPN
      CSPNI = SPN
      IF (MOD(INU,2).EQ.0) GO TO 20
      CSPNR = -CSPNR
      CSPNI = -CSPNI
   20 CONTINUE
      IUF = 0
      C1R = S1R
      C1I = S1I
      C2R = YR(1)
      C2I = YI(1)
      ASCLE = 1.0D+3*D1MACH(1)/TOL
      IF (KODE.EQ.1) GO TO 30
      CALL ZS1S2(ZNR, ZNI, C1R, C1I, C2R, C2I, NW, ASCLE, ALIM, IUF)
      NZ = NZ + NW
      SC1R = C1R
      SC1I = C1I
   30 CONTINUE
      CALL ZMLT(CSPNR, CSPNI, C1R, C1I, STR, STI)
      CALL ZMLT(CSGNR, CSGNI, C2R, C2I, PTR, PTI)
      YR(1) = STR + PTR
      YI(1) = STI + PTI
      IF (N.EQ.1) RETURN
      CSPNR = -CSPNR
      CSPNI = -CSPNI
      S2R = CYR(2)
      S2I = CYI(2)
      C1R = S2R
      C1I = S2I
      C2R = YR(2)
      C2I = YI(2)
      IF (KODE.EQ.1) GO TO 40
      CALL ZS1S2(ZNR, ZNI, C1R, C1I, C2R, C2I, NW, ASCLE, ALIM, IUF)
      NZ = NZ + NW
      SC2R = C1R
      SC2I = C1I
   40 CONTINUE
      CALL ZMLT(CSPNR, CSPNI, C1R, C1I, STR, STI)
      CALL ZMLT(CSGNR, CSGNI, C2R, C2I, PTR, PTI)
      YR(2) = STR + PTR
      YI(2) = STI + PTI
      IF (N.EQ.2) RETURN
      CSPNR = -CSPNR
      CSPNI = -CSPNI
      AZN = ZABS(ZNR,ZNI)
      RAZN = 1.0D0/AZN
      STR = ZNR*RAZN
      STI = -ZNI*RAZN
      RZR = (STR+STR)*RAZN
      RZI = (STI+STI)*RAZN
      FN = FNU + 1.0D0
      CKR = FN*RZR
      CKI = FN*RZI
C-----------------------------------------------------------------------
C     SCALE NEAR EXPONENT EXTREMES DURING RECURRENCE ON K FUNCTIONS
C-----------------------------------------------------------------------
      CSCL = 1.0D0/TOL
      CSCR = TOL
      CSSR(1) = CSCL
      CSSR(2) = CONER
      CSSR(3) = CSCR
      CSRR(1) = CSCR
      CSRR(2) = CONER
      CSRR(3) = CSCL
      BRY(1) = ASCLE
      BRY(2) = 1.0D0/ASCLE
      BRY(3) = D1MACH(2)
      AS2 = ZABS(S2R,S2I)
      KFLAG = 2
      IF (AS2.GT.BRY(1)) GO TO 50
      KFLAG = 1
      GO TO 60
   50 CONTINUE
      IF (AS2.LT.BRY(2)) GO TO 60
      KFLAG = 3
   60 CONTINUE
      BSCLE = BRY(KFLAG)
      S1R = S1R*CSSR(KFLAG)
      S1I = S1I*CSSR(KFLAG)
      S2R = S2R*CSSR(KFLAG)
      S2I = S2I*CSSR(KFLAG)
      CSR = CSRR(KFLAG)
      DO 80 I=3,N
        STR = S2R
        STI = S2I
        S2R = CKR*STR - CKI*STI + S1R
        S2I = CKR*STI + CKI*STR + S1I
        S1R = STR
        S1I = STI
        C1R = S2R*CSR
        C1I = S2I*CSR
        STR = C1R
        STI = C1I
        C2R = YR(I)
        C2I = YI(I)
        IF (KODE.EQ.1) GO TO 70
        IF (IUF.LT.0) GO TO 70
        CALL ZS1S2(ZNR, ZNI, C1R, C1I, C2R, C2I, NW, ASCLE, ALIM, IUF)
        NZ = NZ + NW
        SC1R = SC2R
        SC1I = SC2I
        SC2R = C1R
        SC2I = C1I
        IF (IUF.NE.3) GO TO 70
        IUF = -4
        S1R = SC1R*CSSR(KFLAG)
        S1I = SC1I*CSSR(KFLAG)
        S2R = SC2R*CSSR(KFLAG)
        S2I = SC2I*CSSR(KFLAG)
        STR = SC2R
        STI = SC2I
   70   CONTINUE
        PTR = CSPNR*C1R - CSPNI*C1I
        PTI = CSPNR*C1I + CSPNI*C1R
        YR(I) = PTR + CSGNR*C2R - CSGNI*C2I
        YI(I) = PTI + CSGNR*C2I + CSGNI*C2R
        CKR = CKR + RZR
        CKI = CKI + RZI
        CSPNR = -CSPNR
        CSPNI = -CSPNI
        IF (KFLAG.GE.3) GO TO 80
        PTR = DABS(C1R)
        PTI = DABS(C1I)
        C1M = DMAX1(PTR,PTI)
        IF (C1M.LE.BSCLE) GO TO 80
        KFLAG = KFLAG + 1
        BSCLE = BRY(KFLAG)
        S1R = S1R*CSR
        S1I = S1I*CSR
        S2R = STR
        S2I = STI
        S1R = S1R*CSSR(KFLAG)
        S1I = S1I*CSSR(KFLAG)
        S2R = S2R*CSSR(KFLAG)
        S2I = S2I*CSSR(KFLAG)
        CSR = CSRR(KFLAG)
   80 CONTINUE
      RETURN
   90 CONTINUE
      NZ = -1
      IF(NW.EQ.(-2)) NZ=-2
      RETURN
      END
      SUBROUTINE ZBINU(ZR, ZI, FNU, KODE, N, CYR, CYI, NZ, RL, FNUL,
     * TOL, ELIM, ALIM)
C***BEGIN PROLOGUE  ZBINU
C***REFER TO  ZBESH,ZBESI,ZBESJ,ZBESK,ZAIRY,ZBIRY
C
C     ZBINU COMPUTES THE I FUNCTION IN THE RIGHT HALF Z PLANE
C
C***ROUTINES CALLED  ZABS,ZASYI,ZBUNI,ZMLRI,ZSERI,ZUOIK,ZWRSK
C***END PROLOGUE  ZBINU
      EXTERNAL ZABS
      DOUBLE PRECISION ALIM, AZ, CWI, CWR, CYI, CYR, DFNU, ELIM, FNU,
     * FNUL, RL, TOL, ZEROI, ZEROR, ZI, ZR, ZABS
      INTEGER I, INW, KODE, N, NLAST, NN, NUI, NW, NZ
      DIMENSION CYR(N), CYI(N), CWR(2), CWI(2)
      DATA ZEROR,ZEROI / 0.0D0, 0.0D0 /
C
      NZ = 0
      AZ = ZABS(ZR,ZI)
      NN = N
      DFNU = FNU + DBLE(FLOAT(N-1))
      IF (AZ.LE.2.0D0) GO TO 10
      IF (AZ*AZ*0.25D0.GT.DFNU+1.0D0) GO TO 20
   10 CONTINUE
C-----------------------------------------------------------------------
C     POWER SERIES
C-----------------------------------------------------------------------
      CALL ZSERI(ZR, ZI, FNU, KODE, NN, CYR, CYI, NW, TOL, ELIM, ALIM)
      INW = IABS(NW)
      NZ = NZ + INW
      NN = NN - INW
      IF (NN.EQ.0) RETURN
      IF (NW.GE.0) GO TO 120
      DFNU = FNU + DBLE(FLOAT(NN-1))
   20 CONTINUE
      IF (AZ.LT.RL) GO TO 40
      IF (DFNU.LE.1.0D0) GO TO 30
      IF (AZ+AZ.LT.DFNU*DFNU) GO TO 50
C-----------------------------------------------------------------------
C     ASYMPTOTIC EXPANSION FOR LARGE Z
C-----------------------------------------------------------------------
   30 CONTINUE
      CALL ZASYI(ZR, ZI, FNU, KODE, NN, CYR, CYI, NW, RL, TOL, ELIM,
     * ALIM)
      IF (NW.LT.0) GO TO 130
      GO TO 120
   40 CONTINUE
      IF (DFNU.LE.1.0D0) GO TO 70
   50 CONTINUE
C-----------------------------------------------------------------------
C     OVERFLOW AND UNDERFLOW TEST ON I SEQUENCE FOR MILLER ALGORITHM
C-----------------------------------------------------------------------
      CALL ZUOIK(ZR, ZI, FNU, KODE, 1, NN, CYR, CYI, NW, TOL, ELIM,
     * ALIM)
      IF (NW.LT.0) GO TO 130
      NZ = NZ + NW
      NN = NN - NW
      IF (NN.EQ.0) RETURN
      DFNU = FNU+DBLE(FLOAT(NN-1))
      IF (DFNU.GT.FNUL) GO TO 110
      IF (AZ.GT.FNUL) GO TO 110
   60 CONTINUE
      IF (AZ.GT.RL) GO TO 80
   70 CONTINUE
C-----------------------------------------------------------------------
C     MILLER ALGORITHM NORMALIZED BY THE SERIES
C-----------------------------------------------------------------------
      CALL ZMLRI(ZR, ZI, FNU, KODE, NN, CYR, CYI, NW, TOL)
      IF(NW.LT.0) GO TO 130
      GO TO 120
   80 CONTINUE
C-----------------------------------------------------------------------
C     MILLER ALGORITHM NORMALIZED BY THE WRONSKIAN
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C     OVERFLOW TEST ON K FUNCTIONS USED IN WRONSKIAN
C-----------------------------------------------------------------------
      CALL ZUOIK(ZR, ZI, FNU, KODE, 2, 2, CWR, CWI, NW, TOL, ELIM,
     * ALIM)
      IF (NW.GE.0) GO TO 100
      NZ = NN
      DO 90 I=1,NN
        CYR(I) = ZEROR
        CYI(I) = ZEROI
   90 CONTINUE
      RETURN
  100 CONTINUE
      IF (NW.GT.0) GO TO 130
      CALL ZWRSK(ZR, ZI, FNU, KODE, NN, CYR, CYI, NW, CWR, CWI, TOL,
     * ELIM, ALIM)
      IF (NW.LT.0) GO TO 130
      GO TO 120
  110 CONTINUE
C-----------------------------------------------------------------------
C     INCREMENT FNU+NN-1 UP TO FNUL, COMPUTE AND RECUR BACKWARD
C-----------------------------------------------------------------------
      NUI = INT(SNGL(FNUL-DFNU)) + 1
      NUI = MAX0(NUI,0)
      CALL ZBUNI(ZR, ZI, FNU, KODE, NN, CYR, CYI, NW, NUI, NLAST, FNUL,
     * TOL, ELIM, ALIM)
      IF (NW.LT.0) GO TO 130
      NZ = NZ + NW
      IF (NLAST.EQ.0) GO TO 120
      NN = NLAST
      GO TO 60
  120 CONTINUE
      RETURN
  130 CONTINUE
      NZ = -1
      IF(NW.EQ.(-2)) NZ=-2
      RETURN
      END
      DOUBLE PRECISION FUNCTION DGAMLN(Z,IERR)
C***BEGIN PROLOGUE  DGAMLN
C***DATE WRITTEN   830501   (YYMMDD)
C***REVISION DATE  830501   (YYMMDD)
C***CATEGORY NO.  B5F
C***KEYWORDS  GAMMA FUNCTION,LOGARITHM OF GAMMA FUNCTION
C***AUTHOR  AMOS, DONALD E., SANDIA NATIONAL LABORATORIES
C***PURPOSE  TO COMPUTE THE LOGARITHM OF THE GAMMA FUNCTION
C***DESCRIPTION
C
C               **** A DOUBLE PRECISION ROUTINE ****
C         DGAMLN COMPUTES THE NATURAL LOG OF THE GAMMA FUNCTION FOR
C         Z.GT.0.  THE ASYMPTOTIC EXPANSION IS USED TO GENERATE VALUES
C         GREATER THAN ZMIN WHICH ARE ADJUSTED BY THE RECURSION
C         G(Z+1)=Z*G(Z) FOR Z.LE.ZMIN.  THE FUNCTION WAS MADE AS
C         PORTABLE AS POSSIBLE BY COMPUTIMG ZMIN FROM THE NUMBER OF BASE
C         10 DIGITS IN A WORD, RLN=AMAX1(-ALOG10(R1MACH(4)),0.5E-18)
C         LIMITED TO 18 DIGITS OF (RELATIVE) ACCURACY.
C
C         SINCE INTEGER ARGUMENTS ARE COMMON, A TABLE LOOK UP ON 100
C         VALUES IS USED FOR SPEED OF EXECUTION.
C
C     DESCRIPTION OF ARGUMENTS
C
C         INPUT      Z IS D0UBLE PRECISION
C           Z      - ARGUMENT, Z.GT.0.0D0
C
C         OUTPUT      DGAMLN IS DOUBLE PRECISION
C           DGAMLN  - NATURAL LOG OF THE GAMMA FUNCTION AT Z.NE.0.0D0
C           IERR    - ERROR FLAG
C                     IERR=0, NORMAL RETURN, COMPUTATION COMPLETED
C                     IERR=1, Z.LE.0.0D0,    NO COMPUTATION
C
C
C***REFERENCES  COMPUTATION OF BESSEL FUNCTIONS OF COMPLEX ARGUMENT
C                 BY D. E. AMOS, SAND83-0083, MAY, 1983.
C***ROUTINES CALLED  I1MACH,D1MACH
C***END PROLOGUE  DGAMLN
      DOUBLE PRECISION CF, CON, FLN, FZ, GLN, RLN, S, TLG, TRM, TST,
     * T1, WDTOL, Z, ZDMY, ZINC, ZM, ZMIN, ZP, ZSQ, D1MACH
      INTEGER I, IERR, I1M, K, MZ, NZ, I1MACH
      DIMENSION CF(22), GLN(100)
C           LNGAMMA(N), N=1,100
      DATA GLN(1), GLN(2), GLN(3), GLN(4), GLN(5), GLN(6), GLN(7),
     1     GLN(8), GLN(9), GLN(10), GLN(11), GLN(12), GLN(13), GLN(14),
     2     GLN(15), GLN(16), GLN(17), GLN(18), GLN(19), GLN(20),
     3     GLN(21), GLN(22)/
     4     0.00000000000000000D+00,     0.00000000000000000D+00,
     5     6.93147180559945309D-01,     1.79175946922805500D+00,
     6     3.17805383034794562D+00,     4.78749174278204599D+00,
     7     6.57925121201010100D+00,     8.52516136106541430D+00,
     8     1.06046029027452502D+01,     1.28018274800814696D+01,
     9     1.51044125730755153D+01,     1.75023078458738858D+01,
     A     1.99872144956618861D+01,     2.25521638531234229D+01,
     B     2.51912211827386815D+01,     2.78992713838408916D+01,
     C     3.06718601060806728D+01,     3.35050734501368889D+01,
     D     3.63954452080330536D+01,     3.93398841871994940D+01,
     E     4.23356164607534850D+01,     4.53801388984769080D+01/
      DATA GLN(23), GLN(24), GLN(25), GLN(26), GLN(27), GLN(28),
     1     GLN(29), GLN(30), GLN(31), GLN(32), GLN(33), GLN(34),
     2     GLN(35), GLN(36), GLN(37), GLN(38), GLN(39), GLN(40),
     3     GLN(41), GLN(42), GLN(43), GLN(44)/
     4     4.84711813518352239D+01,     5.16066755677643736D+01,
     5     5.47847293981123192D+01,     5.80036052229805199D+01,
     6     6.12617017610020020D+01,     6.45575386270063311D+01,
     7     6.78897431371815350D+01,     7.12570389671680090D+01,
     8     7.46582363488301644D+01,     7.80922235533153106D+01,
     9     8.15579594561150372D+01,     8.50544670175815174D+01,
     A     8.85808275421976788D+01,     9.21361756036870925D+01,
     B     9.57196945421432025D+01,     9.93306124547874269D+01,
     C     1.02968198614513813D+02,     1.06631760260643459D+02,
     D     1.10320639714757395D+02,     1.14034211781461703D+02,
     E     1.17771881399745072D+02,     1.21533081515438634D+02/
      DATA GLN(45), GLN(46), GLN(47), GLN(48), GLN(49), GLN(50),
     1     GLN(51), GLN(52), GLN(53), GLN(54), GLN(55), GLN(56),
     2     GLN(57), GLN(58), GLN(59), GLN(60), GLN(61), GLN(62),
     3     GLN(63), GLN(64), GLN(65), GLN(66)/
     4     1.25317271149356895D+02,     1.29123933639127215D+02,
     5     1.32952575035616310D+02,     1.36802722637326368D+02,
     6     1.40673923648234259D+02,     1.44565743946344886D+02,
     7     1.48477766951773032D+02,     1.52409592584497358D+02,
     8     1.56360836303078785D+02,     1.60331128216630907D+02,
     9     1.64320112263195181D+02,     1.68327445448427652D+02,
     A     1.72352797139162802D+02,     1.76395848406997352D+02,
     B     1.80456291417543771D+02,     1.84533828861449491D+02,
     C     1.88628173423671591D+02,     1.92739047287844902D+02,
     D     1.96866181672889994D+02,     2.01009316399281527D+02,
     E     2.05168199482641199D+02,     2.09342586752536836D+02/
      DATA GLN(67), GLN(68), GLN(69), GLN(70), GLN(71), GLN(72),
     1     GLN(73), GLN(74), GLN(75), GLN(76), GLN(77), GLN(78),
     2     GLN(79), GLN(80), GLN(81), GLN(82), GLN(83), GLN(84),
     3     GLN(85), GLN(86), GLN(87), GLN(88)/
     4     2.13532241494563261D+02,     2.17736934113954227D+02,
     5     2.21956441819130334D+02,     2.26190548323727593D+02,
     6     2.30439043565776952D+02,     2.34701723442818268D+02,
     7     2.38978389561834323D+02,     2.43268849002982714D+02,
     8     2.47572914096186884D+02,     2.51890402209723194D+02,
     9     2.56221135550009525D+02,     2.60564940971863209D+02,
     A     2.64921649798552801D+02,     2.69291097651019823D+02,
     B     2.73673124285693704D+02,     2.78067573440366143D+02,
     C     2.82474292687630396D+02,     2.86893133295426994D+02,
     D     2.91323950094270308D+02,     2.95766601350760624D+02,
     E     3.00220948647014132D+02,     3.04686856765668715D+02/
      DATA GLN(89), GLN(90), GLN(91), GLN(92), GLN(93), GLN(94),
     1     GLN(95), GLN(96), GLN(97), GLN(98), GLN(99), GLN(100)/
     2     3.09164193580146922D+02,     3.13652829949879062D+02,
     3     3.18152639620209327D+02,     3.22663499126726177D+02,
     4     3.27185287703775217D+02,     3.31717887196928473D+02,
     5     3.36261181979198477D+02,     3.40815058870799018D+02,
     6     3.45379407062266854D+02,     3.49954118040770237D+02,
     7     3.54539085519440809D+02,     3.59134205369575399D+02/
C             COEFFICIENTS OF ASYMPTOTIC EXPANSION
      DATA CF(1), CF(2), CF(3), CF(4), CF(5), CF(6), CF(7), CF(8),
     1     CF(9), CF(10), CF(11), CF(12), CF(13), CF(14), CF(15),
     2     CF(16), CF(17), CF(18), CF(19), CF(20), CF(21), CF(22)/
     3     8.33333333333333333D-02,    -2.77777777777777778D-03,
     4     7.93650793650793651D-04,    -5.95238095238095238D-04,
     5     8.41750841750841751D-04,    -1.91752691752691753D-03,
     6     6.41025641025641026D-03,    -2.95506535947712418D-02,
     7     1.79644372368830573D-01,    -1.39243221690590112D+00,
     8     1.34028640441683920D+01,    -1.56848284626002017D+02,
     9     2.19310333333333333D+03,    -3.61087712537249894D+04,
     A     6.91472268851313067D+05,    -1.52382215394074162D+07,
     B     3.82900751391414141D+08,    -1.08822660357843911D+10,
     C     3.47320283765002252D+11,    -1.23696021422692745D+13,
     D     4.88788064793079335D+14,    -2.13203339609193739D+16/
C
C             LN(2*PI)
      DATA CON                    /     1.83787706640934548D+00/
C
C***FIRST EXECUTABLE STATEMENT  DGAMLN
      IERR=0
      IF (Z.LE.0.0D0) GO TO 70
      IF (Z.GT.101.0D0) GO TO 10
      NZ = INT(Z)
      FZ = Z - FLOAT(NZ)
      IF (FZ.GT.0.0D0) GO TO 10
      IF (NZ.GT.100) GO TO 10
      DGAMLN = GLN(NZ)
      RETURN
   10 CONTINUE
      WDTOL = D1MACH(4)
      WDTOL = DMAX1(WDTOL,0.5D-18)
      I1M = I1MACH(14)
      RLN = D1MACH(5)*FLOAT(I1M)
      FLN = DMIN1(RLN,20.0D0)
      FLN = DMAX1(FLN,3.0D0)
      FLN = FLN - 3.0D0
      ZM = 1.8000D0 + 0.3875D0*FLN
      MZ = INT(SNGL(ZM)) + 1
      ZMIN = FLOAT(MZ)
      ZDMY = Z
      ZINC = 0.0D0
      IF (Z.GE.ZMIN) GO TO 20
      ZINC = ZMIN - FLOAT(NZ)
      ZDMY = Z + ZINC
   20 CONTINUE
      ZP = 1.0D0/ZDMY
      T1 = CF(1)*ZP
      S = T1
      IF (ZP.LT.WDTOL) GO TO 40
      ZSQ = ZP*ZP
      TST = T1*WDTOL
      DO 30 K=2,22
        ZP = ZP*ZSQ
        TRM = CF(K)*ZP
        IF (DABS(TRM).LT.TST) GO TO 40
        S = S + TRM
   30 CONTINUE
   40 CONTINUE
      IF (ZINC.NE.0.0D0) GO TO 50
      TLG = DLOG(Z)
      DGAMLN = Z*(TLG-1.0D0) + 0.5D0*(CON-TLG) + S
      RETURN
   50 CONTINUE
      ZP = 1.0D0
      NZ = INT(SNGL(ZINC))
      DO 60 I=1,NZ
        ZP = ZP*(Z+FLOAT(I-1))
   60 CONTINUE
      TLG = DLOG(ZDMY)
      DGAMLN = ZDMY*(TLG-1.0D0) - DLOG(ZP) + 0.5D0*(CON-TLG) + S
      RETURN
C
C
   70 CONTINUE
      IERR=1
      RETURN
      END
      SUBROUTINE ZACAI(ZR, ZI, FNU, KODE, MR, N, YR, YI, NZ, RL, TOL,
     * ELIM, ALIM)
C***BEGIN PROLOGUE  ZACAI
C***REFER TO  ZAIRY
C
C     ZACAI APPLIES THE ANALYTIC CONTINUATION FORMULA
C
C         K(FNU,ZN*EXP(MP))=K(FNU,ZN)*EXP(-MP*FNU) - MP*I(FNU,ZN)
C                 MP=PI*MR*CMPLX(0.0,1.0)
C
C     TO CONTINUE THE K FUNCTION FROM THE RIGHT HALF TO THE LEFT
C     HALF Z PLANE FOR USE WITH ZAIRY WHERE FNU=1/3 OR 2/3 AND N=1.
C     ZACAI IS THE SAME AS ZACON WITH THE PARTS FOR LARGER ORDERS AND
C     RECURRENCE REMOVED. A RECURSIVE CALL TO ZACON CAN RESULT IF ZACON
C     IS CALLED FROM ZAIRY.
C
C***ROUTINES CALLED  ZASYI,ZBKNU,ZMLRI,ZSERI,ZS1S2,D1MACH,ZABS
C***END PROLOGUE  ZACAI
C     COMPLEX CSGN,CSPN,C1,C2,Y,Z,ZN,CY
      EXTERNAL ZABS
      DOUBLE PRECISION ALIM, ARG, ASCLE, AZ, CSGNR, CSGNI, CSPNR,
     * CSPNI, C1R, C1I, C2R, C2I, CYR, CYI, DFNU, ELIM, FMR, FNU, PI,
     * RL, SGN, TOL, YY, YR, YI, ZR, ZI, ZNR, ZNI, D1MACH, ZABS
      INTEGER INU, IUF, KODE, MR, N, NN, NW, NZ
      DIMENSION YR(N), YI(N), CYR(2), CYI(2)
      DATA PI / 3.14159265358979324D0 /
      NZ = 0
      ZNR = -ZR
      ZNI = -ZI
      AZ = ZABS(ZR,ZI)
      NN = N
      DFNU = FNU + DBLE(FLOAT(N-1))
      IF (AZ.LE.2.0D0) GO TO 10
      IF (AZ*AZ*0.25D0.GT.DFNU+1.0D0) GO TO 20
   10 CONTINUE
C-----------------------------------------------------------------------
C     POWER SERIES FOR THE I FUNCTION
C-----------------------------------------------------------------------
      CALL ZSERI(ZNR, ZNI, FNU, KODE, NN, YR, YI, NW, TOL, ELIM, ALIM)
      GO TO 40
   20 CONTINUE
      IF (AZ.LT.RL) GO TO 30
C-----------------------------------------------------------------------
C     ASYMPTOTIC EXPANSION FOR LARGE Z FOR THE I FUNCTION
C-----------------------------------------------------------------------
      CALL ZASYI(ZNR, ZNI, FNU, KODE, NN, YR, YI, NW, RL, TOL, ELIM,
     * ALIM)
      IF (NW.LT.0) GO TO 80
      GO TO 40
   30 CONTINUE
C-----------------------------------------------------------------------
C     MILLER ALGORITHM NORMALIZED BY THE SERIES FOR THE I FUNCTION
C-----------------------------------------------------------------------
      CALL ZMLRI(ZNR, ZNI, FNU, KODE, NN, YR, YI, NW, TOL)
      IF(NW.LT.0) GO TO 80
   40 CONTINUE
C-----------------------------------------------------------------------
C     ANALYTIC CONTINUATION TO THE LEFT HALF PLANE FOR THE K FUNCTION
C-----------------------------------------------------------------------
      CALL ZBKNU(ZNR, ZNI, FNU, KODE, 1, CYR, CYI, NW, TOL, ELIM, ALIM)
      IF (NW.NE.0) GO TO 80
      FMR = DBLE(FLOAT(MR))
      SGN = -DSIGN(PI,FMR)
      CSGNR = 0.0D0
      CSGNI = SGN
      IF (KODE.EQ.1) GO TO 50
      YY = -ZNI
      CSGNR = -CSGNI*DSIN(YY)
      CSGNI = CSGNI*DCOS(YY)
   50 CONTINUE
C-----------------------------------------------------------------------
C     CALCULATE CSPN=EXP(FNU*PI*I) TO MINIMIZE LOSSES OF SIGNIFICANCE
C     WHEN FNU IS LARGE
C-----------------------------------------------------------------------
      INU = INT(SNGL(FNU))
      ARG = (FNU-DBLE(FLOAT(INU)))*SGN
      CSPNR = DCOS(ARG)
      CSPNI = DSIN(ARG)
      IF (MOD(INU,2).EQ.0) GO TO 60
      CSPNR = -CSPNR
      CSPNI = -CSPNI
   60 CONTINUE
      C1R = CYR(1)
      C1I = CYI(1)
      C2R = YR(1)
      C2I = YI(1)
      IF (KODE.EQ.1) GO TO 70
      IUF = 0
      ASCLE = 1.0D+3*D1MACH(1)/TOL
      CALL ZS1S2(ZNR, ZNI, C1R, C1I, C2R, C2I, NW, ASCLE, ALIM, IUF)
      NZ = NZ + NW
   70 CONTINUE
      YR(1) = CSPNR*C1R - CSPNI*C1I + CSGNR*C2R - CSGNI*C2I
      YI(1) = CSPNR*C1I + CSPNI*C1R + CSGNR*C2I + CSGNI*C2R
      RETURN
   80 CONTINUE
      NZ = -1
      IF(NW.EQ.(-2)) NZ=-2
      RETURN
      END
      SUBROUTINE ZUCHK(YR, YI, NZ, ASCLE, TOL)
C***BEGIN PROLOGUE  ZUCHK
C***REFER TO ZSERI,ZUOIK,ZUNK1,ZUNK2,ZUNI1,ZUNI2,ZKSCL
C
C      Y ENTERS AS A SCALED QUANTITY WHOSE MAGNITUDE IS GREATER THAN
C      EXP(-ALIM)=ASCLE=1.0E+3*D1MACH(1)/TOL. THE TEST IS MADE TO SEE
C      IF THE MAGNITUDE OF THE REAL OR IMAGINARY PART WOULD UNDERFLOW
C      WHEN Y IS SCALED (BY TOL) TO ITS PROPER VALUE. Y IS ACCEPTED
C      IF THE UNDERFLOW IS AT LEAST ONE PRECISION BELOW THE MAGNITUDE
C      OF THE LARGEST COMPONENT; OTHERWISE THE PHASE ANGLE DOES NOT HAVE
C      ABSOLUTE ACCURACY AND AN UNDERFLOW IS ASSUMED.
C
C***ROUTINES CALLED  (NONE)
C***END PROLOGUE  ZUCHK
C
C     COMPLEX Y
      DOUBLE PRECISION ASCLE, SS, ST, TOL, WR, WI, YR, YI
      INTEGER NZ
      NZ = 0
      WR = DABS(YR)
      WI = DABS(YI)
      ST = DMIN1(WR,WI)
      IF (ST.GT.ASCLE) RETURN
      SS = DMAX1(WR,WI)
      ST = ST/TOL
      IF (SS.LT.ST) NZ = 1
      RETURN
      END
      SUBROUTINE ZUNIK(ZRR, ZRI, FNU, IKFLG, IPMTR, TOL, INIT, PHIR,
     * PHII, ZETA1R, ZETA1I, ZETA2R, ZETA2I, SUMR, SUMI, CWRKR, CWRKI)
C***BEGIN PROLOGUE  ZUNIK
C***REFER TO  ZBESI,ZBESK
C
C        ZUNIK COMPUTES PARAMETERS FOR THE UNIFORM ASYMPTOTIC
C        EXPANSIONS OF THE I AND K FUNCTIONS ON IKFLG= 1 OR 2
C        RESPECTIVELY BY
C
C        W(FNU,ZR) = PHI*EXP(ZETA)*SUM
C
C        WHERE       ZETA=-ZETA1 + ZETA2       OR
C                          ZETA1 - ZETA2
C
C        THE FIRST CALL MUST HAVE INIT=0. SUBSEQUENT CALLS WITH THE
C        SAME ZR AND FNU WILL RETURN THE I OR K FUNCTION ON IKFLG=
C        1 OR 2 WITH NO CHANGE IN INIT. CWRK IS A COMPLEX WORK
C        ARRAY. IPMTR=0 COMPUTES ALL PARAMETERS. IPMTR=1 COMPUTES PHI,
C        ZETA1,ZETA2.
C
C***ROUTINES CALLED  ZDIV,ZLOG,ZSQRT,D1MACH
C***END PROLOGUE  ZUNIK
C     COMPLEX CFN,CON,CONE,CRFN,CWRK,CZERO,PHI,S,SR,SUM,T,T2,ZETA1,
C    *ZETA2,ZN,ZR
      DOUBLE PRECISION AC, C, CON, CONEI, CONER, CRFNI, CRFNR, CWRKI,
     * CWRKR, FNU, PHII, PHIR, RFN, SI, SR, SRI, SRR, STI, STR, SUMI,
     * SUMR, TEST, TI, TOL, TR, T2I, T2R, ZEROI, ZEROR, ZETA1I, ZETA1R,
     * ZETA2I, ZETA2R, ZNI, ZNR, ZRI, ZRR, D1MACH
      INTEGER I, IDUM, IKFLG, INIT, IPMTR, J, K, L
      DIMENSION C(120), CWRKR(16), CWRKI(16), CON(2)
      DATA ZEROR,ZEROI,CONER,CONEI / 0.0D0, 0.0D0, 1.0D0, 0.0D0 /
      DATA CON(1), CON(2)  /
     1 3.98942280401432678D-01,  1.25331413731550025D+00 /
      DATA C(1), C(2), C(3), C(4), C(5), C(6), C(7), C(8), C(9), C(10),
     1     C(11), C(12), C(13), C(14), C(15), C(16), C(17), C(18),
     2     C(19), C(20), C(21), C(22), C(23), C(24)/
     3     1.00000000000000000D+00,    -2.08333333333333333D-01,
     4     1.25000000000000000D-01,     3.34201388888888889D-01,
     5    -4.01041666666666667D-01,     7.03125000000000000D-02,
     6    -1.02581259645061728D+00,     1.84646267361111111D+00,
     7    -8.91210937500000000D-01,     7.32421875000000000D-02,
     8     4.66958442342624743D+00,    -1.12070026162229938D+01,
     9     8.78912353515625000D+00,    -2.36408691406250000D+00,
     A     1.12152099609375000D-01,    -2.82120725582002449D+01,
     B     8.46362176746007346D+01,    -9.18182415432400174D+01,
     C     4.25349987453884549D+01,    -7.36879435947963170D+00,
     D     2.27108001708984375D-01,     2.12570130039217123D+02,
     E    -7.65252468141181642D+02,     1.05999045252799988D+03/
      DATA C(25), C(26), C(27), C(28), C(29), C(30), C(31), C(32),
     1     C(33), C(34), C(35), C(36), C(37), C(38), C(39), C(40),
     2     C(41), C(42), C(43), C(44), C(45), C(46), C(47), C(48)/
     3    -6.99579627376132541D+02,     2.18190511744211590D+02,
     4    -2.64914304869515555D+01,     5.72501420974731445D-01,
     5    -1.91945766231840700D+03,     8.06172218173730938D+03,
     6    -1.35865500064341374D+04,     1.16553933368645332D+04,
     7    -5.30564697861340311D+03,     1.20090291321635246D+03,
     8    -1.08090919788394656D+02,     1.72772750258445740D+00,
     9     2.02042913309661486D+04,    -9.69805983886375135D+04,
     A     1.92547001232531532D+05,    -2.03400177280415534D+05,
     B     1.22200464983017460D+05,    -4.11926549688975513D+04,
     C     7.10951430248936372D+03,    -4.93915304773088012D+02,
     D     6.07404200127348304D+00,    -2.42919187900551333D+05,
     E     1.31176361466297720D+06,    -2.99801591853810675D+06/
      DATA C(49), C(50), C(51), C(52), C(53), C(54), C(55), C(56),
     1     C(57), C(58), C(59), C(60), C(61), C(62), C(63), C(64),
     2     C(65), C(66), C(67), C(68), C(69), C(70), C(71), C(72)/
     3     3.76327129765640400D+06,    -2.81356322658653411D+06,
     4     1.26836527332162478D+06,    -3.31645172484563578D+05,
     5     4.52187689813627263D+04,    -2.49983048181120962D+03,
     6     2.43805296995560639D+01,     3.28446985307203782D+06,
     7    -1.97068191184322269D+07,     5.09526024926646422D+07,
     8    -7.41051482115326577D+07,     6.63445122747290267D+07,
     9    -3.75671766607633513D+07,     1.32887671664218183D+07,
     A    -2.78561812808645469D+06,     3.08186404612662398D+05,
     B    -1.38860897537170405D+04,     1.10017140269246738D+02,
     C    -4.93292536645099620D+07,     3.25573074185765749D+08,
     D    -9.39462359681578403D+08,     1.55359689957058006D+09,
     E    -1.62108055210833708D+09,     1.10684281682301447D+09/
      DATA C(73), C(74), C(75), C(76), C(77), C(78), C(79), C(80),
     1     C(81), C(82), C(83), C(84), C(85), C(86), C(87), C(88),
     2     C(89), C(90), C(91), C(92), C(93), C(94), C(95), C(96)/
     3    -4.95889784275030309D+08,     1.42062907797533095D+08,
     4    -2.44740627257387285D+07,     2.24376817792244943D+06,
     5    -8.40054336030240853D+04,     5.51335896122020586D+02,
     6     8.14789096118312115D+08,    -5.86648149205184723D+09,
     7     1.86882075092958249D+10,    -3.46320433881587779D+10,
     8     4.12801855797539740D+10,    -3.30265997498007231D+10,
     9     1.79542137311556001D+10,    -6.56329379261928433D+09,
     A     1.55927986487925751D+09,    -2.25105661889415278D+08,
     B     1.73951075539781645D+07,    -5.49842327572288687D+05,
     C     3.03809051092238427D+03,    -1.46792612476956167D+10,
     D     1.14498237732025810D+11,    -3.99096175224466498D+11,
     E     8.19218669548577329D+11,    -1.09837515608122331D+12/
      DATA C(97), C(98), C(99), C(100), C(101), C(102), C(103), C(104),
     1     C(105), C(106), C(107), C(108), C(109), C(110), C(111),
     2     C(112), C(113), C(114), C(115), C(116), C(117), C(118)/
     3     1.00815810686538209D+12,    -6.45364869245376503D+11,
     4     2.87900649906150589D+11,    -8.78670721780232657D+10,
     5     1.76347306068349694D+10,    -2.16716498322379509D+09,
     6     1.43157876718888981D+08,    -3.87183344257261262D+06,
     7     1.82577554742931747D+04,     2.86464035717679043D+11,
     8    -2.40629790002850396D+12,     9.10934118523989896D+12,
     9    -2.05168994109344374D+13,     3.05651255199353206D+13,
     A    -3.16670885847851584D+13,     2.33483640445818409D+13,
     B    -1.23204913055982872D+13,     4.61272578084913197D+12,
     C    -1.19655288019618160D+12,     2.05914503232410016D+11,
     D    -2.18229277575292237D+10,     1.24700929351271032D+09/
      DATA C(119), C(120)/
     1    -2.91883881222208134D+07,     1.18838426256783253D+05/
C
      IF (INIT.NE.0) GO TO 40
C-----------------------------------------------------------------------
C     INITIALIZE ALL VARIABLES
C-----------------------------------------------------------------------
      RFN = 1.0D0/FNU
C-----------------------------------------------------------------------
C     OVERFLOW TEST (ZR/FNU TOO SMALL)
C-----------------------------------------------------------------------
      TEST = D1MACH(1)*1.0D+3
      AC = FNU*TEST
      IF (DABS(ZRR).GT.AC .OR. DABS(ZRI).GT.AC) GO TO 15
      ZETA1R = 2.0D0*DABS(DLOG(TEST))+FNU
      ZETA1I = 0.0D0
      ZETA2R = FNU
      ZETA2I = 0.0D0
      PHIR = 1.0D0
      PHII = 0.0D0
      RETURN
   15 CONTINUE
      TR = ZRR*RFN
      TI = ZRI*RFN
      SR = CONER + (TR*TR-TI*TI)
      SI = CONEI + (TR*TI+TI*TR)
      CALL ZSQRT(SR, SI, SRR, SRI)
      STR = CONER + SRR
      STI = CONEI + SRI
      CALL ZDIV(STR, STI, TR, TI, ZNR, ZNI)
      CALL ZLOG(ZNR, ZNI, STR, STI, IDUM)
      ZETA1R = FNU*STR
      ZETA1I = FNU*STI
      ZETA2R = FNU*SRR
      ZETA2I = FNU*SRI
      CALL ZDIV(CONER, CONEI, SRR, SRI, TR, TI)
      SRR = TR*RFN
      SRI = TI*RFN
      CALL ZSQRT(SRR, SRI, CWRKR(16), CWRKI(16))
      PHIR = CWRKR(16)*CON(IKFLG)
      PHII = CWRKI(16)*CON(IKFLG)
      IF (IPMTR.NE.0) RETURN
      CALL ZDIV(CONER, CONEI, SR, SI, T2R, T2I)
      CWRKR(1) = CONER
      CWRKI(1) = CONEI
      CRFNR = CONER
      CRFNI = CONEI
      AC = 1.0D0
      L = 1
      DO 20 K=2,15
        SR = ZEROR
        SI = ZEROI
        DO 10 J=1,K
          L = L + 1
          STR = SR*T2R - SI*T2I + C(L)
          SI = SR*T2I + SI*T2R
          SR = STR
   10   CONTINUE
        STR = CRFNR*SRR - CRFNI*SRI
        CRFNI = CRFNR*SRI + CRFNI*SRR
        CRFNR = STR
        CWRKR(K) = CRFNR*SR - CRFNI*SI
        CWRKI(K) = CRFNR*SI + CRFNI*SR
        AC = AC*RFN
        TEST = DABS(CWRKR(K)) + DABS(CWRKI(K))
        IF (AC.LT.TOL .AND. TEST.LT.TOL) GO TO 30
   20 CONTINUE
      K = 15
   30 CONTINUE
      INIT = K
   40 CONTINUE
      IF (IKFLG.EQ.2) GO TO 60
C-----------------------------------------------------------------------
C     COMPUTE SUM FOR THE I FUNCTION
C-----------------------------------------------------------------------
      SR = ZEROR
      SI = ZEROI
      DO 50 I=1,INIT
        SR = SR + CWRKR(I)
        SI = SI + CWRKI(I)
   50 CONTINUE
      SUMR = SR
      SUMI = SI
      PHIR = CWRKR(16)*CON(1)
      PHII = CWRKI(16)*CON(1)
      RETURN
   60 CONTINUE
C-----------------------------------------------------------------------
C     COMPUTE SUM FOR THE K FUNCTION
C-----------------------------------------------------------------------
      SR = ZEROR
      SI = ZEROI
      TR = CONER
      DO 70 I=1,INIT
        SR = SR + TR*CWRKR(I)
        SI = SI + TR*CWRKI(I)
        TR = -TR
   70 CONTINUE
      SUMR = SR
      SUMI = SI
      PHIR = CWRKR(16)*CON(2)
      PHII = CWRKI(16)*CON(2)
      RETURN
      END
      SUBROUTINE ZUNHJ(ZR, ZI, FNU, IPMTR, TOL, PHIR, PHII, ARGR, ARGI,
     * ZETA1R, ZETA1I, ZETA2R, ZETA2I, ASUMR, ASUMI, BSUMR, BSUMI)
C***BEGIN PROLOGUE  ZUNHJ
C***REFER TO  ZBESI,ZBESK
C
C     REFERENCES
C         HANDBOOK OF MATHEMATICAL FUNCTIONS BY M. ABRAMOWITZ AND I.A.
C         STEGUN, AMS55, NATIONAL BUREAU OF STANDARDS, 1965, CHAPTER 9.
C
C         ASYMPTOTICS AND SPECIAL FUNCTIONS BY F.W.J. OLVER, ACADEMIC
C         PRESS, N.Y., 1974, PAGE 420
C
C     ABSTRACT
C         ZUNHJ COMPUTES PARAMETERS FOR BESSEL FUNCTIONS C(FNU,Z) =
C         J(FNU,Z), Y(FNU,Z) OR H(I,FNU,Z) I=1,2 FOR LARGE ORDERS FNU
C         BY MEANS OF THE UNIFORM ASYMPTOTIC EXPANSION
C
C         C(FNU,Z)=C1*PHI*( ASUM*AIRY(ARG) + C2*BSUM*DAIRY(ARG) )
C
C         FOR PROPER CHOICES OF C1, C2, AIRY AND DAIRY WHERE AIRY IS
C         AN AIRY FUNCTION AND DAIRY IS ITS DERIVATIVE.
C
C               (2/3)*FNU*ZETA**1.5 = ZETA1-ZETA2,
C
C         ZETA1=0.5*FNU*CLOG((1+W)/(1-W)), ZETA2=FNU*W FOR SCALING
C         PURPOSES IN AIRY FUNCTIONS FROM CAIRY OR CBIRY.
C
C         MCONJ=SIGN OF AIMAG(Z), BUT IS AMBIGUOUS WHEN Z IS REAL AND
C         MUST BE SPECIFIED. IPMTR=0 RETURNS ALL PARAMETERS. IPMTR=
C         1 COMPUTES ALL EXCEPT ASUM AND BSUM.
C
C***ROUTINES CALLED  ZABS,ZDIV,ZLOG,ZSQRT,D1MACH
C***END PROLOGUE  ZUNHJ
C     COMPLEX ARG,ASUM,BSUM,CFNU,CONE,CR,CZERO,DR,P,PHI,PRZTH,PTFN,
C    *RFN13,RTZTA,RZTH,SUMA,SUMB,TFN,T2,UP,W,W2,Z,ZA,ZB,ZC,ZETA,ZETA1,
C    *ZETA2,ZTH
      EXTERNAL ZABS
      DOUBLE PRECISION ALFA, ANG, AP, AR, ARGI, ARGR, ASUMI, ASUMR,
     * ATOL, AW2, AZTH, BETA, BR, BSUMI, BSUMR, BTOL, C, CONEI, CONER,
     * CRI, CRR, DRI, DRR, EX1, EX2, FNU, FN13, FN23, GAMA, GPI, HPI,
     * PHII, PHIR, PI, PP, PR, PRZTHI, PRZTHR, PTFNI, PTFNR, RAW, RAW2,
     * RAZTH, RFNU, RFNU2, RFN13, RTZTI, RTZTR, RZTHI, RZTHR, STI, STR,
     * SUMAI, SUMAR, SUMBI, SUMBR, TEST, TFNI, TFNR, THPI, TOL, TZAI,
     * TZAR, T2I, T2R, UPI, UPR, WI, WR, W2I, W2R, ZAI, ZAR, ZBI, ZBR,
     * ZCI, ZCR, ZEROI, ZEROR, ZETAI, ZETAR, ZETA1I, ZETA1R, ZETA2I,
     * ZETA2R, ZI, ZR, ZTHI, ZTHR, ZABS, AC, D1MACH
      INTEGER IAS, IBS, IPMTR, IS, J, JR, JU, K, KMAX, KP1, KS, L, LR,
     * LRP1, L1, L2, M, IDUM
      DIMENSION AR(14), BR(14), C(105), ALFA(180), BETA(210), GAMA(30),
     * AP(30), PR(30), PI(30), UPR(14), UPI(14), CRR(14), CRI(14),
     * DRR(14), DRI(14)
      DATA AR(1), AR(2), AR(3), AR(4), AR(5), AR(6), AR(7), AR(8),
     1     AR(9), AR(10), AR(11), AR(12), AR(13), AR(14)/
     2     1.00000000000000000D+00,     1.04166666666666667D-01,
     3     8.35503472222222222D-02,     1.28226574556327160D-01,
     4     2.91849026464140464D-01,     8.81627267443757652D-01,
     5     3.32140828186276754D+00,     1.49957629868625547D+01,
     6     7.89230130115865181D+01,     4.74451538868264323D+02,
     7     3.20749009089066193D+03,     2.40865496408740049D+04,
     8     1.98923119169509794D+05,     1.79190200777534383D+06/
      DATA BR(1), BR(2), BR(3), BR(4), BR(5), BR(6), BR(7), BR(8),
     1     BR(9), BR(10), BR(11), BR(12), BR(13), BR(14)/
     2     1.00000000000000000D+00,    -1.45833333333333333D-01,
     3    -9.87413194444444444D-02,    -1.43312053915895062D-01,
     4    -3.17227202678413548D-01,    -9.42429147957120249D-01,
     5    -3.51120304082635426D+00,    -1.57272636203680451D+01,
     6    -8.22814390971859444D+01,    -4.92355370523670524D+02,
     7    -3.31621856854797251D+03,    -2.48276742452085896D+04,
     8    -2.04526587315129788D+05,    -1.83844491706820990D+06/
      DATA C(1), C(2), C(3), C(4), C(5), C(6), C(7), C(8), C(9), C(10),
     1     C(11), C(12), C(13), C(14), C(15), C(16), C(17), C(18),
     2     C(19), C(20), C(21), C(22), C(23), C(24)/
     3     1.00000000000000000D+00,    -2.08333333333333333D-01,
     4     1.25000000000000000D-01,     3.34201388888888889D-01,
     5    -4.01041666666666667D-01,     7.03125000000000000D-02,
     6    -1.02581259645061728D+00,     1.84646267361111111D+00,
     7    -8.91210937500000000D-01,     7.32421875000000000D-02,
     8     4.66958442342624743D+00,    -1.12070026162229938D+01,
     9     8.78912353515625000D+00,    -2.36408691406250000D+00,
     A     1.12152099609375000D-01,    -2.82120725582002449D+01,
     B     8.46362176746007346D+01,    -9.18182415432400174D+01,
     C     4.25349987453884549D+01,    -7.36879435947963170D+00,
     D     2.27108001708984375D-01,     2.12570130039217123D+02,
     E    -7.65252468141181642D+02,     1.05999045252799988D+03/
      DATA C(25), C(26), C(27), C(28), C(29), C(30), C(31), C(32),
     1     C(33), C(34), C(35), C(36), C(37), C(38), C(39), C(40),
     2     C(41), C(42), C(43), C(44), C(45), C(46), C(47), C(48)/
     3    -6.99579627376132541D+02,     2.18190511744211590D+02,
     4    -2.64914304869515555D+01,     5.72501420974731445D-01,
     5    -1.91945766231840700D+03,     8.06172218173730938D+03,
     6    -1.35865500064341374D+04,     1.16553933368645332D+04,
     7    -5.30564697861340311D+03,     1.20090291321635246D+03,
     8    -1.08090919788394656D+02,     1.72772750258445740D+00,
     9     2.02042913309661486D+04,    -9.69805983886375135D+04,
     A     1.92547001232531532D+05,    -2.03400177280415534D+05,
     B     1.22200464983017460D+05,    -4.11926549688975513D+04,
     C     7.10951430248936372D+03,    -4.93915304773088012D+02,
     D     6.07404200127348304D+00,    -2.42919187900551333D+05,
     E     1.31176361466297720D+06,    -2.99801591853810675D+06/
      DATA C(49), C(50), C(51), C(52), C(53), C(54), C(55), C(56),
     1     C(57), C(58), C(59), C(60), C(61), C(62), C(63), C(64),
     2     C(65), C(66), C(67), C(68), C(69), C(70), C(71), C(72)/
     3     3.76327129765640400D+06,    -2.81356322658653411D+06,
     4     1.26836527332162478D+06,    -3.31645172484563578D+05,
     5     4.52187689813627263D+04,    -2.49983048181120962D+03,
     6     2.43805296995560639D+01,     3.28446985307203782D+06,
     7    -1.97068191184322269D+07,     5.09526024926646422D+07,
     8    -7.41051482115326577D+07,     6.63445122747290267D+07,
     9    -3.75671766607633513D+07,     1.32887671664218183D+07,
     A    -2.78561812808645469D+06,     3.08186404612662398D+05,
     B    -1.38860897537170405D+04,     1.10017140269246738D+02,
     C    -4.93292536645099620D+07,     3.25573074185765749D+08,
     D    -9.39462359681578403D+08,     1.55359689957058006D+09,
     E    -1.62108055210833708D+09,     1.10684281682301447D+09/
      DATA C(73), C(74), C(75), C(76), C(77), C(78), C(79), C(80),
     1     C(81), C(82), C(83), C(84), C(85), C(86), C(87), C(88),
     2     C(89), C(90), C(91), C(92), C(93), C(94), C(95), C(96)/
     3    -4.95889784275030309D+08,     1.42062907797533095D+08,
     4    -2.44740627257387285D+07,     2.24376817792244943D+06,
     5    -8.40054336030240853D+04,     5.51335896122020586D+02,
     6     8.14789096118312115D+08,    -5.86648149205184723D+09,
     7     1.86882075092958249D+10,    -3.46320433881587779D+10,
     8     4.12801855797539740D+10,    -3.30265997498007231D+10,
     9     1.79542137311556001D+10,    -6.56329379261928433D+09,
     A     1.55927986487925751D+09,    -2.25105661889415278D+08,
     B     1.73951075539781645D+07,    -5.49842327572288687D+05,
     C     3.03809051092238427D+03,    -1.46792612476956167D+10,
     D     1.14498237732025810D+11,    -3.99096175224466498D+11,
     E     8.19218669548577329D+11,    -1.09837515608122331D+12/
      DATA C(97), C(98), C(99), C(100), C(101), C(102), C(103), C(104),
     1     C(105)/
     2     1.00815810686538209D+12,    -6.45364869245376503D+11,
     3     2.87900649906150589D+11,    -8.78670721780232657D+10,
     4     1.76347306068349694D+10,    -2.16716498322379509D+09,
     5     1.43157876718888981D+08,    -3.87183344257261262D+06,
     6     1.82577554742931747D+04/
      DATA ALFA(1), ALFA(2), ALFA(3), ALFA(4), ALFA(5), ALFA(6),
     1     ALFA(7), ALFA(8), ALFA(9), ALFA(10), ALFA(11), ALFA(12),
     2     ALFA(13), ALFA(14), ALFA(15), ALFA(16), ALFA(17), ALFA(18),
     3     ALFA(19), ALFA(20), ALFA(21), ALFA(22)/
     4    -4.44444444444444444D-03,    -9.22077922077922078D-04,
     5    -8.84892884892884893D-05,     1.65927687832449737D-04,
     6     2.46691372741792910D-04,     2.65995589346254780D-04,
     7     2.61824297061500945D-04,     2.48730437344655609D-04,
     8     2.32721040083232098D-04,     2.16362485712365082D-04,
     9     2.00738858762752355D-04,     1.86267636637545172D-04,
     A     1.73060775917876493D-04,     1.61091705929015752D-04,
     B     1.50274774160908134D-04,     1.40503497391269794D-04,
     C     1.31668816545922806D-04,     1.23667445598253261D-04,
     D     1.16405271474737902D-04,     1.09798298372713369D-04,
     E     1.03772410422992823D-04,     9.82626078369363448D-05/
      DATA ALFA(23), ALFA(24), ALFA(25), ALFA(26), ALFA(27), ALFA(28),
     1     ALFA(29), ALFA(30), ALFA(31), ALFA(32), ALFA(33), ALFA(34),
     2     ALFA(35), ALFA(36), ALFA(37), ALFA(38), ALFA(39), ALFA(40),
     3     ALFA(41), ALFA(42), ALFA(43), ALFA(44)/
     4     9.32120517249503256D-05,     8.85710852478711718D-05,
     5     8.42963105715700223D-05,     8.03497548407791151D-05,
     6     7.66981345359207388D-05,     7.33122157481777809D-05,
     7     7.01662625163141333D-05,     6.72375633790160292D-05,
     8     6.93735541354588974D-04,     2.32241745182921654D-04,
     9    -1.41986273556691197D-05,    -1.16444931672048640D-04,
     A    -1.50803558053048762D-04,    -1.55121924918096223D-04,
     B    -1.46809756646465549D-04,    -1.33815503867491367D-04,
     C    -1.19744975684254051D-04,    -1.06184319207974020D-04,
     D    -9.37699549891194492D-05,    -8.26923045588193274D-05,
     E    -7.29374348155221211D-05,    -6.44042357721016283D-05/
      DATA ALFA(45), ALFA(46), ALFA(47), ALFA(48), ALFA(49), ALFA(50),
     1     ALFA(51), ALFA(52), ALFA(53), ALFA(54), ALFA(55), ALFA(56),
     2     ALFA(57), ALFA(58), ALFA(59), ALFA(60), ALFA(61), ALFA(62),
     3     ALFA(63), ALFA(64), ALFA(65), ALFA(66)/
     4    -5.69611566009369048D-05,    -5.04731044303561628D-05,
     5    -4.48134868008882786D-05,    -3.98688727717598864D-05,
     6    -3.55400532972042498D-05,    -3.17414256609022480D-05,
     7    -2.83996793904174811D-05,    -2.54522720634870566D-05,
     8    -2.28459297164724555D-05,    -2.05352753106480604D-05,
     9    -1.84816217627666085D-05,    -1.66519330021393806D-05,
     A    -1.50179412980119482D-05,    -1.35554031379040526D-05,
     B    -1.22434746473858131D-05,    -1.10641884811308169D-05,
     C    -3.54211971457743841D-04,    -1.56161263945159416D-04,
     D     3.04465503594936410D-05,     1.30198655773242693D-04,
     E     1.67471106699712269D-04,     1.70222587683592569D-04/
      DATA ALFA(67), ALFA(68), ALFA(69), ALFA(70), ALFA(71), ALFA(72),
     1     ALFA(73), ALFA(74), ALFA(75), ALFA(76), ALFA(77), ALFA(78),
     2     ALFA(79), ALFA(80), ALFA(81), ALFA(82), ALFA(83), ALFA(84),
     3     ALFA(85), ALFA(86), ALFA(87), ALFA(88)/
     4     1.56501427608594704D-04,     1.36339170977445120D-04,
     5     1.14886692029825128D-04,     9.45869093034688111D-05,
     6     7.64498419250898258D-05,     6.07570334965197354D-05,
     7     4.74394299290508799D-05,     3.62757512005344297D-05,
     8     2.69939714979224901D-05,     1.93210938247939253D-05,
     9     1.30056674793963203D-05,     7.82620866744496661D-06,
     A     3.59257485819351583D-06,     1.44040049814251817D-07,
     B    -2.65396769697939116D-06,    -4.91346867098485910D-06,
     C    -6.72739296091248287D-06,    -8.17269379678657923D-06,
     D    -9.31304715093561232D-06,    -1.02011418798016441D-05,
     E    -1.08805962510592880D-05,    -1.13875481509603555D-05/
      DATA ALFA(89), ALFA(90), ALFA(91), ALFA(92), ALFA(93), ALFA(94),
     1     ALFA(95), ALFA(96), ALFA(97), ALFA(98), ALFA(99), ALFA(100),
     2     ALFA(101), ALFA(102), ALFA(103), ALFA(104), ALFA(105),
     3     ALFA(106), ALFA(107), ALFA(108), ALFA(109), ALFA(110)/
     4    -1.17519675674556414D-05,    -1.19987364870944141D-05,
     5     3.78194199201772914D-04,     2.02471952761816167D-04,
     6    -6.37938506318862408D-05,    -2.38598230603005903D-04,
     7    -3.10916256027361568D-04,    -3.13680115247576316D-04,
     8    -2.78950273791323387D-04,    -2.28564082619141374D-04,
     9    -1.75245280340846749D-04,    -1.25544063060690348D-04,
     A    -8.22982872820208365D-05,    -4.62860730588116458D-05,
     B    -1.72334302366962267D-05,     5.60690482304602267D-06,
     C     2.31395443148286800D-05,     3.62642745856793957D-05,
     D     4.58006124490188752D-05,     5.24595294959114050D-05,
     E     5.68396208545815266D-05,     5.94349820393104052D-05/
      DATA ALFA(111), ALFA(112), ALFA(113), ALFA(114), ALFA(115),
     1     ALFA(116), ALFA(117), ALFA(118), ALFA(119), ALFA(120),
     2     ALFA(121), ALFA(122), ALFA(123), ALFA(124), ALFA(125),
     3     ALFA(126), ALFA(127), ALFA(128), ALFA(129), ALFA(130)/
     4     6.06478527578421742D-05,     6.08023907788436497D-05,
     5     6.01577894539460388D-05,     5.89199657344698500D-05,
     6     5.72515823777593053D-05,     5.52804375585852577D-05,
     7     5.31063773802880170D-05,     5.08069302012325706D-05,
     8     4.84418647620094842D-05,     4.60568581607475370D-05,
     9    -6.91141397288294174D-04,    -4.29976633058871912D-04,
     A     1.83067735980039018D-04,     6.60088147542014144D-04,
     B     8.75964969951185931D-04,     8.77335235958235514D-04,
     C     7.49369585378990637D-04,     5.63832329756980918D-04,
     D     3.68059319971443156D-04,     1.88464535514455599D-04/
      DATA ALFA(131), ALFA(132), ALFA(133), ALFA(134), ALFA(135),
     1     ALFA(136), ALFA(137), ALFA(138), ALFA(139), ALFA(140),
     2     ALFA(141), ALFA(142), ALFA(143), ALFA(144), ALFA(145),
     3     ALFA(146), ALFA(147), ALFA(148), ALFA(149), ALFA(150)/
     4     3.70663057664904149D-05,    -8.28520220232137023D-05,
     5    -1.72751952869172998D-04,    -2.36314873605872983D-04,
     6    -2.77966150694906658D-04,    -3.02079514155456919D-04,
     7    -3.12594712643820127D-04,    -3.12872558758067163D-04,
     8    -3.05678038466324377D-04,    -2.93226470614557331D-04,
     9    -2.77255655582934777D-04,    -2.59103928467031709D-04,
     A    -2.39784014396480342D-04,    -2.20048260045422848D-04,
     B    -2.00443911094971498D-04,    -1.81358692210970687D-04,
     C    -1.63057674478657464D-04,    -1.45712672175205844D-04,
     D    -1.29425421983924587D-04,    -1.14245691942445952D-04/
      DATA ALFA(151), ALFA(152), ALFA(153), ALFA(154), ALFA(155),
     1     ALFA(156), ALFA(157), ALFA(158), ALFA(159), ALFA(160),
     2     ALFA(161), ALFA(162), ALFA(163), ALFA(164), ALFA(165),
     3     ALFA(166), ALFA(167), ALFA(168), ALFA(169), ALFA(170)/
     4     1.92821964248775885D-03,     1.35592576302022234D-03,
     5    -7.17858090421302995D-04,    -2.58084802575270346D-03,
     6    -3.49271130826168475D-03,    -3.46986299340960628D-03,
     7    -2.82285233351310182D-03,    -1.88103076404891354D-03,
     8    -8.89531718383947600D-04,     3.87912102631035228D-06,
     9     7.28688540119691412D-04,     1.26566373053457758D-03,
     A     1.62518158372674427D-03,     1.83203153216373172D-03,
     B     1.91588388990527909D-03,     1.90588846755546138D-03,
     C     1.82798982421825727D-03,     1.70389506421121530D-03,
     D     1.55097127171097686D-03,     1.38261421852276159D-03/
      DATA ALFA(171), ALFA(172), ALFA(173), ALFA(174), ALFA(175),
     1     ALFA(176), ALFA(177), ALFA(178), ALFA(179), ALFA(180)/
     2     1.20881424230064774D-03,     1.03676532638344962D-03,
     3     8.71437918068619115D-04,     7.16080155297701002D-04,
     4     5.72637002558129372D-04,     4.42089819465802277D-04,
     5     3.24724948503090564D-04,     2.20342042730246599D-04,
     6     1.28412898401353882D-04,     4.82005924552095464D-05/
      DATA BETA(1), BETA(2), BETA(3), BETA(4), BETA(5), BETA(6),
     1     BETA(7), BETA(8), BETA(9), BETA(10), BETA(11), BETA(12),
     2     BETA(13), BETA(14), BETA(15), BETA(16), BETA(17), BETA(18),
     3     BETA(19), BETA(20), BETA(21), BETA(22)/
     4     1.79988721413553309D-02,     5.59964911064388073D-03,
     5     2.88501402231132779D-03,     1.80096606761053941D-03,
     6     1.24753110589199202D-03,     9.22878876572938311D-04,
     7     7.14430421727287357D-04,     5.71787281789704872D-04,
     8     4.69431007606481533D-04,     3.93232835462916638D-04,
     9     3.34818889318297664D-04,     2.88952148495751517D-04,
     A     2.52211615549573284D-04,     2.22280580798883327D-04,
     B     1.97541838033062524D-04,     1.76836855019718004D-04,
     C     1.59316899661821081D-04,     1.44347930197333986D-04,
     D     1.31448068119965379D-04,     1.20245444949302884D-04,
     E     1.10449144504599392D-04,     1.01828770740567258D-04/
      DATA BETA(23), BETA(24), BETA(25), BETA(26), BETA(27), BETA(28),
     1     BETA(29), BETA(30), BETA(31), BETA(32), BETA(33), BETA(34),
     2     BETA(35), BETA(36), BETA(37), BETA(38), BETA(39), BETA(40),
     3     BETA(41), BETA(42), BETA(43), BETA(44)/
     4     9.41998224204237509D-05,     8.74130545753834437D-05,
     5     8.13466262162801467D-05,     7.59002269646219339D-05,
     6     7.09906300634153481D-05,     6.65482874842468183D-05,
     7     6.25146958969275078D-05,     5.88403394426251749D-05,
     8    -1.49282953213429172D-03,    -8.78204709546389328D-04,
     9    -5.02916549572034614D-04,    -2.94822138512746025D-04,
     A    -1.75463996970782828D-04,    -1.04008550460816434D-04,
     B    -5.96141953046457895D-05,    -3.12038929076098340D-05,
     C    -1.26089735980230047D-05,    -2.42892608575730389D-07,
     D     8.05996165414273571D-06,     1.36507009262147391D-05,
     E     1.73964125472926261D-05,     1.98672978842133780D-05/
      DATA BETA(45), BETA(46), BETA(47), BETA(48), BETA(49), BETA(50),
     1     BETA(51), BETA(52), BETA(53), BETA(54), BETA(55), BETA(56),
     2     BETA(57), BETA(58), BETA(59), BETA(60), BETA(61), BETA(62),
     3     BETA(63), BETA(64), BETA(65), BETA(66)/
     4     2.14463263790822639D-05,     2.23954659232456514D-05,
     5     2.28967783814712629D-05,     2.30785389811177817D-05,
     6     2.30321976080909144D-05,     2.28236073720348722D-05,
     7     2.25005881105292418D-05,     2.20981015361991429D-05,
     8     2.16418427448103905D-05,     2.11507649256220843D-05,
     9     2.06388749782170737D-05,     2.01165241997081666D-05,
     A     1.95913450141179244D-05,     1.90689367910436740D-05,
     B     1.85533719641636667D-05,     1.80475722259674218D-05,
     C     5.52213076721292790D-04,     4.47932581552384646D-04,
     D     2.79520653992020589D-04,     1.52468156198446602D-04,
     E     6.93271105657043598D-05,     1.76258683069991397D-05/
      DATA BETA(67), BETA(68), BETA(69), BETA(70), BETA(71), BETA(72),
     1     BETA(73), BETA(74), BETA(75), BETA(76), BETA(77), BETA(78),
     2     BETA(79), BETA(80), BETA(81), BETA(82), BETA(83), BETA(84),
     3     BETA(85), BETA(86), BETA(87), BETA(88)/
     4    -1.35744996343269136D-05,    -3.17972413350427135D-05,
     5    -4.18861861696693365D-05,    -4.69004889379141029D-05,
     6    -4.87665447413787352D-05,    -4.87010031186735069D-05,
     7    -4.74755620890086638D-05,    -4.55813058138628452D-05,
     8    -4.33309644511266036D-05,    -4.09230193157750364D-05,
     9    -3.84822638603221274D-05,    -3.60857167535410501D-05,
     A    -3.37793306123367417D-05,    -3.15888560772109621D-05,
     B    -2.95269561750807315D-05,    -2.75978914828335759D-05,
     C    -2.58006174666883713D-05,    -2.41308356761280200D-05,
     D    -2.25823509518346033D-05,    -2.11479656768912971D-05,
     E    -1.98200638885294927D-05,    -1.85909870801065077D-05/
      DATA BETA(89), BETA(90), BETA(91), BETA(92), BETA(93), BETA(94),
     1     BETA(95), BETA(96), BETA(97), BETA(98), BETA(99), BETA(100),
     2     BETA(101), BETA(102), BETA(103), BETA(104), BETA(105),
     3     BETA(106), BETA(107), BETA(108), BETA(109), BETA(110)/
     4    -1.74532699844210224D-05,    -1.63997823854497997D-05,
     5    -4.74617796559959808D-04,    -4.77864567147321487D-04,
     6    -3.20390228067037603D-04,    -1.61105016119962282D-04,
     7    -4.25778101285435204D-05,     3.44571294294967503D-05,
     8     7.97092684075674924D-05,     1.03138236708272200D-04,
     9     1.12466775262204158D-04,     1.13103642108481389D-04,
     A     1.08651634848774268D-04,     1.01437951597661973D-04,
     B     9.29298396593363896D-05,     8.40293133016089978D-05,
     C     7.52727991349134062D-05,     6.69632521975730872D-05,
     D     5.92564547323194704D-05,     5.22169308826975567D-05,
     E     4.58539485165360646D-05,     4.01445513891486808D-05/
      DATA BETA(111), BETA(112), BETA(113), BETA(114), BETA(115),
     1     BETA(116), BETA(117), BETA(118), BETA(119), BETA(120),
     2     BETA(121), BETA(122), BETA(123), BETA(124), BETA(125),
     3     BETA(126), BETA(127), BETA(128), BETA(129), BETA(130)/
     4     3.50481730031328081D-05,     3.05157995034346659D-05,
     5     2.64956119950516039D-05,     2.29363633690998152D-05,
     6     1.97893056664021636D-05,     1.70091984636412623D-05,
     7     1.45547428261524004D-05,     1.23886640995878413D-05,
     8     1.04775876076583236D-05,     8.79179954978479373D-06,
     9     7.36465810572578444D-04,     8.72790805146193976D-04,
     A     6.22614862573135066D-04,     2.85998154194304147D-04,
     B     3.84737672879366102D-06,    -1.87906003636971558D-04,
     C    -2.97603646594554535D-04,    -3.45998126832656348D-04,
     D    -3.53382470916037712D-04,    -3.35715635775048757D-04/
      DATA BETA(131), BETA(132), BETA(133), BETA(134), BETA(135),
     1     BETA(136), BETA(137), BETA(138), BETA(139), BETA(140),
     2     BETA(141), BETA(142), BETA(143), BETA(144), BETA(145),
     3     BETA(146), BETA(147), BETA(148), BETA(149), BETA(150)/
     4    -3.04321124789039809D-04,    -2.66722723047612821D-04,
     5    -2.27654214122819527D-04,    -1.89922611854562356D-04,
     6    -1.55058918599093870D-04,    -1.23778240761873630D-04,
     7    -9.62926147717644187D-05,    -7.25178327714425337D-05,
     8    -5.22070028895633801D-05,    -3.50347750511900522D-05,
     9    -2.06489761035551757D-05,    -8.70106096849767054D-06,
     A     1.13698686675100290D-06,     9.16426474122778849D-06,
     B     1.56477785428872620D-05,     2.08223629482466847D-05,
     C     2.48923381004595156D-05,     2.80340509574146325D-05,
     D     3.03987774629861915D-05,     3.21156731406700616D-05/
      DATA BETA(151), BETA(152), BETA(153), BETA(154), BETA(155),
     1     BETA(156), BETA(157), BETA(158), BETA(159), BETA(160),
     2     BETA(161), BETA(162), BETA(163), BETA(164), BETA(165),
     3     BETA(166), BETA(167), BETA(168), BETA(169), BETA(170)/
     4    -1.80182191963885708D-03,    -2.43402962938042533D-03,
     5    -1.83422663549856802D-03,    -7.62204596354009765D-04,
     6     2.39079475256927218D-04,     9.49266117176881141D-04,
     7     1.34467449701540359D-03,     1.48457495259449178D-03,
     8     1.44732339830617591D-03,     1.30268261285657186D-03,
     9     1.10351597375642682D-03,     8.86047440419791759D-04,
     A     6.73073208165665473D-04,     4.77603872856582378D-04,
     B     3.05991926358789362D-04,     1.60315694594721630D-04,
     C     4.00749555270613286D-05,    -5.66607461635251611D-05,
     D    -1.32506186772982638D-04,    -1.90296187989614057D-04/
      DATA BETA(171), BETA(172), BETA(173), BETA(174), BETA(175),
     1     BETA(176), BETA(177), BETA(178), BETA(179), BETA(180),
     2     BETA(181), BETA(182), BETA(183), BETA(184), BETA(185),
     3     BETA(186), BETA(187), BETA(188), BETA(189), BETA(190)/
     4    -2.32811450376937408D-04,    -2.62628811464668841D-04,
     5    -2.82050469867598672D-04,    -2.93081563192861167D-04,
     6    -2.97435962176316616D-04,    -2.96557334239348078D-04,
     7    -2.91647363312090861D-04,    -2.83696203837734166D-04,
     8    -2.73512317095673346D-04,    -2.61750155806768580D-04,
     9     6.38585891212050914D-03,     9.62374215806377941D-03,
     A     7.61878061207001043D-03,     2.83219055545628054D-03,
     B    -2.09841352012720090D-03,    -5.73826764216626498D-03,
     C    -7.70804244495414620D-03,    -8.21011692264844401D-03,
     D    -7.65824520346905413D-03,    -6.47209729391045177D-03/
      DATA BETA(191), BETA(192), BETA(193), BETA(194), BETA(195),
     1     BETA(196), BETA(197), BETA(198), BETA(199), BETA(200),
     2     BETA(201), BETA(202), BETA(203), BETA(204), BETA(205),
     3     BETA(206), BETA(207), BETA(208), BETA(209), BETA(210)/
     4    -4.99132412004966473D-03,    -3.45612289713133280D-03,
     5    -2.01785580014170775D-03,    -7.59430686781961401D-04,
     6     2.84173631523859138D-04,     1.10891667586337403D-03,
     7     1.72901493872728771D-03,     2.16812590802684701D-03,
     8     2.45357710494539735D-03,     2.61281821058334862D-03,
     9     2.67141039656276912D-03,     2.65203073395980430D-03,
     A     2.57411652877287315D-03,     2.45389126236094427D-03,
     B     2.30460058071795494D-03,     2.13684837686712662D-03,
     C     1.95896528478870911D-03,     1.77737008679454412D-03,
     D     1.59690280765839059D-03,     1.42111975664438546D-03/
      DATA GAMA(1), GAMA(2), GAMA(3), GAMA(4), GAMA(5), GAMA(6),
     1     GAMA(7), GAMA(8), GAMA(9), GAMA(10), GAMA(11), GAMA(12),
     2     GAMA(13), GAMA(14), GAMA(15), GAMA(16), GAMA(17), GAMA(18),
     3     GAMA(19), GAMA(20), GAMA(21), GAMA(22)/
     4     6.29960524947436582D-01,     2.51984209978974633D-01,
     5     1.54790300415655846D-01,     1.10713062416159013D-01,
     6     8.57309395527394825D-02,     6.97161316958684292D-02,
     7     5.86085671893713576D-02,     5.04698873536310685D-02,
     8     4.42600580689154809D-02,     3.93720661543509966D-02,
     9     3.54283195924455368D-02,     3.21818857502098231D-02,
     A     2.94646240791157679D-02,     2.71581677112934479D-02,
     B     2.51768272973861779D-02,     2.34570755306078891D-02,
     C     2.19508390134907203D-02,     2.06210828235646240D-02,
     D     1.94388240897880846D-02,     1.83810633800683158D-02,
     E     1.74293213231963172D-02,     1.65685837786612353D-02/
      DATA GAMA(23), GAMA(24), GAMA(25), GAMA(26), GAMA(27), GAMA(28),
     1     GAMA(29), GAMA(30)/
     2     1.57865285987918445D-02,     1.50729501494095594D-02,
     3     1.44193250839954639D-02,     1.38184805735341786D-02,
     4     1.32643378994276568D-02,     1.27517121970498651D-02,
     5     1.22761545318762767D-02,     1.18338262398482403D-02/
      DATA EX1, EX2, HPI, GPI, THPI /
     1     3.33333333333333333D-01,     6.66666666666666667D-01,
     2     1.57079632679489662D+00,     3.14159265358979324D+00,
     3     4.71238898038468986D+00/
      DATA ZEROR,ZEROI,CONER,CONEI / 0.0D0, 0.0D0, 1.0D0, 0.0D0 /
C
      RFNU = 1.0D0/FNU
C-----------------------------------------------------------------------
C     OVERFLOW TEST (Z/FNU TOO SMALL)
C-----------------------------------------------------------------------
      TEST = D1MACH(1)*1.0D+3
      AC = FNU*TEST
      IF (DABS(ZR).GT.AC .OR. DABS(ZI).GT.AC) GO TO 15
      ZETA1R = 2.0D0*DABS(DLOG(TEST))+FNU
      ZETA1I = 0.0D0
      ZETA2R = FNU
      ZETA2I = 0.0D0
      PHIR = 1.0D0
      PHII = 0.0D0
      ARGR = 1.0D0
      ARGI = 0.0D0
      RETURN
   15 CONTINUE
      ZBR = ZR*RFNU
      ZBI = ZI*RFNU
      RFNU2 = RFNU*RFNU
C-----------------------------------------------------------------------
C     COMPUTE IN THE FOURTH QUADRANT
C-----------------------------------------------------------------------
      FN13 = FNU**EX1
      FN23 = FN13*FN13
      RFN13 = 1.0D0/FN13
      W2R = CONER - ZBR*ZBR + ZBI*ZBI
      W2I = CONEI - ZBR*ZBI - ZBR*ZBI
      AW2 = ZABS(W2R,W2I)
      IF (AW2.GT.0.25D0) GO TO 130
C-----------------------------------------------------------------------
C     POWER SERIES FOR CABS(W2).LE.0.25D0
C-----------------------------------------------------------------------
      K = 1
      PR(1) = CONER
      PI(1) = CONEI
      SUMAR = GAMA(1)
      SUMAI = ZEROI
      AP(1) = 1.0D0
      IF (AW2.LT.TOL) GO TO 20
      DO 10 K=2,30
        PR(K) = PR(K-1)*W2R - PI(K-1)*W2I
        PI(K) = PR(K-1)*W2I + PI(K-1)*W2R
        SUMAR = SUMAR + PR(K)*GAMA(K)
        SUMAI = SUMAI + PI(K)*GAMA(K)
        AP(K) = AP(K-1)*AW2
        IF (AP(K).LT.TOL) GO TO 20
   10 CONTINUE
      K = 30
   20 CONTINUE
      KMAX = K
      ZETAR = W2R*SUMAR - W2I*SUMAI
      ZETAI = W2R*SUMAI + W2I*SUMAR
      ARGR = ZETAR*FN23
      ARGI = ZETAI*FN23
      CALL ZSQRT(SUMAR, SUMAI, ZAR, ZAI)
      CALL ZSQRT(W2R, W2I, STR, STI)
      ZETA2R = STR*FNU
      ZETA2I = STI*FNU
      STR = CONER + EX2*(ZETAR*ZAR-ZETAI*ZAI)
      STI = CONEI + EX2*(ZETAR*ZAI+ZETAI*ZAR)
      ZETA1R = STR*ZETA2R - STI*ZETA2I
      ZETA1I = STR*ZETA2I + STI*ZETA2R
      ZAR = ZAR + ZAR
      ZAI = ZAI + ZAI
      CALL ZSQRT(ZAR, ZAI, STR, STI)
      PHIR = STR*RFN13
      PHII = STI*RFN13
      IF (IPMTR.EQ.1) GO TO 120
C-----------------------------------------------------------------------
C     SUM SERIES FOR ASUM AND BSUM
C-----------------------------------------------------------------------
      SUMBR = ZEROR
      SUMBI = ZEROI
      DO 30 K=1,KMAX
        SUMBR = SUMBR + PR(K)*BETA(K)
        SUMBI = SUMBI + PI(K)*BETA(K)
   30 CONTINUE
      ASUMR = ZEROR
      ASUMI = ZEROI
      BSUMR = SUMBR
      BSUMI = SUMBI
      L1 = 0
      L2 = 30
      BTOL = TOL*(DABS(BSUMR)+DABS(BSUMI))
      ATOL = TOL
      PP = 1.0D0
      IAS = 0
      IBS = 0
      IF (RFNU2.LT.TOL) GO TO 110
      DO 100 IS=2,7
        ATOL = ATOL/RFNU2
        PP = PP*RFNU2
        IF (IAS.EQ.1) GO TO 60
        SUMAR = ZEROR
        SUMAI = ZEROI
        DO 40 K=1,KMAX
          M = L1 + K
          SUMAR = SUMAR + PR(K)*ALFA(M)
          SUMAI = SUMAI + PI(K)*ALFA(M)
          IF (AP(K).LT.ATOL) GO TO 50
   40   CONTINUE
   50   CONTINUE
        ASUMR = ASUMR + SUMAR*PP
        ASUMI = ASUMI + SUMAI*PP
        IF (PP.LT.TOL) IAS = 1
   60   CONTINUE
        IF (IBS.EQ.1) GO TO 90
        SUMBR = ZEROR
        SUMBI = ZEROI
        DO 70 K=1,KMAX
          M = L2 + K
          SUMBR = SUMBR + PR(K)*BETA(M)
          SUMBI = SUMBI + PI(K)*BETA(M)
          IF (AP(K).LT.ATOL) GO TO 80
   70   CONTINUE
   80   CONTINUE
        BSUMR = BSUMR + SUMBR*PP
        BSUMI = BSUMI + SUMBI*PP
        IF (PP.LT.BTOL) IBS = 1
   90   CONTINUE
        IF (IAS.EQ.1 .AND. IBS.EQ.1) GO TO 110
        L1 = L1 + 30
        L2 = L2 + 30
  100 CONTINUE
  110 CONTINUE
      ASUMR = ASUMR + CONER
      PP = RFNU*RFN13
      BSUMR = BSUMR*PP
      BSUMI = BSUMI*PP
  120 CONTINUE
      RETURN
C-----------------------------------------------------------------------
C     CABS(W2).GT.0.25D0
C-----------------------------------------------------------------------
  130 CONTINUE
      CALL ZSQRT(W2R, W2I, WR, WI)
      IF (WR.LT.0.0D0) WR = 0.0D0
      IF (WI.LT.0.0D0) WI = 0.0D0
      STR = CONER + WR
      STI = WI
      CALL ZDIV(STR, STI, ZBR, ZBI, ZAR, ZAI)
      CALL ZLOG(ZAR, ZAI, ZCR, ZCI, IDUM)
      IF (ZCI.LT.0.0D0) ZCI = 0.0D0
      IF (ZCI.GT.HPI) ZCI = HPI
      IF (ZCR.LT.0.0D0) ZCR = 0.0D0
      ZTHR = (ZCR-WR)*1.5D0
      ZTHI = (ZCI-WI)*1.5D0
      ZETA1R = ZCR*FNU
      ZETA1I = ZCI*FNU
      ZETA2R = WR*FNU
      ZETA2I = WI*FNU
      AZTH = ZABS(ZTHR,ZTHI)
      ANG = THPI
      IF (ZTHR.GE.0.0D0 .AND. ZTHI.LT.0.0D0) GO TO 140
      ANG = HPI
      IF (ZTHR.EQ.0.0D0) GO TO 140
      ANG = DATAN(ZTHI/ZTHR)
      IF (ZTHR.LT.0.0D0) ANG = ANG + GPI
  140 CONTINUE
      PP = AZTH**EX2
      ANG = ANG*EX2
      ZETAR = PP*DCOS(ANG)
      ZETAI = PP*DSIN(ANG)
      IF (ZETAI.LT.0.0D0) ZETAI = 0.0D0
      ARGR = ZETAR*FN23
      ARGI = ZETAI*FN23
      CALL ZDIV(ZTHR, ZTHI, ZETAR, ZETAI, RTZTR, RTZTI)
      CALL ZDIV(RTZTR, RTZTI, WR, WI, ZAR, ZAI)
      TZAR = ZAR + ZAR
      TZAI = ZAI + ZAI
      CALL ZSQRT(TZAR, TZAI, STR, STI)
      PHIR = STR*RFN13
      PHII = STI*RFN13
      IF (IPMTR.EQ.1) GO TO 120
      RAW = 1.0D0/DSQRT(AW2)
      STR = WR*RAW
      STI = -WI*RAW
      TFNR = STR*RFNU*RAW
      TFNI = STI*RFNU*RAW
      RAZTH = 1.0D0/AZTH
      STR = ZTHR*RAZTH
      STI = -ZTHI*RAZTH
      RZTHR = STR*RAZTH*RFNU
      RZTHI = STI*RAZTH*RFNU
      ZCR = RZTHR*AR(2)
      ZCI = RZTHI*AR(2)
      RAW2 = 1.0D0/AW2
      STR = W2R*RAW2
      STI = -W2I*RAW2
      T2R = STR*RAW2
      T2I = STI*RAW2
      STR = T2R*C(2) + C(3)
      STI = T2I*C(2)
      UPR(2) = STR*TFNR - STI*TFNI
      UPI(2) = STR*TFNI + STI*TFNR
      BSUMR = UPR(2) + ZCR
      BSUMI = UPI(2) + ZCI
      ASUMR = ZEROR
      ASUMI = ZEROI
      IF (RFNU.LT.TOL) GO TO 220
      PRZTHR = RZTHR
      PRZTHI = RZTHI
      PTFNR = TFNR
      PTFNI = TFNI
      UPR(1) = CONER
      UPI(1) = CONEI
      PP = 1.0D0
      BTOL = TOL*(DABS(BSUMR)+DABS(BSUMI))
      KS = 0
      KP1 = 2
      L = 3
      IAS = 0
      IBS = 0
      DO 210 LR=2,12,2
        LRP1 = LR + 1
C-----------------------------------------------------------------------
C     COMPUTE TWO ADDITIONAL CR, DR, AND UP FOR TWO MORE TERMS IN
C     NEXT SUMA AND SUMB
C-----------------------------------------------------------------------
        DO 160 K=LR,LRP1
          KS = KS + 1
          KP1 = KP1 + 1
          L = L + 1
          ZAR = C(L)
          ZAI = ZEROI
          DO 150 J=2,KP1
            L = L + 1
            STR = ZAR*T2R - T2I*ZAI + C(L)
            ZAI = ZAR*T2I + ZAI*T2R
            ZAR = STR
  150     CONTINUE
          STR = PTFNR*TFNR - PTFNI*TFNI
          PTFNI = PTFNR*TFNI + PTFNI*TFNR
          PTFNR = STR
          UPR(KP1) = PTFNR*ZAR - PTFNI*ZAI
          UPI(KP1) = PTFNI*ZAR + PTFNR*ZAI
          CRR(KS) = PRZTHR*BR(KS+1)
          CRI(KS) = PRZTHI*BR(KS+1)
          STR = PRZTHR*RZTHR - PRZTHI*RZTHI
          PRZTHI = PRZTHR*RZTHI + PRZTHI*RZTHR
          PRZTHR = STR
          DRR(KS) = PRZTHR*AR(KS+2)
          DRI(KS) = PRZTHI*AR(KS+2)
  160   CONTINUE
        PP = PP*RFNU2
        IF (IAS.EQ.1) GO TO 180
        SUMAR = UPR(LRP1)
        SUMAI = UPI(LRP1)
        JU = LRP1
        DO 170 JR=1,LR
          JU = JU - 1
          SUMAR = SUMAR + CRR(JR)*UPR(JU) - CRI(JR)*UPI(JU)
          SUMAI = SUMAI + CRR(JR)*UPI(JU) + CRI(JR)*UPR(JU)
  170   CONTINUE
        ASUMR = ASUMR + SUMAR
        ASUMI = ASUMI + SUMAI
        TEST = DABS(SUMAR) + DABS(SUMAI)
        IF (PP.LT.TOL .AND. TEST.LT.TOL) IAS = 1
  180   CONTINUE
        IF (IBS.EQ.1) GO TO 200
        SUMBR = UPR(LR+2) + UPR(LRP1)*ZCR - UPI(LRP1)*ZCI
        SUMBI = UPI(LR+2) + UPR(LRP1)*ZCI + UPI(LRP1)*ZCR
        JU = LRP1
        DO 190 JR=1,LR
          JU = JU - 1
          SUMBR = SUMBR + DRR(JR)*UPR(JU) - DRI(JR)*UPI(JU)
          SUMBI = SUMBI + DRR(JR)*UPI(JU) + DRI(JR)*UPR(JU)
  190   CONTINUE
        BSUMR = BSUMR + SUMBR
        BSUMI = BSUMI + SUMBI
        TEST = DABS(SUMBR) + DABS(SUMBI)
        IF (PP.LT.BTOL .AND. TEST.LT.BTOL) IBS = 1
  200   CONTINUE
        IF (IAS.EQ.1 .AND. IBS.EQ.1) GO TO 220
  210 CONTINUE
  220 CONTINUE
      ASUMR = ASUMR + CONER
      STR = -BSUMR*RFN13
      STI = -BSUMI*RFN13
      CALL ZDIV(STR, STI, RTZTR, RTZTI, BSUMR, BSUMI)
      GO TO 120
      END
      SUBROUTINE ZUNK1(ZR, ZI, FNU, KODE, MR, N, YR, YI, NZ, TOL, ELIM,
     * ALIM)
C***BEGIN PROLOGUE  ZUNK1
C***REFER TO  ZBESK
C
C     ZUNK1 COMPUTES K(FNU,Z) AND ITS ANALYTIC CONTINUATION FROM THE
C     RIGHT HALF PLANE TO THE LEFT HALF PLANE BY MEANS OF THE
C     UNIFORM ASYMPTOTIC EXPANSION.
C     MR INDICATES THE DIRECTION OF ROTATION FOR ANALYTIC CONTINUATION.
C     NZ=-1 MEANS AN OVERFLOW WILL OCCUR
C
C***ROUTINES CALLED  ZKSCL,ZS1S2,ZUCHK,ZUNIK,D1MACH,ZABS
C***END PROLOGUE  ZUNK1
C     COMPLEX CFN,CK,CONE,CRSC,CS,CSCL,CSGN,CSPN,CSR,CSS,CWRK,CY,CZERO,
C    *C1,C2,PHI,PHID,RZ,SUM,SUMD,S1,S2,Y,Z,ZETA1,ZETA1D,ZETA2,ZETA2D,ZR
      EXTERNAL ZABS
      DOUBLE PRECISION ALIM, ANG, APHI, ASC, ASCLE, BRY, CKI, CKR,
     * CONER, CRSC, CSCL, CSGNI, CSPNI, CSPNR, CSR, CSRR, CSSR,
     * CWRKI, CWRKR, CYI, CYR, C1I, C1R, C2I, C2M, C2R, ELIM, FMR, FN,
     * FNF, FNU, PHIDI, PHIDR, PHII, PHIR, PI, RAST, RAZR, RS1, RZI,
     * RZR, SGN, STI, STR, SUMDI, SUMDR, SUMI, SUMR, S1I, S1R, S2I,
     * S2R, TOL, YI, YR, ZEROI, ZEROR, ZETA1I, ZETA1R, ZETA2I, ZETA2R,
     * ZET1DI, ZET1DR, ZET2DI, ZET2DR, ZI, ZR, ZRI, ZRR, D1MACH, ZABS
      INTEGER I, IB, IFLAG, IFN, IL, INIT, INU, IUF, K, KDFLG, KFLAG,
     * KK, KODE, MR, N, NW, NZ, INITD, IC, IPARD, J, M
      DIMENSION BRY(3), INIT(2), YR(N), YI(N), SUMR(2), SUMI(2),
     * ZETA1R(2), ZETA1I(2), ZETA2R(2), ZETA2I(2), CYR(2), CYI(2),
     * CWRKR(16,3), CWRKI(16,3), CSSR(3), CSRR(3), PHIR(2), PHII(2)
      DATA ZEROR,ZEROI,CONER / 0.0D0, 0.0D0, 1.0D0 /
      DATA PI / 3.14159265358979324D0 /
C
      KDFLG = 1
      NZ = 0
C-----------------------------------------------------------------------
C     EXP(-ALIM)=EXP(-ELIM)/TOL=APPROX. ONE PRECISION GREATER THAN
C     THE UNDERFLOW LIMIT
C-----------------------------------------------------------------------
      CSCL = 1.0D0/TOL
      CRSC = TOL
      CSSR(1) = CSCL
      CSSR(2) = CONER
      CSSR(3) = CRSC
      CSRR(1) = CRSC
      CSRR(2) = CONER
      CSRR(3) = CSCL
      BRY(1) = 1.0D+3*D1MACH(1)/TOL
      BRY(2) = 1.0D0/BRY(1)
      BRY(3) = D1MACH(2)
      ZRR = ZR
      ZRI = ZI
      IF (ZR.GE.0.0D0) GO TO 10
      ZRR = -ZR
      ZRI = -ZI
   10 CONTINUE
      J = 2
      DO 70 I=1,N
C-----------------------------------------------------------------------
C     J FLIP FLOPS BETWEEN 1 AND 2 IN J = 3 - J
C-----------------------------------------------------------------------
        J = 3 - J
        FN = FNU + DBLE(FLOAT(I-1))
        INIT(J) = 0
        CALL ZUNIK(ZRR, ZRI, FN, 2, 0, TOL, INIT(J), PHIR(J), PHII(J),
     *   ZETA1R(J), ZETA1I(J), ZETA2R(J), ZETA2I(J), SUMR(J), SUMI(J),
     *   CWRKR(1,J), CWRKI(1,J))
        IF (KODE.EQ.1) GO TO 20
        STR = ZRR + ZETA2R(J)
        STI = ZRI + ZETA2I(J)
        RAST = FN/ZABS(STR,STI)
        STR = STR*RAST*RAST
        STI = -STI*RAST*RAST
        S1R = ZETA1R(J) - STR
        S1I = ZETA1I(J) - STI
        GO TO 30
   20   CONTINUE
        S1R = ZETA1R(J) - ZETA2R(J)
        S1I = ZETA1I(J) - ZETA2I(J)
   30   CONTINUE
        RS1 = S1R
C-----------------------------------------------------------------------
C     TEST FOR UNDERFLOW AND OVERFLOW
C-----------------------------------------------------------------------
        IF (DABS(RS1).GT.ELIM) GO TO 60
        IF (KDFLG.EQ.1) KFLAG = 2
        IF (DABS(RS1).LT.ALIM) GO TO 40
C-----------------------------------------------------------------------
C     REFINE  TEST AND SCALE
C-----------------------------------------------------------------------
        APHI = ZABS(PHIR(J),PHII(J))
        RS1 = RS1 + DLOG(APHI)
        IF (DABS(RS1).GT.ELIM) GO TO 60
        IF (KDFLG.EQ.1) KFLAG = 1
        IF (RS1.LT.0.0D0) GO TO 40
        IF (KDFLG.EQ.1) KFLAG = 3
   40   CONTINUE
C-----------------------------------------------------------------------
C     SCALE S1 TO KEEP INTERMEDIATE ARITHMETIC ON SCALE NEAR
C     EXPONENT EXTREMES
C-----------------------------------------------------------------------
        S2R = PHIR(J)*SUMR(J) - PHII(J)*SUMI(J)
        S2I = PHIR(J)*SUMI(J) + PHII(J)*SUMR(J)
        STR = DEXP(S1R)*CSSR(KFLAG)
        S1R = STR*DCOS(S1I)
        S1I = STR*DSIN(S1I)
        STR = S2R*S1R - S2I*S1I
        S2I = S1R*S2I + S2R*S1I
        S2R = STR
        IF (KFLAG.NE.1) GO TO 50
        CALL ZUCHK(S2R, S2I, NW, BRY(1), TOL)
        IF (NW.NE.0) GO TO 60
   50   CONTINUE
        CYR(KDFLG) = S2R
        CYI(KDFLG) = S2I
        YR(I) = S2R*CSRR(KFLAG)
        YI(I) = S2I*CSRR(KFLAG)
        IF (KDFLG.EQ.2) GO TO 75
        KDFLG = 2
        GO TO 70
   60   CONTINUE
        IF (RS1.GT.0.0D0) GO TO 300
C-----------------------------------------------------------------------
C     FOR ZR.LT.0.0, THE I FUNCTION TO BE ADDED WILL OVERFLOW
C-----------------------------------------------------------------------
        IF (ZR.LT.0.0D0) GO TO 300
        KDFLG = 1
        YR(I)=ZEROR
        YI(I)=ZEROI
        NZ=NZ+1
        IF (I.EQ.1) GO TO 70
        IF ((YR(I-1).EQ.ZEROR).AND.(YI(I-1).EQ.ZEROI)) GO TO 70
        YR(I-1)=ZEROR
        YI(I-1)=ZEROI
        NZ=NZ+1
   70 CONTINUE
      I = N
   75 CONTINUE
      RAZR = 1.0D0/ZABS(ZRR,ZRI)
      STR = ZRR*RAZR
      STI = -ZRI*RAZR
      RZR = (STR+STR)*RAZR
      RZI = (STI+STI)*RAZR
      CKR = FN*RZR
      CKI = FN*RZI
      IB = I + 1
      IF (N.LT.IB) GO TO 160
C-----------------------------------------------------------------------
C     TEST LAST MEMBER FOR UNDERFLOW AND OVERFLOW. SET SEQUENCE TO ZERO
C     ON UNDERFLOW.
C-----------------------------------------------------------------------
      FN = FNU + DBLE(FLOAT(N-1))
      IPARD = 1
      IF (MR.NE.0) IPARD = 0
      INITD = 0
      CALL ZUNIK(ZRR, ZRI, FN, 2, IPARD, TOL, INITD, PHIDR, PHIDI,
     * ZET1DR, ZET1DI, ZET2DR, ZET2DI, SUMDR, SUMDI, CWRKR(1,3),
     * CWRKI(1,3))
      IF (KODE.EQ.1) GO TO 80
      STR = ZRR + ZET2DR
      STI = ZRI + ZET2DI
      RAST = FN/ZABS(STR,STI)
      STR = STR*RAST*RAST
      STI = -STI*RAST*RAST
      S1R = ZET1DR - STR
      S1I = ZET1DI - STI
      GO TO 90
   80 CONTINUE
      S1R = ZET1DR - ZET2DR
      S1I = ZET1DI - ZET2DI
   90 CONTINUE
      RS1 = S1R
      IF (DABS(RS1).GT.ELIM) GO TO 95
      IF (DABS(RS1).LT.ALIM) GO TO 100
C----------------------------------------------------------------------------
C     REFINE ESTIMATE AND TEST
C-------------------------------------------------------------------------
      APHI = ZABS(PHIDR,PHIDI)
      RS1 = RS1+DLOG(APHI)
      IF (DABS(RS1).LT.ELIM) GO TO 100
   95 CONTINUE
      IF (DABS(RS1).GT.0.0D0) GO TO 300
C-----------------------------------------------------------------------
C     FOR ZR.LT.0.0, THE I FUNCTION TO BE ADDED WILL OVERFLOW
C-----------------------------------------------------------------------
      IF (ZR.LT.0.0D0) GO TO 300
      NZ = N
      DO 96 I=1,N
        YR(I) = ZEROR
        YI(I) = ZEROI
   96 CONTINUE
      RETURN
C---------------------------------------------------------------------------
C     FORWARD RECUR FOR REMAINDER OF THE SEQUENCE
C----------------------------------------------------------------------------
  100 CONTINUE
      S1R = CYR(1)
      S1I = CYI(1)
      S2R = CYR(2)
      S2I = CYI(2)
      C1R = CSRR(KFLAG)
      ASCLE = BRY(KFLAG)
      DO 120 I=IB,N
        C2R = S2R
        C2I = S2I
        S2R = CKR*C2R - CKI*C2I + S1R
        S2I = CKR*C2I + CKI*C2R + S1I
        S1R = C2R
        S1I = C2I
        CKR = CKR + RZR
        CKI = CKI + RZI
        C2R = S2R*C1R
        C2I = S2I*C1R
        YR(I) = C2R
        YI(I) = C2I
        IF (KFLAG.GE.3) GO TO 120
        STR = DABS(C2R)
        STI = DABS(C2I)
        C2M = DMAX1(STR,STI)
        IF (C2M.LE.ASCLE) GO TO 120
        KFLAG = KFLAG + 1
        ASCLE = BRY(KFLAG)
        S1R = S1R*C1R
        S1I = S1I*C1R
        S2R = C2R
        S2I = C2I
        S1R = S1R*CSSR(KFLAG)
        S1I = S1I*CSSR(KFLAG)
        S2R = S2R*CSSR(KFLAG)
        S2I = S2I*CSSR(KFLAG)
        C1R = CSRR(KFLAG)
  120 CONTINUE
  160 CONTINUE
      IF (MR.EQ.0) RETURN
C-----------------------------------------------------------------------
C     ANALYTIC CONTINUATION FOR RE(Z).LT.0.0D0
C-----------------------------------------------------------------------
      NZ = 0
      FMR = DBLE(FLOAT(MR))
      SGN = -DSIGN(PI,FMR)
C-----------------------------------------------------------------------
C     CSPN AND CSGN ARE COEFF OF K AND I FUNCTIONS RESP.
C-----------------------------------------------------------------------
      CSGNI = SGN
      INU = INT(SNGL(FNU))
      FNF = FNU - DBLE(FLOAT(INU))
      IFN = INU + N - 1
      ANG = FNF*SGN
      CSPNR = DCOS(ANG)
      CSPNI = DSIN(ANG)
      IF (MOD(IFN,2).EQ.0) GO TO 170
      CSPNR = -CSPNR
      CSPNI = -CSPNI
  170 CONTINUE
      ASC = BRY(1)
      IUF = 0
      KK = N
      KDFLG = 1
      IB = IB - 1
      IC = IB - 1
      DO 270 K=1,N
        FN = FNU + DBLE(FLOAT(KK-1))
C-----------------------------------------------------------------------
C     LOGIC TO SORT OUT CASES WHOSE PARAMETERS WERE SET FOR THE K
C     FUNCTION ABOVE
C-----------------------------------------------------------------------
        M=3
        IF (N.GT.2) GO TO 175
  172   CONTINUE
        INITD = INIT(J)
        PHIDR = PHIR(J)
        PHIDI = PHII(J)
        ZET1DR = ZETA1R(J)
        ZET1DI = ZETA1I(J)
        ZET2DR = ZETA2R(J)
        ZET2DI = ZETA2I(J)
        SUMDR = SUMR(J)
        SUMDI = SUMI(J)
        M = J
        J = 3 - J
        GO TO 180
  175   CONTINUE
        IF ((KK.EQ.N).AND.(IB.LT.N)) GO TO 180
        IF ((KK.EQ.IB).OR.(KK.EQ.IC)) GO TO 172
        INITD = 0
  180   CONTINUE
        CALL ZUNIK(ZRR, ZRI, FN, 1, 0, TOL, INITD, PHIDR, PHIDI,
     *   ZET1DR, ZET1DI, ZET2DR, ZET2DI, SUMDR, SUMDI,
     *   CWRKR(1,M), CWRKI(1,M))
        IF (KODE.EQ.1) GO TO 200
        STR = ZRR + ZET2DR
        STI = ZRI + ZET2DI
        RAST = FN/ZABS(STR,STI)
        STR = STR*RAST*RAST
        STI = -STI*RAST*RAST
        S1R = -ZET1DR + STR
        S1I = -ZET1DI + STI
        GO TO 210
  200   CONTINUE
        S1R = -ZET1DR + ZET2DR
        S1I = -ZET1DI + ZET2DI
  210   CONTINUE
C-----------------------------------------------------------------------
C     TEST FOR UNDERFLOW AND OVERFLOW
C-----------------------------------------------------------------------
        RS1 = S1R
        IF (DABS(RS1).GT.ELIM) GO TO 260
        IF (KDFLG.EQ.1) IFLAG = 2
        IF (DABS(RS1).LT.ALIM) GO TO 220
C-----------------------------------------------------------------------
C     REFINE  TEST AND SCALE
C-----------------------------------------------------------------------
        APHI = ZABS(PHIDR,PHIDI)
        RS1 = RS1 + DLOG(APHI)
        IF (DABS(RS1).GT.ELIM) GO TO 260
        IF (KDFLG.EQ.1) IFLAG = 1
        IF (RS1.LT.0.0D0) GO TO 220
        IF (KDFLG.EQ.1) IFLAG = 3
  220   CONTINUE
        STR = PHIDR*SUMDR - PHIDI*SUMDI
        STI = PHIDR*SUMDI + PHIDI*SUMDR
        S2R = -CSGNI*STI
        S2I = CSGNI*STR
        STR = DEXP(S1R)*CSSR(IFLAG)
        S1R = STR*DCOS(S1I)
        S1I = STR*DSIN(S1I)
        STR = S2R*S1R - S2I*S1I
        S2I = S2R*S1I + S2I*S1R
        S2R = STR
        IF (IFLAG.NE.1) GO TO 230
        CALL ZUCHK(S2R, S2I, NW, BRY(1), TOL)
        IF (NW.EQ.0) GO TO 230
        S2R = ZEROR
        S2I = ZEROI
  230   CONTINUE
        CYR(KDFLG) = S2R
        CYI(KDFLG) = S2I
        C2R = S2R
        C2I = S2I
        S2R = S2R*CSRR(IFLAG)
        S2I = S2I*CSRR(IFLAG)
C-----------------------------------------------------------------------
C     ADD I AND K FUNCTIONS, K SEQUENCE IN Y(I), I=1,N
C-----------------------------------------------------------------------
        S1R = YR(KK)
        S1I = YI(KK)
        IF (KODE.EQ.1) GO TO 250
        CALL ZS1S2(ZRR, ZRI, S1R, S1I, S2R, S2I, NW, ASC, ALIM, IUF)
        NZ = NZ + NW
  250   CONTINUE
        YR(KK) = S1R*CSPNR - S1I*CSPNI + S2R
        YI(KK) = CSPNR*S1I + CSPNI*S1R + S2I
        KK = KK - 1
        CSPNR = -CSPNR
        CSPNI = -CSPNI
        IF (C2R.NE.0.0D0 .OR. C2I.NE.0.0D0) GO TO 255
        KDFLG = 1
        GO TO 270
  255   CONTINUE
        IF (KDFLG.EQ.2) GO TO 275
        KDFLG = 2
        GO TO 270
  260   CONTINUE
        IF (RS1.GT.0.0D0) GO TO 300
        S2R = ZEROR
        S2I = ZEROI
        GO TO 230
  270 CONTINUE
      K = N
  275 CONTINUE
      IL = N - K
      IF (IL.EQ.0) RETURN
C-----------------------------------------------------------------------
C     RECUR BACKWARD FOR REMAINDER OF I SEQUENCE AND ADD IN THE
C     K FUNCTIONS, SCALING THE I SEQUENCE DURING RECURRENCE TO KEEP
C     INTERMEDIATE ARITHMETIC ON SCALE NEAR EXPONENT EXTREMES.
C-----------------------------------------------------------------------
      S1R = CYR(1)
      S1I = CYI(1)
      S2R = CYR(2)
      S2I = CYI(2)
      CSR = CSRR(IFLAG)
      ASCLE = BRY(IFLAG)
      FN = DBLE(FLOAT(INU+IL))
      DO 290 I=1,IL
        C2R = S2R
        C2I = S2I
        S2R = S1R + (FN+FNF)*(RZR*C2R-RZI*C2I)
        S2I = S1I + (FN+FNF)*(RZR*C2I+RZI*C2R)
        S1R = C2R
        S1I = C2I
        FN = FN - 1.0D0
        C2R = S2R*CSR
        C2I = S2I*CSR
        CKR = C2R
        CKI = C2I
        C1R = YR(KK)
        C1I = YI(KK)
        IF (KODE.EQ.1) GO TO 280
        CALL ZS1S2(ZRR, ZRI, C1R, C1I, C2R, C2I, NW, ASC, ALIM, IUF)
        NZ = NZ + NW
  280   CONTINUE
        YR(KK) = C1R*CSPNR - C1I*CSPNI + C2R
        YI(KK) = C1R*CSPNI + C1I*CSPNR + C2I
        KK = KK - 1
        CSPNR = -CSPNR
        CSPNI = -CSPNI
        IF (IFLAG.GE.3) GO TO 290
        C2R = DABS(CKR)
        C2I = DABS(CKI)
        C2M = DMAX1(C2R,C2I)
        IF (C2M.LE.ASCLE) GO TO 290
        IFLAG = IFLAG + 1
        ASCLE = BRY(IFLAG)
        S1R = S1R*CSR
        S1I = S1I*CSR
        S2R = CKR
        S2I = CKI
        S1R = S1R*CSSR(IFLAG)
        S1I = S1I*CSSR(IFLAG)
        S2R = S2R*CSSR(IFLAG)
        S2I = S2I*CSSR(IFLAG)
        CSR = CSRR(IFLAG)
  290 CONTINUE
      RETURN
  300 CONTINUE
      NZ = -1
      RETURN
      END
      SUBROUTINE ZUNK2(ZR, ZI, FNU, KODE, MR, N, YR, YI, NZ, TOL, ELIM,
     * ALIM)
C***BEGIN PROLOGUE  ZUNK2
C***REFER TO  ZBESK
C
C     ZUNK2 COMPUTES K(FNU,Z) AND ITS ANALYTIC CONTINUATION FROM THE
C     RIGHT HALF PLANE TO THE LEFT HALF PLANE BY MEANS OF THE
C     UNIFORM ASYMPTOTIC EXPANSIONS FOR H(KIND,FNU,ZN) AND J(FNU,ZN)
C     WHERE ZN IS IN THE RIGHT HALF PLANE, KIND=(3-MR)/2, MR=+1 OR
C     -1. HERE ZN=ZR*I OR -ZR*I WHERE ZR=Z IF Z IS IN THE RIGHT
C     HALF PLANE OR ZR=-Z IF Z IS IN THE LEFT HALF PLANE. MR INDIC-
C     ATES THE DIRECTION OF ROTATION FOR ANALYTIC CONTINUATION.
C     NZ=-1 MEANS AN OVERFLOW WILL OCCUR
C
C***ROUTINES CALLED  ZAIRY,ZKSCL,ZS1S2,ZUCHK,ZUNHJ,D1MACH,ZABS
C***END PROLOGUE  ZUNK2
C     COMPLEX AI,ARG,ARGD,ASUM,ASUMD,BSUM,BSUMD,CFN,CI,CIP,CK,CONE,CRSC,
C    *CR1,CR2,CS,CSCL,CSGN,CSPN,CSR,CSS,CY,CZERO,C1,C2,DAI,PHI,PHID,RZ,
C    *S1,S2,Y,Z,ZB,ZETA1,ZETA1D,ZETA2,ZETA2D,ZN,ZR
      EXTERNAL ZABS
      DOUBLE PRECISION AARG, AIC, AII, AIR, ALIM, ANG, APHI, ARGDI,
     * ARGDR, ARGI, ARGR, ASC, ASCLE, ASUMDI, ASUMDR, ASUMI, ASUMR,
     * BRY, BSUMDI, BSUMDR, BSUMI, BSUMR, CAR, CIPI, CIPR, CKI, CKR,
     * CONER, CRSC, CR1I, CR1R, CR2I, CR2R, CSCL, CSGNI, CSI,
     * CSPNI, CSPNR, CSR, CSRR, CSSR, CYI, CYR, C1I, C1R, C2I, C2M,
     * C2R, DAII, DAIR, ELIM, FMR, FN, FNF, FNU, HPI, PHIDI, PHIDR,
     * PHII, PHIR, PI, PTI, PTR, RAST, RAZR, RS1, RZI, RZR, SAR, SGN,
     * STI, STR, S1I, S1R, S2I, S2R, TOL, YI, YR, YY, ZBI, ZBR, ZEROI,
     * ZEROR, ZETA1I, ZETA1R, ZETA2I, ZETA2R, ZET1DI, ZET1DR, ZET2DI,
     * ZET2DR, ZI, ZNI, ZNR, ZR, ZRI, ZRR, D1MACH, ZABS
      INTEGER I, IB, IFLAG, IFN, IL, IN, INU, IUF, K, KDFLG, KFLAG, KK,
     * KODE, MR, N, NAI, NDAI, NW, NZ, IDUM, J, IPARD, IC
      DIMENSION BRY(3), YR(N), YI(N), ASUMR(2), ASUMI(2), BSUMR(2),
     * BSUMI(2), PHIR(2), PHII(2), ARGR(2), ARGI(2), ZETA1R(2),
     * ZETA1I(2), ZETA2R(2), ZETA2I(2), CYR(2), CYI(2), CIPR(4),
     * CIPI(4), CSSR(3), CSRR(3)
      DATA ZEROR,ZEROI,CONER,CR1R,CR1I,CR2R,CR2I /
     1         0.0D0, 0.0D0, 1.0D0,
     1 1.0D0,1.73205080756887729D0 , -0.5D0,-8.66025403784438647D-01 /
      DATA HPI, PI, AIC /
     1     1.57079632679489662D+00,     3.14159265358979324D+00,
     1     1.26551212348464539D+00/
      DATA CIPR(1),CIPI(1),CIPR(2),CIPI(2),CIPR(3),CIPI(3),CIPR(4),
     * CIPI(4) /
     1  1.0D0,0.0D0 ,  0.0D0,-1.0D0 ,  -1.0D0,0.0D0 ,  0.0D0,1.0D0 /
C
      KDFLG = 1
      NZ = 0
C-----------------------------------------------------------------------
C     EXP(-ALIM)=EXP(-ELIM)/TOL=APPROX. ONE PRECISION GREATER THAN
C     THE UNDERFLOW LIMIT
C-----------------------------------------------------------------------
      CSCL = 1.0D0/TOL
      CRSC = TOL
      CSSR(1) = CSCL
      CSSR(2) = CONER
      CSSR(3) = CRSC
      CSRR(1) = CRSC
      CSRR(2) = CONER
      CSRR(3) = CSCL
      BRY(1) = 1.0D+3*D1MACH(1)/TOL
      BRY(2) = 1.0D0/BRY(1)
      BRY(3) = D1MACH(2)
      ZRR = ZR
      ZRI = ZI
      IF (ZR.GE.0.0D0) GO TO 10
      ZRR = -ZR
      ZRI = -ZI
   10 CONTINUE
      YY = ZRI
      ZNR = ZRI
      ZNI = -ZRR
      ZBR = ZRR
      ZBI = ZRI
      INU = INT(SNGL(FNU))
      FNF = FNU - DBLE(FLOAT(INU))
      ANG = -HPI*FNF
      CAR = DCOS(ANG)
      SAR = DSIN(ANG)
      C2R = HPI*SAR
      C2I = -HPI*CAR
      KK = MOD(INU,4) + 1
      STR = C2R*CIPR(KK) - C2I*CIPI(KK)
      STI = C2R*CIPI(KK) + C2I*CIPR(KK)
      CSR = CR1R*STR - CR1I*STI
      CSI = CR1R*STI + CR1I*STR
      IF (YY.GT.0.0D0) GO TO 20
      ZNR = -ZNR
      ZBI = -ZBI
   20 CONTINUE
C-----------------------------------------------------------------------
C     K(FNU,Z) IS COMPUTED FROM H(2,FNU,-I*Z) WHERE Z IS IN THE FIRST
C     QUADRANT. FOURTH QUADRANT VALUES (YY.LE.0.0E0) ARE COMPUTED BY
C     CONJUGATION SINCE THE K FUNCTION IS REAL ON THE POSITIVE REAL AXIS
C-----------------------------------------------------------------------
      J = 2
      DO 80 I=1,N
C-----------------------------------------------------------------------
C     J FLIP FLOPS BETWEEN 1 AND 2 IN J = 3 - J
C-----------------------------------------------------------------------
        J = 3 - J
        FN = FNU + DBLE(FLOAT(I-1))
        CALL ZUNHJ(ZNR, ZNI, FN, 0, TOL, PHIR(J), PHII(J), ARGR(J),
     *   ARGI(J), ZETA1R(J), ZETA1I(J), ZETA2R(J), ZETA2I(J), ASUMR(J),
     *   ASUMI(J), BSUMR(J), BSUMI(J))
        IF (KODE.EQ.1) GO TO 30
        STR = ZBR + ZETA2R(J)
        STI = ZBI + ZETA2I(J)
        RAST = FN/ZABS(STR,STI)
        STR = STR*RAST*RAST
        STI = -STI*RAST*RAST
        S1R = ZETA1R(J) - STR
        S1I = ZETA1I(J) - STI
        GO TO 40
   30   CONTINUE
        S1R = ZETA1R(J) - ZETA2R(J)
        S1I = ZETA1I(J) - ZETA2I(J)
   40   CONTINUE
C-----------------------------------------------------------------------
C     TEST FOR UNDERFLOW AND OVERFLOW
C-----------------------------------------------------------------------
        RS1 = S1R
        IF (DABS(RS1).GT.ELIM) GO TO 70
        IF (KDFLG.EQ.1) KFLAG = 2
        IF (DABS(RS1).LT.ALIM) GO TO 50
C-----------------------------------------------------------------------
C     REFINE  TEST AND SCALE
C-----------------------------------------------------------------------
        APHI = ZABS(PHIR(J),PHII(J))
        AARG = ZABS(ARGR(J),ARGI(J))
        RS1 = RS1 + DLOG(APHI) - 0.25D0*DLOG(AARG) - AIC
        IF (DABS(RS1).GT.ELIM) GO TO 70
        IF (KDFLG.EQ.1) KFLAG = 1
        IF (RS1.LT.0.0D0) GO TO 50
        IF (KDFLG.EQ.1) KFLAG = 3
   50   CONTINUE
C-----------------------------------------------------------------------
C     SCALE S1 TO KEEP INTERMEDIATE ARITHMETIC ON SCALE NEAR
C     EXPONENT EXTREMES
C-----------------------------------------------------------------------
        C2R = ARGR(J)*CR2R - ARGI(J)*CR2I
        C2I = ARGR(J)*CR2I + ARGI(J)*CR2R
        CALL ZAIRY(C2R, C2I, 0, 2, AIR, AII, NAI, IDUM)
        CALL ZAIRY(C2R, C2I, 1, 2, DAIR, DAII, NDAI, IDUM)
        STR = DAIR*BSUMR(J) - DAII*BSUMI(J)
        STI = DAIR*BSUMI(J) + DAII*BSUMR(J)
        PTR = STR*CR2R - STI*CR2I
        PTI = STR*CR2I + STI*CR2R
        STR = PTR + (AIR*ASUMR(J)-AII*ASUMI(J))
        STI = PTI + (AIR*ASUMI(J)+AII*ASUMR(J))
        PTR = STR*PHIR(J) - STI*PHII(J)
        PTI = STR*PHII(J) + STI*PHIR(J)
        S2R = PTR*CSR - PTI*CSI
        S2I = PTR*CSI + PTI*CSR
        STR = DEXP(S1R)*CSSR(KFLAG)
        S1R = STR*DCOS(S1I)
        S1I = STR*DSIN(S1I)
        STR = S2R*S1R - S2I*S1I
        S2I = S1R*S2I + S2R*S1I
        S2R = STR
        IF (KFLAG.NE.1) GO TO 60
        CALL ZUCHK(S2R, S2I, NW, BRY(1), TOL)
        IF (NW.NE.0) GO TO 70
   60   CONTINUE
        IF (YY.LE.0.0D0) S2I = -S2I
        CYR(KDFLG) = S2R
        CYI(KDFLG) = S2I
        YR(I) = S2R*CSRR(KFLAG)
        YI(I) = S2I*CSRR(KFLAG)
        STR = CSI
        CSI = -CSR
        CSR = STR
        IF (KDFLG.EQ.2) GO TO 85
        KDFLG = 2
        GO TO 80
   70   CONTINUE
        IF (RS1.GT.0.0D0) GO TO 320
C-----------------------------------------------------------------------
C     FOR ZR.LT.0.0, THE I FUNCTION TO BE ADDED WILL OVERFLOW
C-----------------------------------------------------------------------
        IF (ZR.LT.0.0D0) GO TO 320
        KDFLG = 1
        YR(I)=ZEROR
        YI(I)=ZEROI
        NZ=NZ+1
        STR = CSI
        CSI =-CSR
        CSR = STR
        IF (I.EQ.1) GO TO 80
        IF ((YR(I-1).EQ.ZEROR).AND.(YI(I-1).EQ.ZEROI)) GO TO 80
        YR(I-1)=ZEROR
        YI(I-1)=ZEROI
        NZ=NZ+1
   80 CONTINUE
      I = N
   85 CONTINUE
      RAZR = 1.0D0/ZABS(ZRR,ZRI)
      STR = ZRR*RAZR
      STI = -ZRI*RAZR
      RZR = (STR+STR)*RAZR
      RZI = (STI+STI)*RAZR
      CKR = FN*RZR
      CKI = FN*RZI
      IB = I + 1
      IF (N.LT.IB) GO TO 180
C-----------------------------------------------------------------------
C     TEST LAST MEMBER FOR UNDERFLOW AND OVERFLOW. SET SEQUENCE TO ZERO
C     ON UNDERFLOW.
C-----------------------------------------------------------------------
      FN = FNU + DBLE(FLOAT(N-1))
      IPARD = 1
      IF (MR.NE.0) IPARD = 0
      CALL ZUNHJ(ZNR, ZNI, FN, IPARD, TOL, PHIDR, PHIDI, ARGDR, ARGDI,
     * ZET1DR, ZET1DI, ZET2DR, ZET2DI, ASUMDR, ASUMDI, BSUMDR, BSUMDI)
      IF (KODE.EQ.1) GO TO 90
      STR = ZBR + ZET2DR
      STI = ZBI + ZET2DI
      RAST = FN/ZABS(STR,STI)
      STR = STR*RAST*RAST
      STI = -STI*RAST*RAST
      S1R = ZET1DR - STR
      S1I = ZET1DI - STI
      GO TO 100
   90 CONTINUE
      S1R = ZET1DR - ZET2DR
      S1I = ZET1DI - ZET2DI
  100 CONTINUE
      RS1 = S1R
      IF (DABS(RS1).GT.ELIM) GO TO 105
      IF (DABS(RS1).LT.ALIM) GO TO 120
C----------------------------------------------------------------------------
C     REFINE ESTIMATE AND TEST
C-------------------------------------------------------------------------
      APHI = ZABS(PHIDR,PHIDI)
      RS1 = RS1+DLOG(APHI)
      IF (DABS(RS1).LT.ELIM) GO TO 120
  105 CONTINUE
      IF (RS1.GT.0.0D0) GO TO 320
C-----------------------------------------------------------------------
C     FOR ZR.LT.0.0, THE I FUNCTION TO BE ADDED WILL OVERFLOW
C-----------------------------------------------------------------------
      IF (ZR.LT.0.0D0) GO TO 320
      NZ = N
      DO 106 I=1,N
        YR(I) = ZEROR
        YI(I) = ZEROI
  106 CONTINUE
      RETURN
  120 CONTINUE
      S1R = CYR(1)
      S1I = CYI(1)
      S2R = CYR(2)
      S2I = CYI(2)
      C1R = CSRR(KFLAG)
      ASCLE = BRY(KFLAG)
      DO 130 I=IB,N
        C2R = S2R
        C2I = S2I
        S2R = CKR*C2R - CKI*C2I + S1R
        S2I = CKR*C2I + CKI*C2R + S1I
        S1R = C2R
        S1I = C2I
        CKR = CKR + RZR
        CKI = CKI + RZI
        C2R = S2R*C1R
        C2I = S2I*C1R
        YR(I) = C2R
        YI(I) = C2I
        IF (KFLAG.GE.3) GO TO 130
        STR = DABS(C2R)
        STI = DABS(C2I)
        C2M = DMAX1(STR,STI)
        IF (C2M.LE.ASCLE) GO TO 130
        KFLAG = KFLAG + 1
        ASCLE = BRY(KFLAG)
        S1R = S1R*C1R
        S1I = S1I*C1R
        S2R = C2R
        S2I = C2I
        S1R = S1R*CSSR(KFLAG)
        S1I = S1I*CSSR(KFLAG)
        S2R = S2R*CSSR(KFLAG)
        S2I = S2I*CSSR(KFLAG)
        C1R = CSRR(KFLAG)
  130 CONTINUE
  180 CONTINUE
      IF (MR.EQ.0) RETURN
C-----------------------------------------------------------------------
C     ANALYTIC CONTINUATION FOR RE(Z).LT.0.0D0
C-----------------------------------------------------------------------
      NZ = 0
      FMR = DBLE(FLOAT(MR))
      SGN = -DSIGN(PI,FMR)
C-----------------------------------------------------------------------
C     CSPN AND CSGN ARE COEFF OF K AND I FUNCIONS RESP.
C-----------------------------------------------------------------------
      CSGNI = SGN
      IF (YY.LE.0.0D0) CSGNI = -CSGNI
      IFN = INU + N - 1
      ANG = FNF*SGN
      CSPNR = DCOS(ANG)
      CSPNI = DSIN(ANG)
      IF (MOD(IFN,2).EQ.0) GO TO 190
      CSPNR = -CSPNR
      CSPNI = -CSPNI
  190 CONTINUE
C-----------------------------------------------------------------------
C     CS=COEFF OF THE J FUNCTION TO GET THE I FUNCTION. I(FNU,Z) IS
C     COMPUTED FROM EXP(I*FNU*HPI)*J(FNU,-I*Z) WHERE Z IS IN THE FIRST
C     QUADRANT. FOURTH QUADRANT VALUES (YY.LE.0.0E0) ARE COMPUTED BY
C     CONJUGATION SINCE THE I FUNCTION IS REAL ON THE POSITIVE REAL AXIS
C-----------------------------------------------------------------------
      CSR = SAR*CSGNI
      CSI = CAR*CSGNI
      IN = MOD(IFN,4) + 1
      C2R = CIPR(IN)
      C2I = CIPI(IN)
      STR = CSR*C2R + CSI*C2I
      CSI = -CSR*C2I + CSI*C2R
      CSR = STR
      ASC = BRY(1)
      IUF = 0
      KK = N
      KDFLG = 1
      IB = IB - 1
      IC = IB - 1
      DO 290 K=1,N
        FN = FNU + DBLE(FLOAT(KK-1))
C-----------------------------------------------------------------------
C     LOGIC TO SORT OUT CASES WHOSE PARAMETERS WERE SET FOR THE K
C     FUNCTION ABOVE
C-----------------------------------------------------------------------
        IF (N.GT.2) GO TO 175
  172   CONTINUE
        PHIDR = PHIR(J)
        PHIDI = PHII(J)
        ARGDR = ARGR(J)
        ARGDI = ARGI(J)
        ZET1DR = ZETA1R(J)
        ZET1DI = ZETA1I(J)
        ZET2DR = ZETA2R(J)
        ZET2DI = ZETA2I(J)
        ASUMDR = ASUMR(J)
        ASUMDI = ASUMI(J)
        BSUMDR = BSUMR(J)
        BSUMDI = BSUMI(J)
        J = 3 - J
        GO TO 210
  175   CONTINUE
        IF ((KK.EQ.N).AND.(IB.LT.N)) GO TO 210
        IF ((KK.EQ.IB).OR.(KK.EQ.IC)) GO TO 172
        CALL ZUNHJ(ZNR, ZNI, FN, 0, TOL, PHIDR, PHIDI, ARGDR,
     *   ARGDI, ZET1DR, ZET1DI, ZET2DR, ZET2DI, ASUMDR,
     *   ASUMDI, BSUMDR, BSUMDI)
  210   CONTINUE
        IF (KODE.EQ.1) GO TO 220
        STR = ZBR + ZET2DR
        STI = ZBI + ZET2DI
        RAST = FN/ZABS(STR,STI)
        STR = STR*RAST*RAST
        STI = -STI*RAST*RAST
        S1R = -ZET1DR + STR
        S1I = -ZET1DI + STI
        GO TO 230
  220   CONTINUE
        S1R = -ZET1DR + ZET2DR
        S1I = -ZET1DI + ZET2DI
  230   CONTINUE
C-----------------------------------------------------------------------
C     TEST FOR UNDERFLOW AND OVERFLOW
C-----------------------------------------------------------------------
        RS1 = S1R
        IF (DABS(RS1).GT.ELIM) GO TO 280
        IF (KDFLG.EQ.1) IFLAG = 2
        IF (DABS(RS1).LT.ALIM) GO TO 240
C-----------------------------------------------------------------------
C     REFINE  TEST AND SCALE
C-----------------------------------------------------------------------
        APHI = ZABS(PHIDR,PHIDI)
        AARG = ZABS(ARGDR,ARGDI)
        RS1 = RS1 + DLOG(APHI) - 0.25D0*DLOG(AARG) - AIC
        IF (DABS(RS1).GT.ELIM) GO TO 280
        IF (KDFLG.EQ.1) IFLAG = 1
        IF (RS1.LT.0.0D0) GO TO 240
        IF (KDFLG.EQ.1) IFLAG = 3
  240   CONTINUE
        CALL ZAIRY(ARGDR, ARGDI, 0, 2, AIR, AII, NAI, IDUM)
        CALL ZAIRY(ARGDR, ARGDI, 1, 2, DAIR, DAII, NDAI, IDUM)
        STR = DAIR*BSUMDR - DAII*BSUMDI
        STI = DAIR*BSUMDI + DAII*BSUMDR
        STR = STR + (AIR*ASUMDR-AII*ASUMDI)
        STI = STI + (AIR*ASUMDI+AII*ASUMDR)
        PTR = STR*PHIDR - STI*PHIDI
        PTI = STR*PHIDI + STI*PHIDR
        S2R = PTR*CSR - PTI*CSI
        S2I = PTR*CSI + PTI*CSR
        STR = DEXP(S1R)*CSSR(IFLAG)
        S1R = STR*DCOS(S1I)
        S1I = STR*DSIN(S1I)
        STR = S2R*S1R - S2I*S1I
        S2I = S2R*S1I + S2I*S1R
        S2R = STR
        IF (IFLAG.NE.1) GO TO 250
        CALL ZUCHK(S2R, S2I, NW, BRY(1), TOL)
        IF (NW.EQ.0) GO TO 250
        S2R = ZEROR
        S2I = ZEROI
  250   CONTINUE
        IF (YY.LE.0.0D0) S2I = -S2I
        CYR(KDFLG) = S2R
        CYI(KDFLG) = S2I
        C2R = S2R
        C2I = S2I
        S2R = S2R*CSRR(IFLAG)
        S2I = S2I*CSRR(IFLAG)
C-----------------------------------------------------------------------
C     ADD I AND K FUNCTIONS, K SEQUENCE IN Y(I), I=1,N
C-----------------------------------------------------------------------
        S1R = YR(KK)
        S1I = YI(KK)
        IF (KODE.EQ.1) GO TO 270
        CALL ZS1S2(ZRR, ZRI, S1R, S1I, S2R, S2I, NW, ASC, ALIM, IUF)
        NZ = NZ + NW
  270   CONTINUE
        YR(KK) = S1R*CSPNR - S1I*CSPNI + S2R
        YI(KK) = S1R*CSPNI + S1I*CSPNR + S2I
        KK = KK - 1
        CSPNR = -CSPNR
        CSPNI = -CSPNI
        STR = CSI
        CSI = -CSR
        CSR = STR
        IF (C2R.NE.0.0D0 .OR. C2I.NE.0.0D0) GO TO 255
        KDFLG = 1
        GO TO 290
  255   CONTINUE
        IF (KDFLG.EQ.2) GO TO 295
        KDFLG = 2
        GO TO 290
  280   CONTINUE
        IF (RS1.GT.0.0D0) GO TO 320
        S2R = ZEROR
        S2I = ZEROI
        GO TO 250
  290 CONTINUE
      K = N
  295 CONTINUE
      IL = N - K
      IF (IL.EQ.0) RETURN
C-----------------------------------------------------------------------
C     RECUR BACKWARD FOR REMAINDER OF I SEQUENCE AND ADD IN THE
C     K FUNCTIONS, SCALING THE I SEQUENCE DURING RECURRENCE TO KEEP
C     INTERMEDIATE ARITHMETIC ON SCALE NEAR EXPONENT EXTREMES.
C-----------------------------------------------------------------------
      S1R = CYR(1)
      S1I = CYI(1)
      S2R = CYR(2)
      S2I = CYI(2)
      CSR = CSRR(IFLAG)
      ASCLE = BRY(IFLAG)
      FN = DBLE(FLOAT(INU+IL))
      DO 310 I=1,IL
        C2R = S2R
        C2I = S2I
        S2R = S1R + (FN+FNF)*(RZR*C2R-RZI*C2I)
        S2I = S1I + (FN+FNF)*(RZR*C2I+RZI*C2R)
        S1R = C2R
        S1I = C2I
        FN = FN - 1.0D0
        C2R = S2R*CSR
        C2I = S2I*CSR
        CKR = C2R
        CKI = C2I
        C1R = YR(KK)
        C1I = YI(KK)
        IF (KODE.EQ.1) GO TO 300
        CALL ZS1S2(ZRR, ZRI, C1R, C1I, C2R, C2I, NW, ASC, ALIM, IUF)
        NZ = NZ + NW
  300   CONTINUE
        YR(KK) = C1R*CSPNR - C1I*CSPNI + C2R
        YI(KK) = C1R*CSPNI + C1I*CSPNR + C2I
        KK = KK - 1
        CSPNR = -CSPNR
        CSPNI = -CSPNI
        IF (IFLAG.GE.3) GO TO 310
        C2R = DABS(CKR)
        C2I = DABS(CKI)
        C2M = DMAX1(C2R,C2I)
        IF (C2M.LE.ASCLE) GO TO 310
        IFLAG = IFLAG + 1
        ASCLE = BRY(IFLAG)
        S1R = S1R*CSR
        S1I = S1I*CSR
        S2R = CKR
        S2I = CKI
        S1R = S1R*CSSR(IFLAG)
        S1I = S1I*CSSR(IFLAG)
        S2R = S2R*CSSR(IFLAG)
        S2I = S2I*CSSR(IFLAG)
        CSR = CSRR(IFLAG)
  310 CONTINUE
      RETURN
  320 CONTINUE
      NZ = -1
      RETURN
      END
      SUBROUTINE ZBUNI(ZR, ZI, FNU, KODE, N, YR, YI, NZ, NUI, NLAST,
     * FNUL, TOL, ELIM, ALIM)
C***BEGIN PROLOGUE  ZBUNI
C***REFER TO  ZBESI,ZBESK
C
C     ZBUNI COMPUTES THE I BESSEL FUNCTION FOR LARGE CABS(Z).GT.
C     FNUL AND FNU+N-1.LT.FNUL. THE ORDER IS INCREASED FROM
C     FNU+N-1 GREATER THAN FNUL BY ADDING NUI AND COMPUTING
C     ACCORDING TO THE UNIFORM ASYMPTOTIC EXPANSION FOR I(FNU,Z)
C     ON IFORM=1 AND THE EXPANSION FOR J(FNU,Z) ON IFORM=2
C
C***ROUTINES CALLED  ZUNI1,ZUNI2,ZABS,D1MACH
C***END PROLOGUE  ZBUNI
C     COMPLEX CSCL,CSCR,CY,RZ,ST,S1,S2,Y,Z
      EXTERNAL ZABS
      DOUBLE PRECISION ALIM, AX, AY, CSCLR, CSCRR, CYI, CYR, DFNU,
     * ELIM, FNU, FNUI, FNUL, GNU, RAZ, RZI, RZR, STI, STR, S1I, S1R,
     * S2I, S2R, TOL, YI, YR, ZI, ZR, ZABS, ASCLE, BRY, C1R, C1I, C1M,
     * D1MACH
      INTEGER I, IFLAG, IFORM, K, KODE, N, NL, NLAST, NUI, NW, NZ
      DIMENSION YR(N), YI(N), CYR(2), CYI(2), BRY(3)
      NZ = 0
      AX = DABS(ZR)*1.7321D0
      AY = DABS(ZI)
      IFORM = 1
      IF (AY.GT.AX) IFORM = 2
      IF (NUI.EQ.0) GO TO 60
      FNUI = DBLE(FLOAT(NUI))
      DFNU = FNU + DBLE(FLOAT(N-1))
      GNU = DFNU + FNUI
      IF (IFORM.EQ.2) GO TO 10
C-----------------------------------------------------------------------
C     ASYMPTOTIC EXPANSION FOR I(FNU,Z) FOR LARGE FNU APPLIED IN
C     -PI/3.LE.ARG(Z).LE.PI/3
C-----------------------------------------------------------------------
      CALL ZUNI1(ZR, ZI, GNU, KODE, 2, CYR, CYI, NW, NLAST, FNUL, TOL,
     * ELIM, ALIM)
      GO TO 20
   10 CONTINUE
C-----------------------------------------------------------------------
C     ASYMPTOTIC EXPANSION FOR J(FNU,Z*EXP(M*HPI)) FOR LARGE FNU
C     APPLIED IN PI/3.LT.ABS(ARG(Z)).LE.PI/2 WHERE M=+I OR -I
C     AND HPI=PI/2
C-----------------------------------------------------------------------
      CALL ZUNI2(ZR, ZI, GNU, KODE, 2, CYR, CYI, NW, NLAST, FNUL, TOL,
     * ELIM, ALIM)
   20 CONTINUE
      IF (NW.LT.0) GO TO 50
      IF (NW.NE.0) GO TO 90
      STR = ZABS(CYR(1),CYI(1))
C----------------------------------------------------------------------
C     SCALE BACKWARD RECURRENCE, BRY(3) IS DEFINED BUT NEVER USED
C----------------------------------------------------------------------
      BRY(1)=1.0D+3*D1MACH(1)/TOL
      BRY(2) = 1.0D0/BRY(1)
      BRY(3) = BRY(2)
      IFLAG = 2
      ASCLE = BRY(2)
      CSCLR = 1.0D0
      IF (STR.GT.BRY(1)) GO TO 21
      IFLAG = 1
      ASCLE = BRY(1)
      CSCLR = 1.0D0/TOL
      GO TO 25
   21 CONTINUE
      IF (STR.LT.BRY(2)) GO TO 25
      IFLAG = 3
      ASCLE=BRY(3)
      CSCLR = TOL
   25 CONTINUE
      CSCRR = 1.0D0/CSCLR
      S1R = CYR(2)*CSCLR
      S1I = CYI(2)*CSCLR
      S2R = CYR(1)*CSCLR
      S2I = CYI(1)*CSCLR
      RAZ = 1.0D0/ZABS(ZR,ZI)
      STR = ZR*RAZ
      STI = -ZI*RAZ
      RZR = (STR+STR)*RAZ
      RZI = (STI+STI)*RAZ
      DO 30 I=1,NUI
        STR = S2R
        STI = S2I
        S2R = (DFNU+FNUI)*(RZR*STR-RZI*STI) + S1R
        S2I = (DFNU+FNUI)*(RZR*STI+RZI*STR) + S1I
        S1R = STR
        S1I = STI
        FNUI = FNUI - 1.0D0
        IF (IFLAG.GE.3) GO TO 30
        STR = S2R*CSCRR
        STI = S2I*CSCRR
        C1R = DABS(STR)
        C1I = DABS(STI)
        C1M = DMAX1(C1R,C1I)
        IF (C1M.LE.ASCLE) GO TO 30
        IFLAG = IFLAG+1
        ASCLE = BRY(IFLAG)
        S1R = S1R*CSCRR
        S1I = S1I*CSCRR
        S2R = STR
        S2I = STI
        CSCLR = CSCLR*TOL
        CSCRR = 1.0D0/CSCLR
        S1R = S1R*CSCLR
        S1I = S1I*CSCLR
        S2R = S2R*CSCLR
        S2I = S2I*CSCLR
   30 CONTINUE
      YR(N) = S2R*CSCRR
      YI(N) = S2I*CSCRR
      IF (N.EQ.1) RETURN
      NL = N - 1
      FNUI = DBLE(FLOAT(NL))
      K = NL
      DO 40 I=1,NL
        STR = S2R
        STI = S2I
        S2R = (FNU+FNUI)*(RZR*STR-RZI*STI) + S1R
        S2I = (FNU+FNUI)*(RZR*STI+RZI*STR) + S1I
        S1R = STR
        S1I = STI
        STR = S2R*CSCRR
        STI = S2I*CSCRR
        YR(K) = STR
        YI(K) = STI
        FNUI = FNUI - 1.0D0
        K = K - 1
        IF (IFLAG.GE.3) GO TO 40
        C1R = DABS(STR)
        C1I = DABS(STI)
        C1M = DMAX1(C1R,C1I)
        IF (C1M.LE.ASCLE) GO TO 40
        IFLAG = IFLAG+1
        ASCLE = BRY(IFLAG)
        S1R = S1R*CSCRR
        S1I = S1I*CSCRR
        S2R = STR
        S2I = STI
        CSCLR = CSCLR*TOL
        CSCRR = 1.0D0/CSCLR
        S1R = S1R*CSCLR
        S1I = S1I*CSCLR
        S2R = S2R*CSCLR
        S2I = S2I*CSCLR
   40 CONTINUE
      RETURN
   50 CONTINUE
      NZ = -1
      IF(NW.EQ.(-2)) NZ=-2
      RETURN
   60 CONTINUE
      IF (IFORM.EQ.2) GO TO 70
C-----------------------------------------------------------------------
C     ASYMPTOTIC EXPANSION FOR I(FNU,Z) FOR LARGE FNU APPLIED IN
C     -PI/3.LE.ARG(Z).LE.PI/3
C-----------------------------------------------------------------------
      CALL ZUNI1(ZR, ZI, FNU, KODE, N, YR, YI, NW, NLAST, FNUL, TOL,
     * ELIM, ALIM)
      GO TO 80
   70 CONTINUE
C-----------------------------------------------------------------------
C     ASYMPTOTIC EXPANSION FOR J(FNU,Z*EXP(M*HPI)) FOR LARGE FNU
C     APPLIED IN PI/3.LT.ABS(ARG(Z)).LE.PI/2 WHERE M=+I OR -I
C     AND HPI=PI/2
C-----------------------------------------------------------------------
      CALL ZUNI2(ZR, ZI, FNU, KODE, N, YR, YI, NW, NLAST, FNUL, TOL,
     * ELIM, ALIM)
   80 CONTINUE
      IF (NW.LT.0) GO TO 50
      NZ = NW
      RETURN
   90 CONTINUE
      NLAST = N
      RETURN
      END
      SUBROUTINE ZUNI1(ZR, ZI, FNU, KODE, N, YR, YI, NZ, NLAST, FNUL,
     * TOL, ELIM, ALIM)
C***BEGIN PROLOGUE  ZUNI1
C***REFER TO  ZBESI,ZBESK
C
C     ZUNI1 COMPUTES I(FNU,Z)  BY MEANS OF THE UNIFORM ASYMPTOTIC
C     EXPANSION FOR I(FNU,Z) IN -PI/3.LE.ARG Z.LE.PI/3.
C
C     FNUL IS THE SMALLEST ORDER PERMITTED FOR THE ASYMPTOTIC
C     EXPANSION. NLAST=0 MEANS ALL OF THE Y VALUES WERE SET.
C     NLAST.NE.0 IS THE NUMBER LEFT TO BE COMPUTED BY ANOTHER
C     FORMULA FOR ORDERS FNU TO FNU+NLAST-1 BECAUSE FNU+NLAST-1.LT.FNUL.
C     Y(I)=CZERO FOR I=NLAST+1,N
C
C***ROUTINES CALLED  ZUCHK,ZUNIK,ZUOIK,D1MACH,ZABS
C***END PROLOGUE  ZUNI1
C     COMPLEX CFN,CONE,CRSC,CSCL,CSR,CSS,CWRK,CZERO,C1,C2,PHI,RZ,SUM,S1,
C    *S2,Y,Z,ZETA1,ZETA2
      EXTERNAL ZABS
      DOUBLE PRECISION ALIM, APHI, ASCLE, BRY, CONER, CRSC,
     * CSCL, CSRR, CSSR, CWRKI, CWRKR, C1R, C2I, C2M, C2R, ELIM, FN,
     * FNU, FNUL, PHII, PHIR, RAST, RS1, RZI, RZR, STI, STR, SUMI,
     * SUMR, S1I, S1R, S2I, S2R, TOL, YI, YR, ZEROI, ZEROR, ZETA1I,
     * ZETA1R, ZETA2I, ZETA2R, ZI, ZR, CYR, CYI, D1MACH, ZABS
      INTEGER I, IFLAG, INIT, K, KODE, M, N, ND, NLAST, NN, NUF, NW, NZ
      DIMENSION BRY(3), YR(N), YI(N), CWRKR(16), CWRKI(16), CSSR(3),
     * CSRR(3), CYR(2), CYI(2)
      DATA ZEROR,ZEROI,CONER / 0.0D0, 0.0D0, 1.0D0 /
C
      NZ = 0
      ND = N
      NLAST = 0
C-----------------------------------------------------------------------
C     COMPUTED VALUES WITH EXPONENTS BETWEEN ALIM AND ELIM IN MAG-
C     NITUDE ARE SCALED TO KEEP INTERMEDIATE ARITHMETIC ON SCALE,
C     EXP(ALIM)=EXP(ELIM)*TOL
C-----------------------------------------------------------------------
      CSCL = 1.0D0/TOL
      CRSC = TOL
      CSSR(1) = CSCL
      CSSR(2) = CONER
      CSSR(3) = CRSC
      CSRR(1) = CRSC
      CSRR(2) = CONER
      CSRR(3) = CSCL
      BRY(1) = 1.0D+3*D1MACH(1)/TOL
C-----------------------------------------------------------------------
C     CHECK FOR UNDERFLOW AND OVERFLOW ON FIRST MEMBER
C-----------------------------------------------------------------------
      FN = DMAX1(FNU,1.0D0)
      INIT = 0
      CALL ZUNIK(ZR, ZI, FN, 1, 1, TOL, INIT, PHIR, PHII, ZETA1R,
     * ZETA1I, ZETA2R, ZETA2I, SUMR, SUMI, CWRKR, CWRKI)
      IF (KODE.EQ.1) GO TO 10
      STR = ZR + ZETA2R
      STI = ZI + ZETA2I
      RAST = FN/ZABS(STR,STI)
      STR = STR*RAST*RAST
      STI = -STI*RAST*RAST
      S1R = -ZETA1R + STR
      S1I = -ZETA1I + STI
      GO TO 20
   10 CONTINUE
      S1R = -ZETA1R + ZETA2R
      S1I = -ZETA1I + ZETA2I
   20 CONTINUE
      RS1 = S1R
      IF (DABS(RS1).GT.ELIM) GO TO 130
   30 CONTINUE
      NN = MIN0(2,ND)
      DO 80 I=1,NN
        FN = FNU + DBLE(FLOAT(ND-I))
        INIT = 0
        CALL ZUNIK(ZR, ZI, FN, 1, 0, TOL, INIT, PHIR, PHII, ZETA1R,
     *   ZETA1I, ZETA2R, ZETA2I, SUMR, SUMI, CWRKR, CWRKI)
        IF (KODE.EQ.1) GO TO 40
        STR = ZR + ZETA2R
        STI = ZI + ZETA2I
        RAST = FN/ZABS(STR,STI)
        STR = STR*RAST*RAST
        STI = -STI*RAST*RAST
        S1R = -ZETA1R + STR
        S1I = -ZETA1I + STI + ZI
        GO TO 50
   40   CONTINUE
        S1R = -ZETA1R + ZETA2R
        S1I = -ZETA1I + ZETA2I
   50   CONTINUE
C-----------------------------------------------------------------------
C     TEST FOR UNDERFLOW AND OVERFLOW
C-----------------------------------------------------------------------
        RS1 = S1R
        IF (DABS(RS1).GT.ELIM) GO TO 110
        IF (I.EQ.1) IFLAG = 2
        IF (DABS(RS1).LT.ALIM) GO TO 60
C-----------------------------------------------------------------------
C     REFINE  TEST AND SCALE
C-----------------------------------------------------------------------
        APHI = ZABS(PHIR,PHII)
        RS1 = RS1 + DLOG(APHI)
        IF (DABS(RS1).GT.ELIM) GO TO 110
        IF (I.EQ.1) IFLAG = 1
        IF (RS1.LT.0.0D0) GO TO 60
        IF (I.EQ.1) IFLAG = 3
   60   CONTINUE
C-----------------------------------------------------------------------
C     SCALE S1 IF CABS(S1).LT.ASCLE
C-----------------------------------------------------------------------
        S2R = PHIR*SUMR - PHII*SUMI
        S2I = PHIR*SUMI + PHII*SUMR
        STR = DEXP(S1R)*CSSR(IFLAG)
        S1R = STR*DCOS(S1I)
        S1I = STR*DSIN(S1I)
        STR = S2R*S1R - S2I*S1I
        S2I = S2R*S1I + S2I*S1R
        S2R = STR
        IF (IFLAG.NE.1) GO TO 70
        CALL ZUCHK(S2R, S2I, NW, BRY(1), TOL)
        IF (NW.NE.0) GO TO 110
   70   CONTINUE
        CYR(I) = S2R
        CYI(I) = S2I
        M = ND - I + 1
        YR(M) = S2R*CSRR(IFLAG)
        YI(M) = S2I*CSRR(IFLAG)
   80 CONTINUE
      IF (ND.LE.2) GO TO 100
      RAST = 1.0D0/ZABS(ZR,ZI)
      STR = ZR*RAST
      STI = -ZI*RAST
      RZR = (STR+STR)*RAST
      RZI = (STI+STI)*RAST
      BRY(2) = 1.0D0/BRY(1)
      BRY(3) = D1MACH(2)
      S1R = CYR(1)
      S1I = CYI(1)
      S2R = CYR(2)
      S2I = CYI(2)
      C1R = CSRR(IFLAG)
      ASCLE = BRY(IFLAG)
      K = ND - 2
      FN = DBLE(FLOAT(K))
      DO 90 I=3,ND
        C2R = S2R
        C2I = S2I
        S2R = S1R + (FNU+FN)*(RZR*C2R-RZI*C2I)
        S2I = S1I + (FNU+FN)*(RZR*C2I+RZI*C2R)
        S1R = C2R
        S1I = C2I
        C2R = S2R*C1R
        C2I = S2I*C1R
        YR(K) = C2R
        YI(K) = C2I
        K = K - 1
        FN = FN - 1.0D0
        IF (IFLAG.GE.3) GO TO 90
        STR = DABS(C2R)
        STI = DABS(C2I)
        C2M = DMAX1(STR,STI)
        IF (C2M.LE.ASCLE) GO TO 90
        IFLAG = IFLAG + 1
        ASCLE = BRY(IFLAG)
        S1R = S1R*C1R
        S1I = S1I*C1R
        S2R = C2R
        S2I = C2I
        S1R = S1R*CSSR(IFLAG)
        S1I = S1I*CSSR(IFLAG)
        S2R = S2R*CSSR(IFLAG)
        S2I = S2I*CSSR(IFLAG)
        C1R = CSRR(IFLAG)
   90 CONTINUE
  100 CONTINUE
      RETURN
C-----------------------------------------------------------------------
C     SET UNDERFLOW AND UPDATE PARAMETERS
C-----------------------------------------------------------------------
  110 CONTINUE
      IF (RS1.GT.0.0D0) GO TO 120
      YR(ND) = ZEROR
      YI(ND) = ZEROI
      NZ = NZ + 1
      ND = ND - 1
      IF (ND.EQ.0) GO TO 100
      CALL ZUOIK(ZR, ZI, FNU, KODE, 1, ND, YR, YI, NUF, TOL, ELIM, ALIM)
      IF (NUF.LT.0) GO TO 120
      ND = ND - NUF
      NZ = NZ + NUF
      IF (ND.EQ.0) GO TO 100
      FN = FNU + DBLE(FLOAT(ND-1))
      IF (FN.GE.FNUL) GO TO 30
      NLAST = ND
      RETURN
  120 CONTINUE
      NZ = -1
      RETURN
  130 CONTINUE
      IF (RS1.GT.0.0D0) GO TO 120
      NZ = N
      DO 140 I=1,N
        YR(I) = ZEROR
        YI(I) = ZEROI
  140 CONTINUE
      RETURN
      END
      SUBROUTINE ZUNI2(ZR, ZI, FNU, KODE, N, YR, YI, NZ, NLAST, FNUL,
     * TOL, ELIM, ALIM)
C***BEGIN PROLOGUE  ZUNI2
C***REFER TO  ZBESI,ZBESK
C
C     ZUNI2 COMPUTES I(FNU,Z) IN THE RIGHT HALF PLANE BY MEANS OF
C     UNIFORM ASYMPTOTIC EXPANSION FOR J(FNU,ZN) WHERE ZN IS Z*I
C     OR -Z*I AND ZN IS IN THE RIGHT HALF PLANE ALSO.
C
C     FNUL IS THE SMALLEST ORDER PERMITTED FOR THE ASYMPTOTIC
C     EXPANSION. NLAST=0 MEANS ALL OF THE Y VALUES WERE SET.
C     NLAST.NE.0 IS THE NUMBER LEFT TO BE COMPUTED BY ANOTHER
C     FORMULA FOR ORDERS FNU TO FNU+NLAST-1 BECAUSE FNU+NLAST-1.LT.FNUL.
C     Y(I)=CZERO FOR I=NLAST+1,N
C
C***ROUTINES CALLED  ZAIRY,ZUCHK,ZUNHJ,ZUOIK,D1MACH,ZABS
C***END PROLOGUE  ZUNI2
C     COMPLEX AI,ARG,ASUM,BSUM,CFN,CI,CID,CIP,CONE,CRSC,CSCL,CSR,CSS,
C    *CZERO,C1,C2,DAI,PHI,RZ,S1,S2,Y,Z,ZB,ZETA1,ZETA2,ZN
      EXTERNAL ZABS
      DOUBLE PRECISION AARG, AIC, AII, AIR, ALIM, ANG, APHI, ARGI,
     * ARGR, ASCLE, ASUMI, ASUMR, BRY, BSUMI, BSUMR, CIDI, CIPI, CIPR,
     * CONER, CRSC, CSCL, CSRR, CSSR, C1R, C2I, C2M, C2R, DAII,
     * DAIR, ELIM, FN, FNU, FNUL, HPI, PHII, PHIR, RAST, RAZ, RS1, RZI,
     * RZR, STI, STR, S1I, S1R, S2I, S2R, TOL, YI, YR, ZBI, ZBR, ZEROI,
     * ZEROR, ZETA1I, ZETA1R, ZETA2I, ZETA2R, ZI, ZNI, ZNR, ZR, CYR,
     * CYI, D1MACH, ZABS, CAR, SAR
      INTEGER I, IFLAG, IN, INU, J, K, KODE, N, NAI, ND, NDAI, NLAST,
     * NN, NUF, NW, NZ, IDUM
      DIMENSION BRY(3), YR(N), YI(N), CIPR(4), CIPI(4), CSSR(3),
     * CSRR(3), CYR(2), CYI(2)
      DATA ZEROR,ZEROI,CONER / 0.0D0, 0.0D0, 1.0D0 /
      DATA CIPR(1),CIPI(1),CIPR(2),CIPI(2),CIPR(3),CIPI(3),CIPR(4),
     * CIPI(4)/ 1.0D0,0.0D0, 0.0D0,1.0D0, -1.0D0,0.0D0, 0.0D0,-1.0D0/
      DATA HPI, AIC  /
     1      1.57079632679489662D+00,     1.265512123484645396D+00/
C
      NZ = 0
      ND = N
      NLAST = 0
C-----------------------------------------------------------------------
C     COMPUTED VALUES WITH EXPONENTS BETWEEN ALIM AND ELIM IN MAG-
C     NITUDE ARE SCALED TO KEEP INTERMEDIATE ARITHMETIC ON SCALE,
C     EXP(ALIM)=EXP(ELIM)*TOL
C-----------------------------------------------------------------------
      CSCL = 1.0D0/TOL
      CRSC = TOL
      CSSR(1) = CSCL
      CSSR(2) = CONER
      CSSR(3) = CRSC
      CSRR(1) = CRSC
      CSRR(2) = CONER
      CSRR(3) = CSCL
      BRY(1) = 1.0D+3*D1MACH(1)/TOL
C-----------------------------------------------------------------------
C     ZN IS IN THE RIGHT HALF PLANE AFTER ROTATION BY CI OR -CI
C-----------------------------------------------------------------------
      ZNR = ZI
      ZNI = -ZR
      ZBR = ZR
      ZBI = ZI
      CIDI = -CONER
      INU = INT(SNGL(FNU))
      ANG = HPI*(FNU-DBLE(FLOAT(INU)))
      C2R = DCOS(ANG)
      C2I = DSIN(ANG)
      CAR = C2R
      SAR = C2I
      IN = INU + N - 1
      IN = MOD(IN,4) + 1
      STR = C2R*CIPR(IN) - C2I*CIPI(IN)
      C2I = C2R*CIPI(IN) + C2I*CIPR(IN)
      C2R = STR
      IF (ZI.GT.0.0D0) GO TO 10
      ZNR = -ZNR
      ZBI = -ZBI
      CIDI = -CIDI
      C2I = -C2I
   10 CONTINUE
C-----------------------------------------------------------------------
C     CHECK FOR UNDERFLOW AND OVERFLOW ON FIRST MEMBER
C-----------------------------------------------------------------------
      FN = DMAX1(FNU,1.0D0)
      CALL ZUNHJ(ZNR, ZNI, FN, 1, TOL, PHIR, PHII, ARGR, ARGI, ZETA1R,
     * ZETA1I, ZETA2R, ZETA2I, ASUMR, ASUMI, BSUMR, BSUMI)
      IF (KODE.EQ.1) GO TO 20
      STR = ZBR + ZETA2R
      STI = ZBI + ZETA2I
      RAST = FN/ZABS(STR,STI)
      STR = STR*RAST*RAST
      STI = -STI*RAST*RAST
      S1R = -ZETA1R + STR
      S1I = -ZETA1I + STI
      GO TO 30
   20 CONTINUE
      S1R = -ZETA1R + ZETA2R
      S1I = -ZETA1I + ZETA2I
   30 CONTINUE
      RS1 = S1R
      IF (DABS(RS1).GT.ELIM) GO TO 150
   40 CONTINUE
      NN = MIN0(2,ND)
      DO 90 I=1,NN
        FN = FNU + DBLE(FLOAT(ND-I))
        CALL ZUNHJ(ZNR, ZNI, FN, 0, TOL, PHIR, PHII, ARGR, ARGI,
     *   ZETA1R, ZETA1I, ZETA2R, ZETA2I, ASUMR, ASUMI, BSUMR, BSUMI)
        IF (KODE.EQ.1) GO TO 50
        STR = ZBR + ZETA2R
        STI = ZBI + ZETA2I
        RAST = FN/ZABS(STR,STI)
        STR = STR*RAST*RAST
        STI = -STI*RAST*RAST
        S1R = -ZETA1R + STR
        S1I = -ZETA1I + STI + DABS(ZI)
        GO TO 60
   50   CONTINUE
        S1R = -ZETA1R + ZETA2R
        S1I = -ZETA1I + ZETA2I
   60   CONTINUE
C-----------------------------------------------------------------------
C     TEST FOR UNDERFLOW AND OVERFLOW
C-----------------------------------------------------------------------
        RS1 = S1R
        IF (DABS(RS1).GT.ELIM) GO TO 120
        IF (I.EQ.1) IFLAG = 2
        IF (DABS(RS1).LT.ALIM) GO TO 70
C-----------------------------------------------------------------------
C     REFINE  TEST AND SCALE
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
        APHI = ZABS(PHIR,PHII)
        AARG = ZABS(ARGR,ARGI)
        RS1 = RS1 + DLOG(APHI) - 0.25D0*DLOG(AARG) - AIC
        IF (DABS(RS1).GT.ELIM) GO TO 120
        IF (I.EQ.1) IFLAG = 1
        IF (RS1.LT.0.0D0) GO TO 70
        IF (I.EQ.1) IFLAG = 3
   70   CONTINUE
C-----------------------------------------------------------------------
C     SCALE S1 TO KEEP INTERMEDIATE ARITHMETIC ON SCALE NEAR
C     EXPONENT EXTREMES
C-----------------------------------------------------------------------
        CALL ZAIRY(ARGR, ARGI, 0, 2, AIR, AII, NAI, IDUM)
        CALL ZAIRY(ARGR, ARGI, 1, 2, DAIR, DAII, NDAI, IDUM)
        STR = DAIR*BSUMR - DAII*BSUMI
        STI = DAIR*BSUMI + DAII*BSUMR
        STR = STR + (AIR*ASUMR-AII*ASUMI)
        STI = STI + (AIR*ASUMI+AII*ASUMR)
        S2R = PHIR*STR - PHII*STI
        S2I = PHIR*STI + PHII*STR
        STR = DEXP(S1R)*CSSR(IFLAG)
        S1R = STR*DCOS(S1I)
        S1I = STR*DSIN(S1I)
        STR = S2R*S1R - S2I*S1I
        S2I = S2R*S1I + S2I*S1R
        S2R = STR
        IF (IFLAG.NE.1) GO TO 80
        CALL ZUCHK(S2R, S2I, NW, BRY(1), TOL)
        IF (NW.NE.0) GO TO 120
   80   CONTINUE
        IF (ZI.LE.0.0D0) S2I = -S2I
        STR = S2R*C2R - S2I*C2I
        S2I = S2R*C2I + S2I*C2R
        S2R = STR
        CYR(I) = S2R
        CYI(I) = S2I
        J = ND - I + 1
        YR(J) = S2R*CSRR(IFLAG)
        YI(J) = S2I*CSRR(IFLAG)
        STR = -C2I*CIDI
        C2I = C2R*CIDI
        C2R = STR
   90 CONTINUE
      IF (ND.LE.2) GO TO 110
      RAZ = 1.0D0/ZABS(ZR,ZI)
      STR = ZR*RAZ
      STI = -ZI*RAZ
      RZR = (STR+STR)*RAZ
      RZI = (STI+STI)*RAZ
      BRY(2) = 1.0D0/BRY(1)
      BRY(3) = D1MACH(2)
      S1R = CYR(1)
      S1I = CYI(1)
      S2R = CYR(2)
      S2I = CYI(2)
      C1R = CSRR(IFLAG)
      ASCLE = BRY(IFLAG)
      K = ND - 2
      FN = DBLE(FLOAT(K))
      DO 100 I=3,ND
        C2R = S2R
        C2I = S2I
        S2R = S1R + (FNU+FN)*(RZR*C2R-RZI*C2I)
        S2I = S1I + (FNU+FN)*(RZR*C2I+RZI*C2R)
        S1R = C2R
        S1I = C2I
        C2R = S2R*C1R
        C2I = S2I*C1R
        YR(K) = C2R
        YI(K) = C2I
        K = K - 1
        FN = FN - 1.0D0
        IF (IFLAG.GE.3) GO TO 100
        STR = DABS(C2R)
        STI = DABS(C2I)
        C2M = DMAX1(STR,STI)
        IF (C2M.LE.ASCLE) GO TO 100
        IFLAG = IFLAG + 1
        ASCLE = BRY(IFLAG)
        S1R = S1R*C1R
        S1I = S1I*C1R
        S2R = C2R
        S2I = C2I
        S1R = S1R*CSSR(IFLAG)
        S1I = S1I*CSSR(IFLAG)
        S2R = S2R*CSSR(IFLAG)
        S2I = S2I*CSSR(IFLAG)
        C1R = CSRR(IFLAG)
  100 CONTINUE
  110 CONTINUE
      RETURN
  120 CONTINUE
      IF (RS1.GT.0.0D0) GO TO 140
C-----------------------------------------------------------------------
C     SET UNDERFLOW AND UPDATE PARAMETERS
C-----------------------------------------------------------------------
      YR(ND) = ZEROR
      YI(ND) = ZEROI
      NZ = NZ + 1
      ND = ND - 1
      IF (ND.EQ.0) GO TO 110
      CALL ZUOIK(ZR, ZI, FNU, KODE, 1, ND, YR, YI, NUF, TOL, ELIM, ALIM)
      IF (NUF.LT.0) GO TO 140
      ND = ND - NUF
      NZ = NZ + NUF
      IF (ND.EQ.0) GO TO 110
      FN = FNU + DBLE(FLOAT(ND-1))
      IF (FN.LT.FNUL) GO TO 130
C      FN = CIDI
C      J = NUF + 1
C      K = MOD(J,4) + 1
C      S1R = CIPR(K)
C      S1I = CIPI(K)
C      IF (FN.LT.0.0D0) S1I = -S1I
C      STR = C2R*S1R - C2I*S1I
C      C2I = C2R*S1I + C2I*S1R
C      C2R = STR
      IN = INU + ND - 1
      IN = MOD(IN,4) + 1
      C2R = CAR*CIPR(IN) - SAR*CIPI(IN)
      C2I = CAR*CIPI(IN) + SAR*CIPR(IN)
      IF (ZI.LE.0.0D0) C2I = -C2I
      GO TO 40
  130 CONTINUE
      NLAST = ND
      RETURN
  140 CONTINUE
      NZ = -1
      RETURN
  150 CONTINUE
      IF (RS1.GT.0.0D0) GO TO 140
      NZ = N
      DO 160 I=1,N
        YR(I) = ZEROR
        YI(I) = ZEROI
  160 CONTINUE
      RETURN
      END
      SUBROUTINE XERROR(MESS,NMESS,L1,L2)
C
C     THIS IS A DUMMY XERROR ROUTINE TO PRINT ERROR MESSAGES WITH NMESS
C     CHARACTERS. L1 AND L2 ARE DUMMY PARAMETERS TO MAKE THIS CALL
C     COMPATIBLE WITH THE SLATEC XERROR ROUTINE. THIS IS A FORTRAN 77
C     ROUTINE.
C
      INTEGER NMESS, L1, L2, NN, NR, K, I, KMIN
      CHARACTER*(*) MESS
      NN=NMESS/70
      NR=NMESS-70*NN
      IF(NR.NE.0) NN=NN+1
      K=1
      PRINT 900
  900 FORMAT(/)
      DO 10 I=1,NN
        KMIN=MIN0(K+69,NMESS)
        PRINT *, MESS(K:KMIN)
        K=K+70
   10 CONTINUE
      PRINT 900
      RETURN
      END
C*** zqcbes.f
-----------------------------------------------------------------
C>>>  ZQCBES.FOR:  Double precision quick check programs
-----------------------------------------------------------------

      PROGRAM ZQCBH
C
C***DATE WRITTEN   830501   (YYMMDD)
C***REVISION DATE  890801, 930101   (YYMMDD)
C
C                *** A DOUBLE PRECISION ROUTINE ***
C
C     ZQCBH IS A QUICK CHECK ROUTINE FOR THE COMPLEX H BESSEL FUNCTIONS
C     GENERATED BY SUBROUTINE ZBESH.
C
C     ZQCBH GENERATES SEQUENCES OF H BESSEL FUNCTIONS FOR KIND 2 FROM
C     ZBESH AND CHECKS THEM AGAINST ANALYTIC CONTINUATION FORMULAS
C     IN THE (Z,FNU) SPACE:
C
C     KODE = 1 TESTS (ANALYTIC CONTINUATION FORMULAE, I**2 = -1):
C
C     H(FNU,2,Z)=-EXP(I*PI*FNU)*H(FNU,1,-Z),       -PI.LT.ARG(Z).LE.0
C
C               = 2*COS(PI*FNU)*H(FNU,2,-Z) + EXP(I*PI*FNU)*H(FNU,1,-Z),
C
C                                                   0.LT.ARG(Z).LE.PI
C
C     KODE = 2 TESTS FOR KINDS 1 AND 2:
C
C            EXP(-I*Z)*H(FNU,1,Z) = [EXP(-I*Z)*H(FNU,1,Z)]
C
C            EXP( I*Z)*H(FNU,2,Z) = [EXP( I*Z)*H(FNU,2,Z)]
C
C     WHERE THE LEFT SIDE IS COMPUTED WITH KODE = 1 AND THE RIGHT SIDE
C     WITH KODE = 2.
C
C     THE PARAMETER MQC CAN HAVE VALUES 1 (THE DEFAULT) FOR A FASTER,
C     LESS DEFINITIVE TEST OR 2 FOR A SLOWER, MORE DEFINITIVE TEST.
C
C     MACHINE CONSTANTS ARE DEFINED IN FUNCTIONS I1MACH, R1MACH, AND
C     D1MACH. THESE MUST BE SELECTED BY THE USER OR SET ACCORDING TO
C     PROLOGUE INSTRUCTIONS.
C
C     COMPLEX CW, CI, U, V, W, Y, Z, ZN, CSGN
      EXTERNAL ZABS
      DOUBLE PRECISION AA, AB, AER, ALIM, ATOL, AV, CT, DIG, ERR, ELIM,
     * EPS, ER, ERTOL, FNU, FNUL, PI, R, RL, RM, D1M4, D1M5, R2, ST,
     * T, TOL, TS, XNU, D1MACH, SLAK, FILM, STR, STI, UR, UI, VR, VI,
     * WR, WI, YR, YI, CWR, CWI, CSGNR, CSGNI, ZR, ZI, ZNR, ZNI, ZABS
      INTEGER I, ICASE, IERR, IHP, IL, IR, IRB, IT, ITL, K, KODE, KK,
     *K1, K2, LFLG, LUN, MFLG, M, N, NU, NZ1, NZ2, NZ3, I1MACH, KEPS,
     *MQC, NL, NUL, KDO
      DIMENSION T(20), AER(20), XNU(20), UR(20), UI(20), VR(20), VI(20),
     * WR(20), WI(20), YR(20), YI(20), KEPS(20), KDO(20)
      DATA LUN /7/
      PARAMETER (MQC=1)
      OPEN(LUN,FILE='ZQCBH.OUT')
C-----------------------------------------------------------------------
C     SET PARAMETERS RELATED TO MACHINE CONSTANTS.
C     TOL IS THE APPROXIMATE UNIT ROUNDOFF LIMITED TO 1.0E-18.
C     ELIM IS THE APPROXIMATE EXPONENTIAL OVER- AND UNDERFLOW LIMIT.
C     EXP(-ELIM).LT.EXP(-ALIM)=EXP(-ELIM)/TOL    AND
C     EXP(ELIM).GT.EXP(ALIM)=EXP(ELIM)*TOL       ARE INTERVALS NEAR
C     UNDERFLOW AND OVERFLOW LIMITS WHERE SCALED ARITHMETIC IS DONE.
C     RL IS THE LOWER BOUNDARY OF THE ASYMPTOTIC EXPANSION FOR LARGE Z.
C     DIG = NUMBER OF BASE 10 DIGITS IN TOL = 10**(-DIG).
C     FNUL IS THE LOWER BOUNDARY OF THE ASYMPTOTIC SERIES FOR LARGE FNU
C-----------------------------------------------------------------------
      D1M4 = D1MACH(4)
      TOL = MAX(D1M4,1.0D-18)
      AA = -DLOG10(D1M4)
      K1 = I1MACH(15)
      K2 = I1MACH(16)
      D1M5 = D1MACH(5)
      K = MIN(IABS(K1),IABS(K2))
      ELIM = 2.303D0*(DBLE(FLOAT(K))*D1M5-3.0D0)
      AB = AA*2.303D0
      ALIM = ELIM + MAX(-AB,-41.45D0)
      DIG = MIN(AA,18.0D0)
      FNUL = 10.0D0 + 6.0D0*(DIG-3.0D0)
      RL = 1.2D0*DIG + 3.0D0
      SLAK = 3.0D0+4.0D0*(-DLOG10(TOL)-7.0D0)/11.0D0
      SLAK = MAX(SLAK,3.0D0)
      ERTOL = TOL*10.0D0**SLAK
      RM = 0.5D0*(ALIM + ELIM)
      RM = MIN(RM,200.0D0)
      RM = MAX(RM,RL+10.0D0)
      R2 = MIN(FNUL,RM)
C-----------------------------------------------------------------------
      WRITE (LUN,99999)
99999 FORMAT (' QUICK CHECK ROUTINE FOR THE H BESSEL FUNCTIONS FROM ZBES
     *H'/)
      WRITE (LUN,99998)
99998 FORMAT (' PARAMETERS TOL,ELIM,ALIM,RL,FNUL,DIG')
      WRITE (LUN,99997) TOL, ELIM, ALIM, RL, FNUL, DIG
99997 FORMAT (6D12.4/)
      ATOL = 100.0D0*TOL
      PI = 4.0D0*DATAN(1.0D0)
      WRITE (LUN,99996) MQC
99996 FORMAT (/' CHECKS IN THE (Z,FNU) SPACE WITH MQC = ',I2/)
C-----------------------------------------------------------------------
C     TEST VALUES OF Z IN -PI.LT.ARG(Z).LE.PI NEAR FORMULA BOUNDARIES
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C     KDO(K), K=1,IL  DETERMINES WHICH OF THE IL ANGLES IN -PI TO PI
C     ARE USE TO COMPUTE VALUES OF Z
C       KDO(K) = 0  MEANS THAT THE INDEX K WILL BE USED FOR ONE OR TWO
C                   VALUES OF Z, DEPENDING ON THE CHOICE OF KEPS(K)
C              = 1  MEANS THAT THE INDEX K AND THE CORRESPONDING ANGLE
C                   WILL BE SKIPPED
C     KEPS(K), K=1,IL DETERMINES WHICH OF THE ANGLES GET INCREMENTED
C     UP AND DOWN TO PUT VALUES OF Z IN REGIONS WHERE DIFFERENT
C     FORMULAE ARE USED.
C       KEPS(K) =0  MEANS THAT THE ANGLE WILL BE USED WITHOUT CHANGE
C               =1  MEANS THAT THE ANGLE WILL BE INCREMENTED UP AND
C                   DOWN BY EPS
C     THE ANGLES TO BE USED ARE STORED IN THE T(I) ARRAY, I=1,ITL
C-----------------------------------------------------------------------
      IF (MQC.NE.2) THEN
        NL=2
        IL=5
        DO 5 I=1,IL
          KEPS(I)=0
          KDO(I)=0
    5   CONTINUE
        NUL=5
        XNU(1) = 0.0D0
        XNU(2) = 1.0D0
        XNU(3) = 2.0D0
        XNU(4) = 0.5D0*FNUL
        XNU(5) = FNUL + 1.1D0
      ELSE
        NL=4
        IL=13
        DO 6 I=1,IL
          KDO(I)=0
          KEPS(I)=0
    6   CONTINUE
        KDO(2)=1
        KDO(6)=1
        KDO(8)=1
        KDO(12)=1
        KEPS(3)=1
        KEPS(4)=1
        KEPS(5)=1
        KEPS(9)=1
        KEPS(10)=1
        KEPS(11)=1
        NUL=6
        XNU(1) = 0.0D0
        XNU(2) = 0.6D0
        XNU(3) = 1.3D0
        XNU(4) = 2.0D0
        XNU(5) = 0.5D0*FNUL
        XNU(6) = FNUL + 1.1D0
      ENDIF
      I = 2
      EPS = 0.01D0
      FILM=DBLE(FLOAT(IL-1))
      T(1) = -PI + EPS
      DO 30 K=2,IL
        IF (KDO(K).EQ.0) THEN
          T(I) = PI*DBLE(FLOAT(-IL+2*K-1))/FILM
          IF (KEPS(K).EQ.0) GO TO 20
          TS=T(I)
          T(I) = TS - EPS
          I = I + 1
          T(I) = TS + EPS
   20     CONTINUE
          I = I + 1
        ENDIF
   30 CONTINUE
      ITL = I - 1
      LFLG = 0
      DO 170 KODE=1,2
        DO 160 N=1,NL
          DO 150 NU=1,NUL
            FNU = XNU(NU)
            DO 140 ICASE=1,3
              IRB = MIN(ICASE,2)
              DO 130 IR=IRB,3
                GO TO (50, 60, 70), ICASE
   50           CONTINUE
                R =(EPS*DBLE(FLOAT(3-IR))+2.0D0*DBLE(FLOAT(IR-1)))/2.0D0
                GO TO 80
   60           CONTINUE
                R = (2.0D0*DBLE(FLOAT(3-IR))+R2*DBLE(FLOAT(IR-1)))/2.0D0
                GO TO 80
   70           CONTINUE
                IF (R2.GE.RM) GO TO 140
                R = (R2*DBLE(FLOAT(3-IR))+RM*DBLE(FLOAT(IR-1)))/2.0D0
   80           CONTINUE
                DO 120 IT=1,ITL
                  CT = COS(T(IT))
                  ST = SIN(T(IT))
                  IF (ABS(CT).LT.ATOL) CT = 0.0D0
                  IF (ABS(ST).LT.ATOL) ST = 0.0D0
                  ZR = R*CT
                  ZI = R*ST
                  IF (KODE.EQ.1) THEN
                    M=2
                    CALL ZBESH(ZR,ZI,FNU,KODE,M,N,YR,YI,NZ1,IERR)
                    IF (IERR.NE.0.OR.NZ1.NE.0) GO TO 120
                    IF (ST.LT.0.0D0 .OR. (ST.EQ.0.0D0.AND.CT.GT.0.0D0))
     *              THEN
                      IHP = 1
                      ZNR = -ZR
                      ZNI = -ZI
                      M=1
                      CALL ZBESH(ZNR,ZNI,FNU,KODE,M,N,WR,WI,NZ2,IERR)
                      IF (IERR.NE.0.OR.NZ2.NE.0) GO TO 120
                    ELSE
                      IHP = 2
                      ZNR = -ZR
                      ZNI = -ZI
                      M=2
                      CALL ZBESH(ZNR,ZNI,FNU,KODE,M,N,WR,WI,NZ3,IERR)
                      IF (IERR.NE.0.OR.NZ3.NE.0) GO TO 120
                      M=1
                      CALL ZBESH(ZNR,ZNI,FNU,KODE,M,N,VR,VI,NZ2,IERR)
                      IF (IERR.NE.0.OR.NZ2.NE.0) GO TO 120
                    ENDIF
                    AB=MOD(FNU,2.0D0)*PI
                    CSGNR = COS(AB)
                    CSGNI = SIN(AB)
                    MFLG = 0
                    DO 100 I=1,N
                      AB = FNU+DBLE(FLOAT(I-1))
                      AA = MAX(0.5D0,AB)
                      IF (IHP.EQ.1) THEN
                        VR(I) = -(CSGNR*WR(I)-CSGNI*WI(I))
                        VI(I) = -(CSGNR*WI(I)+CSGNI*WR(I))
                        CWR = YR(I) - VR(I)
                        CWI = YI(I) - VI(I)
                      ELSE
                        CWR = CSGNR+CSGNR
                        STR =   CWR*WR(I) + CSGNR*VR(I)-CSGNI*VI(I)
                        VI(I) = CWR*WI(I) + CSGNR*VI(I)+CSGNI*VR(I)
                        VR(I) = STR
                        CWR = YR(I) - VR(I)
                        CWI = YI(I) - VI(I)
                      ENDIF
                      AV = ZABS(YR(I),YI(I))
                      ER = ZABS(CWR,CWI)
                      IF(ZI.EQ.0.0D0) THEN
                        IF(ABS(ZR).LT.AA) ER = ER/AV
                      ELSE
                        ER = ER/AV
                      ENDIF
                      AER(I) = ER
                      IF (ER.GT.ERTOL) MFLG = 1
                      CSGNR = -CSGNR
                      CSGNI = -CSGNI
  100               CONTINUE
                  ELSE
                    M=1
                    KK=1
                    CALL ZBESH(ZR,ZI,FNU,KK,M,N,UR,UI,NZ1,IERR)
                    IF (IERR.NE.0.OR.NZ1.NE.0) GO TO 120
                    CALL ZBESH(ZR,ZI,FNU,KODE,M,N,VR,VI,NZ2,IERR)
                    IF (IERR.NE.0.OR.NZ2.NE.0) GO TO 120
                    M=2
                    KK=1
                    CALL ZBESH(ZR,ZI,FNU,KK,M,N,WR,WI,NZ1,IERR)
                    IF (IERR.NE.0.OR.NZ1.NE.0) GO TO 120
                    CALL ZBESH(ZR,ZI,FNU,KODE,M,N,YR,YI,NZ2,IERR)
                    IF (IERR.NE.0.OR.NZ2.NE.0) GO TO 120
                    ZNR = -ZI
                    ZNI =  ZR
                    CALL ZEXP(ZNR,ZNI,ZNR,ZNI)
                    MFLG = 0
                    DO 105 I=1,N
                      AB = FNU+DBLE(FLOAT(I-1))
                      AA = MAX(0.5D0,AB)
                      CALL ZDIV(UR(I),UI(I),ZNR,ZNI,STR,STI)
                      CWR = STR - VR(I)
                      CWI = STI - VI(I)
                      AV = ZABS(VR(I),VI(I))
                      ER = ZABS(CWR,CWI)
                      IF(ZI.EQ.0.0D0) THEN
                        IF(ABS(ZR).LT.AA) ER = ER/AV
                      ELSE
                        ER = ER/AV
                      ENDIF
                      ERR = ER
                      IF (ER.GT.ERTOL) MFLG = 1
                      CWR = ZNR*WR(I) - ZNI*WI(I) - YR(I)
                      CWI = ZNR*WI(I) + ZNI*WR(I) - YI(I)
                      AV = ZABS(YR(I),YI(I))
                      ER = ZABS(CWR,CWI)
                      IF(ZI.EQ.0.0D0) THEN
                        IF(ABS(ZR).LT.AA) ER = ER/AV
                      ELSE
                        ER = ER/AV
                      ENDIF
                      IF (ER.GT.ERTOL) MFLG = 1
                      AER(I) = ER+ERR
  105               CONTINUE
                  ENDIF
                  IF (MFLG.EQ.0) GO TO 120
                  IF (LFLG.EQ.1) GO TO 110
                  WRITE (LUN,99995) ERTOL
99995             FORMAT (/' CASES WHICH VIOLATE THE RELATIVE ERROR TEST
     * WITH ERTOL =', D12.4/)
                  WRITE (LUN,99994)
99994             FORMAT (/' OUTPUT FORMAT'/' KODE,N,IR,IT,ICASE')
                  WRITE (LUN,99993)
99993             FORMAT (' ER(K),K=1,N'/' Z,FNU,V(1),Y(1)')
                  LFLG = 1
  110             CONTINUE
                  WRITE (LUN,99992) KODE, N, IR, IT, ICASE
99992             FORMAT (5I5)
                  WRITE (LUN,99991) (AER(K),K=1,N)
                  WRITE (LUN,99991) ZR,ZI,FNU,VR(1),VI(1),YR(1),YI(1)
99991             FORMAT (7D12.4)
  120           CONTINUE
  130         CONTINUE
  140       CONTINUE
  150     CONTINUE
  160   CONTINUE
  170 CONTINUE
      IF (LFLG.EQ.0) WRITE (LUN,99990)
99990 FORMAT (/' QUICK CHECKS OK'/)
      STOP
      END
      PROGRAM ZQCBI
C
C***DATE WRITTEN   830501   (YYMMDD)
C***REVISION DATE  890801, 930101   (YYMMDD)
C
C                *** A DOUBLE PRECISION ROUTINE ***
C
C     ZQCBI IS A QUICK CHECK ROUTINE FOR THE COMPLEX I BESSEL FUNCTION
C     GENERATED BY SUBROUTINE ZBESI.
C
C     ZQCBK GENERATES SEQUENCES OF I AND K BESSEL FUNCTIONS FROM
C     ZBESI AND CBESK AND CHECKS THE WRONSKIAN EVALUATION
C
C           I(FNU,Z)*K(FNU+1,Z) + I(FNU+1,Z)*K(FNU,Z) = 1/Z
C
C     IN THE RIGHT HALF PLANE AND A MODIFIED FORM
C
C          I(FNU+1,Z)*K(FNU,ZR) - I(FNU,Z)*K(FNU+1,ZR) = C/Z
C
C     IN THE LEFT HALF PLANE WHERE ZR=-Z AND C=EXP(I*FNU*SGN) WITH
C     SGN=+1 FOR IM(Z).GE.0 AND SGN=-1 FOR IM(Z).LT.0.
C
C     THE PARAMETER MQC CAN HAVE VALUES 1 (THE DEFAULT) FOR A FASTER,
C     LESS DEFINITIVE TEST OR 2 FOR A SLOWER, MORE DEFINITIVE TEST.
C
C     MACHINE CONSTANTS ARE DEFINED IN FUNCTIONS I1MACH, R1MACH, AND
C     D1MACH. THESE MUST BE SELECTED BY THE USER OR SET ACCORDING TO
C     PROLOGUE INSTRUCTIONS.
C
C     COMPLEX  CONE, CSGN, CC, CV, CW, CY, W, Y, Z, ZR
      EXTERNAL ZABS
      DOUBLE PRECISION AA, AB, AER, ALIM, ARG, ATOL, CCR, CONER, CONEI,
     * CSGNI, CSGNR, CWI, CWR, CYI, CYR, CVR, CVI, DIG, ELIM, EPS, ER,
     * ERTOL, FNU, FNUL, GNU, HPI, PI, R, RL, R1M4, R1M5, R2, RM,
     * STI, STR, T, TOL, WI, WR, YI, YR, ZI, ZR, ZRI, ZRR, D1MACH, ZABS,
     * FILM, SLAK, TS, ST, CT, FFNU, XNU
      INTEGER I, ICASE, IL, IFNU, IPRNT, IR, IT, ITL, K, KK, KODE, K1,
     * K2, LFLG, LUN, MFLG, N, NZ, N1, NL, NU, NUL, MQC, IERR, I1MACH,
     * KEPS, KDO, IRB
      DIMENSION T(20), AER(20), YR(20), YI(20), WR(20), WI(20), XNU(20),
     *  KEPS(20), KDO(20)
      DATA LUN /7/
      PARAMETER (MQC=1)
      OPEN(LUN, FILE='ZQCBI.OUT')
C-----------------------------------------------------------------------
C     SET PARAMETERS RELATED TO MACHINE CONSTANTS.
C     TOL IS THE APPROXIMATE UNIT ROUNDOFF LIMITED TO 1.0E-18.
C     ELIM IS THE APPROXIMATE EXPONENTIAL OVER- AND UNDERFLOW LIMIT.
C     EXP(-ELIM).LT.EXP(-ALIM)=EXP(-ELIM)/TOL    AND
C     EXP(ELIM).GT.EXP(ALIM)=EXP(ELIM)*TOL       ARE INTERVALS NEAR
C     UNDERFLOW AND OVERFLOW LIMITS WHERE SCALED ARITHMETIC IS DONE.
C     RL IS THE LOWER BOUNDARY OF THE ASYMPTOTIC EXPANSION FOR LARGE Z.
C     DIG = NUMBER OF BASE 10 DIGITS IN TOL = 10**(-DIG).
C     FNUL IS THE LOWER BOUNDARY OF THE ASYMPTOTIC SERIES FOR LARGE FNU
C-----------------------------------------------------------------------
      R1M4 = D1MACH(4)
      TOL = MAX(R1M4,1.0D-18)
      AA = -DLOG10(R1M4)
      K1 = I1MACH(15)
      K2 = I1MACH(16)
      R1M5 = D1MACH(5)
      K = MIN(IABS(K1),IABS(K2))
      ELIM = 2.303D0*(DBLE(FLOAT(K))*R1M5-3.0D0)
      AB = AA*2.303D0
      ALIM = ELIM + MAX(-AB,-41.45D0)
      DIG = MIN(AA,18.0D0)
      FNUL = 10.0D0 + 6.0D0*(DIG-3.0D0)
      RL = 1.2D0*DIG + 3.0D0
      SLAK = 3.0D0+4.0D0*(-DLOG10(TOL)-7.0D0)/11.0D0
      SLAK = MAX(SLAK,3.0D0)
      ERTOL = TOL*10.0D0**SLAK
      RM = 0.5D0*(ALIM + ELIM)
      RM = MIN(RM,200.0D0)
      RM = MAX(RM,RL+10.0D0)
      R2 = MIN(RM,FNUL)
C-----------------------------------------------------------------------
      WRITE (LUN,99999)
99999 FORMAT (' QUICK CHECK ROUTINE FOR THE I BESSEL FUNCTION FROM ZBESI
     *'/)
      WRITE (LUN,99998)
99998 FORMAT (' PARAMETERS TOL,ELIM,ALIM,RL,FNUL,DIG')
      WRITE (LUN,99997) TOL, ELIM, ALIM, RL, FNUL, DIG
99997 FORMAT (6D12.4/)
      CONER = 1.0D0
      CONEI = 0.0D0
      ATOL = 100.0D0*TOL
      HPI = 2.0D0*DATAN(1.0D0)
      PI = HPI + HPI
      WRITE (LUN,99996) MQC
99996 FORMAT (/' CHECKS IN THE (Z,FNU) SPACE WITH MQC = ',I2/)
C-----------------------------------------------------------------------
C     TEST VALUES OF Z IN -PI.LT.ARG(Z).LE.PI NEAR FORMULA BOUNDARIES
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C     KDO(K), K=1,IL  DETERMINES WHICH OF THE IL ANGLES IN -PI TO PI
C     ARE USE TO COMPUTE VALUES OF Z
C       KDO(K) = 0  MEANS THAT THE INDEX K WILL BE USED FOR ONE OR TWO
C                   VALUES OF Z, DEPENDING ON THE CHOICE OF KEPS(K)
C              = 1  MEANS THAT THE INDEX K AND THE CORRESPONDING ANGLE
C                   WILL BE SKIPPED
C     KEPS(K), K=1,IL DETERMINES WHICH OF THE ANGLES GET INCREMENTED
C     UP AND DOWN TO PUT VALUES OF Z IN REGIONS WHERE DIFFERENT
C     FORMULAE ARE USED.
C       KEPS(K) =0  MEANS THAT THE ANGLE WILL BE USED WITHOUT CHANGE
C               =1  MEANS THAT THE ANGLE WILL BE INCREMENTED UP AND
C                   DOWN BY EPS
C     THE ANGLES TO BE USED ARE STORED IN THE T(I) ARRAY, I=1,ITL
C-----------------------------------------------------------------------
      IF (MQC.NE.2) THEN
        NL=2
        IL=5
        DO 5 I=1,IL
          KEPS(I)=0
          KDO(I)=0
    5   CONTINUE
        NUL=5
        XNU(1) = 0.0E0
        XNU(2) = 1.0E0
        XNU(3) = 2.0E0
        XNU(4) = 0.5E0*FNUL
        XNU(5) = FNUL + 1.1E0
      ELSE
        NL=4
        IL=13
        DO 6 I=1,IL
          KDO(I)=0
          KEPS(I)=0
    6   CONTINUE
        KDO(2)=1
        KDO(6)=1
        KDO(8)=1
        KDO(12)=1
        KEPS(3)=1
        KEPS(4)=1
        KEPS(5)=1
        KEPS(9)=1
        KEPS(10)=1
        KEPS(11)=1
        NUL=6
        XNU(1) = 0.0E0
        XNU(2) = 0.6E0
        XNU(3) = 1.3E0
        XNU(4) = 2.0E0
        XNU(5) = 0.5E0*FNUL
        XNU(6) = FNUL + 1.1E0
      ENDIF
      I = 2
      EPS = 0.01D0
      FILM=DBLE(FLOAT(IL-1))
      T(1) = -PI + EPS
      DO 30 K=2,IL
        IF (KDO(K).EQ.0) THEN
          T(I) = PI*DBLE(FLOAT(-IL+2*K-1))/FILM
          IF (KEPS(K).EQ.0) GO TO 20
          TS=T(I)
          T(I) = TS - EPS
          I = I + 1
          T(I) = TS + EPS
   20     CONTINUE
          I = I + 1
        ENDIF
   30 CONTINUE
      ITL = I - 1
      LFLG = 0
      DO 200 KODE=1,2
        DO 190 N=1,NL
          N1 = N + 1
          DO 180 NU=1,NUL
            FNU = XNU(NU)
            IFNU = INT(SNGL(FNU))
            FFNU = FNU - DBLE(FLOAT(IFNU))
            ARG = PI*FFNU
            CSGNR = DCOS(ARG)
            CSGNI = DSIN(ARG)
            IF (MOD(IFNU,2).EQ.0) GO TO 50
            CSGNR = -CSGNR
            CSGNI = -CSGNI
   50       CONTINUE
            DO 170 ICASE=1,3
              IRB = MIN(2,ICASE)
              DO 160 IR=IRB,4
                GO TO (60, 70, 80), ICASE
   60           CONTINUE
                R = (0.2D0*DBLE(FLOAT(4-IR))+2.0D0*DBLE(FLOAT(IR-1)))/
     *           3.0D0
                GO TO 90
   70           CONTINUE
                R = (2.0D0*DBLE(FLOAT(4-IR))+R2*DBLE(FLOAT(IR-1)))/3.0D0
                GO TO 90
   80           CONTINUE
                IF (R2.EQ.RM) GO TO 180
                R = (R2*DBLE(FLOAT(4-IR))+RM*DBLE(FLOAT(IR-1)))/3.0D0
   90           CONTINUE
                DO 150 IT=1,ITL
                  CT = COS(T(IT))
                  ST = SIN(T(IT))
                  IF (DABS(CT).LT.ATOL) CT = 0.0D0
                  IF (DABS(ST).LT.ATOL) ST = 0.0D0
                  ZR = R*CT
                  ZI = R*ST
                  IF (CT.GE.0.0E0) THEN
C-----------------------------------------------------------------------
C     WRONSKIAN CHECKS IN THE RIGHT HALF PLANE
C-----------------------------------------------------------------------
                    CALL ZBESI(ZR, ZI, FNU, KODE, N1, WR, WI, NZ, IERR)
                    IF (NZ.NE.0.OR.IERR.NE.0) GO TO 150
                    CALL ZBESK(ZR, ZI, FNU, KODE, N1, YR, YI, NZ, IERR)
                    IF (NZ.NE.0.OR.IERR.NE.0) GO TO 150
C-----------------------------------------------------------------------
C     ADJUSTMENTS TO WRONSKIAN DUE TO SCALING OF I AND K FUNCTIONS
C     ON KODE=2
C-----------------------------------------------------------------------
                    CALL ZDIV(CONER,CONEI,ZR,ZI,CVR,CVI)
                    IF (KODE.EQ.2) THEN
                      STR = COS(ZI)
                      STI = SIN(ZI)
                      AA = STR*CVR - STI*CVI
                      CVI = STR*CVI + STI*CVR
                      CVR = AA
                    ENDIF
                    CCR = 1.0D0
                  ELSE
C-----------------------------------------------------------------------
C     WRONSKIAN CHECKS IN THE LEFT HALF PLANE
C-----------------------------------------------------------------------
                    ZRR = -ZR
                    ZRI = -ZI
                    CALL ZBESI(ZR, ZI, FNU, KODE, N1, YR, YI, NZ, IERR)
                    IF (NZ.NE.0.OR.IERR.NE.0) GO TO 150
                    CALL ZBESK(ZRR,ZRI, FNU, KODE, N1, WR, WI, NZ, IERR)
                    IF (NZ.NE.0.OR.IERR.NE.0) GO TO 150
                    CVR = CSGNR
                    CVI = CSGNI
                    IF (ZI.LT.0.0E0) THEN
                      CVI = -CVI
                    ENDIF
                    CALL ZDIV(CVR,CVI,ZR,ZI,CVR,CVI)
                    IF (KODE.EQ.2) THEN
C-----------------------------------------------------------------------
C     ADJUSTMENTS TO WRONSKIAN DUE TO SCALING OF I AND K FUNCTIONS
C     ON KODE=2. SCALE FACTOR = EXP(-I*YY) FOR RE(Z).LT.0
C-----------------------------------------------------------------------
                      CWR = COS(ZI)
                      CWI = -SIN(ZI)
                      STR = CVR*CWR-CVI*CWI
                      CVI = CVR*CWI+CVI*CWR
                      CVR = STR
                    ENDIF
                    CCR = -1.0D0
                  ENDIF
                  MFLG = 0
                  KK = 0
                  DO 130 I=1,N
                    CWR = WR(I)*YR(I+1)-WI(I)*YI(I+1)
                    CWI = WR(I)*YI(I+1)+WI(I)*YR(I+1)
                    CYR = YR(I)*WR(I+1)-YI(I)*WI(I+1)
                    CYI = YR(I)*WI(I+1)+YI(I)*WR(I+1)
                    CYR = CCR*CYR
                    CYI = CCR*CYI
                    CYR = CYR + CWR - CVR
                    CYI = CYI + CWI - CVI
                    ER = ZABS(CYR,CYI)/ZABS(CVR,CVI)
                    AER(I) = ER
                    IF (ER.GT.ERTOL) THEN
                      IF(KK.EQ.0) THEN
                        MFLG = 1
                        KK=I
                      ENDIF
                    ENDIF
                    IF (CT.LT.0.0E0) THEN
                      CVR = -CVR
                      CVI = -CVI
                    ENDIF
  130             CONTINUE
                  IF (MFLG.EQ.0) GO TO 150
                  IF (LFLG.EQ.1) GO TO 140
                  WRITE (LUN,99995) ERTOL
99995             FORMAT (/' CASES WHICH VIOLATE THE RELATIVE ERROR TEST
     * WITH ERTOL =', E12.4/)
                  WRITE (LUN,99994)
99994             FORMAT (/' OUTPUT FORMAT'/' KODE,N,IR,IT,ICASE,KK')
                  WRITE (LUN,99993)
99993             FORMAT (' ER(K),K=1,N'/' Z,FNU,Y(KK)        KK=INDEX O
     *F FIRST NON-ZERO PAIR'/)
                  LFLG = 1
  140             CONTINUE
                  WRITE (LUN,99992) KODE, N, IR, IT, ICASE, KK
99992             FORMAT (6I5)
                  WRITE (LUN,99991) (AER(K),K=1,N)
                  WRITE (LUN,99991) ZR, ZI, FNU, YR(KK), YI(KK)
99991             FORMAT (6D12.4)
  150           CONTINUE
  160         CONTINUE
  170       CONTINUE
  180     CONTINUE
  190   CONTINUE
  200 CONTINUE
      IF (LFLG.EQ.0) WRITE (LUN,99990)
99990 FORMAT (/' QUICK CHECKS OK'/)
      IF (MQC.EQ.1) STOP
C-----------------------------------------------------------------------
C     CHECKS NEAR UNDERFLOW LIMITS ON SERIES(I=1) AND UNIFORM
C     ASYMPTOTIC EXPANSION(I=2)
C-----------------------------------------------------------------------
      WRITE (LUN,99989)
99989 FORMAT (/' CHECKS NEAR UNDERFLOW AND OVERFLOW LIMITS'/)
      ZR = 1.4D0
      ZI = 1.4D0
      IPRNT = 0
      DO 280 I=1,2
        FNU = 10.2D0
        KODE = 1
        N = 20
  230   CONTINUE
        CALL ZBESI(ZR, ZI, FNU, KODE, N, YR, YI, NZ, IERR)
        IF (NZ.NE.0) GO TO 240
        FNU = FNU + 5.0D0
        GO TO 230
  240   CONTINUE
        IF (NZ.LT.10) GO TO 250
        FNU = FNU - 1.0D0
        GO TO 230
  250   CONTINUE
        CALL ZBESK(ZR, ZI, FNU, KODE, 2, WR, WI, NZ, IERR)
        CALL ZDIV(CONER, CONEI, ZR, ZI, STR, STI)
        CYR = WR(1)*YR(2) - WI(1)*YI(2)
        CYI = WR(1)*YI(2) + WI(1)*YR(2)
        CWR = WR(2)*YR(1) - WI(2)*YI(1)
        CWI = WR(2)*YI(1) + WI(2)*YR(1)
        CWR = CWR + CYR - STR
        CWI = CWI + CYI - STI
        ER = ZABS(CWR,CWI)/ZABS(STR,STI)
        IF (ER.LT.ERTOL) GO TO 270
        IF (IPRNT.EQ.1) GO TO 260
        WRITE (LUN,99988)
99988   FORMAT (/' OUTPUT FORMAT/19H ERROR,Z,FNU,KODE,N'/)
        IPRNT = 1
  260   CONTINUE
        WRITE (LUN,99987) ER, ZR, ZI, FNU, KODE, N
99987   FORMAT (4D12.4, 2I5)
  270   CONTINUE
        ZR = RL +RL
        ZI = 0.0D0
  280 CONTINUE
C-----------------------------------------------------------------------
C     CHECK NEAR OVERFLOW LIMITS
C-----------------------------------------------------------------------
      ZR = ELIM
      ZI = 0.0D0
      FNU = 0.0D0
  290 CONTINUE
      CALL ZBESK(ZR, ZI, FNU, KODE, N, YR, YI, NZ, IERR)
      IF (NZ.LT.10) GO TO 300
      IF (NZ.EQ.N) FNU = FNU + 3.0D0
      FNU = FNU + 2.0D0
      GO TO 290
  300 CONTINUE
      GNU = FNU + DBLE(FLOAT(N-2))
      CALL ZBESI(ZR, ZI, GNU, KODE, 2, WR, WI, NZ, IERR)
      CALL ZDIV(CONER, CONEI, ZR, ZI, STR, STI)
      CYR = YR(N-1)*WR(2) - YI(N-1)*WI(2)
      CYI = YR(N-1)*WI(2) + YI(N-1)*WR(2)
      CWR = YR(N)*WR(1) - YI(N)*WI(1)
      CWI = YR(N)*WI(1) + YI(N)*WR(1)
      CWR = CWR + CYR - STR
      CWI = CWI + CYI - STI
      ER = ZABS(CWR,CWI)/ZABS(STR,STI)
      IF (ER.LT.ERTOL) GO TO 320
      IF (IPRNT.EQ.1) GO TO 310
      WRITE (LUN,99988)
      IPRNT = 1
  310 CONTINUE
      WRITE (LUN,99987) ER, ZR, ZI, FNU, KODE, N
  320 CONTINUE
      IF (IPRNT.EQ.0) WRITE (LUN,99990)
      STOP
      END
      PROGRAM ZQCBJ
C
C***DATE WRITTEN   830501   (YYMMDD)
C***REVISION DATE  890801, 930101   (YYMMDD)
C
C                *** A DOUBLE PRECISION ROUTINE ***
C
C     ZQCBJ IS A QUICK CHECK ROUTINE FOR THE COMPLEX J BESSEL FUNCTION
C     GENERATED BY SUBROUTINE ZBESJ.
C
C     ZQCBJ GENERATES SEQUENCES OF J AND H BESSEL FUNCTIONS FROM ZBESJ
C     AND ZBESH AND CHECKS THE WRONSKIANS
C
C     J(FNU,Z)*H(FNU+1,1,Z)-J(FNU+1,Z)*H(FNU,1,Z)=2/(PI*I*Z)   Y.GE.0
C
C     J(FNU,Z)*H(FNU+1,2,Z)-J(FNU+1,Z)*H(FNU,2,Z)=-2/(PI*I*Z)  Y.LT.0
C
C     IN THEIR RESPECTIVE HALF PLANES.
C
C     THE PARAMETER MQC CAN HAVE VALUES 1 (THE DEFAULT) FOR A FASTER,
C     LESS DEFINITIVE TEST OR 2 FOR A SLOWER, MORE DEFINITIVE TEST.
C
C     MACHINE CONSTANTS ARE DEFINED IN FUNCTIONS I1MACH, R1MACH, AND
C     D1MACH. THESE MUST BE SELECTED BY THE USER OR SET ACCORDING TO
C     PROLOGUE INSTRUCTIONS.
C
C     COMPLEX Z, WR, CJ, CH, CON, T1, T2, CER
      EXTERNAL ZABS
      DOUBLE PRECISION AA, AB, AER, ALIM, ATOL, CONR, CONI, WRR, WRI,
     * CERR, CERI, CT, DIG, ELIM, EPS, ER, ERTOL, FNU, FNUL, GNU, HPI,
     * PI, R, RL, RM, R1M4, R1M5, R2, ST, STR, STI, TOL, T1R, T1I, T2R,
     * T2I, T, XNU, ZI, ZR, D1MACH, ZABS, FILM, SLAK, TS, SGN, CJR, CJI,
     * CHR, CHI
      INTEGER I, ICASE, IL, IR, IRB, IT, ITL, K, KK, KODE, K1, K2,
     * LFLG, LUN, M, MFLG, N, NU, NZJ, NZH, I1MACH, KEPS, KDO, NL,
     *NUL, MQC, IERRJ, IERRH
      DIMENSION T(20), AER(25), XNU(20), CJR(20), CJI(20), CHR(20),
     * CHI(20), KEPS(20), KDO(20)
      DATA LUN /7/
      PARAMETER (MQC=1)
      OPEN(LUN, FILE='ZQCBJ.OUT')
C-----------------------------------------------------------------------
C     SET PARAMETERS RELATED TO MACHINE CONSTANTS.
C     TOL IS THE APPROXIMATE UNIT ROUNDOFF LIMITED TO 1.0E-18.
C     ELIM IS THE APPROXIMATE EXPONENTIAL OVER- AND UNDERFLOW LIMIT.
C     EXP(-ELIM).LT.EXP(-ALIM)=EXP(-ELIM)/TOL    AND
C     EXP(ELIM).GT.EXP(ALIM)=EXP(ELIM)*TOL       ARE INTERVALS NEAR
C     UNDERFLOW AND OVERFLOW LIMITS WHERE SCALED ARITHMETIC IS DONE.
C     RL IS THE LOWER BOUNDARY OF THE ASYMPTOTIC EXPANSION FOR LARGE Z.
C     DIG = NUMBER OF BASE 10 DIGITS IN TOL = 10**(-DIG).
C     FNUL IS THE LOWER BOUNDARY OF THE ASYMPTOTIC SERIES FOR LARGE FNU
C-----------------------------------------------------------------------
      R1M4 = D1MACH(4)
      TOL = DMAX1(R1M4,1.0D-18)
      AA = -DLOG10(R1M4)
      K1 = I1MACH(15)
      K2 = I1MACH(16)
      R1M5 = D1MACH(5)
      K = MIN0(ABS(K1),ABS(K2))
      ELIM = 2.303D0*(DBLE(FLOAT(K))*R1M5-3.0D0)
      AB = AA*2.303D0
      ALIM = ELIM + DMAX1(-AB,-41.45D0)
      DIG = DMIN1(AA,18.0D0)
      FNUL = 10.0D0 + 6.0D0*(DIG-3.0D0)
      RL = 1.2D0*DIG + 3.0D0
      SLAK = 3.0D0+4.0D0*(-DLOG10(TOL)-7.0D0)/11.0D0
      SLAK = DMAX1(SLAK,3.0D0)
      ERTOL = TOL*10.0D0**SLAK
      RM = 0.5D0*(ALIM + ELIM)
      RM = DMIN1(RM,200.0D0)
      RM = DMAX1(RM,RL+10.0D0)
      R2 = DMIN1(RM,FNUL)
C-----------------------------------------------------------------------
      WRITE (LUN,99999)
99999 FORMAT (' QUICK CHECK ROUTINE FOR THE J BESSEL FUNCTION FROM ZBESJ
     *'/)
      WRITE (LUN,99998)
99998 FORMAT (' PARAMETERS TOL,ELIM,ALIM,RL,FNUL,DIG')
      WRITE (LUN,99997) TOL, ELIM, ALIM, RL, FNUL, DIG
99997 FORMAT (6D12.4/)
      ATOL = 100.0D0*TOL
      HPI = 2.0D0*DATAN(1.0D0)
      PI = HPI + HPI
      CONR = 0.0D0
      CONI = -1.0D0/HPI
      WRITE (LUN,99996) MQC
99996 FORMAT (/' CHECKS IN THE (Z,FNU) SPACE WITH MQC = ',I2/)
C-----------------------------------------------------------------------
C     TEST VALUES OF Z IN -PI.LT.ARG(Z).LE.PI
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C     KDO(K), K=1,IL  DETERMINES WHICH OF THE IL ANGLES IN -PI TO PI
C     ARE USE TO COMPUTE VALUES OF Z
C       KDO(K) = 0  MEANS THAT THE INDEX K WILL BE USED FOR ONE OR TWO
C                   VALUES OF Z, DEPENDING ON THE CHOICE OF KEPS(K)
C              = 1  MEANS THAT THE INDEX K AND THE CORRESPONDING ANGLE
C                   WILL BE SKIPPED
C     KEPS(K), K=1,IL DETERMINES WHICH OF THE ANGLES GET INCREMENTED
C     UP AND DOWN TO PUT VALUES OF Z IN REGIONS WHERE DIFFERENT
C     FORMULAE ARE USED.
C       KEPS(K) =0  MEANS THAT THE ANGLE WILL BE USED WITHOUT CHANGE
C               =1  MEANS THAT THE ANGLE WILL BE INCREMENTED UP AND
C                   DOWN BY EPS
C     THE ANGLES TO BE USED ARE STORED IN THE T(I) ARRAY, I=1,ITL
C-----------------------------------------------------------------------
      IF (MQC.NE.2) THEN
        NL=2
        IL=5
        DO 5 I=1,IL
          KEPS(I)=0
          KDO(I)=0
    5   CONTINUE
        NUL=5
        XNU(1) = 0.0D0
        XNU(2) = 1.0D0
        XNU(3) = 2.0D0
        XNU(4) = 0.5D0*FNUL
        XNU(5) = FNUL + 1.1D0
      ELSE
        NL=4
        IL=13
        DO 6 I=1,IL
          KDO(I)=0
          KEPS(I)=0
    6   CONTINUE
        KDO(2)=1
        KDO(6)=1
        KDO(8)=1
        KDO(12)=1
        KEPS(3)=1
        KEPS(4)=1
        KEPS(5)=1
        KEPS(9)=1
        KEPS(10)=1
        KEPS(11)=1
        NUL=6
        XNU(1) = 0.0D0
        XNU(2) = 0.6D0
        XNU(3) = 1.3D0
        XNU(4) = 2.0D0
        XNU(5) = 0.5D0*FNUL
        XNU(6) = FNUL + 1.1D0
      ENDIF
      I = 2
      EPS = 0.01D0
      FILM=DBLE(FLOAT(IL-1))
      T(1) = -PI + EPS
      DO 30 K=2,IL
        IF (KDO(K).EQ.0) THEN
          T(I) = PI*DBLE(FLOAT(-IL+2*K-1))/FILM
          IF (KEPS(K).EQ.0) GO TO 20
          TS=T(I)
          T(I) = TS - EPS
          I = I + 1
          T(I) = TS + EPS
   20     CONTINUE
          I = I + 1
        ENDIF
   30 CONTINUE
      ITL = I - 1
      LFLG = 0
      DO 260 KODE=1,2
        DO 250 N=1,NL
          NP = N+1
          DO 240 NU=1,NUL
            FNU = XNU(NU)
            GNU = FNU + DBLE(FLOAT(N-1)) + 1.0D0
            GNU = SQRT(GNU)
            GNU = MIN(GNU, 0.5D0*RL)
            DO 230 ICASE=1,3
              IRB = MIN0(2,ICASE)
              DO 220 IR=IRB,4
                GO TO (50, 60, 70), ICASE
   50           CONTINUE
                R = (GNU*DBLE(FLOAT(4-IR))+2.0D0*DBLE(FLOAT(IR-1)))/
     *           3.0D0
                GO TO 80
   60           CONTINUE
                R = (2.0D0*DBLE(FLOAT(4-IR))+R2*DBLE(FLOAT(IR-1)))/3.0D0
                GO TO 80
   70           CONTINUE
                IF (R2.GE.RM) GO TO 230
                R = (R2*DBLE(FLOAT(4-IR))+RM*DBLE(FLOAT(IR-1)))/3.0D0
   80           CONTINUE
                DO 210 IT=1,ITL
                  CT = COS(T(IT))
                  ST = SIN(T(IT))
                  IF (ABS(CT).LT.ATOL) CT = 0.0D0
                  IF (ABS(ST).LT.ATOL) ST = 0.0D0
                  ZR = R*CT
                  ZI = R*ST
                  IF(ZR.EQ.0.0.AND.ZI.EQ.0.0) GO TO 210
                  CALL ZDIV(CONR,CONI,ZR,ZI,WRR,WRI)
                  M=1
                  IF(ZI.LT.0.0) THEN
                    M=2
                    WRR = -WRR
                    WRI = -WRI
                  ENDIF
                  CALL ZBESJ(ZR,ZI,FNU,KODE,NP,CJR,CJI,NZJ,IERRJ)
                  CALL ZBESH(ZR,ZI,FNU,KODE,M,NP,CHR,CHI,NZH,IERRH)
                  IF(NZJ.NE.0.OR.NZH.NE.0) GO TO 210
                  IF(IERRJ.NE.0.OR.IERRH.NE.0) GO TO 210
                  IF(KODE.EQ.2) THEN
                    SGN = 3.0D0-2.0D0*DBLE(FLOAT(M))
                    STR = COS(ZR)
                    STI = -SGN*SIN(ZR)
                    T1R = WRR*STR - WRI*STI
                    WRI = WRR*STI + WRI*STR
                    WRR = T1R
                  ENDIF
                  KK = 0
                  MFLG = 0
                  DO 190 I=1,N
                    STR = CJR(I)*CHR(I+1) - CJI(I)*CHI(I+1)
                    T1I = CJR(I)*CHI(I+1) + CJI(I)*CHR(I+1)
                    T1R = STR
                    STR = CJR(I+1)*CHR(I) - CJI(I+1)*CHI(I)
                    T2I = CJR(I+1)*CHI(I) + CJI(I+1)*CHR(I)
                    T2R = STR
                    CERR = T1R - T2R -WRR
                    CERI = T1I - T2I -WRI
                    ER=ZABS(CERR,CERI)/ZABS(WRR,WRI)
                    IF (ER.GT.ERTOL) THEN
                      IF(MFLG.EQ.0) THEN
                        MFLG = 1
                        KK=I
                      ENDIF
                    ENDIF
                    AER(I)=ER
  190             CONTINUE
                  IF (MFLG.EQ.0) GO TO 210
                  IF (LFLG.EQ.1) GO TO 200
                  WRITE (LUN,99995) ERTOL
99995             FORMAT (/' CASES WHICH VIOLATE THE RELATIVE ERROR TEST
     * WITH ERTOL =', D12.4/)
                  WRITE (LUN,99994)
99994             FORMAT (/' OUTPUT FORMAT'/' KODE,N,IR,IT,NZJ,NZH,ICASE
     *')
                  WRITE (LUN,99993)
99993             FORMAT (' ER(K),K=1,N'/' Z,FNU,CJ(KK),CH(KK), KK=INDEX
     * OF FIRST NON-ZERO Y,W PAIR'/)
                  LFLG = 1
  200             CONTINUE
                  WRITE (LUN,99992) KODE, N, IR, IT, NZJ, NZH, ICASE
99992             FORMAT (8I5)
                  WRITE (LUN,99991) (AER(K),K=1,N)
                  WRITE (LUN,99991) ZR, ZI, FNU, CJR(KK), CJI(KK),
     *             CHR(KK), CHI(KK)
99991             FORMAT (9D12.4)
  210           CONTINUE
  220         CONTINUE
  230       CONTINUE
  240     CONTINUE
  250   CONTINUE
  260 CONTINUE
      IF (LFLG.EQ.0) WRITE (LUN,99990)
99990 FORMAT (/' QUICK CHECKS OK'/)
      STOP
      END
      PROGRAM ZQCBK
C
C***DATE WRITTEN   830501   (YYMMDD)
C***REVISION DATE  890801, 930101   (YYMMDD)
C
C                *** A DOUBLE PRECISION ROUTINE ***
C
C     ZQCBK IS A QUICK CHECK ROUTINE FOR THE COMPLEX K BESSEL FUNCTION
C     GENERATED BY SUBROUTINE ZBESK.
C
C     ZQCBK GENERATES SEQUENCES OF I AND K BESSEL FUNCTIONS FROM
C     ZBESI AND ZBESK AND CHECKS THEM AGAINST THE WRONSKIAN EVALUATION
C
C           I(FNU,Z)*K(FNU+1,Z) + I(FNU+1,Z)*K(FNU,Z) = 1/Z
C
C     IN THE RIGHT HALF PLANE AND THE ANALYTIC CONTINUATION FORMULA
C     FOR H(FNU,2,Z) IN TERMS OF THE K FUNCTION
C
C           K(FNU,Z) = C3*H(FNU,2,ZR) + C4*H(FNU,1,ZR)    IM(Z).GE.0
C
C                    = CONJG(K(FNU,CONJG(Z)))             IM(Z).LT.0
C
C     IN THE LEFT HALF PLANE WHERE C3=C1*CONJG(C2)*C5, C4 = C2*C5
C     C1=2*COS(PI*FNU),   C2=EXP(PI*FNU*I/2),   C5 =-PI*I/2   AND
C     ZR = Z*EXP(-3*PI*I/2) = Z*I
C
C     THE PARAMETER MQC CAN HAVE VALUES 1 (THE DEFAULT) FOR A FASTER,
C     LESS DEFINITIVE TEST OR 2 FOR A SLOWER, MORE DEFINITIVE TEST.
C
C     MACHINE CONSTANTS ARE DEFINED IN FUNCTIONS I1MACH, R1MACH, AND
C     D1MACH. THESE MUST BE SELECTED BY THE USER OR SET ACCORDING TO
C     PROLOGUE INSTRUCTIONS.
C
C     COMPLEX CONE, CSGN, CV, CW, CY, C1, C2, C3, C4, V, W, Y, Z, ZR,
C    * CIP
      EXTERNAL ZABS
      DOUBLE PRECISION AA, AB, AER, ALIM, ARG, ATOL, CONEI, CONER,
     * CSGNI, CSGNR, CVI, CVR, CWI, CWR, CYI, CYR, DIG, ELIM, EPS,
     * ER, ERTOL, FFNU, FNU, FNUL, HPI, PI, R, RL, RM, R1M4, R1M5, R2,
     * STI, STR, T, TOL, WI, WR, XNU, YI, YR, ZI, C1R, C1I, C2R, C2I,
     * ZRI, ZRR, ZR, D1MACH, ZABS, FILM, CT, ST, TS, SLAK, C3R, C3I,
     * C4R, C4I, VR, VI, CIPR, CIPI, COEI, ZZR, ZZI
      INTEGER I, ICASE, IFNU, IL, IR, IRB, IT, ITL, I4, K, KK, KODE,
     * K1, K2, LFLG, LUN, M, MFLG, N, NU, N1, IERR, I1MACH,
     *KEPS, KDO, NL, NUL, NZ, MQC
      DIMENSION T(20), AER(25), XNU(20), VR(20), VI(20), YR(20), YI(20),
     * WR(20), WI(20), KEPS(20), KDO(20), CIPR(4), CIPI(4)
      DATA LUN /7/
      DATA CIPR(1), CIPI(1), CIPR(2), CIPI(2), CIPR(3),CIPI(3),CIPR(4),
     * CIPI(4) / 1.0D0,0.0D0,0.0D0,1.0D0,-1.0D0,0.0D0,0.0D0,-1.0D0/
      PARAMETER (MQC=1)
      OPEN(LUN, FILE='ZQCBK.OUT')
C-----------------------------------------------------------------------
C     SET PARAMETERS RELATED TO MACHINE CONSTANTS.
C     TOL IS THE APPROXIMATE UNIT ROUNDOFF LIMITED TO 1.0E-18.
C     ELIM IS THE APPROXIMATE EXPONENTIAL OVER- AND UNDERFLOW LIMIT.
C     EXP(-ELIM).LT.EXP(-ALIM)=EXP(-ELIM)/TOL    AND
C     EXP(ELIM).GT.EXP(ALIM)=EXP(ELIM)*TOL       ARE INTERVALS NEAR
C     UNDERFLOW AND OVERFLOW LIMITS WHERE SCALED ARITHMETIC IS DONE.
C     RL IS THE LOWER BOUNDARY OF THE ASYMPTOTIC EXPANSION FOR LARGE Z.
C     DIG = NUMBER OF BASE 10 DIGITS IN TOL = 10**(-DIG).
C     FNUL IS THE LOWER BOUNDARY OF THE ASYMPTOTIC SERIES FOR LARGE FNU
C-----------------------------------------------------------------------
      R1M4 = D1MACH(4)
      TOL = MAX(R1M4,1.0D-18)
      AA = -DLOG10(R1M4)
      K1 = I1MACH(15)
      K2 = I1MACH(16)
      R1M5 = D1MACH(5)
      K = MIN(IABS(K1),IABS(K2))
      ELIM = 2.303D0*(DBLE(FLOAT(K))*R1M5-3.0D0)
      AB = AA*2.303D0
      ALIM = ELIM + MAX(-AB,-41.45D0)
      DIG = MIN(AA,18.0D0)
      FNUL = 10.0D0 + 6.0D0*(DIG-3.0D0)
      RL = 1.2D0*DIG + 3.0D0
      SLAK = 3.0D0+4.0D0*(-DLOG10(TOL)-7.0D0)/11.0D0
      SLAK = MAX(SLAK,3.0D0)
      ERTOL = TOL*10.0D0**SLAK
      RM = 0.5D0*(ALIM + ELIM)
      RM = MIN(RM,200.0D0)
      RM = MAX(RM,RL+10.0D0)
      R2 = MIN(RM,FNUL)
C-----------------------------------------------------------------------
      WRITE (LUN,99999)
99999 FORMAT (54H QUICK CHECK ROUTINE FOR THE K BESSEL FUNCTION FROM ZB,
     * 3HESK/)
      WRITE (LUN,99998)
99998 FORMAT (37H PARAMETERS TOL,ELIM,ALIM,RL,FNUL,DIG)
      WRITE (LUN,99997) TOL, ELIM, ALIM, RL, FNUL, DIG
99997 FORMAT (6D12.4/)
      CONER = 1.0D0
      CONEI = 0.0D0
      ATOL = 100.0D0*TOL
      HPI = 2.0D0*DATAN(1.0D0)
      PI = HPI + HPI
      WRITE (LUN,99996) MQC
99996 FORMAT (/' CHECKS IN THE (Z,FNU) SPACE WITH MQC = ',I2/)
C-----------------------------------------------------------------------
C     TEST VALUES OF Z IN -PI.LT.ARG(Z).LE.PI NEAR FORMULA BOUNDARIES
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C     KDO(K), K=1,IL  DETERMINES WHICH OF THE IL ANGLES IN -PI TO PI
C     ARE USE TO COMPUTE VALUES OF Z
C       KDO(K) = 0  MEANS THAT THE INDEX K WILL BE USED FOR ONE OR TWO
C                   VALUES OF Z, DEPENDING ON THE CHOICE OF KEPS(K)
C              = 1  MEANS THAT THE INDEX K AND THE CORRESPONDING ANGLE
C                   WILL BE SKIPPED
C     KEPS(K), K=1,IL DETERMINES WHICH OF THE ANGLES GET INCREMENTED
C     UP AND DOWN TO PUT VALUES OF Z IN REGIONS WHERE DIFFERENT
C     FORMULAE ARE USED.
C       KEPS(K) =0  MEANS THAT THE ANGLE WILL BE USED WITHOUT CHANGE
C               =1  MEANS THAT THE ANGLE WILL BE INCREMENTED UP AND
C                   DOWN BY EPS
C     THE ANGLES TO BE USED ARE STORED IN THE T(I) ARRAY, I=1,ITL
C-----------------------------------------------------------------------
      IF (MQC.NE.2) THEN
        NL=2
        IL=5
        DO 5 I=1,IL
          KEPS(I)=0
          KDO(I)=0
    5   CONTINUE
        NUL=5
        XNU(1) = 0.0D0
        XNU(2) = 1.0D0
        XNU(3) = 2.0D0
        XNU(4) = 0.5D0*FNUL
        XNU(5) = FNUL + 1.1D0
      ELSE
        NL=4
        IL=13
        DO 6 I=1,IL
          KDO(I)=0
          KEPS(I)=0
    6   CONTINUE
        KDO(2)=1
        KDO(6)=1
        KDO(8)=1
        KDO(12)=1
        KEPS(3)=1
        KEPS(4)=1
        KEPS(5)=1
        KEPS(9)=1
        KEPS(10)=1
        KEPS(11)=1
        NUL=6
        XNU(1) = 0.0D0
        XNU(2) = 0.6D0
        XNU(3) = 1.3D0
        XNU(4) = 2.0D0
        XNU(5) = 0.5D0*FNUL
        XNU(6) = FNUL + 1.1D0
      ENDIF
      I = 2
      EPS = 0.01D0
      FILM=DBLE(FLOAT(IL-1))
      T(1) = -PI + EPS
      DO 30 K=2,IL
        IF (KDO(K).EQ.0) THEN
          T(I) = PI*DBLE(FLOAT(-IL+2*K-1))/FILM
          IF (KEPS(K).EQ.0) GO TO 20
          TS=T(I)
          T(I) = TS - EPS
          I = I + 1
          T(I) = TS + EPS
   20     CONTINUE
          I = I + 1
        ENDIF
   30 CONTINUE
      ITL = I - 1
      LFLG = 0
      DO 200 KODE=1,2
        DO 190 N=1,NL
          N1 = N + 1
          DO 180 NU=1,NUL
            FNU = XNU(NU)
            IFNU = INT(SNGL(FNU))
            FFNU = FNU - DBLE(FLOAT(IFNU))
            ARG = HPI*FFNU
            CSGNR = DCOS(ARG)
            CSGNI = DSIN(ARG)
            I4 = MOD(IFNU,4)+1
            STR = CSGNR*CIPR(I4) - CSGNI*CIPI(I4)
            CSGNI = CSGNR*CIPI(I4) + CSGNI*CIPR(I4)
            CSGNR = STR
   50       CONTINUE
            DO 170 ICASE=1,3
              IRB = MIN(2,ICASE)
              DO 160 IR=IRB,4
                GO TO (60, 70, 80), ICASE
   60           CONTINUE
                R = (0.2D0*DBLE(FLOAT(4-IR))+2.0D0*DBLE(FLOAT(IR-1)))/
     *           3.0D0
                GO TO 90
   70           CONTINUE
                R = (2.0D0*DBLE(FLOAT(4-IR))+R2*DBLE(FLOAT(IR-1)))/3.0D0
                GO TO 90
   80           CONTINUE
                IF (R2.EQ.RM) GO TO 180
                R = (R2*DBLE(FLOAT(4-IR))+RM*DBLE(FLOAT(IR-1)))/3.0D0
   90           CONTINUE
                DO 150 IT=1,ITL
                  CT = COS(T(IT))
                  ST = SIN(T(IT))
                  IF (DABS(CT).LT.ATOL) CT = 0.0D0
                  IF (DABS(ST).LT.ATOL) ST = 0.0D0
                  ZR = R*CT
                  ZI = R*ST
                  IF (ZR.GE.0.0E0) THEN
C-----------------------------------------------------------------------
C     WRONSKIAN CHECKS IN THE RIGHT HALF PLANE
C-----------------------------------------------------------------------
                    CALL ZBESI(ZR, ZI, FNU, KODE, N1, WR, WI, NZ, IERR)
                    IF (NZ.NE.0.OR.IERR.NE.0) GO TO 150
                    CALL ZBESK(ZR, ZI, FNU, KODE, N1, YR, YI, NZ, IERR)
                    IF (NZ.NE.0.OR.IERR.NE.0) GO TO 150
C-----------------------------------------------------------------------
C     ADJUSTMENTS TO WRONSKIAN DUE TO SCALING OF I AND K FUNCTIONS
C     ON KODE=2
C-----------------------------------------------------------------------
                    CALL ZDIV(CONER,CONEI,ZR,ZI,CVR,CVI)
                    IF (KODE.EQ.2) THEN
                      STR = COS(ZI)
                      STI = SIN(ZI)
                      AA = STR*CVR - STI*CVI
                      CVI = STR*CVI + STI*CVR
                      CVR = AA
                    ENDIF
                    MFLG = 0
                    KK = 0
                    DO 130 I=1,N
                      CWR = WR(I)*YR(I+1)-WI(I)*YI(I+1)
                      CWI = WR(I)*YI(I+1)+WI(I)*YR(I+1)
                      CYR = YR(I)*WR(I+1)-YI(I)*WI(I+1)
                      CYI = YR(I)*WI(I+1)+YI(I)*WR(I+1)
                      CYR = CYR + CWR - CVR
                      CYI = CYI + CWI - CVI
                      ER = ZABS(CYR,CYI)/ZABS(CVR,CVI)
                      AER(I) = ER
                      IF (ER.GT.ERTOL) THEN
                        IF(KK.EQ.0) THEN
                          MFLG = 1
                          KK=I
                        ENDIF
                      ENDIF
  130               CONTINUE
                  ELSE
C-----------------------------------------------------------------------
C     ANALYTIC CONTINUATION FORMULA CHECKS FOR LEFT HALF PLANE IN TERMS
C     OF H(FNU,1,Z) AND H(FNU,2,Z)
C-----------------------------------------------------------------------
                    ZZR = ZR
                    ZZI = ZI
                    IF (ZI.LT.0.0E0) THEN
                      ZZI = -ZZI
                    ENDIF
                    ZRR = -ZZI
                    ZRI = ZZR
                    M=1
                    CALL ZBESH(ZRR,ZRI,FNU,KODE,M,N,WR,WI,NZ,IERR)
                    IF (IERR.NE.0) GO TO 150
                    M=2
                    CALL ZBESH(ZRR,ZRI,FNU,KODE,M,N,VR,VI,NZ,IERR)
                    IF (IERR.NE.0) GO TO 150
                    CALL ZBESK(ZR, ZI, FNU, KODE, N, YR, YI, NZ, IERR)
                    IF (NZ.NE.0.OR.IERR.NE.0) GO TO 150
                    COEI = -HPI
                    MFLG = 0
                    KK = 0
                    AA = 2.0D0*COS(PI*FFNU)
                    IF(MOD(IFNU,2).NE.0) AA = -AA
                    C1R = AA
                    C1I = 0.0D0
                    C2R = CSGNR
                    C2I = CSGNI
                    DO 135 I=1,N
                      C3R = C1R
                      C3I = C1I
                      C4R = C2R
                      C4I = C2I
                      IF (KODE.EQ.2) THEN
C-----------------------------------------------------------------------
C     ADJUSTMENTS TO COEFICIENTS DUE TO SCALING OF H(FNU,1,Z) AND
C     H(FNU,2,Z) FUNCTIONS ON KODE = 2.
C-----------------------------------------------------------------------
                        AB=ZABS(VR(I),VI(I))
                        AA=LOG(AB)+ZR+ZR
                        IF (AA.GT.ELIM) GO TO 150
                        IF (AA.LT.-ELIM) THEN
                          C3R = 0.0D0
                          C3I = 0.0D0
                        ELSE
                          STR = ZZR+ZZR
                          STI = ZZI+ZZI
                          CALL ZEXP(STR,STI,CVR,CVI)
                          STR = C3R*CVR - C3I*CVI
                          C3I = C3R*CVI + C3I*CVR
                          C3R = STR
                        ENDIF
                      ENDIF
C                     CY = (C3*CONJG(C2)*V(I)+C4*W(I))*COE
                      STR = C3R*C2R + C3I*C2I
                      STI = -C3R*C2I + C3I*C2R
                      CYR = STR*VR(I) - STI*VI(I)
                      CYI = STR*VI(I) + STI*VR(I)
                      CYR = CYR + (C4R*WR(I) - C4I*WI(I))
                      CYI = CYI + (C4R*WI(I) + C4I*WR(I))
                      STR = -CYI*COEI
                      CYI =  CYR*COEI
                      CYR = STR
                      IF (ZI.LT.0.0D0) THEN
                        CYI = -CYI
                      ENDIF
                      STR = CYR - YR(I)
                      STI = CYI - YI(I)
                      ER = ZABS(STR,STI)/ZABS(YR(I),YI(I))
                      AER(I) = ER
                      IF (ER.GT.ERTOL) THEN
                        IF(KK.EQ.0) THEN
                          MFLG = 1
                          KK=I
                        ENDIF
                      ENDIF
                      STR = -C2I
                      C2I = C2R
                      C2R = STR
                      C1R = -C1R
                      C1I = -C1I
  135               CONTINUE
                  ENDIF
                  IF (MFLG.EQ.0) GO TO 150
                  IF (LFLG.EQ.1) GO TO 140
                  WRITE (LUN,99995) ERTOL
99995             FORMAT (/' CASES WHICH VIOLATE THE RELATIVE ERROR TEST
     * WITH ERTOL =', E12.4/)
                  WRITE (LUN,99994)
99994             FORMAT (/' OUTPUT FORMAT'/' KODE,N,IR,IT,ICASE,KK')
                  WRITE (LUN,99993)
99993             FORMAT (' ER(K),K=1,N'/' Z,FNU,Y(KK)        KK=INDEX O
     *F FIRST NON-ZERO PAIR'/)
                  LFLG = 1
  140             CONTINUE
                  WRITE (LUN,99992) KODE, N, IR, IT, ICASE, KK
99992             FORMAT (6I5)
                  WRITE (LUN,99991) (AER(K),K=1,N)
                  WRITE (LUN,99991) ZR, ZI, FNU, YR(KK), YI(KK)
99991             FORMAT (6D12.4)
  150           CONTINUE
  160         CONTINUE
  170       CONTINUE
  180     CONTINUE
  190   CONTINUE
  200 CONTINUE
      IF (LFLG.EQ.0) WRITE (LUN,99990)
99990 FORMAT (/' QUICK CHECKS OK'/)
      STOP
      END
      PROGRAM ZQCBY
C
C***DATE WRITTEN   830501   (YYMMDD)
C***REVISION DATE  890801, 930101   (YYMMDD)
C
C                *** A DOUBLE PRECISION ROUTINE ***
C
C     ZQCBY IS A QUICK CHECK ROUTINE FOR THE COMPLEX Y BESSEL FUNCTION
C     GENERATED BY SUBROUTINE ZBESY.
C
C     ZQCBY GENERATES SEQUENCES OF Y BESSEL FUNCTIONS FROM ZBESY AND
C     ZBESYH AND COMPARES THEM FOR A VARIETY OF VALUES IN THE (Z,FNU)
C     SPACE. ZBESYH IS AN OLD VERSION OF ZBESY WHICH COMPUTES THE Y
C     FUNCTION FROM THE H FUNCTIONS OF KINDS 1 AND 2.
C
C     THE PARAMETER MQC CAN HAVE VALUES 1 (THE DEFAULT) FOR A FASTER,
C     LESS DEFINITIVE TEST OR 2 FOR A SLOWER, MORE DEFINITIVE TEST.
C
C     MACHINE CONSTANTS ARE DEFINED IN FUNCTIONS I1MACH, R1MACH, AND
C     D1MACH. THESE MUST BE SELECTED BY THE USER OR SET ACCORDING TO
C     PROLOGUE INSTRUCTIONS.
C
C     COMPLEX CW,CWRK,V,W,Z
      EXTERNAL ZABS
      DOUBLE PRECISION AA, AB, AER, ALIM, ATOL, AV, CWI, CWR, CWRKI,
     * CWRKR, DIG, ELIM, EPS, ER, ERTOL, FFNU, FNU, FNUL, HPI, PI, R,
     * RL, RM, R1M4, R1M5, R2, T, TOL, VI, VR, WI, WR, XNU, ZI, ZR,
     * D1MACH, ZABS, ST, CT, TS, SLAK, FILM
      INTEGER I, ICASE, IFNU, IL, IR, IRB, IT, ITL, K, KK,
     * KODE, K1, K2, LFLG, LUN, MFLG, N, NU, NZ1, NZ2, IERR, I1MACH,
     * MQC, KDO, KEPS, NL, NUL
      DIMENSION T(20), AER(20), XNU(20), WR(20), WI(20), VR(20), VI(20),
     * CWRKR(20), CWRKI(20), KDO(20), KEPS(20)
      DATA LUN /7/
      PARAMETER (MQC=1)
      OPEN(LUN, FILE='ZQCBY.OUT')
C-----------------------------------------------------------------------
C     SET PARAMETERS RELATED TO MACHINE CONSTANTS.
C     TOL IS THE APPROXIMATE UNIT ROUNDOFF LIMITED TO 1.0E-18.
C     ELIM IS THE APPROXIMATE EXPONENTIAL OVER- AND UNDERFLOW LIMIT.
C     EXP(-ELIM).LT.EXP(-ALIM)=EXP(-ELIM)/TOL    AND
C     EXP(ELIM).GT.EXP(ALIM)=EXP(ELIM)*TOL       ARE INTERVALS NEAR
C     UNDERFLOW AND OVERFLOW LIMITS WHERE SCALED ARITHMETIC IS DONE.
C     RL IS THE LOWER BOUNDARY OF THE ASYMPTOTIC EXPANSION FOR LARGE Z.
C     DIG = NUMBER OF BASE 10 DIGITS IN TOL = 10**(-DIG).
C     FNUL IS THE LOWER BOUNDARY OF THE ASYMPTOTIC SERIES FOR LARGE FNU
C-----------------------------------------------------------------------
      R1M4 = D1MACH(4)
      TOL = MAX(R1M4,1.0D-18)
      AA = -DLOG10(R1M4)
      K1 = I1MACH(15)
      K2 = I1MACH(16)
      R1M5 = D1MACH(5)
      K = MIN(IABS(K1),IABS(K2))
      ELIM = 2.303D0*(DBLE(FLOAT(K))*R1M5-3.0D0)
      AB = AA*2.303D0
      ALIM = ELIM + MAX(-AB,-41.45D0)
      DIG = MIN(AA,18.0D0)
      FNUL = 10.0D0 + 6.0D0*(DIG-3.0D0)
      RL = 1.2D0*DIG + 3.0D0
      SLAK = 3.0D0+4.0D0*(-DLOG10(TOL)-7.0D0)/11.0D0
      SLAK = MAX(SLAK,3.0D0)
      ERTOL = TOL*10.0D0**SLAK
      RM = 0.5D0*(ALIM + ELIM)
      RM = MIN(RM,200.0D0)
      RM = MAX(RM,RL+10.0D0)
      R2 = MIN(RM,FNUL)
C-----------------------------------------------------------------------
      WRITE (LUN,99999)
99999 FORMAT (' QUICK CHECK ROUTINE FOR THE Y BESSEL FUNCTION FROM ZBESY
     *'/)
      WRITE (LUN,99998)
99998 FORMAT (' PARAMETERS TOL,ELIM,ALIM,RL,FNUL,DIG')
      WRITE (LUN,99997) TOL, ELIM, ALIM, RL, FNUL, DIG
99997 FORMAT (6D12.4/)
      ATOL = 100.0D0*TOL
      HPI = 2.0D0*DATAN(1.0D0)
      PI = HPI + HPI
      WRITE (LUN,99996) MQC
99996 FORMAT (/' CHECKS IN THE (Z,FNU) SPACE WITH MQC = ',I2/)
C-----------------------------------------------------------------------
C     TEST VALUES OF Z IN -PI/2.LT.ARG(Z).LE.PI NEAR BOUNDARIES
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C     KDO(K), K=1,IL  DETERMINES WHICH OF THE IL ANGLES IN -PI TO PI
C     ARE USE TO COMPUTE VALUES OF Z
C       KDO(K) = 0  MEANS THAT THE INDEX K WILL BE USED FOR ONE OR TWO
C                   VALUES OF Z, DEPENDING ON THE CHOICE OF KEPS(K)
C              = 1  MEANS THAT THE INDEX K AND THE CORRESPONDING ANGLE
C                   WILL BE SKIPPED
C     KEPS(K), K=1,IL DETERMINES WHICH OF THE ANGLES GET INCREMENTED
C     UP AND DOWN TO PUT VALUES OF Z IN REGIONS WHERE DIFFERENT
C     FORMULAE ARE USED.
C       KEPS(K) =0  MEANS THAT THE ANGLE WILL BE USED WITHOUT CHANGE
C               =1  MEANS THAT THE ANGLE WILL BE INCREMENTED UP AND
C                   DOWN BY EPS
C     THE ANGLES TO BE USED ARE STORED IN THE T(I) ARRAY, I=1,ITL
C-----------------------------------------------------------------------
      IF (MQC.NE.2) THEN
        NL=2
        IL=5
        DO 5 I=1,IL
          KEPS(I)=0
          KDO(I)=0
    5   CONTINUE
        KDO(5)=1
        NUL=5
        XNU(1) = 0.0D0
        XNU(2) = 1.0D0
        XNU(3) = 2.0D0
        XNU(4) = 0.5D0*FNUL
        XNU(5) = FNUL + 1.2D0
      ELSE
        NL=4
        IL=13
        DO 6 I=1,IL
          KDO(I)=0
          KEPS(I)=0
    6   CONTINUE
        KDO(2)=1
        KDO(6)=1
        KDO(8)=1
        KDO(11)=1
        KDO(12)=1
        KDO(13)=1
        KEPS(3)=1
        KEPS(4)=1
        KEPS(5)=1
        KEPS(9)=1
        NUL=6
        XNU(1) = 0.0D0
        XNU(2) = 0.6D0
        XNU(3) = 1.3D0
        XNU(4) = 2.0D0
        XNU(5) = 0.5D0*FNUL
        XNU(6) = FNUL + 1.2D0
      ENDIF
      I = 2
      EPS = 0.01D0
      FILM=DBLE(FLOAT(IL-1))
      T(1) = -PI + EPS
      DO 30 K=2,IL
        IF (KDO(K).EQ.0) THEN
          T(I) = PI*DBLE(FLOAT(-IL+2*K-1))/FILM
          IF (KEPS(K).EQ.0) GO TO 20
          TS=T(I)
          T(I) = TS - EPS
          I = I + 1
          T(I) = TS + EPS
   20     CONTINUE
          I = I + 1
        ENDIF
   30 CONTINUE
      ITL = I - 1
      LFLG = 0
      DO 190 KODE=1,2
        DO 180 N=1,NL
          DO 170 NU=1,NUL
            FNU = XNU(NU)
            IFNU = INT(SNGL(FNU))
            FFNU = FNU - DBLE(FLOAT(IFNU))
            DO 160 ICASE=1,3
              IRB = MIN(2,ICASE)
              DO 150 IR=IRB,4
                GO TO (50, 60, 70), ICASE
   50           CONTINUE
                R = (EPS*DBLE(FLOAT(4-IR))+2.0D0*DBLE(FLOAT(IR-1)))/
     *           3.0D0
                GO TO 80
   60           CONTINUE
                R = (2.0D0*DBLE(FLOAT(4-IR))+R2*DBLE(FLOAT(IR-1)))/3.0D0
                GO TO 80
   70           CONTINUE
                IF (RM.EQ.R2) GO TO 160
                R = (R2*DBLE(FLOAT(4-IR))+RM*DBLE(FLOAT(IR-1)))/3.0D0
   80           CONTINUE
                DO 140 IT=1,ITL
                  CT = COS(T(IT))
                  ST = SIN(T(IT))
                  IF (DABS(CT).LT.ATOL) CT = 0.0D0
                  IF (DABS(ST).LT.ATOL) ST = 0.0D0
                  ZR = R*CT
                  ZI = R*ST
                  CALL ZBESY(ZR, ZI, FNU, KODE, N, VR, VI, NZ2, CWRKR,
     &                       CWRKI, IERR)
                  IF (NZ2.NE.0.OR.IERR.NE.0) GO TO 140
                  CALL ZBESYH(ZR, ZI, FNU, KODE, N, WR, WI, NZ1, CWRKR,
     &                        CWRKI, IERR)
                  IF (NZ1.NE.0.OR.IERR.NE.0) GO TO 140
                  MFLG = 0
                  DO 120 I=1,N
                    AB = FNU+DBLE(FLOAT(I-1))
                    AA = MAX(0.5D0,AB)
                    CWR = WR(I) - VR(I)
                    CWI = WI(I) - VI(I)
                    AV = ZABS(VR(I),VI(I))
                    ER = ZABS(CWR,CWI)
                    IF (AV.NE.0.0D0) THEN
                      IF (ZI.EQ.0.0D0) THEN
                        IF (ZR.GT.0.0D0) THEN
                          IF (DABS(ZR).LT.AA) ER = ER/AV
                        ELSE
                          IF (DABS(FFNU-0.5D0).LT.0.125D0) THEN
                            IF (DABS(ZR).LT.AA) ER = ER/AV
                          ELSE
                            ER = ER/AV
                          ENDIF
                        ENDIF
                      ELSE
                        ER = ER/AV
                      ENDIF
                    ENDIF
                    AER(I) = ER
                    IF (ER.GT.ERTOL) MFLG = 1
  120             CONTINUE
                  IF (MFLG.EQ.0) GO TO 140
                  IF (LFLG.EQ.1) GO TO 130
                  WRITE (LUN,99995) ERTOL
99995             FORMAT (/' CASES WHICH VIOLATE THE RELATIVE ERROR TEST
     * WITH ERTOL = ', D12.4/)
                  WRITE (LUN,99994)
99994             FORMAT (/' OUTPUT FORMAT'/' KODE,N,IR,IT,NZ1,NZ2,ICASE
     *')
                  WRITE (LUN,99993)
99993             FORMAT (' ER(K),K=1,N'/' Z,FNU,W(KK),V(KK), KK=INDEX O
     *F FIRST NON-ZERO W,V PAIR'/)
                  LFLG = 1
  130             CONTINUE
                  KK = MAX(NZ1,NZ2) + 1
                  KK = MIN(N,KK)
                  WRITE (LUN,99992) KODE, N, IR, IT, NZ1, NZ2, ICASE
99992             FORMAT (8I5)
                  WRITE (LUN,99991) (AER(K),K=1,N)
                  WRITE (LUN,99991) ZR, ZI, FNU, WR(KK), WI(KK),
     *             VR(KK), VI(KK)
99991             FORMAT (7D12.4)
  140           CONTINUE
  150         CONTINUE
  160       CONTINUE
  170     CONTINUE
  180   CONTINUE
  190 CONTINUE
      IF (LFLG.EQ.0) WRITE (LUN,99990)
99990 FORMAT (/' QUICK CHECKS OK'/)
      STOP
      END
      SUBROUTINE ZBESYH(ZR, ZI, FNU, KODE, N, CYR, CYI, NZ, CWRKR,
     *           CWRKI, IERR)
C***BEGIN PROLOGUE  ZBESYH
C***DATE WRITTEN   830501   (YYMMDD)
C***REVISION DATE  890801, 930101   (YYMMDD)
C***CATEGORY NO.  B5K
C***KEYWORDS  Y-BESSEL FUNCTION,BESSEL FUNCTION OF COMPLEX ARGUMENT,
C             BESSEL FUNCTION OF SECOND KIND
C***AUTHOR  AMOS, DONALD E., SANDIA NATIONAL LABORATORIES
C***PURPOSE  TO COMPUTE THE Y-BESSEL FUNCTION OF A COMPLEX ARGUMENT
C***DESCRIPTION
C
C                ***A DOUBLE PRECISION ROUTINE***
C
C         ON KODE=1, ZBESYH COMPUTES AN N MEMBER SEQUENCE OF COMPLEX
C         BESSEL FUNCTIONS CY(I)=Y(FNU+I-1,Z) FOR REAL, NONNEGATIVE
C         ORDERS FNU+I-1, I=1,...,N AND COMPLEX Z IN THE CUT PLANE
C         -PI.LT.ARG(Z).LE.PI. ON KODE=2, ZBESYH RETURNS THE SCALED
C         FUNCTIONS
C
C         CY(I)=EXP(-ABS(Y))*Y(FNU+I-1,Z)   I = 1,...,N , Y=AIMAG(Z)
C
C         WHICH REMOVE THE EXPONENTIAL GROWTH IN BOTH THE UPPER AND
C         LOWER HALF PLANES FOR Z TO INFINITY. DEFINITIONS AND NOTATION
C         ARE FOUND IN THE NBS HANDBOOK OF MATHEMATICAL FUNCTIONS
C         (REF. 1).
C
C         INPUT      ZR,ZI,FNU ARE DOUBLE PRECISION
C           ZR,ZI  - Z=CMPLX(ZR,ZI), Z.NE.CMPLX(0.0D0,0.0D0),
C                    -PI.LT.ARG(Z).LE.PI
C           FNU    - ORDER OF INITIAL Y FUNCTION, FNU.GE.0.0D0
C           KODE   - A PARAMETER TO INDICATE THE SCALING OPTION
C                    KODE= 1  RETURNS
C                             CY(I)=Y(FNU+I-1,Z), I=1,...,N
C                        = 2  RETURNS
C                             CY(I)=Y(FNU+I-1,Z)*EXP(-ABS(Y)), I=1,...,N
C                             WHERE Y=AIMAG(Z)
C           N      - NUMBER OF MEMBERS OF THE SEQUENCE, N.GE.1
C           CWRKR, - DOUBLE PRECISION WORK VECTORS OF DIMENSION AT
C           CWRKI    AT LEAST N
C
C         OUTPUT     CYR,CYI ARE DOUBLE PRECISION
C           CYR,CYI- DOUBLE PRECISION VECTORS WHOSE FIRST N COMPONENTS
C                    CONTAIN REAL AND IMAGINARY PARTS FOR THE SEQUENCE
C                    CY(I)=Y(FNU+I-1,Z)  OR
C                    CY(I)=Y(FNU+I-1,Z)*EXP(-ABS(Y))  I=1,...,N
C                    DEPENDING ON KODE.
C           NZ     - NZ=0 , A NORMAL RETURN
C                    NZ.GT.0 , NZ COMPONENTS OF CY SET TO ZERO DUE TO
C                    UNDERFLOW (GENERALLY ON KODE=2)
C           IERR   - ERROR FLAG
C                    IERR=0, NORMAL RETURN - COMPUTATION COMPLETED
C                    IERR=1, INPUT ERROR   - NO COMPUTATION
C                    IERR=2, OVERFLOW      - NO COMPUTATION, FNU IS
C                            TOO LARGE OR CABS(Z) IS TOO SMALL OR BOTH
C                    IERR=3, CABS(Z) OR FNU+N-1 LARGE - COMPUTATION DONE
C                            BUT LOSSES OF SIGNIFCANCE BY ARGUMENT
C                            REDUCTION PRODUCE LESS THAN HALF OF MACHINE
C                            ACCURACY
C                    IERR=4, CABS(Z) OR FNU+N-1 TOO LARGE - NO COMPUTA-
C                            TION BECAUSE OF COMPLETE LOSSES OF SIGNIFI-
C                            CANCE BY ARGUMENT REDUCTION
C                    IERR=5, ERROR              - NO COMPUTATION,
C                            ALGORITHM TERMINATION CONDITION NOT MET
C
C***LONG DESCRIPTION
C
C         THE COMPUTATION IS CARRIED OUT BY THE FORMULA
C
C              Y(FNU,Z)=0.5*(H(1,FNU,Z)-H(2,FNU,Z))/I
C
C         WHERE I**2 = -1 AND THE HANKEL BESSEL FUNCTIONS H(1,FNU,Z)
C         AND H(2,FNU,Z) ARE CALCULATED IN ZBESH.
C
C         FOR NEGATIVE ORDERS,THE FORMULA
C
C              Y(-FNU,Z) = Y(FNU,Z)*COS(PI*FNU) + J(FNU,Z)*SIN(PI*FNU)
C
C         CAN BE USED. HOWEVER,FOR LARGE ORDERS CLOSE TO HALF ODD
C         INTEGERS THE FUNCTION CHANGES RADICALLY. WHEN FNU IS A LARGE
C         POSITIVE HALF ODD INTEGER,THE MAGNITUDE OF Y(-FNU,Z)=J(FNU,Z)*
C         SIN(PI*FNU) IS A LARGE NEGATIVE POWER OF TEN. BUT WHEN FNU IS
C         NOT A HALF ODD INTEGER, Y(FNU,Z) DOMINATES IN MAGNITUDE WITH A
C         LARGE POSITIVE POWER OF TEN AND THE MOST THAT THE SECOND TERM
C         CAN BE REDUCED IS BY UNIT ROUNDOFF FROM THE COEFFICIENT. THUS,
C         WIDE CHANGES CAN OCCUR WITHIN UNIT ROUNDOFF OF A LARGE HALF
C         ODD INTEGER. HERE, LARGE MEANS FNU.GT.CABS(Z).
C
C         IN MOST COMPLEX VARIABLE COMPUTATION, ONE MUST EVALUATE ELE-
C         MENTARY FUNCTIONS. WHEN THE MAGNITUDE OF Z OR FNU+N-1 IS
C         LARGE, LOSSES OF SIGNIFICANCE BY ARGUMENT REDUCTION OCCUR.
C         CONSEQUENTLY, IF EITHER ONE EXCEEDS U1=SQRT(0.5/UR), THEN
C         LOSSES EXCEEDING HALF PRECISION ARE LIKELY AND AN ERROR FLAG
C         IERR=3 IS TRIGGERED WHERE UR=DMAX1(D1MACH(4),1.0D-18) IS
C         DOUBLE PRECISION UNIT ROUNDOFF LIMITED TO 18 DIGITS PRECISION.
C         IF EITHER IS LARGER THAN U2=0.5/UR, THEN ALL SIGNIFICANCE IS
C         LOST AND IERR=4. IN ORDER TO USE THE INT FUNCTION, ARGUMENTS
C         MUST BE FURTHER RESTRICTED NOT TO EXCEED THE LARGEST MACHINE
C         INTEGER, U3=I1MACH(9). THUS, THE MAGNITUDE OF Z AND FNU+N-1 IS
C         RESTRICTED BY MIN(U2,U3). ON 32 BIT MACHINES, U1,U2, AND U3
C         ARE APPROXIMATELY 2.0E+3, 4.2E+6, 2.1E+9 IN SINGLE PRECISION
C         ARITHMETIC AND 1.3E+8, 1.8E+16, 2.1E+9 IN DOUBLE PRECISION
C         ARITHMETIC RESPECTIVELY. THIS MAKES U2 AND U3 LIMITING IN
C         THEIR RESPECTIVE ARITHMETICS. THIS MEANS THAT ONE CAN EXPECT
C         TO RETAIN, IN THE WORST CASES ON 32 BIT MACHINES, NO DIGITS
C         IN SINGLE AND ONLY 7 DIGITS IN DOUBLE PRECISION ARITHMETIC.
C         SIMILAR CONSIDERATIONS HOLD FOR OTHER MACHINES.
C
C         THE APPROXIMATE RELATIVE ERROR IN THE MAGNITUDE OF A COMPLEX
C         BESSEL FUNCTION CAN BE EXPRESSED BY P*10**S WHERE P=MAX(UNIT
C         ROUNDOFF,1.0E-18) IS THE NOMINAL PRECISION AND 10**S REPRE-
C         SENTS THE INCREASE IN ERROR DUE TO ARGUMENT REDUCTION IN THE
C         ELEMENTARY FUNCTIONS. HERE, S=MAX(1,ABS(LOG10(CABS(Z))),
C         ABS(LOG10(FNU))) APPROXIMATELY (I.E. S=MAX(1,ABS(EXPONENT OF
C         CABS(Z),ABS(EXPONENT OF FNU)) ). HOWEVER, THE PHASE ANGLE MAY
C         HAVE ONLY ABSOLUTE ACCURACY. THIS IS MOST LIKELY TO OCCUR WHEN
C         ONE COMPONENT (IN ABSOLUTE VALUE) IS LARGER THAN THE OTHER BY
C         SEVERAL ORDERS OF MAGNITUDE. IF ONE COMPONENT IS 10**K LARGER
C         THAN THE OTHER, THEN ONE CAN EXPECT ONLY MAX(ABS(LOG10(P))-K,
C         0) SIGNIFICANT DIGITS; OR, STATED ANOTHER WAY, WHEN K EXCEEDS
C         THE EXPONENT OF P, NO SIGNIFICANT DIGITS REMAIN IN THE SMALLER
C         COMPONENT. HOWEVER, THE PHASE ANGLE RETAINS ABSOLUTE ACCURACY
C         BECAUSE, IN COMPLEX ARITHMETIC WITH PRECISION P, THE SMALLER
C         COMPONENT WILL NOT (AS A RULE) DECREASE BELOW P TIMES THE
C         MAGNITUDE OF THE LARGER COMPONENT. IN THESE EXTREME CASES,
C         THE PRINCIPAL PHASE ANGLE IS ON THE ORDER OF +P, -P, PI/2-P,
C         OR -PI/2+P.
C
C***REFERENCES  HANDBOOK OF MATHEMATICAL FUNCTIONS BY M. ABRAMOWITZ
C                 AND I. A. STEGUN, NBS AMS SERIES 55, U.S. DEPT. OF
C                 COMMERCE, 1955.
C
C               COMPUTATION OF BESSEL FUNCTIONS OF COMPLEX ARGUMENT
C                 BY D. E. AMOS, SAND83-0083, MAY, 1983.
C
C               COMPUTATION OF BESSEL FUNCTIONS OF COMPLEX ARGUMENT
C                 AND LARGE ORDER BY D. E. AMOS, SAND83-0643, MAY, 1983
C
C               A SUBROUTINE PACKAGE FOR BESSEL FUNCTIONS OF A COMPLEX
C                 ARGUMENT AND NONNEGATIVE ORDER BY D. E. AMOS, SAND85-
C                 1018, MAY, 1985
C
C               A PORTABLE PACKAGE FOR BESSEL FUNCTIONS OF A COMPLEX
C                 ARGUMENT AND NONNEGATIVE ORDER BY D. E. AMOS, TRANS.
C                 MATH. SOFTWARE, 12, NO. 3, SEPTEMBER 1986, PP 265-273.
C
C***ROUTINES CALLED  ZBESH,I1MACH,D1MACH
C***END PROLOGUE  ZBESYH
C
C     COMPLEX CWRK,CY,C1,C2,EX,HCI,Z,ZU,ZV
      DOUBLE PRECISION CWRKI, CWRKR, CYI, CYR, C1I, C1R, C2I, C2R,
     * ELIM, EXI, EXR, EY, FNU, HCII, STI, STR, TAY, ZI, ZR, DEXP,
     * D1MACH, ASCLE, RTOL, ATOL, AA, BB, TOL, R1M5
      INTEGER I, IERR, K, KODE, K1, K2, N, NZ, NZ1, NZ2, I1MACH
      DIMENSION CYR(N), CYI(N), CWRKR(N), CWRKI(N)
C***FIRST EXECUTABLE STATEMENT  ZBESYH
      IERR = 0
      NZ=0
      IF (ZR.EQ.0.0D0 .AND. ZI.EQ.0.0D0) IERR=1
      IF (FNU.LT.0.0D0) IERR=1
      IF (KODE.LT.1 .OR. KODE.GT.2) IERR=1
      IF (N.LT.1) IERR=1
      IF (IERR.NE.0) RETURN
      HCII = 0.5D0
      CALL ZBESH(ZR, ZI, FNU, KODE, 1, N, CYR, CYI, NZ1, IERR)
      IF (IERR.NE.0.AND.IERR.NE.3) GO TO 170
      CALL ZBESH(ZR, ZI, FNU, KODE, 2, N, CWRKR, CWRKI, NZ2, IERR)
      IF (IERR.NE.0.AND.IERR.NE.3) GO TO 170
      NZ = MIN0(NZ1,NZ2)
      IF (KODE.EQ.2) GO TO 60
      DO 50 I=1,N
        STR = CWRKR(I) - CYR(I)
        STI = CWRKI(I) - CYI(I)
        CYR(I) = -STI*HCII
        CYI(I) = STR*HCII
   50 CONTINUE
      RETURN
   60 CONTINUE
      TOL = DMAX1(D1MACH(4),1.0D-18)
      K1 = I1MACH(15)
      K2 = I1MACH(16)
      K = MIN0(IABS(K1),IABS(K2))
      R1M5 = D1MACH(5)
C-----------------------------------------------------------------------
C     ELIM IS THE APPROXIMATE EXPONENTIAL UNDER- AND OVERFLOW LIMIT
C-----------------------------------------------------------------------
      ELIM = 2.303D0*(DBLE(FLOAT(K))*R1M5-3.0D0)
      EXR = DCOS(ZR)
      EXI = DSIN(ZR)
      EY = 0.0D0
      TAY = DABS(ZI+ZI)
      IF (TAY.LT.ELIM) EY = DEXP(-TAY)
      IF (ZI.LT.0.0D0) GO TO 90
      C1R = EXR*EY
      C1I = EXI*EY
      C2R = EXR
      C2I = -EXI
   70 CONTINUE
      NZ = 0
      RTOL = 1.0D0/TOL
      ASCLE = D1MACH(1)*RTOL*1.0D+3
      DO 80 I=1,N
C       CY(I) = HCI*(C2*CWRK(I)-C1*CY(I))
C       STR = C1R*CYR(I) - C1I*CYI(I)
C       STI = C1R*CYI(I) + C1I*CYR(I)
C       STR = -STR + C2R*CWRKR(I) - C2I*CWRKI(I)
C       STI = -STI + C2R*CWRKI(I) + C2I*CWRKR(I)
C       CYR(I) = -STI*HCII
C       CYI(I) = STR*HCII
        AA = CWRKR(I)
        BB = CWRKI(I)
        ATOL = 1.0D0
        IF (DMAX1(DABS(AA),DABS(BB)).GT.ASCLE) GO TO 75
          AA = AA*RTOL
          BB = BB*RTOL
          ATOL = TOL
   75   CONTINUE
        STR = (AA*C2R - BB*C2I)*ATOL
        STI = (AA*C2I + BB*C2R)*ATOL
        AA = CYR(I)
        BB = CYI(I)
        ATOL = 1.0D0
        IF (DMAX1(DABS(AA),DABS(BB)).GT.ASCLE) GO TO 85
          AA = AA*RTOL
          BB = BB*RTOL
          ATOL = TOL
   85   CONTINUE
        STR = STR - (AA*C1R - BB*C1I)*ATOL
        STI = STI - (AA*C1I + BB*C1R)*ATOL
        CYR(I) = -STI*HCII
        CYI(I) =  STR*HCII
        IF (STR.EQ.0.0D0 .AND. STI.EQ.0.0D0 .AND. EY.EQ.0.0D0) NZ = NZ
     *   + 1
   80 CONTINUE
      RETURN
   90 CONTINUE
      C1R = EXR
      C1I = EXI
      C2R = EXR*EY
      C2I = -EXI*EY
      GO TO 70
  170 CONTINUE
      NZ = 0
      RETURN
      END
      PROGRAM ZQCAI
C
C***DATE WRITTEN   830501   (YYMMDD)
C***REVISION DATE  890801, 930101   (YYMMDD)
C
C                *** A DOUBLE PRECISION ROUTINE ***
C
C     ZQCAI IS A QUICK CHECK ROUTINE FOR THE COMPLEX AIRY FUNCTIONS
C     GENERATED BY SUBROUTINES ZAIRY AND ZBIRY.
C
C     ZQCAI GENERATES AIRY FUNCTIONS AND THEIR DERIVATIVES FROM ZAIRY
C     AND ZBIRY AND CHECKS THEM AGAINST THE WRONSKIAN EVALUATION IN THE
C     REGION -PI/3 .LE. ARG(Z) .LE. PI/3:
C
C                 AI(Z)*BI'(Z)-AI'(Z)*BI(Z)=1/PI.
C
C     IN THE REMAINDER OF THE CUT PLANE, THE IDENTITIES
C
C              AI(Z)  = SQRT(-Z)*( J(-1/3,ZR) + J(1/3,ZR) )/3
C
C              AI'(Z) =        Z*( J(-2/3,ZR) - J(2/3,ZR) )/3
C
C       BI(Z)  =   I*SQRT(-Z/3)*( C1*H(1/3,1,ZR) - C2*H(1/3,2,ZR) )/2
C
C       BI'(Z) = I*(-Z)/SQRT(3)*( C2*H(2/3,1,ZR) - C1*H(2/3,2,ZR) )/2
C
C     ARE CHECKED WHERE ZR = (2/3)(-Z)**(3/2) WITH C1 = EXP(PI*I/6),
C     C2 = CONJG(C1) AND I**2 = -1.
C
C     THE PARAMETER MQC CAN HAVE VALUES 1 (THE DEFAULT) FOR A FASTER,
C     LESS DEFINITIVE TEST OR 2 FOR A SLOWER, MORE DEFINITIVE TEST.
C
C     MACHINE CONSTANTS ARE DEFINED IN FUNCTIONS I1MACH, R1MACH, AND
C     D1MACH. THESE MUST BE SELECTED BY THE USER OR SET ACCORDING TO
C     PROLOGUE INSTRUCTIONS.
C
C     COMPLEX CA, CAV, CHI, CI, CONA, CONB, CONC, COND, CON1, CON2, CV,
C    *  CW, CY, W, Y, YY, Z, ZR, ZW, SC
      EXTERNAL ZABS
      DOUBLE PRECISION AA, ALIM, ATOL, CAR, CAI, CVR, CVI,
     * CON1I, CON1R, CON2I, CON2R, CONAI, CONAR, CONBI, CONBR, CONCI,
     * CONCR, CONDI, CONDR, C13, C23, C43, CAVR, CAVI, CIR, CII, CHIR,
     * CHII, CWI, CWR, CYI, CYR, CT, DIG, ELIM, EPS, ER, ERTOL, FNUL,
     * FPI, HPI, PI, PTR, R, RL, RPI, R1M4, R1M5, SPI, STI, STR, ST,
     * T, TOL, RM, WI, WR, YI, YR, YYR, YYI, ZI, ZR, ZRI, ZRR, ZWR, ZWI,
     * D1MACH, ZABS, PI3, SLAK, FILM, AB, TS, SCR, SCI
      INTEGER I, ICASE, IL, IR, IRB, IRSET, IT, ITL, K, KEPS, KODE, J,
     * JB, JL, K1, K2, LFLG, LUN, NZ, IERR, I1MACH, MQC, ICL, KDO
      DIMENSION ER(5), T(20), YR(20), YI(20), YYR(20), YYI(20), WR(20),
     * WI(20), KEPS(20), KDO(20)
      DATA LUN /7/
      PARAMETER (MQC=1)
      OPEN(LUN,FILE='ZQCAI.OUT')
C-----------------------------------------------------------------------
C     SET PARAMETERS RELATED TO MACHINE CONSTANTS.
C     TOL IS THE APPROXIMATE UNIT ROUNDOFF LIMITED TO 1.0E-18.
C     ELIM IS THE APPROXIMATE EXPONENTIAL OVER- AND UNDERFLOW LIMIT.
C     EXP(-ELIM).LT.EXP(-ALIM)=EXP(-ELIM)/TOL    AND
C     EXP(ELIM).GT.EXP(ALIM)=EXP(ELIM)*TOL       ARE INTERVALS NEAR
C     UNDERFLOW AND OVERFLOW LIMITS WHERE SCALED ARITHMETIC IS DONE.
C     RL IS THE LOWER BOUNDARY OF THE ASYMPTOTIC EXPANSION FOR LARGE Z.
C     DIG = NUMBER OF BASE 10 DIGITS IN TOL = 10**(-DIG).
C     FNUL IS THE LOWER BOUNDARY OF THE ASYMPTOTIC SERIES FOR LARGE FNU
C-----------------------------------------------------------------------
      R1M4 = D1MACH(4)
      TOL = MAX(R1M4,1.0D-18)
      AA = -DLOG10(R1M4)
      K1 = I1MACH(15)
      K2 = I1MACH(16)
      R1M5 = D1MACH(5)
      K = MIN(ABS(K1),ABS(K2))
      ELIM = 2.303D0*(DBLE(FLOAT(K))*R1M5-3.0D0)
      AB = AA*2.303D0
      ALIM = ELIM + MAX(-AB,-41.45D0)
      DIG = MIN(AA,18.0D0)
      FNUL = 10.0D0 + 6.0D0*(DIG-3.0D0)
      RL = 1.2D0*DIG + 3.0D0
      SLAK = 3.0D0+4.0D0*(-DLOG10(TOL)-7.0D0)/11.0D0
      SLAK = MAX(SLAK,3.0D0)
      ERTOL = TOL*10.0D0**SLAK
      RM = 0.5D0*(ALIM + ELIM)
      RM=MIN(RM,200.0D0)
      RM=MAX(RM,RL+10.0D0)
C-----------------------------------------------------------------------
      WRITE (LUN,99999)
99999 FORMAT (' QUICK CHECK ROUTINE FOR THE AIRY FUNCTIONS FROM ZAIRY AN
     *D ZBIRY'/)
      WRITE (LUN,99998)
99998 FORMAT (' PARAMETERS TOL,ELIM,ALIM,RL,FNUL,DIG')
      WRITE (LUN,99997) TOL, ELIM, ALIM, RL, FNUL, DIG
99997 FORMAT (6D12.4/)
      ATOL = 100.0D0*TOL
      FPI = DATAN(1.0D0)
      HPI = FPI + FPI
      PI = HPI + HPI
      RPI = 1.0D0/PI
      SPI = PI/6.0D0
      CON1R = COS(SPI)
      CON1I = SIN(SPI)
      CON2R = CON1R
      CON2I = -CON1I
      PI3 = SPI+SPI
      C13 = 1.0D0/3.0D0
      C23 = C13+C13
      C43 = C23+C23
      CAVR = SQRT(C13)
      CAVI = 0.0D0
      CHIR = 0.0D0
      CHII = 0.5D0
      CIR = 0.0D0
      CII = 1.0D0
      WRITE (LUN,99996) MQC
99996 FORMAT (/' CHECKS IN THE (Z,FNU) SPACE WITH MQC = ',I2/)
C-----------------------------------------------------------------------
C     TEST VALUES OF Z IN -PI.LT.ARG(Z).LE.PI
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C     KDO(K), K=1,IL  DETERMINES WHICH OF THE IL ANGLES IN -PI TO PI
C     ARE USE TO COMPUTE VALUES OF Z
C       KDO(K) = 0  MEANS THAT THE INDEX K WILL BE USED FOR ONE OR TWO
C                   VALUES OF Z, DEPENDING ON THE CHOICE OF KEPS(K)
C              = 1  MEANS THAT THE INDEX K AND THE CORRESPONDING ANGLE
C                   WILL BE SKIPPED
C     KEPS(K), K=1,IL DETERMINES WHICH OF THE ANGLES GET INCREMENTED
C     UP AND DOWN TO PUT VALUES OF Z IN REGIONS WHERE DIFFERENT
C     FORMULAE ARE USED.
C       KEPS(K) =0  MEANS THAT THE ANGLE WILL BE USED WITHOUT CHANGE
C               =1  MEANS THAT THE ANGLE WILL BE INCREMENTED UP AND
C                   DOWN BY EPS
C     THE ANGLES TO BE USED ARE STORED IN THE T(I) ARRAY, I=1,ITL
C-----------------------------------------------------------------------
      IF (MQC.NE.2) THEN
        ICL=1
        IL=5
        DO 5 I=1,IL
          KDO(I)=0
          KEPS(I)=0
    5   CONTINUE
      ELSE
        ICL=2
        IL=7
        DO 6 I=1,IL
          KDO(I)=0
          KEPS(I)=0
    6   CONTINUE
        KEPS(2)=1
        KEPS(3)=1
        KEPS(5)=1
        KEPS(6)=1
      ENDIF
      I = 2
      EPS=0.01D0
      FILM=DBLE(FLOAT(IL-1))
      T(1) = -PI + EPS
      DO 30 K=2,IL
        IF(KDO(K).EQ.0) THEN
          T(I) = PI*DBLE(FLOAT(-IL+2*K-1))/FILM
          IF (KEPS(K).EQ.0) GO TO 20
          TS=T(I)
          T(I) = TS - EPS
          I = I + 1
          T(I) = TS + EPS
   20     CONTINUE
          I = I + 1
        ENDIF
   30 CONTINUE
      ITL = I - 1
      LFLG = 0
      DO 180 ICASE=1,ICL
        DO 170 KODE=1,2
          DO 160 IRSET=1,3
            IRB = MIN(IRSET,2)
            DO 150 IR=IRB,4
              GO TO (40, 50, 60), IRSET
   40         CONTINUE
              R =(0.2D0*DBLE(FLOAT(4-IR))+2.0D0*DBLE(FLOAT(IR-1)))/3.0D0
              GO TO 70
   50         CONTINUE
              R = (2.0D0*DBLE(FLOAT(4-IR))+RL*DBLE(FLOAT(IR-1)))/3.0D0
              GO TO 70
   60         CONTINUE
              R = (RL*DBLE(FLOAT(4-IR))+RM*DBLE(FLOAT(IR-1)))/3.0D0
   70         CONTINUE
              DO 140 IT=1,ITL
                CT = COS(T(IT))
                ST = SIN(T(IT))
                IF (DABS(CT).LT.ATOL) CT = 0.0D0
                IF (DABS(ST).LT.ATOL) ST = 0.0D0
                ZR = R*CT
                ZI = R*ST
                IF(ABS(T(IT)).LE.PI3) THEN
C-----------------------------------------------------------------------
C     WRONSKIAN CHECK IN -PI/3.LT.ARG(Z).LT.PI/3, TEST #1
C-----------------------------------------------------------------------
                  CALL ZAIRY(ZR, ZI, 0, KODE, YR(1), YI(1), NZ, IERR)
                  IF (NZ.NE.0.OR.IERR.NE.0) GO TO 140
                  CALL ZAIRY(ZR, ZI, 1, KODE, YR(2), YI(2), NZ, IERR)
                  IF (NZ.NE.0.OR.IERR.NE.0) GO TO 140
                  CALL ZBIRY(ZR, ZI, 0, KODE, WR(1), WI(1), IERR)
                  CALL ZBIRY(ZR, ZI, 1, KODE, WR(2), WI(2), IERR)
                  CWR = YR(1)*WR(2)-YI(1)*WI(2)
                  CWI = YR(1)*WI(2)+YI(1)*WR(2)
                  CYR = YR(2)*WR(1)-YI(2)*WI(1)
                  CYI = YR(2)*WI(1)+YI(2)*WR(1)
                  CVR = RPI
                  CVI = 0.0D0
                  IF (KODE.EQ.2) THEN
                    CALL ZSQRT(ZR,ZI,CAR,CAI)
                    ZRR = (ZR*CAR-ZI*CAI)*C23
                    ZRI = (ZR*CAI+ZI*CAR)*C23
                    AA=ABS(ZRR)
                    CAR = ZRR - AA
                    CAI = ZRI
                    CALL ZEXP(CAR,CAI,STR,STI)
                    PTR = STR*CVR-STR*CVI
                    CVI = STR*CVI+STI*CVR
                    CVR = PTR
                  ENDIF
                  CYR = CWR - CYR - CVR
                  CYI = CWI - CYI - CVI
                  ER(1) = ZABS(CYR,CYI)/ZABS(CVR,CVI)
                  JB = 1
                  JL = 1
                ELSE
C-----------------------------------------------------------------------
C     CHECKS IN -PI.LT.ARG(Z).LT.-PI/3 AND PI/3.LT.ARG(Z).LE.PI
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C     CHECK AI    TEST #2
C-----------------------------------------------------------------------
                  CALL ZAIRY(ZR, ZI, 0, KODE, YR(2), YI(2), NZ, IERR)
                  IF (NZ.NE.0.OR.IERR.NE.0) GO TO 140
                  ZRR = -ZR
                  ZRI = -ZI
                  CALL ZSQRT(ZRR,ZRI,CVR,CVI)
                  PTR = (ZRR*CVR-ZRI*CVI)*C23
                  ZRI = (ZRR*CVI+ZRI*CVR)*C23
                  ZRR = PTR
                  CALL ZBESJ(ZRR,ZRI,C23,KODE,2,YYR,YYI,NZ,IERR)
                  IF (NZ.NE.0.OR.IERR.NE.0) GO TO 140
                  STR = YYR(1)*C43
                  STI = YYI(1)*C43
                  CALL ZDIV(STR,STI,ZRR,ZRI,STR,STI)
                  CYR = STR - YYR(2)
                  CYI = STI - YYI(2)
                  CAR = YYR(1)
                  CAI = YYI(1)
                  CALL ZBESJ(ZRR,ZRI,C13,KODE,2,YYR,YYI,NZ,IERR)
                  IF (NZ.NE.0.OR.IERR.NE.0) GO TO 140
                  IF (KODE.EQ.2) THEN
                    AB = ABS(ZRI)
                    CALL ZSQRT(ZR,ZI,CWR,CWI)
                    ZWR = (ZR*CWR-ZI*CWI)*C23
                    ZWI = (ZR*CWI+ZI*CWR)*C23
                    CWR = ZWR+AB
                    CWI = ZWI
                    CALL ZEXP(CWR,CWI,STR,STI)
                    CWR = STR
                    CWI = STI
                    STR = YYR(1)*CWR-YYI(1)*CWI
                    YYI(1) = YYR(1)*CWI+YYI(1)*CWR
                    YYR(1) = STR
                    STR = YYR(2)*CWR-YYI(2)*CWI
                    YYI(2) = YYR(2)*CWI+YYI(2)*CWR
                    YYR(2) = STR
                    STR = CYR*CWR-CYI*CWI
                    CYI = CYR*CWI+CYI*CWR
                    CYR = STR
                    STR = CAR*CWR-CAI*CWI
                    CAI = CAR*CWI+CAI*CWR
                    CAR = STR
                    SCR = CWR
                    SCI = CWI
                  ENDIF
                  CWR = CVR*C13
                  CWI = CVI*C13
                  WR(2) = CWR*(YYR(1)+CYR)-CWI*(YYI(1)+CYI)
                  WI(2) = CWR*(YYI(1)+CYI)+CWI*(YYR(1)+CYR)
                  STR = YR(2)-WR(2)
                  STI = YI(2)-WI(2)
                  ER(2) = ZABS(STR,STI)
                  IF (ZI.NE.0.0D0.OR.ZR.GE.0.0D0) THEN
                    ER(2) = ER(2)/ZABS(YR(2),YI(2))
                  ELSE
                    IF (KODE.EQ.2) THEN
                      ER(2) = ER(2)/ZABS(SCR,SCI)
                    ENDIF
                  ENDIF
C-----------------------------------------------------------------------
C     CHECK AI'   TEST #3
C-----------------------------------------------------------------------
                  STR = YYR(1)*C23
                  STI = YYI(1)*C23
                  CALL ZDIV(STR,STI,ZRR,ZRI,CYR,CYI)
                  CYR = CYR-YYR(2)-CAR
                  CYI = CYI-YYI(2)-CAI
                  WR(3) = (ZR*CYR-ZI*CYI)*C13
                  WI(3) = (ZR*CYI+ZI*CYR)*C13
                  CALL ZAIRY(ZR, ZI, 1, KODE, YR(3), YI(3), NZ, IERR)
                  STR = YR(3)-WR(3)
                  STI = YI(3)-WI(3)
                  ER(3) = ZABS(STR,STI)
                  IF (ZI.NE.0.0D0.OR.ZR.GE.0.0D0) THEN
                    ER(3) = ER(3)/ZABS(YR(3),YI(3))
                  ELSE
                    IF (KODE.EQ.2) THEN
                      ER(3) = ER(3)/ZABS(SCR,SCI)
                    ENDIF
                  ENDIF
C-----------------------------------------------------------------------
C     CHECK BI    TEST #4
C-----------------------------------------------------------------------
                  CALL ZBESH(ZRR,ZRI,C13,KODE,1,1,YR,YI,NZ,IERR)
                  IF (NZ.NE.0.OR.IERR.NE.0) GO TO 140
                  CALL ZBESH(ZRR,ZRI,C13,KODE,2,1,YYR,YYI,NZ,IERR)
                  IF (NZ.NE.0.OR.IERR.NE.0) GO TO 140
                  CONAR = CON1R
                  CONAI = CON1I
                  CONBR = CON2R
                  CONBI = CON2I
                  CONCR = CON2R
                  CONCI = CON2I
                  CONDR = CON1R
                  CONDI = CON1I
                  IF (KODE.EQ.2) THEN
                    AA = ABS(ZWR)
                    ZWR = CIR*ZRR-CII*ZRI-AA
                    ZWI = CIR*ZRI+CII*ZRR
                    CALL ZEXP(ZWR,ZWI,CWR,CWI)
                    STR = CONAR*CWR-CONAI*CWI
                    CONAI =CONAR*CWI+CONAI*CWR
                    CONAR = STR
                    STR = CONCR*CWR-CONCI*CWI
                    CONCI =CONCR*CWI+CONCI*CWR
                    CONCR = STR
                    ZWR = -(CIR*ZRR-CII*ZRI)-AA
                    ZWI = -(CIR*ZRI+CII*ZRR)
                    CALL ZEXP(ZWR,ZWI,CWR,CWI)
                    STR = CONBR*CWR-CONBI*CWI
                    CONBI =CONBR*CWI+CONBI*CWR
                    CONBR = STR
                    STR = CONDR*CWR-CONDI*CWI
                    CONDI =CONDR*CWI+CONDI*CWR
                    CONDR = STR
                    SCR = CWR
                    SCI = CWI
                  ENDIF
                  CWR = CONAR*YR(1)-CONAI*YI(1)
                  CWI = CONAR*YI(1)+CONAI*YR(1)
                  CWR = CWR - (CONBR*YYR(1)-CONBI*YYI(1))
                  CWI = CWI - (CONBR*YYI(1)+CONBI*YYR(1))
                  STR = CVR*CAVR-CVI*CAVI
                  STI = CVR*CAVI+CVI*CAVR
                  PTR = STR*CWR-STI*CWI
                  CWI = STR*CWI+STI*CWR
                  CWR = PTR
                  WR(4) = CWR*CHIR-CWI*CHII
                  WI(4) = CWR*CHII+CWI*CHIR
                  CALL ZBIRY(ZR,ZI,0,KODE,YR(4),YI(4),IERR)
                  STR = YR(4)-WR(4)
                  STI = YI(4)-WI(4)
                  ER(4) = ZABS(STR,STI)
                  IF (ZI.NE.0.0D0.OR.ZR.GE.0.0D0) THEN
                    ER(4) = ER(4)/ZABS(YR(4),YI(4))
                  ELSE
                    IF (KODE.EQ.2) THEN
                      ER(4) = ER(4)/ZABS(SCR,SCI)
                    ENDIF
                  ENDIF
C-----------------------------------------------------------------------
C     CHECK BI'   TEST #5
C-----------------------------------------------------------------------
                  CALL ZBESH(ZRR,ZRI,C23,KODE,1,1,YR,YI,NZ,IERR)
                  IF (NZ.NE.0.OR.IERR.NE.0) GO TO 140
                  CALL ZBESH(ZRR,ZRI,C23,KODE,2,1,YYR,YYI,NZ,IERR)
                  IF (NZ.NE.0.OR.IERR.NE.0) GO TO 140
                  CWR = CONCR*YR(1)-CONCI*YI(1)
                  CWI = CONCR*YI(1)+CONCI*YR(1)
                  CWR = CWR - (CONDR*YYR(1)-CONDI*YYI(1))
                  CWI = CWI - (CONDR*YYI(1)+CONDI*YYR(1))
                  STR = -(ZR*CAVR-ZI*CAVI)
                  STI = -(ZR*CAVI+ZI*CAVR)
                  PTR = STR*CWR-STI*CWI
                  CWI = STR*CWI+STI*CWR
                  CWR = PTR
                  WR(5) = CWR*CHIR-CWI*CHII
                  WI(5) = CWR*CHII+CWI*CHIR
                  CALL ZBIRY(ZR,ZI,1,KODE,YR(5),YI(5),IERR)
                  STR = YR(5)-WR(5)
                  STI = YI(5)-WI(5)
                  ER(5) = ZABS(STR,STI)
                  IF (ZI.NE.0.0D0.OR.ZR.GE.0.0D0) THEN
                    ER(5) = ER(5)/ZABS(YR(5),YI(5))
                  ELSE
                    IF (KODE.EQ.2) THEN
                      ER(5) = ER(5)/ZABS(SCR,SCI)
                    ENDIF
                  ENDIF
                  JB = 2
                  JL = 5
                ENDIF
                DO 190 J=JB,JL
                IF (ER(J).LT.ERTOL) GO TO 190
                IF (LFLG.EQ.1) GO TO 130
                WRITE (LUN,99995) ERTOL
99995           FORMAT (/' CASES WHICH VIOLATE THE RELATIVE ERROR TEST W
     *ITH ERTOL =', E12.4/)
                WRITE (LUN,99994)
99994           FORMAT (/' OUTPUT FORMAT'/' KODE,IR,IT,IRSET,ICASE')
                WRITE (LUN,99993)
99993           FORMAT (' ER'/' I, Z, Y(I), W(I), ON THE I-TH TEST, I=1,
     *5'/)
                LFLG = 1
  130           CONTINUE
                WRITE (LUN,99992) KODE, IR, IT, IRSET, ICASE
99992           FORMAT (5I5)
                WRITE (LUN,99991) ER(J)
                WRITE (LUN,99989) J, ZR, ZI, YR(J), YI(J), WR(J), WI(J)
99991           FORMAT (D12.4)
99989           FORMAT (I5,6D12.4)
  190           CONTINUE
  140         CONTINUE
  150       CONTINUE
  160     CONTINUE
  170   CONTINUE
  180 CONTINUE
      IF (LFLG.EQ.0) WRITE (LUN,99990)
99990 FORMAT (/' QUICK CHECKS OK'/)
      STOP
      END
