C
C     TITLE: TEST DRIVER PROGRAM FOR NEUMAN
C     -------------------------------------
C
C     This program is used to test the subroutine "NEUMAN" which
C     solves interior and exterior Neumann problems for Laplace's
C     equation on a planar region D.  
C
C     The subroutine NEUMAN solves exterior Neumann problems
C     by means of a standard indirect boundary integral equation
C     reformulation, one based on representing the solution as a
C     single layer potential.
C
C     For the interior Neumann problem, the problem is first
C     reformulated as an exterior Neumann problem by means of the
C     Kelvin transformation, and this new problem is solved as for
C     the exterior case.
C
C     This program is limited to problems on simply connected
C     regions D with a smooth boundary curve C.
C
C     * THE USER SUPPLIED SUBROUTINES:
C
C     FUNCTION BDYFCN  : This defines the boundary data, the exterior
C                        normal derivative of the unknown potential.
C                        For the present test program, BDYFCN is
C                        generated from GVFCN.
C     SUBROUTINE GVFCN : This contains a set of test functions, all
C                        potential functions.  This also includes their
C                        first order partial derivatives. The subroutine
C                        BDYFCN calls this routine.
C     SUBROUTINE CURVE : This defines a set of test curves.
C
C     ***** THE SUBROUTINES AND FUNCTIONS CALLED BY SUBROUTINE NEUMAN
C
C     SUBROUTINE EVALU : This subroutine evaluates the single layer
C                        potential by using the density function
C                        evaluated in subroutine INTEQN.
C     SUBROUTINE FUCOEF: This subroutine evaluates the Fourier
C                        coefficients of the density function by
C                        using the FFT subroutines.  This is used
C                        only if the solution is needed on the
C                        boundary.
C     FUNCTION  FUEVAL : This subroutine evaluates the singular part
C                        of the single layer potential when it is
C                        evaluated at boundary points.
C     SUBROUTINE INTEQN: This subroutine generates the density
C                        function.
C     SUBROUTINE NEWCUR: This subroutine defines a transformed curve
C                        when we solve the interior Neumann problem
C                        by calling subroutine KVTRNF.
C     SUBROUTINE KVTRNF: This defines the Kelvin transformation.
C     SUBROUTINES DGESV & DGECON: LAPACK subroutines for LU
C                        decomposition  and condition number
C                        evaluation.
C     FFT SUBROUTINES  : Modified version of selected FFTPACK
C                        subroutines for the Fast Fourier Transform.
C                        The package FFTPACK was written by Paul
C                        Swarztrauber of NCAR.
C     FUNCTION D1MACH  : Subroutine which defines the machine unit round.
C
C    ****************************************************************
C     FOLLOWING ARE THE PARAMETERS THAT MUST BE SUPPLIED BY USERS.
C     IN THIS TEST PROGRAM, THEY ARE REQUESTED INTERACTIVELY FROM
C     THE USER.
C
C     IE          This parameter should be 0 or 1.
C                 IE=0 for the interior Neumann problem.
C                 IE=1 for the exterior Neumann problem.
C     IDBG        The debug parameter. DBG=Y produces an output 
C                 with the debugging information. DBG=N produces a 
C                 standard answer file.
C     NUMCUR      The index of the boundary curve C in subroutine CURVE.
C                 NUMCUR=1 for an ellipse
C                 NUMCUR=2 for a limacon
C                 NUMCUR=3 for the ovals of Cassini
C     A,B         Parameters used in defining the curve C.
C     NUMF        The index of the test functions in subroutine GVFCN.
C                 Indices 1, 2, and 3 are test cases for solving the
C                 interior Neumann problem; and index 4 is for the
C                 test case for solving the exterior Neumann problem.
C     EPS         The absolute error tolerance on the solution of the
C                 Neumann problem.
C     DPTS        The two dimensional double precision array in which
C                 are stored the points at which the harmonic function
C                 is to be evaluated.  Each column denotes a new point
C                 at which the solution is to be evaluated.  For the
C                 definition, see below in the definition of DPTS
C                 given in this main program.  Also, see the
C                 introductory statements of NEUMAN.  See the example
C                 below in this test driver program for setting up
C                 DPTS, especially the distinction between specifying
C                 points on the boundary and points off the boundary.
C                 *** NOTE *** In this test program, our method of
C                 creating DPTS assumes that the region inside of
C                 the curve is starlike with respect to the origin.
C     MAXDPTS     This is the maximum number of points that can be
C                 stored in DPTS.  It is to be set in the PARAMETER
C                 statement given below.
C     NTHETA, NR  The parameters used in defining the test points in
C                 the array DPTS, at which the solution is evaluated.
C                 NTHETA will be the number of angular subdivisions,
C                 and NR is the number of radial subdivisions.
C                 Boundary points are included in this test case.
C     NEUMAN.ANS  The name of the output file.
C
C     The subroutine NEUMAN requires two working storage arrays,
C     WORK and IWORK of respective dimensions NWORK and MWORK.
C     These are set by users, as follows.
C
C     WORK        A double precision work array.  It dimension
C                 should be at least 10,000.  If the array DPTS
C                 contains points close to the boundary, then
C                 the dimension should be increased accordingly,
C                 to obtain accurate numerical integrations for
C                 the potential approximations at such points.
C                 A dimension of 300,000 will allow for very
C                 accurate potential evaluations, even near to
C                 to the boundary.  See the discussion of NWORK,
C                 the dimension of WORK, below in the introductory
C                 comments of SUBROUTINE NEUMAN.
C     IWORK       The integer array.  Its dimension should be at
C                 least
C                    MWORK = SQRT(NWORK + 49)+ 8
C                 where NWORK is the dimension of WORK.
C
      DOUBLE PRECISION A,B,TRUE,CLNGTH,D2X,D2Y,DX,DY,EPS,ERR,F,FNR,
     +        FX,FY,HTHETA,R,THETA,X,Y,FR
      INTEGER MWORK,NWORK,NUMCUR,NUMF,I,IBD,IBEG,IE,IER,II,ILOW,
     +        IOUT,IUP,J,NP,NPTS,NR,NTHETA,TML_IN,TML_OUT,I_FILE,
     +        INP, IDBG
      INTEGER MAXDPTS, MAXFFT
      DOUBLE PRECISION ONE
      PARAMETER (MAXDPTS=200,MAXFFT=1024)
      PARAMETER (MWORK=600,NWORK=300000)
      PARAMETER (ONE=1.0D0)
      DOUBLE PRECISION DPTS(4,MAXDPTS),ERROR(MAXDPTS),U(MAXDPTS),
     *                 WORK(NWORK)
      INTEGER IWORK(MWORK)
      CHARACTER TYN
      EXTERNAL BDYFCN,CURVE,GVFCN,NEUMAN 
      INTRINSIC ABS,ATAN,FLOAT
      COMMON/BLKCUR/A,B,NUMCUR
      COMMON/BLKF/NUMF
      COMMON/DUMMY/WORK
C
C     ************************************************************
C     The standard input and output unit numbers:
      DATA TML_IN/5/,TML_OUT/6/
C     The output unit number for results from this program, to be
C     stored in the file "neuman.ans":
      DATA I_FILE/8/
C     ************************************************************
C
C     INITIALIZATION
      CLNGTH = 8.0D0*ATAN(1.0D0)
      IBEG = 0
      INP = TML_IN
      PRINT *, 'DO YOU WANT THE OUTPUT DIRECTED TO THE TERMINAL? (Y/N)'
      READ(*,'(A1)') TYN
      IF((TYN .EQ. 'Y') .OR. (TYN .EQ. 'y')) THEN
          IOUT = TML_OUT
      ELSE
          IOUT = I_FILE
          OPEN(IOUT,FILE='res1')
      END IF
C
      PRINT *, 'DO YOU WANT A DEBUGGING OUTPUT? (Y/N)'
      READ(*,'(A1)') TYN
      IF((TYN .EQ. 'Y') .OR. (TYN .EQ. 'y')) THEN
          IDBG = 1
      ELSE
          IDBG = 0
      END IF
C
C     INPUT PROBLEM PARAMETERS.
   20 PRINT *,'CURVE DEFINING PARAMETERS; NUMCUR, A, B ?'
      PRINT *,'CHOOSE 1, WITH A,B > 0 FOR AN ELLIPSE'
      PRINT *,'CHOOSE 2, WITH 0 < A < 1, B > 0 FOR A LIMACON'
      PRINT *,'CHOOSE 3, WITH A > 1, B > 0 FOR A CASSINI'
      PRINT *,'CHOOSE 4, WITH ARBITRARY A, B FOR AN AMEOBA'
      READ (INP,*) NUMCUR,A,B
C
C     CHECK THE INPUT PARAMETERS.
      IF (NUMCUR .EQ. 1)  THEN
          IF (A .GT. 0. .AND. B .GT. 0.)  THEN
              GO TO 30
          ELSE
              GO TO 25
          END IF
      ELSE IF (NUMCUR .EQ. 2) THEN
          IF (0. .LT. A .AND. A .LT. 1. .AND. B .GT. 0.)  THEN 
              GO TO 30
          ELSE
              GO TO 25
          END IF
      ELSE IF (NUMCUR .EQ. 3) THEN
          IF (A .GT. 1. .AND. B .GT. 0.) THEN 
              GO TO 30
          ELSE
              GO TO 25
          END IF
      ELSE IF (NUMCUR .EQ. 4)  THEN
          GO TO 30
      ELSE
          GO TO 25
      END IF
   25 PRINT *, 'GIVE THE CURVE-DEFINING PARAMETERS AGAIN.'
      PRINT *, 'SOMETHING WAS WRONG WITH THE PREVIOUS VALUES.'
      GO TO 20
   30 PRINT *,'GIVE THE ERROR TOLERANCE EPS?'
      READ (INP,*) EPS
C
   40 PRINT *,'GIVE IE=0 FOR AN INTERIOR PROBLEM.'
      PRINT *,'GIVE IE=1 FOR AN EXTERIOR PROBLEM.'
      PRINT *,'WHAT IS IE?'
      READ (INP,*) IE
C
      PRINT *,'SPECIFY THE TEST FUNCTION INDEX NUMBF?,'
      PRINT *,'FOR IE=0, CHOOSE NUMF = 1, 2, OR 3'
      PRINT *,'FOR IE=1, CHOOSE NUMF = 4'
      READ (INP,*) NUMF
C
      IF(IE .EQ. 0 .AND. NUMF .GE. 1 .AND. NUMF .LE. 3) GO TO 50
      IF(IE .EQ. 1 .AND. NUMF .EQ. 4) GO TO 50
      PRINT *, 'GIVE IE AND NUMF, AGAIN.'
      GO TO 40
C
   50 PRINT *,'PARAMETERS N_THETA AND N_R ARE USED TO DEFINE A SET'
      PRINT *,'OF MESH POINTS INTERIOR TO D, AT WHICH THE POTENTIAL'
      PRINT *,'U IS TO BE EVALUATED.'
      PRINT *,'GIVE N_THETA AND N_R' 
      READ (INP,*) NTHETA,NR
C
      WRITE (IOUT,9000) A,B,EPS,NUMF,NUMCUR,NTHETA,NR,IE
C
C     SET UP POINTS (X,Y) AT WHICH POTENTIAL SOLUTION U IS
C     TO BE EVALUATED.
C     IE=0 ASSIGNS THE POINTS IN THE INTERIOR OF THE CURVE C.
C     IE=1 ASSIGNS THE POINTS IN THE EXTERIOR OF THE CURVE C.
C
C     IF (X,Y) IS A BOUNDARY POINT, THEN
C        DPTS(1,I)=X, DPTS(2,I)=Y, DPTS(3,I)=S, AND DPTS(4,I)=1.
C     HERE S IS A PARAMETER WHICH IS USED FOR THE PARAMETRIZATION
C     OF THE CURVE C.
C     IF (X,Y) IS A NON BOUNDARY POINT, THEN
C        DPTS(1,I)=X, DPTS(2,I)=Y, DPTS(4,I)=0
C     AND DPTS(3,I) NEED NOT BE SET.
C
      IF(IE .EQ. 0) THEN
C         IT IS AN INTERIOR PROBLEM, AND CHOOSE (0,0) AS THE
C         FIRST POINT.
          DPTS(1,1) = 0.D0
          DPTS(2,1) = 0.D0
          DPTS(4,1) = 0.D0
      ELSE
C         IT IS AN EXTERIOR PROBLEM, AND CHOOSE A POINT SUFFICIENTLY
C         AWAY FROM THE BOUNDARY FOR THE FIRST POINT.
          DPTS(1,1) = 10000.D0
          DPTS(2,1) = 10000.D0
          DPTS(4,1) = 0.D0
      END IF
C
C     THE FOLLOWING SETUP OF DPTS ASSUMES THAT THE INTERIOR OF C IS
C     STARLIKE WITH RESPECT TO THE ORIGIN (0,0).
      HTHETA = CLNGTH/NTHETA
      NP = 1
      FNR = NR
      ILOW = 1
      IUP = NTHETA
      DO II = ILOW,IUP
        I = II - 1
        THETA = I*HTHETA
        CALL CURVE(THETA,X,Y,DX,DY,D2X,D2Y)
        DO J = 1,NR
          NP = NP + 1
C         CHECK WHETHER A SUFFICIENT WORK STORAGE IS GIVEN
C         FOR THE ARRAY DPTS.
          IF(NP .GT. MAXDPTS) THEN
              PRINT *,' MORE WORK SPACE IS NEEDED FOR THE ARRAY DTPS.'
              PRINT *,' INCREASE MAXDPTS AND RECOMPILE THE PROGRAM.'
              STOP
          END IF
          FR = FLOAT(J)/FNR
          R=2*FR - FR*FR 
          IF(IE .EQ. 1) R = ONE/R
          DPTS(1,NP) = R*X
          DPTS(2,NP) = R*Y
          DPTS(3,NP) = THETA
          IF(J .EQ. NR) THEN
              DPTS(4,NP) = 1
          ELSE
              DPTS(4,NP) = 0
          END IF
        END DO
      END DO
      NPTS = NP
C
      CALL NEUMAN(IOUT,IDBG,IE,BDYFCN,CURVE,CLNGTH,DPTS,NPTS,IBEG,
     +            EPS,WORK,IWORK,NWORK,MWORK,MAXFFT,U,ERROR,IER)
C
C     PRINT RESULTS.
C
      WRITE (IOUT,9010) IER
      WRITE (IOUT,9020)
      WRITE (IOUT,9030)
C     COMPARE THE GIVEN TEST FUNCTION AND THE APPROXIMATING SOLUTION.
      DO I = 1,NPTS
          X = DPTS(1,I)
          Y = DPTS(2,I)
          IBD = DPTS(4,I)
          CALL GVFCN(X,Y,F,FX,FY)
          TRUE = F
          ERR = ABS(TRUE-U(I))
          WRITE (IOUT,9050) IBD,X,Y,U(I),ERR,ERROR(I)
      END DO
C
      PRINT *, ' '
      PRINT *, 'THE CALCULATION IS COMPLETE.  DO YOU WISH TO'
      PRINT *, 'CONTINUE WITH ANOTHER CALCULATION? (Y/N)'
      READ(INP,'(A1)') TYN
      IF((TYN .EQ. 'Y') .OR. (TYN .EQ. 'y')) THEN
          PRINT *, 'DO WANT TO USE THE SAME CURVE AND BOUNDARY'
          PRINT *, 'FUNCTION AGAIN, BUT WITH DIFFERENT VALUES OF'
          PRINT *, 'NR AND NTHETA? (Y/N)'
          READ(INP,'(A1)') TYN
          IF((TYN .EQ. 'Y') .OR. (TYN .EQ. 'y')) THEN
              IBEG = 1
              GO TO 50
          ELSE
              IBEG = 0
              GO TO 20
          END IF
      ELSE
          IF (IOUT .EQ. I_FILE) CLOSE(IOUT)
          STOP
      END IF
C
 9000 FORMAT (/,/,' BOUNDARY PARAMETERS: A=',F7.4,3X,'B=',F7.4,/,
     +       ' EPS=',1PD11.3,5X,'BOUNDARY FUNCTION ',I1,5X,'CURVE ',
     +       I1,/,' NTHETA=',I3,3X,'NR=',I3,3X,'IE=',I3,/)
 9010 FORMAT (' IER=',I2)
 9020 FORMAT (/,'In the following table, ERROR is the true error;',/,
     +    'The magnitude of PREDERR is the predicted error bound;',/,
     +    'and if PREDERR is negative, the desired error tolerance',/,
     +    'is predicted to have not been attained.  The variable',/,
     +    'IBD denotes if (X,Y) is inside D [IBD=0] or on the',/,
     +    'boundary [IBD=1].',/)
 9030 FORMAT(/,' IBD',7X,'X',11X,'Y',10X,'U(X,Y)',11X,'ERROR',6X,
     +       'PREDERR',/)
 9050 FORMAT (I3,F11.4,F12.4,1PD20.10,D12.2,D12.2)
      END
C
      DOUBLE PRECISION FUNCTION BDYFCN(S)
C     --------------------------------
C
C     Define the Neumann data on the boundary, using subroutine GVFCN,
C     which defines the test harmonic functions.
C
      DOUBLE PRECISION D2X,D2Y,DX,DY,F,FX,FY,S,X,Y
      EXTERNAL CURVE,GVFCN
      INTRINSIC SQRT
C
      CALL CURVE(S,X,Y,DX,DY,D2X,D2Y)
      CALL GVFCN(X,Y,F,FX,FY)
      BDYFCN = (FX*DY-FY*DX)/SQRT(DX*DX+DY*DY)
      RETURN
      END
C
      SUBROUTINE GVFCN(X,Y,F,FX,FY)
C     ----------------
C
C     This program defines the test functions and their first order
C     partial derivatives.
C
      DOUBLE PRECISION F,FX,FY,X,Y
      INTEGER NUMF
      INTRINSIC COS,EXP,SIN
      COMMON /BLKF/NUMF
C
      GO TO (10,20,30,40) NUMF
   10 F = X**2 - Y**2
      FX = 2.D0*X
      FY = -2.D0*Y
      RETURN
   20 F = X
      FX = 1.D0
      FY = 0.D0
      RETURN
   30 F = EXP(X)*COS(Y) - 1.D0
      FX = EXP(X)*COS(Y)
      FY = -EXP(X)*SIN(Y)
      RETURN
   40 F = X/ (X*X+Y*Y)
      FX = (Y*Y-X*X)/ (X*X+Y*Y)**2
      FY = -2.D0*X*Y/ (X*X+Y*Y)**2
      RETURN
      END

      SUBROUTINE CURVE(S,X,Y,DX,DY,D2X,D2Y)
C     ----------------
C
C     This program defines the boundary of the region C, along
C     with its first and second derivatives with respect to the
C     parameterization variables.  The curves are an ellipse,
C     limacon, and the Ovals of Cassini.
C
      DOUBLE PRECISION A,B,D2X,D2Y,DX,DY,S,X,Y
      INTEGER NUMCUR
      DOUBLE PRECISION CS,CS2,D2R,DR,EC,ES,R,SN,SN2
      INTRINSIC COS,EXP,SIN,SQRT
      COMMON /BLKCUR/A,B,NUMCUR
      GO TO (10,20,40,50) NUMCUR
C
C     DEFINE AN ELLIPSE.
   10 CS = COS(S)
      SN = SIN(S)
      X = A*CS
      Y = B*SN
      DX = -A*SN
      DY = B*CS
      D2X = -X
      D2Y = -Y
      RETURN
C
C     DEFINE A LIMACON.
C     CHOOSE 0 .LE. A .LT. 1,  0 .LT. B.
C     GRAPH CENTERED AT ORIGIN, BETWEEN X=-1 AND X=1,
C     SYMMETRIC ABOUT THE X-AXIS.
   20 CS = COS(S)
      SN = SIN(S)
      R = 1.0 + A*CS
      DR = -A*SN
      D2R = -A*CS
      X = R*CS - A
      Y = B*R*SN
   30 DY = B* (DR*SN+R*CS)
      D2Y = B* (SN* (D2R-R)+2.0D0*DR*CS)
      DX = DR*CS - R*SN
      D2X = (D2R-R)*CS - 2.0D0*DR*SN
      RETURN
C
C     DEFINE THE OVALS OF CASSINI.
C     CHOOSE A .GT. 1, B .GT. 0.
   40 CS = COS(S)
      SN = SIN(S)
      CS2 = CS*CS - SN*SN
      SN2 = 2.0D0*SN*CS
      R = SQRT(CS2+SQRT(A-SN2*SN2))
      DR = -R*SN2/ (R*R-CS2)
      D2R = - (2.0D0*CS2*R+ (2.0D0*R*DR+3.0D0*SN2)*DR)/ (R*R-CS2)
      X = R*CS
      Y = B*R*SN
      GO TO 30
C
C     DEFINE AN "AMOEBA" BOUNDARY.
   50 CS = COS(S)
      SN = SIN(S)
      CS2 = CS*CS - SN*SN
      SN2 = 2.0D0*SN*CS
      EC = EXP(CS)
      ES = EXP(SN)
      R = EC*CS2*CS2 + ES*SN2*SN2
      X = R*CS
      Y = R*SN
      DR = -EC*(SN*CS2*CS2 + 4.0D0*CS2*SN2)
     *     +ES*(4.0D0*CS2*SN2 + CS*SN2*SN2)
      D2R = EC*(-8*CS2*CS2 - CS*CS2*CS2 + CS2*CS2*SN*SN
     *        + 8*CS2*SN*SN2 + 8*SN2*SN2)
     *     +ES*(8*CS2*CS2 +8*CS*CS2*SN2 - 8*SN2*SN2
     *        + CS*CS*SN2*SN2 - SN*SN2*SN2)
      DY = DR*SN+R*CS
      D2Y = SN* (D2R-R)+2.0D0*DR*CS
      DX = DR*CS - R*SN
      D2X = (D2R-R)*CS - 2.0D0*DR*SN
      RETURN
      END
C
C**************************************************************
C     THIS IS THE END OF THE USER SUPPLIED INFORMATION.
C**************************************************************
