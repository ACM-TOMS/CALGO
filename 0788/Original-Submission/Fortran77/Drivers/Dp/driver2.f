C     TITLE: TEST DRIVER PROGRAM FOR DRCHLT
C     -------------------------------------
C
C     This program is used to test the subroutine "DRCHLT" which
C     solves interior and exterior Dirichlet problems for Laplace's
C     equation on a planar region D.  
C
C     The subroutine DRCHLT solves interior Dirichlet problems
C     by means of a standard indirect boundary integral equation
C     reformulation, one based on representing the solution as a
C     double layer potential.
C
C     For the exterior Dirichlet problem, the problem is first
C     reformulated as an interior Dirichlet problem by means of the
C     Kelvin transformation, and this new problem is solved as for
C     the interior case.
C
C     This program is limited to problems on simply connected
C     regions D with a smooth boundary curve C.
C
C     ***** THE USER SUPPLIED SUBROUTINES:
C
C     FUNCTION BDYFCN  : This defines the boundary data.
C     SUBROUTINE CURVE : This defines a set of test curves.
C
C     ***** THE SUBROUTINES INCLUDED IN THE SUBROUTINE DRCHLT
C
C     SUBROUTINE EVALU : This subroutine evaluates the double layer
C                        potential by using the density function
C                        evaluated in subroutine INTEQN.
C     SUBROUTINE INTEQN: This subroutine generates the density
C                        function.
C     SUBROUTINE NEWCUR: This subroutine defines a transformed curve
C                        when we solve the exterior Dirichlet problem
C                        by calling subroutine KVTRNF.
C     SUBROUTINE KVTRNF: This defines the Kelvin transformation.
C     SUBROUTINES DGESV & DGECON: LAPACK subroutines for LU
C                        decomposition  and condition number
C                        evaluation.
C     FUNCTION D1MACH  : Subroutine which defines the machine unit
C                        round.
C
C     *************************************************************
C     FOLLOWING ARE THE PARAMETERS THAT MUST BE SUPPLIED BY USERS.
C     IN THIS TEST PROGRAM, THEY ARE REQUESTED INTERACTIVELY FROM
C     THE USER.
C
c     IE          This parameter should be 0 or 1.
C                 IE=0 for the interior Dirichlet problem.
C                 IE=1 for the exterior Dirichlet problem.
C     IDBG        The debug parameter. DBG=Y produces an output 
C                 with the debugging information. DBG=N produces a 
C                 shortcut answer file.   
C     NUMCUR      The index of the boundary curve C in subroutine CURVE.
C                 NUMCUR=1 for an ellipse
C                 NUMCUR=2 for a limacon
C                 NUMCUR=3 for the ovals of Cassini
C     A,B         Parameters used in defining the curve C.
C     NUMBF       The index of the test functions in BDYFCN.
C                 Indices 1, 2, and 3 are test cases for solving the
C                 interior Dirichlet problem; and indices 4 AND 5
C                 are test cases for solving the exterior Dirichlet
C                 problem.
C     EPS         The absolute error tolerance on the solution of the
C                 Dirichlet problem.
C     R_FORM      This specifies (indirectly) the form of error test to
C                 to be used in subroutines INTEQN and EVALU.
C                 If R_FORM=0, then we use a "normal" error test.  This
C                 attempts to measure the rate of convergence of the 
C                 approximates and to use that in predicting the error.
C                 If R_FORM=1, then we assume the approximates are 
C                 converging at a very slow rate.  Use this with more
C                 ill-behaved solution functions and boundaries.
C     NTHETA, NR  The parameters used in defining the test points in
C                 the array DPTS, at which the solution is evaluated.
C                 NTHETA will be the number of angular subdivisions,
C                 and NR is the number of radial subdivisions
C     DRCHLT.ANS  The name of the output file.      
C
C     The subroutine DRCHLT requires two working storage arrays,
C     WORK and IWORK of respective dimensions NWORK and MWORK.
C     These are set by users, as follows.
C
C     WORK        A double precision work array.  It dimension
C                 should be at least 5,000.  If the array DPTS
C                 contains points close to the boundary, then
C                 the dimension should be increased accordingly,
C                 to obtain accurate numerical integrations for
C                 the potential approximations at such points.
C                 A dimension of 300,000 will allow for very
C                 accurate potential evaluations, even near to
C                 to the boundary.  See the discussion of NWORK,
C                 the dimension of WORK, below in the introductory
C                 comments of SUBROUTINE DRCHLT.
C     IWORK       The integer array.  Its dimension should be at
C                 least
C                    MWORK = SQRT(NWORK + 36)-6
C                 where NWORK is the dimension of WORK.
C
C
      INTEGER NUMBF,NUMCUR,NWORK,MWORK,TML_IN,TML_OUT,IDBG
      INTEGER I,IBEG,IE,IER,II,ILOW,IOUT,IUP,J,NP,NPTS,NR,NTHETA,
     *    R_FORM,I_FILE,MAXDPTS
      DOUBLE PRECISION A,ANS,B,BDYFCN,CLNGTH,D2X,D2Y,DX,DY,EPS,ERR,
     +       FNR,HTHETA,R,THETA,X,Y
      CHARACTER TYN, RFC*12
      PARAMETER(MAXDPTS=200)
      PARAMETER(NWORK=300000, MWORK=600)
      DOUBLE PRECISION DPTS(2,MAXDPTS),ERROR(MAXDPTS),U(MAXDPTS),
     +      WORK(NWORK)
      INTEGER IWORK(MWORK)
      EXTERNAL BDYFCN,CURVE,DRCHLT
      INTRINSIC ATAN,FLOAT
      COMMON /BLKBF/NUMBF
      COMMON /BLKCUR/A,B,NUMCUR
      COMMON /DUMMY/WORK
C     ************************************************************
C     The standard input and output unit numbers:
      DATA TML_IN/5/,TML_OUT/6/
C     The output unit number for results from this program, to be
C     stored in the file "drchlt.ans":
      DATA I_FILE/8/
C     ************************************************************
C
C     INITIALIZATION
      CLNGTH = 8.0D0*ATAN(1.0D0)
      IBEG = 0
C
      PRINT *, 'DO YOU WANT THE OUTPUT DIRECTED TO THE TERMINAL? (Y/N)'
      READ(*,'(A1)') TYN
      IF((TYN .EQ. 'Y') .OR. (TYN .EQ. 'y')) THEN
          IOUT = TML_OUT
      ELSE
          IOUT = I_FILE
          OPEN(IOUT,FILE='res2')
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
C
C     READ PROBLEM PARAMETERS FOR A NEW CURVE AND BOUNDARY FUNCTION.
   20 PRINT *,'CURVE DEFINING PARAMETERS; NUMCUR, A, B ?'
      PRINT *,'CHOOSE 1, WITH A,B > 0 FOR AN ELLIPSE'
      PRINT *,'CHOOSE 2, WITH 0 < A < 1, B > 0 FOR A LIMACON'
      PRINT *,'CHOOSE 3, WITH A > 1, B > 0 FOR A CASSINI'
      PRINT *,'CHOOSE 4, WITH ANY A AND B, FOR THE AMOEBA'
      READ (TML_IN,*) NUMCUR,A,B
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
      GO TO 20
C
   30 PRINT *,'GIVE THE ERROR TOLERANCE EPS?'
      READ (TML_IN,*) EPS
C
      PRINT *,'IS THE ERROR TEST TO BE "NORMAL" (GIVE 0) OR'
      PRINT *,'"CONSERVATIVE" (GIVE 1)?'
      READ (TML_IN,*) R_FORM
C
   40 PRINT *,'GIVE IE=0 FOR AN INTERIOR PROBLEM.'
      PRINT *,'GIVE IE=1 FOR AN EXTERIOR PROBLEM.'
      PRINT *,'WHAT IS IE?'
      READ (TML_IN,*) IE
C
      PRINT *,'SPECIFY THE TEST FUNCTION INDEX NUMBF?,'
      PRINT *,'FOR IE=0, CHOOSE NUMBF = 1, 2, OR 3'
      PRINT *,'FOR IE=1, CHOOSE NUMBF = 4 OR 5'
      READ (TML_IN,*) NUMBF
C
C     CHECK THE INPUT PARAMETERS.
      IF(IE .EQ. 0 .AND. 1 .LE. NUMBF .AND. NUMBF .LE. 3) GO TO 50
      IF(IE .EQ. 1 .AND. 4 .LE. NUMBF .AND. NUMBF .LE. 5) GO TO 50
      PRINT *, 'GIVE IE AND NUMBF, AGAIN'
      GO TO 40
C
   50 PRINT *,'PARAMETERS N_THETA AND N_R ARE USED TO DEFINE A SET'
      PRINT *,'OF MESH POINTS INTERIOR TO D, AT WHICH THE POTENTIAL'
      PRINT *,'U IS TO BE EVALUATED.'
      PRINT *,'GIVE N_THETA AND N_R' 
      READ (TML_IN,*) NTHETA,NR
C
      IF(R_FORM .EQ. 0) THEN
          RFC = 'NORMAL'
      ELSE
          RFC = 'CONSERVATIVE'
      END IF
      WRITE (IOUT,FMT=9000) A,B,EPS,NUMBF,NUMCUR,NTHETA,NR,IE,RFC
C
C     SET UP POINTS (X,Y) AT WHICH POTENTIAL SOLUTION U IS
C     TO BE EVALUATED.
C     IE=0 ASSIGNS THE POINTS IN THE INTERIOR OF THE CURVE C.
C     IE=1 ASSIGNS THE POINTS IN THE EXTERIOR OF THE CURVE C.
C
      IF (IE .EQ. 0) THEN
C     IF IT IS AN INTERIOR PROBLEM, CHOOSE (0,0) AS THE FIRST POINT.
          DPTS(1,1) = 0.0D0
          DPTS(2,1) = 0.0D0
      ELSE
C     IF IT IS AN EXTERIOR PROBLEM, CHOOSE A POINT SUFFICIENTLY AWAY
C     FROM THE BOUNDARY FOR THE FIRST POINT.
          DPTS(1,1) = 100.D0
          DPTS(2,1) = 100.D0
      END IF
      HTHETA = CLNGTH/NTHETA
      NP = 1
      FNR = NR+1
      ILOW = 1
      IUP = NTHETA
      DO II = ILOW,IUP
           I = II - 1
           THETA = I*HTHETA
           CALL CURVE(THETA,X,Y,DX,DY,D2X,D2Y)
           DO  J = 1,NR
                NP = NP + 1
                R = FLOAT(J)/FNR
                R = R* (2.0-R)
                IF (IE .EQ. 1) R = 1/R
                    DPTS(1,NP) = R*X
                    DPTS(2,NP) = R*Y
           END DO
      END DO
C
      NPTS = NP
C     CHECK WHETHER A SUFFICIENT WORK STORAGE IS GIVEN FOR
C     THE ARRAY DPTS.
      IF (NPTS .GT. MAXDPTS) THEN
          PRINT *,'MORE WORK SPACE IS NEEDED FOR THE ARRAY DPTS.'
          PRINT *,'INCREASE MAXDPTS TO AT LEAST NPTS=',NPTS
          STOP
      END IF
C
      CALL DRCHLT(IOUT,IDBG,IBEG,IE,BDYFCN,CURVE,CLNGTH,DPTS,NPTS,EPS,
     +            R_FORM,WORK,IWORK,NWORK,MWORK,U,ERROR,IER)
C
C     PRINT RESULTS.
C
      WRITE (IOUT,FMT=9010) IER
      IF (IER. LT. 0) STOP
      WRITE (IOUT,FMT=9020)
      WRITE (IOUT,FMT=9030)
      DO I = 1,NPTS
          X = DPTS(1,I)
          Y = DPTS(2,I)
          ANS = BDYFCN(X,Y)
          ERR = ANS - U(I)
          WRITE (IOUT,FMT=9050) X,Y, U(I),ERR,ERROR(I)
      END DO
      PRINT *, ' '
      PRINT *, 'THE CALCULATION IS COMPLETE.  DO YOU WISH TO'
      PRINT *, 'CONTINUE WITH ANOTHER CALCULATION? (Y/N)'
      READ(*,'(A1)') TYN
      IF((TYN .EQ. 'Y') .OR. (TYN .EQ. 'y')) THEN
          PRINT *, 'DO WANT TO USE THE SAME CURVE AND BOUNDARY'
          PRINT *, 'FUNCTION AGAIN, BUT WITH DIFFERENT VALUES OF'
          PRINT *, 'NR AND NTHETA? (Y/N)'
          READ(*,'(A1)') TYN
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
 9000 FORMAT (/,/,' BOUNDARY PARAMETERS: A=',F7.4,3X,'B=',F7.4,/,
     +       ' EPS=',1PD11.3,5X,'BOUNDARY FUNCTION ',I1,5X,'CURVE ',
     +       I1,/,' NTHETA=',I3,3X,'NR=',I3,3X,'IE=',I3,/,
     +       ' THE ERROR TEST IS ',A12,/)
 9010 FORMAT (//,' IER=',I2)
 9020 FORMAT (/,'In the following table, ERROR is the true error;',/,
     +    'The magnitude of PREDERR is the predicted error bound;',/,
     +    'and if PREDERR is negative, the desired error tolerance',/,
     +    'is predicted to have not been attained.',/)
 9030 FORMAT (/,7X,'X',11X,'Y',10X,'U(X,Y)',11X,'ERROR',6X,'PREDERR',/)
 9050 FORMAT (F10.4,F12.4,1PD20.10,D12.2,D12.2)
      END

      DOUBLE PRECISION FUNCTION BDYFCN(X,Y)
C     --------------------------------
C
C     This function defines the test Dirichlet boundary data.  It
C     also defines the true solutions to the associated Dirichlet
C     boundary value problem.
C
      DOUBLE PRECISION X,Y,RSQ
      INTEGER NUMBF
      INTRINSIC COS,EXP
      COMMON /BLKBF/NUMBF
C
C     DEFINE THE BOUNDARY FUNCTIONS.
C
      GO TO (10,20,30,40,50) NUMBF
C
   10 BDYFCN = 1.0D0
      RETURN
   20 BDYFCN = X
      RETURN
   30 BDYFCN = EXP(X)*COS(Y)
      RETURN
   40 RSQ = X*X+Y*Y
      BDYFCN = X/RSQ
      RETURN
   50 RSQ = X*X+Y*Y
      BDYFCN = EXP(X/RSQ)*COS(Y/RSQ)
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
      DOUBLE PRECISION CS,CS2,D2R,DR,R,SN,SN2,EC,ES
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
     +     +ES*(4.0D0*CS2*SN2 + CS*SN2*SN2)
      D2R = EC*(-8*CS2*CS2 - CS*CS2*CS2 + CS2*CS2*SN*SN
     +        + 8*CS2*SN*SN2 + 8*SN2*SN2)
     +     +ES*(8*CS2*CS2 +8*CS*CS2*SN2 - 8*SN2*SN2
     +        + CS*CS*SN2*SN2 - SN*SN2*SN2)
      DY = DR*SN+R*CS
      D2Y = SN* (D2R-R)+2.0D0*DR*CS
      DX = DR*CS - R*SN
      D2X = (D2R-R)*CS - 2.0D0*DR*SN
      RETURN
      END


C*******************************************************************
C     THIS IS THE END OF USER SUPPLIED INFORMATION
C*******************************************************************

