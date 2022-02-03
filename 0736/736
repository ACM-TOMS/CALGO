C      ALGORITHM 736, COLLECTED ALGORITHMS FROM ACM.
C      THIS WORK PUBLISHED IN TRANSACTIONS ON MATHEMATICAL SOFTWARE,
C      VOL. 20, NO. 4, DECEMBER, 1994, P. 427-435.
C
C*** ellpti.f
      SUBROUTINE ELLPTI(NDIM, MAXDIM, GAMMA, ERRTOL, RESULT, NX, WORK,
     +                  IER)
C
C     PROGRAM FINDS INTEGRAL OVER THE UNIT SPHERE IN R**N
C     OF SQRT(GAMMA(1)*X(1)**2+...+GAMMA(N)*X(N)**2)
C     THE EXPECTED VALUE OF SQRT(X''*A*X)
C     WHERE A HAS EIGENVALUES GAMMA(1),...,GAMMA(N) AND
C     X IS UNIFORMLY DISTRIBUTED ON THE SPHERE
C
C     EQUIVALENTLY, PROGRAM FINDS THE EXPECTED RADIUS
C     OF THE ELLIPSOID X'*B**(-1)*X = 1
C     WHERE B IS DIAG(DELTA(1)**2,...,DELTA(N)**2)
C     WITH AXES DELTA(1),...,DELTA(N)'
C     WHERE GAMMA(I)=1/DELTA(I)**2, 1<=I<=N
C
C     VERSION - 10/15/91
C
C     FORMAL PARAMETERS
C        NDIM    INTEGER         input:  the number of values in GAMMA.
C        MAXDIM  INTEGER         input:  the dimension of GAMMA in the
C                                        main program.
C
C        GAMMA   REAL array(*)   input:  the values of GAMMA.
C
C        ERRTOL  REAL            input:  on input the user's requested
C                                  &     relative error tolerance, and
C                                output: on output the error tolerance
C                                        used by the program.
C
C        RESULT  REAL array(6)   output: RESULT(1) is the computed
C                                        multivariate elliptic integral,
C                                        RESULT(2) is the lower bound
C                                        estimate for the integral,
C                                        RESULT(3) is the upper bound
C                                        estimate for the integral,
C                                        RESULT(4) is the error estimate
C                                        for the integral,
C                                        RESULT(5) is the computed
C                                        surface measure, and
C                                        RESULT(6) is the error estimate
C                                        for the surface measure.
C
C        NX      INTEGER         output: number of function evaluations
C                                        used.
C
C        WORK    REAL array      input:  work space
C                (2*MAXDIM)
C
C        IER     INTEGER         output: the error flag.  See failure
C                                        indications below.
C
C     FAILURE INDICATIONS
C        If IER = 0, then no error was detected and successful
C           convergence was obtained.
C
C        If IER = 1, then the program did not converge to the
C           required tolerance.  The last value for the integral along
C           with the estimated error are reported, allowing the user
C           to evaluate the utility of the results.
C
C        If IER = 2, then at least one value in GAMMA is
C           nonpositive.  All nonpositive values are set equal to zero.
C           This is a warning that the surface measure is undefined and
C           RESULT(5) and RESULT(6) are set to zero.
C
C        If IER = 3, then both of the above conditions for IER = 1
C           and IER = 2 hold.
C
C        If IER = 4, then there are no positive values in GAMMA and
C           the program terminates.
C
C        If IER = 5, then the dimension of GAMMA is incorrect and
C           the program terminates.
C
C
C     CONSTANTS
C        SMALL is set to be 2**(-18) in REAL mode.  This should be
C           changed to 2**(-40) for DOUBLE PRECISION mode.  These
C           values are about 64 times machine epsilon.

C
C     GLOBAL VARIABLES
      INTEGER IER, MAXDIM, NDIM, NX
C
C     $$$ CHOOSE MODE $$$
C     DOUBLE PRECISION ERRTOL, GAMMA(*), RESULT(6), WORK(2,*)
      REAL             ERRTOL, GAMMA(*), RESULT(6), WORK(2,*)
C     LOCAL VARIABLES
      INTEGER MAXIT
C     CONTROLS NUMBER OF FUNCTION EVALUATIONS
      PARAMETER (MAXIT = 14)
C
C     $$$ CHOOSE MODE $$$
C     DOUBLE PRECISION BETAIN, C10, C20, ERRQF, FACTOR, FACVAL,
      REAL             BETAIN, C10, C20, ERRQF, FACTOR, FACVAL,
     +                 FIVE, FOUR, FX, F4, H, INTGRL, LOWQF, MAXGAM,
     +                 ONE, PROD, Q(MAXIT), Q0, SMALL,
     +                 SUM, THREE, TRAP, TWO, UPQF, V, X,
     +                 XSUM, X3, ZERO
      INTEGER I, IMAX, IP1, J, N
      INTRINSIC ABS, MAX, SQRT
      EXTERNAL BETAIN, BOUNDS, FACTVL
C
C     USE SMALL = 2**(-40) FOR DOUBLE PRECISION MODE
C     USE SMALL = 2**(-18) FOR REAL MODE
C
C     $$$ CHOOSE MODE $$$
C     DATA SMALL / 9.094947017729D-13 /
C     DATA ZERO, ONE, TWO, THREE, FOUR, FIVE, C10, C20 / 0.0D0, 1.0D0,
C    +     2.0D0, 3.0D0, 4.0D0, 5.0D0, 10.0D0, 20.0D0 /
      DATA SMALL / 3.8146973E-06 /
      DATA ZERO, ONE, TWO, THREE, FOUR, FIVE, C10, C20 / 0.0E0, 1.0E0,
     +     2.0E0, 3.0E0, 4.0E0, 5.0E0, 10.0E0, 20.0E0 /
C
C     START OF EXECUTABLE CODE
C
C     INITIALIZE IER = 0
      IER = 0
C
C     CHECK BOUNDS OF NDIM - ABORT IF OUT OF BOUNDS
      IF ((NDIM.LT.2).OR.(NDIM.GT.MAXDIM)) THEN
         IER = 5
         RETURN
      ENDIF
      MAXGAM = ZERO
      IMAX = 0
      DO 100, I=1,NDIM
C
C     CHECK VALUES OF GAMMA TO BE POSITIVE - WARNING ONLY
         IF (GAMMA(I).LE.ZERO) THEN
            IER = 2
            GAMMA(I) = ZERO
         ENDIF
C
C     MAKE GAMMA(1) THE MAXIMUM OF THE VALUES OF GAMMA
         IF (GAMMA(I).GT.MAXGAM) THEN
            MAXGAM = GAMMA(I)
            IMAX = I
         ENDIF
 100  CONTINUE
C     CHECK THAT AT LEAST ONE GAMMA VALUE IS POSITIVE - ABORT IF NOT
      IF (IMAX.EQ.0) THEN
         IER = 4
         RETURN
      ENDIF
      GAMMA(IMAX) = GAMMA(1)

      GAMMA(1) = MAXGAM
      DO 110, I=1,NDIM
         WORK(1,I) = (MAXGAM-GAMMA(I))/MAXGAM
C     WORK(2,I) = 1 - WORK(1,I)
         WORK(2,I) = GAMMA(I)/MAXGAM
 110  CONTINUE
C
C     NOTE: WORK(1,1) = ZERO
C     NOTE: WORK(2,1) = ONE
      FACTOR = C20*SQRT(MAXGAM)*BETAIN(NDIM)/NDIM
      ERRTOL = MAX(ERRTOL,SMALL)
      H = ONE
C
C     EVALUATE FUNCTION AT ONE
      FX = ONE
C     START I = 2 SINCE WORK(2,1) = ONE
      DO 120, I=2,NDIM
         FX = FX+WORK(2,I)
 120  CONTINUE
      FX = FX/SQRT(C10)
C
C     FUNCTION IS ZERO AT ZERO
      TRAP = FX/TWO
      NX = 1
C
C     LOOP FOR ROMBERG INTEGRATION
C     REFERENCE IS
C     DUNKL, C. F. (1962), ROMBERG QUADRATURE TO PRESCRIBED ACCURACY,
C        SHARE FILE NUMBER 7090-1481
      DO 130, N=1,MAXIT
         H = H/TWO
         SUM = ZERO
         NX = NX*2
         DO 140, J=1,NX-1,2
            X = J*H
C
C     IMPLEMENT MODIFICATION OF KAHAN SUBSTITUTION
C     W. M. KAHAN, "HANDHELD CALCULATOR EVALUATES INTEGRALS,"
C        HEWLETT-PACKARD JOURNAL,AUGUST 1980.
C        SET V = 5*X**4 - 4*X**5 TO REMOVE SINGULARITIES AT ZERO AND ONE
            X3 = X**3
            V = (FIVE-FOUR*X)*X*X3
            PROD = ONE/(X*(X*(FOUR*X+THREE)+TWO)+ONE)
            XSUM = ONE
            DO 150, I=2,NDIM
               FX = WORK(2,I)+V*WORK(1,I)
               XSUM = XSUM+(WORK(2,I)/FX)
               PROD = PROD*(V/FX)
 150        CONTINUE
            FX = X3*XSUM*SQRT(PROD)
            SUM = SUM+FX
 140     CONTINUE
         SUM = SUM*H
         TRAP = SUM+TRAP/TWO
         Q(N) = TWO*(TRAP+SUM)/THREE
         Q0 = Q(1)
         IF (N.GT.1) THEN
            F4 = FOUR
            DO 160, I=N-1,1,-1
               F4 = F4*FOUR
               IP1 = I+1
               Q(I) = Q(IP1)+(Q(IP1)-Q(I))/(F4-ONE)
 160        CONTINUE
            ERRQF = ABS(Q(1)-Q0)
            IF (ERRQF.LE.(ERRTOL*Q0)) GOTO 200
         ENDIF
 130  CONTINUE
C
C     PROGRAM DID NOT CONVERGE TO REQUIRED TOLERANCE - WARNING ONLY
      IER = IER+1
C
C     SUCCESSFUL CONVERGENCE
 200  CONTINUE
      INTGRL = FACTOR*Q(1)
      ERRQF = FACTOR*ERRQF
C
C     COMPUTE BOUNDS FOR INTEGRAL
      CALL BOUNDS(NDIM, GAMMA, LOWQF, UPQF)
C
C     PASS RESULTS BACK IN VECTOR RESULT(6)
      RESULT(1) = INTGRL
      RESULT(2) = LOWQF
      RESULT(3) = UPQF
      RESULT(4) = ERRQF
C
C     COMPUTE FACVAL = (PRODUCT OF DELTA'S)*(SURFACE MEASURE OF SPHERE)
C     AND SURFACE MEASURE ONLY IF ALL GAMMA VALUES ARE POSITIVE
      IF (IER.LE.1) THEN
         CALL FACTVL(NDIM, GAMMA, FACVAL)
         RESULT(5) = FACVAL*INTGRL
C
C     COMPUTE ERROR ESTIMATE FOR SURFACE MEASURE
         RESULT(6) = FACVAL*ERRQF
      ELSE
         RESULT(5) = ZERO
         RESULT(6) = ZERO
      ENDIF
C     EXIT
      RETURN
      END
***********************************************************************
      FUNCTION BETAIN(N)
C     COMPUTES BETA(1/2,(N+1)/2)**(-1)
C
C     GLOBAL VARIABLES
C
C     $$$ CHOOSE MODE $$$
C     DOUBLE PRECISION BETAIN
      REAL             BETAIN
      INTEGER N
C
C     LOCAL VARIABLES
      INTEGER I
C
C     $$$ CHOOSE MODE $$$
C     DOUBLE PRECISION HALF, PI, ONE
      REAL             HALF, PI, ONE
      INTRINSIC MOD
C
C     $$$ CHOOSE MODE $$$
C     DATA HALF, ONE, PI / 0.5D0, 1.0D0, 3.1415926535897932D0 /
      DATA HALF, ONE, PI / 0.5E0, 1.0E0, 3.1415926535897932E0 /
C
C     START OF EXECUTABLE CODE
C
C     N ODD
      IF (MOD(N,2).NE.0) THEN
         BETAIN = HALF
         DO 100, I=2,N-1,2
            BETAIN = (BETAIN*(I+1))/I
 100     CONTINUE
C
C     N EVEN
      ELSE
         BETAIN = ONE/PI
         DO 110, I=2,N,2
            BETAIN = (BETAIN*I)/(I-1)
 110     CONTINUE
      ENDIF
      RETURN
      END
***********************************************************************
      SUBROUTINE BOUNDS(NDIM, GAMMA, LOWQF, UPQF)
      INTEGER NDIM, I
C
C     $$$ CHOOSE MODE $$$
C     DOUBLE PRECISION GAMMA(*), LOWQF, UPQF, ZERO
      REAL             GAMMA(*), LOWQF, UPQF, ZERO
      INTRINSIC SQRT
C
C     $$$ CHOOSE MODE $$$
C     DATA ZERO / 0.0D0 /
      DATA ZERO / 0.0E0 /
C
C     START OF EXECUTABLE CODE
C
      LOWQF = ZERO
      UPQF = ZERO
      DO 100, I=1,NDIM
         LOWQF = LOWQF+SQRT(GAMMA(I))
         UPQF = UPQF+GAMMA(I)
 100  CONTINUE
      LOWQF = LOWQF/NDIM
      UPQF = SQRT(UPQF/NDIM)
      RETURN
      END
***********************************************************************
      SUBROUTINE FACTVL(NDIM, GAMMA, FACVAL)
      INTEGER I, NDIM
C
C     $$$ CHOOSE MODE $$$
C     DOUBLE PRECISION GAMMA(*), FACVAL, PI, PROD, TWO
      REAL             GAMMA(*), FACVAL, PI, PROD, TWO
      INTRINSIC MOD, SQRT
C
C     $$$ CHOOSE MODE $$$
C     DATA PI, TWO / 3.1415926535897932D0, 2.0D0 /
      DATA PI, TWO / 3.1415926535897932E0, 2.0E0 /
C
C     START OF EXECUTABLE CODE
C
C     NDIM ODD
      IF (MOD(NDIM,2).NE.0) THEN
         PROD = TWO
         DO 100, I=1,NDIM-1,2
             PROD = PROD*TWO*PI/I
 100     CONTINUE
C     NDIM EVEN
      ELSE
         PROD = TWO*PI
         DO 110, I=1,NDIM/2-1
            PROD = PROD*PI/I
 110     CONTINUE
      ENDIF
      FACVAL = GAMMA(1)
      DO 120, I=2,NDIM
         FACVAL = FACVAL*GAMMA(I)
 120  CONTINUE
      FACVAL = PROD/SQRT(FACVAL)
      RETURN
      END

C*** dpmain.f
      PROGRAM MAIN
C
C     PROGRAMMED IN C BY C. F. DUNKL - 11/27/90
C     CONVERTED TO FORTRAN BY D. E. RAMIREZ
C     COPY OF PROGRAM IS AVAILABLE WITHOUT CHARGE FROM D.E.R AT
C     DER@VIRGINIA.EDU
C     ADDRESS FOR BOTH C.F.D. AND D.E.R. IS
C     MATHEMATICS DEPARTMENT
C     MATH/ASTRO BUILDING
C     UNIVERSITY OF VIRGINIA
C     CHARLOTTESVILLE, VA 22903 USA
C
C     VERSION - 01/09/91
C
      INTEGER MAXDIM, WORKSZ
C     SUGGESTED VALUE FOR MAXDIM IS 20
      PARAMETER (MAXDIM = 20, WORKSZ = 2*MAXDIM)
C
C     $$$ CHOOSE MODE $$$
      DOUBLE PRECISION DELTA(MAXDIM), GAMMA(MAXDIM), RESULT(6),
C     REAL             DELTA(MAXDIM), GAMMA(MAXDIM), RESULT(6),
     +                 INTGRL, LOWQF, UPQF, ERRQF, SRFACE,
     +                 ERRSF, ONE, ERRTOL, WORK(WORKSZ), ZERO
      INTEGER I, IER, NDIM, NX
      CHARACTER ANS*1
      INTRINSIC SQRT
      EXTERNAL ELLPTI
C
C     $$$ CHOOSE MODE $$$
      DATA ONE, ZERO / 1.0D0, 0.0D0 /
C     DATA ONE, ZERO / 1.0E0, 0.0E0 /
C
C     START OF EXECUTABLE CODE
C
 700  FORMAT(1X,A)
 710  FORMAT(1X,A,I5)
 720  FORMAT(1X,A,G25.15,2X,G25.15)
 730  FORMAT(1X,A,E12.4)
 740  FORMAT(A1)
 750  FORMAT(5(G14.7,2X))
C
C     DESCRIPTION OF PROGRAM
C
      PRINT 700, 'PROGRAM FINDS INTEGRAL OVER THE UNIT SPHERE IN R**N'
      PRINT 700, 'OF SQRT(GAMMA(1)*X(1)**2+...+GAMMA(N)*X(N)**2)'
      PRINT 700, 'THE EXPECTED VALUE OF THE SQRT(X''*A*X)'
      PRINT 700, 'WHERE A HAS EIGENVALUES GAMMA(1),...,GAMMA(N)'
      PRINT 700, 'AND X IS UNIFORMLY DISTRIBUTED ON THE SPHERE'
      PRINT 700, ' '
C
      PRINT 700, 'EQUIVALENTLY, PROGRAM FINDS THE EXPECTED RADIUS'
      PRINT 700, 'OF THE ELLIPSOID X''*B**(-1)*X = 1'
      PRINT 700, 'WHERE B IS DIAG(DELTA(1)**2,...,DELTA(N)**2)'
      PRINT 700, 'WITH AXES DELTA(1),...,DELTA(N)'
      PRINT 700, 'WHERE GAMMA(I)=1/DELTA(I)**2, 1<=I<=N'
      PRINT 700, ' '
C
C     ENTRY POINT FOR LOOP
 10   CONTINUE
C
C     ENTER FORM OF DATA INPUT - Q OR E
      PRINT 700, 'ENTER Q FOR QUADRATIC FORM, E FOR ELLIPSOIDAL DATA'
      READ 740, ANS
      IF (ANS.EQ.'q') ANS = 'Q'
      IF (ANS.EQ.'e') ANS = 'E'
C     ERROR IN INPUT - REDO, BACK TO 10
      IF ((ANS.NE.'Q').AND.(ANS.NE.'E')) THEN
         PRINT 700, 'ERROR IN INPUT: RE-ENTER'
         GOTO 10
      ENDIF
C
C     VALUE FOR DIMENSION - BETWEEN 2 AND MAXDIM
C     DIMENSION CHECK IS IN SUBROUTINE
      PRINT 710, 'ENTER 2 <= DIMENSION <= ', MAXDIM
      READ *, NDIM
C
C     ENTER VALUES FOR GAMMA OR DELTA
      IF (ANS.EQ.'Q') THEN
C     VALUES FOR COEFFICIENTS OF QUADRATIC FORM
         PRINT 700, 'ENTER POSITIVE GAMMA VALUE FOR QUADRATIC FORM'
         READ *, (GAMMA(I), I=1,NDIM)
         DO 90, I=1,NDIM
            IF (GAMMA(I).GT.ZERO) THEN
               DELTA(I) = ONE/SQRT(GAMMA(I))
            ELSE
               DELTA(I) = ZERO
            ENDIF
 90      CONTINUE
         PRINT 700, 'DELTA VALUES ARE'
         PRINT 750, (DELTA(I), I=1,NDIM)
      ELSE IF (ANS.EQ.'E') THEN
C     VALUES FOR LENGTHS FOR AXES
         PRINT 700, 'ENTER POSITIVE DELTA VALUES FOR ELLIPSOID'
         READ *, (DELTA(I), I=1,NDIM)
         DO 100, I=1,NDIM
            IF (DELTA(I).GT.ZERO) THEN
               GAMMA(I) = ONE/DELTA(I)**2
            ELSE
               GAMMA(I) = ZERO
            ENDIF
 100     CONTINUE
         PRINT 700, 'GAMMA VALUES ARE'
         PRINT 750, (GAMMA(I), I=1,NDIM)
      ENDIF
C
C     ERROR TOLERANCE - SUGGESTED VALUE IS 1D-10 FOR DOUBLE PRECISION MODE
C     ERROR TOLERANCE - SUGGESTED VALUE IS 1E-06 FOR REAL MODE
C     ERROR TOLERANCE IS ABOUT 64*(MACHINE EPSILON)
C     SUBROUTINE MAY INCREASE TOLERANCE
C     SUBROUTINE WILL RETURN TOLERANCE USED
C
C     $$$ CHOOSE MODE $$$
      PRINT 700, 'ENTER ERROR TOLERANCE (E.G.; 1D-10)'
C     PRINT 700, 'ENTER ERROR TOLERANCE (E.G.; 1E-06)'
      READ *, ERRTOL
C
C     COMPUTE INTEGRAL
      CALL ELLPTI(NDIM, MAXDIM, GAMMA, ERRTOL, RESULT, NX, WORK, IER)
C
C     CHECK ERROR CODES
      PRINT 700, ' '
      PRINT 710, 'ERROR CODE = ', IER
C     ERROR IN DIMENSION - PROGRAM ABORTED
      IF (IER.EQ.5) THEN
         PRINT 700, 'ERROR IN DIMENSION INPUT'
         GOTO 20
      ENDIF
C     NO POSITIVE GAMMA VALUES - PROGRAM ABORTED
      IF (IER.EQ.4) THEN
         PRINT 700, 'NO POSITIVE GAMMA VALUES'
         GOTO 20
      ENDIF
C     SOME NONPOSITIVE VALUES IN GAMMA - THESE VALUES SET EQUAL TO ZERO
C     WARNING ONLY - NOTE: SURFACE MEASURE IS NOT DEFINED
      IF ((IER.EQ.2).OR.(IER.EQ.3)) THEN
         PRINT 700, 'WARNING: SOME VALUES OF GAMMA NOT POSITIVE'
      ENDIF
C     PROGRAM DID NOT MEET ERROR TOLERANCE - WARNING ONLY
      IF ((IER.EQ.1).OR.(IER.EQ.3)) THEN
         PRINT 700, 'WARNING: PROGRAM DID NOT MEET ERROR TOLERANCE'
      ENDIF
C
C     SET VALUES FROM SUBROUTINE
      INTGRL = RESULT(1)
      LOWQF  = RESULT(2)
      UPQF   = RESULT(3)
      ERRQF  = RESULT(4)
      SRFACE = RESULT(5)
      ERRSF  = RESULT(6)
C
C     PRINT FUNCTION EVALUATIONS
      PRINT 710, 'NUMBER OF FUNCTION EVALUATIONS = ', NX
      PRINT 700, ' '
C
C     PRINT ERROR TOLERANCE USED
      PRINT 730, 'ERROR TOLERANCE USED = ', ERRTOL
      PRINT 700, ' '
C
C     PRINT INTEGRAL
      PRINT 720, 'INTEGRAL = ', INTGRL
      PRINT 730, 'ESTIMATED ERROR IN INTEGRAL = ', ERRQF
C
C     PRINT BOUNDS
      PRINT 720, 'BOUNDS  = ', LOWQF, UPQF
      PRINT 700, ' '
C
C     PRINT SURFACE MEASURE
      PRINT 720, 'SURFACE MEASURE = ', SRFACE
      PRINT 730, 'ESTIMATED ERROR IN SURFACE MEASURE = ', ERRSF
C
C     ALLOW FOR LOOPING
 20   CONTINUE
      PRINT 700, 'REPEAT PROGRAM? (N=NO, Y=YES)'
      READ 740, ANS
      IF (ANS.EQ.'y') ANS = 'Y'
      IF (ANS.EQ.'n') ANS = 'N'
C     CHECK FOR VALID INPUT
      IF ((ANS.NE.'N').AND.(ANS.NE.'Y')) GOTO 20
      IF (ANS.EQ.'Y') GOTO 10
C
C     NORMAL EXIT
      PRINT 700, 'NORMAL EXIT'
      END
***********************************************************************
      SUBROUTINE ELLPTI(NDIM, MAXDIM, GAMMA, ERRTOL, RESULT, NX, WORK,
     +                  IER)
C
C     PROGRAM FINDS INTEGRAL OVER THE UNIT SPHERE IN R**N
C     OF SQRT(GAMMA(1)*X(1)**2+...+GAMMA(N)*X(N)**2)
C     THE EXPECTED VALUE OF SQRT(X''*A*X)
C     WHERE A HAS EIGENVALUES GAMMA(1),...,GAMMA(N) AND
C     X IS UNIFORMLY DISTRIBUTED ON THE SPHERE
C
C     EQUIVALENTLY, PROGRAM FINDS THE EXPECTED RADIUS
C     OF THE ELLIPSOID X'*B**(-1)*X = 1
C     WHERE B IS DIAG(DELTA(1)**2,...,DELTA(N)**2)
C     WITH AXES DELTA(1),...,DELTA(N)'
C     WHERE GAMMA(I)=1/DELTA(I)**2, 1<=I<=N
C
C     VERSION - 10/15/91
C
C     FORMAL PARAMETERS
C        NDIM    INTEGER         input:  the number of values in GAMMA.
C        MAXDIM  INTEGER         input:  the dimension of GAMMA in the
C                                        main program.
C
C        GAMMA   REAL array(*)   input:  the values of GAMMA.
C
C        ERRTOL  REAL            input:  on input the user's requested
C                                  &     relative error tolerance, and
C                                output: on output the error tolerance
C                                        used by the program.
C
C        RESULT  REAL array(6)   output: RESULT(1) is the computed
C                                        multivariate elliptic integral,
C                                        RESULT(2) is the lower bound
C                                        estimate for the integral,
C                                        RESULT(3) is the upper bound
C                                        estimate for the integral,
C                                        RESULT(4) is the error estimate
C                                        for the integral,
C                                        RESULT(5) is the computed
C                                        surface measure, and
C                                        RESULT(6) is the error estimate
C                                        for the surface measure.
C
C        NX      INTEGER         output: number of function evaluations
C                                        used.
C
C        WORK    REAL array      input:  work space
C                (2*MAXDIM)
C
C        IER     INTEGER         output: the error flag.  See failure
C                                        indications below.
C
C     FAILURE INDICATIONS
C        If IER = 0, then no error was detected and successful
C           convergence was obtained.
C
C        If IER = 1, then the program did not converge to the
C           required tolerance.  The last value for the integral along
C           with the estimated error are reported, allowing the user
C           to evaluate the utility of the results.
C
C        If IER = 2, then at least one value in GAMMA is
C           nonpositive.  All nonpositive values are set equal to zero.
C           This is a warning that the surface measure is undefined and
C           RESULT(5) and RESULT(6) are set to zero.
C
C        If IER = 3, then both of the above conditions for IER = 1
C           and IER = 2 hold.
C
C        If IER = 4, then there are no positive values in GAMMA and
C           the program terminates.
C
C        If IER = 5, then the dimension of GAMMA is incorrect and
C           the program terminates.
C
C
C     CONSTANTS
C        SMALL is set to be 2**(-18) in REAL mode.  This should be
C           changed to 2**(-40) for DOUBLE PRECISION mode.  These
C           values are about 64 times machine epsilon.

C
C     GLOBAL VARIABLES
      INTEGER IER, MAXDIM, NDIM, NX
C
C     $$$ CHOOSE MODE $$$
      DOUBLE PRECISION ERRTOL, GAMMA(*), RESULT(6), WORK(2,*)
C     REAL             ERRTOL, GAMMA(*), RESULT(6), WORK(2,*)
C     LOCAL VARIABLES
      INTEGER MAXIT
C     CONTROLS NUMBER OF FUNCTION EVALUATIONS
      PARAMETER (MAXIT = 14)
C
C     $$$ CHOOSE MODE $$$
      DOUBLE PRECISION BETAIN, C10, C20, ERRQF, FACTOR, FACVAL,
C     REAL             BETAIN, C10, C20, ERRQF, FACTOR, FACVAL,
     +                 FIVE, FOUR, FX, F4, H, INTGRL, LOWQF, MAXGAM,
     +                 ONE, PROD, Q(MAXIT), Q0, SMALL,
     +                 SUM, THREE, TRAP, TWO, UPQF, V, X,
     +                 XSUM, X3, ZERO
      INTEGER I, IMAX, IP1, J, N
      INTRINSIC ABS, MAX, SQRT
      EXTERNAL BETAIN, BOUNDS, FACTVL
C
C     USE SMALL = 2**(-40) FOR DOUBLE PRECISION MODE
C     USE SMALL = 2**(-18) FOR REAL MODE
C
C     $$$ CHOOSE MODE $$$
      DATA SMALL / 9.094947017729D-13 /
      DATA ZERO, ONE, TWO, THREE, FOUR, FIVE, C10, C20 / 0.0D0, 1.0D0,
     +     2.0D0, 3.0D0, 4.0D0, 5.0D0, 10.0D0, 20.0D0 /
C     DATA SMALL / 3.8146973E-06 /
C     DATA ZERO, ONE, TWO, THREE, FOUR, FIVE, C10, C20 / 0.0E0, 1.0E0,
C    +     2.0E0, 3.0E0, 4.0E0, 5.0E0, 10.0E0, 20.0E0 /
C
C     START OF EXECUTABLE CODE
C
C     INITIALIZE IER = 0
      IER = 0
C
C     CHECK BOUNDS OF NDIM - ABORT IF OUT OF BOUNDS
      IF ((NDIM.LT.2).OR.(NDIM.GT.MAXDIM)) THEN
         IER = 5
         RETURN
      ENDIF
      MAXGAM = ZERO
      IMAX = 0
      DO 100, I=1,NDIM
C
C     CHECK VALUES OF GAMMA TO BE POSITIVE - WARNING ONLY
         IF (GAMMA(I).LE.ZERO) THEN
            IER = 2
            GAMMA(I) = ZERO
         ENDIF
C
C     MAKE GAMMA(1) THE MAXIMUM OF THE VALUES OF GAMMA
         IF (GAMMA(I).GT.MAXGAM) THEN
            MAXGAM = GAMMA(I)
            IMAX = I
         ENDIF
 100  CONTINUE
C     CHECK THAT AT LEAST ONE GAMMA VALUE IS POSITIVE - ABORT IF NOT
      IF (IMAX.EQ.0) THEN
         IER = 4
         RETURN
      ENDIF
      GAMMA(IMAX) = GAMMA(1)

      GAMMA(1) = MAXGAM
      DO 110, I=1,NDIM
         WORK(1,I) = (MAXGAM-GAMMA(I))/MAXGAM
C     WORK(2,I) = 1 - WORK(1,I)
         WORK(2,I) = GAMMA(I)/MAXGAM
 110  CONTINUE
C
C     NOTE: WORK(1,1) = ZERO
C     NOTE: WORK(2,1) = ONE
      FACTOR = C20*SQRT(MAXGAM)*BETAIN(NDIM)/NDIM
      ERRTOL = MAX(ERRTOL,SMALL)
      H = ONE
C
C     EVALUATE FUNCTION AT ONE
      FX = ONE
C     START I = 2 SINCE WORK(2,1) = ONE
      DO 120, I=2,NDIM
         FX = FX+WORK(2,I)
 120  CONTINUE
      FX = FX/SQRT(C10)
C
C     FUNCTION IS ZERO AT ZERO
      TRAP = FX/TWO
      NX = 1
C
C     LOOP FOR ROMBERG INTEGRATION
C     REFERENCE IS
C     DUNKL, C. F. (1962), ROMBERG QUADRATURE TO PRESCRIBED ACCURACY,
C        SHARE FILE NUMBER 7090-1481
      DO 130, N=1,MAXIT
         H = H/TWO
         SUM = ZERO
         NX = NX*2
         DO 140, J=1,NX-1,2
            X = J*H
C
C     IMPLEMENT MODIFICATION OF KAHAN SUBSTITUTION
C     W. M. KAHAN, "HANDHELD CALCULATOR EVALUATES INTEGRALS,"
C        HEWLETT-PACKARD JOURNAL,AUGUST 1980.
C        SET V = 5*X**4 - 4*X**5 TO REMOVE SINGULARITIES AT ZERO AND ONE
            X3 = X**3
            V = (FIVE-FOUR*X)*X*X3
            PROD = ONE/(X*(X*(FOUR*X+THREE)+TWO)+ONE)
            XSUM = ONE
            DO 150, I=2,NDIM
               FX = WORK(2,I)+V*WORK(1,I)
               XSUM = XSUM+(WORK(2,I)/FX)
               PROD = PROD*(V/FX)
 150        CONTINUE
            FX = X3*XSUM*SQRT(PROD)
            SUM = SUM+FX
 140     CONTINUE
         SUM = SUM*H
         TRAP = SUM+TRAP/TWO
         Q(N) = TWO*(TRAP+SUM)/THREE
         Q0 = Q(1)
         IF (N.GT.1) THEN
            F4 = FOUR
            DO 160, I=N-1,1,-1
               F4 = F4*FOUR
               IP1 = I+1
               Q(I) = Q(IP1)+(Q(IP1)-Q(I))/(F4-ONE)
 160        CONTINUE
            ERRQF = ABS(Q(1)-Q0)
            IF (ERRQF.LE.(ERRTOL*Q0)) GOTO 200
         ENDIF
 130  CONTINUE
C
C     PROGRAM DID NOT CONVERGE TO REQUIRED TOLERANCE - WARNING ONLY
      IER = IER+1
C
C     SUCCESSFUL CONVERGENCE
 200  CONTINUE
      INTGRL = FACTOR*Q(1)
      ERRQF = FACTOR*ERRQF
C
C     COMPUTE BOUNDS FOR INTEGRAL
      CALL BOUNDS(NDIM, GAMMA, LOWQF, UPQF)
C
C     PASS RESULTS BACK IN VECTOR RESULT(6)
      RESULT(1) = INTGRL
      RESULT(2) = LOWQF
      RESULT(3) = UPQF
      RESULT(4) = ERRQF
C
C     COMPUTE FACVAL = (PRODUCT OF DELTA'S)*(SURFACE MEASURE OF SPHERE)
C     AND SURFACE MEASURE ONLY IF ALL GAMMA VALUES ARE POSITIVE
      IF (IER.LE.1) THEN
         CALL FACTVL(NDIM, GAMMA, FACVAL)
         RESULT(5) = FACVAL*INTGRL
C
C     COMPUTE ERROR ESTIMATE FOR SURFACE MEASURE
         RESULT(6) = FACVAL*ERRQF
      ELSE
         RESULT(5) = ZERO
         RESULT(6) = ZERO
      ENDIF
C     EXIT
      RETURN
      END
***********************************************************************
      FUNCTION BETAIN(N)
C     COMPUTES BETA(1/2,(N+1)/2)**(-1)
C
C     GLOBAL VARIABLES
C
C     $$$ CHOOSE MODE $$$
      DOUBLE PRECISION BETAIN
C     REAL             BETAIN
      INTEGER N
C
C     LOCAL VARIABLES
      INTEGER I
C
C     $$$ CHOOSE MODE $$$
      DOUBLE PRECISION HALF, PI, ONE
C     REAL             HALF, PI, ONE
      INTRINSIC MOD
C
C     $$$ CHOOSE MODE $$$
      DATA HALF, ONE, PI / 0.5D0, 1.0D0, 3.1415926535897932D0 /
C     DATA HALF, ONE, PI / 0.5E0, 1.0E0, 3.1415926535897932E0 /
C
C     START OF EXECUTABLE CODE
C
C     N ODD
      IF (MOD(N,2).NE.0) THEN
         BETAIN = HALF
         DO 100, I=2,N-1,2
            BETAIN = (BETAIN*(I+1))/I
 100     CONTINUE
C
C     N EVEN
      ELSE
         BETAIN = ONE/PI
         DO 110, I=2,N,2
            BETAIN = (BETAIN*I)/(I-1)
 110     CONTINUE
      ENDIF
      RETURN
      END
***********************************************************************
      SUBROUTINE BOUNDS(NDIM, GAMMA, LOWQF, UPQF)
      INTEGER NDIM, I
C
C     $$$ CHOOSE MODE $$$
      DOUBLE PRECISION GAMMA(*), LOWQF, UPQF, ZERO
C     REAL             GAMMA(*), LOWQF, UPQF, ZERO
      INTRINSIC SQRT
C
C     $$$ CHOOSE MODE $$$
      DATA ZERO / 0.0D0 /
C     DATA ZERO / 0.0E0 /
C
C     START OF EXECUTABLE CODE
C
      LOWQF = ZERO
      UPQF = ZERO
      DO 100, I=1,NDIM
         LOWQF = LOWQF+SQRT(GAMMA(I))
         UPQF = UPQF+GAMMA(I)
 100  CONTINUE
      LOWQF = LOWQF/NDIM
      UPQF = SQRT(UPQF/NDIM)
      RETURN
      END
***********************************************************************
      SUBROUTINE FACTVL(NDIM, GAMMA, FACVAL)
      INTEGER I, NDIM
C
C     $$$ CHOOSE MODE $$$
      DOUBLE PRECISION GAMMA(*), FACVAL, PI, PROD, TWO
C     REAL             GAMMA(*), FACVAL, PI, PROD, TWO
      INTRINSIC MOD, SQRT
C
C     $$$ CHOOSE MODE $$$
      DATA PI, TWO / 3.1415926535897932D0, 2.0D0 /
C     DATA PI, TWO / 3.1415926535897932E0, 2.0E0 /
C
C     START OF EXECUTABLE CODE
C
C     NDIM ODD
      IF (MOD(NDIM,2).NE.0) THEN
         PROD = TWO
         DO 100, I=1,NDIM-1,2
             PROD = PROD*TWO*PI/I
 100     CONTINUE
C     NDIM EVEN
      ELSE
         PROD = TWO*PI
         DO 110, I=1,NDIM/2-1
            PROD = PROD*PI/I
 110     CONTINUE
      ENDIF
      FACVAL = GAMMA(1)
      DO 120, I=2,NDIM
         FACVAL = FACVAL*GAMMA(I)
 120  CONTINUE
      FACVAL = PROD/SQRT(FACVAL)
      RETURN
      END

C*** spmain.f
      PROGRAM MAIN
C
C     PROGRAMMED IN C BY C. F. DUNKL - 11/27/90
C     CONVERTED TO FORTRAN BY D. E. RAMIREZ
C     COPY OF PROGRAM IS AVAILABLE WITHOUT CHARGE FROM D.E.R AT
C     DER@VIRGINIA.EDU
C     ADDRESS FOR BOTH C.F.D. AND D.E.R. IS
C     MATHEMATICS DEPARTMENT
C     MATH/ASTRO BUILDING
C     UNIVERSITY OF VIRGINIA
C     CHARLOTTESVILLE, VA 22903 USA
C
C     VERSION - 01/09/91
C
      INTEGER MAXDIM, WORKSZ
C     SUGGESTED VALUE FOR MAXDIM IS 20
      PARAMETER (MAXDIM = 20, WORKSZ = 2*MAXDIM)
C
C     $$$ CHOOSE MODE $$$
C     DOUBLE PRECISION DELTA(MAXDIM), GAMMA(MAXDIM), RESULT(6),
      REAL             DELTA(MAXDIM), GAMMA(MAXDIM), RESULT(6),
     +                 INTGRL, LOWQF, UPQF, ERRQF, SRFACE,
     +                 ERRSF, ONE, ERRTOL, WORK(WORKSZ), ZERO
      INTEGER I, IER, NDIM, NX
      CHARACTER ANS*1
      INTRINSIC SQRT
      EXTERNAL ELLPTI
C
C     $$$ CHOOSE MODE $$$
C     DATA ONE, ZERO / 1.0D0, 0.0D0 /
      DATA ONE, ZERO / 1.0E0, 0.0E0 /
C
C     START OF EXECUTABLE CODE
C
 700  FORMAT(1X,A)
 710  FORMAT(1X,A,I5)
 720  FORMAT(1X,A,G25.15,2X,G25.15)
 730  FORMAT(1X,A,E12.4)
 740  FORMAT(A1)
 750  FORMAT(5(G14.7,2X))
C
C     DESCRIPTION OF PROGRAM
C
      PRINT 700, 'PROGRAM FINDS INTEGRAL OVER THE UNIT SPHERE IN R**N'
      PRINT 700, 'OF SQRT(GAMMA(1)*X(1)**2+...+GAMMA(N)*X(N)**2)'
      PRINT 700, 'THE EXPECTED VALUE OF THE SQRT(X''*A*X)'
      PRINT 700, 'WHERE A HAS EIGENVALUES GAMMA(1),...,GAMMA(N)'
      PRINT 700, 'AND X IS UNIFORMLY DISTRIBUTED ON THE SPHERE'
      PRINT 700, ' '
C
      PRINT 700, 'EQUIVALENTLY, PROGRAM FINDS THE EXPECTED RADIUS'
      PRINT 700, 'OF THE ELLIPSOID X''*B**(-1)*X = 1'
      PRINT 700, 'WHERE B IS DIAG(DELTA(1)**2,...,DELTA(N)**2)'
      PRINT 700, 'WITH AXES DELTA(1),...,DELTA(N)'
      PRINT 700, 'WHERE GAMMA(I)=1/DELTA(I)**2, 1<=I<=N'
      PRINT 700, ' '
C
C     ENTRY POINT FOR LOOP
 10   CONTINUE
C
C     ENTER FORM OF DATA INPUT - Q OR E
      PRINT 700, 'ENTER Q FOR QUADRATIC FORM, E FOR ELLIPSOIDAL DATA'
      READ 740, ANS
      IF (ANS.EQ.'q') ANS = 'Q'
      IF (ANS.EQ.'e') ANS = 'E'
C     ERROR IN INPUT - REDO, BACK TO 10
      IF ((ANS.NE.'Q').AND.(ANS.NE.'E')) THEN
         PRINT 700, 'ERROR IN INPUT: RE-ENTER'
         GOTO 10
      ENDIF
C
C     VALUE FOR DIMENSION - BETWEEN 2 AND MAXDIM
C     DIMENSION CHECK IS IN SUBROUTINE
      PRINT 710, 'ENTER 2 <= DIMENSION <= ', MAXDIM
      READ *, NDIM
C
C     ENTER VALUES FOR GAMMA OR DELTA
      IF (ANS.EQ.'Q') THEN
C     VALUES FOR COEFFICIENTS OF QUADRATIC FORM
         PRINT 700, 'ENTER POSITIVE GAMMA VALUE FOR QUADRATIC FORM'
         READ *, (GAMMA(I), I=1,NDIM)
         DO 90, I=1,NDIM
            IF (GAMMA(I).GT.ZERO) THEN
               DELTA(I) = ONE/SQRT(GAMMA(I))
            ELSE
               DELTA(I) = ZERO
            ENDIF
 90      CONTINUE
         PRINT 700, 'DELTA VALUES ARE'
         PRINT 750, (DELTA(I), I=1,NDIM)
      ELSE IF (ANS.EQ.'E') THEN
C     VALUES FOR LENGTHS FOR AXES
         PRINT 700, 'ENTER POSITIVE DELTA VALUES FOR ELLIPSOID'
         READ *, (DELTA(I), I=1,NDIM)
         DO 100, I=1,NDIM
            IF (DELTA(I).GT.ZERO) THEN
               GAMMA(I) = ONE/DELTA(I)**2
            ELSE
               GAMMA(I) = ZERO
            ENDIF
 100     CONTINUE
         PRINT 700, 'GAMMA VALUES ARE'
         PRINT 750, (GAMMA(I), I=1,NDIM)
      ENDIF
C
C     ERROR TOLERANCE - SUGGESTED VALUE IS 1D-10 FOR DOUBLE PRECISION MODE
C     ERROR TOLERANCE - SUGGESTED VALUE IS 1E-06 FOR REAL MODE
C     ERROR TOLERANCE IS ABOUT 64*(MACHINE EPSILON)
C     SUBROUTINE MAY INCREASE TOLERANCE
C     SUBROUTINE WILL RETURN TOLERANCE USED
C
C     $$$ CHOOSE MODE $$$
C     PRINT 700, 'ENTER ERROR TOLERANCE (E.G.; 1D-10)'
      PRINT 700, 'ENTER ERROR TOLERANCE (E.G.; 1E-06)'
      READ *, ERRTOL
C
C     COMPUTE INTEGRAL
      CALL ELLPTI(NDIM, MAXDIM, GAMMA, ERRTOL, RESULT, NX, WORK, IER)
C
C     CHECK ERROR CODES
      PRINT 700, ' '
      PRINT 710, 'ERROR CODE = ', IER
C     ERROR IN DIMENSION - PROGRAM ABORTED
      IF (IER.EQ.5) THEN
         PRINT 700, 'ERROR IN DIMENSION INPUT'
         GOTO 20
      ENDIF
C     NO POSITIVE GAMMA VALUES - PROGRAM ABORTED
      IF (IER.EQ.4) THEN
         PRINT 700, 'NO POSITIVE GAMMA VALUES'
         GOTO 20
      ENDIF
C     SOME NONPOSITIVE VALUES IN GAMMA - THESE VALUES SET EQUAL TO ZERO
C     WARNING ONLY - NOTE: SURFACE MEASURE IS NOT DEFINED
      IF ((IER.EQ.2).OR.(IER.EQ.3)) THEN
         PRINT 700, 'WARNING: SOME VALUES OF GAMMA NOT POSITIVE'
      ENDIF
C     PROGRAM DID NOT MEET ERROR TOLERANCE - WARNING ONLY
      IF ((IER.EQ.1).OR.(IER.EQ.3)) THEN
         PRINT 700, 'WARNING: PROGRAM DID NOT MEET ERROR TOLERANCE'
      ENDIF
C
C     SET VALUES FROM SUBROUTINE
      INTGRL = RESULT(1)
      LOWQF  = RESULT(2)
      UPQF   = RESULT(3)
      ERRQF  = RESULT(4)
      SRFACE = RESULT(5)
      ERRSF  = RESULT(6)
C
C     PRINT FUNCTION EVALUATIONS
      PRINT 710, 'NUMBER OF FUNCTION EVALUATIONS = ', NX
      PRINT 700, ' '
C
C     PRINT ERROR TOLERANCE USED
      PRINT 730, 'ERROR TOLERANCE USED = ', ERRTOL
      PRINT 700, ' '
C
C     PRINT INTEGRAL
      PRINT 720, 'INTEGRAL = ', INTGRL
      PRINT 730, 'ESTIMATED ERROR IN INTEGRAL = ', ERRQF
C
C     PRINT BOUNDS
      PRINT 720, 'BOUNDS  = ', LOWQF, UPQF
      PRINT 700, ' '
C
C     PRINT SURFACE MEASURE
      PRINT 720, 'SURFACE MEASURE = ', SRFACE
      PRINT 730, 'ESTIMATED ERROR IN SURFACE MEASURE = ', ERRSF
C
C     ALLOW FOR LOOPING
 20   CONTINUE
      PRINT 700, 'REPEAT PROGRAM? (N=NO, Y=YES)'
      READ 740, ANS
      IF (ANS.EQ.'y') ANS = 'Y'
      IF (ANS.EQ.'n') ANS = 'N'
C     CHECK FOR VALID INPUT
      IF ((ANS.NE.'N').AND.(ANS.NE.'Y')) GOTO 20
      IF (ANS.EQ.'Y') GOTO 10
C
C     NORMAL EXIT
      PRINT 700, 'NORMAL EXIT'
      END
***********************************************************************
      SUBROUTINE ELLPTI(NDIM, MAXDIM, GAMMA, ERRTOL, RESULT, NX, WORK,
     +                  IER)
C
C     PROGRAM FINDS INTEGRAL OVER THE UNIT SPHERE IN R**N
C     OF SQRT(GAMMA(1)*X(1)**2+...+GAMMA(N)*X(N)**2)
C     THE EXPECTED VALUE OF SQRT(X''*A*X)
C     WHERE A HAS EIGENVALUES GAMMA(1),...,GAMMA(N) AND
C     X IS UNIFORMLY DISTRIBUTED ON THE SPHERE
C
C     EQUIVALENTLY, PROGRAM FINDS THE EXPECTED RADIUS
C     OF THE ELLIPSOID X'*B**(-1)*X = 1
C     WHERE B IS DIAG(DELTA(1)**2,...,DELTA(N)**2)
C     WITH AXES DELTA(1),...,DELTA(N)'
C     WHERE GAMMA(I)=1/DELTA(I)**2, 1<=I<=N
C
C     VERSION - 10/15/91
C
C     FORMAL PARAMETERS
C        NDIM    INTEGER         input:  the number of values in GAMMA.
C        MAXDIM  INTEGER         input:  the dimension of GAMMA in the
C                                        main program.
C
C        GAMMA   REAL array(*)   input:  the values of GAMMA.
C
C        ERRTOL  REAL            input:  on input the user's requested
C                                  &     relative error tolerance, and
C                                output: on output the error tolerance
C                                        used by the program.
C
C        RESULT  REAL array(6)   output: RESULT(1) is the computed
C                                        multivariate elliptic integral,
C                                        RESULT(2) is the lower bound
C                                        estimate for the integral,
C                                        RESULT(3) is the upper bound
C                                        estimate for the integral,
C                                        RESULT(4) is the error estimate
C                                        for the integral,
C                                        RESULT(5) is the computed
C                                        surface measure, and
C                                        RESULT(6) is the error estimate
C                                        for the surface measure.
C
C        NX      INTEGER         output: number of function evaluations
C                                        used.
C
C        WORK    REAL array      input:  work space
C                (2*MAXDIM)
C
C        IER     INTEGER         output: the error flag.  See failure
C                                        indications below.
C
C     FAILURE INDICATIONS
C        If IER = 0, then no error was detected and successful
C           convergence was obtained.
C
C        If IER = 1, then the program did not converge to the
C           required tolerance.  The last value for the integral along
C           with the estimated error are reported, allowing the user
C           to evaluate the utility of the results.
C
C        If IER = 2, then at least one value in GAMMA is
C           nonpositive.  All nonpositive values are set equal to zero.
C           This is a warning that the surface measure is undefined and
C           RESULT(5) and RESULT(6) are set to zero.
C
C        If IER = 3, then both of the above conditions for IER = 1
C           and IER = 2 hold.
C
C        If IER = 4, then there are no positive values in GAMMA and
C           the program terminates.
C
C        If IER = 5, then the dimension of GAMMA is incorrect and
C           the program terminates.
C
C
C     CONSTANTS
C        SMALL is set to be 2**(-18) in REAL mode.  This should be
C           changed to 2**(-40) for DOUBLE PRECISION mode.  These
C           values are about 64 times machine epsilon.
C
C     GLOBAL VARIABLES
      INTEGER IER, MAXDIM, NDIM, NX
C
C     $$$ CHOOSE MODE $$$
C     DOUBLE PRECISION ERRTOL, GAMMA(*), RESULT(6), WORK(2,*)
      REAL             ERRTOL, GAMMA(*), RESULT(6), WORK(2,*)
C     LOCAL VARIABLES
      INTEGER MAXIT
C     CONTROLS NUMBER OF FUNCTION EVALUATIONS
      PARAMETER (MAXIT = 14)
C
C     $$$ CHOOSE MODE $$$
C     DOUBLE PRECISION BETAIN, C10, C20, ERRQF, FACTOR, FACVAL,
      REAL             BETAIN, C10, C20, ERRQF, FACTOR, FACVAL,
     +                 FIVE, FOUR, FX, F4, H, INTGRL, LOWQF, MAXGAM,
     +                 ONE, PROD, Q(MAXIT), Q0, SMALL,
     +                 SUM, THREE, TRAP, TWO, UPQF, V, X,
     +                 XSUM, X3, ZERO
      INTEGER I, IMAX, IP1, J, N
      INTRINSIC ABS, MAX, SQRT
      EXTERNAL BETAIN, BOUNDS, FACTVL
C
C     USE SMALL = 2**(-40) FOR DOUBLE PRECISION MODE
C     USE SMALL = 2**(-18) FOR REAL MODE
C
C     $$$ CHOOSE MODE $$$
C     DATA SMALL / 9.094947017729D-13 /
C     DATA ZERO, ONE, TWO, THREE, FOUR, FIVE, C10, C20 / 0.0D0, 1.0D0,
C    +     2.0D0, 3.0D0, 4.0D0, 5.0D0, 10.0D0, 20.0D0 /
      DATA SMALL / 3.8146973E-06 /
      DATA ZERO, ONE, TWO, THREE, FOUR, FIVE, C10, C20 / 0.0E0, 1.0E0,
     +     2.0E0, 3.0E0, 4.0E0, 5.0E0, 10.0E0, 20.0E0 /
C
C     START OF EXECUTABLE CODE
C
C     INITIALIZE IER = 0
      IER = 0
C
C     CHECK BOUNDS OF NDIM - ABORT IF OUT OF BOUNDS
      IF ((NDIM.LT.2).OR.(NDIM.GT.MAXDIM)) THEN
         IER = 5
         RETURN
      ENDIF
      MAXGAM = ZERO
      IMAX = 0
      DO 100, I=1,NDIM
C
C     CHECK VALUES OF GAMMA TO BE POSITIVE - WARNING ONLY
         IF (GAMMA(I).LE.ZERO) THEN
            IER = 2
            GAMMA(I) = ZERO
         ENDIF
C
C     MAKE GAMMA(1) THE MAXIMUM OF THE VALUES OF GAMMA
         IF (GAMMA(I).GT.MAXGAM) THEN
            MAXGAM = GAMMA(I)
            IMAX = I
         ENDIF
 100  CONTINUE
C     CHECK THAT AT LEAST ONE GAMMA VALUE IS POSITIVE - ABORT IF NOT
      IF (IMAX.EQ.0) THEN
         IER = 4
         RETURN
      ENDIF
      GAMMA(IMAX) = GAMMA(1)

      GAMMA(1) = MAXGAM
      DO 110, I=1,NDIM
         WORK(1,I) = (MAXGAM-GAMMA(I))/MAXGAM
C     WORK(2,I) = 1 - WORK(1,I)
         WORK(2,I) = GAMMA(I)/MAXGAM
 110  CONTINUE
C
C     NOTE: WORK(1,1) = ZERO
C     NOTE: WORK(2,1) = ONE
      FACTOR = C20*SQRT(MAXGAM)*BETAIN(NDIM)/NDIM
      ERRTOL = MAX(ERRTOL,SMALL)
      H = ONE
C
C     EVALUATE FUNCTION AT ONE
      FX = ONE
C     START I = 2 SINCE WORK(2,1) = ONE
      DO 120, I=2,NDIM
         FX = FX+WORK(2,I)
 120  CONTINUE
      FX = FX/SQRT(C10)
C
C     FUNCTION IS ZERO AT ZERO
      TRAP = FX/TWO
      NX = 1
C
C     LOOP FOR ROMBERG INTEGRATION
C     REFERENCE IS
C     DUNKL, C. F. (1962), ROMBERG QUADRATURE TO PRESCRIBED ACCURACY,
C        SHARE FILE NUMBER 7090-1481
      DO 130, N=1,MAXIT
         H = H/TWO
         SUM = ZERO
         NX = NX*2
         DO 140, J=1,NX-1,2
            X = J*H
C
C     IMPLEMENT MODIFICATION OF KAHAN SUBSTITUTION
C     W. M. KAHAN, "HANDHELD CALCULATOR EVALUATES INTEGRALS,"
C        HEWLETT-PACKARD JOURNAL,AUGUST 1980.
C        SET V = 5*X**4 - 4*X**5 TO REMOVE SINGULARITIES AT ZERO AND ONE
            X3 = X**3
            V = (FIVE-FOUR*X)*X*X3
            PROD = ONE/(X*(X*(FOUR*X+THREE)+TWO)+ONE)
            XSUM = ONE
            DO 150, I=2,NDIM
               FX = WORK(2,I)+V*WORK(1,I)
               XSUM = XSUM+(WORK(2,I)/FX)
               PROD = PROD*(V/FX)
 150        CONTINUE
            FX = X3*XSUM*SQRT(PROD)
            SUM = SUM+FX
 140     CONTINUE
         SUM = SUM*H
         TRAP = SUM+TRAP/TWO
         Q(N) = TWO*(TRAP+SUM)/THREE
         Q0 = Q(1)
         IF (N.GT.1) THEN
            F4 = FOUR
            DO 160, I=N-1,1,-1
               F4 = F4*FOUR
               IP1 = I+1
               Q(I) = Q(IP1)+(Q(IP1)-Q(I))/(F4-ONE)
 160        CONTINUE
            ERRQF = ABS(Q(1)-Q0)
            IF (ERRQF.LE.(ERRTOL*Q0)) GOTO 200
         ENDIF
 130  CONTINUE
C
C     PROGRAM DID NOT CONVERGE TO REQUIRED TOLERANCE - WARNING ONLY
      IER = IER+1
C
C     SUCCESSFUL CONVERGENCE
 200  CONTINUE
      INTGRL = FACTOR*Q(1)
      ERRQF = FACTOR*ERRQF
C
C     COMPUTE BOUNDS FOR INTEGRAL
      CALL BOUNDS(NDIM, GAMMA, LOWQF, UPQF)
C
C     PASS RESULTS BACK IN VECTOR RESULT(6)
      RESULT(1) = INTGRL
      RESULT(2) = LOWQF
      RESULT(3) = UPQF
      RESULT(4) = ERRQF
C
C     COMPUTE FACVAL = (PRODUCT OF DELTA'S)*(SURFACE MEASURE OF SPHERE)
C     AND SURFACE MEASURE ONLY IF ALL GAMMA VALUES ARE POSITIVE
      IF (IER.LE.1) THEN
         CALL FACTVL(NDIM, GAMMA, FACVAL)
         RESULT(5) = FACVAL*INTGRL
C
C     COMPUTE ERROR ESTIMATE FOR SURFACE MEASURE
         RESULT(6) = FACVAL*ERRQF
      ELSE
         RESULT(5) = ZERO
         RESULT(6) = ZERO
      ENDIF
C     EXIT
      RETURN
      END
***********************************************************************
      FUNCTION BETAIN(N)
C     COMPUTES BETA(1/2,(N+1)/2)**(-1)
C
C     GLOBAL VARIABLES
C
C     $$$ CHOOSE MODE $$$
C     DOUBLE PRECISION BETAIN
      REAL             BETAIN
      INTEGER N
C
C     LOCAL VARIABLES
      INTEGER I
C
C     $$$ CHOOSE MODE $$$
C     DOUBLE PRECISION HALF, PI, ONE
      REAL             HALF, PI, ONE
      INTRINSIC MOD
C
C     $$$ CHOOSE MODE $$$
C     DATA HALF, ONE, PI / 0.5D0, 1.0D0, 3.1415926535897932D0 /
      DATA HALF, ONE, PI / 0.5E0, 1.0E0, 3.1415926535897932E0 /
C
C     START OF EXECUTABLE CODE
C
C     N ODD
      IF (MOD(N,2).NE.0) THEN
         BETAIN = HALF
         DO 100, I=2,N-1,2
            BETAIN = (BETAIN*(I+1))/I
 100     CONTINUE
C
C     N EVEN
      ELSE
         BETAIN = ONE/PI
         DO 110, I=2,N,2
            BETAIN = (BETAIN*I)/(I-1)
 110     CONTINUE
      ENDIF
      RETURN
      END
***********************************************************************
      SUBROUTINE BOUNDS(NDIM, GAMMA, LOWQF, UPQF)
      INTEGER NDIM, I
C
C     $$$ CHOOSE MODE $$$
C     DOUBLE PRECISION GAMMA(*), LOWQF, UPQF, ZERO
      REAL             GAMMA(*), LOWQF, UPQF, ZERO
      INTRINSIC SQRT
C
C     $$$ CHOOSE MODE $$$
C     DATA ZERO / 0.0D0 /
      DATA ZERO / 0.0E0 /
C
C     START OF EXECUTABLE CODE
C
      LOWQF = ZERO
      UPQF = ZERO
      DO 100, I=1,NDIM
         LOWQF = LOWQF+SQRT(GAMMA(I))
         UPQF = UPQF+GAMMA(I)
 100  CONTINUE
      LOWQF = LOWQF/NDIM
      UPQF = SQRT(UPQF/NDIM)
      RETURN
      END
***********************************************************************
      SUBROUTINE FACTVL(NDIM, GAMMA, FACVAL)
      INTEGER I, NDIM
C
C     $$$ CHOOSE MODE $$$
C     DOUBLE PRECISION GAMMA(*), FACVAL, PI, PROD, TWO
      REAL             GAMMA(*), FACVAL, PI, PROD, TWO
      INTRINSIC MOD, SQRT
C
C     $$$ CHOOSE MODE $$$
C     DATA PI, TWO / 3.1415926535897932D0, 2.0D0 /
      DATA PI, TWO / 3.1415926535897932E0, 2.0E0 /
C
C     START OF EXECUTABLE CODE
C
C     NDIM ODD
      IF (MOD(NDIM,2).NE.0) THEN
         PROD = TWO
         DO 100, I=1,NDIM-1,2
             PROD = PROD*TWO*PI/I
 100     CONTINUE
C     NDIM EVEN
      ELSE
         PROD = TWO*PI
         DO 110, I=1,NDIM/2-1
            PROD = PROD*PI/I
 110     CONTINUE
      ENDIF
      FACVAL = GAMMA(1)
      DO 120, I=2,NDIM
         FACVAL = FACVAL*GAMMA(I)
 120  CONTINUE
      FACVAL = PROD/SQRT(FACVAL)
      RETURN
      END
