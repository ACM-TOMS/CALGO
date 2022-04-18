C      ALGORITHM 566, COLLECTED ALGORITHMS FROM ACM.
C      THIS WORK PUBLISHED IN TRANSACTIONS ON MATHEMATICAL SOFTWARE,
C      VOL. 20, NO. 3, SEPTEMBER, 1994, PP. 282-285.
C ----------------------------------------------------------
C THIS FILE CONTAINS THE PROGRAMS ASSOCIATED WITH A
C REMARK SUBMITTED BY V. AVERBUKH, S. FIGUEROA & T. SCHLICK
C TO ALGORITHM 566 (J. MORE', B. GARBOW & K. HILLSTROM,
C ACM TOMS, VOL. 7, PAGES 14-41 AND 136-140, 1981).
C OUR SUPPLEMENTARY PROGRAM, HESFCN, COMPUTES THE SECOND
C DERIVATIVES OF THE 18 TEST FUNCTIONS IN ALGORITHM 566
C FOR UNCONSTRAINED NONLINEAR OPTIMIZATION.
C INCLUDED IN THIS FILE ARE THE FORTRAN PROGRAM
C SEGMENTS OF HESFCN (DOUBLE AND SINGLE PRECISION),
C TESTING PROGRAMS, AND INPUT DATA FILES (SEE BELOW).
C
C A FULL DESCRIPTION OF THE DERIVATIVE FORMULAS PROGRAMMED
C IN HESFCN IS AVAILABLE IN TECHNICAL REPORT 610, COURANT
C INSTITUTE OF MATHEMATICAL SCIENCES, COMPUTER SCIENCE
C DEPARTMENT, NEW YORK UNIVERSITY, 1992,  ENTITLED:
C
C     "HESFCN --- A FORTRAN PACKAGE OF HESSIAN
C      SUBROUTINES FOR TESTING NONLINEAR OPTIMIZATION
C      SOFTWARE".
C
C BY   VICTORIA AVERBUKH, SAMUEL FIGUEROA, AND TAMAR SCHLICK,
C
C     COURANT INSTITUE OF MATHEMATICAL SCIENCES
C     251 MERCER STREET
C     NEW YORK UNIVERSITY
C     NEW YORK,  NEW YORK   10012.
C
C     --------------------------------
C
C COMMENTS CAN BE ADDRESSED TO T. SCHLICK AT THE ADDRESS ABOVE 
C OR BY:
C
C     E-MAIL:    SCHLICK@ACFCLU.NYU.EDU,
C     PHONE:     (212) 998 - 3116, OR
C     FAX:       (212) 995 - 4121.
C
C ----------------------------------------------------------
C THERE ARE FIVE PROGRAM SEGMENTS IN THIS FILE
C (FOLLOWING THESE COMMENTS):
C
C 1. THE HESFCN ROUTINE, DOUBLE PRECISION (A SUPPLEMENT
C    TO SECTION 4 OF ALGORITHM 566)
C
C 2. THE HESFCN ROUTINES, SINGLE PRECISION (A SUPPLEMENT
C    TO SECTION 7 OF ALGORITHM 566)
C
C 3. DRIVER AND ROUTINES FOR TESTING THE SECOND DERIVATIVES
C    OF HESFCN USING TAYLOR EXPANSIONS, DOUBLE PRECISION
C
C 4. DRIVER AND ROUTINES FOR TESTING THE SECOND DERIVATIVES
C    OF HESFCN USING TAYLOR EXPANSIONS, SINGLE PRECISION.
C    NOTE: THE SINGLE PRECISION VERSION OF THE TESTING PROGRAM
C    ----  WILL PROBABLY PERFORM SATISFACTORILY ONLY ON COMPUTERS
C    SUCH AS CRAY SUPERCOMPUTERS, IN WHICH SINGLE PRECISION IS
C    CLOSER TO THE WIDTH OF MANY COMPUTERS' DOUBLE PRECISION.
C
C 5. INPUT FILE FOR TESTING HESFCN (COMMENTED).
C    TO TEST HESFCN, USE AN UNCOMMENTED VERSION OF THIS FILE
C    (INPUT UNIT 5) WITH THE TESTING PROGRAM (SEGMENTS 3 OR
C    4) AND ALGORITHM 566.
C
C ----------------------------------------------------------
C SEGMENT 1: HESFCN, DOUBLE PRECISION
C ----------------------------------------------------------
      SUBROUTINE HESFCN (N,X,HESD,HESL,NPROB)
      INTEGER N,NPROB
      DOUBLE PRECISION X(N),HESD(N),HESL(*)
C     **********
C
C     SUBROUTINE HESFCN
C
C     THIS SUBROUTINE DEFINES THE HESSIAN MATRICES OF 18
C     NONLINEAR UNCONSTRAINED MINIMIZATION PROBLEMS.  THE PROBLEM
C     DIMENSIONS ARE AS DESCRIBED IN OBJFCN.
C
C     THE SUBROUTINE STATEMENT IS
C
C       SUBROUTINE HESFCN (N, X, HESD, HESL, NPROB)
C
C     WHERE
C
C       N IS A POSITIVE INTEGER INPUT VARIABLE.
C
C       X IS AN INPUT ARRAY OF LENGTH N.
C
C       HESD IS AN OUTPUT ARRAY OF LENGTH N WHICH CONTAINS THE
C         DIAGONAL COMPONENTS OF THE HESSIAN MATRIX OF THE NPROB
C         OBJECTIVE FUNCTION EVALUATED AT X.
C
C       HESL IS AN OUTPUT ARRAY OF LENGTH N*(N-1)/2 WHICH CONTAINS
C         THE LOWER TRIANGULAR PART OF THE HESSIAN MATRIX OF THE
C         NPROB OBJECTIVE FUNCTION EVALUATED AT X.
C
C       NPROB IS A POSITIVE INTEGER INPUT VARIABLE WHICH DEFINES THE
C         NUMBER OF THE PROBLEM.  NPROB MUST NOT EXCEED 18.
C
C     SUBPROGRAMS CALLED
C
C       FORTRAN-SUPPLIED ... ABS, ATAN, COS, EXP, LOG, SIGN, SIN,
C                            SQRT
C
C       INTEGER INLINE FUNCTION IX GIVES THE LOCATION OF A HESSIAN
C       ELEMENT (I,J), I>J, IN HESL
C
C     VICTORIA Z. AVERBUKH, SAMUEL A. FIGUEROA, AND
C     TAMAR SCHLICK, 1993.
C     **********
      INTEGER I, J, K, M, II, JJ, IX, IVAR
      DOUBLE PRECISION ZERO, ONE, TWO, THREE, FOUR, FIVE, SIX, EIGHT,
     1      NINE, TEN, FIFTY, CP0001, CP1, CP2, CP25, CP5, C1P5, C2P25,
     2      C2P625, C3P5, C12, C19P8, C25, C29, C50, C90, C100, C120,
     3      C180, C200, C200P2, C202, C220P2, C360, C400, C1000, C1080,
     4      C1200, C2000, C20000, C2E8, C4E8, AP, BP, PI
      PARAMETER (ZERO=0.0D0, ONE=1.0D0, TWO=2.0D0, THREE=3.0D0,
     1      FOUR=4.0D0, FIVE=5.0D0, SIX=6.0D0, EIGHT=8.0D0, NINE=9.0D0,
     2      TEN=1.0D1, FIFTY=5.0D1, CP0001=1.0D-4, CP1=1.0D-1,
     3      CP2=2.0D-1, CP25=2.5D-1, CP5=5.0D-1, C1P5=1.5D0,
     4      C2P25=2.25D0, C2P625=2.625D0, C3P5=3.5D0, C12=1.2D1,
     5      C19P8=1.98D1, C25=2.5D1, C29=2.9D1, C50=5.0D1, C90=9.0D1,
     6      C100=1.0D2, C120=1.2D2, C180=1.8D2, C200=2.0D2,
     7      C200P2=2.002D2, C202=2.02D2, C220P2=2.202D2, C360=3.6D2,
     8      C400=4.0D2, C1000=1.0D3, C1080=1.08D3, C1200=1.2D3,
     9      C2000=2.0D3, C20000=2.0D4, C2E8=2.0D8, C4E8=4.0D8,
     1      AP=1.0D-5, BP=ONE, PI=3.141592653589793D0)
      DOUBLE PRECISION ARG, D1, D2, D3, LOGR, P1, P2, PIARG, PIARG2,
     1      R, R3INV, S1, S2, S3, S1S2, S1S3, S2S3, SS1, SS2,
     2      T, T1, T2, T3, TH, TT, TT1, TT2, TTH
      DOUBLE PRECISION FVEC(50), GVEC(50), Y(15)
      LOGICAL IEV
      DOUBLE PRECISION DFLOAT
      IX(II,JJ)=(II-1)*(II-2)/2+JJ

      DFLOAT(IVAR) = IVAR
      DATA Y /9.0D-4, 4.4D-3, 1.75D-2, 5.4D-2, 1.295D-1, 2.42D-1,
     1      3.521D-1, 3.989D-1, 3.521D-1, 2.42D-1, 1.295D-1, 5.4D-2,
     2      1.75D-2, 4.4D-3, 9.0D-4/

C
C     HESSIAN ROUTINE SELECTOR.
C
      GO TO (100, 200, 300, 400, 500, 600, 700, 800, 900, 1000,
     1      1100, 1200, 1300, 1400, 1500, 1600, 1700, 1800), NPROB
C
C     HELICAL VALLEY FUNCTION.
C
  100 CONTINUE
C
      IF (X(1) .EQ. ZERO) THEN
         TH = SIGN(CP25,X(2))
      ELSE
         TH = ATAN(X(2)/X(1)) / (TWO*PI)
         IF (X(1) .LT. ZERO) TH = TH + CP5
      END IF
      ARG = X(1)**2 + X(2)**2
      PIARG = PI * ARG
      PIARG2 = PIARG * ARG
      R3INV = ONE / SQRT(ARG)**3
      T = X(3) - TEN*TH
      S1 = FIVE*T / PIARG
      P1 = C2000*X(1)*X(2)*T / PIARG2
      P2 = (FIVE/PIARG)**2
      HESD(1) = C200 - C200*(R3INV-P2)*X(2)**2 - P1
      HESD(2) = C200 - C200*(R3INV-P2)*X(1)**2 + P1
      HESD(3) = C202
      HESL(1) = C200*X(1)*X(2)*R3INV +
     1      C1000/PIARG2 * ( T*(X(1)**2-X(2)**2) - FIVE*X(1)*X(2)/PI )
      HESL(2) =  C1000*X(2) / PIARG
      HESL(3) = -C1000*X(1) / PIARG
      RETURN
C
C     BIGGS EXP6 FUNCTION.
C
  200 CONTINUE
      DO 210 I = 1, 6
         HESD(I) = ZERO
  210 CONTINUE
      DO 220 I = 1, 15
         HESL(I) = ZERO
  220 CONTINUE
      DO 230 I = 1, 13
         D1 = DFLOAT(I)/TEN
         D2 = EXP(-D1) - FIVE*EXP(-TEN*D1) + THREE*EXP(-FOUR*D1)
         S1 = EXP(-D1*X(1))
         S2 = EXP(-D1*X(2))
         S3 = EXP(-D1*X(5))
         T = X(3)*S1 - X(4)*S2 + X(6)*S3 - D2
         D2 = D1**2
         S1S2 = S1 * S2
         S1S3 = S1 * S3
         S2S3 = S2 * S3
         HESD(1) = HESD(1) + D2*S1*(T+X(3)*S1)
         HESD(2) = HESD(2) - D2*S2*(T-X(4)*S2)
         HESD(3) = HESD(3) + S1**2
         HESD(4) = HESD(4) + S2**2
         HESD(5) = HESD(5) + D2*S3*(T+X(6)*S3)
         HESD(6) = HESD(6) + S3**2
         HESL(1) = HESL(1) - D2*S1S2
         HESL(2) = HESL(2) - D1*S1*(T+X(3)*S1)
         HESL(3) = HESL(3) + D1*S1S2
         HESL(4) = HESL(4) + D1*S1S2
         HESL(5) = HESL(5) + D1*S2*(T-X(4)*S2)
         HESL(6) = HESL(6) - S1S2
         HESL(7) = HESL(7) + D2*S1S3
         HESL(8) = HESL(8) - D2*S2S3
         HESL(9) = HESL(9) - D1*S1S3
         HESL(10) = HESL(10) + D1*S2S3
         HESL(11) = HESL(11) - D1*S1S3
         HESL(12) = HESL(12) + D1*S2S3
         HESL(13) = HESL(13) + S1S3
         HESL(14) = HESL(14) - S2S3
         HESL(15) = HESL(15) - D1*S3*(T+X(6)*S3)
  230 CONTINUE
      HESD(1) = X(3)*HESD(1)
      HESD(2) = X(4)*HESD(2)
      HESD(5) = X(6)*HESD(5)
      HESL(1) = X(3)*X(4)*HESL(1)
      HESL(3) = X(4)*HESL(3)
      HESL(4) = X(3)*HESL(4)
      HESL(7) = X(3)*X(6)*HESL(7)
      HESL(8) = X(4)*X(6)*HESL(8)
      HESL(9) = X(6)*HESL(9)
      HESL(10) = X(6)*HESL(10)
      HESL(11) = X(3)*HESL(11)
      HESL(12) = X(4)*HESL(12)
      DO 240 I = 1, 6
         HESD(I) = TWO*HESD(I)
  240 CONTINUE
      DO 250 I = 1, 15
         HESL(I) = TWO*HESL(I)
  250 CONTINUE
      RETURN
C
C     GAUSSIAN FUNCTION.
C
  300 CONTINUE
      HESD(1) = ZERO
      HESD(2) = ZERO
      HESD(3) = ZERO
      HESL(1) = ZERO
      HESL(2) = ZERO
      HESL(3) = ZERO
      DO 310 I = 1, 15
         D1 = CP5*DFLOAT(I-1)
         D2 = C3P5 - D1 - X(3)
         ARG = CP5*X(2)*D2**2
         R = EXP(-ARG)
         T = X(1)*R - Y(I)
         T1 = TWO*X(1)*R - Y(I)
         HESD(1) = HESD(1) + R**2
         HESD(2) = HESD(2) + R*T1*D2**4
         HESD(3) = HESD(3) + R*(X(2)*T1*D2**2-T)
         HESL(1) = HESL(1) - R*T1*D2**2
         HESL(2) = HESL(2) + D2*R*T1
         HESL(3) = HESL(3) + D2*R*(T-ARG*T1)
  310 CONTINUE
      HESD(1) = TWO*HESD(1)
      HESD(2) = CP5*X(1)*HESD(2)
      HESD(3) = TWO*X(1)*X(2)*HESD(3)
      HESL(2) = TWO*X(2)*HESL(2)
      HESL(3) = TWO*X(1)*HESL(3)
      RETURN
C
C     POWELL BADLY SCALED FUNCTION.
C
  400 CONTINUE
      S1 = EXP(-X(1))
      S2 = EXP(-X(2))
      T2 = S1 + S2 - ONE - CP0001
      HESD(1) = C2E8*X(2)**2 + TWO*S1*(S1+T2)
      HESD(2) = C2E8*X(1)**2 + TWO*S2*(S2+T2)
      HESL(1) = C4E8*X(1)*X(2) + TWO*S1*S2 - C20000
      RETURN
C
C     BOX 3-DIMENSIONAL FUNCTION.
C
  500 CONTINUE
      HESD(1) = ZERO
      HESD(2) = ZERO
      HESD(3) = ZERO
      HESL(1) = ZERO
      HESL(2) = ZERO
      HESL(3) = ZERO
      DO 510 I = 1, 10
         D1 = DFLOAT(I)
         D2 = D1/TEN
         S1 = EXP(-D2*X(1))
         S2 = EXP(-D2*X(2))
         S3 = EXP(-D2) - EXP(-D1)
         T = S1 - S2 - S3*X(3)
         TH = T*D2**2
         HESD(1) = HESD(1) + TH*S1 + (D2*S1)**2
         HESD(2) = HESD(2) - TH*S2 + (D2*S2)**2
         HESD(3) = HESD(3) + S3**2
         HESL(1) = HESL(1) - S1*S2*D2**2
         HESL(2) = HESL(2) + D2*S1*S3
         HESL(3) = HESL(3) - D2*S2*S3
  510 CONTINUE
      HESD(1) = TWO*HESD(1)
      HESD(2) = TWO*HESD(2)
      HESD(3) = TWO*HESD(3)
      HESL(1) = TWO*HESL(1)
      HESL(2) = TWO*HESL(2)
      HESL(3) = TWO*HESL(3)
      RETURN
C
C     VARIABLY DIMENSIONED FUNCTION.
C
  600 CONTINUE
      T1 = ZERO
      DO 610 J = 1, N
         T1 = T1 + DFLOAT(J)*(X(J)-ONE)
  610 CONTINUE
      T = ONE + SIX*T1**2
      M = 0
      DO 630 J = 1, N
         HESD(J) = TWO + TWO*T*DFLOAT(J)**2
         DO 620 K = 1, J-1
            M = M + 1
            HESL(M) = TWO*T*DFLOAT(J*K)
  620    CONTINUE
  630 CONTINUE
      RETURN
C
C     WATSON FUNCTION.
C
  700 CONTINUE
      DO 710 J = 1, N
         HESD(J) = ZERO
  710 CONTINUE
      DO 720 J = 1, N*(N-1)/2
         HESL(J) = ZERO
  720 CONTINUE
      DO 760 I = 1, 29
         D1 = DFLOAT(I)/C29
         D2 = ONE
         S1 = ZERO
         S2 = X(1)
         DO 730 J = 2, N
            S1 = S1 + DFLOAT(J-1)*D2*X(J)
            D2 = D1*D2
            S2 = S2 + D2*X(J)
  730    CONTINUE
         T = TWO * (S1-S2**2-ONE) * D1**2
         S3 = TWO*D1*S2
         D2 = ONE/D1
         M = 0
         DO 750 J = 1, N
            T1 = DFLOAT(J-1) - S3
            HESD(J) = HESD(J) + (T1**2-T)*D2**2
            D3 = ONE/D1
            DO 740 K = 1, J-1
               M = M + 1
               HESL(M) = HESL(M) + (T1*(DFLOAT(K-1)-S3) - T) * D2*D3
               D3 = D1*D3
  740       CONTINUE
            D2 = D1*D2
  750    CONTINUE
  760 CONTINUE
      T3 = X(2) - X(1)**2 - ONE
      HESD(1) = HESD(1) + ONE - TWO*(T3-TWO*X(1)**2)
      HESD(2) = HESD(2) + ONE
      HESL(1) = HESL(1) - TWO*X(1)
      DO 770 J = 1, N
         HESD(J) = TWO * HESD(J)
  770 CONTINUE
      DO 780 J = 1, N*(N-1)/2
         HESL(J) = TWO * HESL(J)
  780 CONTINUE
      RETURN
C
C     PENALTY FUNCTION I.
C
  800 CONTINUE
      T1 = -CP25
      DO 810 J = 1, N
         T1 = T1 + X(J)**2
  810 CONTINUE
      D1 = TWO*AP
      TH = FOUR*BP*T1
      M = 0
      DO 830 J = 1, N
         HESD(J) = D1 + TH + EIGHT*X(J)**2
         DO 820 K = 1, J-1
            M = M + 1
            HESL(M) = EIGHT*X(J)*X(K)
  820    CONTINUE
  830 CONTINUE
      RETURN
C
C     PENALTY FUNCTION II.
C
  900 CONTINUE
      T1 = -ONE
      DO 910 J = 1, N
         T1 = T1 + DFLOAT(N-J+1)*X(J)**2
  910 CONTINUE
      D1 = EXP(CP1)
      D2 = ONE
      TH = FOUR*BP*T1
      M = 0
      DO 930 J = 1, N
         HESD(J) = EIGHT*BP*(DFLOAT(N-J+1)*X(J))**2 + DFLOAT(N-J+1)*TH
         S1 = EXP(X(J)/TEN)
         IF (J .GT. 1) THEN
            S3 = S1 + S2 - D2*(D1 + ONE)
            HESD(J) = HESD(J) + AP*S1*(S3 + S1 - ONE/D1 + TWO*S1)/C50
            HESD(J-1) = HESD(J-1) + AP*S2*(S2+S3)/C50
            DO 920 K = 1, J-1
               M = M + 1
               T1 = EXP(DFLOAT(K)/TEN)
               HESL(M) = EIGHT*DFLOAT(N-J+1)*DFLOAT(N-K+1)*X(J)*X(K)
  920       CONTINUE
            HESL(M) = HESL(M) + AP*S1*S2/C50
         END IF
         S2 = S1
         D2 = D1*D2
  930 CONTINUE
      HESD(1) = HESD(1) + TWO*BP
      RETURN
C
C     BROWN BADLY SCALED FUNCTION.
C
 1000 CONTINUE
      HESD(1) = TWO + TWO*X(2)**2
      HESD(2) = TWO + TWO*X(1)**2
      HESL(1) = FOUR*X(1)*X(2) - FOUR
      RETURN
C
C     BROWN AND DENNIS FUNCTION.
C
 1100 CONTINUE
      DO 1110 I = 1, 4
         HESD(I) = ZERO
 1110 CONTINUE
      DO 1120 I = 1, 6
         HESL(I) = ZERO
 1120 CONTINUE
      DO 1130 I = 1, 20
         D1 = DFLOAT(I)/FIVE
         D2 = SIN(D1)
         T1 = X(1) + D1*X(2) - EXP(D1)
         T2 = X(3) + D2*X(4) - COS(D1)
         T = EIGHT * T1 * T2
         S1 = C12*T1**2 + FOUR*T2**2
         S2 = C12*T2**2 + FOUR*T1**2
         HESD(1) = HESD(1) + S1
         HESD(2) = HESD(2) + S1*D1**2
         HESD(3) = HESD(3) + S2
         HESD(4) = HESD(4) + S2*D2**2
         HESL(1) = HESL(1) + S1*D1
         HESL(2) = HESL(2) + T
         HESL(4) = HESL(4) + T*D2
         HESL(3) = HESL(3) + T*D1
         HESL(5) = HESL(5) + T*D1*D2
         HESL(6) = HESL(6) + S2*D2
 1130 CONTINUE
      RETURN
C
C     GULF RESEARCH AND DEVELOPMENT FUNCTION.
C
 1200 CONTINUE
      DO 1210 I = 1, 3
         HESD(I) = ZERO
         HESL(I) = ZERO
 1210 CONTINUE
      D1 = TWO/THREE
      DO 1220 I = 1, 99
         ARG = DFLOAT(I)/C100
         R = (-FIFTY*LOG(ARG))**D1+C25-X(2)
         T1 = ABS(R)**X(3)/X(1)
         T2 = EXP(-T1)
         T3 = T1 * T2 * (T1*T2+(T1-ONE)*(T2-ARG))
         T = T1 * T2 * (T2-ARG)
         LOGR = LOG(ABS(R))
         HESD(1) = HESD(1) + T3 - T
         HESD(2) = HESD(2) + (T+X(3)*T3)/R**2
         HESD(3) = HESD(3) + T3*LOGR**2
         HESL(1) = HESL(1) + T3/R
         HESL(2) = HESL(2) - T3*LOGR
         HESL(3) = HESL(3) + (T-X(3)*T3*LOGR)/R
 1220 CONTINUE
      HESD(1) = HESD(1) / X(1)**2
      HESD(2) = HESD(2) * X(3)
      HESL(1) = HESL(1) * X(3)/X(1)
      HESL(2) = HESL(2) / X(1)
      DO 1230 I = 1, 3
         HESD(I) = TWO * HESD(I)
         HESL(I) = TWO * HESL(I)
 1230 CONTINUE
      RETURN
C
C     TRIGONOMETRIC FUNCTION.
C
 1300 CONTINUE
      S1 = ZERO
      DO 1310 J = 1, N
         HESD(J) = SIN(X(J))
         S1 = S1 + COS(X(J))
 1310 CONTINUE
      S2 = ZERO
      M = 0
      DO 1330 J = 1, N
         TH = COS(X(J))
         T = DFLOAT(N+J) - HESD(J) - S1 - DFLOAT(J)*TH
         S2 = S2 + T
         DO 1320 K = 1, J-1
            M = M + 1
            HESL(M) = SIN(X(K))*(DFLOAT(N+J+K)*HESD(J)-TH) -
     *            HESD(J)*COS(X(K))
            HESL(M) = TWO*HESL(M)
 1320    CONTINUE
         HESD(J) = DFLOAT(J*(J+2)+N)*HESD(J)**2 +
     *         TH*(TH-DFLOAT(2*J+2)*HESD(J)) + T*(DFLOAT(J)*TH+HESD(J))
 1330 CONTINUE
      DO 1340 J = 1, N
         HESD(J) = TWO*(HESD(J) + COS(X(J))*S2)
 1340 CONTINUE
      RETURN
C
C     EXTENDED ROSENBROCK FUNCTION.
C
 1400 CONTINUE
      DO 1410 J = 1, N*(N-1)/2
         HESL(J) = ZERO
 1410 CONTINUE
      DO 1420 J = 1, N, 2
         HESD(J+1) = C200
         HESD(J) = C1200*X(J)**2 - C400*X(J+1) + TWO
         HESL(IX(J+1,J)) = -C400*X(J)
 1420 CONTINUE
      RETURN
C
C     EXTENDED POWELL FUNCTION.
C
 1500 CONTINUE
      DO 1510 J = 1, N*(N-1)/2
         HESL(J) = ZERO
 1510 CONTINUE
      DO 1520 J = 1, N, 4
         T2 = X(J+1) - TWO*X(J+2)
         T3 = X(J) - X(J+3)
         S1 = C12 * T2**2
         S2 = C120 * T3**2
         HESD(J) = TWO + S2
         HESD(J+1) = C200 + S1
         HESD(J+2) = TEN + FOUR*S1
         HESD(J+3) = TEN + S2
         HESL(IX(J+1,J)) = TWO*TEN
         HESL(IX(J+2,J)) = ZERO
         HESL(IX(J+2,J+1)) = -TWO*S1
         HESL(IX(J+3,J)) = -S2
         HESL(IX(J+3,J+1)) = ZERO
         HESL(IX(J+3,J+2)) = -TEN
 1520 CONTINUE
      RETURN
C
C     BEALE FUNCTION.
C
 1600 CONTINUE
      S1 = ONE - X(2)
      T1 = C1P5 - X(1)*S1
      S2 = ONE - X(2)**2
      T2 = C2P25 - X(1)*S2
      S3 = ONE - X(2)**3
      T3 = C2P625 - X(1)*S3
      HESD(1) = TWO * (S1**2 + S2**2 + S3**2)
      HESD(2) = TWO*X(1) * (X(1) + TWO*T2 + FOUR*X(1)*X(2)**2 +
     1      SIX*X(2)*T3 + NINE*X(1)*X(2)**4)
      HESL(1) = TWO*(T1-X(1)*S1) + FOUR*X(2)*(T2-X(1)*S2) +
     2      SIX*(T3-X(1)*S3)*X(2)**2
      RETURN
C
C     WOOD FUNCTION.
C
 1700 CONTINUE
      HESD(1) = C1200*X(1)**2 - C400*X(2) + TWO
      HESD(2) = C220P2
      HESD(3) = C1080*X(3)**2 - C360*X(4) + TWO
      HESD(4) = C200P2
      HESL(1) = -C400*X(1)
      HESL(2) = ZERO
      HESL(3) = ZERO
      HESL(4) = ZERO
      HESL(5) = C19P8
      HESL(6) = -C360*X(3)
      RETURN
C
C     CHEBYQUAD FUNCTION.
C
 1800 CONTINUE
      DO 1810 I = 1, N
         FVEC(I) = ZERO
 1810 CONTINUE
      DO 1830 J = 1, N
         T1 = ONE
         T2 = TWO*X(J) - ONE
         T = TWO*T2
         DO 1820 I = 1, N
            FVEC(I) = FVEC(I) + T2
            TH = T*T2 - T1
            T1 = T2
            T2 = TH
 1820    CONTINUE
 1830 CONTINUE
      D1 = ONE/FLOAT(N)
      IEV = .FALSE.
      DO 1840 I = 1, N
         FVEC(I) = D1*FVEC(I)
         IF (IEV) FVEC(I) = FVEC(I) + ONE/(DFLOAT(I)**2 - ONE)
         IEV = .NOT. IEV
 1840 CONTINUE
      D2 = TWO*D1
      M = 0
      DO 1880 J = 1, N
         HESD(J) = FOUR*D1
         T1 = ONE
         T2 = TWO*X(J) - ONE
         T = TWO*T2
         S1 = ZERO
         S2 = TWO
         P1 = ZERO
         P2 = ZERO
         GVEC(1) = S2
         DO 1850 I = 2, N
            TH = FOUR*T2 + T*S2 - S1
            S1 = S2
            S2 = TH
            TH = T*T2 - T1
            T1 = T2
            T2 = TH
            TH = EIGHT*S1 + T*P2 - P1
            P1 = P2
            P2 = TH
            GVEC(I) = S2
            HESD(J) = HESD(J) + FVEC(I)*TH + D1*S2**2
 1850    CONTINUE
         HESD(J) = D2*HESD(J)
         DO 1870 K = 1, J-1
            M = M + 1
            HESL(M) = ZERO
            TT1 = ONE
            TT2 = TWO*X(K) - ONE
            TT = TWO*TT2
            SS1 = ZERO
            SS2 = TWO
            DO 1860 I = 1, N
               HESL(M) = HESL(M) + SS2*GVEC(I)
               TTH = FOUR*TT2 + TT*SS2 - SS1
               SS1 = SS2
               SS2 = TTH
               TTH = TT*TT2 - TT1
               TT1 = TT2
               TT2 = TTH
 1860       CONTINUE
            HESL(M) = D2*D1*HESL(M)
 1870    CONTINUE
 1880 CONTINUE
      RETURN
C
C     LAST CARD OF SUBROUTINE HESFCN.
C
      END
C ----------------------------------------------------------
C SEGMENT 2: HESFCN, SINGLE PRECISION
C ----------------------------------------------------------
      SUBROUTINE HESFCN (N,X,HESD,HESL,NPROB)
      INTEGER N,NPROB
      REAL X(N),HESD(N),HESL(*)
C     **********
C
C     SUBROUTINE HESFCN
C
C     THIS SUBROUTINE DEFINES THE HESSIAN MATRICES OF 18
C     NONLINEAR UNCONSTRAINED MINIMIZATION PROBLEMS.  THE PROBLEM
C     DIMENSIONS ARE AS DESCRIBED IN OBJFCN.
C
C     THE SUBROUTINE STATEMENT IS
C
C       SUBROUTINE HESFCN (N, X, HESD, HESL, NPROB)
C
C     WHERE
C
C       N IS A POSITIVE INTEGER INPUT VARIABLE.
C
C       X IS AN INPUT ARRAY OF LENGTH N.
C
C       HESD IS AN OUTPUT ARRAY OF LENGTH N WHICH CONTAINS THE
C         DIAGONAL COMPONENTS OF THE HESSIAN MATRIX OF THE NPROB
C         OBJECTIVE FUNCTION EVALUATED AT X.
C
C       HESL IS AN OUTPUT ARRAY OF LENGTH N*(N-1)/2 WHICH CONTAINS
C         THE LOWER TRIANGULAR PART OF THE HESSIAN MATRIX OF THE
C         NPROB OBJECTIVE FUNCTION EVALUATED AT X.
C
C       NPROB IS A POSITIVE INTEGER INPUT VARIABLE WHICH DEFINES THE
C         NUMBER OF THE PROBLEM.  NPROB MUST NOT EXCEED 18.
C
C     SUBPROGRAMS CALLED
C
C       FORTRAN-SUPPLIED ... ABS, ATAN, COS, EXP, LOG, SIGN, SIN,
C                            SQRT
C
C       INTEGER INLINE FUNCTION IX GIVES THE LOCATION OF A HESSIAN
C       ELEMENT (I,J), I>J, IN HESL
C
C     VICTORIA Z. AVERBUKH, SAMUEL A. FIGUEROA, AND
C     TAMAR SCHLICK, 1993.
C     **********
      INTEGER I, J, K, M, II, JJ, IX, IVAR
      REAL ZERO, ONE, TWO, THREE, FOUR, FIVE, SIX, EIGHT,
     1      NINE, TEN, FIFTY, CP0001, CP1, CP2, CP25, CP5, C1P5, C2P25,
     2      C2P625, C3P5, C12, C19P8, C25, C29, C50, C90, C100, C120,
     3      C180, C200, C200P2, C202, C220P2, C360, C400, C1000, C1080,
     4      C1200, C2000, C20000, C2E8, C4E8, AP, BP, PI
      PARAMETER (ZERO=0.0E0, ONE=1.0E0, TWO=2.0E0, THREE=3.0E0,
     1      FOUR=4.0E0, FIVE=5.0E0, SIX=6.0E0, EIGHT=8.0E0, NINE=9.0E0,
     2      TEN=1.0E1, FIFTY=5.0E1, CP0001=1.0E-4, CP1=1.0E-1,
     3      CP2=2.0E-1, CP25=2.5E-1, CP5=5.0E-1, C1P5=1.5E0,
     4      C2P25=2.25E0, C2P625=2.625E0, C3P5=3.5E0, C12=1.2E1,
     5      C19P8=1.98E1, C25=2.5E1, C29=2.9E1, C50=5.0E1, C90=9.0E1,
     6      C100=1.0E2, C120=1.2E2, C180=1.8E2, C200=2.0E2,
     7      C200P2=2.002E2, C202=2.02E2, C220P2=2.202E2, C360=3.6E2,
     8      C400=4.0E2, C1000=1.0E3, C1080=1.08E3, C1200=1.2E3,
     9      C2000=2.0E3, C20000=2.0E4, C2E8=2.0E8, C4E8=4.0E8,
     1      AP=1.0E-5, BP=ONE, PI=3.141592653589793E0)
      REAL ARG, D1, D2, D3, LOGR, P1, P2, PIARG, PIARG2,
     1      R, R3INV, S1, S2, S3, S1S2, S1S3, S2S3, SS1, SS2,
     2      T, T1, T2, T3, TH, TT, TT1, TT2, TTH
      REAL FVEC(50), GVEC(50), Y(15)
      LOGICAL IEV
      REAL DFLOAT
      IX(II,JJ)=(II-1)*(II-2)/2+JJ

      DFLOAT(IVAR) = IVAR
      DATA Y /9.0E-4, 4.4E-3, 1.75E-2, 5.4E-2, 1.295E-1, 2.42E-1,
     1      3.521E-1, 3.989E-1, 3.521E-1, 2.42E-1, 1.295E-1, 5.4E-2,
     2      1.75E-2, 4.4E-3, 9.0E-4/

C
C     HESSIAN ROUTINE SELECTOR.
C
      GO TO (100, 200, 300, 400, 500, 600, 700, 800, 900, 1000,
     1      1100, 1200, 1300, 1400, 1500, 1600, 1700, 1800), NPROB
C
C     HELICAL VALLEY FUNCTION.
C
  100 CONTINUE
C
      IF (X(1) .EQ. ZERO) THEN
         TH = SIGN(CP25,X(2))
      ELSE
         TH = ATAN(X(2)/X(1)) / (TWO*PI)
         IF (X(1) .LT. ZERO) TH = TH + CP5
      END IF
      ARG = X(1)**2 + X(2)**2
      PIARG = PI * ARG
      PIARG2 = PIARG * ARG
      R3INV = ONE / SQRT(ARG)**3
      T = X(3) - TEN*TH
      S1 = FIVE*T / PIARG
      P1 = C2000*X(1)*X(2)*T / PIARG2
      P2 = (FIVE/PIARG)**2
      HESD(1) = C200 - C200*(R3INV-P2)*X(2)**2 - P1
      HESD(2) = C200 - C200*(R3INV-P2)*X(1)**2 + P1
      HESD(3) = C202
      HESL(1) = C200*X(1)*X(2)*R3INV +
     1      C1000/PIARG2 * ( T*(X(1)**2-X(2)**2) - FIVE*X(1)*X(2)/PI )
      HESL(2) =  C1000*X(2) / PIARG
      HESL(3) = -C1000*X(1) / PIARG
      RETURN
C
C     BIGGS EXP6 FUNCTION.
C
  200 CONTINUE
      DO 210 I = 1, 6
         HESD(I) = ZERO
  210 CONTINUE
      DO 220 I = 1, 15
         HESL(I) = ZERO
  220 CONTINUE
      DO 230 I = 1, 13
         D1 = DFLOAT(I)/TEN
         D2 = EXP(-D1) - FIVE*EXP(-TEN*D1) + THREE*EXP(-FOUR*D1)
         S1 = EXP(-D1*X(1))
         S2 = EXP(-D1*X(2))
         S3 = EXP(-D1*X(5))
         T = X(3)*S1 - X(4)*S2 + X(6)*S3 - D2
         D2 = D1**2
         S1S2 = S1 * S2
         S1S3 = S1 * S3
         S2S3 = S2 * S3
         HESD(1) = HESD(1) + D2*S1*(T+X(3)*S1)
         HESD(2) = HESD(2) - D2*S2*(T-X(4)*S2)
         HESD(3) = HESD(3) + S1**2
         HESD(4) = HESD(4) + S2**2
         HESD(5) = HESD(5) + D2*S3*(T+X(6)*S3)
         HESD(6) = HESD(6) + S3**2
         HESL(1) = HESL(1) - D2*S1S2
         HESL(2) = HESL(2) - D1*S1*(T+X(3)*S1)
         HESL(3) = HESL(3) + D1*S1S2
         HESL(4) = HESL(4) + D1*S1S2
         HESL(5) = HESL(5) + D1*S2*(T-X(4)*S2)
         HESL(6) = HESL(6) - S1S2
         HESL(7) = HESL(7) + D2*S1S3
         HESL(8) = HESL(8) - D2*S2S3
         HESL(9) = HESL(9) - D1*S1S3
         HESL(10) = HESL(10) + D1*S2S3
         HESL(11) = HESL(11) - D1*S1S3
         HESL(12) = HESL(12) + D1*S2S3
         HESL(13) = HESL(13) + S1S3
         HESL(14) = HESL(14) - S2S3
         HESL(15) = HESL(15) - D1*S3*(T+X(6)*S3)
  230 CONTINUE
      HESD(1) = X(3)*HESD(1)
      HESD(2) = X(4)*HESD(2)
      HESD(5) = X(6)*HESD(5)
      HESL(1) = X(3)*X(4)*HESL(1)
      HESL(3) = X(4)*HESL(3)
      HESL(4) = X(3)*HESL(4)
      HESL(7) = X(3)*X(6)*HESL(7)
      HESL(8) = X(4)*X(6)*HESL(8)
      HESL(9) = X(6)*HESL(9)
      HESL(10) = X(6)*HESL(10)
      HESL(11) = X(3)*HESL(11)
      HESL(12) = X(4)*HESL(12)
      DO 240 I = 1, 6
         HESD(I) = TWO*HESD(I)
  240 CONTINUE
      DO 250 I = 1, 15
         HESL(I) = TWO*HESL(I)
  250 CONTINUE
      RETURN
C
C     GAUSSIAN FUNCTION.
C
  300 CONTINUE
      HESD(1) = ZERO
      HESD(2) = ZERO
      HESD(3) = ZERO
      HESL(1) = ZERO
      HESL(2) = ZERO
      HESL(3) = ZERO
      DO 310 I = 1, 15
         D1 = CP5*DFLOAT(I-1)
         D2 = C3P5 - D1 - X(3)
         ARG = CP5*X(2)*D2**2
         R = EXP(-ARG)
         T = X(1)*R - Y(I)
         T1 = TWO*X(1)*R - Y(I)
         HESD(1) = HESD(1) + R**2
         HESD(2) = HESD(2) + R*T1*D2**4
         HESD(3) = HESD(3) + R*(X(2)*T1*D2**2-T)
         HESL(1) = HESL(1) - R*T1*D2**2
         HESL(2) = HESL(2) + D2*R*T1
         HESL(3) = HESL(3) + D2*R*(T-ARG*T1)
  310 CONTINUE
      HESD(1) = TWO*HESD(1)
      HESD(2) = CP5*X(1)*HESD(2)
      HESD(3) = TWO*X(1)*X(2)*HESD(3)
      HESL(2) = TWO*X(2)*HESL(2)
      HESL(3) = TWO*X(1)*HESL(3)
      RETURN
C
C     POWELL BADLY SCALED FUNCTION.
C
  400 CONTINUE
      S1 = EXP(-X(1))
      S2 = EXP(-X(2))
      T2 = S1 + S2 - ONE - CP0001
      HESD(1) = C2E8*X(2)**2 + TWO*S1*(S1+T2)
      HESD(2) = C2E8*X(1)**2 + TWO*S2*(S2+T2)
      HESL(1) = C4E8*X(1)*X(2) + TWO*S1*S2 - C20000
      RETURN
C
C     BOX 3-DIMENSIONAL FUNCTION.
C
  500 CONTINUE
      HESD(1) = ZERO
      HESD(2) = ZERO
      HESD(3) = ZERO
      HESL(1) = ZERO
      HESL(2) = ZERO
      HESL(3) = ZERO
      DO 510 I = 1, 10
         D1 = DFLOAT(I)
         D2 = D1/TEN
         S1 = EXP(-D2*X(1))
         S2 = EXP(-D2*X(2))
         S3 = EXP(-D2) - EXP(-D1)
         T = S1 - S2 - S3*X(3)
         TH = T*D2**2
         HESD(1) = HESD(1) + TH*S1 + (D2*S1)**2
         HESD(2) = HESD(2) - TH*S2 + (D2*S2)**2
         HESD(3) = HESD(3) + S3**2
         HESL(1) = HESL(1) - S1*S2*D2**2
         HESL(2) = HESL(2) + D2*S1*S3
         HESL(3) = HESL(3) - D2*S2*S3
  510 CONTINUE
      HESD(1) = TWO*HESD(1)
      HESD(2) = TWO*HESD(2)
      HESD(3) = TWO*HESD(3)
      HESL(1) = TWO*HESL(1)
      HESL(2) = TWO*HESL(2)
      HESL(3) = TWO*HESL(3)
      RETURN
C
C     VARIABLY DIMENSIONED FUNCTION.
C
  600 CONTINUE
      T1 = ZERO
      DO 610 J = 1, N
         T1 = T1 + DFLOAT(J)*(X(J)-ONE)
  610 CONTINUE
      T = ONE + SIX*T1**2
      M = 0
      DO 630 J = 1, N
         HESD(J) = TWO + TWO*T*DFLOAT(J)**2
         DO 620 K = 1, J-1
            M = M + 1
            HESL(M) = TWO*T*DFLOAT(J*K)
  620    CONTINUE
  630 CONTINUE
      RETURN
C
C     WATSON FUNCTION.
C
  700 CONTINUE
      DO 710 J = 1, N
         HESD(J) = ZERO
  710 CONTINUE
      DO 720 J = 1, N*(N-1)/2
         HESL(J) = ZERO
  720 CONTINUE
      DO 760 I = 1, 29
         D1 = DFLOAT(I)/C29
         D2 = ONE
         S1 = ZERO
         S2 = X(1)
         DO 730 J = 2, N
            S1 = S1 + DFLOAT(J-1)*D2*X(J)
            D2 = D1*D2
            S2 = S2 + D2*X(J)
  730    CONTINUE
         T = TWO * (S1-S2**2-ONE) * D1**2
         S3 = TWO*D1*S2
         D2 = ONE/D1
         M = 0
         DO 750 J = 1, N
            T1 = DFLOAT(J-1) - S3
            HESD(J) = HESD(J) + (T1**2-T)*D2**2
            D3 = ONE/D1
            DO 740 K = 1, J-1
               M = M + 1
               HESL(M) = HESL(M) + (T1*(DFLOAT(K-1)-S3) - T) * D2*D3
               D3 = D1*D3
  740       CONTINUE
            D2 = D1*D2
  750    CONTINUE
  760 CONTINUE
      T3 = X(2) - X(1)**2 - ONE
      HESD(1) = HESD(1) + ONE - TWO*(T3-TWO*X(1)**2)
      HESD(2) = HESD(2) + ONE
      HESL(1) = HESL(1) - TWO*X(1)
      DO 770 J = 1, N
         HESD(J) = TWO * HESD(J)
  770 CONTINUE
      DO 780 J = 1, N*(N-1)/2
         HESL(J) = TWO * HESL(J)
  780 CONTINUE
      RETURN
C
C     PENALTY FUNCTION I.
C
  800 CONTINUE
      T1 = -CP25
      DO 810 J = 1, N
         T1 = T1 + X(J)**2
  810 CONTINUE
      D1 = TWO*AP
      TH = FOUR*BP*T1
      M = 0
      DO 830 J = 1, N
         HESD(J) = D1 + TH + EIGHT*X(J)**2
         DO 820 K = 1, J-1
            M = M + 1
            HESL(M) = EIGHT*X(J)*X(K)
  820    CONTINUE
  830 CONTINUE
      RETURN
C
C     PENALTY FUNCTION II.
C
  900 CONTINUE
      T1 = -ONE
      DO 910 J = 1, N
         T1 = T1 + DFLOAT(N-J+1)*X(J)**2
  910 CONTINUE
      D1 = EXP(CP1)
      D2 = ONE
      TH = FOUR*BP*T1
      M = 0
      DO 930 J = 1, N
         HESD(J) = EIGHT*BP*(DFLOAT(N-J+1)*X(J))**2 + DFLOAT(N-J+1)*TH
         S1 = EXP(X(J)/TEN)
         IF (J .GT. 1) THEN
            S3 = S1 + S2 - D2*(D1 + ONE)
            HESD(J) = HESD(J) + AP*S1*(S3 + S1 - ONE/D1 + TWO*S1)/C50
            HESD(J-1) = HESD(J-1) + AP*S2*(S2+S3)/C50
            DO 920 K = 1, J-1
               M = M + 1
               T1 = EXP(DFLOAT(K)/TEN)
               HESL(M) = EIGHT*DFLOAT(N-J+1)*DFLOAT(N-K+1)*X(J)*X(K)
  920       CONTINUE
            HESL(M) = HESL(M) + AP*S1*S2/C50
         END IF
         S2 = S1
         D2 = D1*D2
  930 CONTINUE
      HESD(1) = HESD(1) + TWO*BP
      RETURN
C
C     BROWN BADLY SCALED FUNCTION.
C
 1000 CONTINUE
      HESD(1) = TWO + TWO*X(2)**2
      HESD(2) = TWO + TWO*X(1)**2
      HESL(1) = FOUR*X(1)*X(2) - FOUR
      RETURN
C
C     BROWN AND DENNIS FUNCTION.
C
 1100 CONTINUE
      DO 1110 I = 1, 4
         HESD(I) = ZERO
 1110 CONTINUE
      DO 1120 I = 1, 6
         HESL(I) = ZERO
 1120 CONTINUE
      DO 1130 I = 1, 20
         D1 = DFLOAT(I)/FIVE
         D2 = SIN(D1)
         T1 = X(1) + D1*X(2) - EXP(D1)
         T2 = X(3) + D2*X(4) - COS(D1)
         T = EIGHT * T1 * T2
         S1 = C12*T1**2 + FOUR*T2**2
         S2 = C12*T2**2 + FOUR*T1**2
         HESD(1) = HESD(1) + S1
         HESD(2) = HESD(2) + S1*D1**2
         HESD(3) = HESD(3) + S2
         HESD(4) = HESD(4) + S2*D2**2
         HESL(1) = HESL(1) + S1*D1
         HESL(2) = HESL(2) + T
         HESL(4) = HESL(4) + T*D2
         HESL(3) = HESL(3) + T*D1
         HESL(5) = HESL(5) + T*D1*D2
         HESL(6) = HESL(6) + S2*D2
 1130 CONTINUE
      RETURN
C
C     GULF RESEARCH AND DEVELOPMENT FUNCTION.
C
 1200 CONTINUE
      DO 1210 I = 1, 3
         HESD(I) = ZERO
         HESL(I) = ZERO
 1210 CONTINUE
      D1 = TWO/THREE
      DO 1220 I = 1, 99
         ARG = DFLOAT(I)/C100
         R = (-FIFTY*LOG(ARG))**D1+C25-X(2)
         T1 = ABS(R)**X(3)/X(1)
         T2 = EXP(-T1)
         T3 = T1 * T2 * (T1*T2+(T1-ONE)*(T2-ARG))
         T = T1 * T2 * (T2-ARG)
         LOGR = LOG(ABS(R))
         HESD(1) = HESD(1) + T3 - T
         HESD(2) = HESD(2) + (T+X(3)*T3)/R**2
         HESD(3) = HESD(3) + T3*LOGR**2
         HESL(1) = HESL(1) + T3/R
         HESL(2) = HESL(2) - T3*LOGR
         HESL(3) = HESL(3) + (T-X(3)*T3*LOGR)/R
 1220 CONTINUE
      HESD(1) = HESD(1) / X(1)**2
      HESD(2) = HESD(2) * X(3)
      HESL(1) = HESL(1) * X(3)/X(1)
      HESL(2) = HESL(2) / X(1)
      DO 1230 I = 1, 3
         HESD(I) = TWO * HESD(I)
         HESL(I) = TWO * HESL(I)
 1230 CONTINUE
      RETURN
C
C     TRIGONOMETRIC FUNCTION.
C
 1300 CONTINUE
      S1 = ZERO
      DO 1310 J = 1, N
         HESD(J) = SIN(X(J))
         S1 = S1 + COS(X(J))
 1310 CONTINUE
      S2 = ZERO
      M = 0
      DO 1330 J = 1, N
         TH = COS(X(J))
         T = DFLOAT(N+J) - HESD(J) - S1 - DFLOAT(J)*TH
         S2 = S2 + T
         DO 1320 K = 1, J-1
            M = M + 1
            HESL(M) = SIN(X(K))*(DFLOAT(N+J+K)*HESD(J)-TH) -
     *            HESD(J)*COS(X(K))
            HESL(M) = TWO*HESL(M)
 1320    CONTINUE
         HESD(J) = DFLOAT(J*(J+2)+N)*HESD(J)**2 +
     *         TH*(TH-DFLOAT(2*J+2)*HESD(J)) + T*(DFLOAT(J)*TH+HESD(J))
 1330 CONTINUE
      DO 1340 J = 1, N
         HESD(J) = TWO*(HESD(J) + COS(X(J))*S2)
 1340 CONTINUE
      RETURN
C
C     EXTENDED ROSENBROCK FUNCTION.
C
 1400 CONTINUE
      DO 1410 J = 1, N*(N-1)/2
         HESL(J) = ZERO
 1410 CONTINUE
      DO 1420 J = 1, N, 2
         HESD(J+1) = C200
         HESD(J) = C1200*X(J)**2 - C400*X(J+1) + TWO
         HESL(IX(J+1,J)) = -C400*X(J)
 1420 CONTINUE
      RETURN
C
C     EXTENDED POWELL FUNCTION.
C
 1500 CONTINUE
      DO 1510 J = 1, N*(N-1)/2
         HESL(J) = ZERO
 1510 CONTINUE
      DO 1520 J = 1, N, 4
         T2 = X(J+1) - TWO*X(J+2)
         T3 = X(J) - X(J+3)
         S1 = C12 * T2**2
         S2 = C120 * T3**2
         HESD(J) = TWO + S2
         HESD(J+1) = C200 + S1
         HESD(J+2) = TEN + FOUR*S1
         HESD(J+3) = TEN + S2
         HESL(IX(J+1,J)) = TWO*TEN
         HESL(IX(J+2,J)) = ZERO
         HESL(IX(J+2,J+1)) = -TWO*S1
         HESL(IX(J+3,J)) = -S2
         HESL(IX(J+3,J+1)) = ZERO
         HESL(IX(J+3,J+2)) = -TEN
 1520 CONTINUE
      RETURN
C
C     BEALE FUNCTION.
C
 1600 CONTINUE
      S1 = ONE - X(2)
      T1 = C1P5 - X(1)*S1
      S2 = ONE - X(2)**2
      T2 = C2P25 - X(1)*S2
      S3 = ONE - X(2)**3
      T3 = C2P625 - X(1)*S3
      HESD(1) = TWO * (S1**2 + S2**2 + S3**2)
      HESD(2) = TWO*X(1) * (X(1) + TWO*T2 + FOUR*X(1)*X(2)**2 +
     1      SIX*X(2)*T3 + NINE*X(1)*X(2)**4)
      HESL(1) = TWO*(T1-X(1)*S1) + FOUR*X(2)*(T2-X(1)*S2) +
     2      SIX*(T3-X(1)*S3)*X(2)**2
      RETURN
C
C     WOOD FUNCTION.
C
 1700 CONTINUE
      HESD(1) = C1200*X(1)**2 - C400*X(2) + TWO
      HESD(2) = C220P2
      HESD(3) = C1080*X(3)**2 - C360*X(4) + TWO
      HESD(4) = C200P2
      HESL(1) = -C400*X(1)
      HESL(2) = ZERO
      HESL(3) = ZERO
      HESL(4) = ZERO
      HESL(5) = C19P8
      HESL(6) = -C360*X(3)
      RETURN
C
C     CHEBYQUAD FUNCTION.
C
 1800 CONTINUE
      DO 1810 I = 1, N
         FVEC(I) = ZERO
 1810 CONTINUE
      DO 1830 J = 1, N
         T1 = ONE
         T2 = TWO*X(J) - ONE
         T = TWO*T2
         DO 1820 I = 1, N
            FVEC(I) = FVEC(I) + T2
            TH = T*T2 - T1
            T1 = T2
            T2 = TH
 1820    CONTINUE
 1830 CONTINUE
      D1 = ONE/FLOAT(N)
      IEV = .FALSE.
      DO 1840 I = 1, N
         FVEC(I) = D1*FVEC(I)
         IF (IEV) FVEC(I) = FVEC(I) + ONE/(DFLOAT(I)**2 - ONE)
         IEV = .NOT. IEV
 1840 CONTINUE
      D2 = TWO*D1
      M = 0
      DO 1880 J = 1, N
         HESD(J) = FOUR*D1
         T1 = ONE
         T2 = TWO*X(J) - ONE
         T = TWO*T2
         S1 = ZERO
         S2 = TWO
         P1 = ZERO
         P2 = ZERO
         GVEC(1) = S2
         DO 1850 I = 2, N
            TH = FOUR*T2 + T*S2 - S1
            S1 = S2
            S2 = TH
            TH = T*T2 - T1
            T1 = T2
            T2 = TH
            TH = EIGHT*S1 + T*P2 - P1
            P1 = P2
            P2 = TH
            GVEC(I) = S2
            HESD(J) = HESD(J) + FVEC(I)*TH + D1*S2**2
 1850    CONTINUE
         HESD(J) = D2*HESD(J)
         DO 1870 K = 1, J-1
            M = M + 1
            HESL(M) = ZERO
            TT1 = ONE
            TT2 = TWO*X(K) - ONE
            TT = TWO*TT2
            SS1 = ZERO
            SS2 = TWO
            DO 1860 I = 1, N
               HESL(M) = HESL(M) + SS2*GVEC(I)
               TTH = FOUR*TT2 + TT*SS2 - SS1
               SS1 = SS2
               SS2 = TTH
               TTH = TT*TT2 - TT1
               TT1 = TT2
               TT2 = TTH
 1860       CONTINUE
            HESL(M) = D2*D1*HESL(M)
 1870    CONTINUE
 1880 CONTINUE
      RETURN
C
C     LAST CARD OF SUBROUTINE HESFCN.
C
      END
C ----------------------------------------------------------
C SEGMENT 3: DRIVER AND ROUTINES FOR TESTING HESFCN (DOUBLE PRECISION)
C ----------------------------------------------------------
C
      PROGRAM TESTH
C
C TESTH IS THE DRIVER PROGRAM FOR EXERCISING THE VARIOUS COMPONENTS
C OF ALGORITHM 566 WITH THE NEW HESSIAN SEGMENT, ROUTINE HESFCN.
C THE TESTING OF HESFCN IS ACCOMPLISHED THROUGH TAYLOR EXPANSIONS
C (SUBROUTINE TESTGH), WHERE THE RESULTING ERROR FROM THE SECOND-
C ORDER EXPANSION INDICATES WHETHER THE GRADIENT ONLY, OR BOTH THE
C GRADIENT AND HESSIAN, ARE CORRECT.
C
      INTEGER MAXFCN, MAXN
      DOUBLE PRECISION ZERO, ONE, FIVE
      PARAMETER (MAXFCN=18, MAXN=100, ZERO=0.D0, ONE=1.0D0, FIVE=5.0D0)
      INTEGER N, NPROB, NTRIES, NREAD, NWRITE, NVARS(MAXFCN)
      DOUBLE PRECISION FACTOR, F0, X0(MAXN), X(MAXN), Y(MAXN),
     1      H0Y(MAXN), G0(MAXN), H0D(MAXN), H0L(MAXN*(MAXN-1)/2), YHY,
     2      RANVEC(MAXN)
      EXTERNAL INITPT, OBJFCN, GRDFCN, HESFCN, MVPROD
      DATA NREAD, NWRITE /5,6/
      DATA NVARS /3,6,3,2,3,0,0,0,0,2,4,3,0,0,0,2,4,0/
      DATA RANVEC/ 0.908D0, 0.769D0, 0.734D0, 0.644D0,-0.589D0, 0.577D0,
     1   -0.786D0,-0.901D0, 0.517D0, 0.767D0, 0.749D0,-0.978D0, 0.874D0,
     2    0.777D0,-0.945D0,-0.812D0, 0.921D0, 0.580D0,-0.606D0,-0.857D0,
     3   -0.565D0,-0.545D0, 0.637D0, 0.501D0, 0.707D0, 0.513D0, 0.855D0,
     4   -0.969D0, 0.620D0,-0.590D0, 0.659D0, 0.943D0, 0.826D0,-0.575D0,
     5   -0.841D0, 0.693D0,-0.694D0, 0.750D0, 0.574D0, 0.794D0,-0.923D0,
     6   -0.795D0, 0.978D0, 0.778D0, 0.574D0, 0.992D0,-0.704D0, 0.571D0,
     7    0.782D0,-0.626D0,-0.744D0, 0.732D0,-0.981D0,-0.563D0,-0.600D0,
     8   -0.660D0,-0.815D0,-0.563D0, 0.826D0, 0.811D0,-0.902D0, 0.624D0,
     9    0.738D0, 0.695D0,-0.602D0, 0.514D0,-0.951D0,-0.713D0,-0.571D0,
     A    0.974D0,-0.705D0, 0.566D0,-0.943D0,-0.546D0, 0.581D0, 0.536D0,
     B   -0.683D0, 0.627D0,-0.568D0, 0.892D0, 0.728D0, 0.675D0,-0.726D0,
     C   -0.904D0, 0.966D0, 0.826D0, 0.608D0, 0.840D0, 0.954D0, 0.625D0,
     D    0.930D0,-0.736D0,-0.753D0,-0.800D0, 0.909D0, 0.878D0, 0.731D0,
     E   -0.976D0, 0.816D0,-0.720D0/

   10 CONTINUE
C
C READ THE NUMBER ASSOCIATED WITH THE FUNCTION TO BE USED (NPROB), THE
C NUMBER OF VARIABLES (N), AND THE NUMBER OF STARTING POINTS (NTRIES).
C (EACH LINE OF INPUT FILE MUST CONSIST OF 3 INTEGERS, ARBITRARY
C FORMAT, WITH NPROB RANGING BETWEEN 1 AND 18; N MUST BE APPROPRIATE
C FOR THE PROBLEM  -- SEE INPUT FILE -- AND NTRIES MUST BE A POSITIVE 
C INTEGER. FOR FUNCTION NUMBER 18 (CHEBYQUAD), SUBROUTINES GRDFCN 
C AND HESFCN, AS CURRENTLY IMPLEMENTED, CANNOT HANDLE N GREATER THAN 50.
C ALSO, PLEASE NOTE THAT FOR ALL FUNCTIONS, N IS RESTRICTED TO BE NO
C MORE THAN 100 IN ORDER TO ALLOW EIGENVALUE CALCULATIONS RELATED TO THE
C SPECTRUM OF THE HESSIAN TO BE PERFORMED WITHOUT MAJOR DEMANDS ON
C STORAGE AND CPU TIME.  THE INPUT IS CHECKED BELOW.
C
         READ (NREAD, *) NPROB, N, NTRIES
         IF (NPROB .EQ. 0) STOP
         IF (NPROB .LT. 1 .OR. NPROB .GT. 18 .OR. N .LT. 1
     1         .OR. NTRIES .LT. 1) THEN
            WRITE (NWRITE, 850)
            STOP
         ENDIF
         IF (NVARS(NPROB) .NE. 0)
     1      N = NVARS(NPROB)
         IF (NPROB .EQ. 7)
     1      N = MAX(2, MIN(N,31))
         IF (NPROB .EQ. 14)
     1      N = MAX(2, N - MOD(N,2))
         IF (NPROB .EQ. 15)
     1      N = MAX(4, N - MOD(N,4))
         IF (NPROB .EQ. 18)
     1      N = MAX(2, MIN(N,50))
         IF (N .GT. MAXN) THEN
            WRITE (NWRITE, 875) N, MAXN
            GO TO 10
         END IF
C
C OBTAIN THE INITIAL POINT X0, AND COMPUTE THE CORRESPONDING
C FUNCTION VALUE, GRADIENT VECTOR, AND HESSIAN MATRIX
C
         FACTOR = ONE
         DO 40 K = 1, NTRIES
            CALL INITPT(N, X0, NPROB, FACTOR)
            CALL OBJFCN(N, X0, F0, NPROB)
            CALL GRDFCN(N, X0, G0, NPROB)
            CALL HESFCN(N, X0, H0D, H0L, NPROB)
            WRITE (NWRITE, 900) NPROB, N, FACTOR
C
C OBTAIN A PERTURBATION VECTOR Y
C
         DO 20 I = 1, N
            IF (X0(I) .NE. ZERO) THEN
               Y(I) = X0(I) * RANVEC(I)
            ELSE
               Y(I) = RANVEC(I)
            ENDIF
   20    CONTINUE
         WRITE (NWRITE,925)
         WRITE (NWRITE,950) (X0(I), I = 1, N)
         WRITE (NWRITE,975)
         WRITE (NWRITE,950) (Y(I), I = 1, N)
C
C COMPUTE THE INNER PRODUCT Y*HY AT X0 AND
C CALL THE DERIVATIVE TESTING FUNCTION
C
            CALL MVPROD(N, H0D, H0L, Y, H0Y)
            YHY = ZERO
            DO 30 I = 1, N
               YHY = YHY + Y(I)*H0Y(I)
   30       CONTINUE
            CALL TESTGH (N, X0, F0, G0, Y, YHY, X, NPROB)
            FACTOR = FIVE*FACTOR
   40    CONTINUE
         GO TO 10

  850 FORMAT (/4X, 'ERROR IN INPUT FILE'/)
  875 FORMAT (/4X, 'N > MAXN:  N =', I6, ', MAXN =', I6,
     +        '.  PLEASE INCREASE PARAMETER MAXN.'/)
  900 FORMAT (/4X,'TESTING FUNCTION ', I2/4X,'WITH', I5,
     *        ' VARIABLES AT THE STANDARD STARTING POINT'/
     *        4X, 'SCALED BY', 1PE16.2/)
  925 FORMAT (/4X, 'X0 VECTOR:')
  950 FORMAT (4(F14.2,2X))
  975 FORMAT (/4X, 'Y VECTOR:')
      END
C***************************************************************
      SUBROUTINE TESTGH(N,XC,FC,GC,Y,YHY,VEC,NPROB)
C
C TESTGH TESTS USER-SUPPLIED GRADIENT (G) AND HESSIAN (H)
C ROUTINES CORRESPONDING TO A GIVEN FUNCTION F.
C TESTING H IS OPTIONAL.
C
C DERIVATIVES ARE TESTED USING A TAYLOR EXPANSION OF F
C AROUND A GIVEN POINT XC. THE TAYLOR SERIES IS EXPANDED
C AT XC + EPS*Y WHERE Y IS A RANDOM  PERTURBATION VECTOR
C AND EPS IS A SCALAR. IF WE DENOTE THE DOT PRODUCT OF 2
C VECTORS A AND B AS (A,B), WE CAN WRITE OUR EXPANSION AS
C
C F(XC+EPS*Y) = F(XC) + EPS * (G,Y) + (1/2)*(EPS**2) * (Y,HY)
C                     + O(EPS**3),
C
C WHERE G AND H ARE BOTH EVALUATED AT XC, AND HY DENOTES A
C HESSIAN/VECTOR PRODUCT. IF ONLY G ROUTINES ARE TESTED, THE
C SECOND-ORDER TAYLOR TERM IS ZERO, AND THE TRUNCATION ERROR
C IS O(EPS**2).
C
C OUR TEST IS PERFORMED BY COMPUTING THIS TAYLOR APPROX. AT
C SMALLER AND SMALLER VALUES OF EPS AND CHECKING TO SEE
C WHETHER CORRECT TRUNCATION ERRORS ARE OBTAINED --
C O(EPS**2) AND  O(EPS**3) IF THE APPROX. IS CORRECT UP TO
C THE G AND H TERMS, RESPECTIVELY.
C
C WE DIVIDE EPS BY 2 AT EVERY STEP AND TEST IF INDEED THE
C TRUNCATION ERRORS DECREASE AS THEY SHOULD.
C (I.E., IF THE ERROR CORRESPONDING TO EPS IS E1,
C THE ERROR FOR EPS/2 SHOULD BE E1/4 IF THE GRADIENT
C IS CORRECT, AND E1/8 IF THE HESSIAN IS ALSO CORRECT).
C OUR VALUE "RATIO" COMPUTES THIS FACTOR OF THE OLD/NEW
C ERRORS.
C
C THE OUTPUT IS A SERIES OF VALUES FOR RATIO PRINTED FOR
C EACH EPS UNTIL THE TRUNCATION ERROR AND/OR EPS IS VERY SMALL.
C IF RATIO TENDS TO 4 OR 8 AS EPS IS DECREASED (AND THE
C ERROR IS RELATIVELY SMALL) G  IS CORRECT OR G&H ARE CORRECT,
C RESPECTIVELY. IF RATIO TENDS TO 2, WHICH IS O(EPS), NEITHER G
C NOR H ARE CORRECT. IF THE RATIO TENDS TO 1, THE ERRORS MAY
C BE TOO LARGE GIVEN THE PERTURBATION VECTOR Y.
C
C THUS IN GENERAL, RELIABLE VALUES OF RATIO SHOULD
C OCCUR WHEN: (1) EPS IS NOT TOO LARGE AND NOT TOO SMALL,
C AND (2) THE DIFFERENCE BETWEEN F(XC+EPS*Y) AND THE
C TAYLOR SERIES APPROXIMATION IS OF REASONABLE MAGNITUDE.
C (THE VALUES OF EPS AND THE ERRORS APPEAR IN THE OUTPUT).
C IN OTHER WORDS, AN ACCURATE VALUE OF RATIO SHOULD
C APPEAR AROUND THE MIDDLE OF OUR SERIES IF Y IS APPROPRIATE.
C DIFFERENT STARTING POINT AND/OR PERTURBATION VECTORS
C CAN BE TRIED.
C
C USAGE: THE USER MUST SUPPLY THE FOLLOWING INPUT
C ------ VARIABLES IN THE FUNCTION CALL:
C
C N      - DIMENSION (NUMBER OF VARIABLES FOR F)
C XC(N)  - OUR CURRENT VECTOR
C FC     - THE FUNCTION VALUE AT XC
C GC(N)  - THE GRADIENT VECTOR AT XC, ON INPUT
C          ON OUTPUT, GC MAY BE CHANGED IF IT IS USED IN
C          THE FUNCTION CALL TO OBTAIN A NEW GRADIENT IN
C          ADDITION TO A NEW FUNCTION VALUE (SEE BELOW).
C Y(N)   - A RANDOM PERTURBATION VECTOR (Y SHOULD BE CHOSEN
C          SO THAT F(XC+Y) IS IN A REASONABLE RANGE FOR THE
C          PROBLEM)
C YHY    - THE MATRIX INNER PRODUCT -- (Y,HY) -- REPRESENTING
C          THE DOT PRODUCT OF Y WITH THE HESSIAN/VECTOR
C          PRODUCT, HY, WHERE H IS EVALUATED AT XC
C          (IF ONLY THE GRADIENT IS TESTED, SET YHY TO ZERO).
C VEC(N) - A WORK VECTOR
C NPROB  - AN INTEGER VARIABLE THAT MAY BE USED IN THE 
C          USER'S FUNCTION CALL
C
C NOTE: THE USER MAY MODIFY THIS ROUTINE TO TEST OTHER
C ----  FUNCTIONS BY REPLACING THE SAMPLE FUNCTION CALL GIVEN 
C ABOVE THE '40 CONTINUE' STATEMENT WITH THE APPROPRIATE
C CALL FOR HIS/HER PROBLEM. THE INSERTED ROUTINE CALL
C SHOULD PRODUCE A NEW FUNCTION VALUE, FVEC, FOR EACH 
C NEW VECTOR VEC=XC+EPS*Y.
C
      INTEGER N, MP
      DOUBLE PRECISION ZERO,ONE,HALF,TWOP23,EPSMCH,EPSLIM,EPS,FC,
     1                 GY,YHY,TAYLOR,DIFF,FVEC,FOLD,RATIO,TEMP,
     2                 XC(N),GC(N),Y(N),VEC(N)
      PARAMETER (ZERO=0.D0, ONE=1.0D0, HALF=0.5D0, TWOP23=8388608.0D0,
     1           MP=6)
      EXTERNAL OBJFCN

      EPSMCH = ONE / (TWOP23 * TWOP23)
C
C  NOTE:  THE LINE ABOVE MAY BE REPLACED WITH:
C  -----  EPSMCH = D1MACH(3)
C  WHERE D1MACH IS A FUNCTION WHICH MAY BE OBTAINED BY SENDING AN
C  ELECTRONIC MAIL MESSAGE TO NETLIB@ORNL.GOV WITH A SUBJECT LINE OR
C  BODY OF "SEND D1MACH FROM CORE".  EPSMCH WOULD THEN HAVE THE VALUE
C  R**(-P), WHERE R IS THE RADIX OF DOUBLE PRECISION NUMBERS AND P IS
C  THE NUMBER OF RADIX-R DIGITS IN THE MANTISSA OR SIGNIFICAND.  IN
C  OTHER WORDS, EPSMCH WOULD HAVE THE VALUE OF THE SMALLEST RELATIVE
C  SPACING.  HERE, WE HAVE SIMPLY SET EPSMCH TO THE VALUE 2**(-46),
C  WHICH SHOULD GENERALLY BE LARGER THAN THE ACTUAL VALUE OF
C  R**(-P).
C
      EPSLIM = EPSMCH * FLOAT(N*N) * 1.D+2
      EPS    = HALF

      WRITE (MP,900)
      GY = ZERO
      DO 10 I = 1, N
         GY = GY + GC(I)*Y(I)
   10 CONTINUE
      WRITE (MP,910) FC,GY,YHY,EPSMCH
      WRITE (MP,940)

      DIFF = ZERO
      FVEC = FC

   20 CONTINUE
      TEMP = DIFF
      FOLD = FVEC

      DO 30 I = 1, N
         VEC(I) = XC(I) + EPS*Y(I)
   30 CONTINUE

      NOUT = 0
      CALL OBJFCN(N, VEC, FVEC, NPROB)
   40 CONTINUE

      TAYLOR = FC + (EPS*GY) + ( (EPS**2) * HALF  * YHY )
      DIFF   = FVEC - TAYLOR

      IF (ABS(DIFF) .LT. ABS(EPSLIM*FVEC)) THEN
         WRITE (MP,920) ABS(EPSLIM*FVEC)
         GOTO 50
      ENDIF

      IF (ABS(FVEC-FOLD) .LT. ABS(EPSLIM*FOLD)) THEN
         WRITE (MP,930) ABS(EPSLIM*FOLD)
         GOTO 50
      ENDIF

      IF (TEMP .EQ. ZERO .OR. DIFF .EQ. ZERO) THEN
          WRITE (MP, 950) EPS,FVEC,TAYLOR,DIFF
      ELSE
          RATIO = TEMP / DIFF
          WRITE (MP, 950) EPS,FVEC,TAYLOR,DIFF,RATIO
      ENDIF

      EPS = EPS * HALF
      IF (EPS .GT. EPSMCH) GOTO 20

   50 RETURN

  900 FORMAT(/T10,'ENTERING TESTGH ROUTINE:'/)
  910 FORMAT(T5,  'THE FUNCTION VALUE AT X               = ',
     + 1PE16.8/T5,'THE FIRST-ORDER TAYLOR TERM,  (G, Y)  = ',
     + 1PE16.8/T5,'THE SECOND-ORDER TAYLOR TERM, (Y,HY)  = ',
     + 1PE16.8//T5,'THE COMPUTED MACHINE PRECISION        = ',
     + 1PE16.8//)
  920 FORMAT(/T5,'DIFF IS SMALL (LESS THAN ', 1PE16.8,
     + ' IN ABSOLUTE VALUE)'/)
  930 FORMAT(/T5,'CHANGE IN FUNCTION VALUE IS VERY SMALL (LESS THAN ',
     + 1PE16.8,' IN ABSOLUTE VALUE)'/)
  940 FORMAT(4X,'EPS',10X,' F   ',10X,' TAYLOR',9X,
     +  ' DIFF.',9X,'RATIO'/)
  950 FORMAT(1X,1PE10.4,1PE16.8,1PE16.8,1PE16.8,1PE16.8)

      END
C***********************************************************************
      SUBROUTINE MVPROD (N, DIAGA, LOWERA, X, Y)
C
C MVPROD PERFORMS THE MATRIX-VECTOR PRODUCT A*X, AND
C STORES THE RESULT IN THE VECTOR Y.  A IS A SYMMETRIC
C NXN MATRIX, WITH DIAGONAL ELEMENTS STORED IN DIAGA AND
C THE STRICT LOWER TRIANGULAR PART STORED BY ROWS IN LOWERA.
C BOTH X AND Y ARE VECTORS OF LENGTH N.
C THE FUNCTION IX (BELOW) GIVES THE LOCATION OF A MATRIX
C ELEMENT (I,J), I>J, IN THE ONE-DIMENSIONAL ARRAY LOWERA.
C
      INTEGER N
      DOUBLE PRECISION DIAGA(N), LOWERA(N*(N-1)/2), X(N), Y(N)
      INTEGER IX, II, JJ
      IX(II,JJ)=(II-1)*(II-2)/2+JJ

      DO 10 I = 1, N
         Y(I) = DIAGA(I) * X(I)
   10 CONTINUE
      DO 40 I = 1, N
         DO 20 J = 1, I-1
            Y(I) = Y(I) + LOWERA(IX(I,J))*X(J)
   20    CONTINUE
         DO 30 J = I+1, N
            Y(I) = Y(I) + LOWERA(IX(J,I))*X(J)
   30    CONTINUE
   40 CONTINUE
      RETURN
      END
C
C ----------------------------------------------------------
C SEGMENT 4: DRIVER AND ROUTINES FOR TESTING HESFCN (SINGLE PRECISION)
C ----------------------------------------------------------
C
      PROGRAM TESTH
C
C TESTH IS THE DRIVER PROGRAM FOR EXERCISING THE VARIOUS COMPONENTS
C OF ALGORITHM 566 WITH THE NEW HESSIAN SEGMENT, ROUTINE HESFCN.
C THE TESTING OF HESFCN IS ACCOMPLISHED THROUGH TAYLOR EXPANSIONS
C (SUBROUTINE TESTGH), WHERE THE RESULTING ERROR FROM THE SECOND-
C ORDER EXPANSION INDICATES WHETHER THE GRADIENT ONLY, OR BOTH THE
C GRADIENT AND HESSIAN, ARE CORRECT.
C
      INTEGER MAXFCN, MAXN
      REAL ZERO, ONE, FIVE
      PARAMETER (MAXFCN=18, MAXN=100, ZERO=0.E0, ONE=1.0E0, FIVE=5.0E0)
      INTEGER N, NPROB, NTRIES, NREAD, NWRITE, NVARS(MAXFCN)
      REAL FACTOR, F0, X0(MAXN), X(MAXN), Y(MAXN),
     1      H0Y(MAXN), G0(MAXN), H0D(MAXN), H0L(MAXN*(MAXN-1)/2), YHY,
     2      RANVEC(MAXN)
      EXTERNAL INITPT, OBJFCN, GRDFCN, HESFCN, MVPROD
      DATA NREAD, NWRITE /5,6/
      DATA NVARS /3,6,3,2,3,0,0,0,0,2,4,3,0,0,0,2,4,0/
      DATA RANVEC/ 0.908E0, 0.769E0, 0.734E0, 0.644E0,-0.589E0, 0.577E0,
     1   -0.786E0,-0.901E0, 0.517E0, 0.767E0, 0.749E0,-0.978E0, 0.874E0,
     2    0.777E0,-0.945E0,-0.812E0, 0.921E0, 0.580E0,-0.606E0,-0.857E0,
     3   -0.565E0,-0.545E0, 0.637E0, 0.501E0, 0.707E0, 0.513E0, 0.855E0,
     4   -0.969E0, 0.620E0,-0.590E0, 0.659E0, 0.943E0, 0.826E0,-0.575E0,
     5   -0.841E0, 0.693E0,-0.694E0, 0.750E0, 0.574E0, 0.794E0,-0.923E0,
     6   -0.795E0, 0.978E0, 0.778E0, 0.574E0, 0.992E0,-0.704E0, 0.571E0,
     7    0.782E0,-0.626E0,-0.744E0, 0.732E0,-0.981E0,-0.563E0,-0.600E0,
     8   -0.660E0,-0.815E0,-0.563E0, 0.826E0, 0.811E0,-0.902E0, 0.624E0,
     9    0.738E0, 0.695E0,-0.602E0, 0.514E0,-0.951E0,-0.713E0,-0.571E0,
     A    0.974E0,-0.705E0, 0.566E0,-0.943E0,-0.546E0, 0.581E0, 0.536E0,
     B   -0.683E0, 0.627E0,-0.568E0, 0.892E0, 0.728E0, 0.675E0,-0.726E0,
     C   -0.904E0, 0.966E0, 0.826E0, 0.608E0, 0.840E0, 0.954E0, 0.625E0,
     D    0.930E0,-0.736E0,-0.753E0,-0.800E0, 0.909E0, 0.878E0, 0.731E0,
     E   -0.976E0, 0.816E0,-0.720E0/

   10 CONTINUE
C
C READ THE NUMBER ASSOCIATED WITH THE FUNCTION TO BE USED (NPROB), THE
C NUMBER OF VARIABLES (N), AND THE NUMBER OF STARTING POINTS (NTRIES).
C (EACH LINE OF INPUT FILE MUST CONSIST OF 3 INTEGERS, ARBITRARY
C FORMAT, WITH NPROB RANGING BETWEEN 1 AND 18; N MUST BE APPROPRIATE
C FOR THE PROBLEM -- SEE INPUT FILE -- AND NTRIES MUST BE A POSITIVE 
C INTEGER. FOR FUNCTION NUMBER 18 (CHEBYQUAD), SUBROUTINES GRDFCN 
C AND HESFCN, AS CURRENTLY IMPLEMENTED, CANNOT HANDLE N GREATER THAN 50.
C ALSO, PLEASE NOTE THAT FOR ALL FUNCTIONS, N IS RESTRICTED TO BE NO
C MORE THAN 100 IN ORDER TO ALLOW EIGENVALUE CALCULATIONS RELATED TO THE
C SPECTRUM OF THE HESSIAN TO BE PERFORMED WITHOUT MAJOR DEMANDS ON
C STORAGE AND CPU TIME.  THE INPUT IS CHECKED BELOW.
C
         READ (NREAD, *) NPROB, N, NTRIES
         IF (NPROB .EQ. 0) STOP
         IF (NPROB .LT. 1 .OR. NPROB .GT. 18 .OR. N .LT. 1
     1         .OR. NTRIES .LT. 1) THEN
            WRITE (NWRITE, 850)
            STOP
         ENDIF
         IF (NVARS(NPROB) .NE. 0)
     1      N = NVARS(NPROB)
         IF (NPROB .EQ. 7)
     1      N = MAX(2, MIN(N,31))
         IF (NPROB .EQ. 14)
     1      N = MAX(2, N - MOD(N,2))
         IF (NPROB .EQ. 15)
     1      N = MAX(4, N - MOD(N,4))
         IF (NPROB .EQ. 18)
     1      N = MAX(2, MIN(N,50))
         IF (N .GT. MAXN) THEN
            WRITE (NWRITE, 875) N, MAXN
            GO TO 10
         END IF
C
C OBTAIN THE INITIAL POINT X0, AND COMPUTE THE CORRESPONDING
C FUNCTION VALUE, GRADIENT VECTOR, AND HESSIAN MATRIX
C
         FACTOR = ONE
         DO 40 K = 1, NTRIES
            CALL INITPT(N, X0, NPROB, FACTOR)
            CALL OBJFCN(N, X0, F0, NPROB)
            CALL GRDFCN(N, X0, G0, NPROB)
            CALL HESFCN(N, X0, H0D, H0L, NPROB)
            WRITE (NWRITE, 900) NPROB, N, FACTOR
C
C OBTAIN A PERTURBATION VECTOR Y
C
         DO 20 I = 1, N
            IF (X0(I) .NE. ZERO) THEN
               Y(I) = X0(I) * RANVEC(I)
            ELSE
               Y(I) = RANVEC(I)
            ENDIF
   20    CONTINUE
         WRITE (NWRITE,925)
         WRITE (NWRITE,950) (X0(I), I = 1, N)
         WRITE (NWRITE,975)
         WRITE (NWRITE,950) (Y(I), I = 1, N)
C
C COMPUTE THE INNER PRODUCT Y*HY AT X0 AND
C CALL THE DERIVATIVE TESTING FUNCTION
C
            CALL MVPROD(N, H0D, H0L, Y, H0Y)
            YHY = ZERO
            DO 30 I = 1, N
               YHY = YHY + Y(I)*H0Y(I)
   30       CONTINUE
            CALL TESTGH (N, X0, F0, G0, Y, YHY, X, NPROB)
            FACTOR = FIVE*FACTOR
   40    CONTINUE
         GO TO 10

  850 FORMAT (/4X, 'ERROR IN INPUT FILE'/)
  875 FORMAT (/4X, 'N > MAXN:  N =', I6, ', MAXN =', I6,
     +        '.  PLEASE INCREASE PARAMETER MAXN.'/)
  900 FORMAT (/4X,'TESTING FUNCTION ', I2/4X,'WITH', I5,
     *        ' VARIABLES AT THE STANDARD STARTING POINT'/
     *        4X, 'SCALED BY', 1PE16.2/)
  925 FORMAT (/4X, 'X0 VECTOR:')
  950 FORMAT (4(F14.2,2X))
  975 FORMAT (/4X, 'Y VECTOR:')
      END
C***************************************************************
      SUBROUTINE TESTGH(N,XC,FC,GC,Y,YHY,VEC,NPROB)
C
C TESTGH TESTS USER-SUPPLIED GRADIENT (G) AND HESSIAN (H)
C ROUTINES CORRESPONDING TO A GIVEN FUNCTION F.
C TESTING H IS OPTIONAL.
C
C DERIVATIVES ARE TESTED USING A TAYLOR EXPANSION OF F
C AROUND A GIVEN POINT XC. THE TAYLOR SERIES IS EXPANDED
C AT XC + EPS*Y WHERE Y IS A RANDOM  PERTURBATION VECTOR
C AND EPS IS A SCALAR. IF WE DENOTE THE DOT PRODUCT OF 2
C VECTORS A AND B AS (A,B), WE CAN WRITE OUR EXPANSION AS
C
C F(XC+EPS*Y) = F(XC) + EPS * (G,Y) + (1/2)*(EPS**2) * (Y,HY)
C                     + O(EPS**3),
C
C WHERE G AND H ARE BOTH EVALUATED AT XC, AND HY DENOTES A
C HESSIAN/VECTOR PRODUCT. IF ONLY G ROUTINES ARE TESTED, THE
C SECOND-ORDER TAYLOR TERM IS ZERO, AND THE TRUNCATION ERROR
C IS O(EPS**2).
C
C OUR TEST IS PERFORMED BY COMPUTING THIS TAYLOR APPROX. AT
C SMALLER AND SMALLER VALUES OF EPS AND CHECKING TO SEE
C WHETHER CORRECT TRUNCATION ERRORS ARE OBTAINED --
C O(EPS**2) AND  O(EPS**3) IF THE APPROX. IS CORRECT UP TO
C THE G AND H TERMS, RESPECTIVELY.
C
C WE DIVIDE EPS BY 2 AT EVERY STEP AND TEST IF INDEED THE
C TRUNCATION ERRORS DECREASE AS THEY SHOULD.
C (I.E., IF THE ERROR CORRESPONDING TO EPS IS E1,
C THE ERROR FOR EPS/2 SHOULD BE E1/4 IF THE GRADIENT
C IS CORRECT, AND E1/8 IF THE HESSIAN IS ALSO CORRECT).
C OUR VALUE "RATIO" COMPUTES THIS FACTOR OF THE OLD/NEW
C ERRORS.
C
C THE OUTPUT IS A SERIES OF VALUES FOR RATIO PRINTED FOR 
C EACH EPS UNTIL THE TRUNCATION ERROR AND/OR EPS IS VERY SMALL.
C IF RATIO TENDS TO 4 OR 8 AS EPS IS DECREASED (AND THE
C ERROR IS RELATIVELY SMALL) G  IS CORRECT OR G&H ARE CORRECT,
C RESPECTIVELY. IF RATIO TENDS TO 2, WHICH IS O(EPS), NEITHER G
C NOR H ARE CORRECT. IF THE RATIO TENDS TO 1, THE ERRORS MAY
C BE TOO LARGE GIVEN THE PERTURBATION VECTOR Y.
C
C THUS IN GENERAL, RELIABLE VALUES OF RATIO SHOULD
C OCCUR WHEN: (1) EPS IS NOT TOO LARGE AND NOT TOO SMALL,
C AND (2) THE DIFFERENCE BETWEEN F(XC+EPS*Y) AND THE
C TAYLOR SERIES APPROXIMATION IS OF REASONABLE MAGNITUDE.
C (THE VALUES OF EPS AND THE ERRORS APPEAR IN THE OUTPUT).
C IN OTHER WORDS, AN ACCURATE VALUE OF RATIO SHOULD
C APPEAR AROUND THE MIDDLE OF OUR SERIES IF Y IS APPROPRIATE.
C DIFFERENT STARTING POINT AND/OR PERTURBATION VECTORS
C CAN BE TRIED.
C
C USAGE: THE USER MUST SUPPLY THE FOLLOWING INPUT
C ------ VARIABLES IN THE FUNCTION CALL:
C
C N      - DIMENSION (NUMBER OF VARIABLES FOR F)
C XC(N)  - OUR CURRENT VECTOR
C FC     - THE FUNCTION VALUE AT XC
C GC(N)  - THE GRADIENT VECTOR AT XC, ON INPUT
C          ON OUTPUT, GC MAY BE CHANGED IF IT IS USED IN
C          THE FUNCTION CALL TO OBTAIN A NEW GRADIENT IN
C          ADDITION TO A NEW FUNCTION VALUE (SEE BELOW).
C Y(N)   - A RANDOM PERTURBATION VECTOR (Y SHOULD BE CHOSEN
C          SO THAT F(XC+Y) IS IN A REASONABLE RANGE FOR THE
C          PROBLEM)
C YHY    - THE MATRIX INNER PRODUCT -- (Y,HY) -- REPRESENTING
C          THE DOT PRODUCT OF Y WITH THE HESSIAN/VECTOR
C          PRODUCT, HY, WHERE H IS EVALUATED AT XC
C          (IF ONLY THE GRADIENT IS TESTED, SET YHY TO ZERO).
C VEC(N) - A WORK VECTOR
C NPROB  - AN INTEGER VARIABLE THAT MAY BE USED IN THE 
C          USER'S FUNCTION CALL
C
C NOTE: THE USER MAY MODIFY THIS ROUTINE TO TEST OTHER
C ----  FUNCTIONS BY REPLACING THE SAMPLE FUNCTION CALL GIVEN 
C ABOVE THE '40 CONTINUE' STATEMENT WITH THE APPROPRIATE
C CALL FOR HIS/HER PROBLEM. THE INSERTED ROUTINE CALL
C SHOULD PRODUCE A NEW FUNCTION VALUE, FVEC, FOR EACH 
C NEW VECTOR VEC=XC+EPS*Y.
C
      INTEGER N, MP
      REAL ZERO,ONE,HALF,TWOP23,EPSMCH,EPSLIM,EPS,FC,
     1                 GY,YHY,TAYLOR,DIFF,FVEC,FOLD,RATIO,TEMP,
     2                 XC(N),GC(N),Y(N),VEC(N)
      PARAMETER (ZERO=0.E0, ONE=1.0E0, HALF=0.5E0, TWOP23=8388608.0E0,
     1           MP=6)
      EXTERNAL OBJFCN

      EPSMCH = ONE / TWOP23
C
C  NOTE:  THE LINE ABOVE MAY BE REPLACED WITH:
C  -----  EPSMCH = R1MACH(3)
C  WHERE R1MACH IS A FUNCTION WHICH MAY BE OBTAINED BY SENDING AN
C  ELECTRONIC MAIL MESSAGE TO NETLIB@ORNL.GOV WITH A SUBJECT LINE OR
C  BODY OF "SEND R1MACH FROM CORE".  EPSMCH WOULD THEN HAVE THE VALUE
C  R**(-P), WHERE R IS THE RADIX OF DOUBLE PRECISION NUMBERS AND P IS
C  THE NUMBER OF RADIX-R DIGITS IN THE MANTISSA OR SIGNIFICAND.  IN
C  OTHER WORDS, EPSMCH WOULD HAVE THE VALUE OF THE SMALLEST RELATIVE
C  SPACING.  HERE, WE HAVE SIMPLY SET EPSMCH TO THE VALUE 2**(-23),
C  WHICH SHOULD GENERALLY BE LARGER THAN THE ACTUAL VALUE OF
C  R**(-P).
C
      EPSLIM = EPSMCH * FLOAT(N*N) * 1.E+2
      EPS    = HALF

      WRITE (MP,900)
      GY = ZERO
      DO 10 I = 1, N
         GY = GY + GC(I)*Y(I)
   10 CONTINUE
      WRITE (MP,910) FC,GY,YHY,EPSMCH
      WRITE (MP,940)

      DIFF = ZERO
      FVEC = FC

   20 CONTINUE
      TEMP = DIFF
      FOLD = FVEC

      DO 30 I = 1, N
         VEC(I) = XC(I) + EPS*Y(I)
   30 CONTINUE

      NOUT = 0
      CALL OBJFCN(N, VEC, FVEC, NPROB)
   40 CONTINUE

      TAYLOR = FC + (EPS*GY) + ( (EPS**2) * HALF  * YHY )
      DIFF   = FVEC - TAYLOR

      IF (ABS(DIFF) .LT. ABS(EPSLIM*FVEC)) THEN
         WRITE (MP,920) ABS(EPSLIM*FVEC)
         GOTO 50
      ENDIF

      IF (ABS(FVEC-FOLD) .LT. ABS(EPSLIM*FOLD)) THEN
         WRITE (MP,930) ABS(EPSLIM*FOLD)
         GOTO 50
      ENDIF

      IF (TEMP .EQ. ZERO .OR. DIFF .EQ. ZERO) THEN
          WRITE (MP, 950) EPS,FVEC,TAYLOR,DIFF
      ELSE
          RATIO = TEMP / DIFF
          WRITE (MP, 950) EPS,FVEC,TAYLOR,DIFF,RATIO
      ENDIF

      EPS = EPS * HALF
      IF (EPS .GT. EPSMCH) GOTO 20

   50 RETURN

  900 FORMAT(/T10,'ENTERING TESTGH ROUTINE:'/)
  910 FORMAT(T5,  'THE FUNCTION VALUE AT X               = ',
     + 1PE16.8/T5,'THE FIRST-ORDER TAYLOR TERM,  (G, Y)  = ',
     + 1PE16.8/T5,'THE SECOND-ORDER TAYLOR TERM, (Y,HY)  = ',
     + 1PE16.8//T5,'THE COMPUTED MACHINE PRECISION        = ',
     + 1PE16.8//)
  920 FORMAT(/T5,'DIFF IS SMALL (LESS THAN ', 1PE16.8,
     + ' IN ABSOLUTE VALUE)'/)
  930 FORMAT(/T5,'CHANGE IN FUNCTION VALUE IS VERY SMALL (LESS THAN ',
     + 1PE16.8,' IN ABSOLUTE VALUE)'/)
  940 FORMAT(4X,'EPS',10X,' F   ',10X,' TAYLOR',9X,
     +  ' DIFF.',9X,'RATIO'/)
  950 FORMAT(1X,1PE10.4,1PE16.8,1PE16.8,1PE16.8,1PE16.8)

      END
C***********************************************************************
      SUBROUTINE MVPROD (N, DIAGA, LOWERA, X, Y)
C
C MVPROD PERFORMS THE MATRIX-VECTOR PRODUCT A*X, AND
C STORES THE RESULT IN THE VECTOR Y.  A IS A SYMMETRIC
C NXN MATRIX, WITH DIAGONAL ELEMENTS STORED IN DIAGA AND
C THE STRICT LOWER TRIANGULAR PART STORED BY ROWS IN LOWERA.
C BOTH X AND Y ARE VECTORS OF LENGTH N.
C THE FUNCTION IX (BELOW) GIVES THE LOCATION OF A MATRIX
C ELEMENT (I,J), I>J, IN THE ONE-DIMENSIONAL ARRAY LOWERA.
C
      INTEGER N
      REAL DIAGA(N), LOWERA(N*(N-1)/2), X(N), Y(N)
      INTEGER IX, II, JJ
      IX(II,JJ)=(II-1)*(II-2)/2+JJ

      DO 10 I = 1, N
         Y(I) = DIAGA(I) * X(I)
   10 CONTINUE
      DO 40 I = 1, N
         DO 20 J = 1, I-1
            Y(I) = Y(I) + LOWERA(IX(I,J))*X(J)
   20    CONTINUE
         DO 30 J = I+1, N
            Y(I) = Y(I) + LOWERA(IX(J,I))*X(J)
   30    CONTINUE
   40 CONTINUE
      RETURN
      END
C************************************************************
C ----------------------------------------------------------
C SEGMENT 5: INPUT FILE (COMMENTED)
C ----------------------------------------------------------
C EACH LINE CONTAINS 3 INTEGERS THAT ARE READ BY THE
C TESTING PROBLEM: {NPROB,N,NTRIES}, WHERE NPROB IS
C THE PROBLEM NUMBER (1 TO 18), N IS THE NUMBER OF
C VARIABLES (SEE BELOW THE PERMITTED RANGE FOR EACH
C PROBLEM), AND NTRIES IS THE NUMBER OF TIMES
C TESTS WILL BE MADE FOR THE FUNCTION BY SCALING X0.
C FIRST LINE OF FILE SHOULD BEGIN WITH THE FIRST TRIPLET,
C AND LAST LINE SHOULD BE {0,0,0}. THE FORMAT IS FLEXIBLE.
C 
C   1    3    2      Helical Valley, N=3
C   2    6    2      Biggs EXP6, N=6
C   3    3    2      Gaussian, N=3
C   4    2    3      Powell Badly Scaled, N=2
C   5    3    2      Box 3D, N=3
C   6    4    4      Variably Dimensioned, N variable
C   7   10    4      Watson, 2<=N<=31
C   8    2    2      Penalty I, N variable
C   9    2    2      Penalty II, N variable
C  10    2    5      Brown Badly Scaled, N=2
C  11    4    2      Brown and Dennis, N=4
C  12    3    2      Gulf Research and Development, N=3
C  13    2    2      Trigonometric, N variable
C  14   12    5      Extended Rosebrock, N variable, even
C  15    4    2      Extended Powell Singular, N multiple of 4
C  16    2    2      Beale, N=2
C  17    4    4      Wood, N=4
C  18    2    2      Chebyquad, N variable
C   0    0    0

---------------------------------------------------------------------

C ALGORITHM 566
C
C FORTRAN SUBROUTINES FOR TESTING UNCONSTRAINED OPTIMIZATION
C SOFTWARE
C
C BY J.J. MORE, B.S. GARBOW AND K.E. HILLSTROM
C
C ACM TRANSACTIONS ON MATHEMATICAL SOFTWARE 7,1 (MARCH 1981)
C
C ===== THERE ARE 16 PARTS TO THIS FILE
C ===== 1. DOCUMENTATION.
C ===== 2. DOUBLE PRECISION TESTING AIDS FOR NONLINEAR EQUATIONS.
C ===== 3. DOUBLE PRECISION TESTING AIDS FOR NONLINEAR LEAST-SQUARES.
C ===== 4. DOUBLE PRECISION TESTING AIDS FOR UNCONSTRAINED NONLINEAR
C =====     OPTIMIZATION.
C ===== 5. SINGLE PRECISION TESTING AIDS FOR NONLINEAR EQUATIONS.
C ===== 6. SINGLE PRECISION TESTING AIDS FOR NONLINEAR LEAST-SQUARES.
C ===== 7. SINGLE PRECISION TESTING AIDS FOR UNCONSTRAINED NONLINEAR
C =====     OPTIMIZATION.
C ===== 8. SAMPLE DRIVER FOR DOUBLE PRECISION NONLINEAR EQUATIONS.
C ===== 9. SAMPLE DRIVER FOR SINGLE PRECISION NONLINEAR EQUATIONS.
C ===== 10. SAMPLE DRIVER FOR DOUBLE PRECISION NONLINEAR LEAST-SQUARES.
C ===== 11. SAMPLE DRIVER FOR SINGLE PRECISION NONLINEAR LEAST-SQUARES.
C ===== 12. SAMPLE DRIVER FOR DOUBLE PRECISION UNCONSTRAINED NONLINEAR
C =====     MINIMIZATION.
C ===== 13. SAMPLE DRIVER FOR SINGLE PRECISION UNCONSTRAINED NONLINEAR
C =====     MINIMIZATION.
C ===== 14. DATA (NONLINEAR EQUATIONS).
C ===== 15. DATA (NONLINEARR LEAST SQUARES).
C ===== 16. DATA (UNCONSTRAINED NONLINEAR OPTIMIZATION).
C =====
C =====
C =====
C ===== 1. DOCUMENTATION.
 DESCRIPTION

      This is the Fortran package of subroutines described in (1)
 for testing unconstrained optimization software.  The following
 three problem areas are considered.

      1.  Zeros of systems of N nonlinear functions in N variables.

      2.  Least Squares minimization of M nonlinear functions in
          N variables.

      3.  Unconstrained minimization of an objective function with
          N variables.

      The subroutines which define the test functions and starting
 points depend on the dimension parameters M and N and on the
 problem number NPROB.  We first describe the subroutines for the
 test functions.

      For systems of nonlinear functions,

                VECFCN(N,X,FVEC,NPROB)

 returns the function values in the N-vector FVEC, and

                VECJAC(N,X,FJAC,LDFJAC,NPROB)

 returns the Jacobian matrix in the N by N array FJAC.  (The parameter
 LDFJAC is the leading dimension of the array FJAC as defined in the
 main program.)  In order to prevent gross inefficiencies with solvers
 which only require one function value at a time,

                COMFCN(N,K,X,FCNK,NPROB)

 returns the K-th function value in FCNK.

      For nonlinear least squares,

                SSQFCN(M,N,X,FVEC,NPROB)

 returns the function values in the M-vector FVEC, and

                SSQJAC(M,N,X,FJAC,LDFJAC,NPROB)

 returns the Jacobian matrix in the M by N array FJAC.

      For unconstrained minimization,

                OBJFCN(N,X,F,NPROB)

 returns the objective function value in F, and

                GRDFCN(N,X,G,NPROB)

 returns the gradient components in the N-vector G.

      For each problem area, the starting points are generated by

                INITPT(N,X,NPROB,FACTOR)

 which returns in X the starting point corresponding to the
 parameters NPROB and FACTOR.  If XS denotes the standard starting
 point, then X will contain FACTOR*XS, except that if XS is the
 zero vector and FACTOR is not unity, then all the components of X
 will be set to FACTOR.

      To test a code in any of the three problem areas, the user
 must provide a driver and interface routine.  The driver reads in
 the data which defines the dimensions, the problem number, and
 FACTOR, calls INITPT, and then calls the code of interest and
 prints out results.  The interface routine provides a link between
 the code with its particular function routine calling sequences
 and the subroutines for the test functions.

      The package includes example drivers and interface routines
 for each of the problem areas.  Sample data is also provided.


 REFERENCES

 1. More, J.J., Garbow, B.S., and Hillstrom, K.E., Testing
    Unconstrained Optimization Software, ACM Trans. Math. Software
    (this issue).
C ===== 2. DOUBLE PRECISION TESTING AIDS FOR NONLINEAR EQUATIONS.
      SUBROUTINE INITPT(N,X,NPROB,FACTOR)                               00000010
      INTEGER N,NPROB
      DOUBLE PRECISION FACTOR
      DOUBLE PRECISION X(N)
C     **********
C
C     SUBROUTINE INITPT
C
C     THIS SUBROUTINE SPECIFIES THE STANDARD STARTING POINTS FOR
C     THE FUNCTIONS DEFINED BY SUBROUTINES COMFCN AND VECFCN. THE
C     SUBROUTINE RETURNS IN X A MULTIPLE (FACTOR) OF THE STANDARD
C     STARTING POINT. FOR THE SIXTH FUNCTION THE STANDARD STARTING
C     POINT IS ZERO, SO IN THIS CASE, IF FACTOR IS NOT UNITY, THEN THE
C     SUBROUTINE RETURNS THE VECTOR  X(J) = FACTOR, J=1,...,N.
C
C     THE SUBROUTINE STATEMENT IS
C
C       SUBROUTINE INITPT(N,X,NPROB,FACTOR)
C
C     WHERE
C
C       N IS A POSITIVE INTEGER INPUT VARIABLE.
C
C       X IS AN OUTPUT ARRAY OF LENGTH N WHICH CONTAINS THE STANDARD
C         STARTING POINT FOR PROBLEM NPROB MULTIPLIED BY FACTOR.
C
C       NPROB IS A POSITIVE INTEGER INPUT VARIABLE WHICH DEFINES THE
C         NUMBER OF THE PROBLEM. NPROB MUST NOT EXCEED 14.
C
C       FACTOR IS AN INPUT VARIABLE WHICH SPECIFIES THE MULTIPLE OF
C         THE STANDARD STARTING POINT. IF FACTOR IS UNITY, NO
C         MULTIPLICATION IS PERFORMED.
C
C     MINPACK. VERSION OF JULY 1978.
C     BURTON S. GARBOW, KENNETH E. HILLSTROM, JORGE J. MORE
C
C     **********
      INTEGER IVAR,J
      DOUBLE PRECISION C1,H,HALF,ONE,THREE,TJ,ZERO
      DOUBLE PRECISION DFLOAT
      DATA ZERO,HALF,ONE,THREE,C1 /0.0D0,5.0D-1,1.0D0,3.0D0,1.2D0/
      DFLOAT(IVAR) = IVAR
C
C     SELECTION OF INITIAL POINT.
C
      GO TO (10,20,30,40,50,60,80,100,120,120,140,160,180,180), NPROB
C
C     ROSENBROCK FUNCTION.
C
   10 CONTINUE
      X(1) = -C1
      X(2) = ONE
      GO TO 200
C
C     POWELL SINGULAR FUNCTION.
C
   20 CONTINUE
      X(1) = THREE
      X(2) = -ONE
      X(3) = ZERO
      X(4) = ONE
      GO TO 200
C
C     POWELL BADLY SCALED FUNCTION.
C
   30 CONTINUE
      X(1) = ZERO
      X(2) = ONE
      GO TO 200
C
C     WOOD FUNCTION.
C
   40 CONTINUE
      X(1) = -THREE
      X(2) = -ONE
      X(3) = -THREE
      X(4) = -ONE
      GO TO 200
C
C     HELICAL VALLEY FUNCTION.
C
   50 CONTINUE
      X(1) = -ONE
      X(2) = ZERO
      X(3) = ZERO
      GO TO 200
C
C     WATSON FUNCTION.
C
   60 CONTINUE
      DO 70 J = 1, N
         X(J) = ZERO
   70    CONTINUE
      GO TO 200
C
C     CHEBYQUAD FUNCTION.
C
   80 CONTINUE
      H = ONE/DFLOAT(N+1)
      DO 90 J = 1, N
         X(J) = DFLOAT(J)*H
   90    CONTINUE
      GO TO 200
C
C     BROWN ALMOST-LINEAR FUNCTION.
C
  100 CONTINUE
      DO 110 J = 1, N
         X(J) = HALF
  110    CONTINUE
      GO TO 200
C
C     DISCRETE BOUNDARY VALUE AND INTEGRAL EQUATION FUNCTIONS.
C
  120 CONTINUE
      H = ONE/DFLOAT(N+1)
      DO 130 J = 1, N
         TJ = DFLOAT(J)*H
         X(J) = TJ*(TJ - ONE)
  130    CONTINUE
      GO TO 200
C
C     TRIGONOMETRIC FUNCTION.
C
  140 CONTINUE
      H = ONE/DFLOAT(N)
      DO 150 J = 1, N
         X(J) = H
  150    CONTINUE
      GO TO 200
C
C     VARIABLY DIMENSIONED FUNCTION.
C
  160 CONTINUE
      H = ONE/DFLOAT(N)
      DO 170 J = 1, N
         X(J) = ONE - DFLOAT(J)*H
  170    CONTINUE
      GO TO 200
C
C     BROYDEN TRIDIAGONAL AND BANDED FUNCTIONS.
C
  180 CONTINUE
      DO 190 J = 1, N
         X(J) = -ONE
  190    CONTINUE
  200 CONTINUE
C
C     COMPUTE MULTIPLE OF INITIAL POINT.
C
      IF (FACTOR .EQ. ONE) GO TO 250
      IF (NPROB .EQ. 6) GO TO 220
         DO 210 J = 1, N
            X(J) = FACTOR*X(J)
  210       CONTINUE
         GO TO 240
  220 CONTINUE
         DO 230 J = 1, N
            X(J) = FACTOR
  230       CONTINUE
  240 CONTINUE
  250 CONTINUE
      RETURN
C
C     LAST CARD OF SUBROUTINE INITPT.
C
      END
      SUBROUTINE COMFCN(N,K,X,FCNK,NPROB)                               00000010
      INTEGER N,K,NPROB
      DOUBLE PRECISION FCNK
      DOUBLE PRECISION X(N)
C     **********
C
C     SUBROUTINE COMFCN
C
C     THIS SUBROUTINE DEFINES FOURTEEN TEST FUNCTIONS. THE FIRST
C     FIVE TEST FUNCTIONS ARE OF DIMENSIONS 2,4,2,4,3, RESPECTIVELY,
C     WHILE THE REMAINING TEST FUNCTIONS ARE OF VARIABLE DIMENSION
C     N FOR ANY N GREATER THAN OR EQUAL TO 1 (PROBLEM 6 IS AN
C     EXCEPTION TO THIS, SINCE IT DOES NOT ALLOW N = 1).
C
C     THE SUBROUTINE STATEMENT IS
C
C       SUBROUTINE COMFCN(N,K,X,FCNK,NPROB)
C
C     WHERE
C
C       N IS A POSITIVE INTEGER INPUT VARIABLE.
C
C       K IS A POSITIVE INTEGER INPUT VARIABLE NOT GREATER THAN N.
C
C       X IS AN INPUT ARRAY OF LENGTH N.
C
C       FCNK IS AN OUTPUT VARIABLE WHICH CONTAINS THE VALUE OF
C         THE K-TH COMPONENT OF THE NPROB FUNCTION EVALUATED AT X.
C
C       NPROB IS A POSITIVE INTEGER INPUT VARIABLE WHICH DEFINES THE
C         NUMBER OF THE PROBLEM. NPROB MUST NOT EXCEED 14.
C
C     SUBPROGRAMS REQUIRED
C
C       FORTRAN-SUPPLIED ... DATAN,DCOS,DEXP,DSIGN,DSIN,DSQRT,
C                            MAX0,MIN0,MOD
C
C     MINPACK. VERSION OF JULY 1978.
C     BURTON S. GARBOW, KENNETH E. HILLSTROM, JORGE J. MORE
C
C     **********
      INTEGER I,IVAR,J,K1,K2,KP1,ML,MU
      DOUBLE PRECISION C1,C2,C3,C4,C5,C6,C7,C8,C9,EIGHT,FIVE,H,ONE,
     1                 PROD,SUM,SUM1,SUM2,TEMP,TEMP1,TEMP2,TEN,THREE,
     2                 TI,TJ,TK,TPI,TWO,ZERO
      DOUBLE PRECISION DFLOAT
      DATA ZERO,ONE,TWO,THREE,FIVE,EIGHT,TEN
     1     /0.0D0,1.0D0,2.0D0,3.0D0,5.0D0,8.0D0,1.0D1/
      DATA C1,C2,C3,C4,C5,C6,C7,C8,C9
     1     /1.0D4,1.0001D0,2.0D2,2.02D1,1.98D1,1.8D2,2.5D-1,5.0D-1,
     2      2.9D1/
      DFLOAT(IVAR) = IVAR
C
C     PROBLEM SELECTOR.
C
      GO TO (10,20,30,40,50,70,110,150,200,210,250,270,290,300), NPROB
C
C     ROSENBROCK FUNCTION.
C
   10 CONTINUE
      IF (K .EQ. 1) FCNK = ONE - X(1)
      IF (K .EQ. 2) FCNK = TEN*(X(2) - X(1)**2)
      GO TO 320
C
C     POWELL SINGULAR FUNCTION.
C
   20 CONTINUE
      IF (K .EQ. 1) FCNK = X(1) + TEN*X(2)
      IF (K .EQ. 2) FCNK = DSQRT(FIVE)*(X(3) - X(4))
      IF (K .EQ. 3) FCNK = (X(2) - TWO*X(3))**2
      IF (K .EQ. 4) FCNK = DSQRT(TEN)*(X(1) - X(4))**2
      GO TO 320
C
C     POWELL BADLY SCALED FUNCTION.
C
   30 CONTINUE
      IF (K .EQ. 1) FCNK = C1*X(1)*X(2) - ONE
      IF (K .EQ. 2) FCNK = DEXP(-X(1)) + DEXP(-X(2)) - C2
      GO TO 320
C
C     WOOD FUNCTION.
C
   40 CONTINUE
      TEMP1 = X(2) - X(1)**2
      TEMP2 = X(4) - X(3)**2
      IF (K .EQ. 1) FCNK = -C3*X(1)*TEMP1 - (ONE - X(1))
      IF (K .EQ. 2)
     1   FCNK = C3*TEMP1 + C4*(X(2) - ONE) + C5*(X(4) - ONE)
      IF (K .EQ. 3) FCNK = -C6*X(3)*TEMP2 - (ONE - X(3))
      IF (K .EQ. 4)
     1   FCNK = C6*TEMP2 + C4*(X(4) - ONE) + C5*(X(2) - ONE)
      GO TO 320
C
C     HELICAL VALLEY FUNCTION.
C
   50 CONTINUE
      IF (K .NE. 1) GO TO 60
      TPI = EIGHT*DATAN(ONE)
      TEMP1 = DSIGN(C7,X(2))
      IF (X(1) .GT. ZERO) TEMP1 = DATAN(X(2)/X(1))/TPI
      IF (X(1) .LT. ZERO) TEMP1 = DATAN(X(2)/X(1))/TPI + C8
      FCNK = TEN*(X(3) - TEN*TEMP1)
   60 CONTINUE
      IF (K .EQ. 2) FCNK = TEN*(DSQRT(X(1)**2+X(2)**2) - ONE)
      IF (K .EQ. 3) FCNK = X(3)
      GO TO 320
C
C     WATSON FUNCTION.
C
   70 CONTINUE
      FCNK = ZERO
      DO 100 I = 1, 29
         TI = DFLOAT(I)/C9
         SUM1 = ZERO
         TEMP = ONE
         DO 80 J = 2, N
            SUM1 = SUM1 + DFLOAT(J-1)*TEMP*X(J)
            TEMP = TI*TEMP
   80       CONTINUE
         SUM2 = ZERO
         TEMP = ONE
         DO 90 J = 1, N
            SUM2 = SUM2 + TEMP*X(J)
            TEMP = TI*TEMP
   90       CONTINUE
         TEMP1 = SUM1 - SUM2**2 - ONE
         TEMP2 = TWO*TI*SUM2
         FCNK = FCNK + TI**(K - 2)*(DFLOAT(K-1) - TEMP2)*TEMP1
  100    CONTINUE
      TEMP = X(2) - X(1)**2 - ONE
      IF (K .EQ. 1) FCNK = FCNK + X(1)*(ONE - TWO*TEMP)
      IF (K .EQ. 2) FCNK = FCNK + TEMP
      GO TO 320
C
C     CHEBYQUAD FUNCTION.
C
  110 CONTINUE
      SUM = ZERO
      DO 140 J = 1, N
         TEMP1 = ONE
         TEMP2 = TWO*X(J) - ONE
         TEMP = TWO*TEMP2
         IF (K .LT. 2) GO TO 130
         DO 120 I = 2, K
            TI = TEMP*TEMP2 - TEMP1
            TEMP1 = TEMP2
            TEMP2 = TI
  120       CONTINUE
  130    CONTINUE
         SUM = SUM + TEMP2
  140    CONTINUE
      FCNK = SUM/DFLOAT(N)
      IF (MOD(K,2) .EQ. 0) FCNK = FCNK + ONE/(DFLOAT(K)**2 - ONE)
      GO TO 320
C
C     BROWN ALMOST-LINEAR FUNCTION.
C
  150 CONTINUE
      IF (K .EQ. N) GO TO 170
         SUM = -DFLOAT(N+1)
         DO 160 J = 1, N
            SUM = SUM + X(J)
  160       CONTINUE
         FCNK = X(K) + SUM
         GO TO 190
  170 CONTINUE
         PROD = ONE
         DO 180 J = 1, N
            PROD = X(J)*PROD
  180       CONTINUE
         FCNK = PROD - ONE
  190 CONTINUE
      GO TO 320
C
C     DISCRETE BOUNDARY VALUE FUNCTION.
C
  200 CONTINUE
      H = ONE/DFLOAT(N+1)
      TEMP = (X(K) + DFLOAT(K)*H + ONE)**3
      TEMP1 = ZERO
      IF (K .NE. 1) TEMP1 = X(K-1)
      TEMP2 = ZERO
      IF (K .NE. N) TEMP2 = X(K+1)
      FCNK = TWO*X(K) - TEMP1 - TEMP2 + TEMP*H**2/TWO
      GO TO 320
C
C     DISCRETE INTEGRAL EQUATION FUNCTION.
C
  210 CONTINUE
      H = ONE/DFLOAT(N+1)
      TK = DFLOAT(K)*H
      SUM1 = ZERO
      DO 220 J = 1, K
         TJ = DFLOAT(J)*H
         TEMP = (X(J) + TJ + ONE)**3
         SUM1 = SUM1 + TJ*TEMP
  220    CONTINUE
      SUM2 = ZERO
      KP1 = K + 1
      IF (N .LT. KP1) GO TO 240
      DO 230 J = KP1, N
         TJ = DFLOAT(J)*H
         TEMP = (X(J) + TJ + ONE)**3
         SUM2 = SUM2 + (ONE - TJ)*TEMP
  230    CONTINUE
  240 CONTINUE
      FCNK = X(K) + H*((ONE - TK)*SUM1 + TK*SUM2)/TWO
      GO TO 320
C
C     TRIGONOMETRIC FUNCTION.
C
  250 CONTINUE
      SUM = ZERO
      DO 260 J = 1, N
         SUM = SUM + DCOS(X(J))
  260    CONTINUE
      FCNK = DFLOAT(N+K) - DSIN(X(K)) - SUM - DFLOAT(K)*DCOS(X(K))
      GO TO 320
C
C     VARIABLY DIMENSIONED FUNCTION.
C
  270 CONTINUE
      SUM = ZERO
      DO 280 J = 1, N
         SUM = SUM + DFLOAT(J)*(X(J) - ONE)
  280    CONTINUE
      TEMP = SUM*(ONE + TWO*SUM**2)
      FCNK = X(K) - ONE + DFLOAT(K)*TEMP
      GO TO 320
C
C     BROYDEN TRIDIAGONAL FUNCTION.
C
  290 CONTINUE
      TEMP = (THREE - TWO*X(K))*X(K)
      TEMP1 = ZERO
      IF (K .NE. 1) TEMP1 = X(K-1)
      TEMP2 = ZERO
      IF (K .NE. N) TEMP2 = X(K+1)
      FCNK = TEMP - TEMP1 - TWO*TEMP2 + ONE
      GO TO 320
C
C     BROYDEN BANDED FUNCTION.
C
  300 CONTINUE
      ML = 5
      MU = 1
      K1 = MAX0(1,K-ML)
      K2 = MIN0(K+MU,N)
      TEMP = ZERO
      DO 310 J = K1, K2
         IF (J .NE. K) TEMP = TEMP + X(J)*(ONE + X(J))
  310    CONTINUE
      FCNK = X(K)*(TWO + FIVE*X(K)**2) + ONE - TEMP
  320 CONTINUE
      RETURN
C
C     LAST CARD OF SUBROUTINE COMFCN.
C
      END
      SUBROUTINE VECFCN(N,X,FVEC,NPROB)                                 00000010
      INTEGER N,NPROB
      DOUBLE PRECISION X(N),FVEC(N)
C     **********
C
C     SUBROUTINE VECFCN
C
C     THIS SUBROUTINE DEFINES FOURTEEN TEST FUNCTIONS. THE FIRST
C     FIVE TEST FUNCTIONS ARE OF DIMENSIONS 2,4,2,4,3, RESPECTIVELY,
C     WHILE THE REMAINING TEST FUNCTIONS ARE OF VARIABLE DIMENSION
C     N FOR ANY N GREATER THAN OR EQUAL TO 1 (PROBLEM 6 IS AN
C     EXCEPTION TO THIS, SINCE IT DOES NOT ALLOW N = 1).
C
C     THE SUBROUTINE STATEMENT IS
C
C       SUBROUTINE VECFCN(N,X,FVEC,NPROB)
C
C     WHERE
C
C       N IS A POSITIVE INTEGER INPUT VARIABLE.
C
C       X IS AN INPUT ARRAY OF LENGTH N.
C
C       FVEC IS AN OUTPUT ARRAY OF LENGTH N WHICH CONTAINS THE NPROB
C         FUNCTION VECTOR EVALUATED AT X.
C
C       NPROB IS A POSITIVE INTEGER INPUT VARIABLE WHICH DEFINES THE
C         NUMBER OF THE PROBLEM. NPROB MUST NOT EXCEED 14.
C
C     SUBPROGRAMS CALLED
C
C       FORTRAN-SUPPLIED ... DATAN,DCOS,DEXP,DSIGN,DSIN,DSQRT,
C                            MAX0,MIN0
C
C     MINPACK. VERSION OF JULY 1978.
C     BURTON S. GARBOW, KENNETH E. HILLSTROM, JORGE J. MORE
C
C     **********
      INTEGER I,IEV,IVAR,J,K,K1,K2,KP1,ML,MU
      DOUBLE PRECISION C1,C2,C3,C4,C5,C6,C7,C8,C9,EIGHT,FIVE,H,ONE,
     1                 PROD,SUM,SUM1,SUM2,TEMP,TEMP1,TEMP2,TEN,THREE,
     2                 TI,TJ,TK,TPI,TWO,ZERO
      DOUBLE PRECISION DFLOAT
      DATA ZERO,ONE,TWO,THREE,FIVE,EIGHT,TEN
     1     /0.0D0,1.0D0,2.0D0,3.0D0,5.0D0,8.0D0,1.0D1/
      DATA C1,C2,C3,C4,C5,C6,C7,C8,C9
     1     /1.0D4,1.0001D0,2.0D2,2.02D1,1.98D1,1.8D2,2.5D-1,5.0D-1,
     2      2.9D1/
      DFLOAT(IVAR) = IVAR
C
C     PROBLEM SELECTOR.
C
      GO TO (10,20,30,40,50,60,120,170,200,220,270,300,330,350), NPROB
C
C     ROSENBROCK FUNCTION.
C
   10 CONTINUE
      FVEC(1) = ONE - X(1)
      FVEC(2) = TEN*(X(2) - X(1)**2)
      GO TO 380
C
C     POWELL SINGULAR FUNCTION.
C
   20 CONTINUE
      FVEC(1) = X(1) + TEN*X(2)
      FVEC(2) = DSQRT(FIVE)*(X(3) - X(4))
      FVEC(3) = (X(2) - TWO*X(3))**2
      FVEC(4) = DSQRT(TEN)*(X(1) - X(4))**2
      GO TO 380
C
C     POWELL BADLY SCALED FUNCTION.
C
   30 CONTINUE
      FVEC(1) = C1*X(1)*X(2) - ONE
      FVEC(2) = DEXP(-X(1)) + DEXP(-X(2)) - C2
      GO TO 380
C
C     WOOD FUNCTION.
C
   40 CONTINUE
      TEMP1 = X(2) - X(1)**2
      TEMP2 = X(4) - X(3)**2
      FVEC(1) = -C3*X(1)*TEMP1 - (ONE - X(1))
      FVEC(2) = C3*TEMP1 + C4*(X(2) - ONE) + C5*(X(4) - ONE)
      FVEC(3) = -C6*X(3)*TEMP2 - (ONE - X(3))
      FVEC(4) = C6*TEMP2 + C4*(X(4) - ONE) + C5*(X(2) - ONE)
      GO TO 380
C
C     HELICAL VALLEY FUNCTION.
C
   50 CONTINUE
      TPI = EIGHT*DATAN(ONE)
      TEMP1 = DSIGN(C7,X(2))
      IF (X(1) .GT. ZERO) TEMP1 = DATAN(X(2)/X(1))/TPI
      IF (X(1) .LT. ZERO) TEMP1 = DATAN(X(2)/X(1))/TPI + C8
      TEMP2 = DSQRT(X(1)**2+X(2)**2)
      FVEC(1) = TEN*(X(3) - TEN*TEMP1)
      FVEC(2) = TEN*(TEMP2 - ONE)
      FVEC(3) = X(3)
      GO TO 380
C
C     WATSON FUNCTION.
C
   60 CONTINUE
      DO 70 K = 1, N
         FVEC(K) = ZERO
   70    CONTINUE
      DO 110 I = 1, 29
         TI = DFLOAT(I)/C9
         SUM1 = ZERO
         TEMP = ONE
         DO 80 J = 2, N
            SUM1 = SUM1 + DFLOAT(J-1)*TEMP*X(J)
            TEMP = TI*TEMP
   80       CONTINUE
         SUM2 = ZERO
         TEMP = ONE
         DO 90 J = 1, N
            SUM2 = SUM2 + TEMP*X(J)
            TEMP = TI*TEMP
   90       CONTINUE
         TEMP1 = SUM1 - SUM2**2 - ONE
         TEMP2 = TWO*TI*SUM2
         TEMP = ONE/TI
         DO 100 K = 1, N
            FVEC(K) = FVEC(K) + TEMP*(DFLOAT(K-1) - TEMP2)*TEMP1
            TEMP = TI*TEMP
  100       CONTINUE
  110    CONTINUE
      TEMP = X(2) - X(1)**2 - ONE
      FVEC(1) = FVEC(1) + X(1)*(ONE - TWO*TEMP)
      FVEC(2) = FVEC(2) + TEMP
      GO TO 380
C
C     CHEBYQUAD FUNCTION.
C
  120 CONTINUE
      DO 130 K = 1, N
         FVEC(K) = ZERO
  130    CONTINUE
      DO 150 J = 1, N
         TEMP1 = ONE
         TEMP2 = TWO*X(J) - ONE
         TEMP = TWO*TEMP2
         DO 140 I = 1, N
            FVEC(I) = FVEC(I) + TEMP2
            TI = TEMP*TEMP2 - TEMP1
            TEMP1 = TEMP2
            TEMP2 = TI
  140       CONTINUE
  150    CONTINUE
      TK = ONE/DFLOAT(N)
      IEV = -1
      DO 160 K = 1, N
         FVEC(K) = TK*FVEC(K)
         IF (IEV .GT. 0) FVEC(K) = FVEC(K) + ONE/(DFLOAT(K)**2 - ONE)
         IEV = -IEV
  160    CONTINUE
      GO TO 380
C
C     BROWN ALMOST-LINEAR FUNCTION.
C
  170 CONTINUE
      SUM = -DFLOAT(N+1)
      PROD = ONE
      DO 180 J = 1, N
         SUM = SUM + X(J)
         PROD = X(J)*PROD
  180    CONTINUE
      DO 190 K = 1, N
         FVEC(K) = X(K) + SUM
  190    CONTINUE
      FVEC(N) = PROD - ONE
      GO TO 380
C
C     DISCRETE BOUNDARY VALUE FUNCTION.
C
  200 CONTINUE
      H = ONE/DFLOAT(N+1)
      DO 210 K = 1, N
         TEMP = (X(K) + DFLOAT(K)*H + ONE)**3
         TEMP1 = ZERO
         IF (K .NE. 1) TEMP1 = X(K-1)
         TEMP2 = ZERO
         IF (K .NE. N) TEMP2 = X(K+1)
         FVEC(K) = TWO*X(K) - TEMP1 - TEMP2 + TEMP*H**2/TWO
  210    CONTINUE
      GO TO 380
C
C     DISCRETE INTEGRAL EQUATION FUNCTION.
C
  220 CONTINUE
      H = ONE/DFLOAT(N+1)
      DO 260 K = 1, N
         TK = DFLOAT(K)*H
         SUM1 = ZERO
         DO 230 J = 1, K
            TJ = DFLOAT(J)*H
            TEMP = (X(J) + TJ + ONE)**3
            SUM1 = SUM1 + TJ*TEMP
  230       CONTINUE
         SUM2 = ZERO
         KP1 = K + 1
         IF (N .LT. KP1) GO TO 250
         DO 240 J = KP1, N
            TJ = DFLOAT(J)*H
            TEMP = (X(J) + TJ + ONE)**3
            SUM2 = SUM2 + (ONE - TJ)*TEMP
  240       CONTINUE
  250    CONTINUE
         FVEC(K) = X(K) + H*((ONE - TK)*SUM1 + TK*SUM2)/TWO
  260    CONTINUE
      GO TO 380
C
C     TRIGONOMETRIC FUNCTION.
C
  270 CONTINUE
      SUM = ZERO
      DO 280 J = 1, N
         FVEC(J) = DCOS(X(J))
         SUM = SUM + FVEC(J)
  280    CONTINUE
      DO 290 K = 1, N
         FVEC(K) = DFLOAT(N+K) - DSIN(X(K)) - SUM - DFLOAT(K)*FVEC(K)
  290    CONTINUE
      GO TO 380
C
C     VARIABLY DIMENSIONED FUNCTION.
C
  300 CONTINUE
      SUM = ZERO
      DO 310 J = 1, N
         SUM = SUM + DFLOAT(J)*(X(J) - ONE)
  310    CONTINUE
      TEMP = SUM*(ONE + TWO*SUM**2)
      DO 320 K = 1, N
         FVEC(K) = X(K) - ONE + DFLOAT(K)*TEMP
  320    CONTINUE
      GO TO 380
C
C     BROYDEN TRIDIAGONAL FUNCTION.
C
  330 CONTINUE
      DO 340 K = 1, N
         TEMP = (THREE - TWO*X(K))*X(K)
         TEMP1 = ZERO
         IF (K .NE. 1) TEMP1 = X(K-1)
         TEMP2 = ZERO
         IF (K .NE. N) TEMP2 = X(K+1)
         FVEC(K) = TEMP - TEMP1 - TWO*TEMP2 + ONE
  340    CONTINUE
      GO TO 380
C
C     BROYDEN BANDED FUNCTION.
C
  350 CONTINUE
      ML = 5
      MU = 1
      DO 370 K = 1, N
         K1 = MAX0(1,K-ML)
         K2 = MIN0(K+MU,N)
         TEMP = ZERO
         DO 360 J = K1, K2
            IF (J .NE. K) TEMP = TEMP + X(J)*(ONE + X(J))
  360       CONTINUE
         FVEC(K) = X(K)*(TWO + FIVE*X(K)**2) + ONE - TEMP
  370    CONTINUE
  380 CONTINUE
      RETURN
C
C     LAST CARD OF SUBROUTINE VECFCN.
C
      END
      SUBROUTINE VECJAC(N,X,FJAC,LDFJAC,NPROB)                          00000010
      INTEGER N,LDFJAC,NPROB
      DOUBLE PRECISION X(N),FJAC(LDFJAC,N)
C     **********
C
C     SUBROUTINE VECJAC
C
C     THIS SUBROUTINE DEFINES THE JACOBIAN MATRICES OF FOURTEEN
C     TEST FUNCTIONS. THE PROBLEM DIMENSIONS ARE AS DESCRIBED
C     IN THE PROLOGUE COMMENTS OF VECFCN.
C
C     THE SUBROUTINE STATEMENT IS
C
C       SUBROUTINE VECJAC(N,X,FJAC,LDFJAC,NPROB)
C
C     WHERE
C
C       N IS A POSITIVE INTEGER VARIABLE.
C
C       X IS AN ARRAY OF LENGTH N.
C
C       FJAC IS AN N BY N ARRAY. ON OUTPUT FJAC CONTAINS THE
C         JACOBIAN MATRIX OF THE NPROB FUNCTION EVALUATED AT X.
C
C       LDFJAC IS A POSITIVE INTEGER VARIABLE NOT LESS THAN N
C         WHICH SPECIFIES THE LEADING DIMENSION OF THE ARRAY FJAC.
C
C       NPROB IS A POSITIVE INTEGER VARIABLE WHICH DEFINES THE
C         NUMBER OF THE PROBLEM. NPROB MUST NOT EXCEED 14.
C
C     SUBPROGRAMS CALLED
C
C       FORTRAN-SUPPLIED ... DATAN,DCOS,DEXP,DMIN1,DSIN,DSQRT,
C                            MAX0,MIN0
C
C     MINPACK. VERSION OF AUGUST 1978.
C     BURTON S. GARBOW, KENNETH E. HILLSTROM, JORGE J. MORE
C
C     **********
      INTEGER I,IVAR,J,K,K1,K2,ML,MU
      DOUBLE PRECISION C1,C3,C4,C5,C6,C9,EIGHT,FIFTN,FIVE,FOUR,H,
     1                 HUNDRD,ONE,PROD,SIX,SUM,SUM1,SUM2,TEMP,TEMP1,
     2                 TEMP2,TEMP3,TEMP4,TEN,THREE,TI,TJ,TK,TPI,
     3                 TWENTY,TWO,ZERO
      DOUBLE PRECISION DFLOAT
      DATA ZERO,ONE,TWO,THREE,FOUR,FIVE,SIX,EIGHT,TEN,FIFTN,TWENTY,
     1     HUNDRD
     2     /0.0D0,1.0D0,2.0D0,3.0D0,4.0D0,5.0D0,6.0D0,8.0D0,1.0D1,
     3      1.5D1,2.0D1,1.0D2/
      DATA C1,C3,C4,C5,C6,C9 /1.0D4,2.0D2,2.02D1,1.98D1,1.8D2,2.9D1/
      DFLOAT(IVAR) = IVAR
C
C     JACOBIAN ROUTINE SELECTOR.
C
      GO TO (10,20,50,60,90,100,200,230,290,320,350,380,420,450),
     1      NPROB
C
C     ROSENBROCK FUNCTION.
C
   10 CONTINUE
      FJAC(1,1) = -ONE
      FJAC(1,2) = ZERO
      FJAC(2,1) = -TWENTY*X(1)
      FJAC(2,2) = TEN
      GO TO 490
C
C     POWELL SINGULAR FUNCTION.
C
   20 CONTINUE
      DO 40 K = 1, 4
         DO 30 J = 1, 4
            FJAC(K,J) = ZERO
   30       CONTINUE
   40    CONTINUE
      FJAC(1,1) = ONE
      FJAC(1,2) = TEN
      FJAC(2,3) = DSQRT(FIVE)
      FJAC(2,4) = -FJAC(2,3)
      FJAC(3,2) = TWO*(X(2) - TWO*X(3))
      FJAC(3,3) = -TWO*FJAC(3,2)
      FJAC(4,1) = TWO*DSQRT(TEN)*(X(1) - X(4))
      FJAC(4,4) = -FJAC(4,1)
      GO TO 490
C
C     POWELL BADLY SCALED FUNCTION.
C
   50 CONTINUE
      FJAC(1,1) = C1*X(2)
      FJAC(1,2) = C1*X(1)
      FJAC(2,1) = -DEXP(-X(1))
      FJAC(2,2) = -DEXP(-X(2))
      GO TO 490
C
C     WOOD FUNCTION.
C
   60 CONTINUE
      DO 80 K = 1, 4
         DO 70 J = 1, 4
            FJAC(K,J) = ZERO
   70       CONTINUE
   80    CONTINUE
      TEMP1 = X(2) - THREE*X(1)**2
      TEMP2 = X(4) - THREE*X(3)**2
      FJAC(1,1) = -C3*TEMP1 + ONE
      FJAC(1,2) = -C3*X(1)
      FJAC(2,1) = -TWO*C3*X(1)
      FJAC(2,2) = C3 + C4
      FJAC(2,4) = C5
      FJAC(3,3) = -C6*TEMP2 + ONE
      FJAC(3,4) = -C6*X(3)
      FJAC(4,2) = C5
      FJAC(4,3) = -TWO*C6*X(3)
      FJAC(4,4) = C6 + C4
      GO TO 490
C
C     HELICAL VALLEY FUNCTION.
C
   90 CONTINUE
      TPI = EIGHT*DATAN(ONE)
      TEMP = X(1)**2 + X(2)**2
      TEMP1 = TPI*TEMP
      TEMP2 = DSQRT(TEMP)
      FJAC(1,1) = HUNDRD*X(2)/TEMP1
      FJAC(1,2) = -HUNDRD*X(1)/TEMP1
      FJAC(1,3) = TEN
      FJAC(2,1) = TEN*X(1)/TEMP2
      FJAC(2,2) = TEN*X(2)/TEMP2
      FJAC(2,3) = ZERO
      FJAC(3,1) = ZERO
      FJAC(3,2) = ZERO
      FJAC(3,3) = ONE
      GO TO 490
C
C     WATSON FUNCTION.
C
  100 CONTINUE
      DO 120 K = 1, N
         DO 110 J = K, N
            FJAC(K,J) = ZERO
  110       CONTINUE
  120    CONTINUE
      DO 170 I = 1, 29
         TI = DFLOAT(I)/C9
         SUM1 = ZERO
         TEMP = ONE
         DO 130 J = 2, N
            SUM1 = SUM1 + DFLOAT(J-1)*TEMP*X(J)
            TEMP = TI*TEMP
  130       CONTINUE
         SUM2 = ZERO
         TEMP = ONE
         DO 140 J = 1, N
            SUM2 = SUM2 + TEMP*X(J)
            TEMP = TI*TEMP
  140       CONTINUE
         TEMP1 = TWO*(SUM1 - SUM2**2 - ONE)
         TEMP2 = TWO*SUM2
         TEMP = TI**2
         TK = ONE
         DO 160 K = 1, N
            TJ = TK
            DO 150 J = K, N
               FJAC(K,J) = FJAC(K,J)
     1                     + TJ
     2                       *((DFLOAT(K-1)/TI - TEMP2)
     3                         *(DFLOAT(J-1)/TI - TEMP2) - TEMP1)
               TJ = TI*TJ
  150          CONTINUE
            TK = TEMP*TK
  160       CONTINUE
  170    CONTINUE
      FJAC(1,1) = FJAC(1,1) + SIX*X(1)**2 - TWO*X(2) + THREE
      FJAC(1,2) = FJAC(1,2) - TWO*X(1)
      FJAC(2,2) = FJAC(2,2) + ONE
      DO 190 K = 1, N
         DO 180 J = K, N
            FJAC(J,K) = FJAC(K,J)
  180       CONTINUE
  190    CONTINUE
      GO TO 490
C
C     CHEBYQUAD FUNCTION.
C
  200 CONTINUE
      TK = ONE/DFLOAT(N)
      DO 220 J = 1, N
         TEMP1 = ONE
         TEMP2 = TWO*X(J) - ONE
         TEMP = TWO*TEMP2
         TEMP3 = ZERO
         TEMP4 = TWO
         DO 210 K = 1, N
            FJAC(K,J) = TK*TEMP4
            TI = FOUR*TEMP2 + TEMP*TEMP4 - TEMP3
            TEMP3 = TEMP4
            TEMP4 = TI
            TI = TEMP*TEMP2 - TEMP1
            TEMP1 = TEMP2
            TEMP2 = TI
  210       CONTINUE
  220    CONTINUE
      GO TO 490
C
C     BROWN ALMOST-LINEAR FUNCTION.
C
  230 CONTINUE
      PROD = ONE
      DO 250 J = 1, N
         PROD = X(J)*PROD
         DO 240 K = 1, N
            FJAC(K,J) = ONE
  240       CONTINUE
         FJAC(J,J) = TWO
  250    CONTINUE
      DO 280 J = 1, N
         TEMP = X(J)
         IF (TEMP .NE. ZERO) GO TO 270
         TEMP = ONE
         PROD = ONE
         DO 260 K = 1, N
            IF (K .NE. J) PROD = X(K)*PROD
  260       CONTINUE
  270    CONTINUE
         FJAC(N,J) = PROD/TEMP
  280    CONTINUE
      GO TO 490
C
C     DISCRETE BOUNDARY VALUE FUNCTION.
C
  290 CONTINUE
      H = ONE/DFLOAT(N+1)
      DO 310 K = 1, N
         TEMP = THREE*(X(K) + DFLOAT(K)*H + ONE)**2
         DO 300 J = 1, N
            FJAC(K,J) = ZERO
  300       CONTINUE
         FJAC(K,K) = TWO + TEMP*H**2/TWO
         IF (K .NE. 1) FJAC(K,K-1) = -ONE
         IF (K .NE. N) FJAC(K,K+1) = -ONE
  310    CONTINUE
      GO TO 490
C
C     DISCRETE INTEGRAL EQUATION FUNCTION.
C
  320 CONTINUE
      H = ONE/DFLOAT(N+1)
      DO 340 K = 1, N
         TK = DFLOAT(K)*H
         DO 330 J = 1, N
            TJ = DFLOAT(J)*H
            TEMP = THREE*(X(J) + TJ + ONE)**2
            FJAC(K,J) = H*DMIN1(TJ*(ONE-TK),TK*(ONE-TJ))*TEMP/TWO
  330       CONTINUE
         FJAC(K,K) = FJAC(K,K) + ONE
  340    CONTINUE
      GO TO 490
C
C     TRIGONOMETRIC FUNCTION.
C
  350 CONTINUE
      DO 370 J = 1, N
         TEMP = DSIN(X(J))
         DO 360 K = 1, N
            FJAC(K,J) = TEMP
  360       CONTINUE
         FJAC(J,J) = DFLOAT(J+1)*TEMP - DCOS(X(J))
  370    CONTINUE
      GO TO 490
C
C     VARIABLY DIMENSIONED FUNCTION.
C
  380 CONTINUE
      SUM = ZERO
      DO 390 J = 1, N
         SUM = SUM + DFLOAT(J)*(X(J) - ONE)
  390    CONTINUE
      TEMP = ONE + SIX*SUM**2
      DO 410 K = 1, N
         DO 400 J = K, N
            FJAC(K,J) = DFLOAT(K*J)*TEMP
            FJAC(J,K) = FJAC(K,J)
  400       CONTINUE
         FJAC(K,K) = FJAC(K,K) + ONE
  410    CONTINUE
      GO TO 490
C
C     BROYDEN TRIDIAGONAL FUNCTION.
C
  420 CONTINUE
      DO 440 K = 1, N
         DO 430 J = 1, N
            FJAC(K,J) = ZERO
  430       CONTINUE
         FJAC(K,K) = THREE - FOUR*X(K)
         IF (K .NE. 1) FJAC(K,K-1) = -ONE
         IF (K .NE. N) FJAC(K,K+1) = -TWO
  440    CONTINUE
      GO TO 490
C
C     BROYDEN BANDED FUNCTION.
C
  450 CONTINUE
      ML = 5
      MU = 1
      DO 480 K = 1, N
         DO 460 J = 1, N
            FJAC(K,J) = ZERO
  460       CONTINUE
         K1 = MAX0(1,K-ML)
         K2 = MIN0(K+MU,N)
         DO 470 J = K1, K2
            IF (J .NE. K) FJAC(K,J) = -(ONE + TWO*X(J))
  470       CONTINUE
         FJAC(K,K) = TWO + FIFTN*X(K)**2
  480    CONTINUE
  490 CONTINUE
      RETURN
C
C     LAST CARD OF SUBROUTINE VECJAC.
C
      END
C ===== 3. DOUBLE PRECISION TESTING AIDS FOR NONLINEAR LEAST-SQUARES.
      SUBROUTINE INITPT(N,X,NPROB,FACTOR)                               00000010
      INTEGER N,NPROB
      DOUBLE PRECISION FACTOR
      DOUBLE PRECISION X(N)
C     **********
C
C     SUBROUTINE INITPT
C
C     THIS SUBROUTINE SPECIFIES THE STANDARD STARTING POINTS FOR THE
C     FUNCTIONS DEFINED BY SUBROUTINE SSQFCN. THE SUBROUTINE RETURNS
C     IN X A MULTIPLE (FACTOR) OF THE STANDARD STARTING POINT. FOR
C     THE 11TH FUNCTION THE STANDARD STARTING POINT IS ZERO, SO IN
C     THIS CASE, IF FACTOR IS NOT UNITY, THEN THE SUBROUTINE RETURNS
C     THE VECTOR  X(J) = FACTOR, J=1,...,N.
C
C     THE SUBROUTINE STATEMENT IS
C
C       SUBROUTINE INITPT(N,X,NPROB,FACTOR)
C
C     WHERE
C
C       N IS A POSITIVE INTEGER INPUT VARIABLE.
C
C       X IS AN OUTPUT ARRAY OF LENGTH N WHICH CONTAINS THE STANDARD
C         STARTING POINT FOR PROBLEM NPROB MULTIPLIED BY FACTOR.
C
C       NPROB IS A POSITIVE INTEGER INPUT VARIABLE WHICH DEFINES THE
C         NUMBER OF THE PROBLEM. NPROB MUST NOT EXCEED 18.
C
C       FACTOR IS AN INPUT VARIABLE WHICH SPECIFIES THE MULTIPLE OF
C         THE STANDARD STARTING POINT. IF FACTOR IS UNITY, NO
C         MULTIPLICATION IS PERFORMED.
C
C     MINPACK. VERSION OF JULY 1978.
C     BURTON S. GARBOW, KENNETH E. HILLSTROM, JORGE J. MORE
C
C     **********
      INTEGER IVAR,J
      DOUBLE PRECISION C1,C2,C3,C4,C5,C6,C7,C8,C9,C10,C11,C12,C13,C14,
     1                 C15,C16,C17,FIVE,H,HALF,ONE,SEVEN,TEN,THREE,
     2                 TWENTY,TWNTF,TWO,ZERO
      DOUBLE PRECISION DFLOAT
      DATA ZERO,HALF,ONE,TWO,THREE,FIVE,SEVEN,TEN,TWENTY,TWNTF
     1     /0.0D0,5.0D-1,1.0D0,2.0D0,3.0D0,5.0D0,7.0D0,1.0D1,2.0D1,
     2      2.5D1/
      DATA C1,C2,C3,C4,C5,C6,C7,C8,C9,C10,C11,C12,C13,C14,C15,C16,C17
     1     /1.2D0,2.5D-1,3.9D-1,4.15D-1,2.0D-2,4.0D3,2.5D2,3.0D-1,
     2      4.0D-1,1.5D0,1.0D-2,1.3D0,6.5D-1,7.0D-1,6.0D-1,4.5D0,
     3      5.5D0/
      DFLOAT(IVAR) = IVAR
C
C     SELECTION OF INITIAL POINT.
C
      GO TO (10,10,10,30,40,50,60,70,80,90,100,120,130,140,150,170,
     1       190,200), NPROB
C
C     LINEAR FUNCTION - FULL RANK OR RANK 1.
C
   10 CONTINUE
      DO 20 J = 1, N
         X(J) = ONE
   20    CONTINUE
      GO TO 210
C
C     ROSENBROCK FUNCTION.
C
   30 CONTINUE
      X(1) = -C1
      X(2) = ONE
      GO TO 210
C
C     HELICAL VALLEY FUNCTION.
C
   40 CONTINUE
      X(1) = -ONE
      X(2) = ZERO
      X(3) = ZERO
      GO TO 210
C
C     POWELL SINGULAR FUNCTION.
C
   50 CONTINUE
      X(1) = THREE
      X(2) = -ONE
      X(3) = ZERO
      X(4) = ONE
      GO TO 210
C
C     FREUDENSTEIN AND ROTH FUNCTION.
C
   60 CONTINUE
      X(1) = HALF
      X(2) = -TWO
      GO TO 210
C
C     BARD FUNCTION.
C
   70 CONTINUE
      X(1) = ONE
      X(2) = ONE
      X(3) = ONE
      GO TO 210
C
C     KOWALIK AND OSBORNE FUNCTION.
C
   80 CONTINUE
      X(1) = C2
      X(2) = C3
      X(3) = C4
      X(4) = C3
      GO TO 210
C
C     MEYER FUNCTION.
C
   90 CONTINUE
      X(1) = C5
      X(2) = C6
      X(3) = C7
      GO TO 210
C
C     WATSON FUNCTION.
C
  100 CONTINUE
      DO 110 J = 1, N
         X(J) = ZERO
  110    CONTINUE
      GO TO 210
C
C     BOX 3-DIMENSIONAL FUNCTION.
C
  120 CONTINUE
      X(1) = ZERO
      X(2) = TEN
      X(3) = TWENTY
      GO TO 210
C
C     JENNRICH AND SAMPSON FUNCTION.
C
  130 CONTINUE
      X(1) = C8
      X(2) = C9
      GO TO 210
C
C     BROWN AND DENNIS FUNCTION.
C
  140 CONTINUE
      X(1) = TWNTF
      X(2) = FIVE
      X(3) = -FIVE
      X(4) = -ONE
      GO TO 210
C
C     CHEBYQUAD FUNCTION.
C
  150 CONTINUE
      H = ONE/DFLOAT(N+1)
      DO 160 J = 1, N
         X(J) = DFLOAT(J)*H
  160    CONTINUE
      GO TO 210
C
C     BROWN ALMOST-LINEAR FUNCTION.
C
  170 CONTINUE
      DO 180 J = 1, N
         X(J) = HALF
  180    CONTINUE
      GO TO 210
C
C     OSBORNE 1 FUNCTION.
C
  190 CONTINUE
      X(1) = HALF
      X(2) = C10
      X(3) = -ONE
      X(4) = C11
      X(5) = C5
      GO TO 210
C
C     OSBORNE 2 FUNCTION.
C
  200 CONTINUE
      X(1) = C12
      X(2) = C13
      X(3) = C13
      X(4) = C14
      X(5) = C15
      X(6) = THREE
      X(7) = FIVE
      X(8) = SEVEN
      X(9) = TWO
      X(10) = C16
      X(11) = C17
  210 CONTINUE
C
C     COMPUTE MULTIPLE OF INITIAL POINT.
C
      IF (FACTOR .EQ. ONE) GO TO 260
      IF (NPROB .EQ. 11) GO TO 230
         DO 220 J = 1, N
            X(J) = FACTOR*X(J)
  220       CONTINUE
         GO TO 250
  230 CONTINUE
         DO 240 J = 1, N
            X(J) = FACTOR
  240       CONTINUE
  250 CONTINUE
  260 CONTINUE
      RETURN
C
C     LAST CARD OF SUBROUTINE INITPT.
C
      END
      SUBROUTINE SSQFCN(M,N,X,FVEC,NPROB)                               00000010
      INTEGER M,N,NPROB
      DOUBLE PRECISION X(N),FVEC(M)
C     **********
C
C     SUBROUTINE SSQFCN
C
C     THIS SUBROUTINE DEFINES THE FUNCTIONS OF EIGHTEEN NONLINEAR
C     LEAST SQUARES PROBLEMS. THE ALLOWABLE VALUES OF (M,N) FOR
C     FUNCTIONS 1,2 AND 3 ARE VARIABLE BUT WITH M .GE. N.
C     FOR FUNCTIONS 4,5,6,7,8,9 AND 10 THE VALUES OF (M,N) ARE
C     (2,2),(3,3),(4,4),(2,2),(15,3),(11,4) AND (16,3), RESPECTIVELY.
C     FUNCTION 11 (WATSON) HAS M = 31 WITH N USUALLY 6 OR 9.
C     HOWEVER, ANY N, N = 2,...,31, IS PERMITTED.
C     FUNCTIONS 12,13 AND 14 HAVE N = 3,2 AND 4, RESPECTIVELY, BUT
C     ALLOW ANY M .GE. N, WITH THE USUAL CHOICES BEING 10,10 AND 20.
C     FUNCTION 15 (CHEBYQUAD) ALLOWS M AND N VARIABLE WITH M .GE. N.
C     FUNCTION 16 (BROWN) ALLOWS N VARIABLE WITH M = N.
C     FOR FUNCTIONS 17 AND 18, THE VALUES OF (M,N) ARE
C     (33,5) AND (65,11), RESPECTIVELY.
C
C     THE SUBROUTINE STATEMENT IS
C
C       SUBROUTINE SSQFCN(M,N,X,FVEC,NPROB)
C
C     WHERE
C
C       M AND N ARE POSITIVE INTEGER INPUT VARIABLES. N MUST NOT
C         EXCEED M.
C
C       X IS AN INPUT ARRAY OF LENGTH N.
C
C       FVEC IS AN OUTPUT ARRAY OF LENGTH M WHICH CONTAINS THE NPROB
C         FUNCTION EVALUATED AT X.
C
C       NPROB IS A POSITIVE INTEGER INPUT VARIABLE WHICH DEFINES THE
C         NUMBER OF THE PROBLEM. NPROB MUST NOT EXCEED 18.
C
C     SUBPROGRAMS CALLED
C
C       FORTRAN-SUPPLIED ... DATAN,DCOS,DEXP,DSIN,DSQRT,DSIGN
C
C     MINPACK. VERSION OF JULY 1978.
C     BURTON S. GARBOW, KENNETH E. HILLSTROM, JORGE J. MORE
C
C     **********
      INTEGER I,IEV,IVAR,J,NM1
      DOUBLE PRECISION C13,C14,C29,C45,DIV,DX,EIGHT,FIVE,ONE,PROD,SUM,
     1                 S1,S2,TEMP,TEN,TI,TMP1,TMP2,TMP3,TMP4,TPI,TWO,
     2                 ZERO,ZP25,ZP5
      DOUBLE PRECISION V(11),Y1(15),Y2(11),Y3(16),Y4(33),Y5(65)
      DOUBLE PRECISION DFLOAT
      DATA ZERO,ZP25,ZP5,ONE,TWO,FIVE,EIGHT,TEN,C13,C14,C29,C45
     1     /0.0D0,2.5D-1,5.0D-1,1.0D0,2.0D0,5.0D0,8.0D0,1.0D1,1.3D1,
     2      1.4D1,2.9D1,4.5D1/
      DATA V(1),V(2),V(3),V(4),V(5),V(6),V(7),V(8),V(9),V(10),V(11)
     1     /4.0D0,2.0D0,1.0D0,5.0D-1,2.5D-1,1.67D-1,1.25D-1,1.0D-1,
     2      8.33D-2,7.14D-2,6.25D-2/
      DATA Y1(1),Y1(2),Y1(3),Y1(4),Y1(5),Y1(6),Y1(7),Y1(8),Y1(9),
     1     Y1(10),Y1(11),Y1(12),Y1(13),Y1(14),Y1(15)
     2     /1.4D-1,1.8D-1,2.2D-1,2.5D-1,2.9D-1,3.2D-1,3.5D-1,3.9D-1,
     3      3.7D-1,5.8D-1,7.3D-1,9.6D-1,1.34D0,2.1D0,4.39D0/
      DATA Y2(1),Y2(2),Y2(3),Y2(4),Y2(5),Y2(6),Y2(7),Y2(8),Y2(9),
     1     Y2(10),Y2(11)
     2     /1.957D-1,1.947D-1,1.735D-1,1.6D-1,8.44D-2,6.27D-2,4.56D-2,
     3      3.42D-2,3.23D-2,2.35D-2,2.46D-2/
      DATA Y3(1),Y3(2),Y3(3),Y3(4),Y3(5),Y3(6),Y3(7),Y3(8),Y3(9),
     1     Y3(10),Y3(11),Y3(12),Y3(13),Y3(14),Y3(15),Y3(16)
     2     /3.478D4,2.861D4,2.365D4,1.963D4,1.637D4,1.372D4,1.154D4,
     3      9.744D3,8.261D3,7.03D3,6.005D3,5.147D3,4.427D3,3.82D3,
     4      3.307D3,2.872D3/
      DATA Y4(1),Y4(2),Y4(3),Y4(4),Y4(5),Y4(6),Y4(7),Y4(8),Y4(9),
     1     Y4(10),Y4(11),Y4(12),Y4(13),Y4(14),Y4(15),Y4(16),Y4(17),
     2     Y4(18),Y4(19),Y4(20),Y4(21),Y4(22),Y4(23),Y4(24),Y4(25),
     3     Y4(26),Y4(27),Y4(28),Y4(29),Y4(30),Y4(31),Y4(32),Y4(33)
     4     /8.44D-1,9.08D-1,9.32D-1,9.36D-1,9.25D-1,9.08D-1,8.81D-1,
     5      8.5D-1,8.18D-1,7.84D-1,7.51D-1,7.18D-1,6.85D-1,6.58D-1,
     6      6.28D-1,6.03D-1,5.8D-1,5.58D-1,5.38D-1,5.22D-1,5.06D-1,
     7      4.9D-1,4.78D-1,4.67D-1,4.57D-1,4.48D-1,4.38D-1,4.31D-1,
     8      4.24D-1,4.2D-1,4.14D-1,4.11D-1,4.06D-1/
      DATA Y5(1),Y5(2),Y5(3),Y5(4),Y5(5),Y5(6),Y5(7),Y5(8),Y5(9),
     1     Y5(10),Y5(11),Y5(12),Y5(13),Y5(14),Y5(15),Y5(16),Y5(17),
     2     Y5(18),Y5(19),Y5(20),Y5(21),Y5(22),Y5(23),Y5(24),Y5(25),
     3     Y5(26),Y5(27),Y5(28),Y5(29),Y5(30),Y5(31),Y5(32),Y5(33),
     4     Y5(34),Y5(35),Y5(36),Y5(37),Y5(38),Y5(39),Y5(40),Y5(41),
     5     Y5(42),Y5(43),Y5(44),Y5(45),Y5(46),Y5(47),Y5(48),Y5(49),
     6     Y5(50),Y5(51),Y5(52),Y5(53),Y5(54),Y5(55),Y5(56),Y5(57),
     7     Y5(58),Y5(59),Y5(60),Y5(61),Y5(62),Y5(63),Y5(64),Y5(65)
     8     /1.366D0,1.191D0,1.112D0,1.013D0,9.91D-1,8.85D-1,8.31D-1,
     9      8.47D-1,7.86D-1,7.25D-1,7.46D-1,6.79D-1,6.08D-1,6.55D-1,
     A      6.16D-1,6.06D-1,6.02D-1,6.26D-1,6.51D-1,7.24D-1,6.49D-1,
     B      6.49D-1,6.94D-1,6.44D-1,6.24D-1,6.61D-1,6.12D-1,5.58D-1,
     C      5.33D-1,4.95D-1,5.0D-1,4.23D-1,3.95D-1,3.75D-1,3.72D-1,
     D      3.91D-1,3.96D-1,4.05D-1,4.28D-1,4.29D-1,5.23D-1,5.62D-1,
     E      6.07D-1,6.53D-1,6.72D-1,7.08D-1,6.33D-1,6.68D-1,6.45D-1,
     F      6.32D-1,5.91D-1,5.59D-1,5.97D-1,6.25D-1,7.39D-1,7.1D-1,
     G      7.29D-1,7.2D-1,6.36D-1,5.81D-1,4.28D-1,2.92D-1,1.62D-1,
     H      9.8D-2,5.4D-2/
      DFLOAT(IVAR) = IVAR
C
C     FUNCTION ROUTINE SELECTOR.
C
      GO TO (10,40,70,110,120,130,140,150,170,190,210,250,270,290,310,
     1       360,390,410), NPROB
C
C     LINEAR FUNCTION - FULL RANK.
C
   10 CONTINUE
      SUM = ZERO
      DO 20 J = 1, N
         SUM = SUM + X(J)
   20    CONTINUE
      TEMP = TWO*SUM/DFLOAT(M) + ONE
      DO 30 I = 1, M
         FVEC(I) = -TEMP
         IF (I .LE. N) FVEC(I) = FVEC(I) + X(I)
   30    CONTINUE
      GO TO 430
C
C     LINEAR FUNCTION - RANK 1.
C
   40 CONTINUE
      SUM = ZERO
      DO 50 J = 1, N
         SUM = SUM + DFLOAT(J)*X(J)
   50    CONTINUE
      DO 60 I = 1, M
         FVEC(I) = DFLOAT(I)*SUM - ONE
   60    CONTINUE
      GO TO 430
C
C     LINEAR FUNCTION - RANK 1 WITH ZERO COLUMNS AND ROWS.
C
   70 CONTINUE
      SUM = ZERO
      NM1 = N - 1
      IF (NM1 .LT. 2) GO TO 90
      DO 80 J = 2, NM1
         SUM = SUM + DFLOAT(J)*X(J)
   80    CONTINUE
   90 CONTINUE
      DO 100 I = 1, M
         FVEC(I) = DFLOAT(I-1)*SUM - ONE
  100    CONTINUE
      FVEC(M) = -ONE
      GO TO 430
C
C     ROSENBROCK FUNCTION.
C
  110 CONTINUE
      FVEC(1) = TEN*(X(2) - X(1)**2)
      FVEC(2) = ONE - X(1)
      GO TO 430
C
C     HELICAL VALLEY FUNCTION.
C
  120 CONTINUE
      TPI = EIGHT*DATAN(ONE)
      TMP1 = DSIGN(ZP25,X(2))
      IF (X(1) .GT. ZERO) TMP1 = DATAN(X(2)/X(1))/TPI
      IF (X(1) .LT. ZERO) TMP1 = DATAN(X(2)/X(1))/TPI + ZP5
      TMP2 = DSQRT(X(1)**2+X(2)**2)
      FVEC(1) = TEN*(X(3) - TEN*TMP1)
      FVEC(2) = TEN*(TMP2 - ONE)
      FVEC(3) = X(3)
      GO TO 430
C
C     POWELL SINGULAR FUNCTION.
C
  130 CONTINUE
      FVEC(1) = X(1) + TEN*X(2)
      FVEC(2) = DSQRT(FIVE)*(X(3) - X(4))
      FVEC(3) = (X(2) - TWO*X(3))**2
      FVEC(4) = DSQRT(TEN)*(X(1) - X(4))**2
      GO TO 430
C
C     FREUDENSTEIN AND ROTH FUNCTION.
C
  140 CONTINUE
      FVEC(1) = -C13 + X(1) + ((FIVE - X(2))*X(2) - TWO)*X(2)
      FVEC(2) = -C29 + X(1) + ((ONE + X(2))*X(2) - C14)*X(2)
      GO TO 430
C
C     BARD FUNCTION.
C
  150 CONTINUE
      DO 160 I = 1, 15
         TMP1 = DFLOAT(I)
         TMP2 = DFLOAT(16-I)
         TMP3 = TMP1
         IF (I .GT. 8) TMP3 = TMP2
         FVEC(I) = Y1(I) - (X(1) + TMP1/(X(2)*TMP2 + X(3)*TMP3))
  160    CONTINUE
      GO TO 430
C
C     KOWALIK AND OSBORNE FUNCTION.
C
  170 CONTINUE
      DO 180 I = 1, 11
         TMP1 = V(I)*(V(I) + X(2))
         TMP2 = V(I)*(V(I) + X(3)) + X(4)
         FVEC(I) = Y2(I) - X(1)*TMP1/TMP2
  180    CONTINUE
      GO TO 430
C
C     MEYER FUNCTION.
C
  190 CONTINUE
      DO 200 I = 1, 16
         TEMP = FIVE*DFLOAT(I) + C45 + X(3)
         TMP1 = X(2)/TEMP
         TMP2 = DEXP(TMP1)
         FVEC(I) = X(1)*TMP2 - Y3(I)
  200    CONTINUE
      GO TO 430
C
C     WATSON FUNCTION.
C
  210 CONTINUE
      DO 240 I = 1, 29
         DIV = DFLOAT(I)/C29
         S1 = ZERO
         DX = ONE
         DO 220 J = 2, N
            S1 = S1 + DFLOAT(J-1)*DX*X(J)
            DX = DIV*DX
  220       CONTINUE
         S2 = ZERO
         DX = ONE
         DO 230 J = 1, N
            S2 = S2 + DX*X(J)
            DX = DIV*DX
  230       CONTINUE
         FVEC(I) = S1 - S2**2 - ONE
  240    CONTINUE
      FVEC(30) = X(1)
      FVEC(31) = X(2) - X(1)**2 - ONE
      GO TO 430
C
C     BOX 3-DIMENSIONAL FUNCTION.
C
  250 CONTINUE
      DO 260 I = 1, M
         TEMP = DFLOAT(I)
         TMP1 = TEMP/TEN
         FVEC(I) = DEXP(-TMP1*X(1)) - DEXP(-TMP1*X(2))
     1             + (DEXP(-TEMP) - DEXP(-TMP1))*X(3)
  260    CONTINUE
      GO TO 430
C
C     JENNRICH AND SAMPSON FUNCTION.
C
  270 CONTINUE
      DO 280 I = 1, M
         TEMP = DFLOAT(I)
         FVEC(I) = TWO + TWO*TEMP - DEXP(TEMP*X(1)) - DEXP(TEMP*X(2))
  280    CONTINUE
      GO TO 430
C
C     BROWN AND DENNIS FUNCTION.
C
  290 CONTINUE
      DO 300 I = 1, M
         TEMP = DFLOAT(I)/FIVE
         TMP1 = X(1) + TEMP*X(2) - DEXP(TEMP)
         TMP2 = X(3) + DSIN(TEMP)*X(4) - DCOS(TEMP)
         FVEC(I) = TMP1**2 + TMP2**2
  300    CONTINUE
      GO TO 430
C
C     CHEBYQUAD FUNCTION.
C
  310 CONTINUE
      DO 320 I = 1, M
         FVEC(I) = ZERO
  320    CONTINUE
      DO 340 J = 1, N
         TMP1 = ONE
         TMP2 = TWO*X(J) - ONE
         TEMP = TWO*TMP2
         DO 330 I = 1, M
            FVEC(I) = FVEC(I) + TMP2
            TI = TEMP*TMP2 - TMP1
            TMP1 = TMP2
            TMP2 = TI
  330       CONTINUE
  340    CONTINUE
      DX = ONE/DFLOAT(N)
      IEV = -1
      DO 350 I = 1, M
         FVEC(I) = DX*FVEC(I)
         IF (IEV .GT. 0) FVEC(I) = FVEC(I) + ONE/(DFLOAT(I)**2 - ONE)
         IEV = -IEV
  350    CONTINUE
      GO TO 430
C
C     BROWN ALMOST-LINEAR FUNCTION.
C
  360 CONTINUE
      SUM = -DFLOAT(N+1)
      PROD = ONE
      DO 370 J = 1, N
         SUM = SUM + X(J)
         PROD = X(J)*PROD
  370    CONTINUE
      DO 380 I = 1, N
         FVEC(I) = X(I) + SUM
  380    CONTINUE
      FVEC(N) = PROD - ONE
      GO TO 430
C
C     OSBORNE 1 FUNCTION.
C
  390 CONTINUE
      DO 400 I = 1, 33
         TEMP = TEN*DFLOAT(I-1)
         TMP1 = DEXP(-X(4)*TEMP)
         TMP2 = DEXP(-X(5)*TEMP)
         FVEC(I) = Y4(I) - (X(1) + X(2)*TMP1 + X(3)*TMP2)
  400    CONTINUE
      GO TO 430
C
C     OSBORNE 2 FUNCTION.
C
  410 CONTINUE
      DO 420 I = 1, 65
         TEMP = DFLOAT(I-1)/TEN
         TMP1 = DEXP(-X(5)*TEMP)
         TMP2 = DEXP(-X(6)*(TEMP-X(9))**2)
         TMP3 = DEXP(-X(7)*(TEMP-X(10))**2)
         TMP4 = DEXP(-X(8)*(TEMP-X(11))**2)
         FVEC(I) = Y5(I)
     1             - (X(1)*TMP1 + X(2)*TMP2 + X(3)*TMP3 + X(4)*TMP4)
  420    CONTINUE
  430 CONTINUE
      RETURN
C
C     LAST CARD OF SUBROUTINE SSQFCN.
C
      END
      SUBROUTINE SSQJAC(M,N,X,FJAC,LDFJAC,NPROB)                        00000010
      INTEGER M,N,LDFJAC,NPROB
      DOUBLE PRECISION X(N),FJAC(LDFJAC,N)
C     **********
C
C     SUBROUTINE SSQJAC
C
C     THIS SUBROUTINE DEFINES THE JACOBIAN MATRICES OF EIGHTEEN
C     NONLINEAR LEAST SQUARES PROBLEMS. THE PROBLEM DIMENSIONS ARE
C     AS DESCRIBED IN THE PROLOGUE COMMENTS OF SSQFCN.
C
C     THE SUBROUTINE STATEMENT IS
C
C       SUBROUTINE SSQJAC(M,N,X,FJAC,LDFJAC,NPROB)
C
C     WHERE
C
C       M AND N ARE POSITIVE INTEGER INPUT VARIABLES. N MUST NOT
C         EXCEED M.
C
C       X IS AN INPUT ARRAY OF LENGTH N.
C
C       FJAC IS AN M BY N OUTPUT ARRAY WHICH CONTAINS THE JACOBIAN
C         MATRIX OF THE NPROB FUNCTION EVALUATED AT X.
C
C       LDFJAC IS A POSITIVE INTEGER INPUT VARIABLE NOT LESS THAN M
C         WHICH SPECIFIES THE LEADING DIMENSION OF THE ARRAY FJAC.
C
C       NPROB IS A POSITIVE INTEGER VARIABLE WHICH DEFINES THE
C         NUMBER OF THE PROBLEM. NPROB MUST NOT EXCEED 18.
C
C     SUBPROGRAMS CALLED
C
C       FORTRAN-SUPPLIED ... DATAN,DCOS,DEXP,DSIN,DSQRT
C
C     MINPACK. VERSION OF JULY 1978.
C     BURTON S. GARBOW, KENNETH E. HILLSTROM, JORGE J. MORE
C
C     **********
      INTEGER I,IVAR,J,K,MM1,NM1
      DOUBLE PRECISION C14,C20,C29,C45,C100,DIV,DX,EIGHT,FIVE,FOUR,
     1                 ONE,PROD,S2,TEMP,TEN,THREE,TI,TMP1,TMP2,TMP3,
     2                 TMP4,TPI,TWO,ZERO
      DOUBLE PRECISION V(11)
      DOUBLE PRECISION DFLOAT
      DATA ZERO,ONE,TWO,THREE,FOUR,FIVE,EIGHT,TEN,C14,C20,C29,C45,C100
     1     /0.0D0,1.0D0,2.0D0,3.0D0,4.0D0,5.0D0,8.0D0,1.0D1,1.4D1,
     2      2.0D1,2.9D1,4.5D1,1.0D2/
      DATA V(1),V(2),V(3),V(4),V(5),V(6),V(7),V(8),V(9),V(10),V(11)
     1     /4.0D0,2.0D0,1.0D0,5.0D-1,2.5D-1,1.67D-1,1.25D-1,1.0D-1,
     2      8.33D-2,7.14D-2,6.25D-2/
      DFLOAT(IVAR) = IVAR
C
C     JACOBIAN ROUTINE SELECTOR.
C
      GO TO (10,40,70,130,140,150,180,190,210,230,250,310,330,350,370,
     1       400,460,480), NPROB
C
C     LINEAR FUNCTION - FULL RANK.
C
   10 CONTINUE
      TEMP = TWO/DFLOAT(M)
      DO 30 J = 1, N
         DO 20 I = 1, M
            FJAC(I,J) = -TEMP
   20       CONTINUE
         FJAC(J,J) = FJAC(J,J) + ONE
   30    CONTINUE
      GO TO 500
C
C     LINEAR FUNCTION - RANK 1.
C
   40 CONTINUE
      DO 60 J = 1, N
         DO 50 I = 1, M
            FJAC(I,J) = DFLOAT(I)*DFLOAT(J)
   50       CONTINUE
   60    CONTINUE
      GO TO 500
C
C     LINEAR FUNCTION - RANK 1 WITH ZERO COLUMNS AND ROWS.
C
   70 CONTINUE
      DO 90 J = 1, N
         DO 80 I = 1, M
            FJAC(I,J) = ZERO
   80       CONTINUE
   90    CONTINUE
      NM1 = N - 1
      MM1 = M - 1
      IF (NM1 .LT. 2) GO TO 120
      DO 110 J = 2, NM1
         DO 100 I = 2, MM1
            FJAC(I,J) = DFLOAT(I-1)*DFLOAT(J)
  100       CONTINUE
  110    CONTINUE
  120 CONTINUE
      GO TO 500
C
C     ROSENBROCK FUNCTION.
C
  130 CONTINUE
      FJAC(1,1) = -C20*X(1)
      FJAC(1,2) = TEN
      FJAC(2,1) = -ONE
      FJAC(2,2) = ZERO
      GO TO 500
C
C     HELICAL VALLEY FUNCTION.
C
  140 CONTINUE
      TPI = EIGHT*DATAN(ONE)
      TEMP = X(1)**2 + X(2)**2
      TMP1 = TPI*TEMP
      TMP2 = DSQRT(TEMP)
      FJAC(1,1) = C100*X(2)/TMP1
      FJAC(1,2) = -C100*X(1)/TMP1
      FJAC(1,3) = TEN
      FJAC(2,1) = TEN*X(1)/TMP2
      FJAC(2,2) = TEN*X(2)/TMP2
      FJAC(2,3) = ZERO
      FJAC(3,1) = ZERO
      FJAC(3,2) = ZERO
      FJAC(3,3) = ONE
      GO TO 500
C
C     POWELL SINGULAR FUNCTION.
C
  150 CONTINUE
      DO 170 J = 1, 4
         DO 160 I = 1, 4
            FJAC(I,J) = ZERO
  160       CONTINUE
  170    CONTINUE
      FJAC(1,1) = ONE
      FJAC(1,2) = TEN
      FJAC(2,3) = DSQRT(FIVE)
      FJAC(2,4) = -FJAC(2,3)
      FJAC(3,2) = TWO*(X(2) - TWO*X(3))
      FJAC(3,3) = -TWO*FJAC(3,2)
      FJAC(4,1) = TWO*DSQRT(TEN)*(X(1) - X(4))
      FJAC(4,4) = -FJAC(4,1)
      GO TO 500
C
C     FREUDENSTEIN AND ROTH FUNCTION.
C
  180 CONTINUE
      FJAC(1,1) = ONE
      FJAC(1,2) = X(2)*(TEN - THREE*X(2)) - TWO
      FJAC(2,1) = ONE
      FJAC(2,2) = X(2)*(TWO + THREE*X(2)) - C14
      GO TO 500
C
C     BARD FUNCTION.
C
  190 CONTINUE
      DO 200 I = 1, 15
         TMP1 = DFLOAT(I)
         TMP2 = DFLOAT(16-I)
         TMP3 = TMP1
         IF (I .GT. 8) TMP3 = TMP2
         TMP4 = (X(2)*TMP2 + X(3)*TMP3)**2
         FJAC(I,1) = -ONE
         FJAC(I,2) = TMP1*TMP2/TMP4
         FJAC(I,3) = TMP1*TMP3/TMP4
  200    CONTINUE
      GO TO 500
C
C     KOWALIK AND OSBORNE FUNCTION.
C
  210 CONTINUE
      DO 220 I = 1, 11
         TMP1 = V(I)*(V(I) + X(2))
         TMP2 = V(I)*(V(I) + X(3)) + X(4)
         FJAC(I,1) = -TMP1/TMP2
         FJAC(I,2) = -V(I)*X(1)/TMP2
         FJAC(I,3) = FJAC(I,1)*FJAC(I,2)
         FJAC(I,4) = FJAC(I,3)/V(I)
  220    CONTINUE
      GO TO 500
C
C     MEYER FUNCTION.
C
  230 CONTINUE
      DO 240 I = 1, 16
         TEMP = FIVE*DFLOAT(I) + C45 + X(3)
         TMP1 = X(2)/TEMP
         TMP2 = DEXP(TMP1)
         FJAC(I,1) = TMP2
         FJAC(I,2) = X(1)*TMP2/TEMP
         FJAC(I,3) = -TMP1*FJAC(I,2)
  240    CONTINUE
      GO TO 500
C
C     WATSON FUNCTION.
C
  250 CONTINUE
      DO 280 I = 1, 29
         DIV = DFLOAT(I)/C29
         S2 = ZERO
         DX = ONE
         DO 260 J = 1, N
            S2 = S2 + DX*X(J)
            DX = DIV*DX
  260       CONTINUE
         TEMP = TWO*DIV*S2
         DX = ONE/DIV
         DO 270 J = 1, N
            FJAC(I,J) = DX*(DFLOAT(J-1) - TEMP)
            DX = DIV*DX
  270       CONTINUE
  280    CONTINUE
      DO 300 J = 1, N
         DO 290 I = 30, 31
            FJAC(I,J) = ZERO
  290       CONTINUE
  300    CONTINUE
      FJAC(30,1) = ONE
      FJAC(31,1) = -TWO*X(1)
      FJAC(31,2) = ONE
      GO TO 500
C
C     BOX 3-DIMENSIONAL FUNCTION.
C
  310 CONTINUE
      DO 320 I = 1, M
         TEMP = DFLOAT(I)
         TMP1 = TEMP/TEN
         FJAC(I,1) = -TMP1*DEXP(-TMP1*X(1))
         FJAC(I,2) = TMP1*DEXP(-TMP1*X(2))
         FJAC(I,3) = DEXP(-TEMP) - DEXP(-TMP1)
  320    CONTINUE
      GO TO 500
C
C     JENNRICH AND SAMPSON FUNCTION.
C
  330 CONTINUE
      DO 340 I = 1, M
         TEMP = DFLOAT(I)
         FJAC(I,1) = -TEMP*DEXP(TEMP*X(1))
         FJAC(I,2) = -TEMP*DEXP(TEMP*X(2))
  340    CONTINUE
      GO TO 500
C
C     BROWN AND DENNIS FUNCTION.
C
  350 CONTINUE
      DO 360 I = 1, M
         TEMP = DFLOAT(I)/FIVE
         TI = DSIN(TEMP)
         TMP1 = X(1) + TEMP*X(2) - DEXP(TEMP)
         TMP2 = X(3) + TI*X(4) - DCOS(TEMP)
         FJAC(I,1) = TWO*TMP1
         FJAC(I,2) = TEMP*FJAC(I,1)
         FJAC(I,3) = TWO*TMP2
         FJAC(I,4) = TI*FJAC(I,3)
  360    CONTINUE
      GO TO 500
C
C     CHEBYQUAD FUNCTION.
C
  370 CONTINUE
      DX = ONE/DFLOAT(N)
      DO 390 J = 1, N
         TMP1 = ONE
         TMP2 = TWO*X(J) - ONE
         TEMP = TWO*TMP2
         TMP3 = ZERO
         TMP4 = TWO
         DO 380 I = 1, M
            FJAC(I,J) = DX*TMP4
            TI = FOUR*TMP2 + TEMP*TMP4 - TMP3
            TMP3 = TMP4
            TMP4 = TI
            TI = TEMP*TMP2 - TMP1
            TMP1 = TMP2
            TMP2 = TI
  380       CONTINUE
  390    CONTINUE
      GO TO 500
C
C     BROWN ALMOST-LINEAR FUNCTION.
C
  400 CONTINUE
      PROD = ONE
      DO 420 J = 1, N
         PROD = X(J)*PROD
         DO 410 I = 1, N
            FJAC(I,J) = ONE
  410       CONTINUE
         FJAC(J,J) = TWO
  420    CONTINUE
      DO 450 J = 1, N
         TEMP = X(J)
         IF (TEMP .NE. ZERO) GO TO 440
         TEMP = ONE
         PROD = ONE
         DO 430 K = 1, N
            IF (K .NE. J) PROD = X(K)*PROD
  430       CONTINUE
  440    CONTINUE
         FJAC(N,J) = PROD/TEMP
  450    CONTINUE
      GO TO 500
C
C     OSBORNE 1 FUNCTION.
C
  460 CONTINUE
      DO 470 I = 1, 33
         TEMP = TEN*DFLOAT(I-1)
         TMP1 = DEXP(-X(4)*TEMP)
         TMP2 = DEXP(-X(5)*TEMP)
         FJAC(I,1) = -ONE
         FJAC(I,2) = -TMP1
         FJAC(I,3) = -TMP2
         FJAC(I,4) = TEMP*X(2)*TMP1
         FJAC(I,5) = TEMP*X(3)*TMP2
  470    CONTINUE
      GO TO 500
C
C     OSBORNE 2 FUNCTION.
C
  480 CONTINUE
      DO 490 I = 1, 65
         TEMP = DFLOAT(I-1)/TEN
         TMP1 = DEXP(-X(5)*TEMP)
         TMP2 = DEXP(-X(6)*(TEMP-X(9))**2)
         TMP3 = DEXP(-X(7)*(TEMP-X(10))**2)
         TMP4 = DEXP(-X(8)*(TEMP-X(11))**2)
         FJAC(I,1) = -TMP1
         FJAC(I,2) = -TMP2
         FJAC(I,3) = -TMP3
         FJAC(I,4) = -TMP4
         FJAC(I,5) = TEMP*X(1)*TMP1
         FJAC(I,6) = X(2)*(TEMP - X(9))**2*TMP2
         FJAC(I,7) = X(3)*(TEMP - X(10))**2*TMP3
         FJAC(I,8) = X(4)*(TEMP - X(11))**2*TMP4
         FJAC(I,9) = -TWO*X(2)*X(6)*(TEMP - X(9))*TMP2
         FJAC(I,10) = -TWO*X(3)*X(7)*(TEMP - X(10))*TMP3
         FJAC(I,11) = -TWO*X(4)*X(8)*(TEMP - X(11))*TMP4
  490    CONTINUE
  500 CONTINUE
      RETURN
C
C     LAST CARD OF SUBROUTINE SSQJAC.
C
      END
C ===== 4. DOUBLE PRECISION TESTING AIDS FOR UNCONSTRAINED NONLINEAR
C =====     OPTIMIZATION.
      SUBROUTINE INITPT(N,X,NPROB,FACTOR)                               00000010
      INTEGER N,NPROB
      DOUBLE PRECISION FACTOR
      DOUBLE PRECISION X(N)
C     **********
C
C     SUBROUTINE INITPT
C
C     THIS SUBROUTINE SPECIFIES THE STANDARD STARTING POINTS FOR THE
C     FUNCTIONS DEFINED BY SUBROUTINE OBJFCN. THE SUBROUTINE RETURNS
C     IN X A MULTIPLE (FACTOR) OF THE STANDARD STARTING POINT. FOR
C     THE SEVENTH FUNCTION THE STANDARD STARTING POINT IS ZERO, SO IN
C     THIS CASE, IF FACTOR IS NOT UNITY, THEN THE SUBROUTINE RETURNS
C     THE VECTOR  X(J) = FACTOR, J=1,...,N.
C
C     THE SUBROUTINE STATEMENT IS
C
C       SUBROUTINE INITPT(N,X,NPROB,FACTOR)
C
C     WHERE
C
C       N IS A POSITIVE INTEGER INPUT VARIABLE.
C
C       X IS AN OUTPUT ARRAY OF LENGTH N WHICH CONTAINS THE STANDARD
C         STARTING POINT FOR PROBLEM NPROB MULTIPLIED BY FACTOR.
C
C       NPROB IS A POSITIVE INTEGER INPUT VARIABLE WHICH DEFINES THE
C         NUMBER OF THE PROBLEM. NPROB MUST NOT EXCEED 18.
C
C       FACTOR IS AN INPUT VARIABLE WHICH SPECIFIES THE MULTIPLE OF
C         THE STANDARD STARTING POINT. IF FACTOR IS UNITY, NO
C         MULTIPLICATION IS PERFORMED.
C
C     MINPACK. VERSION OF JULY 1978.
C     BURTON S. GARBOW, KENNETH E. HILLSTROM, JORGE J. MORE
C
C     **********
      INTEGER IVAR,J
      DOUBLE PRECISION C1,C2,C3,C4,FIVE,H,HALF,ONE,TEN,THREE,TWENTY,
     1                 TWNTF,TWO,ZERO
      DOUBLE PRECISION DFLOAT
      DATA ZERO,HALF,ONE,TWO,THREE,FIVE,TEN,TWENTY,TWNTF
     1     /0.0D0,0.5D0,1.0D0,2.0D0,3.0D0,5.0D0,1.0D1,2.0D1,2.5D1/
      DATA C1,C2,C3,C4 /4.0D-1,2.5D0,1.5D-1,1.2D0/
      DFLOAT(IVAR) = IVAR
C
C     SELECTION OF INITIAL POINT.
C
      GO TO (10,20,30,40,50,60,80,100,120,140,150,160,170,190,210,230,
     1       240,250), NPROB
C
C     HELICAL VALLEY FUNCTION.
C
   10 CONTINUE
      X(1) = -ONE
      X(2) = ZERO
      X(3) = ZERO
      GO TO 270
C
C     BIGGS EXP6 FUNCTION.
C
   20 CONTINUE
      X(1) = ONE
      X(2) = TWO
      X(3) = ONE
      X(4) = ONE
      X(5) = ONE
      X(6) = ONE
      GO TO 270
C
C     GAUSSIAN FUNCTION.
C
   30 CONTINUE
      X(1) = C1
      X(2) = ONE
      X(3) = ZERO
      GO TO 270
C
C     POWELL BADLY SCALED FUNCTION.
C
   40 CONTINUE
      X(1) = ZERO
      X(2) = ONE
      GO TO 270
C
C     BOX 3-DIMENSIONAL FUNCTION.
C
   50 CONTINUE
      X(1) = ZERO
      X(2) = TEN
      X(3) = TWENTY
      GO TO 270
C
C     VARIABLY DIMENSIONED FUNCTION.
C
   60 CONTINUE
      H = ONE/DFLOAT(N)
      DO 70 J = 1, N
         X(J) = ONE - DFLOAT(J)*H
   70    CONTINUE
      GO TO 270
C
C     WATSON FUNCTION.
C
   80 CONTINUE
      DO 90 J = 1, N
         X(J) = ZERO
   90    CONTINUE
      GO TO 270
C
C     PENALTY FUNCTION I.
C
  100 CONTINUE
      DO 110 J = 1, N
         X(J) = DFLOAT(J)
  110    CONTINUE
      GO TO 270
C
C     PENALTY FUNCTION II.
C
  120 CONTINUE
      DO 130 J = 1, N
         X(J) = HALF
  130    CONTINUE
      GO TO 270
C
C     BROWN BADLY SCALED FUNCTION.
C
  140 CONTINUE
      X(1) = ONE
      X(2) = ONE
      GO TO 270
C
C     BROWN AND DENNIS FUNCTION.
C
  150 CONTINUE
      X(1) = TWNTF
      X(2) = FIVE
      X(3) = -FIVE
      X(4) = -ONE
      GO TO 270
C
C     GULF RESEARCH AND DEVELOPMENT FUNCTION.
C
  160 CONTINUE
      X(1) = FIVE
      X(2) = C2
      X(3) = C3
      GO TO 270
C
C     TRIGONOMETRIC FUNCTION.
C
  170 CONTINUE
      H = ONE/DFLOAT(N)
      DO 180 J = 1, N
         X(J) = H
  180    CONTINUE
      GO TO 270
C
C     EXTENDED ROSENBROCK FUNCTION.
C
  190 CONTINUE
      DO 200 J = 1, N, 2
         X(J) = -C4
         X(J+1) = ONE
  200    CONTINUE
      GO TO 270
C
C     EXTENDED POWELL SINGULAR FUNCTION.
C
  210 CONTINUE
      DO 220 J = 1, N, 4
         X(J) = THREE
         X(J+1) = -ONE
         X(J+2) = ZERO
         X(J+3) = ONE
  220    CONTINUE
      GO TO 270
C
C     BEALE FUNCTION.
C
  230 CONTINUE
      X(1) = ONE
      X(2) = ONE
      GO TO 270
C
C     WOOD FUNCTION.
C
  240 CONTINUE
      X(1) = -THREE
      X(2) = -ONE
      X(3) = -THREE
      X(4) = -ONE
      GO TO 270
C
C     CHEBYQUAD FUNCTION.
C
  250 CONTINUE
      H = ONE/DFLOAT(N+1)
      DO 260 J = 1, N
         X(J) = DFLOAT(J)*H
  260    CONTINUE
  270 CONTINUE
C
C     COMPUTE MULTIPLE OF INITIAL POINT.
C
      IF (FACTOR .EQ. ONE) GO TO 320
      IF (NPROB .EQ. 7) GO TO 290
         DO 280 J = 1, N
            X(J) = FACTOR*X(J)
  280       CONTINUE
         GO TO 310
  290 CONTINUE
         DO 300 J = 1, N
            X(J) = FACTOR
  300       CONTINUE
  310 CONTINUE
  320 CONTINUE
      RETURN
C
C     LAST CARD OF SUBROUTINE INITPT.
C
      END
      SUBROUTINE OBJFCN(N,X,F,NPROB)                                    00000010
      INTEGER N,NPROB
      DOUBLE PRECISION F
      DOUBLE PRECISION X(N)
C     **********
C
C     SUBROUTINE OBJFCN
C
C     THIS SUBROUTINE DEFINES THE OBJECTIVE FUNCTIONS OF EIGHTEEN
C     NONLINEAR UNCONSTRAINED MINIMIZATION PROBLEMS. THE VALUES
C     OF N FOR FUNCTIONS 1,2,3,4,5,10,11,12,16 AND 17 ARE
C     3,6,3,2,3,2,4,3,2 AND 4, RESPECTIVELY.
C     FOR FUNCTION 7, N MAY BE 2 OR GREATER BUT IS USUALLY 6 OR 9.
C     FOR FUNCTIONS 6,8,9,13,14,15 AND 18 N MAY BE VARIABLE,
C     HOWEVER IT MUST BE EVEN FOR FUNCTION 14, A MULTIPLE OF 4 FOR
C     FUNCTION 15, AND NOT GREATER THAN 50 FOR FUNCTION 18.
C
C     THE SUBROUTINE STATEMENT IS
C
C       SUBROUTINE OBJFCN(N,X,F,NPROB)
C
C     WHERE
C
C       N IS A POSITIVE INTEGER INPUT VARIABLE.
C
C       X IS AN INPUT ARRAY OF LENGTH N.
C
C       F IS AN OUTPUT VARIABLE WHICH CONTAINS THE VALUE OF
C         THE NPROB OBJECTIVE FUNCTION EVALUATED AT X.
C
C       NPROB IS A POSITIVE INTEGER INPUT VARIABLE WHICH DEFINES THE
C         NUMBER OF THE PROBLEM. NPROB MUST NOT EXCEED 18.
C
C     SUBPROGRAMS CALLED
C
C       FORTRAN-SUPPLIED ... DABS,DATAN,DCOS,DEXP,DLOG,DSIGN,DSIN,
C                            DSQRT
C
C     MINPACK. VERSION OF JULY 1978.
C     BURTON S. GARBOW, KENNETH E. HILLSTROM, JORGE J. MORE
C
C     **********
      INTEGER I,IEV,IVAR,J
      DOUBLE PRECISION AP,ARG,BP,C2PDM6,CP0001,CP1,CP2,CP25,CP5,C1P5,
     1                 C2P25,C2P625,C3P5,C25,C29,C90,C100,C10000,
     2                 C1PD6,D1,D2,EIGHT,FIFTY,FIVE,FOUR,ONE,R,S1,S2,
     3                 S3,T,T1,T2,T3,TEN,TH,THREE,TPI,TWO,ZERO
      DOUBLE PRECISION FVEC(50),Y(15)
      DOUBLE PRECISION DFLOAT
      DATA ZERO,ONE,TWO,THREE,FOUR,FIVE,EIGHT,TEN,FIFTY
     1     /0.0D0,1.0D0,2.0D0,3.0D0,4.0D0,5.0D0,8.0D0,1.0D1,5.0D1/
      DATA C2PDM6,CP0001,CP1,CP2,CP25,CP5,C1P5,C2P25,C2P625,C3P5,C25,
     1     C29,C90,C100,C10000,C1PD6
     2     /2.0D-6,1.0D-4,1.0D-1,2.0D-1,2.5D-1,5.0D-1,1.5D0,2.25D0,
     3      2.625D0,3.5D0,2.5D1,2.9D1,9.0D1,1.0D2,1.0D4,1.0D6/
      DATA AP,BP /1.0D-5,1.0D0/
      DATA Y(1),Y(2),Y(3),Y(4),Y(5),Y(6),Y(7),Y(8),Y(9),Y(10),Y(11),
     1     Y(12),Y(13),Y(14),Y(15)
     2     /9.0D-4,4.4D-3,1.75D-2,5.4D-2,1.295D-1,2.42D-1,3.521D-1,
     3      3.989D-1,3.521D-1,2.42D-1,1.295D-1,5.4D-2,1.75D-2,4.4D-3,
     4      9.0D-4/
      DFLOAT(IVAR) = IVAR
C
C     FUNCTION ROUTINE SELECTOR.
C
      GO TO (10,20,40,60,70,90,110,150,170,200,210,230,250,280,300,
     1       320,330,340), NPROB
C
C     HELICAL VALLEY FUNCTION.
C
   10 CONTINUE
      TPI = EIGHT*DATAN(ONE)
      TH = DSIGN(CP25,X(2))
      IF (X(1) .GT. ZERO) TH = DATAN(X(2)/X(1))/TPI
      IF (X(1) .LT. ZERO) TH = DATAN(X(2)/X(1))/TPI + CP5
      ARG = X(1)**2 + X(2)**2
      R = DSQRT(ARG)
      T = X(3) - TEN*TH
      F = C100*(T**2 + (R - ONE)**2) + X(3)**2
      GO TO 390
C
C     BIGGS EXP6 FUNCTION.
C
   20 CONTINUE
      F = ZERO
      DO 30 I = 1, 13
         D1 = DFLOAT(I)/TEN
         D2 = DEXP(-D1) - FIVE*DEXP(-TEN*D1) + THREE*DEXP(-FOUR*D1)
         S1 = DEXP(-D1*X(1))
         S2 = DEXP(-D1*X(2))
         S3 = DEXP(-D1*X(5))
         T = X(3)*S1 - X(4)*S2 + X(6)*S3 - D2
         F = F + T**2
   30    CONTINUE
      GO TO 390
C
C     GAUSSIAN FUNCTION.
C
   40 CONTINUE
      F = ZERO
      DO 50 I = 1, 15
         D1 = CP5*DFLOAT(I-1)
         D2 = C3P5 - D1 - X(3)
         ARG = -CP5*X(2)*D2**2
         R = DEXP(ARG)
         T = X(1)*R - Y(I)
         F = F + T**2
   50    CONTINUE
      GO TO 390
C
C     POWELL BADLY SCALED FUNCTION.
C
   60 CONTINUE
      T1 = C10000*X(1)*X(2) - ONE
      S1 = DEXP(-X(1))
      S2 = DEXP(-X(2))
      T2 = S1 + S2 - ONE - CP0001
      F = T1**2 + T2**2
      GO TO 390
C
C     BOX 3-DIMENSIONAL FUNCTION.
C
   70 CONTINUE
      F = ZERO
      DO 80 I = 1, 10
         D1 = DFLOAT(I)
         D2 = D1/TEN
         S1 = DEXP(-D2*X(1))
         S2 = DEXP(-D2*X(2))
         S3 = DEXP(-D2) - DEXP(-D1)
         T = S1 - S2 - S3*X(3)
         F = F + T**2
   80    CONTINUE
      GO TO 390
C
C     VARIABLY DIMENSIONED FUNCTION.
C
   90 CONTINUE
      T1 = ZERO
      T2 = ZERO
      DO 100 J = 1, N
         T1 = T1 + DFLOAT(J)*(X(J) - ONE)
         T2 = T2 + (X(J) - ONE)**2
  100    CONTINUE
      F = T2 + T1**2*(ONE + T1**2)
      GO TO 390
C
C     WATSON FUNCTION.
C
  110 CONTINUE
      F = ZERO
      DO 140 I = 1, 29
         D1 = DFLOAT(I)/C29
         S1 = ZERO
         D2 = ONE
         DO 120 J = 2, N
            S1 = S1 + DFLOAT(J-1)*D2*X(J)
            D2 = D1*D2
  120       CONTINUE
         S2 = ZERO
         D2 = ONE
         DO 130 J = 1, N
            S2 = S2 + D2*X(J)
            D2 = D1*D2
  130       CONTINUE
         T = S1 - S2**2 - ONE
         F = F + T**2
  140    CONTINUE
      T1 = X(2) - X(1)**2 - ONE
      F = F + X(1)**2 + T1**2
      GO TO 390
C
C     PENALTY FUNCTION I.
C
  150 CONTINUE
      T1 = -CP25
      T2 = ZERO
      DO 160 J = 1, N
         T1 = T1 + X(J)**2
         T2 = T2 + (X(J) - ONE)**2
  160    CONTINUE
      F = AP*T2 + BP*T1**2
      GO TO 390
C
C     PENALTY FUNCTION II.
C
  170 CONTINUE
      T1 = -ONE
      T2 = ZERO
      T3 = ZERO
      D1 = DEXP(CP1)
      D2 = ONE
      DO 190 J = 1, N
         T1 = T1 + DFLOAT(N-J+1)*X(J)**2
         S1 = DEXP(X(J)/TEN)
         IF (J .EQ. 1) GO TO 180
         S3 = S1 + S2 - D2*(D1 + ONE)
         T2 = T2 + S3**2
         T3 = T3 + (S1 - ONE/D1)**2
  180    CONTINUE
         S2 = S1
         D2 = D1*D2
  190    CONTINUE
      F = AP*(T2 + T3) + BP*(T1**2 + (X(1) - CP2)**2)
      GO TO 390
C
C     BROWN BADLY SCALED FUNCTION.
C
  200 CONTINUE
      T1 = X(1) - C1PD6
      T2 = X(2) - C2PDM6
      T3 = X(1)*X(2) - TWO
      F = T1**2 + T2**2 + T3**2
      GO TO 390
C
C     BROWN AND DENNIS FUNCTION.
C
  210 CONTINUE
      F = ZERO
      DO 220 I = 1, 20
         D1 = DFLOAT(I)/FIVE
         D2 = DSIN(D1)
         T1 = X(1) + D1*X(2) - DEXP(D1)
         T2 = X(3) + D2*X(4) - DCOS(D1)
         T = T1**2 + T2**2
         F = F + T**2
  220    CONTINUE
      GO TO 390
C
C     GULF RESEARCH AND DEVELOPMENT FUNCTION.
C
  230 CONTINUE
      F = ZERO
      D1 = TWO/THREE
      DO 240 I = 1, 99
         ARG = DFLOAT(I)/C100
         R = DABS((-FIFTY*DLOG(ARG))**D1+C25-X(2))
         T1 = R**X(3)/X(1)
         T2 = DEXP(-T1)
         T = T2 - ARG
         F = F + T**2
  240    CONTINUE
      GO TO 390
C
C     TRIGONOMETRIC FUNCTION.
C
  250 CONTINUE
      S1 = ZERO
      DO 260 J = 1, N
         S1 = S1 + DCOS(X(J))
  260    CONTINUE
      F = ZERO
      DO 270 J = 1, N
         T = DFLOAT(N+J) - DSIN(X(J)) - S1 - DFLOAT(J)*DCOS(X(J))
         F = F + T**2
  270    CONTINUE
      GO TO 390
C
C     EXTENDED ROSENBROCK FUNCTION.
C
  280 CONTINUE
      F = ZERO
      DO 290 J = 1, N, 2
         T1 = ONE - X(J)
         T2 = TEN*(X(J+1) - X(J)**2)
         F = F + T1**2 + T2**2
  290    CONTINUE
      GO TO 390
C
C     EXTENDED POWELL FUNCTION.
C
  300 CONTINUE
      F = ZERO
      DO 310 J = 1, N, 4
         T = X(J) + TEN*X(J+1)
         T1 = X(J+2) - X(J+3)
         S1 = FIVE*T1
         T2 = X(J+1) - TWO*X(J+2)
         S2 = T2**3
         T3 = X(J) - X(J+3)
         S3 = TEN*T3**3
         F = F + T**2 + S1*T1 + S2*T2 + S3*T3
  310    CONTINUE
      GO TO 390
C
C     BEALE FUNCTION.
C
  320 CONTINUE
      S1 = ONE - X(2)
      T1 = C1P5 - X(1)*S1
      S2 = ONE - X(2)**2
      T2 = C2P25 - X(1)*S2
      S3 = ONE - X(2)**3
      T3 = C2P625 - X(1)*S3
      F = T1**2 + T2**2 + T3**2
      GO TO 390
C
C     WOOD FUNCTION.
C
  330 CONTINUE
      S1 = X(2) - X(1)**2
      S2 = ONE - X(1)
      S3 = X(2) - ONE
      T1 = X(4) - X(3)**2
      T2 = ONE - X(3)
      T3 = X(4) - ONE
      F = C100*S1**2 + S2**2 + C90*T1**2 + T2**2 + TEN*(S3 + T3)**2
     1    + (S3 - T3)**2/TEN
      GO TO 390
C
C     CHEBYQUAD FUNCTION.
C
  340 CONTINUE
      DO 350 I = 1, N
         FVEC(I) = ZERO
  350    CONTINUE
      DO 370 J = 1, N
         T1 = ONE
         T2 = TWO*X(J) - ONE
         T = TWO*T2
         DO 360 I = 1, N
            FVEC(I) = FVEC(I) + T2
            TH = T*T2 - T1
            T1 = T2
            T2 = TH
  360       CONTINUE
  370    CONTINUE
      F = ZERO
      D1 = ONE/DFLOAT(N)
      IEV = -1
      DO 380 I = 1, N
         T = D1*FVEC(I)
         IF (IEV .GT. 0) T = T + ONE/(DFLOAT(I)**2 - ONE)
         F = F + T**2
         IEV = -IEV
  380    CONTINUE
  390 CONTINUE
      RETURN
C
C     LAST CARD OF SUBROUTINE OBJFCN.
C
      END
      SUBROUTINE GRDFCN(N,X,G,NPROB)                                    00000010
      INTEGER N,NPROB
      DOUBLE PRECISION X(N),G(N)
C     **********
C
C     SUBROUTINE GRDFCN
C
C     THIS SUBROUTINE DEFINES THE GRADIENT VECTORS OF EIGHTEEN
C     NONLINEAR UNCONSTRAINED MINIMIZATION PROBLEMS. THE PROBLEM
C     DIMENSIONS ARE AS DESCRIBED IN THE PROLOGUE COMMENTS OF OBJFCN.
C
C     THE SUBROUTINE STATEMENT IS
C
C       SUBROUTINE GRDFCN(N,X,G,NPROB)
C
C     WHERE
C
C       N IS A POSITIVE INTEGER INPUT VARIABLE.
C
C       X IS AN INPUT ARRAY OF LENGTH N.
C
C       G IS AN OUTPUT ARRAY OF LENGTH N WHICH CONTAINS THE COMPONENTS
C         OF THE GRADIENT VECTOR OF THE NPROB OBJECTIVE FUNCTION
C         EVALUATED AT X.
C
C       NPROB IS A POSITIVE INTEGER INPUT VARIABLE WHICH DEFINES THE
C         NUMBER OF THE PROBLEM. NPROB MUST NOT EXCEED 18.
C
C     SUBPROGRAMS CALLED
C
C       FORTRAN-SUPPLIED ... DABS,DATAN,DCOS,DEXP,DLOG,DSIGN,DSIN,
C                            DSQRT
C
C     MINPACK. VERSION OF JULY 1978.
C     BURTON S. GARBOW, KENNETH E. HILLSTROM, JORGE J. MORE
C
C     **********
      INTEGER I,IEV,IVAR,J
      DOUBLE PRECISION AP,ARG,BP,C2PDM6,CP0001,CP1,CP2,CP25,CP5,C1P5,
     1                 C2P25,C2P625,C3P5,C19P8,C20P2,C25,C29,C100,
     2                 C180,C200,C10000,C1PD6,D1,D2,EIGHT,FIFTY,FIVE,
     3                 FOUR,ONE,R,S1,S2,S3,T,T1,T2,T3,TEN,TH,THREE,
     4                 TPI,TWENTY,TWO,ZERO
      DOUBLE PRECISION FVEC(50),Y(15)
      DOUBLE PRECISION DFLOAT
      DATA ZERO,ONE,TWO,THREE,FOUR,FIVE,EIGHT,TEN,TWENTY,FIFTY
     1     /0.0D0,1.0D0,2.0D0,3.0D0,4.0D0,5.0D0,8.0D0,1.0D1,2.0D1,
     2      5.0D1/
      DATA C2PDM6,CP0001,CP1,CP2,CP25,CP5,C1P5,C2P25,C2P625,C3P5,
     1     C19P8,C20P2,C25,C29,C100,C180,C200,C10000,C1PD6
     2     /2.0D-6,1.0D-4,1.0D-1,2.0D-1,2.5D-1,5.0D-1,1.5D0,2.25D0,
     3      2.625D0,3.5D0,1.98D1,2.02D1,2.5D1,2.9D1,1.0D2,1.8D2,2.0D2,
     4      1.0D4,1.0D6/
      DATA AP,BP /1.0D-5,1.0D0/
      DATA Y(1),Y(2),Y(3),Y(4),Y(5),Y(6),Y(7),Y(8),Y(9),Y(10),Y(11),
     1     Y(12),Y(13),Y(14),Y(15)
     2     /9.0D-4,4.4D-3,1.75D-2,5.4D-2,1.295D-1,2.42D-1,3.521D-1,
     3      3.989D-1,3.521D-1,2.42D-1,1.295D-1,5.4D-2,1.75D-2,4.4D-3,
     4      9.0D-4/
      DFLOAT(IVAR) = IVAR
C
C     GRADIENT ROUTINE SELECTOR.
C
      GO TO (10,20,50,70,80,100,130,190,220,260,270,290,310,350,370,
     1       390,400,410), NPROB
C
C     HELICAL VALLEY FUNCTION.
C
   10 CONTINUE
      TPI = EIGHT*DATAN(ONE)
      TH = DSIGN(CP25,X(2))
      IF (X(1) .GT. ZERO) TH = DATAN(X(2)/X(1))/TPI
      IF (X(1) .LT. ZERO) TH = DATAN(X(2)/X(1))/TPI + CP5
      ARG = X(1)**2 + X(2)**2
      R = DSQRT(ARG)
      T = X(3) - TEN*TH
      S1 = TEN*T/(TPI*ARG)
      G(1) = C200*(X(1) - X(1)/R + X(2)*S1)
      G(2) = C200*(X(2) - X(2)/R - X(1)*S1)
      G(3) = TWO*(C100*T + X(3))
      GO TO 490
C
C     BIGGS EXP6 FUNCTION.
C
   20 CONTINUE
      DO 30 J = 1, N
         G(J) = ZERO
   30    CONTINUE
      DO 40 I = 1, 13
         D1 = DFLOAT(I)/TEN
         D2 = DEXP(-D1) - FIVE*DEXP(-TEN*D1) + THREE*DEXP(-FOUR*D1)
         S1 = DEXP(-D1*X(1))
         S2 = DEXP(-D1*X(2))
         S3 = DEXP(-D1*X(5))
         T = X(3)*S1 - X(4)*S2 + X(6)*S3 - D2
         TH = D1*T
         G(1) = G(1) - S1*TH
         G(2) = G(2) + S2*TH
         G(3) = G(3) + S1*T
         G(4) = G(4) - S2*T
         G(5) = G(5) - S3*TH
         G(6) = G(6) + S3*T
   40    CONTINUE
      G(1) = TWO*X(3)*G(1)
      G(2) = TWO*X(4)*G(2)
      G(3) = TWO*G(3)
      G(4) = TWO*G(4)
      G(5) = TWO*X(6)*G(5)
      G(6) = TWO*G(6)
      GO TO 490
C
C     GAUSSIAN FUNCTION.
C
   50 CONTINUE
      G(1) = ZERO
      G(2) = ZERO
      G(3) = ZERO
      DO 60 I = 1, 15
         D1 = CP5*DFLOAT(I-1)
         D2 = C3P5 - D1 - X(3)
         ARG = -CP5*X(2)*D2**2
         R = DEXP(ARG)
         T = X(1)*R - Y(I)
         S1 = R*T
         S2 = D2*S1
         G(1) = G(1) + S1
         G(2) = G(2) - D2*S2
         G(3) = G(3) + S2
   60    CONTINUE
      G(1) = TWO*G(1)
      G(2) = X(1)*G(2)
      G(3) = TWO*X(1)*X(2)*G(3)
      GO TO 490
C
C     POWELL BADLY SCALED FUNCTION.
C
   70 CONTINUE
      T1 = C10000*X(1)*X(2) - ONE
      S1 = DEXP(-X(1))
      S2 = DEXP(-X(2))
      T2 = S1 + S2 - ONE - CP0001
      G(1) = TWO*(C10000*X(2)*T1 - S1*T2)
      G(2) = TWO*(C10000*X(1)*T1 - S2*T2)
      GO TO 490
C
C     BOX 3-DIMENSIONAL FUNCTION.
C
   80 CONTINUE
      G(1) = ZERO
      G(2) = ZERO
      G(3) = ZERO
      DO 90 I = 1, 10
         D1 = DFLOAT(I)
         D2 = D1/TEN
         S1 = DEXP(-D2*X(1))
         S2 = DEXP(-D2*X(2))
         S3 = DEXP(-D2) - DEXP(-D1)
         T = S1 - S2 - S3*X(3)
         TH = D2*T
         G(1) = G(1) - S1*TH
         G(2) = G(2) + S2*TH
         G(3) = G(3) - S3*T
   90    CONTINUE
      G(1) = TWO*G(1)
      G(2) = TWO*G(2)
      G(3) = TWO*G(3)
      GO TO 490
C
C     VARIABLY DIMENSIONED FUNCTION.
C
  100 CONTINUE
      T1 = ZERO
      DO 110 J = 1, N
         T1 = T1 + DFLOAT(J)*(X(J) - ONE)
  110    CONTINUE
      T = T1*(ONE + TWO*T1**2)
      DO 120 J = 1, N
         G(J) = TWO*(X(J) - ONE + DFLOAT(J)*T)
  120    CONTINUE
      GO TO 490
C
C     WATSON FUNCTION.
C
  130 CONTINUE
      DO 140 J = 1, N
         G(J) = ZERO
  140    CONTINUE
      DO 180 I = 1, 29
         D1 = DFLOAT(I)/C29
         S1 = ZERO
         D2 = ONE
         DO 150 J = 2, N
            S1 = S1 + DFLOAT(J-1)*D2*X(J)
            D2 = D1*D2
  150       CONTINUE
         S2 = ZERO
         D2 = ONE
         DO 160 J = 1, N
            S2 = S2 + D2*X(J)
            D2 = D1*D2
  160       CONTINUE
         T = S1 - S2**2 - ONE
         S3 = TWO*D1*S2
         D2 = TWO/D1
         DO 170 J = 1, N
            G(J) = G(J) + D2*(DFLOAT(J-1) - S3)*T
            D2 = D1*D2
  170       CONTINUE
  180    CONTINUE
      T1 = X(2) - X(1)**2 - ONE
      G(1) = G(1) + X(1)*(TWO - FOUR*T1)
      G(2) = G(2) + TWO*T1
      GO TO 490
C
C     PENALTY FUNCTION I.
C
  190 CONTINUE
      T1 = -CP25
      DO 200 J = 1, N
         T1 = T1 + X(J)**2
  200    CONTINUE
      D1 = TWO*AP
      TH = FOUR*BP*T1
      DO 210 J = 1, N
         G(J) = D1*(X(J) - ONE) + X(J)*TH
  210    CONTINUE
      GO TO 490
C
C     PENALTY FUNCTION II.
C
  220 CONTINUE
      T1 = -ONE
      DO 230 J = 1, N
         T1 = T1 + DFLOAT(N-J+1)*X(J)**2
  230    CONTINUE
      D1 = DEXP(CP1)
      D2 = ONE
      TH = FOUR*BP*T1
      DO 250 J = 1, N
         G(J) = DFLOAT(N-J+1)*X(J)*TH
         S1 = DEXP(X(J)/TEN)
         IF (J .EQ. 1) GO TO 240
         S3 = S1 + S2 - D2*(D1 + ONE)
         G(J) = G(J) + AP*S1*(S3 + S1 - ONE/D1)/FIVE
         G(J-1) = G(J-1) + AP*S2*S3/FIVE
  240    CONTINUE
         S2 = S1
         D2 = D1*D2
  250    CONTINUE
      G(1) = G(1) + TWO*BP*(X(1) - CP2)
      GO TO 490
C
C     BROWN BADLY SCALED FUNCTION.
C
  260 CONTINUE
      T1 = X(1) - C1PD6
      T2 = X(2) - C2PDM6
      T3 = X(1)*X(2) - TWO
      G(1) = TWO*(T1 + X(2)*T3)
      G(2) = TWO*(T2 + X(1)*T3)
      GO TO 490
C
C     BROWN AND DENNIS FUNCTION.
C
  270 CONTINUE
      G(1) = ZERO
      G(2) = ZERO
      G(3) = ZERO
      G(4) = ZERO
      DO 280 I = 1, 20
         D1 = DFLOAT(I)/FIVE
         D2 = DSIN(D1)
         T1 = X(1) + D1*X(2) - DEXP(D1)
         T2 = X(3) + D2*X(4) - DCOS(D1)
         T = T1**2 + T2**2
         S1 = T1*T
         S2 = T2*T
         G(1) = G(1) + S1
         G(2) = G(2) + D1*S1
         G(3) = G(3) + S2
         G(4) = G(4) + D2*S2
  280    CONTINUE
      G(1) = FOUR*G(1)
      G(2) = FOUR*G(2)
      G(3) = FOUR*G(3)
      G(4) = FOUR*G(4)
      GO TO 490
C
C     GULF RESEARCH AND DEVELOPMENT FUNCTION.
C
  290 CONTINUE
      G(1) = ZERO
      G(2) = ZERO
      G(3) = ZERO
      D1 = TWO/THREE
      DO 300 I = 1, 99
         ARG = DFLOAT(I)/C100
         R = DABS((-FIFTY*DLOG(ARG))**D1+C25-X(2))
         T1 = R**X(3)/X(1)
         T2 = DEXP(-T1)
         T = T2 - ARG
         S1 = T1*T2*T
         G(1) = G(1) + S1
         G(2) = G(2) + S1/R
         G(3) = G(3) - S1*DLOG(R)
  300    CONTINUE
      G(1) = TWO*G(1)/X(1)
      G(2) = TWO*X(3)*G(2)
      G(3) = TWO*G(3)
      GO TO 490
C
C     TRIGONOMETRIC FUNCTION.
C
  310 CONTINUE
      S1 = ZERO
      DO 320 J = 1, N
         G(J) = DCOS(X(J))
         S1 = S1 + G(J)
  320    CONTINUE
      S2 = ZERO
      DO 330 J = 1, N
         TH = DSIN(X(J))
         T = DFLOAT(N+J) - TH - S1 - DFLOAT(J)*G(J)
         S2 = S2 + T
         G(J) = (DFLOAT(J)*TH - G(J))*T
  330    CONTINUE
      DO 340 J = 1, N
         G(J) = TWO*(G(J) + DSIN(X(J))*S2)
  340    CONTINUE
      GO TO 490
C
C     EXTENDED ROSENBROCK FUNCTION.
C
  350 CONTINUE
      DO 360 J = 1, N, 2
         T1 = ONE - X(J)
         G(J+1) = C200*(X(J+1) - X(J)**2)
         G(J) = -TWO*(X(J)*G(J+1) + T1)
  360    CONTINUE
      GO TO 490
C
C     EXTENDED POWELL FUNCTION.
C
  370 CONTINUE
      DO 380 J = 1, N, 4
         T = X(J) + TEN*X(J+1)
         T1 = X(J+2) - X(J+3)
         S1 = FIVE*T1
         T2 = X(J+1) - TWO*X(J+2)
         S2 = FOUR*T2**3
         T3 = X(J) - X(J+3)
         S3 = TWENTY*T3**3
         G(J) = TWO*(T + S3)
         G(J+1) = TWENTY*T + S2
         G(J+2) = TWO*(S1 - S2)
         G(J+3) = -TWO*(S1 + S3)
  380    CONTINUE
      GO TO 490
C
C     BEALE FUNCTION.
C
  390 CONTINUE
      S1 = ONE - X(2)
      T1 = C1P5 - X(1)*S1
      S2 = ONE - X(2)**2
      T2 = C2P25 - X(1)*S2
      S3 = ONE - X(2)**3
      T3 = C2P625 - X(1)*S3
      G(1) = -TWO*(S1*T1 + S2*T2 + S3*T3)
      G(2) = TWO*X(1)*(T1 + X(2)*(TWO*T2 + THREE*X(2)*T3))
      GO TO 490
C
C     WOOD FUNCTION.
C
  400 CONTINUE
      S1 = X(2) - X(1)**2
      S2 = ONE - X(1)
      S3 = X(2) - ONE
      T1 = X(4) - X(3)**2
      T2 = ONE - X(3)
      T3 = X(4) - ONE
      G(1) = -TWO*(C200*X(1)*S1 + S2)
      G(2) = C200*S1 + C20P2*S3 + C19P8*T3
      G(3) = -TWO*(C180*X(3)*T1 + T2)
      G(4) = C180*T1 + C20P2*T3 + C19P8*S3
      GO TO 490
C
C     CHEBYQUAD FUNCTION.
C
  410 CONTINUE
      DO 420 I = 1, N
         FVEC(I) = ZERO
  420    CONTINUE
      DO 440 J = 1, N
         T1 = ONE
         T2 = TWO*X(J) - ONE
         T = TWO*T2
         DO 430 I = 1, N
            FVEC(I) = FVEC(I) + T2
            TH = T*T2 - T1
            T1 = T2
            T2 = TH
  430       CONTINUE
  440    CONTINUE
      D1 = ONE/DFLOAT(N)
      IEV = -1
      DO 450 I = 1, N
         FVEC(I) = D1*FVEC(I)
         IF (IEV .GT. 0) FVEC(I) = FVEC(I) + ONE/(DFLOAT(I)**2 - ONE)
         IEV = -IEV
  450    CONTINUE
      DO 470 J = 1, N
         G(J) = ZERO
         T1 = ONE
         T2 = TWO*X(J) - ONE
         T = TWO*T2
         S1 = ZERO
         S2 = TWO
         DO 460 I = 1, N
            G(J) = G(J) + FVEC(I)*S2
            TH = FOUR*T2 + T*S2 - S1
            S1 = S2
            S2 = TH
            TH = T*T2 - T1
            T1 = T2
            T2 = TH
  460       CONTINUE
  470    CONTINUE
      D2 = TWO*D1
      DO 480 J = 1, N
         G(J) = D2*G(J)
  480    CONTINUE
  490 CONTINUE
      RETURN
C
C     LAST CARD OF SUBROUTINE GRDFCN.
C
      END
C ===== 5. SINGLE PRECISION TESTING AIDS FOR NONLINEAR EQUATIONS.
      SUBROUTINE INITPT(N,X,NPROB,FACTOR)                               00000010
      INTEGER N,NPROB
      REAL FACTOR
      REAL X(N)
C     **********
C
C     SUBROUTINE INITPT
C
C     THIS SUBROUTINE SPECIFIES THE STANDARD STARTING POINTS FOR
C     THE FUNCTIONS DEFINED BY SUBROUTINES COMFCN AND VECFCN. THE
C     SUBROUTINE RETURNS IN X A MULTIPLE (FACTOR) OF THE STANDARD
C     STARTING POINT. FOR THE SIXTH FUNCTION THE STANDARD STARTING
C     POINT IS ZERO, SO IN THIS CASE, IF FACTOR IS NOT UNITY, THEN
C     THE SUBROUTINE RETURNS THE VECTOR  X(J) = FACTOR, J=1,...,N.
C
C     THE SUBROUTINE STATEMENT IS
C
C       SUBROUTINE INITPT(N,X,NPROB,FACTOR)
C
C     WHERE
C
C       N IS A POSITIVE INTEGER VARIABLE.
C
C       X IS AN OUTPUT ARRAY OF LENGTH N WHICH CONTAINS THE
C         STANDARD STARTING POINT FOR PROBLEM NPROB MULTIPLIED BY
C         FACTOR.
C
C       NPROB IS A POSITIVE INTEGER VARIABLE WHICH DEFINES THE
C         NUMBER OF THE PROBLEM. NPROB MUST NOT EXCEED 14.
C
C       FACTOR SPECIFIES THE MULTIPLE OF THE STANDARD STARTING
C         POINT. IF FACTOR IS UNITY, NO MULTIPLICATION IS PERFORMED.
C
C     MINPACK. VERSION OF JULY 1978.
C     BURTON S. GARBOW, KENNETH E. HILLSTROM, JORGE J. MORE
C
C     **********
      INTEGER IVAR,J
      REAL C1,H,HALF,ONE,THREE,TJ,ZERO
      REAL FLOAT
      DATA ZERO,HALF,ONE,THREE,C1 /0.0E0,5.0E-1,1.0E0,3.0E0,1.2E0/
      FLOAT(IVAR) = IVAR
C
C     SELECTION OF INITIAL POINT.
C
      GO TO (10,20,30,40,50,60,80,100,120,120,140,160,180,180), NPROB
C
C     ROSENBROCK FUNCTION.
C
   10 CONTINUE
      X(1) = -C1
      X(2) = ONE
      GO TO 200
C
C     POWELL SINGULAR FUNCTION.
C
   20 CONTINUE
      X(1) = THREE
      X(2) = -ONE
      X(3) = ZERO
      X(4) = ONE
      GO TO 200
C
C     POWELL BADLY SCALED FUNCTION.
C
   30 CONTINUE
      X(1) = ZERO
      X(2) = ONE
      GO TO 200
C
C     WOOD FUNCTION.
C
   40 CONTINUE
      X(1) = -THREE
      X(2) = -ONE
      X(3) = -THREE
      X(4) = -ONE
      GO TO 200
C
C     HELICAL VALLEY FUNCTION.
C
   50 CONTINUE
      X(1) = -ONE
      X(2) = ZERO
      X(3) = ZERO
      GO TO 200
C
C     WATSON FUNCTION.
C
   60 CONTINUE
      DO 70 J = 1, N
         X(J) = ZERO
   70    CONTINUE
      GO TO 200
C
C     CHEBYQUAD FUNCTION.
C
   80 CONTINUE
      H = ONE/FLOAT(N+1)
      DO 90 J = 1, N
         X(J) = FLOAT(J)*H
   90    CONTINUE
      GO TO 200
C
C     BROWN ALMOST-LINEAR FUNCTION.
C
  100 CONTINUE
      DO 110 J = 1, N
         X(J) = HALF
  110    CONTINUE
      GO TO 200
C
C     DISCRETE BOUNDARY VALUE AND INTEGRAL EQUATION FUNCTIONS.
C
  120 CONTINUE
      H = ONE/FLOAT(N+1)
      DO 130 J = 1, N
         TJ = FLOAT(J)*H
         X(J) = TJ*(TJ - ONE)
  130    CONTINUE
      GO TO 200
C
C     TRIGONOMETRIC FUNCTION.
C
  140 CONTINUE
      H = ONE/FLOAT(N)
      DO 150 J = 1, N
         X(J) = H
  150    CONTINUE
      GO TO 200
C
C     VARIABLY DIMENSIONED FUNCTION.
C
  160 CONTINUE
      H = ONE/FLOAT(N)
      DO 170 J = 1, N
         X(J) = ONE - FLOAT(J)*H
  170    CONTINUE
      GO TO 200
C
C     BROYDEN TRIDIAGONAL AND BANDED FUNCTIONS.
C
  180 CONTINUE
      DO 190 J = 1, N
         X(J) = -ONE
  190    CONTINUE
  200 CONTINUE
C
C     COMPUTE MULTIPLE OF INITIAL POINT.
C
      IF (FACTOR .EQ. ONE) GO TO 250
      IF (NPROB .EQ. 6) GO TO 220
         DO 210 J = 1, N
            X(J) = FACTOR*X(J)
  210       CONTINUE
         GO TO 240
  220 CONTINUE
         DO 230 J = 1, N
            X(J) = FACTOR
  230       CONTINUE
  240 CONTINUE
  250 CONTINUE
      RETURN
C
C     LAST CARD OF SUBROUTINE INITPT.
C
      END
      SUBROUTINE COMFCN(N,K,X,FCNK,NPROB)                               00000010
      INTEGER N,K,NPROB
      REAL FCNK
      REAL X(N)
C     **********
C
C     SUBROUTINE COMFCN
C
C     THIS SUBROUTINE DEFINES FOURTEEN TEST FUNCTIONS. THE FIRST
C     FIVE TEST FUNCTIONS ARE OF DIMENSIONS 2,4,2,4,3, RESPECTIVELY,
C     WHILE THE REMAINING TEST FUNCTIONS ARE OF VARIABLE DIMENSION
C     N FOR ANY N GREATER THAN OR EQUAL TO 1 (PROBLEM 6 IS AN
C     EXCEPTION TO THIS, SINCE IT DOES NOT ALLOW N = 1).
C
C     THE SUBROUTINE STATEMENT IS
C
C       SUBROUTINE COMFCN(N,K,X,FCNK,NPROB)
C
C     WHERE
C
C       N IS A POSITIVE INTEGER INPUT VARIABLE.
C
C       K IS A POSITIVE INTEGER INPUT VARIABLE NOT GREATER THAN N.
C
C       X IS AN INPUT ARRAY OF LENGTH N.
C
C       FCNK IS AN OUTPUT VARIABLE WHICH CONTAINS THE VALUE OF
C         THE K-TH COMPONENT OF THE NPROB FUNCTION EVALUATED AT X.
C
C       NPROB IS A POSITIVE INTEGER INPUT VARIABLE WHICH DEFINES THE
C         NUMBER OF THE PROBLEM. NPROB MUST NOT EXCEED 14.
C
C     SUBPROGRAMS REQUIRED
C
C       FORTRAN-SUPPLIED ... ATAN,COS,EXP,SIGN,SIN,SQRT,
C                            MAX0,MIN0,MOD
C
C     MINPACK. VERSION OF JULY 1978.
C     BURTON S. GARBOW, KENNETH E. HILLSTROM, JORGE J. MORE
C
C     **********
      INTEGER I,IVAR,J,K1,K2,KP1,ML,MU
      REAL C1,C2,C3,C4,C5,C6,C7,C8,C9,EIGHT,FIVE,H,ONE,PROD,SUM,SUM1,
     *     SUM2,TEMP,TEMP1,TEMP2,TEN,THREE,TI,TJ,TK,TPI,TWO,ZERO
      REAL FLOAT
      DATA ZERO,ONE,TWO,THREE,FIVE,EIGHT,TEN
     *     /0.0E0,1.0E0,2.0E0,3.0E0,5.0E0,8.0E0,1.0E1/
      DATA C1,C2,C3,C4,C5,C6,C7,C8,C9
     *     /1.0E4,1.0001E0,2.0E2,2.02E1,1.98E1,1.8E2,2.5E-1,5.0E-1,
     *      2.9E1/
      FLOAT(IVAR) = IVAR
C
C     PROBLEM SELECTOR.
C
      GO TO (10,20,30,40,50,70,110,150,200,210,250,270,290,300), NPROB
C
C     ROSENBROCK FUNCTION.
C
   10 CONTINUE
      IF (K .EQ. 1) FCNK = ONE - X(1)
      IF (K .EQ. 2) FCNK = TEN*(X(2) - X(1)**2)
      GO TO 320
C
C     POWELL SINGULAR FUNCTION.
C
   20 CONTINUE
      IF (K .EQ. 1) FCNK = X(1) + TEN*X(2)
      IF (K .EQ. 2) FCNK = SQRT(FIVE)*(X(3) - X(4))
      IF (K .EQ. 3) FCNK = (X(2) - TWO*X(3))**2
      IF (K .EQ. 4) FCNK = SQRT(TEN)*(X(1) - X(4))**2
      GO TO 320
C
C     POWELL BADLY SCALED FUNCTION.
C
   30 CONTINUE
      IF (K .EQ. 1) FCNK = C1*X(1)*X(2) - ONE
      IF (K .EQ. 2) FCNK = EXP(-X(1)) + EXP(-X(2)) - C2
      GO TO 320
C
C     WOOD FUNCTION.
C
   40 CONTINUE
      TEMP1 = X(2) - X(1)**2
      TEMP2 = X(4) - X(3)**2
      IF (K .EQ. 1) FCNK = -C3*X(1)*TEMP1 - (ONE - X(1))
      IF (K .EQ. 2)
     *   FCNK = C3*TEMP1 + C4*(X(2) - ONE) + C5*(X(4) - ONE)
      IF (K .EQ. 3) FCNK = -C6*X(3)*TEMP2 - (ONE - X(3))
      IF (K .EQ. 4)
     *   FCNK = C6*TEMP2 + C4*(X(4) - ONE) + C5*(X(2) - ONE)
      GO TO 320
C
C     HELICAL VALLEY FUNCTION.
C
   50 CONTINUE
      IF (K .NE. 1) GO TO 60
      TPI = EIGHT*ATAN(ONE)
      TEMP1 = SIGN(C7,X(2))
      IF (X(1) .GT. ZERO) TEMP1 = ATAN(X(2)/X(1))/TPI
      IF (X(1) .LT. ZERO) TEMP1 = ATAN(X(2)/X(1))/TPI + C8
      FCNK = TEN*(X(3) - TEN*TEMP1)
   60 CONTINUE
      IF (K .EQ. 2) FCNK = TEN*(SQRT(X(1)**2+X(2)**2) - ONE)
      IF (K .EQ. 3) FCNK = X(3)
      GO TO 320
C
C     WATSON FUNCTION.
C
   70 CONTINUE
      FCNK = ZERO
      DO 100 I = 1, 29
         TI = FLOAT(I)/C9
         SUM1 = ZERO
         TEMP = ONE
         DO 80 J = 2, N
            SUM1 = SUM1 + FLOAT(J-1)*TEMP*X(J)
            TEMP = TI*TEMP
   80       CONTINUE
         SUM2 = ZERO
         TEMP = ONE
         DO 90 J = 1, N
            SUM2 = SUM2 + TEMP*X(J)
            TEMP = TI*TEMP
   90       CONTINUE
         TEMP1 = SUM1 - SUM2**2 - ONE
         TEMP2 = TWO*TI*SUM2
         FCNK = FCNK + TI**(K - 2)*(FLOAT(K-1) - TEMP2)*TEMP1
  100    CONTINUE
      TEMP = X(2) - X(1)**2 - ONE
      IF (K .EQ. 1) FCNK = FCNK + X(1)*(ONE - TWO*TEMP)
      IF (K .EQ. 2) FCNK = FCNK + TEMP
      GO TO 320
C
C     CHEBYQUAD FUNCTION.
C
  110 CONTINUE
      SUM = ZERO
      DO 140 J = 1, N
         TEMP1 = ONE
         TEMP2 = TWO*X(J) - ONE
         TEMP = TWO*TEMP2
         IF (K .LT. 2) GO TO 130
         DO 120 I = 2, K
            TI = TEMP*TEMP2 - TEMP1
            TEMP1 = TEMP2
            TEMP2 = TI
  120       CONTINUE
  130    CONTINUE
         SUM = SUM + TEMP2
  140    CONTINUE
      FCNK = SUM/FLOAT(N)
      IF (MOD(K,2) .EQ. 0) FCNK = FCNK + ONE/(FLOAT(K)**2 - ONE)
      GO TO 320
C
C     BROWN ALMOST-LINEAR FUNCTION.
C
  150 CONTINUE
      IF (K .EQ. N) GO TO 170
         SUM = -FLOAT(N+1)
         DO 160 J = 1, N
            SUM = SUM + X(J)
  160       CONTINUE
         FCNK = X(K) + SUM
         GO TO 190
  170 CONTINUE
         PROD = ONE
         DO 180 J = 1, N
            PROD = X(J)*PROD
  180       CONTINUE
         FCNK = PROD - ONE
  190 CONTINUE
      GO TO 320
C
C     DISCRETE BOUNDARY VALUE FUNCTION.
C
  200 CONTINUE
      H = ONE/FLOAT(N+1)
      TEMP = (X(K) + FLOAT(K)*H + ONE)**3
      TEMP1 = ZERO
      IF (K .NE. 1) TEMP1 = X(K-1)
      TEMP2 = ZERO
      IF (K .NE. N) TEMP2 = X(K+1)
      FCNK = TWO*X(K) - TEMP1 - TEMP2 + TEMP*H**2/TWO
      GO TO 320
C
C     DISCRETE INTEGRAL EQUATION FUNCTION.
C
  210 CONTINUE
      H = ONE/FLOAT(N+1)
      TK = FLOAT(K)*H
      SUM1 = ZERO
      DO 220 J = 1, K
         TJ = FLOAT(J)*H
         TEMP = (X(J) + TJ + ONE)**3
         SUM1 = SUM1 + TJ*TEMP
  220    CONTINUE
      SUM2 = ZERO
      KP1 = K + 1
      IF (N .LT. KP1) GO TO 240
      DO 230 J = KP1, N
         TJ = FLOAT(J)*H
         TEMP = (X(J) + TJ + ONE)**3
         SUM2 = SUM2 + (ONE - TJ)*TEMP
  230    CONTINUE
  240 CONTINUE
      FCNK = X(K) + H*((ONE - TK)*SUM1 + TK*SUM2)/TWO
      GO TO 320
C
C     TRIGONOMETRIC FUNCTION.
C
  250 CONTINUE
      SUM = ZERO
      DO 260 J = 1, N
         SUM = SUM + COS(X(J))
  260    CONTINUE
      FCNK = FLOAT(N+K) - SIN(X(K)) - SUM - FLOAT(K)*COS(X(K))
      GO TO 320
C
C     VARIABLY DIMENSIONED FUNCTION.
C
  270 CONTINUE
      SUM = ZERO
      DO 280 J = 1, N
         SUM = SUM + FLOAT(J)*(X(J) - ONE)
  280    CONTINUE
      TEMP = SUM*(ONE + TWO*SUM**2)
      FCNK = X(K) - ONE + FLOAT(K)*TEMP
      GO TO 320
C
C     BROYDEN TRIDIAGONAL FUNCTION.
C
  290 CONTINUE
      TEMP = (THREE - TWO*X(K))*X(K)
      TEMP1 = ZERO
      IF (K .NE. 1) TEMP1 = X(K-1)
      TEMP2 = ZERO
      IF (K .NE. N) TEMP2 = X(K+1)
      FCNK = TEMP - TEMP1 - TWO*TEMP2 + ONE
      GO TO 320
C
C     BROYDEN BANDED FUNCTION.
C
  300 CONTINUE
      ML = 5
      MU = 1
      K1 = MAX0(1,K-ML)
      K2 = MIN0(K+MU,N)
      TEMP = ZERO
      DO 310 J = K1, K2
         IF (J .NE. K) TEMP = TEMP + X(J)*(ONE + X(J))
  310    CONTINUE
      FCNK = X(K)*(TWO + FIVE*X(K)**2) + ONE - TEMP
  320 CONTINUE
      RETURN
C
C     LAST CARD OF SUBROUTINE COMFCN.
C
      END
      SUBROUTINE VECFCN(N,X,FVEC,NPROB)                                 00000010
      INTEGER N,NPROB
      REAL X(N),FVEC(N)
C     **********
C
C     SUBROUTINE VECFCN
C
C     THIS SUBROUTINE DEFINES FOURTEEN TEST FUNCTIONS. THE FIRST
C     FIVE TEST FUNCTIONS ARE OF DIMENSIONS 2,4,2,4,3, RESPECTIVELY,
C     WHILE THE REMAINING TEST FUNCTIONS ARE OF VARIABLE DIMENSION
C     N FOR ANY N GREATER THAN OR EQUAL TO 1 (PROBLEM 6 IS AN
C     EXCEPTION TO THIS, SINCE IT DOES NOT ALLOW N = 1).
C
C     THE SUBROUTINE STATEMENT IS
C
C       SUBROUTINE VECFCN(N,X,FVEC,NPROB)
C
C     WHERE
C
C       N IS A POSITIVE INTEGER VARIABLE.
C
C       X IS AN ARRAY OF LENGTH N.
C
C       FVEC IS AN OUTPUT ARRAY OF LENGTH N WHICH
C         CONTAINS THE NPROB FUNCTION VECTOR EVALUATED AT X.
C
C       NPROB IS A POSITIVE INTEGER VARIABLE WHICH DEFINES THE
C         NUMBER OF THE PROBLEM. NPROB MUST NOT EXCEED 14.
C
C     SUBPROGRAMS REQUIRED
C
C       FORTRAN-SUPPLIED ... ATAN,COS,EXP,SIGN,SIN,SQRT,
C                            MAX0,MIN0
C
C     MINPACK. VERSION OF JULY 1978.
C     BURTON S. GARBOW, KENNETH E. HILLSTROM, JORGE J. MORE
C
C     **********
      INTEGER I,IEV,IVAR,J,K,K1,K2,KP1,ML,MU
      REAL C1,C2,C3,C4,C5,C6,C7,C8,C9,EIGHT,FIVE,H,ONE,PROD,SUM,SUM1,
     *     SUM2,TEMP,TEMP1,TEMP2,TEN,THREE,TI,TJ,TK,TPI,TWO,ZERO
      REAL FLOAT
      DATA ZERO,ONE,TWO,THREE,FIVE,EIGHT,TEN
     *     /0.0E0,1.0E0,2.0E0,3.0E0,5.0E0,8.0E0,1.0E1/
      DATA C1,C2,C3,C4,C5,C6,C7,C8,C9
     *     /1.0E4,1.0001E0,2.0E2,2.02E1,1.98E1,1.8E2,2.5E-1,5.0E-1,
     *      2.9E1/
      FLOAT(IVAR) = IVAR
C
C     PROBLEM SELECTOR.
C
      GO TO (10,20,30,40,50,60,120,170,200,220,270,300,330,350), NPROB
C
C     ROSENBROCK FUNCTION.
C
   10 CONTINUE
      FVEC(1) = ONE - X(1)
      FVEC(2) = TEN*(X(2) - X(1)**2)
      GO TO 380
C
C     POWELL SINGULAR FUNCTION.
C
   20 CONTINUE
      FVEC(1) = X(1) + TEN*X(2)
      FVEC(2) = SQRT(FIVE)*(X(3) - X(4))
      FVEC(3) = (X(2) - TWO*X(3))**2
      FVEC(4) = SQRT(TEN)*(X(1) - X(4))**2
      GO TO 380
C
C     POWELL BADLY SCALED FUNCTION.
C
   30 CONTINUE
      FVEC(1) = C1*X(1)*X(2) - ONE
      FVEC(2) = EXP(-X(1)) + EXP(-X(2)) - C2
      GO TO 380
C
C     WOOD FUNCTION.
C
   40 CONTINUE
      TEMP1 = X(2) - X(1)**2
      TEMP2 = X(4) - X(3)**2
      FVEC(1) = -C3*X(1)*TEMP1 - (ONE - X(1))
      FVEC(2) = C3*TEMP1 + C4*(X(2) - ONE) + C5*(X(4) - ONE)
      FVEC(3) = -C6*X(3)*TEMP2 - (ONE - X(3))
      FVEC(4) = C6*TEMP2 + C4*(X(4) - ONE) + C5*(X(2) - ONE)
      GO TO 380
C
C     HELICAL VALLEY FUNCTION.
C
   50 CONTINUE
      TPI = EIGHT*ATAN(ONE)
      TEMP1 = SIGN(C7,X(2))
      IF (X(1) .GT. ZERO) TEMP1 = ATAN(X(2)/X(1))/TPI
      IF (X(1) .LT. ZERO) TEMP1 = ATAN(X(2)/X(1))/TPI + C8
      TEMP2 = SQRT(X(1)**2+X(2)**2)
      FVEC(1) = TEN*(X(3) - TEN*TEMP1)
      FVEC(2) = TEN*(TEMP2 - ONE)
      FVEC(3) = X(3)
      GO TO 380
C
C     WATSON FUNCTION.
C
   60 CONTINUE
      DO 70 K = 1, N
         FVEC(K) = ZERO
   70    CONTINUE
      DO 110 I = 1, 29
         TI = FLOAT(I)/C9
         SUM1 = ZERO
         TEMP = ONE
         DO 80 J = 2, N
            SUM1 = SUM1 + FLOAT(J-1)*TEMP*X(J)
            TEMP = TI*TEMP
   80       CONTINUE
         SUM2 = ZERO
         TEMP = ONE
         DO 90 J = 1, N
            SUM2 = SUM2 + TEMP*X(J)
            TEMP = TI*TEMP
   90       CONTINUE
         TEMP1 = SUM1 - SUM2**2 - ONE
         TEMP2 = TWO*TI*SUM2
         TEMP = ONE/TI
         DO 100 K = 1, N
            FVEC(K) = FVEC(K) + TEMP*(FLOAT(K-1) - TEMP2)*TEMP1
            TEMP = TI*TEMP
  100       CONTINUE
  110    CONTINUE
      TEMP = X(2) - X(1)**2 - ONE
      FVEC(1) = FVEC(1) + X(1)*(ONE - TWO*TEMP)
      FVEC(2) = FVEC(2) + TEMP
      GO TO 380
C
C     CHEBYQUAD FUNCTION.
C
  120 CONTINUE
      DO 130 K = 1, N
         FVEC(K) = ZERO
  130    CONTINUE
      DO 150 J = 1, N
         TEMP1 = ONE
         TEMP2 = TWO*X(J) - ONE
         TEMP = TWO*TEMP2
         DO 140 I = 1, N
            FVEC(I) = FVEC(I) + TEMP2
            TI = TEMP*TEMP2 - TEMP1
            TEMP1 = TEMP2
            TEMP2 = TI
  140       CONTINUE
  150    CONTINUE
      TK = ONE/FLOAT(N)
      IEV = -1
      DO 160 K = 1, N
         FVEC(K) = TK*FVEC(K)
         IF (IEV .GT. 0) FVEC(K) = FVEC(K) + ONE/(FLOAT(K)**2 - ONE)
         IEV = -IEV
  160    CONTINUE
      GO TO 380
C
C     BROWN ALMOST-LINEAR FUNCTION.
C
  170 CONTINUE
      SUM = -FLOAT(N+1)
      PROD = ONE
      DO 180 J = 1, N
         SUM = SUM + X(J)
         PROD = X(J)*PROD
  180    CONTINUE
      DO 190 K = 1, N
         FVEC(K) = X(K) + SUM
  190    CONTINUE
      FVEC(N) = PROD - ONE
      GO TO 380
C
C     DISCRETE BOUNDARY VALUE FUNCTION.
C
  200 CONTINUE
      H = ONE/FLOAT(N+1)
      DO 210 K = 1, N
         TEMP = (X(K) + FLOAT(K)*H + ONE)**3
         TEMP1 = ZERO
         IF (K .NE. 1) TEMP1 = X(K-1)
         TEMP2 = ZERO
         IF (K .NE. N) TEMP2 = X(K+1)
         FVEC(K) = TWO*X(K) - TEMP1 - TEMP2 + TEMP*H**2/TWO
  210    CONTINUE
      GO TO 380
C
C     DISCRETE INTEGRAL EQUATION FUNCTION.
C
  220 CONTINUE
      H = ONE/FLOAT(N+1)
      DO 260 K = 1, N
         TK = FLOAT(K)*H
         SUM1 = ZERO
         DO 230 J = 1, K
            TJ = FLOAT(J)*H
            TEMP = (X(J) + TJ + ONE)**3
            SUM1 = SUM1 + TJ*TEMP
  230       CONTINUE
         SUM2 = ZERO
         KP1 = K + 1
         IF (N .LT. KP1) GO TO 250
         DO 240 J = KP1, N
            TJ = FLOAT(J)*H
            TEMP = (X(J) + TJ + ONE)**3
            SUM2 = SUM2 + (ONE - TJ)*TEMP
  240       CONTINUE
  250    CONTINUE
         FVEC(K) = X(K) + H*((ONE - TK)*SUM1 + TK*SUM2)/TWO
  260    CONTINUE
      GO TO 380
C
C     TRIGONOMETRIC FUNCTION.
C
  270 CONTINUE
      SUM = ZERO
      DO 280 J = 1, N
         FVEC(J) = COS(X(J))
         SUM = SUM + FVEC(J)
  280    CONTINUE
      DO 290 K = 1, N
         FVEC(K) = FLOAT(N+K) - SIN(X(K)) - SUM - FLOAT(K)*FVEC(K)
  290    CONTINUE
      GO TO 380
C
C     VARIABLY DIMENSIONED FUNCTION.
C
  300 CONTINUE
      SUM = ZERO
      DO 310 J = 1, N
         SUM = SUM + FLOAT(J)*(X(J) - ONE)
  310    CONTINUE
      TEMP = SUM*(ONE + TWO*SUM**2)
      DO 320 K = 1, N
         FVEC(K) = X(K) - ONE + FLOAT(K)*TEMP
  320    CONTINUE
      GO TO 380
C
C     BROYDEN TRIDIAGONAL FUNCTION.
C
  330 CONTINUE
      DO 340 K = 1, N
         TEMP = (THREE - TWO*X(K))*X(K)
         TEMP1 = ZERO
         IF (K .NE. 1) TEMP1 = X(K-1)
         TEMP2 = ZERO
         IF (K .NE. N) TEMP2 = X(K+1)
         FVEC(K) = TEMP - TEMP1 - TWO*TEMP2 + ONE
  340    CONTINUE
      GO TO 380
C
C     BROYDEN BANDED FUNCTION.
C
  350 CONTINUE
      ML = 5
      MU = 1
      DO 370 K = 1, N
         K1 = MAX0(1,K-ML)
         K2 = MIN0(K+MU,N)
         TEMP = ZERO
         DO 360 J = K1, K2
            IF (J .NE. K) TEMP = TEMP + X(J)*(ONE + X(J))
  360       CONTINUE
         FVEC(K) = X(K)*(TWO + FIVE*X(K)**2) + ONE - TEMP
  370    CONTINUE
  380 CONTINUE
      RETURN
C
C     LAST CARD OF SUBROUTINE VECFCN.
C
      END
      SUBROUTINE VECJAC(N,X,FJAC,LDFJAC,NPROB)                          00000010
      INTEGER N,LDFJAC,NPROB
      REAL X(N),FJAC(LDFJAC,N)
C     **********
C
C     SUBROUTINE VECJAC
C
C     THIS SUBROUTINE DEFINES THE JACOBIAN MATRICES OF FOURTEEN
C     TEST FUNCTIONS. THE PROBLEM DIMENSIONS ARE AS DESCRIBED
C     IN THE PROLOGUE COMMENTS OF VECFCN.
C
C     THE SUBROUTINE STATEMENT IS
C
C       SUBROUTINE VECJAC(N,X,FJAC,LDFJAC,NPROB)
C
C     WHERE
C
C       N IS A POSITIVE INTEGER VARIABLE.
C
C       X IS A LINEAR ARRAY OF LENGTH N.
C
C       FJAC IS AN N BY N ARRAY. ON OUTPUT FJAC CONTAINS THE
C         JACOBIAN MATRIX OF THE NPROB FUNCTION EVALUATED AT X.
C
C       LDFJAC IS A POSITIVE INTEGER VARIABLE NOT LESS THAN N
C         WHICH SPECIFIES THE LEADING DIMENSION OF THE ARRAY FJAC.
C
C       NPROB IS A POSITIVE INTEGER VARIABLE WHICH DEFINES THE
C         NUMBER OF THE PROBLEM. NPROB MUST NOT EXCEED 14.
C
C     SUBPROGRAMS REQUIRED
C
C       FORTRAN-SUPPLIED ... ATAN,COS,EXP,AMIN1,SIN,SQRT,
C                            MAX0,MIN0
C
C     MINPACK. VERSION OF AUGUST 1978.
C     BURTON S. GARBOW, KENNETH E. HILLSTROM, JORGE J. MORE
C
C     **********
      INTEGER I,IVAR,J,K,K1,K2,ML,MU
      REAL C1,C3,C4,C5,C6,C9,EIGHT,FIFTN,FIVE,FOUR,H,HUNDRD,ONE,PROD,
     *     SIX,SUM,SUM1,SUM2,TEMP,TEMP1,TEMP2,TEMP3,TEMP4,TEN,THREE,
     *     TI,TJ,TK,TPI,TWENTY,TWO,ZERO
      REAL FLOAT
      DATA ZERO,ONE,TWO,THREE,FOUR,FIVE,SIX,EIGHT,TEN,FIFTN,TWENTY,
     *     HUNDRD
     *     /0.0E0,1.0E0,2.0E0,3.0E0,4.0E0,5.0E0,6.0E0,8.0E0,1.0E1,
     *      1.5E1,2.0E1,1.0E2/
      DATA C1,C3,C4,C5,C6,C9 /1.0E4,2.0E2,2.02E1,1.98E1,1.8E2,2.9E1/
      FLOAT(IVAR) = IVAR
C
C     JACOBIAN ROUTINE SELECTOR.
C
      GO TO (10,20,50,60,90,100,200,230,290,320,350,380,420,450),
     *      NPROB
C
C     ROSENBROCK FUNCTION.
C
   10 CONTINUE
      FJAC(1,1) = -ONE
      FJAC(1,2) = ZERO
      FJAC(2,1) = -TWENTY*X(1)
      FJAC(2,2) = TEN
      GO TO 490
C
C     POWELL SINGULAR FUNCTION.
C
   20 CONTINUE
      DO 40 K = 1, 4
         DO 30 J = 1, 4
            FJAC(K,J) = ZERO
   30       CONTINUE
   40    CONTINUE
      FJAC(1,1) = ONE
      FJAC(1,2) = TEN
      FJAC(2,3) = SQRT(FIVE)
      FJAC(2,4) = -FJAC(2,3)
      FJAC(3,2) = TWO*(X(2) - TWO*X(3))
      FJAC(3,3) = -TWO*FJAC(3,2)
      FJAC(4,1) = TWO*SQRT(TEN)*(X(1) - X(4))
      FJAC(4,4) = -FJAC(4,1)
      GO TO 490
C
C     POWELL BADLY SCALED FUNCTION.
C
   50 CONTINUE
      FJAC(1,1) = C1*X(2)
      FJAC(1,2) = C1*X(1)
      FJAC(2,1) = -EXP(-X(1))
      FJAC(2,2) = -EXP(-X(2))
      GO TO 490
C
C     WOOD FUNCTION.
C
   60 CONTINUE
      DO 80 K = 1, 4
         DO 70 J = 1, 4
            FJAC(K,J) = ZERO
   70       CONTINUE
   80    CONTINUE
      TEMP1 = X(2) - THREE*X(1)**2
      TEMP2 = X(4) - THREE*X(3)**2
      FJAC(1,1) = -C3*TEMP1 + ONE
      FJAC(1,2) = -C3*X(1)
      FJAC(2,1) = -TWO*C3*X(1)
      FJAC(2,2) = C3 + C4
      FJAC(2,4) = C5
      FJAC(3,3) = -C6*TEMP2 + ONE
      FJAC(3,4) = -C6*X(3)
      FJAC(4,2) = C5
      FJAC(4,3) = -TWO*C6*X(3)
      FJAC(4,4) = C6 + C4
      GO TO 490
C
C     HELICAL VALLEY FUNCTION.
C
   90 CONTINUE
      TPI = EIGHT*ATAN(ONE)
      TEMP = X(1)**2 + X(2)**2
      TEMP1 = TPI*TEMP
      TEMP2 = SQRT(TEMP)
      FJAC(1,1) = HUNDRD*X(2)/TEMP1
      FJAC(1,2) = -HUNDRD*X(1)/TEMP1
      FJAC(1,3) = TEN
      FJAC(2,1) = TEN*X(1)/TEMP2
      FJAC(2,2) = TEN*X(2)/TEMP2
      FJAC(2,3) = ZERO
      FJAC(3,1) = ZERO
      FJAC(3,2) = ZERO
      FJAC(3,3) = ONE
      GO TO 490
C
C     WATSON FUNCTION.
C
  100 CONTINUE
      DO 120 K = 1, N
         DO 110 J = K, N
            FJAC(K,J) = ZERO
  110       CONTINUE
  120    CONTINUE
      DO 170 I = 1, 29
         TI = FLOAT(I)/C9
         SUM1 = ZERO
         TEMP = ONE
         DO 130 J = 2, N
            SUM1 = SUM1 + FLOAT(J-1)*TEMP*X(J)
            TEMP = TI*TEMP
  130       CONTINUE
         SUM2 = ZERO
         TEMP = ONE
         DO 140 J = 1, N
            SUM2 = SUM2 + TEMP*X(J)
            TEMP = TI*TEMP
  140       CONTINUE
         TEMP1 = TWO*(SUM1 - SUM2**2 - ONE)
         TEMP2 = TWO*SUM2
         TEMP = TI**2
         TK = ONE
         DO 160 K = 1, N
            TJ = TK
            DO 150 J = K, N
               FJAC(K,J) = FJAC(K,J)
     *                     + TJ
     *                       *((FLOAT(K-1)/TI - TEMP2)
     *                         *(FLOAT(J-1)/TI - TEMP2) - TEMP1)
               TJ = TI*TJ
  150          CONTINUE
            TK = TEMP*TK
  160       CONTINUE
  170    CONTINUE
      FJAC(1,1) = FJAC(1,1) + SIX*X(1)**2 - TWO*X(2) + THREE
      FJAC(1,2) = FJAC(1,2) - TWO*X(1)
      FJAC(2,2) = FJAC(2,2) + ONE
      DO 190 K = 1, N
         DO 180 J = K, N
            FJAC(J,K) = FJAC(K,J)
  180       CONTINUE
  190    CONTINUE
      GO TO 490
C
C     CHEBYQUAD FUNCTION.
C
  200 CONTINUE
      TK = ONE/FLOAT(N)
      DO 220 J = 1, N
         TEMP1 = ONE
         TEMP2 = TWO*X(J) - ONE
         TEMP = TWO*TEMP2
         TEMP3 = ZERO
         TEMP4 = TWO
         DO 210 K = 1, N
            FJAC(K,J) = TK*TEMP4
            TI = FOUR*TEMP2 + TEMP*TEMP4 - TEMP3
            TEMP3 = TEMP4
            TEMP4 = TI
            TI = TEMP*TEMP2 - TEMP1
            TEMP1 = TEMP2
            TEMP2 = TI
  210       CONTINUE
  220    CONTINUE
      GO TO 490
C
C     BROWN ALMOST-LINEAR FUNCTION.
C
  230 CONTINUE
      PROD = ONE
      DO 250 J = 1, N
         PROD = X(J)*PROD
         DO 240 K = 1, N
            FJAC(K,J) = ONE
  240       CONTINUE
         FJAC(J,J) = TWO
  250    CONTINUE
      DO 280 J = 1, N
         TEMP = X(J)
         IF (TEMP .NE. ZERO) GO TO 270
         TEMP = ONE
         PROD = ONE
         DO 260 K = 1, N
            IF (K .NE. J) PROD = X(K)*PROD
  260       CONTINUE
  270    CONTINUE
         FJAC(N,J) = PROD/TEMP
  280    CONTINUE
      GO TO 490
C
C     DISCRETE BOUNDARY VALUE FUNCTION.
C
  290 CONTINUE
      H = ONE/FLOAT(N+1)
      DO 310 K = 1, N
         TEMP = THREE*(X(K) + FLOAT(K)*H + ONE)**2
         DO 300 J = 1, N
            FJAC(K,J) = ZERO
  300       CONTINUE
         FJAC(K,K) = TWO + TEMP*H**2/TWO
         IF (K .NE. 1) FJAC(K,K-1) = -ONE
         IF (K .NE. N) FJAC(K,K+1) = -ONE
  310    CONTINUE
      GO TO 490
C
C     DISCRETE INTEGRAL EQUATION FUNCTION.
C
  320 CONTINUE
      H = ONE/FLOAT(N+1)
      DO 340 K = 1, N
         TK = FLOAT(K)*H
         DO 330 J = 1, N
            TJ = FLOAT(J)*H
            TEMP = THREE*(X(J) + TJ + ONE)**2
            FJAC(K,J) = H*AMIN1(TJ*(ONE-TK),TK*(ONE-TJ))*TEMP/TWO
  330       CONTINUE
         FJAC(K,K) = FJAC(K,K) + ONE
  340    CONTINUE
      GO TO 490
C
C     TRIGONOMETRIC FUNCTION.
C
  350 CONTINUE
      DO 370 J = 1, N
         TEMP = SIN(X(J))
         DO 360 K = 1, N
            FJAC(K,J) = TEMP
  360       CONTINUE
         FJAC(J,J) = FLOAT(J+1)*TEMP - COS(X(J))
  370    CONTINUE
      GO TO 490
C
C     VARIABLY DIMENSIONED FUNCTION.
C
  380 CONTINUE
      SUM = ZERO
      DO 390 J = 1, N
         SUM = SUM + FLOAT(J)*(X(J) - ONE)
  390    CONTINUE
      TEMP = ONE + SIX*SUM**2
      DO 410 K = 1, N
         DO 400 J = K, N
            FJAC(K,J) = FLOAT(K*J)*TEMP
            FJAC(J,K) = FJAC(K,J)
  400       CONTINUE
         FJAC(K,K) = FJAC(K,K) + ONE
  410    CONTINUE
      GO TO 490
C
C     BROYDEN TRIDIAGONAL FUNCTION.
C
  420 CONTINUE
      DO 440 K = 1, N
         DO 430 J = 1, N
            FJAC(K,J) = ZERO
  430       CONTINUE
         FJAC(K,K) = THREE - FOUR*X(K)
         IF (K .NE. 1) FJAC(K,K-1) = -ONE
         IF (K .NE. N) FJAC(K,K+1) = -TWO
  440    CONTINUE
      GO TO 490
C
C     BROYDEN BANDED FUNCTION.
C
  450 CONTINUE
      ML = 5
      MU = 1
      DO 480 K = 1, N
         DO 460 J = 1, N
            FJAC(K,J) = ZERO
  460       CONTINUE
         K1 = MAX0(1,K-ML)
         K2 = MIN0(K+MU,N)
         DO 470 J = K1, K2
            IF (J .NE. K) FJAC(K,J) = -(ONE + TWO*X(J))
  470       CONTINUE
         FJAC(K,K) = TWO + FIFTN*X(K)**2
  480    CONTINUE
  490 CONTINUE
      RETURN
C
C     LAST CARD OF SUBROUTINE VECJAC.
C
      END
C ===== 6. SINGLE PRECISION TESTING AIDS FOR NONLINEAR EQUATIONS.
      SUBROUTINE INITPT(N,X,NPROB,FACTOR)                               00000010
      INTEGER N,NPROB
      REAL FACTOR
      REAL X(N)
C     **********
C
C     SUBROUTINE INITPT
C
C     THIS SUBROUTINE SPECIFIES THE STANDARD STARTING POINTS FOR THE
C     FUNCTIONS DEFINED BY SUBROUTINE SSQFCN. THE SUBROUTINE RETURNS
C     IN X A MULTIPLE (FACTOR) OF THE STANDARD STARTING POINT. FOR
C     THE 11TH FUNCTION THE STANDARD STARTING POINT IS ZERO, SO IN
C     THIS CASE, IF FACTOR IS NOT UNITY, THEN THE SUBROUTINE RETURNS
C     THE VECTOR  X(J) = FACTOR, J=1,...,N.
C
C     THE SUBROUTINE STATEMENT IS
C
C       SUBROUTINE INITPT(N,X,NPROB,FACTOR)
C
C     WHERE
C
C       N IS A POSITIVE INTEGER VARIABLE.
C
C       X IS AN OUTPUT ARRAY OF LENGTH N THAT CONTAINS THE
C         STANDARD STARTING POINT FOR PROBLEM NPROB MULTIPLIED BY
C         FACTOR.
C
C       NPROB IS A POSITIVE INTEGER VARIABLE WHICH DEFINES THE
C         NUMBER OF THE PROBLEM. NPROB MUST NOT EXCEED 18.
C
C       FACTOR SPECIFIES THE MULTIPLE OF THE STANDARD STARTING
C         POINT. IF FACTOR IS UNITY, NO MULTIPLICATION IS PERFORMED.
C
C     MINPACK. VERSION OF JULY 1978.
C     BURTON S. GARBOW, KENNETH E. HILLSTROM, JORGE J. MORE
C
C     **********
      INTEGER IVAR,J
      REAL C1,C2,C3,C4,C5,C6,C7,C8,C9,C10,C11,C12,C13,C14,C15,C16,C17,
     *     FIVE,H,HALF,ONE,SEVEN,TEN,THREE,TWENTY,TWNTF,TWO,ZERO
      REAL FLOAT
      DATA ZERO,HALF,ONE,TWO,THREE,FIVE,SEVEN,TEN,TWENTY,TWNTF
     *     /0.0E0,5.0E-1,1.0E0,2.0E0,3.0E0,5.0E0,7.0E0,1.0E1,2.0E1,
     *      2.5E1/
      DATA C1,C2,C3,C4,C5,C6,C7,C8,C9,C10,C11,C12,C13,C14,C15,C16,C17
     *     /1.2E0,2.5E-1,3.9E-1,4.15E-1,2.0E-2,4.0E3,2.5E2,3.0E-1,
     *      4.0E-1,1.5E0,1.0E-2,1.3E0,6.5E-1,7.0E-1,6.0E-1,4.5E0,
     *      5.5E0/
      FLOAT(IVAR) = IVAR
C
C     SELECTION OF INITIAL POINT.
C
      GO TO (10,10,10,30,40,50,60,70,80,90,100,120,130,140,150,170,
     *       190,200), NPROB
C
C     LINEAR FUNCTION - FULL RANK OR RANK 1.
C
   10 CONTINUE
      DO 20 J = 1, N
         X(J) = ONE
   20    CONTINUE
      GO TO 210
C
C     ROSENBROCK FUNCTION.
C
   30 CONTINUE
      X(1) = -C1
      X(2) = ONE
      GO TO 210
C
C     HELICAL VALLEY FUNCTION.
C
   40 CONTINUE
      X(1) = -ONE
      X(2) = ZERO
      X(3) = ZERO
      GO TO 210
C
C     POWELL SINGULAR FUNCTION.
C
   50 CONTINUE
      X(1) = THREE
      X(2) = -ONE
      X(3) = ZERO
      X(4) = ONE
      GO TO 210
C
C     FREUDENSTEIN AND ROTH FUNCTION.
C
   60 CONTINUE
      X(1) = HALF
      X(2) = -TWO
      GO TO 210
C
C     BARD FUNCTION.
C
   70 CONTINUE
      X(1) = ONE
      X(2) = ONE
      X(3) = ONE
      GO TO 210
C
C     KOWALIK AND OSBORNE FUNCTION.
C
   80 CONTINUE
      X(1) = C2
      X(2) = C3
      X(3) = C4
      X(4) = C3
      GO TO 210
C
C     MEYER FUNCTION.
C
   90 CONTINUE
      X(1) = C5
      X(2) = C6
      X(3) = C7
      GO TO 210
C
C     WATSON FUNCTION.
C
  100 CONTINUE
      DO 110 J = 1, N
         X(J) = ZERO
  110    CONTINUE
      GO TO 210
C
C     BOX 3-DIMENSIONAL FUNCTION.
C
  120 CONTINUE
      X(1) = ZERO
      X(2) = TEN
      X(3) = TWENTY
      GO TO 210
C
C     JENNRICH AND SAMPSON FUNCTION.
C
  130 CONTINUE
      X(1) = C8
      X(2) = C9
      GO TO 210
C
C     BROWN AND DENNIS FUNCTION.
C
  140 CONTINUE
      X(1) = TWNTF
      X(2) = FIVE
      X(3) = -FIVE
      X(4) = -ONE
      GO TO 210
C
C     CHEBYQUAD FUNCTION.
C
  150 CONTINUE
      H = ONE/FLOAT(N+1)
      DO 160 J = 1, N
         X(J) = FLOAT(J)*H
  160    CONTINUE
      GO TO 210
C
C     BROWN ALMOST-LINEAR FUNCTION.
C
  170 CONTINUE
      DO 180 J = 1, N
         X(J) = HALF
  180    CONTINUE
      GO TO 210
C
C     OSBORNE 1 FUNCTION.
C
  190 CONTINUE
      X(1) = HALF
      X(2) = C10
      X(3) = -ONE
      X(4) = C11
      X(5) = C5
      GO TO 210
C
C     OSBORNE 2 FUNCTION.
C
  200 CONTINUE
      X(1) = C12
      X(2) = C13
      X(3) = C13
      X(4) = C14
      X(5) = C15
      X(6) = THREE
      X(7) = FIVE
      X(8) = SEVEN
      X(9) = TWO
      X(10) = C16
      X(11) = C17
  210 CONTINUE
C
C     COMPUTE MULTIPLE OF INITIAL POINT.
C
      IF (FACTOR .EQ. ONE) GO TO 260
      IF (NPROB .EQ. 11) GO TO 230
         DO 220 J = 1, N
            X(J) = FACTOR*X(J)
  220       CONTINUE
         GO TO 250
  230 CONTINUE
         DO 240 J = 1, N
            X(J) = FACTOR
  240       CONTINUE
  250 CONTINUE
  260 CONTINUE
      RETURN
C
C     LAST CARD OF SUBROUTINE INITPT.
C
      END
      SUBROUTINE SSQFCN(M,N,X,FVEC,NPROB)                               00000010
      INTEGER M,N,NPROB
      REAL X(N),FVEC(M)
C     **********
C
C     SUBROUTINE SSQFCN
C
C     THIS SUBROUTINE DEFINES THE FUNCTIONS OF EIGHTEEN NONLINEAR
C     LEAST SQUARES PROBLEMS. THE ALLOWABLE VALUES OF (M,N) FOR
C     FUNCTIONS 1,2 AND 3 ARE VARIABLE BUT WITH M .GE. N.
C     FOR FUNCTIONS 4,5,6,7,8,9 AND 10 THE VALUES OF (M,N) ARE
C     (2,2),(3,3),(4,4),(2,2),(15,3),(11,4) AND (16,3), RESPECTIVELY.
C     FUNCTION 11 (WATSON) HAS M = 31 WITH N USUALLY 6 OR 9.
C     HOWEVER, ANY N, N = 2,...,31, IS PERMITTED.
C     FUNCTIONS 12,13 AND 14 HAVE N = 3,2 AND 4, RESPECTIVELY, BUT
C     ALLOW ANY M .GE. N, WITH THE USUAL CHOICES BEING 10,10 AND 20.
C     FUNCTION 15 (CHEBYQUAD) ALLOWS M AND N VARIABLE WITH M .GE. N.
C     FUNCTION 16 (BROWN) ALLOWS N VARIABLE WITH M = N.
C     FOR FUNCTIONS 17 AND 18, THE VALUES OF (M,N) ARE
C     (33,5) AND (65,11), RESPECTIVELY.
C
C     THE SUBROUTINE STATEMENT IS
C
C       SUBROUTINE SSQFCN(M,N,X,FVEC,NPROB)
C
C     WHERE
C
C       M AND N ARE POSITIVE INTEGER VARIABLES. N MUST NOT EXCEED M.
C
C       X IS AN ARRAY OF LENGTH N.
C
C       FVEC IS AN OUTPUT ARRAY OF LENGTH M WHICH
C         CONTAINS THE NPROB FUNCTION EVALUATED AT X.
C
C       NPROB IS A POSITIVE INTEGER VARIABLE WHICH DEFINES THE
C         NUMBER OF THE PROBLEM. NPROB MUST NOT EXCEED 18.
C
C     SUBPROGRAMS REQUIRED
C
C       FORTRAN-SUPPLIED ... ATAN,COS,EXP,SIN,SQRT,SIGN
C
C     MINPACK. VERSION OF JULY 1978.
C     BURTON S. GARBOW, KENNETH E. HILLSTROM, JORGE J. MORE
C
C     **********
      INTEGER I,IEV,IVAR,J,NM1
      REAL C13,C14,C29,C45,DIV,DX,EIGHT,FIVE,ONE,PROD,SUM,S1,S2,TEMP,
     *     TEN,TI,TMP1,TMP2,TMP3,TMP4,TPI,TWO,ZERO,ZP25,ZP5
      REAL V(11),Y1(15),Y2(11),Y3(16),Y4(33),Y5(65)
      REAL FLOAT
      DATA ZERO,ZP25,ZP5,ONE,TWO,FIVE,EIGHT,TEN,C13,C14,C29,C45
     *     /0.0E0,2.5E-1,5.0E-1,1.0E0,2.0E0,5.0E0,8.0E0,1.0E1,1.3E1,
     *      1.4E1,2.9E1,4.5E1/
      DATA V(1),V(2),V(3),V(4),V(5),V(6),V(7),V(8),V(9),V(10),V(11)
     *     /4.0E0,2.0E0,1.0E0,5.0E-1,2.5E-1,1.67E-1,1.25E-1,1.0E-1,
     *      8.33E-2,7.14E-2,6.25E-2/
      DATA Y1(1),Y1(2),Y1(3),Y1(4),Y1(5),Y1(6),Y1(7),Y1(8),Y1(9),
     *     Y1(10),Y1(11),Y1(12),Y1(13),Y1(14),Y1(15)
     *     /1.4E-1,1.8E-1,2.2E-1,2.5E-1,2.9E-1,3.2E-1,3.5E-1,3.9E-1,
     *      3.7E-1,5.8E-1,7.3E-1,9.6E-1,1.34E0,2.1E0,4.39E0/
      DATA Y2(1),Y2(2),Y2(3),Y2(4),Y2(5),Y2(6),Y2(7),Y2(8),Y2(9),
     *     Y2(10),Y2(11)
     *     /1.957E-1,1.947E-1,1.735E-1,1.6E-1,8.44E-2,6.27E-2,4.56E-2,
     *      3.42E-2,3.23E-2,2.35E-2,2.46E-2/
      DATA Y3(1),Y3(2),Y3(3),Y3(4),Y3(5),Y3(6),Y3(7),Y3(8),Y3(9),
     *     Y3(10),Y3(11),Y3(12),Y3(13),Y3(14),Y3(15),Y3(16)
     *     /3.478E4,2.861E4,2.365E4,1.963E4,1.637E4,1.372E4,1.154E4,
     *      9.744E3,8.261E3,7.03E3,6.005E3,5.147E3,4.427E3,3.82E3,
     *      3.307E3,2.872E3/
      DATA Y4(1),Y4(2),Y4(3),Y4(4),Y4(5),Y4(6),Y4(7),Y4(8),Y4(9),
     *     Y4(10),Y4(11),Y4(12),Y4(13),Y4(14),Y4(15),Y4(16),Y4(17),
     *     Y4(18),Y4(19),Y4(20),Y4(21),Y4(22),Y4(23),Y4(24),Y4(25),
     *     Y4(26),Y4(27),Y4(28),Y4(29),Y4(30),Y4(31),Y4(32),Y4(33)
     *     /8.44E-1,9.08E-1,9.32E-1,9.36E-1,9.25E-1,9.08E-1,8.81E-1,
     *      8.5E-1,8.18E-1,7.84E-1,7.51E-1,7.18E-1,6.85E-1,6.58E-1,
     *      6.28E-1,6.03E-1,5.8E-1,5.58E-1,5.38E-1,5.22E-1,5.06E-1,
     *      4.9E-1,4.78E-1,4.67E-1,4.57E-1,4.48E-1,4.38E-1,4.31E-1,
     *      4.24E-1,4.2E-1,4.14E-1,4.11E-1,4.06E-1/
      DATA Y5(1),Y5(2),Y5(3),Y5(4),Y5(5),Y5(6),Y5(7),Y5(8),Y5(9),
     *     Y5(10),Y5(11),Y5(12),Y5(13),Y5(14),Y5(15),Y5(16),Y5(17),
     *     Y5(18),Y5(19),Y5(20),Y5(21),Y5(22),Y5(23),Y5(24),Y5(25),
     *     Y5(26),Y5(27),Y5(28),Y5(29),Y5(30),Y5(31),Y5(32),Y5(33),
     *     Y5(34),Y5(35),Y5(36),Y5(37),Y5(38),Y5(39),Y5(40),Y5(41),
     *     Y5(42),Y5(43),Y5(44),Y5(45),Y5(46),Y5(47),Y5(48),Y5(49),
     *     Y5(50),Y5(51),Y5(52),Y5(53),Y5(54),Y5(55),Y5(56),Y5(57),
     *     Y5(58),Y5(59),Y5(60),Y5(61),Y5(62),Y5(63),Y5(64),Y5(65)
     *     /1.366E0,1.191E0,1.112E0,1.013E0,9.91E-1,8.85E-1,8.31E-1,
     *      8.47E-1,7.86E-1,7.25E-1,7.46E-1,6.79E-1,6.08E-1,6.55E-1,
     *      6.16E-1,6.06E-1,6.02E-1,6.26E-1,6.51E-1,7.24E-1,6.49E-1,
     *      6.49E-1,6.94E-1,6.44E-1,6.24E-1,6.61E-1,6.12E-1,5.58E-1,
     *      5.33E-1,4.95E-1,5.0E-1,4.23E-1,3.95E-1,3.75E-1,3.72E-1,
     *      3.91E-1,3.96E-1,4.05E-1,4.28E-1,4.29E-1,5.23E-1,5.62E-1,
     *      6.07E-1,6.53E-1,6.72E-1,7.08E-1,6.33E-1,6.68E-1,6.45E-1,
     *      6.32E-1,5.91E-1,5.59E-1,5.97E-1,6.25E-1,7.39E-1,7.1E-1,
     *      7.29E-1,7.2E-1,6.36E-1,5.81E-1,4.28E-1,2.92E-1,1.62E-1,
     *      9.8E-2,5.4E-2/
      FLOAT(IVAR) = IVAR
C
C     FUNCTION ROUTINE SELECTOR.
C
      GO TO (10,40,70,110,120,130,140,150,170,190,210,250,270,290,310,
     *       360,390,410), NPROB
C
C     LINEAR FUNCTION - FULL RANK.
C
   10 CONTINUE
      SUM = ZERO
      DO 20 J = 1, N
         SUM = SUM + X(J)
   20    CONTINUE
      TEMP = TWO*SUM/FLOAT(M) + ONE
      DO 30 I = 1, M
         FVEC(I) = -TEMP
         IF (I .LE. N) FVEC(I) = FVEC(I) + X(I)
   30    CONTINUE
      GO TO 430
C
C     LINEAR FUNCTION - RANK 1.
C
   40 CONTINUE
      SUM = ZERO
      DO 50 J = 1, N
         SUM = SUM + FLOAT(J)*X(J)
   50    CONTINUE
      DO 60 I = 1, M
         FVEC(I) = FLOAT(I)*SUM - ONE
   60    CONTINUE
      GO TO 430
C
C     LINEAR FUNCTION - RANK 1 WITH ZERO COLUMNS AND ROWS.
C
   70 CONTINUE
      SUM = ZERO
      NM1 = N - 1
      IF (NM1 .LT. 2) GO TO 90
      DO 80 J = 2, NM1
         SUM = SUM + FLOAT(J)*X(J)
   80    CONTINUE
   90 CONTINUE
      DO 100 I = 1, M
         FVEC(I) = FLOAT(I-1)*SUM - ONE
  100    CONTINUE
      FVEC(M) = -ONE
      GO TO 430
C
C     ROSENBROCK FUNCTION.
C
  110 CONTINUE
      FVEC(1) = TEN*(X(2) - X(1)**2)
      FVEC(2) = ONE - X(1)
      GO TO 430
C
C     HELICAL VALLEY FUNCTION.
C
  120 CONTINUE
      TPI = EIGHT*ATAN(ONE)
      TMP1 = SIGN(ZP25,X(2))
      IF (X(1) .GT. ZERO) TMP1 = ATAN(X(2)/X(1))/TPI
      IF (X(1) .LT. ZERO) TMP1 = ATAN(X(2)/X(1))/TPI + ZP5
      TMP2 = SQRT(X(1)**2+X(2)**2)
      FVEC(1) = TEN*(X(3) - TEN*TMP1)
      FVEC(2) = TEN*(TMP2 - ONE)
      FVEC(3) = X(3)
      GO TO 430
C
C     POWELL SINGULAR FUNCTION.
C
  130 CONTINUE
      FVEC(1) = X(1) + TEN*X(2)
      FVEC(2) = SQRT(FIVE)*(X(3) - X(4))
      FVEC(3) = (X(2) - TWO*X(3))**2
      FVEC(4) = SQRT(TEN)*(X(1) - X(4))**2
      GO TO 430
C
C     FREUDENSTEIN AND ROTH FUNCTION.
C
  140 CONTINUE
      FVEC(1) = -C13 + X(1) + ((FIVE - X(2))*X(2) - TWO)*X(2)
      FVEC(2) = -C29 + X(1) + ((ONE + X(2))*X(2) - C14)*X(2)
      GO TO 430
C
C     BARD FUNCTION.
C
  150 CONTINUE
      DO 160 I = 1, 15
         TMP1 = FLOAT(I)
         TMP2 = FLOAT(16-I)
         TMP3 = TMP1
         IF (I .GT. 8) TMP3 = TMP2
         FVEC(I) = Y1(I) - (X(1) + TMP1/(X(2)*TMP2 + X(3)*TMP3))
  160    CONTINUE
      GO TO 430
C
C     KOWALIK AND OSBORNE FUNCTION.
C
  170 CONTINUE
      DO 180 I = 1, 11
         TMP1 = V(I)*(V(I) + X(2))
         TMP2 = V(I)*(V(I) + X(3)) + X(4)
         FVEC(I) = Y2(I) - X(1)*TMP1/TMP2
  180    CONTINUE
      GO TO 430
C
C     MEYER FUNCTION.
C
  190 CONTINUE
      DO 200 I = 1, 16
         TEMP = FIVE*FLOAT(I) + C45 + X(3)
         TMP1 = X(2)/TEMP
         TMP2 = EXP(TMP1)
         FVEC(I) = X(1)*TMP2 - Y3(I)
  200    CONTINUE
      GO TO 430
C
C     WATSON FUNCTION.
C
  210 CONTINUE
      DO 240 I = 1, 29
         DIV = FLOAT(I)/C29
         S1 = ZERO
         DX = ONE
         DO 220 J = 2, N
            S1 = S1 + FLOAT(J-1)*DX*X(J)
            DX = DIV*DX
  220       CONTINUE
         S2 = ZERO
         DX = ONE
         DO 230 J = 1, N
            S2 = S2 + DX*X(J)
            DX = DIV*DX
  230       CONTINUE
         FVEC(I) = S1 - S2**2 - ONE
  240    CONTINUE
      FVEC(30) = X(1)
      FVEC(31) = X(2) - X(1)**2 - ONE
      GO TO 430
C
C     BOX 3-DIMENSIONAL FUNCTION.
C
  250 CONTINUE
      DO 260 I = 1, M
         TEMP = FLOAT(I)
         TMP1 = TEMP/TEN
         FVEC(I) = EXP(-TMP1*X(1)) - EXP(-TMP1*X(2))
     *             + (EXP(-TEMP) - EXP(-TMP1))*X(3)
  260    CONTINUE
      GO TO 430
C
C     JENNRICH AND SAMPSON FUNCTION.
C
  270 CONTINUE
      DO 280 I = 1, M
         TEMP = FLOAT(I)
         FVEC(I) = TWO + TWO*TEMP - EXP(TEMP*X(1)) - EXP(TEMP*X(2))
  280    CONTINUE
      GO TO 430
C
C     BROWN AND DENNIS FUNCTION.
C
  290 CONTINUE
      DO 300 I = 1, M
         TEMP = FLOAT(I)/FIVE
         TMP1 = X(1) + TEMP*X(2) - EXP(TEMP)
         TMP2 = X(3) + SIN(TEMP)*X(4) - COS(TEMP)
         FVEC(I) = TMP1**2 + TMP2**2
  300    CONTINUE
      GO TO 430
C
C     CHEBYQUAD FUNCTION.
C
  310 CONTINUE
      DO 320 I = 1, M
         FVEC(I) = ZERO
  320    CONTINUE
      DO 340 J = 1, N
         TMP1 = ONE
         TMP2 = TWO*X(J) - ONE
         TEMP = TWO*TMP2
         DO 330 I = 1, M
            FVEC(I) = FVEC(I) + TMP2
            TI = TEMP*TMP2 - TMP1
            TMP1 = TMP2
            TMP2 = TI
  330       CONTINUE
  340    CONTINUE
      DX = ONE/FLOAT(N)
      IEV = -1
      DO 350 I = 1, M
         FVEC(I) = DX*FVEC(I)
         IF (IEV .GT. 0) FVEC(I) = FVEC(I) + ONE/(FLOAT(I)**2 - ONE)
         IEV = -IEV
  350    CONTINUE
      GO TO 430
C
C     BROWN ALMOST-LINEAR FUNCTION.
C
  360 CONTINUE
      SUM = -FLOAT(N+1)
      PROD = ONE
      DO 370 J = 1, N
         SUM = SUM + X(J)
         PROD = X(J)*PROD
  370    CONTINUE
      DO 380 I = 1, N
         FVEC(I) = X(I) + SUM
  380    CONTINUE
      FVEC(N) = PROD - ONE
      GO TO 430
C
C     OSBORNE 1 FUNCTION.
C
  390 CONTINUE
      DO 400 I = 1, 33
         TEMP = TEN*FLOAT(I-1)
         TMP1 = EXP(-X(4)*TEMP)
         TMP2 = EXP(-X(5)*TEMP)
         FVEC(I) = Y4(I) - (X(1) + X(2)*TMP1 + X(3)*TMP2)
  400    CONTINUE
      GO TO 430
C
C     OSBORNE 2 FUNCTION.
C
  410 CONTINUE
      DO 420 I = 1, 65
         TEMP = FLOAT(I-1)/TEN
         TMP1 = EXP(-X(5)*TEMP)
         TMP2 = EXP(-X(6)*(TEMP-X(9))**2)
         TMP3 = EXP(-X(7)*(TEMP-X(10))**2)
         TMP4 = EXP(-X(8)*(TEMP-X(11))**2)
         FVEC(I) = Y5(I)
     *             - (X(1)*TMP1 + X(2)*TMP2 + X(3)*TMP3 + X(4)*TMP4)
  420    CONTINUE
  430 CONTINUE
      RETURN
C
C     LAST CARD OF SUBROUTINE SSQFCN.
C
      END
      SUBROUTINE SSQJAC(M,N,X,FJAC,LDFJAC,NPROB)                        00000010
      INTEGER M,N,LDFJAC,NPROB
      REAL X(N),FJAC(LDFJAC,N)
C     **********
C
C     SUBROUTINE SSQJAC
C
C     THIS SUBROUTINE DEFINES THE JACOBIAN MATRICES OF EIGHTEEN
C     NONLINEAR LEAST SQUARES PROBLEMS. THE PROBLEM DIMENSIONS ARE
C     AS DESCRIBED IN THE PROLOGUE COMMENTS OF SSQFCN.
C
C     THE SUBROUTINE STATEMENT IS
C
C       SUBROUTINE SSQJAC(M,N,X,FJAC,LDFJAC,NPROB)
C
C     WHERE
C
C       M AND N ARE POSITIVE INTEGER VARIABLES. N MUST NOT EXCEED M.
C
C       X IS AN ARRAY OF LENGTH N.
C
C       FJAC IS AN M BY N OUTPUT ARRAY WHICH CONTAINS THE
C         JACOBIAN MATRIX OF THE NPROB FUNCTION EVALUATED AT X.
C
C       LDFJAC IS A POSITIVE INTEGER VARIABLE NOT LESS THAN M
C         WHICH SPECIFIES THE LEADING DIMENSION OF THE ARRAY FJAC.
C
C       NPROB IS A POSITIVE INTEGER VARIABLE WHICH DEFINES THE
C         NUMBER OF THE PROBLEM. NPROB MUST NOT EXCEED 18.
C
C     SUBPROGRAMS REQUIRED
C
C       FORTRAN-SUPPLIED ... ATAN,COS,EXP,SIN,SQRT
C
C     MINPACK. VERSION OF JULY 1978.
C     BURTON S. GARBOW, KENNETH E. HILLSTROM, JORGE J. MORE
C
C     **********
      INTEGER I,IVAR,J,K,MM1,NM1
      REAL C14,C20,C29,C45,C100,DIV,DX,EIGHT,FIVE,FOUR,ONE,PROD,S2,
     *     TEMP,TEN,THREE,TI,TMP1,TMP2,TMP3,TMP4,TPI,TWO,ZERO
      REAL V(11)
      REAL FLOAT
      DATA ZERO,ONE,TWO,THREE,FOUR,FIVE,EIGHT,TEN,C14,C20,C29,C45,C100
     *     /0.0E0,1.0E0,2.0E0,3.0E0,4.0E0,5.0E0,8.0E0,1.0E1,1.4E1,
     *      2.0E1,2.9E1,4.5E1,1.0E2/
      DATA V(1),V(2),V(3),V(4),V(5),V(6),V(7),V(8),V(9),V(10),V(11)
     *     /4.0E0,2.0E0,1.0E0,5.0E-1,2.5E-1,1.67E-1,1.25E-1,1.0E-1,
     *      8.33E-2,7.14E-2,6.25E-2/
      FLOAT(IVAR) = IVAR
C
C     JACOBIAN ROUTINE SELECTOR.
C
      GO TO (10,40,70,130,140,150,180,190,210,230,250,310,330,350,370,
     *       400,460,480), NPROB
C
C     LINEAR FUNCTION - FULL RANK.
C
   10 CONTINUE
      TEMP = TWO/FLOAT(M)
      DO 30 J = 1, N
         DO 20 I = 1, M
            FJAC(I,J) = -TEMP
   20       CONTINUE
         FJAC(J,J) = FJAC(J,J) + ONE
   30    CONTINUE
      GO TO 500
C
C     LINEAR FUNCTION - RANK 1.
C
   40 CONTINUE
      DO 60 J = 1, N
         DO 50 I = 1, M
            FJAC(I,J) = FLOAT(I)*FLOAT(J)
   50       CONTINUE
   60    CONTINUE
      GO TO 500
C
C     LINEAR FUNCTION - RANK 1 WITH ZERO COLUMNS AND ROWS.
C
   70 CONTINUE
      DO 90 J = 1, N
         DO 80 I = 1, M
            FJAC(I,J) = ZERO
   80       CONTINUE
   90    CONTINUE
      NM1 = N - 1
      MM1 = M - 1
      IF (NM1 .LT. 2) GO TO 120
      DO 110 J = 2, NM1
         DO 100 I = 2, MM1
            FJAC(I,J) = FLOAT(I-1)*FLOAT(J)
  100       CONTINUE
  110    CONTINUE
  120 CONTINUE
      GO TO 500
C
C     ROSENBROCK FUNCTION.
C
  130 CONTINUE
      FJAC(1,1) = -C20*X(1)
      FJAC(1,2) = TEN
      FJAC(2,1) = -ONE
      FJAC(2,2) = ZERO
      GO TO 500
C
C     HELICAL VALLEY FUNCTION.
C
  140 CONTINUE
      TPI = EIGHT*ATAN(ONE)
      TEMP = X(1)**2 + X(2)**2
      TMP1 = TPI*TEMP
      TMP2 = SQRT(TEMP)
      FJAC(1,1) = C100*X(2)/TMP1
      FJAC(1,2) = -C100*X(1)/TMP1
      FJAC(1,3) = TEN
      FJAC(2,1) = TEN*X(1)/TMP2
      FJAC(2,2) = TEN*X(2)/TMP2
      FJAC(2,3) = ZERO
      FJAC(3,1) = ZERO
      FJAC(3,2) = ZERO
      FJAC(3,3) = ONE
      GO TO 500
C
C     POWELL SINGULAR FUNCTION.
C
  150 CONTINUE
      DO 170 J = 1, 4
         DO 160 I = 1, 4
            FJAC(I,J) = ZERO
  160       CONTINUE
  170    CONTINUE
      FJAC(1,1) = ONE
      FJAC(1,2) = TEN
      FJAC(2,3) = SQRT(FIVE)
      FJAC(2,4) = -FJAC(2,3)
      FJAC(3,2) = TWO*(X(2) - TWO*X(3))
      FJAC(3,3) = -TWO*FJAC(3,2)
      FJAC(4,1) = TWO*SQRT(TEN)*(X(1) - X(4))
      FJAC(4,4) = -FJAC(4,1)
      GO TO 500
C
C     FREUDENSTEIN AND ROTH FUNCTION.
C
  180 CONTINUE
      FJAC(1,1) = ONE
      FJAC(1,2) = X(2)*(TEN - THREE*X(2)) - TWO
      FJAC(2,1) = ONE
      FJAC(2,2) = X(2)*(TWO + THREE*X(2)) - C14
      GO TO 500
C
C     BARD FUNCTION.
C
  190 CONTINUE
      DO 200 I = 1, 15
         TMP1 = FLOAT(I)
         TMP2 = FLOAT(16-I)
         TMP3 = TMP1
         IF (I .GT. 8) TMP3 = TMP2
         TMP4 = (X(2)*TMP2 + X(3)*TMP3)**2
         FJAC(I,1) = -ONE
         FJAC(I,2) = TMP1*TMP2/TMP4
         FJAC(I,3) = TMP1*TMP3/TMP4
  200    CONTINUE
      GO TO 500
C
C     KOWALIK AND OSBORNE FUNCTION.
C
  210 CONTINUE
      DO 220 I = 1, 11
         TMP1 = V(I)*(V(I) + X(2))
         TMP2 = V(I)*(V(I) + X(3)) + X(4)
         FJAC(I,1) = -TMP1/TMP2
         FJAC(I,2) = -V(I)*X(1)/TMP2
         FJAC(I,3) = FJAC(I,1)*FJAC(I,2)
         FJAC(I,4) = FJAC(I,3)/V(I)
  220    CONTINUE
      GO TO 500
C
C     MEYER FUNCTION.
C
  230 CONTINUE
      DO 240 I = 1, 16
         TEMP = FIVE*FLOAT(I) + C45 + X(3)
         TMP1 = X(2)/TEMP
         TMP2 = EXP(TMP1)
         FJAC(I,1) = TMP2
         FJAC(I,2) = X(1)*TMP2/TEMP
         FJAC(I,3) = -TMP1*FJAC(I,2)
  240    CONTINUE
      GO TO 500
C
C     WATSON FUNCTION.
C
  250 CONTINUE
      DO 280 I = 1, 29
         DIV = FLOAT(I)/C29
         S2 = ZERO
         DX = ONE
         DO 260 J = 1, N
            S2 = S2 + DX*X(J)
            DX = DIV*DX
  260       CONTINUE
         TEMP = TWO*DIV*S2
         DX = ONE/DIV
         DO 270 J = 1, N
            FJAC(I,J) = DX*(FLOAT(J-1) - TEMP)
            DX = DIV*DX
  270       CONTINUE
  280    CONTINUE
      DO 300 J = 1, N
         DO 290 I = 30, 31
            FJAC(I,J) = ZERO
  290       CONTINUE
  300    CONTINUE
      FJAC(30,1) = ONE
      FJAC(31,1) = -TWO*X(1)
      FJAC(31,2) = ONE
      GO TO 500
C
C     BOX 3-DIMENSIONAL FUNCTION.
C
  310 CONTINUE
      DO 320 I = 1, M
         TEMP = FLOAT(I)
         TMP1 = TEMP/TEN
         FJAC(I,1) = -TMP1*EXP(-TMP1*X(1))
         FJAC(I,2) = TMP1*EXP(-TMP1*X(2))
         FJAC(I,3) = EXP(-TEMP) - EXP(-TMP1)
  320    CONTINUE
      GO TO 500
C
C     JENNRICH AND SAMPSON FUNCTION.
C
  330 CONTINUE
      DO 340 I = 1, M
         TEMP = FLOAT(I)
         FJAC(I,1) = -TEMP*EXP(TEMP*X(1))
         FJAC(I,2) = -TEMP*EXP(TEMP*X(2))
  340    CONTINUE
      GO TO 500
C
C     BROWN AND DENNIS FUNCTION.
C
  350 CONTINUE
      DO 360 I = 1, M
         TEMP = FLOAT(I)/FIVE
         TI = SIN(TEMP)
         TMP1 = X(1) + TEMP*X(2) - EXP(TEMP)
         TMP2 = X(3) + TI*X(4) - COS(TEMP)
         FJAC(I,1) = TWO*TMP1
         FJAC(I,2) = TEMP*FJAC(I,1)
         FJAC(I,3) = TWO*TMP2
         FJAC(I,4) = TI*FJAC(I,3)
  360    CONTINUE
      GO TO 500
C
C     CHEBYQUAD FUNCTION.
C
  370 CONTINUE
      DX = ONE/FLOAT(N)
      DO 390 J = 1, N
         TMP1 = ONE
         TMP2 = TWO*X(J) - ONE
         TEMP = TWO*TMP2
         TMP3 = ZERO
         TMP4 = TWO
         DO 380 I = 1, M
            FJAC(I,J) = DX*TMP4
            TI = FOUR*TMP2 + TEMP*TMP4 - TMP3
            TMP3 = TMP4
            TMP4 = TI
            TI = TEMP*TMP2 - TMP1
            TMP1 = TMP2
            TMP2 = TI
  380       CONTINUE
  390    CONTINUE
      GO TO 500
C
C     BROWN ALMOST-LINEAR FUNCTION.
C
  400 CONTINUE
      PROD = ONE
      DO 420 J = 1, N
         PROD = X(J)*PROD
         DO 410 I = 1, N
            FJAC(I,J) = ONE
  410       CONTINUE
         FJAC(J,J) = TWO
  420    CONTINUE
      DO 450 J = 1, N
         TEMP = X(J)
         IF (TEMP .NE. ZERO) GO TO 440
         TEMP = ONE
         PROD = ONE
         DO 430 K = 1, N
            IF (K .NE. J) PROD = X(K)*PROD
  430       CONTINUE
  440    CONTINUE
         FJAC(N,J) = PROD/TEMP
  450    CONTINUE
      GO TO 500
C
C     OSBORNE 1 FUNCTION.
C
  460 CONTINUE
      DO 470 I = 1, 33
         TEMP = TEN*FLOAT(I-1)
         TMP1 = EXP(-X(4)*TEMP)
         TMP2 = EXP(-X(5)*TEMP)
         FJAC(I,1) = -ONE
         FJAC(I,2) = -TMP1
         FJAC(I,3) = -TMP2
         FJAC(I,4) = TEMP*X(2)*TMP1
         FJAC(I,5) = TEMP*X(3)*TMP2
  470    CONTINUE
      GO TO 500
C
C     OSBORNE 2 FUNCTION.
C
  480 CONTINUE
      DO 490 I = 1, 65
         TEMP = FLOAT(I-1)/TEN
         TMP1 = EXP(-X(5)*TEMP)
         TMP2 = EXP(-X(6)*(TEMP-X(9))**2)
         TMP3 = EXP(-X(7)*(TEMP-X(10))**2)
         TMP4 = EXP(-X(8)*(TEMP-X(11))**2)
         FJAC(I,1) = -TMP1
         FJAC(I,2) = -TMP2
         FJAC(I,3) = -TMP3
         FJAC(I,4) = -TMP4
         FJAC(I,5) = TEMP*X(1)*TMP1
         FJAC(I,6) = X(2)*(TEMP - X(9))**2*TMP2
         FJAC(I,7) = X(3)*(TEMP - X(10))**2*TMP3
         FJAC(I,8) = X(4)*(TEMP - X(11))**2*TMP4
         FJAC(I,9) = -TWO*X(2)*X(6)*(TEMP - X(9))*TMP2
         FJAC(I,10) = -TWO*X(3)*X(7)*(TEMP - X(10))*TMP3
         FJAC(I,11) = -TWO*X(4)*X(8)*(TEMP - X(11))*TMP4
  490    CONTINUE
  500 CONTINUE
      RETURN
C
C     LAST CARD OF SUBROUTINE SSQJAC.
C
      END
C ===== 7. SINGLE PRECISION TESTING AIDS FOR UNCONSTRAINED NONLINEAR
C =====     OPTIMIZATION.
      SUBROUTINE INITPT(N,X,NPROB,FACTOR)                               00000010
      INTEGER N,NPROB
      REAL FACTOR
      REAL X(N)
C     **********
C
C     SUBROUTINE INITPT
C
C     THIS SUBROUTINE SPECIFIES THE STANDARD STARTING POINTS FOR THE
C     FUNCTIONS DEFINED BY SUBROUTINE OBJFCN. THE SUBROUTINE RETURNS
C     IN X A MULTIPLE (FACTOR) OF THE STANDARD STARTING POINT. FOR
C     THE SEVENTH FUNCTION THE STANDARD STARTING POINT IS ZERO, SO IN
C     THIS CASE, IF FACTOR IS NOT UNITY, THEN THE SUBROUTINE RETURNS
C     THE VECTOR  X(J) = FACTOR, J=1,...,N.
C
C     THE SUBROUTINE STATEMENT IS
C
C       SUBROUTINE INITPT(N,X,NPROB,FACTOR)
C
C     WHERE
C
C       N IS A POSITIVE INTEGER VARIABLE.
C
C       X IS AN OUTPUT ARRAY OF LENGTH N WHICH CONTAINS THE
C         STANDARD STARTING POINT FOR PROBLEM NPROB MULTIPLIED BY
C         FACTOR.
C
C       NPROB IS A POSITIVE INTEGER VARIABLE WHICH DEFINES THE
C         NUMBER OF THE PROBLEM. NPROB MUST NOT EXCEED 18.
C
C       FACTOR SPECIFIES THE MULTIPLE OF THE STANDARD STARTING
C         POINT. IF FACTOR IS UNITY, NO MULTIPLICATION IS PERFORMED.
C
C     MINPACK. VERSION OF JULY 1978.
C     BURTON S. GARBOW, KENNETH E. HILLSTROM, JORGE J. MORE
C
C     **********
      INTEGER IVAR,J
      REAL C1,C2,C3,C4,FIVE,H,HALF,ONE,TEN,THREE,TWENTY,TWNTF,TWO,ZERO
      REAL FLOAT
      DATA ZERO,HALF,ONE,TWO,THREE,FIVE,TEN,TWENTY,TWNTF
     *     /0.0E0,0.5E0,1.0E0,2.0E0,3.0E0,5.0E0,1.0E1,2.0E1,2.5E1/
      DATA C1,C2,C3,C4 /4.0E-1,2.5E0,1.5E-1,1.2E0/
      FLOAT(IVAR) = IVAR
C
C     SELECTION OF INITIAL POINT.
C
      GO TO (10,20,30,40,50,60,80,100,120,140,150,160,170,190,210,230,
     *       240,250), NPROB
C
C     HELICAL VALLEY FUNCTION.
C
   10 CONTINUE
      X(1) = -ONE
      X(2) = ZERO
      X(3) = ZERO
      GO TO 270
C
C     BIGGS EXP6 FUNCTION.
C
   20 CONTINUE
      X(1) = ONE
      X(2) = TWO
      X(3) = ONE
      X(4) = ONE
      X(5) = ONE
      X(6) = ONE
      GO TO 270
C
C     GAUSSIAN FUNCTION.
C
   30 CONTINUE
      X(1) = C1
      X(2) = ONE
      X(3) = ZERO
      GO TO 270
C
C     POWELL BADLY SCALED FUNCTION.
C
   40 CONTINUE
      X(1) = ZERO
      X(2) = ONE
      GO TO 270
C
C     BOX 3-DIMENSIONAL FUNCTION.
C
   50 CONTINUE
      X(1) = ZERO
      X(2) = TEN
      X(3) = TWENTY
      GO TO 270
C
C     VARIABLY DIMENSIONED FUNCTION.
C
   60 CONTINUE
      H = ONE/FLOAT(N)
      DO 70 J = 1, N
         X(J) = ONE - FLOAT(J)*H
   70    CONTINUE
      GO TO 270
C
C     WATSON FUNCTION.
C
   80 CONTINUE
      DO 90 J = 1, N
         X(J) = ZERO
   90    CONTINUE
      GO TO 270
C
C     PENALTY FUNCTION I.
C
  100 CONTINUE
      DO 110 J = 1, N
         X(J) = FLOAT(J)
  110    CONTINUE
      GO TO 270
C
C     PENALTY FUNCTION II.
C
  120 CONTINUE
      DO 130 J = 1, N
         X(J) = HALF
  130    CONTINUE
      GO TO 270
C
C     BROWN BADLY SCALED FUNCTION.
C
  140 CONTINUE
      X(1) = ONE
      X(2) = ONE
      GO TO 270
C
C     BROWN AND DENNIS FUNCTION.
C
  150 CONTINUE
      X(1) = TWNTF
      X(2) = FIVE
      X(3) = -FIVE
      X(4) = -ONE
      GO TO 270
C
C     GULF RESEARCH AND DEVELOPMENT FUNCTION.
C
  160 CONTINUE
      X(1) = FIVE
      X(2) = C2
      X(3) = C3
      GO TO 270
C
C     TRIGONOMETRIC FUNCTION.
C
  170 CONTINUE
      H = ONE/FLOAT(N)
      DO 180 J = 1, N
         X(J) = H
  180    CONTINUE
      GO TO 270
C
C     EXTENDED ROSENBROCK FUNCTION.
C
  190 CONTINUE
      DO 200 J = 1, N, 2
         X(J) = -C4
         X(J+1) = ONE
  200    CONTINUE
      GO TO 270
C
C     EXTENDED POWELL SINGULAR FUNCTION.
C
  210 CONTINUE
      DO 220 J = 1, N, 4
         X(J) = THREE
         X(J+1) = -ONE
         X(J+2) = ZERO
         X(J+3) = ONE
  220    CONTINUE
      GO TO 270
C
C     BEALE FUNCTION.
C
  230 CONTINUE
      X(1) = ONE
      X(2) = ONE
      GO TO 270
C
C     WOOD FUNCTION.
C
  240 CONTINUE
      X(1) = -THREE
      X(2) = -ONE
      X(3) = -THREE
      X(4) = -ONE
      GO TO 270
C
C     CHEBYQUAD FUNCTION.
C
  250 CONTINUE
      H = ONE/FLOAT(N+1)
      DO 260 J = 1, N
         X(J) = FLOAT(J)*H
  260    CONTINUE
  270 CONTINUE
C
C     COMPUTE MULTIPLE OF INITIAL POINT.
C
      IF (FACTOR .EQ. ONE) GO TO 320
      IF (NPROB .EQ. 7) GO TO 290
         DO 280 J = 1, N
            X(J) = FACTOR*X(J)
  280       CONTINUE
         GO TO 310
  290 CONTINUE
         DO 300 J = 1, N
            X(J) = FACTOR
  300       CONTINUE
  310 CONTINUE
  320 CONTINUE
      RETURN
C
C     LAST CARD OF SUBROUTINE INITPT.
C
      END
      SUBROUTINE OBJFCN(N,X,F,NPROB)                                    00000010
      INTEGER N,NPROB
      REAL F
      REAL X(N)
C     **********
C
C     SUBROUTINE OBJFCN
C
C     THIS SUBROUTINE DEFINES THE OBJECTIVE FUNCTIONS OF EIGHTEEN
C     NONLINEAR UNCONSTRAINED MINIMIZATION PROBLEMS. THE VALUES
C     OF N FOR FUNCTIONS 1,2,3,4,5,10,11,12,16 AND 17 ARE
C     3,6,3,2,3,2,4,3,2 AND 4, RESPECTIVELY.
C     FOR FUNCTION 7, N MAY BE 2 OR GREATER BUT IS USUALLY 6 OR 9.
C     FOR FUNCTIONS 6,8,9,13,14,15 AND 18 N MAY BE VARIABLE,
C     HOWEVER IT MUST BE EVEN FOR FUNCTION 14, A MULTIPLE OF 4 FOR
C     FUNCTION 15, AND NOT GREATER THAN 50 FOR FUNCTION 18.
C
C     THE SUBROUTINE STATEMENT IS
C
C       SUBROUTINE OBJFCN(N,X,F,NPROB)
C
C     WHERE
C
C       N IS A POSITIVE INTEGER VARIABLE.
C
C       X IS AN ARRAY OF LENGTH N.
C
C       F IS AN OUTPUT VARIABLE WHICH CONTAINS THE VALUE OF
C         THE NPROB OBJECTIVE FUNCTION EVALUATED AT X.
C
C       NPROB IS A POSITIVE INTEGER INPUT VARIABLE WHICH DEFINES THE
C         NUMBER OF THE PROBLEM. NPROB MUST NOT EXCEED 18.
C
C     SUBPROGRAMS REQUIRED
C
C       FORTRAN-SUPPLIED ... ABS,ATAN,COS,EXP,ALOG,SIGN,SIN,
C                            SQRT
C
C     MINPACK. VERSION OF JULY 1978.
C     BURTON S. GARBOW, KENNETH E. HILLSTROM, JORGE J. MORE
C
C     **********
      INTEGER I,IEV,IVAR,J
      REAL AP,ARG,BP,C2PDM6,CP0001,CP1,CP2,CP25,CP5,C1P5,C2P25,C2P625,
     *     C3P5,C25,C29,C90,C100,C10000,C1PD6,D1,D2,EIGHT,FIFTY,FIVE,
     *     FOUR,ONE,R,S1,S2,S3,T,T1,T2,T3,TEN,TH,THREE,TPI,TWO,ZERO
      REAL FVEC(50),Y(15)
      REAL FLOAT
      DATA ZERO,ONE,TWO,THREE,FOUR,FIVE,EIGHT,TEN,FIFTY
     *     /0.0E0,1.0E0,2.0E0,3.0E0,4.0E0,5.0E0,8.0E0,1.0E1,5.0E1/
      DATA C2PDM6,CP0001,CP1,CP2,CP25,CP5,C1P5,C2P25,C2P625,C3P5,C25,
     *     C29,C90,C100,C10000,C1PD6
     *     /2.0E-6,1.0E-4,1.0E-1,2.0E-1,2.5E-1,5.0E-1,1.5E0,2.25E0,
     *      2.625E0,3.5E0,2.5E1,2.9E1,9.0E1,1.0E2,1.0E4,1.0E6/
      DATA AP,BP /1.0E-5,1.0E0/
      DATA Y(1),Y(2),Y(3),Y(4),Y(5),Y(6),Y(7),Y(8),Y(9),Y(10),Y(11),
     *     Y(12),Y(13),Y(14),Y(15)
     *     /9.0E-4,4.4E-3,1.75E-2,5.4E-2,1.295E-1,2.42E-1,3.521E-1,
     *      3.989E-1,3.521E-1,2.42E-1,1.295E-1,5.4E-2,1.75E-2,4.4E-3,
     *      9.0E-4/
      FLOAT(IVAR) = IVAR
C
C     FUNCTION ROUTINE SELECTOR.
C
      GO TO (10,20,40,60,70,90,110,150,170,200,210,230,250,280,300,
     *       320,330,340), NPROB
C
C     HELICAL VALLEY FUNCTION.
C
   10 CONTINUE
      TPI = EIGHT*ATAN(ONE)
      TH = SIGN(CP25,X(2))
      IF (X(1) .GT. ZERO) TH = ATAN(X(2)/X(1))/TPI
      IF (X(1) .LT. ZERO) TH = ATAN(X(2)/X(1))/TPI + CP5
      ARG = X(1)**2 + X(2)**2
      R = SQRT(ARG)
      T = X(3) - TEN*TH
      F = C100*(T**2 + (R - ONE)**2) + X(3)**2
      GO TO 390
C
C     BIGGS EXP6 FUNCTION.
C
   20 CONTINUE
      F = ZERO
      DO 30 I = 1, 13
         D1 = FLOAT(I)/TEN
         D2 = EXP(-D1) - FIVE*EXP(-TEN*D1) + THREE*EXP(-FOUR*D1)
         S1 = EXP(-D1*X(1))
         S2 = EXP(-D1*X(2))
         S3 = EXP(-D1*X(5))
         T = X(3)*S1 - X(4)*S2 + X(6)*S3 - D2
         F = F + T**2
   30    CONTINUE
      GO TO 390
C
C     GAUSSIAN FUNCTION.
C
   40 CONTINUE
      F = ZERO
      DO 50 I = 1, 15
         D1 = CP5*FLOAT(I-1)
         D2 = C3P5 - D1 - X(3)
         ARG = -CP5*X(2)*D2**2
         R = EXP(ARG)
         T = X(1)*R - Y(I)
         F = F + T**2
   50    CONTINUE
      GO TO 390
C
C     POWELL BADLY SCALED FUNCTION.
C
   60 CONTINUE
      T1 = C10000*X(1)*X(2) - ONE
      S1 = EXP(-X(1))
      S2 = EXP(-X(2))
      T2 = S1 + S2 - ONE - CP0001
      F = T1**2 + T2**2
      GO TO 390
C
C     BOX 3-DIMENSIONAL FUNCTION.
C
   70 CONTINUE
      F = ZERO
      DO 80 I = 1, 10
         D1 = FLOAT(I)
         D2 = D1/TEN
         S1 = EXP(-D2*X(1))
         S2 = EXP(-D2*X(2))
         S3 = EXP(-D2) - EXP(-D1)
         T = S1 - S2 - S3*X(3)
         F = F + T**2
   80    CONTINUE
      GO TO 390
C
C     VARIABLY DIMENSIONED FUNCTION.
C
   90 CONTINUE
      T1 = ZERO
      T2 = ZERO
      DO 100 J = 1, N
         T1 = T1 + FLOAT(J)*(X(J) - ONE)
         T2 = T2 + (X(J) - ONE)**2
  100    CONTINUE
      F = T2 + T1**2*(ONE + T1**2)
      GO TO 390
C
C     WATSON FUNCTION.
C
  110 CONTINUE
      F = ZERO
      DO 140 I = 1, 29
         D1 = FLOAT(I)/C29
         S1 = ZERO
         D2 = ONE
         DO 120 J = 2, N
            S1 = S1 + FLOAT(J-1)*D2*X(J)
            D2 = D1*D2
  120       CONTINUE
         S2 = ZERO
         D2 = ONE
         DO 130 J = 1, N
            S2 = S2 + D2*X(J)
            D2 = D1*D2
  130       CONTINUE
         T = S1 - S2**2 - ONE
         F = F + T**2
  140    CONTINUE
      T1 = X(2) - X(1)**2 - ONE
      F = F + X(1)**2 + T1**2
      GO TO 390
C
C     PENALTY FUNCTION I.
C
  150 CONTINUE
      T1 = -CP25
      T2 = ZERO
      DO 160 J = 1, N
         T1 = T1 + X(J)**2
         T2 = T2 + (X(J) - ONE)**2
  160    CONTINUE
      F = AP*T2 + BP*T1**2
      GO TO 390
C
C     PENALTY FUNCTION II.
C
  170 CONTINUE
      T1 = -ONE
      T2 = ZERO
      T3 = ZERO
      D1 = EXP(CP1)
      D2 = ONE
      DO 190 J = 1, N
         T1 = T1 + FLOAT(N-J+1)*X(J)**2
         S1 = EXP(X(J)/TEN)
         IF (J .EQ. 1) GO TO 180
         S3 = S1 + S2 - D2*(D1 + ONE)
         T2 = T2 + S3**2
         T3 = T3 + (S1 - ONE/D1)**2
  180    CONTINUE
         S2 = S1
         D2 = D1*D2
  190    CONTINUE
      F = AP*(T2 + T3) + BP*(T1**2 + (X(1) - CP2)**2)
      GO TO 390
C
C     BROWN BADLY SCALED FUNCTION.
C
  200 CONTINUE
      T1 = X(1) - C1PD6
      T2 = X(2) - C2PDM6
      T3 = X(1)*X(2) - TWO
      F = T1**2 + T2**2 + T3**2
      GO TO 390
C
C     BROWN AND DENNIS FUNCTION.
C
  210 CONTINUE
      F = ZERO
      DO 220 I = 1, 20
         D1 = FLOAT(I)/FIVE
         D2 = SIN(D1)
         T1 = X(1) + D1*X(2) - EXP(D1)
         T2 = X(3) + D2*X(4) - COS(D1)
         T = T1**2 + T2**2
         F = F + T**2
  220    CONTINUE
      GO TO 390
C
C     GULF RESEARCH AND DEVELOPMENT FUNCTION.
C
  230 CONTINUE
      F = ZERO
      D1 = TWO/THREE
      DO 240 I = 1, 99
         ARG = FLOAT(I)/C100
         R = ABS((-FIFTY*ALOG(ARG))**D1+C25-X(2))
         T1 = R**X(3)/X(1)
         T2 = EXP(-T1)
         T = T2 - ARG
         F = F + T**2
  240    CONTINUE
      GO TO 390
C
C     TRIGONOMETRIC FUNCTION.
C
  250 CONTINUE
      S1 = ZERO
      DO 260 J = 1, N
         S1 = S1 + COS(X(J))
  260    CONTINUE
      F = ZERO
      DO 270 J = 1, N
         T = FLOAT(N+J) - SIN(X(J)) - S1 - FLOAT(J)*COS(X(J))
         F = F + T**2
  270    CONTINUE
      GO TO 390
C
C     EXTENDED ROSENBROCK FUNCTION.
C
  280 CONTINUE
      F = ZERO
      DO 290 J = 1, N, 2
         T1 = ONE - X(J)
         T2 = TEN*(X(J+1) - X(J)**2)
         F = F + T1**2 + T2**2
  290    CONTINUE
      GO TO 390
C
C     EXTENDED POWELL FUNCTION.
C
  300 CONTINUE
      F = ZERO
      DO 310 J = 1, N, 4
         T = X(J) + TEN*X(J+1)
         T1 = X(J+2) - X(J+3)
         S1 = FIVE*T1
         T2 = X(J+1) - TWO*X(J+2)
         S2 = T2**3
         T3 = X(J) - X(J+3)
         S3 = TEN*T3**3
         F = F + T**2 + S1*T1 + S2*T2 + S3*T3
  310    CONTINUE
      GO TO 390
C
C     BEALE FUNCTION.
C
  320 CONTINUE
      S1 = ONE - X(2)
      T1 = C1P5 - X(1)*S1
      S2 = ONE - X(2)**2
      T2 = C2P25 - X(1)*S2
      S3 = ONE - X(2)**3
      T3 = C2P625 - X(1)*S3
      F = T1**2 + T2**2 + T3**2
      GO TO 390
C
C     WOOD FUNCTION.
C
  330 CONTINUE
      S1 = X(2) - X(1)**2
      S2 = ONE - X(1)
      S3 = X(2) - ONE
      T1 = X(4) - X(3)**2
      T2 = ONE - X(3)
      T3 = X(4) - ONE
      F = C100*S1**2 + S2**2 + C90*T1**2 + T2**2 + TEN*(S3 + T3)**2
     *    + (S3 - T3)**2/TEN
      GO TO 390
C
C     CHEBYQUAD FUNCTION.
C
  340 CONTINUE
      DO 350 I = 1, N
         FVEC(I) = ZERO
  350    CONTINUE
      DO 370 J = 1, N
         T1 = ONE
         T2 = TWO*X(J) - ONE
         T = TWO*T2
         DO 360 I = 1, N
            FVEC(I) = FVEC(I) + T2
            TH = T*T2 - T1
            T1 = T2
            T2 = TH
  360       CONTINUE
  370    CONTINUE
      F = ZERO
      D1 = ONE/FLOAT(N)
      IEV = -1
      DO 380 I = 1, N
         T = D1*FVEC(I)
         IF (IEV .GT. 0) T = T + ONE/(FLOAT(I)**2 - ONE)
         F = F + T**2
         IEV = -IEV
  380    CONTINUE
  390 CONTINUE
      RETURN
C
C     LAST CARD OF SUBROUTINE OBJFCN.
C
      END
      SUBROUTINE GRDFCN(N,X,G,NPROB)                                    00000010
      INTEGER N,NPROB
      REAL X(N),G(N)
C     **********
C
C     SUBROUTINE GRDFCN
C
C     THIS SUBROUTINE DEFINES THE GRADIENT VECTORS OF EIGHTEEN
C     NONLINEAR UNCONSTRAINED MINIMIZATION PROBLEMS. THE PROBLEM
C     DIMENSIONS ARE AS DESCRIBED IN THE PROLOGUE COMMENTS OF OBJFCN.
C
C     THE SUBROUTINE STATEMENT IS
C
C       SUBROUTINE GRDFCN(N,X,G,NPROB)
C
C     WHERE
C
C       N IS A POSITIVE INTEGER VARIABLE.
C
C       X IS AN ARRAY OF LENGTH N.
C
C       G IS AN OUTPUT ARRAY OF LENGTH N WHICH CONTAINS
C         THE COMPONENTS OF THE GRADIENT VECTOR OF THE NPROB
C         OBJECTIVE FUNCTION EVALUATED AT X.
C
C       NPROB IS A POSITIVE INTEGER INPUT VARIABLE WHICH DEFINES THE
C         NUMBER OF THE PROBLEM. NPROB MUST NOT EXCEED 18.
C
C     SUBPROGRAMS REQUIRED
C
C       FORTRAN-SUPPLIED ... ABS,ATAN,COS,EXP,ALOG,SIGN,SIN,
C                            SQRT
C
C     MINPACK. VERSION OF JULY 1978.
C     BURTON S. GARBOW, KENNETH E. HILLSTROM, JORGE J. MORE
C
C     **********
      INTEGER I,IEV,IVAR,J
      REAL AP,ARG,BP,C2PDM6,CP0001,CP1,CP2,CP25,CP5,C1P5,C2P25,C2P625,
     *     C3P5,C19P8,C20P2,C25,C29,C100,C180,C200,C10000,C1PD6,D1,D2,
     *     EIGHT,FIFTY,FIVE,FOUR,ONE,R,S1,S2,S3,T,T1,T2,T3,TEN,TH,
     *     THREE,TPI,TWENTY,TWO,ZERO
      REAL FVEC(50),Y(15)
      REAL FLOAT
      DATA ZERO,ONE,TWO,THREE,FOUR,FIVE,EIGHT,TEN,TWENTY,FIFTY
     *     /0.0E0,1.0E0,2.0E0,3.0E0,4.0E0,5.0E0,8.0E0,1.0E1,2.0E1,
     *      5.0E1/
      DATA C2PDM6,CP0001,CP1,CP2,CP25,CP5,C1P5,C2P25,C2P625,C3P5,
     *     C19P8,C20P2,C25,C29,C100,C180,C200,C10000,C1PD6
     *     /2.0E-6,1.0E-4,1.0E-1,2.0E-1,2.5E-1,5.0E-1,1.5E0,2.25E0,
     *      2.625E0,3.5E0,1.98E1,2.02E1,2.5E1,2.9E1,1.0E2,1.8E2,2.0E2,
     *      1.0E4,1.0E6/
      DATA AP,BP /1.0E-5,1.0E0/
      DATA Y(1),Y(2),Y(3),Y(4),Y(5),Y(6),Y(7),Y(8),Y(9),Y(10),Y(11),
     *     Y(12),Y(13),Y(14),Y(15)
     *     /9.0E-4,4.4E-3,1.75E-2,5.4E-2,1.295E-1,2.42E-1,3.521E-1,
     *      3.989E-1,3.521E-1,2.42E-1,1.295E-1,5.4E-2,1.75E-2,4.4E-3,
     *      9.0E-4/
      FLOAT(IVAR) = IVAR
C
C     GRADIENT ROUTINE SELECTOR.
C
      GO TO (10,20,50,70,80,100,130,190,220,260,270,290,310,350,370,
     *       390,400,410), NPROB
C
C     HELICAL VALLEY FUNCTION.
C
   10 CONTINUE
      TPI = EIGHT*ATAN(ONE)
      TH = SIGN(CP25,X(2))
      IF (X(1) .GT. ZERO) TH = ATAN(X(2)/X(1))/TPI
      IF (X(1) .LT. ZERO) TH = ATAN(X(2)/X(1))/TPI + CP5
      ARG = X(1)**2 + X(2)**2
      R = SQRT(ARG)
      T = X(3) - TEN*TH
      S1 = TEN*T/(TPI*ARG)
      G(1) = C200*(X(1) - X(1)/R + X(2)*S1)
      G(2) = C200*(X(2) - X(2)/R - X(1)*S1)
      G(3) = TWO*(C100*T + X(3))
      GO TO 490
C
C     BIGGS EXP6 FUNCTION.
C
   20 CONTINUE
      DO 30 J = 1, N
         G(J) = ZERO
   30    CONTINUE
      DO 40 I = 1, 13
         D1 = FLOAT(I)/TEN
         D2 = EXP(-D1) - FIVE*EXP(-TEN*D1) + THREE*EXP(-FOUR*D1)
         S1 = EXP(-D1*X(1))
         S2 = EXP(-D1*X(2))
         S3 = EXP(-D1*X(5))
         T = X(3)*S1 - X(4)*S2 + X(6)*S3 - D2
         TH = D1*T
         G(1) = G(1) - S1*TH
         G(2) = G(2) + S2*TH
         G(3) = G(3) + S1*T
         G(4) = G(4) - S2*T
         G(5) = G(5) - S3*TH
         G(6) = G(6) + S3*T
   40    CONTINUE
      G(1) = TWO*X(3)*G(1)
      G(2) = TWO*X(4)*G(2)
      G(3) = TWO*G(3)
      G(4) = TWO*G(4)
      G(5) = TWO*X(6)*G(5)
      G(6) = TWO*G(6)
      GO TO 490
C
C     GAUSSIAN FUNCTION.
C
   50 CONTINUE
      G(1) = ZERO
      G(2) = ZERO
      G(3) = ZERO
      DO 60 I = 1, 15
         D1 = CP5*FLOAT(I-1)
         D2 = C3P5 - D1 - X(3)
         ARG = -CP5*X(2)*D2**2
         R = EXP(ARG)
         T = X(1)*R - Y(I)
         S1 = R*T
         S2 = D2*S1
         G(1) = G(1) + S1
         G(2) = G(2) - D2*S2
         G(3) = G(3) + S2
   60    CONTINUE
      G(1) = TWO*G(1)
      G(2) = X(1)*G(2)
      G(3) = TWO*X(1)*X(2)*G(3)
      GO TO 490
C
C     POWELL BADLY SCALED FUNCTION.
C
   70 CONTINUE
      T1 = C10000*X(1)*X(2) - ONE
      S1 = EXP(-X(1))
      S2 = EXP(-X(2))
      T2 = S1 + S2 - ONE - CP0001
      G(1) = TWO*(C10000*X(2)*T1 - S1*T2)
      G(2) = TWO*(C10000*X(1)*T1 - S2*T2)
      GO TO 490
C
C     BOX 3-DIMENSIONAL FUNCTION.
C
   80 CONTINUE
      G(1) = ZERO
      G(2) = ZERO
      G(3) = ZERO
      DO 90 I = 1, 10
         D1 = FLOAT(I)
         D2 = D1/TEN
         S1 = EXP(-D2*X(1))
         S2 = EXP(-D2*X(2))
         S3 = EXP(-D2) - EXP(-D1)
         T = S1 - S2 - S3*X(3)
         TH = D2*T
         G(1) = G(1) - S1*TH
         G(2) = G(2) + S2*TH
         G(3) = G(3) - S3*T
   90    CONTINUE
      G(1) = TWO*G(1)
      G(2) = TWO*G(2)
      G(3) = TWO*G(3)
      GO TO 490
C
C     VARIABLY DIMENSIONED FUNCTION.
C
  100 CONTINUE
      T1 = ZERO
      DO 110 J = 1, N
         T1 = T1 + FLOAT(J)*(X(J) - ONE)
  110    CONTINUE
      T = T1*(ONE + TWO*T1**2)
      DO 120 J = 1, N
         G(J) = TWO*(X(J) - ONE + FLOAT(J)*T)
  120    CONTINUE
      GO TO 490
C
C     WATSON FUNCTION.
C
  130 CONTINUE
      DO 140 J = 1, N
         G(J) = ZERO
  140    CONTINUE
      DO 180 I = 1, 29
         D1 = FLOAT(I)/C29
         S1 = ZERO
         D2 = ONE
         DO 150 J = 2, N
            S1 = S1 + FLOAT(J-1)*D2*X(J)
            D2 = D1*D2
  150       CONTINUE
         S2 = ZERO
         D2 = ONE
         DO 160 J = 1, N
            S2 = S2 + D2*X(J)
            D2 = D1*D2
  160       CONTINUE
         T = S1 - S2**2 - ONE
         S3 = TWO*D1*S2
         D2 = TWO/D1
         DO 170 J = 1, N
            G(J) = G(J) + D2*(FLOAT(J-1) - S3)*T
            D2 = D1*D2
  170       CONTINUE
  180    CONTINUE
      T1 = X(2) - X(1)**2 - ONE
      G(1) = G(1) + X(1)*(TWO - FOUR*T1)
      G(2) = G(2) + TWO*T1
      GO TO 490
C
C     PENALTY FUNCTION I.
C
  190 CONTINUE
      T1 = -CP25
      DO 200 J = 1, N
         T1 = T1 + X(J)**2
  200    CONTINUE
      D1 = TWO*AP
      TH = FOUR*BP*T1
      DO 210 J = 1, N
         G(J) = D1*(X(J) - ONE) + X(J)*TH
  210    CONTINUE
      GO TO 490
C
C     PENALTY FUNCTION II.
C
  220 CONTINUE
      T1 = -ONE
      DO 230 J = 1, N
         T1 = T1 + FLOAT(N-J+1)*X(J)**2
  230    CONTINUE
      D1 = EXP(CP1)
      D2 = ONE
      TH = FOUR*BP*T1
      DO 250 J = 1, N
         G(J) = FLOAT(N-J+1)*X(J)*TH
         S1 = EXP(X(J)/TEN)
         IF (J .EQ. 1) GO TO 240
         S3 = S1 + S2 - D2*(D1 + ONE)
         G(J) = G(J) + AP*S1*(S3 + S1 - ONE/D1)/FIVE
         G(J-1) = G(J-1) + AP*S2*S3/FIVE
  240    CONTINUE
         S2 = S1
         D2 = D1*D2
  250    CONTINUE
      G(1) = G(1) + TWO*BP*(X(1) - CP2)
      GO TO 490
C
C     BROWN BADLY SCALED FUNCTION.
C
  260 CONTINUE
      T1 = X(1) - C1PD6
      T2 = X(2) - C2PDM6
      T3 = X(1)*X(2) - TWO
      G(1) = TWO*(T1 + X(2)*T3)
      G(2) = TWO*(T2 + X(1)*T3)
      GO TO 490
C
C     BROWN AND DENNIS FUNCTION.
C
  270 CONTINUE
      G(1) = ZERO
      G(2) = ZERO
      G(3) = ZERO
      G(4) = ZERO
      DO 280 I = 1, 20
         D1 = FLOAT(I)/FIVE
         D2 = SIN(D1)
         T1 = X(1) + D1*X(2) - EXP(D1)
         T2 = X(3) + D2*X(4) - COS(D1)
         T = T1**2 + T2**2
         S1 = T1*T
         S2 = T2*T
         G(1) = G(1) + S1
         G(2) = G(2) + D1*S1
         G(3) = G(3) + S2
         G(4) = G(4) + D2*S2
  280    CONTINUE
      G(1) = FOUR*G(1)
      G(2) = FOUR*G(2)
      G(3) = FOUR*G(3)
      G(4) = FOUR*G(4)
      GO TO 490
C
C     GULF RESEARCH AND DEVELOPMENT FUNCTION.
C
  290 CONTINUE
      G(1) = ZERO
      G(2) = ZERO
      G(3) = ZERO
      D1 = TWO/THREE
      DO 300 I = 1, 99
         ARG = FLOAT(I)/C100
         R = ABS((-FIFTY*ALOG(ARG))**D1+C25-X(2))
         T1 = R**X(3)/X(1)
         T2 = EXP(-T1)
         T = T2 - ARG
         S1 = T1*T2*T
         G(1) = G(1) + S1
         G(2) = G(2) + S1/R
         G(3) = G(3) - S1*ALOG(R)
  300    CONTINUE
      G(1) = TWO*G(1)/X(1)
      G(2) = TWO*X(3)*G(2)
      G(3) = TWO*G(3)
      GO TO 490
C
C     TRIGONOMETRIC FUNCTION.
C
  310 CONTINUE
      S1 = ZERO
      DO 320 J = 1, N
         G(J) = COS(X(J))
         S1 = S1 + G(J)
  320    CONTINUE
      S2 = ZERO
      DO 330 J = 1, N
         TH = SIN(X(J))
         T = FLOAT(N+J) - TH - S1 - FLOAT(J)*G(J)
         S2 = S2 + T
         G(J) = (FLOAT(J)*TH - G(J))*T
  330    CONTINUE
      DO 340 J = 1, N
         G(J) = TWO*(G(J) + SIN(X(J))*S2)
  340    CONTINUE
      GO TO 490
C
C     EXTENDED ROSENBROCK FUNCTION.
C
  350 CONTINUE
      DO 360 J = 1, N, 2
         T1 = ONE - X(J)
         G(J+1) = C200*(X(J+1) - X(J)**2)
         G(J) = -TWO*(X(J)*G(J+1) + T1)
  360    CONTINUE
      GO TO 490
C
C     EXTENDED POWELL FUNCTION.
C
  370 CONTINUE
      DO 380 J = 1, N, 4
         T = X(J) + TEN*X(J+1)
         T1 = X(J+2) - X(J+3)
         S1 = FIVE*T1
         T2 = X(J+1) - TWO*X(J+2)
         S2 = FOUR*T2**3
         T3 = X(J) - X(J+3)
         S3 = TWENTY*T3**3
         G(J) = TWO*(T + S3)
         G(J+1) = TWENTY*T + S2
         G(J+2) = TWO*(S1 - S2)
         G(J+3) = -TWO*(S1 + S3)
  380    CONTINUE
      GO TO 490
C
C     BEALE FUNCTION.
C
  390 CONTINUE
      S1 = ONE - X(2)
      T1 = C1P5 - X(1)*S1
      S2 = ONE - X(2)**2
      T2 = C2P25 - X(1)*S2
      S3 = ONE - X(2)**3
      T3 = C2P625 - X(1)*S3
      G(1) = -TWO*(S1*T1 + S2*T2 + S3*T3)
      G(2) = TWO*X(1)*(T1 + X(2)*(TWO*T2 + THREE*X(2)*T3))
      GO TO 490
C
C     WOOD FUNCTION.
C
  400 CONTINUE
      S1 = X(2) - X(1)**2
      S2 = ONE - X(1)
      S3 = X(2) - ONE
      T1 = X(4) - X(3)**2
      T2 = ONE - X(3)
      T3 = X(4) - ONE
      G(1) = -TWO*(C200*X(1)*S1 + S2)
      G(2) = C200*S1 + C20P2*S3 + C19P8*T3
      G(3) = -TWO*(C180*X(3)*T1 + T2)
      G(4) = C180*T1 + C20P2*T3 + C19P8*S3
      GO TO 490
C
C     CHEBYQUAD FUNCTION.
C
  410 CONTINUE
      DO 420 I = 1, N
         FVEC(I) = ZERO
  420    CONTINUE
      DO 440 J = 1, N
         T1 = ONE
         T2 = TWO*X(J) - ONE
         T = TWO*T2
         DO 430 I = 1, N
            FVEC(I) = FVEC(I) + T2
            TH = T*T2 - T1
            T1 = T2
            T2 = TH
  430       CONTINUE
  440    CONTINUE
      D1 = ONE/FLOAT(N)
      IEV = -1
      DO 450 I = 1, N
         FVEC(I) = D1*FVEC(I)
         IF (IEV .GT. 0) FVEC(I) = FVEC(I) + ONE/(FLOAT(I)**2 - ONE)
         IEV = -IEV
  450    CONTINUE
      DO 470 J = 1, N
         G(J) = ZERO
         T1 = ONE
         T2 = TWO*X(J) - ONE
         T = TWO*T2
         S1 = ZERO
         S2 = TWO
         DO 460 I = 1, N
            G(J) = G(J) + FVEC(I)*S2
            TH = FOUR*T2 + T*S2 - S1
            S1 = S2
            S2 = TH
            TH = T*T2 - T1
            T1 = T2
            T2 = TH
  460       CONTINUE
  470    CONTINUE
      D2 = TWO*D1
      DO 480 J = 1, N
         G(J) = D2*G(J)
  480    CONTINUE
  490 CONTINUE
      RETURN
C
C     LAST CARD OF SUBROUTINE GRDFCN.
C
      END
C ===== 8. SAMPLE DRIVER FOR DOUBLE PRECISION NONLINEAR EQUATIONS.
C     **********                                                        00000010
C                                                                       00000020
C     THIS PROGRAM TESTS CODES FOR THE SOLUTION OF N NONLINEAR          00000030
C     EQUATIONS IN N VARIABLES. IT CONSISTS OF A DRIVER AND AN          00000040
C     INTERFACE SUBROUTINE FCN. THE DRIVER READS IN DATA, CALLS THE     00000050
C     NONLINEAR EQUATION SOLVER, AND FINALLY PRINTS OUT INFORMATION     00000060
C     ON THE PERFORMANCE OF THE SOLVER. THIS IS ONLY A SAMPLE DRIVER,   00000070
C     MANY OTHER DRIVERS ARE POSSIBLE. THE INTERFACE SUBROUTINE FCN     00000080
C     IS NECESSARY TO TAKE INTO ACCOUNT THE FORMS OF CALLING            00000090
C     SEQUENCES USED BY THE FUNCTION SUBROUTINES IN THE VARIOUS         00000100
C     NONLINEAR EQUATION SOLVERS.                                       00000110
C                                                                       00000120
C     SUBPROGRAMS CALLED                                                00000130
C                                                                       00000140
C       USER-SUPPLIED ...... ENORM,FCN,SOLVER                           00000150
C                                                                       00000160
C       MINPACK-SUPPLIED ... INITPT,VECFCN                              00000170
C                                                                       00000180
C       FORTRAN-SUPPLIED ... DSQRT                                      00000190
C                                                                       00000200
C     MINPACK. VERSION OF NOVEMBER 1978.                                00000210
C     BURTON S. GARBOW, KENNETH E. HILLSTROM, JORGE J. MORE             00000220
C                                                                       00000230
C     **********                                                        00000240
      INTEGER I,IC,INFO,K,LWA,N,NFEV,NPROB,NREAD,NTRIES,NWRITE          00000250
      INTEGER NA(60),NF(60),NP(60),NX(60)                               00000260
      DOUBLE PRECISION FACTOR,FNORM1,FNORM2,ONE,TEN,TOL                 00000270
      DOUBLE PRECISION FNM(60),FVEC(40),WA(2660),X(40)                  00000280
      DOUBLE PRECISION ENORM                                            00000290
      EXTERNAL FCN                                                      00000300
      COMMON /REFNUM/ NPROB,NFEV                                        00000310
C                                                                       00000320
C     LOGICAL INPUT UNIT IS ASSUMED TO BE NUMBER 5.                     00000330
C     LOGICAL OUTPUT UNIT IS ASSUMED TO BE NUMBER 6.                    00000340
C                                                                       00000350
      DATA NREAD,NWRITE /5,6/                                           00000360
C                                                                       00000370
      DATA ONE,TEN,TOL /1.0D0,1.0D1,1.D-10/                             00000380
      LWA = 2660                                                        00000390
      IC = 0                                                            00000400
   10 CONTINUE                                                          00000410
         READ (NREAD,50) NPROB,N,NTRIES                                 00000420
         IF (NPROB .LE. 0) GO TO 30                                     00000430
         FACTOR = ONE                                                   00000440
         DO 20 K = 1, NTRIES                                            00000450
            IC = IC + 1                                                 00000460
            CALL INITPT(N,X,NPROB,FACTOR)                               00000470
            CALL VECFCN(N,X,FVEC,NPROB)                                 00000480
            FNORM1 = ENORM(N,FVEC)                                      00000490
            WRITE (NWRITE,60) NPROB,N                                   00000500
            NFEV = 0                                                    00000510
            CALL SOLVER(FCN,N,X,FVEC,TOL,INFO,WA,LWA)                   00000520
            FNORM2 = ENORM(N,FVEC)                                      00000530
            NP(IC) = NPROB                                              00000540
            NA(IC) = N                                                  00000550
            NF(IC) = NFEV                                               00000560
            NX(IC) = INFO                                               00000570
            FNM(IC) = FNORM2                                            00000580
            WRITE (NWRITE,70) FNORM1,FNORM2,NFEV,INFO,(X(I), I = 1, N)  00000590
            FACTOR = TEN*FACTOR                                         00000600
   20       CONTINUE                                                    00000610
         GO TO 10                                                       00000620
   30 CONTINUE                                                          00000630
      WRITE (NWRITE,80) IC                                              00000640
      WRITE (NWRITE,90)                                                 00000650
      DO 40 I = 1, IC                                                   00000660
         WRITE (NWRITE,100) NP(I),NA(I),NF(I),NX(I),FNM(I)              00000670
   40    CONTINUE                                                       00000680
      STOP                                                              00000690
   50 FORMAT (3I5)                                                      00000700
   60 FORMAT ( //// 5X, 8H PROBLEM, I5, 5X, 10H DIMENSION, I5, 5X //)   00000710
   70 FORMAT (5X, 33H INITIAL L2 NORM OF THE RESIDUALS, D15.7 // 5X,    00000720
     *        33H FINAL L2 NORM OF THE RESIDUALS  , D15.7 // 5X,        00000730
     *        33H NUMBER OF FUNCTION EVALUATIONS  , I10 // 5X,          00000740
     *        15H EXIT PARAMETER, 18X, I10 // 5X,                       00000750
     *        27H FINAL APPROXIMATE SOLUTION // (5X, 5D15.7))           00000760
   80 FORMAT (12H1SUMMARY OF , I3, 16H CALLS TO HYBRD1 /)               00000770
   90 FORMAT (39H NPROB   N    NFEV  INFO  FINAL L2 NORM /)             00000780
  100 FORMAT (I4, I6, I7, I6, 1X, D15.7)                                00000790
C                                                                       00000800
C     LAST CARD OF DRIVER.                                              00000810
C                                                                       00000820
      END                                                               00000830
      SUBROUTINE FCN(N,X,FVEC,IFLAG)                                    00000840
      INTEGER N,IFLAG
      DOUBLE PRECISION X(N),FVEC(N)
C     **********
C
C     THE CALLING SEQUENCE OF FCN SHOULD BE IDENTICAL TO THE
C     CALLING SEQUENCE OF THE FUNCTION SUBROUTINE IN THE NONLINEAR
C     EQUATION SOLVER. FCN SHOULD ONLY CALL THE TESTING FUNCTION
C     SUBROUTINE VECFCN WITH THE APPROPRIATE VALUE OF PROBLEM
C     NUMBER (NPROB).
C
C     SUBPROGRAMS CALLED
C
C       MINPACK-SUPPLIED ... VECFCN
C
C     MINPACK. VERSION OF JULY 1978.
C     BURTON S. GARBOW, KENNETH E. HILLSTROM, JORGE J. MORE
C
C     **********
      INTEGER NPROB,NFEV
      COMMON /REFNUM/ NPROB,NFEV
      CALL VECFCN(N,X,FVEC,NPROB)
      NFEV = NFEV + 1
      RETURN
C
C     LAST CARD OF INTERFACE SUBROUTINE FCN.
C
      END
C ===== 9. SAMPLE DRIVER FOR SINGLE PRECISION NONLINEAR EQUATIONS.
C     **********                                                        00000010
C                                                                       00000020
C     THIS PROGRAM TESTS CODES FOR THE SOLUTION OF N NONLINEAR          00000030
C     EQUATIONS IN N VARIABLES. IT CONSISTS OF A DRIVER AND AN          00000040
C     INTERFACE SUBROUTINE FCN. THE DRIVER READS IN DATA, CALLS THE     00000050
C     NONLINEAR EQUATION SOLVER, AND FINALLY PRINTS OUT INFORMATION     00000060
C     ON THE PERFORMANCE OF THE SOLVER. THIS IS ONLY A SAMPLE DRIVER,   00000070
C     MANY OTHER DRIVERS ARE POSSIBLE. THE INTERFACE SUBROUTINE FCN     00000080
C     IS NECESSARY TO TAKE INTO ACCOUNT THE FORMS OF CALLING            00000090
C     SEQUENCES USED BY THE FUNCTION SUBROUTINES IN THE VARIOUS         00000100
C     NONLINEAR EQUATION SOLVERS.                                       00000110
C                                                                       00000120
C     SUBPROGRAMS CALLED                                                00000130
C                                                                       00000140
C       USER-SUPPLIED ...... ENORM,FCN,SOLVER                           00000150
C                                                                       00000160
C       MINPACK-SUPPLIED ... INITPT,VECFCN                              00000170
C                                                                       00000180
C       FORTRAN-SUPPLIED ... SQRT                                       00000190
C                                                                       00000200
C     MINPACK. VERSION OF NOVEMBER 1978.                                00000210
C     BURTON S. GARBOW, KENNETH E. HILLSTROM, JORGE J. MORE             00000220
C                                                                       00000230
C     **********                                                        00000240
      INTEGER I,IC,INFO,K,LWA,N,NFEV,NPROB,NREAD,NTRIES,NWRITE          00000250
      INTEGER NA(60),NF(60),NP(60),NX(60)                               00000260
      REAL FACTOR,FNORM1,FNORM2,ONE,TEN,TOL                             00000270
      REAL FNM(60),FVEC(40),WA(2660),X(40)                              00000280
      REAL ENORM                                                        00000290
      EXTERNAL FCN                                                      00000300
      COMMON /REFNUM/ NPROB,NFEV                                        00000310
C                                                                       00000320
C     LOGICAL INPUT UNIT IS ASSUMED TO BE NUMBER 5.                     00000330
C     LOGICAL OUTPUT UNIT IS ASSUMED TO BE NUMBER 6.                    00000340
C                                                                       00000350
      DATA NREAD,NWRITE /5,6/                                           00000360
C                                                                       00000370
      DATA ONE,TEN,TOL /1.0E0,1.0E1,1.E-5/                              00000380
      LWA = 2660                                                        00000390
      IC = 0                                                            00000400
   10 CONTINUE                                                          00000410
         READ (NREAD,50) NPROB,N,NTRIES                                 00000420
         IF (NPROB .LE. 0) GO TO 30                                     00000430
         FACTOR = ONE                                                   00000440
         DO 20 K = 1, NTRIES                                            00000450
            IC = IC + 1                                                 00000460
            CALL INITPT(N,X,NPROB,FACTOR)                               00000470
            CALL VECFCN(N,X,FVEC,NPROB)                                 00000480
            FNORM1 = ENORM(N,FVEC)                                      00000490
            WRITE (NWRITE,60) NPROB,N                                   00000500
            NFEV = 0                                                    00000510
            CALL SOLVER(FCN,N,X,FVEC,TOL,INFO,WA,LWA)                   00000520
            FNORM2 = ENORM(N,FVEC)                                      00000530
            NP(IC) = NPROB                                              00000540
            NA(IC) = N                                                  00000550
            NF(IC) = NFEV                                               00000560
            NX(IC) = INFO                                               00000570
            FNM(IC) = FNORM2                                            00000580
            WRITE (NWRITE,70) FNORM1,FNORM2,NFEV,INFO,(X(I), I = 1, N)  00000590
            FACTOR = TEN*FACTOR                                         00000600
   20       CONTINUE                                                    00000610
         GO TO 10                                                       00000620
   30 CONTINUE                                                          00000630
      WRITE (NWRITE,80) IC                                              00000640
      WRITE (NWRITE,90)                                                 00000650
      DO 40 I = 1, IC                                                   00000660
         WRITE (NWRITE,100) NP(I),NA(I),NF(I),NX(I),FNM(I)              00000670
   40    CONTINUE                                                       00000680
      STOP                                                              00000690
   50 FORMAT (3I5)                                                      00000700
   60 FORMAT ( //// 5X, 8H PROBLEM, I5, 5X, 10H DIMENSION, I5, 5X //)   00000710
   70 FORMAT (5X, 33H INITIAL L2 NORM OF THE RESIDUALS, E15.7 // 5X,    00000720
     *        33H FINAL L2 NORM OF THE RESIDUALS  , E15.7 // 5X,        00000730
     *        33H NUMBER OF FUNCTION EVALUATIONS  , I10 // 5X,          00000740
     *        15H EXIT PARAMETER, 18X, I10 // 5X,                       00000750
     *        27H FINAL APPROXIMATE SOLUTION // (5X, 5E15.7))           00000760
   80 FORMAT (12H1SUMMARY OF , I3, 16H CALLS TO HYBRD1 /)               00000770
   90 FORMAT (39H NPROB   N    NFEV  INFO  FINAL L2 NORM /)             00000780
  100 FORMAT (I4, I6, I7, I6, 1X, E15.7)                                00000790
C                                                                       00000800
C     LAST CARD OF DRIVER.                                              00000810
C                                                                       00000820
      END                                                               00000830
      SUBROUTINE FCN(N,X,FVEC,IFLAG)                                    00000840
      INTEGER N,IFLAG
      REAL X(N),FVEC(N)
C     **********
C
C     THE CALLING SEQUENCE OF FCN SHOULD BE IDENTICAL TO THE
C     CALLING SEQUENCE OF THE FUNCTION SUBROUTINE IN THE NONLINEAR
C     EQUATION SOLVER. FCN SHOULD ONLY CALL THE TESTING FUNCTION
C     SUBROUTINE VECFCN WITH THE APPROPRIATE VALUE OF PROBLEM
C     NUMBER (NPROB).
C
C     SUBPROGRAMS CALLED
C
C       MINPACK-SUPPLIED ... VECFCN
C
C     MINPACK. VERSION OF JULY 1978.
C     BURTON S. GARBOW, KENNETH E. HILLSTROM, JORGE J. MORE
C
C     **********
      INTEGER NPROB,NFEV
      COMMON /REFNUM/ NPROB,NFEV
      CALL VECFCN(N,X,FVEC,NPROB)
      NFEV = NFEV + 1
      RETURN
C
C     LAST CARD OF INTERFACE SUBROUTINE FCN.
C
      END
C ===== 10. SAMPLE DRIVER FOR DOUBLE PRECISION NONLINEAR LEAST-SQUARES.
C     **********                                                        00000010
C                                                                       00000020
C     THIS PROGRAM TESTS CODES FOR THE LEAST-SQUARES SOLUTION OF        00000030
C     M NONLINEAR EQUATIONS IN N VARIABLES. IT CONSISTS OF A DRIVER     00000040
C     AND AN INTERFACE SUBROUTINE FCN. THE DRIVER READS IN DATA,        00000050
C     CALLS THE NONLINEAR LEAST-SQUARES SOLVER, AND FINALLY PRINTS      00000060
C     OUT INFORMATION ON THE PERFORMANCE OF THE SOLVER. THIS IS         00000070
C     ONLY A SAMPLE DRIVER, MANY OTHER DRIVERS ARE POSSIBLE. THE        00000080
C     INTERFACE SUBROUTINE FCN IS NECESSARY TO TAKE INTO ACCOUNT THE    00000090
C     FORMS OF CALLING SEQUENCES USED BY THE FUNCTION AND JACOBIAN      00000100
C     SUBROUTINES IN THE VARIOUS NONLINEAR LEAST-SQUARES SOLVERS.       00000110
C                                                                       00000120
C     SUBPROGRAMS CALLED                                                00000130
C                                                                       00000140
C       USER-SUPPLIED ...... ENORM,FCN,SOLVER                           00000150
C                                                                       00000160
C       MINPACK-SUPPLIED ... INITPT,SSQFCN                              00000170
C                                                                       00000180
C       FORTRAN-SUPPLIED ... DSQRT                                      00000190
C                                                                       00000200
C     MINPACK. VERSION OF NOVEMBER 1978.                                00000210
C     BURTON S. GARBOW, KENNETH E. HILLSTROM, JORGE J. MORE             00000220
C                                                                       00000230
C     **********                                                        00000240
      INTEGER I,IC,INFO,K,LDFJAC,LWA,M,N,NFEV,NJEV,NPROB,NREAD,NTRIES,  00000250
     *        NWRITE                                                    00000260
      INTEGER IWA(40),MA(60),NA(60),NF(60),NJ(60),NP(60),NX(60)         00000270
      DOUBLE PRECISION FACTOR,FNORM1,FNORM2,ONE,TEN,TOL                 00000280
      DOUBLE PRECISION FJAC(65,40),FNM(60),FVEC(65),WA(265),X(40)       00000290
      DOUBLE PRECISION ENORM                                            00000300
      EXTERNAL FCN                                                      00000310
      COMMON /REFNUM/ NPROB,NFEV,NJEV                                   00000320
C                                                                       00000330
C     LOGICAL INPUT UNIT IS ASSUMED TO BE NUMBER 5.                     00000340
C     LOGICAL OUTPUT UNIT IS ASSUMED TO BE NUMBER 6.                    00000350
C                                                                       00000360
      DATA NREAD,NWRITE /5,6/                                           00000370
C                                                                       00000380
      DATA ONE,TEN,TOL /1.0D0,1.0D1,1.D-10/                             00000390
      LDFJAC = 65                                                       00000400
      LWA = 265                                                         00000410
      IC = 0                                                            00000420
   10 CONTINUE                                                          00000430
         READ (NREAD,50) NPROB,N,M,NTRIES                               00000440
         IF (NPROB .LE. 0) GO TO 30                                     00000450
         FACTOR = ONE                                                   00000460
         DO 20 K = 1, NTRIES                                            00000470
            IC = IC + 1                                                 00000480
            CALL INITPT(N,X,NPROB,FACTOR)                               00000490
            CALL SSQFCN(M,N,X,FVEC,NPROB)                               00000500
            FNORM1 = ENORM(M,FVEC)                                      00000510
            WRITE (NWRITE,60) NPROB,N,M                                 00000520
            NFEV = 0                                                    00000530
            NJEV = 0                                                    00000540
            CALL SOLVER(FCN,M,N,X,FVEC,FJAC,LDFJAC,TOL,INFO,IWA,WA,     00000550
     *                  LWA)                                            00000560
            CALL SSQFCN(M,N,X,FVEC,NPROB)                               00000570
            FNORM2 = ENORM(M,FVEC)                                      00000580
            NP(IC) = NPROB                                              00000590
            NA(IC) = N                                                  00000600
            MA(IC) = M                                                  00000610
            NF(IC) = NFEV                                               00000620
            NJ(IC) = NJEV                                               00000630
            NX(IC) = INFO                                               00000640
            FNM(IC) = FNORM2                                            00000650
            WRITE (NWRITE,70)                                           00000660
     *            FNORM1,FNORM2,NFEV,NJEV,INFO,(X(I), I = 1, N)         00000670
            FACTOR = TEN*FACTOR                                         00000680
   20       CONTINUE                                                    00000690
         GO TO 10                                                       00000700
   30 CONTINUE                                                          00000710
      WRITE (NWRITE,80) IC                                              00000720
      WRITE (NWRITE,90)                                                 00000730
      DO 40 I = 1, IC                                                   00000740
         WRITE (NWRITE,100) NP(I),NA(I),MA(I),NF(I),NJ(I),NX(I),FNM(I)  00000750
   40    CONTINUE                                                       00000760
      STOP                                                              00000770
   50 FORMAT (4I5)                                                      00000780
   60 FORMAT ( //// 5X, 8H PROBLEM, I5, 5X, 11H DIMENSIONS, 2I5, 5X //  00000790
     *         )                                                        00000800
   70 FORMAT (5X, 33H INITIAL L2 NORM OF THE RESIDUALS, D15.7 // 5X,    00000810
     *        33H FINAL L2 NORM OF THE RESIDUALS  , D15.7 // 5X,        00000820
     *        33H NUMBER OF FUNCTION EVALUATIONS  , I10 // 5X,          00000830
     *        33H NUMBER OF JACOBIAN EVALUATIONS  , I10 // 5X,          00000840
     *        15H EXIT PARAMETER, 18X, I10 // 5X,                       00000850
     *        27H FINAL APPROXIMATE SOLUTION // (5X, 5D15.7))           00000860
   80 FORMAT (12H1SUMMARY OF , I3, 16H CALLS TO LMDER1 /)               00000870
   90 FORMAT (49H NPROB   N    M   NFEV  NJEV  INFO  FINAL L2 NORM /)   00000880
  100 FORMAT (3I5, 3I6, 2X, D15.7)                                      00000890
C                                                                       00000900
C     LAST CARD OF DRIVER.                                              00000910
C                                                                       00000920
      END                                                               00000930
      SUBROUTINE FCN(M,N,X,FVEC,FJAC,LDFJAC,IFLAG)                      00000940
      INTEGER M,N,LDFJAC,IFLAG
      DOUBLE PRECISION X(N),FVEC(M),FJAC(LDFJAC,N)
C     **********
C
C     THE CALLING SEQUENCE OF FCN SHOULD BE IDENTICAL TO THE
C     CALLING SEQUENCE OF THE FUNCTION SUBROUTINE IN THE NONLINEAR
C     LEAST-SQUARES SOLVER. FCN SHOULD ONLY CALL THE TESTING
C     FUNCTION AND JACOBIAN SUBROUTINES SSQFCN AND SSQJAC WITH
C     THE APPROPRIATE VALUE OF PROBLEM NUMBER (NPROB).
C
C     SUBPROGRAMS CALLED
C
C       MINPACK-SUPPLIED ... SSQFCN,SSQJAC
C
C     MINPACK. VERSION OF JULY 1978.
C     BURTON S. GARBOW, KENNETH E. HILLSTROM, JORGE J. MORE
C
C     **********
      INTEGER NPROB,NFEV,NJEV
      COMMON /REFNUM/ NPROB,NFEV,NJEV
      IF (IFLAG .EQ. 1) CALL SSQFCN(M,N,X,FVEC,NPROB)
      IF (IFLAG .EQ. 2) CALL SSQJAC(M,N,X,FJAC,LDFJAC,NPROB)
      IF (IFLAG .EQ. 1) NFEV = NFEV + 1
      IF (IFLAG .EQ. 2) NJEV = NJEV + 1
      RETURN
C
C     LAST CARD OF INTERFACE SUBROUTINE FCN.
C
      END
C ===== 11. SAMPLE DRIVER FOR SINGLE PRECISION NONLINEAR LEAST-SQUARES.
C     **********                                                        00000010
C                                                                       00000020
C     THIS PROGRAM TESTS CODES FOR THE LEAST-SQUARES SOLUTION OF        00000030
C     M NONLINEAR EQUATIONS IN N VARIABLES. IT CONSISTS OF A DRIVER     00000040
C     AND AN INTERFACE SUBROUTINE FCN. THE DRIVER READS IN DATA,        00000050
C     CALLS THE NONLINEAR LEAST-SQUARES SOLVER, AND FINALLY PRINTS      00000060
C     OUT INFORMATION ON THE PERFORMANCE OF THE SOLVER. THIS IS         00000070
C     ONLY A SAMPLE DRIVER, MANY OTHER DRIVERS ARE POSSIBLE. THE        00000080
C     INTERFACE SUBROUTINE FCN IS NECESSARY TO TAKE INTO ACCOUNT THE    00000090
C     FORMS OF CALLING SEQUENCES USED BY THE FUNCTION AND JACOBIAN      00000100
C     SUBROUTINES IN THE VARIOUS NONLINEAR LEAST-SQUARES SOLVERS.       00000110
C                                                                       00000120
C     SUBPROGRAMS CALLED                                                00000130
C                                                                       00000140
C       USER-SUPPLIED ...... ENORM,FCN,SOLVER                           00000150
C                                                                       00000160
C       MINPACK-SUPPLIED ... INITPT,SSQFCN                              00000170
C                                                                       00000180
C       FORTRAN-SUPPLIED ... SQRT                                       00000190
C                                                                       00000200
C     MINPACK. VERSION OF NOVEMBER 1978.                                00000210
C     BURTON S. GARBOW, KENNETH E. HILLSTROM, JORGE J. MORE             00000220
C                                                                       00000230
C     **********                                                        00000240
      INTEGER I,IC,INFO,K,LDFJAC,LWA,M,N,NFEV,NJEV,NPROB,NREAD,NTRIES,  00000250
     *        NWRITE                                                    00000260
      INTEGER IWA(40),MA(60),NA(60),NF(60),NJ(60),NP(60),NX(60)         00000270
      REAL FACTOR,FNORM1,FNORM2,ONE,TEN,TOL                             00000280
      REAL FJAC(65,40),FNM(60),FVEC(65),WA(265),X(40)                   00000290
      REAL ENORM                                                        00000300
      EXTERNAL FCN                                                      00000310
      COMMON /REFNUM/ NPROB,NFEV,NJEV                                   00000320
C                                                                       00000330
C     LOGICAL INPUT UNIT IS ASSUMED TO BE NUMBER 5.                     00000340
C     LOGICAL OUTPUT UNIT IS ASSUMED TO BE NUMBER 6.                    00000350
C                                                                       00000360
      DATA NREAD,NWRITE /5,6/                                           00000370
C                                                                       00000380
      DATA ONE,TEN,TOL /1.0E0,1.0E1,1.E-5/                              00000390
      LDFJAC = 65                                                       00000400
      LWA = 265                                                         00000410
      IC = 0                                                            00000420
   10 CONTINUE                                                          00000430
         READ (NREAD,50) NPROB,N,M,NTRIES                               00000440
         IF (NPROB .LE. 0) GO TO 30                                     00000450
         FACTOR = ONE                                                   00000460
         DO 20 K = 1, NTRIES                                            00000470
            IC = IC + 1                                                 00000480
            CALL INITPT(N,X,NPROB,FACTOR)                               00000490
            CALL SSQFCN(M,N,X,FVEC,NPROB)                               00000500
            FNORM1 = ENORM(M,FVEC)                                      00000510
            WRITE (NWRITE,60) NPROB,N,M                                 00000520
            NFEV = 0                                                    00000530
            NJEV = 0                                                    00000540
            CALL SOLVER(FCN,M,N,X,FVEC,FJAC,LDFJAC,TOL,INFO,IWA,WA,     00000550
     *                  LWA)                                            00000560
            CALL SSQFCN(M,N,X,FVEC,NPROB)                               00000570
            FNORM2 = ENORM(M,FVEC)                                      00000580
            NP(IC) = NPROB                                              00000590
            NA(IC) = N                                                  00000600
            MA(IC) = M                                                  00000610
            NF(IC) = NFEV                                               00000620
            NJ(IC) = NJEV                                               00000630
            NX(IC) = INFO                                               00000640
            FNM(IC) = FNORM2                                            00000650
            WRITE (NWRITE,70)                                           00000660
     *            FNORM1,FNORM2,NFEV,NJEV,INFO,(X(I), I = 1, N)         00000670
            FACTOR = TEN*FACTOR                                         00000680
   20       CONTINUE                                                    00000690
         GO TO 10                                                       00000700
   30 CONTINUE                                                          00000710
      WRITE (NWRITE,80) IC                                              00000720
      WRITE (NWRITE,90)                                                 00000730
      DO 40 I = 1, IC                                                   00000740
         WRITE (NWRITE,100) NP(I),NA(I),MA(I),NF(I),NJ(I),NX(I),FNM(I)  00000750
   40    CONTINUE                                                       00000760
      STOP                                                              00000770
   50 FORMAT (4I5)                                                      00000780
   60 FORMAT ( //// 5X, 8H PROBLEM, I5, 5X, 11H DIMENSIONS, 2I5, 5X //  00000790
     *         )                                                        00000800
   70 FORMAT (5X, 33H INITIAL L2 NORM OF THE RESIDUALS, E15.7 // 5X,    00000810
     *        33H FINAL L2 NORM OF THE RESIDUALS  , E15.7 // 5X,        00000820
     *        33H NUMBER OF FUNCTION EVALUATIONS  , I10 // 5X,          00000830
     *        33H NUMBER OF JACOBIAN EVALUATIONS  , I10 // 5X,          00000840
     *        15H EXIT PARAMETER, 18X, I10 // 5X,                       00000850
     *        27H FINAL APPROXIMATE SOLUTION // (5X, 5E15.7))           00000860
   80 FORMAT (12H1SUMMARY OF , I3, 16H CALLS TO LMDER1 /)               00000870
   90 FORMAT (49H NPROB   N    M   NFEV  NJEV  INFO  FINAL L2 NORM /)   00000880
  100 FORMAT (3I5, 3I6, 2X, E15.7)                                      00000890
C                                                                       00000900
C     LAST CARD OF DRIVER.                                              00000910
C                                                                       00000920
      END                                                               00000930
      SUBROUTINE FCN(M,N,X,FVEC,FJAC,LDFJAC,IFLAG)                      00000940
      INTEGER M,N,LDFJAC,IFLAG
      REAL X(N),FVEC(M),FJAC(LDFJAC,N)
C     **********
C
C     THE CALLING SEQUENCE OF FCN SHOULD BE IDENTICAL TO THE
C     CALLING SEQUENCE OF THE FUNCTION SUBROUTINE IN THE NONLINEAR
C     LEAST-SQUARES SOLVER. FCN SHOULD ONLY CALL THE TESTING
C     FUNCTION AND JACOBIAN SUBROUTINES SSQFCN AND SSQJAC WITH
C     THE APPROPRIATE VALUE OF PROBLEM NUMBER (NPROB).
C
C     SUBPROGRAMS CALLED
C
C       MINPACK-SUPPLIED ... SSQFCN,SSQJAC
C
C     MINPACK. VERSION OF JULY 1978.
C     BURTON S. GARBOW, KENNETH E. HILLSTROM, JORGE J. MORE
C
C     **********
      INTEGER NPROB,NFEV,NJEV
      COMMON /REFNUM/ NPROB,NFEV,NJEV
      IF (IFLAG .EQ. 1) CALL SSQFCN(M,N,X,FVEC,NPROB)
      IF (IFLAG .EQ. 2) CALL SSQJAC(M,N,X,FJAC,LDFJAC,NPROB)
      IF (IFLAG .EQ. 1) NFEV = NFEV + 1
      IF (IFLAG .EQ. 2) NJEV = NJEV + 1
      RETURN
C
C     LAST CARD OF INTERFACE SUBROUTINE FCN.
C
      END
C ===== 12. SAMPLE DRIVER FOR DOUBLE PRECISION UNCONSTRAINED NONLINEAR
C =====     MINIMIZATION.
C     **********                                                        00000010
C                                                                       00000020
C     THIS PROGRAM TESTS CODES FOR THE UNCONSTRAINED OPTIMIZATION OF    00000030
C     A NONLINEAR FUNCTION OF N VARIABLES. IT CONSISTS OF A DRIVER      00000040
C     AND AN INTERFACE SUBROUTINE FCN. THE DRIVER READS IN DATA,        00000050
C     CALLS THE UNCONSTRAINED OPTIMIZER, AND FINALLY PRINTS OUT         00000060
C     INFORMATION ON THE PERFORMANCE OF THE OPTIMIZER. THIS IS          00000070
C     ONLY A SAMPLE DRIVER, MANY OTHER DRIVERS ARE POSSIBLE. THE        00000080
C     INTERFACE SUBROUTINE FCN IS NECESSARY TO TAKE INTO ACCOUNT THE    00000090
C     FORMS OF CALLING SEQUENCES USED BY THE FUNCTION SUBROUTINES       00000100
C     IN THE VARIOUS UNCONSTRAINED OPTIMIZERS.                          00000110
C                                                                       00000120
C     SUBPROGRAMS CALLED                                                00000130
C                                                                       00000140
C       USER-SUPPLIED ...... ENORM,FCN,SOLVER                           00000150
C                                                                       00000160
C       MINPACK-SUPPLIED ... GRDFCN,INITPT,OBJFCN                       00000170
C                                                                       00000180
C       FORTRAN-SUPPLIED ... DSQRT                                      00000190
C                                                                       00000200
C     MINPACK. VERSION OF NOVEMBER 1978.                                00000210
C     BURTON S. GARBOW, KENNETH E. HILLSTROM, JORGE J. MORE             00000220
C                                                                       00000230
C     **********                                                        00000240
      INTEGER I,IC,INFO,K,LWA,N,NFEV,NPROB,NREAD,NTRIES,NWRITE          00000250
      INTEGER NA(120),NF(120),NP(120),NX(120)                           00000260
      DOUBLE PRECISION FACTOR,F1,F2,GNORM1,GNORM2,ONE,TEN,TOL           00000270
      DOUBLE PRECISION FVAL(120),GVEC(100),GNM(120),WA(6130),X(100)     00000280
      DOUBLE PRECISION ENORM                                            00000290
      EXTERNAL FCN                                                      00000300
      COMMON /REFNUM/ NPROB,NFEV                                        00000310
C                                                                       00000320
C     LOGICAL INPUT UNIT IS ASSUMED TO BE NUMBER 5.                     00000330
C     LOGICAL OUTPUT UNIT IS ASSUMED TO BE NUMBER 6.                    00000340
C                                                                       00000350
      DATA NREAD,NWRITE /5,6/                                           00000360
C                                                                       00000370
      DATA ONE,TEN,TOL /1.0D0,1.0D1,1.D-10/                             00000380
      LWA = 6130                                                        00000390
      IC = 0                                                            00000400
   10 CONTINUE                                                          00000410
         READ (NREAD,50) NPROB,N,NTRIES                                 00000420
         IF (NPROB .LE. 0) GO TO 30                                     00000430
         FACTOR = ONE                                                   00000440
         DO 20 K = 1, NTRIES                                            00000450
            IC = IC + 1                                                 00000460
            CALL INITPT(N,X,NPROB,FACTOR)                               00000470
            CALL OBJFCN(N,X,F1,NPROB)                                   00000480
            CALL GRDFCN(N,X,GVEC,NPROB)                                 00000490
            GNORM1 = ENORM(N,GVEC)                                      00000500
            WRITE (NWRITE,60) NPROB,N                                   00000510
            NFEV = 0                                                    00000520
            CALL SOLVER(FCN,N,X,F2,GVEC,TOL,INFO,WA,LWA)                00000530
            CALL OBJFCN(N,X,F2,NPROB)                                   00000540
            CALL GRDFCN(N,X,GVEC,NPROB)                                 00000550
            GNORM2 = ENORM(N,GVEC)                                      00000560
            NP(IC) = NPROB                                              00000570
            NA(IC) = N                                                  00000580
            NF(IC) = NFEV                                               00000590
            NX(IC) = INFO                                               00000600
            FVAL(IC) = F2                                               00000610
            GNM(IC) = GNORM2                                            00000620
            WRITE (NWRITE,70)                                           00000630
     *            F1,F2,GNORM1,GNORM2,NFEV,INFO,(X(I), I = 1, N)        00000640
            FACTOR = TEN*FACTOR                                         00000650
   20       CONTINUE                                                    00000660
         GO TO 10                                                       00000670
   30 CONTINUE                                                          00000680
      WRITE (NWRITE,80) IC                                              00000690
      WRITE (NWRITE,90)                                                 00000700
      DO 40 I = 1, IC                                                   00000710
         WRITE (NWRITE,100) NP(I),NA(I),NF(I),NX(I),FVAL(I),GNM(I)      00000720
   40    CONTINUE                                                       00000730
      STOP                                                              00000740
   50 FORMAT (3I5)                                                      00000750
   60 FORMAT ( //// 5X, 8H PROBLEM, I5, 5X, 10H DIMENSION, I5, 5X //)   00000760
   70 FORMAT (5X, 23H INITIAL FUNCTION VALUE, D15.7 // 5X,              00000770
     *        23H FINAL FUNCTION VALUE  , D15.7 // 5X,                  00000780
     *        23H INITIAL GRADIENT NORM , D15.7 // 5X,                  00000790
     *        23H FINAL GRADIENT NORM   , D15.7 // 5X,                  00000800
     *        33H NUMBER OF FUNCTION EVALUATIONS  , I10 // 5X,          00000810
     *        15H EXIT PARAMETER, 18X, I10 // 5X,                       00000820
     *        27H FINAL APPROXIMATE SOLUTION // (5X, 5D15.7))           00000830
   80 FORMAT (12H1SUMMARY OF , I3, 16H CALLS TO DRVCR1 /)               00000840
   90 FORMAT (25H NPROB   N    NFEV  INFO ,                             00000850
     *        42H FINAL FUNCTION VALUE  FINAL GRADIENT NORM /)          00000860
  100 FORMAT (I4, I6, I7, I6, 5X, D15.7, 6X, D15.7)                     00000870
C                                                                       00000880
C     LAST CARD OF DRIVER.                                              00000890
C                                                                       00000900
      END                                                               00000910
      SUBROUTINE FCN(N,X,F,GVEC,IFLAG)                                  00000920
      INTEGER N,IFLAG
      DOUBLE PRECISION F
      DOUBLE PRECISION X(N),GVEC(N)
C     **********
C
C     THE CALLING SEQUENCE OF FCN SHOULD BE IDENTICAL TO THE
C     CALLING SEQUENCE OF THE FUNCTION SUBROUTINE IN THE
C     UNCONSTRAINED OPTIMIZER. FCN SHOULD ONLY CALL THE TESTING
C     FUNCTION AND GRADIENT SUBROUTINES OBJFCN AND GRDFCN WITH
C     THE APPROPRIATE VALUE OF PROBLEM NUMBER (NPROB).
C
C     SUBPROGRAMS CALLED
C
C       MINPACK-SUPPLIED ... GRDFCN,OBJFCN
C
C     MINPACK. VERSION OF JULY 1978.
C     BURTON S. GARBOW, KENNETH E. HILLSTROM, JORGE J. MORE
C
C     **********
      INTEGER NPROB,NFEV
      COMMON /REFNUM/ NPROB,NFEV
      CALL OBJFCN(N,X,F,NPROB)
      CALL GRDFCN(N,X,GVEC,NPROB)
      NFEV = NFEV + 1
      RETURN
C
C     LAST CARD OF INTERFACE SUBROUTINE FCN.
C
      END
C ===== 13. SAMPLE DRIVER FOR SINGLE PRECISION UNCONSTRAINED NONLINEAR
C =====     MINIMIZATION.
C     **********                                                        00000010
C                                                                       00000020
C     THIS PROGRAM TESTS CODES FOR THE UNCONSTRAINED OPTIMIZATION OF    00000030
C     A NONLINEAR FUNCTION OF N VARIABLES. IT CONSISTS OF A DRIVER      00000040
C     AND AN INTERFACE SUBROUTINE FCN. THE DRIVER READS IN DATA,        00000050
C     CALLS THE UNCONSTRAINED OPTIMIZER, AND FINALLY PRINTS OUT         00000060
C     INFORMATION ON THE PERFORMANCE OF THE OPTIMIZER. THIS IS          00000070
C     ONLY A SAMPLE DRIVER, MANY OTHER DRIVERS ARE POSSIBLE. THE        00000080
C     INTERFACE SUBROUTINE FCN IS NECESSARY TO TAKE INTO ACCOUNT THE    00000090
C     FORMS OF CALLING SEQUENCES USED BY THE FUNCTION SUBROUTINES       00000100
C     IN THE VARIOUS UNCONSTRAINED OPTIMIZERS.                          00000110
C                                                                       00000120
C     SUBPROGRAMS CALLED                                                00000130
C                                                                       00000140
C       USER-SUPPLIED ...... ENORM,FCN,SOLVER                           00000150
C                                                                       00000160
C       MINPACK-SUPPLIED ... GRDFCN,INITPT,OBJFCN                       00000170
C                                                                       00000180
C       FORTRAN-SUPPLIED ... SQRT                                       00000190
C                                                                       00000200
C     MINPACK. VERSION OF NOVEMBER 1978.                                00000210
C     BURTON S. GARBOW, KENNETH E. HILLSTROM, JORGE J. MORE             00000220
C                                                                       00000230
C     **********                                                        00000240
      INTEGER I,IC,INFO,K,LWA,N,NFEV,NPROB,NREAD,NTRIES,NWRITE          00000250
      INTEGER NA(120),NF(120),NP(120),NX(120)                           00000260
      REAL FACTOR,F1,F2,GNORM1,GNORM2,ONE,TEN,TOL                       00000270
      REAL FVAL(120),GVEC(100),GNM(120),WA(6130),X(100)                 00000280
      REAL ENORM                                                        00000290
      EXTERNAL FCN                                                      00000300
      COMMON /REFNUM/ NPROB,NFEV                                        00000310
C                                                                       00000320
C     LOGICAL INPUT UNIT IS ASSUMED TO BE NUMBER 5.                     00000330
C     LOGICAL OUTPUT UNIT IS ASSUMED TO BE NUMBER 6.                    00000340
C                                                                       00000350
      DATA NREAD,NWRITE /5,6/                                           00000360
C                                                                       00000370
      DATA ONE,TEN,TOL /1.0E0,1.0E1,1.E-5/                              00000380
      LWA = 6130                                                        00000390
      IC = 0                                                            00000400
   10 CONTINUE                                                          00000410
         READ (NREAD,50) NPROB,N,NTRIES                                 00000420
         IF (NPROB .LE. 0) GO TO 30                                     00000430
         FACTOR = ONE                                                   00000440
         DO 20 K = 1, NTRIES                                            00000450
            IC = IC + 1                                                 00000460
            CALL INITPT(N,X,NPROB,FACTOR)                               00000470
            CALL OBJFCN(N,X,F1,NPROB)                                   00000480
            CALL GRDFCN(N,X,GVEC,NPROB)                                 00000490
            GNORM1 = ENORM(N,GVEC)                                      00000500
            WRITE (NWRITE,60) NPROB,N                                   00000510
            NFEV = 0                                                    00000520
            CALL SOLVER(FCN,N,X,F2,GVEC,TOL,INFO,WA,LWA)                00000530
            CALL OBJFCN(N,X,F2,NPROB)                                   00000540
            CALL GRDFCN(N,X,GVEC,NPROB)                                 00000550
            GNORM2 = ENORM(N,GVEC)                                      00000560
            NP(IC) = NPROB                                              00000570
            NA(IC) = N                                                  00000580
            NF(IC) = NFEV                                               00000590
            NX(IC) = INFO                                               00000600
            FVAL(IC) = F2                                               00000610
            GNM(IC) = GNORM2                                            00000620
            WRITE (NWRITE,70)                                           00000630
     *            F1,F2,GNORM1,GNORM2,NFEV,INFO,(X(I), I = 1, N)        00000640
            FACTOR = TEN*FACTOR                                         00000650
   20       CONTINUE                                                    00000660
         GO TO 10                                                       00000670
   30 CONTINUE                                                          00000680
      WRITE (NWRITE,80) IC                                              00000690
      WRITE (NWRITE,90)                                                 00000700
      DO 40 I = 1, IC                                                   00000710
         WRITE (NWRITE,100) NP(I),NA(I),NF(I),NX(I),FVAL(I),GNM(I)      00000720
   40    CONTINUE                                                       00000730
      STOP                                                              00000740
   50 FORMAT (3I5)                                                      00000750
   60 FORMAT ( //// 5X, 8H PROBLEM, I5, 5X, 10H DIMENSION, I5, 5X //)   00000760
   70 FORMAT (5X, 23H INITIAL FUNCTION VALUE, E15.7 // 5X,              00000770
     *        23H FINAL FUNCTION VALUE  , E15.7 // 5X,                  00000780
     *        23H INITIAL GRADIENT NORM , E15.7 // 5X,                  00000790
     *        23H FINAL GRADIENT NORM   , E15.7 // 5X,                  00000800
     *        33H NUMBER OF FUNCTION EVALUATIONS  , I10 // 5X,          00000810
     *        15H EXIT PARAMETER, 18X, I10 // 5X,                       00000820
     *        27H FINAL APPROXIMATE SOLUTION // (5X, 5E15.7))           00000830
   80 FORMAT (12H1SUMMARY OF , I3, 16H CALLS TO DRVCR1 /)               00000840
   90 FORMAT (25H NPROB   N    NFEV  INFO ,                             00000850
     *        42H FINAL FUNCTION VALUE  FINAL GRADIENT NORM /)          00000860
  100 FORMAT (I4, I6, I7, I6, 5X, E15.7, 6X, E15.7)                     00000870
C                                                                       00000880
C     LAST CARD OF DRIVER.                                              00000890
C                                                                       00000900
      END                                                               00000910
      SUBROUTINE FCN(N,X,F,GVEC,IFLAG)                                  00000920
      INTEGER N,IFLAG
      REAL F
      REAL X(N),GVEC(N)
C     **********
C
C     THE CALLING SEQUENCE OF FCN SHOULD BE IDENTICAL TO THE
C     CALLING SEQUENCE OF THE FUNCTION SUBROUTINE IN THE
C     UNCONSTRAINED OPTIMIZER. FCN SHOULD ONLY CALL THE TESTING
C     FUNCTION AND GRADIENT SUBROUTINES OBJFCN AND GRDFCN WITH
C     THE APPROPRIATE VALUE OF PROBLEM NUMBER (NPROB).
C
C     SUBPROGRAMS CALLED
C
C       MINPACK-SUPPLIED ... GRDFCN,OBJFCN
C
C     MINPACK. VERSION OF JULY 1978.
C     BURTON S. GARBOW, KENNETH E. HILLSTROM, JORGE J. MORE
C
C     **********
      INTEGER NPROB,NFEV
      COMMON /REFNUM/ NPROB,NFEV
      CALL OBJFCN(N,X,F,NPROB)
      CALL GRDFCN(N,X,GVEC,NPROB)
      NFEV = NFEV + 1
      RETURN
C
C     LAST CARD OF INTERFACE SUBROUTINE FCN.
C
      END
C ===== 14. DATA (NONLINEAR EQUATIONS).
    1    2    3                                                         00000010
    2    4    3                                                         00000020
    3    2    2                                                         00000030
    4    4    3                                                         00000040
    5    3    3                                                         00000050
    6    6    2                                                         00000060
    6    9    2                                                         00000070
    7    5    3                                                         00000080
    7    6    3                                                         00000090
    7    7    3                                                         00000100
    7    8    1                                                         00000110
    7    9    1                                                         00000120
    8   10    3                                                         00000130
    8   30    1                                                         00000140
    8   40    1                                                         00000150
    9   10    3                                                         00000160
   10    1    3                                                         00000170
   10   10    3                                                         00000180
   11   10    3                                                         00000190
   12   10    3                                                         00000200
   13   10    3                                                         00000210
   14   10    3                                                         00000220
    0    0    0                                                         00000230
C ===== 15. DATA (NONLINEAR LEAST SQUARES).
    1    5   10    1                                                    00000010
    1    5   50    1                                                    00000020
    2    5   10    1                                                    00000030
    2    5   50    1                                                    00000040
    3    5   10    1                                                    00000050
    3    5   50    1                                                    00000060
    4    2    2    3                                                    00000070
    5    3    3    3                                                    00000080
    6    4    4    3                                                    00000090
    7    2    2    3                                                    00000100
    8    3   15    3                                                    00000110
    9    4   11    3                                                    00000120
   10    3   16    3                                                    00000130
   11    6   31    3                                                    00000140
   11    9   31    3                                                    00000150
   11   12   31    3                                                    00000160
   12    3   10    1                                                    00000170
   13    2   10    1                                                    00000180
   14    4   20    3                                                    00000190
   15    1    8    3                                                    00000200
   15    8    8    1                                                    00000210
   15    9    9    1                                                    00000220
   15   10   10    1                                                    00000230
   16   10   10    3                                                    00000240
   16   30   30    1                                                    00000250
   16   40   40    1                                                    00000260
   17    5   33    1                                                    00000270
   18   11   65    1                                                    00000280
    0    0    0    0                                                    00000290
C ===== 16. DATA (UNCONSTRAINED NONLINEAR OPTIMIZATION).
    1    3    3                                                         00000010
    2    6    1                                                         00000020
    3    3    1                                                         00000030
    4    2    1                                                         00000040
    5    3    1                                                         00000050
    6   10    3                                                         00000060
    7    9    3                                                         00000070
    7   12    3                                                         00000080
    8   10    3                                                         00000090
    9    1    3                                                         00000100
    9    4    3                                                         00000110
    9   10    3                                                         00000120
   10    2    3                                                         00000130
   11    4    3                                                         00000140
   12    3    2                                                         00000150
   13   10    3                                                         00000160
   14    2    3                                                         00000170
   15    4    3                                                         00000180
   16    2    3                                                         00000190
   17    4    3                                                         00000200
   18    7    1                                                         00000210
   18    8    1                                                         00000220
   18    9    1                                                         00000230
   18   10    1                                                         00000240
    0    0    0                                                         00000250

