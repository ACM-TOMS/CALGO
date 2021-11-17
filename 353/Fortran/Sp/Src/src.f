      SUBROUTINE FSER1 ( F, EPS, MAX, M, C, S, LC, LS )

      IMPLICIT NONE

      REAL A
      REAL B
      REAL B1
      REAL B2
      REAL B3
      REAL B4
      REAL B5
      REAL C
      REAL C1
      REAL COSIP
      REAL EPS
      REAL F
      REAL F0
      REAL F1
      REAL FRAC
      REAL G
      REAL H
      INTEGER I
      INTEGER IN
      INTEGER LC
      INTEGER LLC
      INTEGER LLS
      INTEGER LN
      INTEGER LS
      INTEGER M
      INTEGER MAX
      INTEGER MSWTCH
      INTEGER N
      INTEGER NST
      INTEGER NSTOP
      REAL ODCOS
      REAL ODSIN
      REAL P
      REAL PI
      REAL PVT1
      REAL PVT2
      REAL S
      REAL S1
      REAL S1SQ
      REAL SUMCOS
      REAL SUMSIN
      REAL T
      REAL T1
      REAL T2
      REAL TEMP1
      REAL TEMP2
      REAL THA
      REAL TMAX
      REAL TP
      REAL TSQ
      REAL XI
      REAL XM

      PI = 3.1415926535898E+00
      XM = M
C
C  F1 = COS ( M * PI ) TEMPORARY.
C
      F1 = 1.0E+00 - 2.0E+00 * ( M - ( M / 2 ) * 2 )
      F0 = F ( 0.0E+00 )
      F1 = F ( 1.0E+00 ) * F1
C
C  'CIR' WILL BE USED THROUGHOUT THESE COMMENTS TO STAND FOR 'SIN' OR
C  'COS' WHENEVER THOSE TWO SYMBOLS MAY OCCUR.
C  NOW DEFINE SUMCIR OF THE ENDPOINTS.
C
      SUMCOS = ( F1 + F0 ) * 0.5E+00
      SUMSIN = 0.0E+00
      B1 = 2.0E+00 / 3.0E+00
C
C  TMAX IS THE SWITCH-OVER POINT IN THE ANGLE T.
C  OUR ANALYSIS INDICATES THAT TMAX = 1/6 IS THE BEST FOR THE ILIAC II,
C  WHICH HAS A 44 BIT FLOATING POINT MANTISSA.
C
      TMAX = 0.166E+00
C
C  N IS THE NUMBER OF ITERATIONS.  
C
      N = LOG ( 2.0E+00 * XM ) / 0.693E+00
C
C  BOTH TMAX AND N MAY BE CHANGED IF THE MACHINE FOR WHICH THIS
C  ROUTINE IS INTENDED HAS GREATER OR LESS ACCURACY THAN ILIAC II.
C  IF N IS CHANGED, THEN THE CORRESPONDING CHANGES MUST BE MADE
C  IN THE ASSIGNMENTS OF H AND NSTOP.
C
      H = 1.0E+00 / REAL( 2**N )
C
C  H = 1 / 2**N.
C
      NSTOP = 2**N - 1
C
C  NSTOP = 2**N - 1.
C
      T = H * XM
      TP = T * PI
      NST = 1
C ***** trh ******
C Initialized T2 to zero following Unassigned trap @ line259
      T2 = 0.0E0
C ***** end trh ******
      ASSIGN 67 TO MSWTCH
C
C  LLC AND LLS ARE USED BY THE ROUTINE IN COMPUTED-GO-TO STATEMENTS.
C  AS SOON AS LLS AND LLC HAVE BEEN DEFINED, WE CAN USE LS AND LC
C  AS RETURN PARAMETERS (SEE ABOVE).
C
      IF ( LS ) 1, 1, 2
1     LLS = 2
      GO TO 3
2     LLS = 1
      LS = MAX
3     IF ( LC ) 4, 4, 5
4     LLC = 2
      GO TO 7
5     LLC = 1
      LC = MAX
7     LN = 1
C
C  ALL OF THE ABOVE IS EXECUTED ONLY ONCE PER CALL.
C  NOW THE ITERATION BEGINS.
C
10    ODCOS = 0.0E+00
      ODSIN = 0.0E+00
C
C  BEGIN SUMMATION FOR ODCOS AND ODSIN.
C
      DO 65 I = 1, NSTOP, NST
        XI = I
        THA = XI * T
C
C  THA * PI IS THE ANGLE USED IN THIS I-TH TERM.
C
C  CIR(I*T*PI) IS CALCULATED HERE USING THE IDENTITY
C  CIR ( INTEGER MULTIPLE OF PI + FRACTIONAL MULTIPLE OF PI )
C  = COS ( INTEGER * PI ) * CIR ( FRAC * PI )
C  = ( + 0R - ) * CIR ( FRAC * PI ).
C
        FRAC = THA
        IN = THA
        THA = IN
        FRAC = ( FRAC - THA ) * PI
C
C  THA IS A FLOATING POINT INTEGER, FRAC IS THE FRACTIONAL PART * PI.
C
        COSIP = 1 - 2 * ( IN - 2 * ( IN / 2 ) )
        TEMP1 = COSIP * F ( XI * H )
C
C  TEMP1 = COS ( INTEGER PART ) * F ( I * H ).
C
        GO TO ( 50, 55 ) LLS
50      ODSIN = TEMP1 * SIN ( FRAC ) + ODSIN
55      GO TO ( 60, 65 ) LLC
60      ODCOS = TEMP1 * COS ( FRAC ) + ODCOS
65      CONTINUE

      GO TO MSWTCH, ( 67, 70 )
67    NST = 2
C
C  NOW HAVE MADE UP FOR THE FIRST 4 ITERATION STEPS, SO RESET THESE
C  TWO NUMBERS TO LOOK LIKE THE GENERAL CASE.
C
      NSTOP = 2**N
C
C  NSTOP = 2**N, IN CASE YOU CHANGE STARTING VALUE OF N.
C
      ASSIGN 70 TO MSWTCH
      GO TO 92
70    TSQ = TP * TP
      IF ( T - TMAX ) 74, 74, 75
C
C  74 IS THE POWER SERIES FOR SMALL T, 75 IS THE CLOSED FORM USED WITH
C  LARGER VALUES OF T.
C  THE POWER SERIES (WITH 'TN' = TP**N)
C  A = (2/45)*T3 - (2/315)*T5 + (2/4725)*T7
C  B = (2/3) + (2/15)*T2 - (4/105)*T4 + (2/567)*T6 - (4/22275)*T8
C  G = (4/3) - (2/15)*T2 - (1/210)*T4 + (1/11340)*T6
C  THE NEXT TERM IN G IS TOO SMALL.  IT IS (1/997920)*T8.
C
74    A = TP * TSQ * ( 1.0E+00 - TSQ * ( 1.0E+00 - TSQ / 15.0E+00 ) 
     &  / 7.0E+00 ) / 22.5E+00
      B2 = B1 * TSQ * 0.2E+00
      B3 = B2 * TSQ * 2.0E+00 / 7.0E+00
      B4 = B3 * TSQ / 10.8E+00
      B5 = B4 * TSQ * 14.0E+00 / 275.0E+00
      B = B1 + B2 - B3 + B4 - B5
      G = 2.0E+00 * B1 - B2 + B3 / 8.0E+00 - B4 / 40.0E+00
C
C  G = 2*B1 - B2 + B3/8 - B4/40 + 5*B5/896.  IF YOU WANT THE T8
C  TERM INCLUDED IN G.
C
      GO TO 80
C
C  CLOSED FORM OF THE COEFFICIENTS, WHERE AGAIN 'TN' MEANS TP**N.
C  A = 1/TP + COS(TP)*SIN(TP)/T2 - 2*(SIN(TP))**2/T3
C  B = 2*((1+(COS(TP))**2)/T2 - 2*SIN(TP)*COS(TP)/T3).
C  G = 4*( SIN(TP)/T3 - COS(TP)/T2 ).
C
75    IN = T
      TEMP1 = 1 - 2 * ( IN - 2 * ( IN / 2 ) )
      TEMP2 = IN
C
C  TEMP1 IS COS ( INTEGER PART OF TP), TEMP2 IS FRACTIONAL PART OF TP.
C
      TEMP2 = ( T - TEMP2 ) * PI
      S1 = TEMP1 * SIN ( TEMP2 )
C
C  S1 = SIN(TP).
C
      C1 = TEMP1 * COS ( TEMP2 )
C
C  C1 = COS(TP).
C
      P = S1 * C1
      S1SQ = S1 * S1
      A = ((( -2.0E+00 * S1SQ / TP ) + P ) / TP + 1.0E+00 ) / TP
      B = 2.0E+00 * ( ( -2.0E+00 * P / TP ) + 2.0E+00 - S1SQ ) / TSQ
      G = 4.0E+00 * ( S1 / TP - C1 ) / TSQ
80    GO TO ( 81, 85 ) LLS
C
C  HAVE CALCULATED THE COEFFICIENTS.  NOW READY FOR THE INTEGRATION
C  FORMULAS.
C
81    T2 = H * ( A * ( F0 - F1 ) + B * SUMSIN + G * ODSIN )
C
C  ENDT1 IS A SUBROUTINE WHICH CHECKS FOR THE CONVERGENCE OF THE
C  ITERATIONS.  ENDT1 REQUIRES THE PRESENT VALUE TO AGREE WITH THE
C  PREVIOUS VALUE TO WITHIN EPS2, WHERE
C  EPS2 = ( 1.0 + ABS ( PRESENT VALUE ) * EPS.
C  EPS IS SUPPLIED BY THE USER.
C
      CALL ENDT1 ( PVT2, T2, EPS, S, LLS, LN )
      GO TO ( 85, 84 ), LLS
84    LS = N
85    GO TO ( 86, 90 ), LLC
C
C  THIS IS THE COSINE INTEGRAL.
C
86     T1 = H * ( B * SUMCOS + G * ODCOS )
       CALL ENDT1 ( PVT1, T1, EPS, C, LLC, LN )
       GO TO ( 90, 89 ), LLC
89     LC = N
90     LN = 2
C
C  NOW TEST TO SEE IF DONE.
C
      IF ( LLC + LLS - 3 ) 92, 92, 100
92    N = N + 1
C
C  THIS IS THE BEGINNING OF THE ITERATION.
C
      IF ( N - MAX ) 95, 95, 100
95    H = 0.5E+00 * H
      T = 0.5E+00 * T
      TP = 0.5E+00 * TP
      NSTOP = 2 * NSTOP
      SUMSIN = SUMSIN + ODSIN
      SUMCOS = SUMCOS + ODCOS
      GO TO 10
100   S = T2
      C = T1
      RETURN
      END
      SUBROUTINE ENDT1 ( PREVQT, QUANT, EPS, VALUE, L1, L2 )

      IMPLICIT NONE

      REAL EPS
      INTEGER L1
      INTEGER L2
      REAL PREVQT
      REAL QUANT
      REAL REPS
      REAL VALUE

      GO TO ( 29, 20 ), L2
20    REPS = EPS * ( 1.0E+00 + ABS ( QUANT ) )
23    IF ( ABS ( PREVQT - QUANT ) - REPS ) 25, 25, 29
25    VALUE = QUANT
      L1 = 2
      GO TO 30
29    PREVQT = QUANT
30    RETURN
      END
