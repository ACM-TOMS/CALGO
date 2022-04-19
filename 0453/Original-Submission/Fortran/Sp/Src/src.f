      SUBROUTINE BROMIN ( N, S, TOL, XR, XI, WR, WI, EPS, IER )
      DOUBLE PRECISION AK, AN, ARG, CI, CR, D, D1, D2, E, EPS,
     &  FAC, FACTI, FACTR, PI, PR, QI, QR, RI, RR, S, T1, T2,
     &  TOL, U, V, WI, WR, XI, XR, YI, YR, Z
      DOUBLE PRECISION DGAMMA
      INTEGER IER, J, K, L, N, N1, NUM, NUP, IGNAL
      DIMENSION XR(N), XI(N), WR(N), WI(N)
C  THIS SUBROUTINE CALCULATES ABSCISSAS AND WEIGHTS OF THE
C  GAUSSIAN QUADRATURE FORMULA OF ORDER N FOR THE BROMWICH
C  INTEGRAL.  ONLY THE ABSCISSAS OF THE FIRST QUADRANT OF
C  THE COMPLEX PLANE, THE REAL ABSCISSA (IF N IS ODD) AND
C  THE CORRESPONDING WEIGHTS ARE CALCULATED.  THE OTHER
C  ABSCISSAS AND WEIGHTS ARE COMPLEX CONJUGATES.
C  INPUT PARAMETERS
C    N, ORDER OF THE QUADRATURE FORMULA.
C      N MUST BE GREATER THAN 2.
C    TOL, REQUESTED RELATIVE ACCURACY OF THE ABSCISSAS.
C    S, PARAMETERS OF THE WEIGHT FUNCTION.
C  OUTPUT PARAMETERS
C    XR AND XI CONTAIN THE REAL AND IMAGINARY PARTS OF
C      THE ABSCISSAS.  IF N IS ODD, THE REAL ABSCISSA
C      IS XR(1).
C    WR AND WI CONTAIN THE REAL AND IMAGINARY PARTS OF
C      THE CORRESPONDING WEIGHTS.
C    EPS IS A CRUDE ESTIMATION OF THE OBTAINED RELATIVE
C      ACCURACY OF THE ABSCISSAS.
C    IER IS AN ERROR CODE.
C      IF IER = 0, THE COMPUTATION IS CARRIED OUT TO
C        THE REQUESTED ACCURACY.
C      IF IER.GT.0, THE IER-TH ABSCISSA IS NOT FOUND.
C      IF IER = -1, THE COMPUTATIONS ARE CARRIED OUT,
C        BUT THE REQUESTED ACCURACY IS NOT
C        ACHIEVED.
C      IF IER = -2, N IS LESS THAN 3.
C  FUNCTION PROGRAMS REQUIRED
C    FUNCTION GAMMA(X), WHICH EVALUATES THE GAMMA
C      FUNCTION FOR POSITIVE X.
      IER = -2
      IF ( N .LT. 3 ) RETURN
      N1 = ( N + 1 ) / 2
      L = N - 1
      AN = N
      IER = 0
      EPS = TOL
      ARG = 0.034D0 * ( 30.D0 + AN + AN ) / ( AN - 1.D0 )
      FACTR = DCOS ( ARG )
      FACTI = DSIN ( ARG )
      FAC = 1.D0
      AK = 0.D0
      DO 10 K = 1, L
        AK = AK + 1.D0
        FAC = -FAC * AK
10    CONTINUE
      FAC = FAC * ( AN+AN+S-2.D0)**2 / ( AN * DGAMMA ( AN+S-1.D0 ) )
C  CALCULATION OF AN APPROXIMATION OF THE FIRST ABSCISSA.
      YR = 1.333D0 * AN + S - 1.5D0
      YI = 0.0D0
      IF ( N - N1 - N1 ) 30, 20, 20
20    YI = YI + 1.6D0 + 0.07D0 * S
C  START MAIN LOOP
30    DO 140 K = 1, N1
        E = TOL
        IGNAL = 0
        NUM = 0
        NUP = 0
C  NEWTON-RAPHSON METHOD.
        D = YR * YR + YI * YI
        YR = YR / D
        YI = -YI / D
        GO TO 50
40      IGNAL = 1
50      QR = S * YR - 1.D0
        QI = S * YI
        PR = (S+1.D0) * ( (S+2.D0) * ( YR*YR-YI*YI ) - 2.D0*YR ) + 1.D0
        PI = 2.D0 * ( S + 1.D0 ) * YI * ( ( S + 2.D0 ) * YR - 1.D0 )
        Z = 2.D0
        DO 60 J = 3, N
          RR = QR
          RI = QI
          QR = PR
          QI = PI
          Z = Z + 1.D0
          U = Z + S - 2.D0
          V = U + Z
          D = ( V * YR + ( 2.D0 - S ) / ( V - 2.D0 ) ) / U
          D1 = ( Z - 1.D0 ) * V / ( U * ( V - 2.D0 ) )
          D2 = V * YI / U
          PR = ( V - 1.D0 ) * ( QR * D - QI * D2 ) + D1 * RR
          PI = ( V - 1.D0 ) * ( QI * D + QR * D2 ) + D1 * RI
60      CONTINUE
        IF ( IGNAL .EQ. 1 ) GO TO 100
        D = ( YR * YR + YI * YI ) * V
        D1 = ( ( PR + QR ) * YR + ( PI + QI ) * YI ) / D + PR
        D2 = ( ( PI + QI ) * YR - ( PR + QR ) * YI ) / D + PI
        D = ( D1 * D1 + D2 * D2 ) * AN
        T1 = PR * YR - PI * YI
        T2 = PI * YR + PR * YI
        CR = ( T1 * D1 + T2 * D2 ) / D
        CI = ( T2 * D1 - T1 * D2 ) / D
        YR = YR - CR
        YI = YI - CI
        NUM = NUM + 1
C  TEST OF CONVERGENCE OF ITERATION PROCESS.
        IF ( CR*CR + CI*CI - E*E * ( YR*YR + YI*YI ) ) 40, 40, 70
C  TEST OF NUMBER OF ITERATION STEPS.
70      IF ( NUM - 10 ) 50, 50, 80
80      E = E * 10.D0
        IER = -1
        NUP = NUP + 1
        IF ( NUP - 5 ) 50, 50, 90
90      IER = K
        RETURN
C  CALCULATION OF WEIGHTS.
100     IF ( EPS .GE. E ) GO TO 110
        EPS = E
110     D = ( QR * QR + QI * QI )**2
        D1 = YR * QR + YI * QI
        D2 = YI * QR - YR * QI
        WR(K) = FAC * ( D1 * D1 - D2 * D2 ) / D
        WI(K) = 2.D0 * FAC * D2 * D1 / D
        D = YR * YR + YI * YI
        XR(K) = YR / D
        XI(K) = -YI / D
        IF ( K + 1 - N1 ) 130, 120, 150
120     FACTR = DCOS ( 1.5D0 * ARG )
        FACTI = DSIN ( 1.5D0 * ARG )
C  CALCULATION OF AN APPROXIMATION OF THE (K+1)-TH ABSCISSA.
130     YR = ( XR(K) + 0.67D0*AN ) * FACTR - XI(K) * FACTI - 0.67D0*AN
        YI = ( XR(K) + 0.67D0*AN ) * FACTI + XI(K) * FACTR
140   CONTINUE
150   RETURN
      END
