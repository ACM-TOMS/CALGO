C     ALGORITHM 600, COLLECTED ALGORITHMS FROM ACM.
C     ALGORITHM APPEARED IN ACM-TRANS. MATH. SOFTWARE, VOL.9, NO. 2,
C     JUN., 1983, P. 258-259.
      SUBROUTINE QUINAT(N, X, Y, B, C, D, E, F)                         QUI   10
C
      INTEGER N
      REAL X(N), Y(N), B(N), C(N), D(N), E(N), F(N)
C
C
C
C     QUINAT COMPUTES THE COEFFICIENTS OF A QUINTIC NATURAL QUINTIC SPLI
C     S(X) WITH KNOTS X(I) INTERPOLATING THERE TO GIVEN FUNCTION VALUES:
C               S(X(I)) = Y(I)  FOR I = 1,2, ..., N.
C     IN EACH INTERVAL (X(I),X(I+1)) THE SPLINE FUNCTION S(XX) IS A
C     POLYNOMIAL OF FIFTH DEGREE:
C     S(XX) = ((((F(I)*P+E(I))*P+D(I))*P+C(I))*P+B(I))*P+Y(I)    (*)
C           = ((((-F(I)*Q+E(I+1))*Q-D(I+1))*Q+C(I+1))*Q-B(I+1))*Q+Y(I+1)
C     WHERE  P = XX - X(I)  AND  Q = X(I+1) - XX.
C     (NOTE THE FIRST SUBSCRIPT IN THE SECOND EXPRESSION.)
C     THE DIFFERENT POLYNOMIALS ARE PIECED TOGETHER SO THAT S(X) AND
C     ITS DERIVATIVES UP TO S"" ARE CONTINUOUS.
C
C        INPUT:
C
C     N          NUMBER OF DATA POINTS, (AT LEAST THREE, I.E. N > 2)
C     X(1:N)     THE STRICTLY INCREASING OR DECREASING SEQUENCE OF
C                KNOTS.  THE SPACING MUST BE SUCH THAT THE FIFTH POWER
C                OF X(I+1) - X(I) CAN BE FORMED WITHOUT OVERFLOW OR
C                UNDERFLOW OF EXPONENTS.
C     Y(1:N)     THE PRESCRIBED FUNCTION VALUES AT THE KNOTS.
C
C        OUTPUT:
C
C     B,C,D,E,F  THE COMPUTED SPLINE COEFFICIENTS AS IN (*).
C         (1:N)  SPECIFICALLY
C                B(I) = S'(X(I)), C(I) = S"(X(I))/2, D(I) = S"'(X(I))/6,
C                E(I) = S""(X(I))/24,  F(I) = S""'(X(I))/120.
C                F(N) IS NEITHER USED NOR ALTERED.  THE FIVE ARRAYS
C                B,C,D,E,F MUST ALWAYS BE DISTINCT.
C
C        OPTION:
C
C     IT IS POSSIBLE TO SPECIFY VALUES FOR THE FIRST AND SECOND
C     DERIVATIVES OF THE SPLINE FUNCTION AT ARBITRARILY MANY KNOTS.
C     THIS IS DONE BY RELAXING THE REQUIREMENT THAT THE SEQUENCE OF
C     KNOTS BE STRICTLY INCREASING OR DECREASING.  SPECIFICALLY:
C
C     IF X(J) = X(J+1) THEN S(X(J)) = Y(J) AND S'(X(J)) = Y(J+1),
C     IF X(J) = X(J+1) = X(J+2) THEN IN ADDITION S"(X(J)) = Y(J+2).
C
C     NOTE THAT S""(X) IS DISCONTINUOUS AT A DOUBLE KNOT AND, IN
C     ADDITION, S"'(X) IS DISCONTINUOUS AT A TRIPLE KNOT.  THE
C     SUBROUTINE ASSIGNS Y(I) TO Y(I+1) IN THESE CASES AND ALSO TO
C     Y(I+2) AT A TRIPLE KNOT.  THE REPRESENTATION (*) REMAINS
C     VALID IN EACH OPEN INTERVAL (X(I),X(I+1)).  AT A DOUBLE KNOT,
C     X(J) = X(J+1), THE OUTPUT COEFFICIENTS HAVE THE FOLLOWING VALUES:
C       Y(J) = S(X(J))          = Y(J+1)
C       B(J) = S'(X(J))         = B(J+1)
C       C(J) = S"(X(J))/2       = C(J+1)
C       D(J) = S"'(X(J))/6      = D(J+1)
C       E(J) = S""(X(J)-0)/24     E(J+1) = S""(X(J)+0)/24
C       F(J) = S""'(X(J)-0)/120   F(J+1) = S""'(X(J)+0)/120
C     AT A TRIPLE KNOT, X(J) = X(J+1) = X(J+2), THE OUTPUT
C     COEFFICIENTS HAVE THE FOLLOWING VALUES:
C       Y(J) = S(X(J))         = Y(J+1)    = Y(J+2)
C       B(J) = S'(X(J))        = B(J+1)    = B(J+2)
C       C(J) = S"(X(J))/2      = C(J+1)    = C(J+2)
C       D(J) = S"'((X(J)-0)/6    D(J+1) = 0  D(J+2) = S"'(X(J)+0)/6
C       E(J) = S""(X(J)-0)/24    E(J+1) = 0  E(J+2) = S""(X(J)+0)/24
C       F(J) = S""'(X(J)-0)/120  F(J+1) = 0  F(J+2) = S""'(X(J)+0)/120
C
      INTEGER I, M
      REAL B1, P, PQ, PQQR, PR, P2, P3, Q, QR, Q2, Q3, R, R2, S, T, U, V
C
      IF (N.LE.2) GO TO 190
C
C     COEFFICIENTS OF A POSITIVE DEFINITE, PENTADIAGONAL MATRIX,
C     STORED IN D,E,F FROM 2 TO N-2.
C
      M = N - 2
      Q = X(2) - X(1)
      R = X(3) - X(2)
      Q2 = Q*Q
      R2 = R*R
      QR = Q + R
      D(1) = 0.
      E(1) = 0.
      D(2) = 0.
      IF (Q.NE.0.) D(2) = 6.*Q*Q2/(QR*QR)
C
      IF (M.LT.2) GO TO 40
      DO 30 I=2,M
        P = Q
        Q = R
        R = X(I+2) - X(I+1)
        P2 = Q2
        Q2 = R2
        R2 = R*R
        PQ = QR
        QR = Q + R
        IF (Q) 20, 10, 20
   10   D(I+1) = 0.
        E(I) = 0.
        F(I-1) = 0.
        GO TO 30
   20   Q3 = Q2*Q
        PR = P*R
        PQQR = PQ*QR
        D(I+1) = 6.*Q3/(QR*QR)
        D(I) = D(I) + (Q+Q)*(15.*PR*PR+(P+R)*Q*(20.*PR+7.*Q2)+Q2*(8.*
     *   (P2+R2)+21.*PR+Q2+Q2))/(PQQR*PQQR)
        D(I-1) = D(I-1) + 6.*Q3/(PQ*PQ)
        E(I) = Q2*(P*QR+3.*PQ*(QR+R+R))/(PQQR*QR)
        E(I-1) = E(I-1) + Q2*(R*PQ+3.*QR*(PQ+P+P))/(PQQR*PQ)
        F(I-1) = Q3/PQQR
   30 CONTINUE
C
   40 IF (R.NE.0.) D(M) = D(M) + 6.*R*R2/(QR*QR)
C
C     FIRST AND SECOND ORDER DIVIDED DIFFERENCES OF THE GIVEN FUNCTION
C     VALUES, STORED IN B FROM 2 TO N AND IN C FROM 3 TO N
C     RESPECTIVELY. CARE IS TAKEN OF DOUBLE AND TRIPLE KNOTS.
C
      DO 60 I=2,N
        IF (X(I).NE.X(I-1)) GO TO 50
        B(I) = Y(I)
        Y(I) = Y(I-1)
        GO TO 60
   50   B(I) = (Y(I)-Y(I-1))/(X(I)-X(I-1))
   60 CONTINUE
      DO 80 I=3,N
        IF (X(I).NE.X(I-2)) GO TO 70
        C(I) = B(I)*0.5
        B(I) = B(I-1)
        GO TO 80
   70   C(I) = (B(I)-B(I-1))/(X(I)-X(I-2))
   80 CONTINUE
C
C     SOLVE THE LINEAR SYSTEM WITH C(I+2) - C(I+1) AS RIGHT-HAND SIDE.
C
      IF (M.LT.2) GO TO 100
      P = 0.
      C(1) = 0.
      E(M) = 0.
      F(1) = 0.
      F(M-1) = 0.
      F(M) = 0.
      C(2) = C(4) - C(3)
      D(2) = 1./D(2)
C
      IF (M.LT.3) GO TO 100
      DO 90 I=3,M
        Q = D(I-1)*E(I-1)
        D(I) = 1./(D(I)-P*F(I-2)-Q*E(I-1))
        E(I) = E(I) - Q*F(I-1)
        C(I) = C(I+2) - C(I+1) - P*C(I-2) - Q*C(I-1)
        P = D(I-1)*F(I-1)
   90 CONTINUE
C
  100 I = N - 1
      C(N-1) = 0.
      C(N) = 0.
      IF (N.LT.4) GO TO 120
      DO 110 M=4,N
C        I = N-2, ..., 2
        I = I - 1
        C(I) = (C(I)-E(I)*C(I+1)-F(I)*C(I+2))*D(I)
  110 CONTINUE
C
C     INTEGRATE THE THIRD DERIVATIVE OF S(X).
C
  120 M = N - 1
      Q = X(2) - X(1)
      R = X(3) - X(2)
      B1 = B(2)
      Q3 = Q*Q*Q
      QR = Q + R
      IF (QR) 140, 130, 140
  130 V = 0.
      T = 0.
      GO TO 150
  140 V = C(2)/QR
      T = V
  150 F(1) = 0.
      IF (Q.NE.0.) F(1) = V/Q
      DO 180 I=2,M
        P = Q
        Q = R
        R = 0.
        IF (I.NE.M) R = X(I+2) - X(I+1)
        P3 = Q3
        Q3 = Q*Q*Q
        PQ = QR
        QR = Q + R
        S = T
        T = 0.
        IF (QR.NE.0.) T = (C(I+1)-C(I))/QR
        U = V
        V = T - S
        IF (PQ) 170, 160, 170
  160   C(I) = C(I-1)
        D(I) = 0.
        E(I) = 0.
        F(I) = 0.
        GO TO 180
  170   F(I) = F(I-1)
        IF (Q.NE.0.) F(I) = V/Q
        E(I) = 5.*S
        D(I) = 10.*(C(I)-Q*S)
        C(I) = D(I)*(P-Q) + (B(I+1)-B(I)+(U-E(I))*P3-(V+E(I))*Q3)/PQ
        B(I) = (P*(B(I+1)-V*Q3)+Q*(B(I)-U*P3))/PQ -
     *   P*Q*(D(I)+E(I)*(Q-P))
  180 CONTINUE
C
C     END POINTS X(1) AND X(N).
C
      P = X(2) - X(1)
      S = F(1)*P*P*P
      E(1) = 0.
      D(1) = 0.
      C(1) = C(2) - 10.*S
      B(1) = B1 - (C(1)+S)*P
C
      Q = X(N) - X(N-1)
      T = F(N-1)*Q*Q*Q
      E(N) = 0.
      D(N) = 0.
      C(N) = C(N-1) + 10.*T
      B(N) = B(N) + (C(N)-T)*Q
  190 RETURN
      END
C                                                                       QUI   10
C                                                                       QUI   20
      SUBROUTINE QUINEQ(N, Y, B, C, D, E, F)                            QUI   30
C
      INTEGER N
      REAL Y(N), B(N), C(N), D(N), E(N), F(N)
C
C
C
C     QUINEQ COMPUTES THE COEFFICIENTS OF QUINTIC NATURAL QUINTIC SPLINE
C     S(X) WITH EQUIDISTANT KNOTS X(I) INTERPOLATING THERE TO GIVEN
C     FUNCTION VALUES:
C               S(X(I)) = Y(I)  FOR I = 1,2, ..., N.
C     IN EACH INTERVAL (X(I),X(I+1)) THE SPLINE FUNCTION S(XX) IS
C     A POLYNOMIAL OF FIFTH DEGREE:
C     S(XX)=((((F(I)*P+E(I))*P+D(I))*P+C(I))*P+B(I))*P+Y(I)    (*)
C          =((((-F(I)*Q+E(I+1))*Q-D(I+1))*Q+C(I+1))*Q-B(I+1))*Q+Y(I+1)
C     WHERE  P = (XX - X(I))/X(I+1) - X(I))
C     AND    Q = (X(I+1) - XX)/(X(I+1) - X(I)).
C     (NOTE THE FIRST SUBSCRIPT IN THE SECOND EXPRESSION.)
C     THE DIFFERENT POLYNOMIALS ARE PIECED TOGETHER SO THAT S(X) AND
C     ITS DERIVATIVES UP TO S"" ARE CONTINUOUS.
C
C        INPUT:
C
C     N          NUMBER OF DATA POINTS, (AT LEAST THREE, I.E. N > 2)
C     Y(1:N)     THE PRESCRIBED FUNCTION VALUES AT THE KNOTS
C
C        OUTPUT:
C
C     B,C,D,E,F  THE COMPUTED SPLINE COEFFICIENTS AS IN (*).
C         (1:N)  IF X(I+1) - X(I) = 1., THEN SPECIFICALLY:
C                B(I) = S'X(I)), C(I) = S"(X(I))/2, D(I) = S"'(X(I))/6,
C                E(I) = S""(X(I))/24,  F(I) = S""'(X(I)+0)/120.
C                F(N) IS NEITHER USED NOR ALTERED.  THE ARRAYS
C                Y,B,C,D MUST ALWAYS BE DISTINCT.  IF E AND F ARE
C                NOT WANTED, THE CALL QUINEQ(N,Y,B,C,D,D,D) MAY
C                BE USED TO SAVE STORAGE LOCATIONS.
C
      INTEGER I, M
      REAL P, Q, R, S, T, U, V
C
      IF (N.LE.2) GO TO 50
C
      M = N - 3
      P = 0.
      Q = 0.
      R = 0.
      S = 0.
      T = 0.
      D(M+1) = 0.
      D(M+2) = 0.
      IF (M.LE.0) GO TO 30
      DO 10 I=1,M
        U = P*R
        B(I) = 1./(66.-U*R-Q)
        R = 26. - U
        C(I) = R
        D(I) = Y(I+3) - 3.*(Y(I+2)-Y(I+1)) - Y(I) - U*S - Q*T
        Q = P
        P = B(I)
        T = S
        S = D(I)
   10 CONTINUE
C
      I = N - 2
      DO 20 M=4,N
C        I    = N-3, ..., 1
        I = I - 1
        D(I) = (D(I)-C(I)*D(I+1)-D(I+2))*B(I)
   20 CONTINUE
C
   30 M = N - 1
      Q = 0.
      R = D(1)
      T = R
      V = R
      DO 40 I=2,M
        P = Q
        Q = R
        R = D(I)
        S = T
        T = P - Q - Q + R
        F(I) = T
        U = 5.*(-P+Q)
        E(I) = U
        D(I) = 10.*(P+Q)
        C(I) = 0.5*(Y(I+1)+Y(I-1)+S-T) - Y(I) - U
        B(I) = 0.5*(Y(I+1)-Y(I-1)-S-T) - D(I)
   40 CONTINUE
C
      F(1) = V
      E(1) = 0.
      D(1) = 0.
      C(1) = C(2) - 10.*V
      B(1) = Y(2) - Y(1) - C(1) - V
      E(N) = 0.
      D(N) = 0.
      C(N) = C(N-1) + 10.*T
      B(N) = Y(N) - Y(N-1) + C(N) - T
   50 RETURN
      END
C                                                                       QUI   10
C                                                                       QUI   20
C                                                                       QUI   30
C                                                                       QUI   40
      SUBROUTINE QUINDF(N, X, Y, B, C, D, E, F)                         QUI   50
C
      INTEGER N
      REAL X(N), Y(N), B(N), C(N), D(N), E(N), F(N)
C
C
C
C     QUINDF COMPUTES THE COEFFICIENTS OF A QUINTIC NATURAL QUINTIC
C     SPLINE S(X) WITH KNOTS X(I) FOR WHICH THE FUNCTION VALUES Y(I) AND
C     THE FIRST DERIVATIVES B(I) ARE SPECIFIED AT X(I), I = 1,2, ...,N.
C     IN EACH INTERVAL (X(I),X(I+1)) THE SPLINE FUNCTION S(XX) IS
C     A POLYNOMIAL OF FIFTH DEGREE:
C     S(XX)=((((F(I)*P+E(I))*P+D(I))*P+C(I))*P+B(I))*P+Y(I)    (*)
C     WHERE  P = XX - X(I).
C
C        INPUT:
C
C     N          NUMBER OF DATA POINTS, (AT LEAST TWO, I.E. N > 1)
C     X(1:N)     THE STRICTLY INCREASING OR DECREASING SEQUENCE OF
C                KNOTS.  THE SPACING MUST BE SUCH THAT THE FIFTH
C                POWER OF X(I+1) - X(I) CAN BE FORMED WITHOUT
C                OVERFLOW OR UNDERFLOW OF EXPONENTS
C     Y(1:N)     THE PRESCRIBED FUNCTION VALUES AT THE KNOTS
C     B(1:N)     THE PRESCRIBED DERIVATIVE VALUES AT THE KNOTS
C
C        OUTPUT:
C
C     C,D,E,F    THE COMPUTED SPLINE COEFFICIENTS AS IN (*).
C         (1:N)  E(N) AND F(N) ARE NEITHER USED NOR ALTERED.
C                THE ARRAYS C,D,E,F MUST ALWAYS BE DISTINCT.
C
      INTEGER I, M, N1
      REAL CC, G, H, HH, H2, P, PP, Q, QQ, R, RR
C
      IF (N.LE.1) GO TO 40
      N1 = N - 1
      CC = 0.
      HH = 0.
      PP = 0.
      QQ = 0.
      RR = 0.
      G = 0.
      DO 10 I=1,N1
        H = 1./(X(I+1)-X(I))
        H2 = H*H
        D(I) = 3.*(HH+H) - G*HH
        P = (Y(I+1)-Y(I))*H2*H
        Q = (B(I+1)+B(I))*H2
        R = (B(I+1)-B(I))*H2
        CC = 10.*(P-PP) - 5.*(Q-QQ) + R + RR + G*CC
        C(I) = CC
        G = H/D(I)
        HH = H
        PP = P
        QQ = Q
        RR = R
   10 CONTINUE
C
      C(N) = (-10.*PP+5.*QQ+RR+G*CC)/(3.*HH-G*HH)
      I = N
      DO 20 M=1,N1
C        I      = N-1, ..., 1
        I = I - 1
        D(I+1) = 1./(X(I+1)-X(I))
        C(I) = (C(I)+C(I+1)*D(I+1))/D(I)
   20 CONTINUE
C
      DO 30 I=1,N1
        H = D(I+1)
        P = (((Y(I+1)-Y(I))*H-B(I))*H-C(I))*H
        Q = ((B(I+1)-B(I))*H-C(I)-C(I))*H
        R = (C(I+1)-C(I))*H
        G = Q - 3.*P
        RR = R - 3.*(P+G)
        QQ = -RR - RR + G
        F(I) = RR*H*H
        E(I) = QQ*H
        D(I) = -RR - QQ + P
   30 CONTINUE
C
      D(N) = 0.
      E(N) = 0.
      F(N) = 0.
   40 RETURN
      END
C    DRIVER PROGRAM FOR TEST OF QUINAT, QUINEQ AND QUINDF               MAN   10
C       FOLLOWS.                                                        MAN   20
C                                                                       MAN   30
      INTEGER N,NM1,M,MM,MM1,I,K,J,JJ                                   MAN   40
      REAL Z                                                            MAN   50
      REAL X(200),Y(200),B(200),BB(200),CC(200),DD(200),EE(200),        MAN   60
     *                 FF(200),A(200,6),C(6),DIFF(5),COM(5)             MAN   70
C                                                                       MAN   80
C     N          NUMBER OF DATA POINTS.                                 MAN   90
C     M          2*M-1 IS ORDER OF SPLINE.                              MAN  100
C                   M = 3 ALWAYS FOR QUINTIC SPLINE.                    MAN  110
C     NN,NM1,MM,                                                        MAN  120
C     MM1,I,K,                                                          MAN  130
C     J,JJ       TEMPORARY INTEGER VARIABLES.                           MAN  140
C     Z,P        TEMPORARY REAL VARIABLES.                              MAN  150
C     X(1:N)     THE SEQUENCE OF KNOTS.                                 MAN  160
C     Y(1:N)     THE PRESCRIBED FUNCTION VALUES AT THE KNOTS.           MAN  170
C     B(1:N)     THE PRESCRIBED DERIVATIVE VALUES AT THE KNOTS.         MAN  180
C     BB,CC,DD,                                                         MAN  190
C     EE,FF(1:N) THE COMPUTED SPLINE COEFFICIENTS                       MAN  200
C     A(1:N,1:6) TWO DIMENSIONAL ARRAY WHOSE COLUMNS ARE                MAN  210
C                   Y, BB, CC, DD, EE, FF.                              MAN  220
C     DIFF(1:5)  MAXIMUM VALUES OF DIFFERENCES OF VALUES AND            MAN  230
C                   DERIVATIVES TO RIGHT AND LEFT OF KNOTS.             MAN  240
C     COM(1:5)   MAXIMUM VALUES OF COEFFICIENTS.                        MAN  250
C                                                                       MAN  260
C                                                                       MAN  270
C     TEST OF QUINAT WITH NONEQUIDISTANT KNOTS AND                      MAN  280
C        EQUIDISTANT KNOTS FOLLOWS.                                     MAN  290
C                                                                       MAN  300
      NOUT=6                                                            MAN  310
      WRITE(NOUT,10)                                                    MAN  320
   10 FORMAT(50H1         TEST OF QUINAT WITH NONEQUIDISTANT KNOTS)     MAN  330
      N = 5                                                             MAN  340
      X(1) = -3.0                                                       MAN  350
      X(2) = -1.0                                                       MAN  360
      X(3) =  0.0                                                       MAN  370
      X(4) =  3.0                                                       MAN  380
      X(5) =  4.0                                                       MAN  390
      Y(1) =  7.0                                                       MAN  400
      Y(2) = 11.0                                                       MAN  410
      Y(3) = 26.0                                                       MAN  420
      Y(4) = 56.0                                                       MAN  430
      Y(5) = 29.0                                                       MAN  440
      M = 3                                                             MAN  450
      MM = 2*M                                                          MAN  460
      MM1 = MM - 1                                                      MAN  470
      WRITE(NOUT,15) N,M                                                MAN  480
   15 FORMAT(5H-N = ,I3,7H    M =,I2)                                   MAN  490
      CALL QUINAT(N,X,Y,BB,CC,DD,EE,FF)                                 MAN  500
      DO 20 I = 1,N                                                     MAN  510
         A(I,1) = Y(I)                                                  MAN  520
         A(I,2) = BB(I)                                                 MAN  530
         A(I,3) = CC(I)                                                 MAN  540
         A(I,4) = DD(I)                                                 MAN  550
         A(I,5) = EE(I)                                                 MAN  560
         A(I,6) = FF(I)                                                 MAN  570
   20 CONTINUE                                                          MAN  580
      DO 30 I = 1,MM1                                                   MAN  590
         DIFF(I) = 0.0                                                  MAN  600
         COM(I) = 0.0                                                   MAN  610
   30 CONTINUE                                                          MAN  620
      DO 70 K = 1,N                                                     MAN  630
         DO 35 I = 1,MM                                                 MAN  640
            C(I) = A(K,I)                                               MAN  650
   35 CONTINUE                                                          MAN  660
      WRITE(NOUT,40) K                                                  MAN  670
   40 FORMAT(40H ---------------------------------------,I3,            MAN  680
     *  45H --------------------------------------------)               MAN  690
      WRITE(NOUT,45) X(K)                                               MAN  700
   45 FORMAT(F12.8)                                                     MAN  710
      IF (K .EQ. N) WRITE(NOUT,55) C(1)                                 MAN  720
      IF (K .EQ. N) GO TO 75                                            MAN  730
      WRITE(NOUT,55) (C(I),I=1,MM)                                      MAN  740
      DO 50 I = 1,MM1                                                   MAN  750
         IF (ABS(A(K,I)) .GT. COM(I)) COM(I) = ABS(A(K,I))              MAN  760
   50 CONTINUE                                                          MAN  770
   55 FORMAT(6F16.8)                                                    MAN  780
      Z = X(K+1) - X(K)                                                 MAN  790
      DO 60 I = 2,MM                                                    MAN  800
         DO 60 JJ = I,MM                                                MAN  810
            J = MM + I - JJ                                             MAN  820
            C(J-1) = C(J)*Z + C(J-1)                                    MAN  830
   60 CONTINUE                                                          MAN  840
      WRITE(NOUT,55) (C(I),I=1,MM)                                      MAN  850
      DO 65 I = 1,MM1                                                   MAN  860
         IF ((K .GE. N-1) .AND. (I .NE. 1)) GO TO 65                    MAN  870
         Z = ABS(C(I) - A(K+1,I))                                       MAN  880
         IF (Z .GT. DIFF(I)) DIFF(I) = Z                                MAN  890
   65 CONTINUE                                                          MAN  900
   70 CONTINUE                                                          MAN  910
   75 WRITE(NOUT,80)                                                    MAN  920
   80 FORMAT(41H  MAXIMUM ABSOLUTE VALUES OF DIFFERENCES )              MAN  930
      WRITE(NOUT,85) (DIFF(I),I=1,MM1)                                  MAN  940
   85 FORMAT(5E18.9)                                                    MAN  950
      WRITE(NOUT,90)                                                    MAN  960
   90 FORMAT(42H  MAXIMUM ABSOLUTE VALUES OF COEFFICIENTS )             MAN  970
      IF (ABS(C(1)) .GT. COM(1)) COM(1) =ABS(C(1))                      MAN  980
      WRITE(NOUT,95) (COM(I),I=1,MM1)                                   MAN  990
   95 FORMAT(5F16.8)                                                    MAN 1000
      M = 3                                                             MAN 1010
      DO 200 N = 10,100,10                                              MAN 1020
      MM = 2*M                                                          MAN 1030
      MM1 = MM - 1                                                      MAN 1040
      NM1 = N -1                                                        MAN 1050
      DO 100 I = 1,NM1,2                                                MAN 1060
         X(I)   = I                                                     MAN 1070
         X(I+1) = I + 1                                                 MAN 1080
         Y(I)   = 1.                                                    MAN 1090
         Y(I+1) = 0.                                                    MAN 1100
  100 CONTINUE                                                          MAN 1110
      IF (MOD(N,2) .EQ. 0) GOTO 105                                     MAN 1120
      X(N) = N                                                          MAN 1130
      Y(N) = 1.                                                         MAN 1140
  105 CONTINUE                                                          MAN 1150
      WRITE(NOUT,110) N,M                                               MAN 1160
  110 FORMAT(5H-N = ,I3,7H    M =,I2)                                   MAN 1170
      CALL QUINAT(N,X,Y,BB,CC,DD,EE,FF)                                 MAN 1180
      DO 115 I = 1,N                                                    MAN 1190
         A(I,1) = Y(I)                                                  MAN 1200
         A(I,2) = BB(I)                                                 MAN 1210
         A(I,3) = CC(I)                                                 MAN 1220
         A(I,4) = DD(I)                                                 MAN 1230
         A(I,5) = EE(I)                                                 MAN 1240
         A(I,6) = FF(I)                                                 MAN 1250
  115 CONTINUE                                                          MAN 1260
      DO 120 I = 1, MM1                                                 MAN 1270
         DIFF(I) = 0.0                                                  MAN 1280
         COM(I) = 0.0                                                   MAN 1290
  120 CONTINUE                                                          MAN 1300
      DO 165 K = 1,N                                                    MAN 1310
         DO 125 I = 1,MM                                                MAN 1320
            C(I) = A(K,I)                                               MAN 1330
  125    CONTINUE                                                       MAN 1340
         IF (N .GT. 10) GOTO 140                                        MAN 1350
         WRITE(NOUT,130) K                                              MAN 1360
  130 FORMAT(40H ---------------------------------------,I3,            MAN 1370
     *       45H --------------------------------------------)          MAN 1380
         WRITE(NOUT,135) X(K)                                           MAN 1390
  135 FORMAT(F12.8)                                                     MAN 1400
         IF (K .EQ. N) WRITE(NOUT,150) C(1)                             MAN 1410
  140    CONTINUE                                                       MAN 1420
         IF (K .EQ. N) GO TO 170                                        MAN 1430
         IF (N .LE. 10) WRITE(NOUT,150) (C(I), I=1,MM)                  MAN 1440
         DO 145 I = 1,MM1                                               MAN 1450
            IF (ABS(A(K,I)) .GT.  COM(I)) COM(I) = ABS(A(K,I))          MAN 1460
  145    CONTINUE                                                       MAN 1470
  150 FORMAT(6F16.8)                                                    MAN 1480
         Z = X(K+1) - X(K)                                              MAN 1490
         DO 155 I = 2,MM                                                MAN 1500
            DO 155 JJ = I,MM                                            MAN 1510
               J = MM + I - JJ                                          MAN 1520
               C(J-1) = C(J)*Z + C(J-1)                                 MAN 1530
  155    CONTINUE                                                       MAN 1540
         IF (N .LE. 10) WRITE(NOUT,150) (C(I), I=1,MM)                  MAN 1550
         DO 160 I = 1,MM1                                               MAN 1560
            IF ((K .GE. N-1) .AND. (I .NE. 1)) GO TO 160                MAN 1570
            Z = ABS(C(I) - A(K+1,I))                                    MAN 1580
            IF (Z .GT. DIFF(I)) DIFF(I) = Z                             MAN 1590
  160    CONTINUE                                                       MAN 1600
  165 CONTINUE                                                          MAN 1610
  170 WRITE(NOUT,175)                                                   MAN 1620
  175 FORMAT(41H  MAXIMUM ABSOLUTE VALUES OF DIFFERENCES )              MAN 1630
      WRITE(NOUT,180) (DIFF(I),I=1,MM1)                                 MAN 1640
  180 FORMAT(5E18.9)                                                    MAN 1650
      WRITE(NOUT,185)                                                   MAN 1660
  185 FORMAT(42H  MAXIMUM ABSOLUTE VALUES OF COEFFICIENTS )             MAN 1670
      IF (ABS(C(1)) .GT. COM(1)) COM(1) = ABS(C(1))                     MAN 1680
      WRITE(NOUT,190) (COM(I),I=1,MM1)                                  MAN 1690
  190 FORMAT(5E18.9)                                                    MAN 1700
  200 CONTINUE                                                          MAN 1710
C                                                                       MAN 1720
C                                                                       MAN 1730
C     TEST OF QUINEQ FOLLOWS.                                           MAN 1740
C                                                                       MAN 1750
      WRITE(NOUT,210)                                                   MAN 1760
  210 FORMAT(18H1   TEST OF QUINEQ)                                     MAN 1770
      M = 3                                                             MAN 1780
      DO 400 N = 10,100,10                                              MAN 1790
      MM = 2*M                                                          MAN 1800
      MM1 = MM - 1                                                      MAN 1810
      NM1 = N -1                                                        MAN 1820
      DO 300 I = 1,NM1,2                                                MAN 1830
         X(I)   = I                                                     MAN 1840
         X(I+1) = I + 1                                                 MAN 1850
         Y(I)   = 1.                                                    MAN 1860
         Y(I+1) = 0.                                                    MAN 1870
  300 CONTINUE                                                          MAN 1880
      IF (MOD(N,2) .EQ. 0) GOTO 305                                     MAN 1890
      X(N) = FLOAT(N)                                                   MAN 1900
      Y(N) = 1.                                                         MAN 1910
  305 CONTINUE                                                          MAN 1920
      WRITE(NOUT,310) N,M                                               MAN 1930
  310 FORMAT(5H-N = ,I3,7H    M =,I2)                                   MAN 1940
      CALL QUINEQ(N,Y,BB,CC,DD,EE,FF)                                   MAN 1950
      DO 315 I = 1,N                                                    MAN 1960
         A(I,1) = Y(I)                                                  MAN 1970
         A(I,2) = BB(I)                                                 MAN 1980
         A(I,3) = CC(I)                                                 MAN 1990
         A(I,4) = DD(I)                                                 MAN 2000
         A(I,5) = EE(I)                                                 MAN 2010
         A(I,6) = FF(I)                                                 MAN 2020
  315 CONTINUE                                                          MAN 2030
      DO 320 I = 1, MM1                                                 MAN 2040
         DIFF(I) = 0.0                                                  MAN 2050
         COM(I) = 0.0                                                   MAN 2060
  320 CONTINUE                                                          MAN 2070
      DO 365 K = 1,N                                                    MAN 2080
         DO 325 I = 1,MM                                                MAN 2090
            C(I) = A(K,I)                                               MAN 2100
  325    CONTINUE                                                       MAN 2110
         IF (N .GT. 10) GOTO 340                                        MAN 2120
         WRITE(NOUT,330) K                                              MAN 2130
  330 FORMAT(40H ---------------------------------------,I3,            MAN 2140
     *       45H --------------------------------------------)          MAN 2150
         WRITE(NOUT,335) X(K)                                           MAN 2160
  335 FORMAT(F12.8)                                                     MAN 2170
         IF (K .EQ. N) WRITE(NOUT,350) C(1)                             MAN 2180
  340    CONTINUE                                                       MAN 2190
         IF (K .EQ. N) GO TO 370                                        MAN 2200
         IF (N .LE. 10) WRITE(NOUT,350) (C(I), I=1,MM)                  MAN 2210
         DO 345 I = 1,MM1                                               MAN 2220
            IF (ABS(A(K,I)) .GT.  COM(I)) COM(I) = ABS(A(K,I))          MAN 2230
  345    CONTINUE                                                       MAN 2240
  350 FORMAT(6F16.8)                                                    MAN 2250
         Z = 1.                                                         MAN 2260
         DO 355 I = 2,MM                                                MAN 2270
            DO 355 JJ = I,MM                                            MAN 2280
               J = MM + I - JJ                                          MAN 2290
               C(J-1) = C(J)*Z + C(J-1)                                 MAN 2300
  355    CONTINUE                                                       MAN 2310
         IF (N .LE. 10) WRITE(NOUT,350) (C(I), I=1,MM)                  MAN 2320
         DO 360 I = 1,MM1                                               MAN 2330
            IF ((K .GE. N-1) .AND. (I .NE. 1)) GO TO 360                MAN 2340
            Z = ABS(C(I) - A(K+1,I))                                    MAN 2350
            IF (Z .GT. DIFF(I)) DIFF(I) = Z                             MAN 2360
  360    CONTINUE                                                       MAN 2370
  365 CONTINUE                                                          MAN 2380
  370 WRITE(NOUT,375)                                                   MAN 2390
  375 FORMAT(41H  MAXIMUM ABSOLUTE VALUES OF DIFFERENCES )              MAN 2400
      WRITE(NOUT,380) (DIFF(I),I=1,MM1)                                 MAN 2410
  380 FORMAT(5E18.9)                                                    MAN 2420
      WRITE(NOUT,385)                                                   MAN 2430
  385 FORMAT(42H  MAXIMUM ABSOLUTE VALUES OF COEFFICIENTS )             MAN 2440
      IF (ABS(C(1)) .GT. COM(1)) COM(1) = ABS(C(1))                     MAN 2450
      WRITE(NOUT,390) (COM(I),I=1,MM1)                                  MAN 2460
  390 FORMAT(5E18.9)                                                    MAN 2470
  400 CONTINUE                                                          MAN 2480
C                                                                       MAN 2490
C                                                                       MAN 2500
C     TEST OF QUINDF WITH NONEQUIDISTANT KNOTS FOLLOWS.                 MAN 2510
C                                                                       MAN 2520
      WRITE(NOUT,410)                                                   MAN 2530
  410 FORMAT(50H1         TEST OF QUINDF WITH NONEQUIDISTANT KNOTS)     MAN 2540
      N = 5                                                             MAN 2550
      X(1) = -3.0                                                       MAN 2560
      X(2) = -1.0                                                       MAN 2570
      X(3) =  0.0                                                       MAN 2580
      X(4) =  3.0                                                       MAN 2590
      X(5) =  4.0                                                       MAN 2600
      Y(1) =  7.0                                                       MAN 2610
      Y(2) = 11.0                                                       MAN 2620
      Y(3) = 26.0                                                       MAN 2630
      Y(4) = 56.0                                                       MAN 2640
      Y(5) = 29.0                                                       MAN 2650
      B(1) =  2.0                                                       MAN 2660
      B(2) = 15.0                                                       MAN 2670
      B(3) = 10.0                                                       MAN 2680
      B(4) =-27.0                                                       MAN 2690
      B(5) =-30.0                                                       MAN 2700
      M = 3                                                             MAN 2710
      MM = 2*M                                                          MAN 2720
      MM1 = MM - 1                                                      MAN 2730
      WRITE(NOUT,415) N,M                                               MAN 2740
  415 FORMAT(5H-N = ,I3,7H    M =,I2)                                   MAN 2750
      CALL QUINDF(N,X,Y,B,CC,DD,EE,FF)                                  MAN 2760
      DO 420 I = 1,N                                                    MAN 2770
         A(I,1) = Y(I)                                                  MAN 2780
         A(I,2) = B(I)                                                  MAN 2790
         A(I,3) = CC(I)                                                 MAN 2800
         A(I,4) = DD(I)                                                 MAN 2810
         A(I,5) = EE(I)                                                 MAN 2820
         A(I,6) = FF(I)                                                 MAN 2830
  420 CONTINUE                                                          MAN 2840
      DO 430 I = 1,MM1                                                  MAN 2850
         DIFF(I) = 0.0                                                  MAN 2860
         COM(I) = 0.0                                                   MAN 2870
  430 CONTINUE                                                          MAN 2880
      DO 470 K = 1,N                                                    MAN 2890
         DO 435 I = 1,MM                                                MAN 2900
            C(I) = A(K,I)                                               MAN 2910
  435 CONTINUE                                                          MAN 2920
      WRITE(NOUT,440) K                                                 MAN 2930
  440 FORMAT(40H ---------------------------------------,I3,            MAN 2940
     *  45H --------------------------------------------)               MAN 2950
      WRITE(NOUT,445) X(K)                                              MAN 2960
  445 FORMAT(F12.8)                                                     MAN 2970
      IF (K .EQ. N) WRITE(NOUT,455) C(1)                                MAN 2980
      IF (K .EQ. N) GO TO 475                                           MAN 2990
      WRITE(NOUT,455) (C(I),I=1,MM)                                     MAN 3000
      DO 450 I = 1,MM1                                                  MAN 3010
         IF (ABS(A(K,I)) .GT. COM(I)) COM(I) = ABS(A(K,I))              MAN 3020
  450 CONTINUE                                                          MAN 3030
  455 FORMAT(6F16.8)                                                    MAN 3040
      Z = X(K+1) - X(K)                                                 MAN 3050
      DO 460 I = 2,MM                                                   MAN 3060
         DO 460 JJ = I,MM                                               MAN 3070
            J = MM + I - JJ                                             MAN 3080
            C(J-1) = C(J)*Z + C(J-1)                                    MAN 3090
  460 CONTINUE                                                          MAN 3100
      WRITE(NOUT,455) (C(I),I=1,MM)                                     MAN 3110
      DO 465 I = 1,MM1                                                  MAN 3120
         IF ((K .GE. N-1) .AND. (I .NE. 1)) GO TO 465                   MAN 3130
         Z = ABS(C(I) - A(K+1,I))                                       MAN 3140
         IF (Z .GT. DIFF(I)) DIFF(I) = Z                                MAN 3150
  465 CONTINUE                                                          MAN 3160
  470 CONTINUE                                                          MAN 3170
  475 WRITE(NOUT,480)                                                   MAN 3180
  480 FORMAT(41H  MAXIMUM ABSOLUTE VALUES OF DIFFERENCES )              MAN 3190
      WRITE(NOUT,485) (DIFF(I),I=1,MM1)                                 MAN 3200
  485 FORMAT(5E18.9)                                                    MAN 3210
      WRITE(NOUT,490)                                                   MAN 3220
  490 FORMAT(42H  MAXIMUM ABSOLUTE VALUES OF COEFFICIENTS )             MAN 3230
      IF (ABS(C(1)) .GT. COM(1)) COM(1) =ABS(C(1))                      MAN 3240
      WRITE(NOUT,495) (COM(I),I=1,MM1)                                  MAN 3250
  495 FORMAT(5F16.8)                                                    MAN 3260
      M = 3                                                             MAN 3270
      DO 600 N = 10,100,10                                              MAN 3280
      MM = 2*M                                                          MAN 3290
      MM1 = MM - 1                                                      MAN 3300
      P = 0.0                                                           MAN 3310
      DO 500 I = 1,N                                                    MAN 3320
         P = P + ABS(SIN((FLOAT(I)))) + 0.001*FLOAT(I)                  MAN 3330
         X(I) = P                                                       MAN 3340
         Y(I) = COS((FLOAT(I))) - 0.5                                   MAN 3350
         B(I) = COS((FLOAT(2*I))) - 0.5                                 MAN 3360
  500 CONTINUE                                                          MAN 3370
      WRITE(NOUT,510) N,M                                               MAN 3380
  510 FORMAT(5H-N = ,I3,7H    M =,I2)                                   MAN 3390
      CALL QUINDF(N,X,Y,B,CC,DD,EE,FF)                                  MAN 3400
      DO 515 I = 1,N                                                    MAN 3410
         A(I,1) = Y(I)                                                  MAN 3420
         A(I,2) = B(I)                                                  MAN 3430
         A(I,3) = CC(I)                                                 MAN 3440
         A(I,4) = DD(I)                                                 MAN 3450
         A(I,5) = EE(I)                                                 MAN 3460
         A(I,6) = FF(I)                                                 MAN 3470
  515 CONTINUE                                                          MAN 3480
      DO 520 I = 1, MM1                                                 MAN 3490
         DIFF(I) = 0.0                                                  MAN 3500
         COM(I) = 0.0                                                   MAN 3510
  520 CONTINUE                                                          MAN 3520
      DO 565 K = 1,N                                                    MAN 3530
         DO 525 I = 1,MM                                                MAN 3540
            C(I) = A(K,I)                                               MAN 3550
  525    CONTINUE                                                       MAN 3560
         IF (N .GT. 10) GOTO 540                                        MAN 3570
         WRITE(NOUT,530) K                                              MAN 3580
  530 FORMAT(40H ---------------------------------------,I3,            MAN 3590
     *       45H --------------------------------------------)          MAN 3600
         WRITE(NOUT,535) X(K)                                           MAN 3610
  535 FORMAT(F12.8)                                                     MAN 3620
         IF (K .EQ. N) WRITE(NOUT,550) C(1)                             MAN 3630
  540    CONTINUE                                                       MAN 3640
         IF (K .EQ. N) GO TO 570                                        MAN 3650
         IF (N .LE. 10) WRITE(NOUT,550) (C(I), I=1,MM)                  MAN 3660
         DO 545 I = 1,MM1                                               MAN 3670
            IF (ABS(A(K,I)) .GT.  COM(I)) COM(I) = ABS(A(K,I))          MAN 3680
  545    CONTINUE                                                       MAN 3690
  550 FORMAT(6F16.8)                                                    MAN 3700
         Z = X(K+1) - X(K)                                              MAN 3710
         DO 555 I = 2,MM                                                MAN 3720
            DO 555 JJ = I,MM                                            MAN 3730
               J = MM + I - JJ                                          MAN 3740
               C(J-1) = C(J)*Z + C(J-1)                                 MAN 3750
  555    CONTINUE                                                       MAN 3760
         IF (N .LE. 10) WRITE(NOUT,550) (C(I), I=1,MM)                  MAN 3770
         DO 560 I = 1,MM1                                               MAN 3780
            IF ((K .GE. N-1) .AND. (I .NE. 1)) GO TO 560                MAN 3790
            Z = ABS(C(I) - A(K+1,I))                                    MAN 3800
            IF (Z .GT. DIFF(I)) DIFF(I) = Z                             MAN 3810
  560    CONTINUE                                                       MAN 3820
  565 CONTINUE                                                          MAN 3830
  570 WRITE(NOUT,575)                                                   MAN 3840
  575 FORMAT(41H  MAXIMUM ABSOLUTE VALUES OF DIFFERENCES )              MAN 3850
      WRITE(NOUT,580) (DIFF(I),I=1,MM1)                                 MAN 3860
  580 FORMAT(5E18.9)                                                    MAN 3870
      WRITE(NOUT,585)                                                   MAN 3880
  585 FORMAT(42H  MAXIMUM ABSOLUTE VALUES OF COEFFICIENTS )             MAN 3890
      IF (ABS(C(1)) .GT. COM(1)) COM(1) = ABS(C(1))                     MAN 3900
      WRITE(NOUT,590) (COM(I),I=1,MM1)                                  MAN 3910
  590 FORMAT(5E18.9)                                                    MAN 3920
  600 CONTINUE                                                          MAN 3930
C                                                                       MAN 3940
C                                                                       MAN 3950
C     TEST OF QUINAT WITH NONEQUIDISTANT DOUBLE KNOTS FOLLOWS.          MAN 3960
C                                                                       MAN 3970
      WRITE(NOUT,610)                                                   MAN 3980
  610 FORMAT(50H1  TEST OF QUINAT WITH NONEQUIDISTANT DOUBLE KNOTS)     MAN 3990
      N = 5                                                             MAN 4000
      X(1) = -3.                                                        MAN 4010
      X(2) = -3.                                                        MAN 4020
      X(3) = -1.                                                        MAN 4030
      X(4) = -1.                                                        MAN 4040
      X(5) = 0.                                                         MAN 4050
      X(6) = 0.                                                         MAN 4060
      X(7) = 3.                                                         MAN 4070
      X(8) = 3.                                                         MAN 4080
      X(9) = 4.                                                         MAN 4090
      X(10) = 4.                                                        MAN 4100
      Y(1) = 7.                                                         MAN 4110
      Y(2) = 2.                                                         MAN 4120
      Y(3) = 11.                                                        MAN 4130
      Y(4) = 15.                                                        MAN 4140
      Y(5) = 26.                                                        MAN 4150
      Y(6) = 10.                                                        MAN 4160
      Y(7) = 56.                                                        MAN 4170
      Y(8) = -27.                                                       MAN 4180
      Y(9) = 29.                                                        MAN 4190
      Y(10) = -30.                                                      MAN 4200
      M = 3                                                             MAN 4210
      NN = 2*N                                                          MAN 4220
      MM = 2*M                                                          MAN 4230
      MM1 = MM - 1                                                      MAN 4240
      WRITE(NOUT,615) N,M                                               MAN 4250
  615 FORMAT(5H-N = ,I3,7H    M =,I2)                                   MAN 4260
      CALL QUINAT(NN,X,Y,BB,CC,DD,EE,FF)                                MAN 4270
      DO 620 I = 1,NN                                                   MAN 4280
         A(I,1) = Y(I)                                                  MAN 4290
         A(I,2) = BB(I)                                                 MAN 4300
         A(I,3) = CC(I)                                                 MAN 4310
         A(I,4) = DD(I)                                                 MAN 4320
         A(I,5) = EE(I)                                                 MAN 4330
         A(I,6) = FF(I)                                                 MAN 4340
  620 CONTINUE                                                          MAN 4350
      DO 630 I = 1,MM1                                                  MAN 4360
         DIFF(I) = 0.0                                                  MAN 4370
         COM(I) = 0.0                                                   MAN 4380
  630 CONTINUE                                                          MAN 4390
      DO 670 K = 1,NN                                                   MAN 4400
         DO 635 I = 1,MM                                                MAN 4410
            C(I) = A(K,I)                                               MAN 4420
  635 CONTINUE                                                          MAN 4430
      WRITE(NOUT,640) K                                                 MAN 4440
  640 FORMAT(40H ---------------------------------------,I3,            MAN 4450
     *  45H --------------------------------------------)               MAN 4460
      WRITE(NOUT,645) X(K)                                              MAN 4470
  645 FORMAT(F12.8)                                                     MAN 4480
      IF (K .EQ. NN) WRITE(NOUT,655) C(1)                               MAN 4490
      IF (K .EQ. NN) GO TO 675                                          MAN 4500
      WRITE(NOUT,655) (C(I),I=1,MM)                                     MAN 4510
      DO 650 I = 1,MM1                                                  MAN 4520
         IF (ABS(A(K,I)) .GT. COM(I)) COM(I) = ABS(A(K,I))              MAN 4530
  650 CONTINUE                                                          MAN 4540
  655 FORMAT(6F16.8)                                                    MAN 4550
      Z = X(K+1) - X(K)                                                 MAN 4560
      DO 660 I = 2,MM                                                   MAN 4570
         DO 660 JJ = I,MM                                               MAN 4580
            J = MM + I - JJ                                             MAN 4590
            C(J-1) = C(J)*Z + C(J-1)                                    MAN 4600
  660 CONTINUE                                                          MAN 4610
      WRITE(NOUT,655) (C(I),I=1,MM)                                     MAN 4620
      DO 665 I = 1,MM1                                                  MAN 4630
         IF ((K .GE. NN-1) .AND. (I .NE. 1)) GO TO 665                  MAN 4640
         Z = ABS(C(I) - A(K+1,I))                                       MAN 4650
         IF (Z .GT. DIFF(I)) DIFF(I) = Z                                MAN 4660
  665 CONTINUE                                                          MAN 4670
  670 CONTINUE                                                          MAN 4680
  675 WRITE(NOUT,680)                                                   MAN 4690
  680 FORMAT(41H  MAXIMUM ABSOLUTE VALUES OF DIFFERENCES )              MAN 4700
      WRITE(NOUT,685) (DIFF(I),I=1,MM1)                                 MAN 4710
  685 FORMAT(5E18.9)                                                    MAN 4720
      WRITE(NOUT,690)                                                   MAN 4730
  690 FORMAT(42H  MAXIMUM ABSOLUTE VALUES OF COEFFICIENTS )             MAN 4740
      IF (ABS(C(1)) .GT. COM(1)) COM(1) =ABS(C(1))                      MAN 4750
      WRITE(NOUT,695) (COM(I),I=1,MM1)                                  MAN 4760
  695 FORMAT(5F16.8)                                                    MAN 4770
      M = 3                                                             MAN 4780
      DO 800 N = 10,100,10                                              MAN 4790
      NN = 2*N                                                          MAN 4800
      MM = 2*M                                                          MAN 4810
      MM1 = MM - 1                                                      MAN 4820
      P = 0.0                                                           MAN 4830
      DO 700 I = 1,N                                                    MAN 4840
         P = P + ABS(SIN(FLOAT(I)))                                     MAN 4850
         X(2*I-1) = P                                                   MAN 4860
         X(2*I)   = P                                                   MAN 4870
         Y(2*I-1) = COS(FLOAT(I)) - 0.5                                 MAN 4880
         Y(2*I)   = COS(FLOAT(2*I)) - 0.5                               MAN 4890
  700 CONTINUE                                                          MAN 4900
      WRITE(NOUT,710) N,M                                               MAN 4910
  710 FORMAT(5H-N = ,I3,7H    M =,I2)                                   MAN 4920
      CALL QUINAT(NN,X,Y,BB,CC,DD,EE,FF)                                MAN 4930
      DO 715 I = 1,NN                                                   MAN 4940
         A(I,1) = Y(I)                                                  MAN 4950
         A(I,2) = BB(I)                                                 MAN 4960
         A(I,3) = CC(I)                                                 MAN 4970
         A(I,4) = DD(I)                                                 MAN 4980
         A(I,5) = EE(I)                                                 MAN 4990
         A(I,6) = FF(I)                                                 MAN 5000
  715 CONTINUE                                                          MAN 5010
      DO 720 I = 1, MM1                                                 MAN 5020
         DIFF(I) = 0.0                                                  MAN 5030
         COM(I) = 0.0                                                   MAN 5040
  720 CONTINUE                                                          MAN 5050
      DO 765 K = 1,NN                                                   MAN 5060
         DO 725 I = 1,MM                                                MAN 5070
            C(I) = A(K,I)                                               MAN 5080
  725    CONTINUE                                                       MAN 5090
         IF (N .GT. 10)  GOTO 740                                       MAN 5100
         WRITE(NOUT,730) K                                              MAN 5110
  730 FORMAT(40H ---------------------------------------,I3,            MAN 5120
     *       45H --------------------------------------------)          MAN 5130
         WRITE(NOUT,735) X(K)                                           MAN 5140
  735 FORMAT(F12.8)                                                     MAN 5150
         IF (K .EQ. NN) WRITE(NOUT,750) C(1)                            MAN 5160
  740    CONTINUE                                                       MAN 5170
         IF (K .EQ. NN) GO TO 770                                       MAN 5180
         IF (N .LE. 10) WRITE(NOUT,750) (C(I), I=1,MM)                  MAN 5190
         DO 745 I = 1,MM1                                               MAN 5200
            IF (ABS(A(K,I)) .GT.  COM(I)) COM(I) = ABS(A(K,I))          MAN 5210
  745    CONTINUE                                                       MAN 5220
  750 FORMAT(6F16.8)                                                    MAN 5230
         Z = X(K+1) - X(K)                                              MAN 5240
         DO 755 I = 2,MM                                                MAN 5250
            DO 755 JJ = I,MM                                            MAN 5260
               J = MM + I - JJ                                          MAN 5270
               C(J-1) = C(J)*Z + C(J-1)                                 MAN 5280
  755    CONTINUE                                                       MAN 5290
         IF (N .LE. 10) WRITE(NOUT,750) (C(I), I=1,MM)                  MAN 5300
         DO 760 I = 1,MM1                                               MAN 5310
            IF ((K .GE. NN-1) .AND. (I .NE. 1)) GO TO 760               MAN 5320
            Z = ABS(C(I) - A(K+1,I))                                    MAN 5330
            IF (Z .GT. DIFF(I)) DIFF(I) = Z                             MAN 5340
  760    CONTINUE                                                       MAN 5350
  765 CONTINUE                                                          MAN 5360
  770 WRITE(NOUT,775)                                                   MAN 5370
  775 FORMAT(41H  MAXIMUM ABSOLUTE VALUES OF DIFFERENCES )              MAN 5380
      WRITE(NOUT,780) (DIFF(I),I=1,MM1)                                 MAN 5390
  780 FORMAT(5E18.9)                                                    MAN 5400
      WRITE(NOUT,785)                                                   MAN 5410
  785 FORMAT(42H  MAXIMUM ABSOLUTE VALUES OF COEFFICIENTS )             MAN 5420
      IF (ABS(C(1)) .GT. COM(1)) COM(1) = ABS(C(1))                     MAN 5430
      WRITE(NOUT,790) (COM(I),I=1,MM1)                                  MAN 5440
  790 FORMAT(5E18.9)                                                    MAN 5450
  800 CONTINUE                                                          MAN 5460
C                                                                       MAN 5470
C                                                                       MAN 5480
C     TEST OF QUINAT WITH NONEQUIDISTANT KNOTS, ONE DOUBLE KNOT,        MAN 5490
C        ONE TRIPLE KNOT, FOLLOWS.                                      MAN 5500
C                                                                       MAN 5510
      WRITE(NOUT,805)                                                   MAN 5520
      WRITE(NOUT,810)                                                   MAN 5530
  805 FORMAT(51H1         TEST OF QUINAT WITH NONEQUIDISTANT KNOTS,)    MAN 5540
  810 FORMAT(40H             ONE DOUBLE, ONE TRIPLE KNOT)               MAN 5550
      N = 8                                                             MAN 5560
      X(1) = -3.                                                        MAN 5570
      X(2) = -1.                                                        MAN 5580
      X(3) =  -1.                                                       MAN 5590
      X(4) =  0.                                                        MAN 5600
      X(5) =  3.                                                        MAN 5610
      X(6) = 3.                                                         MAN 5620
      X(7) = 3.                                                         MAN 5630
      X(8) = 4.                                                         MAN 5640
      Y(1) =  7.                                                        MAN 5650
      Y(2) = 11.                                                        MAN 5660
      Y(3) = 15.                                                        MAN 5670
      Y(4) = 26.                                                        MAN 5680
      Y(5) = 56.                                                        MAN 5690
      Y(6) = -30.                                                       MAN 5700
      Y(7) =  -7.                                                       MAN 5710
      Y(8) =  29.                                                       MAN 5720
      M = 3                                                             MAN 5730
      MM = 2*M                                                          MAN 5740
      MM1 = MM - 1                                                      MAN 5750
      WRITE(NOUT,815) N,M                                               MAN 5760
  815 FORMAT(5H-N = ,I3,7H    M =,I2)                                   MAN 5770
      CALL QUINAT(N,X,Y,BB,CC,DD,EE,FF)                                 MAN 5780
      DO 820 I = 1,N                                                    MAN 5790
         A(I,1) = Y(I)                                                  MAN 5800
         A(I,2) = BB(I)                                                 MAN 5810
         A(I,3) = CC(I)                                                 MAN 5820
         A(I,4) = DD(I)                                                 MAN 5830
         A(I,5) = EE(I)                                                 MAN 5840
         A(I,6) = FF(I)                                                 MAN 5850
  820 CONTINUE                                                          MAN 5860
      DO 830 I = 1,MM1                                                  MAN 5870
         DIFF(I) = 0.0                                                  MAN 5880
         COM(I) = 0.0                                                   MAN 5890
  830 CONTINUE                                                          MAN 5900
      DO 870 K = 1,N                                                    MAN 5910
         DO 835 I = 1,MM                                                MAN 5920
            C(I) = A(K,I)                                               MAN 5930
  835 CONTINUE                                                          MAN 5940
      WRITE(NOUT,840) K                                                 MAN 5950
  840 FORMAT(40H ---------------------------------------,I3,            MAN 5960
     *  45H --------------------------------------------)               MAN 5970
      WRITE(NOUT,845) X(K)                                              MAN 5980
  845 FORMAT(F12.8)                                                     MAN 5990
      IF (K .EQ. N) WRITE(NOUT,855) C(1)                                MAN 6000
      IF (K .EQ. N) GO TO 875                                           MAN 6010
      WRITE(NOUT,855) (C(I),I=1,MM)                                     MAN 6020
      DO 850 I = 1,MM1                                                  MAN 6030
         IF (ABS(A(K,I)) .GT. COM(I)) COM(I) = ABS(A(K,I))              MAN 6040
  850 CONTINUE                                                          MAN 6050
  855 FORMAT(6F16.8)                                                    MAN 6060
      Z = X(K+1) - X(K)                                                 MAN 6070
      DO 860 I = 2,MM                                                   MAN 6080
         DO 860 JJ = I,MM                                               MAN 6090
            J = MM + I - JJ                                             MAN 6100
            C(J-1) = C(J)*Z + C(J-1)                                    MAN 6110
  860 CONTINUE                                                          MAN 6120
      WRITE(NOUT,855) (C(I),I=1,MM)                                     MAN 6130
      DO 865 I = 1,MM1                                                  MAN 6140
         IF ((K .GE. N-1) .AND. (I .NE. 1)) GO TO 865                   MAN 6150
         Z = ABS(C(I) - A(K+1,I))                                       MAN 6160
         IF (Z .GT. DIFF(I)) DIFF(I) = Z                                MAN 6170
  865 CONTINUE                                                          MAN 6180
  870 CONTINUE                                                          MAN 6190
  875 WRITE(NOUT,880)                                                   MAN 6200
  880 FORMAT(41H  MAXIMUM ABSOLUTE VALUES OF DIFFERENCES )              MAN 6210
      WRITE(NOUT,885) (DIFF(I),I=1,MM1)                                 MAN 6220
  885 FORMAT(5E18.9)                                                    MAN 6230
      WRITE(NOUT,890)                                                   MAN 6240
  890 FORMAT(42H  MAXIMUM ABSOLUTE VALUES OF COEFFICIENTS )             MAN 6250
      IF (ABS(C(1)) .GT. COM(1)) COM(1) =ABS(C(1))                      MAN 6260
      WRITE(NOUT,895) (COM(I),I=1,MM1)                                  MAN 6270
  895 FORMAT(5F16.8)                                                    MAN 6280
C                                                                       MAN 6290
C                                                                       MAN 6300
C     TEST OF QUINAT WITH NONEQUIDISTANT KNOTS, TWO DOUBLE KNOTS,       MAN 6310
C        ONE TRIPLE KNOT,FOLLOWS.                                       MAN 6320
C                                                                       MAN 6330
      WRITE(NOUT,905)                                                   MAN 6340
      WRITE(NOUT,910)                                                   MAN 6350
  905 FORMAT(51H1         TEST OF QUINAT WITH NONEQUIDISTANT KNOTS,)    MAN 6360
  910 FORMAT(40H             TWO DOUBLE, ONE TRIPLE KNOT)               MAN 6370
      N = 10                                                            MAN 6380
      X(1) = 0.                                                         MAN 6390
      X(2) = 2.                                                         MAN 6400
      X(3) = 2.                                                         MAN 6410
      X(4) = 3.                                                         MAN 6420
      X(5) = 3.                                                         MAN 6430
      X(6) = 3.                                                         MAN 6440
      X(7) = 5.                                                         MAN 6450
      X(8) = 8.                                                         MAN 6460
      X(9) = 9.                                                         MAN 6470
      X(10)= 9.                                                         MAN 6480
      Y(1) = 163.                                                       MAN 6490
      Y(2) = 237.                                                       MAN 6500
      Y(3) = -127.                                                      MAN 6510
      Y(4) = 119.                                                       MAN 6520
      Y(5) = -65.                                                       MAN 6530
      Y(6) = 192.                                                       MAN 6540
      Y(7) = 293.                                                       MAN 6550
      Y(8) =  326.                                                      MAN 6560
      Y(9) = 0.                                                         MAN 6570
      Y(10)= -414.0                                                     MAN 6580
      M = 3                                                             MAN 6590
      MM = 2*M                                                          MAN 6600
      MM1 = MM - 1                                                      MAN 6610
      WRITE(NOUT,915) N,M                                               MAN 6620
  915 FORMAT(5H-N = ,I3,7H    M =,I2)                                   MAN 6630
      CALL QUINAT(N,X,Y,BB,CC,DD,EE,FF)                                 MAN 6640
      DO 920 I = 1,N                                                    MAN 6650
         A(I,1) = Y(I)                                                  MAN 6660
         A(I,2) = BB(I)                                                 MAN 6670
         A(I,3) = CC(I)                                                 MAN 6680
         A(I,4) = DD(I)                                                 MAN 6690
         A(I,5) = EE(I)                                                 MAN 6700
         A(I,6) = FF(I)                                                 MAN 6710
  920 CONTINUE                                                          MAN 6720
      DO 930 I = 1,MM1                                                  MAN 6730
         DIFF(I) = 0.0                                                  MAN 6740
         COM(I) = 0.0                                                   MAN 6750
  930 CONTINUE                                                          MAN 6760
      DO 970 K = 1,N                                                    MAN 6770
         DO 935 I = 1,MM                                                MAN 6780
            C(I) = A(K,I)                                               MAN 6790
  935 CONTINUE                                                          MAN 6800
      WRITE(NOUT,940) K                                                 MAN 6810
  940 FORMAT(40H ---------------------------------------,I3,            MAN 6820
     *  45H --------------------------------------------)               MAN 6830
      WRITE(NOUT,945) X(K)                                              MAN 6840
  945 FORMAT(F12.8)                                                     MAN 6850
      IF (K .EQ. N) WRITE(NOUT,955) C(1)                                MAN 6860
      IF (K .EQ. N) GO TO 975                                           MAN 6870
      WRITE(NOUT,955) (C(I),I=1,MM)                                     MAN 6880
      DO 950 I = 1,MM1                                                  MAN 6890
         IF (ABS(A(K,I)) .GT. COM(I)) COM(I) = ABS(A(K,I))              MAN 6900
  950 CONTINUE                                                          MAN 6910
  955 FORMAT(6F16.8)                                                    MAN 6920
      Z = X(K+1) - X(K)                                                 MAN 6930
      DO 960 I = 2,MM                                                   MAN 6940
         DO 960 JJ = I,MM                                               MAN 6950
            J = MM + I - JJ                                             MAN 6960
            C(J-1) = C(J)*Z + C(J-1)                                    MAN 6970
  960 CONTINUE                                                          MAN 6980
      WRITE(NOUT,955) (C(I),I=1,MM)                                     MAN 6990
      DO 965 I = 1,MM1                                                  MAN 7000
         IF ((K .GE. N-1) .AND. (I .NE. 1)) GO TO 965                   MAN 7010
         Z = ABS(C(I) - A(K+1,I))                                       MAN 7020
         IF (Z .GT. DIFF(I)) DIFF(I) = Z                                MAN 7030
  965 CONTINUE                                                          MAN 7040
  970 CONTINUE                                                          MAN 7050
  975 WRITE(NOUT,980)                                                   MAN 7060
  980 FORMAT(41H  MAXIMUM ABSOLUTE VALUES OF DIFFERENCES )              MAN 7070
      WRITE(NOUT,985) (DIFF(I),I=1,MM1)                                 MAN 7080
  985 FORMAT(5E18.9)                                                    MAN 7090
      WRITE(NOUT,990)                                                   MAN 7100
  990 FORMAT(42H  MAXIMUM ABSOLUTE VALUES OF COEFFICIENTS )             MAN 7110
      IF (ABS(C(1)) .GT. COM(1)) COM(1) =ABS(C(1))                      MAN 7120
      WRITE(NOUT,995) (COM(I),I=1,MM1)                                  MAN 7130
  995 FORMAT(5F16.8)                                                    MAN 7140
      STOP                                                              MAN 7150
      END                                                               MAN 7160
