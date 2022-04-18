C
C     ALGORITHM AS 311: APPLIED STATISTICS (1997), VOL. 46 NO. 1
C
      SUBROUTINE ELF1( M, IM, P, Q, N, W, PHI, THETA, QQ, ISMU, MU,
     *                 ATF, A, SIGMA2, XITOL, LOGELF, F1, F2,
     *                 WS, NWS, IWS, NIWS, IFAULT )
C
      INTEGER    M, IM, P, Q, N, NWS, NIWS, IWS(NIWS), IFAULT,
     *           G, G2, M2, I1, I2, I3, I4, I5, I6, I7, I8,
     *           I9, I10, I11, BIG1, BIG2, BIG3, TOTAL
      LOGICAL    ISMU, ATF
      DOUBLE PRECISION
     *           W(IM,N), PHI(IM,P*M+1), THETA(IM,Q*M+1), QQ(IM,M),
     *           MU(M), A(IM,N), SIGMA2, XITOL, LOGELF, F1, F2,
     *           WS(NWS), EPSIL
C
      INTRINSIC  MAX
      EXTERNAL   MACHEP, ELF2
C
      IFAULT = 0
      IF ( M .LT. 1 ) THEN
         IFAULT = 1
         RETURN
      ELSEIF ( N .LT. 1 ) THEN
         IFAULT = 2
         RETURN
      ELSEIF ( P .LT. 0 ) THEN
         IFAULT = 3
         RETURN
      ELSEIF ( Q .LT. 0 ) THEN
         IFAULT = 4
         RETURN
      ELSEIF ( P .EQ. 0 .AND. Q .EQ. 0 ) THEN
         IFAULT = 5
         RETURN
      ELSEIF (IM .LT. M) THEN
         IFAULT = 6
         RETURN
      ENDIF
C
C     [2]: Check that workspace is big enough
C
      G  = MAX( P, Q )
      G2 = G * G
      M2 = M * M
      IF ( P .GT. 0 ) THEN
         BIG1 = M * (M + 1) / 2 + M2 * (P - 1)
      ELSE
         BIG1 = 1
      ENDIF
      BIG2 = MAX( BIG1, G * M )
      BIG3 = MAX( N, Q )
C
      I1  = 1
      I2  = I1  + M2
      I3  = I2  + M2
      I4  = I3  + M2
      I5  = I4  + BIG2
      I6  = I5  + M
      I7  = I6  + BIG1 * BIG1
      I8  = I7  + M2 * G2
      I9  = I8  + M2 * G2
      I10 = I9  + M2 * G2
      I11 = I10 + M2 * (P + Q) * G
C
      TOTAL  = I11 + M2 * BIG3 - 1
C
      IF ( NWS .LT. TOTAL ) THEN
         IFAULT = 7
         RETURN
      ELSEIF ( NIWS .LT. BIG1 ) THEN
         IFAULT = 8
         RETURN
      ENDIF
C
C     [3]: Calculate machine epsilon (store as EPSIL)
C
      CALL MACHEP( EPSIL )
C
C     [4]: Call ELF2 to evaluate the exact log likelihood
C
      CALL ELF2( M, IM, P, Q, N, G, W, PHI, THETA, QQ, ISMU, MU,
     *           ATF, A, SIGMA2, XITOL, LOGELF, F1, F2, EPSIL,
     *           BIG1, BIG2, BIG3, WS(I1), WS(I2), WS(I3),
     *           WS(I4), WS(I5), WS(I6), WS(I7), WS(I8), WS(I9),
     *           WS(I10), WS(I11), IWS, IFAULT )
C
      END
C
C
C
      SUBROUTINE ELF2( M, IM, P, Q, N, G, W, PHI, THETA, QQ, ISMU, MU,
     *                 ATF, A, SIGMA2, XITOL, LOGELF, F1, F2, EPS,
     *                 BIG1, BIG2, BIG3, Q1, Q1INV, MTMP4,
     *                 VTMP1, VTMP2, MATPHI, MTMP0, MTMP2, MTMP3,
     *                 MTMP1, GAMXI, INDX, IFAULT )
C
      INTEGER    M, IM, P, Q, G, N, BIG1, BIG2, BIG3, INDX(BIG1),
     *           IFAULT, I, II, J, JJ, JL, K, KK, NLIM
      LOGICAL    ISMU, ATF
      DOUBLE PRECISION
     *           W(IM,N), PHI(IM,P*M+1), THETA(IM,Q*M+1), QQ(IM,M),
     *           MU(M), A(IM,N), SIGMA2, XITOL, LOGELF, F1, F2, EPS,
     *           Q1(M,M), Q1INV(M,M), MTMP4(M,M), VTMP1(BIG2),
     *           VTMP2(M), MATPHI(BIG1,BIG1), MTMP0(G*M,G*M),
     *           MTMP2(G*M,G*M), MTMP3(G*M,G*M), MTMP1(M*(P+Q),M*G),
     *           GAMXI(M,M*BIG3), D1, D2, S1, S2, DETQ, DETOM,
     *           ZERO, ONE, TWO, LG2PI, DBLE
C The value of log(2*pi) is corrected here compared to the version of jam197
C available on Jose A. Mauricio's home page, http://www.ucm.es/info/ecocuan/jam
C (it is only single precision there).
      PARAMETER  ( ZERO = 0.0, ONE = 1.0, TWO = 2.0, LG2PI = 
     *                                             1.83787706640935D0)
      INTRINSIC  EXP, LOG, MAX, DBLE
      EXTERNAL   CHOLDC, CHOLFR, CGAMMA, CXI, CRES
C
C     Copy lower triangle of QQ into upper triangles of QQ and Q1
C
      DO 2 I = 1, M
      DO 1 J = I, M
      QQ(I,J) = QQ(J,I)
      Q1(I,J) = QQ(J,I)
    1 CONTINUE
    2 CONTINUE
C
C     [1]: Calculate the inverse of the Cholesky factor of QQ
C          (store as Q1INV) and the determinant of QQ (store as DETQ)
C
      CALL CHOLDC( Q1, M, M, D1, D2, EPS, IFAULT )
C
      IF ( IFAULT .GT. 0 ) THEN
         IFAULT = 9
         RETURN
      ENDIF
C
      DO 5 I = 1, M
      DO 3 J = 1, M
      VTMP1(J) = ZERO
      Q1INV(J,I) = ZERO
    3 CONTINUE
      VTMP1(I) = ONE
      CALL CHOLFR( Q1, M, M, VTMP1 )
      DO 4 J = I, M
      Q1INV(J,I) = VTMP1(J)
    4 CONTINUE
    5 CONTINUE
C
      DETQ = D1 * TWO ** D2
C
C     [2]: Calculate the theoretical autocovariance and
C          cross-covariance matrices (store as GAMXI and VTMP1)
C
      IF ( P .GT. 0 ) THEN
C
         CALL CGAMMA( M, IM, P, Q, PHI, THETA, QQ, BIG1, MATPHI,
     *                MTMP2, MTMP3, GAMXI, VTMP1, INDX, IFAULT )
         IF ( IFAULT .GT. 0 ) THEN
            IFAULT = 10
            RETURN
         ENDIF
C
      ENDIF
C
C     [3]: Calculate M: Cholesky factor of V1*OMEGA*V1' (store as MTMP0)
C
      DO 20 I = 1, M * G
      DO 10 J = 1, M * G
      MTMP0(I,J) = ZERO
   10 CONTINUE
   20 CONTINUE
C
C     [3.1]: Calculate OMEGA*V1' (store as MTMP1)
C
      DO 40 I = 1, M * (P + Q)
      DO 30 J = 1, M * G
      MTMP1(I,J) = ZERO
   30 CONTINUE
   40 CONTINUE
C
      DO 140 I = 1, P
      DO 130 J = 1, G
C
      DO 80 K = J - I, P - I
      DO 70 II = 1, M
      DO 60 JJ = 1, M
      S1 = ZERO
      DO 50 KK = 1, M
C
      IF ( K .GT. 0 ) THEN
         JL = M * (M + 1) / 2 + M * M * (K - 1) + M * (KK - 1) + II
      ELSEIF ( K .LT. 0 ) THEN
         JL = M * (M + 1) / 2 + M * M * (- K - 1) + M * (II - 1) + KK
      ELSE
         IF ( KK .GE. II ) THEN
            JL = KK * (KK - 1) / 2 + II
         ELSE
            JL = II * (II - 1) / 2 + KK
         ENDIF
      ENDIF
C
      S1 = S1 + VTMP1(JL) * PHI(JJ, (P - K - I + J - 1) * M + KK)
   50 CONTINUE
      MTMP1(II + (I - 1) * M, JJ + (J - 1) * M) =
     *      MTMP1(II + (I - 1) * M, JJ + (J - 1) * M) + S1
   60 CONTINUE
   70 CONTINUE
   80 CONTINUE
C
      DO 120 K = J - I, Q - I
      IF ( P + K .LE. Q ) THEN
         DO 110 II = 1, M
         DO 100 JJ = 1, M
         S1 = ZERO
         DO 90 KK = 1, M
         S1 = S1 + GAMXI(II, (Q - P - K) * M + KK)
     *           * THETA(JJ, (Q - K - I + J - 1) * M + KK)
   90    CONTINUE
         MTMP1(II + (I - 1) * M, JJ + (J - 1) * M) =
     *         MTMP1(II + (I - 1) * M, JJ + (J - 1) * M) - S1
  100    CONTINUE
  110    CONTINUE
      ENDIF
  120 CONTINUE
C
  130 CONTINUE
  140 CONTINUE
C
      DO 230 I = P + 1, P + Q
      DO 220 J = 1, G
C
      DO 180 K = P + J - I, P + P - I
      IF ( P - K .LE. Q ) THEN
         DO 170 II = 1, M
         DO 160 JJ = 1, M
         S1 = ZERO
         DO 150 KK = 1, M
         S1 = S1 + GAMXI(KK, (Q - P + K) * M + II)
     *           * PHI(JJ, (P + P - K - I + J - 1) * M + KK)
  150    CONTINUE
         MTMP1(II + (I - 1) * M, JJ + (J - 1) * M) =
     *         MTMP1(II + (I - 1) * M, JJ + (J - 1) * M) + S1
  160    CONTINUE
  170    CONTINUE
      ENDIF
  180 CONTINUE
C
      IF ( P - I + J .LE. 0 ) THEN
         DO 210 II = 1, M
         DO 200 JJ = 1, M
         S1 = ZERO
         DO 190 KK = 1, M
         S1 = S1 + QQ(II,KK)
     *           * THETA(JJ, (Q + P - I + J - 1) * M + KK)
  190    CONTINUE
         MTMP1(II + (I - 1) * M, JJ + (J - 1) * M) =
     *         MTMP1(II + (I - 1) * M, JJ + (J - 1) * M) - S1
  200    CONTINUE
  210    CONTINUE
      ENDIF
C
  220 CONTINUE
  230 CONTINUE
C
C     [3.2]: Calculate V1*OMEGA*V1' (store as MTMP0)
C
      DO 330 I = 1, G
      DO 320 J = I, G
C
      DO 270 K = 0, P - I
      DO 260 II = 1, M
      IF ( I .EQ. J ) THEN
         JL = II
      ELSE
         JL = 1
      ENDIF
      DO 250 JJ = JL, M
      S1 = ZERO
      DO 240 KK = 1, M
      S1 = S1 + PHI(II, (P - K - 1) * M + KK)
     *        * MTMP1(KK + (K + I - 1) * M, JJ + (J - 1) * M)
  240 CONTINUE
      MTMP0(II + (I - 1) * M, JJ + (J - 1) * M) =
     *      MTMP0(II + (I - 1) * M, JJ + (J - 1) * M) + S1
  250 CONTINUE
  260 CONTINUE
  270 CONTINUE
C
      DO 310 K = 0, Q - I
      DO 300 II = 1, M
      IF ( I .EQ. J ) THEN
         JL = II
      ELSE
         JL = 1
      ENDIF
      DO 290 JJ = JL, M
      S1 = ZERO
      DO 280 KK = 1, M
      S1 = S1 + THETA(II, (Q - K - 1) * M + KK)
     *        * MTMP1(KK + (K + P + I - 1) * M, JJ + (J - 1) * M)
  280 CONTINUE
      MTMP0(II + (I - 1) * M, JJ + (J - 1) * M) =
     *      MTMP0(II + (I - 1) * M, JJ + (J - 1) * M) - S1
  290 CONTINUE
  300 CONTINUE
  310 CONTINUE
C
  320 CONTINUE
  330 CONTINUE
C
C     [3.3]: Calculate M (overwrite MTMP0)
C
      CALL CHOLDC( MTMP0, M * G, M * G, D1, D2, EPS, IFAULT )
C
      IF ( IFAULT .GT. 0 ) THEN
         IFAULT = 11
         RETURN
      ENDIF
C
C     [4]: Calculate matrix polynomial R*XI(k) (overwrite GAMXI)
C
      CALL CXI( M, IM, N, Q, THETA, XITOL,
     *           Q1INV, NLIM, GAMXI, MTMP4, IFAULT )
C
      IF ( IFAULT .GT. 0 ) THEN
         IFAULT = 12
         RETURN
      ENDIF
C
C     [5]: Calculate vector eta (store as A)
C
      DO 340 I = 1, M
      DO 335 J = 1, N
      A(I,J) = ZERO
  335 CONTINUE
  340 CONTINUE
C
C     [5.1]: Calculate conditional residuals recursively (store as A)
C
      DO 440 I = 1, N
C
      DO 350 J = 1, M
      VTMP1(J) = ZERO
  350 CONTINUE
      DO 380 J = 1, P
      IF ( I - J .GE. 1 ) THEN
         DO 370 II = 1, M
         S1 = ZERO
         DO 360 K = 1, M
         IF ( ISMU ) THEN
            S2 = W(K, I - J) - MU(K)
         ELSE
            S2 = W(K, I - J)
         ENDIF
         S1 = S1 + PHI(II, (J - 1) * M + K) * S2
  360    CONTINUE
         VTMP1(II) = VTMP1(II) + S1
  370    CONTINUE
      ENDIF
  380 CONTINUE
C
      DO 390 J = 1, M
      VTMP2(J) = ZERO
  390 CONTINUE
      DO 420 J = 1, Q
      IF ( I - J .GE. 1 ) THEN
         DO 410 II = 1, M
         S1 = ZERO
         DO 400 K = 1, M
         S1 = S1 + THETA(II, (J - 1) * M + K) * A(K, I - J)
  400    CONTINUE
         VTMP2(II) = VTMP2(II) + S1
  410    CONTINUE
      ENDIF
  420 CONTINUE
C
      DO 430 II = 1, M
      IF ( ISMU ) THEN
         S2 = W(II,I) - MU(II)
      ELSE
         S2 = W(II,I)
      ENDIF
      A(II,I) = S2 - VTMP1(II) + VTMP2(II)
  430 CONTINUE
C
  440 CONTINUE
C
C     [5.2]: Premultiply each M-block of A by Q1INV (overwrite A)
C
      DO 480 I = 1, N
      DO 460 J = 1, M
      S1 = ZERO
      DO 450 K = 1, J
      S1 = S1 + Q1INV(J,K) * A(K,I)
  450 CONTINUE
      VTMP1(J) = S1
  460 CONTINUE
      DO 470 J = 1, M
      A(J,I) = VTMP1(J)
  470 CONTINUE
  480 CONTINUE
C
C     [6]: Calculate vector M'h (store as VTMP1)
C
      DO 490 I = 1, G * M
      VTMP1(I) = ZERO
  490 CONTINUE
C
C     [6.1]: Calculate vector h (overwrite VTMP1)
C
      DO 530 J = 1, G
      DO 520 I = 0, N - J
      IF ( I .LE. NLIM ) THEN
         DO 510 JJ = 1, M
         S1 = ZERO
         DO 500 K = 1, M
         S1 = S1 + GAMXI(K, I * M + JJ) * A(K,I+J)
  500    CONTINUE
         VTMP1(JJ + (J - 1) * M) = VTMP1( JJ + (J - 1) * M) + S1
  510    CONTINUE
      ENDIF
  520 CONTINUE
  530 CONTINUE
C
C     [6.2]: Premultiply VTMP1 by MTMP0' (overwrite VTMP1)
C
      DO 550 I = 1, M * G
      S1 = ZERO
      DO 540 K = I, M * G
      S1 = S1 + MTMP0(K,I) * VTMP1(K)
  540 CONTINUE
      VTMP1(I) = S1
  550 CONTINUE
C
C     Store M as MTMP3 if residuals have been requested
C
      IF ( ATF ) THEN
         DO 570 I = 1, M * G
         DO 560 J = 1, I
         MTMP3(I,J) = MTMP0(I,J)
  560    CONTINUE
  570    CONTINUE
      ENDIF
C
C     [7]: Calculate H'H (store as MTMP2)
C
      DO 590 I = 1, G * M
      DO 580 J = 1, G * M
      MTMP2(I,J) = ZERO
  580 CONTINUE
  590 CONTINUE
C
      DO 640 I = 1, G
      DO 630 K = 0, N - I
      IF ( K + I - 1 .LE. NLIM ) THEN
C
         DO 620 II = 1, M
         IF ( I .EQ. 1 ) THEN
            JL = II
         ELSE
            JL = M
         ENDIF
         DO 610 JJ = 1, JL
         S1 = ZERO
         DO 600 KK = 1, M
         S1 = S1 + GAMXI(KK, K * M + II)
     *           * GAMXI(KK, (K + I - 1) * M + JJ)
  600    CONTINUE
         MTMP2(II + (I - 1) * M, JJ) = MTMP2(II + (I - 1) * M, JJ) + S1
  610    CONTINUE
  620    CONTINUE
C
      ENDIF
  630 CONTINUE
  640 CONTINUE
C
      DO 690 I = 2, G
      DO 680 J = 2, I
      DO 670 II = 1, M
      IF ( I .EQ. J ) THEN
         JL = II
      ELSE
         JL = M
      ENDIF
      DO 660 JJ = 1, JL
      S1 = ZERO
      IF ( (N - I + 1 .LE. NLIM) .AND. (N - J + 1 .LE. NLIM) ) THEN
C
         DO 650 KK = 1, M
         S1 = S1 + GAMXI(KK, (N - I + 1) * M + II)
     *           * GAMXI(KK, (N - J + 1) * M + JJ)
  650    CONTINUE
C
      ENDIF
      MTMP2(II + (I - 1) * M, JJ + (J - 1) * M) =
     *      MTMP2(II + (I - 2) * M, JJ + (J - 2) * M) - S1
  660 CONTINUE
  670 CONTINUE
  680 CONTINUE
  690 CONTINUE
C
      DO 710 I = 1, G * M
      DO 700 J = I + 1, G * M
      MTMP2(I,J) = MTMP2(J,I)
  700 CONTINUE
  710 CONTINUE
C
C     [8]: Calculate I+M'H'HM and its Cholesky factor (overwrite MTMP0)
C
C     [8.1]: Calculate M'H'H (store as MTMP1)
C
      DO 740 I = 1, G * M
      DO 730 J = 1, G * M
      S1 = ZERO
      DO 720 K = I, G * M
      S1 = S1 + MTMP0(K,I) * MTMP2(K,J)
  720 CONTINUE
      MTMP1(I,J) = S1
  730 CONTINUE
  740 CONTINUE
C
C     [8.2]: Store M as MTMP2 (overwrite) and initialize MTMP0
C
      DO 760 I = 1, G * M
      DO 750 J = 1, G * M
      MTMP2(I,J) = MTMP0(I,J)
      MTMP0(I,J) = ZERO
  750 CONTINUE
  760 CONTINUE
C
C     [8.3]: Calculate I + M'H'HM (store as MTMP0)
C
      DO 790 I = 1, G * M
      DO 780 J = I, G * M
      S1 = ZERO
      DO 770 K = J, G * M
      S1 = S1 + MTMP1(I,K) * MTMP2(K,J)
  770 CONTINUE
      MTMP0(I,J) = S1
  780 CONTINUE
      MTMP0(I,I) = ONE + MTMP0(I,I)
  790 CONTINUE
C
C     [8.4]: Compute the Cholesky factor of I+M'H'HM
C            (overwrite MTMP0) and its determinant (store as DETOM)
C
      CALL CHOLDC( MTMP0, M * G, M * G, D1, D2, EPS, IFAULT )
C
      IF ( IFAULT .GT. 0 ) THEN
         IFAULT = 13
         RETURN
      ENDIF
C
      DETOM = D1 * TWO ** D2
C
C     [9]: Calculate LAMBDA using forward substitution (store as VTMP1)
C
      CALL CHOLFR( MTMP0, M * G, M * G, VTMP1 )
C
C     [10]: Calculate the sum of squares (return as F1)
C
      S1 = ZERO
      DO 810 I = 1, M
      DO 800 J = 1, N
      S1 = S1 + A(I,J) * A(I,J)
  800 CONTINUE
  810 CONTINUE
C
      S2 = ZERO
      DO 820 I = 1, M * G
      S2 = S2 + VTMP1(I) * VTMP1(I)
  820 CONTINUE
C
      F1 = S1 - S2
C
C     Calculate the determinant (return as F2)
C
      F2 = DETOM * (DETQ ** DBLE( N ))
C
C     Calculate the exact log likelihood (return as LOGELF)
C
      LOGELF = - ( DBLE( N * M ) * (LG2PI + LOG( SIGMA2 )) +
     *             DBLE( N ) * LOG( DETQ ) + LOG ( DETOM ) +
     *             F1 / SIGMA2 ) / TWO
C
C     [11]: Calculate residual vector if requested (return as A)
C
      IF ( ATF )
     *   CALL CRES( M, IM, N, G, NLIM, GAMXI,
     *              Q1, MTMP3, MTMP0, VTMP1, A )
C
      END
C
C
C
      SUBROUTINE CGAMMA( M, IM, P, Q, PHI, THETA, QQ, BIG1, MAT,
     *                   WZERO, MZERO, GAMWA, RHS, INDX, IFAULT )
C
      INTEGER    M, IM, P, Q, BIG1, INDX(BIG1), IFAULT, ROW, COL,
     *           H, I, II, J, JJ, K, L, R, S
      DOUBLE PRECISION
     *           PHI(IM,P*M+1), THETA(IM,Q*M+1), QQ(IM,M),
     *           MAT(BIG1,BIG1), WZERO(M,M), MZERO(M,M),
     *           GAMWA(M,Q*M+1), RHS(BIG1), SUM, ZERO, ONE
      PARAMETER  ( ZERO = 0.0, ONE = 1.0 )
      EXTERNAL   DECOMP, SOLVE
C
      IFAULT = 0
C
C     [1]: Compute the Q - 1 cross-covariance matrices (return as GAMWA)
C
      IF ( Q .GE. 1 ) THEN
         DO 20 I = 1, M
         DO 10 J = 1, M
         GAMWA(I,J) = QQ(I,J)
   10    CONTINUE
   20    CONTINUE
      ENDIF
C
      DO 80 K = 1, Q - 1
      DO 70 I = 1, M
      DO 60 J = 1, M
      SUM = ZERO
      DO 30 H = 1, M
      SUM = SUM - THETA(I, (K - 1) * M + H) * QQ(H,J)
   30 CONTINUE
      DO 50 L = 1, K
      IF ( L .LE. P ) THEN
         DO 40 H = 1, M
         SUM = SUM + PHI(I, (L - 1) * M + H)
     *             * GAMWA(H, (K - L) * M + J)
   40    CONTINUE
      ENDIF
   50 CONTINUE
      GAMWA(I, K * M + J) = SUM
   60 CONTINUE
   70 CONTINUE
   80 CONTINUE
C
C     [2]: Compute diagonal and upper triangle of W(0) (store as WZERO)
C
      DO 100 I = 1, M
      DO  90 J = 1, M
      WZERO(I,J) = ZERO
   90 CONTINUE
  100 CONTINUE
C
      DO 180 I = 1, P
      DO 170 J = I, Q
      IF ( J - I .GE. 0 ) THEN
C
         DO 130 II = 1, M
         DO 120 JJ = 1, M
         SUM = ZERO
         DO 110 K = 1, M
         SUM = SUM + PHI(II, (I - 1) * M + K)
     *             * GAMWA(K, (J - I) * M + JJ)
  110    CONTINUE
         MZERO(II,JJ) = SUM
  120    CONTINUE
  130    CONTINUE
C
         DO 160 II = 1, M
         DO 150 JJ = 1, M
         SUM = ZERO
         DO 140 K = 1, M
         SUM = SUM + MZERO(II,K) * THETA(JJ, (J - 1) * M + K)
  140    CONTINUE
         WZERO(II,JJ) = WZERO(II,JJ) + SUM
  150    CONTINUE
  160    CONTINUE
C
      ENDIF
  170 CONTINUE
  180 CONTINUE
C
      DO 200 I = 1, M
      DO 190 J = I, M
      WZERO(I,J) = QQ(I,J) - WZERO(I,J) - WZERO(J,I)
  190 CONTINUE
  200 CONTINUE
C
      DO 270 J = 1, Q
C
      DO 230 II = 1, M
      DO 220 JJ = 1, M
      SUM = ZERO
      DO 210 K = 1, M
      SUM = SUM + THETA(II, (J - 1) * M + K) * QQ(K,JJ)
  210 CONTINUE
      MZERO(II,JJ) = SUM
  220 CONTINUE
  230 CONTINUE
C
      DO 260 II =  1, M
      DO 250 JJ = II, M
      SUM = ZERO
      DO 240 K = 1, M
      SUM = SUM + MZERO(II,K) * THETA(JJ, (J - 1) * M + K)
  240 CONTINUE
      WZERO(II,JJ) = WZERO(II,JJ) + SUM
  250 CONTINUE
  260 CONTINUE
C
  270 CONTINUE
C
C     [3]: Set up system of equations (store as MAT and RHS)
C
      DO 290 I = 1, BIG1
      DO 280 J = 1, BIG1
      MAT(I,J) = ZERO
  280 CONTINUE
      RHS(I) = ZERO
  290 CONTINUE
C
C     [3.1]: Compute the first M * (M + 1) / 2 rows
C
      DO 390 J = 1, M
      DO 380 I = 1, J
      ROW = J * (J - 1) / 2 + I
C
C     [3.1.1]: Compute the first M * (M + 1) / 2 columns
C
      DO 330 L = 1, M
      DO 320 K = 1, L
      COL = L * (L - 1) / 2 + K
      SUM = ZERO
      IF ( K .EQ. L ) THEN
         DO 300 R = 1, P
         SUM = SUM - PHI(I, (R - 1) * M + K) * PHI(J, (R - 1) * M + L)
  300    CONTINUE
      ELSE
         DO 310 R = 1, P
         SUM = SUM - PHI(I, (R - 1) * M + K) * PHI(J, (R - 1) * M + L)
     *             - PHI(I, (R - 1) * M + L) * PHI(J, (R - 1) * M + K)
  310    CONTINUE
      ENDIF
      MAT(ROW,COL) = SUM
  320 CONTINUE
  330 CONTINUE
C
C     [3.1.2]: Compute the remaining M * M * (P - 1) columns
C
      DO 370 S = 1, P - 1
      DO 360 L = 1, M
      DO 350 K = 1, M
      COL = M * (M + 1) / 2 + M * M * (S - 1) + M * (L - 1) + K
      SUM = ZERO
      DO 340 R = 1, P - S
      SUM = SUM - PHI(I, (R + S - 1) * M + K) * PHI(J, (R - 1) * M + L)
     *          - PHI(J, (R + S - 1) * M + K) * PHI(I, (R - 1) * M + L)
  340 CONTINUE
      MAT(ROW,COL) = SUM
  350 CONTINUE
  360 CONTINUE
  370 CONTINUE
C
C     [3.1.3]: Set up RHS and diagonal of MAT
C
      RHS(ROW) = WZERO(I,J)
      MAT(ROW,ROW) = ONE + MAT(ROW,ROW)
C
  380 CONTINUE
  390 CONTINUE
C
C     [3.2]: Compute the remaining M * M * (P - 1) rows
C
      DO 470 S = 1, P - 1
C
      DO 460 I = 1, M
      DO 450 J = 1, M
      ROW = M * (M + 1) / 2 + M * M * (S - 1) + M * (I - 1) + J
C
C     [3.2.1]: Compute the first M * (M + 1) / 2 columns
C
      DO 400 L = 1, M
      IF ( L .LE. J ) THEN
         COL = J * (J - 1) / 2 + L
      ELSE
         COL = L * (L - 1) / 2 + J
      ENDIF
      MAT(ROW,COL) = - PHI(I, (S - 1) * M + L)
  400 CONTINUE
C
C     [3.2.2]: Compute the remaining M * M * (P - 1) columns
C
      DO 420 R = 1, P - 1
      DO 410 L = 1, M
      COL = M * (M + 1) / 2 + (R - 1) * M * M + (J - 1) * M + L
      IF ( R + S .LE. P ) MAT(ROW,COL) = - PHI(I, (R + S - 1) * M + L)
      IF ( S .GT. R ) THEN
         COL = M * (M + 1) / 2 + (R - 1) * M * M + (L - 1) * M + J
         MAT(ROW,COL) = MAT(ROW,COL) - PHI(I, (S - R - 1) * M + L)
      ENDIF
  410 CONTINUE
  420 CONTINUE
C
C     [3.2.3]: Set up RHS and diagonal of MAT
C
      RHS(ROW) = ZERO
      DO 440 II = S, Q
      DO 430  K = 1, M
      RHS(ROW) = RHS(ROW) - GAMWA(J, (II - S) * M + K)
     *                    * THETA(I, (II - 1) * M + K)
  430 CONTINUE
  440 CONTINUE
C
      MAT(ROW,ROW) = ONE + MAT(ROW,ROW)
C
  450 CONTINUE
  460 CONTINUE
  470 CONTINUE
C
C     [4]: Solve for autocovariance matrices (return as RHS)
C
      CALL DECOMP( BIG1, BIG1, MAT, INDX )
      IF ( INDX(BIG1) .EQ. 0 ) THEN
         IFAULT = 1
         RETURN
      ELSE
         CALL SOLVE( BIG1, BIG1, MAT, RHS, INDX )
      ENDIF
C
      END
C
C
C
      SUBROUTINE CXI( M, IM, N, Q, THETA, XITOL,
     *                R, NLIM, XI, MTMP, IFAULT )
C
      INTEGER    M, IM, N, Q, NLIM, IFAULT, H, I, II, J, JJ, K, NQ
      DOUBLE PRECISION
     *           THETA(IM,M*Q+1), XITOL, R(M,M), XI(M,M*N), MTMP(M,M),
     *           S1, S2, MX
      LOGICAL    DELTA
      PARAMETER  ( ZERO = 0.0, ONE = 1.0 )
      INTRINSIC  ABS
C
      IFAULT = 0
      DELTA = .FALSE.
      MX = ZERO
C
      DO 20 I = 1, M
      DO 10 J = 1, M * N
      XI(I,J) = ZERO
   10 CONTINUE
      XI(I,I) = ONE
   20 CONTINUE
C
C     [1]: Update index NLIM and calculate matrix sequence (store as XI)
C
      NLIM = 0
C
   30 IF ( (.NOT. DELTA) .AND. (NLIM .LT. N - 1) ) THEN
C
         NLIM = NLIM + 1
C
C        [1.1]: Calculate the XI matrix for this round
C
         DO 70 J = 1, Q
         IF ( NLIM .GE. J ) THEN
            DO 60 II = 1, M
            DO 50 JJ = 1, M
            S1 = ZERO
            DO 40 H = 1, M
            S1 = S1 + THETA(II, (J - 1) * M + H)
     *              * XI(H, (NLIM - J) * M + JJ)
   40       CONTINUE
            XI(II, NLIM * M + JJ) = XI(II, NLIM * M + JJ) + S1
   50       CONTINUE
   60       CONTINUE
         ENDIF
   70    CONTINUE
C
         S2 = ZERO
         DO 90 II = 1, M
         DO 80 JJ = 1, M
         S2 = S2 + ABS( XI(II, NLIM * M + JJ) )
   80    CONTINUE
   90    CONTINUE
C
         IF ( NLIM .LE. Q ) MX = MX + S2
         IF ( S2 .GT. MX ) THEN
            IFAULT = 1
            RETURN
         ENDIF
C
C        [1.2]: Check for effective convergence
C
         IF ( S2 .LT. XITOL ) THEN
            NQ = 1
            DELTA = .TRUE.
C
  100       IF ((NQ .LE. Q) .AND. (NLIM .LT. N - 1) .AND. (DELTA)) THEN
C
               NQ = NQ + 1
               NLIM = NLIM + 1
               DO 140 J = 1, Q
               IF ( NLIM .GE. J ) THEN
                  DO 130 II = 1, M
                  DO 120 JJ = 1, M
                  S1 = ZERO
                  DO 110 H = 1, M
                  S1 = S1 + THETA(II, (J - 1) * M + H)
     *                    * XI(H, (NLIM - J) * M + JJ)
  110             CONTINUE
                  XI(II, NLIM * M + JJ) = XI(II, NLIM * M + JJ) + S1
  120             CONTINUE
  130             CONTINUE
               ENDIF
  140          CONTINUE
C
               S2 = ZERO
               DO 160 II = 1, M
               DO 150 JJ = 1, M
               S2 = S2 + ABS( XI(II, NLIM * M + JJ) )
  150          CONTINUE
  160          CONTINUE
C
               IF ( NLIM .LE. Q ) MX = MX + S2
               IF ( S2 .GT. MX ) THEN
                  IFAULT = 1
                  RETURN
               ENDIF
               IF ( S2 .GT. XITOL ) DELTA = .FALSE.
               GOTO 100
            ENDIF
C
            IF ( DELTA ) NLIM = NLIM - NQ
         ENDIF
C
         GOTO 30
      ENDIF
C
C     [2]: Premultiply every XI by R (overwrite XI)
C
      DO 220 K = 0, NLIM
C
      DO 190 I = 1, M
      DO 180 J = 1, M
      S1 = ZERO
      DO 170 H = 1, I
      S1 = S1 + R(I,H) * XI(H, K * M + J)
  170 CONTINUE
      MTMP(I,J) = S1
  180 CONTINUE
  190 CONTINUE
C
      DO 210 I = 1, M
      DO 200 J = 1, M
      XI(I, K * M + J) = MTMP(I,J)
  200 CONTINUE
  210 CONTINUE
C
  220 CONTINUE
C
      END
C
C
C
      SUBROUTINE CRES( M, IM, N, G, NLIM, XI,
     *                 Q1, MATM, MATL, LAMBDA, RES )
C
      INTEGER    M, IM, N, G, NLIM, H, I, J, JJ
      DOUBLE PRECISION
     *           XI(M,M*N), Q1(M,M), MATM(G*M,G*M), MATL(G*M,G*M),
     *           LAMBDA(G*M), RES(IM,N), SUM, ZERO
      PARAMETER  ( ZERO = 0.0 )
      EXTERNAL   CHOLBK
C
C     [1]: Solve for C in the system L'C = LAMBDA (overwrite LAMBDA)
C
      CALL CHOLBK( MATL, M * G, M * G, LAMBDA )
C
C     [2]: Calculate D = MC (overwrite LAMBDA)
C
      DO 20 I = M * G, 1, -1
      SUM = ZERO
      DO 10 J = 1, I
      SUM = SUM + MATM(I,J) * LAMBDA(J)
   10 CONTINUE
      LAMBDA(I) = SUM
   20 CONTINUE
C
C     [3]: Calculate residuals (return as RES)
C
      DO 60 I = 1, N
      DO 50 J = 1, I
C
      IF ( (I - J .LE. NLIM) .AND. (J .LE. G) ) THEN
         DO 40 JJ = 1, M
         SUM = ZERO
         DO 30 H = 1, M
         SUM = SUM + XI(JJ, (I - J) * M + H) * LAMBDA(H + (J - 1) * M)
   30    CONTINUE
         RES(JJ,I) = RES(JJ,I) - SUM
   40    CONTINUE
      ENDIF
C
   50 CONTINUE
   60 CONTINUE
C
      DO 90 J = 1, N
      DO 80 I = M, 1, -1
      SUM = ZERO
      DO 70 H = 1, I
      SUM = SUM + Q1(I,H) * RES(H,J)
   70 CONTINUE
      RES(I,J) = SUM
   80 CONTINUE
   90 CONTINUE
C
      END
C
C
C
      SUBROUTINE CHOLDC( M, NDIM, N, D1, D2, EPS, IFAULT )
C
      INTEGER    NDIM, N, IFAULT, I, J, K
      DOUBLE PRECISION
     *           M(NDIM,NDIM), D1, D2, EPS, SUM, ML1, ML2, MLJ,
     *           MXO, MXA, ONE, ZERO, FOUR, P0625, SIXTEN
      PARAMETER  ( ZERO = 0.0, ONE = 1.0, FOUR = 4.0,
     *             P0625 = 0.0625D0, SIXTEN = 16.0 )
      INTRINSIC  ABS, SQRT
C
      IFAULT = 0
      D1 = ONE
      D2 = ZERO
C
C     [1]: Initialize finite arithmetic parameters
C
      ML1 = ZERO
      MXO = SQRT( ABS( M(1,1) ) )
      DO 10 J = 2, N
      SUM = SQRT( ABS( M(J,J) ) )
      IF ( SUM .GT. MXO ) MXO = SUM
   10 CONTINUE
C
      IF ( (MXO * MXO) .LE. SQRT( EPS ) ) THEN
         DO 30 I = 1, N
         DO 20 J = 1, I
         M(I,J) = ZERO
   20    CONTINUE
   30    CONTINUE
         RETURN
      ENDIF
C
      ML2 = SQRT( EPS ) * MXO
      MXA = ZERO
C
C     [2]: Calculate modified Cholesky decomposition
C
      DO 100 J = 1, N
      SUM = M(J,J)
      DO 40 I = 1, J - 1
      SUM = SUM - M(J,I) * M(J,I)
   40 CONTINUE
C
      IF ( ( SUM .NE. ABS( SUM ) ) .AND. ( ABS( SUM ) .GT. ML2 ) ) THEN
         IFAULT = 1
         RETURN
      ELSE
         M(J,J) = SUM
      ENDIF
C
      MLJ = ZERO
C
      DO 60 I = J + 1, N
      SUM = M(J,I)
      DO 50 K = 1, J - 1
      SUM = SUM - M(I,K) * M(J,K)
   50 CONTINUE
      M(I,J) = SUM
      IF ( ABS( M(I,J) ) .GT. MLJ ) MLJ = ABS( M(I,J) )
   60 CONTINUE
C
      IF ( (MLJ / MXO) .GT. ML1 ) THEN
         MLJ = MLJ / MXO
      ELSE
         MLJ = ML1
      ENDIF
C
      IF ( M(J,J) .GT. (MLJ * MLJ) ) THEN
         M(J,J) = SQRT( M(J,J) )
      ELSE
         IF ( MLJ .LT. ML2 ) MLJ = ML2
         IF ( MXA .LT. (MLJ * MLJ - M(J,J)) ) MXA = MLJ * MLJ - M(J,J)
         M(J,J) = MLJ
      ENDIF
C
      D1 = D1 * M(J,J) * M(J,J)
   70 IF ( D1 .GE. ONE ) THEN
         D1 = D1 * P0625
         D2 = D2 + FOUR
         GOTO 70
      ENDIF
   80 IF ( D1 .LT. P0625 ) THEN
         D1 = D1 * SIXTEN
         D2 = D2 - FOUR
         GOTO 80
      ENDIF
C
      DO 90 I = J + 1, N
      M(I,J) = M(I,J) / M(J,J)
   90 CONTINUE
  100 CONTINUE
C
      DO 120 J = 2, N
      DO 110 I = 1, J - 1
      M(I,J) = ZERO
  110 CONTINUE
  120 CONTINUE
C
      END
C
C
C
      SUBROUTINE CHOLFR( MATL, NDIM, N, RHSOL )
C
      INTEGER    NDIM, N, I, J
      DOUBLE PRECISION
     *           MATL(NDIM,NDIM), RHSOL(NDIM), SUM, ZERO
      PARAMETER  ( ZERO = 0.0 )
C
      RHSOL(1) = RHSOL(1) / MATL(1,1)
C
      DO 20 I = 2, N
      SUM = ZERO
      DO 10 J = 1, I - 1
      SUM = SUM + MATL(I,J) * RHSOL(J)
   10 CONTINUE
      RHSOL(I) = (RHSOL(I) - SUM) / MATL(I,I)
   20 CONTINUE
C
      END
C
C
C
      SUBROUTINE CHOLBK( MATL, NDIM, N, RHSOL )
C
      INTEGER    NDIM, N, I, J
      DOUBLE PRECISION
     *           MATL(NDIM,NDIM), RHSOL(NDIM), SUM, ZERO
      PARAMETER  ( ZERO = 0.0 )
C
      RHSOL(N) = RHSOL(N) / MATL(N,N)
C
      DO 20 I = N - 1, 1, -1
      SUM = ZERO
      DO 10 J = I + 1, N
      SUM = SUM + MATL(J,I) * RHSOL(J)
   10 CONTINUE
      RHSOL(I) = (RHSOL(I) - SUM) / MATL(I,I)
   20 CONTINUE
C
      END
C
C
C
      SUBROUTINE MACHEP( EPSIL )
C
      DOUBLE PRECISION
     *           EPSIL, ONE, TWO
      PARAMETER  ( ONE = 1.0, TWO = 2. 0 )
C
      EPSIL = ONE
C
   10 IF ( EPSIL + ONE .GT. ONE ) THEN
         EPSIL = EPSIL / TWO
         GOTO 10
      ENDIF
      EPSIL = TWO * EPSIL
C
      END
C
C     Subroutines DECOMP and SOLVE are those of Moler (1972)
C
      SUBROUTINE DECOMP( N, NDIM, A, IP )
C
      DOUBLE PRECISION
     *           A(NDIM,NDIM), T
      INTEGER    IP(NDIM)
C
      IP(N) = 1
      DO 6 K = 1, N
         IF ( K .EQ. N ) GOTO 5
         KP1 = K + 1
         M = K
         DO 1 I = KP1, N
            IF ( ABS( A(I,K) ) .GT. ABS( A(M,K) ) ) M = I
    1    CONTINUE
         IP(K) = M
         IF ( M .NE. K ) IP(N) = - IP(N)
         T = A(M,K)
         A(M,K) = A(K,K)
         A(K,K) = T
         IF ( T .EQ. 0 ) GOTO 5
         DO 2 I = KP1, N
    2       A(I,K) = - A(I,K) / T
         DO 4 J = KP1, N
            T = A(M,J)
            A(M,J) = A(K,J)
            A(K,J) = T
            IF ( T .EQ. 0 ) GOTO 4
            DO 3 I = KP1, N
    3          A(I,J) = A(I,J) + A(I,K) * T
    4    CONTINUE
    5    IF ( A(K,K) .EQ. 0 ) IP(N) = 0
    6 CONTINUE
      RETURN
      END
C
C
C
      SUBROUTINE SOLVE( N, NDIM, A, B, IP )
C
      DOUBLE PRECISION
     *           A(NDIM,NDIM), B(NDIM), T
      INTEGER    IP(NDIM)
C
      IF ( N .EQ. 1 ) GOTO 9
      NM1 = N - 1
      DO 7 K = 1, NM1
         KP1 = K + 1
         M = IP(K)
         T = B(M)
         B(M) = B(K)
         B(K) = T
         DO 7 I = KP1, N
    7 B(I) = B(I) + A(I,K) * T
      DO 8 K8 = 1, NM1
         KM1 = N - K8
         K = KM1 + 1
         B(K) = B(K) / A(K,K)
         T = - B(K)
         DO 8 I = 1, KM1
    8 B(I) = B(I) + A(I,K) * T
    9 B(1) = B(1) / A(1,1)
      RETURN
      END
C
C
