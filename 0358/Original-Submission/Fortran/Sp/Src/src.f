      SUBROUTINE CSVD ( A, MMAX, NMAX, M, N, P, NU, NV, S, U, V )
C
CC CSVD computes the singular value decomposition of a complex matrix.
C
C  Discussion:
C
C    The singular value decomposition of a complex M by N matrix A
C    has the form
C
C      A = U S V*
C
C    where 
C
C      U is an M by M unitary matrix,
C      S is an M by N diagonal matrix,
C      V is an N by N unitary matrix.
C
C    Moreover, the entries of S are nonnegative and occur on the diagonal
C    in descending order.
C
C    Several internal arrays are dimensioned under the assumption
C    that N <= 100.
C
C  Modified:
C
C    30 October 2006
C
C  Reference:
C
C    Peter Businger, Gene Golub,
C    Algorithm 358:
C    Singular Value Decomposition of a Complex Matrix,
C    Communications of the ACM,
C    Volume 12, Number 10, October 1969, pages 564-565.
C
C  Parameters:
C
C    Input/output, complex A(MMAX,*), the M by N matrix, which may be
C    augmented by P extra columns to which the transformation U*
C    is to be applied.  On output, A has been overwritten, and
C    if 0 < P, columns N+1 through N+P have been premultiplied by U*.
C
C    Input, integer MMAX, the leading dimension of the arrays A
C    and U.
C
C    Input, integer NMAX, the leading dimension of V, and perhaps
C    the second dimension of A and U.
C
C    Input, integer M, N, the number of rows and columns in A.
C    It must be the case that M <= N.  Several internal arrays are
C    dimensioned under the assumption that N <= 100.
C
C    Input, integer P, the number of vectors, stored in A(*,N+1:N+P),
C    to which the transformation U* should be applied.
C
C    Input, integer NU, the number of columns of U to compute.
C
C    Input, integer NV, the number of columns of V to compute.
C
C    Output, real S(N), the computed singular values.
C
C    Output, complex U(MMAX,NU), the first NU columns of U.
C
C    Output, complex V(NMAX,NV), the first NV columns of V.
C
C  Local Parameters:
C
C    Local, real ETA, the relative machine precision.
C    The original text uses ETA = 1.5E-8.
C
C    Local, real TOL, the smallest normalized positive number, divided by ETA.
C    The original test uses TOL = 1.E-31.
C
      IMPLICIT NONE

      INTEGER MMAX
      INTEGER NMAX

      COMPLEX A(MMAX,*)
      REAL B(100)
      REAL C(100)
      REAL CS
      REAL EPS
      REAL ETA
      REAL F
      REAL G
      REAL H
      INTEGER I
      INTEGER J
      INTEGER K
      INTEGER KK
      INTEGER K1
      INTEGER L
      INTEGER L1
      INTEGER LL
      INTEGER M
      INTEGER N
      INTEGER N1
      INTEGER NP
      INTEGER NU
      INTEGER NV
      INTEGER P
      COMPLEX Q
      COMPLEX R
      REAL S(*)
      REAL SN
      REAL T(100)
      REAL TOL
      COMPLEX U(MMAX,*)
      COMPLEX V(NMAX,*)
      REAL W
      REAL X
      REAL Y
      REAL Z

      SAVE ETA
      SAVE TOL

      DATA ETA / 1.5E-07 /
      DATA TOL / 1.5E-31 /

      NP = N + P
      N1 = N + 1
C
C  HOUSEHOLDER REDUCTION.
C
      C(1) = 0.0E+00
      K = 1

10    CONTINUE

      K1 = K + 1
C
C  ELIMINATION OF A(I,K), I = K+1, ..., M.
C
      Z = 0.0E+00
      DO 20 I = K, M
        Z = Z + REAL ( A(I,K) )**2 + AIMAG ( A(I,K) )**2
20    CONTINUE

      B(K) = 0.0E+00

      IF ( Z .LE. TOL ) GO TO 70

      Z = SQRT ( Z )
      B(K) = Z
      W = CABS ( A(K,K) )
      Q = CMPLX ( 1.0E+00, 0.0E+00 )
      IF ( W .NE. 0.0E+00 ) Q = A(K,K) / W
      A(K,K) = Q * ( Z + W )

      IF ( K .EQ. NP ) GO TO 70

      DO 50 J = K1, NP

        Q = CMPLX ( 0.0E+00, 0.0E+00 )
        DO 30 I = K, M
          Q = Q + CONJG ( A(I,K) ) * A(I,J)
30      CONTINUE
        Q = Q / ( Z * ( Z + W ) )

        DO 40 I = K, M
          A(I,J) = A(I,J) - Q * A(I,K)
40      CONTINUE

50    CONTINUE
C
C  PHASE TRANSFORMATION.
C
      Q = -CONJG ( A(K,K) ) / CABS ( A(K,K) )
      DO 60 J = K1, NP
        A(K,J) = Q * A(K,J)
60    CONTINUE
C
C  ELIMINATION OF A(K,J), J = K+2, ..., N
C
70    CONTINUE

      IF ( K .EQ. N ) GO TO 140

      Z = 0.0E+00
      DO 80 J = K1, N
        Z = Z + REAL ( A(K,J) )**2 + AIMAG ( A(K,J) )**2
80    CONTINUE

      C(K1) = 0.0E+00

      IF ( Z .LE. TOL ) GO TO 130

      Z = SQRT ( Z )
      C(K1) = Z
      W = CABS ( A(K,K1) )
      Q = CMPLX ( 1.0E+00, 0.0E+00 )
      IF ( W .NE. 0.0E+00 ) Q = A(K,K1) / W
      A(K,K1) = Q * ( Z + W )

      DO 110 I = K1, M
        Q = CMPLX ( 0.0E+00, 0.0E+00 )
        DO 90 J = K1, N
          Q = Q + CONJG ( A(K,J) ) * A(I,J)
90      CONTINUE
        Q = Q / ( Z * ( Z + W ) )
        DO 100 J = K1, N
          A(I,J) = A(I,J) - Q * A(K,J)
100     CONTINUE
110   CONTINUE
C
C  PHASE TRANSFORMATION.
C
      Q = -CONJG ( A(K,K1) ) / CABS ( A(K,K1) )
      DO 120 I = K1, M
        A(I,K1) = A(I,K1) * Q
120   CONTINUE

130   CONTINUE

      K = K1
      GO TO 10
C
C  TOLERANCE FOR NEGLIGIBLE ELEMENTS.
C
140   CONTINUE

      EPS = 0.0E+00
      DO 150 K = 1, N
        S(K) = B(K)
        T(K) = C(K)
        EPS = AMAX1 ( EPS, S(K) + T(K) )
150   CONTINUE

      EPS = EPS * ETA
C
C  INITIALIZATION OF U AND V.
C
      IF ( NU .EQ. 0 ) GO TO 180

      DO 170 J = 1, NU
        DO 160 I = 1, M
          U(I,J) = CMPLX ( 0.0E+00, 0.0E+00 )
160     CONTINUE
        U(J,J) = CMPLX ( 1.0E+00, 0.0E+00 )
170   CONTINUE

180   CONTINUE

      IF ( NV .EQ. 0 )GO TO 210

      DO 200 J = 1, NV
        DO 190 I = 1, N
          V(I,J) = CMPLX ( 0.0E+00, 0.0E+00 )
190     CONTINUE
        V(J,J) = CMPLX ( 1.0E+00, 0.0E+00 )
200   CONTINUE
C
C  QR DIAGONALIZATION.
C
210   CONTINUE

      DO 380 KK = 1, N

        K = N1 - KK
C
C  TEST FOR SPLIT.
C
220     CONTINUE

        DO 230 LL = 1, K
          L = K + 1 - LL
          IF ( ABS ( T(L) ) .LE. EPS ) GO TO 290
          IF ( ABS ( S(L-1) ) .LE. EPS ) GO TO 240
230     CONTINUE
C
C  CANCELLATION OF E(L).
C
240     CONTINUE

        CS = 0.0E+00
        SN = 1.0E+00
        L1 = L - 1

        DO 280 I = L, K

          F = SN * T(I)
          T(I) = CS * T(I)

          IF ( ABS ( F ) .LE. EPS ) GO TO 290

          H = S(I)
          W = SQRT ( F * F + H * H )
          S(I) = W
          CS = H / W
          SN = - F / W

          IF ( NU .EQ. 0 ) GO TO 260

          DO 250 J = 1, N
            X = REAL ( U(J,L1) )
            Y = REAL ( U(J,I) )
            U(J,L1) = CMPLX ( X * CS + Y * SN, 0.0E+00 )
            U(J,I)  = CMPLX ( Y * CS - X * SN, 0.0E+00 )
250       CONTINUE

260       CONTINUE

          IF ( NP .EQ. N ) GO TO 280

          DO 270 J = N1, NP
            Q = A(L1,J)
            R = A(I,J)
            A(L1,J) = Q * CS + R * SN
            A(I,J)  = R * CS - Q * SN
270       CONTINUE

280     CONTINUE
C
C  TEST FOR CONVERGENCE.
C
290     CONTINUE

        W = S(K)
        IF ( L .EQ. K ) GO TO 360
C
C  ORIGIN SHIFT.
C
        X = S(L)
        Y = S(K-1)
        G = T(K-1)
        H = T(K)
        F = ( ( Y - W ) * ( Y + W ) + ( G - H ) * ( G + H ) ) 
     &    / ( 2.0E+00 * H * Y )
        G = SQRT ( F * F + 1.0E+00 )
        IF ( F .LT. 0.0E+00 ) G = -G
        F = ( ( X - W ) * ( X + W ) + ( Y / ( F + G ) - H ) * H ) / X
C
C  QR STEP.
C
        CS = 1.0E+00
        SN = 1.0E+00
        L1 = L + 1

        DO 350 I = L1, K

          G = T(I)
          Y = S(I)
          H = SN * G
          G = CS * G
          W = SQRT ( H * H + F * F )
          T(I-1) = W
          CS = F / W
          SN = H / W
          F = X * CS + G * SN
          G = G * CS - X * SN
          H = Y * SN
          Y = Y * CS

          IF ( NV .EQ. 0 ) GO TO 310

          DO 300 J = 1, N
            X = REAL ( V(J,I-1) )
            W = REAL ( V(J,I) )
            V(J,I-1) = CMPLX ( X * CS + W * SN, 0.0E+00 )
            V(J,I)   = CMPLX ( W * CS - X * SN, 0.0E+00 )
300       CONTINUE

310       CONTINUE

          W = SQRT ( H * H + F * F )
          S(I-1) = W
          CS = F / W
          SN = H / W
          F = CS * G + SN * Y
          X = CS * Y - SN * G

          IF ( NU .EQ. 0 ) GO TO 330

          DO 320 J = 1, N
            Y = REAL ( U(J,I-1) )
            W = REAL ( U(J,I) )
            U(J,I-1) = CMPLX ( Y * CS + W * SN, 0.0E+00 )
            U(J,I)   = CMPLX ( W * CS - Y * SN, 0.0E+00 )
320       CONTINUE

330       CONTINUE

          IF ( N .EQ. NP ) GO TO 350

          DO 340 J = N1, NP
            Q = A(I-1,J)
            R = A(I,J)
            A(I-1,J) = Q * CS + R * SN
            A(I,J)   = R * CS - Q * SN
340       CONTINUE

350     CONTINUE

        T(L) = 0.0E+00
        T(K) = F
        S(K) = X
        GO TO 220
C
C  CONVERGENCE.
C
360     CONTINUE

        IF ( W .GE. 0.0E+00 ) GO TO 380
        S(K) = -W
        IF ( NV .EQ. 0 ) GO TO 380

        DO 370 J = 1, N
          V(J,K) = -V(J,K)
370     CONTINUE

380   CONTINUE
C
C  SORT SINGULAR VALUES.
C
      DO 450 K = 1, N

        G = -1.0E+00
        J = K

        DO 390 I = K, N
          IF ( S(I) .LE. G ) GO TO 390
          G = S(I)
          J = I
390     CONTINUE

        IF ( J .EQ. K ) GO TO 450

        S(J) = S(K)
        S(K) = G
C
C  Interchange V(1:N,J) and V(1:N,K).
C
        IF ( NV .EQ. 0 ) GO TO 410

        DO 400 I = 1, N
          Q      = V(I,J)
          V(I,J) = V(I,K)
          V(I,K) = Q
400     CONTINUE

410     CONTINUE
C
C  Interchange U(1:N,J) and U(1:N,K).
C
        IF ( NU .EQ. 0 ) GO TO 430

        DO 420 I = 1, N
          Q      = U(I,J)
          U(I,J) = U(I,K)
          U(I,K) = Q
420     CONTINUE

430     CONTINUE
C
C  Interchange A(J,N1:NP) and A(K,N1:NP).
C
        IF ( N .EQ. NP ) GO TO 450

        DO 440 I = N1, NP
          Q      = A(J,I)
          A(J,I) = A(K,I)
          A(K,I) = Q
440     CONTINUE

450   CONTINUE
C
C  BACK TRANSFORMATION.
C
      IF ( NU .EQ. 0 ) GO TO 510

      DO 500 KK = 1, N

        K = N1 - KK
        IF ( B(K) .EQ. 0.0E+00 ) GO TO 500
        Q = -A(K,K) / CABS ( A(K,K) )

        DO 460 J = 1, NU
          U(K,J) = Q * U(K,J)
460     CONTINUE

        DO 490 J = 1, NU

          Q = CMPLX ( 0.0E+00, 0.0E+00 )

          DO 470 I = K, M
            Q = Q + CONJG ( A(I,K) ) * U(I,J)
470       CONTINUE

          Q = Q / ( CABS ( A(K,K) ) * B(K) )

          DO 480 I = K, M
            U(I,J) = U(I,J) - Q * A(I,K)
480       CONTINUE

490     CONTINUE

500   CONTINUE

510   CONTINUE

      IF ( NV .EQ. 0 ) GO TO 570
      IF ( N .LT. 2 ) GO TO 570

      DO 560 KK = 2, N

        K = N1 - KK
        K1 = K + 1
        IF ( C(K1) .EQ. 0.0E+00 ) GO TO 560
        Q = -CONJG ( A(K,K1) ) / CABS ( A(K,K1) )

        DO 520 J = 1, NV
          V(K1,J) = Q * V(K1,J)
520     CONTINUE

        DO 550 J = 1, NV

          Q = CMPLX ( 0.0E+00, 0.0E+00 )

          DO 530 I = K1, N
            Q = Q + A(K,I) * V(I,J)
530       CONTINUE

          Q = Q / ( CABS ( A(K,K1) ) * C(K1) )

          DO 540 I = K1, N
            V(I,J) = V(I,J) - Q * CONJG ( A(K,I) )
540       CONTINUE

550     CONTINUE

560   CONTINUE

570   RETURN
      END
