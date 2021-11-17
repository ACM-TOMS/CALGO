      SUBROUTINE DECOMP(N, NDIM, A, IP)
      INTEGER N, NDIM, K, KP1, M, I, J
      REAL A(NDIM,NDIM),T
      INTEGER IP(NDIM)
C
C  MATRIX TRIANGULARIZATION BY GAUSSIAN ELIMINATION.
C  INPUT..
C    N = ORDER OF MATRIX.
C    NDIM = DECLARED DIMENSION OF ARRAY A.
C    A = MATRIX TO BE TRIANGULARIZED.
C  OUTPUT..
C    A(I,J), I.LE.J = UPPER TRIANGULAR FACTOR U.
C    A(I,J), I.GT.J = MULTIPLIERS = LOWER TRIANGULAR
C                     FACTOR, I-L.
C    IP(K), K.LT.N = INDEX OF K-TH PIVOT ROW.
C    IP(N) = (-1)**(NUMBER OF INTERCHANGES) OR 0.
C  USE 'SOLVE' TO OBTRAIN SOLUTION OF LINEAR SYSTEM.
C  DETERM(A) = IP(N)*A(1,1)*A(2,2)*...*A(N,N).
C  IF IP(N) = 0, A IS SINGULAR, SOLVE WILL DIVIDE BY ZERO.
C  INTERCHANGES FINISHED IN U, ONLY PARTLY IN L.
C
      IP(N) = 1
      DO 6 K = 1,N
        IF (K.EQ.N) GO TO 5
        KP1 = K+1
        M = K
        DO 1 I = KP1,N
          IF(ABS(A(I,K)).GT.ABS(A(M,K)))M=I
1       CONTINUE
        IP(K) = M
        IF(M.NE.K) IP(N) = -IP(N)
        T = A(M,K)
        A(M,K) = A(K,K)
        A(K,K) = T
        IF(T.EQ.0.) GO TO 5
        DO 2 I = KP1,N
          A(I,K) = -A(I,K)/T
2       CONTINUE
        DO 4 J = KP1,N
          T = A(M,J)
          A(M,J) = A(K,J)
          A(K,J) = T
          IF(T.EQ.0.) GO TO 4
          DO 3 I = KP1,N
            A(I,J) = A(I,J) + A(I,K)*T
3         CONTINUE
4       CONTINUE
5       IF(A(K,K).EQ.0.) IP(N) = 0
6     CONTINUE
      RETURN
      END

      SUBROUTINE SOLVE(N, NDIM, A, B, IP)
      INTEGER N, NDIM, NM1, K, KP1, M, I, KB, KM1
      REAL A(NDIM,NDIM),B(NDIM),T
      INTEGER IP(NDIM)
C
C  SOLUTION OF LINEAR SYSTEM, A*X = B.
C  INPUT..
C    N = ORDER 0F MATRIX.
C    NDIM = DECLARED DIMENSION OF ARRAY A.
C    A = TRIANGULARIZED MATRIX OBTAINED FROM 'DECOMP'.
C    B = RIGHT HAND SIDE VECTOR.
C    IP = PIVOT VECTOR OBTAINED FROM 'DECOMP'.
C  DO NOT USE IF DECOMP HAS SET IP(N) = 0.
C  OUTPUT..
C    B = SOLUTION VECTOR, X.
C
      IF(N.EQ.1) GO TO 9
      NM1 = N-1
      DO 75 K = 1,NM1
        KP1 = K+1
        M = IP(K)
        T = B(M)
        B(M) = B(K)
        B(K) = T
        DO 7 I = KP1, N
          B(I) = B(I) + A(I,K)*T
7       CONTINUE
75    CONTINUE
      DO 85 KB = 1,NM1
        KM1 = N-KB
        K = KM1+1
        B(K) = B(K)/A(K,K)
        T = -B(K)
        DO 8 I = 1,KM1
          B(I) = B(I) + A(I,K)*T
8       CONTINUE
85    CONTINUE
9     B(1) = B(1)/A(1,1)
      RETURN
      END
