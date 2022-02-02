      SUBROUTINE PTRAP ( A, B, N, FN, GN, VINT )
C
C  THIS SUBROUTINE USES THE PRODUCT TYPE TRAPEZOIDAL RULE
C  COMPOUNDED N TIMES TO APPROXIMATE THE INTEGRAL FROM A TO B
C  OF THE FUNCTION FN(X) * GN(X).  FN AND GN ARE FUNCTION
C  SUBPROGRAMS WHICH MUST BE SUPPLIED BY THE USER.  THE
C  RESULT IS STORED IN VINT.
C
      DOUBLE PRECISION A, AG, AM(2,2), B, F(2), FN, G(2),
     &  GN, H, VINT, X
      INTEGER I, N, J, K
      DATA AM(1,1), AM(2,2) /2 * 2.D0 /, AM(1,2), AM(2,1)
     &  / 2 * 1.D0/
      H = ( B - A ) / DBLE ( FLOAT ( N ) )
      VINT = 0.D0
      X = A
      F(2) = FN ( A )
      G(2) = GN ( A )
      DO 3 I  =  1, N
        F(1) = F(2)
        G(1) = G(2)
        X = X + H
        F(2) = FN ( X )
        G(2) = GN ( X )
        DO 2 J = 1, 2
          AG = 0.D0
          DO 1 K = 1, 2
            AG = AG + AM(J,K) * G(K)
1     CONTINUE
          VINT = VINT + F(J) * AG
2     CONTINUE
3     CONTINUE
      VINT = H * VINT / 6.D0
      RETURN
      END

