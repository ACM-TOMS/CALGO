      SUBROUTINE P2PGS ( A, B, N, FN, GN, VINT )
C
C  THIS SUBROUTINE USES THE PRODUCT TYPE TWO-POINT GAUSS-
C  LEGENDRE-SIMPSON RULE COMPOUNDED N TIMES TO APPROXIMATE
C  THE INTEGRAL FROM A TO B OF THE FUNCTION FN(X) * GN(X).
C  FN AND GN ARE FUNCTION SUBPROGRAMS WHICH MUST BE SUPPLIED
C  BY THE USER.  THE RESULT IS STORED IN VINT.
C
      INTEGER I, N, J, K
      DOUBLE PRECISION A, AG, AM(2,3), B, F(2), FN, G(3),
     &  GN, H, VINT, X(2), Y(2)
      DATA AM(1,1), AM(2,3) / 2 * 2.732050807568877D0 /,
     &  AM(1,2), AM(2,2) / 2 * 4.D0 /, AM(1,3), AM(2,1)
     &  / 2 * -.7320508075688773D0 /
      H = ( B - A ) / DBLE ( FLOAT ( N ) )
      X(1) = A + .2113248654051871D0 * H
      X(2) = A + .7886751345948129D0 * H
      Y(1) = A + H / 2.D0
      Y(2) = A + H
      VINT = 0.D0
      G(3) = GN ( A )
      DO 4 I  =  1, N
        G(1) = G(3)
        DO 1 J = 1, 2
          F(J) = FN ( X(J) )
          G(J+1) = GN ( Y(J) )
          X(J) = X(J) + H
          Y(J) = Y(J) + H
1     CONTINUE
        DO 3 J = 1, 2
          AG = 0.D0
          DO 2 K = 1, 3
            AG = AG + AM(J,K) * G(K)
2     CONTINUE
          VINT = VINT + F(J) * AG
3     CONTINUE
4     CONTINUE
      VINT = H * VINT / 12.D0
      RETURN
      END
