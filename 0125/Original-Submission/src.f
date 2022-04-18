C Original Fortran translation donated by
C
C M.Dow@anu.edu.au
C ANUSF,  Australian National University
C Canberra Australia
C 
C Tidied up to use workspace arrays, dimension arrays to 
C required lengths, use 0: for dimensions of A, add comments,
C and NAG tools to layout source
C
C trh (20/07/97)

      SUBROUTINE WEIGHTCOEFF(N,Q,E,EPS,W,X,WORK)
C
C This is just a wrapper routine to split up the workspace array
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION EPS
      INTEGER N
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION E(N-1),Q(N),W(N),WORK(9*N+8),X(N)
C     ..
C     .. External Subroutines ..
      EXTERNAL WEIGHTC
C     ..
      CALL WEIGHTC(N,Q,E,EPS,W,X,WORK,WORK(N+1))
      END

      SUBROUTINE WEIGHTC(N,Q,E,EPS,W,X,G,A)
C
C Computes the abscissae x(i) and the weight coefficients w(i) for a
C Gaussian quadrature method
C      \int_0^b w(x)f(x) dx \approx \sum_{i=1}^{n}w_if(x_i) where
C \int_0^b w(x) dx = 1 and w(x)>= 0. The method requires the order n, a
C tolerance eps and the 2n-1 first coefficients of the continued fraction
C      \int_0^b {w(x) \over z-x} = { 1| \over |z} - {q_1 | \over |1} -
C 				 {e_1 | \over |z} - {q_2 | \over |1} -
C 				 {e_2 | \over |z} - \cdots
C to be given, the latter in the two arrays q(n) and e(n-1) all
C components of which are automatically positive by virtue of the
C condition w(x)>= 0.  The method works as well if the upper bound b is
C actually infinity (note that b does not appear directly as an
C argument!) or if the density function w(x) dx is replaced by da(x) with
C a monotonically increasing a(x) with at least n points of of variation.
C The tolerance eps should be given in accordance to the machines
C accuracy (preferably by using the value of d1mach(4)). The result is
C delivered as two arrays w(n) (the weight coefficients) and x(n) (the
C abscissae). For a description of the method see H Rutishauser, ``On a
C modification of the QD-algorithm with Graeffe-type convergence''
C [Proceedings of the IFIPS Congress, Munich, 1962].

C     .. Scalar Arguments ..
      DOUBLE PRECISION EPS
      INTEGER N
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION A(0:N,0:7),E(N-1),G(N),Q(N),W(N),X(N)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION M,P
      INTEGER K
C     ..
C     .. External Subroutines ..
      EXTERNAL QDGRAEFFE
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS,EXP,LOG
C     ..
      X(1) = Q(1) + E(1)
      DO K = 2,N
          G(K-1) = E(K-1)*Q(K)/X(K-1)
          IF (K.EQ.N) THEN
              X(K) = Q(K) - G(K-1)

          ELSE
              X(K) = Q(K) + E(K) - G(K-1)
          END IF

          G(K-1) = G(K-1)/X(K)
          W(K-1) = X(K)/X(K-1)
          X(K-1) = LOG(X(K-1))
      END DO
      X(N) = LOG(X(N))
      P = 1
  30  DO K = 1,N - 1
          IF (ABS(G(K)*W(K)).GT.EPS) GO TO 40
      END DO
      GO TO 50

   40 CALL QDGRAEFFE(N,X,G,W,A)
      P = 2*P
      GO TO 30
C
C What follows is a peculiar method to compute the w(k) from 
C the given ratios g_k = w_{k+1}/w_k such that 
C  \sum_{k=1}^n w_k = 1, but the straightforward formulae to do
C this might well produce overflow of exponent
C

   50 W(1) = 1
      M = 0
      DO K = 1,N - 1
          W(K+1) = W(K)*G(K)
          IF (W(K).GT.M) M = W(K)
      END DO
      WRITE (6,FMT=*) 'm=',M
C   /*do k=1,n w(k)=exp(w(k)-m) */
      M = 0
      DO K = 1,N
          M = M + W(K)
      END DO
      DO K = 1,N
          W(K) = W(K)/M
          X(K) = EXP(X(K)/P)
      END DO
      RETURN

      END

      SUBROUTINE RED(A,F,N)
C
C  Subroutine RED reduces a heptadiagonal matrix a to tridiagonal form as
C  described in the paper referenced above.
C  
C     .. Scalar Arguments ..
      INTEGER N
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION A(0:N,0:7),F(N)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION C
      INTEGER J,K
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS
C     ..
      DO K = 1,N - 1
          DO J = K,N - 1

              C = -F(J)*A(J,7)/A(J,2)
              A(J,7) = 0
              A(J+1,2) = A(J+1,2) + C*A(J,5)
              A(J,1) = A(J,1) - C*F(J)*A(J,4)
              A(J,6) = A(J,6) - C*A(J+1,1)
              A(J+1,3) = A(J+1,3) - C*A(J+1,6)
          END DO
          DO J = K,N - 1
              C = -F(J)*A(J,4)/A(J,1)
              A(J,4) = 0
              A(J+1,1) = A(J+1,1) + C*A(J,6)
              A(J+1,6) = A(J+1,6) + C*A(J+1,3)
              A(J,5) = A(J,5) - C*A(J+1,2)
              A(J+1,0) = A(J+1,0) - C*A(J+1,5)
          END DO
          DO J = K + 1,N - 1
              C = -A(J,3)/A(J-1,6)
              A(J,3) = 0
              A(J,6) = A(J,6) + C*A(J,1)
              A(J-1,5) = A(J-1,5) - C*F(J)*F(J)*A(J,0)
              A(J,2) = A(J,2) - C*F(J)*F(J)*A(J,5)
              A(J,7) = A(J,7) - C*F(J)*A(J+1,2)
          END DO
          DO J = K + 1,N - 1
              C = -A(J,0)/A(J-1,5)
              A(J,0) = 0
              A(J+1,2) = A(J+1,2) + C*F(J)*A(J,7)
              A(J,5) = A(J,5) + C*A(J,2)
              A(J,1) = A(J,1) - C*F(J)*F(J)*A(J,6)
              A(J,4) = A(J,4) - C*F(J)*A(J+1,1)
          END DO
      END DO
      RETURN

      END

      SUBROUTINE QDGRAEFFE(N,H,G,F,A)
C
C Subroutine QDGRAEFFE computes for a given continued fraction
C     f(z) = { 1| \over |z} - {q_1 | \over |1} -
C	     {e_1 | \over |z} - {q_2 | \over |1} -
C	     {e_2 | \over |z} - \cdots - {q_n | \over |1}
C another one, the poles of which are the squares of the poles of f(z)
C However QDGRAEFFE uses not the coefficients q_1 ... q_n and
C e_1 ... e_{n-1} of f(z) but the quotients f_k = q_{k+1}/q_k and
C g_k = e_k/q_{k+1} for k=1,2,...,n-1 and the h_k = ln(abs(q_k)) for
C k=1,2,...,n, and the results are delivered in the same form. Routine 
C QDGRAEFFE can be used independently, but requires subroutine RED
C
C     .. Scalar Arguments ..
      INTEGER N
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION A(0:N,0:7),F(N),G(N),H(N)
C     ..
C     .. Local Scalars ..
      INTEGER K
C     ..
C     .. External Subroutines ..
      EXTERNAL RED
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS,LOG
C     ..
      G(N) = 0
      F(N) = 0
      DO K = 1,N
          A(K-1,4) = 1
          A(K-1,5) = 1
          A(K,2) = 1 + G(K)*F(K)
          A(K,1) = 1 + G(K)*F(K)
          A(K,6) = G(K)
          A(K,7) = G(K)
          A(K,0) = 0
          A(K,3) = 0
      END DO
      A(N,5) = 0
C
C The array a represents the heptadiagonal matrix Q of the paper 
C cited above, but with the modifications needed to avoid large numbers
C and with a peculiar arrangement.
C
      CALL RED(A,F,N)
      DO K = 1,N
          H(K) = 2*H(K) + LOG(ABS(A(K,1)*A(K,2)))
      END DO
      DO K = 1,N - 1
          F(K) = F(K)*F(K)*A(K+1,2)*A(K+1,1)/ (A(K,1)*A(K,2))
          G(K) = A(K,5)*A(K,6)/ (A(K+1,1)*A(K+1,2))
      END DO
      RETURN

      END
