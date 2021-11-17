      SUBROUTINE GREGORY(N,R,T,W,B,C)
C
C
C  Computes the abscissas and weights of the Gregory quadrature rule with r
C  differences:
C  
C    \int_{t_0}^{t_n} f(t) dt  \approx
C  	  h \left( \frac{1}{2} f_0 + f_1 + \cdots +
C            f_{n-1} + \frac{1}{2} f_n \right ) -
C  	  \frac{h}{12}( \nabla f_n - \delta f_0) -
C  	  \frac{h}{24}( \nabla^2 f_n + \delta^2 f_0) - \cdots -
C  	  h c_{r+1}^{*} (\nabla^r f_n + \delta^r f_0)
C  = \sum_{j=0}^{n} w_j f(t_j),
C  
C  where h = (t_n - t_0)/n, and the c_j^* are given in Henrici (1964). The
C  number r must be an integer from 0 to n, the number of subdivisions. The
C  left and right endpoints must be in t(0) and t(n) respectively. The
C  abscissas are returned in t(0) to t(n) and the corresponding weights in
C  w(0) to w(n).
C  
C  If r=0 the Gregory rule is the same as the repeated trapezoid rule, and if
C  r=n the same as the Newton-Cotes rule (closed type). The order p of the
C  quadrature rule is p = r+1 for r odd and p = r+2 for r even. For n >= 9
C  and large r some of the weights can be negative.
C  
C  For n<= 32 and r<= 24, the numerical integration of powers (less than r)
C  of x on the interval [0,1] gave 9 significant digits correct in an 11
C  digit mantissa.
C  
C  Refs:
C    Hildebrand, F. B. Introduction to Numerical Analysis. McGraw-Hill, New
C    York, 1956, p. 155.
C  
C    Henrici, Peter. Elements of Numerical Analysis. Wiley, New York, 1964,
C    p. 252.
C
C
C     .. Parameters ..
      DOUBLE PRECISION HALF,ONE,ZERO
      PARAMETER (HALF=0.5D0,ONE=1.0D0,ZERO=0.0D0)
C     ..
C     .. Scalar Arguments ..
      INTEGER N,R
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION B(0:N), C(0:N+1),T(0:N),W(0:N)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION CJ,H
      INTEGER I,J
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC DBLE,REAL
C     ..
      B(0) = ONE
      C(0) = ONE
      C(1) = -HALF
      B(N) = ZERO
      H = (T(N)-T(0))/N
      W(0) = HALF
      W(N) = HALF
      DO I = 1,N - 1
          W(I) = ONE
          T(I) = DBLE(REAL(I))*H + T(0)
          B(I) = ZERO
      END DO
      IF (R.GT.N) THEN
          R = N
      END IF

      DO J = 1,R
          CJ = HALF*C(J)
          DO I = J,1,-1
              B(I) = B(I) - B(I-1)
          END DO
          DO I = 3,J + 2
              CJ = CJ + C(J+2-I)/I
          END DO
          C(J+1) = -CJ
          DO I = 0,N
              W(I) = W(I) - CJ* (B(N-I)+B(I))
          END DO
      END DO
      DO I = 0,N
          W(I) = W(I)*H
      END DO
      END
