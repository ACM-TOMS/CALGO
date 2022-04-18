C
C     HELMHOLTZ ENERGY FUNCTION
C
*     PARAMETER
      N=100
C
*     SET OF INDICES
      INDEX=1..N
C
*     REAL CONSTANT
      R=8.314
      T=273.0
      C1=1.0 + DSQRT(2.0)
      C2=1.0 - DSQRT(2.0)
      C3=DSQRT(8.0)
      A(I,J)=1/(I+J-1), I IN INDEX, J IN INDEX
      B(I)=0.00001, I IN INDEX
C
*     VARIABLE
      X(I), I IN INDEX
C
*     FUNCTION F
      HBX=SUM(B(I)*X(I), I IN INDEX)
      XAX=SUM(X(I)*SUM(A(I,J)*X(J), J IN INDEX), I IN INDEX)
      F=R*T*SUM(X(I)*DLOG(X(I)/(1-HBX)),I IN INDEX)
     1  -XAX*DLOG((1+C1*HBX)/(1+C2*HBX))/(C3*HBX)
C
*     END
C
