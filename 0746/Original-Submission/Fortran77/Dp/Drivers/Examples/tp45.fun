C
C    Problem:  TP45
C
*     INTEGER CONSTANT
      N=100
C
*     REAL CONSTANT
      FAC = 1.0/120.0
C
*     SET OF INDICES
      INDEXN=1..N
C
*     VARIABLE
      X(I), I in INDEXN
C
*     FUNCTION F
      F=2.0 - FAC*PROD(X(I), I IN INDEXN)
C
*     END
C
