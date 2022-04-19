C
C     Problem:  TP295
C
*     INTEGER CONSTANT
      N=60
      N1=N-1
C
*     SET OF INDICES
      IND=1..N
      IND1=1..N1
C
*     VARIABLE
      X(I), I IN IND
C
*     FUNCTION F
      F=SUM(100*(X(I+1) - X(I)**2)**2 + (1 - X(I))**2, I IN IND1)
C
*     END
C
