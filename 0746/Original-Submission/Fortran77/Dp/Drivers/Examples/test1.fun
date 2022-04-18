C
C    Problem:  1
C
*     SET OF INDICES
      INDN=1..20
      INDNM1=1..19
C
*     VARIABLE
      X(I), I in INDN
C
*     FUNCTION F
      F=SUM(100*(X(I+1)-X(I)**2)**2 + (1 - X(I))**2 , I IN INDNM1)
C
C
*     END
C
