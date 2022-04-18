C
C     Problem:  TP370
C
*     SET OF INDICES
      INDN=1..6
      INDN2=2..6
      INDF=1..29
C
*     VARIABLE
      X(I), I IN INDN
C
*     FUNCTION F
      F=X(1)**2 + (X(2) - X(1)**2 - 1.0)**2
     /  + SUM( (SUM ((J-1)*X(J)*(I/29)**(J-2), J IN INDN2)
     /         - (SUM (X(J)*(I/29)**(J-1), J IN INDN) )**2 - 1.0)**2,
     /         I IN INDF)
C
*     END
C
