C
C     PROBLEM: TP268
C
*     SET OF INDICES
      IND5=1..5
      IND6=1..6
C
*     TABLE D(I,J), I IN IND6, J IN IND5
      1  1  -74.0
      1  2   80.0
      1  3   18.0
      1  4  -11.0
      1  5   -4.0
      2  1   14.0
      2  2  -69.0
      2  3   21.0
      2  4   28.0
      2  5    0.0
      3  1   66.0
      3  2  -72.0
      3  3   -5.0
      3  4    7.0
      3  5    1.0
      4  1  -12.0
      4  2   66.0
      4  3  -30.0
      4  4  -23.0
      4  5    3.0
      5  1    3.0
      5  2    8.0
      5  3   -7.0
      5  4   -4.0
      5  5    1.0
      6  1    4.0
      6  2  -12.0
      6  3    4.0
      6  4    4.0
      6  5    0.0
C
*     TABLE B(I), I IN IND6
      1  51.0
      2 -61.0
      3 -56.0
      4  69.0
      5  10.0
      6 -12.0
C  
*     VARIABLE
      X(I),I IN IND5
C
*     FUNCTION F
      F = SUM(SUM(D(I,J)*X(J),J IN IND5)
     1     *SUM(D(I,J)*X(J),J IN IND5),I IN IND6)
     2     -2*SUM(B(I)*SUM(D(I,J)*X(J),J IN IND5),I IN IND6)
     3     +SUM(B(I)*B(I),I IN IND6)
C
*     FUNCTION G1
      G1=-X(1)-X(2)-X(3)-X(4)-X(5)+5
C
*     FUNCTION G2
      G2=10*X(1)+10*X(2)-3*X(3)+5*X(4)+4*X(5)-20
C
*     FUNCTION G3
      G3=-8*X(1)+X(2)-2*X(3)-5*X(4)+3*X(5)+40
C
*     FUNCTION G4
      G4=8*X(1)-X(2)+2*X(3)+5*X(4)-3*X(5)-11
C
*     FUNCTION G5
      G5=-4*X(1)-2*X(2)+3*X(3)-5*X(4)+X(5)+30
C
*     END
