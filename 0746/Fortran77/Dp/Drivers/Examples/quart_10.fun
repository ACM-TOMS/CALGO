C
C     PROBLEM: QUARTIC (TP291)
C
*     PARAMETER
      N=10
C
*     SET OF INDICES
      IND=1..N
C
*     TABLE A(I,J), I IN IND, J IN IND
      1  1  1.0
      2  2  2.0
      3  3  3.0
      4  4  4.0
      5  5  5.0
      6  6  6.0
      7  7  7.0
      8  8  8.0
      9  9  9.0
      10 10 10.0
C
*     VARIABLE
      X(I),I IN IND
C
*     FUNCTION F
      F=(SUM( SUM(A(I,J)*X(I),I IN IND) *X(J),J IN IND))**2
C
*     END
