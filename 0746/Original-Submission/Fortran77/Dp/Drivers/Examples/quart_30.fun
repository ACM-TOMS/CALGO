C
C     PROBLEM: QUARTIC (TP292)
C
*     PARAMETER
      N=30
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
      11 11 11.0
      12 12 12.0
      13 13 13.0
      14 14 14.0
      15 15 15.0
      16 16 16.0
      17 17 17.0
      18 18 18.0
      19 19 19.0
      20 20 20.0
      21 21 21.0
      22 22 22.0
      23 23 23.0
      24 24 24.0
      25 25 25.0
      26 26 26.0
      27 27 27.0
      28 28 28.0
      29 29 29.0
      30 30 30.0
C
*     VARIABLE
      X(I),I IN IND
C
*     FUNCTION F
      F=(SUM( SUM(A(I,J)*X(I),I IN IND) *X(J),J IN IND))**2
C
*     END
