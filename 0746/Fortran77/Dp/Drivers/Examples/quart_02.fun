C
C     PROBLEM: QUARTIC (TP290)
C
*     PARAMETER
      N=2
C
*     SET OF INDICES
      IND=1..N
C
*     TABLE A(I,J), I IN IND, J IN IND
      1 1 1.0
      2 2 2.0
C
*     VARIABLE
      X(I),I IN IND
C
*     FUNCTION F
      F=(SUM( SUM(A(I,J)*X(I),I IN IND) *X(J),J IN IND))**2
C
*     END
