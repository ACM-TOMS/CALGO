      SUBROUTINE FNOVA ( Y, Z, ROW, MSIZE, NCLS, NFCTR )
C
      IMPLICIT NONE

      INTEGER NCLS
      INTEGER NFCTR

      INTEGER I
      INTEGER J
      INTEGER K
      INTEGER KL1
      INTEGER L
      INTEGER MSIZE(NFCTR)
      INTEGER NF
      INTEGER NRNC
      REAL ROW(*)
      REAL Y(NCLS)
      REAL Z(NCLS)
C
C  LOOP FOR NFCTR CONTRAST MATRICES.
C
      DO 5 NF = 1, NFCTR

        I = 1
C
C  GET SIZE OF THE MATRIX.
C
        K = NFCTR - NF + 1
        NRNC = MSIZE(K)

        DO 3 J = 1, NRNC
C
C  ROW OF A CONTRAST MATRIX.
C
          CALL AROW ( ROW, NRNC, J )
C
C  PERFORM THE 'PAPER STRIP' OPERATION FOR A MATRIX ROW.
C
          DO 2 K = 1, NCLS, NRNC

            Z(I) = 0.0E+00

            DO 1 L = 1, NRNC
              KL1 = K + L - 1
              Z(I) = Z(I) + ROW(L) * Y(KL1)
1     CONTINUE

            I = I + 1
2     CONTINUE

3       CONTINUE
C
C  MOVE Z INTO Y.
C
        DO 4 J = 1, NCLS
          Y(J) = Z(J)
4     CONTINUE

5     CONTINUE

      DO 6 J = 1, NCLS
        Y(J) = Y(J) * Y(J)
6     CONTINUE

      RETURN
      END
      SUBROUTINE AROW ( ROW, NRNC, J )

      IMPLICIT NONE

      INTEGER NRNC

      REAL A
      REAL EL
      INTEGER I
      INTEGER J
      INTEGER JM1
      REAL RJ
      REAL ROW(NRNC)
C
C  IF ROW ONE:
C
      IF ( J - 1 ) 3, 1, 3
1     A = NRNC
      EL = 1.0E+00 / SQRT ( A )

      DO 2 I = 1, NRNC
        ROW(I) = EL
2     CONTINUE
C
C  AND
C
      RETURN
C
C  ELSE
C
3     JM1 = J - 1
      RJ = J
      A = SQRT ( RJ * RJ - RJ )
      EL = 1.0E+00 / A

      DO 4 I = 1, JM1
        ROW(I) = EL
4     CONTINUE

      DO 5 I = J, NRNC
        ROW(I) = 0.0E+00
5     CONTINUE

      ROW(J) = ( 1.0E+00 - RJ ) / A

      RETURN
      END
