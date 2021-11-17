C
      SUBROUTINE DTRSYUP( N, A, LDA )
C
      IMPLICIT NONE
C
C     .. Scalar arguments ..
      INTEGER              N, LDA
C     ..
C     .. Array arguments ..
      DOUBLE PRECISION     A( * )
C     ..
C     
C     Purpose and description
C     =======================
C     To simplify the usage of the SYR2K-update in PTRLYCTD we copy
C     the strictly upper triangular part of the NXN matrix A into 
C     the strictly lower triangular part of the matrix A. We call 
C     this a Double precision TRiangular SYmmetric UPdate.
C
C     .. Local Scalars ..
      INTEGER I, J
C     ..
C     .. Executable Statements ..
C
C     Copy now all elements A(I,J) above the main diagonal to the
C     transposed locations A(J,I) below the diagonal.
C
      DO 10 J = 1, N
         DO 20 I = 1, J-1
            A((I-1)*LDA + J) = A((J-1)*LDA + I)
 20      CONTINUE
 10   CONTINUE
C
      END
C
C     End of DTRSYUP
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
