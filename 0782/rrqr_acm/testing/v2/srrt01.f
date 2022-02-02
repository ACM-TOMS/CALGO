      REAL             FUNCTION SRRT01( TRANS, M, N, A, LDA, JPVT,
     $                 Q, LDQ, R, LDR, WORK, LWORK )
*
*
*
*     .. Scalar Arguments ..
      CHARACTER*1        TRANS
      INTEGER            M, N, LDA, LDQ, LDR, LWORK
*     ..
*     .. Array Arguments ..
      INTEGER            JPVT( * )
      REAL               A( LDA, * ), Q( LDQ, * ), R( LDR, * ),
     $                   WORK( LWORK )
*     ..
*
*  Purpose
*  =======
*
*  SRRT01 tests the Rank-Revealing QR-factorization of matrix A.
*  Array A contains the original matrix being factorized.
*  Argument TRANS says if array Q contains matrix Q or Q' (its
*  transpose).
*  Array R contains matrix R.
*
*  This function returns ||A*P - Q*R||/(||norm(A)||*eps*M)
*
*  Arguments
*  =========
*
*  TRANS   (input) CHARACTER*1
*          If TRANS='N', then array Q contains matrix Q.
*          If TRANS='T', then array Q contains the transpose of
*          matrix Q.
*
*  M       (input) INTEGER
*          Number of rows of matrices A, R and Q,
*          and number of columns of Q.
*
*  N       (input) INTEGER
*          Number of columns of matrices A and R.
*
*  A       (input) REAL array, dimension (LDA, N)
*          Original m by n matrix A.
*
*  JPVT    (input) INTEGER array, dimension (N)
*          Pivot information as returned by SGERQR.
*
*  Q       (input) REAL array, dimension (LDQ,M)
*          Array which contains m by m orthogonal matrix Q
*          or its tranpose, according argument TRANS.
*
*  R       (input) REAL array, dimension (LDR,N)
*          Upper triangular matrix.
*          The lower part of matrix R must contain zeroes.
*
*  LDR     (input) INTEGER
*          Leading dimension of arrays A and R.
*
*  WORK    (workspace) REAL array, dimension (LWORK)
*
*  LWORK   (input) INTEGER
*          Length of array WORK.  LWORK >= M*N.
*
*     .. Parameters ..
      REAL               ZERO, ONE

      PARAMETER          ( ZERO = 0.0E+0, ONE = 1.0E+0 )



*     ..
*     .. Local Scalars ..
      INTEGER            LDWORK, J
      REAL               NORMA
*     ..
*     .. Local Arrays ..
      REAL               RDUMMY( 1 )
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      REAL               SLAMCH, SLANGE
      EXTERNAL           SLAMCH, SLANGE, LSAME
*     ..
*     .. External Subroutines ..
      EXTERNAL           SCOPY, XERBLA, SGEMM
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, REAL
*     ..
*     .. Executable Statements ..
*
      LDWORK = MAX( 1, M )
*
*     Test if there is enough workspace
*
      IF( LWORK.LT.M*N ) THEN
         CALL XERBLA( 'SRRT01', 10 )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( M.LE.0 .OR. N.LE.0 ) THEN
         SRRT01 = ZERO
         RETURN
      END IF
*
      NORMA = SLANGE( 'One-norm', M, N, A, LDA, RDUMMY )
*
*     Compute WORK := A*P.
*
      DO J = 1, N
         CALL SCOPY( M, A( 1, JPVT( J ) ), 1,
     $                WORK( (J-1)*M+1 ), 1 )
      END DO
*
*     Compute WORK := WORK - Q*R.
*
      IF( LSAME( TRANS, 'N' ) ) THEN
         CALL SGEMM( 'No transpose', 'No transpose', M, N, M, -ONE,
     $               Q, LDQ, R, LDR, ONE, WORK, LDWORK )
      ELSE
         CALL SGEMM( 'Transpose', 'No transpose', M, N, M, -ONE,
     $               Q, LDQ, R, LDR, ONE, WORK, LDWORK )
      END IF
*
*     Compute the 1-norm of WORK divided by (max(m,n)*eps).
*
      SRRT01 = SLANGE( 'One-norm', M, N, WORK, LDWORK,
     $                         RDUMMY ) /
     $             ( REAL( MAX( M, N ) )*SLAMCH( 'Epsilon' ) )
      IF( NORMA.NE.ZERO )
     $   SRRT01 = SRRT01 / NORMA
*
      RETURN
*
*     End of SRRT01
*
      END
