      DOUBLE PRECISION FUNCTION DRRT01( TRANS, M, N, A, LDA, JPVT,
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
      DOUBLE PRECISION   A( LDA, * ), Q( LDQ, * ), R( LDR, * ),
     $                   WORK( LWORK )
*     ..
*
*  Purpose
*  =======
*
*  DRRT01 tests the Rank-Revealing QR-factorization of matrix A.
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
*  A       (input) DOUBLE PRECISION array, dimension (LDA, N)
*          Original m by n matrix A.
*
*  JPVT    (input) INTEGER array, dimension (N)
*          Pivot information as returned by DGERQR.
*
*  Q       (input) DOUBLE PRECISION array, dimension (LDQ,M)
*          Array which contains m by m orthogonal matrix Q
*          or its tranpose, according argument TRANS.
*
*  R       (input) DOUBLE PRECISION array, dimension (LDR,N)
*          Upper triangular matrix.
*          The lower part of matrix R must contain zeroes.
*
*  LDR     (input) INTEGER
*          Leading dimension of arrays A and R.
*
*  WORK    (workspace) DOUBLE PRECISION array, dimension (LWORK)
*
*  LWORK   (input) INTEGER
*          Length of array WORK.  LWORK >= M*N.
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE



      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )

*     ..
*     .. Local Scalars ..
      INTEGER            LDWORK, J
      DOUBLE PRECISION   NORMA
*     ..
*     .. Local Arrays ..
      DOUBLE PRECISION   RDUMMY( 1 )
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      DOUBLE PRECISION   DLAMCH, DLANGE
      EXTERNAL           DLAMCH, DLANGE, LSAME
*     ..
*     .. External Subroutines ..
      EXTERNAL           DCOPY, XERBLA, DGEMM
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, DBLE
*     ..
*     .. Executable Statements ..
*
      LDWORK = MAX( 1, M )
*
*     Test if there is enough workspace
*
      IF( LWORK.LT.M*N ) THEN
         CALL XERBLA( 'DRRT01', 10 )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( M.LE.0 .OR. N.LE.0 ) THEN
         DRRT01 = ZERO
         RETURN
      END IF
*
      NORMA = DLANGE( 'One-norm', M, N, A, LDA, RDUMMY )
*
*     Compute WORK := A*P.
*
      DO J = 1, N
         CALL DCOPY( M, A( 1, JPVT( J ) ), 1,
     $                WORK( (J-1)*M+1 ), 1 )
      END DO
*
*     Compute WORK := WORK - Q*R.
*
      IF( LSAME( TRANS, 'N' ) ) THEN
         CALL DGEMM( 'No transpose', 'No transpose', M, N, M, -ONE,
     $               Q, LDQ, R, LDR, ONE, WORK, LDWORK )
      ELSE
         CALL DGEMM( 'Transpose', 'No transpose', M, N, M, -ONE,
     $               Q, LDQ, R, LDR, ONE, WORK, LDWORK )
      END IF
*
*     Compute the 1-norm of WORK divided by (max(m,n)*eps).
*
      DRRT01 = DLANGE( 'One-norm', M, N, WORK, LDWORK,
     $                         RDUMMY ) /
     $             ( DBLE( MAX( M, N ) )*DLAMCH( 'Epsilon' ) )
      IF( NORMA.NE.ZERO )
     $   DRRT01 = DRRT01 / NORMA
*
      RETURN
*
*     End of DRRT01
*
      END
