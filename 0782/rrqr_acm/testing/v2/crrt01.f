      REAL FUNCTION CRRT01( CNJTRN, M, N, A, LDA, JPVT,
     $                 Q, LDQ, R, LDR, WORK, LWORK )
*
*
*
*     .. Scalar Arguments ..
      CHARACTER*1        CNJTRN
      INTEGER            M, N, LDA, LDQ, LDR, LWORK
*     ..
*     .. Array Arguments ..
      INTEGER            JPVT( * )
      COMPLEX            A( LDA, * ), Q( LDQ, * ), R( LDR, * ),
     $                   WORK( LWORK )
*     ..
*
*  Purpose
*  =======
*
*  CRRT01 tests the Rank-Revealing QR-factorization of matrix A.
*  Array A contains the original matrix being factorized.
*  Argument TRANS says if array Q contains matrix Q or its conjugate
*  transpose.
*  Array R contains matrix R.
*
*  This function returns ||A*P - Q*R||/(||norm(A)||*eps*M)
*
*  Arguments
*  =========
*
*  TRANS   (input) CHARACTER*1
*          If TRANS='N', then array Q contains matrix Q.
*          If TRANS='T', then array Q contains the conjugate transpose
*          of matrix Q.
*
*  M       (input) INTEGER
*          Number of rows of matrices A, R and Q,
*          and number of columns of Q.
*
*  N       (input) INTEGER
*          Number of columns of matrices A and R.
*
*  A       (input) COMPLEX array, dimension (LDA, N)
*          Original m by n matrix A.
*
*  JPVT    (input) INTEGER array, dimension (N)
*          Pivot information as returned by CGERQR.
*
*  Q       (input) COMPLEX array, dimension (LDQ,M)
*          Array which contains m by m orthogonal matrix Q
*          or its tranpose, according argument TRANS.
*
*  R       (input) COMPLEX array, dimension (LDR,N)
*          Upper triangular matrix.
*          The lower part of matrix R must contain zeroes.
*
*  LDR     (input) INTEGER
*          Leading dimension of arrays A and R.
*
*  WORK    (workspace) COMPLEX array, dimension (LWORK)
*
*  LWORK   (input) INTEGER
*          Length of array WORK.  LWORK >= M*N.
*
*     .. Parameters ..
      REAL               ZERO, ONE
      COMPLEX            CZERO, CONE

      PARAMETER          ( ZERO = 0.0E+0, ONE = 1.0E+0,
     $                   CZERO = ( 0.0E+0, 0.0E+0 ),
     $                   CONE = ( 1.0E+0, 0.0E+0 ) )





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
      REAL               SLAMCH, CLANGE
      EXTERNAL           LSAME, SLAMCH, SLANGE
*     ..
*     .. External Subroutines ..
      EXTERNAL           CCOPY, XERBLA, CGEMM
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
         CALL XERBLA( 'CRRT01', 10 )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( M.LE.0 .OR. N.LE.0 ) THEN
         CRRT01 = ZERO
         RETURN
      END IF
*
      NORMA = CLANGE( 'One-norm', M, N, A, LDA, RDUMMY )
*
*     Compute WORK := A*P.
*
      DO J = 1, N
         CALL CCOPY( M, A( 1, JPVT( J ) ), 1,
     $                 WORK( (J-1)*M+1 ), 1 )
      END DO
*
*     Compute WORK := WORK - Q*R.
*
      IF( LSAME( CNJTRN, 'N' ) ) THEN
         CALL CGEMM( 'No transpose', 'No transpose', M, N, M,
     $               -CONE, Q, LDQ, R, LDR, CONE, WORK, LDWORK )
      ELSE
         CALL CGEMM( 'Conjugate transpose', 'No transpose',
     $               M, N, M, -CONE, Q, LDQ, R, LDR, CONE,
     $               WORK, LDWORK )
      END IF
*
*     Compute the 1-norm of WORK divided by (max(m,n)*eps).
*
      CRRT01 = CLANGE( 'One-norm', M, N, WORK, LDWORK,
     $                           RDUMMY ) /
     $              ( REAL( MAX( M, N ) )*SLAMCH( 'Epsilon' ) )
      IF( NORMA.NE.ZERO )
     $   CRRT01 = CRRT01 / NORMA
*
      RETURN
*
*     End of CRRT01
*
      END
