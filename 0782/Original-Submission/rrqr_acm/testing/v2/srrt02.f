      REAL             FUNCTION SRRT02( M, QT, LDQT, WORK, LWORK )
*
*
*
*     .. Scalar Arguments ..
      INTEGER            M, LDQT, LWORK
*     ..
*     .. Array Arguments ..
      REAL               QT( LDQT, * ), WORK( LWORK )
*     ..
*
*  Purpose
*  =======
*
*  SRRT02 computes the test ratio
*
*        || Q'*Q - I || / (eps * m)
*
*  where the transpose of matrix Q (Q') is stored in array QT.
*
*
*  Arguments
*  =========
*
*  M       (input) INTEGER
*          Number of rows and columns of matrix QT.
*
*  QT      (input) REAL array, dimension (LDQT,N)
*          Transpose of m by m matrix Q.
*
*  LDQT    (input) INTEGER
*          Leading dimension of array QT.
*
*  WORK    (workspace) REAL array, dimension (LWORK)
*
*  LWORK   (input) INTEGER
*          The length of the array WORK.  LWORK >= M*M.
*
*     .. Parameters ..
      REAL               ZERO, ONE

      PARAMETER          ( ZERO = 0.0E+0, ONE = 1.0E+0 )



*     ..
*     .. Local Scalars ..
      INTEGER            LDWORK
*     ..
*     .. External Functions ..
      REAL               SLAMCH, SLANGE
      EXTERNAL           SLAMCH, SLANGE
*     ..
*     .. External Subroutines ..
      EXTERNAL           SLASET, XERBLA, SGEMM
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, REAL
*     ..
*     .. Local Arrays ..
      REAL               RDUMMY( 1 )
*     ..
*     .. Executable Statements ..
*
      LDWORK = MAX( 1, M )
*
*     Test for sufficient workspace
*
      IF( LWORK.LT.M*M ) THEN
         CALL XERBLA( 'SRRT02', 7 )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( M.LE.0 ) THEN
         SRRT02 = ZERO
         RETURN
      END IF
*
*     Set WORK to the identity.
*
      CALL SLASET( 'All',M, M, ZERO, ONE, WORK, LDWORK )
*
*     Compute WORK := WORK - QT * QT'. That is, WORK := WORK - Q'*Q.
*
      CALL SGEMM( 'No transpose', 'Transpose', M, M, M, -ONE,
     $             QT, LDQT, QT, LDQT, ONE, WORK, LDWORK )
*
*     Compute || WORK || / ( m * eps ).
*
      SRRT02 = SLANGE( 'One-norm', M, M, WORK, LDWORK, RDUMMY )/
     $             ( REAL( M )*SLAMCH( 'Epsilon' ) )
*
      RETURN
*
*     End of SRRT02
*
      END
