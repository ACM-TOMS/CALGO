      DOUBLE PRECISION FUNCTION DRRT02( M, QT, LDQT, WORK, LWORK )
*
*
*
*     .. Scalar Arguments ..
      INTEGER            M, LDQT, LWORK
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   QT( LDQT, * ), WORK( LWORK )
*     ..
*
*  Purpose
*  =======
*
*  DRRT02 computes the test ratio
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
*  QT      (input) DOUBLE PRECISION array, dimension (LDQT,N)
*          Transpose of m by m matrix Q.
*
*  LDQT    (input) INTEGER
*          Leading dimension of array QT.
*
*  WORK    (workspace) DOUBLE PRECISION array, dimension (LWORK)
*
*  LWORK   (input) INTEGER
*          The length of the array WORK.  LWORK >= M*M.
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE



      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )

*     ..
*     .. Local Scalars ..
      INTEGER            LDWORK
*     ..
*     .. External Functions ..
      DOUBLE PRECISION   DLAMCH, DLANGE
      EXTERNAL           DLAMCH, DLANGE
*     ..
*     .. External Subroutines ..
      EXTERNAL           DLASET, XERBLA, DGEMM
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, DBLE
*     ..
*     .. Local Arrays ..
      DOUBLE PRECISION   RDUMMY( 1 )
*     ..
*     .. Executable Statements ..
*
      LDWORK = MAX( 1, M )
*
*     Test for sufficient workspace
*
      IF( LWORK.LT.M*M ) THEN
         CALL XERBLA( 'DRRT02', 7 )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( M.LE.0 ) THEN
         DRRT02 = ZERO
         RETURN
      END IF
*
*     Set WORK to the identity.
*
      CALL DLASET( 'All',M, M, ZERO, ONE, WORK, LDWORK )
*
*     Compute WORK := WORK - QT * QT'. That is, WORK := WORK - Q'*Q.
*
      CALL DGEMM( 'No transpose', 'Transpose', M, M, M, -ONE,
     $             QT, LDQT, QT, LDQT, ONE, WORK, LDWORK )
*
*     Compute || WORK || / ( m * eps ).
*
      DRRT02 = DLANGE( 'One-norm', M, M, WORK, LDWORK, RDUMMY )/
     $             ( DBLE( M )*DLAMCH( 'Epsilon' ) )
*
      RETURN
*
*     End of DRRT02
*
      END
