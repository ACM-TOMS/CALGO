      DOUBLE PRECISION FUNCTION ZRRT02( M, QT, LDQT, WORK, LWORK )
*
*
*
*     .. Scalar Arguments ..
      INTEGER            M, LDQT, LWORK
*     ..
*     .. Array Arguments ..
      COMPLEX*16         QT( LDQT, * ), WORK( LWORK )
*     ..
*
*  Purpose
*  =======
*
*  ZRRT02 computes the test ratio
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
*  QT      (input) COMPLEX*16 array, dimension (LDQT,N)
*          Transpose of m by m matrix Q.
*
*  LDQT    (input) INTEGER
*          Leading dimension of array QT.
*
*  WORK    (workspace) COMPLEX*16 array, dimension (LWORK)
*
*  LWORK   (input) INTEGER
*          The length of the array WORK.  LWORK >= M*M.
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      COMPLEX*16         CZERO, CONE





      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0,
     $                   CZERO = ( 0.0D+0, 0.0D+0 ),
     $                   CONE = ( 1.0D+0, 0.0D+0 ) )

*     ..
*     .. Local Scalars ..
      INTEGER            LDWORK
*     ..
*     .. External Functions ..
      DOUBLE PRECISION   DLAMCH, ZLANGE
      EXTERNAL           DLAMCH, ZLANGE
*     ..
*     .. External Subroutines ..
      EXTERNAL           ZLASET, XERBLA, ZGEMM
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
*     Test for sufficient workspace.
*
      IF( LWORK.LT.M*M ) THEN
         CALL XERBLA( 'ZRRT02', 7 )
         RETURN
      END IF
*
*     Quick return if possible.
*
      IF( M.LE.0 ) THEN
         ZRRT02 = ZERO
         RETURN
      END IF
*
*     Set WORK to the identity.
*
      CALL ZLASET( 'All', M, M, CZERO, CONE, WORK, LDWORK )
*
*     Compute WORK := WORK - QT * QT'. That is, WORK := WORK - Q'*Q.
*
      CALL ZGEMM( 'No transpose', 'Conjugate transpose',
     $             M, M, M, -CONE, QT, LDQT, QT, LDQT, CONE,
     $             WORK, LDWORK )
*
*     Compute || WORK || / ( m * eps ).
*
      ZRRT02 = ZLANGE( 'One-norm', M, M, WORK, LDWORK,
     $                           RDUMMY ) /
     $              ( DBLE( M )*DLAMCH( 'Epsilon' ) )
*
      RETURN
*
*     End of ZRRT02
*
      END
