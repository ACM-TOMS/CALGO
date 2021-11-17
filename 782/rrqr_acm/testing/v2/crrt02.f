      REAL FUNCTION CRRT02( M, QT, LDQT, WORK, LWORK )
*
*
*
*     .. Scalar Arguments ..
      INTEGER            M, LDQT, LWORK
*     ..
*     .. Array Arguments ..
      COMPLEX            QT( LDQT, * ), WORK( LWORK )
*     ..
*
*  Purpose
*  =======
*
*  CRRT02 computes the test ratio
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
*  QT      (input) COMPLEX array, dimension (LDQT,N)
*          Transpose of m by m matrix Q.
*
*  LDQT    (input) INTEGER
*          Leading dimension of array QT.
*
*  WORK    (workspace) COMPLEX array, dimension (LWORK)
*
*  LWORK   (input) INTEGER
*          The length of the array WORK.  LWORK >= M*M.
*
*     .. Parameters ..
      REAL               ZERO, ONE
      COMPLEX            CZERO, CONE

      PARAMETER          ( ZERO = 0.0E+0, ONE = 1.0E+0,
     $                   CZERO = ( 0.0E+0, 0.0E+0 ),
     $                   CONE = ( 1.0E+0, 0.0E+0 ) )





*     ..
*     .. Local Scalars ..
      INTEGER            LDWORK
*     ..
*     .. External Functions ..
      REAL               SLAMCH, CLANGE
      EXTERNAL           SLAMCH, CLANGE
*     ..
*     .. External Subroutines ..
      EXTERNAL           CLASET, XERBLA, CGEMM
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
*     Test for sufficient workspace.
*
      IF( LWORK.LT.M*M ) THEN
         CALL XERBLA( 'CRRT02', 7 )
         RETURN
      END IF
*
*     Quick return if possible.
*
      IF( M.LE.0 ) THEN
         CRRT02 = ZERO
         RETURN
      END IF
*
*     Set WORK to the identity.
*
      CALL CLASET( 'All', M, M, CZERO, CONE, WORK, LDWORK )
*
*     Compute WORK := WORK - QT * QT'. That is, WORK := WORK - Q'*Q.
*
      CALL CGEMM( 'No transpose', 'Conjugate transpose',
     $             M, M, M, -CONE, QT, LDQT, QT, LDQT, CONE,
     $             WORK, LDWORK )
*
*     Compute || WORK || / ( m * eps ).
*
      CRRT02 = CLANGE( 'One-norm', M, M, WORK, LDWORK,
     $                           RDUMMY ) /
     $              ( REAL( M )*SLAMCH( 'Epsilon' ) )
*
      RETURN
*
*     End of CRRT02
*
      END
