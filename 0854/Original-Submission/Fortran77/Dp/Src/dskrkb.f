      SUBROUTINE DSKRKB( UPLO, TRANS, N, K, ALPHA, BETA, A, LDA, B, LDB,
     $                   C, LDC, DWORK, LDWORK, INFO )
      IMPLICIT NONE
C
C     PURPOSE
C
C     To compute the matrix formula
C
C        C := alpha*C + beta*op( A )*B*op( A )',
C
C     where alpha and beta are scalars, B and C are skew-symmetric 
C     matrices, A is a general matrix, and op( A ) is one of
C     
C        op( A ) = A   or   op( A ) = A'.
C     
C     The result is overwritten on C.
C
C     This is a skew-symmetric version of the SLICOT [1] routine MB01RU
C     written by Vasile Sima.
C
C     ARGUMENTS
C
C     Mode Parameters
C
C     UPLO    CHARACTER*1
C             Specifies whether the upper or lower triangular parts of
C             the arrays B and C are to be referenced as follows:
C             = 'U':  only the strictly upper triangular part of B and
C                     C are to be referenced;
C             = 'L':  only the strictly lower triangular part of B and
C                     C are to be referenced.
C            
C     TRANS   CHARACTER*1
C             Specifies the form of op( A ) to be used in the matrix
C             multiplication as follows:
C             = 'N':  op( A ) = A;
C             = 'T':  op( A ) = A';
C             = 'C':  op( A ) = A'.
C
C     Input/Output Parameters
C
C     N       (input) INTEGER
C             The order of the matrix C. N >= 0.
C
C     K       (input) INTEGER
C             The order of the matrix B. K >= 0.
C
C     ALPHA   (input) DOUBLE PRECISION
C             The scalar alpha. When alpha is zero then C need not be  
C             set before entry, except when C is identified with B in
C             the call.
C            
C     BETA    (input) DOUBLE PRECISION
C             The scalar beta. When beta is zero then A and B are not
C             referenced.
C
C     A       (input) DOUBLE PRECISION array, dimension (LDA,m)
C             where m is K when TRANS = 'N' and is N when TRANS = 'T' or
C             TRANS = 'C'.
C             On entry with TRANS = 'N', the leading N-by-K part of this
C             array must contain the matrix A.
C             On entry with TRANS = 'T' or TRANS = 'C', the leading
C             K-by-N part of this array must contain the matrix A.
C
C     LDA     INTEGER
C             The leading dimension of the array A.  LDA >= MAX(1,m),
C             where m is N when TRANS = 'N' and is K when TRANS = 'T' or
C             TRANS = 'C'.
C
C     B       (input) DOUBLE PRECISION array, dimension (LDB,N)
C             On entry, if UPLO = 'U', the leading K-by-K part of this
C             array must contain the strictly upper triangular part of
C             the skew-symmetric matrix B. The lower triangular part
C             of the array is not referenced.
C             On entry, if UPLO = 'L', the leading K-by-K  part of this
C             array must contain the strictly lower triangular part of
C             the symmetric matrix B. The upper triangular part of the
C             array is not referenced.
C
C     LDB     INTEGER
C             The leading dimension of the array B. LDX >= MAX(1,K).
C
C     C       (input/output)  DOUBLE PRECISION array, dimension (LDA,N)
C             On entry with UPLO = 'U', the leading N-by-N part of this
C             array must contain the strictly upper triangular part of
C             the matrix C. The lower triangular part of this array is
C             not referenced.
C             On entry with UPLO = 'L', the leading N-by-N part of this
C             array must contain the strictly lower triangular part of
C             the matrix C. The upper triangular part of this array is
C             not referenced.
C             On exit with UPLO = 'U', the leading N-by-N part of this
C             array contains the strictly upper triangular part of the
C             updated matrix C.
C             On exit with UPLO = 'L', the leading N-by-N part of this
C             array contains the strictly lower triangular part of the
C             updated matrix C.
C
C     LDC     INTEGER
C             The leading dimension of the array C.  LDC >= MAX(1,N)
C
C     Workspace
C
C     DWORK   DOUBLE PRECISION array, dimension (LDWORK,m)
C             where m is K-1 when TRANS = 'N' and is N when TRANS = 'T'
C             or TRANS = 'C'.
C             This array is not referenced when beta = 0.
C
C     LDWORK  INTEGER
C             The leading dimension of the array DWORK.
C             If beta = 0, LDWORK >= MAX(1,m), where m is N when
C             TRANS = 'N' and is K-1 when TRANS = 'T' or TRANS = 'C'.
C             Otherwise, LDWORK >= 1.
C
C     Error Indicator
C
C     INFO    INTEGER
C             = 0:  successful exit;
C             < 0:  if INFO = -i, the i-th argument had an illegal
C                   value.
C
C     METHOD
C
C     The matrix expression is efficiently evaluated taking the
C     skew-symmetry into account. Specifically, let B = T - T', with T
C     an upper or lower triangular matrix, defined by 
C
C        T = triu( B ),  if UPLO = 'U',
C        T = tril( B ),  if UPLO = 'L',
C     
C     where triu and tril denote the upper triangular part and lower
C     triangular part of X, respectively. Then,
C
C        A*B*A' = ( A*T )*A' - A*( A*T )',  for TRANS = 'N',
C        A'*B*A = A'*( T*A ) - ( T*A )'*A,  for TRANS = 'T', or 'C',
C     
C     involving routines DSKR2K and DTRMM.
C
C     CONTRIBUTORS
C
C     D. Kressner (Technical Univ. Berlin, Germany) and
C     P. Benner (Technical Univ. Chemnitz, Germany), December 2003.
C
C     KEYWORDS
C
C     Elementary matrix operations, skew-symmetric matrix.
C
C     ******************************************************************
C
C     .. Parameters ..
      DOUBLE PRECISION  ZERO, ONE
      PARAMETER         ( ZERO = 0.0D0, ONE = 1.0D0 )
C     .. Scalar Arguments ..
      CHARACTER         TRANS, UPLO
      INTEGER           INFO, I, J, K, LDA, LDB, LDC, LDWORK, N
      DOUBLE PRECISION  ALPHA, BETA
C     .. Array Arguments ..
      DOUBLE PRECISION  A(LDA,*), DWORK(LDWORK,*), B(LDB,*), C(LDC,*)
C     .. Local Scalars ..
      LOGICAL           LTRANS, LUPLO
      INTEGER           IERR
C     .. External Functions ..
      LOGICAL           LSAME
      EXTERNAL          LSAME
C     .. External Subroutines ..
      EXTERNAL          DLACPY, DLASCL, DLASET, DSKR2K, DTRMM, XERBLA
C     .. Intrinsic Functions ..
      INTRINSIC         MAX
C     .. Executable Statements ..
C
C     Test the input scalar arguments.
C
      INFO = 0
      LUPLO  = LSAME( UPLO,  'U' )
      LTRANS = LSAME( TRANS, 'T' ) .OR. LSAME( TRANS, 'C' )
C
      IF ( ( .NOT.LUPLO  ).AND.( .NOT.LSAME( UPLO,  'L' ) ) ) THEN
         INFO = -1
      ELSE IF ( ( .NOT.LTRANS ).AND.( .NOT.LSAME( TRANS, 'N' ) ) )THEN
         INFO = -2
      ELSE IF ( N.LT.0 ) THEN
         INFO = -3
      ELSE IF ( K.LT.0 ) THEN
         INFO = -4
      ELSE IF ( LDA.LT.1 .OR. ( LTRANS .AND. LDA.LT.K ) .OR.
     $                  ( .NOT.LTRANS .AND. LDA.LT.N ) ) THEN
         INFO = -8
      ELSE IF ( LDB.LT.MAX( 1, K ) ) THEN
         INFO = -10
      ELSE IF ( LDC.LT.MAX( 1, N ) ) THEN
         INFO = -12
      ELSE IF ( LDWORK.LT.1 .OR. ( BETA.NE.ZERO .AND.
     $          ( ( LTRANS .AND. LDWORK.LT.K-1 ) .OR.
     $            ( .NOT.LTRANS .AND. LDWORK.LT.N ) ) ) ) THEN
         INFO = -14
      END IF
C
      IF ( INFO.NE.0 ) THEN
C
C        Error return.
C
         CALL XERBLA( 'DSKRKB', -INFO )
         RETURN
      END IF
C
C     Quick return if possible.
C
      IF ( N.LE.1 )
     $   RETURN
C
      IF ( LUPLO ) THEN
         I = 1
         J = 2
      ELSE
         I = 2
         J = 1
      END IF
C
      IF ( BETA.EQ.ZERO ) THEN
         IF ( ALPHA.EQ.ZERO ) THEN
C
C           Special case when both alpha = 0 and beta = 0.
C
            CALL DLASET( UPLO, N-1, N-1, ZERO, ZERO, C(I,J), LDC )
         ELSE
C
C           Special case beta = 0.
C
            IF ( ALPHA.NE.ONE ) 
     $         CALL DLASCL( UPLO, 0, 0, ONE, ALPHA, N-1, N-1, C(I,J),
     $                      LDC, INFO )
         END IF
         RETURN
      END IF
C
      IF ( K.LE.1 ) 
     $   RETURN
C     
C     General case: beta <> 0.
C     Compute W = op( A )*T or W = T*op( A ) in DWORK, and apply the
C     updating formula (see METHOD section).
C
      IF ( LTRANS ) THEN
         CALL DLACPY( 'Full', K-1, N, A(J,1), LDA, DWORK, LDWORK )
         CALL DTRMM( 'Left', UPLO, 'NoTranspose', 'Non-unit', K-1, N,
     $               ONE, B(I,J), LDB, DWORK, LDWORK )
         CALL DSKR2K( UPLO, TRANS, N, K-1, -BETA, DWORK, LDWORK, A(I,1),
     $                LDA, ALPHA, C, LDC, IERR )
      ELSE
         CALL DLACPY( 'Full', N, K-1, A(1,I), LDA, DWORK, LDWORK )
         CALL DTRMM(  'Right', UPLO, 'NoTranspose', 'Non-unit', N, K-1,
     $                ONE, B(I,J), LDB, DWORK, LDWORK )
         CALL DSKR2K( UPLO, TRANS, N, K-1, -BETA, DWORK, LDWORK, A(1,J),
     $                LDA, ALPHA, C, LDC, IERR )
      END IF
C
      RETURN
C *** Last line of DSKRKB ***
      END
