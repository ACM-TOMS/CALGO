      SUBROUTINE DSKUPD( UPLO, TRANS, M, N, A, LDA, Z, LDZ, DWORK,
     $                   INFO )
      IMPLICIT NONE
C
C     PURPOSE
C
C     To compute the transformation of the symmetric matrix A by the
C     matrix Z in the form
C
C        A := op(Z)*A*op(Z)',
C
C     where op(Z) is either Z or its transpose, Z'.
C
C     This is a skew-symmetric version of the SLICOT[1] routine MB01RW
C     written by Andras Varga, modified by Thilo Penzl and Vasile
C     Sima. It is also a simpler version of DSKR2B. The routine
C     DSKR2B requires more workspace but less floating point operations.
C
C     ARGUMENTS
C
C     Mode Parameters
C
C     UPLO    CHARACTER*1
C             Specifies whether the upper or lower triangular parts of
C             the array A is to be referenced as follows:
C             = 'U':  only the strictly upper triangular part of A is
C                     to be referenced;
C             = 'L':  only the strictly lower triangular part of A is
C                     to be referenced.
C            
C     TRANS   CHARACTER*1
C             Specifies the form of op( Z ) to be used in the matrix
C             multiplication as follows:
C             = 'N':  op( Z ) = Z;
C             = 'T':  op( Z ) = Z';
C             = 'C':  op( Z ) = Z'.
C
C     Input/Output Parameters
C
C     M       (input) INTEGER
C             The order of the resulting skew-symmetric matrix
C             op(Z)*A*op(Z)'.  M >= 0.
C
C     N       (input) INTEGER
C             The order of the matrix A.  N >= 0.
C
C     A       (input/output)  DOUBLE PRECISION array, dimension
C             ( LDA,MAX( M, N ) )
C             On entry with UPLO = 'U', the leading N-by-N part of this
C             array must contain the strictly upper triangular part of
C             the matrix A.
C             On entry with UPLO = 'L', the leading N-by-N part of this
C             array must contain the strictly lower triangular part of
C             the matrix A.
C             On exit with UPLO = 'U', the leading M-by-M part of this
C             array contains the strictly upper triangular part of the
C             updated matrix A.
C             On exit with UPLO = 'L', the leading M-by-M part of this
C             array contains the strictly lower triangular part of the
C             updated matrix A.
C
C     LDA     INTEGER
C             The leading dimension of the array A.  LDA >= MAX(1,M,N)
C
C     Z       (input) DOUBLE PRECISION array, dimension (LDQ,K)
C             where K = N if TRANS = 'N' and K = M if TRANS = 'T' or
C             TRANS = 'C'.
C             On entry with TRANS = 'N', the leading M-by-N part of this
C             array must contain the matrix Z.
C             On entry with TRANS = 'T' or TRANS = 'C' the leading
C             N-by-M part of this array must contain the matrix Z.
C
C     LDZ     INTEGER
C             The leading dimension of the array Z.
C             LDZ >= MAX(1,M) if TRANS = 'N' and
C             LDZ >= MAX(1,N) if TRANS = 'T' or TRANS = 'C'.
C
C     Workspace
C
C     DWORK   DOUBLE PRECISION array, dimension (N)
C
C     Error Indicator
C
C     INFO    INTEGER
C             = 0:  successful exit;
C             < 0:  if INFO = -i, the i-th argument had an illegal
C                   value.
C
C     REFERENCES
C
C     [1] P. BENNER, V. MEHRMANN, V. SIMA, S. VAN HUFFEL, and A. VARGA:
C     SLICOT - a subroutine library in systems and control theory. In:
C     Applied and computational control, signals, and circuits, Vol. 1,
C     499--539, Birkhaeuser Boston, 1999.
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
      CHARACTER*1       TRANS, UPLO
      INTEGER           INFO, LDA, LDZ, M, N
C     .. Array Arguments ..
      DOUBLE PRECISION  A(LDA,*), DWORK(*), Z(LDZ,*)
C     .. Local Scalars ..
      LOGICAL           NOTTRA, UPPER
      INTEGER           I, J
C     .. External Functions ..
      LOGICAL           LSAME
      EXTERNAL          LSAME
C     .. External Subroutines ..
      EXTERNAL          DCOPY, DGEMV, DSCAL, XERBLA
C     .. Intrinsic Functions ..
      INTRINSIC         MAX
C
C     .. Executable Statements
C
      NOTTRA = LSAME( TRANS, 'N' )
      UPPER  = LSAME( UPLO,  'U' )
C
      INFO = 0
      IF( .NOT.( UPPER .OR. LSAME( UPLO, 'L') ) ) THEN
         INFO = -1
      ELSE IF ( .NOT.( NOTTRA .OR. LSAME( TRANS, 'T')
     $         .OR. LSAME( TRANS, 'C') ) ) THEN
         INFO = -2
      ELSE IF ( M.LT.0 ) THEN
         INFO = -3
      ELSE IF ( N.LT.0 ) THEN
         INFO = -4
      ELSE IF ( LDA.LT.MAX( 1, M, N ) ) THEN
         INFO = -6
      ELSE IF ( (      NOTTRA .AND. LDZ.LT.MAX( 1, M ) )   .OR.
     $          ( .NOT.NOTTRA .AND. LDZ.LT.MAX( 1, N ) ) ) THEN
         INFO = -8
      END IF
C
      IF ( INFO.NE.0 ) THEN
         CALL XERBLA( 'DSKUPD', -INFO )
         RETURN
      END IF
C
C     Quick return if possible.
C
      IF ( N.EQ.0 .OR. M.EQ.0 ) 
     $   RETURN
C
      IF ( NOTTRA ) THEN
C
C        Compute Z*A*Z'.
C
         IF ( UPPER ) THEN
C
C           Compute Z*A in A (M-by-N).
C
            DO 10 J = 1, N
               CALL DCOPY( J-1, A(1,J), 1, DWORK, 1 )
               DWORK(J) = ZERO
               CALL DCOPY( N-J, A(J,J+1), LDA, DWORK(J+1), 1 )
                  CALL DSCAL( N-J, -ONE, DWORK(J+1), 1 )
               CALL DGEMV( TRANS, M, N, ONE, Z, LDZ, DWORK, 1, ZERO,
     $                     A(1,J), 1 )
   10       CONTINUE
C
C           Compute A*Z' in the upper triangular part of A.
C
            DO 20 I = 1, M
               CALL DCOPY( N, A(I,1), LDA, DWORK, 1 )
               CALL DGEMV( TRANS, M-I, N, ONE, Z(I+1,1), LDZ, DWORK, 1,
     $                     ZERO, A(I,I+1), LDA )
   20       CONTINUE
C
         ELSE
C
C           Compute A*Z' in A (N-by-M).
C
            DO 30 I = 1, N
               CALL DCOPY( I-1, A(I,1), LDA, DWORK, 1 )
               DWORK(I) = ZERO
               CALL DCOPY( N-I, A(I+1,I), 1, DWORK(I+1), 1 )
               CALL DSCAL( N-I, -ONE, DWORK(I+1), 1 )
               CALL DGEMV( TRANS, M, N, ONE, Z, LDZ, DWORK, 1, ZERO,
     $                     A(I,1), LDA )
   30       CONTINUE
C
C           Compute Z*A in the lower triangular part of A.
C
            DO 40 J = 1, M
               CALL DCOPY( N, A(1,J), 1, DWORK, 1 )
               CALL DGEMV( TRANS, M-J, N, ONE, Z(J+1,1), LDZ, DWORK, 1, 
     $                     ZERO, A(J+1,J), 1 )
   40       CONTINUE
C
         END IF
      ELSE
C
C        Compute Z'*A*Z.
C
         IF ( UPPER ) THEN
C
C           Compute Z'*A in A (M-by-N).
C
            DO 50 J = 1, N
               CALL DCOPY( J-1, A(1,J), 1, DWORK, 1 )
               DWORK(J) = ZERO
               CALL DCOPY( N-J, A(J,J+1), LDA, DWORK(J+1), 1 )
               CALL DSCAL( N-J, -ONE, DWORK(J+1), 1 )
               CALL DGEMV( TRANS, N, M, ONE, Z, LDZ, DWORK, 1, ZERO,
     $                     A(1,J), 1 )
   50       CONTINUE
C
C           Compute A*Z in the upper triangular part of A.
C
            DO 60 I = 1, M
               CALL DCOPY( N, A(I,1), LDA, DWORK, 1 )
               CALL DGEMV( TRANS, N, M-I, ONE, Z(1,I+1), LDZ, DWORK, 1, 
     $                     ZERO, A(I,I+1), LDA )
   60       CONTINUE
C
         ELSE
C
C           Compute A*Z in A (N-by-M).
C
            DO 70 I = 1, N
               CALL DCOPY( I-1, A(I,1), LDA, DWORK, 1 )
               DWORK(I) = ZERO
               CALL DCOPY( N-I, A(I+1,I), 1, DWORK(I+1), 1 )
               CALL DSCAL( N-I, -ONE, DWORK(I+1), 1 )
               CALL DGEMV( TRANS, N, M, ONE, Z, LDZ, DWORK, 1, ZERO,
     $                     A(I,1), LDA )
   70       CONTINUE
C
C           Compute Z'*A in the lower triangular part of A.
C
            DO 80 J = 1, M
               CALL DCOPY( N, A(1,J), 1, DWORK, 1 )
               CALL DGEMV( TRANS, N, M-J, ONE, Z(1,J+1), LDZ, DWORK, 1, 
     $                     ZERO, A(J+1,J), 1 )
   80       CONTINUE
C
         END IF
      END IF
C     
      RETURN
C *** Last line of DSKUPD ***
      END
