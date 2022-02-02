      SUBROUTINE DPBSV( UPLO, N, KD, NRHS, AB, LDAB, B, LDB, INFO )
*
*  -- LAPACK driver routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     March 31, 1993
*
*     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            INFO, KD, LDAB, LDB, N, NRHS
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   AB( LDAB, * ), B( LDB, * )
*     ..
*
*  Purpose
*  =======
*
*  DPBSV computes the solution to a real system of linear equations
*     A * X = B,
*  where A is an N-by-N symmetric positive definite band matrix and X
*  and B are N-by-NRHS matrices.
*
*  The Cholesky decomposition is used to factor A as
*     A = U**T * U,  if UPLO = 'U', or
*     A = L * L**T,  if UPLO = 'L',
*  where U is an upper triangular band matrix, and L is a lower
*  triangular band matrix, with the same number of superdiagonals or
*  subdiagonals as A.  The factored form of A is then used to solve the
*  system of equations A * X = B.
*
*  Arguments
*  =========
*
*  UPLO    (input) CHARACTER*1
*          = 'U':  Upper triangle of A is stored;
*          = 'L':  Lower triangle of A is stored.
*
*  N       (input) INTEGER
*          The number of linear equations, i.e., the order of the
*          matrix A.  N >= 0.
*
*  KD      (input) INTEGER
*          The number of superdiagonals of the matrix A if UPLO = 'U',
*          or the number of subdiagonals if UPLO = 'L'.  KD >= 0.
*
*  NRHS    (input) INTEGER
*          The number of right hand sides, i.e., the number of columns
*          of the matrix B.  NRHS >= 0.
*
*  AB      (input/output) DOUBLE PRECISION array, dimension (LDAB,N)
*          On entry, the upper or lower triangle of the symmetric band
*          matrix A, stored in the first KD+1 rows of the array.  The
*          j-th column of A is stored in the j-th column of the array AB
*          as follows:
*          if UPLO = 'U', AB(KD+1+i-j,j) = A(i,j) for max(1,j-KD)<=i<=j;
*          if UPLO = 'L', AB(1+i-j,j)    = A(i,j) for j<=i<=min(N,j+KD).
*          See below for further details.
*
*          On exit, if INFO = 0, the triangular factor U or L from the
*          Cholesky factorization A = U**T*U or A = L*L**T of the band
*          matrix A, in the same storage format as A.
*
*  LDAB    (input) INTEGER
*          The leading dimension of the array AB.  LDAB >= KD+1.
*
*  B       (input/output) DOUBLE PRECISION array, dimension (LDB,NRHS)
*          On entry, the N-by-NRHS right hand side matrix B.
*          On exit, if INFO = 0, the N-by-NRHS solution matrix X.
*
*  LDB     (input) INTEGER
*          The leading dimension of the array B.  LDB >= max(1,N).
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value
*          > 0:  if INFO = i, the leading minor of order i of A is not
*                positive definite, so the factorization could not be
*                completed, and the solution has not been computed.
*
*  Further Details
*  ===============
*
*  The band storage scheme is illustrated by the following example, when
*  N = 6, KD = 2, and UPLO = 'U':
*
*  On entry:                       On exit:
*
*      *    *   a13  a24  a35  a46      *    *   u13  u24  u35  u46
*      *   a12  a23  a34  a45  a56      *   u12  u23  u34  u45  u56
*     a11  a22  a33  a44  a55  a66     u11  u22  u33  u44  u55  u66
*
*  Similarly, if UPLO = 'L' the format of A is as follows:
*
*  On entry:                       On exit:
*
*     a11  a22  a33  a44  a55  a66     l11  l22  l33  l44  l55  l66
*     a21  a32  a43  a54  a65   *      l21  l32  l43  l54  l65   *
*     a31  a42  a53  a64   *    *      l31  l42  l53  l64   *    *
*
*  Array elements marked * are not used by the routine.
*
*  =====================================================================
*
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. External Subroutines ..
      EXTERNAL           DPBTRF, DPBTRS, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
      IF( .NOT.LSAME( UPLO, 'U' ) .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( KD.LT.0 ) THEN
         INFO = -3
      ELSE IF( NRHS.LT.0 ) THEN
         INFO = -4
      ELSE IF( LDAB.LT.KD+1 ) THEN
         INFO = -6
      ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
         INFO = -8
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DPBSV ', -INFO )
         RETURN
      END IF
*
*     Compute the Cholesky factorization A = U'*U or A = L*L'.
*
      CALL DPBTRF( UPLO, N, KD, AB, LDAB, INFO )
      IF( INFO.EQ.0 ) THEN
*
*        Solve the system A*X = B, overwriting B with X.
*
         CALL DPBTRS( UPLO, N, KD, NRHS, AB, LDAB, B, LDB, INFO )
*
      END IF
      RETURN
*
*     End of DPBSV
*
      END
      SUBROUTINE DPBTRF( UPLO, N, KD, AB, LDAB, INFO )
*
*  -- LAPACK routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     March 31, 1993
*
*     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            INFO, KD, LDAB, N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   AB( LDAB, * )
*     ..
*
*  Purpose
*  =======
*
*  DPBTRF computes the Cholesky factorization of a real symmetric
*  positive definite band matrix A.
*
*  The factorization has the form
*     A = U**T * U,  if UPLO = 'U', or
*     A = L  * L**T,  if UPLO = 'L',
*  where U is an upper triangular matrix and L is lower triangular.
*
*  Arguments
*  =========
*
*  UPLO    (input) CHARACTER*1
*          = 'U':  Upper triangle of A is stored;
*          = 'L':  Lower triangle of A is stored.
*
*  N       (input) INTEGER
*          The order of the matrix A.  N >= 0.
*
*  KD      (input) INTEGER
*          The number of superdiagonals of the matrix A if UPLO = 'U',
*          or the number of subdiagonals if UPLO = 'L'.  KD >= 0.
*
*  AB      (input/output) DOUBLE PRECISION array, dimension (LDAB,N)
*          On entry, the upper or lower triangle of the symmetric band
*          matrix A, stored in the first KD+1 rows of the array.  The
*          j-th column of A is stored in the j-th column of the array AB
*          as follows:
*          if UPLO = 'U', AB(kd+1+i-j,j) = A(i,j) for max(1,j-kd)<=i<=j;
*          if UPLO = 'L', AB(1+i-j,j)    = A(i,j) for j<=i<=min(n,j+kd).
*
*          On exit, if INFO = 0, the triangular factor U or L from the
*          Cholesky factorization A = U**T*U or A = L*L**T of the band
*          matrix A, in the same storage format as A.
*
*  LDAB    (input) INTEGER
*          The leading dimension of the array AB.  LDAB >= KD+1.
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value
*          > 0:  if INFO = i, the leading minor of order i is not
*                positive definite, and the factorization could not be
*                completed.
*
*  Further Details
*  ===============
*
*  The band storage scheme is illustrated by the following example, when
*  N = 6, KD = 2, and UPLO = 'U':
*
*  On entry:                       On exit:
*
*      *    *   a13  a24  a35  a46      *    *   u13  u24  u35  u46
*      *   a12  a23  a34  a45  a56      *   u12  u23  u34  u45  u56
*     a11  a22  a33  a44  a55  a66     u11  u22  u33  u44  u55  u66
*
*  Similarly, if UPLO = 'L' the format of A is as follows:
*
*  On entry:                       On exit:
*
*     a11  a22  a33  a44  a55  a66     l11  l22  l33  l44  l55  l66
*     a21  a32  a43  a54  a65   *      l21  l32  l43  l54  l65   *
*     a31  a42  a53  a64   *    *      l31  l42  l53  l64   *    *
*
*  Array elements marked * are not used by the routine.
*
*  Contributed by
*  Peter Mayes and Giuseppe Radicati, IBM ECSEC, Rome, March 23, 1989
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
      INTEGER            NBMAX, LDWORK
      PARAMETER          ( NBMAX = 32, LDWORK = NBMAX+1 )
*     ..
*     .. Local Scalars ..
      INTEGER            I, I2, I3, IB, II, J, JJ, NB
*     ..
*     .. Local Arrays ..
      DOUBLE PRECISION   WORK( LDWORK, NBMAX )
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            ILAENV
      EXTERNAL           LSAME, ILAENV
*     ..
*     .. External Subroutines ..
      EXTERNAL           DGEMM, DPBTF2, DPOTF2, DSYRK, DTRSM, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MIN
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
      IF( ( .NOT.LSAME( UPLO, 'U' ) ) .AND.
     $    ( .NOT.LSAME( UPLO, 'L' ) ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( KD.LT.0 ) THEN
         INFO = -3
      ELSE IF( LDAB.LT.KD+1 ) THEN
         INFO = -5
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DPBTRF', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 )
     $   RETURN
*
*     Determine the block size for this environment
*
      NB = ILAENV( 1, 'DPBTRF', UPLO, N, KD, -1, -1 )
*
*     The block size must not exceed the semi-bandwidth KD, and must not
*     exceed the limit set by the size of the local array WORK.
*
      NB = MIN( NB, NBMAX )
*
      IF( NB.LE.1 .OR. NB.GT.KD ) THEN
*
*        Use unblocked code
*
         CALL DPBTF2( UPLO, N, KD, AB, LDAB, INFO )
      ELSE
*
*        Use blocked code
*
         IF( LSAME( UPLO, 'U' ) ) THEN
*
*           Compute the Cholesky factorization of a symmetric band
*           matrix, given the upper triangle of the matrix in band
*           storage.
*
*           Zero the upper triangle of the work array.
*
            DO 20 J = 1, NB
               DO 10 I = 1, J - 1
                  WORK( I, J ) = ZERO
   10          CONTINUE
   20       CONTINUE
*
*           Process the band matrix one diagonal block at a time.
*
            DO 70 I = 1, N, NB
               IB = MIN( NB, N-I+1 )
*
*              Factorize the diagonal block
*
               CALL DPOTF2( UPLO, IB, AB( KD+1, I ), LDAB-1, II )
               IF( II.NE.0 ) THEN
                  INFO = I + II - 1
                  GO TO 150
               END IF
               IF( I+IB.LE.N ) THEN
*
*                 Update the relevant part of the trailing submatrix.
*                 If A11 denotes the diagonal block which has just been
*                 factorized, then we need to update the remaining
*                 blocks in the diagram:
*
*                    A11   A12   A13
*                          A22   A23
*                                A33
*
*                 The numbers of rows and columns in the partitioning
*                 are IB, I2, I3 respectively. The blocks A12, A22 and
*                 A23 are empty if IB = KD. The upper triangle of A13
*                 lies outside the band.
*
                  I2 = MIN( KD-IB, N-I-IB+1 )
                  I3 = MIN( IB, N-I-KD+1 )
*
                  IF( I2.GT.0 ) THEN
*
*                    Update A12
*
                     CALL DTRSM( 'Left', 'Upper', 'Transpose',
     $                           'Non-unit', IB, I2, ONE, AB( KD+1, I ),
     $                           LDAB-1, AB( KD+1-IB, I+IB ), LDAB-1 )
*
*                    Update A22
*
                     CALL DSYRK( 'Upper', 'Transpose', I2, IB, -ONE,
     $                           AB( KD+1-IB, I+IB ), LDAB-1, ONE,
     $                           AB( KD+1, I+IB ), LDAB-1 )
                  END IF
*
                  IF( I3.GT.0 ) THEN
*
*                    Copy the lower triangle of A13 into the work array.
*
                     DO 40 JJ = 1, I3
                        DO 30 II = JJ, IB
                           WORK( II, JJ ) = AB( II-JJ+1, JJ+I+KD-1 )
   30                   CONTINUE
   40                CONTINUE
*
*                    Update A13 (in the work array).
*
                     CALL DTRSM( 'Left', 'Upper', 'Transpose',
     $                           'Non-unit', IB, I3, ONE, AB( KD+1, I ),
     $                           LDAB-1, WORK, LDWORK )
*
*                    Update A23
*
                     IF( I2.GT.0 )
     $                  CALL DGEMM( 'Transpose', 'No Transpose', I2, I3,
     $                              IB, -ONE, AB( KD+1-IB, I+IB ),
     $                              LDAB-1, WORK, LDWORK, ONE,
     $                              AB( 1+IB, I+KD ), LDAB-1 )
*
*                    Update A33
*
                     CALL DSYRK( 'Upper', 'Transpose', I3, IB, -ONE,
     $                           WORK, LDWORK, ONE, AB( KD+1, I+KD ),
     $                           LDAB-1 )
*
*                    Copy the lower triangle of A13 back into place.
*
                     DO 60 JJ = 1, I3
                        DO 50 II = JJ, IB
                           AB( II-JJ+1, JJ+I+KD-1 ) = WORK( II, JJ )
   50                   CONTINUE
   60                CONTINUE
                  END IF
               END IF
   70       CONTINUE
         ELSE
*
*           Compute the Cholesky factorization of a symmetric band
*           matrix, given the lower triangle of the matrix in band
*           storage.
*
*           Zero the lower triangle of the work array.
*
            DO 90 J = 1, NB
               DO 80 I = J + 1, NB
                  WORK( I, J ) = ZERO
   80          CONTINUE
   90       CONTINUE
*
*           Process the band matrix one diagonal block at a time.
*
            DO 140 I = 1, N, NB
               IB = MIN( NB, N-I+1 )
*
*              Factorize the diagonal block
*
               CALL DPOTF2( UPLO, IB, AB( 1, I ), LDAB-1, II )
               IF( II.NE.0 ) THEN
                  INFO = I + II - 1
                  GO TO 150
               END IF
               IF( I+IB.LE.N ) THEN
*
*                 Update the relevant part of the trailing submatrix.
*                 If A11 denotes the diagonal block which has just been
*                 factorized, then we need to update the remaining
*                 blocks in the diagram:
*
*                    A11
*                    A21   A22
*                    A31   A32   A33
*
*                 The numbers of rows and columns in the partitioning
*                 are IB, I2, I3 respectively. The blocks A21, A22 and
*                 A32 are empty if IB = KD. The lower triangle of A31
*                 lies outside the band.
*
                  I2 = MIN( KD-IB, N-I-IB+1 )
                  I3 = MIN( IB, N-I-KD+1 )
*
                  IF( I2.GT.0 ) THEN
*
*                    Update A21
*
                     CALL DTRSM( 'Right', 'Lower', 'Transpose',
     $                           'Non-unit', I2, IB, ONE, AB( 1, I ),
     $                           LDAB-1, AB( 1+IB, I ), LDAB-1 )
*
*                    Update A22
*
                     CALL DSYRK( 'Lower', 'No Transpose', I2, IB, -ONE,
     $                           AB( 1+IB, I ), LDAB-1, ONE,
     $                           AB( 1, I+IB ), LDAB-1 )
                  END IF
*
                  IF( I3.GT.0 ) THEN
*
*                    Copy the upper triangle of A31 into the work array.
*
                     DO 110 JJ = 1, IB
                        DO 100 II = 1, MIN( JJ, I3 )
                           WORK( II, JJ ) = AB( KD+1-JJ+II, JJ+I-1 )
  100                   CONTINUE
  110                CONTINUE
*
*                    Update A31 (in the work array).
*
                     CALL DTRSM( 'Right', 'Lower', 'Transpose',
     $                           'Non-unit', I3, IB, ONE, AB( 1, I ),
     $                           LDAB-1, WORK, LDWORK )
*
*                    Update A32
*
                     IF( I2.GT.0 )
     $                  CALL DGEMM( 'No transpose', 'Transpose', I3, I2,
     $                              IB, -ONE, WORK, LDWORK,
     $                              AB( 1+IB, I ), LDAB-1, ONE,
     $                              AB( 1+KD-IB, I+IB ), LDAB-1 )
*
*                    Update A33
*
                     CALL DSYRK( 'Lower', 'No Transpose', I3, IB, -ONE,
     $                           WORK, LDWORK, ONE, AB( 1, I+KD ),
     $                           LDAB-1 )
*
*                    Copy the upper triangle of A31 back into place.
*
                     DO 130 JJ = 1, IB
                        DO 120 II = 1, MIN( JJ, I3 )
                           AB( KD+1-JJ+II, JJ+I-1 ) = WORK( II, JJ )
  120                   CONTINUE
  130                CONTINUE
                  END IF
               END IF
  140       CONTINUE
         END IF
      END IF
      RETURN
*
  150 CONTINUE
      RETURN
*
*     End of DPBTRF
*
      END
      SUBROUTINE DPBTF2( UPLO, N, KD, AB, LDAB, INFO )
*
*  -- LAPACK routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     February 29, 1992
*
*     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            INFO, KD, LDAB, N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   AB( LDAB, * )
*     ..
*
*  Purpose
*  =======
*
*  DPBTF2 computes the Cholesky factorization of a real symmetric
*  positive definite band matrix A.
*
*  The factorization has the form
*     A = U' * U ,  if UPLO = 'U', or
*     A = L  * L',  if UPLO = 'L',
*  where U is an upper triangular matrix, U' is the transpose of U, and
*  L is lower triangular.
*
*  This is the unblocked version of the algorithm, calling Level 2 BLAS.
*
*  Arguments
*  =========
*
*  UPLO    (input) CHARACTER*1
*          Specifies whether the upper or lower triangular part of the
*          symmetric matrix A is stored:
*          = 'U':  Upper triangular
*          = 'L':  Lower triangular
*
*  N       (input) INTEGER
*          The order of the matrix A.  N >= 0.
*
*  KD      (input) INTEGER
*          The number of super-diagonals of the matrix A if UPLO = 'U',
*          or the number of sub-diagonals if UPLO = 'L'.  KD >= 0.
*
*  AB      (input/output) DOUBLE PRECISION array, dimension (LDAB,N)
*          On entry, the upper or lower triangle of the symmetric band
*          matrix A, stored in the first KD+1 rows of the array.  The
*          j-th column of A is stored in the j-th column of the array AB
*          as follows:
*          if UPLO = 'U', AB(kd+1+i-j,j) = A(i,j) for max(1,j-kd)<=i<=j;
*          if UPLO = 'L', AB(1+i-j,j)    = A(i,j) for j<=i<=min(n,j+kd).
*
*          On exit, if INFO = 0, the triangular factor U or L from the
*          Cholesky factorization A = U'*U or A = L*L' of the band
*          matrix A, in the same storage format as A.
*
*  LDAB    (input) INTEGER
*          The leading dimension of the array AB.  LDAB >= KD+1.
*
*  INFO    (output) INTEGER
*          = 0: successful exit
*          < 0: if INFO = -k, the k-th argument had an illegal value
*          > 0: if INFO = k, the leading minor of order k is not
*               positive definite, and the factorization could not be
*               completed.
*
*  Further Details
*  ===============
*
*  The band storage scheme is illustrated by the following example, when
*  N = 6, KD = 2, and UPLO = 'U':
*
*  On entry:                       On exit:
*
*      *    *   a13  a24  a35  a46      *    *   u13  u24  u35  u46
*      *   a12  a23  a34  a45  a56      *   u12  u23  u34  u45  u56
*     a11  a22  a33  a44  a55  a66     u11  u22  u33  u44  u55  u66
*
*  Similarly, if UPLO = 'L' the format of A is as follows:
*
*  On entry:                       On exit:
*
*     a11  a22  a33  a44  a55  a66     l11  l22  l33  l44  l55  l66
*     a21  a32  a43  a54  a65   *      l21  l32  l43  l54  l65   *
*     a31  a42  a53  a64   *    *      l31  l42  l53  l64   *    *
*
*  Array elements marked * are not used by the routine.
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            UPPER
      INTEGER            J, KLD, KN
      DOUBLE PRECISION   AJJ
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. External Subroutines ..
      EXTERNAL           DSCAL, DSYR, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN, SQRT
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( KD.LT.0 ) THEN
         INFO = -3
      ELSE IF( LDAB.LT.KD+1 ) THEN
         INFO = -5
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DPBTF2', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 )
     $   RETURN
*
      KLD = MAX( 1, LDAB-1 )
*
      IF( UPPER ) THEN
*
*        Compute the Cholesky factorization A = U'*U.
*
         DO 10 J = 1, N
*
*           Compute U(J,J) and test for non-positive-definiteness.
*
            AJJ = AB( KD+1, J )
            IF( AJJ.LE.ZERO )
     $         GO TO 30
            AJJ = SQRT( AJJ )
            AB( KD+1, J ) = AJJ
*
*           Compute elements J+1:J+KN of row J and update the
*           trailing submatrix within the band.
*
            KN = MIN( KD, N-J )
            IF( KN.GT.0 ) THEN
               CALL DSCAL( KN, ONE / AJJ, AB( KD, J+1 ), KLD )
               CALL DSYR( 'Upper', KN, -ONE, AB( KD, J+1 ), KLD,
     $                    AB( KD+1, J+1 ), KLD )
            END IF
   10    CONTINUE
      ELSE
*
*        Compute the Cholesky factorization A = L*L'.
*
         DO 20 J = 1, N
*
*           Compute L(J,J) and test for non-positive-definiteness.
*
            AJJ = AB( 1, J )
            IF( AJJ.LE.ZERO )
     $         GO TO 30
            AJJ = SQRT( AJJ )
            AB( 1, J ) = AJJ
*
*           Compute elements J+1:J+KN of column J and update the
*           trailing submatrix within the band.
*
            KN = MIN( KD, N-J )
            IF( KN.GT.0 ) THEN
               CALL DSCAL( KN, ONE / AJJ, AB( 2, J ), 1 )
               CALL DSYR( 'Lower', KN, -ONE, AB( 2, J ), 1,
     $                    AB( 1, J+1 ), KLD )
            END IF
   20    CONTINUE
      END IF
      RETURN
*
   30 CONTINUE
      INFO = J
      RETURN
*
*     End of DPBTF2
*
      END
      SUBROUTINE DPOTF2( UPLO, N, A, LDA, INFO )
*
*  -- LAPACK routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     February 29, 1992
*
*     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            INFO, LDA, N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * )
*     ..
*
*  Purpose
*  =======
*
*  DPOTF2 computes the Cholesky factorization of a real symmetric
*  positive definite matrix A.
*
*  The factorization has the form
*     A = U' * U ,  if UPLO = 'U', or
*     A = L  * L',  if UPLO = 'L',
*  where U is an upper triangular matrix and L is lower triangular.
*
*  This is the unblocked version of the algorithm, calling Level 2 BLAS.
*
*  Arguments
*  =========
*
*  UPLO    (input) CHARACTER*1
*          Specifies whether the upper or lower triangular part of the
*          symmetric matrix A is stored.
*          = 'U':  Upper triangular
*          = 'L':  Lower triangular
*
*  N       (input) INTEGER
*          The order of the matrix A.  N >= 0.
*
*  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
*          On entry, the symmetric matrix A.  If UPLO = 'U', the leading
*          n by n upper triangular part of A contains the upper
*          triangular part of the matrix A, and the strictly lower
*          triangular part of A is not referenced.  If UPLO = 'L', the
*          leading n by n lower triangular part of A contains the lower
*          triangular part of the matrix A, and the strictly upper
*          triangular part of A is not referenced.
*
*          On exit, if INFO = 0, the factor U or L from the Cholesky
*          factorization A = U'*U  or A = L*L'.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,N).
*
*  INFO    (output) INTEGER
*          = 0: successful exit
*          < 0: if INFO = -k, the k-th argument had an illegal value
*          > 0: if INFO = k, the leading minor of order k is not
*               positive definite, and the factorization could not be
*               completed.
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            UPPER
      INTEGER            J
      DOUBLE PRECISION   AJJ
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      DOUBLE PRECISION   DDOT
      EXTERNAL           LSAME, DDOT
*     ..
*     .. External Subroutines ..
      EXTERNAL           DGEMV, DSCAL, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, SQRT
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -4
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DPOTF2', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 )
     $   RETURN
*
      IF( UPPER ) THEN
*
*        Compute the Cholesky factorization A = U'*U.
*
         DO 10 J = 1, N
*
*           Compute U(J,J) and test for non-positive-definiteness.
*
            AJJ = A( J, J ) - DDOT( J-1, A( 1, J ), 1, A( 1, J ), 1 )
            IF( AJJ.LE.ZERO ) THEN
               A( J, J ) = AJJ
               GO TO 30
            END IF
            AJJ = SQRT( AJJ )
            A( J, J ) = AJJ
*
*           Compute elements J+1:N of row J.
*
            IF( J.LT.N ) THEN
               CALL DGEMV( 'Transpose', J-1, N-J, -ONE, A( 1, J+1 ),
     $                     LDA, A( 1, J ), 1, ONE, A( J, J+1 ), LDA )
               CALL DSCAL( N-J, ONE / AJJ, A( J, J+1 ), LDA )
            END IF
   10    CONTINUE
      ELSE
*
*        Compute the Cholesky factorization A = L*L'.
*
         DO 20 J = 1, N
*
*           Compute L(J,J) and test for non-positive-definiteness.
*
            AJJ = A( J, J ) - DDOT( J-1, A( J, 1 ), LDA, A( J, 1 ),
     $            LDA )
            IF( AJJ.LE.ZERO ) THEN
               A( J, J ) = AJJ
               GO TO 30
            END IF
            AJJ = SQRT( AJJ )
            A( J, J ) = AJJ
*
*           Compute elements J+1:N of column J.
*
            IF( J.LT.N ) THEN
               CALL DGEMV( 'No transpose', N-J, J-1, -ONE, A( J+1, 1 ),
     $                     LDA, A( J, 1 ), LDA, ONE, A( J+1, J ), 1 )
               CALL DSCAL( N-J, ONE / AJJ, A( J+1, J ), 1 )
            END IF
   20    CONTINUE
      END IF
      GO TO 40
*
   30 CONTINUE
      INFO = J
*
   40 CONTINUE
      RETURN
*
*     End of DPOTF2
*
      END
      SUBROUTINE DPBTRS( UPLO, N, KD, NRHS, AB, LDAB, B, LDB, INFO )
*
*  -- LAPACK routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     September 30, 1994
*
*     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            INFO, KD, LDAB, LDB, N, NRHS
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   AB( LDAB, * ), B( LDB, * )
*     ..
*
*  Purpose
*  =======
*
*  DPBTRS solves a system of linear equations A*X = B with a symmetric
*  positive definite band matrix A using the Cholesky factorization
*  A = U**T*U or A = L*L**T computed by DPBTRF.
*
*  Arguments
*  =========
*
*  UPLO    (input) CHARACTER*1
*          = 'U':  Upper triangular factor stored in AB;
*          = 'L':  Lower triangular factor stored in AB.
*
*  N       (input) INTEGER
*          The order of the matrix A.  N >= 0.
*
*  KD      (input) INTEGER
*          The number of superdiagonals of the matrix A if UPLO = 'U',
*          or the number of subdiagonals if UPLO = 'L'.  KD >= 0.
*
*  NRHS    (input) INTEGER
*          The number of right hand sides, i.e., the number of columns
*          of the matrix B.  NRHS >= 0.
*
*  AB      (input) DOUBLE PRECISION array, dimension (LDAB,N)
*          The triangular factor U or L from the Cholesky factorization
*          A = U**T*U or A = L*L**T of the band matrix A, stored in the
*          first KD+1 rows of the array.  The j-th column of U or L is
*          stored in the j-th column of the array AB as follows:
*          if UPLO ='U', AB(kd+1+i-j,j) = U(i,j) for max(1,j-kd)<=i<=j;
*          if UPLO ='L', AB(1+i-j,j)    = L(i,j) for j<=i<=min(n,j+kd).
*
*  LDAB    (input) INTEGER
*          The leading dimension of the array AB.  LDAB >= KD+1.
*
*  B       (input/output) DOUBLE PRECISION array, dimension (LDB,NRHS)
*          On entry, the right hand side matrix B.
*          On exit, the solution matrix X.
*
*  LDB     (input) INTEGER
*          The leading dimension of the array B.  LDB >= max(1,N).
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value
*
*  =====================================================================
*
*     .. Local Scalars ..
      LOGICAL            UPPER
      INTEGER            J
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. External Subroutines ..
      EXTERNAL           DTBSV, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( KD.LT.0 ) THEN
         INFO = -3
      ELSE IF( NRHS.LT.0 ) THEN
         INFO = -4
      ELSE IF( LDAB.LT.KD+1 ) THEN
         INFO = -6
      ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
         INFO = -8
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DPBTRS', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 .OR. NRHS.EQ.0 )
     $   RETURN
*
      IF( UPPER ) THEN
*
*        Solve A*X = B where A = U'*U.
*
         DO 10 J = 1, NRHS
*
*           Solve U'*X = B, overwriting B with X.
*
            CALL DTBSV( 'Upper', 'Transpose', 'Non-unit', N, KD, AB,
     $                  LDAB, B( 1, J ), 1 )
*
*           Solve U*X = B, overwriting B with X.
*
            CALL DTBSV( 'Upper', 'No transpose', 'Non-unit', N, KD, AB,
     $                  LDAB, B( 1, J ), 1 )
   10    CONTINUE
      ELSE
*
*        Solve A*X = B where A = L*L'.
*
         DO 20 J = 1, NRHS
*
*           Solve L*X = B, overwriting B with X.
*
            CALL DTBSV( 'Lower', 'No transpose', 'Non-unit', N, KD, AB,
     $                  LDAB, B( 1, J ), 1 )
*
*           Solve L'*X = B, overwriting B with X.
*
            CALL DTBSV( 'Lower', 'Transpose', 'Non-unit', N, KD, AB,
     $                  LDAB, B( 1, J ), 1 )
   20    CONTINUE
      END IF
*
      RETURN
*
*     End of DPBTRS
*
      END
      SUBROUTINE DPBCON( UPLO, N, KD, AB, LDAB, ANORM, RCOND, WORK,
     $                   IWORK, INFO )
*
*  -- LAPACK routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     September 30, 1994
*
*     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            INFO, KD, LDAB, N
      DOUBLE PRECISION   ANORM, RCOND
*     ..
*     .. Array Arguments ..
      INTEGER            IWORK( * )
      DOUBLE PRECISION   AB( LDAB, * ), WORK( * )
*     ..
*
*  Purpose
*  =======
*
*  DPBCON estimates the reciprocal of the condition number (in the
*  1-norm) of a real symmetric positive definite band matrix using the
*  Cholesky factorization A = U**T*U or A = L*L**T computed by DPBTRF.
*
*  An estimate is obtained for norm(inv(A)), and the reciprocal of the
*  condition number is computed as RCOND = 1 / (ANORM * norm(inv(A))).
*
*  Arguments
*  =========
*
*  UPLO    (input) CHARACTER*1
*          = 'U':  Upper triangular factor stored in AB;
*          = 'L':  Lower triangular factor stored in AB.
*
*  N       (input) INTEGER
*          The order of the matrix A.  N >= 0.
*
*  KD      (input) INTEGER
*          The number of superdiagonals of the matrix A if UPLO = 'U',
*          or the number of subdiagonals if UPLO = 'L'.  KD >= 0.
*
*  AB      (input) DOUBLE PRECISION array, dimension (LDAB,N)
*          The triangular factor U or L from the Cholesky factorization
*          A = U**T*U or A = L*L**T of the band matrix A, stored in the
*          first KD+1 rows of the array.  The j-th column of U or L is
*          stored in the j-th column of the array AB as follows:
*          if UPLO ='U', AB(kd+1+i-j,j) = U(i,j) for max(1,j-kd)<=i<=j;
*          if UPLO ='L', AB(1+i-j,j)    = L(i,j) for j<=i<=min(n,j+kd).
*
*  LDAB    (input) INTEGER
*          The leading dimension of the array AB.  LDAB >= KD+1.
*
*  ANORM   (input) DOUBLE PRECISION
*          The 1-norm (or infinity-norm) of the symmetric band matrix A.
*
*  RCOND   (output) DOUBLE PRECISION
*          The reciprocal of the condition number of the matrix A,
*          computed as RCOND = 1/(ANORM * AINVNM), where AINVNM is an
*          estimate of the 1-norm of inv(A) computed in this routine.
*
*  WORK    (workspace) DOUBLE PRECISION array, dimension (3*N)
*
*  IWORK   (workspace) INTEGER array, dimension (N)
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            UPPER
      CHARACTER          NORMIN
      INTEGER            IX, KASE
      DOUBLE PRECISION   AINVNM, SCALE, SCALEL, SCALEU, SMLNUM
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            IDAMAX
      DOUBLE PRECISION   DLAMCH
      EXTERNAL           LSAME, IDAMAX, DLAMCH
*     ..
*     .. External Subroutines ..
      EXTERNAL           DLACON, DLATBS, DRSCL, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( KD.LT.0 ) THEN
         INFO = -3
      ELSE IF( LDAB.LT.KD+1 ) THEN
         INFO = -5
      ELSE IF( ANORM.LT.ZERO ) THEN
         INFO = -6
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DPBCON', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      RCOND = ZERO
      IF( N.EQ.0 ) THEN
         RCOND = ONE
         RETURN
      ELSE IF( ANORM.EQ.ZERO ) THEN
         RETURN
      END IF
*
      SMLNUM = DLAMCH( 'Safe minimum' )
*
*     Estimate the 1-norm of the inverse.
*
      KASE = 0
      NORMIN = 'N'
   10 CONTINUE
      CALL DLACON( N, WORK( N+1 ), WORK, IWORK, AINVNM, KASE )
      IF( KASE.NE.0 ) THEN
         IF( UPPER ) THEN
*
*           Multiply by inv(U').
*
            CALL DLATBS( 'Upper', 'Transpose', 'Non-unit', NORMIN, N,
     $                   KD, AB, LDAB, WORK, SCALEL, WORK( 2*N+1 ),
     $                   INFO )
            NORMIN = 'Y'
*
*           Multiply by inv(U).
*
            CALL DLATBS( 'Upper', 'No transpose', 'Non-unit', NORMIN, N,
     $                   KD, AB, LDAB, WORK, SCALEU, WORK( 2*N+1 ),
     $                   INFO )
         ELSE
*
*           Multiply by inv(L).
*
            CALL DLATBS( 'Lower', 'No transpose', 'Non-unit', NORMIN, N,
     $                   KD, AB, LDAB, WORK, SCALEL, WORK( 2*N+1 ),
     $                   INFO )
            NORMIN = 'Y'
*
*           Multiply by inv(L').
*
            CALL DLATBS( 'Lower', 'Transpose', 'Non-unit', NORMIN, N,
     $                   KD, AB, LDAB, WORK, SCALEU, WORK( 2*N+1 ),
     $                   INFO )
         END IF
*
*        Multiply by 1/SCALE if doing so will not cause overflow.
*
         SCALE = SCALEL*SCALEU
         IF( SCALE.NE.ONE ) THEN
            IX = IDAMAX( N, WORK, 1 )
            IF( SCALE.LT.ABS( WORK( IX ) )*SMLNUM .OR. SCALE.EQ.ZERO )
     $         GO TO 20
            CALL DRSCL( N, SCALE, WORK, 1 )
         END IF
         GO TO 10
      END IF
*
*     Compute the estimate of the reciprocal condition number.
*
      IF( AINVNM.NE.ZERO )
     $   RCOND = ( ONE / AINVNM ) / ANORM
*
   20 CONTINUE
*
      RETURN
*
*     End of DPBCON
*
      END
      SUBROUTINE DLABAD( SMALL, LARGE )
*
*  -- LAPACK auxiliary routine (version 2.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     October 31, 1992
*
*     .. Scalar Arguments ..
      DOUBLE PRECISION   LARGE, SMALL
*     ..
*
*  Purpose
*  =======
*
*  DLABAD takes as input the values computed by SLAMCH for underflow and
*  overflow, and returns the square root of each of these values if the
*  log of LARGE is sufficiently large.  This subroutine is intended to
*  identify machines with a large exponent range, such as the Crays, and
*  redefine the underflow and overflow limits to be the square roots of
*  the values computed by DLAMCH.  This subroutine is needed because
*  DLAMCH does not compensate for poor arithmetic in the upper half of
*  the exponent range, as is found on a Cray.
*
*  Arguments
*  =========
*
*  SMALL   (input/output) DOUBLE PRECISION
*          On entry, the underflow threshold as computed by DLAMCH.
*          On exit, if LOG10(LARGE) is sufficiently large, the square
*          root of SMALL, otherwise unchanged.
*
*  LARGE   (input/output) DOUBLE PRECISION
*          On entry, the overflow threshold as computed by DLAMCH.
*          On exit, if LOG10(LARGE) is sufficiently large, the square
*          root of LARGE, otherwise unchanged.
*
*  =====================================================================
*
*     .. Intrinsic Functions ..
      INTRINSIC          LOG10, SQRT
*     ..
*     .. Executable Statements ..
*
*     If it looks like we're on a Cray, take the square root of
*     SMALL and LARGE to avoid overflow and underflow problems.
*
      IF( LOG10( LARGE ).GT.2000.D0 ) THEN
         SMALL = SQRT( SMALL )
         LARGE = SQRT( LARGE )
      END IF
*
      RETURN
*
*     End of DLABAD
*
      END
      SUBROUTINE DLACON( N, V, X, ISGN, EST, KASE )
*
*  -- LAPACK auxiliary routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     February 29, 1992
*
*     .. Scalar Arguments ..
      INTEGER            KASE, N
      DOUBLE PRECISION   EST
*     ..
*     .. Array Arguments ..
      INTEGER            ISGN( * )
      DOUBLE PRECISION   V( * ), X( * )
*     ..
*
*  Purpose
*  =======
*
*  DLACON estimates the 1-norm of a square, real matrix A.
*  Reverse communication is used for evaluating matrix-vector products.
*
*  Arguments
*  =========
*
*  N      (input) INTEGER
*         The order of the matrix.  N >= 1.
*
*  V      (workspace) DOUBLE PRECISION array, dimension (N)
*         On the final return, V = A*W,  where  EST = norm(V)/norm(W)
*         (W is not returned).
*
*  X      (input/output) DOUBLE PRECISION array, dimension (N)
*         On an intermediate return, X should be overwritten by
*               A * X,   if KASE=1,
*               A' * X,  if KASE=2,
*         and DLACON must be re-called with all the other parameters
*         unchanged.
*
*  ISGN   (workspace) INTEGER array, dimension (N)
*
*  EST    (output) DOUBLE PRECISION
*         An estimate (a lower bound) for norm(A).
*
*  KASE   (input/output) INTEGER
*         On the initial call to DLACON, KASE should be 0.
*         On an intermediate return, KASE will be 1 or 2, indicating
*         whether X should be overwritten by A * X  or A' * X.
*         On the final return from DLACON, KASE will again be 0.
*
*  Further Details
*  ======= =======
*
*  Contributed by Nick Higham, University of Manchester.
*  Originally named SONEST, dated March 16, 1988.
*
*  Reference: N.J. Higham, "FORTRAN codes for estimating the one-norm of
*  a real or complex matrix, with applications to condition estimation",
*  ACM Trans. Math. Soft., vol. 14, no. 4, pp. 381-396, December 1988.
*
*  =====================================================================
*
*     .. Parameters ..
      INTEGER            ITMAX
      PARAMETER          ( ITMAX = 5 )
      DOUBLE PRECISION   ZERO, ONE, TWO
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0, TWO = 2.0D+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            I, ITER, J, JLAST, JUMP
      DOUBLE PRECISION   ALTSGN, ESTOLD, TEMP
*     ..
*     .. External Functions ..
      INTEGER            IDAMAX
      DOUBLE PRECISION   DASUM
      EXTERNAL           IDAMAX, DASUM
*     ..
*     .. External Subroutines ..
      EXTERNAL           DCOPY
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, DBLE, NINT, SIGN
*     ..
*     .. Save statement ..
      SAVE
*     ..
*     .. Executable Statements ..
*
      IF( KASE.EQ.0 ) THEN
         DO 10 I = 1, N
            X( I ) = ONE / DBLE( N )
   10    CONTINUE
         KASE = 1
         JUMP = 1
         RETURN
      END IF
*
      GO TO ( 20, 40, 70, 110, 140 )JUMP
*
*     ................ ENTRY   (JUMP = 1)
*     FIRST ITERATION.  X HAS BEEN OVERWRITTEN BY A*X.
*
   20 CONTINUE
      IF( N.EQ.1 ) THEN
         V( 1 ) = X( 1 )
         EST = ABS( V( 1 ) )
*        ... QUIT
         GO TO 150
      END IF
      EST = DASUM( N, X, 1 )
*
      DO 30 I = 1, N
         X( I ) = SIGN( ONE, X( I ) )
         ISGN( I ) = NINT( X( I ) )
   30 CONTINUE
      KASE = 2
      JUMP = 2
      RETURN
*
*     ................ ENTRY   (JUMP = 2)
*     FIRST ITERATION.  X HAS BEEN OVERWRITTEN BY TRANDPOSE(A)*X.
*
   40 CONTINUE
      J = IDAMAX( N, X, 1 )
      ITER = 2
*
*     MAIN LOOP - ITERATIONS 2,3,...,ITMAX.
*
   50 CONTINUE
      DO 60 I = 1, N
         X( I ) = ZERO
   60 CONTINUE
      X( J ) = ONE
      KASE = 1
      JUMP = 3
      RETURN
*
*     ................ ENTRY   (JUMP = 3)
*     X HAS BEEN OVERWRITTEN BY A*X.
*
   70 CONTINUE
      CALL DCOPY( N, X, 1, V, 1 )
      ESTOLD = EST
      EST = DASUM( N, V, 1 )
      DO 80 I = 1, N
         IF( NINT( SIGN( ONE, X( I ) ) ).NE.ISGN( I ) )
     $      GO TO 90
   80 CONTINUE
*     REPEATED SIGN VECTOR DETECTED, HENCE ALGORITHM HAS CONVERGED.
      GO TO 120
*
   90 CONTINUE
*     TEST FOR CYCLING.
      IF( EST.LE.ESTOLD )
     $   GO TO 120
*
      DO 100 I = 1, N
         X( I ) = SIGN( ONE, X( I ) )
         ISGN( I ) = NINT( X( I ) )
  100 CONTINUE
      KASE = 2
      JUMP = 4
      RETURN
*
*     ................ ENTRY   (JUMP = 4)
*     X HAS BEEN OVERWRITTEN BY TRANDPOSE(A)*X.
*
  110 CONTINUE
      JLAST = J
      J = IDAMAX( N, X, 1 )
      IF( ( X( JLAST ).NE.ABS( X( J ) ) ) .AND. ( ITER.LT.ITMAX ) ) THEN
         ITER = ITER + 1
         GO TO 50
      END IF
*
*     ITERATION COMPLETE.  FINAL STAGE.
*
  120 CONTINUE
      ALTSGN = ONE
      DO 130 I = 1, N
         X( I ) = ALTSGN*( ONE+DBLE( I-1 ) / DBLE( N-1 ) )
         ALTSGN = -ALTSGN
  130 CONTINUE
      KASE = 1
      JUMP = 5
      RETURN
*
*     ................ ENTRY   (JUMP = 5)
*     X HAS BEEN OVERWRITTEN BY A*X.
*
  140 CONTINUE
      TEMP = TWO*( DASUM( N, X, 1 ) / DBLE( 3*N ) )
      IF( TEMP.GT.EST ) THEN
         CALL DCOPY( N, X, 1, V, 1 )
         EST = TEMP
      END IF
*
  150 CONTINUE
      KASE = 0
      RETURN
*
*     End of DLACON
*
      END
      DOUBLE PRECISION FUNCTION DLAMCH( CMACH )
*
*  -- LAPACK auxiliary routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     October 31, 1992
*
*     .. Scalar Arguments ..
      CHARACTER          CMACH
*     ..
*
*  Purpose
*  =======
*
*  DLAMCH determines double precision machine parameters.
*
*  Arguments
*  =========
*
*  CMACH   (input) CHARACTER*1
*          Specifies the value to be returned by DLAMCH:
*          = 'E' or 'e',   DLAMCH := eps
*          = 'S' or 's ,   DLAMCH := sfmin
*          = 'B' or 'b',   DLAMCH := base
*          = 'P' or 'p',   DLAMCH := eps*base
*          = 'N' or 'n',   DLAMCH := t
*          = 'R' or 'r',   DLAMCH := rnd
*          = 'M' or 'm',   DLAMCH := emin
*          = 'U' or 'u',   DLAMCH := rmin
*          = 'L' or 'l',   DLAMCH := emax
*          = 'O' or 'o',   DLAMCH := rmax
*
*          where
*
*          eps   = relative machine precision
*          sfmin = safe minimum, such that 1/sfmin does not overflow
*          base  = base of the machine
*          prec  = eps*base
*          t     = number of (base) digits in the mantissa
*          rnd   = 1.0 when rounding occurs in addition, 0.0 otherwise
*          emin  = minimum exponent before (gradual) underflow
*          rmin  = underflow threshold - base**(emin-1)
*          emax  = largest exponent before overflow
*          rmax  = overflow threshold  - (base**emax)*(1-eps)
*
* =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            FIRST, LRND
      INTEGER            BETA, IMAX, IMIN, IT
      DOUBLE PRECISION   BASE, EMAX, EMIN, EPS, PREC, RMACH, RMAX, RMIN,
     $                   RND, SFMIN, SMALL, T
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. External Subroutines ..
      EXTERNAL           DLAMC2
*     ..
*     .. Save statement ..
      SAVE               FIRST, EPS, SFMIN, BASE, T, RND, EMIN, RMIN,
     $                   EMAX, RMAX, PREC
*     ..
*     .. Data statements ..
      DATA               FIRST / .TRUE. /
*     ..
*     .. Executable Statements ..
*
      IF( FIRST ) THEN
         FIRST = .FALSE.
         CALL DLAMC2( BETA, IT, LRND, EPS, IMIN, RMIN, IMAX, RMAX )
         BASE = BETA
         T = IT
         IF( LRND ) THEN
            RND = ONE
            EPS = ( BASE**( 1-IT ) ) / 2
         ELSE
            RND = ZERO
            EPS = BASE**( 1-IT )
         END IF
         PREC = EPS*BASE
         EMIN = IMIN
         EMAX = IMAX
         SFMIN = RMIN
         SMALL = ONE / RMAX
         IF( SMALL.GE.SFMIN ) THEN
*
*           Use SMALL plus a bit, to avoid the possibility of rounding
*           causing overflow when computing  1/sfmin.
*
            SFMIN = SMALL*( ONE+EPS )
         END IF
      END IF
*
      IF( LSAME( CMACH, 'E' ) ) THEN
         RMACH = EPS
      ELSE IF( LSAME( CMACH, 'S' ) ) THEN
         RMACH = SFMIN
      ELSE IF( LSAME( CMACH, 'B' ) ) THEN
         RMACH = BASE
      ELSE IF( LSAME( CMACH, 'P' ) ) THEN
         RMACH = PREC
      ELSE IF( LSAME( CMACH, 'N' ) ) THEN
         RMACH = T
      ELSE IF( LSAME( CMACH, 'R' ) ) THEN
         RMACH = RND
      ELSE IF( LSAME( CMACH, 'M' ) ) THEN
         RMACH = EMIN
      ELSE IF( LSAME( CMACH, 'U' ) ) THEN
         RMACH = RMIN
      ELSE IF( LSAME( CMACH, 'L' ) ) THEN
         RMACH = EMAX
      ELSE IF( LSAME( CMACH, 'O' ) ) THEN
         RMACH = RMAX
      END IF
*
      DLAMCH = RMACH
      RETURN
*
*     End of DLAMCH
*
      END
*
************************************************************************
*
      SUBROUTINE DLAMC1( BETA, T, RND, IEEE1 )
*
*  -- LAPACK auxiliary routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     October 31, 1992
*
*     .. Scalar Arguments ..
      LOGICAL            IEEE1, RND
      INTEGER            BETA, T
*     ..
*
*  Purpose
*  =======
*
*  DLAMC1 determines the machine parameters given by BETA, T, RND, and
*  IEEE1.
*
*  Arguments
*  =========
*
*  BETA    (output) INTEGER
*          The base of the machine.
*
*  T       (output) INTEGER
*          The number of ( BETA ) digits in the mantissa.
*
*  RND     (output) LOGICAL
*          Specifies whether proper rounding  ( RND = .TRUE. )  or
*          chopping  ( RND = .FALSE. )  occurs in addition. This may not
*          be a reliable guide to the way in which the machine performs
*          its arithmetic.
*
*  IEEE1   (output) LOGICAL
*          Specifies whether rounding appears to be done in the IEEE
*          'round to nearest' style.
*
*  Further Details
*  ===============
*
*  The routine is based on the routine  ENVRON  by Malcolm and
*  incorporates suggestions by Gentleman and Marovich. See
*
*     Malcolm M. A. (1972) Algorithms to reveal properties of
*        floating-point arithmetic. Comms. of the ACM, 15, 949-951.
*
*     Gentleman W. M. and Marovich S. B. (1974) More on algorithms
*        that reveal properties of floating point arithmetic units.
*        Comms. of the ACM, 17, 276-277.
*
* =====================================================================
*
*     .. Local Scalars ..
      LOGICAL            FIRST, LIEEE1, LRND
      INTEGER            LBETA, LT
      DOUBLE PRECISION   A, B, C, F, ONE, QTR, SAVEC, T1, T2
*     ..
*     .. External Functions ..
      DOUBLE PRECISION   DLAMC3
      EXTERNAL           DLAMC3
*     ..
*     .. Save statement ..
      SAVE               FIRST, LIEEE1, LBETA, LRND, LT
*     ..
*     .. Data statements ..
      DATA               FIRST / .TRUE. /
*     ..
*     .. Executable Statements ..
*
      IF( FIRST ) THEN
         FIRST = .FALSE.
         ONE = 1
*
*        LBETA,  LIEEE1,  LT and  LRND  are the  local values  of  BETA,
*        IEEE1, T and RND.
*
*        Throughout this routine  we use the function  DLAMC3  to ensure
*        that relevant values are  stored and not held in registers,  or
*        are not affected by optimizers.
*
*        Compute  a = 2.0**m  with the  smallest positive integer m such
*        that
*
*           fl( a + 1.0 ) = a.
*
         A = 1
         C = 1
*
*+       WHILE( C.EQ.ONE )LOOP
   10    CONTINUE
         IF( C.EQ.ONE ) THEN
            A = 2*A
            C = DLAMC3( A, ONE )
            C = DLAMC3( C, -A )
            GO TO 10
         END IF
*+       END WHILE
*
*        Now compute  b = 2.0**m  with the smallest positive integer m
*        such that
*
*           fl( a + b ) .gt. a.
*
         B = 1
         C = DLAMC3( A, B )
*
*+       WHILE( C.EQ.A )LOOP
   20    CONTINUE
         IF( C.EQ.A ) THEN
            B = 2*B
            C = DLAMC3( A, B )
            GO TO 20
         END IF
*+       END WHILE
*
*        Now compute the base.  a and c  are neighbouring floating point
*        numbers  in the  interval  ( beta**t, beta**( t + 1 ) )  and so
*        their difference is beta. Adding 0.25 to c is to ensure that it
*        is truncated to beta and not ( beta - 1 ).
*
         QTR = ONE / 4
         SAVEC = C
         C = DLAMC3( C, -A )
         LBETA = C + QTR
*
*        Now determine whether rounding or chopping occurs,  by adding a
*        bit  less  than  beta/2  and a  bit  more  than  beta/2  to  a.
*
         B = LBETA
         F = DLAMC3( B / 2, -B / 100 )
         C = DLAMC3( F, A )
         IF( C.EQ.A ) THEN
            LRND = .TRUE.
         ELSE
            LRND = .FALSE.
         END IF
         F = DLAMC3( B / 2, B / 100 )
         C = DLAMC3( F, A )
         IF( ( LRND ) .AND. ( C.EQ.A ) )
     $      LRND = .FALSE.
*
*        Try and decide whether rounding is done in the  IEEE  'round to
*        nearest' style. B/2 is half a unit in the last place of the two
*        numbers A and SAVEC. Furthermore, A is even, i.e. has last  bit
*        zero, and SAVEC is odd. Thus adding B/2 to A should not  change
*        A, but adding B/2 to SAVEC should change SAVEC.
*
         T1 = DLAMC3( B / 2, A )
         T2 = DLAMC3( B / 2, SAVEC )
         LIEEE1 = ( T1.EQ.A ) .AND. ( T2.GT.SAVEC ) .AND. LRND
*
*        Now find  the  mantissa, t.  It should  be the  integer part of
*        log to the base beta of a,  however it is safer to determine  t
*        by powering.  So we find t as the smallest positive integer for
*        which
*
*           fl( beta**t + 1.0 ) = 1.0.
*
         LT = 0
         A = 1
         C = 1
*
*+       WHILE( C.EQ.ONE )LOOP
   30    CONTINUE
         IF( C.EQ.ONE ) THEN
            LT = LT + 1
            A = A*LBETA
            C = DLAMC3( A, ONE )
            C = DLAMC3( C, -A )
            GO TO 30
         END IF
*+       END WHILE
*
      END IF
*
      BETA = LBETA
      T = LT
      RND = LRND
      IEEE1 = LIEEE1
      RETURN
*
*     End of DLAMC1
*
      END
*
************************************************************************
*
      SUBROUTINE DLAMC2( BETA, T, RND, EPS, EMIN, RMIN, EMAX, RMAX )
*
*  -- LAPACK auxiliary routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     October 31, 1992
*
*     .. Scalar Arguments ..
      LOGICAL            RND
      INTEGER            BETA, EMAX, EMIN, T
      DOUBLE PRECISION   EPS, RMAX, RMIN
*     ..
*
*  Purpose
*  =======
*
*  DLAMC2 determines the machine parameters specified in its argument
*  list.
*
*  Arguments
*  =========
*
*  BETA    (output) INTEGER
*          The base of the machine.
*
*  T       (output) INTEGER
*          The number of ( BETA ) digits in the mantissa.
*
*  RND     (output) LOGICAL
*          Specifies whether proper rounding  ( RND = .TRUE. )  or
*          chopping  ( RND = .FALSE. )  occurs in addition. This may not
*          be a reliable guide to the way in which the machine performs
*          its arithmetic.
*
*  EPS     (output) DOUBLE PRECISION
*          The smallest positive number such that
*
*             fl( 1.0 - EPS ) .LT. 1.0,
*
*          where fl denotes the computed value.
*
*  EMIN    (output) INTEGER
*          The minimum exponent before (gradual) underflow occurs.
*
*  RMIN    (output) DOUBLE PRECISION
*          The smallest normalized number for the machine, given by
*          BASE**( EMIN - 1 ), where  BASE  is the floating point value
*          of BETA.
*
*  EMAX    (output) INTEGER
*          The maximum exponent before overflow occurs.
*
*  RMAX    (output) DOUBLE PRECISION
*          The largest positive number for the machine, given by
*          BASE**EMAX * ( 1 - EPS ), where  BASE  is the floating point
*          value of BETA.
*
*  Further Details
*  ===============
*
*  The computation of  EPS  is based on a routine PARANOIA by
*  W. Kahan of the University of California at Berkeley.
*
* =====================================================================
*
*     .. Local Scalars ..
      LOGICAL            FIRST, IEEE, IWARN, LIEEE1, LRND
      INTEGER            GNMIN, GPMIN, I, LBETA, LEMAX, LEMIN, LT,
     $                   NGNMIN, NGPMIN
      DOUBLE PRECISION   A, B, C, HALF, LEPS, LRMAX, LRMIN, ONE, RBASE,
     $                   SIXTH, SMALL, THIRD, TWO, ZERO
*     ..
*     .. External Functions ..
      DOUBLE PRECISION   DLAMC3
      EXTERNAL           DLAMC3
*     ..
*     .. External Subroutines ..
      EXTERNAL           DLAMC1, DLAMC4, DLAMC5
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, MIN
*     ..
*     .. Save statement ..
      SAVE               FIRST, IWARN, LBETA, LEMAX, LEMIN, LEPS, LRMAX,
     $                   LRMIN, LT
*     ..
*     .. Data statements ..
      DATA               FIRST / .TRUE. / , IWARN / .FALSE. /
*     ..
*     .. Executable Statements ..
*
      IF( FIRST ) THEN
         FIRST = .FALSE.
         ZERO = 0
         ONE = 1
         TWO = 2
*
*        LBETA, LT, LRND, LEPS, LEMIN and LRMIN  are the local values of
*        BETA, T, RND, EPS, EMIN and RMIN.
*
*        Throughout this routine  we use the function  DLAMC3  to ensure
*        that relevant values are stored  and not held in registers,  or
*        are not affected by optimizers.
*
*        DLAMC1 returns the parameters  LBETA, LT, LRND and LIEEE1.
*
         CALL DLAMC1( LBETA, LT, LRND, LIEEE1 )
*
*        Start to find EPS.
*
         B = LBETA
         A = B**( -LT )
         LEPS = A
*
*        Try some tricks to see whether or not this is the correct  EPS.
*
         B = TWO / 3
         HALF = ONE / 2
         SIXTH = DLAMC3( B, -HALF )
         THIRD = DLAMC3( SIXTH, SIXTH )
         B = DLAMC3( THIRD, -HALF )
         B = DLAMC3( B, SIXTH )
         B = ABS( B )
         IF( B.LT.LEPS )
     $      B = LEPS
*
         LEPS = 1
*
*+       WHILE( ( LEPS.GT.B ).AND.( B.GT.ZERO ) )LOOP
   10    CONTINUE
         IF( ( LEPS.GT.B ) .AND. ( B.GT.ZERO ) ) THEN
            LEPS = B
            C = DLAMC3( HALF*LEPS, ( TWO**5 )*( LEPS**2 ) )
            C = DLAMC3( HALF, -C )
            B = DLAMC3( HALF, C )
            C = DLAMC3( HALF, -B )
            B = DLAMC3( HALF, C )
            GO TO 10
         END IF
*+       END WHILE
*
         IF( A.LT.LEPS )
     $      LEPS = A
*
*        Computation of EPS complete.
*
*        Now find  EMIN.  Let A = + or - 1, and + or - (1 + BASE**(-3)).
*        Keep dividing  A by BETA until (gradual) underflow occurs. This
*        is detected when we cannot recover the previous A.
*
         RBASE = ONE / LBETA
         SMALL = ONE
         DO 20 I = 1, 3
            SMALL = DLAMC3( SMALL*RBASE, ZERO )
   20    CONTINUE
         A = DLAMC3( ONE, SMALL )
         CALL DLAMC4( NGPMIN, ONE, LBETA )
         CALL DLAMC4( NGNMIN, -ONE, LBETA )
         CALL DLAMC4( GPMIN, A, LBETA )
         CALL DLAMC4( GNMIN, -A, LBETA )
         IEEE = .FALSE.
*
         IF( ( NGPMIN.EQ.NGNMIN ) .AND. ( GPMIN.EQ.GNMIN ) ) THEN
            IF( NGPMIN.EQ.GPMIN ) THEN
               LEMIN = NGPMIN
*            ( Non twos-complement machines, no gradual underflow;
*              e.g.,  VAX )
            ELSE IF( ( GPMIN-NGPMIN ).EQ.3 ) THEN
               LEMIN = NGPMIN - 1 + LT
               IEEE = .TRUE.
*            ( Non twos-complement machines, with gradual underflow;
*              e.g., IEEE standard followers )
            ELSE
               LEMIN = MIN( NGPMIN, GPMIN )
*            ( A guess; no known machine )
               IWARN = .TRUE.
            END IF
*
         ELSE IF( ( NGPMIN.EQ.GPMIN ) .AND. ( NGNMIN.EQ.GNMIN ) ) THEN
            IF( ABS( NGPMIN-NGNMIN ).EQ.1 ) THEN
               LEMIN = MAX( NGPMIN, NGNMIN )
*            ( Twos-complement machines, no gradual underflow;
*              e.g., CYBER 205 )
            ELSE
               LEMIN = MIN( NGPMIN, NGNMIN )
*            ( A guess; no known machine )
               IWARN = .TRUE.
            END IF
*
         ELSE IF( ( ABS( NGPMIN-NGNMIN ).EQ.1 ) .AND.
     $            ( GPMIN.EQ.GNMIN ) ) THEN
            IF( ( GPMIN-MIN( NGPMIN, NGNMIN ) ).EQ.3 ) THEN
               LEMIN = MAX( NGPMIN, NGNMIN ) - 1 + LT
*            ( Twos-complement machines with gradual underflow;
*              no known machine )
            ELSE
               LEMIN = MIN( NGPMIN, NGNMIN )
*            ( A guess; no known machine )
               IWARN = .TRUE.
            END IF
*
         ELSE
            LEMIN = MIN( NGPMIN, NGNMIN, GPMIN, GNMIN )
*         ( A guess; no known machine )
            IWARN = .TRUE.
         END IF
***
* Comment out this if block if EMIN is ok
         IF( IWARN ) THEN
            FIRST = .TRUE.
            WRITE( 6, FMT = 9999 )LEMIN
         END IF
***
*
*        Assume IEEE arithmetic if we found denormalised  numbers above,
*        or if arithmetic seems to round in the  IEEE style,  determined
*        in routine DLAMC1. A true IEEE machine should have both  things
*        true; however, faulty machines may have one or the other.
*
         IEEE = IEEE .OR. LIEEE1
*
*        Compute  RMIN by successive division by  BETA. We could compute
*        RMIN as BASE**( EMIN - 1 ),  but some machines underflow during
*        this computation.
*
         LRMIN = 1
         DO 30 I = 1, 1 - LEMIN
            LRMIN = DLAMC3( LRMIN*RBASE, ZERO )
   30    CONTINUE
*
*        Finally, call DLAMC5 to compute EMAX and RMAX.
*
         CALL DLAMC5( LBETA, LT, LEMIN, IEEE, LEMAX, LRMAX )
      END IF
*
      BETA = LBETA
      T = LT
      RND = LRND
      EPS = LEPS
      EMIN = LEMIN
      RMIN = LRMIN
      EMAX = LEMAX
      RMAX = LRMAX
*
      RETURN
*
 9999 FORMAT( / / ' WARNING. The value EMIN may be incorrect:-',
     $      '  EMIN = ', I8, /
     $      ' If, after inspection, the value EMIN looks',
     $      ' acceptable please comment out ',
     $      / ' the IF block as marked within the code of routine',
     $      ' DLAMC2,', / ' otherwise supply EMIN explicitly.', / )
*
*     End of DLAMC2
*
      END
*
************************************************************************
*
      DOUBLE PRECISION FUNCTION DLAMC3( A, B )
*
*  -- LAPACK auxiliary routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     October 31, 1992
*
*     .. Scalar Arguments ..
      DOUBLE PRECISION   A, B
*     ..
*
*  Purpose
*  =======
*
*  DLAMC3  is intended to force  A  and  B  to be stored prior to doing
*  the addition of  A  and  B ,  for use in situations where optimizers
*  might hold one of these in a register.
*
*  Arguments
*  =========
*
*  A, B    (input) DOUBLE PRECISION
*          The values A and B.
*
* =====================================================================
*
*     .. Executable Statements ..
*
      DLAMC3 = A + B
*
      RETURN
*
*     End of DLAMC3
*
      END
*
************************************************************************
*
      SUBROUTINE DLAMC4( EMIN, START, BASE )
*
*  -- LAPACK auxiliary routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     October 31, 1992
*
*     .. Scalar Arguments ..
      INTEGER            BASE, EMIN
      DOUBLE PRECISION   START
*     ..
*
*  Purpose
*  =======
*
*  DLAMC4 is a service routine for DLAMC2.
*
*  Arguments
*  =========
*
*  EMIN    (output) EMIN
*          The minimum exponent before (gradual) underflow, computed by
*          setting A = START and dividing by BASE until the previous A
*          can not be recovered.
*
*  START   (input) DOUBLE PRECISION
*          The starting point for determining EMIN.
*
*  BASE    (input) INTEGER
*          The base of the machine.
*
* =====================================================================
*
*     .. Local Scalars ..
      INTEGER            I
      DOUBLE PRECISION   A, B1, B2, C1, C2, D1, D2, ONE, RBASE, ZERO
*     ..
*     .. External Functions ..
      DOUBLE PRECISION   DLAMC3
      EXTERNAL           DLAMC3
*     ..
*     .. Executable Statements ..
*
      A = START
      ONE = 1
      RBASE = ONE / BASE
      ZERO = 0
      EMIN = 1
      B1 = DLAMC3( A*RBASE, ZERO )
      C1 = A
      C2 = A
      D1 = A
      D2 = A
*+    WHILE( ( C1.EQ.A ).AND.( C2.EQ.A ).AND.
*    $       ( D1.EQ.A ).AND.( D2.EQ.A )      )LOOP
   10 CONTINUE
      IF( ( C1.EQ.A ) .AND. ( C2.EQ.A ) .AND. ( D1.EQ.A ) .AND.
     $    ( D2.EQ.A ) ) THEN
         EMIN = EMIN - 1
         A = B1
         B1 = DLAMC3( A / BASE, ZERO )
         C1 = DLAMC3( B1*BASE, ZERO )
         D1 = ZERO
         DO 20 I = 1, BASE
            D1 = D1 + B1
   20    CONTINUE
         B2 = DLAMC3( A*RBASE, ZERO )
         C2 = DLAMC3( B2 / RBASE, ZERO )
         D2 = ZERO
         DO 30 I = 1, BASE
            D2 = D2 + B2
   30    CONTINUE
         GO TO 10
      END IF
*+    END WHILE
*
      RETURN
*
*     End of DLAMC4
*
      END
*
************************************************************************
*
      SUBROUTINE DLAMC5( BETA, P, EMIN, IEEE, EMAX, RMAX )
*
*  -- LAPACK auxiliary routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     October 31, 1992
*
*     .. Scalar Arguments ..
      LOGICAL            IEEE
      INTEGER            BETA, EMAX, EMIN, P
      DOUBLE PRECISION   RMAX
*     ..
*
*  Purpose
*  =======
*
*  DLAMC5 attempts to compute RMAX, the largest machine floating-point
*  number, without overflow.  It assumes that EMAX + abs(EMIN) sum
*  approximately to a power of 2.  It will fail on machines where this
*  assumption does not hold, for example, the Cyber 205 (EMIN = -28625,
*  EMAX = 28718).  It will also fail if the value supplied for EMIN is
*  too large (i.e. too close to zero), probably with overflow.
*
*  Arguments
*  =========
*
*  BETA    (input) INTEGER
*          The base of floating-point arithmetic.
*
*  P       (input) INTEGER
*          The number of base BETA digits in the mantissa of a
*          floating-point value.
*
*  EMIN    (input) INTEGER
*          The minimum exponent before (gradual) underflow.
*
*  IEEE    (input) LOGICAL
*          A logical flag specifying whether or not the arithmetic
*          system is thought to comply with the IEEE standard.
*
*  EMAX    (output) INTEGER
*          The largest exponent before overflow
*
*  RMAX    (output) DOUBLE PRECISION
*          The largest machine floating-point number.
*
* =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D0, ONE = 1.0D0 )
*     ..
*     .. Local Scalars ..
      INTEGER            EXBITS, EXPSUM, I, LEXP, NBITS, TRY, UEXP
      DOUBLE PRECISION   OLDY, RECBAS, Y, Z
*     ..
*     .. External Functions ..
      DOUBLE PRECISION   DLAMC3
      EXTERNAL           DLAMC3
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MOD
*     ..
*     .. Executable Statements ..
*
*     First compute LEXP and UEXP, two powers of 2 that bound
*     abs(EMIN). We then assume that EMAX + abs(EMIN) will sum
*     approximately to the bound that is closest to abs(EMIN).
*     (EMAX is the exponent of the required number RMAX).
*
      LEXP = 1
      EXBITS = 1
   10 CONTINUE
      TRY = LEXP*2
      IF( TRY.LE.( -EMIN ) ) THEN
         LEXP = TRY
         EXBITS = EXBITS + 1
         GO TO 10
      END IF
      IF( LEXP.EQ.-EMIN ) THEN
         UEXP = LEXP
      ELSE
         UEXP = TRY
         EXBITS = EXBITS + 1
      END IF
*
*     Now -LEXP is less than or equal to EMIN, and -UEXP is greater
*     than or equal to EMIN. EXBITS is the number of bits needed to
*     store the exponent.
*
      IF( ( UEXP+EMIN ).GT.( -LEXP-EMIN ) ) THEN
         EXPSUM = 2*LEXP
      ELSE
         EXPSUM = 2*UEXP
      END IF
*
*     EXPSUM is the exponent range, approximately equal to
*     EMAX - EMIN + 1 .
*
      EMAX = EXPSUM + EMIN - 1
      NBITS = 1 + EXBITS + P
*
*     NBITS is the total number of bits needed to store a
*     floating-point number.
*
      IF( ( MOD( NBITS, 2 ).EQ.1 ) .AND. ( BETA.EQ.2 ) ) THEN
*
*        Either there are an odd number of bits used to store a
*        floating-point number, which is unlikely, or some bits are
*        not used in the representation of numbers, which is possible,
*        (e.g. Cray machines) or the mantissa has an implicit bit,
*        (e.g. IEEE machines, Dec Vax machines), which is perhaps the
*        most likely. We have to assume the last alternative.
*        If this is true, then we need to reduce EMAX by one because
*        there must be some way of representing zero in an implicit-bit
*        system. On machines like Cray, we are reducing EMAX by one
*        unnecessarily.
*
         EMAX = EMAX - 1
      END IF
*
      IF( IEEE ) THEN
*
*        Assume we are on an IEEE machine which reserves one exponent
*        for infinity and NaN.
*
         EMAX = EMAX - 1
      END IF
*
*     Now create RMAX, the largest machine number, which should
*     be equal to (1.0 - BETA**(-P)) * BETA**EMAX .
*
*     First compute 1.0 - BETA**(-P), being careful that the
*     result is less than 1.0 .
*
      RECBAS = ONE / BETA
      Z = BETA - ONE
      Y = ZERO
      DO 20 I = 1, P
         Z = Z*RECBAS
         IF( Y.LT.ONE )
     $      OLDY = Y
         Y = DLAMC3( Y, Z )
   20 CONTINUE
      IF( Y.GE.ONE )
     $   Y = OLDY
*
*     Now multiply by BETA**EMAX to get RMAX.
*
      DO 30 I = 1, EMAX
         Y = DLAMC3( Y*BETA, ZERO )
   30 CONTINUE
*
      RMAX = Y
      RETURN
*
*     End of DLAMC5
*
      END
      SUBROUTINE DLATBS( UPLO, TRANS, DIAG, NORMIN, N, KD, AB, LDAB, X,
     $                   SCALE, CNORM, INFO )
*
*  -- LAPACK auxiliary routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     June 30, 1992
*
*     .. Scalar Arguments ..
      CHARACTER          DIAG, NORMIN, TRANS, UPLO
      INTEGER            INFO, KD, LDAB, N
      DOUBLE PRECISION   SCALE
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   AB( LDAB, * ), CNORM( * ), X( * )
*     ..
*
*  Purpose
*  =======
*
*  DLATBS solves one of the triangular systems
*
*     A *x = s*b  or  A'*x = s*b
*
*  with scaling to prevent overflow, where A is an upper or lower
*  triangular band matrix.  Here A' denotes the transpose of A, x and b
*  are n-element vectors, and s is a scaling factor, usually less than
*  or equal to 1, chosen so that the components of x will be less than
*  the overflow threshold.  If the unscaled problem will not cause
*  overflow, the Level 2 BLAS routine DTBSV is called.  If the matrix A
*  is singular (A(j,j) = 0 for some j), then s is set to 0 and a
*  non-trivial solution to A*x = 0 is returned.
*
*  Arguments
*  =========
*
*  UPLO    (input) CHARACTER*1
*          Specifies whether the matrix A is upper or lower triangular.
*          = 'U':  Upper triangular
*          = 'L':  Lower triangular
*
*  TRANS   (input) CHARACTER*1
*          Specifies the operation applied to A.
*          = 'N':  Solve A * x = s*b  (No transpose)
*          = 'T':  Solve A'* x = s*b  (Transpose)
*          = 'C':  Solve A'* x = s*b  (Conjugate transpose = Transpose)
*
*  DIAG    (input) CHARACTER*1
*          Specifies whether or not the matrix A is unit triangular.
*          = 'N':  Non-unit triangular
*          = 'U':  Unit triangular
*
*  NORMIN  (input) CHARACTER*1
*          Specifies whether CNORM has been set or not.
*          = 'Y':  CNORM contains the column norms on entry
*          = 'N':  CNORM is not set on entry.  On exit, the norms will
*                  be computed and stored in CNORM.
*
*  N       (input) INTEGER
*          The order of the matrix A.  N >= 0.
*
*  KD      (input) INTEGER
*          The number of subdiagonals or superdiagonals in the
*          triangular matrix A.  KD >= 0.
*
*  AB      (input) DOUBLE PRECISION array, dimension (LDAB,N)
*          The upper or lower triangular band matrix A, stored in the
*          first KD+1 rows of the array. The j-th column of A is stored
*          in the j-th column of the array AB as follows:
*          if UPLO = 'U', AB(kd+1+i-j,j) = A(i,j) for max(1,j-kd)<=i<=j;
*          if UPLO = 'L', AB(1+i-j,j)    = A(i,j) for j<=i<=min(n,j+kd).
*
*  LDAB    (input) INTEGER
*          The leading dimension of the array AB.  LDAB >= KD+1.
*
*  X       (input/output) DOUBLE PRECISION array, dimension (N)
*          On entry, the right hand side b of the triangular system.
*          On exit, X is overwritten by the solution vector x.
*
*  SCALE   (output) DOUBLE PRECISION
*          The scaling factor s for the triangular system
*             A * x = s*b  or  A'* x = s*b.
*          If SCALE = 0, the matrix A is singular or badly scaled, and
*          the vector x is an exact or approximate solution to A*x = 0.
*
*  CNORM   (input or output) DOUBLE PRECISION array, dimension (N)
*
*          If NORMIN = 'Y', CNORM is an input argument and CNORM(j)
*          contains the norm of the off-diagonal part of the j-th column
*          of A.  If TRANS = 'N', CNORM(j) must be greater than or equal
*          to the infinity-norm, and if TRANS = 'T' or 'C', CNORM(j)
*          must be greater than or equal to the 1-norm.
*
*          If NORMIN = 'N', CNORM is an output argument and CNORM(j)
*          returns the 1-norm of the offdiagonal part of the j-th column
*          of A.
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -k, the k-th argument had an illegal value
*
*  Further Details
*  ======= =======
*
*  A rough bound on x is computed; if that is less than overflow, DTBSV
*  is called, otherwise, specific code is used which checks for possible
*  overflow or divide-by-zero at every operation.
*
*  A columnwise scheme is used for solving A*x = b.  The basic algorithm
*  if A is lower triangular is
*
*       x[1:n] := b[1:n]
*       for j = 1, ..., n
*            x(j) := x(j) / A(j,j)
*            x[j+1:n] := x[j+1:n] - x(j) * A[j+1:n,j]
*       end
*
*  Define bounds on the components of x after j iterations of the loop:
*     M(j) = bound on x[1:j]
*     G(j) = bound on x[j+1:n]
*  Initially, let M(0) = 0 and G(0) = max{x(i), i=1,...,n}.
*
*  Then for iteration j+1 we have
*     M(j+1) <= G(j) / | A(j+1,j+1) |
*     G(j+1) <= G(j) + M(j+1) * | A[j+2:n,j+1] |
*            <= G(j) ( 1 + CNORM(j+1) / | A(j+1,j+1) | )
*
*  where CNORM(j+1) is greater than or equal to the infinity-norm of
*  column j+1 of A, not counting the diagonal.  Hence
*
*     G(j) <= G(0) product ( 1 + CNORM(i) / | A(i,i) | )
*                  1<=i<=j
*  and
*
*     |x(j)| <= ( G(0) / |A(j,j)| ) product ( 1 + CNORM(i) / |A(i,i)| )
*                                   1<=i< j
*
*  Since |x(j)| <= M(j), we use the Level 2 BLAS routine DTBSV if the
*  reciprocal of the largest M(j), j=1,..,n, is larger than
*  max(underflow, 1/overflow).
*
*  The bound on x(j) is also used to determine when a step in the
*  columnwise method can be performed without fear of overflow.  If
*  the computed bound is greater than a large constant, x is scaled to
*  prevent overflow, but if the bound overflows, x is set to 0, x(j) to
*  1, and scale to 0, and a non-trivial solution to A*x = 0 is found.
*
*  Similarly, a row-wise scheme is used to solve A'*x = b.  The basic
*  algorithm for A upper triangular is
*
*       for j = 1, ..., n
*            x(j) := ( b(j) - A[1:j-1,j]' * x[1:j-1] ) / A(j,j)
*       end
*
*  We simultaneously compute two bounds
*       G(j) = bound on ( b(i) - A[1:i-1,i]' * x[1:i-1] ), 1<=i<=j
*       M(j) = bound on x(i), 1<=i<=j
*
*  The initial values are G(0) = 0, M(0) = max{b(i), i=1,..,n}, and we
*  add the constraint G(j) >= G(j-1) and M(j) >= M(j-1) for j >= 1.
*  Then the bound on x(j) is
*
*       M(j) <= M(j-1) * ( 1 + CNORM(j) ) / | A(j,j) |
*
*            <= M(0) * product ( ( 1 + CNORM(i) ) / |A(i,i)| )
*                      1<=i<=j
*
*  and we can safely call DTBSV if 1/M(n) and 1/G(n) are both greater
*  than max(underflow, 1/overflow).
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, HALF, ONE
      PARAMETER          ( ZERO = 0.0D+0, HALF = 0.5D+0, ONE = 1.0D+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            NOTRAN, NOUNIT, UPPER
      INTEGER            I, IMAX, J, JFIRST, JINC, JLAST, JLEN, MAIND
      DOUBLE PRECISION   BIGNUM, GROW, REC, SMLNUM, SUMJ, TJJ, TJJS,
     $                   TMAX, TSCAL, USCAL, XBND, XJ, XMAX
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            IDAMAX
      DOUBLE PRECISION   DASUM, DDOT, DLAMCH
      EXTERNAL           LSAME, IDAMAX, DASUM, DDOT, DLAMCH
*     ..
*     .. External Subroutines ..
      EXTERNAL           DAXPY, DSCAL, DTBSV, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, MIN
*     ..
*     .. Executable Statements ..
*
      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      NOTRAN = LSAME( TRANS, 'N' )
      NOUNIT = LSAME( DIAG, 'N' )
*
*     Test the input parameters.
*
      IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -1
      ELSE IF( .NOT.NOTRAN .AND. .NOT.LSAME( TRANS, 'T' ) .AND. .NOT.
     $         LSAME( TRANS, 'C' ) ) THEN
         INFO = -2
      ELSE IF( .NOT.NOUNIT .AND. .NOT.LSAME( DIAG, 'U' ) ) THEN
         INFO = -3
      ELSE IF( .NOT.LSAME( NORMIN, 'Y' ) .AND. .NOT.
     $         LSAME( NORMIN, 'N' ) ) THEN
         INFO = -4
      ELSE IF( N.LT.0 ) THEN
         INFO = -5
      ELSE IF( KD.LT.0 ) THEN
         INFO = -6
      ELSE IF( LDAB.LT.KD+1 ) THEN
         INFO = -8
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DLATBS', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 )
     $   RETURN
*
*     Determine machine dependent parameters to control overflow.
*
      SMLNUM = DLAMCH( 'Safe minimum' ) / DLAMCH( 'Precision' )
      BIGNUM = ONE / SMLNUM
      SCALE = ONE
*
      IF( LSAME( NORMIN, 'N' ) ) THEN
*
*        Compute the 1-norm of each column, not including the diagonal.
*
         IF( UPPER ) THEN
*
*           A is upper triangular.
*
            DO 10 J = 1, N
               JLEN = MIN( KD, J-1 )
               CNORM( J ) = DASUM( JLEN, AB( KD+1-JLEN, J ), 1 )
   10       CONTINUE
         ELSE
*
*           A is lower triangular.
*
            DO 20 J = 1, N
               JLEN = MIN( KD, N-J )
               IF( JLEN.GT.0 ) THEN
                  CNORM( J ) = DASUM( JLEN, AB( 2, J ), 1 )
               ELSE
                  CNORM( J ) = ZERO
               END IF
   20       CONTINUE
         END IF
      END IF
*
*     Scale the column norms by TSCAL if the maximum element in CNORM is
*     greater than BIGNUM.
*
      IMAX = IDAMAX( N, CNORM, 1 )
      TMAX = CNORM( IMAX )
      IF( TMAX.LE.BIGNUM ) THEN
         TSCAL = ONE
      ELSE
         TSCAL = ONE / ( SMLNUM*TMAX )
         CALL DSCAL( N, TSCAL, CNORM, 1 )
      END IF
*
*     Compute a bound on the computed solution vector to see if the
*     Level 2 BLAS routine DTBSV can be used.
*
      J = IDAMAX( N, X, 1 )
      XMAX = ABS( X( J ) )
      XBND = XMAX
      IF( NOTRAN ) THEN
*
*        Compute the growth in A * x = b.
*
         IF( UPPER ) THEN
            JFIRST = N
            JLAST = 1
            JINC = -1
            MAIND = KD + 1
         ELSE
            JFIRST = 1
            JLAST = N
            JINC = 1
            MAIND = 1
         END IF
*
         IF( TSCAL.NE.ONE ) THEN
            GROW = ZERO
            GO TO 50
         END IF
*
         IF( NOUNIT ) THEN
*
*           A is non-unit triangular.
*
*           Compute GROW = 1/G(j) and XBND = 1/M(j).
*           Initially, G(0) = max{x(i), i=1,...,n}.
*
            GROW = ONE / MAX( XBND, SMLNUM )
            XBND = GROW
            DO 30 J = JFIRST, JLAST, JINC
*
*              Exit the loop if the growth factor is too small.
*
               IF( GROW.LE.SMLNUM )
     $            GO TO 50
*
*              M(j) = G(j-1) / abs(A(j,j))
*
               TJJ = ABS( AB( MAIND, J ) )
               XBND = MIN( XBND, MIN( ONE, TJJ )*GROW )
               IF( TJJ+CNORM( J ).GE.SMLNUM ) THEN
*
*                 G(j) = G(j-1)*( 1 + CNORM(j) / abs(A(j,j)) )
*
                  GROW = GROW*( TJJ / ( TJJ+CNORM( J ) ) )
               ELSE
*
*                 G(j) could overflow, set GROW to 0.
*
                  GROW = ZERO
               END IF
   30       CONTINUE
            GROW = XBND
         ELSE
*
*           A is unit triangular.
*
*           Compute GROW = 1/G(j), where G(0) = max{x(i), i=1,...,n}.
*
            GROW = MIN( ONE, ONE / MAX( XBND, SMLNUM ) )
            DO 40 J = JFIRST, JLAST, JINC
*
*              Exit the loop if the growth factor is too small.
*
               IF( GROW.LE.SMLNUM )
     $            GO TO 50
*
*              G(j) = G(j-1)*( 1 + CNORM(j) )
*
               GROW = GROW*( ONE / ( ONE+CNORM( J ) ) )
   40       CONTINUE
         END IF
   50    CONTINUE
*
      ELSE
*
*        Compute the growth in A' * x = b.
*
         IF( UPPER ) THEN
            JFIRST = 1
            JLAST = N
            JINC = 1
            MAIND = KD + 1
         ELSE
            JFIRST = N
            JLAST = 1
            JINC = -1
            MAIND = 1
         END IF
*
         IF( TSCAL.NE.ONE ) THEN
            GROW = ZERO
            GO TO 80
         END IF
*
         IF( NOUNIT ) THEN
*
*           A is non-unit triangular.
*
*           Compute GROW = 1/G(j) and XBND = 1/M(j).
*           Initially, M(0) = max{x(i), i=1,...,n}.
*
            GROW = ONE / MAX( XBND, SMLNUM )
            XBND = GROW
            DO 60 J = JFIRST, JLAST, JINC
*
*              Exit the loop if the growth factor is too small.
*
               IF( GROW.LE.SMLNUM )
     $            GO TO 80
*
*              G(j) = max( G(j-1), M(j-1)*( 1 + CNORM(j) ) )
*
               XJ = ONE + CNORM( J )
               GROW = MIN( GROW, XBND / XJ )
*
*              M(j) = M(j-1)*( 1 + CNORM(j) ) / abs(A(j,j))
*
               TJJ = ABS( AB( MAIND, J ) )
               IF( XJ.GT.TJJ )
     $            XBND = XBND*( TJJ / XJ )
   60       CONTINUE
            GROW = MIN( GROW, XBND )
         ELSE
*
*           A is unit triangular.
*
*           Compute GROW = 1/G(j), where G(0) = max{x(i), i=1,...,n}.
*
            GROW = MIN( ONE, ONE / MAX( XBND, SMLNUM ) )
            DO 70 J = JFIRST, JLAST, JINC
*
*              Exit the loop if the growth factor is too small.
*
               IF( GROW.LE.SMLNUM )
     $            GO TO 80
*
*              G(j) = ( 1 + CNORM(j) )*G(j-1)
*
               XJ = ONE + CNORM( J )
               GROW = GROW / XJ
   70       CONTINUE
         END IF
   80    CONTINUE
      END IF
*
      IF( ( GROW*TSCAL ).GT.SMLNUM ) THEN
*
*        Use the Level 2 BLAS solve if the reciprocal of the bound on
*        elements of X is not too small.
*
         CALL DTBSV( UPLO, TRANS, DIAG, N, KD, AB, LDAB, X, 1 )
      ELSE
*
*        Use a Level 1 BLAS solve, scaling intermediate results.
*
         IF( XMAX.GT.BIGNUM ) THEN
*
*           Scale X so that its components are less than or equal to
*           BIGNUM in absolute value.
*
            SCALE = BIGNUM / XMAX
            CALL DSCAL( N, SCALE, X, 1 )
            XMAX = BIGNUM
         END IF
*
         IF( NOTRAN ) THEN
*
*           Solve A * x = b
*
            DO 110 J = JFIRST, JLAST, JINC
*
*              Compute x(j) = b(j) / A(j,j), scaling x if necessary.
*
               XJ = ABS( X( J ) )
               IF( NOUNIT ) THEN
                  TJJS = AB( MAIND, J )*TSCAL
               ELSE
                  TJJS = TSCAL
                  IF( TSCAL.EQ.ONE )
     $               GO TO 100
               END IF
               TJJ = ABS( TJJS )
               IF( TJJ.GT.SMLNUM ) THEN
*
*                    abs(A(j,j)) > SMLNUM:
*
                  IF( TJJ.LT.ONE ) THEN
                     IF( XJ.GT.TJJ*BIGNUM ) THEN
*
*                          Scale x by 1/b(j).
*
                        REC = ONE / XJ
                        CALL DSCAL( N, REC, X, 1 )
                        SCALE = SCALE*REC
                        XMAX = XMAX*REC
                     END IF
                  END IF
                  X( J ) = X( J ) / TJJS
                  XJ = ABS( X( J ) )
               ELSE IF( TJJ.GT.ZERO ) THEN
*
*                    0 < abs(A(j,j)) <= SMLNUM:
*
                  IF( XJ.GT.TJJ*BIGNUM ) THEN
*
*                       Scale x by (1/abs(x(j)))*abs(A(j,j))*BIGNUM
*                       to avoid overflow when dividing by A(j,j).
*
                     REC = ( TJJ*BIGNUM ) / XJ
                     IF( CNORM( J ).GT.ONE ) THEN
*
*                          Scale by 1/CNORM(j) to avoid overflow when
*                          multiplying x(j) times column j.
*
                        REC = REC / CNORM( J )
                     END IF
                     CALL DSCAL( N, REC, X, 1 )
                     SCALE = SCALE*REC
                     XMAX = XMAX*REC
                  END IF
                  X( J ) = X( J ) / TJJS
                  XJ = ABS( X( J ) )
               ELSE
*
*                    A(j,j) = 0:  Set x(1:n) = 0, x(j) = 1, and
*                    scale = 0, and compute a solution to A*x = 0.
*
                  DO 90 I = 1, N
                     X( I ) = ZERO
   90             CONTINUE
                  X( J ) = ONE
                  XJ = ONE
                  SCALE = ZERO
                  XMAX = ZERO
               END IF
  100          CONTINUE
*
*              Scale x if necessary to avoid overflow when adding a
*              multiple of column j of A.
*
               IF( XJ.GT.ONE ) THEN
                  REC = ONE / XJ
                  IF( CNORM( J ).GT.( BIGNUM-XMAX )*REC ) THEN
*
*                    Scale x by 1/(2*abs(x(j))).
*
                     REC = REC*HALF
                     CALL DSCAL( N, REC, X, 1 )
                     SCALE = SCALE*REC
                  END IF
               ELSE IF( XJ*CNORM( J ).GT.( BIGNUM-XMAX ) ) THEN
*
*                 Scale x by 1/2.
*
                  CALL DSCAL( N, HALF, X, 1 )
                  SCALE = SCALE*HALF
               END IF
*
               IF( UPPER ) THEN
                  IF( J.GT.1 ) THEN
*
*                    Compute the update
*                       x(max(1,j-kd):j-1) := x(max(1,j-kd):j-1) -
*                                             x(j)* A(max(1,j-kd):j-1,j)
*
                     JLEN = MIN( KD, J-1 )
                     CALL DAXPY( JLEN, -X( J )*TSCAL,
     $                           AB( KD+1-JLEN, J ), 1, X( J-JLEN ), 1 )
                     I = IDAMAX( J-1, X, 1 )
                     XMAX = ABS( X( I ) )
                  END IF
               ELSE IF( J.LT.N ) THEN
*
*                 Compute the update
*                    x(j+1:min(j+kd,n)) := x(j+1:min(j+kd,n)) -
*                                          x(j) * A(j+1:min(j+kd,n),j)
*
                  JLEN = MIN( KD, N-J )
                  IF( JLEN.GT.0 )
     $               CALL DAXPY( JLEN, -X( J )*TSCAL, AB( 2, J ), 1,
     $                           X( J+1 ), 1 )
                  I = J + IDAMAX( N-J, X( J+1 ), 1 )
                  XMAX = ABS( X( I ) )
               END IF
  110       CONTINUE
*
         ELSE
*
*           Solve A' * x = b
*
            DO 160 J = JFIRST, JLAST, JINC
*
*              Compute x(j) = b(j) - sum A(k,j)*x(k).
*                                    k<>j
*
               XJ = ABS( X( J ) )
               USCAL = TSCAL
               REC = ONE / MAX( XMAX, ONE )
               IF( CNORM( J ).GT.( BIGNUM-XJ )*REC ) THEN
*
*                 If x(j) could overflow, scale x by 1/(2*XMAX).
*
                  REC = REC*HALF
                  IF( NOUNIT ) THEN
                     TJJS = AB( MAIND, J )*TSCAL
                  ELSE
                     TJJS = TSCAL
                  END IF
                  TJJ = ABS( TJJS )
                  IF( TJJ.GT.ONE ) THEN
*
*                       Divide by A(j,j) when scaling x if A(j,j) > 1.
*
                     REC = MIN( ONE, REC*TJJ )
                     USCAL = USCAL / TJJS
                  END IF
                  IF( REC.LT.ONE ) THEN
                     CALL DSCAL( N, REC, X, 1 )
                     SCALE = SCALE*REC
                     XMAX = XMAX*REC
                  END IF
               END IF
*
               SUMJ = ZERO
               IF( USCAL.EQ.ONE ) THEN
*
*                 If the scaling needed for A in the dot product is 1,
*                 call DDOT to perform the dot product.
*
                  IF( UPPER ) THEN
                     JLEN = MIN( KD, J-1 )
                     SUMJ = DDOT( JLEN, AB( KD+1-JLEN, J ), 1,
     $                      X( J-JLEN ), 1 )
                  ELSE
                     JLEN = MIN( KD, N-J )
                     IF( JLEN.GT.0 )
     $                  SUMJ = DDOT( JLEN, AB( 2, J ), 1, X( J+1 ), 1 )
                  END IF
               ELSE
*
*                 Otherwise, use in-line code for the dot product.
*
                  IF( UPPER ) THEN
                     JLEN = MIN( KD, J-1 )
                     DO 120 I = 1, JLEN
                        SUMJ = SUMJ + ( AB( KD+I-JLEN, J )*USCAL )*
     $                         X( J-JLEN-1+I )
  120                CONTINUE
                  ELSE
                     JLEN = MIN( KD, N-J )
                     DO 130 I = 1, JLEN
                        SUMJ = SUMJ + ( AB( I+1, J )*USCAL )*X( J+I )
  130                CONTINUE
                  END IF
               END IF
*
               IF( USCAL.EQ.TSCAL ) THEN
*
*                 Compute x(j) := ( x(j) - sumj ) / A(j,j) if 1/A(j,j)
*                 was not used to scale the dotproduct.
*
                  X( J ) = X( J ) - SUMJ
                  XJ = ABS( X( J ) )
                  IF( NOUNIT ) THEN
*
*                    Compute x(j) = x(j) / A(j,j), scaling if necessary.
*
                     TJJS = AB( MAIND, J )*TSCAL
                  ELSE
                     TJJS = TSCAL
                     IF( TSCAL.EQ.ONE )
     $                  GO TO 150
                  END IF
                  TJJ = ABS( TJJS )
                  IF( TJJ.GT.SMLNUM ) THEN
*
*                       abs(A(j,j)) > SMLNUM:
*
                     IF( TJJ.LT.ONE ) THEN
                        IF( XJ.GT.TJJ*BIGNUM ) THEN
*
*                             Scale X by 1/abs(x(j)).
*
                           REC = ONE / XJ
                           CALL DSCAL( N, REC, X, 1 )
                           SCALE = SCALE*REC
                           XMAX = XMAX*REC
                        END IF
                     END IF
                     X( J ) = X( J ) / TJJS
                  ELSE IF( TJJ.GT.ZERO ) THEN
*
*                       0 < abs(A(j,j)) <= SMLNUM:
*
                     IF( XJ.GT.TJJ*BIGNUM ) THEN
*
*                          Scale x by (1/abs(x(j)))*abs(A(j,j))*BIGNUM.
*
                        REC = ( TJJ*BIGNUM ) / XJ
                        CALL DSCAL( N, REC, X, 1 )
                        SCALE = SCALE*REC
                        XMAX = XMAX*REC
                     END IF
                     X( J ) = X( J ) / TJJS
                  ELSE
*
*                       A(j,j) = 0:  Set x(1:n) = 0, x(j) = 1, and
*                       scale = 0, and compute a solution to A'*x = 0.
*
                     DO 140 I = 1, N
                        X( I ) = ZERO
  140                CONTINUE
                     X( J ) = ONE
                     SCALE = ZERO
                     XMAX = ZERO
                  END IF
  150             CONTINUE
               ELSE
*
*                 Compute x(j) := x(j) / A(j,j) - sumj if the dot
*                 product has already been divided by 1/A(j,j).
*
                  X( J ) = X( J ) / TJJS - SUMJ
               END IF
               XMAX = MAX( XMAX, ABS( X( J ) ) )
  160       CONTINUE
         END IF
         SCALE = SCALE / TSCAL
      END IF
*
*     Scale the column norms by 1/TSCAL for return.
*
      IF( TSCAL.NE.ONE ) THEN
         CALL DSCAL( N, ONE / TSCAL, CNORM, 1 )
      END IF
*
      RETURN
*
*     End of DLATBS
*
      END
      SUBROUTINE DRSCL( N, SA, SX, INCX )
*
*  -- LAPACK auxiliary routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     September 30, 1994
*
*     .. Scalar Arguments ..
      INTEGER            INCX, N
      DOUBLE PRECISION   SA
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   SX( * )
*     ..
*
*  Purpose
*  =======
*
*  DRSCL multiplies an n-element real vector x by the real scalar 1/a.
*  This is done without overflow or underflow as long as
*  the final result x/a does not overflow or underflow.
*
*  Arguments
*  =========
*
*  N       (input) INTEGER
*          The number of components of the vector x.
*
*  SA      (input) DOUBLE PRECISION
*          The scalar a which is used to divide each component of x.
*          SA must be >= 0, or the subroutine will divide by zero.
*
*  SX      (input/output) DOUBLE PRECISION array, dimension
*                         (1+(N-1)*abs(INCX))
*          The n-element vector x.
*
*  INCX    (input) INTEGER
*          The increment between successive values of the vector SX.
*          > 0:  SX(1) = X(1) and SX(1+(i-1)*INCX) = x(i),     1< i<= n
*
* =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            DONE
      DOUBLE PRECISION   BIGNUM, CDEN, CDEN1, CNUM, CNUM1, MUL, SMLNUM
*     ..
*     .. External Functions ..
      DOUBLE PRECISION   DLAMCH
      EXTERNAL           DLAMCH
*     ..
*     .. External Subroutines ..
      EXTERNAL           DSCAL
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS
*     ..
*     .. Executable Statements ..
*
*     Quick return if possible
*
      IF( N.LE.0 )
     $   RETURN
*
*     Get machine parameters
*
      SMLNUM = DLAMCH( 'S' )
      BIGNUM = ONE / SMLNUM
      CALL DLABAD( SMLNUM, BIGNUM )
*
*     Initialize the denominator to SA and the numerator to 1.
*
      CDEN = SA
      CNUM = ONE
*
   10 CONTINUE
      CDEN1 = CDEN*SMLNUM
      CNUM1 = CNUM / BIGNUM
      IF( ABS( CDEN1 ).GT.ABS( CNUM ) .AND. CNUM.NE.ZERO ) THEN
*
*        Pre-multiply X by SMLNUM if CDEN is large compared to CNUM.
*
         MUL = SMLNUM
         DONE = .FALSE.
         CDEN = CDEN1
      ELSE IF( ABS( CNUM1 ).GT.ABS( CDEN ) ) THEN
*
*        Pre-multiply X by BIGNUM if CDEN is small compared to CNUM.
*
         MUL = BIGNUM
         DONE = .FALSE.
         CNUM = CNUM1
      ELSE
*
*        Multiply X by CNUM / CDEN and return.
*
         MUL = CNUM / CDEN
         DONE = .TRUE.
      END IF
*
*     Scale the vector X by MUL
*
      CALL DSCAL( N, MUL, SX, INCX )
*
      IF( .NOT.DONE )
     $   GO TO 10
*
      RETURN
*
*     End of DRSCL
*
      END
      INTEGER          FUNCTION IEEECK( ISPEC, ZERO, ONE )
*
*  -- LAPACK auxiliary routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     June 30, 1998
*
*     .. Scalar Arguments ..
      INTEGER            ISPEC
      REAL               ONE, ZERO
*     ..
*
*  Purpose
*  =======
*
*  IEEECK is called from the ILAENV to verify that Infinity and
*  possibly NaN arithmetic is safe (i.e. will not trap).
*
*  Arguments
*  =========
*
*  ISPEC   (input) INTEGER
*          Specifies whether to test just for inifinity arithmetic
*          or whether to test for infinity and NaN arithmetic.
*          = 0: Verify infinity arithmetic only.
*          = 1: Verify infinity and NaN arithmetic.
*
*  ZERO    (input) REAL
*          Must contain the value 0.0
*          This is passed to prevent the compiler from optimizing
*          away this code.
*
*  ONE     (input) REAL
*          Must contain the value 1.0
*          This is passed to prevent the compiler from optimizing
*          away this code.
*
*  RETURN VALUE:  INTEGER
*          = 0:  Arithmetic failed to produce the correct answers
*          = 1:  Arithmetic produced the correct answers
*
*     .. Local Scalars ..
      REAL               NAN1, NAN2, NAN3, NAN4, NAN5, NAN6, NEGINF,
     $                   NEGZRO, NEWZRO, POSINF
*     ..
*     .. Executable Statements ..
      IEEECK = 1
*
      POSINF = ONE / ZERO
      IF( POSINF.LE.ONE ) THEN
         IEEECK = 0
         RETURN
      END IF
*
      NEGINF = -ONE / ZERO
      IF( NEGINF.GE.ZERO ) THEN
         IEEECK = 0
         RETURN
      END IF
*
      NEGZRO = ONE / ( NEGINF+ONE )
      IF( NEGZRO.NE.ZERO ) THEN
         IEEECK = 0
         RETURN
      END IF
*
      NEGINF = ONE / NEGZRO
      IF( NEGINF.GE.ZERO ) THEN
         IEEECK = 0
         RETURN
      END IF
*
      NEWZRO = NEGZRO + ZERO
      IF( NEWZRO.NE.ZERO ) THEN
         IEEECK = 0
         RETURN
      END IF
*
      POSINF = ONE / NEWZRO
      IF( POSINF.LE.ONE ) THEN
         IEEECK = 0
         RETURN
      END IF
*
      NEGINF = NEGINF*POSINF
      IF( NEGINF.GE.ZERO ) THEN
         IEEECK = 0
         RETURN
      END IF
*
      POSINF = POSINF*POSINF
      IF( POSINF.LE.ONE ) THEN
         IEEECK = 0
         RETURN
      END IF
*
*
*
*
*     Return if we were only asked to check infinity arithmetic
*
      IF( ISPEC.EQ.0 )
     $   RETURN
*
      NAN1 = POSINF + NEGINF
*
      NAN2 = POSINF / NEGINF
*
      NAN3 = POSINF / POSINF
*
      NAN4 = POSINF*ZERO
*
      NAN5 = NEGINF*NEGZRO
*
      NAN6 = NAN5*0.0
*
      IF( NAN1.EQ.NAN1 ) THEN
         IEEECK = 0
         RETURN
      END IF
*
      IF( NAN2.EQ.NAN2 ) THEN
         IEEECK = 0
         RETURN
      END IF
*
      IF( NAN3.EQ.NAN3 ) THEN
         IEEECK = 0
         RETURN
      END IF
*
      IF( NAN4.EQ.NAN4 ) THEN
         IEEECK = 0
         RETURN
      END IF
*
      IF( NAN5.EQ.NAN5 ) THEN
         IEEECK = 0
         RETURN
      END IF
*
      IF( NAN6.EQ.NAN6 ) THEN
         IEEECK = 0
         RETURN
      END IF
*
      RETURN
      END
      INTEGER          FUNCTION ILAENV( ISPEC, NAME, OPTS, N1, N2, N3,
     $                 N4 )
*
*  -- LAPACK auxiliary routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     June 30, 1999
*
*     .. Scalar Arguments ..
      CHARACTER*( * )    NAME, OPTS
      INTEGER            ISPEC, N1, N2, N3, N4
*     ..
*
*  Purpose
*  =======
*
*  ILAENV is called from the LAPACK routines to choose problem-dependent
*  parameters for the local environment.  See ISPEC for a description of
*  the parameters.
*
*  This version provides a set of parameters which should give good,
*  but not optimal, performance on many of the currently available
*  computers.  Users are encouraged to modify this subroutine to set
*  the tuning parameters for their particular machine using the option
*  and problem size information in the arguments.
*
*  This routine will not function correctly if it is converted to all
*  lower case.  Converting it to all upper case is allowed.
*
*  Arguments
*  =========
*
*  ISPEC   (input) INTEGER
*          Specifies the parameter to be returned as the value of
*          ILAENV.
*          = 1: the optimal blocksize; if this value is 1, an unblocked
*               algorithm will give the best performance.
*          = 2: the minimum block size for which the block routine
*               should be used; if the usable block size is less than
*               this value, an unblocked routine should be used.
*          = 3: the crossover point (in a block routine, for N less
*               than this value, an unblocked routine should be used)
*          = 4: the number of shifts, used in the nonsymmetric
*               eigenvalue routines
*          = 5: the minimum column dimension for blocking to be used;
*               rectangular blocks must have dimension at least k by m,
*               where k is given by ILAENV(2,...) and m by ILAENV(5,...)
*          = 6: the crossover point for the SVD (when reducing an m by n
*               matrix to bidiagonal form, if max(m,n)/min(m,n) exceeds
*               this value, a QR factorization is used first to reduce
*               the matrix to a triangular form.)
*          = 7: the number of processors
*          = 8: the crossover point for the multishift QR and QZ methods
*               for nonsymmetric eigenvalue problems.
*          = 9: maximum size of the subproblems at the bottom of the
*               computation tree in the divide-and-conquer algorithm
*               (used by xGELSD and xGESDD)
*          =10: ieee NaN arithmetic can be trusted not to trap
*          =11: infinity arithmetic can be trusted not to trap
*
*  NAME    (input) CHARACTER*(*)
*          The name of the calling subroutine, in either upper case or
*          lower case.
*
*  OPTS    (input) CHARACTER*(*)
*          The character options to the subroutine NAME, concatenated
*          into a single character string.  For example, UPLO = 'U',
*          TRANS = 'T', and DIAG = 'N' for a triangular routine would
*          be specified as OPTS = 'UTN'.
*
*  N1      (input) INTEGER
*  N2      (input) INTEGER
*  N3      (input) INTEGER
*  N4      (input) INTEGER
*          Problem dimensions for the subroutine NAME; these may not all
*          be required.
*
* (ILAENV) (output) INTEGER
*          >= 0: the value of the parameter specified by ISPEC
*          < 0:  if ILAENV = -k, the k-th argument had an illegal value.
*
*  Further Details
*  ===============
*
*  The following conventions have been used when calling ILAENV from the
*  LAPACK routines:
*  1)  OPTS is a concatenation of all of the character options to
*      subroutine NAME, in the same order that they appear in the
*      argument list for NAME, even if they are not used in determining
*      the value of the parameter specified by ISPEC.
*  2)  The problem dimensions N1, N2, N3, N4 are specified in the order
*      that they appear in the argument list for NAME.  N1 is used
*      first, N2 second, and so on, and unused problem dimensions are
*      passed a value of -1.
*  3)  The parameter value returned by ILAENV is checked for validity in
*      the calling subroutine.  For example, ILAENV is used to retrieve
*      the optimal blocksize for STRTRI as follows:
*
*      NB = ILAENV( 1, 'STRTRI', UPLO // DIAG, N, -1, -1, -1 )
*      IF( NB.LE.1 ) NB = MAX( 1, N )
*
*  =====================================================================
*
*     .. Local Scalars ..
      LOGICAL            CNAME, SNAME
      CHARACTER*1        C1
      CHARACTER*2        C2, C4
      CHARACTER*3        C3
      CHARACTER*6        SUBNAM
      INTEGER            I, IC, IZ, NB, NBMIN, NX
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          CHAR, ICHAR, INT, MIN, REAL
*     ..
*     .. External Functions ..
      INTEGER            IEEECK
      EXTERNAL           IEEECK
*     ..
*     .. Executable Statements ..
*
      GO TO ( 100, 100, 100, 400, 500, 600, 700, 800, 900, 1000,
     $        1100 ) ISPEC
*
*     Invalid value for ISPEC
*
      ILAENV = -1
      RETURN
*
  100 CONTINUE
*
*     Convert NAME to upper case if the first character is lower case.
*
      ILAENV = 1
      SUBNAM = NAME
      IC = ICHAR( SUBNAM( 1:1 ) )
      IZ = ICHAR( 'Z' )
      IF( IZ.EQ.90 .OR. IZ.EQ.122 ) THEN
*
*        ASCII character set
*
         IF( IC.GE.97 .AND. IC.LE.122 ) THEN
            SUBNAM( 1:1 ) = CHAR( IC-32 )
            DO 10 I = 2, 6
               IC = ICHAR( SUBNAM( I:I ) )
               IF( IC.GE.97 .AND. IC.LE.122 )
     $            SUBNAM( I:I ) = CHAR( IC-32 )
   10       CONTINUE
         END IF
*
      ELSE IF( IZ.EQ.233 .OR. IZ.EQ.169 ) THEN
*
*        EBCDIC character set
*
         IF( ( IC.GE.129 .AND. IC.LE.137 ) .OR.
     $       ( IC.GE.145 .AND. IC.LE.153 ) .OR.
     $       ( IC.GE.162 .AND. IC.LE.169 ) ) THEN
            SUBNAM( 1:1 ) = CHAR( IC+64 )
            DO 20 I = 2, 6
               IC = ICHAR( SUBNAM( I:I ) )
               IF( ( IC.GE.129 .AND. IC.LE.137 ) .OR.
     $             ( IC.GE.145 .AND. IC.LE.153 ) .OR.
     $             ( IC.GE.162 .AND. IC.LE.169 ) )
     $            SUBNAM( I:I ) = CHAR( IC+64 )
   20       CONTINUE
         END IF
*
      ELSE IF( IZ.EQ.218 .OR. IZ.EQ.250 ) THEN
*
*        Prime machines:  ASCII+128
*
         IF( IC.GE.225 .AND. IC.LE.250 ) THEN
            SUBNAM( 1:1 ) = CHAR( IC-32 )
            DO 30 I = 2, 6
               IC = ICHAR( SUBNAM( I:I ) )
               IF( IC.GE.225 .AND. IC.LE.250 )
     $            SUBNAM( I:I ) = CHAR( IC-32 )
   30       CONTINUE
         END IF
      END IF
*
      C1 = SUBNAM( 1:1 )
      SNAME = C1.EQ.'S' .OR. C1.EQ.'D'
      CNAME = C1.EQ.'C' .OR. C1.EQ.'Z'
      IF( .NOT.( CNAME .OR. SNAME ) )
     $   RETURN
      C2 = SUBNAM( 2:3 )
      C3 = SUBNAM( 4:6 )
      C4 = C3( 2:3 )
*
      GO TO ( 110, 200, 300 ) ISPEC
*
  110 CONTINUE
*
*     ISPEC = 1:  block size
*
*     In these examples, separate code is provided for setting NB for
*     real and complex.  We assume that NB will take the same value in
*     single or double precision.
*
      NB = 1
*
      IF( C2.EQ.'GE' ) THEN
         IF( C3.EQ.'TRF' ) THEN
            IF( SNAME ) THEN
               NB = 64
            ELSE
               NB = 64
            END IF
         ELSE IF( C3.EQ.'QRF' .OR. C3.EQ.'RQF' .OR. C3.EQ.'LQF' .OR.
     $            C3.EQ.'QLF' ) THEN
            IF( SNAME ) THEN
               NB = 32
            ELSE
               NB = 32
            END IF
         ELSE IF( C3.EQ.'HRD' ) THEN
            IF( SNAME ) THEN
               NB = 32
            ELSE
               NB = 32
            END IF
         ELSE IF( C3.EQ.'BRD' ) THEN
            IF( SNAME ) THEN
               NB = 32
            ELSE
               NB = 32
            END IF
         ELSE IF( C3.EQ.'TRI' ) THEN
            IF( SNAME ) THEN
               NB = 64
            ELSE
               NB = 64
            END IF
         END IF
      ELSE IF( C2.EQ.'PO' ) THEN
         IF( C3.EQ.'TRF' ) THEN
            IF( SNAME ) THEN
               NB = 64
            ELSE
               NB = 64
            END IF
         END IF
      ELSE IF( C2.EQ.'SY' ) THEN
         IF( C3.EQ.'TRF' ) THEN
            IF( SNAME ) THEN
               NB = 64
            ELSE
               NB = 64
            END IF
         ELSE IF( SNAME .AND. C3.EQ.'TRD' ) THEN
            NB = 32
         ELSE IF( SNAME .AND. C3.EQ.'GST' ) THEN
            NB = 64
         END IF
      ELSE IF( CNAME .AND. C2.EQ.'HE' ) THEN
         IF( C3.EQ.'TRF' ) THEN
            NB = 64
         ELSE IF( C3.EQ.'TRD' ) THEN
            NB = 32
         ELSE IF( C3.EQ.'GST' ) THEN
            NB = 64
         END IF
      ELSE IF( SNAME .AND. C2.EQ.'OR' ) THEN
         IF( C3( 1:1 ).EQ.'G' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR.
     $          C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR.
     $          C4.EQ.'BR' ) THEN
               NB = 32
            END IF
         ELSE IF( C3( 1:1 ).EQ.'M' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR.
     $          C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR.
     $          C4.EQ.'BR' ) THEN
               NB = 32
            END IF
         END IF
      ELSE IF( CNAME .AND. C2.EQ.'UN' ) THEN
         IF( C3( 1:1 ).EQ.'G' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR.
     $          C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR.
     $          C4.EQ.'BR' ) THEN
               NB = 32
            END IF
         ELSE IF( C3( 1:1 ).EQ.'M' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR.
     $          C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR.
     $          C4.EQ.'BR' ) THEN
               NB = 32
            END IF
         END IF
      ELSE IF( C2.EQ.'GB' ) THEN
         IF( C3.EQ.'TRF' ) THEN
            IF( SNAME ) THEN
               IF( N4.LE.64 ) THEN
                  NB = 1
               ELSE
                  NB = 32
               END IF
            ELSE
               IF( N4.LE.64 ) THEN
                  NB = 1
               ELSE
                  NB = 32
               END IF
            END IF
         END IF
      ELSE IF( C2.EQ.'PB' ) THEN
         IF( C3.EQ.'TRF' ) THEN
            IF( SNAME ) THEN
               IF( N2.LE.64 ) THEN
                  NB = 1
               ELSE
                  NB = 32
               END IF
            ELSE
               IF( N2.LE.64 ) THEN
                  NB = 1
               ELSE
                  NB = 32
               END IF
            END IF
         END IF
      ELSE IF( C2.EQ.'TR' ) THEN
         IF( C3.EQ.'TRI' ) THEN
            IF( SNAME ) THEN
               NB = 64
            ELSE
               NB = 64
            END IF
         END IF
      ELSE IF( C2.EQ.'LA' ) THEN
         IF( C3.EQ.'UUM' ) THEN
            IF( SNAME ) THEN
               NB = 64
            ELSE
               NB = 64
            END IF
         END IF
      ELSE IF( SNAME .AND. C2.EQ.'ST' ) THEN
         IF( C3.EQ.'EBZ' ) THEN
            NB = 1
         END IF
      END IF
      ILAENV = NB
      RETURN
*
  200 CONTINUE
*
*     ISPEC = 2:  minimum block size
*
      NBMIN = 2
      IF( C2.EQ.'GE' ) THEN
         IF( C3.EQ.'QRF' .OR. C3.EQ.'RQF' .OR. C3.EQ.'LQF' .OR.
     $       C3.EQ.'QLF' ) THEN
            IF( SNAME ) THEN
               NBMIN = 2
            ELSE
               NBMIN = 2
            END IF
         ELSE IF( C3.EQ.'HRD' ) THEN
            IF( SNAME ) THEN
               NBMIN = 2
            ELSE
               NBMIN = 2
            END IF
         ELSE IF( C3.EQ.'BRD' ) THEN
            IF( SNAME ) THEN
               NBMIN = 2
            ELSE
               NBMIN = 2
            END IF
         ELSE IF( C3.EQ.'TRI' ) THEN
            IF( SNAME ) THEN
               NBMIN = 2
            ELSE
               NBMIN = 2
            END IF
         END IF
      ELSE IF( C2.EQ.'SY' ) THEN
         IF( C3.EQ.'TRF' ) THEN
            IF( SNAME ) THEN
               NBMIN = 8
            ELSE
               NBMIN = 8
            END IF
         ELSE IF( SNAME .AND. C3.EQ.'TRD' ) THEN
            NBMIN = 2
         END IF
      ELSE IF( CNAME .AND. C2.EQ.'HE' ) THEN
         IF( C3.EQ.'TRD' ) THEN
            NBMIN = 2
         END IF
      ELSE IF( SNAME .AND. C2.EQ.'OR' ) THEN
         IF( C3( 1:1 ).EQ.'G' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR.
     $          C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR.
     $          C4.EQ.'BR' ) THEN
               NBMIN = 2
            END IF
         ELSE IF( C3( 1:1 ).EQ.'M' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR.
     $          C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR.
     $          C4.EQ.'BR' ) THEN
               NBMIN = 2
            END IF
         END IF
      ELSE IF( CNAME .AND. C2.EQ.'UN' ) THEN
         IF( C3( 1:1 ).EQ.'G' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR.
     $          C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR.
     $          C4.EQ.'BR' ) THEN
               NBMIN = 2
            END IF
         ELSE IF( C3( 1:1 ).EQ.'M' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR.
     $          C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR.
     $          C4.EQ.'BR' ) THEN
               NBMIN = 2
            END IF
         END IF
      END IF
      ILAENV = NBMIN
      RETURN
*
  300 CONTINUE
*
*     ISPEC = 3:  crossover point
*
      NX = 0
      IF( C2.EQ.'GE' ) THEN
         IF( C3.EQ.'QRF' .OR. C3.EQ.'RQF' .OR. C3.EQ.'LQF' .OR.
     $       C3.EQ.'QLF' ) THEN
            IF( SNAME ) THEN
               NX = 128
            ELSE
               NX = 128
            END IF
         ELSE IF( C3.EQ.'HRD' ) THEN
            IF( SNAME ) THEN
               NX = 128
            ELSE
               NX = 128
            END IF
         ELSE IF( C3.EQ.'BRD' ) THEN
            IF( SNAME ) THEN
               NX = 128
            ELSE
               NX = 128
            END IF
         END IF
      ELSE IF( C2.EQ.'SY' ) THEN
         IF( SNAME .AND. C3.EQ.'TRD' ) THEN
            NX = 32
         END IF
      ELSE IF( CNAME .AND. C2.EQ.'HE' ) THEN
         IF( C3.EQ.'TRD' ) THEN
            NX = 32
         END IF
      ELSE IF( SNAME .AND. C2.EQ.'OR' ) THEN
         IF( C3( 1:1 ).EQ.'G' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR.
     $          C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR.
     $          C4.EQ.'BR' ) THEN
               NX = 128
            END IF
         END IF
      ELSE IF( CNAME .AND. C2.EQ.'UN' ) THEN
         IF( C3( 1:1 ).EQ.'G' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR.
     $          C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR.
     $          C4.EQ.'BR' ) THEN
               NX = 128
            END IF
         END IF
      END IF
      ILAENV = NX
      RETURN
*
  400 CONTINUE
*
*     ISPEC = 4:  number of shifts (used by xHSEQR)
*
      ILAENV = 6
      RETURN
*
  500 CONTINUE
*
*     ISPEC = 5:  minimum column dimension (not used)
*
      ILAENV = 2
      RETURN
*
  600 CONTINUE 
*
*     ISPEC = 6:  crossover point for SVD (used by xGELSS and xGESVD)
*
      ILAENV = INT( REAL( MIN( N1, N2 ) )*1.6E0 )
      RETURN
*
  700 CONTINUE
*
*     ISPEC = 7:  number of processors (not used)
*
      ILAENV = 1
      RETURN
*
  800 CONTINUE
*
*     ISPEC = 8:  crossover point for multishift (used by xHSEQR)
*
      ILAENV = 50
      RETURN
*
  900 CONTINUE
*
*     ISPEC = 9:  maximum size of the subproblems at the bottom of the
*                 computation tree in the divide-and-conquer algorithm
*                 (used by xGELSD and xGESDD)
*
      ILAENV = 25
      RETURN
*
 1000 CONTINUE
*
*     ISPEC = 10: ieee NaN arithmetic can be trusted not to trap
*
C     ILAENV = 0
      ILAENV = 1
      IF( ILAENV.EQ.1 ) THEN
         ILAENV = IEEECK( 0, 0.0, 1.0 ) 
      END IF
      RETURN
*
 1100 CONTINUE
*
*     ISPEC = 11: infinity arithmetic can be trusted not to trap
*
C     ILAENV = 0
      ILAENV = 1
      IF( ILAENV.EQ.1 ) THEN
         ILAENV = IEEECK( 1, 0.0, 1.0 ) 
      END IF
      RETURN
*
*     End of ILAENV
*
      END
      LOGICAL          FUNCTION LSAME( CA, CB )
*
*  -- LAPACK auxiliary routine (version 2.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     January 31, 1994
*
*     .. Scalar Arguments ..
      CHARACTER          CA, CB
*     ..
*
*  Purpose
*  =======
*
*  LSAME returns .TRUE. if CA is the same letter as CB regardless of
*  case.
*
*  Arguments
*  =========
*
*  CA      (input) CHARACTER*1
*  CB      (input) CHARACTER*1
*          CA and CB specify the single characters to be compared.
*
* =====================================================================
*
*     .. Intrinsic Functions ..
      INTRINSIC          ICHAR
*     ..
*     .. Local Scalars ..
      INTEGER            INTA, INTB, ZCODE
*     ..
*     .. Executable Statements ..
*
*     Test if the characters are equal
*
      LSAME = CA.EQ.CB
      IF( LSAME )
     $   RETURN
*
*     Now test for equivalence if both characters are alphabetic.
*
      ZCODE = ICHAR( 'Z' )
*
*     Use 'Z' rather than 'A' so that ASCII can be detected on Prime
*     machines, on which ICHAR returns a value with bit 8 set.
*     ICHAR('A') on Prime machines returns 193 which is the same as
*     ICHAR('A') on an EBCDIC machine.
*
      INTA = ICHAR( CA )
      INTB = ICHAR( CB )
*
      IF( ZCODE.EQ.90 .OR. ZCODE.EQ.122 ) THEN
*
*        ASCII is assumed - ZCODE is the ASCII code of either lower or
*        upper case 'Z'.
*
         IF( INTA.GE.97 .AND. INTA.LE.122 ) INTA = INTA - 32
         IF( INTB.GE.97 .AND. INTB.LE.122 ) INTB = INTB - 32
*
      ELSE IF( ZCODE.EQ.233 .OR. ZCODE.EQ.169 ) THEN
*
*        EBCDIC is assumed - ZCODE is the EBCDIC code of either lower or
*        upper case 'Z'.
*
         IF( INTA.GE.129 .AND. INTA.LE.137 .OR.
     $       INTA.GE.145 .AND. INTA.LE.153 .OR.
     $       INTA.GE.162 .AND. INTA.LE.169 ) INTA = INTA + 64
         IF( INTB.GE.129 .AND. INTB.LE.137 .OR.
     $       INTB.GE.145 .AND. INTB.LE.153 .OR.
     $       INTB.GE.162 .AND. INTB.LE.169 ) INTB = INTB + 64
*
      ELSE IF( ZCODE.EQ.218 .OR. ZCODE.EQ.250 ) THEN
*
*        ASCII is assumed, on Prime machines - ZCODE is the ASCII code
*        plus 128 of either lower or upper case 'Z'.
*
         IF( INTA.GE.225 .AND. INTA.LE.250 ) INTA = INTA - 32
         IF( INTB.GE.225 .AND. INTB.LE.250 ) INTB = INTB - 32
      END IF
      LSAME = INTA.EQ.INTB
*
*     RETURN
*
*     End of LSAME
*
      END
      SUBROUTINE XERBLA( SRNAME, INFO )
*
*  -- LAPACK auxiliary routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     September 30, 1994
*
*     .. Scalar Arguments ..
      CHARACTER*6        SRNAME
      INTEGER            INFO
*     ..
*
*  Purpose
*  =======
*
*  XERBLA  is an error handler for the LAPACK routines.
*  It is called by an LAPACK routine if an input parameter has an
*  invalid value.  A message is printed and execution stops.
*
*  Installers may consider modifying the STOP statement in order to
*  call system-specific exception-handling facilities.
*
*  Arguments
*  =========
*
*  SRNAME  (input) CHARACTER*6
*          The name of the routine which called XERBLA.
*
*  INFO    (input) INTEGER
*          The position of the invalid parameter in the parameter list
*          of the calling routine.
*
* =====================================================================
*
*     .. Executable Statements ..
*
      WRITE( *, FMT = 9999 )SRNAME, INFO
*
      STOP
*
 9999 FORMAT( ' ** On entry to ', A6, ' parameter number ', I2, ' had ',
     $      'an illegal value' )
*
*     End of XERBLA
*
      END
