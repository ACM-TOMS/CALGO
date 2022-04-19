      SUBROUTINE SORQLX( NB0, NX0,
     >                   M, N, K, A, LDA, TAU, WORK, LWORK, INFO )
C      SUBROUTINE DORGQL( M, N, K, A, LDA, TAU, WORK, LWORK, INFO )
*
*  -- LAPACK routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     June 30, 1999
*
*     .. Scalar Arguments ..
      INTEGER            NB0, NX0, INFO, K, LDA, LWORK, M, N
*     ..
*     .. Array Arguments ..
      REAL               A( LDA, * ), TAU( * ), WORK( * )
C      DOUBLE PRECISION   A( LDA, * ), TAU( * ), WORK( * )
*     ..
*
*  Purpose
*  =======
*
*  sORQLx generates an M-by-N real matrix Q with orthonormal columns,
C*  DORGQL generates an M-by-N real matrix Q with orthonormal columns,
*  which is defined as the last N columns of a product of K elementary
*  reflectors of order M
*
*        Q  =  H(k) . . . H(2) H(1)
*
*  as returned by DGEQLF.
*
*  Arguments
*  =========
*
*  nb0     (input) integer
*          The blocking factor.
*
*  nx0     (input) integer
*          The cross-over point for the nonblocked algorithm.
*
*  M       (input) INTEGER
*          The number of rows of the matrix Q. M >= 0.
*
*  N       (input) INTEGER
*          The number of columns of the matrix Q. M >= N >= 0.
*
*  K       (input) INTEGER
*          The number of elementary reflectors whose product defines the
*          matrix Q. N >= K >= 0.
*
*  A       (input/output) real array, dimension (LDA,N)
C*  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
*          On entry, the (n-k+i)-th column must contain the vector which
*          defines the elementary reflector H(i), for i = 1,2,...,k, as
*          returned by DGEQLF in the last k columns of its array
*          argument A.
*          On exit, the M-by-N matrix Q.
*
*  LDA     (input) INTEGER
*          The first dimension of the array A. LDA >= max(1,M).
*
*  TAU     (input) real array, dimension (K)
C*  TAU     (input) DOUBLE PRECISION array, dimension (K)
*          TAU(i) must contain the scalar factor of the elementary
*          reflector H(i), as returned by DGEQLF.
*
*  WORK    (workspace/output) real array, dimension (LWORK)
C*  WORK    (workspace/output) DOUBLE PRECISION array, dimension (LWORK)
*          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
*
*  LWORK   (input) INTEGER
*          The dimension of the array WORK. LWORK >= max(1,N).
*          For optimum performance LWORK >= N*NB, where NB is the
*          optimal blocksize.
*
*          If LWORK = -1, then a workspace query is assumed; the routine
*          only calculates the optimal size of the WORK array, returns
*          this value as the first entry of the WORK array, and no error
*          message related to LWORK is issued by XERBLA.
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument has an illegal value
*
*  =====================================================================
*
*     .. Parameters ..
      REAL               ZERO
C      DOUBLE PRECISION   ZERO
      PARAMETER          ( ZERO = 0.0E+0 )
C      PARAMETER          ( ZERO = 0.0D+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            LQUERY
      INTEGER            I, IB, IINFO, IWS, J, KK, L, LDWORK, LWKOPT,
     $                   NB, NBMIN, NX
*     ..
*     .. External Subroutines ..
      EXTERNAL           SLARFB, SLARFT, SORG2L, XERBLA
C      EXTERNAL           DLARFB, DLARFT, DORG2L, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
*     ..
*     .. External Functions ..
      INTEGER            ILAENV
      EXTERNAL           ILAENV
*     ..
*     .. Executable Statements ..
*
*     Test the input arguments
*
      INFO = 0
      IF ( NB0 .GT. 0 ) THEN
        NB = NB0
      ELSE
        NB = ILAENV( 1, 'sORGQL', ' ', M, N, K, -1 )
      ENDIF
C      NB = ILAENV( 1, 'DORGQL', ' ', M, N, K, -1 )
      LWKOPT = MAX( 1, N )*NB
      WORK( 1 ) = LWKOPT
      LQUERY = ( LWORK.EQ.-1 )
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 .OR. N.GT.M ) THEN
         INFO = -2
      ELSE IF( K.LT.0 .OR. K.GT.N ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -5
      ELSE IF( LWORK.LT.MAX( 1, N ) .AND. .NOT.LQUERY ) THEN
         INFO = -8
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'sORGQL', -INFO )
C         CALL XERBLA( 'DORGQL', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.LE.0 ) THEN
         WORK( 1 ) = 1
         RETURN
      END IF
*
      NBMIN = 2
      NX = 0
      IWS = N
      IF( NB.GT.1 .AND. NB.LT.K ) THEN
*
*        Determine when to cross over from blocked to unblocked code.
*
         IF ( NX0 .GT. 0 ) THEN
           NX = NX0
         ELSE
           NX = MAX( 0, ILAENV( 3, 'sORGQL', ' ', M, N, K, -1 ) )
         ENDIF
C         NX = MAX( 0, ILAENV( 3, 'DORGQL', ' ', M, N, K, -1 ) )
         IF( NX.LT.K ) THEN
*
*           Determine if workspace is large enough for blocked code.
*
            LDWORK = N
            IWS = LDWORK*NB
            IF( LWORK.LT.IWS ) THEN
*
*              Not enough workspace to use optimal NB:  reduce NB and
*              determine the minimum value of NB.
*
               NB = LWORK / LDWORK
               NBMIN = 2
C               NBMIN = MAX( 2, ILAENV( 2, 'DORGQL', ' ', M, N, K, -1 ) )
            END IF
         END IF
      END IF
*
      IF( NB.GE.NBMIN .AND. NB.LT.K .AND. NX.LT.K ) THEN
*
*        Use blocked code after the first block.
*        The last kk columns are handled by the block method.
*
         KK = MIN( K, ( ( K-NX+NB-1 ) / NB )*NB )
*
*        Set A(m-kk+1:m,1:n-kk) to zero.
*
         DO 20 J = 1, N - KK
            DO 10 I = M - KK + 1, M
               A( I, J ) = ZERO
   10       CONTINUE
   20    CONTINUE
      ELSE
         KK = 0
      END IF
*
*     Use unblocked code for the first or only block.
*
      CALL SORG2L( M-KK, N-KK, K-KK, A, LDA, TAU, WORK, IINFO )
C      CALL DORG2L( M-KK, N-KK, K-KK, A, LDA, TAU, WORK, IINFO )
*
      IF( KK.GT.0 ) THEN
*
*        Use blocked code
*
         DO 50 I = K - KK + 1, K, NB
            IB = MIN( NB, K-I+1 )
            IF( N-K+I.GT.1 ) THEN
*
*              Form the triangular factor of the block reflector
*              H = H(i+ib-1) . . . H(i+1) H(i)
*
               CALL SLARFT( 'Backward', 'Columnwise', M-K+I+IB-1, IB,
     $                      A( 1, N-K+I ), LDA, TAU( I ), WORK, LDWORK )
C               CALL DLARFT( 'Backward', 'Columnwise', M-K+I+IB-1, IB,
C     $                      A( 1, N-K+I ), LDA, TAU( I ), WORK, LDWORK )
*
*              Apply H to A(1:m-k+i+ib-1,1:n-k+i-1) from the left
*
               CALL SLARFB( 'Left', 'No transpose', 'Backward',
     $                      'Columnwise', M-K+I+IB-1, N-K+I-1, IB,
     $                      A( 1, N-K+I ), LDA, WORK, LDWORK, A, LDA,
     $                      WORK( IB+1 ), LDWORK )
C               CALL DLARFB( 'Left', 'No transpose', 'Backward',
C     $                      'Columnwise', M-K+I+IB-1, N-K+I-1, IB,
C     $                      A( 1, N-K+I ), LDA, WORK, LDWORK, A, LDA,
C     $                      WORK( IB+1 ), LDWORK )
            END IF
*
*           Apply H to rows 1:m-k+i+ib-1 of current block
*
            CALL SORG2L( M-K+I+IB-1, IB, IB, A( 1, N-K+I ), LDA,
     $                   TAU( I ), WORK, IINFO )
C            CALL DORG2L( M-K+I+IB-1, IB, IB, A( 1, N-K+I ), LDA,
C     $                   TAU( I ), WORK, IINFO )
*
*           Set rows m-k+i+ib:m of current block to zero
*
            DO 40 J = N - K + I, N - K + I + IB - 1
               DO 30 L = M - K + I + IB, M
                  A( L, J ) = ZERO
   30          CONTINUE
   40       CONTINUE
   50    CONTINUE
      END IF
*
      WORK( 1 ) = IWS
      RETURN
*
*     End of sORQLx
C*     End of DORGQL
*
      END
      SUBROUTINE SORQRX( NB0, NX0,
     >                   M, N, K, A, LDA, TAU, WORK, LWORK, INFO )
C      SUBROUTINE DORGQR( M, N, K, A, LDA, TAU, WORK, LWORK, INFO )
*
*  -- LAPACK routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     June 30, 1999
*
*     .. Scalar Arguments ..
      INTEGER            NB0, NX0, INFO, K, LDA, LWORK, M, N
C      INTEGER            INFO, K, LDA, LWORK, M, N
*     ..
*     .. Array Arguments ..
      REAL               A( LDA, * ), TAU( * ), WORK( * )
C      DOUBLE PRECISION   A( LDA, * ), TAU( * ), WORK( * )
*     ..
*
*  Purpose
*  =======
*
*  sORQRx generates an M-by-N real matrix Q with orthonormal columns,
C*  DORGQR generates an M-by-N real matrix Q with orthonormal columns,
*  which is defined as the first N columns of a product of K elementary
*  reflectors of order M
*
*        Q  =  H(1) H(2) . . . H(k)
*
*  as returned by sGEQRF.
C*  as returned by DGEQRF.
*
*  Arguments
*  =========
*
*  nb0     (input) integer
*          The blocking factor.
*
*  nx0     (input) integer
*          The cross-over point for the non-blocked algorithm.
*
*  M       (input) INTEGER
*          The number of rows of the matrix Q. M >= 0.
*
*  N       (input) INTEGER
*          The number of columns of the matrix Q. M >= N >= 0.
*
*  K       (input) INTEGER
*          The number of elementary reflectors whose product defines the
*          matrix Q. N >= K >= 0.
*
*  A       (input/output) real array, dimension (LDA,N)
C*  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
*          On entry, the i-th column must contain the vector which
*          defines the elementary reflector H(i), for i = 1,2,...,k, as
*          returned by sGEQRF in the first k columns of its array
C*          returned by DGEQRF in the first k columns of its array
*          argument A.
*          On exit, the M-by-N matrix Q.
*
*  LDA     (input) INTEGER
*          The first dimension of the array A. LDA >= max(1,M).
*
*  TAU     (input) real array, dimension (K)
C*  TAU     (input) DOUBLE PRECISION array, dimension (K)
*          TAU(i) must contain the scalar factor of the elementary
*          reflector H(i), as returned by sGEQRF.
C*          reflector H(i), as returned by DGEQRF.
*
*  WORK    (workspace/output) real array, dimension (LWORK)
C*  WORK    (workspace/output) DOUBLE PRECISION array, dimension (LWORK)
*          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
*
*  LWORK   (input) INTEGER
*          The dimension of the array WORK. LWORK >= max(1,N).
*          For optimum performance LWORK >= N*NB, where NB is the
*          optimal blocksize.
*
*          If LWORK = -1, then a workspace query is assumed; the routine
*          only calculates the optimal size of the WORK array, returns
*          this value as the first entry of the WORK array, and no error
*          message related to LWORK is issued by XERBLA.
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument has an illegal value
*
*  =====================================================================
*
*     .. Parameters ..
      REAL               ZERO
C      DOUBLE PRECISION   ZERO
      PARAMETER          ( ZERO = 0.0E+0 )
C      PARAMETER          ( ZERO = 0.0D+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            LQUERY
      INTEGER            I, IB, IINFO, IWS, J, KI, KK, L, LDWORK,
     $                   LWKOPT, NB, NBMIN, NX
*     ..
*     .. External Subroutines ..
      EXTERNAL           SLARFB, SLARFT, SORG2R, XERBLA
C      EXTERNAL           DLARFB, DLARFT, DORG2R, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
*     ..
*     .. External Functions ..
      INTEGER            ILAENV
      EXTERNAL           ILAENV
*     ..
*     .. Executable Statements ..
*
*     Test the input arguments
*
      INFO = 0
      IF ( NB0 .GT. 0 ) THEN
        NB = NB0
      ELSE
        NB = ILAENV( 1, 'sORGQR', ' ', M, N, K, -1 )
      ENDIF
C      NB = ILAENV( 1, 'DORGQR', ' ', M, N, K, -1 )
      LWKOPT = MAX( 1, N )*NB
      WORK( 1 ) = LWKOPT
      LQUERY = ( LWORK.EQ.-1 )
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 .OR. N.GT.M ) THEN
         INFO = -2
      ELSE IF( K.LT.0 .OR. K.GT.N ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -5
      ELSE IF( LWORK.LT.MAX( 1, N ) .AND. .NOT.LQUERY ) THEN
         INFO = -8
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'sORGQR', -INFO )
C         CALL XERBLA( 'DORGQR', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.LE.0 ) THEN
         WORK( 1 ) = 1
         RETURN
      END IF
*
      NBMIN = 2
      NX = 0
      IWS = N
      IF( NB.GT.1 .AND. NB.LT.K ) THEN
*
*        Determine when to cross over from blocked to unblocked code.
*
         IF ( NX0 .GT. 0 ) THEN
           NX = NX0
         ELSE
           NX = MAX( 0, ILAENV( 3, 'sORGQR', ' ', M, N, K, -1 ) )
         ENDIF
C         NX = MAX( 0, ILAENV( 3, 'DORGQR', ' ', M, N, K, -1 ) )
         IF( NX.LT.K ) THEN
*
*           Determine if workspace is large enough for blocked code.
*
            LDWORK = N
            IWS = LDWORK*NB
            IF( LWORK.LT.IWS ) THEN
*
*              Not enough workspace to use optimal NB:  reduce NB and
*              determine the minimum value of NB.
*
               NB = LWORK / LDWORK
               NBMIN = 2
C               NBMIN = MAX( 2, ILAENV( 2, 'DORGQR', ' ', M, N, K, -1 ) )
            END IF
         END IF
      END IF
*
      IF( NB.GE.NBMIN .AND. NB.LT.K .AND. NX.LT.K ) THEN
*
*        Use blocked code after the last block.
*        The first kk columns are handled by the block method.
*
         KI = ( ( K-NX-1 ) / NB )*NB
         KK = MIN( K, KI+NB )
*
*        Set A(1:kk,kk+1:n) to zero.
*
         DO 20 J = KK + 1, N
            DO 10 I = 1, KK
               A( I, J ) = ZERO
   10       CONTINUE
   20    CONTINUE
      ELSE
         KK = 0
      END IF
*
*     Use unblocked code for the last or only block.
*
      IF( KK.LT.N )
     $   CALL SORG2R( M-KK, N-KK, K-KK, A( KK+1, KK+1 ), LDA,
     $                TAU( KK+1 ), WORK, IINFO )
C     $   CALL DORG2R( M-KK, N-KK, K-KK, A( KK+1, KK+1 ), LDA,
C     $                TAU( KK+1 ), WORK, IINFO )
*
      IF( KK.GT.0 ) THEN
*
*        Use blocked code
*
         DO 50 I = KI + 1, 1, -NB
            IB = MIN( NB, K-I+1 )
            IF( I+IB.LE.N ) THEN
*
*              Form the triangular factor of the block reflector
*              H = H(i) H(i+1) . . . H(i+ib-1)
*
               CALL SLARFT( 'Forward', 'Columnwise', M-I+1, IB,
     $                      A( I, I ), LDA, TAU( I ), WORK, LDWORK )
C               CALL DLARFT( 'Forward', 'Columnwise', M-I+1, IB,
C     $                      A( I, I ), LDA, TAU( I ), WORK, LDWORK )
*
*              Apply H to A(i:m,i+ib:n) from the left
*
               CALL SLARFB( 'Left', 'No transpose', 'Forward',
     $                      'Columnwise', M-I+1, N-I-IB+1, IB,
     $                      A( I, I ), LDA, WORK, LDWORK, A( I, I+IB ),
     $                      LDA, WORK( IB+1 ), LDWORK )
C               CALL DLARFB( 'Left', 'No transpose', 'Forward',
C     $                      'Columnwise', M-I+1, N-I-IB+1, IB,
C     $                      A( I, I ), LDA, WORK, LDWORK, A( I, I+IB ),
C     $                      LDA, WORK( IB+1 ), LDWORK )
            END IF
*
*           Apply H to rows i:m of current block
*
            CALL SORG2R( M-I+1, IB, IB, A( I, I ), LDA, TAU( I ), WORK,
     $                   IINFO )
C            CALL DORG2R( M-I+1, IB, IB, A( I, I ), LDA, TAU( I ), WORK,
C     $                   IINFO )
*
*           Set rows 1:i-1 of current block to zero
*
            DO 40 J = I, I + IB - 1
               DO 30 L = 1, I - 1
                  A( L, J ) = ZERO
   30          CONTINUE
   40       CONTINUE
   50    CONTINUE
      END IF
*
      WORK( 1 ) = IWS
      RETURN
*
*     End of sORQRx
C*     End of DORGQR
*
      END
      SUBROUTINE SORGTX( NB0, NX0,
     >                   UPLO, N, A, LDA, TAU, WORK, LWORK, INFO )
C      SUBROUTINE DORGTR( UPLO, N, A, LDA, TAU, WORK, LWORK, INFO )
*
*  -- LAPACK routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     June 30, 1999
*
*     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            NB0, NX0, INFO, LDA, LWORK, N
C      INTEGER            INFO, LDA, LWORK, N
*     ..
*     .. Array Arguments ..
      REAL               A( LDA, * ), TAU( * ), WORK( * )
C      DOUBLE PRECISION   A( LDA, * ), TAU( * ), WORK( * )
*     ..
*
*  Purpose
*  =======
*
*  sORGTx generates a real orthogonal matrix Q which is defined as the
C*  DORGTR generates a real orthogonal matrix Q which is defined as the
*  product of n-1 elementary reflectors of order N, as returned by
*  sSYTRD:
C*  DSYTRD:
*
*  if UPLO = 'U', Q = H(n-1) . . . H(2) H(1),
*
*  if UPLO = 'L', Q = H(1) H(2) . . . H(n-1).
*
*  Arguments
*  =========
*
*  nb0     (input) integer
*          The blocking factor.
*
*  nx0     (input) integer
*          The cross-over point for the nonblocked algorithm.
*
*  UPLO    (input) CHARACTER*1
*          = 'U': Upper triangle of A contains elementary reflectors
*                 from sSYTRD;
C*                 from DSYTRD;
*          = 'L': Lower triangle of A contains elementary reflectors
*                 from sSYTRD.
C*                 from DSYTRD.
*
*  N       (input) INTEGER
*          The order of the matrix Q. N >= 0.
*
*  A       (input/output) sPreicision array, dimension (LDA,N)
C*  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
*          On entry, the vectors which define the elementary reflectors,
*          as returned by sSYTRD.
C*          as returned by DSYTRD.
*          On exit, the N-by-N orthogonal matrix Q.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A. LDA >= max(1,N).
*
*  TAU     (input) sPreicision array, dimension (N-1)
C*  TAU     (input) DOUBLE PRECISION array, dimension (N-1)
*          TAU(i) must contain the scalar factor of the elementary
*          reflector H(i), as returned by sSYTRD.
C*          reflector H(i), as returned by DSYTRD.
*
*  WORK    (workspace/output) real array, dimension (LWORK)
C*  WORK    (workspace/output) DOUBLE PRECISION array, dimension (LWORK)
*          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
*
*  LWORK   (input) INTEGER
*          The dimension of the array WORK. LWORK >= max(1,N-1).
*          For optimum performance LWORK >= (N-1)*NB, where NB is
*          the optimal blocksize.
*
*          If LWORK = -1, then a workspace query is assumed; the routine
*          only calculates the optimal size of the WORK array, returns
*          this value as the first entry of the WORK array, and no error
*          message related to LWORK is issued by XERBLA.
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value
*
*  =====================================================================
*
*     .. Parameters ..
      REAL               ZERO, ONE
C      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0E+0, ONE = 1.0E+0 )
C      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            LQUERY, UPPER
      INTEGER            I, IINFO, J, LWKOPT, NB
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            ILAENV
      EXTERNAL           LSAME, ILAENV
*     ..
*     .. External Subroutines ..
      EXTERNAL           SORQLX, SORQRX, XERBLA
C      EXTERNAL           DORGQL, DORGQR, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX
*     ..
*     .. Executable Statements ..
*
*     Test the input arguments
*
      INFO = 0
      LQUERY = ( LWORK.EQ.-1 )
      UPPER = LSAME( UPLO, 'U' )
      IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -4
      ELSE IF( LWORK.LT.MAX( 1, N-1 ) .AND. .NOT.LQUERY ) THEN
         INFO = -7
      END IF
*
      IF( INFO.EQ.0 ) THEN
         IF( UPPER ) THEN
            IF ( NB0 .GT. 0 ) THEN
              NB = NB0
            ELSE
              NB = ILAENV( 1, 'sORGQL', ' ', N-1, N-1, N-1, -1 )
            ENDIF
C            NB = ILAENV( 1, 'DORGQL', ' ', N-1, N-1, N-1, -1 )
         ELSE
            IF ( NB0 .GT. 0 ) THEN
              NB = NB0
            ELSE
              NB = ILAENV( 1, 'sORGQR', ' ', N-1, N-1, N-1, -1 )
            ENDIF
C            NB = ILAENV( 1, 'DORGQR', ' ', N-1, N-1, N-1, -1 )
         END IF
         LWKOPT = MAX( 1, N-1 )*NB
         WORK( 1 ) = LWKOPT
      END IF
*
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'sORGTR', -INFO )
C         CALL XERBLA( 'DORGTR', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 ) THEN
         WORK( 1 ) = 1
         RETURN
      END IF
*
      IF( UPPER ) THEN
*
*        Q was determined by a call to sSYTRD with UPLO = 'U'
C*        Q was determined by a call to DSYTRD with UPLO = 'U'
*
*        Shift the vectors which define the elementary reflectors one
*        column to the left, and set the last row and column of Q to
*        those of the unit matrix
*
         DO 20 J = 1, N - 1
            DO 10 I = 1, J - 1
               A( I, J ) = A( I, J+1 )
   10       CONTINUE
            A( N, J ) = ZERO
   20    CONTINUE
         DO 30 I = 1, N - 1
            A( I, N ) = ZERO
   30    CONTINUE
         A( N, N ) = ONE
*
*        Generate Q(1:n-1,1:n-1)
*
         CALL SORQLX( NB0, NX0,
     >                N-1, N-1, N-1, A, LDA, TAU, WORK, LWORK, IINFO )
C         CALL DORGQL( N-1, N-1, N-1, A, LDA, TAU, WORK, LWORK, IINFO )
*
      ELSE
*
*        Q was determined by a call to sSYTRD with UPLO = 'L'.
C*        Q was determined by a call to DSYTRD with UPLO = 'L'.
*
*        Shift the vectors which define the elementary reflectors one
*        column to the right, and set the first row and column of Q to
*        those of the unit matrix
*
         DO 50 J = N, 2, -1
            A( 1, J ) = ZERO
            DO 40 I = J + 1, N
               A( I, J ) = A( I, J-1 )
   40       CONTINUE
   50    CONTINUE
         A( 1, 1 ) = ONE
         DO 60 I = 2, N
            A( I, 1 ) = ZERO
   60    CONTINUE
         IF( N.GT.1 ) THEN
*
*           Generate Q(2:n,2:n)
*
            CALL SORQRX( NB0, NX0,
     >                   N-1, N-1, N-1, A( 2, 2 ), LDA, TAU, WORK,
     $                   LWORK, IINFO )
C            CALL DORGQR( N-1, N-1, N-1, A( 2, 2 ), LDA, TAU, WORK,
C     $                   LWORK, IINFO )
         END IF
      END IF
      WORK( 1 ) = LWKOPT
      RETURN
*
*     End of sORGTR
C*     End of DORGTR
*
      END
C      SUBROUTINE DSBTRD( VECT, UPLO, N, KD, AB, LDAB, D, E, Q, LDQ,
      SUBROUTINE SSBTRX( VECT, UPLO, N, KD, AB, LDAB, D, E, Q, LDQ,
     $                   WORK, INFO )
*
*  -- LAPACK routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     June 30, 1999
*
*     .. Scalar Arguments ..
      CHARACTER          UPLO, VECT
      INTEGER            INFO, KD, LDAB, LDQ, N
*     ..
*     .. Array Arguments ..
C      DOUBLE PRECISION   AB( LDAB, * ), D( * ), E( * ), Q( LDQ, * ),
      REAL   AB( LDAB, * ), D( * ), E( * ), Q( LDQ, * ),
     $                   WORK( * )
*     ..
*
*  Purpose
*  =======
*
C*  DSBTRD reduces a real symmetric band matrix A to symmetric
*  sSBTRD reduces a real symmetric band matrix A to symmetric
*  tridiagonal form T by an orthogonal similarity transformation:
*  Q**T * A * Q = T.
*
*  Arguments
*  =========
*
*  VECT    (input) CHARACTER*1
*          = 'N':  do not form Q;
*          = 'V':  form Q;
*          = 'U':  update a matrix X, by forming X*Q.
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
C*  AB      (input/output) DOUBLE PRECISION array, dimension (LDAB,N)
*  AB      (input/output) real array, dimension (LDAB,N)
*          On entry, the upper or lower triangle of the symmetric band
*          matrix A, stored in the first KD+1 rows of the array.  The
*          j-th column of A is stored in the j-th column of the array AB
*          as follows:
*          if UPLO = 'U', AB(kd+1+i-j,j) = A(i,j) for max(1,j-kd)<=i<=j;
*          if UPLO = 'L', AB(1+i-j,j)    = A(i,j) for j<=i<=min(n,j+kd).
*          On exit, the diagonal elements of AB are overwritten by the
*          diagonal elements of the tridiagonal matrix T; if KD > 0, the
*          elements on the first superdiagonal (if UPLO = 'U') or the
*          first subdiagonal (if UPLO = 'L') are overwritten by the
*          off-diagonal elements of T; the rest of AB is overwritten by
*          values generated during the reduction.
*
*  LDAB    (input) INTEGER
*          The leading dimension of the array AB.  LDAB >= KD+1.
*
C*  D       (output) DOUBLE PRECISION array, dimension (N)
*  D       (output) real array, dimension (N)
*          The diagonal elements of the tridiagonal matrix T.
*
C*  E       (output) DOUBLE PRECISION array, dimension (N-1)
*  E       (output) real array, dimension (N-1)
*          The off-diagonal elements of the tridiagonal matrix T:
*          E(i) = T(i,i+1) if UPLO = 'U'; E(i) = T(i+1,i) if UPLO = 'L'.
*
C*  Q       (input/output) DOUBLE PRECISION array, dimension (LDQ,N)
*  Q       (input/output) real array, dimension (LDQ,N)
*          On entry, if VECT = 'U', then Q must contain an N-by-N
*          matrix X; if VECT = 'N' or 'V', then Q need not be set.
*
*          On exit:
*          if VECT = 'V', Q contains the N-by-N orthogonal matrix Q;
*          if VECT = 'U', Q contains the product X*Q;
*          if VECT = 'N', the array Q is not referenced.
*
*  LDQ     (input) INTEGER
*          The leading dimension of the array Q.
*          LDQ >= 1, and LDQ >= N if VECT = 'V' or 'U'.
*
C*  WORK    (workspace) DOUBLE PRECISION array, dimension (N)
*  WORK    (workspace) real array, dimension (N)
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value
*
*  Further Details
*  ===============
*
*  Modified by Linda Kaufman, Bell Labs.
*
*  =====================================================================
*
*     .. Parameters ..
C      DOUBLE PRECISION   ZERO, ONE
      REAL   ZERO, ONE
C      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
      PARAMETER          ( ZERO = 0.0E+0, ONE = 1.0E+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            INITQ, UPPER, WANTQ
      INTEGER            I, I2, IBL, INCA, INCX, IQAEND, IQB, IQEND, J,
     $                   J1, J1END, J1INC, J2, JEND, JIN, JINC, K, KD1,
     $                   KDM1, KDN, L, LAST, LEND, NQ, NR, NRT
C      DOUBLE PRECISION   TEMP
      REAL   TEMP
*     ..
*     .. External Subroutines ..
C      EXTERNAL           DLAR2V, DLARGV, DLARTG, DLARTV, DLASET, DROT,
      EXTERNAL           SLAR2V, SLARGV, SLARTG, SLARTV, SLASET, SROT,
     $                   XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters
*
      INITQ = LSAME( VECT, 'V' )
      WANTQ = INITQ .OR. LSAME( VECT, 'U' )
      UPPER = LSAME( UPLO, 'U' )
      KD1 = KD + 1
      KDM1 = KD - 1
      INCX = LDAB - 1
      IQEND = 1
*
      INFO = 0
      IF( .NOT.WANTQ .AND. .NOT.LSAME( VECT, 'N' ) ) THEN
         INFO = -1
      ELSE IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -2
      ELSE IF( N.LT.0 ) THEN
         INFO = -3
      ELSE IF( KD.LT.0 ) THEN
         INFO = -4
      ELSE IF( LDAB.LT.KD1 ) THEN
         INFO = -6
      ELSE IF( LDQ.LT.MAX( 1, N ) .AND. WANTQ ) THEN
         INFO = -10
      END IF
      IF( INFO.NE.0 ) THEN
C         CALL XERBLA( 'DSBTRD', -INFO )
         CALL XERBLA( 'sSBTRD', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 )
     $   RETURN
*
*     Initialize Q to the unit matrix, if needed
*
      IF( INITQ )
C     $   CALL DLASET( 'Full', N, N, ZERO, ONE, Q, LDQ )
     $   CALL SLASET( 'Full', N, N, ZERO, ONE, Q, LDQ )
*
*     Wherever possible, plane rotations are generated and applied in
*     vector operations of length NR over the index set J1:J2:KD1.
*
*     The cosines and sines of the plane rotations are stored in the
*     arrays D and WORK.
*
      INCA = KD1*LDAB
      KDN = MIN( N-1, KD )
      IF( UPPER ) THEN
*
         IF( KD.GT.1 ) THEN
*
*           Reduce to tridiagonal form, working with upper triangle
*
            NR = 0
            J1 = KDN + 2
            J2 = 1
*
            DO 90 I = 1, N - 2
*
*              Reduce i-th row of matrix to tridiagonal form
*
               DO 80 K = KDN + 1, 2, -1
                  J1 = J1 + KDN
                  J2 = J2 + KDN
*
                  IF( NR.GT.0 ) THEN
*
*                    generate plane rotations to annihilate nonzero
*                    elements which have been created outside the band
*
C                     CALL DLARGV( NR, AB( 1, J1-1 ), INCA, WORK( J1 ),
                     CALL SLARGV( NR, AB( 1, J1-1 ), INCA, WORK( J1 ),
     $                            KD1, D( J1 ), KD1 )
*
*                    apply rotations from the right
*
*
*                    Dependent on the the number of diagonals either
C*                    DLARTV or DROT is used
*                    sLARTV or sROT is used
*
                     IF( NR.GE.2*KD-1 ) THEN
                        DO 10 L = 1, KD - 1
C                           CALL DLARTV( NR, AB( L+1, J1-1 ), INCA,
                           CALL SLARTV( NR, AB( L+1, J1-1 ), INCA,
     $                                  AB( L, J1 ), INCA, D( J1 ),
     $                                  WORK( J1 ), KD1 )
   10                   CONTINUE
*
                     ELSE
                        JEND = J1 + ( NR-1 )*KD1
                        DO 20 JINC = J1, JEND, KD1
C                           CALL DROT( KDM1, AB( 2, JINC-1 ), 1,
                           CALL SROT( KDM1, AB( 2, JINC-1 ), 1,
     $                                AB( 1, JINC ), 1, D( JINC ),
     $                                WORK( JINC ) )
   20                   CONTINUE
                     END IF
                  END IF
*
*
                  IF( K.GT.2 ) THEN
                     IF( K.LE.N-I+1 ) THEN
*
*                       generate plane rotation to annihilate a(i,i+k-1)
*                       within the band
*
C                        CALL DLARTG( AB( KD-K+3, I+K-2 ),
                        CALL SLARTG( AB( KD-K+3, I+K-2 ),
     $                               AB( KD-K+2, I+K-1 ), D( I+K-1 ),
     $                               WORK( I+K-1 ), TEMP )
                        AB( KD-K+3, I+K-2 ) = TEMP
*
*                       apply rotation from the right
*
C                        CALL DROT( K-3, AB( KD-K+4, I+K-2 ), 1,
                        CALL SROT( K-3, AB( KD-K+4, I+K-2 ), 1,
     $                             AB( KD-K+3, I+K-1 ), 1, D( I+K-1 ),
     $                             WORK( I+K-1 ) )
                     END IF
                     NR = NR + 1
                     J1 = J1 - KDN - 1
                  END IF
*
*                 apply plane rotations from both sides to diagonal
*                 blocks
*
                  IF( NR.GT.0 )
C     $               CALL DLAR2V( NR, AB( KD1, J1-1 ), AB( KD1, J1 ),
     $               CALL SLAR2V( NR, AB( KD1, J1-1 ), AB( KD1, J1 ),
     $                            AB( KD, J1 ), INCA, D( J1 ),
     $                            WORK( J1 ), KD1 )
*
*                 apply plane rotations from the left
*
                  IF( NR.GT.0 ) THEN
                     IF( 2*KD-1.LT.NR ) THEN
*
*                    Dependent on the the number of diagonals either
C*                    DLARTV or DROT is used
*                    sLARTV or sROT is used
*
                        DO 30 L = 1, KD - 1
                           IF( J2+L.GT.N ) THEN
                              NRT = NR - 1
                           ELSE
                              NRT = NR
                           END IF
                           IF( NRT.GT.0 )
C     $                        CALL DLARTV( NRT, AB( KD-L, J1+L ), INCA,
     $                        CALL SLARTV( NRT, AB( KD-L, J1+L ), INCA,
     $                                     AB( KD-L+1, J1+L ), INCA,
     $                                     D( J1 ), WORK( J1 ), KD1 )
   30                   CONTINUE
                     ELSE
                        J1END = J1 + KD1*( NR-2 )
                        IF( J1END.GE.J1 ) THEN
                           DO 40 JIN = J1, J1END, KD1
C                              CALL DROT( KD-1, AB( KD-1, JIN+1 ), INCX,
                              CALL SROT( KD-1, AB( KD-1, JIN+1 ), INCX,
     $                                   AB( KD, JIN+1 ), INCX,
     $                                   D( JIN ), WORK( JIN ) )
   40                      CONTINUE
                        END IF
                        LEND = MIN( KDM1, N-J2 )
                        LAST = J1END + KD1
                        IF( LEND.GT.0 )
C     $                     CALL DROT( LEND, AB( KD-1, LAST+1 ), INCX,
     $                     CALL SROT( LEND, AB( KD-1, LAST+1 ), INCX,
     $                                AB( KD, LAST+1 ), INCX, D( LAST ),
     $                                WORK( LAST ) )
                     END IF
                  END IF
*
                  IF( WANTQ ) THEN
*
*                    accumulate product of plane rotations in Q
*
                     IF( INITQ ) THEN
*
*                 take advantage of the fact that Q was
*                 initially the Identity matrix
*
                        IQEND = MAX( IQEND, J2 )
                        I2 = MAX( 0, K-3 )
                        IQAEND = 1 + I*KD
                        IF( K.EQ.2 )
     $                     IQAEND = IQAEND + KD
                        IQAEND = MIN( IQAEND, IQEND )
                        DO 50 J = J1, J2, KD1
                           IBL = I - I2 / KDM1
                           I2 = I2 + 1
                           IQB = MAX( 1, J-IBL )
                           NQ = 1 + IQAEND - IQB
                           IQAEND = MIN( IQAEND+KD, IQEND )
C                           CALL DROT( NQ, Q( IQB, J-1 ), 1, Q( IQB, J ),
                           CALL SROT( NQ, Q( IQB, J-1 ), 1, Q( IQB, J ),
     $                                1, D( J ), WORK( J ) )
   50                   CONTINUE
                     ELSE
*
                        DO 60 J = J1, J2, KD1
C                           CALL DROT( N, Q( 1, J-1 ), 1, Q( 1, J ), 1,
                           CALL SROT( N, Q( 1, J-1 ), 1, Q( 1, J ), 1,
     $                                D( J ), WORK( J ) )
   60                   CONTINUE
                     END IF
*
                  END IF
*
                  IF( J2+KDN.GT.N ) THEN
*
*                    adjust J2 to keep within the bounds of the matrix
*
                     NR = NR - 1
                     J2 = J2 - KDN - 1
                  END IF
*
                  DO 70 J = J1, J2, KD1
*
*                    create nonzero element a(j-1,j+kd) outside the band
*                    and store it in WORK
*
                     WORK( J+KD ) = WORK( J )*AB( 1, J+KD )
                     AB( 1, J+KD ) = D( J )*AB( 1, J+KD )
   70             CONTINUE
   80          CONTINUE
   90       CONTINUE
         END IF
*
         IF( KD.GT.0 ) THEN
*
*           copy off-diagonal elements to E
*
            DO 100 I = 1, N - 1
               E( I ) = AB( KD, I+1 )
  100       CONTINUE
         ELSE
*
*           set E to zero if original matrix was diagonal
*
            DO 110 I = 1, N - 1
               E( I ) = ZERO
  110       CONTINUE
         END IF
*
*        copy diagonal elements to D
*
         DO 120 I = 1, N
            D( I ) = AB( KD1, I )
  120    CONTINUE
*
      ELSE
*
         IF( KD.GT.1 ) THEN
*
*           Reduce to tridiagonal form, working with lower triangle
*
            NR = 0
            J1 = KDN + 2
            J2 = 1
*
            DO 210 I = 1, N - 2
*
*              Reduce i-th column of matrix to tridiagonal form
*
               DO 200 K = KDN + 1, 2, -1
                  J1 = J1 + KDN
                  J2 = J2 + KDN
*
                  IF( NR.GT.0 ) THEN
*
*                    generate plane rotations to annihilate nonzero
*                    elements which have been created outside the band
*
C                     CALL DLARGV( NR, AB( KD1, J1-KD1 ), INCA,
                     CALL SLARGV( NR, AB( KD1, J1-KD1 ), INCA,
     $                            WORK( J1 ), KD1, D( J1 ), KD1 )
*
*                    apply plane rotations from one side
*
*
*                    Dependent on the the number of diagonals either
C*                    DLARTV or DROT is used
*                    sLARTV or sROT is used
*
                     IF( NR.GT.2*KD-1 ) THEN
                        DO 130 L = 1, KD - 1
C                           CALL DLARTV( NR, AB( KD1-L, J1-KD1+L ), INCA,
                           CALL SLARTV( NR, AB( KD1-L, J1-KD1+L ), INCA,
     $                                  AB( KD1-L+1, J1-KD1+L ), INCA,
     $                                  D( J1 ), WORK( J1 ), KD1 )
  130                   CONTINUE
                     ELSE
                        JEND = J1 + KD1*( NR-1 )
                        DO 140 JINC = J1, JEND, KD1
C                           CALL DROT( KDM1, AB( KD, JINC-KD ), INCX,
                           CALL SROT( KDM1, AB( KD, JINC-KD ), INCX,
     $                                AB( KD1, JINC-KD ), INCX,
     $                                D( JINC ), WORK( JINC ) )
  140                   CONTINUE
                     END IF
*
                  END IF
*
                  IF( K.GT.2 ) THEN
                     IF( K.LE.N-I+1 ) THEN
*
*                       generate plane rotation to annihilate a(i+k-1,i)
*                       within the band
*
C                        CALL DLARTG( AB( K-1, I ), AB( K, I ),
                        CALL SLARTG( AB( K-1, I ), AB( K, I ),
     $                               D( I+K-1 ), WORK( I+K-1 ), TEMP )
                        AB( K-1, I ) = TEMP
*
*                       apply rotation from the left
*
C                        CALL DROT( K-3, AB( K-2, I+1 ), LDAB-1,
                        CALL SROT( K-3, AB( K-2, I+1 ), LDAB-1,
     $                             AB( K-1, I+1 ), LDAB-1, D( I+K-1 ),
     $                             WORK( I+K-1 ) )
                     END IF
                     NR = NR + 1
                     J1 = J1 - KDN - 1
                  END IF
*
*                 apply plane rotations from both sides to diagonal
*                 blocks
*
                  IF( NR.GT.0 )
C     $               CALL DLAR2V( NR, AB( 1, J1-1 ), AB( 1, J1 ),
     $               CALL SLAR2V( NR, AB( 1, J1-1 ), AB( 1, J1 ),
     $                            AB( 2, J1-1 ), INCA, D( J1 ),
     $                            WORK( J1 ), KD1 )
*
*                 apply plane rotations from the right
*
*
*                    Dependent on the the number of diagonals either
C*                    DLARTV or DROT is used
*                    sLARTV or sROT is used
*
                  IF( NR.GT.0 ) THEN
                     IF( NR.GT.2*KD-1 ) THEN
                        DO 150 L = 1, KD - 1
                           IF( J2+L.GT.N ) THEN
                              NRT = NR - 1
                           ELSE
                              NRT = NR
                           END IF
                           IF( NRT.GT.0 )
C     $                        CALL DLARTV( NRT, AB( L+2, J1-1 ), INCA,
     $                        CALL SLARTV( NRT, AB( L+2, J1-1 ), INCA,
     $                                     AB( L+1, J1 ), INCA, D( J1 ),
     $                                     WORK( J1 ), KD1 )
  150                   CONTINUE
                     ELSE
                        J1END = J1 + KD1*( NR-2 )
                        IF( J1END.GE.J1 ) THEN
                           DO 160 J1INC = J1, J1END, KD1
C                              CALL DROT( KDM1, AB( 3, J1INC-1 ), 1,
                              CALL SROT( KDM1, AB( 3, J1INC-1 ), 1,
     $                                   AB( 2, J1INC ), 1, D( J1INC ),
     $                                   WORK( J1INC ) )
  160                      CONTINUE
                        END IF
                        LEND = MIN( KDM1, N-J2 )
                        LAST = J1END + KD1
                        IF( LEND.GT.0 )
C     $                     CALL DROT( LEND, AB( 3, LAST-1 ), 1,
     $                     CALL SROT( LEND, AB( 3, LAST-1 ), 1,
     $                                AB( 2, LAST ), 1, D( LAST ),
     $                                WORK( LAST ) )
                     END IF
                  END IF
*
*
*
                  IF( WANTQ ) THEN
*
*                    accumulate product of plane rotations in Q
*
                     IF( INITQ ) THEN
*
*                 take advantage of the fact that Q was
*                 initially the Identity matrix
*
                        IQEND = MAX( IQEND, J2 )
                        I2 = MAX( 0, K-3 )
                        IQAEND = 1 + I*KD
                        IF( K.EQ.2 )
     $                     IQAEND = IQAEND + KD
                        IQAEND = MIN( IQAEND, IQEND )
                        DO 170 J = J1, J2, KD1
                           IBL = I - I2 / KDM1
                           I2 = I2 + 1
                           IQB = MAX( 1, J-IBL )
                           NQ = 1 + IQAEND - IQB
                           IQAEND = MIN( IQAEND+KD, IQEND )
C                           CALL DROT( NQ, Q( IQB, J-1 ), 1, Q( IQB, J ),
                           CALL SROT( NQ, Q( IQB, J-1 ), 1, Q( IQB, J ),
     $                                1, D( J ), WORK( J ) )
  170                   CONTINUE
                     ELSE
*
                        DO 180 J = J1, J2, KD1
C                           CALL DROT( N, Q( 1, J-1 ), 1, Q( 1, J ), 1,
                           CALL SROT( N, Q( 1, J-1 ), 1, Q( 1, J ), 1,
     $                                D( J ), WORK( J ) )
  180                   CONTINUE
                     END IF
                  END IF
*
                  IF( J2+KDN.GT.N ) THEN
*
*                    adjust J2 to keep within the bounds of the matrix
*
                     NR = NR - 1
                     J2 = J2 - KDN - 1
                  END IF
*
                  DO 190 J = J1, J2, KD1
*
*                    create nonzero element a(j+kd,j-1) outside the
*                    band and store it in WORK
*
                     WORK( J+KD ) = WORK( J )*AB( KD1, J )
                     AB( KD1, J ) = D( J )*AB( KD1, J )
  190             CONTINUE
  200          CONTINUE
  210       CONTINUE
         END IF
*
         IF( KD.GT.0 ) THEN
*
*           copy off-diagonal elements to E
*
            DO 220 I = 1, N - 1
               E( I ) = AB( 2, I )
  220       CONTINUE
         ELSE
*
*           set E to zero if original matrix was diagonal
*
            DO 230 I = 1, N - 1
               E( I ) = ZERO
  230       CONTINUE
         END IF
*
*        copy diagonal elements to D
*
         DO 240 I = 1, N
            D( I ) = AB( 1, I )
  240    CONTINUE
      END IF
*
      RETURN
*
C*     End of DSBTRD
*     End of sSBTRD
*
      END
      SUBROUTINE SSYTRX( NB0, NX0,
     >                   UPLO, N, A, LDA, D, E, TAU, WORK, LWORK, INFO )
C      SUBROUTINE DSYTRD( UPLO, N, A, LDA, D, E, TAU, WORK, LWORK, INFO )
*
*  -- LAPACK routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     June 30, 1999
*
*     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            NB0, NX0, INFO, LDA, LWORK, N
C      INTEGER            INFO, LDA, LWORK, N
*     ..
*     .. Array Arguments ..
      REAL               A( LDA, * ), D( * ), E( * ), TAU( * ),
     $                   WORK( * )
C      DOUBLE PRECISION   A( LDA, * ), D( * ), E( * ), TAU( * ),
C     $                   WORK( * )
*     ..
*
*  Purpose
*  =======
*
*  sSYTRx reduces a real symmetric matrix A to real symmetric
C*  DSYTRD reduces a real symmetric matrix A to real symmetric
*  tridiagonal form T by an orthogonal similarity transformation:
*  Q**T * A * Q = T.
*
*  Arguments
*  =========
*
*  nb0     (input) integer
*          The blocking factor.
*
*  nx0     (input) integer
*          The cross-over point to the nonblocked algorithm.
C
*
*  UPLO    (input) CHARACTER*1
*          = 'U':  Upper triangle of A is stored;
*          = 'L':  Lower triangle of A is stored.
*
*  N       (input) INTEGER
*          The order of the matrix A.  N >= 0.
*
*  A       (input/output) real array, dimension (LDA,N)
C*  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
*          On entry, the symmetric matrix A.  If UPLO = 'U', the leading
*          N-by-N upper triangular part of A contains the upper
*          triangular part of the matrix A, and the strictly lower
*          triangular part of A is not referenced.  If UPLO = 'L', the
*          leading N-by-N lower triangular part of A contains the lower
*          triangular part of the matrix A, and the strictly upper
*          triangular part of A is not referenced.
*          On exit, if UPLO = 'U', the diagonal and first superdiagonal
*          of A are overwritten by the corresponding elements of the
*          tridiagonal matrix T, and the elements above the first
*          superdiagonal, with the array TAU, represent the orthogonal
*          matrix Q as a product of elementary reflectors; if UPLO
*          = 'L', the diagonal and first subdiagonal of A are over-
*          written by the corresponding elements of the tridiagonal
*          matrix T, and the elements below the first subdiagonal, with
*          the array TAU, represent the orthogonal matrix Q as a product
*          of elementary reflectors. See Further Details.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,N).
*
*  D       (output) real array, dimension (N)
C*  D       (output) DOUBLE PRECISION array, dimension (N)
*          The diagonal elements of the tridiagonal matrix T:
*          D(i) = A(i,i).
*
*  E       (output) real array, dimension (N-1)
C*  E       (output) DOUBLE PRECISION array, dimension (N-1)
*          The off-diagonal elements of the tridiagonal matrix T:
*          E(i) = A(i,i+1) if UPLO = 'U', E(i) = A(i+1,i) if UPLO = 'L'.
*
*  TAU     (output) real array, dimension (N-1)
C*  TAU     (output) DOUBLE PRECISION array, dimension (N-1)
*          The scalar factors of the elementary reflectors (see Further
*          Details).
*
*  WORK    (workspace/output) real array, dimension (LWORK)
C*  WORK    (workspace/output) DOUBLE PRECISION array, dimension (LWORK)
*          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
*
*  LWORK   (input) INTEGER
*          The dimension of the array WORK.  LWORK >= 1.
*          For optimum performance LWORK >= N*NB, where NB is the
*          optimal blocksize.
*
*          If LWORK = -1, then a workspace query is assumed; the routine
*          only calculates the optimal size of the WORK array, returns
*          this value as the first entry of the WORK array, and no error
*          message related to LWORK is issued by XERBLA.
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value
*
*  Further Details
*  ===============
*
*  If UPLO = 'U', the matrix Q is represented as a product of elementary
*  reflectors
*
*     Q = H(n-1) . . . H(2) H(1).
*
*  Each H(i) has the form
*
*     H(i) = I - tau * v * v'
*
*  where tau is a real scalar, and v is a real vector with
*  v(i+1:n) = 0 and v(i) = 1; v(1:i-1) is stored on exit in
*  A(1:i-1,i+1), and tau in TAU(i).
*
*  If UPLO = 'L', the matrix Q is represented as a product of elementary
*  reflectors
*
*     Q = H(1) H(2) . . . H(n-1).
*
*  Each H(i) has the form
*
*     H(i) = I - tau * v * v'
*
*  where tau is a real scalar, and v is a real vector with
*  v(1:i) = 0 and v(i+1) = 1; v(i+2:n) is stored on exit in A(i+2:n,i),
*  and tau in TAU(i).
*
*  The contents of A on exit are illustrated by the following examples
*  with n = 5:
*
*  if UPLO = 'U':                       if UPLO = 'L':
*
*    (  d   e   v2  v3  v4 )              (  d                  )
*    (      d   e   v3  v4 )              (  e   d              )
*    (          d   e   v4 )              (  v1  e   d          )
*    (              d   e  )              (  v1  v2  e   d      )
*    (                  d  )              (  v1  v2  v3  e   d  )
*
*  where d and e denote diagonal and off-diagonal elements of T, and vi
*  denotes an element of the vector defining H(i).
*
*  =====================================================================
*
*     .. Parameters ..
      REAL                ONE
C      DOUBLE PRECISION   ONE
      PARAMETER          ( ONE = 1.0E+0 )
C      PARAMETER          ( ONE = 1.0D+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            LQUERY, UPPER
      INTEGER            I, IINFO, IWS, J, KK, LDWORK, LWKOPT, NB,
     $                   NBMIN, NX
*     ..
*     .. External Subroutines ..
      EXTERNAL           SLATRD, SSYR2K, SSYTD2, XERBLA
C      EXTERNAL           DLATRD, DSYR2K, DSYTD2, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            ILAENV
      EXTERNAL           LSAME, ILAENV
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters
*
      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      LQUERY = ( LWORK.EQ.-1 )
      IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -4
      ELSE IF( LWORK.LT.1 .AND. .NOT.LQUERY ) THEN
         INFO = -9
      END IF
*
      IF( INFO.EQ.0 ) THEN
*
*        Determine the block size.
*
         IF ( NB0 .GT. 0 ) THEN
           NB = NB0
         ELSE
           NB = ILAENV( 1, 'DSYTRD', UPLO, N, -1, -1, -1 )
         ENDIF
C         NB = ILAENV( 1, 'DSYTRD', UPLO, N, -1, -1, -1 )
         LWKOPT = N*NB
         WORK( 1 ) = LWKOPT
      END IF
*
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'sSYTRD', -INFO )
C         CALL XERBLA( 'DSYTRD', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 ) THEN
         WORK( 1 ) = 1
         RETURN
      END IF
*
      NX = N
      IWS = 1
      IF( NB.GT.1 .AND. NB.LT.N ) THEN
*
*        Determine when to cross over from blocked to unblocked code
*        (last block is always handled by unblocked code).
*
         IF ( NX0 .GT. 0 ) THEN
           NX = MAX( NB, NX0 )
         ELSE
           NX = MAX( NB, ILAENV( 3, 'DSYTRD', UPLO, N, -1, -1, -1 ) )
         ENDIF
C         NX = MAX( NB, ILAENV( 3, 'DSYTRD', UPLO, N, -1, -1, -1 ) )
         IF( NX.LT.N ) THEN
*
*           Determine if workspace is large enough for blocked code.
*
            LDWORK = N
            IWS = LDWORK*NB
            IF( LWORK.LT.IWS ) THEN
*
*              Not enough workspace to use optimal NB:  determine the
*              minimum value of NB, and reduce NB or force use of
*              unblocked code by setting NX = N.
*
               NB = MAX( LWORK / LDWORK, 1 )
               NBMIN = 2
C               NBMIN = ILAENV( 2, 'DSYTRD', UPLO, N, -1, -1, -1 )
               IF( NB.LT.NBMIN )
     $            NX = N
            END IF
         ELSE
            NX = N
         END IF
      ELSE
         NB = 1
      END IF
*
      IF( UPPER ) THEN
*
*        Reduce the upper triangle of A.
*        Columns 1:kk are handled by the unblocked method.
*
         KK = N - ( ( N-NX+NB-1 ) / NB )*NB
         DO 20 I = N - NB + 1, KK + 1, -NB
*
*           Reduce columns i:i+nb-1 to tridiagonal form and form the
*           matrix W which is needed to update the unreduced part of
*           the matrix
*
            CALL SLATRD( UPLO, I+NB-1, NB, A, LDA, E, TAU, WORK,
     $                   LDWORK )
C            CALL DLATRD( UPLO, I+NB-1, NB, A, LDA, E, TAU, WORK,
C     $                   LDWORK )
*
*           Update the unreduced submatrix A(1:i-1,1:i-1), using an
*           update of the form:  A := A - V*W' - W*V'
*
            CALL SSYR2K( UPLO, 'No transpose', I-1, NB, -ONE, A( 1, I ),
     $                   LDA, WORK, LDWORK, ONE, A, LDA )
C            CALL DSYR2K( UPLO, 'No transpose', I-1, NB, -ONE, A( 1, I ),
C     $                   LDA, WORK, LDWORK, ONE, A, LDA )
*
*           Copy superdiagonal elements back into A, and diagonal
*           elements into D
*
            DO 10 J = I, I + NB - 1
               A( J-1, J ) = E( J-1 )
               D( J ) = A( J, J )
   10       CONTINUE
   20    CONTINUE
*
*        Use unblocked code to reduce the last or only block
*
         CALL SSYTD2( UPLO, KK, A, LDA, D, E, TAU, IINFO )
C         CALL DSYTD2( UPLO, KK, A, LDA, D, E, TAU, IINFO )
      ELSE
*
*        Reduce the lower triangle of A
*
         DO 40 I = 1, N - NX, NB
*
*           Reduce columns i:i+nb-1 to tridiagonal form and form the
*           matrix W which is needed to update the unreduced part of
*           the matrix
*
            CALL SLATRD( UPLO, N-I+1, NB, A( I, I ), LDA, E( I ),
     $                   TAU( I ), WORK, LDWORK )
C            CALL DLATRD( UPLO, N-I+1, NB, A( I, I ), LDA, E( I ),
C     $                   TAU( I ), WORK, LDWORK )
*
*           Update the unreduced submatrix A(i+ib:n,i+ib:n), using
*           an update of the form:  A := A - V*W' - W*V'
*
            CALL SSYR2K( UPLO, 'No transpose', N-I-NB+1, NB, -ONE,
     $                   A( I+NB, I ), LDA, WORK( NB+1 ), LDWORK, ONE,
     $                   A( I+NB, I+NB ), LDA )
C            CALL DSYR2K( UPLO, 'No transpose', N-I-NB+1, NB, -ONE,
C     $                   A( I+NB, I ), LDA, WORK( NB+1 ), LDWORK, ONE,
C     $                   A( I+NB, I+NB ), LDA )
*
*           Copy subdiagonal elements back into A, and diagonal
*           elements into D
*
            DO 30 J = I, I + NB - 1
               A( J+1, J ) = E( J )
               D( J ) = A( J, J )
   30       CONTINUE
   40    CONTINUE
*
*        Use unblocked code to reduce the last or only block
*
         CALL SSYTD2( UPLO, N-I+1, A( I, I ), LDA, D( I ), E( I ),
     $                TAU( I ), IINFO )
C         CALL DSYTD2( UPLO, N-I+1, A( I, I ), LDA, D( I ), E( I ),
C     $                TAU( I ), IINFO )
      END IF
*
      WORK( 1 ) = LWKOPT
      RETURN
*
*     End of sSYTRD
C*     End of DSYTRD
*
      END
      SUBROUTINE SLAGSY( N, K, D, A, LDA, ISEED, WORK, INFO )
*
*  -- LAPACK auxiliary test routine (version 3.0)
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     February 29, 1992
*
*     .. Scalar Arguments ..
      INTEGER            INFO, K, LDA, N
*     ..
*     .. Array Arguments ..
      INTEGER            ISEED( 4 )
      REAL               A( LDA, * ), D( * ), WORK( * )
*     ..
*
*  Purpose
*  =======
*
*  SLAGSY generates a real symmetric matrix A, by pre- and post-
*  multiplying a real diagonal matrix D with a random orthogonal matrix:
*  A = U*D*U'. The semi-bandwidth may then be reduced to k by additional
*  orthogonal transformations.
*
*  Arguments
*  =========
*
*  N       (input) INTEGER
*          The order of the matrix A.  N >= 0.
*
*  K       (input) INTEGER
*          The number of nonzero subdiagonals within the band of A.
*          0 <= K <= N-1.
*
*  D       (input) REAL array, dimension (N)
*          The diagonal elements of the diagonal matrix D.
*
*  A       (output) REAL array, dimension (LDA,N)
*          The generated n by n symmetric matrix A (the full matrix is
*          stored).
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= N.
*
*  ISEED   (input/output) INTEGER array, dimension (4)
*          On entry, the seed of the random number generator; the array
*          elements must be between 0 and 4095, and ISEED(4) must be
*          odd.
*          On exit, the seed is updated.
*
*  WORK    (workspace) REAL array, dimension (2*N)
*
*  INFO    (output) INTEGER
*          = 0: successful exit
*          < 0: if INFO = -i, the i-th argument had an illegal value
*
*  =====================================================================
*
*     .. Parameters ..
      REAL               ZERO, ONE, HALF
      PARAMETER          ( ZERO = 0.0E+0, ONE = 1.0E+0, HALF = 0.5E+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            I, J
      REAL               ALPHA, TAU, WA, WB, WN
*     ..
*     .. External Subroutines ..
      EXTERNAL           SAXPY, SGEMV, SGER, SLARNV, SSCAL, SSYMV,
     $                   SSYR2, XERBLA
*     ..
*     .. External Functions ..
      REAL               SDOT, SNRM2
      EXTERNAL           SDOT, SNRM2
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, SIGN
*     ..
*     .. Executable Statements ..
*
*     Test the input arguments
*
      INFO = 0
      IF( N.LT.0 ) THEN
         INFO = -1
      ELSE IF( K.LT.0 .OR. K.GT.N-1 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -5
      END IF
      IF( INFO.LT.0 ) THEN
         CALL XERBLA( 'SLAGSY', -INFO )
         RETURN
      END IF
*
*     initialize lower triangle of A to diagonal matrix
*
      DO 20 J = 1, N
         DO 10 I = J + 1, N
            A( I, J ) = ZERO
   10    CONTINUE
   20 CONTINUE
      DO 30 I = 1, N
         A( I, I ) = D( I )
   30 CONTINUE
*
*     Generate lower triangle of symmetric matrix
*
      DO 40 I = N - 1, 1, -1
*
*        generate random reflection
*
         CALL SLARNV( 3, ISEED, N-I+1, WORK )
         WN = SNRM2( N-I+1, WORK, 1 )
         WA = SIGN( WN, WORK( 1 ) )
         IF( WN.EQ.ZERO ) THEN
            TAU = ZERO
         ELSE
            WB = WORK( 1 ) + WA
            CALL SSCAL( N-I, ONE / WB, WORK( 2 ), 1 )
            WORK( 1 ) = ONE
            TAU = WB / WA
         END IF
*
*        apply random reflection to A(i:n,i:n) from the left
*        and the right
*
*        compute  y := tau * A * u
*
         CALL SSYMV( 'Lower', N-I+1, TAU, A( I, I ), LDA, WORK, 1, ZERO,
     $               WORK( N+1 ), 1 )
*
*        compute  v := y - 1/2 * tau * ( y, u ) * u
*
         ALPHA = -HALF*TAU*SDOT( N-I+1, WORK( N+1 ), 1, WORK, 1 )
         CALL SAXPY( N-I+1, ALPHA, WORK, 1, WORK( N+1 ), 1 )
*
*        apply the transformation as a rank-2 update to A(i:n,i:n)
*
         CALL SSYR2( 'Lower', N-I+1, -ONE, WORK, 1, WORK( N+1 ), 1,
     $               A( I, I ), LDA )
   40 CONTINUE
*
*     Reduce number of subdiagonals to K
*
      DO 60 I = 1, N - 1 - K
*
*        generate reflection to annihilate A(k+i+1:n,i)
*
         WN = SNRM2( N-K-I+1, A( K+I, I ), 1 )
         WA = SIGN( WN, A( K+I, I ) )
         IF( WN.EQ.ZERO ) THEN
            TAU = ZERO
         ELSE
            WB = A( K+I, I ) + WA
            CALL SSCAL( N-K-I, ONE / WB, A( K+I+1, I ), 1 )
            A( K+I, I ) = ONE
            TAU = WB / WA
         END IF
*
*        apply reflection to A(k+i:n,i+1:k+i-1) from the left
*
         CALL SGEMV( 'Transpose', N-K-I+1, K-1, ONE, A( K+I, I+1 ), LDA,
     $               A( K+I, I ), 1, ZERO, WORK, 1 )
         CALL SGER( N-K-I+1, K-1, -TAU, A( K+I, I ), 1, WORK, 1,
     $              A( K+I, I+1 ), LDA )
*
*        apply reflection to A(k+i:n,k+i:n) from the left and the right
*
*        compute  y := tau * A * u
*
         CALL SSYMV( 'Lower', N-K-I+1, TAU, A( K+I, K+I ), LDA,
     $               A( K+I, I ), 1, ZERO, WORK, 1 )
*
*        compute  v := y - 1/2 * tau * ( y, u ) * u
*
         ALPHA = -HALF*TAU*SDOT( N-K-I+1, WORK, 1, A( K+I, I ), 1 )
         CALL SAXPY( N-K-I+1, ALPHA, A( K+I, I ), 1, WORK, 1 )
*
*        apply symmetric rank-2 update to A(k+i:n,k+i:n)
*
         CALL SSYR2( 'Lower', N-K-I+1, -ONE, A( K+I, I ), 1, WORK, 1,
     $               A( K+I, K+I ), LDA )
*
         A( K+I, I ) = -WA
         DO 50 J = K + I + 1, N
            A( J, I ) = ZERO
   50    CONTINUE
   60 CONTINUE
*
*     Store full symmetric matrix
*
      DO 80 J = 1, N
         DO 70 I = J + 1, N
            A( J, I ) = A( I, J )
   70    CONTINUE
   80 CONTINUE
      RETURN
*
*     End of SLAGSY
*
      END
      REAL FUNCTION SLARAN( ISEED )
*
*  -- LAPACK auxiliary routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     February 29, 1992
*
*     .. Array Arguments ..
      INTEGER            ISEED( 4 )
*     ..
*
*  Purpose
*  =======
*
*  SLARAN returns a random real number from a uniform (0,1)
*  distribution.
*
*  Arguments
*  =========
*
*  ISEED   (input/output) INTEGER array, dimension (4)
*          On entry, the seed of the random number generator; the array
*          elements must be between 0 and 4095, and ISEED(4) must be
*          odd.
*          On exit, the seed is updated.
*
*  Further Details
*  ===============
*
*  This routine uses a multiplicative congruential method with modulus
*  2**48 and multiplier 33952834046453 (see G.S.Fishman,
*  'Multiplicative congruential random number generators with modulus
*  2**b: an exhaustive analysis for b = 32 and a partial analysis for
*  b = 48', Math. Comp. 189, pp 331-344, 1990).
*
*  48-bit integers are stored in 4 integer array elements with 12 bits
*  per element. Hence the routine is portable across machines with
*  integers of 32 bits or more.
*
*  =====================================================================
*
*     .. Parameters ..
      INTEGER            M1, M2, M3, M4
      PARAMETER          ( M1 = 494, M2 = 322, M3 = 2508, M4 = 2549 )
      REAL               ONE
      PARAMETER          ( ONE = 1.0E+0 )
      INTEGER            IPW2
      REAL               R
      PARAMETER          ( IPW2 = 4096, R = ONE / IPW2 )
*     ..
*     .. Local Scalars ..
      INTEGER            IT1, IT2, IT3, IT4
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MOD, REAL
*     ..
*     .. Executable Statements ..
*
*     multiply the seed by the multiplier modulo 2**48
*
      IT4 = ISEED( 4 )*M4
      IT3 = IT4 / IPW2
      IT4 = IT4 - IPW2*IT3
      IT3 = IT3 + ISEED( 3 )*M4 + ISEED( 4 )*M3
      IT2 = IT3 / IPW2
      IT3 = IT3 - IPW2*IT2
      IT2 = IT2 + ISEED( 2 )*M4 + ISEED( 3 )*M3 + ISEED( 4 )*M2
      IT1 = IT2 / IPW2
      IT2 = IT2 - IPW2*IT1
      IT1 = IT1 + ISEED( 1 )*M4 + ISEED( 2 )*M3 + ISEED( 3 )*M2 +
     $      ISEED( 4 )*M1
      IT1 = MOD( IT1, IPW2 )
*
*     return updated seed
*
      ISEED( 1 ) = IT1
      ISEED( 2 ) = IT2
      ISEED( 3 ) = IT3
      ISEED( 4 ) = IT4
*
*     convert 48-bit integer to a real number in the interval (0,1)
*
      SLARAN = R*( REAL( IT1 )+R*( REAL( IT2 )+R*( REAL( IT3 )+R*
     $         ( REAL( IT4 ) ) ) ) )
      RETURN
*
*     End of SLARAN
*
      END
* **********************************************************************
*
        PROGRAM SRUN
*
* Description:
*
*   This is the testing and timing driver for the SBR toolbox.
*
* Author: Bruno Lang
*         Aachen University of Technology
*         na.blang@na-net.ornl.gov
* Date: May 17, 2000
* Version: SBR Toolbox, Rev. 1.4.1
*
* Local constants and variables:
*
        REAL                  ONE
        PARAMETER             ( ONE = 1.0E0 )
*
*          --- array dimensions ---
*
        INTEGER               NMAX, SZA, SZACPY, SZU, SZWORK
        PARAMETER             ( NMAX = 1400 ,
     >                          SZA = NMAX * NMAX ,
     >                          SZACPY = 500 * NMAX ,
     >                          SZU = NMAX * NMAX ,
     >                          SZWORK = NMAX * NMAX  )
        REAL                  A( SZA ), ACPY( SZACPY ), U( SZU ),
     >                        EVLS( NMAX ), D( NMAX ), E( NMAX ),
     >                        TAU( NMAX ), WORK( SZWORK )
        COMMON    /ARRAYS/    A, ACPY, U, EVLS, D, E, TAU, WORK
*
*   nmax    maximum matrix size
*   a       holds the matrix A to be reduced (full or banded)
*   sza     size of the array a.
*           sza >= nmax * nmax,       if full matrices are handled
*               >= ( bmax+1 ) * nmax, if only banded matrices with
*                                     semibandwidth <= bmax are handled
*   acpy    holds the reduced matrix (banded)
*   szacpy  size of the array acpy.
*           szacpy >= ( bmax+1 ) * nmax, where bmax is the maximum
*                                        semibandwith
*   u       holds the accumulated transformation matrix U, if required
*   szu     size of the array u.
*           szu >= nmax * nmax, if U is accumulated
*               >= 0,           otherwise
*   evls    the eigenvalues of the matrix A (accessed during
*           initialization)
*   d       main diagonal of the reduced tridiagonal matrix
*   e       sub(super)diagonal of the reduced tridiagonal matrix
*   tau     scaling factors in the reduction full -> banded/tridiagonal
*   work    workspace
*   szwork  size of the array work.
*           szwork depends on the reduction algorithms that are run, on
*           the matrix sizes, on the blocking factors, and on the
*           leading dimensions of the matrices.
*           Recomended: szwork >= 100 * nmax.
*   arrays  common block is used only because some machines cannot
*           handle large arrays on the stack
*
*          --- I/O control ---
*
        INTEGER               IUNIT, OUNIT, OUNITF
	INTEGER               I1MACH
*       PARAMETER             ( IUNIT = 5, OUNITF = 8 )
*       CHARACTER*8           IFILE, OFILE
        CHARACTER*8           OFILE
*       PARAMETER             ( IFILE = 'INFILE' )
        INTEGER               OLEVEL, OUTSUM, OUTTUN, OUTRUN, OUTDET,
     >                        LINENO, LINE0
        PARAMETER             ( OUTSUM = 1, OUTTUN = 2, OUTRUN = 3,
     >                          OUTDET = 4 )
        CHARACTER*78          MBUF
*
*   iunit   unit number for input file
*   ounit   unit number for output
*   ounitf  unit number for output file
*   ifile   name of input file
*   ofile   name of output file
*   olevel  level of output
*           olevel = outsum : only summary for all the runs
*                  = outtun : short summaries for tuning
*                  = outrun : summary for each run
*                  = outdet : detailed output for each run
*   lineno  number of current input line
*   line0   first input line corresponding to the current run
*   mbuf    buffer for messages
*
*          --- for statistics ---
*
        INTEGER               SYTRD, ORGTR, SBTRD, SYRDB, SYGTR, SBRDB,
     >                        SBRDT, SY2BC, SY2BI, SB2BC, SB2BI, SYRDD,
     >                        SBRDD, TOTAL
        PARAMETER             ( SYTRD = 1, ORGTR = 2, SBTRD = 3,
     >                          SYRDB = 4, SYGTR = 5, SBRDB = 6,
     >                          SBRDT = 7, SY2BC = 8, SY2BI = 9,
     >                          SB2BC = 10, SB2BI = 11, SYRDD = 12,
     >                          SBRDD = 13, TOTAL = 14 )
        INTEGER               TRY( TOTAL ), SKIP( TOTAL ), FAIL( TOTAL )
        CHARACTER*6           PRBLMN( TOTAL )
        CHARACTER*40          COMMNT( TOTAL )
*
*   sytrd   tests for ssytrd (full -> tridiagonal, LAPACK)
*   orgtr   tests for sorgtr (backward accumulation, LAPACK)
*   sbtrd   tests for ssbtrd (banded -> tridiagonal, LAPACK)
*   syrdb   tests for ssyrdb (full -> banded, SBR)
*   sygtr   tests for ssygtr (backward accumulation, SBR)
*   sbrdb   tests for ssbrdb (banded -> banded, SBR)
*   sbrdt   tests for ssbrdt (banded -> tridiagonal, SBR)
*   sy2bc   tests for ssy2bc (full -> banded copy, SBR)
*   sy2bi   tests for ssy2bi (full -> banded repacking, SBR)
*   sb2bc   tests for ssb2bc (banded -> banded copy, SBR)
*   sb2bi   tests for ssb2bi (banded -> banded repacking, SBR)
*   syrdd   tests for ssyrdd (reduction driver for full matrices, SBR)
*   sbrdd   tests for ssbrdd (reduction driver for banded matrices, SBR)
*   total   all tests
*   count   counters for all tests
*   fail    counters for failed tests
*   skip    counters for skipped tests
*   prblmn  names of the reduction paths (for output)
*   commnt  comments to the reduction paths (for output)
*
        DOUBLE PRECISION      STIME, BTIME, CTIME, RTIME
*
*   stime   time before all the tests
*   btime   time before the procedure call
*   ctime   time for the procedure call
*   rtime   time for a full reduction/accumulation path
*
*        --- controlling the test runs ---
*
        REAL                  EPS, NOTUSD
*
*   eps     machine epsilon
*   notusd  (large) value for marking memory that should not be accessed
*
        CHARACTER*5           PNAME
        INTEGER               MAXDIM, PRBLM, N, B1, LDA1, B2, LDA2, LDU,
     >                        NSTEPS, STEP, LWORK1, LWORK2, NX1, NX2,
     >                        XINFO1, XINFO2, INFO
        CHARACTER             MTYPE, UPLO, FRMT1, FRMT2, JOBU, JOBU0
        REAL                  DIAM, DRPTOL
        INTEGER               B( 100 ), NB( 100 )
*
*   pname   name of the routine to test (without precision specifier)
*   maxdim  maximum matrix size for the runs (runs with larger matrices
*           are skipped)
*   prblm   number of the routine to test (numbering sytrd=1, etc. as
*           above)
*   mtype   matrix type
*           mtype = 'R' : Eigenvalues are uniformly random in [ -1, 1 ).
*                 = 'I' : ISDA-type matrix: eigenvalues are clustered
*                         around 0 and 1.
*   diam    diameter of the eigenvalue clusters in multiples of the
*           machine precision (used only for ISDA-type matrices)
*   uplo    use upper or lower triangle of the matrix ?
*   n       matrix dimension
*   frmt1   format of the matrix before the reduction ('F'ull/'B'anded)
*   b1      semibandwidth before the reduction
*   lda1    leading dimension before the reductino
*   frmt2   format after the reduction ('B'anded/'T'ridiagonal)
*   b2      semibandwidth after the reduction
*   lda2    leading dimension after the reduction
*   jobu0   accumulate the transformations in an orthogonal matrix ?
*   jobu    accumulate in the current routine ?
*   ldu     leading dimension of the array u
*   nsteps  number of reduction steps (used only in the driver routines)
*   step    number of the last step performed or tried
*   b       intermediate bandwidths for multi-step reduction
*   nb      blocking factors for the reduction/accumulation steps
*   drptol  threshold for skipping ransformations
*   lwork1  length of workspace for the reduction
*   lwork2  length of workspace for backward accumulation
*   nx1     cross-over point for LAPACK tridiagonalization ssytrd
*   nx2     cross-over point for LAPACK backward accumulation sorgtr
*   xinfo1  error code expected from the reduction routine
*   xinfo2  error code expected from backward accumulation
*   info    actual error code returned from the routines
*
        CHARACTER             CHECK
        REAL                  ORTHU, SORTHU, TORTHU, RES, SRES, TRES,
     >                        DEVLS, SDEVLS, TDEVLS
        LOGICAL               CHECKS, ERRCHK, OK
        INTEGER               NUSED
*
*   check   is result checking required ? ('C'check/'N'oCheck)
*   checks  are checks enabled ?
*   orthu   ortogonality error, i.e., (Frobenius) norm of U' * U - I
*   sorthu  = orthu / ( n * eps ), the scaled orthogonality error
*   torthu  threshold on sorthu for accepting the run
*   res     residual, i.e., (Frobenius) norm of U * B - A * U, where
*           A is the original matrix and B is the reduced matrix
*   sres    = res / ( n * eps * norm( A ) ), the scaled residual
*   tres    threshold on sres
*   devls   deviation of the eigenvalues, i.e., (Frobenius) norm of
*           spec( A ) - spec( B ), where spec( A ) and spec( B ) are
*           the (ascendingly sorted) eigenvalues of the original matrix
*           and the reduced matrix, respectively
*   sdevls  = devls / ( n * eps * norm( A ) ), the scaled deviation
*   tdevls  threshold on devls
*   errchk  was it a check for error exit ?
*   ok      did the numerical checks meet the thresholds ?
*   nused   number of workspace elements that were modified
*
      INTEGER                 I
*
* Routines called:
*
	REAL                  SECOND
        REAL                  R1MACH
        LOGICAL               LSAME, LXSAME
        EXTERNAL              SECOND, R1MACH, LSAME, LXSAME
*
*   dsecnd  timer (LAPACK)
*   slamch  determine machine parameters (LAPACK)
*   lsame   case-insensitive character-matching (BLAS)
*   lxsame  case-insensitive string matching (this file)
*
        EXTERNAL              SLACPY, SORGTX, SSB2BC, SSB2BI, SSBRDB,
     >                        SSBRDD, SSBRDT, SSBTRX, SSY2BC, SSY2BI,
     >                        SSYGTR, SSYINI, SSYNCK, SSYRDB, SSYRDD,
     >                        SSYTRX
*
*   slacpy  copy matrix (LAPACK)
*   sorgtx  backward accumulation of transformations (LAPACK, modified)
*   ssb2bc  banded -> banded copy (SBR)
*   ssb2bi  banded -> banded repacking (SBR)
*   ssbrdb  banded -> banded reduction (SBR)
*   ssbrdd  reduction driver for banded matrices (SBR)
*   ssbrdt  banded -> tridiagonal reduction (SBR)
*   ssbtrx  banded -> tridiagonal reduction (LAPACK, modified)
*   ssy2bc  full -> banded copy (SBR)
*   ssy2bi  full -> banded repacking (SBR)
*   ssygtr  backward accumulation of transformations (SBR)
*   ssyini  initialize symmetric (full or banded) matrix
*   ssynck  numerical checks for the reduction
*   ssyrdb  full -> banded reduction (SBR)
*   ssyrdd  reduction driver for full matrices (SBR)
*   ssytrx  full -> tridiagonal reduction (LAPACK, modified)
*
        INTRINSIC             MAX, MIN
*
*          --- names and explanations for the test paths ---
*
        DATA
     >      PRBLMN( SYTRD )/'ssytrd'/,
     >      COMMNT( SYTRD )/'full -> tridiagonal, one-step (LAPACK)'/,
     >      PRBLMN( ORGTR )/'sorgtr'/,
     >      COMMNT( ORGTR )/'accumulation of transformations (LAPACK)'/,
     >      PRBLMN( SBTRD )/'ssbtrd'/,
     >      COMMNT( SBTRD )/'banded -> tridiagonal, one-step (LAPACK)'/
        DATA
     >      PRBLMN( SYRDB )/'ssyrdb'/,
     >      COMMNT( SYRDB )/'full -> banded, one-step (SBR)'/,
     >      PRBLMN( SYGTR )/'ssygtr'/,
     >      COMMNT( SYGTR )/'accumulation of transformations (SBR)'/,
     >      PRBLMN( SBRDB )/'ssbrdb'/,
     >      COMMNT( SBRDB )/'banded -> banded, one-step (SBR)'/,
     >      PRBLMN( SBRDT )/'ssbrdt'/,
     >      COMMNT( SBRDT )/'banded -> tridiagonal, one-step (SBR)'/
        DATA
     >      PRBLMN( SY2BC )/'ssy2bc'/,
     >      COMMNT( SY2BC )/'copy full -> banded (SBR)'/,
     >      PRBLMN( SY2BI )/'ssy2bi'/,
     >      COMMNT( SY2BI )/'in-place copy full -> banded (SBR)'/,
     >      PRBLMN( SB2BC )/'ssb2bc'/,
     >      COMMNT( SB2BC )/'copy banded -> banded (SBR)'/,
     >      PRBLMN( SB2BI )/'ssb2bi'/,
     >      COMMNT( SB2BI )/'in-place copy banded -> banded (SBR)'/
        DATA
     >      PRBLMN( SYRDD )/'ssyrdd'/,
     >      COMMNT( SYRDD )/'reduction driver for full matrices (SBR)'/,
     >      PRBLMN( SBRDD )/'ssbrdd'/,
     >      COMMNT( SBRDD )/'reduction driver for banded matrcs (SBR)'/,
     >      PRBLMN( TOTAL )/'Total '/,
     >      COMMNT( TOTAL )/' '/
*
* ----------------------------------------------------------------------
*
        STIME = SECOND()
*
*          --- initialize ---
*
        DO 1000 I = 1, TOTAL
          TRY( I ) = 0
          SKIP( I ) = 0
          FAIL( I ) = 0
 1000   CONTINUE
*
*       EPS = SLAMCH( 'Epsilon' )
*       NOTUSD = ONE / SLAMCH( 'SafeMinimum' )
        EPS = R1MACH(4)
	NOTUSD = ONE/(R1MACH(1)*(ONE + EPS))
*
*          --- open input and output files ---
*
*       OPEN( UNIT=IUNIT, FILE=IFILE )
        LINENO = 0
*
	IUNIT = I1MACH(1)
        OUNIT = I1MACH(2)
        OLEVEL = OUTSUM
        MAXDIM = NMAX
        CHECKS = .TRUE.
*
 2000   CONTINUE
*
*            --- read next line from input file ---
*
          LINENO = LINENO + 1
          LINE0 = LINENO
          READ( IUNIT, * ) PNAME
*         WRITE( *, * ) LINE0
*
*            --- if comment line then just skip it ---
*
          IF ( LSAME( PNAME( 1:1 ), 'Comment' ) ) THEN
            GOTO 2000
*
*            --- filename for output file ? ---
*
          ELSEIF ( LSAME( PNAME( 1:1 ), 'FileName' ) ) THEN
            LINENO = LINENO + 1
            READ( IUNIT, * ) OFILE( 2:8 )
            OFILE( 1:1 ) = 'SinglePrec'
            OUNIT = OUNITF
            OPEN( UNIT=OUNIT, FILE=OFILE )
            GOTO 2000
*
*            --- set output level ? ---
*
          ELSEIF ( LSAME( PNAME( 1:1 ), 'OutputLevel' ) ) THEN
            LINENO = LINENO + 1
            READ( IUNIT, * ) I
            IF ( ( I .LT. OUTSUM ) .OR. ( I .GT. OUTDET ) ) THEN
              WRITE( OUNIT, 9300 ) LINENO, 'Output level out of range'
            ELSE
              OLEVEL = I
            ENDIF
            GOTO 2000
*
*            --- message to print ---
*
          ELSEIF ( LSAME( PNAME( 1:1 ), 'Print' ) ) THEN
            LINENO = LINENO + 1
            READ( IUNIT, * ) MBUF
            WRITE( OUNIT, 9100 ) MBUF
            GOTO 2000
*
*            --- maximum matrix size ---
*
          ELSEIF ( LSAME( PNAME( 1:1 ), 'MaxDim' ) ) THEN
            LINENO = LINENO + 1
            READ( IUNIT, * ) MAXDIM
            GOTO 2000
*
*            --- disable checks ? ---
*
          ELSEIF ( LSAME( PNAME( 1:1 ), 'NoChecks' ) ) THEN
            CHECKS = .FALSE.
            GOTO 2000
*
*            --- test run ?
*                If so then read and set parameters for the run ---
*
          ELSEIF ( LXSAME( PNAME, 'SYTRD' ) ) THEN
*
*              --- LAPACK reduction full -> tridiagonal ---
*
            PRBLM = SYTRD
            READ( IUNIT, * ) MTYPE, DIAM, UPLO, N, LDA1, NB( 1 ), NX1,
     >                       LWORK1, XINFO1, JOBU, NB( 2 ), NX2, LWORK2,
     >                       XINFO2
            FRMT1 = 'Full'
            B1 = MAX( N-1, 0 )
            FRMT2 = 'Tridiagonal'
            B2 = 1
            LDA2 = B2 + 1
            LDU = N
          ELSEIF ( LXSAME( PNAME, 'SBTRD' ) ) THEN
*
*              --- LAPACK reduction banded -> tridiagonal ---
*
            PRBLM = SBTRD
            READ( IUNIT, * ) MTYPE, DIAM, UPLO, N, B1, LDA1, JOBU, LDU,
     >                       XINFO1
            FRMT1 = 'Banded'
            FRMT2 = 'Tridiagonal'
            B2 = 1
            LDA2 = B2 + 1
            LWORK1 = N
          ELSEIF ( LXSAME( PNAME, 'SYRDB' ) ) THEN
*
*              --- SBR reduction full -> banded ---
*
            PRBLM = SYRDB
            READ( IUNIT, * ) MTYPE, DIAM, UPLO, N, LDA1, B2, DRPTOL,
     >                       NB( 1 ), LWORK1, XINFO1
            LINENO = LINENO + 1
            READ( IUNIT, * ) JOBU0, LDU, LWORK2, XINFO2
            FRMT1 = 'Full'
            B1 = MAX( N-1, 0 )
            DRPTOL = DRPTOL * EPS
            FRMT2 = 'Banded'
            LDA2 = B2 + 1
            IF ( LSAME( JOBU0, 'OnFly' ) ) THEN
              JOBU = 'Update'
            ELSEIF ( LSAME( JOBU0, 'Update' ) ) THEN
              JOBU = 'NoU'
            ELSE
              JOBU = JOBU0
            ENDIF
          ELSEIF ( LXSAME( PNAME, 'SBRDB' ) ) THEN
*
*              --- SBR reduction banded -> banded ---
*
            PRBLM = SBRDB
            READ( IUNIT, * ) MTYPE, DIAM, N, B1, LDA1, B2, DRPTOL,
     >                       NB( 1 ), JOBU, LDU, LWORK1, XINFO1
            FRMT1 = 'Banded'
            UPLO = 'Lower'
            DRPTOL = DRPTOL * EPS
            FRMT2 = 'Banded'
            LDA2 = B2 + 1
          ELSEIF ( LXSAME( PNAME, 'SBRDT' ) ) THEN
*
*              --- SBR reduction banded -> tridiagonal ---
*
            PRBLM = SBRDT
            READ( IUNIT, * ) MTYPE, DIAM, N, B1, LDA1, DRPTOL, NB( 1 ),
     >                       JOBU, LDU, LWORK1, XINFO1
            FRMT1 = 'Banded'
            UPLO = 'Lower'
            DRPTOL = DRPTOL * EPS
            FRMT2 = 'Tridiagonal'
            B2 = 1
            LDA2 = B2 + 1
          ELSEIF ( LXSAME( PNAME, 'SY2BC' ) ) THEN
*
*              --- SBR copy full -> banded ---
*
            PRBLM = SY2BC
            READ( IUNIT, * ) MTYPE, DIAM, UPLO, N, B1, LDA1, LDA2,
     >                       XINFO1
            FRMT1 = 'Full'
            FRMT2 = 'Banded'
            B2 = B1
            LWORK1 = 0
            JOBU = 'NoU'
            LDU = N
          ELSEIF ( LXSAME( PNAME, 'SY2BI' ) ) THEN
*
*              --- SBR in-place copy full -> banded ---
*
            PRBLM = SY2BI
            READ( IUNIT, * ) MTYPE, DIAM, UPLO, N, B1, LDA1, LDA2,
     >                       XINFO1
            FRMT1 = 'Full'
            FRMT2 = 'Banded'
            B2 = B1
            LWORK1 = 0
            JOBU = 'NoU'
            LDU = N
          ELSEIF ( LXSAME( PNAME, 'SB2BC' ) ) THEN
*
*              --- SBR copy banded -> banded ---
*
            PRBLM = SB2BC
            READ( IUNIT, * ) MTYPE, DIAM, UPLO, N, B1, LDA1, LDA2,
     >                       XINFO1
            FRMT1 = 'Banded'
            FRMT2 = 'Banded'
            B2 = B1
            LWORK1 = 0
            JOBU = 'NoU'
            LDU = N
          ELSEIF ( LXSAME( PNAME, 'SB2BI' ) ) THEN
*
*              --- SBR in-place copy full -> banded ---
*
            PRBLM = SB2BI
            READ( IUNIT, * ) MTYPE, DIAM, UPLO, N, B1, LDA1, LDA2,
     >                       XINFO1
            FRMT1 = 'Banded'
            FRMT2 = 'Banded'
            B2 = B1
            LWORK1 = 0
            JOBU = 'NoU'
            LDU = N
          ELSEIF ( LXSAME( PNAME, 'SYRDD' ) ) THEN
*
*              --- SBR driver for full matrices ---
*
            PRBLM = SYRDD
            READ( IUNIT, * ) MTYPE, DIAM, UPLO, N, LDA1, B2, LDA2,
     >                       NSTEPS, DRPTOL, JOBU, LWORK1, XINFO1
            FRMT1 = 'Full'
            B1 = MAX( N-1, 0 )
            DRPTOL = DRPTOL * EPS
            IF ( B2 .EQ. 1 ) THEN
              FRMT2 = 'Tridiagonal'
            ELSE
              FRMT2 = 'Banded'
            ENDIF
            LDU = N
            IF ( NSTEPS .GT. 1 ) THEN
              LINENO = LINENO + 1
              READ( IUNIT, * ) ( B( I ), I = 1, NSTEPS-1 )
            ENDIF
            IF ( NSTEPS .GT. 0 ) THEN
              LINENO = LINENO + 1
              READ( IUNIT, * ) ( NB( I ), I = 1, NSTEPS )
            ENDIF
          ELSEIF ( LXSAME( PNAME, 'SBRDD' ) ) THEN
*
*              --- SBR driver for banded matrices ---
*
            PRBLM = SBRDD
            READ( IUNIT, * ) MTYPE, DIAM, UPLO, N, B1, LDA1, B2, NSTEPS,
     >                       DRPTOL, JOBU, LDU, LWORK1, XINFO1
            FRMT1 = 'Banded'
            DRPTOL = DRPTOL * EPS
            IF ( B2 .EQ. 1 ) THEN
              FRMT2 = 'Tridiagonal'
            ELSE
              FRMT2 = 'Banded'
            ENDIF
            LDA2 = B2 + 1
            IF ( NSTEPS .GT. 1 ) THEN
              LINENO = LINENO + 1
              READ( IUNIT, * ) ( B( I ), I = 1, NSTEPS-1 )
            ENDIF
            IF ( NSTEPS .GT. 0 ) THEN
              LINENO = LINENO + 1
              READ( IUNIT, * ) ( NB( I ), I = 1, NSTEPS )
            ENDIF
*
*            --- done ? ---
*
          ELSEIF ( LSAME( PNAME( 1:1 ), 'Quit' ) ) THEN
            GOTO 8000
          ELSE
            WRITE( OUNIT, 9300 ) LINENO, 'Command not recognized'
            GOTO 2000
          ENDIF
*
*  - - - - - - - - if we get here then a test run is required - - - - -
*
          LINENO = LINENO + 2
          READ ( IUNIT, * ) CHECK, TORTHU, TRES, TDEVLS
          ERRCHK = .FALSE.
*
          TRY( PRBLM ) = TRY( PRBLM ) + 1
*
*            --- is the problem too large ? ---
*
          IF ( N .GT. MAXDIM ) THEN
            SKIP( PRBLM ) = SKIP( PRBLM ) + 1
            WRITE( OUNIT, 9310 ) LINE0
            GOTO 2000
          ENDIF
*
*            ---  initialize matrices ---
*
          INFO = 0
          BTIME = SECOND()
          CALL SSYINI( MTYPE, EVLS, NMAX, DIAM,
     >                 FRMT1, UPLO, N, B1, A, LDA1, SZA,
     >                 JOBU, U, LDU, SZU,
     >                 NOTUSD, WORK, SZWORK, INFO )
          CTIME = SECOND() - BTIME
*
*            --- protocol output, if required ---
*
          IF ( OLEVEL .GE. OUTDET ) THEN
            WRITE( OUNIT, 9500 ) MTYPE, NMAX, DIAM, FRMT1, UPLO, N
            WRITE( OUNIT, 9501 ) B1, LDA1, SZA, JOBU, LDU, SZU
            WRITE( OUNIT, 9502 ) NOTUSD, SZWORK
            WRITE( OUNIT, 9400 ) INFO, CTIME
          ENDIF
*
          IF ( ( INFO .EQ. -3 ) .OR. ( INFO .EQ. -11 ) .OR.
     >         ( INFO .EQ. -15 ) .OR. ( INFO .EQ. -18 ) ) THEN
*
*              --- test problem is too large : skip ---
*
            SKIP( PRBLM ) = SKIP( PRBLM ) + 1
            WRITE( OUNIT, 9310 ) LINE0
            GOTO 3000
          ENDIF
*
*            --- mark workspace "unused" ---
*
          DO 2100 I = 1, MIN( MAX( 2*LWORK1, 10000 ), SZWORK )
            WORK( I ) = NOTUSD
 2100     CONTINUE
*
*  - - - - - - - - do the reduction / repacking - - - - - - - - - - - -
*
          INFO = 0
          BTIME = SECOND()
          IF ( PRBLM .EQ. SYTRD ) THEN
            CALL SSYTRX( NB( 1 ), NX1, UPLO, N, A, LDA1, D, E,
     >                   TAU, WORK, LWORK1, INFO )
          ELSEIF ( PRBLM .EQ. SBTRD ) THEN
            CALL SSBTRX( JOBU, UPLO, N, B1, A, LDA1, D, E,
     >                   U, LDU, WORK, INFO )
          ELSEIF ( PRBLM .EQ. SYRDB ) THEN
            CALL SSYRDB( UPLO, JOBU, N, B2, A, LDA1, DRPTOL,
     >                   U, LDU, NB( 1 ), TAU, WORK, LWORK1, INFO )
          ELSEIF ( PRBLM .EQ. SBRDB ) THEN
            CALL SSBRDB( JOBU, N, B1, B2, A, LDA1, DRPTOL,
     >                   U, LDU, NB( 1 ), WORK, LWORK1, INFO )
          ELSEIF ( PRBLM .EQ. SBRDT ) THEN
            CALL SSBRDT( JOBU, N, B1, A, LDA1, DRPTOL, D, E,
     >                   U, LDU, NB( 1 ), WORK, LWORK1, INFO )
          ELSEIF ( PRBLM .EQ. SY2BC ) THEN
            CALL SSY2BC( UPLO, N, B1, A, LDA1, ACPY, LDA2, INFO )
          ELSEIF ( PRBLM .EQ. SY2BI ) THEN
            CALL SSY2BI( UPLO, N, B1, A, LDA1, LDA2, INFO )
          ELSEIF ( PRBLM .EQ. SB2BC ) THEN
            CALL SSB2BC( UPLO, N, B1, A, LDA1, ACPY, LDA2, INFO )
          ELSEIF ( PRBLM .EQ. SB2BI ) THEN
            CALL SSB2BI( UPLO, N, B1, A, LDA1, LDA2, INFO )
          ELSEIF ( PRBLM .EQ. SYRDD ) THEN
            CALL SSYRDD( JOBU, UPLO, N, B2, A, LDA1, DRPTOL,
     >                   ACPY, LDA2, D, E, NSTEPS, B, NB,
     >                   WORK, LWORK1, STEP, INFO )
          ELSEIF ( PRBLM .EQ. SBRDD ) THEN
            CALL SSBRDD( JOBU, UPLO, N, B1, B2, A, LDA1, DRPTOL,
     >                   D, E, U, LDU, NSTEPS, B, NB,
     >                   WORK, LWORK1, STEP, INFO )
          ENDIF
          CTIME = SECOND() - BTIME
          RTIME = CTIME
*
*            --- protocol output, if required ---
*
          IF ( OLEVEL .GE. OUTDET ) THEN
            IF ( PRBLM .EQ. SYTRD ) THEN
              WRITE( OUNIT, 9510 ) UPLO, N, LDA1, NB( 1 ), NX1, LWORK1
              WRITE( OUNIT, 9400 ) INFO, CTIME
            ELSEIF ( PRBLM .EQ. SBTRD ) THEN
              WRITE( OUNIT, 9520 ) JOBU, UPLO, N, B1, LDA1, LDU
              WRITE( OUNIT, 9400 ) INFO, CTIME
            ELSEIF ( PRBLM .EQ. SYRDB ) THEN
              WRITE( OUNIT, 9530 ) UPLO, JOBU, N, B2, LDA1, DRPTOL
              WRITE( OUNIT, 9531 ) LDU, NB( 1 ), LWORK1
              WRITE( OUNIT, 9400 ) INFO, CTIME
            ELSEIF ( PRBLM .EQ. SBRDB ) THEN
              WRITE( OUNIT, 9540 ) JOBU, N, B1, B2, LDA1, DRPTOL
              WRITE( OUNIT, 9531 ) LDU, NB( 1 ), LWORK1
              WRITE( OUNIT, 9400 ) INFO, CTIME
            ELSEIF ( PRBLM .EQ. SBRDT ) THEN
              WRITE( OUNIT, 9550 ) JOBU, N, B1, LDA1, DRPTOL
              WRITE( OUNIT, 9531 ) LDU, NB( 1 ), LWORK1
              WRITE( OUNIT, 9400 ) INFO, CTIME
            ELSEIF ( PRBLM .EQ. SY2BC ) THEN
              WRITE( OUNIT, 9560 ) UPLO, N, B1, LDA1, LDA2
              WRITE( OUNIT, 9400 ) INFO, CTIME
            ELSEIF ( PRBLM .EQ. SY2BI ) THEN
              WRITE( OUNIT, 9570 ) UPLO, N, B1, LDA1, LDA2
              WRITE( OUNIT, 9400 ) INFO, CTIME
            ELSEIF ( PRBLM .EQ. SB2BC ) THEN
              WRITE( OUNIT, 9580 ) UPLO, N, B1, LDA1, LDA2
              WRITE( OUNIT, 9400 ) INFO, CTIME
            ELSEIF ( PRBLM .EQ. SB2BI ) THEN
              WRITE( OUNIT, 9590 ) UPLO, N, B1, LDA1, LDA2
              WRITE( OUNIT, 9400 ) INFO, CTIME
            ELSEIF ( PRBLM .EQ. SYRDD ) THEN
              WRITE( OUNIT, 9600 ) JOBU, UPLO, N, B2, LDA1
              WRITE( OUNIT, 9601 ) DRPTOL, LDA2, NSTEPS, LWORK1
              WRITE( OUNIT, 9401 ) STEP, INFO, CTIME
            ELSEIF ( PRBLM .EQ. SBRDD ) THEN
              WRITE( OUNIT, 9610 ) JOBU, UPLO, N, B1, B2, LDA1
              WRITE( OUNIT, 9611 ) DRPTOL, LDU, NSTEPS, LWORK1
              WRITE( OUNIT, 9401 ) STEP, INFO, CTIME
            ENDIF
          ELSEIF ( OLEVEL .EQ. OUTTUN ) THEN
            IF ( PRBLM .EQ. SYTRD ) THEN
              WRITE( OUNIT, 9515 ) UPLO, N, NB( 1 ), INFO, CTIME
            ELSEIF ( PRBLM .EQ. SBTRD ) THEN
              WRITE( OUNIT, 9525 ) UPLO, JOBU, N, B1, INFO, CTIME
            ELSEIF ( PRBLM .EQ. SYRDB ) THEN
              WRITE( OUNIT, 9535 ) UPLO, JOBU, N, B2, NB( 1 ), INFO,
     >                             CTIME
            ELSEIF ( PRBLM .EQ. SBRDB ) THEN
              WRITE( OUNIT, 9545 ) JOBU, N, B1, B2, NB( 1 ), INFO, CTIME
            ELSEIF ( PRBLM .EQ. SBRDT ) THEN
              WRITE( OUNIT, 9555 ) JOBU, N, B1, NB( 1 ), INFO, CTIME
            ELSEIF ( PRBLM .EQ. SYRDD ) THEN
              WRITE( OUNIT, 9605 ) UPLO, JOBU, N, B2, NSTEPS, STEP,
     >                             INFO, CTIME
            ELSEIF ( PRBLM .EQ. SBRDD ) THEN
              WRITE( OUNIT, 9615 ) UPLO, JOBU, N, B1, B2, NSTEPS, STEP,
     >                             INFO, CTIME
            ENDIF
          ENDIF
*
*  - - - - - - - - check if the error code is as expected - - - - - - -
*
          IF ( INFO .NE. XINFO1 ) THEN
*
*              --- wrong error code ---
*
            FAIL( PRBLM ) = FAIL( PRBLM ) + 1
            IF ( ( PRBLM .EQ. SYRDD ) .OR. ( PRBLM .EQ. SBRDD ) ) THEN
              WRITE( OUNIT, 9900 ) LINE0, PRBLMN( PRBLM ), INFO, STEP,
     >                             XINFO1
            ELSE
              WRITE( OUNIT, 9901 ) LINE0, PRBLMN( PRBLM ), INFO, XINFO1
            ENDIF
            GOTO 3000
          ELSEIF ( ( ( PRBLM .EQ. SYRDD ) .OR. ( PRBLM .EQ. SBRDD ) )
     >             .AND. ( INFO .EQ. 1 ) .AND. ( NSTEPS .GT. 0 )
     >             .AND. ( STEP .NE. NSTEPS ) ) THEN
*
*              --- error code was ok, but the driver took a wrong
*                  number of steps                                ---
*
            FAIL( PRBLM ) = FAIL( PRBLM ) + 1
            WRITE( OUNIT, 9910 ) LINE0, PRBLMN( PRBLM ), STEP,
     >                           NSTEPS
            GOTO 3000
          ELSEIF ( INFO .LT. 0 ) THEN
*
*              --- successful error exit ---
*
            ERRCHK = .TRUE.
            GOTO 3000
          ENDIF
*
*  - - - - - - - - check workspace usage - - - - - - - - - - - - - - - -
*
          NUSED = 0
          DO 2200 I = 1, MIN( MAX( 2*LWORK1, 10000 ), SZWORK )
            IF ( .NOT. ( WORK( I ) .EQ. NOTUSD ) )     NUSED = I
 2200     CONTINUE
*
          IF ( NUSED .GT. LWORK1 ) THEN
*
*              --- ran over the workspace bounds ---
*
            FAIL( PRBLM ) = FAIL( PRBLM ) + 1
            WRITE( OUNIT, 9920 ) LINE0, PRBLMN( PRBLM ), NUSED, LWORK1
            GOTO 3000
          ELSEIF ( OLEVEL .GE. OUTDET ) THEN
*
*            --- write protocol ---
*
            WRITE( OUNIT, 9410 ) NUSED, LWORK1
          ENDIF
*
*  - - - - - - - - accumulate transformations, if required - - - - - - -
*
          IF ( ( PRBLM .EQ. SYTRD ) .AND. LSAME( JOBU, 'Update' ) ) THEN
*
*              --- copy the Householder vectors into the array u ---
*
            CALL SLACPY( UPLO, N, N, A, LDA1, U, LDU )
*
*              --- backward accumulation ---
*
            TRY( ORGTR ) = TRY( ORGTR ) + 1
            INFO = 0
            BTIME = SECOND()
            CALL SORGTX( NB( 2 ), NX2, UPLO, N, U, LDU, TAU,
     >                   WORK, LWORK2, INFO )
            CTIME = SECOND() - BTIME
            RTIME = RTIME + CTIME
*
*              --- protocol output, if required ---
*
            IF ( OLEVEL .GE. OUTDET ) THEN
              WRITE( OUNIT, 9620 ) UPLO, N, LDU, NB( 2 ), NX2, LWORK2
              WRITE( OUNIT, 9400 ) INFO, CTIME
            ENDIF
*
*              --- is the error code as expected ? ---
*
            IF ( INFO .NE. XINFO2 ) THEN
              FAIL( ORGTR ) = FAIL( ORGTR ) + 1
              WRITE( OUNIT, 9901 ) LINE0, 'sorgtr', INFO, XINFO2
              GOTO 3000
            ELSEIF ( INFO .NE. 0 ) THEN
              ERRCHK = .TRUE.
              GOTO 3000
            ENDIF
*
          ELSEIF ( ( PRBLM .EQ. SYRDB ) .AND. LSAME( JOBU0, 'Update' ) )
     >    THEN
*
*              --- copy the Householder vectors into the array u ---
*
            CALL SLACPY( UPLO, N, N, A, LDA1, U, LDU )
*
*              --- backward accumulation ---
*
            TRY( SYGTR ) = TRY( SYGTR ) + 1
            INFO = 0
            BTIME = SECOND()
            CALL SSYGTR( UPLO, N, B2, U, LDU, TAU, WORK, LWORK2, INFO )
            CTIME = SECOND() - BTIME
            RTIME = RTIME + CTIME
*
*              --- protocol output, if required ---
*
            IF ( OLEVEL .GE. OUTDET ) THEN
              WRITE( OUNIT, 9630 ) UPLO, N, B2, LDU, LWORK2
              WRITE( OUNIT, 9400 ) INFO, CTIME
            ENDIF
*
*              --- is the error code as expected ? ---
*
            IF ( INFO .NE. XINFO2 ) THEN
              FAIL( SYGTR ) = FAIL( SYGTR ) + 1
              WRITE( OUNIT, 9901 ) LINE0, 'ssygtr', INFO, XINFO2
              GOTO 3000
            ELSEIF ( INFO .NE. 1 ) THEN
              GOTO 3000
            ENDIF
          ENDIF
*
*  - - - - - - - - check results - - - - - - - - - - - - - - - - - - - -
*
          ORTHU = -ONE
          SORTHU = -ONE
          RES = -ONE
          SRES = -ONE
          DEVLS = -ONE
          SDEVLS = -ONE
*
          IF ( ( .NOT. ERRCHK ) .AND.
     >         CHECKS .AND. LSAME( CHECK, 'Check' ) ) THEN
            OK = .TRUE.
*
*              --- if the reduced matrix is banded then copy it
*                  into the array acpy                          ---
*
            IF ( PRBLM .EQ. SYRDB ) THEN
              INFO = 0
              BTIME = SECOND()
              CALL SSY2BC( UPLO, N, B2, A, LDA1, ACPY, LDA2, INFO )
              CTIME = SECOND() - BTIME
              IF ( OLEVEL .GE. OUTDET ) THEN
                WRITE( OUNIT, 9560 ) UPLO, N, B2, LDA1, LDA2
                WRITE( OUNIT, 9400 ) INFO, CTIME
              ENDIF
              IF ( LSAME( JOBU0, 'Update' ) )      JOBU = 'Update'
            ELSEIF ( ( PRBLM .EQ. SBRDB ) .OR.
     >               ( ( PRBLM .EQ. SBRDD ) .AND. ( B2 .GT. 1 ) ) ) THEN
              INFO = 0
              BTIME = SECOND()
              CALL SSB2BC( 'Lower', N, B2, A, LDA1, ACPY, LDA2, INFO )
              CTIME = SECOND() - BTIME
              IF ( OLEVEL .GE. OUTDET ) THEN
                WRITE( OUNIT, 9580 ) 'Lower', N, B2, LDA1, LDA2
                WRITE( OUNIT, 9400 ) INFO, CTIME
              ENDIF
            ELSEIF ( ( PRBLM .EQ. SY2BI ) .OR. ( PRBLM .EQ. SB2BI ) )
     >      THEN
              INFO = 0
              BTIME = SECOND()
              CALL SSB2BC( 'Lower', N, B2, A, LDA2, ACPY, LDA2, INFO )
              CTIME = SECOND() - BTIME
              IF ( OLEVEL .GE. OUTDET ) THEN
                WRITE( OUNIT, 9580 ) UPLO, N, B2, LDA2, LDA2
                WRITE( OUNIT, 9400 ) INFO, CTIME
              ENDIF
            ELSEIF ( PRBLM .EQ. SYRDD ) THEN
*
*                --- copy orthogonal matrix into u ---
*
              CALL SLACPY( 'Full', N, N, A, LDA1, U, LDU )
            ENDIF
*
*              --- perform numerical checks ---
*
            INFO = 0
            BTIME = SECOND()
            CALL SSYNCK( EVLS, FRMT1, UPLO, N, B1, A, LDA1, SZA,
     >                   FRMT2, B2, ACPY, LDA2, D, E, JOBU, U, LDU,
     >                   ORTHU, SORTHU, RES, SRES, DEVLS, SDEVLS,
     >                   NOTUSD, WORK, SZWORK, INFO )
            CTIME = SECOND() - BTIME
*
*              --- protocol output, if required ---
*
            IF ( OLEVEL .GE. OUTDET ) THEN
              WRITE( OUNIT, 9640 ) FRMT1, UPLO, N, B1, LDA1, SZA
              WRITE( OUNIT, 9641 ) FRMT2, B2, LDA2, JOBU, LDU
              WRITE( OUNIT, 9642 ) ORTHU, SORTHU
              WRITE( OUNIT, 9643 ) RES, SRES
              WRITE( OUNIT, 9644 ) DEVLS, SDEVLS
              WRITE( OUNIT, 9645 ) NOTUSD, SZWORK
              WRITE( OUNIT, 9400 ) INFO, CTIME
            ENDIF
*
*              --- failed any computational check ? ---
*
            IF ( .NOT. ( SORTHU .LE. TORTHU ) ) THEN
              WRITE( OUNIT, 9930 ) LINE0, 'orhogonality', SORTHU, TORTHU
              OK = .FALSE.
            ENDIF
            IF ( .NOT. ( SRES .LE. TRES ) ) THEN
              WRITE( OUNIT, 9930 ) LINE0, 'residual', SRES, TRES
              OK = .FALSE.
            ENDIF
            IF ( .NOT. ( SDEVLS .LE. TDEVLS ) ) THEN
              WRITE( OUNIT, 9930 ) LINE0, 'eigenvalues', SDEVLS, TDEVLS
              OK = .FALSE.
            ENDIF
            IF ( .NOT. OK )     FAIL( PRBLM ) = FAIL( PRBLM ) + 1
          ENDIF
*
*  - - - - - - - - print summary for this run - - - - - - - - - - - - -
*
 3000     IF ( OLEVEL .GE. OUTRUN ) THEN
            WRITE( OUNIT, 9700 ) LINE0, PRBLMN( PRBLM ), COMMNT( PRBLM )
*
            IF ( PRBLM .EQ. SYTRD ) THEN
              WRITE( OUNIT, 9710 ) MTYPE, DIAM, UPLO, N, LDA1
              WRITE( OUNIT, 9711 ) NB( 1 ), NX1, LWORK1, XINFO1
              WRITE( OUNIT, 9712 ) JOBU, NB( 2 ), NX2, LWORK2, XINFO2
            ELSEIF ( PRBLM .EQ. SBTRD ) THEN
              WRITE( OUNIT, 9720 ) MTYPE, DIAM, UPLO, N, B1, LDA1
              WRITE( OUNIT, 9721 ) JOBU, LDU, XINFO1
            ELSEIF ( PRBLM .EQ. SYRDB ) THEN
              WRITE( OUNIT, 9730 ) MTYPE, DIAM, UPLO, N, LDA1, B2
              WRITE( OUNIT, 9731 ) DRPTOL, NB( 1 ), LWORK1, XINFO1
              WRITE( OUNIT, 9732 ) JOBU, LDU, LWORK2, XINFO2
            ELSEIF ( PRBLM .EQ. SBRDB ) THEN
              WRITE( OUNIT, 9740 ) MTYPE, DIAM, N, B1, LDA1, B2
              WRITE( OUNIT, 9741 ) DRPTOL, NB( 1 ), JOBU, LDU, LWORK1
              WRITE( OUNIT, 9742 ) XINFO1
            ELSEIF ( PRBLM .EQ. SBRDT ) THEN
              WRITE( OUNIT, 9750 ) MTYPE, DIAM, N, B1, LDA1
              WRITE( OUNIT, 9741 ) DRPTOL, NB( 1 ), JOBU, LDU, LWORK1
              WRITE( OUNIT, 9742 ) XINFO1
            ELSEIF ( ( PRBLM .EQ. SY2BC )
     >               .OR. ( PRBLM .EQ. SY2BI )
     >               .OR. ( PRBLM .EQ. SB2BC )
     >               .OR. ( PRBLM .EQ. SB2BI ) ) THEN
              WRITE( OUNIT, 9760 ) MTYPE, DIAM, UPLO, N, B1, LDA1
              WRITE( OUNIT, 9761 ) LDA2, XINFO1
            ELSEIF ( PRBLM .EQ. SYRDD ) THEN
              WRITE( OUNIT, 9770 ) MTYPE, DIAM, UPLO, N, LDA1, B2
              WRITE( OUNIT, 9771 ) LDA2, NSTEPS, DRPTOL, JOBU, LWORK1
              WRITE( OUNIT, 9742 ) XINFO1
              IF ( NSTEPS .GT. 1 ) THEN
                WRITE( OUNIT, 9772 ) ( B( I ), I = 1, NSTEPS-1 )
              ENDIF
              IF ( NSTEPS .GT. 0 ) THEN
                WRITE( OUNIT, 9773 ) ( NB( I ), I = 1, NSTEPS )
              ENDIF
            ELSEIF ( PRBLM .EQ. SBRDD ) THEN
              WRITE( OUNIT, 9780 ) MTYPE, DIAM, UPLO, N, B1, LDA1, B2
              WRITE( OUNIT, 9781 ) NSTEPS, DRPTOL, JOBU, LDU
              WRITE( OUNIT, 9782 ) LWORK1, XINFO1
              IF ( NSTEPS .GT. 1 ) THEN
                WRITE( OUNIT, 9772 ) ( B( I ), I = 1, NSTEPS-1 )
              ENDIF
              IF ( NSTEPS .GT. 0 ) THEN
                WRITE( OUNIT, 9773 ) ( NB( I ), I = 1, NSTEPS )
              ENDIF
            ENDIF
*
*              --- print available information on numerical checks ---
*
            IF ( ERRCHK ) THEN
              WRITE( OUNIT, 9420 )
            ELSE
              IF ( .NOT. ( SORTHU .EQ. -ONE ) ) THEN
                WRITE( OUNIT, 9430 ) ORTHU, SORTHU, TORTHU
              ENDIF
              IF ( .NOT. ( SRES .EQ. -ONE ) ) THEN
                WRITE( OUNIT, 9431 ) RES, SRES, TRES
              ENDIF
              IF ( .NOT. ( SDEVLS .EQ. -ONE ) ) THEN
                WRITE( OUNIT, 9432 ) DEVLS, SDEVLS, TDEVLS
              ENDIF
            ENDIF
*
            WRITE( OUNIT, 9440 ) RTIME
          ENDIF
*
          IF ( OLEVEL .GE. OUTDET )     WRITE( OUNIT, 9000 )
*
*            --- end of the main loop ---
*
        GOTO 2000
*
*  - - - - - - - - print "global" summary - - - - - - - - - - - - - - -
*
 8000   WRITE( OUNIT, * )
        WRITE( OUNIT, 9010 )
        WRITE( OUNIT, 9020 )
        WRITE( OUNIT, 9000 )
        DO 8100 PRBLM = SYTRD, TOTAL-1
          WRITE( OUNIT, 9030 ) PRBLMN( PRBLM ), COMMNT( PRBLM ),
     >                         TRY( PRBLM ), SKIP( PRBLM ),
     >                         FAIL( PRBLM )
          TRY( TOTAL ) = TRY( TOTAL ) + TRY( PRBLM )
          SKIP( TOTAL ) = SKIP( TOTAL ) + SKIP( PRBLM )
          FAIL( TOTAL ) = FAIL( TOTAL ) + FAIL( PRBLM )
 8100   CONTINUE
        WRITE( OUNIT, 9000 )
        WRITE( OUNIT, 9030 ) PRBLMN( TOTAL ), COMMNT( TOTAL ),
     >                       TRY( TOTAL ), SKIP( TOTAL ), FAIL( TOTAL )
        WRITE( OUNIT, 9000 )
        WRITE( OUNIT, 9040 ) SECOND() - STIME
*
*  - - - - - - - - formats - - - - - - - - - - - - - - - - - - - - - - -
*
*          --- "global" summary" ---
*
 9000   FORMAT( 1X, '----------------------------------------',
     >          '---------------------------------------' )
 9010   FORMAT( 1X, '========================================',
     >          '=======================================' )
*
 9020   FORMAT( 1X, 'Tests for                               ',
     >          '              Total Skipped  Failed' )
 9030   FORMAT( 1X, A6, ' : ', A40, ' :', I8, I8, I8 )
*
 9040   FORMAT( 1X, 'Total time : ', F8.2, ' seconds' )
*
*          --- messages ---
*
 9100   FORMAT( 1X, A )
*
*          --- problems with the input file ---
*
 9300   FORMAT( 1X, '*** Error in line ', I5, ' : ', A40 )
 9310   FORMAT( 1X, '+++ Line ', I5, ' : Test too large - skipped' )
*
*          --- general status reporting ---
*
 9400   FORMAT( 1X, '          info=', I5, ', time=', F12.6 )
 9401   FORMAT( 1X, '          step=', I5, ', info=', I5,
     >          ', time=', F12.6 )
*
 9410   FORMAT( 1X, '          used ', I8, ' out of ',
     >          I8, ' workspace' )
*
 9420   FORMAT( 1X, '    ((( check for error exit ))),' )
*
 9430   FORMAT( 1X, '    Norm(U''*U-I) =', E10.3, '=',
     >          E10.3, '*n*macheps         (thresh=', E10.3, '),' )
 9431   FORMAT( 1X, '    Norm(U*B-A*U)=', E10.3, '=',
     >          E10.3, '*n*macheps*Norm(A) (thresh=', E10.3, '),' )
 9432   FORMAT( 1X, '    Norm(Dlambda)=', E10.3, '=',
     >          E10.3, '*n*macheps*Norm(A) (thresh=', E10.3, '),' )
*
 9440   FORMAT( 1X, '    reduction time=', F12.6 )
*
*          --- formats for tracing single routine calls ---
*
 9500   FORMAT( 1X, '  ssyini: mtype=', A1, ', szevls=', I5,
     >          ', diam=', E11.4, ', frmt=', A1, ', uplo=', A1,
     >          ', n=', I5, ',' )
 9501   FORMAT( 1X, '          b=', I5, ', lda=', I5, ', sza=', I8,
     >          ', jobu=', A1, ', ldu=', I5, ', szu=', I8, ',' )
 9502   FORMAT( 1X, '          notusd=', E11.4, ', szwork=', I8, ',' )
*
 9510   FORMAT( 1X, '  ssytrd: uplo=', A1, ', n=', I5, ', lda=', I5,
     >          ', nb0=', I5, ', nx0=', I5, ', lwork=', I8, ',' )
 9515   FORMAT( 1X, 'ssytrd: uplo=', A1, ', n=', I5, ', nb=', I3,
     >          ', info=', I3, ', time=', F9.3 )
*
 9520   FORMAT( 1X, '  ssbtrd: vect=', A1, ', uplo=', A1, ', n=', I5,
     >          ', kd=', I5, ', ldab=', I5, ', ldu=', I5, ',' )
 9525   FORMAT( 1X, 'ssbtrd: uplo=', A1, ', jobu=', A1, ', n=', I5,
     >          ', b1=', I3, ', info=', I3, ', time=', F9.3 )
*
 9530   FORMAT( 1X, '  ssyrdb: uplo=', A1, ', job=', A1, ', n=', I5,
     >          ', b=', I5, ', lda=', I5, ', drptol=', E11.4, ',' )
 9531   FORMAT( 1X, '          ldu=', I5, ', nb=', I5, ', lwork=', I8,
     >          ',' )
 9535   FORMAT( 1X, 'ssyrdb: uplo=', A1, ', jobu=', A1, ', n=', I5,
     >          ', b2=', I3, ', nb=', I3, ', info=', I3,
     >          ', time=', F9.3 )
*
 9540   FORMAT( 1X, '  ssbrdb: job=', A1, ', n=', I5, ', b1=', I5,
     >          ', b2=', I5, ', lda=', I5, ', drptol=', E11.4, ',' )
 9545   FORMAT( 1X, 'ssbrdb: jobu=', A1, ', n=', I5, ', b1=', I3,
     >          ', b2=', I3, ', nb=',I3, ', info=', I3,
     >          ', time=', F9.3 )
*
 9550   FORMAT( 1X, '  ssbrdt: job=', A1, ', n=', I5, ', b=', I5,
     >          ', lda=', I5, ', drptol=', E11.4, ',' )
 9555   FORMAT( 1X, 'ssbrdt: jobu=', A1, ', n=', I5, ', b1=', I3,
     >          ', nb=', I3, ', info=', I3, ', time=', F9.3 )
*
 9560   FORMAT( 1X, '  ssy2bc: uplo=', A1, ', n=', I5, ', b=', I5,
     >          ', ldfull=', I5, ', ldband=', I5, ',' )
*
 9570   FORMAT( 1X, '  ssy2bi: uplo=', A1, ', n=', I5, ', b=', I5,
     >          ', ldfull=', I5, ', ldband=', I5, ',' )
*
 9580   FORMAT( 1X, '  ssb2bc: uplo=', A1, ', n=', I5, ', b=', I5,
     >          ', ldsrc=', I5, ', lddst=', I5, ',' )
*
 9590   FORMAT( 1X, '  ssb2bi: uplo=', A1, ', n=', I5, ', b=', I5,
     >          ', ldsrc=', I5, ', lddst=', I5, ',' )
*
 9600   FORMAT( 1X, '  ssyrdd: job=', A1, ', uplo=', A1, ', n=', I5,
     >          ', b2=', I5, ', lda=', I5, ',' )
 9601   FORMAT( 1X, '          drptol=', E11.4, ', ldband=', I5,
     >          ', nsteps=', I5, ', lwork=', I8, ',' )
 9605   FORMAT( 1X, 'ssyrdd: uplo=', A1, ', job=', A1, ', n=', I5,
     >          ', b2=', I3, ', nsteps=', I2, ', step=', I2,
     >          ', info=', I3, ', time=', F9.3 )
*
 9610   FORMAT( 1X, '  ssbrdd: jobu=', A1, ', uplo=', A1, ', n=', I5,
     >          ', b1=', I5, ', b2=', I5, ', lda=', I5, ',' )
 9611   FORMAT( 1X, '          drptol=', E11.4, ', ldu=', I5,
     >          ', nsteps=', I5, ', lwork=', I8, ',' )
 9615   FORMAT( 1X, 'ssbrdd: uplo=', A1, ', jobu=', A1, ', n=', I5,
     >          ', b1=', I3, ', b2=', I3, ', nsteps=', I2,
     >          ', step=', I2, ', info=', I3, ', time=', F9.3 )
*
 9620   FORMAT( 1X, '  sorgtr: uplo=', A1, ', n=', I5, ', lda=', I5,
     >          ', nb=', I5, ', nx=', I5, ', lwork=', I8, ',' )
*
 9630   FORMAT( 1X, '  ssygtr: uplo=', A1, ', n=', I5, ', b=', I5,
     >          ', lda=', I5, ', lwork=', I8, ',' )
*
 9640   FORMAT( 1X, '  ssynck: frmt1=', A1, ', uplo1=', A1, ', n=', I5,
     >          ', b1=', I5, ', lda1=', I5, ', sza1=', I8, ',' )
 9641   FORMAT( 1X, '          frmt2=', A1, ', b2=', I5, ', lda2=', I5,
     >          ', jobu=', A1, ', ldu=', I5, ',' )
 9642   FORMAT( 1X, '          Norm( U''*U - I ) =', E11.4,
     >          ' = ', E11.4, '*n*macheps,' )
 9643   FORMAT( 1X, '          Norm( U*B - A*U )=', E11.4,
     >          ' = ', E11.4, '*n*macheps*Norm( A ),' )
 9644   FORMAT( 1X, '          Norm( Dlambda )  =', E11.4,
     >          ' = ', E11.4, '*n*macheps*Norm( A ),' )
 9645   FORMAT( 1X, '          notusd=', E11.4, ', szwork=', I8, ',' )
*
*          --- formats for tracing single runs ---
*
 9700   FORMAT( 1X, 'Line ', I5, ' / ', A, ' : ', A, ',' )
*
 9710   FORMAT( 1X, '    mtype=', A1, ', diam=', E11.4,
     >          ', uplo=', A1, ', n=', I5, ', lda=', I5, ',' )
 9711   FORMAT( 1X, '    nb1=', I5, ', nx1=', I5, ', lwork1=', I8,
     >          ', xinfo1=', I5, ',' )
 9712   FORMAT( 1X, '    jobu=', A1, ', nb2=', I5, ', nx2=', I5,
     >          ', lwork2=', I8, ', xinfo2=', I5, ',' )
*
 9720   FORMAT( 1X, '    mtype=', A1, ', diam=', E11.4, ', uplo=', A1,
     >          ', n=', I5, ', b1=', I5, ', lda=', I5, ',' )
 9721   FORMAT( 1X, '    jobu=', A1, ', ldu=', I5, ', xinfo=', I5, ',' )
*
 9730   FORMAT( 1X, '    mtype=', A1, ', diam=', E11.4, ', uplo=', A1,
     >          ', n=', I5, ', lda=', I5, ', b2=', I5, ',' )
 9731   FORMAT( 1X, '    drptol=', E11.4, ', nb=', I5, ', lwork1=', I8,
     >          ', xinfo1=', I5, ',' )
 9732   FORMAT( 1X, '    jobu=', A1, ', ldu=', I5, ', lwork2=', I8,
     >          ', xinfo2=', I5, ',' )
*
 9740   FORMAT( 1X, '    mtype=', A1, ', diam=', E11.4, ', n=', I5,
     >          ', b1=', I5, ', lda=', I5, ', b2=', I5, ',' )
 9741   FORMAT( 1X, '    drptol=', E11.4, ', nb=', I5, ', jobu=', A1,
     >          ', ldu=', I5, ', lwork=', I8, ',' )
 9742   FORMAT( 1X, '    xinfo=', I5, ',' )
*
 9750   FORMAT( 1X, '    mtype=', A1, ', diam=', E11.4,
     >          ', n=', I5, ', b1=', I5, ', lda=', I5, ',' )
*
 9760   FORMAT( 1X, '    mtype=', A1, ', diam=', E11.4, ', uplo=', A1,
     >          ', n=', I5, ', b1=', I5, ', lda1=', I5, ',' )
 9761   FORMAT( 1X, '    lda2=', I5, ', xinfo=', I5, ',' )
*
 9770   FORMAT( 1X, '    mtype=', A1, ', diam=', E11.4, ', uplo=', A1,
     >          ', n=', I5, ', lda1=', I5, ', b2=', I5, ',' )
 9771   FORMAT( 1X, '    lda2=', I5, ', nsteps=', I5,
     >          ', drptol=', E11.4, ', jobu=', A1, ', lwork=', I8, ',' )
 9772   FORMAT( 1X, '    b ( 2: )=', 100(I5, ',') )
 9773   FORMAT( 1X, '    nb( 1: )=', 100(I5, ',') )
*
 9780   FORMAT( 1X, '    mtype=', A1, ', diam=', E11.4, ', uplo=', A1,
     >          ', n=', I5, ', b1=', I5, ', lda=', I5, ', b2=', I5,
     >          ',' )
 9781   FORMAT( 1X, '    nsteps=', I5, ', drptol=', E11.4,
     >          ', jobu=', A1, ', ldu=', I5, ',' )
 9782   FORMAT( 1X, '    lwork=', I8, ', xinfo=', I5, ',' )
*
*          --- unexpected behaviour of the routines ---
*
 9900   FORMAT( 1X, '*** Test in line ', I5, ' failed : ', A6,
     >          ' returned info=', I5, ' in step ', I5,
     >          ' (expected ', I5, ')' )
 9901   FORMAT( 1X, '*** Test in line ', I5, ' failed : ', A6,
     >          ' returned info=', I5, ' (expected ', I5, ')' )
*
 9910   FORMAT( 1X, '*** Test in line ', I5, ' failed : ', A6,
     >          ' returned step=', I5, ' (expected ', I5, ')' )
*
 9920   FORMAT( 1X, '*** Test in line ', I5, ' failed : ', A6,
     >          ' used ', I8, ' workspace (allowed ', I8, ')' )
*
 9930   FORMAT( 1X, '*** Test in line ', I5, ' failed ', A,
     >          ' test: ', E10.3, ' (thresh=', E10.3, ')' )
*
*  - - - - - - - - done - - - - - - - - - - - - - - - - - - - - - - - -
*
        STOP
        END
*
* **********************************************************************
*
        SUBROUTINE SSYINI( MTYPE, EVLS, SZEVLS, DIAM,
     >                     FRMT, UPLO, N, B, A, LDA, SZA,
     >                     JOBU, U, LDU, SZU,
     >                     NOTUSD, WORK, SZWORK, INFO )
*
* Description:
*
*   Generates a symmetric (full or banded) matrix A with either
*   prescribed eigenvalues, random eigenvalues, or two clusters of
*   eigenvalues. The matrix is returned either in (upper or lower)
*   symmetric full storage or in (upper or lower) symmetric banded
*   storage. In addition, U can be set to the identity matrix.
*
*   Note: In contrast to the LAPACK routine slagsy, this routine
*         generates the banded matrix within an n-by-(b+1) array.
*
* Author: Bruno Lang
*         Aachen University of Technology
*         na.blang@na-net.ornl.gov
* Date: May 17, 2000
* Version: SBR Toolbox, Rev. 1.4.1
*
* ----------------------------------------------------------------------
*
* Parameters:
*
        CHARACTER*1           MTYPE, FRMT, UPLO, JOBU
        INTEGER               N, B, LDA, SZA, LDU, SZU, SZWORK, SZEVLS,
     >                        INFO
        REAL                  DIAM, NOTUSD
        REAL                  A( LDA, * ), U( LDU, * ), EVLS( * ),
     >                        WORK( * )
*
*   mtype   (in) character*1
*           The matrix type.
*           mtype = 'S' : Generate a matrix with prescribed eigenvalues.
*                 = 'R' : Eigenvalues are uniformly random in [ -1, 1 ).
*                 = 'I' : ISDA-type matrix: eigenvalues are clustered
*                         around 0 and 1.
*
*   evls    (in OR out) real array, dimension ( szevls )
*           The eigenvalues of the matrix A.
*           If mtype = 'S' then evls is an input parameter, and the
*           first n elements of evls are taken as eigenvalues for the
*           matrix A.
*           If mtype = 'R' or 'I', then evls is an output parameter and
*           contains the eigenvalues of the matrix A.
*
*   szevls  (in) integer
*           The size of the array evls.
*           szevls >= n.
*
*   diam    (in) real
*           The diameter of the eigenvalue clusters (in multiples of
*           machine_epsilon).
*           Accessed only if mtype = 'I'.
*
*   frmt    (in) character*1
*           The storage format for the matrix A.
*           frmt = 'F' : (Symmetric upper or lower) full storage.
*                = 'B' : (Symmetric upper or lower) banded storage.
*
*   uplo    (in) character*1
*           Which triangle of the matrix is used?
*           uplo = 'U' : Initialize the upper triangle.
*                = 'L' : Initialize the lower triangle.
*           In banded storage (frmt = 'B'), uplo also determines
*           whether the upper or lower banded storage scheme is used.
*
*   n       (in) integer
*           The dimension of the matrix.
*           n >= 0.
*
*   b       (in) integer
*           The semibandwidth (i.e., the number of subdiagonals) of the
*           matrix A.
*           0 <= b <= max( n-1, 0 ).
*
*   a       (out) real array, dimension ( lda, n )
*           The n-by-n symmetric matrix with semibandwidth b.
*           If frmt = 'F' then A is returned in the leading n-by-n
*           subarray - that is, matrix element A( i, j ) is stored at
*           array position a( i, j ) -, with the strict lower
*           (uplo = 'U') or the strictly upper (uplo = 'L') triangle set
*           to some large value (NOT_USED).
*           If frmt = 'B' then A is returned in banded format in the
*           leading (b+1)-times-n subarray - that is, matrix element
*           A( i, j ) is stored at array position a( b+1+i-j, j ) if
*           uplo = 'U', and at position a( 1+i-j, j ) if uplo = 'L'.
*           In any case, the remainder of the first n columns is set
*           to the same large value mentioned above.
*
*   lda     (in) integer
*           The leading dimension of the array a.
*           lda >= n,   if frmt = 'F'
*               >= b+1, if frmt = 'B'.
*
*   sza     (in) integer
*           The total size (in elements) of the array (not matrix!) a.
*           sza >= lda * n.
*
*   jobu    (in) character*1
*           Set the matrix U to identity?
*           jobu = 'U' : Yes.
*                = 'N' : No.
*
*   u       (out) real array, dimension ( ldu, n )
*           If jobu = 'U' then the leading n-by-n subarray is set to
*           an n-by-n identity matrix.
*           If jobu = 'N' then u is not accessed.
*
*   ldu     (in) integer
*           The leading dimension of the array u.
*           ldu >= n.
*
*   szu     (in) integer
*           The total size (in elements) of the array (not matrix!) u.
*           szu >= ldu * n, if jobu = 'U'
*               >= 0      , if jobu = 'N'.
*
*   notusd  (in) real
*           Value for marking elements that should not be used.
*           Typically a rather large number.
*
*   work    (workspace) real array, dimension ( szwork )
*
*   szwork  (in) integer
*           The size of the array work.
*           szwork >= 2 * n.
*
*   info    (out) integer
*           On exit, info indicates the consistence of the arguments.
*           =   0 : All arguments are OK.
*           =  -1 : mtype is none of 'S', 'R', 'I' (upper/lower case).
*           =  -3 : szevls is too small ( < n ).
*           =  -4 : diam is out of range (negative).
*           =  -5 : frmt is none of 'F', 'B' (upper or lower case).
*           =  -6 : uplo is none of 'U', 'L' (upper or lower case).
*           =  -7 : n is out of range ( < 0 ).
*           =  -8 : b is out of range ( < 0 or > max( n-1, 0 ) ).
*           = -10 : lda is too small ( < n, if full, < b+1, if banded).
*           = -11 : sza is too small ( < n*lda ).
*           = -12 : jobu is none of 'U', 'N' (upper or lower case).
*           = -14 : ldu is too small ( < n ).
*           = -15 : szu is too small ( < n*ldu, if U is required).
*           = -18 : szwork is too small (see above).
*
* ----------------------------------------------------------------------
*
* Local constants and variables:
*
        REAL                  ZERO, HALF, ONE, TWO
        PARAMETER             ( ZERO = 0.0E0, HALF = 0.5E0,
     >                          ONE = 1.0E0, TWO = 2.0E0 )
*
        LOGICAL               SPEC, ISDA, FULL, UPPER, NEEDU
        INTEGER               ISEED( 4 )
        REAL                  SDIAM
*
*   spec    matrix with prescribed spectrum ?
*   isda    ISDA-type matrix or random matrix ?
*   full    store the matrix in full or banded format ?
*   upper   use upper or lower triangle ?
*   needu   update I or not ?
*   iseed   seed for the random number generator
*   sdiam   scaled diameter
*
        INTEGER               D, I, J, K, Q
        REAL                  C, S, FILL, A11, A12, A21, A22, TMP
*
*   d       the number of the new diagonal
*   c       cosine of the rotation angle
*   s       sine of the rotation angle
*   fill    new fill-in element
*
* Routines called:
*
        LOGICAL               LSAME
        EXTERNAL              LSAME
*
*   lsame   case-insensitive character matching (BLAS)
*
        REAL                  R1MACH, SLARAN
        EXTERNAL              SCOPY, SLAGSY, R1MACH, SLARAN, SLARTG,
     >                        SLASET
*
*   scopy   vector copy (BLAS)
*   slagsy  generate symmetric matrix (LAPACK TMG)
*   slamch  inquire machine parameters (LAPACK)
*   slaran  random number generator (LAPACK TMG)
*   slartg  generate rotation (LAPACK)
*   slaset  initialize diagonal matrix (LAPACK)
*
        INTRINSIC             MAX, MIN
*
* ----------------------------------------------------------------------
*
*          --- sanity check for the arguments ---
*
        SPEC = LSAME( MTYPE, 'Spectrum' )
        ISDA = LSAME( MTYPE, 'Isda' )
        FULL = LSAME( FRMT, 'Full' )
        UPPER = LSAME( UPLO, 'Upper' )
        NEEDU = LSAME( JOBU, 'Update' )
*
        IF ( .NOT. ( ISDA .OR. SPEC .OR. LSAME( MTYPE, 'Random' ) ) )
     >  THEN
          INFO = -1
        ELSEIF ( SZEVLS .LT. N ) THEN
          INFO = -3
        ELSEIF ( DIAM .LT. ZERO ) THEN
          INFO = -4
        ELSEIF ( .NOT. ( FULL .OR. LSAME( FRMT, 'Banded' ) ) ) THEN
          INFO = -5
        ELSEIF ( .NOT. ( UPPER .OR. LSAME( UPLO, 'Lower' ) ) ) THEN
          INFO = -6
        ELSEIF ( N .LT. 0 ) THEN
          INFO = -7
        ELSEIF ( ( B .LT. 0 ) .OR. ( B .GT. MAX( N-1, 0 ) ) ) THEN
          INFO = -8
        ELSEIF ( ( FULL .AND. ( LDA .LT. N ) )
     >           .OR. ( LDA .LT. B+1 ) ) THEN
          INFO = -10
        ELSEIF ( SZA .LT. N*LDA ) THEN
          INFO = -11
        ELSEIF ( .NOT. ( NEEDU .OR. LSAME( JOBU, 'NoUpdate' ) ) ) THEN
          INFO = -12
        ELSEIF ( LDU .LT. N ) THEN
          INFO = -14
        ELSEIF ( NEEDU .AND. ( SZU .LT. N*LDU ) ) THEN
          INFO = -15
        ELSEIF ( SZWORK .LT. 2*N ) THEN
          INFO = -18
        ELSE
          INFO = 0
        ENDIF
        IF ( INFO .NE. 0 )     GOTO 9999
*
*          --- check for quick return ---
*
        IF ( N .EQ. 0 )     GOTO 9999
*
        IF ( .NOT. SPEC ) THEN
*
*            --- set up eigenvalue distribution ---
*
          ISEED( 1 ) = 4001
          ISEED( 2 ) = 3001
          ISEED( 3 ) = 2001
          ISEED( 4 ) = 1001
*
          SDIAM = TWO * DIAM * R1MACH(4)
          DO 1000 I = 1, N
            IF ( ISDA ) THEN
              EVLS( I ) = SDIAM * ( SLARAN( ISEED ) - HALF )
              IF ( SLARAN( ISEED ) .LT. HALF ) THEN
                EVLS( I ) = EVLS( I ) + ONE
              ENDIF
            ELSE
              EVLS( I ) = TWO * ( SLARAN( ISEED ) - HALF )
            ENDIF
 1000     CONTINUE
        ENDIF
*
*          --- initialize random number generator ---
*
        ISEED( 1 ) = 1001
        ISEED( 2 ) = 2001
        ISEED( 3 ) = 3001
        ISEED( 4 ) = 4001
*
        IF ( FULL ) THEN
*
*  - - - - - - - initialize matrix in full storage - - - - - - - - - - -
*
          CALL SLAGSY( N, B, EVLS, A, LDA, ISEED, WORK, INFO )
*
          IF ( UPPER ) THEN
*
*              --- mark lower triangle "not used" ---
*
            DO 2100 J = 1, N-1
              DO 2010 I = J+1, N
                A( I, J ) = NOTUSD
 2010         CONTINUE
 2100       CONTINUE
          ELSE
*
*              --- mark upper triangle "not used" ---
*
            DO 2300 J = 2, N
              DO 2210 I = 1, J-1
                A( I, J ) = NOTUSD
 2210         CONTINUE
 2300       CONTINUE
          ENDIF
*
        ELSE
*
* - - - - - - - - initialize matrix in banded storage - - - - - - - - -
*
* - - - - - - - - first initialize in lower packed storage - - - - - - -
*
          CALL SCOPY( N, EVLS, 1, A( 1, 1 ), LDA )
*
          DO 3500 D = 1, B
*
*              --- add d-th diagonal ---
*
            DO 3300 I = N-1, 1, -1
*
*                --- to generate an new element at position
*                    ( i+1, i-d ), ...                      ---
*
*                --- ... first determine a random rotation, ... ---
*
              A11 = SLARAN( ISEED )
              A12 = SLARAN( ISEED )
              CALL SLARTG( A11, A12, C, S, TMP )
*
*                --- ... apply it to rows i and i+1, ... ---
*
              DO 3010 Q = MAX( I-D+1, 1 ), I-1
                IF ( Q .EQ. ( I-D+1 ) ) THEN
                  A( D+1, Q ) = -S * A( D, Q )
                  A( D, Q ) = C * A( D, Q )
                ELSE
                  TMP = C * A( I-Q+2, Q ) - S * A( I-Q+1, Q )
                  A( I-Q+1, Q ) = C * A( I-Q+1, Q ) + S * A( I-Q+2, Q )
                  A( I-Q+2, Q ) = TMP
                ENDIF
 3010         CONTINUE
*
*                --- ... apply it to A( i:i+1, i:i+1 ), ... ---
*
              IF ( D .EQ. 1 ) THEN
                A11 = C * A( 1, I )
                A21 = -S * A( 1, I )
                A12 = S * A( 1, I+1 )
                A22 = C * A( 1, I+1 )
              ELSE
                A11 = C * A( 1, I ) + S * A( 2, I )
                A21 = C * A( 2, I ) - S * A( 1, I )
                A12 = C * A( 2, I ) + S * A( 1, I+1 )
                A22 = C * A( 1, I+1 ) - S * A( 2, I )
              ENDIF
              A( 1, I ) = C * A11 + S * A12
              A( 2, I ) = C * A21 + S * A22
              A( 1, I+1 ) = C * A22 - S * A21
*
*                --- ... apply it to columns i and i+1, ... ---
*
              DO 3020 Q = I+2, MIN( I+D, N )
                TMP = C * A( Q-I, I+1 ) - S * A( Q-I+1, I )
                A( Q-I+1, I ) = C * A( Q-I+1, I ) + S * A( Q-I, I+1 )
                A( Q-I, I+1 ) = TMP
 3020         CONTINUE
              IF ( ( I+D+1 ) .LE. N ) THEN
                FILL = S * A( D+1, I+1 )
                A( D+1, I+1 ) = C * A( D+1, I+1 )
              ENDIF
*
*                --- ... and chase the fill-in element down the band ---
*
              DO 3200 K = I+D, N-1, D
*
                CALL SLARTG( A( D+1, K-D ), FILL, C, S, TMP )
                A( D+1, K-D ) = TMP
*
                DO 3110 Q = K-D+1, K-1
                  TMP = C * A( K-Q+2, Q ) - S * A( K-Q+1, Q )
                  A( K-Q+1, Q ) = C * A( K-Q+1, Q ) + S * A( K-Q+2, Q )
                  A( K-Q+2, Q ) = TMP
 3110           CONTINUE
*
                A11 = C * A( 1, K ) + S * A( 2, K )
                A21 = C * A( 2, K ) - S * A( 1, K )
                A12 = C * A( 2, K ) + S * A( 1, K+1 )
                A22 = C * A( 1, K+1 ) - S * A( 2, K )
                A( 1, K ) = C * A11 + S * A12
                A( 2, K ) = C * A21 + S * A22
                A( 1, K+1 ) = C * A22 - S * A21
*
                DO 3120 Q = K+2, MIN( K+D, N )
                  TMP = C * A( Q-K, K+1 ) - S * A( Q-K+1, K )
                  A( Q-K+1, K ) = C * A( Q-K+1, K ) + S * A( Q-K, K+1 )
                  A( Q-K, K+1 ) = TMP
 3120           CONTINUE
                IF ( ( K+D+1 ) .LE. N ) THEN
                  FILL = S * A( D+1, K+1 )
                  A( D+1, K+1 ) = C * A( D+1, K+1 )
                ENDIF
 3200         CONTINUE
*
 3300       CONTINUE
*
*              --- mark the last d entries of the row "not used" ---
*
            DO 3410 I = N-D+1, N
              A( D+1, I ) = NOTUSD
 3410       CONTINUE
*
 3500     CONTINUE
*
* - - - - - - - - end : first initialize in lower packed storage - - - -
*
          IF ( UPPER ) THEN
*
*              --- convert from lower to upper banded storage ---
*
            DO 3700 I = 0, B/2
              CALL SCOPY( N-B+I, A( B+1-I, 1 ), LDA, WORK, 1 )
              DO 3610 J = 1, N-I
                A( B+1-I, I+J ) = A( I+1, J )
 3610         CONTINUE
              DO 3620 K = 1, I
                A( B+1-I, K ) = NOTUSD
 3620         CONTINUE
              CALL SCOPY( N-B+I, WORK, 1, A( I+1, B+1-I ), LDA )
              DO 3630 K = 1, B-I
                A( I+1, K ) = NOTUSD
 3630         CONTINUE
 3700       CONTINUE
          ENDIF
*
* - - - - - - - - end : initialize matrix in banded storage - - - - - -
*
        ENDIF
*
*          --- generate identity matrix in U ---
*
        IF ( NEEDU ) THEN
          CALL SLASET( 'All', N, N, ZERO, ONE, U, LDU )
        ENDIF
*
 9999   RETURN
        END
*
* **********************************************************************
*
        SUBROUTINE SSYNCK( EVLS, FRMT1, UPLO1, N, B1, A1, LDA1, SZA1,
     >                     FRMT2, B2, A2, LDA2, D, E,
     >                     JOBU, U, LDU,
     >                     ORTHU, SORTHU, RES, SRES, DEVLS, SDEVLS,
     >                     NOTUSD, WORK, SZWORK, INFO )
*
* Description:
*
*   This routine performs numerical tests after an orthogonal
*   transformation
*                 T
*      A --> B = U  * A * U.
*   A is rebuilt by calling ssyini, B may be given in banded or
*   tridiagonal format.
*   Depending on whether U was build explicitly during the reduction,
*   either
*                           T
*      the orthogonality   U  * U - I      of U
*      and the residual    U * B - A * U
*   or
*      the difference between A's and B's eigenvalues
*   is checked.
*
* Author: Bruno Lang
*         Aachen University of Technology
*         na.blang@na-net.ornl.gov
* Date: May 17, 2000
* Version: SBR Toolbox, Rev. 1.4.1
*
* ----------------------------------------------------------------------
*
* Parameters:
*
        CHARACTER*1           FRMT1, UPLO1, FRMT2, JOBU
        INTEGER               N, B1, LDA1, SZA1, B2, LDA2, LDU,
     >                        SZWORK, INFO
        REAL                  ORTHU, SORTHU, RES, SRES, DEVLS, SDEVLS,
     >                        NOTUSD
        REAL                  EVLS( * ), A1( LDA1, * ), A2( LDA2, * ),
     >                        D( * ), E( * ), U( LDU, * ), WORK( * )
*
*   evls    (in) real array, dimension ( n )
*           The eigenvalues of the matrix A, as returned from ssyini.
*
*   frmt1   (in) character*1
*           The storage format of the original matrix A.
*           frmt1 = 'F' : (Symmetric) full storage.
*                 = 'B' : (Symmetric) banded storage.
*
*   uplo1   (in) character*1
*           Which triangle of the matrix A (and B) is used?
*           uplo1 = 'U' : Use the upper triangle.
*                 = 'L' : Use the lower triangle.
*           In banded storage (frmt1 = 'B'), uplo1 also determines
*           whether the upper or lower banded storage scheme is used.
*
*   n       (in) integer
*           The dimension of the matrix.
*           n >= 0.
*
*   b1      (in) integer
*           The semibandwidth (i.e., the number of subdiagonals) of the
*           original matrix A.
*           0 <= b1 <= max( n-1, 0 ).
*
*   a1      (out) real array, dimension ( lda1, n )
*           a1 is used to build a copy of the original n-by-n symmetric
*           matrix A with semibandwidth b1.
*           If frmt = 'F' then A is returned in the leading n-by-n
*           subarray - that is, matrix element A( i, j ) is stored at
*           array position a1( i, j ) -, with the strict lower
*           (uplo = 'U') or the strictly upper (uplo = 'L') triangle set
*           to some large value (NOT_USED).
*           If frmt = 'B' then A is returned in banded format in the
*           leading (b1+1)-times-n subarray - that is, matrix element
*           A( i, j ) is stored at array position a1( b1+1+i-j, j ) if
*           uplo = 'U', and at position a1( 1+i-j, j ) if uplo = 'L'.
*
*   lda1    (in) integer
*           The leading dimension of the array a1.
*           lda >= n,    if frmt1 = 'F'
*               >= b1+1, if frmt1 = 'B'.
*
*   sza1    (in) integer
*           The total size (in elements) of the array (not matrix!) a1.
*           sza1 >= lda1 * n.
*
*   frmt2   (in) character*1
*           The storage format of the reduced matrix B.
*           frmt1 = 'B' : (Symmetric) banded storage.
*                 = 'T' : (Symmetric) tridiagonal in arrays d and e.
*
*   b2      (in) integer
*           The semibandwidth (i.e., the number of subdiagonals) of the
*           reduced matrix B.
*           0 <= b2 <= max( n-1, 0 ).
*
*   a2      (in) real array, dimension ( lda2, n )
*           If frmt2 = 'B' then the leading (b2+1)-times-n subarray of
*           a2 holds the reduced symmetric matrix B with semibandwidth
*           b2 in banded format - that is, matrix element A( i, j ) is
*           stored at array position a2( b2+1+i-j, j ) if uplo = 'U',
*           and at position a2( 1+i-j, j ) if uplo = 'L'.
*           If frmt2 = 'T' then a2 is not referenced.
*
*   lda2    (in) integer
*           The leading dimension of the array a2.
*           lda2 >= b2+1.
*
*   d       (in) real array, dimension ( n )
*           If frmt2 = 'T' then the first n elements of d hold the
*           main diagonal of the tridiagonal matrix B.
*           If frmt2 = 'B' then d is not referenced.
*
*   e       (in) real array, dimension ( n-1 )
*           If frmt2 = 'T' then the first n-1 elements of e hold the
*           subdiagonal (superdiagonal) of the tridiagonal matrix B.
*           If frmt2 = 'B' then e is not referenced.
*
*   jobu    (in) character*1
*           Does U contain the accumulated transformations ?
*           jobu = 'U' : Yes.
*                = 'N' : No.
*
*   u       (in) real array, dimension ( ldu, n )
*           If jobu = 'U' then the leading n-by-n subarray contains
*           the accumulated transformations from the reduction (an
*           orthogonal matrix).
*           If jobu = 'N' then u is not accessed.
*
*   ldu     (in) integer
*           The leading dimension of the array u.
*           ldu >= n.
*
*   orthu   (out) real
*           The orthogonality error : F-norm( U''*U - I ).
*           If jobu = 'N' then orthu is set to 0.
*
*   sorthu  (out) real
*           The scaled orthogonality error : orthu / ( n * macheps).
*           If jobu = 'N' then sorthu is set to 0.
*
*   res     (out) real
*           The residual : F-norm( U*B - A*U ).
*           If jobu = 'N' then res is set to 0.
*
*   sres    (out) real
*           The scaled residual : res / ( n * macheps * F-norm( A ) ).
*           If jobu = 'N' then sres is set to 0.
*
*   devls   (out) real
*           The change of the eigenvalues :
*              2-norm( sort(evls) - sort(spec(B)) ).
*           If jobu = 'U' then devls is set to 0.
*
*   sdevls  (out) real
*           The scaled eigenvalue changes :
*              devls / ( n * macheps * F-norm( A ) ).
*           If jobu = 'U' then sdevls is set to 0.
*
*   notusd  (in) real
*           Value for marking elements that should not be used.
*           Typically a rather large number.
*
*   work    (workspace) real array, dimension ( szwork )
*
*   szwork  (in) integer
*           The length of the array work.
*           szwork >= max( 3*n-2, 2*n ).
*
*   info    (out) integer
*           On exit, info indicates the consistence of the arguments.
*           =   0 : All arguments are OK.
*           =  -2 : frmt1 is none of 'F', 'B' (upper/lower case).
*           =  -3 : uplo1 is none of 'U', 'L' (upper/lower case).
*           =  -4 : n is out of range ( < 0 ).
*           =  -5 : b1 is out of range ( < 0 or > max( n-1, 0 ) ).
*           =  -7 : lda1 too small ( < n, if full, < b1+1, if banded).
*           =  -8 : sza1 is too small ( < n * lda1 ).
*           =  -9 : frmt2 is none of 'B', 'T' (upper/lower case).
*           = -10 : b2 is out of range ( < 0 or > max( n-1, 0 ),
*                   or > 1 if frmt2 = 'T' ).
*           = -12 : lda2 is too small ( < b2+1 ).
*           = -15 : jobu is none of 'U', 'N' (upper/lower case).
*           = -17 : ldu is too small ( < n ).
*           = -26 : szwork is too small (see above).
*
* ----------------------------------------------------------------------
*
* Local constants and variables:
*
        REAL                  ZERO, ONE
        PARAMETER             ( ZERO = 0.0E0, ONE = 1.0E0 )
*
        LOGICAL               FULL1, UPPER, TRI2, HAVEU
        REAL                  EPS, NORMA
        INTEGER               NCOLS, COLS, J0, J, IX
*
*   full1   build A in full storage ?
*   upper   use upper triangle ?
*   tri2    is B given in tridiagonal storage ?
*   haveu   is the accumulated matrix U available ?
*   eps     machine epsilon
*   norma   the Frobenius norm of the matrix A
*   ncols   maximum width of a block column that fits into workspace
*   cols    width of current block column
*   j0      first column of the current block column
*   j       current column
*   ix      points to the position of the current column in workspace
*
* Routines called:
*
        LOGICAL               LSAME
        EXTERNAL              LSAME
*
*   lsame   case-insensitive character matching (BLAS)
*
        REAL                  R1MACH, SLANGE, SNRM2
        EXTERNAL              SAXPY, SCOPY, SGEMM, SGEMV, R1MACH,
     >                        SLANGE, SLASRT, SNRM2, SSBEV, SSBMV,
     >                        SSCAL, SSTEV, SSYINI, SSYMM
*
*   saxpy   add multiple of a vector to another one (BLAS)
*   scopy   vector copy (BLAS)
*   sgemm   matrix-matrix product (BLAS)
*   sgemv   matrix-vector product (BLAS)
*   slamch  determine machine parameters (LAPACK)
*   slange  matrix norm (LAPACK)
*   slasrt  sort a vector (LAPACK)
*   snrm2   2-norm of a vector (BLAS)
*   ssbev   compute eigensystem of symmetric banded matrix (LAPACK)
*   ssbmv   banded symmetri matrix-vector product (BLAS)
*   sscal   scale a vector (BLAS)
*   sstev   compute eigensystem of symmetric tridiagonal matrix (LAPACK)
*   ssyini  initialize symmetric matrix (this file)
*   ssymm   symmetric matrix-matrix product (BLAS)
*
        INTRINSIC             MAX, MIN, SQRT
*
* ----------------------------------------------------------------------
*
*          --- sanity check for the arguments ---
*
        FULL1 = LSAME( FRMT1, 'Full' )
        UPPER = LSAME( UPLO1, 'Upper' )
        TRI2 = LSAME( FRMT2, 'Tridiagonal' )
        HAVEU = LSAME( JOBU, 'U' )
*
        IF ( .NOT. ( FULL1 .OR. LSAME( FRMT1, 'Banded' ) ) ) THEN
          INFO = -2
        ELSEIF ( .NOT. ( UPPER .OR. LSAME( UPLO1, 'Lower' ) ) ) THEN
          INFO = -3
        ELSEIF ( N .LT. 0 ) THEN
          INFO = -4
        ELSEIF ( ( B1 .LT. 0 ) .OR. ( B1 .GT. MAX( N-1, 0 ) ) ) THEN
          INFO = -5
        ELSEIF ( ( FULL1 .AND. ( LDA1 .LT. N ) ) .OR.
     >           ( LDA1 .LT. B1+1 ) ) THEN
          INFO = -7
        ELSEIF ( SZA1 .LT. N*LDA1 ) THEN
          INFO = -8
        ELSEIF ( .NOT. ( TRI2 .OR. LSAME( FRMT2, 'Banded' ) ) ) THEN
          INFO = -9
        ELSEIF ( ( B2 .LT. 0 ) .OR.
     >           ( ( N .GT. 0 ) .AND.
     >             ( ( B2 .GT. N-1 ) .OR.
     >               ( TRI2 .AND. ( B2 .GT. 1 ) ) ) ) ) THEN
          INFO = -10
        ELSEIF ( LDA2 .LT. B2+1 ) THEN
          INFO = -12
        ELSEIF ( .NOT. ( HAVEU .OR. LSAME( JOBU, 'NoU' ) ) ) THEN
          INFO = -15
        ELSEIF ( LDU .LT. N ) THEN
          INFO = -17
        ELSEIF ( SZWORK .LT. MAX( 3*N-2, 2*N ) ) THEN
          INFO = -26
        ELSE
          INFO = 0
        ENDIF
        IF ( INFO .NE. 0 )     GOTO 9999
*
*          --- check for quick return ---
*
        IF ( N .EQ. 0 )     GOTO 9999
*
*          --- determine machine precision ---
*
        EPS = R1MACH(4)
*                                                     T
*  - - - - - - compute orthogonality error = F-Norm( U  * U - I ) - - -
*
        IF ( LSAME( JOBU, 'Update' ) ) THEN
          ORTHU = ZERO
*                                     T
*            --- how many columns of U  * U - I fit into workspace ? ---
*
          NCOLS = SZWORK / N
          DO 1100 J0 = 1, N, NCOLS
            COLS = MIN( NCOLS, N-J0+1 )
*                                                     T
*              --- compute n-by-cols block column of U  * U ---
*
            CALL SGEMM( 'Transpose', 'NoTranspose', N, COLS, N,
     >                  ONE, U, LDU, U( 1, J0 ), LDU, ZERO, WORK, N )
*
*              --- subtract corresponding parts of I ---
*
            DO 1010 J = J0, J0+COLS-1
              WORK( (J-J0)*N+J ) = WORK( (J-J0)*N+J ) - ONE
 1010       CONTINUE
            ORTHU = ORTHU + SLANGE( 'Frobenius', N, COLS, WORK, N,
     >                              WORK )**2
 1100     CONTINUE
          ORTHU = SQRT( ORTHU )
          SORTHU = ORTHU / ( N * EPS )
        ELSE
*
*            --- no orthogonality information available ---
*
          ORTHU = - ONE
          SORTHU = - ONE
        ENDIF
*                                                           T
*  - - - - - - end : compute orthogonality error = F-Norm( U  * U - I )
*
        NORMA = SNRM2( N, EVLS, 1 )
*
        IF ( LSAME( JOBU, 'Update' ) ) THEN
*
*            --- do not compute changes in eigenvalues ---
*
          DEVLS = - ONE
          SDEVLS = - ONE
*
*  - - - - - - - compute residual error = F-Norm( U * B - A * U ) - - -
*
*            --- build another copy of the original matrix A ---
*
          CALL SSYINI( 'Spectrum', EVLS, N, ZERO,
     >                 FRMT1, UPLO1, N, B1, A1, LDA1, SZA1,
     >                 'NoU', WORK, N, 0,
     >                 NOTUSD, WORK, SZWORK, INFO )
*
          RES = ZERO
*
*            --- how many cols of U * B - A * U fit into workspace ? ---
*
          NCOLS = SZWORK / N
*
          DO 2000 J0 = 1, N, NCOLS
            COLS = MIN( NCOLS, N-J0+1 )
*
*              --- build the next cols columns of U * B ---
*
            IF ( LSAME( FRMT2, 'Tridiagonal' ) ) THEN
              DO 2010 J = J0, J0+COLS-1
                IX = ( J - J0 ) * N + 1
                CALL SCOPY( N, U( 1, J ), 1, WORK( IX ), 1 )
                CALL SSCAL( N, D( J ), WORK( IX ), 1 )
                IF ( J .GT. 1 ) THEN
                  CALL SAXPY( N, E( J-1 ), U( 1, J-1 ), 1,
     >                        WORK( IX ), 1 )
                ENDIF
                IF ( J .LT. N ) THEN
                  CALL SAXPY( N, E( J ), U( 1, J+1 ), 1, WORK( IX ), 1 )
                ENDIF
 2010         CONTINUE
            ELSE
              DO 2020 J = J0, J0+COLS-1
                IX = ( J - J0 ) * N + 1
                CALL SGEMV( 'NoTranspose', N, MIN( B2+1, N-J+1 ),
     >                      ONE, U( 1, J ), LDU, A2( 1, J ), 1,
     >                      ZERO, WORK( IX ), 1 )
                IF ( B2 .GT. 0 ) THEN
                  IF ( J .GT. B2 ) THEN
                    CALL SGEMV( 'NoTranspose', N, B2,
     >                          ONE, U( 1, J-B2 ), LDU,
     >                          A2( B2+1, J-B2 ), LDA2 - 1,
     >                          ONE, WORK( IX ), 1 )
                  ELSEIF ( J .GT. 1 ) THEN
                    CALL SGEMV( 'NoTranspose', N, J - 1,
     >                          ONE, U( 1, 1 ), LDU,
     >                          A2( J, 1 ), LDA2 - 1,
     >                          ONE, WORK( IX ), 1 )
                  ENDIF
                ENDIF
 2020         CONTINUE
            ENDIF
*
*              --- subtract corresponding columns of A * U ---
*
            IF ( LSAME( FRMT1, 'Full' ) ) THEN
              CALL SSYMM( 'Left', UPLO1, N, COLS,
     >                   -ONE, A1, LDA1, U( 1, J0 ), LDU,
     >                   ONE, WORK, N )
            ELSE
              DO 2110 J = J0, J0+COLS-1
                CALL SSBMV( UPLO1, N, B1, -ONE, A1, LDA1, U( 1, J ), 1,
     >                      ONE, WORK( (J-J0)*N+1 ), 1 )
 2110         CONTINUE
            ENDIF
*
            RES = RES + SLANGE( 'Frobenius', N, COLS, WORK, N, WORK )**2
 2000     CONTINUE
*
          RES = SQRT( RES )
          SRES = RES / ( N * NORMA * EPS )
*
*  - - - - - - - end : compute residual error = F-Norm( U * B - A * U )
*
        ELSE
*
*            --- no residual information available ---
*
          RES = - ONE
          SRES = - ONE
*
*  - - - - - - - check if eigenvalues have changed - - - - - - - - - - -
*
*            --- compute eigenvalues of the tridiagonal or banded
*                matrix                                           ---
*
          IF ( LSAME( FRMT2, 'Tridiagonal' ) ) THEN
            CALL SSTEV( 'NoEigenvectors', N, D, E, WORK, N, WORK, INFO )
          ELSE
            CALL SSBEV( 'NoEigenvectors', 'Lower', N, B2, A2, LDA2,
     >                  D, WORK, N, WORK, INFO )
          ENDIF
*
*            --- compare with (sorted) original eigenvalues ---
*
          IF ( INFO .EQ. 0 ) THEN
            CALL SLASRT( 'Increasing', N, EVLS, INFO )
            CALL SAXPY( N, -ONE, EVLS, 1, D, 1 )
            DEVLS = SNRM2( N, D, 1 )
            SDEVLS = DEVLS / ( N * NORMA * EPS )
          ELSE
            DEVLS = - ONE
            SDEVLS = - ONE
          ENDIF
*
*  - - - - - - - end : check if eigenvalues have changed - - - - - - - -
*
        ENDIF
*
 9999   RETURN
        END
*
* **********************************************************************
*
        LOGICAL FUNCTION LXSAME( S1, S2 )
*
* Description:
*
*   Checks if two strings agree (ignoring case).
*
* Author: Bruno Lang
*         Aachen University of Technology
*         na.blang@na-net.ornl.gov
* Date: May 17, 2000
* Version: SBR Toolbox, Rev. 1.4.1
*
* Parameters:
*
        CHARACTER*( * )       S1, S2
*
*   s1      (in) character*( * )
*           The first string.
*
*   s2      (in) character*( * )
*           The second string.
*
* Local constants and variables:
*
        INTEGER               I
*
* Routines called:
*
        LOGICAL               LSAME
        EXTERNAL              LSAME
*
*   lsame   case-insensitive character matching (BLAS)
*
        INTRINSIC             LEN
*
* ----------------------------------------------------------------------
*
        LXSAME = ( LEN( S1 ) .EQ. LEN( S2 ) )
        IF ( LXSAME ) THEN
          DO 100 I = 1, LEN( S1 )
            LXSAME = LXSAME .AND. LSAME( S1( I:I ), S2( I:I ) )
  100     CONTINUE
        ENDIF
*
        RETURN
        END
*
