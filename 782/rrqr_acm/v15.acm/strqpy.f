      SUBROUTINE STRQPY( JOB, M, N, K, A, LDA, C, LDC, JPVT, IRCOND,
     $                   ORCOND, RANK, SVLUES, WORK, LWORK, INFO )
*
*     This code is part of a release of the package for computing
*     rank-revealing QR Factorizations written by:
*     ==================================================================
*     Christian H. Bischof        and   Gregorio Quintana-Orti
*     Math. and Comp. Sci. Div.         Departamento de Informatica
*     Argonne National Lab.             Universidad Jaime I
*     Argonne, IL 60439                 Campus P. Roja, 12071 Castellon
*     USA                               Spain
*     bischof@mcs.anl.gov               gquintan@inf.uji.es
*     ==================================================================
*     $Revision: 1.84 $
*     $Date: 96/12/30 16:59:19 $
*
*     .. Scalar Arguments ..
      INTEGER            JOB, M, N, K, LDA, LDC, RANK, LWORK, INFO
      REAL               IRCOND, ORCOND
*     ..
*     .. Array Arguments ..
      REAL               A( LDA, * ), C( LDC, * ), SVLUES( 4 )
      REAL               WORK( LWORK )
      INTEGER            JPVT( * )
*     ..
*
*  Purpose
*  =======
*
*  STRQPY detects the right rank for upper triangular matrix A.
*  The algorithm used here is an version of Pan and Tang's RRQR
*  algorithm number 3.
*  This algorithm is applied to matrix A until the right rank is
*  obtained. If the input ordering of matrix A is not accepted, the
*  matrix will be permuted and retriangularized until the rank is
*  revealed.
*
*  Arguments
*  =========
*
*  JOB     (input) INTEGER
*          The job to do:
*          = 1: The orthogonal transformations needed in the
*               triangularization are only applied to matrix A.
*               Thus, matrix C is not updated.
*          = 2: The same orthogonal transformations needed in the
*               triangularization of matrix A are applied to
*               matrix C from the left.
*               That is, if Q'*A*P=R, then C := Q'*C.
*               In this case, matrix C is m-by-k.
*          = 3: The transpose of the orthogonal transformations needed
*               in the triangularization of matrix A are applied
*               to matrix C from the right.
*               That is, if Q'*A*P=R, then C := C*Q.
*               In this case, matrix C is k-by-m.
*          In these three cases, the permutations are always stored
*          in vector JPVT.
*
*  M       (input) INTEGER
*          The number of rows of matrices A. M >= 0.
*          If JOB=2, M is the number of rows of matrix C.
*          If JOB=3, M is the number of columns of matrix C.
*
*  N       (input) INTEGER
*          The number of columns of matrix A.  N >= 0.
*
*  K       (input) INTEGER
*          It defines the dimension of matrix C. K >= 0.
*          If JOB=2, K is the number of columns of matrix C.
*          If JOB=3, K is the number of rows of matrix C.
*
*  A       (input/output) REAL array, dimension (LDA,N)
*          On entry, the m by n matrix A.
*          On exit, the upper triangle of the array contains the
*          min(m,n) by n upper trapezoidal matrix R; the lower triangle
*          array is filled with zeros.
*
*  LDA     (input) INTEGER
*          The leading dimension of array A. LDA >= max(1,M).
*
*  C       (input/output) REAL array, dimension
*                ( LDC, K ) if JOB=2.
*                ( LDC, M ) if JOB=3.
*          If argument JOB asks, all the orthogonal transformations
*          applied to matrix A are also applied to matrix C.
*
*  LDC     (input) INTEGER
*          The leading dimension of array C.
*          If JOB=2, then LDC >= MAX(1,M).
*          If JOB=3, then LDC >= MAX(1,K).
*
*  JPVT    (input/output) INTEGER array, dimension ( N )
*          If JPVT(I) = K, then the Ith column of the permuted
*          A was the Kth column of the original A (just before
*          the preprocessing). If a permutation occurs, JPVT will
*          be updated correctly.
*          JPVT(1:RANK) contains the indices of the columns considered
*          linearly independent.
*          JPVT(RANK+1:N) contains the indices of the columns considered
*          linearly dependent from the previous ones.
*
*  IRCOND  (input) REAL
*          1/IRCOND specifies an upper bound on the condition number
*          of R11. If IRCOND == 0, IRCOND = machine precision is chosen
*          as default. IRCOND must be >= 0.
*
*  ORCOND  (output) REAL
*          1/ORCOND is an estimate for the condition number of R11.
*
*  RANK    (output) INTEGER
*          An estimate of the rank offered by this algorithm.
*          0 <= RANK <= MIN(M,N).
*
*  SVLUES  (output) REAL array, dimension (4)
*          On exit, SVLUES contains estimates of some of the
*          singular values of the triangular factor R.
*          SVLUES(1): largest singular value of R(1:RANK,1:RANK)
*          SVLUES(2): smallest singular value of R(1:RANK,1:RANK)
*          SVLUES(3): smallest singular value of R(1:RANK+1,1:RANK+1)
*          SVLUES(4): smallest singular value of R
*          If the triangular factorization is a rank-revealing one
*          (which will be the case if the leading columns were well-
*          conditioned), then SVLUES(1) will also be an estimate
*          for the largest singular value of A, SVLUES(2) and SVLUES(3)
*          will be estimates for the RANK-th and (RANK+1)-st singular
*          value of A, and SVLUES(4) will be an estimate for the
*          smallest singular value of A.
*          By examining these values, one can confirm that the rank is
*          well defined with respect to the threshold chosen.
*
*  WORK    (workspace) REAL array, dimension ( LWORK )
*
*  LWORK   (input) INTEGER
*          The dimension of array WORK. LWORK >= N+3*MN, where
*          MN = min(M,N).
*
*  INFO    (output) INTEGER
*          = 0: Successful exit.
*          < 0: If INFO = -i, the i-th argument had an illegal value
*          > 0: Problems in the computation of the rank.
*                   1: Exceeded the allowed maximum number of steps.
*                   2: Rank not well defined.
*          In adition, vector SVLUES tell if rank is not well defined.
*          When INFO.NE.0, the contents of ORCOND may be not the right
*          one.
*
*
*  ===================================================================
*
*     .. Parameters ..
      INTEGER            INB
      REAL               ZERO
      PARAMETER          ( INB = 1, ZERO = 0.0E+0 )
*
*     Indices into the 'svlues' array.
*
      INTEGER            IMAX, IBEFOR, IAFTER, IMIN
      PARAMETER          ( IMAX = 1, IBEFOR = 2, IAFTER = 3, IMIN = 4 )
*     ..
*     ..
*     .. Common Block ..
      INTEGER            NB
      COMMON             /BSPRQR/ NB
*     ..
*     .. Local Scalars ..
      REAL               RCNR, RCNRP1, RCOND
      LOGICAL            GOLEFT, RNKDTD
      INTEGER            MN, OINFO
*     ..
*     .. External Functions ..
      INTEGER            ILAENV
      REAL               SLAMCH
      EXTERNAL           ILAENV, SLAMCH
*     ..
*     .. External Subroutines ..
      EXTERNAL           XERBLA, STRQYC, STRRNK
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MIN, MAX
*     ..
*     .. Executable Statements ..
*
      MN = MIN( M, N )
      NB = ILAENV( INB, 'SGEQRF', ' ', M, N, 0, 0 )
*
*     Test input arguments
*     ====================
*
      INFO = 0
      IF( ( JOB.LT.1 ).OR.( JOB.GT.3 ) ) THEN
         INFO = -1
      ELSE IF( M.LT.0 ) THEN
         INFO = -2
      ELSE IF( N.LT.0 ) THEN
         INFO = -3
      ELSE IF( K.LT.0 ) THEN
         INFO = -4
      ELSE IF( LDA.LT.MAX(1,M) ) THEN
         INFO = -6
      ELSE IF( ( ( JOB.EQ.1 ).AND.( LDC.LT.1 ) ).OR.
     $         ( ( JOB.EQ.2 ).AND.( LDC.LT.MAX( 1, M ) ) ).OR.
     $         ( ( JOB.EQ.3 ).AND.( LDC.LT.MAX( 1, K ) ) ) ) THEN
         INFO = -8
      ELSE IF( IRCOND.LT.ZERO ) THEN
         INFO = -10
      ELSE IF( LWORK.LT.MAX(1,N+3*MN) ) THEN
         INFO = -15
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'STRQPY', -INFO )
         RETURN
      END IF
*
*     Quick return if possible.
*
      IF( MN.EQ.0 ) THEN
         RANK = 0
         ORCOND = ZERO
         SVLUES( IMAX )   = ZERO
         SVLUES( IBEFOR ) = ZERO
         SVLUES( IAFTER ) = ZERO
         SVLUES( IMIN )   = ZERO
         RETURN
      END IF
*
*     Check whether Threshold for condition number was supplied.
*     If not, choose machine precision as default for RCOND.
*
      IF( IRCOND.GT.ZERO ) THEN
         RCOND = IRCOND
      ELSE
         RCOND = SLAMCH( 'Epsilon' )
      END IF
*
*     Compute the initial estimate for the rank.
*
      CALL STRRNK( MN, A, LDA, RCOND, RANK, WORK, INFO )
*
*     ************************
*     * First call to xTRQYC *
*     ************************
*
*     Get tighter bounds for the value RANK.
*
      CALL STRQYC( JOB, M, N, K, A, LDA, C, LDC, JPVT,
     $            RANK, SVLUES, RCNR, RCNRP1, WORK, INFO )
      OINFO = 0
      IF( INFO.NE.0 )
     $   OINFO = INFO
*
*     Check if the numerical rank is larger, equal or smaller than
*     the contents of RANK.
*
      IF( ( ( RCNR.GE.RCOND ).AND.( RANK.EQ.MN ) ).OR.
     $     ( ( RCNR.GE.RCOND ).AND.( RCNRP1.LT.RCOND ) ) ) THEN
         RNKDTD = .TRUE.
      ELSE IF( ( RCNR.GE.RCOND ).AND.( RCNRP1.GE.RCOND ) ) THEN
         RNKDTD = .FALSE.
         GOLEFT = .FALSE.
         RANK = RANK + 1
      ELSE IF( ( RCNR.LT.RCOND ).AND.( RCNRP1.LT.RCOND ) ) THEN
         IF( RANK.EQ.1 ) THEN
            RNKDTD = .TRUE.
            IF( ABS( A( 1, 1 ) ).EQ.ZERO ) THEN
               RANK = 0
               ORCOND = ZERO
               SVLUES( IMAX )   = ZERO
               SVLUES( IBEFOR ) = ZERO
               SVLUES( IAFTER ) = ZERO
               SVLUES( IMIN )   = ZERO
            ELSE
               RANK = 1
            END IF
         ELSE
            RNKDTD = .FALSE.
            GOLEFT = .TRUE.
            RANK = RANK - 1
         END IF
      ELSE
         RNKDTD = .TRUE.
         INFO = 2
      END IF
*
*     *****************
*     * Start of Loop *
*     *****************
*
*     Loop for the detection of the actual rank. The variable RANK is
*     updated until the rank is found. To avoid infinite loops, the
*     variable RANK either increases or decreases.
*
 10   CONTINUE
      IF( .NOT. RNKDTD ) THEN
*
*        Call to xTRQYC to get tighter bounds for the value RANK.
*
         CALL STRQYC( JOB, M, N, K, A, LDA, C, LDC, JPVT,
     $               RANK, SVLUES, RCNR, RCNRP1, WORK, INFO )
         IF( INFO.NE.0 )
     $      OINFO = INFO
*
*        Check if the numerical rank is larger, equal or smaller than
*        the contents of RANK.
*
         IF( ( ( RCNR.GE.RCOND ).AND.( RANK.EQ.MN ) ).OR.
     $        ( ( RCNR.GE.RCOND ).AND.( RCNRP1.LT.RCOND ) ) ) THEN
            RNKDTD = .TRUE.
         ELSE IF( ( RCNR.GE.RCOND ).AND.( RCNRP1.GE.RCOND ) ) THEN
            IF( .NOT. GOLEFT ) THEN
               RANK = RANK + 1
            ELSE
               RNKDTD = .TRUE.
               INFO = 2
            END IF
         ELSE IF( ( RCNR.LT.RCOND ).AND.( RCNRP1.LT.RCOND ) ) THEN
            IF( RANK.EQ.1 ) THEN
               RNKDTD = .TRUE.
               IF( ABS( A( 1, 1 ) ).EQ.ZERO ) THEN
                  RANK = 0
                  ORCOND = ZERO
                  SVLUES( IMAX )   = ZERO
                  SVLUES( IBEFOR ) = ZERO
                  SVLUES( IAFTER ) = ZERO
                  SVLUES( IMIN )   = ZERO
               ELSE
                  RANK = 1
               END IF
            ELSE
               GOLEFT = .TRUE.
               RANK = RANK - 1
            END IF
         ELSE
            RNKDTD = .TRUE.
            INFO = 2
         END IF
*
*        Jump to the beginning of the loop.
*
         GOTO 10
      END IF
*
*     ***************
*     * end of loop *
*     ***************
*
*     Give back the obtained value of RCOND and check the value of INFO.
*
      ORCOND = RCNR
      IF( OINFO.NE.0 )
     $   INFO = OINFO
*
      RETURN
*
*     End of STRQPY
*
      END
