      SUBROUTINE CGEQPB( JOB, M, N, K, A, LDA, C, LDC, JPVT,
     $                   IRCOND, ORCOND, RANK, SVLUES, WORK, LWORK,
     $                   RWORK, INFO )
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
*     $Revision: 1.42 $
*     $Date: 96/12/30 16:59:36 $
*
*     .. Scalar Arguments ..
      INTEGER            JOB, M, N, K, LDA, LDC, RANK, LWORK, INFO
      REAL               IRCOND, ORCOND
*     ..
*     .. Array Arguments ..
      INTEGER            JPVT( * )
      COMPLEX            A( LDA, * ), C( LDC, * ), WORK( * )
      REAL               SVLUES( 4 ), RWORK( * )
*     ..
*
*  Purpose
*  =======
*
*  CGEQPB computes a QR factorization
*       A*P = Q*[ R11 R12 ]
*               [  0  R22 ]
*  of a real m by n matrix A. The permutation P is
*  chosen with the goal to reveal the rank of A by a
*  suitably dimensioned trailing submatrix R22 with norm(R22)
*  being small.
*
*  Arguments
*  =========
*
*  JOB     (input) INTEGER
*          The job to do:
*          = 1: The transformations needed in the
*               triangularization are only applied to matrix A.
*               Thus, matrix C is not updated.
*          = 2: The same transformations needed in the
*               triangularization of matrix A are applied to
*               matrix C from the left.
*               That is, if Q'*A*P=R, then C := Q'*C.
*               In this case, matrix C is m-by-k.
*          = 3: The transpose of the transformations needed
*               in the triangularization of matrix A are applied
*               to matrix C from the right.
*               That is, if Q'*A*P=R, then C := C*Q.
*               In this case, matrix C is k-by-m.
*          In these three cases, the permutations are always saved
*          into vector JPVT.
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
*  A       (input/output) COMPLEX array, dimension (LDA,N)
*          On entry, the m by n matrix A.
*          On exit, the upper triangle of the array contains the
*          min(m,n) by n upper trapezoidal matrix R; the lower triangle
*          array is filled with zeros.
*
*  LDA     (input) INTEGER
*          The leading dimension of array A. LDA >= MAX(1,M).
*
*  C       (input/output) COMPLEX array, dimension
*                ( LDC, K ) if JOB=2.
*                ( LDC, M ) if JOB=3.
*          If argument JOB asks, all the transformations
*          applied to matrix A are also applied to matrix C.
*
*  LDC     (input) INTEGER
*          The leading dimension of array C.
*          If JOB=2, then LDC >= MAX(1,M).
*          If JOB=3, then LDC >= MAX(1,K).
*
*  JPVT    (output) INTEGER array, dimension (N)
*          JPVT(I) = K <==> Column K of A has been permuted into
*                           position I in AP.
*          JPVT(1:RANK) contains the indices of the columns considered
*          linearly independent.
*          JPVT(RANK+1:N) contains the indices of the columns considered
*          linearly dependent from the previous ones.
*
*  IRCOND  (input) FLOATING_DECLARE
*          1/IRCOND specifies an upper bound on the condition number
*          of R11. If IRCOND == 0, IRCOND = machine precision is chosen
*          as default. IRCOND must be >= 0.
*
*  ORCOND  (output) FLOATING_DECLARE
*          1/ORCOND is an estimate for the condition number of R11.
*
*  RANK    (output) INTEGER
*          RANK is an estimate for the numerical rank of A with respect
*          to the threshold 1/IRCOND in the sense that
*               RANK = arg_max(cond_no(R(1:r,1:r))<1/IRCOND)
*          This may be an underestimate of the rank if the leading
*          columns were not well-conditioned.
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
*  WORK    (workspace) COMPLEX array, dimension (LWORK)
*          On exit: WORK(1) is the size of the storage array needed
*                   for optimal performance
*
*  LWORK   (input) INTEGER
*          The dimension of array WORK.
*          If JOB=1:
*             The unblocked strategy requires that:
*                 LWORK >= 2*MN+N.
*             The block algorithm requires that:
*                 LWORK >= 2*MN+N*NB.
*          If JOB<>1:
*             The unblocked strategy requires that:
*                 LWORK >= 2*MN+MAX(K,N).
*             The block algorithm requires that:
*                 LWORK >= 2*MN+NB*NB+NB*MAX(K,N).
*          Where MN = min(M,N) and NB is the maximum of blocksize
*          used within xGEQRF and blocksize used within xUNMQR.
*          In both cases, the minimum required workspace is the
*          one for the unblocked strategy.
*
*  RWORK   (workspace) REAL array, dimension (2*N)
*
*  INFO    (output) INTEGER
*          = 0: Successful exit
*          < 0: If INFO = -i, the i-th argument had an illegal value
*
*  =====================================================================
*
*     .. Parameters ..
      REAL               ZERO
      COMPLEX            CZERO
      PARAMETER          ( ZERO = 0.0E+0, CZERO = ( 0.0E+0, 0.0E+0 ) )
      INTEGER            INB, INBMIN, IXOVER
      PARAMETER          ( INB = 1, INBMIN = 2, IXOVER = 3 )
*
*     Indices into the 'svlues' array.
*
      INTEGER            IMAX, IBEFOR, IAFTER, IMIN
      PARAMETER          ( IMAX = 1, IBEFOR = 2, IAFTER = 3, IMIN = 4 )
*     ..
*     .. Local Scalars ..
      INTEGER            I, J, MN, ITEMP, KK, LACPTD, MVIDX, STREJ,
     $                   ACCPTD, NB, LWSIZE, NLLITY, KB, WSIZE, WKMIN
      REAL               SMIN, MXNM, RCOND
      LOGICAL            BLOCK
*     ..
*     .. External Subroutines ..
      EXTERNAL           XERBLA, CGEQPC, CGEQPW,
     $                   CLARFT, CLARFB
*     ..
*     .. External Functions ..
      INTEGER            ILAENV
      REAL               SLAMCH, SLASMX
      EXTERNAL           ILAENV, SLAMCH, SLASMX
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN, REAL, INT
*     ..
*     .. Executable Statements ..
*
      MN = MIN( M, N )
*
*     Compute the minimum required complex workspace.
*
      IF( JOB.EQ.1 ) THEN
         WKMIN = 2*MN + N
      ELSE
         WKMIN = 2*MN + MAX( N, K )
      END IF
*
*     Test input arguments
*     ====================
*
      INFO = 0
      IF( ( JOB.LT.1 ).OR.( JOB.GT.3 ) ) THEN
         INFO = -1
      ELSE IF( M.LT.0 ) THEN
         INFO = -2
      ELSE IF( M.LT.0 ) THEN
         INFO = -3
      ELSE IF( K.LT.0 ) THEN
         INFO = -4
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -6
      ELSE IF( ( ( JOB.EQ.1 ).AND.( LDC.LT.1 ) ).OR.
     $         ( ( JOB.EQ.2 ).AND.( LDC.LT.MAX( 1, M ) ) ).OR.
     $         ( ( JOB.EQ.3 ).AND.( LDC.LT.MAX( 1, K ) ) ) ) THEN
         INFO = -8
      ELSE IF( IRCOND.LT.ZERO ) THEN
         INFO = -10
      ELSE IF( LWORK.LT.MAX( 1, WKMIN ) ) THEN
         INFO = -15
      END IF
      IF( ( INFO .EQ. 0 .OR. INFO .EQ. -15 ).AND. LWORK.GE.1 ) THEN
*
*        Compute the optimal complex workspace.
*
         IF( JOB.EQ.1 ) THEN
            NB = ILAENV( INB, 'CGEQRF', ' ', M, N, 0, 0 )
            WSIZE = 2*MN + MAX( 3*N, N*NB )
         ELSE
            NB = MAX( ILAENV( INB, 'CGEQRF', ' ', M, N, 0, 0 ),
     $                ILAENV( INB, 'CUNMQR', ' ', M, N, 0, 0 ) )
            WSIZE = MAX( 2*MN + MAX( N, K ),
     $                   2*MN + NB*NB + NB*MAX( N, K ) )
         END IF
         WORK( 1 ) = REAL( WSIZE )
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'CGEQPB', -INFO )
         RETURN
      END IF
*
*     Initialization of vector JPVT.
*
      DO 70 J = 1, N
         JPVT( J ) = J
 70   CONTINUE
*
*     Quick return if possible
*
      IF( MN.EQ.0 ) THEN
         RANK = 0
         ORCOND = ZERO
         SVLUES( IMAX ) = ZERO
         SVLUES( IBEFOR ) = ZERO
         SVLUES( IAFTER ) = ZERO
         SVLUES( IMIN ) = ZERO
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
*     Determine block size and whether to use blocked code at all
*
      IF( LWORK.LT.WSIZE ) THEN
         IF( JOB.EQ.1 ) THEN
            NB = ( LWORK-2*MN )/N
         ELSE
            ITEMP = INT( SQRT( REAL(
     $               MAX( K, N )**2+4*LWORK-8*MN ) ) )
            NB = ( ITEMP-MAX( K, N ) )/2
         END IF
      END IF
*
      BLOCK = ( ( NB.GT.1 ).AND.
     $   ( NB.GE.ILAENV( INBMIN, 'CGEQRF', ' ', M, N, 0, 0 ) ).AND.
     $   ( MN.GE.ILAENV( IXOVER, 'CGEQRF', ' ', M, N, 0, 0 ) ) )
*
*     The size of the pivot window is chosen to be NB + NLLITY
*     for the blocked algorithm.
*
      NLLITY = MIN( MN, MAX( 10, NB/2+(N*5)/100 ) )
*
*     ***************************************************
*     * Move column with largest residual norm up front *
*     ***************************************************
*
      CALL CGEQPC( JOB, M, N, K, A, LDA, C, LDC, 1, 0,
     $             RCOND, LACPTD, JPVT, WORK( 1 ), WORK( MN+1 ),
     $             SVLUES, MXNM, WORK( 2*MN+1 ), LWORK-2*MN, RWORK )
      IF( LACPTD.EQ.1 ) THEN
         IF( LACPTD.EQ.MN ) THEN
            RANK = 1
            ORCOND = SVLUES( IBEFOR )/SVLUES( IMAX )
            GOTO 30
         ELSE
            SMIN = SVLUES( IBEFOR )
         END IF
      ELSE
         RANK = 0
         ORCOND = ZERO
         SVLUES( IMAX ) = ZERO
         SVLUES( IBEFOR ) = ZERO
         SVLUES( IAFTER ) = ZERO
         SVLUES( IMIN ) = ZERO
         GOTO 30
      END IF
*
*     ****************************
*     * Factor remaining columns *
*     ****************************
*
      IF( BLOCK ) THEN
*
*        *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
*        * Using blocked code with restricted pivoting strategy  *
*        *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
*
         STREJ = N+1
         KK = 2
*
 10      IF( ( KK.GE.STREJ ).OR.( KK.GT.MN ) ) GOTO 20
*
*           invariant: A(:,KK) is the first column in currently
*           considered block column.
*
            KB = MIN( NB, MIN( MN+1, STREJ )-KK )
*
*           The goal now is to find "KB" independent columns
*           among the remaining STREJ-KK not yet rejected columns.
*
            LWSIZE = MIN( STREJ-KK, KB+NLLITY )
            CALL CGEQPW( M, LWSIZE, KB, KK-1, LACPTD, A, LDA, JPVT,
     $                   RCOND, WORK( MN+1 ), SMIN, MXNM,
     $                   WORK( 1 ), WORK( 2*MN+1 ), RWORK )
            IF( LACPTD.GT.0 ) THEN
*
*              Accumulate Householder vectors in a block reflector.
*
               CALL CLARFT( 'Forward', 'Columnwise', M-KK+1,
     $                      LACPTD, A( KK, KK ), LDA, WORK( KK ),
     $                      WORK( 2*MN+1 ), LACPTD )
*
*              Apply block reflector to A(KK:M,KK+LWSIZE:N).
*
               IF( ( KK+LWSIZE ).LE.N ) THEN
                  CALL CLARFB( 'Left', 'Conjugate Transpose',
     $                         'Forward', 'Columnwise',
     $                         M-KK+1, N-KK-LWSIZE+1, LACPTD,
     $                         A( KK, KK ), LDA,
     $                         WORK( 2*MN+1 ), LACPTD,
     $                         A( KK, KK+LWSIZE ), LDA,
     $                         WORK( 2*MN+LACPTD*LACPTD+1 ),
     $                            N-KK-LWSIZE+1 )
               END IF
*
*              Apply block reflector to the corresponding part
*              of matrix C.
*
               IF( ( JOB.EQ.2 ).AND.( K.GT.0 ) ) THEN
*
*                 Apply it to matrix C(KK:M,1:K) from the left.
*
                  CALL CLARFB( 'Left', 'Conjugate Transpose',
     $                         'Forward', 'Columnwise', M-KK+1, K,
     $                         LACPTD, A( KK, KK ), LDA,
     $                         WORK( 2*MN+1 ), LACPTD,
     $                         C( KK, 1 ), LDC,
     $                         WORK( 2*MN+LACPTD*LACPTD+1 ), K )
               ELSE IF( ( JOB.EQ.3 ).AND.( K.GT.0 ) ) THEN
*
*                 Apply the transpose of it to matrix C(1:K,KK:M)
*                 from the right.
*
                  CALL CLARFB( 'Right', 'No Transpose', 'Forward',
     $                         'Columnwise', K, M-KK+1, LACPTD,
     $                         A( KK, KK ), LDA,
     $                         WORK( 2*MN+1 ), LACPTD,
     $                         C( 1, KK ), LDC,
     $                         WORK( 2*MN+LACPTD*LACPTD+1 ), K )
               END IF
            END IF
*
*           Move rejected columns to the end if there is space.
*
            IF( LACPTD.LT.KB ) THEN
               IF( STREJ.LE.( KK+LWSIZE ) ) THEN
                  STREJ = KK + LACPTD
               ELSE
                  MVIDX = STREJ
                  DO 40 I = KK+LACPTD,
     $                      MIN( KK+LWSIZE-1, STREJ-LWSIZE+LACPTD-1 )
                     MVIDX = MVIDX - 1
                     CALL CSWAP( M, A( 1, I ),1, A( 1, MVIDX ),1 )
                     ITEMP = JPVT( I )
                     JPVT( I ) = JPVT( MVIDX )
                     JPVT( MVIDX ) = ITEMP
 40               CONTINUE
                  STREJ = MVIDX
               END IF
            END IF
            KK = KK + LACPTD
         GOTO 10
 20      CONTINUE
         ACCPTD = KK-1
         SVLUES( IMAX ) = SLASMX( ACCPTD )*MXNM
         SVLUES( IBEFOR ) = SMIN
         IF( ACCPTD.LT.MN ) THEN
*
*           Process rejected columns.
*
            CALL CGEQPC( JOB, M, N, K, A, LDA, C, LDC, MN-KK+1,
     $                   KK-1, RCOND, LACPTD, JPVT, WORK( 1 ),
     $                   WORK( MN+1 ), SVLUES, MXNM, WORK( 2*MN+1 ),
     $                   LWORK-2*MN, RWORK )
            RANK = ACCPTD + LACPTD
         ELSE
            RANK = ACCPTD
            SVLUES( IAFTER ) = SMIN
            SVLUES( IMIN ) = SMIN
         END IF
      ELSE
*
*        *-*-*-*-*-*-*-*-*-*-*-*-*
*        * using unblocked code  *
*        *-*-*-*-*-*-*-*-*-*-*-*-*
*
         ACCPTD = 1
         CALL CGEQPC( JOB, M, N, K, A, LDA, C, LDC, MN-ACCPTD,
     $                ACCPTD, RCOND, LACPTD, JPVT, WORK( 1 ),
     $                WORK( MN+1 ), SVLUES, MXNM, WORK( 2*MN+1 ),
     $                LWORK-2*MN, RWORK )
         RANK = ACCPTD+LACPTD
*
      END IF
      ORCOND = SVLUES( IBEFOR )/SVLUES( IMAX )
*
*     Nullify the lower part of matrix A.
*
 30   CONTINUE
      DO 50 J = 1, MN
         DO 60 I = J+1, M
            A( I, J ) = CZERO
 60      CONTINUE
 50   CONTINUE
*
      WORK( 1 ) = REAL( WSIZE )
      RETURN
*
*     End of CGEQPB
*
      END
