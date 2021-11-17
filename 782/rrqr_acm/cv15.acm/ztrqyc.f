      SUBROUTINE ZTRQYC( JOB, M, N, K, A, LDA, C, LDC, JPVT, RANK,
     $                   SVLUES, RCNR, RCNRP1, WORK, RWORK, INFO )
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
*     $Date: 96/12/30 16:59:44 $
*
*     .. Scalars Arguments ..
      INTEGER            JOB, M, N, K, LDA, LDC, RANK, INFO
      DOUBLE PRECISION   RCNR, RCNRP1
*     ..
*     .. Array Arguments ..
      COMPLEX*16         A( LDA, * ), C( LDC, * ), WORK( * )
      DOUBLE PRECISION   SVLUES( 4 ), RWORK( * )
      INTEGER            JPVT( * )
*     ..
*
*  Purpose
*  =======
*
*  PREC_TRQYC carries out Pan-Tang Algorithm 3 for the stage RANK.
*  This is a mofified version of the original algorithm. The improved
*  features are the following:
*  o Use of Bischof's ICE to reduce the computational cost.
*  o Reorganization of the main loop to save computations.
*  o No permutation is carried out if not strictly needed.
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
*          The number of rows of matrices A and C. M >= 0.
*
*  N       (input) INTEGER
*          The number of columns of matrix A. N >= 0.
*
*  K       (input) INTEGER
*          The number of columns of matrix C. K >= 0.
*
*  A       (input/output) COMPLEX*16 array (LDA,N)
*          Upper triangular m by n matrix.
*
*  LDA     (input) INTEGER
*          The leading dimension of array A. LDA >= MAX( 1, M ).
*
*  C       (input/output) COMPLEX*16 array (LDC,K)
*          Matrix of dimension m x k where to accumulate
*          orthogonal transformations from the left.
*
*  LDC     (input) INTEGER
*          The leading dimension of array C. LDC >= MAX( 1, M ).
*
*  JPVT    (input/output) INTEGER array (N)
*          Vector with the actual permutation of matrix A.
*
*  RANK    (input) INTEGER
*          The estimate for the rank. 1 <= RANK <= MIN(M,N).
*
*  SVLUES  (output) DOUBLE PRECISION array, dimension (4)
*          On exit, SVLUES contains estimates of some of the singular
*          values of the triangular factor R.
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
*  RCNR    (output) DOUBLE PRECISION
*          The estimate for the inverse of the condition number of
*          block R(1:RANK,1:RANK).
*
*  RCNRP1  (output) DOUBLE PRECISION
*          The estimate for the inverse of the condition number of
*          block R(1:RANK+1,1:RANK+1).
*
*  WORK    (workspace) COMPLEX*16 array, dimension (2*MIN(M,N))
*
*  RWORK   (workspace) DOUBLE PRECISION array, dimension (N+MIN(M,N))
*
*  INFO    (output) INTEGER
*          = 0: Successful exit.
*          < 0: If info = -i, the i-th argument had an illegal value.
*          = 4: Exceeded the allowed maximum number of steps. That is,
*               the matrix presents a slow convergence.
*
*
*  Further Details
*  ===============
*
*    If the leading block of R is singular or near singular, there will
*  be no permutation because in that case the right (and left) singular
*  vectors are the canonical ones ((0,0,...0,1)^T).
*    That is, there will not be permutation if
*  RCOND <= SF * DLAMCH('Safe Minimum'), where SF (Safe Factor) is
*  a security factor to avoid arithmetic exceptions.
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   FP, SF
      COMPLEX*16         CONE
      PARAMETER          ( FP = 0.9D+0, SF = 1.0D+2,
     $                   CONE = ( 1.0D+0, 0.0D+0 ) )
*
*     Indices into the 'svlues' array.
*
      INTEGER            IMAX, IBEFOR, IAFTER, IMIN
      PARAMETER          ( IMAX = 1, IBEFOR = 2, IAFTER = 3, IMIN = 4 )
*     ..
*     .. Local Scalars ..
      INTEGER            I, II, ITEMP, J, JJ, MN, MXSTPS, NCA, NCTBA
      DOUBLE PRECISION   F, RCOS, SMAX, SMAXPR, SMIN, SMINPR, SMXRP1,
     $                   SMNRP1, TEMP
      COMPLEX*16         COSINE, CTEMP, DIAG, SINE
      INTEGER            NS
*     .. Local Arrays ..
      COMPLEX*16         CDUMMY( 1, 1 )
*     ..
*     .. External Subroutines ..
      EXTERNAL           ZCOPY, ZGRET, ZHESS, ZLAIC1,
     $                   ZLARTG, ZSCAL, ZSWAP, ZTRSV
*     ..
*     .. External Functions ..
      EXTERNAL           IDAMAX, IZAMAX, DLAMCH,
     $                   DLASMX, DZNRM2
      INTEGER            IDAMAX, IZAMAX
      DOUBLE PRECISION   DLAMCH, DLASMX, DZNRM2
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, MIN, SQRT, DBLE
*     ..
*     .. Executable Statements ..
      MN = MIN( M, N )
      MXSTPS = N+25
      NS = 0
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
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -6
      ELSE IF( ( ( JOB.EQ.1 ).AND.( LDC.LT.1 ) ).OR.
     $         ( ( JOB.EQ.2 ).AND.( LDC.LT.MAX( 1, M ) ) ).OR.
     $         ( ( JOB.EQ.3 ).AND.( LDC.LT.MAX( 1, K ) ) ) ) THEN
         INFO = -8
      ELSE IF( ( RANK.LT.1 ).OR.( RANK.GT.MN ) ) THEN
         INFO = -10
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZTRQYC', -INFO )
         RETURN
      END IF
*
*     Quick return if possible.
*
      IF( MN.EQ.0 )
     $   RETURN
*
      IF( RANK.EQ.MN ) THEN
*
*        ************************
*        ************************
*        * Apply Chan Algorithm *
*        ************************
*        ************************
*
         F = FP
*
*        Move the best column of A(1:M,M:N) to position M-th.
*
         JJ = MN - 1 + IZAMAX( N-MN+1, A( MN, MN ), LDA )
         IF( JJ.GT.MN ) THEN
            CALL ZSWAP( M, A( 1, MN ), 1, A( 1, JJ ), 1 )
            ITEMP = JPVT( MN )
            JPVT( MN ) = JPVT( JJ )
            JPVT( JJ ) = ITEMP
         END IF
*
*        Estimation of the largest singular value, the smallest
*        singular value, and its corresponding left singular vector.
*
         SMAX = ABS( A( 1, 1 ) )
         WORK( 1 ) = CONE
         SMIN = SMAX
         WORK( MN+1 ) = CONE
         DO 10 J = 2, RANK
            CALL ZLAIC1( 1, J-1, WORK( 1 ), SMAX, A( 1, J ),
     $                   A( J, J ), SMAXPR, SINE, COSINE )
            CALL ZSCAL( J-1, SINE, WORK( 1 ), 1 )
            WORK( J ) = COSINE
            SMAX = SMAXPR
*
            CALL ZLAIC1( 2, J-1, WORK( MN+1 ), SMIN, A( 1, J ),
     $                   A( J, J ), SMINPR, SINE, COSINE )
            CALL ZSCAL( J-1, SINE, WORK( MN+1 ), 1 )
            WORK( MN+J ) = COSINE
            SMIN = SMINPR
10       CONTINUE
*
*        Determine if matrix A is singular or nearly singular.
*        SF (Safe Factor) is used to say whether or not it is.
*
         IF( SMIN.GT.( SMAX*SF*DLAMCH( 'Safe Minimum' ) ) ) THEN
*
*           Matrix is not singular or not nearly singular.
*           Follow usual method: Estimate the right singular vector
*           corresponding to the smallest singular value of upper
*           triangular block A(1:RANK,1:RANK).
*
            CALL ZTRSV( 'Upper', 'No transpose', 'No unit',
     $                  RANK, A, LDA, WORK( MN+1 ), 1 )
*
*           Find the index with largest absolute value in vector
*           WORK( MN+1:2*MN ).
*
            JJ = IZAMAX( RANK, WORK( MN+1 ), 1 )
*
*           Permut if necessary.
*
            IF( ( JJ.LT.RANK ).AND.( ( ABS( WORK( MN+JJ ) )*F )
     $         .GT.ABS( WORK( MN+RANK ) ) ) ) THEN
*
               NS = 1
*
*              Exchange cyclically to the left the columns of A between
*              JJ and RANK, that is: JJ->RANK, JJ+1->JJ, JJ+2->JJ+1,...,
*              RANK->RANK-1.
*
               CALL ZCOPY( RANK, A( 1, JJ ), 1, WORK( 1 ), 1 )
               DO 20 J = JJ+1, RANK
                  CALL ZCOPY( J, A( 1, J ), 1, A( 1, J-1 ), 1 )
20             CONTINUE
               CALL ZCOPY( RANK, WORK( 1 ), 1, A( 1, RANK ), 1 )
*
*              Exchange in the same way vector JPVT.
*
               ITEMP = JPVT( JJ )
               DO 30 J = JJ+1, RANK
                  JPVT( J-1 ) = JPVT( J )
30             CONTINUE
               JPVT( RANK ) = ITEMP
*
*              Retriangularization of matrix A after the permutation.
*
               IF( JOB.EQ.1 ) THEN
                  CALL ZHESS( JOB, RANK-JJ+1, N-JJ+1, K,
     $                        A( JJ, JJ ), LDA, CDUMMY, 1,
     $                        WORK( 1 ), RWORK( N+1 ), INFO )
               ELSE IF( JOB.EQ.2 ) THEN
                  CALL ZHESS( JOB, RANK-JJ+1, N-JJ+1, K,
     $                        A( JJ, JJ ), LDA, C( JJ, 1 ), LDC,
     $                        WORK( 1 ), RWORK( N+1 ), INFO )
               ELSE IF( JOB.EQ.3 ) THEN
                  CALL ZHESS( JOB, RANK-JJ+1, N-JJ+1, K,
     $                        A( JJ, JJ ), LDA, C( 1, JJ ), LDC,
     $                        WORK( 1 ), RWORK( N+1 ), INFO )
               END IF
            END IF
         END IF
*
*        Computation of the contents of vector SVLUES, RCNR and RCNRP1.
*
         SVLUES( IMAX ) = SMAX
         SVLUES( IBEFOR ) = SMIN
         SVLUES( IAFTER ) = SMIN
         SVLUES( IMIN ) = SMIN
         RCNR = SVLUES( IBEFOR ) / SVLUES( IMAX )
         RCNRP1 = RCNR
      ELSE
*
*        ***************************************
*        ***************************************
*        * Apply Modified Pan&Tang Algorithm 3 *
*        ***************************************
*        ***************************************
*
*        Adjust the value of f.
*
         F = FP / SQRT( DBLE( RANK+1 ) )
*
*        Compute the norms of columns of matrix A. Store them into
*        vector RWORK(1:N).
*
         DO 100 J = 1, N
            RWORK( J ) = DZNRM2( MIN( M, J ), A( 1, J ), 1 )
100      CONTINUE
*
*        Estimate the smallest singular value of A(1:RANK,1:RANK) and
*        its corresponding left singular vector.
*        SMIN will contain the smallest singular value and
*        WORK(1:MN) will contain the left singular vector.
*
         SMIN = ABS( A( 1, 1 ) )
         WORK( 1 ) = CONE
         DO 110 J = 2, RANK
            CALL ZLAIC1( 2, J-1, WORK( 1 ), SMIN, A( 1, J ),
     $                   A( J , J ), SMINPR, SINE, COSINE )
            CALL ZSCAL( J-1, SINE, WORK( 1 ), 1 )
            WORK( J ) = COSINE
            SMIN = SMINPR
110      CONTINUE
*
*        Initialize loop variables.
*
         NCA = 0
         NCTBA = N-RANK
         II = RANK+1
*
*        ***********************
*        * Start of Loop WHILE *
*        ***********************
*
 1000    IF( ( NCA.LT.NCTBA ).AND.( NS.LT.MXSTPS ) ) THEN
*
*           Estimate the smallest singular value of A(1:RANK+1,1:RANK+1)
*           and its corresponding left singular vector as if column II
*           of matrix A were on column RANK+1.
*
            DIAG = A( MIN( MN, II ), II )
            DO 120 I = MIN( MN, II )-1, RANK+1, -1
               CALL ZLARTG( A( I, II ), DIAG, RCOS, SINE, CTEMP )
               DIAG = CTEMP
120         CONTINUE
*
            CALL ZLAIC1( 2, RANK, WORK( 1 ), SMIN, A( 1, II ),
     $                   DIAG, SMNRP1, SINE, COSINE )
            IF( SMNRP1.GE.( F*ABS( DIAG ) ) ) THEN
*
*              Column II accepted on the right part of matrix A.
*
               NCA = NCA+1
               IF( II.EQ.N ) THEN
                  II = RANK+1
               ELSE
                  II = II+1
               END IF
            ELSE
*
*              Column II not accepted on the right part of matrix A.
*
*              *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
*              * Permut column II to position RANK+1 *
*              *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
*
*              Exchange cyclically to the right the columns of A between
*              RANK+1 and II, that is, RANK+1->RANK+2, RANK+2->RANK+3,
*              ...,II-1->II,II->RANK+1.
*
               CALL ZCOPY( MIN( MN, II ), A( 1, II ), 1,
     $                      WORK( MN+1 ), 1 )
               DO 130 J = II-1, RANK+1, -1
                  CALL ZCOPY( MIN( MN, J+1 ), A( 1, J ), 1 ,
     $                         A( 1, J+1 ), 1 )
130            CONTINUE
               CALL ZCOPY( MIN( MN, II ), WORK( MN+1 ), 1,
     $                      A( 1, RANK+1 ), 1 )
*
*              Exchange in the same way vector JPVT.
*
               ITEMP = JPVT( II )
               DO 140 J = II-1, RANK+1, -1
                  JPVT( J+1 ) = JPVT( J )
140            CONTINUE
               JPVT( RANK+1 ) = ITEMP
*
*              Exchange in the same way vector RWORK(1:N).
*
               TEMP = RWORK( II )
               DO 150 J = II-1, RANK+1, -1
                  RWORK( J+1 ) = RWORK( J )
150            CONTINUE
               RWORK( RANK+1 ) = TEMP
*
*              Retriangularize matrix A after permutation.
*
               IF( JOB.EQ.1 ) THEN
                  CALL ZGRET( JOB, MIN( M, II )-RANK, N-RANK, K,
     $                        A( RANK+1, RANK+1 ), LDA, CDUMMY, 1,
     $                        WORK( MN+1 ), RWORK( N+1 ), INFO )
               ELSE IF( JOB.EQ.2 ) THEN
                  CALL ZGRET( JOB, MIN( M, II )-RANK, N-RANK, K,
     $                        A( RANK+1, RANK+1 ), LDA, C( RANK+1, 1 ),
     $                        LDC, WORK( MN+1 ), RWORK( N+1 ), INFO )
               ELSE IF( JOB.EQ.3 ) THEN
                  CALL ZGRET( JOB, MIN( M, II )-RANK, N-RANK, K,
     $                        A( RANK+1, RANK+1 ), LDA, C( 1, RANK+1 ),
     $                        LDC, WORK( MN+1 ), RWORK( N+1 ), INFO )
               END IF
*
*              Estimate the largest singular value.
*
               ITEMP = IDAMAX( RANK+1, RWORK, 1 )
               SMXRP1 = DLASMX( RANK+1 )*RWORK( ITEMP )
*
*              *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
*              * Estimate the right singular vector  *
*              *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
*
               IF( SMNRP1.GT.
     $            ( SMXRP1*SF*DLAMCH( 'Safe minimum' ) ) ) THEN
*
*                 Matrix is not singular or not nearly singular.
*
*                 First, end the estimation of the left singular vector.
*                 No problem to access WORK(MN+RANK+1) since RANK<MN.
*
                  CALL ZCOPY( RANK, WORK( 1 ), 1, WORK( MN+1 ), 1 )
                  CALL ZSCAL( RANK, SINE, WORK( MN+1 ), 1 )
                  WORK( MN+RANK+1 ) = COSINE
*
*                 Obtain the right singular vector from the left one.
*
                  CALL ZTRSV( 'Upper', 'No transpose', 'No unit',
     $                        RANK+1, A, LDA, WORK( MN+1 ), 1 )
*
                  JJ = IZAMAX( RANK+1, WORK( MN+1 ), 1 )
*
*                 *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
*                 * Permut column JJ to position RANK+1 *
*                 *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
*
                  IF( JJ.LT.( RANK+1 ) ) THEN
*
*                    Exchange cyclically to the left the columns of A
*                    between JJ and RANK+1, that is, JJ->RANK+1,JJ+1->JJ,
*                    JJ+2->JJ+1,...,RANK+1->RANK.
*
                     CALL ZCOPY( RANK+1, A( 1, JJ ), 1,
     $                            WORK( MN+1 ), 1 )
                     DO 160 J = JJ+1, RANK+1
                        CALL ZCOPY( J, A( 1, J ), 1,
     $                               A( 1, J-1 ), 1 )
160                  CONTINUE
                     CALL ZCOPY( RANK+1, WORK( MN+1 ), 1,
     $                            A( 1, RANK+1 ), 1 )
*
*                    Exchange in the same way vector JPVT.
*
                     ITEMP = JPVT( JJ )
                     DO 170 J = JJ+1, RANK+1
                        JPVT( J-1 ) = JPVT( J )
170                  CONTINUE
                     JPVT( RANK+1 ) = ITEMP
*
*                    Exchange in the same way vector RWORK.
*
                     TEMP = RWORK( JJ )
                     DO 180 J = JJ+1, RANK+1
                        RWORK( J-1 ) = RWORK( J )
180                  CONTINUE
                     RWORK( RANK+1 ) = TEMP
*
*                    Retriangularize matrix A after the permutation.
*
                     IF( JOB.EQ.1 ) THEN
                        CALL ZHESS( JOB, RANK-JJ+2, N-JJ+1, K,
     $                              A( JJ, JJ ), LDA, CDUMMY, 1,
     $                              WORK( MN+1 ), RWORK( N+1 ), INFO )
                     ELSE IF( JOB.EQ.2 ) THEN
                        CALL ZHESS( JOB, RANK-JJ+2, N-JJ+1, K,
     $                              A( JJ, JJ ), LDA, C( JJ, 1 ), LDC,
     $                              WORK( MN+1 ), RWORK( N+1 ), INFO )
                     ELSE IF( JOB.EQ.3 ) THEN
                        CALL ZHESS( JOB, RANK-JJ+2, N-JJ+1, K,
     $                              A( JJ, JJ ), LDA, C( 1, JJ ), LDC,
     $                              WORK( MN+1 ), RWORK( N+1 ), INFO )
                     END IF
*
*                    Estimate the smallest singular value of
*                    A(1:RANK,1:RANK) and its corresponding left
*                    singular vector.
*                    SMIN will contain the smallest singular value and
*                    WORK(1:MN) will contain the left singular
*                    vector.
*
                     SMIN = ABS( A( 1, 1 ) )
                     WORK( 1 ) = CONE
                     DO 190 J = 2, RANK
                        CALL ZLAIC1( 2, J-1, WORK( 1 ), SMIN,
     $                               A( 1, J ), A( J , J ), SMINPR,
     $                               SINE, COSINE )
                        CALL ZSCAL( J-1, SINE, WORK( 1 ), 1 )
                        WORK( J ) = COSINE
                        SMIN = SMINPR
190                  CONTINUE
                  END IF
               END IF
*
*              Update loop variables.
*
               NCA = 0
               NS = NS+1
               IF( II.EQ.N ) THEN
                  II = RANK+1
               ELSE
                  II = II+1
               END IF
            END IF
            GOTO 1000
         END IF
*
*        *********************
*        * End of Loop WHILE *
*        *********************
*
*        ******************
*        * Final Pivoting *
*        ******************
*
*        Exchange column in R(RANK+1:M,RANK+1:N) with largest norm to
*        position RANK+1.
*
         JJ = RANK+IDAMAX( N-RANK, RWORK( RANK+1 ), 1 )
         IF( ( JJ.GT.( RANK+1 ) ).AND.
     $       ( F*ABS( RWORK( JJ ) ).GT.ABS( RWORK( RANK+1 ) ) ) ) THEN
*
*           Exchange column JJ to position RANK+1.
*
            CALL ZCOPY( MIN( MN, JJ ), A( 1, JJ ), 1,
     $                   WORK( MN+1 ), 1 )
            DO 200 J = JJ-1, RANK+1, -1
               CALL ZCOPY( MIN( MN, J+1 ), A( 1, J ), 1 ,
     $                      A( 1, J+1 ), 1 )
200         CONTINUE
            CALL ZCOPY( MIN( MN, JJ ), WORK( MN+1 ), 1,
     $                   A( 1, RANK+1 ), 1 )
*
*           Exchange in the same way vector JPVT.
*
            ITEMP = JPVT( JJ )
            DO 210 J = JJ-1, RANK+1, -1
               JPVT( J+1 ) = JPVT( J )
210         CONTINUE
            JPVT( RANK+1 ) = ITEMP
*
*           Retriangularize matrix A after permutation.
*
            IF( JOB.EQ.1 ) THEN
               CALL ZGRET( JOB, MIN( M, JJ )-RANK, N-RANK, K,
     $                     A( RANK+1, RANK+1 ), LDA, CDUMMY, 1,
     $                     WORK( MN+1 ), RWORK( N+1 ), INFO )
            ELSE IF( JOB.EQ.2 ) THEN
               CALL ZGRET( JOB, MIN( M, JJ )-RANK, N-RANK, K,
     $                     A( RANK+1, RANK+1 ), LDA, C( RANK+1, 1 ),
     $                     LDC, WORK( MN+1 ), RWORK( N+1 ), INFO )
            ELSE IF( JOB.EQ.3 ) THEN
               CALL ZGRET( JOB, MIN( M, JJ )-RANK, N-RANK, K,
     $                     A( RANK+1, RANK+1 ), LDA, C( 1, RANK+1 ),
     $                     LDC, WORK( MN+1 ), RWORK( N+1 ), INFO )
            END IF
         END IF
*
*        **************************************************************
*        * Computation of vector SVLUES and variables RCNR and RCNRP1 *
*        **************************************************************
*
*        Computation of the largest singular value of A(1:RANK,1:RANK).
*
         SMAX = ABS( A( 1, 1 ) )
         WORK( MN+1 ) = CONE
*
         DO 220 J = 2, RANK
            CALL ZLAIC1( 1, J-1, WORK( MN+1 ), SMAX,
     $                   A( 1, J ), A( J, J ), SMAXPR, SINE, COSINE )
            CALL ZSCAL( J-1, SINE, WORK( MN+1 ), 1 )
            WORK( MN+J ) = COSINE
            SMAX = SMAXPR
220      CONTINUE
         SVLUES( IMAX ) = SMAX
         SVLUES( IBEFOR ) = SMIN
*
*        Computation of the largest singular value and the smallest
*        singular value of A(1:RANK+1,1:RANK+1).
*
         CALL ZLAIC1( 1, RANK, WORK( MN+1 ), SMAX,
     $                A( 1, RANK+1 ), A( RANK+1, RANK+1 ), SMXRP1,
     $                SINE, COSINE )
         CALL ZLAIC1( 2, RANK, WORK( 1 ), SMIN,
     $                A( 1, RANK+1 ), A( RANK+1, RANK+1 ), SMINPR,
     $                SINE, COSINE )
         CALL ZSCAL( RANK, SINE, WORK( 1 ), 1 )
         WORK( RANK+1 ) = COSINE
         SMIN = SMINPR
         SVLUES( IAFTER ) = SMIN
*
*        Computation of the smallest singular value of A(1:MN,1:MN).
*
         DO 230 J = RANK+2, MN
            CALL ZLAIC1( 2, J-1, WORK( 1 ), SMIN,
     $                   A( 1, J ), A( J, J ), SMINPR, SINE, COSINE )
            CALL ZSCAL( J-1, SINE, WORK( 1 ), 1 )
            WORK( J ) = COSINE
            SMIN = SMINPR
230      CONTINUE
         SVLUES( IMIN ) = SMIN
*
*        Computation of RCNR and RCNRP1.
*
         RCNR = SVLUES( IBEFOR ) / SVLUES( IMAX )
         RCNRP1 = SVLUES( IAFTER ) / SMXRP1
*
         IF( NS.GE.MXSTPS ) THEN
            INFO = 1
         END IF
      END IF
      RETURN
*
*     End of ZTRQYC
*
      END
