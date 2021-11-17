      SUBROUTINE CTRQXC( JOB, M, N, K, A, LDA, C, LDC, JPVT,
     $                  RANK, SVLUES, RCNR, RCNRP1, WORK, RWORK, INFO )
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
*     $Date: 96/12/30 16:59:43 $
*
*     .. Scalar Arguments ..
      INTEGER            JOB, M, N, K, LDA, LDC, RANK, INFO
      REAL               RCNR, RCNRP1
*     ..
*     .. Array Arguments ..
      COMPLEX            A( LDA, * ), C( LDC, * ), WORK( * )
      REAL               SVLUES( 4 ), RWORK( * )
      INTEGER            JPVT( * )
*     ..
*
*  Purpose
*  =======
*
*  CTRQXC carries out an algorithm related to algorithm Hybrid-III
*  by Chandrasekaran and Ipsen for the stage RANK. The algorithm used
*  here offers the following advantages:
*  o It is faster since it is based on Chan-II instead of Stewart-II.
*  o This algorithm uses the F factor technique to reduce the number of
*    cycling problems due to roundoff errors.
*  o The final steps that do not improve the ordering are saved.
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
*  A       (input/output) COMPLEX array, dimension (LDA,N)
*          On entry, the m by n matrix A.
*          On exit, the upper triangle of the array contains the
*          min(m,n) by n upper trapezoidal matrix R; the lower triangle
*          array is filled with zeros.
*
*  LDA     (input) INTEGER
*          The leading dimension of array A. LDA >= max(1,M).
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
*  JPVT    (input/output) INTEGER array, dimension ( N )
*          If JPVT(I) = K, then the Ith column of the permuted
*          A was the Kth column of the original A (just before
*          the preprocessing). If a permutation occurs, JPVT will
*          be updated correctly.
*
*  RANK    (input) INTEGER
*          The estimate of the rank. 1 <= RANK <= MIN(M,N).
*
*  SVLUES  (output) REAL array, dimension (4)
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
*          value of A, and SVLUES(4) wil be an estimate for the
*          smallest singular value of A.
*          By examining these values, one can confirm that the rank is
*          well defined with respect to the threshold chosen.
*
*  RCNR    (output) REAL
*          The estimate for the inverse of the condition number of
*          block R(1:RANK,1:RANK).
*
*  RCNRP1  (output) REAL
*          The estimate for the inverse of the condition number of
*          block R(1:RANK+1,1:RANK+1).
*
*  WORK    (workspace) COMPLEX array, dimension ( 2*MIN(M,N) ).
*
*  RWORK   (workspace) REAL array, dimension ( MIN(M,N)+N ).
*
*  INFO    (output) INTEGER
*          = 0: Successful exit.
*          < 0: If info = -i, the i-th argument had an illegal value.
*          = 4: Exceeded the allowed maximum number of steps. That is,
*               the matrix presents a slow convergence.
*
*
*  ===================================================================
*
*     .. Parameters ..
      REAL               F
      COMPLEX            CONE
      PARAMETER          (  F = 0.5E+0, CONE = ( 1.0E+0, 0.0E+0 ) )
*
*     Indices into the 'svlues' array.
*
      INTEGER            IMAX, IBEFOR, IAFTER, IMIN
      PARAMETER          ( IMAX = 1, IBEFOR = 2, IAFTER = 3, IMIN = 4 )
*     ..
*     .. Local Scalars ..
      REAL               SMAX, SMAXPR, SMIN, SMINPR, SMXRP1
      COMPLEX            COSINE, SINE
      LOGICAL            PERMUT
      INTEGER            J, MN, MXSTPS, NACPTD
      INTEGER            NS
*     ..
*     .. External Subroutines ..
      EXTERNAL           CSCAL
*     ..
*     .. External Functions ..
      INTEGER            ISAMAX
      REAL               SLASMX, SCNRM2
      EXTERNAL           ISAMAX, SLASMX, SCNRM2
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
*     ..
*     .. Executable Statements ..
*
      MN = MIN( M, N )
      NS = 0
      MXSTPS = N + 25
      INFO = 0
*
*     Quick return if possible.
*
      IF( MN.EQ.0 )
     $   RETURN
*
*     Inicialization of variable NACPTD, which controls main loop.
*
      NACPTD = 0
*
*     Compute the norms of block A(1:RANK,1:RANK) and store them
*     in vector RWORK(1:RANK). It is computed only once at the
*     beginning and updated every iteration. It is used to estimate
*     the largest singular value in order to pass it to Chan-II.
*
      DO 10 J = 1, RANK
         RWORK( J ) = SCNRM2( J, A( 1, J ), 1 )
 10   CONTINUE
*
*     *****************
*     * start of loop *
*     *****************
*
 20   CONTINUE
*
*     *-*-*-*-*-*-*-*-*-*-*-*-*
*     * call to Golub-I(rank) *
*     *-*-*-*-*-*-*-*-*-*-*-*-*
*
      IF( NACPTD.LT.4 ) THEN
*
*        Apply Golub-I for the stage RANK.
*
         CALL CGLBIF( JOB, M, N, K, A, LDA, C, LDC, JPVT,
     $                F, RANK, PERMUT, WORK, RWORK( MN+1 ), INFO )
*
*        If necessary, update the contents of WORK(RANK).
*
         IF( PERMUT )
     $      RWORK( RANK ) = SCNRM2( RANK, A( 1, RANK ), 1 )
*
*        Update variables NACPTD and NS.
*
         IF( PERMUT ) THEN
            NACPTD = 1
         ELSE
            NACPTD = NACPTD+1
         END IF
         NS = NS + 1
      END IF
*
*     *-*-*-*-*-*-*-*-*-*-*-*-*-*
*     * call to Golub-I(rank+1) *
*     *-*-*-*-*-*-*-*-*-*-*-*-*-*
*
      IF( NACPTD.LT.4 ) THEN
*
*        Determine if the application of Golub-I(rank+1) is necessary.
*
         IF( RANK.EQ.MN ) THEN
*
*           Not necessary. Therefore, no permutation occurs.
*
            PERMUT = .FALSE.
         ELSE
*
*           Apply Golub-I for the stage RANK+1.
*
            CALL CGLBIF( JOB, M, N, K, A, LDA, C, LDC, JPVT,
     $                   F, RANK+1, PERMUT, WORK, RWORK( MN+1 ), INFO )
*
*           Update variable NS.
*
            NS = NS+1
         END IF
*
*        Update variable NACPTD.
*
         IF( PERMUT ) THEN
            NACPTD = 1
         ELSE
            NACPTD = NACPTD+1
         END IF
      END IF
*
*     *-*-*-*-*-*-*-*-*-*-*-*-*-*
*     * call to Chan-II (rank+1)*
*     *-*-*-*-*-*-*-*-*-*-*-*-*-*
*
      IF( NACPTD.LT.4 ) THEN
*
*        Determine if the application of Chan-II(rank+1) is necessary.
*
         IF( RANK.EQ.MN ) THEN
*
*           Not necessary. Therefore, no permutation occurs.
*
            PERMUT = .FALSE.
         ELSE
*
*           Extend vector WORK(1:RANK) to vector WORK(1:RANK+1).
*           So, pivoting vector WORK(1:N) inside Chan-II will be
*           easier.
*
            RWORK( RANK+1 ) = SCNRM2( RANK+1, A( 1, RANK+1 ), 1 )
*
*           Apply Chan-II for the stage RANK+1
*           on block A(1:RANK+1,1:RANK+1).
*
            CALL CCNIIF( JOB, M, N, K, A, LDA, C, LDC, JPVT,
     $                   RWORK, F, RANK+1, PERMUT, WORK,
     $                   RWORK( MN+1 ), INFO )
*
*           Update variable NS.
*
            NS = NS+1
         END IF
*
*        Update variable NACPTD.
*
         IF( PERMUT ) THEN
            NACPTD = 1
         ELSE
            NACPTD = NACPTD+1
         END IF
      END IF
*
*     *-*-*-*-*-*-*-*-*-*-*-*-*
*     * call to Chan-II(rank) *
*     *-*-*-*-*-*-*-*-*-*-*-*-*
*
      IF( NACPTD.LT.4 ) THEN
*
*        Apply Chan-II for the stage RANK on block A(1:RANK,1:RANK).
*
         CALL CCNIIF( JOB, M, N, K, A, LDA, C, LDC, JPVT,
     $                RWORK, F, RANK, PERMUT, WORK,
     $                RWORK( MN+1 ), INFO )
*
*        Update variables NACPTD and NS.
*
         IF( PERMUT ) THEN
            NACPTD = 1
         ELSE
            NACPTD = NACPTD+1
         END IF
         NS = NS + 1
      END IF
*
*     Check if loop must finish.
*
      IF( NS.GE.MXSTPS ) THEN
         INFO = 1
      ELSE IF( NACPTD.LT.4 ) THEN
         GOTO 20
      END IF
*
*     ***************
*     * end of loop *
*     ***************
*
*     Computation of the largest singular value of A(1:RANK,1:RANK).
*
      SMAX = ABS( A( 1, 1 ) )
      WORK( 1 ) = CONE
      SMIN = SMAX
      WORK( MN+1 ) = CONE
*
      DO 30 J = 2, RANK
         CALL CLAIC1( 1, J-1, WORK( 1 ), SMAX,
     $                A( 1, J ), A( J, J ), SMAXPR, SINE, COSINE )
         CALL CSCAL( J-1, SINE, WORK( 1 ), 1 )
         WORK( J ) = COSINE
         SMAX = SMAXPR
         CALL CLAIC1( 2, J-1, WORK( MN+1 ), SMIN,
     $                A( 1, J ), A( J, J ), SMINPR, SINE, COSINE )
         CALL CSCAL( J-1, SINE, WORK( MN+1 ), 1 )
         WORK( MN+J ) = COSINE
         SMIN = SMINPR
 30   CONTINUE
      SVLUES( IMAX ) = SMAX
      SVLUES( IBEFOR ) = SMIN
*
*     Computation of the largest singular value and the smallest
*     singular value of A(1:RANK+1,1:RANK+1).
*
      IF( RANK.LT.MN ) THEN
         CALL CLAIC1( 1, RANK, WORK( 1 ), SMAX,
     $                A( 1, RANK+1 ), A( RANK+1, RANK+1 ), SMAXPR,
     $                SINE, COSINE )
         SMAX = SMAXPR
         CALL CLAIC1( 2, RANK, WORK( MN+1 ), SMIN,
     $                A( 1, RANK+1 ), A( RANK+1, RANK+1 ), SMINPR,
     $                SINE, COSINE )
         CALL CSCAL( RANK, SINE, WORK( MN+1 ), 1 )
         WORK( MN+RANK+1 ) = COSINE
         SMIN = SMINPR
      END IF
      SMXRP1 = SMAX
      SVLUES( IAFTER ) = SMIN
*
*     Computation of the smallest singular value of A(1:MN,1:MN).
*
      DO 40 J = RANK+2, MN
         CALL CLAIC1( 2, J-1, WORK( MN+1 ), SMIN,
     $                A( 1, J ), A( J, J ), SMINPR, SINE, COSINE )
         CALL CSCAL( J-1, SINE, WORK( MN+1 ), 1 )
         WORK( MN+J ) = COSINE
         SMIN = SMINPR
 40   CONTINUE
      SVLUES( IMIN ) = SMIN
*
*     Computation of RCNR and RCNRP1.
*
      RCNR = SVLUES( IBEFOR ) / SVLUES( IMAX )
      RCNRP1 = SVLUES( IAFTER ) / SMXRP1
      RETURN
*
*     End of CTRQXC
*
      END
      SUBROUTINE CGLBIF( JOB, M, N, K, A, LDA, C, LDC, JPVT,
     $                    F, RANK, PERMUT, WORK, RWORK, INFO )
*
*
*     .. Scalar Arguments ..
      INTEGER            JOB, M, N, K, LDA, LDC, RANK, INFO
      REAL               F
      LOGICAL            PERMUT
*     ..
*     .. Array Arguments ..
      COMPLEX            A( LDA, * ), C( LDC, * ), WORK( * )
      REAL               RWORK( * )
      INTEGER            JPVT( * )
*     ..
*
*  Purpose
*  =======
*
*  CGLBIF computes the column index of A(RANK:M,RANK:N) with largest
*  norm and determines if pivoting is necessary. If so, it pivots it
*  into column RANK, permuts vector JPVT, adjusts vector VNORM and
*  permuts and retriangularizes matrix A. It does only one permutation.
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
*  A       (input/output) COMPLEX array, dimension (LDA,N)
*          On entry, the m by n matrix A.
*          On exit, the upper triangle of the array contains the
*          min(m,n) by n upper trapezoidal matrix R; the lower triangle
*          array is filled with zeros.
*
*  LDA     (input) INTEGER
*          The leading dimension of array A. LDA >= max(1,M).
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
*  JPVT    (input/output) INTEGER array, dimension ( N )
*          If JPVT(I) = K, then the Ith column of the permuted
*          A was the Kth column of the original A (just before the
*          preprocessing). If a permutation occurs, it will be
*          updated correctly.
*
*  F       (input) REAL
*          F factor for the pivoting. It must be always 0 < f <= 1.
*
*  RANK    (input) INTEGER
*          The estimate for the rank. 1 <= RANK <= MIN(M,N).
*
*  PERMUT  (output) LOGICAL
*          Tells if a permutation occurred.
*
*  WORK    (workspace) COMPLEX array, dimension ( MIN(M,N) )
*
*  RWORK   (workspace) REAL array, dimension ( N )
*
*  INFO    (output) INTEGER
*          = 0: Successful exit.
*          < 0: If info = -i, the i-th argument had an illegal value.
*
*  =====================================================================
*
*     .. Parameters ..
*     ..
*     .. Local Scalars ..
      INTEGER             MN, JJ, J, ITEMP
*     ..
*     .. Local Scalars ..
      COMPLEX            CDUMMY( 1, 1 )
*     ..
*     .. External Subroutines ..
      EXTERNAL           CCOPY, CGRET
*     ..
*     .. External Functions ..
      INTEGER            ISAMAX
      REAL               SCNRM2
      EXTERNAL           ISAMAX, SCNRM2
*     ..
*     .. Executable Statements ..
*
      MN = MIN( M, N )
      INFO = 0
*
*     Quick return if possible.
*
      IF( ( MN.EQ.0 ).OR.( RANK.EQ.N ) ) THEN
         PERMUT = .FALSE.
         RETURN
      END IF
*
*     Compute the norms of the columns of block A(RANK:M,RANK:N)
*     and store them in vector RWORK(RANK:N).
*
      DO 10 J = RANK, N
         RWORK( J ) =
     $       SCNRM2( MIN( M, J )-RANK+1, A( RANK, J ), 1 )
 10   CONTINUE
*
*     Find column with largest two-norm of upper triangular block
*     A(RANK:M,RANK:N). Use the data stored in vector RWORK(RANK:N).
*
      JJ = RANK - 1 + ISAMAX( N-RANK+1, RWORK( RANK ), 1)
*
*     Determine if a permutation must occur.
*
      PERMUT = ( ( JJ.GT.RANK ).AND.
     $         ( ( ABS( RWORK( JJ ) )*F ).GT.ABS( RWORK( RANK ) ) ) )
*
      IF( PERMUT ) THEN
*
*        Exchage cyclically to the right the columns of matrix A
*        between RANK and JJ. That is, RANK->RANK+1,
*        RANK+1->RANK+2,...,JJ-1->JJ,JJ->K. Use vector WORK(1:MN)
*        to store temporal data.
*
         CALL CCOPY( MIN( MN, JJ ), A( 1, JJ ), 1, WORK, 1 )
         DO 20 J = JJ-1, RANK, -1
            CALL CCOPY( MIN( MN, J+1 ), A( 1, J ), 1,
     $                  A( 1, J+1 ), 1 )
 20      CONTINUE
         CALL CCOPY( MIN( MN, JJ ), WORK, 1, A( 1, RANK ), 1 )
*
*        Exchange in the same way vector JPVT.
*
         ITEMP = JPVT( JJ )
         DO 30 J = JJ-1, RANK, -1
            JPVT( J+1 ) = JPVT( J )
 30      CONTINUE
         JPVT( RANK ) = ITEMP
*
*        Retriangularize matrix A after the permutation.
*
         IF( JOB.EQ.1 ) THEN
            CALL CGRET( JOB, MIN( M, JJ )-RANK+1, N-RANK+1, K,
     $                  A( RANK, RANK ), LDA, CDUMMY, 1,
     $                  WORK, RWORK, INFO )
         ELSE IF( JOB.EQ.2 ) THEN
            CALL CGRET( JOB, MIN( M, JJ )-RANK+1, N-RANK+1, K,
     $                  A( RANK, RANK ), LDA, C( RANK, 1 ), LDC,
     $                  WORK, RWORK, INFO )
         ELSE IF( JOB.EQ.3 ) THEN
            CALL CGRET( JOB, MIN( M, JJ )-RANK+1, N-RANK+1, K,
     $                  A( RANK, RANK ), LDA, C( 1, RANK ), LDC,
     $                  WORK, RWORK, INFO )
         END IF
      END IF
      RETURN
*
*     End of CGLBIF
*
      END
      SUBROUTINE CCNIIF( JOB, M, N, K, A, LDA, C, LDC, JPVT, VNORM,
     $                   F, RANK, PERMUT, WORK, RWORK, INFO )
*
*
*     .. Scalar Arguments ..
      INTEGER            JOB, M, N, K, LDA, LDC, RANK, INFO
      REAL               F
      LOGICAL            PERMUT
*     ..
*     .. Array Arguments ..
      COMPLEX            A( LDA, * ), C( LDC, * ), WORK( * )
      REAL               VNORM( * ), RWORK( * )
      INTEGER            JPVT( * )
*     ..
*
*  Purpose
*  =======
*
*  CCNIIF computes the "worst" column in A(1:RANK,1:RANK) and
*  determines if pivoting is necessary. If so, it pivots it into column
*  RANK, permuts vector JPVT, adjusts vector VNORM and permuts and
*  retriangularizes matrix A. It does only one permutation.
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
*  JPVT    (input/output) INTEGER array, dimension (N)
*          If JPVT(I) = K, then the Ith column of the permuted
*          A was the Kth column of the original A (just before the
*          preprocessing). If a permutation occurs, this vector will
*          be updated correctly.
*
*  VNORM   (input/output) REAL array, dimension ( N )
*          VNORM(1:RANK) contains the norms of the columns of upper
*          triangular block A(1:RANK,1:RANK). If a permutation occurs,
*          this vector will be updated correctly.
*
*  F       (input) REAL
*          F factor for the pivoting. It must be always 0 < f <= 1.
*
*  RANK    (input) INTEGER
*          The estimate for the rank. 1 <= RANK <= MIN(M,N).
*
*  PERMUT  (output) LOGICAL
*          Tells if a permutation occurred.
*
*  WORK    (workspace) COMPLEX array, dimension ( MIN(M,N) )
*
*  RWORK   (workspace) REAL array, dimension ( MIN(M,N) )
*
*  INFO    (output) INTEGER
*          = 0: Successful exit.
*          < 0: If info = -i, the i-th argument had an illegal value.
*
*  Further Details
*  ===============
*
*    If block R(1:RANK,1:RANK) is singular or near singular, there will
*  be no permutation because in that case the right (and left) singular
*  vectors are the canonical ones ((0,0,...0,1)^T).
*    That is, there will not be permutation if
*  RCOND <= SF * SLAMCH('Safe Minimum'), where SF (Safe Factor) is
*  a security factor to avoid arithmetic exceptions.
*
*  =====================================================================
*
*     .. Parameters ..
      REAL               SF
      PARAMETER          ( SF = 1.0E+2 )
      COMPLEX            CONE
      PARAMETER          ( CONE = ( 1.0E+0, 0.0E+0 ) )
*     ..
*     .. Local Scalars ..
      INTEGER            MN, JJ, J, ITEMP
      REAL               SMAX, SMIN, SMINPR, TEMP
      COMPLEX            SINE, COSINE
*     ..
*     .. Local Arrays ..
      COMPLEX            CDUMMY( 1, 1 )
*     ..
*     .. External Subroutines ..
      EXTERNAL           CCOPY, CTRSV, CHESS
*     ..
*     .. External Functions ..
      INTEGER            ISAMAX, ICAMAX
      REAL               SCNRM2, SLAMCH, SLASMX
      EXTERNAL           ISAMAX, ICAMAX, SCNRM2,
     $                   SLAMCH, SLASMX
*     ..
*     .. Executable Statements ..
*
      MN = MIN( M, N )
      INFO = 0
*
*     Quick return if possible.
*
      IF( ( MN.EQ.0 ).OR.( RANK.EQ.1 ) ) THEN
         PERMUT = .FALSE.
         RETURN
      END IF
*
*     Estimation of the largest singular value of block
*     A(1:RANK,1:RANK) by using the contents of vector
*     VNORM.
*
      ITEMP = ISAMAX( RANK, VNORM, 1 )
      SMAX = SLASMX( RANK ) * VNORM( ITEMP )
*
*     Estimation of the smallest singular value of block
*     A(1:RANK,1:RANK) and its corresponding left singular vector.
*     Save left singular vector in vector WORK(1:RANK).
*
      SMIN = ABS( A( 1, 1 ) )
      WORK( 1 ) = CONE
      DO 10 J = 2, RANK
         CALL CLAIC1( 2, J-1, WORK( 1 ), SMIN, A( 1, J ),
     $                A( J , J ), SMINPR, SINE, COSINE )
         CALL CSCAL( J-1, SINE, WORK( 1 ), 1 )
         WORK( J ) = COSINE
         SMIN = SMINPR
 10   CONTINUE
*
*     Determine if matrix A(1:RANK,1:RANK) is singular or nearly
*     singular. SF (Safe Factor) is used to say if it is singular or not.
*
      IF( SMIN.LE.( SMAX*SF*SLAMCH( 'Safe minimum' ) ) ) THEN
*
*        Singular or nearly singular matrix. Its right singular
*        vector is (0,0,...,0,1)^T. So, no pivoting is needed.
*
         PERMUT = .FALSE.
      ELSE
*
*        Follow usual method: Estimate the right singular vector
*        corresponding to the smallest singular value of upper
*        triangular block A(1:RANK,1:RANK) and store in vector
*        WORK(1:RANK).
*
         CALL CTRSV( 'Upper', 'No transpose', 'No unit',
     $                RANK, A, LDA, WORK, 1)
*
*        Find the index with largest absolute value in vector
*        WORK(1:RANK).
*
         JJ = ICAMAX( RANK, WORK, 1 )
*
*        Determine if a permutation must occur.
*
         PERMUT = ( ( JJ.LT.RANK ).AND.
     $          ( ( ABS( WORK( JJ ) )*F ).GT.ABS( WORK( RANK ) ) ) )
*
         IF( PERMUT ) THEN
*
*           Exchange cyclically to the left the colums of matrix A
*           between JJ and RANK. That is, JJ->RANK,JJ+1->JJ,...,
*           RANK->RANK-1. Use vector WORK to store temporal data.
*
            CALL CCOPY( RANK, A( 1, JJ ), 1, WORK, 1 )
            DO 20 J = JJ+1, RANK
               CALL CCOPY( J, A( 1, J ), 1, A( 1, J-1 ), 1 )
 20         CONTINUE
            CALL CCOPY( RANK, WORK, 1, A( 1, RANK ), 1 )
*
*           Exchange in the same way vector JPVT.
*
            ITEMP = JPVT( JJ )
            DO 30 J = JJ+1, RANK
                JPVT( J-1 ) = JPVT( J )
 30         CONTINUE
            JPVT( RANK ) = ITEMP
*
*           Adjust the contents of VNORM.
*
            TEMP = VNORM( JJ )
            DO 40 J = JJ+1, RANK
               VNORM( J-1 ) = VNORM( J )
 40         CONTINUE
            VNORM( RANK ) = TEMP
*
*           Retriangularize matrix A after the permutation.
*
            IF( JOB.EQ.1 ) THEN
               CALL CHESS( JOB, RANK-JJ+1, N-JJ+1, K,
     $                     A( JJ, JJ ), LDA, CDUMMY, 1,
     $                     WORK, RWORK, INFO )
            ELSE IF( JOB.EQ.2 ) THEN
               CALL CHESS( JOB, RANK-JJ+1, N-JJ+1, K,
     $                     A( JJ, JJ ), LDA, C( JJ, 1 ), LDC,
     $                     WORK, RWORK, INFO )
            ELSE IF( JOB.EQ.3 ) THEN
               CALL CHESS( JOB, RANK-JJ+1, N-JJ+1, K,
     $                     A( JJ, JJ ), LDA, C( 1, JJ ), LDC,
     $                     WORK, RWORK, INFO )
            END IF
         END IF
      END IF
      RETURN
*
*     End of CCNIIF
*
      END
      SUBROUTINE CGRET( JOB, M, N, K, A, LDA, C, LDC,
     $                       WORK, RWORK, INFO )
*
*
*     .. Scalar Arguments ..
      INTEGER            JOB, M, N, K, LDA, LDC, INFO
*     ..
*     .. Array Arguments ..
      COMPLEX            A( LDA, * ), C( LDC, * ), WORK( * )
      REAL               RWORK( * )
*     ..
*
*  Purpose
*  =======
*
*  CGRET retriangularizes a special matrix. This has the following
*  features: its first column is non-zero and its main diagonal is zero
*  except the first element. Now it is showed a 4 by 8 small example:
*                           x x x x x x x x
*                           x 0 x x x x x x
*                           x 0 0 x x x x x
*                           x 0 0 0 x x x x
*  This subroutine assumes that in all cases N>=M.
*  The transformations applied to matrix A can be also
*  applied to matrix C.
*
*  Parameters
*  ==========
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
*  A       (input/output) COMPLEX array, dimension (LDA,N)
*          On entry, the m by n matrix A.
*          On exit, the upper triangle of the array contains the
*          min(m,n) by n upper trapezoidal matrix R; the lower triangle
*          array is filled with zeros.
*
*  LDA     (input) INTEGER
*          The leading dimension of array A. LDA >= max(1,M).
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
*  WORK    (workspace) COMPLEX array, dimension ( M )
*
*  RWORK   (workspace) REAL array, dimension ( M )
*
*  INFO    (output) INTEGER
*          = 0: Successful exit.
*          < 0: If info = -i, the i-th argument had an illegal value.
*
*  =====================================================================
*
*     .. Parameters ..
      COMPLEX            CZERO
      PARAMETER          ( CZERO = ( 0.0E+0, 0.0E+0 ) )
*     ..
*     .. Common Blocks ..
      INTEGER            NB
      COMMON             /BSPRQR/ NB
*     ..
*     .. Local Scalars ..
      INTEGER            I, J, JB, ITEMP
      REAL               COSINE
      COMPLEX            R, SINE
*     ..
*     .. External Subroutines ..
      EXTERNAL           CLARTG, CROT
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MIN, CONJG
*     ..
*     .. Executable Statements ..
*
      INFO = 0
*
*     Quick return if possible.
*
      IF( ( M.EQ.0 ).OR.( M.EQ.1 ).OR.( N.EQ.0 ) )
     $   RETURN
      IF( NB.GT.1 ) THEN
*
*        Block Algorithm
*        ===============
*
*        Compute Givens Rotations needed to nullify the first column
*        of matrix A and apply on the fly to that column. Store
*        temporally the sine and cosine of the angles of the applied
*        Givens Rotations in vectors WORK and RWORK.
*
         DO 10 I = M, 2, -1
            CALL CLARTG( A( I-1, 1 ), A( I, 1 ),
     $                   RWORK( I ), WORK( I ), R )
            A( I-1, 1 ) = R
            A( I, 1 ) = CZERO
 10      CONTINUE
*
*        Apply the previously computed Givens Rotations to the rest
*        of matrix A.
*
         DO 20 J = 2, N, NB
            JB = MIN( NB, N-J+1 )
            DO 30 I = MIN( M, J+JB-1 ), J, -1
               CALL CROT( J+JB-I, A( I-1, I ), LDA, A( I, I ), LDA,
     $                    RWORK( I ), WORK( I ) )
 30         CONTINUE
            DO 40 I = MIN( M, J-1 ), 2, -1
               CALL CROT( JB, A( I-1, J ), LDA, A( I, J ), LDA,
     $                    RWORK( I ), WORK( I ) )
 40         CONTINUE
 20      CONTINUE
*
*        Update the corresponding part of matrix C.
*
         IF( ( JOB.EQ.2 ).AND.( K.GT.0 ) ) THEN
*
*           Apply the previously computed rotations from the left.
*
            DO 50 J = 1, K, NB
               JB = MIN( NB, K-J+1 )
               DO 60 I = M, 2, -1
                  CALL CROT( JB, C( I-1, J ), LDC, C( I, J ), LDC,
     $                      RWORK( I ), WORK( I ) )
 60            CONTINUE
 50         CONTINUE
         ELSE IF( ( JOB.EQ.3 ).AND.( K.GT.0 ) ) THEN
*
*           Apply the transpose of the previously computed rotations
*           from the right.
*
            DO 70 I = M, 2, -1
               CALL CROT( K, C( 1, I-1 ), 1, C( 1, I ), 1,
     $                   RWORK( I ), CONJG( WORK( I ) ) )
 70         CONTINUE
         END IF
      ELSE
*
*        Non-Block Algorithm
*        ===================
*
         DO 90 I = M, 2, -1
            ITEMP = I - 1
*
*           Compute the rotation parameters and update column 1 of A.
*
            CALL CLARTG( A( ITEMP, 1 ), A( I , 1 ), COSINE, SINE, R )
            A( ITEMP, 1 ) = R
            A( I, 1 ) = CZERO
*
*           Update columns I:N of matrix A.
*
            CALL CROT( N-I+1, A( ITEMP, I ), LDA, A( I, I ), LDA,
     $                 COSINE, SINE )
*
*           Update the corresponding part of matrix C.
*
            IF( ( JOB.EQ.2 ).AND.( K.GT.0 ) ) THEN
*
*              Apply the previously computed rotations from the left.
*
               CALL CROT( K, C( ITEMP, 1 ), LDC, C( I, 1 ), LDC,
     $                     COSINE, SINE )
            ELSE IF( ( JOB.EQ.3 ).AND.( K.GT.0 ) ) THEN
*
*              Apply the transpose of the previously computed rotations
*              from the right.
*
               CALL CROT( K, C( 1, ITEMP ), 1, C( 1, I ), 1,
     $                     COSINE, CONJG( SINE ) )
            END IF
 90      CONTINUE
      END IF
      RETURN
*
*     End of CGRET
*
      END
      SUBROUTINE CHESS( JOB, M, N, K, A, LDA, C, LDC,
     $                       WORK, RWORK, INFO )
*
*
*     .. Scalar Arguments ..
      INTEGER            JOB, M, N, K, LDA, LDC, INFO
*     ..
*     .. Array Arguments ..
      COMPLEX            A( LDA, * ), C( LDC, * ), WORK( * )
      REAL               RWORK( * )
*     ..
*
*  Purpose
*  =======
*
*  CHESS reduces the upper hessemberg matrix A to upper triangular form.
*  applied to matrix C if argument JOB asks.
*  This subroutine assumes that in all cases N>=M.
*
*  Parameters
*  ==========
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
*  A       (input/output) COMPLEX array, dimension (LDA,N)
*          On entry, the m by n matrix A.
*          On exit, the upper triangle of the array contains the
*          min(m,n) by n upper trapezoidal matrix R; the lower triangle
*          array is filled with zeros.
*
*  LDA     (input) INTEGER
*          The leading dimension of array A. LDA >= max(1,M).
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
*  WORK    (workspace) COMPLEX array, dimension ( M )
*
*  RWORK   (workspace) REAL array, dimension ( M )
*
*  INFO    (output) INTEGER
*          = 0: Successful exit.
*          < 0: If info = -i, the i-th argument had an illegal value.
*
*  =====================================================================
*
*     .. Parameters ..
      COMPLEX            CZERO
      PARAMETER          ( CZERO = ( 0.0E+0, 0.0E+0 ) )
*     ..
*     .. Common Blocks ..
      INTEGER            NB
      COMMON             /BSPRQR/ NB
*     ..
*     .. Local Scalars ..
      INTEGER            I, J, ITEMP, JB
      REAL               COSINE
      COMPLEX            R, SINE
*     ..
*     .. External Subroutines ..
      EXTERNAL           CLARTG, CROT
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MIN, CONJG
*     ..
*     .. Executable Statements ..
*
      INFO = 0
*
*     Quick return if possible.
*
      IF( ( M.EQ.0 ).OR.( M.EQ.1 ).OR.( N.EQ.0 ) )
     $   RETURN
      IF( NB.GT.1 ) THEN
*
*        Block Algorithm
*        ===============
*
*        Compute Givens Rotations needed to reduce upper hessenberg
*        matrix A to triangular form, and apply on the fly those
*        rotations to matrix. Store temporally the sine and cosine
*        of the angles of the applied Givens Rotations in
*        vectors WORK and RWORK.
*
         DO 10 J = 1, N, NB
            JB = MIN( NB, N-J+1 )
            DO 20 I = 2, MIN( M, J )
               CALL CROT( JB, A( I-1, J ), LDA, A( I, J ), LDA,
     $                    RWORK( I ), WORK( I ) )
 20         CONTINUE
            DO 30 I = J+1, MIN( M, J+JB )
               ITEMP = I-1
               CALL CLARTG( A( ITEMP, ITEMP ), A( I, ITEMP ),
     $                      RWORK( I ), WORK( I ), R )
               A( ITEMP, ITEMP ) = R
               A( I, ITEMP ) = CZERO
               CALL CROT( J+JB-I, A( ITEMP, I ), LDA,
     $                    A( I, I ), LDA, RWORK( I ), WORK( I ) )
 30         CONTINUE
 10      CONTINUE
*
*        Update the corresponding part of matrix C.
*
         IF( ( JOB.EQ.2 ).AND.( K.GT.0 ) ) THEN
*
*           Apply the previously computed rotations from the left.
*
            DO 40 J = 1, K, NB
               JB = MIN( NB, K-J+1 )
               DO 50 I = 2, M
                  CALL CROT( JB, C( I-1, J ), LDC, C( I, J ), LDC,
     $                       RWORK( I ), WORK( I ) )
 50            CONTINUE
 40         CONTINUE
         ELSE IF( ( JOB.EQ.3 ).AND.( K.GT.0 ) ) THEN
*
*           Apply the transpose of the previously computed rotations
*           from the right.
*
            DO 60 I = 2, M
               CALL CROT( K, C( 1, I-1 ), 1, C( 1, I ), 1,
     $                    RWORK( I ), CONJG( WORK( I ) ) )
 60         CONTINUE
         END IF
      ELSE
*
*        Non-Block Algorithm
*        ===================
*
         DO 80 I = 2, M
            ITEMP = I - 1
*
*           Compute the rotation parameters.
*
            CALL CLARTG( A( ITEMP, ITEMP ), A( I, ITEMP ),
     $                   COSINE, SINE, R )
*
*           Update columns I-1:N of matrix A.
*
            A( ITEMP, ITEMP ) = R
            A( I, ITEMP ) = CZERO
            CALL CROT( N-I+1, A( ITEMP, I ), LDA, A( I, I ), LDA,
     $                 COSINE, SINE )
*
*           Update the corresponding part of matrix C.
*
            IF( ( JOB.EQ.2 ).AND.( K.GT.0 ) ) THEN
*
*              Apply the previously computed rotations from the left.
*
               CALL CROT( K, C( ITEMP, 1 ), LDC, C( I, 1 ), LDC,
     $                    COSINE, SINE )
            ELSE IF( ( JOB.EQ.3 ).AND.( K.GT.0 ) ) THEN
*
*              Apply the transpose of the previously computed rotations
*              from the right.
*
               CALL CROT( K, C( 1, ITEMP ), 1, C( 1, I ), 1,
     $                    COSINE, CONJG( SINE ) )
            END IF
 80      CONTINUE
      END IF
      RETURN
*
*     End of CHESS
*
      END
