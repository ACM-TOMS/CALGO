      SUBROUTINE DLAQP3( M, N, A, LDA, JPVT, RCOND, RANK, TAU, WORK,
     $                   LWORK, INFO )
*
*     This code is part of a package for solving rank deficient least
*     squares problems, written by:
*     ==================================================================
*     L. Foster                   and   R. Kommu
*     Department of Mathematics         Department of Physics
*     San Jose State University         San Jose State University
*     San Jose, CA 95192                San Jose, CA 95192
*     foster@math.sjsu.edu              rkommu@email.sjsu.edu
*     ==================================================================
*     03/05/2004
*
*     .. Scalar Arguments ..
      INTEGER            INFO, LDA, LWORK, M, N, RANK
      DOUBLE PRECISION   RCOND
*     ..
*     .. Array Arguments ..
      INTEGER            JPVT( * )
      DOUBLE PRECISION   A( LDA, * ), TAU( * ), WORK( * )
*     ..
*
*  Purpose
*  =======
*
*  DLAQP3 uses Level 3 BLAS to compute enough of a QR factorization
*  with column pivoting so that P, Q, R11 and R12 can be determined in:
*      A * P = Q * [ R11 R12 ]
*                  [  0  R22 ]
*  with R11 defined as the largest leading submatrix whose estimated
*  condition number is less than 1/RCOND.  The order of R11, RANK,
*  is the effective rank of A.
*
*  Arguments
*  =========
*
*  M       (input) INTEGER
*          The number of rows of the matrix A. M >= 0.
*
*  N       (input) INTEGER
*          The number of columns of the matrix A.  N >= 0.
*
*  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
*          On entry, the M-by-N matrix A.
*          On exit, in the first RANK columns, the upper triangular
*          part contains the upper triangular matrix R11 and the
*          elements of A below the diagonal, together with TAU,
*          represent the orthogonal matrix Q as a product of
*          RANK elementary reflections.  The submatrix R12 is in
*          the rows 1 to RANK and columns RANK + 1 to N of A.
*          The portion of A in rows RANK + 1 to N and columns
*          RANK + 1 to M contains information related to R22 that
*          is used in the calculation of P, Q, R11 and R12.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A. LDA >= max(1,M).
*
*  JPVT    (input/output) INTEGER array, dimension (N)
*          On entry, if JPVT(J).ne.0, the J-th column of A is permuted
*          to the front of A*P (a leading column); if JPVT(J)=0,
*          the J-th column of A is a free column.  The determination
*          of the effective rank, RANK, is more reliable if JPVT and
*          RCOND are chosen so that the estimated condition number of
*          the fixed columns is less than 1 / RCOND.
*          On exit, if JPVT(J)=K, then the J-th column of A*P was the
*          the K-th column of A.
*
*  RCOND   (input) DOUBLE PRECISION
*          RCOND is used to determine the effective rank of A, which
*          is defined as the order of the largest leading triangular
*          submatrix R11 in the QR factorization with pivoting of A,
*          whose estimated condition number < 1/RCOND.
*
*  RANK    (output) INTEGER
*          The effective rank of A, i.e., the order of the submatrix
*          R11.
*
*  TAU     (output) DOUBLE PRECISION array, dimension (min(M,N))
*          The scalar factors of the elementary reflectors.
*
*  WORK    (workspace/output) DOUBLE PRECISION array, dimension (LWORK)
*          On exit, if INFO=0, WORK(1) returns the optimal LWORK.
*
*  LWORK   (input) INTEGER
*          The dimension of the array WORK. LWORK >= 2*MN + 3*N + 1,
*          where MN = min( M, N ).
*          For optimal performance LWORK >= 2*MN + 2*N +( N+1 )*NB,
*          where NB is the optimal blocksize.
*
*          If LWORK = -1, then a workspace query is assumed; the routine
*          only calculates the optimal size of the WORK array, returns
*          this value as the first entry of the WORK array, and no error
*          message related to LWORK is issued by XERBLA.
*
*  INFO    (output) INTEGER
*          = 0: successful exit.
*          < 0: if INFO = -i, the i-th argument had an illegal value.
*
*  Further Details
*  ===============
*
*  The matrix Q is represented as a product of elementary reflectors
*
*     Q = H(1) H(2) . . . H(k), where k = RANK.
*
*  Each H(i) has the form
*
*     H(i) = I - tau * v * v'
*
*  where tau is a real/complex scalar, and v is a real/complex vector
*  with v(1:i-1) = 0 and v(i) = 1; v(i+1:m) is stored on exit in
*  A(i+1:m,i), and tau in TAU(i).
*
*  This is a modification of LAPACK routine xGEQP3 written by
*    G. Quintana-Orti, Depto. de Informatica, Universidad Jaime I, Spain
*    X. Sun, Computer Science Dept., Duke University, USA
*  Modified by L. Foster, Department of Mathematics, San Jose State
*    University, San Jose, CA.
*
*  =====================================================================
*
*     .. Parameters ..
      INTEGER            INB, INBMIN, IXOVER
      PARAMETER          ( INB = 1, INBMIN = 2, IXOVER = 3 )
      INTEGER            IMAX, IMIN
      PARAMETER          ( IMAX = 1, IMIN = 2 )
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            LQUERY
      INTEGER            FJB, I, ISMIN, ISMAX, IWS, J, JB, LWKOPT, MN,
     $                   MINWS, NA, NB, NBMIN, NFXD, NX, SM, SMINMN, SN,
     $                   TOPBMN
      DOUBLE PRECISION   C1, C2, S1, S2, SMAX, SMAXPR, SMIN, SMINPR
*
*     ..
*     .. External Subroutines ..
      EXTERNAL           DGEQRF, DLAIC1, DLAQPP, DLAQPQ, DORMQR, DSWAP,
     $                   XERBLA
*     ..
*     .. External Functions ..
      INTEGER            ILAENV
      DOUBLE PRECISION   DNRM2
      EXTERNAL           ILAENV, DNRM2
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, DBLE, INT, MAX, MIN
*     ..
*     .. Executable Statements ..
*
      MN = MIN( M, N )
      IWS = 2*MN + 3*N + 1
      RANK = 0
*
*     Test input arguments
*     ====================
*
      INFO = 0
      NB = ILAENV( INB, 'DGEQRF', ' ', M, N, -1, -1 )
      LWKOPT = 2*N + ( N+1 )*NB + 2*MN
      WORK( 1 ) = LWKOPT
      LQUERY = ( LWORK.EQ.-1 )
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -4
      ELSE IF( ( LWORK.LT.IWS ) .AND. .NOT.LQUERY ) THEN
         INFO = -10
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DLAQP3', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF
*
*     Quick return if possible.
*
      IF( MN.EQ.0 ) THEN
         WORK( 1 ) = 1
         RETURN
      END IF
*
*     Move initial columns up front.
*
      NFXD = 1
      DO 10 J = 1, N
         IF( JPVT( J ).NE.0 ) THEN
            IF( J.NE.NFXD ) THEN
               CALL DSWAP( M, A( 1, J ), 1, A( 1, NFXD ), 1 )
               JPVT( J ) = JPVT( NFXD )
               JPVT( NFXD ) = J
            ELSE
               JPVT( J ) = J
            END IF
            NFXD = NFXD + 1
         ELSE
            JPVT( J ) = J
         END IF
   10 CONTINUE
      NFXD = NFXD - 1
*
*     Factorize fixed columns
*     =======================
*
*     Compute the QR factorization of fixed columns and update
*     remaining columns.
*
      IF( NFXD.GT.0 ) THEN
         NA = MIN( M, NFXD )
*CC      CALL DGEQR2( M, NA, A, LDA, TAU, WORK, INFO )
         CALL DGEQRF( M, NA, A, LDA, TAU, WORK, LWORK, INFO )
         IWS = MAX( IWS, INT( WORK( 1 ) ) )
         IF( NA.LT.N ) THEN
*CC         CALL DORM2R( 'Left', 'Transpose', M, N-NA, NA, A, LDA,
*CC  $                   TAU, A( 1, NA+1 ), LDA, WORK, INFO )
            CALL DORMQR( 'Left', 'Transpose', M, N-NA, NA, A, LDA, TAU,
     $                   A( 1, NA+1 ), LDA, WORK, LWORK, INFO )
            IWS = MAX( IWS, INT( WORK( 1 ) ) )
         END IF
*
*        Determine RANK of the fixed columns using incremental
*        condition estimation.  The determination of the effective
*        rank, RANK, is more reliable if JPVT and RCOND are chosen
*        so that the estimated condition number of the fixed columns
*        is less than 1 / RCOND.
*
         ISMIN = 1
         ISMAX = MN + 1
         WORK( ISMIN ) = ONE
         WORK( ISMAX ) = ONE
         SMAX = ABS( A( 1, 1 ) )
         SMIN = SMAX
         IF( ABS( A( 1, 1 ) ).EQ.ZERO ) THEN
*           Can exit early, partial factorization done since the
*           calculated rank is zero
            RANK = 0
            WORK( 1 ) = DBLE( LWKOPT )
            RETURN
         ELSE
            RANK = 1
         END IF
*
   20    CONTINUE
         IF( RANK.LT.NA ) THEN
            I = RANK + 1
            CALL DLAIC1( IMIN, RANK, WORK( ISMIN ), SMIN, A( 1, I ),
     $                   A( I, I ), SMINPR, S1, C1 )
            CALL DLAIC1( IMAX, RANK, WORK( ISMAX ), SMAX, A( 1, I ),
     $                   A( I, I ), SMAXPR, S2, C2 )
*
            IF( SMAXPR*RCOND.LE.SMINPR ) THEN
               DO 30 I = 1, RANK
                  WORK( ISMIN+I-1 ) = S1*WORK( ISMIN+I-1 )
                  WORK( ISMAX+I-1 ) = S2*WORK( ISMAX+I-1 )
   30          CONTINUE
               WORK( ISMIN+RANK ) = C1
               WORK( ISMAX+RANK ) = C2
               SMIN = SMINPR
               SMAX = SMAXPR
               RANK = RANK + 1
               GO TO 20
            END IF
*           Can exit early, due to the RCOND test the partial
*           factorization is done
            WORK( 1 ) = DBLE( LWKOPT )
            RETURN
         END IF
*
      END IF
*
*
*     Factorize free columns
*     ======================
*
      IF( NFXD.LT.MN ) THEN
*
         SM = M - NFXD
         SN = N - NFXD
         SMINMN = MN - NFXD
*
*        Determine the block size.
*
         NB = ILAENV( INB, 'DGEQRF', ' ', SM, SN, -1, -1 )
         NBMIN = 2
         NX = 0
*
         IF( ( NB.GT.1 ) .AND. ( NB.LT.SMINMN ) ) THEN
*
*           Determine when to cross over from blocked to unblocked code.
*
            NX = MAX( 0, ILAENV( IXOVER, 'DGEQRF', ' ', SM, SN, -1,
     $           -1 ) )
*
*
            IF( NX.LT.SMINMN ) THEN
*
*              Determine if workspace is large enough for blocked code.
*
*              Based on the do 40 loop, below, MINWS should be:
               MINWS = 2*N + ( SN+1 )*NB + 2*MN
               IWS = MAX( IWS, MINWS )
               IF( LWORK.LT.MINWS ) THEN
*
*                 Not enough workspace to use optimal NB: Reduce NB and
*                 determine the minimum value of NB.
*
                  NB = ( LWORK-2*N-2*MN ) / ( SN+1 )
                  NBMIN = MAX( 2, ILAENV( INBMIN, 'DGEQRF', ' ', SM, SN,
     $                    -1, -1 ) )
*
*
               END IF
            END IF
         END IF
*
*        Initialize partial column norms. The elements 2*MN + 1 through
*        2*MN + N of work store the exact column norms.
*
         DO 40 J = NFXD + 1, N
            WORK( 2*MN+J ) = DNRM2( SM, A( NFXD+1, J ), 1 )
            WORK( 2*MN+N+J ) = WORK( 2*MN+J )
   40    CONTINUE
*
         IF( ( NB.GE.NBMIN ) .AND. ( NB.LT.SMINMN ) .AND.
     $       ( NX.LT.SMINMN ) ) THEN
*
*           Use blocked code initially.
*
            J = NFXD + 1
*
*           Compute factorization: while loop.
*
            TOPBMN = MN - NX
   50       CONTINUE
            IF( J.LE.TOPBMN ) THEN
               JB = MIN( NB, TOPBMN-J+1 )
*
*              Factorize JB columns among columns J:N.
*
               CALL DLAQPP( M, N-J+1, J-1, JB, FJB, A( 1, J ), LDA,
     $                      JPVT( J ), RCOND, RANK, TAU( J ), SMIN,
     $                      WORK( 1 ), SMAX, WORK( MN+1 ),
     $                      WORK( 2*MN+J ), WORK( 2*MN+N+J ),
     $                      WORK( 2*MN+2*N+1 ), WORK( 2*MN+2*N+JB+1 ),
     $                      N-J+1 )
*
               J = J + FJB
               IF( RANK.LT.( J-1 ) ) THEN
*                 Can exit early, due to the RCOND test the partial
*                 factorization is done
                  WORK( 1 ) = IWS
                  RETURN
               END IF
               GO TO 50
            END IF
         ELSE
            J = NFXD + 1
         END IF
*
*        Use unblocked code to factor the last or only block.
*
         IF( J.LE.MN )
     $      CALL DLAQPQ( M, N-J+1, J-1, A( 1, J ), LDA, JPVT( J ),
     $                   RCOND, RANK, TAU( J ), SMIN, WORK( 1 ), SMAX,
     $                   WORK( MN+1 ), WORK( 2*MN+J ), WORK( 2*MN+N+J ),
     $                   WORK( 2*MN+2*N+1 ) )
*
      END IF
      J = MN
*
      WORK( 1 ) = IWS
      RETURN
*
*     End of DLAQP3
*
      END
