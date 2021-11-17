      SUBROUTINE ZGEQPW( M, LWSIZE, NB, OFFSET, LACPTD, A, LDA, JPVT,
     $                    IRCOND, X, SMIN, MXNM, TAU, WORK, RWORK )
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
*     $Date: 96/12/30 16:59:38 $
*
*     .. Scalar Arguments ..
      INTEGER            M, LWSIZE, NB, OFFSET, LACPTD, LDA
      DOUBLE PRECISION   IRCOND, SMIN, MXNM
*     ..
*     .. Array Arguments ..
      INTEGER            JPVT( * )
      COMPLEX*16         A( LDA, * ), TAU( * ), X( * ), WORK( * )
      DOUBLE PRECISION   RWORK( * )
*
*
*  Purpose
*  =======
*
*  ZGEQPW applies one block step of the Householder QR
*  factorization algorithm with restricted pivoting. It is called
*  by ZGEQPB to factorize a window of the matrix.
*
*  Let A be the partial QR factorization of an M by (OFFSET+LWSIZE)
*  matrix C, i.e. we have computed an orthogonal matrix Q1 and a
*  permutation matrix P1 such that
*            C * P1 = Q1 * A
*  and A(:,1:OFFSET) is upper triangular. Let us denote A(:,1:OFFSET)
*  by B. Then in addition let
*  X be an approximate smallest left singular vector of B in the sense
*  that
*       sigma_min(B) ~ twonorm(B'*X) = SMIN
*  and
*       sigma_max(B) ~ ((offset)**(1./3.))*MXNM = SMAX
*  with
*       cond_no(B) ~ SMAX/SMIN <= 1/IRCOND
*
*  Then ZGEQP2 tries to identify NB columns in
*   A(:,OFFSET+1:OFFSET+LWSIZE) such that
*       cond_no([B,D]) < 1/IRCOND
*  where D are the KB columns of A(:,OFFSET+1:OFFSET+LWSIZE) that were
*  considered independent with respect to the threshold 1/IRCOND.
*
*  On exit,
*       C * P2 = Q2 * A
*  is again a partial QR factorization of C, but columns
*  OFFSET+1:OFFSET+LACPTD of A have been reduced via
*  a series of elementary reflectors to upper
*  trapezoidal form. Further
*       sigma_min(A(:,1:OFFSET+LACPTD))
*                ~ twonorm(A(:,1:OFFSET+LACPTD)'*x) = SMIN
*  and
*       sigma_max(A(:,1:OFFSET+LACPTD)) ~ sqrt(OFFSET+LACPTD)*MXNM = SMAX
*  with
*       cond_no(A(:,1:OFFSET+LACPTD))
*                   ~ SMAX/SMIN <= 1/IRCOND.
*
*  In the ideal case, LACPTD = NB, that is,
*  we found NB independent columns in the window consisting of
*  the first LWSIZE columns of A.
*
*
*  Arguments
*  =========
*
*  M       (input) INTEGER
*          The number of rows of the matrix A.
*
*  LWSIZE  (input) INTEGER
*          The size of the pivot window in A.
*
*  NB      (input) INTEGER
*          The number of independent columns one would like to identify.
*          This equals the desired blocksize in ZGEQPB.
*
*  OFFSET  (input) INTEGER
*          The number of rows and columns of A that need not be updated.
*
*  LACPTD  (output) INTEGER
*          The number of columns in A(:,OFFSET+LWSIZE) that were
*          accepted as linearly independent.
*
*  A       (input/output) COMPLEX*16 array, dimension (LDA,OFFSET+LWSIZE)
*          On entry, the upper triangle of A(:,1:OFFSET) contains the
*          partially completed triangular factor R; the elements below
*          the diagonal, with the array TAU, represent the matrix Q1 as
*          a product of elementary reflectors.
*          On exit, the upper triangle of A(:,OFFSET+LACPTD) contains
*          the partially completed upper triangular factor R; the
*          elements below the diagonal, with the array TAU, represent
*          the matrix Q2 as a product of elementary reflectors.
*          A(OFFSET:M,LACPTD+1:LWSIZE) has been updated by the product
*          of these elementary reflectors.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A. LDA >= M.
*
*  JPVT    (input/output) INTEGER array, dimension (OFFSET+LWSIZE)
*          On entry and exit, jpvt(i) = k if the i-th column
*          of A was the k-th column of C.
*
*  IRCOND (input) DOUBLE PRECISION
*          1/IRCOND is the threshold for the condition number.
*
*  X       (input/output) COMPLEX*16 array, dimension (OFFSET+NB)
*          On entry, X(1:OFFSET) is an approximate left nullvector of
*          the upper triangle of A(1:OFFSET,1:OFFSET).
*          On exit, X(1:OFFSET+LACPTD) is an approximate left
*          nullvector of the matrix in the upper triangle of
*          A(1:OFFSET+LACPTD,1:OFFSET+LACPTD).
*
*  SMIN    (input/output) DOUBLE PRECISION
*          On entry, SMIN is an estimate for the smallest singular
*          value of the upper triangle of A(1:OFFSET,1:OFFSET).
*          On exit, SMIN is an estimate for the smallest singular
*          value of the matrix in the upper triangle of
*          A(1:OFFSET+LACPTD,1:OFFSET+LACPTD).
*
*  MXNM    (input) FLOATING_DECLARE
*          The norm of the largest column in matrix A.
*
*  TAU     (output) COMPLEX*16 array, dimension (OFFSET+LWSIZE)
*          On exit, TAU(1:OFFSET+LACPTD) contains details of
*          the matrix Q2.
*
*  WORK    (workspace) COMPLEX*16 array, dimension (LWSIZE)
*
*  RWORK   (workspace) DOUBLE PRECISION array, dimension (2*LWSIZE)
*
*  ================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            I, K, I1, PVTIDX, LASTK
      DOUBLE PRECISION   TEMP, TEMP2, SMAX
      COMPLEX*16         GAMMA, AKK
*     ..
*     .. External Subroutines ..
      EXTERNAL           DZNRM2, ZSCAL, ZSWAP, ZLARFG,
     $                   ZLARF, IDAMAX, ZLAUC1, DLASMX
      INTEGER            IDAMAX
      DOUBLE PRECISION   DZNRM2, DLASMX
      LOGICAL            ZLAUC1
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, MIN, SQRT, DBLE, DCMPLX, DCONJG
*     ..
*     .. Executable Statements ..
*
*     Initialize partial column norms (stored in the first LWSIZE
*     entries of WORK) and exact column norms (stored in the second
*     LWSIZE entries of WORK) for the first batch of columns.
*
      DO 10 I = 1,LWSIZE
         RWORK( I ) = DZNRM2( M-OFFSET, A( OFFSET+1, OFFSET+I ), 1 )
         RWORK( LWSIZE+I ) = RWORK( I )
10    CONTINUE
*
*     *************
*     * Main loop *
*     *************
*
      LASTK = MIN( M, OFFSET+LWSIZE )
      LACPTD = 0
1000  IF( LACPTD.EQ.NB ) GOTO 2000
*
*        Determine pivot candidate.
*        =========================
         PVTIDX = OFFSET + LACPTD +
     $            IDAMAX( LWSIZE-LACPTD, RWORK( LACPTD+1 ), 1 )
         K = OFFSET + LACPTD + 1
*
*        Exchange current column and pivot column.
*
         IF( PVTIDX.NE.K ) THEN
            CALL ZSWAP( M, A( 1, PVTIDX ), 1, A( 1, K ), 1 )
            I1 = JPVT( PVTIDX )
            JPVT( PVTIDX ) = JPVT( K )
            JPVT( K ) = I1
            TEMP = RWORK( PVTIDX-OFFSET )
            RWORK( PVTIDX-OFFSET ) = RWORK( K-OFFSET )
            RWORK( K-OFFSET ) = TEMP
            TEMP = RWORK( PVTIDX-OFFSET+LWSIZE )
            RWORK( PVTIDX-OFFSET+LWSIZE ) = RWORK( K+LWSIZE-OFFSET )
            RWORK( K+LWSIZE-OFFSET ) = TEMP
         END IF
*
*        Determine (offset+lacptd+1)st diagonal element
*        GAMMA of matrix A should elementary reflector be applied.
*
         TEMP = DBLE( A( K, K ) )
         IF( TEMP.EQ.ZERO ) THEN
            GAMMA = -RWORK( K-OFFSET )
         ELSE
            GAMMA = -SIGN( RWORK( K-OFFSET ), TEMP )
         END IF
*
*        Update estimate for largest singular value.
*
         SMAX = DLASMX( K )*MXNM
*
*        Is candidate pivot column acceptable ?
*        =====================================
         IF( ZLAUC1( K, X, SMIN, A( 1, K ), GAMMA, SMAX*IRCOND ) )
     $      THEN
*
*           Pivot candidate was accepted.
*           ============================
*
            LACPTD = LACPTD + 1
*
*           Generate Householder vector.
*
            IF( K.LT.M ) THEN
               CALL ZLARFG( M-K+1, A( K, K ), A( K+1, K ), 1,
     $                      TAU( K ) )
            ELSE
               CALL ZLARFG( 1, A( M, K), A( M, K ), 1, TAU ( K ) )
            END IF
*
*           Apply Householder reflection to A(k:m,k+1:lwsize).
*
            IF( LACPTD.LT.LWSIZE ) THEN
               AKK = A( K, K )
               A( K, K ) = DCMPLX( ONE )
               CALL ZLARF( 'Left', M-K+1, LWSIZE-LACPTD,
     $                     A( K, K ), 1, DCONJG( TAU( K ) ),
     $                     A( K, K+1 ), LDA, WORK )
               A( K, K ) = AKK
            END IF
*
*           Update partial column norms.
*
            IF( K.LT.LASTK ) THEN
               DO 20 I = LACPTD+1,LWSIZE
                  IF( RWORK( I ).NE.ZERO ) THEN
                     TEMP = ONE-
     $                         ( ABS( A( K, OFFSET+I ) )/RWORK( I ) )**2
                     TEMP = MAX( TEMP, ZERO )
                     TEMP2 = ONE+ 0.05*TEMP*
     $                               ( RWORK( I )/RWORK( I+LWSIZE ) )**2
                     IF( TEMP2.EQ.ONE ) THEN
                        RWORK( I ) = DZNRM2( M-K,
     $                                         A( K+1, OFFSET+I ), 1 )
                        RWORK( I+LWSIZE ) = RWORK( I )
                     ELSE
                        RWORK( I ) = RWORK( I )*SQRT( TEMP )
                     END IF
                  END IF
20             CONTINUE
            END IF
         ELSE
*
*           Reject all remaining columns in pivot window.
*           ============================================
*
            GOTO 2000
         END IF
*
*     End while.
*
      GOTO 1000
2000  CONTINUE
      RETURN
*
*     End of ZGEQPW
*
      END
