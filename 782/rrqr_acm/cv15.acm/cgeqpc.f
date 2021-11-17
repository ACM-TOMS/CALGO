      SUBROUTINE CGEQPC( JOB, M, N, K, A, LDA, C, LDC, DSRD, OFFSET,
     $                  IRCOND, LACPTD, JPVT, TAU, X, SVLUES, MXNM,
     $                  WORK, LWORK, RWORK )
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
*     $Date: 96/12/30 16:59:37 $
*
*     .. Scalar Arguments ..
      INTEGER            JOB, M, N, K, LDA, LDC, DSRD, OFFSET, LACPTD,
     $                   LWORK
      REAL               IRCOND, MXNM
*     ..
*     .. Array Arguments ..
      INTEGER            JPVT( * )
      COMPLEX            A( LDA, * ), C( LDC, * ), TAU( * ),
     $                   WORK( * ), X( * )
      REAL               SVLUES( 4 ), RWORK( * )
*     ..
*
*  Purpose:
*  =======
*
*  CGEQPC continues a partial QR factorization of A. If
*  A(1:OFFSET,1:OFFSET) has been reduced to upper triangular
*  form, then SGQPC applies the traditional column pivoting
*  strategy to identify DSRD more independent columns of A with
*  the restriction that the condition number of the leading
*  triangle of A should not be larger than 1/IRCOND.  If
*  LACPTD ( <= DSRD) such columns are found, then the condition
*  number of
*     A(1:OFFSET+LACPTD,1:OFFSET+LACPTD) is less than 1/IRCOND.
*  If LACPTD < DSRD, then the QR factorization of A is completed,
*  otherwise only DSRD new steps were performed.
*
*  Arguments:
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
*  DSRD    (input) INTEGER
*          The number of independent columns one would like to
*          extract.
*
*  OFFSET  (input) INTEGER
*          A(1:OFFSET,1:OFFSET) has already been factored.
*          OFFSET >= 0.
*
*  IRCOND  (input) REAL
*          1/IRCOND is threshold for condition number.
*
*  LACPTD  (output) INTEGER
*          The number of additional columns that were identified
*          as independent.
*
*  JPVT    (input/output) INTEGER array, dimension (N)
*          If JPVT(I) = K, then the Ith column of the permuted
*          A was the Kth column of the original A.
*
*  TAU     (input/output) COMPLEX array, dimension (MIN(M,N))
*          Further details of the matrix Q (see A).
*
*  X       (input/output) COMPLEX array, dimension (MIN(M,N))
*          On entry: X(1:OFFSET) contains an approximate smallest
*          left singular vector of A(1:OFFSET,1:OFFSET)
*          On exit: X(1:OFFSET+LACPTD) contains an approximate
*          smallest left singular vector of
*          A(1:OFFSET+LACPTD,1:OFFSET+LACPTD).
*
*  SVLUES  (input/output) REAL array, dimension(4)
*          estimates of singular values.
*          On entry: SVLUES(1) = sigma_max(A(1:M,1:N))
*                    SVLUES(2) = sigma_min(A(1:OFFSET,1:OFFSET))
*          On exit: SVLUES(1) = sigma_max(A(1:M,1:N))
*                   SVLUES(2) = sigma_r(B)
*                   SVLUES(3) = sigma_(min(r+1,min(m,n)))(B)
*                   SVLUES(4) = sigma_min(A)
*          where r = OFFSET+LACPTD and B = A(1:r,1:r)
*
*  MXNM    (input/output) FLOATING_DECLARE
*          On entry: norm of largest column in A(1:OFFSET,1:OFFSET)
*          On exit: norm of largest column in
*                   A(1:J,1:J) where J = OFFSET+LACPTD
*
*  WORK    (workspace) FLOATING_DECLARE array, dimension (LWORK)
*
*  LWORK   (input) INTEGER
*             MAX( 1, N*NB )              if JOB=1, or
*             MAX( 1, MAX( N, K )*NB )    otherwise.
*          where NB is the maximum of blocksize used within xGEQRF and
*          blocksize used within xUNMQR.
*
*  RWORK   (workspace) REAL array, dimension ( 2*N ).
*
*  Further Details
*  ===============
*
*  The matrix Q is represented as a product of elementary reflectors
*
*     Q = H(1) H(2) . . . H(n)
*
*  Each H(i) has the form
*
*     H = I - tau * v * v'
*
*  where tau is a complex scalar, and v is a complex vector with
*  v(1:i-1) = 0 and v(i) = 1; v(i+1:m) is stored on exit in A(i+1:m,i).
*
*  The matrix P is represented in jpvt as follows: If
*     jpvt(j) = i
*  then the jth column of P is the ith canonical unit vector.
*
*  =====================================================================
*
*     .. Parameters ..
      REAL               ZERO, ONE
      COMPLEX            CONE
      PARAMETER          ( ZERO = 0.0E+0, ONE = 1.0E+0,
     $                   CONE = ( 1.0E+0, 0.0E+0 ) )
*
*     Indices into the 'svlues' array.
*
      INTEGER            IMAX, IBEFOR, IAFTER, IMIN
      PARAMETER          ( IMAX = 1, IBEFOR = 2, IAFTER = 3, IMIN = 4 )
*     ..
*     .. Local Scalars ..
      INTEGER            I, J, PVT, MN, ITEMP, INFO, LASTI
      REAL               TEMP, TEMP2, SMIN, SMINPR, SMAX, SMAXPR
      COMPLEX            AII, SINE, COSINE
*     ..
*     .. External Subroutines ..
      EXTERNAL           CLARFG, CLARF, CSWAP, CSCAL,
     $                   CLAIC1, CGEQRF, CUNMQR
*     ..
*     .. External Functions ..
      INTEGER            ISAMAX
      REAL               SCNRM2, SLASMX
      LOGICAL            CLAUC1
      EXTERNAL           ISAMAX, SCNRM2, SLASMX, CLAUC1
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, MIN, CONJG
*     ..
*     .. Executable Statements ..
*
      MN = MIN( M, N )
      LACPTD = 0
      IF( OFFSET.GT.0 ) THEN
         SMAX = SVLUES( IMAX )
         SMIN = SVLUES( IBEFOR )
      END IF
*
*     Initialize partial column norms. The first n entries of
*     work store the exact column norms.
*
      DO 10 I = OFFSET+1,N
         RWORK( I ) = SCNRM2( M-OFFSET, A( OFFSET+1, I ), 1 )
         RWORK( N+I ) = RWORK( I )
10    CONTINUE
*
*     Compute factorization.
*
      LASTI = MIN( MN, OFFSET+DSRD )
      DO 20 I = OFFSET+1, LASTI
*
*        Determine ith pivot column and swap if necessary.
*
         PVT = ( I-1 )+ISAMAX( N-I+1, RWORK( I ), 1 )
         IF( PVT.NE.I ) THEN
            CALL CSWAP( M, A( 1, PVT ), 1, A( 1, I ), 1 )
            ITEMP = JPVT( PVT )
            JPVT( PVT ) = JPVT( I )
            JPVT( I ) = ITEMP
            RWORK( PVT ) = RWORK( I )
            RWORK( N+PVT ) = RWORK( N+I )
         END IF
*
*        Generate elementary reflector H(i).
*
         IF( I.LT.M ) THEN
            CALL CLARFG( M-I+1, A( I, I ), A( I+1, I ), 1,
     $                   TAU( I ) )
         ELSE
            CALL CLARFG( 1, A( M, M ), A( M, M ), 1, TAU( M ) )
         END IF
*
*        Apply elementary reflector H(I) to the corresponding blocks
*        of matrices A and C.
*
         AII = A( I, I )
         A( I, I ) = CONE
         IF( I.LT.N ) THEN
*
*           Apply H(I) to A(I:M,I+1:N) from the left.
*
            CALL CLARF( 'Left', M-I+1, N-I, A( I, I ), 1,
     $                  CONJG( TAU( I ) ), A( I, I+1 ), LDA, WORK )
         END IF
         IF( ( JOB.EQ.2 ).AND.( K.GT.0 ) ) THEN
*
*           Apply H(I) to C(I:M,1:K) from the left.
*
            CALL CLARF( 'Left', M-I+1, K, A( I, I ), 1,
     $                   CONJG( TAU( I ) ), C( I, 1 ), LDC, WORK )
         ELSE IF( ( JOB.EQ.3 ).AND.( K.GT.0 ) ) THEN
*
*           Apply the transpose of H(I) to C(1:K,I:M) from the right.
*
            CALL CLARF( 'Right', K, M-I+1, A( I, I ), 1,
     $                   TAU( I ), C( 1, I ), LDC, WORK )
         END IF
         A( I, I ) = AII
*
*        Update partial column norms.
*
         IF( I.LT.LASTI ) THEN
            DO 30 J = I+1, N
               IF( RWORK( J ).NE.ZERO ) THEN
                  TEMP = ONE-( ABS( A( I, J ) )/RWORK( J ) )**2
                  TEMP = MAX( TEMP, ZERO )
                  TEMP2 = ONE+0.05*TEMP*( RWORK( J )/RWORK( N+J ) )**2
                  IF( TEMP2.EQ.ONE ) THEN
                     RWORK( J ) = SCNRM2( M-I, A( I+1, J ), 1 )
                     RWORK( N+J ) = RWORK( J )
                  ELSE
                     RWORK( J ) = RWORK( J )*SQRT( TEMP )
                  END IF
               END IF
 30         CONTINUE
         END IF
*
*        Check new column for independence.
*
         IF( I.EQ.1 ) THEN
            MXNM = ABS( A( 1, 1 ) )
            SMIN = MXNM
            SMAX = MXNM
            X( 1 ) = CONE
            IF( MXNM.GT.ZERO ) THEN
               LACPTD = 1
            ELSE
               SVLUES( IAFTER ) = SMIN
               GOTO 50
            END IF
         ELSE
            SMAXPR = SLASMX( I )*MXNM
            IF( CLAUC1( I, X, SMIN, A( 1, I ), A( I, I ),
     $          SMAXPR*IRCOND ) ) THEN
*
*              Column accepted.
*
               SMAX = SMAXPR
               LACPTD = LACPTD + 1
            ELSE
*
*              Column rejected.
*
               GOTO 50
            END IF
         END IF
 20   CONTINUE
 50   SVLUES( IMAX ) = SMAX
      SVLUES( IBEFOR ) = SMIN
      IF( LACPTD.EQ.DSRD ) THEN
*
*        DSRD independent columns have been found.
*
         SVLUES( IAFTER ) = SMIN
         SVLUES( IMIN ) = SMIN
      ELSE
*
*        All remaining columns rejected.
*
         I = OFFSET + LACPTD + 1
         IF( I.LT.MN ) THEN
*
*           Factor remaining columns.
*
            CALL CGEQRF( M-I, N-I, A( I+1, I+1 ), LDA, TAU( I+1 ),
     $                   WORK, LWORK, INFO )
*
*           Apply the transformations computed in CGEQRF to the
*           corresponding part of matrix C.
*
            IF( ( JOB.EQ.2 ).AND.( K.GT.0 ) ) THEN
*
*              Apply them to C(I+1:M,1:K) from the left.
*
               CALL CUNMQR( 'Left', 'Conjugate Transpose',
     $                     M-I, K, MN-I,
     $                     A( I+1, I+1 ), LDA, TAU( I+1 ),
     $                     C( I+1, 1 ), LDC, WORK, LWORK, INFO )
            ELSE IF( ( JOB.EQ.3 ).AND.( K.GT.0 ) ) THEN
*
*              Apply the transpose of them to C(1:K,I+1:M) from the
*              right.
*
               CALL CUNMQR( 'Right', 'No Transpose', K, M-I, MN-I,
     $                     A( I+1, I+1 ), LDA, TAU( I+1 ),
     $                     C( 1, I+1 ), LDC, WORK, LWORK, INFO )
            END IF
         END IF
*
*        Use incremental condition estimation to get an estimate
*        of the smallest singular value.
*
         DO 60 I = MAX( 2, OFFSET+LACPTD+1 ), MN
            CALL CLAIC1( 2, I-1, X, SMIN, A( 1, I ), A( I, I ),
     $                  SMINPR, SINE, COSINE )
            CALL CSCAL( I-1, SINE, X, 1 )
            X( I ) = COSINE
            SMIN = SMINPR
            IF( I.EQ.OFFSET+LACPTD+1 ) THEN
                SVLUES( IAFTER ) = SMIN
            END IF
 60      CONTINUE
         SVLUES( IMIN ) = SMIN
      END IF
      RETURN
*
*     End of CGEQPC
*
      END
