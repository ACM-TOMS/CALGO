      SUBROUTINE DGELSZ( M, N, NRHS, A, LDA, B, LDB, JPVT, RCOND, RANK,
     $                   WORK, LWORK, INFO )
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
*     instrumented version of DGELSZ to count ops
*
*     .. Scalar Arguments ..
      INTEGER            INFO, LDA, LDB, LWORK, M, N, NRHS, RANK
      DOUBLE PRECISION   RCOND
*     ..
*     .. Array Arguments ..
      INTEGER            JPVT( * )
      DOUBLE PRECISION   A( LDA, * ), B( LDB, * ), WORK( * )
*     ..
*     Common block to return operation counts and timings
*     .. Common blocks ..
      COMMON             / LATIME / OPS, ITCNT
      COMMON             / LSTIME / OPCNT, TIMNG
*     ..
*     .. Scalars in Common ..
      DOUBLE PRECISION   ITCNT, OPS
*     ..
*     .. Arrays in Common ..
      DOUBLE PRECISION   OPCNT( 6 ), TIMNG( 6 )
*     ..
*
*  Purpose
*  =======
*
*  DGELSZ computes the minimum-norm solution to a real linear least
*  squares problem:
*      minimize || A * X - B ||
*  using a complete orthogonal factorization of A.  A is an M-by-N
*  matrix which may be rank-deficient.
*
*  Several right hand side vectors b and solution vectors x can be
*  handled in a single call; they are stored as the columns of the
*  M-by-NRHS right hand side matrix B and the N-by-NRHS solution
*  matrix X.
*
*  The routine first computes enough of a QR factorization with column
*  pivoting so that P, Q, R11 and R12 can be determined in:
*      A * P = Q * [ R11 R12 ]
*                  [  0  R22 ]
*  with R11 defined as the largest leading submatrix whose estimated
*  condition number is less than 1/RCOND.  The order of R11, RANK,
*  is the effective rank of A.
*
*  Then, R22 is considered to be negligible, and R12 is annihilated
*  by orthogonal transformations from the right, arriving at the
*  complete orthogonal factorization:
*     A * P = Q * [ T11 0 ] * Z
*                 [  0  0 ]
*  The minimum-norm solution is then
*     X = P * Z' [ inv(T11)*Q1'*B ]
*                [        0       ]
*  where Q1 consists of the first RANK columns of Q.
*
*  This routine is basically identical to the original xGELSX except
*  three differences:
*    o The call to the subroutine xGEQPF has been substituted by the
*      the call to the subroutine xLAQP3. The subroutine xLAQP3 uses
*      Level 3 BLAS to compute a partial QR factorization of A. Only
*      enough of the factorization is completed so that P, Q, R11
*      and R12 can be determined.  For low-rank problems the
*      calculation of X is faster than if xGELSZ calls xGEQPF or
*      xGEQP3, which factor A completely.
*    o Matrix B (the right hand side) is updated with Blas-3.
*    o The permutation of matrix B (the right hand side) is faster and
*      more simple.
*
*  Arguments
*  =========
*
*  M       (input) INTEGER
*          The number of rows of the matrix A.  M >= 0.
*
*  N       (input) INTEGER
*          The number of columns of the matrix A.  N >= 0.
*
*  NRHS    (input) INTEGER
*          The number of right hand sides, i.e., the number of
*          columns of matrices B and X. NRHS >= 0.
*
*  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
*          On entry, the M-by-N matrix A.
*          On exit, A has been overwritten by details of its
*          complete orthogonal factorization.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,M).
*
*  B       (input/output) DOUBLE PRECISION array, dimension (LDB,NRHS)
*          On entry, the M-by-NRHS right hand side matrix B.
*          On exit, the N-by-NRHS solution matrix X.
*
*  LDB     (input) INTEGER
*          The leading dimension of the array B. LDB >= max(1,M,N).
*
*  JPVT    (input/output) INTEGER array, dimension (N)
*          On entry, if JPVT(i) .ne. 0, the i-th column of A is permuted
*          to the front of AP, otherwise column i is a free column.
*          The determination of the effective rank, RANK, is more
*          reliable if JPVT and RCOND are chosen so that the estimated
*          condition number of the fixed columns is less than 1 / RCOND.
*          On exit, if JPVT(i) = k, then the i-th column of AP
*          was the k-th column of A.
*
*  RCOND   (input) DOUBLE PRECISION
*          RCOND is used to determine the effective rank of A, which
*          is defined as the order of the largest leading triangular
*          submatrix R11 in the QR factorization with pivoting of A,
*          whose estimated condition number < 1/RCOND.
*
*  RANK    (output) INTEGER
*          The effective rank of A, i.e., the order of the submatrix
*          R11.  This is the same as the order of the submatrix T11
*          in the complete orthogonal factorization of A.
*
*  WORK    (workspace/output) DOUBLE PRECISION array, dimension (LWORK)
*          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
*
*  LWORK   (input) INTEGER
*          The dimension of the array WORK.
*          The unblocked strategy requires that:
*             LWORK >= MAX( 3*MN+3*N+1, 2*MN+NRHS ),
*          where MN = min( M, N ).
*          The block algorithm requires that:
*             LWORK >= MAX( 3*MN+2*N+NB*(N+1), 2*MN+NB*NRHS ),
*          where NB is an upper bound on the blocksize returned
*          by ILAENV for the routines DGEQRF, DGERQF, DORMQR,
*          and DORMRQ.
*
*          If LWORK = -1, then a workspace query is assumed; the routine
*          only calculates the optimal size of the WORK array, returns
*          this value as the first entry of the WORK array, and no error
*          message related to LWORK is issued by XERBLA.
*
*  INFO    (output) INTEGER
*          = 0: successful exit
*          < 0: If INFO = -i, the i-th argument had an illegal value.
*
*  Further Details
*  ===============
*
*  This is a modification of LAPACK routine xGELSY written by
*    A. Petitet, Computer Science Dept., Univ. of Tenn., Knoxville, USA
*    E. Quintana-Orti, Depto. de Informatica, Universidad Jaime I, Spain
*    G. Quintana-Orti, Depto. de Informatica, Universidad Jaime I, Spain
*  Modified by L. Foster, Department of Mathematics, San Jose State
*    University, San Jose, CA  and R. Kommu, Department of Physics, San
*    Jose State University, San Jose, CA.
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            LQUERY
      INTEGER            GELSZ, I, IASCL, IBSCL, J, LAQP3, LWKOPT, MN,
     $                   NB, NB1, NB2, NB3, NB4, ORMQR, ORMRZ, TRSM,
     $                   TZRZF
      DOUBLE PRECISION   ADDS, ANRM, BIGNUM, BNRM, EM, EN, ER, MULTS,
     $                   SMLNUM, T1, T2, WSIZE
*     ..
*     .. External Functions ..
      INTEGER            ILAENV
      DOUBLE PRECISION   DLAMCH, DLANGE, DOPBL3, DOPLA, DSECND
      EXTERNAL           ILAENV, DLAMCH, DLANGE, DOPBL3, DOPLA, DSECND
*     ..
*     .. External Subroutines ..
      EXTERNAL           DCOPY, DLABAD, DLAQP3, DLASCL, DLASET, DORMQR,
     $                   DORMRZ, DTRSM, DTZRZF, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, DBLE, MAX, MIN
*     ..
*     .. Data statements ..
      DATA               GELSZ / 1 / , LAQP3 / 2 / , ORMQR / 4 / ,
     $                   ORMRZ / 6 / , TRSM / 5 / , TZRZF / 3 /
*     ..
*     .. Executable Statements ..
      MN = MIN( M, N )
*
*     Test the input arguments.
*
      INFO = 0
      NB1 = ILAENV( 1, 'DGEQRF', ' ', M, N, -1, -1 )
      NB2 = ILAENV( 1, 'DGERQF', ' ', M, N, -1, -1 )
      NB3 = ILAENV( 1, 'DORMQR', ' ', M, N, NRHS, -1 )
      NB4 = ILAENV( 1, 'DORMRQ', ' ', M, N, NRHS, -1 )
      NB = MAX( NB1, NB2, NB3, NB4 )
      LWKOPT = MAX( 1, 3*MN+2*N+NB*( N+1 ), 2*MN+NB*NRHS )
      WORK( 1 ) = DBLE( LWKOPT )
      LQUERY = ( LWORK.EQ.-1 )
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( NRHS.LT.0 ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -5
      ELSE IF( LDB.LT.MAX( 1, M, N ) ) THEN
         INFO = -7
      ELSE IF( LWORK.LT.MAX( 1, 3*MN+3*N+1, 2*MN+NRHS ) .AND. .NOT.
     $         LQUERY ) THEN
         INFO = -12
      END IF
*
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DGELSZ', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( MIN( M, N, NRHS ).EQ.0 ) THEN
         RANK = 0
         RETURN
      END IF
*
*     Get machine parameters
*
      OPCNT( GELSZ ) = OPCNT( GELSZ ) + DBLE( 2 )
      SMLNUM = DLAMCH( 'S' ) / DLAMCH( 'P' )
      BIGNUM = ONE / SMLNUM
      CALL DLABAD( SMLNUM, BIGNUM )
*
*     Scale A, B if max entries outside range [SMLNUM,BIGNUM]
*
      ANRM = DLANGE( 'M', M, N, A, LDA, WORK )
      IASCL = 0
      IF( ANRM.GT.ZERO .AND. ANRM.LT.SMLNUM ) THEN
*
*        Scale matrix norm up to SMLNUM
*
         OPCNT( GELSZ ) = OPCNT( GELSZ ) + DBLE( M*N )
         CALL DLASCL( 'G', 0, 0, ANRM, SMLNUM, M, N, A, LDA, INFO )
         IASCL = 1
      ELSE IF( ANRM.GT.BIGNUM ) THEN
*
*        Scale matrix norm down to BIGNUM
*
         OPCNT( GELSZ ) = OPCNT( GELSZ ) + DBLE( M*N )
         CALL DLASCL( 'G', 0, 0, ANRM, BIGNUM, M, N, A, LDA, INFO )
         IASCL = 2
      ELSE IF( ANRM.EQ.ZERO ) THEN
*
*        Matrix all zero. Return zero solution.
*
         CALL DLASET( 'F', MAX( M, N ), NRHS, ZERO, ZERO, B, LDB )
         RANK = 0
         GO TO 50
      END IF
*
      BNRM = DLANGE( 'M', M, NRHS, B, LDB, WORK )
      IBSCL = 0
      IF( BNRM.GT.ZERO .AND. BNRM.LT.SMLNUM ) THEN
*
*        Scale matrix norm up to SMLNUM
*
         OPCNT( GELSZ ) = OPCNT( GELSZ ) + DBLE( M*NRHS )
         CALL DLASCL( 'G', 0, 0, BNRM, SMLNUM, M, NRHS, B, LDB, INFO )
         IBSCL = 1
      ELSE IF( BNRM.GT.BIGNUM ) THEN
*
*        Scale matrix norm down to BIGNUM
*
         OPCNT( GELSZ ) = OPCNT( GELSZ ) + DBLE( M*NRHS )
         CALL DLASCL( 'G', 0, 0, BNRM, BIGNUM, M, NRHS, B, LDB, INFO )
         IBSCL = 2
      END IF
*
*     Compute enough of a QR factorization with column pivoting
*     so that P, Q, R11 and R12 can be determined in:
*         A * P = Q * [ R11 R12 ]
*                     [  0  R22 ]
*     with R11 defined as the largest leading submatrix whose estimated
*     condition number is less than 1/RCOND.  The order of R11, RANK,
*     is the effective rank of A.
*
      OPS = 0
      T1 = DSECND( )
      CALL DLAQP3( M, N, A, LDA, JPVT, RCOND, RANK, WORK( 1 ),
     $             WORK( MN+1 ), LWORK-MN, INFO )
      T2 = DSECND( )
      TIMNG( LAQP3 ) = TIMNG( LAQP3 ) + ( T2-T1 )
*
*     Calculate approximate flop count in partial factorization.
*
      EM = M
      EN = N
      ER = RANK
      MULTS = 2*EN*EM + ER*( 3*EM+5*EN+2*EM*EN-( ER+1 )*
     $        ( 4+EN+EM-( 2*ER+1 ) / 3 ) )
      ADDS = EN*EM + ER*( 2*EM+EN+2*EM*EN-( ER+1 )*
     $       ( 2+EN+EM-( 2*ER+1 ) / 3 ) )
      OPCNT( LAQP3 ) = OPCNT( LAQP3 ) + ADDS + MULTS
*
*     Calculate approximate flop count in condition estimation.
*
      OPCNT( GELSZ ) = OPCNT( GELSZ ) + DBLE( RANK*RANK-1 ) + OPS
*
      WSIZE = MN + WORK( MN+1 )
*
*     workspace: 3*MN+2*N+NB*(N+1).
*     Details of Householder rotations stored in WORK(1:RANK).
*
      IF( ABS( A( 1, 1 ) ).EQ.ZERO ) THEN
         CALL DLASET( 'F', MAX( M, N ), NRHS, ZERO, ZERO, B, LDB )
         GO TO 50
      END IF
*
*     workspace: 3*MN.
*
*     Logically partition R = [ R11 R12 ]
*                             [  0  R22 ]
*     where R11 = R(1:RANK,1:RANK)
*
*     [R11,R12] = [ T11, 0 ] * Y
*
      IF( RANK.LT.N ) THEN
         OPCNT( TZRZF ) = OPCNT( TZRZF ) +
     $                    DOPLA( 'DTZRQF', RANK, N, 0, 0, 0 )
         T1 = DSECND( )
         CALL DTZRZF( RANK, N, A, LDA, WORK( MN+1 ), WORK( 2*MN+1 ),
     $                LWORK-2*MN, INFO )
         T2 = DSECND( )
         TIMNG( TZRZF ) = TIMNG( TZRZF ) + ( T2-T1 )
      END IF
*
*     workspace: 2*MN.
*     Details of Householder rotations stored in WORK(MN+1:2*MN)
*
*     B(1:M,1:NRHS) := Q' * B(1:M,1:NRHS)
*
      OPCNT( ORMQR ) = OPCNT( ORMQR ) +
     $                 DOPLA( 'DORMQR', M, NRHS, RANK, 0, 0 )
      T1 = DSECND( )
      CALL DORMQR( 'Left', 'Transpose', M, NRHS, RANK, A, LDA,
     $             WORK( 1 ), B, LDB, WORK( 2*MN+1 ), LWORK-2*MN, INFO )
      T2 = DSECND( )
      TIMNG( ORMQR ) = TIMNG( ORMQR ) + ( T2-T1 )
      WSIZE = MAX( WSIZE, 2*MN+WORK( 2*MN+1 ) )
*
*     workspace: 2*MN+NB*NRHS.
*
*     B(1:RANK,1:NRHS) := inv(T11) * B(1:RANK,1:NRHS)
*
      OPCNT( TRSM ) = OPCNT( TRSM ) + DOPBL3( 'DTRSM ', RANK, NRHS, 0 )
      T1 = DSECND( )
      CALL DTRSM( 'Left', 'Upper', 'No transpose', 'Non-unit', RANK,
     $            NRHS, ONE, A, LDA, B, LDB )
      T2 = DSECND( )
      TIMNG( TRSM ) = TIMNG( TRSM ) + ( T2-T1 )
*
      DO 20 J = 1, NRHS
         DO 10 I = RANK + 1, N
            B( I, J ) = ZERO
   10    CONTINUE
   20 CONTINUE
*
*     B(1:N,1:NRHS) := Y' * B(1:N,1:NRHS)
*
      IF( RANK.LT.N ) THEN
         NB = ILAENV( 1, 'DORMRQ', 'LT', N, NRHS, RANK, -1 )
         OPCNT( ORMRZ ) = OPCNT( ORMRZ ) +
     $                    DOPLA( 'DORMRQ', N, NRHS, RANK, 0, NB )
         T1 = DSECND( )
         CALL DORMRZ( 'Left', 'Transpose', N, NRHS, RANK, N-RANK, A,
     $                LDA, WORK( MN+1 ), B, LDB, WORK( 2*MN+1 ),
     $                LWORK-2*MN, INFO )
         T2 = DSECND( )
         TIMNG( ORMRZ ) = TIMNG( ORMRZ ) + ( T2-T1 )
      END IF
*
*     workspace: 2*MN+NRHS.
*
*     B(1:N,1:NRHS) := P * B(1:N,1:NRHS)
*
      DO 40 J = 1, NRHS
         DO 30 I = 1, N
            WORK( JPVT( I ) ) = B( I, J )
   30    CONTINUE
         CALL DCOPY( N, WORK( 1 ), 1, B( 1, J ), 1 )
   40 CONTINUE
*
*     workspace: N.
*
*     Undo scaling
*
      IF( IASCL.EQ.1 ) THEN
         OPCNT( GELSZ ) = OPCNT( GELSZ ) + DBLE( N*NRHS+RANK*RANK )
         CALL DLASCL( 'G', 0, 0, ANRM, SMLNUM, N, NRHS, B, LDB, INFO )
         CALL DLASCL( 'U', 0, 0, SMLNUM, ANRM, RANK, RANK, A, LDA,
     $                INFO )
      ELSE IF( IASCL.EQ.2 ) THEN
         OPCNT( GELSZ ) = OPCNT( GELSZ ) + DBLE( N*NRHS+RANK*RANK )
         CALL DLASCL( 'G', 0, 0, ANRM, BIGNUM, N, NRHS, B, LDB, INFO )
         CALL DLASCL( 'U', 0, 0, BIGNUM, ANRM, RANK, RANK, A, LDA,
     $                INFO )
      END IF
      IF( IBSCL.EQ.1 ) THEN
         OPCNT( GELSZ ) = OPCNT( GELSZ ) + DBLE( N*NRHS )
         CALL DLASCL( 'G', 0, 0, SMLNUM, BNRM, N, NRHS, B, LDB, INFO )
      ELSE IF( IBSCL.EQ.2 ) THEN
         OPCNT( GELSZ ) = OPCNT( GELSZ ) + DBLE( N*NRHS )
         CALL DLASCL( 'G', 0, 0, BIGNUM, BNRM, N, NRHS, B, LDB, INFO )
      END IF
*
   50 CONTINUE
      WORK( 1 ) = DBLE( LWKOPT )
*
      RETURN
*
*     End of DGELSZ
*
      END