      SUBROUTINE DGNTST( WHICH, M, N, RTHRESH, GAP, STRIP, ISEED,
     $                   RANK, S, A, LDA, WORK )
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
*     $Date: 96/12/30 16:59:14 $
*
*     .. Scalar Arguments ..
      INTEGER            WHICH, M, N, STRIP,  RANK, LDA
      DOUBLE PRECISION   RTHRESH, GAP
*     ..
*     .. Array Arguments ..
      INTEGER            ISEED( 4 )
      DOUBLE PRECISION   A( LDA, * ), S( * ), WORK( * )
*     ..
*
*  Purpose:
*  =======
*
*  DGNTST forms a test matrix for DGEQPF, the QR factorization with
*  restricted column pivoting.
*
*  Arguments:
*  =========
*
*  WHICH   (input) INTEGER
*          Determines what kind of matrix is generated. Setting
*          MN = min(M,N), we have for WHICH =
*          1: The last MN/2 columns are generated with MN/2-1
*             singular values equal 1 and one equal to RTHRESH/GAP.
*             The remaining columns are random linear combinations of
*             those MN/2 columns, scaled by sqrt(sqrt(eps)).
*             Argument STRIP is not referenced.
*          2: columns 2:MN are generated with singular values between
*             1 and GAP*RTHRESH in arithmetic progression. Column 1 is
*             a random multiple of column2. Columns MN+1:N are random
*             linear combinations of previous columns, scaled by
*             sqrt(sqrt(eps)).
*             Argument STRIP is not referenced.
*          3: Generates an m-by-n matrix with singular values between
*             1 and GAP*RTHRESH in geometric sequence.
*             Argument STRIP is not referenced.
*          4: Generates a matrix which has STRIP
*             columns with norms in the order of sqrt(sqrt(eps))
*             up front. The rest of the columns is generated with a
*             geometric distribution of singular values between 1 and
*             GAP*RTHRESH.
*          5: Generates a matrix which has STRIP
*             columns with an arithmetic distribution of singular
*             values between 1 and GAP*RTHRESH up front. The remaining
*             columns are random linear combinations of these columns
*             with permutations of the order of sqrt(epsilon).
*          6: Matrix with random singular values between 1 and
*             RTHRESH*GAP, except for six small singular values, which
*             are all small around RTHRESH*GAP.
*          any other value or when the matrix sizes are too small for
*          a selected option to make sense:
*             generate null matrix.
*
*  M       (input) INTEGER
*          The number or rows of the matrix.
*
*  N       (input) INTEGER
*          The number of columns of the matrix A.
*
*  RTHRESH (input) DOUBLE PRECISION
*          1/RTHRESH is the acceptance threshold for the condition
*          number of a matrix.
*
*  GAP     (input) DOUBLE PRECISION
*          GAP (.gt. ONE) determines singular values around threshold.
*          The smallest singular value above RTHRESH will be
*          GAP*RTHRESH, the next singular value below RTHRESH
*          will be RTHRESH/GAP.
*
*  STRIP   (input) INTEGER
*          The width of dependent strips.
*
*  ISEED   (input/output) INTEGER array, dimension (4)
*          seed for the random number generator. ISEED(4) must be odd.
*
*  RANK    (output) INTEGER
*          The rank of the matrix generated with respect to the
*          threshold 1/RTHRESH.
*
*  S       (output) DOUBLE PRECISION array, dimension min(M,N)
*          singular values of A
*
*  A       (output) DOUBLE PRECISION array, dimension (LDA,N)
*          The m-by-n matrix being generated
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A. LDA >= M.
*
*  WORK    (workspace) DOUBLE PRECISION array,
*                      dimension M*N+3*max(M,min(M,N))+max(M,N)
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
*
*        input arguments for SLATMS:
*           BREAK:  all sing. values 1 except for last one
*           ARITH:  arithmetic sequence
*           GEOM:   geometric sequence
*
      INTEGER            BREAK, ARITH, GEOM
      PARAMETER          ( BREAK = 2, GEOM = 3, ARITH = 4 )
*     ..
*     .. Local Scalars ..
      INTEGER            MN, WIDTH1, WIDTH2, INFO, I, J
      DOUBLE PRECISION   DUMMY, RTEPS, EPS, TEMP
*     ..
*     .. External Subroutines ..
      EXTERNAL           DGEMM, DCOPY, DSCAL, DLARNV,
     $                   DGEBD2, DBDSQR
*     ..
*     .. External Functions ..
      EXTERNAL           DLARAN, SFRANK, DNRM2, DLAMCH
      INTEGER            SFRANK
      DOUBLE PRECISION   DLARAN, DLAMCH, DNRM2
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MIN, MAX, MOD, SQRT, DBLE
*     ..
*     .. Executable Statements ..
*
      MN = MIN( M, N )
      EPS = DLAMCH( 'Epsilon' )
      RTEPS = SQRT( SQRT( EPS ) )
*
      IF( WHICH.EQ.1 ) THEN
*
*        columns MN/2+1:MN of A are linearly dependent with condition
*        number GAP/RTHRESH.
*
         IF( MN.LE.1 ) GOTO 1111
         IF( MOD( MN, 2 ).EQ.0 ) THEN
            WIDTH1 = MN/2
            WIDTH2 = MN/2
         ELSE
            WIDTH1 = ( MN-1 )/2
            WIDTH2 = ( MN+1 )/2
         END IF
*
*        generate A(:,1+WIDTH1:MN) such that all singular values
*        are 1 except for last one which is RTHRESH/GAP
*
         CALL DLATMS( M, WIDTH2, 'Uniform Distribution', ISEED,
     $               'Nonsymmetric', S, -BREAK, GAP/RTHRESH, ONE, M,
     $               WIDTH2, 'No Packing', A( 1, WIDTH1+1 ), LDA,
     $               WORK, INFO )
*
*        multiply A(:,1+WIDTH1:MN) with a random
*        WIDTH2-by-WIDTH1 matrix to generate A(1,1:WIDTH1).
*
         IF( WIDTH1.GT.0 ) THEN
            CALL DLARNV(1,ISEED,WIDTH1*WIDTH2,WORK(1))
            CALL DGEMM( 'no transpose', 'no transpose', M, WIDTH1,
     $                 WIDTH2, RTEPS, A( 1, WIDTH1+1 ), LDA, WORK,
     $                 WIDTH2, ZERO, A( 1, 1 ), LDA )
         END IF
*
*        multiply A(:,1+WIDTH1:MN) with a random
*        WIDTH2-by-(N-MN) matrix to generate A(:,MN+1:N).
*
         IF( MN.LT.N ) THEN
           CALL DLARNV(1,ISEED,WIDTH2*(N-MN),WORK(1))
           CALL DGEMM( 'no transpose', 'no transpose', M, N-MN,
     $                WIDTH2, RTEPS, A( 1, WIDTH1+1 ), LDA, WORK,
     $                WIDTH2, ZERO, A( 1, MN+1 ), LDA )
         END IF
*
*        compute SVD of A
*
         CALL DLACPY( 'full', M, N, A, LDA, WORK( MN+1 ), M )
         CALL DGEBD2( M, N, WORK( MN+1 ), M, S, WORK( 1 ),
     $             WORK( MN+M*N+1 ),WORK( 2*MN+M*N+1 ),
     $             WORK( 3*MN+M*N+1 ),INFO )
         CALL DBDSQR( 'upper', MN, 0, 0, 0, S, WORK(1),
     $             DUMMY, MN, DUMMY, 1, DUMMY, MN, WORK( MN+1 ), INFO )
*
*        initialize RANK
*
         RANK = SFRANK( S, MN, RTHRESH )
         RETURN
      ELSEIF( WHICH.EQ.2 ) THEN
*
*        columns 2:MN are linearly independent. column 1 is dependent.
*        The singular values of A are similar to an arithmetic sequence
*        from 1 to GAP*RTHRESH.
*
*        generate A(:,2:MN) such that singular values decline in
*        arithmetic progression from 1 to GAP*RTHRESH.
*
         IF( MN.LT.2 ) GOTO 1111
         CALL DLATMS( M, MN-1, 'Uniform Distribution', ISEED,
     $               'Nonsymmetric', S, -ARITH, ONE/(GAP*RTHRESH), ONE,
     $               M,MN-1, 'No Packing', A( 1, 2 ), LDA, WORK, INFO )
*
*        first column is random multiple of second column
*
         CALL DCOPY( M, A( 1, 2 ), 1, A( 1, 1 ), 1 )
         CALL DSCAL( M, DLARAN( ISEED ), A( 1, 1 ), 1 )
*
*        multiply A(:,2:MN) with a random
*        (MN-1)-by-(N-MN) matrix to generate A(:,MN+1:N).
*
         IF( MN.LT.N ) THEN
           CALL DLARNV( 1, ISEED, (MN-1)*(N-MN), WORK( 1 ) )
           CALL DGEMM( 'no transpose' , 'no transpose', M, N-MN,
     $                MN-1, RTEPS, A( 1, 2 ), LDA, WORK, MN-1, ZERO,
     $                A( 1, MN+1 ), LDA )
         END IF
*
*        compute SVD of A
*
         CALL DLACPY( 'full matrix', M, N, A, LDA, WORK( MN+1 ), M )
         CALL DGEBD2( M, N, WORK( MN+1 ), M, S, WORK( 1 ),
     $             WORK( MN+M*N+1 ),WORK( 2*MN+M*N+1 ),
     $             WORK( 3*MN+M*N+1 ),INFO )
         CALL DBDSQR( 'upper',MN, 0, 0, 0, S, WORK( 1 ),
     $             DUMMY, MN, DUMMY, 1, DUMMY, MN, WORK( MN+1 ), INFO )
*
*        initialize RANK
*
         RANK = SFRANK( S, MN, RTHRESH )
         RETURN
      ELSEIF( WHICH.EQ.3 ) THEN
*
*        generate a matrix with full rank and fix the first (MN-1)
*        columns. The singular values are generated with a geometric
*        distribution.
*
         CALL DLATMS( M, N, 'Uniform', ISEED, 'Nonsymmetric', S,
     $              GEOM, ONE/(GAP*RTHRESH), ONE, M, N, 'No Packing',
     $              A, LDA, WORK, INFO )
         RANK = MN
         RETURN
      ELSEIF( WHICH.EQ.4 ) THEN
*
*        generate a matrix which has min(STRIP,N-1) small columns up front,
*        the rest of the columns is independent and generated with
*        a geometric distribution of singular values.
*
         WIDTH1 = MAX( 1, MIN( STRIP, N-1 ) )
         DO 80 J = 1, WIDTH1
            CALL DLARNV( 1, ISEED, M, A( 1, J ) )
****            CALL SSCAL(M,RTEPS,A(1,J),1)
80       CONTINUE
         IF( N.EQ.1 ) THEN
            S( 1 ) = DNRM2( M, A( 1, 1 ),1 )
            RANK = 0
            RETURN
         ELSEIF( M.EQ.1 ) THEN
            IF( N.GT.1 ) THEN
               DO 85 I = 2, N
                  A( 1, I ) = ONE
 85            CONTINUE
               S( 1 ) = DNRM2( N, A( 1, 1 ), LDA )
               RANK = 1
            ELSE
               S( 1 ) = ABS( A( 1, 1 ) )
               RANK = 0
            END IF
            RETURN
         END IF
         CALL DLATMS( M, N-WIDTH1, 'Uniform', ISEED, 'Nonsymmetric',
     $              S, GEOM, ONE/( GAP*RTHRESH ), ONE, M, N-WIDTH1,
     $              'No Packing', A( 1, WIDTH1+1 ), LDA, WORK, INFO )
*
*        compute SVD
*
         CALL DLACPY( 'full matrix', M, N, A, LDA, WORK( MN+1 ), M )
         CALL DGEBD2( M, N, WORK( MN+1 ), M, S, WORK( 1 ),
     $             WORK( MN+M*N+1 ),WORK( 2*MN+M*N+1 ),
     $             WORK( 3*MN+M*N+1 ), INFO )
         CALL DBDSQR( 'upper', MN, 0, 0, 0, S, WORK( 1 ),
     $             DUMMY, MN, DUMMY, 1, DUMMY, MN, WORK( MN+1 ), INFO )
         RANK = SFRANK( S, MN, RTHRESH )
         RETURN
      ELSEIF( WHICH.EQ.5 ) THEN
*
*        generate a matrix which has STRIP independent columns up front,
*        using an arithmetic sequence of singular values.
*        The rest of the columns is generated as linear combinations of
*        the previous ones with a perturbation of order epsilon.
*
         WIDTH1 = MIN( STRIP, MN )
         CALL DLATMS( M, WIDTH1, 'Uniform', ISEED, 'Nonsymmetric',
     $               S, ARITH, ONE/( GAP*RTHRESH ), ONE, M, WIDTH1,
     $               'No Packing', A, LDA, WORK, INFO )
        IF( N.GT.WIDTH1 ) THEN
           CALL DLARNV( 1, ISEED, WIDTH1*( N-WIDTH1 ), WORK( 1 ) )
           DO 110 J = WIDTH1+1, N
              CALL DLARNV( 1, ISEED, M, A( 1, J ) )
110        CONTINUE
           CALL DGEMM( 'no transpose', 'no transpose', M, N-WIDTH1,
     $                WIDTH1, RTEPS, A( 1, 1 ), LDA, WORK, WIDTH1, EPS,
     $                A( 1, WIDTH1+1 ), LDA )
        END IF
*
*       compute SVD of A
*
        CALL DLACPY( 'full matrix', M, N, A, LDA, WORK( MN+1 ), M )
        CALL DGEBD2( M, N, WORK( MN+1 ), M, S, WORK( 1 ),
     $            WORK( MN+M*N+1 ), WORK( 2*MN+M*N+1 ),
     $            WORK( 3*MN+M*N+1 ), INFO )
        CALL DBDSQR( 'upper', MN, 0, 0, 0, S, WORK( 1 ),
     $            DUMMY, MN, DUMMY, 1, DUMMY, MN, WORK( MN+1 ), INFO )
        RANK = SFRANK( S, MN, RTHRESH )
        RETURN
      ELSEIF( WHICH.EQ.6 ) THEN
*
*       Peter's suggestion: Matrix with random singular values,
*       between 1 and RTHRESH*GAP, and six very close singular
*       values around RTHRESH*GAP.
*
        S( 1 ) = ONE
        DO 160 I = 2, MN-6
170        TEMP = DLARAN( ISEED )
           IF( TEMP.GE.RTHRESH*GAP ) THEN
              S( I ) = TEMP
           ELSE
              GOTO 170
           END IF
160     CONTINUE
        DO 180 I = MAX( 2, MN-5 ), MN
           S( I ) = RTHRESH*GAP*( ONE+3.0*DLARAN( ISEED ) )
180     CONTINUE
        CALL DLATMS( M, N, 'Uniform', ISEED, 'Nonsymmetric', S, 0,
     $              DUMMY, ONE, M, N, 'No Packing', A, LDA, WORK, INFO )
        CALL DSORT( MN, S, 1, 'decreasing' )
        RANK = MN
        RETURN
      ELSEIF( WHICH.EQ.7 ) THEN
*
*       Generate Null matrix
*
        DO 130 J = 1, N
           DO 140 I = 1, M
              A( I, J ) = ZERO
140        CONTINUE
130     CONTINUE
        DO 150 I = 1, MN
           S( I ) = ZERO
150     CONTINUE
        RANK = 0
        RETURN
      END IF
 1111 CONTINUE
*
*     Default: Generate matrix of all ones
*
        DO 190 J = 1,N
           DO 200 I = 1, M
              A( I, J ) = ONE
 200       CONTINUE
190     CONTINUE
        S( 1 ) = SQRT( DBLE( M*N ) )
        DO 210 I = 2, MN
           S( I ) = ZERO
 210    CONTINUE
        RANK = 1
      RETURN
*
*     End of DGNTST
*
      END
