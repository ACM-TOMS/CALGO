* **********************************************************************
*
        INTEGER FUNCTION NBDFLT( ALGO, PREC, REQEST, N, B1, B2, NEEDU )
*
* ----------------------------------------------------------------------
*
* Description:
*
*   Determine an appropriate default blocking factor or crossover point
*   or intermediate bandwidth for the reduction routine specified in
*   algo.
*
* Author: Bruno Lang
*         Aachen University of Technology
*         na.blang@na-net.ornl.gov
* Date: May 17, 2000
* Version: SBR Toolbox, Rev. 1.4.1
*
* Note: >>> This routine performs no sanity checks on its arguments. <<<
*
* ----------------------------------------------------------------------
*
* Parameters:
*
        CHARACTER*5           ALGO
        CHARACTER*1           PREC, REQEST
        INTEGER               N, B1, B2
        LOGICAL               NEEDU
*
*   algo    (in) character*5
*           For which SBR routine do we need a control parameter ?
*           algo = 'SYRDB' : Reduction full -> banded.
*                = 'SBRDB' : Reduction banded -> banded.
*                = 'SBRDT' : Reduction banded -> tridiagonal.
*                = 'SYRDD' : Driver for reducing full matrices.
*                = 'SBRDD' : Driver for reducing banded matrices.
*
*   prec    (in) character*1
*           Is it the single or double precision routine ?
*           prec = 'D' : Double precision.
*                = 'S' : Single precision.
*
*   reqest  (in) character*1
*           Which control parameter is required ?
*           reqest = 'B' : Blocking factor or intermediate bandwidth.
*                  = 'C' : Crossover point to LAPACK code.
*
*   n       (in) integer
*           The size of the matrix to reduce.
*           n >= 0.
*
*   b1      (in) integer
*           The bandwidth of the matrix before the reduction.
*           0 <= b1 < n.
*
*   b2      (in) integer
*           The bandwidth of the matrix after the reduction.
*           1 <= b2 <= b1 or b1 = b2 = 0.
*
*   needu   (in) logical
*           Do we update the transformation matrix ?
*           needu = .true.  : Yes.
*                 = .false. : No.
*
* ----------------------------------------------------------------------
*
* Local variables:
*
        LOGICAL               SINGLE
        INTEGER               DELTAB
*
*   single  single or double precision routine ?
*   deltab  number of diagonals to reduce
*
* Routines called:
*
        LOGICAL               LSAME
        EXTERNAL              LSAME
*
*   lsame   case-insensitive character matching (BLAS)
*
* ----------------------------------------------------------------------
*
        SINGLE = LSAME( PREC, 'SinglePrecision' )
        DELTAB = B1 - B2
*
        IF ( ( ALGO .EQ. 'SYRDB' ) .AND.
     >       ( LSAME( REQEST, 'BlockingFactor' ) ) ) THEN
*
*            --- block size for reduction from full to banded ---
*
          IF ( NEEDU ) THEN
            IF ( SINGLE ) THEN
              NBDFLT = 32
            ELSE
              NBDFLT = 32
            ENDIF
          ELSE
            IF ( SINGLE ) THEN
              NBDFLT = 32
            ELSE
              NBDFLT = 32
            ENDIF
          ENDIF
*
        ELSEIF ( ( ALGO .EQ. 'SBRDB' ) .AND.
     >           ( LSAME( REQEST, 'BlockingFactor' ) ) ) THEN
*
*            --- block size for bandwidth reduction ---
*
          IF ( NEEDU ) THEN
            IF ( SINGLE ) THEN
              NBDFLT = ( ( DELTAB + 7 ) / 8 ) * 2
            ELSE
              NBDFLT = ( ( DELTAB + 7 ) / 8 ) * 2
            ENDIF
            IF ( NBDFLT .LT. 4 )     NBDFLT = 4
            IF ( NBDFLT .GT. DELTAB )     NBDFLT = DELTAB
          ELSE
            IF ( SINGLE ) THEN
              NBDFLT = ( ( DELTAB + 7 ) / 8 ) * 2
            ELSE
              NBDFLT = ( ( DELTAB + 7 ) / 8 ) * 2
            ENDIF
            IF ( NBDFLT .LT. 4 )     NBDFLT = 4
            IF ( NBDFLT .GT. DELTAB )     NBDFLT = DELTAB
          ENDIF
          IF ( NBDFLT .LT. 1 )     NBDFLT = 1
*
        ELSEIF ( ( ALGO .EQ. 'SBRDT' ) .AND.
     >           ( LSAME( REQEST, 'BlockingFactor' ) ) ) THEN
*
*            --- block size for reduction to tridiagonal form ---
*
          IF ( NEEDU ) THEN
            IF ( SINGLE ) THEN
              NBDFLT = ( ( B1 + 5 ) / 6 ) * 2
            ELSE
              NBDFLT = ( ( B1 + 5 ) / 6 ) * 2
            ENDIF
            IF ( NBDFLT .LT. 4 )     NBDFLT = 4
            IF ( NBDFLT .GT. DELTAB )     NBDFLT = DELTAB
          ELSE
            NBDFLT = 1
          ENDIF
*
        ELSEIF ( ( ALGO .EQ. 'SYRDD' ) .AND.
     >           ( LSAME( REQEST, 'Bandwidth' ) ) ) THEN
*
*            --- intermediate bandwidth in the reduction driver ---
*
          IF ( NEEDU ) THEN
            NBDFLT = 1
          ELSE
            IF ( ( N .GE. 300 ) .OR. ( B2 .GT. 1 ) ) THEN
              IF ( SINGLE ) THEN
                NBDFLT = 32
              ELSE
                NBDFLT = 32
              ENDIF
            ELSE
              NBDFLT = 1
            ENDIF
          ENDIF
*
        ELSEIF ( ( ALGO .EQ. 'SBRDD' ) .AND.
     >           ( LSAME( REQEST, 'Bandwidth' ) ) ) THEN
*
*            --- intermediate bandwidth in the reduction driver ---
*
          IF ( NEEDU ) THEN
            NBDFLT = 1
          ELSE
            IF ( SINGLE ) THEN
              NBDFLT = 32
            ELSE
              NBDFLT = 32
            ENDIF
          ENDIF
*
        ELSEIF ( ( ALGO .EQ. 'SBRDD' ) .AND.
     >           ( LSAME( REQEST, 'CrossoverPoint' ) ) ) THEN
*
*            --- crossover point from SBRDT to SBTRD ---
*
          IF ( NEEDU ) THEN
            IF ( SINGLE ) THEN
              NBDFLT = 8
            ELSE
              NBDFLT = 8
            ENDIF
          ELSE
            IF ( SINGLE ) THEN
              NBDFLT = 16
            ELSE
              NBDFLT = 16
            ENDIF
          ENDIF
*
        ELSE
          NBDFLT = 1
        ENDIF
*
        RETURN
        END
*
* **********************************************************************
*
        SUBROUTINE DSYRDB( UPLO, JOB, N, B, A, LDA, DRPTOL,
     >                     U, LDU, NB, TAU,
     >                     WORK, LWORK, INFO )
*
* ----------------------------------------------------------------------
*
* Description:
*
*   Given a symmetric n-times-n matrix A with either its lower or upper
*   triangle explicitly stored, dsyrdb reduces A to a symmetric banded
*   matrix B with b sub(super)diagonals by using a sequence of
*   orthogonal similarity transformations. If desired, the same
*   transformations are also applied from the right to the matrix U.
*   That is, if initially U was the identity, then the resulting banded
*   matrix B and the updated U fulfill the equation
*
*       T
*      U  * A * U   =   B  .
*
*   A must be given in full storage (with either its lower or upper
*   triangle explicitly stored) and will be overwritten with B and
*   the "mantissae" of the Householder vectors used for the reduction.
*
* Author: Bruno Lang
*         Aachen University of Technology
*         na.blang@na-net.ornl.gov
* Date: May 17, 2000
* Version: SBR Toolbox, Rev. 1.4.1
*
* Parameters:
*
        CHARACTER*1           UPLO, JOB
        INTEGER               N, B, LDA, LDU, NB, LWORK, INFO
        DOUBLE PRECISION      DRPTOL
        DOUBLE PRECISION      A( LDA, * ), U( LDU, * ), TAU( * ),
     >                        WORK( * )
*
*   uplo    (in) character*1
*           Reduce the upper or lower triangle of A ?
*           uplo = 'U' : Upper triangle.
*                = 'L' : Lower triangle.
*
*   job     (in) character*1
*           Specifies if U is required or not.
*           job = 'U' : Update U.
*               = 'N' : Do not update U.
*
*   n       (in) integer
*           The order of the matrix A.
*           n >= 0.
*
*   b       (in) integer
*           The number of sub(super-)diagonals of the reduced matrix B.
*           1 <= b < n, if n > 1.
*           b = 0     , if n <= 1.
*
*   a       (in/out) double precision array, dimension ( lda, n )
*           On entry, this array contains the symmetric matrix A with
*           either its lower (uplo = 'L') or upper (uplo = 'U') triangle
*           stored.
*           On exit, the main diagonal and the first b sub-(super-)
*           diagonals contain the banded matrix B, and the Householder
*           vectors that were used in the reduction are stored in the
*           remaining parts of the lower (upper) triangle.
*
*   lda     (in) integer
*           The leading dimension of the array a.
*           lda >= max( n, 1 ).
*
*   drptol  (in) double precision
*           Threshold for dropping the Householder transformations: if
*           the norm of the vector to be eliminated is already smaller
*           than drptol then the transform is skipped.
*           drptol >= 0.0.
*           If you do not know, use 0.0.
*
*   u       (in/out) double precision array, dimension ( ldu, n )
*           On entry, u contains some n-by-n matrix U.
*           On exit, U is postmultiplied with all orthogonal transforms
*           that were used to reduce A to banded form.
*           Accessed only if the update of U is required.
*
*   ldu     (in) integer
*           The leading dimension of the array u.
*           ldu >= max( n, 1 ).
*
*   nb      (in) integer
*           The blocking factor, as suggested by the user, for the
*           reduction and the update of U.
*           0 <= nb <= b.
*           nb = 0 : User does not know how many Householder transforms
*                    should be blocked for reasonable performance.
*                    Let the routine figure it out (depending on the
*                    workspace available).
*           nb = 1 : Do not block.
*           nb > 1 : Try to block nb Housholder transforms
*                    (if the workspace is not sufficient then a smaller
*                    blocking factor will be used).
*
*   tau     (out) double precision array, dimension ( n )
*           The scaling factors for the Householder transforms used in
*           the reduction.
*
*   work    (workspace) double precision array, dimension ( lwork )
*
*   lwork   (in) integer
*           The size of the workspace (must be provided on entry).
*           lwork >= 3 * n - 2 * b.
*
*           Workspace that small prevents blocking and may degrade
*           performance. For using nb-blocked transformations, provide
*           lwork >= nb * ( 3 * n - 2 * b ).
*           As a rule, lwork >= 50 * n should be fine.
*
*   info    (out) integer
*           On exit, info indicates the consistence of the arguments.
*           info = - 1 : uplo is none of 'U', 'L' (upper/lower case).
*                = - 2 : job is none of 'U', 'N' (upper/lower case).
*                = - 3 : n is out of range (negative).
*                = - 4 : b is out of range ( < 0, or = 0 while n > 1,
*                        or >= n ).
*                = - 6 : lda is out of range ( < n or < 1 ).
*                = - 7 : drptol is out of range (negative).
*                = - 9 : ldu is out of range ( < n or < 1 ).
*                = -10 : nb is out of range (negative).
*                = -13 : lwork is too small (see above).
*           info >=  1 : All arguments are OK.
*                        In this case, info returns the blocking factor
*                        that was used in the reduction and in the
*                        update of U.
*
* ----------------------------------------------------------------------
*
* Local variables:
*
        DOUBLE PRECISION      ZERO
        PARAMETER             ( ZERO = 0.0D0 )
*
        LOGICAL               UPPER, NEEDU
        INTEGER               NB0, NB1, NB2, J, J0,
     >                        LEN, IW, IY, IWORK
*
*   upper   reduce the upper triangle of A ?
*   needu   update U ?
*   nb0     the final blocking factor for "full" block columns
*   nb1     width of the current block column to reduce (= nb0, except
*           for the last block column)
*   nb2     number of nonzero Householder vectors in the current
*           WY transform
*   j       first column of the current block column
*   j0      last column of the current block column
*   len     length of the current block column
*   iw      points to the portion of the workspace where the block
*           transform W is stored
*   iy      points to the location of Y in the workspace
*   iwork   points to the portion of work that is used as workspace for
*           other routines
*
* The workspace is used as follows:
*
*   - The first nb1 * len <= nb0 * ( n - b ) elements hold the block W
*     of the current block WY transform.
*     Here, nb1 <= nb0 denotes the width of the current block column and
*     len <= n - b denotes its height; nb0 is the final blocking factor.
*   - The next nb1 * len elements hold the block Y.
*   - Another ( len + nb1 ) * nb1 <= n * nb0 elements are used as
*     workspace in the routines dgeqrl, dgewyg, dgewy, and dsywy.
*
* Routines called:
*
        LOGICAL               LSAME
        INTEGER               NBDFLT
        EXTERNAL              LSAME, NBDFLT
*
*   lsame   case-insensitive character matching (BLAS)
*   nbdflt  determine default blocking factor (SBR)
*
        EXTERNAL              DGEQRL, DGEWY, DGEWYG, DSYWY
*
*   dgeqrl  QR or QL factorization (SBR)
*   dgewy   apply WY transform from one side (SBR)
*   dgewyg  generate WY factors (SBRtoolbox)
*   dsywy   apply two-sided WY transform to a symmetric matrix (SBR)
*
        INTRINSIC             MAX, MIN
*
* ----------------------------------------------------------------------
*
*          --- check for errors in input ---
*
        UPPER = LSAME( UPLO, 'U' )
        NEEDU = LSAME( JOB, 'U' )
*
        IF ( .NOT. ( UPPER .OR. LSAME( UPLO, 'L' ) ) ) THEN
          INFO = - 1
        ELSEIF ( .NOT. ( NEEDU .OR. LSAME( JOB, 'N' ) ) ) THEN
          INFO = - 2
        ELSEIF ( N .LT. 0 ) THEN
          INFO = - 3
        ELSEIF ( ( ( N .LE. 1 ) .AND. ( B .NE. 0 ) ) .OR.
     >           ( ( N .GT. 1 ) .AND.
     >             ( ( B .LT. 1 ) .OR. ( B .GE. N ) ) ) ) THEN
          INFO = - 4
        ELSEIF ( LDA .LT. MAX( N, 1 ) ) THEN
          INFO = - 6
        ELSEIF ( DRPTOL .LT. ZERO ) THEN
          INFO = - 7
        ELSEIF ( LDU .LT. MAX( N, 1 ) ) THEN
          INFO = - 9
        ELSEIF ( NB .LT. 0 ) THEN
          INFO = -10
        ELSEIF ( LWORK .LT. ( 3*N-2*B ) ) THEN
          INFO = -13
        ELSE
          INFO = 1
        ENDIF
*
        IF ( INFO .NE. 1 )     GOTO 999
*
*          --- check for quick return ---
*
        IF ( ( N .LE. 1 ) .OR. ( B .EQ. ( N-1 ) ) )      GOTO 999
*
*          --- compute the maximum blocking factor for the
*                given workspace                             ---
*
        NB0 = LWORK / ( 3 * N - 2 * B )
*
*          --- adjust to the user-supplied or default value ---
*
        IF ( NB .GT. 0 ) THEN
          NB0 = MIN( NB, B, NB0 )
        ELSE
          NB0 = MIN( NBDFLT( 'SYRDB', 'DoublePrec', 'BlockingFactor',
     >                       N, N-1, B, NEEDU ),
     >               B, NB0 )
        ENDIF
*
* ......................................................................
*
        IW = 1
*
        IF ( UPPER ) THEN
*
*            --- reduce the upper triangle ---
*
          DO 100 J0 = N, B+2, - NB0
*
            LEN = J0 - B
            NB1 = MIN( NB0, LEN-1 )
            J = J0 - NB1 + 1
            IY = IW + NB1 * LEN
            IWORK = IY + NB1 * LEN
*
            CALL DGEQRL( 'L', LEN, NB1, A( 1, J ), LDA, DRPTOL,
     >                   TAU( J ), WORK( IWORK ) )
            CALL DGEWYG( 'Upper', LEN, NB1, A( 1, J ), LDA, TAU( J ),
     >                   WORK( IW ), LEN, WORK( IY ), LEN, NB2,
     >                   WORK( IWORK ) )
            CALL DGEWY( 'Left', LEN, B - NB1, NB2, A( 1, LEN+1 ), LDA,
     >                  WORK( IW ), LEN, WORK( IY ), LEN,
     >                  WORK( IWORK ) )
            CALL DSYWY( 'Upper', LEN, NB2, A( 1, 1 ), LDA,
     >                  WORK( IW ), LEN, WORK( IY ), LEN,
     >                  WORK( IWORK ) )
*
            IF ( NEEDU ) THEN
              CALL DGEWY( 'Right', N, LEN, NB2, U( 1, 1 ), LDU,
     >                    WORK( IW ), LEN, WORK( IY ), LEN,
     >                    WORK( IWORK ) )
            ENDIF
*
  100     CONTINUE
*
          TAU( B+1 ) = ZERO
*
        ELSE
*
*            --- reduce the lower triangle ---
*
          DO 200 J = 1, N-B-1, NB0
*
            LEN = N - ( J + B ) + 1
            NB1 = MIN( NB0, LEN - 1 )
            IY = IW + NB1 * LEN
            IWORK = IY + NB1 * LEN
*
            CALL DGEQRL( 'R', LEN, NB1, A( J+B, J ), LDA, DRPTOL,
     >                   TAU( J ), WORK( IWORK ) )
            CALL DGEWYG( 'Lower', LEN, NB1, A( J+B, J ), LDA, TAU( J ),
     >                   WORK( IW ), LEN, WORK( IY ), LEN, NB2,
     >                   WORK( IWORK ) )
            CALL DGEWY( 'Left', LEN, B - NB1, NB2, A( J+B, J+NB1 ), LDA,
     >                  WORK( IW ), LEN, WORK( IY ), LEN,
     >                  WORK( IWORK ) )
            CALL DSYWY( 'Lower', LEN, NB2, A( J+B, J+B ), LDA,
     >                  WORK( IW ), LEN, WORK( IY ), LEN,
     >                  WORK( IWORK ) )
*
            IF ( NEEDU ) THEN
              CALL DGEWY( 'Right', N, LEN, NB2, U( 1, J+B ), LDU,
     >                    WORK( IW ), LEN, WORK( IY ), LEN,
     >                    WORK( IWORK ) )
            ENDIF
     >
  200     CONTINUE
*
          TAU( N-B ) = ZERO
*
        ENDIF
*
*          --- return the blocking factor ---
*
        INFO = NB0
*
  999   RETURN
        END
*
* **********************************************************************
*
        SUBROUTINE DSYGTR( UPLO, N, B, A, LDA, TAU, WORK, LWORK, INFO )
*
* ----------------------------------------------------------------------
*
* Description:
*
*   dsygtr generates a real orthogonal matrix U of size n, defined as
*   the product of n - b Householder transformations, as returned by
*   dsyrdb.
*
*   If uplo = 'U' then
*
*      U  =  H( n ) * ... * H( b + 2 ),
*
*   where the "mantissa" of the length-( k - b ) Householder
*   transformation H( k ) is stored in a( 1:k-b-1, k ) and its scaling
*   factor is stored in tau( k ),
*
*   if uplo = 'L' then
*
*      U  =  H( 1 ) * ... * H( n - b - 1 ),
*
*   where the "mantiassa" of the length-(n - b - k - 1) Householder
*   transformation H( k ) is stored in a( b+k+1:n, k ) and its scaling
*   factor is stored in tau( k ).
*
*   On exit, a is overwritten with the orthogonal matrix U.
*
* Author: Bruno Lang
*         Aachen University of Technology
*         na.blang@na-net.ornl.gov
* Date: May 17, 2000
* Version: SBR Toolbox, Rev. 1.4.1
*
* Parameters:
*
        CHARACTER*1           UPLO
        INTEGER               N, B, LDA, LWORK, INFO
        DOUBLE PRECISION      A( LDA, * ), TAU( * ), WORK( * )
*
*   uplo    (in) character*1
*           Are the Householder vectors stored in the upper or lower
*           triangle of A ?
*           uplo = 'U' : Upper triangle.
*                = 'L' : Lower triangle.
*
*   n       (in) integer
*           The dimension of the reduced matrix A (and of the requested
*           orthogonal matrix).
*           n >= 0.
*
*   b       (in) integer
*           The number of sub(super-)diagonals of the reduced matrix.
*           1 <= b < n, if n > 1.
*           b = 0     , if n <= 1.
*
*   a       (in/out) double precision array, dimension ( lda, n )
*           On entry, the columns of a hold the "mantissae" of the
*           Householder transformations, as returned from dsyrdb.
*           On exit, a is overwritten with the orthogonal matrix U.
*
*   lda     (in) integer
*           The leading dimension of the array a.
*           lda >= max( n, 1 ).
*
*   tau     (in) double precision array, dimension ( n ).
*           The scaling factors for the Householder transformations, as
*           returned by dsyrdb.
*
*   work    (workspace) double precision array, dimension ( lwork )
*
*   lwork   (in) integer
*           The size of the workspace (must be provided on entry).
*           lwork >= n - b.
*
*           Workspace that small prevents blocking and may degrade
*           performance. For using nb-blocked transformations, provide
*           lwork >= nb * ( n - b ).
*           As a rule, lwork >= 20 * n should be fine.
*
*   info    (out) integer
*           On exit, info indocates the consistence of the arguments.
*           info =  1 : All arguments were OK.
*           info = -1 : uplo is none of 'U', 'L' (upper/lower case).
*                = -2 : n is out of range (negative).
*                = -3 : b is out of range ( < 0, or = 0 while n > 1,
*                       or >= n ).
*                = -5 : lda is out of range ( < n or < 1 ).
*                = -8 : lwork is out of range ( < n - b ).
*
* ----------------------------------------------------------------------
*
* Local variables:
*
        DOUBLE PRECISION      ZERO, ONE
        PARAMETER             ( ZERO = 0.0D0, ONE = 1.0D0 )
*
        LOGICAL               UPPER
        INTEGER               B1, J
*
*   upper   are the Householder vectors stored in the upper triangle ?
*   b1      = b + 1
*
* Routines called:
*
        LOGICAL               LSAME
        EXTERNAL              LSAME
*
*   lsame   case-insensitive character matching (BLAS)
*
        EXTERNAL              DCOPY, DLASET, DORGQL, DORGQR
*
*   dcopy   vector copy (BLAS)
*   dlaset  initialize matrix (LAPACK)
*   dorgql  build orthogonal matrix after QL factorization (LAPACK)
*   dorgqr  build orthogonal matrix after QR factorization (LAPACK)
*
        INTRINSIC             MAX
*
* ----------------------------------------------------------------------
*
*          --- check for errors in the input ---
*
        UPPER = LSAME( UPLO, 'U' )
        B1 = B + 1
*
        IF ( .NOT. ( UPPER .OR. ( LSAME( UPLO, 'L' ) ) ) ) THEN
          INFO = - 1
        ELSEIF ( N .LT. 0 ) THEN
          INFO = - 2
        ELSEIF ( ( ( N .LE. 1 ) .AND. ( B .NE. 0 ) ) .OR.
     >           ( ( N .GT. 1 ) .AND.
     >             ( ( B .LT. 1 ) .OR. ( B .GE. N ) ) ) ) THEN
          INFO = - 3
        ELSEIF ( LDA .LT. MAX( N, 1 ) ) THEN
          INFO = - 5
        ELSEIF ( LWORK .LT. ( N-B ) ) THEN
          INFO = - 8
        ELSE
          INFO = 1
        ENDIF
*
        IF ( INFO .NE. 1 )     GOTO 999
*
*          --- check for quick return ---
*
        IF ( B1 .GE. N ) THEN
          CALL DLASET( 'All', N, N, ZERO, ONE, A, LDA )
          GOTO 999
        ENDIF
*
* ......................................................................
*
        IF ( UPPER ) THEN
*
*            --- move the Householder vectors by b columns to
*                the left                                     ---
*
          DO 100 J = B+2, N
            CALL DCOPY( J-B1, A( 1, J ), 1, A( 1, J-B ), 1 )
  100     CONTINUE
*
*            --- set the last b columns and rows to identity ---
*
          CALL DLASET( 'All', B, N-B, ZERO, ZERO, A( N-B+1, 1 ), LDA )
          CALL DLASET( 'All', N-B, B, ZERO, ZERO, A( 1, N-B+1 ), LDA )
          CALL DLASET( 'All', B, B, ZERO, ONE, A( N-B+1, N-B+1 ), LDA )
*
*            --- generate the leading ( n-b )-by-( n-b )
*                orthogonal matrix                       ---
*
          CALL DORGQL( N-B, N-B, N-B, A, LDA, TAU( B+1 ),
     >                 WORK, LWORK, INFO )
*
        ELSE
*
*            --- move the Householder vectors by b columns to
*                the right                                    ---
*
          DO 200 J = N-1, B1, -1
            CALL DCOPY( N-J, A( J+1, J-B ), 1, A( J+1, J ), 1 )
  200     CONTINUE
*
*            --- set the first b columns and rows to identity ---
*
          CALL DLASET( 'All', N, B, ZERO, ONE, A, LDA )
          CALL DLASET( 'All', B, N-B, ZERO, ZERO, A( 1, B1 ), LDA )
*
*            --- generate the trailing ( n-b )-by-( n-b )
*                orthogonal matrix                        ---
*
          CALL DORGQR( N-B, N-B, N-B, A( B1, B1 ), LDA, TAU,
     >                 WORK, LWORK, INFO )
*
        ENDIF
*
        IF ( INFO .EQ. 0 )     INFO = 1
*
  999   RETURN
        END
*
* **********************************************************************
*
        SUBROUTINE DSBRDB( JOB, N, B1, B2, A, LDA, DRPTOL,
     >                     U, LDU, NB,
     >                     WORK, LWORK, INFO )
*
* ----------------------------------------------------------------------
*
* Description:
*
*   Given a symmetric n-times-n matrix A with b1 sub(or super-)
*   diagonals, dsbrdb reduces A to a narrower band with b2 sub(super)
*   diagonals by using a sequence of orthogonal similarity
*   transformations. If desired, the same transformations are also
*   applied from the right to the matrix U. That is, if initially U was
*   the identity, then the resulting banded matrix B and the updated U
*   fulfill the equation
*
*       T
*      U  * A * U  =  B  .
*
*   A must be given in the LAPACK packed storage scheme for lower
*   triangular banded matrices and will be overwritten with B.
*
* Author: Bruno Lang
*         Aachen University of Technology
*         na.blang@na-net.ornl.gov
* Date: May 17, 2000
* Version: SBR Toolbox, Rev. 1.4.1
*
* Parameters:
*
        CHARACTER*1           JOB
        INTEGER               N, B1, B2, LDA, LDU, NB, LWORK, INFO
        DOUBLE PRECISION      DRPTOL
        DOUBLE PRECISION      A( LDA, * ), U( LDU, * ), WORK( * )
*
*   job     (in) character*1
*           Specifies whether U is required or not.
*           job = 'U' : Update U.
*               = 'N' : Do not update U.
*
*   n       (in) integer
*           The order of the matrix.
*           n >= 0.
*
*   b1      (in) integer
*           The number of sub(super-)diagonals of A.
*           b1 >= 0 and b1 < n (if n > 0).
*
*   b2      (in) integer
*           The number of sub(super-)diagonals of B.
*           1 <= b2 <= b1 or b2 = b1 = 0.
*
*   a       (in/out) double precision array, dimension ( lda, n )
*           On entry, this array contains the symmetric matrix A. The
*           leading ( b1+1 )-by-n part of a contains the lower triangle
*           of A with A( i, j ) stored in a( 1+i-j, j ), and the
*           strictly upper triangle is not referenced.
*           On exit, for i = 1, ..., b2+1, the first n - i + 1 elements
*           in the i-th row of a contain the reduced matrix B,
*           for i = b2 + 2, ..., b1+1, the first n - i + 1 elements in
*           the i-th row of a are zeroed, and
*           for i = b1 + 2, ..., min( b1+(b1-b2)+nb, n, lda ), the first
*           n - i + 1 elements in the i-th row of a are destroyed.
*
*   lda     (in) integer
*           The leading dimension of the array a.
*           lda >= min( b1 + 1, n ).
*
*           To avoid unnecessary copying and therefore improve the
*           performance in the blocked case, we recommend
*           lda >= nb + b1 + ( b1 - b2 ).
*
*   drptol  (in) double precision
*           Threshold for dropping the Householder transformations: if
*           the norm of the vector to be eliminated is already smaller
*           than drptol then the transform is skipped.
*           drptol >= 0.0.
*           If you do not know, use 0.0.
*
*   u       (in/out) double precision array, dimension ( ldu, n )
*           On entry, u contains some n-by-n matrix U.
*           On exit, U is postmultiplied by all orthogonal transforms
*           that were used to reduce the bandwidth of A.
*           Accessed only if the update of U is required.
*
*   ldu     (in) integer
*           The leading dimension of the array u.
*           ldu >= max( n, 1 ).
*
*   nb      (in) integer
*           The blocking factor nb, as suggested by the user, for the
*           reduction of A and the update of U.
*           0 <= nb <= b2.
*           nb = 0 : User does not know how many Householder transforms
*                    should be blocked for reasonable performance. Let
*                    the routine figure it out (depending on the
*                    workspace available).
*              = 1 : Do not block
*              > 1 : Try to block nb Householder transforms (if the
*                    workspace is not sufficient then a smaller blocking
*                    factor will be used).
*
*   work    (workspace) double precision array, dimension ( lwork )
*
*   lwork   (in) integer
*           The size of the workspace (must be provided on entry).
*           lwork >= 3 * b1 + b2,    if U is not required, and
*                 >= n + 2 * b1 + 1, if U is required.
*
*           Workspace that small enforces non-blocked transforms.
*           For using nb-blocked transforms, we need
*           lwork >= ( 3 * b1 + b2 ) * nb,    if U is not required, and
*                 >= ( n + 2 * b1 + 1 ) * nb, if U is required.
*
*           If lda < b1 + ( b1 - b2 ) + nb, then ADDITIONAL workspace
*           is needed to buffer some of the intermediate fill-in during
*           the reduction. More precisely, we need
*           another ndiag * ( ndiag+1 ) / 2 elements to buffer one
*              ndiag-by-ndiag triangle of each diagonal block, where
*              ndiag = min( b1+(b1-b2)+nb, n ) - lda, and
*           another ( n / b1 ) * nsub * ( nsub+1 ) / 2  elements to
*              buffer an nsub-by-nsub triangle from each of the n / b1
*              subdiagonal blocks, where
*              nsub = min( b1+(b1-b2), n ) - lda
*           Note that roughly n * ( b1 + (b1-b2) - lda ) / 2 elements
*           are sufficient, and that no buffer space is required if
*           lda >= b1 + ( b1 - b2 ) + nb.
*
*   info    (out) integer
*           On exit, info indicates the consistence of the arguments.
*           info = - 1 : job is none of 'U', 'N' (upper/lower case).
*                = - 2 : n is out of range (negative).
*                = - 3 : b1 is out of range (negative or >= n ).
*                = - 4 : b2 is out of range ( < 0, or = 0 while b1 >= 1,
*                        or > b1 ).
*                = - 6 : lda is out of range ( < b1 + ( b1 - b2 ) + 1
*                        and < n ).
*                = - 7 : drptol is out of range ( < 0.0 ).
*                = - 9 : ldu is out of range ( < n ).
*                = -10 : nb is out of range (negative).
*                = -12 : lwork is too small (see above).
*           info >=  1 : All arguments are OK.
*                        In this case, info returns the blocking factor
*                        that was used in the reduction and in the
*                        update of U.
*
* ----------------------------------------------------------------------
*
* Local variables:
*
*            --- constants ---
*
        DOUBLE PRECISION      ZERO
        PARAMETER             ( ZERO = 0.0D0 )
*
*            --- for the reduction part ---
*
        LOGICAL               NEEDU
        INTEGER               DELTAB, NB0, LDH, PREW, K,
     >                        L, NBLK, NB1, LDWY
*
*   needu   is U required ?
*   deltab  = b1 - b2, the number of diagonals to remove
*   nb0     blocking factor that is ultimately used in the reduction
*           and in the update
*   ldh     = min( deltab + nb0, n - b2 ), the size of the blocks
*           affected by each block transform
*   prew    = b2 - nb0, the width of the block to be premultiplied
*   k       numbers the "elimination steps"
*   l       numbers the bulge chasing operations associated with each
*           elimination step
*   nblk    number of the current block column
*   nb1     number of nonzero Householder vectors in the current block
*           transform
*   ldwy    leading dimension of the arrays containing W and Y
*
*          --- workspace management ---
*
        INTEGER               NSUB, SZSUB, LSUB, NSUB1, NDIAG, LDIAG,
     >                        NDIAG1, ISUB, IDIAG, IC, ID, ID1, IX,
     >                        IW, IY, IWORK, J
*
*   nsub    the dimension of the triangles from the subdiagonal blocks
*           that must be buffered in the workspace
*           (nsub <= b1 + ( b1 - b2 ) - lda)
*   szsub   = nsub * ( nsub+1 ) / 2, the size of each such triangle
*   lsub    = ( n / b ) * szsub, the total size of all these triangles
*   nsub1   the dimension of the current triangle, plus 1
*   ndiag   the dimension of the triangle from the main diagonal block
*           that must be buffered in the workspace
*           (ndiag <= b1 + ( b1 - b2 ) + nb0 - lda )
*   ldiag   = ndiag * ( ndiag+1 ) / 2, the size of this triangle
*   ndiag1  the dimension of the current triangle, plus 1
*   isub    points to the portion of the workspace that will hold the
*           triangles from all the subdiagonal blocks
*   idiag   points to the portion of the workspace that will hold the
*           triangle from the current diagonal block
*   ic      column position of the next element to be buffered
*   id      diagonal position of the next element to be buffered
*   id1     = id + 1
*   ix      position in the workspace for this element
*   iw      points to the portion of work that will hold the block
*           W of the WY tranform
*   iy      points to the portion of work that will hold Y
*   iwork   points to the portion of work that is used as workspace
*           in the routines dqrupd, dsywy, and dgewy
*
* The workspace is used as follows:
*
*   - The first ( n / b1 ) * nsub * ( nsub+1 ) / 2 elements are used to
*     buffer nsub-by-nsub triangles of fill-in from the subdiagonal
*     blocks (see the discussion on lwork).
*     Here, nsub = min( b1+(b1-b2), n ) - lda <= b1 - b2 - 1.
*     These elements are needed only if nsub > 0.
*   - The next ndiag * ( ndiag+1 ) / 2 elements are used to buffer one
*     ndiag-by-ndiag triangle from the current diagonal block.
*     Here, ndiag = min( b1+(b1-b2)+nb0, n ) - lda <= b1 - 1.
*     These elements are needed only if ndiag > 0.
*   - The following ldh * nb0 elements hold the block W of the block WY
*     transform.
*     Here, ldh = min( b1-b2+nb0, n-b2 ) <= b1 - b2 + nb0 <= b1.
*   - The following ldh * nb0 elements hold the block Y.
*   - The remaining workspace is alternatively used by the routines
*     dqrupd (which needs ( n+1 )*nb0 elements if U is required and
*     ( ldh+1 )*nb0 elements otherwise), dsywy (( b1+1 )*nb0 elements),
*     and dgewy (b1*nb0 elements).
*
* Routines called:
*
        LOGICAL               LSAME
        INTEGER               NBDFLT
        EXTERNAL              LSAME, NBDFLT
*
*   lsame   case-insensitive character matching (BLAS)
*   nbdflt  determine default blocking factor
*
        EXTERNAL              DGEWY, DQRUPD, DSYWY
*
*   dgewy   apply WY block transform (SBR)
*   dqrupd  QR decomposition of a block of A, update of U, and premulti-
*           plication of another block of A (SBR)
*   dsywy   apply two-sided block transform to a symmetric matrix (SBR)
*
        INTRINSIC             MAX, MIN
*
* ----------------------------------------------------------------------
*
*          --- check for errors in input ---
*
        NEEDU = LSAME( JOB, 'U' )
*
        IF ( .NOT. ( NEEDU .OR. LSAME( JOB, 'N' ) ) ) THEN
          INFO = - 1
        ELSEIF ( N .LT. 0 ) THEN
          INFO = - 2
        ELSEIF ( ( B1 .LT. 0 ) .OR.
     >           ( ( N .GT. 0 ) .AND. ( B1 .GE. N ) ) ) THEN
          INFO = - 3
        ELSEIF ( ( B2 .LT. 0 ) .OR.
     >           ( ( B1 .GE. 1 ) .AND. ( B2 .EQ. 0 ) ) .OR.
     >           ( B2 .GT. B1 ) ) THEN
          INFO = - 4
        ELSEIF ( LDA .LT. MIN( B1+1, N ) ) THEN
          INFO = - 6
        ELSEIF ( DRPTOL .LT. ZERO ) THEN
          INFO = - 7
        ELSEIF ( NEEDU .AND. ( LDU .LT. MAX( N, 1 ) ) ) THEN
          INFO = - 9
        ELSEIF ( NB .LT. 0 ) THEN
          INFO = - 10
        ELSE
          INFO = 1
        ENDIF
*
        IF ( INFO .NE. 1 )     GOTO 999
*
*          --- check for quick return ---
*
        IF ( ( N .EQ. 0 ) .OR. ( B2 .EQ. B1 ) )     GOTO 999
*
*          --- determine blocking factor ---
*
        IF ( NB .GT. 0 ) THEN
          NB0 = NB
        ELSE
          NB0 = NBDFLT( 'SBRDB', 'PrecString', 'BlockingFactor',
     >                  N, B1, B2, NEEDU )
        ENDIF
        NB0 = MIN( NB0, B2 )
*
*          --- how many diagonals are eliminated, and how many must
*              be buffered ?                                        ---
*
        DELTAB = B1 - B2
        NSUB = MAX( MIN( B1+DELTAB, N ) - LDA, 0 )
        SZSUB = ( NSUB * ( NSUB + 1 ) ) / 2
        LSUB = ( N / B1 ) * SZSUB
*
*          --- reduce nb0 until workspace is sufficient ---
*
   10   CONTINUE
*
*            --- how much of the diagonal blocks must be buffered ? ---
*
          NDIAG = MAX( MIN( B1+DELTAB+NB0, N ) - LDA, 0 )
          LDIAG = ( NDIAG * ( NDIAG + 1 ) ) / 2
          IF ( NEEDU ) THEN
            IWORK = ( N + 2 * B1 + 1 ) * NB0
          ELSE
            IWORK = ( 3 * B1 + B2 ) * NB0
          ENDIF
          IF ( LWORK .LT. ( LSUB+LDIAG+IWORK ) ) THEN
            IF ( NB0 .GT. 1 ) THEN
              NB0 = NB0 - 1
              GOTO 10
            ELSE
              INFO = -12
              GOTO 999
            ENDIF
          ENDIF
*
*          --- OK, now that we have nb0, set other constants ---
*
        LDH = MIN( DELTAB+NB0, N-B2 )
        PREW = MAX( B2-NB0, 0 )
*
*          --- set up pointers to the submatrices stored in work
*              (the data layout is pointed out above)            ---
*
        ISUB = 1
        IDIAG = ISUB + LSUB
        IW = IDIAG + LDIAG
        IY = IW + LDH * NB0
        IWORK = IY + LDH * NB0
*
*          --- make sure that the "bulge buffer" is empty ---
*
        DO 30 K = B1+2, MIN( B1+NB0+DELTAB, LDA, N )
          DO 20 J = 1, N-K+1
            A( K, J ) = ZERO
   20     CONTINUE
   30   CONTINUE
*
        DO 40 K = ISUB, IDIAG+LDIAG-1
          WORK( K ) = ZERO
   40   CONTINUE
*
* ......................................................................
*
        DO 400 K = 1, N - B2 - 1, NB0
*
*            --- QR decomposition to eliminate the leading nb0
*                columns of the band                           ---
*
          CALL DQRUPD( MIN( LDH, N-K-B2+1 ), NB0, DELTAB,
     >                 A( B2+1, K ), LDA - 1, DRPTOL,
     >                 NEEDU, N, U( 1, K+B2 ), LDU,
     >                 PREW, A( PREW+1, K+NB0 ),
     >                 NB1, WORK( IY ), WORK( IW ), LDWY,
     >                 WORK( IWORK ) )
*
*              --- bulge chasing ---
*
          DO 300 L = K + B2, N - 1, B1
*
            NBLK = 1 + ( L - K - B2 ) / B1
            IF ( NB1 .GT. 0 ) THEN
*
*                --- symmetric update of diagonal block ---
*
              CALL DSYWY( 'Lower', MIN( LDH, N-L+1 ), NB1,
     >                    A( 1, L ), LDA - 1,
     >                    WORK( IW ), LDWY, WORK( IY ), LDWY,
     >                    WORK( IWORK ) )
            ENDIF
*
*              --- if lda < b1 + ( b1 - b2 ) + nb0, then before
*                  transforming the subdiagonal block, an
*                  ndiag-by-ndiag triangle of the diagonal block
*                  must be copied into the workspace ...         ---
*
            NDIAG = MIN( B1+LDH, N-L+1 ) - LDA
            NDIAG1 = NDIAG + 1
            IX = IDIAG
            DO 120 J = 1, NDIAG
              IC = L + J
              DO 110 ID = 1, NDIAG1-J
                WORK( IX ) = A( ID, IC )
                IX = IX + 1
  110         CONTINUE
  120       CONTINUE
*
*              --- ... the bottommost rows of the subdiagonal block
*                  must be set to zero                              ---
*
            NSUB = MIN( MAX( N-L+1-LDA, 0 ), B1+(B1-B2)-LDA )
            NSUB1 = NSUB + 1
            DO 140 ID = NSUB1, NDIAG
              ID1 = ID + 1
              DO 130 J = 1, ID
                A( ID1-J, L+J ) = ZERO
  130         CONTINUE
  140       CONTINUE
*
*              --- ... and an nsub-by-nsub triangle of the
*                  subdiagonal block must be restored from workspace ---
*
            IX = ISUB + ( NBLK - 1 ) * SZSUB
            DO 160 J = 1, NSUB
              IC = L + J
              DO 150 ID = 1, NSUB1-J
                A( ID, IC ) = WORK( IX )
                IX = IX + 1
  150         CONTINUE
  160       CONTINUE
*
            IF ( NB1 .GT. 0 ) THEN
*
*                  --- postmultiply ---
*
              IF ( N .GE. ( L+LDH ) ) THEN
                CALL DGEWY( 'Right', MIN( B1, N-L-LDH+1 ),
     >                      LDH, NB1, A( LDH+1, L ), LDA - 1,
     >                      WORK( IW ), LDWY, WORK( IY ), LDWY,
     >                      WORK( IWORK ) )
              ENDIF
            ENDIF
*
            IF ( N .GT. ( L+B1 ) ) THEN
*
*                --- QR decomposition to eliminate the leading
*                    nb0 columns of the bulge                  ---
*
              CALL DQRUPD( MIN( LDH, N-L-B1+1 ), NB0, DELTAB,
     >                     A( B1+1, L ), LDA - 1, DRPTOL,
     >                     NEEDU, N, U( 1, L+B1 ), LDU,
     >                     B1 - NB0, A( B1-NB0+1, L+NB0 ),
     >                     NB1, WORK( IY ), WORK( IW ), LDWY,
     >                     WORK( IWORK ) )
            ENDIF
*
*              --- if lda < b1 + ( b1 - b2 ) + nb0, then an
*                  nsub-by-nsub triangle of the subdiagonal block
*                  must be copied into the workspace ...          ---
*
            NSUB = MIN( MAX( N-L-NB0+1-LDA, 0 ), B1+(B1-B2)-LDA )
            NSUB1 = NSUB + 1
            IX = ISUB + ( NBLK - 1 ) * SZSUB
            DO 220 J = 1, NSUB
              IC = L + J + NB0
              DO 210 ID = 1, NSUB1-J
                WORK( IX ) = A( ID, IC )
                IX = IX + 1
  210         CONTINUE
  220       CONTINUE
*
*              --- ... and the ndiag-by-ndiag triangle of the
*                  diagonal block must be restored from workspace ---
*
            IX = IDIAG
            DO 240 J = 1, NDIAG
              IC = L + J
              DO 230 ID = 1, NDIAG1-J
                A( ID, IC ) = WORK( IX )
                IX = IX + 1
  230         CONTINUE
  240       CONTINUE
*
  300     CONTINUE
*
  400   CONTINUE
*
        INFO = NB0
*
  999   RETURN
        END
*
* **********************************************************************
*
        SUBROUTINE DSBRDT( JOB, N, B, A, LDA, DRPTOL,
     >                     D, E,
     >                     U, LDU, NB,
     >                     WORK, LWORK, INFO )
*
* ----------------------------------------------------------------------
*
* Description:
*
*   Given a symmetric n-by-n matrix A with b sub(or super-) diagonals,
*   dsbrdt reduces A to tridiagonal form using a sequence of orthogonal
*   similarity transformations. If desired, the same transformations are
*   also applied from the right to the matrix U. That is, if initially U
*   was the identity, then the resulting tridiagonal matrix T and the
*   updated U fulfill the equation
*
*       T
*      U  * A * U = T  .
*
*   A must be given in the LAPACK packed storage scheme for lower
*   triangular banded matrices.
*
* Author: Bruno Lang
*         Aachen University of technology
*         na.blang@na-net.ornl.gov
* Date: May 17, 2000
* Version: SBR Toolbox, Rev. 1.4.1
*
* Parameters:
*
        CHARACTER*1           JOB
        INTEGER               N, B, LDA, LDU, NB, LWORK, INFO
        DOUBLE PRECISION      DRPTOL
        DOUBLE PRECISION      A( LDA, * ), D( * ), E( * ), U( LDU, * ),
     >                        WORK( * )
*
*   job     (in) character*1
*           Specifies whether U is required or not.
*           job = 'U' : Update U.
*               = 'N' : Do not update U.
*
*   n       (in) integer
*           The order of the matrix A.
*           n >= 0.
*
*   b       (in) integer
*           The number of sub(super-)diagonals of A.
*           0 <= b < n, if n >= 1.
*           b = 0     , if n = 0.
*
*   a       (in/out) double precision array, dimension ( lda, n )
*           On entry, this array contains the symmetric matrix A. The
*           leading ( b+1 )-by-n part of a contains the lower triangle
*           of A with A( i, j ) stored in a( 1+i-j, j ), and the
*           strictly upper triangle is not referenced.
*           On exit, for i = 1, ..., min( 2*b, n, lda ), the first
*           n - i + 1 elements in the i-th row of a are destroyed.
*
*   lda     (in) integer
*           The leading dimension of the array a.
*           lda >= b+1.
*           For optimum performance, we recommend lda >= 2*b.
*
*   drptol  (in) double precision
*           The threshold for dropping Householder transformations: if
*           the norm of the vector to be eliminated is already smaller
*           than drptol then the transform is skipped.
*           drptol >= 0.0.
*           If you do not know, use 0.0.
*
*   d       (out) double precision array, dimension ( n )
*           On exit, d( 1 .. n ) contains the main diagonal of the
*           reduced tridiagonal matrix.
*
*   e       (out) double precision array, dimension ( n-1 )
*           On exit, e( 1 .. n-1 ) contains the sub(super-)diagonal of
*           the reduced tridiagonal matrix.
*
*   u       (in/out) double precision array, dimension ( ldu, n )
*           On entry, u contains some n-by-n matrix U.
*           On exit, U is postmultiplied with all orthogonal transforms
*           that were used to reduce the bandwidth of A.
*           Accessed only if the update of U is required.
*
*   ldu     (in) integer
*           The leading dimension of the array u.
*           ldu >= max( n, 1 ).
*
*   nb      (in) integer
*           The blocking factor nb, as suggested by the user, for the
*           update of U (nb is accessed only if U is required).
*           nb >= 0.
*               = 0 : User does not know how many Householder transforms
*                     should be blocked for reasonable performance. Let
*                     the routine figure it out (depending on the
*                     workspace available).
*               = 1 : Do not block the update of U.
*               > 1 : Try to block nb Householder transforms (if the
*                     workspace is not sufficient then a smaller
*                     blocking factor will be used).
*
*   work    (workspace) double precision array, dimension ( lwork )
*
*   lwork   (in) integer
*           The size of the workspace (must be provided on entry).
*           If lda >= 2*b, then
*             lwork >= 2 * b, if U is not required.
*             lwork >= b + n, if U is required and is updated with
*                             non-blocked Householder transforms.
*             lwork >= 2 * nb * ( n + b + nb - 1 ), if U is updated
*                             with nb-blocked transforms.
*           If lda < 2*b, then ADDITIONAL ( n/b + 1 )*szsub + nsave
*                         elements of workspace are needed to buffer
*                         intermediate elements. Here, nsave = 2*b-lda,
*                         and szsub = (nsave*(nsave-1))/2.
*                         Note that this is roughly n*(2*b-lda)/2
*                         elements.
*
*   info    (out) integer
*           On exit, info indicates the consistence of the arguments.
*           info = - 1 : job is none of 'U', 'N' (upper/lower case).
*                = - 2 : n is out of range (negative).
*                = - 3 : b is out of range (negative or
*                        >= n, while n > 0 ).
*                = - 5 : lda is too small ( < min( b+1, n ) ).
*                = - 6 : drptol is out of range (negative).
*                = -10 : ldu is too small ( < n or < 1 ).
*                = -11 : nb is out of range (negative).
*                = -13 : workspace is too small (see above).
*           info >=  1 : All arguments are OK.
*                        In this case, info returns the blocking factor
*                        that was used in the update of U (if job='U'),
*                        or info = 1 (if job='N').
*
* ----------------------------------------------------------------------
*
* Local variables:
*
*          --- constants ---
*
        DOUBLE PRECISION      ZERO, ONE
        PARAMETER             ( ZERO = 0.0D0, ONE = 1.0D0 )
*
*          --- for the reduction part ---
*
        INTEGER               NSAVE, NSAVE1, IDIAG, ISUB, SZSUB, ID, IC,
     >                        IX, NBLK, HEIGHT, WIDTH, K, J, I, J0
        DOUBLE PRECISION      TAU
*
*   nsave   the number of diagonals that must be saved
*   nsave1  = nsave + 1
*   idiag   points to the portion of the workspace that holds the part
*           of the current diagonal block of the band
*   isub    points to the portion of the workspace that holds the part
*           of the subdiagonal blocks of the band
*   szsub   number of elements that must be saved from each subdiagonal
*           block
*   id      diagonal position of the element to save
*   ic      column position of the element to save
*   ix      position in the workspace for saving the element
*   nblk    number of the block to process next
*   height  height of the block to process next
*   width   width of the subdiagonal block
*   k       numbers the columns to be eliminated
*   j       counter for the "progress" of the bulge chasing
*   tau     scaling factor of the current Householder transform
*
*          --- for updating U ---
*
        LOGICAL               NEEDU, BLOCK, ONFLY
        INTEGER               NB0, IV, IY, IW, ITAU, IWORK, LDWY0,
     >                        K0, IFIRST, NB1, LDWY, NB2, V
        DOUBLE PRECISION      BETA
*
*   needu   is U needed ?
*   block   is U needed, and is there enough workspace for blocking the
*           updates ?
*   onfly   is U needed, and is blocking prohibited by small workspace ?
*   nb0     blocking factor that is ultimately used in the update
*   iv      points to the portion of work that will hold the Householder
*           vectors for nb consecutive "eliminate & chase" sweeps
*           (required size for this buffer : n * nb0)
*   iy      points to the portion of work that will hold the block
*           transform Y in the WY representation (required size for this
*           buffer : ( b + nb0 - 1 ) * nb0)
*   iw      points to the portion of work that will hold the block
*           transform W in the WY representation (required size for this
*           buffer : ( b + nb0 - 1 ) * nb0)
*   itau    points to the portion of work that will hold the scaling
*           factors (required size for this buffer: nb0)
*   iwork   points to the protion of work that may be used as working
*           space by the routines dlarfx, dsyrf, and dgewy.
*           (required size for this buffer : n * nb0, if U is needed,
*           b otherwise)
*   ldwy0   = b + nb0 - 1, the leading dimension of "full" block
*           Householder reflectors W and Y, resp.
*   k0      number of the first sweep that contributed to the the
*           current block transform
*   ifirst  first row affected by the current block transform
*   nb1     number of transforms in the current block update (usually
*           = nb0, except for the one at the bottom or the last sweeps)
*   ldwy    size of the current block transform (usually = ldwy0,
*           except for the last block transform of the sweep)
*   nb2     number of nonzero transforms in the current block update
*           (usually = nb1, except for special matrices)
*   v       points to the portion of work that holds the current
*           Householder vector
*   beta    coefficient in the quadratic equation for determining the
*           maximum possible blocking factor
*
* The workspace is used as follows:
*
*   If lda < 2 * b then some of the intermediate fill-in is held in the
*   workspace. More precisely, nsave = 2 * b - lda is the size of the
*   (triangular) blocks that extend beyond the array a. Then
*   nsave-by-nsave blocks of the band must be swapped between the array
*   a and the workspace. Note that actually only (nsave-1)-by-(nsave-1)
*   blocks of size szsub = ( ( nsave - 1 ) * nsave ) / 2 must be saved
*   as the remaining nsave elements are made zero by the algorithm.
*   - The first szsub + nsave elements of the workspace are used to
*     buffer a triangular portion of the current diagonal block, and
*   - the next ( n / b ) * szsub elements hold n / b triangular
*     (nsave-1)-by-(nsave-1) portions of the subdiagonal blocks.
*   If lda >= 2 * b then the whole fill-in fits into the array a, and
*   the above-mentioned buffers are not needed.
*
*   If U is not needed (or U is needed, but its update is non-blocked)
*   then the next b places of work hold the current Householder vector
*   and the following b (or n, resp.) places are used as workspace for
*   applying the transformation.
*
*   If U is updated via blocked Householder transforms (comprising nb0
*   Householder transforms)
*                                                  T
*      (affected part of)   U  =  U  *  ( I + W * Y  )
*   then in the current implementation
*   - the next n * nb0 places hold the Householder vectors for
*     nb0 successive "eliminate & chase" sweeps down the band,
*   - the next ( b + nb0 - 1 ) * nb0 places hold the block Y,
*   - the next ( b + nb0 - 1 ) * nb0 places hold W,
*   - the next n * nb0 places are alternatively used as workspace for
*     the routines dsyrf (b elements), dlarfx (b elements), and dgewyr
*     (n * nb0 elements), and
*   - the buffer for nb0 scaling factors is overlapped with the first
*     nb0 elements of the "multiuse" workspace.
*
* Routines called:
*
        LOGICAL               LSAME
        INTEGER               NBDFLT
        EXTERNAL              LSAME, NBDFLT
*
*   lsame   case-insensitive character matching (BLAS)
*   nbdflt  determine default blocking factor
*
        EXTERNAL              DCOPY, DGBWYG, DGEWY, DLARFX, DLBRFG,
     >                        DSYRF
*
*   dcopy   vector copy (BLAS)
*   dgbwyg  generate WY block reflector (SBR)
*   dgewy   apply WY block reflector (SBR)
*   dlbrfg  determine Householder vector (SBR)
*   dsyrf   apply two-sided transform to a symmetric matrix (SBR)
*   dlarfx  apply Householder transform (LAPACK)
*
        INTRINSIC             DBLE, INT, MAX, MIN, MOD, SQRT
*
* ----------------------------------------------------------------------
*
*          --- check for errors in input ---
*
        NEEDU = LSAME( JOB, 'U' )
*
        IF ( .NOT. ( NEEDU .OR. LSAME( JOB, 'N' ) ) ) THEN
          INFO = -1
        ELSEIF ( N .LT. 0 ) THEN
          INFO = -2
        ELSEIF ( ( B .LT. 0 ) .OR. ( ( N .GT. 0 ) .AND. ( B .GE. N ) ) )
     >  THEN
          INFO = -3
        ELSEIF ( LDA .LT. ( B+1 ) ) THEN
          INFO = -5
        ELSEIF ( DRPTOL .LT. ZERO ) THEN
          INFO = -6
        ELSEIF ( LDU .LT. MAX( N, 1 ) ) THEN
          INFO = -10
        ELSEIF ( NEEDU .AND. ( NB .LT. 0 ) ) THEN
          INFO = -11
        ELSE
          INFO = 1
        ENDIF
*
        IF ( INFO .NE. 1 )     GOTO 999
*
*          --- check for quick return ---
*
        IF ( N .EQ. 0 ) THEN
          GOTO 999
        ELSEIF ( B .EQ. 0 ) THEN
          CALL DCOPY( N, A( 1, 1 ), LDA, D, 1 )
          DO 10 J = 1, N-1
            E( J ) = ZERO
   10     CONTINUE
          GOTO 999
        ELSEIF ( B .EQ. 1 ) THEN
          CALL DCOPY( N, A( 1, 1 ), LDA, D, 1 )
          CALL DCOPY( N-1, A( 2, 1 ), LDA, E, 1 )
          GOTO 999
        ENDIF
*
*          --- check if workspace is sufficient ---
*
        NSAVE = 2 * B - LDA
        IF ( NSAVE .GT. 0 ) THEN
          SZSUB = ( ( NSAVE - 1 ) * NSAVE ) / 2
          IX = SZSUB + NSAVE + ( N / B ) * SZSUB
        ELSE
          SZSUB = 0
          IX = 0
        ENDIF
        IF ( ( LWORK .LT. ( IX+2*B ) ) .OR.
     >       ( NEEDU .AND. ( LWORK .LT. ( IX+B+N ) ) ) ) THEN
          INFO = -13
          GOTO 999
        ENDIF
*
*          --- compute the final blocking factor for updating U ---
*
        IF ( NEEDU ) THEN
*
*            --- first determine the maximum possible value nb0
*                (for blocking nb0 updates, we need
*                2 * nb0 * ( n + b + nb0 - 1 ) workspace).      ---
*
          BETA = DBLE( N + B - 1 )
          NB0 = INT( ( SQRT( BETA**2 + 2.0D0*( LWORK-IX ) ) - BETA )
     >               / 2.0D0 )
*
*            --- adjust to the user-supplied or default value ---
*
          IF ( NB .GT. 0 ) THEN
            NB0 = MIN( NB, NB0, N-2 )
          ELSE
            NB0 = MIN( NBDFLT( 'SBRDT', 'DoublePrec', 'BlockingFactor',
     >                         N, B, 1, .TRUE. ),
     >                 NB0 )
          ENDIF
*
          IF ( NB0 .LT. 1 )     NB0 = 1
*
*            --- OK, now we are done ---
*
          BLOCK = ( NB0 .GT. 1 )
          ONFLY = ( .NOT. BLOCK )
          INFO = NB0
        ELSE
          NB0 = 1
          BLOCK = .FALSE.
          ONFLY = .FALSE.
        ENDIF
*
*          --- set up pointers to the submatrices stored in work
*              (the data layout is pointed out above)            ---
*
        IF ( NSAVE .GT. 0 ) THEN
          IDIAG = 1
          ISUB = IDIAG + ( SZSUB + NSAVE )
          IV = ISUB + ( N / B ) * SZSUB
        ELSE
          IDIAG = -1
          ISUB = -1
          IV = 1
        ENDIF
*
        LDWY0 = B + NB0 - 1
        IF ( BLOCK ) THEN
          IY = IV + N * NB0
          IW = IY + LDWY0 * NB0
          IWORK = IW + LDWY0 * NB0
          ITAU = IWORK
        ELSE
          IY = -1
          IW = -1
          IWORK = IV + B
          ITAU = -1
        ENDIF
        V = IV
*
*          --- make sure that the "bulge buffer" is empty ---
*
        DO 30 J = 2, N-B
          DO 20 K = B+2, MIN( 2*B, LDA, N-J+1 )
            A( K, J ) = ZERO
   20     CONTINUE
   30   CONTINUE
*
        DO 40 IX = 1, IV-1
          WORK( IX ) = ZERO
   40   CONTINUE
*
* ......................................................................
*
        K0 = 1
*
*          --- in the k-th pass of the following main loop,
*              the k-th column of the band is zeroed and the
*              resulting bulge is chased down the band       ---
*
        DO 700 K = 1, N-2
*
*            --- determine the length of the column and the
*                Householder vector                         ---
*
          HEIGHT = MIN( B, N-K )
          CALL DLBRFG( HEIGHT, A( 2, K ), A( 3, K ), 1, TAU, DRPTOL )
*
*            --- save the vector and zero column of A ---
*
          IF ( BLOCK ) THEN
            V = IV + ( K - K0 ) * ( N - K0 + 1 )
          ENDIF
          WORK( V ) = ONE
          CALL DCOPY( HEIGHT-1, A( 3, K ), 1, WORK( V+1 ), 1 )
          DO 110 I = 3, HEIGHT+1
            A( I, K ) = ZERO
  110     CONTINUE
*
*            --- save the new tridiagonal elements ---
*
          D( K ) = A( 1, K )
          E( K ) = A( 2, K )
*
          IF ( TAU .NE. ZERO ) THEN
*
*              --- apply the Householder transform to both sides
*                  of the first diagonal block                   ---
*
            CALL DSYRF( 'Lower', HEIGHT, WORK( V ), 1, TAU,
     >                  A( 1, K+1 ), LDA - 1, WORK( IWORK ) )
*
*              --- if the transformations of U are not blocked
*                  then apply them as they are generated       ---
*
            IF ( ONFLY ) THEN
              CALL DLARFX( 'Right', N, HEIGHT, WORK( V ), TAU,
     >                     U( 1, K+1 ), LDU, WORK( IWORK ) )
            ENDIF
          ENDIF
*
*            --- start the bulge chasing ---
*
          DO 500 J = K+1, N-B, B
            NBLK = 1 + ( J - K - 1 ) / B
*
*              --- determine the size of the subdiagonal block ---
*
            WIDTH = HEIGHT
            HEIGHT = MIN( B, N-J-(B-1) )
*
*              --- if lda is too small then parts of the diagonal
*                  block must be saved ...                        ---
*
            NSAVE = B + HEIGHT - LDA
            NSAVE1 = NSAVE + 1
            IX = IDIAG
            DO 220 J0 = 1, NSAVE
              IC = J + J0
              DO 210 ID = 1, NSAVE1-J0
                WORK( IX ) = A( ID, IC )
                IX = IX + 1
  210         CONTINUE
  220       CONTINUE
*
*              --- ... and the subdiagonal block must be restored ---
*
            IX = ISUB + ( NBLK - 1 ) * SZSUB
            IF ( HEIGHT .EQ. B ) THEN
              DO 230 J0 = 1, NSAVE
                A( NSAVE1-J0, J+J0 ) = ZERO
  230         CONTINUE
              DO 250 J0 = 2, NSAVE
                IC = J + J0 - 1
                DO 240 ID = 1, NSAVE1-J0
                  A( ID, IC ) = WORK( IX )
                  IX = IX + 1
  240           CONTINUE
  250         CONTINUE
            ELSE
              DO 270 J0 = 1, NSAVE
                IC = J + J0
                DO 260 ID = 1, NSAVE1-J0
                  A( ID, IC ) = WORK( IX )
                  IX = IX + 1
  260           CONTINUE
  270         CONTINUE
            ENDIF
*
            IF ( TAU .NE. ZERO ) THEN
*
*                --- apply the previous Householder transform
*                    to the subdiagonal block (from the right) ---
*
              CALL DLARFX( 'Right', HEIGHT, WIDTH, WORK( V ), TAU,
     >                     A( B+1, J ), LDA-1, WORK( IWORK ) )
            ENDIF
*
*              --- save the previous scaling factor ---
*
            WORK( V ) = TAU
*
*              --- generate a new Householder transform to
*                  zero the first column of the subdiagonal block ---
*
            CALL DLBRFG( HEIGHT, A( B+1, J ), A( B+2, J ), 1,
     >                   TAU, DRPTOL )
*
            IF ( BLOCK ) THEN
              V = IV + ( K - K0 ) * ( N - K0 + 1 ) + NBLK * B
            ENDIF
            WORK( V ) = ONE
            CALL DCOPY( HEIGHT-1, A( B+2, J ), 1, WORK( V+1 ), 1 )
            DO 310 I = B+2, B+HEIGHT
              A( I, J ) = ZERO
  310       CONTINUE
*
            IF ( TAU .NE. ZERO ) THEN
*
*                --- apply the transform to the subdiagonal block
*                    (from the left)                              ---
*
              CALL DLARFX( 'Left', HEIGHT, WIDTH-1, WORK( V ), TAU,
     >                     A( B, J+1 ), LDA-1, WORK( IWORK ) )
            ENDIF
*
*              --- if lda is too small then the subdiagonal block
*                  must be saved ...                              ---
*
            IX = ISUB + ( NBLK - 1 ) * SZSUB
            DO 420 J0 = 2, NSAVE
              IC = J + J0
              DO 410 ID = 1, NSAVE1-J0
                WORK( IX ) = A( ID, IC )
                IX = IX + 1
  410         CONTINUE
  420       CONTINUE
*
*              --- ... and the diagonal block must be restored ---
*
            IX = IDIAG
            DO 440 J0 = 1, NSAVE
              IC = J + J0
              DO 430 ID = 1, NSAVE1-J0
                A( ID, IC ) = WORK( IX )
                IX = IX + 1
  430         CONTINUE
  440       CONTINUE
            IF ( TAU .NE. ZERO ) THEN
*
*                --- apply the transform to the diagonal block
*                    (two-sided)                               ---
*
              CALL DSYRF( 'Lower', HEIGHT, WORK( V ), 1, TAU,
     >                    A( 1, J+B ), LDA-1, WORK( IWORK ) )
*
*                --- if the transformations of U are not blocked
*                    then apply them as they are generated       ---
*
              IF ( ONFLY ) THEN
                CALL DLARFX( 'Right', N, HEIGHT, WORK( V ), TAU,
     >                      U( 1, J+B ), LDU, WORK( IWORK ) )
              ENDIF
            ENDIF
*
  500     CONTINUE
*
*            --- save the last scaling factor ---
*
          WORK( V ) = TAU
*
* ......................................................................
*
          IF ( BLOCK ) THEN
*
*              --- (delayed) update of the transformation matrix U.
*                  To allow blocked operations, the update is only
*                  done in every nb0-th pass through the k loop     ---
*
            IF ( ( MOD( K, NB0 ) .EQ. 0 ) .OR. ( K .EQ. ( N-2 ) ) )
     >      THEN
*
*                --- the following loop applies block transforms
*                    reversely to their generation, e.g., starting
*                    at the bottom and working up.
*                      First, determine which columns ifirst ...
*                    are affected by the block transform and
*                    how many (nb1) transforms it contains         ---
*
              DO 600 IFIRST = N-MOD( N-K0-1, B ), K0+1, -B
*
*                  --- how many Householder transforms contribute
*                      to the current block transform ?           ---
*
                NB1 = MIN( K-K0+1, N-IFIRST+1 )
                LDWY = MIN( LDWY0, N-IFIRST+1 )
*
*                  --- determine the block Householder reflector
*                      corresponding to these nb1 vectors        ---
*
                V = IFIRST - K0 - 1 + IV
                CALL DCOPY( NB1, WORK( V ), N-K0+1, WORK( ITAU ), 1 )
                CALL DGBWYG( LDWY, NB1, B-1,
     >                       WORK( V ), N-K0, WORK( ITAU ),
     >                       WORK( IW ), LDWY, WORK( IY ), LDWY, NB2 )
*
*                  --- update U by postmultiplying ---
*
                CALL DGEWY( 'Right', N, LDWY, NB2,
     >                      U( 1, IFIRST ), LDU,
     >                      WORK( IW ), LDWY, WORK( IY ), LDWY,
     >                      WORK( IWORK ) )
*
  600         CONTINUE
*
              K0 = K0 + NB0
*
            ENDIF
          ENDIF
*
* ......................................................................
*
  700   CONTINUE
*
*          --- end of the main loop: cleanup for the very tail
*              of the band                                     ---
*
        D( N-1 ) = A( 1, N-1 )
        E( N-1 ) = A( 2, N-1 )
*
        D( N ) = A( 1, N )
*
  999   RETURN
        END
*
* **********************************************************************
*
        SUBROUTINE DSY2BC( UPLO, N, B, AFULL, LDFULL, ABAND, LDBAND,
     >                     INFO )
*
* ----------------------------------------------------------------------
*
* Description:
*
*   dsy2bc copies a symmetric banded A matrix from (upper or lower)
*   symmetric (full) storage to lower banded storage.
*
*   Note that afull and aband must refer to different memory locations,
*   i.e., A may NOT be repacked within the same array.
*
* Author: Bruno Lang
*         Aachen University of Technology
*         na.blang@na-net.ornl.gov
* Date: May 17, 2000
* Version: SBR Toolbox, Rev. 1.4.1
*
* Parameters:
*
        CHARACTER*1           UPLO
        INTEGER               N, B, LDFULL, LDBAND, INFO
        DOUBLE PRECISION      AFULL( LDFULL, * ), ABAND( LDBAND, * )
*
*   uplo    (in) character*1
*           Is the matrix A stored in the upper or lower triangle of the
*           array afull ?
*           uplo = 'U' : Upper triangle.
*                = 'L' : Lower triangle.
*
*   n       (in) integer
*           The size of the matrix A.
*           n >= 0.
*
*   b       (in) integer
*           The (semi-)bandwidth of the matrix A.
*           0 <= b < n, if n >= 1.
*           b = 0     , if n = 0.
*
*   afull   (in) double precision array, dimension ( ldfull, n )
*           The symmetric banded matrix A in upper or lower symmetric
*           storage.
*
*   ldfull  (in) integer
*           The leading dimension of the array afull.
*           ldfull >= max( n, 1 ).
*
*   aband   (out) double precision array, dimension ( ldband, n )
*           The symmetric banded matrix A in lower banded storage (upper
*           banded storage is not supported).
*
*   ldband  (in)
*           The leading dimension of the array aband.
*           ldband >= b + 1.
*
*   info    (out) integer
*           On exit, info indicates the consistency of the arguments.
*           info = -1 : uplo was none of 'U', 'L' (upper/lower case).
*                = -2 : n was out of range (negative).
*                = -3 : b was out of range (negative or >= n ).
*                = -5 : ldfull was out of range ( < n or < 1 ).
*                = -7 : ldband was out of range ( < b + 1 ).
*           info =  1 : All arguments were ok.
*
* ----------------------------------------------------------------------
*
* Local variables:
*
        LOGICAL               UPPER
        INTEGER               J, I
*
*   upper   is A given in upper symmetric storage ?
*
* Routines called:
*
        LOGICAL               LSAME
        EXTERNAL              LSAME
*
*   lsame   case-insensitive character matching (BLAS)
*
        INTRINSIC             MAX, MIN
*
* ----------------------------------------------------------------------
*
*          --- check for errors in the parameters ---
*
        UPPER = LSAME( UPLO, 'U' )
*
        IF ( .NOT. ( UPPER .OR. ( LSAME( UPLO, 'L' ) ) ) ) THEN
          INFO = - 1
        ELSEIF ( N .LT. 0 ) THEN
          INFO = - 2
        ELSEIF ( ( B .LT. 0 ) .OR.
     >           ( ( N .GT. 0 ) .AND. ( B .GE. N ) ) ) THEN
          INFO = - 3
        ELSEIF ( LDFULL .LT. MAX( N, 1 ) ) THEN
          INFO = - 5
        ELSEIF ( LDBAND .LT. ( B+1 ) ) THEN
          INFO = - 7
        ELSE
          INFO = 1
        ENDIF
*
        IF ( INFO .NE. 1 )     GOTO 999
*
*          --- check for quick return ---
*
        IF ( N .EQ. 0 )     GOTO 999
*
*          --- non-trivial case ---
*
        IF ( UPPER ) THEN
*
*            --- "upper to lower" ---
*
          DO 110 J = 1, N
            DO 100 I = MAX( J-B, 1 ), J
              ABAND( 1+J-I, I ) = AFULL( I, J )
  100       CONTINUE
  110     CONTINUE
        ELSE
*
*            --- "lower to lower" ---
*
          DO 210 J = 1, N
            DO 200 I = J, MIN( J+B, N )
              ABAND( 1+I-J, J ) = AFULL( I, J )
  200       CONTINUE
  210     CONTINUE
        ENDIF
*
  999   RETURN
        END
*
* **********************************************************************
*
        SUBROUTINE DSY2BI( UPLO, N, B, A, LDFULL, LDBAND, INFO )
*
* ----------------------------------------------------------------------
*
* Description:
*
*   dsy2bi copies a symmetric banded A matrix from (upper or lower)
*   symmetric (full) storage to lower banded storage within the same
*   array a (i.e., in-place).
*
* Author: Bruno Lang
*         Aachen University of Technology
*         na.blang@na-net.ornl.gov
* Date: May 17, 2000
* Version: SBR Toolbox, Rev. 1.4.1
*
* Parameters:
*
        CHARACTER*1           UPLO
        INTEGER               N, B, LDFULL, LDBAND, INFO
        DOUBLE PRECISION      A( * )
*
*   uplo    (in) character*1
*           Is the matrix A stored in the upper or lower triangle of the
*           array a ?
*           uplo = 'U' : Upper triangle.
*                = 'L' : Lower triangle.
*
*   n       (in) integer
*           The size of the matrix A.
*           n >= 0.
*
*   b       (in) integer
*           The (semi-)bandwidth of the matrix A.
*           0 <= b < n, if n >= 1.
*           b = 0     , if n = 0.
*
*   a       (in/out) double precision array,
*           dimension ( max( ldfull, ldband), n )
*           On entry, the symmetric banded matrix A in upper or lower
*           symmetric storage.
*           On exit, the symmetric banded matrix A in lower banded
*           storage (upper banded storage is not supported).
*           Note that this routine treats a as a one-dimensional array.
*
*   ldfull  (in) integer
*           The leading dimension of the array a (in full storage).
*           ldfull >= max( n, 1 ).
*
*   ldband  (in)
*           The leading dimension of the array a (in banded storage).
*           ldband >= b + 1.
*
*   info    (out) integer
*           On exit, info indicates the consistency of the arguments.
*           info = -1 : uplo was none of 'U', 'L' (upper/lower case).
*                = -2 : n was out of range (negative).
*                = -3 : b was out of range (negative or >= n ).
*                = -5 : ldfull was out of range ( < n or < 1 ).
*                = -6 : ldband was out of range ( < b + 1 ).
*           info =  1 : All arguments were ok.
*
* ----------------------------------------------------------------------
*
* Local variables:
*
        LOGICAL               UPPER
        INTEGER               B1, IBAND, IFULL, J, I
*
*   upper   is A given in upper symmetric storage ?
*   b1      = b + 1
*   iband   iband + i points to the next entry in a (banded storage)
*   ifull   ifull + i points to the next entry in a (full storage)
*
* Routines called:
*
        LOGICAL               LSAME
        EXTERNAL              LSAME
*
*   lsame   case-insensitive character matching (BLAS)
*
        EXTERNAL              DSB2BI
*
*   dsb2bi  in-place repacking from banded to banded storage (SBR)
*
        INTRINSIC             MAX, MIN
*
* ----------------------------------------------------------------------
*
*          --- check for errors in the parameters ---
*
        UPPER = LSAME( UPLO, 'U' )
        B1 = B + 1
*
        IF ( .NOT. ( UPPER .OR. ( LSAME( UPLO, 'L' ) ) ) ) THEN
          INFO = - 1
        ELSEIF ( N .LT. 0 ) THEN
          INFO = - 2
        ELSEIF ( ( B .LT. 0 ) .OR.
     >           ( ( N .GT. 0 ) .AND. ( B .GE. N ) ) ) THEN
          INFO = - 3
        ELSEIF ( LDFULL .LT. MAX( N, 1 ) ) THEN
          INFO = - 5
        ELSEIF ( LDBAND .LT. B1 ) THEN
          INFO = - 6
        ELSE
          INFO = 1
        ENDIF
*
        IF ( INFO .NE. 1 )     GOTO 999
*
*          --- check for quick return ---
*
        IF ( N .EQ. 0 )     GOTO 999
*
*          --- non-trivial case ---
*
        IF ( UPPER ) THEN
*
*            --- "upper to lower"; avoid overwriting by first
*                packing the band as tight as possible (i.e.,
*                leading dimension b + 1) and then "stretching"
*                to leading dimension ldband                    ---
*
          DO 110 J = 1, N
              IBAND = ( J - 1 ) * B1 + 1 - J
              DO 100 I = J, MIN( J+B, N )
                A( IBAND+I ) = A( (I-1)*LDFULL+J )
  100         CONTINUE
  110     CONTINUE
*
          IF ( LDBAND .GT. B1 ) THEN
            CALL DSB2BI( 'Lower', N, B, A, B1, LDBAND, INFO )
          ENDIF
*
        ELSE
*
*            --- "lower to lower"; avoid overwriting by
*                 copying the columns in the "safe" order ---
*
          IF ( LDBAND .LE. ( LDFULL+1 ) ) THEN
            DO 210 J = 1, N
              IBAND = ( J - 1 ) * LDBAND + 1 - J
              IFULL = ( J - 1 ) * LDFULL
              DO 200 I = J, MIN( J+B, N )
                A( IBAND+I ) = A( IFULL+I )
  200         CONTINUE
  210       CONTINUE
          ELSE
            DO 310 J = N, 1, -1
              IBAND = ( J - 1 ) * LDBAND + 1 - J
              IFULL = ( J - 1 ) * LDFULL
              DO 300 I = MIN( J+B, N ), J, -1
                A( IBAND+I ) = A( IFULL+I )
  300         CONTINUE
  310       CONTINUE
          ENDIF
*
        ENDIF
*
  999   RETURN
        END
*
* **********************************************************************
*
        SUBROUTINE DSB2BC( UPLO, N, B, ASRC, LDSRC, ADST, LDDST, INFO )
*
* ----------------------------------------------------------------------
*
* Description:
*
*   dsb2bc copies a symmetric banded A matrix from upper or lower
*   banded storage to lower banded storage. It may be used to pack the
*   matrix tightly by setting lddst accordingly.
*
*   Note that asrc and adst must refer to different memory locations,
*   i.e., A may NOT be repacked within the same array.
*
* Author: Bruno Lang
*         Aachen University of Technology
*         na.blang@na-net.ornl.gov
* Date: May 17, 2000
* Version: SBR Toolbox, Rev. 1.4.1
*
* Parameters:
*
        CHARACTER*1           UPLO
        INTEGER               N, B, LDSRC, LDDST, INFO
        DOUBLE PRECISION      ASRC( LDSRC, * ), ADST( LDDST, * )
*
*   uplo    (in) character*1
*           Is the matrix A given in upper or lower banded storage ?
*           uplo = 'U' : Upper banded storage.
*                = 'L' : Lower banded storage.
*
*   n       (in) integer
*           The size of the matrix A.
*           n >= 0.
*
*   b       (in) integer
*           The (semi-)bandwidth of the matrix A.
*           0 <= b < n, if n >= 1.
*           b = 0     , if n = 0.
*
*   asrc    (in) double precision array, dimension ( ldsrc, n )
*           The symmetric banded matrix A in upper or lower banded
*           storage.
*
*   ldsrc   (in) integer
*           The leading dimension of the array asrc.
*           ldsrc >= b + 1.
*
*   adst    (out) double precision array, dimension ( lddst, n )
*           The symmetric banded matrix A in lower banded storage.
*
*   lddst   (in) integer
*           The leading dimension of the array adst.
*           lddst >= b + 1.
*
*   info    (out) integer
*           On exit, info indicates the consistency of the arguments.
*           info = -1 : uplo was none of 'U', 'L' (upper/lower case).
*                = -2 : n was out of range (negative).
*                = -3 : b was out of range (negative or >= n ).
*                = -5 : ldsrc was out of range ( < b + 1 ).
*                = -7 : lddst was out of range ( < b + 1 ).
*           info =  1 : All arguments were ok.
*
* ----------------------------------------------------------------------
*
* Local variables:
*
        LOGICAL               UPPER
        INTEGER               B1, J, I
*
*   upper   is A given in upper banded storage ?
*   b1      = b + 1
*
* Routines called:
*
        LOGICAL               LSAME
        EXTERNAL              LSAME
*
*   lsame   case-insensitive character matching (BLAS)
*
        INTRINSIC             MIN
*
* ----------------------------------------------------------------------
*
*          --- check for errors in the parameters ---
*
        UPPER = LSAME( UPLO, 'U' )
        B1 = B + 1
*
        IF ( .NOT. ( UPPER .OR. ( LSAME( UPLO, 'L' ) ) ) ) THEN
          INFO = - 1
        ELSEIF ( N .LT. 0 ) THEN
          INFO = - 2
        ELSEIF ( ( B .LT. 0 ) .OR.
     >           ( ( N .GT. 0 ) .AND. ( B .GE. N ) ) ) THEN
          INFO = - 3
        ELSEIF ( LDSRC .LT. B1 ) THEN
          INFO = - 5
        ELSEIF ( LDDST .LT. B1 ) THEN
          INFO = - 7
        ELSE
          INFO = 1
        ENDIF
*
        IF ( INFO .NE. 1 )     GOTO 999
*
*          --- check for quick return ---
*
        IF ( N .EQ. 0 )     GOTO 999
*
*          --- non-trivial case ---
*
        IF ( UPPER ) THEN
*
*              --- "upper to lower" ---
*
          DO 110 J = 1, N
            DO 100 I = 1, MIN( B1, N-J+1 )
              ADST( I, J ) = ASRC( B+2-I, J+I-1 )
  100       CONTINUE
  110     CONTINUE
        ELSE
*
*            --- "lower to lower" ---
*
          DO 210 J = 1, N
            DO 200 I = 1, MIN( B1, N-J+1 )
              ADST( I, J ) = ASRC( I, J )
  200       CONTINUE
  210     CONTINUE
*
        ENDIF
*
  999   RETURN
        END
*
* **********************************************************************
*
        SUBROUTINE DSB2BI( UPLO, N, B, A, LDSRC, LDDST, INFO )
*
* ----------------------------------------------------------------------
*
* Description:
*
*   dsb2bi copies a symmetric banded A matrix from upper or lower banded
*   storage to lower banded storage within the same array a (i.e.,
*   in-place). It may be used to pack the matrix tightly by setting
*   lddst accordingly.
*
* Author: Bruno Lang
*         Aachen University of Technology
*         na.blang@na-net.ornl.gov
* Date: May 17, 2000
* Version: SBR Toolbox, Rev. 1.4.1
*
* Parameters:
*
        CHARACTER*1           UPLO
        INTEGER               N, B, LDSRC, LDDST, INFO
        DOUBLE PRECISION      A( * )
*
*   uplo    (in) character*1
*           Is the matrix A given in the upper or lower banded storage ?
*           uplo = 'U' : Upper banded storage.
*                = 'L' : Lower banded storage.
*
*   n       (in) integer
*           The size of the matrix A.
*           n >= 0.
*
*   b       (in) integer
*           The (semi-)bandwidth of the matrix A.
*           0 <= b < n.
*
*   a       (in/out) double precision array,
*           dimension ( max( ldsrc, lddst ), n )
*           On entry, the symmetric banded matrix A in upper or lower
*           banded storage.
*           On exit, the symmetric banded matrix A in lower banded
*           storage (upper banded storage is not supported).
*           Note that this routine treats a as a one-dimensional array.
*
*   ldsrc   (in) integer
*           The leading dimension of the array a (input storage scheme).
*           ldsrc >= b + 1.
*
*   lddst   (in) integer
*           The leading dimension of the array a (output storage
*           scheme).
*           lddst >= b + 1.
*
*   info    (out) integer
*           On exit, info indicates the consistency of the arguments.
*           info = -1 : uplo was none of 'U', 'L' (upper/lower case).
*                = -2 : n was out of range (negative).
*                = -3 : b was out of range (negative or >= n ).
*                = -5 : ldsrc was out of range ( < b + 1 ).
*                = -6 : lddst was out of range ( < b + 1 ).
*           info =  1 : All arguments were ok.
*
* ----------------------------------------------------------------------
*
* Local variables:
*
        LOGICAL               UPPER
        INTEGER               B1, LDSRC1, ISRC, IDST, J, I
*
*   upper   is A given in upper banded storage ?
*   b1      = b + 1
*   ldsrc1  leading dimension of the current source, minus 1
*   isrc    isrc + i points to the next source element in a
*   idst    idst + something points to the next destination element in a
*
* Routines called:
*
        LOGICAL               LSAME
        EXTERNAL              LSAME
*
*   lsame   case-insensitive character matching (BLAS)
*
        INTRINSIC             MAX, MIN
*
* ----------------------------------------------------------------------
*
*          --- check for errors in the parameters ---
*
        UPPER = LSAME( UPLO, 'U' )
        B1 = B + 1
*
        IF ( .NOT. ( UPPER .OR. ( LSAME( UPLO, 'L' ) ) ) ) THEN
          INFO = - 1
        ELSEIF ( N .LT. 0 ) THEN
          INFO = - 2
        ELSEIF ( ( B .LT. 0 ) .OR.
     >           ( ( N .GT. 0 ) .AND. ( B .GE. N ) ) ) THEN
          INFO = - 3
        ELSEIF ( LDSRC .LT. B1 ) THEN
          INFO = - 5
        ELSEIF ( LDDST .LT. B1 ) THEN
          INFO = - 6
        ELSE
          INFO = 1
        ENDIF
*
        IF ( INFO .NE. 1 )     GOTO 999
*
*          --- check for quick return ---
*
        IF ( N .EQ. 0 )     GOTO 999
*
*          --- non-trivial case ---
*
        IF ( UPPER ) THEN
*
*            --- "upper to lower" ---
*
          IF ( LDDST .LE. LDSRC ) THEN
*
*              --- no danger of overwriting ---
*
            LDSRC1 = LDSRC - 1
            DO 110 J = 1, N
              IDST = ( J - 1 ) * LDDST
              ISRC = ( J - 2 ) * LDSRC + B + 2
              DO 100 I = 1, MIN( B1, N-J+1 )
                A( IDST+I ) = A( ISRC+I*LDSRC1 )
  100         CONTINUE
  110       CONTINUE
          ELSE
*
*              --- first copy to upper storage, then transpose ---
*
            DO 210 J = N, 1, -1
              IDST = ( J - 1 ) * LDDST
              ISRC = ( J - 1 ) * LDSRC
              DO 200 I = B1, MAX( B1-J+1, 1 ), -1
                A( IDST+I ) = A( ISRC+I )
  200         CONTINUE
  210       CONTINUE
*
            LDSRC1 = LDDST - 1
            DO 230 J = 1, N
              IDST = ( J - 1 ) * LDDST
              ISRC = ( J - 2 ) * LDDST + B + 2
              DO 220 I = 1, MIN( B1, N-J+1 )
                A( IDST+I ) = A( ISRC+I*LDSRC1 )
  220         CONTINUE
  230       CONTINUE
          ENDIF
*
        ELSE
*
*            --- "lower to lower"; avoid overwriting by
*                 copying the columns in the "safe" order ---
*
          IF ( LDDST .LE. LDSRC )
     >    THEN
            DO 310 J = 1, N
              IDST = ( J - 1 ) * LDDST
              ISRC = ( J - 1 ) * LDSRC
              DO 300 I = 1, MIN( B1, N-J+1 )
                A( IDST+I ) = A( ISRC+I )
  300         CONTINUE
  310       CONTINUE
          ELSE
            DO 410 J = N, 1, -1
              IDST = ( J - 1 ) * LDDST
              ISRC = ( J - 1 ) * LDSRC
              DO 400 I = MIN( B1, N-J+1 ), 1, -1
                A( IDST+I ) = A( ISRC+I )
  400         CONTINUE
  410       CONTINUE
          ENDIF
*
        ENDIF
*
  999   RETURN
        END
*
* **********************************************************************
*
        SUBROUTINE DSBRDD( JOBU, UPLO, N, B1, B2, A, LDA, DRPTOL, D, E,
     >                     U, LDU, NSTEPS, B, NB,
     >                     WORK, LWORK, STEP, INFO )
*
* Description:
*
*   Driver routine for the multi-step reduction of symmetric banded
*   matrices to narrower banded form.
*   Given a symmetric n-by-n matrix A with b1 sub(or super-) diagonals,
*   dsbrdd reduces A to a symmetric banded matrix B with b2 sub(or
*   super-) diagonals using a sequence of orthogonal transformations.
*   If desired, the same transformations are also applied from the right
*   to the matrix U. That is, if initially U was the identity, then the
*   resulting banded matrix B and the updated U fulfill the equation
*
*       T
*      U  * A * U  =  B  .
*
*   A must be given in the LAPACK packed storage scheme for upper or
*   lower banded matrices and is returned in lower banded storage (upper
*   banded is not supported for output).
*
* Author: Bruno Lang
*         Aachen University of Technology
*         na.blang@na-net.ornl.gov
* Date: May 17, 2000
* Version: SBR Toolbox, Rev. 1.4.1
*
* Parameters:
*
        CHARACTER*1           JOBU, UPLO
        INTEGER               N, B1, B2, LDA, LDU, NSTEPS, LWORK, STEP,
     >                        INFO
        INTEGER               B( * ), NB( * )
        DOUBLE PRECISION      DRPTOL
        DOUBLE PRECISION      A( LDA, * ), D( * ), E( * ), U( LDU, * ),
     >                        WORK( * )
*
*   jobu    (in) character*1
*           Specifies whether U is required or not.
*           jobu = 'U' : Update U.
*                = 'N' : Do not update U.
*
*   uplo    (in) character*1
*           Is the matrix A given in the upper or lower banded storage ?
*           uplo = 'U' : Upper banded storage.
*                = 'L' : Lower banded storage.
*
*   n       (in) integer
*           The size of the matrix A.
*           n >= 0.
*
*   b1      (in) integer
*           The number of sub(super-)diagonals of A.
*           0 <= b1 < n, if n > 0.
*           b1 = 0,      if n = 0.
*
*   b2      (in) integer
*           The number of sub(super-)diagonals of B.
*           1 <= b2 <= b1 or b2 = b1 = 0.
*
*   a       (in/out) double precision array, dimension ( lda, n )
*           On entry, the leading (b1+1)-by-n portion of a contains the
*           symmetric banded matrix A in upper or lower banded storage.
*           If uplo='U' then for i = 1, ..., n, the matrix elements
*           A( i, j ), j = i, ..., min( i+b1, n ), are stored at
*           position a( b1+1+i-j, j ), whereas for uplo='L', the matrix
*           elements A( i, j ), j = max( i-b, 1 ), ..., i, are stored at
*           position a( 1+i-j, j ).
*           On exit, the leading (b2+1)-by-n portion of a contains the
*           reduced matrix B in lower (!) banded storage.
*
*   lda     (in) integer
*           The leading dimension of the array a.
*           lda >= min( b1+1, n ) and lda > 0.
*
*   drptol  (in) double precision
*           Threshold for dropping the Householder transformations: if
*           the norm of the vector to be eliminated is already smaller
*           than drptol then the transform is skipped.
*           drptol >= 0.0.
*           If you do not know, use 0.0.
*
*   d       (out) double precision array, dimension ( n )
*           If b2 = 1 then, on exit, d( 1 .. n ) contains the main
*           diagonal of the resulting tridiagonal matrix B.
*
*   e       (out) double precision array, dimension ( n-1 )
*           If b2 = 1 then, on exit, e( 1 .. n-1 ) contains the
*           subdiagonal of the resulting tridiagonal matrix B.
*
*   u       (in/out) double precision array, dimension ( ldu, n )
*           On entry, u contains some n-by-n matrix U.
*           On exit, U is postmultiplied with all orthogonal transforms
*           that were used to reduce the bandwidth of A.
*           Accessed only if the update of U is required.
*
*   ldu     (in) integer
*           The leading dimension of the array u.
*           ldu >= max( n, 1 ).
*
*   nsteps  (in) integer
*           The number of reduction steps.
*           nsteps >= 0.
*
*           If nsteps > 0 then the reduction from b1 to b2 is done in
*           exactly nsteps steps with user-specified intermediate
*           bandwidths and blocking factors (see paramaters b and nb,
*           below).
*           If nsteps = 0 then the routine tries to determine these
*           parameters on its own.
*           If you do not know or if you do not care, use 0.
*
*   b       (in) integer array, dimension ( nsteps-1 )
*           The sequence of intermediate bandwidths. The nsteps
*           bandwidth reduction steps are
*              b1 --> b( 1 ),
*              b( 1 ) --> b( 2 ),
*              ...
*              b( nsteps-1 ) --> b2.
*           This array is accessed only if nsteps > 1.
*
*   nb      (in) integer array, dimension ( nsteps )
*           The blocking factors nb, as suggested by the user, for the
*           reduction and/or the accumulation of the transformations
*           in each of the nb reduction steps.
*           Note that the automatic selection of suitable blocking
*           factors is enabled by setting nb( i ) = 0 for some steps i.
*           This array is accessed only if nsteps > 0.
*
*   work    (workspace) double precision array, dimension ( lwork )
*
*   lwork   (in) integer
*           The length of the workspace array.
*           The exact amount of workspace required depends on the
*           sequence of intermediate bandwidths (i.e., b), on the
*           leading dimension of the array a (i.e., lda), and on the
*           blocking factors for the reduction steps (i.e., nb).
*           lwork >= maximum workspace requirement for any of the
*                    reduction steps.
*           As a rule, lwork >= b1 * n should work fine.
*
*   step    (out) integer
*           On exit, step returns the number of the last reduction step
*           that was performed (or tried).
*
*   info    (out) integer
*           On exit, info indicates the status of the reduction.
*           info =   1 : All arguments were OK, and the reduction went
*                        without problems.
*                =  -1 : jobu is none of 'U', 'N' (upper/lower case).
*                =  -2 : uplo is none of 'U', 'L' (upper/lower case).
*                =  -3 : n is out of range (negative).
*                =  -4 : b1 is out of range ( < 0 or >= n ).
*                =  -5 : b2 is out of range ( < 0, or = 0 while b1 >= 1,
*                        or > b1 ).
*                =  -7 : lda is too small (see above).
*                =  -8 : drptol is out of range (negative).
*                = -12 : ldu is out of range ( < n ).
*                = -13 : nsteps is out of range (negative).
*                = -14 : the sequence of bandwidths b1, b( 1 ), b( 2 ),
*                        ..., b( nsteps-1 ), b2 is not monotonically
*                        decreasing
*                = -15 : nb contains negative entries
*                = -17 : lwork is too small (see above).
*           If info < 0 and step > 0 then the parameter for the step-th
*           reduction step was not OK.
*
* Local constants and variables:
*
        DOUBLE PRECISION      ZERO
        PARAMETER             ( ZERO = 0.0D0 )
*
        LOGICAL               NEEDU, UPPER
        INTEGER               BCURR, LDCURR, NBCURR, BNEXT, LDNEXT, NX
        CHARACTER*1           ULCURR
*
*   needu   is U needed ?
*   upper   is A given in upper or lower banded storage ?
*   bcurr   current bandwidth
*   ldcurr  current leading dimension of the band
*   nbcurr  blocking factor for the current step
*   bnext   bandwidth after the next reduction step
*   ldnext  leading dimension of the band after repacking
*   nx      cross-over point from LAPACK to SBR code, and block size
*   ulcurr  is the current band in upper or lower storage ?
*
* Routines called:
*
        LOGICAL               LSAME
        INTEGER               NBDFLT
        EXTERNAL              LSAME, NBDFLT
*
*   lsame   case-insensitive character matching (BLAS)
*   nbdflt  determine default crossover point from LAPACK to SBR code
*           or default intermediate bandwidth (SBR)
*
        EXTERNAL              DSB2BI, DSBRDB, DSBRDT, DSBTRD
*
*   dsb2bi  in-place repacking of banded matrices (SBR)
*   dsbrdb  one-step reduction banded -> banded (SBR)
*   dsbrdt  one-step reduction banded -> tridiagonal (SBR)
*   dsbtrd  one-step reduction banded -> tridiagonal (LAPACK)
*
        INTRINSIC             MAX, MIN
*
* ----------------------------------------------------------------------
*
*          --- check for errors in input ---
*
        NEEDU = LSAME( JOBU, 'U' )
        UPPER = LSAME( UPLO, 'U' )
*
        IF ( .NOT. ( NEEDU .OR. LSAME( JOBU, 'N' ) ) ) THEN
          INFO = -1
        ELSEIF ( .NOT. ( UPPER .OR. LSAME( UPLO, 'Lower' ) ) ) THEN
          INFO = -2
        ELSEIF ( N .LT. 0 ) THEN
          INFO = -3
        ELSEIF ( ( B1 .LT. 0 ) .OR.
     >           ( ( N .GT. 0 ) .AND. ( B1 .GE. N ) ) ) THEN
          INFO = -4
        ELSEIF ( ( B2 .LT. 0 ) .OR.
     >           ( ( B2 .EQ. 0 ) .AND. ( B1 .GE. 1 ) ) .OR.
     >           ( B2 .GT. B1 ) ) THEN
          INFO = -5
        ELSEIF ( LDA .LT. ( B1+1 ) ) THEN
          INFO = -7
        ELSEIF ( DRPTOL .LT. ZERO ) THEN
          INFO = -8
        ELSEIF ( NEEDU .AND. ( LDU .LT. MAX( N, 1 ) ) ) THEN
          INFO = -12
        ELSEIF ( NSTEPS .LT. 0 ) THEN
          INFO = -13
        ELSE
          INFO = 0
        ENDIF
        IF ( INFO .LT. 0 )     GOTO 999
*
        IF ( NSTEPS .GT. 1 ) THEN
          IF ( ( B( 1 ) .GT. B1 ) .OR. ( B2 .GT. B( NSTEPS-1 ) ) ) THEN
            INFO = -14
          ELSE
            DO 10 STEP = 1, NSTEPS-2
              IF ( B( STEP+1 ) .GT. B( STEP ) )     INFO = -14
   10       CONTINUE
          ENDIF
        ENDIF
        IF ( INFO .LT. 0 )     GOTO 999
*
        IF ( NSTEPS .GT. 0 ) THEN
          DO 20 STEP = 1, NSTEPS
            IF ( NB( STEP ) .LT. 0 )     INFO = -15
   20     CONTINUE
        ENDIF
        IF ( INFO .LT. 0 )     GOTO 999
*
*          --- set current bandwidth, leading dimension, etc. ---
*
        ULCURR = UPLO
        BCURR = B1
        LDCURR = LDA
        STEP = 0
*
*          --- while bandwidth is too large do next reduction step ---
*
  100   IF ( BCURR .GT. B2 ) THEN
          STEP = STEP + 1
*
*            --- determine target bandwidth and blocking factor
*                for the next reduction step                    ---
*
          IF ( STEP .LT. NSTEPS ) THEN
            BNEXT = B( STEP )
            NBCURR = NB( STEP )
          ELSEIF ( STEP .EQ. NSTEPS ) THEN
            BNEXT = B2
            NBCURR = NB( STEP )
          ELSE
*
*              --- automatic bandwidth selection ---
*
            BNEXT = NBDFLT( 'SBRDD', 'DoublePrec', 'Bandwidth',
     >                      N, BCURR, B2, NEEDU )
*
*              --- if the proposed intermediate bandwidth is out
*                  of range then do a one-step reduction         ---
*
            IF ( ( BNEXT .GE. BCURR ) .OR. ( BNEXT .LT. B2 ) ) THEN
              BNEXT = B2
            ENDIF
            NBCURR = 0
          ENDIF
*
          IF ( BNEXT .EQ. 1 ) THEN
*
*              --- reduction to tridiagonal form;
*                  use SBR or LAPACK code ?       ---
*
            NX = NBDFLT( 'SBRDD', 'DoublePrec', 'CrossoverPoint',
     >                   N, BCURR, B2, NEEDU )
            IF ( BCURR .LT. NX ) THEN
*
*                --- LAPACK tridiagonalization; first repack the
*                    band to improve performance                 ---
*
              LDNEXT = BCURR + 1
              IF ( LDCURR .NE. LDNEXT ) THEN
                CALL DSB2BI( ULCURR, N, BCURR, A, LDCURR, LDNEXT, INFO )
                LDCURR = LDNEXT
                ULCURR = 'Lower'
              ENDIF
*
              IF ( LWORK .LT. N ) THEN
                INFO = -17
              ELSE
                CALL DSBTRD( JOBU, ULCURR, N, BCURR, A, LDCURR, D, E,
     >                       U, LDU, WORK, INFO )
              ENDIF
            ELSE
*
*                --- SBR tridiagonalization; first repack the band
*                    to improve performance                        ---
*
              LDNEXT = MIN( 2*BCURR, LDA )
              IF ( ( LDCURR .NE. LDNEXT ) .OR.
     >             LSAME( ULCURR, 'Upper' ) ) THEN
                CALL DSB2BI( ULCURR, N, BCURR, A, LDCURR, LDNEXT, INFO )
                LDCURR = LDNEXT
                ULCURR = 'Lower'
              ENDIF
*
              CALL DSBRDT( JOBU, N, BCURR, A, LDCURR, DRPTOL, D, E,
     >                     U, LDU, NBCURR, WORK, LWORK, INFO )
              IF ( INFO .EQ. -13 )     INFO = -17
            ENDIF
*
          ELSE
*
*              --- SBR reduction to narrower banded form; first
*                  repack the band to improve performance (optimum
*                  leading dimension depends on the blocking
*                  factor - so determine that one first)           ---
*
            IF ( NBCURR .EQ. 0 ) THEN
              NX = NBDFLT( 'SBRDB', 'DoublePrec', 'BlockingFactor',
     >                     N, BCURR, BNEXT, NEEDU )
            ELSE
              NX = NBCURR
            ENDIF
            IF ( NX .GT. BNEXT )     NX = BNEXT
            IF ( NX .LT. 1 )     NX = 1
            LDNEXT = MIN( BCURR+(BCURR-BNEXT)+NX, LDA )
            IF ( ( LDCURR .NE. LDNEXT ) .OR.
     >           LSAME( ULCURR, 'Upper' ) ) THEN
              CALL DSB2BI( ULCURR, N, BCURR, A, LDCURR, LDNEXT, INFO )
              LDCURR = LDNEXT
              ULCURR = 'Lower'
            ENDIF
*
            CALL DSBRDB( JOBU, N, BCURR, BNEXT, A, LDCURR, DRPTOL,
     >                   U, LDU, NBCURR, WORK, LWORK, INFO )
            IF ( INFO .EQ. -12 )     INFO = -17
          ENDIF
*
          BCURR = BNEXT
          IF ( INFO .GE. 0 )     GOTO 100
*
        ENDIF
*
*          --- after the reduction, repack the band to leading
*              dimension lda and lower banded storage scheme   ---
*
        IF ( ( INFO .GE. 0 ) .AND.
     >       ( ( LDCURR .NE. LDA ) .OR. LSAME( ULCURR, 'Upper' ) ) )
     >  THEN
          CALL DSB2BI( ULCURR, N, BCURR, A, LDCURR, LDA, INFO )
        ENDIF
*
        IF ( INFO .GE. 0 )     INFO = 1
*
  999   RETURN
        END
*
* **********************************************************************
*
        SUBROUTINE DSYRDD( JOBU, UPLO, N, B2, A, LDA, DRPTOL,
     >                     ABAND, LDBAND, D, E, NSTEPS, B, NB,
     >                     WORK, LWORK, STEP, INFO )
*
* Description:
*
*   Driver routine for the multi-step reduction of symmetric matrices to
*   banded form.
*   Given a symmetric n-by-n matrix A, dsyrdd reduces A to a symmetric
*   banded matrix B with b2 sub(or super-) diagonals using a sequence of
*   orthogonal transformations.
*   If desired, the same transformations are also accumulated in an
*   orthogonal matrix U. That is, the resulting banded matrix B and the
*   accumulated matrix U fulfill the equation
*
*       T
*      U  * A * U  =  B  .
*
*   A must be given in the LAPACK symmetric storage scheme with either
*   the upper or lower triangle stored explicitly, and B is returned in
*   lower banded storage (upper banded is not supported for output).
*
* Author: Bruno Lang
*         Aachen University of Technology
*         na.blang@na-net.ornl.gov
* Date: May 17, 2000
* Version: SBR Toolbox, Rev. 1.4.1
*
* Parameters:
*
        CHARACTER*1           JOBU, UPLO
        INTEGER               N, B2, LDA, LDBAND, NSTEPS, LWORK, STEP,
     >                        INFO
        INTEGER               B( * ), NB( * )
        DOUBLE PRECISION      DRPTOL
        DOUBLE PRECISION      A( LDA, * ), ABAND( LDBAND, * ),
     >                        D( * ), E( * ), WORK( * )
*
*   jobu    (in) character*1
*           Specifies whether U is required or not.
*           jobu = 'U' : Build U.
*                = 'N' : Do not build U.
*
*   uplo    (in) character*1
*           Is the matrix A stored in the upper or lower triangle of a ?
*           uplo = 'U' : Upper triangle.
*                = 'L' : Lower triangle.
*
*   n       (in) integer
*           The size of the matrix A.
*           n >= 0.
*
*   b2      (in) integer
*           The number of sub(super-)diagonals of B.
*           1 <= b2 < n, if n > 1.
*           b2 = 0     , if n <= 1.
*
*   a       (in/out) double precision array, dimension ( lda, n )
*           On entry, the leading n-by-n portion of a contains the
*           symmetric matrix A in upper or lower symmetric storage.
*           If uplo='U' then the upper triangle of A is stored in
*           the upper triangle of the leading n-by-n portion of a;
*           if uplo='L' then the lower triangle of A is stored in
*           the lower triangle of the leading n-by-n portion of a.
*           On exit, A is destroyed.
*           If jobu='U' then the leading n-by-n portion of a holds
*           the orthogonal matrix U.
*
*   lda     (in) integer
*           The leading dimension of the array a.
*           lda >= max( n, 1 ).
*
*   drptol  (in) double precision
*           Threshold for dropping the Householder transformations: if
*           the norm of the vector to be eliminated is already smaller
*           than drptol then the transform is skipped.
*           drptol >= 0.0.
*           If you do not know, use 0.0.
*
*   aband   (out) double precision array, dimension ( ldband, n )
*           On exit, the leading (b2+1)-by-n portion of aband holds the
*           symmetric banded matrix B in lower banded storage (i.e., for
*           j = 1, ..., n and i = j, ..., min( j+b2, n ), the matrix
*           element B( i, j ) is stored in location aband( i-j+1, j ).
*
*   ldband  (in) integer
*           The leading dimension of the array aband.
*           ldband >= b( 1 ) + 1, if nsteps > 1.
*                  >= b2 + 1    , otherwise.
*           For reasonable performance, set ldband >= 2 * b( 1 ) if
*           nsteps > 1, and ldband >= 50 if nsteps = 0.
*
*   d       (out) double precision array, dimension ( n )
*           If b2 = 1 then, on exit, d( 1 .. n ) contains the main
*           diagonal of the resulting tridiagonal matrix B.
*
*   e       (out) double precision array, dimension ( n-1 )
*           If b2 = 1 then, on exit, e( 1 .. n-1 ) contains the
*           subdiagonal of the resulting tridiagonal matrix B.
*
*   nsteps  (in) integer
*           The number of reduction steps.
*           nsteps >= 0.
*
*           If nsteps > 0 then the reduction from full to b2 diagonals
*           is done in exactly nsteps steps with user-specified
*           intermediate bandwidths and blocking factors (see paramaters
*           b and nb, below).
*           If nsteps = 0 then the routine tries to determine these
*           parameters on its own.
*           If you do not know or if you do not care, use 0.
*
*   b       (in) integer array, dimension ( nsteps-1 )
*           The sequence of intermediate bandwidths. The nsteps
*           bandwidth reduction steps are
*              n-1 --> b( 1 ),
*              b( 1 ) --> b( 2 ),
*              ...
*              b( nsteps-1 ) --> b2.
*           This array is accessed only if nsteps > 1.
*
*   nb      (in) integer array, dimension ( nsteps )
*           The blocking factors nb, as suggested by the user, for the
*           reduction and/or the accumulation of the transformations
*           in each of the nb reduction steps.
*           Note that the automatic selection of suitable blocking
*           factors is enabled by setting nb( i ) = 0 for some steps i.
*           This array is accessed only if nsteps > 0.
*
*   work    (workspace) double precision array, dimension ( lwork )
*
*   lwork   (in) integer
*           The length of the workspace array.
*           The exact amount of workspace required depends on the
*           sequence of intermediate bandwidths (i.e., b), on the
*           leading dimension of the array a (i.e., lda), and on the
*           blocking factors for the reduction steps (i.e., nb).
*           lwork >= maximum workspace requirement for any of the
*                    reduction steps.
*           Minimum requirement is 4 * n, but setting lwork this low may
*           severely degrade performance.
*           As a rule, lwork >= 50 * n should work fine; if you can
*           spare it, use lwork >= 100 * n.
*
*   step    (out) integer
*           On exit, step returns the number of the last reduction step
*           that was performed (or tried).
*
*   info    (out) integer
*           On exit, info indicates the status of the reduction.
*           info =   1 : All arguments were OK, and the reduction went
*                        without problems.
*                =  -1 : jobu is none of 'U', 'N' (upper/lower case).
*                =  -2 : uplo is none of 'U', 'L' (upper/lower case).
*                =  -3 : n is out of range (negative).
*                =  -4 : b2 is out of range ( < 0,
*                        or < 1 or >= n while n > 1 ).
*                =  -6 : lda is too small (see above).
*                =  -7 : drptol is out of range (negative).
*                =  -9 : ldband is out of range ( < b(1)+1 ).
*                = -12 : nsteps is out of range (negative).
*                = -13 : the sequence of bandwidths n-1, b( 1 ), b( 2 ),
*                        ..., b( nsteps-1 ), b2 is not monotonically
*                        decreasing
*                = -14 : nb contains negative entries
*                = -16 : lwork is too small (see above).
*           If info < 0 and step > 0 then the parameter for the step-th
*           reduction step was not OK.
*
* Local constants and variables:
*
        DOUBLE PRECISION      ZERO
        PARAMETER             ( ZERO = 0.0D0 )
*
        LOGICAL               NEEDU, UPPER
        INTEGER               BNEXT, NBCURR
*
*   needu   is U needed ?
*   upper   is the upper triangle of A stored ?
*   nbcurr  blocking factor for the current step
*   bnext   bandwidth after the next reduction step
*
* Routines called:
*
        LOGICAL               LSAME
        INTEGER               NBDFLT
        EXTERNAL              LSAME, NBDFLT
*
*   lsame   case-insensitive character matching (BLAS)
*   nbdflt  determine intermediate bandwidth (SBR)
*
        EXTERNAL              DORGTR, DSBRDD, DSY2BC, DSYGTR, DSYRDB,
     >                        DSYTRD
*
*   dorgtr  accumulation of the transformations (LAPACK)
*   dsbrdd  bandwidth reduction for banded matrices (SBR)
*   dsy2bc  repacking from full to banded storage (SBR)
*   dsygtr  accumulation of the transformations (SBR)
*   dsyrdb  reduction full -> banded (SBR)
*   dsytrd  one-step reduction full -> tridiagonal (LAPACK)
*
        INTRINSIC             MAX
*
* ----------------------------------------------------------------------
*
        STEP = 1
*
*          --- check for errors in input ---
*
        NEEDU = LSAME( JOBU, 'U' )
        UPPER = LSAME( UPLO, 'U' )
*
        IF ( .NOT. ( NEEDU .OR. LSAME( JOBU, 'N' ) ) ) THEN
          INFO = -1
        ELSEIF ( .NOT. ( UPPER .OR. LSAME( UPLO, 'Lower' ) ) ) THEN
          INFO = -2
        ELSEIF ( N .LT. 0 ) THEN
          INFO = -3
        ELSEIF ( ( B2 .LT. 0 ) .OR.
     >           ( ( N .GT. 1 ) .AND.
     >             ( ( B2 .LT. 1 ) .OR. ( B2 .GT. N-1 ) ) ) ) THEN
          INFO = -4
        ELSEIF ( LDA .LT. MAX( N, 1 ) ) THEN
          INFO = -6
        ELSEIF ( DRPTOL .LT. ZERO ) THEN
          INFO = -7
        ELSEIF ( ( LDBAND .LE. B2 ) .OR.
     >           ( ( NSTEPS .GT. 1 ) .AND. ( LDBAND .LE. B( 1 ) ) ) )
     >  THEN
          INFO = -9
        ELSEIF ( NSTEPS .LT. 0 ) THEN
          INFO = -12
        ELSEIF ( LWORK .LT. 4*N ) THEN
          INFO = -16
        ELSE
          INFO = 0
        ENDIF
        IF ( INFO .LT. 0 )     GOTO 999
*
        IF ( NSTEPS .GT. 1 ) THEN
          IF ( ( B( 1 ) .GT. ( N-1 ) ) .OR.
     >         ( B2 .GT. B( NSTEPS-1 ) ) ) THEN
            INFO = -13
          ELSE
            DO 10 STEP = 1, NSTEPS-2
              IF ( B( STEP+1 ) .GT. B( STEP ) )     INFO = -13
   10       CONTINUE
          ENDIF
        ENDIF
        IF ( INFO .LT. 0 )     GOTO 999
*
        IF ( NSTEPS .GT. 0 ) THEN
          DO 20 STEP = 1, NSTEPS
            IF ( NB( STEP ) .LT. 0 )     INFO = -14
   20     CONTINUE
        ENDIF
        IF ( INFO .LT. 0 )     GOTO 999
*
*          --- check for quick return ---
*
        IF ( N .EQ. 0 )     GOTO 999
*
        STEP = 1
*
*          --- determine first intermediate bandwidth and blocking
*              factor                                              ---
*
        IF ( NSTEPS .GT. 0 ) THEN
          IF ( NSTEPS .EQ. 1 ) THEN
            BNEXT = B2
          ELSE
            BNEXT = B( 1 )
          ENDIF
          NBCURR = NB( 1 )
        ELSE
*
*            --- automatic bandwidth selection ---
*
          BNEXT = NBDFLT( 'SYRDD', 'DoublePrec', 'Bandwidth',
     >                    N, N - 1, B2, NEEDU )
*
*              --- if the proposed intermediate bandwidth is out
*                  of range then do a one-step reduction         ---
*
          IF ( ( BNEXT .GT. N-1 ) .OR. ( BNEXT .LT. B2 ) .OR.
     >         ( LDBAND .LT. 2*BNEXT ) ) THEN
            BNEXT = B2
          ENDIF
          NBCURR = 0
        ENDIF
*
        IF ( BNEXT .EQ. 1 ) THEN
*
*            --- use LAPACK code for direct tridiagonalization ---
*
          CALL DSYTRD( UPLO, N, A, LDA, D, E, WORK( 1 ),
     >                 WORK( N ), LWORK-N+1, INFO )
          IF ( NEEDU ) THEN
            CALL DORGTR( UPLO, N, A, LDA, WORK( 1 ),
     >                   WORK( N ), LWORK-N+1, INFO )
          ENDIF
        ELSE
*
*            --- reduction to bandwidth > 1 width dsyrdb ---
*
          CALL DSYRDB( UPLO, 'NoU', N, BNEXT, A, LDA, DRPTOL,
     >                 WORK( 1 ), N, NBCURR, WORK( 1 ),
     >                 WORK( N+1 ), LWORK-N, INFO )
*
*            --- repack the reduced matrix to banded storage ---
*
          CALL DSY2BC( UPLO, N, BNEXT, A, LDA, ABAND, LDBAND, INFO )
*
*            --- accumulation of the transformations, if required ---
*
          IF ( NEEDU ) THEN
            CALL DSYGTR( UPLO, N, BNEXT, A, LDA, WORK( 1 ),
     >                   WORK( N+1 ), LWORK-N, INFO )
          ENDIF
*
*            --- further reduce the banded matrix, if required ---
*
          IF ( ( BNEXT .GT. B2 ) .AND. ( INFO .GE. 0 ) ) THEN
            IF ( NSTEPS .EQ. 0 ) THEN
              CALL DSBRDD( JOBU, 'Lower', N, BNEXT, B2, ABAND, LDBAND,
     >                     DRPTOL, D, E, A, LDA, 0, B, NB,
     >                     WORK, LWORK, STEP, INFO )
            ELSE
              CALL DSBRDD( JOBU, 'Lower', N, BNEXT, B2, ABAND, LDBAND,
     >                     DRPTOL, D, E, A, LDA, NSTEPS-1, B( 2 ),
     >                     NB( 2 ), WORK, LWORK, STEP, INFO )
            ENDIF
            STEP = STEP + 1
          ENDIF
        ENDIF
*
  999   IF ( INFO .GE. 0 )     INFO = 1
*
        RETURN
        END
*
* **********************************************************************
*
        SUBROUTINE DLBRFG( N, ALPHA, X, INCX, TAU, DRPTOL )
*
* ----------------------------------------------------------------------
*
* Description:
*
*   dlbrfg generates an elementary reflector H of order n, such that
*
*                                    T
*      H * ( alpha ) = ( beta ) ,   H  * H = I ,
*          (   x   )   (   0  )
*
*   where alpha and beta are scalars, and x is an ( n - 1 )-element
*   vector. H is represented in the form
*
*                                 T
*      H = I - tau * ( 1 ) * ( 1 v  ) ,
*                    ( v )
*
*   where tau is a scalar and v is an ( n - 1 )-element vector.
*
*   If the elements of x are all zero, or the 2-norm of x does not
*   exceed drptol, then x is simply zeroed out, tau is set to 0,
*   indicating that H is the identity matrix. In this case, the
*   transformation is effectively skipped.
*
*   Otherwise  1 <= tau <= 2.
*
*   This routine is a modified version of dlbrfg by Bischof/Sun which
*   in turn is a modified version of the LAPACK routine dlarfg.
*
* Author: Bruno Lang
*         Aachen University of Technology
*         na.blang@na-net.ornl.gov
* Date: May 17, 2000
* Version: SBR Toolbox, Rev. 1.4.1
*
* Note: >>> This routine performs no sanity checks on its arguments. <<<
*
* Parameters:
*
        INTEGER               N, INCX
        DOUBLE PRECISION      ALPHA, TAU, DRPTOL
        DOUBLE PRECISION      X( * )
*
*   n       (in) integer
*           The order of the elementary reflector.
*           n >= 0.
*
*   alpha   (in/out) double precision
*           On entry, the value alpha.
*           On exit, it is overwritten with the value beta.
*
*   x       (in/out) double precision array, dimension ( 1+(n-2)*incx )
*           On entry, the vector x.
*           On exit, the "x elements" are overwritten with the vector v.
*
*   incx    (in) integer
*           The increment between elements of x.
*           incx > 0.
*
*   tau     (out) double precision
*           The scaling factor tau.
*
*   drptol  (in) double precision
*           Threshold for dropping Householder transforms.
*           If the norm of the vector to eliminate is already <= drptol
*           then the transform is skipped.
*           drptol >= 0.0.
*           If you do not know, use 0.0.
*
* ----------------------------------------------------------------------
*
* Local variables:
*
        DOUBLE PRECISION      ZERO, ONE
        PARAMETER             ( ZERO = 0.0D0, ONE = 1.0D0 )
*
        DOUBLE PRECISION      XNORM, BETA, SAFMIN, RSAFMN
        INTEGER               CNT, J
*
*   xnorm   2-norm of x
*   beta    the computed value beta
*   safmin  safe minimum (such that relative accuracy can be guaranteed)
*   rsafmn  1 / safmin
*   cnt     counts the rescaling steps
*
* Routines called:
*
        DOUBLE PRECISION      DLAMCH, DLAPY2, DNRM2
        EXTERNAL              DLAMCH, DLAPY2, DNRM2, DSCAL
*
*   dlamch        determine machine parameters (LAPACK)
*   dlapy2        stable 2-norm of 2-vector (LAPACK)
*   dnrm2         2-norm of a vector (BLAS)
*   dscal         vector scaling (BLAS)
*
        INTRINSIC             SIGN, ABS
*
* ----------------------------------------------------------------------
*
*            --- check for quick exit ---
*
        IF ( N .LE. 1 ) THEN
          TAU = ZERO
        ELSE
*
          XNORM = DNRM2( N-1, X, INCX )
*
          IF ( XNORM .LE. DRPTOL ) THEN
*
*              --- H = I ---
*
            TAU = ZERO
            DO 10 J = 1, N-1
              X( 1+(J-1)*INCX ) = ZERO
   10       CONTINUE
*
          ELSE
*
*              --- general case ---
*
            BETA = - SIGN( DLAPY2( ALPHA, XNORM ), ALPHA )
            SAFMIN = DLAMCH( 'S' ) / DLAMCH( 'E' )
*
            IF ( ABS( BETA ) .LT. SAFMIN ) THEN
*
*                --- xnorm, beta may be inaccurate; scale x
*                    and recompute them                     ---
*
              RSAFMN = ONE / SAFMIN
              CNT = 0
*
  100         CONTINUE
                CNT = CNT + 1
                ALPHA = ALPHA * RSAFMN
                CALL DSCAL( N-1, RSAFMN, X, INCX )
                BETA = BETA * RSAFMN
              IF ( ABS( BETA ) .LT. SAFMIN )     GOTO 100
*
*                --- now abs( beta ) is at least safmin ---
*
              XNORM = DNRM2( N-1, X, INCX )
              BETA = - SIGN( DLAPY2( ALPHA, XNORM ), ALPHA )
              TAU = ( BETA - ALPHA ) / BETA
              CALL DSCAL( N-1, ONE / ( ALPHA-BETA ), X, INCX )
*
*                --- if alpha is subnormal, it may lose
*                    relative accuracy                  ---
*
              ALPHA = BETA
*
              DO 110 J = 1, CNT
                ALPHA = ALPHA * SAFMIN
  110         CONTINUE
*
            ELSE
*
              TAU = ( BETA - ALPHA ) / BETA
              CALL DSCAL( N-1, ONE / ( ALPHA-BETA ), X, INCX )
              ALPHA = BETA
*
            ENDIF
*
          ENDIF
*
        ENDIF
*
        RETURN
        END
*
* **********************************************************************
*
        SUBROUTINE DSYRF( UPLO, N, V, INCV, TAU,
     >                    A, LDA, WORK )
*
* ----------------------------------------------------------------------
*
* Description:
*
*                                                      T
*   Apply the Householder transform H = I - tau * v * v  from both
*   sides to the symmetric matrix A, i.e., update
*
*                            T
*      A := H * A * H   ( = H  * A * H )  .
*
* Author: Bruno Lang
*         Aachen University of Technology
*         na.blang@na-net.ornl.gov
* Date: May 17, 2000
* Version: SBR Toolbox, Rev. 1.4.1
*
* Note: >>> This routine performs no sanity checks on its arguments. <<<
*
* ----------------------------------------------------------------------
*
* Parameters:
*
        CHARACTER*1           UPLO
        INTEGER               N, INCV, LDA
        DOUBLE PRECISION      TAU
        DOUBLE PRECISION      V( * ), A( LDA, * ), WORK( * )
*
*  uplo     (in) character*1
*           Is the upper or lower triangle of A stored ?
*           uplo = 'U' : Upper triangle.
*                = 'L' : Lower triangle.
*
*  n        (in) integer
*           The dimension of the matrix A.
*           n >= 0.
*
*  v        (in) double precision array, dimension ( 1+(n-1)*incv )
*           The Householder vector associated with H.
*           Usually, v( 1 ) = 1, and v( 2:n ) contains the "mantissa" of
*           the transform.
*
*  incv     (in) integer
*           The increment between the elements of v.
*           incv > 0.
*
*  tau      (in) double precision
*           The scaling factor for the transform.
*
*  a        (in/out) double precision array, dimension ( lda, n )
*           On entry, an n-by-n symmetric matrix A with either the lower
*           or the upper triangle stored.
*           On exit, the leading A is overwritten with H * A * H.
*
*  lda      (in) integer
*           The leading dimension of the array a.
*           lda >= max( n, 1 ).
*
*  work     (workspace) double precision array, dimension ( n )
*
* ----------------------------------------------------------------------
*
* Local variables:
*
        DOUBLE PRECISION      ZERO, ONE
        PARAMETER             ( ZERO = 0.0D0, ONE = 1.0D0 )
*
        DOUBLE PRECISION      ALPHA
*
*   alpha   auxiliary factor
*
* Routines called:
*
        DOUBLE PRECISION      DDOT
        EXTERNAL              DAXPY, DDOT, DSYMV, DSYR2
*
*   daxpy   add multiple of a vector to another vector (BLAS)
*   ddot    dot product of two vectors (BLAS)
*   dsymv   symmetric matrix-vector product (BLAS)
*   dsyr2   symmetric rank-2 update (BLAS)
*
* ----------------------------------------------------------------------
*
        IF ( ( N .GT. 0 ) .AND. ( TAU .NE. ZERO ) ) THEN
*
*              --- work( 1:n ) := x = tau * A * v ---
*
          CALL DSYMV( UPLO, N, TAU, A, LDA, V, INCV, ZERO, WORK, 1 )
*
*                                                        T
*              --- work( 1:n ) := w = x - 1/2 * tau * ( x  * v ) * v ---
*
          ALPHA = -0.5D0 * TAU * DDOT( N, WORK, 1, V, INCV )
          CALL DAXPY( N, ALPHA, V, INCV, WORK, 1 )
*
*                                T        T
*              --- A := A - v * w  - w * v  ---
*
          CALL DSYR2( UPLO, N, -ONE, V, INCV, WORK, 1, A, LDA )
*
        ENDIF
*
        RETURN
        END
*
* **********************************************************************
*
        SUBROUTINE DGBQR2( M, N, D, A, LDA, TAU, DRPTOL, WORK )
*
* ----------------------------------------------------------------------
*
* Description:
*
*   Given an m-by-n matrix A with d subdiagonals, dgeqr2 computes
*   a QR factorization of A.
*
*   The matrix Q is represented as a product of elementary reflectors
*
*      Q  =  H( 1 ) * ... * H( k )   ,   where k = min( m, n ).
*
*   Each H( j ) has the form
*
*                                  T
*      H( j )  =  I  -  tau * v * v
*
*   where tau is a scalar, and v is a vector of length m with
*   v( 1:j-1 ) = 0   and   v( j ) = 1.
*   On exit, v( j+1:m ) is stored in a( j+1:m, j ) and tau in tau( j ).
*
* Author: Bruno Lang
*         Aachen University of Technology
*         na.blang@na-net.ornl.gov
* Date: May 17, 2000
* Version: SBR Toolbox, Rev. 1.4.1
*
* Note: >>> This routine performs no sanity checks on its arguments. <<<
*
* Parameters:
*
        INTEGER               M, N, D, LDA
        DOUBLE PRECISION      DRPTOL
        DOUBLE PRECISION      A( LDA, * ), TAU( * ), WORK( * )
*
*   m       (in) integer
*           The number of rows in the matrix A.
*           m >= 0.
*
*   n       (in) integer
*           The number of columns in the matrix A.
*           n >= 0.
*
*   d       (in) integer
*           The number of subdiagonals of the matrix A.
*           0 <= d <= m - 1, if m > 0.
*           d = 0          , if m = 0.
*
*   a       (in/out) double precision array, dimension ( lda, n )
*           On entry, the m-by-n matrix A with lower bandwidth b.
*           On exit, the elements on and above the diagonal contain the
*           min( m, n )-by-n upper trapezoidal matrix R (R is upper
*           triangular if m >= n); the elements below the diagonal,
*           together with the array tau, represent the orthogonal
*           matrix Q as a product of Householder transformations.
*
*   lda     (in) integer
*           The leading dimension of the array a.
*           lda >= min( m, 1 ).
*
*   tau     (out) double precision array, dimension ( min( m-1, n ) )
*           The scaling factors of the Householder transformations.
*
*   drptol  (in) double precision
*           Threshold for dropping Householder transforms.
*           If the norm of the vector to be zeroed is already less than
*           drptol then the transform is skipped.
*           drptol >= 0.0.
*
*   work    (workspace) double precision array, dimension ( n )
*
* ----------------------------------------------------------------------
*
* Local variables:
*
        DOUBLE PRECISION      ZERO, ONE
        PARAMETER             ( ZERO = 0.0D0, ONE = 1.0D0 )
*
        INTEGER               J, LEN
        DOUBLE PRECISION      AJJ
*
*   len     length of the current transform
*   ajj     temporary buffer for the element a( j, j )
*
* Routines called:
*
        EXTERNAL              DLARFX, DLBRFG
*
*   dlarfx  apply Householder transformation (LAPACK)
*   dlbrfg  generate Householder vector (SBR)
*
        INTRINSIC             MIN
*
* ----------------------------------------------------------------------
*
        DO 100 J = 1, MIN( M-1, N )
*
          LEN = MIN( D+1, M-J+1 )
*
*            --- generate Householder reflector H( j ) to zero
*                A( j+1:j+d, j )                               ---
*
          CALL DLBRFG( LEN, A( J, J ), A( J+1, J ), 1,
     >                 TAU( J ), DRPTOL )
*
*            --- apply the transform ---
*
          IF ( TAU( J ) .NE. ZERO ) THEN
*
            AJJ = A( J, J )
            A( J, J ) = ONE
            CALL DLARFX( 'Left', LEN, N-J, A( J, J ), TAU( J ),
     >                  A( J, J+1 ), LDA, WORK )
            A( J, J ) = AJJ
*
          ENDIF
*
  100   CONTINUE
*
        RETURN
        END
*
* **********************************************************************
*
        SUBROUTINE DGEQRL( DECOMP, M, N, A, LDA, DRPTOL,
     >                     TAU, WORK )
*
* ----------------------------------------------------------------------
*
* Description:
*
*   Compute the QR or QL factorization of the m-by-n matrix A,
*   where m >= n.
*
*   The matrix Q is represented as a product of elementary reflectors
*
*      Q  =  H( 1 ) * ... * H( n )        (QR decomposition)
*
*   or
*
*      Q  =  H( n ) * ... * H( 1 )        (QL decomposition).
*
*   Each H( j ) has the form
*
*                                  T
*      H( j )  =  I  -  tau * v * v
*
*   where tau is a scalar, and v is a vector of length m.
*
*   In the QR case,
*     v( 1:j-1 ) = 0   and   v( j ) = 1   and, on exit, v( j+1:m ) is
*     stored in a( j+1:m, j ) and tau in tau( j ).
*
*   In the QL case,
*     v( m-n+j+1:m ) = 0   and   v( m-n+j ) = 1   and, on exit,
*     v( 1:m-n+j-1 ) is stored in a( 1:m-n+j-1, j ) and tau in tau( j ).
*
* Author: Bruno Lang
*         Aachen University of Technology
*         na.blang@na-net.ornl.gov
* Date: May 17, 2000
* Version: SBR Toolbox, Rev. 1.4.1
*
* Note: >>> This routine performs no sanity checks on its arguments. <<<
*
* Parameters:
*
        CHARACTER*1           DECOMP
        INTEGER               M, N, LDA
        DOUBLE PRECISION      DRPTOL
        DOUBLE PRECISION      A( LDA, * ), TAU( * ), WORK( * )
*
*   decomp  (in) character*1
*           Which decomposition is required ?
*           decomp = 'L' : QL decomposition.
*                  = 'R' : QR decomposition.
*
*   m       (in) integer
*           The number of rows of the matrix A.
*           m >= 0.
*
*   n       (in) integer
*           The number of columns of the matrix A.
*           0 <= n <= m.
*
*   a       (in/out) double precision array, dimension ( lda, n )
*           On entry, the m-by-n matrix A.
*           If decomp = 'R' then, on exit, the elements on and above the
*           diagonal contain the n-by-n upper triangular matrix R; the
*           elements below the diagonal, together with the array tau,
*           represent the orthogonal matrix Q as a product of
*           Householder transformations.
*           If decomp = 'L' then, on exit, the n-by-n matrix L is stored
*           in the lower triangle at the "bottom" of A, and the
*           Householder vectors are stored in the remaining trapezoid.
*
*   lda     (in) integer
*           The leading dimension of the array a.
*           lda >= min( m, 1 ).
*
*   tau     (out) double precision array, dimension ( n )
*           The scaling factors of the Householder transformations.
*
*   drptol  (in) double precision
*           Threshold for dropping Householder transforms.
*           If the norm of the vector to be zeroed is already less than
*           drptol then the transform is skipped.
*           drptol >= 0.0.
*
*   work    (workspace) double precision array, dimension ( n )
*
* ----------------------------------------------------------------------
*
* Local variables:
*
        DOUBLE PRECISION      ONE
        PARAMETER             ( ONE = 1.0D0 )
*
        INTEGER               I, J
        DOUBLE PRECISION      TMP
*
* Routines called:
*
        LOGICAL               LSAME
        EXTERNAL              LSAME
*
*   lsame   case-insensitive character matching (BLAS)
*
        EXTERNAL              DLARFX, DLBRFG
*
*   dlarfx  apply Householder transform (LAPACK)
*   dlbrfg  generate Householder vector (SBR)
*
        INTRINSIC             MIN
*
* ----------------------------------------------------------------------
*
        IF ( N .GT. 0 ) THEN
*
          IF ( LSAME( DECOMP, 'L' ) ) THEN
*
*              --- QL factorization ---
*
            DO 100 J = N, 1, -1
*
*                --- generate Householder vector to zero out
*                    A( 1:m-n+j, j)                          ---
*
              I = M - N + J
              CALL DLBRFG( I, A( I, J ), A( 1, J ), 1, TAU( J ),
     >                     DRPTOL )
*
*                --- apply transformation from the left ---
*
              TMP = A( I, J )
              A( I, J ) = ONE
              CALL DLARFX( 'Left', I, J-1, A( 1, J ), TAU( J ),
     >                     A, LDA, WORK )
              A( I, J ) = TMP
*
  100       CONTINUE
*
          ELSE
*
*              --- QR factorization ---
*
            DO 200 J = 1, N
*
*                --- generate Householder vector to zero out
*                    A( j+1:m, j )                           ---
*
              CALL DLBRFG( M-J+1, A( J, J ),
     >                     A( MIN( J+1, M ), J ), 1, TAU( J ),
     >                     DRPTOL )
*
*                --- apply transformation from the left ---
*
              TMP = A( J, J )
              A( J, J ) = ONE
              CALL DLARFX( 'Left', M-J+1, N-J, A( J, J ), TAU( J ),
     >                     A( J, J+1 ), LDA, WORK )
              A( J, J ) = TMP
*
  200       CONTINUE
*
          ENDIF
*
        ENDIF
*
        RETURN
        END
*
* **********************************************************************
*
        SUBROUTINE DGBWYG( M, N, D, A, LDA, TAU,
     >                     W, LDW, Y, LDY, K )
*
* ----------------------------------------------------------------------
*
* Description:
*
*   Given a set of Householder vectors as returned by dgbqr2, this
*   routine returns W and Y such that
*
*                       T
*      Q  =  I  +  W * Y
*
*   where Q is the product of the Householder matrices.
*
*   Each Householder transformation is represented by a vector vj and
*   a scalar tau( j ) such that
*                                         T  T          T
*      H( j )  =  I  -  tau( j ) * ( 1, vj  )  * ( 1, vj  )  .
*
*   The matrix Q is accumulated as
*                                                  T
*      Q( j + 1 )  =  I  +  W( j + 1 ) * Y( j + 1 )   =  Q( j ) * H( j )
*
*   where
*                                            T          T
*      Y( j + 1 )  =  [ Y( j ), y ]   with  y  = [ 0, vj  ]   and
*
*      W( j + 1 )  =  [ W( j ), w ]   with  w = - tau( j ) * Q( j ) * y.
*
* Author: Bruno Lang
*         Aachen University of Technology
*         na.blang@na-net.ornl.gov
* Date: May 17, 2000
* Version: SBR Toolbox, Rev. 1.4.1
*
* Note: >>> This routine performs no sanity checks on its arguments. <<<
*
* Parameters:
*
        INTEGER               M, N, D, LDA, LDW, LDY, K
        DOUBLE PRECISION      A( LDA, * ), TAU( * ),
     >                        W( LDW, * ), Y( LDY, * )
*
*   m       (in) integer
*           The number of rows of the matrix A.
*           m >= 0.
*
*   n       (in) integer
*           The number of columns of the matrix A.
*           n >= 0.
*
*   d       (in) integer
*           The number of subdiagonals of the matrix A (as the vectors
*           vj are stored in these subdiagonals, d is the length of
*           these transforms, minus 1).
*           0 <= d <= m - 1, if m > 0.
*           d = 0          , if m = 0.
*
*   a       (in/out) double precision array, dimension ( lda, n )
*           On entry, the d subdiagonals of the m-by-n matrix A contain
*           the Householder vectors vj.
*           On exit, A is zero below the main diagonal.
*
*   lda     (in) integer
*           The leading dimension of the array A.
*           lda >= max( m, 1 ).
*
*   tau     (in/out) double precision array, dimension ( n )
*           On entry, tau contains the scaling factors for the
*           Householder transforms.
*           On exit, tau is destroyed.
*
*   w       (out) double precision array, dimension ( ldw, n )
*           On exit, the first k columns of w contain the block
*           reflector W (zero below the d-th subdiagonal).
*
*   ldw     (in) integer
*           The leading dimension of the array w.
*           ldw >= max( m, 1 ).
*
*   y       (out) double precision array, dimension ( ldy, n )
*           On exit, the first k columns of y contain the block
*           reflector Y (zero above the main diagonal and below the d-th
*           subdiagonal).
*
*   ldy     (in) integer
*           The leading dimension of the array y.
*           ldy >= max( m, 1 ).
*
*   k       (out) integer
*           On exit, k is the "active blocksize"
*           = the number of valid columns in W and Y
*           = the number of nontrivial Householder transforms contained
*             in W and Y
*           = the number of nonzero tau entries in the input.
*           0 <= k <= n.
*
* ----------------------------------------------------------------------
*
* Local variables:
*
        DOUBLE PRECISION      ZERO, ONE
        PARAMETER             ( ZERO = 0.0D0, ONE = 1.0D0 )
*
        INTEGER               I, J, LEN
*
*   len     length of the vector v( j )
*
* Routines called:
*
        EXTERNAL              DGEMV, DSWAP
*
*   dgemv   matrix-vector product (BLAS)
*   dswap   swap two vectors (BLAS)
*
        INTRINSIC             MIN
*
* ----------------------------------------------------------------------
*
        K = 0
*
        DO 200 J = 1, MIN( M-1, N )
*
          IF ( TAU( J ) .NE. ZERO ) THEN
*
*              --- W and Y grow ---
*
            K = K + 1
            LEN = MIN( M-J, D )
*
*              --- copy vj to Y( :, k ) and zero A( j+1:j+d, j ) ---
*
            DO 110 I = 1, M
              Y( I, K ) = ZERO
  110       CONTINUE
            Y( J, K ) = ONE
            CALL DSWAP( LEN, A( J+1, J ), 1, Y( J+1, K ), 1 )
*
*                --- W( :, k ) = - tau( j ) * Y( :, k ) ---
*
            DO 120 I = 1, J-1
              W( I, K ) = ZERO
  120       CONTINUE
            DO 130 I = J, J+LEN
              W( I, K ) = -TAU( J ) * Y( I, K )
  130       CONTINUE
            DO 140 I = J+LEN+1, M
              W( I, K ) = ZERO
  140       CONTINUE
*
*              --- W( :, k ) = - tau( j ) * Q( j-1 ) * y( j )
*                                                    T
*                         = ( I + W( j-1 ) * Y( k-1 )  ) * W( :, k ) ---
*
            CALL DGEMV( 'Transpose', LEN+J, K-1,
     >                  ONE, Y, LDY, W( 1, K ), 1, ZERO, TAU, 1 )
            CALL DGEMV( 'NoTranspose', LEN+J, K-1,
     >                  ONE, W, LDW, TAU, 1, ONE, W( 1, K ), 1 )
*
          ENDIF
*
  200   CONTINUE
*
        RETURN
        END
*
* **********************************************************************
*
        SUBROUTINE DGEWYG( UPLO, M, N, A, LDA, TAU,
     >                     W, LDW, Y, LDY, K,
     >                     WORK )
*
* ----------------------------------------------------------------------
*
* Description:
*
*   Given a set of Householder vectors as returned by dgeqrl, this
*   routine returns W and Y such that
*
*                       T
*      Q  =  I  +  W * Y
*
*   where Q is the product of the Householder matrices.
*
*   Each Householder transformation is represented by a vector vj and
*   a scalar tau( j ) such that
*                                         T  T          T
*      H( j )  =  I  -  tau( j ) * ( 1, vj  )  * ( 1, vj  )  .
*
*   The matrix Q is accumulated as
*                                                  T
*      Q( j + 1 )  =  I  +  W( j + 1 ) * Y( j + 1 )   =  Q( j ) * H( j )
*
*   where
*                                            T          T
*      Y( j + 1 )  =  [ Y( j ), y ]   with  y  = [ 0, vj  ]   and
*
*      W( j + 1 )  =  [ W( j ), w ]   with  w = - tau( j ) * Q( j ) * y.
*
* Author: Bruno Lang
*         Aachen University of Technology
*         na.blang@na-net.ornl.gov
* Date: May 17, 2000
* Version: SBR Toolbox, Rev. 1.4.1
*
* Note: >>> This routine performs no sanity checks on its arguments. <<<
*
* Parameters:
*
        CHARACTER*1           UPLO
        INTEGER               M, N, LDA, LDW, LDY, K
        DOUBLE PRECISION      A( LDA, * ), TAU( * ),
     >                        W( LDW, * ), Y( LDY, * ), WORK( * )
*
*   uplo    (in) character*1
*           Are the Housholder vectors stored in the upper or lower
*           triangle of A ?
*           uplo = 'U' : Upper triangle.
*                = 'L' : Lower triangle.
*
*   m       (in) integer
*           The row dimension of the matrix A.
*           m >= 0.
*
*   n       (in) integer
*           The column dimension of A.
*           0 <= n <= m.
*
*   a       (in) double precision array, dimension ( lda, n )
*           The matrix containing the Householder vectors.
*
*   lda     (in) integer
*           The leading dimension of the array a.
*           lda >= max( m, 1 ).
*
*   tau     (in) double precision array, dimension ( n )
*           The scaling factors for the Householder transformations.
*
*   w       (out) double precision array, dimension ( ldw, n )
*           On exit, the first k columns of w contain the block
*           reflector W.
*
*   ldw     (in) integer
*           The leading dimension of the array w.
*           ldw >= max( m, 1 ).
*
*   y       (out) double precision array, dimension ( ldy, n )
*           On exit, the first k columns of y contain the block
*           reflector Y.
*
*   ldy     (in) integer
*           The leading dimension of the array y.
*           ldy >= max( m, 1 ).
*
*   k       (out) integer
*           On exit, k is the "active blocksize"
*           = the number of valid columns in W and Y
*           = the number of nontrivial Householder transforms contained
*             in W and Y
*           = the number of nonzero tau entries in the input.
*           0 <= k <= n.
*
*   work    (workspace) double precision array, dimension ( n-1 )
*
* ----------------------------------------------------------------------
*
* Local variables:
*
        DOUBLE PRECISION      ZERO, ONE
        PARAMETER             ( ZERO = 0.0D0, ONE = 1.0D0 )
*
        INTEGER               J, I
        DOUBLE PRECISION      TMP
*
* Routines called:
*
        LOGICAL               LSAME
        EXTERNAL              LSAME
*
*   lsame   case-insensitive character matching (BLAS)
*
        EXTERNAL              DCOPY, DGEMV
*
*   dcopy   BLAS
*   dgemv   BLAS
*
* ----------------------------------------------------------------------
*
        K = 0
*
        IF ( ( M .GT. 0 ) .AND. ( N .GT. 0 ) ) THEN
*
          IF ( LSAME( UPLO, 'Lower' ) ) THEN
*
*              --- vectors come from a QR decomposition ---
*
            DO 200 J = 1, N
*
              TMP = - TAU( J )
*
              IF ( TMP .NE. ZERO ) THEN
                K = K + 1
*
*                  --- Y( :, k ) := vj ---
*
                DO 110 I = 1, J-1
                  Y( I, K ) = ZERO
  110           CONTINUE
                Y( J, K ) = ONE
                CALL DCOPY( M-J, A( J+1, J ), 1, Y( J+1, K ), 1 )
*
*                  --- W( :, k ) := - tau( j ) * vj ---
*
                DO 120 I = 1, J-1
                  W( I, K ) = ZERO
  120           CONTINUE
                DO 130 I = J, M
                  W( I, K ) = TMP * Y( I, K )
  130           CONTINUE
*
*                     --- apply previous transformations ---
*
                IF ( K .GT. 1 ) THEN
                  CALL DGEMV( 'Transpose', M, K-1,
     >                        ONE, Y, LDY, W( 1, K ), 1, ZERO, WORK, 1 )
                  CALL DGEMV( 'NoTranspose', M, K-1,
     >                        ONE, W, LDW, WORK, 1, ONE, W( 1, K ), 1 )
                ENDIF
              ENDIF
*
  200       CONTINUE
*
          ELSE
*
*                --- vectors come from a QL decomposition ---
*
            DO 400 J = N, 1, -1
*
              TMP = - TAU( J )
*
              IF ( TMP .NE. ZERO ) THEN
                K = K + 1
*
*                    --- Y( :, k ) := vj ---
*
                DO 310 I = M-N+J+1, M
                  Y( I, K ) = ZERO
  310           CONTINUE
                Y( M-N+J, K ) = ONE
                CALL DCOPY( M-N+J-1, A( 1, J ), 1, Y( 1, K ), 1 )
*
*                    --- W( :, k ) := - tau( j ) * vj ---
*
                DO 320 I = M-N+J+1, M
                  W( I, K ) = ZERO
  320           CONTINUE
                DO 330 I = 1, M-N+J
                  W( I, K ) = TMP * Y( I, K )
  330           CONTINUE
*
*                     --- apply previous transformations ---
*
                IF ( K .GT. 1 ) THEN
                  CALL DGEMV( 'Transpose', M, K-1,
     >                        ONE, Y, LDY, W( 1, K ), 1, ZERO, WORK, 1 )
                  CALL DGEMV( 'NoTranspose', M, K-1,
     >                        ONE, W, LDW, WORK, 1, ONE, W( 1, K ), 1 )
                ENDIF
              ENDIF
*
  400       CONTINUE
*
          ENDIF
*
        ENDIF
*
        RETURN
        END
*
* **********************************************************************
*
        SUBROUTINE DGEWY( SIDE, M, N, K, A, LDA,
     >                    W, LDW, Y, LDY,
     >                    WORK )
*
* ----------------------------------------------------------------------
*
* Description:
*
*   Apply a WY block transform to the matrix A (from the left or from
*   the right), i.e., update A as
*
*                        T  T
*      A  :=  ( I + W * Y  )  * A   , or
*
*                            T
*      A  :=  A * ( I + W * Y  )    , resp.
*
* Author: Bruno Lang
*         Aachen University of Technology
*         na.blang@na-net.ornl.gov
* Date: May 17, 2000
* Version: SBR Toolbox, Rev. 1.4.1
*
* Note: >>> This routine performs no sanity checks on its arguments. <<<
*
* Parameters:
*
        CHARACTER*1           SIDE
        INTEGER               M, N, K, LDA, LDW, LDY
        DOUBLE PRECISION      A( LDA, * ), W( LDW, * ), Y( LDY, * ),
     >                        WORK( * )
*
*   side    (in) character*1
*           Apply the transforms from the left or from the right ?
*           side = 'L' : From the left.
*                = 'R' : From the right.
*
*   m       (in) integer
*           The number of rows of the matrix A.
*           m >= 0.
*
*   n       (in) integer
*           The number of columns of the matrix A.
*           n >= 0.
*
*   k       (in) integer
*           The number of columns of the matrices W and Y.
*           k >= 0.
*
*   a       (in/out) double precision array, dimension ( lda, n )
*           On entry, the matrix A.
*                               T  T
*           On exit, ( I + W * Y  )  * A, if side = 'L', or
*                                   T
*                    A * ( I + W * Y  ),  if side = 'R'.
*
*   lda     (in) integer
*           The leading dimension of the array a.
*           lda >= max( m, 1 ).
*
*   w       (in) double precision array, dimension ( ldw, k )
*           The matrix W.
*
*   ldw     (in) integer
*           The leading dimension of the array w.
*           ldw >= max( m, 1 ), if side = 'L'.
*               >= max( n, 1 ), if side = 'R'.
*           (W has row dimension m or n, resp.)
*
*   y       (in) double precision array, dimension ( ldy, k )
*           The matrix Y.
*
*   ldy     (in) integer
*           The leading dimension of the array y.
*           ldy >= max( m, 1 ), if side = 'L'.
*               >= max( n, 1 ), if side = 'R'.
*           (Y has row dimension m or n, resp.)
*
*   work    (workspace) double precision array,
*           dimension ( n*k ), if side = 'L'.
*                     ( m*k ), if side = 'R'.
*
* ----------------------------------------------------------------------
*
* Local variables:
*
        DOUBLE PRECISION      ZERO, ONE
        PARAMETER             ( ZERO = 0.0D0, ONE = 1.0D0 )
*
* Routines called:
*
        LOGICAL               LSAME
        EXTERNAL              LSAME
*
*   lsame   case-insensitive character matching (BLAS)
*
        EXTERNAL              DGEMM
*
*   dgemm   matrix-matrix product (BLAS)
*
* ----------------------------------------------------------------------
*
        IF ( ( K .GT. 0 ) .AND. ( M .GT. 0 ) .AND. ( N .GT. 0 ) ) THEN
*
          IF ( LSAME( SIDE, 'L' ) ) THEN
*
*                            T
*                --- work = W  * A ---
*
            CALL DGEMM( 'Transpose', 'NoTranspose', K, N, M,
     >                  ONE, W, LDW, A, LDA, ZERO, WORK, K )
*
*                --- A = Y * work + A ---
*
            CALL DGEMM( 'NoTranspose', 'NoTranspose', M, N, K,
     >                  ONE, Y, LDY, WORK, K, ONE, A, LDA )
*
          ELSE
*
*                --- work = A * W ---
*
            CALL DGEMM( 'NoTranspose', 'NoTranspose', M, K, N,
     >                  ONE, A, LDA, W, LDW, ZERO, WORK, M )
*
*                                    T
*                --- A = A + work * Y  ---
*
            CALL DGEMM( 'NoTranspose', 'Transpose', M, N, K,
     >                  ONE, WORK, M, Y, LDY, ONE, A, LDA )
*
          ENDIF
        ENDIF
*
        RETURN
        END
*
* **********************************************************************
*
        SUBROUTINE DQRUPD( M, N, B, A, LDA, DRPTOL,
     >                     NEEDA1, M1, A1, LDA1, N2, A2,
     >                     NB, Y, W, LDWY, WORK )
*
* ----------------------------------------------------------------------
*
* Description:
*
*   dqrupd forms the QR factorization
*
*      A  =  Q * R
*
*   of the m-by-n block A and overwrites A by R.
*   In addition, the matrices W and Y in the block representation
*
*                     T
*      Q  =  I + W * Y
*
*   are built from the Householder transforms that were used in the
*   QR decomposition, a matrix A1 may be updated via
*
*      A1 := A1 * Q  ,
*
*   and another matrix A2 is updated via
*
*             T
*      A2 := Q  * A2  .
*
* Author: Bruno Lang
*         Aachen University of Technology
*         na.blang@na-net.ornl.gov
* Date: May 17, 2000
* Version: SBR Toolbox, Rev. 1.4.1
*
* Note: >>> This routine performs no sanity checks on its arguments. <<<
* Parameters:
*
        INTEGER               M, N, B, LDA, M1, LDA1, N2, NB, LDWY
        DOUBLE PRECISION      DRPTOL
        LOGICAL               NEEDA1
        DOUBLE PRECISION      A( LDA, * ), A1( LDA1, * ), A2( LDA, * ),
     >                        Y( * ), W( * ), WORK( * )
*
*   m       (in) integer
*           The number of rows in the block A, for which a QR
*           decomposition is sought.
*           m > 0.
*
*   n       (in) integer
*           The number of columns in that block.
*           n > 0.
*
*   b       (in) integer
*           The number of nonzero subdiagonals in the block (the
*           subdiagonals below the b-th are supposed to be zero and are
*           not accessed.)
*           b > 0.
*
*   a       (in/out) double precision array, dimension ( lda, n )
*           On entry, A is a m-by-n matrix with lower bandwidth <= b.
*           Because of the banded storage scheme used in the calling
*           routine, the zero elements below the b-th subdiagonal may
*           not be explicitly stored. Therefore, they are not accessed
*           in this routine.
*           On exit, the subdiagonals 1, ..., b of A are zeroed out, and
*           the upper triangle of A is overwritten with the matrix R of
*           the factorization A = Q * R.
*
*   lda     (in) integer
*           The leading dimension of the array a.
*           lda >= m.
*
*   drptol  (in) double precision
*           Threshold for dropping the Householder transformations: if
*           the norm of the vector to be eliminated is already smaller
*           than drptol then the transformation is skipped.
*           drptol >= 0.0.
*
*   needa1  (in) logical
*           Is A1 required or not ?
*           needa1 = .true.  : Update A1.
*                  = .false. : Do not update A1.
*
*   m1      (in) integer
*           The number of rows of the matrix A1 (the number of columns
*           must be m). Accessed only if A1 is required.
*           m1 >= 0.
*
*   a1      (in/out) double precision array, dimension ( lda1, m )
*           On entry, the m1-by-m matrix A1.
*           On exit, A1 may be overwritten with A1 * Q.
*
*   lda1    (in) integer
*           The leading dimension of the array a1. Accessed only if A1
*           is required.
*           lda1 >= max( m1, 1 ).
*
*   n2      (in) integer
*           The number of columns in the block A2.
*           0 <= n2 <= lda.
*
*   a2      (in/out) double precision array, dimension ( lda, n2 )
*           On entry, the m-by-n2 matrix A2.
*                                            T
*           On exit, A2 is overwritten with Q  * A2.
*
*   nb      (out) integer
*           The number of non-zero Householder vectors contributing to
*           Q.
*           0 <= nb <= min( m-1, n ).
*
*   y       (out) double precision array, dimension ( ldwy, n )           T
*           The m-by-nb matrix Y in the representation Q = I + W * Y .
*
*   w       (out) double precision array, dimension ( ldwy, n )           T
*           The m-by-nb matrix W in the representation Q = I + W * Y .
*
*   ldwy    (out) integer
*           The leading dimension of the arrays w and y.
*           ldwy = m.
*
*   work    (workspace) double precision array,
*           dimension ( max( m1, m )*n + n ),   if A1 is required.
*                     ( m*n + n )           ,   otherwise.
*
* ----------------------------------------------------------------------
*
* Local variables :
*
        INTEGER               ITAU, IWORK
*
*  itau     points to the portion of work that will hold the scaling
*           values tau for the blocked transforms (required size for
*           this buffer : n elements)
*  iwork    points to the portion of work that will be used as workspace
*           for the routines dgbqr2 (n elements) and dgewy
*           ( max( m1, m ) * nb elements).
*
* Routines called:
*
        EXTERNAL              DGBQR2, DGBWYG, DGEWY
*
*   dgbqr2  QR factorization of a banded matrix (SBR)
*   dgbwyg  generate WY factors (SBR)
*   dgewy   apply WY transform (SBR)
*
* ----------------------------------------------------------------------
*
*            --- set pointers to subdivide workspace ---
*
        ITAU = 1
        IWORK = ITAU + N
*
*            --- factorize A = Q * R ---
*
        CALL DGBQR2( M, N, B, A, LDA, WORK( ITAU ), DRPTOL,
     >               WORK( IWORK ) )
*
*            --- build W and Y factors ---
*
        LDWY = M
        CALL DGBWYG( M, N, B, A, LDA, WORK( ITAU ),
     >               W, LDWY, Y, LDWY, NB )
*
        IF ( NB .GT. 0 ) THEN
*
*              --- update A1 ---
*
          IF ( NEEDA1 ) THEN
            CALL DGEWY( 'Right', M1, M, NB, A1, LDA1,
     >                  W, LDWY, Y, LDWY, WORK( IWORK ) )
          ENDIF
*
*              --- premultiply A2 ---
*
          IF ( N2 .GT. 0 ) THEN
            CALL DGEWY( 'Left', M, N2, NB, A2, LDA,
     >                  W, LDWY, Y, LDWY, WORK( IWORK ) )
          ENDIF
        ENDIF
*
        RETURN
        END
*
* **********************************************************************
*
        SUBROUTINE DSYWY( UPLO, N, K, A, LDA, W, LDW, Y, LDY, WORK )
*
* ----------------------------------------------------------------------
*
* Description:
*
*   Apply a two-sided WY block transform to the symmetric matrix A,
*                                   T  T                  T
*   i.e., replace A with ( I + W * Y  )  * A * ( I + W * Y  ).
*
* Author: Bruno Lang
*         Aachen University of Technology
*         na.blang@na-net.ornl.gov
* Date: May 17, 2000
* Version: SBR Toolbox, Rev. 1.4.1
*
* Note: >>> This routine performs no sanity checks on its arguments. <<<
*
* Parameters:
*
        CHARACTER*1           UPLO
        INTEGER               N, K, LDA, LDW, LDY
        DOUBLE PRECISION      A( LDA, * ), W( LDW, * ), Y( LDY, * ),
     >                        WORK( * )
*
*   uplo    (in) character*1
*           Is the upper or lower triangle of A stored ?
*           uplo = 'U' : Upper triangle.
*                = 'L' : Lower triangle.
*
*   n       (in) integer
*           The size of the block to transform.
*           n >= 0.
*
*   k       (in) integer
*           The number of nontrivial Householder transforms in the
*           block WY transform.
*           k >= 0.
*
*   a       (in/out) double precision array, dimension ( lda, n )
*           On entry, the n-by-n symmetric matrix A.
*                                                T  T             T
*           On exit, A is overwritten with ( I+WY  )  * A * ( I+WY  ).
*
*   lda     (in) integer
*           The leading dimension of the array a.
*           lda >= max( n, 1 ).
*
*   w       (in) double precision array, dimension ( ldw, k )
*           Contains the block reflector W.
*
*   ldw     (in) integer
*           The leading dimension of the array w.
*           ldw >= max( n, 1 ).
*
*   y       (in) double precision array, dimension ( ldy, k )
*           Contains the block reflector Y.
*
*   ldy     (in) integer
*           The leading dimension of the array y.
*           ldy >= max( n, 1 ).
*
*   work    (workspace) double precision array, dimension ( (n+k)*k )
*
* ----------------------------------------------------------------------
*
* Local variables:
*
        DOUBLE PRECISION      ZERO, ONE
        PARAMETER             ( ZERO = 0.0D0, ONE = 1.0D0 )
*
* Routines called:
*
        EXTERNAL              DGEMM, DSYMM, DSYR2K
*
*   dgemm   matrix-matrix product (BLAS)
*   dsymm   symmetric matrix-matrix product (BLAS)
*   dsyr2k  symmetric rank-2 update (BLAS)
*
* ----------------------------------------------------------------------
*
        IF ( ( N .GT. 0 ) .AND. ( K .GT. 0 ) ) THEN
*
*              --- work( 1 .. n*k ) := X1 = A * W ---
*
          CALL DSYMM( 'Left', UPLO, N, K,
     >                ONE, A, LDA, W, LDW, ZERO, WORK, N )
*
*                                                           T
*              --- work( n*k+1 .. n*k+k*k ) := X2 = 1/2 * X1  * W ---
*
          CALL DGEMM( 'Transpose', 'NoTranspose', K, K, N,
     >                0.5D0, WORK, N, W, LDW, ZERO, WORK( N*K+1 ), K )
*
*              --- work( 1 .. n*k ) := X3 = X1 + Y * X2 ---
*
          CALL DGEMM( 'NoTranspose', 'NoTranspose', N, K, K,
     >                ONE, Y, LDY, WORK( N*K+1 ), K, ONE, WORK, N )
*
*                                 T        T
*              --- A := A + X3 * Y + Y * X3  ---
*
          CALL DSYR2K( UPLO, 'NoTranspose', N, K,
     >                 ONE, WORK, N, Y, LDY, ONE, A, LDA )
*
        ENDIF
*
        RETURN
        END
