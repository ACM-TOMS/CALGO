*     =================================================================
      SUBROUTINE DGEMG( MATTYP, M, N, ISEED, RTHRESH, A, LDA, VS, WORK )
*
      IMPLICIT NONE
*
*     .. Scalar Arguments ..
      INTEGER            M, MATTYP, N
      DOUBLE PRECISION   RTHRESH
*     ..
*     .. Array Arguments ..
      INTEGER            LDA, ISEED( * )
      DOUBLE PRECISION   A( LDA, * ), VS( * ), WORK( * )
*     ..
*
*  Purpose
*  =======
*
*  DGEMG generates a distributed random double-precision
*  m-by-n matrix A.
*
*  Arguments
*  =========
*
*  MATTYP  (input) INTEGER
*          Determines what kind of matrix is generated.
*          The first 7 are generated with driver xGNTST.
*          The following 6 are generated with driver xQRMTX.
*          The last 6 are generated with driver xLATMS.
*          These matrix types are the same as the ones used in the RRQR
*          papers by C.H.Bischof and G.Quintana-Orti.
*
*  M       (input) INTEGER
*          The number of rows of A.
*
*  N       (input) INTEGER
*          The number of columns of A.
*
*  ISEED   (input/output) INTEGER array, dimension ( 4 )
*          Seed for the random number generator. ISEED(4) must be odd.
*
*  RTHRESH (input) DOUBLE PRECISION
*          The threshold for the inverse of the condition number. Used
*          for determining the number of independent columns of the
*          matrix.
*
*  A       (output) DOUBLE PRECISION array, dimension (LDA,N)
*          The random m-by-n matrix.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A. LDA >= max(1,M).
*
*  VS      (output) DOUBLE PRECISION array, dimension ( MIN( M, N ) )
*          The singular values of matrix A.
*
*  WORK    (workspace) DOUBLE PRECISION array, dimension
*            ( MAX( M*N + 3*M + MAX( M, N ),
*                   MAX( 3*MIN( M, N ), 8 ),
*                   3 * MAX( M, N ) ).
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            CRANK, INFO, IPW, MODE, STRIP
      DOUBLE PRECISION   EPS, GAP, ORTHRESH, SCALE
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MIN
*     ..
*     .. External Functions ..
      INTEGER            IRANK
      DOUBLE PRECISION   DLAMCH
      EXTERNAL           IRANK, DLAMCH
*     ..
*     .. External Subroutines ..
      EXTERNAL           DGNTST, DLATMS, DQRMTX, DSORT
*     ..
*     .. Executable Statements ..
*
*     Get grid parameters.
*
*     Get matrix sizes.
*
      EPS  = DLAMCH( 'Epsilon' )
      IPW  = 1
*
*     Generate test matrix of size n-by-n using
*     test matrix generators indicated by `mattyp'.
*     =============================================
*
      ORTHRESH = RTHRESH
      GAP     = 50.0D+0
      STRIP   = 2
      SCALE   = EPS**1
      IF( ( MATTYP.GE.1 ).AND.( MATTYP.LE.7 ) ) THEN
*
        CALL DGNTST( MATTYP, M, N, ORTHRESH, GAP, STRIP,
     $               ISEED, CRANK, VS, A, LDA,
     $               WORK( IPW ) )
*
      ELSE IF( ( MATTYP.GE.8 ).AND.( MATTYP.LE.13 ) ) THEN
*
        IF( MATTYP.EQ. 8 ) MODE =  2
        IF( MATTYP.EQ. 9 ) MODE = -2
        IF( MATTYP.EQ.10 ) MODE =  3
        IF( MATTYP.EQ.11 ) MODE = -3
        IF( MATTYP.EQ.12 ) MODE =  4
        IF( MATTYP.EQ.13 ) MODE = -4
C       Increased the scale by a great factor to avoid very 
C       small values in the matrix once factorized.
        CALL DQRMTX( 'All', 1.0D+6 * SCALE, M, N, ORTHRESH*GAP,
     $               STRIP, MODE, ISEED, CRANK, VS,
     $               A, LDA, WORK( IPW ) )
        CRANK = IRANK( VS, N, ORTHRESH )
*
      ELSE IF( ( MATTYP.GE.14 ).AND.( MATTYP.LE.19 ) ) THEN
*
        IF( MATTYP.EQ.14 ) MODE =  2
        IF( MATTYP.EQ.15 ) MODE = -2
        IF( MATTYP.EQ.16 ) MODE =  3
        IF( MATTYP.EQ.17 ) MODE = -3
        IF( MATTYP.EQ.18 ) MODE =  4
        IF( MATTYP.EQ.19 ) MODE = -4
        CALL DLATMS( M, N, 'Uniform', ISEED, 'Nonsymmetric',
     $               VS, MODE, GAP/ORTHRESH, ONE, M, N,
     $               'No packing', A, LDA, WORK( IPW ), INFO )
        IF( MODE.LT.0 ) THEN
          CALL DSORT( N, VS, 1, 'Decreasing' )
        END IF
        CRANK = IRANK( VS, N, ORTHRESH )
*
      ELSE 
        PRINT *, ''
        PRINT *, 'ERROR in dgemg: Unknown matrix type: ', mattyp
        PRINT *, ''
      END IF
*
      RETURN
*
*     End of DGEMG
*
      END
*     =================================================================
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
*     $Revision: 1.49 $
*     $Date: 1999/02/17 14:04:19 $
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
      EXTERNAL           DLARAN, IRANK, DNRM2, DLAMCH
      INTEGER            IRANK
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
         RANK = IRANK( S, MN, RTHRESH )
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
         RANK = IRANK( S, MN, RTHRESH )
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
         RANK = IRANK( S, MN, RTHRESH )
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

        RANK = IRANK( S, MN, RTHRESH )
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
*     =================================================================
      SUBROUTINE DQRMTX( OPT, SCALE, M, N, RCOND, WIDTH,
     $                   MODE, ISEED, RANK, S, A, LDA, WORK )
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
*     $Revision: 1.49 $
*     $Date: 1999/02/17 14:04:19 $
*
      CHARACTER*1        OPT
      INTEGER            M, N, WIDTH, LDA, RANK, MODE
      INTEGER            ISEED( 4 )
      DOUBLE PRECISION   SCALE, RCOND, A( LDA, * ), S( * ), WORK( * )
*
*  Purpose: 
*  =======  
*
*  generates a matrix for testing DGEQPF. The independent and 
*  dependent columns of A are arranged in a zebra-like fashion.
*  That is, if m = 5, n = 12, and width = 2, 
*          columns 1:2 are independent
*          columns 3:4 are a linear combination of columns 1:2
*          columns 5:6 are independent
*          columns 7:8 are a linear combination of columns 5:6 or 
*                      [1:6], depending on the value of 'opt'.
*          column  9   is independent (there can't be more than
*                      min(m,n) independent columns)
*          columns 10:12 are again linear combinations of previous 
*                      columns
*
*  Arguments:
*  =========
*
*  OPT     (input) CHARACTER*1
*               OPT == 'l' or 'L': dependent columns are linear
*                                  combinations of the last set of 
*                                  independent columns
*               any other value  : dependent columns are linear 
*                                  combinations of all previous 
*                                  independent columns
*  SCALE   (input) DOUBLE PRECISION
*          dependent columns are a random linear combination of
*          previous ones multiplied by SCALE.
*
*  M       (input) INTEGER
*          The number of rows of the matrix A.
* 
*  N       (input) INTEGER
*          The number of columns of the matrix A.
*
*  RCOND   (input) DOUBLE PRECISION
*          1/RCOND is the condition number of the matrix to be
*          generated. Singular values for the submatrix consisting
*          of independent columns are generated between
*          1 and RCOND dependent on MODE.
*
*  WIDTH   (input) INTEGER
*          The width of a strip of dependent or independent columns.
*
*  MODE    (input) INTEGER
*          is passed to DLATMS to determine how diagonal entries
*          are generated between 1 and RCOND.
*              MODE = {-,+}1 : all diagonal entries are RCOND except for
*                              {last,first} one.
*              MODE = {-,+}2 : all diagonal entries are 1 except for
*                              {first,last} one.
*              MODE = {-,+}3 : exponentially {declining,increasing}
*              MODE = {-,+}4 : arithmetically {decl.,incr.}
*
*  ISEED   (input/output) INTEGER array, dimension(4)
*          Seed for random number generator. ISEED(4) must be odd.
*
*  RANK    (output) INTEGER
*          The number of independent columns generated. Note that
*          this need not necessarily be the numerical rank of A
*          as determined by the SVD due to the permutation generated
*          by adding the columns which are linear combinations of
*          previous ones.
*
*  S       (output) DOUBLE PRECISION array (min(M,N))
*          The singular values of A
*
*  A       (output) DOUBLE PRECISION array, dimension (M,N)
*          matrix with singular value distribution given in S
*          and pattern of dependent/independent columns determined
*          by WIDTH.
*
*  LDA     (input) INTEGER
*          leading dimension of A.
*
*  WORK    (workspace) DOUBLE PRECISION array, 
*          dimension max(3*min(m,n),width*width) if OPT == 'L' or 'l'
*          dimension max(3*min(m,n),width*width*2) otherwise
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE

      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )

*     ..
*     .. Local Scalars ..
      INTEGER            MN, WLAST, NSTRPS, INFO, OFFSET, 
     $                   NCOLS, CLSLFT, I
      DOUBLE PRECISION   DUMMY
*     ..
*     ..
*     .. External Subroutines
      EXTERNAL           DLATMS, DLACPY, LSAME, 
     $                   DGEBD2, DBDSQR, DLARNV
      LOGICAL            LSAME
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MIN, MOD, SQRT
*     .. 
*     .. Executable Statements ..

      MN = MIN( M, N )
*
*     How many strips do fit and what is width of last strip?
*
      WLAST = MOD( N, WIDTH )
      IF( WLAST.EQ.0 ) THEN
         NSTRPS = N/WIDTH
         WLAST = WIDTH
      ELSE
         NSTRPS = N/WIDTH + 1 
      END IF
*
*     What is the rank of A?
*
      IF( MOD( NSTRPS, 2 ).EQ.0 ) THEN
         RANK = MIN( MN, NSTRPS/2*WIDTH )
      ELSE
         RANK = MIN( MN, (NSTRPS-1)/2*WIDTH + WLAST )
      END IF
*
*     How many strips is the matrix of size m-by-rank partitioned into?
*
      WLAST = MOD( RANK, WIDTH )
      IF( WLAST.EQ.0 ) THEN
         NSTRPS = RANK/WIDTH
         WLAST = WIDTH
      ELSE
         NSTRPS = RANK/WIDTH + 1
      END IF
*
*     Generate 'rank' independent columns in
*     A(:,(nstrips-1)*width+1:(nstrips-1)*width+rank))
*
      OFFSET = ( NSTRPS-1 )*WIDTH

      CALL DLATMS( M, RANK, 'Uniform', ISEED, 'Nonsymmetric',
     $             S, MODE, ONE/RCOND, ONE, M, RANK, 'No Packing',
     $             A( 1, OFFSET+1 ), LDA, WORK, INFO )
      IF( INFO.GT.0 ) THEN
         WRITE(*,999) INFO
         STOP
      END IF
*
*     Redistribute independent columns and generate dependent
*     ones in columns 1 through offset+rank
*
      DO 10 I = 1,NSTRPS-1

         CALL DLACPY( 'full matrix', M, WIDTH,
     $               A( 1, OFFSET+( I-1 )*WIDTH+1 ), LDA,
     $               A( 1, 2*( I-1 )*WIDTH+1 ),LDA )

         IF( LSAME( OPT, 'L' ) ) THEN
             NCOLS = WIDTH
         ELSE
             NCOLS = MIN( 2, I )*WIDTH
         END IF
         CALL DLARNV( 1, ISEED, NCOLS*WIDTH, WORK( 1 ) )
         CALL DGEMM( 'no transpose', 'no transpose', M, WIDTH, NCOLS,
     $              SCALE, A( 1, ( 2*I-1 )*WIDTH-NCOLS+1 ), LDA,
     $              WORK, NCOLS, ZERO, A( 1,( 2*I-1 )*WIDTH+1 ), LDA )

10    CONTINUE
*
*     generate dependent columns offset+rank+1 through n
*
      CLSLFT = N-( OFFSET+RANK )
      IF( CLSLFT.GT.0 ) THEN
     
         IF( LSAME( OPT, 'L' ) ) THEN
            NCOLS = WLAST
         ELSE
            NCOLS = MIN( OFFSET+RANK, WLAST+WIDTH )
         END IF
         CALL DLARNV( 1, ISEED, NCOLS*CLSLFT, WORK( 1 ) )
         CALL DGEMM( 'no transpose', 'no transpose', M, CLSLFT, NCOLS,
     $              SCALE, A( 1, OFFSET+RANK+1-NCOLS ), LDA, WORK, 
     $              NCOLS, ZERO,A( 1, OFFSET+RANK+1 ), LDA )
      END IF
*
*     compute singular value decomposition of A
*
      CALL DLACPY( 'full matrix', M, N, A, LDA, WORK( MN+1 ), M )
      CALL DGEBD2( M, N, WORK( MN+1 ), M, S, WORK( 1 ),
     $            WORK( MN+M*N+1 ), WORK( 2*MN+M*N+1 ),
     $            WORK( 3*MN+M*N+1 ), INFO )
      CALL DBDSQR( 'upper', MN, 0, 0, 0, S, WORK( 1 ), 
     $            DUMMY, MN, DUMMY, 1, DUMMY, MN, WORK( MN+1 ), INFO )

      RETURN
999   FORMAT( '** ERROR in sqrmtx: DLATMS returns INFO = ',i2 )
*
*     End of DQRMTX
*
      END
*     =================================================================
      SUBROUTINE DSORT( n, x, incrx, job )
      INTEGER     n, incrx
      CHARACTER*1 job
      DOUBLE PRECISION x(*)
*
*     Subroutine to sort a vector
*
*     on entry:
*     ========
*      
*     x    vector of length n to be sorted
*     n    length of vector
*     incrx     element spacing in x
*     job  = 'i' or 'i' sorts in increasing order
*          = 'd' or 'd' sorts in decreasing order
*          otherwise the routine returns without performing
*          any computation
*
*     on exit:
*     ========
*
*     x    sorted in the prescribed order
*
*
*     EXTERNAL entries
*     ================
*
      LOGICAL lsame
      EXTERNAL lsame
*
*     internal variables
*     ==================
*
      INTEGER i, curelt, nextelt, switch,k
      DOUBLE PRECISION temp
       
      switch = 0
      IF( lsame(job,'i')) switch = 1
      IF( lsame(job,'d')) switch = 2
       
      IF( switch .eq. 0) RETURN
       
      GOTO (100,200) switch
*
*     sort in increasing order
*
100   DO 10 i = n-1,1,-1
         k = i
 20      IF( k .eq. n) GOTO 10
            curelt = 1+(k-1)*incrx
            nextelt = 1 + k*incrx
            IF( x(curelt) .le. x(nextelt)) then
               GOTO 10
            ELSE
               temp = x(curelt)
               x(curelt) = x(nextelt)
               x(nextelt) = temp
            END IF
            k = k+1
            GOTO 20
 10   CONTINUE
      RETURN
*
*     sort in decreasing order
*
200   DO 30 i = n-1,1,-1
         k = i
 40      IF( k .eq. n) GOTO 30
            curelt = 1+(k-1)*incrx
            nextelt = 1 + k*incrx
            IF( x(curelt) .ge. x(nextelt)) then
               GOTO 30
            ELSE
               temp = x(curelt)
               x(curelt) = x(nextelt)
               x(nextelt) = temp
            END IF
            k = k+1
            GOTO 40
 30   CONTINUE
      RETURN    
*
*     End of subroutine DSORT
*
      END
*     =================================================================
      INTEGER FUNCTION IRANK( S, N, RCOND )
*
      IMPLICIT NONE
*
*     Returns MAX { 1 <= i <= n | s(1)/s(i) < 1/RCOND }.
*     The entries of S are assumed to be nonnegative and
*     monotonically decreasing.
*
      INTEGER          N
      DOUBLE PRECISION RCOND, S( * )
*
      INTEGER I
      IRANK = 1
      DO I = N, 2, -1
        IF( ( S( 1 )*RCOND ).LT.S( I ) ) THEN
          IRANK = I
          GOTO 20
        END IF
      END DO
20    RETURN
*
*     END OF IRANK
*
      END
*     =================================================================
