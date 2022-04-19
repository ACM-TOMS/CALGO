      SUBROUTINE SQRMTX( OPT, SCALE, M, N, RCOND, WIDTH,
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
*     $Revision: 1.84 $
*     $Date: 96/12/30 16:59:18 $
*
      CHARACTER*1        OPT
      INTEGER            M, N, WIDTH, LDA, RANK, MODE
      INTEGER            ISEED( 4 )
      REAL               SCALE, RCOND, A( LDA, * ), S( * ), WORK( * )
*
*  Purpose:
*  =======
*
*  generates a matrix for testing SGEQPF. The independent and
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
*  SCALE   (input) REAL
*          dependent columns are a random linear combination of
*          previous ones multiplied by SCALE.
*
*  M       (input) INTEGER
*          The number of rows of the matrix A.
*
*  N       (input) INTEGER
*          The number of columns of the matrix A.
*
*  RCOND   (input) REAL
*          1/RCOND is the condition number of the matrix to be
*          generated. Singular values for the submatrix consisting
*          of independent columns are generated between
*          1 and RCOND dependent on MODE.
*
*  WIDTH   (input) INTEGER
*          The width of a strip of dependent or independent columns.
*
*  MODE    (input) INTEGER
*          is passed to SLATMS to determine how diagonal entries
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
*  S       (output) REAL array (min(M,N))
*          The singular values of A
*
*  A       (output) REAL array, dimension (M,N)
*          matrix with singular value distribution given in S
*          and pattern of dependent/independent columns determined
*          by WIDTH.
*
*  LDA     (input) INTEGER
*          leading dimension of A.
*
*  WORK    (workspace) REAL array,
*          dimension max(3*min(m,n),width*width) if OPT == 'L' or 'l'
*          dimension max(3*min(m,n),width*width*2) otherwise
*
*  =====================================================================
*
*     .. Parameters ..
      REAL               ZERO, ONE
      PARAMETER          ( ZERO = 0.0E+0, ONE = 1.0E+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            MN, WLAST, NSTRPS, INFO, OFFSET,
     $                   NCOLS, CLSLFT, I
      REAL               DUMMY
*     ..
*     ..
*     .. External Subroutines
      EXTERNAL           SLATMS, SLACPY, LSAME,
     $                   SGEBD2, SBDSQR, SLARNV
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
*     How many strips is the matrix of size m -by- rank partitioned into?
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
      CALL SLATMS( M, RANK, 'Uniform', ISEED, 'Nonsymmetric',
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
         CALL SLACPY( 'full matrix', M, WIDTH,
     $               A( 1, OFFSET+( I-1 )*WIDTH+1 ), LDA,
     $               A( 1, 2*( I-1 )*WIDTH+1 ),LDA )
         IF( LSAME( OPT, 'L' ) ) THEN
             NCOLS = WIDTH
         ELSE
             NCOLS = MIN( 2, I )*WIDTH
         END IF
         CALL SLARNV( 1, ISEED, NCOLS*WIDTH, WORK( 1 ) )
         CALL SGEMM( 'no transpose', 'no transpose', M, WIDTH, NCOLS,
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
         CALL SLARNV( 1, ISEED, NCOLS*CLSLFT, WORK( 1 ) )
         CALL SGEMM( 'no transpose', 'no transpose', M, CLSLFT, NCOLS,
     $              SCALE, A( 1, OFFSET+RANK+1-NCOLS ), LDA, WORK,
     $              NCOLS, ZERO,A( 1, OFFSET+RANK+1 ), LDA )
      END IF
*
*     compute singular value decomposition of A
*
      CALL SLACPY( 'full matrix', M, N, A, LDA, WORK( MN+1 ), M )
      CALL SGEBD2( M, N, WORK( MN+1 ), M, S, WORK( 1 ),
     $            WORK( MN+M*N+1 ), WORK( 2*MN+M*N+1 ),
     $            WORK( 3*MN+M*N+1 ), INFO )
      CALL SBDSQR( 'upper', MN, 0, 0, 0, S, WORK( 1 ),
     $            DUMMY, MN, DUMMY, 1, DUMMY, MN, WORK( MN+1 ), INFO )
      RETURN
999   FORMAT( '** ERROR in sqrmtx: SLATMS returns INFO = ',i2 )
*
*     End of SQRMTX
*
      END
