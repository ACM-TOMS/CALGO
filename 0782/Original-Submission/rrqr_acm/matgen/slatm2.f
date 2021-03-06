      REAL             FUNCTION SLATM2( M, N, I, J, KL, KU, IDIST,
     $                 ISEED, D, IGRADE, DL, DR, IPVTNG, IWORK, SPARSE )
*
*  -- LAPACK auxiliary test routine (version 2.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     February 29, 1992
*
*     .. Scalar Arguments ..
*
      INTEGER            I, IDIST, IGRADE, IPVTNG, J, KL, KU, M, N
      REAL               SPARSE
*     ..
*
*     .. Array Arguments ..
*
      INTEGER            ISEED( 4 ), IWORK( * )
      REAL               D( * ), DL( * ), DR( * )
*     ..
*
*  Purpose
*  =======
*
*     SLATM2 returns the (I,J) entry of a random matrix of dimension
*     (M, N) described by the other paramters. It is called by the
*     SLATMR routine in order to build random test matrices. No error
*     checking on parameters is done, because this routine is called in
*     a tight loop by SLATMR which has already checked the parameters.
*
*     Use of SLATM2 differs from SLATM3 in the order in which the random
*     number generator is called to fill in random matrix entries.
*     With SLATM2, the generator is called to fill in the pivoted matrix
*     columnwise. With SLATM3, the generator is called to fill in the
*     matrix columnwise, after which it is pivoted. Thus, SLATM3 can
*     be used to construct random matrices which differ only in their
*     order of rows and/or columns. SLATM2 is used to construct band
*     matrices while avoiding calling the random number generator for
*     entries outside the band (and therefore generating random numbers
*
*     The matrix whose (I,J) entry is returned is constructed as
*     follows (this routine only computes one entry):
*
*       If I is outside (1..M) or J is outside (1..N), return zero
*          (this is convenient for generating matrices in band format).
*
*       Generate a matrix A with random entries of distribution IDIST.
*
*       Set the diagonal to D.
*
*       Grade the matrix, if desired, from the left (by DL) and/or
*          from the right (by DR or DL) as specified by IGRADE.
*
*       Permute, if desired, the rows and/or columns as specified by
*          IPVTNG and IWORK.
*
*       Band the matrix to have lower bandwidth KL and upper
*          bandwidth KU.
*
*       Set random entries to zero as specified by SPARSE.
*
*  Arguments
*  =========
*
*  M      - INTEGER
*           Number of rows of matrix. Not modified.
*
*  N      - INTEGER
*           Number of columns of matrix. Not modified.
*
*  I      - INTEGER
*           Row of entry to be returned. Not modified.
*
*  J      - INTEGER
*           Column of entry to be returned. Not modified.
*
*  KL     - INTEGER
*           Lower bandwidth. Not modified.
*
*  KU     - INTEGER
*           Upper bandwidth. Not modified.
*
*  IDIST  - INTEGER
*           On entry, IDIST specifies the type of distribution to be
*           used to generate a random matrix .
*           1 => UNIFORM( 0, 1 )
*           2 => UNIFORM( -1, 1 )
*           3 => NORMAL( 0, 1 )
*           Not modified.
*
*  ISEED  - INTEGER array of dimension ( 4 )
*           Seed for random number generator.
*           Changed on exit.
*
*  D      - REAL array of dimension ( MIN( I , J ) )
*           Diagonal entries of matrix. Not modified.
*
*  IGRADE - INTEGER
*           Specifies grading of matrix as follows:
*           0  => no grading
*           1  => matrix premultiplied by diag( DL )
*           2  => matrix postmultiplied by diag( DR )
*           3  => matrix premultiplied by diag( DL ) and
*                         postmultiplied by diag( DR )
*           4  => matrix premultiplied by diag( DL ) and
*                         postmultiplied by inv( diag( DL ) )
*           5  => matrix premultiplied by diag( DL ) and
*                         postmultiplied by diag( DL )
*           Not modified.
*
*  DL     - REAL array ( I or J, as appropriate )
*           Left scale factors for grading matrix.  Not modified.
*
*  DR     - REAL array ( I or J, as appropriate )
*           Right scale factors for grading matrix.  Not modified.
*
*  IPVTNG - INTEGER
*           On entry specifies pivoting permutations as follows:
*           0 => none.
*           1 => row pivoting.
*           2 => column pivoting.
*           3 => full pivoting, i.e., on both sides.
*           Not modified.
*
*  IWORK  - INTEGER array ( I or J, as appropriate )
*           This array specifies the permutation used. The
*           row (or column) in position K was originally in
*           position IWORK( K ).
*           This differs from IWORK for SLATM3. Not modified.
*
*  SPARSE - REAL    between 0. and 1.
*           On entry specifies the sparsity of the matrix
*           if sparse matix is to be generated.
*           SPARSE should lie between 0 and 1.
*           A uniform ( 0, 1 ) random number x is generated and
*           compared to SPARSE; if x is larger the matrix entry
*           is unchanged and if x is smaller the entry is set
*           to zero. Thus on the average a fraction SPARSE of the
*           entries will be set to zero.
*           Not modified.
*
*  =====================================================================
*
*     .. Parameters ..
*
      REAL               ZERO
      PARAMETER          ( ZERO = 0.0E0 )
*     ..
*
*     .. Local Scalars ..
*
      INTEGER            ISUB, JSUB
      REAL               TEMP
*     ..
*
*     .. External Functions ..
*
      REAL               SLARAN, SLARND
      EXTERNAL           SLARAN, SLARND
*     ..
*
*-----------------------------------------------------------------------
*
*     .. Executable Statements ..
*
*
*     Check for I and J in range
*
      IF( I.LT.1 .OR. I.GT.M .OR. J.LT.1 .OR. J.GT.N ) THEN
         SLATM2 = ZERO
         RETURN
      END IF
*
*     Check for banding
*
      IF( J.GT.I+KU .OR. J.LT.I-KL ) THEN
         SLATM2 = ZERO
         RETURN
      END IF
*
*     Check for sparsity
*
      IF( SPARSE.GT.ZERO ) THEN
         IF( SLARAN( ISEED ).LT.SPARSE ) THEN
            SLATM2 = ZERO
            RETURN
         END IF
      END IF
*
*     Compute subscripts depending on IPVTNG
*
      IF( IPVTNG.EQ.0 ) THEN
         ISUB = I
         JSUB = J
      ELSE IF( IPVTNG.EQ.1 ) THEN
         ISUB = IWORK( I )
         JSUB = J
      ELSE IF( IPVTNG.EQ.2 ) THEN
         ISUB = I
         JSUB = IWORK( J )
      ELSE IF( IPVTNG.EQ.3 ) THEN
         ISUB = IWORK( I )
         JSUB = IWORK( J )
      END IF
*
*     Compute entry and grade it according to IGRADE
*
      IF( ISUB.EQ.JSUB ) THEN
         TEMP = D( ISUB )
      ELSE
         TEMP = SLARND( IDIST, ISEED )
      END IF
      IF( IGRADE.EQ.1 ) THEN
         TEMP = TEMP*DL( ISUB )
      ELSE IF( IGRADE.EQ.2 ) THEN
         TEMP = TEMP*DR( JSUB )
      ELSE IF( IGRADE.EQ.3 ) THEN
         TEMP = TEMP*DL( ISUB )*DR( JSUB )
      ELSE IF( IGRADE.EQ.4 .AND. ISUB.NE.JSUB ) THEN
         TEMP = TEMP*DL( ISUB ) / DL( JSUB )
      ELSE IF( IGRADE.EQ.5 ) THEN
         TEMP = TEMP*DL( ISUB )*DL( JSUB )
      END IF
      SLATM2 = TEMP
      RETURN
*
*     End of SLATM2
*
      END
