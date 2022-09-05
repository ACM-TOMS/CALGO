* =====================================================================
      SUBROUTINE DLATMS( M, N, DIST, ISEED, SYM, D, MODE, COND, DMAX,
     $                   KL, KU, PACK, A, LDA, WORK, INFO )
*
*  -- LAPACK test routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     September 30, 1994
*
*     .. Scalar Arguments ..
      CHARACTER          DIST, PACK, SYM
      INTEGER            INFO, KL, KU, LDA, M, MODE, N
      DOUBLE PRECISION   COND, DMAX
*     ..
*     .. Array Arguments ..
      INTEGER            ISEED( 4 )
      DOUBLE PRECISION   A( LDA, * ), D( * ), WORK( * )
*     ..
*
*  Purpose
*  =======
*
*     DLATMS generates random matrices with specified singular values
*     (or symmetric/hermitian with specified eigenvalues)
*     for testing LAPACK programs.
*
*     DLATMS operates by applying the following sequence of
*     operations:
*
*       Set the diagonal to D, where D may be input or
*          computed according to MODE, COND, DMAX, and SYM
*          as described below.
*
*       Generate a matrix with the appropriate band structure, by one
*          of two methods:
*
*       Method A:
*           Generate a dense M x N matrix by multiplying D on the left
*               and the right by random unitary matrices, then:
*
*           Reduce the bandwidth according to KL and KU, using
*           Householder transformations.
*
*       Method B:
*           Convert the bandwidth-0 (i.e., diagonal) matrix to a
*               bandwidth-1 matrix using Givens rotations, "chasing"
*               out-of-band elements back, much as in QR; then
*               convert the bandwidth-1 to a bandwidth-2 matrix, etc.
*               Note that for reasonably small bandwidths (relative to
*               M and N) this requires less storage, as a dense matrix
*               is not generated.  Also, for symmetric matrices, only
*               one triangle is generated.
*
*       Method A is chosen if the bandwidth is a large fraction of the
*           order of the matrix, and LDA is at least M (so a dense
*           matrix can be stored.)  Method B is chosen if the bandwidth
*           is small (< 1/2 N for symmetric, < .3 N+M for
*           non-symmetric), or LDA is less than M and not less than the
*           bandwidth.
*
*       Pack the matrix if desired. Options specified by PACK are:
*          no packing
*          zero out upper half (if symmetric)
*          zero out lower half (if symmetric)
*          store the upper half columnwise (if symmetric or upper
*                triangular)
*          store the lower half columnwise (if symmetric or lower
*                triangular)
*          store the lower triangle in banded format (if symmetric
*                or lower triangular)
*          store the upper triangle in banded format (if symmetric
*                or upper triangular)
*          store the entire matrix in banded format
*       If Method B is chosen, and band format is specified, then the
*          matrix will be generated in the band format, so no repacking
*          will be necessary.
*
*  Arguments
*  =========
*
*  M      - INTEGER
*           The number of rows of A. Not modified.
*
*  N      - INTEGER
*           The number of columns of A. Not modified.
*
*  DIST   - CHARACTER*1
*           On entry, DIST specifies the type of distribution to be used
*           to generate the random eigen-/singular values.
*           'U' => UNIFORM( 0, 1 )  ( 'U' for uniform )
*           'S' => UNIFORM( -1, 1 ) ( 'S' for symmetric )
*           'N' => NORMAL( 0, 1 )   ( 'N' for normal )
*           Not modified.
*
*  ISEED  - INTEGER array, dimension ( 4 )
*           On entry ISEED specifies the seed of the random number
*           generator. They should lie between 0 and 4095 inclusive,
*           and ISEED(4) should be odd. The random number generator
*           uses a linear congruential sequence limited to small
*           integers, and so should produce machine independent
*           random numbers. The values of ISEED are changed on
*           exit, and can be used in the next call to DLATMS
*           to continue the same random number sequence.
*           Changed on exit.
*
*  SYM    - CHARACTER*1
*           If SYM='S' or 'H', the generated matrix is symmetric, with
*             eigenvalues specified by D, COND, MODE, and DMAX; they
*             may be positive, negative, or zero.
*           If SYM='P', the generated matrix is symmetric, with
*             eigenvalues (= singular values) specified by D, COND,
*             MODE, and DMAX; they will not be negative.
*           If SYM='N', the generated matrix is nonsymmetric, with
*             singular values specified by D, COND, MODE, and DMAX;
*             they will not be negative.
*           Not modified.
*
*  D      - DOUBLE PRECISION array, dimension ( MIN( M , N ) )
*           This array is used to specify the singular values or
*           eigenvalues of A (see SYM, above.)  If MODE=0, then D is
*           assumed to contain the singular/eigenvalues, otherwise
*           they will be computed according to MODE, COND, and DMAX,
*           and placed in D.
*           Modified if MODE is nonzero.
*
*  MODE   - INTEGER
*           On entry this describes how the singular/eigenvalues are to
*           be specified:
*           MODE = 0 means use D as input
*           MODE = 1 sets D(1)=1 and D(2:N)=1.0/COND
*           MODE = 2 sets D(1:N-1)=1 and D(N)=1.0/COND
*           MODE = 3 sets D(I)=COND**(-(I-1)/(N-1))
*           MODE = 4 sets D(i)=1 - (i-1)/(N-1)*(1 - 1/COND)
*           MODE = 5 sets D to random numbers in the range
*                    ( 1/COND , 1 ) such that their logarithms
*                    are uniformly distributed.
*           MODE = 6 set D to random numbers from same distribution
*                    as the rest of the matrix.
*           MODE < 0 has the same meaning as ABS(MODE), except that
*              the order of the elements of D is reversed.
*           Thus if MODE is positive, D has entries ranging from
*              1 to 1/COND, if negative, from 1/COND to 1,
*           If SYM='S' or 'H', and MODE is neither 0, 6, nor -6, then
*              the elements of D will also be multiplied by a random
*              sign (i.e., +1 or -1.)
*           Not modified.
*
*  COND   - DOUBLE PRECISION
*           On entry, this is used as described under MODE above.
*           If used, it must be >= 1. Not modified.
*
*  DMAX   - DOUBLE PRECISION
*           If MODE is neither -6, 0 nor 6, the contents of D, as
*           computed according to MODE and COND, will be scaled by
*           DMAX / max(abs(D(i))); thus, the maximum absolute eigen- or
*           singular value (which is to say the norm) will be abs(DMAX).
*           Note that DMAX need not be positive: if DMAX is negative
*           (or zero), D will be scaled by a negative number (or zero).
*           Not modified.
*
*  KL     - INTEGER
*           This specifies the lower bandwidth of the  matrix. For
*           example, KL=0 implies upper triangular, KL=1 implies upper
*           Hessenberg, and KL being at least M-1 means that the matrix
*           has full lower bandwidth.  KL must equal KU if the matrix
*           is symmetric.
*           Not modified.
*
*  KU     - INTEGER
*           This specifies the upper bandwidth of the  matrix. For
*           example, KU=0 implies lower triangular, KU=1 implies lower
*           Hessenberg, and KU being at least N-1 means that the matrix
*           has full upper bandwidth.  KL must equal KU if the matrix
*           is symmetric.
*           Not modified.
*
*  PACK   - CHARACTER*1
*           This specifies packing of matrix as follows:
*           'N' => no packing
*           'U' => zero out all subdiagonal entries (if symmetric)
*           'L' => zero out all superdiagonal entries (if symmetric)
*           'C' => store the upper triangle columnwise
*                  (only if the matrix is symmetric or upper triangular)
*           'R' => store the lower triangle columnwise
*                  (only if the matrix is symmetric or lower triangular)
*           'B' => store the lower triangle in band storage scheme
*                  (only if matrix symmetric or lower triangular)
*           'Q' => store the upper triangle in band storage scheme
*                  (only if matrix symmetric or upper triangular)
*           'Z' => store the entire matrix in band storage scheme
*                      (pivoting can be provided for by using this
*                      option to store A in the trailing rows of
*                      the allocated storage)
*
*           Using these options, the various LAPACK packed and banded
*           storage schemes can be obtained:
*           GB               - use 'Z'
*           PB, SB or TB     - use 'B' or 'Q'
*           PP, SP or TP     - use 'C' or 'R'
*
*           If two calls to DLATMS differ only in the PACK parameter,
*           they will generate mathematically equivalent matrices.
*           Not modified.
*
*  A      - DOUBLE PRECISION array, dimension ( LDA, N )
*           On exit A is the desired test matrix.  A is first generated
*           in full (unpacked) form, and then packed, if so specified
*           by PACK.  Thus, the first M elements of the first N
*           columns will always be modified.  If PACK specifies a
*           packed or banded storage scheme, all LDA elements of the
*           first N columns will be modified; the elements of the
*           array which do not correspond to elements of the generated
*           matrix are set to zero.
*           Modified.
*
*  LDA    - INTEGER
*           LDA specifies the first dimension of A as declared in the
*           calling program.  If PACK='N', 'U', 'L', 'C', or 'R', then
*           LDA must be at least M.  If PACK='B' or 'Q', then LDA must
*           be at least MIN( KL, M-1) (which is equal to MIN(KU,N-1)).
*           If PACK='Z', LDA must be large enough to hold the packed
*           array: MIN( KU, N-1) + MIN( KL, M-1) + 1.
*           Not modified.
*
*  WORK   - DOUBLE PRECISION array, dimension ( 3*MAX( N , M ) )
*           Workspace.
*           Modified.
*
*  INFO   - INTEGER
*           Error code.  On exit, INFO will be set to one of the
*           following values:
*             0 => normal return
*            -1 => M negative or unequal to N and SYM='S', 'H', or 'P'
*            -2 => N negative
*            -3 => DIST illegal string
*            -5 => SYM illegal string
*            -7 => MODE not in range -6 to 6
*            -8 => COND less than 1.0, and MODE neither -6, 0 nor 6
*           -10 => KL negative
*           -11 => KU negative, or SYM='S' or 'H' and KU not equal to KL
*           -12 => PACK illegal string, or PACK='U' or 'L', and SYM='N';
*                  or PACK='C' or 'Q' and SYM='N' and KL is not zero;
*                  or PACK='R' or 'B' and SYM='N' and KU is not zero;
*                  or PACK='U', 'L', 'C', 'R', 'B', or 'Q', and M is not
*                  N.
*           -14 => LDA is less than M, or PACK='Z' and LDA is less than
*                  MIN(KU,N-1) + MIN(KL,M-1) + 1.
*            1  => Error return from DLATM1
*            2  => Cannot scale to DMAX (max. sing. value is 0)
*            3  => Error return from DLAGGE or SLAGSY
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO
      PARAMETER          ( ZERO = 0.0D0 )
      DOUBLE PRECISION   ONE
      PARAMETER          ( ONE = 1.0D0 )
      DOUBLE PRECISION   TWOPI
      PARAMETER          ( TWOPI = 6.2831853071795864769252867663D+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            GIVENS, ILEXTR, ILTEMP, TOPDWN
      INTEGER            I, IC, ICOL, IDIST, IENDCH, IINFO, IL, ILDA,
     $                   IOFFG, IOFFST, IPACK, IPACKG, IR, IR1, IR2,
     $                   IROW, IRSIGN, ISKEW, ISYM, ISYMPK, J, JC, JCH,
     $                   JKL, JKU, JR, K, LLB, MINLDA, MNMIN, MR, NC,
     $                   UUB
      DOUBLE PRECISION   ALPHA, ANGLE, C, DUMMY, EXTRA, S, TEMP
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      DOUBLE PRECISION   DLARND
      EXTERNAL           LSAME, DLARND
*     ..
*     .. External Subroutines ..
      EXTERNAL           DCOPY, DLAGGE, DLAGSY, DLAROT, DLARTG, DLASET,
     $                   DLATM1, DSCAL, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, COS, DBLE, MAX, MIN, MOD, SIN
*     ..
*     .. Executable Statements ..
*
*     1)      Decode and Test the input parameters.
*             Initialize flags & seed.
*
      INFO = 0
*
*     Quick return if possible
*
      IF( M.EQ.0 .OR. N.EQ.0 )
     $   RETURN
*
*     Decode DIST
*
      IF( LSAME( DIST, 'U' ) ) THEN
         IDIST = 1
      ELSE IF( LSAME( DIST, 'S' ) ) THEN
         IDIST = 2
      ELSE IF( LSAME( DIST, 'N' ) ) THEN
         IDIST = 3
      ELSE
         IDIST = -1
      END IF
*
*     Decode SYM
*
      IF( LSAME( SYM, 'N' ) ) THEN
         ISYM = 1
         IRSIGN = 0
      ELSE IF( LSAME( SYM, 'P' ) ) THEN
         ISYM = 2
         IRSIGN = 0
      ELSE IF( LSAME( SYM, 'S' ) ) THEN
         ISYM = 2
         IRSIGN = 1
      ELSE IF( LSAME( SYM, 'H' ) ) THEN
         ISYM = 2
         IRSIGN = 1
      ELSE
         ISYM = -1
      END IF
*
*     Decode PACK
*
      ISYMPK = 0
      IF( LSAME( PACK, 'N' ) ) THEN
         IPACK = 0
      ELSE IF( LSAME( PACK, 'U' ) ) THEN
         IPACK = 1
         ISYMPK = 1
      ELSE IF( LSAME( PACK, 'L' ) ) THEN
         IPACK = 2
         ISYMPK = 1
      ELSE IF( LSAME( PACK, 'C' ) ) THEN
         IPACK = 3
         ISYMPK = 2
      ELSE IF( LSAME( PACK, 'R' ) ) THEN
         IPACK = 4
         ISYMPK = 3
      ELSE IF( LSAME( PACK, 'B' ) ) THEN
         IPACK = 5
         ISYMPK = 3
      ELSE IF( LSAME( PACK, 'Q' ) ) THEN
         IPACK = 6
         ISYMPK = 2
      ELSE IF( LSAME( PACK, 'Z' ) ) THEN
         IPACK = 7
      ELSE
         IPACK = -1
      END IF
*
*     Set certain internal parameters
*
      MNMIN = MIN( M, N )
      LLB = MIN( KL, M-1 )
      UUB = MIN( KU, N-1 )
      MR = MIN( M, N+LLB )
      NC = MIN( N, M+UUB )
*
      IF( IPACK.EQ.5 .OR. IPACK.EQ.6 ) THEN
         MINLDA = UUB + 1
      ELSE IF( IPACK.EQ.7 ) THEN
         MINLDA = LLB + UUB + 1
      ELSE
         MINLDA = M
      END IF
*
*     Use Givens rotation method if bandwidth small enough,
*     or if LDA is too small to store the matrix unpacked.
*
      GIVENS = .FALSE.
      IF( ISYM.EQ.1 ) THEN
         IF( DBLE( LLB+UUB ).LT.0.3D0*DBLE( MAX( 1, MR+NC ) ) )
     $      GIVENS = .TRUE.
      ELSE
         IF( 2*LLB.LT.M )
     $      GIVENS = .TRUE.
      END IF
      IF( LDA.LT.M .AND. LDA.GE.MINLDA )
     $   GIVENS = .TRUE.
*
*     Set INFO if an error
*
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( M.NE.N .AND. ISYM.NE.1 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( IDIST.EQ.-1 ) THEN
         INFO = -3
      ELSE IF( ISYM.EQ.-1 ) THEN
         INFO = -5
      ELSE IF( ABS( MODE ).GT.6 ) THEN
         INFO = -7
      ELSE IF( ( MODE.NE.0 .AND. ABS( MODE ).NE.6 ) .AND. COND.LT.ONE )
     $          THEN
         INFO = -8
      ELSE IF( KL.LT.0 ) THEN
         INFO = -10
      ELSE IF( KU.LT.0 .OR. ( ISYM.NE.1 .AND. KL.NE.KU ) ) THEN
         INFO = -11
      ELSE IF( IPACK.EQ.-1 .OR. ( ISYMPK.EQ.1 .AND. ISYM.EQ.1 ) .OR.
     $         ( ISYMPK.EQ.2 .AND. ISYM.EQ.1 .AND. KL.GT.0 ) .OR.
     $         ( ISYMPK.EQ.3 .AND. ISYM.EQ.1 .AND. KU.GT.0 ) .OR.
     $         ( ISYMPK.NE.0 .AND. M.NE.N ) ) THEN
         INFO = -12
      ELSE IF( LDA.LT.MAX( 1, MINLDA ) ) THEN
         INFO = -14
      END IF
*
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DLATMS', -INFO )
         RETURN
      END IF
*
*     Initialize random number generator
*
      DO 10 I = 1, 4
         ISEED( I ) = MOD( ABS( ISEED( I ) ), 4096 )
   10 CONTINUE
*
      IF( MOD( ISEED( 4 ), 2 ).NE.1 )
     $   ISEED( 4 ) = ISEED( 4 ) + 1
*
*     2)      Set up D  if indicated.
*
*             Compute D according to COND and MODE
*
      CALL DLATM1( MODE, COND, IRSIGN, IDIST, ISEED, D, MNMIN, IINFO )
      IF( IINFO.NE.0 ) THEN
         INFO = 1
         RETURN
      END IF
*
*     Choose Top-Down if D is (apparently) increasing,
*     Bottom-Up if D is (apparently) decreasing.
*
      IF( ABS( D( 1 ) ).LE.ABS( D( MNMIN ) ) ) THEN
         TOPDWN = .TRUE.
      ELSE
         TOPDWN = .FALSE.
      END IF
*
      IF( MODE.NE.0 .AND. ABS( MODE ).NE.6 ) THEN
*
*        Scale by DMAX
*
         TEMP = ABS( D( 1 ) )
         DO 20 I = 2, MNMIN
            TEMP = MAX( TEMP, ABS( D( I ) ) )
   20    CONTINUE
*
         IF( TEMP.GT.ZERO ) THEN
            ALPHA = DMAX / TEMP
         ELSE
            INFO = 2
            RETURN
         END IF
*
         CALL DSCAL( MNMIN, ALPHA, D, 1 )
*
      END IF
*
*     3)      Generate Banded Matrix using Givens rotations.
*             Also the special case of UUB=LLB=0
*
*               Compute Addressing constants to cover all
*               storage formats.  Whether GE, SY, GB, or SB,
*               upper or lower triangle or both,
*               the (i,j)-th element is in
*               A( i - ISKEW*j + IOFFST, j )
*
      IF( IPACK.GT.4 ) THEN
         ILDA = LDA - 1
         ISKEW = 1
         IF( IPACK.GT.5 ) THEN
            IOFFST = UUB + 1
         ELSE
            IOFFST = 1
         END IF
      ELSE
         ILDA = LDA
         ISKEW = 0
         IOFFST = 0
      END IF
*
*     IPACKG is the format that the matrix is generated in. If this is
*     different from IPACK, then the matrix must be repacked at the
*     end.  It also signals how to compute the norm, for scaling.
*
      IPACKG = 0
      CALL DLASET( 'Full', LDA, N, ZERO, ZERO, A, LDA )
*
*     Diagonal Matrix -- We are done, unless it
*     is to be stored SP/PP/TP (PACK='R' or 'C')
*
      IF( LLB.EQ.0 .AND. UUB.EQ.0 ) THEN
         CALL DCOPY( MNMIN, D, 1, A( 1-ISKEW+IOFFST, 1 ), ILDA+1 )
         IF( IPACK.LE.2 .OR. IPACK.GE.5 )
     $      IPACKG = IPACK
*
      ELSE IF( GIVENS ) THEN
*
*        Check whether to use Givens rotations,
*        Householder transformations, or nothing.
*
         IF( ISYM.EQ.1 ) THEN
*
*           Non-symmetric -- A = U D V
*
            IF( IPACK.GT.4 ) THEN
               IPACKG = IPACK
            ELSE
               IPACKG = 0
            END IF
*
            CALL DCOPY( MNMIN, D, 1, A( 1-ISKEW+IOFFST, 1 ), ILDA+1 )
*
            IF( TOPDWN ) THEN
               JKL = 0
               DO 50 JKU = 1, UUB
*
*                 Transform from bandwidth JKL, JKU-1 to JKL, JKU
*
*                 Last row actually rotated is M
*                 Last column actually rotated is MIN( M+JKU, N )
*
                  DO 40 JR = 1, MIN( M+JKU, N ) + JKL - 1
                     EXTRA = ZERO
                     ANGLE = TWOPI*DLARND( 1, ISEED )
                     C = COS( ANGLE )
                     S = SIN( ANGLE )
                     ICOL = MAX( 1, JR-JKL )
                     IF( JR.LT.M ) THEN
                        IL = MIN( N, JR+JKU ) + 1 - ICOL
                        CALL DLAROT( .TRUE., JR.GT.JKL, .FALSE., IL, C,
     $                               S, A( JR-ISKEW*ICOL+IOFFST, ICOL ),
     $                               ILDA, EXTRA, DUMMY )
                     END IF
*
*                    Chase "EXTRA" back up
*
                     IR = JR
                     IC = ICOL
                     DO 30 JCH = JR - JKL, 1, -JKL - JKU
                        IF( IR.LT.M ) THEN
                           CALL DLARTG( A( IR+1-ISKEW*( IC+1 )+IOFFST,
     $                                  IC+1 ), EXTRA, C, S, DUMMY )
                        END IF
                        IROW = MAX( 1, JCH-JKU )
                        IL = IR + 2 - IROW
                        TEMP = ZERO
                        ILTEMP = JCH.GT.JKU
                        CALL DLAROT( .FALSE., ILTEMP, .TRUE., IL, C, -S,
     $                               A( IROW-ISKEW*IC+IOFFST, IC ),
     $                               ILDA, TEMP, EXTRA )
                        IF( ILTEMP ) THEN
                           CALL DLARTG( A( IROW+1-ISKEW*( IC+1 )+IOFFST,
     $                                  IC+1 ), TEMP, C, S, DUMMY )
                           ICOL = MAX( 1, JCH-JKU-JKL )
                           IL = IC + 2 - ICOL
                           EXTRA = ZERO
                           CALL DLAROT( .TRUE., JCH.GT.JKU+JKL, .TRUE.,
     $                                  IL, C, -S, A( IROW-ISKEW*ICOL+
     $                                  IOFFST, ICOL ), ILDA, EXTRA,
     $                                  TEMP )
                           IC = ICOL
                           IR = IROW
                        END IF
   30                CONTINUE
   40             CONTINUE
   50          CONTINUE
*
               JKU = UUB
               DO 80 JKL = 1, LLB
*
*                 Transform from bandwidth JKL-1, JKU to JKL, JKU
*
                  DO 70 JC = 1, MIN( N+JKL, M ) + JKU - 1
                     EXTRA = ZERO
                     ANGLE = TWOPI*DLARND( 1, ISEED )
                     C = COS( ANGLE )
                     S = SIN( ANGLE )
                     IROW = MAX( 1, JC-JKU )
                     IF( JC.LT.N ) THEN
                        IL = MIN( M, JC+JKL ) + 1 - IROW
                        CALL DLAROT( .FALSE., JC.GT.JKU, .FALSE., IL, C,
     $                               S, A( IROW-ISKEW*JC+IOFFST, JC ),
     $                               ILDA, EXTRA, DUMMY )
                     END IF
*
*                    Chase "EXTRA" back up
*
                     IC = JC
                     IR = IROW
                     DO 60 JCH = JC - JKU, 1, -JKL - JKU
                        IF( IC.LT.N ) THEN
                           CALL DLARTG( A( IR+1-ISKEW*( IC+1 )+IOFFST,
     $                                  IC+1 ), EXTRA, C, S, DUMMY )
                        END IF
                        ICOL = MAX( 1, JCH-JKL )
                        IL = IC + 2 - ICOL
                        TEMP = ZERO
                        ILTEMP = JCH.GT.JKL
                        CALL DLAROT( .TRUE., ILTEMP, .TRUE., IL, C, -S,
     $                               A( IR-ISKEW*ICOL+IOFFST, ICOL ),
     $                               ILDA, TEMP, EXTRA )
                        IF( ILTEMP ) THEN
                           CALL DLARTG( A( IR+1-ISKEW*( ICOL+1 )+IOFFST,
     $                                  ICOL+1 ), TEMP, C, S, DUMMY )
                           IROW = MAX( 1, JCH-JKL-JKU )
                           IL = IR + 2 - IROW
                           EXTRA = ZERO
                           CALL DLAROT( .FALSE., JCH.GT.JKL+JKU, .TRUE.,
     $                                  IL, C, -S, A( IROW-ISKEW*ICOL+
     $                                  IOFFST, ICOL ), ILDA, EXTRA,
     $                                  TEMP )
                           IC = ICOL
                           IR = IROW
                        END IF
   60                CONTINUE
   70             CONTINUE
   80          CONTINUE
*
            ELSE
*
*              Bottom-Up -- Start at the bottom right.
*
               JKL = 0
               DO 110 JKU = 1, UUB
*
*                 Transform from bandwidth JKL, JKU-1 to JKL, JKU
*
*                 First row actually rotated is M
*                 First column actually rotated is MIN( M+JKU, N )
*
                  IENDCH = MIN( M, N+JKL ) - 1
                  DO 100 JC = MIN( M+JKU, N ) - 1, 1 - JKL, -1
                     EXTRA = ZERO
                     ANGLE = TWOPI*DLARND( 1, ISEED )
                     C = COS( ANGLE )
                     S = SIN( ANGLE )
                     IROW = MAX( 1, JC-JKU+1 )
                     IF( JC.GT.0 ) THEN
                        IL = MIN( M, JC+JKL+1 ) + 1 - IROW
                        CALL DLAROT( .FALSE., .FALSE., JC+JKL.LT.M, IL,
     $                               C, S, A( IROW-ISKEW*JC+IOFFST,
     $                               JC ), ILDA, DUMMY, EXTRA )
                     END IF
*
*                    Chase "EXTRA" back down
*
                     IC = JC
                     DO 90 JCH = JC + JKL, IENDCH, JKL + JKU
                        ILEXTR = IC.GT.0
                        IF( ILEXTR ) THEN
                           CALL DLARTG( A( JCH-ISKEW*IC+IOFFST, IC ),
     $                                  EXTRA, C, S, DUMMY )
                        END IF
                        IC = MAX( 1, IC )
                        ICOL = MIN( N-1, JCH+JKU )
                        ILTEMP = JCH + JKU.LT.N
                        TEMP = ZERO
                        CALL DLAROT( .TRUE., ILEXTR, ILTEMP, ICOL+2-IC,
     $                               C, S, A( JCH-ISKEW*IC+IOFFST, IC ),
     $                               ILDA, EXTRA, TEMP )
                        IF( ILTEMP ) THEN
                           CALL DLARTG( A( JCH-ISKEW*ICOL+IOFFST,
     $                                  ICOL ), TEMP, C, S, DUMMY )
                           IL = MIN( IENDCH, JCH+JKL+JKU ) + 2 - JCH
                           EXTRA = ZERO
                           CALL DLAROT( .FALSE., .TRUE.,
     $                                  JCH+JKL+JKU.LE.IENDCH, IL, C, S,
     $                                  A( JCH-ISKEW*ICOL+IOFFST,
     $                                  ICOL ), ILDA, TEMP, EXTRA )
                           IC = ICOL
                        END IF
   90                CONTINUE
  100             CONTINUE
  110          CONTINUE
*
               JKU = UUB
               DO 140 JKL = 1, LLB
*
*                 Transform from bandwidth JKL-1, JKU to JKL, JKU
*
*                 First row actually rotated is MIN( N+JKL, M )
*                 First column actually rotated is N
*
                  IENDCH = MIN( N, M+JKU ) - 1
                  DO 130 JR = MIN( N+JKL, M ) - 1, 1 - JKU, -1
                     EXTRA = ZERO
                     ANGLE = TWOPI*DLARND( 1, ISEED )
                     C = COS( ANGLE )
                     S = SIN( ANGLE )
                     ICOL = MAX( 1, JR-JKL+1 )
                     IF( JR.GT.0 ) THEN
                        IL = MIN( N, JR+JKU+1 ) + 1 - ICOL
                        CALL DLAROT( .TRUE., .FALSE., JR+JKU.LT.N, IL,
     $                               C, S, A( JR-ISKEW*ICOL+IOFFST,
     $                               ICOL ), ILDA, DUMMY, EXTRA )
                     END IF
*
*                    Chase "EXTRA" back down
*
                     IR = JR
                     DO 120 JCH = JR + JKU, IENDCH, JKL + JKU
                        ILEXTR = IR.GT.0
                        IF( ILEXTR ) THEN
                           CALL DLARTG( A( IR-ISKEW*JCH+IOFFST, JCH ),
     $                                  EXTRA, C, S, DUMMY )
                        END IF
                        IR = MAX( 1, IR )
                        IROW = MIN( M-1, JCH+JKL )
                        ILTEMP = JCH + JKL.LT.M
                        TEMP = ZERO
                        CALL DLAROT( .FALSE., ILEXTR, ILTEMP, IROW+2-IR,
     $                               C, S, A( IR-ISKEW*JCH+IOFFST,
     $                               JCH ), ILDA, EXTRA, TEMP )
                        IF( ILTEMP ) THEN
                           CALL DLARTG( A( IROW-ISKEW*JCH+IOFFST, JCH ),
     $                                  TEMP, C, S, DUMMY )
                           IL = MIN( IENDCH, JCH+JKL+JKU ) + 2 - JCH
                           EXTRA = ZERO
                           CALL DLAROT( .TRUE., .TRUE.,
     $                                  JCH+JKL+JKU.LE.IENDCH, IL, C, S,
     $                                  A( IROW-ISKEW*JCH+IOFFST, JCH ),
     $                                  ILDA, TEMP, EXTRA )
                           IR = IROW
                        END IF
  120                CONTINUE
  130             CONTINUE
  140          CONTINUE
            END IF
*
         ELSE
*
*           Symmetric -- A = U D U'
*
            IPACKG = IPACK
            IOFFG = IOFFST
*
            IF( TOPDWN ) THEN
*
*              Top-Down -- Generate Upper triangle only
*
               IF( IPACK.GE.5 ) THEN
                  IPACKG = 6
                  IOFFG = UUB + 1
               ELSE
                  IPACKG = 1
               END IF
               CALL DCOPY( MNMIN, D, 1, A( 1-ISKEW+IOFFG, 1 ), ILDA+1 )
*
               DO 170 K = 1, UUB
                  DO 160 JC = 1, N - 1
                     IROW = MAX( 1, JC-K )
                     IL = MIN( JC+1, K+2 )
                     EXTRA = ZERO
                     TEMP = A( JC-ISKEW*( JC+1 )+IOFFG, JC+1 )
                     ANGLE = TWOPI*DLARND( 1, ISEED )
                     C = COS( ANGLE )
                     S = SIN( ANGLE )
                     CALL DLAROT( .FALSE., JC.GT.K, .TRUE., IL, C, S,
     $                            A( IROW-ISKEW*JC+IOFFG, JC ), ILDA,
     $                            EXTRA, TEMP )
                     CALL DLAROT( .TRUE., .TRUE., .FALSE.,
     $                            MIN( K, N-JC )+1, C, S,
     $                            A( ( 1-ISKEW )*JC+IOFFG, JC ), ILDA,
     $                            TEMP, DUMMY )
*
*                    Chase EXTRA back up the matrix
*
                     ICOL = JC
                     DO 150 JCH = JC - K, 1, -K
                        CALL DLARTG( A( JCH+1-ISKEW*( ICOL+1 )+IOFFG,
     $                               ICOL+1 ), EXTRA, C, S, DUMMY )
                        TEMP = A( JCH-ISKEW*( JCH+1 )+IOFFG, JCH+1 )
                        CALL DLAROT( .TRUE., .TRUE., .TRUE., K+2, C, -S,
     $                               A( ( 1-ISKEW )*JCH+IOFFG, JCH ),
     $                               ILDA, TEMP, EXTRA )
                        IROW = MAX( 1, JCH-K )
                        IL = MIN( JCH+1, K+2 )
                        EXTRA = ZERO
                        CALL DLAROT( .FALSE., JCH.GT.K, .TRUE., IL, C,
     $                               -S, A( IROW-ISKEW*JCH+IOFFG, JCH ),
     $                               ILDA, EXTRA, TEMP )
                        ICOL = JCH
  150                CONTINUE
  160             CONTINUE
  170          CONTINUE
*
*              If we need lower triangle, copy from upper. Note that
*              the order of copying is chosen to work for 'q' -> 'b'
*
               IF( IPACK.NE.IPACKG .AND. IPACK.NE.3 ) THEN
                  DO 190 JC = 1, N
                     IROW = IOFFST - ISKEW*JC
                     DO 180 JR = JC, MIN( N, JC+UUB )
                        A( JR+IROW, JC ) = A( JC-ISKEW*JR+IOFFG, JR )
  180                CONTINUE
  190             CONTINUE
                  IF( IPACK.EQ.5 ) THEN
                     DO 210 JC = N - UUB + 1, N
                        DO 200 JR = N + 2 - JC, UUB + 1
                           A( JR, JC ) = ZERO
  200                   CONTINUE
  210                CONTINUE
                  END IF
                  IF( IPACKG.EQ.6 ) THEN
                     IPACKG = IPACK
                  ELSE
                     IPACKG = 0
                  END IF
               END IF
            ELSE
*
*              Bottom-Up -- Generate Lower triangle only
*
               IF( IPACK.GE.5 ) THEN
                  IPACKG = 5
                  IF( IPACK.EQ.6 )
     $               IOFFG = 1
               ELSE
                  IPACKG = 2
               END IF
               CALL DCOPY( MNMIN, D, 1, A( 1-ISKEW+IOFFG, 1 ), ILDA+1 )
*
               DO 240 K = 1, UUB
                  DO 230 JC = N - 1, 1, -1
                     IL = MIN( N+1-JC, K+2 )
                     EXTRA = ZERO
                     TEMP = A( 1+( 1-ISKEW )*JC+IOFFG, JC )
                     ANGLE = TWOPI*DLARND( 1, ISEED )
                     C = COS( ANGLE )
                     S = -SIN( ANGLE )
                     CALL DLAROT( .FALSE., .TRUE., N-JC.GT.K, IL, C, S,
     $                            A( ( 1-ISKEW )*JC+IOFFG, JC ), ILDA,
     $                            TEMP, EXTRA )
                     ICOL = MAX( 1, JC-K+1 )
                     CALL DLAROT( .TRUE., .FALSE., .TRUE., JC+2-ICOL, C,
     $                            S, A( JC-ISKEW*ICOL+IOFFG, ICOL ),
     $                            ILDA, DUMMY, TEMP )
*
*                    Chase EXTRA back down the matrix
*
                     ICOL = JC
                     DO 220 JCH = JC + K, N - 1, K
                        CALL DLARTG( A( JCH-ISKEW*ICOL+IOFFG, ICOL ),
     $                               EXTRA, C, S, DUMMY )
                        TEMP = A( 1+( 1-ISKEW )*JCH+IOFFG, JCH )
                        CALL DLAROT( .TRUE., .TRUE., .TRUE., K+2, C, S,
     $                               A( JCH-ISKEW*ICOL+IOFFG, ICOL ),
     $                               ILDA, EXTRA, TEMP )
                        IL = MIN( N+1-JCH, K+2 )
                        EXTRA = ZERO
                        CALL DLAROT( .FALSE., .TRUE., N-JCH.GT.K, IL, C,
     $                               S, A( ( 1-ISKEW )*JCH+IOFFG, JCH ),
     $                               ILDA, TEMP, EXTRA )
                        ICOL = JCH
  220                CONTINUE
  230             CONTINUE
  240          CONTINUE
*
*              If we need upper triangle, copy from lower. Note that
*              the order of copying is chosen to work for 'b' -> 'q'
*
               IF( IPACK.NE.IPACKG .AND. IPACK.NE.4 ) THEN
                  DO 260 JC = N, 1, -1
                     IROW = IOFFST - ISKEW*JC
                     DO 250 JR = JC, MAX( 1, JC-UUB ), -1
                        A( JR+IROW, JC ) = A( JC-ISKEW*JR+IOFFG, JR )
  250                CONTINUE
  260             CONTINUE
                  IF( IPACK.EQ.6 ) THEN
                     DO 280 JC = 1, UUB
                        DO 270 JR = 1, UUB + 1 - JC
                           A( JR, JC ) = ZERO
  270                   CONTINUE
  280                CONTINUE
                  END IF
                  IF( IPACKG.EQ.5 ) THEN
                     IPACKG = IPACK
                  ELSE
                     IPACKG = 0
                  END IF
               END IF
            END IF
         END IF
*
      ELSE
*
*        4)      Generate Banded Matrix by first
*                Rotating by random Unitary matrices,
*                then reducing the bandwidth using Householder
*                transformations.
*
*                Note: we should get here only if LDA .ge. N
*
         IF( ISYM.EQ.1 ) THEN
*
*           Non-symmetric -- A = U D V
*
            CALL DLAGGE( MR, NC, LLB, UUB, D, A, LDA, ISEED, WORK,
     $                   IINFO )
         ELSE
*
*           Symmetric -- A = U D U'
*
            CALL DLAGSY( M, LLB, D, A, LDA, ISEED, WORK, IINFO )
*
         END IF
         IF( IINFO.NE.0 ) THEN
            INFO = 3
            RETURN
         END IF
      END IF
*
*     5)      Pack the matrix
*
      IF( IPACK.NE.IPACKG ) THEN
         IF( IPACK.EQ.1 ) THEN
*
*           'U' -- Upper triangular, not packed
*
            DO 300 J = 1, M
               DO 290 I = J + 1, M
                  A( I, J ) = ZERO
  290          CONTINUE
  300       CONTINUE
*
         ELSE IF( IPACK.EQ.2 ) THEN
*
*           'L' -- Lower triangular, not packed
*
            DO 320 J = 2, M
               DO 310 I = 1, J - 1
                  A( I, J ) = ZERO
  310          CONTINUE
  320       CONTINUE
*
         ELSE IF( IPACK.EQ.3 ) THEN
*
*           'C' -- Upper triangle packed Columnwise.
*
            ICOL = 1
            IROW = 0
            DO 340 J = 1, M
               DO 330 I = 1, J
                  IROW = IROW + 1
                  IF( IROW.GT.LDA ) THEN
                     IROW = 1
                     ICOL = ICOL + 1
                  END IF
                  A( IROW, ICOL ) = A( I, J )
  330          CONTINUE
  340       CONTINUE
*
         ELSE IF( IPACK.EQ.4 ) THEN
*
*           'R' -- Lower triangle packed Columnwise.
*
            ICOL = 1
            IROW = 0
            DO 360 J = 1, M
               DO 350 I = J, M
                  IROW = IROW + 1
                  IF( IROW.GT.LDA ) THEN
                     IROW = 1
                     ICOL = ICOL + 1
                  END IF
                  A( IROW, ICOL ) = A( I, J )
  350          CONTINUE
  360       CONTINUE
*
         ELSE IF( IPACK.GE.5 ) THEN
*
*           'B' -- The lower triangle is packed as a band matrix.
*           'Q' -- The upper triangle is packed as a band matrix.
*           'Z' -- The whole matrix is packed as a band matrix.
*
            IF( IPACK.EQ.5 )
     $         UUB = 0
            IF( IPACK.EQ.6 )
     $         LLB = 0
*
            DO 380 J = 1, UUB
               DO 370 I = MIN( J+LLB, M ), 1, -1
                  A( I-J+UUB+1, J ) = A( I, J )
  370          CONTINUE
  380       CONTINUE
*
            DO 400 J = UUB + 2, N
               DO 390 I = J - UUB, MIN( J+LLB, M )
                  A( I-J+UUB+1, J ) = A( I, J )
  390          CONTINUE
  400       CONTINUE
         END IF
*
*        If packed, zero out extraneous elements.
*
*        Symmetric/Triangular Packed --
*        zero out everything after A(IROW,ICOL)
*
         IF( IPACK.EQ.3 .OR. IPACK.EQ.4 ) THEN
            DO 420 JC = ICOL, M
               DO 410 JR = IROW + 1, LDA
                  A( JR, JC ) = ZERO
  410          CONTINUE
               IROW = 0
  420       CONTINUE
*
         ELSE IF( IPACK.GE.5 ) THEN
*
*           Packed Band --
*              1st row is now in A( UUB+2-j, j), zero above it
*              m-th row is now in A( M+UUB-j,j), zero below it
*              last non-zero diagonal is now in A( UUB+LLB+1,j ),
*                 zero below it, too.
*
            IR1 = UUB + LLB + 2
            IR2 = UUB + M + 2
            DO 450 JC = 1, N
               DO 430 JR = 1, UUB + 1 - JC
                  A( JR, JC ) = ZERO
  430          CONTINUE
               DO 440 JR = MAX( 1, MIN( IR1, IR2-JC ) ), LDA
                  A( JR, JC ) = ZERO
  440          CONTINUE
  450       CONTINUE
         END IF
      END IF
*
      RETURN
*
*     End of DLATMS
*
      END
* =====================================================================
      SUBROUTINE DLATM1( MODE, COND, IRSIGN, IDIST, ISEED, D, N, INFO )
*
*  -- LAPACK auxiliary test routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     September 30, 1994
*
*     .. Scalar Arguments ..
      INTEGER            IDIST, INFO, IRSIGN, MODE, N
      DOUBLE PRECISION   COND
*     ..
*     .. Array Arguments ..
      INTEGER            ISEED( 4 )
      DOUBLE PRECISION   D( * )
*     ..
*
*  Purpose
*  =======
*
*     DLATM1 computes the entries of D(1..N) as specified by
*     MODE, COND and IRSIGN. IDIST and ISEED determine the generation
*     of random numbers. DLATM1 is called by SLATMR to generate
*     random test matrices for LAPACK programs.
*
*  Arguments
*  =========
*
*  MODE   - INTEGER
*           On entry describes how D is to be computed:
*           MODE = 0 means do not change D.
*           MODE = 1 sets D(1)=1 and D(2:N)=1.0/COND
*           MODE = 2 sets D(1:N-1)=1 and D(N)=1.0/COND
*           MODE = 3 sets D(I)=COND**(-(I-1)/(N-1))
*           MODE = 4 sets D(i)=1 - (i-1)/(N-1)*(1 - 1/COND)
*           MODE = 5 sets D to random numbers in the range
*                    ( 1/COND , 1 ) such that their logarithms
*                    are uniformly distributed.
*           MODE = 6 set D to random numbers from same distribution
*                    as the rest of the matrix.
*           MODE < 0 has the same meaning as ABS(MODE), except that
*              the order of the elements of D is reversed.
*           Thus if MODE is positive, D has entries ranging from
*              1 to 1/COND, if negative, from 1/COND to 1,
*           Not modified.
*
*  COND   - DOUBLE PRECISION
*           On entry, used as described under MODE above.
*           If used, it must be >= 1. Not modified.
*
*  IRSIGN - INTEGER
*           On entry, if MODE neither -6, 0 nor 6, determines sign of
*           entries of D
*           0 => leave entries of D unchanged
*           1 => multiply each entry of D by 1 or -1 with probability .5
*
*  IDIST  - CHARACTER*1
*           On entry, IDIST specifies the type of distribution to be
*           used to generate a random matrix .
*           1 => UNIFORM( 0, 1 )
*           2 => UNIFORM( -1, 1 )
*           3 => NORMAL( 0, 1 )
*           Not modified.
*
*  ISEED  - INTEGER array, dimension ( 4 )
*           On entry ISEED specifies the seed of the random number
*           generator. The random number generator uses a
*           linear congruential sequence limited to small
*           integers, and so should produce machine independent
*           random numbers. The values of ISEED are changed on
*           exit, and can be used in the next call to DLATM1
*           to continue the same random number sequence.
*           Changed on exit.
*
*  D      - DOUBLE PRECISION array, dimension ( MIN( M , N ) )
*           Array to be computed according to MODE, COND and IRSIGN.
*           May be changed on exit if MODE is nonzero.
*
*  N      - INTEGER
*           Number of entries of D. Not modified.
*
*  INFO   - INTEGER
*            0  => normal termination
*           -1  => if MODE not in range -6 to 6
*           -2  => if MODE neither -6, 0 nor 6, and
*                  IRSIGN neither 0 nor 1
*           -3  => if MODE neither -6, 0 nor 6 and COND less than 1
*           -4  => if MODE equals 6 or -6 and IDIST not in range 1 to 3
*           -7  => if N negative
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ONE
      PARAMETER          ( ONE = 1.0D0 )
      DOUBLE PRECISION   HALF
      PARAMETER          ( HALF = 0.5D0 )
*     ..
*     .. Local Scalars ..
      INTEGER            I
      DOUBLE PRECISION   ALPHA, TEMP
*     ..
*     .. External Functions ..
      DOUBLE PRECISION   DLARAN
      EXTERNAL           DLARAN
*     ..
*     .. External Subroutines ..
      EXTERNAL           DLARNV, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, DBLE, EXP, LOG
*     ..
*     .. Executable Statements ..
*
*     Decode and Test the input parameters. Initialize flags & seed.
*
      INFO = 0
*
*     Quick return if possible
*
      IF( N.EQ.0 )
     $   RETURN
*
*     Set INFO if an error
*
      IF( MODE.LT.-6 .OR. MODE.GT.6 ) THEN
         INFO = -1
      ELSE IF( ( MODE.NE.-6 .AND. MODE.NE.0 .AND. MODE.NE.6 ) .AND.
     $         ( IRSIGN.NE.0 .AND. IRSIGN.NE.1 ) ) THEN
         INFO = -2
      ELSE IF( ( MODE.NE.-6 .AND. MODE.NE.0 .AND. MODE.NE.6 ) .AND.
     $         COND.LT.ONE ) THEN
         INFO = -3
      ELSE IF( ( MODE.EQ.6 .OR. MODE.EQ.-6 ) .AND.
     $         ( IDIST.LT.1 .OR. IDIST.GT.3 ) ) THEN
         INFO = -4
      ELSE IF( N.LT.0 ) THEN
         INFO = -7
      END IF
*
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DLATM1', -INFO )
         RETURN
      END IF
*
*     Compute D according to COND and MODE
*
      IF( MODE.NE.0 ) THEN
         GO TO ( 10, 30, 50, 70, 90, 110 )ABS( MODE )
*
*        One large D value:
*
   10    CONTINUE
         DO 20 I = 1, N
            D( I ) = ONE / COND
   20    CONTINUE
         D( 1 ) = ONE
         GO TO 120
*
*        One small D value:
*
   30    CONTINUE
         DO 40 I = 1, N
            D( I ) = ONE
   40    CONTINUE
         D( N ) = ONE / COND
         GO TO 120
*
*        Exponentially distributed D values:
*
   50    CONTINUE
         D( 1 ) = ONE
         IF( N.GT.1 ) THEN
            ALPHA = COND**( -ONE / DBLE( N-1 ) )
            DO 60 I = 2, N
               D( I ) = ALPHA**( I-1 )
   60       CONTINUE
         END IF
         GO TO 120
*
*        Arithmetically distributed D values:
*
   70    CONTINUE
         D( 1 ) = ONE
         IF( N.GT.1 ) THEN
            TEMP = ONE / COND
            ALPHA = ( ONE-TEMP ) / DBLE( N-1 )
            DO 80 I = 2, N
               D( I ) = DBLE( N-I )*ALPHA + TEMP
   80       CONTINUE
         END IF
         GO TO 120
*
*        Randomly distributed D values on ( 1/COND , 1):
*
   90    CONTINUE
         ALPHA = LOG( ONE / COND )
         DO 100 I = 1, N
            D( I ) = EXP( ALPHA*DLARAN( ISEED ) )
  100    CONTINUE
         GO TO 120
*
*        Randomly distributed D values from IDIST
*
  110    CONTINUE
         CALL DLARNV( IDIST, ISEED, N, D )
*
  120    CONTINUE
*
*        If MODE neither -6 nor 0 nor 6, and IRSIGN = 1, assign
*        random signs to D
*
         IF( ( MODE.NE.-6 .AND. MODE.NE.0 .AND. MODE.NE.6 ) .AND.
     $       IRSIGN.EQ.1 ) THEN
            DO 130 I = 1, N
               TEMP = DLARAN( ISEED )
               IF( TEMP.GT.HALF )
     $            D( I ) = -D( I )
  130       CONTINUE
         END IF
*
*        Reverse if MODE < 0
*
         IF( MODE.LT.0 ) THEN
            DO 140 I = 1, N / 2
               TEMP = D( I )
               D( I ) = D( N+1-I )
               D( N+1-I ) = TEMP
  140       CONTINUE
         END IF
*
      END IF
*
      RETURN
*
*     End of DLATM1
*
      END
* =====================================================================
      DOUBLE PRECISION FUNCTION DLARND( IDIST, ISEED )
*
*  -- LAPACK auxiliary routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     September 30, 1994
*
*     .. Scalar Arguments ..
      INTEGER            IDIST
*     ..
*     .. Array Arguments ..
      INTEGER            ISEED( 4 )
*     ..
*
*  Purpose
*  =======
*
*  DLARND returns a random real number from a uniform or normal
*  distribution.
*
*  Arguments
*  =========
*
*  IDIST   (input) INTEGER
*          Specifies the distribution of the random numbers:
*          = 1:  uniform (0,1)
*          = 2:  uniform (-1,1)
*          = 3:  normal (0,1)
*
*  ISEED   (input/output) INTEGER array, dimension (4)
*          On entry, the seed of the random number generator; the array
*          elements must be between 0 and 4095, and ISEED(4) must be
*          odd.
*          On exit, the seed is updated.
*
*  Further Details
*  ===============
*
*  This routine calls the auxiliary routine DLARAN to generate a random
*  real number from a uniform (0,1) distribution. The Box-Muller method
*  is used to transform numbers from a uniform to a normal distribution.
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ONE, TWO
      PARAMETER          ( ONE = 1.0D+0, TWO = 2.0D+0 )
      DOUBLE PRECISION   TWOPI
      PARAMETER          ( TWOPI = 6.2831853071795864769252867663D+0 )
*     ..
*     .. Local Scalars ..
      DOUBLE PRECISION   T1, T2
*     ..
*     .. External Functions ..
      DOUBLE PRECISION   DLARAN
      EXTERNAL           DLARAN
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          COS, LOG, SQRT
*     ..
*     .. Executable Statements ..
*
*     Generate a real random number from a uniform (0,1) distribution
*
      T1 = DLARAN( ISEED )
*
      IF( IDIST.EQ.1 ) THEN
*
*        uniform (0,1)
*
         DLARND = T1
      ELSE IF( IDIST.EQ.2 ) THEN
*
*        uniform (-1,1)
*
         DLARND = TWO*T1 - ONE
      ELSE IF( IDIST.EQ.3 ) THEN
*
*        normal (0,1)
*
         T2 = DLARAN( ISEED )
         DLARND = SQRT( -TWO*LOG( T1 ) )*COS( TWOPI*T2 )
      END IF
      RETURN
*
*     End of DLARND
*
      END
* =====================================================================
      SUBROUTINE DLAROT( LROWS, LLEFT, LRIGHT, NL, C, S, A, LDA, XLEFT,
     $                   XRIGHT )
*
*  -- LAPACK auxiliary test routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     February 29, 1992
*
*     .. Scalar Arguments ..
      LOGICAL            LLEFT, LRIGHT, LROWS
      INTEGER            LDA, NL
      DOUBLE PRECISION   C, S, XLEFT, XRIGHT
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   A( * )
*     ..
*
*  Purpose
*  =======
*
*     DLAROT applies a (Givens) rotation to two adjacent rows or
*     columns, where one element of the first and/or last column/row
*     may be a separate variable.  This is specifically indended
*     for use on matrices stored in some format other than GE, so
*     that elements of the matrix may be used or modified for which
*     no array element is provided.
*
*     One example is a symmetric matrix in SB format (bandwidth=4), for
*     which UPLO='L':  Two adjacent rows will have the format:
*
*     row j:     *  *  *  *  *  .  .  .  .
*     row j+1:      *  *  *  *  *  .  .  .  .
*
*     '*' indicates elements for which storage is provided,
*     '.' indicates elements for which no storage is provided, but
*     are not necessarily zero; their values are determined by
*     symmetry.  ' ' indicates elements which are necessarily zero,
*      and have no storage provided.
*
*     Those columns which have two '*'s can be handled by DROT.
*     Those columns which have no '*'s can be ignored, since as long
*     as the Givens rotations are carefully applied to preserve
*     symmetry, their values are determined.
*     Those columns which have one '*' have to be handled separately,
*     by using separate variables "p" and "q":
*
*     row j:     *  *  *  *  *  p  .  .  .
*     row j+1:   q  *  *  *  *  *  .  .  .  .
*
*     The element p would have to be set correctly, then that column
*     is rotated, setting p to its new value.  The next call to
*     DLAROT would rotate columns j and j+1, using p, and restore
*     symmetry.  The element q would start out being zero, and be
*     made non-zero by the rotation.  Later, rotations would presumably
*     be chosen to zero q out.
*
*     Typical Calling Sequences: rotating the i-th and (i+1)-st rows.
*     ------- ------- ---------
*
*       General dense matrix:
*
*               CALL DLAROT(.TRUE.,.FALSE.,.FALSE., N, C,S,
*                       A(i,1),LDA, DUMMY, DUMMY)
*
*       General banded matrix in GB format:
*
*               j = MAX(1, i-KL )
*               NL = MIN( N, i+KU+1 ) + 1-j
*               CALL DLAROT( .TRUE., i-KL.GE.1, i+KU.LT.N, NL, C,S,
*                       A(KU+i+1-j,j),LDA-1, XLEFT, XRIGHT )
*
*               [ note that i+1-j is just MIN(i,KL+1) ]
*
*       Symmetric banded matrix in SY format, bandwidth K,
*       lower triangle only:
*
*               j = MAX(1, i-K )
*               NL = MIN( K+1, i ) + 1
*               CALL DLAROT( .TRUE., i-K.GE.1, .TRUE., NL, C,S,
*                       A(i,j), LDA, XLEFT, XRIGHT )
*
*       Same, but upper triangle only:
*
*               NL = MIN( K+1, N-i ) + 1
*               CALL DLAROT( .TRUE., .TRUE., i+K.LT.N, NL, C,S,
*                       A(i,i), LDA, XLEFT, XRIGHT )
*
*       Symmetric banded matrix in SB format, bandwidth K,
*       lower triangle only:
*
*               [ same as for SY, except:]
*                   . . . .
*                       A(i+1-j,j), LDA-1, XLEFT, XRIGHT )
*
*               [ note that i+1-j is just MIN(i,K+1) ]
*
*       Same, but upper triangle only:
*                    . . .
*                       A(K+1,i), LDA-1, XLEFT, XRIGHT )
*
*       Rotating columns is just the transpose of rotating rows, except
*       for GB and SB: (rotating columns i and i+1)
*
*       GB:
*               j = MAX(1, i-KU )
*               NL = MIN( N, i+KL+1 ) + 1-j
*               CALL DLAROT( .TRUE., i-KU.GE.1, i+KL.LT.N, NL, C,S,
*                       A(KU+j+1-i,i),LDA-1, XTOP, XBOTTM )
*
*               [note that KU+j+1-i is just MAX(1,KU+2-i)]
*
*       SB: (upper triangle)
*
*                    . . . . . .
*                       A(K+j+1-i,i),LDA-1, XTOP, XBOTTM )
*
*       SB: (lower triangle)
*
*                    . . . . . .
*                       A(1,i),LDA-1, XTOP, XBOTTM )
*
*  Arguments
*  =========
*
*  LROWS  - LOGICAL
*           If .TRUE., then DLAROT will rotate two rows.  If .FALSE.,
*           then it will rotate two columns.
*           Not modified.
*
*  LLEFT  - LOGICAL
*           If .TRUE., then XLEFT will be used instead of the
*           corresponding element of A for the first element in the
*           second row (if LROWS=.FALSE.) or column (if LROWS=.TRUE.)
*           If .FALSE., then the corresponding element of A will be
*           used.
*           Not modified.
*
*  LRIGHT - LOGICAL
*           If .TRUE., then XRIGHT will be used instead of the
*           corresponding element of A for the last element in the
*           first row (if LROWS=.FALSE.) or column (if LROWS=.TRUE.) If
*           .FALSE., then the corresponding element of A will be used.
*           Not modified.
*
*  NL     - INTEGER
*           The length of the rows (if LROWS=.TRUE.) or columns (if
*           LROWS=.FALSE.) to be rotated.  If XLEFT and/or XRIGHT are
*           used, the columns/rows they are in should be included in
*           NL, e.g., if LLEFT = LRIGHT = .TRUE., then NL must be at
*           least 2.  The number of rows/columns to be rotated
*           exclusive of those involving XLEFT and/or XRIGHT may
*           not be negative, i.e., NL minus how many of LLEFT and
*           LRIGHT are .TRUE. must be at least zero; if not, XERBLA
*           will be called.
*           Not modified.
*
*  C, S   - DOUBLE PRECISION
*           Specify the Givens rotation to be applied.  If LROWS is
*           true, then the matrix ( c  s )
*                                 (-s  c )  is applied from the left;
*           if false, then the transpose thereof is applied from the
*           right.  For a Givens rotation, C**2 + S**2 should be 1,
*           but this is not checked.
*           Not modified.
*
*  A      - DOUBLE PRECISION array.
*           The array containing the rows/columns to be rotated.  The
*           first element of A should be the upper left element to
*           be rotated.
*           Read and modified.
*
*  LDA    - INTEGER
*           The "effective" leading dimension of A.  If A contains
*           a matrix stored in GE or SY format, then this is just
*           the leading dimension of A as dimensioned in the calling
*           routine.  If A contains a matrix stored in band (GB or SB)
*           format, then this should be *one less* than the leading
*           dimension used in the calling routine.  Thus, if
*           A were dimensioned A(LDA,*) in DLAROT, then A(1,j) would
*           be the j-th element in the first of the two rows
*           to be rotated, and A(2,j) would be the j-th in the second,
*           regardless of how the array may be stored in the calling
*           routine.  [A cannot, however, actually be dimensioned thus,
*           since for band format, the row number may exceed LDA, which
*           is not legal FORTRAN.]
*           If LROWS=.TRUE., then LDA must be at least 1, otherwise
*           it must be at least NL minus the number of .TRUE. values
*           in XLEFT and XRIGHT.
*           Not modified.
*
*  XLEFT  - DOUBLE PRECISION
*           If LLEFT is .TRUE., then XLEFT will be used and modified
*           instead of A(2,1) (if LROWS=.TRUE.) or A(1,2)
*           (if LROWS=.FALSE.).
*           Read and modified.
*
*  XRIGHT - DOUBLE PRECISION
*           If LRIGHT is .TRUE., then XRIGHT will be used and modified
*           instead of A(1,NL) (if LROWS=.TRUE.) or A(NL,1)
*           (if LROWS=.FALSE.).
*           Read and modified.
*
*  =====================================================================
*
*     .. Local Scalars ..
      INTEGER            IINC, INEXT, IX, IY, IYT, NT
*     ..
*     .. Local Arrays ..
      DOUBLE PRECISION   XT( 2 ), YT( 2 )
*     ..
*     .. External Subroutines ..
      EXTERNAL           DROT, XERBLA
*     ..
*     .. Executable Statements ..
*
*     Set up indices, arrays for ends
*
      IF( LROWS ) THEN
         IINC = LDA
         INEXT = 1
      ELSE
         IINC = 1
         INEXT = LDA
      END IF
*
      IF( LLEFT ) THEN
         NT = 1
         IX = 1 + IINC
         IY = 2 + LDA
         XT( 1 ) = A( 1 )
         YT( 1 ) = XLEFT
      ELSE
         NT = 0
         IX = 1
         IY = 1 + INEXT
      END IF
*
      IF( LRIGHT ) THEN
         IYT = 1 + INEXT + ( NL-1 )*IINC
         NT = NT + 1
         XT( NT ) = XRIGHT
         YT( NT ) = A( IYT )
      END IF
*
*     Check for errors
*
      IF( NL.LT.NT ) THEN
         CALL XERBLA( 'DLAROT', 4 )
         RETURN
      END IF
      IF( LDA.LE.0 .OR. ( .NOT.LROWS .AND. LDA.LT.NL-NT ) ) THEN
         CALL XERBLA( 'DLAROT', 8 )
         RETURN
      END IF
*
*     Rotate
*
      CALL DROT( NL-NT, A( IX ), IINC, A( IY ), IINC, C, S )
      CALL DROT( NT, XT, 1, YT, 1, C, S )
*
*     Stuff values back into XLEFT, XRIGHT, etc.
*
      IF( LLEFT ) THEN
         A( 1 ) = XT( 1 )
         XLEFT = YT( 1 )
      END IF
*
      IF( LRIGHT ) THEN
         XRIGHT = XT( NT )
         A( IYT ) = YT( NT )
      END IF
*
      RETURN
*
*     End of DLAROT
*
      END
* =====================================================================
      SUBROUTINE DLAGGE( M, N, KL, KU, D, A, LDA, ISEED, WORK, INFO )
*
*  -- LAPACK auxiliary test routine (version 3.0)
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     February 29, 1992
*
*     .. Scalar Arguments ..
      INTEGER            INFO, KL, KU, LDA, M, N
*     ..
*     .. Array Arguments ..
      INTEGER            ISEED( 4 )
      DOUBLE PRECISION   A( LDA, * ), D( * ), WORK( * )
*     ..
*
*  Purpose
*  =======
*
*  DLAGGE generates a real general m by n matrix A, by pre- and post-
*  multiplying a real diagonal matrix D with random orthogonal matrices:
*  A = U*D*V. The lower and upper bandwidths may then be reduced to
*  kl and ku by additional orthogonal transformations.
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
*  KL      (input) INTEGER
*          The number of nonzero subdiagonals within the band of A.
*          0 <= KL <= M-1.
*
*  KU      (input) INTEGER
*          The number of nonzero superdiagonals within the band of A.
*          0 <= KU <= N-1.
*
*  D       (input) DOUBLE PRECISION array, dimension (min(M,N))
*          The diagonal elements of the diagonal matrix D.
*
*  A       (output) DOUBLE PRECISION array, dimension (LDA,N)
*          The generated m by n matrix A.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= M.
*
*  ISEED   (input/output) INTEGER array, dimension (4)
*          On entry, the seed of the random number generator; the array
*          elements must be between 0 and 4095, and ISEED(4) must be
*          odd.
*          On exit, the seed is updated.
*
*  WORK    (workspace) DOUBLE PRECISION array, dimension (M+N)
*
*  INFO    (output) INTEGER
*          = 0: successful exit
*          < 0: if INFO = -i, the i-th argument had an illegal value
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            I, J
      DOUBLE PRECISION   TAU, WA, WB, WN
*     ..
*     .. External Subroutines ..
      EXTERNAL           DGEMV, DGER, DLARNV, DSCAL, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN, SIGN
*     ..
*     .. External Functions ..
      DOUBLE PRECISION   DNRM2
      EXTERNAL           DNRM2
*     ..
*     .. Executable Statements ..
*
*     Test the input arguments
*
      INFO = 0
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( KL.LT.0 .OR. KL.GT.M-1 ) THEN
         INFO = -3
      ELSE IF( KU.LT.0 .OR. KU.GT.N-1 ) THEN
         INFO = -4
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -7
      END IF
      IF( INFO.LT.0 ) THEN
         CALL XERBLA( 'DLAGGE', -INFO )
         RETURN
      END IF
*
*     initialize A to diagonal matrix
*
      DO 20 J = 1, N
         DO 10 I = 1, M
            A( I, J ) = ZERO
   10    CONTINUE
   20 CONTINUE
      DO 30 I = 1, MIN( M, N )
         A( I, I ) = D( I )
   30 CONTINUE
*
*     pre- and post-multiply A by random orthogonal matrices
*
      DO 40 I = MIN( M, N ), 1, -1
         IF( I.LT.M ) THEN
*
*           generate random reflection
*
            CALL DLARNV( 3, ISEED, M-I+1, WORK )
            WN = DNRM2( M-I+1, WORK, 1 )
            WA = SIGN( WN, WORK( 1 ) )
            IF( WN.EQ.ZERO ) THEN
               TAU = ZERO
            ELSE
               WB = WORK( 1 ) + WA
               CALL DSCAL( M-I, ONE / WB, WORK( 2 ), 1 )
               WORK( 1 ) = ONE
               TAU = WB / WA
            END IF
*
*           multiply A(i:m,i:n) by random reflection from the left
*
            CALL DGEMV( 'Transpose', M-I+1, N-I+1, ONE, A( I, I ), LDA,
     $                  WORK, 1, ZERO, WORK( M+1 ), 1 )
            CALL DGER( M-I+1, N-I+1, -TAU, WORK, 1, WORK( M+1 ), 1,
     $                 A( I, I ), LDA )
         END IF
         IF( I.LT.N ) THEN
*
*           generate random reflection
*
            CALL DLARNV( 3, ISEED, N-I+1, WORK )
            WN = DNRM2( N-I+1, WORK, 1 )
            WA = SIGN( WN, WORK( 1 ) )
            IF( WN.EQ.ZERO ) THEN
               TAU = ZERO
            ELSE
               WB = WORK( 1 ) + WA
               CALL DSCAL( N-I, ONE / WB, WORK( 2 ), 1 )
               WORK( 1 ) = ONE
               TAU = WB / WA
            END IF
*
*           multiply A(i:m,i:n) by random reflection from the right
*
            CALL DGEMV( 'No transpose', M-I+1, N-I+1, ONE, A( I, I ),
     $                  LDA, WORK, 1, ZERO, WORK( N+1 ), 1 )
            CALL DGER( M-I+1, N-I+1, -TAU, WORK( N+1 ), 1, WORK, 1,
     $                 A( I, I ), LDA )
         END IF
   40 CONTINUE
*
*     Reduce number of subdiagonals to KL and number of superdiagonals
*     to KU
*
      DO 70 I = 1, MAX( M-1-KL, N-1-KU )
         IF( KL.LE.KU ) THEN
*
*           annihilate subdiagonal elements first (necessary if KL = 0)
*
            IF( I.LE.MIN( M-1-KL, N ) ) THEN
*
*              generate reflection to annihilate A(kl+i+1:m,i)
*
               WN = DNRM2( M-KL-I+1, A( KL+I, I ), 1 )
               WA = SIGN( WN, A( KL+I, I ) )
               IF( WN.EQ.ZERO ) THEN
                  TAU = ZERO
               ELSE
                  WB = A( KL+I, I ) + WA
                  CALL DSCAL( M-KL-I, ONE / WB, A( KL+I+1, I ), 1 )
                  A( KL+I, I ) = ONE
                  TAU = WB / WA
               END IF
*
*              apply reflection to A(kl+i:m,i+1:n) from the left
*
               CALL DGEMV( 'Transpose', M-KL-I+1, N-I, ONE,
     $                     A( KL+I, I+1 ), LDA, A( KL+I, I ), 1, ZERO,
     $                     WORK, 1 )
               CALL DGER( M-KL-I+1, N-I, -TAU, A( KL+I, I ), 1, WORK, 1,
     $                    A( KL+I, I+1 ), LDA )
               A( KL+I, I ) = -WA
            END IF
*
            IF( I.LE.MIN( N-1-KU, M ) ) THEN
*
*              generate reflection to annihilate A(i,ku+i+1:n)
*
               WN = DNRM2( N-KU-I+1, A( I, KU+I ), LDA )
               WA = SIGN( WN, A( I, KU+I ) )
               IF( WN.EQ.ZERO ) THEN
                  TAU = ZERO
               ELSE
                  WB = A( I, KU+I ) + WA
                  CALL DSCAL( N-KU-I, ONE / WB, A( I, KU+I+1 ), LDA )
                  A( I, KU+I ) = ONE
                  TAU = WB / WA
               END IF
*
*              apply reflection to A(i+1:m,ku+i:n) from the right
*
               CALL DGEMV( 'No transpose', M-I, N-KU-I+1, ONE,
     $                     A( I+1, KU+I ), LDA, A( I, KU+I ), LDA, ZERO,
     $                     WORK, 1 )
               CALL DGER( M-I, N-KU-I+1, -TAU, WORK, 1, A( I, KU+I ),
     $                    LDA, A( I+1, KU+I ), LDA )
               A( I, KU+I ) = -WA
            END IF
         ELSE
*
*           annihilate superdiagonal elements first (necessary if
*           KU = 0)
*
            IF( I.LE.MIN( N-1-KU, M ) ) THEN
*
*              generate reflection to annihilate A(i,ku+i+1:n)
*
               WN = DNRM2( N-KU-I+1, A( I, KU+I ), LDA )
               WA = SIGN( WN, A( I, KU+I ) )
               IF( WN.EQ.ZERO ) THEN
                  TAU = ZERO
               ELSE
                  WB = A( I, KU+I ) + WA
                  CALL DSCAL( N-KU-I, ONE / WB, A( I, KU+I+1 ), LDA )
                  A( I, KU+I ) = ONE
                  TAU = WB / WA
               END IF
*
*              apply reflection to A(i+1:m,ku+i:n) from the right
*
               CALL DGEMV( 'No transpose', M-I, N-KU-I+1, ONE,
     $                     A( I+1, KU+I ), LDA, A( I, KU+I ), LDA, ZERO,
     $                     WORK, 1 )
               CALL DGER( M-I, N-KU-I+1, -TAU, WORK, 1, A( I, KU+I ),
     $                    LDA, A( I+1, KU+I ), LDA )
               A( I, KU+I ) = -WA
            END IF
*
            IF( I.LE.MIN( M-1-KL, N ) ) THEN
*
*              generate reflection to annihilate A(kl+i+1:m,i)
*
               WN = DNRM2( M-KL-I+1, A( KL+I, I ), 1 )
               WA = SIGN( WN, A( KL+I, I ) )
               IF( WN.EQ.ZERO ) THEN
                  TAU = ZERO
               ELSE
                  WB = A( KL+I, I ) + WA
                  CALL DSCAL( M-KL-I, ONE / WB, A( KL+I+1, I ), 1 )
                  A( KL+I, I ) = ONE
                  TAU = WB / WA
               END IF
*
*              apply reflection to A(kl+i:m,i+1:n) from the left
*
               CALL DGEMV( 'Transpose', M-KL-I+1, N-I, ONE,
     $                     A( KL+I, I+1 ), LDA, A( KL+I, I ), 1, ZERO,
     $                     WORK, 1 )
               CALL DGER( M-KL-I+1, N-I, -TAU, A( KL+I, I ), 1, WORK, 1,
     $                    A( KL+I, I+1 ), LDA )
               A( KL+I, I ) = -WA
            END IF
         END IF
*
         DO 50 J = KL + I + 1, M
            A( J, I ) = ZERO
   50    CONTINUE
*
         DO 60 J = KU + I + 1, N
            A( I, J ) = ZERO
   60    CONTINUE
   70 CONTINUE
      RETURN
*
*     End of DLAGGE
*
      END
* =====================================================================
      SUBROUTINE DLAGSY( N, K, D, A, LDA, ISEED, WORK, INFO )
*
*  -- LAPACK auxiliary test routine (version 3.0)
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     February 29, 1992
*
*     .. Scalar Arguments ..
      INTEGER            INFO, K, LDA, N
*     ..
*     .. Array Arguments ..
      INTEGER            ISEED( 4 )
      DOUBLE PRECISION   A( LDA, * ), D( * ), WORK( * )
*     ..
*
*  Purpose
*  =======
*
*  DLAGSY generates a real symmetric matrix A, by pre- and post-
*  multiplying a real diagonal matrix D with a random orthogonal matrix:
*  A = U*D*U'. The semi-bandwidth may then be reduced to k by additional
*  orthogonal transformations.
*
*  Arguments
*  =========
*
*  N       (input) INTEGER
*          The order of the matrix A.  N >= 0.
*
*  K       (input) INTEGER
*          The number of nonzero subdiagonals within the band of A.
*          0 <= K <= N-1.
*
*  D       (input) DOUBLE PRECISION array, dimension (N)
*          The diagonal elements of the diagonal matrix D.
*
*  A       (output) DOUBLE PRECISION array, dimension (LDA,N)
*          The generated n by n symmetric matrix A (the full matrix is
*          stored).
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= N.
*
*  ISEED   (input/output) INTEGER array, dimension (4)
*          On entry, the seed of the random number generator; the array
*          elements must be between 0 and 4095, and ISEED(4) must be
*          odd.
*          On exit, the seed is updated.
*
*  WORK    (workspace) DOUBLE PRECISION array, dimension (2*N)
*
*  INFO    (output) INTEGER
*          = 0: successful exit
*          < 0: if INFO = -i, the i-th argument had an illegal value
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE, HALF
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0, HALF = 0.5D+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            I, J
      DOUBLE PRECISION   ALPHA, TAU, WA, WB, WN
*     ..
*     .. External Subroutines ..
      EXTERNAL           DAXPY, DGEMV, DGER, DLARNV, DSCAL, DSYMV,
     $                   DSYR2, XERBLA
*     ..
*     .. External Functions ..
      DOUBLE PRECISION   DDOT, DNRM2
      EXTERNAL           DDOT, DNRM2
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, SIGN
*     ..
*     .. Executable Statements ..
*
*     Test the input arguments
*
      INFO = 0
      IF( N.LT.0 ) THEN
         INFO = -1
      ELSE IF( K.LT.0 .OR. K.GT.N-1 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -5
      END IF
      IF( INFO.LT.0 ) THEN
         CALL XERBLA( 'DLAGSY', -INFO )
         RETURN
      END IF
*
*     initialize lower triangle of A to diagonal matrix
*
      DO 20 J = 1, N
         DO 10 I = J + 1, N
            A( I, J ) = ZERO
   10    CONTINUE
   20 CONTINUE
      DO 30 I = 1, N
         A( I, I ) = D( I )
   30 CONTINUE
*
*     Generate lower triangle of symmetric matrix
*
      DO 40 I = N - 1, 1, -1
*
*        generate random reflection
*
         CALL DLARNV( 3, ISEED, N-I+1, WORK )
         WN = DNRM2( N-I+1, WORK, 1 )
         WA = SIGN( WN, WORK( 1 ) )
         IF( WN.EQ.ZERO ) THEN
            TAU = ZERO
         ELSE
            WB = WORK( 1 ) + WA
            CALL DSCAL( N-I, ONE / WB, WORK( 2 ), 1 )
            WORK( 1 ) = ONE
            TAU = WB / WA
         END IF
*
*        apply random reflection to A(i:n,i:n) from the left
*        and the right
*
*        compute  y := tau * A * u
*
         CALL DSYMV( 'Lower', N-I+1, TAU, A( I, I ), LDA, WORK, 1, ZERO,
     $               WORK( N+1 ), 1 )
*
*        compute  v := y - 1/2 * tau * ( y, u ) * u
*
         ALPHA = -HALF*TAU*DDOT( N-I+1, WORK( N+1 ), 1, WORK, 1 )
         CALL DAXPY( N-I+1, ALPHA, WORK, 1, WORK( N+1 ), 1 )
*
*        apply the transformation as a rank-2 update to A(i:n,i:n)
*
         CALL DSYR2( 'Lower', N-I+1, -ONE, WORK, 1, WORK( N+1 ), 1,
     $               A( I, I ), LDA )
   40 CONTINUE
*
*     Reduce number of subdiagonals to K
*
      DO 60 I = 1, N - 1 - K
*
*        generate reflection to annihilate A(k+i+1:n,i)
*
         WN = DNRM2( N-K-I+1, A( K+I, I ), 1 )
         WA = SIGN( WN, A( K+I, I ) )
         IF( WN.EQ.ZERO ) THEN
            TAU = ZERO
         ELSE
            WB = A( K+I, I ) + WA
            CALL DSCAL( N-K-I, ONE / WB, A( K+I+1, I ), 1 )
            A( K+I, I ) = ONE
            TAU = WB / WA
         END IF
*
*        apply reflection to A(k+i:n,i+1:k+i-1) from the left
*
         CALL DGEMV( 'Transpose', N-K-I+1, K-1, ONE, A( K+I, I+1 ), LDA,
     $               A( K+I, I ), 1, ZERO, WORK, 1 )
         CALL DGER( N-K-I+1, K-1, -TAU, A( K+I, I ), 1, WORK, 1,
     $              A( K+I, I+1 ), LDA )
*
*        apply reflection to A(k+i:n,k+i:n) from the left and the right
*
*        compute  y := tau * A * u
*
         CALL DSYMV( 'Lower', N-K-I+1, TAU, A( K+I, K+I ), LDA,
     $               A( K+I, I ), 1, ZERO, WORK, 1 )
*
*        compute  v := y - 1/2 * tau * ( y, u ) * u
*
         ALPHA = -HALF*TAU*DDOT( N-K-I+1, WORK, 1, A( K+I, I ), 1 )
         CALL DAXPY( N-K-I+1, ALPHA, A( K+I, I ), 1, WORK, 1 )
*
*        apply symmetric rank-2 update to A(k+i:n,k+i:n)
*
         CALL DSYR2( 'Lower', N-K-I+1, -ONE, A( K+I, I ), 1, WORK, 1,
     $               A( K+I, K+I ), LDA )
*
         A( K+I, I ) = -WA
         DO 50 J = K + I + 1, N
            A( J, I ) = ZERO
   50    CONTINUE
   60 CONTINUE
*
*     Store full symmetric matrix
*
      DO 80 J = 1, N
         DO 70 I = J + 1, N
            A( J, I ) = A( I, J )
   70    CONTINUE
   80 CONTINUE
      RETURN
*
*     End of DLAGSY
*
      END
* =====================================================================
      DOUBLE PRECISION FUNCTION DLARAN( ISEED )
*
*  -- LAPACK auxiliary routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     February 29, 1992
*
*     .. Array Arguments ..
      INTEGER            ISEED( 4 )
*     ..
*
*  Purpose
*  =======
*
*  DLARAN returns a random real number from a uniform (0,1)
*  distribution.
*
*  Arguments
*  =========
*
*  ISEED   (input/output) INTEGER array, dimension (4)
*          On entry, the seed of the random number generator; the array
*          elements must be between 0 and 4095, and ISEED(4) must be
*          odd.
*          On exit, the seed is updated.
*
*  Further Details
*  ===============
*
*  This routine uses a multiplicative congruential method with modulus
*  2**48 and multiplier 33952834046453 (see G.S.Fishman,
*  'Multiplicative congruential random number generators with modulus
*  2**b: an exhaustive analysis for b = 32 and a partial analysis for
*  b = 48', Math. Comp. 189, pp 331-344, 1990).
*
*  48-bit integers are stored in 4 integer array elements with 12 bits
*  per element. Hence the routine is portable across machines with
*  integers of 32 bits or more.
*
*  =====================================================================
*
*     .. Parameters ..
      INTEGER            M1, M2, M3, M4
      PARAMETER          ( M1 = 494, M2 = 322, M3 = 2508, M4 = 2549 )
      DOUBLE PRECISION   ONE
      PARAMETER          ( ONE = 1.0D+0 )
      INTEGER            IPW2
      DOUBLE PRECISION   R
      PARAMETER          ( IPW2 = 4096, R = ONE / IPW2 )
*     ..
*     .. Local Scalars ..
      INTEGER            IT1, IT2, IT3, IT4
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          DBLE, MOD
*     ..
*     .. Executable Statements ..
*
*     multiply the seed by the multiplier modulo 2**48
*
      IT4 = ISEED( 4 )*M4
      IT3 = IT4 / IPW2
      IT4 = IT4 - IPW2*IT3
      IT3 = IT3 + ISEED( 3 )*M4 + ISEED( 4 )*M3
      IT2 = IT3 / IPW2
      IT3 = IT3 - IPW2*IT2
      IT2 = IT2 + ISEED( 2 )*M4 + ISEED( 3 )*M3 + ISEED( 4 )*M2
      IT1 = IT2 / IPW2
      IT2 = IT2 - IPW2*IT1
      IT1 = IT1 + ISEED( 1 )*M4 + ISEED( 2 )*M3 + ISEED( 3 )*M2 +
     $      ISEED( 4 )*M1
      IT1 = MOD( IT1, IPW2 )
*
*     return updated seed
*
      ISEED( 1 ) = IT1
      ISEED( 2 ) = IT2
      ISEED( 3 ) = IT3
      ISEED( 4 ) = IT4
*
*     convert 48-bit integer to a real number in the interval (0,1)
*
      DLARAN = R*( DBLE( IT1 )+R*( DBLE( IT2 )+R*( DBLE( IT3 )+R*
     $         ( DBLE( IT4 ) ) ) ) )
      RETURN
*
*     End of DLARAN
*
      END
* =====================================================================

