      SUBROUTINE DGHUTP( JOB, COMPQ1, COMPQ2, N, A, LDA, DE, LDDE, C1,
     $                   LDC1, VW, LDVW, Q1, LDQ1, Q2, LDQ2, B, LDB, F,
     $                   LDF, C2, LDC2, ALPHAR, ALPHAI, BETA, IWORK,
     $                   LIWORK, DWORK, LDWORK, INFO )
C
C     PURPOSE
C
C     To compute the eigenvalues of a real N-by-N skew-Hamiltonian/
C     Hamiltonian pencil aS - bH with
C
C           (  A  D  )         (  C  V  )
C       S = (        ) and H = (        ).                           (1)
C           (  E  A' )         (  W -C' )
C
C     Optionally, if JOB = 'T', decompositions of S and H will be
C     computed via orthogonal transformations Q1 and Q2 as follows:
C
C                       (  Aout  Dout  )
C       Q1' S J Q1 J' = (              ),
C                       (   0    Aout' )
C
C                       (  Bout  Fout  )
C       J' Q2' J S Q2 = (              ) =: T,                       (2)
C                       (   0    Bout' )
C
C                  (  C1out  Vout  )            (  0  I  )
C       Q1' H Q2 = (               ), where J = (        )
C                  (  0     C2out' )            ( -I  0  )
C
C     and Aout, Bout, C1out are upper triangular, C2out is upper quasi-
C     triangular and Dout and Fout are skew-symmetric. The notation M'
C     denotes the transpose of the matrix M.
C     Optionally, if COMPQ1 = 'I' or COMPQ1 = 'U', then the orthogonal
C     transformation matrix Q1 will be computed.
C     Optionally, if COMPQ2 = 'I' or COMPQ2 = 'U', then the orthogonal
C     transformation matrix Q2 will be computed.
C
C     ARGUMENTS
C
C     Mode Parameters
C
C     JOB     CHARACTER*1
C             Specifies the computation to be performed, as follows:
C             = 'E': compute the eigenvalues only; S and H will not
C                    necessarily be transformed as in (2).
C             = 'T': put S and H into the forms in (2) and return the
C                    eigenvalues in ALPHAR, ALPHAI and BETA.
C
C     COMPQ1  CHARACTER*1
C             Specifies whether to compute the orthogonal transformation
C             matrix Q1, as follows:
C             = 'N':  Q1 is not computed;
C             = 'I':  the array Q1 is initialized internally to the unit
C                     matrix, and the orthogonal matrix Q1 is returned;
C             = 'U':  the array Q1 contains an orthogonal matrix Q on
C                     entry, and the product Q*Q1 is returned, where Q1
C                     is the product of the orthogonal transformations
C                     that are applied to the pencil aS - bH to reduce
C                     S and H to the forms in (2), for COMPQ1 = 'I'.
C
C     COMPQ2  CHARACTER*1
C             Specifies whether to compute the orthogonal transformation
C             matrix Q2, as follows:
C             = 'N':  Q2 is not computed;
C             = 'I':  on exit, the array Q2 contains the orthogonal
C                     matrix Q2;
C             = 'U':  on exit, the array Q2 contains the matrix product
C                     J*Q*J'*Q2, where Q2 is the product of the
C                     orthogonal transformations that are applied to
C                     the pencil aS - bH to reduce S and H to the forms
C                     in (2), for COMPQ2 = 'I'.
C                     Setting COMPQ2 <> 'N' assumes COMPQ2 = COMPQ1.
C
C     Input/Output Parameters
C
C     N       (input) INTEGER
C             The order of the pencil aS - bH.  N >= 0, even.
C
C     A       (input/output) DOUBLE PRECISION array, dimension
C                            (LDA, N/2)
C             On entry, the leading N/2-by-N/2 part of this array must
C             contain the matrix A.
C             On exit, if JOB = 'T', the leading N/2-by-N/2 part of this
C             array contains the matrix Aout; otherwise, it contains the
C             upper triangular matrix A obtained just before the
C             application of the periodic QZ algorithm.
C
C     LDA     INTEGER
C             The leading dimension of the array A.  LDA >= MAX(1, N/2).
C
C     DE      (input/output) DOUBLE PRECISION array, dimension
C                            (LDDE, N/2+1)
C             On entry, the leading N/2-by-N/2 strictly lower triangular
C             part of this array must contain the strictly lower
C             triangular part of the skew-symmetric matrix E, and the
C             N/2-by-N/2 strictly upper triangular part of the submatrix
C             in the columns 2 to N/2+1 of this array must contain the
C             strictly upper triangular part of the skew-symmetric
C             matrix D.
C             The entries on the diagonal and the first superdiagonal of
C             this array need not be set, but are assumed to be zero.
C             On exit, if JOB = 'T', the leading N/2-by-N/2 strictly
C             upper triangular part of the submatrix in the columns
C             2 to N/2+1 of this array contains the strictly upper
C             triangular part of the skew-symmetric matrix Dout.
C             If JOB = 'E', the leading N/2-by-N/2 strictly upper
C             triangular part of the submatrix in the columns 2 to N/2+1
C             of this array contains the strictly upper triangular part
C             of the skew-symmetric matrix D just before the application
C             of the periodic QZ algorithm. The remaining entries are
C             meaningless.
C
C     LDDE    INTEGER
C             The leading dimension of the array DE.
C             LDDE >= MAX(1, N/2).
C
C     C1      (input/output) DOUBLE PRECISION array, dimension
C                            (LDC1, N/2)
C             On entry, the leading N/2-by-N/2 part of this array must
C             contain the matrix C1 = C.
C             On exit, if JOB = 'T', the leading N/2-by-N/2 part of this
C             array contains the matrix C1out; otherwise, it contains the
C             upper triangular matrix C1 obtained just before the
C             application of the periodic QZ algorithm.
C
C     LDC1    INTEGER
C             The leading dimension of the array C1.
C             LDC1 >= MAX(1, N/2).
C
C     VW      (input/output) DOUBLE PRECISION array, dimension
C                            (LDVW, N/2+1)
C             On entry, the leading N/2-by-N/2 lower triangular part of
C             this array must contain the lower triangular part of the
C             symmetric matrix W, and the N/2-by-N/2 upper triangular
C             part of the submatrix in the columns 2 to N/2+1 of this
C             array must contain the upper triangular part of the
C             symmetric matrix V.
C             On exit, if JOB = 'T', the N/2-by-N/2 part in the columns
C             2 to N/2+1 of this array contains the matrix Vout.
C             If JOB = 'E', the N/2-by-N/2 part in the columns 2 to
C             N/2+1 of this array contains the matrix V just before the
C             application of the periodic QZ algorithm.
C
C     LDVW    INTEGER
C             The leading dimension of the array VW.
C             LDVW >= MAX(1, N/2).
C
C     Q1      (input/output) DOUBLE PRECISION array, dimension (LDQ1, N)
C             On entry, if COMPQ1 = 'U', then the leading N-by-N part of
C             this array must contain a given matrix Q, and on exit,
C             the leading N-by-N part of this array contains the product
C             of the input matrix Q and the transformation matrix Q1
C             used to transform the matrices S and H.
C             On exit, if COMPQ1 = 'I', then the leading N-by-N part of
C             this array contains the orthogonal transformation matrix
C             Q1.
C             If COMPQ1 = 'N', this array is not referenced.
C
C     LDQ1    INTEGER
C             The leading dimension of the array Q1.
C             LDQ1 >= 1,         if COMPQ1 = 'N';
C             LDQ1 >= MAX(1, N), if COMPQ1 = 'I' or COMPQ1 = 'U'.
C
C     Q2      (output) DOUBLE PRECISION array, dimension (LDQ2, N)
C             On exit, if COMPQ2 = 'U', then the leading N-by-N part of
C             this array contains the product of the matrix J*Q*J' and
C             the transformation matrix Q2 used to transform the
C             matrices S and H.
C             On exit, if COMPQ2 = 'I', then the leading N-by-N part of
C             this array contains the orthogonal transformation matrix
C             Q2.
C             If COMPQ2 = 'N', this array is not referenced.
C
C     LDQ2    INTEGER
C             The leading dimension of the array Q2.
C             LDQ2 >= 1,         if COMPQ2 = 'N';
C             LDQ2 >= MAX(1, N), if COMPQ2 = 'I' or COMPQ2 = 'U'.
C
C     B       (output) DOUBLE PRECISION array, dimension (LDB, N/2)
C             On exit, if JOB = 'T', the leading N/2-by-N/2 part of this
C             array contains the matrix Bout; otherwise, it contains the
C             upper triangular matrix B obtained just before the
C             application of the periodic QZ algorithm.
C
C     LDB     INTEGER
C             The leading dimension of the array B.  LDB >= MAX(1, N/2).
C
C     F       (output) DOUBLE PRECISION array, dimension (LDF, N/2)
C             On exit, if JOB = 'T', the leading N/2-by-N/2 strictly
C             upper triangular part of this array contains the strictly
C             upper triangular part of the skew-symmetric matrix Fout.
C             If JOB = 'E', the leading N/2-by-N/2 strictly upper
C             triangular part of this array contains the strictly upper
C             triangular part of the skew-symmetric matrix F just before
C             the application of the periodic QZ algorithm.
C             The entries on the leading N/2-by-N/2 lower triangular
C             part of this array are not referenced.
C
C     LDF     INTEGER
C             The leading dimension of the array F.  LDF >= MAX(1, N/2).
C
C     C2      (output) DOUBLE PRECISION array, dimension (LDC2, N/2)
C             On exit, if JOB = 'T', the leading N/2-by-N/2 part of this
C             array contains the matrix C2out; otherwise, it contains
C             the upper Hessenberg matrix C2 obtained just before the
C             application of the periodic QZ algorithm.
C
C     LDC2    INTEGER
C             The leading dimension of the array C2.
C             LDC2 >= MAX(1, N/2).
C
C     ALPHAR  (output) DOUBLE PRECISION array, dimension (N/2)
C             The real parts of each scalar alpha defining an eigenvalue
C             of the pencil aS - bH.
C
C     ALPHAI  (output) DOUBLE PRECISION array, dimension (N/2)
C             The imaginary parts of each scalar alpha defining an
C             eigenvalue of the pencil aS - bH.
C             If ALPHAI(j) is zero, then the j-th eigenvalue is real.
C
C     BETA    (output) DOUBLE PRECISION array, dimension (N/2)
C             The scalars beta that define the eigenvalues of the pencil
C             aS - bH.
C             Together, the quantities alpha = (ALPHAR(j),ALPHAI(j)) and
C             beta = BETA(j) represent the j-th eigenvalue of the pencil
C             aS - bH, in the form lambda = alpha/beta. Since lambda may
C             overflow, the ratios should not, in general, be computed.
C             Due to the skew-Hamiltonian/Hamiltonian structure of the
C             pencil, for every eigenvalue lambda, -lambda is also an
C             eigenvalue, and thus it has only to be saved once in
C             ALPHAR, ALPHAI and BETA.
C             Specifically, only eigenvalues with imaginary parts
C             greater than or equal to zero are stored; their conjugate
C             eigenvalues are not stored. If imaginary parts are zero
C             (i.e., for real eigenvalues), only positive eigenvalues
C             are stored. The remaining eigenvalues have opposite signs.
C             As a consequence, pairs of complex eigenvalues, stored in
C             consecutive locations, are not complex conjugate.
C
C     Workspace
C
C     IWORK   INTEGER array, dimension (LIWORK)
C             On exit, if INFO = 3, IWORK(1) contains the number of
C             possibly inaccurate eigenvalues, q <= N/2, and IWORK(2),
C             ..., IWORK(q+1) indicate their indices. Specifically, a
C             positive value is an index of a real or purely imaginary
C             eigenvalue, corresponding to the 1-by-1 blocks, while the
C             absolute value of a negative entry in IWORK is an index to
C             the first eigenvalue in a pair of consecutively stored
C             eigenvalues, corresponding to the 2-by-2 blocks.
C             IWORK(q+2) contains the number s of finite eigenvalues
C             corresponding to the 1-by-1 blocks, and IWORK(q+3)
C             contains the number t of the 2-by-2 blocks. A 2-by-2 block
C             may have two complex, two real, two purely imaginary, or
C             one real and one purely imaginary eigenvalues.
C             If INFO = 0, then q = 0.
C
C     LIWORK  INTEGER
C             The dimension of the array IWORK.
C             LIWORK >= N+12.
C
C     DWORK   DOUBLE PRECISION array, dimension (LDWORK)
C             On exit, if INFO = 0 or INFO = 3, DWORK(1) returns the
C             optimal LDWORK, and DWORK(2), ..., DWORK(5) contain the
C             Frobenius norms of the factors of the formal matrix
C             product used by the algorithm; in addition, DWORK(6), ...,
C             DWORK(5+4*s) contain the s quadruple values corresponding
C             to the used 1-by-1 blocks. Their eigenvalues are finite
C             real and purely imaginary. (Such an eigenvalue is obtained
C             from -i*sqrt(a1*a3/a2/a4), but always taking a positive
C             sign, where a1, ..., a4 are the corresponding quadruple
C             values.)
C             Moreover, DWORK(6+4*s), ..., DWORK(5+4*s+16*t) contain the
C             t groups of quadruple 2-by-2 matrices corresponding to the
C             used 2-by-2 blocks. Their eigenvalue pairs are finite,
C             either complex, or placed on the real and imaginary axes.
C             (Such an eigenvalue pair is obtained from the eigenvalues
C             of the matrix product A1*inv(A2)*A3*inv(A4), where A1,
C             ..., A4 define the corresponding 2-by-2 matrix quadruple;
C             such matrix products are not evaluated.)
C             On exit, if INFO = -27, DWORK(1) returns the minimum
C             value of LDWORK.
C
C     LDWORK  INTEGER
C             The dimension of the array DWORK.
C             If JOB = 'E' and COMPQ1 = 'N' and COMPQ2 = 'N',
C                LDWORK >= N**2 + MAX(L, 36);
C             if JOB = 'T' or COMPQ1 <> 'N' or COMPQ2 <> 'N',
C                LDWORK >= 2*N**2 + MAX(L, 36);
C             where
C                L = 4*N + 4, if N/2 is even, and
C                L = 4*N    , if N/2 is odd.
C             For good performance LDWORK should generally be larger.
C             The dimension of the array DWORK.
C
C     Error Indicator
C
C     INFO    INTEGER
C             = 0: succesful exit;
C             < 0: if INFO = -i, the i-th argument had an illegal value;
C             = 1: problem during computation of the eigenvalues;
C             = 2: periodic QZ algorithm did not converge in the SLICOT
C                  Library subroutine MB03BD;
C             = 3: some eigenvalues might be inaccurate, and details can
C                  be found in IWORK and DWORK. This is a warning.
C
C     METHOD
C
C     The algorithm uses Givens rotations and Householder reflections to
C     annihilate elements in S, T, and H such that A, B, and C1 are
C     upper triangular and C2 is upper Hessenberg. Finally, the periodic
C     QZ algorithm is applied to transform C2 to upper quasi-triangular
C     form while A, B, and C1 stay in upper triangular form.
C     See also page 27 in [1] for more details.
C
C     REFERENCES
C
C     [1] Benner, P., Byers, R., Losse, P., Mehrmann, V. and Xu, H.
C         Numerical Solution of Real Skew-Hamiltonian/Hamiltonian
C         Eigenproblems.
C         Tech. Rep., Technical University Chemnitz, Germany,
C         Nov. 2007.
C
C     NUMERICAL ASPECTS
C                                                               3
C     The algorithm is numerically backward stable and needs O(N ) real
C     floating point operations.
C
C     FURTHER COMMENTS
C
C     For large values of N, the routine applies the transformations
C     for reducing T on panels of columns. The user may specify in INFO
C     the desired number of columns. If on entry INFO <= 0, then the
C     routine estimates a suitable value of this number.
C
C     CONTRIBUTOR
C
C     V. Sima, Research Institute for Informatics, Bucharest, Jul. 2011
C
C     REVISIONS
C
C     M. Voigt, Jan. 2012, July 2013, June 2014, July 2014.
C     V. Sima, Oct. 2012, Jan. 2013, Feb. 2013, July 2013, July 2014.
C
C     KEYWORDS
C
C     periodic QZ algorithm, upper (quasi-)triangular matrix,
C     skew-Hamiltonian/Hamiltonian pencil.
C
C     ******************************************************************
C
C     .. Parameters ..
      DOUBLE PRECISION   ZERO, HALF, ONE, TWO
      PARAMETER          ( ZERO = 0.0D+0, HALF = 0.5D+0, ONE = 1.0D+0,
     $                     TWO  = 2.0D+0 )
      INTEGER            DKIND
      PARAMETER          ( DKIND = KIND( ZERO ) )
C
C     .. Scalar Arguments ..
      CHARACTER          COMPQ1, COMPQ2, JOB
      INTEGER            INFO, LDA, LDB, LDC1, LDC2, LDDE, LDF, LDQ1,
     $                   LDQ2, LDVW, LDWORK, LIWORK, N
C
C     .. Array Arguments ..
      INTEGER            IWORK( * )
      DOUBLE PRECISION   A( LDA, * ), ALPHAI( * ), ALPHAR( * ),
     $                   B( LDB, * ), BETA( * ), C1( LDC1, * ),
     $                   C2( LDC2, * ), DE( LDDE, * ), DWORK( * ),
     $                   F( LDF, * ), Q1( LDQ1, * ), Q2( LDQ2, * ),
     $                   VW( LDVW, * )
C
C     .. Local Scalars ..
      LOGICAL            I2R, LCMPQ1, LCMPQ2, LINIQ1, LINIQ2, LTRI,
     $                   LUPDQ1, LUPDQ2
      CHARACTER*16       CMPQ, CMPSC
      INTEGER            I, I11, I22, IMAT, IWARN, IWRK, J, K, M, M1,
     $                   MJ2, MJ3, MK1, MK2, MK3, MM, OPTDW
      DOUBLE PRECISION   BASE, CO, EMAX, EMIN, MU, NU, PREC, SI, TMP1,
     $                   TMP2
      COMPLEX(DKIND) ::  EIG
C
C     .. Local Arrays ..
      INTEGER            IDUM( 1 )
      DOUBLE PRECISION   DUM( 1 )
C
C     .. External Functions ..
      LOGICAL            LSAME
      DOUBLE PRECISION   DDOT, DLAMCH, DLAPY2
      EXTERNAL           DDOT, DLAMCH, DLAPY2, LSAME
C
C     .. External Subroutines ..
      EXTERNAL           DAXPY, DCOPY, DGEMM, DGEQRF, DLACPY, DLARF,
     $                   DLARFG, DLARTG, DLASET, DROT, DSYMV, DSYR2,
     $                   MA02AD, MB01LD, MB01MD, MB01ND, MB03BD, XERBLA
C
C     .. Intrinsic Functions ..
      INTRINSIC          ABS, AIMAG, CMPLX, DBLE, INT, MAX, MIN, MOD,
     $                   SQRT
C
C *** Variables for blocked code.
C
      INTEGER            IC, ICS, JA, JB, JC, JE, JS, N1, N2, N3, NB, NC
      DOUBLE PRECISION,  ALLOCATABLE :: R( : )
C
C ***
C     .. Executable Statements ..
C
C     Decode the input arguments.
C
C ***
      NB = INFO
C ***
      M  = N/2
      MM = M*M
      M1 = MAX( 1, M )
      LTRI   = LSAME( JOB,    'T' )
      LINIQ1 = LSAME( COMPQ1, 'I' )
      LINIQ2 = LSAME( COMPQ2, 'I' )
      LUPDQ1 = LSAME( COMPQ1, 'U' )
      LUPDQ2 = LSAME( COMPQ2, 'U' )
      LCMPQ1 = LUPDQ1 .OR. LINIQ1
      LCMPQ2 = LUPDQ2 .OR. LINIQ2
C
C     Test the input arguments.
C
      INFO = 0
      IF( .NOT.( LSAME( JOB, 'E' ) .OR. LTRI ) ) THEN
         INFO = -1
      ELSE IF( .NOT.( LSAME( COMPQ1, 'N' ) .OR. LCMPQ1 ) ) THEN
         INFO = -2
      ELSE IF( .NOT.( LSAME( COMPQ2, 'N' ) .OR. LCMPQ2 ) ) THEN
         INFO = -3
      ELSE IF( ( LINIQ2 .AND. .NOT.LINIQ1 ) .OR.
     $         ( LUPDQ2 .AND. .NOT.LUPDQ1 ) ) THEN
         INFO = -3
      ELSE IF( N.LT.0 .OR. MOD( N, 2 ).NE.0 ) THEN
         INFO = -4
      ELSE IF(  LDA.LT.M1 ) THEN
         INFO = -6
      ELSE IF( LDDE.LT.M1 ) THEN
         INFO = -8
      ELSE IF( LDC1.LT.M1 ) THEN
         INFO = -10
      ELSE IF( LDVW.LT.M1 ) THEN
         INFO = -12
      ELSE IF( LDQ1.LT.1 .OR. ( LCMPQ1 .AND. LDQ1.LT.N ) ) THEN
         INFO = -14
      ELSE IF( LDQ2.LT.1 .OR. ( LCMPQ2 .AND. LDQ2.LT.N ) ) THEN
         INFO = -16
      ELSE IF(  LDB.LT.M1 ) THEN
         INFO = -18
      ELSE IF(  LDF.LT.M1 ) THEN
         INFO = -20
      ELSE IF( LDC2.LT.M1 ) THEN
         INFO = -22
      ELSE IF( LIWORK.LT.N+12 ) THEN
         INFO = -27
      ELSE
         IF( MOD( M, 2 ).EQ.0 ) THEN
            I = MAX( 4*N, 32 ) + 4
         ELSE
            I = MAX( 4*N, 36 )
         END IF
         IF( LTRI .OR. LCMPQ1 .OR. LCMPQ2 ) THEN
            OPTDW = 8*MM + I
         ELSE
            OPTDW = 4*MM + I
         END IF
         IF( LDWORK.LT.OPTDW )
     $      INFO = -29
      END IF
C
      IF( INFO.NE.0) THEN
         CALL XERBLA( 'DGHUTP', -INFO )
         RETURN
      END IF
C
C     Quick return if possible.
C
      IF( N.EQ.0 ) THEN
         DWORK( 1 ) = ONE
         RETURN
      END IF
C
C     A block algorithm is used for large M.
C
      IF( NB.LE.0 ) THEN
         CALL DGEQRF( M, M, A, LDA, DWORK, DWORK, -1, INFO )
         NB = MIN( MAX( INT( DWORK( 1 ) )/M1, 2 ), M )
      END IF
C
C     Determine machine constants.
C
      BASE = DLAMCH( 'Base' )
      EMIN = DLAMCH( 'Minimum Exponent' )
      EMAX = DLAMCH( 'Largest Exponent' )
      PREC = DLAMCH( 'Precision' )*DBLE( N )
C
C     STEP 1: Reduce S to skew-Hamiltonian triangular form.
C
      IF( LINIQ1 )
     $   CALL DLASET( 'Full', N, N, ZERO, ONE, Q1, LDQ1 )
C
      DUM( 1 ) = ZERO
C
      DO 10 K = 1, M - 1
C
C        Generate elementary reflector H(k) = I - nu * v * v' to
C        annihilate E(k+2:m,k).
C
         MK2  = MIN( K+2, M )
         MK3  = MK2 + 1
         TMP1 = DE( K+1, K )
         CALL DLARFG( M-K, TMP1, DE( MK2, K ), 1, NU )
         IF( NU.NE.ZERO ) THEN
            DE( K+1, K ) = ONE
C
C           Apply H(k) from both sides to E(k+1:m,k+1:m).
C           Compute  x := nu * E(k+1:m,k+1:m) * v.
C
            CALL MB01MD( 'Lower', M-K, NU, DE( K+1, K+1 ), LDDE,
     $                   DE( K+1, K ), 1, ZERO, DWORK, 1 )
C
C           Compute  w := x - 1/2 * nu * (x'*v) * v in x.
C
            MU = -HALF*NU*DDOT( M-K, DWORK, 1, DE( K+1, K ), 1 )
            CALL DAXPY( M-K, MU, DE( K+1, K ), 1, DWORK, 1 )
C
C           Apply the transformation as a skew-symmetric rank-2 update:
C                E := E + v * w' - w * v'.
C
            CALL MB01ND( 'Lower', M-K, ONE, DE( K+1, K ), 1, DWORK, 1,
     $                   DE( K+1, K+1 ), LDDE )
C
C           Apply H(k) to W(k+1:m,1:k) from the left (and implicitly to
C           W(1:k,k+1:m) from the right).
C
            CALL DLARF( 'Left', M-K, K, DE( K+1, K ), 1, NU,
     $                  VW( K+1, 1 ), LDVW, DWORK )
C
C           Apply H(k) from both sides to W(k+1:m,k+1:m).
C           Compute  x := nu * W(k+1:m,k+1:m) * v.
C
            CALL DSYMV( 'Lower', M-K, NU, VW( K+1, K+1 ), LDVW,
     $                  DE( K+1, K ), 1, ZERO, DWORK, 1 )
C
C           Compute  w := x - 1/2 * nu * (x'*v) * v.
C
            MU = -HALF*NU*DDOT( M-K, DWORK, 1, DE( K+1, K ), 1 )
            CALL DAXPY( M-K, MU, DE( K+1, K ), 1, DWORK, 1 )
C
C           Apply the transformation as a rank-2 update:
C                W := W - v * w' - w * v'.
C
            CALL DSYR2( 'Lower', M-K, -ONE, DE( K+1, K ), 1, DWORK, 1,
     $                  VW( K+1, K+1 ), LDVW )
C
C           Apply H(k) from the right hand side to A(1:m,k+1:m) and
C           C1(1:m,k+1:m).
C
            CALL DLARF( 'Right', M, M-K, DE( K+1, K ), 1, NU,
     $                  A( 1, K+1 ), LDA, DWORK )
            CALL DLARF( 'Right', M, M-K, DE( K+1, K ), 1, NU,
     $                  C1( 1, K+1 ), LDC1, DWORK )
C
            IF( LCMPQ1 ) THEN
C
C              Apply H(k) from the right hand side to Q1(1:n,m+k+1:n).
C
               CALL DLARF( 'Right', N, M-K, DE( K+1, K ), 1, NU,
     $                     Q1( 1, M+K+1 ), LDQ1, DWORK )
            END IF
            DE( K+1, K ) = TMP1
         END IF
C
C        Determine a Givens rotation to annihilate E(k+1,k) from the
C        left.
C
         TMP2 = A( K+1, K )
         CALL DLARTG( TMP2, DE( K+1, K ), CO, SI, A( K+1, K ) )
C
C        Update A, E and D.
C
         CALL DROT( M-K-1, DE( MK2, K+1 ), 1, A( K+1, MK2 ), LDA, CO,
     $              SI )
         CALL DROT( K, A( 1, K+1 ), 1, DE( 1, K+2 ), 1, CO, SI )
         CALL DROT( M-K-1, DE( K+1, MK3 ), LDDE, A( MK2, K+1 ), 1, CO,
     $              SI )
C
C        Update C1, W and V.
C
         CALL DROT( K, VW( K+1, 1 ), LDVW, C1( K+1, 1 ), LDC1, CO, -SI )
         CALL DROT( M-K-1, VW( MK2, K+1 ), 1, C1( K+1, MK2 ), LDC1, CO,
     $              -SI )
         CALL DROT( K, C1( 1, K+1 ), 1, VW( 1, K+2 ), 1, CO, SI )
         CALL DROT( M-K-1, VW( K+1, MK3 ), LDVW, C1( MK2, K+1 ), 1, CO,
     $              -SI )
C
C        Fix the diagonal part.
C
         TMP1 = C1( K+1, K+1 )
         TMP2 = VW( K+1, K+2 )
         C1( K+1, K+1 ) = ( CO - SI )*( CO + SI )*TMP1 +
     $                    CO*SI*( VW( K+1, K+1 ) + TMP2 )
         TMP1 = TWO*CO*SI*TMP1
         VW( K+1, K+2 ) = CO**2*TMP2 - SI**2*VW( K+1, K+1 ) - TMP1
         VW( K+1, K+1 ) = CO**2*VW( K+1, K+1 ) - SI**2*TMP2 - TMP1
C
         IF( LCMPQ1 ) THEN
C
C           Update Q1.
C
            CALL DROT( N, Q1( 1, K+1 ), 1, Q1( 1, M+K+1 ), 1, CO, SI )
         END IF
C
C        Generate elementary reflector P(k) to annihilate A(k+1:m,k).
C
         TMP1 = A( K, K )
         CALL DLARFG( M-K+1, TMP1, A( K+1, K ), 1, NU )
         IF( NU.NE.ZERO ) THEN
            A( K, K ) = ONE
C
C           Apply P(k) from the left hand side to A(k:m,k+1:m).
C
            CALL DLARF( 'Left', M-K+1, M-K, A( K, K ), 1, NU,
     $                  A( K, K+1 ), LDA, DWORK )
C
C           Apply P(k) to D(1:k-1,k:m) from the right (and implicitly
C           to D(k:m,1:k-1) from the left).
C
            CALL DLARF( 'Right', K-1, M-K+1, A( K, K ), 1, NU,
     $                  DE( 1, K+1 ), LDDE, DWORK )
C
C           Apply P(k) from both sides to D(k:m,k:m).
C           Compute  x := nu * D(k:m,k:m) * v.
C
            CALL MB01MD( 'Upper', M-K+1, NU, DE( K, K+1 ), LDDE,
     $                   A( K, K ), 1, ZERO, DWORK, 1 )
C
C           Compute  w := x - 1/2 * nu * (x'*v) * v in x.
C
            MU = -HALF*NU*DDOT( M-K+1, DWORK, 1, A( K, K ), 1 )
            CALL DAXPY( M-K+1, MU, A( K, K ), 1, DWORK, 1 )
C
C           Apply the transformation as a skew-symmetric rank-2 update:
C                D := D + v * w' - w * v'.
C
            CALL MB01ND( 'Upper', M-K+1, ONE, A( K, K ), 1, DWORK, 1,
     $                   DE( K, K+1 ), LDDE )
C
C           Apply P(k) from the left hand side to C1(k:m,1:m).
C
            CALL DLARF( 'Left', M-K+1, M, A( K, K ), 1, NU, C1( K, 1 ),
     $                  LDC1, DWORK )
C
C           Apply P(k) to V(1:k-1,k:m) from the right (and implicitly
C           to V(k:m,1:k-1) from the left).
C
            CALL DLARF( 'Right', K-1, M-K+1, A( K, K ), 1, NU,
     $                  VW( 1, K+1 ), LDVW, DWORK )
C
C           Apply P(k) from both sides to V(k:m,k:m).
C           Compute  x := nu * V(k:m,k:m) * v.
C
            CALL DSYMV( 'Upper', M-K+1, NU, VW( K, K+1 ), LDVW,
     $                  A( K, K ), 1, ZERO, DWORK, 1 )
C
C           Compute  w := x - 1/2 * nu * (x'*v) * v.
C
            MU = -HALF*NU*DDOT( M-K+1, DWORK, 1, A( K, K ), 1 )
            CALL DAXPY( M-K+1, MU, A( K, K ), 1, DWORK, 1 )
C
C           Apply the transformation as a rank-2 update:
C                V := V - v * w' - w * v'.
C
            CALL DSYR2( 'Upper', M-K+1, -ONE, A( K, K ), 1, DWORK, 1,
     $                  VW( K, K+1 ), LDVW )
C
            IF( LCMPQ1 ) THEN
C
C              Apply P(k) from the right hand side to Q1(1:n,k:m).
C
               CALL DLARF( 'Right', N, M-K+1, A( K, K ), 1, NU,
     $                     Q1( 1, K ), LDQ1, DWORK )
            END IF
            A( K, K ) = TMP1
         END IF
C
C        Set A(k+1:m,k) to zero in order to be able to apply MB03BD.
C
         CALL DCOPY( M-K, DUM, 0, A( K+1, K ), 1 )
   10 CONTINUE
C
C     The following operations do not preserve the Hamiltonian structure
C     of H. -C1 is copied to C2. The lower triangular part of W(1:m,1:m)
C     and its transpose are stored in DWORK. Then, the transpose of the
C     upper triangular part of V(1:m,1:m) is saved in the lower
C     triangular part of VW(1:m,2:m+1).
C
      CALL DLACPY( 'Full', M, M, A, LDA, B, LDB )
      CALL DLACPY( 'Upper', M, M, DE( 1, 2 ), LDDE, F, LDF )
C
      DO 30 J = 1, M
         DO 20 I = 1, M
            C2( I, J ) = -C1( I, J )
   20    CONTINUE
   30 CONTINUE
C
      CALL DLACPY( 'Lower', M, M, VW, LDVW, DWORK, M )
      CALL MA02AD( 'Lower', M, M, VW, LDVW, DWORK, M )
C
      CALL MA02AD( 'Upper', M, M, VW( 1, 2 ), LDVW, VW( 1, 2 ), LDVW )
C
      IF ( LCMPQ2 ) THEN
         CALL DLACPY( 'Full', M, M, Q1( M+1, M+1 ), LDQ1, Q2, LDQ2 )
C
         DO 50 J = 1, M
            DO 40 I = M + 1, N
               Q2( I, J ) = -Q1( I-M, J+M )
   40       CONTINUE
   50    CONTINUE
C
         DO 70 J = M + 1, N
            DO 60 I = 1, M
               Q2( I, J ) = -Q1( I+M, J-M )
   60       CONTINUE
   70    CONTINUE
C
         CALL DLACPY( 'Full', M, M, Q1, LDQ1, Q2( M+1, M+1 ), LDQ2 )
      END IF
C      
C     Allocate space for storing rotations.
C
      ALLOCATE( R( 2*N ) )
C
C     STEP 2: Eliminations in H.
C
      DO 390 K = 1, M
         MK1 = MIN( K+1, M )
C
C        I. Annihilate W(k:m-1,k).
C
         JS = K
         JE = MIN( M, JS+NB-1 )
         JB = JE
         IC = 1
         JC = 2*( M - K ) + 1
         N1 = MIN( M-K,   JE-JS+1 )
         N2 = MIN( M-K+1, JE-JS+1 )
         N3 = MIN( M, NB )
C
         DO 100 J = K, M - 1
            MJ3  = MIN( J+3, M+1 )
C
C           Determine a Givens rotation to annihilate W(j,k) from the
C           left.
C
            CALL DLARTG( DWORK( ( K-1 )*M+J+1 ), DWORK( ( K-1 )*M+J ),
     $                   CO, SI, TMP1 )
            R( IC   ) = CO
            R( IC+1 ) = SI
            IC = IC + 2
C
C           Update C2 and W.
C
            CALL DROT( M, C2( 1, J+1 ), 1, C2( 1, J ), 1, CO, SI )
            DWORK( ( K-1 )*M+J+1 ) = TMP1
            DWORK( ( K-1 )*M+J )   = ZERO
            CALL DROT( N1, DWORK( K*M+J+1 ), M, DWORK( K*M+J ), M, CO,
     $                 SI )
C
C           Update A.
C
            IF( J.EQ.JE ) THEN
C
C              Update the next panel (the columns JE+1 to JE+NB) of A
C              and D for previous row transformations.
C
               JS = JE + 1
               JE = MIN( M, JE+NB )
               NC = JE - JS + 1
               JA = 2*( M - K ) + 1
               DO 80 I = K, J - 1
                  CALL DROT( NC, A( I, JS ), LDA, A( I+1, JS ), LDA,
     $                       R( JA ), R( JA+1 ) )
                  JA = JA + 2
   80          CONTINUE
               JA = 2*( M - K ) + 1
               DO 90 I = K, J - 1
                  CALL DROT( NC, DE( I, JS+1 ), LDDE, DE( I+1, JS+1 ),
     $                       LDDE, R( JA ), R( JA+1 ) )
                  JA = JA + 2
   90          CONTINUE
            END IF
C
            CALL DROT( J, A( 1, J+1 ), 1, A( 1, J ), 1, CO, SI )
            TMP1          = -SI*A( J+1, J+1 )
            A( J+1, J+1 ) =  CO*A( J+1, J+1 )
C
            IF( LCMPQ1 ) THEN
C
C              Update Q1.
C
               CALL DROT( N, Q1( 1, M+J+1 ), 1, Q1( 1, M+J ), 1, CO, SI
     $                     )
            END IF
C
C           Determine a Givens rotation to annihilate A(j+1,j) from the
C           left.
C
            CALL DLARTG( A( J, J ), TMP1, CO, SI, TMP2 )
            R( JC   ) = CO
            R( JC+1 ) = SI
            JC = JC + 2
C
C           Update A and D.
C
            NC = MIN( JE-J, JE-JS+1  )
            A( J, J ) = TMP2
            CALL DROT( NC, A( J, J+1 ), LDA, A( J+1, J+1 ), LDA, CO, SI
     $                  )
            CALL DROT( J-1, DE( 1, J+1 ), 1, DE( 1, J+2 ), 1, CO, SI )
            CALL DROT( NC-1, DE( J, MJ3 ), LDDE, DE( J+1, MJ3 ), LDDE,
     $                 CO, SI )
C
C           Update C1 and V.
C
            CALL DROT( N2, C1( J, K ), LDC1, C1( J+1, K ), LDC1, CO,
     $                 SI )
            CALL DROT( N3, VW( J, 2 ), LDVW, VW( J+1, 2 ), LDVW, CO,
     $                 SI )
C
            IF( LCMPQ1 ) THEN
C
C              Update Q1.
C
               CALL DROT( N, Q1( 1, J ), 1, Q1( 1, J+1 ), 1, CO, SI )
            END IF
  100    CONTINUE
C
C        Update the remaining panels of columns of W, C1, and V for
C        previous row transformations.
C
         DO 130 JS = JB + 1, M, NB
            JE = MIN( M, JS+NB-1 )
            NC = JE - JS + 1
            IC = 1
            DO 110 J = K, M - 1
               CALL DROT( NC, DWORK( JS*M+J+1 ), M, DWORK( JS*M+J ), M,
     $                    R( IC ), R( IC+1 ) )
               IC = IC + 2
  110       CONTINUE
            JC = 2*( M - K ) + 1
            DO 120 J = K, M - 1
               CALL DROT( NC, C1( J, JS ), LDC1, C1( J+1, JS ), LDC1,
     $                    R( JC ), R( JC+1 ) )
               JC = JC + 2
  120       CONTINUE
  130    CONTINUE
C
         DO 150 JS = 2 + NB, M + 1, NB
            JE = MIN( M+1, JS+NB-1 )
            NC = JE - JS + 1
            JC = 2*( M - K ) + 1
            DO 140 J = K, M - 1
               CALL DROT( NC, VW( J, JS ), LDVW, VW( J+1, JS ), LDVW,
     $                    R( JC ), R( JC+1 ) )
               JC = JC + 2
  140       CONTINUE
  150    CONTINUE
C
C        II. Annihilate W(m,k).
C
C        Determine a Givens rotation to annihilate W(m,k) from the left.
C
         CALL DLARTG( C1( M, K ), DWORK( M*K ), CO, SI, TMP1 )
C
C        Update C1 and W.
C
         C1( M, K )   = TMP1
         DWORK( M*K ) = ZERO
         CALL DROT( M-K, C1( M, MK1 ), LDC1, DWORK( M*MK1 ), M, CO, SI )
         CALL DROT( M, VW( M, 2 ), LDVW, C2( 1, M ), 1, CO, SI )
C
C        Update A and D.
C
         CALL DROT( M-1, A( 1, M ), 1, DE( 1, M+1 ), 1, CO, SI )
C
         IF( LCMPQ1 ) THEN
C
C           Update Q1.
C
            CALL DROT( N, Q1( 1, M ), 1, Q1( 1, N ), 1, CO, SI )
         END IF
C
C        III. Annihilate C1(k+1:m,k).
C
         JS = K
         JE = MIN( M, JS+NB-1 )
         IC = 1
         JC = 2*( M - K ) + 1
         N1 = MIN( M-K,   JE-JS+1 )
         N2 = MIN( M-K+1, JE-JS+1 )
         N3 = MIN( M, NB )
C
         DO 160 J = M, K + 1, -1
C
C           Determine a Givens rotation to annihilate C1(j,k) from the
C           left.
C
            CALL DLARTG( C1( J-1, K ), C1( J, K ), CO, SI, TMP1 )
            R( IC   ) = CO
            R( IC+1 ) = SI
            IC = IC + 2
C
C           Update C1 and V.
C
            C1( J-1, K ) = TMP1
            C1( J, K )   = ZERO
            CALL DROT( N1, C1( J-1, MK1 ), LDC1, C1( J, MK1 ), LDC1,
     $                 CO, SI )
            CALL DROT( N3, VW( J-1, 2 ), LDVW, VW( J, 2 ), LDVW , CO, SI
     $                  )
C
C           Update A and D.
C
            TMP1          = -SI*A( J-1, J-1 )
            A( J-1, J-1 ) =  CO*A( J-1, J-1 )
            CALL DROT( 1, A( J-1, J ), LDA, A( J, J ), LDA, CO, SI )
            CALL DROT( J-2, DE( 1, J ), 1, DE( 1, J+1 ), 1, CO, SI )
C
            IF( LCMPQ1 ) THEN
C
C              Update Q1.
C
               CALL DROT( N, Q1( 1, J-1 ), 1, Q1( 1, J ), 1, CO, SI )
            END IF
C
C           Determine a Givens rotation to annihilate A(j,j-1) from the
C           right.
C
            CALL DLARTG( A( J, J ), TMP1, CO, SI, TMP2 )
            R( JC   ) = CO
            R( JC+1 ) = SI
            JC = JC + 2
C
C           Update A.
C
            A( J, J ) = TMP2
            CALL DROT( J-1, A( 1, J ), 1, A( 1, J-1 ), 1, CO, SI )
C
C           Update C2 and W.
C
            CALL DROT( M, C2( 1, J ), 1, C2( 1, J-1 ), 1, CO, SI )
            CALL DROT( N2, DWORK( ( K-1 )*M+J ), M,
     $                 DWORK( ( K-1)*M+J-1 ), M, CO, SI )
C
            IF( LCMPQ1 ) THEN
C
C              Update Q1.
C
               CALL DROT( N, Q1( 1, M+J ), 1, Q1( 1, M+J-1 ), 1, CO, SI
     $                     )
            END IF
  160    CONTINUE
C
C        Update the remaining panels of columns of C1, W, V, A, and D
C        for previous row transformations.
C
         DO 180 JS = MK1 + N1, M, NB
            JE = MIN( M, JS+NB-1 )
            NC = JE - JS + 1
            IC = 1
            DO 170 J = M, K + 1, -1
               CALL DROT( NC, C1( J-1, JS ), LDC1, C1( J, JS ), LDC1,
     $                    R( IC ), R( IC+1 ) )
               IC = IC + 2
  170       CONTINUE
  180    CONTINUE
C
         DO 200 JS = K - 1 + N2, M, NB
            JE = MIN( M, JS+NB-1 )
            NC = JE - JS + 1
            JC = 2*( M - K ) + 1
            DO 190 J = M, K + 1, -1
               CALL DROT( NC, DWORK( JS*M+J ), M, DWORK(JS*M+J-1 ), M,
     $                    R( JC ), R( JC+1 ) )
               JC = JC + 2
  190       CONTINUE
  200    CONTINUE
C
         DO 220 JS = 2 + NB, M + 1, NB
            JE = MIN( M+1, JS+NB-1 )
            NC = JE - JS + 1
            IC = 1
            DO 210 J = M, K + 1, -1
               CALL DROT( NC, VW( J-1, JS ), LDVW, VW( J, JS ), LDVW,
     $                    R( IC ), R( IC+1 ) )
               IC = IC + 2
  210       CONTINUE
  220    CONTINUE
C
         ICS = 3
         JE  = M
C
C        WHILE( JE.GT.2 ) DO
  230    CONTINUE
         IF( JE.GT.2 ) THEN
            NC  = 0
            IC  = ICS
            ICS = ICS + 2*NB
            DO 240 J = JE - 1, K + 1, -1
               NC = MIN( NC+1, NB )
               JS = JE - NC + 1
               CALL DROT( NC, A( J-1, JS ), LDA, A( J, JS ), LDA,
     $                    R( IC ), R( IC+1 ) )
               IC = IC + 2
  240       CONTINUE
            JE = JE - NB
            GO TO 230
         END IF
C        END WHILE 230
C
         ICS = 3
         JE  = M + 1
C
C        WHILE( JE.GT.3 ) DO
  250    CONTINUE
         IF( JE.GT.3 ) THEN
            NC  = 0
            IC  = ICS
            ICS = ICS + 2*NB
            DO 260 J = JE - 2, K + 1, -1
               NC = MIN( NC+1, NB )
               JS = JE - NC + 1
               CALL DROT( NC, DE( J-1, JS ), LDDE, DE( J, JS ), LDDE,
     $                    R( IC ), R( IC+1 ) )
               IC = IC + 2
  260       CONTINUE
            JE = JE - NB
            GO TO 250
         END IF
C        END WHILE 250
C
C        IV. Annihilate W(k,k+1:m-1).
C
         JS = K + 1
         JE = MIN( M, JS+NB-1 )
         JC = 1
         N2 = MIN( M-K+1, JE-JS+1 )
C
         DO 290 J = K + 1, M - 1
            MJ2  = MIN( J+2, M )
C
C           Determine a Givens rotation to annihilate W(k,j) from the
C           right.
C
            CALL DLARTG( DWORK( J*M+K ), DWORK( ( J-1 )*M+K ), CO,
     $                   SI, TMP1 )
C
C           Update C1 and W.
C
            CALL DROT( M, C1( 1, J+1 ), 1, C1( 1, J ), 1, CO, SI )
            DWORK( ( J-1 )*M+K ) = ZERO
            DWORK( J*M+K )       = TMP1
            CALL DROT( M-K, DWORK( J*M+MK1 ), 1, DWORK( ( J-1 )*M+K+1 ),
     $                 1, CO, SI )
C
C           Update B.
C
            IF( J.EQ.JE ) THEN
C
C              Update the columns JE+1 to JE+NB of B and F for previous
C              row transformations.
C
               JS = JE + 1
               JE = MIN( M, JE+NB )
               NC = JE - JS + 1
               JA = 1
               DO 270 I = K + 1, J - 1
                  CALL DROT( NC, B( I, JS ), LDB, B( I+1, JS ), LDB,
     $                       R( JA ), R( JA+1 ) )
                  JA = JA + 2
  270          CONTINUE
               JA = 1
               DO 280 I = K + 1, J - 1
                  CALL DROT( NC, F( I, JS ), LDF, F( I+1, JS ), LDF,
     $                       R( JA ), R( JA+1 ) )
                  JA = JA + 2
  280          CONTINUE
            END IF
C
            CALL DROT( J, B( 1, J+1 ), 1, B( 1, J ), 1, CO, SI )
            TMP1          = -SI*B( J+1, J+1 )
            B( J+1, J+1 ) =  CO*B( J+1, J+1 )
C
            IF( LCMPQ2 ) THEN
C
C              Update Q2.
C
               CALL DROT( N, Q2( 1, J+1 ), 1, Q2( 1, J ), 1, CO, SI )
            END IF
C
C           Determine a Givens rotation to annihilate B(j+1,j) from the
C           left.
C
            CALL DLARTG( B( J, J ), TMP1, CO, SI, TMP2 )
            R( JC   ) = CO
            R( JC+1 ) = SI
            JC = JC + 2
C
C           Update B and F.
C
            NC = MIN( JE-J, JE-JS+1  )
            B( J, J ) = TMP2
            CALL DROT( NC, B( J, J+1 ), LDB, B( J+1, J+1 ), LDB, CO, SI
     $                  )
            CALL DROT( J-1, F( 1, J ), 1, F( 1, J+1 ), 1, CO, SI )
            CALL DROT( NC-1, F( J, MJ2 ), LDF, F( J+1, MJ2 ), LDF, CO,
     $                 SI )
C
C           Update C2 and V.
C
            CALL DROT( N2, C2( J, K ), LDC2, C2( J+1, K ), LDC2, CO,
     $                 SI )
            CALL DROT( M, VW( 1, J+1 ), 1, VW( 1, J+2 ), 1, CO, SI )
C
            IF( LCMPQ2 ) THEN
C
C              Update Q2.
C
               CALL DROT( N, Q2( 1, M+J ), 1, Q2( 1, M+J+1 ), 1, CO, SI
     $                     )
            END IF
  290    CONTINUE
C
C        Update the remaining panels of columns of C2 for previous row
C        transformations.
C
         DO 310 JS = K + N2, M, NB
            JE = MIN( M, JS+NB-1 )
            NC = JE - JS + 1
            JC = 1
            DO 300 J = K + 1, M - 1
               CALL DROT( NC, C2( J, JS ), LDC2, C2( J+1, JS ), LDC2,
     $                    R( JC ), R( JC+1 ) )
               JC = JC + 2
  300       CONTINUE
  310    CONTINUE
C
C        V. Annihilate W(k,m).
C
         IF( K.LT.M ) THEN
C
C           Determine a Givens rotation to annihilate W(k,m) from the
C           right.
C
            CALL DLARTG( C2( M, K ), DWORK( ( M-1 )*M+K ), CO, SI, TMP1
     $                    )
C
C           Update C1, C2, W and V.
C
            CALL DROT( M, VW( 1, M+1 ), 1, C1( 1, M ), 1, CO, SI )
            C2( M, K )           = TMP1
            DWORK( ( M-1 )*M+K ) = ZERO
            CALL DROT( M-K, C2( M, K+1 ), LDC2, DWORK( ( M-1 )*M+K+1 ),
     $                 1, CO, SI )
C
C           Update B and F.
C
            CALL DROT( M-1, F( 1, M ), 1, B( 1, M ), 1, CO, SI )
C
            IF( LCMPQ2 ) THEN
C
C               Update Q2.
C
                CALL DROT( N, Q2( 1, N ), 1, Q2( 1, M ), 1, CO, SI )
            END IF
         ELSE
C
C           Determine a Givens rotation to annihilate W(m,m) from the
C           left.
C
            CALL DLARTG( C1( M, M ), DWORK( MM ), CO, SI, TMP1 )
C
C           Update C1, C2, W and V.
C
            C1( M, M )  = TMP1
            DWORK( MM ) = ZERO
            CALL DROT( M, VW( M, 2 ), LDVW, C2( 1, M ), 1, CO, SI )
C
C           Update A and D.
C
            CALL DROT( M-1, A( 1, M ), 1, DE( 1, M+1 ), 1, CO, SI )
C
            IF( LCMPQ1 ) THEN
C
C              Update Q1.
C
               CALL DROT( N, Q1( 1, M ), 1, Q1( 1, N ), 1, CO, SI )
            END IF
         END IF
C
C        VI. Annihilate C2(k+2:m,k).
C
         JS = K
         JE = MIN( M, JS+NB-1 )
         IC = 1
         N1 = MIN( M-K, JE-JS+1 )
C
         DO 320 J = M, K + 2, -1
C
C           Determine a Givens rotation to annihilate C2(j,k) from the
C           left.
C
            CALL DLARTG( C2( J-1, K ), C2( J, K ), CO, SI, TMP1 )
            R( IC   ) = CO
            R( IC+1 ) = SI
            IC = IC + 2
C
C           Update C2 and V.
C
            C2( J-1, K ) = TMP1
            C2( J, K )   = ZERO
            CALL DROT( N1, C2( J-1, MK1 ), LDC2, C2( J, MK1 ), LDC2,
     $                 CO, SI )
            CALL DROT( M, VW( 1, J ), 1, VW( 1, J+1 ), 1, CO, SI )
C
C           Update B and F.
C
            CALL DROT( 1, B( J-1, J ), LDB, B( J, J ), LDB, CO, SI )
            TMP1          = -SI*B( J-1, J-1 )
            B( J-1, J-1 ) =  CO*B( J-1, J-1 )
            CALL DROT( J-2, F( 1, J-1 ), 1, F( 1, J ), 1, CO, SI )
C
            IF( LCMPQ2 ) THEN
C
C              Update Q2.
C
               CALL DROT( N, Q2( 1, M+J-1 ), 1, Q2( 1, M+J ), 1, CO, SI
     $                     )
            END IF
C
C           Determine a Givens rotation to annihilate B(j,j-1) from the
C           right.
C
            CALL DLARTG( B( J, J ), TMP1, CO, SI, TMP2 )
            B( J, J ) = TMP2
C
C           Update B.
C
            CALL DROT( J-1, B( 1, J ), 1, B( 1, J-1 ), 1, CO, SI )
C
C           Update C1 and W.
C
            CALL DROT( M, C1( 1, J ), 1, C1( 1, J-1 ), 1, CO, SI )
            CALL DROT( M-K+1, DWORK( ( J-1 )*M+K ), 1,
     $                 DWORK( ( J-2 )*M+K ), 1, CO, SI )
C
            IF( LCMPQ2 ) THEN
C
C              Update Q2.
C
               CALL DROT( N, Q2( 1, J ), 1, Q2( 1, J-1 ), 1, CO, SI )
            END IF
  320    CONTINUE
C
C        Update the remaining panels of columns of C2, B, and F for
C        previous row transformations.
C
         DO 340 JS = MK1 + N1, M, NB
            JE = MIN( M, JS+NB-1 )
            NC = JE - JS + 1
            IC = 1
            DO 330 J = M, K + 2, -1
               CALL DROT( NC, C2( J-1, JS ), LDC2, C2( J, JS ), LDC2,
     $                    R( IC ), R( IC+1 ) )
               IC = IC + 2
  330       CONTINUE
  340    CONTINUE
C
         ICS = 3
         JE  = M
C
C        WHILE( JE.GT.2 ) DO
  350    CONTINUE
         IF( JE.GT.2 ) THEN
            NC  = 0
            IC  = ICS
            ICS = ICS + 2*NB
            DO 360 J = JE - 1, K + 2, -1
               NC = MIN( NC+1, NB )
               JS = JE - NC + 1
               CALL DROT( NC, B( J-1, JS ), LDB, B( J, JS ), LDB,
     $                    R( IC ), R( IC+1 ) )
               IC = IC + 2
  360       CONTINUE
            JE = JE - NB
            GO TO 350
         END IF
C        END WHILE 350
C
         ICS = 3
         JE  = M
C
C        WHILE( JE.GT.3 ) DO
  370    CONTINUE
         IF( JE.GT.3 ) THEN
            NC  = 0
            IC  = ICS
            ICS = ICS + 2*NB
            DO 380 J = JE - 1, K + 2, -1
               NC = MIN( NC+1, NB )
               JS = JE - NC + 1
               CALL DROT( NC, F( J-1, JS ), LDF, F( J, JS ), LDF,
     $                    R( IC ), R( IC+1 ) )
               IC = IC + 2
  380       CONTINUE
            JE = JE - NB
            GO TO 370
         END IF
C        END WHILE 370
C
  390 CONTINUE
C
C     Deallocate the array storing rotations.
C
      DEALLOCATE( R )
C
C                     (  A1  D1  )      (  B1  F1  )      (  C11  V1  )
C     Now we have S = (          ), T = (          ), H = (           ),
C                     (   0  A1' )      (   0  B1' )      (   0  C21' )
C
C     where A1, B1, and C11 are upper triangular, C21 is upper
C     Hessenberg, and D1 and F1 are skew-symmetric.
C
C     STEP 3: Apply the periodic QZ algorithm to the generalized matrix
C
C                           -1       -1
C             product C21 A1   C11 B1   in order to make C21 upper
C             quasi-triangular.
C
C     Determine the mode of computations.
C
      IF( LTRI .OR. LCMPQ1 .OR. LCMPQ2 ) THEN
         CMPQ = 'Initialize'
         IMAT = 4*MM + 1
         IWRK = 8*MM + 1
      ELSE
         CMPQ = 'No Computation'
         IMAT = 1
         IWRK = 4*MM + 1
      END IF
C
      IF( LTRI ) THEN
         CMPSC = 'Schur Form'
      ELSE
         CMPSC = 'Eigenvalues Only'
      END IF
C
C     Save matrices in the form that is required by MB03BD.
C
      CALL DLACPY( 'Full', M, M, C2, LDC2, DWORK( IMAT ),      M )
      CALL DLACPY( 'Full', M, M, A,  LDA,  DWORK( IMAT+MM ),   M )
      CALL DLACPY( 'Full', M, M, C1, LDC1, DWORK( IMAT+2*MM ), M )
      CALL DLACPY( 'Full', M, M, B,  LDB,  DWORK( IMAT+3*MM ), M )
      IWORK( 1 ) =  1
      IWORK( 2 ) = -1
      IWORK( 3 ) =  1
      IWORK( 4 ) = -1
C
C     Apply periodic QZ algorithm.
C     Workspace:    need   IWRK + MAX( N, 32 ) + 3.
C
      CALL MB03BD( CMPSC, 'Careful', CMPQ, IDUM, 4, M, 1, 1, M, IWORK,
     $             DWORK( IMAT ), M, M, DWORK, M, M, ALPHAR, ALPHAI,
     $             BETA, IWORK( 5 ), IWORK( M+5 ), LIWORK-( M+4 ),
     $             DWORK( IWRK ), LDWORK-IWRK+1, IWARN, INFO )
      IF( IWARN.GT.0 .AND. IWARN.LT.M ) THEN
         INFO = 1
         RETURN
      ELSE IF( IWARN.EQ.M+1 ) THEN
         INFO = 3
      ELSE IF( INFO.GT.0 ) THEN
         INFO = 2
         RETURN
      END IF
      OPTDW = MAX( OPTDW, INT( DWORK( IWRK ) ) + IWRK - 1 )
      I11   = 0
      I22   = 0
      I2R   = .TRUE.
C
C     Compute the eigenvalues with nonnegative imaginary parts of the
C     pencil aS - bH. Also, count the number of 1-by-1 or 2-by-2
C     diagonal blocks, corresponding to finite eigenvalues.
C
      DO 400 I = 1, M
         IF( IWORK( I+4 ).GE.2*EMIN .AND. IWORK( I+4 ).LE.2*EMAX ) THEN
C
C           B = SQRT(BASE**IWORK(i+4)) is between underflow and overflow
C           threshold, BETA(i) is divided by B.
C           Set to zero negligible real and imaginary parts.
C
            BETA( I ) = BETA( I )/BASE**( HALF*IWORK( I+4 ) )
            IF( BETA( I ).NE.ZERO ) THEN
               IF( INFO.EQ.0 .AND. IWORK( M+I+5 ).LT.0 ) THEN
                  IF( I2R ) THEN
                     I2R = .FALSE.
                     I22 = I22 + 2
                     IWORK( M+I+6 ) = -1
                     IWORK( I+4   ) =  2
                  ELSE
                     I2R = .TRUE.
                  END IF
               ELSE
                  IF( ALPHAI( I ).EQ.ZERO ) THEN
                     I11 = I11 + 1
                     IWORK( I+4 ) = 1
                  ELSE
                     I22 = I22 + 1
                     IWORK( I+4 ) = 2
                  END IF
               END IF
            ELSE
               IWORK( I+4 ) = 0
            END IF
            EIG = SQRT( CMPLX( ALPHAR( I ), ALPHAI( I ) ) )
            ALPHAR( I ) = AIMAG( EIG )
            ALPHAI( I ) = DBLE(  EIG )
            TMP2 = PREC*DLAPY2( ALPHAR( I ), ALPHAI( I ) )
            IF( ABS( ALPHAR( I ) ).LT.TMP2 ) THEN
               ALPHAR( I ) = ZERO
               IF( ALPHAI( I ).LT.ZERO )
     $             ALPHAI( I ) = -ALPHAI( I )
            END IF
            IF( ABS( ALPHAI( I ) ).LT.TMP2 ) THEN
               ALPHAI( I ) = ZERO
               IF( ALPHAR( I ).LT.ZERO )
     $             ALPHAR( I ) = -ALPHAR( I )
            END IF
         ELSE
C
C           Error.
C
            INFO = 1
            RETURN
         END IF
  400 CONTINUE
C
      I22 = I22/2
C
      DWORK( IWRK )   = DWORK( IWRK+1 )
      DWORK( IWRK+1 ) = DWORK( IWRK+2 )
      DWORK( IWRK+2 ) = DWORK( IWRK+3 )
      DWORK( IWRK+3 ) = DWORK( IWRK+4 )
      K = IWRK + 4
C
      IF( I11.GT.0 ) THEN
C
C        Save the I11 quadruples corresponding to the 1-by-1
C        diagonal blocks.
C
         DO 410  I = 1, M
            IF( IWORK( I+4 ).EQ.1 ) THEN
               DWORK( K )   = DWORK( IMAT+(M+1)*(I-1) )
               DWORK( K+1 ) = DWORK( IMAT+(M+1)*(I-1)+MM )
               DWORK( K+2 ) = DWORK( IMAT+(M+1)*(I-1)+2*MM )
               DWORK( K+3 ) = DWORK( IMAT+(M+1)*(I-1)+3*MM )
               K = K + 4
            END IF
  410    CONTINUE
C
      END IF
C
      IF( I22.GT.0 ) THEN
C
C        Save the I22 quadruple of 2-by-2 matrices corresponding to
C        the 2-by-2 diagonal blocks.
C
C        WHILE( I < M ) DO
         I = 1
  420    CONTINUE
            IF( IWORK( I+4 ).EQ.2 ) THEN
               CALL DLACPY( 'Full', 2, 2, DWORK( IMAT+(M+1)*(I-1) ), M,
     $                      DWORK( K ), 2 )
               CALL DLACPY( 'Full', 2, 2, DWORK( IMAT+(M+1)*(I-1)+MM ),
     $                      M, DWORK( K+4 ), 2 )
               CALL DLACPY( 'Full', 2, 2,
     $                      DWORK( IMAT+(M+1)*(I-1)+2*MM ), M,
     $                      DWORK( K+8 ), 2 )
               CALL DLACPY( 'Full', 2, 2,
     $                      DWORK( IMAT+(M+1)*(I-1)+3*MM ), M,
     $                      DWORK( K+12 ), 2 )
               K = K + 16
               I = I + 2
            ELSE
               I = I + 1
            END IF
            IF( I.LT.M )
     $         GO TO 420
C        END WHILE 420
C
      END IF
C
C     Check for possible loss of accuracy.
C
      IF( IWARN.GT.M ) THEN
         I = 2
C
C        Compress IWORK.
C
         DO 430  J = M + 6, 2*M + 5
            K = IWORK( J )
            IF( K.NE.0 ) THEN
               IWORK( I ) = K
               I = I + 1
            END IF
  430    CONTINUE
         IWORK( 1 ) = MAX( 0, I - 2 )
      ELSE
         IWORK( 1 ) = 0
      END IF
      I = IWORK( 1 )
      IWORK( I+2 ) = I11
      IWORK( I+3 ) = I22
C
      IF( LTRI ) THEN
C
C        Update C1 and C2.
C
         CALL DLACPY( 'Upper', M, M, DWORK( IMAT+2*MM ), M, C1, LDC1 )
         CALL DLACPY( 'Full',  M, M, DWORK( IMAT ), M, C2, LDC2 )
C
C        Update V.
C
         CALL DGEMM( 'Transpose', 'No Transpose', M, M, M, ONE,
     $               DWORK( 2*MM+1 ), M, VW( 1, 2 ), LDVW, ZERO,
     $               DWORK( IMAT ), M )
         CALL DGEMM( 'No Transpose', 'No Transpose', M, M, M, ONE,
     $               DWORK( IMAT ), M, DWORK, M, ZERO, VW( 1, 2 ),
     $               LDVW )
C
C        Update A.
C
         CALL DLACPY( 'Upper', M, M, DWORK( IMAT+MM ), M, A, LDA )
C
C        Skew-symmetric update of D.
C
         CALL MB01LD( 'Upper', 'Transpose', M, M, ZERO, ONE, DE( 1, 2 ),
     $                LDDE, DWORK( 2*MM+1 ), M, DE( 1, 2 ), LDDE,
     $                DWORK( IMAT ), MM, INFO )
C
C        Update B.
C
         CALL DLACPY( 'Upper', M, M, DWORK( IMAT+3*MM ), M, B, LDB )
C
C        Skew-symmetric update of F.
C
         CALL MB01LD( 'Upper', 'Transpose', M, M, ZERO, ONE, F, LDF,
     $                DWORK, M, F, LDF, DWORK( IMAT ), LDWORK-IMAT+1,
     $                INFO )
C
         IF( LCMPQ1 ) THEN
C
C           Update Q1.
C
            CALL DGEMM( 'No Transpose', 'No Transpose', N, M, M, ONE,
     $                  Q1, LDQ1, DWORK( 2*MM+1 ), M, ZERO,
     $                  DWORK( IMAT ), N )
            CALL DLACPY( 'Full', N, M, DWORK( IMAT ), N, Q1, LDQ1 )
            CALL DGEMM( 'No Transpose', 'No Transpose', N, M, M, ONE,
     $                  Q1( 1, M+1 ), LDQ1, DWORK( MM+1 ), M, ZERO,
     $                  DWORK( IMAT ), N )
            CALL DLACPY( 'Full', N, M, DWORK( IMAT ), N, Q1( 1, M+1 ),
     $                   LDQ1 )
         END IF
C
         IF( LCMPQ2 ) THEN
C
C           Update Q2.
C
            CALL DGEMM( 'No Transpose', 'No Transpose', N, M, M, ONE,
     $                  Q2, LDQ2, DWORK( 3*MM+1 ), M, ZERO,
     $                  DWORK( IMAT ), N )
            CALL DLACPY( 'Full', N, M, DWORK( IMAT ), N, Q2, LDQ2 )
            CALL DGEMM( 'No Transpose', 'No Transpose', N, M, M, ONE,
     $                  Q2( 1, M+1 ), LDQ2, DWORK, M, ZERO,
     $                  DWORK( IMAT ), N )
            CALL DLACPY( 'Full', N, M, DWORK( IMAT ), N, Q2( 1, M+1 ),
     $                   LDQ2 )
         END IF
      END IF
C
C     Move the norms, the I11 quadruples, and the I22 quadruple
C     2-by-2 matrices in front.
C
      K = IWRK + 4*I11 + 16*I22 + 3
C
      DO 440  J = 4*I11 + 16*I22 + 5, 2, -1
         DWORK( J ) = DWORK( K )
         K = K - 1
  440 CONTINUE
C
      DWORK( 1 ) = OPTDW
      RETURN
C *** Last line of DGHUTP ***
      END
