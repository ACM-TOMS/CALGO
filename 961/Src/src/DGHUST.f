      SUBROUTINE DGHUST( JOB, COMPQ, N, A, LDA, DE, LDDE, B, LDB,
     $                   FG, LDFG, Q, LDQ, ALPHAR, ALPHAI, BETA, DWORK,
     $                   LDWORK, INFO )
C
C     PURPOSE
C
C     To compute the eigenvalues of a real N-by-N skew-Hamiltonian/
C     skew-Hamiltonian pencil aS - bT with
C
C           (  A  D  )         (  B  F  )
C       S = (        ) and T = (        ).                           (1)
C           (  E  A' )         (  G  B' )
C
C     Optionally, if JOB = 'T', the pencil aS - bT will be transformed
C     to the structured Schur form: an orthogonal transformation matrix
C     Q is computed such that
C
C                     (  Aout  Dout  )
C       J Q' J' S Q = (              ), and
C                     (   0    Aout' )
C                                                                    (2)
C                     (  Bout  Fout  )            (  0  I  )
C       J Q' J' T Q = (              ), where J = (        ),
C                     (   0    Bout' )            ( -I  0  )
C
C     Aout is upper triangular, and Bout is upper quasi-triangular. The
C     notation M' denotes the transpose of the matrix M.
C     Optionally, if COMPQ = 'I' or COMPQ = 'U', the orthogonal
C     transformation matrix Q will be computed.
C
C     ARGUMENTS
C
C     Mode Parameters
C
C     JOB     CHARACTER*1
C             Specifies the computation to be performed, as follows:
C             = 'E':  compute the eigenvalues only; S and T will not
C                     necessarily be put into skew-Hamiltonian
C                     triangular form (2);
C             = 'T':  put S and T into skew-Hamiltonian triangular form
C                     (2), and return the eigenvalues in ALPHAR, ALPHAI
C                     and BETA.
C
C     COMPQ   CHARACTER*1
C             Specifies whether to compute the orthogonal transformation
C             matrix Q as follows:
C             = 'N':  Q is not computed;
C             = 'I':  the array Q is initialized internally to the unit
C                     matrix, and the orthogonal matrix Q is returned;
C             = 'U':  the array Q contains an orthogonal matrix Q0 on
C                     entry, and the product Q0*Q is returned, where Q
C                     is the product of the orthogonal transformations
C                     that are applied to the pencil aS - bT to reduce
C                     S and T to the forms in (2), for COMPQ = 'I'.
C
C     Input/Output Parameters
C
C     N       (input) INTEGER
C             The order of the pencil aS - bT.  N >= 0, even.
C
C     A       (input/output) DOUBLE PRECISION array, dimension
C                            (LDA, N/2)
C             On entry, the leading N/2-by-N/2 part of this array must
C             contain the matrix A.
C             On exit, if JOB = 'T', the leading N/2-by-N/2 part of this
C             array contains the matrix Aout; otherwise, it contains
C             meaningless elements.
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
C             this array are not referenced, but are assumed to be zero.
C             On exit, if JOB = 'T', the leading N/2-by-N/2 strictly
C             upper triangular part of the submatrix in the columns
C             2 to N/2+1 of this array contains the strictly upper
C             triangular part of the skew-symmetric matrix Dout.
C             If JOB = 'E', the leading N/2-by-N/2 strictly upper
C             triangular part of the submatrix in the columns 2 to N/2+1
C             of this array contains the strictly upper triangular part
C             of the skew-symmetric matrix D just before the application
C             of the QZ algorithm. The remaining entries are
C             meaningless.
C
C     LDDE    INTEGER
C             The leading dimension of the array DE.
C             LDDE >= MAX(1, N/2).
C
C     B       (input/output) DOUBLE PRECISION array, dimension
C                            (LDB, N/2)
C             On entry, the leading N/2-by-N/2 part of this array must
C             contain the matrix B.
C             On exit, if JOB = 'T', the leading N/2-by-N/2 part of this
C             array contains the matrix Bout; otherwise, it contains
C             meaningless elements.
C
C     LDB     INTEGER
C             The leading dimension of the array B.  LDB >= MAX(1, N/2).
C
C     FG      (input/output) DOUBLE PRECISION array, dimension
C                            (LDFG, N/2+1)
C             On entry, the leading N/2-by-N/2 strictly lower triangular
C             part of this array must contain the strictly lower
C             triangular part of the skew-symmetric matrix G, and the
C             N/2-by-N/2 strictly upper triangular part of the submatrix
C             in the columns 2 to N/2+1 of this array must contain the
C             strictly upper triangular part of the skew-symmetric
C             matrix F.
C             The entries on the diagonal and the first superdiagonal of
C             this array are not referenced, but are assumed to be zero.
C             On exit, if JOB = 'T', the leading N/2-by-N/2 strictly
C             upper triangular part of the submatrix in the columns
C             2 to N/2+1 of this array contains the strictly upper
C             triangular part of the skew-symmetric matrix Fout.
C             If JOB = 'E', the leading N/2-by-N/2 strictly upper
C             triangular part of the submatrix in the columns 2 to N/2+1
C             of this array contains the strictly upper triangular part
C             of the skew-symmetric matrix F just before the application
C             of the QZ algorithm. The remaining entries are
C             meaningless.
C
C     LDFG    INTEGER
C             The leading dimension of the array FG.
C             LDFG >= MAX(1, N/2).
C
C     Q       (input/output) DOUBLE PRECISION array, dimension (LDQ, N)
C             On entry, if COMPQ = 'U', then the leading N-by-N part of
C             this array must contain a given matrix Q0, and on exit,
C             the leading N-by-N part of this array contains the product
C             of the input matrix Q0 and the transformation matrix Q
C             used to transform the matrices S and T.
C             On exit, if COMPQ = 'I', then the leading N-by-N part of
C             this array contains the orthogonal transformation matrix
C             Q.
C             If COMPQ = 'N', this array is not referenced.
C
C     LDQ     INTEGER
C             The leading dimension of the array Q.
C             LDQ >= 1,         if COMPQ = 'N';
C             LDQ >= MAX(1, N), if COMPQ = 'I' or COMPQ = 'U'.
C
C     ALPHAR  (output) DOUBLE PRECISION array, dimension (N/2)
C             The real parts of each scalar alpha defining an eigenvalue
C             of the pencil aS - bT.
C
C     ALPHAI  (output) DOUBLE PRECISION array, dimension (N/2)
C             The imaginary parts of each scalar alpha defining an
C             eigenvalue of the pencil aS - bT.
C             If ALPHAI(j) is zero, then the j-th eigenvalue is real; if
C             positive, then the j-th and (j+1)-st eigenvalues are a
C             complex conjugate pair.
C
C     BETA    (output) DOUBLE PRECISION array, dimension (N/2)
C             The scalars beta that define the eigenvalues of the pencil
C             aS - bT.
C             Together, the quantities alpha = (ALPHAR(j),ALPHAI(j)) and
C             beta = BETA(j) represent the j-th eigenvalue of the pencil
C             aS - bT, in the form lambda = alpha/beta. Since lambda may
C             overflow, the ratios should not, in general, be computed.
C             Due to the skew-Hamiltonian/skew-Hamiltonian structure of
C             the pencil, every eigenvalue occurs twice and thus it has
C             only to be saved once in ALPHAR, ALPHAI and BETA.
C
C     Workspace
C
C     DWORK   DOUBLE PRECISION array, dimension (LDWORK)
C             On exit, if INFO = 0, DWORK(1) returns the optimal LDWORK.
C             On exit, if INFO = -18, DWORK(1) returns the minimum
C             value of LDWORK.
C
C     LDWORK  INTEGER
C             The dimension of the array DWORK.
C             LDWORK >= MAX(1,N/2),        if JOB = 'E' and COMPQ = 'N';
C             LDWORK >= MAX(1,N**2/4+N/2), if JOB = 'T' and COMPQ = 'N';
C             LDWORK >= MAX(1,3*N**2/4),   if               COMPQ<> 'N'.
C             For good performance LDWORK should generally be larger.
C
C             If LDWORK = -1, then a workspace query is assumed;
C             the routine only calculates the optimal size of the
C             DWORK array, returns this value as the first entry of
C             the DWORK array, and no error message related to LDWORK
C             is issued by XERBLA.
C
C     Error Indicator
C
C     INFO    INTEGER
C             = 0: succesful exit;
C             < 0: if INFO = -i, the i-th argument had an illegal value;
C             = 1: QZ iteration failed in the LAPACK Library routine
C                  DHGEQZ. (QZ iteration did not converge or computation
C                  of the shifts failed.)
C
C     METHOD
C
C     The algorithm uses Givens rotations and Householder reflections to
C     annihilate elements in S and T such that S is in skew-Hamiltonian
C     triangular form and T is in skew-Hamiltonian Hessenberg form:
C
C         (  A1  D1  )      (  B1  F1  )
C     S = (          ), T = (          ),
C         (   0  A1' )      (   0  B1' )
C
C     where A1 is upper triangular and B1 is upper Hessenberg.
C     Subsequently, the QZ algorithm is applied to the pencil aA1 - bB1
C     to determine orthogonal matrices Q1 and Q2 such that
C     Q2' A1 Q1 is upper triangular and Q2' B1 Q1 is upper quasi-
C     triangular.
C     See also page 40 in [1] for more details.
C
C     REFERENCES
C
C     [1] Benner, P., Byers, R., Mehrmann, V. and Xu, H.
C         Numerical Computation of Deflating Subspaces of Embedded
C         Hamiltonian Pencils.
C         Tech. Rep. SFB393/99-15, Technical University Chemnitz,
C         Germany, June 1999.
C
C     NUMERICAL ASPECTS
C                                                               3
C     The algorithm is numerically backward stable and needs O(N )
C     real floating point operations.
C
C     CONTRIBUTOR
C
C     M. Voigt, Technische Universitaet Chemnitz, Mar. 2009.
C     V. Sima, Research Institute for Informatics, Bucharest, Nov. 2009.
C
C     REVISIONS
C
C     V. Sima, Dec. 2010, Jan. 2011, Aug. 2011.
C     M. Voigt, Jan. 2012, Jul. 2013. 
C
C     KEYWORDS
C
C     QZ algorithm, upper (quasi-)triangular matrix,
C     skew-Hamiltonian/skew-Hamiltonian pencil.
C
C     ******************************************************************
C
C     .. Parameters ..
      DOUBLE PRECISION   ZERO, HALF, ONE, TWO
      PARAMETER          ( ZERO = 0.0D+0, HALF = 0.5D+0, ONE = 1.0D+0,
     $                     TWO  = 2.0D+0 )
C
C     .. Scalar Arguments ..
      CHARACTER          COMPQ, JOB
      INTEGER            INFO, LDA, LDB, LDDE, LDFG, LDQ, LDWORK, N
C
C     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), ALPHAI( * ), ALPHAR( * ),
     $                   B( LDB, * ), BETA( * ), DE( LDDE, * ),
     $                   DWORK( * ), FG( LDFG, * ), Q( LDQ, * )
C
C     .. Local Scalars ..
      LOGICAL            LCMPQ, LINIQ, LQUERY, LTRI, LUPDQ
      CHARACTER*16       CMPQ, CMPSC, CMPZ
      INTEGER            IQ1, IQ2, IWRK, J, K, M, M1, MINDW, MJ1, MJ2,
     $                   MJ3, MK2, MK3, MM, OPTDW
      DOUBLE PRECISION   CO, MU, NU, SI, TMP1, TMP2
C
C     .. Local Arrays ..
      DOUBLE PRECISION   DUM( 1 )
C
C     .. External Functions ..
      LOGICAL            LSAME
      DOUBLE PRECISION   DDOT
      EXTERNAL           DDOT, LSAME
C
C     .. External Subroutines ..
      EXTERNAL           DAXPY, DCOPY, DGEMM, DHGEQZ, DLACPY, DLARF,
     $                   DLARFG, DLARTG, DLASET, DROT, MB01LD, MB01MD,
     $                   MB01ND, XERBLA
C
C     .. Intrinsic Functions ..
      INTRINSIC          INT, MAX, MIN, MOD
C
C     .. Executable Statements ..
C
C     Decode the input arguments.
C
      M  = N/2
      MM = M*M
      M1 = MAX( 1, M )
      LTRI   = LSAME( JOB,   'T' )
      LINIQ  = LSAME( COMPQ, 'I' )
      LUPDQ  = LSAME( COMPQ, 'U' )
      LCMPQ  = LINIQ .OR. LUPDQ
      LQUERY = LDWORK.EQ.-1
C
C     Determine the mode of computations.
C
      IQ1 = 1
      IF( LCMPQ ) THEN
         CMPQ  = 'Initialize'
         CMPZ  = CMPQ
         IQ2   = IQ1 + MM
         IWRK  = IQ2 + MM
         MINDW = MAX( 1, IWRK - 1 + MM )
      ELSE IF( LTRI ) THEN
         CMPQ  = 'Initialize'
         CMPZ  = 'No Computation'
         IQ2   = 1
         IWRK  = IQ2 + MM
         MINDW = MAX( 1, IWRK - 1 + M )
      ELSE
         CMPQ  = 'No Computation'
         CMPZ  = CMPQ
         IQ2   = 1
         IWRK  = 1
         MINDW = M1
      END IF
      K = IWRK - 1
C
      IF( LTRI ) THEN
         CMPSC = 'Schur Form'
      ELSE
         CMPSC = 'Eigenvalues Only'
      END IF
C
C     Test the input arguments.
C
      INFO = 0
      IF( .NOT.( LSAME( JOB, 'E' ) .OR. LTRI ) ) THEN
         INFO = -1
      ELSE IF( .NOT.( LSAME( COMPQ, 'N' ) .OR. LCMPQ ) ) THEN
         INFO = -2
      ELSE IF( N.LT.0 .OR. MOD( N, 2 ).NE.0 ) THEN
         INFO = -3
      ELSE IF(  LDA.LT.M1 ) THEN
         INFO = -5
      ELSE IF( LDDE.LT.M1 ) THEN
         INFO = -7
      ELSE IF(  LDB.LT.M1 ) THEN
         INFO = -9
      ELSE IF( LDFG.LT.M1 ) THEN
         INFO = -11
      ELSE IF(  LDQ.LT.1 .OR. ( LCMPQ .AND. LDQ.LT.N ) ) THEN
         INFO = -13
      ELSE IF( LDWORK.LT.MINDW .AND. .NOT.LQUERY ) THEN
         DWORK( 1 ) = MINDW
         INFO = -18
      END IF
C
      IF( INFO.NE.0) THEN
         CALL XERBLA( 'DGHUST', -INFO )
         RETURN
      ELSE IF( N.GT.0 .AND. LQUERY ) THEN
         CALL DHGEQZ( CMPSC, CMPQ, CMPZ, M, 1, M, B, LDB, A, LDA,
     $                ALPHAR, ALPHAI, BETA, DWORK, M1, DWORK, M1,
     $                DWORK, -1, INFO )
C
         IF( LCMPQ ) THEN
            OPTDW = K + MAX( K, INT( DWORK( 1 ) ) )
         ELSE IF( LTRI ) THEN
            OPTDW = K + MAX( K - M, INT( DWORK( 1 ) ) )
         ELSE
            OPTDW = INT( DWORK( 1 ) )
         END IF
         DWORK( 1 ) = OPTDW
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
C     STEP 1: Reduce S to skew-Hamiltonian triangular form.
C
C     Workspace:    need   N, if COMPQ <> 'N';
C                          M, if COMPQ  = 'N'.
C
      IF( LINIQ )
     $   CALL DLASET( 'Full', N, N, ZERO, ONE, Q, LDQ )
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
C           Apply H(k) to G(k+1:m,1:k) from the left (and implicitly to
C           G(1:k,k+1:m) from the right).
C
            CALL DLARF( 'Left', M-K, K, DE( K+1, K ), 1, NU,
     $                  FG( K+1, 1 ), LDFG, DWORK )
C
C           Apply H(k) from both sides to G(k+1:m,k+1:m).
C           Compute  x := nu * G(k+1:m,k+1:m) * v.
C
            CALL MB01MD( 'Lower', M-K, NU, FG( K+1, K+1 ), LDFG,
     $                   DE( K+1, K ), 1, ZERO, DWORK, 1 )
C
C           Compute  w := x - 1/2 * nu * (x'*v) * v.
C
            MU = -HALF*NU*DDOT( M-K, DWORK, 1, DE( K+1, K ), 1 )
            CALL DAXPY( M-K, MU, DE( K+1, K ), 1, DWORK, 1 )
C
C           Apply the transformation as a skew-symmetric rank-2 update:
C                G := G + v * w' - w * v'.
C
            CALL MB01ND( 'Lower', M-K, ONE, DE( K+1, K ), 1, DWORK, 1,
     $                   FG( K+1, K+1 ), LDFG )
C
C           Apply H(k) from the right hand side to A(1:m,k+1:m) and
C           B(1:m,k+1:m).
C
            CALL DLARF( 'Right', M, M-K, DE( K+1, K ), 1, NU,
     $                  A( 1, K+1 ), LDA, DWORK )
            CALL DLARF( 'Right', M, M-K, DE( K+1, K ), 1, NU,
     $                  B( 1, K+1 ), LDB, DWORK )
C
            IF( LCMPQ ) THEN
C
C              Apply H(k) from the right hand side to Q(1:n,k+1:m).
C
               CALL DLARF( 'Right', N, M-K, DE( K+1, K ), 1, NU,
     $                     Q( 1, K+1 ), LDQ, DWORK )
            END IF
            DE( K+1, K ) = TMP1
         END IF
C
C        Determine a Givens rotation to annihilate E(k+1,k) from the
C        left.
C
         TMP2 = A( K+1, K )
         CALL DLARTG( TMP2, TMP1, CO, SI, A( K+1, K ) )
C
C        Update A, D and E.
C
         CALL DROT( M-K-1, DE( MK2, K+1 ), 1, A( K+1, MK2 ), LDA, CO,
     $              SI )
         CALL DROT( K, A( 1, K+1 ), 1, DE( 1, K+2 ), 1, CO, SI )
         CALL DROT( M-K-1, DE( K+1, MK3 ), LDDE, A( MK2, K+1 ), 1, CO,
     $              SI )
C
C        Update B, G and F.
C
         CALL DROT( K, FG( K+1, 1 ), LDFG, B( K+1, 1 ), LDB, CO, -SI )
         CALL DROT( M-K-1, FG( MK2, K+1 ), 1, B( K+1, MK2 ), LDB, CO,
     $              SI )
         CALL DROT( K, B( 1, K+1 ), 1, FG( 1, K+2 ), 1, CO, SI )
         CALL DROT( M-K-1, FG( K+1, MK3 ), LDFG, B( MK2, K+1 ), 1, CO,
     $              SI )
C
         IF( LCMPQ ) THEN
C
C           Update Q.
C
            CALL DROT( N, Q( 1, M+K+1 ), 1, Q( 1, K+1 ), 1, CO, -SI )
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
C           Apply P(k) to D(1:k-1,k:m) from the right (and implicitly to
C           D(k:m,1:k-1) from the left).
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
C           Apply P(k) from the left hand side to B(k:m,1:m).
C
            CALL DLARF( 'Left', M-K+1, M, A( K, K ), 1, NU, B( K, 1 ),
     $                  LDB, DWORK )
C
C           Apply P(k) to F(1:k-1,k:m) from the right (and implicitly to
C           F(k:m,1:k-1) from the left).
C
            CALL DLARF( 'Right', K-1, M-K+1, A( K, K ), 1, NU,
     $                  FG( 1, K+1 ), LDFG, DWORK )
C
C           Apply P(k) from both sides to F(k:m,k:m).
C           Compute  x := nu * F(k:m,k:m) * v.
C
            CALL MB01MD( 'Upper', M-K+1, NU, FG( K, K+1 ), LDFG,
     $                   A( K, K ), 1, ZERO, DWORK, 1 )
C
C           Compute  w := x - 1/2 * nu * (x'*v) * v.
C
            MU = -HALF*NU*DDOT( M-K+1, DWORK, 1, A( K, K ), 1 )
            CALL DAXPY( M-K+1, MU, A( K, K ), 1, DWORK, 1 )
C
C           Apply the transformation as a skew-symmetric rank-2 update:
C                F := F + v * w' - w * v'.
C
            CALL MB01ND( 'Upper', M-K+1, ONE, A( K, K ), 1, DWORK, 1,
     $                   FG( K, K+1 ), LDFG )
C
            IF( LCMPQ ) THEN
C
C              Apply P(k) from the right hand side to Q(1:n,m+k:n).
C
               CALL DLARF( 'Right', N, M-K+1, A( K, K ), 1, NU,
     $                     Q( 1, M+K ), LDQ, DWORK )
            END IF
            A( K, K ) = TMP1
         END IF
C
C        Set A(k+1:m,k) to zero in order to be able to apply DHGEQZ.
C
         CALL DCOPY( M-K, DUM, 0, A( K+1, K ), 1 )
   10 CONTINUE
C
C     STEP 2: Reduce T to skew-Hamiltonian triangular form.
C
      DO 40 K = 1, M - 1
C
C        I. Annihilate G(k+1:m-1,k).
C
         DO 20 J = K + 1, M - 1
            MJ2 = MIN( J+2, M )
            MJ3 = MJ2 + 1
C
C           Determine a Givens rotation to annihilate G(j,k) from the
C           left.
C
            CALL DLARTG( FG( J+1, K ), FG( J, K ), CO, SI, TMP1 )
C
C           Update B and G.
C
            CALL DROT( M, B( 1, J+1 ), 1, B( 1, J ), 1, CO, SI )
            FG( J+1, K ) = TMP1
            CALL DROT( M-J-1, FG( MJ2, J+1 ), 1, FG( MJ2, J ), 1, CO, SI
     $                  )
            CALL DROT( J-K-1, FG( J+1, K+1 ), LDFG, FG( J, K+1 ), LDFG,
     $                 CO, SI )
C
C           Update A.
C
            CALL DROT( J, A( 1, J+1 ), 1, A( 1, J ), 1, CO, SI )
            TMP1          = -SI*A( J+1, J+1 )
            A( J+1, J+1 ) =  CO*A( J+1, J+1 )
C
            IF( LCMPQ ) THEN
C
C              Update Q.
C
               CALL DROT( N, Q( 1, J+1 ), 1, Q( 1, J ), 1, CO, SI )
            END IF
C
C           Determine a Givens rotation to annihilate A(j+1,j) from the
C           left.
C
            CALL DLARTG( A( J, J ), TMP1, CO, SI, TMP2 )
C
C           Update A and D.
C
            A( J, J ) = TMP2
            CALL DROT( M-J, A( J, J+1 ), LDA, A( J+1, J+1 ), LDA, CO, SI
     $                  )
            CALL DROT( J-1, DE( 1, J+1 ), 1, DE( 1, J+2 ), 1, CO, SI )
            CALL DROT( M-J-1, DE( J, MJ3 ), LDDE, DE( J+1, MJ3 ), LDDE,
     $                 CO, SI )
C
C           Update B and F.
C
            CALL DROT( M-K+1, B( J, K ), LDB, B( J+1, K ), LDB, CO, SI )
            CALL DROT( J-1, FG( 1, J+1 ), 1, FG( 1, J+2 ), 1, CO, SI )
            CALL DROT( M-J-1, FG( J, MJ3 ), LDFG, FG( J+1, MJ3 ), LDFG,
     $                 CO, SI )
C
            IF( LCMPQ ) THEN
C
C              Update Q.
C
               CALL DROT( N, Q( 1, M+J ), 1, Q( 1, M+J+1 ), 1, CO, SI )
            END IF
   20    CONTINUE
C
C        II. Annihilate G(m,k).
C
C        Determine a Givens rotation to annihilate G(m,k) from the
C        left.
C
         CALL DLARTG( B( M, K ), -FG( M, K ), CO, SI, TMP1 )
C
C        Update B and G.
C
         B(  M, K ) = TMP1
         FG( M, K ) = ZERO
         CALL DROT( M-1, FG( 1, M+1 ), 1, B( 1, M ), 1, CO, SI )
         CALL DROT( M-K-1, FG( M, K+1 ), LDFG, B( M, K+1 ), LDB, CO, SI
     $               )
C
C        Update A and D.
C
         CALL DROT( M-1, DE( 1, M+1 ), 1, A( 1, M ), 1, CO, SI )
C
         IF( LCMPQ ) THEN
C
C           Update Q.
C
            CALL DROT( N, Q( 1, N ), 1, Q( 1, M ), 1, CO, SI )
         END IF
C
C        III. Annihilate B(k+2:m,k).
C
         DO 30 J = M, K + 2, -1
            MJ1 = MIN( J+1, M )
            MJ2 = MJ1 + 1
C
C           Determine a Givens rotation to annihilate B(j,k) from the
C           left.
C
            CALL DLARTG( B( J-1, K ), B( J, K ), CO, SI, TMP1 )
C
C           Update B and F.
C
            B( J-1, K ) = TMP1
            B(   J, K ) = ZERO
            CALL DROT( J-2, FG( 1, J ), 1, FG( 1, J+1 ), 1, CO, SI )
            CALL DROT( M-K, B( J-1, K+1 ), LDB, B( J, K+1 ), LDB, CO, SI
     $                  )
            CALL DROT( M-J, FG( J-1, MJ2 ), LDFG, FG( J, MJ2 ), LDFG,
     $                 CO, SI )
C
C           Update A and D.
C
            TMP1          = -SI*A( J-1, J-1 )
            A( J-1, J-1 ) =  CO*A( J-1, J-1 )
            CALL DROT( M-J+1, A( J-1, J ), LDA, A( J, J ), LDA, CO, SI )
            CALL DROT( J-2, DE( 1, J ), 1, DE( 1, J+1 ), 1, CO, SI )
            CALL DROT( M-J, DE( J-1, MJ2 ), LDDE, DE( J, MJ2 ), LDDE,
     $                 CO, SI )
C
            IF( LCMPQ ) THEN
C
C              Update Q.
C
               CALL DROT( N, Q( 1, M+J-1 ), 1, Q( 1, M+J ), 1, CO, SI )
            END IF
C
C           Determine a Givens rotation to annihilate A(j,j-1) from the
C           right.
C
            CALL DLARTG( A( J, J ), TMP1, CO, SI, TMP2 )
C
C           Update A.
C
            A( J, J ) = TMP2
            CALL DROT( J-1, A( 1, J ), 1, A( 1, J-1 ), 1, CO, SI )
C
C           Update B and G.
C
            CALL DROT( M, B( 1, J ), 1, B( 1, J-1 ), 1, CO, SI )
            CALL DROT( J-K-1, FG( J, K ), LDFG, FG( J-1, K ), LDFG, CO,
     $                 SI )
            CALL DROT( M-J, FG( MJ1, J ), 1, FG( MJ1, J-1 ), 1, CO, SI )
C
            IF( LCMPQ ) THEN
C
C              Update Q.
C
               CALL DROT( N, Q( 1, J ), 1, Q( 1, J-1 ), 1, CO, SI )
            END IF
   30    CONTINUE
   40 CONTINUE
C
C                     (  A1  D1  )         (  B1  F1  )
C     Now we have S = (          ) and T = (          ),
C                     (   0  A1' )         (   0  B1' )
C
C     where A1 is upper triangular and B1 is upper Hessenberg.
C
C     STEP 3: Apply the QZ algorithm to the pencil aA1 - bB1 to
C             determine orthogonal matrices Q1 and Q2 such that
C             Q2' A1 Q1 is upper triangular and Q2' B1 Q1 is upper
C             quasi-triangular.
C
C     Workspace:    need   w + M, where
C                          w = 2*M**2, if COMPQ <> 'N';
C                          w =   M**2, if COMPQ  = 'N' and JOB = 'T';
C                          w =      0, if COMPQ  = 'N' and JOB = 'E';
C                   prefer larger.
C
      CALL DHGEQZ( CMPSC, CMPQ, CMPZ, M, 1, M, B, LDB, A, LDA, ALPHAR,
     $             ALPHAI, BETA, DWORK( IQ1 ), M, DWORK( IQ2 ), M,
     $             DWORK( IWRK ), LDWORK-IWRK+1, INFO )
      IF( INFO.GT.0 ) THEN
         INFO = 1
         RETURN
      END IF
      OPTDW = MAX( MINDW, INT( DWORK( IWRK ) ) + IWRK - 1 )
C
C     Enforce the needed equalities in complex eigenvalues.
C
      J = 1
C     WHILE( J.LT.M ) DO
   50 CONTINUE
      IF( J.LT.M ) THEN
         IF( ALPHAI( J ).NE.ZERO ) THEN
            IF( BETA( J ).GE. BETA( J + 1 ) ) THEN
               TMP2 = BETA( J + 1 )/BETA( J )
               TMP1 = ( ALPHAR( J )*TMP2 + ALPHAR( J+1 ) )/TWO
               TMP2 = ( ALPHAI( J )*TMP2 - ALPHAI( J+1 ) )/TWO
               BETA( J ) = BETA( J+1 )
            ELSE
               TMP2 = BETA( J )/BETA( J + 1 )
               TMP1 = ( ALPHAR( J+1 )*TMP2 + ALPHAR( J ) )/TWO
               TMP2 = ( ALPHAI( J+1 )*TMP2 - ALPHAI( J ) )/TWO
               BETA( J+1 ) = BETA( J )
            END IF
            ALPHAR( J )   =  TMP1
            ALPHAR( J+1 ) =  TMP1
            ALPHAI( J )   =  TMP2
            ALPHAI( J+1 ) = -TMP2
            J = J + 2
         ELSE
            J = J + 1
         END IF
         GO TO 50
      END IF
C     END WHILE 50
C
      IF( LTRI ) THEN
C
C        Skew-symmetric update of D.
C
C        Workspace:    need   w + M;
C                      prefer w + M*(M-1).
C
         CALL MB01LD( 'Upper', 'Transpose', M, M, ZERO, ONE, DE( 1, 2 ),
     $                LDDE, DWORK( IQ1 ), M, DE( 1, 2 ), LDDE,
     $                DWORK( IWRK ), LDWORK-IWRK+1, INFO )
C
C        Skew-symmetric update of F.
C
         CALL MB01LD( 'Upper', 'Transpose', M, M, ZERO, ONE, FG( 1, 2 ),
     $                LDFG, DWORK( IQ1 ), M, FG( 1, 2 ), LDFG,
     $                DWORK( IWRK ), LDWORK-IWRK+1, INFO )
      END IF
C
      IF( LCMPQ ) THEN
C
C        Update Q.
C        Workspace:    need   3*M*M;
C                      prefer 4*M*M.
C
         IF( LDWORK.GE.N*N ) THEN
            CALL DGEMM(  'No Transpose', 'No Transpose', N, M, M, ONE,
     $                   Q, LDQ, DWORK( IQ2 ), M, ZERO, DWORK( IWRK ),
     $                   N )
            CALL DLACPY( 'Full', N, M, DWORK( IWRK ), N, Q, LDQ )
            CALL DGEMM(  'No Transpose', 'No Transpose', N, M, M, ONE,
     $                   Q( 1, M+1 ), LDQ, DWORK( IQ1 ), M, ZERO,
     $                   DWORK( IWRK ), N )
            CALL DLACPY( 'Full', N, M, DWORK( IWRK ), N, Q( 1, M+1 ),
     $                   LDQ )
         ELSE
            CALL DGEMM(  'No Transpose', 'No Transpose', M, M, M, ONE,
     $                   Q, LDQ, DWORK( IQ2 ), M, ZERO, DWORK( IWRK ),
     $                   M )
            CALL DLACPY( 'Full', M, M, DWORK( IWRK ), M, Q, LDQ )
            CALL DGEMM(  'No Transpose', 'No Transpose', M, M, M, ONE,
     $                   Q( M+1, 1 ), LDQ, DWORK( IQ2 ), M, ZERO,
     $                   DWORK( IWRK ), M )
            CALL DLACPY( 'Full', M, M, DWORK( IWRK ), M, Q( M+1, 1 ),
     $                   LDQ )
            CALL DGEMM(  'No Transpose', 'No Transpose', M, M, M, ONE,
     $                   Q( 1, M+1 ), LDQ, DWORK( IQ1 ), M, ZERO,
     $                   DWORK( IWRK ), M )
            CALL DLACPY( 'Full', M, M, DWORK( IWRK ), M, Q( 1, M+1 ),
     $                   LDQ )
            CALL DGEMM(  'No Transpose', 'No Transpose', M, M, M, ONE,
     $                   Q( M+1, M+1 ), LDQ, DWORK( IQ1 ), M, ZERO,
     $                   DWORK( IWRK ), M )
            CALL DLACPY( 'Full', M, M, DWORK( IWRK ), M, Q( M+1, M+1 ),
     $                   LDQ )
         END IF
      END IF
C
      DWORK( 1 ) = OPTDW
      RETURN
C *** Last line of DGHUST ***
      END
