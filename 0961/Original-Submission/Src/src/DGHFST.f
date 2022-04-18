      SUBROUTINE DGHFST( JOB, COMPQ, COMPU, N, Z, LDZ, B, LDB, FG, LDFG,
     $                   Q, LDQ, U1, LDU1, U2, LDU2, ALPHAR, ALPHAI,
     $                   BETA, IWORK, LIWORK, DWORK, LDWORK, INFO )
C
C     PURPOSE
C
C     To compute the eigenvalues of a real N-by-N skew-Hamiltonian/
C     skew-Hamiltonian pencil aS - bT with
C
C                             (  B  F  )            (  0  I  )
C       S = J Z' J' Z and T = (        ), where J = (        ).      (1)
C                             (  G  B' )            ( -I  0  )
C
C     Optionally, if JOB = 'T', the pencil aS - bT will be transformed
C     to the structured Schur form: an orthogonal transformation matrix
C     Q and an orthogonal symplectic transformation matrix U are
C     computed, such that
C
C                (  Z11  Z12  )
C       U' Z Q = (            ) = Zout, and
C                (   0   Z22  )
C                                                                    (2)
C                     (  Bout  Fout  )
C       J Q' J' T Q = (              ),
C                     (   0    Bout' )
C
C     where Z11 and Z22' are upper triangular and Bout is upper quasi-
C     triangular. The notation M' denotes the transpose of the matrix M.
C     Optionally, if COMPQ = 'I', the orthogonal transformation matrix Q
C     will be computed.
C     Optionally, if COMPU = 'I' or COMPU = 'U', the orthogonal
C     symplectic transformation matrix
C
C           (  U1  U2  )
C       U = (          )
C           ( -U2  U1  )
C
C     will be computed.
C
C     ARGUMENTS
C
C     Mode Parameters
C
C     JOB     CHARACTER*1
C             Specifies the computation to be performed, as follows:
C             = 'E':  compute the eigenvalues only; Z and T will not
C                     necessarily be put into the forms in (2);
C             = 'T':  put Z and T into the forms in (2), and return the
C                     eigenvalues in ALPHAR, ALPHAI and BETA.
C
C     COMPQ   CHARACTER*1
C             Specifies whether to compute the orthogonal transformation
C             matrix Q as follows:
C             = 'N':  Q is not computed;
C             = 'I':  the array Q is initialized internally to the unit
C                     matrix, and the orthogonal matrix Q is returned.
C
C     COMPU   CHARACTER*1
C             Specifies whether to compute the orthogonal symplectic
C             transformation matrix U as follows:
C             = 'N':  U is not computed;
C             = 'I':  the array U is initialized internally to the unit
C                     matrix, and the orthogonal matrix U is returned;
C             = 'U':  the arrays U1 and U2 contain the corresponding
C                     submatrices of an orthogonal symplectic matrix U0
C                     on entry, and the updated submatrices U1 and U2
C                     of the matrix product U0*U are returned, where U
C                     is the product of the orthogonal symplectic
C                     transformations that are applied to the pencil
C                     aS - bT to reduce Z and T to the forms in (2), for
C                     COMPU = 'I'.
C
C     Input/Output Parameters
C
C     N       (input) INTEGER
C             The order of the pencil aS - bT.  N >= 0, even.
C
C     Z       (input/output) DOUBLE PRECISION array, dimension (LDZ, N)
C             On entry, the leading N-by-N part of this array must
C             contain the matrix Z.
C             On exit, if JOB = 'T', the leading N-by-N part of this
C             array contains the matrix Zout; otherwise, it contains the
C             matrix Z just before the application of the periodic QZ
C             algorithm. The entries in the rows N/2+1 to N and the
C             first N/2 columns are unchanged.
C
C     LDZ     INTEGER
C             The leading dimension of the array Z.  LDZ >= MAX(1, N).
C
C     B       (input/output) DOUBLE PRECISION array, dimension
C                            (LDB, N/2)
C             On entry, the leading N/2-by-N/2 part of this array must
C             contain the matrix B.
C             On exit, if JOB = 'T', the leading N/2-by-N/2 part of this
C             array contains the matrix Bout; otherwise, it contains the
C             matrix B just before the application of the periodic QZ
C             algorithm.
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
C             On exit, if JOB = 'T', the leading N/2-by-N/2 strictly
C             upper triangular part of the submatrix in the columns 2 to
C             N/2+1 of this array contains the strictly upper triangular
C             part of the skew-symmetric matrix Fout.
C             If JOB = 'E', the leading N/2-by-N/2 strictly upper
C             triangular part of the submatrix in the columns 2 to N/2+1
C             of this array contains the strictly upper triangular part
C             of the skew-symmetric matrix F just before the application
C             of the QZ algorithm.
C             The entries on the diagonal and the first superdiagonal of
C             this array are not referenced, but are assumed to be zero.
C             Moreover, the diagonal and the first subdiagonal of this
C             array on exit coincide to the corresponding diagonals of
C             this array on entry.
C
C     LDFG    INTEGER
C             The leading dimension of the array FG.
C             LDFG >= MAX(1, N/2).
C
C     Q       (output) DOUBLE PRECISION array, dimension (LDQ, N)
C             On exit, if COMPQ = 'I', the leading N-by-N part of this
C             array contains the orthogonal transformation matrix Q.
C             On exit, if COMPQ = 'N', the leading N-by-N part of this
C             array contains the orthogonal matrix Q1, such that
C
C                      (  Z11  Z12  )
C               Z*Q1 = (            ),
C                      (   0   Z22  )
C
C             where Z11 and Z22' are upper triangular (the first step
C             of the algorithm).
C
C     LDQ     INTEGER
C             The leading dimension of the array Q.  LDQ >= MAX(1, N).
C
C     U1      (input/output) DOUBLE PRECISION array, dimension
C                            (LDU1, N/2)
C             On entry, if COMPU = 'U', then the leading N/2-by-N/2 part
C             of this array must contain the upper left block of a
C             given matrix U0, and on exit, the leading N/2-by-N/2 part
C             of this array contains the updated upper left block U1 of
C             the product of the input matrix U0 and the transformation
C             matrix U used to transform the matrices Z and T.
C             On exit, if COMPU = 'I', then the leading N/2-by-N/2 part
C             of this array contains the upper left block U1 of the
C             orthogonal symplectic transformation matrix U.
C             If COMPU = 'N' this array is not referenced.
C
C     LDU1    INTEGER
C             The leading dimension of the array U1.
C             LDU1 >= 1,           if COMPU = 'N';
C             LDU1 >= MAX(1, N/2), if COMPU = 'I' or COMPU = 'U'.
C
C     U2      (input/output) DOUBLE PRECISION array, dimension
C                            (LDU2, N/2)
C             On entry, if COMPU = 'U', then the leading N/2-by-N/2 part
C             of this array must contain the upper right block of a
C             given matrix U0, and on exit, the leading N/2-by-N/2 part
C             of this array contains the updated upper right block U2 of
C             the product of the input matrix U0 and the transformation
C             matrix U used to transform the matrices Z and T.
C             On exit, if COMPU = 'I', then the leading N/2-by-N/2 part
C             of this array contains the upper right block U2 of the
C             orthogonal symplectic transformation matrix U.
C             If COMPU = 'N' this array is not referenced.
C
C     LDU2    INTEGER
C             The leading dimension of the array U2.
C             LDU2 >= 1,           if COMPU = 'N';
C             LDU2 >= MAX(1, N/2), if COMPU = 'I' or COMPU = 'U'.
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
C             complex conjugate pair, with ALPHAI(j+1) = -ALPHAI(j).
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
C     IWORK   INTEGER array, dimension (LIWORK)
C
C     LIWORK  INTEGER
C             The dimension of the array IWORK.  LIWORK >= N+9.
C
C     DWORK   DOUBLE PRECISION array, dimension (LDWORK)
C             On exit, if INFO = 0, DWORK(1) returns the optimal LDWORK.
C             On exit, if INFO = -23, DWORK(1) returns the minimum
C             value of LDWORK.
C
C     LDWORK  INTEGER
C             The dimension of the array DWORK.
C             If JOB = 'E' and COMPQ = 'N' and COMPU = 'N',
C                   LDWORK >= 3/4*N**2+MAX(N, 24)+3;
C             else, LDWORK >= 3/2*N**2+MAX(N, 24)+3.
C             For good performance LDWORK should generally be larger.
C
C             If LDWORK = -1, then a workspace query is assumed; the
C             routine only calculates the optimal size of the DWORK
C             array, returns this value as the first entry of the DWORK
C             array, and no error message related to LDWORK is issued by
C             XERBLA.
C
C     Error Indicator
C
C     INFO    INTEGER
C             = 0: succesful exit;
C             < 0: if INFO = -i, the i-th argument had an illegal value;
C             = 1: problem during computation of the eigenvalues;
C             = 2: periodic QZ algorithm did not converge in the SLICOT
C                  Library routine MB03BD.
C
C     METHOD
C
C     The algorithm uses Givens rotations and Householder reflections to
C     annihilate elements in Z and T such that Z is in a special block
C     triangular form and T is in skew-Hamiltonian Hessenberg form:
C
C         (  Z11  Z12  )      (  B1  F1  )                  
C     Z = (            ), T = (          ),
C         (   0   Z22  )      (   0  B1' )
C
C     with Z11 and Z22' upper triangular and B1 upper Hessenberg.
C     Subsequently, the periodic QZ algorithm is applied to the pencil
C     aZ22' Z11 - bB1 to determine orthogonal matrices Q1, Q2 and U such
C     that U' Z11 Q1, Q2' Z22' U are upper triangular and Q2' B1 Q1 is
C     upper quasi-triangular. See also page 35 in [1] for more details.
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
C     M. Voigt, Technische Universitaet Chemnitz, Apr. 2009.
C     V. Sima, Research Institute for Informatics, Bucharest, Aug. 2009.
C
C     REVISIONS
C
C     V. Sima, Dec. 2010, Jan. 2011, Aug. 2011, Nov. 2011, Jul. 2012,
C     Jul. 2013.
C     M. Voigt, Jan. 2012, Jul. 2013.
C
C     KEYWORDS
C
C     Periodic QZ algorithm, upper (quasi-)triangular matrix,
C     skew-Hamiltonian/skew-Hamiltonian pencil.
C
C     ******************************************************************
C
C     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
C
C     .. Scalar Arguments ..
      CHARACTER          COMPQ, COMPU, JOB
      INTEGER            INFO, LDB, LDFG, LDQ, LDU1, LDU2, LDWORK, LDZ,
     $                   LIWORK, N
C
C     .. Array Arguments ..
      INTEGER            IWORK( * )
      DOUBLE PRECISION   ALPHAI( * ), ALPHAR( * ), B( LDB, * ),
     $                   BETA( * ), DWORK( * ), FG( LDFG, * ),
     $                   Q( LDQ, * ), U1( LDU1, * ), U2( LDU2, * ),
     $                   Z( LDZ, * )
C
C     .. Local Scalars ..
      LOGICAL            LCMPQ, LCMPU, LINIU, LQUERY, LTRI, LUPDU
      CHARACTER*120      CMPQ, CMPSC
      INTEGER            I, IB1, ICF, ICG, IQ1, IQ2, ITAU, IU, IWARN,
     $                   IWRK, IZ11, IZ22, J, K, M, MINDW, MJ1, MJ2,
     $                   MJ3, MM, OPTDW
      DOUBLE PRECISION   BASE, CO, EMAX, EMIN, SI, TMP1, TMP2
C
C     .. Local Arrays ..
      INTEGER            IDUM( 1 )
      DOUBLE PRECISION   DUM(  4 )
C
C     .. External Functions ..
      LOGICAL            LSAME
      DOUBLE PRECISION   DLAMCH
      EXTERNAL           DLAMCH, LSAME
C
C     .. External Subroutines ..
      EXTERNAL           DCOPY, DGEMM, DGEQRF, DGERQF, DLACPY, DLARTG,
     $                   DLASET, DORMQR, DORMRQ, DROT, DSWAP, MA02AD,
     $                   MB01KD, MB01LD, MB01MD, MB03BD, XERBLA
C
C     .. Intrinsic Functions ..
      INTRINSIC          INT, MAX, MIN, MOD
C
C     .. Executable Statements ..
C
C     Decode the input arguments
C
      M  = N/2
      MM = M*M
      LTRI  = LSAME( JOB,   'T' )
      LCMPQ = LSAME( COMPQ, 'I' )
      LINIU = LSAME( COMPU, 'I' )
      LUPDU = LSAME( COMPU, 'U' )
      LCMPU = LINIU .OR. LUPDU
      IF( N.EQ.0 ) THEN
         MINDW = 1
      ELSE IF( LTRI .OR. LCMPQ .OR. LCMPU ) THEN
         MINDW = 6*MM + MAX( N, 24 ) + 3
      ELSE
         MINDW = 3*MM + MAX( N, 24 ) + 3
      END IF
      LQUERY = LDWORK.EQ.-1
C
C     Test the input arguments.
C
      INFO = 0
      IF( .NOT.( LSAME( JOB, 'E' ) .OR. LTRI ) ) THEN
         INFO = -1
      ELSE IF( .NOT.( LSAME( COMPQ, 'N' ) .OR. LCMPQ ) ) THEN
         INFO = -2
      ELSE IF( .NOT.( LSAME( COMPU, 'N' ) .OR. LCMPU ) ) THEN
         INFO = -3
      ELSE IF( N.LT.0 .OR. MOD( N, 2 ).NE.0 ) THEN
         INFO = -4
      ELSE IF(  LDZ.LT.MAX( 1, N ) ) THEN
         INFO = -6
      ELSE IF(  LDB.LT.MAX( 1, M ) ) THEN
         INFO = -8
      ELSE IF( LDFG.LT.MAX( 1, M ) ) THEN
         INFO = -10
      ELSE IF(  LDQ.LT.MAX( 1, N ) ) THEN
         INFO = -12
      ELSE IF( LDU1.LT.1 .OR. ( LCMPU .AND. LDU1.LT.M ) ) THEN
         INFO = -14
      ELSE IF( LDU2.LT.1 .OR. ( LCMPU .AND. LDU2.LT.M ) ) THEN
         INFO = -16
      ELSE IF( LIWORK.LT.N+9 ) THEN
         INFO = -21
      ELSE IF( .NOT. LQUERY .AND. LDWORK.LT.MINDW ) THEN
         DWORK( 1 ) = MINDW
         INFO = -23
      END IF
C
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DGHFST', -INFO )
         RETURN
      ELSE IF( N.GT.0 ) THEN
         IF( LQUERY ) THEN
C
C           Compute optimal workspace.
C
            CALL DGEQRF( N, M, DWORK, N, DWORK, DUM, -1, INFO )
            CALL DORMQR( 'Right', 'No Transpose', N, N, M, DWORK, N,
     $                   DWORK, Q, LDQ, DUM( 2 ), -1, INFO )
            CALL DGERQF( M, M, Z, LDZ, DWORK, DUM( 3 ), -1, INFO )
            CALL DORMRQ( 'Right', 'Transpose', N, M, M, Z, LDZ, DWORK,
     $                   Q, LDQ, DUM( 4 ), -1, INFO )
            J = MAX( MAX( N*M + MAX( INT( DUM( 1 ) ), INT( DUM( 2 ) ) ),
     $                    INT( DUM( 3 ) ), INT( DUM( 4 ) ) ) + M, 3*MM )
            DWORK( 1 ) = MAX( MINDW, J )
            RETURN
         END IF
      END IF
C
C     Quick return if possible.
C
      IF( N.EQ.0 ) THEN
         DWORK( 1 ) = ONE
         RETURN
      END IF
C
C     Determine machine constants.
C
      BASE = DLAMCH( 'Base' )
      EMIN = DLAMCH( 'Minimum Exponent' )
      EMAX = DLAMCH( 'Largest Exponent' )
C
C     Initializations.
C
      CALL DLASET( 'Full', N, N, ZERO, ONE, Q, LDQ )
C
      IF( LINIU ) THEN
         CALL DLASET( 'Full', M, M, ZERO, ONE,  U1, LDU1 )
         CALL DLASET( 'Full', M, M, ZERO, ZERO, U2, LDU2 )
      END IF
C
C     STEP 1: By changing the elimination order in the classical RQ
C             decomposition, determine an orthogonal matrix Q1 such that
C
C                 (  Z11  Z12  )
C             Z = (            ) Q1',
C                 (   0   Z22  )
C
C             where Z11 and Z22' are upper triangular.
C             Update Q and T subsequently.
C
      ITAU = N*M  + 1
      IWRK = ITAU + M
      CALL MA02AD( 'Full', M, N, Z( M+1, 1 ), LDZ, DWORK, N )
C
C     Perform a QR decomposition ( Z21  Z22 )' = Q1*R1, and
C     update ( Z11  Z12 ).
C
C     Workspace:     need   N*M+N;
C                    prefer N*M+M+M*NB, where NB is the optimal
C                                       blocksize.
C
      CALL DGEQRF( N, M, DWORK, N, DWORK( ITAU ), DWORK( IWRK ),
     $             LDWORK-IWRK+1, INFO )
      OPTDW = MAX( MINDW, INT( DWORK( IWRK ) ) + IWRK - 1 )
      CALL DORMQR( 'Right', 'No Transpose', M, N, M, DWORK, N,
     $             DWORK( ITAU ), Z, LDZ, DWORK( IWRK ), LDWORK-IWRK+1,
     $             INFO )
C
C     Copy R1' to Z22 and set the strictly upper triangular part of Z22
C     to zero.
C
      CALL MA02AD( 'Upper', M, M, DWORK, N, Z( M+1, M+1 ), LDZ )
      IF( M.GT.1 )
     $   CALL DLASET( 'Upper', M-1, M-1, ZERO, ZERO, Z( M+1, M+2 ),
     $                LDZ )
C
C     Update Q.
C
      CALL DORMQR( 'Right', 'No Transpose', N, N, M, DWORK, N,
     $             DWORK( ITAU ), Q, LDQ, DWORK( IWRK ), LDWORK-IWRK+1,
     $             INFO )
      OPTDW = MAX( OPTDW, INT( DWORK( IWRK ) ) + IWRK - 1 )
C
      DO 10 I = 1, M
         CALL DSWAP( N, Q( 1, I ), 1, Q( 1, M+I ), 1 )
   10 CONTINUE     
C
      ITAU = 1
      IWRK = ITAU + M
C
C     Perform an RQ decomposition Z12 = R2*Q2.
C
C     Workspace:     need   N;
C                    prefer M+M*NB, where NB is the optimal blocksize.
C
      CALL DGERQF( M, M, Z( 1, M+1 ), LDZ, DWORK( ITAU ), DWORK( IWRK ),
     $             LDWORK-IWRK+1, INFO )
      OPTDW = MAX( OPTDW, INT( DWORK( IWRK ) ) + IWRK - 1 )
C
C     Update Q.
C
      CALL DORMRQ( 'Right', 'Transpose', N, M, M, Z( 1, M+1 ), LDZ,
     $             DWORK( ITAU ), Q, LDQ, DWORK( IWRK ), LDWORK-IWRK+1,
     $             INFO )
      OPTDW = MAX( OPTDW, INT( DWORK( IWRK ) ) + IWRK - 1 )
C
C     Exchange Z11 and Z12 and set the strictly lower triangular part
C     of Z11 to zero.
C
      DUM( 1 ) = ZERO
      DO 20 J = 1, M - 1
         CALL DSWAP(  M, Z( 1, J ), 1, Z( 1, M+J ), 1 )
         CALL DCOPY( M-J, DUM( 1 ), 0, Z( J+1, J ), 1 )
   20 CONTINUE
C
      CALL DSWAP( M, Z( 1, M ), 1, Z( 1, N ), 1 )
C
C     Apply the transformations to B and FG.
C
C     Workspace:     need   3*M*M.
C
C     Copy the strictly upper triangular part of F and the transpose of
C     the strictly lower triangular part of G to appropriate locations
C     in DWORK.
C
      ICF  = 1
      ICG  = ICF + MM
      IWRK = ICG + MM
      IF( M.GT.1 ) THEN
         CALL DLACPY( 'Upper', M-1, M-1, FG( 1, 3 ), LDFG,
     $                DWORK( ICF+M ), M )
         CALL MA02AD( 'Lower', M-1, M-1, FG( 2, 1 ), LDFG,
     $                DWORK( ICG+M ), M )
      END IF
C
C     Skew-symmetric updates to determine the new F.
C
C       Fnew :=   Q(m+1:n,m+1:n)'*B *Q(1:m,m+1:n)
C               - Q(1:m,m+1:n)'  *B'*Q(m+1:n,m+1:n)
C               + Q(m+1:n,m+1:n)'*F *Q(m+1:n,m+1:n)
C               - Q(1:m,m+1:n)'  *G *Q(1:m,m+1:n).
C
      CALL MB01LD( 'Upper', 'Transpose', M, M, ZERO, ONE, DWORK( ICF ),
     $             M, Q( M+1, M+1 ), LDQ, DWORK( ICF ), M,
     $             DWORK( IWRK ), LDWORK-IWRK+1, INFO )
      CALL MB01LD( 'Upper', 'Transpose', M, M, ONE, ONE, DWORK( ICF ),
     $             M, Q( 1, M+1 ), LDQ, DWORK( ICG ), M, DWORK( IWRK ),
     $             LDWORK-IWRK+1, INFO )
      CALL DGEMM(  'No Transpose', 'No Transpose', M, M, M, ONE, B, LDB,
     $             Q( 1, M+1 ), LDQ, ZERO, DWORK( IWRK ), M )
      CALL MB01KD( 'Upper', 'Transpose', M, M, ONE, Q( M+1, M+1 ), LDQ,
     $             DWORK( IWRK ), M, ONE, DWORK( ICF ), M, INFO )
C
C     Copy the strictly lower triangular part of G and the transpose of
C     the strictly upper triangular part of F to appropriate locations
C     in DWORK.
C
      IF( M.GT.1 ) THEN
         CALL DLACPY( 'Lower', M-1, M-1, FG( 2, 1 ), LDFG,
     $                DWORK( ICG+1 ), M )
         CALL MA02AD( 'Upper', M-1, M-1, FG( 1, 3 ), LDFG,
     $                DWORK( ICF+1 ), M )
      END IF
C
C     Skew-symmetric updates to determine the new G.
C
C       Gnew :=   Q(1:m,1:m)'  *B'*Q(m+1:n,1:m)
C               - Q(m+1:n,1:m)'*B *Q(1:m,1:m)
C               + Q(1:m,1:m)'  *G *Q(1:m,1:m)
C               - Q(m+1:n,1:m)'*F *Q(m+1:n,1:m).
C
      CALL MB01LD( 'Lower', 'Transpose', M, M, ZERO, ONE, DWORK( ICG ),
     $             M, Q, LDQ, DWORK( ICG ), M, DWORK( IWRK ),
     $             LDWORK-IWRK+1, INFO )
      CALL MB01LD( 'Lower', 'Transpose', M, M, ONE, ONE, DWORK( ICG ),
     $             M, Q( M+1, 1 ), LDQ, DWORK( ICF ), M, DWORK( IWRK ),
     $             LDWORK-IWRK+1, INFO )
      CALL DGEMM(  'Transpose', 'No Transpose', M, M, M, ONE, B, LDB,
     $             Q( M+1, 1 ), LDQ, ZERO, DWORK( IWRK ), M )
      CALL MB01KD( 'Lower', 'Transpose', M, M, ONE, Q, LDQ,
     $             DWORK( IWRK ), M, ONE, DWORK( ICG ), M, INFO )
C
C     Determine the new B.
C
C       Bnew :=   Q(m+1:n,m+1:n)'*B *Q(1:m,1:m)
C               - Q(1:m,m+1:n)'  *B'*Q(m+1:n,1:m)
C               + Q(m+1:n,m+1:n)'*F *Q(m+1:n,1:m)
C               - Q(1:m,m+1:n)'  *G *Q(1:m,1:m).
C
      DO 30 I = 1, M
         CALL MB01MD( 'Upper', M, ONE, FG( 1, 2 ), LDFG, Q( M+1, I ), 1,
     $                ZERO, DWORK( IWRK+( I-1 )*M ), 1 )
   30 CONTINUE
      IF( M.GT.1 )
     $   CALL DLACPY( 'Upper', M-1, M-1, DWORK( ICF+M ), M, FG( 1, 3 ),
     $                LDFG )
      CALL DGEMM(  'No Transpose', 'No Transpose', M, M, M, ONE, B, LDB,
     $             Q, LDQ, ONE, DWORK( IWRK ), M )
      CALL DGEMM(  'Transpose', 'No Transpose', M, M, M, ONE,
     $             Q( M+1, M+1 ), LDQ, DWORK( IWRK ), M, ZERO,
     $             DWORK( ICF ), M )
      DO 40 I = 1, M
         CALL MB01MD( 'Lower', M, ONE, FG, LDFG, Q( 1, I ), 1, ZERO,
     $                DWORK( IWRK+( I-1 )*M ), 1 )
   40 CONTINUE
      CALL DGEMM(  'Transpose', 'No Transpose', M, M, M, ONE, B, LDB,
     $             Q( M+1, 1 ), LDQ, ONE, DWORK( IWRK ), M )
      CALL DGEMM(  'Transpose', 'No Transpose', M, M, M, -ONE,
     $             Q( 1, M+1 ), LDQ, DWORK( IWRK ), M, ONE,
     $             DWORK( ICF ), M )
      IF( M.GT.1 )
     $   CALL DLACPY( 'Lower', M-1, M-1, DWORK( ICG+1 ), M, FG( 2, 1 ),
     $                LDFG )
      CALL DLACPY( 'Full', M, M, DWORK( ICF ), M, B, LDB )
C
C     STEP 2: Reduce T to skew-Hamiltonian Hessenberg form.
C
      DO 70 K = 1, M - 1
C
C        I. Annihilate T(m+k, k+1:m-1) as well as T(m+k+1:n-1, k),
C           i.e., G(k+1:m-1, k).
C
         DO 50 J = K + 1, M - 1
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
            CALL DROT( M-J-1, FG( MJ2, J+1 ), 1, FG( MJ2, J ), 1, CO, SI
     $                  )
            FG( J+1, K ) = TMP1
            CALL DROT( J-K-1, FG( J+1, K+1 ), LDFG, FG( J, K+1 ), LDFG,
     $                 CO, SI )
C
C           Update Z.
C
            CALL DROT( J, Z( 1, J+1 ), 1, Z( 1, J ), 1, CO, SI )
            TMP1          = -SI*Z( J+1, J+1 )
            Z( J+1, J+1 ) =  CO*Z( J+1, J+1 )
C
            IF( LCMPQ ) THEN
C
C              Update Q.
C
               CALL DROT( N, Q( 1, J+1 ), 1, Q( 1, J ), 1, CO, SI )
            END IF
C
C           Determine a Givens rotation to annihilate Z(j+1,j) from the
C           left.
C
            CALL DLARTG( Z( J, J ), TMP1, CO, SI, TMP2 )
C
C           Update Z.
C
            Z(   J, J ) = TMP2
            Z( J+1, J ) = ZERO
            CALL DROT( N-J, Z( J, J+1 ), LDZ, Z( J+1, J+1 ), LDZ, CO,
     $                 SI )
            CALL DROT( J, Z( M+J, M+1 ), LDZ, Z( M+J+1, M+1 ), LDZ,
     $                 CO, SI )
            TMP1              = SI*Z( M+J+1, M+J+1 )
            Z( M+J+1, M+J+1 ) = CO*Z( M+J+1, M+J+1 )
C
            IF( LCMPU ) THEN
C
C              Update U1 and U2.
C
               CALL DROT( M, U1( 1, J ), 1, U1( 1, J+1 ), 1, CO, SI )
               CALL DROT( M, U2( 1, J ), 1, U2( 1, J+1 ), 1, CO, SI )
            END IF
C
C           Determine a Givens rotation to annihilate Z(m+j,m+j+1) from
C           the right.
C
            CALL DLARTG( Z( M+J, M+J ), TMP1, CO, SI, TMP2 )
C
C           Update Z.
C
            CALL DROT( M, Z( 1, M+J ), 1, Z( 1, M+J+1 ), 1, CO, SI )
            Z( M+J, M+J ) = TMP2
            CALL DROT( M-J, Z( M+J+1, M+J ), 1, Z( M+J+1, M+J+1 ), 1,
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
   50    CONTINUE
C
C        II. Annihilate G(k,m) (and also G(m,k)).
C
C        Determine a Givens rotation to annihilate G(m,k) from the
C        left.
C
         CALL DLARTG( B( M, K ), -FG( M, K ), CO, SI, TMP1 )
C
C        Update B, F and G.
C
         CALL DROT( M-1, FG( 1, M+1 ), 1, B( 1, M ), 1, CO, SI )
         B( M, K ) = TMP1
         CALL DROT( M-K-1, FG( M, K+1 ), LDFG, B( M, K+1 ), LDB, CO,
     $              SI )
C
C        Update Z.
C
         CALL DROT( M, Z( 1, N ), 1, Z( 1, M ), 1, CO, SI )
         TMP1      = -SI*Z( N, N )
         Z( N, N ) =  CO*Z( N, N )
C
         IF( LCMPQ ) THEN
C
C           Update Q.
C
            CALL DROT( N, Q( 1, N ), 1, Q( 1, M ), 1, CO, SI )
         END IF
C
C        Determine a Givens rotation to annihilate Z(n,m) from the left.
C
         CALL DLARTG( Z( M, M ), TMP1, CO, SI, TMP2 )
C
C        Update Z.
C
         Z( M, M ) = TMP2
         CALL DROT( M, Z( M, M+1 ), LDZ, Z( N, M+1 ), LDZ, CO, SI )
C
         IF( LCMPU ) THEN
C
C           Update U1 and U2.
C           
            CALL DROT( M, U1( 1, M ), 1, U2( 1, M ), 1, CO, SI )
         END IF
C
C        III. Annihilate B(k+2:m,k).
C
         DO 60 J = M, K + 2, -1
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
            CALL DROT( J-2, FG( 1, J ), 1, FG( 1, J+1 ), 1, CO, SI )
            B( J-1, K ) = TMP1
            B(   J, K ) = ZERO
            CALL DROT( M-K, B( J-1, K+1 ), LDB, B( J, K+1 ), LDB, CO,
     $                 SI )
            CALL DROT( M-J, FG( J-1, MJ2 ), LDFG, FG( J, MJ2 ), LDFG,
     $                 CO, SI )
C
C           Update Z.
C
            CALL DROT( M, Z( 1, M+J-1 ), 1, Z( 1, M+J ), 1, CO, SI )
            TMP1              = -SI*Z( M+J-1, M+J-1 )
            Z( M+J-1, M+J-1 ) =  CO*Z( M+J-1, M+J-1 )
            CALL DROT( M-J+1, Z( M+J, M+J-1 ), 1, Z( M+J, M+J ), 1, CO,
     $                 SI )
C
            IF( LCMPQ ) THEN
C
C              Update Q.
C
               CALL DROT( N, Q( 1, M+J-1 ), 1, Q( 1, M+J ), 1, CO, SI )
            END IF
C
C           Determine a Givens rotation to annihilate Z(m+j-1,m+j) from
C           the left.
C
            CALL DLARTG( Z( M+J, M+J ), TMP1, CO, SI, TMP2 )
C
C           Update Z.
C
            TMP1          = SI*Z( J-1, J-1 )
            Z( J-1, J-1 ) = CO*Z( J-1, J-1 )
            CALL DROT( N-J+1, Z( J, J ), LDZ, Z( J-1, J ), LDZ, CO, SI )
            CALL DROT( J-1, Z( M+J, M+1 ), LDZ, Z( M+J-1, M+1 ), LDZ,
     $                 CO, SI )
            Z( M+J, M+J ) = TMP2
C
            IF( LCMPU ) THEN
C
C              Update U1 and U2.
C
               CALL DROT( M, U1( 1, J ), 1, U1( 1, J-1 ), 1, CO, SI )
               CALL DROT( M, U2( 1, J ), 1, U2( 1, J-1 ), 1, CO, SI )
            END IF
C
C           Determine a Givens rotation to annihilate Z(j,j-1) from the
C           right.
C
            CALL DLARTG( Z( J, J ), TMP1, CO, SI, TMP2 )
C
C           Update Z.
C
            Z( J, J ) = TMP2
            CALL DROT( J-1, Z( 1, J ), 1, Z( 1, J-1 ), 1, CO, SI )
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
   60    CONTINUE
   70 CONTINUE
C
C             (  Z11  Z12  )      (  B1  F1  )
C     Now Z = (            ), T = (          ),
C             (   0   Z22  )      (   0  B1' )
C
C     where Z11 and Z22' are upper triangular and B1 upper Hessenberg.
C
C     STEP 3: Apply the periodic QZ algorithm to the pencil
C             aZ22' Z11 - bB1 to determine orthogonal matrices
C             Q1, Q2 and U such that U' Z11 Q1 and Q2' Z22' U are upper
C             triangular, and Q2' B1 Q1 is upper quasi-triangular.
C
C     Determine mode of computations.
C
      IQ2  = 1
      IF( LTRI .OR. LCMPQ .OR. LCMPU ) THEN
         CMPQ = 'Initialize'
         IQ1  = IQ2 + MM
         IU   = IQ1 + MM
         IB1  = IU  + MM
      ELSE
         CMPQ = 'No Computation'
         IB1  = 1
      END IF
      IZ11 = IB1  + MM
      IZ22 = IZ11 + MM
      IWRK = IZ22 + MM
C
      IF( LTRI ) THEN
         CMPSC = 'Schur Form'
      ELSE
         CMPSC = 'Eigenvalues Only'
      END IF
C
      IWORK( 1 ) =  1
      IWORK( 2 ) = -1
      IWORK( 3 ) = -1
C
      CALL DLACPY( 'Full', M, M, B, LDB, DWORK( IB1 ),  M )
      CALL DLACPY( 'Full', M, M, Z, LDZ, DWORK( IZ11 ), M )
      CALL MA02AD( 'Full', M, M, Z( M+1, M+1 ), LDZ, DWORK( IZ22 ), M )
C
C     Periodic QZ iteration.
C  
C     Real workspace:    need   w1 + MAX( N,24 ) + 3, where
C                               w1 = 6*M**2, if JOB = 'T', or
C                                            COMPQ = 'I' or COMPU <>'N';
C                               w1 = 3*M**2, otherwise.
C     Integer workspace: need   N + 9.
C
      CALL MB03BD( CMPSC, 'Careful', CMPQ, IDUM, 3, M, 1, 1, M, IWORK,
     $             DWORK( IB1 ), M, M, DWORK( IQ2 ), M, M, ALPHAR,
     $             ALPHAI, BETA, IWORK( 4 ), IWORK( M+4 ),
     $             LIWORK-( M+3 ), DWORK( IWRK ), LDWORK-IWRK+1, IWARN,
     $             INFO )
      IF( IWARN.GT.0 .AND. IWARN.LT.M ) THEN
         INFO = 1
         RETURN
      ELSE IF( INFO.GT.0 ) THEN
         INFO = 2
         RETURN
      END IF
C
C     Compute the eigenvalues of the pencil aS - bT.
C
      DO 80 I = 1, M
         IF( IWORK( I+3 ).GE.EMIN .AND. IWORK( I+3 ).LE.EMAX ) THEN
C
C           b = BASE**IWORK(i+3) is between underflow and overflow
C           threshold, BETA(i) is divided by b.
C
            BETA( I ) = BETA( I )/BASE**IWORK( I+3 )
         ELSE
C
C           Error.
C
            INFO = 1
            RETURN
         END IF
   80 CONTINUE
C
      IF( LTRI ) THEN
C
C        Update Z.
C
C        Workspace:     need   w2 = 5*M**2, if JOB = 'T';
C                              w2 = 2*M**2, otherwise.
C
         IWRK = IZ11
         CALL DLACPY( 'Upper', M, M, DWORK( IZ11 ), M, Z, LDZ )
         CALL DGEMM(  'Transpose', 'No Transpose', M, M, M, ONE,
     $                DWORK( IU ), M, Z( 1, M+1 ), LDZ, ZERO,
     $                DWORK( IWRK ), M )
         CALL DGEMM(  'No Transpose', 'No Transpose', M, M, M, ONE,
     $                DWORK( IWRK ), M, DWORK( IQ2 ), M, ZERO,
     $                Z( 1, M+1 ), LDZ )
         CALL MA02AD( 'Upper', M, M, DWORK( IZ22 ), M, Z( M+1, M+1 ),
     $                LDZ )
C
C        Update B.
C
         CALL DLACPY( 'Full', M, M, DWORK( IB1 ), M, B, LDB )
         IWRK = IB1
C
C        Skew-symmetric update of F.
C
C        Workspace:     need   w3 + M, where
C                              w3 = 3*M**2, if JOB = 'T', or
C                                           COMPQ = 'I' or COMPU <>'N';
C                              w3 = 0,      otherwise;
C                       prefer w3 + M*(M-1).
C
         CALL MB01LD( 'Upper', 'Transpose', M, M, ZERO, ONE, FG( 1, 2 ),
     $                LDFG, DWORK( IQ2 ), M, FG( 1, 2 ), LDFG,
     $                DWORK( IWRK ), LDWORK-IWRK+1, INFO )
C
         IF( LCMPQ ) THEN
C
C           Update Q.
C
C           Workspace:     need   w3 + N*M.
C
            CALL DGEMM(  'No Transpose', 'No Transpose', N, M, M, ONE,
     $                   Q, LDQ, DWORK( IQ1 ), M, ZERO, DWORK( IWRK ),
     $                   N )
            CALL DLACPY( 'Full', N, M, DWORK( IWRK ), N, Q, LDQ )
            CALL DGEMM(  'No Transpose', 'No Transpose', N, M, M, ONE,
     $                   Q( 1, M+1 ), LDQ, DWORK( IQ2 ), M, ZERO,
     $                   DWORK( IWRK ), N )
            CALL DLACPY( 'Full', N, M, DWORK( IWRK ), N, Q( 1, M+1 ),
     $                   LDQ )
         END IF
C
         IF( LCMPU ) THEN
C
C           Update U.
C
            CALL DGEMM(  'No Transpose', 'No Transpose', M, M, M, ONE,
     $                   U1, LDU1, DWORK( IU ), M, ZERO, DWORK( IWRK ),
     $                   M )
            CALL DLACPY( 'Full', M, M, DWORK( IWRK ), M, U1, LDU1 )
            CALL DGEMM(  'No Transpose', 'No Transpose', M, M, M, ONE,
     $                   U2, LDU2, DWORK( IU ), M, ZERO, DWORK( IWRK ),
     $                   M )
            CALL DLACPY( 'Full', M, M, DWORK( IWRK ), M, U2, LDU2 )
         END IF
      END IF
C
      DWORK( 1 ) = OPTDW
      RETURN
C *** Last line of DGHFST ***
      END 
