      SUBROUTINE PSB02MD( DICO, HINV, UPLO, SCAL, SORT, N, A, IA, JA,
     $                    DESCA, G, IG, JG, DESCG, Q, IQ, JQ, DESCQ, 
     $                    RCOND, WR, WI, S, IS, JS, DESCS, U, IU, JU,
     $                    DESCU, DWORK, LDWORK, IWORK, LIWORK, INFO,
     $                    TIMES )
C   
C  -- ScaLAPACK/PSLICOT-style routine --
C     Preliminary version.
C     Dept. Computing Science and HPC2N, Univ. of Umeå, Sweden
C     December 7, 2007.
C
      IMPLICIT NONE
C
C     .. Scalar Arguments ..
      CHARACTER         DICO, HINV, SCAL, SORT, UPLO
      INTEGER           INFO, IA, JA, IG, JG, IQ, JQ, IS, JS, IU, JU,
     $                  LDWORK, LIWORK, N
      DOUBLE PRECISION  RCOND
C     .. Array Arguments ..
      INTEGER           DESCA(*), DESCG(*), DESCQ(*), DESCS(*), 
     $                  DESCU(*), IWORK(*)
      DOUBLE PRECISION  A(*), DWORK(*), G(*), Q(*), S(*), U(*), WR(*), 
     $                  WI(*) 
C
C     PURPOSE
C
C     To solve for X either the continuous-time algebraic Riccati
C     equation
C                              -1
C        Q + A'*X + X*A - X*B*R  B'*X = 0                            (1)
C
C     or the discrete-time algebraic Riccati equation
C                                        -1
C        X = A'*X*A - A'*X*B*(R + B'*X*B)  B'*X*A + Q                (2)
C
C     where A, B, Q and R are N-by-N, N-by-M, N-by-N and M-by-M matrices
C     respectively, with Q symmetric and R symmetric nonsingular; X is
C     an N-by-N symmetric matrix.
C                       -1
C     The matrix G = B*R  B' must be provided on input, instead of B and
C     R, that is, for instance, the continuous-time equation
C
C        Q + A'*X + X*A - X*G*X = 0                                  (3)
C
C     is solved, where G is an N-by-N symmetric matrix. The routine 
C     PSB02MT should be used to compute G, given B and R.
C
C     The routine also returns the computed values of the closed-loop
C     spectrum of the optimal system, i.e., the stable eigenvalues
C     lambda(1),...,lambda(N) of the corresponding Hamiltonian or
C     symplectic matrix associated to the optimal problem.
C
C     ARGUMENTS
C
C     Mode Parameters
C
C     DICO    CHARACTER*1
C             Specifies the type of Riccati equation to be solved as
C             follows:
C             = 'C':  Equation (3), continuous-time case;
C             = 'D':  Equation (2), discrete-time case.
C
C     HINV    CHARACTER*1
C             If DICO = 'D', specifies which symplectic matrix is to be
C             constructed, as follows:
C             = 'D':  The matrix H in (5) (see METHOD) is constructed;
C             = 'I':  The inverse of the matrix H in (5) is constructed.
C             HINV is not used if DICO = 'C'.
C
C     UPLO    CHARACTER*1
C             Specifies which triangle of the matrices G and Q is
C             stored, as follows:
C             = 'U':  Upper triangle is stored;
C             = 'L':  Lower triangle is stored.
C
C     SCAL    CHARACTER*1
C             Specifies whether or not a scaling strategy should be
C             used, as follows:
C             = 'G':  General scaling should be used;
C             = 'N':  No scaling should be used.
C
C     SORT    CHARACTER*1
C             Specifies which eigenvalues should be obtained in the top
C             of the Schur form, as follows:
C             = 'S':  Stable   eigenvalues come first;
C             = 'U':  Unstable eigenvalues come first.
C
C     Input/Output Parameters
C
C     N       (global input) INTEGER
C             The order of the matrices A, Q, G and X.  N >= 0.
C
C     A       (local input/output) DOUBLE PRECISION pointer into the
C             local memory to an array, dimension (LLD_A, LOCq(N)).
C             This array must contain the local pieces of the N-by-N
C             distributed matrix A.
C             On exit, if DICO = 'D', and INFO = 0 or INFO > 1, this 
C             array contains the matrix inv(A).
C             Otherwise, the array A is unchanged on exit.
C
C     IA      (global input) INTEGER
C     JA      (global input) INTEGER
C             The row and column index in the global array A indicating 
C             the first column of sub( A ). IA = JA = 1 must hold for
C             this version.
C
C     DESCA   (global and local input) INTEGER array, dimension DLEN_
C             The array descriptor for the distributed matrix A.
C
C     G       (local input) DOUBLE PRECISION pointer into the
C             local memory to an array, dimension (LLD_G, LOCq(N)).
C             This array must contain the local pieces of the N-by-N
C             distributed matrix G.
C             Matrix G must be symmetric positive definite (SPD).
C
C     IG      (global input) INTEGER
C     JG      (global input) INTEGER
C             The row and column index in the global array G indicating 
C             the first column of sub( G ). IG = JG = 1 must hold for
C             this version.
C
C     DESCG   (global and local input) INTEGER array, dimension DLEN_
C             The array descriptor for the distributed matrix G.
C
C     Q       (local input/output) DOUBLE PRECISION pointer into the
C             local memory to an array, dimension (LLD_Q,LOCq(N))
C             This array must contain the local pieces of the N-by-N
C             distributed matrix Q.
C             Matrix Q must be symmetric positive definite (SPD).
C             On exit, if INFO = 0, this array contains the local pieces
C             solution of the globally distributed solution matrix X.
C
C     IQ      (global input) INTEGER
C     JQ      (global input) INTEGER
C             The row and column index in the global array Q indicating 
C             the first column of sub( Q ). IQ = JQ = 1 must hold for
C             this version.
C
C     DESCQ   (global and local input) INTEGER array, dimension DLEN_
C             The array descriptor for the distributed matrix Q.
C
C     RCOND   (output) DOUBLE PRECISION
C             An estimate of the reciprocal of the condition number (in
C             the 1-norm) of the N-th order system of algebraic
C             equations from which the solution matrix X is obtained.
C
C     WR      (output) DOUBLE PRECISION array, dimension (2*N)
C     WI      (output) DOUBLE PRECISION array, dimension (2*N)
C             If INFO = 0 or INFO = 5, these arrays contain the real and
C             imaginary parts, respectively, of the eigenvalues of the
C             2N-by-2N matrix S, ordered as specified by SORT (except
C             for the case HINV = 'D', when the order is opposite to
C             that specified by SORT). The leading N elements of these
C             arrays contain the closed-loop spectrum of the system
C                           -1                         
C             matrix A - B*R  *B'*X, if DICO = 'C', or of the matrix
C                               -1      
C             A - B*(R + B'*X*B)  B'*X*A, if DICO = 'D'. Specifically,
C                lambda(k) = WR(k) + j*WI(k), for k = 1,2,...,N.
C
C     S       (local output) DOUBLE PRECISION pointer into the
C             local memory to an array, dimension (LLD_S,LOCq(2*N))
C             If INFO = 0 or INFO = 5, the this array contains the 
C             ordered real Schur form S of the Hamiltonian or 
C             symplectic matrix H. That is,
C
C                    (S   S  )
C                    ( 11  12)
C                S = (       ),
C                    (0   S  )
C                    (     22)
C
C             where S  , S   and S   are N-by-N matrices.
C                    11   12      22
C
C     IS      (global input) INTEGER
C     JS      (global input) INTEGER
C             The row and column index in the global array S indicating 
C             the first column of sub( S ). IS = JS = 1 must hold for
C             this version.
C
C     DESCS   (global and local input) INTEGER array, dimension DLEN_
C             The array descriptor for the distributed matrix S.
C
C     U       (global output) DOUBLE PRECISION pointer into the
C             local memory to an array, dimension (LLD_S,LOCq(2*N))
C             If INFO = 0 or INFO = 5, this array contains the 
C             transformation matrix U which reduces the Hamiltonian or 
C             symplectic matrix H to the ordered real Schur form S. That 
C             is,
C
C                    (U   U  )
C                    ( 11  12)
C                U = (       ),
C                    (U   U  )
C                    ( 21  22)
C
C             where U  , U  , U   and U   are N-by-N matrices.
C                    11   12   21      22
C
C     IU      (global input) INTEGER
C     JU      (global input) INTEGER
C             The row and column index in the global array U indicating 
C             the first column of sub( U ). IU = JU = 1 must hold for
C             this version.
C
C     DESCU   (global and local input) INTEGER array, dimension DLEN_
C             The array descriptor for the distributed matrix U.
C
C     Workspace
C
C     DWORK   (local workspace) DOUBLE PRECISION array, dimension (LDWORK)
C             On exit, if INFO = 0, DWORK(1) returns the optimal value
C             of LDWORK and DWORK(2) returns the scaling factor used
C             (set to 1 if SCAL = 'N'), also set if INFO = 5;
C             if DICO = 'D', DWORK(3) returns the reciprocal condition
C             number of the given matrix  A.
C
C     LDWORK  (local input) INTEGER
C             The dimension of the real array used for workspace.
C             LDWORK >= if DICO = 'C';
C             LDWORK >= if DICO = 'D'.
C
C             If LDWORK = -1, then LDWORK is global input and a
C             workspace query is assumed; the routine only calculates
C             the minimum and optimal size for all work arrays. Each of
C             these values is returned in the first entry of the
C             corresponding work array, and no error message is issued
C             by PXERBLA.
C
C     IWORK   (local workspace) INTEGER array, dimension (LIWORK)
C             On exit, IWORK(1) returns the minimal and optimal LIWORK.
C
C     LIWORK  (local input) INTEGER
C             The dimension of the integer array used for workspace.
C             LIWORK >=  if DICO = 'C';
C             LIWORK >=  if DICO = 'D'.             
C
C             If LIWORK = -1, then LIWORK is global input and a
C             workspace query is assumed; the routine only calculates
C             the minimum and optimal size for all work arrays. Each
C             of these values is returned in the first entry of the
C             corresponding work array, and no error message is issued
C             by PXERBLA.
C
C     Error Indicator
C
C     INFO    (global output) INTEGER
C             = 0: Successful exit.
C
C     Error Indicator
C
C     INFO    INTEGER
C             = 0:  successful exit;
C             < 0:  If the i-th argument is an array and the j-entry
C                   had an illegal value, then INFO = -(i*100+j), if the
C                   i-th argument is a scalar and had an illegal value,
C                   then INFO = -i.
C             = 1:  if matrix A is (numerically) singular in discrete-
C                   time case;
C             = 2:  if the Hamiltonian or symplectic matrix H cannot be
C                   reduced to real Schur form;
C             = 3:  if the real Schur form of the Hamiltonian or
C                   symplectic matrix H cannot be appropriately ordered;
C             = 4:  if the Hamiltonian or symplectic matrix H has less
C                   than N stable eigenvalues;
C             = 5:  if the N-th order system of linear algebraic
C                   equations, from which the solution matrix X would
C                   be obtained, is singular to working precision.
C
C     METHOD
C
C     The method used is the Schur vector approach proposed by Laub.
C     It is assumed that [A,B] is a stabilizable pair (where for (3) B
C     is any matrix such that B*B' = G with rank(B) = rank(G)), and
C     [E,A] is a detectable pair, where E is any matrix such that
C     E*E' = Q with rank(E) = rank(Q). Under these assumptions, any of
C     the algebraic Riccati equations (1)-(3) is known to have a unique
C     non-negative definite solution. See [2].
C     Now consider the 2N-by-2N Hamiltonian or symplectic matrix
C
C                 ( A   -G )
C            H =  (        ),                                    (4)
C                 (-Q   -A'),
C
C     for continuous-time equation, and
C                    -1        -1
C                 ( A         A  *G   )
C            H =  (   -1          -1  ),                         (5)
C                 (Q*A    A' + Q*A  *G)
C                                                            -1
C     for discrete-time equation, respectively, where G = B*R  *B'.
C     The assumptions guarantee that H in (4) has no pure imaginary
C     eigenvalues, and H in (5) has no eigenvalues on the unit circle.
C     If Y is an N-by-N matrix then there exists an orthogonal matrix U
C     such that U'*Y*U is an upper quasi-triangular matrix. Moreover, U
C     can be chosen so that the 2-by-2 and 1-by-1 diagonal blocks
C     (corresponding to the complex conjugate eigenvalues and real
C     eigenvalues respectively) appear in any desired order. This is the
C     ordered real Schur form. Thus, we can find an orthogonal
C     similarity transformation U which puts (4) or (5) in ordered real
C     Schur form
C
C            U'*H*U = S = (S(1,1)  S(1,2))
C                         (  0     S(2,2))
C
C     where S(i,j) is an N-by-N matrix and the eigenvalues of S(1,1)
C     have negative real parts in case of (4), or moduli greater than
C     one in case of (5). If U is conformably partitioned into four
C     N-by-N blocks
C
C               U = (U(1,1)  U(1,2))
C                   (U(2,1)  U(2,2))
C
C     with respect to the assumptions we then have
C     (a) U(1,1) is invertible and X = U(2,1)*inv(U(1,1)) solves (1),
C         (2), or (3) with X = X' and non-negative definite;
C     (b) the eigenvalues of S(1,1) (if DICO = 'C') or S(2,2) (if
C         DICO = 'D') are equal to the eigenvalues of optimal system
C         (the 'closed-loop' spectrum).
C
C     [A,B] is stabilizable if there exists a matrix F such that (A-BF)
C     is stable. [E,A] is detectable if [A',E'] is stabilizable.
C
C     REFERENCES
C
C     [1] Laub, A.J.
C         A Schur Method for Solving Algebraic Riccati equations.
C         IEEE Trans. Auto. Contr., AC-24, pp. 913-921, 1979.
C
C     [2] Wonham, W.M.
C         On a matrix Riccati equation of stochastic control.
C         SIAM J. Contr., 6, pp. 681-697, 1968.
C
C     [3] Sima, V.
C         Algorithms for Linear-Quadratic Optimization.
C         Pure and Applied Mathematics: A Series of Monographs and
C         Textbooks, vol. 200, Marcel Dekker, Inc., New York, 1996.
C
C     NUMERICAL ASPECTS
C                               3
C     The algorithm requires 0(N ) operations.
C
C     PARALLEL EXECUTION RECOMMENDATIONS
C
C     Restrictions:
C      o All matrices should be be completely aligned, i.e.,
C        having the same MB, NB, RSRC, CSRC, etc. 
C
C     FURTHER COMMENTS
C
C     To obtain a stabilizing solution of the algebraic Riccati
C     equation for DICO = 'D', set SORT = 'U', if HINV = 'D', or set
C     SORT = 'S', if HINV = 'I'.
C
C     The routine can also compute the anti-stabilizing solutions of
C     the algebraic Riccati equations, by specifying 
C         SORT = 'U' if DICO = 'D' and HINV = 'I', or DICO = 'C', or 
C         SORT = 'S' if DICO = 'D' and HINV = 'D'.
C
C     Usually, the combinations HINV = 'D' and SORT = 'U', or HINV = 'I'
C     and SORT = 'U', will be faster then the other combinations [3].
C
C     CONTRIBUTOR
C
C     Robert Granat, Dept Computing Science and HPC2N, Dec. 2007.
C
C     REVISIONS
C
C     None.
C
C     KEYWORDS
C
C     Algebraic Riccati equation, closed loop system, continuous-time
C     system, discrete-time system, optimal regulator, Schur form.
C
C     ******************************************************************
C
C     .. Parameters ..
      LOGICAL           DEBUG, PRINT
      INTEGER           BLOCK_CYCLIC_2D, CSRC_, CTXT_, DLEN_, DT_,
     $                  LLD_, MB_, M_, NB_, N_, RSRC_
      PARAMETER         ( BLOCK_CYCLIC_2D = 1, DLEN_ = 9, DT_ = 1,
     $                    CTXT_ = 2, M_ = 3, N_ = 4, MB_ = 5, NB_ = 6,
     $                    RSRC_ = 7, CSRC_ = 8, LLD_ = 9,
     $                    DEBUG = .FALSE., PRINT = .FALSE. )
      DOUBLE PRECISION  ZERO, HALF, ONE
      PARAMETER         ( ZERO = 0.0D0, HALF = 0.5D0, ONE = 1.0D0 )
C     .. Local Scalars ..
      LOGICAL           DISCR, LHINV, LSCAL, LSORT, LUPLO, LQUERY
      INTEGER           I, IERR, ISCL, N2, NP1, NROT, ICTXT, NPROW,
     $                  NPCOL, NIPIV, MYROW, MYCOL
      DOUBLE PRECISION  GNORM, QNORM, RCONDA, UNORM, WRKOPT, IWRKOPT
C     .. Local Arrays ..
      DOUBLE PRECISION  TIMES(3)
C     .. External Functions ..
      LOGICAL           LSAME, SB02MR, SB02MS, SB02MV, SB02MW
      INTEGER           NUMROC
      DOUBLE PRECISION  DLAMCH, PDLANGE, PDLANSY
      EXTERNAL          DLAMCH, PDLANGE, PDLANSY, LSAME, SB02MR, SB02MS, 
     $                  SB02MV, SB02MW, NUMROC
C     .. External Subroutines ..
      EXTERNAL          PDAXPY, PDCOPY, PDGECON, PDGEES, PDGETRF, 
     $                  PDGETRS, PDLAMVE, PDLASCL, PDLASET, PDSCAL, 
     $                  PDSWAP, PSB02MU, PXERBLA, PCHK2MAT
C     .. Intrinsic Functions ..
      INTRINSIC         INT, MAX
C     .. Executable Statements ..
C
C     Get grid parameters and initialize some variables.
C
      IF( DEBUG ) WRITE(*,*) 'PSB02MD'
      INFO = 0
      LQUERY = ( ( LDWORK.EQ.-1 ).OR.( LIWORK.EQ.-1 ) )
      ICTXT = DESCA( CTXT_ )
      CALL BLACS_GRIDINFO( ICTXT, NPROW, NPCOL, MYROW, MYCOL )
      IF( NPROW.EQ.-1 ) INFO = -(600+CTXT_)
      IF( INFO.EQ.0 ) THEN
         N2  = N + N
         NP1 = N + 1
         DISCR = LSAME( DICO, 'D' )
         LSCAL = LSAME( SCAL, 'G' )
         LSORT = LSAME( SORT, 'S' )
         LUPLO = LSAME( UPLO, 'U' )
         IF ( DISCR ) LHINV = LSAME( HINV, 'D' )
C     
C        Test the input scalar arguments.
C     
         IF( .NOT.DISCR .AND. .NOT.LSAME( DICO, 'C' ) ) THEN
            INFO = -1
         ELSE IF( DISCR ) THEN
            IF( .NOT.LHINV .AND. .NOT.LSAME( HINV, 'I' ) )
     $           INFO = -2
         END IF
         IF( .NOT.LUPLO .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
            INFO = -3
         ELSE IF( .NOT.LSCAL .AND. .NOT.LSAME( SCAL, 'N' ) ) THEN
            INFO = -4
         ELSE IF( .NOT.LSORT .AND. .NOT.LSAME( SORT, 'U' ) ) THEN
            INFO = -5
         ELSE IF( N.LT.0 ) THEN
            INFO = -6
         END IF
C     
C        Test input matrices
C     
         IF( INFO.EQ.0 ) THEN
            CALL CHK1MAT( N, 6, N, 6, IA, JA, DESCA, 10, INFO )
            CALL CHK1MAT( N, 6, N, 6, IG, JG, DESCG, 14, INFO )
            CALL CHK1MAT( N, 6, N, 6, IQ, JQ, DESCQ, 18, INFO )
c            CALL CHK1MAT( N2, 6, N2, 6, IS, JS, DESCS, 25, INFO )
c            CALL CHK1MAT( N2, 6, N2, 6, IU, JU, DESCU, 29, INFO )
c            CALL PCHK2MAT( N2, 6, N2, 6, IS, JS, DESCS, 25, N2, 6, N2, 
c     $                     6, IU, JU, DESCU, 29, 0, IWORK, IWORK, INFO )
         END IF
C
C        Test workspace - TODO!
C
         IF( INFO.EQ.0 ) THEN
            WRKOPT = N
            IWRKOPT = N
            IF( LDWORK.LT.WRKOPT ) THEN
               INFO = -31
            ELSEIF( LIWORK.LT.IWRKOPT ) THEN
               INFO = -33
            END IF
         END IF
C
C        Error or workspace check return.
C
         IF( .NOT.LQUERY .AND. INFO.NE.0 ) THEN
            CALL PXERBLA( 'PSB02MD', -INFO )
            RETURN
         ELSEIF( LQUERY ) THEN
            DWORK( 1 ) = WRKOPT
            IWORK( 1 ) = IWRKOPT
            RETURN
         END IF
      END IF
C
C     Quick return if possible.
C
      IF ( N.EQ.0 ) THEN
         RCOND = ONE
         DWORK(1) = ONE
         DWORK(2) = ONE
         IF ( DISCR ) DWORK(3) = ONE
         RETURN
      END IF
C
      IF ( LSCAL ) THEN
C
C        Compute the norms of the matrices Q and G.     
C
         QNORM = PDLANSY( '1-norm', UPLO, N, Q, IQ, JQ, DESCQ, DWORK )
         GNORM = PDLANSY( '1-norm', UPLO, N, G, IG, JG, DESCG, DWORK ) 
      END IF
C
C     Initialize the Hamiltonian or symplectic matrix associated with
C     the problem.
C
      IF( DEBUG ) WRITE(*,*) 'PSB02MU in...'
      CALL PSB02MU( DICO, HINV, UPLO, N, A, IA, JA, DESCA, G, IG, JG, 
     $              DESCG, Q, IQ, JQ, DESCQ, S, IS, JS, DESCS, DWORK, 
     $              LDWORK, IWORK, LIWORK, INFO )
      IF ( INFO.NE.0 ) THEN
         INFO = 1
         RETURN
      END IF
C
      IF( PRINT ) THEN
         CALL PDLAPRNT( N2, N2, S, IS, JS, DESCS, 0, 0,
     $                 'S0', 6, DWORK )
      END IF

C
      WRKOPT = DWORK(1)
      IF ( DISCR ) RCONDA = DWORK(2)
C
      ISCL = 0
      IF ( LSCAL ) THEN
C
C        Scale the Hamiltonian or symplectic matrix.     
C
         IF( QNORM.GT.GNORM .AND. GNORM.GT.ZERO ) THEN
            CALL PDLASCL( 'G', QNORM, GNORM, N, N, S, IS+NP1-1, JS, 
     $                    DESCS, IERR ) 
            CALL PDLASCL( 'G', GNORM, QNORM, N, N, S, IS, JS+NP1-1, 
     $                    DESCS, IERR )
            ISCL = 1
         END IF
      END IF
C
C     Find the ordered Schur factorization of S = U*H*U'.
C
      IF( DEBUG ) WRITE(*,*) 'PDGEES in...'
      IF ( .NOT.DISCR ) THEN
         IF ( LSORT ) THEN
            CALL PDGEES( 'Vectors', 'Sorted', SB02MV, N2, S, IS, JS,
     $                   DESCS, NROT, WR, WI, U, IU, JU, DESCU, DWORK, 
     $                   LDWORK, IWORK, LIWORK, INFO )
         ELSE
            CALL PDGEES( 'Vectors', 'Sorted', SB02MR, N2, S, IS, JS,
     $                   DESCS, NROT, WR, WI, U, IU, JU, DESCU, DWORK, 
     $                   LDWORK, IWORK, LIWORK, INFO )
         END IF
      ELSE 
         IF ( LSORT ) THEN
            CALL PDGEES( 'Vectors', 'Sorted', SB02MW, N2, S, IS, JS,
     $                   DESCS, NROT, WR, WI, U, IU, JU, DESCU, DWORK, 
     $                   LDWORK, IWORK, LIWORK, INFO )
         ELSE
            CALL PDGEES( 'Vectors', 'Sorted', SB02MS, N2, S, IS, JS,
     $                   DESCS, NROT, WR, WI, U, IU, JU, DESCU, LDWORK, 
     $                   DWORK, IWORK, LIWORK, INFO )
         END IF
         IF ( LHINV ) THEN
            CALL DSWAP( N, WR, 1, WR(NP1), 1 )
            CALL DSWAP( N, WI, 1, WI(NP1), 1 )
         END IF
      END IF
      TIMES( 1 ) = DWORK( 2 )
      TIMES( 2 ) = DWORK( 3 )
      TIMES( 3 ) = DWORK( 4 )
C
      IF( PRINT ) THEN
         CALL PDLAPRNT( N2, N2, S, IS, JS, DESCS, 0, 0,
     $                 'S1', 6, DWORK )
         CALL PDLAPRNT( N2, N2, U, IU, JU, DESCU, 0, 0,
     $                 'U1', 6, DWORK )
      END IF  
      IF ( INFO.GT.N2 ) THEN
         INFO = 3
      ELSE IF ( INFO.GT.0 ) THEN
         INFO = 2
      ELSE IF ( NROT.NE.N ) THEN
         INFO = 4
      END IF
      IF( DEBUG ) WRITE(*,*) 'INFO from PDGEES =',INFO
      IF ( INFO.NE.0 )
     $     RETURN
C
      IF( DEBUG ) WRITE(*,*) 'WRKOPT'
      WRKOPT = MAX( WRKOPT, DWORK(1) )
C
C     Check if U(1,1) is singular.  Use the (2,1) block of S as a
C     workspace for storing the original U(1,1).
C
      UNORM = PDLANGE( '1-norm', N, N, U, IU, JU, DESCU, DWORK )
C
      IF( DEBUG ) WRITE(*,*) 'PDLAMVE 1'
      CALL PDLAMVE( 'Full', N, N, U, IU, JU, DESCU, S, IS+NP1-1, JS, 
     $              DESCS, DWORK )
      IF( DEBUG ) WRITE(*,*) 'PDGETRF 1'
      CALL PDGETRF( N, N, U, IU, JU, DESCU, IWORK, INFO )
      IF( PRINT ) THEN
         CALL PDLAPRNT( N2, N2, S, IS, JS, DESCS, 0, 0,
     $                 'S2', 6, DWORK )
         CALL PDLAPRNT( N2, N2, U, IU, JU, DESCU, 0, 0,
     $                 'U2', 6, DWORK )
      END IF  
C
      IF ( INFO.GT.0 ) THEN
C
C        Singular matrix.  Set INFO and RCOND for error return.
C
         INFO  = 5
         RCOND = ZERO
         GO TO 100
      END IF
C
C     Estimate the reciprocal condition of U(1,1).
C
      NIPIV = NUMROC( N, DESCS(MB_), MYROW, DESCS(RSRC_), NPROW ) +
     $                DESCS(MB_)
      IF( DEBUG ) WRITE(*,*) 'PDGECON 1'
      CALL PDGECON( '1-norm', N, U, IU, JU, DESCU, UNORM, RCOND,
     $              DWORK, LDWORK, IWORK(NIPIV+1), LIWORK-NIPIV, INFO )
C
      IF ( RCOND.LT.DLAMCH( 'Epsilon' ) ) THEN
C
C        Nearly singular matrix.  Set INFO for error return.
C
         INFO = 5
         RETURN
      END IF
C
C     Transpose U(2,1) in Q and compute the solution.
C
      IF( DEBUG ) WRITE(*,*) 'PDLAMVE 2'
      CALL PDTRAN( N, N, ONE, U, IU+NP1-1, JU, DESCU, ZERO, Q, IQ, JQ, 
     $             DESCQ )
C
      IF( DEBUG ) WRITE(*,*) 'PDGETRS 1',IS,JS,IQ,JQ,NP1
      CALL PDGETRS( 'Transpose', N, N, U, IU, JU, DESCU, IWORK, Q, 
     $              IQ, JQ, DESCQ, INFO )
      IF( PRINT ) THEN
         CALL PDLAPRNT( N, N, Q, IQ, JQ, DESCQ, 0, 0,
     $                 'QQ', 6, DWORK )
      END IF 
C
C     Copy back the original U(1,1) and set S(2,1) to zero.
C
      CALL PDLAMVE( 'Full', N, N, S, IS+NP1-1, JS, DESCS, U, IU, JU, 
     $              DESCU, DWORK )
      IF( DEBUG ) WRITE(*,*) 'PDLASET 1'
      CALL PDLASET( 'Full', N, N, ZERO, ZERO, S, IS+NP1-1, JS, DESCS )
      IF( PRINT ) THEN
         CALL PDLAPRNT( N2, N2, S, IS, JS, DESCS, 0, 0,
     $                 'S3', 6, DWORK )
         CALL PDLAPRNT( N2, N2, U, IU, JU, DESCU, 0, 0,
     $                 'U3', 6, DWORK )
      END IF  
C
C     Make sure the solution matrix X is symmetric.
C
      IF( DEBUG ) WRITE(*,*) 'Symmetric check'
      DO 80 I = 1, N - 1
         CALL PDAXPY( N-I, ONE, Q, IQ+I-1, JQ+I, DESCQ, DESCQ(M_), 
     $                Q, IQ+I, JQ+I-1, DESCQ, 1 )
         CALL PDSCAL( N-I, HALF, Q, IQ+I, JQ+I-1, DESCQ, 1 )
         CALL PDCOPY( N-I, Q, IQ+I, JQ+I-1, DESCQ, 1, Q, IQ+I-1, JQ+I, 
     $                DESCQ, DESCQ(M_) )
   80 CONTINUE
      IF( PRINT ) THEN
         CALL PDLAPRNT( N, N, Q, IQ, JQ, DESCQ, 0, 0,
     $                 'QQQ', 6, DWORK )
      END IF 
      IF( DEBUG ) WRITE(*,*) 'Finished Symmetric check'
C
      IF( LSCAL ) THEN
C
C        Undo scaling for the solution matrix.
C
         IF( ISCL.EQ.1 )
     $      CALL PDLASCL( 'G', GNORM, QNORM, N, N, Q, IQ, JQ, DESCQ, 
     $                    IERR )
      END IF
C
C     Set the optimal workspace, the scaling factor, and reciprocal
C     condition number (if any).
C
      DWORK( 1 ) = WRKOPT
      IWORK( 1 ) = IWRKOPT
  100 CONTINUE
      IF( ISCL.EQ.1 ) THEN
         DWORK(2) = QNORM / GNORM
      ELSE
         DWORK(2) = ONE
      END IF
      IF ( DISCR ) DWORK(3) = RCONDA
C
      RETURN
C *** Last line of PSB02MD ***
      END
