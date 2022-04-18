      SUBROUTINE PSB02MU( DICO, HINV, UPLO, N, A, IA, JA, DESCA, G, IG, 
     $                    JG, DESCG, Q, IQ, JQ, DESCQ, S, IS, JS, DESCS, 
     $                    DWORK, LDWORK, IWORK, LIWORK, INFO )
C   
C  -- ScaLAPACK/PSLICOT-style routine --
C     Preliminary version.
C     Dept. Computing Science and HPC2N, Univ. of Umeå, Sweden
C     December 7, 2007.
C
      IMPLICIT NONE
C
C     .. Scalar Arguments ..
      CHARACTER         DICO, HINV, UPLO
      INTEGER           INFO, IA, JA, IG, JG, IQ, JQ, IS, JS, LDWORK, 
     $                  LIWORK, N
C     .. Array Arguments ..
      INTEGER           DESCA(*), DESCG(*), DESCQ(*), DESCS(*), IWORK(*)
      DOUBLE PRECISION  A(*), DWORK(*), G(*), Q(*), S(*)
C
C     PURPOSE
C
C     To construct the 2n-by-2n Hamiltonian or symplectic matrix S
C     associated to the linear-quadratic optimization problem, used to
C     solve the continuous- or discrete-time algebraic Riccati equation,
C     respectively.
C
C     For a continuous-time problem, S is defined by
C
C             (  A  -G )
C         S = (        ),                                       (1)
C             ( -Q  -A')
C
C     and for a discrete-time problem by
C
C                 -1       -1
C             (  A        A  *G     )
C         S = (   -1           -1   ),                          (2)
C             ( QA     A' + Q*A  *G )
C
C     or
C
C                       -T         -T
C             (  A + G*A  *Q   -G*A   )
C         S = (      -T            -T ),                        (3)
C             (    -A  *Q         A   )
C
C     where A, G, and Q are N-by-N matrices, with G and Q symmetric.
C     Matrix A must be nonsingular in the discrete-time case.
C
C     ARGUMENTS
C
C     Mode Parameters
C
C     DICO    CHARACTER*1
C             Specifies the type of the system as follows:
C             = 'C':  Continuous-time system;
C             = 'D':  Discrete-time system.
C
C     HINV    CHARACTER*1
C             If DICO = 'D', specifies which of the matrices (2) or (3)
C             is constructed, as follows:
C             = 'D':  The matrix S in (2) is constructed;
C             = 'I':  The (inverse) matrix S in (3) is constructed.
C             HINV is not referenced if DICO = 'C'.
C
C     UPLO    CHARACTER*1
C             Specifies which triangle of the matrices G and Q is
C             stored, as follows:
C             = 'U':  Upper triangle is stored;
C             = 'L':  Lower triangle is stored.
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
C             The array descriptor for the distributed matrix Q
C
C     S       (local output) DOUBLE PRECISION pointer into the
C             local memory to an array, dimension (LLD_S,LOCq(2*N))
C             If INFO = 0 or INFO = 5, the this array contains the 
C             ordered real Schur form S of the Hamiltonian or 
C             symplectic matrix H. 
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
C     Workspace
C
C     DWORK   (local workspace) DOUBLE PRECISION array, dimension (LDWORK)
C             On exit, if INFO = 0, DWORK(1) returns the optimal value
C             of LDWORK; if DICO = 'D', DWORK(2) returns the reciprocal
C             condition number of the given matrix  A.
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
C     INFO    INTEGER
C             = 0:   successful exit;
C             < 0:   If the i-th argument is an array and the j-entry
C                    had an illegal value, then INFO = -(i*100+j), if the
C                    i-th argument is a scalar and had an illegal value,
C                    then INFO = -i.
C             = i:   if the leading i-by-i (1 <= i <= N) upper triangular
C                    submatrix of A is singular in discrete-time case;
C             = N+1: if matrix A is numerically singular in discrete-
C                    time case.
C
C     METHOD
C
C     For a continuous-time problem, the 2n-by-2n Hamiltonian matrix (1)
C     is constructed.
C     For a discrete-time problem, the 2n-by-2n symplectic matrix (2) or
C     (3) - the inverse of the matrix in (2) - is constructed.
C
C     NUMERICAL ASPECTS
C
C     The discrete-time case needs the inverse of the matrix A, hence
C     the routine should not be used when A is ill-conditioned.
C                               3
C     The algorithm requires 0(n ) floating point operations in the
C     discrete-time case.
C
C     PARALLEL EXECUTION RECOMMENDATIONS
C
C     Restrictions:
C      o All matrices should be be completely aligned, i.e.,
C        having the same MB, NB, RSRC, CSRC. 
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
      DOUBLE PRECISION  ZERO, ONE
      PARAMETER         ( ZERO = 0.0D0, ONE = 1.0D0 )
C     .. Local Scalars ..
      LOGICAL           DISCR, LHINV, LUPLO, LQUERY
      INTEGER           I, J, MAXWRK, N2, NJ, NP1, ICTXT, NPROW, NPCOL,
     $                  MYROW, MYCOL, WRKOPT, IWRKOPT, NIPIV, NB,
     $                  IPW1
      DOUBLE PRECISION  ANORM, RCOND
C     .. External Functions ..
      LOGICAL           LSAME
      INTEGER           NUMROC, ICEIL
      DOUBLE PRECISION  DLAMCH, PDLANGE
      EXTERNAL          DLAMCH, PDLANGE, LSAME, NUMROC, ICEIL
C     .. External Subroutines ..
      EXTERNAL          PDLACPY, PDGECON, PDGEMM, PDGETRF, PDGETRI, 
     $                  PDGETRS, PDLAMVE, PXERBLA
C     .. Intrinsic Functions ..
      INTRINSIC         MAX
C     .. Executable Statements ..
C
      IF( DEBUG ) WRITE(*,*) 'PSB02MU'
      INFO  = 0
      LQUERY = ( ( LDWORK.EQ.-1 ).OR.( LIWORK.EQ.-1 ) )
      ICTXT = DESCA( CTXT_ )
      NB = DESCA(MB_)
      CALL BLACS_GRIDINFO( ICTXT, NPROW, NPCOL, MYROW, MYCOL )
      IF( NPROW.EQ.-1 ) INFO = -(600+CTXT_)
      IF( INFO.EQ.0 ) THEN
         N2 = N + N
         DISCR = LSAME( DICO,  'D' )
         LUPLO = LSAME( UPLO,  'U' )
         IF( DISCR ) THEN
            LHINV = LSAME( HINV, 'D' )
         ELSE
            LHINV = .FALSE.
         END IF
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
         ELSE IF( N.LT.0 ) THEN
            INFO = -4
         END IF
C     
C        Test input matrices
C     
         IF( INFO.EQ.0 ) THEN
            CALL CHK1MAT( N, 4, N, 4, IA, JA, DESCA, 8, INFO )
            CALL CHK1MAT( N, 4, N, 4, IG, JG, DESCG, 12, INFO )
            CALL CHK1MAT( N, 4, N, 4, IQ, JQ, DESCQ, 16, INFO )
c            CALL CHK1MAT( N2, 4, N2, 4, IS, JS, DESCS, 20, INFO )
c            CALL PCHK2MAT( N, 6, N, 6, IA, JA, DESCA, 8, N2, 6, N2, 
c     $                     6, IS, JS, DESCS, 20, 0, IWORK, IWORK, INFO )
         END IF
C
C        Test workspace - TODO!
C
         IF( INFO.EQ.0 ) THEN
            WRKOPT = N
            IWRKOPT = N
            IF( LDWORK.LT.WRKOPT ) THEN
               INFO = -22
            ELSEIF( LIWORK.LT.IWRKOPT ) THEN
               INFO = -24
            END IF
         END IF
C
C        Error or workspace check return.
C
         IF( .NOT.LQUERY .AND. INFO.NE.0 ) THEN
            CALL PXERBLA( 'PSB02MU', -INFO )
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
         DWORK(1) = ONE
         IF ( DISCR ) DWORK(2) = ONE
         RETURN
      END IF
C
C     Compute upper bound on buffer size used in PDLAMVE
C
      IPW1 = ICEIL(ICEIL(N,NB),NPROW)*NB *
     $     ICEIL(ICEIL(N,NB),NPCOL)*NB
C
C     The code tries to exploit data locality as much as possible.
C
      IF ( .NOT.LHINV ) THEN
         CALL PDLACPY( 'Full', N, N, A, IA, JA, DESCA, S, IS, JS, 
     $                 DESCS )
C
C        Construct Hamiltonian matrix in the continuous-time case, or
C        prepare symplectic matrix in (3) in the discrete-time case:
C
C        Construct full Q in S(N+1:2*N,1:N) and change the sign, and
C        construct full G in S(1:N,N+1:2*N) and change the sign.
C
         IF( LUPLO ) THEN
            IF( DEBUG ) WRITE(*,*) MYROW,MYCOL,'PDLAMVE 1'
            CALL PDLAMVE( 'Upper', N, N, Q, IQ, JQ, DESCQ, S, IS+N, JS, 
     $                    DESCS, DWORK )
            CALL PDTRAN( N, N, ONE, Q, IQ, JQ, DESCQ, ZERO, DWORK,
     $                   IQ, JQ, DESCQ )
            IF( DEBUG ) WRITE(*,*) MYROW,MYCOL,'PDLAMVE 2'
            CALL PDLAMVE( 'Lower', N-1, N-1, DWORK, IQ+1, JQ, DESCQ, S, 
     $                    IS+N+1, JS, DESCS, DWORK(IPW1) )
            DO 20 J = 1, N
               CALL PDSCAL( N, -ONE, S, IS+N, JS+J-1, DESCS, 1 )
 20         CONTINUE 
C
            IF( DEBUG ) WRITE(*,*) MYROW,MYCOL,'PDLAMVE 3'
            CALL PDLAMVE( 'Upper', N, N, G, IG, JG, DESCG, S, IS, JS+N, 
     $                    DESCS, DWORK )
            CALL PDTRAN( N, N, ONE, G, IG, JG, DESCG, ZERO, DWORK,
     $                   IG, JG, DESCG )
            IF( DEBUG ) WRITE(*,*) MYROW,MYCOL,'PDLAMVE 4'
            CALL PDLAMVE( 'Lower', N-1, N-1, DWORK, IG+1, JG, DESCG, S, 
     $                    IS+1, JS+N, DESCS, DWORK(IPW1) )
            DO 30 J = 1, N
               CALL PDSCAL( N, -ONE, S, IS, JS+N+J-1, DESCS, 1 )
 30         CONTINUE  
         ELSE
            IF( DEBUG ) WRITE(*,*) MYROW,MYCOL,'PDLAMVE 5'
            CALL PDLAMVE( 'Lower', N, N, Q, IQ, JQ, DESCQ, S, IS+N, JS, 
     $                    DESCS, DWORK )
            CALL PDTRAN( N, N, ONE, Q, IQ, JQ, DESCQ, ZERO, DWORK,
     $                   IQ, JQ, DESCQ )
           IF( DEBUG ) WRITE(*,*) MYROW,MYCOL,'PDLAMVE 6' 
            CALL PDLAMVE( 'Upper', N-1, N-1, DWORK, IQ, JQ+1, DESCQ, S, 
     $                    IS+N, JS+1, DESCS, DWORK(IPW1) )
            DO 40 J = 1, N
               CALL PDSCAL( N, -ONE, S, IS+N, JS+J-1, DESCS, 1 )
 40         CONTINUE 
C
            IF( DEBUG ) WRITE(*,*) MYROW,MYCOL,'PDLAMVE 7'
            CALL PDLAMVE( 'Lower', N, N, G, IG, JG, DESCG, S, IS, JS+N, 
     $                    DESCS, DWORK )
            CALL PDTRAN( N, N, ONE, G, IG, JG, DESCG, ZERO, DWORK,
     $                   IG, JG, DESCG )
            IF( DEBUG ) WRITE(*,*) MYROW,MYCOL,'PDLAMVE 8'
            CALL PDLAMVE( 'Upper', N-1, N-1, DWORK, IG, JG+1, DESCG, S, 
     $                    IS, JS+N+1, DESCS, DWORK(IPW1) )
            DO 50 J = 1, N
               CALL PDSCAL( N, -ONE, S, IS, JS+N+J-1, DESCS, 1 )
 50         CONTINUE  
         END IF
C
         IF ( .NOT.DISCR ) THEN
            CALL PDTRAN( N, N, -ONE, A, IA, JA, DESCA, ZERO, S, IS+N, 
     $                   JS+N, DESCS )
C
            DWORK(1) = ONE
         END IF
      END IF
C
      IF ( DISCR ) THEN
C
C        Construct the symplectic matrix (2) or (3) in the discrete-time
C        case.
C
         IF ( LHINV ) THEN
C
C           Put  A'  in  S(N+1:2*N,N+1:2*N).
C
            CALL PDTRAN( N, N, ONE, A, IA, JA, DESCA, ZERO, S, IS+N, 
     $                  JS+N, DESCS )
C
         END IF
C
C        Compute the norm of the matrix A.
C
         ANORM = PDLANGE( '1-norm', N, N, A, IA, JA, DESCA, DWORK )
C
C        Compute the LU factorization of A.
C
         CALL PDGETRF( N, N, A, IA, JA, DESCA, IWORK, INFO )
C
C        Return if INFO is non-zero.
C
         IF( INFO.GT.0 ) THEN
            DWORK(2) = ZERO
            RETURN
         END IF
C
C        Compute the reciprocal of the condition number of A.
C
         NIPIV = NUMROC( N, DESCA(MB_), MYROW, DESCA(RSRC_), NPROW ) +
     $           DESCA(MB_)
         CALL PDGECON( '1-norm', N, A, IA, JA, DESCA, ANORM, RCOND, 
     $                DWORK, IWORK(NIPIV+1), INFO )
C
C        Return if the matrix is singular to working precision.
C
         IF( RCOND.LT.DLAMCH( 'Epsilon' ) ) THEN
            INFO = N + 1
            DWORK(2) = RCOND
            RETURN
         END IF
C
         IF ( LHINV ) THEN
C
C           Compute S in (2).
C
C           Construct full Q in S(N+1:2*N,1:N).
C
            IF( LUPLO ) THEN
               IF( DEBUG ) WRITE(*,*) MYROW,MYCOL,'PDLAMVE 9'
               CALL PDLAMVE( 'Upper', N, N, Q, IQ, JQ, DESCQ, S, IS+N, 
     $                       JS, DESCS, DWORK )
               CALL PDTRAN( N, N, ONE, Q, IQ, JQ, DESCQ, ZERO, DWORK,
     $                      IQ, JQ, DESCQ )
               IF( DEBUG ) WRITE(*,*) MYROW,MYCOL,'PDLAMVE 10'
               CALL PDLAMVE( 'Lower', N-1, N-1, DWORK, IQ+1, JQ, DESCQ, 
     $                       S, IS+N+1, JS, DESCS, DWORK(IPW1) )
            ELSE
               IF( DEBUG ) WRITE(*,*) MYROW,MYCOL,'PDLAMVE 11'
               CALL PDLAMVE( 'Lower', N, N, Q, IQ, JQ, DESCQ, S, IS+N, 
     $                       JS, DESCS, DWORK )
               CALL PDTRAN( N, N, ONE, Q, IQ, JQ, DESCQ, ZERO, DWORK,
     $                      IQ, JQ, DESCQ )
               IF( DEBUG ) WRITE(*,*) MYROW,MYCOL,'PDLAMVE 12'
               CALL PDLAMVE( 'Upper', N-1, N-1, DWORK, IQ, JQ+1, DESCQ, 
     $                       S, IS+N, JS+1, DESCS, DWORK(IPW1) )
            END IF
C
C           Compute the solution matrix  X  of the system  X*A = Q  by
C                                                                    -1
C           solving  A'*X' = Q and transposing the result to get  Q*A  .
C
            CALL PDGETRS( 'Transpose', N, N, A, IA, JA, DESCA, IWORK, 
     $                    S, IS+N, JS, DESCS, INFO )
            CALL PDTRAN( N, N, ONE, S, IS+N, JS, DESCS, ZERO, DWORK,
     $                   IQ, JQ, DESCQ )
            IF( DEBUG ) WRITE(*,*) MYROW,MYCOL,'PDLAMVE 13'
            CALL PDLAMVE( 'All', N, N, DWORK, IQ, JQ, DESCQ, S, IS+N, 
     $                    JS, DESCS, DWORK(IPW1) )
C
C           Construct full G in S(1:N,N+1:2*N).
C
            IF( LUPLO ) THEN
               IF( DEBUG ) WRITE(*,*) MYROW,MYCOL,'PDLAMVE 14'
               CALL PDLAMVE( 'Upper', N, N, G, IG, JG, DESCG, S, IS, 
     $                       JS+N, DESCS, DWORK )
               CALL PDTRAN( N, N, ONE, G, IG, JG, DESCG, ZERO, DWORK,
     $                      IG, JG, DESCG )
               IF( DEBUG ) WRITE(*,*) MYROW,MYCOL,'PDLAMVE 15'
               CALL PDLAMVE( 'Lower', N-1, N-1, DWORK, IG+1, JG, DESCG, 
     $                       S, IS+1, JS+N, DESCS, DWORK(IPW1) )  
            ELSE
               IF( DEBUG ) WRITE(*,*) MYROW,MYCOL,'PDLAMVE 16'
               CALL PDLAMVE( 'Lower', N, N, G, IG, JG, DESCG, S, IS, 
     $                       JS+N, DESCS, DWORK )
               CALL PDTRAN( N, N, ONE, G, IG, JG, DESCG, ZERO, DWORK,
     $                      IG, JG, DESCG )
               IF( DEBUG ) WRITE(*,*) MYROW,MYCOL,'PDLAMVE 17'
               CALL PDLAMVE( 'Upper', N-1, N-1, DWORK, IG, JG+1, DESCG, 
     $                       S, IS, JS+N+1, DESCS, DWORK(IPW1) )
            END IF
C                            -1
C           Compute  A' + Q*A  *G  in  S(N+1:2N,N+1:2N).
C
            CALL PDGEMM( 'No transpose', 'No transpose', N, N, N, ONE,
     $                   S, IS+N, JS, DESCS, S, IS, JS+N, DESCS, ONE, 
     $                   S, IS+N, JS+N, DESCS )
C
C           Compute the solution matrix  Y  of the system  A*Y = G.
C
            CALL PDGETRS( 'No transpose', N, N, A, IA, JA, DESCA, IWORK, 
     $                    S, IS, JS+N, DESCS, INFO )
C
C           Compute the inverse of  A  in situ.
C
            CALL PDGETRI( N, A, IA, JA, DESCA, IWORK, DWORK, LDWORK, 
     $                    INFO )
C                  -1
C           Copy  A    in  S(1:N,1:N).
C
            IF( DEBUG ) WRITE(*,*) MYROW,MYCOL,'PDLAMVE 18'
            CALL PDLAMVE( 'Full', N, N, A, IA, JA, DESCA, S, IS, JS,
     $                    DESCS, DWORK )
         ELSE
C
C           Compute S in (3) using the already prepared part.
C
C           Compute the solution matrix  X'  of the system  A*X' = -G
C                                                       -T
C           and transpose the result to obtain  X = -G*A  .
C
            CALL PDGETRS( 'No transpose', N, N, A, IA, JA, DESCA, IWORK, 
     $                    S, IS, JS+N, DESCS, INFO )
            CALL PDTRAN( N, N, ONE, S, IS, JS+N, DESCS, ZERO, DWORK,
     $                   IA, JA, DESCA )
            IF( DEBUG ) WRITE(*,*) MYROW,MYCOL,'PDLAMVE 19'
            CALL PDLAMVE( 'All', N, N, DWORK, IA, JA, DESCA, S, IS, 
     $                    JS+N, DESCS, DWORK(IPW1) )
C
C                           -T
C           Compute  A + G*A  *Q  in  S(1:N,1:N).
C
            CALL PDGEMM( 'No transpose', 'No transpose', N, N, N, ONE,
     $                   S, IS, JS+N, DESCS, S, IS+N, JS, DESCS, ONE, S, 
     $                   IS, JS, DESCS )
C
C           Compute the solution matrix  Y  of the system  A'*Y = -Q.
C
            CALL PDGETRS( 'Transpose', N, N, A, IA, JA, DESCA, IWORK, 
     $                    S, IS+N, JS, DESCS, INFO )
C
C           Compute the inverse of  A  in situ.
C
            CALL PDGETRI( N, A, IA, JA, DESCA, IWORK, DWORK, LDWORK, 
     $                    INFO )
C
C                  -T
C           Copy  A    in  S(N+1:2N,N+1:2N).
C
            CALL PDTRAN( N, N, ONE, A, IA, JA, DESCA, ZERO, S, IS+N,
     $                   JS+N, DESCS )
C
         END IF
         DWORK(1) = WRKOPT
         IWORK(1) = IWRKOPT
         DWORK(2) = RCOND
      END IF
C
      RETURN
C
C *** Last line of SB02MU ***
      END
