      PROGRAM TEST
*
* ======================================================================
*
*  Test routine of SRRIT - a subroutine to calculate the dominant
*  invariant subspace of a nonsymmetric matrix. The SRRIT user guide is
*  documented in
*
*  Z. Bai and G. W. Stewart, SRRIT - A FORTRAN Subroutine to
*    calculate the dominant invariant subspace of a nonsymmetric
*    matrix, submitted to ACM TOMS.
*
*  With the file ``input'' as input, this test driver will exercise 
*  major features of SRRIT. The output is going to be in the file 
*  named ``output''.
*
*  The current array sizes for testing SRRIT are set to handle the order
*  of a test matrix less than 1000, and the dimension of subspace less
*  than 100. For running a larger test or a larger subspace, user needs 
*  to change the dimension parameters LDA and MAXSPC.
*  In addition, see the comments of SCHECK for the limit of the test 
*  problem sizes for the testing subroutine SLAQR3. 
*
*  ---------------- 
*
*  Copyright (C) 1994, Zhaojun Bai and G. W. Stewart
*  All rights reserved.  Use at your own risk.
*  Please send comments to bai@ms.uky.edu or stewart@cs.umd.edu
*
* ======================================================================
*
*     .. Parameters ..
      INTEGER            LDA, LDQ, MAXSPC, LDT
      PARAMETER          ( LDA = 1000, LDQ = LDA, MAXSPC = 100,
     $                   LDT = MAXSPC )
      INTEGER            LWORK
      PARAMETER          ( LWORK = MAXSPC*MAXSPC+5*MAXSPC )
      REAL               ZERO, HALF, ONE
      PARAMETER          ( ZERO = 0.0E+0, HALF = 0.5E+0, ONE = 1.0E+0 )
      INTEGER            NIN, NOUT
      PARAMETER          ( NIN = 5, NOUT = 6 )
*     ..
*     .. Local Scalars ..
      CHARACTER*3        PTYPE
      INTEGER            I, INFO, ISTART, J, K, M, MAXIT, N, NV
      REAL               QNORM, TOL
*     ..
*     .. Local Arrays ..
      INTEGER            ITRSD( MAXSPC ), IWORK( LWORK ), ISEED( 4 ), 
     $                   NMAX( MAXSPC ) 
      REAL               AQ( LDA, MAXSPC ), DUMMY( 1 ),
     $                   Q( LDQ, MAXSPC ), RSD( MAXSPC ), T( LDT, LDT ),
     $                   WI( LDT ), WORK( LWORK ), WR( LDT )
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      REAL               SLANGE, SLARAN
      EXTERNAL           LSAME, SLANGE, SLARAN
*     ..
*     .. External Subroutines ..
      EXTERNAL           SGEMM, SLASET, SRRIT
      EXTERNAL           ATQODE, ATQPDE, ATQRWK, ATQBWM, ATQ4X4, ATQF44
*
*     .. Intrinsic functions ..
      INTRINSIC          SIGN
*     ..
*     .. Executable Statements ..
*
      OPEN( NOUT, FILE = 'soutput' )
*
   10 CONTINUE
*
*     Input the Type of Test Problem
*
      READ( NIN, FMT = * )PTYPE
*
      IF( .NOT.LSAME( PTYPE, 'END' ) )
     $   WRITE( NOUT, FMT = 9999 )PTYPE
*
      IF( LSAME( PTYPE, 'END' ) ) THEN
         WRITE( NOUT, FMT = 9998 )
         GO TO 80
      END IF
*
*     Input Parameter K, which decides the size of test problem
*
      READ( NIN, FMT = * )K
*
      IF( LSAME( PTYPE, 'RWK' ) .OR. LSAME( PTYPE, 'LRW' ) .OR. 
     $    LSAME( PTYPE, 'DWK' ) ) THEN
         N = ( K+1 )*( K+2 ) / 2
      ELSE IF( LSAME( PTYPE, 'ODE' ) ) THEN
         N = K + 1
      ELSE IF( LSAME( PTYPE, 'PDE' ) ) THEN
         N = K*K
      ELSE IF( LSAME( PTYPE, 'BWM' ) ) THEN
         N = K
      ELSE IF( LSAME( PTYPE, 'M44' ) .OR. LSAME( PTYPE, 'N44' ) ) THEN 
         N = K 
      ELSE IF( LSAME( PTYPE, 'F44' ) ) THEN 
         N = K 
      ELSE IF( LSAME( PTYPE, 'SUB' ) ) THEN 
         GO TO 90 
      END IF
*
      WRITE( NOUT, FMT = 9996 )N
*
      IF( N.GT.LDA ) THEN
         WRITE( NOUT, FMT = 9995 )
         GO TO 80
      END IF
*
*     Input NV, the number of the desired dominant eigenvalues.
*           M,  the subspace size  
*     Note that it is required that M >= NV+2.
*
      READ( NIN, FMT = * )NV, M
*
      WRITE( NOUT, FMT = 1111 )NV, M
 1111 FORMAT( ' NV and M = ', 2I5) 
*
*     Set stopping criterion TOL and maximal number of iterations. 
*     They may be altered if desired.
*
      TOL = 1E-5
      MAXIT = 10*N
*
*     Call main SRRIT routine
*
      IF( LSAME( PTYPE, 'RWK' ) ) THEN
         INFO = 1 
         ISTART = -1
         CALL SRRIT( N, NV, M, MAXIT, ISTART, Q, LDQ, AQ, LDA, T, LDT,
     $               WR, WI, RSD, ITRSD, IWORK, WORK, LWORK, INFO, TOL,
     $               ATQRWK )
      ELSE IF( LSAME( PTYPE, 'LRW' ) ) THEN
         INFO = 0
         ISTART = 0  
*
*        generate initial vectors (istart = 0)
* 
         ISEED( 1 ) = 3192 
         ISEED( 2 ) = 1221 
         ISEED( 3 ) = 1322 
         ISEED( 4 ) = 8097 
         DO 30 J = 1, M
            DO 20 I = 1, N 
               Q( I,J ) = SIGN( ONE, SLARAN( ISEED ) - HALF )  
   20       CONTINUE
   30    CONTINUE
*
         CALL SRRIT( N, NV, M, MAXIT, ISTART, Q, LDQ, AQ, LDA, T, LDT,
     $               WR, WI, RSD, ITRSD, IWORK, WORK, LWORK, INFO, TOL,
     $               ATQRWK )
      ELSE IF( LSAME( PTYPE, 'DRW' ) ) THEN
         INFO = 0
         ISTART = 1 
*
*        generate and orthogonalize initial vectors (istart = 1)
* 
         ISEED( 1 ) = 3192
         ISEED( 2 ) = 1221
         ISEED( 3 ) = 1322
         ISEED( 4 ) = 8097
         DO 50 J = 1, M
            DO 40 I = 1, N
               Q( I,J ) = SLARAN( ISEED )
   40       CONTINUE
   50    CONTINUE
         CALL ORTH( N, 1, M, Q, LDQ, INFO )
         IF( INFO.NE.0 )THEN
            WRITE( NOUT, 1112 )INFO
 1112       FORMAT( 'Given starting vectors are linearly dependent' )
            GO TO 80
         END IF 
*
         CALL SRRIT( N, NV, M, MAXIT, ISTART, Q, LDQ, AQ, LDA, T, LDT,
     $               WR, WI, RSD, ITRSD, IWORK, WORK, LWORK, INFO, TOL,
     $               ATQRWK )
      ELSE IF( LSAME( PTYPE, 'ODE' ) ) THEN
         INFO = 1 
         ISTART = -1
         CALL SRRIT( N, NV, M, MAXIT, ISTART, Q, LDQ, AQ, LDA, T, LDT,
     $               WR, WI, RSD, ITRSD, IWORK, WORK, LWORK, INFO, TOL,
     $               ATQODE )
      ELSE IF( LSAME( PTYPE, 'PDE' ) ) THEN
         INFO = 0 
         ISTART = -1
         CALL SRRIT( N, NV, M, MAXIT, ISTART, Q, LDQ, AQ, LDA, T, LDT,
     $               WR, WI, RSD, ITRSD, IWORK, WORK, LWORK, INFO, TOL,
     $               ATQPDE )
      ELSE IF( LSAME( PTYPE, 'BWM' ) ) THEN
         INFO = 0
         ISTART = -1
         CALL SRRIT( N, NV, M, MAXIT, ISTART, Q, LDQ, AQ, LDA, T, LDT,
     $               WR, WI, RSD, ITRSD, IWORK, WORK, LWORK, INFO, TOL,
     $               ATQBWM )
      ELSE IF( LSAME( PTYPE, 'M44' ) ) THEN
         INFO = 0
         ISTART = -1
         CALL SRRIT( N, NV, M, MAXIT, ISTART, Q, LDQ, AQ, LDA, T, LDT,
     $               WR, WI, RSD, ITRSD, IWORK, WORK, LWORK, INFO, TOL,
     $               ATQ4X4 )
      ELSE IF( LSAME( PTYPE, 'N44' ) ) THEN
         INFO = 0
         ISTART = 0 
         DO 55 J = 1, M
            DO 45 I = 1, N
               Q( I,J ) = ONE
   45       CONTINUE
   55    CONTINUE
         CALL SRRIT( N, NV, M, MAXIT, ISTART, Q, LDQ, AQ, LDA, T, LDT,
     $               WR, WI, RSD, ITRSD, IWORK, WORK, LWORK, INFO, TOL,
     $               ATQ4X4 )
      ELSE IF( LSAME( PTYPE, 'F44' ) ) THEN
         INFO = 0
         ISTART = 1
         DO 75 J = 1, M
            DO 65 I = 1, N
               Q( I,J ) = ZERO 
   65       CONTINUE
            Q( J, J ) = ONE
   75    CONTINUE
         CALL SRRIT( N, NV, M, MAXIT, ISTART, Q, LDQ, AQ, LDA, T, LDT,
     $               WR, WI, RSD, ITRSD, IWORK, WORK, LWORK, INFO, TOL,
     $               ATQF44 )
      END IF
*
      IF( INFO.EQ.1 )THEN
         WRITE( NOUT, FMT = 1113 )INFO 
 1113    FORMAT( ' ORTH subroutine fails, INFO = ', I2 )
      ELSE IF( INFO.EQ.2 )THEN
         WRITE( NOUT, FMT = 1114 )INFO 
 1114    FORMAT( ' SRRSTP subroutine fails, INFO = ', I2 )
      ELSE IF( INFO.EQ.3 )THEN
         WRITE( NOUT, FMT = 1115 )INFO 
 1115    FORMAT( ' Reduced matrix T is singular, INFO = ', I2 )
      ELSE IF( INFO.EQ.4 )THEN
         WRITE( NOUT, FMT = 1116 )INFO 
 1116    FORMAT( ' SRRIT fails to converge after MAXIT, INF0 = ', I2) 
         WRITE( NOUT, FMT = 1118 )
 1118    FORMAT( ' 2-norms of Residual Vectors ' ) 
         WRITE( NOUT, FMT = 9989 )( RSD( I ), I = 1, M )
         WRITE( NOUT, FMT = 1119 )
 1119    FORMAT( ' Iteration numbers where the res. are computed ' )
         WRITE( NOUT, FMT = 9985 )( ITRSD( I ), I = 1, M )
         DO 61 K = 1, M 
            WRITE( NOUT, FMT = 9991 )WR( K ), WI( K )
   61    CONTINUE
      ELSE 
         WRITE( NOUT, FMT = 9993 )INFO
         WRITE( NOUT, FMT = 1117 )MAXIT
 1117    FORMAT( ' Total number of SRRIT iterations = ', I5 ) 
      END IF 
*
      IF( INFO.EQ.0 .AND. NV.GT.0 )THEN 
         WRITE( NOUT, FMT = 9992 )
         DO 60 K = 1, NV
            WRITE( NOUT, FMT = 9991 )WR( K ), WI( K )
   60    CONTINUE
*
*        Test residual: R = A*Q - Q*T
*
         CALL SGEMM( 'No transpose', 'No transpose', N, NV, NV, -ONE, 
     $               Q, LDQ, T, LDT, ONE, AQ, LDA )
*
         DO 70 I = 1, NV
            WORK( I ) = SLANGE( 'MAX', N, 1, AQ( 1, I ), 1, DUMMY )
   70    CONTINUE
*
         WRITE( NOUT, FMT = 9990 )
         WRITE( NOUT, FMT = 9989 )( WORK( I ), I = 1, NV )
*
         CALL SLASET( 'Full', NV, NV, ZERO, ONE, WORK, NV )
         CALL SGEMM( 'Transpose', 'No transpose', NV, NV, N, -ONE, Q, 
     $               LDQ, Q, LDQ, ONE, WORK, NV )
         QNORM = SLANGE( 'Max', NV, NV, WORK, NV, DUMMY )
*
         WRITE( NOUT, FMT = 9988 )QNORM
*
      END IF
*
      GO TO 10
*
 9999 FORMAT( / 10X, ' ...... Type of Test problem: ', A3 , ' ......' )
 9998 FORMAT( / ' ..... END OF SRRIT TEST .....' )
 9996 FORMAT( / ' The test matrix size = ', I5 )
 9995 FORMAT( ' N is larger than array leading dimension LDA' )
 9993 FORMAT( / ' The INFO from SRRIT = ', I2 )
 9992 FORMAT( ' .. Converged Eigenvalues ..' )
 9991 FORMAT( ' WR = ', 1P,E15.6, ' WI = ', 1P,E15.6 )
 9990 FORMAT( ' Final residual norms ',
     $      'of (A*Q - Q*T) of converged RR pairs = ' )
 9989 FORMAT( 1P,6E10.3 )
 9985 FORMAT( 6I10 )
 9988 FORMAT( ' norm(Q*Q - I)   = ', 1P,E15.6 )
*
*     Test subroutines SLAQR3
*
   90 CONTINUE 
*
      READ( NIN, FMT = * )( NMAX( I ), I = 1, K ) 
*
      CALL SCHECK( K, NMAX, INFO )
      IF( INFO .EQ. 0 )THEN
         WRITE( NOUT, FMT = 9987 ) 
 9987    FORMAT( ' Passed all SLAQR3 tests ' ) 
      ELSE 
         WRITE( NOUT, FMT = 9986 )INFO
 9986    FORMAT( ' Test routine SCHECK INFO = ', I5 ) 
      END IF
*
      GO TO 10 
*
   80 CONTINUE
*
      CLOSE ( NOUT )
*
      END
*
****** begin of scheck.f ********************************************
*
      SUBROUTINE SCHECK( MAX, NMAX, INFO ) 
*
*     .. Scalar Arguments .. 
      INTEGER         MAX, INFO 
*
*     .. Arrary Arguments .. 
      INTEGER         NMAX( * ) 
*
*  Purpose
*  =======
*
*  SCHECK tests the subroutines SLAQR3. 
*  On normal return, INFO is to be zero, otherwise, if 
*
*   INFO = 1, SLAQR3 has illegal parameter in the calling sequence
*             or it fails to compute the complete Schur decomposition.    
*        = 2, the backward errors of the Schur decomposition fails 
*             test. 
*        = 3, the order of the computed eigenvalues has not been 
*             sorted correctly with option WANTT = .TRUE. and
*             WANTZ = .TRUE. in SLAQR3
*        = 4, the order of the computed eigenvalues has not been 
*             sorted correctly with option WANTT = .FALSE. and
*             WANTZ = .FALSE. in SLAQR3
*
*  Note: the internal parameter LDA is set so that the order of any 
*  test matrix is smaller than 50. 
*
*  ===================================================================
*
*     .. Parameters .. 
      INTEGER         LDA, LDB, LDH, LDU, LWORK
      PARAMETER       ( LDA = 50, LDB = LDA, LDH = LDA, LDU = LDA,
     $                  LWORK = 2*LDA*LDA )
*
      REAL            ZERO, TOL
      PARAMETER       ( ZERO = 0.0E+0, TOL = 2.0E+01 )
*
*     .. Local Scalars ..
      INTEGER         I, J, K, N, IERR
*
*     .. Local Arrays ..
      INTEGER         ISEED( 4 )
      REAL            A( LDA,LDA ), B( LDB,LDB ), U( LDU,LDU ),
     $                H( LDH, LDH ), WR( LDA ), WI(LDA), 
     $                WORK( LWORK ), RESULT( 2 )
*
*     .. External functions .. 
      REAL            SLAPY2, SLARAN
      EXTERNAL        SLAPY2, SLARAN
*
*     .. External Subroutines ..
      EXTERNAL        SLACPY, SGEHRD, SORGHR, SLAQR3, SGET21
*
      INFO = 0  
      ISEED( 1 ) = 8321
      ISEED( 2 ) = 9735
      ISEED( 3 ) = 8673
      ISEED( 4 ) = 2437
*
      DO 10 K = 1, MAX
         N = NMAX( K ) 
*
*        Generate the random test matrices 
*
         DO 30 J = 1,N
            DO 20  I = 1, N
               A( I,J ) = SLARAN( ISEED )
  20        CONTINUE 
  30     CONTINUE
*
*        Save a copy of the matrix
*
         CALL SLACPY( 'FULL', N, N, A, LDA, B, LDB )
*
*        Hessenberg reduction
*
         CALL SGEHRD( N, 1, N, A, LDA, WR, WORK, LWORK, IERR )
         DO 50 J = 1, N - 1
            DO 40 I = J + 2, N
               U( I, J ) = A( I, J )
               A( I, J ) = ZERO
   40       CONTINUE
   50    CONTINUE
         CALL SORGHR( N, 1, N, U, LDU, WR, WORK, LWORK, IERR )
*
         CALL SLACPY( 'FULL', N, N, A, LDA, H, LDH )
*
*        Compute the Schur decomposition of upper Hessenberg matrix
*
         CALL SLAQR3( .TRUE., .TRUE., N, 1, N, A, LDA, WR, WI, 1, N, 
     $                U, LDU, WORK, IERR )
         IF( IERR .NE. 0 )THEN
             INFO =  1
             RETURN
         END IF 
*
         CALL SGET21( N, B, LDB, A, LDA, U, LDU, WORK, RESULT )
         IF( RESULT( 1 ) .GT. TOL .OR. RESULT( 2 ) .GT. TOL )THEN
             INFO  = 2  
             RETURN
         END IF  
*
*        Check the ordering 
*
         DO 60 I = 1,N
            WORK( I ) = SLAPY2( WR( I ), WI( I ) )
  60     CONTINUE
         DO 80 I = N-1,-1,1
            DO 70 J = I+1,-1,1
               IF( WORK(I) .LT. WORK(J) )THEN
                  INFO = 3 
                  RETURN 
               END IF 
  70        CONTINUE          
  80     CONTINUE 
*
*        Compute the eigenvalues only. 
*
         CALL SLAQR3( .FALSE., .FALSE., N, 1, N, H, LDH, WR, WI, 
     $                1, N, U, LDU, WORK, IERR )
*
*        Check the ordering
*
         DO 90 I = 1,N
            WORK( I ) = SLAPY2( WR( I ), WI( I ) )
  90     CONTINUE
         DO 110 I = N-1,-1,1
            DO 100 J = I+1,-1,1
               IF( WORK(I) .LT. WORK(J) )THEN
                  INFO = 4
                  RETURN
               END IF
 100        CONTINUE
 110     CONTINUE

  10  CONTINUE 
*
      RETURN
      END  
*
********* begin of sget21.f ***************************************** 
*
      SUBROUTINE SGET21( N, A, LDA, B, LDB, U, LDU, WORK, RESULT )
*     ..
*     .. Scalar Arguments ..
      INTEGER            LDA, LDB, LDU, N
*     ..
*     .. Array Arguments ..
      REAL               A( LDA, * ), B( LDB, * ), RESULT( * ),
     $                   U( LDU, * ), WORK( * )
*
*  Purpose
*  =======
*
*  SGET21  generally checks a decomposition of the form
*
*               A = U B U'
*
*  where ' means transpose and U is orthogonal. Specifically, 
*
*               RESULT(1) = | A - U B U' | / ( |A| n ulp )
*  and 
*               RESULT(2) = | I - UU' | / ( n ulp )
*
*  Arguments
*  ==========
*
*  N      - INTEGER
*           The size of the matrix.  If it is zero, SGET21 does nothing.
*           It must be at least zero.
*           Not modified.
*
*  A      - REAL array of dimension ( LDA , N )
*           The original (unfactored) matrix.
*           Not modified.
*
*  LDA    - INTEGER
*           The leading dimension of A.  It must be at least 1
*           and at least N.
*           Not modified.
*
*  B      - REAL array of dimension ( LDB , N )
*           The factored matrix.
*           Not modified.
*
*  LDB    - INTEGER
*           The leading dimension of B.  It must be at least 1
*           and at least N.
*           Not modified.
*
*  U      - REAL array of dimension ( LDU, N ).
*           The orthogonal matrix in the decomposition.  
*           Not modified.
*
*  LDU    - INTEGER
*           The leading dimension of U.  LDU must be at least N and
*           at least 1.
*           Not modified.
*
*  WORK   - REAL array of dimension ( 2*N**2 )
*           Workspace.
*           Modified.
*
*  RESULT - REAL array of dimension ( 2 )
*           The values computed by the two tests described above.  The
*           values are currently limited to 1/ulp, to avoid overflow.
*           Errors are flagged by RESULT(1)=10/ulp.
*
*-----------------------------------------------------------------------
*
*     .. Parameters ..
*
      REAL               ZERO, ONE
      PARAMETER          ( ZERO = 0.0E+0, ONE = 1.0E+0 )
*     ..
*     .. Local Scalars ..
*
      INTEGER            JDIAG
      REAL               ANORM, ULP, UNFL, WNORM
*     ..
*     .. External Functions ..
*
      REAL               SLAMCH, SLANGE
      EXTERNAL           SLAMCH, SLANGE
*     ..
*     .. External Subroutines ..
*
      EXTERNAL           SGEMM, SLACPY
*     ..
*     .. Intrinsic Functions ..
*
      INTRINSIC          MAX, MIN, REAL
*     ..
*     .. Executable Statements ..
*
      RESULT( 1 ) = ZERO
      RESULT( 2 ) = ZERO
      IF( N.LE.0 )
     $   RETURN
*
*     Constants
*
      UNFL = SLAMCH( 'Safe minimum' )
      ULP = SLAMCH( 'Epsilon' )*SLAMCH( 'Base' )
*
*     Do Test 1
*
*     Norm of A:
*
      ANORM = MAX( SLANGE( '1', N, N, A, LDA, WORK ), UNFL )
*
*     Norm of A - UBU'
*
      CALL SLACPY( ' ', N, N, A, LDA, WORK, N )
      CALL SGEMM( 'N', 'N', N, N, N, ONE, U, LDU, B, LDB, ZERO,
     $           WORK( N**2+1 ), N )
*
      CALL SGEMM( 'N', 'C', N, N, N, -ONE, WORK( N**2+1 ), N, U, LDU,
     $           ONE, WORK, N )
*
      WNORM = SLANGE( '1', N, N, WORK, N, WORK( N**2+1 ) )
*
      IF( ANORM.GT.WNORM ) THEN
         RESULT( 1 ) = ( WNORM / ANORM ) / ( N*ULP )
      ELSE
         IF( ANORM.LT.ONE ) THEN
            RESULT( 1 ) = ( MIN( WNORM, N*ANORM ) / ANORM ) / ( N*ULP )
         ELSE
            RESULT( 1 ) = MIN( WNORM / ANORM, REAL( N ) ) / ( N*ULP )
         END IF
      END IF
*
*     Do Test 2
*
*     Compute  UU' - I
*
      CALL SGEMM( 'N', 'C', N, N, N, ONE, U, LDU, U, LDU, ZERO, WORK,
     $           N )
*
      DO 40 JDIAG = 1, N
         WORK( ( N+1 )*( JDIAG-1 )+1 ) = WORK( ( N+1 )*( JDIAG-1 )+
     $         1 ) - ONE
   40 CONTINUE
*
      RESULT( 2 ) = MIN( SLANGE( '1', N, N, WORK, N,
     $              WORK( N**2+1 ) ), REAL( N ) ) / ( N*ULP )
*
      RETURN
*
*     End of SGET21
*
      END
*
* 
****** begin of atqrwk.f *********************************************
*
      SUBROUTINE ATQRWK( N, L, M, Q, LDQ, AQ, LDA )
*     .. 
*     .. Scalar Arguments ..
      INTEGER            L, LDA, LDQ, M, N
*     ..
*     .. Array Arguments ..
      REAL               AQ( LDA, * ), Q( LDQ, * )
*     ..
*
*  Purpose
*  =======
*
*  Compute
*
*               AQ(:,L:M) = A*Q(:,L:M)
*
*  Note that we only need to access the matrix A in question via
*  matrix-vector product form. Therefore, any special structure
*  and/or sparsity of matrix A can be exploited in forming such
*  matrix-vector products.
*
*  The matrix A in this subroutine is for the random walk example in
*  the manuscript ``SRRIT - A Fortran Subroutine to Calculate the
*  Dominant Invariant subspace of a Nonsymmetric matrix'' by Z. Bai
*  and G. W. Stewart.
*
*  The order of the random walk test matrix is
*
*        N = (K+1)*(K+2)/2
*
*  where K is the size of triangular random walk grid.
*
*  Arguments
*  =========
*
*  N    (input) INTEGER
*       The order of the matrix A.
*
*  L, M (input) INTEGERS
*       The first and last columns of Q to multiply by the matrix Q.
*
*  Q    (input) REAL array, dimension ( LDQ, M )
*       Q contains the matrix Q.
*
*  AQ   (output) REAL array, dimension (LDQ, M )
*       On return, columns L through M of AQ should contain the product
*       of the matrix A with columns L through M of Q.
*
*  =====================================================================
*
*     .. Parameters ..
      REAL               ZERO, HALF
      PARAMETER          ( ZERO = 0.0E+0, HALF = 0.5E+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            I, J, K, KK, LK, NN
      REAL               UPPER
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          REAL, SQRT
*     ..
*     .. Statement Functions ..
      INTEGER            IINDEX
*     ..
*     .. Statement Function definitions ..
      IINDEX( K, I, J ) = I*K - ( ( I-1 )*( I-2 ) ) / 2 + 1 + J + 1
*     ..
*     .. Executable Statements ..
*
      NN = ( -3+SQRT( REAL( 8*N+1 ) ) ) / 2
      DO 30 KK = L, M
*
*        Compute the vector AQ(1:N,KK) = A*Q(1:N,KK)
*
         DO 20 I = 0, NN
            DO 10 J = 0, NN - I
*
               K = IINDEX( NN, I, J )
               UPPER = REAL( I+J ) / REAL( 2*NN )
               AQ( K, KK ) = ZERO
*
               IF( I.EQ.0 .AND. J.EQ.0 ) THEN
                  LK = IINDEX( NN, I, J+1 )
                  AQ( K, KK ) = AQ( K, KK ) + ( HALF-UPPER )*Q( LK, KK )
                  LK = IINDEX( NN, I+1, J )
                  AQ( K, KK ) = AQ( K, KK ) + ( HALF-UPPER )*Q( LK, KK )
*
               ELSE IF( I.EQ.0 .AND. J.LE.NN-1 ) THEN
                  LK = IINDEX( NN, I, J-1 )
                  AQ( K, KK ) = AQ( K, KK ) + 2*UPPER*Q( LK, KK )
                  LK = IINDEX( NN, I+1, J )
                  AQ( K, KK ) = AQ( K, KK ) + ( HALF-UPPER )*Q( LK, KK )
                  LK = IINDEX( NN, I, J+1 )
                  AQ( K, KK ) = AQ( K, KK ) + ( HALF-UPPER )*Q( LK, KK )
*
               ELSE IF( I.EQ.0 .AND. J.EQ.NN ) THEN
                  LK = IINDEX( NN, I, J-1 )
                  AQ( K, KK ) = AQ( K, KK ) + 2*UPPER*Q( LK, KK )
*
               ELSE IF( I.LE.NN-1 .AND. J.EQ.0 ) THEN
                  LK = IINDEX( NN, I-1, J )
                  AQ( K, KK ) = AQ( K, KK ) + 2*UPPER*Q( LK, KK )
                  LK = IINDEX( NN, I+1, J )
                  AQ( K, KK ) = AQ( K, KK ) + ( HALF-UPPER )*Q( LK, KK )
                  LK = IINDEX( NN, I, J+1 )
                  AQ( K, KK ) = AQ( K, KK ) + ( HALF-UPPER )*Q( LK, KK )
*
               ELSE IF( I.EQ.NN .AND. J.EQ.0 ) THEN
                  LK = IINDEX( NN, I-1, J )
                  AQ( K, KK ) = AQ( K, KK ) + 2*UPPER*Q( LK, KK )
*
               ELSE IF( ( I+J ).EQ.NN ) THEN
                  LK = IINDEX( NN, I-1, J )
                  AQ( K, KK ) = AQ( K, KK ) + UPPER*Q( LK, KK )
                  LK = IINDEX( NN, I, J-1 )
                  AQ( K, KK ) = AQ( K, KK ) + UPPER*Q( LK, KK )
*
               ELSE
                  LK = IINDEX( NN, I-1, J )
                  AQ( K, KK ) = AQ( K, KK ) + UPPER*Q( LK, KK )
                  LK = IINDEX( NN, I, J-1 )
                  AQ( K, KK ) = AQ( K, KK ) + UPPER*Q( LK, KK )
                  LK = IINDEX( NN, I, J+1 )
                  AQ( K, KK ) = AQ( K, KK ) + ( HALF-UPPER )*Q( LK, KK )
                  LK = IINDEX( NN, I+1, J )
                  AQ( K, KK ) = AQ( K, KK ) + ( HALF-UPPER )*Q( LK, KK )
*
               END IF
   10       CONTINUE
   20    CONTINUE
*
   30 CONTINUE
*
      RETURN
*
*     End of ATQRWK
      END
*
****** begin of atqpde.f **********************************************
      SUBROUTINE ATQPDE( N, L, M, Q, LDQ, AQ, LDA )
*     ..
*     .. Scalar Arguments ..
      INTEGER            L, LDA, LDQ, M, N
*     ..
*     .. Array Arguments ..
      REAL               AQ( LDA, * ), Q( LDQ, * )
*     ..
*
*  Purpose
*  =======
*
*  Compute
*
*               AQ(:,L:M) = A*Q(:,L:M)
*
*  Note that we only need to access the matrix A in question via
*  matrix-vector product form. Therefore, any special structure
*  and/or sparsity of matrix A can be exploited in forming such
*  matrix-vector products.
*
*  The matrix A in this subroutine is for the convection-diffusion
*  PDE example in the manuscript ``SRRIT - A Fortran Subroutine to
*  Calculate the Dominant Invariant subspace of a Nonsymmetric
*  matrix'' by Z. Bai and G. W. Stewart.
*
*  The convection-diffusion operator is discretized on a K x K square 
*  grid, and resulted a PDE test matrix of order 
*
*        N  = K^2.
*
*  The constant parameters p1, p2, and p3 in the convection-diffusion
*  operator can changed in this subroutine.
*
*  Arguments
*  =========
*
*  N    (input) INTEGER
*       The order of the matrix A.
*
*  L, M (input) INTEGERS
*       The first and last columns of Q to multiply by the matrix Q.
*
*  Q    (input) REAL array, dimension ( LDQ, M )
*       Q contains the matrix Q.
*
*  AQ   (output) REAL array, dimension (LDQ, M )
*       On return, columns L through M of AQ should contain the product
*       of the matrix A with columns L through M of Q.
*
*  ===================================================================
*
*     .. Parameters ..
      REAL               ONE, FOUR
      PARAMETER          ( ONE = 1.0E+0, FOUR = 4.0E+0 )
      REAL               P1, P2, P3
      PARAMETER          ( P1 = 1.0E+0, P2 = 1.0E+0, P3 = 1.0E+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            I, II, K, NS
      REAL               BETA, GAMMA, HSTEP, SIGMA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          REAL, SQRT
*     ..
*     .. Executable Statements ..
*
      NS = SQRT( REAL( N ) )
      HSTEP = ONE / ( REAL( NS )+ONE )
      SIGMA = P3*HSTEP*HSTEP
      GAMMA = P2*HSTEP
      BETA = P1*HSTEP
*
      DO 70 K = L, M
*
*        Compute the vector AQ(1:N,K) = A*Q(1:N,K)
*
         DO 60 I = 1, NS
*
            DO 10 II = ( I-1 )*NS + 1, I*NS
               AQ( II, K ) = ( FOUR-SIGMA )*Q( II, K )
   10       CONTINUE
*
            DO 20 II = ( I-1 )*NS + 1, I*NS - 1
               AQ( II, K ) = AQ( II, K ) + ( GAMMA-ONE )*Q( II+1, K )
   20       CONTINUE
*
            DO 30 II = ( I-1 )*NS + 1 + 1, I*NS
               AQ( II, K ) = AQ( II, K ) + ( -GAMMA-ONE )*Q( II-1, K )
   30       CONTINUE
*
            IF( I.GT.1 ) THEN
               DO 40 II = ( I-1 )*NS + 1, I*NS
                  AQ( II, K ) = AQ( II, K ) - ( BETA+ONE )*Q( II-NS, K )
   40          CONTINUE
            END IF
            IF( I.LT.NS ) THEN
               DO 50 II = ( I-1 )*NS + 1, I*NS
                  AQ( II, K ) = AQ( II, K ) + ( BETA-ONE )*Q( II+NS, K )
   50          CONTINUE
            END IF
*
   60    CONTINUE
*
   70 CONTINUE
*
      RETURN
*
*     End of ATQPDE
      END
*
********* begin of atqode.f *******************************************
      SUBROUTINE ATQODE( N, L, M, Q, LDQ, AQ, LDA )
*     ..
*     .. Scalar Arguments ..
      INTEGER            L, LDA, LDQ, M, N
*     ..
*     .. Array Arguments ..
      REAL               AQ( LDA, * ), Q( LDQ, * )
*     ..
*
*  Purpose
*  =======
*
*  Compute
*
*               AQ(:,L:M) = A*Q(:,L:M)
*
*  Note that we only need to access the matrix A in question via
*  matrix-vector product form. Therefore, any special structure
*  and/or sparsity of matrix A can be exploited in forming such
*  matrix-vector products.
*
*  The matrix A in this subroutine is for the ODE boundary value
*  problem example in the manuscript ``SRRIT - A Fortran Subroutine
*  to Calculate the Dominant Invariant subspace of a Nonsymmetric
*  matrix'' by Z. Bai and G. W. Stewart.
*
*  The order of the ODE test matrix is
*
*        N = K + 1,
*
*  where K is the number of subintervals of [0,1].
*
*  NOTE: N should be smaller than 1001, if one wants to test a larger
*  size of N, the size of working arrays A and U in this subroutine 
*  should be increased. 
*
*  Arguments
*  =========
*
*  N    (input) INTEGER
*       The order of the matrix A.
*
*  L, M (input) INTEGERS
*       The first and last columns of Q to multiply by the matrix Q.
*
*  Q    (input) REAL array, dimension ( LDQ, M )
*       Q contains the matrix Q.
*
*  AQ   (output) REAL array, dimension (LDQ, M )
*       On return, columns L through M of AQ should contain the product
*       of the matrix A with columns L through M of Q.
*
*  ===================================================================
*
*     .. Parameters ..
      REAL               ZERO, ONE, TWO
      PARAMETER          ( ZERO = 0.0E+0, ONE = 1.0E+0, TWO = 2.0E+0 )
      REAL               THREE, FOUR
      PARAMETER          ( THREE = 3.0E+0, FOUR = 4.0E+0 )
      REAL               GAMMA
      PARAMETER          ( GAMMA = 1.0E-2 )
      INTEGER            LDU
      PARAMETER          ( LDU = 1001 )
*     ..
*     .. Local Scalars ..
      INTEGER            I, INFO, J, K
      REAL               DD, DELTA, HSTEP, HSTEP2
*     ..
*     .. Local Arrays ..
      REAL               A( 2, LDU ), U( LDU, 1 )
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          REAL
*     ..
*     .. External Subroutines ..
      EXTERNAL           SPBTF2, SPBTRS
*     ..
*     .. Executable Statements ..
*
      HSTEP = ONE / REAL( N )
      HSTEP2 = HSTEP*HSTEP
*
      A( 1, 1 ) = ZERO
      DO 10 J = 2, N - 1
         A( 1, J ) = -ONE
   10 CONTINUE
*
      DO 20 J = 1, N - 1
         A( 2, J ) = TWO
   20 CONTINUE
*
*     Cholesky decomposition of Symmetric tridiagonal matrix
*
      CALL SPBTF2( 'Upper', N-1, 1, A, 2, INFO )
*
      DO 30 J = 1, N - 2
         U( J, 1 ) = ZERO
   30 CONTINUE
      U( N-1, 1 ) = ONE
*
      CALL SPBTRS( 'Upper', N-1, 1, 1, A, 2, U, 301, INFO )
*
      DELTA = THREE*GAMMA + ( FOUR*U( 1, 1 )-U( 2, 1 )+GAMMA*
     $        U( N-2, 1 )-FOUR*GAMMA*U( N-1, 1 ) )
*
      DO 60 K = L, M
*
*        Compute the vector AQ(1:N,K) = A*Q(1:N,K)
*
         DO 40 I = 1, N - 1
            AQ( I, K ) = HSTEP2*Q( I, K )
   40    CONTINUE
         AQ( N, K ) = ZERO
*
         CALL SPBTRS( 'Upper', N-1, 1, 1, A, 2, AQ( 1, K ), LDA, INFO )
*
         DD = FOUR*AQ( 1, K ) - AQ( 2, K ) + GAMMA*AQ( N-2, K ) -
     $        FOUR*GAMMA*AQ( N-1, K )
         AQ( N, K ) = DD / DELTA
*
         DO 50 I = 1, N - 1
            AQ( I, K ) = -AQ( I, K ) + U( I, 1 )*AQ( N, K )
   50    CONTINUE
*
   60 CONTINUE
*
      DO 80 K = L, M 
         DO 70 I = 1, N
               AQ( I, K ) = -AQ( I, K )
   70    CONTINUE
   80 CONTINUE
*
      RETURN
*
*     End of ATQODE
      END
********* begin of atqbwm.f ******************************************
      SUBROUTINE ATQBWM( N, L, M, X, LDX, Y, LDY )
*     ..
*     .. Scalar Arguments ..
      INTEGER            LDY, LDX, L, M, N
*     ..
*     .. Array Arguments ..
      REAL               Y( LDY, * ), X( LDX, * )
*     ..
*
*  Purpose
*  =======
*
*  Compute
*
*               Y(:,L:M) = A*X(:,L:M)
*
*  The matrix A is a 2 by 2 block matrix resulted from the finite 
*  difference discretization of the Brusselator wave model. 
*
*  Arguments
*  =========
*
*  N       (input) INTEGER
*          The order of the matrix A. N has to be an even number.
*
*  L,M     (input) INTEGERS
*          The number of the first and the last column of Q to 
*          multiply by A.
*
*  X       (input) REAL             array, dimension ( LDX, M )
*          contains the vectors X.
*
*  LDX     (input) INTEGER
*          The leading dimension of the array X, LDX >= max( 1, N )
*
*  Y       (output) REAL             array, dimension (LDY, M )
*          contains the product of the matrix op(A) with X.
*
*  LDY     (input) INTEGER
*          The leading dimension of the array Y, LDY >= max( 1, N )
*
*  ===================================================================
*
*     .. Parameters ..
      REAL               ONE, TWO
      PARAMETER          ( ONE = 1.0E+0, TWO = 2.0E+0 )
      REAL               ALPHA, BETA, DELT1, DELT2, LL
      PARAMETER          ( ALPHA = 2.0E+0, BETA = 5.45E+0, 
     $                     DELT1 = 8.0E-3, DELT2 = 4.0E-3,
     $                     LL = 0.51302E+0 )
*
*     .. Local Variables .. 
      INTEGER            I, K, NHALF
      REAL               HSTEP, TAU1, TAU2, DIAG1, DIAG2, ALPHA2 
*
*     set up middle parameters
*
      NHALF = N / 2 
      HSTEP = ONE / ( REAL( NHALF ) + 1 ) 
      TAU1 = DELT1 / (HSTEP*LL)**2  
      TAU2 = DELT2 / (HSTEP*LL)**2  
      DIAG1 = -TWO*TAU1 + BETA - ONE 
      ALPHA2 = ALPHA**2 
      DIAG2 = -TWO*TAU2 - ALPHA2 
* 
*     Compute Y(:,L:M) = A*X(:,L:M)  
*
      DO 30 K = L, M 
*
*        the 1st row 
*
         Y( 1, K ) = DIAG1*X( 1, K ) + TAU1*X( 2, K ) + 
     $               ALPHA2*X( NHALF+1, K )    
*
*        the 2nd row to (NHALF-1)th row              
*
         DO 10 I = 2, NHALF - 1 
            Y( I, K ) = TAU1*X( I-1, K ) + DIAG1*X( I, K ) +
     $               TAU1*X( I+1, K ) + ALPHA2*X( NHALF+I, K )   
   10    CONTINUE  
*
*        the NHALFth row
*
         Y( NHALF, K ) = TAU1*X( NHALF-1, K ) + DIAG1*X( NHALF, K ) 
     $                   + ALPHA2*X( N,K )  
*          
*        the (NHALF+1)th row
*
         Y( NHALF+1, K ) = -BETA*X( 1, K ) + DIAG2*X( NHALF+1, K ) 
     $                     + TAU2*X( NHALF+2, K )  
*
*        the (NHALF+2)th row to (N-1)th row
*
         DO 20 I = NHALF+2, N - 1 
            Y( I, K ) = -BETA*X( I-NHALF, K ) + TAU2*X( I-1, K ) + 
     $                  DIAG2*X( I, K ) + TAU2*X( I+1, K ) 
   20    CONTINUE     
*
*        the last (Nth) row  
*
         Y( N, K ) = -BETA*X( N-NHALF, K ) + TAU2*X( N-1,K ) + 
     $               DIAG2*X( N, K ) 
*
   30 CONTINUE 
*
      RETURN
*
*     End of ATQBWM
*
      END
********* begin of atq4x4.f ******************************************
      SUBROUTINE ATQ4X4( N, L, M, Q, LDQ, AQ, LDA )
*     ..
*     .. Scalar Arguments ..
      INTEGER            L, LDA, LDQ, M, N
*     ..
*     .. Array Arguments ..
      REAL               AQ( LDA, * ), Q( LDQ, * )
*
*  Purpose
*  =======
*
*  Compute
*
*               AQ(:,L:M) = A*Q(:,L:M)
*
*  Note that we only need to access the matrix A in question via
*  matrix-vector product form. Therefore, any special structure
*  and/or sparsity of matrix A can be exploited in forming such
*  matrix-vector products.
* 
*  The matrix A in this subroutine is a 3 by 3 matrix: 
*
*            A = [ 2  1  1 1 ] 
*                [ 0  2  1 1 ]
*                [ 0  0  2 1 ]
*                [ 0  0    1 ]
*
*  which has exact eigenvalues 2, 2, 2, 1. 
* 
*  Arguments
*  =========
*
*  N    (input) INTEGER
*       The order of the matrix A.
*
*  L, M (input) INTEGERS
*       The first and last columns of Q to multiply by the matrix Q.
*
*  Q    (input) REAL array, dimension ( LDQ, M )
*       Q contains the matrix Q.
*
*  AQ   (output) REAL array, dimension (LDQ, M )
*       On return, columns L through M of AQ should contain the product
*       of the matrix A with columns L through M of Q.
*
*  ===================================================================
*
*     .. Parameters ..
      REAL               TWO 
      PARAMETER          ( TWO = 2.0E+0 )
*
*     ..
*     .. Local Scalars ..
      INTEGER            K
*
      N = 4 
      DO 10 K = L, M 
*
*        Compute the vector AQ(1:N,K) = A*Q(1:N,K)
*
         AQ( 1, K ) = TWO*Q( 1, K ) + Q( 2, K ) + Q( 3, K ) + Q( 4, K )
         AQ( 2, K ) = TWO*Q( 2, K ) + Q( 3, K ) + Q( 4, K )
         AQ( 3, K ) = TWO*Q( 3, K ) + Q( 4, K )
         AQ( 4, K ) = Q( 4, K )
*
10    CONTINUE 
*
      RETURN
      END 
      SUBROUTINE ATQF44( N, L, M, Q, LDQ, AQ, LDA )
*     ..
*     .. Scalar Arguments ..
      INTEGER            L, LDA, LDQ, M, N
*     ..
*     .. Array Arguments ..
      REAL               AQ( LDA, * ), Q( LDQ, * )
*
*  Purpose
*  =======
*
*  Compute
*
*               AQ(:,L:M) = A*Q(:,L:M)
*
*  Note that we only need to access the matrix A in question via
*  matrix-vector product form. Therefore, any special structure
*  and/or sparsity of matrix A can be exploited in forming such
*  matrix-vector products.
* 
*  The matrix A in this subroutine is a 3 by 3 matrix: 
*
*            A = [ 2  0  1 2 ] 
*                [ 0  0  1 1 ]
*                [ 0  0  2 1 ]
*                [ 0  0  1 1 ]
* 
*  Arguments
*  =========
*
*  N    (input) INTEGER
*       The order of the matrix A.
*
*  L, M (input) INTEGERS
*       The first and last columns of Q to multiply by the matrix Q.
*
*  Q    (input) REAL array, dimension ( LDQ, M )
*       Q contains the matrix Q.
*
*  AQ   (output) REAL array, dimension (LDQ, M )
*       On return, columns L through M of AQ should contain the product
*       of the matrix A with columns L through M of Q.
*
*  ===================================================================
*
*     .. Parameters ..
      REAL               TWO 
      PARAMETER          ( TWO = 2.0E+0 )
*
*     ..
*     .. Local Scalars ..
      INTEGER            K
*
      N = 4
      DO 10 K = L, M 
*
*        Compute the vector AQ(1:N,K) = A*Q(1:N,K)
*
         AQ( 1, K ) = TWO*Q( 1, K ) + Q( 3, K ) + TWO*Q( 4, K )
         AQ( 2, K ) = Q( 3, K ) + Q( 4, K )
         AQ( 3, K ) = TWO*Q( 3, K ) + Q( 4, K )
         AQ( 4, K ) = Q( 3, K ) + Q( 4, K )
*
10    CONTINUE 
*
      RETURN
      END 
