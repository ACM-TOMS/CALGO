      SUBROUTINE DSRRIT( N, NV, M, MAXIT, ISTART, Q, LDQ, AQ, LDA, T,
     $                   LDT, WR, WI, RSD, ITRSD, IWORK, WORK, LWORK,
     $                   INFO, EPS, ATQ )
*
* ====================================================================
*
*     Copyright (C) 1994, Zhaojun Bai and G. W. Stewart
*     All rights reserved.  Use at your own risk.
*     Please send comments to bai@ms.uky.edu or stewart@cs.umd.edu
*
* ====================================================================
*     ..
*     .. Scalar Arguments ..
      INTEGER            INFO, ISTART, LDA, LDQ, LDT, LWORK, M, MAXIT,
     $                   N, NV
      DOUBLE PRECISION   EPS
*     ..
*     .. Array Arguments ..
      INTEGER            ITRSD( * ), IWORK( * )
      DOUBLE PRECISION   AQ( LDA, * ), Q( LDQ, * ), RSD( * ),
     $                   T( LDT, * ), WI( * ), WORK( * ), WR( * )
*     ..
*     .. Subroutine Arguments ..
      EXTERNAL           ATQ
*     ..
*
*  Purpose
*  =======
*
*  SRRIT computes a nested sequence of orthonormal bases for the
*  dominant invariant subspaces of a real matrix A of order N.
*  Specifically, the program returns an N x NV matrix Q with
*  orthonormal columns and an NV x NV quasi-triangular matrix T
*  satisfying
*
*            A*Q = Q*T + O(eps)
*
*  The eigenvalues of T in 1 x 1 diagonal blocks are real, in the 2 x 2
*  diagonal blocks are complex conjugate pairs. The  eigenvalues of T:
*  W(J) = WR(J) + i*WI(J), J = 1, 2, ... , NV, are ordered so that
*
*         |W(1)| >= |W(2)| >= ... >= |W(NV)|.
*
*  These eigenvalues approximate the dominant eigenvalues of A, i.e.,
*  the eigenvalues of largest absolute magnitudes. These facts have the
*  following consequences
*
*     1. If W(L).ne.W(L+1) and W(L).ne.conj(W(L+1)), then columns
*        1,2,...,L of Q form an approximate basis for the invariant
*        subspace corresponding to the L dominant eigenvalues of A.
*        The L by L leading principal submatrix of T is a representation       
*        of A in that subspace with respect to the basis Q.
*
*     2. If z is an eigenvector of T corresponding to W(J), then Q*z is
*        an approximate eigenvector of A corresponding to W(J). z may
*        be computed by calling LAPACK subroutine STREVC.
*
*  The program actually iterates with an N x M matrix Q and an M x M
*  matrix T. The rate of convergence of the L-th column of Q is
*  essentially linear with ratio |W(M+1)/W(L)|. The program may fail to
*  converge after the user specified the maximal number of iterations, 
*  MAXIT, see the arguments NV, MAXIT and INFO. In general, it pays that
*  the user sets M larger than the desired number, NV, of eigenvalues, 
*  say M = NV + 5 or larger.
*
*  The user must furnish a subroutine to compute the product A*Q.
*  The calling sequence is
*
*        ATQ( N, L, M, Q, LDQ, AQ, LDA )
*
*  For J = L, L+1, ..., M, the program must place the product A*Q(:,J)
*  in AQ(:,J).
*
*  The method is based on the paper by G. W. Stewart 'Simultaneous
*  iteration for computing invariant subspaces of non-hermitian
*  matrices', Numer, Math. 25(1976), pp.123-136.
*
*  For the further details of algorithm implementation and usage, see 
*
*  [1] Z. Bai and G. W. Stewart, SRRIT - A Fortran subroutine to 
*  calculate the dominant invariant subspace of a nonsymmetric matrix, 
*  Technical Report TR-2908, Department of Computer Science, 
*  University of Maryland, 1992 (Submitted to ACM TOMS for 
*  publication).
*
*  Arguments
*  =========
*
*  N      (input) INTEGER
*         N is the order of the matrix A.
*
*  NV     (input/output) INTEGER
*         One entry, NV is the size of the dominant invariant subspace 
*         of A that the user desired. On exit, NV is the actual size of 
*         converged invariant subspace. 
*         Note that if input NV is larger than output NV, say, on exit, 
*         NV = 0, then it means that the method fails to converge to 
*         the number of desired eigenvalues after the MAXIT number of 
*         iterations.
*
*  M      (input) INTEGER
*         M is the size of iteration space. NV <= M <= N.
*
*  MAXIT  (input/output) INTEGER
*         On entry, MAXIT is an upper bound on the number of iterations
*         the program is to execute. On exit, MAXIT is the actual number       
*         of iteration the program executed (also see NV). 
*
*  ISTART (input) INTEGER
*         ISTART specifies whether user supply initial Q.
*           < 0, Q is initialized by the program, random vectors with
*                uniform (0,1) distribution. 
*           = 0, starting Q has been set in the input, but
*                it is not orthonormal.
*           > 0, starting Q has been set in the input, it is also
*                orthonormal.
*
*  Q      (input/output) DOUBLE PRECISION array, dimension( LDQ,M )
*         On entry, if ISTART > 0, Q contains the starting Q which will
*         be used in the simultaneous iteration.
*         On exit, Q contains the orthonormal vectors described above.
*
*  LDQ    (input) INTEGER
*         The leading dimension of Q, LDQ >= max(1,N).
*
*  AQ     (output) DOUBLE PRECISION array, dimension( LDA, M)
*         On exit, AQ contains the product A*Q.
*
*  LDA    (input) INTEGER
*         The leading dimension of AQ, LDA >= max(1,N).
*
*  T      (output) DOUBLE PRECISION array, dimension( LDT,M )
*         On exit, T contains of representation of A described above.
*
*  LDT    (input) INTEGER
*         The leading dimension of T, LDT >= max(1,M).
*
*  WR,WI  (output) DOUBLE PRECISION arrays, dimension (M)
*         On exit, WR and WI contain the real and imaginary parts, resp.
*         of the eigenvalues of T, which are the approximations of the
*         dominant eigenvalues of matrix A. The eigenvalues are ordered
*         in decreasing.
*
*  RSD    (output) DOUBLE PRECISION arrays, dimension( M )
*         On exit, RSD contains the 2-norm of the residual vectors
*         R(:,J) = A*Q(:,J) - Q*T(:,J) for J = 1,...,M.
*
*  ITRSD  (output) INTEGER array, dimension( M )
*         On exit, ITRSD contains the iteration numbers at which
*         the residuals were computed.
*
*  IWORK  (workspace) INTEGER array, dimension( 2*M )
*
*  WORK   (workspace) DOUBLE PRECISION array, WORK( LWORK )
*
*  LWORK  (input) INTEGER
*         The length of workspace WORK, LWORK >= M*M + 5*M
*         Note: the first M*M workspace is used for the matrix U
*         in DSRRSP, next 2*M is the workspace for DSRRSP, additional
*         2*M for storing old WR and WI, and finally, M for storing
*         old RSD in this program.
*
*  INFO   (input/output) INTEGER
*         On entry, if INFO is set to be a nonzero integer, then
*         the information about the course of the iteration will be
*         printed, otherwise, it is omitted.
*         On exit, if INFO is set to
*            0,   normal return. 
*           -k,   if input argument number k is illegal.
*            1,   Orthonormalization subroutine DORTH fails. 
*            2,   error from subroutine DSRRSP.
*            3,   matrix T is singular, see DCOND.
*            4,   method fails to converge the number of desired 
*                 eigenvalues (NV) after MAXIT iterations. RSD 
*                 contains the 2-norm of residual vectors. 
*
*  EPS    (Input) DOUBLE PRECISION
*         A convergence criterion supplied by user.
*
*  ATQ    External procedure name, which passes the matrix-vectors
*         product operation.
*
*
*  Further details: 
*  ================
*
*  The internal control parameters INIT, STPFAC, ALPHA, BETA, GRPTOL, 
*  CNVTOL and ORTTOL are described in reference [1]. Currently, they 
*  are set as default initial values.  
*
* =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE, TWO
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0, TWO = 2.0D+0 )
      DOUBLE PRECISION   STPFAC, ALPHA, BETA
      PARAMETER          ( STPFAC = 1.5, ALPHA = ONE, BETA = 1.1D+0 )
      DOUBLE PRECISION   ORTTOL
      PARAMETER          ( ORTTOL = TWO )
      DOUBLE PRECISION   GRPTOL, CNVTOL
      PARAMETER          ( GRPTOL = 1.0D-8, CNVTOL = 1.0D-6 )
*
      INTEGER            INIT, NOUT
      PARAMETER          ( INIT = 5, NOUT = 6 )
*     ..
*     .. Local Scalars ..
      LOGICAL            PINFO
      INTEGER            I, IDORT, IDSRR, IERR, IT, J, L, M2, NGRP,
     $                   NOGRP, NV0, NXTORT, NXTSRR
      DOUBLE PRECISION   AE, AQMAX, ARSD, CTR, OAE, OARSD, OCTR, REC,
     $                   TCOND
*     ..
*     .. Local Arrays ..
      INTEGER            ISEED( 4 )
      DOUBLE PRECISION   DUMMY( 1 )
*     ..
*     .. External Functions ..
      DOUBLE PRECISION   DCOND, DLANGE, DLARAN
      EXTERNAL           DCOND, DLANGE, DLARAN
*     ..
*     .. External Subroutines ..
      EXTERNAL           DGROUP, DORTH, DRESID, DSRRSP, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, INT, LOG, LOG10, MAX, MIN, DBLE
*     ..
*     .. Executable Statements ..
*
*     Decide whether print the middle computed results.
*
      PINFO = .TRUE.
      IF( INFO.EQ.0 )
     $   PINFO = .FALSE.
*
*     Test the input arguments
*
      INFO = 0
      IF( N.LT.0 ) THEN
         INFO = -1
      ELSE IF( NV.GT.N .OR. NV.GT.M ) THEN
         INFO = -2
      ELSE IF( M.GT.N ) THEN
         INFO = -3
      ELSE IF( LDQ.LT.MAX( 1, N ) ) THEN
         INFO = -7
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -9
      ELSE IF( LDT.LT.MAX( 1, M ) ) THEN
         INFO = -11
      ELSE IF( LWORK.LT.( M*M+5*M ) ) THEN
         INFO = -18
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DSRRIT', -INFO )
         RETURN
      END IF
*
*     Initialize
*
      IT = 0
      L = 1
      M2 = M*M + 2*M
      DO 10 J = 1, M
         WR( J ) = ZERO
         WI( J ) = ZERO
         RSD( J ) = ZERO
         ITRSD( J ) = -1
   10 CONTINUE
      NV0 = NV
*
      IF( ISTART.GE.0 )
     $   GO TO 40
*
*     Set initial random vectors with uniform (0,1) distribution.
*
      ISEED( 1 ) = 1978
      ISEED( 2 ) = 1989
      ISEED( 3 ) = 1992
      ISEED( 4 ) = 1993
      DO 30 J = 1, M
         DO 20 I = 1, N
            Q( I, J ) = DLARAN( ISEED )
   20    CONTINUE
   30 CONTINUE
*
*     Use the input vectors as initial vectors. 
*
   40 CONTINUE
      IF( ISTART.GT.0 )
     $   GO TO 50
*
*     Orthogonalize Q if necessary
*
      CALL DORTH( N, 1, M, Q, LDQ, IERR )
      IF( IERR.NE.0 )THEN 
         INFO = 1
         RETURN
      END IF
*
   50 CONTINUE
*
*     The main Schur-Rayleigh-Ritz (SRR) loop
*
   60 CONTINUE
*
      DO 70 I = L, M
         WORK( M2+I ) = WR( I )
         WORK( M2+M+I ) = WI( I )
   70 CONTINUE
*
      IF( PINFO ) THEN
         WRITE( NOUT, FMT = 9999 )IT
 9999    FORMAT( 30X, 'IT = ', I4 )
      END IF
*
      CALL DSRRSP( N, L, M, Q, LDQ, AQ, LDA, T, LDT, WR, WI, WORK, M,
     $             WORK( M*M+1 ), IERR, ATQ )
      IF( IERR.NE.0 )THEN 
         INFO = 2
         RETURN
      END IF 
*
      DO 80 I = L, M
         WORK( M2+2*M+I ) = RSD( I )
         IWORK( M+I ) = ITRSD( I )
         ITRSD( I ) = IT
   80 CONTINUE
*
*     Compute the 2-norms of residual vectors: 
*      || R(1:M,L:M) || = || A*Q(1:M,L:M) - Q*T(1:M,L:M) ||
*
      CALL DRESID( N, L, M, Q, LDQ, AQ, LDA, T, LDT, RSD )
*
      IF( PINFO ) THEN
         WRITE( NOUT, FMT = 9998 )( WR( I ), I = 1, M )
 9998    FORMAT( ' WR   =', 1P,6D12.3 )
         WRITE( NOUT, FMT = 9997 )( WI( I ), I = 1, M )
 9997    FORMAT( ' WI   =', 1P,6D12.3 )
         WRITE( NOUT, FMT = 9996 )( RSD( I ), I = 1, M )
 9996    FORMAT( ' RSD  =', 1P,6D12.3 )
      END IF
*
   90 CONTINUE
*
*     Find clusters of computed eigenvalues
*
      CALL DGROUP( L, M, WR, WI, RSD, NGRP, CTR, AE, ARSD, GRPTOL )
      CALL DGROUP( L, M, WORK( M2+1 ), WORK( M2+M+1 ), WORK( M2+2*M+1 ),
     $             NOGRP, OCTR, OAE, OARSD, GRPTOL ) 
*
      IF( PINFO ) THEN
         WRITE( NOUT, FMT = 9995 )NGRP
 9995    FORMAT( ' NGRP =', I5 )
         WRITE( NOUT, FMT = 9994 )CTR, AE, ARSD
 9994    FORMAT( ' CTR  =', 1P,D12.3, '  AE =', 1P,D12.3, 
     $'  ARSD = ', 1P,D12.3 )
      END IF
*
      IF( NGRP.NE.NOGRP )
     $   GO TO 100
      IF( NGRP.EQ.0 )
     $   GO TO 100
*
      IF( ABS( AE-OAE ).GT.CTR*CNVTOL*DBLE( ITRSD( L )-IWORK( M+L ) ) )
     $   GO TO 100
*
*     Test for the convergence.
*
      IF( ARSD.GT.CTR*EPS )
     $   GO TO 100
      L = L + NGRP
      IF( L.GT.M )
     $   GO TO 100
*
      GO TO 90
  100 CONTINUE
*
*     Exit if the required number of eigenvalues have converged.
*
      IF( L.GT.NV )
     $   GO TO 210
*
*     Exit if iteration count exceeds the maximum number of iterations
*
      IF( IT.GE.MAXIT )
     $   GO TO 210
*
*     Determine when the next SRR step (NXTSRR) is to be taken.
*
      NXTSRR = MIN( MAXIT, MAX( INT( STPFAC*DBLE( IT ) ), INIT ) )
      IDSRR = NXTSRR - IT
      IF( NGRP.NE.NOGRP )
     $   GO TO 110
      IF( NGRP.EQ.0 )
     $   GO TO 110
      IF( ARSD.GE.OARSD )
     $   GO TO 110
      REC = ALPHA + BETA*DBLE( IWORK( M+L )-ITRSD( L ) )*
     $      LOG( ARSD / EPS ) / LOG( ARSD / OARSD )
      IDSRR = MAX( 1, INT( REC ) )
*
  110 CONTINUE
      NXTSRR = MIN( NXTSRR, IT+IDSRR )
*
*     Determine the interval between orthogonalizations.
*     The next orthogonalization step at NXTORT
*
      DO 130 J = 1, M
         DO 120 I = 1, M
            WORK( I+( J-1 )*M ) = T( I, J )
  120    CONTINUE
  130 CONTINUE
      TCOND = DCOND( M, WORK, M, IWORK, WORK( M*M+1 ), IERR )
      IF( IERR.NE.0 )THEN 
         INFO = 3
         RETURN
      END IF
      IDORT = MAX( 1, INT( ORTTOL / MAX( ONE, LOG10( TCOND ) ) ) )
      NXTORT = MIN( IT+IDORT, NXTSRR )
*
      IF( PINFO ) THEN
         WRITE( NOUT, FMT = 9993 )NXTSRR, IDORT
 9993    FORMAT( ' NXTSRR =', I4, '  IDORT  =', I4 )
      END IF
*
      DO 150 J = L, M
         DO 140 I = 1, N
            Q( I, J ) = AQ( I, J )
  140    CONTINUE
  150 CONTINUE
      IT = IT + 1
*
  160 CONTINUE
*
*     Orthogonalization loop
*
*     Power loop. After each matrix-vector step, we scale the
*     product to prevent the possible overflow.
*
  170 CONTINUE
      IF( IT.EQ.NXTORT )
     $   GO TO 200
*
      CALL ATQ( N, L, M, Q, LDQ, AQ, LDA )
      AQMAX = DLANGE( 'M', N, M-L+1, AQ( 1, L ), LDA, DUMMY )
      REC = ONE / AQMAX
      DO 190 J = L, M
         DO 180 I = 1, N
            Q( I, J ) = AQ( I, J )*REC
  180    CONTINUE
  190 CONTINUE
      IT = IT + 1
      GO TO 170
*
  200 CONTINUE
      CALL DORTH( N, L, M, Q, LDQ, IERR )
      IF( IERR.NE.0 )THEN
         INFO = 1
         RETURN
      END IF 
      NXTORT = MIN( IT+IDORT, NXTSRR )
      IF( IT.LT.NXTSRR )
     $   GO TO 160
*
      GO TO 60
*
  210 CONTINUE
*
*     Final number of converged approx. eigenvalues, and number of 
*     total iterations (= the number of calls to ATQ).  
*
      NV = L - 1
      MAXIT = IT
*
      IF( NV.LT.NV0 )
     $   INFO = 4
*
      RETURN
*
*     End of DSRRIT
      END
*
******** begin of dsrrsp.f ********************************************        
      SUBROUTINE DSRRSP( N, L, M, Q, LDQ, AQ, LDA, T, LDT, WR, WI, U,
     $                   LDU, WORK, INFO, ATQ )
*     .. 
*     .. Scalar Arguments ..
      INTEGER            INFO, L, LDA, LDQ, LDT, LDU, M, N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   AQ( LDA, * ), Q( LDQ, * ), T( LDT, * ),
     $                   U( LDU, * ), WI( * ), WORK( * ), WR( * )
*     ..
*     .. Subroutine Arguments ..
      EXTERNAL           ATQ
*     ..
*
*  Purpose
*  =======
*
*  DRRSP performs a Schur-Rayleigh-Ritz refinement on the set of
*  M-L+1 orthogonal N-vectors contained in the array Q(:,L:M) in
*  the following three steps:
*
*       1). Compute T = Q(:,L:M)'*A*Q(:,L:M);
*       2). Compute Schur decomposition of T: T := U'*T*U;
*       3). Update Q(:,L:M)  = Q(:,L:M)*U;
*                  AQ(:,L:M) = AQ(:,L:M)*U.
*
*  Note that in the Schur decomposition of T, the eigenvalues in the 
*  Schur form are ordered in decreasing in terms of absolute magnitude.
*
*  Arguments
*  =========
*
*  N       (input) INTEGER
*          N is the order of the matrix A.
*
*  L, M    (input) INTEGER
*          The first and last column indices of Q for SRR iteration.
*
*  Q       (input/output) DOUBLE PRECISION array, dimension( LDQ,M )
*          On entry, Q contains approximate invariant subpace matrix.
*          On exit, Q is updated by Schur vectors of matrix T.
*
*  LDQ     (input) INTEGER
*          The leading dimension of Q, LDQ >= max(1,N).
*
*  AQ      (output) DOUBLE PRECISION array, dimension( LDA, M)
*          On exit, AQ contains the product A*Q.
*
*  LDA     (input) INTEGER
*          The leading dimension of AQ, LDA >= max(1,N).
*
*  T       (output) DOUBLE PRECISION array, dimension( LDT,M )
*          On exit, T contains of Schur form computed by SRR iteration.
*
*  LDT     (input) INTEGER
*          The leading dimension of T, LDT >= max(1,M).
*
*  WR,WI   (output) DOUBLE PRECISION arrays, dimension (M)
*          On exit, WR and WI contain the real and imaginary parts of 
*          the eigenvalues of T, respectively. The eigenvalues are 
*          ordered in decreasing.
*
*  U       (output) DOUBLE PRECISION array, dimension(M,M)
*          On exit, U contains the orthogonal matrix of the Schur
*          decomposition.
*
*  LDU     (input) INTEGER
*          The leading dimension of U, LDU >= max(1,M).
*
*  WORK    (workspace) DOUBLE PRECISION array, dimension( 2*M )
*
*  INFO    (output) INTEGER
*          On exit, if INFO is set to 
*            0,  normal return
*            1,  error from reduction T to Hessenberg form (DGEHD2).
*            2,  error from forming orthogonal matrix U (DORGN2).
*            3,  error from Schur decomposition routine (DLAQR3).
*
*  ATQ    External procedure name, which passes the matrix-vectors 
*         product operation.
*
* =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            I, IERR, J, K, NC
*     ..
*     .. External Subroutines ..
      EXTERNAL           DGEHRD, DGEMM, DLAQR3, DLASET, DORGHR
*     ..
*     .. Executable Statements ..
*
      INFO = 0
*
*     Calculate T(1:M,L:M) = Q(:,1:M)'*A*Q(:,L:M)
*
      NC = M - L + 1
      CALL ATQ( N, L, M, Q, LDQ, AQ, LDA )
      CALL DGEMM( 'Transpose', 'No transpose', M, NC, N, ONE, Q, LDQ,
     $            AQ( 1, L ), LDA, ZERO, T( 1, L ), LDT )
*
*     Hessenberg reduction 
*
      CALL DGEHRD( M, L, M, T, LDT, WORK, WORK( M+1 ), M, IERR )
      IF( IERR.NE.0 )THEN
         INFO = 1
         RETURN
      END IF
      DO 20 J = 1, M - 1
         DO 10 I = J + 2, M
            U( I, J ) = T( I, J )
            T( I, J ) = ZERO
   10    CONTINUE
   20 CONTINUE
      CALL DORGHR( M, L, M, U, LDU, WORK, WORK( M+1 ), M, IERR )   
*
*     Compute Schur decomposition
*
      CALL DLAQR3( .TRUE., .TRUE., M, L, M, T, LDT, WR, WI, 1, M, 
     $             U, LDU, WORK, IERR )
      IF( IERR.NE.0 )THEN
         INFO = 3
         RETURN
      END IF
*
*     Update: Q(:,L:M)  := Q(:,L:M)*U
*             AQ(:,L:M) := AQ(:,L:M)*U
*
      DO 60 I = 1, N
         DO 40 J = L, M
            WORK( J ) = ZERO
            WORK( J+M ) = ZERO
            DO 30 K = 1, M
               WORK( J ) = WORK( J ) + Q( I, K )*U( K, J )
               WORK( J+M ) = WORK( J+M ) + AQ( I, K )*U( K, J )
   30       CONTINUE
   40    CONTINUE
         DO 50 J = L, M
            Q( I, J ) = WORK( J )
            AQ( I, J ) = WORK( J+M )
   50    CONTINUE
   60 CONTINUE
*
      RETURN
*
*     End of DSRRSP
      END
*
******** begin of dorth.f **********************************************
      SUBROUTINE DORTH( N, L, M, Q, LDQ, INFO )
*     .. 
*     .. Scalar Arguments ..
      INTEGER            INFO, L, LDQ, M, N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   Q( LDQ, * )
*     ..
*
*  Purpose
*  =======
*
*  DORTH orthonormalizes columns L through M of the array Q with respect
*  to columns 1 to M. Column 1 through L-1 are assumed to be 
*  orthonormal. The method is the modified Gram-Schmidt method with 
*  reorthogonalization. A column is accepted when an orthogonalization 
*  does not reduce its 2-norm by a factor of more than TOL (0.5). If 
*  this is not done in MAXTRY (5) attempts the program stops. DORTH also
*  stops if it encounters a zero vector.
*
*  Arguments
*  =========
*
*  N       (input) INTEGER
*          The number of rows of Q.
*
*  L,M     (input) INTEGERs
*          The first and last column indices of Q which is
*          going to be orthonormalized.
*
*  Q       (input/output) DOUBLE PRECISION array, dimension(LDQ,M)
*          On entry, Q contains the vectors which is going to be
*          orthonormalized.
*          On exit, Q is overwritten by the orthonormalized vectors.
*
*  LDQ     (input) INTEGER
*          The leading dimension of array Q, LDQ >=max(1,N)
*
*  INFO    (output) INTEGER
*          On exit, if INFO is set to
*            0,  successful return.
*            1,  zero vectors in Q.
*            2,  reorthogonalization iteration fails after MAXTRY 
*                attempts. 
*
* =====================================================================
*
*     .. Parameters ..
      INTEGER            MAXTRY
      PARAMETER          ( MAXTRY = 5 ) 
      DOUBLE PRECISION   ZERO, ONE, TOL
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0, TOL = 0.5D+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            ITRY, J, JM1, K
      DOUBLE PRECISION   NORM, QQ, REC
*     ..
*     .. External Functions ..
      DOUBLE PRECISION   DDOT, DNRM2
      EXTERNAL           DDOT, DNRM2
*     ..
*     .. External Subroutines ..
      EXTERNAL           DAXPY, DSCAL
*     ..
*     .. Executable Statements ..
*
      INFO = 0
*
      DO 30 J = L, M
         JM1 = J - 1
         ITRY = 0
   10    CONTINUE
*
*        Compute the 2-norm
*
         NORM = DNRM2( N, Q( 1, J ), 1 )
         IF( NORM.EQ.ZERO )THEN
            INFO = 1 
            RETURN
         END IF 
*
*        Scale the vector
*
         REC = ONE / NORM
         CALL DSCAL( N, REC, Q( 1, J ), 1 )
*
*        Test to see if the J-th vector is orthogonal
*
         IF( J.EQ.1 )
     $      GO TO 30
         IF( ITRY.EQ.0 )
     $      NORM = ZERO
         IF( NORM.GT.TOL )
     $      GO TO 30
         ITRY = ITRY + 1
         IF( ITRY.GT.MAXTRY )THEN
            INFO = 2
            RETURN 
         END IF 
*
*        Perform one step of the modified Gram-Schmidt process:
*           Q(:,J) = ( I - Q(:,1:J-1)*Q(:,1:J-1)')*Q(:,J)  
*
         DO 20 K = 1, JM1
            QQ = DDOT( N, Q( 1, K ), 1, Q( 1, J ), 1 )
            CALL DAXPY( N, -QQ, Q( 1, K ), 1, Q( 1, J ), 1 )
   20    CONTINUE
         GO TO 10
*
   30 CONTINUE
*
      RETURN
*
*     End of DORTH
      END
*
******** begin of dresid.f *******************************************
      SUBROUTINE DRESID( N, L, M, Q, LDQ, AQ, LDA, T, LDT, RSD )
*     .. 
*     .. Scalar Arguments ..
      INTEGER            L, LDA, LDQ, LDT, M, N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   AQ( LDA, * ), Q( LDQ, * ), RSD( * ),
     $                   T( LDT, * )
*     ..
*
*  Purpose
*  =======
*
*  Computes the column norms of residual vectors
*
*             A*Q(1:N,L:M) - Q*T(1:M,L:M)
*
*  where on entry, A*Q has been computed and stored in the array AQ.
*
*  Arguments
*  =========
*
*  N       (input) INTEGER
*          N is the order of the matrix A.
*
*  L, M    (input) INTEGER
*          The first and last column indices of Q.
*
*  Q       (input) DOUBLE PRECISION array, dimension( LDQ,M )
*          On entry, Q contains approximate subspace matrix.
*
*  LDQ     (input) INTEGER
*          The leading dimension of Q, LDQ >= max(1,N).
*
*  AQ      (input) DOUBLE PRECISION array, dimension(LDA, M)
*          On entry, AQ contains the product A*Q.
*
*  LDA     (input) INTEGER
*          The leading dimension of AQ, LDA >= max(1,N).
*
*  T       (input) DOUBLE PRECISION array, dimension( LDT,M )
*          On entry, T contains of Schur form.
*
*  LDT     (input) INTEGER
*          The leading dimension of T, LDT >= max(1,M).
*
*  RSD     (output) DOUBLE PRECISION array, dimension(M)
*          On exit, RSD(L) to RSD(M) contains the 2-norm of the residual
*          vectors.
*
* ======================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, TWO
      PARAMETER          ( ZERO = 0.0D+0, TWO = 2.0D+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            I, J, JNEXT, K, KU
      DOUBLE PRECISION   S
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MIN, SQRT
*     ..
*     .. Executable Statements ..
*
*     Quick return if necessary
*
      IF( L.GT.M )
     $   RETURN
*
*     Compute the residual R = A*Q - Q*T, where A*Q has been 
*     computed and stored in the array AQ. 
*
      DO 30 J = L, M
         KU = MIN( J+1, M )
         IF( J.LT.M )THEN
            IF( T( J+1, J ).EQ.ZERO )
     $         KU = J
         END IF 
*
         RSD( J ) = ZERO
         DO 20 I = 1, N
            S = ZERO
            DO 10 K = 1, KU
               S = S + Q( I, K )*T( K, J )
   10       CONTINUE
            RSD( J ) = RSD( J ) + ( AQ( I, J )-S )*( AQ( I, J )-S )
   20    CONTINUE
*
   30 CONTINUE
*
*     Compute the norms of residual vectors 
*
      JNEXT = L
      DO 40 J = L, M
         IF( J.LT.JNEXT )
     $      GO TO 40
*
         IF( J.EQ.M ) THEN
            RSD( J ) = SQRT( RSD( J ) )
            JNEXT = JNEXT + 1
            GO TO 40
         END IF
*
         IF( T( J+1, J ).EQ.ZERO ) THEN
            RSD( J ) = SQRT( RSD( J ) )
            JNEXT = JNEXT + 1
         ELSE
            RSD( J ) = ( RSD( J ) + RSD( J+1 ) ) / TWO
            RSD( J ) = SQRT( RSD( J ) )
            RSD( J+1 ) = RSD( J )
            JNEXT = JNEXT + 2
         END IF
   40 CONTINUE
*
      RETURN
*
*     End of DRESID
      END
*
******** begin of dgroup.f *********************************************
      SUBROUTINE DGROUP( L, M, WR, WI, RSD, NGRP, CTR, AE, ARSD, GRPTOL)      
*     ..
*     .. Scalar Arguments ..
      INTEGER            L, M, NGRP
      DOUBLE PRECISION   AE, ARSD, CTR, GRPTOL
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   RSD( * ), WI( * ), WR( * )
*     ..
*
*  Purpose
*  =======
*
*  Find a cluster of complex numbers whose real parts are contained
*  in the array WR and imaginary parts are contained in the array WI.
*  These numbers are assumed to be stored in descending order of
*  magnitude. NGRP is determined as the largest integer less than or
*  equal to M-L+1 for which the absolute value W(J) of the number
*  WR(J)+i*WI(J) satisfies
*
*       |W(L) - W(L+J-1)| <= GRPTOL*(W(L)+W(L+J-1))               (1)
*
*  for J = L, L+1, ..., NGRP. CTR is set to (W(L)+W(L+NGRP-1))/2, and
*  AE to the average of the numbers WR(L), ..., WR(L+NGRP-1), and ARSD 
*  to the average of RSD(L), ..., RSD(L+NGRP-1).
*
*  Arguments
*  =========
*
*  L, M    (input) INTEGER
*          The first and the last indices of array WR and WI that will
*          be checked for a possible cluster of complex numbers.
*
*  WR,WI   (input) DOUBLE PRECISION array, dimension(M)
*          On entry, WR and WI contains the real and imaginary parts of        
*          the complex array that is going to be processed, it is 
*          assumed that they are stored in descending order of 
*          magnitude.      
*
*  RSD     (input) DOUBLE PRECISION array, dimension(M)
*          RSD contains the 2-norms of residual vectors computed 
*          by subroutine RESID.       
*
*  NGRP    (output) INTEGER
*          On exit, NGRP is the number of entries that are grouped
*          together to satisfy the condition (1).
*
*  CTR     (output) DOUBLE PRECISION
*          On exit, CTR is the average of the magnitude of the absolute
*          sum of of the the grouped complex numbers.
*
*  AE      (output) DOUBLE PRECISION
*          On exit, AE is the average of the real part of the grouped
*          complex numbers.
*
*  ARSD    (output) DOUBLE PRECISION
*          On exit, ARSD is the average of the residual norm of the
*          grouped complex numbers, see argument RSD.
*
*  GRPTOL  (input) DOUBLE PRECISION
*          On entry, GRPTOL is the threshold to decide whether the two
*          complex numbers are considered to be close and put into the
*          cluster. 
*
* =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, TWO
      PARAMETER          ( ZERO = 0.0D+0, TWO = 2.0D+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            ITYPE, J, L1
      DOUBLE PRECISION   RMOD, RMOD1
*     ..
*     .. External Functions ..
      DOUBLE PRECISION   DLAPY2
      EXTERNAL           DLAPY2
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, DBLE, SQRT
*     ..
*     .. Executable Statements ..
*
      NGRP = 0
      RMOD = DLAPY2( WR( L ), WI( L ) )
      CTR = ZERO
*
*     Find the cluster group 
*
   10 CONTINUE
      L1 = L + NGRP
      IF( L1.GT.M )
     $   GO TO 20
      RMOD1 = DLAPY2( WR( L1 ), WI( L1 ) )
      IF( ABS( RMOD-RMOD1 ).GT.GRPTOL*( RMOD+RMOD1 ) )
     $   GO TO 20
      CTR = ( RMOD+RMOD1 ) / TWO
      ITYPE = 0
      IF( WI( L1 ).NE.ZERO )
     $   ITYPE = 1
      NGRP = NGRP + ITYPE + 1
      GO TO 10
*
*     Compute the averages of the real part (AE) and the residual 
*     norms (ARSD) of the grouped complex numbers.
*
   20 CONTINUE
      AE = ZERO
      ARSD = ZERO
      IF( NGRP.EQ.0 )
     $   GO TO 40
      L1 = L + NGRP - 1
      DO 30 J = L, L1
         AE = AE + WR( J )
         ARSD = ARSD + RSD( J )*RSD( J )
   30 CONTINUE
      AE = AE / DBLE( NGRP )
      ARSD = SQRT( ARSD / DBLE( NGRP ) )
*
   40 CONTINUE
*
      RETURN
*
*     End of DGROUP
      END
******** begin of cond.f **********************************************
      DOUBLE PRECISION FUNCTION DCOND( M, H, LDH, IPVT, MULT, INFO )
*     .. 
*     .. Scalar Arguments ..
      INTEGER            INFO, LDH, M
*     ..
*     .. Array Arguments ..
      INTEGER            IPVT( * )
      DOUBLE PRECISION   H( LDH, * ), MULT( * )
*     ..
*
*  Purpose
*  =======
*
*  Compute the inf-norm condition number of an upper Hessenberg
*  matrix H:
*
*             COND = norm(H)*norm(inv(H))
*
*  In this application, the code uses Gaussian elimination with partial
*  pivoting to compute the inverse explicitly. 
*
*  Arguments
*  =========
*
*  M       (input) INTEGER
*          M is the order of matrix H.
*
*  H       (input/output) DOUBLE PRECISION array, dimension( LDH,M )
*          On entry, H contains an M by M upper Hessenberg matrix.
*          On exit, H contains the inverse of H.
*
*  LDH     (input) INTEGER
*          The leading dimension of the array H, LDH >= max(1,M)
*
*  IPVT    (workspace) INTEGER array, dimension(M)
*
*  MULT    (workspace) DOUBLE PRECISION array, dimension(2*M)
*
*  INFO    (output) INTEGER
*          On exit, if INFO is set to
*             0,  normal return.
*             1,  the matrix is singular, condition number is infinity.
*
*  ====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE, UPPER
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0, UPPER = 1.0D+8 )       
*     ..
*     .. Local Scalars ..
      INTEGER            I, J, JJ, JM1, K
      DOUBLE PRECISION   HIN, HN, S
*     ..
*     .. External Functions ..
      DOUBLE PRECISION   DLANHS, DLANGE
      EXTERNAL           DLANHS, DLANGE
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS
*     ..
*     .. Executable Statements ..
*
      INFO = 0
*
*     Compute inf-norm of the upper Hessenberg matrix H
*
      HN = DLANHS( 'Inf', M, H, LDH, MULT ) 
*
*     Compute the inverse of matrix H
*
*     Step 1: Gaussian elimination with partial pivoting: 
*           (M_{m-1} P_{m-1} .... M_1 P_1 )H = H1 
*     where P_i and M_i are permutation and elementray transformation
*     matrices, resp. (stored in MULT and IPVT), H1 is an upper 
*     triangular 
*
      DO 60 I = 1, M - 1
         IPVT( I ) = 0
         MULT( I ) = ZERO
         IF( H( I+1, I ).EQ.ZERO )
     $      GO TO 60
*
         IF( ABS( H( I+1, I ) ).LE.ABS( H( I, I ) ) )
     $      GO TO 40
*
         IPVT( I ) = 1
         DO 30 J = I, M
            S = H( I, J )
            H( I, J ) = H( I+1, J )
            H( I+1, J ) = S
   30    CONTINUE
*
   40    CONTINUE
*
         MULT( I ) = H( I+1, I ) / H( I, I )
         H( I+1, I ) = ZERO
         DO 50 J = I + 1, M
            H( I+1, J ) = H( I+1, J ) - MULT( I )*H( I, J )
   50    CONTINUE
   60 CONTINUE
*
*     Step 2: compute the inverse of H1 (triangular matrix)
*
      DO 110 J = 1, M
         IF( H( J, J ).EQ.ZERO )THEN 
            DCOND = UPPER
            INFO = 1
            RETURN
         END IF 
*
         H( J, J ) = ONE / H( J, J )
         IF( J.EQ.1 )
     $      GO TO 100
         JM1 = J - 1
         DO 90 I = 1, JM1
            S = ZERO
            DO 80 K = I, JM1
               S = S + H( I, K )*H( K, J )
   80       CONTINUE
            H( I, J ) = -S*H( J, J )
   90    CONTINUE
*
  100    CONTINUE
*
  110 CONTINUE
*
*     Step 3: Compute the inverse of H:
*         inv(H) = inv(H1)*(M_{m-1} P_{m-1} ... M_1 P_1) 
*
      DO 160 JJ = 1, M - 1
         J = M - JJ
         IF( MULT( J ).EQ.ZERO )
     $      GO TO 130
         DO 120 I = 1, M 
            H( I, J ) = H( I, J ) - MULT( J )*H( I, J+1 )
  120    CONTINUE
*
  130    CONTINUE
         IF( IPVT( J ).EQ.0 )
     $      GO TO 150
         DO 140 I = 1, M 
            S = H( I, J )
            H( I, J ) = H( I, J+1 )
            H( I, J+1 ) = S
  140    CONTINUE
*
  150    CONTINUE
  160 CONTINUE
*
*     Compute the inf-norm of the inverse matrix.
*
      HIN = DLANGE( 'Inf', M, M, H, LDH, MULT )
*
*     Compute the condition number 
*
      DCOND = HN*HIN
*
      RETURN
*
*     End of DCOND
      END
*
******** begin of dlaran.f *********************************************
      DOUBLE PRECISION FUNCTION DLARAN( ISEED )
*
*     .. Array Arguments ..
      INTEGER            ISEED( 4 )
*     ..
*
*  Purpose
*  =======
*
*  DLARAN returns a random real number from a uniform (0,1) 
*  distribution. This subroutine is taken from the LAPACK test matrix 
*  suite.
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
*  ===================================================================
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
      INTRINSIC          MOD, DBLE
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
*
******** begin of dlaqr3.f ********************************************
      SUBROUTINE DLAQR3( WANTT, WANTZ, N, ILO, IHI, H, LDH, WR, WI,
     $                   ILOZ, IHIZ, Z, LDZ, WORK, INFO )
*     ..
*     .. Scalar Arguments ..
      LOGICAL            WANTT, WANTZ
      INTEGER            IHI, IHIZ, ILO, ILOZ, INFO, LDH, LDZ, N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   H( LDH, * ), WI( * ), WR( * ), WORK( * ), 
     $                   Z( LDZ, * )
*     ..
*
*  Purpose
*  =======
*
*  This subroutine computes the Schur factorization:
*
*                       H = Z T Z'
*
*  of a real upper Hessenberg matrix H, by dealing with the Hessenberg 
*  submatrix in rows and columns ILO to IHI. T is a matrix in Schur 
*  canonical form, Z is orthogonal, and Z' denotes the transpose of Z. 
*  The eigenvalues of Schur form appear in descending order of magnitude
*  along its diagonal.
*
*  This subroutine is based on the LAPACK auxiliary routines SLAQR3 for 
*  computing the Schur decomposition by the double shift QR algorithm 
*  and the LAPACK subroutine SLAEXC for swapping adjacent diagonal 1 by 
*  1 or 2 by 2 blocks in the upper quasi-triangular matrix T by an 
*  orthogonal similarity transformation.
*
*  Arguments
*  =========
*
*  WANTT   (input) LOGICAL
*          = .TRUE. : the full Schur form T is required;
*          = .FALSE.: only eigenvalues are required.
*
*  WANTZ   (input) LOGICAL
*          = .TRUE. : the matrix of Schur vectors Z is required;
*          = .FALSE.: Schur vectors are not required.
*
*  N       (input) INTEGER
*          The order of the matrix H.  N >= 0.
*
*  ILO     (input) INTEGER
*  IHI     (input) INTEGER
*          It is assumed that H is already upper quasi-triangular in
*          rows and columns IHI+1:N, and that H(ILO,ILO-1) = 0 (unless
*          ILO = 1). SLAQR3 works primarily with the Hessenberg
*          submatrix in rows and columns ILO to IHI, but applies
*          transformations to all of H if WANTT is .TRUE..
*          1 <= ILO <= max(1,IHI); IHI <= N.
*
*  H       (input/output) DOUBLE PRECISION array, dimension (LDH,N)
*          On entry, the upper Hessenberg matrix H.
*          On exit, if WANTT is .TRUE., H is upper quasi-triangular in
*          rows and columns ILO:IHI, with any 2-by-2 diagonal blocks in
*          standard form. If WANTT is .FALSE., the contents of H are
*          unspecified on exit.
*
*  LDH     (input) INTEGER
*          The leading dimension of the array H. LDH >= max(1,N).
*
*  WR      (output) DOUBLE PRECISION array, dimension (N)
*  WI      (output) DOUBLE PRECISION array, dimension (N)
*          The real and imaginary parts, respectively, of the computed
*          eigenvalues ILO to IHI are stored in the corresponding
*          elements of WR and WI. If two eigenvalues are computed as a
*          complex conjugate pair, they are stored in consecutive
*          elements of WR and WI, say the i-th and (i+1)th, with
*          WI(i) > 0 and WI(i+1) < 0. If WANTT is .TRUE., the
*          eigenvalues are stored in the same order as on the diagonal
*          of the Schur form returned in H, with WR(i) = H(i,i), and, 
*          if H(i:i+1,i:i+1) is a 2-by-2 diagonal block,
*          WI(i) = sqrt(H(i+1,i)*H(i,i+1)) and WI(i+1) = -WI(i).
*
*  ILOZ    (input) INTEGER
*  IHIZ    (input) INTEGER
*          Specify the rows of Z to which transformations must be
*          applied if WANTZ is .TRUE..
*          1 <= ILOZ <= ILO; IHI <= IHIZ <= N.
*
*  Z       (input/output) DOUBLE PRECISION array, dimension (LDZ,N)
*          If WANTZ is .TRUE., on entry Z must contain the current
*          matrix Z of transformations accumulated by SHSEQR, and on
*          exit Z has been updated; transformations are applied only to
*          the submatrix Z(ILOZ:IHIZ,ILO:IHI).
*          If WANTZ is .FALSE., Z is not referenced.
*
*  LDZ     (input) INTEGER
*          The leading dimension of the array Z. LDZ >= max(1,N).
*
*  WORK    (workspace) DOUBLE PRECISION array, dimension (N)
*
*  INFO    (output) INTEGER
*          = 0: successful exit
*          = 1,...,N: SLAQR3 failed to compute all the eigenvalues ILO 
*               to IHI in a total of 30*(IHI-ILO+1) iterations; 
*               if INFO = i, elements i+1:ihi of WR and WI contain those
*               eigenvalues which have been successfully computed.
*          = N+i: SLAEXC failed to swap the computed eigenvalues in the
*               desired order for the ith computed eigenvalues 
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
      DOUBLE PRECISION   DAT1, DAT2
      PARAMETER          ( DAT1 = 0.75D+0, DAT2 = -0.4375D+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            I, I1, I2, ITN, ITS, J, K, L, M, NH, NR, NZ,
     $                   IC, IN, IERR
      DOUBLE PRECISION   CS, H00, H10, H11, H12, H21, H22, H33, H33S,
     $                   H43H34, H44, H44S, OVFL, S, SMLNUM, SN, SUM,
     $                   T1, T2, T3, TEMP, TST1, ULP, UNFL, V1, V2, V3,
     $                   WR1
*     ..
*     .. Local Arrays ..
      DOUBLE PRECISION   V( 3 )
*     ..
*     .. External Functions ..
      DOUBLE PRECISION   DLAMCH, DLANHS, DLAPY2
      EXTERNAL           DLAMCH, DLANHS, DLAPY2
*     ..
*     .. External Subroutines ..
      EXTERNAL           DCOPY, DLABAD, DLANV2, DLARFG, DROT, DLAEXC
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, MIN
*     ..
*     .. Executable Statements ..
*
      INFO = 0
*
*     Quick return if possible
*
      IF( N.EQ.0 )
     $   RETURN
      IF( ILO.EQ.IHI ) THEN
         WR( ILO ) = H( ILO, ILO )
         WI( ILO ) = ZERO
         RETURN
      END IF
*
      NH = IHI - ILO + 1
      NZ = IHIZ - ILOZ + 1
*
*     Set machine-dependent constants for the stopping criterion.
*     If norm(H) <= sqrt(OVFL), overflow should not occur.
*
      UNFL = DLAMCH( 'Safe minimum' )
      OVFL = ONE / UNFL
      CALL DLABAD( UNFL, OVFL )
      ULP = DLAMCH( 'Precision' )
      SMLNUM = UNFL*( NH / ULP )
*
*     I1 and I2 are the indices of the first row and last column of H
*     to which transformations must be applied. If eigenvalues only are
*     being computed, I1 and I2 are set inside the main loop.
*
      IF( WANTT ) THEN
         I1 = 1
         I2 = N
      END IF
*
*     ITN is the total number of QR iterations allowed.
*
      ITN = 30*NH
*
*     The main loop begins here. I is the loop index and decreases from       
*     IHI to ILO in steps of 1 or 2. Each iteration of the loop works
*     with the active submatrix in rows and columns L to I.
*     Eigenvalues I+1 to IHI have already converged. Either L = ILO or
*     H(L,L-1) is negligible so that the matrix splits.
*
      I = IHI
   10 CONTINUE
      L = ILO
      IF( I.LT.ILO )
     $   GO TO 150
*
*     Perform QR iterations on rows and columns ILO to I until a
*     submatrix of order 1 or 2 splits off at the bottom because a
*     subdiagonal element has become negligible.
*
      DO 130 ITS = 0, ITN
*
*        Look for a single small subdiagonal element.
*
         DO 20 K = I, L + 1, -1
            TST1 = ABS( H( K-1, K-1 ) ) + ABS( H( K, K ) )
            IF( TST1.EQ.ZERO )
     $         TST1 = DLANHS( '1', I-L+1, H( L, L ), LDH, WORK )
            IF( ABS( H( K, K-1 ) ).LE.MAX( ULP*TST1, SMLNUM ) )
     $         GO TO 30
   20    CONTINUE
   30    CONTINUE
         L = K
         IF( L.GT.ILO ) THEN
*
*           H(L,L-1) is negligible
*
            H( L, L-1 ) = ZERO
         END IF
*
*        Exit from loop if a submatrix of order 1 or 2 has split off.
*
         IF( L.GE.I-1 )
     $      GO TO 140
*
*        Now the active submatrix is in rows and columns L to I. If
*        eigenvalues only are being computed, only the active submatrix
*        need be transformed.
*
         IF( .NOT.WANTT ) THEN
            I1 = L
            I2 = I
         END IF
*
         IF( ITS.EQ.10 .OR. ITS.EQ.20 ) THEN
*
*           Exceptional shift.
*
            S = ABS( H( I, I-1 ) ) + ABS( H( I-1, I-2 ) )
            H44 = DAT1*S
            H33 = H44
            H43H34 = DAT2*S*S
         ELSE
*
*           Prepare to use Wilkinson's double shift
*
            H44 = H( I, I )
            H33 = H( I-1, I-1 )
            H43H34 = H( I, I-1 )*H( I-1, I )
         END IF
*
*        Look for two consecutive small subdiagonal elements.
*
         DO 40 M = I - 2, L, -1
*
*           Determine the effect of starting the double-shift QR
*           iteration at row M, and see if this would make H(M,M-1)
*           negligible.
*
            H11 = H( M, M )
            H22 = H( M+1, M+1 )
            H21 = H( M+1, M )
            H12 = H( M, M+1 )
            H44S = H44 - H11
            H33S = H33 - H11
            V1 = ( H33S*H44S-H43H34 ) / H21 + H12
            V2 = H22 - H11 - H33S - H44S
            V3 = H( M+2, M+1 )
            S = ABS( V1 ) + ABS( V2 ) + ABS( V3 )
            V1 = V1 / S
            V2 = V2 / S
            V3 = V3 / S
            V( 1 ) = V1
            V( 2 ) = V2
            V( 3 ) = V3
            IF( M.EQ.L )
     $         GO TO 50
            H00 = H( M-1, M-1 )
            H10 = H( M, M-1 )
            TST1 = ABS( V1 )*( ABS( H00 )+ABS( H11 )+ABS( H22 ) )
            IF( ABS( H10 )*( ABS( V2 )+ABS( V3 ) ).LE.ULP*TST1 )
     $         GO TO 50
   40    CONTINUE
   50    CONTINUE
*
*        Double-shift QR step
*
         DO 120 K = M, I - 1
*
*           The first iteration of this loop determines a reflection G
*           from the vector V and applies it from left and right to H,
*           thus creating a nonzero bulge below the subdiagonal.
*
*           Each subsequent iteration determines a reflection G to
*           restore the Hessenberg form in the (K-1)th column, and thus
*           chases the bulge one step toward the bottom of the active
*           submatrix. NR is the order of G.
*
            NR = MIN( 3, I-K+1 )
            IF( K.GT.M )
     $         CALL DCOPY( NR, H( K, K-1 ), 1, V, 1 )
            CALL DLARFG( NR, V( 1 ), V( 2 ), 1, T1 )
            IF( K.GT.M ) THEN
               H( K, K-1 ) = V( 1 )
               H( K+1, K-1 ) = ZERO
               IF( K.LT.I-1 )
     $            H( K+2, K-1 ) = ZERO
            ELSE IF( M.GT.L ) THEN
               H( K, K-1 ) = -H( K, K-1 )
            END IF
            V2 = V( 2 )
            T2 = T1*V2
            IF( NR.EQ.3 ) THEN
               V3 = V( 3 )
               T3 = T1*V3
*
*              Apply G from the left to transform the rows of the matrix
*              in columns K to I2.
*
               DO 60 J = K, I2
                  SUM = H( K, J ) + V2*H( K+1, J ) + V3*H( K+2, J )
                  H( K, J ) = H( K, J ) - SUM*T1
                  H( K+1, J ) = H( K+1, J ) - SUM*T2
                  H( K+2, J ) = H( K+2, J ) - SUM*T3
   60          CONTINUE
*
*              Apply G from the right to transform the columns of the
*              matrix in rows I1 to min(K+3,I).
*
               DO 70 J = I1, MIN( K+3, I )
                  SUM = H( J, K ) + V2*H( J, K+1 ) + V3*H( J, K+2 )
                  H( J, K ) = H( J, K ) - SUM*T1
                  H( J, K+1 ) = H( J, K+1 ) - SUM*T2
                  H( J, K+2 ) = H( J, K+2 ) - SUM*T3
   70          CONTINUE
*
               IF( WANTZ ) THEN
*
*                 Accumulate transformations in the matrix Z
*
                  DO 80 J = ILOZ, IHIZ
                     SUM = Z( J, K ) + V2*Z( J, K+1 ) + V3*Z( J, K+2 )
                     Z( J, K ) = Z( J, K ) - SUM*T1
                     Z( J, K+1 ) = Z( J, K+1 ) - SUM*T2
                     Z( J, K+2 ) = Z( J, K+2 ) - SUM*T3
   80             CONTINUE
               END IF
            ELSE IF( NR.EQ.2 ) THEN
*
*              Apply G from the left to transform the rows of the matrix
*              in columns K to I2.
*
               DO 90 J = K, I2
                  SUM = H( K, J ) + V2*H( K+1, J )
                  H( K, J ) = H( K, J ) - SUM*T1
                  H( K+1, J ) = H( K+1, J ) - SUM*T2
   90          CONTINUE
*
*              Apply G from the right to transform the columns of the
*              matrix in rows I1 to min(K+3,I).
*
               DO 100 J = I1, I
                  SUM = H( J, K ) + V2*H( J, K+1 )
                  H( J, K ) = H( J, K ) - SUM*T1
                  H( J, K+1 ) = H( J, K+1 ) - SUM*T2
  100          CONTINUE
*
               IF( WANTZ ) THEN
*
*                 Accumulate transformations in the matrix Z
*
                  DO 110 J = ILOZ, IHIZ
                     SUM = Z( J, K ) + V2*Z( J, K+1 )
                     Z( J, K ) = Z( J, K ) - SUM*T1
                     Z( J, K+1 ) = Z( J, K+1 ) - SUM*T2
  110             CONTINUE
               END IF
            END IF
  120    CONTINUE
*
  130 CONTINUE
*
*     Failure to converge in remaining number of iterations
*
      INFO = I
      RETURN
*
  140 CONTINUE
*
      IF( L.EQ.I ) THEN
*
*        H(I,I-1) is negligible: one eigenvalue has converged.
*        loop to position 1 x 1 block
*
         IF( I.EQ.IHI )
     $       GO TO 270
*
         WR1 = ABS( H( I, I ) )
         IC = I 
         IN = I + 1
*
  250    CONTINUE
         IF( IN.EQ.IHI .OR. H( IN+1, IN ).EQ.ZERO )THEN
            IF( ABS( H( IN, IN ) ).GT.WR1 )THEN
                CALL DLAEXC( WANTZ, N, H, LDH, Z, LDZ, IC, 1, 1, WORK,
     $                       IERR )
                IF( IERR.NE.0 )THEN
                   INFO = N + I
                   RETURN
                END IF 
                IC = IC + 1
                IN = IC + 1
            ELSE
                GO TO 270
            END IF
         ELSE
            TEMP = SQRT(ABS(H( IN+1, IN )))*SQRT(ABS( H( IN, IN+1 )))
            IF( DLAPY2( H( IN, IN ), TEMP ).GT.WR1 )THEN
               CALL DLAEXC( WANTZ, N, H, LDH, Z, LDZ, IC, 1, 2, WORK,
     $                      IERR )
               IF( IERR.NE.0 )THEN
                   INFO = N + I
                   RETURN
               END IF 
               IC = IC + 2
               IN = IC + 1
            ELSE
               GO TO 270
            END IF
         END IF
*
         IF( IN.GT.IHI )
     $      GO TO 270
*
         GO TO 250
*
      ELSE IF( L.EQ.I-1 ) THEN
*
*        H(I-1,I-2) is negligible: a pair of eigenvalues have converged.
*
*        Transform the 2-by-2 submatrix to standard Schur form,
*        and compute and store the eigenvalues.
*
         CALL DLANV2( H( I-1, I-1 ), H( I-1, I ), H( I, I-1 ),
     $                H( I, I ), WR( I-1 ), WI( I-1 ), WR( I ), WI( I ),
     $                CS, SN )
*
         IF( WANTT ) THEN
*
*           Apply the transformation to the rest of H.
*
            IF( I2.GT.I )
     $         CALL DROT( I2-I, H( I-1, I+1 ), LDH, H( I, I+1 ), LDH,
     $                    CS, SN )
            CALL DROT( I-I1-1, H( I1, I-1 ), 1, H( I1, I ), 1, CS, SN )
         END IF
         IF( WANTZ ) THEN
*
*           Apply the transformation to Z.
*
            CALL DROT( NZ, Z( ILOZ, I-1 ), 1, Z( ILOZ, I ), 1, CS, SN )
         END IF
*        
         IF( H( I, I-1 ).EQ.ZERO ) THEN
*
*          two real eigenvalues, loop to position the real eigenvalues 
*
           IF( I.EQ.IHI )
     $        GO TO 240
*
           WR1 = ABS( H( I, I ) )
           IC = I 
           IN = I + 1
*
  230      CONTINUE 
           IF( IN.EQ.IHI .OR. H( IN+1, IN ).EQ.ZERO ) THEN
               IF( ABS( H( IN, IN ) ).GT.WR1 )THEN 
                  CALL DLAEXC( WANTZ, N, H, LDH, Z, LDZ, IC, 1, 1, 
     $                         WORK, IERR )
                  IF( IERR.NE.0 )THEN
                     INFO = N + I
                     RETURN
                  END IF 
                  IC = IC+1 
                  IN = IC+1
               ELSE 
                  GO TO 240 
               END IF 
            ELSE
               TEMP = SQRT(ABS(H( IN+1, IN )))*SQRT(ABS(H( IN, IN+1 )))
               IF( DLAPY2( H( IN, IN ), TEMP ).GT.WR1 )THEN
                  CALL DLAEXC( WANTZ, N, H, LDH, Z, LDZ, IC, 1, 2,  
     $                         WORK, IERR )
                  IF( IERR.NE.0 )THEN
                     INFO = N + I
                     RETURN
                  END IF 
                  IC = IC+2  
                  IN = IC+1 
               ELSE 
                  GO TO 240 
               END IF 
            END IF
            IF( IN.GT.IHI )
     $          GO TO 240
*
            GO TO 230
*
*           loop to position the first real eigenvalue
*
  240       CONTINUE
            WR1 = ABS( H( I-1, I-1 ) )
            IC = I - 1 
            IN = IC + 1 
*
  245       CONTINUE
            IF( IN.EQ.IHI .OR. H( IN+1, IN ).EQ.ZERO ) THEN
               IF( ABS( H( IN, IN ) ).GT.WR1 )THEN
                  CALL DLAEXC( WANTZ, N, H, LDH, Z, LDZ, IC, 1, 1, 
     $                         WORK, IERR )
                  IF( IERR.NE.0 )THEN
                     INFO = N + I - 1
                     RETURN
                  END IF 
                  IC = IC+1
                  IN = IC+1
               ELSE
                  GO TO 270
               END IF
            ELSE
               TEMP = SQRT(ABS(H( IN+1, IN )))*SQRT(ABS(H( IN, IN+1 )))
               IF( DLAPY2( H( IN, IN ), TEMP ).GT.WR1 )THEN
                  CALL DLAEXC( WANTZ, N, H, LDH, Z, LDZ, IC, 1, 2, 
     $                         WORK, IERR )
                  IF( IERR.NE.0 )THEN
                     INFO = N + I - 1 
                     RETURN
                  END IF 
                  IC = IC+2 
                  IN = IC+1 
               ELSE
                  GO TO 270
               END IF
            END IF
            IF( IN.GT.IHI )
     $          GO TO 270 
*
            GO TO 245
*
         ELSE
*
*           Complex conjugate eigenvalues, loop to position 2 x 2 block
*
            IF( I.EQ.IHI )
     $         GO TO 270
*
            WR1 = DLAPY2( H( I-1, I-1 ), 
     $                SQRT(ABS(H(I, I-1 )))*SQRT(ABS(H( I-1, I ))) )
            IC = I-1
            IN = IC + 2 
*
  260       CONTINUE 
            IF( IN.EQ.IHI .OR. H( IN+1, IN ).EQ.ZERO ) THEN
               IF( ABS( H( IN, IN ) ).GT.WR1 )THEN 
                  CALL DLAEXC( WANTZ, N, H, LDH, Z, LDZ, IC, 2, 1,  
     $                         WORK, IERR )
                  IF( IERR.NE.0 )THEN
                     INFO = N + I - 1
                     RETURN
                  END IF 
                  IC = IC + 1  
                  IN = IC + 2
               ELSE
                  GO TO 270 
               END IF 
            ELSE
               TEMP = SQRT(ABS(H( IN+1, IN )))*SQRT(ABS(H( IN, IN+1 ))) 
               IF( DLAPY2( H( IN, IN ), TEMP) .GT.WR1 )THEN
                  CALL DLAEXC( WANTZ, N, H, LDH, Z, LDZ, IC, 2, 2,  
     $                          WORK, IERR )
                  IF( IERR.NE.0 )THEN
                     INFO = N + I - 1
                     RETURN
                  END IF 
                  IC = IC + 2  
                  IN = IC + 2
               ELSE
                  GO TO 270 
               END IF 
            END IF
            IF( IN.GT.IHI )
     $         GO TO 270 
*
            GO TO 260 
*
         END IF
      END IF
*
  270 CONTINUE 
*
*     Decrement number of remaining iterations, and return to start of
*     the main loop with new value of I.
*
      ITN = ITN - ITS
      I = L - 1
      GO TO 10
*
  150 CONTINUE
*
*     All eigenvalues have been found and ordered, set the final
*     eigenvalues to WR and WI
*
      IN = ILO
      DO 300 I = ILO, IHI
         IF( I.LT.IN )
     $      GO TO 300
         IF( I.EQ.IHI ) THEN
            WR( I ) = H( I, I )
            WI( I ) = ZERO
            GO TO 300
         END IF
*
         IF( H( I+1, I ).EQ.ZERO ) THEN
            WR( I ) = H( I, I )
            WI( I ) = ZERO
            IN = I + 1
         ELSE
            WR( I ) = H( I, I )
            WR( I+1 ) = H( I, I )
            WI( I ) = SQRT( ABS( H( I+1, I ) ) )*
     $                SQRT( ABS( H( I, I+1 ) ) )
            WI( I+1 ) = -WI( I )
            IN = I + 2
         END IF
  300 CONTINUE
*
      RETURN
*
*     End of DLAQR3
*
      END
