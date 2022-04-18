#include "fintrf.h"
C SKEWHAMIL2FEIG.F - Gateway function for computing the eigenvalues of
C                    a real skew-Hamiltonian/skew-Hamiltonian pencil 
C                    using the routine DGHFST.
C
C Matlab call:
C   [ALPHAR,ALPHAI,BETA(,Zo,Bo,Fo(,Q,U1,U2))] =
C           skewHamil2feig(Z,B,FG(,job,compq,compu(,U01,U02)))
C
C   [ALPHAR,ALPHAI,BETA]          = skewHamil2feig(Z,B,FG)
C   [ALPHAR,ALPHAI,BETA,Zo,Bo,Fo] = skewHamil2feig(Z,B,FG,1)
C
C Purpose:
C   To compute the eigenvalues of a real n-by-n skew-Hamiltonian/
C   skew-Hamiltonian pencil aS - bT, with n = 2m, and
C
C                           (  B  F  )            (  0  I  )
C     S = J Z' J' Z and T = (        ), where J = (        ).        (1)
C                           (  G  B' )            ( -I  0  )
C
C   Optionally, if job = 1, aS - bT will be transformed to
C   structured Schur form: an orthogonal transformation matrix Q and
C   an orthogonal symplectic transformation matrix U are computed,
C   such that
C
C              (  Z11  Z12  )
C     U' Z Q = (            ) = Zo, and
C              (   0   Z22' )
C                                                                    (2)
C                   (  Bo  Fo  )
C     J Q' J' T Q = (          ),
C                   (   0  Bo' )
C
C   where Z11 and Z22 are upper triangular and Bo is upper quasi-
C   triangular. The notation M' denotes the transpose of the matrix M.
C   Optionally, if compq <> 0, the orthogonal transformation matrix Q
C   will be computed.
C   Optionally, if compu > 0, the orthogonal symplectic transformation
C   matrix
C
C         (  U1  U2  )
C     U = (          )
C         ( -U2  U1  )
C
C   will be computed.
C
C Input parameters:
C   Z      - the n-by-n matrix Z.
C   B      - the m-by-m matrix B.
C   FG     - an  m-by-(m+1) matrix containing the strict triangles of
C            the skew-symmetric matrices F and G, as follows:
C            the leading m-by-m strictly lower triangular part contains
C            the strictly lower triangle of the matrix G, and the
C            m-by-m strictly upper triangular part of the submatrix in
C            the columns 2 to m+1 contains the strictly upper triangle
C            of the matrix F of T in (1).
C            So, if i > j, then G(i,j) = -G(j,i) is stored in FG(i,j)
C            and F(j,i) = -F(i,j) is stored in FG(j,i+1).
C            The entries on the diagonal and the first superdiagonal of
C            FG need not be set, but are assumed to be zero.
C            FG is an empty matrix if m = 0.
C   job    - (optional) scalar indicating the computation to be
C            performed, as follows:
C            = 0 :  compute the eigenvalues only (default);
C            = 1 :  compute the eigenvalues and the matrices of the
C                   transformed pencil in (2).
C   compq  - (optional) scalar indicating whether the orthogonal
C            transformation matrix Q is returned, as follows:
C            = 0 :  Q is not required (default);
C            = 1 :  on exit, Q contains the orthogonal matrix Q;
C            =-1 :  on exit, for job = 0, Q contains the orthogonal
C                   matrix Q1 which reduced Z to 
C
C                            (  Z11  Z12  )
C                     Z*Q1 = (            ),
C                            (   0   Z22  )
C
C                   where Z11 and Z22' are upper triangular (the first
C                   step of the algorithm).
C   compu  - (optional) scalar indicating whether the orthogonal
C            symplectic transformation matrix U is returned, as follows:
C            = 0 :  U is not required (default);
C            = 1 :  on exit, U1 and U2 contain the submatrices 
C                   of the orthogonal symplectic matrix U;
C            = 2 :  the orthogonal transformations are accumulated
C                   into U;
C                   on input, U01 and U02 must contain the corresponding
C                   submatrices of an orthogonal symplectic matrix U0;
C                   on exit, U1 and U2 contain the updated submatrices
C                   U1 and U2 of the matrix product U0*U, with U
C                   returned for compu = 1.
C   U01,   - if compu = 2, the m-by-m submatrices U01 and U02 of U0,
C   U02      respectively.
C
C Output parameters:
C   ALPHAR,- the m-vectors of real parts and imaginary parts,
C   ALPHAI   respectively, of each scalar alpha defining an eigenvalue
C            of the pencil aS - bT.
C            If ALPHAI(j) is zero, then the j-th eigenvalue is real.
C   BETA     the m-vector of the scalars beta that define the eigenvalues
C            of the pencil aS - bT.
C            Together, the quantities alpha = (ALPHAR(j),ALPHAI(j)) and
C            beta = BETA(j) represent the j-th eigenvalue of the pencil
C            aS - bT, in the form lambda = alpha/beta. Since lambda may
C            overflow, the ratios should not, in general, be computed.
C            Due to the skew-Hamiltonian/skew-Hamiltonian structure of
C            the pencil, every eigenvalue occurs twice and thus it has
C            only to be saved once in ALPHAR, ALPHAI and BETA.
C   Zo     - if job = 1, the computed n-by-n matrix Zo in (2).
C   Bo     - if job = 1, the computed m-by-m strictly upper triangular
C            part of the skew-symmetric matrix Bo in (2).
C   Fo     - if job = 1, the computed m-by-m strictly upper triangular
C            part of the skew-symmetric matrix Fo in (2).
C   Q      - if compq <> 0, an n-by-n matrix containing the computed
C            orthogonal matrix Q or Q1.
C   U1     - if compu > 0, an m-by-m matrix containing the computed
C            matrix U1.
C   U2     - if compu > 0, an m-by-m matrix containing the computed
C            matrix U2.
C            If job = 0, Zo, Bo, Fo, Q, U1, and U2 contain the
C            corresponding matrices just before the application of the
C            periodic QZ algorithm.
C
C Contributors:
C   V. Sima, Research Institute for Informatics, Bucharest, Oct. 2009.
C
C Revisions:
C   V. Sima, Dec. 2010, July 2012, Nov. 2012, Jul. 2013.
C   M. Voigt, Jul. 2013.
C
C     ******************************************************************
C
      SUBROUTINE MEXFUNCTION( NLHS, PLHS, NRHS, PRHS )
C
C     .. Mex-file interface parameters ..
      mwPointer         PLHS( * ), PRHS( * )
      INTEGER*4         NLHS, NRHS
C
C     .. Mex-file integer functions ..
      mwPointer         mxCreateDoubleMatrix, mxGetPr
      INTEGER*4         mxGetM, mxGetN, mxIsNumeric, mxIsComplex
C
C     .. Scalar parameters used by subroutines ..
      CHARACTER         COMPQ, COMPU, JOB
      INTEGER           INFO, LDB, LDFG, LDQ, LDU1, LDU2, LDWORK, LDZ,
     $                  LIWORK, N
C
C     .. Allocatable arrays ..
C     !Fortran 90/95 (Fixed dimensions should be used with Fortran 77.)
      INTEGER,          ALLOCATABLE :: IWORK( : )
      DOUBLE PRECISION, ALLOCATABLE :: DW( : )
      mwPointer         ALPHAI, ALPHAR, B, BETA, DWORK, FG, Q, U1, U2, Z
C
C     .. Local variables and constant dimension arrays ..
      CHARACTER*120     TEXT
      INTEGER           ICMPQ, ICMPU, IDW, IJOB, IP, M, M1, N1
      DOUBLE PRECISION  DUM( 1 ), TEMP
C
C     .. External Subroutines ..
      EXTERNAL          DLACPY, DGHFST
C
C     ..Intrinsic Functions..
      INTRINSIC         INT, MAX, MIN, MOD
C
C     Check for proper number of arguments.
C 
      IF( NRHS.LT.3 ) THEN
         CALL mexErrMsgTxt
     $        ( 'SKEWHAMIL2FEIG requires at least 3 input arguments.' )
      ELSE IF ( NLHS.LT.3 ) THEN
         CALL mexErrMsgTxt
     $       ( 'SKEWHAMIL2FEIG requires at least 3 output arguments.' )
      END IF
C
C   Z(nxn), B(mxm), FG(mx(m+1))(, job, compq, compu(, U01(mxm), U02(mxm))).
C
      N  = mxGetM( PRHS( 1 ) )
      N1 = mxGetN( PRHS( 1 ) )
      IF ( mxIsNumeric( PRHS(1) ).EQ.0 .OR.
     $     mxIsComplex( PRHS(1) ).EQ.1 ) THEN
         CALL mexErrMsgTxt( 'Z must be a real matrix' )
      END IF
      IF ( N.NE.N1 ) THEN
         CALL mexErrMsgTxt( 'Z must be a square matrix' )
      END IF
      IF ( MOD( N, 2 ).NE.0 ) THEN
         CALL mexErrMsgTxt( 'N must be even' )
      END IF
      M = N/2
C
      M1 = mxGetM( PRHS(2) )
      N1 = mxGetN( PRHS(2) )
      IF ( mxIsNumeric( PRHS(2) ).EQ.0 .OR.
     $     mxIsComplex( PRHS(2) ).EQ.1 ) THEN
         CALL mexErrMsgTxt( 'B must be a real matrix' )
      END IF
      IF ( M1.NE.M ) THEN
         CALL mexErrMsgTxt( 'B must have half as many rows as Z' )
      END IF
      IF ( N1.NE.M ) THEN
         CALL mexErrMsgTxt( 'B must be a square matrix' )
      END IF
C
      M1 = mxGetM( PRHS(3) )
      N1 = mxGetN( PRHS(3) )
      IF ( mxIsNumeric( PRHS(3) ).EQ.0 .OR.
     $     mxIsComplex( PRHS(3) ).EQ.1 ) THEN
         CALL mexErrMsgTxt( 'FG must be a real matrix' )
      END IF
      IF ( M1.NE.M ) THEN
         CALL mexErrMsgTxt( 'FG must have half as many rows as Z' )
      END IF
      IF ( N1.NE.M+1 ) THEN
         CALL mexErrMsgTxt( 'FG must have one more columns than rows' )
      END IF
C
      IP = 4
      IF ( NRHS.GE.IP ) THEN
         IF ( mxGetM( PRHS(IP) ).NE.1 .OR.
     $        mxGetN( PRHS(IP) ).NE.1 ) THEN
            CALL mexErrMsgTxt( 'JOB must be a scalar' )
         END IF
         IF ( mxIsNumeric( PRHS(IP) ).EQ.0 .OR.
     $        mxIsComplex( PRHS(IP) ).EQ.1 ) THEN
            CALL mexErrMsgTxt( 'JOB must be an integer scalar' )
         END IF
         CALL mxCopyPtrToReal8( mxGetPr( PRHS(IP) ), TEMP, 1 )
         IJOB = TEMP
         IF ( IJOB.LT.0 .OR. IJOB.GT.1 ) THEN
            CALL mexErrMsgTxt
     $           ( 'JOB has 0 or 1 the only admissible values' )
         END IF
         IP = IP + 1
      ELSE
         IJOB = 0
      END IF
C
      IF ( IJOB.EQ.0 ) THEN
         JOB = 'E'
      ELSE
         JOB = 'T'
      END IF
C
      IF ( NRHS.GE.IP ) THEN
         IF ( mxGetM( PRHS(IP) ).NE.1 .OR.
     $        mxGetN( PRHS(IP) ).NE.1 ) THEN
            CALL mexErrMsgTxt( 'COMPQ must be a scalar' )
         END IF
         IF ( mxIsNumeric( PRHS(IP) ).EQ.0 .OR.
     $        mxIsComplex( PRHS(IP) ).EQ.1 ) THEN
            CALL mexErrMsgTxt( 'COMPQ must be an integer scalar' )
         END IF
         CALL mxCopyPtrToReal8( mxGetPr( PRHS(IP) ), TEMP, 1 )
         ICMPQ = TEMP
         IF ( ICMPQ.LT.-1 .OR. ICMPQ.GT.1 ) THEN
            CALL mexErrMsgTxt
     $          ( 'COMPQ has -1, 0, or 1 the only admissible values' )
         END IF
         IP = IP + 1
C
      ELSE
         ICMPQ = 0
      END IF
C
      IF ( ICMPQ.LE.0 ) THEN
         COMPQ = 'N'
      ELSE
         COMPQ = 'I'
      END IF
C
      IF ( NRHS.GE.IP ) THEN
         IF ( mxGetM( PRHS(IP) ).NE.1 .OR.
     $        mxGetN( PRHS(IP) ).NE.1 ) THEN
            CALL mexErrMsgTxt( 'COMPU must be a scalar' )
         END IF
         IF ( mxIsNumeric( PRHS(IP) ).EQ.0 .OR.
     $        mxIsComplex( PRHS(IP) ).EQ.1 ) THEN
            CALL mexErrMsgTxt( 'COMPU must be an integer scalar' )
         END IF
         CALL mxCopyPtrToReal8( mxGetPr( PRHS(IP) ), TEMP, 1 )
         ICMPU = TEMP
         IF ( ICMPU.LT.0 .OR. ICMPU.GT.2 ) THEN
            CALL mexErrMsgTxt
     $           ( 'COMPU has 0, 1, or 2 the only admissible values' )
         END IF
C
      ELSE
         ICMPU = 0
      END IF
C
      IF ( ICMPU.EQ.0 ) THEN
         COMPU = 'N'
      ELSE IF ( ICMPU.EQ.1 ) THEN
         COMPU = 'I'
      ELSE
         COMPU = 'U'
      END IF
C
C Determine dimensions of the arrays and the workspace.
C
      LDZ  = MAX( 1, N )
      LDB  = MAX( 1, M )
      LDFG = LDB
      LDQ  = LDZ
      IF ( ICMPU.EQ.0 ) THEN
         LDU1 = 1
         LDU2 = 1
         N1   = 0
      ELSE
         LDU1 = LDB
         LDU2 = LDB
         N1   = M
      END IF
      LIWORK = N + 9
C
C Allocate variable dimension local arrays.
C !Fortran 90/95
C
      ALLOCATE( IWORK( LIWORK ) )
      ALPHAI = mxCalloc(      M,     8 )
      ALPHAR = mxCalloc(      M,     8 )
      B      = mxCalloc(  LDB*M,     8 )
      BETA   = mxCalloc(      M,     8 )
      FG     = mxCalloc( LDFG*(M+1), 8 )
      Q      = mxCalloc(  LDQ*N,     8 )
      U1     = mxCalloc( LDU1*N1,    8 )
      U2     = mxCalloc( LDU2*N1,    8 )
      Z      = mxCalloc(  LDZ*N,     8 )
C
      CALL DGHFST( JOB, COMPQ, COMPU, N, Z, LDZ, B, LDB, FG, LDFG, Q,
     $             LDQ, U1, LDU1, U2, LDU2, ALPHAR, ALPHAI, BETA,
     $             IWORK, LIWORK, DUM, -1, INFO )
      LDWORK = INT( DUM( 1 ) )
      DWORK  = mxCalloc( LDWORK, 8 )
C
C Copy inputs from MATLAB workspace to locally allocated arrays.
C
      CALL mxCopyPtrToReal8( mxGetPr( PRHS( 1 ) ), %VAL( Z ),  N*N )
      CALL mxCopyPtrToReal8( mxGetPr( PRHS( 2 ) ), %VAL( B ),  M*M )
      CALL mxCopyPtrToReal8( mxGetPr( PRHS( 3 ) ), %VAL( FG ),
     $                       M*( M + 1 ) )
      IP = 7
      IF ( NRHS.GE.IP .AND. ICMPU.EQ.2 ) THEN
         CALL mxCopyPtrToReal8( mxGetPr( PRHS( IP ) ), %VAL( U1 ),
     $                          M*M )
         IP = IP + 1
         CALL mxCopyPtrToReal8( mxGetPr( PRHS( IP ) ), %VAL( U2 ),
     $                          M*M )
      END IF
C
C Do the actual computations.
C
      CALL DGHFST( JOB, COMPQ, COMPU, N, %VAL( Z ), LDZ, %VAL( B ), LDB,
     $             %VAL( FG ), LDFG, %VAL( Q ), LDQ, %VAL( U1 ), LDU1,
     $             %VAL( U2 ), LDU2, %VAL( ALPHAR ), %VAL( ALPHAI ),
     $             %VAL( BETA ), IWORK, LIWORK, %VAL( DWORK ), LDWORK,
     $             INFO )
C
C Copy output to MATLAB workspace.
C
      IF ( INFO.EQ.0 ) THEN
         PLHS( 1 ) = mxCreateDoubleMatrix( M, 1, 0 )
         CALL mxCopyReal8ToPtr( %VAL( ALPHAR ),
     $                          mxGetPr( PLHS( 1 ) ), M )
C
         PLHS( 2 ) = mxCreateDoubleMatrix( M, 1, 0 )
         CALL mxCopyReal8ToPtr( %VAL( ALPHAI ),
     $                          mxGetPr( PLHS( 2 ) ), M )
C
         PLHS( 3 ) = mxCreateDoubleMatrix( M, 1, 0 )
         CALL mxCopyReal8ToPtr( %VAL( BETA ), mxGetPr( PLHS( 3 ) ),
     $                          M )
C
         IP = 4
         IF ( NLHS.GE.IP ) THEN
            PLHS( IP ) = mxCreateDoubleMatrix( N, N, 0 )
            CALL mxCopyReal8ToPtr( %VAL( Z ), mxGetPr( PLHS( IP ) ),
     $                             N*N )
            IP = IP + 1
         END IF
C
         IF ( NLHS.GE.IP ) THEN
            PLHS( IP ) = mxCreateDoubleMatrix( M, M, 0 )
            CALL mxCopyReal8ToPtr( %VAL( B ), mxGetPr( PLHS( IP ) ),
     $                             M*M )
            IP = IP + 1
         END IF
C
         IF ( NLHS.GE.IP ) THEN
            ALLOCATE( DW( M*(M+1) ) )
            CALL DLACPY( 'Full', M, M+1, %VAL( FG ), LDFG, DW, LDFG )
            PLHS( IP ) = mxCreateDoubleMatrix( M, M, 0 )
            IDW = MIN( M*( M+1 ), M+1 )
            CALL mxCopyReal8ToPtr( DW( IDW ), mxGetPr( PLHS( IP ) ),
     $                             M*M )
            DEALLOCATE( DW )
            IP = IP + 1
         END IF
C
         IF ( NLHS.GE.IP .AND. ICMPQ.NE.0 ) THEN
            PLHS( IP ) = mxCreateDoubleMatrix( N, N, 0 )
            CALL mxCopyReal8ToPtr( %VAL( Q ), mxGetPr( PLHS( IP ) ),
     $                             N*N )
            IP = IP + 1
         END IF
C
         IF ( NLHS.GE.IP+1 .AND. ICMPU.GT.0 ) THEN
            PLHS( IP ) = mxCreateDoubleMatrix( M, M, 0 )
            CALL mxCopyReal8ToPtr( %VAL( U1 ),
     $                             mxGetPr( PLHS( IP ) ), M*M )
            IP = IP + 1
            PLHS( IP ) = mxCreateDoubleMatrix( M, M, 0 )
            CALL mxCopyReal8ToPtr( %VAL( U2 ),
     $                             mxGetPr( PLHS( IP ) ), M*M )
         END IF
      END IF
C
C Deallocate variable dimension arrays.
C !Fortran 90/95
C
      DEALLOCATE( IWORK )
      CALL mxFree( ALPHAI )
      CALL mxFree( ALPHAR )
      CALL mxFree( B      )
      CALL mxFree( BETA   )
      CALL mxFree( DWORK  )
      CALL mxFree( FG     )
      CALL mxFree( Q      )
      CALL mxFree( U1     )
      CALL mxFree( U2     )
      CALL mxFree( Z      )
C
C Error and warning handling.
C
      IF ( INFO.NE.0 ) THEN
         WRITE( TEXT, '('' INFO = '',I4,'' ON EXIT FROM DGHFST'')'
     $        ) INFO
         CALL mexErrMsgTxt( TEXT )
      END IF
C
      RETURN
C *** Last line of skewHamil2feig ***
      END
 
