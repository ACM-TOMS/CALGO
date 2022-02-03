#include "fintrf.h"
C SKEWHAMIL2EIG_NB.F - Gateway function for computing the eigenvalues 
C                      of a real skew-Hamiltonian/skew-Hamiltonian
C                      pencil using the routine DGHUSP ("block" 
C                      routine).
C                      The block size, NB, can be specified.
C
C Matlab call:
C   [ALPHAR,ALPHAI,BETA(,Ao,Do,Bo,Fo(,Q))] =
C                      skewHamil2eig_nb(A,DE,B,FG(,nb,job,compq(,Q)))
C
C   [ALPHAR,ALPHAI,BETA]             = skewHamil2eig_nb(A,DE,B,FG)
C   [ALPHAR,ALPHAI,BETA,Ao,Do,Bo,Fo] = skewHamil2eig_nb(A,DE,B,FG,nb,1)
C
C Purpose:
C   To compute the eigenvalues of a real n-by-n skew-Hamiltonian/
C   skew-Hamiltonian pencil aS - bT, with n = 2m,
C
C         (  A  D  )         (  B  F  )
C     S = (        ) and T = (        ).                             (1)
C         (  E  A' )         (  G  B' )
C
C   Optionally, if job = 1, decompositions of S and T will be
C   computed via an orthogonal transformation Q as follows:
C
C                   (  Ao  Do  )
C     J Q' J' S Q = (          ),
C                   (   0  Ao' )
C
C                   (  Bo  Fo  )            (  0  I  )
C     J Q' J' T Q = (          ), where J = (        ),              (2)
C                   (   0  Bo' )            ( -I  0  )
C
C   and Ao is upper triangular, Bo is upper quasi-triangular, and Do and
C   Fo are skew-symmetric. The notation M' denotes the transpose of the
C   matrix M.
C   Optionally, if compq = 1, the orthogonal transformation matrix Q
C   will be computed.
C
C Input parameters:
C   A      - the m-by-m matrix A, with m = n/2.
C   DE     - an  m-by-(m+1) matrix containing the strict triangles of
C            the skew-symmetric matrices D and E, as follows:
C            the leading m-by-m strictly lower triangular part contains
C            the strictly lower triangle of the matrix E, and the
C            m-by-m strictly upper triangular part of the submatrix in
C            the columns 2 to m+1 contains the strictly upper triangle
C            of the matrix D of S in (1).
C            So, if i > j, then E(i,j) = -E(j,i) is stored in DE(i,j)
C            and D(j,i) = -D(i,j) is stored in DE(j,i+1).
C            The entries on the diagonal and the first superdiagonal of
C            DE need not be set, but are assumed to be zero.
C            DE is an empty matrix if m = 0.
C   B      - the m-by-m matrix B.
C   FG     - an  m-by-(m+1) matrix containing the triangles of the
C            skew-symmetric matrices F and G, as follows:
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
C   nb     - (optional) scalar indicating the block size. An internally
C            computed block size will be used for a non-positive value.
C            Default:  nb = 32.
C   job    - (optional) scalar indicating the computation to be
C            performed, as follows:
C            = 0 :  compute the eigenvalues only (default);
C            = 1 :  compute the eigenvalues and the matrices of the
C                   transformed pencil in (2).
C   compq  - (optional) scalar indicating whether the orthogonal
C            transformation matrix Q is returned or accumulated
C            into an orthogonal matrix, or if Q is not required,
C            as follows:
C            = 0 :  Q is not required (default);
C            = 1 :  on entry, Q need not be set;
C                   on exit, Q contains the orthogonal matrix Q;
C            = 2 :  the orthogonal transformations are accumulated
C                   into Q;
C                   on input, Q must contain an orthogonal matrix Q0;
C                   on exit, Q contains Q0*Q.
C   Q      - if compq = 2, the n-by-n orthogonal matrix Q0.
C
C Output parameters:
C   ALPHAR,- the m-vectors of real parts and imaginary parts,
C   ALPHAI   respectively, of each scalar alpha defining an eigenvalue
C            of the pencil aS - bT.
C            If ALPHAI(j) is zero, then the j-th eigenvalue is real; if
C            positive, then the j-th and (j+1)-st eigenvalues are a
C            complex conjugate pair.
C   BETA     the m-vector of the scalars beta that define the
C            eigenvalues of the pencil aS - bT.
C            Together, the quantities alpha = (ALPHAR(j),ALPHAI(j)) and
C            beta = BETA(j) represent the j-th eigenvalue of the pencil
C            aS - bT, in the form lambda = alpha/beta. Since lambda may
C            overflow, the ratios should not, in general, be computed.
C            Due to the skew-Hamiltonian/skew-Hamiltonian structure of
C            the pencil, every eigenvalue occurs twice and thus it has
C            only to be saved once in ALPHAR, ALPHAI and BETA.
C   Ao     - if job = 1, the computed m-by-m submatrix Ao in (2).
C   Do     - if job = 1, the computed m-by-m strictly upper triangular
C            part of the skew-symmetric matrix Do in (2). The lower
C            triangle is not set.
C   Bo     - if job = 1, the computed m-by-m submatrix Bo in (2).
C   Fo     - if job = 1, the computed m-by-m strictly upper triangular
C            part of the skew-symmetric matrix Fo in (2). The lower
C            triangle is not set.
C   Q      - if compq > 0, an n-by-n matrix containing the computed
C            orthogonal matrix Q (if compq = 1) or Q0*Q (if compq = 2).
C            If job = 0, Do and Fo contain the corresponding matrices
C            just before the application of the QZ algorithm, but Ao and
C            Bo contain meaningless elements.
C
C Contributors:
C   V. Sima, Research Institute for Informatics, Bucharest, Nov. 2009.
C
C Revisions:
C   V. Sima, Dec. 2010, July 2013, July 2014.
C   M. Voigt, July 2013.
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
      mwPointer         mxCalloc, mxCreateDoubleMatrix,
     $                  mxGetPr
      INTEGER*4         mxGetM, mxGetN, mxIsNumeric, mxIsComplex
C
C     .. Scalar parameters used by subroutines ..
      CHARACTER         COMPQ, JOB
      INTEGER           INFO, LDA, LDB, LDDE, LDFG, LDQ, LDWORK, N
C
C     .. Allocatable arrays ..
C     !Fortran 90/95 (Fixed dimensions should be used with Fortran 77.)
      DOUBLE PRECISION, ALLOCATABLE :: DW( : )
      mwPointer         A, ALPHAI, ALPHAR, B, BETA, DE, DWORK, FG, Q
C
C     .. Local variables and constant dimension arrays ..
      CHARACTER*120     TEXT
      INTEGER           ICMPQ, IDW, IJOB, IP, M, M1, N1, NB
      DOUBLE PRECISION  DUM( 1 ), TEMP
C
C     .. External Subroutines ..
      EXTERNAL          DGHUSP
C
C     ..Intrinsic Functions..
      INTRINSIC         INT, MAX, MIN
C
C     Check for proper number of arguments.
C
      IF( NRHS.LT.4 ) THEN
         CALL mexErrMsgTxt
     $        ( 'SKEWHAMIL2EIG requires at least 4 input arguments.' )
      ELSE IF ( NLHS.LT.3 ) THEN
         CALL mexErrMsgTxt
     $        ( 'SKEWHAMIL2EIG requires at least 3 output arguments.' )
      END IF
C
C   A(mxm), DE(mx(m+1)), B(mxm), FG(mx(m+1))(, nb, job, compq(, Q(nxn)) ).
C
      M  = mxGetM( PRHS( 1 ) )
      N1 = mxGetN( PRHS( 1 ) )
      IF ( mxIsNumeric( PRHS(1) ).EQ.0 .OR.
     $     mxIsComplex( PRHS(1) ).EQ.1 ) THEN
         CALL mexErrMsgTxt( 'A must be a real matrix' )
      END IF
      IF ( M.NE.N1 ) THEN
         CALL mexErrMsgTxt( 'A must be a square matrix' )
      END IF
      N  = 2*M
C
      M1 = mxGetM( PRHS(2) )
      N1 = mxGetN( PRHS(2) )
      IF ( mxIsNumeric( PRHS(2) ).EQ.0 .OR.
     $     mxIsComplex( PRHS(2) ).EQ.1 ) THEN
         CALL mexErrMsgTxt( 'DE must be a real matrix' )
      END IF
      IF ( M1.NE.M ) THEN
         CALL mexErrMsgTxt( 'DE must have the same number of rows as A'
     $                    )
      END IF
      IF ( N1.NE.M+1 ) THEN
         CALL mexErrMsgTxt( 'DE must have one more columns than rows' )
      END IF
C
      M1 = mxGetM( PRHS(3) )
      N1 = mxGetN( PRHS(3) )
      IF ( mxIsNumeric( PRHS(3) ).EQ.0 .OR.
     $     mxIsComplex( PRHS(3) ).EQ.1 ) THEN
         CALL mexErrMsgTxt( 'B must be a real matrix' )
      END IF
      IF ( M1.NE.M ) THEN
         CALL mexErrMsgTxt( 'B must have the same number of rows as A' )
      END IF
      IF ( N1.NE.M ) THEN
         CALL mexErrMsgTxt( 'B must be a square matrix' )
      END IF
C
      M1 = mxGetM( PRHS(4) )
      N1 = mxGetN( PRHS(4) )
      IF ( mxIsNumeric( PRHS(4) ).EQ.0 .OR.
     $     mxIsComplex( PRHS(4) ).EQ.1 ) THEN
         CALL mexErrMsgTxt( 'FG must be a real matrix' )
      END IF
      IF ( M1.NE.M ) THEN
         CALL mexErrMsgTxt( 'FG must have the same number of rows as A'
     $                    )
      END IF
      IF ( N1.NE.M+1 ) THEN
         CALL mexErrMsgTxt( 'FG must have one more columns than rows' )
      END IF
C
      IP = 5
      IF ( NRHS.GE.IP ) THEN
         IF ( mxGetM( PRHS(IP) ).NE.1 .OR.
     $        mxGetN( PRHS(IP) ).NE.1 ) THEN
            CALL mexErrMsgTxt( 'NB must be a scalar' )
         END IF
         IF ( mxIsNumeric( PRHS(IP) ).EQ.0 .OR.
     $        mxIsComplex( PRHS(IP) ).EQ.1 ) THEN
            CALL mexErrMsgTxt( 'NB must be an integer scalar' )
         END IF
         CALL mxCopyPtrToReal8( mxGetPr( PRHS(IP) ), TEMP, 1 )
         NB = TEMP
         IP = IP + 1
C
      ELSE
         NB = 32
      END IF
C
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
         IF ( ICMPQ.LT.0 .OR. ICMPQ.GT.2 ) THEN
            CALL mexErrMsgTxt
     $           ( 'COMPQ has 0, 1 or 2 the only admissible values' )
         END IF
         IP = IP + 1
C
         IF ( NRHS.GE.IP .AND. ICMPQ.EQ.2 ) THEN
            M1 = mxGetM( PRHS(IP) )
            N1 = mxGetN( PRHS(IP) )
            IF ( mxIsNumeric( PRHS(IP) ).EQ.0 .OR.
     $           mxIsComplex( PRHS(IP) ).EQ.1 ) THEN
               CALL mexErrMsgTxt( 'Q must be a real matrix' )
            END IF
            IF ( M1.NE.N ) THEN
               CALL mexErrMsgTxt( 'Q must have twice as many rows as A'
     $                          )
            END IF
            IF ( N1.NE.N ) THEN
               CALL mexErrMsgTxt( 'Q must be a square matrix' )
            END IF
         END IF
C
      ELSE
         ICMPQ = 0
      END IF
C
      IF ( ICMPQ.EQ.0 ) THEN
         COMPQ = 'N'
      ELSE IF ( ICMPQ.EQ.1 ) THEN
         COMPQ = 'I'
      ELSE
         COMPQ = 'U'
      END IF
C
C Determine dimensions of the arrays and the workspace.
C
      LDA  = MAX( 1, M )
      LDDE = LDA
      LDB  = LDA
      LDFG = LDA
      IF ( ICMPQ.EQ.0 ) THEN
         LDQ = 1
         M1  = 0
      ELSE
         LDQ = MAX( 1, N )
         M1  = N
      END IF
C
C Allocate variable dimension local arrays.
C !Fortran 90/95
C
      ALLOCATE( DW( LDDE*(M+1) ) )
      A      = mxCalloc(  LDA*M,     8 )
      ALPHAI = mxCalloc(      M,     8 )
      ALPHAR = mxCalloc(      M,     8 )
      B      = mxCalloc(  LDB*M,     8 )
      BETA   = mxCalloc(      M,     8 )
      DE     = mxCalloc( LDDE*(M+1), 8 )
      FG     = mxCalloc( LDFG*(M+1), 8 )
      Q      = mxCalloc(  LDQ*M1,    8 )
C
      CALL DGHUSP( JOB, COMPQ, N, %VAL( A ), LDA, %VAL( DE ), LDDE,
     $             %VAL( B ), LDB, %VAL( FG ), LDFG, %VAL( Q ), LDQ,
     $             %VAL( ALPHAR ), %VAL( ALPHAI ), %VAL( BETA ), DUM,
     $             -1, INFO )
      LDWORK = INT( DUM(1) )
      DWORK  = mxCalloc( LDWORK, 8 )
C
C Copy inputs from MATLAB workspace to locally allocated arrays.
C
      CALL mxCopyPtrToReal8( mxGetPr( PRHS( 1 ) ), %VAL( A ),  M*M )
      CALL mxCopyPtrToReal8( mxGetPr( PRHS( 2 ) ), %VAL( DE ),
     $                       M*( M + 1 ) )
      CALL mxCopyPtrToReal8( mxGetPr( PRHS( 3 ) ), %VAL( B ),  M*M )
      CALL mxCopyPtrToReal8( mxGetPr( PRHS( 4 ) ), %VAL( FG ),
     $                       M*( M + 1 ) )
      IF ( NRHS.GE.7 .AND. ICMPQ.EQ.2 )
     $   CALL mxCopyPtrToReal8( mxGetPr( PRHS( 7 ) ), %VAL( Q ),
     $                          N*N )
C
C Do the actual computations.
C
      INFO = NB
      CALL DGHUSP( JOB, COMPQ, N, %VAL( A ), LDA, %VAL( DE ), LDDE,
     $             %VAL( B ), LDB, %VAL( FG ), LDFG, %VAL( Q ), LDQ,
     $             %VAL( ALPHAR ), %VAL( ALPHAI ), %VAL( BETA ),
     $             %VAL( DWORK ), LDWORK, INFO )
C
C Copy output to MATLAB workspace.
C
      IF ( INFO.EQ.0 ) THEN
         PLHS( 1 ) = mxCreateDoubleMatrix( M, 1, 0 )
         CALL mxCopyReal8ToPtr( %VAL( ALPHAR ), mxGetPr( PLHS( 1 ) ),
     $                          M )
C
         PLHS( 2 ) = mxCreateDoubleMatrix( M, 1, 0 )
         CALL mxCopyReal8ToPtr( %VAL( ALPHAI ), mxGetPr( PLHS( 2 ) ),
     $                          M )
C
         PLHS( 3 ) = mxCreateDoubleMatrix( M, 1, 0 )
         CALL mxCopyReal8ToPtr( %VAL( BETA ), mxGetPr( PLHS( 3 ) ),
     $                          M )
C
         IP = 4
         IF ( NLHS.GE.IP ) THEN
            PLHS( IP ) = mxCreateDoubleMatrix( M, M, 0 )
            CALL mxCopyReal8ToPtr( %VAL( A ), mxGetPr( PLHS( IP ) ),
     $                             M*M )
            IP = IP + 1
         END IF
C
         IF ( NLHS.GE.IP ) THEN
            CALL DLACPY( 'Full', M, M+1, %VAL( DE ), LDDE, DW, LDA )
            PLHS( IP ) = mxCreateDoubleMatrix( M, M, 0 )
            IDW = MIN( M*( M+1 ), M+1 )
            CALL mxCopyReal8ToPtr( DW( IDW ), mxGetPr( PLHS( IP ) ),
     $                             M*M )
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
            CALL DLACPY( 'Full', M, M+1, %VAL( FG ), LDFG, DW, LDA )
            PLHS( IP ) = mxCreateDoubleMatrix( M, M, 0 )
            CALL mxCopyReal8ToPtr( DW( IDW ), mxGetPr( PLHS( IP ) ),
     $                             M*M )
            IP = IP + 1
         END IF
C
         IF ( NLHS.GE.IP .AND. ICMPQ.GE.1 ) THEN
            PLHS( IP ) = mxCreateDoubleMatrix( N, N, 0 )
            CALL mxCopyReal8ToPtr( %VAL( Q ), mxGetPr( PLHS( IP ) ),
     $                             N*N )
         END IF
      END IF
C
C Deallocate variable dimension arrays.
C !Fortran 90/95
C
      DEALLOCATE( DW )
      CALL mxFree( A      )
      CALL mxFree( ALPHAI )
      CALL mxFree( ALPHAR )
      CALL mxFree( B      )
      CALL mxFree( BETA   )
      CALL mxFree( DE     )
      CALL mxFree( DWORK  )
      CALL mxFree( FG     )
      CALL mxFree( Q      )
C
C Error and warning handling.
C
      IF ( INFO.NE.0 ) THEN
         WRITE( TEXT, '('' INFO = '',I4,'' ON EXIT FROM DGHUSP'')'
     $        ) INFO
         CALL mexErrMsgTxt( TEXT )
      END IF
C
      RETURN
C *** Last line of skewHamil2eig_nb ***
      END
