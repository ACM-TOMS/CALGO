#include "fintrf.h"
C SKEWHAMILEIG.F - Gateway function for computing the eigenvalues of a
C                  skew-Hamiltonian/Hamiltonian pencil, using the routine
C                  DGHUTR.
C
C Matlab call:
C   [ALPHAR,ALPHAI,BETA(,Ao,Do,Bo,Fo,C1o,Vo,C2o(,Q1,Q2))] =
C                      skewHamileig(A,DE,C,VW(,job,compq1(,Q),compq2))
C
C   [ALPHAR,ALPHAI,BETA]                 = skewHamileig(A,DE,C,VW)
C   [ALPHAR,ALPHAI,BETA,Ao,Do,Bo,Fo,C1o,Vo,C2o] ...
C                                        = skewHamileig(A,DE,C,VW,1)
C
C Purpose:
C   To compute the eigenvalues of a real n-by-n skew-Hamiltonian/
C   Hamiltonian pencil aS - bH with n = 2m,
C
C         (  A  D  )         (  C  V  )
C     S = (        ) and H = (        ).                             (1)
C         (  E  A' )         (  W -C' )
C
C   Optionally, if job = 1, decompositions of S and H will be
C   computed via orthogonal transformations Q1 and Q2 as follows:
C
C                     (  Ao  Do  )
C     Q1' S J Q1 J' = (          ),
C                     (   0  Ao' )
C
C                     (  Bo  Fo  )
C     J' Q2' J S Q2 = (          ) =: T,                             (2)
C                     (   0  Bo' )
C
C                (  C1o  Vo   )            (  0  I  )
C     Q1' H Q2 = (            ), where J = (        ),
C                (  0    C2o' )            ( -I  0  )
C
C   and Ao, Bo, C1o are upper triangular, C2o is upper quasi-triangular
C   and Do and Fo are skew-symmetric. The notation M' denotes the
C   transpose of the matrix M.
C   Optionally, if compq1 = 1, the orthogonal transformation matrix
C   Q1 will be computed.
C   Optionally, if compq2 = 1, the orthogonal transformation matrix
C   Q2 will be computed.
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
C   C      - the m-by-m matrix C.
C   VW     - an  m-by-(m+1) matrix containing the triangles of the
C            symmetric matrices V and W, as follows:
C            the leading m-by-m lower triangular part contains the lower
C            triangle of the matrix W, and the m-by-m upper triangular
C            part of the submatrix in the columns 2 to m+1 contains the
C            upper triangle of the matrix V of H in (1).
C            So, if i >= j, then W(i,j) = W(j,i) is stored in VW(i,j)
C            and V(j,i) = V(i,j) is stored in VW(j,i+1).
C            VW is an empty matrix if m = 0.
C   job    - (optional) scalar indicating the computation to be
C            performed, as follows:
C            = 0 :  compute the eigenvalues only (default);
C            = 1 :  compute the eigenvalues and the matrices of the
C                   transformed pencil in (2).
C   compq1 - (optional) scalar indicating whether the orthogonal
C            transformation matrix Q1 is returned or accumulated
C            into an orthogonal matrix, or if Q1 is not required,
C            as follows:
C            = 0 :  Q1 is not required (default);
C            = 1 :  on entry, Q1 need not be set;
C                   on exit, Q1 contains the orthogonal matrix Q1;
C            = 2 :  the orthogonal transformations are accumulated
C                   into Q1;
C                   on input, Q1 must contain an orthogonal matrix Q;
C                   on exit, Q1 contains Q*Q1.
C   Q      - if compq1 = 2, the n-by-n orthogonal matrix Q.
C   compq2 - (optional) scalar indicating whether the orthogonal
C            transformation matrix Q2 is returned or accumulated
C            into an orthogonal matrix, or if Q2 is not required,
C            as follows:
C            = 0 :  Q2 is not required (default);
C            = 1 :  on exit, Q2 contains the orthogonal matrix Q2;
C            = 2 :  on exit, Q2 contains the matrix product J*Q*J'*Q2,
C                   where Q2 is the product of the orthogonal
C                   transformations that are applied to the pencil
C                   aS - bH to reduce S and H to the forms in (2),
C                   for compq2 = 1.
C                   If nonzero, compq2 should coincide with compq1.
C
C Output parameters:
C   ALPHAR,- the m-vectors of real parts and imaginary parts,
C   ALPHAI   respectively, of each scalar alpha defining an eigenvalue
C            of the pencil aS - bH.
C            If ALPHAI(j) is zero, then the j-th eigenvalue is real.
C   BETA     the m-vector of the scalars beta that define the eigenvalues
C            of the pencil aS - bH.
C            Together, the quantities alpha = (ALPHAR(j),ALPHAI(j)) and
C            beta = BETA(j) represent the j-th eigenvalue of the pencil
C            aS - bH, in the form lambda = alpha/beta. Since lambda may
C            overflow, the ratios should not, in general, be computed.
C            Due to the skew-Hamiltonian/Hamiltonian structure of the
C            pencil, for every eigenvalue lambda, -lambda is also an
C            eigenvalue, and thus it has only to be saved once in
C            ALPHAR, ALPHAI and BETA.
C            Specifically, only eigenvalues with imaginary parts greater
C            than or equal to zero are stored; their conjugate 
C            eigenvalues are not stored. If imaginary parts are zero
C            (i.e., for real eigenvalues), only positive eigenvalues
C            are stored.
C   Ao     - if job = 1, the computed m-by-m submatrix Ao in (2).
C   Do     - if job = 1, the computed m-by-m strictly upper triangular
C            part of the skew-symmetric matrix Do in (2).
C   Bo     - if job = 1, the computed m-by-m submatrix Bo in (2).
C   Fo     - if job = 1, the computed m-by-m strictly upper triangular
C            part of the skew-symmetric matrix Fo in (2).
C   C1o    - if job = 1, the computed m-by-m submatrix C1o in (2).
C   Vo     - if job = 1, the computed m-by-m matrix Vo in (2).
C   C2o    - if job = 1, the computed m-by-m submatrix C2o in (2).
C   Q1     - if compq1 > 0, an n-by-n matrix containing the computed
C            orthogonal matrix Q1 (if compq1 = 1) or Q*Q1
C            (if compq1 = 2).
C   Q2     - if compq2 > 0, an n-by-n matrix containing the computed
C            orthogonal matrix Q2 (if compq2 = compq1 = 1) or J*Q*J'*Q2
C            (if compq2 = compq1 = 2).
C            If job = 0, Ao, Do, Bo, Fo, C1o, Vo, and C2o contain the
C            corresponding matrices just before the application of the
C            periodic QZ algorithm.
C
C Contributors:
C   V. Sima, Research Institute for Informatics, Bucharest, Oct. 2009.
C
C Revisions:
C   V. Sima, Oct. 2010, Nov. 2010, Dec. 2010, July 2011, Sep. 2012,
C   Jan. 2013, Feb. 2013, July 2014.
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
C     .. Scalar parameters used by the subroutines ..
      CHARACTER         COMPQ1, COMPQ2, JOB
      INTEGER           INFO, LDA, LDB, LDC1, LDC2, LDDE, LDF, LDQ1,
     $                  LDQ2, LDVW, LDWORK, LIWORK, N 
C
C     .. Allocatable arrays ..
C     !Fortran 90/95 (Fixed dimensions should be used with Fortran 77.)
      INTEGER,          ALLOCATABLE :: IWORK( : )
      mwPointer         A, ALPHAI, ALPHAR, B, BETA, C1, C2, DE, DW,
     $                  DWORK, F, Q1, Q2, VW
C
C     .. Local variables and constant dimension arrays ..
      CHARACTER*120     TEXT
      INTEGER           I, ICMPQ1, ICMPQ2, IDW, IJOB, IP, LDQ, M, M1, N1
      DOUBLE PRECISION  TEMP
C
C     .. External Subroutines ..
      EXTERNAL          DLACPY, DLACPV, DGHUTR
C
C     ..Intrinsic Functions..
      INTRINSIC         MAX, MIN, MOD
C
C     Check for proper number of arguments.
C 
      IF( NRHS.LT.4 ) THEN
         CALL mexErrMsgTxt
     $        ( 'SKEWHAMILEIG requires at least 4 input arguments.' )
      ELSE IF ( NLHS.LT.3 ) THEN
         CALL mexErrMsgTxt
     $        ( 'SKEWHAMILEIG requires at least 3 output arguments.' )
      END IF
C
C   A(mxm), DE(mx(m+1)), C1(mxm), VW(mx(m+1))(, job, compq1(, Q(nxn)),
C   compq2).
C
      M  = mxGetM( PRHS(1) )
      N1 = mxGetN( PRHS(1) )
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
         CALL mexErrMsgTxt( 'C must be a real matrix' )
      END IF
      IF ( M1.NE.M ) THEN
         CALL mexErrMsgTxt( 'C must have the same number of rows as A' )
      END IF
      IF ( N1.NE.M ) THEN
         CALL mexErrMsgTxt( 'C must be a square matrix' )
      END IF
C
      M1 = mxGetM( PRHS(4) )
      N1 = mxGetN( PRHS(4) )
      IF ( mxIsNumeric( PRHS(4) ).EQ.0 .OR.
     $     mxIsComplex( PRHS(4) ).EQ.1 ) THEN
         CALL mexErrMsgTxt( 'VW must be a real matrix' )
      END IF
      IF ( M1.NE.M ) THEN
         CALL mexErrMsgTxt( 'VW must have the same number of rows as A'
     $                    )
      END IF
      IF ( N1.NE.M+1 ) THEN
         CALL mexErrMsgTxt( 'VW must have one more columns than rows' )
      END IF
C
      IP = 5
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
            CALL mexErrMsgTxt( 'COMPQ1 must be a scalar' )
         END IF
         IF ( mxIsNumeric( PRHS(IP) ).EQ.0 .OR.
     $        mxIsComplex( PRHS(IP) ).EQ.1 ) THEN
            CALL mexErrMsgTxt( 'COMPQ1 must be an integer scalar' )
         END IF
         CALL mxCopyPtrToReal8( mxGetPr( PRHS(IP) ), TEMP, 1 )
         ICMPQ1 = TEMP
         IF ( ICMPQ1.LT.0 .OR. ICMPQ1.GT.2 ) THEN
            CALL mexErrMsgTxt
     $           ( 'COMPQ1 has 0, 1 or 2 the only admissible values' )
         END IF
         IP = IP + 1
C
         IF ( NRHS.GE.IP .AND. ICMPQ1.EQ.2 ) THEN
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
            IP = IP + 1
         END IF
C
      ELSE
         ICMPQ1 = 0
      END IF
C
      IF ( ICMPQ1.EQ.0 ) THEN
         COMPQ1 = 'N'
      ELSE IF ( ICMPQ1.EQ.1 ) THEN
         COMPQ1 = 'I'
      ELSE
         COMPQ1 = 'U'
      END IF
C
      IF ( NRHS.GE.IP ) THEN
         IF ( mxGetM( PRHS(IP) ).NE.1 .OR.
     $        mxGetN( PRHS(IP) ).NE.1 ) THEN
            CALL mexErrMsgTxt( 'COMPQ2 must be a scalar' )
         END IF
         IF ( mxIsNumeric( PRHS(IP) ).EQ.0 .OR.
     $        mxIsComplex( PRHS(IP) ).EQ.1 ) THEN
            CALL mexErrMsgTxt( 'COMPQ2 must be an integer scalar' )
         END IF
         CALL mxCopyPtrToReal8( mxGetPr( PRHS(IP) ), TEMP, 1 )
         ICMPQ2 = TEMP
         IF ( ICMPQ2.LT.0 .OR. ICMPQ2.GT.2 ) THEN
            CALL mexErrMsgTxt
     $           ( 'COMPQ2 has 0, 1 or 2 the only admissible values' )
         END IF
         IF ( ICMPQ2.NE.0 .AND. ICMPQ2.NE.ICMPQ1 ) THEN
            CALL mexErrMsgTxt( 'COMPQ2 must coincide to COMPQ1' )
         END IF
C
      ELSE
         ICMPQ2 = 0
      END IF
C
      IF ( ICMPQ2.EQ.0 ) THEN
         COMPQ2 = 'N'
      ELSE IF ( ICMPQ2.EQ.1 ) THEN
         COMPQ2 = 'I'
      ELSE
         COMPQ2 = 'U'
      END IF
C
C Determine dimensions of the arrays and the workspace.
C
      LDA  = MAX( 1, M )
      LDDE = LDA
      LDC1 = LDA
      LDVW = LDA
      LDB  = LDA
      LDF  = LDA
      LDC2 = LDA
      LDQ  = MAX( 1, N )
      IF ( ICMPQ1.EQ.0 ) THEN
         LDQ1 = 1
         M1   = 0
      ELSE
         LDQ1 = LDQ
         M1   = N
      END IF
      IF ( ICMPQ2.EQ.0 ) THEN
         LDQ2 = 1
         N1   = 0
      ELSE
         LDQ2 = LDQ
         N1   = N
      END IF
      LIWORK = N + 12
      IF( MOD( M, 2 ).EQ.0 ) THEN
         I = MAX( 4*N, 32 ) + 4
      ELSE
         I = MAX( 4*N, 36 )
      END IF
      IF( IJOB.EQ.1 .OR. ICMPQ1.GT.0 .OR. ICMPQ2.GT.0 ) THEN
         LDWORK = 2*N*N + I
      ELSE
         LDWORK = N*N + I
      END IF
C
C Allocate variable dimension local arrays.
C !Fortran 90/95
C
      ALLOCATE( IWORK( LIWORK ) )
      A      = mxCalloc(  LDA*M,     8 )
      ALPHAI = mxCalloc(      M,     8 )
      ALPHAR = mxCalloc(      M,     8 )
      B      = mxCalloc(  LDB*M,     8 )
      BETA   = mxCalloc(      M,     8 )
      C1     = mxCalloc( LDC1*M,     8 )
      C2     = mxCalloc( LDC2*M,     8 )
      DE     = mxCalloc( LDDE*(M+1), 8 )
      DW     = mxCalloc(    M*M,     8 )
      DWORK  = mxCalloc( LDWORK,     8 )
      F      = mxCalloc(  LDF*M,     8 )
      Q1     = mxCalloc( LDQ1*M1,    8 )
      Q2     = mxCalloc( LDQ2*N1,    8 )
      VW     = mxCalloc( LDVW*(M+1), 8 )
C
C Copy inputs from MATLAB workspace to locally allocated arrays.
C
      CALL mxCopyPtrToReal8( mxGetPr( PRHS( 1 ) ), %VAL( A ),  M*M )
      CALL mxCopyPtrToReal8( mxGetPr( PRHS( 2 ) ), %VAL( DE ),
     $                       M*( M + 1 ) )
      CALL mxCopyPtrToReal8( mxGetPr( PRHS( 3 ) ), %VAL( C1 ), M*M )
      CALL mxCopyPtrToReal8( mxGetPr( PRHS( 4 ) ), %VAL( VW ),
     $                       M*( M + 1 ) )
      IF ( NRHS.GE.7 .AND. ICMPQ1.EQ.2 )
     $   CALL mxCopyPtrToReal8( mxGetPr( PRHS( 7 ) ), %VAL( Q1 ),
     $                          N*N )
C
C Do the actual computations.
C
      CALL DGHUTR( JOB, COMPQ1, COMPQ2, N, %VAL( A ), LDA, %VAL( DE ),
     $             LDDE, %VAL( C1 ), LDC1, %VAL( VW ), LDVW, %VAL( Q1 ),
     $             LDQ1, %VAL( Q2 ), LDQ2, %VAL( B ), LDB, %VAL( F ),
     $             LDF, %VAL( C2 ), LDC2, %VAL( ALPHAR ),
     $             %VAL( ALPHAI ), %VAL( BETA ), IWORK, LIWORK,
     $             %VAL( DWORK ), LDWORK, INFO )
C
C Copy output to MATLAB workspace.
C
      IF ( INFO.EQ.0 .OR. INFO.EQ.3 ) THEN
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
            IDW = MIN( M*( M+1 ), M+1 )
            CALL DLACPV( M*M, IDW, 1, %VAL( DE ), %VAL( DW ) )
            PLHS( IP ) = mxCreateDoubleMatrix( M, M, 0 )
            CALL mxCopyReal8ToPtr( %VAL( DW ),
     $                             mxGetPr( PLHS( IP ) ), M*M )
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
            PLHS( IP ) = mxCreateDoubleMatrix( M, M, 0 )
            CALL mxCopyReal8ToPtr( %VAL( F ), mxGetPr( PLHS( IP ) ),
     $                             M*M )
            IP = IP + 1
         END IF
C
         IF ( NLHS.GE.IP ) THEN
            PLHS( IP ) = mxCreateDoubleMatrix( M, M, 0 )
            CALL mxCopyReal8ToPtr( %VAL( C1 ), mxGetPr( PLHS( IP ) ),
     $                             M*M )
            IP = IP + 1
         END IF
C
         IF ( NLHS.GE.IP ) THEN
            CALL DLACPV( M*M, IDW, 1, %VAL( VW ), %VAL( DW ) )
            PLHS( IP ) = mxCreateDoubleMatrix( M, M, 0 )
            CALL mxCopyReal8ToPtr( %VAL( DW ),
     $                             mxGetPr( PLHS( IP ) ), M*M )
            IP = IP + 1
         END IF
C
         IF ( NLHS.GE.IP ) THEN
            PLHS( IP ) = mxCreateDoubleMatrix( M, M, 0 )
            CALL mxCopyReal8ToPtr( %VAL( C2 ), mxGetPr( PLHS( IP ) ),
     $                             M*M )
            IP = IP + 1
         END IF
C
         IF ( NLHS.GE.IP .AND. ICMPQ1.GE.1 ) THEN
            PLHS( IP ) = mxCreateDoubleMatrix( N, N, 0 )
            CALL mxCopyReal8ToPtr( %VAL( Q1 ), mxGetPr( PLHS( IP ) ), 
     $                             N*N )
            IP = IP + 1
         END IF
C
         IF ( NLHS.GE.IP .AND. ICMPQ2.GE.1 ) THEN
            PLHS( IP ) = mxCreateDoubleMatrix( N, N, 0 )
            CALL mxCopyReal8ToPtr( %VAL( Q2 ), mxGetPr( PLHS( IP ) ),
     $                             N*N )
         END IF
      END IF
C
C Deallocate variable dimension arrays.
C !Fortran 90/95
C
      DEALLOCATE(  IWORK  )
      CALL mxFree( A      )
      CALL mxFree( ALPHAI )
      CALL mxFree( ALPHAR )
      CALL mxFree( B      )
      CALL mxFree( BETA   )
      CALL mxFree( C1     )
      CALL mxFree( C2     )
      CALL mxFree( DE     )
      CALL mxFree( DW     )
      CALL mxFree( DWORK  )
      CALL mxFree( F      )
      CALL mxFree( Q1     )
      CALL mxFree( Q2     )
      CALL mxFree( VW     )
C
C Error and warning handling.
C
      IF ( INFO.NE.0 ) THEN
         WRITE( TEXT, '('' INFO = '',I4,'' ON EXIT FROM DGHUTR'')'
     $        ) INFO
         IF ( INFO.EQ.3 ) THEN
            CALL mexPrintf( TEXT )
         ELSE
            CALL mexErrMsgTxt( TEXT )
         END IF
      END IF
C
      RETURN
C *** Last line of skewHamileig ***
      END
