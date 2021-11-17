#include "fintrf.h"
C SKEWHAMILDEFLZ.F - Gateway function for computing the eigenvalues of a
C                    complex skew-Hamiltonian/Hamiltonian pencil and the
C                    right deflating subspace corresponding to the
C                    eigenvalues with strictly negative real part, using
C                    the routine ZGHUDF.
C                    The gateway accepts real or complex input matrices.
C
C Matlab call:
C   [ALPHA,BETA(,Q,neig)] = skewHamildeflZ(A,DE,B,FG(,compq(,orthm)))
C
C Purpose:
C   To compute the eigenvalues of a complex n-by-n skew-Hamiltonian/
C   Hamiltonian pencil aS - bH, with n = 2m,
C
C         (  A  D  )         (  B  F  )
C     S = (        ),    H = (        ),                             (1)
C         (  E  A' )         (  G -B' )
C
C   where M' denotes the conjugate transpose of the matrix M.
C   The structured Schur form of the embedded real skew-Hamiltonian/
C   skew-Hamiltonian pencil aB_S - bB_T, defined as
C
C
C           (  Re(A)  -Im(A)  |  Re(D)  -Im(D)  )
C           (                 |                 )
C           (  Im(A)   Re(A)  |  Im(D)   Re(D)  )
C           (                 |                 )
C     B_S = (-----------------+-----------------) , and
C           (                 |                 )
C           (  Re(E)  -Im(E)  |  Re(A')  Im(A') )
C           (                 |                 )
C           (  Im(E)   Re(E)  | -Im(A')  Re(A') )
C                                                                    (2)
C           ( -Im(B)  -Re(B)  | -Im(F)  -Re(F)  )
C           (                 |                 )
C           (  Re(B)  -Im(B)  |  Re(F)  -Im(F)  )
C           (                 |                 )
C     B_T = (-----------------+-----------------) ,  T = i*H,
C           (                 |                 )
C           ( -Im(G)  -Re(G)  | -Im(B')  Re(B') )
C           (                 |                 )
C           (  Re(G)  -Im(G)  | -Re(B') -Im(B') )
C
C   is determined and used to compute the eigenvalues. Optionally,
C   if compq = 1, an orthonormal basis of the right deflating
C   subspace of the pencil aS - bH, corresponding to the eigenvalues
C   with strictly negative real part, is computed. Namely, after
C   transforming aB_S - bB_H by unitary matrices, we have
C
C              ( Ao  Do  )              ( Bo  Fo  )
C     B_Sout = (         ) and B_Hout = (         ),                 (3)
C              (  0  Ao' )              (  0 -Bo' )
C
C   and the eigenvalues with strictly negative real part of the
C   complex pencil aB_Sout - bB_Hout are moved to the top.
C
C Input parameters:
C   A      - the m-by-m matrix A, with m = n/2.
C   DE     - an  m-by-(m+1) matrix containing the strict triangles of
C            the skew-Hermitian matrices D and E, as follows:
C            the leading m-by-m lower triangular part contains the lower
C            triangle of the matrix E, and the m-by-m upper triangular
C            part of the submatrix in the columns 2 to m+1 contains the
C            upper triangle of the matrix D of S in (1).
C            So, if i >= j, then E(i,j) = -E(j,i) is stored in DE(i,j)
C            and D(j,i) = -D(i,j) is stored in DE(j,i+1).
C            DE is an empty matrix if m = 0.
C   B      - the m-by-m matrix B.
C   FG     - an  m-by-(m+1) matrix containing the triangles of the
C            Hermitian matrices F and G, as follows:
C            the leading m-by-m lower triangular part contains the lower
C            triangle of the matrix G, and the m-by-m upper triangular
C            part of the submatrix in the columns 2 to m+1 contains the
C            upper triangle of the matrix F of H in (1).
C            So, if i >= j, then G(i,j) = -G(j,i) is stored in FG(i,j)
C            and F(j,i) = -F(i,j) is stored in FG(j,i+1).
C            FG is an empty matrix if m = 0.
C   compq  - (optional) scalar indicating whether an orthonormal basis Q
C            of the right deflating subspace is returned, or if Q is not
C            required, as follows:
C            = 0 :  Q is not required (default);
C            = 1 :  on exit, Q contains a matrix with orthonormal
C                   columns.
C   orthm  - (optional) if compq = 1, scalar indicating the technique
C            for computing the unitary basis of the deflating subspace,
C            as follows:
C            = 0 :  QR factorization, the fastest technique (default);
C            = 1 :  QR factorization with column pivoting;
C            = 2 :  singular value decomposition.
C            If compq = 0, the orthm value is not used.
C            Usually, orthm = 0 gives acceptable results, but badly
C            scaled or ill-conditioned problems might need to set
C            orthm = 1 or even orthm = 2.
C
C Output parameters:
C   ALPHA  - the n-vector of the numerators alpha defining the
C            eigenvalues of the pencil aS - bT.
C   BETA     the n-vector of the denominators beta defining the
C            eigenvalues of the pencil aS - bT.
C            Together, the quantities alpha = ALPHA(j) and
C            beta = BETA(j) represent the j-th eigenvalue of the pencil
C            aS - bT, in the form lambda = alpha/beta. Since lambda may
C            overflow, the ratios should not, in general, be computed.
C   Q      - if compq = 1, an n-by-neig matrix containing an orthonormal
C            basis of the right deflating subspace.
C   neig   - if compq = 1, the number of eigenvalues with negative real
C            parts.
C
C Contributors:
C   V. Sima, Research Institute for Informatics, Bucharest, Dec. 2010.
C
C Revisions:
C   V. Sima, Jan. 2011, Jul. 2011, Dec. 2011, Aug. 2012, Sep. 2012,
C   Nov. 2012, Jul. 2013, June 2014.
C   M. Voigt, Jul. 2013.
C
C     ******************************************************************
C
      SUBROUTINE MEXFUNCTION( NLHS, PLHS, NRHS, PRHS )
C
      INTEGER           DKIND
      PARAMETER         ( DKIND = KIND( 0.0D0 ) )
C
C     .. Mex-file interface parameters ..
      mwPointer         PLHS( * ), PRHS( * )
      INTEGER*4         NLHS, NRHS
C
C     .. Mex-file integer functions ..
      mwPointer         mxCalloc, mxCreateDoubleMatrix,
     $                  mxGetPi, mxGetPr
      INTEGER*4         mxGetM, mxGetN, mxIsNumeric, mxIsComplex
C
C     .. Scalar parameters used by subroutines ..
      CHARACTER         COMPQ, ORTH
      INTEGER           INFO, LDA, LDB, LDDE, LDFG, LDQ, LBWORK, LDWORK,
     $                  LZWORK, N, NEIG
C
C     .. Allocatable arrays ..
C     !Fortran 90/95 (Fixed dimensions should be used with Fortran 77.)
      LOGICAL,          ALLOCATABLE :: BWORK( : )
      INTEGER,          ALLOCATABLE :: IWORK( : )
      mwPointer         A, ALPHAI, ALPHAR, B, BETA, DE, DW, DWORK, FG,
     $                  Q, TMP, ZWORK
C
C     .. Local variables and constant dimension arrays ..
      CHARACTER*120     TEXT
      INTEGER           ICMPQ, IORTH, IP, LIWORK, LZWTMP, M, M1, MP1, N1
      DOUBLE PRECISION  DUM(  1 ), TEMP
      COMPLEX(DKIND) :: ZDUM( 2 )
C
C     .. External Subroutines ..
      EXTERNAL          ZLACPV, ZGHUDF, ZLACPY
C
C     ..Intrinsic Functions..
      INTRINSIC         INT, MAX, MIN
C
C     Check for proper number of arguments.
C
      IF( NRHS.LT.4 ) THEN
         CALL mexErrMsgTxt
     $        ( 'SKEWHAMILDEFLZ requires at least 4 input arguments.' )
      ELSE IF( NLHS.LT.2 ) THEN
         CALL mexErrMsgTxt
     $        ( 'SKEWHAMILDEFLZ requires at least 2 output arguments.' )
      END IF
C
C     Check input arguments.
C
C   A(mxm), DE(mx(m+1)), B(mxm), FG(mx(m+1))(, compq, orthm)) ).
C
      M  = mxGetM( PRHS( 1 ) )
      N1 = mxGetN( PRHS( 1 ) )
      IF ( MIN( M, N1 ).GT.0 ) THEN
         IF ( mxIsNumeric( PRHS(1) ).EQ.0 ) THEN
            CALL mexErrMsgTxt( 'A must be a real or complex matrix' )
         END IF
      END IF
      IF ( M.NE.N1 ) THEN
         CALL mexErrMsgTxt( 'A must be a square matrix' )
      END IF
      N   = 2*M
      MP1 = M + 1
C
      M1 = mxGetM( PRHS(2) )
      N1 = mxGetN( PRHS(2) )
      IF ( MIN( M1, N1 ).GT.1 ) THEN
         IF ( mxIsNumeric( PRHS(2) ).EQ.0 ) THEN
            CALL mexErrMsgTxt( 'DE must be a real or complex matrix' )
         END IF
      END IF
      IF ( M1.NE.M ) THEN
         CALL mexErrMsgTxt( 'DE must have the same number of rows as A'
     $                    )
      END IF
      IF ( N1.NE.MP1 ) THEN
         CALL mexErrMsgTxt( 'DE must have one more columns than rows' )
      END IF
C
      M1 = mxGetM( PRHS(3) )
      N1 = mxGetN( PRHS(3) )
      IF ( MIN( M1, N1 ).GT.0 ) THEN
         IF ( mxIsNumeric( PRHS(3) ).EQ.0 ) THEN
            CALL mexErrMsgTxt( 'B must be a real or complex matrix' )
         END IF
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
      IF ( MIN( M1, N1 ).GT.1 ) THEN
         IF ( mxIsNumeric( PRHS(4) ).EQ.0 ) THEN
            CALL mexErrMsgTxt( 'FG must be a real or complex matrix' )
         END IF
      END IF
      IF ( M1.NE.M ) THEN
         CALL mexErrMsgTxt( 'FG must have the same number of rows as A'
     $                    )
      END IF
      IF ( N1.NE.MP1 ) THEN
         CALL mexErrMsgTxt( 'FG must have one more columns than rows' )
      END IF
C
      IP = 5
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
         IF ( ICMPQ.LT.0 .OR. ICMPQ.GT.1 ) THEN
            CALL mexErrMsgTxt
     $           ( 'COMPQ has 0 or 1 the only admissible values' )
         END IF
         IP = IP + 1
C
      ELSE
         ICMPQ = 0
      END IF
C
      IF ( ICMPQ.EQ.0 ) THEN
         COMPQ = 'No Computation'
      ELSE
         COMPQ = 'Computation'
      END IF
C
      IF ( ICMPQ.EQ.1 .AND. NRHS.GE.IP ) THEN
         IF ( mxGetM( PRHS(IP) ).NE.1 .OR.
     $        mxGetN( PRHS(IP) ).NE.1 ) THEN
            CALL mexErrMsgTxt( 'ORTHM must be a scalar' )
         END IF
         IF ( mxIsNumeric( PRHS(IP) ).EQ.0 .OR.
     $        mxIsComplex( PRHS(IP) ).EQ.1 ) THEN
            CALL mexErrMsgTxt( 'ORTHM must be an integer scalar' )
         END IF
         CALL mxCopyPtrToReal8( mxGetPr( PRHS(IP) ), TEMP, 1 )
         IORTH = TEMP
         IF ( IORTH.LT.0 .OR. IORTH.GT.2 ) THEN
            CALL mexErrMsgTxt
     $           ( 'ORTHM has 0, 1, or 2 the only admissible values' )
         END IF
C
      ELSE
         IORTH = 0
      END IF
C
      IF ( IORTH.EQ.0 ) THEN
         ORTH = 'Q'
      ELSE IF ( IORTH.EQ.1 ) THEN
         ORTH = 'P'
      ELSE
         ORTH = 'S'
      END IF
C
C     Determine dimensions of the arrays and the workspace.
C
      LDA  = MAX( 1, N )
      LDB  = LDA
      LDDE = LDA
      LDFG = LDA
C
      IF ( ICMPQ.EQ.1 ) THEN
         N1  = 2*N
         LDQ = MAX( 1, N1 )
         LBWORK = N
      ELSE
         N1  = 0
         LDQ = 1
         LBWORK = 0
      END IF
C     
      IF ( IORTH.EQ.1 ) THEN
         LIWORK = N
      ELSE
         LIWORK = 0
      END IF
C
      LZWTMP = MAX( 1, M*MP1 )
C
C     Allocate part of the variable dimension local arrays.
C     !Fortran 90/95
C
      ALLOCATE( BWORK( LBWORK ), IWORK( LIWORK ) )
      A      = mxCalloc(  LDA*N, 16 )
      B      = mxCalloc(  LDB*N, 16 )
      DE     = mxCalloc( LDDE*N, 16 )
      FG     = mxCalloc( LDFG*N, 16 )
      ZWORK  = mxCalloc( LZWTMP, 16 )
C
C Copy inputs from MATLAB workspace to locally allocated arrays.
C
      M1 = MAX( 1, M )
      IF ( mxIsComplex( PRHS(1) ).EQ.0 ) THEN
         DW = mxCalloc( M*M, 8 )
         CALL mxCopyPtrToReal8( mxGetPr( PRHS( 1 ) ), %VAL( DW ),
     $                          M*M )
         CALL ZLACPV( M*M, 1, 1, %VAL( DW ), %VAL( ZWORK ) )
         CALL mxFree( DW )
      ELSE
         CALL mxCopyPtrToComplex16( mxGetPr( PRHS( 1 ) ),
     $                              mxGetPi( PRHS( 1 ) ),
     $                              %VAL( ZWORK ), M*M )
      END IF
      CALL ZLACPY( 'Full', M, M, %VAL( ZWORK ), M1, %VAL( A ), LDA )
C
      IF ( mxIsComplex( PRHS(2) ).EQ.0 ) THEN
         DW = mxCalloc( M*MP1, 8 )
         CALL mxCopyPtrToReal8( mxGetPr( PRHS( 2 ) ), %VAL( DW ),
     $                          M*MP1 )
         CALL ZLACPV( M*MP1, 1, 1, %VAL( DW ), %VAL( ZWORK ) )
         CALL mxFree( DW )
      ELSE
         CALL mxCopyPtrToComplex16( mxGetPr( PRHS( 2 ) ),
     $                              mxGetPi( PRHS( 2 ) ),
     $                              %VAL( ZWORK ), M*MP1 )
      END IF
      CALL ZLACPY( 'Full', M, MP1, %VAL( ZWORK ), M1, %VAL( DE ), LDDE )
C
      IF ( mxIsComplex( PRHS(3) ).EQ.0 ) THEN
         DW = mxCalloc( M*M, 8 )
         CALL mxCopyPtrToReal8( mxGetPr( PRHS( 3 ) ), %VAL( DW ),
     $                          M*M )
         CALL ZLACPV( M*M, 1, 1, %VAL( DW ), %VAL( ZWORK ) )
         CALL mxFree( DW )
      ELSE
         CALL mxCopyPtrToComplex16( mxGetPr( PRHS( 3 ) ),
     $                              mxGetPi( PRHS( 3 ) ),
     $                              %VAL( ZWORK ), M*M )
      END IF
      CALL ZLACPY( 'Full', M, M, %VAL( ZWORK ), M1, %VAL( B ), LDB )
C
      IF ( mxIsComplex( PRHS(4) ).EQ.0 ) THEN
         DW = mxCalloc( M*MP1, 8 )
         CALL mxCopyPtrToReal8( mxGetPr( PRHS( 4 ) ), %VAL( DW ),
     $                          M*MP1 )
         CALL ZLACPV( M*MP1, 1, 1, %VAL( DW ), %VAL( ZWORK ) )
         CALL mxFree( DW )
      ELSE
         CALL mxCopyPtrToComplex16( mxGetPr( PRHS( 4 ) ),
     $                              mxGetPi( PRHS( 4 ) ),
     $                              %VAL( ZWORK ), M*MP1 )
      END IF
      CALL ZLACPY( 'Full', M, MP1, %VAL( ZWORK ), M1, %VAL( FG ), LDFG )
      CALL mxFree( ZWORK )
C
C   ldwork, lzwork.
C
      CALL ZGHUDF( COMPQ, ORTH, N, %VAL( A ), LDA, %VAL( DE ), LDDE,
     $             %VAL( B ), LDB, %VAL( FG ), LDFG, NEIG, %VAL( A ),
     $             LDQ, DUM, DUM, DUM, IWORK, DUM, -1, ZDUM, -1, BWORK,
     $             INFO )
      LDWORK = INT(  DUM( 1 ) )
      LZWORK = INT( ZDUM( 1 ) )
C
C     Allocate the remaining variable dimension local arrays.
C     !Fortran 90/95
C
      ALPHAI = mxCalloc(      N,  8 )
      ALPHAR = mxCalloc(      N,  8 )
      BETA   = mxCalloc(      N,  8 )
      DWORK  = mxCalloc( LDWORK,  8 )
      Q      = mxCalloc( LDQ*N1, 16 )
      ZWORK  = mxCalloc( LZWORK, 16 )
C
C Do the actual computations.
C
      CALL ZGHUDF( COMPQ, ORTH, N, %VAL( A ), LDA, %VAL( DE ), LDDE,
     $             %VAL( B ), LDB, %VAL( FG ), LDFG, NEIG, %VAL( Q ),
     $             LDQ, %VAL( ALPHAR ), %VAL( ALPHAI ), %VAL( BETA ),
     $             IWORK, %VAL( DWORK ), LDWORK, %VAL( ZWORK ), LZWORK,
     $             BWORK, INFO )
C
C Copy output to MATLAB workspace.
C
      IF ( INFO.EQ.0 ) THEN
         PLHS( 1 ) = mxCreateDoubleMatrix( N, 1, 1 )
         CALL mxCopyReal8ToPtr( %VAL( ALPHAR ),
     $                          mxGetPr( PLHS( 1 ) ), N )
         CALL mxCopyReal8ToPtr( %VAL( ALPHAI ),
     $                          mxGetPi( PLHS( 1 ) ), N )
C
         PLHS( 2 ) = mxCreateDoubleMatrix( N, 1, 0 )
         CALL mxCopyReal8ToPtr( %VAL( BETA ),
     $                          mxGetPr( PLHS( 2 ) ), N )
C
         IP = 3
         IF( ICMPQ.GT.0 .AND. NLHS.GE.IP ) THEN
            TMP = mxCalloc( N*NEIG, 16 )
            CALL ZLACPY( 'Full', N, NEIG, %VAL( Q ), LDQ, %VAL( TMP ),
     $                   LDA )
            PLHS( IP ) = mxCreateDoubleMatrix( N, NEIG, 1 )
            CALL mxCopyComplex16ToPtr( %VAL( TMP ),
     $                                 mxGetPr( PLHS( IP ) ),
     $                                 mxGetPi( PLHS( IP ) ),
     $                                 N*NEIG )
            CALL mxFree( TMP )
            IP = IP + 1
C
            IF( NLHS.GE.IP ) THEN
               PLHS( IP ) = mxCreateDoubleMatrix( 1, 1, 0 )
               TEMP       = NEIG
               CALL mxCopyReal8ToPtr( TEMP, mxGetPr( PLHS( IP ) ),
     $                                1 )
            END IF
         END IF
      END IF
C
C     Deallocate variable dimension arrays.
C     !Fortran 90/95
C
      DEALLOCATE(  BWORK, IWORK )
      CALL mxFree( A      )
      CALL mxFree( ALPHAI )
      CALL mxFree( ALPHAR )
      CALL mxFree( B      )
      CALL mxFree( BETA   )
      CALL mxFree( DE     )
      CALL mxFree( DWORK  )
      CALL mxFree( FG     )
      CALL mxFree( Q      )
      CALL mxFree( ZWORK  )
C
C Error and warning handling.
C
      IF ( INFO.NE.0 ) THEN
         WRITE( TEXT, '('' INFO = '',I4,'' ON EXIT FROM ZGHUDF'')'
     $        ) INFO
         CALL mexErrMsgTxt( TEXT )
      END IF
C
      RETURN
C *** Last line of skewHamildeflZ ***
      END
