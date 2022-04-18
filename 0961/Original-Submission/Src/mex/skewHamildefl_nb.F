#include "fintrf.h"
C SKEWHAMILDEFL_NB.f - Gateway function for computing the eigenvalues of
C                   a skew-Hamiltonian/Hamiltonian pencil and the right
C                   deflating subspace corresponding to the eigenvalues
C                   with strictly negative real part, using the routine
C                   DGHUDP.
C                   The block size, NB, can be specified.
C
C Matlab call:
C   [ALPHAR,ALPHAI,BETA(,Q,neig)] = ...
C                  skewHamildefl_nb(A,DE,B,FG(,nb,compq(,orthm)))
C
C Purpose:
C   To compute the eigenvalues of an n-by-n skew-Hamiltonian/Hamiltonian
C   pencil aS - bH, with
C
C         (  A  D  )         (  B  F  )
C     S = (        ) and H = (        ),                              (1)
C         (  E  A' )         (  G -B' )
C
C   using the structured Schur form of an embedded pencil. Here, M'
C   denotes the transpose of the matrix M. Optionally, the right
C   deflating subspace of aS - bH corresponding to the eigenvalues with
C   strictly negative real part is computed.
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
C            symmetric matrices F and G, as follows:
C            the leading m-by-m lower triangular part contains the lower
C            triangle of the matrix G, and the m-by-m upper triangular
C            part of the submatrix in the columns 2 to m+1 contains the
C            upper triangle of the matrix F of H in (1).
C            So, if i >= j, then G(i,j) = G(j,i) is stored in FG(i,j)
C            and F(j,i) = F(i,j) is stored in FG(j,i+1).
C            FG is an empty matrix if m = 0.
C   nb     - (optional) scalar indicating the block size. An internally
C            computed block size will be used for a non-positive value.
C            Default:  nb = 32.
C   compq  - (optional) scalar indicating whether the orthogonal
C            transformation matrix Q is returned, or if Q is not
C            required, as follows:
C            = 0 :  Q is not required (default);
C            = 1 :  on exit, Q contains a matrix Q with orthogonal
C                   columns.
C   orthm  - (optional) if compq = 1, scalar indicating the technique
C            for computing the orthogonal basis of the deflating
C            subspace, as follows:
C            = 0 :  QR factorization, the fastest technique (default);
C            = 1 :  QR factorization with column pivoting;
C            = 2 :  singular value decomposition.
C            If compq = 0, the orthm value is not used.
C            Usually, orthm = 0 gives acceptable results, but badly
C            scaled or ill-conditioned problems might need to set
C            orthm = 1 or even orthm = 2.
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
C            eigenvalue. Only eigenvalues with imaginary parts greater
C            than or equal to zero are stored; their conjugate
C            eigenvalues are not stored. If imaginary parts are zero
C            (i.e., for real eigenvalues), only positive eigenvalues
C            are stored.
C   Q      - if compq = 1, an n-by-neig matrix containing the computed
C            basis of the right deflating subspace.
C   neig   - if compq = 1, the number of eigenvalues with negative real
C            parts.
C
C Contributors:
C   V. Sima, Research Institute for Informatics, Bucharest, Oct. 2010.
C
C Revisions:
C   V. Sima, Nov. 2010, Dec. 2010, Mar. 2011, July 2013, July 2014.
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
      CHARACTER         COMPQ, ORTH
      INTEGER           INFO, LDA, LDB, LDDE, LDFG, LDQ, LDWORK, LIWORK,
     $                  N, NEIG
C
C     .. Allocatable arrays ..
C     !Fortran 90/95 (Fixed dimensions should be used with Fortran 77.)
      LOGICAL,          ALLOCATABLE :: BWORK( : )
      INTEGER,          ALLOCATABLE :: IWORK( : )
      mwPointer         A, ALPHAI, ALPHAR, B, BETA, DE, DWORK, FG, Q
C
C     .. Local variables and constant dimension arrays ..
      CHARACTER*120     TEXT
      INTEGER           ICMPQ, IORTH, IP, M, M1, N1, N2, NB
      DOUBLE PRECISION  DUM( 1 ), TEMP
      REAL              TIME
      INTEGER           TIC, TOC, RATE
C
C     .. External Subroutines ..
      EXTERNAL          DLACPY, DGHUDP
C
C     ..Intrinsic Functions..
      INTRINSIC         INT, MAX
C
C     Check for proper number of arguments.
C 
      IF( NRHS.LT.4 ) THEN
         CALL mexErrMsgTxt
     $        ( 'SKEWHAMILDEFL requires at least 4 input arguments.' )
      ELSE IF ( NLHS.LT.3 ) THEN
         CALL mexErrMsgTxt
     $        ( 'SKEWHAMILDEFL requires at least 3 output arguments.' )
      END IF
C
C   A(mxm), DE(mx(m+1)), B(mxm), FG(mx(m+1))(, nb, compq, orthm).
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
      N2 = 2*N
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
         COMPQ = 'N'
      ELSE
         COMPQ = 'C'
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
C Determine dimensions of the arrays and the workspace.
C
      LDA  = MAX( 1, M )
      LDDE = LDA
      LDB  = LDA
      LDFG = LDA
      IF ( ICMPQ.EQ.0 ) THEN
         LDQ = 1
         N1  = 0
      ELSE
         LDQ = MAX( 1, N2 )
         N1  = N2
      END IF
      LIWORK = MAX( 32, N + 12, N2 + 1 )
C
C Allocate variable dimension local arrays.
C !Fortran 90/95
C
      ALLOCATE( BWORK( M ), IWORK( LIWORK ) )
      A      = mxCalloc(  LDA*M,     8 )
      ALPHAI = mxCalloc(      M,     8 )
      ALPHAR = mxCalloc(      M,     8 )
      B      = mxCalloc(  LDB*M,     8 )
      BETA   = mxCalloc(      M,     8 )
      DE     = mxCalloc( LDDE*(M+1), 8 )
      FG     = mxCalloc( LDFG*(M+1), 8 )
      Q      = mxCalloc( LDQ*N1,     8 )
C
C   ldwork
C
      CALL DGHUDP( COMPQ, ORTH, N, %VAL( A ), LDA, %VAL( DE ), LDDE,
     $             %VAL( B ), LDB, %VAL( FG ), LDFG, NEIG, %VAL( Q ),
     $             LDQ, %VAL( ALPHAR ), %VAL( ALPHAI ), %VAL( BETA ),
     $             IWORK, LIWORK, DUM, -1, BWORK, INFO )
      LDWORK = INT( DUM( 1 ) )
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
C
C Do the actual computations.
C
      CALL SYSTEM_CLOCK( TIC, RATE )
      INFO = NB
      CALL DGHUDP( COMPQ, ORTH, N, %VAL( A ), LDA, %VAL( DE ), LDDE,
     $             %VAL( B ), LDB, %VAL( FG ), LDFG, NEIG, %VAL( Q ),
     $             LDQ, %VAL( ALPHAR ), %VAL( ALPHAI ), %VAL( BETA ),
     $             IWORK, LIWORK, %VAL( DWORK ), LDWORK, BWORK, INFO )
      CALL SYSTEM_CLOCK( TOC )
      TIME = REAL( TOC-TIC )/RATE
**    WRITE( TEXT, '('' TIME = '',F20.12,'' ON EXIT FROM BLK SOLVER'')'
**   $      ) TIME
**    CALL mexPrintf( TEXT )
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
         CALL mxCopyReal8ToPtr( %VAL( BETA ),
     $                          mxGetPr( PLHS( 3 ) ), M )
C
         IP = 4
         IF ( NLHS.GE.IP .AND. ICMPQ.GE.1 ) THEN
            CALL DLACPY( 'Full', N, NEIG, %VAL( Q ), LDQ, %VAL( DWORK ),
     $                    MAX( 1, N ) )
            PLHS( IP ) = mxCreateDoubleMatrix( N, NEIG, 0 )
            CALL mxCopyReal8ToPtr( %VAL( DWORK ),
     $                             mxGetPr( PLHS( IP ) ), N*NEIG )
            IP = IP + 1
C
            IF ( NLHS.GE.IP ) THEN
               PLHS( IP ) = mxCreateDoubleMatrix( 1, 1, 0 )
               TEMP = NEIG
               CALL mxCopyReal8ToPtr( TEMP,
     $                                mxGetPr( PLHS( IP ) ), 1 )
            END IF
         END IF
      END IF
C
C Deallocate variable dimension arrays.
C !Fortran 90/95
C
      DEALLOCATE( BWORK, IWORK )
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
         WRITE( TEXT, '('' INFO = '',I4,'' ON EXIT FROM DGHUDP'')'
     $        ) INFO
         CALL mexErrMsgTxt( TEXT )
      END IF
C
      RETURN
C *** Last line of skewHamildefl_nb ***
      END
