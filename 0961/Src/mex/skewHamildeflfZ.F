#include "fintrf.h"
C SKEWHAMILDEFLFZ.F - Gateway function for computing the eigenvalues of
C                     a complex skew-Hamiltonian/Hamiltonian pencil and
C                     the right deflating subspace and/or the companion
C                     subspace corresponding to the eigenvalues with
C                     strictly negative real part, using the routine
C                     ZGHFDF.
C                     The gateway accepts real or complex input
C                     matrices.
C
C Matlab call:
C   [ALPHA,BETA(,Q,U,neig)] = ...
C                skewHamildeflfZ(Z,B,FG(,compq,compu(,orthm)))
C
C Purpose:
C   To compute the eigenvalues of a complex n-by-n skew-Hamiltonian/
C   Hamiltonian pencil aS - bH, with n = 2m,
C
C                             (  B  F  )       (  Z11  Z12  )
C       S = J Z' J' Z and H = (        ), Z =: (            ),       (1)
C                             (  G -B' )       (  Z21  Z22  )
C
C   where M' denotes the conjugate transpose of the matrix M.
C   The structured Schur form of the embedded real skew-Hamiltonian/
C   skew-Hamiltonian pencil, aB_S - bB_T, with B_S = J B_Z' J' B_Z,
C
C             (  Re(Z11)  -Im(Z11)  |  Re(Z12)  -Im(Z12)  )
C             (                     |                     )
C             (  Im(Z11)   Re(Z11)  |  Im(Z12)   Re(Z12)  )
C             (                     |                     )
C       B_Z = (---------------------+---------------------) ,
C             (                     |                     )
C             (  Re(Z21)  -Im(Z21)  |  Re(Z22)  -Im(Z22)  )
C             (                     |                     )
C             (  Im(Z21)   Re(Z21)  |  Im(Z22)   Re(Z22)  )
C                                                                    (2)
C             ( -Im(B)  -Re(B)  | -Im(F)  -Re(F)  )
C             (                 |                 )
C             (  Re(B)  -Im(B)  |  Re(F)  -Im(F)  )
C             (                 |                 )
C       B_T = (-----------------+-----------------) ,  T = i*H,
C             (                 |                 )
C             ( -Im(G)  -Re(G)  | -Im(B')  Re(B') )
C             (                 |                 )
C             (  Re(G)  -Im(G)  | -Re(B') -Im(B') )
C
C   is determined and used to compute the eigenvalues. Optionally, if
C   compq = 1, an orthonormal basis of the right deflating subspace,
C   Def_-(S, H), of the pencil aS - bH in (1), corresponding to the
C   eigenvalues with strictly negative real part, is computed. Namely,
C   after transforming aB_S - bB_H, in the factored form, by unitary
C   matrices, we have B_Sout = J B_Zout' J' B_Zout,
C
C                ( Ao  Do  )              ( Bo  Fo  )
C       B_Zout = (         ) and B_Hout = (         ),               (3)
C                (  0  Co  )              (  0 -Bo' )
C
C   and the eigenvalues with strictly negative real part of the
C   complex pencil aB_Sout - bB_Hout are moved to the top. Optionally,
C   if compu = 1, an orthonormal basis of the companion subspace,
C   range(P_U), which corresponds to the eigenvalues with negative
C   real part, is computed.
C
C Input parameters:
C   Z      - the n-by-n matrix Z.
C   B      - the m-by-m matrix B, with m = n/2.
C   FG     - an  m-by-(m+1) matrix containing the triangles of the
C            Hermitian matrices F and G, as follows:
C            the leading m-by-m lower triangular part contains the lower
C            triangle of the matrix G, and the m-by-m upper triangular
C            part of the submatrix in the columns 2 to m+1 contains the
C            upper triangle of the matrix F of H in (1).
C            So, if i >= j, then G(i,j) = -G(j,i) is stored in FG(i,j)
C            and F(j,i) = -F(i,j) is stored in FG(j,i+1).
C            FG is an empty matrix if m = 0.
C   compq  - (optional) scalar indicating whether to compute a basis Q
C            of the deflating subspace corresponding to the eigenvalues
C            of the pencil aS - bH with strictly negative real part.
C            = 0 :  Q is not required (default);
C            = 1 :  on exit, Q contains a matrix with orthonormal
C                   columns.
C   compu  - (optional) scalar indicating whether to compute a basis U
C            of the companion subspace corresponding to the eigenvalues
C            of the pencil aS - bH with strictly negative real part.
C            = 0 :  U is not required (default);
C            = 1 :  on exit, U contains a matrix with orthonormal
C                   columns.
C   orthm  - (optional) if compq = 1 or compu = 1, scalar indicating the
C            technique for computing the unitary basis of the deflating
C            subspace, or the companion subspace, as follows:
C            = 0 :  QR factorization, the fastest technique (default);
C            = 1 :  QR factorization with column pivoting;
C            = 2 :  singular value decomposition.
C            If compq = 0 and compu = 0, the orthm value is not used.
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
C   U      - if compu = 1, an n-by-neig matrix containing an orthonormal
C            basis of the companion subspace.
C   neig   - if compq = 1 or compu = 1, the number of eigenvalues with
C            negative real parts.
C
C Contributors:
C   V. Sima, Research Institute for Informatics, Bucharest, Jan. 2011.
C
C Revisions:
C   V. Sima, July 2011, Jul 2012, Sep. 2012, Nov. 2012, Jul 2013,
C   June 2014.
C   M. Voigt, Jan. 2012, Jul. 2013.
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
     $                  mxGetPr, mxGetPi
      INTEGER*4         mxGetM, mxGetN, mxIsNumeric, mxIsComplex
C
C     .. Scalar parameters used by subroutines ..
      CHARACTER         COMPQ, COMPU, ORTH
      INTEGER           INFO, LDB, LDC, LDD, LDFG, LDQ, LDU, LBWORK,
     $                  LDWORK, LDZ, LIWORK, LZWORK, N, NEIG
C
C     .. Allocatable arrays ..
C     !Fortran 90/95 (Fixed dimensions should be used with Fortran 77.)
      LOGICAL,          ALLOCATABLE :: BWORK( : )
      INTEGER,          ALLOCATABLE :: IWORK( : )
      mwPointer         ALPHAI, ALPHAR, B, BETA, C, D, DW, DWORK, FG, Q,
     $                  U, TMP, Z, ZWORK
C
C     .. Local variables and constant dimension arrays ..
      CHARACTER*120     TEXT
      INTEGER           ICMPQ, ICMPU, IORTH, IP, M, M1, MP1, N1, N2, N3
      DOUBLE PRECISION  DUM(  1 ), TEMP
      COMPLEX(DKIND) :: ZDUM( 2 )
C
C     .. External Subroutines ..
      EXTERNAL          ZLACPV, ZGHFDF, ZLACPY
C
C     ..Intrinsic Functions..
      INTRINSIC         INT, MAX, MIN, MOD
C
C     Check for proper number of arguments.
C 
      IF( NRHS.LT.3 ) THEN
         CALL mexErrMsgTxt
     $      ( 'SKEWHAMILDEFLFZ requires at least 3 input arguments.' )
      ELSE IF( NLHS.LT.2 ) THEN
         CALL mexErrMsgTxt
     $      ( 'SKEWHAMILDEFLFZ requires at least 2 output arguments.' )
      END IF
C
C     Check input arguments.
C
C   Z(nxn), B(mxm), FG(mx(m+1))(, compq, compu, orthm)) ).
C
      N  = mxGetM( PRHS( 1 ) )
      N1 = mxGetN( PRHS( 1 ) )
      IF ( MIN( N, N1 ).GT.0 ) THEN
         IF ( mxIsNumeric( PRHS(1) ).EQ.0 ) THEN
            CALL mexErrMsgTxt( 'Z must be a complex matrix' )
         END IF
      END IF
      IF ( N.NE.N1 ) THEN
         CALL mexErrMsgTxt( 'Z must be a square matrix' )
      END IF
      IF ( MOD( N, 2 ).NE.0 ) THEN
         CALL mexErrMsgTxt( 'Z must have an even order' )
      END IF
      M   = N/2
      MP1 = M + 1
C
      M1 = mxGetM( PRHS(2) )
      N1 = mxGetN( PRHS(2) )
      IF ( MIN( M1, N1 ).GT.0 ) THEN
         IF ( mxIsNumeric( PRHS(2) ).EQ.0 ) THEN
            CALL mexErrMsgTxt( 'B must be a complex matrix' )
         END IF
      END IF
      IF ( M1.NE.M ) THEN
         CALL mexErrMsgTxt( 'B must have half the number of rows of Z'
     $                    )
      END IF
      IF ( N1.NE.M ) THEN
         CALL mexErrMsgTxt( 'B must be a square matrix' )
      END IF
C
      M1 = mxGetM( PRHS(3) )
      N1 = mxGetN( PRHS(3) )
      IF ( MIN( M1, N1 ).GT.1 ) THEN
         IF ( mxIsNumeric( PRHS(3) ).EQ.0 ) THEN
            CALL mexErrMsgTxt( 'FG must be a complex matrix' )
         END IF
      END IF
      IF ( M1.NE.M ) THEN
         CALL mexErrMsgTxt( 'FG must have the same number of rows as B'
     $                    )
      END IF
      IF ( N1.NE.MP1 ) THEN
         CALL mexErrMsgTxt( 'FG must have one more columns than rows' )
      END IF
C
      IP = 4
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
         IF ( ICMPU.LT.0 .OR. ICMPU.GT.1 ) THEN
            CALL mexErrMsgTxt
     $           ( 'COMPU has 0 or 1 the only admissible values' )
         END IF
         IP = IP + 1
C
      ELSE
         ICMPU = 0
      END IF
C
      IF ( ICMPU.EQ.0 ) THEN
         COMPU = 'No Computation'
      ELSE
         COMPU = 'Computation'
      END IF
C
      IF ( MAX( ICMPQ, ICMPU ).GT.0 .AND. NRHS.GE.IP ) THEN
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
      LDZ  = MAX( 1, N )
      LDB  = LDZ
      LDFG = LDZ
      IF ( ICMPQ.GT.0 ) THEN
         N1  = 2*N
         LDQ = MAX( 1, N1 )
      ELSE
         N1  = 0
         LDQ = 1
      END IF
C
      IF ( ICMPU.GT.0 ) THEN
         LDU = MAX( 1, N )
         N2  = 2*N
      ELSE
         LDU = 1
         N2  = 0
      END IF
C
      IF ( MAX( ICMPQ, ICMPU ).GT.0 ) THEN
         LDC = LDZ
         LDD = LDZ
         N3  = N
         LBWORK = N
      ELSE
         LDC = 1
         LDD = 1
         N3  = 0
         LBWORK = 0
      END IF
C
      LIWORK = 2*N + 9
C
C     Allocate part of the variable dimension local arrays.
C     !Fortran 90/95
C
      ALLOCATE( BWORK( LBWORK ), IWORK( LIWORK ) )
      B     = mxCalloc(  LDB*N, 16 )
      FG    = mxCalloc( LDFG*N, 16 )
      Z     = mxCalloc(  LDZ*N, 16 )
C
C Copy inputs from MATLAB workspace to locally allocated arrays.
C
      M1 = MAX( 1, M )
      IF ( mxIsComplex( PRHS( 1 ) ).EQ.0 ) THEN
         DW = mxCalloc( N*N, 8 )
         CALL mxCopyPtrToReal8( mxGetPr( PRHS( 1 ) ), %VAL( DW ),
     $                          N*N )
         CALL ZLACPV( N*N, 1, 1, %VAL( DW ), %VAL( Z ) )
         CALL mxFree( DW )
      ELSE
         CALL mxCopyPtrToComplex16( mxGetPr( PRHS( 1 ) ),
     $                              mxGetPi( PRHS( 1 ) ),
     $                              %VAL( Z ), N*N )
      END IF
C
      TMP = mxCalloc( M*MP1, 16 )
      IF ( mxIsComplex( PRHS( 2 ) ).EQ.0 ) THEN
         DW = mxCalloc( M*M, 8 )
         CALL mxCopyPtrToReal8( mxGetPr( PRHS( 2 ) ), %VAL( DW ),
     $                          M*M )
         CALL ZLACPV( M*M, 1, 1, %VAL( DW ), %VAL( TMP ) )
         CALL mxFree( DW )
      ELSE
         CALL mxCopyPtrToComplex16( mxGetPr( PRHS( 2 ) ),
     $                              mxGetPi( PRHS( 2 ) ),
     $                              %VAL( TMP ), M*M )
      END IF
      CALL ZLACPY( 'Full', M, M, %VAL( TMP ), M1, %VAL( B ), LDB )
C
      IF ( mxIsComplex( PRHS( 3 ) ).EQ.0 ) THEN
         DW = mxCalloc( M*MP1, 8 )
         CALL mxCopyPtrToReal8( mxGetPr( PRHS( 3 ) ), %VAL( DW ),
     $                          M*MP1 )
         CALL ZLACPV( M*MP1, 1, 1, %VAL( DW ), %VAL( TMP ) )
         CALL mxFree( DW )
      ELSE
         CALL mxCopyPtrToComplex16( mxGetPr( PRHS( 3 ) ),
     $                              mxGetPi( PRHS( 3 ) ),
     $                              %VAL( TMP ), M*MP1 )
      END IF
      CALL ZLACPY( 'Full', M, MP1, %VAL( TMP ), M1, %VAL( FG ), LDFG )
      CALL mxFree( TMP )
C
C   ldwork, lzwork.
C
      CALL ZGHFDF( COMPQ, COMPU, ORTH, N, %VAL( Z ), LDZ, %VAL( B ),
     $             LDB, %VAL( FG ), LDFG, NEIG, %VAL( Z ), LDD,
     $             %VAL( Z ), LDC, %VAL( Z ), LDQ, %VAL( Z ), LDU, DUM,
     $             DUM, DUM, IWORK, LIWORK, DUM, -1, ZDUM, -1, BWORK,
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
      C      = mxCalloc( LDC*N3, 16 )
      D      = mxCalloc( LDD*N3, 16 )
      DWORK  = mxCalloc( LDWORK,  8 )
      Q      = mxCalloc( LDQ*N1, 16 )
      U      = mxCalloc( LDU*N2, 16 )
      ZWORK  = mxCalloc( LZWORK, 16 )
C
C Do the actual computations.
C
      CALL ZGHFDF( COMPQ, COMPU, ORTH, N, %VAL( Z ), LDZ, %VAL( B ),
     $             LDB, %VAL( FG ), LDFG, NEIG, %VAL( D ), LDD,
     $             %VAL( C ), LDC, %VAL( Q ), LDQ, %VAL( U ), LDU,
     $             %VAL( ALPHAR ), %VAL( ALPHAI ), %VAL( BETA ), IWORK,
     $             LIWORK, %VAL( DWORK ), LDWORK, %VAL( ZWORK ), LZWORK,
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
         IF( MAX( ICMPQ, ICMPU ).GT.0 .AND. NLHS.GE.IP ) THEN
            TMP = mxCalloc( N*NEIG, 16 )
            IF( ICMPQ.GT.0 ) THEN
               CALL ZLACPY( 'Full', N, NEIG, %VAL( Q ), LDQ,
     $                      %VAL( TMP ), LDZ )
               PLHS( IP ) = mxCreateDoubleMatrix( N, NEIG, 1 )
               CALL mxCopyComplex16ToPtr( %VAL( TMP ),
     $                                    mxGetPr( PLHS( IP ) ),
     $                                    mxGetPi( PLHS( IP ) ),
     $                                    N*NEIG )
               IP = IP + 1
               IF( ICMPU.EQ.0 )
     $            CALL mxFree( TMP )
            END IF
C
            IF( ICMPU.GT.0 ) THEN
               CALL ZLACPY( 'Full', N, NEIG, %VAL( U ), LDU,
     $                      %VAL( TMP ), LDZ )
               PLHS( IP ) = mxCreateDoubleMatrix( N, NEIG, 1 )
               CALL mxCopyComplex16ToPtr( %VAL( TMP ),
     $                                    mxGetPr( PLHS( IP ) ),
     $                                    mxGetPi( PLHS( IP ) ),
     $                                    N*NEIG )
               CALL mxFree( TMP )
               IP = IP + 1
            END IF
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
      DEALLOCATE(  BWORK, IWORK   )
      CALL mxFree( ALPHAI )
      CALL mxFree( ALPHAR )
      CALL mxFree( B      )
      CALL mxFree( BETA   )
      CALL mxFree( C      )
      CALL mxFree( D      )
      CALL mxFree( DWORK  )
      CALL mxFree( FG     )
      CALL mxFree( Q      )
      CALL mxFree( U      )
      CALL mxFree( Z      )
      CALL mxFree( ZWORK  )
C
C Error and warning handling.
C
      IF ( INFO.NE.0 ) THEN
         WRITE( TEXT, '('' INFO = '',I4,'' ON EXIT FROM ZGHFDF'')'
     $        ) INFO
         CALL mexErrMsgTxt( TEXT )
      END IF
C
      RETURN
C *** Last line of skewHamildeflfZ ***
      END
