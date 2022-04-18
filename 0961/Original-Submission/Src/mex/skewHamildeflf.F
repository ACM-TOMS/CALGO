#include "fintrf.h"
C SKEWHAMILDEFLF.F - Gateway function for computing the eigenvalues of a
C                    skew-Hamiltonian/Hamiltonian pencil, the right
C                    deflating subspace and the companion subspace
C                    corresponding to the eigenvalues with strictly
C                    negative real part, using the routine DGHFDF.
C
C Matlab call:
C   [ALPHAR,ALPHAI,BETA(,Q,U,neig)] = skewHamildeflf(Z,B,FG(,...
C                                              compq,compu(,orthm)))
C
C Purpose:
C   To compute the relevant eigenvalues of a real n-by-n skew-
C   Hamiltonian/Hamiltonian pencil aS - bH, with
C
C                                 (  B  F  )      (  0  I  )
C     S = T Z = J Z' J' Z and H = (        ), J = (        ),        (1)
C                                 (  G -B' )      ( -I  0  )
C
C   where the notation M' denotes the transpose of the matrix M.
C   Optionally, if compq = 1, an orthogonal basis of the right deflating
C   subspace of aS - bH corresponding to the eigenvalues with strictly
C   negative real part is computed. Optionally, if compu = 1, an
C   orthogonal basis of the companion subspace, range(P_U), which
C   corresponds to the eigenvalues with strictly negative real part, is
C   computed.
C
C Input parameters:
C   Z      - the n-by-n matrix Z.
C   B      - the m-by-m matrix B, with m = n/2.
C   FG     - an  m-by-(m+1) matrix containing the triangles of the
C            symmetric matrices F and G, as follows:
C            the leading m-by-m lower triangular part contains the lower
C            triangle of the matrix G, and the m-by-m upper triangular
C            part of the submatrix in the columns 2 to m+1 contains the
C            upper triangle of the matrix F of H in (1).
C            So, if i >= j, then G(i,j) = G(j,i) is stored in FG(i,j)
C            and F(j,i) = F(i,j) is stored in FG(j,i+1).
C            FG is an empty matrix if m = 0.
C   compq  - (optional) scalar indicating whether the orthogonal basis Q
C            of the right deflating subspace corresponding to the
C            eigenvalues with strictly negative real part is returned,
C            or if Q is not required, as follows:
C            = 0 :  Q is not required (default);
C            = 1 :  on exit, Q contains a matrix Q with orthogonal
C                   columns.
C   compu  - (optional) scalar indicating whether the orthogonal basis U
C            of the companion subspace corresponding to the eigenvalues
C            with strictly negative real part is returned, or if U is
C            not required, as follows:
C            = 0 :  U is not required (default);
C            = 1 :  on exit, U contains a matrix U with orthogonal
C                   columns.
C   orthm  - (optional) if compq = 1 or compu = 1, scalar indicating the
C            technique for computing the orthogonal basis of the
C            deflating subspace, or the companion subspace, as follows:
C            = 0 :  QR factorization, the fastest technique (default);
C            = 1 :  QR factorization with column pivoting;
C            = 2 :  singular value decomposition.
C            If compq = 0 and compu = 0, the orthm value is not used.
C            Usually, orthm = 0 gives acceptable results, but badly
C            scaled or ill-conditioned problems might need to set
C            orthm = 1 or even orthm = 2.
C
C Output parameters:
C   ALPHAR,- the m-vectors of real parts and imaginary parts,
C   ALPHAI   respectively, of each scalar alpha defining an eigenvalue
C            of the pencil aS - bH.
C            If ALPHAI(j) is zero, then the j-th eigenvalue is real.
C   BETA     the m-vector of the scalars beta that define the
C            eigenvalues of the pencil aS - bH.
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
C   U      - if compu = 1, an n-by-neig matrix containing the computed
C            basis of the companion subspace.
C   neig   - if compq = 1 or compu = 1, the number of eigenvalues with
C            negative real parts.
C
C Contributors:
C   V. Sima, Research Institute for Informatics, Bucharest, Feb. 2011.
C
C Revisions:
C   V. Sima, July 2013, July 2014.
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
      INTEGER*4         mxGetM, mxGetN, mxIsComplex, mxIsNumeric
C
C     .. Scalar parameters used by SLICOT subroutines ..
      CHARACTER         COMPQ, COMPU, ORTH
      INTEGER           INFO, IWARN, LDB, LDFG, LDQ, LDU, LDWORK, LDZ,
     $                  LIWORK, N, NEIG
C
C     .. Allocatable arrays ..
C     !Fortran 90/95 (Fixed dimensions should be used with Fortran 77.)
      LOGICAL,          ALLOCATABLE :: BWORK( : )
      INTEGER,          ALLOCATABLE :: IWORK( : )
      DOUBLE PRECISION, ALLOCATABLE :: DW( : )
      mwPointer         ALPHAI, ALPHAR, B, BETA, DWORK, FG, Q, U, Z
C
C     .. Local variables and constant dimension arrays ..
      CHARACTER*120     TEXT
      INTEGER           I, ICMPQ, ICMPU, IORTH, IP, M, M1, N1, N2, NU
      DOUBLE PRECISION  DUM( 2 ), TEMP
C
C     .. External Subroutines ..
      EXTERNAL          DCOPY, DLACPY, DGHFDF
C
C     ..Intrinsic Functions..
      INTRINSIC         INT, MAX, MOD
C
C     Check for proper number of arguments.
C 
      IF( NRHS.LT.3 ) THEN
         CALL mexErrMsgTxt
     $        ( 'SKEWHAMILDEFLF requires at least 3 input arguments.' )
      ELSE IF ( NLHS.LT.3 ) THEN
         CALL mexErrMsgTxt
     $        ( 'SKEWHAMILDEFLF requires at least 3 output arguments.' )
      END IF
C
C   Z(nxn), B(mxm), FG(mx(m+1))(, compq, compu, orthm).
C
      N  = mxGetM( PRHS(1) )
      N1 = mxGetN( PRHS(1) )
      IF ( N.LT.0 ) THEN
         CALL mexErrMsgTxt( 'N must be non-negative' )
      ELSE IF ( MOD( N, 2 ).NE.0 ) THEN
         CALL mexErrMsgTxt( 'N must be even' )
      END IF
C
      IF ( mxIsNumeric( PRHS(1) ).EQ.0 .OR.
     $     mxIsComplex( PRHS(1) ).EQ.1 ) THEN
         CALL mexErrMsgTxt( 'Z must be a real matrix' )
      END IF
      IF ( N.NE.N1 ) THEN
         CALL mexErrMsgTxt( 'Z must be a square matrix' )
      END IF
      M  = N/2
      N2 = 2*N
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
         CALL mexErrMsgTxt( 'FG must have the same number of rows as B'
     $                    )
      END IF
      IF ( N1.NE.M+1 ) THEN
         CALL mexErrMsgTxt( 'FG must have one more columns than rows' )
      END IF
C
      IP = 4
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
         COMPU = 'N'
      ELSE
         COMPU = 'C'
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
C Determine dimensions of the arrays and the workspace.
C
      LDZ  = MAX( 1, N )
      LDB  = MAX( 1, M )
      LDFG = LDB
      IF ( ICMPQ.EQ.0 ) THEN
         LDQ = 1
         N1  = 0
      ELSE
         LDQ = MAX( 1, N2 )
         N1  = N2
      END IF
      IF ( ICMPU.EQ.0 ) THEN
         LDU = 1
         NU  = 0
      ELSE
         LDU = MAX( 1, N )
         NU  = N2
      END IF
      IF ( MAX( ICMPQ, ICMPU ).EQ.0 ) THEN
         LIWORK = N + 18
      ELSE
         LIWORK = MAX( N + 18, M + 48, 5*M + 1 )
      END IF
C
C Allocate variable dimension local arrays.
C !Fortran 90/95
C
      ALLOCATE(  BWORK( M ), IWORK( LIWORK ) )
      ALPHAI = mxCalloc(      M,     8 )
      ALPHAR = mxCalloc(      M,     8 )
      B      = mxCalloc(  LDB*M,     8 )
      BETA   = mxCalloc(      M,     8 )
      FG     = mxCalloc( LDFG*(M+1), 8 )
      Q      = mxCalloc( LDQ*N1,     8 )
      U      = mxCalloc( LDU*NU,     8 )
      Z      = mxCalloc(  LDZ*N,     8 )
C
C   ldwork
C
      CALL DGHFDF( COMPQ, COMPU, ORTH, N, %VAL( Z ), LDZ, %VAL( B ),
     $             LDB, %VAL( FG ), LDFG, NEIG, %VAL( Q ), LDQ,
     $             %VAL( U ), LDU, %VAL( ALPHAR ), %VAL( ALPHAI ),
     $             %VAL( BETA ), IWORK, LIWORK, DUM, -1, BWORK, IWARN,
     $             INFO )
      LDWORK = INT( DUM( 1 ) )
      DWORK  = mxCalloc( LDWORK, 8 )
C
C Copy inputs from MATLAB workspace to locally allocated arrays.
C
      CALL mxCopyPtrToReal8( mxGetPr( PRHS( 1 ) ), %VAL( Z ),  N*N )
      CALL mxCopyPtrToReal8( mxGetPr( PRHS( 2 ) ), %VAL( B ),  M*M )
      CALL mxCopyPtrToReal8( mxGetPr( PRHS( 3 ) ), %VAL( FG ),
     $                       M*( M + 1 ) )
C
C Do the actual computations.
C
      CALL DGHFDF( COMPQ, COMPU, ORTH, N, %VAL( Z ), LDZ, %VAL( B ),
     $             LDB, %VAL( FG ), LDFG, NEIG, %VAL( Q ), LDQ,
     $             %VAL( U ), LDU, %VAL( ALPHAR ), %VAL( ALPHAI ),
     $             %VAL( BETA ), IWORK, LIWORK, %VAL( DWORK ), LDWORK,
     $             BWORK, IWARN, INFO )
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
         IF ( IWARN.GT.0 ) THEN
            ALLOCATE( DW( M ) )
            CALL DCOPY( M, %VAL( BETA ),  1, DW,  1 )
            CALL DCOPY( 2, %VAL( DWORK ), 1, DUM, 1 )
            DO 10 I = 1, M
               DW( I ) = DW( I )/DUM( 2 )**IWORK( I )
   10       CONTINUE
            CALL DCOPY( M, DW, 1, %VAL( BETA ), 1 )
            DEALLOCATE( DW )
         END IF
         PLHS( 3 ) = mxCreateDoubleMatrix( M, 1, 0 )
         CALL mxCopyReal8ToPtr( %VAL( BETA ), mxGetPr( PLHS( 3 ) ),
     $                          M )
C
         IP = 4
         IF ( NLHS.GE.IP .AND. ICMPQ.GT.0 ) THEN
            CALL DLACPY( 'Full', N, NEIG, %VAL( Q ), LDQ, %VAL( DWORK ),
     $                   LDZ )
            PLHS( IP ) = mxCreateDoubleMatrix( N, NEIG, 0 )
            CALL mxCopyReal8ToPtr( %VAL( DWORK ),
     $                             mxGetPr( PLHS( IP ) ), N*NEIG )
            IP = IP + 1
         END IF
C
         IF ( NLHS.GE.IP .AND. ICMPU.GT.0 ) THEN
            CALL DLACPY( 'Full', N, NEIG, %VAL( U ), LDU, %VAL( DWORK ),
     $                   LDZ )
            PLHS( IP ) = mxCreateDoubleMatrix( N, NEIG, 0 )
            CALL mxCopyReal8ToPtr( %VAL( DWORK ),
     $                             mxGetPr( PLHS( IP ) ), N*NEIG )
            IP = IP + 1
         END IF
C
         IF ( NLHS.GE.IP .AND. MAX( ICMPQ, ICMPU ).GT.0 ) THEN
            PLHS( IP ) = mxCreateDoubleMatrix( 1, 1, 0 )
            TEMP = NEIG
            CALL mxCopyReal8ToPtr( TEMP, mxGetPr( PLHS( IP ) ), 1 )
         END IF
      END IF
C
C Deallocate variable dimension arrays.
C !Fortran 90/95
C
      DEALLOCATE( BWORK, IWORK   )
      CALL mxFree( ALPHAI )
      CALL mxFree( ALPHAR )
      CALL mxFree( B      )
      CALL mxFree( BETA   )
      CALL mxFree( DWORK  )
      CALL mxFree( FG     )
      CALL mxFree( Q      )
      CALL mxFree( U      )
      CALL mxFree( Z      )
C
C Error and warning handling.
C
      IF ( IWARN.NE.0 ) THEN
         WRITE( TEXT, '('' Warning: Eigenvalue denominators BETA'',
     $                  '' might not be representable'')' )
         CALL mexPrintf( TEXT )
      END IF
C
      IF ( INFO.NE.0 ) THEN
         WRITE( TEXT, '('' INFO = '',I4,'' ON EXIT FROM DGHFDF'')'
     $        ) INFO
         CALL mexErrMsgTxt( TEXT )
      END IF
C
      RETURN
C *** Last line of skewHamildeflf ***
      END
