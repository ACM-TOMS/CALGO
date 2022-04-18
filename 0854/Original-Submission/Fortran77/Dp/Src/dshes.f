      SUBROUTINE DSHES( JOBU, N, A, LDA, QG, LDQG, U1, LDU1, U2, LDU2,
     $                  WR, WI, DWORK, LDWORK, INFO )
      IMPLICIT NONE
C
C     PURPOSE
C
C     To compute the real skew-Hamiltonian Schur form of a
C     skew-Hamiltonian matrix,
C
C                   [  A   G  ]
C             W  =  [       T ],
C                   [  Q   A  ]
C
C     where A is an N-by-N matrix and G,Q are N-by-N skew-symmetric
C     matrices. That is, an orthogonal symplectic U is computed so that
C
C               T       [  Aout  Gout  ]
C              U W U =  [            T ] ,
C                       [    0   Aout  ]
C
C     where Aout is in Schur canonical form (as returned by the LAPACK
C     routine DHSEQR). That is, A is block upper triangular with 1-by-1
C     and 2-by-2 diagonal blocks; each 2-by-2 diagonal block has its
C     diagonal elements equal and its off-diagonal elements of opposite
C     sign.
C
C     Optionally, the matrix U is returned in terms of its first N/2
C     rows
C
C                      [  U1   U2 ]
C                  U = [          ].
C                      [ -U2   U1 ]
C
C     ARGUMENTS
C
C     Mode Parameters
C
C     JOBU    (input) CHARACTER*1
C             = 'N': transformation matrix U is not computed;
C             = 'U': transformation matrix U is computed.
C
C     Input/Output Parameters
C
C     N       (input) INTEGER
C             The order of the matrix A. N >= 0.
C
C     A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
C             On entry, the leading N-by-N part of this array must
C             contain the matrix A.
C             On exit, the leading N-by-N part of this array contains
C             the matrix Aout in Schur canonical form.
C
C     LDA     INTEGER
C             The leading dimension of the array A.  LDA >= MAX(1,N).
C
C     QG      (input/output) DOUBLE PRECISION array, dimension
C                            (LDQG,N+1)
C             On entry, the leading N-by-N+1 part of this array must
C             contain in columns 1:N the strictly lower triangular part
C             of the matrix Q and in columns 2:N+1 the strictly upper
C             triangular part of the matrix G.
C             On exit, the leading N-by-N+1 part of this array must
C             contain in columns 2:N+1 the strictly upper triangular
C             part of the skew-symmetric matrix G. The part which
C             contained the matrix Q is set to zero.
C             Note that the parts containing the diagonal and the first
C             supdiagonal of this array are not overwritten by zeros
C             only if JOBU = 'U' or LDWORK >= 2*N*N - N.
C
C     LDQG    INTEGER
C             The leading dimension of the array QG.  LDQG >= MAX(1,N).
C
C     U1      (output) DOUBLE PRECISION array, dimension (LDU1,N)
C             On exit, the leading N-by-N part of this array contains
C             the matrix U1.
C
C     LDU1    INTEGER
C             The leading dimension of the array U1.
C             LDU1 >= MAX(1,N),  if JOBU = 'U'; 
C             LDU1 >= 1,         if JOBU = 'N'.
C
C
C     U2      (output) DOUBLE PRECISION array, dimension (LDU2,N)
C             On exit, the leading N-by-N part of this array contains
C             the matrix U2.
C
C     LDU2    INTEGER
C             The leading dimension of the array U1.
C             LDU2 >= MAX(1,N),  if JOBU = 'U'; 
C             LDU2 >= 1,         if JOBU = 'N'.
C
C     WR      (output) DOUBLE PRECISION array, dimension (N)
C     WI      (output) DOUBLE PRECISION array, dimension (N)
C             The real and imaginary parts, respectively, of the
C             eigenvalues of Aout, which are half of the eigenvalues
C             of S. The eigenvalues are stored in the same order as on
C             the diagonal of Aout, with WR(i) = Aout(i,i) and, if
C             Aout(i:i+1,i:i+1) is a 2-by-2 diagonal block, WI(i) > 0
C             and WI(i+1) = -WI(i).
C
C     Workspace
C
C     DWORK   DOUBLE PRECISION array, dimension (LDWORK)
C             On exit, if INFO = 0,  DWORK(1)  returns the optimal
C             value of LDWORK.
C             On exit, if  INFO = -14,  DWORK(1)  returns the minimum
C             value of LDWORK.
C
C     LDWORK  INTEGER
C             The length of the array DWORK.
C             LDWORK >= MAX(1,(N+5)*N),      if JOBU = 'U';
C             LDWORK >= MAX(1,5*N,(N+1)*N),  if JOBU = 'N'.
C
C     Error Indicator
C
C     INFO    INTEGER
C             = 0:  successful exit;
C             < 0:  if INFO = -i, the i-th argument had an illegal 
C                   value.
C
C     METHOD
C
C     First, using routine DSHPVB, an orthogonal symplectic matrix UP
C     is computed so that
C
C                T      [  AP  GP  ]
C            UP W UP =  [        T ]
C                       [  0   AP  ]
C
C     is in Paige/Van Loan form. Next, the LAPACK routine DHSEQR is
C     applied to the matrix AP to compute an orthogonal matrix V
C     so that Aout = V'*AP*V is in Schur canonical form.
C     Finally, the transformations
C
C                       [ V  0 ]
C             U =  UP * [      ],     Gout = V'*G*V,
C                       [ 0  V ]
C
C     using the routines DSKRKB/DSKUPD for the latter, are performed.
C
C     REFERENCES
C
C     [1] C. F. VAN LOAN:
C     A symplectic method for approximating all the eigenvalues of
C     a Hamiltonian matrix.
C     Linear Algebra and its Applications 61 (1984), pp. 233-251.
C
C     [2] D. KRESSNER:
C     Block algorithms for orthogonal symplectic factorizations.
C     To appear in BIT, 2002.
C
C     CONTRIBUTORS
C
C     D. Kressner (Technical Univ. Berlin, Germany) and
C     P. Benner (Technical Univ. Chemnitz, Germany), December 2003.
C
C     KEYWORDS
C
C     Schur form, eigenvalues, skew-Hamiltonian matrix.
C
C     ******************************************************************
C
C     .. Parameters ..
      DOUBLE PRECISION  ZERO, ONE
      PARAMETER         ( ZERO = 0.0D0, ONE = 1.0D0 )
C     .. Scalar Arguments ..
      INTEGER           ILO, INFO, LDA, LDQG, LDU1, LDU2, LDWORK, N
      CHARACTER*1       JOBU
C     .. Array Arguments ..
      DOUBLE PRECISION  A(LDA,*), DWORK(*), QG(LDQG,*), U1(LDU1,*),
     $                  U2(LDU2,*), WR(*), WI(*)
C     .. Local Scalars ..
      LOGICAL           COMPU, SCALEW
      INTEGER           I, I1, I2, IERR, INXT, PBAL, PCS, PDV, PDW, PHO,
     $                  PTAU, WRKMIN, WRKOPT
      DOUBLE PRECISION  BIGNUM, CSCALE, EPS, SMLNUM, WNRM
C     .. External Function ..
      EXTERNAL          DLAMCH, DLANHA, LSAME
      LOGICAL           LSAME
      DOUBLE PRECISION  DLANHA, DLAMCH
C     .. External Subroutines ..
      EXTERNAL          DCOPY, DHABAK, DHSEQR, DLABAD, DLACPY, DLASCL,
     $                  DLASET, DOSMPV, DSCAL, DSHBAL, DSHPVB, DSKRKB,
     $                  DSKUPD, DSWAP, XERBLA
C     .. Intrinsic Functions ..
      INTRINSIC         DBLE, MAX
C
C     .. Executable Statements ..
C
C     Check the scalar input parameters.
C
      INFO  = 0
      COMPU = LSAME( JOBU,  'U' )
C
      IF ( COMPU ) THEN
         WRKMIN = MAX( 1, (N+5)*N )
         WRKOPT = WRKMIN
      ELSE
         WRKMIN = MAX( 1, 5*N,(N+1)*N )
         WRKOPT = MAX( WRKMIN, 2*N*(N-1) )
      END IF
C
      IF ( .NOT.COMPU .AND. .NOT.LSAME( JOBU, 'N' ) ) THEN
         INFO = -1
      ELSE IF ( N.LT.0 ) THEN
         INFO = -2
      ELSE IF ( LDA.LT.MAX( 1,N ) ) THEN
         INFO = -4
      ELSE IF ( LDQG.LT.MAX( 1,N ) ) THEN
         INFO = -6
      ELSE IF ( LDU1.LT.1 .OR. ( COMPU.AND.LDU1.LT.N ) ) THEN
         INFO = -8
      ELSE IF ( LDU2.LT.1 .OR. ( COMPU.AND.LDU2.LT.N ) ) THEN
         INFO = -10
      ELSE IF ( LDWORK.LT.WRKMIN ) THEN
         DWORK(1) = DBLE( WRKMIN )
         INFO = -14
      END IF
C
C     Return if there were illegal values.
C
      IF ( INFO.NE.0 ) THEN
         CALL XERBLA( 'DSHES', -INFO )
         RETURN
      END IF
C
C     Quick return if possible.
C
      IF ( N.EQ.0 ) THEN
         DWORK(1) = ONE
         RETURN
      END IF
C
C     Get machine constants
C
      EPS = DLAMCH( 'P' )
      SMLNUM = DLAMCH( 'S' )
      BIGNUM = ONE / SMLNUM
      CALL DLABAD( SMLNUM, BIGNUM )
      SMLNUM = SQRT( SMLNUM ) / EPS
      BIGNUM = ONE / SMLNUM
C
C     Scale W if max element outside range [SMLNUM,BIGNUM]
C
      WNRM = DLANHA( 'Skew-Hamiltonian', 'Max-Norm', N, A, LDA, QG,
     $               LDQG, DWORK )
      SCALEW = .FALSE.
      IF ( WNRM.GT.ZERO .AND. WNRM.LT.SMLNUM ) THEN
         SCALEW = .TRUE.
         CSCALE = SMLNUM
      ELSE IF ( WNRM.GT.BIGNUM ) THEN
         SCALEW = .TRUE.
         CSCALE = BIGNUM
      END IF
      IF( SCALEW ) THEN
         CALL DLASCL( 'General', 0, 0, WNRM, CSCALE, N, N, A, LDA,
     $                IERR )
         IF ( N.GT.1 ) THEN
            CALL DLASCL( 'Lower', 0, 0, WNRM, CSCALE, N-1, N-1, QG(2,1),
     $                   LDQG, IERR )
            CALL DLASCL( 'Upper', 0, 0, WNRM, CSCALE, N-1, N-1, QG(1,3),
     $                   LDQG, IERR )
         END IF
      END IF
C
C     Permute to make W closer to skew-Hamiltonian Schur form.
C
      PBAL = 1
      CALL DSHBAL( 'Permute', N, A, LDA, QG, LDQG, ILO, DWORK(PBAL),
     $             IERR )
C
C     Reduce to Paige/Van Loan form
C
      PCS  = N + PBAL
      PTAU = 2*N + PCS
      PDW = N + PTAU
      CALL DSHPVB( N, ILO, A, LDA, QG, LDQG, DWORK(PCS), DWORK(PTAU),
     $             DWORK(PDW), LDWORK-PDW+1, IERR )
      WRKOPT = MAX( WRKOPT, INT( DWORK(PDW) ) + PDW - 1 )
      IF ( COMPU ) THEN
C
C        Copy information about Householder vectors to workspace.
C
         PHO = PDW
         PDW = PDW + N*N
         CALL DLACPY( 'L', N, N, A, LDA, DWORK(PHO), N )
C
C        Perform QR iteration, accumulating Schur vectors in U1.
C
         CALL DHSEQR( 'Schur', 'Initialize', N, MIN( ILO, N ), N, A,
     $                LDA, WR, WI, U1, LDU1, DWORK(PDW), LDWORK-PDW+1,
     $                INFO )
C
C        Update G = V'*G*V.
C
         CALL DSKRKB( 'Upper', 'Transpose', N, N, ZERO, ONE, U1, LDU1,
     $                 QG(1,2), LDQG, QG(1,2), LDQG, U2, LDU2, INFO )
C
C        Apply orthogonal symplectic matrix from PVL reduction to [V;0].
C
         CALL DSCAL( N-1, -ONE, DWORK(PCS+1), 2 )
         CALL DLASET( 'All', N, N, ZERO, ZERO, U2, LDU2 )
         CALL DOSMPV( 'No Transpose', 'No Transpose', 'No Transpose', N,
     $                N, ILO, DWORK(PHO), N, QG, LDQG, U1, LDU1, U2,
     $                LDU2, DWORK(PCS), DWORK(PTAU), DWORK(PDW),
     $                LDWORK-PDW+1, IERR )
         WRKOPT = MAX( WRKOPT, INT( DWORK(PDW) ) + PDW - 1 )
C
C        Annihilate Q.
C
         IF ( N.GT.1 )
     $      CALL DLASET( 'Lower', N-1, N-1, ZERO, ZERO, QG(2,1), LDQG )
C
C        Undo balancing
C
         CALL DHABAK( 'Permute', 'Positive', N, ILO, DWORK(PBAL), N,
     $                U1, LDU1, U2, LDU2, IERR )
      ELSE
C
C        Perform QR iteration, accumulating Schur vectors in U1.
C
         PDV = 1
         PDW = N*N+PDV
         CALL DHSEQR( 'Schur', 'Initialize', N, MIN( ILO, N ), N, A,
     $                LDA, WR, WI, DWORK(PDV), N, DWORK(PDW),
     $                LDWORK-PDW+1, INFO )
C
C        Update G = V'*G*V
C
         IF ( LDWORK-PDW+1.GE.N*(N-1) ) THEN
            CALL DSKRKB( 'Upper', 'Transpose', N, N, ZERO, ONE,
     $                   DWORK(PDV), N, QG(1,2), LDQG, QG(1,2), LDQG,
     $                   DWORK(PDW), N-1, INFO )
C
C           Annihilate Q.
C
            IF ( N.GT.1 )
     $         CALL DLASET( 'Lower', N-1, N-1, ZERO, ZERO, QG(2,1),
     $                      LDQG )
         ELSE
            CALL DSKUPD( 'Upper', 'Transpose', N, N, QG(1,2), LDQG,
     $                   DWORK(PDV), N, DWORK(PDW), INFO )
            CALL DLASET( 'All', N, 1, ZERO, ZERO, QG, LDQG )
            IF ( N.GT.1 )
     $         CALL DLASET( 'Lower', N, N, ZERO, ZERO, QG(1,2), LDQG )
         END IF
      END IF
C
      IF ( SCALEW ) THEN
C
C        Undo scaling for the skew-Hamiltonian Schur form.
C
         CALL DLASCL( 'Hessenberg', 0, 0, CSCALE, WNRM, N, N, A, LDA,
     $                IERR )
         IF ( N.GT.1 ) THEN
            CALL DLASCL( 'Upper', 0, 0, CSCALE, WNRM, N-1, N-1, QG(1,3),
     $                   LDQG, IERR )
         END IF
         CALL DCOPY( N, A, LDA+1, WR, 1 )
C
         IF ( CSCALE.EQ.SMLNUM ) THEN
C
C           If scaling back towards underflow, adjust WI if an
C           offdiagonal element of a 2-by-2 block in the Schur form
C           underflows.
C
            IF( INFO.GT.0 ) THEN
               I1 = INFO + 1
               I2 = N-1
               CALL DLASCL( 'General', 0, 0, CSCALE, WNRM, ILO-1, 1, WI,
     $                      MAX( ILO-1, 1 ), IERR )
            ELSE
               I1 = ILO
               I2 = N - 1
            END IF
            INXT = I1 - 1
            DO 10 I = I1, I2
               IF ( I.LT.INXT )
     $            GO TO 10
               IF ( WI(I).EQ.ZERO ) THEN
                  INXT = I + 1
               ELSE
                  IF ( A(I+1,I).EQ.ZERO ) THEN
                     WI(I) = ZERO
                     WI(I+1) = ZERO
                  ELSE IF ( A(I+1,I).NE.ZERO .AND. A(I,I+1).EQ.
     $                     ZERO ) THEN
                     WI(I) = ZERO
                     WI(I+1) = ZERO
                     IF ( I.GT.1 )
     $                  CALL DSWAP( I-1, A(1,I), 1, A(1,I+1), 1 )
                     IF( N.GT.I+1 )
     $                  CALL DSWAP( N-I-1, A(I,I+2), LDA,
     $                              A(I+1,I+2), LDA )
                     A(I,I+1) = A(I+1,I)
                     A(I+1,I ) = ZERO
C
                     CALL DSWAP( I-1, QG(1,I+2), 1, QG(1,I+1), 1 )

                     IF ( N.GT.I+1 )
     $                  CALL DSWAP( N-I-1, QG(I+1,I+3), LDQG, QG(I,I+3),
     $                              LDQG )
                     QG(I,I+2) = -QG(I,I+2)
                     IF ( COMPU ) THEN
                        CALL DSWAP( N, U1(1,I), 1, U1(1,I+1), 1 )
                        CALL DSWAP( N, U2(1,I), 1, U2(1,I+1), 1 )
                     END IF
                  END IF
                  INXT = I + 2
               END IF
   10       CONTINUE
         END IF
C
C        Undo scaling for imaginary parts of the eigenvalues.
C
         CALL DLASCL( 'General', 0, 0, CSCALE, WNRM, N-INFO, 1,
     $                WI(INFO+1), MAX(N-INFO,1), IERR )
      END IF
C
      DWORK(1) = DBLE( WRKOPT )
C
      RETURN
C *** Last line of DSHES ***
      END
