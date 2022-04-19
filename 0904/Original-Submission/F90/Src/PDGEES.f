      SUBROUTINE PDGEES( JOBVS, SORT, SELECT, N, A, IA, JA, DESCA, SDIM, 
     $                   WR, WI, VS, IVS, JVS, DESCVS, DWORK, LDWORK, 
     $                   IWORK, LIWORK, INFO )
C   
C  -- ScaLAPACK/PSLICOT-style routine --
C     Preliminary version.
C     Dept. Computing Science and HPC2N, Univ. of Umeå, Sweden
C     December 10, 2007.
C
      IMPLICIT NONE
C
C     .. Scalar Arguments ..
      CHARACTER          JOBVS, SORT
      INTEGER            INFO, IA, JA, LDVS, LDWORK, LIWORK, N, SDIM,
     $                   IVS, JVS
C     ..
C     .. Array Arguments ..
      INTEGER            DESCA( * ), DESCVS(*), IWORK( * )
      DOUBLE PRECISION   A( * ), VS( * ), WI( * ), DWORK( * ), WR( * )
C     ..
C     .. Function Arguments ..
      LOGICAL            SELECT
      EXTERNAL           SELECT
C     ..
C
C  Purpose
C  =======
C
C  PDGEES computes for an N-by-N real nonsymmetric matrix A, the
C  eigenvalues, the real Schur form T, and, optionally, the matrix of
C  Schur vectors Z.  This gives the Schur factorization A = Z*T*(Z**T).
C
C  Optionally, it also orders the eigenvalues on the diagonal of the
C  real Schur form so that selected eigenvalues are at the top left.
C  The leading columns of Z then form an orthonormal basis for the
C  invariant subspace corresponding to the selected eigenvalues.
C
C  A matrix is in real Schur form if it is upper quasi-triangular with
C  1-by-1 and 2-by-2 blocks. 2-by-2 blocks will be standardized in the
C  form
C          [  a  b  ]
C          [  c  a  ]
C
C  where b*c < 0. The eigenvalues of such a block are a +- sqrt(bc).
C
C  Arguments
C  =========
C
C  JOBVS   (global input) CHARACTER*1
C          = 'N': Schur vectors are not computed;
C          = 'V': Schur vectors are computed.
C
C  SORT    (global input) CHARACTER*1
C          Specifies whether or not to order the eigenvalues on the
C          diagonal of the Schur form.
C          = 'N': Eigenvalues are not ordered;
C          = 'S': Eigenvalues are ordered (see SELECT).
C
C  SELECT  (global input) LOGICAL FUNCTION of two DOUBLE PRECISION 
C          arguments.
C          SELECT must be declared EXTERNAL in the calling subroutine.
C          If SORT = 'S', SELECT is used to select eigenvalues to sort
C          to the top left of the Schur form.
C          If SORT = 'N', SELECT is not referenced.
C          An eigenvalue WR(j)+sqrt(-1)*WI(j) is selected if
C          SELECT(WR(j),WI(j)) is true; i.e., if either one of a complex
C          conjugate pair of eigenvalues is selected, then both complex
C          eigenvalues are selected.
C          Note that a selected complex eigenvalue may no longer
C          satisfy SELECT(WR(j),WI(j)) = .TRUE. after ordering, since
C          ordering may change the value of complex eigenvalues
C          (especially if the eigenvalue is ill-conditioned); in this
C          case INFO is set to N+2 (see INFO below).
C
C  N       (global input) INTEGER
C          The order of the matrix A. N >= 0.
C
C  A       (local input/output) DOUBLE PRECISION pointer into the
C          local memory to an array, dimension (LLD_A, LOCq(N)).
C          This array must contain the local pieces of the N-by-N
C          distributed matrix A. On exit, A has been overwritten by 
C          its real Schur form T.
C
C  IA      (global input) INTEGER
C  JA      (global input) INTEGER
C          The row and column index in the global array A indicating 
C          the first column of sub( A ). IA = JA = 1 must hold for
C          this version.
C
C  DESCA   (global and local input) INTEGER array, dimension DLEN_
C          The array descriptor for the distributed matrix A.
C
C  SDIM    (global output) INTEGER
C          If SORT = 'N', SDIM = 0.
C          If SORT = 'S', SDIM = number of eigenvalues (after sorting)
C                         for which SELECT is true. (Complex conjugate
C                         pairs for which SELECT is true for either
C                         eigenvalue count as 2.)
*
C  WR      (global output) DOUBLE PRECISION array, dimension (N)
C  WI      (global output) DOUBLE PRECISION array, dimension (N)
C          WR and WI contain the real and imaginary parts,
C          respectively, of the computed eigenvalues in the same order
C          that they appear on the diagonal of the output Schur form T.
C          Complex conjugate pairs of eigenvalues will appear
C          consecutively with the eigenvalue having the positive
C          imaginary part first.
C
C  VS      (local output) DOUBLE PRECISION pointer into the
C          local memory to an array, dimension (LLD_VS, LOCq(N)).
C          If JOBVS = 'V', VS contains the orthogonal matrix Z of Schur
C          vectors.
C          If JOBVS = 'N', VS is not referenced.
C
C  IVS     (global input) INTEGER
C  JVS     (global input) INTEGER
C          The row and column index in the global array VS indicating 
C          the first column of sub( VS ). IVS = JVS = 1 must hold for
C          this version.
C
C  DESCVS  (global and local input) INTEGER array, dimension DLEN_
C          The array descriptor for the distributed matrix VS.
C
C  DWORK   (local workspace/output) DOUBLE PRECISION array, dimension (LWORK)
C          On exit, if INFO = 0, WORK(1) contains the optimal LWORK.
C
C  LDWORK  (local input) INTEGER
C          The dimension of the array WORK.  LWORK >= .
C          For good performance, LWORK must generally be larger.
C
C          If LWORK = -1, then a workspace query is assumed; the routine
C          only calculates the optimal size of the WORK array, returns
C          this value as the first entry of the WORK array, and no error
C          message related to LWORK is issued by XERBLA.
C
C  IWORK   (local workspace) INTEGER array, dimension (LIWORK)
C          Not referenced if SORT = 'N'.
C
C  LIWORK  (local input) INTEGER
C          The dimension of the array IWORK.
C          If SORT = 'N', LIWORK >= N.
C
C  INFO    (global output) INTEGER
C          = 0: successful exit
C          < 0: if INFO = -i, the i-th argument had an illegal value.
C          > 0: if INFO = i, and i is
C             <= N: the QR algorithm failed to compute all the
C                   eigenvalues; elements 1:ILO-1 and i+1:N of WR and WI
C                   contain those eigenvalues which have converged; if
C                   JOBVS = 'V', VS contains the matrix which reduces A
C                   to its partially converged Schur form.
C             = N+1: the eigenvalues could not be reordered because some
C                   eigenvalues were too close to separate (the problem
C                   is very ill-conditioned);
C             = N+2: after reordering, roundoff changed values of some
C                   complex eigenvalues so that leading eigenvalues in
C                   the Schur form no longer satisfy SELECT=.TRUE.  This
C                   could also be caused by underflow due to scaling.
C
C  =====================================================================
C
C     .. Parameters ..
      LOGICAL           DEBUG, PRINT
      INTEGER           BLOCK_CYCLIC_2D, CSRC_, CTXT_, DLEN_, DT_,
     $                  LLD_, MB_, M_, NB_, N_, RSRC_
      PARAMETER         ( BLOCK_CYCLIC_2D = 1, DLEN_ = 9, DT_ = 1,
     $                    CTXT_ = 2, M_ = 3, N_ = 4, MB_ = 5, NB_ = 6,
     $                    RSRC_ = 7, CSRC_ = 8, LLD_ = 9,
     $                    DEBUG = .FALSE., PRINT = .FALSE. )
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D0, ONE = 1.0D0 )
C     ..
C     .. Local Scalars ..
      LOGICAL            CURSL, LASTSL, LQUERY, LST2SL, SCALEA, WANTST,
     $                   WANTVS
      INTEGER            HSWORK, I, I1, I2, IBAL, ICOND, IERR, IEVAL,
     $                   IHI, ILO, INXT, IP, ITAU, IWRK, K, MAXB,
     $                   MAXWRK, MINWRK, NPROW, NPCOL, MYROW, MYCOL,
     $                   WRKOPT, IWRKOPT, ICTXT
      DOUBLE PRECISION   ANRM, BIGNUM, CSCALE, EPS, S, SEP, SMLNUM,
     $                   TMP1, TMP2, TMP3, TMP4, COS, SIN, DUM1, DUM2,
     $                   DUM3, DUM4, TIME1, TIME2, TIME3, STAMP
C     ..
C     .. Local Arrays ..
      INTEGER            IDUM( 1 ), PARA(6)
      DOUBLE PRECISION   DUM( 1 )
C     ..
C     .. External Subroutines ..
      EXTERNAL           PDGEHRD, PDLAHQR, PDLACPY, PDLASCL, PDORGHR, 
     $                   PDSWAP, PBDTRSEN, PXERBLA, DLABAD
C     ..
C     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            NUMROC, ILG2NT
      DOUBLE PRECISION   DLAMCH, PDLANGE, MPI_WTIME
      EXTERNAL           LSAME, NUMROC, DLAMCH, PDLANGE, MPI_WTIME,
     $                   ILG2NT
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN, SQRT
C     ..
C     .. Executable Statements ..
C
C     Test the input arguments
C
      IF( DEBUG ) WRITE(*,*) 'PDGEES'
      INFO = 0
      LQUERY = ( LDWORK.EQ.-1.OR.LDWORK.EQ.-1 )
      ICTXT = DESCA( CTXT_ )
      CALL BLACS_GRIDINFO( ICTXT, NPROW, NPCOL, MYROW, MYCOL )
      IF( NPROW.EQ.-1 ) INFO = -(600+CTXT_)
      IF( INFO.EQ.0 ) THEN
         WANTVS = LSAME( JOBVS, 'V' )
         WANTST = LSAME( SORT, 'S' )
         IF( ( .NOT.WANTVS ) .AND. ( .NOT.LSAME( JOBVS, 'N' ) ) ) THEN
            INFO = -1
         ELSE IF( ( .NOT.WANTST ) .AND. ( .NOT.LSAME( SORT, 'N')) ) THEN
            INFO = -2
         ELSE IF( N.LT.0 ) THEN
            INFO = -4
         ELSE 
            CALL CHK1MAT( N, 4, N, 4, IA, JA, DESCA, 8, INFO )
            CALL CHK1MAT( N, 4, N, 4, IVS, JVS, DESCVS, 15, INFO )
            CALL PCHK2MAT( N, 4, N, 4, IA, JA, DESCA, 8, N, 4, N, 4,
     $           IVS, JVS, DESCVS, 15, 0, IWORK, IWORK, INFO )
         END IF
C     
C     Test workspace - TODO!
C     
         IF( INFO.EQ.0 ) THEN
            WRKOPT = N
            IWRKOPT = N
            IF( .NOT.LQUERY.AND.LDWORK.LT.WRKOPT ) THEN
               INFO = -17
            ELSEIF( .NOT.LQUERY.AND.LIWORK.LT.IWRKOPT ) THEN
               INFO = -19
            END IF
         END IF
C     
C     Error or workspace check return.
C     
         IF( .NOT.LQUERY .AND. INFO.NE.0 ) THEN
            CALL PXERBLA( 'PDGEES', -INFO )
            RETURN
         ELSEIF( LQUERY ) THEN
            DWORK( 1 ) = WRKOPT
            IWORK( 1 ) = IWRKOPT
            RETURN
         END IF
      END IF
C     
C     Quick return if possible
C
      IF( N.EQ.0 ) THEN
         SDIM = 0
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
C     Scale A if max element outside range [SMLNUM,BIGNUM]
C
      ANRM = PDLANGE( 'M', N, N, A, IA, JA, DESCA, DWORK )
      SCALEA = .FALSE.
      IF( ANRM.GT.ZERO .AND. ANRM.LT.SMLNUM ) THEN
         SCALEA = .TRUE.
         CSCALE = SMLNUM
      ELSE IF( ANRM.GT.BIGNUM ) THEN
         SCALEA = .TRUE.
         CSCALE = BIGNUM
      END IF
      IF( SCALEA )
     $   CALL PDLASCL( 'G', ANRM, CSCALE, N, N, A, IA, JA, DESCA, IERR )
C
C     Set submatrix indices for reduction to real Schur form
C
C
C     Permute the matrix to make it more nearly triangular
C     (Workspace: need N)
C
      CALL PDGEBAL( 'P', N, A, DESCA, ILO, IHI, CSCALE, INFO )
C      ILO = IA
C      IHI = IA+N-1
C
C     Reduce to upper Hessenberg form
C
      ITAU = 1
      IWRK = ITAU + NUMROC( JA+N-2, DESCA(NB_), MYCOL, DESCA(CSRC_), 
     $                      NPCOL )
      STAMP = MPI_WTIME()
      CALL PDGEHRD( N, ILO, IHI, A, IA, JA, DESCA, DWORK( ITAU ), 
     $              DWORK( IWRK ), LDWORK-IWRK+1, IERR )
      IF( DEBUG ) WRITE(*,*) 'PDGEHRD IERR =',IERR
C
      IF( WANTVS ) THEN
C
C        Generate Hessenberg transformation explicitly in VS and extract
C        explicit Hessenberg matrix from A
C
         CALL PDLASET( 'All', N, N, ZERO, ONE, VS, IVS, JVS, DESCVS )     
         CALL PDORMHR( 'L','N', N, N, ILO, IHI, A, IA, JA, DESCA, 
     $                 DWORK( ITAU ), VS, IVS, JVS, DESCVS, 
     $                 DWORK( IWRK ), LDWORK-IWRK+1, IERR )
         IF( DEBUG ) WRITE(*,*) 'PDORMHR IERR =',IERR
C
         CALL PDLASET( 'Lower triangular', N-2, N-2, ZERO, ZERO, A, 
     $                 IA + 2, JA, DESCA )
      END IF
      TIME1 = MPI_WTIME() - STAMP 
      IF( PRINT ) THEN
         CALL PDLAPRNT( N, N, A, IA, JA, DESCA, 0, 0, 'H1', 6, DWORK )
         CALL PDLAPRNT( N, N, VS, IVS, JVS, DESCVS, 0, 0, 'Q1', 6,DWORK)
      END IF
C
C     Perform QR iterations, accumulating Schur vectors in VS if desired
C
      IWRK = ITAU
      STAMP = MPI_WTIME()
      CALL PDLAHQR( .TRUE., WANTVS, N, ILO, IHI, A, DESCA, WR, WI, ILO, 
     $              IHI, VS, DESCVS, DWORK( IWRK ), LDWORK-IWRK+1,
     $              IWORK, LIWORK, IEVAL )
      IF( DEBUG ) WRITE(*,*) 'PDLAHQR IEVAL =',IEVAL
      IF( IEVAL.GT.0 )
     $   INFO = IEVAL
      IF( PRINT ) THEN
         CALL PDLAPRNT( N, N, A, IA, JA, DESCA, 0, 0, 'H2', 6, DWORK )
         CALL PDLAPRNT( N, N, VS, IVS, JVS, DESCVS, 0, 0, 'Q2', 6,DWORK)
      END IF
C
C     Check for unreduced 2-by-2 blocks corresponding to real eigenvalues
C     and get rid off them to make the upcoming separation possible
C
      DO 7  I = ILO, IHI-1
         CALL PDELGET( 'All', ' ', TMP3, A, I+1, I, DESCA )
         IF( TMP3.NE.0.0D+00 ) THEN
            CALL PDELGET( 'All', ' ', TMP1, A, I, I, DESCA )
            CALL PDELGET( 'All', ' ', TMP2, A, I, I+1, DESCA )
            CALL PDELGET( 'All', ' ', TMP4, A, I+1, I+1,
     $           DESCA )
            CALL DLANV2( TMP1, TMP2, TMP3, TMP4, DUM1, DUM2, DUM3,
     $           DUM4, COS, SIN )
            IF( TMP3.EQ.0.0D+00 ) THEN
               IF( I+2.LE.N )
     $              CALL PDROT( N-I-1, A, I, I+2, DESCA,
     $              N, A, I+1, I+2, DESCA, N, COS, SIN, 
     $              DWORK(IWRK), LDWORK-IWRK+1, IERR )
               CALL PDROT( I-1, A, 1, I, DESCA, 1,
     $              A, 1, I+1, DESCA, 1, COS, SIN,
     $              DWORK(IWRK), LDWORK-IWRK+1, IERR )
               CALL PDROT( N, VS, 1, I, DESCVS, 1,
     $              VS, 1, I+1, DESCVS, 1, COS, SIN,
     $              DWORK(IWRK), LDWORK-IWRK+1, IERR )
               CALL PDELSET( A, I, I, DESCA, TMP1 )
               CALL PDELSET( A, I, I+1, DESCA, TMP2 )
               CALL PDELSET( A, I+1, I, DESCA, TMP3 )
               CALL PDELSET( A, I+1, I+1, DESCA, TMP4 )
            END IF
         END IF
 7    CONTINUE
      IF( PRINT ) THEN
         CALL PDLAPRNT( N, N, A, IA, JA, DESCA, 0, 0, 'T', 6, DWORK )
         CALL PDLAPRNT( N, N, VS, IVS, JVS, DESCVS, 0, 0, 'Z', 6, DWORK)
      END IF
      IF( DEBUG ) WRITE(*,*) 'PDROT IERR =',IERR
C
C     Set dimension of invariant subspace
C
      SDIM = 0
C
C     Set optimal configuration parameters for parallel reordering
C
      PARA(1) = MIN(NPROW,NPCOL)
      PARA(2) = MIN(DESCA(MB_)/2,40)
      PARA(3) = MIN(DESCA(MB_),80)
      PARA(4) = 50
      PARA(5) = MIN(DESCA(NB_),32)
      PARA(6) = MIN(DESCA(MB_)/2,40)
C
C     Read out the eigenvalues
C
      DO 8 I = 1, N
         IWORK( I ) = 0
 8    CONTINUE
      CALL PBDTRSEN( 'N', 'N', IWORK, PARA, N, A, IA, JA, DESCA, VS, 
     $               IVS, JVS, DESCVS, WR, WI, SDIM, S, SEP, 
     $               DWORK(IWRK), LDWORK-IWRK+1, IWORK(N+1), 
     $               LIWORK-N, ICOND )
      IF( DEBUG ) WRITE(*,*) 'PBDTRSEN 1 ICOND =',ICOND
      TIME2 = MPI_WTIME() - STAMP 
C     
C     Sort eigenvalues if desired
C
      STAMP = MPI_WTIME()
      IF( WANTST .AND. INFO.EQ.0 ) THEN
         IF( SCALEA ) THEN
            CALL DLASCL( 'G', 0, 0, CSCALE, ANRM, N, 1, WR, N, IERR )
            CALL DLASCL( 'G', 0, 0, CSCALE, ANRM, N, 1, WI, N, IERR )
         END IF
         DO 10 I = 1, N
            IWORK( I ) = ILG2NT(SELECT( WR( I ), WI( I ) ))
   10    CONTINUE
C
C        Reorder eigenvalues and transform Schur vectors
C
         CALL PBDTRSEN( 'N', 'V', IWORK, PARA, N, A, IA, JA, DESCA, VS, 
     $                  IVS, JVS, DESCVS, WR, WI, SDIM, S, SEP, 
     $                  DWORK(IWRK), LDWORK-IWRK+1, IWORK(N+1), 
     $                  LIWORK-N, ICOND )
         IF( DEBUG ) WRITE(*,*) 'PBDTRSEN 2 ICOND =',ICOND
         IF( ICOND.GT.0 )
     $      INFO = N + ICOND
      END IF
      TIME3 = MPI_WTIME() - STAMP
      IF( PRINT ) THEN
         CALL PDLAPRNT( N, N, A, IA, JA, DESCA, 0, 0, 'TT', 6, DWORK )
         CALL PDLAPRNT( N, N, VS, IVS, JVS, DESCVS, 0, 0, 'ZZ', 6,DWORK)
      END IF
C
      IF( SCALEA ) THEN
C
C        Undo scaling for the Schur form of A and read out eigenvalues
C        again
C
         CALL PDLASCL( 'H', CSCALE, ANRM, N, N, A, IA, JA, DESCA, IERR )
         DO 12 I = 1, N
            IWORK( I ) = 0
 12      CONTINUE
         CALL PBDTRSEN( 'N', 'N', IWORK, PARA, N, A, IA, JA, DESCA, VS, 
     $                  IVS, JVS, DESCVS, WR, WI, SDIM, S, SEP, 
     $                  DWORK(IWRK), LDWORK-IWRK+1, IWORK(N+1), 
     $                  LIWORK-N, ICOND )  
         IF( DEBUG ) WRITE(*,*) 'PBDTRSEN 3 ICOND =',ICOND
         IF( CSCALE.EQ.SMLNUM ) THEN
C
C           If scaling back towards underflow, adjust WI if an
C           offdiagonal element of a 2-by-2 block in the Schur form
C           underflows.
C
            IF( IEVAL.GT.0 ) THEN
               I1 = IEVAL + 1
               I2 = IHI - 1
               CALL DLASCL( 'G', 0, 0, CSCALE, ANRM, ILO-1, 1, WI,
     $                      MAX( ILO-1, 1 ), IERR )
            ELSE IF( WANTST ) THEN
               I1 = 1
               I2 = N - 1
            ELSE
               I1 = ILO
               I2 = IHI - 1
            END IF
            INXT = I1 - 1
            DO 20 I = I1, I2
               IF( I.LT.INXT )
     $            GO TO 20
               IF( WI( I ).EQ.ZERO ) THEN
                  INXT = I + 1
               ELSE
                  CALL PDELGET( 'All', '1-Tree', TMP1, A, I+1, I, DESCA)
                  CALL PDELGET( 'All', '1-Tree', TMP2, A, I, I+1, DESCA)
                  IF( TMP1.EQ.ZERO ) THEN
                     WI( I ) = ZERO
                     WI( I+1 ) = ZERO
                  ELSE IF( TMP1.NE.ZERO .AND. TMP2.EQ.ZERO ) THEN
                     WI( I ) = ZERO
                     WI( I+1 ) = ZERO
                     IF( DEBUG ) WRITE(*,*) 'PDSWAP 1'
                     IF( I.GT.1 )
     $                  CALL PDSWAP( I-1, A, 1, I, DESCA, 1, A, 1, I+1, 
     $                               DESCA, 1 )
                     IF( DEBUG ) WRITE(*,*) 'PDSWAP 2'
                     IF( N.GT.I+1 )
     $                    CALL PDSWAP( N-I-1, A, I, I+2, DESCA, 
     $                    DESCA(M_), A, I+1, I+2, DESCA, DESCA(M_) )
                     IF( DEBUG ) WRITE(*,*) 'PDSWAP 3'
                     CALL PDSWAP( N, VS, 1, I, DESCVS, 1, VS, 1, I+1, 
     $                            DESCA, 1 )
                     CALL PDELSET( A, I, I+1, TMP1, DESCA ) 
                     CALL PDELSET( A, I+1, I, ZERO, DESCA ) 
                  END IF
                  INXT = I + 2
               END IF
   20       CONTINUE
         END IF
C
C        Undo scaling for the imaginary part of the eigenvalues
C
         CALL DLASCL( 'G', 0, 0, CSCALE, ANRM, N-IEVAL, 1,
     $                WI( IEVAL+1 ), MAX( N-IEVAL, 1 ), IERR )
      END IF
C
      IF( WANTST .AND. INFO.EQ.0 ) THEN
C
C        Check if reordering successful
C
         LASTSL = .TRUE.
         LST2SL = .TRUE.
         SDIM = 0
         IP = 0
         DO 30 I = 1, N
            CURSL = SELECT( WR( I ), WI( I ) )
            IF( WI( I ).EQ.ZERO ) THEN
               IF( CURSL )
     $            SDIM = SDIM + 1
               IP = 0
               IF( CURSL .AND. .NOT.LASTSL )
     $            INFO = N + 2
            ELSE
               IF( IP.EQ.1 ) THEN
C
C                 Last eigenvalue of conjugate pair
C
                  CURSL = CURSL .OR. LASTSL
                  LASTSL = CURSL
                  IF( CURSL )
     $               SDIM = SDIM + 2
                  IP = -1
                  IF( CURSL .AND. .NOT.LST2SL )
     $               INFO = N + 2
               ELSE
C
C                 First eigenvalue of conjugate pair
C
                  IP = 1
               END IF
            END IF
            LST2SL = LASTSL
            LASTSL = CURSL
   30    CONTINUE
      END IF
C
      DWORK( 1 ) = WRKOPT
      DWORK( 2 ) = TIME1
      DWORK( 3 ) = TIME2
      DWORK( 4 ) = TIME3
      IWORK( 1 ) = IWRKOPT
C
      IF( DEBUG ) WRITE(*,*) 'End of PDGEES'
      RETURN
C *** Last line of PDGEES ***
      END
