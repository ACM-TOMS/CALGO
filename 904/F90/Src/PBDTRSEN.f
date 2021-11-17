      SUBROUTINE PBDTRSEN( JOB, COMPQ, SELECT, PARA, N, T, IT, JT, 
     $     DESCT, Q, IQ, JQ, DESCQ, WR, WI, M, S, SEP, DWORK, LDWORK, 
     $     IWORK, LIWORK, INFO )
C   
C  -- ScaLAPACK-style routine --
C     Preliminary version.
C     Dept. Computing Science and HPC2N, Univ. of Umeå, Sweden
C     March 7, 2006.
C
      IMPLICIT NONE
C
C     .. Scalar Arguments ..
      CHARACTER          COMPQ, JOB
      INTEGER            INFO, LIWORK, LDWORK, M, N,
     $                   IT, JT, IQ, JQ
      DOUBLE PRECISION   S, SEP
C     ..
C     .. Array Arguments ..
      INTEGER            SELECT( N )
      INTEGER            PARA( 6 ), DESCT( * ), DESCQ( * ), IWORK( * )
      DOUBLE PRECISION   Q( * ), T( * ), WI( * ), DWORK( * ), WR( * )
C     ..
C
C  Purpose
C  =======
C     
C  PBDTRSEN reorders the real Schur factorization of a real matrix
C  A = Q * T * Q**T, so that a selected cluster of eigenvalues appears 
C  in the leading diagonal blocks of the upper quasi-triangular matrix 
C  T, and the leading columns of Q form an orthonormal basis of the
C  corresponding right invariant subspace.
C     
C  Optionally the routine computes the reciprocal condition numbers of
C  the cluster of eigenvalues and/or the invariant subspace.
C     
C  T must be in Schur form (as returned by PDLAHQR), that is, block 
C  upper triangular with 1-by-1 and 2-by-2 diagonal blocks.
C     
C  This subroutine uses a delay and accumulate procedure for performing
C  the off-diagonal updates (see references for details).
C
C  Notes
C  =====
C
C  Each global data object is described by an associated description
C  vector.  This vector stores the information required to establish
C  the mapping between an object element and its corresponding process
C  and memory location.
C
C  Let A be a generic term for any 2D block cyclicly distributed array.
C  Such a global array has an associated description vector DESCA.
C  In the following comments, the character _ should be read as
C  "of the global array".
C
C  NOTATION        STORED IN      EXPLANATION
C  --------------- -------------- --------------------------------------
C  DTYPE_A(global) DESCA( DTYPE_ )The descriptor type.  In this case,
C                                 DTYPE_A = 1.
C  CTXT_A (global) DESCA( CTXT_ ) The BLACS context handle, indicating
C                                 the BLACS process grid A is distribu-
C                                 ted over. The context itself is glo-
C                                 bal, but the handle (the integer
C                                 value) may vary.
C  M_A    (global) DESCA( M_ )    The number of rows in the global
C                                 array A.
C  N_A    (global) DESCA( N_ )    The number of columns in the global
C                                 array A.
C  MB_A   (global) DESCA( MB_ )   The blocking factor used to distribute
C                                 the rows of the array.
C  NB_A   (global) DESCA( NB_ )   The blocking factor used to distribute
C                                 the columns of the array.
C  RSRC_A (global) DESCA( RSRC_ ) The process row over which the first
C                                 row of the array A is distributed.
C  CSRC_A (global) DESCA( CSRC_ ) The process column over which the
C                                 first column of the array A is
C                                 distributed.
C  LLD_A  (local)  DESCA( LLD_ )  The leading dimension of the local
C                                 array.  LLD_A >= MAX(1,LOCr(M_A)).
C
C  Let K be the number of rows or columns of a distributed matrix,
C  and assume that its process grid has dimension p x q.
C  LOCr( K ) denotes the number of elements of K that a process
C  would receive if K were distributed over the p processes of its
C  process column.
C  Similarly, LOCc( K ) denotes the number of elements of K that a
C  process would receive if K were distributed over the q processes of
C  its process row.
C  The values of LOCr() and LOCc() may be determined via a call to the
C  ScaLAPACK tool function, NUMROC:
C          LOCr( M ) = NUMROC( M, MB_A, MYROW, RSRC_A, NPROW ),
C          LOCc( N ) = NUMROC( N, NB_A, MYCOL, CSRC_A, NPCOL ).
C  An upper bound for these quantities may be computed by:
C          LOCr( M ) <= ceil( ceil(M/MB_A)/NPROW )*MB_A
C          LOCc( N ) <= ceil( ceil(N/NB_A)/NPCOL )*NB_A
C     
C  Arguments
C  =========
C     
C  JOB     (global input) CHARACTER*1
C          Specifies whether condition numbers are required for the
C          cluster of eigenvalues (S) or the invariant subspace (SEP):
C          = 'N': none;
C          = 'E': for eigenvalues only (S);
C          = 'V': for invariant subspace only (SEP);
C          = 'B': for both eigenvalues and invariant subspace (S and
C                 SEP).
C     
C  COMPQ   (global input) CHARACTER*1
C          = 'V': update the matrix Q of Schur vectors;
C          = 'N': do not update Q.
C     
C  SELECT  (global input/output) INTEGER/LOGICAL  array, dimension (N)
C          SELECT specifies the eigenvalues in the selected cluster. To
C          select a real eigenvalue w(j), SELECT(j) must be set to
C          .TRUE. (1). To select a complex conjugate pair of eigenvalues
C          w(j) and w(j+1), corresponding to a 2-by-2 diagonal block,
C          either SELECT(j) or SELECT(j+1) or both must be set to
C          .TRUE. (1); a complex conjugate pair of eigenvalues must be
C          either both included in the cluster or both excluded.
C          On output, the (partial) reordering is displayed.
C
C  PARA    (global input) INTEGER*6
C          Block parameters (some should be replaced by calls to
C          PILAENV and others by meaningful default values):
C          PARA(1) = maximum number of concurrent computational windows
C                    allowed in the algorithm; 
C                    0 < PARA(1) <= min(NPROW,NPCOL) must hold;
C          PARA(2) = number of eigenvalues in each window; 
C                    0 < PARA(2) < PARA(3) must hold; 
C          PARA(3) = window size; PARA(2) < PARA(3) < DESCT(MB_) 
C                    must hold;
C          PARA(4) = minimal percentage of flops required for
C                    performing matrix-matrix multiplications instead
C                    of pipelined orthogonal transformations;
C                    0 <= PARA(4) <= 100 must hold;
C          PARA(5) = width of block column slabs for row-wise
C                    application of pipelined orthogonal
C                    transformations in their factorized form; 
C                    0 < PARA(5) <= DESCT(MB_) must hold.
C          PARA(6) = the maximum number of eigenvalues moved together
C                    over a process border; in practice, this will be
C                    approximately half of the cross border window size
C                    0 < PARA(6) <= PARA(2) must hold; 
C
C  N       (global input) INTEGER
C          The order of the globally distributed matrix T. N >= 0.
C     
C  T       (local input/output) DOUBLE PRECISION array, 
C          dimension (LLD_T,LOCc(N)).
C          On entry, the local pieces of the global distributed 
C          upper quasi-triangular matrix T, in Schur form. On exit, T is 
C          overwritten by the local pieces of the reordered matrix T, 
C          again in Schur form, with the selected eigenvalues in the 
C          globally leading diagonal blocks.
C
C  IT      (global input) INTEGER
C  JT      (global input) INTEGER
C          The row and column index in the global array T indicating the
C          first column of sub( T ). IT = JT = 1 must hold.
C     
C  DESCT   (global and local input) INTEGER array of dimension DLEN_.
C          The array descriptor for the global distributed matrix T.
C     
C  Q       (local input/output) DOUBLE PRECISION array, 
C          dimension (LLD_Q,LOCc(N)).
C          On entry, if COMPQ = 'V', the local pieces of the global
C          distributed matrix Q of Schur vectors.
C          On exit, if COMPQ = 'V', Q has been postmultiplied by the
C          global orthogonal transformation matrix which reorders T; the
C          leading M columns of Q form an orthonormal basis for the
C          specified invariant subspace.
C          If COMPQ = 'N', Q is not referenced.
C
C  IQ      (global input) INTEGER
C  JQ      (global input) INTEGER
C          The column index in the global array Q indicating the
C          first column of sub( Q ). IQ = JQ = 1 must hold.
C     
C  DESCQ   (global and local input) INTEGER array of dimension DLEN_.
C          The array descriptor for the global distributed matrix Q.
C     
C  WR      (global output) DOUBLE PRECISION array, dimension (N)
C  WI      (global output) DOUBLE PRECISION array, dimension (N)
C          The real and imaginary parts, respectively, of the reordered
C          eigenvalues of T. The eigenvalues are in principle stored in 
C          the same order as on the diagonal of T, with WR(i) = T(i,i) 
C          and, if T(i:i+1,i:i+1) is a 2-by-2 diagonal block, WI(i) > 0 
C          and WI(i+1) = -WI(i). 
C          Note also that if a complex eigenvalue is sufficiently 
C          ill-conditioned, then its value may differ significantly 
C          from its value before reordering. 
C     
C  M       (global output) INTEGER
C          The dimension of the specified invariant subspace.
C          0 <= M <= N.
C     
C  S       (global output) DOUBLE PRECISION
C          If JOB = 'E' or 'B', S is a lower bound on the reciprocal
C          condition number for the selected cluster of eigenvalues.
C          S cannot underestimate the true reciprocal condition number
C          by more than a factor of sqrt(N). If M = 0 or N, S = 1.
C          If JOB = 'N' or 'V', S is not referenced.
C     
C  SEP     (global output) DOUBLE PRECISION
C          If JOB = 'V' or 'B', SEP is the estimated reciprocal
C          condition number of the specified invariant subspace. If
C          M = 0 or N, SEP = norm(T).
C          If JOB = 'N' or 'E', SEP is not referenced.
C     
C  DWORK   (local workspace/output) DOUBLE PRECISION array, 
C          dimension (LDWORK) 
C          On exit, if INFO = 0, DWORK(1) returns the optimal LDWORK.
C     
C  LDWORK  (local input) INTEGER
C          The dimension of the array DWORK.
C     
C          If LDWORK = -1, then a workspace query is assumed; the routine
C          only calculates the optimal size of the DWORK array, returns
C          this value as the first entry of the DWORK array, and no error
C          message related to LDWORK is issued by PXERBLA.
C     
C  IWORK   (local workspace/output) INTEGER array, dimension (LIWORK)
C     
C  LIWORK  (local input) INTEGER
C          The dimension of the array IWORK.
C     
C          If LIWORK = -1, then a workspace query is assumed; the
C          routine only calculates the optimal size of the IWORK array,
C          returns this value as the first entry of the IWORK array, and
C          no error message related to LIWORK is issued by PXERBLA.
C     
C  INFO    (global output) INTEGER
C          = 0: successful exit
C          < 0: if INFO = -i, the i-th argument had an illegal value.
C          If the i-th argument is an array and the j-entry had
C          an illegal value, then INFO = -(i*1000+j), if the i-th
C          argument is a scalar and had an illegal value, then INFO = -i.
C          > 0: here we have several possibilites
C            *) Reordering of T failed because some eigenvalues are too
C               close to separate (the problem is very ill-conditioned);
C               T may have been partially reordered, and WR and WI
C               contain the eigenvalues in the same order as in T; S and
C               SEP (if requested) are set to zero. The process that 
C               failed in the reordering will return INFO = {the index 
C               of T where the swap failed}; all others will return 
C               INFO = 1.
C            *) A 2-by-2 block to be reordered split into two 1-by-1
C               blocks and the second block failed to swap with an
C               adjacent block. The process that failed in the 
C               reordering will return INFO = {the index of T where the 
C               swap failed}; all others will return INFO = 1.
C            *) If INFO = 2, the routines used in the calculation of the
C               condition numbers raised a positive warning flag (see
C               the documentation for PGESYCTD and PSYCTCON of the
C               SCASY library).
C            *) If INFO = 3, there is no valid BLACS context (see the
C               BLACS documentation for details).
C            *) If INFO = 333, PGESYCTD raised an input error flag;
C               please report this bug to the authors (see below).
C               If INFO = 444, PSYCTCON raised an input error flag;
C               please report this bug to the authors (see below).
C          In a future release this subroutine may distinguish between
C          the case 1 and 2 above.  
C     
C  Method
C  ======
C
C  This routine performs parallel eigenvalue reordering in real Schur
C  form by parallelizing the approach proposed in [3]. The condition
C  number estimation part is performed by using techniques and code
C  from SCASY, see http://www.cs.umu.se/research/parallel/scasy.
C
C  Additional requirements
C  =======================
C
C  The following alignment requirements must hold:
C  (a) DESCT( MB_ ) = DESCT( NB_ ) = DESCQ( MB_ ) = DESCQ( NB_ )
C  (b) DESCT( RSRC_ ) = DESCQ( RSRC_ )
C  (c) DESCT( CSRC_ ) = DESCQ( CSRC_ ) 
C
C  All matrices must be blocked by a block factor larger than or
C  equal to two (3). This to simplify reordering across processor
C  borders in the presence of 2-by-2 blocks.
C
C  Limitations
C  ===========
C
C  This algorithm cannot work on submatrices of T and Q, i.e.,
C  IT = JT = IQ = JQ = 1 must hold. This is however no limitation
C  since PDLAHQR does not compute Schur forms of submatrices anyway.
C
C  References
C  ==========
C
C  [1] Z. Bai and J. W. Demmel; On swapping diagonal blocks in real
C      Schur form, Linear Algebra Appl., 186:73--95, 1993. Also as
C      LAPACK Working Note 54.
C
C  [2] Z. Bai, J. W. Demmel, and A. McKenney; On computing condition
C      numbers for the nonsymmetric eigenvalue problem, ACM Trans.
C      Math. Software, 19(2):202--223, 1993. Also as LAPACK Working
C      Note 13.
C
C  [3] D. Kressner; Block algorithms for reordering standard and
C      generalized Schur forms, ACM TOMS, 32(4):521-532, 2006. 
C      Also LAPACK Working Note 171.
C
C  [4] R. Granat, B. Kågström, D. Kressner; Parallel eigenvalue 
C      reordering in real Schur form, submitted to Concurrency and
C      Computations: Practice & Experience, Septemeber, 2007. Also as 
C      LAPACK Working Note ???.
C
C  Parallel execution recommendations
C  ==================================
C 
C  Use a square grid, if possible, for maximum performance. The block
C  parameters in PARA should be kept well below the data distribution
C  block size. In particular, see [3,4] for recommended settings for
C  these parameters.
C  
C  In general, the parallel algorithm strives to perform as much work 
C  as possible without crossing the block borders on the main block 
C  diagonal.
C
C  Contributors
C  ============
C  
C  Implemented by Robert Granat, Umea University and HPC2N, March 2007,
C  in collaboration with Bo Kågström and Daniel Kressner.
C 
C  Revisions
C  =========
C
C  Please send bug-reports to granat@cs.umu.se
C
C  Keywords
C  ========
C
C  Real Schur form, eigenvalue reordering, Sylvester matrix equation
C
C  =====================================================================
C     ..
C     .. Parameters ..
      CHARACTER          TOP
      INTEGER            BLOCK_CYCLIC_2D, CSRC_, CTXT_, DLEN_, DTYPE_,
     $                   LLD_, MB_, M_, NB_, N_, RSRC_
      LOGICAL            DBG, DBG2, DBG3, DBG4, DBG5
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( TOP = '1-Tree',
     $                     BLOCK_CYCLIC_2D = 1, DLEN_ = 9, DTYPE_ = 1,
     $                     CTXT_ = 2, M_ = 3, N_ = 4, MB_ = 5, NB_ = 6,
     $                     RSRC_ = 7, CSRC_ = 8, LLD_ = 9,
     $                     ZERO = 0.0D+0, ONE = 1.0D+0, DBG = .FALSE.,
     $                     DBG2 = .FALSE., DBG3 = .FALSE., 
     $                     DBG4 = .FALSE., DBG5 = .FALSE. ) 
C     ..
C     .. Local Scalars ..
      LOGICAL            LQUERY, PAIR, SWAP, WANTBH, WANTQ, WANTS,
     $                   WANTSP, UPDATE, ISHH, FIRST, SKIP1CR, BORDER,
     $                   LASTWAIT
      INTEGER            NPROW, NPCOL, MYROW, MYCOL, NB, NPROCS, 
     $                   IERR, KKK, LLL, DB, DIM1, INDX, IDUM1,
     $                   IDUM2, LLDT, TRSRC, TCSRC, ILOC1, JLOC1,
     $                   ICOFFT12, T12ROWS, T12COLS, SPACE, MYIERR,
     $                   NOEXSY, IPIW1, IPIW2, N1N2, ITER, ICTXT, 
     $                   RSRC1, CSRC1, ILOC3, JLOC3, TRSRC3, 
     $                   TCSRC3, ILOC, JLOC, TRSRC4, TCSRC4, IIT,
     $                   FLOPS, I, ILO, IHI, J, K, KASE, KK, KKS, 
     $                   KS, LIWMIN, LWMIN, MMULT, MMWORK, N1, N2,
     $                   NCB, NDTRAF, NITRAF, NN, NWIN, NUMWIN, PDTRAF,
     $                   PITRAF, PDW, SEL, WINEIG, WINSIZ, LLDQ, RINDX,
     $                   CINDX, RSRC, CSRC, ILILO, ILIHI, ILSEL, IRSRC,
     $                   ICSRC, IPIW, IPW1, IPW2, IPW3, TIHI, TILO,
     $                   NWIN2, ISEL, LIHI, WINDOW, LILO, LSEL, BUFFER,
     $                   NMWIN2, BUFFLEN, LROWS, LCOLS, ILOC2, JLOC2,
     $                   WNEICR, WINDOW0, RSRC4, CSRC4, LIHI4, RSRC3,
     $                   CSRC3, RSRC2, CSRC2, LIHIC, LIHI1, ILEN4,
     $                   SELI4, ILEN1, DIM4, IPW4, QROWS, TROWS,
     $                   TCOLS, IPW5, IPW6, IPW7, IPW8, JLOC4,
     $                   EAST, WEST, ILOC4, SOUTH, NORTH, INDXS,
     $                   ITT, JTT, WRK1, IWRK1, WRK2, IWRK2, WRK3,
     $                   IWRK3, ILEN, DLEN, INDXE, TRSRC1, TCSRC1,
     $                   TRSRC2, TCSRC2, ILOS, TMPNWIN, IPAIR, DIR,
     $                   TLIHI, TLILO, TLSEL, ROUND, LAST, WIN0S,
     $                   WIN0E, WINE, MMAX, MMIN
      DOUBLE PRECISION   EST, RNORM, SCALE, ELEM, DPDUM1, ELEM1, ELEM2,
     $                   ELEM3, ELEM4, SN, CS, TMP, ELEM5, TIME1,
     $                   TIME2, TIME, TIME3, TIME4
C     ..
C     .. Local Arrays ..
      INTEGER            DESCT12( DLEN_ ), DESCXM( DLEN_ ),
     $                   MBNB2( 2 ), IBUFF( 8 )
C     ..
C     .. External Functions ..
      LOGICAL            LSAME, INT2LG
      INTEGER            NUMROC, INDXG2P, INDXG2L, ILG2NT
      DOUBLE PRECISION   PDLANGE, MPI_WTIME
      EXTERNAL           LSAME, PDLANGE, NUMROC, INDXG2P, INDXG2L,
     $                   MPI_WTIME, INT2LG, ILG2NT
C     ..
C     .. External Subroutines ..
      EXTERNAL           PSYCTCON, PDLACPY, PGESYCTD, PXERBLA,
     $                   DGEMX, DLACPY, ILACPY, CHK1MAT, CHK2MAT, 
     $                   INFOG2L, DGSUM2D, DGESD2D, DGERV2D, DGEBS2D, 
     $                   DGEBR2D, IGSUM2D, BLACS_GRIDINFO, IGEBS2D, 
     $                   IGEBR2D, IGAMX2D, IGAMN2D
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, SQRT, MIN
C     ..
C     .. Local Functions .. 
      INTEGER            ICEIL
C     .. 
C     .. Executable Statements ..
C
C     Get grid parameters 
C     
      ICTXT = DESCT( CTXT_ )
      CALL BLACS_GRIDINFO( ICTXT, NPROW, NPCOL, MYROW, MYCOL )
      NPROCS = NPROW*NPCOL
      IF(DBG) WRITE(*,*)MYROW,MYCOL, 'Entering PBDTRSEN'
      IF(DBG) WRITE(*,*)MYROW,MYCOL, 'NPROW,NPCOL,NPROCS=',NPROW,NPCOL,
     $     NPROCS
      IF(DBG) WRITE(*,*)MYROW,MYCOL, 'N,NB=',N,DESCT(MB_)
C     
C     Test if grid is O.K., i.e., the context is valid
C     
      INFO = 0
      IF( NPROW.EQ.-1 ) THEN
         INFO = 3
      END IF
C
C     Check if workspace
C
      LQUERY = LDWORK.EQ.-1 .OR. LIWORK.EQ.-1
C     
C     Test dimensions for local sanity
C     
      IF( INFO.EQ.0 ) THEN
         CALL CHK1MAT( N, 5, N, 5, IT, JT, DESCT, 9, INFO )
      END IF
      IF( INFO.EQ.0 ) THEN
         CALL CHK1MAT( N, 5, N, 5, IQ, JQ, DESCQ, 13, INFO )
      END IF
C     
C     Check the blocking sizes for alignment requirements
C     
      IF( INFO.EQ.0 ) THEN
         IF( DESCT( MB_ ).NE.DESCT( NB_ ) ) INFO = -(1000*9 + MB_)
      END IF
      IF( INFO.EQ.0 ) THEN
         IF( DESCQ( MB_ ).NE.DESCQ( NB_ ) ) INFO = -(1000*13 + MB_)
      END IF
      IF( INFO.EQ.0 ) THEN
         IF( DESCT( MB_ ).NE.DESCQ( MB_ ) ) INFO = -(1000*9 + MB_)
      END IF
C     
C     Check the blocking sizes for minimum sizes
C
      IF( INFO.EQ.0 ) THEN
         IF( N.NE.DESCT( MB_ ) .AND. DESCT( MB_ ).LT.3 ) 
     $        INFO = -(1000*9 + MB_)
         IF( N.NE.DESCQ( MB_ ) .AND. DESCQ( MB_ ).LT.3 )
     $        INFO = -(1000*13 + MB_)
      END IF
C     
C     Check parameters in PARA  
C
      NB = DESCT( MB_ )
      IF( INFO.EQ.0 ) THEN
         IF( PARA(1).LT.1 .OR. PARA(1).GT.MIN(NPROW,NPCOL) ) 
     $        INFO = -(1000 * 4 + 1)
         IF( PARA(2).LT.1 .OR. PARA(2).GE.PARA(3) ) 
     $        INFO = -(1000 * 4 + 2)
         IF( PARA(3).LT.1 .OR. PARA(3).GT.NB ) 
     $        INFO = -(1000 * 4 + 3)
         IF( PARA(4).LT.0 .OR. PARA(4).GT.100 ) 
     $        INFO = -(1000 * 4 + 4)
         IF( PARA(5).LT.1 .OR. PARA(5).GT.NB )
     $        INFO = -(1000 * 4 + 5)  
         IF( PARA(6).LT.1 .OR. PARA(6).GT.PARA(2) )
     $        INFO = -(1000 * 4 + 6)  
      END IF
C     
C     Check requirements on IT, JT, IQ and JQ
C     
      IF( INFO.EQ.0 ) THEN
         IF( IT.NE.1 ) INFO = -7
         IF( JT.NE.IT ) INFO = -8
         IF( IQ.NE.1 ) INFO = -11
         IF( JQ.NE.IQ ) INFO = -12
      END IF
      IF(DBG) WRITE(*,*)MYROW,MYCOL,'Test input'
C     
C     Test input parameters for global sanity
C     
      IF( INFO.EQ.0 ) THEN
         IF(DBG) WRITE(*,*)MYROW,MYCOL, 'PCHK1MAT 1'
         CALL PCHK1MAT( N, 5, N, 5, IT, JT, DESCT, 9, 0, IDUM1, 
     $        IDUM2, INFO )
      END IF
      IF( INFO.EQ.0 ) THEN
         IF(DBG) WRITE(*,*)MYROW,MYCOL, 'PCHK1MAT 2'
         CALL PCHK1MAT( N, 5, N, 5, IQ, JQ, DESCQ, 13, 0, IDUM1, 
     $        IDUM2, INFO )
      END IF
      IF( INFO.EQ.0 ) THEN
         IF(DBG) WRITE(*,*)MYROW,MYCOL, 'PCHK2MAT'
         CALL PCHK2MAT( N, 5, N, 5, IT, JT, DESCT, 9, N, 5, N, 5, 
     $        IQ, JQ, DESCQ, 13, 0, IDUM1, IDUM2, INFO )
      END IF
C     
C     Decode and test the input parameters
C     
      IF(DBG) WRITE(*,*)MYROW,MYCOL, 'Decode and test input param.'
      IF( INFO.EQ.0 .OR. LQUERY ) THEN
         WANTBH = LSAME( JOB, 'B' )
         WANTS = LSAME( JOB, 'E' ) .OR. WANTBH
         WANTSP = LSAME( JOB, 'V' ) .OR. WANTBH
         WANTQ = LSAME( COMPQ, 'V' )
C    
         IF( .NOT.LSAME( JOB, 'N' ) .AND. .NOT.WANTS .AND. .NOT.WANTSP )
     $        THEN
            INFO = -1
         ELSEIF( .NOT.LSAME( COMPQ, 'N' ) .AND. .NOT.WANTQ ) THEN
            INFO = -2
         ELSEIF( N.LT.0 ) THEN
            INFO = -4
         ELSE
C     
C     Extract local leading dimension
C     
            LLDT = DESCT( LLD_ )
            LLDQ = DESCQ( LLD_ )
C     
C     Check the SELECT vector for consistency and set M to the dimension 
C     of the specified invariant subspace.
C     
            IF(DBG) WRITE(*,*)MYROW,MYCOL, 'SELECT=',SELECT(1:N)
            IF(DBG) WRITE(*,*)MYROW,MYCOL, 'Calculating M'
            M = 0
            DO 10 K = 1, N
               IF( K.LT.N ) THEN
                  CALL INFOG2L( K+1, K, DESCT, NPROW, NPCOL, 
     $                 MYROW, MYCOL, ITT, JTT, TRSRC, TCSRC )
                  IF( MYROW.EQ.TRSRC .AND. MYCOL.EQ.TCSRC ) THEN
                     ELEM = T( (JTT-1)*LLDT + ITT )
                     IF( ELEM.NE.ZERO ) THEN
                        IF( DBG ) WRITE(*,*) MYROW,MYCOL,'K,ELEM=',
     $                       K,ELEM
                        IF( INT2LG(SELECT(K)).AND.
     $                       .NOT.INT2LG(SELECT(K+1)) ) THEN
                           IF(DBG) WRITE(*,*) MYROW,MYCOL,
     $                          'SELECT ERROR 1 K=',K
C                           INFO = -3
                           SELECT(K+1) = 1
                        ELSEIF( .NOT.INT2LG(SELECT(K)).AND.
     $                          INT2LG(SELECT(K+1)) ) THEN
                           IF(DBG) WRITE(*,*) MYROW,MYCOL,
     $                          'SELECT ERROR 2 K=',K
C                           INFO = -3
                           SELECT(K) = 1
                        END IF
                     END IF
                  END IF
               END IF
               IF( INT2LG(SELECT(K)) ) M = M + 1
               IF(DBG) WRITE(*,*)MYROW,MYCOL,'M=',M 
 10         CONTINUE
            MMAX = M
            MMIN = M
            IF( NPROCS.GT.1 ) 
     $           CALL IGAMX2D( ICTXT, 'All', TOP, 1, 1, MMAX, 1, -1, 
     $           -1, -1, -1, -1 )
            IF( NPROCS.GT.1 ) 
     $           CALL IGAMN2D( ICTXT, 'All', TOP, 1, 1, MMIN, 1, -1, 
     $           -1, -1, -1, -1 )
            IF( MMAX.GT.MMIN ) THEN
               M = MMAX
               IF( NPROCS.GT.1 ) 
     $              CALL IGAMX2D( ICTXT, 'All', TOP, N, 1, SELECT, N, 
     $              -1, -1, -1, -1, -1 )
            END IF
            IF(DBG5) WRITE(*,*)MYROW,MYCOL, 'Done Calculating M=',M
            IF(DBG) WRITE(*,*) 'INFO=',INFO
C     
C     Set parameters for deep pipelining in parallel Sylvester solver
C
            MBNB2( 1 ) = MIN( MAX( PARA( 3 ), PARA( 2 )*2 ), NB )
            MBNB2( 2 ) = MBNB2( 1 )
C     
C     Compute needed workspace
C     
            IF(DBG) WRITE(*,*) MYROW,MYCOL,'computing workspace'
            N1 = M
            N2 = N - M
            IF( WANTS ) THEN
               CALL PGESYCTD( 'Solve', 'Schur', 'Schur', 'Notranspose', 
     $              'Notranspose', -1, 'Demand', N1, N2, T, 1, 1, DESCT, 
     $              T, N1+1, N1+1, DESCT, T, 1, N1+1, DESCT, MBNB2, 
     $              DWORK, -1, IWORK, -1, NOEXSY, SCALE, IERR )
               WRK1 = INT(DWORK(1))
               IWRK1 = IWORK(1)
            ELSE
               WRK1 = 0
               IWRK1 = 0
            END IF
C
            IF( WANTSP ) THEN
               CALL PSYCTCON( 'Notranspose', 'Notranspose', -1, 
     $              'Demand', N1, N2, T, 1, 1, DESCT, T, N1+1, N1+1, 
     $              DESCT, MBNB2, DWORK, -1, IWORK, -1, EST, ITER, 
     $              IERR )
               WRK2 = INT(DWORK(1))
               IWRK2 = IWORK(1)
            ELSE
               WRK2 = 0
               IWRK2 = 0
            END IF
C
            TROWS = NUMROC( N, NB, MYROW, DESCT(RSRC_), NPROW )
            TCOLS = NUMROC( N, NB, MYCOL, DESCT(CSRC_), NPCOL )
            WRK3 = N + 7*NB**2 + 2*TROWS*PARA( 3 ) + TCOLS*PARA( 3 ) +
     $           MAX( TROWS*PARA( 3 ), TCOLS*PARA( 3 ) )
            IWRK3 = 5*PARA( 1 ) + PARA(2)*PARA(3) - 
     $           PARA(2) * (PARA(2) + 1 ) / 2
C
            IF( WANTSP ) THEN
               LWMIN = MAX( 1, MAX( WRK2, WRK3) )
               LIWMIN = MAX( 1, MAX( IWRK2, IWRK3 ) )
            ELSE IF( LSAME( JOB, 'N' ) ) THEN
               LWMIN = MAX( 1, WRK3 )
               LIWMIN = IWRK3
            ELSE IF( LSAME( JOB, 'E' ) ) THEN
               LWMIN = MAX( 1, MAX( WRK1, WRK3) )
               LIWMIN = MAX( 1, MAX( IWRK1, IWRK3 ) )
            END IF
C     
            IF( LDWORK.LT.LWMIN .AND. .NOT.LQUERY ) THEN
               INFO = -20
            ELSE IF( LIWORK.LT.LIWMIN .AND. .NOT.LQUERY ) THEN
               INFO = -22
            END IF
            IF(DBG) WRITE(*,*) MYROW,MYCOL,'done computing workspace'
         END IF
      END IF
C     
C     Global maximum on info
C     
      IF(DBG) WRITE(*,*) MYROW,MYCOL,'Global max on info NPROCS=',
     $     NPROCS
      IF( NPROCS.GT.1 )
     $     CALL IGAMX2D( ICTXT, 'All', TOP, 1, 1, INFO, 1, -1, -1, -1, 
     $     -1, -1 )
C     
C     Return if some argument is incorrect
C     
      IF(DBG) WRITE(*,*) MYROW,MYCOL,'return if argument incorrect'
      IF( INFO.NE.0 .AND. .NOT.LQUERY ) THEN
         M = 0
         S = ONE
         SEP = ZERO
         CALL PXERBLA( ICTXT, 'PBDTRSEN', -INFO )
         RETURN
      ELSEIF( LQUERY ) THEN
         DWORK( 1 ) = DBLE(LWMIN)
         IWORK( 1 ) = LIWMIN
         RETURN
      END IF
C     
C     Quick return if possible.
C    
      IF(DBG) WRITE(*,*) MYROW,MYCOL,'quick return if possible' 
      IF( M.EQ.N .OR. M.EQ.0 ) THEN
         IF( WANTS )
     $        S = ONE
         IF( WANTSP )
     $        SEP = PDLANGE( '1', N, N, T, IT, JT, DESCT, DWORK )
         GO TO 545
      END IF
C     
C     Set parameters
C     
      IF(DBG) WRITE(*,*) MYROW,MYCOL,'set parameters'
      NUMWIN = PARA( 1 )
      WINEIG = MAX( PARA( 2 ), 2 )
      WINSIZ = MIN( MAX( PARA( 3 ), PARA( 2 )*2 ), NB )
      MMULT  = PARA( 4 )
      NCB    = PARA( 5 )
      WNEICR = PARA( 6 )
C     
C     Insert some pointers into INTEGER workspace
C     
      IF(DBG) WRITE(*,*) MYROW,MYCOL,'insert pointers'
      ILILO = 1
      ILIHI = ILILO + NUMWIN
      ILSEL = ILIHI + NUMWIN
      IRSRC = ILSEL + NUMWIN
      ICSRC = IRSRC + NUMWIN
      IPIW  = ICSRC + NUMWIN
C     
C     Insert some pointers into DOUBLE PRECISION workspace - for now we
C     only need two pointers
C     
      IPW1 = 1
      IPW2 = IPW1 + NB
C     
C     Collect the selected blocks at the top-left corner of T.
C     
C     Globally: ignore eigenvalues that are already in order.
C     ILO is a global variable and is kept updated to be consistent
C     throughout the process mesh.
C     
      ILO = 0
 40   CONTINUE
      ILO = ILO + 1
      IF( ILO.LE.N ) THEN
         IF( INT2LG(SELECT(ILO)) ) GO TO 40
      END IF
      IF(DBG) WRITE(*,*)MYROW,MYCOL, 'ILO=',ILO
C     
C     Globally: start the collection at the top of the matrix. Here,
C     IHI is a global variable and is kept updated to be consistent
C     throughout the process mesh.
C     
      IHI = N
C
C     Zero out local timers
C
      TIME  = ZERO
      TIME1 = ZERO
      TIME2 = ZERO
      TIME3 = ZERO
      TIME4 = ZERO
C     
C     Globally:  While ( ILO <= M ) do
 50   CONTINUE
      IF(DBG3) WRITE(*,*)MYROW,MYCOL, 'back at 50'
C     
      IF( ILO.LE.M ) THEN
C
C     Depending on the value of ILO, find the diagonal block index J,
C     such that T(1+(J-1)*NB:1+J*NB,1+(J-1)*NB:1+J*NB) contains the
C     first unsorted eigenvalue. Check that J does not point to a block
C     with only one selected eigenvalue in the last position which
C     belongs to a splitted 2-by-2 block.
C
         ILOS = ILO - 1
 52      CONTINUE
         ILOS = ILOS + 1
         IF( .NOT. INT2LG(SELECT(ILOS)) ) GO TO 52
         IF( ILOS.LT.N ) THEN
            IF( INT2LG(SELECT( ILOS+1 )) .AND. MOD(ILOS,NB).EQ.0 ) THEN
               CALL PDELGET( 'All', TOP, ELEM, T, ILOS+1, ILOS, DESCT )
               IF( ELEM.NE.ZERO ) GO TO 52
            END IF
         END IF
         J = ICEIL(ILOS,NB)
         IF(DBG3) WRITE(*,*)MYROW,MYCOL, 
     $        'first block with non-selected eigenvalues=',J
C     
C     Globally: Set startvalues of LILO and LIHI for all processes. 
C     Choose also the number of selected eigenvalues at top of each
C     diagonal block such that the number of eigenvalues which remain 
C     to be reordered is an integer multiple of WINEIG. 
C     
C     All the information is saved into the INTEGER workspace such that 
C     all processors are aware of each others operations.  
C
C     Compute the number of concurrent windows
C     
         NMWIN2 = (ICEIL(IHI,NB)*NB - (ILO-MOD(ILO,NB)+1)+1) / NB
         NMWIN2 = MIN( MIN( NUMWIN, NMWIN2 ), ICEIL(N,NB) - J + 1 )
         IF(DBG3) WRITE(*,*)MYROW,MYCOL, 
     $        'NMWIN2 =',NMWIN2
C     
C     For all windows, set LSEL = 0 and find a proper start value of
C     LILO such that LILO points at the first non-selected entry in
C     the corresponding diagonal block of T
C
         DO 80 K = 1, NMWIN2
            IWORK( ILSEL+K-1) = 0
            IWORK( ILILO+K-1) = MAX( ILO, (J-1)*NB+(K-1)*NB+1 )
            LILO = IWORK( ILILO+K-1 )
 82         CONTINUE
            IF( INT2LG(SELECT(LILO)) .AND. LILO.LT.(J+K-1)*NB ) THEN
               LILO = LILO + 1
               IF( LILO.LE.N ) GO TO 82
            END IF
            IWORK( ILILO+K-1 ) = LILO
C
C     Fix each LILO to avoid that no 2-by-2 block is cut in top of the 
C     submatrix (LILO:LIHI,LILO:LIHI)
C
            LILO = IWORK(ILILO+K-1)
            IF(DBG3) WRITE(*,*)MYROW,MYCOL, 'LILO=',LILO
            IF( LILO.GT.NB ) THEN
               CALL PDELGET( 'All', TOP, ELEM, T, LILO, LILO-1, DESCT ) 
               IF( ELEM.NE.ZERO ) THEN
                  IF( LILO.LT.(J+K-1)*NB ) THEN
                     IWORK(ILILO+K-1) = IWORK(ILILO+K-1) + 1
                  ELSE
                     IWORK(ILILO+K-1) = IWORK(ILILO+K-1) - 1
                  END IF 
               END IF
            END IF
C
C     Set a proper LIHI value for each window. Also find the processors
C     corresponding to the corresponding windows.
C
            IWORK( ILIHI+K-1 ) =  IWORK( ILILO+K-1 )
            IF(DBG3) WRITE(*,*)MYROW,MYCOL,'K,LIHI1=',
     $           K,IWORK( ILIHI+K-1 )
            IWORK( IRSRC+K-1 ) = INDXG2P( IWORK(ILILO+K-1), NB, MYROW, 
     $           DESCT( RSRC_ ), NPROW )
            IWORK( ICSRC+K-1 ) = INDXG2P( IWORK(ILILO+K-1), NB, MYCOL, 
     $           DESCT( CSRC_ ), NPCOL )
            TILO = IWORK(ILILO+K-1)
            TIHI = MIN( N, ICEIL( TILO, NB ) * NB )
            DO 90 KK = TIHI, TILO, -1
               IF( INT2LG(SELECT( KK )) ) THEN
                  IWORK(ILIHI+K-1) = MAX(IWORK(ILIHI+K-1) , KK )
                  IF(DBG3) WRITE(*,*)MYROW,MYCOL,'K,LIHI2=',
     $                 K,IWORK( ILIHI+K-1 )
                  IWORK(ILSEL+K-1) = IWORK(ILSEL+K-1) + 1
                  IF( IWORK(ILSEL+K-1).GT.WINEIG ) THEN
                     IWORK(ILIHI+K-1) = KK
                     IF(DBG3) WRITE(*,*)MYROW,MYCOL,'K,LIHI3=',
     $                    K,IWORK( ILIHI+K-1 )
                     IWORK(ILSEL+K-1) = 1
                  END IF
               END IF
 90         CONTINUE
C     
C     Fix each LIHI to avoid that bottom of window cuts 2-by-2 block. 
C     We exclude such a block if located on block (process) border and 
C     on window border or if an inclusion would cause violation on the 
C     maximum number of eigenvalues to reorder inside each window. If 
C     only on window border, we include it. The excluded block is 
C     included automatically later when a subcluster is reordered into 
C     the block from South-East.
C     
            LIHI = IWORK(ILIHI+K-1)
            IF( LIHI.LT.N ) THEN
               CALL PDELGET( 'All', TOP, ELEM, T, LIHI+1, LIHI, DESCT ) 
               IF( ELEM.NE.ZERO ) THEN
                  IF( ICEIL( LIHI, NB ) .NE. ICEIL( LIHI+1, NB ) .OR. 
     $                 IWORK( ILSEL+K-1 ).EQ.WINEIG ) THEN
                     IWORK( ILIHI+K-1 ) = IWORK( ILIHI+K-1 ) - 1
                     IF(DBG3) WRITE(*,*)MYROW,MYCOL,'K,LIHI4=',
     $                    K,IWORK( ILIHI+K-1 )
                     IF( IWORK( ILSEL+K-1 ).GT.2 )
     $                    IWORK( ILSEL+K-1 ) = IWORK( ILSEL+K-1 ) - 1
                  ELSE
                     IWORK( ILIHI+K-1 ) = IWORK( ILIHI+K-1 ) + 1
                     IF(DBG3) WRITE(*,*)MYROW,MYCOL,'K,LIHI5=',
     $                    K,IWORK( ILIHI+K-1 )
                     IF( INT2LG(SELECT(LIHI+1)) )
     $                    IWORK( ILSEL+K-1 ) = IWORK( ILSEL+K-1 ) + 1
                  END IF
               END IF
            END IF
 80      CONTINUE
C
C     Fix the special cases of LSEL = 0 and LILO = LIHI for each window 
C     by assuring that the stop-condition for local reordering is 
C     fulfilled directly. Do this by setting LIHI = startposition for the 
C     corresponding block and LILO = LIHI + 1.
C
         DO 85 K = 1, NMWIN2
            LILO = IWORK( ILILO + K - 1 )
            LIHI = IWORK( ILIHI + K - 1 )
            LSEL = IWORK( ILSEL + K - 1 )
            IF( LSEL.EQ.0 .OR. LILO.EQ.LIHI ) THEN
               LIHI = IWORK( ILIHI + K - 1 )
               IWORK( ILIHI + K - 1 ) = (ICEIL(LIHI,NB)-1)*NB + 1
               IF(DBG3) WRITE(*,*)MYROW,MYCOL,'K,LIHI6=',
     $              K,IWORK( ILIHI+K-1 )
               IWORK( ILILO + K - 1 ) = IWORK( ILIHI + K - 1 ) + 1
            END IF
 85      CONTINUE
C
         IF( DBG )CALL BLACS_BARRIER( ICTXT, 'All' )
         IF(DBG) WRITE(*,*) MYROW,MYCOL,
     $        'IWORK(IRSRC:IRSRC+NMWIN2-1)=',
     $        IWORK(IRSRC:IRSRC+NMWIN2-1)
         IF( DBG ) CALL BLACS_BARRIER( ICTXT, 'All' )
         IF(DBG) WRITE(*,*) MYROW,MYCOL,
     $        'IWORK(ICSRC:ICSRC+NMWIN2-1)=',
     $        IWORK(ICSRC:ICSRC+NMWIN2-1)
         IF( DBG ) CALL BLACS_BARRIER( ICTXT, 'All' )
         IF(DBG) WRITE(*,*) MYROW,MYCOL,
     $        'IWORK(ILIHI:ILIHI+NMWIN2-1)=',
     $        IWORK(ILIHI:ILIHI+NMWIN2-1)
         IF( DBG ) CALL BLACS_BARRIER( ICTXT, 'All' )
         IF(DBG) WRITE(*,*) MYROW,MYCOL,
     $        'IWORK(ILILO:ILILO+NMWIN2-1)=',
     $        IWORK(ILILO:ILILO+NMWIN2-1)
         IF( DBG ) CALL BLACS_BARRIER( ICTXT, 'All' )
         IF(DBG) WRITE(*,*) MYROW,MYCOL,
     $        'IWORK(ILSEL:ILSEL+NMWIN2-1)=',
     $        IWORK(ILSEL:ILSEL+NMWIN2-1)
         IF( DBG ) CALL BLACS_BARRIER( ICTXT, 'All' )
C
C     Associate all processors with the first computational window 
C     that should be activated, if possible. 
C
         LILO = IHI
         LIHI = ILO
         LSEL = M
         FIRST = .TRUE.
         DO 95 WINDOW = 1, NMWIN2
            RSRC = IWORK(IRSRC+WINDOW-1)
            CSRC = IWORK(ICSRC+WINDOW-1)
            IF( MYROW.EQ.RSRC .OR. MYCOL.EQ.CSRC ) THEN
               TLILO = IWORK( ILILO + WINDOW - 1 )
               TLIHI = IWORK( ILIHI + WINDOW - 1 )
               TLSEL = IWORK( ILSEL + WINDOW - 1 )
               IF( (.NOT. ( LIHI .GE. LILO + LSEL ) ) .AND.
     $             ( (TLIHI .GE. TLILO + TLSEL) .OR. FIRST ) ) THEN
                  IF( FIRST ) FIRST = .FALSE.
                  LILO = TLILO
                  LIHI = TLIHI
                  LSEL = TLSEL
                  GO TO 97
               END IF
            END IF
 95      CONTINUE
 97      CONTINUE
C
C     Exclude all processors that are not involved in any computational
C     window right now
C
         IF( LILO.EQ.IHI .AND. LIHI.EQ.ILO .AND. LSEL.EQ.M ) THEN
            IF( DBG3 )WRITE(*,*) MYROW,MYCOL,
     $           'excluded from current set of windows'
            GO TO 114
         END IF
C
C     Make sure all processors associated with a compuational window
C     enter the local reordering the first time 
C
         FIRST = .TRUE.
         IF(DBG3) WRITE(*,*)MYROW,MYCOL, 'entering local reordering...'
         TIME = MPI_WTIME()
C     
C     Globally for all computational windows:
C     While ( LIHI >= LILO + LSEL ) do
         ROUND = 1
 130     CONTINUE
         IF(DBG2) WRITE(*,*) MYROW,MYCOL,'LIHI,LILO,LSEL=',LIHI,LILO,
     $        LSEL
         IF(DBG) WRITE(*,*) MYROW,MYCOL,'LIHI .GE. LILO + LSEL =',
     $        LIHI .GE. LILO + LSEL
         IF( FIRST .OR. ( LIHI .GE. LILO + LSEL ) ) THEN
            IF(DBG2) WRITE(*,*) MYROW,MYCOL,' ROUND =',ROUND
            IF( ROUND.GT.20 ) THEN
               IF(DBG5) WRITE(*,*) MYROW,MYCOL,'ROUND>20 -> STOP!'
               IF(DBG5) WRITE(*,*) N, NB, M
               STOP
            END IF
C     
C     Perform computations in parallel: loop thorugh all compuational
C     windows, do local reordering and accumulate transformations,
C     broadcast them in the corresponding block row and columns and
C     compute the corresponding updates
C     
            IF(DBG2) WRITE(*,*) MYROW,MYCOL,'local reordering'
            DO 110 WINDOW = 1, NMWIN2
               RSRC = IWORK(IRSRC+WINDOW-1)
               CSRC = IWORK(ICSRC+WINDOW-1)
C     
C     The process on the block diagonal computes the reordering 
C     
               IF( MYROW.EQ.RSRC .AND. MYCOL.EQ.CSRC ) THEN
                  LILO = IWORK(ILILO+WINDOW-1)
                  LIHI = IWORK(ILIHI+WINDOW-1)
                  LSEL = IWORK(ILSEL+WINDOW-1)
C
C     Compute the local value of I
C 
                  I = MAX( LILO, LIHI - WINSIZ + 1 )
C     
C     Fix my I to avoid that top of window cuts 2-by-2 block.
C     
                  IF( I.GT.LILO ) THEN
                     CALL INFOG2L( I, I-1, DESCT, NPROW, NPCOL, MYROW,
     $                    MYCOL, ILOC, JLOC, RSRC, CSRC )
                     IF( T( LLDT*(JLOC-1) + ILOC ).NE.ZERO )   
     $                    I = I + 1
                  END IF
C
C     Compute local indicies for submatrix to operate on
C
                  CALL INFOG2L( I, I, DESCT, NPROW, NPCOL, 
     $                 MYROW, MYCOL, ILOC1, JLOC1, RSRC, 
     $                 CSRC )
C     
C     The active window is ( I:LIHI, I:LIHI ). Reorder eigenvalues
C     within this window and pipeline transformations.
C     
                  IF(DBG5) WRITE(*,*)MYROW,MYCOL, 'I,LIHI=',I,LIHI
                  IF(DBG5) WRITE(*,*)MYROW,MYCOL, 'SELECT(I:LIHI)=',
     $                 SELECT(I:LIHI)
                  NWIN = LIHI - I + 1
                  KS = 0
                  PITRAF = IPIW
                  PDTRAF = IPW2
C     
                  PAIR = .FALSE.
                  DO 140 K = I, LIHI
                     IF( PAIR ) THEN
                        PAIR = .FALSE.
                     ELSE
                        SWAP = INT2LG(SELECT( K ))
                        IF( K.LT.LIHI ) THEN
                           CALL INFOG2L( K+1, K, DESCT, NPROW, NPCOL, 
     $                          MYROW, MYCOL, ILOC, JLOC, RSRC, CSRC )
                           IF( T( LLDT*(JLOC-1) + ILOC ).NE.ZERO )
     $                          PAIR = .TRUE.
                        END IF
                        IF( SWAP ) THEN
                           KS = KS + 1
C     
C     Swap the K-th block to position I+KS-1.
C     
                           IERR = 0
                           KK  = K - I + 1
                           KKS = KS
                           IF( KK.NE.KS ) THEN
                              NITRAF = LIWORK - PITRAF + 1
                              NDTRAF = LDWORK - PDTRAF + 1
                              IF( DBG4 ) WRITE(*,*) 'BDTREXC 1'
                              CALL BDTREXC( NWIN, 
     $                             T(LLDT*(JLOC1-1) + ILOC1), LLDT, KK, 
     $                             KKS, NITRAF, IWORK( PITRAF ), NDTRAF, 
     $                             DWORK( PDTRAF ), DWORK(IPW1), IERR )
                              IF( DBG5 ) WRITE(*,*)MYROW,MYCOL, 'IERR=',
     $                             IERR
                              PITRAF = PITRAF + NITRAF
                              PDTRAF = PDTRAF + NDTRAF
C     
C     Update array SELECT.
C     
                              IF ( PAIR ) THEN
                                 DO 150 J = I+KK-1, I+KKS, -1
                                    SELECT(J+1) = SELECT(J-1)
 150                             CONTINUE
                                 SELECT(I+KKS-1) = 1
                                 SELECT(I+KKS) = 1
                              ELSE
                                 DO 160 J = I+KK-1, I+KKS, -1
                                    SELECT(J) = SELECT(J-1)
 160                             CONTINUE
                                 SELECT(I+KKS-1) = 1
                              END IF
C     
                              IF ( IERR.EQ.1 .OR. IERR.EQ.2 ) THEN
C     
C     Some blocks are too close to swap: prepare to leave in a clean 
C     fashion. If IERR.EQ.2, we must update SELECT to account for the 
C     fact that the 2 by 2 block to be reordered did split and the first 
C     part of this block is already reordered.
C     
                                 IF ( IERR.EQ.2 ) THEN
                                    SELECT( I+KKS-3 ) = 1
                                    SELECT( I+KKS-1 ) = 0
                                    KKS = KKS + 1
                                 END IF
C     
C     Update off-diagonal blocks immediately.
C     
                                 GO TO 170
                              END IF
                              KS = KKS     
                           END IF
                           IF( PAIR )
     $                          KS = KS + 1
                        END IF
                     END IF
 140              CONTINUE
               END IF
 110        CONTINUE
 170        CONTINUE
C
C     The on-diagonal processes save their information from the local
C     reordering in the integer buffer. This buffer is broadcasted to
C     updating processors, see below.
C
            DO 175 WINDOW = 1, NMWIN2
               RSRC = IWORK(IRSRC+WINDOW-1)
               CSRC = IWORK(ICSRC+WINDOW-1)
               IF( MYROW.EQ.RSRC .AND. MYCOL.EQ.CSRC ) THEN
                  IBUFF( 1 ) = I
                  IBUFF( 2 ) = NWIN
                  IBUFF( 3 ) = PITRAF
                  IBUFF( 4 ) = KS
                  IBUFF( 5 ) = PDTRAF
                  IBUFF( 6 ) = NDTRAF
                  ILEN = PITRAF - IPIW 
                  DLEN = PDTRAF - IPW2
                  IBUFF( 7 ) = ILEN
                  IBUFF( 8 ) = DLEN
               END IF
 175        CONTINUE
C     
C     For the updates with respect to the local reordering, we organize
C     the updates in two phases where the update "direction" (controlled
C     by the DIR variable below) is first chosen to be the corresponding
C     rows, then the corresponding columns
C
            DO 1111 DIR = 1, 2
            IF(DBG) WRITE(*,*) MYROW,MYCOL,
     $              'DIR =',DIR
C     
C     Broadcast information about the reordering and the accumulated
C     transformations: I, NWIN, PITRAF, NITRAF, PDTRAF, NDTRAF. If no
C     broadcast is performed, use a artificial value of KS to prevent
C     updating indicies for windows already finished (use KS = -1)
C     
            IF(DBG2) WRITE(*,*) MYROW,MYCOL,
     $           'broadcast info and acc trans'
            DO 111 WINDOW = 1, NMWIN2
               RSRC = IWORK(IRSRC+WINDOW-1)
               CSRC = IWORK(ICSRC+WINDOW-1)
               IF( MYROW.EQ.RSRC .OR. MYCOL.EQ.CSRC ) THEN
                  LILO = IWORK(ILILO+WINDOW-1)
                  LIHI = IWORK(ILIHI+WINDOW-1)
                  LSEL = IWORK(ILSEL+WINDOW-1)
               END IF
               IF(DBG) WRITE(*,*) MYROW,MYCOL,
     $              'WINDOW,RSRC,CSRC=',
     $              WINDOW,RSRC,CSRC
               IF( MYROW.EQ.RSRC .AND. MYCOL.EQ.CSRC ) THEN
                  IF(DBG2) WRITE(*,*) MYROW,MYCOL,'send info'
                  IF(DBG) WRITE(*,*) MYROW,MYCOL,'IBUFF1',
     $                    IBUFF(1:8)
                  IF( NPCOL.GT.1 .AND. DIR.EQ.1 )
     $                 CALL IGEBS2D( ICTXT, 'Row', TOP, 8, 1, IBUFF, 8 )
                  IF( NPROW.GT.1 .AND. DIR.EQ.2 )
     $                 CALL IGEBS2D( ICTXT, 'Col', TOP, 8, 1, IBUFF, 8 )
                  IF(DBG2) WRITE(*,*) MYROW,MYCOL,'send info done'
               ELSEIF( MYROW.EQ.RSRC .OR. MYCOL.EQ.CSRC ) THEN
                  IF( NPCOL.GT.1 .AND. DIR.EQ.1 .AND. MYROW.EQ.RSRC ) 
     $                 THEN
                     IF( FIRST .OR. (LIHI .GE. LILO + LSEL) ) THEN 
                        IF(DBG2) WRITE(*,*) MYROW,MYCOL,'recv info row'
                        CALL IGEBR2D( ICTXT, 'Row', TOP, 8, 1, IBUFF, 8, 
     $                       RSRC, CSRC )
                        IF(DBG2) WRITE(*,*) MYROW,MYCOL,'recv info done'
                        IF(DBG) WRITE(*,*) MYROW,MYCOL,'IBUFF1',
     $                       IBUFF(1:8)
                        I = IBUFF( 1 ) 
                        NWIN = IBUFF( 2 ) 
                        PITRAF = IBUFF( 3 ) 
                        KS = IBUFF( 4 ) 
                        PDTRAF = IBUFF( 5 ) 
                        NDTRAF = IBUFF( 6 ) 
                        ILEN = IBUFF( 7 )
                        DLEN = IBUFF( 8 )
                     ELSE
                        ILEN = 0
                        DLEN = 0
                        KS = -1
                     END IF
                  END IF
                  IF( NPROW.GT.1 .AND. DIR.EQ.2 .AND. MYCOL.EQ.CSRC ) 
     $                 THEN
                     IF( FIRST .OR. (LIHI .GE. LILO + LSEL) ) THEN 
                        IF(DBG2) WRITE(*,*) MYROW,MYCOL,'recv info col'
                        CALL IGEBR2D( ICTXT, 'Col', TOP, 8, 1, IBUFF, 8, 
     $                       RSRC, CSRC )
                        IF(DBG2) WRITE(*,*) MYROW,MYCOL,'recv info done'
                        IF(DBG) WRITE(*,*) MYROW,MYCOL,'IBUFF1',
     $                       IBUFF(1:8)
                        I = IBUFF( 1 ) 
                        NWIN = IBUFF( 2 ) 
                        PITRAF = IBUFF( 3 ) 
                        KS = IBUFF( 4 ) 
                        PDTRAF = IBUFF( 5 ) 
                        NDTRAF = IBUFF( 6 )
                        ILEN = IBUFF( 7 )
                        DLEN = IBUFF( 8 ) 
                     ELSE
                        ILEN = 0
                        DLEN = 0
                        KS = -1
                     END IF
                  END IF
               END IF
               IF(DBG) WRITE(*,*) MYROW,MYCOL,'bc info done!'
C     
C     Broadcast the accumulated transformations - copy all information
C     from IWORK(IPIW:PITRAF-1) and DWORK(IPW2:PDTRAF-1) to a buffer and 
C     broadcast this buffer in the corresponding block row and column.
C     On arrival, copy the information back to the correct part of the
C     workspace. This step is avoided if no computations were performed
C     at the diagonal processor, i.e., BUFFLEN = 0.
C     
               IF( MYROW.EQ.RSRC .AND. MYCOL.EQ.CSRC ) THEN
                  BUFFER = PDTRAF
                  BUFFLEN = DLEN + ILEN
                  IF( BUFFLEN.NE.0 ) THEN
                     IF(DBG) WRITE(*,*) MYROW,MYCOL,'converting data'
                     DO 180 INDX = 1, ILEN
                        DWORK( BUFFER+INDX-1 ) = 
     $                       DBLE( IWORK(IPIW+INDX-1) )
 180                 CONTINUE
                     CALL DLACPY( 'All', DLEN, 1, DWORK( IPW2 ), 
     $                    DLEN, DWORK(BUFFER+ILEN), DLEN )
                     IF(DBG) WRITE(*,*) MYROW,MYCOL,'broadcasting data'
                     IF( NPCOL.GT.1 .AND. DIR.EQ.1 ) THEN
                        IF(DBG) WRITE(*,*) MYROW,MYCOL,'bc data row'
                        IF(DBG) WRITE(*,*) MYROW,MYCOL,'BUFFLEN=',
     $                       BUFFLEN
                        CALL DGEBS2D( ICTXT, 'Row', TOP, BUFFLEN, 1, 
     $                       DWORK(BUFFER), BUFFLEN )
                        IF(DBG) WRITE(*,*) MYROW,MYCOL,'bc row done'
                     END IF
                     IF( NPROW.GT.1 .AND. DIR.EQ.2 ) THEN
                        IF(DBG) WRITE(*,*) MYROW,MYCOL,'bc data col'
                        IF(DBG) WRITE(*,*) MYROW,MYCOL,'BUFFLEN=',
     $                       BUFFLEN
                        CALL DGEBS2D( ICTXT, 'Col', TOP, BUFFLEN, 1, 
     $                       DWORK(BUFFER), BUFFLEN )
                        IF(DBG) WRITE(*,*) MYROW,MYCOL,'bc col done'
                     END IF
                  END IF
               ELSEIF( MYROW.EQ.RSRC .OR. MYCOL.EQ.CSRC ) THEN
                  IF(DBG) WRITE(*,*) MYROW,MYCOL,'receving data'
                  IF( NPCOL.GT.1 .AND. DIR.EQ.1 .AND. MYROW.EQ.RSRC ) 
     $                 THEN
                     BUFFER = PDTRAF
                     BUFFLEN = DLEN + ILEN
                     IF( BUFFLEN.NE.0 ) THEN 
                        IF(DBG) WRITE(*,*) MYROW,MYCOL,'recv data row'
                        IF(DBG) WRITE(*,*) MYROW,MYCOL,'BUFFER=',BUFFER
                        IF(DBG) WRITE(*,*) MYROW,MYCOL,'BUFFLEN=',
     $                       BUFFLEN
                        CALL DGEBR2D( ICTXT, 'Row', TOP, BUFFLEN, 1, 
     $                       DWORK(BUFFER), BUFFLEN, RSRC, CSRC )
                        IF(DBG) WRITE(*,*) MYROW,MYCOL,
     $                       'recv data row done'
                     END IF
                  END IF
                  IF( NPROW.GT.1 .AND. DIR.EQ.2 .AND. MYCOL.EQ.CSRC ) 
     $                 THEN
                     BUFFER = PDTRAF
                     BUFFLEN = DLEN + ILEN
                     IF( BUFFLEN.NE.0 ) THEN
                        IF(DBG) WRITE(*,*) MYROW,MYCOL,'recv data col'
                        IF(DBG) WRITE(*,*) MYROW,MYCOL,'BUFFER=',BUFFER
                        IF(DBG) WRITE(*,*) MYROW,MYCOL,'BUFFLEN=',
     $                       BUFFLEN
                        CALL DGEBR2D( ICTXT, 'Col', TOP, BUFFLEN, 1, 
     $                       DWORK(BUFFER), BUFFLEN, RSRC, CSRC )
                        IF(DBG) WRITE(*,*) MYROW,MYCOL,
     $                       'recv data col done'
                     END IF
                  END IF
                  IF((NPCOL.GT.1.AND.DIR.EQ.1.AND.MYROW.EQ.RSRC).OR.
     $               (NPROW.GT.1.AND.DIR.EQ.2.AND.MYCOL.EQ.CSRC ) ) 
     $                 THEN
                     IF( BUFFLEN.NE.0 ) THEN
                        IF(DBG) WRITE(*,*) MYROW,MYCOL,
     $                       'converting data back'
                        DO 190 INDX = 1, ILEN
                           IWORK(IPIW+INDX-1) = 
     $                          INT(DWORK( BUFFER+INDX-1 ))
 190                    CONTINUE
                        CALL DLACPY( 'All', DLEN, 1, 
     $                       DWORK( BUFFER+ILEN ), DLEN, 
     $                       DWORK( IPW2 ), DLEN )
                     END IF
                  END IF
               END IF
 111        CONTINUE
C     
C     Now really perform the updates by applying the orthogonal 
C     transformations to the out-of-window parts of T and Q. This
C     step is avoided if no reordering was performed by the on-
C     diagonal processor from the beginning, i.e., BUFFLEN = 0.
C
C     Count number of operations to decide whether to use matrix-matrix 
C     multiplications for updating off-diagonal parts or not. 
C     
            IF(DBG) WRITE(*,*) MYROW,MYCOL,'count operations'
            DO 112 WINDOW = 1, NMWIN2
               IF(DBG) WRITE(*,*) MYROW,MYCOL,'WINDOW =',WINDOW 
               RSRC = IWORK(IRSRC+WINDOW-1)
               CSRC = IWORK(ICSRC+WINDOW-1)
C
               IF(DBG) WRITE(*,*) MYROW,MYCOL,'RSRC,CSRC =',RSRC,CSRC
               IF( (MYROW.EQ.RSRC .AND. DIR.EQ.1 ).OR.
     $             (MYCOL.EQ.CSRC .AND. DIR.EQ.2 ) ) THEN
                  IF(DBG) WRITE(*,*) MYROW,MYCOL,'On the race...'
                  LILO = IWORK(ILILO+WINDOW-1)
                  LIHI = IWORK(ILIHI+WINDOW-1)
                  LSEL = IWORK(ILSEL+WINDOW-1)
                   IF(DBG) WRITE(*,*) MYROW,MYCOL,'LILO,LIHI,LSEL=',
     $                 LILO,LIHI,LSEL
C
C     Skip update part for current WINDOW if BUFFLEN = 0
C
                  IF( BUFFLEN.EQ.0 ) GO TO 295
C
                  NITRAF = PITRAF - IPIW
                  ISHH = .FALSE.
                  FLOPS = 0
                  DO 200 K = 1, NITRAF
                     IF( IWORK( IPIW + K - 1 ).LE.NWIN ) THEN
                        FLOPS = FLOPS + 6
                     ELSE
                        FLOPS = FLOPS + 11
                        ISHH = .TRUE.
                     END IF
 200              CONTINUE
C     
C     Compute amount of work space necessary for performing matrix-
C     matrix multiplications.
C     
                  PDW = BUFFER
                  IPW3 = PDW + NWIN*NWIN
               ELSE
                  FLOPS = 0
               END IF
C     
               IF( FLOPS.NE.0 .AND.
     $              ( FLOPS*100 ) / ( 2*NWIN*NWIN ) .GE. MMULT ) THEN
C     
C     Update off-diagonal blocks and Q using matrix-matrix
C     multiplications; if there are no Householder reflectors
C     it is preferable to take the triangular block structure
C     of the transformation matrix into account.
C     
                  CALL DLASET( 'All', NWIN, NWIN, ZERO, ONE, 
     $                 DWORK( PDW ), NWIN )
                  IF(DBG) WRITE(*,*) MYROW,MYCOL,'bdlaapp local'
                  CALL BDLAAPP( 1, NWIN, NWIN, NCB, DWORK( PDW ), NWIN,
     $                 NITRAF, IWORK(IPIW), DWORK( IPW2 ), DWORK(IPW3) )
                  IF(DBG) WRITE(*,*) MYROW,MYCOL,'bdlaapp local done'
C     
                  IF( ISHH ) THEN
C     
C     Loop through the local blocks of the distributed matrices T and 
C     Q and update them according to the performed reordering
C     
C     Update the columns of T and Q affected by the reordering
C     
                     IF( DIR.EQ.2 ) THEN
                     DO 210 INDX = 1, I-1, NB
                        CALL INFOG2L( INDX, I, DESCT, NPROW, NPCOL,
     $                       MYROW, MYCOL, ILOC, JLOC, RSRC1, CSRC1 )
                        IF( MYROW.EQ.RSRC1 .AND. MYCOL.EQ.CSRC1 ) THEN
                           LROWS = MIN(NB,I-INDX)
                           IF(DBG)WRITE(*,*)MYROW,MYCOL, 'dgemm 1'
                           CALL DGEMX( 'No transpose', 'No transpose', 
     $                          LROWS, NWIN, NWIN, ONE, 
     $                          T((JLOC-1)*LLDT+ILOC), LLDT, 
     $                          DWORK( PDW ), NWIN, ZERO, DWORK(IPW3), 
     $                          LROWS )
                           CALL DLACPY( 'All', LROWS, NWIN, DWORK(IPW3), 
     $                          LROWS, T((JLOC-1)*LLDT+ILOC), LLDT )
                        END IF
 210                 CONTINUE
                     IF( WANTQ ) THEN
                        DO 220 INDX = 1, N, NB
                           CALL INFOG2L( INDX, I, DESCQ, NPROW, NPCOL,
     $                          MYROW, MYCOL, ILOC, JLOC, RSRC1, CSRC1 )
                           IF( MYROW.EQ.RSRC1.AND.MYCOL.EQ.CSRC1 ) THEN
                              LROWS = MIN(NB,N-INDX+1)
                              IF(DBG) WRITE(*,*)MYROW,MYCOL, 'dgemm 2'
                              CALL DGEMX( 'No transpose', 
     $                             'No transpose', LROWS, NWIN, NWIN, 
     $                             ONE, Q((JLOC-1)*LLDQ+ILOC), LLDQ, 
     $                             DWORK( PDW ), NWIN, ZERO, 
     $                             DWORK(IPW3), LROWS )
                              CALL DLACPY( 'All', LROWS, NWIN, 
     $                             DWORK(IPW3), LROWS, 
     $                             Q((JLOC-1)*LLDQ+ILOC), LLDQ )
                           END IF
 220                    CONTINUE
                     END IF
                     END IF
C     
C     Update the rows of T affected by the reordering
C     
                     IF( DIR.EQ.1 ) THEN
                     IF( LIHI.LT.N ) THEN
                        IF( MOD(LIHI,NB).GT.0 ) THEN
                           INDX = LIHI + 1
                           CALL INFOG2L( I, INDX, DESCT, NPROW, NPCOL,
     $                          MYROW, MYCOL, ILOC, JLOC, RSRC1, CSRC1 )
                           IF( MYROW.EQ.RSRC1.AND.MYCOL.EQ.CSRC1 ) THEN
                              LCOLS = MOD( MIN( NB-MOD(LIHI,NB), 
     $                             N-LIHI ), NB )
                              IF(DBG) WRITE(*,*)MYROW,MYCOL, 
     $                             'dgemm 3 LCOLS=',LCOLS
                              CALL DGEMX( 'Transpose', 'No Transpose', 
     $                             NWIN, LCOLS, NWIN, ONE, DWORK( PDW ), 
     $                             NWIN, T((JLOC-1)*LLDT+ILOC), LLDT, 
     $                             ZERO, DWORK(IPW3), NWIN )
                              CALL DLACPY( 'All', NWIN, LCOLS, 
     $                             DWORK(IPW3), NWIN, 
     $                             T((JLOC-1)*LLDT+ILOC), LLDT )
                           END IF
                        END IF
                        INDXS = ICEIL(LIHI,NB)*NB + 1
                        DO 230 INDX = INDXS, N, NB 
                           CALL INFOG2L( I, INDX, DESCT, NPROW, NPCOL,
     $                          MYROW, MYCOL, ILOC, JLOC, RSRC1, CSRC1 )
                           IF( MYROW.EQ.RSRC1.AND.MYCOL.EQ.CSRC1 ) THEN
                              LCOLS = MIN( NB, N-INDX+1 )
                              IF(DBG) WRITE(*,*)MYROW,MYCOL, 
     $                          'dgemm 4 LCOLS=',LCOLS
                              CALL DGEMX( 'Transpose', 'No Transpose', 
     $                             NWIN, LCOLS, NWIN, ONE, DWORK( PDW ), 
     $                             NWIN, T((JLOC-1)*LLDT+ILOC), LLDT, 
     $                             ZERO, DWORK(IPW3), NWIN )
                              CALL DLACPY( 'All', NWIN, LCOLS, 
     $                             DWORK(IPW3), NWIN, 
     $                             T((JLOC-1)*LLDT+ILOC), LLDT )
                           END IF
 230                    CONTINUE
                     END IF
                     END IF
                  ELSE
C     
C     The NWIN-by-NWIN matrix U containing the accumulated orthogonal 
C     transformations has the following structure:
C     
C                   [ U11  U12 ]
C               U = [          ],
C                   [ U21  U22 ]
C     
C     where U21 is KS-by-KS upper triangular and U12 is
C     (NWIN-KS)-by-(NWIN-KS) lower triangular.
C     
C     Update the columns of T and Q affected by the reordering
C     
C     Compute T2*U21 + T1*U11 in workspace.
C     
                     IF( DIR.EQ.2 ) THEN
                     DO 240 INDX = 1, I-1, NB
                        CALL INFOG2L( INDX, I, DESCT, NPROW, NPCOL,
     $                       MYROW, MYCOL, ILOC, JLOC, RSRC1, CSRC1 )
                        IF( MYROW.EQ.RSRC1 .AND. MYCOL.EQ.CSRC1 ) THEN
                           JLOC1 = INDXG2L( I+NWIN-KS, NB, MYCOL, 
     $                          DESCT( CSRC_ ), NPCOL ) 
                           LROWS = MIN(NB,I-INDX)
                           CALL DLACPY( 'All', LROWS, KS, 
     $                          T((JLOC1-1)*LLDT+ILOC ), LLDT,
     $                          DWORK(IPW3), LROWS )
                           CALL DTRMX( 'Right', 'Upper', 'No transpose',
     $                          'No transpose', LROWS, KS, ONE,
     $                          DWORK( PDW+NWIN-KS ), NWIN, DWORK(IPW3),
     $                          LROWS )
                           IF(DBG) WRITE(*,*)MYROW,MYCOL, 'dgemm 5'
                           CALL DGEMX( 'No transpose', 'No transpose', 
     $                          LROWS, KS, NWIN-KS, ONE, 
     $                          T((JLOC-1)*LLDT+ILOC), LLDT,
     $                          DWORK( PDW ), NWIN, ONE, DWORK(IPW3),
     $                          LROWS )
C     
C     Compute T1*U12 + T2*U22 in workspace.
C     
                           CALL DLACPY( 'All', LROWS, NWIN-KS, 
     $                          T((JLOC-1)*LLDT+ILOC), LLDT,
     $                          DWORK( IPW3+KS*LROWS ), LROWS )
                           CALL DTRMX( 'Right', 'Lower', 'No transpose',
     $                          'No transpose', LROWS, NWIN-KS, ONE,
     $                          DWORK( PDW+NWIN*KS ), NWIN,
     $                          DWORK( IPW3+KS*LROWS ), LROWS ) 
                           IF(DBG) WRITE(*,*)MYROW,MYCOL, 'dgemm 6'
                           CALL DGEMX( 'No transpose', 'No transpose', 
     $                          LROWS, NWIN-KS, KS, ONE, 
     $                          T((JLOC1-1)*LLDT+ILOC), LLDT,
     $                          DWORK( PDW+NWIN*KS+NWIN-KS ), NWIN,
     $                          ONE, DWORK( IPW3+KS*LROWS ), 
     $                          LROWS )
C     
C     Copy workspace to T.
C     
                           CALL DLACPY( 'All', LROWS, NWIN, DWORK(IPW3), 
     $                          LROWS, T((JLOC-1)*LLDT+ILOC), LLDT )
                        END IF 
 240                 CONTINUE
                     IF( WANTQ ) THEN
C     
C     Compute Q2*U21 + Q1*U11 in workspace.
C     
                        DO 250 INDX = 1, N, NB
                           CALL INFOG2L( INDX, I, DESCQ, NPROW, NPCOL,
     $                          MYROW, MYCOL, ILOC, JLOC, RSRC1, CSRC1 )
                           IF( MYROW.EQ.RSRC1.AND.MYCOL.EQ.CSRC1 ) THEN
                              JLOC1 = INDXG2L( I+NWIN-KS, NB, 
     $                             MYCOL, DESCQ( CSRC_ ), NPCOL ) 
                              LROWS = MIN(NB,N-INDX+1)
                              CALL DLACPY( 'All', LROWS, KS, 
     $                             Q((JLOC1-1)*LLDQ+ILOC ), LLDQ,
     $                             DWORK(IPW3), LROWS )
                              CALL DTRMX( 'Right', 'Upper', 
     $                             'No transpose', 'No transpose', 
     $                             LROWS, KS, ONE, DWORK( PDW+NWIN-KS ), 
     $                             NWIN, DWORK(IPW3), LROWS )
                              IF(DBG) WRITE(*,*)MYROW,MYCOL, 'dgemm 7'
                              CALL DGEMX( 'No transpose', 
     $                             'No transpose', LROWS, KS, NWIN-KS, 
     $                             ONE, Q((JLOC-1)*LLDQ+ILOC), LLDQ,
     $                             DWORK( PDW ), NWIN, ONE, DWORK(IPW3),
     $                             LROWS )
C     
C     Compute Q1*U12 + Q2*U22 in workspace.
C     
                              CALL DLACPY( 'All', LROWS, NWIN-KS, 
     $                             Q((JLOC-1)*LLDQ+ILOC), LLDQ,
     $                             DWORK( IPW3+KS*LROWS ), LROWS)
                              CALL DTRMX( 'Right', 'Lower', 
     $                             'No transpose', 'No transpose', 
     $                             LROWS, NWIN-KS, ONE,
     $                             DWORK( PDW+NWIN*KS ), NWIN,
     $                             DWORK( IPW3+KS*LROWS ), LROWS) 
                              IF(DBG) WRITE(*,*)MYROW,MYCOL, 'dgemm 8'
                              CALL DGEMX( 'No transpose', 
     $                             'No transpose', LROWS, NWIN-KS, KS, 
     $                             ONE, Q((JLOC1-1)*LLDQ+ILOC), LLDQ,
     $                             DWORK( PDW+NWIN*KS+NWIN-KS ), NWIN,
     $                             ONE, DWORK( IPW3+KS*LROWS ), 
     $                             LROWS )
C     
C     Copy workspace to Q.
C     
                              CALL DLACPY( 'All', LROWS, NWIN, 
     $                             DWORK(IPW3), LROWS, 
     $                             Q((JLOC-1)*LLDQ+ILOC), LLDQ )
                           END IF 
 250                    CONTINUE
                     END IF
                     END IF
C     
                     IF( DIR.EQ.1 ) THEN
                     IF ( LIHI.LT.N ) THEN
C     
C     Compute U21**T*T2 + U11**T*T1 in workspace.
C     
                        IF( MOD(LIHI,NB).GT.0 ) THEN
                           INDX = LIHI + 1
                           CALL INFOG2L( I, INDX, DESCT, NPROW, NPCOL,
     $                          MYROW, MYCOL, ILOC, JLOC, RSRC1, CSRC1 )
                           IF( MYROW.EQ.RSRC1.AND.MYCOL.EQ.CSRC1 ) THEN
                              ILOC1 = INDXG2L( I+NWIN-KS, NB, MYROW, 
     $                             DESCT( RSRC_ ), NPROW )
                              LCOLS = MOD( MIN( NB-MOD(LIHI,NB), 
     $                             N-LIHI ), NB )
                              CALL DLACPY( 'All', KS, LCOLS,
     $                             T((JLOC-1)*LLDT+ILOC1), LLDT, 
     $                             DWORK(IPW3), NWIN )
                              CALL DTRMX( 'Left', 'Upper', 'Transpose',
     $                             'No transpose', KS, LCOLS, ONE,
     $                             DWORK( PDW+NWIN-KS ), NWIN, 
     $                             DWORK(IPW3), NWIN )
                              IF(DBG) WRITE(*,*)MYROW,MYCOL, 'dgemm 9'
                              CALL DGEMX( 'Transpose', 'No transpose', 
     $                             KS, LCOLS, NWIN-KS, ONE, DWORK(PDW), 
     $                             NWIN, T((JLOC-1)*LLDT+ILOC), LLDT, 
     $                             ONE, DWORK(IPW3), NWIN )
C     
C     Compute U12**T*T1 + U22**T*T2 in workspace.
C     
                              CALL DLACPY( 'All', NWIN-KS, LCOLS,
     $                             T((JLOC-1)*LLDT+ILOC), LLDT, 
     $                             DWORK( IPW3+KS ), NWIN )
                              CALL DTRMX( 'Left', 'Lower', 'Transpose',
     $                             'No transpose', NWIN-KS, LCOLS, ONE,
     $                             DWORK( PDW+NWIN*KS ), NWIN,
     $                             DWORK( IPW3+KS ), NWIN )
                              IF(DBG) WRITE(*,*)MYROW,MYCOL, 'dgemm 10'
                              CALL DGEMX( 'Transpose', 'No Transpose', 
     $                             NWIN-KS, LCOLS, KS, ONE,
     $                             DWORK( PDW+NWIN*KS+NWIN-KS ), NWIN,
     $                             T((JLOC-1)*LLDT+ILOC1), LLDT,
     $                             ONE, DWORK( IPW3+KS ), NWIN )
C     
C     Copy workspace to T.
C     
                              CALL DLACPY( 'All', NWIN, LCOLS, 
     $                             DWORK(IPW3), NWIN, 
     $                             T((JLOC-1)*LLDT+ILOC), LLDT )
                           END IF
                        END IF
                        INDXS = ICEIL(LIHI,NB)*NB + 1
                        DO 260 INDX = INDXS, N, NB 
                           CALL INFOG2L( I, INDX, DESCT, NPROW, NPCOL,
     $                          MYROW, MYCOL, ILOC, JLOC, RSRC1, CSRC1 )
                           IF( MYROW.EQ.RSRC1.AND.MYCOL.EQ.CSRC1 ) THEN
C     
C     Compute U21**T*T2 + U11**T*T1 in workspace.
C     
                              ILOC1 = INDXG2L( I+NWIN-KS, NB, 
     $                             MYROW, DESCT( RSRC_ ), NPROW ) 
                              LCOLS = MIN( NB, N-INDX+1 )
                              CALL DLACPY( 'All', KS, LCOLS,
     $                             T((JLOC-1)*LLDT+ILOC1), LLDT, 
     $                             DWORK(IPW3), NWIN )
                              CALL DTRMX( 'Left', 'Upper', 'Transpose',
     $                             'No transpose', KS, LCOLS, ONE,
     $                             DWORK( PDW+NWIN-KS ), NWIN, 
     $                             DWORK(IPW3), NWIN )
                              IF(DBG) WRITE(*,*)MYROW,MYCOL, 'dgemm 11'
                              CALL DGEMX( 'Transpose', 'No transpose', 
     $                             KS, LCOLS, NWIN-KS, ONE, DWORK(PDW), 
     $                             NWIN, T((JLOC-1)*LLDT+ILOC), LLDT, 
     $                             ONE, DWORK(IPW3), NWIN )
C     
C     Compute U12**T*T1 + U22**T*T2 in workspace.
C     
                              CALL DLACPY( 'All', NWIN-KS, LCOLS,
     $                             T((JLOC-1)*LLDT+ILOC), LLDT, 
     $                             DWORK( IPW3+KS ), NWIN )
                              CALL DTRMX( 'Left', 'Lower', 'Transpose',
     $                             'No transpose', NWIN-KS, LCOLS, ONE,
     $                             DWORK( PDW+NWIN*KS ), NWIN,
     $                             DWORK( IPW3+KS ), NWIN )
                              IF(DBG) WRITE(*,*)MYROW,MYCOL, 'dgemm 12'
                              CALL DGEMX( 'Transpose', 'No Transpose', 
     $                             NWIN-KS, LCOLS, KS, ONE,
     $                             DWORK( PDW+NWIN*KS+NWIN-KS ), NWIN,
     $                             T((JLOC-1)*LLDT+ILOC1), LLDT,
     $                             ONE, DWORK( IPW3+KS ), NWIN )
C     
C     Copy workspace to T.
C     
                              CALL DLACPY( 'All', NWIN, LCOLS, 
     $                             DWORK(IPW3), NWIN, 
     $                             T((JLOC-1)*LLDT+ILOC), LLDT )
                           END IF
 260                    CONTINUE
                     END IF
                     END IF 
                  END IF
               ELSEIF( FLOPS.NE.0 ) THEN
C     
C     Update off-diagonal blocks and Q using the pipelined elementary 
C     transformations.
C     
                  IF( DIR.EQ.2 ) THEN
                  DO 270 INDX = 1, I-1, NB
                     CALL INFOG2L( INDX, I, DESCT, NPROW, NPCOL,
     $                    MYROW, MYCOL, ILOC, JLOC, RSRC1, CSRC1 )
                     IF( MYROW.EQ.RSRC1 .AND. MYCOL.EQ.CSRC1 ) THEN
                        LROWS = MIN(NB,I-INDX)
                        CALL BDLAAPP( 1, LROWS, NWIN, NCB, 
     $                       T((JLOC-1)*LLDT+ILOC ), LLDT, NITRAF,
     $                       IWORK(IPIW), DWORK( IPW2 ), DWORK(IPW3) )
                     END IF
 270              CONTINUE
                  IF( WANTQ ) THEN
                     DO 280 INDX = 1, N, NB
                        CALL INFOG2L( INDX, I, DESCQ, NPROW, NPCOL,
     $                       MYROW, MYCOL, ILOC, JLOC, RSRC1, CSRC1 )
                        IF( MYROW.EQ.RSRC1 .AND. MYCOL.EQ.CSRC1 ) THEN
                           LROWS = MIN(NB,N-INDX+1)
                           CALL BDLAAPP( 1, LROWS, NWIN, NCB, 
     $                          Q((JLOC-1)*LLDQ+ILOC), LLDQ, NITRAF,
     $                          IWORK(IPIW), DWORK( IPW2 ), 
     $                          DWORK(IPW3) )
                        END IF
 280                 CONTINUE
                  END IF
                  END IF
                  IF( DIR.EQ.1 ) THEN
                  IF( LIHI.LT.N ) THEN
                     IF( MOD(LIHI,NB).GT.0 ) THEN
                        INDX = LIHI + 1
                        CALL INFOG2L( I, INDX, DESCT, NPROW, NPCOL,
     $                       MYROW, MYCOL, ILOC, JLOC, RSRC1, CSRC1 )
                        IF( MYROW.EQ.RSRC1 .AND. MYCOL.EQ.CSRC1 ) THEN
                           LCOLS = MOD( MIN( NB-MOD(LIHI,NB), N-LIHI ),
     $                          NB )
                           CALL BDLAAPP( 0, NWIN, LCOLS, NCB, 
     $                          T((JLOC-1)*LLDT+ILOC), LLDT, NITRAF, 
     $                          IWORK(IPIW), DWORK( IPW2 ), 
     $                          DWORK(IPW3) )
                        END IF
                     END IF
                     INDXS = ICEIL(LIHI,NB)*NB + 1
                     DO 290 INDX = INDXS, N, NB 
                        CALL INFOG2L( I, INDX, DESCT, NPROW, NPCOL,
     $                       MYROW, MYCOL, ILOC, JLOC, RSRC1, CSRC1 )
                        IF( MYROW.EQ.RSRC1 .AND. MYCOL.EQ.CSRC1 ) THEN
                           LCOLS = MIN( NB, N-INDX+1 )
                           CALL BDLAAPP( 0, NWIN, LCOLS, NCB, 
     $                          T((JLOC-1)*LLDT+ILOC), LLDT, NITRAF, 
     $                          IWORK(IPIW), DWORK( IPW2 ), 
     $                          DWORK(IPW3) )
                        END IF
 290                 CONTINUE
                  END IF
                  END IF
               END IF
C
C     If I was not involved in the updates for the current window or
C     the window was fully processed, I go here and try again for the
C     next window
C
 295           CONTINUE
C     
C     Update LIHI and LIHI depending on the number of eigenvalues really 
C     moved - for on-diagonal processes we do this update only once since
C     each on-diagonal process is only involved with one window at one
C     time. The indicies are updated in three cases:
C       1) When some reordering was really performed - indicated by
C          BUFFLEN > 0.
C       2) When no selected eigenvalues was found in the current window
C          - indicated by KS = 0.
C       3) When some selected eigenvalues was found in the current window 
C          but no one of them was moved (KS > 0 and BUFFLEN = 0)
C     False index updating is avoided by sometimes setting KS = -1. This
C     will affect processors involved in more than one window and where
C     the first one ends up with KS = 0 and for the second one is done
C     already. 
C     
               IF(DBG2) WRITE(*,*) MYROW,MYCOL,'Update LIHI and LIHI'
               IF(DBG2) WRITE(*,*) MYROW,MYCOL,'DIR,WINDOW=',DIR,WINDOW
               IF(DBG2) WRITE(*,*) MYROW,MYCOL,'RSRC,CSRC=',RSRC,CSRC
               IF(DBG2) WRITE(*,*) MYROW,MYCOL,'LILO,LIHI=',LILO,LIHI
               IF(DBG2) WRITE(*,*) MYROW,MYCOL,'I,KS=',I,KS
               IF(DBG2) WRITE(*,*) MYROW,MYCOL,'BUFFLEN=',BUFFLEN
               IF( MYROW.EQ.RSRC.AND.MYCOL.EQ.CSRC ) THEN
                  IF( DIR.EQ.2 ) THEN
                     IF( BUFFLEN.NE.0 .OR. KS.EQ.0 .OR.
     $                   ( BUFFLEN.EQ.0 .AND. KS.GT.0 ) ) 
     $                    LIHI = I + KS - 1
                     IF(DBG2) WRITE(*,*)MYROW,MYCOL, 'UPDATING LIHI=',
     $                    LIHI
                     IWORK( ILIHI+WINDOW-1 ) = LIHI
                     IF( .NOT. LIHI.GE.LILO+LSEL ) THEN
                        LILO = LILO + LSEL
                        IWORK( ILILO+WINDOW-1 ) = LILO
                        IF(DBG2) WRITE(*,*)MYROW,MYCOL, 
     $                       'UPDATING LILO=',LILO
                     END IF
                  END IF
               ELSEIF( MYROW.EQ.RSRC .AND. DIR.EQ.1 ) THEN
                  IF( BUFFLEN.NE.0 .OR. KS.EQ.0 .OR.
     $                ( BUFFLEN.EQ.0 .AND. KS.GT.0 ) )
     $                 LIHI = I + KS - 1
                  IF(DBG2) WRITE(*,*)MYROW,MYCOL, 'UPDATING LIHI=',
     $                 LIHI
                  IWORK( ILIHI+WINDOW-1 ) = LIHI
                  IF( .NOT. LIHI.GE.LILO+LSEL ) THEN
                     LILO = LILO + LSEL
                     IWORK( ILILO+WINDOW-1 ) = LILO
                     IF(DBG2) WRITE(*,*)MYROW,MYCOL, 
     $                    'UPDATING LILO=',LILO
                     IF(DBG) WRITE(*,*)MYROW,MYCOL,'ILILO+WINDOW-1=',
     $                    ILILO+WINDOW-1
                  END IF
               ELSEIF( MYCOL.EQ.CSRC .AND. DIR.EQ.2 ) THEN
                  IF( BUFFLEN.NE.0 .OR. KS.EQ.0 .OR.
     $                ( BUFFLEN.EQ.0 .AND. KS.GT.0 ) )
     $                 LIHI = I + KS - 1
                  IF(DBG2) WRITE(*,*)MYROW,MYCOL, 'UPDATING LIHI=',
     $                 LIHI
                  IWORK( ILIHI+WINDOW-1 ) = LIHI
                  IF( .NOT. LIHI.GE.LILO+LSEL ) THEN
                     LILO = LILO + LSEL
                     IWORK( ILILO+WINDOW-1 ) = LILO
                     IF(DBG2) WRITE(*,*)MYROW,MYCOL, 
     $                    'UPDATING LILO=',LILO
                     IF(DBG2) WRITE(*,*)MYROW,MYCOL,'ILILO+WINDOW-1=',
     $                    ILILO+WINDOW-1
                  END IF
               END IF
C     
 112        CONTINUE
C
C     End of direction loop for updates with respect to local reordering
C
 1111       CONTINUE
C
C     Associate each process with one of the corresponding computational 
C     windows such that the test for another round of local reordering
C     is carried out properly. Since the column updates were computed 
C     after the row updates, it is sufficient to test for changing the
C     association to the window in the corresponding process row
C
            IF(DBG2) WRITE(*,*)MYROW,MYCOL,
     $           'Associated with LILO,LIHI,LSEL='
            IF(DBG) WRITE(*,*)MYROW,MYCOL,LILO,LIHI,LSEL 
            IF(DBG) WRITE(*,*) MYROW,MYCOL,
     $           'IWORK(IRSRC:IRSRC+NMWIN2-1)='
            IF(DBG) WRITE(*,*) MYROW,MYCOL,IWORK(IRSRC:IRSRC+NMWIN2-1)
            IF(DBG) WRITE(*,*) MYROW,MYCOL,
     $           'IWORK(ICSRC:ICSRC+NMWIN2-1)='
            IF(DBG) WRITE(*,*) MYROW,MYCOL,IWORK(ICSRC:ICSRC+NMWIN2-1)
            IF(DBG) WRITE(*,*) MYROW,MYCOL,
     $           'IWORK(ILIHI:ILIHI+NMWIN2-1)='
            IF(DBG) WRITE(*,*) MYROW,MYCOL,IWORK(ILIHI:ILIHI+NMWIN2-1)
            IF(DBG) WRITE(*,*) MYROW,MYCOL,
     $           'IWORK(ILILO:ILILO+NMWIN2-1)='
            IF(DBG) WRITE(*,*) MYROW,MYCOL,IWORK(ILILO:ILILO+NMWIN2-1)
            IF(DBG) WRITE(*,*) MYROW,MYCOL,
     $           'IWORK(ILSEL:ILSEL+NMWIN2-1)='
            IF(DBG) WRITE(*,*) MYROW,MYCOL,IWORK(ILSEL:ILSEL+NMWIN2-1)
            DO 113 WINDOW = 1, NMWIN2
               IF(DBG) WRITE(*,*)MYROW,MYCOL,'WINDOW=',WINDOW
               RSRC = IWORK( IRSRC + WINDOW - 1 )
               IF( MYROW.EQ.RSRC .AND. (.NOT. LIHI.GE.LILO+LSEL ) ) THEN
                  LILO = IWORK( ILILO + WINDOW - 1 )
                  IF(DBG) WRITE(*,*)MYROW,MYCOL,'ILILO+WINDOW-1=',
     $                 ILILO+WINDOW-1
                  LIHI = IWORK( ILIHI + WINDOW - 1 )
                  LSEL = IWORK( ILSEL + WINDOW - 1 )
                  IF(DBG2) WRITE(*,*)MYROW,MYCOL,
     $           'Changed association to LILO,LIHI,LSEL='
                  IF(DBG2) WRITE(*,*)MYROW,MYCOL,LILO,LIHI,LSEL
               END IF
 113        CONTINUE
C     
C     End While ( LIHI >= LILO + LSEL )
            ROUND = ROUND + 1
            IF( FIRST ) FIRST = .FALSE.
            GO TO 130
         END IF
C
C     All processors excluded from the local reordering go here
C
 114     CONTINUE
C     
C     Barrier to collect the processes before proceeding
C     
         IF(DBG3) WRITE(*,*)MYROW,MYCOL,'barrier 1 in'
         CALL BLACS_BARRIER( ICTXT, 'All' )
         IF(DBG3) WRITE(*,*)MYROW,MYCOL,'barrier 1 out'
C     
C     Compute global maximum of IERR so that we know if some process
C     experienced a failure in the reordering
C     
         IF(DBG) WRITE(*,*) MYROW,MYCOL,'IERR=',IERR
         MYIERR = IERR
         IF(DBG) WRITE(*,*) MYROW,MYCOL,'MYIERR=',MYIERR
         IF( NPROCS.GT.1 )
     $        CALL IGAMX2D( ICTXT, 'All', TOP, 1, 1, IERR, 1, -1, 
     $        -1, -1, -1, -1 )
         IF(DBG) WRITE(*,*) MYROW,MYCOL,'IERR=',IERR
C     
         IF( IERR.EQ.1 .OR. IERR.EQ.2 ) THEN
            IF(DBG) WRITE(*,*) MYROW,MYCOL,
     $           'an error has occured - exiting'
C     
C     When calling BDTREXC, the block at position I+KKS-1 failed to swap.
C     
            IF( IERR.EQ.1 .OR. IERR.EQ.2 ) THEN
               INFO = MAX(1,I+KKS-1)
            ELSE 
               INFO = 1
            END IF
            IF(DBG) WRITE(*,*) MYROW,MYCOL,
     $           'an error has occured - resulting INFO=',INFO
            IF( WANTS )
     $           S = ZERO
            IF( WANTSP )
     $           SEP = ZERO
            GO TO 300
         END IF
C
         TIME1 = TIME1 + MPI_WTIME() - TIME
         TIME = MPI_WTIME()
C
         IF(DBG) WRITE(*,*) MYROW,MYCOL,
     $        'going for cross border computations'
C     
C     Now, for each compuational window, move the selected eigenvalues
C     across the process border. Do this by forming the processors into
C     groups of four working together to bring the window over the 
C     border. The processes are numbered as follows
C
C             1 | 2
C             -----
C             3 | 4
C     
C     where '|' and '-' denotes the process (and block) borders.
C     This implies that the cluster to be reordered over the border
C     is held by process 4, process 1 will receive the cluster after
C     the reordering, process 3 holds the local (2,1)th element of a 
C     2-by-2 diagonal block located on the block border and process 2
C     holds the closest off-diagonal part of the windo that is affected 
C     by the cross-border reordering.
C
C     The active window is now ( I : LIHI[4], I : LIHI[4] ), where 
C     I = MAX( ILO, LIHI - 2*MOD(LIHI,NB) ). If this active window is to 
C     large compared to the value of PARA( 6 ), it will be truncated in 
C     both ends such that a maximum of PARA( 6 ) eigenvalues is reordered 
C     across the border this time.
C     
C     The active window will be collected and built in workspace at 
C     process 1 and 4, which both compute the reordering and return the 
C     updated parts to the corresponding processes 2-3. Next, the 
C     accumulated transformations are broadcasted for updates in the 
C     block rows and column that corresponds to the process rows and 
C     columns where process 1 and 4 reside.  
C
C     The off-diagonal blocks are updated by the processes receiving
C     from the broadcasts of the orthogonal transformations. Since the
C     active window is split over the process borders, the updates
C     of T and Q requires that stripes of block rows of columns are
C     exchanged between neighboring processes in the corresponding
C     process rows and columns.
C
C     First, form each group of processors involved in the crossborder
C     reordering. Do this in two (or three) phases: 
C     1) Reorder each odd window over the border.
C     2) Reorder each even window over the border.
C     3) Reorder the last odd window over the border, if it was not
C        processed in the first phase.
C
C     When reordering the odd windows over the border, we must make sure
C     that no process row or column is involved in both the first and 
C     the last window at the same time. This happens when the total 
C     number of windows is odd, greater than one and equal to the 
C     minumum process mesh dimension. Therefore the last odd window may 
C     be reordered over the border at last.
C         
         LASTWAIT = NMWIN2.GT.1 .AND. MOD(NMWIN2,2).EQ.1 .AND. 
     $        NMWIN2.EQ.MIN(NPROW,NPCOL)
         IF(DBG) WRITE(*,*) MYROW,MYCOL,'LASTWAIT=',LASTWAIT
C
         LAST = 0
 308     CONTINUE
         IF( LASTWAIT ) THEN
            IF( LAST.EQ.0 ) THEN
               WIN0S = 1
               WIN0E = 2
               WINE = NMWIN2 - 1
            ELSE
               WIN0S = NMWIN2
               WIN0E = NMWIN2
               WINE = NMWIN2
            END IF
         ELSE
            WIN0S = 1
            WIN0E = 2
            WINE = NMWIN2
         END IF
         IF(DBG) WRITE(*,*) MYROW,MYCOL,'LAST=',LAST
         DO 310 WINDOW0 = WIN0S, WIN0E
            DO 320 WINDOW = WINDOW0, WINE, 2
C
C     Define the process holding the down-right part of the window
C
               RSRC4 = IWORK(IRSRC+WINDOW-1)
               CSRC4 = IWORK(ICSRC+WINDOW-1)
C
C     Define the other processes in the group of four
C
               RSRC3 = RSRC4
               CSRC3 = MOD( CSRC4 - 1 + NPCOL, NPCOL )
               RSRC2 = MOD( RSRC4 - 1 + NPROW, NPROW )
               CSRC2 = CSRC4
               RSRC1 = RSRC2
               CSRC1 = CSRC3
               IF( ( MYROW.EQ.RSRC1 .AND. MYCOL.EQ.CSRC1 ) .OR.
     $             ( MYROW.EQ.RSRC2 .AND. MYCOL.EQ.CSRC2 ) .OR.
     $             ( MYROW.EQ.RSRC3 .AND. MYCOL.EQ.CSRC3 ) .OR.
     $             ( MYROW.EQ.RSRC4 .AND. MYCOL.EQ.CSRC4 ) ) THEN
C     
C     Compute the correct active window - for reordering into a block
C     that has not been active at all before, we try to reorder as
C     many of our eigenvalues over the border as possible without
C     knowing of the situation on the other side - this may cause
C     very few eigenvalues to be reordered over the border this time
C     (perhaps not any) but this should be an initial problem. 
C     Anyway, the bottom-right position of the block will be at
C     position LIHIC.
C
                  IF( MYROW.EQ.RSRC4 .AND. MYCOL.EQ.CSRC4 ) THEN
                     LIHI4 = ( IWORK( ILILO + WINDOW - 1 ) + 
     $                    IWORK( ILIHI + WINDOW - 1 ) ) / 2 
                     IF(DBG) WRITE(*,*)MYROW,MYCOL, 'LIHI4=',LIHI4
                     LIHIC = MIN(LIHI4,(ICEIL(LIHI4,NB)-1)*NB+WNEICR)
                     IF(DBG) WRITE(*,*)MYROW,MYCOL, 'LIHIC=',LIHIC
C
C     Fix LIHIC to avoid that bottom of window cuts 2-by-2 block and
C     make sure all processors in the group knows about the correct
C     value
C
                     IF( (.NOT. LIHIC.LE.NB) .AND. LIHIC.LT.N ) THEN
                        ILOC = INDXG2L( LIHIC+1, NB, MYROW, 
     $                       DESCT( RSRC_ ), NPROW )
                        JLOC = INDXG2L( LIHIC, NB, MYCOL, 
     $                       DESCT( CSRC_ ), NPCOL )
                        IF( T( (JLOC-1)*LLDT+ILOC ).NE.ZERO ) THEN
                           IF( MOD( LIHIC, NB ).EQ.1 .OR.
     $                          ( MOD( LIHIC, NB ).EQ.2 .AND.
     $                          .NOT. INT2LG(SELECT(LIHIC-2) )) ) 
     $                          THEN
                              LIHIC = LIHIC + 1
                           ELSE
                              LIHIC = LIHIC - 1
                           END IF
                        END IF
                     END IF
                     IF(DBG) WRITE(*,*)MYROW,MYCOL, 'LIHIC=',LIHIC
                     IF( RSRC4.NE.RSRC1 .OR. CSRC4.NE.CSRC1 ) 
     $                    CALL IGESD2D( ICTXT, 1, 1, LIHIC, 1, RSRC1,
     $                    CSRC1 )
                     IF( RSRC4.NE.RSRC2 .OR. CSRC4.NE.CSRC2 ) 
     $                    CALL IGESD2D( ICTXT, 1, 1, LIHIC, 1, RSRC2,
     $                    CSRC2 )
                     IF( RSRC4.NE.RSRC3 .OR. CSRC4.NE.CSRC3 ) 
     $                    CALL IGESD2D( ICTXT, 1, 1, LIHIC, 1, RSRC3,
     $                    CSRC3 ) 
                  END IF
                  IF( MYROW.EQ.RSRC1 .AND. MYCOL.EQ.CSRC1 ) THEN
                     IF( RSRC4.NE.RSRC1 .OR. CSRC4.NE.CSRC1 ) 
     $                    CALL IGERV2D( ICTXT, 1, 1, LIHIC, 1, RSRC4,
     $                    CSRC4 )
                  END IF
                  IF( MYROW.EQ.RSRC2 .AND. MYCOL.EQ.CSRC2 ) THEN
                     IF( RSRC4.NE.RSRC2 .OR. CSRC4.NE.CSRC2 ) 
     $                    CALL IGERV2D( ICTXT, 1, 1, LIHIC, 1, RSRC4,
     $                    CSRC4 )
                  END IF
                  IF( MYROW.EQ.RSRC3 .AND. MYCOL.EQ.CSRC3 ) THEN
                     IF( RSRC4.NE.RSRC3 .OR. CSRC4.NE.CSRC3 ) 
     $                    CALL IGERV2D( ICTXT, 1, 1, LIHIC, 1, RSRC4,
     $                    CSRC4 )
                  END IF
C
C     Avoid going over the border with the first window if it resides
C     in the block where the last global position T(ILO,ILO) is or
C     ILO has been updated to point to a position right of 
C     T(LIHIC,LIHIC)
C
                  SKIP1CR = WINDOW.EQ.1 .AND. 
     $                 ICEIL(LIHIC,NB).LE.ICEIL(ILO,NB) 
C     
C     Decide I, where to put top of window, such that top of window 
C     does not cut 2-by-2 block. Make sure that we do not end up
C     in a situation where a 2-by-2 block splitted on the border
C     is left in its original place - this can cause infinite loops.
C     Remedy: make sure that the part of the window that resides left
C     to the border is at least of dimension two (2) in case we have
C     2-by-2 blocks in top of the cross border window.
C
C     Also make sure all processors in the group knows about the correct 
C     value of I. When skipping the crossborder reordering, just set
C     I = LIHIC.
C
                  IF( .NOT. SKIP1CR ) THEN
                     IF( MYROW.EQ.RSRC1 .AND. MYCOL.EQ.CSRC1 ) THEN
                        IF( WINDOW.EQ.1 ) THEN
                           LIHI1 = ILO 
                        ELSE
                           LIHI1 = IWORK( ILIHI + WINDOW - 2 )
                        END IF
                        IF(DBG) WRITE(*,*) MYROW,MYCOL,'LIHI1=',LIHI1
                        I = MAX( LIHI1, 
     $                       MIN( LIHIC-2*MOD(LIHIC,NB) + 1, 
     $                       (ICEIL(LIHIC,NB)-1)*NB - 1  ) )
                        IF(DBG) WRITE(*,*)MYROW,MYCOL, 'I=',I
                        ILOC = INDXG2L( I, NB, MYROW, DESCT( RSRC_ ), 
     $                       NPROW )
                        JLOC = INDXG2L( I-1, NB, MYCOL, DESCT( CSRC_ ), 
     $                       NPCOL )
                        IF( T( (JLOC-1)*LLDT+ILOC ).NE.ZERO ) THEN
                           I = I - 1
                        END IF
                        IF(DBG) WRITE(*,*)MYROW,MYCOL, 'I=',I
                        IF( RSRC1.NE.RSRC4 .OR. CSRC1.NE.CSRC4 ) 
     $                       CALL IGESD2D( ICTXT, 1, 1, I, 1, RSRC4, 
     $                       CSRC4 )
                        IF( RSRC1.NE.RSRC2 .OR. CSRC1.NE.CSRC2 ) 
     $                       CALL IGESD2D( ICTXT, 1, 1, I, 1, RSRC2,
     $                       CSRC2 )
                        IF( RSRC1.NE.RSRC3 .OR. CSRC1.NE.CSRC3 ) 
     $                       CALL IGESD2D( ICTXT, 1, 1, I, 1, RSRC3,
     $                       CSRC3 ) 
                     END IF
                     IF( MYROW.EQ.RSRC2 .AND. MYCOL.EQ.CSRC2 ) THEN
                        IF( RSRC1.NE.RSRC2 .OR. CSRC1.NE.CSRC2 ) 
     $                       CALL IGERV2D( ICTXT, 1, 1, I, 1, RSRC1,
     $                       CSRC1 )
                     END IF
                     IF( MYROW.EQ.RSRC3 .AND. MYCOL.EQ.CSRC3 ) THEN
                        IF( RSRC1.NE.RSRC3 .OR. CSRC1.NE.CSRC3 ) 
     $                       CALL IGERV2D( ICTXT, 1, 1, I, 1, RSRC1,
     $                       CSRC1 )
                     END IF
                     IF( MYROW.EQ.RSRC4 .AND. MYCOL.EQ.CSRC4 ) THEN
                        IF( RSRC1.NE.RSRC4 .OR. CSRC1.NE.CSRC4 ) 
     $                       CALL IGERV2D( ICTXT, 1, 1, I, 1, RSRC1,
     $                       CSRC1 )
                     END IF            
                  ELSE
                     I = LIHIC
                  END IF
C
C     Finalize computation of window size: active window is now
C     (I:LIHIC,I:LIHIC)
C
                  IF(DBG) WRITE(*,*)MYROW,MYCOL, 'WINDOW=',WINDOW
                  IF(DBG) WRITE(*,*)MYROW,MYCOL, 'I,LIHIC=',I,LIHIC
                  NWIN = LIHIC - I + 1
                  KS = 0
C
C     Skip rest of this part if appropriate
C
                  IF( SKIP1CR ) GO TO 360
C
C     Divide workspace - put active window in DWORK(IPW2:IPW2+NWIN**2-1)
C     and orthogonal transformations in DWORK(IPW3:...)
C
                  CALL DLASET( 'All', NWIN, NWIN, ZERO, ZERO, 
     $                 DWORK( IPW2 ), NWIN )  
C     
                  PITRAF = IPIW
                  IPW3 = IPW2 + NWIN*NWIN
                  PDTRAF = IPW3
C     
C     Exchange the current view of SELECT for the active window between 
C     process 1 and 4 to make sure that exactly the same job is performed
C     for both processes.
C
                  IF(DBG) WRITE(*,*) MYROW,MYCOL,'Exchange SELECT'
                  IF( RSRC1.NE.RSRC4 .OR. CSRC1.NE.CSRC4 ) THEN
                     ILEN4 = MOD(LIHIC,NB)
                     SELI4 = ICEIL(I,NB)*NB+1
                     ILEN1 = NWIN - ILEN4 
                     IF(DBG) WRITE(*,*) MYROW,MYCOL,'ILEN1,ILEN4=',
     $                    ILEN1,ILEN4
                     IF( MYROW.EQ.RSRC1 .AND. MYCOL.EQ.CSRC1 ) THEN
                        IF(DBG) WRITE(*,*) MYROW,MYCOL,'send'
                        CALL IGESD2D( ICTXT, ILEN1, 1, SELECT(I), 
     $                       ILEN1, RSRC4, CSRC4 )
                        IF(DBG) WRITE(*,*) MYROW,MYCOL,'revc'
                        CALL IGERV2D( ICTXT, ILEN4, 1, SELECT(SELI4), 
     $                       ILEN4, RSRC4, CSRC4 )
                     END IF
                     IF( MYROW.EQ.RSRC4 .AND. MYCOL.EQ.CSRC4 ) THEN
                        IF(DBG) WRITE(*,*) MYROW,MYCOL,'send'
                        CALL IGESD2D( ICTXT, ILEN4, 1, SELECT(SELI4), 
     $                       ILEN4, RSRC1, CSRC1 )
                        IF(DBG) WRITE(*,*) MYROW,MYCOL,'revc'
                        CALL IGERV2D( ICTXT, ILEN1, 1, SELECT(I), 
     $                       ILEN1, RSRC1, CSRC1 )
                     END IF
                  END IF
C     
C     Form the active window by a series of point-to-point sends
C     and receives.  
C
                  IF(DBG) WRITE(*,*) MYROW,MYCOL,'Form active window'
                  DIM1 = NB - MOD(I-1,NB) 
                  DIM4 = NWIN - DIM1
                  IF(DBG) WRITE(*,*)MYROW,MYCOL, 'DIM1,DIM4=',DIM1,DIM4
                  IF( MYROW.EQ.RSRC1 .AND. MYCOL.EQ.CSRC1 ) THEN
                     ILOC = INDXG2L( I, NB, MYROW, DESCT( RSRC_ ), 
     $                    NPROW )
                     JLOC = INDXG2L( I, NB, MYCOL, DESCT( CSRC_ ), 
     $                    NPCOL ) 
                     CALL DLACPY( 'All', DIM1, DIM1, 
     $                    T((JLOC-1)*LLDT+ILOC), LLDT, DWORK(IPW2),
     $                    NWIN )
                     IF( RSRC1.NE.RSRC4 .OR. CSRC1.NE.CSRC4 ) THEN
                        IF(DBG) WRITE(*,*) MYROW,MYCOL,'send 1->4'
                        CALL DGESD2D( ICTXT, DIM1, DIM1, 
     $                       DWORK(IPW2), NWIN, RSRC4, CSRC4 )
                        IF(DBG) WRITE(*,*) MYROW,MYCOL,'recv 4<-1'
                        CALL DGERV2D( ICTXT, DIM4, DIM4, 
     $                       DWORK(IPW2+DIM1*NWIN+DIM1), NWIN, RSRC4, 
     $                       CSRC4 )
                     END IF
                  END IF
                  IF( MYROW.EQ.RSRC4 .AND. MYCOL.EQ.CSRC4 ) THEN
                     ILOC = INDXG2L( I+DIM1, NB, MYROW, DESCT( RSRC_ ), 
     $                    NPROW )
                     JLOC = INDXG2L( I+DIM1, NB, MYCOL, DESCT( CSRC_ ), 
     $                    NPCOL ) 
                     CALL DLACPY( 'All', DIM4, DIM4, 
     $                    T((JLOC-1)*LLDT+ILOC), LLDT, 
     $                    DWORK(IPW2+DIM1*NWIN+DIM1), NWIN )
                     IF( RSRC4.NE.RSRC1 .OR. CSRC4.NE.CSRC1 ) THEN
                        IF(DBG) WRITE(*,*) MYROW,MYCOL,'send 4->1'
                        CALL DGESD2D( ICTXT, DIM4, DIM4, 
     $                       DWORK(IPW2+DIM1*NWIN+DIM1), NWIN, RSRC1, 
     $                       CSRC1 ) 
                        IF(DBG) WRITE(*,*) MYROW,MYCOL,'recv 1<-4'
                        CALL DGERV2D( ICTXT, DIM1, DIM1, 
     $                       DWORK(IPW2), NWIN, RSRC1, CSRC1 )
                     END IF
                  END IF
                  IF( MYROW.EQ.RSRC2 .AND. MYCOL.EQ.CSRC2 ) THEN
                     ILOC = INDXG2L( I, NB, MYROW, DESCT( RSRC_ ), 
     $                    NPROW )
                     JLOC = INDXG2L( I+DIM1, NB, MYCOL, DESCT( CSRC_ ), 
     $                    NPCOL )
                     IF(DBG) WRITE(*,*) MYROW,MYCOL,' IPW2,DIM1,NWIN=',
     $                    IPW2,DIM1,NWIN     
                     CALL DLACPY( 'All', DIM1, DIM4, 
     $                    T((JLOC-1)*LLDT+ILOC), LLDT, 
     $                    DWORK(IPW2+DIM1*NWIN), NWIN )
                     IF( RSRC2.NE.RSRC1 .OR. CSRC2.NE.CSRC1 ) THEN
                        IF(DBG) WRITE(*,*) MYROW,MYCOL,'send 2->1'
                        CALL DGESD2D( ICTXT, DIM1, DIM4, 
     $                       DWORK(IPW2+DIM1*NWIN), NWIN, RSRC1, CSRC1 )
                     END IF
                  END IF
                  IF( MYROW.EQ.RSRC2 .AND. MYCOL.EQ.CSRC2 ) THEN 
                     IF( RSRC2.NE.RSRC4 .OR. CSRC2.NE.CSRC4 ) THEN
                        IF(DBG) WRITE(*,*) MYROW,MYCOL,'send 2->4'
                        CALL DGESD2D( ICTXT, DIM1, DIM4, 
     $                       DWORK(IPW2+DIM1*NWIN), NWIN, RSRC4, CSRC4 )
                     END IF
                  END IF
                  IF( MYROW.EQ.RSRC3 .AND. MYCOL.EQ.CSRC3 ) THEN
                     ILOC = INDXG2L( I+DIM1, NB, MYROW, DESCT( RSRC_ ), 
     $                    NPROW )
                     JLOC = INDXG2L( I+DIM1-1, NB, MYCOL, 
     $                    DESCT( CSRC_ ), NPCOL )
                     CALL DLACPY( 'All', 1, 1, 
     $                    T((JLOC-1)*LLDT+ILOC), LLDT, 
     $                    DWORK(IPW2+(DIM1-1)*NWIN+DIM1), NWIN )
                     IF( RSRC3.NE.RSRC1 .OR. CSRC3.NE.CSRC1 ) THEN 
                        IF(DBG) WRITE(*,*) MYROW,MYCOL,'send 3->1'
                        CALL DGESD2D( ICTXT, 1, 1, 
     $                       DWORK(IPW2+(DIM1-1)*NWIN+DIM1), NWIN, 
     $                       RSRC1, CSRC1 )
                     END IF
                  END IF
                  IF( MYROW.EQ.RSRC3 .AND. MYCOL.EQ.CSRC3 ) THEN
                     IF( RSRC3.NE.RSRC4 .OR. CSRC3.NE.CSRC4 ) THEN 
                        IF(DBG) WRITE(*,*) MYROW,MYCOL,'send 3->4'
                        CALL DGESD2D( ICTXT, 1, 1, 
     $                       DWORK(IPW2+(DIM1-1)*NWIN+DIM1), NWIN, 
     $                       RSRC4, CSRC4 )
                     END IF
                  END IF
                  IF( MYROW.EQ.RSRC1 .AND. MYCOL.EQ.CSRC1 ) THEN
                     IF( RSRC1.NE.RSRC2 .OR. CSRC1.NE.CSRC2 ) THEN
                        IF(DBG) WRITE(*,*) MYROW,MYCOL,'recv 1<-2'
                        CALL DGERV2D( ICTXT, DIM1, DIM4,
     $                       DWORK(IPW2+DIM1*NWIN), NWIN, RSRC2, 
     $                       CSRC2 )
                     END IF
                     IF( RSRC1.NE.RSRC3 .OR. CSRC1.NE.CSRC3 ) THEN
                        IF(DBG) WRITE(*,*) MYROW,MYCOL,'recv 1<-3'
                        CALL DGERV2D( ICTXT, 1, 1,
     $                       DWORK(IPW2+(DIM1-1)*NWIN+DIM1), NWIN, 
     $                       RSRC3, CSRC3 )
                     END IF
                  END IF
                  IF( MYROW.EQ.RSRC4 .AND. MYCOL.EQ.CSRC4 ) THEN
                     IF( RSRC4.NE.RSRC2 .OR. CSRC4.NE.CSRC2 ) THEN
                        IF(DBG) WRITE(*,*) MYROW,MYCOL,'recv 4<-2'
                        CALL DGERV2D( ICTXT, DIM1, DIM4,
     $                       DWORK(IPW2+DIM1*NWIN), NWIN, RSRC2, 
     $                       CSRC2 )
                     END IF
                     IF( RSRC4.NE.RSRC3 .OR. CSRC4.NE.CSRC3 ) THEN
                        IF(DBG) WRITE(*,*) MYROW,MYCOL,'recv 4<-3'
                        CALL DGERV2D( ICTXT, 1, 1,
     $                       DWORK(IPW2+(DIM1-1)*NWIN+DIM1), NWIN, 
     $                       RSRC3, CSRC3 )
                     END IF
                  END IF
C
C     Compute the reordering (just as in the total local case) and
C     accumulate the transformations (ONLY ON-DIAGONAL PROCESSES)
C     
                  IF( ( MYROW.EQ.RSRC1 .AND. MYCOL.EQ.CSRC1 ) .OR.
     $                ( MYROW.EQ.RSRC4 .AND. MYCOL.EQ.CSRC4 ) ) THEN 
                     IF(DBG) WRITE(*,*)MYROW,MYCOL, 
     $                    'computing crossborder reordering'
                     PAIR = .FALSE.
                     DO 330 K = I, LIHIC
                        IF( PAIR ) THEN
                           PAIR = .FALSE.
                        ELSE
                           SWAP = INT2LG(SELECT( K ))
                           IF( K.LT.LIHIC ) THEN
                              ELEM = DWORK(IPW2+(K-I)*NWIN+K-I+1)
                              IF( ELEM.NE.ZERO )
     $                             PAIR = .TRUE.
                           END IF
                           IF( SWAP ) THEN
                              KS = KS + 1
C     
C     Swap the K-th block to position I+KS-1.
C     
                              IERR = 0
                              KK  = K - I + 1
                              KKS = KS
                              IF( KK.NE.KS ) THEN
                                 NITRAF = LIWORK - PITRAF + 1
                                 NDTRAF = LDWORK - PDTRAF + 1
                                 IF( DBG4 ) WRITE(*,*) 'BDTREXC 2'
                                 CALL BDTREXC( NWIN, DWORK(IPW2), NWIN, 
     $                                KK, KKS, NITRAF, IWORK( PITRAF ), 
     $                                NDTRAF, DWORK( PDTRAF ), 
     $                                DWORK(IPW1), IERR )
                                 PITRAF = PITRAF + NITRAF
                                 PDTRAF = PDTRAF + NDTRAF
C     
C     Update array SELECT.
C     
                                 IF ( PAIR ) THEN
                                    DO 340 J = I+KK-1, I+KKS, -1
                                       SELECT(J+1) = SELECT(J-1)
 340                                CONTINUE
                                    SELECT(I+KKS-1) = 1
                                    SELECT(I+KKS) = 1
                                 ELSE
                                    DO 350 J = I+KK-1, I+KKS, -1
                                       SELECT(J) = SELECT(J-1)
 350                                CONTINUE
                                    SELECT(I+KKS-1) = 1
                                 END IF
C     
                                 IF ( IERR.EQ.1 .OR. IERR.EQ.2 ) THEN
C     
                                    IF ( IERR.EQ.2 ) THEN
                                       SELECT( I+KKS-3 ) = 1
                                       SELECT( I+KKS-1 ) = 0
                                       KKS = KKS + 1
                                    END IF
C     
                                    GO TO 360
                                 END IF
                                 KS = KKS     
                              END IF
                              IF( PAIR )
     $                             KS = KS + 1
                           END IF
                        END IF
 330                 CONTINUE
                  END IF
 360              CONTINUE
C     
C     Save information about the reordering 
C     
                  IF( ( MYROW.EQ.RSRC1 .AND. MYCOL.EQ.CSRC1 ) .OR.
     $                ( MYROW.EQ.RSRC4 .AND. MYCOL.EQ.CSRC4 ) ) THEN 
                     IBUFF( 1 ) = I
                     IBUFF( 2 ) = NWIN
                     IBUFF( 3 ) = PITRAF
                     IBUFF( 4 ) = KS
                     IBUFF( 5 ) = PDTRAF
                     IBUFF( 6 ) = NDTRAF 
                     ILEN = PITRAF - IPIW + 1
                     DLEN = PDTRAF - IPW3 + 1
                     IBUFF( 7 ) = ILEN
                     IBUFF( 8 ) = DLEN
C     
C     Put reordered data back into global matrix if a reordering
C     took place
C    
                     IF( .NOT. SKIP1CR ) THEN 
                        IF( MYROW.EQ.RSRC1 .AND. MYCOL.EQ.CSRC1 ) THEN
                           ILOC = INDXG2L( I, NB, MYROW, DESCT( RSRC_ ), 
     $                          NPROW )
                           JLOC = INDXG2L( I, NB, MYCOL, DESCT( CSRC_ ), 
     $                          NPCOL ) 
                           CALL DLACPY( 'All', DIM1, DIM1, DWORK(IPW2),
     $                          NWIN, T((JLOC-1)*LLDT+ILOC), LLDT )
                        END IF
                        IF( MYROW.EQ.RSRC4 .AND. MYCOL.EQ.CSRC4 ) THEN
                           ILOC = INDXG2L( I+DIM1, NB, MYROW, 
     $                          DESCT( RSRC_ ), NPROW )
                           JLOC = INDXG2L( I+DIM1, NB, MYCOL, 
     $                          DESCT( CSRC_ ), NPCOL ) 
                           CALL DLACPY( 'All', DIM4, DIM4, 
     $                          DWORK(IPW2+DIM1*NWIN+DIM1), NWIN,
     $                          T((JLOC-1)*LLDT+ILOC), LLDT )
                        END IF
                     END IF
                  END IF
C
C     Break if appropriate - IBUFF(3:8) may now contain nonsens, but
C     that's no problem. The processors outside the cross border group
C     only needs to know about I and NWIN to get a correct value of
C     SKIP1CR (see below) and to skip the cross border updates if
C     necessary.
C     
                  IF( WINDOW.EQ.1 .AND. SKIP1CR ) GO TO 325
C     
C     Return reordered data to process 2 and 3
C     
                  IF(DBG) WRITE(*,*)MYROW,MYCOL, 'return reordered data'
                  IF( MYROW.EQ.RSRC1 .AND. MYCOL.EQ.CSRC1 ) THEN
                     IF( RSRC1.NE.RSRC3 .OR. CSRC1.NE.CSRC3 ) THEN
                        IF(DBG) WRITE(*,*) MYROW,MYCOL,'send 1->3'
                        CALL DGESD2D( ICTXT, 1, 1, 
     $                       DWORK( IPW2+(DIM1-1)*NWIN+DIM1 ), NWIN,
     $                       RSRC3, CSRC3 )
                        IF(DBG) WRITE(*,*) MYROW,MYCOL,'send 1->3 done'
                     END IF
                  END IF
                  IF( MYROW.EQ.RSRC4 .AND. MYCOL.EQ.CSRC4 ) THEN
                     IF( RSRC4.NE.RSRC2 .OR. CSRC4.NE.CSRC2 ) THEN
                        IF(DBG) WRITE(*,*) MYROW,MYCOL,'send 4->2'
                        CALL DGESD2D( ICTXT, DIM1, DIM4, 
     $                       DWORK( IPW2+DIM1*NWIN), NWIN, RSRC2, 
     $                       CSRC2 )
                        IF(DBG) WRITE(*,*) MYROW,MYCOL,'send 4->2 done' 
                     END IF
                  END IF
                  IF( MYROW.EQ.RSRC2 .AND. MYCOL.EQ.CSRC2 ) THEN
                     ILOC = INDXG2L( I, NB, MYROW, DESCT( RSRC_ ), 
     $                    NPROW )
                     JLOC = INDXG2L( I+DIM1, NB, MYCOL, 
     $                    DESCT( CSRC_ ), NPCOL )
                     IF( RSRC2.NE.RSRC4 .OR. CSRC2.NE.CSRC4 ) THEN 
                        IF(DBG) WRITE(*,*) MYROW,MYCOL,'recv 2<-4'
                        CALL DGERV2D( ICTXT, DIM1, DIM4, 
     $                       DWORK(IPW2+DIM1*NWIN), NWIN, RSRC4, CSRC4 )
                        IF(DBG) WRITE(*,*) MYROW,MYCOL,'recv 2<-4 done'
                     END IF
                     CALL DLACPY( 'All', DIM1, DIM4, 
     $                    DWORK( IPW2+DIM1*NWIN ), NWIN,
     $                    T((JLOC-1)*LLDT+ILOC), LLDT )
                  END IF
                  IF( MYROW.EQ.RSRC3 .AND. MYCOL.EQ.CSRC3 ) THEN
                     ILOC = INDXG2L( I+DIM1, NB, MYROW, 
     $                    DESCT( RSRC_ ), NPROW )
                     JLOC = INDXG2L( I+DIM1-1, NB, MYCOL, 
     $                    DESCT( CSRC_ ), NPCOL )
                     IF( RSRC3.NE.RSRC1 .OR. CSRC3.NE.CSRC1 ) THEN
                        IF(DBG) WRITE(*,*) MYROW,MYCOL,'recv 3<-1'
                        CALL DGERV2D( ICTXT, 1, 1, 
     $                       DWORK( IPW2+(DIM1-1)*NWIN+DIM1 ), NWIN,
     $                       RSRC1, CSRC1 )
                        IF(DBG) WRITE(*,*) MYROW,MYCOL,'recv 3<-1 done'
                     END IF
                     T((JLOC-1)*LLDT+ILOC) = 
     $                    DWORK( IPW2+(DIM1-1)*NWIN+DIM1 )
                  END IF
               END IF
C
 325           CONTINUE
C
 320        CONTINUE
C
C     For the crossborder updates, we use the same directions as in
C     the local reordering case above
C
            DO 2222 DIR = 1, 2
               IF(DBG) WRITE(*,*) MYROW,MYCOL,'DIR=',DIR
C     
C     Broadcast informations about the reordering
C     
            DO 321 WINDOW = WINDOW0, WINE, 2
               IF(DBG) WRITE(*,*)MYROW,MYCOL, 
     $              'broadcasting info about crb reordering'
               RSRC4 = IWORK(IRSRC+WINDOW-1)
               CSRC4 = IWORK(ICSRC+WINDOW-1)
               RSRC1 = MOD( RSRC4 - 1 + NPROW, NPROW )
               CSRC1 = MOD( CSRC4 - 1 + NPCOL, NPCOL )
               IF(DBG) WRITE(*,*) MYROW,MYCOL,
     $              'RSRC1,CSRC1,RSRC4,CSRC4=',
     $              RSRC1,CSRC1,RSRC4,CSRC4 
               IF( MYROW.EQ.RSRC1 .AND. MYCOL.EQ.CSRC1 ) THEN
                  IF(DBG) WRITE(*,*) MYROW,MYCOL,'bcast send row 1'
                  IF(DBG) WRITE(*,*) MYROW,MYCOL,'IBUFF2',IBUFF(1:8)
                  IF( NPCOL.GT.1 .AND. DIR.EQ.1 )
     $                 CALL IGEBS2D( ICTXT, 'Row', TOP, 8, 1, 
     $                 IBUFF, 8 )
                  IF(DBG) WRITE(*,*) MYROW,MYCOL,'bcast send col 1'
                  IF( NPROW.GT.1 .AND. DIR.EQ.2 )
     $                 CALL IGEBS2D( ICTXT, 'Col', TOP, 8, 1, 
     $                 IBUFF, 8 )
                  SKIP1CR = WINDOW.EQ.1 .AND. 
     $                 ICEIL(LIHIC,NB).LE.ICEIL(ILO,NB)
               ELSEIF( MYROW.EQ.RSRC1 .OR. MYCOL.EQ.CSRC1 ) THEN
                  IF( NPCOL.GT.1 .AND. DIR.EQ.1 .AND. MYROW.EQ.RSRC1 ) 
     $                 THEN
                     IF(DBG) WRITE(*,*) MYROW,MYCOL,'bcast recv row 1'
                     CALL IGEBR2D( ICTXT, 'Row', TOP, 8, 1, 
     $                    IBUFF, 8, RSRC1, CSRC1 )
                     IF(DBG) WRITE(*,*) MYROW,MYCOL,'IBUFF2',IBUFF(1:8)
                     I = IBUFF( 1 ) 
                     NWIN = IBUFF( 2 ) 
                     PITRAF = IBUFF( 3 ) 
                     KS = IBUFF( 4 ) 
                     PDTRAF = IBUFF( 5 ) 
                     NDTRAF = IBUFF( 6 ) 
                     ILEN = IBUFF( 7 )
                     DLEN = IBUFF( 8 )
                     BUFFLEN = ILEN + DLEN
                     IPW3 = IPW2 + NWIN*NWIN
                     DIM1 = NB - MOD(I-1,NB) 
                     DIM4 = NWIN - DIM1
                     LIHIC = NWIN + I - 1
                     SKIP1CR = WINDOW.EQ.1 .AND. 
     $                    ICEIL(LIHIC,NB).LE.ICEIL(ILO,NB)
                  END IF
                  IF( NPROW.GT.1 .AND. DIR.EQ.2 .AND. MYCOL.EQ.CSRC1 ) 
     $                 THEN
                     IF(DBG) WRITE(*,*) MYROW,MYCOL,'bcast recv col 1'
                     CALL IGEBR2D( ICTXT, 'Col', TOP, 8, 1, 
     $                    IBUFF, 8, RSRC1, CSRC1 ) 
                     IF(DBG) WRITE(*,*) MYROW,MYCOL,'IBUFF2',IBUFF(1:8)
                     I = IBUFF( 1 ) 
                     NWIN = IBUFF( 2 ) 
                     PITRAF = IBUFF( 3 ) 
                     KS = IBUFF( 4 ) 
                     PDTRAF = IBUFF( 5 ) 
                     NDTRAF = IBUFF( 6 )
                     ILEN = IBUFF( 7 )
                     DLEN = IBUFF( 8 )
                     BUFFLEN = ILEN + DLEN
                     IPW3 = IPW2 + NWIN*NWIN
                     DIM1 = NB - MOD(I-1,NB) 
                     DIM4 = NWIN - DIM1
                     LIHIC = NWIN + I - 1
                     SKIP1CR = WINDOW.EQ.1 .AND. 
     $                    ICEIL(LIHIC,NB).LE.ICEIL(ILO,NB)
                  END IF
               END IF
               IF( RSRC1.NE.RSRC4 ) THEN
                  IF( MYROW.EQ.RSRC4 .AND. MYCOL.EQ.CSRC4 ) THEN
                     IF(DBG) WRITE(*,*) MYROW,MYCOL,'bcast send row 4'
                     IF(DBG) WRITE(*,*) MYROW,MYCOL,'IBUFF2',IBUFF(1:8)
                     IF( NPCOL.GT.1 .AND. DIR.EQ.1 )
     $                    CALL IGEBS2D( ICTXT, 'Row', TOP, 8, 1, 
     $                    IBUFF, 8 )
                     SKIP1CR = WINDOW.EQ.1 .AND. 
     $                    ICEIL(LIHIC,NB).LE.ICEIL(ILO,NB)
                  ELSEIF( MYROW.EQ.RSRC4 ) THEN
                     IF( NPCOL.GT.1 .AND. DIR.EQ.1 ) THEN
                        IF(DBG) WRITE(*,*) MYROW,MYCOL,
     $                       'bcast recv row 4'
                        CALL IGEBR2D( ICTXT, 'Row', TOP, 8, 1, 
     $                       IBUFF, 8, RSRC4, CSRC4 )
                        IF(DBG) WRITE(*,*) MYROW,MYCOL,'IBUFF2',
     $                       IBUFF(1:8)
                        I = IBUFF( 1 ) 
                        NWIN = IBUFF( 2 ) 
                        PITRAF = IBUFF( 3 ) 
                        KS = IBUFF( 4 ) 
                        PDTRAF = IBUFF( 5 ) 
                        NDTRAF = IBUFF( 6 ) 
                        ILEN = IBUFF( 7 )
                        DLEN = IBUFF( 8 )
                        BUFFLEN = ILEN + DLEN
                        IPW3 = IPW2 + NWIN*NWIN
                        DIM1 = NB - MOD(I-1,NB) 
                        DIM4 = NWIN - DIM1
                        LIHIC = NWIN + I - 1
                        SKIP1CR = WINDOW.EQ.1 .AND. 
     $                       ICEIL(LIHIC,NB).LE.ICEIL(ILO,NB)
                     END IF
                  END IF
               END IF
               IF( CSRC1.NE.CSRC4 ) THEN
                  IF(DBG) WRITE(*,*) MYROW,MYCOL,'RSRC4,CSRC4=',
     $                 RSRC4,CSRC4
                  IF( MYROW.EQ.RSRC4 .AND. MYCOL.EQ.CSRC4 ) THEN
                     IF(DBG) WRITE(*,*) MYROW,MYCOL,'bcast send col 4'
                     IF(DBG) WRITE(*,*) MYROW,MYCOL,'IBUFF2',
     $                    IBUFF(1:8)
                     IF( NPROW.GT.1 .AND. DIR.EQ.2 )
     $                    CALL IGEBS2D( ICTXT, 'Col', TOP, 8, 1, 
     $                    IBUFF, 8 )
                     SKIP1CR = WINDOW.EQ.1 .AND. 
     $                    ICEIL(LIHIC,NB).LE.ICEIL(ILO,NB)
                  ELSEIF( MYCOL.EQ.CSRC4 ) THEN
                     IF( NPROW.GT.1 .AND. DIR.EQ.2 ) THEN
                        IF(DBG) WRITE(*,*) MYROW,MYCOL,
     $                       'bcast recv col 4'
                        CALL IGEBR2D( ICTXT, 'Col', TOP, 8, 1, 
     $                       IBUFF, 8, RSRC4, CSRC4 )
                        IF(DBG) WRITE(*,*) MYROW,MYCOL,'IBUFF2',
     $                       IBUFF(1:8)
                        I = IBUFF( 1 ) 
                        NWIN = IBUFF( 2 ) 
                        PITRAF = IBUFF( 3 ) 
                        KS = IBUFF( 4 ) 
                        PDTRAF = IBUFF( 5 ) 
                        NDTRAF = IBUFF( 6 ) 
                        ILEN = IBUFF( 7 )
                        DLEN = IBUFF( 8 )
                        BUFFLEN = ILEN + DLEN
                        IPW3 = IPW2 + NWIN*NWIN
                        DIM1 = NB - MOD(I-1,NB) 
                        DIM4 = NWIN - DIM1
                        LIHIC = NWIN + I - 1
                        SKIP1CR = WINDOW.EQ.1 .AND. 
     $                       ICEIL(LIHIC,NB).LE.ICEIL(ILO,NB)
                     END IF
                  END IF
               END IF
C
C     Skip rest of broadcasts and updates if appropriate
C               
               IF(DBG) WRITE(*,*) MYROW,MYCOL,
     $              'WINDOW=',WINDOW
               IF(DBG) WRITE(*,*) MYROW,MYCOL,
     $              'LIHIC=',LIHIC
               IF(DBG) WRITE(*,*) MYROW,MYCOL,
     $              'ILO=',ILO
               IF(DBG) WRITE(*,*) MYROW,MYCOL,
     $              'NB=',NB
               IF(DBG) WRITE(*,*) MYROW,MYCOL,
     $              'SKIP1CR=',SKIP1CR
               IF( SKIP1CR ) GO TO 326 
C     
C     Broadcast the orthogonal transformations
C     
               IF(DBG) WRITE(*,*) MYROW,MYCOL,'bcast orthogonal transf.'
               IF(DBG) WRITE(*,*) MYROW,MYCOL,
     $              'BUFFER,ILEN,DLEN,BUFFLEN=',
     $              BUFFER,ILEN,DLEN,BUFFLEN
               IF( MYROW.EQ.RSRC1 .AND. MYCOL.EQ.CSRC1 ) THEN
                  BUFFER = PDTRAF
                  BUFFLEN = DLEN + ILEN
                  IF( (NPROW.GT.1 .AND. DIR.EQ.2) .OR. 
     $                (NPCOL.GT.1 .AND. DIR.EQ.1) ) THEN
                     DO 370 INDX = 1, ILEN
                        DWORK( BUFFER+INDX-1 ) = 
     $                       DBLE( IWORK(IPIW+INDX-1) )
 370                 CONTINUE
                     CALL DLACPY( 'All', DLEN, 1, DWORK( IPW3 ), 
     $                    DLEN, DWORK(BUFFER+ILEN), DLEN )
                  END IF
                  IF( NPCOL.GT.1 .AND. DIR.EQ.1 ) THEN
                     IF(DBG) WRITE(*,*) MYROW,MYCOL,'bcast send row 1'
                     CALL DGEBS2D( ICTXT, 'Row', TOP, BUFFLEN, 1, 
     $                    DWORK(BUFFER), BUFFLEN )
                  END IF
                  IF( NPROW.GT.1 .AND. DIR.EQ.2 ) THEN
                     IF(DBG) WRITE(*,*) MYROW,MYCOL,'bcast send col 1'
                     IF(DBG) WRITE(*,*) MYROW,MYCOL,'BUFFER=',
     $                    DWORK(BUFFER:BUFFER+BUFFLEN-1)
                     CALL DGEBS2D( ICTXT, 'Col', TOP, BUFFLEN, 1, 
     $                    DWORK(BUFFER), BUFFLEN )
                  END IF
               ELSEIF( MYROW.EQ.RSRC1 .OR. MYCOL.EQ.CSRC1 ) THEN
                  IF( NPCOL.GT.1 .AND. DIR.EQ.1 .AND. MYROW.EQ.RSRC1 ) 
     $                 THEN
                     BUFFER = PDTRAF
                     BUFFLEN = DLEN + ILEN
                     IF(DBG) WRITE(*,*) MYROW,MYCOL,'bcast recv row 1'
                     CALL DGEBR2D( ICTXT, 'Row', TOP, BUFFLEN, 1, 
     $                    DWORK(BUFFER), BUFFLEN, RSRC1, CSRC1 )
                  END IF
                  IF( NPROW.GT.1 .AND. DIR.EQ.2 .AND. MYCOL.EQ.CSRC1 ) 
     $                 THEN
                     BUFFER = PDTRAF
                     BUFFLEN = DLEN + ILEN
                     IF(DBG) WRITE(*,*) MYROW,MYCOL,'bcast recv col 1'
                     CALL DGEBR2D( ICTXT, 'Col', TOP, BUFFLEN, 1, 
     $                    DWORK(BUFFER), BUFFLEN, RSRC1, CSRC1 )
                     IF(DBG) WRITE(*,*) MYROW,MYCOL,'BUFFER=',
     $                    DWORK(BUFFER:BUFFER+BUFFLEN-1)
                  END IF
                  IF( (NPCOL.GT.1.AND.DIR.EQ.1.AND.MYROW.EQ.RSRC1) .OR.
     $                (NPROW.GT.1.AND.DIR.EQ.2.AND. MYCOL.EQ.CSRC1) ) 
     $                 THEN
                     DO 380 INDX = 1, ILEN
                        IWORK(IPIW+INDX-1) = 
     $                       INT( DWORK( BUFFER+INDX-1 ) )
 380                 CONTINUE
                     CALL DLACPY( 'All', DLEN, 1, 
     $                    DWORK( BUFFER+ILEN ), DLEN, 
     $                    DWORK( IPW3 ), DLEN )
                  END IF
               END IF
               IF( RSRC1.NE.RSRC4 ) THEN
                  IF( MYROW.EQ.RSRC4 .AND. MYCOL.EQ.CSRC4 ) THEN
                     BUFFER = PDTRAF
                     BUFFLEN = DLEN + ILEN
                     IF( NPCOL.GT.1 .AND. DIR.EQ.1 ) THEN
                        DO 390 INDX = 1, ILEN
                           DWORK( BUFFER+INDX-1 ) = 
     $                          DBLE( IWORK(IPIW+INDX-1) )
 390                    CONTINUE
                        CALL DLACPY( 'All', DLEN, 1, DWORK( IPW3 ), 
     $                       DLEN, DWORK(BUFFER+ILEN), DLEN )
                        IF(DBG) WRITE(*,*) MYROW,MYCOL,
     $                       'bcast send row 4'
                        CALL DGEBS2D( ICTXT, 'Row', TOP, BUFFLEN, 
     $                       1, DWORK(BUFFER), BUFFLEN )
                     END IF
                  ELSEIF( MYROW.EQ.RSRC4.AND.DIR.EQ.1.AND.NPCOL.GT.1 ) 
     $                    THEN
                     IF(DBG) WRITE(*,*) MYROW,MYCOL,
     $                    'bcast recv row 4'
                     BUFFER = PDTRAF
                     BUFFLEN = DLEN + ILEN
                     CALL DGEBR2D( ICTXT, 'Row', TOP, BUFFLEN, 
     $                    1, DWORK(BUFFER), BUFFLEN, RSRC4, CSRC4 )
                     DO 400 INDX = 1, ILEN
                        IWORK(IPIW+INDX-1) = 
     $                       INT( DWORK( BUFFER+INDX-1 ) )
 400                 CONTINUE
                     CALL DLACPY( 'All', DLEN, 1, 
     $                    DWORK( BUFFER+ILEN ), DLEN, 
     $                    DWORK( IPW3 ), DLEN )
                  END IF
               END IF
               IF( CSRC1.NE.CSRC4 ) THEN
                  IF( MYROW.EQ.RSRC4 .AND. MYCOL.EQ.CSRC4 ) THEN
                     BUFFER = PDTRAF
                     BUFFLEN = DLEN + ILEN
                     IF( NPROW.GT.1 .AND. DIR.EQ.2 ) THEN 
                        DO 395 INDX = 1, ILEN
                           DWORK( BUFFER+INDX-1 ) = 
     $                          DBLE( IWORK(IPIW+INDX-1) )
 395                    CONTINUE
                        CALL DLACPY( 'All', DLEN, 1, DWORK( IPW3 ), 
     $                       DLEN, DWORK(BUFFER+ILEN), DLEN )
                        IF(DBG) WRITE(*,*) MYROW,MYCOL,
     $                       'bcast send col 4, BUFFLEN=',
     $                       BUFFLEN
                        CALL DGEBS2D( ICTXT, 'Col', TOP, BUFFLEN, 
     $                       1, DWORK(BUFFER), BUFFLEN )
                        IF(DBG) WRITE(*,*) MYROW,MYCOL,
     $                       'bcast send col 4 done'
                     END IF
                  ELSEIF( MYCOL.EQ.CSRC4.AND.DIR.EQ.2.AND.NPROW.GT.1 ) 
     $                    THEN
                     BUFFER = PDTRAF
                     BUFFLEN = DLEN + ILEN
                     IF(DBG) WRITE(*,*) MYROW,MYCOL,
     $                    'bcast recv col 4, BUFFLEN=',
     $                    BUFFLEN
                     CALL DGEBR2D( ICTXT, 'Col', TOP, BUFFLEN, 1, 
     $                    DWORK(BUFFER), BUFFLEN, RSRC4, CSRC4 )
                     IF(DBG) WRITE(*,*) MYROW,MYCOL,
     $                    'bcast recv col 4 done'
                     DO 402 INDX = 1, ILEN
                        IWORK(IPIW+INDX-1) = 
     $                       INT( DWORK( BUFFER+INDX-1 ) )
 402                 CONTINUE
                     CALL DLACPY( 'All', DLEN, 1, 
     $                    DWORK( BUFFER+ILEN ), DLEN, 
     $                    DWORK( IPW3 ), DLEN )
                  END IF
               END IF
C     
 326           CONTINUE
C     
 321        CONTINUE
C
C     Compute crossborder updates
C
            DO 322 WINDOW = WINDOW0, WINE, 2
               IF(DBG) WRITE(*,*)MYROW,MYCOL, 
     $              'compute crossborder updates, part 2'
               IF( WINDOW.EQ.1 .AND. SKIP1CR ) GO TO 327
               RSRC4 = IWORK(IRSRC+WINDOW-1)
               CSRC4 = IWORK(ICSRC+WINDOW-1)
               RSRC1 = MOD( RSRC4 - 1 + NPROW, NPROW )
               CSRC1 = MOD( CSRC4 - 1 + NPCOL, NPCOL )
C     
C     Prepare workspaces for updates: 
C       IPW3 holds now the orthogonal transformations
C       IPW4 holds the explicit orthogonal matrix, if formed
C       IPW5 holds the crossborder block column of T 
C       IPW6 holds the crossborder block row of T
C       IPW7 holds the crossborder block column of Q (if WANTQ=.TRUE.)
C       IPW8 points to the leftover workspace used as lhs in matmults 
C     
               IF(((MYCOL.EQ.CSRC1.OR.MYCOL.EQ.CSRC4).AND.DIR.EQ.2).OR.
     $            ((MYROW.EQ.RSRC1.OR.MYROW.EQ.RSRC4).AND.DIR.EQ.1)) 
     $              THEN
                  IPW4 = BUFFER 
                  IF( DIR.EQ.2 ) THEN
                     IF( WANTQ ) THEN
                        QROWS = NUMROC( N, NB, MYROW, DESCQ( RSRC_ ), 
     $                       NPROW )
                     ELSE
                        QROWS = 0
                     END IF
                     IF(DBG) WRITE(*,*)MYROW,MYCOL,
     $                    'computing TROWS for I=',I
                     TROWS = NUMROC( I-1, NB, MYROW, DESCT( RSRC_ ), 
     $                    NPROW )
                     
                  ELSE
                     QROWS = 0
                     TROWS = 0
                  END IF
                  IF( DIR.EQ.1 ) THEN
                     TCOLS = NUMROC( N - (I+DIM1-1), NB, MYCOL, CSRC4, 
     $                    NPCOL )
                     IF( MYCOL.EQ.CSRC4 ) TCOLS = TCOLS - DIM4
                  ELSE
                     TCOLS = 0
                  END IF
                  IF(DBG) WRITE(*,*) MYROW,MYCOL,'TROWS,TCOLS,QROWS=',
     $                 TROWS,TCOLS,QROWS
                  IPW5 = IPW4 + NWIN*NWIN
                  IPW6 = IPW5 + TROWS * NWIN
                  IF( WANTQ ) THEN
                     IPW7 = IPW6 + NWIN * TCOLS
                     IPW8 = IPW7 + QROWS * NWIN
                  ELSE
                     IPW8 = IPW6 + NWIN * TCOLS
                  END IF
               END IF
C     
C     Let each process row and column involved in the updates
C     exchange data in T and Q with their neighbours 
C     
               IF( DIR.EQ.2 ) THEN
               IF( MYCOL.EQ.CSRC1 .OR. MYCOL.EQ.CSRC4 ) THEN
                  DO 410 INDX = 1, NPROW
                     IF( MYCOL.EQ.CSRC1 ) THEN
                        CALL INFOG2L( 1+(INDX-1)*NB, I, DESCT, NPROW, 
     $                       NPCOL, MYROW, MYCOL, ILOC, JLOC1, RSRC, 
     $                       CSRC1 )
                        IF( MYROW.EQ.RSRC ) THEN
                           CALL DLACPY( 'All', TROWS, DIM1, 
     $                          T((JLOC1-1)*LLDT+ILOC), LLDT, 
     $                          DWORK(IPW5), TROWS )
                           IF( NPCOL.GT.1 ) THEN
                              EAST = MOD( MYCOL + 1, NPCOL )
                              IF(DBG) WRITE(*,*) MYROW,MYCOL,
     $                             'send east', RSRC, EAST
                              CALL DGESD2D( ICTXT, TROWS, DIM1,
     $                             DWORK(IPW5), TROWS, RSRC, EAST )
                              IF(DBG) WRITE(*,*) MYROW,MYCOL,'recv east'
                              CALL DGERV2D( ICTXT, TROWS, DIM4,
     $                             DWORK(IPW5+TROWS*DIM1), TROWS, RSRC, 
     $                             EAST )
                              IF(DBG) WRITE(*,*) MYROW,MYCOL,
     $                             'recv east done'
                           END IF
                        END IF
                     END IF
                     IF( MYCOL.EQ.CSRC4 ) THEN
                        CALL INFOG2L( 1+(INDX-1)*NB, I+DIM1, DESCT, 
     $                       NPROW, NPCOL, MYROW, MYCOL, ILOC, JLOC4, 
     $                       RSRC, CSRC4 )
                        IF( MYROW.EQ.RSRC ) THEN
                           CALL DLACPY( 'All', TROWS, DIM4, 
     $                          T((JLOC4-1)*LLDT+ILOC), LLDT, 
     $                          DWORK(IPW5+TROWS*DIM1), TROWS )
                           IF( NPCOL.GT.1 ) THEN
                              WEST = MOD( MYCOL - 1 + NPCOL, NPCOL )
                              IF(DBG) WRITE(*,*) MYROW,MYCOL,
     $                             'send west',RSRC, 
     $                             WEST
                              CALL DGESD2D( ICTXT, TROWS, DIM4,
     $                             DWORK(IPW5+TROWS*DIM1), TROWS, RSRC, 
     $                             WEST )
                              IF(DBG) WRITE(*,*) MYROW,MYCOL,
     $                             'send west done'
                              IF(DBG) WRITE(*,*) MYROW,MYCOL,'recv west'
                              CALL DGERV2D( ICTXT, TROWS, DIM1,
     $                             DWORK(IPW5), TROWS, RSRC, WEST )
                              IF(DBG) WRITE(*,*) MYROW,MYCOL,
     $                             'recv west done'
                           END IF
                        END IF
                     END IF
 410              CONTINUE
               END IF
               END IF
C     
               IF( DIR.EQ.1 ) THEN
               IF( MYROW.EQ.RSRC1 .OR. MYROW.EQ.RSRC4 ) THEN
                  IF(DBG) WRITE(*,*) MYROW,MYCOL,'Entering rowexch'
                  IF(DBG) WRITE(*,*) MYROW,MYCOL,'RSRC1,RSRC4=',
     $                 RSRC1,RSRC4 
                  DO 420 INDX = 1, NPCOL
                     IF(DBG) WRITE(*,*) MYROW,MYCOL,'INDX=',INDX
                     IF( MYROW.EQ.RSRC1 ) THEN
                        IF( INDX.EQ.1 ) THEN
                           CALL INFOG2L( I, LIHIC+1, DESCT, NPROW, 
     $                          NPCOL, MYROW, MYCOL, ILOC1, JLOC, 
     $                          RSRC1, CSRC )
                           IF(DBG) WRITE(*,*) MYROW,MYCOL,
     $                          'I,LIHIC,ILOC1,JLOC,RSRC1,CSRC=',
     $                          I,LIHIC,ILOC1,JLOC,RSRC1,CSRC
                        ELSE
                           CALL INFOG2L( I,
     $                          (ICEIL(LIHIC,NB)+(INDX-2))*NB+1, 
     $                          DESCT, NPROW, NPCOL, MYROW, MYCOL, 
     $                          ILOC1, JLOC, RSRC1, CSRC )
                            IF(DBG) WRITE(*,*) MYROW,MYCOL,
     $                          'ICEIL(LIHIC,NB)+(INDX-2))*NB+1=',
     $                          (ICEIL(LIHIC,NB)+(INDX-2))*NB+1
                           IF(DBG) WRITE(*,*) MYROW,MYCOL,
     $                          'I,ILOC1,JLOC,RSRC1,CSRC=',
     $                          I,ILOC1,JLOC,RSRC1,CSRC
                        END IF
                        IF( MYCOL.EQ.CSRC ) THEN
                           CALL DLACPY( 'All', DIM1, TCOLS, 
     $                          T((JLOC-1)*LLDT+ILOC1), LLDT, 
     $                          DWORK(IPW6), NWIN )
                           IF( NPROW.GT.1 ) THEN
                              SOUTH = MOD( MYROW + 1, NPROW )
                              IF(DBG) WRITE(*,*) MYROW,MYCOL,
     $                             'send south'
                              CALL DGESD2D( ICTXT, DIM1, TCOLS,
     $                             DWORK(IPW6), NWIN, SOUTH, CSRC )
                              IF(DBG) WRITE(*,*) MYROW,MYCOL,
     $                             'recv south'
                              CALL DGERV2D( ICTXT, DIM4, TCOLS,
     $                             DWORK(IPW6+DIM1), NWIN, SOUTH, 
     $                             CSRC )
                              IF(DBG) WRITE(*,*) MYROW,MYCOL,
     $                             'recv south done'
                           END IF
                        END IF
                     END IF
                     IF( MYROW.EQ.RSRC4 ) THEN
                        IF( INDX.EQ.1 ) THEN
                           CALL INFOG2L( I+DIM1, LIHIC+1, DESCT, NPROW, 
     $                          NPCOL, MYROW, MYCOL, ILOC4, JLOC, RSRC4, 
     $                          CSRC )
                           IF(DBG) WRITE(*,*) MYROW,MYCOL,
     $                          'I+DIM1,LIHIC,ILOC4,JLOC,RSRC4,CSRC=',
     $                          I+DIM1,LIHIC,ILOC4,JLOC,RSRC4,CSRC
                        ELSE
                           CALL INFOG2L( I+DIM1, 
     $                          (ICEIL(LIHIC,NB)+(INDX-2))*NB+1, DESCT, 
     $                          NPROW, NPCOL, MYROW, MYCOL, ILOC4, JLOC, 
     $                          RSRC4, CSRC )
                           IF(DBG) WRITE(*,*) MYROW,MYCOL,
     $                          'ICEIL(LIHIC,NB)+(INDX-2))*NB+1=',
     $                          (ICEIL(LIHIC,NB)+(INDX-2))*NB+1
                           IF(DBG) WRITE(*,*) MYROW,MYCOL,
     $                          'I+DIM1,ILOC4,JLOC,RSRC4,CSRC=',
     $                          I+DIM1,ILOC4,JLOC,RSRC4,CSRC
                        END IF
                        IF( MYCOL.EQ.CSRC ) THEN
                           CALL DLACPY( 'All', DIM4, TCOLS, 
     $                          T((JLOC-1)*LLDT+ILOC4), LLDT, 
     $                          DWORK(IPW6+DIM1), NWIN )
                           IF( NPROW.GT.1 ) THEN
                              NORTH = MOD( MYROW - 1 + NPROW, NPROW )
                              IF(DBG) WRITE(*,*) MYROW,MYCOL,
     $                             'send north'
                              CALL DGESD2D( ICTXT, DIM4, TCOLS,
     $                             DWORK(IPW6+DIM1), NWIN, NORTH, 
     $                             CSRC )
                              IF(DBG) WRITE(*,*) MYROW,MYCOL,
     $                             'recv north'
                              CALL DGERV2D( ICTXT, DIM1, TCOLS,
     $                             DWORK(IPW6), NWIN, NORTH, CSRC )
                              IF(DBG) WRITE(*,*) MYROW,MYCOL,
     $                             'recv north done'
                           END IF
                        END IF
                     END IF
 420              CONTINUE
               END IF
               END IF
C     
               IF( DIR.EQ.2 ) THEN
               IF( WANTQ ) THEN
                  IF( MYCOL.EQ.CSRC1 .OR. MYCOL.EQ.CSRC4 ) THEN
                     DO 430 INDX = 1, NPROW
                        IF( MYCOL.EQ.CSRC1 ) THEN
                           CALL INFOG2L( 1+(INDX-1)*NB, I, DESCQ, 
     $                          NPROW, NPCOL, MYROW, MYCOL, ILOC, JLOC1, 
     $                          RSRC, CSRC1 )
                           IF( MYROW.EQ.RSRC ) THEN
                              CALL DLACPY( 'All', QROWS, DIM1, 
     $                             Q((JLOC1-1)*LLDQ+ILOC), LLDQ, 
     $                             DWORK(IPW7), QROWS )
                              IF( NPCOL.GT.1 ) THEN
                                 EAST = MOD( MYCOL + 1, NPCOL )
                                 IF(DBG) WRITE(*,*) MYROW,MYCOL,
     $                                'send Q east'
                                 CALL DGESD2D( ICTXT, QROWS, DIM1,
     $                                DWORK(IPW7), QROWS, RSRC, EAST )
                                 IF(DBG) WRITE(*,*) MYROW,MYCOL,
     $                                'recv Q east' 
                                 CALL DGERV2D( ICTXT, QROWS, DIM4,
     $                                DWORK(IPW7+QROWS*DIM1), QROWS, 
     $                                RSRC, EAST )
                              END IF
                           END IF
                        END IF
                        IF( MYCOL.EQ.CSRC4 ) THEN
                           CALL INFOG2L( 1+(INDX-1)*NB, I+DIM1, DESCQ, 
     $                          NPROW, NPCOL, MYROW, MYCOL, ILOC, JLOC4, 
     $                          RSRC, CSRC4 )
                           IF( MYROW.EQ.RSRC ) THEN
                              CALL DLACPY( 'All', QROWS, DIM4, 
     $                             Q((JLOC4-1)*LLDQ+ILOC), LLDQ, 
     $                             DWORK(IPW7+QROWS*DIM1), QROWS )
                              IF( NPCOL.GT.1 ) THEN
                                 WEST = MOD( MYCOL - 1 + NPCOL, NPCOL )
                                 IF(DBG) WRITE(*,*) MYROW,MYCOL,
     $                                'send Q west'
                                 CALL DGESD2D( ICTXT, QROWS, DIM4,
     $                                DWORK(IPW7+QROWS*DIM1), QROWS, 
     $                                RSRC, WEST )
                                 IF(DBG) WRITE(*,*) MYROW,MYCOL,
     $                                'recv Q west'
                                 CALL DGERV2D( ICTXT, QROWS, DIM1,
     $                                DWORK(IPW7), QROWS, RSRC, WEST )
                              END IF
                           END IF
                        END IF
 430                 CONTINUE
                  END IF
               END IF
               END IF
C     
 327           CONTINUE
C     
 322        CONTINUE
C     
            DO 323 WINDOW = WINDOW0, WINE, 2
               IF(DBG) WRITE(*,*)MYROW,MYCOL, 
     $              'compute crossborder updates, part 2'
               IF(DBG) WRITE(*,*)MYROW,MYCOL,'WINDOW=',WINDOW
               IF(DBG) WRITE(*,*)MYROW,MYCOL,'DIR=',DIR
               RSRC4 = IWORK(IRSRC+WINDOW-1)
               CSRC4 = IWORK(ICSRC+WINDOW-1)
               RSRC1 = MOD( RSRC4 - 1 + NPROW, NPROW )
               CSRC1 = MOD( CSRC4 - 1 + NPCOL, NPCOL )
               IF(DBG) WRITE(*,*)MYROW,MYCOL,'RSRC1,CSRC1=',
     $              RSRC1,CSRC1
               IF(DBG) WRITE(*,*)MYROW,MYCOL,'RSRC4,CSRC4=',
     $              RSRC4,CSRC4
               FLOPS = 0
               IF(((MYCOL.EQ.CSRC1.OR.MYCOL.EQ.CSRC4).AND.DIR.EQ.2).OR.
     $            ((MYROW.EQ.RSRC1.OR.MYROW.EQ.RSRC4).AND.DIR.EQ.1)) 
     $              THEN
                  IF(DBG) WRITE(*,*)MYROW,MYCOL,'into crb updates'
C
C     Skip this part of the updates if appropriate
C
                  IF(DBG) WRITE(*,*)MYROW,MYCOL,'SKIP1CR=',
     $                 SKIP1CR   
                  IF( WINDOW.EQ.1 .AND. SKIP1CR ) GO TO 328
C     
C     Count number of operations to decide whether to use matrix-matrix 
C     multiplications for updating off-diagonal parts or not.
C     
                  NITRAF = PITRAF - IPIW
                  ISHH = .FALSE.
                  DO 405 K = 1, NITRAF
                     IF( IWORK( IPIW + K - 1 ).LE.NWIN ) THEN
                        FLOPS = FLOPS + 6
                     ELSE
                        FLOPS = FLOPS + 11
                        ISHH = .TRUE.
                     END IF
 405              CONTINUE
                  IF(DBG) WRITE(*,*)MYROW,MYCOL,'FLOPS=',
     $                 FLOPS
C     
C     Perform updates in parallel
C     
                  IF( FLOPS.NE.0 .AND.
     $                 ( 2*FLOPS*100 )/( 2*NWIN*NWIN ) .GE. MMULT ) THEN
C     
                     CALL DLASET( 'All', NWIN, NWIN, ZERO, ONE, 
     $                    DWORK( IPW4 ), NWIN )
                     IF(DBG) WRITE(*,*)MYROW,MYCOL,'bdlaapp 1 ',
     $                    'NWIN,NCB,IPW4,NITRAF,IPIW,IPW3,IPW8=',
     $                    NWIN,NCB,IPW4,NITRAF,IPIW,IPW3,IPW8 
                     IF(DBG) WRITE(*,*)MYROW,MYCOL,
     $                    'DWORK(IPW3:IPW4+NWIN*NWIN-1)=',
     $                    DWORK(IPW3:IPW4+NWIN*NWIN-1)
                     IF(DBG) WRITE(*,*)MYROW,MYCOL,'IWORK(IPIW)=',
     $                    IWORK(IPIW:IPIW+NITRAF-1)
                     DWORK(IPW8) = DBLE(MYROW)
                     DWORK(IPW8+1) = DBLE(MYCOL) 
                     CALL BDLAAPP( 1, NWIN, NWIN, NCB, DWORK( IPW4 ), 
     $                    NWIN, NITRAF, IWORK(IPIW), DWORK( IPW3 ), 
     $                    DWORK(IPW8) )
                     IF(DBG) WRITE(*,*)MYROW,MYCOL,'bdlaapp 1 done!'
C     
C     Test if sparsity structure of orthogonal matrix can be exploited
C     (see below) 
C     
                     IF( ISHH .OR. DIM1.NE.KS ) THEN
C     
C     Update the columns of T and Q affected by the reordering
C     
                        IF( DIR.EQ.2 ) THEN
                        DO 440 INDX = 1, MIN(I-1,1+(NPROW-1)*NB), NB
                           IF( MYCOL.EQ.CSRC1 ) THEN
                              CALL INFOG2L( INDX, I, DESCT, NPROW, 
     $                             NPCOL, MYROW, MYCOL, ILOC, JLOC, 
     $                             RSRC, CSRC1 )
                              IF( MYROW.EQ.RSRC ) THEN
                                 IF(DBG) WRITE(*,*)MYROW,MYCOL, 
     $                                'dgemm 13'
                                 CALL DGEMX( 'No transpose',
     $                                'No transpose', TROWS, DIM1, NWIN, 
     $                                ONE, DWORK( IPW5 ), TROWS, 
     $                                DWORK( IPW4 ), NWIN, ZERO, 
     $                                DWORK(IPW8), TROWS )
                                 CALL DLACPY( 'All', TROWS, DIM1, 
     $                                DWORK(IPW8), TROWS, 
     $                                T((JLOC-1)*LLDT+ILOC), LLDT )
                              END IF
                           END IF
                           IF( MYCOL.EQ.CSRC4 ) THEN
                              CALL INFOG2L( INDX, I+DIM1, DESCT, NPROW, 
     $                             NPCOL, MYROW, MYCOL, ILOC, JLOC, 
     $                             RSRC, CSRC4 )
                              IF( MYROW.EQ.RSRC ) THEN
                                 IF(DBG) WRITE(*,*)MYROW,MYCOL, 
     $                                'dgemm 14'
                                 CALL DGEMX( 'No transpose',
     $                                'No transpose', TROWS, DIM4, NWIN, 
     $                                ONE, DWORK( IPW5 ), TROWS, 
     $                                DWORK( IPW4+NWIN*DIM1 ), NWIN, 
     $                                ZERO, DWORK(IPW8), TROWS )
                                 CALL DLACPY( 'All', TROWS, DIM4, 
     $                                DWORK(IPW8), TROWS, 
     $                                T((JLOC-1)*LLDT+ILOC), LLDT )
                              END IF
                           END IF
 440                    CONTINUE
C
                        IF( WANTQ ) THEN
                           DO 450 INDX = 1, MIN(N,1+(NPROW-1)*NB), NB
                              IF( MYCOL.EQ.CSRC1 ) THEN
                                 CALL INFOG2L( INDX, I, DESCQ, NPROW, 
     $                                NPCOL, MYROW, MYCOL, ILOC, JLOC, 
     $                                RSRC, CSRC1 )
                                 IF( MYROW.EQ.RSRC ) THEN
                                    IF(DBG) WRITE(*,*)MYROW,MYCOL, 
     $                                   'dgemm 15'
                                    CALL DGEMX( 'No transpose', 
     $                                   'No transpose', QROWS, DIM1, 
     $                                   NWIN, ONE, DWORK( IPW7 ), 
     $                                   QROWS, DWORK( IPW4 ), NWIN, 
     $                                   ZERO, DWORK(IPW8), QROWS )
                                    CALL DLACPY( 'All', QROWS, DIM1, 
     $                                   DWORK(IPW8), QROWS, 
     $                                   Q((JLOC-1)*LLDQ+ILOC), LLDQ )
                                 END IF
                              END IF
                              IF( MYCOL.EQ.CSRC4 ) THEN
                                 CALL INFOG2L( INDX, I+DIM1, DESCQ, 
     $                                NPROW, NPCOL, MYROW, MYCOL, ILOC, 
     $                                JLOC, RSRC, CSRC4 )
                                 IF( MYROW.EQ.RSRC ) THEN
                                    IF(DBG) WRITE(*,*)MYROW,MYCOL, 
     $                                   'dgemm 16'
                                    CALL DGEMX( 'No transpose', 
     $                                   'No transpose', QROWS, DIM4, 
     $                                   NWIN, ONE, DWORK( IPW7 ), 
     $                                   QROWS, 
     $                                   DWORK( IPW4+NWIN*DIM1 ), NWIN, 
     $                                   ZERO, DWORK(IPW8), QROWS )
                                    CALL DLACPY( 'All', QROWS, DIM4, 
     $                                   DWORK(IPW8), QROWS, 
     $                                   Q((JLOC-1)*LLDQ+ILOC), LLDQ )
                                 END IF
                              END IF
 450                       CONTINUE
                        END IF
                        END IF
C     
C     Update the rows of T affected by the reordering
C     
                        IF( DIR.EQ.1 ) THEN
                        IF ( LIHIC.LT.N ) THEN
                           IF( MYROW.EQ.RSRC1.AND.MYCOL.EQ.CSRC4.AND.
     $                         MOD(LIHIC,NB).NE.0 ) THEN
                              INDX = LIHIC + 1
                              CALL INFOG2L( I, INDX, DESCT, NPROW, 
     $                             NPCOL, MYROW, MYCOL, ILOC, JLOC, 
     $                             RSRC1, CSRC4 )
                              IF(DBG) WRITE(*,*)MYROW,MYCOL, 'dgemm 17'
                              CALL DGEMX( 'Transpose', 'No Transpose', 
     $                             DIM1, TCOLS, NWIN, ONE, DWORK(IPW4), 
     $                             NWIN, DWORK( IPW6 ), NWIN, ZERO, 
     $                             DWORK(IPW8), DIM1 )
                              CALL DLACPY( 'All', DIM1, TCOLS, 
     $                             DWORK(IPW8), DIM1, 
     $                             T((JLOC-1)*LLDT+ILOC), LLDT )
                           END IF
                           IF( MYROW.EQ.RSRC4.AND.MYCOL.EQ.CSRC4.AND.
     $                         MOD(LIHIC,NB).NE.0 ) THEN
                              INDX = LIHIC + 1
                              CALL INFOG2L( I+DIM1, INDX, DESCT, NPROW, 
     $                             NPCOL, MYROW, MYCOL, ILOC, JLOC, 
     $                             RSRC4, CSRC4 )
                              IF(DBG) WRITE(*,*)MYROW,MYCOL, 'dgemm 18'
                              CALL DGEMX( 'Transpose', 'No Transpose', 
     $                             DIM4, TCOLS, NWIN, ONE, 
     $                             DWORK( IPW4+DIM1*NWIN ), NWIN, 
     $                             DWORK( IPW6), NWIN, ZERO, 
     $                             DWORK(IPW8), DIM4 )
                              CALL DLACPY( 'All', DIM4, TCOLS, 
     $                             DWORK(IPW8), DIM4, 
     $                             T((JLOC-1)*LLDT+ILOC), LLDT )
                              IF(DBG) WRITE(*,*)MYROW,MYCOL, 
     $                             'dgemm 18 done'
                           END IF
                           INDXS = ICEIL(LIHIC,NB)*NB + 1 
                           INDXE = MIN(N,INDXS+(NPCOL-2)*NB)
                           DO 460 INDX = INDXS, INDXE, NB
                              IF( MYROW.EQ.RSRC1 ) THEN
                                 CALL INFOG2L( I, INDX, DESCT, NPROW, 
     $                                NPCOL, MYROW, MYCOL, ILOC, JLOC, 
     $                                RSRC1, CSRC )
                                 IF( MYCOL.EQ.CSRC ) THEN
                                    IF(DBG) WRITE(*,*)MYROW,MYCOL, 
     $                                   'dgemm 19'
                                    CALL DGEMX( 'Transpose',
     $                                   'No Transpose', DIM1, TCOLS, 
     $                                   NWIN, ONE, DWORK( IPW4 ), NWIN, 
     $                                   DWORK( IPW6 ), NWIN, ZERO, 
     $                                   DWORK(IPW8), DIM1 )
                                    CALL DLACPY( 'All', DIM1, TCOLS, 
     $                                   DWORK(IPW8), DIM1, 
     $                                   T((JLOC-1)*LLDT+ILOC), LLDT )
                                    IF(DBG) WRITE(*,*)MYROW,MYCOL, 
     $                                   'dgemm 19 done'
                                 END IF
                              END IF
                              IF( MYROW.EQ.RSRC4 ) THEN
                                 CALL INFOG2L( I+DIM1, INDX, DESCT, 
     $                                NPROW, NPCOL, MYROW, MYCOL, ILOC, 
     $                                JLOC, RSRC4, CSRC )
                                 IF( MYCOL.EQ.CSRC ) THEN
                                    IF(DBG) WRITE(*,*)MYROW,MYCOL, 
     $                                   'dgemm 20 DIM4=',DIM4
                                    CALL DGEMX( 'Transpose',
     $                                   'No Transpose', DIM4, TCOLS, 
     $                                   NWIN, ONE, 
     $                                   DWORK( IPW4+NWIN*DIM1 ), NWIN, 
     $                                   DWORK( IPW6 ), NWIN, 
     $                                   ZERO, DWORK(IPW8), DIM4 )
                                    CALL DLACPY( 'All', DIM4, TCOLS, 
     $                                   DWORK(IPW8), DIM4, 
     $                                   T((JLOC-1)*LLDT+ILOC), LLDT )
                                 END IF
                              END IF
 460                       CONTINUE
                        END IF
                        END IF
                     ELSE 
C     
C     The NWIN-by-NWIN matrix U containing the accumulated orthogonal 
C     transformations has the following structure:
C     
C                   [ U11  U12 ]
C               U = [          ],
C                   [ U21  U22 ]
C     
C     where U21 is KS-by-KS upper triangular and U12 is
C     (NWIN-KS)-by-(NWIN-KS) lower triangular. For reordering over
C     the border the structure is only exploited when the border cuts
C     the columns of U conformally with the structure itself. This 
C     happens exactly when all eigenvalues in the subcluster was moved 
C     to the other side of the border and fits perfectly in their
C     new positions, i.e., the reordering stops when the last eigenvalue
C     to cross the border is reordered to the position closest to the
C     border. Tested by checking is KS = DIM1 (see above). This should
C     hold quite often...
C     
C     Update the columns of T and Q affected by the reordering
C     
C     Compute T2*U21 + T1*U11 on the left side of the border.
C     
                        IF( DIR.EQ.2 ) THEN
                        IF(DBG) WRITE(*,*) MYROW,MYCOL,
     $                       'I,LIHIC,TROWS,TCOLS,QROWS=',
     $                       I,LIHIC,TROWS,TCOLS,QROWS
                        INDXE = MIN(I-1,1+(NPROW-1)*NB)
                        DO 470 INDX = 1, INDXE, NB
                           IF( MYCOL.EQ.CSRC1 ) THEN
                              CALL INFOG2L( INDX, I, DESCT, NPROW, 
     $                             NPCOL, MYROW, MYCOL, ILOC, JLOC, 
     $                             RSRC, CSRC1 )
                              IF( MYROW.EQ.RSRC ) THEN
                                 CALL DLACPY( 'All', TROWS, KS, 
     $                                DWORK( IPW5+TROWS*DIM4), TROWS,
     $                                DWORK(IPW8), TROWS )
                                 CALL DTRMX( 'Right', 'Upper', 
     $                                'No transpose',
     $                                'No transpose', TROWS, KS, ONE,
     $                                DWORK( IPW4+DIM4 ), NWIN, 
     $                                DWORK(IPW8), TROWS )
                                 IF(DBG) WRITE(*,*)MYROW,MYCOL, 
     $                                'dgemm 21',
     $                                ' INDX=',INDX
                                 CALL DGEMX( 'No transpose', 
     $                                'No transpose', TROWS, KS, DIM4, 
     $                                ONE, DWORK( IPW5 ), TROWS, 
     $                                DWORK( IPW4 ), NWIN, ONE, 
     $                                DWORK(IPW8), TROWS )
                                 CALL DLACPY( 'All', TROWS, KS, 
     $                                DWORK(IPW8), TROWS, 
     $                                T((JLOC-1)*LLDT+ILOC), LLDT )
                              END IF
                           END IF
C     
C     Compute T1*U12 + T2*U22 on the right side of the border.
C     
                           IF( MYCOL.EQ.CSRC4 ) THEN
                              CALL INFOG2L( INDX, I+DIM1, DESCT, NPROW, 
     $                             NPCOL, MYROW, MYCOL, ILOC, JLOC, 
     $                             RSRC, CSRC4 )
                              IF( MYROW.EQ.RSRC ) THEN
                                 CALL DLACPY( 'All', TROWS, DIM4, 
     $                                DWORK(IPW5), TROWS, DWORK( IPW8 ), 
     $                                TROWS )
                                 CALL DTRMX( 'Right', 'Lower', 
     $                                'No transpose',
     $                                'No transpose', TROWS, DIM4, ONE,
     $                                DWORK( IPW4+NWIN*KS ), NWIN,
     $                                DWORK( IPW8 ), TROWS )
                                 IF(DBG) WRITE(*,*)MYROW,MYCOL, 
     $                                'dgemm 22',
     $                                ' INDX=',INDX
                                 CALL DGEMX( 'No transpose', 
     $                                'No transpose', TROWS, DIM4, KS, 
     $                                ONE, DWORK( IPW5+TROWS*DIM4), 
     $                                TROWS,
     $                                DWORK( IPW4+NWIN*KS+DIM4 ), NWIN,
     $                                ONE, DWORK( IPW8 ), TROWS )
                                 CALL DLACPY( 'All', TROWS, DIM4, 
     $                                DWORK(IPW8), TROWS, 
     $                                T((JLOC-1)*LLDT+ILOC), LLDT )
                              END IF 
                           END IF
 470                    CONTINUE
                        IF( WANTQ ) THEN
C     
C     Compute Q2*U21 + Q1*U11 on the left side of border.
C     
                           INDXE = MIN(N,1+(NPROW-1)*NB)
                           DO 480 INDX = 1, INDXE, NB
                              IF( MYCOL.EQ.CSRC1 ) THEN
                                 CALL INFOG2L( INDX, I, DESCQ, NPROW, 
     $                                NPCOL, MYROW, MYCOL, ILOC, JLOC, 
     $                                RSRC, CSRC1 )
                                 IF( MYROW.EQ.RSRC ) THEN
                                    CALL DLACPY( 'All', QROWS, KS, 
     $                                   DWORK( IPW7+QROWS*DIM4), 
     $                                   QROWS, DWORK(IPW8), QROWS )
                                    CALL DTRMX( 'Right', 'Upper', 
     $                                   'No transpose',
     $                                   'No transpose', QROWS, KS, ONE,
     $                                   DWORK( IPW4+DIM4 ), NWIN, 
     $                                   DWORK(IPW8), QROWS )
                                    IF(DBG) WRITE(*,*)MYROW,MYCOL, 
     $                                   'dgemm 23',
     $                                   ' INDX=',INDX
                                    CALL DGEMX( 'No transpose', 
     $                                   'No transpose', QROWS, KS, 
     $                                   DIM4, ONE, DWORK( IPW7 ), 
     $                                   QROWS, DWORK( IPW4 ), NWIN, 
     $                                   ONE, DWORK(IPW8), QROWS )
                                    CALL DLACPY( 'All', QROWS, KS, 
     $                                   DWORK(IPW8), QROWS, 
     $                                   Q((JLOC-1)*LLDQ+ILOC), LLDQ )
                                 END IF
                              END IF
C     
C     Compute Q1*U12 + Q2*U22 on the right side of border.
C     
                              IF( MYCOL.EQ.CSRC4 ) THEN
                                 CALL INFOG2L( INDX, I+DIM1, DESCQ, 
     $                                NPROW, NPCOL, MYROW, MYCOL, ILOC, 
     $                                JLOC, RSRC, CSRC4 )
                                 IF( MYROW.EQ.RSRC ) THEN
                                    CALL DLACPY( 'All', QROWS, DIM4, 
     $                                   DWORK(IPW7), QROWS, 
     $                                   DWORK( IPW8 ), QROWS )
                                    CALL DTRMX( 'Right', 'Lower', 
     $                                   'No transpose',
     $                                   'No transpose', QROWS, DIM4, 
     $                                   ONE, DWORK( IPW4+NWIN*KS ), 
     $                                   NWIN, DWORK( IPW8 ), QROWS )
                                    IF(DBG) WRITE(*,*)MYROW,MYCOL, 
     $                                   'dgemm 24',
     $                                   ' INDX=',INDX
                                    CALL DGEMX( 'No transpose', 
     $                                   'No transpose', QROWS, DIM4, 
     $                                   KS, ONE, 
     $                                   DWORK( IPW7+QROWS*(DIM4)), 
     $                                   QROWS,
     $                                   DWORK( IPW4+NWIN*KS+DIM4 ), 
     $                                   NWIN, ONE, DWORK( IPW8 ), 
     $                                   QROWS )
                                    CALL DLACPY( 'All', QROWS, DIM4, 
     $                                   DWORK(IPW8), QROWS, 
     $                                   Q((JLOC-1)*LLDQ+ILOC), LLDQ )
                                 END IF 
                              END IF
 480                       CONTINUE
                        END IF
                        END IF
C     
                        IF( DIR.EQ.1 ) THEN
                        IF ( LIHIC.LT.N ) THEN
C     
C     Compute U21**T*T2 + U11**T*T1 on the upper side of the border.
C     
                           IF( MYROW.EQ.RSRC1.AND.MYCOL.EQ.CSRC4.AND.
     $                         MOD(LIHIC,NB).NE.0 ) THEN
                              INDX = LIHIC + 1
                              CALL INFOG2L( I, INDX, DESCT, NPROW, 
     $                             NPCOL, MYROW, MYCOL, ILOC, JLOC, 
     $                             RSRC1, CSRC4 ) 
                              CALL DLACPY( 'All', KS, TCOLS,
     $                             DWORK( IPW6+DIM4 ), NWIN, 
     $                             DWORK(IPW8), KS )
                              CALL DTRMX( 'Left', 'Upper', 'Transpose',
     $                             'No transpose', KS, TCOLS, ONE,
     $                             DWORK( IPW4+DIM4 ), NWIN, 
     $                             DWORK(IPW8), KS )
                              IF(DBG) WRITE(*,*)MYROW,MYCOL, 'dgemm 25,
     $                             INDX=',
     $                             INDX
                              CALL DGEMX( 'Transpose', 'No transpose', 
     $                             KS, TCOLS, DIM4, ONE, DWORK(IPW4), 
     $                             NWIN, DWORK(IPW6), NWIN, 
     $                             ONE, DWORK(IPW8), KS )
                              CALL DLACPY( 'All', KS, TCOLS, 
     $                             DWORK(IPW8), KS, 
     $                             T((JLOC-1)*LLDT+ILOC), LLDT )
                           END IF
C     
C     Compute U12**T*T1 + U22**T*T2 one the lower side of the border.
C     
                           IF( MYROW.EQ.RSRC4.AND.MYCOL.EQ.CSRC4.AND.
     $                         MOD(LIHIC,NB).NE.0 ) THEN
                              INDX = LIHIC + 1
                              CALL INFOG2L( I+DIM1, INDX, DESCT, NPROW, 
     $                             NPCOL, MYROW, MYCOL, ILOC, JLOC, 
     $                             RSRC4, CSRC4 )
                              CALL DLACPY( 'All', DIM4, TCOLS,
     $                             DWORK( IPW6 ), NWIN, 
     $                             DWORK( IPW8 ), DIM4 )
                              CALL DTRMX( 'Left', 'Lower', 'Transpose',
     $                             'No transpose', DIM4, TCOLS, ONE,
     $                             DWORK( IPW4+NWIN*KS ), NWIN,
     $                             DWORK( IPW8 ), DIM4 )
                              IF(DBG) WRITE(*,*)MYROW,MYCOL, 'dgemm 26, 
     $                             INDX=',
     $                             INDX
                              CALL DGEMX( 'Transpose', 'No Transpose', 
     $                             DIM4, TCOLS, KS, ONE,
     $                             DWORK( IPW4+NWIN*KS+DIM4 ), NWIN,
     $                             DWORK( IPW6+DIM1 ), NWIN,
     $                             ONE, DWORK( IPW8), DIM4 )
                              CALL DLACPY( 'All', DIM4, TCOLS, 
     $                             DWORK(IPW8), DIM4, 
     $                             T((JLOC-1)*LLDT+ILOC), LLDT )
                           END IF
C     
C     Compute U21**T*T2 + U11**T*T1 on upper side on border.
C     
                           INDXS = ICEIL(LIHIC,NB)*NB+1
                           INDXE = MIN(N,INDXS+(NPCOL-2)*NB)
                           DO 490 INDX = INDXS, INDXE, NB
                              IF( MYROW.EQ.RSRC1 ) THEN
                                 CALL INFOG2L( I, INDX, DESCT, NPROW, 
     $                                NPCOL, MYROW, MYCOL, ILOC, JLOC, 
     $                                RSRC1, CSRC ) 
                                 IF( MYCOL.EQ.CSRC ) THEN
                                    CALL DLACPY( 'All', KS, TCOLS,
     $                                   DWORK( IPW6+DIM4 ), NWIN, 
     $                                   DWORK(IPW8), KS )
                                    CALL DTRMX( 'Left', 'Upper', 
     $                                   'Transpose', 'No transpose', 
     $                                   KS, TCOLS, ONE, 
     $                                   DWORK( IPW4+DIM4 ), NWIN, 
     $                                   DWORK(IPW8), KS )
                                    IF(DBG) WRITE(*,*)MYROW,MYCOL, 
     $                                   'dgemm 27',
     $                                   ' INDX=',INDX
                                    CALL DGEMX( 'Transpose',
     $                                   'No transpose', KS, TCOLS, 
     $                                   DIM4, ONE, DWORK(IPW4), NWIN, 
     $                                   DWORK(IPW6), NWIN, ONE, 
     $                                   DWORK(IPW8), KS )
                                    CALL DLACPY( 'All', KS, TCOLS, 
     $                                   DWORK(IPW8), KS, 
     $                                   T((JLOC-1)*LLDT+ILOC), LLDT )
                                 END IF
                              END IF
C     
C     Compute U12**T*T1 + U22**T*T2 on lower side of border.
C     
                              IF( MYROW.EQ.RSRC4 ) THEN
                                 CALL INFOG2L( I+DIM1, INDX, DESCT, 
     $                                NPROW, NPCOL, MYROW, MYCOL, ILOC, 
     $                                JLOC, RSRC4, CSRC )
                                 IF( MYCOL.EQ.CSRC ) THEN
                                    CALL DLACPY( 'All', DIM4, TCOLS,
     $                                   DWORK( IPW6 ), NWIN, 
     $                                   DWORK( IPW8 ), DIM4 )
                                    CALL DTRMX( 'Left', 'Lower', 
     $                                   'Transpose','No transpose', 
     $                                   DIM4, TCOLS, ONE,
     $                                   DWORK( IPW4+NWIN*KS ), NWIN,
     $                                   DWORK( IPW8 ), DIM4 )
                                    IF(DBG) WRITE(*,*)MYROW,MYCOL, 
     $                                   'dgemm 28',
     $                                   ' INDX=',INDX 
                                    CALL DGEMX( 'Transpose',
     $                                   'No Transpose', DIM4, TCOLS, 
     $                                   KS, ONE,
     $                                   DWORK( IPW4+NWIN*KS+DIM4 ), 
     $                                   NWIN, DWORK( IPW6+DIM1 ), 
     $                                   NWIN, ONE, DWORK( IPW8), 
     $                                   DIM4 )
                                    CALL DLACPY( 'All', DIM4, TCOLS, 
     $                                   DWORK(IPW8), DIM4, 
     $                                   T((JLOC-1)*LLDT+ILOC), LLDT )
                                 END IF
                              END IF
 490                       CONTINUE
                        END IF
                        END IF 
                     END IF
                  ELSEIF( FLOPS.NE.0 ) THEN
C     
C     Update off-diagonal blocks and Q using the pipelined elementary 
C     transformations. Now we have a delicate problem: how to do this
C     without redundant work? For now, we let the processes involved
C     compute the whole crossborder block rows and column saving only
C     the part belonging to the corresponding side of the border. To
C     make this a realistic alternative, we have modified the ratio
C     r_flops (see Kressner paper above) to give more favor to the
C     ordinary matrix multiplication.
C     
                     IF( DIR.EQ.2 ) THEN
                     IF(DBG) WRITE(*,*)MYROW,MYCOL, 
     $                    'factorized crb updates'
                     INDXE =  MIN(I-1,1+(NPROW-1)*NB)
                     DO 500 INDX = 1, INDXE, NB
                        CALL INFOG2L( INDX, I, DESCT, NPROW, NPCOL,
     $                       MYROW, MYCOL, ILOC, JLOC, RSRC, CSRC )
                        IF( MYROW.EQ.RSRC .AND. MYCOL.EQ.CSRC ) THEN
                           CALL BDLAAPP( 1, TROWS, NWIN, NCB, 
     $                          DWORK(IPW5), TROWS, NITRAF,
     $                          IWORK(IPIW), DWORK( IPW3 ), 
     $                          DWORK(IPW8) )
                           CALL DLACPY( 'All', TROWS, DIM1, 
     $                          DWORK(IPW5), TROWS,
     $                          T((JLOC-1)*LLDT+ILOC ), LLDT )
                        END IF
                        CALL INFOG2L( INDX, I+DIM1, DESCT, NPROW, 
     $                       NPCOL, MYROW, MYCOL, ILOC, JLOC, RSRC, 
     $                       CSRC )
                        IF( MYROW.EQ.RSRC .AND. MYCOL.EQ.CSRC ) THEN
                           IF( NPCOL.GT.1 )
     $                          CALL BDLAAPP( 1, TROWS, NWIN, NCB, 
     $                          DWORK(IPW5), TROWS, NITRAF,
     $                          IWORK(IPIW), DWORK( IPW3 ), 
     $                          DWORK(IPW8) )
                           CALL DLACPY( 'All', TROWS, DIM4, 
     $                          DWORK(IPW5+TROWS*DIM1), TROWS,
     $                          T((JLOC-1)*LLDT+ILOC ), LLDT )
                        END IF
 500                 CONTINUE
                     IF( WANTQ ) THEN
                        INDXE = MIN(N,1+(NPROW-1)*NB)
                        DO 510 INDX = 1, INDXE, NB
                           CALL INFOG2L( INDX, I, DESCQ, NPROW, NPCOL,
     $                          MYROW, MYCOL, ILOC, JLOC, RSRC, CSRC )
                           IF( MYROW.EQ.RSRC .AND. MYCOL.EQ.CSRC ) THEN
                              CALL BDLAAPP( 1, QROWS, NWIN, NCB, 
     $                             DWORK(IPW7), QROWS, NITRAF,
     $                             IWORK(IPIW), DWORK( IPW3 ), 
     $                             DWORK(IPW8) )
                              CALL DLACPY( 'All', QROWS, DIM1, 
     $                             DWORK(IPW7), QROWS,
     $                             Q((JLOC-1)*LLDQ+ILOC ), LLDQ )
                           END IF
                           CALL INFOG2L( INDX, I+DIM1, DESCQ, NPROW, 
     $                          NPCOL, MYROW, MYCOL, ILOC, JLOC, RSRC, 
     $                          CSRC )
                           IF( MYROW.EQ.RSRC .AND. MYCOL.EQ.CSRC ) THEN
                              IF( NPCOL.GT.1 ) 
     $                             CALL BDLAAPP( 1, QROWS, NWIN, NCB, 
     $                             DWORK(IPW7), QROWS, NITRAF,
     $                             IWORK(IPIW), DWORK( IPW3 ), 
     $                             DWORK(IPW8) )
                              CALL DLACPY( 'All', QROWS, DIM4, 
     $                             DWORK(IPW7+QROWS*DIM1), QROWS,
     $                             Q((JLOC-1)*LLDQ+ILOC ), LLDQ )
                           END IF  
 510                    CONTINUE
                     END IF
                     END IF
C     
                     IF( DIR.EQ.1 ) THEN
                     IF( LIHIC.LT.N ) THEN
                        INDX = LIHIC + 1
                        CALL INFOG2L( I, INDX, DESCT, NPROW, NPCOL,
     $                       MYROW, MYCOL, ILOC, JLOC, RSRC, CSRC )
                        IF( MYROW.EQ.RSRC .AND. MYCOL.EQ.CSRC.AND.
     $                      MOD(LIHIC,NB).NE.0 ) THEN
                           CALL BDLAAPP( 0, NWIN, TCOLS, NCB, 
     $                          DWORK( IPW6 ), NWIN, NITRAF, 
     $                          IWORK(IPIW), DWORK( IPW3 ), 
     $                          DWORK(IPW8) )
                           CALL DLACPY( 'All', DIM1, TCOLS, 
     $                          DWORK( IPW6 ), NWIN, 
     $                          T((JLOC-1)*LLDT+ILOC), LLDT )
                        END IF
                        CALL INFOG2L( I+DIM1, INDX, DESCT, NPROW, 
     $                       NPCOL, MYROW, MYCOL, ILOC, JLOC, RSRC, 
     $                       CSRC )
                        IF( MYROW.EQ.RSRC .AND. MYCOL.EQ.CSRC.AND.
     $                      MOD(LIHIC,NB).NE.0 ) THEN
                           IF( NPROW.GT.1 )
     $                          CALL BDLAAPP( 0, NWIN, TCOLS, NCB, 
     $                          DWORK( IPW6 ), NWIN, NITRAF, 
     $                          IWORK(IPIW), DWORK( IPW3 ), 
     $                          DWORK(IPW8) )
                           CALL DLACPY( 'All', DIM4, TCOLS, 
     $                          DWORK( IPW6+DIM1 ), NWIN, 
     $                          T((JLOC-1)*LLDT+ILOC), LLDT )
                        END IF
                        INDXS = ICEIL(LIHIC,NB)*NB + 1
                        INDXE = MIN(N,INDXS+(NPCOL-2)*NB)
                        DO 520 INDX = INDXS, INDXE, NB
                           CALL INFOG2L( I, INDX, DESCT, NPROW, 
     $                          NPCOL, MYROW, MYCOL, ILOC, JLOC, 
     $                          RSRC, CSRC )
                           IF( MYROW.EQ.RSRC .AND. MYCOL.EQ.CSRC ) THEN
                              CALL BDLAAPP( 0, NWIN, TCOLS, NCB, 
     $                             DWORK(IPW6), NWIN, NITRAF, 
     $                             IWORK(IPIW), DWORK( IPW3 ), 
     $                             DWORK(IPW8) )
                              CALL DLACPY( 'All', DIM1, TCOLS, 
     $                             DWORK( IPW6 ), NWIN, 
     $                             T((JLOC-1)*LLDT+ILOC), LLDT )
                           END IF
                           CALL INFOG2L( I+DIM1, INDX, DESCT, NPROW, 
     $                          NPCOL, MYROW, MYCOL, ILOC, JLOC, 
     $                          RSRC, CSRC )
                           IF( MYROW.EQ.RSRC .AND. MYCOL.EQ.CSRC ) THEN
                              IF( NPROW.GT.1 )
     $                             CALL BDLAAPP( 0, NWIN, TCOLS, NCB, 
     $                             DWORK(IPW6), NWIN, NITRAF, 
     $                             IWORK(IPIW), DWORK( IPW3 ), 
     $                             DWORK(IPW8) )
                              CALL DLACPY( 'All', DIM4, TCOLS, 
     $                             DWORK( IPW6+DIM1 ), NWIN, 
     $                             T((JLOC-1)*LLDT+ILOC), LLDT )
                           END IF
 520                    CONTINUE
                     END IF
                     END IF
                  END IF
               END IF
C     
 328           CONTINUE
               IF(DBG) WRITE(*,*)MYROW,MYCOL,'RSRC1,CSRC1=',
     $              RSRC1,CSRC1
               IF(DBG) WRITE(*,*)MYROW,MYCOL,'RSRC4,CSRC4=',
     $              RSRC4,CSRC4
C
 323        CONTINUE
C
C     End of loops over directions (DIR)
C
 2222       CONTINUE
C
C     End of loops over diagonal blocks for reordering over the block
C     diagonal
C     
 310     CONTINUE
         LAST = LAST + 1
         IF( LASTWAIT .AND. LAST.LT.2 ) GO TO 308
C     
C     Barrier to collect the processes before proceeding
C     
         IF(DBG3) WRITE(*,*) MYROW,MYCOL,'barrier 2 in'
         CALL BLACS_BARRIER( ICTXT, 'All' )
         IF(DBG3) WRITE(*,*) MYROW,MYCOL,'barrier 2 out'
C
         TIME2 = TIME2 + MPI_WTIME() - TIME
C     
C     Compute global maximum of IERR so that we know if some process
C     experienced a failure in the reordering
C     
         MYIERR = IERR
         IF( NPROCS.GT.1 ) 
     $        CALL IGAMX2D( ICTXT, 'All', TOP, 1, 1, IERR, 1, -1, 
     $        -1, -1, -1, -1 )
C     
         IF( IERR.EQ.1 .OR. IERR.EQ.2 ) THEN
            IF(DBG) WRITE(*,*) MYROW,MYCOL,
     $        'an error has occured '
C     
C     When calling BDTREXC, the block at position I+KKS-1 failed to swap.
C     
            IF( IERR.EQ.1 .OR. IERR.EQ.2 ) THEN
               INFO = MAX(1,I+KKS-1)
            ELSE 
               INFO = 1
            END IF
            IF( WANTS )
     $           S = ZERO
            IF( WANTSP )
     $           SEP = ZERO
            GO TO 300
         END IF
C     
C     Do a global update of the SELECT vector
C     
         DO 530 K = 1, N
            RSRC = INDXG2P( K, NB, MYROW, DESCT( RSRC_ ), NPROW )
            CSRC = INDXG2P( K, NB, MYCOL, DESCT( CSRC_ ), NPCOL )
            IF( MYROW.NE.RSRC .OR. MYCOL.NE.CSRC ) 
     $           SELECT( K ) = 0
 530     CONTINUE 
         IF(DBG3) WRITE(*,*) MYROW,MYCOL,'global select update'
         IF( NPROCS.GT.1 )
     $        CALL IGSUM2D( ICTXT, 'All', TOP, N, 1, SELECT, N, -1, -1 )
C     
C     Find the global minumum of ILO and IHI
C     
         ILO = ILO - 1
 523     CONTINUE
         ILO = ILO + 1
         IF( ILO.LE.N ) THEN
            IF( INT2LG(SELECT(ILO)) ) GO TO 523
         END IF
         IF(DBG3) WRITE(*,*)MYROW,MYCOL, 'new global value of ILO=',ILO
         IHI = IHI + 1
 527     CONTINUE
         IHI = IHI - 1
         IF( IHI.GE.1 ) THEN
            IF( .NOT. INT2LG(SELECT(IHI)) ) GO TO 527
         END IF
         IF(DBG3) WRITE(*,*)MYROW,MYCOL, 'new global value of IHI=',IHI
C     
C     End While ( ILO <= M )
         IF(DBG3) WRITE(*,*)MYROW,MYCOL, 'returning to 50'
         GO TO 50
      END IF
C     
      IF( WANTS ) THEN
          TIME = MPI_WTIME()
C     
C     Solve Sylvester equation T11*R - R*T2 = scale*T12 for R in 
C     parallel
C     
C     Copy T12 to workspace
C     
         IF(DBG) WRITE(*,*) MYROW,MYCOL,'T11*R - R*T21 = scale*T12'
         IF(DBG) WRITE(*,*) MYROW,MYCOL,'N1,N2=',N1,N2
         CALL INFOG2L( 1, N1+1, DESCT, NPROW, NPCOL, MYROW, 
     $        MYCOL, ILOC1, JLOC1, TRSRC, TCSRC )
         ICOFFT12 = MOD( N1, NB )
         T12ROWS = NUMROC( N1, NB, MYROW, TRSRC, NPROW )
         T12COLS = NUMROC( N2+ICOFFT12, NB, MYCOL, TCSRC, NPCOL )
         CALL DESCINIT( DESCT12, N1, N2+ICOFFT12, NB, NB, TRSRC, 
     $        TCSRC, ICTXT, MAX(1,T12ROWS), IERR )
         CALL PDLACPY( 'All', N1, N2, T, 1, N1+1, DESCT, DWORK, 
     $        1, 1+ICOFFT12, DESCT12 )
C     
C     Solve the equation to get the solution in workspace 
C     
         SPACE = DESCT12( LLD_ ) * T12COLS
         IPW1 = 1 + SPACE
         CALL PGESYCTD( 'Solve', 'Schur', 'Schur', 'Notranspose', 
     $        'Notranspose', -1, 'Demand', N1, N2, T, 1, 1, DESCT, T, 
     $        N1+1, N1+1, DESCT, DWORK, 1, 1+ICOFFT12, DESCT12, MBNB2, 
     $        DWORK(IPW1), LDWORK-SPACE+1, IWORK, LIWORK, NOEXSY, 
     $        SCALE, IERR )
         IF( IERR.LT.0 ) THEN
            INFO = 333
         ELSE
            INFO = 2
         END IF 
C     
C     Estimate the reciprocal of the condition number of the cluster
C     of eigenvalues.
C     
         RNORM = PDLANGE( 'Frobenius', N1, N2, DWORK, 1, 1+ICOFFT12, 
     $        DESCT12, DPDUM1 )
         IF( RNORM.EQ.ZERO ) THEN
            S = ONE
         ELSE
            S = SCALE / ( SQRT( SCALE*SCALE / RNORM+RNORM )*
     $           SQRT( RNORM ) )
         END IF
         TIME3 = MPI_WTIME() - TIME
      END IF
C     
      IF( WANTSP ) THEN
         TIME = MPI_WTIME()
C     
C     Estimate sep(T11,T21) in parallel.
C     
         IF(DBG) WRITE(*,*) MYROW,MYCOL,'Estimate sep(T11,T21)'
         CALL PSYCTCON( 'Notranspose', 'Notranspose', -1, 'Demand', N1, 
     $        N2, T, 1, 1, DESCT, T, N1+1, N1+1, DESCT, MBNB2, DWORK, 
     $        LDWORK, IWORK, LIWORK, EST, ITER, IERR )
         EST = EST * SQRT(DBLE(N1*N2))
         SEP = ONE / EST
         IF( IERR.LT.0 ) THEN
            INFO = 444
         ELSE
            INFO = 2
         END IF
         TIME4 = MPI_WTIME() - TIME
      END IF
C     
 300  CONTINUE
C
C     In case an error occured, do an additional global update of SELECT
C
      IF(DBG) WRITE(*,*) MYROW,MYCOL,'Additional SELECT update INFO=',
     $     INFO
      IF( INFO.NE.0 ) THEN
         DO 540 K = 1, N
            RSRC = INDXG2P( K, NB, MYROW, DESCT( RSRC_ ), NPROW )
            CSRC = INDXG2P( K, NB, MYCOL, DESCT( CSRC_ ), NPCOL )
            IF( MYROW.NE.RSRC .OR. MYCOL.NE.CSRC ) 
     $           SELECT( K ) = 0
 540     CONTINUE
         IF(DBG) WRITE(*,*) MYROW,MYCOL,'SELECT=',SELECT(1:N)
         IF(DBG) WRITE(*,*) MYROW,MYCOL,'IGSUM2D on SELECT',NPROCS
         IF( NPROCS.GT.1 )
     $        CALL IGSUM2D( ICTXT, 'All', TOP, N, 1, SELECT, N, -1, -1 )
      END IF   
C
 545  CONTINUE
C     
C     Store the output eigenvalues in WR and WI: first let all the
C     processes compute the eigenvalue inside their diagonal blocks
C     in parallel, except for the eigenvalue located next to a block
C     border. After that, compute all eigenvalues located next to the
C     block borders. Finally, do a global summation over WR and WI
C     so that all processors receive the result. Notice: real eigenvalues
C     extracted from a non-canonical 2-by-2 block are not stored in
C     in any particular order
C     
      IF(DBG) WRITE(*,*) MYROW,MYCOL,'Store output eigenvalues'
      DO 550 K = 1, N
         WR( K ) = ZERO
         WI( K ) = ZERO
 550  CONTINUE
      PAIR = .FALSE.
      DO 560 K = 1, N
         IF( .NOT. PAIR ) THEN
            IF( K.LT.N ) THEN
               BORDER = MOD( K, NB ).EQ.0 .OR. ( K.NE.1 .AND. 
     $              MOD( K, NB ).EQ.1 )
               IF( .NOT. BORDER ) THEN
                  CALL INFOG2L( K, K, DESCT, NPROW, NPCOL, MYROW, MYCOL,
     $                 ILOC1, JLOC1, TRSRC1, TCSRC1 )
                  IF( MYROW.EQ.TRSRC1 .AND. MYCOL.EQ.TCSRC1 ) THEN
                     ELEM1 = T((JLOC1-1)*LLDT+ILOC1)
                     ELEM3 = T((JLOC1-1)*LLDT+ILOC1+1)
                     IF( ELEM3.NE.ZERO ) THEN
                        ELEM2 = T((JLOC1)*LLDT+ILOC1)
                        ELEM4 = T((JLOC1)*LLDT+ILOC1+1)
                        CALL DLANV2( ELEM1, ELEM2, ELEM3, ELEM4, 
     $                       WR( K ), WI( K ), WR( K+1 ), WI( K+1 ), SN, 
     $                       CS )
                        PAIR = .TRUE.
                     ELSE 
                        IF( K.GT.1 ) THEN
                           TMP = T((JLOC1-2)*LLDT+ILOC1)
                           IF( TMP.NE.ZERO ) THEN
                              ELEM1 = T((JLOC1-2)*LLDT+ILOC1-1)
                              ELEM2 = T((JLOC1-1)*LLDT+ILOC1-1)
                              ELEM3 = T((JLOC1-2)*LLDT+ILOC1)
                              ELEM4 = T((JLOC1-1)*LLDT+ILOC1)
                              CALL DLANV2( ELEM1, ELEM2, ELEM3, ELEM4, 
     $                             WR( K-1 ), WI( K-1 ), WR( K ), 
     $                             WI( K ), SN, CS )
                           ELSE
                              WR( K ) = ELEM1 
                           END IF
                        ELSE
                           WR( K ) = ELEM1
                        END IF
                     END IF
                  END IF
               END IF
            ELSE
               CALL INFOG2L( K, K, DESCT, NPROW, NPCOL, MYROW, MYCOL,
     $              ILOC1, JLOC1, TRSRC1, TCSRC1 )
               IF( MYROW.EQ.TRSRC1 .AND. MYCOL.EQ.TCSRC1 ) 
     $              WR( K ) = T((JLOC1-1)*LLDT+ILOC1)
            END IF
         ELSE
            PAIR = .FALSE.
         END IF
 560  CONTINUE
      DO 570 K = NB, N-1, NB
         CALL INFOG2L( K, K, DESCT, NPROW, NPCOL, MYROW, MYCOL,
     $        ILOC1, JLOC1, TRSRC1, TCSRC1 )
         CALL INFOG2L( K, K+1, DESCT, NPROW, NPCOL, MYROW, MYCOL,
     $        ILOC2, JLOC2, TRSRC2, TCSRC2 )
         CALL INFOG2L( K+1, K, DESCT, NPROW, NPCOL, MYROW, MYCOL,
     $        ILOC3, JLOC3, TRSRC3, TCSRC3 )
         CALL INFOG2L( K+1, K+1, DESCT, NPROW, NPCOL, MYROW, MYCOL,
     $        ILOC4, JLOC4, TRSRC4, TCSRC4 )
         IF( MYROW.EQ.TRSRC2 .AND. MYCOL.EQ.TCSRC2 ) THEN
            ELEM2 = T((JLOC2-1)*LLDT+ILOC2)
            IF( TRSRC1.NE.TRSRC2 .OR. TCSRC1.NE.TCSRC2 ) 
     $           CALL DGESD2D( ICTXT, 1, 1, ELEM2, 1, TRSRC1, TCSRC1 ) 
         END IF
         IF( MYROW.EQ.TRSRC3 .AND. MYCOL.EQ.TCSRC3 ) THEN
            ELEM3 = T((JLOC3-1)*LLDT+ILOC3)
            IF( TRSRC1.NE.TRSRC3 .OR. TCSRC1.NE.TCSRC3 ) 
     $           CALL DGESD2D( ICTXT, 1, 1, ELEM3, 1, TRSRC1, TCSRC1 )
         END IF
         IF( MYROW.EQ.TRSRC4 .AND. MYCOL.EQ.TCSRC4 ) THEN
            DWORK(1) = T((JLOC4-1)*LLDT+ILOC4)
            IF( K+1.LT.N ) THEN
               DWORK(2) = T((JLOC4-1)*LLDT+ILOC4+1)
            ELSE
               DWORK(2) = ZERO
            END IF
            IF( TRSRC1.NE.TRSRC4 .OR. TCSRC1.NE.TCSRC4 ) 
     $           CALL DGESD2D( ICTXT, 2, 1, DWORK, 2, TRSRC1, 
     $           TCSRC1 )
         END IF
         IF( MYROW.EQ.TRSRC1 .AND. MYCOL.EQ.TCSRC1 ) THEN
            ELEM1 = T((JLOC1-1)*LLDT+ILOC1)
            IF( TRSRC1.NE.TRSRC2 .OR. TCSRC1.NE.TCSRC2 ) 
     $           CALL DGERV2D( ICTXT, 1, 1, ELEM2, 1, TRSRC2, TCSRC2 ) 
            IF( TRSRC1.NE.TRSRC3 .OR. TCSRC1.NE.TCSRC3 ) 
     $           CALL DGERV2D( ICTXT, 1, 1, ELEM3, 1, TRSRC3, TCSRC3 )
            IF( TRSRC1.NE.TRSRC4 .OR. TCSRC1.NE.TCSRC4 ) 
     $           CALL DGERV2D( ICTXT, 2, 1, DWORK, 2, TRSRC4, 
     $           TCSRC4 )
            ELEM4 = DWORK(1)
            ELEM5 = DWORK(2) 
            IF( ELEM5.EQ.ZERO ) THEN
               IF( WR( K ).EQ.ZERO .AND. WI( K ).EQ.ZERO ) THEN
                  CALL DLANV2( ELEM1, ELEM2, ELEM3, ELEM4, WR( K ), 
     $                 WI( K ), WR( K+1 ), WI( K+1 ), SN, CS )
               ELSEIF( WR( K+1 ).EQ.ZERO .AND. WI( K+1 ).EQ.ZERO ) THEN 
                  WR( K+1 ) = ELEM4
               END IF
            ELSEIF( WR( K ).EQ.ZERO .AND. WI( K ).EQ.ZERO ) THEN
               WR( K ) = ELEM1
            END IF
         END IF
 570  CONTINUE
C     
      IF( NPROCS.GT.1 ) THEN
         CALL DGSUM2D( ICTXT, 'All', TOP, N, 1, WR, N, -1, -1 )
         CALL DGSUM2D( ICTXT, 'All', TOP, N, 1, WI, N, -1, -1 )
      END IF
C     
C     Store storage requirements in workspaces
C     
      IF(DBG) WRITE(*,*) MYROW,MYCOL,'Storage requirements'
      DWORK( 1 ) = DBLE(LWMIN)
      IWORK( 1 ) = LIWMIN
C     
C     Return to calling program
C     
      DWORK(2) = TIME1
      DWORK(3) = TIME2
      DWORK(4) = TIME3
      DWORK(5) = TIME4
      IF(DBG) WRITE(*,*) MYROW,MYCOL,'Returning from PBDTRSEN'
C
      RETURN
C     
C     End of PBDTRSEN
C     
      END
C     
