      RECURSIVE SUBROUTINE PDLAQR3( WANTT, WANTZ, N, KTOP, 
     $     KBOT, NW, H, DESCH, ILOZ, IHIZ, Z, DESCZ, NS, ND, 
     $     SR, SI, V, DESCV, NH, T, DESCT, NV, WV, DESCW, 
     $     WORK, LWORK, IWORK, LIWORK, TIMINGS )
*
*  -- ScaLAPACK auxiliary routine (version 1.7.x) --
*     University of Umea and HPC2N, Umea, Sweden
*     February 2008
*
      IMPLICIT NONE
*
*     .. Scalar Arguments ..
      INTEGER            IHIZ, ILOZ, KBOT, KTOP, LWORK, N, ND, NH, NS, 
     $                   NV, NW, LIWORK
      LOGICAL            WANTT, WANTZ
*     ..
*     .. Array Arguments ..
      INTEGER            DESCH( * ), DESCZ( * ), DESCT( * ), DESCV( * ),
     $                   DESCW( * ), IWORK( * )
      DOUBLE PRECISION   H( * ), SI( KBOT ), SR( KBOT ), T( * ),
     $                   V( * ), WORK( * ), WV( * ),
     $                   Z( * ), TIMINGS( 10 )
*     ..
*
*     ******************************************************************
*     Aggressive early deflation:
*
*     This subroutine accepts as input an upper Hessenberg matrix
*     H and performs an orthogonal similarity transformation
*     designed to detect and deflate fully converged eigenvalues from
*     a trailing principal submatrix.  On output H has been over-
*     written by a new Hessenberg matrix that is a perturbation of
*     an orthogonal similarity transformation of H.  It is to be
*     hoped that the final version of H has many zero subdiagonal
*     entries.
*
*     ******************************************************************
*
*  Notes
*  =====
*
*  Each global data object is described by an associated description
*  vector.  This vector stores the information required to establish
*  the mapping between an object element and its corresponding process
*  and memory location.
*
*  Let A be a generic term for any 2D block cyclicly distributed array.
*  Such a global array has an associated description vector DESCA.
*  In the following comments, the character _ should be read as
*  "of the global array".
*
*  NOTATION        STORED IN      EXPLANATION
*  --------------- -------------- --------------------------------------
*  DTYPE_A(global) DESCA( DTYPE_ )The descriptor type.  In this case,
*                                 DTYPE_A = 1.
*  CTXT_A (global) DESCA( CTXT_ ) The BLACS context handle, indicating
*                                 the BLACS process grid A is distribu-
*                                 ted over. The context itself is glo-
*                                 bal, but the handle (the integer
*                                 value) may vary.
*  M_A    (global) DESCA( M_ )    The number of rows in the global
*                                 array A.
*  N_A    (global) DESCA( N_ )    The number of columns in the global
*                                 array A.
*  MB_A   (global) DESCA( MB_ )   The blocking factor used to distribute
*                                 the rows of the array.
*  NB_A   (global) DESCA( NB_ )   The blocking factor used to distribute
*                                 the columns of the array.
*  RSRC_A (global) DESCA( RSRC_ ) The process row over which the first
*                                 row of the array A is distributed.
*  CSRC_A (global) DESCA( CSRC_ ) The process column over which the
*                                 first column of the array A is
*                                 distributed.
*  LLD_A  (local)  DESCA( LLD_ )  The leading dimension of the local
*                                 array.  LLD_A >= MAX(1,LOCr(M_A)).
*
*  Let K be the number of rows or columns of a distributed matrix,
*  and assume that its process grid has dimension p x q.
*  LOCr( K ) denotes the number of elements of K that a process
*  would receive if K were distributed over the p processes of its
*  process column.
*  Similarly, LOCc( K ) denotes the number of elements of K that a
*  process would receive if K were distributed over the q processes of
*  its process row.
*  The values of LOCr() and LOCc() may be determined via a call to the
*  ScaLAPACK tool function, NUMROC:
*          LOCr( M ) = NUMROC( M, MB_A, MYROW, RSRC_A, NPROW ),
*          LOCc( N ) = NUMROC( N, NB_A, MYCOL, CSRC_A, NPCOL ).
*  An upper bound for these quantities may be computed by:
*          LOCr( M ) <= ceil( ceil(M/MB_A)/NPROW )*MB_A
*          LOCc( N ) <= ceil( ceil(N/NB_A)/NPCOL )*NB_A
*
* ======================================================================

*     WANTT   (global input) LOGICAL
*          If .TRUE., then the Hessenberg matrix H is fully updated
*          so that the quasi-triangular Schur factor may be
*          computed (in cooperation with the calling subroutine).
*          If .FALSE., then only enough of H is updated to preserve
*          the eigenvalues.
*
*     WANTZ   (global input) LOGICAL
*          If .TRUE., then the orthogonal matrix Z is updated so
*          so that the orthogonal Schur factor may be computed
*          (in cooperation with the calling subroutine).
*          If .FALSE., then Z is not referenced.
*
*     N       (global input) INTEGER
*          The order of the matrix H and (if WANTZ is .TRUE.) the
*          order of the orthogonal matrix Z.
*
*     KTOP    (global input) INTEGER
*          It is assumed that either KTOP = 1 or H(KTOP,KTOP-1)=0.
*          KBOT and KTOP together determine an isolated block
*          along the diagonal of the Hessenberg matrix.
*
*     KBOT    (global input) INTEGER
*          It is assumed without a check that either
*          KBOT = N or H(KBOT+1,KBOT)=0.  KBOT and KTOP together
*          determine an isolated block along the diagonal of the
*          Hessenberg matrix.
*
*     NW      (global input) INTEGER
*          Deflation window size.  1 .LE. NW .LE. (KBOT-KTOP+1).
*
*     H       (local input/output) DOUBLE PRECISION array, dimension
*             (DESCH(LLD_),*)
*          On input the initial N-by-N section of H stores the
*          Hessenberg matrix undergoing aggressive early deflation.
*          On output H has been transformed by an orthogonal
*          similarity transformation, perturbed, and the returned
*          to Hessenberg form that (it is to be hoped) has some
*          zero subdiagonal entries.
*
*     DESCH   (global and local input) INTEGER array of dimension DLEN_.
*          The array descriptor for the distributed matrix H.
*
*     ILOZ    (global input) INTEGER
*     IHIZ    (global input) INTEGER
*          Specify the rows of Z to which transformations must be
*          applied if WANTZ is .TRUE.. 1 .LE. ILOZ .LE. IHIZ .LE. N.
*
*     Z       (input/output) DOUBLE PRECISION array, dimension 
*             (DESCH(LLD_),*)
*          IF WANTZ is .TRUE., then on output, the orthogonal
*          similarity transformation mentioned above has been
*          accumulated into Z(ILOZ:IHIZ,ILO:IHI) from the right.
*          If WANTZ is .FALSE., then Z is unreferenced.
*
*     DESCZ   (global and local input) INTEGER array of dimension DLEN_.
*          The array descriptor for the distributed matrix Z.
*
*     NS      (global output) INTEGER
*          The number of unconverged (ie approximate) eigenvalues
*          returned in SR and SI that may be used as shifts by the
*          calling subroutine.
*
*     ND      (global output) INTEGER
*          The number of converged eigenvalues uncovered by this
*          subroutine.
*
*     SR      (global output) DOUBLE PRECISION array, dimension KBOT
*     SI      (global output) DOUBLE PRECISION array, dimension KBOT
*          On output, the real and imaginary parts of approximate
*          eigenvalues that may be used for shifts are stored in
*          SR(KBOT-ND-NS+1) through SR(KBOT-ND) and
*          SI(KBOT-ND-NS+1) through SI(KBOT-ND), respectively.
*          The real and imaginary parts of converged eigenvalues
*          are stored in SR(KBOT-ND+1) through SR(KBOT) and
*          SI(KBOT-ND+1) through SI(KBOT), respectively.
*
*     V       (global workspace) DOUBLE PRECISION array, dimension 
*             (DESCV(LLD_),*)
*          An NW-by-NW distributed work array.
*
*     DESCV   (global and local input) INTEGER array of dimension DLEN_.
*          The array descriptor for the distributed matrix V.
*
*     NH      (input) INTEGER scalar
*          The number of columns of T.  NH.GE.NW.
*
*     T       (global workspace) DOUBLE PRECISION array, dimension 
*             (DESCV(LLD_),*)
*
*     DESCT   (global and local input) INTEGER array of dimension DLEN_.
*          The array descriptor for the distributed matrix T.
*
*     NV      (global input) INTEGER
*          The number of rows of work array WV available for
*          workspace.  NV.GE.NW.
*
*     WV      (global workspace) DOUBLE PRECISION array, dimension 
*             (DESCW(LLD_),*)
*
*     DESCW   (global and local input) INTEGER array of dimension DLEN_.
*          The array descriptor for the distributed matrix WV.
*
*     WORK    (local workspace) DOUBLE PRECISION array, dimension LWORK.
*          On exit, WORK(1) is set to an estimate of the optimal value
*          of LWORK for the given values of N, NW, KTOP and KBOT.
*
*     LWORK   (local input) INTEGER
*          The dimension of the work array WORK.  LWORK = 2*NW
*          suffices, but greater efficiency may result from larger
*          values of LWORK.
*
*          If LWORK = -1, then a workspace query is assumed; PDLAQR3
*          only estimates the optimal workspace size for the given
*          values of N, NW, KTOP and KBOT.  The estimate is returned
*          in WORK(1).  No error message related to LWORK is issued
*          by XERBLA.  Neither H nor Z are accessed.
*
*     IWORK   (local workspace) INTEGER array, dimension (ILWORK)
*
*     ILWORK  (local input) INTEGER
*          The length of the workspace array IWORK
*
*     ================================================================
*     Based on contributions by
*        Robert Granat, Department of Computing Science and HPC2N,
*        University of Umea, Sweden.
*
*     ==================================================================
*     .. Parameters ..
      INTEGER            BLOCK_CYCLIC_2D, CSRC_, CTXT_, DLEN_, DTYPE_,
     $                   LLD_, MB_, M_, NB_, N_, RSRC_
      LOGICAL            DEBUG, PRINT, SORTGRAD, RECURSION, AED_ON
      PARAMETER          ( BLOCK_CYCLIC_2D = 1, DLEN_ = 9, DTYPE_ = 1,
     $                     CTXT_ = 2, M_ = 3, N_ = 4, MB_ = 5, NB_ = 6,
     $                     RSRC_ = 7, CSRC_ = 8, LLD_ = 9,
     $                     DEBUG = .FALSE., PRINT = .FALSE.,
     $                     SORTGRAD = .FALSE., RECURSION = .TRUE.,
     $                     AED_ON = .TRUE. )
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0d0, ONE = 1.0d0 )
*     ..
*     .. Local Scalars ..
      DOUBLE PRECISION   AA, BB, BETA, CC, CS, DD, EVI, EVK, FOO, S,
     $                   SAFMAX, SAFMIN, SMLNUM, SN, TAU, ULP, DDUM,
     $                   ELEM, ELEM1, ELEM2, ELEM3, R1, ANORM, RNORM,
     $                   DPDUM, STAMP, RESAED
      INTEGER            I, IFST, ILST, INFO, INFQR, J, JW, K, KCOL,
     $                   KEND, KLN, KROW, KWTOP, LTOP, LWK1, LWK2, LWK3,
     $                   LWKOPT, NMIN, LLDH, LLDZ, LLDT, LLDV, LLDWV,
     $                   ICTXT, NPROW, NPCOL, MYROW, MYCOL, NB, IROFFH,
     $                   M, RCOLS, TAUROWS, RROWS, TAUCOLS, ITAU, IR,
     $                   IPW, IDUM, NPROCS, MLOC, IROFFHH, ICOFFHH,
     $                   HHRSRC, HHCSRC, HHROWS, HHCOLS, IROFFZZ, 
     $                   ICOFFZZ, ZZRSRC, ZZCSRC, ZZROWS, ZZCOLS,
     $                   ITER, IERR, TZROWS0, TZCOLS0, IERR0, IPT0,
     $                   IPZ0, IPW0, NB2, ROUND, LILST, KK, LILST0,
     $                   IWRK1, RSRC, CSRC, LWK4, LWK5, IWRK2, LWK6,
     $                   LWK7, LWK8, ILWKOPT, TZROWS, TZCOLS
      LOGICAL            BULGE, SORTED, LQUERY
*     ..
*     .. Local Arrays ..
      INTEGER            SELECT( N ), PAR( 6 ), DESCR( DLEN_ ),
     $                   DESCTAU( DLEN_ ), DESCHH( DLEN_ ),
     $                   DESCZZ( DLEN_ ), DESCTZ0( DLEN_ )
*     ..
*     .. External Functions ..
      DOUBLE PRECISION   DLAMCH, PDLANGE, MPI_WTIME, PCHKRESI
      INTEGER            PILAENVX, NUMROC, INDXG2P, ICEIL
      EXTERNAL           DLAMCH, PILAENVX, NUMROC, INDXG2P, PDLANGE,
     $                   MPI_WTIME, ICEIL, PCHKRESI
*     ..
*     .. External Subroutines ..
      EXTERNAL           PDCOPY, PDGEHRD, PDGEMM, DLABAD, PDLACPY, 
     $                   PDLAQR1, DLANV2, PDLAQR0, PDLARF, PDLARFG, 
     $                   PDLASET, PDORGHR, PBDTRSEN, PDELGET, PDELSET,
     $                   PDLAMVE
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, DBLE, INT, MAX, MIN, SQRT
*     ..
*     .. Executable Statements ..
      IF( DEBUG ) WRITE(*,*) MYROW,MYCOL, 
     $     'PDLAQR3 N,NB,KTOP,KBOT=',N,DESCH(MB_),KTOP,KBOT
      ICTXT = DESCH( CTXT_ )
      CALL BLACS_GRIDINFO( ICTXT, NPROW, NPCOL, MYROW, MYCOL )
      NPROCS = NPROW*NPCOL
*
*     Extract local leading dimensions, blockfactors, offset for
*     keeping the alignment requirements and size of defl. window
*
      LLDH  = DESCH( LLD_ )
      LLDZ  = DESCZ( LLD_ )
      LLDT  = DESCT( LLD_ )
      LLDV  = DESCV( LLD_ )
      LLDWV = DESCW( LLD_ )
      NB = DESCH( MB_ )
      IROFFH = MOD( KTOP - 1, NB )
      ITER = IWORK(1)
      IF( DEBUG ) WRITE(*,*) MYROW,MYCOL,'ITER=',ITER
      JW = MIN( NW, KBOT-KTOP+1 )
      IF( DEBUG ) WRITE(*,*) MYROW,MYCOL, 'JW=',JW
*
*     Extract environment variables for parallel eigenvalue reordering
*
      PAR(1) = PILAENVX(ICTXT, 17, 'PDLAQR3', 'SV', JW, NB, IDUM, IDUM)
      PAR(2) = PILAENVX(ICTXT, 18, 'PDLAQR3', 'SV', JW, NB, IDUM, IDUM)
      PAR(3) = PILAENVX(ICTXT, 19, 'PDLAQR3', 'SV', JW, NB, IDUM, IDUM)
      PAR(4) = PILAENVX(ICTXT, 20, 'PDLAQR3', 'SV', JW, NB, IDUM, IDUM)
      PAR(5) = PILAENVX(ICTXT, 21, 'PDLAQR3', 'SV', JW, NB, IDUM, IDUM)
      PAR(6) = PILAENVX(ICTXT, 22, 'PDLAQR3', 'SV', JW, NB, IDUM, IDUM)
*
*     Check if workspace query
*
      LQUERY = LWORK.EQ.-1 .OR. LIWORK.EQ.-1
*
*     ==== Estimate optimal workspace. ====
*
      IF( JW.LE.2 ) THEN
         LWKOPT = 1
      ELSE
*
*        ==== Workspace query calls to PDGEHRD and PDORMHR ====
*
         TAUROWS = NUMROC( 1, 1, MYCOL, DESCV(RSRC_), NPROW )
         TAUCOLS = NUMROC( JW+IROFFH, NB, MYCOL, DESCV(CSRC_), 
     $        NPCOL )
         CALL PDGEHRD( JW, 1, JW, T, 1, 1, DESCT, WORK, WORK, 
     $                 -1, INFO )
         LWK1 = INT( WORK( 1 ) ) + TAUROWS*TAUCOLS
*
*        ==== Workspace query call to PDORMHR ====
*
         CALL PDORMHR( 'Right', 'No', JW, JW, 1, JW, T, 1, 1, 
     $                 DESCT, WORK, V, 1, 1, DESCV, WORK, -1, INFO )
         LWK2 = INT( WORK( 1 ) )
*
*        ==== Workspace query call to PDLAQR0 ====
*
         NMIN = PILAENVX( ICTXT, 12, 'PDLAQR3', 'SV', JW, 1, JW, LWORK )
         IF( JW+IROFFH.GT.NMIN .AND. RECURSION ) THEN
            IF( DEBUG ) WRITE(*,*) MYROW,MYCOL, 
     $           'Workspace query call to PDLAQR0'
            CALL PDLAQR0( .TRUE., .TRUE., JW+IROFFH, 1+IROFFH, 
     $           JW+IROFFH, T, DESCT, SR, SI, 1, JW, V, DESCV, WORK, 
     $           -1, IWORK, LIWORK, INFQR, TIMINGS )
            LWK3 = INT( WORK( 1 ) )
            IWRK1 = IWORK(1)
         ELSE
            RSRC = DESCT( RSRC_ )
            CSRC = DESCT( CSRC_ )
            DESCT( RSRC_ ) = 0
            DESCT( CSRC_ ) = 0
            CALL PDLAQR1( .TRUE., .TRUE., JW+IROFFH, 1, JW+IROFFH, 
     $           T, DESCT, SR, SI, 1, JW+IROFFH, V, DESCV, WORK, 
     $           -1, IWORK, LIWORK, INFQR )
            DESCT( RSRC_ ) = RSRC
            DESCT( CSRC_ ) = CSRC
            LWK3 = INT( WORK( 1 ) )
            IWRK1 = IWORK(1)
         END IF
*
*        ==== Workspace in case of alignment problems ====
*
         TZROWS0 = NUMROC( JW+IROFFH, NB, MYROW, 0, NPROW )
         TZCOLS0 = NUMROC( JW+IROFFH, NB, MYCOL, 0, NPCOL )
         LWK4 = 2 * TZROWS0*TZCOLS0
*
*        ==== Workspace check for reordering ====
*
         CALL PBDTRSEN( 'No condition numbers', 'Vectors', 
     $        SELECT, PAR, JW+IROFFH, T, 1, 1, DESCT, V, 1, 1, 
     $        DESCV, DDUM, DDUM, MLOC, DDUM, DDUM, WORK, -1,
     $        IWORK, LIWORK, INFO )
         LWK5 = INT( WORK( 1 ) )
         IWRK2 = IWORK(1)
*
*        ==== Extra workspace for reflecting back spike ====
*        (workspace for PDLARF approximated for simplicity)
*
         RROWS =  NUMROC( NS+IROFFH, NB, MYROW, DESCV(RSRC_), NPROW )
         RCOLS =  NUMROC( 1, 1, MYCOL, DESCV(CSRC_), NPCOL )
         LWK6 = RROWS*RCOLS + TAUROWS*TAUCOLS +
     $        2*ICEIL(ICEIL(JW+IROFFH,NB),NPROW)*NB
     $         *ICEIL(ICEIL(JW+IROFFH,NB),NPCOL)*NB
*
*        ==== Extra workspace needed by PBLAS update calls
*        (also estimated for simplicity) ====
*
         LWK7 = MAX( ICEIL(ICEIL(JW,NB),NPROW)*NB *
     $               ICEIL(ICEIL(N-KBOT,NB),NPCOL)*NB,
     $               ICEIL(ICEIL(IHIZ-ILOZ+1,NB),NPROW)*NB *
     $               ICEIL(ICEIL(JW,NB),NPCOL)*NB,
     $               ICEIL(ICEIL(KBOT-JW,NB),NPROW)*NB *
     $               ICEIL(ICEIL(JW,NB),NPCOL)*NB ) 
*
*        === Residual check workspace ===
*
#ifdef USE_AED_RES
         TZROWS = NUMROC( JW+IROFFH, NB, MYROW, DESCT(RSRC_), NPROW )
         TZCOLS = NUMROC( JW+IROFFH, NB, MYCOL, DESCT(CSRC_), NPCOL )
         LWK8 = 2*TZROWS*TZCOLS
#else
         LWK8 = 0
#endif
*
*        ==== Optimal workspace ====
*
         LWKOPT = MAX( LWK1, LWK2, LWK3+LWK4, LWK5, LWK6, LWK7, LWK8 )
         ILWKOPT = MAX( IWRK1, IWRK2 )
      END IF
*
*     ==== Quick return in case of workspace query. ====
*
      IF( LQUERY ) THEN
         WORK( 1 ) = DBLE( LWKOPT )
         IWORK( 1 ) = ILWKOPT
         RETURN
      END IF
*
*     ==== Nothing to do ...
*     ... for an empty active block ... ====
      NS = 0
      ND = 0
      IF( KTOP.GT.KBOT )
     $   RETURN
*     ... nor for an empty deflation window. ====
      IF( NW.LT.1 )
     $   RETURN
*
*     ==== Machine constants ====
*
      SAFMIN = DLAMCH( 'SAFE MINIMUM' )
      SAFMAX = ONE / SAFMIN
      CALL DLABAD( SAFMIN, SAFMAX )
      ULP = DLAMCH( 'PRECISION' )
      SMLNUM = SAFMIN*( DBLE( N ) / ULP )
*
*     ==== Setup deflation window ====
*
      JW = MIN( NW, KBOT-KTOP+1 )
      KWTOP = KBOT - JW + 1
      IF( DEBUG) WRITE(*,*) MYROW,MYCOL, 'KWTOP=',KWTOP
      IF( KWTOP.EQ.KTOP ) THEN
         S = ZERO
      ELSE
         CALL PDELGET( 'All', '1-Tree', S, H, KWTOP, KWTOP-1, DESCH )
      END IF
*
      IF( KBOT.EQ.KWTOP ) THEN
*
*        ==== 1-by-1 deflation window: not much to do ====
*
         CALL PDELGET( 'All', '1-Tree', SR( KWTOP ), H, KWTOP, KWTOP, 
     $                 DESCH )
         SI( KWTOP ) = ZERO
         NS = 1
         ND = 0
         IF( ABS( S ).LE.MAX( SMLNUM, ULP*ABS( SR( KWTOP ) ) ) )
     $        THEN
            NS = 0
            ND = 1
            IF( KWTOP.GT.KTOP )
     $      CALL PDELSET( H, KWTOP, KWTOP-1 , DESCH, ZERO )
         END IF
         RETURN
      END IF
*
      IF( KWTOP.EQ.KTOP .AND. KBOT-KWTOP.EQ.1 ) THEN
*
*     ==== 2-by-2 deflation window: a little more to do ====
*
         CALL PDELGET( 'All', '1-Tree', AA, H, KWTOP, KWTOP, 
     $                 DESCH )
         CALL PDELGET( 'All', '1-Tree', BB, H, KWTOP, KWTOP+1, 
     $                 DESCH )
         CALL PDELGET( 'All', '1-Tree', CC, H, KWTOP+1, KWTOP, 
     $                 DESCH )
         CALL PDELGET( 'All', '1-Tree', DD, H, KWTOP+1, KWTOP+1, 
     $                 DESCH )
         CALL DLANV2( AA, BB, CC, DD, SR(KWTOP), SI(KWTOP), 
     $                SR(KWTOP+1), SI(KWTOP+1), CS, SN )
         NS = 0
         ND = 2
         IF( CC.EQ.ZERO ) THEN
            I = KWTOP
            IF( I+2.LE.N )
     $           CALL PDROT( N-I-1, H, I, I+2, DESCH,
     $           DESCH(M_), H, I+1, I+2, DESCH, DESCH(M_), 
     $           CS, SN, WORK, LWORK, INFO )
            CALL PDROT( I-1, H, 1, I, DESCH, 1,
     $           H, 1, I+1, DESCH, 1, CS, SN,
     $           WORK, LWORK, INFO )
            CALL PDROT( N, Z, 1, I, DESCZ, 1,
     $           Z, 1, I+1, DESCZ, 1, CS, SN,
     $           WORK, LWORK, INFO )
            CALL PDELSET( H, I, I, DESCH, AA )
            CALL PDELSET( H, I, I+1, DESCH, BB )
            CALL PDELSET( H, I+1, I, DESCH, CC )
            CALL PDELSET( H, I+1, I+1, DESCH, DD )
         END IF
         RETURN
      END IF
*
*     Calculate new value for IROFFH in case deflation window
*     was adjusted
*
      IROFFH = MOD( KWTOP - 1, NB ) 
      IF( DEBUG ) WRITE(*,*) MYROW,MYCOL, 'IROFFH=',IROFFH
*
*     Adjust number of rows and columns of T matrix descriptor
*     to prepare for call to PDBTRSEN
*
      DESCT( M_ ) = JW+IROFFH
      DESCT( N_ ) = JW+IROFFH
*
*     ==== Convert to spike-triangular form.  (In case of a
*     .    rare QR failure, this routine continues to do
*     .    aggressive early deflation using that part of
*     .    the deflation window that converged using INFQR
*     .    here and there to keep track.) ====
*
      CALL PDLASET( 'All', IROFFH, JW+IROFFH, ZERO, ONE, T, 1, 1, 
     $              DESCT )
      CALL PDLASET( 'All', IROFFH, JW-IROFFH, ZERO, ZERO, T, 1, 
     $              1+IROFFH, DESCT ) 
      CALL PDLACPY( 'Upper', JW-1, JW-1, H, KWTOP+1, KWTOP, DESCH, T, 
     $              1+IROFFH+1, 1+IROFFH, DESCT )
      CALL PDLACPY( 'All', 1, JW, H, KWTOP, KWTOP, DESCH, T, 1+IROFFH, 
     $              1+IROFFH, DESCT )
      CALL PDLACPY( 'All', JW-1, 1, H, KWTOP+1, KWTOP+JW-1, DESCH, T, 
     $              1+IROFFH+1, 1+IROFFH+JW-1, DESCT )
      IF( JW.GT.2 )
     $     CALL PDLASET( 'Lower', JW-2, JW-2, ZERO, ZERO, T, 1+IROFFH+2, 
     $                   1+IROFFH, DESCT )
      CALL PDLASET( 'All', JW, IROFFH, ZERO, ZERO, T, 1+IROFFH, 1, 
     $              DESCT ) 
*
      CALL PDLASET( 'All', JW+IROFFH, JW+IROFFH, ZERO, ONE, V, 1, 
     $              1, DESCV )
      IF( PRINT.AND.ITER.EQ.3 ) THEN
         CALL PDLAPRNT( JW+IROFFH, JW+IROFFH, T, 1, 1, DESCT, 0, 0, 
     $        'H3', 6, WORK )
         CALL PDLAPRNT( JW+IROFFH, JW+IROFFH, V, 1, 1, DESCV, 0, 0, 
     $        'Q3', 6, WORK )
c         IF( DEBUG ) WRITE(*,*) MYROW,MYCOL, 'SR=',SR(1:KBOT)
c         IF( DEBUG ) WRITE(*,*) MYROW,MYCOL, 'SI=',SI(1:KBOT)
      END IF
      NMIN = PILAENVX( ICTXT, 12, 'PDLAQR3', 'SV', JW, 1, JW, LWORK )
      IF( JW+IROFFH.GT.NMIN .AND. RECURSION ) THEN
         IF( DEBUG ) WRITE(*,*) MYROW,MYCOL, 'PDLAQR3 -> PDLAQR0'
         IWORK(1) = ITER
         STAMP = MPI_WTIME()
         CALL PDLAQR0( .TRUE., .TRUE., JW+IROFFH, 1+IROFFH, JW+IROFFH, 
     $                 T, DESCT, SR( KWTOP-IROFFH ), SI( KWTOP-IROFFH ), 
     $                 1+IROFFH, JW+IROFFH, V, DESCV, WORK, LWORK, 
     $                 IWORK, LIWORK, INFQR, TIMINGS )
         TIMINGS( 4 ) = TIMINGS( 4 ) + MPI_WTIME() - STAMP
         TIMINGS( 3 ) = TIMINGS( 3 ) - MPI_WTIME() + STAMP
         IF( DEBUG ) WRITE(*,*) MYROW,MYCOL, 'PDLAQR3 <- PDLAQR0'    
      ELSE
         IF( DEBUG ) WRITE(*,*) MYROW,MYCOL,'MB_,NB_ 0=',DESCT(MB_),
     $        DESCT(NB_)
         STAMP = MPI_WTIME()
         IF( DESCT(RSRC_).EQ.0 .AND. DESCT(CSRC_).EQ.0 ) THEN
            IF( DEBUG ) WRITE(*,*) MYROW,MYCOL,'MB_,NB_ 1=',DESCT(MB_),
     $           DESCT(NB_)
            IF( JW+IROFFH.GT.DESCT( MB_ ) ) THEN
               IF( DEBUG ) WRITE(*,*) MYROW,MYCOL, 'PDLAQR3 -> PDLAQR1'
               CALL PDLAQR1( .TRUE., .TRUE., JW+IROFFH, 1, 
     $                       JW+IROFFH, T, DESCT, SR( KWTOP-IROFFH ), 
     $                       SI( KWTOP-IROFFH ), 1, JW+IROFFH, V, 
     $                       DESCV, WORK, LWORK, IWORK, LIWORK, INFQR )
            ELSE
               IF( MYROW.EQ.0 .AND. MYCOL.EQ.0 ) THEN
                  IF( DEBUG ) WRITE(*,*) MYROW,MYCOL, 
     $                 'PDLAQR3 -> DLAHQR 1'
                  CALL DLAHQR( .TRUE., .TRUE., JW+IROFFH, 1+IROFFH, 
     $                 JW+IROFFH, T, DESCT(LLD_), 
     $                 SR( KWTOP-IROFFH ), SI( KWTOP-IROFFH ), 
     $                 1+IROFFH, JW+IROFFH, V, DESCV(LLD_), 
     $                 INFQR )
               ELSE
                  INFQR = 0
               END IF
               IF( NPROCS.GT.1 )
     $              CALL IGAMN2D( ICTXT, 'All', '1-Tree', 1, 1, INFQR, 
     $              1, -1, -1, -1, -1, -1 )
            END IF
         ELSEIF( JW+IROFFH.LE.DESCT( MB_ ) ) THEN
            IF( MYROW.EQ.DESCT(RSRC_) .AND. MYCOL.EQ.DESCT(CSRC_) ) THEN
               IF( DEBUG ) WRITE(*,*) MYROW,MYCOL, 'PDLAQR3 -> DLAHQR 2'
               CALL DLAHQR( .TRUE., .TRUE., JW+IROFFH, 1+IROFFH, 
     $                      JW+IROFFH, T, DESCT(LLD_), 
     $                      SR( KWTOP-IROFFH ), SI( KWTOP-IROFFH ), 
     $                      1+IROFFH, JW+IROFFH, V, DESCV(LLD_), 
     $                      INFQR )
            ELSE
               INFQR = 0
            END IF
            IF( NPROCS.GT.1 )
     $           CALL IGAMN2D( ICTXT, 'All', '1-Tree', 1, 1, INFQR, 
     $           1, -1, -1, -1, -1, -1 )
         ELSE
            TZROWS0 = NUMROC( JW+IROFFH, NB, MYROW, 0, NPROW )
            TZCOLS0 = NUMROC( JW+IROFFH, NB, MYCOL, 0, NPCOL )
            CALL DESCINIT( DESCTZ0, JW+IROFFH, JW+IROFFH, NB, NB, 0, 
     $                     0, ICTXT, MAX(1,TZROWS0), IERR0 )
            IPT0 = 1
            IPZ0 = IPT0 + TZROWS0*TZCOLS0
            IPW0 = IPZ0 + TZROWS0*TZCOLS0
            CALL PDLAMVE( 'All', JW+IROFFH, JW+IROFFH, T, 1, 1, DESCT, 
     $           WORK(IPT0), 1, 1, DESCTZ0, WORK(IPW0) )
            CALL PDLASET( 'All', JW+IROFFH, JW+IROFFH, ZERO, ONE, 
     $                    WORK(IPZ0), 1, 1, DESCTZ0 )
            IF( DEBUG ) WRITE(*,*) MYROW,MYCOL,'MB_,NB_ 2=',
     $           DESCTZ0(MB_), DESCTZ0(NB_)
            IF( DEBUG ) WRITE(*,*) MYROW,MYCOL, 'PDLAQR3 -> PDLAQR1'
            CALL PDLAQR1( .TRUE., .TRUE., JW+IROFFH, 1, 
     $                    JW+IROFFH, WORK(IPT0), DESCTZ0, 
     $                    SR( KWTOP-IROFFH ), SI( KWTOP-IROFFH ), 
     $                    1, JW+IROFFH, WORK(IPZ0), 
     $                    DESCTZ0, WORK(IPW0), LWORK-IPW0+1, IWORK, 
     $                    LIWORK, INFQR )
            CALL PDLAMVE( 'All', JW+IROFFH, JW+IROFFH, WORK(IPT0), 1, 1, 
     $                    DESCTZ0, T, 1, 1, DESCT, WORK(IPW0) )
            CALL PDLAMVE( 'All', JW+IROFFH, JW+IROFFH, WORK(IPZ0), 1, 1, 
     $                    DESCTZ0, V, 1, 1, DESCV, WORK(IPW0) )
         END IF
         TIMINGS( 1 ) = TIMINGS( 1 ) + MPI_WTIME() - STAMP
         TIMINGS( 3 ) = TIMINGS( 3 ) - MPI_WTIME() + STAMP
         IF( DEBUG ) WRITE(*,*) MYROW,MYCOL, 'PDLAQR3 <- PDLAQR1' 
      END IF
*
*     === Debug output ===
*
      IF( PRINT.AND.ITER.EQ.3 ) THEN
         CALL PDLAPRNT( JW+IROFFH, JW+IROFFH, T, 1, 1, DESCT, 0, 0, 
     $        'S3', 6, WORK )
         CALL PDLAPRNT( JW+IROFFH, JW+IROFFH, V, 1, 1, DESCV, 0, 0, 
     $        'Z3', 6, WORK ) 
c         IF( DEBUG ) WRITE(*,*) MYROW,MYCOL, 'SR=',SR(1:KBOT)
c         IF( DEBUG ) WRITE(*,*) MYROW,MYCOL, 'SI=',SI(1:KBOT)
c         RETURN
      END IF
*
*     Adjust INFQR for offset from block border in submatrices
*
      IF( DEBUG ) WRITE(*,*) MYROW,MYCOL, 'Adjust infqr'
      IF( INFQR.NE.0 )
     $     INFQR = INFQR - IROFFH
*
*     In case of no AED, set ND and NS and return
*
      IF( .NOT. AED_ON ) THEN
         NS = JW
         ND = 0
         RETURN
      END IF
*
*     ==== PBDTRSEN needs a clean margin near the diagonal ====
*
      IF( DEBUG ) WRITE(*,*) MYROW,MYCOL, 'Clear margin'
      DO 10 J = 1, JW - 3
         CALL PDELSET( T, J+2, J, DESCT, ZERO )
         CALL PDELSET( T, J+3, J, DESCT, ZERO )
   10 CONTINUE
      IF( JW.GT.2 )
     $   CALL PDELSET( T, JW, JW-2, DESCT, ZERO )
      IF( DEBUG ) WRITE(*,*) MYROW,MYCOL, 'Clear margin done!'
*
*     === Check local residual for AED Schur decomposition ===
*
#ifdef USE_AED_RES
      RESAED = PCHKRESI( JW, H, KWTOP, KWTOP, DESCH, T, 
     $                   1+IROFFH, 1+IROFFH, DESCT, V, 1+IROFFH, 
     $                   1+IROFFH, DESCV, WORK, LWORK )
      TIMINGS(10) = MAX(TIMINGS(10),RESAED)
      IF( DEBUG ) WRITE(*,*) MYROW,MYCOL,'RESAED=',RESAED
#else
      RESAED = 0.0D+00
      TIMINGS(10) = MAX(TIMINGS(10),RESAED)
#endif
*
*     ==== Clean up SELECT ====
*
       IF( DEBUG ) WRITE(*,*) MYROW,MYCOL, 'Clean select'
      DO 13 J = 1, N
         SELECT( J ) = 0
 13   CONTINUE
*
*     ==== Set local M counter to zero  ====
*
      MLOC = 0
*
*     ==== Outer deflation detection loop ====
*
      IF( DEBUG ) WRITE(*,*) MYROW,MYCOL, 'Defl. det. loop'
      DO 15 J = 1, IROFFH + INFQR
         SELECT( J ) = 1
 15   CONTINUE
 9992 FORMAT('Outer 0 select:',100(1X,I3))
      IF( DEBUG ) WRITE(*,FMT=9992) 
     $     (SELECT(KK), KK=1,JW+IROFFH)
*
      NS = JW
      ILST = INFQR + 1 + IROFFH
      IF( ILST.GT.1 ) THEN
         CALL PDELGET( 'All', '1-Tree', ELEM, T, ILST, ILST-1, 
     $        DESCT )
         BULGE = ELEM.NE.ZERO
         IF( BULGE ) ILST = ILST+1
      END IF
*
 18   CONTINUE
      IF( ILST.LE.NS+IROFFH ) THEN
         LILST = MAX(ILST,NS+IROFFH-NB+1)
         IF( LILST.GT.1 ) THEN
            CALL PDELGET( 'All', '1-Tree', ELEM, T, LILST, LILST-1, 
     $           DESCT )
            BULGE = ELEM.NE.ZERO
            IF( BULGE ) LILST = LILST+1
         END IF
         DO 19 J = IROFFH+1, LILST-1
            SELECT( J ) = 1
 19      CONTINUE
         LILST0 = LILST
*
 9993    FORMAT('Outer 1 select:',100(1X,I3))
         IF( DEBUG ) WRITE(*,FMT=9993) 
     $        (SELECT(KK), KK=1,JW+IROFFH)
*
*     ==== Inner deflation detection loop ====
*
 20      CONTINUE
         IF( LILST.LE.NS+IROFFH ) THEN
            IF( NS.EQ.1 ) THEN
               BULGE = .FALSE.
            ELSE
               CALL PDELGET( 'All', '1-Tree', ELEM, T, NS+IROFFH, 
     $              NS+IROFFH-1, DESCT )
               BULGE = ELEM.NE.ZERO
            END IF
*     
*     ==== Small spike tip test for deflation ====
*     
            IF( .NOT.BULGE ) THEN
*     
*     ==== Real eigenvalue ====
*     
               CALL PDELGET( 'All', '1-Tree', ELEM, T, NS+IROFFH, 
     $              NS+IROFFH, DESCT )
               FOO = ABS( ELEM )
               IF( FOO.EQ.ZERO )
     $              FOO = ABS( S )
               CALL PDELGET( 'All', '1-Tree', ELEM, V, 1+IROFFH, 
     $              NS+IROFFH, DESCV )
               IF( ABS( S*ELEM ).LE.MAX( SMLNUM, ULP*FOO ) ) THEN
*     
*     ==== Deflatable ====
*     
                  NS = NS - 1
               ELSE
*     
*     ==== Undeflatable. Move it up out of the way. ====
*     
                  IFST = NS
                  DO 22 J = LILST, JW+IROFFH
                     SELECT( J ) = 0
 22               CONTINUE
                  SELECT( IFST+IROFFH ) = 1
 9994             FORMAT('Inner 1 select:',100(1X,I3))
                  IF( DEBUG ) WRITE(*,FMT=9994) 
     $                 (SELECT(KK), KK=1,JW+IROFFH)
                  STAMP = MPI_WTIME()
                  CALL PBDTRSEN( 'No condition numbers', 'Vectors', 
     $                 SELECT, PAR, JW+IROFFH, T, 1, 1, DESCT, V, 1, 1, 
     $                 DESCV, WORK, WORK(JW+IROFFH+1), MLOC, DDUM, DDUM, 
     $                 WORK(2*(JW+IROFFH)+1), LWORK-2*(JW+IROFFH), 
     $                 IWORK, LIWORK, INFO )
                  TIMINGS( 7 ) = TIMINGS( 7 ) + MPI_WTIME() - STAMP
                  TIMINGS( 3 ) = TIMINGS( 3 ) - MPI_WTIME() + STAMP
c                  SELECT( LILST ) = 1
 9995             FORMAT('Inner 1 update:',100(1X,I3))
                  IF( DEBUG ) WRITE(*,FMT=9995) 
     $                 (SELECT(KK), KK=1,JW+IROFFH)
                  LILST = LILST + 1
               END IF
            ELSE
*     
*     ==== Complex conjugate pair ====
*     
               CALL PDELGET( 'All', '1-Tree', ELEM1, T, NS+IROFFH, 
     $              NS+IROFFH, DESCT )
               CALL PDELGET( 'All', '1-Tree', ELEM2, T, NS+IROFFH, 
     $              NS+IROFFH-1, DESCT )
               CALL PDELGET( 'All', '1-Tree', ELEM3, T, NS+IROFFH-1, 
     $              NS+IROFFH, DESCT )
               FOO = ABS( ELEM1 ) + SQRT( ABS( ELEM2 ) )*
     $              SQRT( ABS( ELEM3 ) )
               IF( FOO.EQ.ZERO )
     $              FOO = ABS( S )
               CALL PDELGET( 'All', '1-Tree', ELEM1, V, 1+IROFFH, 
     $              NS+IROFFH, DESCV )
               CALL PDELGET( 'All', '1-Tree', ELEM2, V, 1+IROFFH, 
     $              NS+IROFFH-1, DESCV )
               IF( MAX( ABS( S*ELEM1 ), ABS( S*ELEM2 ) ).LE.
     $              MAX( SMLNUM, ULP*FOO ) ) THEN
*     
*     ==== Deflatable ====
*     
                  NS = NS - 2
               ELSE
*     
*     ==== Undeflatable. Move them up out of the way. ====
*     
                  IFST = NS
                  DO 24 J = LILST, JW+IROFFH
                     SELECT( J ) = 0
 24               CONTINUE
                  SELECT( IFST+IROFFH ) = 1
                  SELECT( IFST+IROFFH-1 ) = 1
 9996             FORMAT('Inner 2 select:',100(1X,I3))
                  IF( DEBUG ) WRITE(*,FMT=9996) 
     $                 (SELECT(KK), KK=1,JW+IROFFH)
                  STAMP = MPI_WTIME()
                  CALL PBDTRSEN( 'No condition numbers', 'Vectors', 
     $                 SELECT, PAR, JW+IROFFH, T, 1, 1, DESCT, V, 
     $                 1, 1, DESCV, WORK, WORK(JW+IROFFH+1), MLOC, DDUM, 
     $                 DDUM, WORK(2*(JW+IROFFH)+1), LWORK-2*(JW+IROFFH), 
     $                 IWORK, LIWORK, INFO )
                  TIMINGS( 7 ) = TIMINGS( 7 ) + MPI_WTIME() - STAMP
                  TIMINGS( 3 ) = TIMINGS( 3 ) - MPI_WTIME() + STAMP
c                  SELECT( LILST ) = 1
c                  SELECT( LILST+1 ) = 1
                  IF( DEBUG ) WRITE(*,FMT=9997) 
     $                 (SELECT(KK), KK=1,JW+IROFFH)
 9997             FORMAT('Inner 2 update:',100(1X,I3))                 
                  LILST = LILST + 2
               END IF
            END IF
*     
*     ==== End inner deflation detection loop ====
*     
            GO TO 20
         END IF
         DO 21 J = ILST, LILST0-1
            SELECT( J ) = 0
 21      CONTINUE
 9998    FORMAT('Outer 2 select:',100(1X,I3))
         IF( DEBUG ) WRITE(*,FMT=9998) 
     $        (SELECT(KK), KK=1,JW+IROFFH)
         STAMP = MPI_WTIME()
         CALL PBDTRSEN( 'No condition numbers', 'Vectors', 
     $        SELECT, PAR, JW+IROFFH, T, 1, 1, DESCT, V, 1, 1, 
     $        DESCV, WORK, WORK(JW+IROFFH+1), M, DDUM, DDUM, 
     $        WORK(2*(JW+IROFFH)+1), LWORK-2*(JW+IROFFH), 
     $        IWORK, LIWORK, INFO )
 9999    FORMAT('Outer 2 update:',100(1X,I3)) 
         IF( DEBUG ) WRITE(*,FMT=9999) 
     $        (SELECT(KK), KK=1,JW+IROFFH)
         TIMINGS( 7 ) = TIMINGS( 7 ) + MPI_WTIME() - STAMP
         TIMINGS( 3 ) = TIMINGS( 3 ) - MPI_WTIME() + STAMP
         ILST = M + 1 
*     
*     ==== End outer deflation detection loop ====
*  
         GO TO 18
      END IF

      IF( DEBUG ) WRITE(*,*) MYROW,MYCOL, 'NS=',NS
C      IF( DEBUG ) WRITE(*,*) MYROW,MYCOL, 'SELECT=',
c     $     SELECT(1:JW+IROFFH)
*
*     Post-reordering step: copy output eigenvalues to output
*
      CALL DCOPY( JW, WORK(1+IROFFH), 1, SR( KWTOP ), 1 )
      CALL DCOPY( JW, WORK(JW+2*IROFFH+1), 1, SI( KWTOP ), 1 )
      IF( PRINT.AND.ITER.EQ.3 ) THEN
         CALL PDLAPRNT( JW+IROFFH, JW+IROFFH, T, 1, 1, DESCT, 0, 0, 
     $        'S4', 6, WORK )
         CALL PDLAPRNT( JW+IROFFH, JW+IROFFH, V, 1, 1, DESCV, 0, 0, 
     $        'Z4', 6, WORK ) 
      END IF
c      IF( DEBUG ) WRITE(*,*) MYROW,MYCOL, 'SR=',SR(KWTOP:KWTOP+JW-1)
c      IF( DEBUG ) WRITE(*,*) MYROW,MYCOL, 'SI=',SI(KWTOP:KWTOP+JW-1)
*
*     === Check local residual for reordered AED Schur 
*     ... decomposition ===
*
#ifdef USE_AED_RES
      RESAED = PCHKRESI( JW, H, KWTOP, KWTOP, DESCH, T, 
     $                   1+IROFFH, 1+IROFFH, DESCT, V, 1+IROFFH, 
     $                   1+IROFFH, DESCV, WORK, LWORK )
      TIMINGS(10) = MAX(TIMINGS(10),RESAED)
      IF( DEBUG ) WRITE(*,*) MYROW,MYCOL,'RESAED=',RESAED
#else
      RESAED = 0.0D+00
      TIMINGS(10) = MAX(TIMINGS(10),RESAED)
#endif
*
*        ==== Return to Hessenberg form ====
*
      IF( NS.EQ.0 )
     $   S = ZERO
*
      IF( NS.LT.JW .AND. SORTGRAD ) THEN
*
*        ==== sorting diagonal blocks of T improves accuracy for
*        .    graded matrices.  Bubble sort deals well with
*        .    exchange failures. Eigenvalues/shifts from T are also
*        .    restored.  ====
*
         IF( DEBUG ) WRITE(*,*) MYROW,MYCOL, 
     $        'sorting graded matrices...'
         ROUND = 0
         SORTED = .FALSE.
         I = NS + 1 + IROFFH
   30    CONTINUE
         IF( SORTED )
     $      GO TO 50
         SORTED = .TRUE.
         ROUND = ROUND + 1
         IF( DEBUG ) WRITE(*,*) MYROW,MYCOL, 'graded round=',ROUND
         IF( PRINT.AND.ITER.EQ.3 ) THEN
            CALL PDLAPRNT( JW+IROFFH, JW+IROFFH, T, 1, 1, DESCT, 0, 0, 
     $           'SG', 6, WORK )
            CALL PDLAPRNT( JW+IROFFH, JW+IROFFH, V, 1, 1, DESCV, 0, 0, 
     $           'ZG', 6, WORK ) 
c            WRITE(*,*) MYROW,MYCOL, 'SR=',SR(KWTOP:KWTOP+JW-1)
c            WRITE(*,*) MYROW,MYCOL, 'SI=',SI(KWTOP:KWTOP+JW-1)
         END IF
*
         KEND = I - 1
         I = INFQR + 1 + IROFFH
         IF( I.EQ.NS+IROFFH ) THEN
            K = I + 1
         ELSE IF( SI( KWTOP-IROFFH + I-1 ).EQ.ZERO ) THEN
            K = I + 1
         ELSE
            K = I + 2
         END IF    
   40    CONTINUE
         IF( DEBUG ) WRITE(*,*) MYROW,MYCOL,'swap K,I=',K,I
         IF( K.LE.KEND ) THEN
            IF( K.EQ.I+1 ) THEN
               EVI = ABS( SR( KWTOP-IROFFH+I-1 ) )
               IF( DEBUG ) WRITE(*,*) MYROW,MYCOL,'EVI 1'
            ELSE
               EVI = ABS( SR( KWTOP-IROFFH+I-1 ) ) + 
     $              ABS( SI( KWTOP-IROFFH+I-1 ) ) 
               IF( DEBUG ) WRITE(*,*) MYROW,MYCOL,'EVI 2'
            END IF
            IF( DEBUG ) WRITE(*,*) MYROW,MYCOL,'EVI=',EVI
*
            IF( K.EQ.KEND ) THEN
               EVK = ABS( SR( KWTOP-IROFFH+K-1 ) )
               IF( DEBUG ) WRITE(*,*) MYROW,MYCOL,'EVK 1'
            ELSEIF( SI( KWTOP-IROFFH+K-1 ).EQ.ZERO ) THEN
               EVK = ABS( SR( KWTOP-IROFFH+K-1 ) )
               IF( DEBUG ) WRITE(*,*) MYROW,MYCOL,'EVK 2'
            ELSE
               EVK = ABS( SR( KWTOP-IROFFH+K-1 ) ) + 
     $              ABS( SI( KWTOP-IROFFH+K-1 ) )
               IF( DEBUG ) WRITE(*,*) MYROW,MYCOL,'EVK 3'
            END IF
            IF( DEBUG ) WRITE(*,*) MYROW,MYCOL,'EVK=',EVK
*
            IF( EVI.GE.EVK ) THEN
               I = K
            ELSE
               MLOC = 0
               SORTED = .false.
               IFST = I
               ILST = K
               DO 52 J = 1, I-1
                  SELECT( J ) = 1
                  MLOC = MLOC + 1
 52            CONTINUE
               IF( K.EQ.I+2 ) THEN
                  SELECT( I ) = 0
                  SELECT(I+1) = 0
               ELSE
                  SELECT( I ) = 0
               END IF
               IF( K.NE.KEND .AND. SI( KWTOP-IROFFH+K-1 ).NE.ZERO ) THEN
                  SELECT( K ) = 1
                  SELECT(K+1) = 1
                  MLOC = MLOC + 2
               ELSE
                  SELECT( K ) = 1
                  IF( K.LT.KEND ) SELECT(K+1) = 0
                  MLOC = MLOC + 1
               END IF
               DO 55 J = K+2, JW+IROFFH
                  SELECT( J ) = 0
 55            CONTINUE
               IF( DEBUG ) WRITE(*,*) MYROW,MYCOL,
     $              'PDLAQR3 -> PBDTRSEN 2'
               IF( DEBUG ) WRITE(*,*) MYROW,MYCOL, 
     $                 'SELECT=',SELECT(1:JW+IROFFH)
               CALL PBDTRSEN( 'No condition numbers', 'Vectors', 
     $              SELECT, PAR, JW+IROFFH, T, 1, 1, DESCT, V, 1, 1, 
     $              DESCV, WORK, WORK(JW+IROFFH+1), M, DDUM, DDUM, 
     $              WORK(2*(JW+IROFFH)+1), LWORK-2*(JW+IROFFH), 
     $              IWORK, LIWORK, IERR )
               IF( DEBUG ) WRITE(*,*) MYROW,MYCOL, 
     $              'PDLAQR3 <- PBDTRSEN 2' 
               CALL DCOPY( JW, WORK(1+IROFFH), 1, SR( KWTOP ), 1 )
               CALL DCOPY( JW, WORK(JW+2*IROFFH+1), 1, SI( KWTOP ), 1 )
               IF( IERR.EQ.0 ) THEN
                  I = ILST
               ELSE
                  I = K
               END IF
            END IF
            IF( I.EQ.KEND ) THEN
               K = I + 1
            ELSE IF( SI( KWTOP-IROFFH+I-1 ).EQ.ZERO ) THEN
               K = I + 1
            ELSE
               K = I + 2
            END IF
            GO TO 40
         END IF
         GO TO 30
   50    CONTINUE
      END IF
      IF( PRINT.AND.ITER.EQ.3 ) THEN
         CALL PDLAPRNT( JW+IROFFH, JW+IROFFH, T, 1, 1, DESCT, 0, 0, 
     $        'S5', 6, WORK )
         CALL PDLAPRNT( JW+IROFFH, JW+IROFFH, V, 1, 1, DESCV, 0, 0, 
     $        'Z5', 6, WORK ) 
      END IF
c      IF( DEBUG ) WRITE(*,*) MYROW,MYCOL, 'SR=',SR(1:KBOT)
c      IF( DEBUG ) WRITE(*,*) MYROW,MYCOL, 'SI=',SI(1:KBOT)
*
*     Restore number of rows and columns of T matrix descriptor
*
      DESCT( M_ ) = NW+IROFFH
      DESCT( N_ ) = NH+IROFFH
*
      IF( DEBUG ) WRITE(*,*) MYROW,MYCOL, 'NS,JW,S=',NS,JW,S
      IF( NS.LT.JW .OR. S.EQ.ZERO ) THEN
         IF( NS.GT.1 .AND. S.NE.ZERO ) THEN
*
*           ==== Reflect spike back into lower triangle ====
*
            IF( DEBUG ) WRITE(*,*) MYROW,MYCOL, 'Reflect spike back'   
            RROWS =  NUMROC( NS+IROFFH, NB, MYROW, DESCV(RSRC_), NPROW )
            RCOLS =  NUMROC( 1, 1, MYCOL, DESCV(CSRC_), NPCOL )
            CALL DESCINIT( DESCR, NS+IROFFH, 1, NB, 1, DESCV(RSRC_), 
     $                     DESCV(CSRC_), ICTXT, MAX(1, RROWS), INFO )
            TAUROWS = NUMROC( 1, 1, MYCOL, DESCV(RSRC_), NPROW )
            TAUCOLS = NUMROC( JW+IROFFH, NB, MYCOL, DESCV(CSRC_), 
     $                        NPCOL )
            CALL DESCINIT( DESCTAU, 1, JW+IROFFH, 1, NB, DESCV(RSRC_), 
     $                     DESCV(CSRC_), ICTXT, MAX(1, TAUROWS), INFO )
*
            IR = 1
            ITAU = IR + DESCR( LLD_ ) * RCOLS
            IPW  = ITAU + DESCTAU( LLD_ ) * TAUCOLS
*
            IF( DEBUG ) WRITE(*,*) MYROW,MYCOL, '...pdlaset...'
            CALL PDLASET( 'All', NS+IROFFH, 1, ZERO, ZERO, WORK(ITAU), 
     $                    1, 1, DESCTAU )
*
            IF( DEBUG ) WRITE(*,*) MYROW,MYCOL, '...pdcopy...'
            CALL PDCOPY( NS, V, 1+IROFFH, 1+IROFFH, DESCV, DESCV(M_), 
     $                   WORK(IR), 1+IROFFH, 1, DESCR, 1 )
           IF( DEBUG ) WRITE(*,*) MYROW,MYCOL, '...pdlarfg...' 
            CALL PDLARFG( NS, BETA, 1+IROFFH, 1, WORK(IR), 2+IROFFH, 1, 
     $                    DESCR, 1, WORK(ITAU+IROFFH) ) 
            IF( DEBUG ) WRITE(*,*) MYROW,MYCOL, '...pdelset...'
            CALL PDELSET( WORK(IR), 1+IROFFH, 1, DESCR, ONE )
*
            IF( DEBUG ) WRITE(*,*) MYROW,MYCOL, 
     $           '...pdlaset...again...'
            CALL PDLASET( 'Lower', JW-2, JW-2, ZERO, ZERO, T, 3+IROFFH, 
     $                    1+IROFFH, DESCT )
*
            IF( DEBUG ) WRITE(*,*) MYROW,MYCOL, '...pdlarf...'
            CALL PDLARF( 'Left', NS, JW, WORK(IR), 1+IROFFH, 1, DESCR,
     $                   1, WORK(ITAU+IROFFH), T, 1+IROFFH, 1+IROFFH, 
     $                   DESCT, WORK( IPW ) )
            IF( DEBUG ) WRITE(*,*) MYROW,MYCOL, '...pdlaset...again...'
            CALL PDLARF( 'Right', NS, NS, WORK(IR), 1+IROFFH, 1, DESCR,
     $                   1, WORK(ITAU+IROFFH), T, 1+IROFFH, 1+IROFFH, 
     $                   DESCT, WORK( IPW ) ) 
            IF( DEBUG ) WRITE(*,*) MYROW,MYCOL, 
     $           '...pdlarf...again...third time...'
            CALL PDLARF( 'Right', JW, NS, WORK(IR), 1+IROFFH, 1, DESCR,
     $                   1, WORK(ITAU+IROFFH), V, 1+IROFFH, 1+IROFFH, 
     $                   DESCV, WORK( IPW ) ) 
            IF( PRINT.AND.ITER.EQ.3 ) THEN
               CALL PDLAPRNT( JW+IROFFH, JW+IROFFH, T, 1, 1, DESCT,0,0, 
     $              'S6', 6, WORK(IPW) )
                CALL PDLAPRNT( JW+IROFFH, JW+IROFFH, V, 1, 1, DESCV,0,0, 
     $              'Z6', 6, WORK(IPW) ) 
            END IF   
*
            ITAU = 1
            IPW = ITAU + DESCTAU( LLD_ ) * TAUCOLS
            IF( DEBUG ) WRITE(*,*) MYROW,MYCOL, 'PDLAQR3->PDGEHRD',ITAU
            CALL PDGEHRD( JW+IROFFH, 1+IROFFH, NS+IROFFH, T, 1, 1, 
     $                    DESCT, WORK(ITAU), WORK( IPW ), LWORK-IPW+1, 
     $                    INFO )
            IF( PRINT.AND.ITER.EQ.3 ) THEN
               CALL PDLAPRNT( JW+IROFFH, JW+IROFFH, T, 1, 1, DESCT,0,0, 
     $              'S7', 6, WORK(IPW) )
               CALL PDLAPRNT( 1, JW+IROFFH, WORK(ITAU), 1, 1, DESCTAU,
     $              0,0, 'TAU1', 6, WORK(IPW) )
            END IF
            IF( DEBUG ) WRITE(*,*) MYROW,MYCOL, 'PDLAQR3<-PDGEHRD'
         END IF
*
*        ==== Copy updated reduced window into place ====
*
         IF( DEBUG ) WRITE(*,*) MYROW,MYCOL, 
     $        'Copy updated reduced window into place'
         IF( KWTOP.GT.1 ) THEN
            CALL PDELGET( 'All', '1-Tree', ELEM, V, 1+IROFFH, 
     $                    1+IROFFH, DESCV )
            CALL PDELSET( H, KWTOP, KWTOP-1, DESCH, S*ELEM )
         END IF
         CALL PDLACPY( 'Upper', JW-1, JW-1, T, 1+IROFFH+1, 1+IROFFH, 
     $                 DESCT, H, KWTOP+1, KWTOP, DESCH )
         CALL PDLACPY( 'All', 1, JW, T, 1+IROFFH, 1+IROFFH, DESCT, H,
     $                 KWTOP, KWTOP, DESCH )
         CALL PDLACPY( 'All', JW-1, 1, T, 1+IROFFH+1, 1+IROFFH+JW-1, 
     $                 DESCT, H, KWTOP+1, KWTOP+JW-1, DESCH )
*
*        ==== Accumulate orthogonal matrix in order to update
*        .    H and Z, if requested. ====
*        
         IF( NS.GT.1 .AND. S.NE.ZERO ) THEN
            IF( DEBUG ) WRITE(*,*) MYROW,MYCOL, 
     $           'Accumulate orthogonal matrix'
C            CALL PDORGHR( JW+IROFFH, 1+IROFFH, NS+IROFFH, T, 1, 1, 
C     $                    DESCT, WORK(ITAU), WORK( IPW ), LWORK-IPW+1, 
C     $                    INFO )
            IF( PRINT.AND.ITER.EQ.3 ) THEN
               CALL PDLAPRNT( JW+IROFFH, JW+IROFFH, V, 1, 1, DESCV,0,0, 
     $              'Z71', 6, WORK(IPW) )
               CALL PDLAPRNT( 1, JW+IROFFH, WORK(ITAU), 1, 1, DESCTAU,
     $              0,0, 'TAU2', 6, WORK(IPW) )
            END IF
            IF( DEBUG ) WRITE(*,*) MYROW,MYCOL, 'PDLAQR3->PDORMHR',ITAU
            CALL PDORMHR( 'Right', 'No', JW+IROFFH, NS+IROFFH, 1+IROFFH, 
     $                    NS+IROFFH, T, 1, 1, DESCT, WORK(ITAU), V, 1, 
     $                    1, DESCV, WORK( IPW ), LWORK-IPW+1, INFO )
            IF( PRINT.AND.ITER.EQ.3 ) THEN
               CALL PDLAPRNT( JW+IROFFH, JW+IROFFH, V, 1, 1, DESCV,0,0, 
     $              'Z72', 6, WORK(IPW) ) 
            END IF
C            CALL PDGEMM( 'No', 'No', JW, NS, NS, ONE, V, 1+IROFFH, 
C     $                   1+IROFFH, DESCV, T, 1+IROFFH, 1+IROFFH, DESCT, 
C     $                   ZERO, WV, 1+IROFFH, 1+IROFFH, DESCW )
C            CALL PDLACPY( 'All', JW, NS, WV, 1+IROFFH, 1+IROFFH, DESCW, 
C     $                    V, 1+IROFFH, 1+IROFFH, DESCV )
        END IF
*
*        ==== Update vertical slab in H ====
*
         IF( DEBUG ) WRITE(*,*) MYROW,MYCOL, 
     $       'Update vertical slab in H'
         IF( WANTT ) THEN
            LTOP = 1
         ELSE
            LTOP = KTOP
         END IF
         KLN = MAX( 0, KWTOP-LTOP )
         IROFFHH = MOD( LTOP-1, NB )
         ICOFFHH = MOD( KWTOP-1, NB )
         HHRSRC = INDXG2P( LTOP, NB, MYROW, DESCH(RSRC_), NPROW )
         HHCSRC = INDXG2P( KWTOP, NB, MYCOL, DESCH(CSRC_), NPCOL )
         HHROWS = NUMROC( KLN+IROFFHH, NB, MYROW, HHRSRC, NPROW )
         HHCOLS = NUMROC( JW+ICOFFHH, NB, MYCOL, HHCSRC, NPCOL )
         CALL DESCINIT( DESCHH, KLN+IROFFHH, JW+ICOFFHH, NB, NB, 
     $                  HHRSRC, HHCSRC, ICTXT, MAX(1, HHROWS), 
     $                  IERR )
         CALL PDGEMM( 'No', 'No', KLN, JW, JW, ONE, H, LTOP, 
     $                KWTOP, DESCH, V, 1+IROFFH, 1+IROFFH, DESCV, ZERO, 
     $                WORK, 1+IROFFHH, 1+ICOFFHH, DESCHH )
         CALL PDLACPY( 'All', KLN, JW, WORK, 1+IROFFHH, 1+ICOFFHH, 
     $                 DESCHH, H, LTOP, KWTOP, DESCH )
*
*        ==== Update horizontal slab in H ====
*
         IF( DEBUG )  WRITE(*,*) MYROW,MYCOL, 
     $        'Update horizontal slab in H'
         IF( WANTT ) THEN
            KLN = N-KBOT
            IROFFHH = MOD( KWTOP-1, NB )
            ICOFFHH = MOD( KBOT, NB )
            HHRSRC = INDXG2P( KWTOP, NB, MYROW, DESCH(RSRC_), NPROW )
            HHCSRC = INDXG2P( KBOT+1, NB, MYCOL, DESCH(CSRC_), NPCOL )
            HHROWS = NUMROC( JW+IROFFHH, NB, MYROW, HHRSRC, NPROW )
            HHCOLS = NUMROC( KLN+ICOFFHH, NB, MYCOL, HHCSRC, NPCOL )
            CALL DESCINIT( DESCHH, JW+IROFFHH, KLN+ICOFFHH, NB, NB, 
     $                     HHRSRC, HHCSRC, ICTXT, MAX(1, HHROWS), 
     $                     IERR )
            CALL PDGEMM( 'Tr', 'No', JW, KLN, JW, ONE, V, 
     $                   1+IROFFH, 1+IROFFH, DESCV, H, KWTOP, KBOT+1, 
     $                   DESCH, ZERO, WORK, 1+IROFFHH, 1+ICOFFHH, 
     $                   DESCHH )
            CALL PDLACPY( 'All', JW, KLN, WORK, 1+IROFFHH, 1+ICOFFHH, 
     $                    DESCHH, H, KWTOP, KBOT+1, DESCH )
         END IF
*
*        ==== Update vertical slab in Z ====
*     
         IF( DEBUG ) WRITE(*,*) MYROW,MYCOL, 
     $        'Update vertical slab in Z' 
         IF( WANTZ ) THEN
            KLN = IHIZ-ILOZ+1
            IROFFZZ = MOD( ILOZ-1, NB )
            ICOFFZZ = MOD( KWTOP-1, NB )
            ZZRSRC = INDXG2P( ILOZ, NB, MYROW, DESCZ(RSRC_), NPROW )
            ZZCSRC = INDXG2P( KWTOP, NB, MYCOL, DESCZ(CSRC_), NPCOL )
            ZZROWS = NUMROC( KLN+IROFFZZ, NB, MYROW, ZZRSRC, NPROW )
            ZZCOLS = NUMROC( JW+ICOFFZZ, NB, MYCOL, ZZCSRC, NPCOL )
            CALL DESCINIT( DESCZZ, KLN+IROFFZZ, JW+ICOFFZZ, NB, NB, 
     $                     ZZRSRC, ZZCSRC, ICTXT, MAX(1, ZZROWS), 
     $                     IERR )
            CALL PDGEMM( 'No', 'No', KLN, JW, JW, ONE, Z, ILOZ, 
     $                   KWTOP, DESCZ, V, 1+IROFFH, 1+IROFFH, DESCV, 
     $                   ZERO, WORK, 1+IROFFZZ, 1+ICOFFZZ, DESCZZ )
            CALL PDLACPY( 'All', KLN, JW, WORK, 1+IROFFZZ, 1+ICOFFZZ, 
     $                    DESCZZ, Z, ILOZ, KWTOP, DESCZ )
         END IF
      END IF
*
*     ==== Return the number of deflations ... ====
*
      ND = JW - NS
*
*     ==== ... and the number of shifts. (Subtracting
*     .    INFQR from the spike length takes care
*     .    of the case of a rare QR failure while
*     .    calculating eigenvalues of the deflation
*     .    window.)  ====
*
      NS = NS - INFQR
*
*      ==== Return optimal workspace. ====
*
      WORK( 1 ) = DBLE( LWKOPT )
*
*     ==== End of PDLAQR3 ====
*
      END
