      SUBROUTINE PDLAQR5( WANTT, WANTZ, KACC22, N, KTOP, KBOT, NSHFTS,
     $                    SR, SI, H, DESCH, ILOZ, IHIZ, Z, DESCZ, WORK, 
     $                    LWORK, IWORK, LIWORK, TIMINGS )
*
*  -- ScaLAPACK auxiliary routine (version 1.7.x) --
*     University of Umea and HPC2N, Umea, Sweden
*     February 2008
*
      IMPLICIT NONE
*
*     .. Scalar Arguments ..
      INTEGER            IHIZ, ILOZ, KACC22, KBOT, KTOP, N, NSHFTS, 
     $                   LWORK, LIWORK
      LOGICAL            WANTT, WANTZ
*     ..
*     .. Array Arguments ..
      INTEGER            DESCH( * ), DESCZ( * ), IWORK( * )
      DOUBLE PRECISION   H( * ), SI( * ), SR( * ), Z( * ), WORK( * ),
     $                   TIMINGS( 10 )
*     ..
*
*     This auxiliary subroutine called by PDLAQR0 performs a
*     single small-bulge multi-shift QR sweep by chasing separated
*     groups of bulges along the main block diagonal of H.
*
*      WANTT  (global input) logical scalar
*             WANTT = .TRUE. if the quasi-triangular Schur factor
*             is being computed.  WANTT is set to .FALSE. otherwise.
*
*      WANTZ  (global input) logical scalar
*             WANTZ = .TRUE. if the orthogonal Schur factor is being
*             computed.  WANTZ is set to .FALSE. otherwise.
*
*      KACC22 (global input) integer with value 0, 1, or 2.
*             Specifies the computation mode of far-from-diagonal
*             orthogonal updates.
*        = 0: PDLAQR5 does not accumulate reflections and does not
*             use matrix-matrix multiply to update far-from-diagonal
*             matrix entries.
*        = 1: PDLAQR5 accumulates reflections and uses matrix-matrix
*             multiply to update the far-from-diagonal matrix entries.
*        = 2: PDLAQR5 accumulates reflections, uses matrix-matrix
*             multiply to update the far-from-diagonal matrix entries,
*             and takes advantage of 2-by-2 block structure during
*             matrix multiplies.
*
*      N      (global input) integer scalar
*             N is the order of the Hessenberg matrix H upon which this
*             subroutine operates.
*
*      KTOP   (global input) integer scalar
*      KBOT   (global input) integer scalar
*             These are the first and last rows and columns of an
*             isolated diagonal block upon which the QR sweep is to be
*             applied. It is assumed without a check that
*                       either KTOP = 1  or   H(KTOP,KTOP-1) = 0
*             and
*                       either KBOT = N  or   H(KBOT+1,KBOT) = 0.
*
*      NSHFTS (global input) integer scalar
*             NSHFTS gives the number of simultaneous shifts.  NSHFTS
*             must be positive and even.
*
*      SR     (global input) DOUBLE PRECISION array of size (NSHFTS)
*      SI     (global input) DOUBLE PRECISION array of size (NSHFTS)
*             SR contains the real parts and SI contains the imaginary
*             parts of the NSHFTS shifts of origin that define the
*             multi-shift QR sweep.
*
*      H      (local input/output) DOUBLE PRECISION array of size 
*             (DESCH(LLD_),*)
*             On input H contains a Hessenberg matrix.  On output a
*             multi-shift QR sweep with shifts SR(J)+i*SI(J) is applied
*             to the isolated diagonal block in rows and columns KTOP
*             through KBOT.
*
*      DESCH  (global and local input) INTEGER array of dimension DLEN_.
*             The array descriptor for the distributed matrix H.
*
*      ILOZ   (global input) INTEGER
*      IHIZ   (global input) INTEGER
*             Specify the rows of Z to which transformations must be
*             applied if WANTZ is .TRUE.. 1 .LE. ILOZ .LE. IHIZ .LE. N
*
*      Z      (local input/output) DOUBLE PRECISION array of size
*             (DESCZ(LLD_),*)
*             If WANTZ = .TRUE., then the QR Sweep orthogonal
*             similarity transformation is accumulated into
*             Z(ILOZ:IHIZ,ILO:IHI) from the right.
*             If WANTZ = .FALSE., then Z is unreferenced.
*
*      DESCZ  (global and local input) INTEGER array of dimension DLEN_.
*             The array descriptor for the distributed matrix Z.
*
*      WORK   (local workspace) DOUBLE PRECISION array, dimension(DWORK)
*
*      LWORK  (local input) INTEGER
*             The length of the workspace array WORK.
*
*      IWORK  (local workspace) INTEGER array, dimension (ILWORK)
*
*      ILWORK (local input) INTEGER
*             The length of the workspace array IWORK
*
*     ================================================================
*     Based on contributions by
*        Robert Granat, Department of Computing Science and HPC2N,
*        University of Umea, Sweden.
*
*     ============================================================
*     Reference:
*
*     K. Braman, R. Byers and R. Mathias, The Multi-Shift QR
*     Algorithm Part I: Maintaining Well Focused Shifts, and
*     Level 3 Performance, SIAM Journal of Matrix Analysis,
*     volume 23, pages 929--947, 2002.
*
*     R. Granat and D. Kressner. A new implementation of the 
*     unsymmetric QR-algorithm for ScaLAPACK. Lapack Working
*     Note XYZ, 2008.  
*
*     ============================================================
*     .. Parameters ..
      INTEGER            BLOCK_CYCLIC_2D, CSRC_, CTXT_, DLEN_, DTYPE_,
     $                   LLD_, MB_, M_, NB_, N_, RSRC_
      LOGICAL            DEBUG, PRINT, CHKORTH
      PARAMETER          ( BLOCK_CYCLIC_2D = 1, DLEN_ = 9, DTYPE_ = 1,
     $                     CTXT_ = 2, M_ = 3, N_ = 4, MB_ = 5, NB_ = 6,
     $                     RSRC_ = 7, CSRC_ = 8, LLD_ = 9,
     $                     DEBUG = .FALSE., PRINT = .FALSE.,
     $                     CHKORTH = .FALSE. )
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0d0, ONE = 1.0d0 )
      INTEGER            NTINY
      PARAMETER          ( NTINY = 11 )
*     ..
*     .. Local Scalars ..
      DOUBLE PRECISION   ALPHA, BETA, H11, H12, H21, H22, REFSUM,
     $                   SAFMAX, SAFMIN, SCL, SMLNUM, SWAP, TST1, TST2,
     $                   ULP, TAU, ELEM, STAMP, DDUM, ORTH
      INTEGER            I, I2, I4, INCOL, J, J2, J4, JBOT, JCOL, JLEN,
     $                   JROW, JTOP, K, K1, KDU, KMS, KNZ, KRCOL, KZS,
     $                   M, M22, MBOT, MEND, MSTART, MTOP, NBMPS, NDCOL,
     $                   NS, NU, LLDH, LLDZ, LLDU, LLDV, LLDW, LLDWH,
     $                   INFO, ICTXT, NPROW, NPCOL, NB, IROFFH, ITOP,
     $                   NWIN, MYROW, MYCOL, LNS, NUMWIN, LKACC22,
     $                   LCHAIN, WIN, IDONEJOB, IPNEXT, ANMWIN, LENRBUF,
     $                   LENCBUF, ICHOFF, LRSRC, LCSRC, LKTOP, LKBOT,
     $                   II, JJ, SWIN, EWIN, LNWIN, DIM, LLKTOP, LLKBOT,
     $                   IPV, IPU, IPH, IPW, KU, KWH, KWV, NVE, LKS, 
     $                   IDUM, NHO, DIR, WINID, INDX, ILOC, JLOC, RSRC1,
     $                   CSRC1, RSRC2, CSRC2, RSRC3, CSRC3, RSRC4, IPUU,
     $                   CSRC4, LROWS, LCOLS, INDXS, KS, JLOC1, ILOC1,
     $                   LKTOP1, LKTOP2, WCHUNK, NUMCHUNK, ODDEVEN,
     $                   CHUNKNUM, DIM1, DIM4, IPW3, HROWS, ZROWS,
     $                   HCOLS, IPW1, IPW2, RSRC, EAST, JLOC4, ILOC4,
     $                   WEST, CSRC, SOUTH, NORHT, INDXE, NORTH,
     $                   IHH, IPIW, LKBOT1, NPROCS, LIROFFH,
     $                   WINFIN, SWEEP, RWS3, CLS3, INDX2, HROWS2, 
     $                   ZROWS2, HCOLS2, THREADS, MYTHREAD, MNRBUF,
     $                   MXRBUF, MNCBUF, MXCBUF, LWKOPT
      LOGICAL            ACCUM, BLK22, BMP22, INTRO, DONEJOB, ODDNPROW,
     $                   ODDNPCOL, LQUERY, BCDONE
      CHARACTER          JBCMPZ*2, JOB
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            PILAENVX, ICEIL, INDXG2P, INDXG2L, NUMROC,
     $                   OMP_GET_NUM_THREADS, OMP_GET_THREAD_NUM
      DOUBLE PRECISION   DLAMCH, MPI_WTIME, DLANGE, PCHKORTH
      EXTERNAL           DLAMCH, PILAENVX, ICEIL, INDXG2P, INDXG2L,
     $                   NUMROC, LSAME, MPI_WTIME, DLANGE,
     $                   PCHKORTH, OMP_GET_NUM_THREADS, 
     $                   OMP_GET_THREAD_NUM 
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, DBLE, MAX, MIN, MOD
*     ..
*     .. Local Arrays ..
      DOUBLE PRECISION   VT( 3 )
*     ..
*     .. External Subroutines ..
      EXTERNAL           DGEMM, DLABAD, DLACPY, DLAQR1, DLARFG, DLASET,
     $                   DTRMM, DLAQR6
*     ..
*     .. Executable Statements ..
*
      INFO = 0
      SWEEP = IWORK(1)
      ICTXT = DESCH( CTXT_ )
      CALL BLACS_GRIDINFO( ICTXT, NPROW, NPCOL, MYROW, MYCOL )
      NPROCS = NPROW*NPCOL
      LLDH = DESCH( LLD_ )
      LLDZ = DESCZ( LLD_ )
      NB = DESCH( MB_ )
      IROFFH = MOD( KTOP - 1, NB )
      LQUERY = LWORK.EQ.-1 .OR. LIWORK.EQ.-1
      IF( DEBUG ) WRITE(*,*)MYROW,MYCOL, 
     $     'Entering PDLAQR5: N,KTOP,KBOT=',
     $     N,KTOP,KBOT,NSHFTS
      IF( DEBUG ) WRITE(*,*) MYROW,MYCOL,'SWEEP=',SWEEP
*
*     ==== If there are no shifts, then there is nothing to do. ====
*
      IF( .NOT. LQUERY .AND. NSHFTS.LT.2 )
     $   RETURN
*
*     ==== If the active block is empty or 1-by-1, then there
*     .    is nothing to do. ====
*
      IF( .NOT. LQUERY .AND. KTOP.GE.KBOT )
     $   RETURN
*
*     ==== Shuffle shifts into pairs of real shifts and pairs
*     .    of complex conjugate shifts assuming complex
*     .    conjugate shifts are already adjacent to one
*     .    another. ====
*
      DO 10 I = 1, NSHFTS - 2, 2
         IF( SI( I ).NE.-SI( I+1 ) ) THEN
*
            SWAP = SR( I )
            SR( I ) = SR( I+1 )
            SR( I+1 ) = SR( I+2 )
            SR( I+2 ) = SWAP
*
            SWAP = SI( I )
            SI( I ) = SI( I+1 )
            SI( I+1 ) = SI( I+2 )
            SI( I+2 ) = SWAP
         END IF
   10 CONTINUE
*
*     ==== NSHFTS is supposed to be even, but if is odd,
*     .    then simply reduce it by one.  The shuffle above
*     .    ensures that the dropped shift is real and that
*     .    the remaining shifts are paired. ====
*
      NS = NSHFTS - MOD( NSHFTS, 2 )
*
*     === Extract the size of the computational window ===
*
      NWIN = PILAENVX( ICTXT, 19, 'PDLAQR5', JBCMPZ, N, NB, NB, NB )
      NWIN = MIN( NWIN, KBOT-KTOP+1 )
      IF( DEBUG ) WRITE(*,*)MYROW,MYCOL, 'PDLAQR5: NWIN =',NWIN
*
*     === Adjust number of simultaneous shifts if it exceeds the limit
*     .   set by the number of diagonal blocks in the active submatrix 
*     .   H(KTOP:KBOT,KTOP:KBOT). ===
*
      IF( DEBUG ) WRITE(*,*)MYROW,MYCOL, 'PDLAQR5: ICEIL 1'
      NS = MAX( 2, MIN( NS, ICEIL( KBOT-KTOP+1, NB )*NWIN/3 ) )
      NS = NS - MOD( NS, 2 )

      IF( DEBUG ) WRITE(*,*)MYROW,MYCOL, 'PDLAQR5: NS=',NS
*
*     === Decide the number of simultaneous computational windows 
*     .   from the number of shifts - each window should contain up to
*     .   (NWIN / 3) shifts. Also compute the number of shifts per
*     .   window and make sure that number is even. ===
*
      LNS = MIN( MAX( 2, NWIN / 3 ), MAX( 2, NS / MIN(NPROW,NPCOL) ) )
      IF( DEBUG ) WRITE(*,*)MYROW,MYCOL, 'PDLAQR5: LNS=',LNS
      LNS = LNS - MOD( LNS, 2 )
      IF( DEBUG ) WRITE(*,*)MYROW,MYCOL, 'PDLAQR5: LNS=',LNS
      IF( DEBUG ) WRITE(*,*)MYROW,MYCOL, 'PDLAQR5: ICEIL 2-3'
      IF( DEBUG ) WRITE(*,*)MYROW,MYCOL, 
     $     'PDLAQR5: NS, LNS, KBOT-KTOP+1, NB=',
     $     NS, LNS, KBOT-KTOP+1, NB
      NUMWIN = MAX( 1, MIN( ICEIL( NS, LNS ), 
     $     ICEIL( KBOT-KTOP+1, NB ) - 1 ) )
      IF( NPROW.NE.NPCOL ) THEN
         NUMWIN = MIN( NUMWIN, MIN(NPROW,NPCOL) )
         LNS = MIN( LNS, MAX( 2, NS / MIN(NPROW,NPCOL) ) )
         LNS = LNS - MOD( LNS, 2 )
      END IF
      IF( DEBUG ) WRITE(*,*)MYROW,MYCOL, 'NUMWIN =',NUMWIN
*
*     ==== Machine constants for deflation ====
*
      SAFMIN = DLAMCH( 'SAFE MINIMUM' )
      SAFMAX = ONE / SAFMIN
      CALL DLABAD( SAFMIN, SAFMAX )
      ULP = DLAMCH( 'PRECISION' )
      SMLNUM = SAFMIN*( DBLE( N ) / ULP )
*
*     ==== Use accumulated reflections to update far-from-diagonal
*     .    entries from a global perspective? ====
*
      ACCUM = ( KACC22.EQ.1 ) .OR. ( KACC22.EQ.2 )
      IF( DEBUG ) WRITE(*,*)MYROW,MYCOL, 'ACCUM,KACC22=',ACCUM,KACC22
*
*     ==== Use accumulated reflections to update far-from-diagonal
*     .    entries on a local level? ====
*
      IF( ACCUM ) THEN
         IF( LNS.LT.14 ) THEN
            LKACC22 = 1
         ELSE
            LKACC22 = 2
         END IF
      ELSE
         LKACC22 = 0
      END IF
*
*     ==== If so, exploit the 2-by-2 block structure? ====
*
      BLK22 = ( LNS.GT.2 ) .AND. ( KACC22.EQ.2 )
*
*     !!! TO BE REMOVED AFTER DEBUGGING !!!
*
      ACCUM = .TRUE.
      BLK22 = .FALSE.

      IF( DEBUG ) WRITE(*,*)MYROW,MYCOL, 'ACCUM,BLK22=',ACCUM,BLK22
*
*     ==== clear trash ====
*
      IF( .NOT. LQUERY .AND. KTOP+2.LE.KBOT ) THEN
         IF( DEBUG ) WRITE(*,*)MYROW,MYCOL, 'KTOP,KBOT=',KTOP,KBOT
         IF( DEBUG ) WRITE(*,*)MYROW,MYCOL,'pdelset'
         CALL PDELSET( H, KTOP+2, KTOP, DESCH, ZERO )
         IF( DEBUG ) WRITE(*,*)MYROW,MYCOL,'pdelset done!'
      END IF
*
*     ==== NBMPS = number of 2-shift bulges in each chain ====
*
      NBMPS = LNS / 2
*
*     ==== KDU = width of slab ====
*
      KDU = 6*NBMPS - 3
*
*     ==== LCHAIN = length of each chain 
*
      LCHAIN = 3 * NBMPS + 1
*
*     ==== Check if workspace query ====
*     !! TODO !!
*
      IF( LQUERY ) THEN
         IF( DEBUG ) WRITE(*,*)MYROW,MYCOL,'lquery - return!'
         HROWS = NUMROC( N, NB, MYROW, DESCH(RSRC_), NPROW )
         HCOLS = NUMROC( N, NB, MYCOL, DESCH(CSRC_), NPCOL )
         LWKOPT = (5+2*NUMWIN)*NB**2 + 2*HROWS*NB + HCOLS*NB +
     $        MAX( HROWS*NB, HCOLS*NB )
         WORK(1)  = DBLE(LWKOPT) 
         IWORK(1) = 5*NUMWIN
         RETURN
      END IF
*
*     Check that we do not have something stupid in the call
*
      IF( DEBUG ) WRITE(*,*)MYROW,MYCOL,'chk stupid'
      IF( KTOP.LT.1 .OR. KBOT.GT.N ) STOP
      IF( DEBUG ) WRITE(*,*)MYROW,MYCOL,'chk stupid done!'
*
*     Get number of OpenMP threads and their identities
*
      IF( DEBUG ) WRITE(*,*)MYROW,MYCOL,'OpenMP threads...'
#ifdef USE_OMP
C$OMP PARALLEL SHARED(THREADS), PRIVATE( MYTHREAD )
C$OMP MASTER
      THREADS = OMP_GET_NUM_THREADS()
C$OMP END MASTER
      MYTHREAD = OMP_GET_THREAD_NUM()
C$OMP END PARALLEL
#else
      THREADS = 1
      MYTHREAD = 0
#endif
      IF( DEBUG ) WRITE(*,*)MYROW,MYCOL,'THREADS=',THREADS
*
*     ==== Create and chase NUMWIN chains of NBMPS bulges ====
*
*     ==== Set up window introduction ===
*
      ANMWIN = 0
      INTRO = .TRUE.
      IPIW = 1
*
*     ==== While-loop over the computational windows which is
*     .    terminated when all windows have been introduced,
*     .    chased down to the bottom of the considered submatrix
*     .    and chased off. === 
*
      IF( DEBUG ) WRITE(*,*)MYROW,MYCOL, 
     $     'PDLAQR5: While-loop over comp win'
 20   CONTINUE
*
*     Set up next window as long as we have less than the prescribed
*     number of windows. Each window is described an integer quadruple:
*     1. Local value of KTOP (below denoted by LKTOP)
*     2. Local value of KBOT (below denoted by LKBOT)
*     3-4. Processor indices (LRSRC,LCSRC) associated with the window.
*     (5. Mark that decides if a window is fully processed or not)  
*
*     Notice - the next window is only introduced if the first block
*     in the active submatrix does not contain any other windows.
*
      IF( ANMWIN.GT.0 ) THEN
         LKTOP = IWORK( 1+(ANMWIN-1)*5 )
      ELSE
         LKTOP = KTOP
      END IF
      IF( INTRO .AND. (ANMWIN.EQ.0 .OR. LKTOP.GT.ICEIL(KTOP,NB)*NB) ) 
     $     THEN
         ANMWIN = ANMWIN + 1
         IWORK( 1+(ANMWIN-1)*5 ) = KTOP
         IF( DEBUG ) WRITE(*,*)MYROW,MYCOL, 'NWIN,NB,IROFFH,N,KTOP=',
     $        NWIN,NB,IROFFH,N,KTOP
         IWORK( 2+(ANMWIN-1)*5 ) = KTOP + 
     $                             MIN( NWIN,NB-IROFFH,KBOT-KTOP+1 ) - 1 
         IWORK( 3+(ANMWIN-1)*5 ) = INDXG2P( IWORK(1+(ANMWIN-1)*5), NB, 
     $                             MYROW, DESCH(RSRC_), NPROW )
         IWORK( 4+(ANMWIN-1)*5 ) = INDXG2P( IWORK(2+(ANMWIN-1)*5), NB, 
     $                             MYCOL, DESCH(CSRC_), NPCOL )
         IWORK( 5+(ANMWIN-1)*5 ) = 0
         IPIW = 6+(ANMWIN-1)*5
         IF( ANMWIN.EQ.NUMWIN ) INTRO = .FALSE.
         IF( DEBUG ) WRITE(*,*)MYROW,MYCOL, 
     $        'Set up window n:o',ANMWIN,'with data',
     $        IWORK(1+(ANMWIN-1)*5:5+(ANMWIN-1)*5)
      END IF
*
*     Do-loop over the number of windows
*
      IPNEXT = 1
      DONEJOB = .FALSE.
      IDONEJOB = 0
      LENRBUF = 0
      LENCBUF = 0
      ICHOFF = 0
      IF( DEBUG ) WRITE(*,*)MYROW,MYCOL, 'PDLAQR5: do-loop over num win'
      DO 40 WIN = 1, ANMWIN
         IF( DEBUG ) WRITE(*,*)MYROW,MYCOL, 'WIN=',WIN
*
*     Extract window information to simplify the rest
*
         LRSRC = IWORK( 3+(WIN-1)*5 )
         LCSRC = IWORK( 4+(WIN-1)*5 )
         LKTOP = IWORK( 1+(WIN-1)*5 )
         LKBOT = IWORK( 2+(WIN-1)*5 )
         LNWIN = LKBOT - LKTOP + 1
         IF( DEBUG ) WRITE(*,*)MYROW,MYCOL, 'LKTOP,LKBOT=',
     $        LKTOP,LKBOT
         IF( DEBUG ) WRITE(*,*)MYROW,MYCOL, 'LNWIN,MARK=',
     $        LNWIN,IWORK( 5+(ANMWIN-1)*5 )
*
*     Check if anything to do for current window, i.e., if the local
*     chain of bulges has reached the next block border etc. 
*
         IF( DEBUG ) WRITE(*,*)MYROW,MYCOL,'LNWIN,LCHAIN=',LNWIN,LCHAIN,
     $        LNWIN.GE.LCHAIN
         IF( DEBUG ) WRITE(*,*)MYROW,MYCOL, 'MARK=',IWORK(5+(WIN-1)*5)
         IF( IWORK(5+(WIN-1)*5).LT.2 .AND. LNWIN.GT.1 .AND. 
     $        (LNWIN.GT.LCHAIN .OR. LKBOT.EQ.KBOT ) ) THEN        
            IF( DEBUG ) WRITE(*,*)MYROW,MYCOL, 
     $           'LRSRC,LCSRC,LKTOP,LKBOT,LNWIN=',
     $           LRSRC,LCSRC,LKTOP,LKBOT,LNWIN
            LIROFFH = MOD(LKTOP-1,NB)
            IF( DEBUG ) WRITE(*,*)MYROW,MYCOL, 'LIROFFH=',LIROFFH
            SWIN = LKTOP-LIROFFH
            EWIN = MIN(KBOT,LKTOP-LIROFFH+NB-1)
            DIM = EWIN-SWIN+1
            IF( DEBUG ) WRITE(*,*)MYROW,MYCOL, 'LNWIN,DIM=',LNWIN,DIM
            IF( DIM.LE.NTINY .AND. .NOT.LKBOT.EQ.KBOT ) THEN
               IF( DEBUG ) WRITE(*,*)MYROW,MYCOL, 
     $              'Dim less than NTINY - skip!'
               IWORK( 5+(WIN-1)*5 ) = 2
               GO TO 45
            END IF
            IDONEJOB = 1
            IF( IWORK(5+(WIN-1)*5).EQ.0 ) THEN
               IWORK(5+(WIN-1)*5) = 1
            END IF    
*
*     Let the process that owns the corresponding window do the
*     local bulge chase
*
            IF( MYROW.EQ.LRSRC .AND. MYCOL.EQ.LCSRC ) THEN
               IF( DEBUG ) WRITE(*,*)MYROW,MYCOL, 
     $              'LRSRC,LCSRC=',LRSRC,LCSRC      
*
*     Set the kind of job to do in DLAQR6:
*     1. JOB = 'I': Introduce and chase bulges in window WIN
*     2. JOB = 'C': Chase bulges from top to bottom of window WIN
*     3. JOB = 'O': Chase bulges off window WIN
*     4. JOB = 'A': All of 1-3 above is done - this will for example happen
*                   for very small active submatrices (like 2-by-2) 
*     
               IF( DEBUG ) WRITE(*,*)MYROW,MYCOL, 
     $              'Select job for DLAQR6:'
               IF( DEBUG ) WRITE(*,*)MYROW,MYCOL, 
     $              'KTOP,KBOT,LKTOP,KBOT=',
     $              KTOP,KBOT,LKTOP,KBOT
               LLKBOT = LLKTOP + LNWIN - 1
               IF( LKTOP.EQ.KTOP .AND. LKBOT.EQ.KBOT ) THEN
                  JOB = 'All steps'
                  ICHOFF = 1
               ELSEIF( LKTOP.EQ.KTOP ) THEN
                  JOB = 'Introduce and chase'
               ELSEIF( LKBOT.EQ.KBOT ) THEN
                  JOB = 'Off-chase bulges'
                  ICHOFF = 1
               ELSE
                  JOB = 'Chase bulges'
               END IF
               IF( DEBUG ) WRITE(*,*)MYROW,MYCOL, 'JOB=',JOB
*     
*     Copy submatrix of H corresponding to window WIN into
*     workspace and set ut additional workspace for storing
*     orthogonal transformations. This submatrix must be 
*     at least (NTINY+1)-by-(NTINY+1) to fit into DLAQR6 - if not, abort
*     and go for cross border reordering with this particular
*     window.
*
               II = INDXG2L( SWIN, NB, MYROW, DESCH(RSRC_), NPROW )
               JJ = INDXG2L( SWIN, NB, MYCOL, DESCH(CSRC_), NPCOL )
               LLKTOP = 1 + LIROFFH
               LLKBOT = LLKTOP + LNWIN - 1
               IF( DEBUG ) WRITE(*,*)MYROW,MYCOL, 
     $              'LLKTOP,LLKBOT=',LLKTOP,LLKBOT
*
               IPU = IPNEXT
               IF( DEBUG ) WRITE(*,*)MYROW,MYCOL, 'IPU=',IPU
               IF( .NOT. ACCUM ) THEN
                  IPH = IPU + 3*LNWIN
               ELSE
                  IPH = IPU + LNWIN**2
               END IF
               IPUU = IPH + MAX(NTINY+1,DIM)**2
               IPV = IPUU + MAX(NTINY+1,DIM)**2
               IPNEXT = IPH               
* 
               IF( LSAME( JOB, 'A' ) .OR. LSAME( JOB, 'O' ) .AND.
     $             DIM.LT.NTINY+1 ) THEN
                  CALL DLASET( 'All', NTINY+1, NTINY+1, ZERO, ONE, 
     $                 WORK(IPH), NTINY+1 )
               END IF
               IF( DEBUG ) WRITE(*,*)MYROW,MYCOL, 
     $              'Copy upper part of comp win'
               CALL DLACPY( 'Upper', DIM, DIM, H(II+(JJ-1)*LLDH), LLDH,
     $                      WORK(IPH), MAX(NTINY+1,DIM) ) 
               IF( DEBUG ) WRITE(*,*)MYROW,MYCOL, 
     $              'Copy subd. part of comp win'
               CALL DCOPY(  DIM-1, H(II+(JJ-1)*LLDH+1), LLDH+1, 
     $                      WORK(IPH+1), MAX(NTINY+1,DIM)+1 )
               IF( LSAME( JOB, 'C' ) .OR. LSAME( JOB, 'O') ) THEN
                  CALL DCOPY(  DIM-2, H(II+(JJ-1)*LLDH+2), LLDH+1, 
     $                 WORK(IPH+2), MAX(NTINY+1,DIM)+1 ) 
                  CALL DCOPY(  DIM-3, H(II+(JJ-1)*LLDH+3), LLDH+1, 
     $                 WORK(IPH+3), MAX(NTINY+1,DIM)+1 ) 
                  CALL DLASET( 'Lower', DIM-4, DIM-4, ZERO, 
     $                 ZERO, WORK(IPH+4), MAX(NTINY+1,DIM) )
               ELSE
                  CALL DLASET( 'Lower', DIM-2, DIM-2, ZERO, 
     $                 ZERO, WORK(IPH+2), MAX(NTINY+1,DIM) )
               END IF
*
               KU = MAX(NTINY+1,DIM) - KDU + 1
               IF( DEBUG ) WRITE(*,*)MYROW,MYCOL, 'KU=',KU
               KWH = KDU + 1
               IF( DEBUG ) WRITE(*,*)MYROW,MYCOL, 'KWH=',KWH
               NHO = ( MAX(NTINY+1,DIM)-KDU+1-4 ) - ( KDU+1 ) + 1
               IF( DEBUG ) WRITE(*,*)MYROW,MYCOL, 'NHO=',NHO
               KWV = KDU + 4
               IF( DEBUG ) WRITE(*,*)MYROW,MYCOL, 'KWV=',KWV
               NVE = MAX(NTINY+1,DIM) - KDU - KWV + 1
               IF( DEBUG ) WRITE(*,*)MYROW,MYCOL, 'NVE=',NVE
               IF( DEBUG ) WRITE(*,*)MYROW,MYCOL, 
     $              'Set local U to identity'
               IF( ACCUM )
     $              CALL DLASET( 'All', MAX(NTINY+1,DIM), 
     $              MAX(NTINY+1,DIM), ZERO, ONE, WORK(IPUU), 
     $              MAX(NTINY+1,DIM) )
*
*              ==== Small-bulge multi-shift QR sweep ====
*
               IF( DEBUG ) WRITE(*,*)MYROW,MYCOL, 'LKBOT,LNS=',LKBOT,LNS
               LKS = MAX( 1, NS - WIN*LNS + 1 )
               IF( DEBUG ) WRITE(*,*)MYROW,MYCOL, 
     $              'WIN,NS,LNS,LKS=',WIN,NS,LNS,LKS
               IF( DEBUG ) WRITE(*,*)MYROW,MYCOL, 'PDLAQR5 -> DLAQR6'
               IF( PRINT ) THEN
                  CALL DLAPRNT(DIM,DIM,WORK(IPH),1,1,
     $                 MAX(NTINY+1,DIM),'HH1',6)
                  CALL DLAPRNT(DIM,DIM,WORK(IPUU),1,1,MAX(NTINY+1,DIM),
     $                 'UU1',6)
c                  CALL DLAPRNT(3,LNWIN,WORK(IPV),1,1,3,'VV1',6)
               END IF
               IF( DEBUG ) WRITE(*,*)MYROW,MYCOL,'SR=',SR(LKS:LKS+LNS-1)
               IF( DEBUG ) WRITE(*,*)MYROW,MYCOL,'SI=',SI(LKS:LKS+LNS-1)
               STAMP = MPI_WTIME()
               CALL DLAQR6( JOB, WANTT, ACCUM, LKACC22, 
     $              MAX(NTINY+1,DIM), LLKTOP, LLKBOT, LNS, SR( LKS ), 
     $              SI( LKS ), WORK(IPH), MAX(NTINY+1,DIM), LLKTOP, 
     $              LLKBOT, WORK(IPUU), MAX(NTINY+1,DIM), WORK(IPV), 
     $              3, IWORK(IPIW), 2, WORK(IPU), 3, WORK( IPH+KU-1 ), 
     $              MAX(NTINY+1,DIM), NVE, WORK( IPH+KWV-1 ), 
     $              MAX(NTINY+1,DIM), NHO, WORK( IPH-1+KU+(KWH-1)*
     $              MAX(NTINY+1,DIM) ), MAX(NTINY+1,DIM) )
               TIMINGS( 6 ) = TIMINGS( 6 ) + MPI_WTIME() - STAMP
               IF( PRINT ) THEN
                  CALL DLAPRNT(DIM,DIM,WORK(IPH),1,1,
     $                 MAX(NTINY+1,DIM),'HH2',6)
                  CALL DLAPRNT(DIM,DIM,WORK(IPUU),1,1,
     $                 MAX(NTINY+1,DIM),'UU2',6)
C               CALL DLAPRNT(3,LNWIN,WORK(IPV),1,1,3,'VV2',6)
               END IF
               IF( DEBUG ) WRITE(*,*)MYROW,MYCOL, 'PDLAQR5 <- DLAQR6'
*
*     Copy submatrix of H back
*
               CALL DLACPY( 'Upper', DIM, DIM, WORK(IPH), 
     $              MAX(NTINY+1,DIM), H(II+(JJ-1)*LLDH), LLDH )
               CALL DCOPY( DIM-1, WORK(IPH+1), MAX(NTINY+1,DIM)+1, 
     $                     H(II+(JJ-1)*LLDH+1), LLDH+1 )
               IF( LSAME( JOB, 'I' ) .OR. LSAME( JOB, 'C' ) ) THEN
                  CALL DCOPY( DIM-2, WORK(IPH+2), DIM+1, 
     $                 H(II+(JJ-1)*LLDH+2), LLDH+1 )
                  CALL DCOPY( DIM-3, WORK(IPH+3), DIM+1, 
     $                 H(II+(JJ-1)*LLDH+3), LLDH+1 )
               ELSE
                  CALL DLASET( 'Lower', DIM-2, DIM-2, ZERO, 
     $                 ZERO, H(II+(JJ-1)*LLDH+2), LLDH )
               END IF
*
*     Copy actual submatrix of U to the correct place of the buffer
*
c               IF( ACCUM )
c     $              CALL DLACPY( 'All', LNWIN, LNWIN, 
c     $              WORK(IPUU+(DIM-LNWIN)*MAX(NTINY+1,DIM)+DIM-LNWIN), 
c     $              MAX(NTINY+1,DIM), WORK(IPU), LNWIN  )
               IF( ACCUM )
     $              CALL DLACPY( 'All', LNWIN, LNWIN, 
     $              WORK(IPUU+(MAX(NTINY+1,DIM)*LIROFFH)+LIROFFH), 
     $              MAX(NTINY+1,DIM), WORK(IPU), LNWIN  )
c               CALL DLAPRNT(LNWIN,LNWIN,WORK(IPU),1,1,
c     $              LNWIN,'U2',6)
            END IF
*
*     In case the local submatrix was smaller than 
*     (NTINY+1)-by-(NTINY+1) we go here and proceed.
*
 45         CONTINUE
         ELSE
            IF( DEBUG ) WRITE(*,*)MYROW,MYCOL, 
     $           'Lack of local space  - skip and go for crossborder!'
            IWORK( 5+(WIN-1)*5 ) = 2
         END IF
*
*     Increment counter for buffers of orthogonal transformations
*
         IF( MYROW.EQ.LRSRC .OR. MYCOL.EQ.LCSRC ) THEN
            IF( IDONEJOB.EQ.1 .AND. IWORK(5+(WIN-1)*5).LT.2 ) THEN
               IF( .NOT. ACCUM ) THEN
                  IF( MYROW.EQ.LRSRC ) LENRBUF = LENRBUF + 3*LNWIN
                  IF( MYCOL.EQ.LCSRC ) LENCBUF = LENCBUF + 3*LNWIN
               ELSE
                  IF( MYROW.EQ.LRSRC ) LENRBUF = LENRBUF + LNWIN*LNWIN
                  IF( MYCOL.EQ.LCSRC ) LENCBUF = LENCBUF + LNWIN*LNWIN
               END IF 
            END IF
            IF( DEBUG ) WRITE(*,*)MYROW,MYCOL, 'LENRBUF,LENCBUF=',
     $           LENRBUF,LENCBUF 
         END IF
 40   CONTINUE
*
*     Did some work in the above do-loop? 
*
      IF( DEBUG ) WRITE(*,*)MYROW,MYCOL, 'PDLAQR5: IGSUM2D IDONEJOB'
      CALL IGSUM2D( ICTXT, 'All', '1-Tree', 1, 1, IDONEJOB, 1, -1, -1 )
      DONEJOB = IDONEJOB.GT.0
      IF( DEBUG ) WRITE(*,*)MYROW,MYCOL, 'DONEJOB=',DONEJOB
*
*     Chased off bulges from first window?
*
      IF( DEBUG ) WRITE(*,*)MYROW,MYCOL, 'PDLAQR5: IGAMX2D ICHOFF 1'
      IF( NPROCS.GT.1 )
     $     CALL IGAMX2D( ICTXT, 'All', '1-Tree', 1, 1, ICHOFF, 1, -1, 
     $                   -1, -1, -1, -1 )
      IF( DEBUG ) WRITE(*,*)MYROW,MYCOL, 'ICHOFF=',ICHOFF
*
*     If work was done in the do-loop over local windows, perform
*     updates, otherwise go for cross border reordering and updates
*
      IF( DONEJOB ) THEN  
         IF( DEBUG ) WRITE(*,*)MYROW,MYCOL, 'PDLAQR5: local broadcasts'
*
*     In debugging, check for inconsistency in buffer lengths
*
         IF( DEBUG ) THEN
            MNRBUF = LENRBUF
            MXRBUF = LENRBUF
            MNCBUF = LENCBUF
            MXCBUF = LENCBUF
            IF( NPROCS.GT.1 )
     $           CALL IGAMN2D( ICTXT, 'Row', '1-Tree', 1, 1, MNRBUF, 1, 
     $           -1, -1, -1, -1, -1 )
            IF( NPROCS.GT.1 )
     $           CALL IGAMX2D( ICTXT, 'Row', '1-Tree', 1, 1, MXRBUF, 1, 
     $           -1, -1, -1, -1, -1 )
            IF( NPROCS.GT.1 )
     $           CALL IGAMN2D( ICTXT, 'Col', '1-Tree', 1, 1, MNCBUF, 1, 
     $           -1, -1, -1, -1, -1 )
            IF( NPROCS.GT.1 )
     $           CALL IGAMX2D( ICTXT, 'Col', '1-Tree', 1, 1, MXCBUF, 1, 
     $           -1, -1, -1, -1, -1 )
            IF( MNRBUF.NE.MXRBUF ) THEN
               IF( MYROW.EQ.MYCOL ) THEN
                  WRITE(*,*) '*** Error before local BC ***'
                  WRITE(*,*) '*** Inconsistency in row buffer lengths:'
                  WRITE(*,*) 'Min,Max,Corr,row=',MNRBUF,MXRBUF,
     $                 LENRBUF,MYROW
               END IF
               STOP
            END IF
            IF( MNCBUF.NE.MXCBUF ) THEN
               IF( MYROW.EQ.MYCOL ) THEN
                  WRITE(*,*) '*** Error before local BC ***'
                  WRITE(*,*) 
     $                 '*** Inconsistency in column buffer lengths:'
                  WRITE(*,*) 'Min,Max,Corr,col=',MNCBUF,MXCBUF,
     $                 LENCBUF,MYCOL
               END IF
               STOP
            END IF
         END IF
*     
*     Broadcast orthogonal transformations 
*     
 49      CONTINUE
         IF( LENRBUF.GT.0 .OR. LENCBUF.GT.0 ) THEN
            DO 50 DIR = 1, 2
               BCDONE = .FALSE.
               DO 60 WIN = 1, ANMWIN
                  IF( ( LENRBUF.EQ.0 .AND. LENCBUF.EQ.0 ) .OR. 
     $                 BCDONE ) GOTO 62
                  LRSRC = IWORK( 3+(WIN-1)*5 )
                  LCSRC = IWORK( 4+(WIN-1)*5 )
                  IF( MYROW.EQ.LRSRC .AND. MYCOL.EQ.LCSRC ) THEN 
                     IF( DIR.EQ.1 .AND. LENRBUF.GT.0 .AND.
     $                    NPCOL.GT.1 ) THEN
                        IF( DEBUG ) WRITE(*,*)MYROW,MYCOL,'Loc. BC row',
     $                       LENRBUF
                        CALL DGEBS2D( ICTXT, 'Row', '1-Tree', LENRBUF, 
     $                       1, WORK, LENRBUF )
                        IF( DEBUG ) WRITE(*,*)MYROW,MYCOL, 
     $                       'Loc. BC row done!'
                     ELSEIF( DIR.EQ.2 .AND. LENCBUF.GT.0 .AND.
     $                       NPROW.GT.1 ) THEN
                        IF( DEBUG ) WRITE(*,*)MYROW,MYCOL,'Loc. BC col',
     $                       LENCBUF
                        CALL DGEBS2D( ICTXT, 'Col', '1-Tree', LENCBUF, 
     $                       1, WORK, LENCBUF )
                        IF( DEBUG ) WRITE(*,*)MYROW,MYCOL, 
     $                       'Loc. BC col done!'
                     END IF
                     IF( LENRBUF.GT.0 )
     $                    CALL DLACPY( 'All', LENRBUF, 1, WORK, LENRBUF, 
     $                    WORK(1+LENRBUF), LENCBUF )
                     BCDONE = .TRUE.
                  ELSEIF( MYROW.EQ.LRSRC .AND. DIR.EQ.1 ) THEN
                     IF( LENRBUF.GT.0 .AND. NPCOL.GT.1 ) THEN
                        IF( DEBUG ) WRITE(*,*)MYROW,MYCOL, 
     $                       'Loc. BC row recv',LENRBUF
                        CALL DGEBR2D( ICTXT, 'Row', '1-Tree', LENRBUF, 
     $                       1, WORK, LENRBUF, LRSRC, LCSRC )
                        BCDONE = .TRUE.
                        IF( DEBUG ) WRITE(*,*)MYROW,MYCOL, 
     $                       'Loc. BC row recv done!'
                     END IF
                  ELSEIF( MYCOL.EQ.LCSRC .AND. DIR.EQ.2 ) THEN
                     IF( LENCBUF.GT.0 .AND. NPROW.GT.1 ) THEN
                        IF( DEBUG ) WRITE(*,*)MYROW,MYCOL, 
     $                       'Loc. BC col recv',LENCBUF
                        CALL DGEBR2D( ICTXT, 'Col', '1-Tree', LENCBUF, 
     $                       1, WORK(1+LENRBUF), LENCBUF, LRSRC, LCSRC )
                        BCDONE = .TRUE.
                        IF( DEBUG ) WRITE(*,*)MYROW,MYCOL, 
     $                       'Loc. BC col recv done!'
                     END IF
                  END IF
 62               CONTINUE
 60            CONTINUE
 50         CONTINUE
         END IF
         IF( DEBUG .AND. PRINT ) THEN
            CALL BLACS_BARRIER( ICTXT, 'All' )
            IF( MYROW.EQ.0 .AND. MYCOL.EQ.0 ) THEN
               CALL DLAPRNT( LENRBUF, 1, WORK, 1, 1, LENRBUF, 'R00', 6 )
               CALL DLAPRNT( LENCBUF, 1, WORK(1+LENRBUF), 1, 1, LENCBUF, 
     $              'C00', 6 )
            END IF
            IF( MYROW.EQ.0 .AND. MYCOL.EQ.1 ) THEN
               CALL DLAPRNT( LENRBUF, 1, WORK, 1, 1, LENRBUF, 'R01', 6 )
               CALL DLAPRNT( LENCBUF, 1, WORK(1+LENRBUF), 1, 1, LENCBUF, 
     $              'C01', 6 )
            END IF
            IF( MYROW.EQ.1 .AND. MYCOL.EQ.0 ) THEN
               CALL DLAPRNT( LENRBUF, 1, WORK, 1, 1, LENRBUF, 'R10', 6 )
               CALL DLAPRNT( LENCBUF, 1, WORK(1+LENRBUF), 1, 1, LENCBUF, 
     $              'C10', 6 )
            END IF
            IF( MYROW.EQ.1 .AND. MYCOL.EQ.1 ) THEN
               CALL DLAPRNT( LENRBUF, 1, WORK, 1, 1, LENRBUF, 'R11', 6 )
               CALL DLAPRNT( LENCBUF, 1, WORK(1+LENRBUF), 1, 1, LENCBUF, 
     $              'C11', 6 )
            END IF
            CALL BLACS_BARRIER( ICTXT, 'All' )
         END IF
*     
*     Compute updates - make sure to skip windows that was skipped
*     regarding local reordering
*     
         IF( DEBUG ) WRITE(*,*)MYROW,MYCOL, 'PDLAQR5: Local updates'
         DO 65 DIR = 1, 2
            IF( DEBUG ) WRITE(*,*)MYROW,MYCOL, 'DIR=',DIR
            WINID = 0
            IF( DIR.EQ.1 ) THEN
               IPNEXT = 1
            ELSE
               IPNEXT = 1 + LENRBUF
            END IF
            DO 70 WIN = 1, ANMWIN
               IF( DEBUG ) WRITE(*,*)MYROW,MYCOL, 'WIN=',WIN
               IF( IWORK( 5+(WIN-1)*5 ).EQ.2 ) GOTO 75
               LRSRC = IWORK( 3+(WIN-1)*5 )
               LCSRC = IWORK( 4+(WIN-1)*5 )
               LKTOP = IWORK( 1+(WIN-1)*5 )
               LKBOT = IWORK( 2+(WIN-1)*5 )
               LNWIN = LKBOT - LKTOP + 1
               IF( (MYROW.EQ.LRSRC.AND.LENRBUF.GT.0.AND.DIR.EQ.1) .OR. 
     $             (MYCOL.EQ.LCSRC.AND.LENCBUF.GT.0.AND.DIR.EQ.2 ) ) 
     $              THEN
*
*     Set up workspaces
*
                  IPU = IPNEXT
                  IF( DEBUG ) WRITE(*,*)MYROW,MYCOL, 'IPU=',IPU
                  IF( ACCUM ) THEN
                     IPNEXT = IPU + LNWIN*LNWIN
                  ELSE
                     IPNEXT = IPU + 3*LNWIN
                  END IF
                  IPW = 1 + LENRBUF + LENCBUF
                  LIROFFH = MOD(LKTOP-1,NB)
                  WINID = WINID + 1
*
*     Recompute JOB to see if block structure of U could possibly
*     be exploited or not
*
                  IF( LKTOP.EQ.KTOP .AND. LKBOT.EQ.KBOT ) THEN
                     JOB = 'All steps'
                  ELSEIF( LKTOP.EQ.KTOP ) THEN
                     JOB = 'Introduce and chase'
                  ELSEIF( LKBOT.EQ.KBOT ) THEN
                     JOB = 'Off-chase bulges'
                  ELSE
                     JOB = 'Chase bulges'
                  END IF
                  IF( DEBUG ) WRITE(*,*)MYROW,MYCOL, 'JOB=',JOB
               END IF
*     
*     ==== Use U (if accumulated) to update far-from-diagonal
*     .    entries in H.  If required, use U to update Z as
*     .    well. ====
*     
               IF( ACCUM .AND. (.NOT. BLK22 .OR. LSAME(JOB,'I') 
     $              .OR. LSAME(JOB,'O') .OR. LNS.LE.2 ) ) THEN
*
                  IF( DIR.EQ.2 .AND. LENCBUF.GT.0 .AND. 
     $                 MYCOL.EQ.LCSRC ) THEN
                     IF( WANTT ) THEN
                        IF( DEBUG ) WRITE(*,*)MYROW,MYCOL, 
     $                       'PDLAQR5: ICEIL 4'
#ifdef USE_OMP
*$OMP PARALLEL DEFAULT(NONE), 
*$OMP& PRIVATE(INDX,LROWS,ILOC,JLOC,MYTHREAD,RSRC1,CSRC1),
*$OMP& SHARED(H,WORK,LKTOP,LIROFFH,NB,DESCH,NPROW,NPCOL,MYROW,MYCOL,
*$OMP& LNWIN,LLDH,ONE,IPU,ZERO,IPW)        
                        MYTHREAD = OMP_GET_THREAD_NUM()
*$OMP DO
#endif
                        DO 80 INDX = 1, LKTOP-LIROFFH-1, NB
                           IF( DEBUG ) WRITE(*,*)MYROW,MYCOL, 
     $                          'INDX=',INDX
                           CALL INFOG2L( INDX, LKTOP, DESCH, NPROW, 
     $                          NPCOL, MYROW, MYCOL, ILOC, JLOC, RSRC1, 
     $                          CSRC1 )
                           IF( MYROW.EQ.RSRC1.AND.MYCOL.EQ.CSRC1 ) THEN
                              LROWS = MIN( NB, LKTOP-INDX )
                              IF( DEBUG ) WRITE(*,*)MYROW,MYCOL, 
     $                             'LROWS=',LROWS
                              IF( PRINT )
     $                             CALL DLAPRNT(LNWIN,LNWIN,WORK(IPU),
     $                             1,1,LNWIN,'UUU',6)
                              CALL DGEMM('No transpose', 'No transpose', 
     $                             LROWS, LNWIN, LNWIN, ONE, 
     $                             H((JLOC-1)*LLDH+ILOC), LLDH, 
     $                             WORK( IPU ), LNWIN, ZERO, 
     $                             WORK(IPW+MYTHREAD*NB**2), 
     $                             LROWS )
                              CALL DLACPY( 'All', LROWS, LNWIN, 
     $                             WORK(IPW+MYTHREAD*NB**2), LROWS, 
     $                             H((JLOC-1)*LLDH+ILOC), LLDH )
                           END IF
 80                     CONTINUE
#ifdef USE_OMP
*$OMP END DO 
*$OMP END PARALLEL
#endif
                     END IF
                     IF( WANTZ ) THEN
#ifdef USE_OMP
*$OMP PARALLEL DEFAULT(NONE), 
*$OMP& PRIVATE(INDX,LROWS,ILOC,JLOC,MYTHREAD,RSRC1,CSRC1),
*$OMP& SHARED(Z,WORK,LKTOP,LIROFFH,NB,DESCZ,NPROW,NPCOL,MYROW,MYCOL,
*$OMP& LNWIN,LLDZ,ONE,IPU,ZERO,IPW,N,WIN,DDUM)     
                        MYTHREAD = OMP_GET_THREAD_NUM()
*$OMP DO
#endif
                        DO 90 INDX = 1, N, NB
                           IF( DEBUG ) WRITE(*,*)MYROW,MYCOL, 
     $                          'INDX,LKTOP=',INDX,LKTOP
                           CALL INFOG2L( INDX, LKTOP, DESCZ, NPROW, 
     $                          NPCOL, MYROW, MYCOL, ILOC, JLOC, RSRC1, 
     $                          CSRC1 )
                           IF( DEBUG ) WRITE(*,*) 'ILOC,JLOC=',
     $                          ILOC,JLOC
                           IF( MYROW.EQ.RSRC1.AND.MYCOL.EQ.CSRC1 ) THEN
                              LROWS = MIN(NB,N-INDX+1)
                              IF( DEBUG ) WRITE(*,*)MYROW,MYCOL, 
     $                             'LROWS=',LROWS
                              IF( DEBUG ) WRITE(*,*) MYROW,MYCOL,
     $                             'Z*Uloc - WIN,IPU=',WIN,IPU
                              IF( PRINT )
     $                             CALL DLAPRNT(LNWIN,LNWIN,
     $                             WORK(IPU),1,1,LNWIN,'Uloc',6)
                              CALL DGEMM( 'No transpose', 
     $                             'No transpose', LROWS, LNWIN, LNWIN, 
     $                             ONE, Z((JLOC-1)*LLDZ+ILOC), LLDZ, 
     $                             WORK( IPU ), LNWIN, ZERO, 
     $                             WORK(IPW+MYTHREAD*NB**2), LROWS )
                              CALL DLACPY( 'All', LROWS, LNWIN, 
     $                             WORK(IPW+MYTHREAD*NB**2), LROWS, 
     $                             Z((JLOC-1)*LLDZ+ILOC), LLDZ )
                           END IF
 90                     CONTINUE
#ifdef USE_OMP
*$OMP END DO 
*$OMP END PARALLEL
#endif
                     END IF
                  END IF
*     
*     Update the rows of T affected by the bulge-chase
*     
                  IF( DIR.EQ.1 .AND. LENRBUF.GT.0 .AND. 
     $                 MYROW.EQ.LRSRC ) THEN
                     IF( WANTT ) THEN
                        IF( ICEIL(LKBOT,NB).EQ.ICEIL(KBOT,NB) ) THEN
                           LCOLS = MIN(ICEIL(KBOT,NB)*NB,N) - KBOT
                        ELSE
                           LCOLS = 0
                        END IF
                        IF( LCOLS.GT.0 ) THEN
                           INDX = KBOT + 1
                           IF( DEBUG ) WRITE(*,*)MYROW,MYCOL, 
     $                          'INDX=',INDX
                           CALL INFOG2L( LKTOP, INDX, DESCH, NPROW, 
     $                          NPCOL, MYROW, MYCOL, ILOC, JLOC, 
     $                          RSRC1, CSRC1 )
                           IF( MYROW.EQ.RSRC1.AND.MYCOL.EQ.CSRC1 ) THEN
                              IF( DEBUG ) WRITE(*,*)MYROW,MYCOL, 
     $                             'LCOLS=',LCOLS
                              CALL DGEMM( 'Transpose', 'No Transpose', 
     $                             LNWIN, LCOLS, LNWIN, ONE, WORK(IPU), 
     $                             LNWIN, H((JLOC-1)*LLDH+ILOC), LLDH, 
     $                             ZERO, WORK(IPW), LNWIN )
                              CALL DLACPY( 'All', LNWIN, LCOLS, 
     $                             WORK(IPW), LNWIN, 
     $                             H((JLOC-1)*LLDH+ILOC), LLDH )
                           END IF
                        END IF
 93                     CONTINUE
                        IF( DEBUG ) WRITE(*,*)MYROW,MYCOL, 
     $                       'PDLAQR5: ICEIL 5'
                        INDXS = ICEIL(LKBOT,NB)*NB + 1
#ifdef USE_OMP
*$OMP PARALLEL DEFAULT(NONE), 
*$OMP& PRIVATE(INDX,LCOLS,ILOC,JLOC,MYTHREAD,RSRC1,CSRC1),
*$OMP& SHARED(H,WORK,LKTOP,LIROFFH,NB,DESCH,NPROW,NPCOL,MYROW,MYCOL,
*$OMP& LNWIN,LLDH,ONE,IPU,ZERO,IPW,N,INDXS)     
                        MYTHREAD = OMP_GET_THREAD_NUM()
*$OMP DO
#endif
                        DO 95 INDX = INDXS, N, NB 
                           IF( DEBUG ) WRITE(*,*)MYROW,MYCOL, 
     $                          'INDX=',INDX
                           CALL INFOG2L( LKTOP, INDX, 
     $                          DESCH, NPROW, NPCOL, MYROW, MYCOL, 
     $                          ILOC, JLOC, RSRC1, CSRC1 )
                           IF( MYROW.EQ.RSRC1.AND.MYCOL.EQ.CSRC1 ) THEN
                              LCOLS = MIN( NB, N-INDX+1 )
                              IF( DEBUG ) WRITE(*,*)MYROW,MYCOL, 
     $                             'LCOLS=',LCOLS
                              CALL DGEMM( 'Transpose', 'No Transpose', 
     $                             LNWIN, LCOLS, LNWIN, ONE, WORK(IPU), 
     $                             LNWIN, H((JLOC-1)*LLDH+ILOC), LLDH, 
     $                             ZERO, WORK(IPW+MYTHREAD*NB**2), 
     $                             LNWIN )
                              CALL DLACPY( 'All', LNWIN, LCOLS, 
     $                             WORK(IPW+MYTHREAD*NB**2), LNWIN, 
     $                             H((JLOC-1)*LLDH+ILOC), LLDH )
                           END IF
 95                     CONTINUE
#ifdef USE_OMP
*$OMP END DO 
*$OMP END PARALLEL
#endif
                     END IF
                  END IF
               ELSEIF( ACCUM .AND. BLK22 ) THEN
                  KS = LNS
*     
*     The LNWIN-by-LNWIN matrix U containing the accumulated orthogonal 
*     transformations has the following structure:
*     
*         [ U11  U12 ]
*     U = [          ],
*         [ U21  U22 ]
*     
*     where U21 is KS-by-KS upper triangular and U12 is
*     (LNWIN-KS)-by-(LNWIN-KS) lower triangular. Here, KS = LNS.
*     
*     Update the columns of H and Z affected by the reordering
*     
*     Compute H2*U21 + H1*U11 in workspace.
*     
                  IF( DIR.EQ.2 .AND. LENCBUF.GT.0 .AND. 
     $                 MYCOL.EQ.LCSRC ) THEN
                     IF( WANTT ) THEN
                        IF( DEBUG ) WRITE(*,*)MYROW,MYCOL, 
     $                       'PDLAQR5: ICEIL 6'
                        DO 100 INDX = 1, LKTOP-LIROFFH-1, NB
                           CALL INFOG2L( INDX, LKTOP, DESCH, NPROW, 
     $                          NPCOL, MYROW, MYCOL, ILOC, JLOC, RSRC1, 
     $                          CSRC1 )
                           IF( MYROW.EQ.RSRC1.AND.MYCOL.EQ.CSRC1 ) THEN
                              JLOC1 = INDXG2L( LKTOP+LNWIN-KS, NB, 
     $                             MYCOL, DESCH( CSRC_ ), NPCOL ) 
                              LROWS = MIN( NB, LKTOP-INDX ) 
                              CALL DLACPY( 'All', LROWS, KS, 
     $                             H((JLOC1-1)*LLDH+ILOC ), LLDH,
     $                             WORK(IPW), LROWS )
                              CALL DTRMM( 'Right', 'Upper', 
     $                             'No transpose','No transpose', LROWS, 
     $                             KS, ONE, WORK( IPU+LNWIN-KS ), LNWIN, 
     $                             WORK(IPW), LROWS )
                              CALL DGEMM('No transpose', 'No transpose', 
     $                             LROWS, KS, LNWIN-KS, ONE, 
     $                             H((JLOC-1)*LLDH+ILOC), LLDH,
     $                             WORK( IPU ), LNWIN, ONE, WORK(IPW),
     $                             LROWS )
*     
*     Compute H1*U12 + H2*U22 in workspace.
*     
                              CALL DLACPY( 'All', LROWS, LNWIN-KS, 
     $                             H((JLOC-1)*LLDH+ILOC), LLDH,
     $                             WORK( IPW+KS*LROWS ), LROWS )
                              CALL DTRMM( 'Right', 'Lower', 
     $                             'No transpose', 'No transpose', 
     $                             LROWS, LNWIN-KS, ONE,
     $                             WORK( IPU+LNWIN*KS ), LNWIN,
     $                             WORK( IPW+KS*LROWS ), LROWS ) 
                              CALL DGEMM('No transpose', 'No transpose', 
     $                             LROWS, LNWIN-KS, KS, ONE, 
     $                             H((JLOC1-1)*LLDH+ILOC), LLDH,
     $                             WORK( IPU+LNWIN*KS+LNWIN-KS ), LNWIN,
     $                             ONE, WORK( IPW+KS*LROWS ), 
     $                             LROWS )
*     
*     Copy workspace to H.
*     
                              CALL DLACPY( 'All', LROWS, LNWIN, 
     $                             WORK(IPW), LROWS, 
     $                             H((JLOC-1)*LLDH+ILOC), LLDH )
                           END IF 
 100                    CONTINUE
                     END IF
*     
                     IF( WANTZ ) THEN
*     
*     Compute Z2*U21 + Z1*U11 in workspace.
*     
                        DO 110 INDX = 1, N, NB
                           CALL INFOG2L( INDX, LKTOP, DESCZ, NPROW, 
     $                          NPCOL, MYROW, MYCOL, ILOC, JLOC, RSRC1, 
     $                          CSRC1 )
                           IF( MYROW.EQ.RSRC1.AND.MYCOL.EQ.CSRC1 ) THEN
                              JLOC1 = INDXG2L( LKTOP+LNWIN-KS, NB, 
     $                             MYCOL, DESCZ( CSRC_ ), NPCOL ) 
                              LROWS = MIN(NB,N-INDX+1)
                              CALL DLACPY( 'All', LROWS, KS, 
     $                             Z((JLOC1-1)*LLDZ+ILOC ), LLDZ,
     $                             WORK(IPW), LROWS )
                              CALL DTRMM( 'Right', 'Upper', 
     $                             'No transpose', 'No transpose', 
     $                             LROWS, KS, ONE, WORK( IPU+LNWIN-KS ), 
     $                             LNWIN, WORK(IPW), LROWS )
                              CALL DGEMM( 'No transpose', 
     $                             'No transpose', LROWS, KS, LNWIN-KS, 
     $                             ONE, Z((JLOC-1)*LLDZ+ILOC), LLDZ,
     $                             WORK( IPU ), LNWIN, ONE, WORK(IPW),
     $                             LROWS )
*     
*     Compute Z1*U12 + Z2*U22 in workspace.
*     
                              CALL DLACPY( 'All', LROWS, LNWIN-KS, 
     $                             Z((JLOC-1)*LLDZ+ILOC), LLDZ,
     $                             WORK( IPW+KS*LROWS ), LROWS)
                              CALL DTRMM( 'Right', 'Lower', 
     $                             'No transpose', 'No transpose', 
     $                             LROWS, LNWIN-KS, ONE,
     $                             WORK( IPU+LNWIN*KS ), LNWIN,
     $                             WORK( IPW+KS*LROWS ), LROWS ) 
                              CALL DGEMM( 'No transpose', 
     $                             'No transpose', LROWS, LNWIN-KS, KS, 
     $                             ONE, Z((JLOC1-1)*LLDZ+ILOC), LLDZ,
     $                             WORK( IPU+LNWIN*KS+LNWIN-KS ), LNWIN,
     $                             ONE, WORK( IPW+KS*LROWS ), 
     $                             LROWS )
*     
*     Copy workspace to Z.
*     
                              CALL DLACPY( 'All', LROWS, LNWIN, 
     $                             WORK(IPW), LROWS, 
     $                             Z((JLOC-1)*LLDZ+ILOC), LLDZ )
                           END IF 
 110                    CONTINUE
                     END IF
                  END IF
*     
                  IF( DIR.EQ.1 .AND. LENRBUF.GT.0 .AND. 
     $                 MYROW.EQ.LRSRC ) THEN
                     IF( WANTT ) THEN
                        IF( DEBUG ) WRITE(*,*)MYROW,MYCOL, 
     $                       'PDLAQR5: ICEIL 7'
                        INDXS = ICEIL(LKBOT,NB)*NB + 1
                        DO 120 INDX = INDXS, N, NB  
                           CALL INFOG2L( LKTOP, INDX, 
     $                          DESCH, NPROW, NPCOL, MYROW, MYCOL, ILOC, 
     $                          JLOC, RSRC1, CSRC1 )
                           IF( MYROW.EQ.RSRC1.AND.MYCOL.EQ.CSRC1 ) THEN
*     
*     Compute U21**T*H2 + U11**T*H1 in workspace.
*     
                              ILOC1 = INDXG2L( LKBOT+LNWIN-KS, NB, 
     $                             MYROW, DESCH( RSRC_ ), NPROW ) 
                              LCOLS = MIN( NB, N-INDX+1 )
                              CALL DLACPY( 'All', KS, LCOLS,
     $                             H((JLOC-1)*LLDH+ILOC1), LLDH, 
     $                             WORK(IPW), LNWIN )
                              CALL DTRMM( 'Left', 'Upper', 'Transpose',
     $                             'No transpose', KS, LCOLS, ONE,
     $                             WORK( IPU+LNWIN-KS ), LNWIN, 
     $                             WORK(IPW), LNWIN )
                              CALL DGEMM( 'Transpose', 'No transpose', 
     $                             KS, LCOLS, LNWIN-KS, ONE, WORK(IPU), 
     $                             LNWIN, H((JLOC-1)*LLDH+ILOC), LLDH, 
     $                             ONE, WORK(IPW), LNWIN )
*     
*     Compute U12**T*H1 + U22**T*H2 in workspace.
*     
                              CALL DLACPY( 'All', LNWIN-KS, LCOLS,
     $                             H((JLOC-1)*LLDH+ILOC), LLDH, 
     $                             WORK( IPW+KS ), LNWIN )
                              CALL DTRMM( 'Left', 'Lower', 'Transpose',
     $                             'No transpose', LNWIN-KS, LCOLS, ONE,
     $                             WORK( IPU+LNWIN*KS ), LNWIN,
     $                             WORK( IPW+KS ), LNWIN )
                              CALL DGEMM( 'Transpose', 'No Transpose', 
     $                             LNWIN-KS, LCOLS, KS, ONE,
     $                             WORK( IPU+LNWIN*KS+LNWIN-KS ), LNWIN,
     $                             H((JLOC-1)*LLDH+ILOC1), LLDH,
     $                             ONE, WORK( IPW+KS ), LNWIN )
*     
*     Copy workspace to T.
*     
                              CALL DLACPY( 'All', LNWIN, LCOLS, 
     $                             WORK(IPW), LNWIN, 
     $                             H((JLOC-1)*LLDH+ILOC), LLDH )
                           END IF
 120                    CONTINUE
                     END IF
                  END IF
               ELSE
*     
*     Use unaccumulated transformations to update far-from-diagonal
*     entries
*     
                  IF( DIR.EQ.2 .AND. LENCBUF.GT.0 .AND. 
     $                 MYCOL.EQ.LCSRC ) THEN
                     IF( WANTT ) THEN
                        IF( DEBUG ) WRITE(*,*)MYROW,MYCOL, 
     $                       'PDLAQR5: ICEIL 8'
                        DO 122 INDX = 1, (ICEIL(LKTOP-1,NB)-1)*NB, NB
                           CALL INFOG2L( INDX, LKTOP, DESCH, NPROW, 
     $                          NPCOL, MYROW, MYCOL, ILOC, JLOC, RSRC1, 
     $                          CSRC1 )
                           IF( MYROW.EQ.RSRC1.AND.MYCOL.EQ.CSRC1 ) THEN
                              LROWS = NB
                              DO 123 IHH = 1, LNWIN 
                                 TAU = WORK(IPU+(IHH-1)*3)
                                 WORK(IPU+(IHH-1)*3) = ONE
                                 CALL DLARFX( 'Right', LROWS, 
     $                                MIN(3,LNWIN-IHH+1),
     $                                WORK(IPU+(IHH-1)*3), TAU, 
     $                                H((JLOC-1+IHH-1)*LLDH+ILOC), 
     $                                LLDH, WORK(IPW) )
                                 WORK(IPU+(IHH-1)*3) = TAU
 123                          CONTINUE
                           END IF
 122                    CONTINUE
                     END IF
                     IF( WANTZ ) THEN
                        DO 124 INDX = 1, N, NB
                           CALL INFOG2L( INDX, LKTOP, DESCZ, NPROW, 
     $                          NPCOL, MYROW, MYCOL, ILOC, JLOC, RSRC1, 
     $                          CSRC1 )
                           IF( MYROW.EQ.RSRC1.AND.MYCOL.EQ.CSRC1 ) THEN
                              LROWS = MIN(NB,N-INDX+1)
                              DO 125 IHH = 1, LNWIN 
                                 TAU = WORK(IPU+(IHH-1)*3)
                                 WORK(IPU+(IHH-1)*3) = ONE
                                 CALL DLARFX ( 'Right', LROWS, 
     $                                MIN(3,LNWIN-IHH+1),
     $                                WORK(IPU+(IHH-1)*3), TAU, 
     $                                Z((JLOC-1+IHH-1)*LLDZ+ILOC), 
     $                                LLDZ, WORK(IPW) )
                                 WORK(IPU+(IHH-1)*3) = TAU
 125                          CONTINUE
                           END IF
 124                    CONTINUE
                     END IF
                  END IF
*     
*     Update the rows of T affected by the bulge-chase
*     
                  IF( DIR.EQ.1 .AND. LENRBUF.GT.0 .AND. 
     $                 MYROW.EQ.LRSRC ) THEN
                     IF( WANTT ) THEN
                        IF( DEBUG ) WRITE(*,*)MYROW,MYCOL, 
     $                       'PDLAQR5: ICEIL 9'
                        INDXS = ICEIL(LKBOT,NB)*NB + 1
                        DO 126 INDX = INDXS, N, NB 
                           CALL INFOG2L( LKTOP, INDX, DESCH, NPROW, 
     $                          NPCOL, MYROW, MYCOL, ILOC, JLOC, RSRC1, 
     $                          CSRC1 )
                           IF( MYROW.EQ.RSRC1.AND.MYCOL.EQ.CSRC1 ) THEN
                              LCOLS = MIN( NB, N-INDX+1 )
                              DO 127 IHH = 1, LNWIN 
                                 TAU = WORK(IPU+(IHH-1)*3)
                                 WORK(IPU+(IHH-1)*3) = ONE
                                 CALL DLARFX ( 'Left', 
     $                                MIN(3,LNWIN-IHH+1), LCOLS,
     $                                WORK(IPU+(IHH-1)*3), TAU, 
     $                                H((JLOC-1)*LLDH+ILOC+IHH-1), 
     $                                LLDH, WORK(IPW) )
                                 WORK(IPU+(IHH-1)*3) = TAU
 127                          CONTINUE
                           END IF
 126                    CONTINUE
                     END IF
                  END IF
               END IF
*     
*     Update position information about current window
*     
               IF( DIR.EQ.2 ) THEN
                  IF( DEBUG ) WRITE(*,*)MYROW,MYCOL, 
     $                 'Old value of LKTOP=',LKTOP
                  IF( DEBUG) WRITE(*,*)MYROW,MYCOL, 
     $                 'Old value of LKBOT=',LKBOT
                  IF( DEBUG) WRITE(*,*)MYROW,MYCOL, 
     $                 'LNWIN,LCHAIN,NB,KBOT=',
     $                 LNWIN,LCHAIN,NB,KBOT
                  IF( LKBOT.EQ.KBOT ) THEN
                     LKTOP = KBOT+1
                     LKBOT = KBOT+1
                     IWORK( 1+(WIN-1)*5 ) = LKTOP
                     IWORK( 2+(WIN-1)*5 ) = LKBOT
                     IWORK(5+(WIN-1)*5) = 2
                     IF( DEBUG) WRITE(*,*)MYROW,MYCOL, 
     $                    'New value of LKTOP=',LKTOP
                     IF( DEBUG) WRITE(*,*)MYROW,MYCOL, 
     $                    'New value of LKBOT=',LKBOT
                  ELSE
                     LKTOP = MIN( LKTOP + LNWIN - LCHAIN, 
     $                    ICEIL( LKTOP, NB )*NB - LCHAIN + 1,
     $                    KBOT )
                     IWORK( 1+(WIN-1)*5 ) = LKTOP
                     IF( DEBUG) WRITE(*,*)MYROW,MYCOL, 
     $                    'New value of LKTOP=',LKTOP
                     LKBOT = MIN( LKBOT + LNWIN - LCHAIN, 
     $                    ICEIL( LKBOT, NB )*NB, KBOT )
                     IWORK( 2+(WIN-1)*5 ) = LKBOT 
                     IF( DEBUG) WRITE(*,*)MYROW,MYCOL, 
     $                    'New value of LKBOT=',LKBOT
                     LNWIN = LKBOT-LKTOP+1 
                     IF( LNWIN.EQ.LCHAIN ) IWORK(5+(WIN-1)*5) = 2
                  END IF
                  IF( DEBUG) WRITE(*,*)MYROW,MYCOL,'SWEEP,LKTOP,LKBOT=',
     $                 SWEEP,LKTOP,LKBOT
c                    IF( SWEEP.EQ.5 .AND. WIN.EQ.3 .AND.
c     $                 LKTOP.EQ.21 .AND. LKBOT.EQ.24 ) RETURN
               END IF
 75            CONTINUE
 70         CONTINUE
 65      CONTINUE
*
*     In debugging, check orthogonality
*
         IF( DEBUG .AND. CHKORTH ) THEN
            IPW = MAX(1,IPW)
            ORTH = PCHKORTH( N, Z, DESCZ, WORK(IPW), LWORK-IPW+1 )
            WRITE(*,*) MYROW,MYCOL,'SWEEP,ORTH=',SWEEP,ORTH
         END IF
*     
*     If bulges was chasen off from first window it is removed
*     
         IF( DEBUG ) WRITE(*,*)MYROW,MYCOL, 'ICHOFF=',ICHOFF
         IF( ICHOFF.GT.0 ) THEN
            DO 128 WIN = 2, ANMWIN
               IWORK( 1+(WIN-2)*5 ) = IWORK( 1+(WIN-1)*5 )
               IWORK( 2+(WIN-2)*5 ) = IWORK( 2+(WIN-1)*5 )
               IWORK( 3+(WIN-2)*5 ) = IWORK( 3+(WIN-1)*5 )
               IWORK( 4+(WIN-2)*5 ) = IWORK( 4+(WIN-1)*5 )
               IWORK( 5+(WIN-2)*5 ) = IWORK( 5+(WIN-1)*5 )
 128        CONTINUE
            ANMWIN = ANMWIN - 1
            IF( DEBUG ) WRITE(*,*)MYROW,MYCOL, 'ANMWIN=',ANMWIN
            IPIW = 6+(ANMWIN-1)*5
         END IF
*     
*     If we have no more windows, return
*     
         IF( DEBUG ) WRITE(*,*)MYROW,MYCOL, 'Return if ANMWIN < 1'
         IF( ANMWIN.LT.1 ) RETURN
*
      ELSE
*     
*     Set up windows such that as many bulges as possible can be moved
*     over the border to the next block. Make sure that the cross border
*     window is at least (NTINY+1)-by-(NTINY+1), unless we are chasing 
*     of the bulges from the last window. This is accomplished by setting 
*     the bottom index LKBOT such that the local window has the correct 
*     size. 
*
*     If LKBOT then becomes larger than KBOT, the endpoint of the whole
*     global submatrix, or LKTOP from a window located already residing 
*     at the other side of the border, this is taken care of by some
*     dirty tricks.
*     
         IF( DEBUG ) WRITE(*,*)MYROW,MYCOL, 
     $        'PDLAQR5: Set up windows for cr. border'
         IF( DEBUG ) WRITE(*,*)MYROW,MYCOL, 'ANMWIN=',ANMWIN
         DO 130 WIN = 1, ANMWIN
            IF( DEBUG ) WRITE(*,*)MYROW,MYCOL, 'Setting up window', WIN
            LKTOP1 = IWORK( 1+(WIN-1)*5 )
            IF( DEBUG ) WRITE(*,*)MYROW,MYCOL, 'LKTOP1=',LKTOP1
c            IF( WIN.GT.1 .AND. MOD(WIN,2).EQ.1 ) THEN
c               IF( DEBUG ) WRITE(*,*)MYROW,MYCOL, 'lktop2 from next window'
c               LKTOP2 = IWORK( 1+(WIN-2)*5 )
c            ELSE
c               IF( DEBUG ) WRITE(*,*)MYROW,MYCOL, 'lktop2=n+1'
c            LKTOP2 = KBOT+1
c            END IF
c            IF( DEBUG ) WRITE(*,*)MYROW,MYCOL, 'LKTOP2=',LKTOP2
            LKBOT = IWORK( 2+(WIN-1)*5 )
            IF( DEBUG ) WRITE(*,*)MYROW,MYCOL, 'LKBOT=',LKBOT
            LNWIN = MAX( 6, MIN( LKBOT - LKTOP1 + 1, LCHAIN ) )
            IF( DEBUG ) WRITE(*,*)MYROW,MYCOL, 'LNWIN=',LNWIN
            LKBOT1 = MAX( MIN( KBOT, ICEIL(LKTOP1,NB)*NB+LCHAIN), 
     $           MIN( KBOT, MIN( LKTOP1+2*LNWIN-1, 
     $           (ICEIL(LKTOP1,NB)+1)*NB ) ) )
            IF( DEBUG ) WRITE(*,*)MYROW,MYCOL, 'LKBOT1=',LKBOT1
            IWORK( 2+(WIN-1)*5 ) = LKBOT1
 130     CONTINUE
         ICHOFF = 0
*
*     Keep a record over what windows that were moved over the borders
*     such that we can delay some windows due to lack of space on the
*     other side of the border; we do not want to leave any of the 
*     bulges behind...
*
*     IWORK( 5+(WIN-1)*5 ) = 0: window WIN has not been processed
*     IWORK( 5+(WIN-1)*5 ) = 1: window WIN is being processed (need to
*                               know for updates)
*     IWORK( 5+(WIN-1)*5 ) = 2: window WIN has been fully processed
*
*     So, start by marking all windows as not processed.
*
         DO 135 WIN = 1, ANMWIN
            IWORK( 5+(WIN-1)*5 ) = 0
 135     CONTINUE
*     
*     Do the cross border bulge-chase as follows: Start from the
*     first window (the one that is closest to be chased off the
*     diagonal of H) and take the odd windows first followed by the 
*     even ones. To not get into hang-problems on processor meshes 
*     with at least one odd dimension, the windows will in such a case
*     be processed in chunks of {the minumum odd process dimension}-1
*     windows to avoid overlapping processor scopes in forming the
*     cross border computational windows and the cross border update
*     regions    
*     
*         ODDNPROW = MOD( NPROW, 2 ).EQ.1
*         ODDNPCOL = MOD( NPCOL, 2 ).EQ.1
*         IF( ODDNPROW .AND. .NOT. ODDNPCOL ) THEN
*            WCHUNK = MAX( NPROW-1, 1 )
*         ELSEIF( .NOT. ODDNPROW .AND. ODDNPCOL ) THEN
*            WCHUNK = MAX( NPCOL-1, 1 )
*         ELSEIF( ODDNPROW .AND. ODDNPCOL ) THEN
*            WCHUNK = MAX( MIN( NPROW, NPCOL ) - 1, 1 )
*         ELSE
*            WCHUNK = ANMWIN
*         END IF
         WCHUNK = MAX( 1, MIN( ANMWIN, NPROW-1, NPCOL-1 ) )
         IF( DEBUG ) WRITE(*,*)MYROW,MYCOL, 'PDLAQR5: ANMWIN,WCHUNK=',
     $        ANMWIN,WCHUNK
         NUMCHUNK = ICEIL( ANMWIN, WCHUNK )
         IF( DEBUG ) WRITE(*,*)MYROW,MYCOL, 
     $        'PDLAQR5: NUMCHUNK=',NUMCHUNK
*     
*     Based on the computed chunk of windows, start working with
*     crossborder bulge-chasing. Repeat this as long as there is
*     still work left to do (137 is a kind of do-while statement)
*     
         IF( DEBUG ) CALL BLACS_BARRIER( ICTXT, 'All' )
         IF( DEBUG ) WRITE(*,*) MYROW,MYCOL,'137 while-loop'
         IF( DEBUG ) CALL BLACS_BARRIER( ICTXT, 'All' )
 137     CONTINUE
         IF( DEBUG ) WRITE(*,*) MYROW,MYCOL,
     $        'Next round - 137 while-loop'
*
*     Zero out LENRBUF and LENCBUF each time we restart this loop
*
         LENRBUF = 0
         LENCBUF = 0
*
         DO 140 ODDEVEN = 1, MIN( 2, ANMWIN )
            IF( DEBUG ) WRITE(*,*) MYROW,MYCOL,'Next round - oddeven'
         DO 150 CHUNKNUM = 1, NUMCHUNK
            IF( DEBUG ) WRITE(*,*) MYROW,MYCOL,'Next round - chunknum'
            IPNEXT = 1
            DO 160 WIN = ODDEVEN+(CHUNKNUM-1)*WCHUNK, 
     $           MIN(ANMWIN,MAX(1,ODDEVEN+(CHUNKNUM)*WCHUNK-1)), 2
               IF( DEBUG ) CALL BLACS_BARRIER( ICTXT, 'All' )
               IF( DEBUG ) WRITE(*,*)MYROW,MYCOL, 'Consider window',WIN
               IF( DEBUG ) CALL BLACS_BARRIER( ICTXT, 'All' )
*     
*     Get position and size of the WIN:th active window and make sure
*     that we skip the cross border bulge for this window if the window
*     is not shared between several data layout blocks (and processors).
*
*     Also, delay windows that do not have sufficient size of the other
*     side of the border. Moreover, make sure to skip windows that was
*     already processed in the last round of the do-while loop (137)
*     
               IF( IWORK( 5+(WIN-1)*5 ).EQ.2 ) GO TO 165
               LKTOP = IWORK( 1+(WIN-1)*5 )
               LKBOT = IWORK( 2+(WIN-1)*5 )
               IF( WIN.GT.1 ) THEN
                  LKTOP2 = IWORK( 1+(WIN-2)*5 )
               ELSE
                  LKTOP2 = KBOT+1
               END IF
               IF( ICEIL(LKTOP,NB).EQ.ICEIL(LKBOT,NB) .OR.
     $             LKBOT.GE.LKTOP2 ) GO TO 165
               LNWIN = LKBOT - LKTOP + 1
               IF( LNWIN.LE.NTINY .AND. LKBOT.NE.KBOT .AND.
     $             .NOT. MOD(LKBOT,NB).EQ.0  ) GO TO 165
*
*     If window is going to be processed, mark it as processed
*
               IWORK( 5+(WIN-1)*5 ) = 1
*     
*     Extract processors for current cross border window, as below:
*     
*               1 | 2
*               -----
*               3 | 4
*     
               RSRC1 = IWORK( 3+(WIN-1)*5 )
               CSRC1 = IWORK( 4+(WIN-1)*5 )
               RSRC2 = RSRC1
               CSRC2 = MOD( CSRC1+1, NPCOL )
               RSRC3 = MOD( RSRC1+1, NPROW )
               CSRC3 = CSRC1
               RSRC4 = MOD( RSRC1+1, NPROW )
               CSRC4 = MOD( CSRC1+1, NPCOL )
*     
               IF( DEBUG ) WRITE(*,*)MYROW,MYCOL, 
     $              'WIN,LKTOP,LKBOT,LNWIN=',
     $              WIN,LKTOP,LKBOT,LNWIN
*     
*     Form group of four processors for cross border window
*     
               IF( ( MYROW.EQ.RSRC1 .AND. MYCOL.EQ.CSRC1 ) .OR.
     $              ( MYROW.EQ.RSRC2 .AND. MYCOL.EQ.CSRC2 ) .OR.
     $              ( MYROW.EQ.RSRC3 .AND. MYCOL.EQ.CSRC3 ) .OR.
     $              ( MYROW.EQ.RSRC4 .AND. MYCOL.EQ.CSRC4 ) ) THEN
                  IF( DEBUG ) WRITE(*,*) MYROW,MYCOL,
     $                 'My process is involved in group of four'
*
*     Compute the upper and lower parts of the active window 
*
                  DIM1 = NB - MOD(LKTOP-1,NB) 
                  DIM4 = LNWIN - DIM1
*
*     Temporarily compute a new value of the size of the computational 
*     window that is larger than or equal to NTINY+1; call the *real*
*     value DIM
*
                  DIM = LNWIN
                  LNWIN = MAX(NTINY+1,LNWIN)
*     
*     Divide workspace 
*     
                  IPU = IPNEXT
                  IF( .NOT. ACCUM ) THEN
                     IPH = IPU + 3*DIM
                  ELSE
                     IPH = IPU + DIM**2
                  END IF
                  IPUU = IPH + LNWIN**2
                  IPV = IPUU + LNWIN**2
                  IPNEXT = IPH                  
                  IF( DEBUG ) WRITE(*,*)MYROW,MYCOL, 'IPU,IPH,IPNEXT=',
     $                 IPU,IPH,IPNEXT
                  IF( DIM.LT.LNWIN ) THEN
                     CALL DLASET( 'All', LNWIN, LNWIN, ZERO, 
     $                    ONE, WORK( IPH ), LNWIN )
                  ELSE
                     CALL DLASET( 'All', DIM, DIM, ZERO, 
     $                    ZERO, WORK( IPH ), LNWIN )
                  END IF
*     
*     Form the active window 
*     
                  IF( DEBUG ) WRITE(*,*)MYROW,MYCOL,
     $                 'Form active window'
                  IF( MYROW.EQ.RSRC1 .AND. MYCOL.EQ.CSRC1 ) THEN
                     ILOC = INDXG2L( LKTOP, NB, MYROW, 
     $                    DESCH( RSRC_ ), NPROW )
                     JLOC = INDXG2L( LKTOP, NB, MYCOL, 
     $                    DESCH( CSRC_ ), NPCOL )
                     IF( DEBUG ) WRITE(*,*)MYROW,MYCOL, 'ILOC,JLOC=',
     $                    ILOC,JLOC
                     CALL DLACPY( 'All', DIM1, DIM1, 
     $                    H((JLOC-1)*LLDH+ILOC), LLDH, WORK(IPH),
     $                    LNWIN )
                     IF( RSRC1.NE.RSRC4 .OR. CSRC1.NE.CSRC4 ) THEN
                        IF( DEBUG ) WRITE(*,*)MYROW,MYCOL,
     $                       'Proc. 1 sending to proc. 4'
                        CALL DGESD2D( ICTXT, DIM1, DIM1, 
     $                       WORK(IPH), LNWIN, RSRC4, CSRC4 )
                        IF( DEBUG ) WRITE(*,*)MYROW,MYCOL,
     $                       'Proc. 1 receiving from proc. 4'
                        CALL DGERV2D( ICTXT, DIM4, DIM4, 
     $                       WORK(IPH+DIM1*LNWIN+DIM1), 
     $                       LNWIN, RSRC4, CSRC4 )
                     END IF
                  END IF
                  IF( MYROW.EQ.RSRC4 .AND. MYCOL.EQ.CSRC4 ) THEN
                     ILOC = INDXG2L( LKTOP+DIM1, NB, MYROW, 
     $                    DESCH( RSRC_ ), NPROW )
                     JLOC = INDXG2L( LKTOP+DIM1, NB, MYCOL, 
     $                    DESCH( CSRC_ ), NPCOL )
                     IF( DEBUG ) WRITE(*,*)MYROW,MYCOL, 
     $                    'IPH,DIM1,LNWIN=',
     $                    IPH, DIM1, LNWIN
                     CALL DLACPY( 'All', DIM4, DIM4, 
     $                    H((JLOC-1)*LLDH+ILOC), LLDH, 
     $                    WORK(IPH+DIM1*LNWIN+DIM1), 
     $                    LNWIN )
                     IF( RSRC4.NE.RSRC1 .OR. CSRC4.NE.CSRC1 ) THEN
                        IF( DEBUG ) WRITE(*,*)MYROW,MYCOL,
     $                       'Proc. 4 sending to proc. 1'
                        CALL DGESD2D( ICTXT, DIM4, DIM4, 
     $                       WORK(IPH+DIM1*LNWIN+DIM1), 
     $                       LNWIN, RSRC1, CSRC1 ) 
                        IF( DEBUG ) WRITE(*,*)MYROW,MYCOL,
     $                       'Proc. 4 receiving from proc. 1'
                        CALL DGERV2D( ICTXT, DIM1, DIM1, 
     $                       WORK(IPH), LNWIN, RSRC1, CSRC1 )
                     END IF
                  END IF
                  IF( MYROW.EQ.RSRC2 .AND. MYCOL.EQ.CSRC2 ) THEN
                     ILOC = INDXG2L( LKTOP, NB, MYROW, 
     $                    DESCH( RSRC_ ), NPROW )
                     JLOC = INDXG2L( LKTOP+DIM1, NB, MYCOL, 
     $                    DESCH( CSRC_ ), NPCOL )
                     CALL DLACPY( 'All', DIM1, DIM4, 
     $                    H((JLOC-1)*LLDH+ILOC), LLDH, 
     $                    WORK(IPH+DIM1*LNWIN), 
     $                    LNWIN )
                     IF( RSRC2.NE.RSRC1 .OR. CSRC2.NE.CSRC1 ) THEN
                        IF( DEBUG ) WRITE(*,*)MYROW,MYCOL,
     $                       'Proc. 2 sending to proc. 1'
                        CALL DGESD2D( ICTXT, DIM1, DIM4, 
     $                       WORK(IPH+DIM1*LNWIN), 
     $                       LNWIN, RSRC1, CSRC1 )
                     END IF
                  END IF
                  IF( MYROW.EQ.RSRC2 .AND. MYCOL.EQ.CSRC2 ) THEN 
                     IF( RSRC2.NE.RSRC4 .OR. CSRC2.NE.CSRC4 ) THEN
                        IF( DEBUG ) WRITE(*,*)MYROW,MYCOL,
     $                       'Proc. 2 sending to proc. 4'
                        CALL DGESD2D( ICTXT, DIM1, DIM4, 
     $                       WORK(IPH+DIM1*LNWIN), 
     $                       LNWIN, RSRC4, CSRC4 )
                     END IF
                  END IF
                  IF( MYROW.EQ.RSRC3 .AND. MYCOL.EQ.CSRC3 ) THEN
                     ILOC = INDXG2L( LKTOP+DIM1, NB, MYROW, 
     $                    DESCH( RSRC_ ), NPROW )
                     JLOC = INDXG2L( LKTOP+DIM1-1, NB, MYCOL, 
     $                    DESCH( CSRC_ ), NPCOL )
                     CALL DLACPY( 'All', 1, 1, 
     $                    H((JLOC-1)*LLDH+ILOC), LLDH, 
     $                    WORK(IPH+(DIM1-1)*LNWIN+DIM1), 
     $                    LNWIN )
                     IF( RSRC3.NE.RSRC1 .OR. CSRC3.NE.CSRC1 ) THEN 
                        IF( DEBUG ) WRITE(*,*)MYROW,MYCOL,
     $                       'Proc. 3 sending to proc. 1'
                        CALL DGESD2D( ICTXT, 1, 1, 
     $                       WORK(IPH+(DIM1-1)*LNWIN+DIM1), 
     $                       LNWIN, RSRC1, CSRC1 )
                     END IF
                  END IF
                  IF( MYROW.EQ.RSRC3 .AND. MYCOL.EQ.CSRC3 ) THEN
                     IF( RSRC3.NE.RSRC4 .OR. CSRC3.NE.CSRC4 ) THEN 
                        IF( DEBUG ) WRITE(*,*)MYROW,MYCOL,
     $                       'Proc. 3 sending to proc. 4'
                        CALL DGESD2D( ICTXT, 1, 1, 
     $                       WORK(IPH+(DIM1-1)*LNWIN+DIM1), 
     $                       LNWIN, RSRC4, CSRC4 )
                     END IF
                  END IF
                  IF( MYROW.EQ.RSRC1 .AND. MYCOL.EQ.CSRC1 ) THEN
                     IF( RSRC1.NE.RSRC2 .OR. CSRC1.NE.CSRC2 ) THEN
                        IF( DEBUG ) WRITE(*,*)MYROW,MYCOL,
     $                       'Proc. 1 receiving from proc. 2'
                        CALL DGERV2D( ICTXT, DIM1, DIM4,
     $                       WORK(IPH+DIM1*LNWIN), 
     $                       LNWIN, RSRC2, CSRC2 )
                     END IF
                     IF( RSRC1.NE.RSRC3 .OR. CSRC1.NE.CSRC3 ) THEN
                        IF( DEBUG ) WRITE(*,*)MYROW,MYCOL,
     $                       'Proc. 1 receiving from proc. 3'
                        CALL DGERV2D( ICTXT, 1, 1,
     $                       WORK(IPH+(DIM1-1)*LNWIN+DIM1), 
     $                       LNWIN, RSRC3, CSRC3 )
                     END IF
                  END IF
                  IF( MYROW.EQ.RSRC4 .AND. MYCOL.EQ.CSRC4 ) THEN
                     IF( RSRC4.NE.RSRC2 .OR. CSRC4.NE.CSRC2 ) THEN
                        IF( DEBUG ) WRITE(*,*)MYROW,MYCOL,
     $                       'Proc. 4 receiving from proc. 2'
                        CALL DGERV2D( ICTXT, DIM1, DIM4,
     $                       WORK(IPH+DIM1*LNWIN), 
     $                       LNWIN, RSRC2, CSRC2 )
                     END IF
                     IF( RSRC4.NE.RSRC3 .OR. CSRC4.NE.CSRC3 ) THEN
                        IF( DEBUG ) WRITE(*,*)MYROW,MYCOL,
     $                       'Proc. 4 receiving from proc. 3'
                        CALL DGERV2D( ICTXT, 1, 1,
     $                       WORK(IPH+(DIM1-1)*LNWIN+DIM1), 
     $                       LNWIN, RSRC3, CSRC3 )
                     END IF
                  END IF
                  IF( DEBUG ) WRITE(*,*)MYROW,MYCOL,
     $                 'Done forming active window'
*     
*     Prepare for call to DLAQR6 - it could happen that no bulges
*     where introduced in the pre-cross border step since the chain
*     was to long to fit in the top-left part of the cross border
*     window. In such a case, the bulges are introduced here instead.
*     It could also happen that the bottom-right part is to small to
*     hold the whole chain - in such a case, the bulges are chasen
*     off immediately, as well.
*     
                  IF( (MYROW.EQ.RSRC1 .AND. MYCOL.EQ.CSRC1) .OR.
     $                 (MYROW.EQ.RSRC4 .AND. MYCOL.EQ.CSRC4) ) THEN
                     IF( LKTOP.EQ.KTOP .AND. LKBOT.EQ.KBOT .AND.
     $                   (DIM1.LE.LCHAIN .OR. DIM1.LE.NTINY ) ) THEN
                        JOB = 'All steps'
                        ICHOFF = 1
                     ELSEIF( LKTOP.EQ.KTOP .AND. 
     $                       ( DIM1.LE.LCHAIN .OR. DIM1.LE.NTINY ) ) 
     $                       THEN
                        JOB = 'Introduce and chase'
                     ELSEIF( LKBOT.EQ.KBOT ) THEN
                        JOB = 'Off-chase bulges'
                        ICHOFF = 1
                     ELSE
                        JOB = 'Chase bulges'
                     END IF
                     IF( DEBUG ) WRITE(*,*)MYROW,MYCOL, 'JOB=',JOB
                     KU = LNWIN - KDU + 1
                     KWH = KDU + 1
                     NHO = ( LNWIN-KDU+1-4 ) - ( KDU+1 ) + 1
                     KWV = KDU + 4
                     NVE = LNWIN - KDU - KWV + 1
                     IF( ACCUM )
     $                    CALL DLASET( 'All', LNWIN, LNWIN, 
     $                    ZERO, ONE, WORK(IPUU), LNWIN )
*     
*     ==== Small-bulge multi-shift QR sweep ====
*     
                     LKS = MAX(1, NS - WIN*LNS + 1)
                     IF( DEBUG ) WRITE(*,*)MYROW,MYCOL, 
     $                    'DIM,LNWIN,LKS=',DIM,LNWIN,LKS
                     IF( DEBUG ) WRITE(*,*)MYROW,MYCOL, 
     $                    'WIN,NS,LNS,LKS=',WIN,NS,LNS,LKS
                     IF( DEBUG ) WRITE(*,*)MYROW,MYCOL, 
     $                    'PDLAQR5 -> DLAQR6'
                     IF( PRINT.AND.MYROW.EQ.0.AND.MYCOL.EQ.0 ) THEN
                        CALL DLAPRNT(LNWIN,LNWIN,WORK(IPH),1,1,
     $                       LNWIN,'HH3',6)
                        CALL DLAPRNT(LNWIN,LNWIN,WORK(IPUU),1,1,LNWIN,
     $                       'UU3',6)
C                     CALL DLAPRNT(3,LNWIN,WORK(IPV),1,1,3,'VV3',6)
                     END IF
                     IF( PRINT.AND.MYROW.EQ.1.AND.MYCOL.EQ.1 ) THEN
                        CALL DLAPRNT(LNWIN,LNWIN,WORK(IPH),1,1,
     $                       LNWIN,'HHH3',6)
                        CALL DLAPRNT(LNWIN,LNWIN,WORK(IPUU),1,1,LNWIN,
     $                       'UUU3',6)
C                     CALL DLAPRNT(3,LNWIN,WORK(IPV),1,1,3,'VV3',6)
                     END IF
                     IF( DEBUG ) 
     $                    WRITE(*,*)MYROW,MYCOL,'SR=',SR(LKS:LKS+LNS-1)
                     IF( DEBUG ) 
     $                    WRITE(*,*)MYROW,MYCOL,'SI=',SI(LKS:LKS+LNS-1)
                     STAMP = MPI_WTIME()
                     CALL DLAQR6( JOB, WANTT, ACCUM, LKACC22, LNWIN, 
     $                    1, DIM, LNS, SR( LKS ), SI( LKS ), 
     $                    WORK(IPH), LNWIN, 1, DIM, 
     $                    WORK(IPUU), LNWIN, WORK(IPV), 3, 
     $                    IWORK(IPIW), 2, WORK(IPU), 3, 
     $                    WORK( IPH+KU-1 ), LNWIN, NVE, 
     $                    WORK( IPH+KWV-1 ), LNWIN, NHO, 
     $                    WORK( IPH-1+KU+(KWH-1)*LNWIN ), 
     $                    LNWIN )
                     TIMINGS( 6 ) = TIMINGS( 6 ) + MPI_WTIME() - STAMP
                     IF( PRINT.AND.MYROW.EQ.0.AND.MYCOL.EQ.0 ) THEN
                        CALL DLAPRNT(LNWIN,LNWIN,WORK(IPH),1,1,
     $                       LNWIN,'HH4',6)
                        CALL DLAPRNT(LNWIN,LNWIN,WORK(IPUU),1,1,LNWIN,
     $                       'UU4',6)
C               CALL DLAPRNT(3,LNWIN,WORK(IPV),1,1,3,'VV4',6)      
                     END IF
                     IF( PRINT.AND.MYROW.EQ.1.AND.MYCOL.EQ.1 ) THEN
                        CALL DLAPRNT(LNWIN,LNWIN,WORK(IPH),1,1,
     $                       LNWIN,'HHH4',6)
                        CALL DLAPRNT(LNWIN,LNWIN,WORK(IPUU),1,1,LNWIN,
     $                       'UUU4',6)
C               CALL DLAPRNT(3,LNWIN,WORK(IPV),1,1,3,'VV4',6)      
                     END IF
                     IF( DEBUG ) WRITE(*,*)MYROW,MYCOL, 
     $                    'PDLAQR5 <- DLAQR6' 
*     
*     Copy local submatrices of H back to global matrix 
*     
                     IF( DEBUG ) WRITE(*,*)MYROW,MYCOL, 
     $                    'Copy local submatrices of H back'
                     IF( MYROW.EQ.RSRC1 .AND. MYCOL.EQ.CSRC1 ) THEN
                        IF( DEBUG ) WRITE(*,*)MYROW,MYCOL, 
     $                       'Copy H back 1'
                        ILOC = INDXG2L( LKTOP, NB, MYROW, 
     $                       DESCH( RSRC_ ), NPROW )
                        JLOC = INDXG2L( LKTOP, NB, MYCOL, 
     $                       DESCH( CSRC_ ), NPCOL ) 
                        CALL DLACPY( 'All', DIM1, DIM1, WORK(IPH),
     $                       LNWIN, H((JLOC-1)*LLDH+ILOC), 
     $                       LLDH )
                     END IF
                     IF( MYROW.EQ.RSRC4 .AND. MYCOL.EQ.CSRC4 ) THEN
                        IF( DEBUG ) WRITE(*,*)MYROW,MYCOL, 
     $                       'Copy H back 4'
                        ILOC = INDXG2L( LKTOP+DIM1, NB, MYROW, 
     $                       DESCH( RSRC_ ), NPROW )
                        JLOC = INDXG2L( LKTOP+DIM1, NB, MYCOL, 
     $                       DESCH( CSRC_ ), NPCOL ) 
                        CALL DLACPY( 'All', DIM4, DIM4, 
     $                       WORK(IPH+DIM1*LNWIN+DIM1), 
     $                       LNWIN, H((JLOC-1)*LLDH+ILOC), 
     $                       LLDH )
                     END IF
*
*     Copy actual submatrix of U to the correct place of the buffer
*
                     IF( ACCUM )
     $                    CALL DLACPY( 'All', DIM, DIM, 
     $                    WORK(IPUU), LNWIN, WORK(IPU), DIM )
                  END IF
*     
*     Return data to process 2 and 3
*     
                  IF( DEBUG ) WRITE(*,*)MYROW,MYCOL, 
     $                 'Return data to process 2 and 3'
                  RWS3 = MIN(3,DIM4)
                  CLS3 = MIN(3,DIM1)
                  IF( DEBUG ) WRITE(*,*) 'RWS3,CLS3=',
     $                 RWS3,CLS3
                  IF( MYROW.EQ.RSRC1 .AND. MYCOL.EQ.CSRC1 ) THEN
                     IF( RSRC1.NE.RSRC3 .OR. CSRC1.NE.CSRC3 ) THEN
                        IF( DEBUG ) WRITE(*,*)MYROW,MYCOL,
     $                       'Proc. 1 sending to proc. 3'
                        CALL DGESD2D( ICTXT, RWS3, CLS3, 
     $                       WORK( IPH+(DIM1-CLS3)*LNWIN+DIM1 ), 
     $                       LNWIN, RSRC3, CSRC3 )
                     END IF
                  END IF
                  IF( MYROW.EQ.RSRC4 .AND. MYCOL.EQ.CSRC4 ) THEN
                     IF( RSRC4.NE.RSRC2 .OR. CSRC4.NE.CSRC2 ) THEN
                        IF( DEBUG ) WRITE(*,*)MYROW,MYCOL,
     $                       'Proc. 4 sending to proc. 2'
                        CALL DGESD2D( ICTXT, DIM1, DIM4, 
     $                       WORK( IPH+DIM1*LNWIN), 
     $                       LNWIN, RSRC2, CSRC2 ) 
                     END IF
                  END IF
                  IF( MYROW.EQ.RSRC2 .AND. MYCOL.EQ.CSRC2 ) THEN
                     ILOC = INDXG2L( LKTOP, NB, MYROW, 
     $                    DESCH( RSRC_ ), NPROW )
                     JLOC = INDXG2L( LKTOP+DIM1, NB, MYCOL, 
     $                    DESCH( CSRC_ ), NPCOL )
                     IF( RSRC2.NE.RSRC4 .OR. CSRC2.NE.CSRC4 ) THEN
                        IF( DEBUG ) WRITE(*,*)MYROW,MYCOL,
     $                       'Proc. 2 receiving from proc. 4'
                        CALL DGERV2D( ICTXT, DIM1, DIM4, 
     $                       WORK(IPH+DIM1*LNWIN), 
     $                       LNWIN, RSRC4, CSRC4 )
                        IF( DEBUG ) WRITE(*,*)MYROW,MYCOL,
     $                       'Done!'
                     END IF
                     CALL DLACPY( 'All', DIM1, DIM4, 
     $                    WORK( IPH+DIM1*LNWIN ), 
     $                    LNWIN,
     $                    H((JLOC-1)*LLDH+ILOC), LLDH )
                  END IF
                  IF( MYROW.EQ.RSRC3 .AND. MYCOL.EQ.CSRC3 ) THEN
                     ILOC = INDXG2L( LKTOP+DIM1, NB, MYROW, 
     $                    DESCH( RSRC_ ), NPROW )
                     JLOC = INDXG2L( LKTOP+DIM1-CLS3, NB, MYCOL, 
     $                    DESCH( CSRC_ ), NPCOL )
                     IF( RSRC3.NE.RSRC1 .OR. CSRC3.NE.CSRC1 ) THEN
                        IF( DEBUG ) WRITE(*,*)MYROW,MYCOL,
     $                       'Proc. 3 receiving from proc. 1'
                        CALL DGERV2D( ICTXT, RWS3, CLS3, 
     $                       WORK( IPH+(DIM1-CLS3)*LNWIN+DIM1 ), 
     $                       LNWIN, RSRC1, CSRC1 )
                        IF( DEBUG ) WRITE(*,*)MYROW,MYCOL,
     $                       'Done!'
                     END IF
                     CALL DLACPY( 'Upper', RWS3, CLS3, 
     $                    WORK( IPH+(DIM1-CLS3)*LNWIN+DIM1 ), 
     $                    LNWIN, H((JLOC-1)*LLDH+ILOC), 
     $                    LLDH )
                     IF( RWS3.GT.1 .AND. CLS3.GT.1 ) THEN
                        ELEM = WORK( IPH+(DIM1-CLS3)*LNWIN+DIM1+1 )
                        IF( ELEM.NE.ZERO ) THEN
                           CALL DLACPY( 'Lower', RWS3-1, CLS3-1, 
     $                          WORK( IPH+(DIM1-CLS3)*LNWIN+DIM1+1 ), 
     $                          LNWIN, H((JLOC-1)*LLDH+ILOC+1), 
     $                          LLDH )
                        END IF
                     END IF
                  END IF
*
*     Restore correct value of LNWIN
*
                  LNWIN = DIM
*
               END IF
*     
*     Increment counter for buffers of orthogonal transformations
*     
               IF( DEBUG ) WRITE(*,*)MYROW,MYCOL, 
     $              'Increment counter for buffers'
               IF( MYROW.EQ.RSRC1 .OR. MYCOL.EQ.CSRC1 .OR.
     $              MYROW.EQ.RSRC4 .OR. MYCOL.EQ.CSRC4 ) THEN
                  IF( .NOT. ACCUM ) THEN
                     IF( MYROW.EQ.RSRC1 .OR. MYROW.EQ.RSRC4 ) 
     $                    LENRBUF = LENRBUF + 3*LNWIN
                     IF( MYCOL.EQ.CSRC1 .OR. MYCOL.EQ.CSRC4 ) 
     $                    LENCBUF = LENCBUF + 3*LNWIN
                  ELSE
                     IF( MYROW.EQ.RSRC1 .OR. MYROW.EQ.RSRC4 ) 
     $                    LENRBUF = LENRBUF + LNWIN*LNWIN
                     IF( MYCOL.EQ.CSRC1 .OR. MYCOL.EQ.CSRC4 ) 
     $                    LENCBUF = LENCBUF + LNWIN*LNWIN
                  END IF
                  IF( DEBUG ) WRITE(*,*)MYROW,MYCOL, 
     $                 'LENRBUF,LENCBUF=',LENRBUF,LENCBUF
               END IF
*     
*     If no cross border reordering was performed for the current WIN:th
*     window, the processor jump to this point and consider the next one
*     
 165           CONTINUE
*     
 160        CONTINUE
*
*     In debugging, check for inconsistency in buffer lengths
*
            IF( DEBUG ) THEN
               MNRBUF = LENRBUF
               MXRBUF = LENRBUF
               MNCBUF = LENCBUF
               MXCBUF = LENCBUF
               IF( NPROCS.GT.1 )
     $              CALL IGAMN2D( ICTXT, 'Row', '1-Tree', 1, 1, MNRBUF, 
     $              1, -1, -1, -1, -1, -1 )
               IF( NPROCS.GT.1 )
     $              CALL IGAMX2D( ICTXT, 'Row', '1-Tree', 1, 1, MXRBUF, 
     $              1, -1, -1, -1, -1, -1 )
               IF( NPROCS.GT.1 )
     $              CALL IGAMN2D( ICTXT, 'Col', '1-Tree', 1, 1, MNCBUF, 
     $              1, -1, -1, -1, -1, -1 )
               IF( NPROCS.GT.1 )
     $              CALL IGAMX2D( ICTXT, 'Col', '1-Tree', 1, 1, MXCBUF, 
     $              1, -1, -1, -1, -1, -1 )
               IF( MNRBUF.NE.MXRBUF ) THEN
                  IF( MYROW.EQ.MYCOL ) THEN
                     WRITE(*,*) '*** Error before cross BC ***'
                     WRITE(*,*) 
     $                    '*** Inconsistency in row buffer lengths:'
                     WRITE(*,*) 'Min,Max,Corr,row=',MNRBUF,MXRBUF,
     $                    LENRBUF,MYROW
                  END IF
                  STOP
               END IF
               IF( MNCBUF.NE.MXCBUF ) THEN
                  IF( MYROW.EQ.MYCOL ) THEN
                     WRITE(*,*) '*** Error before cross BC ***'
                     WRITE(*,*) 
     $                    '*** Inconsistency in column buffer lengths:'
                     WRITE(*,*) 'Min,Max,Corr,col=',MNCBUF,MXCBUF,
     $                    LENCBUF,MYCOL
                  END IF
                  STOP
               END IF
            END IF
*     
*     Broadcast orthogonal transformations - this will only happen if
*     the buffer associated with the orthogonal transformations is not
*     empty (controlled by LENRBUF, for row-wise broadcasts, and LENCBUF,
*     for column-wise broadcasts).
*     
            IF( DEBUG ) WRITE(*,*) MYROW,MYCOL,'BC crossborder...'
            DO 170 DIR = 1, 2
               BCDONE = .FALSE.
               IF( DEBUG ) WRITE(*,*) MYROW,MYCOL,'DIR=',DIR
               DO 180 WIN = ODDEVEN+(CHUNKNUM-1)*WCHUNK, 
     $              MIN(ANMWIN,MAX(1,ODDEVEN+(CHUNKNUM)*WCHUNK-1)), 2
                  IF( DEBUG ) WRITE(*,*)MYROW,MYCOL, 
     $                 'Status: WIN,mark=',WIN,
     $                 IWORK( 5+(WIN-1)*5 ) 
                  IF( ( LENRBUF.EQ.0 .AND. LENCBUF.EQ.0 ) .OR.
     $                 BCDONE ) GO TO 185
                  RSRC1 = IWORK( 3+(WIN-1)*5 )
                  CSRC1 = IWORK( 4+(WIN-1)*5 )
                  RSRC4 = MOD( RSRC1+1, NPROW )
                  CSRC4 = MOD( CSRC1+1, NPCOL )
                  IF( DEBUG ) WRITE(*,*) 
     $                 MYROW,MYCOL,'RSRC1,CSRC1,RSRC4,CSRC1=',
     $                 RSRC1,CSRC1,RSRC4,CSRC1
                  IF( ( MYROW.EQ.RSRC1 .AND. MYCOL.EQ.CSRC1 ) .OR.
     $                 ( MYROW.EQ.RSRC4 .AND. MYCOL.EQ.CSRC4 ) ) THEN
                     IF( DIR.EQ.1 .AND. LENRBUF.GT.0 .AND. 
     $                    NPCOL.GT.1 .AND. NPROCS.GT.2 ) THEN
                        IF( DEBUG ) WRITE(*,*) MYROW,MYCOL,'BC row'
                        IF( MYROW.EQ.RSRC1 .OR. ( MYROW.EQ.RSRC4
     $                       .AND. RSRC4.NE.RSRC1 ) ) THEN
                           CALL DGEBS2D( ICTXT, 'Row', '1-Tree', 
     $                          LENRBUF, 1, WORK, LENRBUF )
                        ELSE
                           CALL DGEBR2D( ICTXT, 'Row', '1-Tree', 
     $                          LENRBUF, 1, WORK, LENRBUF, RSRC1, 
     $                          CSRC1 )
                        END IF
                     ELSEIF( DIR.EQ.2 .AND. LENCBUF.GT.0 .AND. 
     $                       NPROW.GT.1 .AND. NPROCS.GT.2 ) THEN
                        IF( DEBUG ) WRITE(*,*) MYROW,MYCOL,'BC col'
                        IF( MYCOL.EQ.CSRC1 .OR. ( MYCOL.EQ.CSRC4
     $                       .AND. CSRC4.NE.CSRC1 ) ) THEN
                           CALL DGEBS2D( ICTXT, 'Col', '1-Tree', 
     $                          LENCBUF, 1, WORK, LENCBUF )
                        ELSE
                           CALL DGEBR2D( ICTXT, 'Col', '1-Tree', 
     $                          LENCBUF, 1, WORK(1+LENRBUF), LENCBUF, 
     $                          RSRC1, CSRC1 )
                        END IF
                     END IF
                     IF( DEBUG ) WRITE(*,*) 
     $                    MYROW,MYCOL,'Copy rbuf into cbuf:',
     $                    LENRBUF
                     IF( LENRBUF.GT.0 .AND. ( MYCOL.EQ.CSRC1 .OR. 
     $                    ( MYCOL.EQ.CSRC4 .AND. CSRC4.NE.CSRC1 ) ) )  
     $                    CALL DLACPY( 'All', LENRBUF, 1, WORK, LENRBUF, 
     $                    WORK(1+LENRBUF), LENCBUF )
                     BCDONE = .TRUE.
                  ELSEIF( MYROW.EQ.RSRC1 .AND. DIR.EQ.1 ) THEN
                     IF( DEBUG ) WRITE(*,*) MYROW,MYCOL,'BC row recv.'
                     IF( LENRBUF.GT.0 .AND. NPCOL.GT.1 )
     $                    CALL DGEBR2D( ICTXT, 'Row', '1-Tree', LENRBUF, 
     $                    1, WORK, LENRBUF, RSRC1, CSRC1 )
                     BCDONE = .TRUE.
                  ELSEIF( MYCOL.EQ.CSRC1 .AND. DIR.EQ.2 ) THEN
                     IF( DEBUG ) WRITE(*,*) MYROW,MYCOL,'BC col recv.'
                     IF( LENCBUF.GT.0 .AND. NPROW.GT.1 )
     $                    CALL DGEBR2D( ICTXT, 'Col', '1-Tree', LENCBUF, 
     $                    1, WORK(1+LENRBUF), LENCBUF, RSRC1, CSRC1 )
                     BCDONE = .TRUE.
                  ELSEIF( MYROW.EQ.RSRC4 .AND. DIR.EQ.1 ) THEN
                     IF( DEBUG ) WRITE(*,*) MYROW,MYCOL,'BC row recv.'
                     IF( LENRBUF.GT.0 .AND. NPCOL.GT.1 )
     $                    CALL DGEBR2D( ICTXT, 'Row', '1-Tree', LENRBUF, 
     $                    1, WORK, LENRBUF, RSRC4, CSRC4 )
                     BCDONE = .TRUE.
                  ELSEIF( MYCOL.EQ.CSRC4 .AND. DIR.EQ.2 ) THEN
                     IF( DEBUG ) WRITE(*,*) MYROW,MYCOL,'BC col recv.'
                     IF( LENCBUF.GT.0 .AND. NPROW.GT.1 )
     $                    CALL DGEBR2D( ICTXT, 'Col', '1-Tree', LENCBUF, 
     $                    1, WORK(1+LENRBUF), LENCBUF, RSRC4, CSRC4 ) 
                     BCDONE = .TRUE.
                  END IF
 185              CONTINUE
 180           CONTINUE
 170        CONTINUE
            IF( DEBUG .AND. PRINT ) THEN
               CALL BLACS_BARRIER( ICTXT, 'All' )
               IF( MYROW.EQ.0 .AND. MYCOL.EQ.0 ) THEN
                  CALL DLAPRNT( LENRBUF, 1, WORK, 1, 1, LENRBUF, 'R00', 
     $                 6 )
                  CALL DLAPRNT( LENCBUF, 1, WORK(1+LENRBUF), 1, 1, 
     $                 LENCBUF, 'C00', 6 )
               END IF
               IF( MYROW.EQ.0 .AND. MYCOL.EQ.1 ) THEN
                  CALL DLAPRNT( LENRBUF, 1, WORK, 1, 1, LENRBUF, 'R01', 
     $                 6 )
                  CALL DLAPRNT( LENCBUF, 1, WORK(1+LENRBUF), 1, 1, 
     $                 LENCBUF, 'C01', 6 )
               END IF
               IF( MYROW.EQ.1 .AND. MYCOL.EQ.0 ) THEN
                  CALL DLAPRNT( LENRBUF, 1, WORK, 1, 1, LENRBUF, 'R10', 
     $                 6 )
                  CALL DLAPRNT( LENCBUF, 1, WORK(1+LENRBUF), 1, 1, 
     $                 LENCBUF, 'C10', 6 )
               END IF
               IF( MYROW.EQ.1 .AND. MYCOL.EQ.1 ) THEN
                  CALL DLAPRNT( LENRBUF, 1, WORK, 1, 1, LENRBUF, 'R11', 
     $                 6 )
                  CALL DLAPRNT( LENCBUF, 1, WORK(1+LENRBUF), 1, 1, 
     $                 LENCBUF, 'C11', 6 )
               END IF
               CALL BLACS_BARRIER( ICTXT, 'All' )
            END IF
*     
*     Prepare for computing cross border updates by exchanging data in 
*     cross border update regions in H and Z 
*     
            DO 190 DIR = 1, 2
               IF( DEBUG ) WRITE(*,*) MYROW,MYCOL,'DIR=',DIR
               WINID = 0
               IPW3 = 1
               DO 200 WIN = ODDEVEN+(CHUNKNUM-1)*WCHUNK, 
     $              MIN(ANMWIN,MAX(1,ODDEVEN+(CHUNKNUM)*WCHUNK-1)), 2
                  IF( IWORK( 5+(WIN-1)*5 ).NE.1 ) GO TO 205
*     
*     Make sure this part of the code is only executed when there has
*     been some work performed on the WIN:th window
*     
                  LKTOP = IWORK( 1+(WIN-1)*5 )
                  LKBOT = IWORK( 2+(WIN-1)*5 )
*     
*     Extract processor indices associated with the current window 
*     
                  RSRC1 = IWORK( 3+(WIN-1)*5 )
                  CSRC1 = IWORK( 4+(WIN-1)*5 )
                  RSRC4 = MOD( RSRC1+1, NPROW )
                  CSRC4 = MOD( CSRC1+1, NPCOL )
*     
*     Compute local number of rows and columns of H and Z to exchange
*     
                  IF(((MYCOL.EQ.CSRC1.OR.MYCOL.EQ.CSRC4).AND.DIR.EQ.2)
     $                 .OR.((MYROW.EQ.RSRC1.OR.MYROW.EQ.RSRC4).AND.
     $                 DIR.EQ.1)) THEN
                     WINID = WINID + 1
                     LNWIN = LKBOT - LKTOP + 1
                     IPU = IPNEXT
                     DIM1 = NB - MOD(LKTOP-1,NB) 
                     DIM4 = LNWIN - DIM1
                     IF( ACCUM ) THEN
                        IPNEXT = IPU + LNWIN*LNWIN
                     ELSE
                        IPNEXT = IPU + 3*LNWIN
                     END IF
                     IF( DIR.EQ.2 ) THEN
                        IF( WANTZ ) THEN
                           ZROWS = NUMROC( N, NB, MYROW, DESCZ( RSRC_ ), 
     $                          NPROW )
                        ELSE
                           ZROWS = 0
                        END IF
                        IF( WANTT ) THEN
                           HROWS = NUMROC( LKTOP-1, NB, MYROW, 
     $                          DESCH( RSRC_ ), NPROW )
                        ELSE
                           HROWS = 0
                        END IF
                     ELSE
                        ZROWS = 0
                        HROWS = 0
                     END IF
                     IF( DIR.EQ.1 ) THEN
                        IF( WANTT ) THEN
                           HCOLS = NUMROC( N - (LKTOP+DIM1-1), NB, 
     $                          MYCOL, CSRC4, NPCOL )
                           IF( MYCOL.EQ.CSRC4 ) HCOLS = HCOLS - DIM4
                        ELSE
                           HCOLS = 0
                        END IF
                     ELSE
                        HCOLS = 0
                     END IF
                     IPW = MAX( 1 + LENRBUF + LENCBUF, IPW3 )   
                     IPW1 = IPW + HROWS * LNWIN
                     IF( WANTZ ) THEN
                        IPW2 = IPW1 + LNWIN * HCOLS
                        IPW3 = IPW2 + ZROWS * LNWIN
                     ELSE
                        IPW3 = IPW1 + LNWIN * HCOLS
                     END IF
                  END IF
*     
*     Let each process row and column involved in the updates
*     exchange data in H and Z with their neighbours 
*     
                  IF( DIR.EQ.2 .AND. WANTT .AND. LENCBUF.GT.0 ) THEN
                     IF( MYCOL.EQ.CSRC1 .OR. MYCOL.EQ.CSRC4 ) THEN
                        DO 210 INDX = 1, NPROW
                           IF( MYCOL.EQ.CSRC1 ) THEN
                              CALL INFOG2L( 1+(INDX-1)*NB, LKTOP, DESCH, 
     $                             NPROW, NPCOL, MYROW, MYCOL, ILOC, 
     $                             JLOC1, RSRC, CSRC1 )
                              IF( MYROW.EQ.RSRC ) THEN
                                 CALL DLACPY( 'All', HROWS, DIM1, 
     $                                H((JLOC1-1)*LLDH+ILOC), LLDH, 
     $                                WORK(IPW), HROWS )
                                 IF( NPCOL.GT.1 ) THEN
                                    EAST = MOD( MYCOL + 1, NPCOL )
                                    CALL DGESD2D( ICTXT, HROWS, DIM1,
     $                                   WORK(IPW), HROWS, RSRC, EAST )
                                    CALL DGERV2D( ICTXT, HROWS, DIM4,
     $                                   WORK(IPW+HROWS*DIM1), HROWS, 
     $                                   RSRC, EAST )
                                 END IF
                              END IF
                           END IF
                           IF( MYCOL.EQ.CSRC4 ) THEN
                              CALL INFOG2L( 1+(INDX-1)*NB, LKTOP+DIM1, 
     $                             DESCH, NPROW, NPCOL, MYROW, MYCOL, 
     $                             ILOC, JLOC4, RSRC, CSRC4 )
                              IF( MYROW.EQ.RSRC ) THEN
                                 CALL DLACPY( 'All', HROWS, DIM4, 
     $                                H((JLOC4-1)*LLDH+ILOC), LLDH, 
     $                                WORK(IPW+HROWS*DIM1), HROWS )
                                 IF( NPCOL.GT.1 ) THEN
                                    WEST = MOD( MYCOL - 1 + NPCOL, 
     $                                   NPCOL )
                                    CALL DGESD2D( ICTXT, HROWS, DIM4,
     $                                   WORK(IPW+HROWS*DIM1), HROWS, 
     $                                   RSRC, WEST )
                                    CALL DGERV2D( ICTXT, HROWS, DIM1,
     $                                   WORK(IPW), HROWS, RSRC, WEST )
                                 END IF
                              END IF
                           END IF
 210                    CONTINUE
                     END IF
                  END IF
*     
                  IF( DIR.EQ.1 .AND. WANTT .AND. LENRBUF.GT.0 ) THEN
                     IF( MYROW.EQ.RSRC1 .OR. MYROW.EQ.RSRC4 ) THEN
                        DO 220 INDX = 1, NPCOL
                           IF( MYROW.EQ.RSRC1 ) THEN
                              IF( INDX.EQ.1 ) THEN
                                 IF( DEBUG ) WRITE(*,*)MYROW,MYCOL, 
     $                                'LKTOP,LKBOT=',LKTOP,LKBOT
                                 IF( LKBOT.LT.N ) THEN
                                    CALL INFOG2L( LKTOP, LKBOT+1, DESCH, 
     $                                   NPROW, NPCOL, MYROW, MYCOL, 
     $                                   ILOC1, JLOC, RSRC1, CSRC )
                                 ELSE
                                    CSRC = -1
                                 END IF
                              ELSEIF( MOD(LKBOT,NB).NE.0 ) THEN
                                 IF( DEBUG ) WRITE(*,*)MYROW,MYCOL, 
     $                                'PDLAQR5: ICEIL 12'
                                 CALL INFOG2L( LKTOP,
     $                                (ICEIL(LKBOT,NB)+(INDX-2))*NB+1, 
     $                                DESCH, NPROW, NPCOL, MYROW, MYCOL, 
     $                                ILOC1, JLOC, RSRC1, CSRC )
                              ELSE
                                 CALL INFOG2L( LKTOP,
     $                                (ICEIL(LKBOT,NB)+(INDX-1))*NB+1, 
     $                                DESCH, NPROW, NPCOL, MYROW, MYCOL, 
     $                                ILOC1, JLOC, RSRC1, CSRC )
                              END IF
                              IF( DEBUG ) WRITE(*,*)MYROW,MYCOL, 
     $                             'ILOC1,JLOC=',
     $                             ILOC1,JLOC   
                              IF( MYCOL.EQ.CSRC ) THEN
                                 CALL DLACPY( 'All', DIM1, HCOLS, 
     $                                H((JLOC-1)*LLDH+ILOC1), LLDH, 
     $                                WORK(IPW1), LNWIN )
                                 IF( NPROW.GT.1 ) THEN
                                    SOUTH = MOD( MYROW + 1, NPROW )
                                    CALL DGESD2D( ICTXT, DIM1, HCOLS,
     $                                   WORK(IPW1), LNWIN, SOUTH, 
     $                                   CSRC )
                                    CALL DGERV2D( ICTXT, DIM4, HCOLS,
     $                                   WORK(IPW1+DIM1), LNWIN, SOUTH, 
     $                                   CSRC )
                                 END IF
                              END IF
                           END IF
                           IF( MYROW.EQ.RSRC4 ) THEN
                              IF( INDX.EQ.1 ) THEN
                                 IF( LKBOT.LT.N ) THEN
                                    CALL INFOG2L( LKTOP+DIM1, LKBOT+1, 
     $                                   DESCH, NPROW, NPCOL, MYROW, 
     $                                   MYCOL, ILOC4, JLOC, RSRC4, 
     $                                   CSRC )
                                 ELSE
                                    CSRC = -1
                                 END IF
                              ELSEIF( MOD(LKBOT,NB).NE.0 ) THEN
                                 CALL INFOG2L( LKTOP+DIM1, 
     $                                (ICEIL(LKBOT,NB)+(INDX-2))*NB+1, 
     $                                DESCH, NPROW, NPCOL, MYROW, MYCOL, 
     $                                ILOC4, JLOC, RSRC4, CSRC )
                              ELSE
                                 CALL INFOG2L( LKTOP+DIM1, 
     $                                (ICEIL(LKBOT,NB)+(INDX-1))*NB+1, 
     $                                DESCH, NPROW, NPCOL, MYROW, MYCOL, 
     $                                ILOC4, JLOC, RSRC4, CSRC )
                              END IF
                              IF( MYCOL.EQ.CSRC ) THEN
                                 CALL DLACPY( 'All', DIM4, HCOLS, 
     $                                H((JLOC-1)*LLDH+ILOC4), LLDH, 
     $                                WORK(IPW1+DIM1), LNWIN )
                                 IF( NPROW.GT.1 ) THEN
                                    NORTH = MOD( MYROW - 1 + NPROW, 
     $                                           NPROW )
                                    CALL DGESD2D( ICTXT, DIM4, HCOLS,
     $                                   WORK(IPW1+DIM1), LNWIN, NORTH, 
     $                                   CSRC )
                                    CALL DGERV2D( ICTXT, DIM1, HCOLS,
     $                                   WORK(IPW1), LNWIN, NORTH, 
     $                                   CSRC )
                                 END IF
                              END IF
                           END IF
 220                    CONTINUE
                     END IF
                  END IF
*     
                  IF( DIR.EQ.2 .AND. WANTZ .AND. LENCBUF.GT.0) THEN
                     IF( MYCOL.EQ.CSRC1 .OR. MYCOL.EQ.CSRC4 ) THEN
                        DO 230 INDX = 1, NPROW
                           IF( MYCOL.EQ.CSRC1 ) THEN
                              CALL INFOG2L( 1+(INDX-1)*NB, LKTOP, 
     $                             DESCZ, NPROW, NPCOL, MYROW, MYCOL, 
     $                             ILOC, JLOC1, RSRC, CSRC1 )
                              IF( MYROW.EQ.RSRC ) THEN
                                 CALL DLACPY( 'All', ZROWS, DIM1, 
     $                                Z((JLOC1-1)*LLDZ+ILOC), LLDZ, 
     $                                WORK(IPW2), ZROWS )
                                 IF( NPCOL.GT.1 ) THEN
                                    EAST = MOD( MYCOL + 1, NPCOL )
                                    CALL DGESD2D( ICTXT, ZROWS, DIM1,
     $                                   WORK(IPW2), ZROWS, RSRC, 
     $                                   EAST )
                                    CALL DGERV2D( ICTXT, ZROWS, DIM4,
     $                                   WORK(IPW2+ZROWS*DIM1), 
     $                                   ZROWS, RSRC, EAST )
                                 END IF
                              END IF
                           END IF
                           IF( MYCOL.EQ.CSRC4 ) THEN
                              CALL INFOG2L( 1+(INDX-1)*NB, 
     $                             LKTOP+DIM1, DESCZ, NPROW, NPCOL, 
     $                             MYROW, MYCOL, ILOC, JLOC4, RSRC, 
     $                             CSRC4 )
                              IF( MYROW.EQ.RSRC ) THEN
                                 CALL DLACPY( 'All', ZROWS, DIM4, 
     $                                Z((JLOC4-1)*LLDZ+ILOC), LLDZ, 
     $                                WORK(IPW2+ZROWS*DIM1), ZROWS )
                                 IF( NPCOL.GT.1 ) THEN
                                    WEST = MOD( MYCOL - 1 + NPCOL, 
     $                                   NPCOL )
                                    CALL DGESD2D( ICTXT, ZROWS, DIM4,
     $                                   WORK(IPW2+ZROWS*DIM1), 
     $                                   ZROWS, RSRC, WEST )
                                    CALL DGERV2D( ICTXT, ZROWS, DIM1,
     $                                   WORK(IPW2), ZROWS, RSRC, 
     $                                   WEST )
                                 END IF
                              END IF
                           END IF
 230                    CONTINUE
                     END IF
                  END IF
*
*     If no exchanges was performed for the current window, all 
*     processors jump to this point and try the next one
*
 205              CONTINUE
*
 200           CONTINUE
*     
*     Compute crossborder bulge-chase updates
*
               WINID = 0
               IF( DIR.EQ.1 ) THEN
                  IPNEXT = 1
               ELSE
                  IPNEXT = 1 + LENRBUF
               END IF
               IPW3 = 1
               DO 240 WIN = ODDEVEN+(CHUNKNUM-1)*WCHUNK, 
     $              MIN(ANMWIN,MAX(1,ODDEVEN+(CHUNKNUM)*WCHUNK-1)), 2
                  IF( DEBUG ) WRITE(*,*)MYROW,MYCOL, 
     $                 'Compute crossborder bulge-chase updates',
     $                 WIN,IWORK( 5+(WIN-1)*5 )
                  IF( IWORK( 5+(WIN-1)*5 ).NE.1 ) GO TO 245
*
*     Only perform this part of the code if there was really some
*     work performed on the WIN:th window
*
                  LKTOP = IWORK( 1+(WIN-1)*5 )
                  LKBOT = IWORK( 2+(WIN-1)*5 )
                  LNWIN = LKBOT - LKTOP + 1
*
*     Extract the processor indices associated with the current window
*
                  RSRC1 = IWORK( 3+(WIN-1)*5 )
                  CSRC1 = IWORK( 4+(WIN-1)*5 )
                  RSRC4 = MOD( RSRC1+1, NPROW )
                  CSRC4 = MOD( CSRC1+1, NPCOL )
*
                  IF(((MYCOL.EQ.CSRC1.OR.MYCOL.EQ.CSRC4).AND.DIR.EQ.2)
     $                 .OR.((MYROW.EQ.RSRC1.OR.MYROW.EQ.RSRC4).AND.
     $                 DIR.EQ.1)) THEN
*     
*     Set up workspaces
*     
                     WINID = WINID + 1
                     LKTOP = IWORK( 1+(WIN-1)*5 )
                     LKBOT = IWORK( 2+(WIN-1)*5 )
                     LNWIN = LKBOT - LKTOP + 1
                     DIM1 = NB - MOD(LKTOP-1,NB) 
                     DIM4 = LNWIN - DIM1
                     IF( ACCUM ) THEN
                        IPU = IPNEXT + (WINID-1)*LNWIN*LNWIN
                     ELSE
                        IPU = IPNEXT + (WINID-1)*3*LNWIN
                     END IF
                     IF( DEBUG ) WRITE(*,*)MYROW,MYCOL, 'WIN,IPU=',
     $                    WIN,IPU
                     IF( DIR.EQ.2 ) THEN
                        IF( WANTZ ) THEN
                           ZROWS = NUMROC( N, NB, MYROW, DESCZ( RSRC_ ), 
     $                          NPROW )
                        ELSE
                           ZROWS = 0
                        END IF
                        IF( WANTT ) THEN
                           HROWS = NUMROC( LKTOP-1, NB, MYROW, 
     $                          DESCH( RSRC_ ), NPROW )
                        ELSE
                           HROWS = 0
                        END IF
                     ELSE
                        ZROWS = 0
                        HROWS = 0
                     END IF
                     IF( DIR.EQ.1 ) THEN
                        IF( WANTT ) THEN
                           HCOLS = NUMROC( N - (LKTOP+DIM1-1), NB, 
     $                          MYCOL, CSRC4, NPCOL )
                           IF( MYCOL.EQ.CSRC4 ) HCOLS = HCOLS - DIM4
                        ELSE
                           HCOLS = 0
                        END IF
                     ELSE
                        HCOLS = 0
                     END IF
*     
*     IPW  = local copy of overlapping column block of H
*     IPW1 = local copy of overlapping row block of H
*     IPW2 = local copy of overlapping column block of Z
*     IPW3 = workspace for right hand side of matrix multiplication
*     
                     IPW = MAX( 1 + LENRBUF + LENCBUF, IPW3 )   
                     IPW1 = IPW + HROWS * LNWIN
                     IF( WANTZ ) THEN
                        IPW2 = IPW1 + LNWIN * HCOLS
                        IPW3 = IPW2 + ZROWS * LNWIN
                     ELSE
                        IPW3 = IPW1 + LNWIN * HCOLS
                     END IF
*     
*     Recompute job to see if special structure of U could possibly
*     be exploited
*     
                     IF( LKTOP.EQ.KTOP .AND. LKBOT.EQ.KBOT ) THEN
                        JOB = 'All steps'
                     ELSEIF( LKTOP.EQ.KTOP .AND. 
     $                       ( DIM1.LT.LCHAIN+1 .OR. DIM1.LE.NTINY ) ) 
     $                       THEN
                        JOB = 'Introduce and chase'
                     ELSEIF( LKBOT.EQ.KBOT ) THEN
                        JOB = 'Off-chase bulges'
                     ELSE
                        JOB = 'Chase bulges'
                     END IF
                  END IF
*     
*     Test if to exploit sparsity structure of orthogonal matrix U 
*     
                  KS = LNS
                  IF( ACCUM .AND. (.NOT. BLK22 .OR. DIM1.NE.KS .OR. 
     $                 LSAME(JOB,'I') .OR. LSAME(JOB,'O') .OR. 
     $                 LNS.LE.2 ) ) THEN
*     
*     Update the columns of H and Z 
*     
                     IF( DEBUG ) WRITE(*,*)MYROW,MYCOL,
     $                    'Cross border column updates...'
                     IF( DIR.EQ.2 .AND. WANTT .AND. LENCBUF.GT.0 ) THEN
                        DO 250 INDX = 1, MIN(LKTOP-1,1+(NPROW-1)*NB), NB
                           IF( MYCOL.EQ.CSRC1 ) THEN
                              CALL INFOG2L( INDX, LKTOP, DESCH, NPROW, 
     $                             NPCOL, MYROW, MYCOL, ILOC, JLOC, 
     $                             RSRC, CSRC1 )
                              IF( DEBUG ) WRITE(*,*)MYROW,MYCOL, 
     $                             'INDX,HROWS,DIM1=',
     $                             INDX, HROWS, DIM1, ILOC, JLOC
c                              CALL DLAPRNT(LNWIN,LNWIN,WORK(IPU),1,1,
c     $                             LNWIN,'U',6)
c                              CALL DLAPRNT(HROWS,DIM1+DIM4,WORK(IPW),
c     $                             1,1,HROWS,'W1',6)
                              IF( MYROW.EQ.RSRC ) THEN
#ifdef USE_OMP
*$OMP PARALLEL DEFAULT(SHARED), PRIVATE( INDX2, HROWS2 )
*$OMP DO
                                 DO 252 INDX2 = 1, HROWS, NB
                                    HROWS2 = MIN( NB, HROWS-INDX2+1 )
                                    CALL DGEMM( 'No transpose',
     $                                   'No transpose', HROWS2, DIM1, 
     $                                   LNWIN, ONE, WORK(IPW+INDX2-1), 
     $                                   HROWS, WORK( IPU ), LNWIN, 
     $                                   ZERO, WORK(IPW3+INDX2-1), 
     $                                   HROWS )
                                    CALL DLACPY( 'All', HROWS2, DIM1, 
     $                                   WORK(IPW3+INDX2-1), HROWS, 
     $                                   H((JLOC-1)*LLDH+ILOC+INDX2-1), 
     $                                   LLDH ) 
 252                             CONTINUE
*$OMP END DO 
*$OMP END PARALLEL
#else
                                 CALL DGEMM( 'No transpose',
     $                                'No transpose', HROWS, DIM1, 
     $                                LNWIN, ONE, WORK( IPW ), HROWS, 
     $                                WORK( IPU ), LNWIN, ZERO, 
     $                                WORK(IPW3), HROWS )
                                 CALL DLACPY( 'All', HROWS, DIM1, 
     $                                WORK(IPW3), HROWS, 
     $                                H((JLOC-1)*LLDH+ILOC), LLDH )
#endif
                              END IF
c                              CALL DLAPRNT(HROWS,DIM1,WORK(IPW3),1,
c     $                             1, HROWS,'W2',6)
                           END IF
                           IF( MYCOL.EQ.CSRC4 ) THEN
                              CALL INFOG2L( INDX, LKTOP+DIM1, DESCH, 
     $                             NPROW, NPCOL, MYROW, MYCOL, ILOC, 
     $                             JLOC, RSRC, CSRC4 )
                              IF( DEBUG ) WRITE(*,*)MYROW,MYCOL, 
     $                             'INDX,HROWS,DIM4=',
     $                             INDX, HROWS, DIM4, ILOC, JLOC
c                               CALL DLAPRNT(LNWIN,LNWIN,WORK(IPU),1,1,
c     $                             LNWIN,'U',6)
c                              CALL DLAPRNT(HROWS,DIM1+DIM4,WORK(IPW),
c     $                              1,1,HROWS,'W3',6)
                              IF( MYROW.EQ.RSRC ) THEN
#ifdef USE_OMP
*$OMP PARALLEL DEFAULT(SHARED), PRIVATE( INDX2, HROWS2 )
*$OMP DO
                                 DO 253 INDX2 = 1, HROWS, NB
                                    HROWS2 = MIN( NB, HROWS-INDX2+1 )
                                    CALL DGEMM( 'No transpose',
     $                                   'No transpose', HROWS2, DIM4, 
     $                                   LNWIN, ONE, WORK(IPW+INDX2-1), 
     $                                   HROWS, WORK( IPU+LNWIN*DIM1 ), 
     $                                   LNWIN, ZERO, 
     $                                   WORK(IPW3+INDX2-1), HROWS )
                                    CALL DLACPY( 'All', HROWS2, DIM4, 
     $                                   WORK(IPW3+INDX2-1), HROWS, 
     $                                   H((JLOC-1)*LLDH+ILOC+INDX2-1), 
     $                                   LLDH ) 
 253                             CONTINUE
*$OMP END DO 
*$OMP END PARALLEL
#else
                                 CALL DGEMM( 'No transpose',
     $                                'No transpose', HROWS, DIM4, 
     $                                LNWIN, ONE, WORK( IPW ), HROWS, 
     $                                WORK( IPU+LNWIN*DIM1 ), LNWIN, 
     $                                ZERO, WORK(IPW3), HROWS )
                                 CALL DLACPY( 'All', HROWS, DIM4, 
     $                                WORK(IPW3), HROWS, 
     $                                H((JLOC-1)*LLDH+ILOC), LLDH )
#endif
                              END IF
c                              CALL DLAPRNT(HROWS,DIM4,WORK(IPW3),
c     $                             1,1,HROWS,'W4',6)
                           END IF
 250                    CONTINUE
                     END IF
*
                     IF( DIR.EQ.2 .AND. WANTZ .AND. LENCBUF.GT.0 ) THEN
                        IF( DEBUG ) WRITE(*,*) MYROW,MYCOL,
     $                       'Updating Z...'
                        DO 260 INDX = 1, MIN(N,1+(NPROW-1)*NB), NB
                           IF( DEBUG ) WRITE(*,*) MYROW,MYCOL,
     $                          'INDX=',INDX
                           IF( MYCOL.EQ.CSRC1 ) THEN
                              CALL INFOG2L( INDX, LKTOP, DESCZ, NPROW, 
     $                             NPCOL, MYROW, MYCOL, ILOC, JLOC, 
     $                             RSRC, CSRC1 )
                              IF( MYROW.EQ.RSRC ) THEN
                                 IF( DEBUG ) WRITE(*,*) MYROW,MYCOL,
     $                                'Z*Uloc - CSRC1, WIN,IPU=',WIN,IPU
                                 IF( PRINT )
     $                                CALL DLAPRNT(LNWIN,LNWIN,
     $                                WORK(IPU),1,1,LNWIN,'Uloc1',6)
                                 IF( DEBUG ) THEN
                                    CALL DLACPY( 'All', LNWIN, LNWIN,
     $                                   WORK(IPU), LNWIN, WORK(IPW3),
     $                                   LNWIN )
                                    CALL DLASET( 'ALL', LNWIN, LNWIN, 
     $                                   ZERO, ONE, 
     $                                   WORK(IPW3+LNWIN*LNWIN), 
     $                                   LNWIN )
                                    CALL DGEMM('T','N',LNWIN,LNWIN,
     $                                   LNWIN,ONE,WORK(IPU),LNWIN,
     $                                   WORK(IPW3),LNWIN,-ONE,
     $                                   WORK(IPW3+LNWIN*LNWIN),
     $                                   LNWIN )
                                    DDUM = DLANGE('F',LNWIN,LNWIN,
     $                                   WORK(IPW3+LNWIN*LNWIN),LNWIN,
     $                                   WORK(IPW3+2*LNWIN*LNWIN))
                                 WRITE(*,*) MYROW,MYCOL, 
     $                                'NORM(Uloc1^T*Uloc1-eye)=',DDUM
                                 END IF
                                 IF( DEBUG ) WRITE(*,*) MYROW,MYCOL,
     $                                'ZROWS,DIM1,INDX=',ZROWS,DIM1,INDX
#ifdef USE_OMP
*$OMP PARALLEL DEFAULT(SHARED), PRIVATE( INDX2, ZROWS2 )
*$OMP DO
                                 DO 262 INDX2 = 1, ZROWS, NB
                                    ZROWS2 = MIN( NB, ZROWS-INDX2+1 )
                                    CALL DGEMM( 'No transpose', 
     $                                   'No transpose', ZROWS2, DIM1, 
     $                                   LNWIN, ONE, 
     $                                   WORK( IPW2+INDX2-1 ), 
     $                                   ZROWS, WORK( IPU ), LNWIN, 
     $                                   ZERO, WORK(IPW3+INDX2-1), 
     $                                   ZROWS )
                                    CALL DLACPY( 'All', ZROWS2, DIM1, 
     $                                   WORK(IPW3+INDX2-1), ZROWS, 
     $                                   Z((JLOC-1)*LLDZ+ILOC+INDX2-1), 
     $                                   LLDZ )
 262                             CONTINUE
*$OMP END DO 
*$OMP END PARALLEL
#else
                                 CALL DGEMM( 'No transpose', 
     $                                'No transpose', ZROWS, DIM1, 
     $                                LNWIN, ONE, WORK( IPW2 ), 
     $                                ZROWS, WORK( IPU ), LNWIN, 
     $                                ZERO, WORK(IPW3), ZROWS )
                                 CALL DLACPY( 'All', ZROWS, DIM1, 
     $                                WORK(IPW3), ZROWS, 
     $                                Z((JLOC-1)*LLDZ+ILOC), LLDZ )
#endif
                              END IF
                           END IF
                           IF( MYCOL.EQ.CSRC4 ) THEN
                              CALL INFOG2L( INDX, LKTOP+DIM1, DESCZ, 
     $                             NPROW, NPCOL, MYROW, MYCOL, ILOC, 
     $                             JLOC, RSRC, CSRC4 )
                              IF( MYROW.EQ.RSRC ) THEN
                                 IF( DEBUG ) WRITE(*,*) MYROW,MYCOL,
     $                                'Z*Uloc - CSRC4, WIN,IPU=',WIN,IPU
                                 IF( PRINT ) 
     $                                CALL DLAPRNT(LNWIN,LNWIN,
     $                                WORK(IPU),1,1,LNWIN,'Uloc4',6)
                                 IF( DEBUG ) THEN
                                    CALL DLACPY( 'All', LNWIN, LNWIN,
     $                                   WORK(IPU), LNWIN, WORK(IPW3),
     $                                   LNWIN )
                                    CALL DLASET( 'ALL', LNWIN, LNWIN, 
     $                                   ZERO, ONE, 
     $                                   WORK(IPW3+LNWIN*LNWIN), LNWIN )
                                    CALL DGEMM('T','N',LNWIN,LNWIN,
     $                                   LNWIN,ONE,WORK(IPU),LNWIN,
     $                                   WORK(IPW3),LNWIN,-ONE,
     $                                   WORK(IPW3+LNWIN*LNWIN),LNWIN )
                                    DDUM = DLANGE('F',LNWIN,LNWIN,
     $                                   WORK(IPW3+LNWIN*LNWIN),LNWIN,
     $                                   WORK(IPW3+2*LNWIN*LNWIN))
                                    WRITE(*,*) MYROW,MYCOL,
     $                                   'NORM(Uloc4^T*Uloc4-eye)=',DDUM
                                 END IF
                                 IF( DEBUG ) WRITE(*,*) MYROW,MYCOL,
     $                                'ZROWS,DIM4,INDX=',ZROWS,DIM4,INDX
#ifdef USE_OMP
*$OMP PARALLEL DEFAULT(SHARED), PRIVATE( INDX2, ZROWS2 )
*$OMP DO
                                 DO 263 INDX2 = 1, ZROWS, NB
                                    ZROWS2 = MIN( NB, ZROWS-INDX2+1 )
                                    CALL DGEMM( 'No transpose', 
     $                                   'No transpose', ZROWS2, DIM4, 
     $                                   LNWIN, ONE, 
     $                                   WORK( IPW2+INDX2-1 ), 
     $                                   ZROWS, WORK( IPU+LNWIN*DIM1 ), 
     $                                   LNWIN, ZERO, 
     $                                   WORK(IPW3+INDX2-1), ZROWS )
                                    CALL DLACPY( 'All', ZROWS2, DIM4, 
     $                                   WORK(IPW3+INDX2-1), ZROWS, 
     $                                   Z((JLOC-1)*LLDZ+ILOC+INDX2-1), 
     $                                   LLDZ )
 263                             CONTINUE
*$OMP END DO 
*$OMP END PARALLEL
#else
                                 CALL DGEMM( 'No transpose', 
     $                                'No transpose', ZROWS, DIM4, 
     $                                LNWIN, ONE, WORK( IPW2 ), 
     $                                ZROWS, 
     $                                WORK( IPU+LNWIN*DIM1 ), LNWIN, 
     $                                ZERO, WORK(IPW3), ZROWS )
                                 CALL DLACPY( 'All', ZROWS, DIM4, 
     $                                WORK(IPW3), ZROWS, 
     $                                Z((JLOC-1)*LLDZ+ILOC), LLDZ )
#endif
                              END IF
                           END IF
 260                    CONTINUE
                     END IF
*     
*     Update the rows of H 
*     
                      IF( DEBUG ) WRITE(*,*)MYROW,MYCOL,
     $                    'Cross border row updates...'
                     IF( DIR.EQ.1 .AND. WANTT .AND. LENRBUF.GT.0 ) THEN
                        IF( LKBOT.LT.N ) THEN
                           IF( MYROW.EQ.RSRC1.AND.MYCOL.EQ.CSRC4 .AND.
     $                          MOD(LKBOT,NB).NE.0 ) THEN
                              INDX = LKBOT + 1
                              CALL INFOG2L( LKTOP, INDX, DESCH, NPROW, 
     $                             NPCOL, MYROW, MYCOL, ILOC, JLOC, 
     $                             RSRC1, CSRC4 )
                              IF( DEBUG ) WRITE(*,*)MYROW,MYCOL, 
     $                             'INDX,HCOLS 1=',
     $                             INDX,HCOLS, LNWIN, DIM1
                              IF( DEBUG ) WRITE(*,*) MYROW,MYCOL,
     $                                   'U14: WIN,IPU',WIN,IPU
c                              CALL DLAPRNT(LNWIN,LNWIN,WORK(IPU),1,1,
c     $                             LNWIN,'U14',6)
c                              CALL DLAPRNT(DIM1,HCOLS,WORK(IPW1),1,1,
c     $                             LNWIN,'W1',6)
                              CALL DGEMM( 'Transpose', 'No Transpose', 
     $                             DIM1, HCOLS, LNWIN, ONE, WORK(IPU), 
     $                             LNWIN, WORK( IPW1 ), LNWIN, ZERO, 
     $                             WORK(IPW3), DIM1 )
c                              CALL DLAPRNT(DIM1,HCOLS,WORK(IPW3),1,1,
c     $                             LNWIN,'W2',6)
                              CALL DLACPY( 'All', DIM1, HCOLS, 
     $                             WORK(IPW3), DIM1, 
     $                             H((JLOC-1)*LLDH+ILOC), LLDH )
                           END IF
                           IF( MYROW.EQ.RSRC4.AND.MYCOL.EQ.CSRC4 .AND.
     $                          MOD(LKBOT,NB).NE.0 ) THEN
                              INDX = LKBOT + 1
                              CALL INFOG2L( LKTOP+DIM1, INDX, DESCH, 
     $                             NPROW, NPCOL, MYROW, MYCOL, ILOC, 
     $                             JLOC, RSRC4, CSRC4 )
                              IF( DEBUG ) WRITE(*,*)MYROW,MYCOL, 
     $                             'INDX,HCOLS 2=',
     $                             INDX,HCOLS, LNWIN, DIM4
                              IF( DEBUG ) WRITE(*,*) MYROW,MYCOL,
     $                                   'U44: WIN,IPU',WIN,IPU
c                              CALL DLAPRNT(LNWIN,LNWIN,WORK(IPU),1,1,
c     $                             LNWIN,'U44',6)
c                              CALL DLAPRNT(DIM4,HCOLS,WORK(IPW1),1,1,
c     $                             LNWIN,'W3',6)
                              CALL DGEMM( 'Transpose', 'No Transpose', 
     $                             DIM4, HCOLS, LNWIN, ONE, 
     $                             WORK( IPU+DIM1*LNWIN ), LNWIN, 
     $                             WORK( IPW1), LNWIN, ZERO, 
     $                             WORK(IPW3), DIM4 )
c                              CALL DLAPRNT(DIM4,HCOLS,WORK(IPW3),1,1,
c     $                             LNWIN,'W4',6)
                              CALL DLACPY( 'All', DIM4, HCOLS, 
     $                             WORK(IPW3), DIM4, 
     $                             H((JLOC-1)*LLDH+ILOC), LLDH )
                           END IF
                           IF( DEBUG ) WRITE(*,*)MYROW,MYCOL, 
     $                          'PDLAQR5: ICEIL 13'
                           INDXS = ICEIL(LKBOT,NB)*NB + 1 
                           IF( MOD(LKBOT,NB).NE.0 ) THEN
                              INDXE = MIN(N,INDXS+(NPCOL-2)*NB)
                           ELSE
                              INDXE = MIN(N,INDXS+(NPCOL-1)*NB)
                           END IF
                           IF( DEBUG ) WRITE(*,*)MYROW,MYCOL, 
     $                          'INDXS,INDXE=',INDXS,INDXE
                           DO 270 INDX = INDXS, INDXE, NB
                              IF( MYROW.EQ.RSRC1 ) THEN
                                 CALL INFOG2L( LKTOP, INDX, DESCH, 
     $                                NPROW, NPCOL, MYROW, MYCOL, ILOC, 
     $                                JLOC, RSRC1, CSRC )
                                 IF( MYCOL.EQ.CSRC ) THEN
                                    IF( DEBUG ) WRITE(*,*)MYROW,MYCOL,
     $                                   'INDX,HCOLS=',
     $                                   INDX,HCOLS
                                    IF( DEBUG ) WRITE(*,*) MYROW,MYCOL,
     $                                   'U1: WIN,IPU',WIN,IPU
c                                    CALL DLAPRNT(LNWIN,LNWIN,WORK(IPU),
c     $                                   1,1,LNWIN,'U1',6)
#ifdef USE_OMP
*$OMP PARALLEL DEFAULT(SHARED), PRIVATE( INDX2, HCOLS2 )
*$OMP DO
                                    DO 272 INDX2 = 1, HCOLS, NB
                                       HCOLS2 = MIN( NB, HCOLS-INDX2+1 )
                                       CALL DGEMM( 'Transpose',
     $                                      'No Transpose', DIM1, 
     $                                      HCOLS2, LNWIN, ONE, 
     $                                      WORK( IPU ), LNWIN, 
     $                                      WORK(IPW1+LNWIN*(INDX2-1)), 
     $                                      LNWIN, ZERO, 
     $                                      WORK(IPW3+LNWIN*(INDX2-1)), 
     $                                      DIM1 )
                                       CALL DLACPY( 'All', DIM1, HCOLS2, 
     $                                      WORK(IPW3+LNWIN*(INDX2-1)), 
     $                                      DIM1, 
     $                                      H((JLOC-1+INDX2-1)*LLDH+
     $                                      ILOC), LLDH )
 272                                CONTINUE
*$OMP END DO 
*$OMP END PARALLEL
#else
                                    CALL DGEMM( 'Transpose',
     $                                   'No Transpose', DIM1, HCOLS, 
     $                                   LNWIN, ONE, WORK( IPU ), LNWIN, 
     $                                   WORK( IPW1 ), LNWIN, ZERO, 
     $                                   WORK(IPW3), DIM1 )
                                    CALL DLACPY( 'All', DIM1, HCOLS, 
     $                                   WORK(IPW3), DIM1, 
     $                                   H((JLOC-1)*LLDH+ILOC), LLDH )
#endif
                                 END IF
                              END IF
                              IF( MYROW.EQ.RSRC4 ) THEN
                                 CALL INFOG2L( LKTOP+DIM1, INDX, DESCH, 
     $                                NPROW, NPCOL, MYROW, MYCOL, ILOC, 
     $                                JLOC, RSRC4, CSRC )
                                 IF( MYCOL.EQ.CSRC ) THEN
                                    IF( DEBUG ) WRITE(*,*)MYROW,MYCOL,
     $                                   'INDX,HCOLS=',
     $                                   INDX,HCOLS
                                    IF( DEBUG ) WRITE(*,*) MYROW,MYCOL,
     $                                   'U4: WIN,IPU',WIN,IPU
c                                    CALL DLAPRNT(LNWIN,LNWIN,WORK(IPU),
c     $                                   1,1,LNWIN,'U4',6)
#ifdef USE_OMP
*$OMP PARALLEL DEFAULT(SHARED), PRIVATE( INDX2, HCOLS2 )
*$OMP DO
                                    DO 273 INDX2 = 1, HCOLS, NB
                                       HCOLS2 = MIN( NB, HCOLS-INDX2+1 )
                                       CALL DGEMM( 'Transpose',
     $                                      'No Transpose', DIM4, 
     $                                      HCOLS2, LNWIN, ONE, 
     $                                      WORK( IPU+LNWIN*DIM1 ), 
     $                                      LNWIN, 
     $                                      WORK(IPW1+LNWIN*(INDX2-1)), 
     $                                      LNWIN, ZERO, 
     $                                      WORK(IPW3+LNWIN*(INDX2-1)), 
     $                                      DIM4 )
                                       CALL DLACPY( 'All', DIM4, HCOLS2, 
     $                                      WORK(IPW3+LNWIN*(INDX2-1)), 
     $                                      DIM4, 
     $                                      H((JLOC-1+INDX2-1)*LLDH+
     $                                      ILOC), LLDH )
 273                                CONTINUE
*$OMP END DO 
*$OMP END PARALLEL
#else
                                    CALL DGEMM( 'Transpose',
     $                                   'No Transpose', DIM4, HCOLS, 
     $                                   LNWIN, ONE, 
     $                                   WORK( IPU+LNWIN*DIM1 ), LNWIN, 
     $                                   WORK( IPW1 ), LNWIN, 
     $                                   ZERO, WORK(IPW3), DIM4 )
                                    CALL DLACPY( 'All', DIM4, HCOLS, 
     $                                   WORK(IPW3), DIM4, 
     $                                   H((JLOC-1)*LLDH+ILOC), LLDH )
#endif
                                 END IF
                              END IF
 270                       CONTINUE
                        END IF
                     END IF
                  ELSEIF( ACCUM .AND. BLK22 ) THEN 
*     
*     Update the columns of H and Z
*     
*     Compute H2*U21 + H1*U11 on the left side of the border.
*     
                     IF( DIR.EQ.2 .AND. WANTT .AND. LENCBUF.GT.0 ) THEN
                        INDXE = MIN(LKTOP-1,1+(NPROW-1)*NB)
                        DO 280 INDX = 1, INDXE, NB
                           IF( MYCOL.EQ.CSRC1 ) THEN
                              CALL INFOG2L( INDX, LKTOP, DESCH, NPROW, 
     $                             NPCOL, MYROW, MYCOL, ILOC, JLOC, 
     $                             RSRC, CSRC1 )
                              IF( MYROW.EQ.RSRC ) THEN
                                 CALL DLACPY( 'All', HROWS, KS, 
     $                                WORK( IPW+HROWS*DIM4), HROWS,
     $                                WORK(IPW3), HROWS )
                                 CALL DTRMM( 'Right', 'Upper', 
     $                                'No transpose',
     $                                'No transpose', HROWS, KS, ONE,
     $                                WORK( IPU+DIM4 ), LNWIN, 
     $                                WORK(IPW3), HROWS )
                                 CALL DGEMM( 'No transpose', 
     $                                'No transpose', HROWS, KS, DIM4, 
     $                                ONE, WORK( IPW ), HROWS, 
     $                                WORK( IPU ), LNWIN, ONE, 
     $                                WORK(IPW3), HROWS )
                                 CALL DLACPY( 'All', HROWS, KS, 
     $                                WORK(IPW3), HROWS, 
     $                                H((JLOC-1)*LLDH+ILOC), LLDH )
                              END IF
                           END IF
*     
*     Compute H1*U12 + H2*U22 on the right side of the border.
*     
                           IF( MYCOL.EQ.CSRC4 ) THEN
                              CALL INFOG2L( INDX, LKTOP+DIM1, DESCH, 
     $                             NPROW, NPCOL, MYROW, MYCOL, ILOC, 
     $                             JLOC, RSRC, CSRC4 )
                              IF( MYROW.EQ.RSRC ) THEN
                                 CALL DLACPY( 'All', HROWS, DIM4, 
     $                                WORK(IPW), HROWS, WORK( IPW3 ), 
     $                                HROWS )
                                 CALL DTRMM( 'Right', 'Lower', 
     $                                'No transpose',
     $                                'No transpose', HROWS, DIM4, ONE,
     $                                WORK( IPU+LNWIN*KS ), LNWIN,
     $                                WORK( IPW3 ), HROWS )
                                 CALL DGEMM( 'No transpose', 
     $                                'No transpose', HROWS, DIM4, KS, 
     $                                ONE, WORK( IPW+HROWS*DIM4), 
     $                                HROWS,
     $                                WORK( IPU+LNWIN*KS+DIM4 ), LNWIN,
     $                                ONE, WORK( IPW3 ), HROWS )
                                 CALL DLACPY( 'All', HROWS, DIM4, 
     $                                WORK(IPW3), HROWS, 
     $                                H((JLOC-1)*LLDH+ILOC), LLDH )
                              END IF 
                           END IF
 280                    CONTINUE
                     END IF
*
                     IF( DIR.EQ.2 .AND. WANTZ .AND. LENCBUF.GT.0 ) THEN
*     
*     Compute Z2*U21 + Z1*U11 on the left side of border.
*     
                        INDXE = MIN(N,1+(NPROW-1)*NB)
                        DO 290 INDX = 1, INDXE, NB
                           IF( MYCOL.EQ.CSRC1 ) THEN
                              CALL INFOG2L( INDX, I, DESCZ, NPROW, 
     $                             NPCOL, MYROW, MYCOL, ILOC, JLOC, 
     $                             RSRC, CSRC1 )
                              IF( MYROW.EQ.RSRC ) THEN
                                 CALL DLACPY( 'All', ZROWS, KS, 
     $                                WORK( IPW2+ZROWS*DIM4), 
     $                                ZROWS, WORK(IPW3), ZROWS )
                                 CALL DTRMM( 'Right', 'Upper', 
     $                                'No transpose',
     $                                'No transpose', ZROWS, KS, ONE,
     $                                WORK( IPU+DIM4 ), LNWIN, 
     $                                WORK(IPW3), ZROWS )
                                 CALL DGEMM( 'No transpose', 
     $                                'No transpose', ZROWS, KS, 
     $                                DIM4, ONE, WORK( IPW2 ), 
     $                                ZROWS, WORK( IPU ), LNWIN, 
     $                                ONE, WORK(IPW3), ZROWS )
                                 CALL DLACPY( 'All', ZROWS, KS, 
     $                                WORK(IPW3), ZROWS, 
     $                                Z((JLOC-1)*LLDZ+ILOC), LLDZ )
                              END IF
                           END IF
*     
*     Compute Z1*U12 + Z2*U22 on the right side of border.
*     
                           IF( MYCOL.EQ.CSRC4 ) THEN
                              CALL INFOG2L( INDX, I+DIM1, DESCZ, 
     $                             NPROW, NPCOL, MYROW, MYCOL, ILOC, 
     $                             JLOC, RSRC, CSRC4 )
                              IF( MYROW.EQ.RSRC ) THEN
                                 CALL DLACPY( 'All', ZROWS, DIM4, 
     $                                WORK(IPW2), ZROWS, 
     $                                WORK( IPW3 ), ZROWS )
                                 CALL DTRMM( 'Right', 'Lower', 
     $                                'No transpose',
     $                                'No transpose', ZROWS, DIM4, 
     $                                ONE, WORK( IPU+LNWIN*KS ), 
     $                                LNWIN, WORK( IPW3 ), ZROWS )
                                 CALL DGEMM( 'No transpose', 
     $                                'No transpose', ZROWS, DIM4, 
     $                                KS, ONE, 
     $                                WORK( IPW2+ZROWS*(DIM4)), 
     $                                ZROWS,
     $                                WORK( IPU+LNWIN*KS+DIM4 ), 
     $                                LNWIN, ONE, WORK( IPW3 ), 
     $                                ZROWS )
                                 CALL DLACPY( 'All', ZROWS, DIM4, 
     $                                WORK(IPW3), ZROWS, 
     $                                Z((JLOC-1)*LLDZ+ILOC), LLDZ )
                              END IF 
                           END IF
 290                    CONTINUE
                     END IF
*     
                     IF( DIR.EQ.1 .AND. WANTT .AND. LENRBUF.GT.0) THEN
                        IF ( LKBOT.LT.N ) THEN
*     
*     Compute U21**T*H2 + U11**T*H1 on the upper side of the border.
*     
                           IF( MYROW.EQ.RSRC1.AND.MYCOL.EQ.CSRC4.AND.
     $                          MOD(LKBOT,NB).NE.0 ) THEN
                              INDX = LKBOT + 1
                              CALL INFOG2L( LKTOP, INDX, DESCH, NPROW, 
     $                             NPCOL, MYROW, MYCOL, ILOC, JLOC, 
     $                             RSRC1, CSRC4 ) 
                              CALL DLACPY( 'All', KS, HCOLS,
     $                             WORK( IPW1+DIM4 ), LNWIN, 
     $                             WORK(IPW3), KS )
                              CALL DTRMM( 'Left', 'Upper', 'Transpose',
     $                             'No transpose', KS, HCOLS, ONE,
     $                             WORK( IPU+DIM4 ), LNWIN, 
     $                             WORK(IPW3), KS )
                              CALL DGEMM( 'Transpose', 'No transpose', 
     $                             KS, HCOLS, DIM4, ONE, WORK(IPU), 
     $                             LNWIN, WORK(IPW1), LNWIN, 
     $                             ONE, WORK(IPW3), KS )
                              CALL DLACPY( 'All', KS, HCOLS, 
     $                             WORK(IPW3), KS, 
     $                             H((JLOC-1)*LLDH+ILOC), LLDH )
                           END IF
*     
*     Compute U12**T*H1 + U22**T*H2 one the lower side of the border.
*     
                           IF( MYROW.EQ.RSRC4.AND.MYCOL.EQ.CSRC4.AND.
     $                          MOD(LKBOT,NB).NE.0 ) THEN
                              INDX = LKBOT + 1
                              CALL INFOG2L( LKTOP+DIM1, INDX, DESCH, 
     $                             NPROW, NPCOL, MYROW, MYCOL, ILOC, 
     $                             JLOC, RSRC4, CSRC4 )
                              CALL DLACPY( 'All', DIM4, HCOLS,
     $                             WORK( IPW1 ), LNWIN, 
     $                             WORK( IPW3 ), DIM4 )
                              CALL DTRMM( 'Left', 'Lower', 'Transpose',
     $                             'No transpose', DIM4, HCOLS, ONE,
     $                             WORK( IPU+LNWIN*KS ), LNWIN,
     $                             WORK( IPW3 ), DIM4 )
                              CALL DGEMM( 'Transpose', 'No Transpose', 
     $                             DIM4, HCOLS, KS, ONE,
     $                             WORK( IPU+LNWIN*KS+DIM4 ), LNWIN,
     $                             WORK( IPW1+DIM1 ), LNWIN,
     $                             ONE, WORK( IPW3), DIM4 )
                              CALL DLACPY( 'All', DIM4, HCOLS, 
     $                             WORK(IPW3), DIM4, 
     $                             H((JLOC-1)*LLDH+ILOC), LLDH )
                           END IF
*     
*     Compute U21**T*H2 + U11**T*H1 on upper side on border.
*     
                           IF( DEBUG ) WRITE(*,*)MYROW,MYCOL, 
     $                          'PDLAQR5: ICEIL 14'
                           INDXS = ICEIL(LKBOT,NB)*NB+1
                           IF( MOD(LKBOT,NB).NE.0 ) THEN
                              INDXE = MIN(N,INDXS+(NPCOL-2)*NB)
                           ELSE
                              INDXE = MIN(N,INDXS+(NPCOL-1)*NB)
                           END IF
                           DO 300 INDX = INDXS, INDXE, NB
                              IF( MYROW.EQ.RSRC1 ) THEN
                                 CALL INFOG2L( LKTOP, INDX, DESCH, 
     $                                NPROW, NPCOL, MYROW, MYCOL, ILOC, 
     $                                JLOC, RSRC1, CSRC ) 
                                 IF( MYCOL.EQ.CSRC ) THEN
                                    CALL DLACPY( 'All', KS, HCOLS,
     $                                   WORK( IPW1+DIM4 ), LNWIN, 
     $                                   WORK(IPW3), KS )
                                    CALL DTRMM( 'Left', 'Upper', 
     $                                   'Transpose', 'No transpose', 
     $                                   KS, HCOLS, ONE, 
     $                                   WORK( IPU+DIM4 ), LNWIN, 
     $                                   WORK(IPW3), KS )
                                    CALL DGEMM( 'Transpose',
     $                                   'No transpose', KS, HCOLS, 
     $                                   DIM4, ONE, WORK(IPU), LNWIN, 
     $                                   WORK(IPW1), LNWIN, ONE, 
     $                                   WORK(IPW3), KS )
                                    CALL DLACPY( 'All', KS, HCOLS, 
     $                                   WORK(IPW3), KS, 
     $                                   H((JLOC-1)*LLDH+ILOC), LLDH )
                                 END IF
                              END IF
*     
*     Compute U12**T*H1 + U22**T*H2 on lower side of border.
*     
                              IF( MYROW.EQ.RSRC4 ) THEN
                                 CALL INFOG2L( LKTOP+DIM1, INDX, DESCH, 
     $                                NPROW, NPCOL, MYROW, MYCOL, ILOC, 
     $                                JLOC, RSRC4, CSRC )
                                 IF( MYCOL.EQ.CSRC ) THEN
                                    CALL DLACPY( 'All', DIM4, HCOLS,
     $                                   WORK( IPW1 ), LNWIN, 
     $                                   WORK( IPW3 ), DIM4 )
                                    CALL DTRMM( 'Left', 'Lower', 
     $                                   'Transpose','No transpose', 
     $                                   DIM4, HCOLS, ONE,
     $                                   WORK( IPU+LNWIN*KS ), LNWIN,
     $                                   WORK( IPW3 ), DIM4 )
                                    CALL DGEMM( 'Transpose',
     $                                   'No Transpose', DIM4, HCOLS, 
     $                                   KS, ONE,
     $                                   WORK( IPU+LNWIN*KS+DIM4 ), 
     $                                   LNWIN, WORK( IPW1+DIM1 ), 
     $                                   LNWIN, ONE, WORK( IPW3), 
     $                                   DIM4 )
                                    CALL DLACPY( 'All', DIM4, HCOLS, 
     $                                   WORK(IPW3), DIM4, 
     $                                   H((JLOC-1)*LLDH+ILOC), LLDH )
                                 END IF
                              END IF
 300                       CONTINUE
                        END IF
                     END IF
                  ELSE 
*
*     Update cross border regions using nonaccumulated transformations
*
                     IF( DIR.EQ.2 .AND. WANTT .AND. LENCBUF.GT.0 ) THEN
                        DO 301 INDX = 1, MIN(LKTOP-1,1+(NPROW-1)*NB), NB
                           IF( MYCOL.EQ.CSRC1 ) THEN
                              CALL INFOG2L( INDX, LKTOP, DESCH, NPROW, 
     $                             NPCOL, MYROW, MYCOL, ILOC, JLOC, 
     $                             RSRC, CSRC1 )
                              IF( MYROW.EQ.RSRC ) THEN
                                 DO 302 IHH = 1, LNWIN 
                                    TAU = WORK(IPU+(IHH-1)*3)
                                    WORK(IPU+(IHH-1)*3) = ONE
                                    CALL DLARFX ( 'Right', HROWS, 
     $                                   MIN(3,LNWIN-IHH+1),
     $                                   WORK(IPU+(IHH-1)*3), TAU, 
     $                                   WORK(IPW+(IHH-1)*HROWS), 
     $                                   HROWS, WORK(IPW3) )
                                    WORK(IPU+(IHH-1)*3) = TAU
 302                             CONTINUE
                                 CALL DLACPY( 'All', HROWS, DIM1, 
     $                                WORK(IPW3), HROWS, 
     $                                H((JLOC-1)*LLDH+ILOC), LLDH )
                              END IF
                           END IF
                           IF( MYCOL.EQ.CSRC4 ) THEN
                              CALL INFOG2L( INDX, LKTOP+DIM1, DESCH, 
     $                             NPROW, NPCOL, MYROW, MYCOL, ILOC, 
     $                             JLOC, RSRC, CSRC4 )
                              IF( MYROW.EQ.RSRC ) THEN
                                 DO 303 IHH = 1, LNWIN 
                                    TAU = WORK(IPU+(IHH-1)*3)
                                    WORK(IPU+(IHH-1)*3) = ONE
                                    CALL DLARFX ( 'Right', HROWS, 
     $                                   MIN(3,LNWIN-IHH+1),
     $                                   WORK(IPU+(IHH-1)*3), TAU, 
     $                                   WORK(IPW+(IHH-1)*HROWS), 
     $                                   HROWS, WORK(IPW3) )
                                    WORK(IPU+(IHH-1)*3) = TAU
 303                             CONTINUE
                                 CALL DLACPY( 'All', HROWS, DIM4, 
     $                                WORK(IPW3+DIM1*HROWS), HROWS, 
     $                                H((JLOC-1)*LLDH+ILOC), LLDH )
                              END IF
                           END IF
 301                    CONTINUE
                     END IF
*
                     IF( DIR.EQ.2 .AND. WANTZ .AND. LENCBUF.GT.0 ) THEN
                        DO 304 INDX = 1, MIN(N,1+(NPROW-1)*NB), NB
                           IF( MYCOL.EQ.CSRC1 ) THEN
                              CALL INFOG2L( INDX, LKTOP, DESCZ, NPROW, 
     $                             NPCOL, MYROW, MYCOL, ILOC, JLOC, 
     $                             RSRC, CSRC1 )
                              IF( MYROW.EQ.RSRC ) THEN
                                 DO 305 IHH = 1, LNWIN 
                                    TAU = WORK(IPU+(IHH-1)*3)
                                    WORK(IPU+(IHH-1)*3) = ONE
                                    CALL DLARFX ( 'Right', ZROWS, 
     $                                   MIN(3,LNWIN-IHH+1),
     $                                   WORK(IPU+(IHH-1)*3), TAU, 
     $                                   WORK(IPW2+(IHH-1)*ZROWS), 
     $                                   ZROWS, WORK(IPW3) )
                                    WORK(IPU+(IHH-1)*3) = TAU
 305                             CONTINUE
                                 CALL DLACPY( 'All', ZROWS, DIM1, 
     $                                WORK(IPW3), ZROWS, 
     $                                Z((JLOC-1)*LLDZ+ILOC), LLDZ )
                              END IF
                           END IF
                           IF( MYCOL.EQ.CSRC4 ) THEN
                              CALL INFOG2L( INDX, LKTOP+DIM1, DESCZ, 
     $                             NPROW, NPCOL, MYROW, MYCOL, ILOC, 
     $                             JLOC, RSRC, CSRC4 )
                              IF( MYROW.EQ.RSRC ) THEN
                                 DO 306 IHH = 1, LNWIN 
                                    TAU = WORK(IPU+(IHH-1)*3)
                                    WORK(IPU+(IHH-1)*3) = ONE
                                    CALL DLARFX ( 'Right', ZROWS, 
     $                                   MIN(3,LNWIN-IHH+1),
     $                                   WORK(IPU+(IHH-1)*3), TAU, 
     $                                   WORK(IPW2+(IHH-1)*ZROWS), 
     $                                   ZROWS, WORK(IPW3) )
                                    WORK(IPU+(IHH-1)*3) = TAU
 306                             CONTINUE
                                 CALL DLACPY( 'All', ZROWS, DIM4, 
     $                                WORK(IPW3+DIM1*ZROWS), ZROWS, 
     $                                Z((JLOC-1)*LLDZ+ILOC), LLDZ )
                              END IF
                           END IF
 304                    CONTINUE
                     END IF
*     
*     Update the rows of H 
*     
                     IF( DIR.EQ.1 .AND. WANTT .AND. LENRBUF.GT.0 ) THEN
                        IF( LKBOT.LT.N ) THEN
                           IF( MYROW.EQ.RSRC1.AND.MYCOL.EQ.CSRC4.AND.
     $                          MOD(LKBOT,NB).NE.0 ) THEN
                              INDX = LKBOT + 1
                              CALL INFOG2L( LKTOP, INDX, DESCH, NPROW, 
     $                             NPCOL, MYROW, MYCOL, ILOC, JLOC, 
     $                             RSRC1, CSRC4 )
                              DO 307 IHH = 1, LNWIN 
                                 TAU = WORK(IPU+(IHH-1)*3)
                                 WORK(IPU+(IHH-1)*3) = ONE
                                 CALL DLARFX( 'Left', 
     $                                MIN(3,LNWIN-IHH+1), HCOLS,
     $                                WORK(IPU+(IHH-1)*3), TAU, 
     $                                H(IPW1+IHH-1), LNWIN, 
     $                                WORK(IPW3) )
                                 WORK(IPU+(IHH-1)*3) = TAU
 307                          CONTINUE
                              CALL DLACPY( 'All', DIM1, HCOLS, 
     $                             WORK(IPW3), LNWIN, 
     $                             H((JLOC-1)*LLDH+ILOC), LLDH )
                           END IF
                           IF( MYROW.EQ.RSRC4.AND.MYCOL.EQ.CSRC4.AND.
     $                          MOD(LKBOT,NB).NE.0 ) THEN
                              INDX = LKBOT + 1
                              CALL INFOG2L( LKTOP+DIM1, INDX, DESCH, 
     $                             NPROW, NPCOL, MYROW, MYCOL, ILOC, 
     $                             JLOC, RSRC4, CSRC4 )
                              DO 308 IHH = 1, LNWIN 
                                 TAU = WORK(IPU+(IHH-1)*3)
                                 WORK(IPU+(IHH-1)*3) = ONE
                                 CALL DLARFX( 'Left', 
     $                                MIN(3,LNWIN-IHH+1), HCOLS,
     $                                WORK(IPU+(IHH-1)*3), TAU, 
     $                                H(IPW1+IHH-1), LNWIN, 
     $                                WORK(IPW3) )
                                 WORK(IPU+(IHH-1)*3) = TAU
 308                          CONTINUE
                              CALL DLACPY( 'All', DIM4, HCOLS, 
     $                             WORK(IPW3+DIM1), LNWIN, 
     $                             H((JLOC-1)*LLDH+ILOC), LLDH )
                           END IF
                           IF( DEBUG ) WRITE(*,*)MYROW,MYCOL, 
     $                          'PDLAQR5: ICEIL 15'
                           INDXS = ICEIL(LKBOT,NB)*NB + 1 
                           IF( MOD(LKBOT,NB).NE.0 ) THEN
                              INDXE = MIN(N,INDXS+(NPCOL-2)*NB)
                           ELSE
                              INDXE = MIN(N,INDXS+(NPCOL-1)*NB)
                           END IF
                           DO 309 INDX = INDXS, INDXE, NB
                              IF( MYROW.EQ.RSRC1 ) THEN
                                 CALL INFOG2L( LKTOP, INDX, DESCH, 
     $                                NPROW, NPCOL, MYROW, MYCOL, ILOC, 
     $                                JLOC, RSRC1, CSRC )
                                 IF( MYCOL.EQ.CSRC ) THEN
                                    DO 311 IHH = 1, LNWIN 
                                       TAU = WORK(IPU+(IHH-1)*3)
                                       WORK(IPU+(IHH-1)*3) = ONE
                                       CALL DLARFX( 'Left', 
     $                                      MIN(3,LNWIN-IHH+1), HCOLS,
     $                                      WORK(IPU+(IHH-1)*3), TAU, 
     $                                      H(IPW1+IHH-1), LNWIN, 
     $                                      WORK(IPW3) )
                                       WORK(IPU+(IHH-1)*3) = TAU
 311                                   CONTINUE
                                    CALL DLACPY( 'All', DIM1, HCOLS, 
     $                                   WORK(IPW3), LNWIN, 
     $                                   H((JLOC-1)*LLDH+ILOC), LLDH )
                                 END IF
                              END IF
                              IF( MYROW.EQ.RSRC4 ) THEN
                                 CALL INFOG2L( LKTOP+DIM1, INDX, DESCH, 
     $                                NPROW, NPCOL, MYROW, MYCOL, ILOC, 
     $                                JLOC, RSRC4, CSRC )
                                 IF( MYCOL.EQ.CSRC ) THEN
                                    DO 312 IHH = 1, LNWIN 
                                       TAU = WORK(IPU+(IHH-1)*3)
                                       WORK(IPU+(IHH-1)*3) = ONE
                                       CALL DLARFX( 'Left', 
     $                                      MIN(3,LNWIN-IHH+1), HCOLS,
     $                                      WORK(IPU+(IHH-1)*3), TAU, 
     $                                      H(IPW1+IHH-1), LNWIN, 
     $                                      WORK(IPW3) )
                                       WORK(IPU+(IHH-1)*3) = TAU
 312                                CONTINUE
                                    CALL DLACPY( 'All', DIM4, HCOLS, 
     $                                   WORK(IPW3+DIM1), LNWIN, 
     $                                   H((JLOC-1)*LLDH+ILOC), LLDH )
                                 END IF
                              END IF
 309                       CONTINUE
                        END IF
                     END IF
                  END IF
*
*     Update window information - mark processed windows are completed
*
                  IF( DIR.EQ.2 ) THEN
                     IF( LKBOT.EQ.KBOT ) THEN
                        LKTOP = KBOT+1
                        LKBOT = KBOT+1
                        IWORK( 1+(WIN-1)*5 ) = LKTOP
                        IWORK( 2+(WIN-1)*5 ) = LKBOT
                     ELSE
                        IF( DEBUG ) WRITE(*,*)MYROW,MYCOL, 
     $                       'PDLAQR5: ICEIL 16'
                        IF( DEBUG ) WRITE(*,*)MYROW,MYCOL, 
     $                       'Lchain=',LCHAIN
                        IF( DEBUG ) WRITE(*,*)MYROW,MYCOL,  
     $                       'Old LKTOP,LKBOT=',
     $                       LKTOP,LKBOT   
                        IF( DEBUG ) WRITE(*,*)MYROW,MYCOL,  
     $                       'Old LNWIN,KBOT=',LNWIN,KBOT
                        LKTOP = MIN( LKTOP + LNWIN - LCHAIN,
     $                       MIN( KBOT, ICEIL( LKBOT, NB )*NB ) -
     $                       LCHAIN + 1 )
                        IWORK( 1+(WIN-1)*5 ) = LKTOP
                        IF( DEBUG ) WRITE(*,*)MYROW,MYCOL, 
     $                       'PDLAQR5: ICEIL 17'
                        LKBOT = MIN( MAX( LKBOT + LNWIN - LCHAIN,
     $                       LKTOP + NWIN - 1), MIN( KBOT, 
     $                       ICEIL( LKBOT, NB )*NB ) )
                        IWORK( 2+(WIN-1)*5 ) = LKBOT
                        IF(DEBUG) WRITE(*,*)MYROW,MYCOL, 
     $                       'New LKTOP,LKBOT=',
     $                       LKTOP,LKBOT 
                     END IF
                     IF( IWORK( 5+(WIN-1)*5 ).EQ.1 )
     $                    IWORK( 5+(WIN-1)*5 ) = 2
                     IWORK( 3+(WIN-1)*5 ) = RSRC4
                     IWORK( 4+(WIN-1)*5 ) = CSRC4
                  END IF
*
*     If nothing was done for the WIN:th window, all processors come
*     here and consider the next one instead
*     
 245              CONTINUE
                  IF( DEBUG .AND. CHKORTH ) THEN
                     WRITE(*,*) MYROW,MYCOL,'Barrier 1 in...'
                     CALL BLACS_BARRIER( ICTXT, 'All' )
                     WRITE(*,*) MYROW,MYCOL,'Barrier 1 out...'
                     IPW3 = MAX(1,IPW3,1+LENCBUF+LENRBUF)
                     WRITE(*,*) MYROW,MYCOL,'IPW3=',IPW3
                     ORTH = PCHKORTH( N, Z, DESCZ, WORK(IPW3), 
     $                    LWORK-IPW3+1 )
                     WRITE(*,*) MYROW,MYCOL,'SWEEP,ORTH=',SWEEP,ORTH
                     WRITE(*,*) MYROW,MYCOL,'Curr win=',WIN
                     WRITE(*,*) MYROW,MYCOL,'Barrier 2 in...'
                     CALL BLACS_BARRIER( ICTXT, 'All' )
                     WRITE(*,*) MYROW,MYCOL,'Barrier 2 out...'
                  END IF
c                  IF( DIR.EQ.2 .AND. SWEEP.EQ.5 .AND. WIN.EQ.3 .AND. 
c     $                 LKTOP.EQ.17 .AND. LKBOT.EQ.24 ) RETURN
 240           CONTINUE
 190        CONTINUE
 150     CONTINUE
 140     CONTINUE
*
*     In debugging, check orthogonality
*
         IF( DEBUG .AND. CHKORTH ) THEN
            WRITE(*,*) MYROW,MYCOL,'Chk orthogonality...'
            IPW3 = MAX(1,IPW3)
            ORTH = PCHKORTH( N, Z, DESCZ, WORK(IPW3), LWORK-IPW3+1 )
            WRITE(*,*) MYROW,MYCOL,'SWEEP,ORTH=',SWEEP,ORTH
         END IF
*     
*     Chased off bulges from first window?
*     
         IF( DEBUG ) WRITE(*,*)MYROW,MYCOL, 'PDLAQR5: IGAMX2D ICHOFF 2'
         IF( NPROCS.GT.1 )
     $        CALL IGAMX2D( ICTXT, 'All', '1-Tree', 1, 1, ICHOFF, 1, 
     $                      -1, -1, -1, -1, -1 )
*     
*     If bulges was chasen off from first window it is removed
*     
         IF( DEBUG ) WRITE(*,*)MYROW,MYCOL, 
     $        'cross border ICHOFF=',ICHOFF
         IF( ICHOFF.GT.0 ) THEN
            DO 198 WIN = 2, ANMWIN
               IWORK( 1+(WIN-2)*5 ) = IWORK( 1+(WIN-1)*5 )
               IWORK( 2+(WIN-2)*5 ) = IWORK( 2+(WIN-1)*5 )
               IWORK( 3+(WIN-2)*5 ) = IWORK( 3+(WIN-1)*5 )
               IWORK( 4+(WIN-2)*5 ) = IWORK( 4+(WIN-1)*5 ) 
 198        CONTINUE
            ANMWIN = ANMWIN - 1
            IF( DEBUG ) WRITE(*,*)MYROW,MYCOL, 'ANMWIN=',ANMWIN
            IPIW = 6+(ANMWIN-1)*5
         END IF
*     
*     If we have no more windows, return
*     
         IF( DEBUG ) WRITE(*,*)MYROW,MYCOL, 'Return if ANMWIN < 1'
         IF( ANMWIN.LT.1 ) RETURN
*     
*     Check for any more windows to bring over the border
*     
         WINFIN = 0
         DO 199 WIN = 1, ANMWIN
            WINFIN = WINFIN+IWORK( 5+(WIN-1)*5 )
 199     CONTINUE
         IF( DEBUG ) WRITE(*,*) MYROW,MYCOL,
     $        'check for more windows to bring over border'
         IF( WINFIN.LT.2*ANMWIN ) THEN
            IF( DEBUG ) CALL BLACS_BARRIER( ICTXT, 'All' )
            IF( DEBUG ) WRITE(*,*) MYROW,MYCOL,'redo 137 while-loop'
            IF( DEBUG ) CALL BLACS_BARRIER( ICTXT, 'All' )
            GO TO 137
         END IF      
*
*     Zero out process mark for each window - this is legal now when
*     the process starts over with local bulge-chasing etc.
*
         IF( DEBUG ) WRITE(*,*)MYROW,MYCOL, 'Zero out process marks...'
         DO 201 WIN = 1, ANMWIN
            IWORK( 5+(WIN-1)*5 ) = 0
 201     CONTINUE
*
      END IF
*
*     Go back to local bulge-chase and see if there is more work to do
*
      IF( DEBUG ) WRITE(*,*)MYROW,MYCOL, 
     $     'PDLAQR5: go back to beginning of while loop'
      GO TO 20
*
*     ==== End of PDLAQR5 ====
*
      END
