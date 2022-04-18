      RECURSIVE SUBROUTINE PDLAQR0( WANTT, WANTZ, N, ILO, IHI, H, 
     $     DESCH, WR, WI, ILOZ, IHIZ, Z, DESCZ, WORK, LWORK, 
     $     IWORK, LIWORK, INFO, TIMINGS )
*
*  -- ScaLAPACK auxiliary routine (version 1.7.x) --
*     University of Umea and HPC2N, Umea, Sweden
*     February 2008
*
      IMPLICIT NONE
*
*     .. Scalar Arguments ..
      INTEGER            IHI, IHIZ, ILO, ILOZ, INFO, LIWORK, LWORK, N
      LOGICAL            WANTT, WANTZ
*     ..
*     .. Array Arguments ..
      INTEGER            DESCH( * ), DESCZ( * ), IWORK( * )
      DOUBLE PRECISION   H( * ), WI( N ), WORK( * ), WR( N ),
     $                   Z( * ), TIMINGS( 10 )
*     ..
*
*  Purpose
*  =======
*
*  PDLAQR0 computes the eigenvalues of a Hessenberg matrix H
*  and, optionally, the matrices T and Z from the Schur decomposition
*  H = Z T Z**T, where T is an upper quasi-triangular matrix (the
*  Schur form), and Z is the orthogonal matrix of Schur vectors.
*
*  Optionally Z may be postmultiplied into an input orthogonal
*  matrix Q so that this routine can give the Schur factorization
*  of a matrix A which has been reduced to the Hessenberg form H
*  by the orthogonal matrix Q:  A = Q*H*Q**T = (QZ)*T*(QZ)**T.
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
*  Arguments
*  =========
*
*  WANTT   (global input) LOGICAL
*          = .TRUE. : the full Schur form T is required;
*          = .FALSE.: only eigenvalues are required.
*
*  WANTZ   (global input) LOGICAL
*          = .TRUE. : the matrix of Schur vectors Z is required;
*          = .FALSE.: Schur vectors are not required.
*
*  N       (global input) INTEGER
*          The order of the Hessenberg matrix H (and Z if WANTZ).
*          N >= 0.
*
*  ILO     (global input) INTEGER
*  IHI     (global input) INTEGER
*          It is assumed that H is already upper triangular in rows
*          and columns 1:ILO-1 and IHI+1:N. ILO and IHI are normally
*          set by a previous call to PDGEBAL, and then passed to PDGEHRD
*          when the matrix output by PDGEBAL is reduced to Hessenberg
*          form. Otherwise ILO and IHI should be set to 1 and N
*          respectively.  If N.GT.0, then 1.LE.ILO.LE.IHI.LE.N.
*          If N = 0, then ILO = 1 and IHI = 0.
*
*  H       (global input/output) DOUBLE PRECISION array, dimension
*          (DESCH(LLD_),*)
*          On entry, the upper Hessenberg matrix H.
*          On exit, if JOB = 'S', H is upper quasi-triangular in
*          rows and columns ILO:IHI, with 1-by-1 and 2-by-2 blocks on
*          the main diagonal.  The 2-by-2 diagonal blocks (corresponding 
*          to complex conjugate pairs of eigenvalues) are returned in
*          standard form, with H(i,i) = H(i+1,i+1) and
*          H(i+1,i)*H(i,i+1).LT.0. If INFO = 0 and JOB = 'E', the
*          contents of H are unspecified on exit. 
*
*  DESCH   (global and local input) INTEGER array of dimension DLEN_.
*          The array descriptor for the distributed matrix H.
*
*  WR      (global output) DOUBLE PRECISION array, dimension (N)
*  WI      (global output) DOUBLE PRECISION array, dimension (N)
*          The real and imaginary parts, respectively, of the computed
*          eigenvalues ILO to IHI are stored in the corresponding
*          elements of WR and WI. If two eigenvalues are computed as a
*          complex conjugate pair, they are stored in consecutive
*          elements of WR and WI, say the i-th and (i+1)th, with
*          WI(i) > 0 and WI(i+1) < 0. If JOB = 'S', the
*          eigenvalues are stored in the same order as on the diagonal
*          of the Schur form returned in H.
*
*  Z       (global input/output) DOUBLE PRECISION array.
*          If COMPZ = 'V', on entry Z must contain the current
*          matrix Z of accumulated transformations from, e.g., PDGEHRD, 
*          and on exit Z has been updated; transformations are applied 
*          only to the submatrix Z(ILO:IHI,ILO:IHI).
*          If COMPZ = 'N', Z is not referenced.
*          If COMPZ = 'I', on entry Z need not be set and on exit,
*          if INFO = 0, Z contains the orthogonal matrix Z of the Schur
*           vectors of H. 
*
*  DESCZ   (global and local input) INTEGER array of dimension DLEN_.
*          The array descriptor for the distributed matrix Z.
*
*  WORK    (local workspace) DOUBLE PRECISION array, dimension(DWORK)
*
*  LWORK   (local input) INTEGER
*          The length of the workspace array WORK.
*
*  IWORK   (local workspace) INTEGER array, dimension (ILWORK)
*
*  ILWORK  (local input) INTEGER
*          The length of the workspace array IWORK
*
*  INFO    (output) INTEGER
*          =    0:  successful exit
*          .LT. 0:  if INFO = -i, the i-th argument had an illegal
*                   value
*          .GT. 0:  if INFO = i, PDLAQR0 failed to compute all of
*                   the eigenvalues.  Elements 1:ilo-1 and i+1:n of WR
*                   and WI contain those eigenvalues which have been
*                   successfully computed.  (Failures are rare.)
*
*                If INFO .GT. 0 and JOB = 'E', then on exit, the
*                remaining unconverged eigenvalues are the eigen-
*                values of the upper Hessenberg matrix rows and
*                columns ILO through INFO of the final, output
*                value of H.
*
*                If INFO .GT. 0 and JOB   = 'S', then on exit
*
*           (*)  (initial value of H)*U  = U*(final value of H)
*
*                where U is an orthogonal matrix.  The final
*                value of H is upper Hessenberg and quasi-triangular
*                in rows and columns INFO+1 through IHI.
*
*                If INFO .GT. 0 and COMPZ = 'V', then on exit
*
*                  (final value of Z)  =  (initial value of Z)*U
*
*                where U is the orthogonal matrix in (*) (regard-
*                less of the value of JOB.)
*
*                If INFO .GT. 0 and COMPZ = 'I', then on exit
*                      (final value of Z)  = U
*                where U is the orthogonal matrix in (*) (regard-
*                less of the value of JOB.)
*
*                If INFO .GT. 0 and COMPZ = 'N', then Z is not
*                accessed.
*
*     ================================================================
*     Based on contributions by
*        Robert Granat, Department of Computing Science and HPC2N,
*        University of Umea, Sweden.
*
*     ================================================================
*     References:
*       K. Braman, R. Byers and R. Mathias, The Multi-Shift QR
*       Algorithm Part I: Maintaining Well Focused Shifts, and Level 3
*       Performance, SIAM Journal of Matrix Analysis, volume 23, pages
*       929--947, 2002.
*
*       K. Braman, R. Byers and R. Mathias, The Multi-Shift QR
*       Algorithm Part II: Aggressive Early Deflation, SIAM Journal
*       of Matrix Analysis, volume 23, pages 948--973, 2002.
*
*       R. Granat and D. Kressner. A new implementation of the 
*       unsymmetric QR-algorithm for ScaLAPACK. Lapack Working
*       Note XYZ, 2008.  
*
*     ================================================================
*
*     .. Parameters ..
*
*     ==== Exceptional deflation windows:  try to cure rare
*     .    slow convergence by increasing the size of the
*     .    deflation window after KEXNW iterations. =====
*
*     ==== Exceptional shifts: try to cure rare slow convergence
*     .    with ad-hoc exceptional shifts every KEXSH iterations.
*     .    The constants WILK1 and WILK2 are used to form the
*     .    exceptional shifts. ====
*
      INTEGER            BLOCK_CYCLIC_2D, CSRC_, CTXT_, DLEN_, DTYPE_,
     $                   LLD_, MB_, M_, NB_, N_, RSRC_
      LOGICAL            DEBUG, USEHQR, RECURSION, PROTOCOL,
     $                   PROFILE, DEBUG_REC
      PARAMETER          ( BLOCK_CYCLIC_2D = 1, DLEN_ = 9, DTYPE_ = 1,
     $                     CTXT_ = 2, M_ = 3, N_ = 4, MB_ = 5, NB_ = 6,
     $                     RSRC_ = 7, CSRC_ = 8, LLD_ = 9,
     $                     DEBUG = .FALSE., RECURSION = .TRUE.,
     $                     PROTOCOL = .FALSE., PROFILE = .FALSE. )
      INTEGER            NTINY
      PARAMETER          ( NTINY = 11 )
      INTEGER            KEXNW, KEXSH
      PARAMETER          ( KEXNW = 5, KEXSH = 6 )
      DOUBLE PRECISION   WILK1, WILK2
      PARAMETER          ( WILK1 = 0.75d0, WILK2 = -0.4375d0 )
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0d0, ONE = 1.0d0 )
*     ..
*     .. Local Scalars ..
      DOUBLE PRECISION   AA, BB, CC, CS, DD, SN, SS, SWAP, ELEM, T0,
     $                   ELEM1, ELEM2, ELEM3, ALPHA, SDSUM, STAMP
      INTEGER            I, INF, IT, ITMAX, K, KACC22, KBOT, KDU, KS,
     $                   KT, KTOP, KU, KV, KWH, KWTOP, KWV, LD, LS,
     $                   LWKOPT, NDFL, NH, NHO, NIBBLE, NMIN, NS, NSMAX,
     $                   NSR, NVE, NW, NWMAX, NWR, LLDH, LLDZ, II, JJ,
     $                   ICTXT, NPROW, NPCOL, MYROW, MYCOL, IPV, IPT,
     $                   IPW, IPWRK, VROWS, VCOLS, TROWS, TCOLS, WROWS,
     $                   WCOLS, HRSRC, HCSRC, NB, IS, IE, NPROCS, KK,
     $                   IROFFH, ICOFFH, HRSRC3, HCSRC3, NWIN, TOTIT,
     $                   SWEEP, JW, TOTNS, ILWKOPT
      LOGICAL            NWINC, SORTED, LQUERY
      CHARACTER          JBCMPZ*2
*     ..
*     .. External Functions ..
      DOUBLE PRECISION   MPI_WTIME
      INTEGER            PILAENVX, NUMROC, INDXG2P, ICEIL
      EXTERNAL           PILAENVX, NUMROC, INDXG2P, ICEIL, MPI_WTIME
*     ..
*     .. Local Arrays ..
      INTEGER            DESCV( DLEN_ ), DESCT( DLEN_ ), DESCW( DLEN_ )
      DOUBLE PRECISION   ZDUM( 1, 1 )
*     ..
*     .. External Subroutines ..
      EXTERNAL           PDLACPY, PDLAQR1, DLANV2, PDLAQR3, PDLAQR5, 
     $                   PDELGET
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, DBLE, INT, MAX, MIN, MOD
*     ..
*     .. Executable Statements ..
      T0 = MPI_WTIME()
      IF( DEBUG ) WRITE(*,*) MYROW,MYCOL, 
     $     'PDLAQR0 N,ILO,IHI=',N,ILO,IHI
      INFO = 0
      ICTXT = DESCH( CTXT_ )
      CALL BLACS_GRIDINFO( ICTXT, NPROW, NPCOL, MYROW, MYCOL )
      NPROCS = NPROW*NPCOL
*
*     ==== Quick return for N = 0: nothing to do. ====
*
      IF( N.EQ.0 ) THEN
         WORK( 1 ) = ONE
         RETURN
      END IF
*
*     ==== Set up job flags for PILAENV. ====
*
      IF( WANTT ) THEN
         JBCMPZ( 1: 1 ) = 'S'
      ELSE
         JBCMPZ( 1: 1 ) = 'E'
      END IF
      IF( WANTZ ) THEN
         JBCMPZ( 2: 2 ) = 'V'
      ELSE
         JBCMPZ( 2: 2 ) = 'N'
      END IF
*
*     Check if workspace query
*
      LQUERY = LWORK.EQ.-1 .OR. LIWORK.EQ.-1
*
*     Extract local leading dimensions and block factors of matrices 
*     H and Z
*
      LLDH = DESCH( LLD_ )
      LLDZ = DESCZ( LLD_ )
      NB = DESCH( MB_ )
*
*     ==== Tiny (sub-) matrices must use PDLAQR1. (Stops recursion) ====
*
      IF( N.LE.NTINY ) THEN
*     
*     ==== Estimate optimal workspace. ====
*     
         CALL PDLAQR1( WANTT, WANTZ, N, ILO, IHI, H, DESCH, WR, WI,
     $        ILOZ, IHIZ, Z, DESCZ, WORK, LWORK, IWORK, 
     $        LIWORK, INFO )
         LWKOPT = INT( WORK( 1 ) )
         TOTIT = IWORK(1)
*
*     === Completely local matrices uses LAPACK. (Stops recursion) ===
*
      ELSEIF( N.LE.NB ) THEN
         IF( MYROW.EQ.DESCH(RSRC_) .AND. MYCOL.EQ.DESCH(CSRC_) ) 
     $        CALL DLAQR0( WANTT, WANTZ, N, ILO, IHI, H, DESCH(LLD_), 
     $        WR, WI, ILOZ, IHIZ, Z, DESCZ(LLD_), WORK, LWORK, INFO ) 
*
*     === Do one more step of recursion ===
*
      ELSE
*
*        ==== Zero out iteration and sweep counters for debugging 
*        .    purposes ====
*
         TOTIT = 0
         SWEEP = 0
         TOTNS = 0
*
*        ==== Use small bulge multi-shift QR with aggressive early
*        .    deflation on larger-than-tiny matrices. ====
*
*        ==== Hope for the best. ====
*
         INFO = 0
*
*        ==== NWR = recommended deflation window size.  At this
*        .    point,  N .GT. NTINY = 11, so there is enough
*        .    subdiagonal workspace for NWR.GE.2 as required.
*        .    (In fact, there is enough subdiagonal space for
*        .    NWR.GE.3.) ====
*
         NWR = PILAENVX( ICTXT, 13, 'PDLAQR0', JBCMPZ, N, ILO, IHI, 
     $                  LWORK )
         NWR = MAX( 2, NWR )
         NWR = MIN( IHI-ILO+1, NWR )
         NW = NWR
*
*        ==== NSR = recommended number of simultaneous shifts.
*        .    At this point N .GT. NTINY = 11, so there is at
*        .    enough subdiagonal workspace for NSR to be even
*        .    and greater than or equal to two as required. ====
*
         NWIN = PILAENVX( ICTXT, 19, 'PDLAQR0', JBCMPZ, N, NB, NB, NB )
         NSR = PILAENVX( ICTXT, 15, 'PDLAQR0', JBCMPZ, N, ILO, IHI, 
     $                   MAX(NWIN,NB) )
         NSR = MIN( NSR, IHI-ILO )
         NSR = MAX( 2, NSR-MOD( NSR, 2 ) )
*
*        ==== Estimate optimal workspace ====
*
         LWKOPT = 3*ICEIL(NWR,NPROW)*ICEIL(NWR,NPCOL)
*
*        ==== Workspace query call to PDLAQR3 ====
*
         IF( DEBUG ) WRITE(*,*) MYROW,MYCOL, 
     $        'Workspace query to PDLAQR3'
         CALL PDLAQR3( WANTT, WANTZ, N, ILO, IHI, NWR+1, H, 
     $        DESCH, ILOZ, IHIZ, Z, DESCZ, LS, LD, WR, WI, H, 
     $        DESCH, N, H, DESCH, N, H, DESCH, WORK, -1, IWORK, 
     $        LIWORK, TIMINGS )
         LWKOPT = LWKOPT + INT( WORK( 1 ) )
         ILWKOPT = IWORK( 1 )
         IF( DEBUG ) WRITE(*,*) MYROW,MYCOL, 
     $        'Workspace query to PDLAQR3 done'
*
*        ==== Workspace query call to PDLAQR5 ====
*
         IF( DEBUG ) WRITE(*,*) MYROW,MYCOL, 
     $        'Workspace query to PDLAQR5'
         CALL PDLAQR5( WANTT, WANTZ, 2, N, 1, N, N, WR, WI, H, 
     $                 DESCH, ILOZ, IHIZ, Z, DESCZ, WORK, -1, IWORK, 
     $                 LIWORK, TIMINGS )
         IF( DEBUG ) WRITE(*,*) MYROW,MYCOL, 
     $        'Workspace query to PDLAQR5 done'
*
*        ==== Optimal workspace = MAX(PDLAQR3, PDLAQR5) ====
*
         LWKOPT = MAX( LWKOPT, INT( WORK( 1 ) ) )
         ILWKOPT = MAX( ILWKOPT, IWORK( 1 ) )
*
*        ==== Quick return in case of workspace query. ====
*
         IF( LQUERY ) THEN
            WORK( 1 ) = DBLE( LWKOPT )
            IWORK( 1 ) = ILWKOPT
            RETURN
         END IF
*
*        ==== DLAHQR/DLAQR0 crossover point ====
*
         NMIN = PILAENVX( ICTXT, 12, 'PDLAQR0', JBCMPZ, N, ILO, IHI, 
     $                   LWORK )
         NMIN = MAX( NTINY, NMIN )
*
*        ==== Nibble crossover point ====
*
         NIBBLE = PILAENVX( ICTXT, 14, 'PDLAQR0', JBCMPZ, N, ILO, IHI, 
     $                     LWORK )
         NIBBLE = MAX( 0, NIBBLE )
*
*        ==== Accumulate reflections during ttswp?  Use block
*        .    2-by-2 structure during matrix-matrix multiply? ====
*
         KACC22 = PILAENVX( ICTXT, 16, 'PDLAQR0', JBCMPZ, N, ILO, IHI, 
     $                     LWORK )
         KACC22 = MAX( 0, KACC22 )
         KACC22 = MIN( 2, KACC22 )
*
*        ==== NWMAX = the largest possible deflation window for
*        .    which there is sufficient workspace. ====
*
         NWMAX = MIN( ( N-1 ) / 3, LWORK / 2 )
*
*        ==== NSMAX = the Largest number of simultaneous shifts
*        .    for which there is sufficient workspace. ====
*
         NSMAX = MIN( ( N+6 ) / 9, 2*LWORK / 3 )
         NSMAX = NSMAX - MOD( NSMAX, 2 )
         IF( DEBUG) WRITE(*,*) MYROW,MYCOL, 'NSMAX =',NSMAX
*
*        ==== NDFL: an iteration count restarted at deflation. ====
*
         NDFL = 1
*
*        ==== ITMAX = iteration limit ====
*
         ITMAX = MAX( 30, 2*KEXSH )*MAX( 10, ( IHI-ILO+1 ) )
*
*        ==== Last row and column in the active block ====
*
         KBOT = IHI
*
*        ==== Main Loop ====
*
         DO 80 IT = 1, ITMAX
            TOTIT = TOTIT + 1
            IF( DEBUG ) WRITE(*,*) MYROW,MYCOL, 'ITER=',IT
*
*           ==== Done when KBOT falls below ILO ====
*
            IF( KBOT.LT.ILO )
     $         GO TO 90
*
*           ==== Locate active block ====
*
            DO 10 K = KBOT, ILO + 1, -1
               CALL INFOG2L( K, K-1, DESCH, NPROW, NPCOL, MYROW, MYCOL, 
     $                       II, JJ, HRSRC, HCSRC )
               IF( MYROW.EQ.HRSRC .AND. MYCOL.EQ.HCSRC ) THEN
                  IF( H( II + (JJ-1)*LLDH ).EQ.ZERO )
     $                 GO TO 20
               END IF
   10       CONTINUE
            K = ILO
   20       CONTINUE
            KTOP = K
            IF( DEBUG ) WRITE(*,*) MYROW,MYCOL, 
     $           'PDLAQR0: IGAMX2D KTOP=',KTOP
            IF( NPROCS.GT.1 )
     $           CALL IGAMX2D( ICTXT, 'All', '1-Tree', 1, 1, KTOP, 1, 
     $                         -1, -1, -1, -1, -1 )
*
*           ==== Select deflation window size ====
*
            NH = KBOT - KTOP + 1
            IF( DEBUG ) WRITE(*,*) MYROW,MYCOL, 'KTOP,KBOT=',KTOP,KBOT
            IF( NH.LE.NTINY ) THEN
               NW = NH
            ELSEIF( NDFL.LT.KEXNW .OR. NH.LT.NW ) THEN
*
*              ==== Typical deflation window.  If possible and
*              .    advisable, nibble the entire active block.
*              .    If not, use size NWR or NWR+1 depending upon
*              .    which has the smaller corresponding subdiagonal
*              .    entry (a heuristic). ====
*
               NWINC = .TRUE.
               IF( NH.LE.MIN( NMIN, NWMAX ) ) THEN
                  NW = NH
               ELSE
                  NW = MIN( NWR, NH, NWMAX )
                  IF( NW.LT.NWMAX ) THEN
                     IF( NW.GE.NH-1 ) THEN
                        NW = NH
                     ELSE
                        KWTOP = KBOT - NW + 1
                        CALL PDELGET( 'All', '1-Tree', ELEM1, H, KWTOP, 
     $                                KWTOP-1, DESCH )
                        CALL PDELGET( 'All', '1-Tree', ELEM2, H, 
     $                                KWTOP-1, KWTOP-2, DESCH )
                        IF( ABS( ELEM1 ).GT.ABS( ELEM2 ) ) NW = NW + 1
                     END IF
                  END IF
               END IF
            ELSE
*
*              ==== Exceptional deflation window.  If there have
*              .    been no deflations in KEXNW or more iterations,
*              .    then vary the deflation window size.   At first,
*              .    because, larger windows are, in general, more
*              .    powerful than smaller ones, rapidly increase the
*              .    window up to the maximum reasonable and possible.
*              .    Then maybe try a slightly smaller window.  ====
*
               IF( NWINC .AND. NW.LT.MIN( NWMAX, NH ) ) THEN
                  NW = MIN( NWMAX, NH, 2*NW )
               ELSE
                  NWINC = .FALSE.
                  IF( NW.EQ.NH .AND. NH.GT.2 )
     $               NW = NH - 1
               END IF
            END IF
            IF( DEBUG ) WRITE(*,*) MYROW,MYCOL, 
     $           'delfation win size NW=',NW
*
*           ==== Aggressive early deflation:
*           .    split workspace into
*           .      - an nw-by-nw work array V for orthogonal matrix
*           .      - an NW-by-at-least-NW-but-more-is-better
*           .        (NW-by-NHO) horizontal work array for Schur factor
*           .      - an at-least-NW-but-more-is-better (NVE-by-NW)
*           .        vertical work array for matrix multiplications
*           .      - align T, V and W with the deflation window
*           .        ====
*
            KV = N - NW + 1
            KT = NW + 1
            NHO = ( N-NW-1 ) - KT + 1
            KWV = NW + 2
            NVE = ( N-NW ) - KWV + 1
*
            JW = MIN( NW, KBOT-KTOP+1 )
            KWTOP = KBOT - JW + 1
            IROFFH = MOD( KWTOP - 1, NB ) 
            ICOFFH = IROFFH
            HRSRC = INDXG2P( KWTOP, NB, MYROW, DESCH(RSRC_), NPROW )
            HCSRC = INDXG2P( KWTOP, NB, MYCOL, DESCH(CSRC_), NPCOL )
            IF( DEBUG ) WRITE(*,*) MYROW,MYCOL, 
     $           'KTOP,HRSRC,HCSRC=',KTOP,HRSRC,HCSRC
            VROWS = NUMROC( JW+IROFFH, NB, MYROW, HRSRC, NPROW )
            VCOLS = NUMROC( JW+ICOFFH, NB, MYCOL, HCSRC, NPCOL )
            CALL DESCINIT( DESCV, JW+IROFFH, JW+ICOFFH, NB, NB, 
     $                     HRSRC, HCSRC, ICTXT, MAX(1, VROWS), 
     $                     INFO )
*            
            TROWS = NUMROC( JW+IROFFH, NB, MYROW, HRSRC, NPROW )
            TCOLS = NUMROC( JW+ICOFFH, NB, MYCOL, HCSRC, NPCOL )
            CALL DESCINIT( DESCT, JW+IROFFH, JW+ICOFFH, NB, NB, 
     $                     HRSRC, HCSRC, ICTXT, MAX(1, TROWS), 
     $                     INFO )
            IF( DEBUG ) WRITE(*,*) MYROW,MYCOL,
     $           'DESCT,TROWS,TCOLS=',DESCT(1:DLEN_),TROWS,TCOLS
            WROWS = NUMROC( JW+IROFFH, NB, MYROW, HRSRC, NPROW )
            WCOLS = NUMROC( JW+ICOFFH, NB, MYCOL, HCSRC, NPCOL )
            CALL DESCINIT( DESCW, JW+IROFFH, JW+ICOFFH, NB, NB, 
     $                     HRSRC, HCSRC, ICTXT, MAX(1, WROWS), 
     $                     INFO )
*
            IPV   = 1
            IPT   = IPV + DESCV( LLD_ ) * VCOLS
            IPW   = IPT + DESCT( LLD_ ) * TCOLS
            IPWRK = IPW + DESCW( LLD_ ) * WCOLS
*
*           ==== Aggressive early deflation ====
*
            IF( DEBUG ) WRITE(*,*) MYROW,MYCOL, 
     $           'PDLAQR0 -> PDLAQR3 KTOP,KBOT,IT=',
     $           KTOP, KBOT, IT
c            IF( IT.EQ.1 ) RETURN
            IF( DEBUG ) WRITE(*,*) MYROW,MYCOL, 'ITER=',IT
            IWORK(1) = IT
            STAMP = MPI_WTIME()
            CALL PDLAQR3( WANTT, WANTZ, N, KTOP, KBOT, NW, H, 
     $           DESCH, ILOZ, IHIZ, Z, DESCZ, LS, LD, WR, WI, 
     $           WORK(IPV), DESCV, NHO, WORK(IPT), DESCT, NVE, 
     $           WORK(IPW), DESCW, WORK(IPWRK), LWORK-IPWRK+1, 
     $           IWORK, LIWORK, TIMINGS )
            TIMINGS( 3 ) = TIMINGS( 3 ) + MPI_WTIME() - STAMP
            IF( DEBUG ) WRITE(*,*) MYROW,MYCOL, 
     $           'PDLAQR0 <- PDLAQR3 LS,LD,SWEEP=',
     $           LS,LD,SWEEP
c            IF( IT.EQ.1 ) RETURN
            IF( DEBUG ) WRITE(*,*) MYROW,MYCOL,'WR=',WR(1:N),'WI=',
     $           WI(1:N)
C            WRITE(*,*) MYROW,MYCOL, 'SWEEP=',SWEEP
            IF( PROTOCOL .AND. MYROW.EQ.0 .AND. MYCOL.EQ.0 )
     $           write(*, FMT='(I5,I5,I8,I8,I8,I8,I8,I8,I8,F8.2)') 
     $           IT,SWEEP,KTOP,KBOT,NW,LS,LD,100*LD,NW*NIBBLE,
     $           MPI_WTIME()-T0
            IF( PROFILE .AND. MYROW.EQ.0 .AND. MYCOL.EQ.0 )
     $           write(*, FMT='(I5,I5,F8.2,F8.2,F8.2,F8.2,
     $           F8.2,F8.2,F8.2,F8.2)') 
     $           IT,SWEEP,TIMINGS(1),TIMINGS(2),TIMINGS(3),
     $           TIMINGS(4),TIMINGS(5),TIMINGS(6),TIMINGS(7),
     $           MPI_WTIME()-T0
            IF( PROTOCOL .OR. PROFILE )
     $           CALL BLACS_BARRIER( ICTXT, 'All' )
*
*           ==== Adjust KBOT accounting for new deflations. ====
*
            KBOT = KBOT - LD
*
*           ==== KS points to the shifts. ====
*
            KS = KBOT - LS + 1
*
*           ==== Skip an expensive QR sweep if there is a (partly
*           .    heuristic) reason to expect that many eigenvalues
*           .    will deflate without it.  Here, the QR sweep is
*           .    skipped if many eigenvalues have just been deflated
*           .    or if the remaining active block is small.
*
            IF( ( LD.EQ.0 ) .OR. ( ( 100*LD.LE.NW*NIBBLE ) .AND. ( KBOT-
     $           KTOP+1.GT.MIN( NMIN, NWMAX ) ) ) ) THEN
               IF( DEBUG ) WRITE(*,*) MYROW,MYCOL, 
     $              'Perform QR sweep, ITER=',IT
*     
*              ==== NS = nominal number of simultaneous shifts.
*              .    This may be lowered (slightly) if PDLAQR3
*              .    did not provide that many shifts. ====
*
               NS = MIN( NSMAX, NSR, MAX( 2, KBOT-KTOP ) )
               IF( DEBUG) WRITE(*,*) MYROW,MYCOL, 'NS 1 =',NS
               NS = NS - MOD( NS, 2 )
               IF( DEBUG) WRITE(*,*) MYROW,MYCOL, 'NS 2 =',NS
*
*              ==== If there have been no deflations
*              .    in a multiple of KEXSH iterations,
*              .    then try exceptional shifts.
*              .    Otherwise use shifts provided by
*              .    PDLAQR3 above or from the eigenvalues
*              .    of a trailing principal submatrix. ====
*
               IF( MOD( NDFL, KEXSH ).EQ.0 ) THEN
                  IF( DEBUG ) WRITE(*,*) MYROW,MYCOL, 
     $                 '...Exceptional shifts...'
                  KS = KBOT - NS + 1
                  DO 30 I = KBOT, MAX( KS+1, KTOP+2 ), -2
                     CALL PDELGET( 'All', '1-Tree', ELEM1, H, I, I-1, 
     $                             DESCH )
                     CALL PDELGET( 'All', '1-Tree', ELEM2, H, I-1, I-2, 
     $                             DESCH )
                     CALL PDELGET( 'All', '1-Tree', ELEM3, H, I, I, 
     $                             DESCH )
                     SS = ABS( ELEM1 ) + ABS( ELEM2 )
                     AA = WILK1*SS + ELEM3
                     BB = SS
                     CC = WILK2*SS
                     DD = AA
                     CALL DLANV2( AA, BB, CC, DD, WR( I-1 ), WI( I-1 ),
     $                            WR( I ), WI( I ), CS, SN )
   30             CONTINUE
                  IF( KS.EQ.KTOP ) THEN
                     CALL PDELGET( 'All', '1-Tree', ELEM1, H, KS+1, 
     $                             KS+1, DESCH )
                     WR( KS+1 ) = ELEM1
                     WI( KS+1 ) = ZERO
                     WR( KS ) = WR( KS+1 )
                     WI( KS ) = WI( KS+1 )
                  END IF
               ELSE
*
*                 ==== Got NS/2 or fewer shifts? Use PDLAQR0 or
*                 .    PDLAQR1 on a trailing principal submatrix to
*                 .    get more. ====
*     
                  IF( DEBUG ) WRITE(*,*) MYROW,MYCOL, 
     $                 '...Check if we need more shifts...'
                  IF( KBOT-KS+1.LE.NS / 2 ) THEN
                     IF( DEBUG ) WRITE(*,*) MYROW,MYCOL, 
     $                    '...Yes, get more shifts...'
                     KS = KBOT - NS + 1
                     KT = N - NS + 1
                     IROFFH = MOD( KS - 1, NB )
                     ICOFFH = IROFFH
                     IF( NS.GT.NMIN ) THEN
                        HRSRC = INDXG2P( KS, NB, MYROW, DESCH(RSRC_), 
     $                                   NPROW )
                        HCSRC = INDXG2P( KS, NB, MYROW, DESCH(CSRC_), 
     $                                   NPCOL )
                     ELSE
                        HRSRC = 0
                        HCSRC = 0
                     END if
                     TROWS = NUMROC( NS+IROFFH, NB, MYROW, HRSRC, 
     $                               NPROW )
                     TCOLS = NUMROC( NS+ICOFFH, NB, MYCOL, HCSRC, 
     $                               NPCOL )
                     CALL DESCINIT( DESCT, NS+IROFFH, NS+ICOFFH, NB, 
     $                              NB, HRSRC, HCSRC, ICTXT, 
     $                              MAX(1, TROWS), INFO )
                     IPT = 1
                     IPWRK = IPT + DESCT(LLD_) * TCOLS
*                        
                     IF( NS.GT.NMIN .AND. RECURSION ) THEN
                        CALL PDLACPY( 'All', NS, NS, H, KS, KS, DESCH,
     $                       WORK(IPT), 1+IROFFH, 1+ICOFFH, DESCT )
                        IF( DEBUG ) WRITE(*,*) MYROW,MYCOL, 
     $                       'PDLAQR0 -> PDLAQR4'
                        STAMP = MPI_WTIME()
                        CALL PDLAQR0( .FALSE., .FALSE., IROFFH+NS, 
     $                                1+IROFFH, IROFFH+NS, WORK(IPT), 
     $                                DESCT, WR( KS-IROFFH ), 
     $                                WI( KS-IROFFH ), 1, 1, ZDUM, 
     $                                DESCZ, WORK( IPWRK ), 
     $                                LWORK-IPWRK+1, IWORK, LIWORK, 
     $                                INF, TIMINGS )
                        TIMINGS( 4 ) = TIMINGS( 4 ) + MPI_WTIME() - 
     $                       STAMP
                        IF( DEBUG ) WRITE(*,*) MYROW,MYCOL, 
     $                       'PDLAQR0 <- PDLAQR4'
                     ELSE
                        CALL PDLAMVE( 'All', NS, NS, H, KS, KS, DESCH, 
     $                       WORK(IPT), 1+IROFFH, 1+ICOFFH, DESCT, 
     $                       WORK(IPWRK) )
                        IF( DEBUG ) WRITE(*,*) MYROW,MYCOL, 
     $                       'PDLAQR0 -> PDLAQR1',
     $                       NS, 1+IROFFH, IROFFH+NS
                        STAMP = MPI_WTIME()
                        CALL PDLAQR1( .FALSE., .FALSE., IROFFH+NS, 
     $                                1+IROFFH, IROFFH+NS, WORK(IPT), 
     $                                DESCT, WR( KS-IROFFH ), 
     $                                WI( KS-IROFFH ), 1, 1, ZDUM, 
     $                                DESCZ, WORK( IPWRK ), 
     $                                LWORK-IPWRK+1, IWORK, LIWORK, 
     $                                INF )
                        TIMINGS( 1 ) = TIMINGS( 1 ) + MPI_WTIME() - 
     $                       STAMP
                        IF( DEBUG ) WRITE(*,*) MYROW,MYCOL, 
     $                       'PDLAQR0 <- PDLAQR1'
                     END IF
                     KS = KS + INF
*
*                    ==== In case of a rare QR failure use
*                    .    eigenvalues of the trailing 2-by-2
*                    .    principal submatrix.  ====
*
                     IF( KS.GE.KBOT ) THEN
                        CALL PDELGET( 'All', '1-Tree', AA, H, KBOT-1, 
     $                                KBOT-1, DESCH )
                        CALL PDELGET( 'All', '1-Tree', CC, H, KBOT, 
     $                                KBOT-1, DESCH )
                        CALL PDELGET( 'All', '1-Tree', BB, H, KBOT-1, 
     $                                KBOT, DESCH )
                        CALL PDELGET( 'All', '1-Tree', DD, H, KBOT, 
     $                                KBOT, DESCH )
                        CALL DLANV2( AA, BB, CC, DD, WR( KBOT-1 ),
     $                               WI( KBOT-1 ), WR( KBOT ),
     $                               WI( KBOT ), CS, SN )
                        KS = KBOT - 1
                     END IF
                  ELSE
                     IF( DEBUG ) WRITE(*,*) MYROW,MYCOL, 
     $                    '...No, we have enough for now...'
                  END IF
*
                  IF( KBOT-KS+1.GT.NS ) THEN
                     IF( DEBUG ) WRITE(*,*) MYROW,MYCOL, 
     $                    '...Sort shifts...'
*
*                    ==== Sort the shifts (Helps a little)
*                    .    Bubble sort keeps complex conjugate
*                    .    pairs together. ====
*
                     SORTED = .FALSE.
                     DO 50 K = KBOT, KS + 1, -1
                        IF( SORTED )
     $                     GO TO 60
                        SORTED = .true.
                        DO 40 I = KS, K - 1
                           IF( ABS( WR( I ) )+ABS( WI( I ) ).LT.
     $                         ABS( WR( I+1 ) )+ABS( WI( I+1 ) ) ) THEN
                              SORTED = .FALSE.
*
                              SWAP = WR( I )
                              WR( I ) = WR( I+1 )
                              WR( I+1 ) = SWAP
*
                              SWAP = WI( I )
                              WI( I ) = WI( I+1 )
                              WI( I+1 ) = SWAP
                           END IF
   40                   CONTINUE
   50                CONTINUE
   60                CONTINUE
                  END IF
                  
C                  IF( SWEEP.EQ.0 ) RETURN
*
*                 ==== Shuffle shifts into pairs of real shifts
*                 .    and pairs of complex conjugate shifts
*                 .    assuming complex conjugate shifts are
*                 .    already adjacent to one another. (Yes,
*                 .    they are.)  ====
*
                  DO 70 I = KBOT, KS + 2, -2
                     IF( WI( I ).NE.-WI( I-1 ) ) THEN
*
                        SWAP = WR( I )
                        WR( I ) = WR( I-1 )
                        WR( I-1 ) = WR( I-2 )
                        WR( I-2 ) = SWAP
*
                        SWAP = WI( I )
                        WI( I ) = WI( I-1 )
                        WI( I-1 ) = WI( I-2 )
                        WI( I-2 ) = SWAP
                     END IF
   70             CONTINUE
               END IF
C               IF( SWEEP.EQ.0 ) RETURN
*
*              ==== If there are only two shifts and both are
*              .    real, then use only one.  ====
*
               IF( KBOT-KS+1.EQ.2 ) THEN
                  IF( WI( KBOT ).EQ.ZERO ) THEN
                     CALL PDELGET( 'All', '1-Tree', ELEM, H, KBOT, 
     $                             KBOT, DESCH )
                     IF( ABS( WR( KBOT )-ELEM ).LT.
     $                   ABS( WR( KBOT-1 )-ELEM ) ) THEN
                        WR( KBOT-1 ) = WR( KBOT )
                     ELSE
                        WR( KBOT ) = WR( KBOT-1 )
                     END IF
                  END IF
               END IF
C               IF( SWEEP.EQ.0 ) RETURN
*
*              ==== Use up to NS of the the smallest magnatiude
*              .    shifts.  If there aren't NS shifts available,
*              .    then use them all, possibly dropping one to
*              .    make the number of shifts even. ====
*
               NS = MIN( NS, KBOT-KS+1 )
               IF( DEBUG) WRITE(*,*) MYROW,MYCOL, 'NS 3 =',NS
               NS = NS - MOD( NS, 2 )
               IF( DEBUG) WRITE(*,*) MYROW,MYCOL, 'NS 4 =',NS
               KS = KBOT - NS + 1
               IF( DEBUG ) WRITE(*,*) MYROW,MYCOL, 'KS =',KS
*
*              ==== Small-bulge multi-shift QR sweep ====
*
               TOTNS = TOTNS + NS
               IF( DEBUG ) THEN
                  SDSUM = 0.0D+00
                  DO 111 KK = 1,N-2
                     CALL PDELGET('All',' ',ALPHA,H,KK+1,KK,DESCH)
                     SDSUM = SDSUM + ABS(ALPHA)
 111              CONTINUE
               END IF
               IF( DEBUG ) WRITE(*,*) MYROW,MYCOL, 
     $              'PDLAQR0 -> PDLAQR5 NS,SDSUM=',NS,
     $              SDSUM
               SWEEP = SWEEP + 1
c               IF( SWEEP.EQ.1 ) RETURN
c               IF( DEBUG ) WRITE(*,*) MYROW,MYCOL, 'SWEEP=',SWEEP
               IWORK(1) = SWEEP
               IF( DEBUG ) WRITE(*,*) MYROW,MYCOL,'WR=',WR(1:N),'WI=',
     $           WI(1:N)
               IF( DEBUG ) WRITE(*,*) MYROW,MYCOL,'WR(KS)=',
     $              WR(KS:KS+NS-1), 'WI(KS)=',WI(KS:KS+NS-1)
               STAMP = MPI_WTIME()
               CALL PDLAQR5( WANTT, WANTZ, KACC22, N, KTOP, KBOT, 
     $              NS, WR( KS ), WI( KS ), H, DESCH, ILOZ, IHIZ, Z, 
     $              DESCZ, WORK, LWORK, IWORK, LIWORK, TIMINGS )
               TIMINGS( 5 ) = TIMINGS( 5 ) + MPI_WTIME() - STAMP
C              IF( IT.EQ.5 ) return
               IF( DEBUG ) THEN
                  SDSUM = 0.0D+00
                  DO 112 KK = 1,N-2
                     CALL PDELGET('All',' ',ALPHA,H,KK+1,KK,DESCH)
                     SDSUM = SDSUM + ABS(ALPHA)
 112              CONTINUE
               END IF
               IF( DEBUG ) WRITE(*,*) MYROW,MYCOL, 
     $              'PDLAQR0 <- PDLAQR5 SDSUM,SWEEP=',
     $              SDSUM,SWEEP
c               IF( SWEEP.EQ.1 ) RETURN
            END IF
*
*           ==== Note progress (or the lack of it). ====
*
            IF( LD.GT.0 ) THEN
               NDFL = 1
            ELSE
               NDFL = NDFL + 1
            END IF
*
*           ==== End of main loop ====
   80    CONTINUE
*
*        ==== Iteration limit exceeded.  Set INFO to show where
*        .    the problem occurred and exit. ====
*
         INFO = KBOT
   90    CONTINUE
      END IF
*
*     ==== Return the optimal value of LWORK ====
*     
      WORK( 1 ) = DBLE( LWKOPT )
      IWORK( 1 ) = ILWKOPT
      IF( .NOT. LQUERY ) THEN
         IWORK( 1 ) = TOTIT
         IWORK( 2 ) = SWEEP
         IWORK( 3 ) = TOTNS
      END IF
      RETURN
*
*     ==== End of PDLAQR0 ====
*
      END
