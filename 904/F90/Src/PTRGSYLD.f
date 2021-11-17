CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C         This file is a part of the SCASY software library.           C
C         See http://www.cs.umu.se/parallel/research/scasy             C
C         Contributors: Robert Granat and Bo Kågström.                 C
C         SCASY Version 0.10,    March 31, 2006.                       C
C         Copyright 2006, Umeå University, Sweden.                     C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE PTRGSYLD( TRANAC, TRANBD, ISGN, COMM, M, N, A, IA, JA,
     $                     DESCA, B, IB, JB, DESCB, C, IC, JC, DESCC, 
     $                     D, ID, JD, DESCD, E, IE, JE, DESCE, MB2,
     $                     DWORK, LDWORK, IWORK, LIWORK, NOEXSY, SCALE, 
     $                     INFO )
C
C  -- ScaLAPACK-style routine (preliminary version ) --
C     HPC2N and Department of Computing Science,
C     University of Umeå, Sweden.
C     Written by Robert Granat, (granat@cs.umu.se)
C     October 20, 2006.
C
      IMPLICIT NONE
C
C     .. Scalar Arguments ..
      CHARACTER*1        TRANAC, TRANBD, COMM
      INTEGER            M, N, IA, JA, IB, JB, IC, JC, ID, JD, IE, JE, 
     $                   ISGN, INFO, LDWORK, LIWORK, NOEXSY, MB2
      DOUBLE PRECISION   SCALE
C     ..
C     .. Array Arguments ..
      INTEGER            DESCA( * ), DESCB( * ), DESCC( * ), 
     $                   DESCD( * ), DESCE( * ), IWORK( * )
      DOUBLE PRECISION   A( * ), B( * ), C( * ), D( * ), E( * ),
     $                   DWORK( * )
C     ..
C
C  Purpose
C  =======
C
C  The subroutine solves the (quasi-)triangular generalized Sylvester 
C  Equation (GSYL)
C
C     op(sub(A)) * X * op(sub(B)) +/- op(sub(C)) * X * op(sub(D)) = E   
C
C  where sub(A) = A(IA:IA+M-1,JA:JA+M-1) and sub(C) = 
C  C(IC:IC+M-1,JC:JC+M-1) are M-by-M matrices, 
C  sub(B) = B(IB:IB+N-1,JB:JB+N-1) and sub(D) = D(ID:ID+N-1,JD:JD+N-1) 
C  are an N-by-N matrices and the solution X is a M-by-N matrix which 
C  overwrites sub(E) = E(IE:IE+M-1,JE:JE+N-1). 
C  
C  The notation op(_) means the transpose or non-transpose of a matrix. 
C
C  This routine should *not* be called directly, but through PGEGSYLD.
C
C  Notes
C  =====
C
C  Each global data object is described by an associated description
C  vector called DESC_.  This vector stores the information required to 
C  establish the mapping between an object element and its corresponding 
C  process and memory location. 
C
C  Let A be a generic term for any 2D block cyclicly distributed array.
C  Such a global array has an associated description vector DESCA.  In
C  the following comments, the character _ should be read as "of the
C  global array".
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
C
C
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
C  Let K be the number of rows or columns of a distributed matrix, and
C  assume that its process grid has dimension p x q.
C  LOCr( K ) denotes the number of elements of K that a process would
C  receive if K were distributed over the p processes of its process
C  column.
C  Similarly, LOCc( K ) denotes the number of elements of K that a
C  process would receive if K were distributed over the q processes of
C  its process row.
C  The values of LOCr() and LOCc() may be determined via a call to the
C  ScaLAPACK tool function, NUMROC:
C          LOCr( M ) = NUMROC( M, MB_A, MYROW, RSRC_A, NPROW ),
C          LOCc( N ) = NUMROC( N, NB_A, MYCOL, CSRC_A, NPCOL ). 
C  
C  An upper bound for these quantities may be computed by:
C          LOCr( M ) <= ceil( ceil(M/MB_A)/NPROW )*MB_A
C          LOCc( N ) <= ceil( ceil(N/NB_A)/NPCOL )*NB_A
C
C
C  Arguments
C  =========
C
C  Mode parameters
C
C  TRANAC    (global input) CHARACTER*1
C            If TRANAC = 'N' then op(A) = A and op(C) = C
C            If TRANAC = 'T' then op(A) = A**T and op(C) = C**T
C
C  TRANBD    (global input) CHARACTER*1
C            If TRANBD = 'N' then op(B) = B and op(D) = D
C            If TRANBD = 'T' then op(B) = B**T and op(D) = D**T
C
C  ISGN      (global input) INTEGER*1
C            If ISGN = 1, we solve the equation with a '+'.
C            If ISGN = -1, we solve the equation with a '-'.
C
C  Input/output arguments
C
C  COMM      (global input/output) CHARACTER*1
C            This subroutine uses two different communications schemes in
C            solving the reduced triangular problem:
C              If COMM = 'S', the "shifts" scheme is used.
C              If COMM = 'D', the "communicate on demand" scheme is used.
C            The choice COMM = 'S' is only valid for TRANAC = TRANBD = 'N' 
C            or TRANAC = TRANBD = 'T'. The scheme used will be output.
C            See the references for details.
C
C  M         (global input) INTEGER
C            The number of rows and columns of the global distributed 
C            matrices A and C. This is also the number of rows of the
C            global distributed matrix E. M >= 0.
C
C  N         (global input) INTEGER
C            The number of rows and columns of the global distributed
C            matrices B and D. This is also the number of columns of the
C            global distributed matrix E. N >= 0.
C
C  A         (local input) DOUBLE PRECISION array 
C            Array of dimension (LLD_A,LOCc(M)). Contains the local
C            pieces of the global distributed matrix A. 
C
C  IA        (global input) INTEGER
C            Row start index for sub(A), i.e., the submatrix to operate
C            on. IA >= 1.
C
C  JA        (global input) INTEGER
C            Column start index for sub(A), i.e., the submatrix to 
C            operate on. JA = IA must hold.
C
C  DESCA     (global and local input) INTEGER array of dimension DLEN_.
C            The array descriptor for the global distributed matrix A.
C
C  B         (local input) DOUBLE PRECISION array 
C            Array of dimension (LLD_B,LOCc(N)). Contains the local
C            pieces of the global distributed matrix B. 
C
C  IB        (global input) INTEGER
C            Row start index for sub(B), i.e., the submatrix to operate
C            on. IB >= 1.
C
C  JB        (global input) INTEGER
C            Column start index for sub(B), i.e., the submatrix to 
C            operate on. JB = IB must hold.
C
C  DESCB     (global and local input) INTEGER array of dimension DLEN_.
C            The array descriptor for the global distributed matrix B.
C
C  C         (local input/local output) DOUBLE PRECISION array 
C            Array of dimension (LLD_C,LOCc(M)). Contains the local 
C            pieces of the global distributed matrix C.
C
C  IC        (global input) INTEGER
C            Row start index for sub(C), i.e., the submatrix to operate
C            on. IC = IA must hold.
C
C  JC        (global input) INTEGER
C            Column start index for sub(C), i.e., the submatrix to 
C            operate on. JC = IC must hold.
C
C  DESCC     (global and local input) INTEGER array of dimension DLEN_.
C            The array descriptor for the global distributed matrix C.
C
C  D         (local input) DOUBLE PRECISION array 
C            Array of dimension (LLD_D,LOCc(N)). Contains the local
C            pieces of the global distributed matrix D.
C
C  ID        (global input) INTEGER
C            Row start index for sub(D), i.e., the submatrix to operate
C            on. ID = IB.
C
C  JD        (global input) INTEGER
C            Column start index for sub(D), i.e., the submatrix to 
C            operate on. JD = ID must hold.
C
C  DESCD     (global and local input) INTEGER array of dimension DLEN_.
C            The array descriptor for the global distributed matrix D.
C
C  E         (local input) DOUBLE PRECISION array 
C            Array of dimension (LLD_E,LOCc(N)). 
C            On entry E contains the local pieces of the global 
C            distributed matrix E. On exit, it contains the local 
C            pieces of the global distributed solution X.
C
C  IE        (global input) INTEGER
C            Row start index for sub(E), i.e., the submatrix to operate 
C            on. MOD(IE,MB_A) = MOD(IA,MB_A) must hold. 
C
C  JE        (global input) INTEGER
C            Column start index for sub(E), i.e., the submatrix to 
C            operate on. MOD(JE,MB_B) = MOD(JB,MB_B) must hold. 
C
C  DESCE     (global and local input) INTEGER array of dimension DLEN_.
C            With JOB = 'S', the array descriptor for the global 
C            distributed matrix E.
C            Otherwise, E is not referenced.
C
C  MB2       (global input) INTEGER a
C            Internal blocking factors for pipelining of subsolutions
C            for updates of the matrix E.
C            1 < = MB2 <= DESCE( MB_ ) must hold.
C
C  Workspace
C
C  DWORK     (local workspace) DOUBLE PRECISION array, dimension
C            LDWORK. 
C
C  LDWORK    (local or global input) INTEGER
C            The dimension of the array DWORK.
C            The optimal value of LDWORK is very complex and cannot
C            be expressed in a simple way here - the optimal value
C            should be calculated via a workspace query, see below.
C 
C            If LDWORK = -1, LDWORK is global input and a workspace 
C            query is assumed. The routine will then calculate the 
C            optimal workspace needed, store it in DWORK(1) and return
C            immediately. No error will then be signaled by PXERBLA.
C
C  IWORK     (global input) INTEGER array
C            Integer workspace of dimension LIWORK.
C
C  LIWORK    (global input) INTEGER
C            The dimension of the array IWORK.
C            LIWORK >= DBA + DBB + 8 * MIN( P_r, P_c ), where 
C            DBA = ICEIL(LOCr(IA+IROFFA),MB_A) and DBB = 
C            = ICEIL(LOCr(IB+IROFFB),MB_B).
C           
C  Output information
C 
C  NOEXSY    (local output) INTEGER
C            When solving the triangular problem in PTRGSYLD it is possible
C            that we have to extend some subsystems to not lose any data
C            from some 2x2 block of conjugate pairs of eigenvalues. NOEXSY 
C            helps us to keep track of the number of such extensions. 
C
C  SCALE     (global output) DOUBLE PRECISION
C            A scale factor, usually 1.0. A scale factor < 1.0 means
C            the solution may have overflowed. See INFO for details.
C
C  Error indicator
C
C  INFO      (global output) INTEGER
C             = 0:  successful exit
C             < 0:  unsuccessful exit. 
C             If the i-th argument is an array and the j-entry had
C             an illegal value, then INFO = -(i*100+j), if the i-th
C             argument is a scalar and had an illegal value, then
C             INFO = -i. If INFO = 1, we had no valid BLACS context.
C             If INFO = 2, (A,C) and (B,D) have common or very close 
C             eigenvalues; perturbed values were used to solve the 
C             equations (but (A,C) and (B,D) are unchanged).
C             If INFO = 3,  the problem is badly scaled - E is scaled
C             a factor SCALE to guarantee an overflow free solution.
C
C  Method
C  ======
C  This subroutine implements a parallel wave-front algorithm for
C  solving the triangular generalized Sylvester equation. See the 
C  references for details.
C
C  Additional requirements
C  =======================
C
C  (A,C) and (B,D) must be distributed using the same blocking factor in 
C  each direction. Moreover, for (C,F) the blocksize in the row direction 
C  must agree with (A,C)'s, and the blocksize in the column direction must
C  agree with (B,D)'s. 
C
C  The three matrix pairs (A,D), (B,E) and (C,F) must be aligned
C  internally, i.e., A must be aligned with D, and so on.
C
C  References
C  ==========
C
C  [1] Robert Granat and Bo Kågström, Parallel ScaLAPACK-style Algorithms 
C      for Standard and Generalized Sylvester-Type Matrix Equations, in 
C      preparation, Department of Computing Science and HPC2N, Umeå 
C      University, 2006.
C 
C  [2] Robert Granat and Bo Kågström, SCASY - A ScaLAPACK-style High 
C      Performance Library for Sylvester-Type Matrix Equations, in 
C      preparation, Department of Computing Science and HPC2N, Umeå 
C      University, 2006. 
C
C  See also: http://www.cs.umu.se/research/parallel/scasy
C  
C  Parallel execution recommendations
C  ==================================
C  Use a squarish process grid, if possible, for best performance.
C
C  Revisions
C  =========
C  Please report bugs to <granat@cs.umu.se>.
C
C  Keywords
C  ========
C  Wave-front algorithm, generalized Schur form, generalized Sylvester 
C  equation, explicit blocking, GEMM-updates, matrix shifts, on demand 
C
C  =====================================================================
C     ..
C     .. Parameters ..
      DOUBLE PRECISION MONE, ONE, ZERO
      PARAMETER ( MONE = -1.0D+0, ONE = 1.0D+0, ZERO = 0.0D+0 )
      INTEGER            BLOCK_CYCLIC_2D, CSRC_, CTXT_, DLEN_, DT_,
     $                   LLD_, MB_, M_, NB_, N_, RSRC_
      PARAMETER          ( BLOCK_CYCLIC_2D = 1, DLEN_ = 9, DT_ = 1,
     $                     CTXT_ = 2, M_ = 3, N_ = 4, MB_ = 5, NB_ = 6,
     $                     RSRC_ = 7, CSRC_ = 8, LLD_ = 9 )
C     ..
C     .. Local Scalars ..
      INTEGER  MB, NB, DBAC, DBBD, NROLL, IS, JS, IEND, K,
     $         MYCOL, MYROW, NPCOL, NPROW, J, NPROCS, I, IDUM, SND,
     $         ACROWS, ACCOLS, BDROWS, BDCOLS, ROWS, COLS, LINFO, JEND,
     $         IX, JX, RSRC, CSRC, LIAC, LJAC, LIBD, LJBD, INDX, ICTXT,
     $         LLDA, LLDB, LLDC, SRSRC, SCSRC, WRK, EXMEME, EXE,
     $         MWORKneeded, XIJ, ACII, BDJJ, EIJ, RRSRC, RCSRC, EXC,GIJ, 
     $         ACRSRC, ACCSRC, BDRSRC, BDCSRC, ERSRC, ECSRC, EXA, EXB,
     $         EXMEMAC, EXMEMBD, EXMEMFG, IWRK, EXACINF, EXBDINF, 
     $         EXBUFF, LEXBUFF, GI, GJ, LBI, LBJ, NBCAC, NBCBD, POS, 
     $         ACROWS2, BDCOLS2, GINDX, ACCOLS2, BDROWS2, DD, SUBS, FIJ,
     $         GFROWS, GFCOLS, SIZEFG, SIZE_SUBS, IIC, JJC, CCOLS, IIA, 
     $         EI, EJ, JSTART, Da, Db, MWORK, LMATR, NORTH, SOUTH, WEST,
     $         EAST, LLDD, LLDFG, F, G, EXD, EXF, EXG, LLDE, JJA,
     $         IROFFAC, IROFFBD, LIE, LJE, ASI, ASJ, BSI, BSJ,
     $         CSI, CSJ, LIAS, LJAS, DSI, DSJ, ESI, ESJ, LIBS, LJBS,
     $         LJCS, LICS, LIAD, LJAD, LIES, LJES, MINELEM, LIDS, LJDS,
     $         GSI, GSJ, INDXS, INDXE, INDXU, GSIND, ECOLS, FGI, FGJ,
     $         PHASE, PHASES, MNPDIM, IROWS, ICOLS, IGSI, IGSJ, IEXRW, 
     $         IEXCL, JJS, IIS, IIE, NIDEEP, IDEEPS, IDEEPE, IDEEPU, 
     $         IDEEP, XRWS, XRIND, IXRWS, IXRIND, BDUPBL, ACUPBL, ACUP, 
     $         KACUP, EKJ, BDUP, KBDUP, SIZE_IPW, ACII2, IPW, KKK 
      DOUBLE PRECISION SCALOC, USIGN, LSCALC
      LOGICAL  LQUERY, SNODE, TRANSAC, TRANSBD, EXROW, EXCOL, CEXROW, 
     $         CEXCOL, EEXCOL, EEXROW, ACEXT, BDEXT, SHIFT, SHALT1,
     $         SHALT2, SHALT3
C     ..
C     .. Local Arrays ..
      INTEGER IDUM1(1), IDUM2(1), DESCFG( DLEN_ ), DESCSA( DLEN_ ),
     $        DESCSB( DLEN_ ), DESCSC( DLEN_ ), DESCSD( DLEN_ ),
     $        DESCSE( DLEN_ ), IBUFF( 2 ) 
C     ..
C     .. External Subroutines ..
      EXTERNAL           BLACS_GRIDINFO, CHK1MAT, SB04PY, 
     $                   PXERBLA, INFOG2L, DGEMM, DLACPY,
     $                   PCHK2MAT, DSCAL, DGAMN2D,DGESD2D,
     $                   DGERV2D, IGAMX2D, DMATADD
C     ..
C     .. External Functions ..
      LOGICAL  LSAME, INT2LG
      INTEGER  NUMROC, ICEIL, ILG2NT
      EXTERNAL LSAME, NUMROC, ICEIL, INT2LG, ILG2NT  
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC          MOD, MAX, MIN 
C     ..
C     .. Executable Statements ..
C
C     Get grid parameters
C
      ICTXT = DESCE( CTXT_ )
      CALL BLACS_GRIDINFO( ICTXT, NPROW, NPCOL, MYROW, MYCOL )
C
C     Test the input parameters
C
      LINFO = 0
      INFO = 0
      IF( NPROW.EQ.-1 ) THEN
         INFO = 1
      ELSE
         CALL CHK1MAT( M, 6, M, 6, IA, JA, DESCA, 11, INFO )
         CALL CHK1MAT( N, 7, N, 7, IB, JB, DESCB, 15, INFO )
         CALL CHK1MAT( M, 6, M, 6, IC, JC, DESCC, 19, INFO )
         CALL CHK1MAT( N, 7, N, 7, ID, JD, DESCD, 23, INFO )
         CALL CHK1MAT( M, 6, N, 7, IE, JE, DESCE, 27, INFO )
         CALL PCHK2MAT( M, 6, M, 6, IA, JA, DESCA, 11, N, 7, N, 7, IB,
     $                 JB, DESCB, 15, 0, IDUM1, IDUM2, INFO )
         CALL PCHK2MAT( M, 6, M, 6, IC, JC, DESCC, 19, N, 7, N, 7, ID,
     $                 JD, DESCD, 23, 0, IDUM1, IDUM2, INFO )
         CALL PCHK2MAT( M, 6, M, 6, IA, JA, DESCA, 10, M, 6, N, 7, IE, 
     $                  JE, DESCE, 27, 0, IDUM1, IDUM2, INFO )
      END IF
C
      IF( INFO.EQ.0 ) THEN
         IF( .NOT.( LSAME( TRANAC, 'T' ) .OR. 
     $        LSAME( TRANAC, 'N' ) ) ) THEN
            INFO = -1
         END IF
      END IF
      IF( INFO.EQ.0 ) THEN
         IF( .NOT.( LSAME( TRANBD, 'T' ) .OR. 
     $        LSAME( TRANBD, 'N' ) ) ) THEN
            INFO = -2
         END IF
      END IF
C
      IF( INFO.EQ.0 ) THEN
         IF( .NOT.(ISGN.EQ.1 .OR. ISGN.EQ.-1) ) THEN
            INFO = -3
         END IF
      END IF  
C       
      IF( INFO.EQ.0 ) THEN
         IF( .NOT.( LSAME( COMM, 'S' ) .OR. 
     $        LSAME( COMM, 'D' ) ) ) THEN
            INFO = -4
         END IF
      END IF      
C
C     Check if workspace query
C
      LQUERY = LDWORK.EQ.-1 .OR. LIWORK.EQ.-1
C
C     Check what communication scheme to use - shifts will only be used
C     when the matrices A and B are both untransposed or transposed.
C
      IF( LSAME( COMM, 'S') ) THEN
         IF( .NOT. ( (LSAME( TRANAC, 'N' ).AND.LSAME( TRANBD, 'N' ) )
     $        .OR.( (LSAME( TRANAC, 'T' ).AND.LSAME( TRANBD, 'T' )))))
     $        COMM = 'D'
      END IF
C
C     Check that the last row and column of the matrices sub(A), sub(B) 
C     and sub(C) can be mapped onto the last process row or column. 
C     If not, change to on-demand communication.
C     
      MB = DESCA( MB_ )
      NB = DESCB( MB_ ) 
      IROFFAC = MOD( IA - 1, MB )
      IROFFBD = MOD( IB - 1, NB ) 
      DBAC = ICEIL( M + IROFFAC, DESCA(MB_) )
      DBBD = ICEIL( N + IROFFBD, DESCB(MB_) )
C
      IF( INFO.EQ.0 .AND. LSAME(COMM,'S')) THEN 
         DD = MAX( NPROW, NPCOL )
         IF( MIN( M + IROFFAC, N + IROFFBD ) .GT. (DD-1)**2 ) THEN
            IF( .NOT.(MOD(DBAC,DD).EQ.0 .OR. MOD(DBBD,DD).EQ.0) ) THEN
               COMM = 'D'
            END IF
         ELSE
            COMM = 'D'
         END IF
      END IF            
C     
C     Set some initial values
C
      SCALOC = ONE
      SCALE = ONE
      NPROCS = NPROW * NPCOL
      MNPDIM = MIN( NPROW, NPCOL )
      TRANSAC = LSAME( TRANAC, 'T' )
      TRANSBD = LSAME( TRANBD, 'T' )
      SHIFT = LSAME( COMM, 'S' )
      ACEXT = .FALSE.
      BDEXT = .FALSE.
      NOEXSY = 0
      LLDA = DESCA(LLD_)
      LLDB = DESCB(LLD_)
      LLDC = DESCC(LLD_)
      LLDD = DESCD(LLD_)
      LLDE = DESCE(LLD_)
      USIGN = -DBLE(ISGN)
C
C     Compute the number of blocks to keep from A and B during the 
C     updates when deep pipelining is turned on
C
      IF( MB2.NE.MB ) THEN
         ACUPBL = MAX( 1, ICEIL( DBAC-1, MNPDIM ) ) + 2
         BDUPBL = MAX( 1, ICEIL( DBBD-1, MNPDIM ) ) + 2
      ELSE
         ACUPBL = 1
         BDUPBL = 1
      END IF     
C
C     Compute the number of rows and columns held by each process and
C     the number of block columns held by each process for the smallest
C     submatrices conforming with ScaLAPACK conventions 
C
      CALL INFOG2L( IA, JA, DESCA, NPROW, NPCOL, MYROW, MYCOL, LIAC, 
     $              LJAC, ACRSRC, ACCSRC )
      CALL INFOG2L( IB, JB, DESCB, NPROW, NPCOL, MYROW, MYCOL, LIBD, 
     $              LJBD, BDRSRC, BDCSRC )
      CALL INFOG2L( IE, JE, DESCE, NPROW, NPCOL, MYROW, MYCOL, LIE, 
     $              LJE, ERSRC, ECSRC )
      ACROWS = NUMROC( M + IROFFAC, MB, MYROW, ACRSRC, NPROW )
      ACCOLS = NUMROC( M + IROFFAC, MB, MYCOL, ACCSRC, NPCOL )
      BDROWS = NUMROC( N + IROFFBD, NB, MYROW, BDRSRC, NPROW )
      BDCOLS = NUMROC( N + IROFFBD, NB, MYCOL, BDCSRC, NPCOL )
      NBCAC = ICEIL( ACCOLS, MB )
      NBCBD = ICEIL( BDCOLS, NB )
C
C     Init a matrix descriptor for matrices F and G
C
      CALL DESCINIT( DESCFG,  M+IROFFAC, N+IROFFBD, DESCE(MB_), 
     $               DESCE(NB_), ERSRC, ECSRC, ICTXT, MAX( 1, ACROWS ), 
     $               INFO )
      LLDFG = DESCFG(LLD_)
C
C     We need six extra matrix descriptors for shifts and for 
C     calculating the correct indicies in sub(A), sub(B)...
C     These descriptors are related to the smallest submatrices  
C     including sub(_) conforming with ScaLAPACK conventions.
C    
      CALL DESCINIT ( DESCSA, M+IROFFAC, M+IROFFAC, DESCA(MB_), 
     $                DESCA(NB_), ACRSRC, ACCSRC, ICTXT, LLDA, INFO )
      CALL DESCINIT ( DESCSB, N+IROFFBD, N+IROFFBD, DESCB(MB_), 
     $                DESCB(NB_), BDRSRC, BDCSRC, ICTXT, LLDB, INFO )
      CALL DESCINIT ( DESCSC, M+IROFFAC, M+IROFFAC, DESCC(MB_), 
     $                DESCC(NB_), ACRSRC, ACCSRC, ICTXT, LLDC, INFO )
      CALL DESCINIT ( DESCSD, N+IROFFBD, N+IROFFBD, DESCD(MB_), 
     $                DESCD(NB_), BDRSRC, BDCSRC, ICTXT, LLDD, INFO )
      CALL DESCINIT ( DESCSE, M+IROFFAC, N+IROFFBD, DESCE(MB_), 
     $                DESCE(NB_), ERSRC, ECSRC, ICTXT, LLDE, INFO )
C
C     Compute the needed workspace for holding the each of the
C     matrices G and F
C
      SIZEFG = ACROWS * BDCOLS
C  
C     Compute the needed memory for holding data for the extended
C     subsystems. Even if all blocks cannot be extended, we add memory
C     for all blocks to maintain indicies simple.
C
      EXMEMAC = ICEIL(ACROWS,MB)*ICEIL(ACCOLS,MB)*(2*MB+1)
      EXMEMBD = ICEIL(BDROWS,NB)*ICEIL(BDCOLS,NB)*(2*NB+1)
      EXMEME  = ICEIL(ACROWS,MB)*ICEIL(BDCOLS,NB)*(MB+NB+1)
      EXMEMFG = EXMEME
C
C     Adjust ACROWS - BDCOLS to sub(A) and sub(B) 
C
      IF( MYROW.EQ.ACRSRC ) ACROWS = ACROWS - IROFFAC
      IF( MYCOL.EQ.ACCSRC ) ACCOLS = ACCOLS - IROFFAC
      IF( MYROW.EQ.BDRSRC ) BDROWS = BDROWS - IROFFBD
      IF( MYCOL.EQ.BDCSRC ) BDCOLS = BDCOLS - IROFFBD
C
C     Compute workspace requirements for subsystem solving
C     routine
C
      SIZE_IPW = (NB+1)*(MB2+1)
C
C     Test work space
C     
      IF( INFO.EQ.0 ) THEN
C
         WRK = 2*SIZEFG + 5*(MB+1)*(NB+1) + 2*(MB+1)**2 +
     $        2*(NB+1)**2 + 2*EXMEMAC + 2*EXMEMBD + EXMEME + 
     $        2*EXMEMFG + MAX(MB, NB) + 2*ACUPBL*(MB+1)**2 +
     $        2*BDUPBL*(NB+1)**2 + SIZE_IPW
      END IF
C     
C     Compute needed integer workspace
C     
      IF( INFO.EQ.0 ) THEN
         IWRK = DBAC + DBBD + 8 * MIN( NPROW, NPCOL )
      END IF
C     
C     Check if the call to PTRGSYLD was a workspace query, if so
C     store the needed workspace in DWORK(1) and IWORK(1) and return
C     If not check if the supplied memory is big enough
C     
      IF( INFO.EQ.0 .OR. LQUERY ) THEN
         IF( WRK.GT.LDWORK .AND. .NOT.LQUERY ) THEN
            INFO = -29
         ELSEIF( IWRK.GT.LIWORK .AND. .NOT. LQUERY ) THEN
            INFO = -31
         ELSEIF( LQUERY ) THEN 
            DWORK( 1 ) = WRK
            IWORK( 1 ) = IWRK
            INFO = 0
            RETURN
         END IF
      END IF 
C     
C     Check if we shall continue of interrupt
C     
      IF( INFO.NE.0 ) THEN
         CALL PXERBLA( ICTXT, 'PTRGSYLD', -INFO )
         RETURN
      END IF
C     
C     Init some local pointers into the DWORK-array
C     
      F     = 1
      G     = F + SIZEFG
      XIJ   = G + SIZEFG
      ACII  = XIJ + (MB+1) * (NB+1)
      BDJJ  = ACII + 2 * (MB+1) ** 2
      EIJ   = BDJJ + 2 * (NB+1) ** 2
      FIJ   = EIJ + (MB+1) * (NB+1)
      GIJ   = FIJ + (MB+1) * (NB+1)
      SND   = GIJ + (MB+1) * (NB+1)
      EXA   = SND + MAX( MB, NB )
      EXB   = EXA + EXMEMAC
      EXC   = EXB + EXMEMBD
      EXD   = EXC + EXMEMAC
      EXE   = EXD + EXMEMBD
      EXF   = EXE + EXMEME
      EXG   = EXF + EXMEMFG
      EKJ   = EXG + EXMEMFG
      ACII2 = EKJ + (MB+1) * (NB+1) 
      ACUP  = ACII2 + 2 * (MB+1) ** 2
      BDUP  = ACUP + 2*ACUPBL * (MB+1) ** 2
      IPW   = BDUP + 2*BDUPBL * (NB+1) ** 2
C
C     Init two local pointers into the IDWORK-array. Those should hold
C     the startindicies of two integer array describing the extensions
C     of the submatrices of (A,C) and (B,D) (and E, since it depends on 
C     the two others). 
C
      EXACINF = 1
      EXBDINF = EXACINF + DBAC
      IROWS = EXBDINF + DBBD
      ICOLS = IROWS + MNPDIM
      IGSI  = ICOLS + MNPDIM
      IGSJ  = IGSI + MNPDIM
      IEXRW = IGSJ + MNPDIM
      IEXCL = IEXRW + MNPDIM
      IXRWS = IEXCL + MNPDIM
      IXRIND = IXRWS + MNPDIM
C
C     Check for 2x2 diagonal-block-split between any blocks of op(A) or
C     op(B) and set the extensions.
C 
      CALL PDEXTCHK( M, A, IA, JA, DESCA, IWORK( EXACINF ), DBAC, ACEXT, 
     $               INFO )
C     
      CALL PDEXTCHK( N, B, IB, JB, DESCB, IWORK( EXBDINF ), DBBD, BDEXT, 
     $               INFO )
C     
C     Do an implicit redistribution of the elements in (A,C), (B,D) and 
C     E based on the extension information just set
C
      IF( ACEXT.OR.BDEXT ) 
     $     CALL PDIMPRED( 'ABC', M, N, A, IA, JA, DESCA, B, IB, JB,
     $                    DESCB, E, IE, JE, DESCE, IWORK( EXACINF ),
     $                    IWORK( EXBDINF ),  DWORK( EXA ), EXMEMAC, 
     $                    DWORK( EXB ), EXMEMBD, DWORK( EXE ), EXMEME, 
     $                    DWORK( SND ), MAX( MB, NB ), INFO )
      IF( ACEXT.OR.BDEXT )
     $     CALL PDIMPRED( 'AB_', M, N, C, IC, JC, DESCC, D, ID, JD,
     $                    DESCD, E, IE, JE, DESCC, IWORK( EXACINF ),
     $                    IWORK( EXBDINF ),  DWORK( EXC ), EXMEMAC, 
     $                    DWORK( EXD ), EXMEMBD, DWORK( EXE ), EXMEME, 
     $                    DWORK( SND ), MAX( MB, NB ), INFO )        
C     
C     Compute the number of block diagonals of C and some communication
C     directions
C     
      NROLL = DBAC + DBBD - 1
      NORTH = MOD(MYROW + NPROW - 1, NPROW)
      EAST =  MOD(MYCOL + 1, NPCOL)
      SOUTH = MOD(MYROW + 1, NPROW)
      WEST =  MOD(MYCOL + NPCOL - 1, NPCOL)
C 
C     Depending on the transposes, set some loop variables
C     
      IF( (.NOT.TRANSAC).AND.(.NOT.TRANSBD) ) THEN
         JS = 1
         IS = DBAC
         IEND = DBAC       
      ELSEIF( TRANSAC.AND.(.NOT.TRANSBD) ) THEN
         JS = 1
         IS = 1
         IEND = 1
      ELSEIF( (.NOT.TRANSAC).AND.TRANSBD ) THEN
         JS = DBBD
         IS = DBAC
         IEND = DBAC
      ELSEIF( TRANSAC.AND.TRANSBD ) THEN
         JS = DBBD
         IS = 1
         IEND = 1
      END IF
C     
      SNODE = .TRUE.
C
C     Compute the number of rounds to do in deep pipelining and
C     set looplimits depending on the transposes
C
      NIDEEP = ICEIL( MB, MB2 )
      IF( .NOT.TRANSAC  ) THEN
         IDEEPS = NIDEEP
         IDEEPE = 1
         IDEEPU = -1
      ELSE
         IDEEPS = 1
         IDEEPE = NIDEEP
         IDEEPU = 1
      END IF
C
C     Compute the local indicies where my part of sub(_) begins
C
      ASI = IA + MB * MOD( NPROW + MYROW - ACRSRC, NPROW ) - IROFFAC 
      ASJ = JA + MB * MOD( NPCOL + MYCOL - ACCSRC, NPCOL ) - IROFFAC
      BSI = IB + NB * MOD( NPROW + MYROW - BDRSRC, NPROW ) - IROFFBD
      BSJ = JB + NB * MOD( NPCOL + MYCOL - BDCSRC, NPCOL ) - IROFFBD
      CSI = IC + MB * MOD( NPROW + MYROW - ACRSRC, NPROW ) - IROFFAC
      CSJ = JC + MB * MOD( NPCOL + MYCOL - ACCSRC, NPCOL ) - IROFFAC
      DSI = ID + NB * MOD( NPROW + MYROW - BDRSRC, NPROW ) - IROFFBD 
      DSJ = JD + NB * MOD( NPCOL + MYCOL - BDCSRC, NPCOL ) - IROFFBD
      ESI = IE + MB * MOD( NPROW + MYROW - ERSRC, NPROW )  - IROFFAC
      ESJ = JE + NB * MOD( NPCOL + MYCOL - ECSRC, NPCOL )  - IROFFBD
      CALL INFOG2L( ASI, ASJ, DESCA, NPROW, NPCOL, MYROW, MYCOL, 
     $              LIAS, LJAS, IDUM1, IDUM2 )
      CALL INFOG2L( BSI, BSJ, DESCB, NPROW, NPCOL, MYROW, MYCOL, 
     $              LIBS, LJBS, IDUM1, IDUM2 )
      CALL INFOG2L( CSI, CSJ, DESCC, NPROW, NPCOL, MYROW, MYCOL, 
     $              LICS, LJCS, IDUM1, IDUM2 )
      CALL INFOG2L( DSI, DSJ, DESCD, NPROW, NPCOL, MYROW, MYCOL, 
     $              LIDS, LJDS, IDUM1, IDUM2 )
      CALL INFOG2L( ESI, ESJ, DESCE, NPROW, NPCOL, MYROW, MYCOL, 
     $              LIES, LJES, IDUM1, IDUM2 )
C
C     Compute what alternative to use in shifts 
C
      MINELEM = MIN(2*M**2 + 2*N**2,MIN(2*M**2 + 3*M*N,2*N**2 + 3*M*N))
      SHALT1  = 2*M**2 + 2*N**2 .EQ. MINELEM 
      SHALT2  = .NOT. SHALT1 .AND. (2*M**2 + 3*M*N .EQ. MINELEM)
      SHALT3  = .NOT. (SHALT1.OR.SHALT2).AND.(2*N**2 + 3*M*N.EQ.MINELEM)
C     
C     Main loop over number of diagonals in C
C     
      DO 10 K = 1, NROLL
C        
         IF( (.NOT.TRANSAC).AND.(.NOT.TRANSBD) ) THEN
            IF ( K.GT.DBBD) IS = IS - 1
         ELSEIF( TRANSAC.AND.(.NOT.TRANSBD) ) THEN
            IF ( K.GT.DBBD) IEND = IEND + 1
         ELSEIF( (.NOT.TRANSAC).AND.TRANSBD ) THEN
            IF ( K.GT.DBBD) IS = IS - 1
         ELSEIF( TRANSAC.AND.TRANSBD ) THEN
            IF ( K.GT.DBBD) IEND = IEND + 1
         END IF
C
C     If shifts are supposed to be used, we do it right here
C
         IF( SHIFT .AND. NPROCS.GT.1 .AND. 
     $     ( SHALT1.OR.N.EQ.NB.OR.M.EQ.MB ) ) THEN
C     
C     Shift (A,C) East if NPCOL > 1 and (A,C)**T West if NPCOL > 1
C     Shift (B,D) North if NPROW > 1 and (B,D)**T South if NPROW > 1
C     When shifting, also send/receive the extension elements
C
            IF ( NPCOL.GT.1.AND.M.NE.MB ) THEN
               IF( (.NOT.TRANSAC) .AND. (.NOT.TRANSBD) ) THEN
                  CALL DGESD2D( ICTXT, ACROWS, ACCOLS, 
     $                 A((LJAS-1)*LLDA+LIAS), LLDA, MYROW, 
     $                 EAST )
                  DESCSA(CSRC_) = MOD(DESCSA(CSRC_) + 1, NPCOL)
                  CALL DGERV2D( ICTXT, ACROWS, ACCOLS, 
     $                 A((LJAS-1)*LLDA+LIAS), LLDA, MYROW,
     $                 WEST )
                  CALL DGESD2D( ICTXT, ACROWS, ACCOLS, 
     $                 C((LJCS-1)*LLDC+LICS), LLDC, MYROW, 
     $                 EAST )
                  DESCSC(CSRC_) = MOD(DESCSC(CSRC_) + 1, NPCOL)
                  CALL DGERV2D( ICTXT, ACROWS, ACCOLS, 
     $                 C((LJCS-1)*LLDC+LICS), LLDC, MYROW,
     $                 WEST )
                  IF( ACEXT ) THEN
                     CALL DGESD2D( ICTXT, EXMEMAC, 1, DWORK(EXA),
     $                    EXMEMAC, MYROW, EAST )
                     CALL DGERV2D( ICTXT, EXMEMAC, 1, DWORK(EXA),
     $                    EXMEMAC, MYROW, WEST )
                     CALL DGESD2D( ICTXT, EXMEMAC, 1, DWORK(EXC),
     $                    EXMEMAC, MYROW, EAST )
                     CALL DGERV2D( ICTXT, EXMEMAC, 1, DWORK(EXC),
     $                    EXMEMAC, MYROW, WEST )
                  END IF
               ELSEIF( TRANSAC .AND. TRANSBD ) THEN
                  CALL DGESD2D( ICTXT, ACROWS, ACCOLS, 
     $                 A((LJAS-1)*LLDA+LIAS), LLDA, MYROW, 
     $                 WEST )   
                  DESCSA(CSRC_) = MOD(DESCSA(CSRC_) + NPCOL - 1, 
     $                 NPCOL )
                  CALL DGERV2D( ICTXT, ACROWS, ACCOLS,
     $                 A((LJAS-1)*LLDA+LIAS), LLDA, MYROW, 
     $                 EAST )
                  CALL DGESD2D( ICTXT, ACROWS, ACCOLS, 
     $                 C((LJCS-1)*LLDC+LICS), LLDC, MYROW, 
     $                 WEST )   
                  DESCSC(CSRC_) = MOD(DESCSC(CSRC_) + NPCOL - 1, 
     $                 NPCOL )
                  CALL DGERV2D( ICTXT, ACROWS, ACCOLS,
     $                 C((LJCS-1)*LLDC+LICS), LLDC, MYROW, 
     $                 EAST )
                  IF( ACEXT ) THEN
                     CALL DGESD2D( ICTXT, EXMEMAC, 1, DWORK(EXA), 
     $                    EXMEMAC, MYROW, WEST )
                     CALL DGERV2D( ICTXT, EXMEMAC, 1, DWORK(EXA), 
     $                    EXMEMAC, MYROW, EAST )
                     CALL DGESD2D( ICTXT, EXMEMAC, 1, DWORK(EXC), 
     $                    EXMEMAC, MYROW, WEST )
                     CALL DGERV2D( ICTXT, EXMEMAC, 1, DWORK(EXC), 
     $                    EXMEMAC, MYROW, EAST )
                  END IF
               END IF
            END IF
            IF ( NPROW.GT.1.AND.N.NE.NB ) THEN
               IF( (.NOT.TRANSAC) .AND. (.NOT.TRANSBD) ) THEN
                  CALL DGESD2D( ICTXT, BDROWS, BDCOLS,
     $                 B((LJBS-1)*LLDB+LIBS), LLDB, NORTH,
     $                 MYCOL )
                  DESCSB(RSRC_) = MOD(DESCSB(RSRC_) + NPROW - 1,
     $                 NPROW )
                  CALL DGERV2D( ICTXT, BDROWS, BDCOLS,
     $                 B((LJBS-1)*LLDB+LIBS), LLDB, SOUTH,
     $                 MYCOL )
                  CALL DGESD2D( ICTXT, BDROWS, BDCOLS,
     $                 D((LJDS-1)*LLDD+LIDS), LLDD, NORTH,
     $                 MYCOL )
                  DESCSD(RSRC_) = MOD(DESCSD(RSRC_) + NPROW - 1,
     $                 NPROW )
                  CALL DGERV2D( ICTXT, BDROWS, BDCOLS,
     $                 D((LJDS-1)*LLDD+LIDS), LLDD, SOUTH,
     $                 MYCOL )
                  IF( BDEXT ) THEN
                     CALL DGESD2D( ICTXT, EXMEMBD, 1, DWORK(EXB), 
     $                    EXMEMBD, NORTH, MYCOL )
                     CALL DGERV2D( ICTXT, EXMEMBD, 1, DWORK(EXB),
     $                    EXMEMBD, SOUTH, MYCOL )
                     CALL DGESD2D( ICTXT, EXMEMBD, 1, DWORK(EXD), 
     $                    EXMEMBD, NORTH, MYCOL )
                     CALL DGERV2D( ICTXT, EXMEMBD, 1, DWORK(EXD),
     $                    EXMEMBD, SOUTH, MYCOL )
                  END IF
               ELSEIF( TRANSAC .AND. TRANSBD ) THEN
                  CALL DGESD2D( ICTXT, BDROWS, BDCOLS,
     $                 B((LJBS-1)*LLDB+LIBS), LLDB, SOUTH,
     $                 MYCOL )
                  DESCSB(RSRC_) = MOD(DESCSB(RSRC_) + 1, NPROW)
                  CALL DGERV2D( ICTXT, BDROWS, BDCOLS,
     $                 B((LJBS-1)*LLDB+LIBS), LLDB, NORTH,
     $                 MYCOL )
                  CALL DGESD2D( ICTXT, BDROWS, BDCOLS,
     $                 D((LJDS-1)*LLDD+LIDS), LLDD, SOUTH,
     $                 MYCOL )
                  DESCSD(RSRC_) = MOD(DESCSD(RSRC_) + 1, NPROW )
                  CALL DGERV2D( ICTXT, BDROWS, BDCOLS,
     $                 D((LJDS-1)*LLDD+LIDS), LLDD, NORTH,
     $                 MYCOL )
                  IF( BDEXT ) THEN
                     CALL DGESD2D( ICTXT, EXMEMBD, 1, DWORK(EXB),
     $                    EXMEMBD, SOUTH, MYCOL )
                     CALL DGERV2D( ICTXT, EXMEMBD, 1, DWORK(EXB),
     $                    EXMEMBD, NORTH, MYCOL )
                     CALL DGESD2D( ICTXT, EXMEMBD, 1, DWORK(EXD),
     $                    EXMEMBD, SOUTH, MYCOL )
                     CALL DGERV2D( ICTXT, EXMEMBD, 1, DWORK(EXD),
     $                    EXMEMBD, NORTH, MYCOL )
                  END IF
               END IF
            END IF
         ELSEIF( SHIFT .AND. NPROCS.GT.1 .AND. SHALT2 ) THEN 
C     
C     Shift (A,C) SouthEast and E and (F,G) South  if NPROW > 1 
C     or Shift (A,C)**T NorthWest and E and (F,G) North if NPROW > 1 
C     
            IF( (.NOT.TRANSAC) .AND. (.NOT.TRANSBD) ) THEN
               CALL DGESD2D( ICTXT, ACROWS, ACCOLS,
     $              A((LJAS-1)*LLDA+LIAS), LLDA, SOUTH, 
     $              EAST )
               DESCSA(RSRC_) = MOD(DESCSA(RSRC_) + 1, NPROW)
               DESCSA(CSRC_) = MOD(DESCSA(CSRC_) + 1, NPCOL)
               CALL DGERV2D( ICTXT, ACROWS, ACCOLS, 
     $              A((LJAS-1)*LLDA+LIAS), LLDA, NORTH, 
     $              WEST )
               CALL DGESD2D( ICTXT, ACROWS, ACCOLS,
     $              C((LJCS-1)*LLDC+LICS), LLDC, SOUTH, 
     $              EAST )
               DESCSC(RSRC_) = MOD(DESCSC(RSRC_) + 1, NPROW )
               DESCSC(CSRC_) = MOD(DESCSC(CSRC_) + 1, NPCOL )
               CALL DGERV2D( ICTXT, ACROWS, ACCOLS, 
     $              C((LJCS-1)*LLDC+LICS), LLDC, NORTH, 
     $              WEST )
               IF( ACEXT ) THEN
                  CALL DGESD2D( ICTXT, EXMEMAC, 1, DWORK(EXA),
     $                 EXMEMAC, SOUTH, EAST )
                  CALL DGERV2D( ICTXT, EXMEMAC, 1, DWORK(EXA),
     $                 EXMEMAC, NORTH, WEST )
                  CALL DGESD2D( ICTXT, EXMEMAC, 1, DWORK(EXC),
     $                 EXMEMAC, SOUTH, EAST )
                  CALL DGERV2D( ICTXT, EXMEMAC, 1, DWORK(EXC),
     $                 EXMEMAC, NORTH, WEST )
               END IF
            ELSEIF( TRANSAC .AND. TRANSBD ) THEN
               CALL DGESD2D( ICTXT, ACROWS, ACCOLS,
     $              A((LJAS-1)*LLDA+LIAS), LLDA, NORTH, 
     $              WEST )
               DESCSA(RSRC_) = MOD(DESCSA(RSRC_) + NPROW - 1, NPROW)
               DESCSA(CSRC_) = MOD(DESCSA(CSRC_) + NPCOL - 1, NPCOL)
               CALL DGERV2D( ICTXT, ACROWS, ACCOLS,
     $              A((LJAS-1)*LLDA+LIAS), LLDA, SOUTH, 
     $              EAST )
               CALL DGESD2D( ICTXT, ACROWS, ACCOLS,
     $              C((LJCS-1)*LLDC+LICS), LLDC, NORTH, 
     $              WEST )
               DESCSC(RSRC_) = MOD(DESCSC(RSRC_) + NPROW - 1, NPROW)
               DESCSC(CSRC_) = MOD(DESCSC(CSRC_) + NPCOL - 1, NPCOL)
               CALL DGERV2D( ICTXT, ACROWS, ACCOLS,
     $              C((LJCS-1)*LLDC+LICS), LLDC, SOUTH, 
     $              EAST )
               IF( ACEXT ) THEN
                  CALL DGESD2D( ICTXT, EXMEMAC, 1, DWORK(EXA), 
     $                 EXMEMAC, NORTH, WEST )
                  CALL DGERV2D( ICTXT, EXMEMAC, 1, DWORK(EXA), 
     $                 EXMEMAC, SOUTH, EAST)
                  CALL DGESD2D( ICTXT, EXMEMAC, 1, DWORK(EXC), 
     $                 EXMEMAC, NORTH, WEST )
                  CALL DGERV2D( ICTXT, EXMEMAC, 1, DWORK(EXC), 
     $                 EXMEMAC, SOUTH, EAST)
               END IF
            END IF
C     
C     Shift E
C     
            IF ( NPROW.GT.1 ) THEN
               IF( (.NOT.TRANSAC) .AND. (.NOT.TRANSBD) ) THEN
                  CALL DGESD2D( ICTXT, ACROWS, BDCOLS, 
     $                 E((LJES-1)*LLDE+LIES), LLDE, SOUTH,
     $                 MYCOL )
                  DESCSE(RSRC_) = MOD(DESCSE(RSRC_) + 1, NPROW)
                  CALL DGERV2D( ICTXT, ACROWS, BDCOLS,
     $                 E((LJES-1)*LLDE+LIES), LLDE, NORTH,
     $                 MYCOL )
                  IF( ACEXT .OR. BDEXT ) THEN
                     CALL DGESD2D( ICTXT, EXMEME, 1, DWORK(EXE),
     $                    EXMEME, SOUTH, MYCOL )
                     CALL DGERV2D( ICTXT, EXMEME, 1, DWORK(EXE), 
     $                    EXMEME, NORTH, MYCOL )
                  END IF
               ELSEIF( TRANSAC .AND. TRANSBD ) THEN
                  CALL DGESD2D( ICTXT, ACROWS, BDCOLS, 
     $                 E((LJES-1)*LLDE+LIES), LLDE, NORTH,
     $                 MYCOL )
                  DESCSE(RSRC_) = MOD(DESCSE(RSRC_) + NPROW-1, NPROW)
                  CALL DGERV2D( ICTXT, ACROWS, BDCOLS,
     $                 E((LJES-1)*LLDE+LIES), LLDE, SOUTH,
     $                 MYCOL )
                  IF( ACEXT .OR. BDEXT ) THEN
                     CALL DGESD2D( ICTXT, EXMEME, 1, DWORK(EXE),
     $                    EXMEME, NORTH, MYCOL )
                     CALL DGERV2D( ICTXT, EXMEME, 1, DWORK(EXE),
     $                    EXMEME, SOUTH, MYCOL )
                  END IF
               END IF
            END IF
C     
C     Shift (F,G)
C     
            IF ( NPROW.GT.1 ) THEN
               IF( (.NOT.TRANSAC) .AND. (.NOT.TRANSBD) ) THEN
                  CALL DGESD2D( ICTXT, ACROWS, 2*BDCOLS, 
     $                 DWORK( F ), LLDFG, SOUTH,
     $                 MYCOL )
                  DESCFG(RSRC_) = MOD(DESCFG(RSRC_) + 1, NPROW)
                  CALL DGERV2D( ICTXT, ACROWS, 2*BDCOLS,
     $                 DWORK( F ), LLDFG, NORTH,
     $                 MYCOL )
                  IF( ACEXT .OR. BDEXT ) THEN
                     CALL DGESD2D( ICTXT, 2*EXMEMFG, 1, DWORK(EXF),
     $                    2*EXMEMFG, SOUTH, MYCOL )
                     CALL DGERV2D( ICTXT, 2*EXMEMFG, 1, DWORK(EXF), 
     $                    2*EXMEMFG, NORTH, MYCOL )
                  END IF
               ELSEIF( TRANSAC .AND. TRANSBD ) THEN
                  CALL DGESD2D( ICTXT, ACROWS, 2*BDCOLS, 
     $                 DWORK( F ), LLDFG, NORTH,
     $                 MYCOL )
                  DESCFG(RSRC_) = MOD(DESCFG(RSRC_) + NPROW-1,NPROW)
                  CALL DGERV2D( ICTXT, ACROWS, 2*BDCOLS,
     $                 DWORK( F ), LLDFG, SOUTH,
     $                 MYCOL )
                  IF( ACEXT .OR. BDEXT ) THEN
                     CALL DGESD2D( ICTXT, 2*EXMEMFG, 1, DWORK(EXF),
     $                    2*EXMEMFG, NORTH, MYCOL )
                     CALL DGERV2D( ICTXT, 2*EXMEMFG, 1, DWORK(EXF),
     $                    2*EXMEMFG, SOUTH, MYCOL )
                  END IF
               END IF
            END IF
         ELSEIF ( SHIFT .AND. NPROCS.GT.1 .AND. SHALT3 ) THEN 
C     
C     Shift (B,D) NorthWest and E and (F,G) West if NPCOL > 1 
C     or shift (B,D)**T SouthEast and E and (F,G) East if NPCOL > 1 
C     
            IF( (.NOT.TRANSAC) .AND. (.NOT.TRANSBD) ) THEN
               CALL DGESD2D( ICTXT, BDROWS, BDCOLS,
     $              B((LJBS-1)*LLDB+LIBS), LLDB, NORTH, 
     $              WEST )
               DESCSB(RSRC_) = MOD(DESCSB(RSRC_) + NPROW - 1, NPROW )
               DESCSB(CSRC_) = MOD(DESCSB(CSRC_) + NPCOL - 1, NPCOL )
               CALL DGERV2D( ICTXT, BDROWS, BDCOLS,
     $              B((LJBS-1)*LLDB+LIBS), LLDB, SOUTH, 
     $              EAST )
               CALL DGESD2D( ICTXT, BDROWS, BDCOLS,
     $              D((LJDS-1)*LLDD+LIDS), LLDD, NORTH, 
     $              WEST )
               DESCSD(RSRC_) = MOD(DESCSD(RSRC_) + NPROW - 1, NPROW )
               DESCSD(CSRC_) = MOD(DESCSD(CSRC_) + NPCOL - 1, NPCOL )
               CALL DGERV2D( ICTXT, BDROWS, BDCOLS,
     $              D((LJDS-1)*LLDD+LIDS), LLDD, SOUTH, 
     $              EAST )
               IF( BDEXT ) THEN
                  CALL DGESD2D( ICTXT, EXMEMBD, 1, DWORK(EXB), 
     $                 EXMEMBD, NORTH, WEST )
                  CALL DGERV2D( ICTXT, EXMEMBD, 1, DWORK(EXB), 
     $                 EXMEMBD, SOUTH, EAST )
                  CALL DGESD2D( ICTXT, EXMEMBD, 1, DWORK(EXD), 
     $                 EXMEMBD, NORTH, WEST )
                  CALL DGERV2D( ICTXT, EXMEMBD, 1, DWORK(EXD), 
     $                 EXMEMBD, SOUTH, EAST )
               END IF
            ELSEIF( TRANSAC .AND. TRANSBD ) THEN
               CALL DGESD2D( ICTXT, BDROWS, BDCOLS,
     $              B((LJBS-1)*LLDB+LIBS), LLDB, SOUTH, 
     $              EAST )
               DESCSB(RSRC_) = MOD( DESCSB(RSRC_) + 1, NPROW )
               DESCSB(CSRC_) = MOD( DESCSB(CSRC_) + 1, NPCOL )
               CALL DGERV2D( ICTXT, BDROWS, BDCOLS, 
     $              B((LJBS-1)*LLDB+LIBS), LLDB, NORTH, 
     $              WEST )
               CALL DGESD2D( ICTXT, BDROWS, BDCOLS,
     $              D((LJDS-1)*LLDD+LIDS), LLDD, SOUTH, 
     $              EAST )
               DESCSD(RSRC_) = MOD( DESCSD(RSRC_) + 1, NPROW )
               DESCSD(CSRC_) = MOD( DESCSD(CSRC_) + 1, NPCOL )
               CALL DGERV2D( ICTXT, BDROWS, BDCOLS, 
     $              D((LJDS-1)*LLDD+LIDS), LLDD, NORTH, 
     $              WEST )
               IF( BDEXT ) THEN
                  CALL DGESD2D( ICTXT, EXMEMBD, 1, DWORK(EXB), 
     $                 EXMEMBD, SOUTH, EAST )
                  CALL DGERV2D( ICTXT, EXMEMBD, 1, DWORK(EXB), 
     $                 EXMEMBD, NORTH, WEST )
                  CALL DGESD2D( ICTXT, EXMEMBD, 1, DWORK(EXD), 
     $                 EXMEMBD, SOUTH, EAST )
                  CALL DGERV2D( ICTXT, EXMEMBD, 1, DWORK(EXD), 
     $                 EXMEMBD, NORTH,WEST)
               END IF
            END IF
C     
C     Shift E
C     
            IF ( NPCOL.GT.1 ) THEN
               IF( (.NOT.TRANSAC) .AND. (.NOT.TRANSBD) ) THEN
                  CALL DGESD2D( ICTXT, ACROWS, BDCOLS,
     $                 E((LJES-1)*LLDE+LIES), LLDE, MYROW,
     $                 WEST )
                  DESCSE(CSRC_) = MOD(DESCSE(CSRC_) + NPCOL-1,NPCOL)
                  CALL DGERV2D( ICTXT, ACROWS, BDCOLS, 
     $                 E((LJES-1)*LLDE+LIES), LLDE, MYROW, 
     $                 EAST )
                  IF( ACEXT .OR. BDEXT ) THEN
                     CALL DGESD2D( ICTXT, EXMEME, 1, DWORK(EXE), 
     $                    EXMEME, MYROW, WEST )
                     CALL DGERV2D( ICTXT, EXMEME, 1, DWORK(EXE),
     $                    EXMEME, MYROW, EAST )
                  END IF
               ELSEIF( TRANSAC .AND. TRANSBD ) THEN
                  CALL DGESD2D( ICTXT, ACROWS, BDCOLS,
     $                 E((LJES-1)*LLDE+LIES), LLDE, MYROW,
     $                 EAST )
                  DESCSE(CSRC_) = MOD(DESCSE(CSRC_) + 1, NPCOL)
                  CALL DGERV2D( ICTXT, ACROWS, BDCOLS,
     $                 E((LJES-1)*LLDE+LIES), LLDE, MYROW, 
     $                 WEST )
                  IF( ACEXT .OR. BDEXT ) THEN
                     CALL DGESD2D( ICTXT, EXMEME, 1, DWORK(EXE),
     $                    EXMEME, MYROW, EAST )
                     CALL DGERV2D( ICTXT, EXMEME, 1, DWORK(EXE),
     $                    EXMEME, MYROW, WEST )
                  END IF
               END IF
            END IF
C     
C     Shift (F,G)
C     
            IF ( NPCOL.GT.1 ) THEN
               IF( (.NOT.TRANSAC) .AND. (.NOT.TRANSBD) ) THEN
                  CALL DGESD2D( ICTXT, ACROWS, 2*BDCOLS,
     $                 DWORK( F ), LLDFG, MYROW,
     $                 WEST )
                  DESCFG(CSRC_) = MOD(DESCFG(CSRC_) + NPCOL-1,NPCOL)
                  CALL DGERV2D( ICTXT, ACROWS, 2*BDCOLS, 
     $                 DWORK( F ), LLDFG, MYROW, 
     $                 EAST )
                  IF( ACEXT .OR. BDEXT ) THEN
                     CALL DGESD2D( ICTXT, 2*EXMEMFG, 1, DWORK(EXF), 
     $                    2*EXMEMFG, MYROW, WEST )
                     CALL DGERV2D( ICTXT, 2*EXMEMFG, 1, DWORK(EXF),
     $                    2*EXMEMFG, MYROW, EAST )
                  END IF
               ELSEIF( TRANSAC .AND. TRANSBD ) THEN
                  CALL DGESD2D( ICTXT, ACROWS, 2*BDCOLS,
     $                 DWORK( F ), LLDFG, MYROW,
     $                 EAST )
                  DESCFG(CSRC_) = MOD(DESCFG(CSRC_) + 1, NPCOL)
                  CALL DGERV2D( ICTXT, ACROWS, 2*BDCOLS,
     $                 DWORK( F ), LLDFG, MYROW, 
     $                 WEST )
                  IF( ACEXT .OR. BDEXT ) THEN
                     CALL DGESD2D( ICTXT, 2*EXMEMFG, 1, DWORK(EXF),
     $                    2*EXMEMFG, MYROW, EAST )
                     CALL DGERV2D( ICTXT, 2*EXMEMFG, 1, DWORK(EXF),
     $                    2*EXMEMFG, MYROW, WEST )
                  END IF
               END IF
            END IF
         END IF
C     
C     Recompute the number of block columns of (A,C) and (B,D) held by 
C     my node  
C     
         IF( SHIFT ) THEN
            ACCOLS2 = NUMROC( M+IROFFAC, MB, MYCOL, DESCA(CSRC_), NPCOL)
            BDCOLS2 = NUMROC( N+IROFFBD, NB, MYCOL, DESCB(CSRC_), NPCOL)
            NBCAC = ICEIL( ACCOLS2, MB )
            NBCBD = ICEIL( BDCOLS2, NB )
         END IF        
C
         JJS = JS 
C     
C     Solve subsystems on current diagonal in parallel - do it in 
C     a number of PHASES depending on the process configuration
C     
         PHASES = ICEIL( IS-IEND+1, MNPDIM )
         DO 20 PHASE = 1, PHASES, 1
            IIS = IS-(PHASE-1)*MNPDIM
            IIE = MAX(IS-PHASE*MNPDIM+1,IEND)
            IF( (.NOT.TRANSAC).AND.(.NOT.TRANSBD) ) THEN
               JJS = JS - (PHASE-1)*MNPDIM
            ELSEIF( TRANSAC.AND.(.NOT.TRANSBD) ) THEN
               JJS = JS + (PHASE-1)*MNPDIM
            ELSEIF( (.NOT.TRANSAC).AND.TRANSBD ) THEN
               JJS = JS + (PHASE-1)*MNPDIM
            ELSEIF( TRANSAC.AND.TRANSBD ) THEN
               JJS = JS - (PHASE-1)*MNPDIM
            END IF
            DO 24 IDEEP = IDEEPS, IDEEPE, IDEEPU
            SCALOC = ONE
            LINFO = 0
            J = JJS
            DO 30 I = IIS, IIE, -1       
C     
C     Here we check if the systems to solve are extended, and extract
C     the necessary data before communicating and calling the solving 
C     routine.
C     
C     Check if (Aii,Cii) is extended and set some variables describing 
C     Eij for this particular solve. To handle submatrices, we distinguish
C     between I = 1 and I > 1.
C     
               IF( I.EQ.1 ) THEN
                  IF( IWORK(EXACINF).EQ.0 ) THEN
                     ROWS = MIN( MB - IROFFAC, M )
                     GSI = 1 + IROFFAC
                     EXROW = .FALSE.
                  ELSEIF( IWORK(EXACINF).EQ.1 ) THEN
                     ROWS = MB - IROFFAC + 1
                     GSI = 1 + IROFFAC
                     EXROW = .TRUE.
                  END IF
               ELSE
                  IF( IWORK(EXACINF+(I-1)).EQ.0 ) THEN
                     ROWS = MIN(MB, (M - MB + IROFFAC) - (I - 2) * MB)
                     GSI = (I - 1) * MB + 1
                     EXROW = .FALSE.
                  ELSEIF( IWORK(EXACINF+(I-1)).EQ.1 ) THEN
                     ROWS = MB + 1
                     GSI = (I - 1) * MB + 1
                     EXROW = .TRUE.
                  ELSEIF( IWORK(EXACINF+(I-1)).EQ.2 ) THEN
                     ROWS = MIN(MB, (M - MB + IROFFAC) - 
     $                    (I - 2) * MB) - 1
                     GSI = (I - 1) * MB + 2
                     EXROW = .FALSE.
                  ELSEIF( IWORK(EXACINF+(I-1)).EQ.3 ) THEN
                     ROWS = MB
                     GSI = (I - 1) * MB + 2
                     EXROW = .TRUE.
                  END IF
               END IF
               IWORK( IROWS + IIS - I ) = ROWS
               IWORK( IGSI + IIS - I ) = GSI
               IWORK( IEXRW + IIS - I ) = ILG2NT(EXROW)
C     
C     Check if (Bjj,Djj) is extended and set some variables describing 
C     Eij
C     
               IF( J.EQ.1 ) THEN
                  IF( IWORK(EXBDINF).EQ.0 ) THEN
                     COLS = MIN( NB - IROFFBD, N )
                     GSJ = 1 + IROFFBD
                     EXCOL = .FALSE.
                  ELSEIF( IWORK(EXBDINF).EQ.1 ) THEN
                     COLS = NB - IROFFBD + 1
                     GSJ = 1 + IROFFBD
                     EXCOL = .TRUE.
                  END IF
               ELSE
                  IF( IWORK(EXBDINF+(J-1)).EQ.0 ) THEN
                     COLS = MIN(NB, (N - NB + IROFFBD) - (J - 2) * NB)
                     GSJ = (J - 1) * NB + 1
                     EXCOL = .FALSE.
                  ELSEIF( IWORK(EXBDINF+(J-1)).EQ.1 ) THEN
                     COLS = NB + 1
                     GSJ = (J - 1) * NB + 1
                     EXCOL = .TRUE.
                  ELSEIF( IWORK(EXBDINF+(J-1)).EQ.2 ) THEN
                     COLS = MIN(NB, (N - NB + IROFFBD) - 
     $                    (J - 2) * NB) - 1
                     GSJ = (J - 1) * NB + 2
                     EXCOL = .FALSE.
                  ELSEIF( IWORK(EXBDINF+(J-1)).EQ.3 ) THEN
                     COLS = NB 
                     GSJ = (J - 1) * NB + 2
                     EXCOL = .TRUE.
                  END IF
               END IF
               IWORK( ICOLS + IIS - I ) = COLS
               IWORK( IGSJ + IIS - I ) = GSJ
               IWORK( IEXCL + IIS - I ) = ILG2NT(EXCOL)
C
C     If some dimension is zero, skip this subsystem and the following
C     updates and go on to the next one
C
               IF( ROWS.EQ.0 .OR. COLS.EQ.0 ) GO TO 35
C     
C     Update the extended system counter
C     
               IF( EXROW .OR. EXCOL ) NOEXSY = NOEXSY + 1
C     
C     Get starting indicies and the process id:s needed etc.
C     
               CALL INFOG2L( GSI, GSI, DESCSA, NPROW, NPCOL, MYROW, 
     $              MYCOL, LIAC, LJAC, ACRSRC, ACCSRC )
               CALL INFOG2L( GSJ, GSJ, DESCSB, NPROW, NPCOL, MYROW, 
     $              MYCOL, LIBD, LJBD, BDRSRC, BDCSRC )
               CALL INFOG2L( GSI, GSJ, DESCSE, NPROW, NPCOL, MYROW, 
     $              MYCOL, IX, JX, ERSRC, ECSRC )
C
C     Submatrix construction is omitted most of the times for
C     deep pipelining
C
               IF( IDEEP.NE.IDEEPS ) GO TO 32
C     
C     Build the extended matrix (Aii,Cii) and send it to the process 
C     holding Eij
C     
               IF( MYROW.EQ.ACRSRC .AND. MYCOL.EQ.ACCSRC ) THEN
                  CALL DBEXMAT( EXROW, EXROW, ROWS, ROWS, LIAC, LJAC, 
     $                 A((LJAS-1)*LLDA+LIAS), LLDA, NBCAC, MB, MB,
     $                 DWORK( EXA ), DWORK( ACII ), MB + 1 )
                  CALL DBEXMAT( EXROW, EXROW, ROWS, ROWS, LIAC, LJAC, 
     $                 C((LJCS-1)*LLDC+LICS), LLDC, NBCAC, MB, MB, 
     $                 DWORK( EXC ), DWORK( ACII+(MB+1)*ROWS ), 
     $                 MB + 1 )
                  IF( (ACRSRC.NE.ERSRC) .OR. (ACCSRC.NE.ECSRC) ) THEN
                     CALL DGESD2D( ICTXT, ROWS, 2*ROWS, DWORK( ACII ), 
     $                    MB + 1, ERSRC, ECSRC )
                  END IF
               END IF
C     
C     Build the extended matrix (Bjj,Djj) and send it to the process 
C     holding Eij
C     
               IF( MYROW.EQ.BDRSRC .AND. MYCOL.EQ.BDCSRC ) THEN
                  CALL DBEXMAT( EXCOL, EXCOL, COLS, COLS, LIBD, LJBD, 
     $                 B((LJBS-1)*LLDB+LIBS), LLDB, NBCBD, NB, NB, 
     $                 DWORK( EXB ), DWORK( BDJJ ), NB + 1 )
                  CALL DBEXMAT( EXCOL, EXCOL, COLS, COLS, LIBD, LJBD, 
     $                 D((LJDS-1)*LLDD+LIDS), LLDD, NBCBD, NB, NB, 
     $                 DWORK( EXD ), DWORK( BDJJ+(NB+1)*COLS ), 
     $                 NB + 1 )
                  IF( (BDRSRC.NE.ERSRC) .OR. (BDCSRC.NE.ECSRC) ) THEN
                     CALL DGESD2D( ICTXT, COLS, 2*COLS, DWORK( BDJJ ),
     $                    NB + 1, ERSRC, ECSRC )
                  END IF
               END IF
C
C     
C     Receive from sends at solving node
C     
               IF( MYROW.EQ.ERSRC .AND. MYCOL.EQ.ECSRC ) THEN
                  IF( ACRSRC.NE.ERSRC .OR. ACCSRC.NE.ECSRC ) THEN
                     CALL DGERV2D( ICTXT, ROWS, 2*ROWS, 
     $                    DWORK( ACII ), MB + 1, ACRSRC, ACCSRC )   
                  END IF
                  IF( BDRSRC.NE.ERSRC .OR. BDCSRC.NE.ECSRC ) THEN
                     CALL DGERV2D( ICTXT, COLS, 2*COLS, 
     $                    DWORK( BDJJ ), NB + 1, BDRSRC, BDCSRC )
                  END IF
               END IF     
C
C     Skipped submatrix construction in deep pipelining? 
C
 32            CONTINUE  
C     
C     Set some solution variables
C     
               SNODE = MYROW.EQ.ERSRC.AND.MYCOL.EQ.ECSRC
               SRSRC = ERSRC
               SCSRC = ECSRC
C     
C     The solvingprocess: build the extended matrix Eij
C     
               IF( SNODE ) THEN
                  IF( IDEEP.EQ.IDEEPS )
     $                 CALL DBEXMAT( EXROW, EXCOL, ROWS, COLS, IX, JX, 
     $                 E((LJES-1)*LLDE+LIES), LLDE, NBCBD, MB, NB, 
     $                 DWORK( EXE ), DWORK( EIJ ), MB + 1 )
C     
C     Solve sub-system op(Aii)*Xij*op(Bjj) +/- op(Cii)*Xij*op(Djj) = Eij 
C     
                  CALL DTRGSYL( 'S', TRANAC, TRANBD, ISGN, IDEEP, 
     $                 ROWS, MB2, COLS, DWORK( ACII ), MB + 1, 
     $                 DWORK( BDJJ ), NB + 1, DWORK(ACII+ROWS*(MB+1)), 
     $                 MB + 1, DWORK(BDJJ+COLS*(NB+1)), NB + 1, 
     $                 DWORK( EIJ ), MB + 1, DWORK( IPW ), SIZE_IPW, 
     $                 XRWS, XRIND, SCALOC, LINFO )
                  IWORK( IXRWS + IIS - I ) = XRWS
                  IWORK( IXRIND + IIS - I ) = XRIND
C     
                  IF ( LINFO.NE.0 ) INFO = LINFO
C     
C     Copy the given solution back to the global matrix and update the
C     extension elements
C     
                  IF( IDEEP.EQ.IDEEPE .AND. SCALOC.EQ. ONE )
     $                 CALL DUBEXMA( EXROW, EXCOL, ROWS, COLS, IX, JX, 
     $                 E((LJES-1)*LLDE+LIES), LLDE, NBCBD, MB, NB, 
     $                 DWORK( EXE ), DWORK( EIJ ), MB + 1 )
C     
C     Warn for owerflow 
C     
                  IF (SCALOC.NE.ONE) THEN
                     INFO = 3
                  END IF
C     
C     Let the solving process copy the solution into the DWORK( XIJ ) 
C     space as a preparation for the updates of the global solution
C     
                  IF( XRWS.NE.0 )
     $                 CALL DLACPY( 'All', XRWS, COLS, 
     $                 DWORK( EIJ+XRIND-1 ), MB + 1, 
     $                 DWORK( XIJ+XRIND-1 ), MB + 1 )
               END IF
 35            CONTINUE
C     
C     Update inner loop variable j (column index)
C     
               IF( (.NOT.TRANSAC).AND.(.NOT.TRANSBD) ) THEN
                  J = J - 1
               ELSEIF( TRANSAC.AND.(.NOT.TRANSBD) ) THEN
                  J = J + 1
               ELSEIF( (.NOT.TRANSAC).AND.TRANSBD ) THEN
                  J = J + 1
               ELSEIF( TRANSAC.AND.TRANSBD ) THEN
                  J = J - 1
               END IF
C     
 30         CONTINUE
C
C     Compute global value of SCALOC by a k-to-all reduction, where 
C     k = MNPDIM. First, do a broadcast of the local SCALOC in the 
C     largest processor mesh dimension
C
            J = JJS
            DO 111 I = IIS, IIE, -1
C     
               GSI = IWORK( IGSI + IIS - I ) 
               GSJ = IWORK( IGSJ + IIS - I )
C
               CALL INFOG2L( GSI, GSJ, DESCSE, NPROW, NPCOL, 
     $              MYROW, MYCOL, IX, JX, SRSRC, SCSRC )
               SNODE = MYROW.EQ.SRSRC.AND.MYCOL.EQ.SCSRC
C
               IF( SNODE ) THEN
                  IF( NPCOL.GT.NPROW ) THEN
                     CALL DGEBS2D( ICTXT, 'Row', ' ', 1, 1, SCALOC, 1 )
                  ELSEIF( NPROW.GT.1 ) THEN
                     CALL DGEBS2D( ICTXT, 'Col', ' ', 1, 1, SCALOC, 1 )
                  END IF
                  LSCALC = SCALOC
               ELSE
                  IF( NPCOL.GT.NPROW ) THEN
                     IF( MYROW.EQ.SRSRC ) THEN
                        CALL DGEBR2D( ICTXT, 'Row', ' ', 1, 1, SCALOC, 
     $                       1, SRSRC, SCSRC )
                     END IF
                  ELSEIF( NPROW.GT.1 ) THEN
                     IF( MYCOL.EQ.SCSRC ) THEN
                        CALL DGEBR2D( ICTXT, 'Col', ' ', 1, 1, SCALOC, 
     $                       1, SRSRC, SCSRC )
                     END IF
                  END IF
               END IF
C     
               IF( (.NOT.TRANSAC).AND.(.NOT.TRANSBD) ) THEN
                  J = J - 1
               ELSEIF( TRANSAC.AND.(.NOT.TRANSBD) ) THEN
                  J = J + 1
               ELSEIF( (.NOT.TRANSAC).AND.TRANSBD ) THEN
                  J = J + 1
               ELSEIF( TRANSAC.AND.TRANSBD ) THEN
                  J = J - 1
               END IF
C
 111        CONTINUE
C
C     Then, do an all-to-one reduction in the smallest dimension,
C     collecting the result to the process number 0 in the scope
C
            IF( NPCOL.GT.NPROW ) THEN
               IF( NPROW.GT.1 ) THEN
                  CALL DGAMN2D( ICTXT, 'Col', ' ', 1, 1, SCALOC, 1, 
     $                 -1, -1, -1, 0, MYCOL )
               END IF
            ELSE
               IF( NPCOL.GT.1 ) THEN
                  CALL DGAMN2D( ICTXT, 'Row', ' ', 1, 1, SCALOC, 1, 
     $                 -1, -1, -1, MYROW, 0 )
               END IF
            END IF
C
C     Finally, do a broadcast of the computed minmum in the smallest
C     dimension using processor 0 as root
C
            IF( NPCOL.GT.NPROW ) THEN
               IF( NPROW.GT.1 ) THEN
                  IF( MYROW.EQ.0 ) THEN
                     CALL DGEBS2D( ICTXT, 'Col', ' ', 1, 1, SCALOC, 1 )
                  ELSE
                     CALL DGEBR2D( ICTXT, 'Col', ' ', 1, 1, SCALOC, 1,
     $                             0, MYCOL )
                  END IF
               END IF
            ELSE
               IF( NPCOL.GT.1 ) THEN
                  IF( MYCOL.EQ.0 ) THEN
                     CALL DGEBS2D( ICTXT, 'Row', ' ', 1, 1, SCALOC, 1 )
                  ELSE
                     CALL DGEBR2D( ICTXT, 'Row', ' ', 1, 1, SCALOC, 1,
     $                             MYROW, 0 )
                  END IF 
               END IF
            END IF
C
C     Now, do global scaling on the right hand side and do local
C     scaling of the current (anti-) diagonal blocks of the right 
C     hand side. If it is time to put the solution back into the 
C     global matrix, we do so.
C
            IF( SCALOC.NE.ONE ) THEN
C
C     Global scaling - right hand side and local extension elements
C
               DO 123 KKK = 1, N
                  CALL PDSCAL( M, SCALOC, E, IE, JE+KKK-1, DESCE, 1 )
                  CALL PDSCAL( M, SCALOC, DWORK( F ), 1, KKK, DESCFG, 
     $                 1 )
                  CALL PDSCAL( M, SCALOC, DWORK( G ), 1, KKK, DESCFG, 
     $                 1 )
 123           CONTINUE
               IF( ACEXT.OR.BDEXT ) THEN
                  CALL DSCAL( EXMEME, SCALOC, DWORK( EXE ), 1 )
                  CALL DSCAL( EXMEMFG, SCALOC, DWORK( EXF ), 1 )
                  CALL DSCAL( EXMEMFG, SCALOC, DWORK( EXG ), 1 )
               END IF
C
C     Local subsystem scaling
C
               J = JJS
               DO 222 I = IIS, IIE, -1
C     
                  ROWS = IWORK( IROWS + IIS - I )
                  GSI = IWORK( IGSI + IIS - I ) 
                  COLS = IWORK( ICOLS + IIS - I ) 
                  GSJ = IWORK( IGSJ + IIS - I )
C     
                  IF( ROWS.EQ.0 .OR. COLS.EQ.0 ) GO TO 37
C     
                  CALL INFOG2L( GSI, GSJ, DESCSE, NPROW, NPCOL, 
     $                 MYROW, MYCOL, IX, JX, SRSRC, SCSRC )     
                  SNODE = MYROW.EQ.SRSRC.AND.MYCOL.EQ.SCSRC
C
                  IF( SNODE ) THEN 
                     XRWS = IWORK( IXRWS + IIS - I )  
                     XRIND = IWORK( IXRIND + IIS - I )
                     IF( SCALOC.NE.LSCALC ) THEN
                        DO 333 KKK = 1, COLS
                           CALL DSCAL( XRWS, SCALOC / LSCALC, 
     $                          DWORK( XIJ + XRIND - 1 + 
     $                          (KKK-1)*(MB+1) ), 1 )
 333                    CONTINUE
                     END IF
                     DO 555 KKK = 1, COLS
                        CALL DSCAL( ROWS, SCALOC, 
     $                       DWORK( EIJ + (KKK-1)*(MB+1) ), 1 )
 555                 CONTINUE
                     CALL DLACPY( 'All', XRWS, COLS, 
     $                    DWORK( XIJ + XRIND - 1), MB+1,  
     $                    DWORK( EIJ + XRIND - 1), MB+1 )
                     IF( IDEEP.EQ.IDEEPE ) THEN
                        EXROW = INT2LG(IWORK( IEXRW + IIS - I )) 
                        EXCOL = INT2LG(IWORK( IEXCL + IIS - I ))
                        CALL DUBEXMA( EXROW, EXCOL, ROWS, COLS, IX, JX, 
     $                       E((LJES-1)*LLDE+LIES), LLDE, NBCBD, MB, NB, 
     $                       DWORK( EXE ), DWORK( EIJ ), MB + 1 )
                     END IF
                  END IF
C
 37               CONTINUE
C
                  IF( (.NOT.TRANSAC).AND.(.NOT.TRANSBD) ) THEN
                     J = J - 1
                  ELSEIF( TRANSAC.AND.(.NOT.TRANSBD) ) THEN
                     J = J + 1
                  ELSEIF( (.NOT.TRANSAC).AND.TRANSBD ) THEN
                     J = J + 1
                  ELSEIF( TRANSAC.AND.TRANSBD ) THEN
                     J = J - 1
                  END IF
C     
 222           CONTINUE
C     
C     Update value of SCALE according to SCALOC
C     
               SCALE = SCALE * SCALOC
            END IF
C     
C     Broadcast the local solution Xij to block row i
C
            IF( K.LT.NROLL ) THEN
               J = JJS
               DO 40 I = IIS, IIE, -1
C     
                  ROWS = IWORK( IROWS + IIS - I )
                  GSI = IWORK( IGSI + IIS - I ) 
                  COLS = IWORK( ICOLS + IIS - I ) 
                  GSJ = IWORK( IGSJ + IIS - I )
C     
                  IF( ROWS.EQ.0 .OR. COLS.EQ.0 ) GO TO 45
C     
                  CALL INFOG2L( GSI, GSJ, DESCSE, NPROW, NPCOL, MYROW, 
     $                 MYCOL, IX, JX, SRSRC, SCSRC )
C     
                  SNODE = MYROW.EQ.SRSRC.AND.MYCOL.EQ.SCSRC
C     
                  IF( SNODE ) THEN
                     XRWS = IWORK( IXRWS + IIS - I )  
                     XRIND = IWORK( IXRIND + IIS - I )
                     IBUFF(1) = XRWS
                     IBUFF(2) = XRIND
C     
                     IF (NPCOL.GT.1) THEN
                        CALL IGEBS2D( ICTXT, 'Row', ' ', 2, 
     $                       1, IBUFF, 2 )
                        IF( XRWS.GT.0 )
     $                       CALL DGEBS2D( ICTXT, 'Row', ' ', 
     $                       XRWS, COLS, DWORK(XIJ+XRIND-1), 
     $                       MB+1 )
                     END IF
                  ELSE
C     
                     IF ( NPCOL.GT.1 .AND. MYROW.EQ.SRSRC ) THEN
                        CALL IGEBR2D( ICTXT, 'Row', ' ', 2, 
     $                       1, IBUFF, 2, SRSRC, SCSRC )
                        IWORK( IXRWS + IIS - I ) = IBUFF(1)
                        IWORK( IXRIND + IIS - I ) = IBUFF(2)
                        XRWS = IWORK( IXRWS + IIS - I ) 
                        XRIND = IWORK( IXRIND + IIS - I )
                        IF( XRWS.GT.0 )
     $                       CALL DGEBR2D( ICTXT, 'Row', ' ', 
     $                       XRWS, COLS, DWORK(XIJ+XRIND-1), 
     $                       MB+1, SRSRC, SCSRC )
                     END IF
                  END IF
 45               CONTINUE
C     
C     Update inner loop variable j (column index)
C     
                  IF( (.NOT.TRANSAC).AND.(.NOT.TRANSBD) ) THEN
                     J = J - 1
                  ELSEIF( TRANSAC.AND.(.NOT.TRANSBD) ) THEN
                     J = J + 1
                  ELSEIF( (.NOT.TRANSAC).AND.TRANSBD ) THEN
                     J = J + 1
                  ELSEIF( TRANSAC.AND.TRANSBD ) THEN
                     J = J - 1
                  END IF
 40            CONTINUE
            END IF
C
C     Update subsolution in deep pipelining
C
            J = JJS
            DO 47 I = IIS, IIE, -1
C     
               ROWS = IWORK( IROWS + IIS - I )
               GSI = IWORK( IGSI + IIS - I ) 
               COLS = IWORK( ICOLS + IIS - I ) 
               GSJ = IWORK( IGSJ + IIS - I )
C     
               IF( ROWS.EQ.0 .OR. COLS.EQ.0 ) GO TO 48
C     
               CALL INFOG2L( GSI, GSJ, DESCSE, NPROW, NPCOL, 
     $              MYROW, MYCOL, IX, JX, SRSRC, SCSRC )
C     
               SNODE = MYROW.EQ.SRSRC.AND.MYCOL.EQ.SCSRC
C     
               IF( SNODE ) THEN
                  XRWS = IWORK( IXRWS + IIS - I )  
                  XRIND = IWORK( IXRIND + IIS - I )
                  IF( XRWS.GT.0 )
     $                 CALL DTRGSYL( 'U', TRANAC, TRANBD, ISGN, IDEEP, 
     $                 ROWS, MB2, COLS, DWORK( ACII ), MB + 1, 
     $                 DWORK( BDJJ ), NB + 1, DWORK(ACII+ROWS*(MB+1)), 
     $                 MB + 1, DWORK(BDJJ+COLS*(NB+1)), NB + 1, 
     $                 DWORK( EIJ ), MB + 1, DWORK( IPW ), SIZE_IPW, 
     $                 XRWS, XRIND, SCALOC, LINFO )
               END IF
 48            CONTINUE
C     
C     Update inner loop variable j (column index)
C     
               IF( (.NOT.TRANSAC).AND.(.NOT.TRANSBD) ) THEN
                  J = J - 1
               ELSEIF( TRANSAC.AND.(.NOT.TRANSBD) ) THEN
                  J = J + 1
               ELSEIF( (.NOT.TRANSAC).AND.TRANSBD ) THEN
                  J = J + 1
               ELSEIF( TRANSAC.AND.TRANSBD ) THEN
                  J = J - 1
               END IF
 47         CONTINUE
C     
C     Update rest of global system wrt to current solution
C     
C     For every update we do we also need to find out the dimensions and
C     possible extensions of the submatrices involved.
C
            IF( .NOT. TRANSBD ) THEN
               J = JJS
               DO 50 I = IIS, IIE, -1
C     
                  ROWS = IWORK( IROWS + IIS - I )
                  GSI = IWORK( IGSI + IIS - I ) 
                  EXROW = INT2LG(IWORK( IEXRW + IIS - I )) 
                  COLS = IWORK( ICOLS + IIS - I ) 
                  GSJ = IWORK( IGSJ + IIS - I )
                  EXCOL = INT2LG(IWORK( IEXCL + IIS - I ))
C     
                  IF( ROWS.EQ.0 .OR. COLS.EQ.0 ) GO TO 55
C
                  IF( MOD( I, 2 ).EQ.0 ) THEN
                     INDXS = J
                     IF( I.EQ.1 .AND. .NOT.TRANSAC ) INDXS = J + 1
                     IF( I.EQ.DBAC .AND. TRANSAC )   INDXS = J + 1
                     INDXE = DBBD
                     INDXU = 1
                  ELSE
                     INDXS = DBBD
                     INDXE = J
                     IF( I.EQ.1 .AND. .NOT.TRANSAC ) INDXE = J + 1
                     IF( I.EQ.DBAC .AND. TRANSAC )   INDXE = J + 1
                     INDXU = -1
                  END IF
                  KBDUP = 1
                  DO 60 INDX = INDXS, INDXE, INDXU
                     IF( IWORK(EXBDINF+(INDX-1)).EQ.0 ) THEN
                        BDCOLS2 = MIN(NB, (N - NB + IROFFBD) - 
     $                       (INDX - 2) * NB)
                        GSIND = (INDX - 1) * NB + 1
                        EEXCOL = .FALSE.
                     ELSEIF( IWORK(EXBDINF+(INDX-1)).EQ.1 ) THEN
                        BDCOLS2 = NB + 1
                        GSIND = (INDX - 1) * NB + 1
                        EEXCOL = .TRUE.
                     ELSEIF( IWORK(EXBDINF+(INDX-1)).EQ.2 ) THEN
                        BDCOLS2 = MIN(NB, (N - NB + IROFFBD) - 
     $                       (INDX - 2) * NB) - 1
                        GSIND = (INDX - 1) * NB + 2
                        EEXCOL = .FALSE.
                     ELSEIF( IWORK(EXBDINF+(INDX-1)).EQ.3 ) THEN
                        BDCOLS2 = NB 
                        GSIND = (INDX - 1) * NB + 2
                        EEXCOL = .TRUE.
                     END IF
C     
                     CALL INFOG2L( GSJ, GSIND, DESCSB, NPROW, NPCOL,
     $                    MYROW, MYCOL, LIBD, LJBD, RSRC, CSRC )
                     CALL INFOG2L( GSI, GSIND, DESCFG, NPROW, NPCOL, 
     $                    MYROW, MYCOL, IX, JX, RRSRC, RCSRC )
C
C     Submatrix constuction is omitted most of the times for 
C     deep pipelining
C
                     IF( IDEEP.NE.IDEEPS ) GO TO 62
                     IF( MYROW.EQ.RSRC.AND.MYCOL.EQ.CSRC ) THEN
                        IF( (RSRC.NE.RRSRC).OR.(CSRC.NE.RCSRC) ) THEN
                           CALL DBEXMAT( EXCOL, EEXCOL, COLS, BDCOLS2, 
     $                          LIBD, LJBD, B((LJBS-1)*LLDB+LIBS), LLDB, 
     $                          NBCBD, NB, NB, DWORK( EXB ), 
     $                          DWORK( BDUP+2*(BDUPBL-1)*(NB+1)**2 ), 
     $                          NB + 1 )
                           CALL DBEXMAT( EXCOL, EEXCOL, COLS, BDCOLS2,
     $                          LIBD, LJBD, D((LJDS-1)*LLDD+LIDS), 
     $                          LLDD, NBCBD, NB, NB, DWORK( EXD ), 
     $                          DWORK( BDUP+2*(BDUPBL-1)*(NB+1)**2 +
     $                          BDCOLS2*(NB+1) ), NB + 1 )
                           CALL DGESD2D( ICTXT, COLS, 2*BDCOLS2, 
     $                          DWORK(BDUP+2*(BDUPBL-1)*(NB+1)**2), 
     $                          NB + 1, RRSRC, RCSRC )
                        ELSE
                           CALL DBEXMAT( EXCOL, EEXCOL, COLS, BDCOLS2, 
     $                          LIBD, LJBD, B((LJBS-1)*LLDB+LIBS), LLDB, 
     $                          NBCBD, NB, NB, DWORK( EXB ), 
     $                          DWORK( BDUP+2*(KBDUP-1)*(NB+1)**2 ), 
     $                          NB + 1 )
                           CALL DBEXMAT( EXCOL, EEXCOL, COLS, BDCOLS2,
     $                          LIBD, LJBD, D((LJDS-1)*LLDD+LIDS), 
     $                          LLDD, NBCBD, NB, NB, DWORK( EXD ), 
     $                          DWORK( BDUP+2*(KBDUP-1)*(NB+1)**2 +
     $                          BDCOLS2*(NB+1) ), NB + 1 )
                        END IF
                     END IF
C
C     Skipped submatrix construction?
C
 62                  CONTINUE
C
                     IF(MYROW.EQ.RRSRC.AND.MYCOL.EQ.RCSRC) THEN
                        IF( J.EQ.1 .AND. IDEEP.EQ.IDEEPS ) THEN
                           CALL DLASET( 'All', MB + 1, NB + 1, ZERO, 
     $                          ZERO, DWORK(FIJ), MB + 1 )
                           CALL DLASET( 'All', MB + 1, NB + 1, ZERO, 
     $                          ZERO, DWORK(GIJ), MB + 1 )
                        ELSE
                           CALL DBEXMAT( EXROW, EEXCOL, ROWS, BDCOLS2, 
     $                          IX, JX, DWORK(F), LLDFG, NBCBD, MB, 
     $                          NB, DWORK( EXF ), DWORK(FIJ), 
     $                          MB + 1 )
                           CALL DBEXMAT( EXROW, EEXCOL, ROWS, BDCOLS2, 
     $                          IX, JX, DWORK(G), LLDFG, NBCBD, MB, 
     $                          NB, DWORK( EXG ), DWORK(GIJ), 
     $                          MB + 1 ) 
                        END IF
                        IF( IDEEP.NE.IDEEPS ) GO TO 64
                        IF( (RSRC.NE.RRSRC) .OR. (CSRC.NE.RCSRC) ) THEN
                           CALL DGERV2D( ICTXT, COLS, 2*BDCOLS2, 
     $                          DWORK( BDUP+2*(KBDUP-1)*(NB+1)**2 ), 
     $                          NB + 1, RSRC, CSRC )
                        END IF
 64                     CONTINUE
C     
C     Perform the update of Fik and Gik
C     
                        XRWS = IWORK( IXRWS + IIS - I )
                        XRIND = IWORK( IXRIND + IIS - I )
                        IF( XRWS.GT.0 .AND. BDCOLS2.GT.0 ) THEN
                           CALL DGEMM( 'N', 'N', XRWS, BDCOLS2, COLS, 
     $                          ONE, DWORK( XIJ+XRIND-1 ), MB + 1, 
     $                          DWORK( BDUP+2*(KBDUP-1)*(NB+1)**2 ),
     $                          NB + 1, ONE, DWORK(FIJ+XRIND-1), 
     $                          MB + 1 )
                           CALL DGEMM( 'N', 'N', XRWS, BDCOLS2, COLS, 
     $                          ONE, DWORK( XIJ+XRIND-1 ), MB + 1, 
     $                          DWORK( BDUP+2*(KBDUP-1)*(NB+1)**2 + 
     $                          BDCOLS2*(NB+1) ), NB + 1, ONE, 
     $                          DWORK(GIJ+XRIND-1), MB + 1 )
                        END IF
C     
C     Save the result in the matrices F and G
C     
                        CALL DUBEXMA( EXROW, EEXCOL, ROWS, BDCOLS2, IX, 
     $                       JX, DWORK(F), LLDFG, NBCBD, MB, NB, 
     $                       DWORK( EXF ), DWORK(FIJ), MB + 1 )
                        CALL DUBEXMA( EXROW, EEXCOL, ROWS, BDCOLS2, IX, 
     $                       JX, DWORK(G), LLDFG, NBCBD, MB, NB, 
     $                       DWORK( EXG ), DWORK(GIJ), MB + 1 ) 
                        IF( BDUPBL.NE.1 ) KBDUP = KBDUP + 1
                     END IF
 60               CONTINUE
C     
 55               CONTINUE
C     
C     Update inner loop variable j (column index)
C     
                  IF( .NOT.TRANSAC ) THEN
                     J = J - 1
                  ELSEIF( TRANSAC ) THEN
                     J = J + 1
                  END IF
 50            CONTINUE
C     
            ELSEIF( TRANSBD ) THEN
C     
               J = JJS
               DO 70 I = IIS, IIE, -1
C     
                  ROWS = IWORK( IROWS + IIS - I )
                  GSI = IWORK( IGSI + IIS - I ) 
                  EXROW = INT2LG(IWORK( IEXRW + IIS - I )) 
                  COLS = IWORK( ICOLS + IIS - I ) 
                  GSJ = IWORK( IGSJ + IIS - I )
                  EXCOL = INT2LG(IWORK( IEXCL + IIS - I ))
C     
                  IF( ROWS.EQ.0 .OR. COLS.EQ.0 ) GO TO 75
C     
                  IF( MOD( I, 2 ).EQ.0 ) THEN
                     INDXS = 1
                     INDXE = J
                     IF( I.EQ.1 .AND. .NOT.TRANSAC ) INDXE = J - 1
                     IF( I.EQ.DBAC .AND. TRANSAC )   INDXE = J - 1
                     INDXU = 1 
                  ELSE
                     INDXS = J
                     IF( I.EQ.1 .AND. .NOT.TRANSAC ) INDXS = J - 1
                     IF( I.EQ.DBAC .AND. TRANSAC )   INDXS = J - 1
                     INDXE = 1
                     INDXU = -1
                  END IF
                  KBDUP = 1
                  DO 80 INDX = INDXS, INDXE, INDXU
C     
C     Set some constants descibing the involved submatrices
C     
                     IF( INDX.EQ.1 ) THEN
                        IF( IWORK(EXBDINF).EQ.0 ) THEN
                           BDROWS2 = NB - IROFFBD
                           GSIND = 1 + IROFFBD
                           EEXCOL = .FALSE.
                        ELSEIF( IWORK(EXBDINF).EQ.1 ) THEN
                           BDROWS2 = NB - IROFFBD + 1
                           GSIND = 1 + IROFFBD
                           EEXCOL = .TRUE.
                        END IF
                     ELSE
                        IF( IWORK(EXBDINF+(INDX-1)).EQ.0 ) THEN
                           BDROWS2 = MIN(NB, (N - NB + IROFFBD) - 
     $                          (INDX - 2) * NB)
                           GSIND = (INDX - 1) * NB + 1
                           EEXCOL = .FALSE.
                        ELSEIF( IWORK(EXBDINF+(INDX-1)).EQ.1 ) THEN
                           BDROWS2 = NB + 1
                           GSIND = (INDX - 1) * NB + 1
                           EEXCOL = .TRUE.
                        ELSEIF( IWORK(EXBDINF+(INDX-1)).EQ.2 ) THEN
                           BDROWS2 = MIN(NB, (N - NB + IROFFBD) - 
     $                          (INDX - 2) * NB) - 1
                           GSIND = (INDX - 1) * NB + 2
                           EEXCOL = .FALSE.
                        ELSEIF( IWORK(EXBDINF+(INDX-1)).EQ.3 ) THEN
                           BDROWS2 = NB 
                           GSIND = (INDX - 1) * NB + 2
                           EEXCOL = .TRUE.
                        END IF
                     END IF
C     
                     CALL INFOG2L( GSIND, GSJ, DESCSB, NPROW, NPCOL,
     $                    MYROW, MYCOL, LIBD, LJBD, RSRC, CSRC )
                     CALL INFOG2L( GSI, GSIND, DESCFG, NPROW, NPCOL, 
     $                    MYROW, MYCOL, IX, JX, RRSRC, RCSRC )
C
C     Submatrix constuction is omitted most of the times for 
C     deep pipelining
C
                     IF( IDEEP.NE.IDEEPS ) GO TO 82
                     IF( MYROW.EQ.RSRC.AND.MYCOL.EQ.CSRC ) THEN
                        IF( (RSRC.NE.RRSRC).OR.(CSRC.NE.RCSRC) ) THEN 
                           CALL DBEXMAT( EEXCOL, EXCOL, BDROWS2, COLS, 
     $                          LIBD, LJBD, B((LJBS-1)*LLDB+LIBS), LLDB, 
     $                          NBCBD, NB, NB, DWORK( EXB ), 
     $                          DWORK( BDUP+2*(BDUPBL-1)*(NB+1)**2 ), 
     $                          NB + 1 )
                           CALL DBEXMAT( EEXCOL, EXCOL, BDROWS2, COLS,
     $                          LIBD, LJBD, D((LJDS-1)*LLDD+LIDS), 
     $                          LLDD, NBCBD, NB, NB, DWORK( EXD ), 
     $                          DWORK( BDUP+2*(BDUPBL-1)*(NB+1)**2 +
     $                          COLS*(NB+1) ), NB + 1 )
                           CALL DGESD2D( ICTXT, BDROWS2, 2*COLS, 
     $                          DWORK(BDUP+2*(BDUPBL-1)*(NB+1)**2), 
     $                          NB + 1, RRSRC, RCSRC )
                        ELSE
                           CALL DBEXMAT( EEXCOL, EXCOL, BDROWS2, COLS, 
     $                          LIBD, LJBD, B((LJBS-1)*LLDB+LIBS), LLDB, 
     $                          NBCBD, NB, NB, DWORK( EXB ), 
     $                          DWORK( BDUP+2*(KBDUP-1)*(NB+1)**2 ), 
     $                          NB + 1 )
                           CALL DBEXMAT( EEXCOL, EXCOL, BDROWS2, COLS,
     $                          LIBD, LJBD, D((LJDS-1)*LLDD+LIDS), 
     $                          LLDD, NBCBD, NB, NB, DWORK( EXD ), 
     $                          DWORK( BDUP+2*(KBDUP-1)*(NB+1)**2 +
     $                          COLS*(NB+1) ), NB + 1 )
                        END IF      
                     END IF
C
 82                  CONTINUE
C
                     IF(MYROW.EQ.RRSRC .AND. MYCOL.EQ.RCSRC) THEN
                        IF( J.EQ.DBBD .AND. IDEEP.EQ.IDEEPS ) THEN
                           CALL DLASET( 'All', MB + 1, NB + 1, ZERO, 
     $                          ZERO, DWORK(FIJ), MB + 1 )
                           CALL DLASET( 'All', MB + 1, NB + 1, ZERO, 
     $                          ZERO, DWORK(GIJ), MB + 1 )
                        ELSE
                           CALL DBEXMAT( EXROW, EEXCOL, ROWS, BDROWS2, 
     $                          IX, JX, DWORK(F), LLDFG, NBCBD, MB, 
     $                          NB, DWORK( EXF ), DWORK(FIJ), 
     $                          MB + 1 )
                           CALL DBEXMAT( EXROW, EEXCOL, ROWS, BDROWS2, 
     $                          IX, JX, DWORK(G), LLDFG, NBCBD, MB, 
     $                          NB, DWORK( EXG ), DWORK(GIJ), 
     $                          MB + 1 ) 
                        END IF     
                        IF( IDEEP.NE.IDEEPS ) GO TO 84
                        IF( (RSRC.NE.RRSRC) .OR. (CSRC.NE.RCSRC) ) THEN
                           CALL DGERV2D( ICTXT, BDROWS2, 2*COLS,
     $                          DWORK( BDUP+2*(KBDUP-1)*(NB+1)**2 ), 
     $                          NB + 1, RSRC, CSRC )
                        END IF        
 84                     CONTINUE
C     
C     Perform the update of Fik and Gik
C     
                        XRWS = IWORK( IXRWS + IIS - I )
                        XRIND = IWORK( IXRIND + IIS - I )
                        IF( XRWS.GT.0 .AND. BDROWS2.GT.0 ) THEN 
                           CALL DGEMM( 'N', 'T', XRWS, BDROWS2, COLS, 
     $                          ONE, DWORK( XIJ+XRIND-1 ), MB + 1, 
     $                          DWORK( BDUP+2*(KBDUP-1)*(NB+1)**2 ),
     $                          NB + 1, ONE, DWORK(FIJ+XRIND-1), 
     $                          MB + 1 )
                           CALL DGEMM( 'N', 'T', XRWS, BDROWS2, COLS, 
     $                          ONE, DWORK( XIJ+XRIND-1 ), MB + 1, 
     $                          DWORK( BDUP+2*(KBDUP-1)*(NB+1)**2 +
     $                          COLS*(NB+1) ), NB + 1, ONE, 
     $                          DWORK(GIJ+XRIND-1), MB + 1 )
                        END IF                       
C     
C     Save the result in the matrices F and G
C     
                        CALL DUBEXMA( EXROW, EEXCOL, ROWS, BDROWS2, IX, 
     $                       JX, DWORK(F), LLDFG, NBCBD, MB, NB, 
     $                       DWORK( EXF ),DWORK(FIJ), MB + 1 )
                        CALL DUBEXMA( EXROW, EEXCOL, ROWS, BDROWS2, IX, 
     $                       JX, DWORK(G), LLDFG, NBCBD, MB, NB, 
     $                       DWORK( EXG ), DWORK(GIJ), MB + 1 )
                        IF( BDUPBL.NE.1 ) KBDUP = KBDUP + 1
                     END IF
 80               CONTINUE
C
 75               CONTINUE
C   
C     Update inner loop variable j (column index)
C     
                  IF( .NOT.TRANSAC ) THEN
                     J = J + 1
                  ELSEIF( TRANSAC ) THEN
                     J = J - 1
                  END IF
 70            CONTINUE
            END IF
C     
C     Broadcast the matrix pair (Fij,Gij) in block column j
C     
            IF( K.LT.NROLL .AND. .NOT. TRANSBD ) THEN
               J = JJS
               DO 90 I = IIS, IIE, -1     
                  ROWS = IWORK( IROWS + IIS - I )
                  GSI = IWORK( IGSI + IIS - I ) 
                  EXROW = INT2LG(IWORK( IEXRW + IIS - I )) 
                  COLS = IWORK( ICOLS + IIS - I ) 
                  GSJ = IWORK( IGSJ + IIS - I )
                  EXCOL = INT2LG(IWORK( IEXCL + IIS - I ))
C     
                  IF( ROWS.EQ.0 .OR. COLS.EQ.0 ) GO TO 95
C
                  CALL INFOG2L( GSI, GSJ, DESCFG, NPROW, NPCOL, MYROW, 
     $                          MYCOL, IX, JX, ERSRC, ECSRC )
C     
                  SNODE = MYROW.EQ.ERSRC.AND.MYCOL.EQ.ECSRC
                  SRSRC = ERSRC
                  SCSRC = ECSRC
C
                  IF( SNODE ) THEN
                     CALL DBEXMAT( EXROW, EXCOL, ROWS, COLS, 
     $                    IX, JX, DWORK(F), LLDFG, NBCBD, MB, 
     $                    NB, DWORK( EXF ), DWORK(FIJ), 
     $                    MB + 1 )
                     CALL DBEXMAT( EXROW, EXCOL, ROWS, COLS, 
     $                    IX, JX, DWORK(G), LLDFG, NBCBD, MB, 
     $                    NB, DWORK( EXG ), DWORK(GIJ), 
     $                    MB + 1 ) 
                     XRWS = IWORK( IXRWS + IIS - I )  
                     XRIND = IWORK( IXRIND + IIS - I )
                     IBUFF(1) = XRWS
                     IBUFF(2) = XRIND
                     IF( NPROW.GT.1 ) THEN
                        CALL IGEBS2D( ICTXT, 'Col', ' ', 2, 1, IBUFF, 2)
                        CALL DGEBS2D( ICTXT, 'Col', ' ', XRWS, 2*(NB+1), 
     $                       DWORK(FIJ+XRIND-1), MB + 1 ) 
                     END IF                     
                  ELSEIF( NPROW.GT.1 .AND. MYCOL.EQ.SCSRC ) THEN
                     CALL IGEBR2D( ICTXT, 'Col', ' ', 2, 1, IBUFF, 2, 
     $                    SRSRC, SCSRC )
                     IWORK( IXRWS + IIS - I ) = IBUFF(1)
                     IWORK( IXRIND + IIS - I ) = IBUFF(2)
                     XRWS = IWORK( IXRWS + IIS - I ) 
                     XRIND = IWORK( IXRIND + IIS - I )
                     CALL DGEBR2D( ICTXT, 'Col', ' ', XRWS, 2*(NB+1),
     $                    DWORK(FIJ+XRIND-1), MB + 1, SRSRC, SCSRC )
                  END IF
 95               CONTINUE
C     
C     Update inner loop variable j (column index)
C     
                  IF( .NOT.TRANSAC ) THEN
                     J = J - 1
                  ELSEIF( TRANSAC ) THEN
                     J = J + 1
                  END IF
 90            CONTINUE
C
            ELSEIF( K.LT.NROLL .AND. TRANSBD ) THEN
               J = JJS
               DO 100 I = IIS, IIE, -1     
                  ROWS = IWORK( IROWS + IIS - I )
                  GSI = IWORK( IGSI + IIS - I ) 
                  EXROW = INT2LG(IWORK( IEXRW + IIS - I )) 
                  COLS = IWORK( ICOLS + IIS - I ) 
                  GSJ = IWORK( IGSJ + IIS - I )
                  EXCOL = INT2LG(IWORK( IEXCL + IIS - I ))
C     
                  IF( ROWS.EQ.0 .OR. COLS.EQ.0 ) GO TO 105
C
                  CALL INFOG2L( GSI, GSJ, DESCFG, NPROW, NPCOL, MYROW, 
     $                          MYCOL, IX, JX, ERSRC, ECSRC )
C     
                  SNODE = MYROW.EQ.ERSRC.AND.MYCOL.EQ.ECSRC
                  SRSRC = ERSRC
                  SCSRC = ECSRC
C
                  IF( SNODE ) THEN
                     CALL DBEXMAT( EXROW, EXCOL, ROWS, COLS, 
     $                    IX, JX, DWORK(F), LLDFG, NBCBD, MB, 
     $                    NB, DWORK( EXF ), DWORK(FIJ), 
     $                    MB + 1 )
                     CALL DBEXMAT( EXROW, EXCOL, ROWS, COLS, 
     $                    IX, JX, DWORK(G), LLDFG, NBCBD, MB, 
     $                    NB, DWORK( EXG ), DWORK(GIJ), 
     $                    MB + 1 ) 
                     XRWS = IWORK( IXRWS + IIS - I )  
                     XRIND = IWORK( IXRIND + IIS - I )
                     IBUFF(1) = XRWS
                     IBUFF(2) = XRIND
                     IF( NPROW.GT.1 ) THEN
                        CALL IGEBS2D( ICTXT, 'Col', ' ', 2, 1, IBUFF, 2)
                        CALL DGEBS2D( ICTXT, 'Col', ' ', XRWS, 2*(NB+1), 
     $                       DWORK(FIJ+XRIND-1), MB + 1 ) 
                     END IF                     
                  ELSEIF( NPROW.GT.1 .AND. MYCOL.EQ.SCSRC ) THEN
                     CALL IGEBR2D( ICTXT, 'Col', ' ', 2, 1, IBUFF, 2, 
     $                    SRSRC, SCSRC )
                     IWORK( IXRWS + IIS - I ) = IBUFF(1)
                     IWORK( IXRIND + IIS - I ) = IBUFF(2)
                     XRWS = IWORK( IXRWS + IIS - I ) 
                     XRIND = IWORK( IXRIND + IIS - I )
                     CALL DGEBR2D( ICTXT, 'Col', ' ', XRWS, 2*(NB+1),
     $                    DWORK(FIJ+XRIND-1), MB + 1, SRSRC, SCSRC )
                  END IF
 105              CONTINUE
C     
C     Update inner loop variable j (column index)
C     
                  IF( .NOT.TRANSAC ) THEN
                     J = J + 1
                  ELSEIF( TRANSAC ) THEN
                     J = J - 1
                  END IF
 100           CONTINUE
            END IF
C     
C     Go on with the updates. Now update E with respect to the values of
C     F and G.
C     
            J = JJS
            DO 110 I = IIS, IIE, -1
C     
               ROWS = IWORK( IROWS + IIS - I )
               GSI = IWORK( IGSI + IIS - I ) 
               EXROW = INT2LG(IWORK( IEXRW + IIS - I )) 
               COLS = IWORK( ICOLS + IIS - I ) 
               GSJ = IWORK( IGSJ + IIS - I )
               EXCOL = INT2LG(IWORK( IEXCL + IIS - I ))
C     
               IF( ROWS.EQ.0 .OR. COLS.EQ.0 ) GO TO 115
C
               IF( .NOT.TRANSAC ) THEN
                  IF( MOD( J, 2 ).EQ.0 ) THEN
                     INDXS = 1
                     INDXE = I - 1
                     INDXU = 1
                  ELSE
                     INDXS = I - 1
                     INDXE = 1
                     INDXU = -1
                  END IF
                  KACUP = 1
                  DO 120 INDX = INDXS, INDXE, INDXU
C     
C     Set some constants descibing the involved submatrices
C     
                     IF( INDX.EQ.1 ) THEN
                        IF( IWORK(EXACINF).EQ.0 ) THEN
                           ACROWS2 = MB - IROFFAC
                           GSIND = 1 + IROFFAC
                           CEXROW = .FALSE.
                        ELSEIF( IWORK(EXACINF).EQ.1 ) THEN
                           ACROWS2 = MB -IROFFAC + 1
                           GSIND = 1 + IROFFAC
                           CEXROW = .TRUE.
                        END IF
                     ELSE
                        IF( IWORK(EXACINF+(INDX-1)).EQ.0 ) THEN
                           ACROWS2 = MIN(MB, (M - MB + IROFFAC) - 
     $                          (INDX - 2) * MB)
                           GSIND = (INDX - 1) * MB + 1
                           CEXROW = .FALSE.
                        ELSEIF( IWORK(EXACINF+(INDX-1)).EQ.1 ) THEN
                           ACROWS2 = MB + 1
                           GSIND = (INDX - 1) * MB + 1
                           CEXROW = .TRUE.
                        ELSEIF( IWORK(EXACINF+(INDX-1)).EQ.2 ) THEN
                           ACROWS2 = MIN(MB, (M - MB + IROFFAC) - 
     $                          (INDX - 2) * MB) - 1
                           GSIND = (INDX - 1) * MB + 2
                           CEXROW = .FALSE.
                        ELSEIF( IWORK(EXACINF+(INDX-1)).EQ.3 ) THEN
                           ACROWS2 = MB 
                           GSIND = (INDX - 1) * MB + 2
                           CEXROW = .TRUE.
                        END IF
                     END IF
C     
                     CALL INFOG2L( GSIND, GSI, DESCSA, NPROW, NPCOL,
     $                    MYROW, MYCOL, LIAC, LJAC, RSRC, CSRC )
                     CALL INFOG2L( GSIND, GSJ, DESCSE, NPROW, NPCOL, 
     $                    MYROW, MYCOL, IX, JX, RRSRC, RCSRC )
C
C     Submatrix construction is omitted most of the times for
C     deep pipelining
C     
                     IF( IDEEP.NE.IDEEPS ) GO TO 122
                     IF( MYROW.EQ.RSRC .AND. MYCOL.EQ.CSRC ) THEN
                        IF( (RSRC.NE.RRSRC) .OR. (CSRC.NE.RCSRC) ) THEN
                           CALL DBEXMAT( CEXROW, EXROW, ACROWS2, ROWS, 
     $                          LIAC, LJAC, A((LJAS-1)*LLDA+LIAS), LLDA, 
     $                          NBCAC, MB, MB, DWORK( EXA ), 
     $                          DWORK( ACUP+2*(ACUPBL-1)*(MB+1)**2 ), 
     $                          MB + 1 )
                           CALL DBEXMAT( CEXROW, EXROW, ACROWS2, ROWS, 
     $                          LIAC, LJAC, C((LJCS-1)*LLDC+LICS), LLDC, 
     $                          NBCAC, MB, MB, DWORK( EXC ), 
     $                          DWORK( ACUP+2*(ACUPBL-1)*(MB+1)**2 +
     $                          ROWS*(MB+1) ), MB + 1 )
                           CALL DGESD2D( ICTXT, ACROWS2, 2*ROWS, 
     $                          DWORK( ACUP+2*(ACUPBL-1)*(MB+1)**2 ), 
     $                          MB + 1, RRSRC, RCSRC )
                        ELSE
                           CALL DBEXMAT( CEXROW, EXROW, ACROWS2, ROWS, 
     $                          LIAC, LJAC, A((LJAS-1)*LLDA+LIAS), LLDA, 
     $                          NBCAC, MB, MB, DWORK( EXA ), 
     $                          DWORK( ACUP+2*(KACUP-1)*(MB+1)**2 ), 
     $                          MB + 1 )
                           CALL DBEXMAT( CEXROW, EXROW, ACROWS2, ROWS, 
     $                          LIAC, LJAC, C((LJCS-1)*LLDC+LICS), LLDC, 
     $                          NBCAC, MB, MB, DWORK( EXC ), 
     $                          DWORK( ACUP+2*(KACUP-1)*(MB+1)**2 +
     $                          ROWS*(MB+1) ), MB + 1 )
                        END IF
                     END IF
C
C     Skipped submatrix construction?
C
 122                 CONTINUE
C
                     IF( MYROW.EQ.RRSRC .AND. MYCOL.EQ.RCSRC ) THEN
                        CALL DBEXMAT( CEXROW, EXCOL, ACROWS2, COLS, IX, 
     $                       JX, E((LJES-1)*LLDE+LIES), LLDE, NBCBD, 
     $                       MB, NB, DWORK( EXE ), DWORK(EKJ), 
     $                       MB + 1 )
                        IF( IDEEP.NE.IDEEPS ) GO TO 124
                        IF( (RSRC.NE.RRSRC) .OR. (CSRC.NE.RCSRC) ) THEN
                           CALL DGERV2D( ICTXT, ACROWS2, 2*ROWS,
     $                          DWORK( ACUP+2*(KACUP-1)*(MB+1)**2 ), 
     $                          MB + 1, RSRC, CSRC )
                        END IF
 124                    CONTINUE
C     
C     Perform the update of Ekj
C     
                        XRWS = IWORK( IXRWS + IIS - I )
                        XRIND = IWORK( IXRIND + IIS - I )
                        IF( XRWS.GT.0 .AND. ACROWS2.GT.0 ) THEN 
                           CALL DGEMM( 'N', 'N', ACROWS2, COLS, XRWS, 
     $                          -ONE, DWORK(ACUP+2*(KACUP-1)*(MB+1)**2+
     $                          (XRIND-1)*(MB+1)), MB + 1, 
     $                          DWORK( FIJ+XRIND-1 ), MB + 1, ONE, 
     $                          DWORK(EKJ), MB + 1 )
                           CALL DGEMM( 'N', 'N', ACROWS2, COLS, XRWS, 
     $                          USIGN, DWORK(ACUP+2*(KACUP-1)*(MB+1)**2+
     $                          (XRIND-1)*(MB+1) + ROWS*(MB+1)), MB + 1, 
     $                          DWORK( GIJ+XRIND-1 ), MB + 1, ONE, 
     $                          DWORK(EKJ), MB + 1 )
                        END IF
C     
C     Save the result in the matrix E
C     
                        CALL DUBEXMA( CEXROW, EXCOL, ACROWS2, COLS, IX, 
     $                       JX, E((LJES-1)*LLDE+LIES), LLDE, NBCBD,
     $                       MB, NB, DWORK( EXE ), DWORK(EKJ), 
     $                       MB + 1 ) 
                        IF( ACUPBL.NE.1 ) KACUP = KACUP + 1
                     END IF
 120              CONTINUE
C     
               ELSEIF( TRANSAC ) THEN
                  IF( MOD( J, 2 ).EQ.0 ) THEN
                     INDXS = I + 1
                     INDXE = DBAC
                     INDXU = 1
                  ELSE
                     INDXS = DBAC
                     INDXE = I + 1
                     INDXU = -1
                  END IF
                  KACUP = 1
                  DO 130 INDX = INDXS, INDXE, INDXU
C     
C     Set some constants descibing the involved submatrices
C     
                     IF( IWORK(EXACINF+(INDX-1)).EQ.0 ) THEN
                        ACCOLS2 = MIN(MB, (M - MB + IROFFAC) - 
     $                       (INDX - 2) * MB)
                        GSIND = (INDX - 1) * MB + 1
                        CEXROW = .FALSE.
                     ELSEIF( IWORK(EXACINF+(INDX-1)).EQ.1 ) THEN
                        ACCOLS2 = MB + 1
                        GSIND = (INDX - 1) * MB + 1
                        CEXROW = .TRUE.
                     ELSEIF( IWORK(EXACINF+(INDX-1)).EQ.2 ) THEN
                        ACCOLS2 = MIN(MB, (M - MB + IROFFAC) - 
     $                       (INDX - 2) * MB) - 1
                        GSIND = (INDX - 1) * MB + 2
                        CEXROW = .FALSE.
                     ELSEIF( IWORK(EXACINF+(INDX-1)).EQ.3 ) THEN
                        ACCOLS2 = MB 
                        GSIND = (INDX - 1) * MB + 2
                        CEXROW = .TRUE.
                     END IF
C     
                     CALL INFOG2L( GSI, GSIND, DESCSA, NPROW, NPCOL,
     $                    MYROW, MYCOL, LIAC, LJAC, RSRC, CSRC )
                     CALL INFOG2L( GSIND, GSJ, DESCSE, NPROW, NPCOL, 
     $                    MYROW, MYCOL, IX, JX, RRSRC, RCSRC )
C
C     Submatrix construction is omitted most of the times for
C     deep pipelining
C     
                     IF( IDEEP.NE.IDEEPS ) GO TO 132
                     IF( MYROW.EQ.RSRC .AND. MYCOL.EQ.CSRC ) THEN
                        IF( (RSRC.NE.RRSRC) .OR. (CSRC.NE.RCSRC) ) THEN
                           CALL DBEXMAT( EXROW, CEXROW, ROWS, ACCOLS2,  
     $                          LIAC, LJAC, A((LJAS-1)*LLDA+LIAS), LLDA, 
     $                          NBCAC, MB, MB, DWORK( EXA ), 
     $                          DWORK( ACUP+2*(ACUPBL-1)*(MB+1)**2 ), 
     $                          MB + 1 )
                           CALL DBEXMAT( EXROW, CEXROW, ROWS, ACCOLS2, 
     $                          LIAC, LJAC, C((LJCS-1)*LLDC+LICS), LLDC, 
     $                          NBCAC, MB, MB, DWORK( EXC ), 
     $                          DWORK( ACUP+2*(ACUPBL-1)*(MB+1)**2 +
     $                          ACCOLS2*(MB+1) ), MB + 1 )
                           CALL DGESD2D( ICTXT, ROWS, 2*ACCOLS2, 
     $                          DWORK( ACUP+2*(ACUPBL-1)*(MB+1)**2 ), 
     $                          MB + 1, RRSRC, RCSRC )
                        ELSE
                           CALL DBEXMAT( EXROW, CEXROW, ROWS, ACCOLS2, 
     $                          LIAC, LJAC, A((LJAS-1)*LLDA+LIAS), LLDA, 
     $                          NBCAC, MB, MB, DWORK( EXA ), 
     $                          DWORK( ACUP+2*(KACUP-1)*(MB+1)**2 ), 
     $                          MB + 1 )
                           CALL DBEXMAT( EXROW, CEXROW, ROWS, ACCOLS2, 
     $                          LIAC, LJAC, C((LJCS-1)*LLDC+LICS), LLDC, 
     $                          NBCAC, MB, MB, DWORK( EXC ), 
     $                          DWORK( ACUP+2*(KACUP-1)*(MB+1)**2 +
     $                          ACCOLS2*(MB+1) ), MB + 1 )
                        END IF
                     END IF
C
C     Skipped submatrix construction?
C
 132                 CONTINUE
C
                     IF( MYROW.EQ.RRSRC .AND. MYCOL.EQ.RCSRC ) THEN
                        CALL DBEXMAT( CEXROW, EXCOL, ACCOLS2, COLS, IX, 
     $                       JX, E((LJES-1)*LLDE+LIES), LLDE, NBCBD, 
     $                       MB, NB, DWORK( EXE ), DWORK(EKJ), 
     $                       MB + 1 )
                        IF( IDEEP.NE.IDEEPS ) GO TO 134
                        IF( (RSRC.NE.RRSRC) .OR. (CSRC.NE.RCSRC) ) THEN
                           CALL DGERV2D( ICTXT, ROWS, 2*ACCOLS2,
     $                          DWORK( ACUP+2*(KACUP-1)*(MB+1)**2 ), 
     $                          MB + 1, RSRC, CSRC )
                        END IF
 134                    CONTINUE
C     
C     Perform the update of Ekj
C     
                        XRWS = IWORK( IXRWS + IIS - I )
                        XRIND = IWORK( IXRIND + IIS - I )
                        IF( XRWS.GT.0 .AND. ACCOLS2.GT.0 ) THEN 
                           CALL DGEMM( 'T', 'N', ACCOLS2, COLS, XRWS, 
     $                          -ONE, DWORK(ACUP+2*(KACUP-1)*(MB+1)**2+
     $                          XRIND-1), MB + 1, DWORK(FIJ+XRIND-1),
     $                          MB + 1, ONE, DWORK(EKJ), MB + 1 )
                           CALL DGEMM( 'T', 'N', ACCOLS2, COLS, XRWS, 
     $                          USIGN, DWORK(ACUP+2*(KACUP-1)*(MB+1)**2+
     $                          XRIND-1 + ACCOLS2*(MB+1) ), MB+1, 
     $                          DWORK(GIJ+XRIND-1), MB + 1, ONE, 
     $                          DWORK(EKJ), MB + 1 ) 
                        END IF
C     
C     Save the result in the matrix E
C     
                        CALL DUBEXMA( CEXROW, EXCOL, ACCOLS2, COLS, IX, 
     $                       JX, E((LJES-1)*LLDE+LIES), LLDE, NBCBD, 
     $                       MB, NB, DWORK( EXE ), DWORK(EKJ), 
     $                       MB + 1 ) 
                        IF( ACUPBL.NE.1 ) KACUP = KACUP + 1
                     END IF
 130              CONTINUE
               END IF
C     
 115           CONTINUE
C     
C     Update inner loop variable j (column index)
C     
               IF( (.NOT.TRANSAC).AND.(.NOT.TRANSBD) ) THEN
                  J = J - 1
               ELSEIF( TRANSAC.AND.(.NOT.TRANSBD) ) THEN
                  J = J + 1
               ELSEIF( (.NOT.TRANSAC).AND.TRANSBD ) THEN
                  J = J + 1
               ELSEIF( TRANSAC.AND.TRANSBD ) THEN
                  J = J - 1
               END IF
C     
 110        CONTINUE
C     
C     Prepare for solving for the next block diagonal of E
C     
            J = JJS
            DO 140 I = IIS, IIE, -1
C     
               ROWS = IWORK( IROWS + IIS - I )
               GSI = IWORK( IGSI + IIS - I ) 
               EXROW = INT2LG(IWORK( IEXRW + IIS - I )) 
               COLS = IWORK( ICOLS + IIS - I ) 
               GSJ = IWORK( IGSJ + IIS - I )
               EXCOL = INT2LG(IWORK( IEXCL + IIS - I ))
C     
               IF( ROWS.EQ.0 .OR. COLS.EQ.0 ) GO TO 145
C
               IF((.NOT.TRANSAC.AND.TRANSBD).OR.(TRANSAC.AND.TRANSBD))
     $              THEN
                  IF( J.GT.1 ) THEN
C     
                     IF( J.EQ.2 ) THEN
                        IF( IWORK(EXBDINF).EQ.0 ) THEN
                           ECOLS = NB - IROFFBD
                           GSIND = 1 + IROFFBD
                           CEXCOL = .FALSE.
                        ELSEIF( IWORK(EXBDINF).EQ.1 ) THEN
                           ECOLS = NB - IROFFBD + 1
                           GSIND = 1 + IROFFBD
                           CEXCOL = .TRUE.
                        END IF
                     ELSE
                        IF( IWORK(EXBDINF+(J-2)).EQ.0 ) THEN
                           ECOLS = MIN(NB, (N - NB + IROFFBD) - 
     $                          (J - 3) * NB)
                           GSIND = (J - 2) * NB + 1
                           CEXCOL = .FALSE.
                        ELSEIF( IWORK(EXBDINF+(J-2)).EQ.1 ) THEN
                           ECOLS = NB + 1
                           GSIND = (J - 2) * NB + 1
                           CEXCOL = .TRUE.
                        ELSEIF( IWORK(EXBDINF+(J-2)).EQ.2 ) THEN
                           ECOLS = MIN(NB, (N - NB + IROFFBD) - 
     $                          (J - 3) * NB) - 1
                           GSIND = (J - 2) * NB + 2
                           CEXCOL = .FALSE.
                        ELSEIF( IWORK(EXBDINF+(J-2)).EQ.3 ) THEN
                           ECOLS = NB 
                           GSIND = (J - 2) * NB + 2
                           CEXCOL = .TRUE.
                        END IF
                     END IF
C     
                     CALL INFOG2L( GSI, GSI, DESCSA, NPROW, NPCOL, 
     $                    MYROW, MYCOL, IIA, JJA, ACRSRC, ACCSRC )
                     CALL INFOG2L( GSI, GSIND, DESCSE, NPROW, NPCOL, 
     $                    MYROW, MYCOL, EI, EJ, ERSRC, ECSRC )
                     CALL INFOG2L( GSI, GSIND, DESCFG, NPROW, NPCOL, 
     $                    MYROW, MYCOL, FGI, FGJ, ERSRC, ECSRC )
C
C     Submatrix construction is omitted most of the times for
C     deep pipelining
C
                     IF( IDEEP.NE.IDEEPS ) GO TO 142
C     
                     IF( MYROW.EQ.ACRSRC .AND. MYCOL.EQ.ACCSRC ) THEN
                        CALL DBEXMAT( EXROW, EXROW, ROWS, ROWS, IIA, 
     $                       JJA, A((LJAS-1)*LLDA+LIAS), LLDA, NBCAC, 
     $                       MB, MB, DWORK( EXA ), DWORK( ACII2 ), 
     $                       MB + 1 )
                        CALL DBEXMAT( EXROW, EXROW, ROWS, ROWS, IIA, 
     $                       JJA, C((LJCS-1)*LLDC+LICS), LLDC, NBCAC, 
     $                       MB, MB, DWORK( EXC ),
     $                       DWORK( ACII2+ROWS*(MB+1) ), MB + 1 )
                        IF( ERSRC.NE.ACRSRC .OR. ECSRC.NE.ACCSRC ) THEN
                           CALL DGESD2D( ICTXT, ROWS, 2*ROWS, 
     $                          DWORK( ACII2 ), MB + 1, ERSRC, ECSRC )
                        END IF
                     END IF
C
C     Skipped submatrix construction?
C
 142                 CONTINUE
C
                     IF( MYROW.EQ.ERSRC .AND. MYCOL.EQ.ECSRC ) THEN     
                        CALL DBEXMAT( EXROW, CEXCOL, ROWS, ECOLS, EI, 
     $                       EJ, E((LJES-1)*LLDE+LIES), LLDE, NBCBD, 
     $                       MB, NB, DWORK( EXE ), DWORK( EKJ ), 
     $                       MB + 1 )
                        CALL DBEXMAT( EXROW, CEXCOL, ROWS, ECOLS, FGI, 
     $                       FGJ, DWORK( F ), LLDFG, NBCBD, MB, NB, 
     $                       DWORK( EXF ), DWORK( FIJ ), MB + 1 )
                        CALL DBEXMAT( EXROW, CEXCOL, ROWS, ECOLS, FGI, 
     $                       FGJ, DWORK( G ), LLDFG, NBCBD, MB, NB, 
     $                       DWORK( EXG ), DWORK( GIJ ), MB + 1 )
                        IF( IDEEP.NE.IDEEPS ) GO TO 144
                        IF( ERSRC.NE.ACRSRC .OR. ECSRC.NE.ACCSRC ) THEN
                           CALL DGERV2D( ICTXT, ROWS, 2*ROWS, 
     $                          DWORK( ACII2 ), MB + 1, ACRSRC, 
     $                          ACCSRC )
                        END IF
 144                    CONTINUE
C     
C     Do update
C     
                        XRWS = IWORK( IXRWS + IIS - I )
                        XRIND = IWORK( IXRIND + IIS - I )
                        IF( .NOT. TRANSAC ) THEN
                           CALL DGEMM( TRANAC, 'N', ROWS, ECOLS, XRWS, 
     $                          -ONE, DWORK( ACII2+(XRIND-1)*(MB+1) ), 
     $                          MB + 1, DWORK( FIJ+XRIND-1 ), MB + 1, 
     $                          ONE, DWORK(EKJ), MB + 1 )
                           CALL DGEMM( TRANAC, 'N', ROWS, ECOLS, XRWS, 
     $                          USIGN, DWORK( ACII2+(XRIND-1)*(MB+1)+
     $                          ROWS*(MB+1) ), MB+1, DWORK(GIJ+XRIND-1), 
     $                          MB + 1, ONE, DWORK(EKJ), MB + 1 ) 
                        ELSE
                           CALL DGEMM( TRANAC, 'N', ROWS, ECOLS, XRWS, 
     $                          -ONE, DWORK( ACII2+XRIND-1), MB + 1, 
     $                          DWORK( FIJ+XRIND-1 ), MB + 1, ONE, 
     $                          DWORK(EKJ), MB + 1 )
                           CALL DGEMM( TRANAC, 'N', ROWS, ECOLS, XRWS, 
     $                          USIGN, DWORK( ACII2+XRIND-1 +
     $                          ROWS*(MB+1) ), MB+1, DWORK(GIJ+XRIND-1), 
     $                          MB + 1, ONE, DWORK(EKJ), MB + 1 ) 
                        END IF
                           
C     
C     Save the result in the matrix E
C     
                        CALL DUBEXMA( EXROW, CEXCOL, ROWS, ECOLS, EI, 
     $                       EJ, E((LJES-1)*LLDE+LIES), LLDE, NBCBD, 
     $                       MB, NB, DWORK( EXE ), DWORK(EKJ), 
     $                       MB + 1 ) 
                     END IF
                  END IF
C     
               ELSEIF((TRANSAC.AND..NOT.TRANSBD).OR.(.NOT.TRANSAC.AND.
     $                 .NOT.TRANSBD )) THEN
                  IF( J.LT.DBBD ) THEN
C     
                     IF( IWORK(EXBDINF+J).EQ.0 ) THEN
                        ECOLS = MIN(NB, (N - NB + IROFFBD) - (J-1) * NB)
                        GSIND = J * NB + 1
                        CEXCOL = .FALSE.
                     ELSEIF( IWORK(EXBDINF+J).EQ.1 ) THEN
                        ECOLS = NB + 1
                        GSIND = J * NB + 1
                        CEXCOL = .TRUE.
                     ELSEIF( IWORK(EXBDINF+J).EQ.2 ) THEN
                        ECOLS = MIN(NB, (N - NB + IROFFBD) - 
     $                       (J-1) * NB) - 1
                        GSIND = J * NB + 2
                        CEXCOL = .FALSE.
                     ELSEIF( IWORK(EXBDINF+J).EQ.3 ) THEN
                        ECOLS = NB 
                        GSIND = J * NB + 2
                        CEXCOL = .TRUE.
                     END IF
C     
                     CALL INFOG2L( GSI, GSI, DESCSA, NPROW, NPCOL, 
     $                    MYROW, MYCOL, IIA, JJA, ACRSRC, ACCSRC )
                     CALL INFOG2L( GSI, GSIND, DESCSE, NPROW, NPCOL, 
     $                    MYROW, MYCOL, EI, EJ, ERSRC, ECSRC )
                     CALL INFOG2L( GSI, GSIND, DESCFG, NPROW, NPCOL, 
     $                    MYROW, MYCOL, FGI, FGJ, ERSRC, ECSRC )
C
C     Submatrix construction is omitted most of the times for
C     deep pipelining
C
                     IF( IDEEP.NE.IDEEPS ) GO TO 146
                     IF( MYROW.EQ.ACRSRC .AND. MYCOL.EQ.ACCSRC ) THEN
                        CALL DBEXMAT( EXROW, EXROW, ROWS, ROWS, IIA, 
     $                       JJA, A((LJAS-1)*LLDA+LIAS), LLDA, NBCAC, 
     $                       MB, MB, DWORK( EXA ), DWORK( ACII2 ), 
     $                       MB + 1 )
                        CALL DBEXMAT( EXROW, EXROW, ROWS, ROWS, IIA, 
     $                       JJA, C((LJCS-1)*LLDC+LICS), LLDC, NBCAC, 
     $                       MB, MB, DWORK( EXC ),
     $                       DWORK( ACII2+ROWS*(MB+1) ), MB + 1 )
                        IF( ERSRC.NE.ACRSRC .OR. ECSRC.NE.ACCSRC ) THEN
                           CALL DGESD2D( ICTXT, ROWS, 2*ROWS, 
     $                          DWORK( ACII2 ), MB + 1, ERSRC, 
     $                          ECSRC )
                        END IF
                     END IF
C     
C     Skipped submatrix construction?
C     
 146                 CONTINUE
C     
                     IF( MYROW.EQ.ERSRC .AND. MYCOL.EQ.ECSRC ) THEN
                        CALL DBEXMAT( EXROW, CEXCOL, ROWS, ECOLS, EI, 
     $                       EJ, E((LJES-1)*LLDE+LIES), LLDE, NBCBD, 
     $                       MB, NB, DWORK( EXE ), DWORK( EKJ ), 
     $                       MB + 1 )
                        CALL DBEXMAT( EXROW, CEXCOL, ROWS, ECOLS, FGI, 
     $                       FGJ, DWORK( F ), LLDFG, NBCBD, MB, NB, 
     $                       DWORK( EXF ), DWORK( FIJ ), MB + 1 )
                        CALL DBEXMAT( EXROW, CEXCOL, ROWS, ECOLS, FGI, 
     $                       FGJ, DWORK( G ), LLDFG, NBCBD, MB, NB, 
     $                       DWORK( EXG ), DWORK( GIJ ), MB + 1 )
                        IF( IDEEP.NE.IDEEPS ) GO TO 148
                        IF( ERSRC.NE.ACRSRC .OR. ECSRC.NE.ACCSRC ) THEN
                           CALL DGERV2D( ICTXT, ROWS, 2*ROWS, 
     $                          DWORK( ACII2 ), MB + 1, ACRSRC, 
     $                          ACCSRC )
                        END IF
 148                    CONTINUE
C     
C     Do update
C     
                        XRWS = IWORK( IXRWS + IIS - I )
                        XRIND = IWORK( IXRIND + IIS - I )
                        IF( TRANSAC ) THEN
                           CALL DGEMM( TRANAC, 'N', ROWS, ECOLS, XRWS, 
     $                          -ONE, DWORK( ACII2+XRIND-1 ), MB + 1, 
     $                          DWORK( FIJ+XRIND-1 ), MB + 1, ONE, 
     $                          DWORK(EKJ), MB + 1 )
                           CALL DGEMM( TRANAC, 'N', ROWS, ECOLS, XRWS, 
     $                          USIGN, DWORK(ACII2+XRIND-1+ROWS*(MB+1)), 
     $                          MB + 1, DWORK( GIJ+XRIND-1 ), MB + 1, 
     $                          ONE, DWORK(EKJ), MB + 1 )
                        ELSE
                           CALL DGEMM( TRANAC, 'N', ROWS, ECOLS, XRWS, 
     $                          -ONE, DWORK( ACII2+(XRIND-1)*(MB+1) ), 
     $                          MB + 1, DWORK( FIJ+XRIND-1 ), MB + 1, 
     $                          ONE, DWORK(EKJ), MB + 1 )
                           CALL DGEMM( TRANAC, 'N', ROWS, ECOLS, XRWS, 
     $                          USIGN, DWORK(ACII2+(XRIND-1)*(MB+1)+
     $                          ROWS*(MB+1)), MB + 1, 
     $                          DWORK( GIJ+XRIND-1 ), MB + 1, ONE, 
     $                          DWORK(EKJ), MB + 1 )
                        END IF
C     
C     Save the result in the matrix E
C     
                        CALL DUBEXMA( EXROW, CEXCOL, ROWS, ECOLS, EI, 
     $                       EJ, E((LJES-1)*LLDE+LIES), LLDE, NBCBD, 
     $                       MB, NB, DWORK( EXE ), DWORK(EKJ), 
     $                       MB + 1 ) 
                     END IF
                  END IF             
               END IF
C     
 145           CONTINUE
C     
C     Update inner loop variable j (column index)
C     
               IF( (.NOT.TRANSAC).AND.(.NOT.TRANSBD) ) THEN
                  J = J - 1
               ELSEIF( TRANSAC.AND.(.NOT.TRANSBD) ) THEN
                  J = J + 1
               ELSEIF( (.NOT.TRANSAC).AND.TRANSBD ) THEN
                  J = J + 1
               ELSEIF( TRANSAC.AND.TRANSBD ) THEN
                  J = J - 1
               END IF
C     
 140        CONTINUE
 24         CONTINUE
 20      CONTINUE
C     
C     Update outer loop variables
C     
         IF( (.NOT.TRANSAC).AND.(.NOT.TRANSBD) ) THEN
            JS = MIN( JS + 1, DBBD )
            IEND = MAX( 1, IEND - 1 )
         ELSEIF( TRANSAC.AND.(.NOT.TRANSBD) ) THEN
            IF( K .GE. DBAC ) JS = JS + 1
            IS = MIN( IS + 1, DBAC )
         ELSEIF( (.NOT.TRANSAC).AND.TRANSBD ) THEN
            JS = MAX( JS - 1, 1 )
            IEND = MAX( IEND - 1, 1 )
         ELSEIF( TRANSAC.AND.TRANSBD ) THEN
            IF( K .GE. DBAC ) JS = JS - 1
            IS = MIN( IS + 1, DBAC )
         END IF
C     
 10   CONTINUE
C     
C     Before we go on we must back redistributed the elements in E
C     
      IF( ACEXT.OR.BDEXT) THEN
         CALL PDBCKRD( M, N, E, IE, JE, DESCE, IWORK( EXACINF ), 
     $                 IWORK( EXBDINF ), DWORK( EXE ), EXMEME,
     $                 DWORK( SND ), MAX( MB, NB ), INFO )
      END IF
C
C     Restore Data in Original Position. Now we don't need to
C     shift the extended elements or intermediate matrices, we shift 
C     just the original matrices. 
C
      IF( SHIFT .AND. NPROCS.GT.1 .AND. 
     $     ( SHALT1.OR.N.EQ.NB.OR.M.EQ.MB ) ) THEN
C     
C     Shift (A,C) East if NPCOL > 1 and (A,C)**T West if NPCOL > 1
C     Shift (B,D) North if NPROW > 1 and (B,D)**T South if NPROW > 1
C     
         IF ( NPCOL.GT.1.AND.M.NE.MB ) THEN
            IF( (.NOT.TRANSAC) .AND. (.NOT.TRANSBD) ) THEN
               CALL DGESD2D( ICTXT, ACROWS, ACCOLS, 
     $              A((LJAS-1)*LLDA+LIAS), LLDA, MYROW, 
     $              EAST )
               DESCSA(CSRC_) = MOD(DESCSA(CSRC_) + 1, NPCOL)
               CALL DGERV2D( ICTXT, ACROWS, ACCOLS, 
     $              A((LJAS-1)*LLDA+LIAS), LLDA, MYROW,
     $              WEST )
               CALL DGESD2D( ICTXT, ACROWS, ACCOLS, 
     $              C((LJCS-1)*LLDC+LICS), LLDC, MYROW, 
     $              EAST )
               DESCSC(CSRC_) = MOD(DESCSC(CSRC_) + 1, NPCOL)
               CALL DGERV2D( ICTXT, ACROWS, ACCOLS, 
     $              C((LJCS-1)*LLDC+LICS), LLDC, MYROW,
     $              WEST )
            ELSEIF( TRANSAC .AND. TRANSBD ) THEN
               CALL DGESD2D( ICTXT, ACROWS, ACCOLS, 
     $              A((LJAS-1)*LLDA+LIAS), LLDA, MYROW, 
     $              WEST )   
               DESCSA(CSRC_) = MOD(DESCSA(CSRC_) + NPCOL - 1, 
     $              NPCOL )
               CALL DGERV2D( ICTXT, ACROWS, ACCOLS,
     $              A((LJAS-1)*LLDA+LIAS), LLDA, MYROW, 
     $              EAST )
               CALL DGESD2D( ICTXT, ACROWS, ACCOLS, 
     $              C((LJCS-1)*LLDC+LICS), LLDC, MYROW, 
     $              WEST )   
               DESCSC(CSRC_) = MOD(DESCSC(CSRC_) + NPCOL - 1, 
     $              NPCOL )
               CALL DGERV2D( ICTXT, ACROWS, ACCOLS,
     $              C((LJCS-1)*LLDC+LICS), LLDC, MYROW, 
     $              EAST )
            END IF
         END IF
         IF ( NPROW.GT.1.AND.N.NE.NB ) THEN
            IF( (.NOT.TRANSAC) .AND. (.NOT.TRANSBD) ) THEN
               CALL DGESD2D( ICTXT, BDROWS, BDCOLS,
     $              B((LJBS-1)*LLDB+LIBS), LLDB, NORTH,
     $              MYCOL )
               DESCSB(RSRC_) = MOD(DESCSB(RSRC_) + NPROW - 1,
     $              NPROW )
               CALL DGERV2D( ICTXT, BDROWS, BDCOLS,
     $              B((LJBS-1)*LLDB+LIBS), LLDB, SOUTH,
     $              MYCOL )
               CALL DGESD2D( ICTXT, BDROWS, BDCOLS,
     $              D((LJDS-1)*LLDD+LIDS), LLDD, NORTH,
     $              MYCOL )
               DESCSD(RSRC_) = MOD(DESCSD(RSRC_) + NPROW - 1,
     $              NPROW )
               CALL DGERV2D( ICTXT, BDROWS, BDCOLS,
     $              D((LJDS-1)*LLDD+LIDS), LLDD, SOUTH,
     $              MYCOL )
            ELSEIF( TRANSAC .AND. TRANSBD ) THEN
               CALL DGESD2D( ICTXT, BDROWS, BDCOLS,
     $              B((LJBS-1)*LLDB+LIBS), LLDB, SOUTH,
     $              MYCOL )
               DESCSB(RSRC_) = MOD(DESCSB(RSRC_) + 1, NPROW)
               CALL DGERV2D( ICTXT, BDROWS, BDCOLS,
     $              B((LJBS-1)*LLDB+LIBS), LLDB, NORTH,
     $              MYCOL )
               CALL DGESD2D( ICTXT, BDROWS, BDCOLS,
     $              D((LJDS-1)*LLDD+LIDS), LLDD, SOUTH,
     $              MYCOL )
               DESCSD(RSRC_) = MOD(DESCSD(RSRC_) + 1, NPROW )
               CALL DGERV2D( ICTXT, BDROWS, BDCOLS,
     $              D((LJDS-1)*LLDD+LIDS), LLDD, NORTH,
     $              MYCOL )
            END IF
         END IF
      ELSEIF( SHIFT .AND. NPROCS.GT.1 .AND. SHALT2 ) THEN 
C     
C     Shift (A,C) SouthEast and E South if NPROW > 1 
C     or Shift (A,C)**T NorthWest and E North if NPROW > 1 
C     
         IF( (.NOT.TRANSAC) .AND. (.NOT.TRANSBD) ) THEN
            CALL DGESD2D( ICTXT, ACROWS, ACCOLS,
     $           A((LJAS-1)*LLDA+LIAS), LLDA, SOUTH, 
     $           EAST )
            DESCSA(RSRC_) = MOD(DESCSA(RSRC_) + 1, NPROW)
            DESCSA(CSRC_) = MOD(DESCSA(CSRC_) + 1, NPCOL)
            CALL DGERV2D( ICTXT, ACROWS, ACCOLS, 
     $           A((LJAS-1)*LLDA+LIAS), LLDA, NORTH, 
     $           WEST )
            CALL DGESD2D( ICTXT, ACROWS, ACCOLS,
     $           C((LJCS-1)*LLDC+LICS), LLDC, SOUTH, 
     $           EAST )
            DESCSC(RSRC_) = MOD(DESCSC(RSRC_) + 1, NPROW )
            DESCSC(CSRC_) = MOD(DESCSC(CSRC_) + 1, NPCOL )
            CALL DGERV2D( ICTXT, ACROWS, ACCOLS, 
     $           C((LJCS-1)*LLDC+LICS), LLDC, NORTH, 
     $           WEST )
         ELSEIF( TRANSAC .AND. TRANSBD ) THEN
            CALL DGESD2D( ICTXT, ACROWS, ACCOLS,
     $           A((LJAS-1)*LLDA+LIAS), LLDA, NORTH, 
     $           WEST )
            DESCSA(RSRC_) = MOD(DESCSA(RSRC_) + NPROW - 1, NPROW)
            DESCSA(CSRC_) = MOD(DESCSA(CSRC_) + NPCOL - 1, NPCOL)
            CALL DGERV2D( ICTXT, ACROWS, ACCOLS,
     $           A((LJAS-1)*LLDA+LIAS), LLDA, SOUTH, 
     $           EAST )
            CALL DGESD2D( ICTXT, ACROWS, ACCOLS,
     $           C((LJCS-1)*LLDC+LICS), LLDC, NORTH, 
     $           WEST )
            DESCSC(RSRC_) = MOD(DESCSC(RSRC_) + NPROW - 1, NPROW)
            DESCSC(CSRC_) = MOD(DESCSC(CSRC_) + NPCOL - 1, NPCOL)
            CALL DGERV2D( ICTXT, ACROWS, ACCOLS,
     $           C((LJCS-1)*LLDC+LICS), LLDC, SOUTH, 
     $           EAST )
         END IF
C     
C     Shift E
C     
         IF ( NPROW.GT.1 ) THEN
            IF( (.NOT.TRANSAC) .AND. (.NOT.TRANSBD) ) THEN
               CALL DGESD2D( ICTXT, ACROWS, BDCOLS, 
     $              E((LJES-1)*LLDE+LIES), LLDE, SOUTH,
     $              MYCOL )
               DESCSE(RSRC_) = MOD(DESCSE(RSRC_) + 1, NPROW)
               CALL DGERV2D( ICTXT, ACROWS, BDCOLS,
     $              E((LJES-1)*LLDE+LIES), LLDE, NORTH,
     $              MYCOL )
            ELSEIF( TRANSAC .AND. TRANSBD ) THEN
               CALL DGESD2D( ICTXT, ACROWS, BDCOLS, 
     $              E((LJES-1)*LLDE+LIES), LLDE, NORTH,
     $              MYCOL )
               DESCSE(RSRC_) = MOD(DESCSE(RSRC_) + NPROW-1, NPROW)
               CALL DGERV2D( ICTXT, ACROWS, BDCOLS,
     $              E((LJES-1)*LLDE+LIES), LLDE, SOUTH,
     $              MYCOL )
            END IF
         END IF
      ELSEIF ( SHIFT .AND. NPROCS.GT.1 .AND. SHALT3 ) THEN 
C     
C     Shift (B,D) NorthWest and E West if NPCOL > 1 
C     or shift (B,D)**T SouthEast and E East if NPCOL > 1 
C     
         IF( (.NOT.TRANSAC) .AND. (.NOT.TRANSBD) ) THEN
            CALL DGESD2D( ICTXT, BDROWS, BDCOLS,
     $           B((LJBS-1)*LLDB+LIBS), LLDB, NORTH, 
     $           WEST )
            DESCSB(RSRC_) = MOD(DESCSB(RSRC_) + NPROW - 1, NPROW )
            DESCSB(CSRC_) = MOD(DESCSB(CSRC_) + NPCOL - 1, NPCOL )
            CALL DGERV2D( ICTXT, BDROWS, BDCOLS,
     $           B((LJBS-1)*LLDB+LIBS), LLDB, SOUTH, 
     $           EAST )
            CALL DGESD2D( ICTXT, BDROWS, BDCOLS,
     $           D((LJDS-1)*LLDD+LIDS), LLDD, NORTH, 
     $           WEST )
            DESCSD(RSRC_) = MOD(DESCSD(RSRC_) + NPROW - 1, NPROW )
            DESCSD(CSRC_) = MOD(DESCSD(CSRC_) + NPCOL - 1, NPCOL )
            CALL DGERV2D( ICTXT, BDROWS, BDCOLS,
     $           D((LJDS-1)*LLDD+LIDS), LLDD, SOUTH, 
     $           EAST )
         ELSEIF( TRANSAC .AND. TRANSBD ) THEN
            CALL DGESD2D( ICTXT, BDROWS, BDCOLS,
     $           B((LJBS-1)*LLDB+LIBS), LLDB, SOUTH, 
     $           EAST )
            DESCSB(RSRC_) = MOD( DESCSB(RSRC_) + 1, NPROW )
            DESCSB(CSRC_) = MOD( DESCSB(CSRC_) + 1, NPCOL )
            CALL DGERV2D( ICTXT, BDROWS, BDCOLS, 
     $           B((LJBS-1)*LLDB+LIBS), LLDB, NORTH, 
     $           WEST )
            CALL DGESD2D( ICTXT, BDROWS, BDCOLS,
     $           D((LJDS-1)*LLDD+LIDS), LLDD, SOUTH, 
     $           EAST )
            DESCSD(RSRC_) = MOD( DESCSD(RSRC_) + 1, NPROW )
            DESCSD(CSRC_) = MOD( DESCSD(CSRC_) + 1, NPCOL )
            CALL DGERV2D( ICTXT, BDROWS, BDCOLS, 
     $           D((LJDS-1)*LLDD+LIDS), LLDD, NORTH, 
     $           WEST )
         END IF
C     
C     Shift E
C     
         IF ( NPCOL.GT.1 ) THEN
            IF( (.NOT.TRANSAC) .AND. (.NOT.TRANSBD) ) THEN
               CALL DGESD2D( ICTXT, ACROWS, BDCOLS,
     $              E((LJES-1)*LLDE+LIES), LLDE, MYROW,
     $              WEST )
               DESCSE(CSRC_) = MOD(DESCSE(CSRC_) + NPCOL-1,NPCOL)
               CALL DGERV2D( ICTXT, ACROWS, BDCOLS, 
     $              E((LJES-1)*LLDE+LIES), LLDE, MYROW, 
     $              EAST )
            ELSEIF( TRANSAC .AND. TRANSBD ) THEN
               CALL DGESD2D( ICTXT, ACROWS, BDCOLS,
     $              E((LJES-1)*LLDE+LIES), LLDE, MYROW,
     $              EAST )
               DESCSE(CSRC_) = MOD(DESCSE(CSRC_) + 1, NPCOL)
               CALL DGERV2D( ICTXT, ACROWS, BDCOLS,
     $              E((LJES-1)*LLDE+LIES), LLDE, MYROW, 
     $              WEST )
            END IF
         END IF
      END IF
C     
C     Global max on INFO
C     
      IF( NPROCS.GT.1 ) 
     $     CALL IGAMX2D( ICTXT, 'All', ' ', 1, 1, INFO, 1, -1, -1, -1, 
     $                   -1, -1 )
C     
C     If INFO = 2, global min on SCALE
C     
      IF( INFO.EQ.2 ) THEN
         IF( NPROCS.GT.1 ) THEN
            CALL DGAMN2D( ICTXT, 'All', ' ', 1, 1, SCALE, 1, -1,- 1, -1,
     $                    -1, -1 )
         END IF
         RETURN
      END IF
C
      END
C     
C     End of PTRGSYLD
C     
C *** Last line of PTRGSYLD ***
