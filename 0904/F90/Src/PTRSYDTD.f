CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C         This file is a part of the SCASY software library.           C
C         See http://www.cs.umu.se/parallel/research/scasy             C
C         Contributors: Robert Granat and Bo Kågström.                 C
C         SCASY Version 0.10,    March 31, 2006.                       C
C         Copyright 2006, Umeå University, Sweden.                     C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE PTRSYDTD( TRANSA, TRANSB, ISGN, COMM, M, N, A, IA,
     $                     JA, DESCA, B, IB, JB, DESCB, C, IC, JC, 
     $                     DESCC, MB2, DWORK, LDWORK, IWORK, 
     $                     LIWORK, NOEXSY, SCALE, INFO)
C
C  -- ScaLAPACK-style routine (preliminary version ) --
C     HPC2N and Department of Computing Science,
C     University of Umeå, Sweden.
C     Written by Robert Granat, (granat@cs.umu.se)
C     October 18, 2006.
C
      IMPLICIT NONE
C
C     .. Scalar Arguments ..
      CHARACTER*1        TRANSA, TRANSB, COMM
      INTEGER            ISGN, IA, IB, IC, INFO, JA, JB, JC, M, N, 
     $                   LDWORK, LIWORK, NOEXSY, MB2
      DOUBLE PRECISION   SCALE
C     ..
C     .. Array Arguments ..
      INTEGER            DESCA( * ), DESCB( * ), DESCC( * ), IWORK( * )
      DOUBLE PRECISION   A( * ), B( * ), C( * ), DWORK( * )
C     ..
C
C  Purpose
C  =======
C
C  This subroutine solves the real (quasi-)triangular discrete time 
C  Sylvester equation
C
C     op( sub( A ) ) * X * op( sub( B ) ) +/- X =  sub( C ),
C
C  where sub( A ) = A(IA:IA+M-1,JA:JA+M-1) is an M-by-M distributed
C  matrix and sub( B ) = B(IB:IB+N-1,JB:JB+N-1) is an N-by-N
C  distributed matrix. sub( C ) =  C(IC:IC+M-1,JC:JC+N-1) is an M-by-N
C  distributed matrix and will be overwritten by the solution sub( X ).
C
C  The notation op(_) means the transpose or non-transpose of a matrix.
C
C  This routine should *not* be called directly, but through PGESYDTD.
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
C  Arguments
C  =========
C
C  Mode parameters
C
C  TRANSA    (global input) CHARACTER*1
C            If TRANSA = 'N' then op(A) = A
C            If TRANSA = 'T' then op(A) = A**T
C
C  TRANSB    (global input) CHARACTER*1
C            If TRANSB = 'N' then op(B) = B
C            If TRANSB = 'T' then op(B) = B**T
C
C  ISGN      (global input) INTEGER
C            If ISGN =  1, we solve the equation 
C              op(A) * X + X * op(B) = C
C            If ISGN = -1, we solve the equation 
C              op(A) * X - X * op(B) = C
C
C  Input/Output parameters
C
C  COMM      (global input/output) CHARACTER*1
C            This subroutine uses two different communications schemes in
C            solving the reduced triangular problem. 
C              If COMM = 'S', the "shifts" scheme is used.
C              If COMM = 'D', the "communicate on demand" scheme is used.
C            The choice COMM = 'S' is only valid for TRANSA = TRANSB = 'N' 
C            or TRANSA = TRANSB = 'T'. The scheme used will be output.
C            See the references for details.
C
C  M         (global input) INTEGER
C            The number of rows and columns of the global distributed 
C            matrix A. This is also the number of rows of the
C            global distributed matrix C. M >= 0.
C
C  N         (global input) INTEGER
C            The number of rows and columns of the global distributed
C            matrix B. This is also the number of columns of the
C            global distributed matrix C. N >= 0.
C
C  A         (local input) DOUBLE PRECISION array 
C            Array of dimension (LLD_A,LOCc(M)). Contains the local
C            pieces of the global distributed matrix A in real
C            Schur form.
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
C            pieces of the global distributed matrix B in real
C            Schur form.
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
C            Array of dimension (LLD_C,LOCc(N)). 
C            On entry C contains the local pieces of the global 
C            distributed matrix C . On exit, it contains the local
C            pieces of the global distributed solution X.
C
C  IC        (global input) INTEGER
C            Row start index for sub(C), i.e., the submatrix to operate 
C            on. MOD(IC,MB_A) = MOD(IA,MB_A) must hold. 
C
C  JC        (global input) INTEGER
C            Column start index for sub(C), i.e., the submatrix to 
C            operate on. MOD(JC,MB_B) = MOD(JB,MB_B) must hold.
C
C  DESCC     (global and local input) INTEGER array of dimension DLEN_.
C            The array descriptor for the global distributed matrix C.
C
C  MB2       (global input) INTEGER 
C            Internal blocking factor for pipelining of subsolutions
C            for updates of the matrix C.
C            1 < = MB2 <= DESCC( MB_ ) must hold.
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
C            LIWORK >= DBA + DBB + 8 * MIN( P_r, P_c ), 
C            where DBA = ICEIL(LOCr(IA+IROFFA),MB_A) and 
C            DBB = ICEIL(LOCr(IB+IROFFB),MB_B).
C
C            If LIWORK = -1, LIWORK is global input and a workspace 
C            query is assumed. The routine will then calculate the 
C            optimal workspace needed, store it in IWORK(1) and return
C            immediately. No error will then be signaled by PXERBLA.
C
C  Output information
C            
C  NOEXSY    (local output) INTEGER
C            When solving the triangular problem in PTRSYDTD it is possible
C            that we have to extend some subsystems to not lose any data
C            from some 2x2 block of conjugate pairs of eigenvalues. NOEXSY 
C            helps us to keep track of the number of such extensions.
C
C  SCALE     (global output) DOUBLE PRECISION
C            A scale factor, usually 1.0. See INFO for details.
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
C             If INFO = 2, A and B have common or very close eigenvalues; 
C             perturbed values were used to solve the equations
C             (but A and B are unchanged).
C             If INFO = 3,  the problem is badly scaled - C is scaled
C             a factor SCALE to guarantee an overflow free solution.
C
C  Method
C  ======
C  This subroutine implements a parallel wave-front algorithm for
C  solving the triangular discrete-time Sylvester equation. See
C  the references for details.
C
C  Additional requirements
C  =======================
C
C  A and B must be distributed using the same blocking factor in 
C  each direction, i.e., MB_A=NB_A, MB_B=NB_B. Moreover, for C the 
C  blocksize in the row direction must agree with A's, i.e. MB_C=MB_A
C  must hold, and the blocksize in the column direction must agree with
C  B's, i.e. NB_C=NB_B must hold. 
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
C  Wave-front algorithm, real Schur form, Sylvester equation,
C  explicit blocking, GEMM-updates, matrix shifts, on demand 
C
C  =====================================================================
C
C     .. Parameters ..
      DOUBLE PRECISION MONE, ONE, ZERO
      PARAMETER        ( MONE = -1.0D+0, ONE = 1.0D+0, ZERO = 0.0D+0 )
      INTEGER          BLOCK_CYCLIC_2D, CSRC_, CTXT_, DLEN_, DT_,
     $                 LLD_, MB_, M_, NB_, N_, RSRC_
      PARAMETER        ( BLOCK_CYCLIC_2D = 1, DLEN_ = 9, DT_ = 1,
     $                   CTXT_ = 2, M_ = 3, N_ = 4, MB_ = 5, NB_ = 6,
     $                   RSRC_ = 7, CSRC_ = 8, LLD_ = 9 )
C     ..
C     .. Local Scalars ..
      INTEGER  MB, NB, DBA, DBB, NROLL, IS, JS, IE, K,
     $         MYCOL, MYROW, NPCOL, NPROW, J, NPROC, I, IDUM,
     $         AROWS, ACOLS, BROWS, BCOLS, ROWS, COLS, LINFO,
     $         IX, JX, RSRC, CSRC, LIA, LJA, LIB, LJB, INDX, ICTXT,
     $         LLDA, LLDB, LLDC, SRSRC, SCSRC, WRK, EXMEME, EXE,
     $         MWORKneeded, XIJ, AII, BJJ, CIJ, EIJ, E, RRSRC, RCSRC, 
     $         ARSRC, ACSRC, BRSRC, BCSRC, CRSRC, CCSRC, EXA, EXB, EXC,
     $         EXMEMA, EXMEMB, EXMEMC, IWRK, EXAINF, EXBINF, SND,
     $         EXBUFF, LEXBUFF, GI, GJ, LBI, LBJ, NBCA, NBCB, POS, JEND,
     $         AROWS2, BCOLS2, GINDX, ACOLS2, BROWS2, D, IPW, LLDE,
     $         EROWS, ECOLS, SIZE_E, SIZE_IPW, IIC, JJC, CCOLS, IIA, 
     $         JJA, EI, EJ, JSTART, Da, Db, MWORK, LMATR, NORTH, SOUTH, 
     $         WEST, EAST, NPROCS, IROFFA, IROFFB, LIC, LJC, ASI, ASJ,
     $         BSI, BSJ, CSI, CSJ, LIAS, LJAS, LIBS, LJBS, LICS, LJCS,
     $         GSI, GSJ, GSIND, INDXS, INDXE, INDXU, MINELEM, PHASE,
     $         PHASES, MNPDIM, IROWS, ICOLS, IGSI, IGSJ, IEXRW, IEXCL,
     $         JJS, IIS, IIE, NIDEEP, IDEEPS, IDEEPE, IDEEPU, IDEEP, 
     $         XRWS, XRIND, IXRWS, IXRIND, BUPBL, AUPBL, AUP, 
     $         KAUP, CKJ, BUP, KBUP, SIZE_IPW2, AII2, IPW2, KKK
      DOUBLE PRECISION SCALOC, LSCALC
      LOGICAL  LQUERY, SNODE, TRANA, TRANB, EXROW, EXCOL, CEXROW, 
     $         CEXCOL, EEXCOL, EEXROW, AEXT, BEXT, SHIFT, SHALT1,
     $         SHALT2, SHALT3
C     ..
C     .. Local Arrays ..
      INTEGER IDUM1(1), IDUM2(1), DESCSC( DLEN_ ), DESCSA( DLEN_ ),
     $        DESCSB( DLEN_ ), DESCE( DLEN_ ), IBUFF( 2 )
C     ..
C     .. External Subroutines ..
      EXTERNAL           BLACS_GRIDINFO, CHK1MAT, SB04PY, 
     $                   PXERBLA, INFOG2L, DGEMM, DLACPY,
     $                   PCHK2MAT, DSCAL, DGAMN2D, DGESD2D,
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
      ICTXT = DESCC( CTXT_ )
      CALL BLACS_GRIDINFO( ICTXT, NPROW, NPCOL, MYROW, MYCOL )
C
C     Test the input parameters
C
      LINFO = 0
      INFO = 0
      IF( NPROW.EQ.-1 ) THEN
         INFO = -807
      ELSE
         CALL CHK1MAT( M, 4, M, 4, IA, JA, DESCA, 10, INFO )
         CALL CHK1MAT( N, 5, N, 5, IB, JB, DESCB, 14, INFO )
         CALL CHK1MAT( M, 4, N, 5, IC, JC, DESCC, 18, INFO )
         CALL PCHK2MAT( M, 4, M, 4, IA, JA, DESCA, 10, N, 5, N, 5, IB,
     $                 JB, DESCB, 14, 0, IDUM1, IDUM2, INFO )
         CALL PCHK2MAT( M, 4, M, 4, IA, JA, DESCA, 10, M, 4, N, 5, IC, 
     $                  JC, DESCC, 18, 0, IDUM1, IDUM2, INFO )
      END IF
C
      IF( INFO.EQ.0 ) THEN
         IF( .NOT.( LSAME( TRANSA, 'T' ) .OR. 
     $        LSAME( TRANSA, 'N' ) ) ) THEN
            INFO = -1
         END IF
      END IF
      IF( INFO.EQ.0 ) THEN
         IF( .NOT.( LSAME( TRANSB, 'T' ) .OR. 
     $        LSAME( TRANSB, 'N' ) ) ) THEN
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
         IF( .NOT. ( (LSAME( TRANSA, 'N' ).AND.LSAME( TRANSB, 'N' ) )
     $        .OR.( (LSAME( TRANSA, 'T' ).AND.LSAME( TRANSB, 'T' )))))
     $        COMM = 'D'
      END IF
C
C     Check that the last row and column of the matrices sub(A), sub(B) 
C     and sub(C) can be mapped onto the last process row or column. 
C     If not, change to on-demand communication.
C     
      MB = DESCA( MB_ )
      NB = DESCB( MB_ ) 
      IROFFA = MOD( IA - 1, MB )
      IROFFB = MOD( IB - 1, NB ) 
      DBA = ICEIL( M + IROFFA, DESCA(MB_) )
      DBB = ICEIL( N + IROFFB, DESCB(MB_) )
C
      IF( INFO.EQ.0 .AND. LSAME(COMM,'S')) THEN 
         D = MAX( NPROW, NPCOL )
         IF( MIN( M + IROFFA, N + IROFFB ) .GT. (D-1)**2 ) THEN
            IF( .NOT.(MOD(DBA,D).EQ.0 .OR. MOD(DBB,D).EQ.0) ) THEN
               COMM = 'D'
            END IF
         ELSE
            COMM = 'D'
         END IF
      END IF            
C
C     Set some initial values
C
      NPROCS = NPROW*NPCOL 
      MNPDIM = MIN( NPROW, NPCOL )
      SCALOC = ONE
      SCALE = ONE
      TRANA = LSAME( TRANSA, 'T' )
      TRANB = LSAME( TRANSB, 'T' )
      SHIFT = LSAME( COMM, 'S' )
      AEXT = .FALSE.
      BEXT = .FALSE.
      NOEXSY = 0
      LLDA = DESCA(LLD_)
      LLDB = DESCB(LLD_)
      LLDC = DESCC(LLD_)
C
C     Compute the number of blocks to keep from A and B during the 
C     updates when deep pipelining is turned on
C
      IF( MB2.NE.MB ) THEN
         AUPBL = MAX( 1, ICEIL( DBA-1, MNPDIM ) ) + 1
         BUPBL = MAX( 1, ICEIL( DBB-1, MNPDIM ) ) + 1
      ELSE
         AUPBL = 1
         BUPBL = 1
      END IF     
C
C     Compute the number of rows and columns held by each process and
C     the number of block columns held by each process for the smallest
C     submatrices of A and B including sub(A) and sub(B) conforming with 
C     ScaLAPACK conventions 
C
      CALL INFOG2L( IA, JA, DESCA, NPROW, NPCOL, MYROW, MYCOL, LIA, LJA,
     $              ARSRC, ACSRC )
      CALL INFOG2L( IB, JB, DESCB, NPROW, NPCOL, MYROW, MYCOL, LIB, LJB,
     $              BRSRC, BCSRC )
      CALL INFOG2L( IC, JC, DESCC, NPROW, NPCOL, MYROW, MYCOL, LIC, LJC,
     $              CRSRC, CCSRC )
      AROWS = NUMROC( M + IROFFA, MB, MYROW, ARSRC, NPROW )
      ACOLS = NUMROC( M + IROFFA, MB, MYCOL, ACSRC, NPCOL )
      BROWS = NUMROC( N + IROFFB, NB, MYROW, BRSRC, NPROW )
      BCOLS = NUMROC( N + IROFFB, NB, MYCOL, BCSRC, NPCOL )
      NBCA = ICEIL( ACOLS, MB )
      NBCB = ICEIL( BCOLS, NB )
C
C     Compute the needed memory for holding data for the extended
C     subsystems. Even if all blocks cannot be extended, we add memory
C     for all blocks to maintain indicies simple.
C
      EXMEMA = ICEIL(AROWS,MB)*ICEIL(ACOLS,MB)*(2*MB+1)
      EXMEMB = ICEIL(BROWS,NB)*ICEIL(BCOLS,NB)*(2*NB+1)
      EXMEMC = ICEIL(AROWS,MB)*ICEIL(BCOLS,NB)*(MB+NB+1)
      EXMEME = EXMEMC
C
C     We need three extra matrix descriptors for shifts and for 
C     calculating the correct indicies in sub(A), sub(B) and sub(C).
C     These descriptors are related to the the smallest submatrices of 
C     A and B including sub(A) and sub(B) conforming with ScaLAPACK 
C     conventions.
C    
      CALL DESCINIT( DESCSA, M+IROFFA, M+IROFFA, DESCA(MB_), 
     $               DESCA(NB_), ARSRC, ACSRC, ICTXT, LLDA, INFO )
      CALL DESCINIT( DESCSB, N+IROFFB, N+IROFFB, DESCB(MB_), 
     $               DESCB(NB_), BRSRC, BCSRC, ICTXT, LLDB, INFO )
      CALL DESCINIT( DESCSC, M+IROFFA, N+IROFFB, DESCC(MB_), 
     $               DESCC(NB_), CRSRC, CCSRC, ICTXT, LLDC, INFO )
      CALL DESCINIT( DESCE, M+IROFFA, N+IROFFB, DESCC(MB_), 
     $               DESCC(NB_), CRSRC, CCSRC, ICTXT, MAX( 1, AROWS ), 
     $               INFO )
C
C     Adjust AROWS - BCOLS to sub(A) and sub(B) 
C
      IF( MYROW.EQ.ARSRC ) AROWS = AROWS - IROFFA
      IF( MYCOL.EQ.ACSRC ) ACOLS = ACOLS - IROFFA
      IF( MYROW.EQ.BRSRC ) BROWS = BROWS - IROFFB
      IF( MYCOL.EQ.BCSRC ) BCOLS = BCOLS - IROFFB
C
C     Test work space
C     
      IF( INFO.EQ.0 ) THEN
C
C     Compute the needed workspace for holding the matrix E
C
         SIZE_E = AROWS * BCOLS
         LLDE = DESCE(LLD_)
C
C     Compute workspace requirements for subsystem solving routine
C
         SIZE_IPW = (NB+1)*(MB2+1) 
C
C     Compute the needed workspace
C
         WRK = SIZE_E + 4*(MB+1)*(NB+1) + (MB+1)**2 +
     $        (NB+1)**2 + EXMEMA + EXMEMB + EXMEMC + EXMEME + 
     $        MAX( MB, NB ) + SIZE_IPW + 
     $        AUPBL * (MB+1) ** 2 + BUPBL * (NB+1) ** 2 
C     
C     Compute needed integer workspace
C     
         IWRK = DBA + DBB + 8 * MNPDIM
      END IF
C     
C     Check if the call to PTRSYDTD was a workspace query, if so
C     store the needed workspace in DWORK(1) and IWORK(1) and return
C     If not check if the supplied memory is big enough
C     
      IF( INFO.EQ.0 .OR. LQUERY ) THEN
         LQUERY = ( LDWORK.EQ.-1 )
         IF( WRK.GT.LDWORK .AND. .NOT.LQUERY ) THEN
            INFO = -19
         ELSEIF( IWRK.GT.LIWORK .AND. .NOT. LQUERY ) THEN
            INFO = -21
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
         CALL PXERBLA( ICTXT, 'PTRSYDTD', -INFO )
         RETURN
      END IF
C     
C     Init some local pointers into the DWORK-array
C     
      E    = 1
      XIJ  = E + SIZE_E
      AII  = XIJ + (MB+1) * (NB+1)
      BJJ  = AII + (MB+1) ** 2
      CIJ  = BJJ + (NB+1) ** 2
      EIJ  = CIJ + (MB+1) * (NB+1)
      SND  = EIJ + (MB+1) * (NB+1)
      EXA  = SND + MAX( MB, NB )
      EXB  = EXA + EXMEMA
      EXC  = EXB + EXMEMB
      EXE  = EXC + EXMEMC
      CKJ  = EXE + EXMEME
      AII2 = CKJ + (MB+1) * (NB+1)
      AUP  = AII2 + (MB+1) ** 2
      BUP  = AUP + AUPBL * (MB+1) ** 2
      IPW  = BUP + BUPBL * (NB+1) ** 2       
C
C     Init two local pointers into the IWORK-array
C
      EXAINF = 1
      EXBINF = EXAINF + DBA
      IROWS = EXBINF + DBB
      ICOLS = IROWS + MNPDIM
      IGSI  = ICOLS + MNPDIM
      IGSJ  = IGSI + MNPDIM
      IEXRW = IGSJ + MNPDIM
      IEXCL = IEXRW + MNPDIM
      IXRWS = IEXCL + MNPDIM
      IXRIND = IXRWS + MNPDIM
C
C     Check for 2x2 diagonal-block-split between any blocks of sub(A) or
C     sub(B) and set the extensions.
C
      CALL PDEXTCHK( M, A, IA, JA, DESCA, IWORK( EXAINF ), DBA, AEXT, 
     $               INFO )
C     
      CALL PDEXTCHK( N, B, IB, JB, DESCB, IWORK( EXBINF ), DBB, BEXT, 
     $               INFO )
C
C     
C     Do an implicit redistribution of the elements in sub(A), 
C     sub(B) and sub(C) based on the extension information just set
C
      IF( AEXT.OR.BEXT ) 
     $     CALL PDIMPRED( 'ABC', M, N, A, IA, JA, DESCA, B, IB, JB,
     $                    DESCB, C, IC, JC, DESCC, IWORK( EXAINF ),
     $                    IWORK( EXBINF ),  DWORK( EXA ), EXMEMA, 
     $                    DWORK( EXB ), EXMEMB, DWORK( EXC ), EXMEMC, 
     $                    DWORK( SND ), MAX( MB, NB ), INFO )
C     
C     Compute the number of block diagonals of C and some communication
C     directions
C     
      NROLL = DBA + DBB - 1
      NORTH = MOD(MYROW + NPROW - 1, NPROW)
      EAST =  MOD(MYCOL + 1, NPCOL)
      SOUTH = MOD(MYROW + 1, NPROW)
      WEST =  MOD(MYCOL + NPCOL - 1, NPCOL)
C 
C     Depending on the transposes, set some loop variables
C     
      IF( (.NOT.TRANA).AND.(.NOT.TRANB) ) THEN
         JS = 1
         IS = DBA
         IE = DBA       
      ELSEIF( TRANA.AND.(.NOT.TRANB) ) THEN
         JS = 1
         IS = 1
         IE = 1
      ELSEIF( (.NOT.TRANA).AND.TRANB ) THEN
         JS = DBB
         IS = DBA
         IE = DBA
      ELSEIF( TRANA.AND.TRANB ) THEN
         JS = DBB
         IS = 1
         IE = 1
      END IF
C     
      SNODE = .TRUE.
C
C     Compute the number of rounds to do in deep pipelining and
C     set looplimits depending on the transposes
C
      NIDEEP = ICEIL( MB, MB2 )
      IF( .NOT.TRANA  ) THEN
         IDEEPS = NIDEEP
         IDEEPE = 1
         IDEEPU = -1
      ELSE
         IDEEPS = 1
         IDEEPE = NIDEEP
         IDEEPU = 1
      END IF
C
C     Compute the local indicies where my part of sub(A), sub(B) 
C     and sub(C) begins.
C
      ASI = IA + MB * MOD( NPROW + MYROW - ARSRC, NPROW ) - IROFFA 
      ASJ = JA + MB * MOD( NPCOL + MYCOL - ACSRC, NPCOL ) - IROFFA
      BSI = IB + NB * MOD( NPROW + MYROW - BRSRC, NPROW ) - IROFFB
      BSJ = JB + NB * MOD( NPCOL + MYCOL - BCSRC, NPCOL ) - IROFFB
      CSI = IC + MB * MOD( NPROW + MYROW - CRSRC, NPROW ) - IROFFA
      CSJ = JC + NB * MOD( NPCOL + MYCOL - CCSRC, NPCOL ) - IROFFB
      CALL INFOG2L( ASI, ASJ, DESCA, NPROW, NPCOL, MYROW, MYCOL, 
     $              LIAS, LJAS, IDUM1, IDUM2 )
      CALL INFOG2L( BSI, BSJ, DESCB, NPROW, NPCOL, MYROW, MYCOL, 
     $              LIBS, LJBS, IDUM1, IDUM2 )
      CALL INFOG2L( CSI, CSJ, DESCC, NPROW, NPCOL, MYROW, MYCOL, 
     $              LICS, LJCS, IDUM1, IDUM2 )
C     
C     Compute what alternative to use in shifts - these alternatives
C     is not used in the restoration of the data in initial positions
C     after the solution has been computed, since then no temporary
C     matrix E is shifted along with the right hand side C.
C
      MINELEM = MIN(M**2 + N**2, MIN(M**2 + 2*M*N, N**2 + 2*M*N))
      SHALT1  = M**2 + N**2 .EQ. MINELEM 
      SHALT2  = .NOT. SHALT1 .AND. (M**2 + 2*M*N .EQ. MINELEM)
      SHALT3  = .NOT. (SHALT1.OR.SHALT2).AND.(N**2 + 2*M*N .EQ. MINELEM)
C     
C     Main loop over number of diagonals in C
C     
      DO 10 K = 1, NROLL
C        
         IF( (.NOT.TRANA).AND.(.NOT.TRANB) ) THEN
            IF ( K.GT.DBB) IS = IS - 1
         ELSEIF( TRANA.AND.(.NOT.TRANB) ) THEN
            IF ( K.GT.DBB) IE = IE + 1
         ELSEIF( (.NOT.TRANA).AND.TRANB ) THEN
            IF ( K.GT.DBB) IS = IS - 1
         ELSEIF( TRANA.AND.TRANB ) THEN
            IF ( K.GT.DBB) IE = IE + 1
         END IF
C
C     If shifts are supposed to be used, we do it right here
C
         IF( SHIFT .AND. NPROCS.GT.1 .AND. 
     $     ( SHALT1.OR.N.EQ.NB.OR.M.EQ.MB ) ) THEN
C     
C     Shift sub(A) East if NPCOL > 1 and A**T West if NPCOL > 1
C     Shift sub(B) North if NPROW > 1 and B**T South if NPROW > 1
C     When shifting, also send/receive the extension elements
C    
            IF ( NPCOL.GT.1.AND.M.NE.MB ) THEN
               IF( (.NOT.TRANA) .AND. (.NOT.TRANB) ) THEN
                  CALL DGESD2D( ICTXT, AROWS, ACOLS, 
     $                 A((LJAS-1)*LLDA+LIAS), LLDA, MYROW, 
     $                 EAST)
                  DESCSA(CSRC_) = MOD(DESCSA(CSRC_) + 1, NPCOL)
                  CALL DGERV2D( ICTXT, AROWS, ACOLS, 
     $                 A((LJAS-1)*LLDA+LIAS), LLDA, MYROW,
     $                 WEST )
                  IF( AEXT ) THEN
                     CALL DGESD2D( ICTXT, EXMEMA, 1, DWORK(EXA),
     $                    EXMEMA, MYROW,EAST )
                     CALL DGERV2D( ICTXT, EXMEMA, 1, DWORK(EXA),
     $                                EXMEMA, MYROW, WEST )
                  END IF
               ELSEIF( TRANA .AND. TRANB ) THEN
                  CALL DGESD2D( ICTXT, AROWS, ACOLS, 
     $                 A((LJAS-1)*LLDA+LIAS), LLDA, MYROW, 
     $                 WEST )   
                  DESCSA(CSRC_) = MOD(DESCSA(CSRC_) + NPCOL - 1, 
     $                 NPCOL )
                  CALL DGERV2D( ICTXT, AROWS, ACOLS,
     $                 A((LJAS-1)*LLDA+LIAS), LLDA, MYROW, 
     $                 EAST)
                  IF( AEXT ) THEN
                     CALL DGESD2D( ICTXT, EXMEMA, 1, DWORK(EXA), 
     $                    EXMEMA, MYROW, WEST )
                     CALL DGERV2D( ICTXT, EXMEMA, 1, DWORK(EXA), 
     $                    EXMEMA, MYROW,EAST )
                  END IF
               END IF
            END IF
            IF ( NPROW.GT.1.AND.N.NE.NB ) THEN
               IF( (.NOT.TRANA) .AND. (.NOT.TRANB) ) THEN
                  CALL DGESD2D( ICTXT, BROWS, BCOLS,
     $                 B((LJBS-1)*LLDB+LIBS), LLDB, NORTH,
     $                 MYCOL )
                  DESCSB(RSRC_) = MOD(DESCSB(RSRC_) + NPROW - 1,
     $                 NPROW )
                  CALL DGERV2D( ICTXT, BROWS, BCOLS,
     $                 B((LJBS-1)*LLDB+LIBS), LLDB, SOUTH,
     $                 MYCOL )
                  IF( BEXT ) THEN
                     CALL DGESD2D( ICTXT, EXMEMB, 1, DWORK(EXB), 
     $                    EXMEMB, NORTH, MYCOL )
                     CALL DGERV2D( ICTXT, EXMEMB, 1, DWORK(EXB),
     $                    EXMEMB, SOUTH, MYCOL )
                  END IF
               ELSEIF( TRANA .AND. TRANB ) THEN
                  CALL DGESD2D( ICTXT, BROWS, BCOLS,
     $                 B((LJBS-1)*LLDB+LIBS), LLDB, SOUTH,
     $                 MYCOL )
                  DESCSB(RSRC_) = MOD(DESCSB(RSRC_) + 1, NPROW)
                  CALL DGERV2D( ICTXT, BROWS, BCOLS,
     $                 B((LJBS-1)*LLDB+LIBS), LLDB, NORTH,
     $                 MYCOL )
                  IF( BEXT ) THEN
                     CALL DGESD2D( ICTXT, EXMEMB, 1, DWORK(EXB),
     $                    EXMEMB, SOUTH, MYCOL )
                     CALL DGERV2D( ICTXT, EXMEMB, 1, DWORK(EXB),
     $                    EXMEMB, NORTH, MYCOL )
                  END IF
               END IF
            END IF
         ELSEIF( SHIFT .AND. NPROCS.GT.1 .AND. SHALT2 ) THEN 
C     
C     Shift A SouthEast and (C,E) South  if NPROW > 1 
C     or Shift A**T NorthWest and (C,E) North if NPROW > 1 
C     
            IF( (.NOT.TRANA) .AND. (.NOT.TRANB) ) THEN
               CALL DGESD2D( ICTXT, AROWS, ACOLS,
     $              A((LJAS-1)*LLDA+LIAS), LLDA, SOUTH, 
     $              EAST )
               DESCSA(RSRC_) = MOD(DESCSA(RSRC_) + 1, NPROW)
               DESCSA(CSRC_) = MOD(DESCSA(CSRC_) + 1, NPCOL)
               CALL DGERV2D( ICTXT, AROWS, ACOLS, 
     $              A((LJAS-1)*LLDA+LIAS), LLDA, NORTH, 
     $              WEST )
               IF( AEXT ) THEN
                  CALL DGESD2D( ICTXT, EXMEMA, 1, DWORK(EXA),
     $                 EXMEMA, SOUTH, EAST )
                  CALL DGERV2D( ICTXT, EXMEMA, 1, DWORK(EXA),
     $                 EXMEMA, NORTH, WEST )
               END IF
            ELSEIF( TRANA .AND. TRANB ) THEN
               CALL DGESD2D( ICTXT, AROWS, ACOLS,
     $              A((LJAS-1)*LLDA+LIAS), LLDA, NORTH, 
     $              WEST )
               DESCSA(RSRC_) = MOD(DESCSA(RSRC_) + NPROW - 1, NPROW)
               DESCSA(CSRC_) = MOD(DESCSA(CSRC_) + NPCOL - 1, NPCOL)
               CALL DGERV2D( ICTXT, AROWS, ACOLS,
     $              A((LJAS-1)*LLDA+LIAS), LLDA, SOUTH, 
     $              EAST )
               IF( AEXT ) THEN
                  CALL DGESD2D( ICTXT, EXMEMA, 1, DWORK(EXA), EXMEMA, 
     $                 NORTH, WEST )
                  CALL DGERV2D( ICTXT, EXMEMA, 1, DWORK(EXA), EXMEMA,
     $                 SOUTH, EAST)
               END IF
            END IF
C     
C     Shift C
C     
            IF ( NPROW.GT.1 ) THEN
               IF( (.NOT.TRANA) .AND. (.NOT.TRANB) ) THEN
                  CALL DGESD2D( ICTXT, AROWS, BCOLS, 
     $                 C((LJCS-1)*LLDC+LICS), LLDC, SOUTH,
     $                 MYCOL )
                  DESCSC(RSRC_) = MOD(DESCSC(RSRC_) + 1, NPROW)
                  CALL DGERV2D( ICTXT, AROWS, BCOLS,
     $                 C((LJCS-1)*LLDC+LICS), LLDC, NORTH,
     $                 MYCOL )
                  IF( AEXT .OR. BEXT ) THEN
                     CALL DGESD2D( ICTXT, EXMEMC, 1, DWORK(EXC),
     $                    EXMEMC, SOUTH, MYCOL )
                     CALL DGERV2D( ICTXT, EXMEMC, 1, DWORK(EXC), 
     $                    EXMEMC, NORTH, MYCOL )
                  END IF
               ELSEIF( TRANA .AND. TRANB ) THEN
                  CALL DGESD2D( ICTXT, AROWS, BCOLS, 
     $                 C((LJCS-1)*LLDC+LICS), LLDC, NORTH,
     $                 MYCOL )
                  DESCSC(RSRC_) = MOD(DESCSC(RSRC_) + NPROW-1,NPROW)
                  CALL DGERV2D( ICTXT, AROWS, BCOLS,
     $                 C((LJCS-1)*LLDC+LICS), LLDC, SOUTH,
     $                 MYCOL )
                  IF( AEXT .OR. BEXT ) THEN
                     CALL DGESD2D( ICTXT, EXMEMC, 1, DWORK(EXC),
     $                    EXMEMC, NORTH, MYCOL )
                     CALL DGERV2D( ICTXT, EXMEMC, 1, DWORK(EXC),
     $                    EXMEMC, SOUTH, MYCOL )
                  END IF
               END IF
            END IF
C     
C     Shift E
C     
            IF ( NPROW.GT.1 ) THEN
               IF( (.NOT.TRANA) .AND. (.NOT.TRANB) ) THEN
                  CALL DGESD2D( ICTXT, AROWS, BCOLS, 
     $                 DWORK( E ), LLDE, SOUTH,
     $                 MYCOL )
                  DESCE(RSRC_) = MOD(DESCE(RSRC_) + 1, NPROW)
                  CALL DGERV2D( ICTXT, AROWS, BCOLS,
     $                 DWORK( E ), LLDE, NORTH,
     $                 MYCOL )
                  IF( AEXT .OR. BEXT ) THEN
                     CALL DGESD2D( ICTXT, EXMEME, 1, DWORK(EXE),
     $                    EXMEME, SOUTH, MYCOL )
                     CALL DGERV2D( ICTXT, EXMEME, 1, DWORK(EXE), 
     $                    EXMEME, NORTH, MYCOL )
                  END IF
               ELSEIF( TRANA .AND. TRANB ) THEN
                  CALL DGESD2D( ICTXT, AROWS, BCOLS, 
     $                 DWORK( E ), LLDE, NORTH,
     $                 MYCOL )
                  DESCE(RSRC_) = MOD(DESCE(RSRC_) + NPROW-1,NPROW)
                  CALL DGERV2D( ICTXT, AROWS, BCOLS,
     $                 DWORK( E ), LLDE, SOUTH,
     $                 MYCOL )
                  IF( AEXT .OR. BEXT ) THEN
                     CALL DGESD2D( ICTXT, EXMEME, 1, DWORK(EXE),
     $                    EXMEME, NORTH, MYCOL )
                     CALL DGERV2D( ICTXT, EXMEME, 1, DWORK(EXE),
     $                    EXMEME, SOUTH, MYCOL )
                  END IF
               END IF
            END IF
         ELSEIF ( SHIFT .AND. NPROCS.GT.1 .AND. SHALT3 ) THEN 
C     
C     Shift B NorthWest and C West if NPCOL > 1 
C     or shift B**T SouthEast and C East if NPCOL > 1 
C     
            IF( (.NOT.TRANA) .AND. (.NOT.TRANB) ) THEN
               CALL DGESD2D( ICTXT, BROWS, BCOLS,
     $              B((LJBS-1)*LLDB+LIBS), LLDB, NORTH, 
     $              WEST )
               DESCSB(RSRC_) = MOD(DESCSB(RSRC_) + NPROW - 1, NPROW )
               DESCSB(CSRC_) = MOD(DESCSB(CSRC_) + NPCOL - 1, NPCOL )
               CALL DGERV2D( ICTXT, BROWS, BCOLS,
     $              B((LJBS-1)*LLDB+LIBS), LLDB, SOUTH, 
     $              EAST )
               IF( BEXT ) THEN
                  CALL DGESD2D( ICTXT, EXMEMB, 1, DWORK(EXB), EXMEMB,
     $                 NORTH,WEST)
                  CALL DGERV2D( ICTXT, EXMEMB, 1, DWORK(EXB), EXMEMB,
     $                 SOUTH, EAST )
               END IF
            ELSEIF( TRANA .AND. TRANB ) THEN
               CALL DGESD2D( ICTXT, BROWS, BCOLS,
     $              B((LJBS-1)*LLDB+LIBS), LLDB, SOUTH, 
     $              EAST )
               DESCSB(RSRC_) = MOD( DESCSB(RSRC_) + 1, NPROW )
               DESCSB(CSRC_) = MOD( DESCSB(CSRC_) + 1, NPCOL )
               CALL DGERV2D( ICTXT, BROWS, BCOLS, 
     $              B((LJBS-1)*LLDB+LIBS), LLDB, NORTH, 
     $              WEST )
               IF( BEXT ) THEN
                  CALL DGESD2D( ICTXT, EXMEMB, 1, DWORK(EXB), EXMEMB,
     $                 SOUTH, EAST )
                  CALL DGERV2D( ICTXT, EXMEMB, 1, DWORK(EXB), EXMEMB,
     $                 NORTH,WEST)
               END IF
            END IF
C     
C     Shift C
C     
            IF ( NPCOL.GT.1 ) THEN
               IF( (.NOT.TRANA) .AND. (.NOT.TRANB) ) THEN
                  CALL DGESD2D( ICTXT, AROWS, BCOLS,
     $                 C((LJCS-1)*LLDC+LICS), LLDC, MYROW,
     $                 WEST )
                  DESCSC(CSRC_) = MOD(DESCSC(CSRC_) + NPCOL-1,NPCOL)
                  CALL DGERV2D( ICTXT, AROWS, BCOLS, 
     $                 C((LJCS-1)*LLDC+LICS), LLDC, MYROW, 
     $                 EAST )
                  IF( AEXT .OR. BEXT ) THEN
                     CALL DGESD2D( ICTXT, EXMEMC, 1, DWORK(EXC), 
     $                    EXMEMC, MYROW, WEST )
                     CALL DGERV2D( ICTXT, EXMEMC, 1, DWORK(EXC),
     $                    EXMEMC, MYROW, EAST )
                  END IF
               ELSEIF( TRANA .AND. TRANB ) THEN
                  CALL DGESD2D( ICTXT, AROWS, BCOLS,
     $                 C((LJCS-1)*LLDC+LICS), LLDC, MYROW,
     $                 EAST )
                  DESCSC(CSRC_) = MOD(DESCSC(CSRC_) + 1, NPCOL)
                  CALL DGERV2D( ICTXT, AROWS, BCOLS,
     $                 C((LJCS-1)*LLDC+LICS), LLDC, MYROW, 
     $                 WEST )
                  IF( AEXT .OR. BEXT ) THEN
                     CALL DGESD2D( ICTXT, EXMEMC, 1, DWORK(EXC),
     $                    EXMEMC, MYROW, EAST )
                     CALL DGERV2D( ICTXT, EXMEMC, 1, DWORK(EXC),
     $                    EXMEMC, MYROW, WEST )
                  END IF
               END IF
            END IF
C     
C     Shift E
C     
            IF ( NPCOL.GT.1 ) THEN
               IF( (.NOT.TRANA) .AND. (.NOT.TRANB) ) THEN
                  CALL DGESD2D( ICTXT, AROWS, BCOLS,
     $                 DWORK( E ), LLDE, MYROW,
     $                 WEST )
                  DESCE(CSRC_) = MOD(DESCE(CSRC_) + NPCOL-1,NPCOL)
                  CALL DGERV2D( ICTXT, AROWS, BCOLS, 
     $                 DWORK( E ), LLDE, MYROW, 
     $                 EAST )
                  IF( AEXT .OR. BEXT ) THEN
                     CALL DGESD2D( ICTXT, EXMEME, 1, DWORK(EXE), 
     $                    EXMEME, MYROW, WEST )
                     CALL DGERV2D( ICTXT, EXMEME, 1, DWORK(EXE),
     $                    EXMEME, MYROW, EAST )
                  END IF
               ELSEIF( TRANA .AND. TRANB ) THEN
                  CALL DGESD2D( ICTXT, AROWS, BCOLS,
     $                 DWORK( E ), LLDE, MYROW,
     $                 EAST )
                  DESCE(CSRC_) = MOD(DESCE(CSRC_) + 1, NPCOL)
                  CALL DGERV2D( ICTXT, AROWS, BCOLS,
     $                 DWORK( E ), LLDE, MYROW, 
     $                 WEST )
                  IF( AEXT .OR. BEXT ) THEN
                     CALL DGESD2D( ICTXT, EXMEME, 1, DWORK(EXE),
     $                    EXMEME, MYROW, EAST )
                     CALL DGERV2D( ICTXT, EXMEME, 1, DWORK(EXE),
     $                    EXMEME, MYROW, WEST )
                  END IF
               END IF
            END IF
         END IF
C     
C     Recompute the number of block columns of sub(A) and sub(B) held by 
C     my node  
C     
         IF( SHIFT ) THEN    
            ACOLS2 = NUMROC( M+IROFFA, MB, MYCOL, DESCSA(CSRC_), NPCOL )
            BCOLS2 = NUMROC( N+IROFFB, NB, MYCOL, DESCSB(CSRC_), NPCOL )
            NBCA = ICEIL( ACOLS2, MB )
            NBCB = ICEIL( BCOLS2, NB )
         END IF
C
         JJS = JS        
C     
C     Solve subsystems on current diagonal in parallel - do it in 
C     a number of PHASES depending on the process configuration
C     
         PHASES = ICEIL( IS-IE+1, MNPDIM )
         DO 20 PHASE = 1, PHASES, 1
            IIS = IS-(PHASE-1)*MNPDIM
            IIE = MAX(IS-PHASE*MNPDIM+1,IE)
            IF( (.NOT.TRANA).AND.(.NOT.TRANB) ) THEN
               JJS = JS - (PHASE-1)*MNPDIM
            ELSEIF( TRANA.AND.(.NOT.TRANB) ) THEN
               JJS = JS + (PHASE-1)*MNPDIM
            ELSEIF( (.NOT.TRANA).AND.TRANB ) THEN
               JJS = JS + (PHASE-1)*MNPDIM
            ELSEIF( TRANA.AND.TRANB ) THEN
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
C     Check if Aii is extended and set some variables describing Cij
C     for this particular solve. To handle submatrices, we distinguish
C     between I = 1 and I > 1.
C     
               IF( I.EQ.1 ) THEN
                  IF( IWORK(EXAINF).EQ.0 ) THEN
                     ROWS = MIN( MB - IROFFA, M )
                     GSI = 1 + IROFFA
                     EXROW = .FALSE.
                  ELSEIF( IWORK(EXAINF).EQ.1 ) THEN
                     ROWS = MB -IROFFA + 1
                     GSI = 1 + IROFFA
                     EXROW = .TRUE.
                  END IF
               ELSE
                  IF( IWORK(EXAINF+(I-1)).EQ.0 ) THEN
                     ROWS = MIN(MB, (M - MB + IROFFA) - (I - 2) * MB)
                     GSI = (I - 1) * MB + 1
                     EXROW = .FALSE.
                  ELSEIF( IWORK(EXAINF+(I-1)).EQ.1 ) THEN
                     ROWS = MB + 1
                     GSI = (I - 1) * MB + 1
                     EXROW = .TRUE.
                  ELSEIF( IWORK(EXAINF+(I-1)).EQ.2 ) THEN
                     ROWS = MIN(MB, (M - MB + IROFFA) - (I - 2) * MB) -1
                     GSI = (I - 1) * MB + 2
                     EXROW = .FALSE.
                  ELSEIF( IWORK(EXAINF+(I-1)).EQ.3 ) THEN
                     ROWS = MB
                     GSI = (I - 1) * MB + 2
                     EXROW = .TRUE.
                  END IF
               END IF
               IWORK( IROWS + IIS - I ) = ROWS
               IWORK( IGSI + IIS - I ) = GSI
               IWORK( IEXRW + IIS - I ) = ILG2NT(EXROW)
C     
C     Check if Bjj is extended and set some variables describing Cij
C     
               IF( J.EQ.1 ) THEN
                  IF( IWORK(EXBINF).EQ.0 ) THEN
                     COLS = MIN( NB - IROFFB, N )
                     GSJ = 1 + IROFFB
                     EXCOL = .FALSE.
                  ELSEIF( IWORK(EXBINF).EQ.1 ) THEN
                     COLS = NB - IROFFB + 1
                     GSJ = 1 + IROFFB
                     EXCOL = .TRUE.
                  END IF
               ELSE
                  IF( IWORK(EXBINF+(J-1)).EQ.0 ) THEN
                     COLS = MIN(NB, (N - NB + IROFFB ) - (J - 2) * NB)
                     GSJ = (J - 1) * NB + 1
                     EXCOL = .FALSE.
                  ELSEIF( IWORK(EXBINF+(J-1)).EQ.1 ) THEN
                     COLS = NB + 1
                     GSJ = (J - 1) * NB + 1
                     EXCOL = .TRUE.
                  ELSEIF( IWORK(EXBINF+(J-1)).EQ.2 ) THEN
                     COLS = MIN(NB, (N - NB + IROFFB ) - (J - 2) * NB)-1
                     GSJ = (J - 1) * NB + 2
                     EXCOL = .FALSE.
                  ELSEIF( IWORK(EXBINF+(J-1)).EQ.3 ) THEN
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
     $                       MYCOL, LIA, LJA, ARSRC, ACSRC )
               CALL INFOG2L( GSJ, GSJ, DESCSB, NPROW, NPCOL, MYROW, 
     $                       MYCOL, LIB, LJB, BRSRC, BCSRC )
               CALL INFOG2L( GSI, GSJ, DESCSC, NPROW, NPCOL, MYROW, 
     $                       MYCOL, IX, JX, CRSRC, CCSRC )
C
C     Submatrix construction is omitted most of the times for
C     deep pipelining
C
               IF( IDEEP.NE.IDEEPS ) GO TO 32
C     
C     Build the extended matrix Aii and send it to the process holding
C     Cij
C     
               IF( MYROW.EQ.ARSRC .AND. MYCOL.EQ.ACSRC ) THEN
                  CALL DBEXMAT( EXROW, EXROW, ROWS, ROWS, LIA, LJA, 
     $                          A((LJAS-1)*LLDA+LIAS), LLDA, NBCA, MB, 
     $                          MB, DWORK( EXA ), DWORK( AII ), 
     $                          MB + 1 )
                  IF( (ARSRC.NE.CRSRC).OR.(ACSRC.NE.CCSRC) ) THEN
                     CALL DGESD2D( ICTXT, ROWS, ROWS, DWORK( AII ), 
     $                             MB+1, CRSRC, CCSRC )
                  END IF
               END IF
C     
C     Build the extended matrix Bjj and send it to the process holding  
C     Cij
C     
               IF( MYROW.EQ.BRSRC .AND. MYCOL.EQ.BCSRC ) THEN
                  CALL DBEXMAT( EXCOL, EXCOL, COLS, COLS, LIB, LJB, 
     $                          B((LJBS-1)*LLDB+LIBS), LLDB, NBCB, NB,
     $                          NB, DWORK( EXB ), DWORK( BJJ ), NB + 1 )
                  IF( (BRSRC.NE.CRSRC).OR.(BCSRC.NE.CCSRC) ) THEN
                     CALL DGESD2D( ICTXT, COLS, COLS, DWORK( BJJ ), 
     $                             NB+1, CRSRC, CCSRC )
                  END IF
               END IF
C     
               IF( MYROW.EQ.CRSRC .AND. MYCOL.EQ.CCSRC ) THEN
                  IF( ARSRC.NE.CRSRC .OR. ACSRC.NE.CCSRC ) THEN
                     CALL DGERV2D( ICTXT, ROWS, ROWS, DWORK( AII ), 
     $                             MB+1, ARSRC, ACSRC )   
                  END IF
                  IF( BRSRC.NE.CRSRC .OR. BCSRC.NE.CCSRC ) THEN
                     CALL DGERV2D( ICTXT, COLS, COLS, DWORK( BJJ ), 
     $                             NB+1, BRSRC, BCSRC )
                  END IF
               END IF
C
C     Skipped submatrix construction in deep pipelining? 
C
 32            CONTINUE  
C     
C     Set some solution variables
C     
               SNODE = MYROW.EQ.CRSRC.AND.MYCOL.EQ.CCSRC
               SRSRC = CRSRC
               SCSRC = CCSRC
C     
C     Build the extended matrix Cij
C     
               IF( SNODE ) THEN
                  IF( IDEEP.EQ.IDEEPS )
     $                 CALL DBEXMAT( EXROW, EXCOL, ROWS, COLS, IX, JX, 
     $                 C((LJCS-1)*LLDC+LICS), LLDC, NBCB, MB, NB, 
     $                 DWORK( EXC ), DWORK( CIJ ), MB + 1 )   
C     
C     Solve sub-system op(Aii)*Xij*op(Bjj) +/- Xij = Cij 
C     
                  CALL DTRSYDT( 'S', TRANSA, TRANSB, ISGN, IDEEP, 
     $                 ROWS, MB2, COLS, DWORK( AII ), MB + 1, 
     $                 DWORK( BJJ ), NB + 1, DWORK( CIJ ), MB + 1, 
     $                 DWORK( IPW ), SIZE_IPW, XRWS, XRIND, SCALOC, 
     $                 LINFO )
                  IWORK( IXRWS + IIS - I ) = XRWS
                  IWORK( IXRIND + IIS - I ) = XRIND
C     
                  IF ( LINFO.NE.0 ) INFO = LINFO
C     
C     Copy the given solution back to the global matrix and update the
C     extension elements for the last round of deep pipelining
C     
                  IF( IDEEP.EQ.IDEEPE .AND. SCALOC.EQ.ONE )     
     $                 CALL DUBEXMA( EXROW, EXCOL, ROWS, COLS, IX, JX, 
     $                 C((LJCS-1)*LLDC+LICS), LLDC, NBCB, MB, 
     $                 NB, DWORK( EXC ), DWORK( CIJ ), MB + 1 )
C     
C     Warn for owerflow 
C     
                  IF( SCALOC.NE.ONE ) THEN
                     INFO = 3
                  END IF
C     
C     Let the solving process copy the solution into the DWORK( XIJ ) 
C     space as a preparation for the updates of the global solution
C     
                  IF( XRWS.NE.0 ) 
     $                 CALL DLACPY( 'All', XRWS, COLS, 
     $                 DWORK( CIJ + XRIND - 1), MB + 1,  
     $                 DWORK( XIJ + XRIND - 1), MB + 1 )
               END IF
 35            CONTINUE
C     
C     Update inner loop variable j (column index)
C     
               IF( (.NOT.TRANA).AND.(.NOT.TRANB) ) THEN
                  J = J - 1
               ELSEIF( TRANA.AND.(.NOT.TRANB) ) THEN
                  J = J + 1
               ELSEIF( (.NOT.TRANA).AND.TRANB ) THEN
                  J = J + 1
               ELSEIF( TRANA.AND.TRANB ) THEN
                  J = J - 1
               END IF   
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
               CALL INFOG2L( GSI, GSJ, DESCSC, NPROW, NPCOL, 
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
               IF( (.NOT.TRANA).AND.(.NOT.TRANB) ) THEN
                  J = J - 1
               ELSEIF( TRANA.AND.(.NOT.TRANB) ) THEN
                  J = J + 1
               ELSEIF( (.NOT.TRANA).AND.TRANB ) THEN
                  J = J + 1
               ELSEIF( TRANA.AND.TRANB ) THEN
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
                  CALL PDSCAL( M, SCALOC, C, IC, JC+KKK-1, DESCC, 1 )
                  CALL PDSCAL( M, SCALOC, DWORK( E ), 1, KKK, DESCE, 
     $                 1 )
 123           CONTINUE
               IF( AEXT.OR.BEXT ) THEN
                  CALL DSCAL( EXMEMC, SCALOC, DWORK( EXC ), 1 )
                  CALL DSCAL( EXMEME, SCALOC, DWORK( EXE ), 1 )
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
                  CALL INFOG2L( GSI, GSJ, DESCSC, NPROW, NPCOL, 
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
     $                       DWORK( CIJ + (KKK-1)*(MB+1) ), 1 )
 555                 CONTINUE
                     CALL DLACPY( 'All', XRWS, COLS, 
     $                    DWORK( XIJ + XRIND - 1), MB+1,  
     $                    DWORK( CIJ + XRIND - 1), MB+1 )
                     IF( IDEEP.EQ.IDEEPE ) THEN
                        EXROW = INT2LG(IWORK( IEXRW + IIS - I )) 
                        EXCOL = INT2LG(IWORK( IEXCL + IIS - I ))
                        CALL DUBEXMA( EXROW, EXCOL, ROWS, COLS, IX, 
     $                       JX, C((LJCS-1)*LLDC+LICS), LLDC, NBCB, MB, 
     $                       NB, DWORK( EXC ), DWORK( CIJ ), MB + 1 )
                     END IF
                  END IF
C
 37               CONTINUE
C
                  IF( (.NOT.TRANA).AND.(.NOT.TRANB) ) THEN
                     J = J - 1
                  ELSEIF( TRANA.AND.(.NOT.TRANB) ) THEN
                     J = J + 1
                  ELSEIF( (.NOT.TRANA).AND.TRANB ) THEN
                     J = J + 1
                  ELSEIF( TRANA.AND.TRANB ) THEN
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
                  CALL INFOG2L( GSI, GSJ, DESCSC, NPROW, NPCOL, MYROW, 
     $                 MYCOL, IX, JX, CRSRC, CCSRC )
C     
                  SNODE = MYROW.EQ.CRSRC.AND.MYCOL.EQ.CCSRC
                  SRSRC = CRSRC
                  SCSRC = CCSRC
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
                  IF( (.NOT.TRANA).AND.(.NOT.TRANB) ) THEN
                     J = J - 1
                  ELSEIF( TRANA.AND.(.NOT.TRANB) ) THEN
                     J = J + 1
                  ELSEIF( (.NOT.TRANA).AND.TRANB ) THEN
                     J = J + 1
                  ELSEIF( TRANA.AND.TRANB ) THEN
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
               CALL INFOG2L( GSI, GSJ, DESCSC, NPROW, NPCOL, 
     $              MYROW, MYCOL, IX, JX, CRSRC, CCSRC )
C     
               SNODE = MYROW.EQ.CRSRC.AND.MYCOL.EQ.CCSRC
C     
               IF( SNODE ) THEN
                  XRWS = IWORK( IXRWS + IIS - I )  
                  XRIND = IWORK( IXRIND + IIS - I )
                  IF( XRWS.GT.0 )
     $                 CALL DTRSYDT( 'U', TRANSA, TRANSB, ISGN, IDEEP, 
     $                 ROWS, MB2, COLS, DWORK( AII ), MB + 1, 
     $                 DWORK( BJJ ), NB + 1, DWORK( CIJ ), MB + 1, 
     $                 DWORK( IPW ), SIZE_IPW, XRWS, XRIND, SCALOC, 
     $                 LINFO )
               END IF
 48            CONTINUE
C     
C     Update inner loop variable j (column index)
C     
               IF( (.NOT.TRANA).AND.(.NOT.TRANB) ) THEN
                  J = J - 1
               ELSEIF( TRANA.AND.(.NOT.TRANB) ) THEN
                  J = J + 1
               ELSEIF( (.NOT.TRANA).AND.TRANB ) THEN
                  J = J + 1
               ELSEIF( TRANA.AND.TRANB ) THEN
                  J = J - 1
               END IF
 47         CONTINUE
C     
C     Update rest of global system wrt to current solution - first
C     compute intermediate values Eik = Eik + Xij*Bkj^T or
C     Eik = Eik + Xij*Bjk . Choose the right branch depending on 
C     the transposes in the equation.
C
            IF( .NOT. TRANB ) THEN
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
C     For every update we do we also need to find out the dimensions and
C     possible extensions of the submatrices involved.
C     
                  IF( MOD( I, 2 ).EQ.0 ) THEN
                     INDXS = J
                     IF( I.EQ.1 .AND. .NOT.TRANA ) INDXS = J + 1
                     IF( I.EQ.DBA .AND. TRANA )    INDXS = J + 1
                     INDXE = DBB
                     INDXU = 1
                  ELSE
                     INDXS = DBB
                     INDXE = J
                     IF( I.EQ.1 .AND. .NOT.TRANA ) INDXE = J + 1
                     IF( I.EQ.DBA .AND. TRANA )    INDXE = J + 1
                     INDXU = -1
                  END IF
                  KBUP = 1
                  DO 60 INDX = INDXS, INDXE, INDXU
                     IF( IWORK(EXBINF+(INDX-1)).EQ.0 ) THEN
                        BCOLS2 = MIN(NB, (N - NB + IROFFB) - 
     $                       (INDX - 2) * NB)
                        GSIND = (INDX - 1) * NB + 1
                        EEXCOL = .FALSE.
                     ELSEIF( IWORK(EXBINF+(INDX-1)).EQ.1 ) THEN
                        BCOLS2 = NB + 1
                        GSIND = (INDX - 1) * NB + 1
                        EEXCOL = .TRUE.
                     ELSEIF( IWORK(EXBINF+(INDX-1)).EQ.2 ) THEN
                        BCOLS2 = MIN(NB, (N - NB + IROFFB) - 
     $                       (INDX - 2) * NB) - 1
                        GSIND = (INDX - 1) * NB + 2
                        EEXCOL = .FALSE.
                     ELSEIF( IWORK(EXBINF+(INDX-1)).EQ.3 ) THEN
                        BCOLS2 = NB 
                        GSIND = (INDX - 1) * NB + 2
                        EEXCOL = .TRUE.
                     END IF
C     
C     Build and, if needed, communicate the submatrices 
C     
                     CALL INFOG2L( GSJ, GSIND, DESCSB, NPROW, NPCOL,
     $                             MYROW, MYCOL, LIB, LJB, RSRC, 
     $                             CSRC )
                     CALL INFOG2L( GSI, GSIND, DESCE, NPROW, NPCOL, 
     $                             MYROW, MYCOL, IX, JX, RRSRC, RCSRC )
C
C     Submatrix constuction is omitted most of the times for 
C     deep pipelining
C
                     IF( IDEEP.NE.IDEEPS ) GO TO 62
                     IF( MYROW.EQ.RSRC.AND.MYCOL.EQ.CSRC ) THEN
                        IF( (RSRC.NE.RRSRC).OR.(CSRC.NE.RCSRC) ) THEN
                           CALL DBEXMAT( EXCOL, EEXCOL, COLS, BCOLS2,
     $                          LIB, LJB, B((LJBS-1)*LLDB+LIBS), 
     $                          LLDB, NBCB, NB, NB, DWORK( EXB ), 
     $                          DWORK( BUP+(BUPBL-1)*(NB+1)**2 ), 
     $                          NB + 1 )
                           CALL DGESD2D( ICTXT, COLS, BCOLS2, 
     $                          DWORK( BUP+(BUPBL-1)*(NB+1)**2  ), 
     $                          NB + 1, RRSRC, RCSRC )
                        ELSE
                           CALL DBEXMAT( EXCOL, EEXCOL, COLS, BCOLS2,
     $                          LIB, LJB, B((LJBS-1)*LLDB+LIBS), 
     $                          LLDB, NBCB, NB, NB, DWORK( EXB ), 
     $                          DWORK( BUP+(KBUP-1)*(NB+1)**2 ), 
     $                          NB + 1 )
                        END IF
                     END IF
C
C     Skipped submatrix construction?
C
 62                  CONTINUE
C
                     IF( MYROW.EQ.RRSRC .AND. MYCOL.EQ.RCSRC ) THEN
C     
C     If we are on the beginning of a new line set Eik=[0] else
C     build the matrix E from previously saved results. For deep
C     pipelining, we do this only once per block.
C     
                        IF( J.EQ.1 .AND. IDEEP.EQ.IDEEPS ) THEN
                           CALL DLASET( 'All', MB + 1, NB + 1, ZERO, 
     $                                  ZERO, DWORK(EIJ), MB + 1 )
                        ELSE
                           CALL DBEXMAT( EXROW, EEXCOL, ROWS, BCOLS2, 
     $                          IX, JX, DWORK(E), LLDE, NBCB, MB, NB, 
     $                          DWORK( EXE ), DWORK(EIJ), MB + 1 ) 
                        END IF
                        IF( IDEEP.NE.IDEEPS ) GO TO 64
                        IF( (RSRC.NE.RRSRC) .OR. (CSRC.NE.RCSRC) ) THEN
                           CALL DGERV2D( ICTXT, COLS, BCOLS2, 
     $                          DWORK( BUP+(KBUP-1)*(NB+1)**2 ), 
     $                          NB + 1, RSRC, CSRC )
                        END IF
 64                     CONTINUE
C     
C     Perform the update of Eik
C     
                        XRWS = IWORK( IXRWS + IIS - I )  
                        XRIND = IWORK( IXRIND + IIS - I )
                        IF( XRWS.GT.0 .AND. BCOLS2.GT.0 ) 
     $                       CALL DGEMM( 'N', 'N', XRWS, BCOLS2, COLS, 
     $                       ONE, DWORK( XIJ+XRIND-1 ), MB + 1, 
     $                       DWORK(BUP+(KBUP-1)*(NB+1)**2),
     $                       NB + 1, ONE, DWORK(EIJ+XRIND-1), MB + 1 )
C
C     Save the result in the matrix E
C     
                        CALL DUBEXMA( EXROW, EEXCOL, ROWS, BCOLS2,
     $                       IX, JX, DWORK(E), LLDE, NBCB, MB,
     $                       NB, DWORK( EXE ), DWORK(EIJ),
     $                       MB + 1 )
                        IF( BUPBL.NE.1 ) KBUP = KBUP + 1
                     END IF
C
 60               CONTINUE
C
 55               CONTINUE
C
C     Update inner loop variable j (column index)
C
                  IF( .NOT.TRANA ) THEN
                     J = J - 1
                  ELSEIF( TRANA ) THEN
                     J = J + 1
                  END IF
C     
 50            CONTINUE
C     
            ELSEIF( TRANB ) THEN
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
                  IF( MOD( I, 2 ).EQ.0 ) THEN
                     INDXS = J
                     IF( I.EQ.1 .AND. .NOT.TRANA ) INDXS = J - 1
                     IF( I.EQ.DBA .AND. TRANA ) INDXS = J - 1
                     INDXE = 1
                     INDXU = -1
                  ELSE
                     INDXS = 1 
                     INDXE = J
                     IF( I.EQ.1 .AND. .NOT.TRANA ) INDXE = J - 1
                     IF( I.EQ.DBA .AND. TRANA ) INDXE = J - 1
                     INDXU = 1
                  END IF
                  KBUP = 1
                  DO 80 INDX =  INDXS, INDXE, INDXU
C     
C     Set some constants descibing the involved submatrices
C     
C     
                     IF( INDX.EQ.1 ) THEN
                        IF( IWORK(EXBINF).EQ.0 ) THEN
                           BROWS2 = NB - IROFFB
                           GSIND = 1 + IROFFB
                           EEXCOL = .FALSE.
                        ELSEIF( IWORK(EXBINF).EQ.1 ) THEN
                           BROWS2 = NB -IROFFB + 1
                           GSIND = 1 + IROFFB
                           EEXCOL = .TRUE.
                        END IF
                     ELSE
                        IF( IWORK(EXBINF+(INDX-1)).EQ.0 ) THEN
                           BROWS2 = MIN(NB, (N - NB + IROFFB) - 
     $                          (INDX - 2) * NB)
                           GSIND = (INDX - 1) * NB + 1
                           EEXCOL = .FALSE.
                        ELSEIF( IWORK(EXBINF+(INDX-1)).EQ.1 ) THEN
                           BROWS2 = NB + 1
                           GSIND = (INDX - 1) * NB + 1
                           EEXCOL = .TRUE.
                        ELSEIF( IWORK(EXBINF+(INDX-1)).EQ.2 ) THEN
                           BROWS2 = MIN(NB, (N - NB + IROFFB) - 
     $                          (INDX - 2) * NB) - 1
                           GSIND = (INDX - 1) * NB + 2
                           EEXCOL = .FALSE.
                        ELSEIF( IWORK(EXBINF+(INDX-1)).EQ.3 ) THEN
                           BROWS2 = NB 
                           GSIND = (INDX - 1) * NB + 2
                           EEXCOL = .TRUE.
                        END IF
                     END IF
C     
C     Build and, if needed, communicate the submatrices 
C     
                     CALL INFOG2L( GSIND, GSJ, DESCSB, NPROW, NPCOL,
     $                             MYROW, MYCOL, LIB, LJB, RSRC, CSRC )
                     CALL INFOG2L( GSI, GSIND, DESCE, NPROW, NPCOL, 
     $                             MYROW, MYCOL, IX, JX, RRSRC, RCSRC )
C
C     Submatrix constuction is omitted most of the times for 
C     deep pipelining
C
                     IF( IDEEP.NE.IDEEPS ) GO TO 82
                     IF( MYROW.EQ.RSRC.AND.MYCOL.EQ.CSRC ) THEN
                        IF( (RSRC.NE.RRSRC).OR.(CSRC.NE.RCSRC) ) THEN 
                           CALL DBEXMAT( EEXCOL, EXCOL, BROWS2, COLS, 
     $                          LIB, LJB, B((LJBS-1)*LLDB+LIBS), LLDB, 
     $                          NBCB, NB, NB, DWORK( EXB ), 
     $                          DWORK( BUP+(BUPBL-1)*(NB+1)**2 ), 
     $                          NB + 1 )
                           CALL DGESD2D( ICTXT, BROWS2, COLS, 
     $                          DWORK( BUP+(BUPBL-1)*(NB+1)**2 ), 
     $                          NB + 1, RRSRC, RCSRC )
                        ELSE
                           CALL DBEXMAT( EEXCOL, EXCOL, BROWS2, COLS, 
     $                          LIB, LJB, B((LJBS-1)*LLDB+LIBS), LLDB, 
     $                          NBCB, NB, NB, DWORK( EXB ), 
     $                          DWORK( BUP+(KBUP-1)*(NB+1)**2 ), 
     $                          NB + 1 )
                        END IF
                     END IF
C     
C    Skipped submatrix construction?
C     
 82                  CONTINUE
C     
                     IF(MYROW.EQ.RRSRC.AND.MYCOL.EQ.RCSRC) THEN
C     
C     If we are on the beginning of a new line set Eik=[0] else
C     build the matrix E from previously saved results. For 
C     
                        IF( J.EQ.DBB .AND. IDEEP.EQ.IDEEPS ) THEN
                           CALL DLASET( 'All', MB + 1, NB + 1, ZERO, 
     $                          ZERO, DWORK(EIJ), MB + 1 )
                        ELSE
                           CALL DBEXMAT( EXROW, EEXCOL, ROWS, BROWS2, 
     $                          IX, JX, DWORK(E), LLDE, NBCB, 
     $                          MB, NB, DWORK( EXE ), DWORK(EIJ), 
     $                          MB + 1 ) 
                        END IF
                        IF( IDEEP.NE.IDEEPS ) GO TO 84
                        IF( (RSRC.NE.RRSRC) .OR. (CSRC.NE.RCSRC) ) THEN
                           CALL DGERV2D( ICTXT, BROWS2, COLS, 
     $                          DWORK( BUP+(KBUP-1)*(NB+1)**2 ), 
     $                          NB + 1, RSRC, CSRC )
                        END IF
 84                     CONTINUE
C     
C     Perform the update of Eik
C     
                        XRWS = IWORK( IXRWS + IIS - I ) 
                        XRIND = IWORK( IXRIND + IIS - I )
                        IF( XRWS.GT.0 .AND. BROWS2.GT.0 )
     $                       CALL DGEMM( 'N','T', XRWS, BROWS2, COLS, 
     $                       ONE, DWORK( XIJ+XRIND-1 ), MB + 1, 
     $                       DWORK(BUP+(KBUP-1)*(NB+1)**2), NB + 1, 
     $                       ONE, DWORK(EIJ+XRIND-1), MB + 1 ) 
C     
C     Save the result in the matrix E
C     
                        CALL DUBEXMA( EXROW, EEXCOL, ROWS, BROWS2,
     $                       IX, JX, DWORK(E), LLDE, NBCB, MB,
     $                       NB, DWORK( EXE ), DWORK(EIJ),
     $                       MB + 1 )
                        IF( BUPBL.NE.1 ) KBUP = KBUP + 1
                     END IF
 80               CONTINUE
C     
 75               CONTINUE
C     
C     Update inner loop variable j (column index)
C     
                  IF( .NOT.TRANA ) THEN
                     J = J + 1
                  ELSEIF( TRANA ) THEN
                     J = J - 1
                  END IF
C     
 70            CONTINUE
            END IF
C
C     Broadcast the submatrix Eij in the corresponding block column
C
            IF( K.LT.NROLL .AND. .NOT.TRANB ) THEN
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
                  CALL INFOG2L( GSI, GSJ, DESCE, NPROW, NPCOL, MYROW, 
     $                          MYCOL, IX, JX, CRSRC, CCSRC )
C     
                  SNODE = MYROW.EQ.CRSRC.AND.MYCOL.EQ.CCSRC
                  SRSRC = CRSRC
                  SCSRC = CCSRC
C
                  IF( SNODE ) THEN
                     CALL DBEXMAT( EXROW, EXCOL, ROWS, COLS, IX, JX, 
     $                             DWORK(E), LLDE, NBCB, MB, NB, 
     $                             DWORK( EXE ), DWORK( EIJ ), MB + 1 )
                     XRWS = IWORK( IXRWS + IIS - I )  
                     XRIND = IWORK( IXRIND + IIS - I )
                     IBUFF(1) = XRWS
                     IBUFF(2) = XRIND
C     
                     IF (NPROW.GT.1) THEN
                        CALL IGEBS2D( ICTXT, 'Col', ' ', 2, 1, IBUFF, 2)
                        CALL DGEBS2D( ICTXT, 'Col', ' ', XRWS, COLS,
     $                       DWORK(EIJ+XRIND-1), MB + 1 )
                     END IF
                  ELSE
                     IF ( NPROW.GT.1 .AND. MYCOL.EQ.SCSRC ) THEN 
                        CALL IGEBR2D( ICTXT, 'Col', ' ', 2, 1, IBUFF, 2, 
     $                       SRSRC, SCSRC )
                        IWORK( IXRWS + IIS - I ) = IBUFF(1)
                        IWORK( IXRIND + IIS - I ) = IBUFF(2)
                        XRWS = IWORK( IXRWS + IIS - I ) 
                        XRIND = IWORK( IXRIND + IIS - I )
                        CALL DGEBR2D( ICTXT, 'Col', ' ', XRWS, COLS,
     $                       DWORK(EIJ+XRIND-1), MB + 1, SRSRC, SCSRC )
                     END IF
                  END IF
 95               CONTINUE
C     
C     Update inner loop variable j (column index)
C     
                  IF( .NOT.TRANA ) THEN
                     J = J - 1
                  ELSEIF( TRANA ) THEN
                     J = J + 1
                  END IF
 90            CONTINUE
C
            ELSEIF( K.LT.NROLL .AND. TRANB ) THEN
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
                  CALL INFOG2L( GSI, GSJ, DESCE, NPROW, NPCOL, MYROW, 
     $                          MYCOL, IX, JX, CRSRC, CCSRC )
C     
                  SNODE = MYROW.EQ.CRSRC.AND.MYCOL.EQ.CCSRC
                  SRSRC = CRSRC
                  SCSRC = CCSRC 
C
                  IF( SNODE ) THEN
                     CALL DBEXMAT( EXROW, EXCOL, ROWS, COLS, IX, JX, 
     $                             DWORK(E), LLDE, NBCB, MB, NB, 
     $                             DWORK( EXE ), DWORK( EIJ ), MB + 1 )
                     XRWS = IWORK( IXRWS + IIS - I ) 
                     XRIND = IWORK( IXRIND + IIS - I )
                     IBUFF(1) = XRWS
                     IBUFF(2) = XRIND
C     
                     IF (NPROW.GT.1) THEN
                        CALL IGEBS2D( ICTXT, 'Col', ' ', 2, 1, IBUFF, 2)
                        CALL DGEBS2D( ICTXT, 'Col', ' ', XRWS, COLS,
     $                       DWORK(EIJ+XRIND-1), MB + 1 )
                     END IF
                  ELSE
                     IF ( NPROW.GT.1 .AND. MYCOL.EQ.SCSRC ) THEN 
                        CALL IGEBR2D( ICTXT, 'Col', ' ', 2, 1, IBUFF, 2, 
     $                       SRSRC, SCSRC )
                        IWORK( IXRWS + IIS - I ) = IBUFF(1)
                        IWORK( IXRIND + IIS - I ) = IBUFF(2)
                        XRWS = IWORK( IXRWS + IIS - I ) 
                        XRIND = IWORK( IXRIND + IIS - I )
                        CALL DGEBR2D( ICTXT, 'Col', ' ', XRWS, COLS,
     $                       DWORK(EIJ+XRIND-1), MB + 1, SRSRC, SCSRC )
                     END IF
                  END IF
 105              CONTINUE
C     
C     Update inner loop variable j (column index)
C     
                  IF( (.NOT.TRANA) ) THEN
                     J = J + 1
                  ELSEIF( TRANA ) THEN
                     J = J - 1
                  END IF  
 100           CONTINUE                  
            END IF
C     
C     Go on with the updates:   
C     Now update block column j of C with respect to Eij in
C     the operation Ckj = Ckj - Aki*Eij or Ckj = Ckj - Aik^T*Eij 
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
               IF( .NOT.TRANA ) THEN
                  IF( MOD( J, 2 ).EQ.0 ) THEN
                     INDXS = 1
                     INDXE = I - 1
                     INDXU = 1
                  ELSE
                     INDXS = I - 1
                     INDXE = 1
                     INDXU = -1
                  END IF
                  KAUP = 1
                  DO 120 INDX = INDXS, INDXE, INDXU
C     
C     Set some constants descibing the involved submatrices
C     
                     IF( INDX.EQ.1 ) THEN
                        IF( IWORK(EXAINF).EQ.0 ) THEN
                           AROWS2 = MB - IROFFA
                           GSIND = 1 + IROFFA
                           CEXROW = .FALSE.
                        ELSEIF( IWORK(EXAINF).EQ.1 ) THEN
                           AROWS2 = MB - IROFFA + 1
                           GSIND = 1 + IROFFA
                           CEXROW = .TRUE.
                        END IF
                     ELSE
                        IF( IWORK(EXAINF+(INDX-1)).EQ.0 ) THEN
                           AROWS2 = MIN(MB, (M - MB + IROFFA) - 
     $                          (INDX - 2) * MB)
                           GSIND = (INDX - 1) * MB + 1
                           CEXROW = .FALSE.
                        ELSEIF( IWORK(EXAINF+(INDX-1)).EQ.1 ) THEN
                           AROWS2 = MB + 1
                           GSIND = (INDX - 1) * MB + 1
                           CEXROW = .TRUE.
                        ELSEIF( IWORK(EXAINF+(INDX-1)).EQ.2 ) THEN
                           AROWS2 = MIN(MB, (M - MB + IROFFA) - 
     $                          (INDX - 2) * MB) - 1
                           GSIND = (INDX - 1) * MB + 2
                           CEXROW = .FALSE.
                        ELSEIF( IWORK(EXAINF+(INDX-1)).EQ.3 ) THEN
                           AROWS2 = MB 
                           GSIND = (INDX - 1) * MB + 2
                           CEXROW = .TRUE.
                        END IF
                     END IF
C     
C     Build and, if needed, communicate the submatrices 
C     
                     CALL INFOG2L( GSIND, GSI, DESCSA, NPROW, NPCOL,
     $                             MYROW, MYCOL, LIA, LJA, RSRC, CSRC )
                     CALL INFOG2L( GSIND, GSJ, DESCSC, NPROW, NPCOL, 
     $                             MYROW, MYCOL, IX, JX, RRSRC, RCSRC )
C
C     Submatrix constuction is omitted most of the times for 
C     deep pipelining
C
                     IF( IDEEP.NE.IDEEPS ) GO TO 122
                     IF( MYROW.EQ.RSRC.AND.MYCOL.EQ.CSRC ) THEN
                        IF( (RSRC.NE.RRSRC).OR.(CSRC.NE.RCSRC) ) THEN
                           CALL DBEXMAT( CEXROW, EXROW, AROWS2, ROWS, 
     $                          LIA, LJA, A((LJAS-1)*LLDA+LIAS), LLDA, 
     $                          NBCA, MB, MB, DWORK( EXA ), 
     $                          DWORK(AUP+(AUPBL-1)*(MB+1)**2), MB + 1 )
                           CALL DGESD2D( ICTXT, AROWS2, ROWS, 
     $                          DWORK(AUP+(AUPBL-1)*(MB+1)**2), MB + 1, 
     $                          RRSRC, RCSRC )
                        ELSE
                           CALL DBEXMAT( CEXROW, EXROW, AROWS2, ROWS, 
     $                          LIA, LJA, A((LJAS-1)*LLDA+LIAS), LLDA, 
     $                          NBCA, MB, MB, DWORK( EXA ), 
     $                          DWORK(AUP+(KAUP-1)*(MB+1)**2), MB + 1 )
                        END IF
                     END IF
C     
C    Skipped submatrix construction?
C     
 122                 CONTINUE
                     IF(MYROW.EQ.RRSRC.AND.MYCOL.EQ.RCSRC) THEN
                        CALL DBEXMAT( CEXROW, EXCOL, AROWS2, COLS, IX, 
     $                                JX, C((LJCS-1)*LLDC+LICS), LLDC, 
     $                                NBCB, MB, NB, DWORK( EXC ), 
     $                                DWORK(CKJ), MB + 1 )
                        IF( IDEEP.NE.IDEEPS ) GO TO 124
                        IF( (RSRC.NE.RRSRC).OR.(CSRC.NE.RCSRC) ) THEN
                           CALL DGERV2D( ICTXT, AROWS2, ROWS, 
     $                          DWORK(AUP+(KAUP-1)*(MB+1)**2), MB + 1, 
     $                          RSRC, CSRC )
                        END IF
 124                    CONTINUE
C     
C     Perform the update of Ckj
C     
                        XRWS = IWORK( IXRWS + IIS - I )  
                        XRIND = IWORK( IXRIND + IIS - I )
                        CALL DGEMM( 'N', 'N', AROWS2, COLS, XRWS, 
     $                       -ONE, DWORK( AUP+(KAUP-1)*(MB+1)**2 + 
     $                       (XRIND-1)*(MB+1) ), MB + 1, 
     $                       DWORK(EIJ+XRIND-1), MB + 1, ONE, 
     $                       DWORK(CKJ), MB + 1 )
C     
C     Save the result in the matrix C
C     
                        CALL DUBEXMA( CEXROW, EXCOL, AROWS2, COLS, IX, 
     $                       JX, C((LJCS-1)*LLDC+LICS), LLDC, NBCB, MB, 
     $                       NB, DWORK( EXC ), DWORK(CKJ), MB + 1 )
                        IF( AUPBL.NE.1 ) KAUP = KAUP + 1
                     END IF
 120              CONTINUE
C     
               ELSEIF( TRANA ) THEN
                  IF( MOD( J, 2 ).EQ.0 ) THEN
                     INDXS = I + 1
                     INDXE = DBA
                     INDXU = 1
                  ELSE
                     INDXS = DBA
                     INDXE = I + 1
                     INDXU = -1
                  END IF
                  KAUP = 1
                  DO 130 INDX =  INDXS, INDXE, INDXU
C     
C     Set some constants descibing the involved submatrices
C     
                     IF( IWORK(EXAINF+(INDX-1)).EQ.0 ) THEN
                        ACOLS2 = MIN(MB, (M - MB + IROFFA) - 
     $                       (INDX - 2) * MB)
                        GSIND = (INDX - 1) * MB + 1
                        CEXROW = .FALSE.
                     ELSEIF( IWORK(EXAINF+(INDX-1)).EQ.1 ) THEN
                        ACOLS2 = MB + 1
                        GSIND = (INDX - 1) * MB + 1
                        CEXROW = .TRUE.
                     ELSEIF( IWORK(EXAINF+(INDX-1)).EQ.2 ) THEN
                        ACOLS2 = MIN(MB, (M - MB + IROFFA) - 
     $                       (INDX - 2) * MB) - 1
                        GSIND = (INDX - 1) * MB + 2
                        CEXROW = .FALSE.
                     ELSEIF( IWORK(EXAINF+(INDX-1)).EQ.3 ) THEN
                        ACOLS2 = MB 
                        GSIND = (INDX - 1) * MB + 2
                        CEXROW = .TRUE.
                     END IF
C     
C     Build and, if needed, communicate the submatrices 
C     
                     CALL INFOG2L( GSI, GSIND, DESCSA, NPROW, NPCOL,
     $                             MYROW, MYCOL, LIA, LJA, RSRC, CSRC )
                     CALL INFOG2L( GSIND, GSJ, DESCSC, NPROW, NPCOL, 
     $                             MYROW, MYCOL, IX, JX, RRSRC, RCSRC )
C
C     Submatrix constuction is omitted most of the times for 
C     deep pipelining
C
                     IF( IDEEP.NE.IDEEPS ) GO TO 132
C
                     IF( MYROW.EQ.RSRC.AND.MYCOL.EQ.CSRC ) THEN
                        IF( (RSRC.NE.RRSRC) .OR. (CSRC.NE.RCSRC) ) THEN
                           CALL DBEXMAT( EXROW, CEXROW, ROWS, ACOLS2, 
     $                          LIA, LJA, A((LJAS-1)*LLDA+LIAS), LLDA, 
     $                          NBCA, MB, MB, DWORK( EXA ), 
     $                          DWORK( AUP+(AUPBL-1)*(MB+1)**2 ), 
     $                          MB + 1 )
                           CALL DGESD2D( ICTXT, ROWS, ACOLS2, 
     $                          DWORK( AUP+(AUPBL-1)*(MB+1)**2 ), 
     $                          MB + 1, RRSRC, RCSRC )
                        ELSE
                           CALL DBEXMAT( EXROW, CEXROW, ROWS, ACOLS2, 
     $                          LIA, LJA, A((LJAS-1)*LLDA+LIAS), LLDA, 
     $                          NBCA, MB, MB, DWORK( EXA ), 
     $                          DWORK( AUP+(KAUP-1)*(MB+1)**2 ), 
     $                          MB + 1 )
                        END IF
                     END IF
C
C     Skipped submatrix construction?
C
 132                 CONTINUE
C
                     IF( MYROW.EQ.RRSRC .AND. MYCOL.EQ.RCSRC ) THEN
                        CALL DBEXMAT( CEXROW, EXCOL, ACOLS2, COLS, IX, 
     $                       JX, C((LJCS-1)*LLDC+LICS), LLDC, NBCB, MB, 
     $                       NB, DWORK( EXC ), DWORK(CKJ), MB + 1 )
                        IF( IDEEP.NE.IDEEPS ) GO TO 134
                        IF( (RSRC.NE.RRSRC) .OR. (CSRC.NE.RCSRC) ) THEN
                           CALL DGERV2D( ICTXT, ROWS, ACOLS2, 
     $                          DWORK( AUP+(KAUP-1)*(MB+1)**2 ), 
     $                          MB + 1, RSRC, CSRC )
                        END IF
 134                    CONTINUE
C     
C     Perform the update of Ckj
C     
                        XRWS = IWORK( IXRWS + IIS - I )  
                        XRIND = IWORK( IXRIND + IIS - I )
                        CALL DGEMM( 'T', 'N', ACOLS2, COLS, XRWS, -ONE,
     $                       DWORK(AUP+(KAUP-1)*(MB+1)**2+XRIND-1), 
     $                       MB + 1, DWORK(EIJ+XRIND-1), MB + 1, ONE, 
     $                       DWORK(CKJ), MB + 1 )
C     
C     Save the result in the matrix C
C     
                        CALL DUBEXMA( CEXROW, EXCOL, ACOLS2, COLS, IX, 
     $                       JX, C((LJCS-1)*LLDC+LICS), LLDC, NBCB, MB, 
     $                       NB, DWORK( EXC ), DWORK(CKJ), MB + 1 )
                        IF( AUPBL.NE.1 ) KAUP = KAUP + 1
                     END IF
 130              CONTINUE
               END IF
C
 115           CONTINUE
C     
C     Update inner loop variable j (column index)
C     
               IF( (.NOT.TRANA).AND.(.NOT.TRANB) ) THEN
                  J = J - 1
               ELSEIF( TRANA.AND.(.NOT.TRANB) ) THEN
                  J = J + 1
               ELSEIF( (.NOT.TRANA).AND.TRANB ) THEN
                  J = J + 1
               ELSEIF( TRANA.AND.TRANB ) THEN
                  J = J - 1
               END IF
C     
 110        CONTINUE               
C     
C     Prepare for solving for the next block diagonal of C:
C     update Ci,j-1 = Ci,j-1 - Aii * Ei,j-1 or
C     do update  Ci,j+1 = Ci,j+1 - Aii^T * Ei,j+1
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
               IF( (.NOT.TRANA .AND. TRANB) .OR. 
     $              (TRANA .AND. TRANB)) THEN
                  IF( J.GT.1 ) THEN
C     
C     Build all submatrices and, if needed, communicate the submatrix 
C     Aii involved.      
C     Check if Bj-1,j-1 is extended and set some variables describing 
C     Ci,j-1 and Ei,j-1
C     
                     IF( J.EQ.2 ) THEN
                        IF( IWORK(EXBINF).EQ.0 ) THEN
                           CCOLS = NB - IROFFB
                           GSIND = 1 + IROFFB
                           CEXCOL = .FALSE.
                        ELSEIF( IWORK(EXBINF).EQ.1 ) THEN
                           CCOLS = NB - IROFFB + 1
                           GSIND = 1 + IROFFB
                           CEXCOL = .TRUE.
                        END IF
                     ELSE
                        IF( IWORK(EXBINF+(J-2)).EQ.0 ) THEN
                           CCOLS = MIN(NB, (N - NB + IROFFB) - 
     $                          (J - 3) * NB)
                           GSIND = (J - 2) * NB + 1
                           CEXCOL = .FALSE.
                        ELSEIF( IWORK(EXBINF+(J-2)).EQ.1 ) THEN
                           CCOLS = NB + 1
                           GSIND = (J - 2) * NB + 1
                           CEXCOL = .TRUE.
                        ELSEIF( IWORK(EXBINF+(J-2)).EQ.2 ) THEN
                           CCOLS = MIN(NB, (N - NB + IROFFB) - 
     $                          (J - 3) * NB) - 1
                           GSIND = (J - 2) * NB + 2
                           CEXCOL = .FALSE.
                        ELSEIF( IWORK(EXBINF+(J-2)).EQ.3 ) THEN
                           CCOLS = NB 
                           GSIND = (J - 2) * NB + 2
                           CEXCOL = .TRUE.
                        END IF
                     END IF
C     
                     CALL INFOG2L( GSI, GSI, DESCSA, NPROW, NPCOL, 
     $                             MYROW, MYCOL, IIA, JJA, ARSRC, 
     $                             ACSRC )
                     CALL INFOG2L( GSI, GSIND, DESCSC, NPROW, NPCOL, 
     $                             MYROW, MYCOL, IIC, JJC, CRSRC, 
     $                             CCSRC )
                     CALL INFOG2L( GSI, GSIND, DESCE, NPROW, NPCOL, 
     $                             MYROW, MYCOL, EI, EJ, CRSRC, CCSRC )
C
C     Submatrix constuction is omitted most of the times for 
C     deep pipelining
C
                     IF( IDEEP.NE.IDEEPS ) GO TO 142
C     
                     IF( MYROW.EQ.ARSRC .AND. MYCOL.EQ.ACSRC ) THEN
                        CALL DBEXMAT( EXROW, EXROW, ROWS, ROWS, IIA, 
     $                       JJA, A((LJAS-1)*LLDA+LIAS), LLDA, NBCA, 
     $                       MB, MB, DWORK( EXA ), DWORK( AII2 ), 
     $                       MB + 1 )
                        IF( CRSRC.NE.ARSRC .OR. CCSRC.NE.ACSRC ) THEN
                           CALL DGESD2D( ICTXT, ROWS, ROWS,
     $                          DWORK( AII2 ), MB + 1, CRSRC, CCSRC )
                        END IF
                     END IF
C
C     Skipped submatrix construction?
C     
 142                 CONTINUE
C     
C     Do update
C     
                     IF( MYROW.EQ.CRSRC .AND. MYCOL.EQ.CCSRC ) THEN
                        CALL DBEXMAT( EXROW, CEXCOL, ROWS, CCOLS, IIC, 
     $                       JJC, C((LJCS-1)*LLDC+LICS), LLDC,
     $                       NBCB,  MB, NB, DWORK( EXC ), 
     $                       DWORK( CKJ ), MB + 1 )
                        CALL DBEXMAT( EXROW, CEXCOL, ROWS, CCOLS, EI, 
     $                       EJ, DWORK( E ), LLDE, NBCB, MB, 
     $                       NB, DWORK( EXE ), DWORK( EIJ ), 
     $                       MB + 1 )
                        IF( IDEEP.NE.IDEEPS ) GO TO 144
                        IF( CRSRC.NE.ARSRC .OR. CCSRC.NE.ACSRC ) THEN
                           CALL DGERV2D( ICTXT, ROWS, ROWS, 
     $                          DWORK( AII2 ), MB + 1, ARSRC, 
     $                          ACSRC )
                        END IF
 144                    CONTINUE
                        XRWS = IWORK( IXRWS + IIS - I ) 
                        XRIND = IWORK( IXRIND + IIS - I )
                        IF( .NOT. TRANA ) THEN
                           CALL DGEMM( 'N', 'N', ROWS, CCOLS, XRWS, 
     $                          -ONE, DWORK( AII2+(XRIND-1)*(MB+1) ), 
     $                          MB + 1, DWORK( EIJ+XRIND-1 ),
     $                          MB + 1, ONE, DWORK(CKJ), MB + 1 )
                        ELSE
                           CALL DGEMM( 'T', 'N', ROWS, CCOLS, XRWS, 
     $                          -ONE, DWORK( AII2 + XRIND - 1), 
     $                          MB + 1, DWORK( EIJ+XRIND-1 ),
     $                          MB + 1, ONE, DWORK(CKJ), MB + 1 )
                        END IF
C     
C     Save the result in the matrix C
C     
                        CALL DUBEXMA( EXROW, CEXCOL, ROWS, CCOLS, IIC, 
     $                                JJC, C((LJCS-1)*LLDC+LICS), LLDC, 
     $                                NBCB, MB, NB, DWORK( EXC ), 
     $                                DWORK(CKJ), MB + 1 ) 
                     END IF
                  END IF
C     
               ELSEIF( (TRANA .AND. .NOT.TRANB) .OR. 
     $                 (.NOT.TRANA .AND. .NOT.TRANB )) THEN
                  IF( J.LT.DBB ) THEN
C     
C     Build all submatrices and, if needed, communicate the submatrix 
C     Aii involved.      
C     Check if Bj+1,j+1 is extended and set some variables describing 
C     Ci,j+1 and Ei,j+1
C     
                     IF( IWORK(EXBINF+J).EQ.0 ) THEN
                        CCOLS = MIN(NB, (N - NB + IROFFB) - 
     $                       (J - 1) * NB)
                        GSIND = J * NB + 1
                        CEXCOL = .FALSE.
                     ELSEIF( IWORK(EXBINF+J).EQ.1 ) THEN
                        CCOLS = NB + 1
                        GSIND = J * NB + 1
                        CEXCOL = .TRUE.
                     ELSEIF( IWORK(EXBINF+J).EQ.2 ) THEN
                        CCOLS = MIN(NB, (N - NB + IROFFB) - 
     $                       (J-1) * NB) - 1
                        GSIND = J * NB + 2
                        CEXCOL = .FALSE.
                     ELSEIF( IWORK(EXBINF+J).EQ.3 ) THEN
                        CCOLS = NB 
                        GSIND = J * NB + 2
                        CEXCOL = .TRUE.
                     END IF
C     
                     CALL INFOG2L( GSI, GSI, DESCSA, NPROW, NPCOL, 
     $                             MYROW, MYCOL, IIA, JJA, ARSRC, 
     $                             ACSRC )
                     CALL INFOG2L( GSI, GSIND, DESCSC, NPROW, NPCOL, 
     $                             MYROW, MYCOL, IIC, JJC, CRSRC, 
     $                             CCSRC )
                     CALL INFOG2L( GSI, GSIND, DESCE, NPROW, NPCOL, 
     $                             MYROW, MYCOL, EI, EJ, CRSRC, CCSRC )
C
C     Submatrix constuction is omitted most of the times for 
C     deep pipelining
C
                     IF( IDEEP.NE.IDEEPS ) GO TO 146
C      
                     IF( MYROW.EQ.ARSRC .AND. MYCOL.EQ.ACSRC ) THEN
                        CALL DBEXMAT( EXROW, EXROW, ROWS, ROWS, IIA, 
     $                       JJA, A((LJAS-1)*LLDA+LICS), LLDA, NBCA, 
     $                       MB, MB, DWORK( EXA ), DWORK( AII2 ), 
     $                       MB + 1 )
                        IF( CRSRC.NE.ARSRC .OR. CCSRC.NE.ACSRC ) THEN
                           CALL DGESD2D( ICTXT, ROWS, ROWS, 
     $                          DWORK( AII2 ), MB + 1, CRSRC, CCSRC )
                        END IF
                     END IF
C
C     Skipped submatrix construction?
C
 146                 CONTINUE
C     
                     IF( MYROW.EQ.CRSRC .AND. MYCOL.EQ.CCSRC ) THEN
                        CALL DBEXMAT( EXROW, CEXCOL, ROWS, CCOLS, IIC, 
     $                                JJC, C((LJCS-1)*LLDC+LICS), LLDC, 
     $                                NBCB, MB, NB, DWORK( EXC ), 
     $                                DWORK( CKJ ), MB + 1 )
                        CALL DBEXMAT( EXROW, CEXCOL, ROWS, CCOLS, EI, 
     $                                EJ, DWORK( E ), LLDE, NBCB, MB, 
     $                                NB, DWORK( EXE ), DWORK( EIJ ), 
     $                                MB + 1 )
                        IF( IDEEP.NE.IDEEPS ) GO TO 148
                        IF( CRSRC.NE.ARSRC .OR. CCSRC.NE.ACSRC ) THEN
                           CALL DGERV2D( ICTXT, ROWS, ROWS, 
     $                          DWORK( AII2 ), MB + 1, ARSRC, 
     $                          ACSRC )
                        END IF     
C     
C     Do update
C     
 148                    CONTINUE
                        XRWS = IWORK( IXRWS + IIS - I ) 
                        XRIND = IWORK( IXRIND + IIS - I )
                        IF( TRANA ) THEN
                           CALL DGEMM( 'T', 'N', ROWS, CCOLS, XRWS, 
     $                          -ONE, DWORK(AII2+XRIND-1), MB + 1, 
     $                          DWORK( EIJ+XRIND-1 ), MB + 1, ONE, 
     $                          DWORK(CKJ), MB + 1 )
                        ELSE
                           CALL DGEMM( 'N', 'N', ROWS, CCOLS, XRWS, 
     $                          -ONE, DWORK(AII2+(MB+1)*(XRIND-1)), 
     $                          MB + 1, DWORK( EIJ+XRIND-1 ), MB + 1, 
     $                          ONE, DWORK(CKJ), MB + 1 )
                        END IF
C     
C     Save the result in the matrix C
C     
                        CALL DUBEXMA( EXROW, CEXCOL, ROWS, CCOLS, IIC, 
     $                       JJC, C((LJCS-1)*LLDC+LICS), LLDC, NBCB, 
     $                       MB, NB, DWORK( EXC ), DWORK(CKJ), MB + 1 ) 
                     END IF
                  END IF             
               END IF
C     
 145           CONTINUE
C     
C     Update inner loop variable j (column index)
C     
               IF( (.NOT.TRANA).AND.(.NOT.TRANB) ) THEN
                  J = J - 1
               ELSEIF( TRANA.AND.(.NOT.TRANB) ) THEN
                  J = J + 1
               ELSEIF( (.NOT.TRANA).AND.TRANB ) THEN
                  J = J + 1
               ELSEIF( TRANA.AND.TRANB ) THEN
                  J = J - 1
               END IF
C     
 140        CONTINUE
 24         CONTINUE
 20      CONTINUE
C     
C     Update outer loop variables
C     
         IF( (.NOT.TRANA).AND.(.NOT.TRANB) ) THEN
            JS = MIN( JS + 1, DBB )
            IE = MAX( 1, IE - 1 )
         ELSEIF( TRANA.AND.(.NOT.TRANB) ) THEN
            IF( K .GE. DBA ) JS = JS + 1
            IS = MIN( IS + 1, DBA )
         ELSEIF( (.NOT.TRANA).AND.TRANB ) THEN
            JS = MAX( JS - 1, 1 )
            IE = MAX( IE - 1, 1 )
         ELSEIF( TRANA.AND.TRANB ) THEN
            IF( K .GE. DBA ) JS = JS - 1
            IS = MIN( IS + 1, DBA )
         END IF
C     
 10   CONTINUE
C     
C     Before we go on we must back redistributed the elements in C
C     
      IF( AEXT.OR.BEXT) THEN
         CALL PDBCKRD( M, N, C, IC, JC, DESCC, IWORK( EXAINF ), 
     $                 IWORK( EXBINF ), DWORK( EXC ), EXMEMC,
     $                 DWORK( SND ), MAX( MB, NB ), INFO )
      END IF
C
C     Restore Data in Original Position. Now we don't need to
C     shift the extended elements or the intermediate matrices, we shift
C     just the original  matrices. 
C
      IF (SHIFT .AND. NPROCS.GT.1 .AND. 
     $     (SHALT1.OR.N.EQ.NB.OR.M.EQ.MB)) THEN
C     
C     Shift sub(A) East if NPCOL > 1 and A**T West if NPCOL > 1
C     Shift sub(B) North if NPROW > 1 and B**T South if NPROW > 1
C    
         IF ( NPCOL.GT.1.AND.M.NE.MB ) THEN
            IF( (.NOT.TRANA) .AND. (.NOT.TRANB) ) THEN
               CALL DGESD2D( ICTXT, AROWS, ACOLS, 
     $              A((LJAS-1)*LLDA+LIAS), LLDA, MYROW, 
     $              EAST )
               DESCSA(CSRC_) = MOD(DESCSA(CSRC_) + 1, NPCOL)
               CALL DGERV2D( ICTXT, AROWS, ACOLS, 
     $              A((LJAS-1)*LLDA+LIAS), LLDA, MYROW,
     $              WEST )
            ELSEIF( TRANA .AND. TRANB ) THEN
               CALL DGESD2D( ICTXT, AROWS, ACOLS, 
     $              A((LJAS-1)*LLDA+LIAS), LLDA, MYROW, 
     $              WEST )   
               DESCSA(CSRC_) = MOD(DESCSA(CSRC_) + NPCOL - 1, NPCOL)
               CALL DGERV2D( ICTXT, AROWS, ACOLS,
     $              A((LJAS-1)*LLDA+LIAS), LLDA, MYROW, 
     $              EAST )
            END IF
         END IF
         IF ( NPROW.GT.1.AND.N.NE.NB ) THEN
            IF( (.NOT.TRANA) .AND. (.NOT.TRANB) ) THEN
               CALL DGESD2D( ICTXT, BROWS, BCOLS,
     $              B((LJBS-1)*LLDB+LIBS), LLDB, NORTH,
     $              MYCOL )
               DESCSB(RSRC_) = MOD(DESCSB(RSRC_) + NPROW - 1,NPROW)
               CALL DGERV2D( ICTXT, BROWS, BCOLS,
     $              B((LJBS-1)*LLDB+LIBS), LLDB, SOUTH,
     $              MYCOL )
            ELSEIF( TRANA .AND. TRANB ) THEN
               CALL DGESD2D( ICTXT, BROWS, BCOLS,
     $              B((LJBS-1)*LLDB+LIBS), LLDB, SOUTH,
     $              MYCOL )
               DESCSB(RSRC_) = MOD(DESCSB(RSRC_) + 1, NPROW)
               CALL DGERV2D( ICTXT, BROWS, BCOLS,
     $              B((LJBS-1)*LLDB+LIBS), LLDB, NORTH,
     $              MYCOL )
            END IF
         END IF
      ELSEIF (SHIFT .AND. NPROCS.GT.1 .AND. SHALT2) THEN 
C     
C     Shift A SouthEast and C South  if NPROW > 1 
C     or Shift A**T NorthWest and C North if NPROW > 1 
C     
         IF( (.NOT.TRANA) .AND. (.NOT.TRANB) ) THEN
            CALL DGESD2D( ICTXT, AROWS, ACOLS,
     $           A((LJAS-1)*LLDA+LIAS), LLDA, SOUTH, EAST )
            DESCSA(RSRC_) = MOD(DESCSA(RSRC_) + 1, NPROW)
            DESCSA(CSRC_) = MOD(DESCSA(CSRC_) + 1, NPCOL)
            CALL DGERV2D( ICTXT, AROWS, ACOLS, 
     $           A((LJAS-1)*LLDA+LIAS), LLDA, NORTH, WEST )
         ELSEIF( TRANA .AND. TRANB ) THEN
            CALL DGESD2D( ICTXT, AROWS, ACOLS,
     $           A((LJAS-1)*LLDA+LIAS), LLDA, NORTH, WEST )
            DESCSA(RSRC_) = MOD(DESCSA(RSRC_) + NPROW - 1, NPROW)
            DESCSA(CSRC_) = MOD(DESCSA(CSRC_) + NPCOL - 1, NPCOL)
            CALL DGERV2D( ICTXT, AROWS, ACOLS,
     $           A((LJAS-1)*LLDA+LIAS), LLDA, SOUTH, EAST )
         END IF
         IF ( NPROW.GT.1 ) THEN
            IF( (.NOT.TRANA) .AND. (.NOT.TRANB) ) THEN
               CALL DGESD2D( ICTXT, AROWS, BCOLS, 
     $              C((LJCS-1)*LLDC+LICS), LLDC, SOUTH,
     $              MYCOL )
               DESCSC(RSRC_) = MOD(DESCSC(RSRC_) + 1, NPROW)
               CALL DGERV2D( ICTXT, AROWS, BCOLS,
     $              C((LJCS-1)*LLDC+LICS), LLDC, NORTH,
     $              MYCOL )
            ELSEIF( TRANA .AND. TRANB ) THEN
               CALL DGESD2D( ICTXT, AROWS, BCOLS, 
     $              C((LJCS-1)*LLDC+LICS), LLDC, NORTH,
     $              MYCOL )
               DESCSC(RSRC_) = MOD(DESCSC(RSRC_) + NPROW-1,NPROW)
               CALL DGERV2D( ICTXT, AROWS, BCOLS,
     $              C((LJCS-1)*LLDC+LICS), LLDC, SOUTH,
     $              MYCOL )
            END IF
         END IF
      ELSEIF (SHIFT .AND. NPROCS.GT.1 .AND. SHALT3) THEN 
C     
C     Shift B NorthWest and C West if NPCOL > 1 
C     or shift B**T SouthEast and C East if NPCOL > 1 
C     
         IF( (.NOT.TRANA) .AND. (.NOT.TRANB) ) THEN
            CALL DGESD2D( ICTXT, BROWS, BCOLS,
     $           B((LJBS-1)*LLDB+LIBS), LLDB, NORTH, WEST )
            DESCSB(RSRC_) = MOD(DESCSB(RSRC_) + NPROW - 1, NPROW )
            DESCSB(CSRC_) = MOD(DESCSB(CSRC_) + NPCOL - 1, NPCOL )
            CALL DGERV2D( ICTXT, BROWS, BCOLS,
     $           B((LJBS-1)*LLDB+LIBS), LLDB, SOUTH, EAST )
         ELSEIF( TRANA .AND. TRANB ) THEN
            CALL DGESD2D( ICTXT, BROWS, BCOLS,
     $           B((LJBS-1)*LLDB+LIBS), LLDB, SOUTH, EAST )
            DESCSB(RSRC_) = MOD( DESCSB(RSRC_) + 1, NPROW )
            DESCSB(CSRC_) = MOD( DESCSB(CSRC_) + 1, NPCOL )
            CALL DGERV2D( ICTXT, BROWS, BCOLS, 
     $           B((LJBS-1)*LLDB+LIBS), LLDB, NORTH, WEST )
         END IF
         IF ( NPCOL.GT.1 ) THEN
            IF( (.NOT.TRANA) .AND. (.NOT.TRANB) ) THEN
               CALL DGESD2D( ICTXT, AROWS, BCOLS,
     $              C((LJCS-1)*LLDC+LICS), LLDC, MYROW,
     $              WEST )
               DESCSC(CSRC_) = MOD(DESCSC(CSRC_) + NPCOL-1,NPCOL)
               CALL DGERV2D( ICTXT, AROWS, BCOLS, 
     $              C((LJCS-1)*LLDC+LICS), LLDC, MYROW, 
     $              EAST )
            ELSEIF( TRANA .AND. TRANB ) THEN
               CALL DGESD2D( ICTXT, AROWS, BCOLS,
     $              C((LJCS-1)*LLDC+LICS), LLDC, MYROW,
     $              EAST )
               DESCSC(CSRC_) = MOD(DESCSC(CSRC_) + 1, NPCOL)
               CALL DGERV2D( ICTXT, AROWS, BCOLS,
     $              C((LJCS-1)*LLDC+LICS), LLDC, MYROW, 
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
      END
C     
C     End of PTRSYDTD
C     
C *** Last line of PTRSYDTD ***
