CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C         This file is a part of the SCASY software library.           C
C         See http://www.cs.umu.se/parallel/research/scasy             C
C         Contributors: Robert Granat and Bo Kågström.                 C
C         SCASY Version 0.10,    March 31, 2006.                       C
C         Copyright 2006, Umeå University, Sweden.                     C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE PTRLYDTD( CSYMM, OP, N, A, IA, JA, DESCA, C, IC, JC,
     $                     DESCC, NB2, DWORK, LDWORK, IWORK, LIWORK, 
     $                     NOEXSY, SCALE, INFO)
C
C  -- ScaLAPACK-style routine (preliminary version)--
C     HPC2N and Department of Computing Science,
C     University of Umeå, Sweden.
C     Written by Robert Granat, (granat@cs.umu.se)
C     March 31, 2006.
C
      IMPLICIT NONE
C
C     .. Scalar Arguments ..
      LOGICAL            CSYMM
      CHARACTER*1        OP
      INTEGER            N, IA, JA, IC, JC, NB2, LDWORK, LIWORK, NOEXSY, 
     $                   INFO
      DOUBLE PRECISION   SCALE
C     ..
C     .. Array Arguments ..
      INTEGER            DESCA( * ), DESCC( * ), IWORK( * )
      DOUBLE PRECISION   A( * ), C( * ), DWORK( * )
C     ..
C
C  Purpose
C  =======
C
C  This subroutine solves the real (quasi-)triangular discrete-time 
C  Lyapunov equation
C
C     op( sub( A ) ) * X * op( sub( A ) )^T - X =  sub( C ),
C
C  where sub(A) = A(IA:IA+N-1,JA:JA+N-1), sub( C ) = 
C  C(IC:IC+N-1,JC:JC+N-1) and X are N-by-N distributed matrices.  
C
C  The matrix C (and X) may be symmetric (see CSYMM below).
C  The notation op(_) means the transpose or non-transpose of a matrix.
C
C  This routine should *not* be called directly, but through PGELYDTD.
C
C  NOTICE: this routine uses on demand communication only (see the 
C  references).
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
C  CSYMM     (global input) LOGICAL
C            If CSYMM = .TRUE., then the matrix C is assumed to be 
C            symmetric, that is, C = C^T. If CSYMM = .FALSE., the matrix 
C            C is not assumed to be symmetric.
C     
C  OP        (global input) CHARACTER*1
C            If OP = 'N', then we solve A * X * A^T - X = C
C            If OP = 'T', then we solve A^T * X * A - X = C.
C
C  Input/Output parameters
C
C  N         (global input) INTEGER
C            The number of rows and columns of the global distributed
C            matrix A. This is also the number of rows and columns of the
C            global distributed matrix C (and X). N >= 0.
C
C  A         (local input) DOUBLE PRECISION array 
C            Array of dimension (LLD_A,LOCc(M)). Contains the local
C            pieces of the global distributed matrix A in real Schur 
C            form. 
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
C  C         (local input/local output) DOUBLE PRECISION array 
C            Array of dimension (LLD_C,LOCc(N)). 
C            On entry, C contains the local pieces of the global 
C            distributed matrix C . On exit, it contains the local 
C            pieces of the global distributed solution X.
C
C  IC        (global input) INTEGER
C            Row start index for sub(C), i.e., the submatrix to operate 
C            on. MOD(IC,MB_A) = MOD(IA,MB_A) must hold. 
C
C  JC        (global input) INTEGER
C            Column start index for sub(C), i.e., the submatrix to 
C            operate on. MOD(JC,NB_A) = MOD(JA,NB_A) must hold.
C
C  DESCC     (global and local input) INTEGER array of dimension DLEN_.
C            The array descriptor for the global distributed matrix C.
C
C
C  NB2       (global input) INTEGER 
C            Internal blocking factor for pipelining of subsolutions
C            for updates of the matrix C (see the references 
C            for details).
C            1 < = NB2 <= DESCC( NB_ ) must hold.
C            Currently not used - placeholder for next release.
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
C            LIWORK >= DBA + 6 * MIN( P_r, P_c ), where 
C            DBA = ICEIL(LOCr(IA+IROFFA),MB_A).
C
C            If LIWORK = -1, LIWORK is global input and a workspace 
C            query is assumed. The routine will then calculate the 
C            optimal workspace needed, store it in IWORK(1) and return
C            immediately. No error will then be signaled by PXERBLA.
C            
C  Output information
C
C  NOEXSY    (local output) INTEGER
C            When solving the triangular problem in it is possible
C            that we have to extend some subsystems to not lose any data
C            from some 2x2 block of conjugate pairs of eigenvalues. 
C            NOEXSY helps us to keep track of the number of such 
C            extensions. 
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
C             If INFO = 2, the conditions for a unique solution was not
C             fulfilled; perturbed values was used to solve the equation 
C             (but the matrix A is unchanged).
C             If INFO = 3 the problem is badly scaled - C should have 
C             been scaled a factor SCALE before calling this routine to 
C             guarantee an overflow free solution, i.e., the solution
C             may well have overflowed.
C
C  Method
C  ======
C  This subroutine implements a parallel wave-front algorithm for
C  solving the triangular discrete-time Lyapunov equation. See
C  the references for details.
C
C  Additional requirements
C  =======================
C
C  A must be distributed using the same blocking factor in 
C  each direction, i.e., MB_A=NB_A. Moreover, C must be blocked 
C  with the same blocksize as A in both direction.
C
C  Limitations
C  ===========
C  In contrary to SLICOTs SB03MD this routine do not scale against 
C  overflow in the solution. See SCALE and INFO.
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
C  Wave-front algorithm, real Schur form, Lyapunov equation,
C  explicit blocking, GEMM-updates, on demand 
C
C  =====================================================================
C
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
      INTEGER  NB, DBA, NROLL, IS, JS, IE, K, RCMAX, MYCOL, MYROW,
     $         NPCOL, NPROW, J, NPROCS, I, IDUM, AROWS, ACOLS, ROWS, 
     $         COLS, LINFO, IX, JX, RSRC, CSRC, LIA, LJA, INDX, ICTXT,
     $         LLDA, LLDC, SRSRC, SCSRC, WRK, EXMEME, EXE, XIJ, AII, 
     $         AJJ, CIJ, EIJ, E, RRSRC, RCSRC, ARSRC, ACSRC, CRSRC, 
     $         CCSRC, EXA, EXC, EXMEMA, EXMEMC, IWRK, EXAINF, IPW,
     $         NBCA, POS, JEND, XIJT, AROWS2, GINDX, ACOLS2, D, SUBS, 
     $         LLDE, EROWS, ECOLS, STEP, ATRSRC, ATCSRC, LIAT, LJAT, 
     $         TJX, TIX, CTRSRC, TMP, CTCSRC, IDUMMY, SIZE_E, SIZE_SUBS, 
     $         LIC, LJC, CCOLS, LEI, LEJ, JSTART, ISTART, IEND, INDXS,
     $         INDXE, INDXU, IROFFA, GSIND, LIAS, LJAS, LJCS, LICS, KKK,
     $         ASI, ASJ, CSI, CSJ, GSI, GSJ, PHASES, PHASE, MNPDIM,
     $         IROWS, ICOLS, IEXRW, IEXCL, IGSI, IGSJ, IIS, JJS, IIE,
     $         THREADS
      DOUBLE PRECISION SCALOC, DUMMY, LSCALC
      LOGICAL  LQUERY, SNODE, RSIDE, EXROW, EXCOL, CEXROW, CEXCOL, 
     $         EEXCOL, EEXROW, AEXT
C     ..
C     .. Local Arrays ..
      INTEGER IDUM1(1), IDUM2(1), DESCSC( DLEN_ ), DESCSA( DLEN_ ),
     $        DESCE( DLEN_ )
      DOUBLE PRECISION DPDUM( 1 ), MACHINE(10)
C     ..
C     .. External Subroutines ..
      EXTERNAL           BLACS_GRIDINFO, CHK1MAT, SB04PY,
     $                   SB03MD, PXERBLA, INFOG2L, DGEMM, DLACPY,
     $                   PCHK2MAT, DSCAL, DGAMN2D, DGESD2D,
     $                   DGERV2D, IGAMX2D, DMATADD
C     ..
C     .. External Functions ..
      LOGICAL  LSAME, INT2LG
      INTEGER  NUMROC, ICEIL, OMP_GET_NUM_THREADS, ILG2NT
      EXTERNAL LSAME, NUMROC, ICEIL, OMP_GET_NUM_THREADS, INT2LG,
     $         ILG2NT
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC          MOD, MAX, MIN, ABS 
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
C
      CALL CHK1MAT( N, 3, N, 3, IA, JA, DESCA, 7, INFO )
      CALL CHK1MAT( N, 3, N, 3, IC, JC, DESCC, 11, INFO )
      CALL PCHK2MAT( N, 3, N, 3, IA, JA, DESCA, 7, N, 3, N, 3, 
     $               IC, JC, DESCC, 11, 0, IDUM1, IDUM2, INFO )
C
C     Check the value of OP
C 
      IF( INFO.EQ.0 ) THEN
         IF( .NOT.( LSAME( OP, 'T' ) .OR. LSAME( OP, 'N' ))) 
     $        THEN
            INFO = -2
         END IF
      END IF        
C
C     Check if workspace query
C
      LQUERY = LDWORK.EQ.-1 .OR. LIWORK.EQ.-1
C
C     Set some initial values
C
      SCALOC = ONE
      SCALE = ONE
      NPROCS = NPROW * NPCOL
      MNPDIM = MIN( NPROW, NPCOL )
      NB = DESCA( MB_ )
      IROFFA = MOD( IA - 1, NB )
      DBA = ICEIL( N + IROFFA, NB )
      AEXT = .FALSE.
      NOEXSY = 0
      LLDA = DESCA(LLD_)
      LLDC = DESCC(LLD_)
C
C     Set the logical value RSIDE which defines the equation to be solved.
C     If RSIDE = TRUE then the transposed matrix is located on the right
C     side of the solution X. Otherwise, RSIDE = FALSE.
C
      IF ( LSAME( OP, 'N') ) THEN
         RSIDE = .TRUE.
      ELSEIF( LSAME( OP, 'T' ) ) THEN
         RSIDE = .FALSE.
      END IF
C
C     Compute the number of rows and columns held by each process and
C     the number of block columns held by each process for the smallest
C     submatrix of A including sub(A) conforming with ScaLAPACK conventions 
C
      CALL INFOG2L( IA, JA, DESCA, NPROW, NPCOL, MYROW, MYCOL, LIA, LJA,
     $              ARSRC, ACSRC )
      CALL INFOG2L( IC, JC, DESCC, NPROW, NPCOL, MYROW, MYCOL, LIC, LJC,
     $              CRSRC, CCSRC )
      AROWS = NUMROC( N + IROFFA, NB, MYROW, ARSRC, NPROW )
      ACOLS = NUMROC( N + IROFFA, NB, MYCOL, ACSRC, NPCOL )
      NBCA = ICEIL( ACOLS, NB )
C
C     We need two extra matrix descriptors for computing indicies
C     in sub(A) and sub(C). These descriptors are related to the 
C     smallest submatrices of A and C including sub(A) and sub(C)
C     that conform with ScaLAPACK conventions
C
      CALL DESCINIT( DESCSA, N+IROFFA, N+IROFFA, DESCA(MB_), 
     $               DESCA(NB_), ARSRC, ACSRC, ICTXT, LLDA, INFO )
      CALL DESCINIT( DESCSC, N+IROFFA, N+IROFFA, DESCC(MB_), 
     $               DESCC(NB_), CRSRC, CCSRC, ICTXT, LLDC, INFO )
      CALL DESCINIT( DESCE, N+IROFFA, N+IROFFA, DESCC(MB_), 
     $               DESCC(NB_), CRSRC, CCSRC, ICTXT, MAX( 1, AROWS ), 
     $               INFO )
C
C     Compute the needed workspace for holding the matrix E
C
      SIZE_E = AROWS * ACOLS
      LLDE = DESCE(LLD_)
C
C     Compute workspace requirements for subsystem solving routines
C     SB04PY and SB03MD
C
      SIZE_SUBS = MAX( 2 * ( NB + 1 ), ( NB + 1 ) ** 2 )
C  
C     Compute the needed memory for holding data for the extended
C     subsystems. Even if all blocks cannot be extended, we add memory
C     for all blocks to maintain indicies simple.
C
      EXMEMA = ICEIL(AROWS,NB)*ICEIL(ACOLS,NB)*(2*NB+1)
      EXMEMC = EXMEMA
      EXMEME = EXMEMC
C
C     Adjust AROWS and ACOLS
C     
      IF( MYROW.EQ.ARSRC ) AROWS = AROWS - IROFFA
      IF( MYCOL.EQ.ACSRC ) ACOLS = ACOLS - IROFFA
C
C     Test work space
C     
      IF( INFO.EQ.0 ) THEN
C
         WRK = SIZE_E + 6 * (NB+1)**2 + EXMEMA + EXMEMC + EXMEME
     $         + NB + SIZE_SUBS 
      END IF
C     
C     Compute needed integer workspace
C     
      IF( INFO.EQ.0 ) THEN
         IWRK = DBA + 6 * MNPDIM
      END IF
C     
C     Check if the call to PTRLYDTD was a workspace query, if so
C     store the needed workspace in DWORK(1) and IWORK(1) and return
C     If not check if the supplied memory is big enough
C     
      IF( INFO.EQ.0 .OR. LQUERY ) THEN
         IF( WRK.GT.LDWORK .AND. .NOT.LQUERY ) THEN
            INFO = -14
         ELSEIF( IWRK.GT.LIWORK .AND. .NOT. LQUERY ) THEN
            INFO = -16
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
         CALL PXERBLA( ICTXT, 'PTRLYDTD', -INFO )
         RETURN
      END IF
C
C     Read out the default number of OpenMP threads
C
#ifdef USE_OMP
C$OMP PARALLEL SHARED(THREADS)
C$OMP MASTER
      THREADS = OMP_GET_NUM_THREADS()
C$OMP END MASTER
C$OMP END PARALLEL 
#else
      THREADS = 1
#endif
C     
C     Init some local pointers into the DWORK-array
C     
      E    = 1
      XIJ  = E + SIZE_E
      AII  = XIJ + (NB+1) ** 2
      AJJ  = AII + (NB+1) ** 2
      CIJ  = AJJ + (NB+1) ** 2
      EIJ  = CIJ + (NB+1) ** 2
      XIJT = EIJ + (NB+1) ** 2
      EXA  = XIJT + (NB+1) ** 2 
      EXC  = EXA + EXMEMA
      EXE  = EXC + EXMEMC
      IPW  = EXE + EXMEME
C
C     Init one local pointers into the IDWORK-array.  
C
      EXAINF = 1
      IROWS = EXAINF + DBA
      ICOLS = IROWS + MNPDIM
      IGSI  = ICOLS + MNPDIM
      IGSJ  = IGSI + MNPDIM
      IEXRW = IGSJ + MNPDIM
      IEXCL = IEXRW + MNPDIM
C
C     Check for 2x2 diagonal-block-split between any blocks of A and set
C     the extensions.
C
      CALL PDEXTCHK( N, A, IA, JA, DESCA, IWORK( EXAINF ), DBA,
     $               AEXT, INFO )
C     
C     
C     Do an implicit redistribution of the elements in A and C
C     based on the extension information just set
C     
      IF( AEXT )
     $     CALL PDIMPRED( 'A_C', N, N, A, IA, JA, DESCA, DPDUM, IDUM, 
     $                    IDUM, IDUM, C, IC, JC, DESCC, IWORK( EXAINF ),
     $                    IWORK( EXAINF ),  DWORK( EXA ), EXMEMA, DPDUM, 
     $                    IDUM, DWORK( EXC ), EXMEMC, DWORK( IPW ), NB, 
     $                    INFO )
C     
C     Compute the number of block diagonals of C
C     
      NROLL = 2 * DBA - 1
C 
C     Depending on OP (RSIDE), set some loop variables
C            
      IF( .NOT. RSIDE ) THEN 
         JS = 1
         IS = 1
         IE = 1
         STEP = + 1
      ELSEIF( RSIDE ) THEN
         JS = DBA
         IS = DBA
         IE = DBA
         STEP = -1
      END IF
C     
      SNODE = .TRUE.
C
C     Compute the local indicies where my part of sub(A) and sub(C) 
C     begins
C
      ASI = IA + NB * MOD( NPROW + MYROW - ARSRC, NPROW ) - IROFFA 
      ASJ = JA + NB * MOD( NPCOL + MYCOL - ACSRC, NPCOL ) - IROFFA
      CSI = IC + NB * MOD( NPROW + MYROW - CRSRC, NPROW ) - IROFFA
      CSJ = JC + NB * MOD( NPCOL + MYCOL - CCSRC, NPCOL ) - IROFFA
      CALL INFOG2L( ASI, ASJ, DESCA, NPROW, NPCOL, MYROW, MYCOL, 
     $              LIAS, LJAS, IDUM1, IDUM2 )
      CALL INFOG2L( CSI, CSJ, DESCC, NPROW, NPCOL, MYROW, MYCOL, 
     $              LICS, LJCS, IDUM1, IDUM2 )
C
C     Init machine
C
      CALL RECSY_MACHINE( MACHINE )
C     
C     Main loop over number of diagonals in C
C     
      DO 10 K = 1, NROLL
C        
         IF( .NOT. RSIDE ) THEN
            IF ( K.GT.DBA) IS = IS + 1
         ELSEIF( RSIDE ) THEN
            IF ( K.GT.DBA) IS = IS - 1
         END IF
C
         JJS = JS        
C     
C     Solve subsystems on current diagonal in parallel
C     
         PHASES = ICEIL( ABS(IS-IE)+1, MNPDIM )
         DO 20 PHASE = 1, PHASES, 1
            IF( .NOT. RSIDE ) THEN
               IIS = IS+(PHASE-1)*MNPDIM
               IIE = MIN(IS+PHASE*MNPDIM-1,IE)
               JJS = JS - (PHASE-1)*MNPDIM
            ELSEIF( RSIDE ) THEN 
               IIS = IS-(PHASE-1)*MNPDIM
               IIE = MAX(IS-PHASE*MNPDIM+1,IE)
               JJS = JS + (PHASE-1)*MNPDIM
            END IF
            J = JJS
            DO 30 I = IIS, IIE, STEP
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
                     ROWS = MIN( NB - IROFFA, N )
                     GSI = 1 + IROFFA
                     EXROW = .FALSE.
                  ELSEIF( IWORK(EXAINF).EQ.1 ) THEN
                     ROWS = NB - IROFFA + 1
                     GSI = 1 + IROFFA
                     EXROW = .TRUE.
                  END IF 
               ELSE
                  IF( IWORK(EXAINF+(I-1)).EQ.0 ) THEN
                     ROWS = MIN(NB, (N - NB + IROFFA) - (I - 2) * NB)
                     GSI = (I - 1) * NB + 1
                     EXROW = .FALSE.
                  ELSEIF( IWORK(EXAINF+(I-1)).EQ.1 ) THEN
                     ROWS = NB + 1
                     GSI = (I - 1) * NB + 1
                     EXROW = .TRUE.
                  ELSEIF( IWORK(EXAINF+(I-1)).EQ.2 ) THEN
                     ROWS = MIN(NB, (N - NB + IROFFA) - 
     $                    (I - 2) * NB) - 1
                     GSI = (I - 1) * NB + 2
                     EXROW = .FALSE.
                  ELSEIF( IWORK(EXAINF+(I-1)).EQ.3 ) THEN
                     ROWS = NB
                     GSI = (I - 1) * NB + 2
                     EXROW = .TRUE.
                  END IF
               END IF
               IWORK( IROWS + ABS(IIS - I) ) = ROWS
               IWORK( IGSI + ABS(IIS - I) ) = GSI
               IWORK( IEXRW + ABS(IIS - I) ) = ILG2NT(EXROW)
C     
C     Check if Ajj is extended and set some variables describing Cij
C     
               IF( J.EQ.1 ) THEN
                  IF( IWORK(EXAINF).EQ.0 ) THEN
                     COLS = MIN( NB - IROFFA, N )
                     GSJ = 1 + IROFFA
                     EXCOL = .FALSE.
                  ELSEIF( IWORK(EXAINF).EQ.1 ) THEN
                     COLS = NB - IROFFA + 1
                     GSJ = 1 + IROFFA
                     EXCOL = .TRUE.
                  END IF
               ELSE
                  IF( IWORK(EXAINF+(J-1)).EQ.0 ) THEN
                     COLS = MIN(NB, (N - NB + IROFFA) - (J - 2) * NB)
                     GSJ = (J - 1) * NB + 1
                     EXCOL = .FALSE.
                  ELSEIF( IWORK(EXAINF+(J-1)).EQ.1 ) THEN
                     COLS = NB + 1
                     GSJ = (J - 1) * NB + 1
                     EXCOL = .TRUE.
                  ELSEIF( IWORK(EXAINF+(J-1)).EQ.2 ) THEN
                     COLS = MIN(NB, (N - NB + IROFFA) - 
     $                    (J - 2) * NB) - 1
                     GSJ = (J - 1) * NB + 2
                     EXCOL = .FALSE.
                  ELSEIF( IWORK(EXAINF+(J-1)).EQ.3 ) THEN
                     COLS = NB 
                     GSJ = (J - 1) * NB + 2
                     EXCOL = .TRUE.
                  END IF
               END IF
               IWORK( ICOLS + ABS(IIS - I) ) = COLS
               IWORK( IGSJ + ABS(IIS - I) ) = GSJ
               IWORK( IEXCL + ABS(IIS - I) ) = ILG2NT(EXCOL)
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
C     If the matrix C is symmetric we solve only for the blocks located
C     on or below/above the main block diagonal of C.
C     
               IF(  (CSYMM .AND. RSIDE .AND. J.LE.I ) .OR. 
     $              (CSYMM. AND. (.NOT. RSIDE) .AND. J.GE.I ) 
     $              .OR. (.NOT. CSYMM) ) THEN
C     
C     Get starting indicies and the process id:s needed etc.
C     NOTICE: here we use the descriptors for the submatrices instead,
C     which means that all indicies below are calculated for the
C     smallest submatrices sub(A) and sub(C) that conform with ScaLAPACK
C     conventions
C     
                  CALL INFOG2L( GSI, GSI, DESCSA, NPROW, NPCOL, MYROW, 
     $                          MYCOL, LIA, LJA, ARSRC, ACSRC )
                  CALL INFOG2L( GSJ, GSJ, DESCSA, NPROW, NPCOL, MYROW, 
     $                          MYCOL, LIAT, LJAT, ATRSRC, ATCSRC )
                  CALL INFOG2L( GSI, GSJ, DESCSC, NPROW, NPCOL, MYROW, 
     $                          MYCOL, IX, JX, CRSRC, CCSRC )
C     
C     Build the extended matrix Aii and send it to the process holding
C     Cij
C     
                  IF( MYROW.EQ.ARSRC .AND. MYCOL.EQ.ACSRC ) THEN
                     CALL DBEXMAT( EXROW, EXROW, ROWS, ROWS, LIA, LJA, 
     $                             A((LJAS-1)*LLDA+LIAS), LLDA, NBCA, 
     $                             NB, NB, DWORK( EXA ), DWORK( AII ), 
     $                             NB + 1 )
                     IF( (ARSRC.NE.CRSRC).OR.(ACSRC.NE.CCSRC) ) THEN
                        CALL DGESD2D( ICTXT, ROWS, ROWS, DWORK( AII ),
     $                                NB + 1, CRSRC, CCSRC )
                     END IF
                  END IF
C     
C     Build the extended matrix Ajj and send it to the process holding  
C     Cij. This is only done if the I /= J.
C     
                  IF( I.NE.J .AND. MYROW.EQ.ATRSRC .AND. 
     $                 MYCOL.EQ.ATCSRC ) THEN
                     CALL DBEXMAT( EXCOL, EXCOL, COLS, COLS, LIAT, LJAT, 
     $                             A((LJAS-1)*LLDA+LIAS), LLDA, NBCA, 
     $                             NB, NB, DWORK( EXA ), DWORK( AJJ ), 
     $                             NB + 1 )
                     IF( (ATRSRC.NE.CRSRC) .OR. (ATCSRC.NE.CCSRC) ) THEN
                        CALL DGESD2D( ICTXT, COLS, COLS, DWORK( AJJ ), 
     $                                NB + 1, CRSRC, CCSRC )
                     END IF
                  END IF
C     
C     Receive Aii and Ajj
C     
                  IF( MYROW.EQ.CRSRC .AND. MYCOL.EQ.CCSRC ) THEN
                     IF( ARSRC.NE.CRSRC .OR. ACSRC.NE.CCSRC ) THEN
                        CALL DGERV2D( ICTXT, ROWS, ROWS, DWORK( AII ), 
     $                                NB + 1, ARSRC, ACSRC )
                     END IF
                     IF( I.EQ.J .AND. .NOT. CSYMM ) THEN
                        CALL DLACPY( 'All', COLS, COLS, DWORK( AII ),
     $                               NB + 1, DWORK( AJJ ), NB + 1 )
                     ELSEIF( I.NE.J ) THEN
                        IF( ATRSRC.NE.CRSRC .OR. ATCSRC.NE.CCSRC ) THEN
                           CALL DGERV2D( ICTXT, COLS, COLS, 
     $                                   DWORK( AJJ ), NB + 1, ATRSRC, 
     $                                   ATCSRC )
                        END IF
                     END IF
                  END IF  
C     
C     Set some solution variables
C     
                  SNODE = MYROW.EQ.CRSRC.AND.MYCOL.EQ.CCSRC
                  SRSRC = CRSRC
                  SCSRC = CCSRC
C     
C     Find out who owns the transposed subsolution
C     
                  CALL INFOG2L( GSJ, GSI, DESCSC, NPROW, NPCOL, MYROW, 
     $                          MYCOL, TIX, TJX, CTRSRC, CTCSRC )
C     
C     Build the extended matrix Cij 
C     
                  IF( SNODE ) THEN    
                     CALL DBEXMAT( EXROW, EXCOL, ROWS, COLS, IX, JX, 
     $                    C((LJCS-1)*LLDC+LICS), LLDC, NBCA, NB, NB, 
     $                    DWORK( EXC ), DWORK( CIJ ), NB + 1 )
                     IF( CSYMM .AND. I.EQ.J ) 
     $                    CALL DLATCPY( 'Lower', ROWS-1, COLS-1, 
     $                    DWORK(CIJ+1), NB+1, DWORK(CIJ+NB+1), NB+1 )
C     
C     Solve sub-system op(Aii)*Xij*op(Ajj^T) +/- Xij = Cij 
C     
                     IF( .NOT. CSYMM .OR. I.NE.J ) THEN
                        IF( I.NE.J ) THEN
#ifdef USE_OMP
                           IF( RSIDE ) THEN
                              CALL RECSYDT_P( THREADS, 3, SCALOC, 
     $                             ROWS, COLS, DWORK( AII ), NB + 1, 
     $                             DWORK(AJJ), NB + 1, DWORK( CIJ ), 
     $                             NB + 1, LINFO, MACHINE, DWORK(IPW),
     $                             LDWORK-IPW+1 )
                           ELSEIF( .NOT. RSIDE ) THEN
                              CALL RECSYDT_P( THREADS, 5, SCALOC, 
     $                             ROWS, COLS, DWORK( AII ), NB + 1, 
     $                             DWORK(AJJ), NB + 1, DWORK( CIJ ), 
     $                             NB + 1, LINFO, MACHINE, DWORK(IPW),
     $                             LDWORK-IPW+1 )
                           END IF
#else
                           IF( RSIDE ) THEN
                              CALL RECSYDT( 3, SCALOC, ROWS, COLS, 
     $                             DWORK( AII ), NB + 1, DWORK(AJJ), 
     $                             NB + 1, DWORK( CIJ ), NB + 1, LINFO,
     $                             MACHINE, DWORK(IPW), LDWORK-IPW+1 )
                           ELSEIF( .NOT. RSIDE ) THEN
                              CALL RECSYDT( 5, SCALOC, ROWS, COLS, 
     $                             DWORK( AII ), NB + 1, DWORK(AJJ), 
     $                             NB + 1, DWORK( CIJ ), NB + 1, LINFO,
     $                             MACHINE, DWORK(IPW), LDWORK-IPW+1 )
                           END IF
#endif
                        ELSEIF( I.EQ.J ) THEN
#ifdef USE_OMP
                           IF( RSIDE ) THEN
                              CALL RECSYDT_P( THREADS, 3, SCALOC, 
     $                             ROWS, COLS, DWORK( AII ), NB + 1, 
     $                             DWORK(AII), NB + 1, DWORK( CIJ ), 
     $                             NB + 1, LINFO, MACHINE, DWORK(IPW),
     $                             LDWORK-IPW+1 )
                           ELSEIF( .NOT. RSIDE ) THEN
                              CALL RECSYDT_P( THREADS, 5, SCALOC, 
     $                             ROWS, COLS, DWORK( AII ), NB + 1, 
     $                             DWORK(AII), NB + 1, DWORK( CIJ ), 
     $                             NB + 1, LINFO, MACHINE, DWORK(IPW),
     $                             LDWORK-IPW+1 )
                           END IF
#else
                           IF( RSIDE ) THEN
                              CALL RECSYDT( 3, SCALOC, ROWS, COLS, 
     $                             DWORK( AII ), NB + 1, DWORK(AII), 
     $                             NB + 1, DWORK( CIJ ), NB + 1, LINFO,
     $                             MACHINE, DWORK(IPW), LDWORK-IPW+1 )
                           ELSEIF( .NOT. RSIDE ) THEN
                              CALL RECSYDT( 5, SCALOC, ROWS, COLS, 
     $                             DWORK( AII ), NB + 1, DWORK(AII), 
     $                             NB + 1, DWORK( CIJ ), NB + 1, LINFO,
     $                             MACHINE, DWORK(IPW), LDWORK-IPW+1 )
                           END IF
#endif
                        END IF
                     ELSEIF( I.EQ.J .AND. CSYMM ) THEN
#ifdef USE_OMP
                        IF( RSIDE ) THEN 
                           CALL RECLYDT_P( THREADS, 0, SCALOC, ROWS, 
     $                          DWORK( AII ), NB + 1, DWORK( CIJ ), 
     $                          NB + 1, LINFO, MACHINE, DWORK(IPW),
     $                          LDWORK-IPW+1 )
                        ELSEIF( .NOT. RSIDE ) THEN
                           CALL RECLYDT_P( THREADS, 1, SCALOC, ROWS, 
     $                          DWORK( AII ), NB + 1, DWORK( CIJ ), 
     $                          NB + 1, LINFO, MACHINE, DWORK(IPW),
     $                          LDWORK-IPW+1 )
                        END IF
#else
                        IF( RSIDE ) THEN 
                           CALL RECLYDT( 0, SCALOC, ROWS, 
     $                          DWORK( AII ), NB + 1, DWORK( CIJ ), 
     $                          NB + 1, LINFO, MACHINE, DWORK(IPW),
     $                          LDWORK-IPW+1 )
                        ELSEIF( .NOT. RSIDE ) THEN
                           CALL RECLYDT( 1, SCALOC, ROWS, 
     $                          DWORK( AII ), NB + 1, DWORK( CIJ ), 
     $                          NB + 1, LINFO, MACHINE, DWORK(IPW),
     $                          LDWORK-IPW+1 )
                        END IF
#endif
                        CALL DLATCPY( 'Upper', ROWS-1, ROWS-1, 
     $                       DWORK(CIJ+NB+1), NB+1, DWORK(CIJ+1), NB+1 )
                     END IF
C     
                     IF (LINFO.NE.0) INFO = LINFO
C     
C     Copy the given solution back to the global matrix and update the
C     extension elements
C     
                     IF( SCALOC.EQ.ONE )
     $                    CALL DUBEXMA( EXROW, EXCOL, ROWS, COLS, IX, 
     $                    JX, C((LJCS-1)*LLDC+LICS), LLDC, NBCA, NB, NB, 
     $                    DWORK( EXC ), DWORK( CIJ ), NB + 1 )
C     
C     Warn for owerflow 
C     
                     IF (SCALOC.NE.ONE) THEN
                        INFO = 3
                     END IF
C     
C     Let the solving process copy the solution into the DWORK( XIJ ) 
C     space as a preparation for the updates of the global solution.
C     
                     CALL DLACPY( 'All', ROWS, COLS, DWORK( CIJ ), 
     $                    NB + 1, DWORK( XIJ ), NB + 1 )
                  END IF
               END IF
 35            CONTINUE
               IF( .NOT. RSIDE ) THEN
                  J = J - 1
               ELSEIF( RSIDE ) THEN
                  J = J + 1
               END IF
 30         CONTINUE
C
C     Compute global value of SCALOC by a k-to-all reduction, where 
C     k = MNPDIM. First, do a broadcast of the local SCALOC in the 
C     largest processor mesh dimension
C
            J = JJS
            DO 111 I = IIS, IIE, STEP
C
               IF(  (CSYMM .AND. RSIDE .AND. J.LE.I ) .OR. 
     $              (CSYMM. AND. (.NOT. RSIDE) .AND. J.GE.I ) 
     $              .OR. (.NOT. CSYMM) ) THEN
C     
                  GSI = IWORK( IGSI + ABS(IIS - I) ) 
                  GSJ = IWORK( IGSJ + ABS(IIS - I) )
C
                  CALL INFOG2L( GSI, GSJ, DESCSC, NPROW, NPCOL, 
     $                 MYROW, MYCOL, IX, JX, SRSRC, SCSRC )
                  SNODE = MYROW.EQ.SRSRC.AND.MYCOL.EQ.SCSRC
C     
                  IF( SNODE ) THEN
                     IF( NPCOL.GT.NPROW ) THEN
                        CALL DGEBS2D( ICTXT, 'Row', ' ', 1, 1, SCALOC, 
     $                       1 )
                     ELSEIF( NPROW.GT.1 ) THEN
                        CALL DGEBS2D( ICTXT, 'Col', ' ', 1, 1, SCALOC, 
     $                       1 )
                     END IF
                     LSCALC = SCALOC
                  ELSE
                     IF( NPCOL.GT.NPROW ) THEN
                        IF( MYROW.EQ.SRSRC ) THEN
                           CALL DGEBR2D( ICTXT, 'Row', ' ', 1, 1, 
     $                          SCALOC, 1, SRSRC, SCSRC )
                        END IF
                     ELSEIF( NPROW.GT.1 ) THEN
                        IF( MYCOL.EQ.SCSRC ) THEN
                           CALL DGEBR2D( ICTXT, 'Col', ' ', 1, 1, 
     $                          SCALOC, 1, SRSRC, SCSRC )
                        END IF
                     END IF
                  END IF
C     
                  IF( .NOT. RSIDE ) THEN
                     J = J - 1
                  ELSEIF( RSIDE ) THEN
                     J = J + 1
                  END IF
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
                  CALL PDSCAL( N, SCALOC, C, IC, JC+KKK-1, DESCC, 1 )
                  CALL PDSCAL( N, SCALOC, DWORK( E ), 1, KKK, DESCE, 
     $                 1 )
 123           CONTINUE
               IF( AEXT ) THEN
                  CALL DSCAL( EXMEMC, SCALOC, DWORK( EXC ), 1 )
                  CALL DSCAL( EXMEME, SCALOC, DWORK( EXE ), 1 )
               END IF
C
C     Local subsystem scaling
C
               J = JJS
               DO 222 I = IIS, IIE, STEP
C     
                  IF(  (CSYMM .AND. RSIDE .AND. J.LE.I ) .OR. 
     $                 (CSYMM. AND. (.NOT. RSIDE) .AND. J.GE.I ) 
     $                 .OR. (.NOT. CSYMM) ) THEN
C
                     ROWS = IWORK( IROWS + ABS(IIS - I) )
                     GSI = IWORK( IGSI + ABS(IIS - I) ) 
                     COLS = IWORK( ICOLS + ABS(IIS - I) ) 
                     GSJ = IWORK( IGSJ + ABS(IIS - I) )
C     
                     IF( ROWS.EQ.0 .OR. COLS.EQ.0 ) GO TO 777
C     
                     CALL INFOG2L( GSI, GSJ, DESCSC, NPROW, NPCOL, 
     $                    MYROW, MYCOL, IX, JX, SRSRC, SCSRC )     
                     SNODE = MYROW.EQ.SRSRC.AND.MYCOL.EQ.SCSRC
C     
                     IF( SNODE ) THEN 
                        IF( SCALOC.NE.LSCALC ) THEN
                           DO 333 KKK = 1, COLS
                              CALL DSCAL( ROWS, SCALOC / LSCALC, 
     $                             DWORK( XIJ + (KKK-1)*(NB+1) ), 1 )
 333                       CONTINUE
                        END IF
                        DO 555 KKK = 1, COLS
                           CALL DSCAL( ROWS, SCALOC, 
     $                          DWORK( CIJ + (KKK-1)*(NB+1) ), 1 )
 555                    CONTINUE
                        CALL DLACPY( 'All', ROWS, COLS, 
     $                       DWORK( XIJ ), NB+1, DWORK( CIJ ), NB+1 )
                        EXROW = INT2LG(IWORK( IEXRW + ABS(IIS - I) )) 
                        EXCOL = INT2LG(IWORK( IEXCL + ABS(IIS - I) ))
                        CALL DUBEXMA( EXROW, EXCOL, ROWS, COLS, IX, 
     $                       JX, C((LJCS-1)*LLDC+LICS), LLDC, NBCA, 
     $                       NB, NB, DWORK( EXC ), DWORK( CIJ ), 
     $                       NB + 1 )
                     END IF
                  END IF
C     
 777              CONTINUE
C     
                  IF( .NOT. RSIDE ) THEN
                     J = J - 1
                  ELSEIF( RSIDE ) THEN
                     J = J + 1
                  END IF
C     
 222           CONTINUE
C     
C     Update value of SCALE according to SCALOC
C     
               SCALE = SCALE * SCALOC
               SCALOC = ONE
            END IF
C     
C     In the case of symmetry (and we are not located on the main block
C     diagonal), we possibly send the current subsolution
C     to the process holding Cji, if that node hasn't received it yet
C     
            IF( CSYMM ) THEN
               J = JJS
               DO 36 I = IIS, IIE, STEP
                  IF( ( RSIDE .AND. J.LE.I ) .OR. 
     $                 (.NOT. RSIDE .AND. J.GE.I ) ) THEN 
C     
                     ROWS = IWORK( IROWS + ABS(IIS - I) )
                     GSI = IWORK( IGSI + ABS(IIS - I) ) 
                     EXROW = INT2LG(IWORK( IEXRW + ABS(IIS - I) )) 
                     COLS = IWORK( ICOLS + ABS(IIS - I) ) 
                     GSJ = IWORK( IGSJ + ABS(IIS - I) )
                     EXCOL = INT2LG(IWORK( IEXCL + ABS(IIS - I) ))
C     
                     IF( ROWS.EQ.0 .OR. COLS.EQ.0 ) GO TO 37
C     
                     CALL INFOG2L( GSI, GSJ, DESCSC, NPROW, NPCOL, 
     $                    MYROW, MYCOL, IX, JX, SRSRC, SCSRC )
                     CALL INFOG2L( GSJ, GSI, DESCSC, NPROW, NPCOL, 
     $                    MYROW, MYCOL, TIX, TJX, CTRSRC, 
     $                    CTCSRC )
C
                     IF( CTRSRC.NE.SRSRC .OR. CTCSRC.NE.SCSRC ) THEN
                        IF( MYROW.EQ.SRSRC .AND. MYCOL.EQ.SCSRC ) THEN
                           CALL DGESD2D( ICTXT, ROWS, COLS, 
     $                          DWORK( XIJ ), NB + 1, CTRSRC, CTCSRC )
                        END IF
                        IF(MYROW.EQ.CTRSRC.AND.MYCOL.EQ.CTCSRC) THEN
                           CALL DGERV2D( ICTXT, ROWS, COLS, 
     $                          DWORK( IPW ), ROWS, SRSRC, SCSRC )
                        END IF   
                     END IF
C     
                     IF( MYROW.EQ.CTRSRC.AND.MYCOL.EQ.CTCSRC ) THEN
                        IF(CTRSRC.NE.SRSRC.OR.CTCSRC.NE.SCSRC) THEN
                           CALL DLATCPY( 'All', ROWS, COLS, 
     $                          DWORK( IPW ), ROWS, DWORK( XIJT ), 
     $                          NB + 1 )
                        ELSE
                           CALL DLATCPY( 'All', ROWS, COLS, 
     $                          DWORK( XIJ ), NB + 1, DWORK( XIJT ), 
     $                          NB + 1 )
                        END IF
C     
C     Put Xij^T back into the local part of the global matrix
C     
                        CALL DUBEXMA( EXCOL, EXROW, COLS, ROWS, 
     $                       TIX, TJX, C((LJCS-1)*LLDC+LICS), LLDC, 
     $                       NBCA, NB, NB, DWORK( EXC ), DWORK( XIJT ), 
     $                       NB + 1 )
                     END IF
                  END IF
 37               CONTINUE
C     
C     Update inner loop variable j (column index)
C     
                  IF( .NOT. RSIDE ) THEN
                     J = J - 1
                  ELSEIF( RSIDE ) THEN
                     J = J + 1
                  END IF
 36            CONTINUE
            END IF
C     
C     Broadcast the local solution Xij in block row i
C     
            IF ( K.LT.NROLL ) THEN
               J = JJS
               DO 40 I = IIS, IIE, STEP
C     
                  ROWS = IWORK( IROWS + ABS(IIS - I) )
                  GSI = IWORK( IGSI + ABS(IIS - I) ) 
                  EXROW = INT2LG(IWORK( IEXRW + ABS(IIS - I) )) 
                  COLS = IWORK( ICOLS + ABS(IIS - I) ) 
                  GSJ = IWORK( IGSJ + ABS(IIS - I) )
                  EXCOL = INT2LG(IWORK( IEXCL + ABS(IIS - I) ))
C     
                  IF( ROWS.EQ.0 .OR. COLS.EQ.0 ) GOTO 45
C     
                  CALL INFOG2L( GSI, GSJ, DESCSC, NPROW, NPCOL, 
     $                          MYROW, MYCOL, IX, JX, CRSRC, CCSRC )
                  SNODE = MYROW.EQ.CRSRC.AND.MYCOL.EQ.CCSRC
                  SRSRC = CRSRC
                  SCSRC = CCSRC
C
C     Build the matrix Xij 
C
                  IF( SNODE ) THEN
                     CALL DBEXMAT( EXROW, EXCOL, ROWS, COLS, IX, JX, 
     $                    C((LJCS-1)*LLDC+LICS), LLDC, NBCA, NB, NB, 
     $                    DWORK( EXC ), DWORK( XIJ ), NB + 1 )
                  END IF
C     
C     Depending on who had the block, PERHAPS broadcast it
C     
                  IF( NPCOL.GT.1 ) THEN
                     IF ( SNODE ) THEN
                        CALL DGEBS2D( ICTXT, 'Row', ' ', ROWS, COLS, 
     $                       DWORK( XIJ ), NB + 1 )
                     ELSEIF( MYROW.EQ.SRSRC ) THEN
                        CALL DGEBR2D( ICTXT, 'Row', ' ', ROWS, COLS, 
     $                       DWORK( XIJ ), NB + 1, SRSRC, SCSRC )
                     END IF
                  END IF 
 45               CONTINUE
                  IF( .NOT. RSIDE ) THEN
                     J = J - 1
                  ELSEIF( RSIDE ) THEN
                     J = J + 1
                  END IF
 40            CONTINUE
            END IF
C     
C     Update rest of global system wrt to current solution - first
C     compute intermediate values Eik = Eik + Xij*Akj^T or
C     Eik = Eik + Xij*Ajk . Choose the right branch depending on 
C     the value of OP (RSIDE) in the equation.
C     
C     For every update we do we also need to find out the dimensions and
C     possible extensions of the submatrices involved.
C     
            IF( .NOT. RSIDE ) THEN
               J = JJS
               DO 50 I = IIS, IIE, STEP
                  ROWS = IWORK( IROWS + ABS(IIS - I) )
                  GSI = IWORK( IGSI + ABS(IIS - I) ) 
                  EXROW = INT2LG(IWORK( IEXRW + ABS(IIS - I) )) 
                  COLS = IWORK( ICOLS + ABS(IIS - I) ) 
                  GSJ = IWORK( IGSJ + ABS(IIS - I) )
                  EXCOL = INT2LG(IWORK( IEXCL + ABS(IIS - I) ))
C     
                  IF( ROWS.EQ.0 .OR. COLS.EQ.0 ) GO TO 55
C     
                  IF( MOD( I, 2 ).EQ.0 ) THEN
                     IF( CSYMM .AND. I.NE.DBA ) THEN
                        IF( I.EQ.J ) THEN
                           INDXS = MAX( I, J + 1 )
                        ELSE
                           INDXS = MAX( I, J )
                        END IF 
                     ELSEIF( CSYMM .AND. I.EQ.DBA ) THEN
                        INDXS = MAX( I, J + 1 )
                     ELSEIF( (.NOT. CSYMM) .AND. I.EQ.DBA ) THEN
                        INDXS = J + 1
                     ELSEIF( (.NOT. CSYMM) .AND. I.NE.DBA ) THEN
                        INDXS = J
                     END IF
                     INDXE = DBA
                     INDXU = 1
                  ELSE
                     INDXS = DBA
                     IF( CSYMM .AND. I.NE.DBA ) THEN
                        IF( I.EQ.J ) THEN
                           INDXE = MAX( I, J + 1 )
                        ELSE
                           INDXE = MAX( I, J )
                        END IF 
                     ELSEIF( CSYMM .AND. I.EQ.DBA ) THEN
                        INDXE = MAX( I, J + 1 )
                     ELSEIF( (.NOT. CSYMM) .AND. I.EQ.DBA ) THEN
                        INDXE = J + 1
                     ELSEIF( (.NOT. CSYMM) .AND. I.NE.DBA ) THEN
                        INDXE = J
                     END IF
                     INDXU = -1
                  END IF
                  DO 60 INDX = INDXS, INDXE, INDXU 
                     IF( INDX.EQ.1 ) THEN
                        IF( IWORK(EXAINF).EQ.0 ) THEN
                           ACOLS2 = NB - IROFFA
                           GSIND = 1 + IROFFA
                           EEXCOL = .FALSE.
                        ELSEIF( IWORK(EXAINF).EQ.1 ) THEN
                           ACOLS2 = NB - IROFFA + 1
                           GSIND = 1 + IROFFA
                           EEXCOL = .TRUE.
                        END IF
                     ELSE
                        IF( IWORK(EXAINF+(INDX-1)).EQ.0 ) THEN
                           ACOLS2 = MIN(NB, (N - NB + IROFFA) - 
     $                          (INDX - 2) * NB)
                           GSIND = (INDX - 1) * NB + 1
                           EEXCOL = .FALSE.
                        ELSEIF( IWORK(EXAINF+(INDX-1)).EQ.1 ) THEN
                           ACOLS2 = NB + 1
                           GSIND = (INDX - 1) * NB + 1
                           EEXCOL = .TRUE.
                        ELSEIF( IWORK(EXAINF+(INDX-1)).EQ.2 ) THEN
                           ACOLS2 = MIN(NB, (N - NB + IROFFA) - 
     $                          (INDX - 2) * NB) - 1
                           GSIND = (INDX - 1) * NB + 2
                           EEXCOL = .FALSE.
                        ELSEIF( IWORK(EXAINF+(INDX-1)).EQ.3 ) THEN
                           ACOLS2 = NB 
                           GSIND = (INDX - 1) * NB + 2
                           EEXCOL = .TRUE.
                        END IF
                     END IF
C     
C     Build and, if needed, communicate the submatrices 
C     
                     CALL INFOG2L( GSJ, GSIND, DESCSA, NPROW, NPCOL,
     $                             MYROW, MYCOL, LIA, LJA, RSRC, CSRC )
                     CALL INFOG2L( GSI, GSIND, DESCE, NPROW, NPCOL, 
     $                             MYROW, MYCOL, IX, JX, RRSRC, RCSRC )
                     IF( MYROW.EQ.RSRC .AND. MYCOL.EQ.CSRC ) THEN
                        CALL DBEXMAT( EXCOL, EEXCOL, COLS, ACOLS2, LIA, 
     $                                LJA, A((LJAS-1)*LLDA+LIAS), LLDA, 
     $                                NBCA, NB, NB, DWORK( EXA ), 
     $                                DWORK( AJJ ), NB + 1 )
                        IF( (RSRC.NE.RRSRC) .OR. (CSRC.NE.RCSRC) ) THEN
                           CALL DGESD2D( ICTXT, COLS, ACOLS2, 
     $                                   DWORK( AJJ ), NB + 1, RRSRC, 
     $                                   RCSRC )
                        END IF
                     END IF
                     IF( MYROW.EQ.RRSRC .AND. MYCOL.EQ.RCSRC ) THEN
C     
C     If we are on the beginning of a new line set Eik=[0] else
C     build the matrix E from previously saved results
C     
                        IF( J.EQ.1 ) THEN
                           CALL DLASET( 'All', NB + 1, NB + 1, ZERO, 
     $                                  ZERO, DWORK(EIJ), NB + 1 )
                        ELSE
                           CALL DBEXMAT( EXROW, EEXCOL, ROWS, ACOLS2, 
     $                          IX, JX, DWORK(E), LLDE, NBCA, NB, NB, 
     $                          DWORK( EXE ), DWORK(EIJ), NB + 1 ) 
                        END IF
                        IF( (RSRC.NE.RRSRC) .OR. (CSRC.NE.RCSRC) ) THEN
                           CALL DGERV2D( ICTXT, COLS, ACOLS2, 
     $                          DWORK( AJJ ), NB + 1, RSRC, CSRC )
                        END IF
C     
C     Perform the update of Eik
C     
                        CALL DGEMM( 'N','N', ROWS, ACOLS2, COLS, ONE,
     $                       DWORK( XIJ ), NB + 1, DWORK( AJJ ),
     $                       NB + 1, ONE, DWORK(EIJ), NB + 1 )
C     
C     Save the result in the matrix E
C     
                        CALL DUBEXMA( EXROW, EEXCOL, ROWS, ACOLS2, IX, 
     $                       JX, DWORK(E), LLDE, NBCA, NB, NB, 
     $                       DWORK( EXE ), DWORK(EIJ), NB + 1 )
                     END IF                     
 60               CONTINUE
C     
 55               CONTINUE
C     
                  J = J - 1
 50            CONTINUE
C     
            ELSEIF( RSIDE ) THEN
               J = JJS
               DO 70 I = IIS, IIE, STEP
                  ROWS = IWORK( IROWS + ABS(IIS - I) )
                  GSI = IWORK( IGSI + ABS(IIS - I) ) 
                  EXROW = INT2LG(IWORK( IEXRW + ABS(IIS - I) )) 
                  COLS = IWORK( ICOLS + ABS(IIS - I) ) 
                  GSJ = IWORK( IGSJ + ABS(IIS - I) )
                  EXCOL = INT2LG(IWORK( IEXCL + ABS(IIS - I) ))
C     
                  IF( ROWS.EQ.0 .OR. COLS.EQ.0 ) GO TO 75
C     
                  IF( MOD( I, 2 ).EQ.0 ) THEN
                     INDXS = 1
                     IF( CSYMM .AND. I.NE.1 ) THEN
                        IF( I.EQ.J ) THEN
                           INDXE = MIN( I, J - 1)
                        ELSE
                           INDXE = MIN( I, J )
                        END IF 
                     ELSEIF( CSYMM .AND. I.EQ.1 ) THEN
                        INDXE = MIN( I, J - 1 )
                     ELSEIF( (.NOT. CSYMM) .AND. I.EQ.1 ) THEN
                        INDXE = J - 1
                     ELSEIF( (.NOT. CSYMM) .AND. I.NE.1 ) THEN
                        INDXE = J
                     END IF
                     INDXU = 1
                  ELSE
                     IF( CSYMM .AND. I.NE.1 ) THEN
                        IF( I.EQ.J ) THEN
                           INDXS = MIN( I, J - 1)
                        ELSE
                           INDXS = MIN( I, J )
                        END IF 
                     ELSEIF( CSYMM .AND. I.EQ.1 ) THEN
                        INDXS = MIN( I, J - 1 )
                     ELSEIF( (.NOT. CSYMM) .AND. I.EQ.1 ) THEN
                        INDXS = J - 1
                     ELSEIF( (.NOT. CSYMM) .AND. I.NE.1 ) THEN
                        INDXS = J
                     END IF
                     INDXE = 1
                     INDXU = -1
                  END IF
                  DO 80 INDX =  INDXS, INDXE, INDXU
C     
C     Set some constants describing the involved submatrices
C     
                     IF( INDX.EQ.1 ) THEN
                        IF( IWORK(EXAINF).EQ.0 ) THEN
                           AROWS2 = NB - IROFFA
                           GSIND = 1 + IROFFA
                           EEXCOL = .FALSE.
                        ELSEIF( IWORK(EXAINF).EQ.1 ) THEN
                           AROWS2 = NB - IROFFA + 1
                           GSIND = 1 + IROFFA
                           EEXCOL = .TRUE.
                        END IF
                     ELSE
                        IF( IWORK(EXAINF+(INDX-1)).EQ.0 ) THEN
                           AROWS2 = MIN(NB, (N - NB + IROFFA) - 
     $                          (INDX - 2) * NB)
                           GSIND = (INDX - 1) * NB + 1
                           EEXCOL = .FALSE.
                        ELSEIF( IWORK(EXAINF+(INDX-1)).EQ.1 ) THEN
                           AROWS2 = NB + 1
                           GSIND = (INDX - 1) * NB + 1
                           EEXCOL = .TRUE.
                        ELSEIF( IWORK(EXAINF+(INDX-1)).EQ.2 ) THEN
                           AROWS2 = MIN(NB, (N - NB + IROFFA) - 
     $                          (INDX - 2) * NB) - 1
                           GSIND = (INDX - 1) * NB + 2
                           EEXCOL = .FALSE.
                        ELSEIF( IWORK(EXAINF+(INDX-1)).EQ.3 ) THEN
                           AROWS2 = NB 
                           GSIND = (INDX - 1) * NB + 2
                           EEXCOL = .TRUE.
                        END IF
                     END IF
C     
C     Build and, if needed, communicate the submatrices 
C     
                     CALL INFOG2L( GSIND, GSJ, DESCSA, NPROW, NPCOL,
     $                    MYROW, MYCOL, LIA, LJA, RSRC, CSRC )
                     CALL INFOG2L( GSI, GSIND, DESCE, NPROW, NPCOL, 
     $                    MYROW, MYCOL, IX, JX, RRSRC, RCSRC )
                     IF( MYROW.EQ.RSRC .AND. MYCOL.EQ.CSRC ) THEN
                        CALL DBEXMAT( EEXCOL, EXCOL, AROWS2, COLS, LIA, 
     $                       LJA, A((LJAS-1)*LLDA+LIAS), LLDA, 
     $                       NBCA, NB, NB, DWORK( EXA ), 
     $                       DWORK( AJJ ), NB + 1 )
                        IF( (RSRC.NE.RRSRC) .OR. (CSRC.NE.RCSRC) ) THEN
                           CALL DGESD2D( ICTXT, AROWS2, COLS, 
     $                          DWORK( AJJ ), NB + 1, RRSRC, RCSRC )
                        END IF
                     END IF
                     IF(MYROW.EQ.RRSRC .AND. MYCOL.EQ.RCSRC) THEN
C     
C     If we are on the beginning of a new line set Eik=[0] else
C     build the matrix E from previously saved results
C     
                        IF( J.EQ.DBA ) THEN
                           CALL DLASET( 'All', NB + 1, NB + 1, ZERO, 
     $                          ZERO, DWORK(EIJ), NB + 1 )
                        ELSE
                           CALL DBEXMAT( EXROW, EEXCOL, ROWS, AROWS2, 
     $                          IX, JX, DWORK(E), LLDE, NBCA, NB, NB, 
     $                          DWORK( EXE ), DWORK(EIJ), NB + 1 ) 
                        END IF
                        IF( (RSRC.NE.RRSRC) .OR. (CSRC.NE.RCSRC) ) THEN
                           CALL DGERV2D( ICTXT, AROWS2, COLS, 
     $                          DWORK( AJJ ), NB + 1, RSRC, CSRC )
                        END IF
C     
C     Perform the update of Eik
C     
                        CALL DGEMM( 'N','T', ROWS, AROWS2, COLS, ONE,
     $                       DWORK( XIJ ), NB + 1, DWORK( AJJ ),
     $                       NB + 1, ONE, DWORK(EIJ), NB + 1 )
C     
C     Save the result in the matrix E
C     
                        CALL DUBEXMA( EXROW, EEXCOL, ROWS, AROWS2, IX, 
     $                       JX, DWORK(E), LLDE, NBCA, NB, NB, 
     $                       DWORK( EXE ), DWORK(EIJ), NB + 1 )
                     END IF                     
 80               CONTINUE
C     
 75               CONTINUE
C     
                  J = J + 1
 70            CONTINUE
            END IF
C     
C     Broadcast the matrix Eij in the corresponding block column 
C     
            IF( K.LT.NROLL .AND. .NOT. RSIDE ) THEN
               J = JJS
               DO 90 I = IIS, IIE, STEP
                  IF( (CSYMM .AND. I.LT.J) .OR. (.NOT. CSYMM) ) THEN
                     ROWS = IWORK( IROWS + ABS(IIS - I) )
                     GSI = IWORK( IGSI + ABS(IIS - I) ) 
                     EXROW = INT2LG(IWORK( IEXRW + ABS(IIS - I) )) 
                     COLS = IWORK( ICOLS + ABS(IIS - I) ) 
                     GSJ = IWORK( IGSJ + ABS(IIS - I) )
                     EXCOL = INT2LG(IWORK( IEXCL + ABS(IIS - I) ))
C     
                     IF( ROWS.EQ.0 .OR. COLS.EQ.0 ) GO TO 95
C     
                     CALL INFOG2L( GSI, GSJ, DESCE, NPROW, NPCOL, 
     $                             MYROW, MYCOL, IX, JX, CRSRC, CCSRC )
C     
                     SNODE = MYROW.EQ.CRSRC.AND.MYCOL.EQ.CCSRC
                     SRSRC = CRSRC
                     SCSRC = CCSRC
C     
                     IF( SNODE ) THEN
                        CALL DBEXMAT( EXROW, EXCOL, ROWS, COLS, IX, JX, 
     $                       DWORK(E), LLDE, NBCA, NB, NB, DWORK( EXE ), 
     $                       DWORK( EIJ ), NB + 1 )   
                        IF( NPROW.GT.1 ) THEN
                           CALL DGEBS2D( ICTXT, 'Col', ' ', ROWS, COLS, 
     $                                   DWORK(EIJ), NB + 1 ) 
                        END IF
                     ELSEIF( NPROW.GT.1 .AND. MYCOL.EQ.SCSRC ) THEN
                        CALL DGEBR2D( ICTXT, 'Col', ' ', ROWS, COLS,
     $                                DWORK(EIJ), NB + 1, SRSRC, SCSRC )
                     END IF
C     
 95                  CONTINUE
                  END IF
C     
                  J = J - 1
 90            CONTINUE
C     
            ELSEIF( K.LT.NROLL .AND. RSIDE ) THEN
               J = JJS
               DO 100 I = IIS, IIE, STEP
                  IF( (CSYMM .AND. I.GT.J) .OR. (.NOT. CSYMM) ) THEN
                     ROWS = IWORK( IROWS + ABS(IIS - I) )
                     GSI = IWORK( IGSI + ABS(IIS - I) ) 
                     EXROW = INT2LG(IWORK( IEXRW + ABS(IIS - I) )) 
                     COLS = IWORK( ICOLS + ABS(IIS - I) ) 
                     GSJ = IWORK( IGSJ + ABS(IIS - I) )
                     EXCOL = INT2LG(IWORK( IEXCL + ABS(IIS - I) ))
C     
                     IF( ROWS.EQ.0 .OR. COLS.EQ.0 ) GO TO 105
C     
                     CALL INFOG2L( GSI, GSJ, DESCE, NPROW, NPCOL, 
     $                    MYROW, MYCOL, IX, JX, CRSRC, CCSRC )
C     
                     SNODE = MYROW.EQ.CRSRC.AND.MYCOL.EQ.CCSRC
                     SRSRC = CRSRC
                     SCSRC = CCSRC
C     
                     IF( SNODE ) THEN
                        CALL DBEXMAT( EXROW, EXCOL, ROWS, COLS, IX, JX, 
     $                       DWORK(E), LLDE, NBCA, NB, NB, DWORK( EXE ), 
     $                       DWORK( EIJ ), NB + 1 )   
                        IF( NPROW.GT.1 ) THEN
                           CALL DGEBS2D( ICTXT, 'Col', ' ', ROWS, COLS, 
     $                                   DWORK(EIJ), NB + 1 ) 
                        END IF
                     ELSEIF( NPROW.GT.1 .AND. MYCOL.EQ.SCSRC ) THEN
                        CALL DGEBR2D( ICTXT, 'Col', ' ', ROWS, COLS,
     $                       DWORK(EIJ), NB + 1, SRSRC, SCSRC )
                     END IF
C     
 105                 CONTINUE
                  END IF
C     
                  J = J + 1
 100           CONTINUE
            END IF
C     
C     Go on with the updates:      
C     Now update block column j of C with respect to Eij in
C     the operation Ckj = Ckj - Aki*Eij or Ckj = Ckj - Aik^T*Eij 
C     
            J = JJS
            DO 110 I = IIS, IIE, STEP
C     
               ROWS = IWORK( IROWS + ABS(IIS - I) )
               GSI = IWORK( IGSI + ABS(IIS - I) ) 
               EXROW = INT2LG(IWORK( IEXRW + ABS(IIS - I) )) 
               COLS = IWORK( ICOLS + ABS(IIS - I) ) 
               GSJ = IWORK( IGSJ + ABS(IIS - I) )
               EXCOL = INT2LG(IWORK( IEXCL + ABS(IIS - I) ))
C     
               IF( ROWS.EQ.0 .OR. COLS.EQ.0 ) GO TO 115
C     
               IF( RSIDE .AND. ((CSYMM.AND.I.GT.J).OR.(.NOT.CSYMM) ) ) 
     $              THEN
                  IF( MOD( J, 2 ).EQ.0 ) THEN
                     INDXS = 1
                     IF( CSYMM ) INDXS = J
                     INDXE = I - 1
                     INDXU = 1
                  ELSE
                     INDXS = I - 1
                     INDXE = 1
                     IF( CSYMM ) INDXE = J
                     INDXU = -1
                  END IF
                  DO 120 INDX =  INDXS, INDXE, INDXU
C     
C     Set some constants describing the involved submatrices
C     
                     IF( INDX.EQ.1 ) THEN
                        IF( IWORK(EXAINF).EQ.0 ) THEN
                           AROWS2 = NB - IROFFA
                           GSIND = 1 + IROFFA
                           CEXROW = .FALSE.
                        ELSEIF( IWORK(EXAINF).EQ.1 ) THEN
                           AROWS2 = NB - IROFFA + 1
                           GSIND = 1 + IROFFA
                           CEXROW = .TRUE.
                        END IF
                     ELSE
                        IF( IWORK(EXAINF+(INDX-1)).EQ.0 ) THEN
                           AROWS2 = MIN(NB, (N - NB + IROFFA) - 
     $                          (INDX - 2) * NB)
                           GSIND = (INDX - 1) * NB + 1
                           CEXROW = .FALSE.
                        ELSEIF( IWORK(EXAINF+(INDX-1)).EQ.1 ) THEN
                           AROWS2 = NB + 1
                           GSIND = (INDX - 1) * NB + 1
                           CEXROW = .TRUE.
                        ELSEIF( IWORK(EXAINF+(INDX-1)).EQ.2 ) THEN
                           AROWS2 = MIN(NB, (N - NB + IROFFA) - 
     $                          (INDX - 2) * NB) - 1
                           GSIND = (INDX - 1) * NB + 2
                           CEXROW = .FALSE.
                        ELSEIF( IWORK(EXAINF+(INDX-1)).EQ.3 ) THEN
                           AROWS2 = NB 
                           GSIND = (INDX - 1) * NB + 2
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
                     IF( MYROW.EQ.RSRC .AND. MYCOL.EQ.CSRC ) THEN
                        CALL DBEXMAT( CEXROW, EXROW, AROWS2, ROWS, LIA, 
     $                                LJA, A((LJAS-1)*LLDA+LIAS), LLDA, 
     $                                NBCA, NB, NB, DWORK( EXA ), 
     $                                DWORK( AII ), NB + 1 )
                        IF( (RSRC.NE.RRSRC) .OR. (CSRC.NE.RCSRC) ) THEN
                           CALL DGESD2D( ICTXT, AROWS2, ROWS, 
     $                                   DWORK( AII ), NB + 1, RRSRC, 
     $                                   RCSRC )
                        END IF
                     END IF
                     IF(MYROW.EQ.RRSRC .AND. MYCOL.EQ.RCSRC) THEN  
                        CALL DBEXMAT( CEXROW, EXCOL, AROWS2, COLS, IX, 
     $                       JX, C((LJCS-1)*LLDC+LICS), LLDC, NBCA, NB, 
     $                       NB, DWORK( EXC ), DWORK(CIJ), NB + 1 ) 
                        IF( (RSRC.NE.RRSRC) .OR. (CSRC.NE.RCSRC) ) THEN
                           CALL DGERV2D( ICTXT, AROWS2, ROWS, 
     $                          DWORK( AII ), NB + 1, RSRC, CSRC )
                        END IF     
C     
C     Perform the update of Ckj
C     
                        CALL DGEMM( 'N','N',AROWS2, COLS, ROWS, -ONE,
     $                       DWORK( AII ), NB + 1, DWORK( EIJ ),
     $                       NB + 1, ONE, DWORK(CIJ), NB + 1 )
C     
C     Save the result in the matrix C
C     
                        CALL DUBEXMA( CEXROW, EXCOL, AROWS2, COLS, IX, 
     $                       JX, C((LJCS-1)*LLDC+LICS), LLDC, NBCA, 
     $                       NB, NB, DWORK( EXC ), DWORK(CIJ), NB + 1 ) 
                     END IF
 120              CONTINUE
C     
               ELSEIF( .NOT. RSIDE .AND. 
     $                 ((CSYMM.AND.I.LT.J).OR.(.NOT.CSYMM) )) THEN
                  IF( MOD( J, 2 ).EQ.0 ) THEN
                     INDXS = I + 1
                     INDXE = DBA
                     IF( CSYMM ) INDXE = J
                     INDXU = 1
                  ELSE
                     INDXS = DBA
                     IF( CSYMM ) INDXS = J
                     INDXE = I + 1
                     INDXU = -1
                  END IF
                  DO 130 INDX =  INDXS, INDXE, INDXU
C     
C     Set some constants descibing the involved submatrices
C     
                     IF( IWORK(EXAINF+(INDX-1)).EQ.0 ) THEN
                        ACOLS2 = MIN(NB, (N - NB + IROFFA) - 
     $                       (INDX - 2) * NB)
                        GSIND = (INDX - 1) * NB + 1
                        CEXROW = .FALSE.
                     ELSEIF( IWORK(EXAINF+(INDX-1)).EQ.1 ) THEN
                        ACOLS2 = NB + 1
                        GSIND = (INDX - 1) * NB + 1
                        CEXROW = .TRUE.
                     ELSEIF( IWORK(EXAINF+(INDX-1)).EQ.2 ) THEN
                        ACOLS2 = MIN(NB, (N - NB + IROFFA) - 
     $                       (INDX - 2) * NB) - 1
                        GSIND = (INDX - 1) * NB + 2
                        CEXROW = .FALSE.
                     ELSEIF( IWORK(EXAINF+(INDX-1)).EQ.3 ) THEN
                        ACOLS2 = NB 
                        GSIND = (INDX - 1) * NB + 2
                        CEXROW = .TRUE.
                     END IF
C     
C     Build and, if needed, communicate the submatrices 
C     
                     CALL INFOG2L( GSI, GSIND, DESCSA, NPROW, NPCOL,
     $                             MYROW, MYCOL, LIA, LJA, RSRC, CSRC )
                     CALL INFOG2L( GSIND, GSJ, DESCSC, NPROW, NPCOL, 
     $                             MYROW, MYCOL, IX, JX, RRSRC, RCSRC )
                     IF( MYROW.EQ.RSRC .AND. MYCOL.EQ.CSRC ) THEN
                        CALL DBEXMAT( EXROW, CEXROW, ROWS, ACOLS2, LIA, 
     $                       LJA, A((LJAS-1)*LLDA+LIAS), LLDA, NBCA, NB, 
     $                       NB, DWORK( EXA ), DWORK( AII ), NB + 1 )
                        IF( (RSRC.NE.RRSRC) .OR. (CSRC.NE.RCSRC) ) THEN
                           CALL DGESD2D( ICTXT, ROWS, ACOLS2, 
     $                          DWORK( AII ), NB + 1, RRSRC, RCSRC )
                        END IF
                     END IF
                     IF(MYROW.EQ.RRSRC .AND. MYCOL.EQ.RCSRC) THEN
                        CALL DBEXMAT( CEXROW, EXCOL, ACOLS2, COLS, 
     $                       IX, JX, C((LJCS-1)*LLDC+LICS), LLDC, 
     $                       NBCA, NB, NB, DWORK( EXC ), DWORK(CIJ), 
     $                       NB + 1 )
                        IF( (RSRC.NE.RRSRC).OR.(CSRC.NE.RCSRC) ) THEN
                           CALL DGERV2D( ICTXT, ROWS, ACOLS2, 
     $                          DWORK( AII ), NB + 1, RSRC, CSRC )
                        END IF 
C     
C     Perform the update of Ckj
C     
                        CALL DGEMM( 'T','N', ACOLS2, COLS, ROWS, 
     $                       -ONE, DWORK( AII ), NB + 1, 
     $                       DWORK( EIJ ), NB + 1, ONE, DWORK(CIJ), 
     $                       NB + 1 )
C     
C     Save the result in the matrix C
C     
                        CALL DUBEXMA( CEXROW, EXCOL, ACOLS2, COLS, 
     $                       IX, JX, C((LJCS-1)*LLDC+LICS), LLDC, 
     $                       NBCA, NB, NB, DWORK( EXC ), DWORK(CIJ), 
     $                       NB + 1 ) 
                     END IF
 130              CONTINUE
               END IF
C     
 115           CONTINUE
C     
C     Update inner loop variable j (column index)
C     
               IF( .NOT. RSIDE ) THEN
                  J = J - 1
               ELSEIF( RSIDE ) THEN
                  J = J + 1
               END IF
C    
 110        CONTINUE
C     
C     Prepare for solving for the next block diagonal of C:
C     update Ci,j-1 = Ci,j-1 - Aii * Ei,j-1 or
C     do update  Ci,j+1 = Ci,j+1 - Aii^T * Ei,j+1
C     
            J = JJS
            DO 140 I = IIS, IIE, STEP
C     
               ROWS = IWORK( IROWS + ABS(IIS - I) )
               GSI = IWORK( IGSI + ABS(IIS - I) ) 
               EXROW = INT2LG(IWORK( IEXRW + ABS(IIS - I) )) 
               COLS = IWORK( ICOLS + ABS(IIS - I) ) 
               GSJ = IWORK( IGSJ + ABS(IIS - I) )
               EXCOL = INT2LG(IWORK( IEXCL + ABS(IIS - I) ))
C     
               IF( ROWS.EQ.0 .OR. COLS.EQ.0 ) GO TO 145
C     
               IF( RSIDE ) THEN
                  IF( J.GT.1 .AND. ( (CSYMM .AND. J.LE.( I+1 )) .OR. 
     $                 ( .NOT. CSYMM ) ) ) THEN
C     
C     Build all submatrices and, if needed, communicate the submatrix 
C     Aii involved.      
C     Check if Aj-1,j-1 is extended and set some variables describing 
C     Ci,j-1 and Ei,j-1
C     
                     IF( J.EQ.2 ) THEN
                        IF( IWORK(EXAINF).EQ.0 ) THEN
                           CCOLS = NB - IROFFA
                           GSIND = 1 + IROFFA
                           CEXCOL = .FALSE.
                        ELSEIF( IWORK(EXAINF).EQ.1 ) THEN
                           CCOLS = NB - IROFFA + 1
                           GSIND = 1 + IROFFA
                           CEXCOL = .TRUE.
                        END IF
                     ELSE
                        IF( IWORK(EXAINF+(J-2)).EQ.0 ) THEN
                           CCOLS = MIN(NB, (N - NB + IROFFA) - 
     $                          (J - 3) * NB)
                           GSIND = (J - 2) * NB + 1
                           CEXCOL = .FALSE.
                        ELSEIF( IWORK(EXAINF+(J-2)).EQ.1 ) THEN
                           CCOLS = NB + 1
                           GSIND = (J - 2) * NB + 1
                           CEXCOL = .TRUE.
                        ELSEIF( IWORK(EXAINF+(J-2)).EQ.2 ) THEN
                           CCOLS = MIN(NB, (N - NB + IROFFA) - 
     $                          (J - 3) * NB) - 1
                           GSIND = (J - 2) * NB + 2
                           CEXCOL = .FALSE.
                        ELSEIF( IWORK(EXAINF+(J-2)).EQ.3 ) THEN
                           CCOLS = NB 
                           GSIND = (J - 2) * NB + 2
                           CEXCOL = .TRUE.
                        END IF
                     END IF
C     
                     CALL INFOG2L( GSI, GSI, DESCSA, NPROW, NPCOL, 
     $                             MYROW, MYCOL, LIA, LJA, ARSRC, 
     $                             ACSRC )
                     CALL INFOG2L( GSI, GSIND, DESCSC, NPROW, NPCOL, 
     $                             MYROW, MYCOL, LIC, LJC, CRSRC, 
     $                             CCSRC )
                     CALL INFOG2L( GSI, GSIND, DESCE, NPROW, NPCOL, 
     $                             MYROW, MYCOL, LEI, LEJ, CRSRC, 
     $                             CCSRC )
C     
                     IF( MYROW.EQ.ARSRC .AND. MYCOL.EQ.ACSRC ) THEN
                        CALL DBEXMAT( EXROW, EXROW, ROWS, ROWS, LIA, 
     $                       LJA, A((LJAS-1)*LLDA+LIAS), LLDA, NBCA, 
     $                       NB, NB, DWORK( EXA ), DWORK( AII ), 
     $                       NB + 1 )
                        IF( CRSRC.NE.ARSRC .OR. CCSRC.NE.ACSRC ) THEN
                           CALL DGESD2D( ICTXT, ROWS, ROWS, 
     $                          DWORK( AII ), NB + 1, CRSRC, CCSRC )
                        END IF
                     END IF
C     
                     IF( MYROW.EQ.CRSRC .AND. MYCOL.EQ.CCSRC ) THEN    
                        CALL DBEXMAT( EXROW, CEXCOL, ROWS, CCOLS, LIC, 
     $                       LJC, C((LJCS-1)*LLDC+LICS), LLDC, NBCA, 
     $                       NB, NB, DWORK( EXC ), DWORK( CIJ ), 
     $                       NB + 1 )
                        CALL DBEXMAT( EXROW, CEXCOL, ROWS, CCOLS, LEI, 
     $                       LEJ, DWORK( E ), LLDE, NBCA, NB, NB, 
     $                       DWORK( EXE ), DWORK( EIJ ), NB + 1 )
                        IF( CRSRC.NE.ARSRC .OR. CCSRC.NE.ACSRC ) THEN
                           CALL DGERV2D( ICTXT, ROWS, ROWS, 
     $                          DWORK( AII ), NB + 1, ARSRC, ACSRC )
                        END IF
C     
C     Do update
C     
                        CALL DGEMM( 'N', 'N', ROWS, CCOLS, ROWS, -ONE,
     $                       DWORK( AII ), NB + 1, DWORK( EIJ ), NB + 1, 
     $                       ONE, DWORK(CIJ), NB + 1 )
C     
C     Save the result in the matrix C
C     
                        CALL DUBEXMA( EXROW, CEXCOL, ROWS, CCOLS, LIC, 
     $                       LJC, C((LJCS-1)*LLDC+LICS), LLDC, NBCA, 
     $                       NB, NB, DWORK( EXC ), DWORK(CIJ), NB + 1 ) 
                     END IF
                  END IF
C     
               ELSEIF( .NOT. RSIDE ) THEN
                  IF( J.LT.DBA .AND. ( ( CSYMM .AND. J.GE.(I - 1) .OR.
     $                 ( .NOT.CSYMM ) ))) THEN
C     
C     Build all submatrices and, if needed, communicate the submatrix 
C     Aii involved.      
C     Check if Aj+1,j+1 is extended and set some variables describing 
C     Ci,j+1 and Ei,j+1
C     
                     IF( IWORK(EXAINF+J).EQ.0 ) THEN
                        CCOLS = MIN(NB, (N - NB + IROFFA) - (J-1) * NB)
                        GSIND = J * NB + 1
                        CEXCOL = .FALSE.
                     ELSEIF( IWORK(EXAINF+J).EQ.1 ) THEN
                        CCOLS = NB + 1
                        GSIND = J * NB + 1
                        CEXCOL = .TRUE.
                     ELSEIF( IWORK(EXAINF+J).EQ.2 ) THEN
                        CCOLS = MIN(NB, (N - NB + IROFFA) - 
     $                       (J-1) * NB) - 1
                        GSIND = J * NB + 2
                        CEXCOL = .FALSE.
                     ELSEIF( IWORK(EXAINF+J).EQ.3 ) THEN
                        CCOLS = NB 
                        GSIND = J * NB + 2
                        CEXCOL = .TRUE.
                     END IF
C     
                     CALL INFOG2L( GSI, GSI, DESCSA, NPROW, NPCOL, 
     $                             MYROW, MYCOL, LIA, LJA, ARSRC, 
     $                             ACSRC )
                     CALL INFOG2L( GSI, GSIND, DESCSC, NPROW, NPCOL, 
     $                             MYROW, MYCOL, LIC, LJC, CRSRC, 
     $                             CCSRC )
                     CALL INFOG2L( GSI, GSIND, DESCE, NPROW, NPCOL, 
     $                             MYROW, MYCOL, LEI, LEJ, CRSRC, 
     $                             CCSRC )
C     
                     IF( MYROW.EQ.ARSRC .AND. MYCOL.EQ.ACSRC ) THEN
                        CALL DBEXMAT( EXROW, EXROW, ROWS, ROWS, LIA, 
     $                       LJA, A((LJAS-1)*LLDA+LIAS), LLDA, NBCA, 
     $                       NB, NB, DWORK( EXA ), DWORK( AII ), 
     $                       NB + 1 )
                        IF( CRSRC.NE.ARSRC .OR. CCSRC.NE.ACSRC ) THEN
                           CALL DGESD2D( ICTXT, ROWS, ROWS, 
     $                          DWORK( AII ), NB + 1, CRSRC, CCSRC )
                        END IF
                     END IF
C     
                     IF( MYROW.EQ.CRSRC .AND. MYCOL.EQ.CCSRC ) THEN
                        CALL DBEXMAT( EXROW, CEXCOL, ROWS, CCOLS, LIC, 
     $                       LJC, C((LJCS-1)*LLDC+LICS), LLDC, NBCA, 
     $                       NB, NB, DWORK( EXC ), DWORK( CIJ ), 
     $                       NB + 1 )
                        CALL DBEXMAT( EXROW, CEXCOL, ROWS, CCOLS, LEI, 
     $                       LEJ, DWORK( E ), LLDE, NBCA, NB, NB, 
     $                       DWORK( EXE ), DWORK( EIJ ), NB + 1 )
                        IF( CRSRC.NE.ARSRC .OR. CCSRC.NE.ACSRC ) THEN
                           CALL DGERV2D( ICTXT, ROWS, ROWS, 
     $                          DWORK( AII ), NB + 1, ARSRC, ACSRC )
                        END IF     
C     
C     Do update
C     
                        CALL DGEMM( 'T','N', ROWS, CCOLS, ROWS, -ONE,
     $                       DWORK( AII ), NB + 1, DWORK( EIJ ),
     $                       NB + 1, ONE, DWORK(CIJ), NB + 1 )
C     
C     Save the result in the matrix C
C     
                        CALL DUBEXMA( EXROW, CEXCOL, ROWS, CCOLS, LIC, 
     $                       LJC, C((LJCS-1)*LLDC+LICS), LLDC, NBCA, 
     $                       NB, NB, DWORK( EXC ), DWORK(CIJ), 
     $                       NB + 1 ) 
                     END IF
                  END IF             
               END IF
C     
 145           CONTINUE
C     
C     Update inner loop variable j (column index)
C     
               IF( .NOT. RSIDE ) THEN
                  J = J - 1
               ELSEIF( RSIDE ) THEN
                  J = J + 1
               END IF
C    
 140        CONTINUE
 20      CONTINUE
C     
C     Update outer loop variables
C     
         IF( .NOT. RSIDE ) THEN
            JS = MIN( JS + 1, DBA )
            IE = MIN( IE + 1, DBA )
         ELSEIF( RSIDE ) THEN
            JS = MAX( JS - 1, 1 )
            IE = MAX( IE - 1, 1 )
         END IF
C     
 10   CONTINUE
C     
C     Before we go on we must back redistributed the elements in C
C     
      IF( AEXT )
     $     CALL PDBCKRD( N, N, C, IC, JC, DESCC, IWORK( EXAINF ), 
     $                   IWORK( EXAINF ), DWORK( EXC ), EXMEMC,
     $                   DWORK( IPW ), NB, INFO )
C     
C     Global max on INFO
C     
      IF( NPROCS.GT.1 ) 
     $     CALL IGAMX2D( ICTXT, 'All', ' ', 1, 1, INFO, 1, -1, -1, -1, 
     $                   -1, -1 )
C     
      END
C     
C     End of PTRLYDTD
C     
C *** Last line of PTRLYDTD ***
