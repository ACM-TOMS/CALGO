CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C         This file is a part of the SCASY software library.           C
C         See http://www.cs.umu.se/parallel/research/scasy             C
C         Contributors: Robert Granat and Bo Kågström.                 C
C         SCASY Version 0.10,    March 31, 2006.                       C
C         Copyright 2006, Umeå University, Sweden.                     C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE PTRLYCTD( CSYMM, OP, N, A, IA, JA, DESCA, C, IC, JC, 
     $                     DESCC, NB2, DWORK, LDWORK, IWORK, LIWORK, 
     $                     NOEXSY, SCALE, INFO)
C
C  -- ScaLAPACK-style routine (preliminary version ) --
C     HPC2N and Department of Computing Science,
C     University of Umeå, Sweden.
C     Written by Robert Granat, (granat@cs.umu.se)
C     March 30, 2006.
C
      IMPLICIT NONE
C
C     .. Scalar Arguments ..
      LOGICAL            CSYMM
      CHARACTER*1        OP
      INTEGER            N, IA, JA, IC, JC, LDWORK, LIWORK, NOEXSY, 
     $                   NB2, INFO
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
C  This subroutine solves the real (quasi-)triangular Lyapunov equation
C
C     op( sub( A ) ) * X + X * op( sub(A)^T ) = sub( C ),
C
C  where sub(A) = A(IA:IA+N-1,JA:JA+N-1), sub( C ) = 
C  C(IC:IC+N-1,JC:JC+N-1) and X are N-by-N distributed matrices.  
C
C  The matrix C (and X) may be symmetric (see CSYMM below).
C  The notation op(_) means the transpose or non-transpose of a matrix.
C
C  This routine should *not* be called directly, but through PGELYCTD.
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
C            If OP = 'N', then we solve A * X + X * A^T = C
C            If OP = 'T', then we solve A^T * X + X * A = C.
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
C  NB2       (global input) INTEGER 
C            Internal blocking factor for pipelining of subsolutions
C            for updates of the matrix C (see the references for 
C            details).
C            1 < = NB2 <= DESCC( NB_ ) must hold.
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
C            LIWORK >= DBA + 10 * MIN( P_r, P_c ), where 
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
C  solving the triangular continuous-time Lyapunov equation. See
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
C  In contrary to LAPACK DTRSYL this routine do not scale against 
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
      PARAMETER        ( MONE = -1.0D+0, ONE = 1.0D+0, ZERO = 0.0D+0 )
      INTEGER          BLOCK_CYCLIC_2D, CSRC_, CTXT_, DLEN_, DT_,
     $                 LLD_, MB_, M_, NB_, N_, RSRC_
      PARAMETER        ( BLOCK_CYCLIC_2D = 1, DLEN_ = 9, DT_ = 1,
     $                   CTXT_ = 2, M_ = 3, N_ = 4, MB_ = 5, NB_ = 6,
     $                   RSRC_ = 7, CSRC_ = 8, LLD_ = 9 )
C     ..
C     .. Local Scalars ..
      INTEGER  NB, DBA, NROLL, IS, JS, IE, K, PHASE, PHASES, MNPDIM,
     $         MYCOL, MYROW, NPCOL, NPROW, J, NPROCS, I, IDUM,
     $         AROWS, ACOLS, ROWS, COLS, LINFO, IX, JX, RSRC, CSRC,
     $         LIA, LJA, INDX, ICTXT, LLDA, LLDC, SRSRC, SCSRC,
     $         WRK, MWORKneeded, LMATR, XIJ1, XIJ2, AII, AJJ, CIJ,XIJT, 
     $         RRSRC, RCSRC, ARSRC, ACSRC, ATRSRC, ATCSRC, CRSRC, CCSRC, 
     $         EXA, EXC, EXMEMA, EXMEMC, IWRK, EXAINF, SND, RCDIR,
     $         GI, GJ, NBCA, POS, AROWS2, ATCOLS2, GINDX, ACOLS2, INDXS, 
     $         ATROWS2, LIAT, LJAT, SYMMRSRC, SYMMCSRC, SYMMIX, SYMMJX,
     $         JEND, INDXEND, IPW, IDUMMY, EYE, IROFFA, LIC, LJC,
     $         CSI, CSJ, LIAS, LJAS, LICS, LJCS, GSI, INDXE, INDXU, 
     $         GSIND, ASI, ASJ, GSJ, IROWS, ICOLS, IGSI, IGSJ, IEXRW, 
     $         IEXCL, IIS, IIE, JJS, NIDEEP, NJDEEP, IDEEPS, JDEEPS, 
     $         IDEEPE, JDEEPE, IDEEPU, JDEEPU, IDEEP, JDEEP, XRWS, 
     $         XCLS, XRIND, XCIND, IXRWS, IXCLS, IXRIND, IXCIND, AUPBL, 
     $         CKJ, AUP1, AUP2, AUP3, KAUP1, KAUP2, KAUP3, KKK
      DOUBLE PRECISION SCALOC, DUMMY, LSCALC
      LOGICAL  LQUERY, SNODE, EXROW, EXCOL, CEXROW, CEXCOL, AEXT, RSIDE
C     ..
C     .. Local Arrays ..
      INTEGER IDUM1(1), IDUM2(1), DESCSA( DLEN_ ), DESCSC( DLEN_ ),
     $        IBUFF(4)
      DOUBLE PRECISION DPDUM(1)
C     ..
C     .. External Subroutines ..
      EXTERNAL           BLACS_GRIDINFO, CHK1MAT, DTRSYL, 
     $                   PXERBLA, INFOG2L, DGEMM, DLACPY,
     $                   PCHK2MAT, DSCAL, DGAMN2D, DGESD2D,
     $                   DGERV2D, IGAMX2D, DLATCPY, DSYR2K,
     $                   SB03MD
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
C     Compute the number of blocks to keep from A during the 
C     updates when deep pipelining is turned on
C
      IF( NB2.NE.NB ) THEN
         AUPBL = MAX( 1, ICEIL( DBA-1, MNPDIM ) ) + 1
      ELSE
         AUPBL = 1
      END IF     
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
C     Compute the needed memory for holding data for the extended
C     subsystems. Even if all blocks cannot be extended, we add memory
C     for all blocks to maintain indicies simple.
C
      EXMEMA = ICEIL(AROWS,NB)*ICEIL(ACOLS,NB)*(2*NB+1)
      EXMEMC = EXMEMA 
C
C     Adjust AROWS and ACOLS
C     
      IF( MYROW.EQ.ARSRC ) AROWS = AROWS - IROFFA
      IF( MYCOL.EQ.ACSRC ) ACOLS = ACOLS - IROFFA
C
C     We need two extra matrix descriptors for computing indicies
C     in sub(A) and sub(C). These descriptors are related to the 
C     smallest submatrices of A and C including sub(A) and sub(C)
C     that conform with ScaLAPACK conventions
C
      CALL DESCINIT ( DESCSA, N+IROFFA, N+IROFFA, DESCA(MB_), 
     $                DESCA(NB_), ARSRC, ACSRC, ICTXT, LLDA, INFO )
      CALL DESCINIT ( DESCSC, N+IROFFA, N+IROFFA, DESCC(MB_), 
     $                DESCC(NB_), CRSRC, CCSRC, ICTXT, LLDC, INFO )
C
C     Test work space
C     
      IF( INFO.EQ.0 ) THEN
         WRK = 10*(NB+1)**2 + EXMEMA + EXMEMC + NB + 
     $        3 * AUPBL * (NB+1)**2 + (NB2+1)**2
      END IF
C     
C     Compute needed integer workspace
C     
      IF( INFO.EQ.0 ) THEN
         IWRK = DBA + 10 * MIN( NPROW, NPCOL )
      END IF
C     
C     Check if the call to PTRLYCTD was a workspace query, if so
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
C     Check if we shall continue or interrupt
C     
      IF( INFO.NE.0 ) THEN
         CALL PXERBLA( ICTXT, 'PTRLYCTD', -INFO )
         RETURN
      END IF
C     
C     Init some local pointers into the DWORK-array
C     
      XIJ1 = 1
      XIJ2 = XIJ1  + (NB+1) ** 2
      AII  = XIJ2  + (NB+1) ** 2
      AJJ  = AII  + (NB+1) ** 2
      CIJ  = AJJ  + (NB+1) ** 2
      XIJT = CIJ  + (NB+1) ** 2
      EYE  = XIJT + (NB+1) ** 2
      SND  = EYE  + (NB+1) ** 2
      EXA  = SND  +  NB
      EXC  = EXA  + EXMEMA
      CKJ  = EXC  + EXMEMC
      AUP1 = CKJ + (NB+1)**2
      AUP2 = AUP1 + AUPBL * (NB+1)**2
      IF( CSYMM ) THEN
         AUP3 = AUP2 + AUPBL * (NB+1)**2
         IPW  = AUP3 + AUPBL * (NB+1)**2
      ELSE
         IPW  = AUP2 + AUPBL * (NB+1)**2
      END IF
C
C     Init a local pointer into the IDWORK-array. 
C
      EXAINF = 1
      IROWS = EXAINF + DBA
      ICOLS = IROWS + MNPDIM
      IGSI  = ICOLS + MNPDIM
      IGSJ  = IGSI + MNPDIM
      IEXRW = IGSJ + MNPDIM
      IEXCL = IEXRW + MNPDIM
      IXRWS = IEXCL + MNPDIM
      IXCLS = IXRWS + MNPDIM
      IXRIND = IXCLS + MNPDIM
      IXCIND = IXRIND + MNPDIM
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
     $                    IWORK( EXAINF ), DWORK( EXA ), EXMEMA, DPDUM, 
     $                    IDUM, DWORK( EXC ), EXMEMC, DWORK( SND ), NB, 
     $                    INFO )
C
C     Init a small identity matrix
C
      CALL DLASET( 'All', NB + 1, NB + 1, ZERO, ONE, DWORK( EYE ), 
     $             NB + 1 )
C     
C     Compute the number of block diagonals of C
C     
      NROLL = 2 * DBA - 1
C 
C     Set some loop variables which tells us where to start solving
C     the problem (in matrix C/X)
C
      IF( RSIDE ) THEN
         JS = DBA
         IS = DBA
         IE = DBA
      ELSEIF( .NOT. RSIDE ) THEN
         JS = 1
         IS = 1
         IE = 1
      END IF
C     
      SNODE = .TRUE.
C
C     Compute the number of rounds to do in deep pipelining and
C     set looplimits depending on the transposes
C
      NIDEEP = ICEIL( NB, NB2 )
      NJDEEP = NIDEEP 
      IF( RSIDE  ) THEN
         IDEEPS = NIDEEP
         IDEEPE = 1
         IDEEPU = -1
      ELSE
         IDEEPS = 1
         IDEEPE = NIDEEP
         IDEEPU = 1
      END IF
      IF( .NOT.RSIDE  ) THEN
         JDEEPS = 1
         JDEEPE = NJDEEP
         JDEEPU = 1
      ELSE
         JDEEPS = NJDEEP
         JDEEPE = 1
         JDEEPU = -1
      END IF 
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
C     Main loop over number of diagonals in C
C     
      DO 10 K = 1, NROLL
C     
         IF( RSIDE ) THEN
            IF ( K.GT.DBA ) IS = IS - 1
         ELSEIF( .NOT. RSIDE ) THEN
            IF ( CSYMM .AND. MOD( K, 2 ).EQ.0 ) THEN
               IE = IE + 1
            ELSEIF ( .NOT. CSYMM .AND. K.GT.DBA ) THEN 
               IE = IE + 1
            END IF
         END IF   
C     
         JJS = JS
C     
C     Solve subsystems on the current block diagonal in parallel
C     
         PHASES = ICEIL( IS-IE+1, MNPDIM )
         DO 20 PHASE = 1, PHASES, 1
            IIS = IS-(PHASE-1)*MNPDIM
            IIE = MAX(IS-PHASE*MNPDIM+1,IE)
            JJS = JS + (PHASE-1)*MNPDIM
            DO 22 JDEEP = JDEEPS, JDEEPE, JDEEPU
            DO 24 IDEEP = IDEEPS, IDEEPE, IDEEPU
            SCALOC = ONE
            LINFO = 0
            J = JJS
            DO 30 I = IIS, IIE, -1
C
C     Skip iteration on diagonalblock (I=J) if on upper triangular 
C     part of the block
C
               IF( CSYMM .AND. I.EQ.J .AND. IDEEP.LT.JDEEP ) GOTO 35
C     
C     Here we check if the systems to solve are extended, and extract
C     the necessary data before communicating and calling the solving 
C     routine.
C     
C     Check if Aii is extended and set some variables describing Cij
C     for this particular solve. To handle submatrices, we dstinguish
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
                     ROWS = MIN(NB, (N - NB + IROFFA) - (I - 2) * NB) -1
                     GSI = (I - 1) * NB + 2 
                     EXROW = .FALSE.
                  ELSEIF( IWORK(EXAINF+(I-1)).EQ.3 ) THEN
                     ROWS = NB
                     GSI = (I - 1) * NB + 2 
                     EXROW = .TRUE.
                  END IF
               END IF
               IWORK( IROWS + IIS - I ) = ROWS
               IWORK( IGSI + IIS - I ) = GSI
               IWORK( IEXRW + IIS - I ) = ILG2NT(EXROW)
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
                     COLS = MIN(NB, (N - NB + IROFFA) - (J - 2) * NB) -1
                     GSJ = (J - 1) * NB + 2 
                     EXCOL = .FALSE.
                  ELSEIF( IWORK(EXAINF+(J-1)).EQ.3 ) THEN
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
C     NOTICE: here we use the descriptors for the submatrices instead,
C     which means that all indicies below are calculated for the
C     smallest submatrices sub(A) and sub(C) that conform with ScaLAPACK
C     conventions
C     
               CALL INFOG2L( GSI, GSI, DESCSA, NPROW, NPCOL, MYROW, 
     $              MYCOL, LIA, LJA, ARSRC, ACSRC )
               CALL INFOG2L( GSJ, GSJ, DESCSA, NPROW, NPCOL, MYROW, 
     $              MYCOL, LIAT, LJAT, ATRSRC, ATCSRC )
               CALL INFOG2L( GSI, GSJ, DESCSC, NPROW, NPCOL, MYROW, 
     $              MYCOL, IX, JX, CRSRC, CCSRC ) 
C     
C     Submatrix construction is omitted most of the times for
C     deep pipelining
C     
               IF( IDEEP.NE.IDEEPS .OR. JDEEP.NE.JDEEPS ) GO TO 32
C     
C     Build the extended matrix Aii and send it (if needed) to the 
C     process holding Cij 
C     
               IF( MYROW.EQ.ARSRC .AND. MYCOL.EQ.ACSRC ) THEN
                  CALL DBEXMAT( EXROW, EXROW, ROWS, ROWS, LIA, 
     $                 LJA, A((LJAS-1)*LLDA+LIAS), LLDA, NBCA,
     $                 NB, NB, DWORK( EXA ), DWORK( AII ), 
     $                 NB + 1 )
                  IF( (ARSRC.NE.CRSRC).OR.(ACSRC.NE.CCSRC) ) THEN
                     CALL DGESD2D( ICTXT, ROWS, ROWS, DWORK( AII ), 
     $                    NB + 1, CRSRC, CCSRC )
                  END IF
               END IF
C     
C     Build the extended matrix Ajj and send it to the process holding  
C     Cij
C     
               IF( I.NE.J .AND. (MYROW.EQ.ATRSRC .AND. 
     $              MYCOL.EQ.ATCSRC) ) THEN
                  CALL DBEXMAT( EXCOL, EXCOL, COLS, COLS, LIAT,
     $                 LJAT, A((LJAS-1)*LLDA+LIAS), LLDA, NBCA, 
     $                 NB, NB, DWORK( EXA ), DWORK( AJJ ), 
     $                 NB + 1 )
                  IF( (ATRSRC.NE.CRSRC).OR.(ATCSRC.NE.CCSRC) ) THEN
                     CALL DGESD2D( ICTXT, COLS, COLS, DWORK( AJJ ), 
     $                    NB+1, CRSRC, CCSRC )
                  END IF
               END IF
C     
C     Receive the submatrices Aii and Ajj
C     
               IF( MYROW.EQ.CRSRC .AND. MYCOL.EQ.CCSRC ) THEN
                  IF( ARSRC.NE.CRSRC .OR. ACSRC.NE.CCSRC ) THEN
                     CALL DGERV2D( ICTXT, ROWS, ROWS, DWORK( AII ), 
     $                    NB+1, ARSRC, ACSRC )           
                  END IF
                  IF( I.NE.J .AND. (ATRSRC.NE.CRSRC .OR. 
     $                 ATCSRC.NE.CCSRC) ) THEN
                     CALL DGERV2D( ICTXT, COLS, COLS, DWORK( AJJ ), 
     $                    NB+1, ATRSRC, ATCSRC )
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
                  IF( IDEEP.EQ.IDEEPS .AND. JDEEP.EQ.JDEEPS ) THEN
                     CALL DBEXMAT( EXROW, EXCOL, ROWS, COLS, IX, JX, 
     $                    C((LJCS-1)*LLDC+LICS), LLDC, NBCA, NB, NB,  
     $                    DWORK( EXC ), DWORK( CIJ ), NB + 1 ) 
                     IF( CSYMM .AND. I.EQ.J )
     $                    CALL DLATCPY( 'Lower', ROWS-1, COLS-1, 
     $                    DWORK(CIJ+1), NB+1, DWORK(CIJ+NB+1), NB+1 )
                  END IF
C     
C     Solve sub-system, choose branch depending on RSIDE 
C     
                  IF( I.EQ.J ) THEN
                     IF( CSYMM ) THEN
                        CALL DTRLYCT( 'S', OP, IDEEP, JDEEP, ROWS, NB2, 
     $                       DWORK(AII), NB+1, DWORK(CIJ), NB+1, 
     $                       DWORK(IPW), (NB2+1)**2, XRWS, XCLS, XRIND, 
     $                       XCIND, SCALOC, LINFO )
                     ELSE
                        IF( RSIDE ) THEN
                           CALL DTRSYCT( 'S', 'N', 'T', +1, IDEEP, 
     $                          JDEEP, ROWS, NB2, COLS, NB2, 
     $                          DWORK(AII), NB+1, DWORK(AII), NB+1, 
     $                          DWORK(CIJ), NB+1, XRWS, XCLS, XRIND, 
     $                          XCIND, SCALOC, LINFO ) 
                        ELSEIF( .NOT. RSIDE ) THEN
                           CALL DTRSYCT( 'S', 'T', 'N', +1, IDEEP, 
     $                          JDEEP, ROWS, NB2, COLS, NB2, 
     $                          DWORK(AII), NB+1, DWORK(AII), NB+1, 
     $                          DWORK(CIJ), NB+1, XRWS, XCLS, XRIND, 
     $                          XCIND, SCALOC, LINFO ) 
                        END IF
                     END IF
                  ELSE
                     IF( RSIDE ) THEN
                        CALL DTRSYCT( 'S', 'N', 'T', +1, IDEEP, 
     $                       JDEEP, ROWS, NB2, COLS, NB2, 
     $                       DWORK(AII), NB+1, DWORK(AJJ), NB+1, 
     $                       DWORK(CIJ), NB+1, XRWS, XCLS, XRIND, 
     $                       XCIND, SCALOC, LINFO ) 
                     ELSEIF( .NOT. RSIDE ) THEN
                        CALL DTRSYCT( 'S', 'T', 'N', +1, IDEEP, 
     $                       JDEEP, ROWS, NB2, COLS, NB2, 
     $                       DWORK(AII), NB+1, DWORK(AJJ), NB+1, 
     $                       DWORK(CIJ), NB+1, XRWS, XCLS, XRIND, 
     $                       XCIND, SCALOC, LINFO ) 
                     END IF
                  END IF
                  IWORK( IXRWS + IIS - I ) = XRWS
                  IWORK( IXCLS + IIS - I ) = XCLS
                  IWORK( IXRIND + IIS - I ) = XRIND
                  IWORK( IXCIND + IIS - I ) = XCIND
C     
                  IF (LINFO.NE.0) THEN
                     INFO = LINFO
                  END IF
C     
C     Copy the given solution back to the global matrix and update the
C     extension elements
C     
                  IF( IDEEP.EQ.IDEEPE .AND. JDEEP.EQ.JDEEPE .AND. 
     $                SCALOC.EQ.ONE ) THEN
                     IF( CSYMM .AND. I.EQ.J )
     $                    CALL DLATCPY( 'Upper', ROWS-1, COLS-1, 
     $                    DWORK(CIJ+NB+1), NB+1, DWORK(CIJ+1), NB+1 )
                     CALL DUBEXMA( EXROW, EXCOL, ROWS, COLS, IX, JX, 
     $                    C((LJCS-1)*LLDC+LICS), LLDC, NBCA, NB, 
     $                    NB, DWORK( EXC ), DWORK( CIJ ), NB + 1 )
                  END IF
C     
C     Warn for owerflow 
C     
                  IF (SCALOC.NE.ONE) THEN
                     INFO = 3
                  END IF
C     
C     Copy the solution to broadcast into two buffers, one for each 
C     direction
C     
                  CALL DLACPY( 'All', XRWS, XCLS, 
     $                 DWORK( CIJ + (XCIND-1)*(NB+1) + XRIND - 1), NB+1,  
     $                 DWORK( XIJ1 + (XCIND-1)*(NB+1) + XRIND - 1), 
     $                 NB+1 )
                  CALL DLACPY( 'All', XRWS, XCLS, 
     $                 DWORK( CIJ + (XCIND-1)*(NB+1) + XRIND - 1), NB+1,  
     $                 DWORK( XIJ2 + (XCIND-1)*(NB+1) + XRIND - 1), 
     $                 NB+1 )
               END IF
 35            CONTINUE
C     
C     Update inner loop variable j (column index)
C     
               J = J + 1
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
              J = J + 1
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
 123           CONTINUE
               IF( AEXT )
     $              CALL DSCAL( EXMEMC, SCALOC, DWORK( EXC ), 1 )
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
                  IF( ROWS.EQ.0 .OR. COLS.EQ.0 ) GO TO 223
C     
                  CALL INFOG2L( GSI, GSJ, DESCSC, NPROW, NPCOL, 
     $                 MYROW, MYCOL, IX, JX, SRSRC, SCSRC )     
                  SNODE = MYROW.EQ.SRSRC.AND.MYCOL.EQ.SCSRC
C
                  IF( SNODE ) THEN 
                     XRWS = IWORK( IXRWS + IIS - I ) 
                     XCLS = IWORK( IXCLS + IIS - I ) 
                     XRIND = IWORK( IXRIND + IIS - I )
                     XCIND = IWORK( IXCIND + IIS - I )
                     IF( SCALOC.NE.LSCALC ) THEN
                        DO 333 KKK = 1, XCLS
                           CALL DSCAL( XRWS, SCALOC / LSCALC, 
     $                          DWORK( XIJ1 + (XCIND-1)*(NB+1) + 
     $                          XRIND - 1 + (KKK-1)*(NB+1) ), 1 )
 333                    CONTINUE
                        DO 444 KKK = 1, XCLS
                           CALL DSCAL( XRWS, SCALOC / LSCALC, 
     $                          DWORK( XIJ2 + (XCIND-1)*(NB+1) + 
     $                          XRIND - 1 + (KKK-1)*(NB+1) ), 1 )
 444                    CONTINUE
                     END IF
                     DO 555 KKK = 1, COLS
                        CALL DSCAL( ROWS, SCALOC, 
     $                       DWORK( CIJ + (KKK-1)*(NB+1) ), 1 )
 555                 CONTINUE
                     CALL DLACPY( 'All', XRWS, XCLS, 
     $                    DWORK( XIJ1 + (XCIND-1)*(NB+1) +
     $                    XRIND - 1), NB+1,  DWORK( CIJ + 
     $                    (XCIND-1)*(NB+1) + XRIND - 1), NB+1 )
                     IF( IDEEP.EQ.IDEEPE .AND. JDEEP.EQ.JDEEPE ) THEN
                        IF( CSYMM )
     $                       CALL DLATCPY( 'Upper', ROWS-1, COLS-1, 
     $                       DWORK(CIJ+NB+1), NB+1, DWORK(CIJ+1), NB+1 )
                        EXROW = INT2LG(IWORK( IEXRW + IIS - I )) 
                        EXCOL = INT2LG(IWORK( IEXCL + IIS - I ))
                        CALL DUBEXMA( EXROW, EXCOL, ROWS, COLS, IX, 
     $                       JX, C((LJCS-1)*LLDC+LICS), LLDC, NBCA, NB, 
     $                       NB, DWORK( EXC ), DWORK( CIJ ), NB + 1 )
                     END IF
                  END IF
C
C     In case of symmetry, scale the local block at the transposed
C     position as well
C
                  IF( CSYMM ) THEN
                     CALL INFOG2L( GSJ, GSI, DESCSC, NPROW, NPCOL, 
     $                    MYROW, MYCOL, SYMMIX, SYMMJX, SYMMRSRC, 
     $                    SYMMCSRC )
                     IF( MYROW.EQ.SYMMRSRC.AND.MYCOL.EQ.SYMMCSRC ) THEN
                        DO 666 KKK = 1, ROWS
                           CALL DSCAL( COLS, SCALOC, 
     $                          DWORK( XIJT + (KKK-1)*(NB+1) ), 1 )
 666                    CONTINUE
                     END IF
                  END IF
C
 223              CONTINUE
C
                J = J + 1
 222           CONTINUE
C     
C     Update value of SCALE according to SCALOC
C     
               SCALE = SCALE * SCALOC
            END IF
C     
C     In the case of symmetry (and we are not located on the main block
C     diagonal), we possibly send the current subsolution
C     to the process holding Cji, if that node hasn't received it yet
C     
            IF( CSYMM ) THEN
               J = JJS
               DO 36 I = IIS, IIE, -1
                  IF( I.GT.J ) THEN
C     
                     ROWS = IWORK( IROWS + IIS - I )
                     GSI = IWORK( IGSI + IIS - I ) 
                     EXROW = INT2LG(IWORK( IEXRW + IIS - I )) 
                     COLS = IWORK( ICOLS + IIS - I ) 
                     GSJ = IWORK( IGSJ + IIS - I )
                     EXCOL = INT2LG(IWORK( IEXCL + IIS - I ))
C     
                     IF( ROWS.EQ.0 .OR. COLS.EQ.0 ) GO TO 37
C
                     CALL INFOG2L( GSI, GSJ, DESCSC, NPROW, NPCOL, 
     $                    MYROW, MYCOL, IX, JX, SRSRC, SCSRC )
                     CALL INFOG2L( GSJ, GSI, DESCSC, NPROW, NPCOL, 
     $                    MYROW, MYCOL, SYMMIX, SYMMJX, SYMMRSRC, 
     $                    SYMMCSRC )
C
                     IF( SYMMRSRC.NE.SRSRC .OR. SYMMCSRC.NE.SCSRC ) THEN
                        IF( MYROW.EQ.SRSRC .AND. MYCOL.EQ.SCSRC ) THEN
                           IBUFF(1) = IWORK( IXRWS + IIS - I )
                           IBUFF(2) = IWORK( IXCLS + IIS - I )
                           IBUFF(3) = IWORK( IXRIND + IIS - I )
                           IBUFF(4) = IWORK( IXCIND + IIS - I )
                           CALL IGESD2D( ICTXT, 4, 1, IBUFF, 4, 
     $                          SYMMRSRC, SYMMCSRC )
                           XRWS = IWORK( IXRWS + IIS - I ) 
                           XCLS = IWORK( IXCLS + IIS - I ) 
                           XRIND = IWORK( IXRIND + IIS - I )
                           XCIND = IWORK( IXCIND + IIS - I )
                           CALL DGESD2D( ICTXT, XRWS, XCLS, 
     $                          DWORK( XIJ1+(XCIND-1)*(NB+1) + XRIND-1),
     $                          NB + 1, SYMMRSRC, SYMMCSRC )
                        END IF
                        IF(MYROW.EQ.SYMMRSRC.AND.MYCOL.EQ.SYMMCSRC) THEN
                           CALL IGERV2D( ICTXT, 4, 1, IBUFF, 4, SRSRC, 
     $                          SCSRC )
                           IWORK( IXRWS + IIS - I ) = IBUFF(1)
                           IWORK( IXCLS + IIS - I ) = IBUFF(2)
                           IWORK( IXRIND + IIS - I ) = IBUFF(3)
                           IWORK( IXCIND + IIS - I ) = IBUFF(4)
                           XRWS = IWORK( IXRWS + IIS - I ) 
                           XCLS = IWORK( IXCLS + IIS - I ) 
                           XRIND = IWORK( IXRIND + IIS - I )
                           XCIND = IWORK( IXCIND + IIS - I )
                           CALL DGERV2D( ICTXT, XRWS, XCLS, 
     $                          DWORK( IPW ), XRWS, SRSRC, SCSRC )
                        END IF   
                     END IF
C     
                     IF( MYROW.EQ.SYMMRSRC.AND.MYCOL.EQ.SYMMCSRC ) THEN
                        IF(SYMMRSRC.NE.SRSRC.OR.SYMMCSRC.NE.SCSRC) THEN
                           CALL DLATCPY( 'All', XRWS, XCLS, 
     $                          DWORK( IPW ), XRWS, 
     $                          DWORK( XIJT+(XRIND-1)*(NB+1) + XCIND-1), 
     $                          NB + 1 )
                        ELSE
                           CALL DLATCPY( 'All', XRWS, XCLS, 
     $                          DWORK( XIJ1+(XCIND-1)*(NB+1) + XRIND-1), 
     $                          NB + 1, 
     $                          DWORK( XIJT+(XRIND-1)*(NB+1) + XCIND-1), 
     $                          NB + 1 )
                        END IF
C     
C     Put Xij^T back into the local part of the global matrix
C     
                        IF( IDEEP.EQ.IDEEPE .AND. JDEEP.EQ.JDEEPE )
     $                       CALL DUBEXMA( EXCOL, EXROW, COLS, ROWS, 
     $                       SYMMIX, SYMMJX, C((LJCS-1)*LLDC+LICS), 
     $                       LLDC, NBCA, NB, NB, DWORK( EXC ), 
     $                       DWORK( XIJT ), NB + 1 )
                     END IF
                  END IF
 37               CONTINUE
C     
C     Update inner loop variable j (column index)
C     
                  J = J + 1
 36            CONTINUE
            END IF
C     
C     Broadcasting of the solutions in corresponding block rows and column
C     
            IF( K.LT.NROLL ) THEN
               DO 39 RCDIR = 1, 2
                  J = JJS
                  DO 40 I = IIS, IIE, -1
C     
C
C     Skip iteration on diagonalblock (I=J) if on upper triangular 
C     part of the block
C
               IF( CSYMM .AND. I.EQ.J .AND. IDEEP.LT.JDEEP ) GOTO 45
C
                     ROWS = IWORK( IROWS + IIS - I )
                     GSI = IWORK( IGSI + IIS - I ) 
                     EXROW = INT2LG(IWORK( IEXRW + IIS - I )) 
                     COLS = IWORK( ICOLS + IIS - I ) 
                     GSJ = IWORK( IGSJ + IIS - I )
                     EXCOL = INT2LG(IWORK( IEXCL + IIS - I ))
C     
                     IF( ROWS.EQ.0 .OR. COLS.EQ.0 ) GO TO 45
C     
                     CALL INFOG2L( GSI, GSJ, DESCSC, NPROW, NPCOL, 
     $                    MYROW, MYCOL, IX, JX, CRSRC, CCSRC )
C     
                     SNODE = MYROW.EQ.CRSRC.AND.MYCOL.EQ.CCSRC
                     SRSRC = CRSRC
                     SCSRC = CCSRC
C     
                     IF( SNODE ) THEN
                        XRWS = IWORK( IXRWS + IIS - I ) 
                        XCLS = IWORK( IXCLS + IIS - I ) 
                        XRIND = IWORK( IXRIND + IIS - I )
                        XCIND = IWORK( IXCIND + IIS - I )
                        IBUFF(1) = XRWS
                        IBUFF(2) = XCLS
                        IBUFF(3) = XRIND
                        IBUFF(4) = XCIND
C     
                        IF (NPCOL.GT.1 .AND. RCDIR.EQ.1 ) THEN
                           CALL IGEBS2D( ICTXT, 'Row', ' ', 4, 
     $                          1, IBUFF, 4 )
                           CALL DGEBS2D( ICTXT, 'Row', ' ', XRWS, 
     $                          XCLS,
     $                          DWORK(XIJ1+(XCIND-1)*(NB+1)+XRIND-1), 
     $                          NB+1 )
                        END IF
                        IF (NPROW.GT.1 .AND. RCDIR.EQ.2 ) THEN
                           CALL IGEBS2D( ICTXT, 'Col', ' ', 4, 
     $                          1, IBUFF, 4 )
                           CALL DGEBS2D( ICTXT, 'Col', ' ', XRWS, 
     $                          XCLS,
     $                          DWORK(XIJ2+(XCIND-1)*(NB+1)+XRIND-1), 
     $                          NB+1 )
                        END IF
                     ELSE
C     
                        IF ( NPCOL.GT.1 .AND. MYROW.EQ.SRSRC
     $                       .AND. RCDIR.EQ.1 ) THEN
                           CALL IGEBR2D( ICTXT, 'Row', ' ', 4, 
     $                          1, IBUFF, 4, SRSRC, SCSRC )
                           IWORK( IXRWS + IIS - I ) = IBUFF(1)
                           IWORK( IXCLS + IIS - I ) = IBUFF(2)
                           IWORK( IXRIND + IIS - I ) = IBUFF(3)
                           IWORK( IXCIND + IIS - I ) = IBUFF(4)
                           XRWS = IWORK( IXRWS + IIS - I ) 
                           XCLS = IWORK( IXCLS + IIS - I ) 
                           XRIND = IWORK( IXRIND + IIS - I )
                           XCIND = IWORK( IXCIND + IIS - I )
                           CALL DGEBR2D( ICTXT, 'Row', ' ', XRWS, 
     $                          XCLS,
     $                          DWORK(XIJ1+(XCIND-1)*(NB+1)+XRIND-1), 
     $                          NB+1, SRSRC, SCSRC )
                        END IF
                        IF ( NPROW.GT.1 .AND. MYCOL.EQ.SCSRC
     $                       .AND. RCDIR.EQ.2 ) THEN 
                           CALL IGEBR2D( ICTXT, 'Col', ' ', 4, 
     $                          1, IBUFF, 4, SRSRC, SCSRC )
                           IWORK( IXRWS + IIS - I ) = IBUFF(1)
                           IWORK( IXCLS + IIS - I ) = IBUFF(2)
                           IWORK( IXRIND + IIS - I ) = IBUFF(3)
                           IWORK( IXCIND + IIS - I ) = IBUFF(4)
                           XRWS = IWORK( IXRWS + IIS - I ) 
                           XCLS = IWORK( IXCLS + IIS - I ) 
                           XRIND = IWORK( IXRIND + IIS - I )
                           XCIND = IWORK( IXCIND + IIS - I )
                           CALL DGEBR2D( ICTXT, 'Col', ' ', XRWS, 
     $                          XCLS,
     $                          DWORK(XIJ2+(XCIND-1)*(NB+1)+XRIND-1), 
     $                          NB+1, SRSRC, SCSRC )
                        END IF
                     END IF
C     
C     Now broadcast the subsolution Xji to processors holding blocks
C     in block row j of C or block column i, depending on RSIDE
C     
                     IF( CSYMM .AND. I.GT.J ) THEN
                        CALL INFOG2L( GSJ, GSI, DESCSC, NPROW, NPCOL, 
     $                       MYROW, MYCOL, SYMMIX, SYMMJX, SYMMRSRC, 
     $                       SYMMCSRC )
                        IF( MYROW.EQ.SYMMRSRC.AND.MYCOL.EQ.SYMMCSRC ) 
     $                       THEN
                           IF ( RSIDE ) THEN
                              IF ( NPCOL.GT.1 .AND. RCDIR.EQ.1 ) THEN
                                 IF( SYMMRSRC.NE.SRSRC ) THEN
                                    XRWS = IWORK( IXRWS + IIS - I ) 
                                    XCLS = IWORK( IXCLS + IIS - I ) 
                                    XRIND = IWORK( IXRIND + IIS - I )
                                    XCIND = IWORK( IXCIND + IIS - I )
                                    IBUFF(1) = XRWS
                                    IBUFF(2) = XCLS
                                    IBUFF(3) = XRIND
                                    IBUFF(4) = XCIND
                                    CALL IGEBS2D( ICTXT, 'Row', ' ', 4, 
     $                                   1, IBUFF, 4 )
                                    CALL DGEBS2D( ICTXT, 'Row', ' ', 
     $                                   XCLS, XRWS, 
     $                                   DWORK( XIJT+(XRIND-1)*(NB+1) + 
     $                                   XCIND - 1 ), NB + 1 )
                                 END IF
                              END IF
                           ELSEIF ( .NOT. RSIDE ) THEN
                              IF ( NPROW.GT.1 .AND. RCDIR.EQ.2 ) THEN
                                 IF( SYMMCSRC.NE.SCSRC ) THEN
                                    XRWS = IWORK( IXRWS + IIS - I ) 
                                    XCLS = IWORK( IXCLS + IIS - I ) 
                                    XRIND = IWORK( IXRIND + IIS - I )
                                    XCIND = IWORK( IXCIND + IIS - I )
                                    IBUFF(1) = XRWS
                                    IBUFF(2) = XCLS
                                    IBUFF(3) = XRIND
                                    IBUFF(4) = XCIND
                                    CALL IGEBS2D( ICTXT, 'Col', ' ', 4, 
     $                                   1, IBUFF, 4 )
                                    CALL DGEBS2D( ICTXT, 'Col', ' ', 
     $                                   XCLS, XRWS, 
     $                                   DWORK( XIJT+(XRIND-1)*(NB+1) + 
     $                                   XCIND - 1 ), NB + 1 )
                                 END IF
                              END IF
                           END IF
C     
C     The processors holding blocks in block row j or block column
C     i receives Xji, again depending on RSIDE
C     
                        ELSE
                           IF ( RSIDE ) THEN
                              IF( NPCOL.GT.1 .AND. MYROW.EQ.SYMMRSRC 
     $                             .AND. RCDIR.EQ.1 ) THEN
                                 IF( SYMMRSRC.NE.SRSRC ) THEN
                                    CALL IGEBR2D( ICTXT, 'Row', ' ', 4, 
     $                                   1, IBUFF, 4, SYMMRSRC, 
     $                                   SYMMCSRC )
                                    IWORK( IXRWS + IIS - I ) = IBUFF(1)
                                    IWORK( IXCLS + IIS - I ) = IBUFF(2)
                                    IWORK( IXRIND + IIS - I ) = IBUFF(3)
                                    IWORK( IXCIND + IIS - I ) = IBUFF(4)
                                    XRWS = IWORK( IXRWS + IIS - I ) 
                                    XCLS = IWORK( IXCLS + IIS - I ) 
                                    XRIND = IWORK( IXRIND + IIS - I )
                                    XCIND = IWORK( IXCIND + IIS - I )
                                    CALL DGEBR2D( ICTXT, 'Row', ' ', 
     $                                   XCLS, XRWS, 
     $                                   DWORK( XIJT+(XRIND-1)*(NB+1) + 
     $                                   XCIND - 1 ), 
     $                                   NB + 1, SYMMRSRC, SYMMCSRC )
                                 ELSE
                                    CALL DLATCPY( 'All', XRWS, XCLS, 
     $                                   DWORK( XIJ1+(XCIND-1)*(NB+1) + 
     $                                   XRIND - 1 ), NB + 1, 
     $                                   DWORK( XIJT+(XRIND-1)*(NB+1) + 
     $                                   XCIND - 1 ), NB + 1 )
                                 END IF
                              END IF
                           ELSEIF( .NOT. RSIDE ) THEN
                              IF( NPROW.GT.1 .AND. MYCOL.EQ.SYMMCSRC 
     $                             .AND. RCDIR.EQ.2 ) THEN
                                 IF( SYMMCSRC.NE.SCSRC ) THEN
                                    CALL IGEBR2D( ICTXT, 'Col', ' ', 4, 
     $                                   1, IBUFF, 4, SYMMRSRC, 
     $                                   SYMMCSRC )
                                    IWORK( IXRWS + IIS - I ) = IBUFF(1)
                                    IWORK( IXCLS + IIS - I ) = IBUFF(2)
                                    IWORK( IXRIND + IIS - I ) = IBUFF(3)
                                    IWORK( IXCIND + IIS - I ) = IBUFF(4)
                                    XRWS = IWORK( IXRWS + IIS - I ) 
                                    XCLS = IWORK( IXCLS + IIS - I ) 
                                    XRIND = IWORK( IXRIND + IIS - I )
                                    XCIND = IWORK( IXCIND + IIS - I )
                                    CALL DGEBR2D( ICTXT, 'Col', ' ', 
     $                                   XCLS, XRWS, 
     $                                   DWORK( XIJT+(XRIND-1)*(NB+1) + 
     $                                   XCIND - 1 ), 
     $                                   NB + 1, SYMMRSRC, SYMMCSRC )
                                 ELSE
                                    CALL DLATCPY( 'All', XRWS, XCLS, 
     $                                   DWORK( XIJ2+(XCIND-1)*(NB+1) + 
     $                                   XRIND - 1 ), NB + 1, 
     $                                   DWORK( XIJT+(XRIND-1)*(NB+1) + 
     $                                   XCIND - 1 ), NB + 1 )
                                 END IF
                              END IF
                           END IF
                        END IF
                     END IF
C     
 45                  CONTINUE
C     
C     Update inner loop variable j (column index)
C     
                     J = J + 1
 40               CONTINUE
 39            CONTINUE
            END IF
C     
C     Update subsolution in deep pipelining
C     
            J = JJS
            DO 47 I = IIS, IIE, -1
C     
C
C     Skip iteration on diagonalblock (I=J) if on upper triangular 
C     part of the block
C
               IF( CSYMM .AND. I.EQ.J .AND. IDEEP.LT.JDEEP ) GOTO 48
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
                  XCLS = IWORK( IXCLS + IIS - I ) 
                  XRIND = IWORK( IXRIND + IIS - I )
                  XCIND = IWORK( IXCIND + IIS - I )
                  IF( I.EQ.J ) THEN
                     IF( CSYMM ) THEN
                        CALL DTRLYCT( 'U', OP, IDEEP, JDEEP, ROWS, NB2, 
     $                       DWORK(AII), NB+1, DWORK(CIJ), NB+1, 
     $                       DWORK(IPW), (NB2+1)**2, XRWS, XCLS, XRIND, 
     $                       XCIND, SCALOC, LINFO )
                     ELSE
                        IF( RSIDE ) THEN
                           CALL DTRSYCT( 'U', 'N', 'T', +1, IDEEP, 
     $                          JDEEP, ROWS, NB2, COLS, NB2, 
     $                          DWORK(AII), NB+1, DWORK(AII), NB+1, 
     $                          DWORK(CIJ), NB+1, XRWS, XCLS, XRIND, 
     $                          XCIND, SCALOC, LINFO ) 
                        ELSEIF( .NOT. RSIDE ) THEN
                           CALL DTRSYCT( 'U', 'T', 'N', +1, IDEEP, 
     $                          JDEEP, ROWS, NB2, COLS, NB2, 
     $                          DWORK(AII), NB+1, DWORK(AII), NB+1, 
     $                          DWORK(CIJ), NB+1, XRWS, XCLS, XRIND, 
     $                          XCIND, SCALOC, LINFO ) 
                        END IF
                     END IF
                  ELSE
                     IF( RSIDE ) THEN
                        CALL DTRSYCT( 'U', 'N', 'T', +1, IDEEP, 
     $                       JDEEP, ROWS, NB2, COLS, NB2, 
     $                       DWORK(AII), NB+1, DWORK(AJJ), NB+1, 
     $                       DWORK(CIJ), NB+1, XRWS, XCLS, XRIND, 
     $                       XCIND, SCALOC, LINFO ) 
                     ELSEIF( .NOT. RSIDE ) THEN
                        CALL DTRSYCT( 'U', 'T', 'N', +1, IDEEP, 
     $                       JDEEP, ROWS, NB2, COLS, NB2, 
     $                       DWORK(AII), NB+1, DWORK(AJJ), NB+1, 
     $                       DWORK(CIJ), NB+1, XRWS, XCLS, XRIND, 
     $                       XCIND, SCALOC, LINFO ) 
                     END IF
                  END IF
               END IF
 48            CONTINUE
C     
C     Update inner loop variable j (column index)
C     
               J = J + 1
 47         CONTINUE
C     
C     Update rest of global system wrt to current subsolution. 
C     Choose branch depending on RSIDE.
C     This first branch corresponds to solving A*X + X*A^T = C
C     
            J = JJS
            DO 50 I = IIS, IIE, -1
C
C     Skip iteration on diagonalblock (I=J) if on upper triangular 
C     part of the block
C
               IF( CSYMM .AND. I.EQ.J .AND. IDEEP.LT.JDEEP ) GOTO 55
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
               IF ( RSIDE ) THEN
C     
C     For every update we do we also need to find out the dimensions and
C     possible extensions of the submatrices involved.
C     
C     Start with block column j, that is, do the operation (ISTART is
C     decided by the symmetric properties)
C     
C     Ckj = Ckj - Aki*Xij, for all k that belongs to [ISTART,i-1]
C     
C     If C is symmetric, and we shall update Ckj, where k=j,
C     we do a SYR2K-update, namely the operation
C     
C     Ckj = Ckj - Aki*Xij - Xji*Aki^T = Ckj- Aki*Xji^T - Xji*Aki^T 
C     
C     To do this as a regular SYR2K-operation the matrix Xij must have 
C     been transposed locally to Xji before the call to the BLAS routine
C     DSYR2K.
C     
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
                  KAUP1 = 1
                  DO 60 INDX =  INDXS, INDXE, INDXU
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
                     CALL INFOG2L( GSIND, GSI, DESCSA, NPROW, NPCOL, 
     $                             MYROW, MYCOL, LIA, LJA, RSRC, CSRC )
                     CALL INFOG2L( GSIND, GSJ, DESCSC, NPROW, NPCOL, 
     $                             MYROW, MYCOL, IX, JX, RRSRC, RCSRC )
C
C     Submatrix constuction is omitted most of the times for 
C     deep pipelining
C
                     IF( IDEEP.NE.IDEEPS .OR. JDEEP.NE.JDEEPS )
     $                    GO TO 62
C
                     IF( MYROW.EQ.RSRC.AND.MYCOL.EQ.CSRC ) THEN
                        IF( (RSRC.NE.RRSRC).OR.(CSRC.NE.RCSRC) ) THEN
                           CALL DBEXMAT( CEXROW, EXROW, AROWS2, ROWS, 
     $                          LIA, LJA, A((LJAS-1)*LLDA+LIAS), LLDA, 
     $                          NBCA, NB, NB, DWORK( EXA ), 
     $                          DWORK( AUP1+(AUPBL-1)*(NB+1)**2 ), 
     $                          NB + 1 )
                           CALL DGESD2D( ICTXT, AROWS2, ROWS, 
     $                          DWORK( AUP1+(AUPBL-1)*(NB+1)**2 ), 
     $                          NB + 1, RRSRC, RCSRC )
                        ELSE
                           CALL DBEXMAT( CEXROW, EXROW, AROWS2, ROWS, 
     $                          LIA, LJA, A((LJAS-1)*LLDA+LIAS), LLDA, 
     $                          NBCA, NB, NB, DWORK( EXA ), 
     $                          DWORK( AUP1+(KAUP1-1)*(NB+1)**2 ), 
     $                          NB + 1 )
                        END IF
                     END IF
C     
C    Skipped submatrix construction?
C     
 62                  CONTINUE
C     
                     IF(MYROW.EQ.RRSRC.AND.MYCOL.EQ.RCSRC) THEN
                        CALL DBEXMAT( CEXROW, EXCOL, AROWS2, COLS, IX, 
     $                       JX, C((LJCS-1)*LLDC+LICS), LLDC, 
     $                       NBCA, NB, NB, DWORK( EXC ), 
     $                       DWORK(CKJ), NB + 1 )
                        IF( IDEEP.NE.IDEEPS .OR. JDEEP.NE.JDEEPS )
     $                       GO TO 64
                        IF( (RSRC.NE.RRSRC) .OR. (CSRC.NE.RCSRC) ) THEN
                           CALL DGERV2D( ICTXT, AROWS2, ROWS, 
     $                          DWORK( AUP1+(KAUP1-1)*(NB+1)**2 ), 
     $                          NB + 1, RSRC, CSRC )
                        END IF
 64                     CONTINUE
                        XRWS = IWORK( IXRWS + IIS - I ) 
                        XCLS = IWORK( IXCLS + IIS - I ) 
                        XRIND = IWORK( IXRIND + IIS - I )
                        XCIND = IWORK( IXCIND + IIS - I )
                        IF( NB2.EQ.NB .AND. CSYMM .AND. INDX.EQ.J ) THEN
                           CALL DSYR2K( 'U', 'N', XCLS, XRWS, 
     $                          -1.0D+00, 
     $                          DWORK(AUP1+(KAUP1-1)*(NB+1)**2+
     $                          (XRIND-1)*(NB+1)+XCIND-1 ), NB + 1, 
     $                          DWORK(XIJT+(XRIND-1)*(NB+1)+XCIND-1), 
     $                          NB + 1, ONE, 
     $                          DWORK(CKJ+(XCIND-1)*(NB+1)+XCIND-1), 
     $                          NB + 1)
                           CALL DLATCPY( 'Upper', XCLS-1, XCLS-1, 
     $                          DWORK(CKJ+XCIND*(NB+1)+XCIND-1), NB + 1,
     $                          DWORK(CKJ+(XCIND-1)*(NB+1)+XCIND), 
     $                          NB + 1 )
                        ELSE
                           CALL DGEMM( 'N', 'N', AROWS2, XCLS, XRWS, 
     $                          -1.0D+00, 
     $                          DWORK(AUP1+(KAUP1-1)*(NB+1)**2+
     $                          (XRIND-1)*(NB+1) ), NB + 1, 
     $                          DWORK(XIJ2+(XCIND-1)*(NB+1)+XRIND-1), 
     $                          NB + 1, ONE, 
     $                          DWORK(CKJ+(XCIND-1)*(NB+1)), NB + 1 )
                           IF( CSYMM .AND. I.EQ.J .AND. IDEEP.GT.JDEEP )
     $                          CALL DGEMM( 'N', 'T', AROWS2, XRWS, 
     $                          XCLS, 
     $                          -1.0D+00, 
     $                          DWORK(AUP1+(KAUP1-1)*(NB+1)**2+
     $                          (XCIND-1)*(NB+1) ), NB + 1, 
     $                          DWORK(XIJ2+(XCIND-1)*(NB+1)+XRIND-1), 
     $                          NB + 1, ONE, 
     $                          DWORK(CKJ+(XRIND-1)*(NB+1)), NB + 1 )   
                        END IF
                        CALL DUBEXMA( CEXROW, EXCOL, AROWS2, COLS, IX, 
     $                       JX, C((LJCS-1)*LLDC+LICS), LLDC, 
     $                       NBCA, NB, NB, DWORK( EXC ), 
     $                       DWORK(CKJ), NB + 1 )
                        IF( AUPBL.NE.1 ) KAUP1 = KAUP1 + 1
                     END IF
 60               CONTINUE   
C     
C     Now update block row i of C, that is do the operation
C     
C     Cik = Cik - Xij * Akj^T, for all k that belongs to [1,j-1]  
C     
                  IF( MOD( I, 2 ).EQ.0 ) THEN
                     INDXS = 1
                     INDXE = J - 1
                     INDXU = 1
                  ELSE
                     INDXS = J - 1
                     INDXE = 1
                     INDXU = -1
                  END IF
                  KAUP2 = 1
                  DO 70 INDX =  INDXS, INDXE, INDXU
                     IF( INDX.EQ.1 ) THEN
                        IF( IWORK(EXAINF).EQ.0 ) THEN
                           ATROWS2 = NB - IROFFA
                           GSIND = 1 + IROFFA
                           CEXCOL = .FALSE.
                        ELSEIF( IWORK(EXAINF).EQ.1 ) THEN
                           ATROWS2 = NB - IROFFA + 1
                           GSIND = 1 + IROFFA
                           CEXCOL = .TRUE.
                        END IF
                     ELSE
                        IF( IWORK(EXAINF+(INDX-1)).EQ.0 ) THEN
                           ATROWS2 = MIN(NB, (N - NB + IROFFA) - 
     $                          (INDX - 2) * NB)
                           GSIND = (INDX - 1) * NB + 1 
                           CEXCOL = .FALSE.
                        ELSEIF( IWORK(EXAINF+(INDX-1)).EQ.1 ) THEN
                           ATROWS2 = NB + 1
                           GSIND = (INDX - 1) * NB + 1 
                           CEXCOL = .TRUE.
                        ELSEIF( IWORK(EXAINF+(INDX-1)).EQ.2 ) THEN
                           ATROWS2 = MIN(NB, (N - NB + IROFFA) - 
     $                          (INDX - 2) * NB) - 1
                           GSIND = (INDX - 1) * NB + 2 
                           CEXCOL = .FALSE.
                        ELSEIF( IWORK(EXAINF+(INDX-1)).EQ.3 ) THEN
                           ATROWS2 = NB 
                           GSIND = (INDX - 1) * NB + 2 
                           CEXCOL = .TRUE.
                        END IF
                     END IF
C     
                     CALL INFOG2L( GSIND, GSJ, DESCSA, NPROW, NPCOL, 
     $                             MYROW, MYCOL, LIAT, LJAT, RSRC, 
     $                             CSRC )
                     CALL INFOG2L( GSI, GSIND, DESCSC, NPROW, NPCOL, 
     $                             MYROW, MYCOL, IX, JX, RRSRC, RCSRC )
C
C     Submatrix constuction is omitted most of the times for 
C     deep pipelining
C
                     IF( IDEEP.NE.IDEEPS .OR. JDEEP.NE.JDEEPS )
     $                    GO TO 72
C
                     IF( MYROW.EQ.RSRC.AND.MYCOL.EQ.CSRC ) THEN
                        IF( (RSRC.NE.RRSRC).OR.(CSRC.NE.RCSRC) ) THEN
                           CALL DBEXMAT( CEXCOL, EXCOL, ATROWS2, COLS, 
     $                          LIAT, LJAT, A((LJAS-1)*LLDA+LIAS), 
     $                          LLDA, NBCA, NB, NB, DWORK( EXA ), 
     $                          DWORK( AUP2+(AUPBL-1)*(NB+1)**2 ), 
     $                          NB + 1 )
                           CALL DGESD2D( ICTXT, ATROWS2, COLS, 
     $                          DWORK( AUP2+(AUPBL-1)*(NB+1)**2  ), 
     $                          NB + 1, RRSRC, RCSRC )
                        ELSE
                           CALL DBEXMAT( CEXCOL, EXCOL, ATROWS2, COLS, 
     $                          LIAT, LJAT, A((LJAS-1)*LLDA+LIAS), 
     $                          LLDA, NBCA, NB, NB, DWORK( EXA ), 
     $                          DWORK( AUP2+(KAUP2-1)*(NB+1)**2 ), 
     $                          NB + 1 )
                        END IF
                     END IF
C     
C    Skipped submatrix construction?
C     
 72                  CONTINUE
C 
                     IF( MYROW.EQ.RRSRC .AND. MYCOL.EQ.RCSRC ) THEN
                        CALL DBEXMAT( EXROW, CEXCOL, ROWS, ATROWS2, IX, 
     $                       JX, C((LJCS-1)*LLDC+LICS), LLDC, 
     $                       NBCA, NB, NB, DWORK( EXC ), 
     $                       DWORK(CKJ), NB + 1 )
                       IF( IDEEP.NE.IDEEPS .OR. JDEEP.NE.JDEEPS )
     $                       GO TO 74 
                        IF( (RSRC.NE.RRSRC) .OR. (CSRC.NE.RCSRC) ) THEN
                           CALL DGERV2D( ICTXT, ATROWS2, COLS, 
     $                          DWORK( AUP2+(KAUP2-1)*(NB+1)**2 ), 
     $                          NB + 1, RSRC, CSRC )
                        END IF
 74                     CONTINUE
                        XRWS = IWORK( IXRWS + IIS - I ) 
                        XCLS = IWORK( IXCLS + IIS - I ) 
                        XRIND = IWORK( IXRIND + IIS - I )
                        XCIND = IWORK( IXCIND + IIS - I )
                        CALL DGEMM( 'N', 'T', XRWS, ATROWS2, XCLS, 
     $                       -1.0D+00, DWORK( XIJ1+
     $                       (XCIND-1)*(NB+1)+XRIND-1 ), NB + 1, 
     $                       DWORK( AUP2+(KAUP2-1)*(NB+1)**2+
     $                       (XCIND-1)*(NB+1)), NB + 1, 1.0D+00, 
     $                       DWORK( CKJ+XRIND-1 ), NB + 1 )
                        IF( CSYMM .AND. I.EQ.J .AND. IDEEP.GT.JDEEP )
     $                       CALL DGEMM( 'T', 'T', XCLS, ATROWS2, XRWS, 
     $                       -1.0D+00, DWORK( XIJ1+
     $                       (XCIND-1)*(NB+1)+XRIND-1 ), NB + 1, 
     $                       DWORK( AUP2+(KAUP2-1)*(NB+1)**2+
     $                       (XRIND-1)*(NB+1)), NB + 1, 1.0D+00, 
     $                       DWORK( CKJ+XCIND-1 ), NB + 1 )
                        CALL DUBEXMA( EXROW, CEXCOL, ROWS, ATROWS2, IX, 
     $                       JX, C((LJCS-1)*LLDC+LICS), LLDC, 
     $                       NBCA, NB, NB, DWORK( EXC ), 
     $                       DWORK(CKJ), NB + 1 )
                        IF( AUPBL.NE.1 ) KAUP2 = KAUP2 + 1
                     END IF
 70               CONTINUE
C     
C     In the case of symmetric C and we are below the main block
C     diagonal: update the jth block row of C, that is, do the operation
C     
C     Cjk = Cjk - Xji*Aki^T, where k belongs to [1,j-1]
C     
                  IF( CSYMM .AND. I.GT.J ) THEN
                     IF( MOD( I, 2 ).EQ.0 ) THEN
                        INDXS = 1
                        INDXE = J
                        IF( NB2.EQ.NB ) INDXE = J - 1
                        INDXU = 1
                     ELSE
                        INDXS = J
                        IF( NB2.EQ.NB ) INDXS = J - 1
                        INDXE = 1
                        INDXU = -1
                     END IF
                     KAUP3 = 1 
                     DO 80 INDX = INDXS, INDXE, INDXU
                        IF( INDX.EQ.1 ) THEN
                           IF( IWORK(EXAINF).EQ.0 ) THEN
                              ATCOLS2 = NB - IROFFA
                              GSIND = 1 + IROFFA
                              CEXCOL = .FALSE.       
                           ELSEIF( IWORK(EXAINF).EQ.1 ) THEN
                              ATCOLS2 = NB - IROFFA + 1
                              GSIND = 1 + IROFFA
                              CEXCOL = .TRUE.
                           END IF
                        ELSE
                           IF( IWORK(EXAINF+(INDX-1)).EQ.0 ) THEN
                              ATCOLS2 = MIN(NB, (N - NB + IROFFA) - 
     $                             (INDX - 2) * NB)
                              GSIND = (INDX - 1) * NB + 1
                              CEXCOL = .FALSE.
                           ELSEIF( IWORK(EXAINF+(INDX-1)).EQ.1 ) THEN
                              ATCOLS2 = NB + 1
                              GSIND = (INDX - 1) * NB + 1
                              CEXCOL = .TRUE.
                           ELSEIF( IWORK(EXAINF+(INDX-1)).EQ.2 ) THEN
                              ATCOLS2 = MIN(NB, (N - NB + IROFFA) - 
     $                             (INDX - 2) * NB) - 1
                              GSIND = (INDX - 1) * NB + 2
                              CEXCOL = .FALSE.
                           ELSEIF( IWORK(EXAINF+(INDX-1)).EQ.3 ) THEN
                              ATCOLS2 = NB 
                              GSIND = (INDX - 1) * NB + 2
                              CEXCOL = .TRUE.
                           END IF
                        END IF
C     
                        CALL INFOG2L( GSIND, GSI, DESCSA, NPROW, NPCOL, 
     $                                MYROW, MYCOL, LIAT, LJAT, RSRC, 
     $                                CSRC )
                        CALL INFOG2L( GSJ, GSIND, DESCSC, NPROW, NPCOL,
     $                                MYROW, MYCOL, IX, JX, RRSRC, 
     $                                RCSRC )
C
C     Submatrix constuction is omitted most of the times for 
C     deep pipelining
C
                     IF( IDEEP.NE.IDEEPS .OR. JDEEP.NE.JDEEPS )
     $                    GO TO 82
C
                        IF( MYROW.EQ.RSRC .AND. MYCOL.EQ.CSRC ) THEN
                           IF( (RSRC.NE.RRSRC).OR.(CSRC.NE.RCSRC) ) THEN
                              CALL DBEXMAT( CEXCOL, EXROW, ATCOLS2, 
     $                             ROWS, LIAT, LJAT, 
     $                             A((LJAS-1)*LLDA+LIAS), LLDA, 
     $                             NBCA, NB, NB, DWORK( EXA ), 
     $                             DWORK( AUP3+(AUPBL-1)*(NB+1)**2 ), 
     $                             NB + 1 )
                              CALL DGESD2D( ICTXT, ATCOLS2, ROWS, 
     $                             DWORK( AUP3+(AUPBL-1)*(NB+1)**2 ), 
     $                             NB + 1, RRSRC, RCSRC )
                           ELSE
                              CALL DBEXMAT( CEXCOL, EXROW, ATCOLS2, 
     $                             ROWS, LIAT, LJAT, 
     $                             A((LJAS-1)*LLDA+LIAS), LLDA, 
     $                             NBCA, NB, NB, DWORK( EXA ), 
     $                             DWORK( AUP3+(KAUP3-1)*(NB+1)**2 ), 
     $                             NB + 1 )
                           END IF
                        END IF
C     
C    Skipped submatrix construction?
C     
 82                     CONTINUE
C
                        IF(MYROW.EQ.RRSRC.AND.MYCOL.EQ.RCSRC) THEN
                           CALL DBEXMAT( EXCOL, CEXCOL, COLS, ATCOLS2, 
     $                          IX, JX, C((LJCS-1)*LLDC+LICS), 
     $                          LLDC, NBCA, NB, NB, DWORK(EXC),
     $                          DWORK(CKJ), NB + 1 )
                           IF( IDEEP.NE.IDEEPS .OR. JDEEP.NE.JDEEPS )
     $                       GO TO 84
                           IF( (RSRC.NE.RRSRC).OR.(CSRC.NE.RCSRC) ) THEN
                              CALL DGERV2D( ICTXT, ATCOLS2, ROWS, 
     $                             DWORK(AUP3+(KAUP3-1)*(NB+1)**2), 
     $                             NB + 1, RSRC, CSRC )
                           END IF
 84                        CONTINUE
                           XRWS = IWORK( IXRWS + IIS - I ) 
                           XCLS = IWORK( IXCLS + IIS - I ) 
                           XRIND = IWORK( IXRIND + IIS - I )
                           XCIND = IWORK( IXCIND + IIS - I )
                           CALL DGEMM( 'N', 'T', XCLS, ATCOLS2, XRWS, 
     $                          -1.0D+00, DWORK( XIJT+(XRIND-1)*(NB+1)+
     $                          XCIND-1 ), NB + 1, 
     $                          DWORK( AUP3+(KAUP3-1)*(NB+1)**2+
     $                          (XRIND-1)*(NB+1) ), NB + 1, ONE, 
     $                          DWORK( CKJ+XCIND-1 ), NB + 1 )
                           CALL DUBEXMA( EXCOL, CEXCOL, COLS, ATCOLS2, 
     $                          IX, JX, C((LJCS-1)*LLDC+LICS), 
     $                          LLDC, NBCA, NB, NB, DWORK(EXC), 
     $                          DWORK(CKJ), NB + 1 )
                           IF( AUPBL.NE.1 ) KAUP3 = KAUP3 + 1
                        END IF
 80                  CONTINUE   
                  END IF
C     
C     This branch corresponds to solving A^T * X + X*A = C
C
C     Start updating block column j, i.e., do the operation
C     Ckj = Ckj - Aik^T * Xij, for all k that belongs to [i+1,DBA]
C
               ELSEIF( .NOT. RSIDE ) THEN
                  IF( MOD( J, 2 ).EQ.0 ) THEN
                     INDXS = I + 1 
                     INDXE = DBA
                     INDXU = 1
                  ELSE
                     INDXS = DBA
                     INDXE = I + 1
                     INDXU = -1
                  END IF
                  KAUP1 = 1
                  DO 90 INDX =  INDXS, INDXE, INDXU
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
                     CALL INFOG2L( GSI, GSIND, DESCSA, NPROW, NPCOL,
     $                    MYROW, MYCOL, LIA, LJA, RSRC, CSRC )
                     CALL INFOG2L( GSIND, GSJ, DESCSC, NPROW, NPCOL, 
     $                    MYROW, MYCOL, IX, JX, RRSRC, RCSRC )
C     
C     Submatrix constuction is omitted most of the times for 
C     deep pipelining
C     
                     IF( IDEEP.NE.IDEEPS .OR. JDEEP.NE.JDEEPS )
     $                    GO TO 92
C     
                     IF( MYROW.EQ.RSRC .AND. MYCOL.EQ.CSRC ) THEN
                        IF( (RSRC.NE.RRSRC) .OR. (CSRC.NE.RCSRC) ) THEN
                           CALL DBEXMAT( EXROW, CEXROW, ROWS, ACOLS2, 
     $                          LIA, LJA, A((LJAS-1)*LLDA+LIAS), 
     $                          LLDA, NBCA, NB, NB, DWORK( EXA ), 
     $                          DWORK( AUP1+(AUPBL-1)*(NB+1)**2 ), 
     $                          NB + 1 )
                           CALL DGESD2D( ICTXT, ROWS, ACOLS2, 
     $                          DWORK( AUP1+(AUPBL-1)*(NB+1)**2 ), 
     $                          NB + 1, RRSRC, RCSRC )
                        ELSE
                           CALL DBEXMAT( EXROW, CEXROW, ROWS, ACOLS2, 
     $                          LIA, LJA, A((LJAS-1)*LLDA+LIAS), 
     $                          LLDA, NBCA, NB, NB, DWORK( EXA ), 
     $                          DWORK( AUP1+(KAUP1-1)*(NB+1)**2 ), 
     $                          NB + 1 )
                        END IF
                     END IF    
C     
C     Skipped submatrix construction?
C     
 92                  CONTINUE
C     
                     IF(MYROW.EQ.RRSRC .AND. MYCOL.EQ.RCSRC) THEN
                        CALL DBEXMAT( CEXROW, EXCOL, ACOLS2, COLS, IX, 
     $                       JX, C((LJCS-1)*LLDC+LICS), LLDC, 
     $                       NBCA, NB, NB, DWORK( EXC ), 
     $                       DWORK(CKJ), NB + 1 )
                        IF( IDEEP.NE.IDEEPS .OR. JDEEP.NE.JDEEPS )
     $                       GO TO 94
                        IF( (RSRC.NE.RRSRC) .OR. (CSRC.NE.RCSRC) ) THEN
                           CALL DGERV2D( ICTXT, ROWS, ACOLS2, 
     $                          DWORK( AUP1+(KAUP1-1)*(NB+1)**2 ), 
     $                          NB + 1, RSRC, CSRC )
                        END IF
 94                     CONTINUE
                        XRWS = IWORK( IXRWS + IIS - I ) 
                        XCLS = IWORK( IXCLS + IIS - I ) 
                        XRIND = IWORK( IXRIND + IIS - I )
                        XCIND = IWORK( IXCIND + IIS - I )
                        CALL DGEMM( 'T', 'N', ACOLS2, XCLS, XRWS, 
     $                       -1.0D+00,
     $                       DWORK( AUP1+(KAUP1-1)*(NB+1)**2+
     $                       XRIND-1 ), NB + 1, 
     $                       DWORK( XIJ2+(XCIND-1)*(NB+1)+XRIND-1 ),
     $                       NB + 1, ONE, DWORK(CKJ+(XCIND-1)*(NB+1)), 
     $                       NB + 1 )
                        IF( CSYMM .AND. I.EQ.J .AND. IDEEP.GT.JDEEP )
     $                       CALL DGEMM( 'T', 'T', ACOLS2, XRWS, XCLS, 
     $                       -1.0D+00,
     $                       DWORK( AUP1+(KAUP1-1)*(NB+1)**2+
     $                       XCIND-1 ), NB + 1, 
     $                       DWORK( XIJ2+(XCIND-1)*(NB+1)+XRIND-1 ),
     $                       NB + 1, ONE, DWORK(CKJ+(XRIND-1)*(NB+1)), 
     $                       NB + 1 )
                        CALL DUBEXMA( CEXROW, EXCOL, ACOLS2, COLS, IX, 
     $                       JX, C((LJCS-1)*LLDC+LICS), LLDC, 
     $                       NBCA, NB, NB, DWORK( EXC ), 
     $                       DWORK(CKJ), NB + 1 )
                        IF( AUPBL.NE.1 ) KAUP1 = KAUP1 + 1
                     END IF
 90               CONTINUE
C     
C     Now update block row i of C, that is do the operation
C     
C     Cik = Cik - Xij * Ajk, for all k that belongs to [j+1,IEND].
C     
C     If the matrix C is symmetric we do the symmetric update
C     and we shall update Cik, where k=i, we do a SYR2K-update, 
C     namely the operation
C     
C     Cik = Cik - Xij * Ajk - Ajk^T * Xji = (Ajk^T) * Xij^T -
C                                            Xij * (Ajk^T)^T
C     
C     To do this as a regular SYR2K-operation we will have to
C     transpose Ajk explicitely before calling the SYR2K Blas
C     routine.
C     
                  IF( MOD( I, 2 ).EQ.0 ) THEN
                     INDXS = J + 1
                     INDXE = DBA
                     IF( CSYMM ) INDXE = I
                     INDXU = 1
                  ELSE
                     INDXS = DBA
                     IF( CSYMM ) INDXS = I
                     INDXE = J + 1
                     INDXU = -1
                  END IF
                  KAUP2 = 1
                  DO 100 INDX =  INDXS, INDXE, INDXU
                     IF( IWORK(EXAINF+(INDX-1)).EQ.0 ) THEN
                        ACOLS2 = MIN(NB, (N - NB + IROFFA) - 
     $                       (INDX - 2) * NB)
                        GSIND = (INDX - 1) * NB + 1
                        CEXCOL = .FALSE.
                     ELSEIF( IWORK(EXAINF+(INDX-1)).EQ.1 ) THEN
                        ACOLS2 = NB + 1
                        GSIND = (INDX - 1) * NB + 1
                        CEXCOL = .TRUE.
                     ELSEIF( IWORK(EXAINF+(INDX-1)).EQ.2 ) THEN
                        ACOLS2 = MIN(NB, (N - NB + IROFFA) - 
     $                       (INDX - 2) * NB) - 1
                        GSIND = (INDX - 1) * NB + 2
                        CEXCOL = .FALSE.
                     ELSEIF( IWORK(EXAINF+(INDX-1)).EQ.3 ) THEN
                        ACOLS2 = NB 
                        GSIND = (INDX - 1) * NB + 2
                        CEXCOL = .TRUE.
                     END IF
C     
                     CALL INFOG2L( GSJ, GSIND, DESCSA, NPROW, NPCOL, 
     $                             MYROW, MYCOL, LIA, LJA, RSRC, CSRC )
                     CALL INFOG2L( GSI, GSIND, DESCSC, NPROW, NPCOL, 
     $                             MYROW, MYCOL, IX, JX, RRSRC, RCSRC )
C
C     Submatrix constuction is omitted most of the times for 
C     deep pipelining
C
                     IF( IDEEP.NE.IDEEPS .OR. JDEEP.NE.JDEEPS )
     $                    GO TO 102                     
C
                     IF( MYROW.EQ.RSRC .AND. MYCOL.EQ.CSRC ) THEN
                        IF( (RSRC.NE.RRSRC) .OR. (CSRC.NE.RCSRC) ) THEN
                           CALL DBEXMAT( EXCOL, CEXCOL, COLS, ACOLS2, 
     $                          LIA, LJA, A((LJAS-1)*LLDA+LIAS), LLDA, 
     $                          NBCA, NB, NB, DWORK( EXA ), 
     $                          DWORK( AUP2+(AUPBL-1)*(NB+1)**2 ), 
     $                          NB + 1 )
                           CALL DGESD2D( ICTXT, COLS, ACOLS2, 
     $                          DWORK( AUP2+(AUPBL-1)*(NB+1)**2 ), 
     $                          NB + 1, RRSRC, RCSRC )
                        ELSE
                           CALL DBEXMAT( EXCOL, CEXCOL, COLS, ACOLS2, 
     $                          LIA, LJA, A((LJAS-1)*LLDA+LIAS), LLDA, 
     $                          NBCA, NB, NB, DWORK( EXA ), 
     $                          DWORK( AUP2+(KAUP2-1)*(NB+1)**2 ), 
     $                          NB + 1 )
                        END IF
                     END IF
C     
C    Skipped submatrix construction?
C     
 102                 CONTINUE
C
                     IF( MYROW.EQ.RRSRC .AND. MYCOL.EQ.RCSRC ) THEN
                        CALL DBEXMAT( EXROW, CEXCOL, ROWS, ACOLS2, IX, 
     $                       JX, C((LJCS-1)*LLDC+LICS), LLDC, 
     $                       NBCA, NB, NB, DWORK( EXC ), 
     $                       DWORK(CKJ), NB + 1 )
                        IF( IDEEP.NE.IDEEPS .OR. JDEEP.NE.JDEEPS )
     $                       GO TO 104
                        IF( (RSRC.NE.RRSRC) .OR. (CSRC.NE.RCSRC) ) THEN
                           CALL DGERV2D( ICTXT, COLS, ACOLS2, 
     $                          DWORK( AUP2+(KAUP2-1)*(NB+1)**2 ), 
     $                          NB + 1, RSRC, CSRC )
                        END IF
 104                    CONTINUE
                        XRWS = IWORK( IXRWS + IIS - I ) 
                        XCLS = IWORK( IXCLS + IIS - I ) 
                        XRIND = IWORK( IXRIND + IIS - I )
                        XCIND = IWORK( IXCIND + IIS - I )
                        IF( NB2.EQ.NB .AND. CSYMM .AND. INDX.EQ.I ) THEN
                           CALL DLATCPY( 'All', XCLS, XRWS, 
     $                          DWORK( AUP2+(KAUP2-1)*(NB+1)**2+
     $                          XCIND-1 ), NB + 1, DWORK( IPW ), 
     $                          XRWS )
                           CALL DSYR2K( 'U', 'N', XRWS, XCLS, -1.0D+00, 
     $                          DWORK( IPW ), XRWS, 
     $                          DWORK(XIJ1+(XCIND-1)*(NB+1)+XRIND-1),
     $                          NB + 1, 1.0D+00, 
     $                          DWORK(CKJ+(XRIND-1)*(NB+1)+XRIND-1), 
     $                          NB+1 )
                           CALL DLATCPY( 'Upper', XRWS-1, XRWS-1, 
     $                          DWORK(CKJ+(XRIND-1)*(NB+1)+XRIND+NB), 
     $                          NB + 1,
     $                          DWORK(CKJ+(XRIND-1)*(NB+1)+XRIND), 
     $                          NB + 1 )
                        ELSE
                           CALL DGEMM( 'N', 'N', XRWS, ACOLS2, XCLS, 
     $                          -1.0D+00, 
     $                          DWORK( XIJ1+(XCIND-1)*(NB+1)+XRIND-1 ), 
     $                          NB+1, DWORK( AUP2+(KAUP2-1)*(NB+1)**2+
     $                          XCIND-1  ), NB + 1, 1.0D+00, 
     $                          DWORK(CKJ+XRIND-1), NB + 1 )
                           IF( CSYMM .AND. I.EQ.J .AND. IDEEP.GT.JDEEP )
     $                          CALL DGEMM( 'T', 'N', XCLS, ACOLS2, 
     $                          XRWS, -1.0D+00, 
     $                          DWORK( XIJ1+(XCIND-1)*(NB+1)+XRIND-1 ), 
     $                          NB+1, DWORK( AUP2+(KAUP2-1)*(NB+1)**2+
     $                          XRIND-1  ), NB + 1, 1.0D+00, 
     $                          DWORK(CKJ+XCIND-1), NB + 1 )  
                        END IF
                        CALL DUBEXMA( EXROW, CEXCOL, ROWS, ACOLS2, IX, 
     $                       JX, C((LJCS-1)*LLDC+LICS), LLDC, NBCA, NB, 
     $                       NB, DWORK( EXC ), DWORK(CKJ), NB + 1 )
                        IF( AUPBL.NE.1 ) KAUP2 = KAUP2 + 1
                     END IF 
 100              CONTINUE
C     
C     In the case of symmetric C and we are below the main block
C     diagonal: update the ith block column of C, that is, do the 
C     operation
C     
C     Cki = Cki - Ajk^T * Xji, where k belongs to [i+1,DBA].
C     
                  IF( CSYMM .AND. I.GT.J ) THEN
                     IF( MOD( J, 2 ).EQ.0 ) THEN
                        INDXS = I
                        IF( NB2.EQ.NB ) INDXS = I + 1
                        INDXE = DBA
                        INDXU = 1
                     ELSE
                        INDXS = DBA
                        INDXE = I 
                        IF( NB2.EQ.NB ) INDXE = I + 1
                        INDXU = -1
                     END IF
                     KAUP3 = 1
                     DO 110 INDX =  INDXS, INDXE, INDXU 
                        IF( IWORK(EXAINF+(INDX-1)).EQ.0 ) THEN
                           ATCOLS2 = MIN(NB, (N - NB + IROFFA) - 
     $                          (INDX - 2) * NB)
                           GSIND = (INDX - 1) * NB + 1
                           CEXROW = .FALSE.
                        ELSEIF( IWORK(EXAINF+(INDX-1)).EQ.1 ) THEN
                           ATCOLS2 = NB + 1
                           GSIND = (INDX - 1) * NB + 1
                           CEXROW = .TRUE.
                        ELSEIF( IWORK(EXAINF+(INDX-1)).EQ.2 ) THEN
                           ATCOLS2 = MIN(NB, (N - NB + IROFFA) - 
     $                          (INDX - 2) * NB) - 1
                           GSIND = (INDX - 1) * NB + 2
                           CEXROW = .FALSE.
                        ELSEIF( IWORK(EXAINF+(INDX-1)).EQ.3 ) THEN
                           ATCOLS2 = NB 
                           GSIND = (INDX - 1) * NB + 2
                           CEXROW = .TRUE.
                        END IF
C     
                        CALL INFOG2L( GSJ, GSIND, DESCSA, NPROW, NPCOL,
     $                                MYROW, MYCOL, LIAT, LJAT, RSRC, 
     $                                CSRC )
                        CALL INFOG2L( GSIND, GSI, DESCSC, NPROW, NPCOL,
     $                                MYROW, MYCOL, IX, JX, RRSRC, 
     $                                RCSRC )
C
C     Submatrix constuction is omitted most of the times for 
C     deep pipelining
C
                        IF( IDEEP.NE.IDEEPS .OR. JDEEP.NE.JDEEPS )
     $                       GO TO 112
                        IF( MYROW.EQ.RSRC .AND. MYCOL.EQ.CSRC ) THEN
                           IF( (RSRC.NE.RRSRC).OR.(CSRC.NE.RCSRC) ) THEN
                              CALL DBEXMAT( EXCOL, CEXROW, COLS, 
     $                             ATCOLS2, LIAT, LJAT, 
     $                             A((LJAS-1)*LLDA+LIAS), LLDA, NBCA, 
     $                             NB, NB, DWORK( EXA ), 
     $                             DWORK( AUP3+(AUPBL-1)*(NB+1)**2 ), 
     $                             NB + 1 )
                              CALL DGESD2D( ICTXT, COLS, ATCOLS2, 
     $                             DWORK( AUP3+(AUPBL-1)*(NB+1)**2 ), 
     $                             NB + 1, RRSRC, RCSRC )
                           ELSE
                              CALL DBEXMAT( EXCOL, CEXROW, COLS, 
     $                             ATCOLS2, LIAT, LJAT, 
     $                             A((LJAS-1)*LLDA+LIAS), LLDA, NBCA, 
     $                             NB, NB, DWORK( EXA ), 
     $                             DWORK( AUP3+(KAUP3-1)*(NB+1)**2 ), 
     $                             NB + 1 )
                           END IF
                        END IF
C     
C    Skipped submatrix construction?
C     
 112                    CONTINUE
C     
                        IF( MYROW.EQ.RRSRC .AND. MYCOL.EQ.RCSRC ) THEN
                           CALL DBEXMAT( CEXROW, EXROW, ATCOLS2, ROWS, 
     $                          IX, JX, C((LJCS-1)*LLDC+LICS), LLDC, 
     $                          NBCA, NB, NB, DWORK( EXC ), DWORK(CKJ), 
     $                          NB + 1 )
                           IF( IDEEP.NE.IDEEPS .OR. JDEEP.NE.JDEEPS )
     $                          GO TO 124
                           IF( (RSRC.NE.RRSRC).OR.(CSRC.NE.RCSRC) ) THEN
                              CALL DGERV2D( ICTXT, COLS, ATCOLS2, 
     $                             DWORK(AUP3+(KAUP3-1)*(NB+1)**2), 
     $                             NB + 1, RSRC, CSRC )
                           END IF
 124                       CONTINUE
                           XRWS = IWORK( IXRWS + IIS - I ) 
                           XCLS = IWORK( IXCLS + IIS - I ) 
                           XRIND = IWORK( IXRIND + IIS - I )
                           XCIND = IWORK( IXCIND + IIS - I )
                           CALL DGEMM( 'T', 'N', ATCOLS2, XRWS, XCLS, 
     $                          -1.0D+00, 
     $                          DWORK( AUP3+(KAUP3-1)*(NB+1)**2+
     $                          XCIND-1 ), NB + 1, 
     $                          DWORK( XIJT+(XRIND-1)*(NB+1)+XCIND-1 ), 
     $                          NB + 1, ONE, 
     $                          DWORK( CKJ+(XRIND-1)*(NB+1) ), NB + 1 )
                           CALL DUBEXMA( CEXROW, EXROW, ATCOLS2, ROWS, 
     $                          IX, JX, C((LJCS-1)*LLDC+LICS), LLDC, 
     $                          NBCA, NB, NB, DWORK( EXC ), DWORK(CKJ), 
     $                          NB + 1 )
                           IF( AUPBL.NE.1 ) KAUP3 = KAUP3 + 1
                        END IF
 110                 CONTINUE   
                  END IF
               END IF 
C     
 55            CONTINUE
C     
C     Update inner loop variable j (block column index)
C     
               J = J + 1 
C     
 50         CONTINUE
 24         CONTINUE
 26         CONTINUE
 22         CONTINUE
 20      CONTINUE
C     
C     Update outer loop variables (where to go the next roll) - this we
C     do differently depending on the symmetry of the problem
C
         IF( RSIDE ) THEN
            JS = MAX( 1, JS - 1 )
            IF( CSYMM .AND. MOD( K, 2 ).EQ.0 ) THEN
               IE = MAX( 1, IE - 1 )
            ELSEIF( .NOT. CSYMM ) THEN
               IE = MAX( 1, IE - 1 )
            END IF
         ELSEIF( .NOT. RSIDE ) THEN
            IF( K .GE. DBA ) JS = JS + 1
            IS = MIN ( IS + 1, DBA )
         END IF
C     
 10   CONTINUE
C     
C     Before we go on we must back redistributed the elements in C
C     
      IF( AEXT )
     $     CALL PDBCKRD( N, N, C, IC, JC, DESCC, IWORK( EXAINF ), 
     $                   IWORK( EXAINF ), DWORK( EXC ), EXMEMC,
     $                   DWORK( SND ), NB, INFO )
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
            CALL DGAMN2D( ICTXT, 'All', ' ', 1, 1, SCALE, 1, -1, -1, -1,
     $                    -1, -1 )
         END IF
         RETURN
      END IF
C     
      END
C     
C     End of PTRLYCTD
C     
C *** Last line of PTRLYCTD ***
