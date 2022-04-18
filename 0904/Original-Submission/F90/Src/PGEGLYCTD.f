CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C         This file is a part of the SCASY software library.           C
C         See http://www.cs.umu.se/parallel/research/scasy             C
C         Contributors: Robert Granat and Bo Kågström.                 C
C         SCASY Version 0.10,    March 31, 2006.                       C
C         Copyright 2006, Umeå University, Sweden.                     C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE PGEGLYCTD( JOB, SYMM, OP, AESCHR, N, A, IA, JA, DESCA, 
     $                      E, IE, JE, DESCE, C, IC, JC, DESCC, NB2,
     $                      DWORK, LDWORK, IWORK, LIWORK, NOEXSY, SCALE, 
     $                      INFO )
C
C  -- ScaLAPACK-style routine (preliminary version ) --
C     HPC2N and Department of Computing Science,
C     University of Umeå, Sweden.
C     Written by Robert Granat, (granat@cs.umu.se)
C     January 22, 2006.
C
      IMPLICIT NONE
C
C     .. Scalar Arguments ..
      CHARACTER*1          JOB, SYMM, OP, AESCHR
      INTEGER              N, IA, JA, IE, JE, IC, JC, LDWORK, LIWORK, 
     $                     NOEXSY, INFO, NB2
      DOUBLE PRECISION     SCALE
C     ..
C     .. Array Arguments ..
      INTEGER              DESCA( * ), DESCE( * ), DESCC( * ), 
     $                     IWORK( * )
      DOUBLE PRECISION     A( * ), E( * ), C( * ), DWORK( * )
C 
C  Purpose
C  =======
C  The subroutine solves the general generalized continuous-time
C  Lyapunov Equation (GLYCT) 
C
C     op(sub(A)) * X * op(sub(E)^T) + op(sub(E)) * X * op(sub(A)^T) = C,
C
C  where sub(A) = A(IA:IA+M-1,JA:JA+M-1), sub(E) = E(IE:IE+M-1,JE:JE+N-1)
C  and the solution X which overwrites sub(C) = C(IC:IC+M-1,JC:JC+M-1) 
C  are N-by-N matrices.
C  
C  The notation op(_) means the transpose or non-transpose of a matrix. 
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
C  JOB       (global input) CHARACTER*1
C            With JOB = 'S', the equation is solved. For JOB = 'R',
C            only the reduction step regarding the matrix pair (A,E)
C            is performed.
C
C  SYMM      (global input) CHARACTER*1
C            With JOB = 'S':
C              If SYMM = 'S', the matrix C is assumed to be symmetric.
C              If SYMM = 'N', the matrix C is not assumed to be symmetric.
C            Otherwise, SYMM is not referenced.
C
C  OP        (global input) CHARACTER*1
C            With JOB = 'S':
C              If OP = 'N', then op(A) = A and op(E) = E.
C              If OP = 'T', then op(A) = A^T and op(E) = E^T.
C            Otherwise, OP is not referenced.
C
C  AESCHR   (global input) CHARACTER*1
C            If AESCHR = 'S', then the matrix pair (A,E) is supposed to 
C            be in generalized Schur form. No reduction to generalized
C            Schur form of this matrix pair will be done.
C            If AESCHR = 'N', then the matrix pair (A,E) is not in 
C            generalized triangular form and the full reduction to 
C            generalized Schur form will be carried out.
C
C  Input/Output parameters
C
C  N         (global input) INTEGER
C            The number of rows and columns of the global distributed
C            matrices A, E and C.
C
C  A         (local input) DOUBLE PRECISION array 
C            Array of dimension (LLD_A,LOCc(N)). Contains the local
C            pieces of the global distributed matrix A. On exit, with
C            AESCHR = 'N', it contains the local pieces of the 
C            generalized Schur form of (A,E) corresponding to the
C            quasi-triangular matrix A.
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
C  E         (local input) DOUBLE PRECISION array 
C            Array of dimension (LLD_E,LOCc(N)). Contains the local
C            pieces of the global distributed matrix E. On exit, with
C            AESCHR = 'N', it contains the local pieces of the 
C            generalized Schur form of (A,E) corresponding to the
C            upper triangular matrix E.
C
C  IE        (global input) INTEGER
C            Row start index for sub(E), i.e., the submatrix to operate 
C            on. IE = IA must hold. 
C
C  JE        (global input) INTEGER
C            Column start index for sub(E), i.e., the submatrix to 
C            operate on. JE = JA must hold. 
C
C  DESCE     (global and local input) INTEGER array of dimension DLEN_.
C            The array descriptor for the global distributed matrix E.
C
C  C         (local input/local output) DOUBLE PRECISION array 
C            Array of dimension (LLD_C,LOCc(N)). 
C            With JOB = 'S': 
C            On entry C contains the local pieces of the global 
C            distributed matrix C . On exit, it contains the local 
C            pieces of the global distributed solution X.
C            Otherwise, C is not referenced. 
C
C  IC        (global input) INTEGER
C            With JOB='S': Row start index for sub(C), i.e., the 
C            submatrix to operate on. MOD(IC,MB_A) = MOD(IA,MB_A) 
C            must hold.
C            Otherwise, IC is not referenced. 
C
C  JC        (global input) INTEGER
C            With JOB = S: Column start index for sub(C), i.e., the 
C            submatrix to operate on. MOD(JC,NB_A) = MOD(JA,NB_A) 
C            must hold.
C            Otherwise, JC is not referenced.  
C
C  DESCC     (global and local input) INTEGER array of dimension DLEN_.
C            With JOB = 'S', the array descriptor for the global 
C            distributed matrix C.
C            Otherwise, DESCC is not referenced.
C
C  NB2       (global input) INTEGER 
C            Internal blocking factor for pipelining of subsolutions
C            for updates of the matrix C in PTRGLYCTD (see the references 
C            for details).
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
C            LIWORK >= DBA + 2*MB_A + 8, where 
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
C            With JOB = 'S':
C            When solving the triangular problem in PTRGLYCTD it is possible
C            that we have to extend some subsystems to not lose any data
C            from some 2x2 block of conjugate pairs of eigenvalues. NOEXSY 
C            helps us to keep track of the number of such extensions. See 
C            the references for details.
C            Otherwise, NOEXSY is not referenced.
C
C  SCALE     (global output) DOUBLE PRECISION
C            With JOB = 'S':
C            A scale factor, usually 1.0. A scale factor < 1.0 means
C            the solution may have overflowed. See INFO for details.
C            Otherwise, SCALE is not referenced.
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
C             If INFO = 2, some eigenvalue of (A,E) is equal to or 
C             nearly equal to the reciprocal of an eigenvalue of (A,E).
C             Perturbed values was used to solve the equation (but 
C             (A,E) is unchanged).
C             If INFO = 3,  the problem is badly scaled - C should have 
C             been scaled a factor SCALE before calling this routine too 
C             guarantee an overflow free solution, i.e., the solution
C             may well have overflowed.
C
C  Method
C  ======
C  This subroutine implements a parallel version of the Bartels-Stewart
C  method for solving the general generalized continuous-time Lyapunov 
C  equation (see the references for details).
C
C  Additional requirements
C  =======================
C  (A,E) and C be distributed using the same blocking factor in each 
C  direction.  
C
C  The matrix pair (A,E) must be aligned internally.
C
C  Limitations
C  ===========
C  In contrary to SLICOTs SG03AY this routine do not scale against 
C  overflow in the solution. See SCALE and INFO.
C
C  For AESCHR = 'N' it is not possible to work with submatrices, i.e., 
C  AESCHR = 'N' implies IA = IE = JA = JE = 1. This limitation comes from 
C  the fact that ScaLAPACKs PDLAHQR does not work on submatrices - this 
C  might be changed in a future release.
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
C  Known issues
C  ============
C  Experienced hung problems outside this routine when PDORGQR was used 
C  with more that one cpu. Temporary fix: use PDORMQR instead. The 
C  choice is controlled by logical parameter GLYCT_PDORGQR, see below.
C
C  Keywords
C  ========
C  QZ-algorithm, generalized Schur form, generalized continuous-time 
C  Lyapunov equation, Hessenberg-triangular transformation, 
C  PDGEMM-updates.
C
C  =====================================================================
C
C     ..Parameters ..
      INTEGER            BLOCK_CYCLIC_2D, CSRC_, CTXT_, DLEN_, DTYPE_,
     $                   LLD_, MB_, M_, NB_, N_, RSRC_, NINARG,
     $                   WRKSPQ, MAXBULG
      LOGICAL            GLYCT_PDORGQR
      DOUBLE PRECISION   ZERO, ONE, MONE
      CHARACTER*1        CPYALL
      PARAMETER          ( BLOCK_CYCLIC_2D = 1, DLEN_ = 9, DTYPE_ = 1,
     $                     CTXT_ = 2, M_ = 3, N_ = 4, MB_ = 5, NB_ = 6,
     $                     RSRC_ = 7, CSRC_ = 8, LLD_ = 9, NINARG = 17,
     $                     WRKSPQ = -1, MAXBULG = 1, 
     $                     GLYCT_PDORGQR = .FALSE., ZERO = 0.0D0, 
     $                     ONE = 1.0D+0, MONE = -1.0D+0, CPYALL = 'A') 
C     .. 
C     .. Local Scalars ..
      LOGICAL            SOLVE, LQUERY, SCHRAE, CSYMM, RSIDE
      DOUBLE PRECISION   T1, T2, T22, T3
      INTEGER            ICTXT, MYCOL, MYROW, NPCOL, NPROW, NROWSAE,
     $                   NCOLSAE, SIZEAE, NB, LW1, LW2, LINFO, LW3,
     $                   LWTRI, LSWORK, WRK, DBAE, EIGR, SIZEE,
     $                   BETA, TAUE, J, GJ, I, GI, EIGI, IROFFC,
     $                   ICOFFC, LIC, LJC, CRSRC, CCSRC, IWRK,
     $                   NCOLSC, NROWSC, BULGES
C     ..
C     .. Local Arrays and Local Pointers ..
      INTEGER             Q, Z, CCOPY, TAUA, WR, WI, SWORK,
     $                    IDUM1( 1 ), IDUM2( 1 ), DESCQZ( DLEN_ ),
     $                    DESCCC( DLEN_ )
      DOUBLE PRECISION    DPDUM( 1 )
C     ..
C     .. External Subroutines ..
C     
      EXTERNAL           BLACS_GRIDINFO, CHK1MAT, PCHK2MAT, PXERBLA,
     $                   PDGEMM, PDLACPY, PDHEXTR, PTRGLYCTD, DESCINIT,
     $                   PDGGHRD, PDHGEQZ, PDGEQRF, PDORMQR, PDORGQR
C     .. 
C     .. External Functions ..
      DOUBLE PRECISION   MPI_WTIME
      INTEGER            NUMROC, ICEIL, INDXL2G
      LOGICAL            LSAME
      EXTERNAL           NUMROC, ICEIL, LSAME, MPI_WTIME, INDXL2G
      
C     .. Intrinsic Functions ..
      INTRINSIC MAX, MIN, MOD, INT, REAL
C     ..
C     .. Executable Statements ..
      ICTXT = DESCA( CTXT_ )
      CALL BLACS_GRIDINFO( ICTXT, NPROW, NPCOL, MYROW, MYCOL )
C     
C     Test the input parameters
C     
C     See first if we have a valid context
C
      INFO = 0
      IF( NPROW.EQ.-1 ) THEN
         INFO = 1
      END IF
C
C     Check JOB
C
      IF( .NOT.( LSAME( JOB, 'R' ) .OR. LSAME( JOB, 'S' ) ) ) THEN
         INFO = -1
      ELSE
         IF( LSAME( JOB, 'S' ) ) THEN
            SOLVE = .TRUE.
         ELSE
            SOLVE = .FALSE.
         END IF
      END IF
C     
C     Check dimensions
C     
      IF( INFO.EQ.0 ) THEN
         CALL CHK1MAT( N, 5, N, 5, IA, JA, DESCA, 9, INFO )
         IF( INFO.EQ.0 )
     $        CALL CHK1MAT( N, 5, N, 5, IE, JE, DESCE, 13, INFO )
         IF( INFO.EQ.0 ) 
     $        CALL CHK1MAT( N, 5, N, 5, IC, JC, DESCC, 17, INFO )
      END IF
C     
C     Check the blocking sizes for equivalence and being not to small
C     
      IF( INFO.EQ.0 ) THEN
         IF( DESCA( MB_ ).NE.DESCA( NB_ ) ) INFO = -(100*9 + MB_)
      END IF
      IF( INFO.EQ.0 ) THEN
         IF( DESCE( MB_ ).NE.DESCE( MB_ ) ) INFO = -(100*13 + MB_)
         IF( DESCE( MB_ ).NE.DESCA( MB_ ) ) INFO = -(100*13 + MB_)
      END IF
      IF( INFO.EQ.0 .AND. SOLVE ) THEN
         IF( DESCC( MB_ ).NE.DESCC( NB_ ) ) INFO = -(100*17 + MB_)
         IF( DESCC( MB_ ).NE.DESCA( MB_ ) ) INFO = -(100*17 + MB_)
      END IF
C
      IF( INFO.EQ.0 ) THEN
         IF( N.NE.DESCA( MB_ ) .AND. DESCA( MB_ ).LT.9 .AND. 
     $        LSAME( AESCHR, 'N' ) )
     $        INFO = -(100*9 + MB_)
      END IF
C
C     Check the value of SYMM
C     
      IF( INFO.EQ.0 .AND. SOLVE ) THEN
         IF( .NOT.( LSAME( SYMM, 'N' ) .OR. LSAME( SYMM, 'S' ))) 
     $        THEN
            INFO = -2
         ELSE
            IF( LSAME( SYMM, 'S' ) ) THEN
               CSYMM = .TRUE.
            ELSE
               CSYMM = .FALSE.
            END IF
         END IF
      END IF
C
C     Check alignment, part 1: check RSRC_ and CSRC_
C
      IF( INFO.EQ.0 ) THEN
         IF( DESCE( RSRC_ ) .NE. DESCA( RSRC_ ) ) 
     $        INFO = -(100*13 + RSRC_)
         IF( DESCE( CSRC_ ) .NE. DESCA( CSRC_ ) ) 
     $        INFO = -(100*17 + CSRC_)
      END IF
C
C     Check the values of SCHRAE and set the logical value SCHRAE
C 
      IF( INFO.EQ.0 ) THEN
         IF( .NOT.( LSAME( AESCHR, 'N' ) .OR. LSAME( AESCHR, 'S' ))) 
     $        THEN
            INFO = -3
         ELSE
            IF( LSAME( AESCHR, 'S' ) ) THEN
               SCHRAE = .TRUE.
            ELSE
               SCHRAE = .FALSE.
            END IF
         END IF
      END IF
C
C     Check the values of OP 
C
      IF( INFO.EQ.0 .AND. SOLVE ) THEN
         IF( .NOT.( LSAME( OP, 'N' ) .OR. LSAME( OP, 'T' ))) 
     $        THEN
            INFO = -4
         ELSE
            RSIDE = LSAME( OP, 'N' )
         END IF
      END IF
C
C     Check alignment, part 2: check IA, JA, IE, JE, IC, JC.
C
      IF( INFO.EQ.0 ) THEN
         IF( IA.LT.1 .OR. ( .NOT. SCHRAE .AND. IA.GT.1 ) ) INFO = -7
         IF( JA.NE.IA ) INFO = -8
         IF( IE .NE. IA ) INFO = -11
         IF( JE .NE. JA ) INFO = -12
      END IF
      IF( INFO.EQ.0 .AND. SOLVE ) THEN
         IF( MOD(IC,DESCC(MB_)) .NE. MOD(IA,DESCA(MB_)) ) INFO = -15
         IF( MOD(JC,DESCC(NB_)) .NE. MOD(JA,DESCA(NB_)) ) INFO = -16
      END IF
C
C     Check if workspace query
C
      LQUERY = LDWORK.EQ.-1 .OR. LIWORK.EQ.-1
C     
C     Test working space
C     
C     Check the number of rows and columns of the smallest submatrix
C     of C including sub(C) that conform with the ScaLAPACK conventions and 
C     compute the local number of rows and columns of (A,E).
C
      IF( INFO.EQ.0 .OR. LQUERY ) THEN 
         IROFFC = MOD( IC - 1, DESCC(MB_) ) 
         ICOFFC = MOD( JC - 1, DESCC(NB_) )
         CALL INFOG2L( IC, JC, DESCC, NPROW, NPCOL, MYROW, MYCOL,
     $                 LIC, LJC, CRSRC, CCSRC )
C         
         NROWSAE = NUMROC( DESCA(M_), DESCA(MB_), MYROW, DESCA(RSRC_), 
     $                     NPROW )
         NROWSC  = NUMROC( N + IROFFC, DESCE(MB_), MYROW, CRSRC, NPROW )
         NCOLSAE = NUMROC( DESCA(N_), DESCA(NB_), MYCOL, DESCA(CSRC_), 
     $                     NPCOL )
         NCOLSC  = NUMROC( N + ICOFFC, DESCE(NB_), MYCOL, CCSRC, NPCOL )
C
C     Initialize matrix descriptor for orthogonal transformations
C     Q1, Z1, and the copy of sub(C) and 
C
         IF( .NOT. SCHRAE .AND. SOLVE ) THEN
            CALL DESCINIT( DESCQZ, DESCA(M_), DESCA(N_), DESCA(MB_), 
     $                     DESCA(NB_), DESCA(RSRC_), DESCA(CSRC_), 
     $                     ICTXT, MAX( 1, NROWSAE ), INFO )
            CALL DESCINIT( DESCCC, N + IROFFC, N + ICOFFC, DESCC(MB_),
     $                     DESCC(NB_), CRSRC, CCSRC, ICTXT,
     $                     MAX( 1, NROWSC ), INFO ) 
         END IF         
C     
C     Compute the sizes of the different needed DP workspaces
C     
         SIZEAE = NROWSAE * NCOLSAE 
C     
C     For the LSWORK we need to check the workspace needed by work-
C     space quiries to each individual routines. First, for PDGGHRD...
C     
         IF( .NOT. SCHRAE ) THEN
            CALL PDGGHRD( 'Vectors', 'Vectors', N, IA, IA+N-1, A, DESCA, 
     $           E, DESCE, A, DESCA, E, DESCE, DPDUM, -1, LINFO )
            LW1 = MAX( LW1, INT( DPDUM(1) ) )
         ELSE
            LW1 = 0
         END IF
C     
C     Then for PDHGEQZ...
C     
         IF( .NOT. SCHRAE ) THEN
            BULGES = MAX( 1, MIN( MAXBULG, ICEIL( N, DESCA(MB_) ) - 1 ))
            CALL PDHGEQZ( 'Schur', 'Vectors', 'Vectors', N, IA, IA+N-1, 
     $           A, DESCA, E, DESCE, DWORK, DWORK, DWORK, A, DESCA, 
     $           E, DESCE, IA, IA+N-1, IA, IA+N-1, BULGES, DPDUM, 
     $           -1, LINFO )
            LW2 = MAX( LW2, INT( DPDUM(1) ) )
         ELSE
            LW2 = 0
         END IF
C
C     Then for PDGEQRF, PDORMQR and PDORGQR
C
         IF( .NOT. SCHRAE ) THEN
            CALL PDGEQRF( N, N, C, IC, JC, DESCC, DWORK, DPDUM, -1, 
     $           LINFO )
            LW3 = INT( DPDUM( 1 ) )
         ELSE
            LW3 = 0
         END IF
         IF( .NOT. SCHRAE .AND. SOLVE ) THEN
            CALL PDORMQR( 'Left', 'NoTranspose', N, N, N, E, IE, JE, 
     $           DESCE, DWORK, A, IA, JA, DESCA, DPDUM, -1, LINFO )
            LW3 = MAX( INT( DPDUM( 1 ) ), LW3 )
            CALL PDORMQR( 'Right', 'Transpose', N, N, N, E, IE, JE, 
     $           DESCE, DWORK, A, IA, JA, DESCA, DPDUM, -1, LINFO )
            LW3 = MAX( INT( DPDUM( 1 ) ), LW3 )
            CALL PDORGQR( N, N, N, E, IE, JE, DESCE, DWORK, DPDUM, -1, 
     $           LINFO )
            LW3 = MAX( INT( DPDUM( 1 ) ), LW3 )
         ELSE
            LW3 = 0
         END IF
C 
C     Now compute the workspace needed in PTRGLYCTD using a workspace
C     query as above
C
         IF( SOLVE ) THEN
            CALL PTRGLYCTD( CSYMM, OP, N, A, IA, JA, DESCA, E, IE, JE, 
     $           DESCE, C, IC, JC, DESCC, NB2, DPDUM, -1, IWORK, 
     $           LIWORK, IDUM1, SCALE, LINFO )
            LWTRI = INT( DPDUM( 1 ) )
            IWRK = IDUM1( 1 )
         ELSE
            IWRK = 0
            LWTRI = 0
         END IF
C     
C     Take now the maximum 
C     
         LSWORK = MAX( LWTRI, MAX( LW1, MAX( LW2, LW3 ) ) )
C     
C     Add together the síze of the whole DP WORK requirements
C     
         WRK = LSWORK + 3 * N
C     
         IF( .NOT. SCHRAE .AND. SOLVE )  WRK = WRK + 2*SIZEAE
C     
C     Compute the integer workspace
C     
         DBAE = ICEIL( N, DESCA(MB_) )
         IWRK = MAX(IWRK, DBAE + 2 * DESCA(MB_) + 8)
C
C     Now check if the call to this routine was a workspace query
C     and check if the workspace supplied is enough. 
C     
         IF( WRK.GT.LDWORK .AND. .NOT.LQUERY ) THEN
            INFO = -20
         ELSEIF ( LIWORK.LT.IWRK .AND. .NOT.LQUERY ) THEN
            INFO = -22
         ELSEIF( LQUERY ) THEN 
            DWORK( 1 ) = WRK
            IWORK( 1 ) = IWRK
            INFO = 0
            RETURN
         END IF
      END IF
C     
C     Make sure all global variables are indeed global and that INFO is
C     really set to the same value at all the processes
C     
      IF( INFO.EQ.0 ) THEN
         CALL PCHK1MAT( N, 5, N, 5, IA, JA, DESCA, 9, 0, IDUM1, IDUM2, 
     $                  INFO )
         IF( INFO.EQ.0 ) THEN
            CALL PCHK1MAT( N, 5, N, 5, IE, JE, DESCE, 13, 0, IDUM1, 
     $                     IDUM2, INFO )
         END IF
      END IF
      IF( INFO.EQ.0 .AND. SOLVE ) THEN
         CALL PCHK1MAT( N, 5, N, 5, IC, JC, DESCC, 17, 0, IDUM1,
     $                  IDUM2, INFO )
      END IF
C     
C     Checking if we may continue or should interrupt
C     
      IF( INFO.NE.0 ) THEN
         CALL PXERBLA( ICTXT, 'PGEGLYCTD', -INFO )
         RETURN
      END IF
C
C     Quick return if possible
C 
      IF( N.EQ.0 ) RETURN
C     
C     Partition workspace between local pointers
C
      IF( SOLVE ) THEN
         IF( SCHRAE ) THEN
            SWORK = 1
         ELSEIF( .NOT. SCHRAE ) THEN
            Q = 1
            Z = Q + SIZEAE
            CCOPY = Z + SIZEAE
            EIGR = CCOPY + SIZEE
            EIGI = EIGR + N
            BETA = EIGI + N
            TAUE = BETA + N
            SWORK = TAUE + NCOLSAE
         END IF
      ELSE
         IF( SCHRAE ) THEN
            SWORK = 1
         ELSEIF( .NOT. SCHRAE ) THEN
            EIGR = 1
            EIGI = EIGR + N
            BETA = EIGI + N
            TAUE = BETA + N
            SWORK = TAUE + NCOLSAE
         END IF
      END IF
      LSWORK = LDWORK-SWORK+1
C
C     Start the solution process - turn (A,E) into
C     Hessenberg-triangular form
C
C     QR-factorize E
C
      T1 = MPI_WTIME()
      IF( .NOT. SCHRAE ) 
     $     CALL PDGEQRF( N, N, E, IE, JE, DESCE, DWORK(TAUE), 
     $                   DWORK(SWORK), LSWORK, LINFO )
C
C     Update matrix A with respect to the QR-factorization
C
      IF( .NOT. SCHRAE )
     $     CALL PDORMQR( 'Left', 'Transpose', N, N, N, E, IE, JE, DESCE, 
     $                   DWORK(TAUE), A, IA, JA, DESCA, DWORK(SWORK),
     $                   LSWORK, LINFO )
C
C     Extract matrix Q from implicit storage in the matrix E.
C     Also initialize the matrix Z as identity of order N. Notice:
C     since we experienced problems with PDORGQR, we do this in
C     two different fashions depending on the logical parameter
C     GLYCT_PDORGQR. See also "Known issues" above.
C
      IF( GLYCT_PDORGQR ) THEN
         IF( .NOT. SCHRAE .AND. SOLVE ) THEN
            CALL PDLACPY( 'All', N, N, E, IE, JE, DESCE, DWORK(Q), IA, 
     $                    JA, DESCQZ )
            CALL PDORGQR( N, N, N, DWORK(Q), IA, JA, DESCQZ, 
     $                    DWORK(TAUE), DWORK(SWORK), LSWORK, LINFO )
            CALL PDLASET( 'All', N, N, ZERO, ONE, DWORK(Z), IA, JA, 
     $                    DESCQZ )
         END IF
      ELSE
         IF( .NOT. SCHRAE .AND. SOLVE ) THEN
            CALL PDLASET( 'All', N, N, ZERO, ONE, DWORK(Q), IA, JA, 
     $                    DESCQZ )
            CALL PDORMQR( 'Left', 'NoTranspose', N, N, N, E, IE, JE, 
     $                    DESCE, DWORK(TAUE), DWORK(Q), IA, JA, 
     $                    DESCQZ, DWORK(SWORK), LSWORK, LINFO ) 
            CALL PDLASET( 'All', N, N, ZERO, ONE, DWORK(Z), IA, JA, 
     $                    DESCQZ )
         END IF
      END IF
C    
C     Extract upper triangular matrix from E
C
      IF( .NOT. SCHRAE ) THEN
         CALL PDLASET( 'Lower triangular', N-1, N-1, ZERO, ZERO, E, 
     $                 IE+1, JE, DESCE )
      END IF
C
C     Call Hessenberg-triangular reduction routine for the 
C     full-triangular pair (A,E)
C      
      IF( .NOT. SCHRAE .AND. SOLVE ) THEN  
         CALL PDGGHRD( 'Vectors', 'Vectors', N, IA, IA+N-1, A, DESCA, E, 
     $                 DESCE, DWORK(Q), DESCQZ, DWORK(Z), DESCQZ, 
     $                 DWORK(SWORK), LSWORK, LINFO )
      ELSEIF( .NOT. SCHRAE ) THEN
         CALL PDGGHRD( 'NoVectors', 'NoVectors', N, IA, IA+N-1, A, 
     $                 DESCA, E, DESCE, DWORK(Q), DESCQZ, DWORK(Z), 
     $                 DESCQZ, DWORK(SWORK), LSWORK, LINFO )
      END IF
C 
C    Apply the QZ-algorithm
C
      IF( .NOT. SCHRAE .AND. SOLVE ) THEN
         BULGES = MAXBULG
         CALL PDHGEQZ( 'Schur', 'Vectors', 'Vectors', N, IA, IA+N-1, A, 
     $                 DESCA, E, DESCE, DWORK(EIGR), DWORK(EIGI), 
     $                 DWORK(BETA), DWORK(Q), DESCQZ, DWORK(Z), 
     $                 DESCQZ, IA, IA+N-1, IA, IA+N-1, BULGES, 
     $                 DWORK(SWORK), LSWORK, LINFO )
      ELSEIF( .NOT. SCHRAE ) THEN
         BULGES = MAXBULG
         CALL PDHGEQZ( 'Schur', 'NoVectors', 'NoVectors', N, IA, IA+N-1,
     $                 A, DESCA, E, DESCE, DWORK(EIGR), DWORK(EIGI), 
     $                 DWORK(BETA), DWORK(Q), DESCQZ, DWORK(Z), 
     $                 DESCQZ, IA, IA+N-1, IA, IA+N-1, BULGES, 
     $                 DWORK(SWORK), LSWORK, LINFO )
      END IF
      T1 = MPI_WTIME() - T1
      IF( .NOT. SOLVE ) RETURN
C     
C     Update C with respect to the Schur decomposition
C    
      T2 = MPI_WTIME()
      IF( .NOT. SCHRAE ) THEN
         IF( RSIDE ) THEN
            CALL PDLACPY( 'All', N, N, C, IC, JC, DESCC, DWORK( CCOPY ),
     $                    1 + IROFFC, 1 + ICOFFC, DESCCC )
            CALL PDGEMM( 'T', 'N', N, N, N, ONE, DWORK( Q ), IA, JA, 
     $                   DESCQZ, DWORK( CCOPY ), 1 + IROFFC, 1 + ICOFFC, 
     $                   DESCCC, ZERO, C, IC, JC, DESCC )
            CALL PDLACPY( 'All', N, N, C, IC, JC, DESCC, DWORK( CCOPY ),
     $                    1 + IROFFC, 1 + ICOFFC, DESCCC ) 
            CALL PDGEMM( 'N', 'N', N, N, N, ONE, DWORK( CCOPY ), 
     $                   1 + IROFFC, 1 + ICOFFC, DESCCC, DWORK( Q ), IA, 
     $                   JA, DESCQZ, ZERO, C, IC, JC, DESCC )
         ELSE
            CALL PDLACPY( 'All', N, N, C, IC, JC, DESCC, DWORK( CCOPY ),
     $                    1 + IROFFC, 1 + ICOFFC, DESCCC )
            CALL PDGEMM( 'T', 'N', N, N, N, ONE, DWORK( Z ), IA, JA, 
     $                   DESCQZ, DWORK( CCOPY ), 1 + IROFFC, 1 + ICOFFC, 
     $                   DESCCC, ZERO, C, IC, JC, DESCC )
            CALL PDLACPY( 'All', N, N, C, IC, JC, DESCC, DWORK( CCOPY ),
     $                    1 + IROFFC, 1 + ICOFFC, DESCCC ) 
            CALL PDGEMM( 'N', 'N', N, N, N, ONE, DWORK( CCOPY ), 
     $                   1 + IROFFC, 1 + ICOFFC, DESCCC, DWORK( Z ), 
     $                   IA, JA, DESCQZ, ZERO, C, IC, JC, DESCC )
         END IF
      END IF
      T2 = MPI_WTIME() - T2
C     
C     Now we have reduced our general equation to a (quasi-)triangular
C     case. Solve now the reduced problem with a call to a triangular 
C     solver - the solution to that reduced problem is stored
C     on C.
C
      T3 = MPI_WTIME()
      CALL PTRGLYCTD( CSYMM, OP, N, A, IA, JA, DESCA, E, IE, JE, DESCE, 
     $                C, IC, JC, DESCC, NB2, DWORK( SWORK ), LDWORK, 
     $                IWORK, LIWORK, NOEXSY, SCALE, INFO)
      IF( INFO.NE.0 ) THEN
         IF( LINFO.EQ.1 ) INFO = 2
         IF( LINFO.EQ.2 ) INFO = 3
      END IF
      T3 = MPI_WTIME() - T3
C     
C     Transform the solution back to the original coordinate system
C
      T22 = MPI_WTIME()
      IF( .NOT. SCHRAE ) THEN
         IF( RSIDE ) THEN
            CALL PDLACPY( 'All', N, N, C, IC, JC, DESCC, DWORK( CCOPY ),
     $                    1 + IROFFC, 1 + ICOFFC, DESCCC )
            CALL PDGEMM( 'N', 'N', N, N, N, ONE, DWORK( Z ), IA, JA, 
     $                   DESCQZ, DWORK( CCOPY ), 1 + IROFFC, 1 + ICOFFC, 
     $                   DESCCC, ZERO, C, IC, JC, DESCC )
            CALL PDLACPY( 'All', N, N, C, IC, JC, DESCC, DWORK( CCOPY ),
     $                    1 + IROFFC, 1 + ICOFFC, DESCCC ) 
            CALL PDGEMM( 'N', 'T', N, N, N, ONE, DWORK( CCOPY ), 
     $                   1 + IROFFC, 1 + ICOFFC, DESCCC, DWORK( Z ), IA, 
     $                   JA, DESCQZ, ZERO, C, IC, JC, DESCC )
         ELSE
            CALL PDLACPY( 'All', N, N, C, IC, JC, DESCC, DWORK( CCOPY ),
     $                    1 + IROFFC, 1 + ICOFFC, DESCC )
            CALL PDGEMM( 'N', 'N', N, N, N, ONE, DWORK( Q ), IA, JA,
     $                   DESCQZ, DWORK( CCOPY ), 1 + IROFFC, 1 + ICOFFC, 
     $                   DESCCC, ZERO, C, IC, JC, DESCC )
            CALL PDLACPY( 'All', N, N, C, IC, JC, DESCC, DWORK( CCOPY ),
     $                    1 + IROFFC, 1 + ICOFFC, DESCCC ) 
            CALL PDGEMM( 'N', 'T', N, N, N, ONE, DWORK( CCOPY ), 
     $                   1 + IROFFC, 1 + ICOFFC, DESCCC, DWORK( Q ), IA, 
     $                   JA, DESCQZ, ZERO, C, IC, JC, DESCC )
         END IF 
      END IF
      T2 = T2 + MPI_WTIME() - T22
C
      DWORK(1) = T1
      DWORK(2) = T2
      DWORK(3) = T3
C     
      END
C
C     END OF PGEGLYCTD 
C
C *** Last line of PGEGLYCTD ***
