CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C         This file is a part of the SCASY software library.           C
C         See http://www.cs.umu.se/parallel/research/scasy             C
C         Contributors: Robert Granat and Bo Kågström.                 C
C         SCASY Version 0.10,    March 31, 2006.                       C
C         Copyright 2006, Umeå University, Sweden.                     C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
       SUBROUTINE PGELYCTD( JOB, SYMM, OP, ASCHUR, N, A, IA, JA, DESCA, 
     $                      C, IC, JC, DESCC, NB2, DWORK, LDWORK, 
     $                      IWORK, LIWORK, NOEXSY, SCALE, INFO )
C
C  -- ScaLAPACK-style routine (preliminary version)--
C     HPC2N and Department of Computing Science,
C     University of Umeå, Sweden.
C     Written by Robert Granat, (granat@cs.umu.se)
C     January 25, 2006.
C
      IMPLICIT NONE
C
C     .. Scalar Arguments ..
      CHARACTER*1          JOB, SYMM, OP, ASCHUR
      INTEGER              N, IA, JA, IC, JC, LDWORK, LIWORK, NOEXSY, 
     $                     NB2, INFO
      DOUBLE PRECISION     SCALE
C     ..
C     .. Array Arguments ..
      INTEGER              DESCA( * ), DESCC( * ), IWORK( * )
      DOUBLE PRECISION     A( * ), C( * ), DWORK( * )
C 
C  Purpose
C  =======
C  The subroutine solves the general Lyapunov Equation 
C
C     op(sub(A)) * X + X * op(sub(A)^T) = C,
C
C  where sub(A), sub(C) and X (which overwrites sub(C) on output) are 
C  N-by-N matrices.
C
C  The notation op(_) denotes the transpose or non-transpose of a 
C  matrix. 
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
C  JOB       (global input) CHARACTER*1
C            If JOB = 'S', the equation is solved.
C            If JOB = 'R', only the reduction step regarding A is 
C            performed.
C
C  SYMM      (global input) CHARACTER*1
C            With JOB = 'S':
C            If SYMM = 'S', then the matrix C is assumed to be symmetric,
C            that is, C = C^T. If SYMM = 'N', the matrix C is not assumed
C            to be symmetric.
C            Otherwise, SYMM is not referenced.
C     
C  OP        (global input) CHARACTER*1
C            With JOB = 'S':
C              If OP = 'N', then we solve A * X + X * A^T = C
C              If OP = 'T', then we solve A^T * X + X * A = C.
C            Otherwise, OP is not referenced.
C
C  ASCHUR    (global input) CHARACTER*1
C            If ASCHUR = 'S', then the matrix A is supposed to be an
C            upper (quasi-)triangular matrix. No reduction to the real
C            Schur form of this matrix will be done.
C            If ASCHUR = 'N', then the matrix A is not in upper (quasi-)
C            triangular form and a reduction to the real Schur form will
C            be done.
C
C  Input/Output parameters
C
C  N         (global input) INTEGER
C            The number of rows and columns of the global distributed
C            matrix A. This is also the number of rows and columns of the
C            global distributed matrix C (and X). N >= 0.
C
C  A         (local input/output) DOUBLE PRECISION array 
C            Array of dimension (LLD_A,LOCc(M)). Contains the local
C            pieces of the global distributed matrix A. On output,
C            the local part of the distributed matrix A in real
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
C  C         (local input/local output) DOUBLE PRECISION array 
C            Array of dimension (LLD_C,LOCc(N)). 
C            With JOB = 'S': On entry, C contains the local pieces of 
C            the global distributed matrix C . On exit, it contains the
C            local pieces of the global distributed solution X.
C            If SYMM = 'S':

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
C  NB2       (global input) INTEGER 
C            Internal blocking factor for pipelining of subsolutions
C            for updates of the matrix C in PTRLYCTD (see the references 
C            for details).
C            1 < = NB2 <= DESCC( NB_ ) must hold.
C
C  DESCC     (global and local input) INTEGER array of dimension DLEN_.
C            With JOB = 'S', the array descriptor for the global 
C            distributed matrix C.
C            Otherwise, DESCC is not referenced.
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
C            LIWORK >= DBA, where DBA = ICEIL(LOCr(IA+IROFFA),MB_A).
C
C            If LIWORK = -1, LIWORK is global input and a workspace 
C            query is assumed. The routine will then calculate the 
C            optimal workspace needed, store it in IWORK(1) and return
C            immediately. No error will then be signaled by PXERBLA.
C            
C  Output information
C
C  NOEXSY    (local output) INTEGER
C            When solving the triangular problem in PTRLYCTD it is possible
C            that we have to extend some subsystems to not lose any data
C            from some 2x2 block of conjugate pairs of eigenvalues. 
C            NOEXSY helps us to keep track of the number of such extensions. 
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
C  This subroutine implements a parallel version of a Bartels-Stewart-like
C  method for solving the general continous-time Lyapunov equation (see
C  the references for details).
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
C  For ASCHUR = 'N' it is not possible to work with submatrices, i.e., 
C  ASCHUR = 'N' implies IA = JA = 1. This limitation comes from the fact 
C  that ScaLAPACKs PDLAHQR does not work on submatrices - this might be 
C  changed in a future release.
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
C  Symmetric right hand side update fixed, 17 December, 2007.
C
C  Report bugs to <granat@cs.umu.se>
C
C  Keywords
C  ========
C  QR-algorithm, real Schur form, continous-time Lyapunov equation, 
C  Hessenberg transformation, PDGEMM-updates, symmetric matrix.
C
C  =====================================================================
C
C     ..Parameters ..
      INTEGER            BLOCK_CYCLIC_2D, CSRC_, CTXT_, DLEN_, DTYPE_,
     $                   LLD_, MB_, M_, NB_, N_, RSRC_, NOINARG,
     $                   WORKSPACE_Q
      DOUBLE PRECISION   ZERO, ONE, MINONE, HALF
      CHARACTER*1        COPY_ALL
      PARAMETER          ( BLOCK_CYCLIC_2D = 1, DLEN_ = 9, DTYPE_ = 1,
     $                     CTXT_ = 2, M_ = 3, N_ = 4, MB_ = 5, NB_ = 6,
     $                     RSRC_ = 7, CSRC_ = 8, LLD_ = 9, NOINARG = 14,
     $                     WORKSPACE_Q = -1, ZERO = 0.0D+0,ONE = 1.0D+0,
     $                     MINONE = -1.0D+0, HALF = 0.5D+00, 
     $                     COPY_ALL = 'A')         
C     .. 
C     .. Local Scalars ..
      LOGICAL            SOLVE, LQUERY, A_SCHUR, CSYMM, SIDE, WANTT, 
     $                   WANTZ
      DOUBLE PRECISION   T1, T2, T3, T22
      INTEGER            ICTXT, MYCOL, MYROW, NPCOL, NPROW, D, Da,
     $                   NROWS_A, NCOLS_A, LSWORK, SIZE_A,
     $                   SIZE_C, SIZE_TAUA, SIZE_WI, SIZE_WR, 
     $                   LWORK1, LWORK2, LWORK3, LWORK4, WRK, IWRK, 
     $                   DUMMY, NPROCS, I, IS, IX, JX, RSRC, CSRC, 
     $                   DBA, IROFFA, ICOFFA, IROFFC, ICOFFC, LIA, LJA,
     $                   ARSRC, ACSRC, LIC, LJC, CRSRC, CCSRC, NROWS_C,
     $                   K, W
C     ..
C     .. Local Arrays and Local Pointers ..
      INTEGER             Z_A, C_COPY, TAUA, WR, WI, SWORK, 
     $                    IDUM1( 1 ), IDUM2( 1 ), DUMMY_ARRAY( 1 ),
     $                    DESCZ_A( DLEN_ ), DESCCC( DLEN_ ),
     $                    DESCW( DLEN_ )
      DOUBLE PRECISION    DPDUM1( 1 ), DPDUM2( 1 ), DPDUM3( 1 ), 
     $                    DPDUM4( 1 )
C     ..
C     .. External Subroutines ..
C     
      EXTERNAL           BLACS_GRIDINFO, CHK1MAT, PCHK1MAT, PXERBLA,
     $                   PDGEHRD, PDORMHR, PDLAHQR, PDGEMM, PDLACPY, 
     $                   PDIDMAT, PDHESS, PTRLYCTD, DLATCPY, PDCOPY
C     .. 
C     .. External Functions ..
      DOUBLE PRECISION   MPI_WTIME
      INTEGER            NUMROC, ILCM, ICEIL
      LOGICAL            LSAME
      EXTERNAL           NUMROC, ILCM, ICEIL, LSAME, MPI_WTIME
      
C     .. Intrinsic Functions ..
      INTRINSIC MAX, MIN, MOD, INT, REAL
C     ..
C     .. Executable Statements ..
C     
C     Get grid parameters
C     
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
      IF(.NOT.( LSAME( JOB, 'R' ) .OR. LSAME( JOB, 'S' ))) 
     $     THEN
         INFO = -1
      ELSE
         IF( LSAME( JOB, 'S' ) ) THEN
            SOLVE = .TRUE.
            WANTT = .TRUE.
            WANTZ = .TRUE.
         ELSE
            SOLVE = .FALSE.
            WANTT = .TRUE.
            WANTZ = .FALSE.
         END IF
      END IF
C     
C     Check dimensions
C     
      IF( INFO.EQ.0 ) THEN
         CALL CHK1MAT( N, 5, N, 5, IA, JA, DESCA, 9, INFO )
         IF( INFO.EQ.0 .AND. SOLVE ) THEN
            CALL CHK1MAT( N, 5, N, 5, IC, JC, DESCC, 13, INFO )
         END IF
      END IF
C     
C     Check the blocking sizes for equivalence and being not to small
C     
      IF( INFO.EQ.0 ) THEN
         IF( DESCA( MB_ ).NE.DESCA( NB_ ) ) INFO = -(100*9 + MB_)
      END IF
      IF( INFO.EQ.0 ) THEN
         IF( DESCC( MB_ ).NE.DESCA( MB_ ) ) INFO = -(100*13 + MB_)
      END IF
      IF( INFO.EQ.0 .AND. SOLVE ) THEN
         IF( DESCC( NB_ ).NE.DESCA( NB_ ) ) INFO = -(100*13 + NB_)
      END IF
C
      IF( INFO.EQ.0 ) THEN
         IF( N.NE.DESCA( MB_ ) .AND. DESCA( MB_ ).LT.6 .AND. 
     $        LSAME( ASCHUR, 'N' ) ) 
     $        INFO = -(100*9 + MB_)
      END IF
C
C     Check the value of SYMM and set the logical value of
C     CSYMM for future reference
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
C     Check the value of OP
C 
      IF( INFO.EQ.0 .AND. SOLVE ) THEN
         IF( .NOT.( LSAME( OP, 'T' ) .OR. LSAME( OP, 'N' ))) 
     $        THEN
            INFO = -3
         END IF
      END IF        
C
C     Check the value of ASCHUR and set the logical value of
C     A_SCHUR for future reference
C 
      IF( INFO.EQ.0 ) THEN
         IF( .NOT.( LSAME( ASCHUR, 'N' ) .OR. LSAME( ASCHUR, 'S' ))) 
     $        THEN
            INFO = -4
         ELSE
            IF( LSAME( ASCHUR, 'S' ) ) THEN
               A_SCHUR = .TRUE.
            ELSE
               A_SCHUR = .FALSE.
            END IF
         END IF
      END IF
C
C     Check the values of IA, JA, IC, JC
C
      IF( INFO.EQ.0 ) THEN
         IF( IA.LT.1 .OR. ( .NOT. A_SCHUR .AND. IA.GT.1 ) ) INFO = -7
         IF( JA.NE.IA ) INFO = -8
      END IF
      IF( INFO.EQ.0 .AND. SOLVE ) THEN
         IF( MOD(IC,DESCC(MB_)) .NE. MOD(IA,DESCA(MB_)) ) INFO = -11
         IF( MOD(JC,DESCC(NB_)) .NE. MOD(JA,DESCA(NB_)) ) INFO = -12
      END IF
C
C     Check if workspace query
C
      LQUERY = LDWORK.EQ.-1 .OR. LIWORK.EQ.-1
C     
C     Test working space
C     
C     Check the number of rows and columns of the smallest submatrix
C     of A including sub(A) that conform with the ScaLAPACK conventions
C     
      IF( INFO.EQ.0 .OR. LQUERY ) THEN
         IROFFA = MOD( IA - 1, DESCA(MB_) )
         ICOFFA = MOD( JA - 1, DESCA(NB_) )
         IROFFC = MOD( IC - 1, DESCC(MB_) )
         ICOFFC = MOD( JC - 1, DESCC(NB_) )
         CALL INFOG2L( IA, JA, DESCA, NPROW, NPCOL, MYROW, MYCOL,
     $                 LIA, LJA, ARSRC, ACSRC )
         CALL INFOG2L( IC, JC, DESCC, NPROW, NPCOL, MYROW, MYCOL,
     $                 LIC, LJC, CRSRC, CCSRC )
         NROWS_A = NUMROC( N, DESCA(MB_), MYROW, DESCA(RSRC_), NPROW )
         NROWS_C = NUMROC( N, DESCC(MB_), MYROW, DESCC(RSRC_), NPROW )
         NCOLS_A = NUMROC( N, DESCA(NB_), MYCOL, DESCA(CSRC_), NPCOL )
         NCOLS_A = NUMROC( N, DESCC(NB_), MYCOL, DESCC(CSRC_), NPCOL )
C
C     Init matrix descriptor for Z_A and W or C_COPY
C
         IF( .NOT. A_SCHUR .AND. SOLVE ) THEN
            CALL DESCINIT( DESCZ_A, N + IROFFA, N + ICOFFA, DESCA(MB_), 
     $                     DESCA(NB_), ARSRC, ACSRC, ICTXT, 
     $                     MAX( 1, NROWS_A ), INFO )
            IF( .NOT. CSYMM ) THEN
               CALL DESCINIT( DESCCC, N + IROFFC, N + ICOFFC, 
     $                        DESCC(MB_), DESCC(NB_), CRSRC, CCSRC, 
     $                        ICTXT, MAX( 1, NROWS_C ), INFO )
            ELSE
               CALL DESCINIT( DESCW, N + IROFFA, N + ICOFFA, 
     $                        DESCA(MB_), DESCA(NB_), ARSRC, ACSRC, 
     $                        ICTXT, MAX( 1, NROWS_A ), INFO ) 
           END IF
         END IF
C     
C     Compute the sizes of the different needed workspaces
C     
         SIZE_A = DESCA( LLD_ ) * NCOLS_A 
         SIZE_C = SIZE_A
         SIZE_TAUA = NUMROC( N-1, DESCA( NB_ ), MYCOL, DESCA( CSRC_ ), 
     $                       NPCOL )
         SIZE_WR = N
         SIZE_WI = N
C     
C     Do some workspace queries to the routines called below
C     
         IF( .NOT. A_SCHUR ) THEN
            CALL PDGEHRD( N, 1, N, A, IA, JA, DESCA, DPDUM1, 
     $                    DPDUM2, WORKSPACE_Q, INFO )
            LWORK1 = INT( DPDUM2( 1 ) )
         ELSE
            LWORK1 = 0
         END IF
C     
         IF( .NOT. A_SCHUR .AND. SOLVE ) THEN
            CALL PDORMHR( 'L','N', N, N, 1, N, A, IA, JA, DESCA, 
     $                    DPDUM1, DWORK, 1+IROFFA, 1+ICOFFA, DESCZ_A, 
     $                    DPDUM2, WORKSPACE_Q, INFO )
            LWORK2 = INT( DPDUM2( 1 ) )
         ELSE
            LWORK2 = 0
         END IF
C
         IF( .NOT. A_SCHUR ) THEN
            CALL PDLAHQR( WANTT, WANTZ, N, IA, IA+N-1, A, DESCA, 
     $                    DWORK, DWORK, IROFFA+1, IROFFA+N, DWORK, 
     $                    DESCZ_A, DPDUM1, WORKSPACE_Q, IDUM1, 
     $                    WORKSPACE_Q, INFO )
            LWORK3 = INT( DPDUM1( 1 ) )
            IWRK = IDUM1(1)
         ELSE
            LWORK3 = 0
            IWRK = 0
         END IF
C
         IF( SOLVE ) THEN
            CALL PTRLYCTD( CSYMM, OP, N, A, IA, JA, DESCA, C, IC, JC, 
     $                     DESCC, NB2, DPDUM1, WORKSPACE_Q, IDUM1, 
     $                     LIWORK, NOEXSY, SCALE, INFO )
            LWORK4 = INT( DPDUM1( 1 ) )
            IWRK = MAX( IWRK, IDUM1( 1 ))
         ELSE
            LWORK4 = 0
            IWRK = MAX( IWRK, 0 )
         END IF
C     
         LSWORK = MAX( LWORK4, MAX( LWORK1,
     $            MAX( LWORK2, LWORK3 ) ) )
C     
C     Add together the síze of the whole DP WORK requirements
C     
         WRK = LSWORK
C     
         IF( .NOT. A_SCHUR )  WRK = WRK + SIZE_TAUA + SIZE_WR +
     $                              SIZE_WI
         IF( SOLVE ) THEN
            IF( .NOT. A_SCHUR ) WRK = WRK + SIZE_A + SIZE_C
         END IF
C
C     Now check if the call to this routine was a workspace query
C     and check if the workspace supplied is enough. 
C     
         IF( WRK.GT.LDWORK .AND. .NOT.LQUERY ) THEN
            INFO = -16
         ELSEIF( IWRK.GT.LIWORK .AND. .NOT. LQUERY ) THEN
            INFO = -18
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
         CALL PCHK1MAT( N, 5, N, 5, 1, 1, DESCA, 9, 0, IDUM1, 
     $                  IDUM2, INFO )
      END IF
      IF( INFO.EQ.0 .AND. SOLVE ) THEN
         CALL PCHK1MAT( N, 5, N, 5, 1, 1, DESCC, 13, 0, IDUM1,
     $                  IDUM2, INFO )
      END IF
C     
C     Checking if we may continue or should interrupt
C     
      IF( INFO.NE.0 ) THEN
         CALL PXERBLA( ICTXT, 'PGELYCTD', -INFO )
         RETURN
      END IF
C
C     Quick return if possible
C 
      IF( N.EQ.0 ) RETURN
C     
C     Divide the workspace between the different local pointers
C     which are needed
C 
      IF( SOLVE ) THEN
         IF( A_SCHUR ) THEN
            SWORK = 1
         ELSEIF( .NOT. A_SCHUR ) THEN
            Z_A = 1
            IF( .NOT. CSYMM ) THEN
               C_COPY = Z_A + SIZE_A
               TAUA = C_COPY + SIZE_C
            ELSE
               W= Z_A + SIZE_A
               TAUA = W + SIZE_A
            END IF
            WR = TAUA + SIZE_TAUA
            WI = WR + SIZE_WR
            SWORK = WI + SIZE_WI
         END IF 
      ELSE
         IF( A_SCHUR ) THEN
            SWORK = 1
         ELSEIF( .NOT. A_SCHUR ) THEN
            TAUA = 1
            WR = TAUA + SIZE_TAUA
            WI = WR + SIZE_WR
            SWORK = WI + SIZE_WI
         END IF 
      END IF
      LSWORK = LDWORK-SWORK+1
C
C     Turn A into upper Hessenberg form 
C 
      T1 = MPI_WTIME()
      IF( .NOT. A_SCHUR ) THEN
         CALL PDGEHRD( N, 1, N, A, IA, JA, DESCA, DWORK( TAUA ),
     $                 DWORK( SWORK ), LSWORK, INFO )
C
C     Form Q_A explicitely in Z_A by using PDORMHR on an identity
C     matrix stored in Z_A. 
C
         IF( SOLVE ) THEN
            CALL PDLASET( 'A', N, N, ZERO, ONE, DWORK(Z_A), 1+IROFFA, 
     $                    1+ICOFFA, DESCZ_A )   
C     
            CALL PDORMHR( 'L','N', N, N, 1, N, A, IA, JA, DESCA, 
     $                    DWORK( TAUA ), DWORK( Z_A ), 1+IROFFA, 
     $                    1+ICOFFA, DESCZ_A, DWORK( SWORK ), LSWORK, 
     $                    INFO )
         END IF
C
C     Extract the upper Hessenberg part of A 
C   
         CALL PDLASET( 'Lower triangular', N-2, N-2, ZERO, ZERO, A, 
     $                 IA + 2, JA, DESCA )
C     
C     Compute the real Schur form of the Hessenberg matrix A
C
         CALL PDLAHQR( WANTT, WANTZ, N, IA, IA+N-1, A, DESCA, 
     $                 DWORK( WR ), DWORK( WI ), 1+IROFFA, IROFFA+N, 
     $                 DWORK( Z_A ), DESCZ_A, DWORK( SWORK ), LSWORK, 
     $                 IWORK, LIWORK, INFO )
C
      END IF
      T1 = MPI_WTIME() - T1
      IF( .NOT. SOLVE ) RETURN
C
C     Update C with respect to the transformations done
C     Bug fix 17.12.07: with U = triu(C) - 0.5*diag(diag(C)),
C     we compute CC = Q'*C*Q as CC = Q'*U*Q + (Q'*U*Q)' 
C     (symmetric right hand side only).
C
      T2 = MPI_WTIME()
      IF( .NOT. A_SCHUR ) THEN
         IF( .NOT. CSYMM ) THEN
            CALL PDLACPY( COPY_ALL, N, N, C, IC, JC, DESCC, 
     $           DWORK( C_COPY ), IROFFC+1, ICOFFC+1, DESCCC )
            CALL PDGEMM( 'Transpose', 'Notranspose', N, N, N, ONE, 
     $           DWORK( Z_A ), IROFFA+1, ICOFFA+1, DESCZ_A, 
     $           DWORK( C_COPY ), IROFFC+1, ICOFFC+1, DESCCC, 
     $           ZERO, C, IC, JC, DESCC )
            CALL PDLACPY( COPY_ALL, N, N, C, IC, JC, DESCC, 
     $           DWORK( C_COPY ), IROFFC+1, ICOFFC+1, DESCCC )
            CALL PDGEMM( 'Notranspose','Notranspose', N, N, N, ONE, 
     $           DWORK( C_COPY ), IROFFC+1, ICOFFC+1, DESCCC, 
     $           DWORK( Z_A ), IROFFA+1, ICOFFA+1, DESCA, 
     $           ZERO, C, IC, JC, DESCC )
         ELSE
            DO 10 K = 1, N
               CALL PDSCAL( 1, HALF, C, IC+K-1, JC+K-1, DESCC, 1 )
 10         CONTINUE
            CALL PDLACPY( COPY_ALL, N, N, DWORK(Z_A), IROFFA+1, 
     $                    ICOFFA+1, DESCZ_A,  DWORK( W ), IROFFA+1, 
     $                    ICOFFA+1, DESCW )
            CALL PDTRMM( 'Left', 'Lower', 'Notranspose', 'Non-Unit', N, 
     $                   N, ONE, C, IC, JC, DESCC, DWORK( W ), IROFFA+1, 
     $                   ICOFFA+1, DESCW )
            CALL PDSYR2K( 'Lower', 'Transpose', N, N, ONE, DWORK( W ), 
     $                    IROFFA+1, ICOFFA+1, DESCW, DWORK( Z_A ), 
     $                    IROFFA+1, ICOFFA+1, DESCZ_A, ZERO, C, IC, JC, 
     $                    DESCC )
            DO 15 K = 1, N - 1
               CALL PDCOPY( N-K, C, IC+K, JC+K-1, DESCC, 1, C, IC+K-1, 
     $                      JC+K, DESCC, DESCC(M_) )
 15         CONTINUE
         END IF
      END IF
      T2 = MPI_WTIME() - T2
C     
C     Solve the reduced triangular problem
C
      T3 = MPI_WTIME()
      CALL PTRLYCTD( CSYMM, OP, N, A, IA, JA, DESCA, C, IC, JC, DESCC, 
     $               NB2, DWORK( SWORK ), LSWORK, IWORK, LIWORK, NOEXSY, 
     $               SCALE, INFO )
C
      IF( INFO.NE.0 ) THEN
         IF( INFO.EQ.1 ) INFO = 2
         IF( INFO.EQ.2 ) INFO = 3
      END IF
      T3 = MPI_WTIME() - T3
C     
C     Transform the solution back to the original coordinate system
C
      T22 = MPI_WTIME()
      IF( .NOT. A_SCHUR ) THEN
         IF( .NOT. CSYMM ) THEN
            CALL PDLACPY( COPY_ALL, N, N, C, IC, JC, DESCC, 
     $           DWORK( C_COPY ), IROFFC+1, ICOFFC+1, DESCCC )
            CALL PDGEMM( 'NoTranspose', 'Notranspose', N, N, N, ONE, 
     $           DWORK( Z_A ), IROFFA+1, ICOFFA+1, DESCZ_A, 
     $           DWORK( C_COPY ), IROFFC+1, ICOFFC+1, DESCCC, 
     $           ZERO, C, IC, JC, DESCC )
            CALL PDLACPY( COPY_ALL, N, N, C, IC, JC, DESCC, 
     $           DWORK( C_COPY ), IROFFC+1, ICOFFC+1, DESCCC )
            CALL PDGEMM( 'Notranspose','Transpose', N, N, N, ONE, 
     $           DWORK( C_COPY ), IROFFC+1, ICOFFC+1, DESCCC, 
     $           DWORK( Z_A ), IROFFA+1, ICOFFA+1, DESCA, 
     $           ZERO, C, IC, JC, DESCC )
         ELSE
            DO 20 K = 1, N
               CALL PDSCAL( 1, HALF, C, IC+K-1, JC+K-1, DESCC, 1 )
 20         CONTINUE
            CALL PDLACPY( COPY_ALL, N, N, DWORK(Z_A), IROFFA+1, 
     $                    ICOFFA+1, DESCZ_A,  DWORK( W ), IROFFA+1, 
     $                    ICOFFA+1, DESCW )
            CALL PDTRMM( 'Right', 'Lower', 'Notranspose', 'Non-Unit', N, 
     $                   N, ONE, C, IC, JC, DESCC, DWORK( W ), IROFFA+1, 
     $                   ICOFFA+1, DESCW )
            CALL PDSYR2K( 'Lower', 'Notranspose', N, N, ONE, DWORK( W ), 
     $                    IROFFA+1, ICOFFA+1, DESCW, DWORK( Z_A ), 
     $                    IROFFA+1, ICOFFA+1, DESCZ_A, ZERO, C, IC, JC, 
     $                    DESCC )
            DO 25 K = 1, N - 1
               CALL PDCOPY( N-K, C, IC+K, JC+K-1, DESCC, 1, C, IC+K-1, 
     $                      JC+K, DESCC, DESCC(M_) )
 25         CONTINUE
         END IF
      END IF
      T2 = T2 + MPI_WTIME() - T22 
C
      DWORK( 1 ) = T1
      DWORK( 2 ) = T2
      DWORK( 3 ) = T3
C
      END
C
C     END OF PGELYCTD 
C
C *** Last line of PGELYCTD ***
