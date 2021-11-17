CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C         This file is a part of the SCASY software library.           C
C         See http://www.cs.umu.se/parallel/research/scasy             C
C         Contributors: Robert Granat and Bo Kågström.                 C
C         SCASY Version 0.10,    March 31, 2006.                       C
C         Copyright 2006, Umeå University, Sweden.                     C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE PGESYDTD( JOB, ASCHUR, BSCHUR, TRANSA, TRANSB, ISGN,
     $                     COMM, M, N, A, IA, JA, DESCA, B, IB, JB, 
     $                     DESCB, C, IC, JC, DESCC, MB2, DWORK, 
     $                     LDWORK, IWORK, LIWORK, NOEXSY, SCALE, INFO )
C
C  -- ScaLAPACK-style routine (preliminary version)--
C     HPC2N and Department of Computing Science,
C     University of Umeå, Sweden.
C     Written by Robert Granat, (granat@cs.umu.se)
C     January 27, 2006.
C
      IMPLICIT NONE
C
C     .. Scalar Arguments ..
      CHARACTER*1          JOB, ASCHUR, BSCHUR, TRANSA, TRANSB, COMM
      INTEGER              ISGN, M, N, IA, JA, IB, JB, IC, JC, LDWORK, 
     $                     LIWORK, NOEXSY, INFO, MB2
      DOUBLE PRECISION     SCALE
C     ..
C     .. Array Arguments ..
      INTEGER              DESCA( * ), DESCB( * ), DESCC( * ), 
     $                     IWORK( * )
      DOUBLE PRECISION     A( * ), B( * ), C( * ), DWORK( * )
C 
C  Purpose
C  =======
C  The subroutine solves the general Discrete Time Sylvester 
C  Equation (SYDT) 
C
C     op(sub(A)) * X * op(sub(B)) +/- X  = sub(C),
C
C  where sub(A) is an M-by-M matrix, sub(B) is an N-by-N matrix and the 
C  solution X is an M-by-N matrix which overwrites sub(C). 
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
C            If JOB = 'S', the equation is solved. If JOB = 'R',
C            only the reduction step regarding A and B is performed.
C
C  ASCHUR    (global input) CHARACTER*1
C            If ASCHUR = 'S', then the matrix A is supposed to be an
C            upper (quasi-)triangular matrix. No reduction to the real
C            Schur form of this matrix will be done.
C            If ASCHUR = 'N', then the matrix A is not in upper (quasi-)
C            triangular form and a reduction to the real Schur form will
C            be done.
C
C  BSCHUR    (global input) CHARACTER*1
C            If BSCHUR = 'S', then the matrix B is supposed to be an
C            upper (quasi-)triangular matrix. No reduction to the real
C            Schur form of this matrix will be done.
C            If BSCHUR = 'N', then the matrix B is not in upper (quasi-)
C            triangular form and a reduction to the real Schur form will
C            be done.
C
C  TRANSA    (global input) CHARACTER*1
C            With JOB = 'S':
C              If TRANSA = 'N' then op(A) = A
C              If TRANSA = 'T' then op(A) = A**T
C            Otherwise, TRANSA is not referenced.
C
C  TRANSB    (global input) CHARACTER*1
C            With JOB = 'S':
C              If TRANSB = 'N' then op(B) = B
C              If TRANSB = 'T' then op(B) = B**T
C            Otherwise, TRANSB is not referenced.
C
C  ISGN      (global input) INTEGER*1
C            With JOB = 'S':
C              If ISGN =  1, we solve the equation 
C                op(A) * X + X * op(B) = C
C              If ISGN = -1, we solve the equation 
C                op(A) * X - X * op(B) = C
C            Otherwise, ISGN is not referenced.
C
C  Input/Output parameters
C
C  COMM      (global input/output) CHARACTER*1
C            This subroutine uses two different communications schemes in
C            solving the reduced triangular problem. 
C            With JOB = 'S':
C              If COMM = 'S', the "shifts" scheme is used.
C              If COMM = 'D', the "communicate on demand" scheme is used.
C            The choice COMM = 'S' is only valid for TRANSA = TRANSB = 'N' 
C            or TRANSA = TRANSB = 'T'. The scheme used will be output.
C            See the references for details.
C            Otherwise, COMM is not referenced.
C
C  M         (global input) INTEGER
C            The number of rows and columns of the global distributed 
C            matrix sub(A). This is also the number of rows of the
C            global distributed matrix sub(C). M >= 0.
C
C  N         (global input) INTEGER
C            The number of rows and columns of the global distributed
C            matrix sub(B). This is also the number of columns of the
C            global distributed matrix sub(C). N >= 0.
C
C  A         (local input/output) DOUBLE PRECISION array 
C            Array of dimension (LLD_A,LOCc(M)). Contains the local
C            pieces of the global distributed matrix A. On output,
C            the local parts of the distributed matrix A in real
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
C  B         (local input/output) DOUBLE PRECISION array 
C            Array of dimension (LLD_B,LOCc(N)). Contains the local
C            pieces of the global distributed matrix B. On output,
C            the local parts of the distributed matrix B in real
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
C            With JOB='S': On entry C contains the local pieces of the
C            global distributed matrix C . On exit, it contains the local
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
C            submatrix to operate on. MOD(JC,MB_B) = MOD(JB,MB_B) 
C            must hold.
C            Otherwise, JC is not referenced. 
C
C  DESCC     (global and local input) INTEGER array of dimension DLEN_.
C            With JOB = 'S', the array descriptor for the global 
C            distributed matrix C.
C            Otherwise, DESCC is not referenced.
C
C  MB2       (global input) INTEGER
C            Internal blocking factors for pipelining of subsolutions
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
C            LIWORK >= DBA + DBB, where DBA = ICEIL(LOCr(IA+IROFFA),MB_A)
C            and DBB = ICEIL(LOCr(IB+IROFFB),MB_B).
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
C            When solving the triangular problem in PTRSYCTD it is possible
C            that we have to extend some subsystems to not lose any data
C            from some 2x2 block of conjugate pairs of eigenvalues. 
C            NOEXSY helps us to keep track of the number of such 
C            extensions. Otherwise, NOEXSY is not referenced.
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
C             If INFO = 2, some eigenvalue of A is equal to or nearly
C             equal to the reciprocal of an eigenvalue of -ISGN*B.
C             Perturbed values was used to solve the equation (but A
C             and B are unchanged).
C             If INFO = 3,  the problem is badly scaled - C should have 
C             been scaled a factor SCALE before calling this routine too 
C             guarantee an overflow free solution, i.e., the solution
C             may well have overflowed.
C
C  Method
C  ======
C  This subroutine implements a parallel version of the Bartels-Stewart
C  method for solving the general discrete-time Sylvester equation (see
C  the references for details).
C
C  Additional requirements
C  =======================
C  A and B must be distributed using the same blocking factor in 
C  each direction, i.e., MB_A=NB_A, MB_B=NB_B. Moreover, for C the 
C  blocksize in the row direction must agree with A's, i.e. MB_C=MB_A
C  must hold, and the blocksize in the column direction must agree with
C  B's, i.e. NB_C=NB_B must hold. The blocksizes must be larger than
C  or equal to six (6).
C
C  Limitations
C  ===========
C  In contrary to SLICOTs SB04PY this routine do not scale against 
C  overflow in the solution. See SCALE and INFO.
C
C  For ASCHUR = 'N' or BSCHUR = 'N' it is not possible to work with
C  submatrices, e.g., ASCHUR = 'N' implies IA = JA = 1. This limitation
C  comes from the fact that ScaLAPACKs PDLAHQR does not work on 
C  submatrices - this might be changed in a future release.
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
C  QR-algorithm, real Schur form, Sylvester equation, 
C  Hessenberg transformation, PDGEMM-updates.
C
C  =====================================================================
C
C     ..Parameters ..
      INTEGER            BLOCK_CYCLIC_2D, CSRC_, CTXT_, DLEN_, DTYPE_,
     $                   LLD_, MB_, M_, NB_, N_, RSRC_, NOINARG,
     $                   WRKSPCQ
      DOUBLE PRECISION   ZERO, ONE, MINONE
      CHARACTER*1        COPY_ALL
      PARAMETER          ( BLOCK_CYCLIC_2D = 1, DLEN_ = 9, DTYPE_ = 1,
     $                     CTXT_ = 2, M_ = 3, N_ = 4, MB_ = 5, NB_ = 6,
     $                     RSRC_ = 7, CSRC_ = 8, LLD_ = 9, NOINARG = 20,
     $                     WRKSPCQ = -1, ZERO = 0.0D0, ONE = 1.0D0, 
     $                     MINONE = -1.0D0, COPY_ALL = 'A')         
C     .. 
C     .. Local Scalars ..
      LOGICAL            SOLVE, LQUERY, A_SCHUR, B_SCHUR, WANTT, WANTZ
      DOUBLE PRECISION   FNORM, ANORM, BNORM, XNORM, T1, T2, T3, T22
      INTEGER            ICTXT, MYCOL, MYROW, NPCOL, NPROW, D, Da, Db,
     $                   NROWS_A, NROWS_B, NCOLS_A, NCOLS_B, NCOLS_C,
     $                   LSWORK, SIZE_A, SIZE_B, SIZE_C, SIZE_TAUA, 
     $                   SIZE_TAUB, SIZE_WI, SIZE_WR, LWORK1, 
     $                   LWORK2, LWORK3, LWORK4, WRK, MAXMN, DUMMY,
     $                   NPROCS, I, IS, IX, JX, RSRC, CSRC, DBA, DBB,
     $                   IROFFA, ICOFFA, IROFFB, ICOFFB, IROFFC, 
     $                   ICOFFC, LIA, LJA, ARSRC, ACSRC, LIB, LJB, 
     $                   BRSRC, BCSRC, LIC, LJC, CRSRC, CCSRC,
     $                   NROWS_C, IWRK
C     ..
C     .. Local Arrays and Local Pointers ..
      INTEGER             Z_A, Z_B, C_COPY, TAUA, TAUB, WR, WI, SWORK, 
     $                    IDUM1( 1 ), IDUM2( 1 ), DUMMY_ARRAY( 1 ),
     $                    DESCZ_A( DLEN_ ), DESCZ_B( DLEN_ ),
     $                    DESCCC( DLEN_ )
      DOUBLE PRECISION    DPDUM1( 1 ), DPDUM2( 1 ), DPDUM3( 1 ), 
     $                    DPDUM4( 1 )
C     ..
C     .. External Subroutines ..
C     
      EXTERNAL           BLACS_GRIDINFO, CHK1MAT, PCHK1MAT, PXERBLA,
     $                   PDGEHRD, PDORMHR, PDLAHQR, PDGEMM, PDLACPY, 
     $                   PDHESS, PTRSYDTD
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
C     Check dimension
C     
      IF( INFO.EQ.0 ) THEN
         CALL CHK1MAT( M, 8, M, 8, IA, JA, DESCA, 13, INFO )
         IF( INFO.EQ.0 ) THEN
            CALL CHK1MAT( N, 9, N, 9, IB, JB, DESCB, 17, INFO )
         END IF
         IF( INFO.EQ.0 .AND. SOLVE ) THEN
            CALL CHK1MAT( M, 8, N, 9, IC, JC, DESCC, 21, INFO )
         END IF
      END IF
C     
C     Check the blocking sizes for equivalence and being not to small
C     
      IF( INFO.EQ.0 ) THEN
         IF( DESCA( MB_ ).NE.DESCA( NB_ ) ) INFO = -(100*13 + MB_)
      END IF
      IF( INFO.EQ.0 ) THEN
         IF( DESCB( MB_ ).NE.DESCB( NB_ ) ) INFO = -(100*17 + MB_)
      END IF
      IF( INFO.EQ.0 .AND. SOLVE ) THEN
         IF( DESCC( MB_ ).NE.DESCA( MB_ ) ) INFO = -(100*21 + MB_)
      END IF
      IF( INFO.EQ.0 .AND. SOLVE ) THEN
         IF( DESCC( NB_ ).NE.DESCB( NB_ ) ) INFO = -(100*21 + NB_)
      END IF
C
      IF( INFO.EQ.0 ) THEN
         IF( M.NE.DESCA( MB_ ) .AND. DESCA( MB_ ).LT.6 .AND. 
     $        LSAME( ASCHUR, 'N' ) ) 
     $        INFO = -(100*13 + MB_)
         IF( N.NE.DESCB( MB_ ) .AND. DESCB( MB_ ).LT.6 .AND. 
     $        LSAME( BSCHUR, 'N' ) )
     $        INFO = -(100*17 + MB_)
      END IF
C
C     Check the values of ASCHUR and BSCHUR and set the logical values
C     A_SCHUR and B_SCHUR for future reference
C 
      IF( INFO.EQ.0 ) THEN
         IF( .NOT.( LSAME( ASCHUR, 'N' ) .OR. LSAME( ASCHUR, 'S' ))) 
     $        THEN
            INFO = -2
         ELSE
            IF( LSAME( ASCHUR, 'S' ) ) THEN
               A_SCHUR = .TRUE.
            ELSE
               A_SCHUR = .FALSE.
            END IF
         END IF
      END IF
C
      IF( INFO.EQ.0 ) THEN
         IF( .NOT.( LSAME( BSCHUR, 'N' ) .OR. LSAME( BSCHUR, 'S' ))) 
     $        THEN
            INFO = -3
         ELSE
            IF( LSAME( BSCHUR, 'S' ) ) THEN
               B_SCHUR = .TRUE.
            ELSE
               B_SCHUR = .FALSE.
            END IF  
         END IF  
      END IF
      
C
C     Check the values of TRANSA and TRANSB
C
      IF( INFO.EQ.0 .AND. SOLVE ) THEN
         IF( .NOT.( LSAME( TRANSA, 'N' ) .OR. LSAME( TRANSA, 'T' ))) 
     $        THEN
            INFO = -4
         END IF
         IF( INFO.EQ.0 ) THEN
            IF( .NOT.( LSAME( TRANSB, 'N' ) .OR. LSAME( TRANSB, 'T' ))) 
     $           THEN
               INFO = -5
            END IF  
         END IF
      END IF
C
C     Check the value of ISGN
C
      IF( INFO.EQ.0 .AND. SOLVE ) THEN
         IF( .NOT.(ISGN.EQ.1 .OR. ISGN.EQ.-1) ) THEN
            INFO = -6
         END IF
      END IF
C
C     Check the value of COMM
C
      IF( INFO.EQ.0 .AND. SOLVE ) THEN
         IF( .NOT.( LSAME( COMM, 'S' ) .OR. LSAME( COMM, 'D' ))) 
     $        THEN
            INFO = -7
         END IF
      END IF
C
C     Check the values of IA, JA, IB, JB, IC, JC
C
      IF( INFO.EQ.0 ) THEN
         IF( IA.LT.1 .OR. ( .NOT. A_SCHUR .AND. IA.GT.1 ) ) INFO = -11
         IF( JA.NE.IA ) INFO = -12
         IF( IB.LT.1 .OR. ( .NOT. B_SCHUR .AND. IB.GT.1 ) ) INFO = -15
         IF( JB.NE.IB ) INFO = -16
      END IF
      IF( INFO.EQ.0 .AND. SOLVE ) THEN
         IF( MOD(IC,DESCC(MB_)) .NE. MOD(IA,DESCA(MB_)) ) INFO = -19
         IF( MOD(JC,DESCC(NB_)) .NE. MOD(JB,DESCB(NB_)) ) INFO = -20
      END IF
C
C     Check if workspace query
C
      LQUERY = LDWORK.EQ.-1 .OR. LIWORK.EQ.-1
C     
C     Test working space
C     
C     Check the number of rows and columns of the smallest submatrices
C     of A and B including sub(A), sub(B) and sub(C) that conform with 
C     the ScaLAPACK conventions 
C     
      IF( INFO.EQ.0 .OR. LQUERY ) THEN
         IROFFA = MOD( IA - 1, DESCA(MB_) )
         ICOFFA = MOD( JA - 1, DESCA(NB_) ) 
         IROFFB = MOD( IB - 1, DESCB(MB_) ) 
         ICOFFB = MOD( JB - 1, DESCB(NB_) ) 
         IROFFC = MOD( IC - 1, DESCC(MB_) ) 
         ICOFFC = MOD( JC - 1, DESCC(NB_) )
         CALL INFOG2L( IA, JA, DESCA, NPROW, NPCOL, MYROW, MYCOL,
     $                 LIA, LJA, ARSRC, ACSRC )
         CALL INFOG2L( IB, JB, DESCB, NPROW, NPCOL, MYROW, MYCOL,
     $                 LIB, LJB, BRSRC, BCSRC )
         CALL INFOG2L( IC, JC, DESCC, NPROW, NPCOL, MYROW, MYCOL,
     $                 LIC, LJC, CRSRC, CCSRC )
         NROWS_A = NUMROC( M + IROFFA, DESCA(MB_), MYROW, ARSRC, NPROW )
         NROWS_B = NUMROC( N + IROFFB, DESCB(MB_), MYROW, BRSRC, NPROW )
         NROWS_C = NUMROC( M + ICOFFC, DESCC(MB_), MYROW, CRSRC, NPROW ) 
         NCOLS_A = NUMROC( M + ICOFFA, DESCA(NB_), MYCOL, ACSRC, NPCOL )
         NCOLS_B = NUMROC( N + ICOFFB, DESCB(NB_), MYCOL, BCSRC, NPCOL )
         NCOLS_C = NUMROC( N + ICOFFC, DESCC(NB_), MYCOL, CCSRC, NPCOL )
C
C     Initialize matrix descriptors for Z_A, Z_B and C_COPY
C
         IF( .NOT. A_SCHUR .AND. SOLVE ) THEN
            CALL DESCINIT( DESCZ_A, M + IROFFA, M + ICOFFA, DESCA(MB_), 
     $                     DESCA(NB_), ARSRC, ACSRC, ICTXT, 
     $                     MAX( 1, NROWS_A ), INFO )
         END IF
         IF( .NOT. B_SCHUR .AND. SOLVE ) THEN
            CALL DESCINIT( DESCZ_B, N + IROFFB, N + ICOFFB, DESCB(MB_), 
     $                     DESCB(NB_), BRSRC, BCSRC, ICTXT, 
     $                     MAX( 1, NROWS_B ), INFO )
         END IF
         IF( SOLVE .AND. (.NOT. A_SCHUR .OR. .NOT. B_SCHUR ) ) THEN
            CALL DESCINIT( DESCCC, M + IROFFC, N + ICOFFC, DESCC(MB_),
     $                     DESCC(NB_), CRSRC, CCSRC, ICTXT,
     $                     MAX( 1, NROWS_C ), INFO ) 
         END IF
C     
C     Check which matrix that has the greatest order
C     
         MAXMN = MAX( M, N )
C     
C     Compute the sizes of the different needed workspaces
C     
         SIZE_A = NROWS_A * NCOLS_A 
         SIZE_B = NROWS_B * NCOLS_B
         SIZE_C = NROWS_C * NCOLS_C
         SIZE_TAUA = NUMROC( M - 1 + IROFFA, DESCA(NB_), MYCOL, ACSRC, 
     $                       NPCOL )
         SIZE_TAUB = NUMROC( N - 1 + ICOFFA, DESCB(NB_), MYCOL, BCSRC, 
     $                       NPCOL )
         SIZE_WR = MAXMN
         SIZE_WI = MAXMN
C     
C     Do some workspace queries to the routines called below.
C     
         IF( .NOT. A_SCHUR ) THEN
            CALL PDGEHRD( M, 1, M, A, IA, JA, DESCA, DPDUM1, 
     $                    DPDUM2, WRKSPCQ, INFO )
            LWORK2 = INT( DPDUM2( 1 ) )
         ELSE
            LWORK2 = 0
         END IF
C     
         IF( .NOT. B_SCHUR ) THEN
            CALL PDGEHRD( N, 1, N, B, IB, JB, DESCB, DPDUM1, 
     $                    DPDUM2, WRKSPCQ, INFO )
            LWORK2 = MAX( LWORK2, INT( DPDUM2(1) ) )
         END IF
C     
         IF( .NOT. A_SCHUR .AND. SOLVE ) THEN
            CALL PDORMHR( 'L','N', M, M, 1, M, A, IA, JA, DESCA, 
     $                    DPDUM1, DWORK, IROFFA+1, ICOFFA+1, DESCZ_A, 
     $                    DPDUM2, WRKSPCQ, INFO )
            LWORK1 = INT( DPDUM2( 1 ) )
         ELSE
            LWORK1 = 0
         END IF
C         
         IF( .NOT. B_SCHUR .AND. SOLVE ) THEN
            CALL PDORMHR( 'L','N', N, N, 1, N, B, IB, JB, DESCB, 
     $                    DPDUM1, DWORK, IROFFB+1, ICOFFB+1, DESCZ_B, 
     $                    DPDUM2, WRKSPCQ, INFO ) 
            LWORK1 = MAX(LWORK1, INT( DPDUM2( 1 ) ) )
         END IF
C
C
         IF( .NOT. A_SCHUR ) THEN
            CALL PDLAHQR( WANTT, WANTZ, M, IA, IA+M-1, A, DESCA, 
     $                    DWORK, DWORK, IROFFA+1, IROFFA+M, DWORK, 
     $                    DESCZ_A, DPDUM1, WRKSPCQ, IDUM1, WRKSPCQ, 
     $                    INFO )
         ELSE
            DPDUM1(1) = ZERO
            IDUM1(1) = 0
         END IF
         IF( .NOT. B_SCHUR ) THEN
            CALL PDLAHQR( WANTT, WANTZ, N, IB, IB+N-1, B, DESCB, 
     $                    DWORK, DWORK, IROFFB+1, IROFFB+N, DWORK, 
     $                    DESCZ_B, DPDUM2, WRKSPCQ, IDUM2, WRKSPCQ, 
     $                    INFO )
         ELSE
            DPDUM2(1) = ZERO
            IDUM2(1) = 0
         END IF
         IF( .NOT. A_SCHUR. OR. .NOT. B_SCHUR ) THEN
            LWORK3 = MAX( INT( DPDUM1( 1 ) ),  INT( DPDUM2( 1 )) )
            IWRK = MAX( IDUM1(1), IDUM2(1) )
         ELSE
            LWORK3 = 0
            IWRK = 0
         END IF
C
         IF( SOLVE ) THEN
            CALL PTRSYDTD( TRANSA, TRANSB, ISGN, COMM, M, N, A, IA, JA, 
     $                     DESCA, B, IB, JB, DESCB, C, IC, JC, DESCC, 
     $                     MB2, DPDUM1, WRKSPCQ, IDUM1, LIWORK, 
     $                     NOEXSY, SCALE, INFO )
            LWORK4 = INT( DPDUM1( 1 ) )
            IWRK = MAX( IWRK, IDUM1( 1 ) )
         ELSE
            LWORK4 = 0
            IWRK = MAX( IWRK, 0 )
         END IF
C     
         LSWORK = MAX( LWORK4, MAX( LWORK2,
     $            MAX( LWORK1, LWORK3 ) ) ) 
C
C     Add together the síze of the whole DP WORK requirements
C     
         WRK = LSWORK
C     
         IF( .NOT. A_SCHUR )  WRK = WRK + SIZE_TAUA
         IF( .NOT. B_SCHUR )  WRK = WRK + SIZE_TAUB
         IF( .NOT. A_SCHUR .OR. .NOT. B_SCHUR ) THEN
            WRK = WRK + SIZE_WR + SIZE_WI
         END IF
         IF( SOLVE ) THEN
            IF( .NOT. A_SCHUR .OR. .NOT. B_SCHUR )
     $           WRK = WRK +  SIZE_C
            IF( .NOT. A_SCHUR )  WRK = WRK + SIZE_A
            IF( .NOT. B_SCHUR )  WRK = WRK + SIZE_B  
         END IF
C
C     Now check if the call to this routine was a workspace query
C     and check if the workspace supplied is enough. 
C     
         IF( WRK.GT.LDWORK .AND. .NOT.LQUERY ) THEN
            INFO = -24
         ELSEIF( IWRK.GT.LIWORK .AND. .NOT. LQUERY ) THEN 
            INFO = -26
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
         CALL PCHK1MAT( M, 8, M, 8, 1, 1, DESCA, 13, 0, IDUM1, 
     $                  IDUM2, INFO )
      END IF
      IF( INFO.EQ.0 ) THEN
         CALL PCHK1MAT( N, 9, N, 9, 1, 1, DESCB, 17, 0, IDUM1, 
     $                  IDUM2, INFO )
      END IF
      IF( INFO.EQ.0 .AND. SOLVE ) THEN
         CALL PCHK1MAT( M, 8, N, 9, 1, 1, DESCC, 21, 0, IDUM1,
     $                  IDUM2, INFO )
      END IF
C     
C     Checking if we may continue or should interrupt
C     
      IF( INFO.NE.0 ) THEN
         CALL PXERBLA( ICTXT, 'PGESYDTD', -INFO )
         RETURN
      END IF
C
C     Quick return if possible
C 
      IF( M.EQ.0 .OR. N.EQ.0 .OR. INFO.EQ.1 ) RETURN
C     
C     Divide the workspace between the different local pointers
C     which are needed
C 
      IF( SOLVE ) THEN
         IF( A_SCHUR .AND. B_SCHUR ) THEN
            SWORK = 1
         ELSEIF( A_SCHUR .AND. .NOT. B_SCHUR ) THEN
            Z_B = 1
            C_COPY = Z_B + SIZE_B
            TAUB = C_COPY + SIZE_C
            WR = TAUB + SIZE_TAUB
            WI = WR + SIZE_WR
            SWORK = WI + SIZE_WI
         ELSEIF( .NOT. A_SCHUR .AND. B_SCHUR ) THEN
            Z_A = 1
            C_COPY = Z_A + SIZE_A
            TAUA = C_COPY + SIZE_C
            WR = TAUA + SIZE_TAUA
            WI = WR + SIZE_WR
            SWORK = WI + SIZE_WI
         ELSEIF( .NOT. A_SCHUR .AND. .NOT. B_SCHUR ) THEN 
            Z_A = 1
            Z_B = Z_A + SIZE_A
            C_COPY = Z_B + SIZE_B
            TAUA = C_COPY + SIZE_C
            TAUB = TAUA + SIZE_TAUA
            WR = TAUB + SIZE_TAUB
            WI = WR + SIZE_WR
            SWORK = WI + SIZE_WI
         END IF 
      ELSE
         IF( A_SCHUR .AND. B_SCHUR ) THEN
            SWORK = 1
         ELSEIF( A_SCHUR .AND. .NOT. B_SCHUR ) THEN
            TAUB = 1
            WR = TAUB + SIZE_TAUB
            WI = WR + SIZE_WR
            SWORK = WI + SIZE_WI
         ELSEIF( .NOT. A_SCHUR .AND. B_SCHUR ) THEN
            TAUA = 1
            WR = TAUA + SIZE_TAUA
            WI = WR + SIZE_WR
            SWORK = WI + SIZE_WI
         ELSEIF( .NOT. A_SCHUR .AND. .NOT. B_SCHUR ) THEN 
            TAUA = 1
            TAUB = TAUA + SIZE_TAUA
            WR = TAUB + SIZE_TAUB
            WI = WR + SIZE_WR
            SWORK = WI + SIZE_WI
         END IF 
      END IF
      LSWORK = LDWORK-SWORK+1
C
C     Turn A and B into upper Hessenberg form 
C
      T1 = MPI_WTIME()
      IF( .NOT. A_SCHUR ) 
     $     CALL PDGEHRD( M, 1, M, A, IA, JA, DESCA, DWORK( TAUA ),
     $                   DWORK( SWORK ), LSWORK, INFO )
C
      IF( .NOT. B_SCHUR )
     $     CALL PDGEHRD( N, 1, N, B, IB, JB, DESCB, DWORK( TAUB ), 
     $                   DWORK( SWORK ), LSWORK, INFO )
C
C     Form Q_A and Q_B explicitely in Z_A and Z_B by using PDORMHR on
C     an identity matrix stored in Z_A and Z_B.
C
      IF( .NOT. A_SCHUR .AND. SOLVE )
     $     CALL PDLASET( 'A', M, M, ZERO, ONE, DWORK(Z_A), IROFFA+1, 
     $                   ICOFFA+1, DESCZ_A )
C
      IF( .NOT. B_SCHUR .AND. SOLVE )
     $     CALL PDLASET( 'A', N, N, ZERO, ONE, DWORK(Z_B), IROFFB+1, 
     $                   ICOFFB+1, DESCZ_B )
C
      IF( .NOT. A_SCHUR .AND. SOLVE )
     $     CALL PDORMHR( 'L','N', M, M, 1, M, A, IA, JA, DESCA, 
     $                   DWORK( TAUA ), DWORK( Z_A ), IROFFA+1, 
     $                   ICOFFA+1, DESCZ_A, DWORK( SWORK ), LSWORK, 
     $                   INFO )
C
      IF( .NOT. B_SCHUR .AND. SOLVE )
     $     CALL PDORMHR( 'L','N', N, N, 1, N, B, IB, JB, DESCB, 
     $                   DWORK( TAUB ), DWORK( Z_B ), IROFFB+1, 
     $                   ICOFFB+1, DESCZ_B, DWORK( SWORK ), LSWORK, 
     $                   INFO )
C
C     Extract the upper Hessenberg parts of A and B
C
      IF( .NOT. A_SCHUR ) THEN
         CALL PDLASET( 'Lower triangular', M-2, M-2, ZERO, ZERO, A, 
     $                 IA + 2, JA, DESCA )
      END IF
C     
      IF( .NOT. B_SCHUR ) THEN
         CALL PDLASET( 'Lower triangular', N-2, N-2, ZERO, ZERO, B, 
     $                 IB + 2, JB, DESCB )
      END IF
C     
C     Compute the real Schur forms of Hessenberg matrices A and B
C
      IF( .NOT. A_SCHUR )
     $     CALL PDLAHQR( WANTT, WANTZ, M, IA, IA+M-1, A, DESCA, 
     $                   DWORK( WR ), DWORK( WI ), IROFFA+1, IROFFA+M, 
     $                   DWORK( Z_A ), DESCZ_A, DWORK( SWORK ), LSWORK, 
     $                   IWORK, LIWORK, INFO )
C
      IF( .NOT. B_SCHUR )
     $     CALL PDLAHQR( WANTT, WANTZ, N, IB, IB+N-1, B, DESCB, 
     $                   DWORK( WR ), DWORK( WI ), IROFFB+1, IROFFB+N, 
     $                   DWORK( Z_B ), DESCZ_B, DWORK( SWORK ), LSWORK, 
     $                   IWORK, LIWORK, INFO )
C      
      T1 = MPI_WTIME() - T1
      IF( .NOT. SOLVE ) RETURN
C     
C     Update C with respect to the transformations done
C     
      T2 = MPI_WTIME()
      IF( .NOT. A_SCHUR ) THEN
         CALL PDLACPY( COPY_ALL, M, N, C, IC, JC, DESCC, 
     $                 DWORK( C_COPY ), IROFFC+1, ICOFFC+1, DESCCC )
         CALL PDGEMM( 'Transpose','No transpose', M, N, M, ONE, 
     $                DWORK( Z_A ), IROFFA+1, ICOFFA+1, DESCZ_A, 
     $                DWORK( C_COPY ), IROFFC+1, ICOFFC+1, DESCCC, 
     $                ZERO, C, IC, JC, DESCC )
      END IF
C     
      IF( .NOT. B_SCHUR ) THEN
         CALL PDLACPY( COPY_ALL, M, N, C, IC, JC, DESCC, 
     $                 DWORK( C_COPY ), IROFFC+1, ICOFFC+1, DESCCC )
         CALL PDGEMM( 'No transpose','No transpose', M, N, N, ONE, 
     $                DWORK( C_COPY ), IROFFC+1, ICOFFC+1, DESCCC, 
     $                DWORK( Z_B ), IROFFA+1, ICOFFA+1, DESCZ_B, 
     $                ZERO, C, IC, JC, DESCC )
      END IF
      T2 = MPI_WTIME() - T2
C     
C     Solve the reduced triangular problem
C     
      T3 = MPI_WTIME()
      CALL PTRSYDTD( TRANSA, TRANSB, ISGN, COMM, M, N, A, IA, JA, DESCA, 
     $               B, IB, JB, DESCB, C, IC, JC, DESCC, MB2, 
     $               DWORK( SWORK ), LSWORK, IWORK, LIWORK, NOEXSY, 
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
         CALL PDLACPY( COPY_ALL, M, N, C, IC, JC, DESCC, 
     $                 DWORK( C_COPY ), IROFFC+1, ICOFFC+1, DESCCC )
         CALL PDGEMM( 'No transpose','No transpose', M, N, M, ONE, 
     $                DWORK( Z_A ), IROFFA+1, ICOFFA+1, DESCZ_A, 
     $                DWORK( C_COPY ), IROFFC+1, ICOFFC+1, DESCCC, 
     $                ZERO, C, IC, JC, DESCC )
      END IF
C
      IF( .NOT. B_SCHUR ) THEN
         CALL PDLACPY( COPY_ALL, M, N, C, IC, JC, DESCC, 
     $                 DWORK( C_COPY ), IROFFC+1, ICOFFC+1, DESCCC ) 
         CALL PDGEMM( 'No transpose','Transpose', M, N, N, ONE, 
     $                DWORK( C_COPY ), IROFFC+1, ICOFFC+1, DESCCC,
     $                DWORK( Z_B ), IROFFA+1, ICOFFA+1, DESCZ_B, 
     $                ZERO, C, IC, JC, DESCC )
      END IF
      T2 = T2 + MPI_WTIME() - T22
C
      DWORK(1) = T1
      DWORK(2) = T2
      DWORK(3) = T3
C
      END
C
C     END OF PGESYDTD 
C
C *** Last line of PGESYDTD ***
