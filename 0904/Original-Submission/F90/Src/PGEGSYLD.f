CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C         This file is a part of the SCASY software library.           C
C         See http://www.cs.umu.se/parallel/research/scasy             C
C         Contributors: Robert Granat and Bo Kågström.                 C
C         SCASY Version 0.10,    March 31, 2006.                       C
C         Copyright 2006, Umeå University, Sweden.                     C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE PGEGSYLD( JOB, ACSCHR, BDSCHR, TRANAC, TRANBD, ISGN, 
     $                     COMM, M, N, A, IA, JA, DESCA, B, IB, JB, 
     $                     DESCB, C, IC, JC, DESCC, D, ID, JD, DESCD, 
     $                     E, IE, JE, DESCE, MB2, DWORK, LDWORK, 
     $                     IWORK, LIWORK, NOEXSY, SCALE, INFO )
C
C  -- ScaLAPACK-style routine (preliminary version ) --
C     HPC2N and Department of Computing Science,
C     University of Umeå, Sweden.
C     Written by Robert Granat, (granat@cs.umu.se)
C     January 14, 2006.
C
      IMPLICIT NONE
C
C     .. Scalar Arguments ..
      CHARACTER*1          JOB, ACSCHR, BDSCHR, TRANAC, TRANBD, 
     $                     COMM
      INTEGER              ISGN, M, N, IA, JA, IB, JB, IC, JC, ID, JD,
     $                     IE, JE, LDWORK, LIWORK, NOEXSY, INFO, MB2
      DOUBLE PRECISION     SCALE
C     ..
C     .. Array Arguments ..
      INTEGER              DESCA( * ), DESCD( * ), DESCB( * ), 
     $                     DESCE( * ), DESCC( * ), IWORK( * )
      DOUBLE PRECISION     A( * ), D( * ), B( * ), E( * ), C( * ), 
     $                     DWORK( * )
C 
C  Purpose
C  =======
C  The subroutine solves the general generalized Sylvester 
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
C            If JOB = 'S', the equation is solved. For JOB = 'R',
C            only the reduction step regarding the matrix pairs
C            (A,C) and (B,D) is performed.
C
C  ACSCHR    (global input) CHARACTER*1
C            If ACSCHR = 'S', then the matrix pair (A,C) is supposed to 
C            be in generalized Schur form. No reduction to the generalized
C            Schur form of this matrix pair will be done.
C            If ACSCHR = 'N', then the matrix pair (A,C) is not in 
C            generalized Schur form and a reduction will be performed.
C
C  BDSCHR    (global input) CHARACTER*1
C            If BDSCHR = 'S', then the matrix pair (B,D) is supposed to 
C            be in generalized Schur form. No reduction to the generalized
C            Schur form of this matrix pair will be done.
C            If BDSCHR = 'N', then the matrix pair (B,D) is not in 
C            generalized Schur form and a reduction will be performed.
C
C  TRANAC    (global input) CHARACTER*1
C            With JOB = 'S':
C              If TRANAC = 'N' then op(A) = A and op(C) = C
C              If TRANAC = 'T' then op(A) = A**T and op(C) = C**T
C            Otherwise, TRANAC is not referenced.
C
C  TRANBD    (global input) CHARACTER*1
C            With JOB = 'S':
C              If TRANBD = 'N' then op(B) = B and op(D) = D
C              If TRANBD = 'T' then op(B) = B**T and op(D) = D**T
C            Otherwise, TRANBD is not referenced.
C
C  ISGN      (global input) INTEGER*1
C            With JOB = 'S':
C              If ISGN = 1, we solve the equation with a '+'.
C              If ISGN = -1, we solve the equation with a '-'.
C            Otherwise, ISGN is not referenced.
C
C  Input/output arguments
C
C  COMM      (global input/output) CHARACTER*1
C            With JOB = 'S':
C            This subroutine uses two different communications schemes in
C            solving the reduced triangular problem:
C              If COMM = 'S', the "shifts" scheme is used.
C              If COMM = 'D', the "communicate on demand" scheme is used.
C            The choice COMM = 'S' is only valid for TRANAC = TRANBD = 'N' 
C            or TRANAC = TRANBD = 'T'. The scheme used will be output.
C            See the references for details.
C            Otherwise, COMM is not referenced.
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
C            pieces of the global distributed matrix A. On output,
C            the local pieces of the A in the generalized Schur form
C            of the matrix pair (A,C).
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
C            pieces of the global distributed matrix B. On output,
C            the local pieces of the B in the generalized Schur form
C            of the matrix pair (B,D).
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
C            pieces of the global distributed matrix C. On output,
C            the local pieces of the C in the generalized Schur form
C            of the matrix pair (A,C).
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
C            pieces of the global distributed matrix D. On output,
C            the local pieces of the D in the generalized Schur form
C            of the matrix pair (B,D).
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
C            With JOB = 'S': 
C            On entry E contains the local pieces of the global 
C            distributed matrix E. On exit, it contains the local 
C            pieces of the global distributed solution X.
C            Otherwise, E is not referenced.
C
C  IE        (global input) INTEGER
C            With JOB = 'S': Row start index for sub(E), i.e., the 
C            submatrix to operate on. MOD(IE,MB_A) = MOD(IA,MB_A) must 
C            hold.
C            Otherwise, IE is not referenced. 
C
C  JE        (global input) INTEGER
C            With JOB = 'S': Column start index for sub(E), i.e., the 
C            submatrix to operate on. MOD(JE,MB_B) = MOD(JB,MB_B) must 
C            hold.
C            Otherwise, IE is not referenced. 
C
C  DESCE     (global and local input) INTEGER array of dimension DLEN_.
C            With JOB = 'S', the array descriptor for the global 
C            distributed matrix E.
C            Otherwise, E is not referenced.
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
C            LIWORK >= DBA + DBB + MB_A + MB_B + 8, where 
C            DBA = ICEIL(LOCr(IA+IROFFA),MB_A) and DBB = 
C            = ICEIL(LOCr(IB+IROFFB),MB_B).
C           
C  Output information
C 
C  NOEXSY    (local output) INTEGER
C            With JOB = 'S':
C            When solving the triangular problem in PTRGSYLD it is possible
C            that we have to extend some subsystems to not lose any data
C            from some 2x2 block of conjugate pairs of eigenvalues. This
C            extension can cause bad performance. NOEXSY helps us to 
C            keep track of the number of such extensions.
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
C             If INFO = 2, (A,C) and (B,D) have common or very close 
C             eigenvalues; perturbed values were used to solve the 
C             equations (but (A,C) and (B,D) are unchanged).
C             If INFO = 3,  the problem is badly scaled - (C,F) should 
C             be scaled a factor SCALE before calling this routine too 
C             guarantee an overflow free solution, i.e., the solution
C             may well have overflowed.
C
C  Method
C  ======
C  This subroutine implements a parallel version of the Bartels-Stewart
C  method for solving the general generalized Sylvester equation (see 
C  the references for details).
C
C  Additional requirements
C  =======================
C
C  (A,C) and (B,D) must be distributed using the same blocking factor in 
C  each direction. Moreover, for E the blocksize in the row direction 
C  must agree with (A,C)'s, and the blocksize in the column direction must
C  agree with (B,D)'s. 
C
C  The matrix pairs (A,C) and  (B,D) must be aligned internally, i.e., A 
C  must be aligned with C, and so on.
C
C  Limitations
C  ===========
C  In contrary to RECSYs RECGSYL this routine do not scale against 
C  overflow in the solution. See SCALE and INFO.
C
C  For ACSCHR = 'N' or BDSCHR = 'N' it is not possible to work with
C  submatrices, e.g., ACSCHR = 'N' implies IA = IC = JA = JC = 1. This 
C  limitation comes from the fact that utilized PDHGEQZ does not work on 
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
C  Known issues
C  ============
C  Experienced hung problems outside this routine when PDORGQR was used 
C  with more that one cpu. Temporary fix: use PDORMQR instead. The 
C  choice is controlled by logical parameter GSYL_PDORGQR, see below.
C
C  Keywords
C  ========
C  QZ-algorithm, generalized Schur form, generalized Sylvester equation, 
C  Hessenberg-triangular transformation, PDGEMM-updates.
C
C  =====================================================================
C     ..
C     ..Parameters ..
      INTEGER            BLOCK_CYCLIC_2D, CSRC_, CTXT_, DLEN_, DTYPE_,
     $                   LLD_, MB_, M_, NB_, N_, RSRC_, NOINARG,
     $                   WRKQUE, MAXBULG
      LOGICAL            GSYL_PDORGQR
      DOUBLE PRECISION   ZERO, ONE, MINONE
      CHARACTER*1        CPYALL
      PARAMETER          ( BLOCK_CYCLIC_2D = 1, DLEN_ = 9, DTYPE_ = 1,
     $                     CTXT_ = 2, M_ = 3, N_ = 4, MB_ = 5, NB_ = 6,
     $                     RSRC_ = 7, CSRC_ = 8, LLD_ = 9, NOINARG = 27,
     $                     WRKQUE = -1, MAXBULG = 8,
     $                     GSYL_PDORGQR = .FALSE., ZERO = 0.0D0, 
     $                     ONE = 1.0D0, MINONE = -1.0D0, CPYALL = 'A')       
C     .. 
C     .. Local Scalars ..
      LOGICAL            SOLVE, LQUERY, SCHRAC, SCHRBD, TRANSAC, TRANSBD
      INTEGER            ICTXT, MYCOL, MYROW, NPCOL, NPROW, LLDF, DD,
     $                   NROWSAC, NROWSBD, NCOLSAC, NCOLSBD, NCOLSE,
     $                   LSWORK, SIZEAC, SIZEBD, SIZEE, MAXMB, LLDD,
     $                   SIZE_WI, SIZE_WR, LW1, LW2, LWTRI, BETA, LLDE,
     $                   WRK, LIWNEED,  MAXmn, DUMMY, NPROCS, I, IS, 
     $                   IX, JX, RSRC, CSRC, DBAD, DBBE, EIGI, EIGR, 
     $                   LW3, GI, GJ, J, LINFO, IROFFE, ICOFFE, LIE,
     $                   LJE, ERSRC, ECSRC, NROWSE, DBAC, DBBD, IWRK,
     $                   BULGES
      DOUBLE PRECISION   FNORM, ANORM, BNORM, XNORM, T1, T2, T22, T3
C     ..
C     .. Local Arrays and Local Pointers ..
      INTEGER             Q1, Q2, Z1, Z2, ECOPY, TAUC, TAUD, WR, 
     $                    WI, SWORK, IDUM1( 1 ), IDUM2( 1 ), DUMARR(1),
     $                    DESCQZ1( DLEN_ ), DESCQZ2( DLEN_ ),
     $                    DESCEC( DLEN_ )
      DOUBLE PRECISION    DPDUM( 1 )
C     ..
C     .. External Subroutines ..
C     
      EXTERNAL           BLACS_GRIDINFO, CHK1MAT, PCHK1MAT, PXERBLA,
     $                   PDGGHRD, PDHGEQZ, PDGEMM, PDLACPY, 
     $                   PDIDMAT, PTRGSYLD, PDGEQRF, PDGERQF,
     $                   PDORMQR, PDORGQR
C     .. 
C     .. External Functions ..
      DOUBLE PRECISION   MPI_WTIME
      INTEGER            NUMROC, ILCM, ICEIL, INDXL2G
      LOGICAL            LSAME
      EXTERNAL           NUMROC, ILCM, ICEIL, LSAME, MPI_WTIME, INDXL2G
      
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
         ELSE
            SOLVE = .FALSE.
         END IF
      END IF
C     
C     Check dimensions.
C     
      IF( INFO.EQ.0 ) THEN
         CALL CHK1MAT( M, 8, M, 8, IA, JA, DESCA, 13, INFO )
         IF( INFO.EQ.0 )
     $         CALL CHK1MAT( M, 8, M, 8, IC, JC, DESCC, 21, INFO )
         IF( INFO.EQ.0 ) THEN
            CALL CHK1MAT( N, 9, N, 9, IB, JB, DESCB, 17, INFO )
            IF( INFO.EQ.0 )
     $           CALL CHK1MAT( N, 9, N, 9, ID, JD, DESCD, 25, INFO )   
         END IF
         IF( INFO.EQ.0 .AND. SOLVE ) THEN
            CALL CHK1MAT( M, 8, N, 9, IE, JE, DESCE, 29, INFO )
         END IF
      END IF
C     
C     Check the blocking sizes for equivalence and being not to small
C     
      IF( INFO.EQ.0 ) THEN
         IF( DESCA( MB_ ).NE.DESCA( NB_ ) ) INFO = -(100*13 + MB_)
      END IF
      IF( INFO.EQ.0 ) THEN
         IF( DESCC( MB_ ).NE.DESCC( NB_ ) ) INFO = -(100*21 + MB_)
         IF( DESCC( MB_ ).NE.DESCA( MB_ ) ) INFO = -(100*21 + MB_)
      END IF
      IF( INFO.EQ.0 ) THEN
         IF( DESCB( MB_ ).NE.DESCB( NB_ ) ) INFO = -(100*17 + MB_)
      END IF
      IF( INFO.EQ.0 ) THEN
         IF( DESCD( MB_ ).NE.DESCD( NB_ ) ) INFO = -(100*25 + MB_)
         IF( DESCD( MB_ ).NE.DESCB( MB_ ) ) INFO = -(100*25 + MB_) 
      END IF
      IF( INFO.EQ.0 .AND. SOLVE ) THEN
         IF( DESCE( MB_ ).NE.DESCA( MB_ ) ) INFO = -(100*29 + MB_)
      END IF
      IF( INFO.EQ.0 .AND. SOLVE ) THEN
         IF( DESCE( NB_ ).NE.DESCB( NB_ ) ) INFO = -(100*29 + NB_)
      END IF
C
      IF( INFO.EQ.0 ) THEN
         IF( M.NE.DESCA( MB_ ) .AND. DESCA( MB_ ).LT.9 .AND. 
     $        LSAME( ACSCHR, 'N' ) ) 
     $        INFO = -(100*13 + MB_)
         IF( N.NE.DESCB( MB_ ) .AND. DESCB( MB_ ).LT.9 .AND. 
     $        LSAME( BDSCHR, 'N' ) )
     $        INFO = -(100*17 + MB_)
      END IF
C
C     Check alignment, part 1: check RSRC_ and CSRC_
C
      IF( INFO.EQ.0 ) THEN
         IF( DESCC( RSRC_ ) .NE. DESCA( RSRC_ ) ) 
     $        INFO = -(100*21 + RSRC_)
         IF( DESCC( CSRC_ ) .NE. DESCA( CSRC_ ) ) 
     $        INFO = -(100*21 + CSRC_)
         IF( DESCD( RSRC_ ) .NE. DESCB( RSRC_ ) ) 
     $        INFO = -(100*25 + RSRC_)
         IF( DESCD( CSRC_ ) .NE. DESCB( CSRC_ ) ) 
     $        INFO = -(100*25 + CSRC_)
      END IF
C
C     Check the values of ACSCHR and BDSCHR and set the logical values
C     SCHRAC and SCHRBD for future reference
C 
      IF( INFO.EQ.0 ) THEN
         IF( .NOT.( LSAME( ACSCHR, 'N' ) .OR. LSAME( ACSCHR, 'S' ))) 
     $        THEN
            INFO = -2
         ELSE
            IF( LSAME( ACSCHR, 'S' ) ) THEN
               SCHRAC = .TRUE.
            ELSE
               SCHRAC = .FALSE.
            END IF
         END IF
      END IF
C
      IF( INFO.EQ.0 ) THEN
         IF( .NOT.( LSAME( BDSCHR, 'N' ) .OR. LSAME( BDSCHR, 'S' ))) 
     $        THEN
            INFO = -3
         ELSE
            IF( LSAME( BDSCHR, 'S' ) ) THEN
               SCHRBD = .TRUE.
            ELSE
               SCHRBD = .FALSE.
            END IF  
         END IF  
      END IF
C
C     Check the values of TRANAC and TRANSB
C
      IF( INFO.EQ.0 .AND. SOLVE ) THEN
         IF( .NOT.( LSAME( TRANAC, 'N' ) .OR. LSAME( TRANAC, 'T' ))) 
     $        THEN
            INFO = -4
         ELSE
            TRANSAC = LSAME( TRANAC, 'T' )
         END IF
         IF( INFO.EQ.0 ) THEN
            IF( .NOT.( LSAME( TRANBD, 'N' ) .OR. LSAME( TRANBD, 'T' ))) 
     $           THEN
               INFO = -5
            ELSE
               TRANSBD = LSAME( TRANBD, 'T' )
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
C     Check alignment, part 2: check IA, JA, IB, JB, IC, JC etc.
C
      IF( INFO.EQ.0 ) THEN
         IF( IA.LT.1 .OR. ( .NOT. SCHRAC .AND. IA.GT.1 ) ) INFO = -11
         IF( JA.NE.IA ) INFO = -12
         IF( IB.LT.1 .OR. ( .NOT. SCHRBD .AND. IB.GT.1 ) ) INFO = -15
         IF( JB.NE.IB ) INFO = -16
         IF( IC .NE. IA ) INFO = -19
         IF( JC .NE. JA ) INFO = -20
         IF( ID .NE. IB ) INFO = -23
         IF( JD .NE. JB ) INFO = -24
      END IF
      IF( INFO.EQ.0 .AND. SOLVE ) THEN
         IF( MOD(IE,DESCE(MB_)) .NE. MOD(IA,DESCA(MB_)) ) INFO = -27
         IF( MOD(JE,DESCE(NB_)) .NE. MOD(JB,DESCB(NB_)) ) INFO = -28
      END IF
C
C     Check if workspace query
C
      LQUERY = LDWORK.EQ.-1 .OR. LIWORK.EQ.-1
C     
C     Test working space
C     
C     Check the number of rows and columns of the smallest submatrix
C     of E including sub(E) that conform with the ScaLAPACK conventions 
C     and compute the local number of rows and columns of (A,C) and (B,D).
C
      IF( INFO.EQ.0 .OR. LQUERY ) THEN 
         IROFFE = MOD( IE - 1, DESCE(MB_) ) 
         ICOFFE = MOD( JE - 1, DESCE(NB_) )
         CALL INFOG2L( IE, JE, DESCE, NPROW, NPCOL, MYROW, MYCOL,
     $                 LIE, LJE, ERSRC, ECSRC )
C         
         NROWSAC = NUMROC( DESCA(M_), DESCA(MB_), MYROW, DESCA(RSRC_), 
     $                     NPROW )
         NROWSBD = NUMROC( DESCB(M_), DESCB(MB_), MYROW, DESCB(RSRC_), 
     $                     NPROW )
         NROWSE  = NUMROC( M + IROFFE, DESCE(MB_), MYROW, ERSRC, NPROW )
         NCOLSAC = NUMROC( DESCA(N_), DESCA(NB_), MYCOL, DESCA(CSRC_), 
     $                     NPCOL )
         NCOLSBD = NUMROC( DESCB(N_), DESCB(NB_), MYCOL, DESCB(CSRC_), 
     $                     NPCOL )
         NCOLSE  = NUMROC( N + ICOFFE, DESCE(NB_), MYCOL, ECSRC, NPCOL )
C
C     Initialize matrix descriptor for orthogonal transformations
C     Q1, Q2, Z1, Z2 and the copies of C and F
C
         IF( .NOT. SCHRAC .AND. SOLVE ) THEN
            CALL DESCINIT( DESCQZ1, DESCA(M_), DESCA(N_), DESCA(MB_), 
     $                     DESCA(NB_), DESCA(RSRC_), DESCA(CSRC_), 
     $                     ICTXT, MAX( 1, NROWSAC ), INFO )
         END IF
         IF( .NOT. SCHRBD .AND. SOLVE ) THEN
            CALL DESCINIT( DESCQZ2, DESCB(M_), DESCB(N_), DESCB(MB_), 
     $                     DESCB(NB_), DESCB(RSRC_), DESCB(CSRC_), 
     $                     ICTXT, MAX( 1, NROWSBD ), INFO )
         END IF
         IF( SOLVE .AND. (.NOT. SCHRAC .OR. .NOT. SCHRBD ) ) THEN
            CALL DESCINIT( DESCEC, M + IROFFE, N + ICOFFE, DESCE(MB_),
     $                     DESCE(NB_), ERSRC, ECSRC, ICTXT,
     $                     MAX( 1, NROWSE ), INFO ) 
         END IF         
C     
C     Check which matrix of A and B that has the greatest order for
C     future reference - also check the greatest blocking size
C     
         MAXMN = MAX( M, N )
         MAXMB = MAX( DESCA(MB_), DESCB(MB_) )
C     
C     Compute the sizes of the different needed DP workspaces
C     
         SIZEAC = NROWSAC * NCOLSAC 
         SIZEBD = NROWSBD * NCOLSBD
         SIZEE  = NROWSE  * NCOLSE
C     
C     For the LSWORK we need to check the workspace needed by work-
C     space quiries to each individual routines. First, for PDGGHRD...
C        
         IF( .NOT. SCHRAC ) THEN
            CALL PDGGHRD( 'Vectors', 'Vectors', M, IA, IA+M-1, A, DESCA, 
     $           C, DESCC, A, DESCA, C, DESCC, DPDUM, -1, LINFO )
            LW1 = MAX( LW1, INT( DPDUM(1) ) )
         ELSE
            LW1 = 0
         END IF
         IF( .NOT. SCHRBD ) THEN
            CALL PDGGHRD( 'Vectors', 'Vectors', N, IB, IB+N-1, B, DESCB, 
     $           D, DESCD, B, DESCB, D, DESCD, DPDUM, -1, LINFO )
            LW1 = MAX( LW1, INT( DPDUM(1) ) )
         END IF
C     
C     Then for PDHGEQZ...
C     
         IF( .NOT. SCHRAC ) THEN
            BULGES = MAX( 1, MIN( MAXBULG, ICEIL( M, DESCA(MB_) ) - 1 ))
            CALL PDHGEQZ( 'Schur', 'Vectors', 'Vectors', M, IA, IA+M-1, 
     $           A, DESCA, C, DESCC, DWORK, DWORK, DWORK, A, DESCA, 
     $           C, DESCC, IA, IA+M-1, IA, IA+M-1, BULGES, DPDUM, 
     $           -1, LINFO )
            LW2 = MAX( LW2, INT( DPDUM(1) ) )
         ELSE
            LW2 = 0
         END IF
         IF( .NOT. SCHRBD ) THEN 
            BULGES = MAX( 1, MIN( MAXBULG, ICEIL( N, DESCB(MB_) ) - 1 ))
            CALL PDHGEQZ( 'Schur', 'Vectors', 'Vectors', N, IB, IB+N-1, 
     $           B, DESCB, D, DESCD, DWORK, DWORK, DWORK, B, DESCB, 
     $           D, DESCD, IB, IB+N-1, IB, IB+N-1, BULGES, DPDUM, 
     $           -1, LINFO )
            LW2 = MAX( LW2, INT( DPDUM(1) ) )
         END IF
C
C     Then for PDGEQRF, PDORMQR and PDORGQR
C
         IF( .NOT. SCHRAC .OR. .NOT. SCHRBD ) THEN
            CALL PDGEQRF( MAXMN, MAXMN, C, IC, JC, DESCC, DWORK(TAUC), 
     $           DPDUM, -1, LINFO )
            LW3 = INT( DPDUM( 1 ) )
         ELSE
            LW3 = 0
         END IF
         IF( .NOT. SCHRAC .OR. .NOT. SCHRBD .AND. SOLVE ) THEN
            CALL PDORMQR( 'Left', 'NoTranspose', MAXMN, MAXMN, MAXMN, 
     $           C, IC, JC, DESCC, DWORK(TAUC), A, IA, JA, DESCA,
     $           DPDUM, -1, LINFO )
            LW3 = MAX( INT( DPDUM( 1 ) ), LW3 )
            CALL PDORMQR( 'Right', 'Transpose', MAXMN, MAXMN, MAXMN, 
     $           C, IC, JC, DESCC, DWORK(TAUC), A, IA, JA, DESCA, 
     $           DPDUM, -1, LINFO )
            LW3 = MAX( INT( DPDUM( 1 ) ), LW3 )
            CALL PDORGQR( MAXMN, MAXMN, MAXMN, C, IC, JC, DESCC, 
     $           DWORK(TAUC), DPDUM, -1, LINFO )
            LW3 = MAX( INT( DPDUM( 1 ) ), LW3 )
         ELSE
            LW3 = 0
         END IF
C 
C     Now compute the workspace needed in PTRGSYLD using a workspace
C     query as above
C
         IF( SOLVE ) THEN
            CALL PTRGSYLD( TRANAC, TRANBD, ISGN, COMM, M, N, A, IA, JA, 
     $           DESCA, B, IB, JB, DESCB, C, IC, JC, DESCC, D, 
     $           ID, JD, DESCD, E, IE, JE, DESCE, MB2, DPDUM, -1, 
     $           IDUM1, LIWORK, NOEXSY, SCALE, LINFO )
            LWTRI = INT( DPDUM( 1 ) )
            IWRK = IDUM1( 1 )
         ELSE
            LWTRI = 0
            IWRK = 0
         END IF 
C     
C     Take now the maximum 
C     
         LSWORK = MAX( LWTRI, MAX( LW1, MAX( LW2, LW3 ) ) )
C     
C     Add together the síze of the whole DP WORK requirements
C     depending on the given forms of (A,C) and (B,D)
C     
         WRK = LSWORK + 3 * MAXMN
C     
         IF( SOLVE ) THEN
            IF( .NOT. SCHRAC )  WRK = WRK + 2*SIZEAC
            IF( .NOT. SCHRBD )  WRK = WRK + 2*SIZEBD
            IF( .NOT. SCHRAC .OR. .NOT. SCHRBD ) THEN
               WRK = WRK + SIZEE 
            END IF
         END IF
C     
C     Check the integer workspace
C     
         DBAC = ICEIL( M, DESCA(MB_) )
         DBBD = ICEIL( N, DESCB(MB_) )
         IWRK = MAX(IWRK,DBAC + DBBD + DESCA(MB_) + DESCB(MB_) + 8) 
C
C     Now check if the call to this routine was a workspace query
C     and check if the workspace supplied is enough. 
C     
         IF( WRK.GT.LDWORK .AND. .NOT.LQUERY ) THEN
            INFO = -32
         ELSEIF( IWRK.GT.LIWORK .AND. .NOT. LQUERY ) THEN 
            INFO = -34
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
         CALL PCHK1MAT( M, 8, M, 8, IA, JA, DESCA, 13, 0, IDUM1, 
     $                  IDUM2, INFO )
         IF( INFO.EQ.0 ) 
     $        CALL PCHK1MAT( M, 8, M, 8, IC, JC, DESCC, 21, 0, IDUM1, 
     $                       IDUM2, INFO )
      END IF
      IF( INFO.EQ.0 ) THEN
         CALL PCHK1MAT( N, 9, N, 9, IB, JB, DESCB, 17, 0, IDUM1, 
     $                  IDUM2, INFO )
         IF( INFO.EQ.0 )
     $        CALL PCHK1MAT( N, 9, N, 9, ID, JD, DESCD, 25, 0, IDUM1, 
     $                       IDUM2, INFO )    
      END IF
      IF( INFO.EQ.0 .AND. SOLVE ) THEN
         CALL PCHK1MAT( M, 8, N, 9, IE, JE, DESCE, 29, 0, IDUM1,
     $                  IDUM2, INFO )
      END IF
C     
C     Checking if we may continue or should interrupt
C     
      IF( INFO.NE.0 ) THEN
         CALL PXERBLA( ICTXT, 'PGEGSYLD', -INFO )
         RETURN
      END IF
C
C     Quick return if possible
C 
      IF( M.EQ.0 .OR. N.EQ.0 ) RETURN
C     
C     Divide the workspace between the different local pointers
C     which are needed - that depends on the given forms of 
C     (A,C) and (B,D).
C
      IF( SOLVE ) THEN
         IF( SCHRAC .AND. SCHRBD ) THEN
            SWORK = 1
         ELSEIF( SCHRAC .AND. (.NOT. SCHRBD) ) THEN
            Q2 = 1
            Z2 = Q2 + SIZEBD
            ECOPY = Z2 + SIZEBD
            EIGR = ECOPY + SIZEE
            EIGI = EIGR + MAXMN
            BETA = EIGI + MAXMN
            TAUD = BETA + MAXMN
            SWORK = TAUD + NCOLSBD
         ELSEIF( (.NOT. SCHRAC) .AND. SCHRBD ) THEN
            Q1 = 1
            Z1 = Q1 + SIZEAC
            ECOPY = Z1 + SIZEAC
            EIGR = ECOPY + SIZEE
            EIGI = EIGR + MAXMN
            BETA = EIGI + MAXMN
            TAUC = BETA + MAXMN
            SWORK = TAUC + NCOLSAC
         ELSEIF( (.NOT. SCHRAC) .AND. (.NOT. SCHRBD) ) THEN 
            Q1 = 1
            Z1 = Q1 + SIZEAC 
            Q2 = Z1 + SIZEAC
            Z2 = Q2 + SIZEBD
            ECOPY = Z2 + SIZEBD
            EIGR = ECOPY + SIZEE
            EIGI = EIGR + MAXMN
            BETA = EIGI + MAXMN
            TAUC = BETA + MAXMN
            TAUD = TAUC + NCOLSAC
            SWORK = TAUD + NCOLSBD
         END IF 
      ELSE
         IF( SCHRAC .AND. SCHRBD ) THEN
            SWORK = 1
         ELSEIF( SCHRAC .AND. (.NOT. SCHRBD) ) THEN
            EIGR = 1
            EIGI = EIGR + MAXMN
            BETA = EIGI + MAXMN
            TAUD = BETA + MAXMN
            SWORK = TAUD + NCOLSBD
         ELSEIF( (.NOT. SCHRAC) .AND. SCHRBD ) THEN
            EIGR = 1
            EIGI = EIGR + MAXMN
            BETA = EIGI + MAXMN
            TAUC = BETA + MAXMN
            SWORK = TAUC + NCOLSAC
         ELSEIF( (.NOT. SCHRAC) .AND. (.NOT. SCHRBD) ) THEN 
            EIGR = 1
            EIGI = EIGR + MAXMN
            BETA = EIGI + MAXMN
            TAUC = BETA + MAXMN
            TAUD = TAUC + NCOLSAC
            SWORK = TAUD + NCOLSBD
         END IF 
      END IF
      LSWORK = LDWORK-SWORK+1
C
C     Start the solution process - turn (A,C) and (B,D) into
C     Hessenberg-triangular form
C
C     QR-factorize C and D
C
      T1 = MPI_WTIME()
      IF( .NOT. SCHRAC ) 
     $     CALL PDGEQRF( M, M, C, IC, JC, DESCC, DWORK(TAUC), 
     $                   DWORK(SWORK), LSWORK, LINFO )
      IF( .NOT. SCHRBD )
     $     CALL PDGEQRF( N, N, D, ID, JD, DESCD, DWORK(TAUD), 
     $                   DWORK(SWORK), LSWORK, LINFO )
C
C     Update matrices A and B with respect to the QR-factorizations
C
      IF( .NOT. SCHRAC )
     $     CALL PDORMQR( 'Left', 'Transpose', M, M, M, C, IC, JC, DESCC, 
     $                   DWORK(TAUC), A, IA, JA, DESCA, DWORK(SWORK),
     $                   LSWORK, LINFO )
      IF( .NOT. SCHRBD )
     $     CALL PDORMQR( 'Left', 'Transpose', N, N, N, D, ID, JD, DESCD, 
     $                   DWORK(TAUD), B, IB, JB, DESCB, DWORK(SWORK),
     $                   LSWORK, LINFO )
C
C     Extract matrices Q1 and Q2 from implicit storage in
C     the matrices D and E. Also initialize the matrice Z1 and Z2 
C     as identity matrices of order M and N, respectively. Notice:
C     since we experienced problems with PDORGQR, we do this in
C     two different fashions depending on the logical parameter
C     GSYL_PDORGQR. See also "Known issues" above.
C
      IF( GSYL_PDORGQR ) THEN
         IF( .NOT. SCHRAC .AND. SOLVE ) THEN
            CALL PDLACPY( 'All', M, M, C, IC, JC, DESCC, DWORK(Q1), 
     $                    IA, JA, DESCQZ1 )
            CALL PDORGQR( M, M, M, DWORK(Q1), IA, JA, DESCQZ1, 
     $                    DWORK(TAUC), DWORK(SWORK), LSWORK, LINFO )
            CALL PDLASET( 'All', M, M, ZERO, ONE, DWORK(Z1), IA, JA, 
     $                    DESCQZ1 )
         END IF
         IF( .NOT. SCHRBD .AND. SOLVE ) THEN
            CALL PDLACPY( 'All', N, N, D, ID, JD, DESCD, DWORK(Q2), 
     $                    IB, JB, DESCQZ2 )
            CALL PDORGQR( N, N, N, DWORK(Q2), IB, JB, DESCQZ2, 
     $                    DWORK(TAUD), DWORK(SWORK), LSWORK, LINFO )
            CALL PDLASET( 'All', N, N, ZERO, ONE, DWORK(Z2), IB, JB, 
     $                    DESCQZ2 )
         END IF
      ELSE
         IF( .NOT. SCHRAC .AND. SOLVE ) THEN
            CALL PDLASET( 'All', M, M, ZERO, ONE, DWORK(Q1), IA, JA, 
     $                    DESCQZ1 )
            CALL PDORMQR( 'Left', 'NoTranspose', M, M, M, C, IC, JC, 
     $                    DESCC, DWORK(TAUC), DWORK(Q1), IA, JA, 
     $                    DESCQZ1, DWORK(SWORK), LSWORK, LINFO ) 
            CALL PDLASET( 'All', M, M, ZERO, ONE, DWORK(Z1), IA, JA, 
     $                    DESCQZ1 )
         END IF
         IF( .NOT. SCHRBD .AND. SOLVE ) THEN
            CALL PDLASET( 'All', N, N, ZERO, ONE, DWORK(Q2), IB, JB, 
     $                    DESCQZ2 )
            CALL PDORMQR( 'Left', 'NoTranspose', N, N, N, D, ID, JD, 
     $                    DESCD, DWORK(TAUD), DWORK(Q2), IB, JB, 
     $                    DESCQZ2, DWORK(SWORK), LSWORK, LINFO ) 
            CALL PDLASET( 'All', N, N, ZERO, ONE, DWORK(Z2), IB, JB, 
     $                    DESCQZ2 )
         END IF
      END IF
C    
C     Extract upper triangular matrices from C and D
C
      IF( .NOT. SCHRAC ) THEN
         CALL PDLASET( 'Lower triangular', M-1, M-1, ZERO, ZERO, C, 
     $                 IC+1, JC, DESCC )
      END IF
      IF( .NOT. SCHRBD ) THEN
         CALL PDLASET( 'Lower triangular', N-1, N-1, ZERO, ZERO, D, 
     $                 ID+1, JD, DESCD )
      END IF
C
C     Call Hessenberg-triangular reduction routine from the 
C     full-triangular pairs (A,C) and (B,D)
C      
      IF( .NOT. SCHRAC .AND. SOLVE ) THEN  
         CALL PDGGHRD( 'Vectors', 'Vectors', M, IA, IA+M-1, A, DESCA, C, 
     $                 DESCC, DWORK(Q1), DESCQZ1, DWORK(Z1), DESCQZ1, 
     $                 DWORK(SWORK), LSWORK, LINFO )
      ELSEIF( .NOT. SCHRAC ) THEN
         CALL PDGGHRD( 'NoVectors', 'NoVectors', M, IA, IA+M-1, A, 
     $                 DESCA, C, DESCC, DWORK(Q1), DESCQZ1, DWORK(Z1), 
     $                 DESCQZ1, DWORK(SWORK), LSWORK, LINFO )
      END IF
C
      IF( .NOT. SCHRBD .AND. SOLVE ) THEN
         CALL PDGGHRD( 'Vectors', 'Vectors', N, IB, IB+N-1, B, DESCB, D, 
     $                 DESCD, DWORK(Q2), DESCQZ2, DWORK(Z2), DESCQZ2, 
     $                 DWORK(SWORK), LSWORK, LINFO )
      ELSEIF( .NOT. SCHRBD ) THEN
         CALL PDGGHRD( 'NoVectors', 'NoVectors', N, IB, IB+N-1, B, 
     $                 DESCB, D, DESCD, DWORK(Q2), DESCQZ2, DWORK(Z2), 
     $                 DESCQZ2, DWORK(SWORK), LSWORK, LINFO )
      END IF
C 
C    Apply the QZ-algorithm to compute the generalized Schur forms
C    (A,C) and (B,E)
C
      IF( .NOT. SCHRAC .AND. SOLVE ) THEN
         BULGES = MAXBULG
         CALL PDHGEQZ( 'Schur', 'Vectors', 'Vectors', M, IA, IA+M-1, A, 
     $                 DESCA, C, DESCC, DWORK(EIGR), DWORK(EIGI), 
     $                 DWORK(BETA), DWORK(Q1), DESCQZ1, DWORK(Z1), 
     $                 DESCQZ1, IA, IA+M-1, IA, IA+M-1, BULGES, 
     $                 DWORK(SWORK), LSWORK, LINFO )
      ELSEIF( .NOT. SCHRAC ) THEN
         BULGES = MAXBULG
         CALL PDHGEQZ( 'Schur', 'NoVectors', 'NoVectors', M, IA, IA+M-1,
     $                 A, DESCA, C, DESCC, DWORK(EIGR), DWORK(EIGI), 
     $                 DWORK(BETA), DWORK(Q1), DESCQZ1, DWORK(Z1), 
     $                 DESCQZ1, IA, IA+M-1, IA, IA+M-1, BULGES, 
     $                 DWORK(SWORK), LSWORK, LINFO )
      END IF
C
      IF( .NOT. SCHRBD .AND. SOLVE ) THEN
         BULGES = MAXBULG
         CALL PDHGEQZ( 'Schur', 'Vectors', 'Vectors', N, IB, IB+N-1, B,
     $                 DESCB, D, DESCD, DWORK(EIGR), DWORK(EIGI), 
     $                 DWORK(BETA), DWORK(Q2), DESCQZ2, DWORK(Z2), 
     $                 DESCQZ2, IB, IB+N-1, IB, IB+N-1, BULGES, 
     $                 DWORK(SWORK), LSWORK, LINFO )
      ELSEIF( .NOT. SCHRBD ) THEN
         BULGES = MAXBULG
         CALL PDHGEQZ( 'Schur', 'NoVectors', 'NoVectors', N, IB, IB+N-1,
     $                 B, DESCB, D, DESCD, DWORK(EIGR), DWORK(EIGI), 
     $                 DWORK(BETA), DWORK(Q2), DESCQZ2, DWORK(Z2), 
     $                 DESCQZ2, IB, IB+N-1, IB, IB+N-1, BULGES, 
     $                 DWORK(SWORK), LSWORK, LINFO )
      END IF
      T1 = MPI_WTIME() - T1
      IF( .NOT. SOLVE ) RETURN
C     
C     Update E with respect to the Schur decompositions
C     
      T2 = MPI_WTIME()
      IF( .NOT. SCHRAC ) THEN
         IF( .NOT. TRANSAC ) THEN
            CALL PDLACPY( 'All', M, N, E, IE, JE, DESCE, DWORK( ECOPY ),
     $                    1 + IROFFE, 1 + ICOFFE, DESCEC )
            CALL PDGEMM( 'T', 'N', M, N, M, ONE, DWORK( Q1 ), IA, JA, 
     $                   DESCQZ1, DWORK( ECOPY ), 1 + IROFFE, 
     $                   1 + ICOFFE, DESCEC, ZERO, E, IE, JE, DESCE )
         ELSE
            CALL PDLACPY( 'All', M, N, E, IE, JE, DESCE, DWORK( ECOPY ),
     $                    1 + IROFFE, 1 + ICOFFE, DESCEC )
            CALL PDGEMM( 'T', 'N', M, N, M, ONE, DWORK( Z1 ), IA, JA, 
     $                   DESCQZ1, DWORK( ECOPY ), 1 + IROFFE, 
     $                   1 + ICOFFE, DESCEC, ZERO, E, IE, JE, DESCE )
         END IF
      END IF
C     
      IF( .NOT. SCHRBD ) THEN
         IF( .NOT. TRANSBD ) THEN
            CALL PDLACPY( 'All', M, N, E, IE, JE, DESCE, DWORK( ECOPY ),
     $                    1 + IROFFE, 1 + ICOFFE, DESCEC )
            CALL PDGEMM( 'N', 'N', M, N, N, ONE, DWORK( ECOPY ), 
     $                   1 + IROFFE, 1 + ICOFFE, DESCEC, DWORK( Z2 ), 
     $                   IB, JB, DESCQZ2, ZERO, E, IE, JE, DESCE )
         ELSE
            CALL PDLACPY( 'All', M, N, E, IE, JE, DESCE, DWORK( ECOPY ),
     $                    1 + IROFFE, 1 + ICOFFE, DESCEC )
            CALL PDGEMM( 'N', 'N', M, N, N, ONE, DWORK( ECOPY ), 
     $                   1 + IROFFE, 1 + ICOFFE, DESCEC, DWORK( Q2 ), 
     $                   IB, JB, DESCQZ2, ZERO, E, IE, JE, DESCE )
         END IF
      END IF
      T2 = MPI_WTIME() - T2
C     
C     Now we have reduced our general equation to a (quasi-)triangular
C     case. Solve now the reduced problem with a call to a triangular 
C     solver - the solution to that reduced problem is stored
C     on E. 
C
      T3 = MPI_WTIME()
      CALL PTRGSYLD( TRANAC, TRANBD, ISGN, COMM, M, N, A, IA, JA, DESCA,
     $               B, IB, JB, DESCB, C, IC, JC, DESCC, D, ID, JD, 
     $               DESCD, E, IE, JE, DESCE, MB2, DWORK( SWORK ), 
     $               LSWORK, IWORK, LIWORK, NOEXSY, SCALE, LINFO )
      IF( INFO.NE.0 ) THEN
         IF( LINFO.EQ.1 ) INFO = 2
         IF( LINFO.EQ.2 ) INFO = 3
      END IF
      T3 = MPI_WTIME() - T3
C     
C     Transform the solution back to the original coordinate system
C
      T22 = MPI_WTIME()
      IF( .NOT. SCHRAC ) THEN
         IF( .NOT. TRANSAC ) THEN
            CALL PDLACPY( 'All', M, N, E, IE, JE, DESCE, DWORK( ECOPY ),
     $                    1 + IROFFE, 1 + ICOFFE, DESCEC )
            CALL PDGEMM( 'N','N', M, N, M, ONE, DWORK( Z1 ), IA, JA, 
     $                   DESCQZ1, DWORK( ECOPY ), 1 + IROFFE, 
     $                   1 + ICOFFE, DESCEC, ZERO, E, IE, JE, DESCE )
         ELSE
            CALL PDLACPY( 'All', M, N, E, IE, JE, DESCE, DWORK( ECOPY ),
     $                    1 + IROFFE, 1 + ICOFFE, DESCEC )
            CALL PDGEMM( 'N', 'N', M, N, M, ONE, DWORK( Q1 ), IA, JA, 
     $                   DESCQZ1, DWORK( ECOPY ), 1 + IROFFE, 
     $                   1 + ICOFFE, DESCEC, ZERO, E, IE, JE, DESCE )
         END IF 
      END IF
C
      IF( .NOT. SCHRBD ) THEN
         IF( .NOT. TRANSBD ) THEN
            CALL PDLACPY( 'All', M, N, E, IE, JE, DESCE, DWORK( ECOPY ),
     $                    1 + IROFFE, 1 + ICOFFE, DESCEC ) 
            CALL PDGEMM( 'N', 'T', M, N, N, ONE, DWORK( ECOPY ), 
     $                   1 + IROFFE, 1 + ICOFFE, DESCEC, DWORK( Q2 ), 
     $                   IA, JA, DESCQZ2, ZERO, E, IE, JE, DESCE )
         ELSE
            CALL PDLACPY( 'All', M, N, E, IE, JE, DESCE, DWORK( ECOPY ),
     $                    1 + IROFFE, 1 + ICOFFE, DESCE ) 
            CALL PDGEMM( 'N', 'T', M, N, N, ONE, DWORK( ECOPY ), 
     $                   1 + IROFFE, 1 + ICOFFE, DESCEC, DWORK( Z2 ), 
     $                   IB, JB, DESCQZ2, ZERO, E, IE, JE, DESCE )
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
C     END OF PGEGSYLD 
C
C *** Last line of PGEGSYLD ***
