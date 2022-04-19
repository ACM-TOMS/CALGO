CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C         This file is a part of the SCASY software library.           C
C         See http://www.cs.umu.se/parallel/research/scasy             C
C         Contributors: Robert Granat and Bo Kågström.                 C
C         SCASY Version 0.10,    March 31, 2006.                       C
C         Copyright 2006, Umeå University, Sweden.                     C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE PGEGCSYD( JOB, TRANZ, ADSCHR, BESCHR, TRANAD, TRANBE, 
     $                     ISGN, COMM, M, N, A, IA, JA, DESCA, B, IB, 
     $                     JB, DESCB, C, IC, JC, DESCC, D, ID, JD, 
     $                     DESCD, E, IE, JE, DESCE, F, IF, JF, DESCF, 
     $                     MBNB2, DWORK, LDWORK, IWORK, LIWORK, NOEXSY, 
     $                     SCALE, INFO )
C
C  -- ScaLAPACK-style routine (preliminary version ) --
C     HPC2N and Department of Computing Science,
C     University of Umeå, Sweden.
C     Written by Robert Granat, (granat@cs.umu.se)
C     January 10, 2006.
C
      IMPLICIT NONE
C
C     .. Scalar Arguments ..
      CHARACTER*1          JOB, TRANZ, ADSCHR, BESCHR, TRANAD, TRANBE,
     $                     COMM
      INTEGER              ISGN, M, N, IA, JA, IB, JB, IC, JC, ID, JD,
     $                     IE, JE, IF, JF, LDWORK, LIWORK, NOEXSY, INFO
      DOUBLE PRECISION     SCALE
C     ..
C     .. Array Arguments ..
      INTEGER              DESCA( * ), DESCD( * ), DESCB( * ), 
     $                     DESCE( * ), DESCC( * ), DESCF( * ), 
     $                     IWORK( * ), MBNB2( 2 )
      DOUBLE PRECISION     A( * ), D( * ), B( * ), E( * ), C( * ), 
     $                     F( * ), DWORK( * )
C 
C  Purpose
C  =======
C  The subroutine solves the general generalized coupled Sylvester 
C  Equation (GCSY)
C
C    ( op( sub( A ) ) * X +/- Y * op( sub( B ) ) = C,
C      op( sub( D ) ) * X +/- Y * op( sub( E ) ) = F )   (1)
C
C  where sub(A) and sub(D) are M-by-M matrices, sub(B) and sub(E) are 
C  N-by-N matrices and the solution matrices X and Y are M-by-N which
C  overwrites C and F. 
C  
C  The notation op(_) means the transpose or non-transpose of a matrix. 
C  Each matrix in the pairs (A,D) and (B,E) are either both transposed
C  or not.
C
C  In matrix notation (1) (without transposes and ISGN = -1) is 
C  equivalent to solve the system Zx = scale b, where Z is defined as
C
C             Z = [ kron(In, A)  -kron(B', Im) ]         (2)
C                 [ kron(In, D)  -kron(E', Im) ].
C
C  Here Ik is the identity matrix of size k and X' is the transpose of
C  X. kron(X, Y) is the Kronecker product between the matrices X and Y.
C
C  If TRANZ = 'T', this routine solves the transposed system 
C  Z'*y = scale*b, which is equivalent to solve for X and Y in
C
C              A' * X  +  D' * Y  = scale *  C           (3)
C              X  * B' +  Y  * E' = scale * ( +/- F )
C
C  This case (TRANZ = 'T') is used to compute a one-norm-based estimate
C  of Dif[(A,D),(B,E)], the separation between the matrix pairs (A,D)
C  and (B,E), using PDLACON in the routine PGCSYCON.
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
C            If JOB = 'S', the equation is solved. For JOB = 'R',
C            only the reduction step regarding the matrix pairs
C            (A,D) and (B,E) is performed.
C
C  TRANZ     (global input) CHARACTER*1
C            With JOB = 'S':
C            If TRANZ = 'T', this routines solves for the equation (2).
C            Notice that TRANZ = 'T' is only allowed for 
C            TRANAD = TRANBE = 'N' or TRANAD = TRANBE = 'T'.
C            If TRANZ = 'N', this routines solves equation (1), and all
C            possible combination of the transposes TRANAD and TRANBE
C            are allowed.
C            Otherwise, TRANZ is not referenced.
C
C  ADSCHR    (global input) CHARACTER*1
C            If ADSCHR = 'S', then the matrix pair (A,D) is supposed to 
C            be in generalized Schur form. No reduction to the generalized
C            Schur form of this matrix pair will be done.
C            If ADSCHR = 'N', then the matrix pair (A,D) is not in 
C            generalized Schur form and a reduction will be performed.
C
C  BESCHR    (global input) CHARACTER*1
C            If BESCHR = 'S', then the matrix pair (B,E) is supposed to 
C            be in generalized Schur form. No reduction to the generalized
C            Schur form of this matrix pair will be done.
C            If BESCHR = 'N', then the matrix pair (B,E) is not in 
C            generalized Schur form and a reduction will be performed.
C
C  TRANAD    (global input) CHARACTER*1
C            With JOB = 'S':
C              If TRANAD = 'N' then op(A) = A and op(D) = D
C              If TRANAD = 'T' then op(A) = A**T and op(D) = D**T
C            Otherwise, TRANAD is not referenced.
C
C  TRANBE    (global input) CHARACTER*1
C            With JOB = 'S':
C              If TRANBE = 'N' then op(B) = B and op(E) = E
C              If TRANBE = 'T' then op(B) = B**T and op(E) = E**T
C            Otherwise, TRANBE is not referenced.
C
C  ISGN      (global input) INTEGER*1
C            With JOB = 'S':
C            If ISGN = 1, we solve the equation (1) with a '+'.
C            If ISGN = -1, we solve the equation (1) with a '-'.
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
C            The choice COMM = 'S' is only valid for TRANAD = TRANBE = 'N' 
C            or TRANAD = TRANBE = 'T'. The scheme used will be output.
C            See the references for details.
C            Otherwise, COMM is not referenced.
C
C  M         (global input) INTEGER
C            The number of rows and columns of the global distributed 
C            matrices A and D. This is also the number of rows of the
C            global distributed matrices C and F. M >= 0.
C
C  N         (global input) INTEGER
C            The number of rows and columns of the global distributed
C            matrices B and E. This is also the number of columns of the
C            global distributed matrices C and F. N >= 0.
C
C  A         (local input/output) DOUBLE PRECISION array 
C            Array of dimension (LLD_A,LOCc(M)). Contains the local
C            pieces of the global distributed matrix A. On output,
C            the local pieces of the A in the generalized Schur form
C            of the matrix pair (A,D).
C
C  IA        (global input) INTEGER
C            Row start index for sub(A), i.e., the submatrix to operate
C            on. IA >= 1.
C
C  JA        (global input) INTEGER
C            Column start index for sub(A), i.e., the submatrix to 
C            operate on. JA = IA must hold.
C
C
C  DESCA     (global and local input) INTEGER array of dimension DLEN_.
C            The array descriptor for the global distributed matrix A.
C
C  B         (local input/output) DOUBLE PRECISION array 
C            Array of dimension (LLD_B,LOCc(N)). Contains the local
C            pieces of the global distributed matrix B. On output,
C            the local pieces of the B in the generalized Schur form
C            of the matrix pair (B,E).
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
C            With JOB = 'S': 
C            On entry C contains the local pieces of the global 
C            distributed matrix C. On exit, it contains the local 
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
C  D         (local input/output) DOUBLE PRECISION array 
C            Array of dimension (LLD_D,LOCc(M)). Contains the local
C            pieces of the global distributed matrix D. On output,
C            the local pieces of the D in the generalized Schur form
C            of the matrix pair (A,D).
C
C  ID        (global input) INTEGER
C            Row start index for sub(D), i.e., the submatrix to operate
C            on. ID = IA must hold.
C
C  JD        (global input) INTEGER
C            Column start index for sub(D), i.e., the submatrix to 
C            operate on. JD = JA must hold. 
C
C  DESCD     (global and local input) INTEGER array of dimension DLEN_.
C            The array descriptor for the global distributed matrix D.
C
C  E         (local input/output) DOUBLE PRECISION array 
C            Array of dimension (LLD_E,LOCc(N)). Contains the local
C            pieces of the global distributed matrix E. On output,
C            the local pieces of the E in the generalized Schur form
C            of the matrix pair (B,E).
C
C  IE        (global input) INTEGER
C            Row start index for sub(E), i.e., the submatrix to operate
C            on. IE = IB must hold.
C
C  JE        (global input) INTEGER
C            Column start index for sub(E), i.e., the submatrix to 
C            operate on. JE = JB must hold.
C
C  DESCE     (global and local input) INTEGER array of dimension DLEN_.
C            The array descriptor for the global distributed matrix E.
C
C  F         (local input/local output) DOUBLE PRECISION array 
C            Array of dimension (LLD_F,LOCc(N)). 
C            With JOB = 'S': 
C            On entry F contains the local pieces of the global 
C            distributed matrix F. On exit, it contains the local 
C            pieces of the global distributed solution Y.
C            Otherwise, F is not referenced.
C 
C  IF        (global input) INTEGER
C            With JOB='S': Row start index for sub(F), i.e., the 
C            submatrix to operate on. IF = IC must hold.
C            Otherwise, IF is not referenced. 
C
C  JF        (global input) INTEGER
C            With JOB = S: Column start index for sub(F), i.e., the 
C            submatrix to operate on. JF = JC must hold.
C            Otherwise, JF is not referenced.  
C
C  DESCF     (global and local input) INTEGER array of dimension DLEN_.
C            With JOB = 'S', the array descriptor for the global 
C            distributed matrix F.
C            Otherwise, DESCF is not referenced.
C
C  MBNB2     (global input) INTEGER array of dimension 2.
C            Internal blocking factors for pipelining of subsolutions
C            for updates of the matrix C in PTRGCSYD (see the references 
C            for details).
C            1 < = MBNB2(1) <= DESCC( MB_ ) and 
C            1 < = MBNB2(2) <= DESCC( NB_ ) must hold.
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
C            When solving the triangular problem in PTRGCSYD it is possible
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
C             If INFO = 2, (A,D) and (B,E) have common or very close 
C             eigenvalues; perturbed values were used to solve the 
C             equations (but (A,D) and (B,E) are unchanged).
C             If INFO = 3,  the problem is badly scaled - C should have 
C             been scaled a factor SCALE before calling this routine to 
C             guarantee an overflow free solution, i.e., the solution
C             may well have overflowed.
C
C  Method
C  ======
C  This subroutine implements a parallel version of the Bartels-Stewart
C  method for solving the general generalized coupled Sylvester equation 
C  (see the references for details).
C
C  Additional requirements
C  =======================
C
C  (A,D) and (B,E) must be distributed using the same blocking factor in 
C  each direction. Moreover, for (C,F) the blocksize in the row direction 
C  must agree with (A,D)'s, and the blocksize in the column direction must
C  agree with (B,E)'s. 
C
C  The three matrix pairs (A,D), (B,E) and (C,F) must be aligned
C  internally, i.e., A must be aligned with D, and so on.
C
C  Limitations
C  ===========
C  In contrary to LAPACK DTGSYL this routine do not scale against 
C  overflow in the solution. See SCALE and INFO.
C
C  For ADSCHR = 'N' or BESCHR = 'N' it is not possible to work with
C  submatrices, e.g., ADSCHR = 'N' implies IA = ID = JA = JD = 1. This 
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
C  choice is controlled by logical parameter GCSY_PDORGQR, see below.
C
C  Keywords
C  ========
C  QZ-algorithm, generalized Schur form, generalized coupled Sylvester 
C  equation, Hessenberg-triangular transformation, PDGEMM-updates.
C
C  =====================================================================
C     ..
C     ..Parameters ..
      INTEGER            BLOCK_CYCLIC_2D, CSRC_, CTXT_, DLEN_, DTYPE_,
     $                   LLD_, MB_, M_, NB_, N_, RSRC_, NOINARG,
     $                   WRKQUE, MAXBULG
      LOGICAL            GCSY_PDORGQR
      DOUBLE PRECISION   ZERO, ONE, MINONE
      CHARACTER*1        CPYALL
      PARAMETER          ( BLOCK_CYCLIC_2D = 1, DLEN_ = 9, DTYPE_ = 1,
     $                     CTXT_ = 2, M_ = 3, N_ = 4, MB_ = 5, NB_ = 6,
     $                     RSRC_ = 7, CSRC_ = 8, LLD_ = 9, NOINARG = 28,
     $                     WRKQUE = -1, MAXBULG = 8, 
     $                     GCSY_PDORGQR = .FALSE., ZERO = 0.0D0, 
     $                     ONE = 1.0D0, MINONE = -1.0D0, CPYALL = 'A')         
C     .. 
C     .. Local Scalars ..
      LOGICAL            SOLVE, LQUERY, SCHRAD, SCHRBE, TRANSAD, 
     $                   TRANSBE, TRANSZ
      INTEGER            ICTXT, MYCOL, MYROW, NPCOL, NPROW, LLDF, DD,
     $                   NROWSAD, NROWSBE, NCOLSAD, NCOLSBE, NCOLSCF,
     $                   LSWORK, SIZEAD, SIZEBE, SIZECF, MAXMB, LLDD,
     $                   SIZE_WI, SIZE_WR, LW1, LW2, LWTRI, BETA, LLDE,
     $                   WRK, LIWNEED,  MAXmn, DUMMY, NPROCS, I, IS, 
     $                   IX, JX, RSRC, CSRC, DBAD, DBBE, EIGI, EIGR, 
     $                   LW3, GI, GJ, J, LINFO, IROFFCF, ICOFFCF, LIC,
     $                   LJC, CFRSRC, CFCSRC, NROWSCF, CRSRC, CCSRC, 
     $                   IWRK, BULGES
      DOUBLE PRECISION   FNORM, ANORM, BNORM, XNORM, T1, T2, T22, T3
C     ..
C     .. Local Arrays and Local Pointers ..
      INTEGER             Q1, Q2, Z1, Z2, CFCOPY, TAUD, TAUE, WR, 
     $                    WI, SWORK, IDUM1( 1 ), IDUM2( 1 ), DUMARR(1),
     $                    DESCQZ1( DLEN_ ), DESCQZ2( DLEN_ ),
     $                    DESCCF( DLEN_ )
      DOUBLE PRECISION    DPDUM( 1 )
C     ..
C     .. External Subroutines ..
C     
      EXTERNAL           BLACS_GRIDINFO, CHK1MAT, PCHK1MAT, PXERBLA,
     $                   PDGGHRD, PDHGEQZ, PDGEMM, PDLACPY, 
     $                   PDIDMAT, PTRGCSYD, PDGEQRF, PDGERQF,
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
      NPROCS = NPROW*NPCOL
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
C     Check dimensions
C     
      IF( INFO.EQ.0 ) THEN
         CALL CHK1MAT( M, 9, M, 9, IA, JA, DESCA, 14, INFO )
         IF( INFO.EQ.0 ) 
     $        CALL CHK1MAT( M, 9, M, 9, ID, JD, DESCD, 26, INFO )
         IF( INFO.EQ.0 ) THEN
            CALL CHK1MAT( N, 10, N, 10, IB, JB, DESCB, 18, INFO )
            IF( INFO.EQ.0 )
     $           CALL CHK1MAT( N, 10, N, 10, IE, JE, DESCE, 30, INFO )     
         END IF
         IF( INFO.EQ.0 .AND. SOLVE ) THEN
            CALL CHK1MAT( M, 9, N, 10, IC, JC, DESCC, 22, INFO )
            IF( INFO.EQ.0 )
     $           CALL CHK1MAT( M, 9, N, 10, IF, JF, DESCF, 34, INFO )   
         END IF
      END IF
C     
C     Check the blocking sizes for equivalence and being not to small
C     
      IF( INFO.EQ.0 ) THEN
         IF( DESCA( MB_ ).NE.DESCA( NB_ ) ) INFO = -(100*14 + MB_)
      END IF
      IF( INFO.EQ.0 ) THEN
         IF( DESCD( MB_ ).NE.DESCD( NB_ ) ) INFO = -(100*26 + MB_)
         IF( DESCD( MB_ ).NE.DESCA( MB_ ) ) INFO = -(100*26 + MB_)
      END IF
      IF( INFO.EQ.0 ) THEN
         IF( DESCB( MB_ ).NE.DESCB( NB_ ) ) INFO = -(100*18 + MB_)
      END IF
      IF( INFO.EQ.0 ) THEN
         IF( DESCE( MB_ ).NE.DESCE( NB_ ) ) INFO = -(100*30 + MB_)
         IF( DESCE( MB_ ).NE.DESCB( MB_ ) ) INFO = -(100*30 + MB_) 
      END IF
      IF( INFO.EQ.0 .AND. SOLVE ) THEN
         IF( DESCC( MB_ ).NE.DESCA( MB_ ) ) INFO = -(100*22 + MB_)
      END IF
      IF( INFO.EQ.0 .AND. SOLVE ) THEN
         IF( DESCC( NB_ ).NE.DESCB( NB_ ) ) INFO = -(100*22 + NB_)
      END IF
      IF( INFO.EQ.0 .AND. SOLVE ) THEN
         IF( DESCF( MB_ ).NE.DESCA( MB_ ) ) INFO = -(100*34 + MB_)
      END IF
      IF( INFO.EQ.0 .AND. SOLVE ) THEN
         IF( DESCF( NB_ ).NE.DESCB( NB_ ) ) INFO = -(100*34 + NB_)
      END IF
C
      IF( INFO.EQ.0 ) THEN
         IF( M.NE.DESCA( MB_ ) .AND. DESCA( MB_ ).LT.9 .AND. 
     $        LSAME( ADSCHR, 'N' ) ) 
     $        INFO = -(100*14 + MB_)
         IF( N.NE.DESCB( MB_ ) .AND. DESCB( MB_ ).LT.9 .AND. 
     $        LSAME( BESCHR, 'N' ) )
     $        INFO = -(100*18 + MB_)
      END IF
C
C     Check alignment, part 1: check RSRC_ and CSRC_
C
      IF( INFO.EQ.0 ) THEN
         IF( DESCD( RSRC_ ) .NE. DESCA( RSRC_ ) ) 
     $        INFO = -(100*26 + RSRC_)
         IF( DESCD( CSRC_ ) .NE. DESCA( CSRC_ ) ) 
     $        INFO = -(100*26 + CSRC_)
         IF( DESCE( RSRC_ ) .NE. DESCB( RSRC_ ) ) 
     $        INFO = -(100*30 + RSRC_)
         IF( DESCE( CSRC_ ) .NE. DESCB( CSRC_ ) ) 
     $        INFO = -(100*30 + CSRC_)
         IF( DESCF( RSRC_ ) .NE. DESCC( RSRC_ ) ) 
     $        INFO = -(100*34 + RSRC_)
         IF( DESCF( CSRC_ ) .NE. DESCC( CSRC_ ) ) 
     $        INFO = -(100*34 + CSRC_)
      END IF
C
C     Check the value of TRANZ
C     
      IF( INFO.EQ.0 .AND. SOLVE ) THEN
         IF( .NOT. ( LSAME( TRANZ, 'N') .OR. LSAME( TRANZ, 'T'))) THEN
            INFO = -2
         ELSE
            TRANSZ = LSAME( TRANZ, 'T' )
         END IF
      END IF
C
C     Check the values of ADSCHR and BESCHR and set the logical values
C     SCHRAD and SCHRBE for future reference
C 
      IF( INFO.EQ.0 ) THEN
         IF( .NOT.( LSAME( ADSCHR, 'N' ) .OR. LSAME( ADSCHR, 'S' ))) 
     $        THEN
            INFO = -3
         ELSE
            IF( LSAME( ADSCHR, 'S' ) ) THEN
               SCHRAD = .TRUE.
            ELSE
               SCHRAD = .FALSE.
            END IF
         END IF
      END IF
C
      IF( INFO.EQ.0 ) THEN
         IF( .NOT.( LSAME( BESCHR, 'N' ) .OR. LSAME( BESCHR, 'S' ))) 
     $        THEN
            INFO = -4
         ELSE
            IF( LSAME( BESCHR, 'S' ) ) THEN
               SCHRBE = .TRUE.
            ELSE
               SCHRBE = .FALSE.
            END IF  
         END IF  
      END IF
C
C     Check the values of TRANAD and TRANSB
C
      IF( INFO.EQ.0 .AND. SOLVE ) THEN
         IF( .NOT.( LSAME( TRANAD, 'N' ) .OR. LSAME( TRANAD, 'T' ))) 
     $        THEN
            INFO = -5
         ELSE
            TRANSAD = LSAME( TRANAD, 'T' )
         END IF
         IF( INFO.EQ.0 ) THEN
            IF( .NOT.( LSAME( TRANBE, 'N' ) .OR. LSAME( TRANBE, 'T' ))) 
     $           THEN
               INFO = -6
            ELSE
               TRANSBE = LSAME( TRANBE, 'T' )
            END IF  
         END IF
      END IF
C
C     Check the value of TRANZ combined with TRANAD and TRANBE
C
      IF( SOLVE ) THEN
         IF( TRANBE.NE.TRANAD .AND. TRANSZ ) THEN
            INFO = -2
         END IF
      END IF
C
C     Check the value of ISGN
C
      IF( INFO.EQ.0 .AND. SOLVE ) THEN
         IF( .NOT.(ISGN.EQ.1 .OR. ISGN.EQ.-1) ) THEN
            INFO = -7
         END IF
      END IF
C
C     Check the value of COMM
C
      IF( INFO.EQ.0 .AND. SOLVE ) THEN
         IF( .NOT.( LSAME( COMM, 'S' ) .OR. LSAME( COMM, 'D' ))) 
     $        THEN
            INFO = -8
         END IF
      END IF
C
C     Check alignment, part 2: check IA, JA, IB, JB, IC, JC etc.
C
      IF( INFO.EQ.0 ) THEN
         IF( IA.LT.1 .OR. ( .NOT. SCHRAD .AND. IA.GT.1 ) ) INFO = -12
         IF( JA.NE.IA ) INFO = -13
         IF( IB.LT.1 .OR. ( .NOT. SCHRBE .AND. IB.GT.1 ) ) INFO = -16
         IF( JB.NE.IB ) INFO = -17
         IF( ID .NE. IA ) INFO = -24
         IF( JD .NE. JA ) INFO = -25
         IF( IE .NE. IB ) INFO = -28
         IF( JE .NE. JB ) INFO = -29
      END IF
      IF( INFO.EQ.0 .AND. SOLVE ) THEN
         IF( MOD(IC,DESCC(MB_)) .NE. MOD(IA,DESCA(MB_)) ) INFO = -20
         IF( MOD(JC,DESCC(NB_)) .NE. MOD(JB,DESCB(NB_)) ) INFO = -21
         IF( IF .NE. IC ) INFO = -32
         IF( JF .NE. JC ) INFO = -33
      END IF
C
C     Check if workspace query
C
      LQUERY = LDWORK.EQ.-1 .OR. LIWORK.EQ.-1
C     
C     Test working space
C     
C     Check the number of rows and columns of the smallest submatrix
C     of (C,F) including sub(C,F) that conform with the ScaLAPACK 
C     conventions and compute the local number of rows and columns of
C     (A,D) and (B,E).
C
      IF( INFO.EQ.0 .OR. LQUERY ) THEN 
         IROFFCF = MOD( IC - 1, DESCC(MB_) ) 
         ICOFFCF = MOD( JC - 1, DESCC(NB_) )
         CALL INFOG2L( IC, JC, DESCC, NPROW, NPCOL, MYROW, MYCOL,
     $                 LIC, LJC, CFRSRC, CFCSRC )
C         
         NROWSAD = NUMROC( DESCA(M_), DESCA(MB_), MYROW, DESCA(RSRC_), 
     $                     NPROW )
         NROWSBE = NUMROC( DESCB(M_), DESCB(MB_), MYROW, DESCB(RSRC_), 
     $                     NPROW )
         NROWSCF = NUMROC( M + IROFFCF, DESCC(MB_), MYROW, CFRSRC, 
     $                     NPROW )
         NCOLSAD = NUMROC( DESCA(N_), DESCA(NB_), MYCOL, DESCA(CSRC_), 
     $                     NPCOL )
         NCOLSBE = NUMROC( DESCB(N_), DESCB(NB_), MYCOL, DESCB(CSRC_), 
     $                     NPCOL )
         NCOLSCF = NUMROC( N + ICOFFCF, DESCC(NB_), MYCOL, CFCSRC, 
     $                     NPCOL )
C
C     Initialize matrix descriptor for orthogonal transformations
C     Q1, Q2, Z1, Z2 and the copies of C and F
C
         IF( .NOT. SCHRAD .AND. SOLVE ) THEN
            CALL DESCINIT( DESCQZ1, DESCA(M_), DESCA(N_), DESCA(MB_), 
     $                     DESCA(NB_), DESCA(RSRC_), DESCA(CSRC_), 
     $                     ICTXT, MAX( 1, NROWSAD ), INFO )
         END IF
         IF( .NOT. SCHRBE .AND. SOLVE ) THEN
            CALL DESCINIT( DESCQZ2, DESCB(M_), DESCB(N_), DESCB(MB_), 
     $                     DESCB(NB_), DESCB(RSRC_), DESCB(CSRC_), 
     $                     ICTXT, MAX( 1, NROWSBE ), INFO )
         END IF
         IF( SOLVE .AND. (.NOT. SCHRAD .OR. .NOT. SCHRBE ) ) THEN
            CALL DESCINIT( DESCCF, M + IROFFCF, N + ICOFFCF, DESCC(MB_),
     $                     DESCC(NB_), CFRSRC, CFCSRC, ICTXT,
     $                     MAX( 1, NROWSCF ), INFO ) 
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
         SIZEAD = NROWSAD * NCOLSAD 
         SIZEBE = NROWSBE * NCOLSBE
         SIZECF = NROWSAD * NCOLSCF
C     
C     For the LSWORK we need to check the workspace needed by work-
C     space quiries to each individual routines. First, for PDGGHRD...
C     
         IF( .NOT. SCHRAD ) THEN
            CALL PDGGHRD( 'Vectors', 'Vectors', M, IA, IA+M-1, A, DESCA, 
     $           D, DESCD, A, DESCA, D, DESCD, DPDUM, -1, LINFO )
            LW1 = INT( DPDUM(1) ) 
         ELSE
            LW1 = 0
         END IF
C     
         IF( .NOT. SCHRBE ) THEN
            CALL PDGGHRD( 'Vectors', 'Vectors', N, IB, IB+N-1, B, DESCB, 
     $           E, DESCE, B, DESCB, E, DESCE, DPDUM, -1, LINFO )
            LW1 = MAX( LW1, INT( DPDUM(1) ) )
         END IF
C     
C     Then for PDHGEQZ...
C     
         IF( .NOT. SCHRAD ) THEN
            BULGES = MAX( 1, MIN( MAXBULG, ICEIL( M, DESCA(MB_) ) - 1 ))
            CALL PDHGEQZ( 'Schur', 'Vectors', 'Vectors', M, IA, IA+M-1, 
     $           A, DESCA, D, DESCD, DWORK, DWORK, DWORK, A, DESCA, 
     $           D, DESCD, IA, IA+M-1, IA, IA+M-1, BULGES, DPDUM, 
     $           -1, LINFO )
            LW2 = INT( DPDUM(1) )
         ELSE
            LW2 = 0
         END IF
C
         IF( .NOT. SCHRBE ) THEN
            BULGES = MAX( 1, MIN( MAXBULG, ICEIL( N, DESCB(MB_) ) - 1 ))
            CALL PDHGEQZ( 'Schur', 'Vectors', 'Vectors', N, IB, IB+N-1, 
     $           B, DESCB, E, DESCE, DWORK, DWORK, DWORK, B, DESCB, 
     $           E, DESCE, IB, IB+N-1, IB, IB+N-1, BULGES, DPDUM, 
     $           -1, LINFO )
            LW2 = MAX( LW2, INT( DPDUM(1) ) )
         END IF
C
C     Then for PDGEQRF, PDORMQR and PDORGQR
C
         IF( .NOT. SCHRAD .OR. .NOT. SCHRBE ) THEN
            CALL PDGEQRF( MAXMN, MAXMN, D, ID, JD, DESCD, DWORK, 
     $           DPDUM, -1, LINFO )
            LW3 = INT( DPDUM( 1 ) )
         ELSE
            LW3 = 0
         END IF
         IF( .NOT. SCHRAD .OR. .NOT. SCHRBE .AND. SOLVE ) THEN
            CALL PDORMQR( 'Left', 'NoTranspose', MAXMN, MAXMN, MAXMN, 
     $           D, ID, JD, DESCD, DWORK, A, IA, JA, DESCA,
     $           DPDUM, -1, LINFO )
            LW3 = MAX( INT( DPDUM( 1 ) ), LW3 )
            CALL PDORMQR( 'Right', 'Transpose', MAXMN, MAXMN, MAXMN, 
     $           D, ID, JD, DESCD, DWORK, A, IA, JA, DESCA, 
     $           DPDUM, -1, LINFO )
            LW3 = MAX( INT( DPDUM( 1 ) ), LW3 )
            CALL PDORGQR( MAXMN, MAXMN, MAXMN, D, ID, JD, DESCD, 
     $           DWORK, DPDUM, -1, LINFO )
            LW3 = MAX( INT( DPDUM( 1 ) ), LW3 )
         END IF
C 
C     Now compute the workspace needed in PTRGCSYD using a workspace
C     query as above
C
         IF( SOLVE ) THEN
            CALL PTRGCSYD( TRANZ, TRANAD, TRANBE, ISGN, COMM, M, N, A, 
     $           IA, JA, DESCA, B, IB, JB, DESCB, C, IC, JC, 
     $           DESCC, D, ID, JD, DESCD, E, IE, JE, DESCE,
     $           F, IF, JF, DESCF, MBNB2, DPDUM, -1, IDUM1, LIWORK, 
     $           NOEXSY, SCALE, LINFO )
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
C     depending on the given forms of (A,D) and (B,E)
C     
         WRK = LSWORK + 3 * MAXMN
C     
         IF( SOLVE ) THEN
            IF( .NOT. SCHRAD )  WRK = WRK + 2*SIZEAD
            IF( .NOT. SCHRBE )  WRK = WRK + 2*SIZEBE
            IF( .NOT. SCHRAD .OR. .NOT. SCHRBE ) THEN
               WRK = WRK + SIZECF 
            END IF
         END IF
C     
C     Check the integer workspace
C     
         DBAD = ICEIL( M, DESCA(MB_) )
         DBBE = ICEIL( N, DESCB(MB_) )
         IWRK = MAX(IWRK, DBAD + DBBE + DESCA(MB_) + DESCB(MB_) + 8) 
C
C     Now check if the call to this routine was a workspace query
C     and check if the workspace supplied is enough. 
C     
         IF( WRK.GT.LDWORK .AND. .NOT.LQUERY ) THEN
            INFO = -37
         ELSEIF( IWRK.GT.LIWORK .AND. .NOT. LQUERY ) THEN 
            INFO = -39
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
         CALL PCHK1MAT( M, 9, M, 9, IA, JA, DESCA, 14, 0, IDUM1, 
     $                  IDUM2, INFO )
         IF( INFO.EQ.0 )
     $        CALL PCHK1MAT( M, 9, M, 9, ID, JD, DESCD, 26, 0, IDUM1, 
     $                       IDUM2, INFO )   
      END IF
      IF( INFO.EQ.0 ) THEN
         CALL PCHK1MAT( N, 10, N, 10, IB, JB, DESCB, 18, 0, IDUM1, 
     $                  IDUM2, INFO )
         IF( INFO.EQ.0 )
     $        CALL PCHK1MAT( N, 10, N, 10, IE, JE, DESCE, 30, 0, IDUM1, 
     $                       IDUM2, INFO ) 
      END IF
      IF( INFO.EQ.0 .AND. SOLVE ) THEN
         CALL PCHK1MAT( M, 9, N, 10, IC, JC, DESCC, 22, 0, IDUM1,
     $                  IDUM2, INFO )
         IF( INFO.EQ.0 )
     $        CALL PCHK1MAT( M, 9, N, 10, IF, JF, DESCF, 34, 0, IDUM1,
     $                       IDUM2, INFO )  
      END IF
C     
C     Checking if we may continue or should interrupt
C     
      IF( INFO.NE.0 ) THEN
         CALL PXERBLA( ICTXT, 'PGEGCSYD', -INFO )
         RETURN
      END IF
C
C     Quick return if possible
C 
      IF( M.EQ.0 .OR. N.EQ.0 ) RETURN
C     
C     Divide the workspace between the different local pointers
C     which are needed - that depends on the given forms of A and B.
C
      IF( SOLVE ) THEN
         IF( SCHRAD .AND. SCHRBE ) THEN
            SWORK = 1
         ELSEIF( SCHRAD .AND. (.NOT. SCHRBE) ) THEN
            Q2 = 1
            Z2 = Q2 + SIZEBE
            CFCOPY = Z2 + SIZEBE
            EIGR = CFCOPY + SIZECF
            EIGI = EIGR + MAXMN
            BETA = EIGI + MAXMN
            TAUE = BETA + MAXMN
            SWORK = TAUE + NCOLSBE
         ELSEIF( (.NOT. SCHRAD) .AND. SCHRBE ) THEN
            Q1 = 1
            Z1 = Q1 + SIZEAD
            CFCOPY = Z1 + SIZEAD
            EIGR = CFCOPY + SIZECF
            EIGI = EIGR + MAXMN
            BETA = EIGI + MAXMN
            TAUD = BETA + MAXMN
            SWORK = TAUD + NCOLSAD
         ELSEIF( (.NOT. SCHRAD) .AND. (.NOT. SCHRBE) ) THEN 
            Q1 = 1
            Z1 = Q1 + SIZEAD 
            Q2 = Z1 + SIZEAD
            Z2 = Q2 + SIZEBE
            CFCOPY = Z2 + SIZEBE
            EIGR = CFCOPY + SIZECF
            EIGI = EIGR + MAXMN
            BETA = EIGI + MAXMN
            TAUD = BETA + MAXMN
            TAUE = TAUD + NCOLSAD
            SWORK = TAUE + NCOLSBE
         END IF
      ELSE
         IF( SCHRAD .AND. SCHRBE ) THEN
            SWORK = 1
         ELSEIF( SCHRAD .AND. (.NOT. SCHRBE) ) THEN
            EIGR = 1
            EIGI = EIGR + MAXMN
            BETA = EIGI + MAXMN
            TAUE = BETA + MAXMN
            SWORK = TAUE + NCOLSBE
         ELSEIF( (.NOT. SCHRAD) .AND. SCHRBE ) THEN
            EIGR = 1
            EIGI = EIGR + MAXMN
            BETA = EIGI + MAXMN
            TAUD = BETA + MAXMN
            SWORK = TAUD + NCOLSAD
         ELSEIF( (.NOT. SCHRAD) .AND. (.NOT. SCHRBE) ) THEN 
            EIGR = 1
            EIGI = EIGR + MAXMN
            BETA = EIGI + MAXMN
            TAUD = BETA + MAXMN
            TAUE = TAUD + NCOLSAD
            SWORK = TAUE + NCOLSBE
         END IF
      END IF 
      LSWORK = LDWORK-SWORK+1
C
C     Start the solution process - turn (A,D) and (B,E) into
C     Hessenberg-triangular form
C
C     QR-factorize D and E
C
      T1 = MPI_WTIME()
      IF( .NOT. SCHRAD ) 
     $     CALL PDGEQRF( M, M, D, ID, JD, DESCD, DWORK(TAUD), 
     $                   DWORK(SWORK), LSWORK, LINFO )
      IF( .NOT. SCHRBE )
     $     CALL PDGEQRF( N, N, E, IE, JE, DESCE, DWORK(TAUE), 
     $                   DWORK(SWORK), LSWORK, LINFO )
C
C     Update matrices A and B with respect to the QR-factorizations
C
      IF( .NOT. SCHRAD )
     $     CALL PDORMQR( 'Left', 'Transpose', M, M, M, D, ID, JD, DESCD, 
     $                   DWORK(TAUD), A, IA, JA, DESCA, DWORK(SWORK),
     $                   LSWORK, LINFO )
      IF( .NOT. SCHRBE )
     $     CALL PDORMQR( 'Left', 'Transpose', N, N, N, E, IE, JE, DESCE, 
     $                   DWORK(TAUE), B, IB, JB, DESCB, DWORK(SWORK),
     $                   LSWORK, LINFO )
C
C     Extract matrices Q1 and Q2 from implicit storage in
C     the matrices D and E. Also initialize the matrice Z1 and Z2 
C     as identity matrices of order M and N, respectively. Notice:
C     since we experienced problems with PDORGQR, we do this in
C     two different fashions depending on the logical parameter
C     GCSY_PDORGQR. See also "Known issues" above.
C
      IF( GCSY_PDORGQR ) THEN
         IF( .NOT. SCHRAD .AND. SOLVE ) THEN
            CALL PDLACPY( 'All', M, M, D, ID, JD, DESCD, DWORK(Q1), 
     $                    IA, JA, DESCQZ1 )
            CALL PDORGQR( M, M, M, DWORK(Q1), IA, JA, DESCQZ1, 
     $                    DWORK(TAUD), DWORK(SWORK), LSWORK, LINFO )
            CALL PDLASET( 'All', M, M, ZERO, ONE, DWORK(Z1), IA, JA, 
     $                    DESCQZ1 )
         END IF
         IF( .NOT. SCHRBE .AND. SOLVE ) THEN
            CALL PDLACPY( 'All', N, N, E, IE, JE, DESCE, DWORK(Q2), 
     $                    IB, JB, DESCQZ2 )
            CALL PDORGQR( N, N, N, DWORK(Q2), IB, JB, DESCQZ2, 
     $                    DWORK(TAUE), DWORK(SWORK), LSWORK, LINFO )
            CALL PDLASET( 'All', N, N, ZERO, ONE, DWORK(Z2), IB, JB, 
     $                    DESCQZ2 )
         END IF
      ELSE
         IF( .NOT. SCHRAD .AND. SOLVE ) THEN
            CALL PDLASET( 'All', M, M, ZERO, ONE, DWORK(Q1), IA, JA, 
     $                    DESCQZ1 )
            CALL PDORMQR( 'Left', 'NoTranspose', M, M, M, D, ID, JD, 
     $                    DESCD, DWORK(TAUD), DWORK(Q1), IA, JA, 
     $                    DESCQZ1, DWORK(SWORK), LSWORK, LINFO ) 
            CALL PDLASET( 'All', M, M, ZERO, ONE, DWORK(Z1), IA, JA, 
     $                    DESCQZ1 )
         END IF
         IF( .NOT. SCHRBE .AND. SOLVE ) THEN
            CALL PDLASET( 'All', N, N, ZERO, ONE, DWORK(Q2), IB, JB, 
     $                    DESCQZ2 )
            CALL PDORMQR( 'Left', 'NoTranspose', N, N, N, E, IE, JE, 
     $                    DESCE, DWORK(TAUE), DWORK(Q2), IB, JB, 
     $                    DESCQZ2, DWORK(SWORK), LSWORK, LINFO ) 
            CALL PDLASET( 'All', N, N, ZERO, ONE, DWORK(Z2), IB, JB, 
     $                    DESCQZ2 )
         END IF
      END IF
C    
C     Extract upper triangular matrices from D and E
C
      IF( .NOT. SCHRAD ) THEN
         CALL PDLASET( 'Lower triangular', M-1, M-1, ZERO, ZERO, D, 
     $                 ID+1, JD, DESCD ) 
      END IF
      IF( .NOT. SCHRBE ) THEN
         CALL PDLASET( 'Lower triangular', N-1, N-1, ZERO, ZERO, E, 
     $                 IE+1, JE, DESCE ) 
      END IF
C
C     Call Hessenberg-triangular reduction routine for the 
C     full-triangular pairs (A,D) and (B,E)
C      
      IF( .NOT. SCHRAD .AND. SOLVE ) THEN   
         CALL PDGGHRD( 'Vectors', 'Vectors', M, IA, IA+M-1, A, DESCA, D, 
     $                 DESCD, DWORK(Q1), DESCQZ1, DWORK(Z1), DESCQZ1, 
     $                 DWORK(SWORK), LSWORK, LINFO )
      ELSEIF( .NOT. SCHRAD ) THEN
         CALL PDGGHRD( 'NoVectors', 'NoVectors', M, IA, IA+M-1, A, 
     $                 DESCA, D, DESCD, DWORK(Q1), DESCQZ1, DWORK(Z1), 
     $                 DESCQZ1, DWORK(SWORK), LSWORK, LINFO )
      END IF   
C
      IF( .NOT. SCHRBE .AND. SOLVE ) THEN
         CALL PDGGHRD( 'Vectors', 'Vectors', N, IB, IB+N-1, B, DESCB, E, 
     $                 DESCE, DWORK(Q2), DESCQZ2, DWORK(Z2), DESCQZ2, 
     $                 DWORK(SWORK), LSWORK, LINFO )
      ELSEIF( .NOT. SCHRBE ) THEN
         CALL PDGGHRD( 'NoVectors', 'NoVectors', N, IB, IB+N-1, B, 
     $                 DESCB, E, DESCE, DWORK(Q2), DESCQZ2, DWORK(Z1), 
     $                 DESCQZ1, DWORK(SWORK), LSWORK, LINFO )
      END IF
C     
C     Apply the QZ-algorithm to compute the generalized Schur forms
C     of (A,D) and (B,E)
C     
      IF( .NOT. SCHRAD .AND. SOLVE ) THEN
         BULGES = MAXBULG
         CALL PDHGEQZ( 'Schur', 'Vectors', 'Vectors', M, IA, IA+M-1, 
     $        A, DESCA, D, DESCD, DWORK(EIGR), DWORK(EIGI), 
     $        DWORK(BETA), DWORK(Q1), DESCQZ1, DWORK(Z1), 
     $        DESCQZ1, IA, IA+M-1, IA, IA+M-1, BULGES, 
     $        DWORK(SWORK), LSWORK, LINFO )
      ELSEIF( .NOT. SCHRAD ) THEN
         BULGES = MAXBULG
         CALL PDHGEQZ( 'Schur', 'NoVectors', 'NoVectors', M, IA, 
     $        IA+M-1, A, DESCA, D, DESCD, DWORK(EIGR), DWORK(EIGI), 
     $        DWORK(BETA), DWORK(Q1), DESCQZ1, DWORK(Z1), 
     $        DESCQZ1, IA, IA+M-1, IA, IA+M-1, BULGES, 
     $        DWORK(SWORK), LSWORK, LINFO )
      END IF
C
      IF( .NOT. SCHRBE .AND. SOLVE ) THEN
         BULGES = MAXBULG
         CALL PDHGEQZ( 'Schur', 'Vectors', 'Vectors', N, IB, IB+N-1, 
     $        B, DESCB, E, DESCE, DWORK(EIGR), DWORK(EIGI), 
     $        DWORK(BETA), DWORK(Q2), DESCQZ2, DWORK(Z2), 
     $        DESCQZ2, IB, IB+N-1, IB, IB+N-1, BULGES, 
     $        DWORK(SWORK), LSWORK, LINFO )
      ELSEIF( .NOT. SCHRBE ) THEN
         BULGES = MAXBULG
         CALL PDHGEQZ( 'Schur', 'NoVectors', 'NoVectors', N, IB, 
     $        IB+N-1, B, DESCB, E, DESCE, DWORK(EIGR), DWORK(EIGI), 
     $        DWORK(BETA), DWORK(Q2), DESCQZ2, DWORK(Z2), 
     $        DESCQZ2, IB, IB+N-1, IB, IB+N-1, BULGES, 
     $        DWORK(SWORK), LSWORK, LINFO )
      END IF
C
      T1 = MPI_WTIME() - T1
      IF( .NOT. SOLVE ) RETURN
C     
C     Update (C,F) with respect to the Schur decompositions
C     
      T2 = MPI_WTIME()
      IF( .NOT. TRANSZ .AND. .NOT. SCHRAD ) THEN
         IF( .NOT. TRANSAD ) THEN
            CALL PDLACPY( 'All', M, N, C, IC, JC, DESCC, 
     $                    DWORK( CFCOPY ), 1 + IROFFCF, 1 + ICOFFCF, 
     $                    DESCCF )
            CALL PDGEMM( 'T', 'N', M, N, M, ONE, DWORK( Q1 ), IA, JA, 
     $                   DESCQZ1, DWORK( CFCOPY ), 1 + IROFFCF, 
     $                   1 + ICOFFCF, DESCCF, ZERO, C, IC, JC, DESCC )
C
            CALL PDLACPY( 'All', M, N, F, IF, JF, DESCF, 
     $                    DWORK( CFCOPY ), 1 + IROFFCF, 1 + ICOFFCF, 
     $                    DESCCF )
            CALL PDGEMM( 'T', 'N', M, N, M, ONE, DWORK( Q1 ), IA, JA, 
     $                   DESCQZ1, DWORK( CFCOPY ), 1 + IROFFCF, 
     $                   1 + ICOFFCF, DESCCF, ZERO, F, IF, JF, DESCF )
         ELSE
            CALL PDLACPY( 'All', M, N, C, IC, JC, DESCC, 
     $                    DWORK( CFCOPY ), 1 + IROFFCF, 1 + ICOFFCF, 
     $                    DESCCF )
            CALL PDGEMM( 'T', 'N', M, N, M, ONE, DWORK( Z1 ), IA, JA, 
     $                   DESCQZ1, DWORK( CFCOPY ), 1 + IROFFCF, 
     $                   1 + ICOFFCF, DESCCF, ZERO, C, IC, JC, DESCC )
C
            CALL PDLACPY( 'All', M, N, F, IF, JF, DESCF, 
     $                    DWORK( CFCOPY ), 1 + IROFFCF, 1 + ICOFFCF, 
     $                    DESCCF )
            CALL PDGEMM( 'T', 'N', M, N, M, ONE, DWORK( Z1 ), IA, JA, 
     $                   DESCQZ1, DWORK( CFCOPY ), 1 + IROFFCF, 
     $                   1 + ICOFFCF, DESCCF, ZERO, F, IF, JF, DESCF )
         END IF
      END IF
C     
      IF( .NOT. TRANSZ .AND. .NOT. SCHRBE ) THEN
         IF( .NOT. TRANSBE ) THEN
            CALL PDLACPY( 'All', M, N, C, IC, JC, DESCC, 
     $                    DWORK( CFCOPY ), 1 + IROFFCF, 1 + ICOFFCF, 
     $                    DESCCF )
            CALL PDGEMM( 'N', 'N', M, N, N, ONE, DWORK( CFCOPY ), 
     $                   1 + IROFFCF, 1 + ICOFFCF, DESCCF, DWORK( Z2 ), 
     $                   IB, JB, DESCQZ2, ZERO, C, IC, JC, DESCC )
C     
            CALL PDLACPY( 'All', M, N, F, IF, JF, DESCF, 
     $                    DWORK( CFCOPY ), 1 + IROFFCF, 1 + ICOFFCF, 
     $                    DESCCF )
            CALL PDGEMM( 'N', 'N', M, N, N, ONE, DWORK( CFCOPY ), 
     $                   1 + IROFFCF, 1 + ICOFFCF, DESCCF, DWORK( Z2 ), 
     $                   IB, JB, DESCQZ2, ZERO, F, IF, JF, DESCF )
         ELSE
            CALL PDLACPY( 'All', M, N, C, IC, JC, DESCC, 
     $                    DWORK( CFCOPY ), 1 + IROFFCF, 1 + ICOFFCF, 
     $                    DESCCF )
            CALL PDGEMM( 'N', 'N', M, N, N, ONE, DWORK( CFCOPY ), 
     $                   1 + IROFFCF, 1 + ICOFFCF, DESCCF, DWORK( Q2 ), 
     $                   IB, JB, DESCQZ2, ZERO, C, IC, JC, DESCC )
C     
            CALL PDLACPY( 'All', M, N, F, IF, JF, DESCF, 
     $                    DWORK( CFCOPY ), 1 + IROFFCF, 1 + ICOFFCF, 
     $                    DESCCF )
            CALL PDGEMM( 'N', 'N', M, N, N, ONE, DWORK( CFCOPY ), 
     $                   1 + IROFFCF, 1 + ICOFFCF, DESCCF, DWORK( Q2 ), 
     $                   IB, JB, DESCQZ2, ZERO, F, IF, JF, DESCF )
         END IF
      END IF 
C
      IF( TRANSZ .AND. .NOT. SCHRAD ) THEN
         IF( .NOT. TRANSAD ) THEN
            CALL PDLACPY( 'All', M, N, C, IC, JC, DESCC, 
     $                    DWORK( CFCOPY ), 1 + IROFFCF, 1 + ICOFFCF, 
     $                    DESCCF )
            CALL PDGEMM( 'T', 'N', M, N, M, ONE, DWORK( Z1 ), IA, JA, 
     $                   DESCQZ1, DWORK( CFCOPY ), 1 + IROFFCF, 
     $                   1 + ICOFFCF, DESCCF, ZERO, C, IC, JC, DESCC )
C         
            CALL PDLACPY( 'All', M, N, F, IF, JF, DESCF, 
     $                    DWORK( CFCOPY ), 1 + IROFFCF, 1 + ICOFFCF, 
     $                    DESCCF )
            CALL PDGEMM( 'T', 'N', M, N, M, ONE, DWORK( Q1 ), IA, JA, 
     $                   DESCQZ1, DWORK( CFCOPY ), 1 + IROFFCF, 
     $                   1 + ICOFFCF, DESCCF, ZERO, F, IF, JF, DESCF )
         ELSE
            CALL PDLACPY( 'All', M, N, C, IC, JC, DESCC, 
     $                    DWORK( CFCOPY ), 1 + IROFFCF, 1 + ICOFFCF, 
     $                    DESCC )
            CALL PDGEMM( 'T', 'N', M, N, M, ONE, DWORK( Q1 ), IA, JA, 
     $                   DESCQZ1, DWORK( CFCOPY ), 1 + IROFFCF, 
     $                   1 + ICOFFCF, DESCCF, ZERO, C, IC, JC, DESCC )
C         
            CALL PDLACPY( 'All', M, N, F, IF, JF, DESCF, 
     $                    DWORK( CFCOPY ), 1 + IROFFCF, ICOFFCF, DESCF )
            CALL PDGEMM( 'T', 'N', M, N, M, ONE, DWORK( Z1 ), IA, JA, 
     $                   DESCQZ1, DWORK( CFCOPY ), 1 + IROFFCF, 
     $                   1 + ICOFFCF, DESCCF, ZERO, F, IF, JF, DESCF )
         END IF 
      END IF
C
      IF( TRANSZ .AND. .NOT. SCHRBE ) THEN
         IF( .NOT. TRANSBE ) THEN
            CALL PDLACPY( 'All', M, N, C, IC, JC, DESCC, 
     $                    DWORK( CFCOPY ), 1 + IROFFCF, 1 + ICOFFCF, 
     $                    DESCCF )
            CALL PDGEMM( 'N', 'N', M, N, N, ONE, DWORK( CFCOPY ), 
     $                   1 + IROFFCF, 1 + ICOFFCF, DESCCF, DWORK( Z2 ), 
     $                   IB, JB, DESCQZ2, ZERO, C, IC, JC, DESCC )
C
            CALL PDLACPY( 'All', M, N, F, IF, JF, DESCF, 
     $                    DWORK( CFCOPY ), 1 + IROFFCF, 1 + ICOFFCF, 
     $                    DESCCF )            
            CALL PDGEMM( 'N', 'N', M, N, N, ONE, DWORK( CFCOPY ), 
     $                   1 + IROFFCF, 1 + ICOFFCF, DESCCF, DWORK( Q2 ), 
     $                   IB, JB, DESCQZ2, ZERO, F, IF, JF, DESCF )
         ELSE
            CALL PDLACPY( 'All', M, N, C, IC, JC, DESCC, 
     $                    DWORK( CFCOPY ), 1 + IROFFCF, 1 + ICOFFCF, 
     $                    DESCCF ) 
            CALL PDGEMM( 'N', 'N', M, N, N, ONE, DWORK( CFCOPY ),
     $                   1 + IROFFCF, 1 + ICOFFCF, DESCCF, DWORK( Q2 ), 
     $                   IB, JB, DESCQZ2, ZERO, C, IC, JC, DESCC )
C
            CALL PDLACPY( 'All', M, N, F, IF, JF, DESCF, 
     $                    DWORK( CFCOPY ), 1 + IROFFCF, 1 + ICOFFCF, 
     $                    DESCF )
            CALL PDGEMM( 'N', 'N', M, N, N, ONE, DWORK( CFCOPY ),
     $                   1 + IROFFCF, 1 + ICOFFCF, DESCCF, DWORK( Z2 ), 
     $                   IB, JB, DESCQZ2, ZERO, F, IF, JF, DESCCF )
         END IF
      END IF
      T2 = MPI_WTIME() - T2
C
C     Now we have reduced our general equation to a (quasi-)triangular
C     case. Solve now the reduced problem with a call to a triangular 
C     solver - the solution to that reduced problem is stored
C     on (C,F). 
C     
      T3 = MPI_WTIME()
      CALL PTRGCSYD( TRANZ, TRANAD, TRANBE, ISGN, COMM, M, N, A, IA,
     $               JA, DESCA, B, IB, JB, DESCB, C, IC, JC, DESCC, D, 
     $               ID, JD, DESCD, E, IE, JE, DESCE, F, IF, JF, DESCF, 
     $               MBNB2, DWORK( SWORK ), LSWORK, IWORK, LIWORK, 
     $               NOEXSY, SCALE, LINFO )

      IF( INFO.NE.0 ) THEN
         IF( LINFO.EQ.1 ) INFO = 2
         IF( LINFO.EQ.2 ) INFO = 3
      END IF
      T3 = MPI_WTIME() - T3
C     
C     Transform the solution back to the original coordinate system
C
      T22 = MPI_WTIME()
      IF( .NOT. TRANSZ .AND. .NOT. SCHRAD ) THEN
         IF( .NOT. TRANSAD ) THEN
            CALL PDLACPY( 'All', M, N, C, IC, JC, DESCC, 
     $                    DWORK( CFCOPY ), 1 + IROFFCF, 1 + ICOFFCF, 
     $                    DESCCF )
            CALL PDGEMM( 'N', 'N', M, N, M, ONE, DWORK( Z1 ), IA, JA, 
     $                   DESCQZ1, DWORK( CFCOPY ), 1 + IROFFCF, 
     $                   1 + ICOFFCF, DESCCF, ZERO, C, IC, JC, DESCC )
C         
            CALL PDLACPY( 'All', M, N, F, IF, JF, DESCF, 
     $                    DWORK( CFCOPY ), 1 + IROFFCF, 1 + ICOFFCF, 
     $                    DESCCF )
            CALL PDGEMM( 'N', 'N', M, N, M, ONE, DWORK( Q1 ), IA, JA, 
     $                   DESCQZ1, DWORK( CFCOPY ), 1 + IROFFCF, 
     $                   1 + ICOFFCF, DESCCF, ZERO, F, IF, JF, DESCF )
         ELSE
            CALL PDLACPY( 'All', M, N, C, IC, JC, DESCC, 
     $                    DWORK( CFCOPY ), 1 + IROFFCF, 1 + ICOFFCF, 
     $                    DESCC )
            CALL PDGEMM( 'N', 'N', M, N, M, ONE, DWORK( Q1 ), IA, JA, 
     $                   DESCQZ1, DWORK( CFCOPY ), 1 + IROFFCF, 
     $                   1 + ICOFFCF, DESCCF, ZERO, C, IC, JC, DESCC )
C         
            CALL PDLACPY( 'All', M, N, F, IF, JF, DESCF, 
     $                    DWORK( CFCOPY ), 1 + IROFFCF, ICOFFCF, DESCF )
            CALL PDGEMM( 'N', 'N', M, N, M, ONE, DWORK( Z1 ), IA, JA, 
     $                   DESCQZ1, DWORK( CFCOPY ), 1 + IROFFCF, 
     $                   1 + ICOFFCF, DESCCF, ZERO, F, IF, JF, DESCF )
         END IF 
      END IF
C
      IF( .NOT. TRANSZ .AND. .NOT. SCHRBE ) THEN
         IF( .NOT. TRANSBE ) THEN
            CALL PDLACPY( 'All', M, N, C, IC, JC, DESCC, 
     $                    DWORK( CFCOPY ), 1 + IROFFCF, 1 + ICOFFCF, 
     $                    DESCCF ) 
            CALL PDGEMM( 'N', 'T', M, N, N, ONE, DWORK( CFCOPY ), 
     $                   1 + IROFFCF, 1 + ICOFFCF, DESCCF, DWORK( Z2 ), 
     $                   IB, JB, DESCQZ2, ZERO, C, IC, JC, DESCC )
C
            CALL PDLACPY( 'All', M, N, F, IF, JF, DESCF, 
     $                    DWORK( CFCOPY ), 1 + IROFFCF, 1 + ICOFFCF, 
     $                    DESCCF )
            CALL PDGEMM( 'N', 'T', M, N, N, ONE, DWORK( CFCOPY ), 
     $                   1 + IROFFCF, 1 + ICOFFCF, DESCCF, DWORK( Q2 ), 
     $                   IB, JB, DESCQZ2, ZERO, F, IF, JF, DESCF )
         ELSE
            CALL PDLACPY( 'All', M, N, C, IC, JC, DESCC, 
     $                    DWORK( CFCOPY ), 1 + IROFFCF, 1 + ICOFFCF, 
     $                    DESCCF ) 
            CALL PDGEMM( 'N', 'T', M, N, N, ONE, DWORK( CFCOPY ),
     $                   1 + IROFFCF, 1 + ICOFFCF, DESCCF, DWORK( Q2 ), 
     $                   IB, JB, DESCQZ2, ZERO, C, IC, JC, DESCC )
C
            CALL PDLACPY( 'All', M, N, F, IF, JF, DESCF, 
     $                    DWORK( CFCOPY ), 1 + IROFFCF, 1 + ICOFFCF, 
     $                    DESCF )
            CALL PDGEMM( 'N', 'T', M, N, N, ONE, DWORK( CFCOPY ),
     $                   1 + IROFFCF, 1 + ICOFFCF, DESCCF, DWORK( Z2 ), 
     $                   IB, JB, DESCQZ2, ZERO, F, IF, JF, DESCCF )
         END IF
      END IF
C
      IF( TRANSZ .AND. .NOT. SCHRAD ) THEN
         IF( .NOT. TRANSAD ) THEN
            CALL PDLACPY( 'All', M, N, C, IC, JC, DESCC, 
     $                    DWORK( CFCOPY ), 1 + IROFFCF, 1 + ICOFFCF, 
     $                    DESCCF )
            CALL PDGEMM( 'N', 'N', M, N, M, ONE, DWORK( Q1 ), IA, JA, 
     $                   DESCQZ1, DWORK( CFCOPY ), 1 + IROFFCF, 
     $                   1 + ICOFFCF, DESCCF, ZERO, C, IC, JC, DESCC )
C
            CALL PDLACPY( 'All', M, N, F, IF, JF, DESCF, 
     $                    DWORK( CFCOPY ), 1 + IROFFCF, 1 + ICOFFCF, 
     $                    DESCCF )
            CALL PDGEMM( 'N', 'N', M, N, M, ONE, DWORK( Q1 ), IA, JA, 
     $                   DESCQZ1, DWORK( CFCOPY ), 1 + IROFFCF, 
     $                   1 + ICOFFCF, DESCCF, ZERO, F, IF, JF, DESCCF )
         ELSE
            CALL PDLACPY( 'All', M, N, C, IC, JC, DESCC, 
     $                    DWORK( CFCOPY ), 1 + IROFFCF, 1 + ICOFFCF, 
     $                    DESCCF )
            CALL PDGEMM( 'N', 'N', M, N, M, ONE, DWORK( Z1 ), IA, JA, 
     $                   DESCQZ1, DWORK( CFCOPY ), 1 + IROFFCF, 
     $                   1 + ICOFFCF, DESCCF, ZERO, C, IC, JC, DESCC )
C
            CALL PDLACPY( 'All', M, N, F, IF, JF, DESCF, 
     $                    DWORK( CFCOPY ), 1 + IROFFCF, 1 + ICOFFCF, 
     $                    DESCCF )
            CALL PDGEMM( 'N', 'N', M, N, M, ONE, DWORK( Z1 ), IA, JA, 
     $                   DESCQZ1, DWORK( CFCOPY ), 1 + IROFFCF, 
     $                   1 + ICOFFCF, DESCCF, ZERO, F, IF, JF, DESCF )
         END IF
      END IF
C     
      IF( TRANSZ .AND. .NOT. SCHRBE ) THEN
         IF( .NOT. TRANSBE ) THEN
            CALL PDLACPY( 'All', M, N, C, IC, JC, DESCC, 
     $                    DWORK( CFCOPY ), 1 + IROFFCF, 1 + ICOFFCF, 
     $                    DESCCF )
            CALL PDGEMM( 'N', 'T', M, N, N, ONE, DWORK( CFCOPY ), 
     $                   1 + IROFFCF, 1 + ICOFFCF, DESCCF, DWORK( Z2 ), 
     $                   IB, JB, DESCQZ2, ZERO, C, IC, JC, DESCC )
C     
            CALL PDLACPY( 'All', M, N, F, IF, JF, DESCF, 
     $                    DWORK( CFCOPY ), 1 + IROFFCF, 1 + ICOFFCF, 
     $                    DESCF )
            CALL PDGEMM( 'N', 'T', M, N, N, ONE, DWORK( CFCOPY ), 
     $                   1 + IROFFCF, 1 + ICOFFCF, DESCCF, DWORK( Z2 ), 
     $                   IB, JB, DESCQZ2, ZERO, F, IF, JF, DESCF )
         ELSE
            CALL PDLACPY( 'All', M, N, C, IC, JC, DESCC, 
     $                    DWORK( CFCOPY ), 1 + IROFFCF, 1 + ICOFFCF, 
     $                    DESCCF )
            CALL PDGEMM( 'N', 'T', M, N, N, ONE, DWORK( CFCOPY ), 
     $                   1 + IROFFCF, 1 + ICOFFCF, DESCCF, DWORK( Q2 ), 
     $                   IB, JB, DESCQZ2, ZERO, C, IC, JC, DESCC )
C     
            CALL PDLACPY( 'All', M, N, F, IF, JF, DESCF, 
     $                    DWORK( CFCOPY ), 1 + IROFFCF, 1 + ICOFFCF, 
     $                    DESCCF )
            CALL PDGEMM( 'N', 'T', M, N, N, ONE, DWORK( CFCOPY ), 
     $                   1 + IROFFCF, 1 + ICOFFCF, DESCCF, DWORK( Q2 ), 
     $                   IB, JB, DESCQZ2, ZERO, F, IF, JF, DESCF )
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
C     END OF PGEGCSYD 
C
C *** Last line of PGEGCSYD ***
