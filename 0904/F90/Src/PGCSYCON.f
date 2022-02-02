CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C         This file is a part of the SCASY software library.           C
C         See http://www.cs.umu.se/parallel/research/scasy             C
C         Contributors: Robert Granat and Bo Kågström.                 C
C         SCASY Version 0.10,    March 31, 2006.                       C
C         Copyright 2006, Umeå University, Sweden.                     C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE PGCSYCON( TRANAD, TRANBE, ISGN, COMM, M, N, A, IA, 
     $                     JA, DESCA, B, IB, JB, DESCB, D, ID, JD, 
     $                     DESCD, E, IE, JE, DESCE, MBNB2, DWORK, 
     $                     LDWORK, IWORK, LIWORK, EST, NOITER, INFO )
C
C  -- ScaLAPACK-style routine (preliminary version)--
C     HPC2N and Department of Computing Science,
C     University of Umeå, Sweden.
C     Written by Robert Granat, (granat@cs.umu.se)
C     January 21, 2006.
C
      IMPLICIT NONE
C
C     .. Scalar Arguments ..
      CHARACTER*1          TRANAD, TRANBE, COMM
      INTEGER              ISGN, M, N, IA, JA, IB, JB, ID, JD, IE, JE,
     $                     LDWORK, LIWORK, 
     $                     NOITER, INFO
      DOUBLE PRECISION     EST
C     ..
C     .. Array Arguments ..
      INTEGER              DESCA( * ), DESCB( * ),  DESCD( * ), 
     $                     DESCE( * ), IWORK( * ), MBNB2( 2 )
      DOUBLE PRECISION     A( * ), B( * ), D( * ), E( * ), DWORK( * )
C     ..
C 
C  Purpose
C  =======
C  This subroutine gives a 1-norm based bound of the spectral norm of the
C  2*M*N x 2*M*N matrix Z^(-1) corresponding to the triangular GCSY 
c  equation 
C
C  (op(A) * X +/- Y * op(B),op(D) * X +/- Y * op(E)) = (C,F)    
C
C  where A and D are M x M, B and E are N x N and C and F are M x N.     
C  The matrices A, B have to be generalized Schur form.
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
C  TRANAD    (global input) CHARACTER*1
C            If TRANAD = 'N', then op(A) = A and op(D) = D.
C            If TRANAD = 'T', then op(A) = A**T and op(D) = D**T.
C            
C  TRANBE    (global input) CHARACTER*1
C            If TRANBE = 'N', then op(B) = B and op(E) = E.
C            If TRANBE = 'T', then op(B) = B**T and op(E) = E**T. 
C
C  ISGN      (global input) INTEGER*1
C            This integer determines the sign in the equation.
C
C  Input/Output parameters
C
C  COMM      (global input/output) CHARACTER*1
C            This subroutine uses two different communications schemes in
C            solving the reduced triangular problem, as follows:
C              If COMM = 'S', the "shifts" scheme is used.
C              If COMM = 'D', the "on demand" scheme is used.
C            The choice COMM = 'S' is only valid for TRANAD = TRANBE = 'N' 
C            or TRANAD = TRANBE = 'T'. The scheme used will be output.
C            See the references for details.
C
C  M         (global input) INTEGER
C            The number of rows and columns of the global distributed
C            matrix A. M >= 0.  
C
C  N         (global input) INTEGER
C            The number of rows and columns of the global distributed
C            matrix B. N >= 0.
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
C  D         (local input) DOUBLE PRECISION array 
C            Array of dimension (LLD_D,LOCc(M)). Contains the local
C            pieces of the global distributed matrix D.
C
C  ID        (global input) INTEGER
C            Row start index for sub(D), i.e., the submatrix to operate
C            on. ID >= 1.
C
C  JD        (global input) INTEGER
C            Column start index for sub(D), i.e., the submatrix to 
C            operate on. JD = ID must hold.
C
C  DESCD     (global and local input) INTEGER array of dimension DLEN_.
C            The array descriptor for the global distributed matrix D.
C
C  E         (local input) DOUBLE PRECISION array 
C            Array of dimension (LLD_E,LOCc(N)). Contains the local
C            pieces of the global distributed matrix E.
C
C  IE        (global input) INTEGER
C            Row start index for sub(E), i.e., the submatrix to operate
C            on. IE >= 1.
C
C  JE        (global input) INTEGER
C            Column start index for sub(E), i.e., the submatrix to 
C            operate on. JE = IE must hold.
C
C  DESCE     (global and local input) INTEGER array of dimension DLEN_.
C            The array descriptor for the global distributed matrix E.
C
C  MBNB2     (global input) INTEGER array of dimension 2.
C            Internal blocking factors for pipelining of subsolutions
C            for updates of the matrix C in PTRGCSYD (see the references 
C            for details).
C            1 < = MBNB2(1) <= DESCA( MB_ ) and 
C            1 < = MBNB2(2) <= DESCB( NB_ ) must hold.
C
C  Workspaces
C
C  DWORK     (local workspace) DOUBLE PRECISION array, dimension
C            LDWORK. This array is supposed to contain all the extra
C            arrays and workspaces we need to use throughout the 
C            solution process.
C
C  LDWORK    (local or global input) INTEGER
C            The dimension of the array DWORK. 
C            LDWORK >= BUFF + 2*LOCr(M*N+IROFFAD) + LOCr(M+IROFFAD) *
C                      LOCc(N+IROFFBE), 
C            where BUFF is the workspace needed by PTRGCSYD (i.e., the 
C            same as for PGEGCSYD with ASCHUR = BSCHUR = 'S').
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
C            LIWORK >= MAX( 2*BUFF, DBA + DBB ) + LOCr(M*N+IROFFAD), 
C            where DBA = ICEIL(LOCr(IA+IROFFAD),MB_A) and
C            DBB = ICEIL(LOCr(IB+IROFFBE),MB_B). For BUFF, see LDWORK.
C            
C            If LIWORK = -1, LIWORK is global input and a workspace 
C            query is assumed. The routine will then calculate the 
C            optimal workspace needed, store it in IWORK(1) and return 
C            immediately. No error will then be signaled by PXERBLA.
C
C Output information
C
C  EST       (global output) DOUBLE PRECISION
C            The resulting 1-norm based estimation of the spectral norm
C            of Z^(-1).
C
C  NOITER    (global output) INTEGER
C            The number of iterations used to compute the result. NOITER
C            corresponds to the number of call to PDLACON and PGEGCSYD.
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
C             (but the matrices A and B is unchanged).
C
C  Method
C  ======
C  See the references.
C
C  Additional requirements
C  =======================
C
C  (A,D) must be distributed using the same blocking factor in 
C  each direction, i.e., MB_A=NB_A=MB_D=NB_D. (B,E) must be distributed 
C  using the same blocking factor in each direction, i.e., MB_B=NB_B=
C  =MB_E=NB_E. The matrix pairs (A,D) and (B,E) must be internally 
C  aligned.
C
C  In this version, the bound can only be computed for TRANAD = TRANBE.
C
C  References
C  ==========
C
C  [1] KÅGSTRÖM, B., POROMAA, P., (1992). Distributed and Shared Memory 
C      Block Algorithms for the TriangularSylvester Equation with sep^{-1} 
C      Estimators, SIAM J. Matrix Anal. Appl., VOL. 13, NO 1, pp. 90-101,
C      January 1992.
C
C  [2] Robert Granat and Bo Kågström, Parallel ScaLAPACK-style Algorithms 
C      for Standard and Generalized Sylvester-Type Matrix Equations, in 
C      preparation, Department of Computing Science and HPC2N, Umeå 
C      University, 2006.
C 
C  [3] Robert Granat and Bo Kågström, SCASY - A ScaLAPACK-style High 
C      Performance Library for Sylvester-Type Matrix Equations, in 
C      preparation, Department of Computing Science and HPC2N, Umeå 
C      University, 2006. 
C
C  See also: http://www.cs.umu.se/research/parallel/scasy
C
C  Parallel execution recommendations
C  ==================================
C  None.
C
C  Revisions
C  =========
C  Please report bugs to <granat@cs.umu.se>
C
C  Keywords
C  ========
C  Z-GCSY, 1-norm based lower bound, continuous-time Sylvester equation, 
C  sep^-1, condition number, condition estimation.
C
C  =====================================================================
C
C     .. Parameters ..
      INTEGER            BLOCK_CYCLIC_2D, CSRC_, CTXT_, DLEN_, DTYPE_,
     $                   LLD_, MB_, M_, NB_, N_, RSRC_
      DOUBLE PRECISION   ZERO
      PARAMETER          ( BLOCK_CYCLIC_2D = 1, DLEN_ = 9, DTYPE_ = 1,
     $                     CTXT_ = 2, M_ = 3, N_ = 4, MB_ = 5, NB_ = 6,
     $                     RSRC_ = 7, CSRC_ = 8, LLD_ = 9, 
     $                     ZERO = 0.0D+0 )
C     .. Local Scalars ..
      LOGICAL            LQUERY, FIRST, SHIFT, TRANSAD, TRANSBE
      CHARACTER*1        TRANZ1, TRANZ2
      INTEGER            ICTXT, NPROW, NPCOL, MYROW, MYCOL, KASE, ITER,
     $                   NOEXSY, WRK, IWRK, CONWRK, WRKNEED, SIZE_CF,
     $                   CFCOLS, IPW, IPV, IPC, IPX, XCSRC, MB, NB,
     $                   LLDWORK, LLDV, MN, NPROCS, IPX2, IPIW1, IPIW2,
     $                   ICF, JCF, IROFFAD, ICOFFAD, IROFFBE, ICOFFBE, 
     $                   LIA, LJA, ADRSRC, ADCSRC, LIB, LJB, BERSRC, 
     $                   BECSRC, ADROWS, ADCOLS, CFROWS, BEROWS, VECDIM, 
     $                   BECOLS, EXMEMAD, EXMEMBE, EXMEMCF, DBAD, DBBE, 
     $                   BUFF, DIMVEC, IPF, MWORK, XVCDIM, CNWRK2, IPY,
     $                   EAST, WEST, INDXX, XSIZE, ROUNDS, XPCOL, XCOLS,
     $                   INDXY, YSIZE, YPCOL, YCOLS, GKASE, GKASE_RST,
     $                   NPCOL0
      DOUBLE PRECISION   SCALE, NORM, T1, T2, T3, T
C     ..
C     .. Local Arrays ..
      INTEGER            DESCXV( DLEN_ ), DESCCF( DLEN_ )
C     ..
C     .. External Subroutines ..
      EXTERNAL           BLACS_GRIDINFO, PDLACON, PGEGCSYD, PMAT2VEC,
     $                   PVEC2MAT, CHK1MAT
C     ..
C     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            NUMROC, ICEIL
      DOUBLE PRECISION   PDLANGE, MPI_WTIME
      EXTERNAL           NUMROC, LSAME, PDLANGE, ICEIL, MPI_WTIME
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC          MAX
C     ..
C     .. Executable Statements ..
C
C     Get grid parameters
C     
      ICTXT = DESCA( CTXT_ )
      CALL BLACS_GRIDINFO( ICTXT, NPROW, NPCOL, MYROW, MYCOL )
      NPROCS = NPROW * NPCOL
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
C     Check if NPROW and NPCOL are O.K. etc.
C     
      IF( INFO.EQ.0 ) THEN
         CALL CHK1MAT( M, 5, M, 5, IA, JA, DESCA, 10, INFO )
         IF( INFO.EQ.0 ) THEN
            CALL CHK1MAT( N, 6, N, 6, IB, JB, DESCB, 14, INFO )
         END IF
         IF( INFO.EQ.0 ) THEN
            CALL CHK1MAT( M, 5, M, 5, ID, JD, DESCD, 18, INFO )
         END IF
         IF( INFO.EQ.0 ) THEN
            CALL CHK1MAT( N, 6, N, 6, IE, JE, DESCE, 22, INFO )
         END IF
      END IF
C     
C     Check the blocking sizes for equivalence
C     
      IF( INFO.EQ.0 ) THEN
         IF( DESCA( MB_ ).NE.DESCA( NB_ ) ) INFO = -(100*10 + MB_)
      END IF
      IF( INFO.EQ.0 ) THEN
         IF( DESCB( MB_ ).NE.DESCB( NB_ ) ) INFO = -(100*14 + MB_)
      END IF
      IF( INFO.EQ.0 ) THEN
         IF( DESCD( MB_ ).NE.DESCD( NB_ ) ) INFO = -(100*18 + MB_)
         IF( INFO.EQ.0 ) THEN
            IF( DESCD( MB_ ).NE.DESCA( MB_ ) ) INFO = -(100*18 + MB_)
         END IF
         IF( INFO.EQ.0 ) THEN
            IF( DESCD( NB_ ).NE.DESCA( NB_ ) ) INFO = -(100*18 + NB_)
         END IF
      END IF
      IF( INFO.EQ.0 ) THEN
         IF( DESCE( MB_ ).NE.DESCE( NB_ ) ) INFO = -(100*22 + MB_)
         IF( INFO.EQ.0 ) THEN
            IF( DESCE( MB_ ).NE.DESCB( MB_ ) ) INFO = -(100*22 + MB_)
         END IF
         IF( INFO.EQ.0 ) THEN
            IF( DESCE( NB_ ).NE.DESCB( NB_ ) ) INFO = -(100*22 + NB_)
         END IF
      END IF
C
C     Check alignment, part 1: check RSRC_ and CSRC_
C
      IF( INFO.EQ.0 ) THEN
         IF( DESCD( RSRC_ ) .NE. DESCA( RSRC_ ) ) 
     $        INFO = -(100*18 + RSRC_)
         IF( DESCD( CSRC_ ) .NE. DESCA( CSRC_ ) ) 
     $        INFO = -(100*18 + CSRC_)
         IF( DESCE( RSRC_ ) .NE. DESCB( RSRC_ ) ) 
     $        INFO = -(100*22 + RSRC_)
         IF( DESCE( CSRC_ ) .NE. DESCB( CSRC_ ) ) 
     $        INFO = -(100*22 + CSRC_)
      END IF
C      
C     Check the values of TRANAD and TRANBE
C
      IF( INFO.EQ.0 ) THEN
         IF( .NOT.( LSAME( TRANAD, 'N' ) .OR. LSAME( TRANAD, 'T' ))) 
     $        THEN
            INFO = -1
         END IF
         IF( INFO.EQ.0 ) THEN
            IF( .NOT.( LSAME( TRANBE, 'N' ) .OR. LSAME( TRANBE, 'T' ))) 
     $           THEN
               INFO = -2
            END IF  
         END IF
      END IF
C
C     Check that TRANAD and TRANBE are really the same
C
      IF( .NOT. LSAME( TRANAD, TRANBE ) ) THEN
         INFO = -1
      END IF
      TRANSAD = LSAME( TRANAD, 'T' )
      TRANSBE = LSAME( TRANBE, 'T' )
C
C     Check the value of ISGN
C
      IF( INFO.EQ.0 ) THEN
         IF( .NOT.(ISGN.EQ.1 .OR. ISGN.EQ.-1) ) THEN
            INFO = -3
         END IF
      END IF
C
C     Check the value of COMM
C
      IF( INFO.EQ.0 ) THEN
         IF( .NOT.( LSAME( COMM, 'S' ) .OR. LSAME( COMM, 'D' ))) 
     $        THEN
            INFO = -4
         END IF
      END IF
      SHIFT = LSAME( COMM, 'S' )
C
C     Check alignment, part 2: check IA, JA, IB, JB, IC, JC etc.
C
      IF( INFO.EQ.0 ) THEN
         IF( IA.LT.1 ) INFO = -8
         IF( JA.NE.IA ) INFO = -9
         IF( IB.LT.1 ) INFO = -12
         IF( JB.NE.IB ) INFO = -13
         IF( ID .NE. IA ) INFO = -16
         IF( JD .NE. JA ) INFO = -17
         IF( IE .NE. IB ) INFO = -20
         IF( JE .NE. JB ) INFO = -21
      END IF
C
C     Check if workspace query
C
      LQUERY = LDWORK.EQ.-1 .OR. LIWORK.EQ.-1
C
C     Test workspace
C
      IF( INFO.EQ.0 .OR. LQUERY ) THEN
C
C     Compute IROFFAD, IROFFBE, etc.
C
         IROFFAD = MOD( IA - 1, DESCA(MB_) )
         IROFFBE = MOD( IB - 1, DESCB(MB_) )
         CALL INFOG2L( IA, JA, DESCA, NPROW, NPCOL, MYROW, MYCOL,
     $                 LIA, LJA, ADRSRC, ADCSRC )
         CALL INFOG2L( IB, JB, DESCB, NPROW, NPCOL, MYROW, MYCOL,
     $                 LIB, LJB, BERSRC, BECSRC )
         ADROWS = NUMROC( M+IROFFAD, DESCA( MB_ ), MYROW, ADRSRC, 
     $                    NPROW )
         ADCOLS = NUMROC( M+IROFFAD, DESCA( NB_ ), MYCOL, ADCSRC, 
     $                    NPCOL )
         CFROWS = ADROWS
         BEROWS = NUMROC( N+IROFFBE, DESCB( MB_ ), MYROW, BERSRC, 
     $                    NPROW )
         BECOLS = NUMROC( N+IROFFBE, DESCB( NB_ ), MYCOL, BECSRC, 
     $                    NPCOL )
         CFCOLS = BECOLS
C
C     Set a descriptor for the global distributed matrix pair (C,F), 
C     which will form our right hand sides in the equation
C
         CALL DESCINIT( DESCCF, M+IROFFAD, N+IROFFBE, DESCA(MB_), 
     $                  DESCB(NB_), ADRSRC, BECSRC, ICTXT,  
     $                  MAX( 1, CFROWS ), INFO )
         ICF = 1 + IROFFAD
         JCF = 1 + IROFFBE
         VECDIM = 2*(M+IROFFAD)*(N+IROFFBE)
         XVCDIM = (M+IROFFAD)*(N+IROFFBE)
*         VECDIM = 2*M*N + IROFFAD
*         XVCDIM = M*N + IROFFAD
C     
C     Compute needed workspace
C     
         MB = DESCA( MB_ )
         NB = DESCB( NB_ )
         EXMEMAD = 2*ICEIL(ADROWS,MB)*ICEIL(ADCOLS,MB)*(2*MB+1)
         EXMEMBE = 2*ICEIL(BEROWS,NB)*ICEIL(BECOLS,NB)*(2*NB+1)
         EXMEMCF = 2*ICEIL(ADROWS,MB)*ICEIL(BECOLS,NB)*(MB+NB+1) 
         DBAD = ICEIL( M + IROFFAD, DESCA(MB_) )
         DBBE = ICEIL( N + IROFFBE, DESCB(MB_) )   
         SIZE_CF = DESCCF( LLD_ ) * CFCOLS
         CONWRK = NUMROC( VECDIM, DESCCF(MB_), MYROW, DESCCF(RSRC_), 
     $                    NPROW )
         CNWRK2 = NUMROC( XVCDIM, DESCCF(MB_), MYROW, DESCCF(RSRC_), 
     $                    NPROW )
C 
         IF( SHIFT ) THEN
            IF( M.EQ.N ) MWORK = MAX( ADROWS*ADCOLS,BEROWS*BECOLS )
            IF( M.GT.N .AND. NPCOL.EQ.1 ) MWORK = ADROWS*ADCOLS
            IF( M.GT.N .AND. NPCOL.GT.1 ) MWORK = 
     $           MAX( ADROWS*ADCOLS, ADROWS*BECOLS )
            IF( M.LE.N .AND. NPROW.EQ.1 ) MWORK = BEROWS*BECOLS
            IF( M.LE.N .AND. NPROW.GT.1 ) MWORK = 
     $           MAX(BEROWS*BECOLS , ADROWS*BECOLS )
            IF( SHIFT .AND. TRANSAD .AND. .NOT. TRANSBE ) 
     $           MWORK = MAX( ADROWS*ADCOLS, MWORK )
            IF( SHIFT .AND. TRANSBE .AND. .NOT. TRANSAD )
     $           MWORK = MAX( BEROWS*BECOLS, MWORK )
         ELSE
            MWORK = 0
         END IF
         BUFF = MAX( 1, EXMEMAD + EXMEMBE + EXMEMCF ) + 
     $          4*(MB+1)*(NB+1) + 2*(MB+1)**2 + 2*(NB+1)**2 + 
     $          MAX( MB, NB ) + MWORK
         WRK = BUFF + 2 * CONWRK + 2 * SIZE_CF
         IWRK = MAX( 2*BUFF, DBAD + DBBE ) + CONWRK
         IF( LQUERY ) THEN
            DWORK(1) = WRK
            IWORK(1) = IWRK
            INFO = 0
            RETURN
         ELSEIF( LDWORK.LT.WRK ) THEN
            DWORK(1) = WRK
            INFO = -25
         ELSEIF( LIWORK.LT.IWRK ) THEN
            IWORK(1) = IWRK
            INFO = -27
         END IF
      END IF
C     
C     If INFO is negative call PXERBLA
C
      IF( INFO.LT.0 ) THEN
         CALL PXERBLA( ICTXT, 'PGCSYCON', -INFO )
         RETURN
      END IF
C
C     Quick return if needed
C
      IF( INFO.GT.0 ) RETURN 
C
C     Divide the workspace among ScaLAPACK arrays
C
      IPC  = 1
      IPF  = IPC + SIZE_CF
      IPV  = IPF + SIZE_CF
      IPX  = IPV + CONWRK
      IPY  = IPX + CNWRK2
      IPW  = IPX + CONWRK
C
C     Divide the integer workspace between PDLACON and PGEGCSYD
C
      IPIW1 = 1
      IPIW2 = IPIW1 + CONWRK
C
C     Start estimation work
C
C     Set some starting values
C
      NOITER = 0
      TRANZ1 = 'N'
      TRANZ2 = 'T'
C
C     Compute communication directions for column-oriented
C     all-to-all broadcast
C
      EAST = MOD( MYCOL + 1, NPCOL )
      WEST = MOD( MYCOL - 1 + NPCOL, NPCOL )
C
C     Set starting values for KASE and ITER used in PDLACON
C
      KASE = 0
      ITER = 0
C
C     Set the descriptor for the global distributed vectors
C     vec(X,Y) and vec(C,F), where vec(X,Y) = Z^(-1) * vec(C,F)
C
      LLDV = CONWRK
      CALL DESCINIT( DESCXV, VECDIM, 1, DESCCF( MB_ ), 1, 
     $               DESCCF( RSRC_ ), DESCCF( CSRC_ ), DESCCF( CTXT_ ), 
     $               MAX( 1, LLDV ), INFO )
      XCSRC = DESCXV( CSRC_ )
C
C     Set timers to zero
C
      T1 = ZERO
      T2 = ZERO
      T3 = ZERO
C
C     Enter do-while loop to compute 1-norm lower bound 
C   
      FIRST = .TRUE.
 10   CONTINUE
C
C     Perform all-to-all broadcast of the vector X rowwise in the grid.
C     Store the messages in the order the come in - each process
C     column will perform their own condition estimation call.
C     
      CALL BLACS_BARRIER( ICTXT , 'All' )
      T = MPI_WTIME()
C
C     Copy my own part of the matrix C to the vector X
C
      IF( .NOT. FIRST ) THEN
         CALL DLACPY( 'All', CFROWS, CFCOLS, DWORK(IPC), DESCCF(LLD_),
     $                DWORK( IPX ), CFROWS )
         INDXX = IPX
         XSIZE = CFROWS*CFCOLS
      END IF
      IF (.NOT. FIRST .AND. NPCOL.GT. 1 ) THEN
         DO 20 ROUNDS = 1, NPCOL-1
            CALL DGESD2D( ICTXT, XSIZE, 1, DWORK(INDXX), XSIZE, MYROW,
     $                    EAST )
            INDXX = INDXX + XSIZE
            XPCOL = MOD( MYCOL - ROUNDS + NPCOL, NPCOL ) 
            XCOLS = NUMROC( N + IROFFBE, NB, XPCOL, BECSRC, NPCOL )
            XSIZE = CFROWS * XCOLS
            CALL DGERV2D( ICTXT, XSIZE, 1, DWORK(INDXX), XSIZE, MYROW,
     $                    WEST )
 20      CONTINUE
      END IF
C
C     Copy my own part of the matrix F to the vector Y
C
      IF( .NOT. FIRST ) THEN
         CALL DLACPY( 'All', CFROWS, CFCOLS, DWORK(IPF), DESCCF(LLD_),
     $                DWORK( IPY ), CFROWS )
         INDXY = IPY
         YSIZE = CFROWS*CFCOLS
      END IF
      IF (.NOT. FIRST .AND. NPCOL.GT. 1 ) THEN
         DO 30 ROUNDS = 1, NPCOL-1
            CALL DGESD2D( ICTXT, YSIZE, 1, DWORK(INDXY), YSIZE, MYROW,
     $                    EAST )
            INDXY = INDXY + YSIZE
            YPCOL = MOD( MYCOL - ROUNDS + NPCOL, NPCOL ) 
            YCOLS = NUMROC( N + IROFFBE, NB, YPCOL, BECSRC, NPCOL )
            YSIZE = CFROWS * YCOLS
            CALL DGERV2D( ICTXT, YSIZE, 1, DWORK(INDXY), YSIZE, MYROW,
     $                    WEST )
 30      CONTINUE
      END IF
      T2 = T2 + MPI_WTIME() - T
C
C     Change FIRST
C
      IF( FIRST ) FIRST = .FALSE.
C
C     Compute new estimation (all processors involved)
C
      T = MPI_WTIME()
      DESCXV( CSRC_ ) = MYCOL
      CALL PDLACON( VECDIM, DWORK(IPV), 1, 1, DESCXV, DWORK(IPX), 
     $              1, 1, DESCXV, IWORK(IPIW1), EST, KASE )
      DESCXV( CSRC_ ) = XCSRC
      T1 = T1 + MPI_WTIME() - T
C
C     Make sure that KASE is set to the same for all processors
C
      IF( NPCOL.GT.1 ) THEN
         IF( KASE.EQ.0 ) THEN
            GKASE = 2*NPCOL+1
         ELSE
            GKASE = KASE
         END IF
         CALL IGSUM2D( ICTXT, 'Row', ' ', 1, 1, GKASE, 1, -1, -1, -1, -1,
     $                 -1 )
         GKASE_RST = MOD( GKASE, 2*NPCOL+1 )
         NPCOL0 = GKASE / (2*NPCOL+1)
         IF( GKASE_RST.EQ.0 ) THEN
            KASE = 0
         ELSEIF( GKASE_RST.LT.3*(NPCOL-NPCOL0)/2 ) THEN
            KASE = 1
         ELSE
            KASE = 2
         END IF
      END IF
C
C     Take maximum of all values of EST in each process row
C     
      CALL DGAMX2D( ICTXT, 'Row', ' ', 1, 1, EST, 1, -1, -1, -1, -1,
     $              -1 ) 
C
C     Check whether we should continue or leave the loop
C
      IF( KASE.NE.0 ) THEN
C
C     Update iteration counter
C
         ITER = ITER + 1
C
C     Build the right hand side C from X coming out of the last 
C     call to PDLACON
C
         CALL BLACS_BARRIER( ICTXT , 'All' )
         T = MPI_WTIME()
         CALL DLACPY( 'All', CFROWS, CFCOLS, DWORK(IPX), CFROWS, 
     $                DWORK( IPC ), DESCCF( LLD_ ) )
         CALL DLACPY( 'All', CFROWS, CFCOLS, DWORK(IPY), CFROWS, 
     $                DWORK( IPF ), DESCCF( LLD_ ) )
         T2 = T2 + MPI_WTIME() - T 
C        
C     Choose branch depending on the value of KASE
C
         CALL BLACS_BARRIER( ICTXT , 'All' )
         T = MPI_WTIME()
         IF( KASE.EQ.1 ) THEN
            CALL PGEGCSYD( 'Solve', TRANZ1, 'Schur', 'Schur', TRANAD, 
     $                     TRANBE, ISGN, COMM, M, N, A, IA, JA, DESCA, 
     $                     B, IB, JB, DESCB, DWORK(IPC), ICF, JCF, 
     $                     DESCCF, D, ID, JD, DESCD, E, IE, JE, DESCE, 
     $                     DWORK(IPF), ICF, JCF, DESCCF, MBNB2,
     $                     DWORK(IPW), LDWORK-IPW+1, IWORK(IPIW2), 
     $                     LIWORK-IPIW2+1, NOEXSY, SCALE, INFO )
         ELSEIF( KASE.EQ.2 ) THEN
            CALL PGEGCSYD( 'Solve', TRANZ2, 'Schur', 'Schur', TRANAD, 
     $                     TRANBE, ISGN, COMM, M, N, A, IA, JA, DESCA, 
     $                     B, IB, JB, DESCB, DWORK(IPC), ICF, JCF, 
     $                     DESCCF, D, ID, JD, DESCD, E, IE, JE, DESCE, 
     $                     DWORK(IPF), ICF, JCF, DESCCF, MBNB2, 
     $                     DWORK(IPW), LDWORK-IPW+1, IWORK(IPIW2), 
     $                     LIWORK-IPIW2+1, NOEXSY, SCALE, INFO )
         END IF
         T3 = T3 + MPI_WTIME() - T 
C
         GO TO 10
      END IF
C     
C     Compute lower bound on spectral (2-) norm 
C
      EST = EST / (SQRT(DBLE(2*M*N))*SCALE)
C
C     Give number of iterations to output
C 
      NOITER = ITER
C     
C     Save timings to DWORK
C
      DWORK( 1 ) = T1
      DWORK( 2 ) = T2
      DWORK( 3 ) = T3
C
      END
C     
C ***Last line of PGCSYCON***
