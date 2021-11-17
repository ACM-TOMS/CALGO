CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C         This file is a part of the SCASY software library.           C
C         See http://www.cs.umu.se/parallel/research/scasy             C
C         Contributors: Robert Granat and Bo Kågström.                 C
C         SCASY Version 0.10,    March 31, 2006.                       C
C         Copyright 2006, Umeå University, Sweden.                     C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE PTRGCSYD( TRANZ, TRANAD, TRANBE, ISGN, COMM, M, N, A, 
     $                     IA, JA, DESCA, B, IB, JB, DESCB, C, IC, JC, 
     $                     DESCC, D, ID, JD, DESCD, E, IE, JE, DESCE, F, 
     $                     IF, JF, DESCF, MBNB2, DWORK, LDWORK, IWORK, 
     $                     LIWORK, NOEXSY, SCALE, INFO )
C
C  -- ScaLAPACK-style routine (preliminary version ) --
C     HPC2N and Department of Computing Science,
C     University of Umeå, Sweden.
C     Written by Robert Granat, (granat@cs.umu.se)
C     September 4, 2006.
C
      IMPLICIT NONE
C
C     .. Scalar Arguments ..
      CHARACTER*1        TRANZ, TRANAD, TRANBE, COMM
      INTEGER            ISGN, M, N, IA, JA, IB, JB, IC, JC, ID, JD, 
     $                   IE, JE, IF, JF, INFO, LDWORK, LIWORK, NOEXSY
      DOUBLE PRECISION   SCALE
C     ..
C     .. Array Arguments ..
      INTEGER            DESCA( * ), DESCD( * ), DESCB( * ), DESCE( * ),
     $                   DESCC( * ), DESCF( * ), IWORK( * ), MBNB2( 2 )
      DOUBLE PRECISION   A( * ), D( * ), B( * ), E( * ), C( * ), F( * ),
     $                   DWORK( * )
C     ..
C
C 
C  Purpose
C  =======
C  This subroutine solves the real (quasi-)triangular generalized 
C  coupled Sylvester Equation (GCSY)
C
C    ( op( sub( A ) ) * X +/- Y * op( sub( B ) ) = C,
C      op( sub( D ) ) * X +/- Y * op( sub( E ) ) = F )   (1)
C
C  where sub(A) = A(IA:IA+M-1,JA:JA+M-1) and 
C  sub(D) = D(ID:ID+M-1,JD:JD+M-1) are M-by-M matrices, 
C  sub(B) = B(IB:IB+N-1,JB:JB+N-1) and sub(E) = E(IE:IE+N-1,JE:JE+N-1) 
C  are N-by-N matrices and the solution matrices X and Y are M-by-N which
C  overwrites C(IC:IC+M-1,JC:JC+N-1) and F(IF:IF+M-1,JF:JF+N-1). 
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
C  Z'*y = scale*b, which is equivalent to solve for R and L in
C
C              A' * X  +  D' * Y  = scale *  C           (3)
C              X  * B' +  Y  * E' = scale * ( +/- F )
C
C  This case (TRANZ = 'T') is used to compute a one-norm-based estimate
C  of Dif[(A,D),(B,E)], the separation between the matrix pairs (A,D)
C  and (B,E), using PDLACON in the routine PGCSYCON.
C
C  This routine should *not* be called directly, but through PGEGCSYD.
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
C  TRANZ     (global input) CHARACTER*1
C            If TRANZ = 'T', this routines solves for the equation (2).
C            Notice that TRANZ = 'T' is only allowed for 
C            TRANAD = TRANBE = 'N' or TRANAD = TRANBE = 'T'.
C            If TRANZ = 'N', this routines solves equation (1), and all
C            possible combination of the transposes TRANAD and TRANBE
C            are allowed.
C
C  TRANAD    (global input) CHARACTER*1
C              If TRANAD = 'N' then op(A) = A and op(D) = D
C              If TRANAD = 'T' then op(A) = A**T and op(D) = D**T
C
C  TRANBE    (global input) CHARACTER*1
C              If TRANBE = 'N' then op(B) = B and op(E) = E
C              If TRANBE = 'T' then op(B) = B**T and op(E) = E**T
C
C  ISGN      (global input) INTEGER*1
C            If ISGN = 1, we solve the equation (1) with a '+'.
C            If ISGN = -1, we solve the equation (1) with a '-'.
C
C  Input/output arguments
C
C  COMM      (global input/output) CHARACTER*1
C            This subroutine uses two different communications schemes in
C            solving the reduced triangular problem:
C              If COMM = 'S', the "shifts" scheme is used.
C              If COMM = 'D', the "communicate on demand" scheme is used.
C            The choice COMM = 'S' is only valid for TRANAD = TRANBE = 'N' 
C            or TRANAD = TRANBE = 'T'. The scheme used will be output.
C            See the references for details.
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
C            Array of dimension (LLD_C,LOCc(N)). 
C            On entry C contains the local pieces of the global 
C            distributed matrix C. On exit, it contains the local 
C            pieces of the global distributed solution X.
C
C  IC        (global input) INTEGER
C            Row start index for sub(C), i.e., the submatrix to operate 
C            on. MOD(IC,MB_A) = MOD(IA,MB_A) must hold. 
C
C  JC        (global input) INTEGER
C            Column start index for sub(C), i.e., the submatrix to operate 
C            on. MOD(JC,MB_B) = MOD(JB,MB_B) must hold.  
C
C  DESCC     (global and local input) INTEGER array of dimension DLEN_.
C            The array descriptor for the global distributed matrix C.
C
C  D         (local input) DOUBLE PRECISION array 
C            Array of dimension (LLD_D,LOCc(M)). Contains the local
C            pieces of the global distributed matrix D. 
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
C  E         (local input) DOUBLE PRECISION array 
C            Array of dimension (LLD_E,LOCc(N)). Contains the local
C            pieces of the global distributed matrix E. 
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
C            On entry F contains the local pieces of the global 
C            distributed matrix F. On exit, it contains the local 
C            pieces of the global distributed solution Y.
C 
C  IF        (global input) INTEGER
C            Row start index for sub(F), i.e., the submatrix to operate 
C            on. IF = IC must hold. 
C
C  JF        (global input) INTEGER
C            Column start index for sub(F), i.e., the submatrix to operate 
C            on. JF = JC must hold.  
C
C  DESCF     (global and local input) INTEGER array of dimension DLEN_.
C            The array descriptor for the global distributed matrix F.
C
C  MBNB2     (global input) INTEGER array of dimension 2.
C            Internal blocking factors for pipelining of subsolutions
C            for updates of the matrix C in PTRSYCTD (see the references 
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
C            LIWORK >= DBA + DBB + MB_A + MB_B + 
C            10 * MIN( P_r, P_C ) + 8, where 
C            DBA = ICEIL(LOCr(IA+IROFFA),MB_A) and DBB = 
C            = ICEIL(LOCr(IB+IROFFB),MB_B).
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
C  This routine implements a parallel wave-front algorithm for solving
C  the triangular generalized coupled Sylvester equation. See the 
C  references for details.
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
C  Wave-front algorithm, generalized Schur form, generalized coupled
C  Sylvester equation, explicit blocking, GEMM-updates, matrix shifts,
C  on demand
C
C  =====================================================================
C     ..
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
      INTEGER          MB, NB, DBAD, DBB, NROLL, IS, JS, IEND, K, RCDIR,
     $                 MYCOL, MYROW, NPCOL, NPROW, J, NPROCS, I, IDUM,
     $                 ADROWS, ADCOLS, BEROWS, BECOLS, ROWS, COLS, LINFO, 
     $                 D, DA, DB, IX, JX, RSRC, CSRC, LIAD, LJAD, LIBE, 
     $                 LJBE, INDX, ICTXT, LLDA, LLDB, LLDC, NORTH, WEST,
     $                 SOUTH, EAST, SRSRC, SCSRC, WRK, LMATR, FIJ, 
     $                 XIJ1, YIJ1, XIJ2, YIJ2, ADII, BEJJ, CIJ, MWORK, 
     $                 RRSRC, IIS, IIE, JJS, PHASES, PHASE, MNPDIM, 
     $                 RCSRC, ADRSRC, ADCSRC, BERSRC, BECSRC, CFRSRC, 
     $                 CFCSRC, EXA, EXD, EXB, EXE, EXC, EXF, EXMEMAD, 
     $                 EXMEMBE, EXMEMCF, IWRK, EXADINF, EXBEINF, SND, 
     $                 LEXBUFF, GI, GJ, LBI, LBJ, NBCAD, NBCBE, 
     $                 POS, ADROWS2, BECOLS2, GINDX, ADCOLS2, BEROWS2,
     $                 IW, LINFO, DD, LLDD, LLDE, LLDF, DBBE,
     $                 IROFFAD, IROFFBE, LICF, LJCF, ASI, ASJ, BSI, 
     $                 BSJ, CSI, CSJ, DSI, DSJ, ESI, ESJ, FSI, FSJ,
     $                 LIAS, LJAS, LIBS, LJBS, LICS, LJCS, LIES, LJES,
     $                 LIFS, LJFS, LIDS, LJDS, GSI, GSJ, GSIND, INDXS,
     $                 INDXE, INDXU, IROWS, ICOLS, IGSI, IGSJ, IEXRW, 
     $                 IEXCL, MB2, NB2, NIDEEP, NJDEEP, IDEEPS, JDEEPS, 
     $                 IDEEPE, JDEEPE, IDEEPU, JDEEPU, IDEEP, JDEEP, 
     $                 XYRWS, XYCLS, XYRIND, XYCIND, IXYRWS, IXYCLS, 
     $                 IXYRIND, IXYCIND, ADUPBL, BEUPBL, CKJ, FKJ,
     $                 ADUP, BEUP, KADUP, KBEUP, IPW, KKK
      DOUBLE PRECISION SCALOC, USIGN, DIF, DUMMY, LSCALC
      LOGICAL          LQUERY, SNODE, TRANSAD, TRANSBE, EXROW, EXCOL, 
     $                 CEXROW, CEXCOL, ADEXT, BEEXT, SHIFT, TRANSZ
C     ..
C     .. Local Arrays ..
      INTEGER          IDUM1(1), IDUM2(1), DESCSA( DLEN_ ),
     $                 DESCSB( DLEN_ ), DESCSC( DLEN_ ), 
     $                 DESCSF( DLEN_ ), DESCSD( DLEN_ ),
     $                 DESCSE( DLEN_ ), IBUFF( 4 )

      DOUBLE PRECISION MACHINE(10)
C     ..
C     .. External Subroutines ..
      EXTERNAL         BLACS_GRIDINFO, CHK1MAT, DTGSYL, PXERBLA, 
     $                 INFOG2L, DGEMM, DLACPY, PCHK2MAT, DSCAL, 
     $                 DGAMN2D, DGESD2D, DGERV2D, IGAMX2D, DMATADD,
     $                 PDLACPY, DLATCPY, RECGCSY, RECSY_MACHINE
C     ..
C     .. External Functions ..
      LOGICAL          LSAME, INT2LG
      INTEGER          NUMROC, ICEIL, ILG2NT
      EXTERNAL         LSAME, NUMROC, ICEIL, INT2LG, ILG2NT
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC        MOD, MAX, MIN 
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
         INFO = 1
      ELSE
         CALL CHK1MAT( M, 6, M, 6, IA, JA, DESCA, 11, INFO )
         CALL CHK1MAT( N, 7, N, 7, IB, JB, DESCB, 15, INFO )
         CALL CHK1MAT( M, 6, N, 7, IC, JC, DESCC, 19, INFO )
         CALL CHK1MAT( M, 6, M, 6, ID, JD, DESCD, 23, INFO )
         CALL CHK1MAT( N, 7, N, 7, IE, JE, DESCB, 27, INFO )
         CALL CHK1MAT( M, 6, N, 7, IC, JC, DESCC, 31, INFO )
         CALL PCHK2MAT( M, 6, M, 6, IA, JA, DESCA, 11, N, 7, N, 7, IB,
     $                 JB, DESCB, 15, 0, IDUM1, IDUM2, INFO )
         CALL PCHK2MAT( M, 6, M, 6, ID, JD, DESCD, 23, N, 7, N, 7, IE,
     $                 JE, DESCE, 27, 0, IDUM1, IDUM2, INFO )
         CALL PCHK2MAT( M, 6, M, 6, IA, JA, DESCA, 10, M, 6, N, 7, IC, 
     $                  JC, DESCC, 19, 0, IDUM1, IDUM2, INFO )
         CALL PCHK2MAT( M, 6, M, 6, ID, JD, DESCD, 23, M, 6, N, 7, IF, 
     $                  JF, DESCF, 31, 0, IDUM1, IDUM2, INFO )
      END IF
C
      IF( INFO.EQ.0 ) THEN
         IF( .NOT.( LSAME( TRANAD, 'T' ) .OR. 
     $        LSAME( TRANAD, 'N' ) ) ) THEN
            INFO = -1
         END IF
      END IF
      IF( INFO.EQ.0 ) THEN
         IF( .NOT.( LSAME( TRANBE, 'T' ) .OR. 
     $        LSAME( TRANBE, 'N' ) ) ) THEN
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
         IF( .NOT. ( (LSAME( TRANAD, 'N' ).AND.LSAME( TRANBE, 'N' ) )
     $        .OR.( (LSAME( TRANAD, 'T' ).AND.LSAME( TRANBE, 'T' )))))
     $        COMM = 'D'
      END IF
C
C     Check that the last row and column of the matrices sub(A), sub(B) 
C     and sub(C) can be mapped onto the last process row or column. 
C     If not, change to on-demand communication.
C     
      MB = DESCA( MB_ )
      NB = DESCB( MB_ ) 
      IROFFAD = MOD( IA - 1, MB )
      IROFFBE = MOD( IB - 1, NB ) 
      DBAD = ICEIL( M + IROFFAD, DESCA(MB_) )
      DBBE = ICEIL( N + IROFFBE, DESCB(MB_) )
C
      IF( INFO.EQ.0 .AND. LSAME(COMM,'S')) THEN 
         DD = MAX( NPROW, NPCOL )
         IF( MIN( M + IROFFAD, N + IROFFBE ) .GT. (DD-1)**2 ) THEN
            IF( .NOT.(MOD(DBAD,DD).EQ.0 .OR. MOD(DBBE,DD).EQ.0) ) THEN
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
      TRANSAD = LSAME( TRANAD, 'T' )
      TRANSBE = LSAME( TRANBE, 'T' )
      TRANSZ = LSAME( TRANZ, 'T' )
      SHIFT = LSAME( COMM, 'S' )
      ADEXT = .FALSE.
      BEEXT = .FALSE.
      NOEXSY = 0
      LLDA = DESCA(LLD_)
      LLDD = DESCD(LLD_)
      LLDB = DESCB(LLD_)
      LLDE = DESCE(LLD_)
      LLDC = DESCC(LLD_)
      LLDF = DESCF(LLD_)
C
C     Compute the number of blocks to keep from A and B during the 
C     updates when deep pipelining is turned on
C
      MB2 = MBNB2(1)
      NB2 = MBNB2(2)
      IF( MB2.NE.MB .OR. NB2.NE.NB ) THEN
         ADUPBL = 2 * MAX( 1, ICEIL( DBAD-1, MNPDIM ) ) + 2
         BEUPBL = 2 * MAX( 1, ICEIL( DBBE-1, MNPDIM ) ) + 2
      ELSE
         ADUPBL = 1
         BEUPBL = 1
      END IF     
C
C     Compute the number of rows and columns held by each process and
C     the number of block columns held by each process for the smallest
C     submatrices of A and B including sub(A) and sub(B) conforming with 
C     ScaLAPACK conventions 
C
      CALL INFOG2L( IA, JA, DESCA, NPROW, NPCOL, MYROW, MYCOL, LIAD, 
     $              LJAD, ADRSRC, ADCSRC )
      CALL INFOG2L( IB, JB, DESCB, NPROW, NPCOL, MYROW, MYCOL, LIBE, 
     $              LJBE, BERSRC, BECSRC )
      CALL INFOG2L( IC, JC, DESCC, NPROW, NPCOL, MYROW, MYCOL, LICF, 
     $              LJCF, CFRSRC, CFCSRC )
      ADROWS = NUMROC( M + IROFFAD, MB, MYROW, ADRSRC, NPROW )
      ADCOLS = NUMROC( M + IROFFAD, MB, MYCOL, ADCSRC, NPCOL )
      BEROWS = NUMROC( N + IROFFBE, NB, MYROW, BERSRC, NPROW )
      BECOLS = NUMROC( N + IROFFBE, NB, MYCOL, BECSRC, NPCOL )
      NBCAD = ICEIL( ADCOLS, MB )
      NBCBE = ICEIL( BECOLS, NB )
C
C     Compute the needed memory for holding data for the extended
C     subsystems. Even if all blocks cannot be extended, we add memory
C     for all blocks to maintain indicies simple.
C
      EXMEMAD = 2*ICEIL(ADROWS,MB)*ICEIL(ADCOLS,MB)*(2*MB+1)
      EXMEMBE = 2*ICEIL(BEROWS,NB)*ICEIL(BECOLS,NB)*(2*NB+1)
      EXMEMCF = 2*ICEIL(ADROWS,MB)*ICEIL(BECOLS,NB)*(MB+NB+1) 
C
C     Adjust AROWS - BCOLS to sub(A) and sub(B) 
C
      IF( MYROW.EQ.ADRSRC ) ADROWS = ADROWS - IROFFAD
      IF( MYCOL.EQ.ADCSRC ) ADCOLS = ADCOLS - IROFFAD
      IF( MYROW.EQ.BERSRC ) BEROWS = BEROWS - IROFFBE
      IF( MYCOL.EQ.BECSRC ) BECOLS = BECOLS - IROFFBE
C
C     We need six extra matrix descriptors for shifts and for 
C     calculating the correct indicies in sub(A), sub(B)...
C     These descriptors are related to the the smallest submatrices of 
C     including sub(_) conforming with ScaLAPACK conventions.
C    
      CALL DESCINIT ( DESCSA, M+IROFFAD, M+IROFFAD, DESCA(MB_), 
     $                DESCA(NB_), ADRSRC, ADCSRC, ICTXT, LLDA, INFO )
      CALL DESCINIT ( DESCSB, N+IROFFBE, N+IROFFBE, DESCB(MB_), 
     $                DESCB(NB_), BERSRC, BECSRC, ICTXT, LLDB, INFO )
      CALL DESCINIT ( DESCSC, M+IROFFAD, N+IROFFBE, DESCC(MB_), 
     $                DESCC(NB_), CFRSRC, CFCSRC, ICTXT, LLDC, INFO )
      CALL DESCINIT ( DESCSD, M+IROFFAD, M+IROFFAD, DESCD(MB_), 
     $                DESCD(NB_), ADRSRC, ADCSRC, ICTXT, LLDD, INFO )
      CALL DESCINIT ( DESCSE, N+IROFFBE, N+IROFFBE, DESCE(MB_), 
     $                DESCE(NB_), BERSRC, BECSRC, ICTXT, LLDE, INFO )
      CALL DESCINIT ( DESCSF, M+IROFFAD, N+IROFFBE, DESCF(MB_), 
     $                DESCF(NB_), CFRSRC, CFCSRC, ICTXT, LLDF, INFO )
C
C     Test work space
C     
C     
      WRK = 8*(MB+1)*(NB+1) + 2*(MB+1)**2 +
     $     2*(NB+1)**2 + EXMEMAD + EXMEMBE + EXMEMCF +
     $     MAX( MB, NB ) + ADUPBL * (MB+1) ** 2 + 
     $     BEUPBL * (NB+1) ** 2 + 2 * (MB2+1) * (NB2+1)
C     
C     Compute needed integer workspace
C     
      IF( INFO.EQ.0 ) THEN
         IWRK = DBAD + DBBE + MB + NB + 10 * MNPDIM + 8 
      END IF
C     
C     Check if the call to PTRGCSYD was a workspace query, if so
C     store the needed workspace in DWORK(1) and IWORK(1) and return
C     If not check if the supplied memory is big enough
C     
      IF( INFO.EQ.0 .OR. LQUERY ) THEN
         IF( WRK.GT.LDWORK .AND. .NOT.LQUERY ) THEN
            INFO = -34
         ELSEIF( IWRK.GT.LIWORK .AND. .NOT. LQUERY ) THEN
            INFO = -35
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
         CALL PXERBLA( ICTXT, 'PTRGCSYD', -INFO )
         RETURN
      END IF
C
C     Set the sign in the BE-updates according to ISGN: if ISGN = 1 we 
C     update by subtraction, and if ISGN = -1 we update by addition.
C     
      USIGN = DBLE(-1*ISGN)
C     
C     Init some local pointers into the DWORK-array
C     
      XIJ1 = 1
      YIJ1 = XIJ1 + (MB+1) * (NB+1)
      XIJ2 = YIJ1 + (MB+1) * (NB+1)
      YIJ2 = XIJ2 + (MB+1) * (NB+1)
      ADII = YIJ2 + (MB+1) * (NB+1)
      BEJJ = ADII + 2 * (MB+1) ** 2
      CIJ = BEJJ + 2 * (NB+1) ** 2
      FIJ = CIJ + (MB+1) * (NB+1)
      SND = FIJ + (MB+1) * (NB+1)
      EXA = SND + MAX( MB, NB )
      EXD = EXA + EXMEMAD
      EXB = EXD + EXMEMAD
      EXE = EXB + EXMEMBE
      EXC = EXE + EXMEMBE
      EXF = EXC + EXMEMCF
      CKJ = EXF + EXMEMCF
      FKJ = CKJ + (MB+1) * (NB+1)
      ADUP = FKJ + (MB+1) * (NB+1)
      BEUP = ADUP + ADUPBL * (MB+1) ** 2
      IPW = BEUP + BEUPBL * (NB+1) ** 2
C
C     Init two local pointers into the IDWORK-array. Those should hold
C     the startindicies of two integer array describing the extensions
C     of the submatrices of (A,D), (B,E) and (C,F), respectively.
C
      EXADINF = 1
      EXBEINF = EXADINF + DBAD
      IROWS = EXBEINF + DBBE
      ICOLS = IROWS + MNPDIM
      IGSI  = ICOLS + MNPDIM
      IGSJ  = IGSI + MNPDIM
      IEXRW = IGSJ + MNPDIM
      IEXCL = IEXRW + MNPDIM
      IXYRWS = IEXCL + MNPDIM
      IXYCLS = IXYRWS + MNPDIM
      IXYRIND = IXYCLS + MNPDIM
      IXYCIND = IXYRIND + MNPDIM
      IW = IXYCIND + MNPDIM
C
C     Check for 2x2 diagonal-block-split between any blocks of op(A) or
C     op(B) and set the extensions.
C 
      CALL PDEXTCHK( M, A, IA, JA, DESCA, IWORK( EXADINF ), DBAD, ADEXT, 
     $               INFO )
C     
      CALL PDEXTCHK( N, B, IB, JB, DESCB, IWORK( EXBEINF ), DBBE, BEEXT, 
     $               INFO )
C     
C     Do an implicit redistribution of the elements in (A,D), (B,E) and 
C     (C,F) based on the extension information just set
C
      IF( ADEXT.OR.BEEXT ) 
     $     CALL PDIMPRED( 'ABC', M, N, A, IA, JA, DESCA, B, IB, JB,
     $                    DESCB, C, IC, JC, DESCC, IWORK( EXADINF ),
     $                    IWORK( EXBEINF ),  DWORK( EXA ), EXMEMAD, 
     $                    DWORK( EXB ), EXMEMBE, DWORK( EXC ), EXMEMCF, 
     $                    DWORK( SND ), MAX( MB, NB ), INFO )
      IF( ADEXT.OR.BEEXT ) 
     $     CALL PDIMPRED( 'ABC', M, N, D, ID, JD, DESCD, E, IE, JE,
     $                    DESCE, F, IF, JF, DESCF, IWORK( EXADINF ),
     $                    IWORK( EXBEINF ),  DWORK( EXD ), EXMEMAD, 
     $                    DWORK( EXE ), EXMEMBE, DWORK( EXF ), EXMEMCF, 
     $                    DWORK( SND ), MAX( MB, NB ), INFO )   
C     
C     Compute the number of block diagonals of (C,F) and some 
C     communication directions
C     
      NROLL = DBAD + DBBE - 1
      NORTH = MOD(MYROW + NPROW - 1, NPROW)
      EAST =  MOD(MYCOL + 1, NPCOL)
      SOUTH = MOD(MYROW + 1, NPROW)
      WEST =  MOD(MYCOL + NPCOL - 1, NPCOL)
C 
C     Depending on the transposes, set some loop variables
C     
      IF( .NOT. TRANSZ ) THEN
         IF( (.NOT.TRANSAD).AND.(.NOT.TRANSBE) ) THEN
            JS = 1
            IS = DBAD
            IEND = DBAD       
         ELSEIF( TRANSAD.AND.(.NOT.TRANSBE) ) THEN
            JS = 1
            IS = 1
            IEND = 1
         ELSEIF( (.NOT.TRANSAD).AND.TRANSBE ) THEN
            JS = DBBE
            IS = DBAD
            IEND = DBAD
         ELSEIF( TRANSAD.AND.TRANSBE ) THEN
            JS = DBBE
            IS = 1
            IEND = 1
         END IF
      ELSE
         IF( (.NOT.TRANSAD).AND.(.NOT.TRANSBE) ) THEN
            JS = DBBE
            IS = 1
            IEND = 1   
         ELSEIF( TRANSAD.AND.TRANSBE ) THEN
            JS = 1
            IS = DBAD
            IEND = DBAD  
         END IF
      END IF
C     
      SNODE = .TRUE.
C
C     Compute the number of rounds to do in deep pipelining and
C     set looplimits depending on the transposes
C
      NIDEEP = ICEIL( MB, MB2 )
      NJDEEP = ICEIL( NB, NB2 ) 
      IF( .NOT.TRANSZ ) THEN
         IF( .NOT.TRANSAD  ) THEN
            IDEEPS = NIDEEP
            IDEEPE = 1
            IDEEPU = -1
         ELSE
            IDEEPS = 1
            IDEEPE = NIDEEP
            IDEEPU = 1
         END IF
         IF( .NOT.TRANSBE  ) THEN
            JDEEPS = 1
            JDEEPE = NJDEEP
            JDEEPU = 1
         ELSE
            JDEEPS = NJDEEP
            JDEEPE = 1
            JDEEPU = -1
         END IF 
      ELSE
         IF( .NOT.TRANSAD .AND. .NOT.TRANSBE ) THEN
            IDEEPS = 1
            IDEEPE = NIDEEP
            IDEEPU = 1
            JDEEPS = NJDEEP
            JDEEPE = 1
            JDEEPU = -1
         ELSEIF( TRANSAD .AND. TRANSBE ) THEN
            IDEEPS = NIDEEP
            IDEEPE = 1
            IDEEPU = -1
            JDEEPS = 1
            JDEEPE = NJDEEP
            JDEEPU = 1
         END IF 
      END IF
C
C     Compute the local indicies where my part of sub(_) begins
C
      ASI = IA + MB * MOD( NPROW + MYROW - ADRSRC, NPROW ) - IROFFAD 
      ASJ = JA + MB * MOD( NPCOL + MYCOL - ADCSRC, NPCOL ) - IROFFAD
      BSI = IB + NB * MOD( NPROW + MYROW - BERSRC, NPROW ) - IROFFBE
      BSJ = JB + NB * MOD( NPCOL + MYCOL - BECSRC, NPCOL ) - IROFFBE
      CSI = IC + MB * MOD( NPROW + MYROW - CFRSRC, NPROW ) - IROFFAD
      CSJ = JC + NB * MOD( NPCOL + MYCOL - CFCSRC, NPCOL ) - IROFFBE
      DSI = ID + MB * MOD( NPROW + MYROW - ADRSRC, NPROW ) - IROFFAD 
      DSJ = JD + MB * MOD( NPCOL + MYCOL - ADCSRC, NPCOL ) - IROFFAD
      ESI = IE + NB * MOD( NPROW + MYROW - BERSRC, NPROW ) - IROFFBE
      ESJ = JE + NB * MOD( NPCOL + MYCOL - BECSRC, NPCOL ) - IROFFBE
      FSI = IF + MB * MOD( NPROW + MYROW - CFRSRC, NPROW ) - IROFFAD
      FSJ = JF + NB * MOD( NPCOL + MYCOL - CFCSRC, NPCOL ) - IROFFBE
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
      CALL INFOG2L( FSI, FSJ, DESCF, NPROW, NPCOL, MYROW, MYCOL, 
     $              LIFS, LJFS, IDUM1, IDUM2 )
C     
C     Main loop over number of diagonals in C
C     
      DO 10 K = 1, NROLL
C        
         IF( .NOT. TRANSZ ) THEN
            IF( (.NOT.TRANSAD).AND.(.NOT.TRANSBE) ) THEN
               IF ( K.GT.DBBE) IS = IS - 1
            ELSEIF( TRANSAD.AND.(.NOT.TRANSBE) ) THEN
               IF ( K.GT.DBBE) IEND = IEND + 1
            ELSEIF( (.NOT.TRANSAD).AND.TRANSBE ) THEN
               IF ( K.GT.DBBE) IS = IS - 1
            ELSEIF( TRANSAD.AND.TRANSBE ) THEN
               IF ( K.GT.DBBE) IEND = IEND + 1
            END IF
         ELSE
            IF( (.NOT.TRANSAD).AND.(.NOT.TRANSBE) ) THEN
               IF ( K.GT.DBBE) IEND = IEND + 1
            ELSEIF( TRANSAD.AND.TRANSBE ) THEN
               IF ( K.GT.DBBE) IS = IS - 1 
            END IF
         END IF
C
C     If shifts are supposed to be used, we do it right here
C
         IF (SHIFT .AND. NPROCS.GT.1 .AND. 
     $        (M.EQ.N.OR.N.EQ.NB.OR.M.EQ.MB)) THEN            
C     
C     Shift (A,D) East if NPCOL > 1 and (A,D)**T West if NPCOL > 1
C     Shift (B,E) North if NPROW > 1 and (B,E)**T South if NPROW > 1
C     When shifting, also send/receive the extension elements
C
            IF( TRANSZ ) THEN
               TRANSAD = .NOT. TRANSAD
               TRANSBE = .NOT. TRANSBE
            END IF
            IF ( NPCOL.GT.1.AND.M.NE.MB ) THEN
               IF( (.NOT.TRANSAD) .AND. (.NOT.TRANSBE) ) THEN
                  CALL DGESD2D( ICTXT, ADROWS, ADCOLS, 
     $                 A((LJAS-1)*LLDA+LIAS), LLDA, MYROW,
     $                 EAST )
                  DESCSA(CSRC_) = MOD(DESCSA(CSRC_) + 1, NPCOL)
                  CALL DGERV2D( ICTXT, ADROWS, ADCOLS,
     $                 A((LJAS-1)*LLDA+LIAS), LLDA, MYROW,
     $                 WEST )
                  CALL DGESD2D( ICTXT, ADROWS, ADCOLS,
     $                 D((LJDS-1)*LLDD+LIDS), LLDD, MYROW,
     $                 EAST )
                  DESCSD(CSRC_) = MOD(DESCSD(CSRC_) + 1, NPCOL)
                  CALL DGERV2D( ICTXT, ADROWS, ADCOLS,
     $                 D((LJDS-1)*LLDD+LIDS), LLDD, MYROW,
     $                 WEST )
                  IF( ADEXT ) THEN
                     CALL DGESD2D( ICTXT, EXMEMAD, 1, DWORK(EXA),
     $                    EXMEMAD, MYROW, EAST )
                     CALL DGERV2D( ICTXT, EXMEMAD, 1, DWORK(EXA),
     $                    EXMEMAD, MYROW, WEST )
                     CALL DGESD2D( ICTXT, EXMEMAD, 1, DWORK(EXD),
     $                    EXMEMAD, MYROW,EAST )
                     CALL DGERV2D( ICTXT, EXMEMAD, 1, DWORK(EXD), 
     $                    EXMEMAD, MYROW, WEST )
                  END IF
               ELSEIF( TRANSAD .AND. TRANSBE ) THEN
                  CALL DGESD2D( ICTXT, ADROWS, ADCOLS,
     $                 A((LJAS-1)*LLDA+LIAS), LLDA, MYROW,
     $                 WEST )  
                  DESCSA(CSRC_) = MOD(DESCSA(CSRC_) + NPCOL - 1, 
     $                 NPCOL)
                  CALL DGERV2D( ICTXT, ADROWS, ADCOLS,
     $                 A((LJAS-1)*LLDA+LIAS), LLDA, MYROW,
     $                 EAST )
                  CALL DGESD2D( ICTXT, ADROWS, ADCOLS,
     $                 D((LJDS-1)*LLDD+LIDS), LLDD, MYROW,
     $                 WEST )  
                  DESCSD(CSRC_) = MOD(DESCSD(CSRC_) + NPCOL - 1, 
     $                 NPCOL)
                  CALL DGERV2D( ICTXT, ADROWS, ADCOLS,
     $                 D((LJDS-1)*LLDD+LIDS), LLDD, MYROW,
     $                 EAST )
                  IF( ADEXT ) THEN
                     CALL DGESD2D( ICTXT, EXMEMAD, 1, DWORK(EXA),
     $                    EXMEMAD, MYROW, WEST )
                     CALL DGERV2D( ICTXT, EXMEMAD, 1, DWORK(EXA),
     $                    EXMEMAD, MYROW, EAST )
                     CALL DGESD2D( ICTXT, EXMEMAD, 1, DWORK(EXD),
     $                    EXMEMAD, MYROW, WEST )
                     CALL DGERV2D( ICTXT, EXMEMAD, 1, DWORK(EXD),
     $                    EXMEMAD, MYROW, EAST )
                  END IF
               END IF
            END IF
            IF ( NPROW.GT.1.AND.N.NE.NB ) THEN
               IF( (.NOT.TRANSAD) .AND. (.NOT.TRANSBE) ) THEN
                  CALL DGESD2D( ICTXT, BEROWS, BECOLS, 
     $                 B((LJBS-1)*LLDB+LIBS), LLDB, NORTH,
     $                 MYCOL )
                  DESCSB(RSRC_) = MOD( DESCSB(RSRC_) + NPROW - 1,
     $                 NPROW )
                  CALL DGERV2D( ICTXT, BEROWS, BECOLS,
     $                 B((LJBS-1)*LLDB+LIBS), LLDB, SOUTH,
     $                 MYCOL )
                  CALL DGESD2D( ICTXT, BEROWS, BECOLS,
     $                 E((LJES-1)*LLDE+LIES), LLDE, NORTH,
     $                 MYCOL )
                  DESCSE(RSRC_) = MOD( DESCSE(RSRC_) + NPROW - 1,
     $                 NPROW )
                  CALL DGERV2D( ICTXT, BEROWS, BECOLS,
     $                 E((LJES-1)*LLDE+LIES), LLDE, SOUTH,
     $                 MYCOL )
                  IF( BEEXT ) THEN
                     CALL DGESD2D( ICTXT, EXMEMBE, 1, DWORK(EXB),
     $                    EXMEMBE, NORTH, MYCOL )
                     CALL DGERV2D( ICTXT, EXMEMBE, 1, DWORK(EXB),
     $                    EXMEMBE, SOUTH, MYCOL )
                     CALL DGESD2D( ICTXT, EXMEMBE, 1, DWORK(EXE),
     $                    EXMEMBE, NORTH, MYCOL )
                     CALL DGERV2D( ICTXT, EXMEMBE, 1, DWORK(EXE),
     $                    EXMEMBE, SOUTH, MYCOL )
                  END IF
               ELSEIF( TRANSAD .AND. TRANSBE ) THEN
                  CALL DGESD2D( ICTXT, BEROWS, BECOLS,
     $                 B((LJBS-1)*LLDB+LIBS), LLDB, SOUTH,
     $                 MYCOL )
                  DESCSB(RSRC_) = MOD(DESCSB(RSRC_) + 1, NPROW)
                  CALL DGERV2D( ICTXT, BEROWS, BECOLS,
     $                 B((LJBS-1)*LLDB+LIBS), LLDB, NORTH,
     $                 MYCOL )
                  CALL DGESD2D( ICTXT, BEROWS, BECOLS, 
     $                 E((LJES-1)*LLDE+LIES), LLDE, SOUTH,
     $                 MYCOL )
                  DESCSE(RSRC_) = MOD(DESCSE(RSRC_) + 1, NPROW)
                  CALL DGERV2D( ICTXT, BEROWS, BECOLS,
     $                 E((LJES-1)*LLDE+LIES), LLDE, NORTH,
     $                 MYCOL )
                  IF( BEEXT ) THEN
                     CALL DGESD2D( ICTXT, EXMEMBE, 1, DWORK(EXB),
     $                    EXMEMBE, SOUTH, MYCOL )
                     CALL DGERV2D( ICTXT, EXMEMBE, 1, DWORK(EXB),
     $                    EXMEMBE, NORTH, MYCOL )
                     CALL DGESD2D( ICTXT, EXMEMBE, 1, DWORK(EXE),
     $                    EXMEMBE, SOUTH, MYCOL )
                     CALL DGERV2D( ICTXT, EXMEMBE, 1, DWORK(EXE),
     $                    EXMEMBE, NORTH, MYCOL )
                  END IF
               END IF
            END IF
         ELSEIF (SHIFT .AND. NPROCS.GT.1 .AND. M.LT.N) THEN 
C     
C     Shift (A,D) SouthEast and (C,F) South  if NPROW > 1 
C     or Shift (A,D)**T NorthWest and (C,F) North if NPROW > 1 
C     
            IF( TRANSZ ) THEN
               TRANSAD = .NOT. TRANSAD
               TRANSBE = .NOT. TRANSBE
            END IF
            IF( (.NOT.TRANSAD) .AND. (.NOT.TRANSBE) ) THEN
               CALL DGESD2D( ICTXT, ADROWS, ADCOLS,
     $              A((LJAS-1)*LLDA+LIAS), LLDA, SOUTH, 
     $              EAST )
               DESCSA(RSRC_) = MOD(DESCSA(RSRC_) + 1, NPROW)
               DESCSA(CSRC_) = MOD(DESCSA(CSRC_) + 1, NPCOL)
               CALL DGERV2D( ICTXT, ADROWS, ADCOLS, 
     $              A((LJAS-1)*LLDA+LIAS), LLDA, NORTH,
     $              WEST )
               CALL DGESD2D( ICTXT, ADROWS, ADCOLS,
     $              D((LJDS-1)*LLDD+LIDS), LLDD, SOUTH,
     $              EAST )
               DESCSD(RSRC_) = MOD(DESCSD(RSRC_) + 1, NPROW)
               DESCSD(CSRC_) = MOD(DESCSD(CSRC_) + 1, NPCOL)
               CALL DGERV2D( ICTXT, ADROWS, ADCOLS,
     $              D((LJDS-1)*LLDD+LIDS), LLDD, NORTH,
     $              WEST )
               IF( ADEXT ) THEN
                  CALL DGESD2D( ICTXT, EXMEMAD, 1, DWORK(EXA),
     $                 EXMEMAD, SOUTH, EAST )
                  CALL DGERV2D( ICTXT, EXMEMAD, 1, DWORK(EXA),
     $                 EXMEMAD, NORTH, WEST )
                  CALL DGESD2D( ICTXT, EXMEMAD, 1, DWORK(EXD),
     $                 EXMEMAD, SOUTH,EAST )
                  CALL DGERV2D( ICTXT, EXMEMAD, 1, DWORK(EXD),
     $                 EXMEMAD, NORTH, WEST )
               END IF
            ELSEIF( TRANSAD .AND. TRANSBE ) THEN
               CALL DGESD2D( ICTXT, ADROWS, ADCOLS,
     $              A((LJAS-1)*LLDA+LIAS), LLDA, NORTH,
     $              WEST )
               DESCSA(RSRC_) = MOD(DESCSA(RSRC_) + NPROW - 1, NPROW)
               DESCSA(CSRC_) = MOD(DESCSA(CSRC_) + NPCOL - 1, NPCOL)
               CALL DGERV2D( ICTXT, ADROWS, ADCOLS, 
     $              A((LJAS-1)*LLDA+LIAS), LLDA, SOUTH,
     $              EAST )
               CALL DGESD2D( ICTXT, ADROWS, ADCOLS,
     $              D((LJDS-1)*LLDD+LIDS), LLDD, NORTH,
     $              WEST )
               DESCSD(RSRC_) = MOD(DESCSD(RSRC_) + NPROW - 1, NPROW)
               DESCSD(CSRC_) = MOD(DESCSD(CSRC_) + NPCOL - 1, NPCOL)
               CALL DGERV2D( ICTXT, ADROWS, ADCOLS,
     $              D((LJDS-1)*LLDD+LIDS), LLDD, SOUTH,
     $              EAST )
               IF( ADEXT ) THEN
                  CALL DGESD2D( ICTXT, EXMEMAD, 1, DWORK(EXA),
     $                 EXMEMAD, NORTH, WEST )
                  CALL DGERV2D( ICTXT, EXMEMAD, 1, DWORK(EXA),
     $                 EXMEMAD, SOUTH, EAST )
                  CALL DGESD2D( ICTXT, EXMEMAD, 1, DWORK(EXD),
     $                 EXMEMAD, NORTH, WEST )
                  CALL DGERV2D( ICTXT, EXMEMAD, 1, DWORK(EXD),
     $                 EXMEMAD, SOUTH, EAST )
               END IF
            END IF
            IF ( NPROW.GT.1 ) THEN
               IF( (.NOT.TRANSAD) .AND. (.NOT.TRANSBE) ) THEN
                  CALL DGESD2D( ICTXT, ADROWS, BECOLS, 
     $                 C((LJCS-1)*LLDC+LICS), LLDC, SOUTH,
     $                 MYCOL )
                  DESCSC(RSRC_) = MOD(DESCSC(RSRC_) + 1, NPROW)
                  CALL DGERV2D( ICTXT, ADROWS, BECOLS, 
     $                 C((LJCS-1)*LLDC+LICS), LLDC, NORTH,
     $                 MYCOL )
                  CALL DGESD2D( ICTXT, ADROWS, BECOLS,
     $                 F((LJFS-1)*LLDF+LIFS), LLDF, SOUTH,
     $                 MYCOL )
                  DESCSF(RSRC_) = MOD(DESCSF(RSRC_) + 1, NPROW)
                  CALL DGERV2D( ICTXT, ADROWS, BECOLS,
     $                 F((LJFS-1)*LLDF+LIFS), LLDF, NORTH,
     $                 MYCOL )
                  IF( ADEXT .OR. BEEXT ) THEN
                     CALL DGESD2D( ICTXT, EXMEMCF, 1, DWORK(EXC),
     $                    EXMEMCF, SOUTH, MYCOL )
                     CALL DGERV2D( ICTXT, EXMEMCF, 1, DWORK(EXC),
     $                    EXMEMCF, NORTH, MYCOL )
                     CALL DGESD2D( ICTXT, EXMEMCF, 1, DWORK(EXF),
     $                    EXMEMCF, SOUTH, MYCOL )
                     CALL DGERV2D( ICTXT, EXMEMCF, 1, DWORK(EXF),
     $                    EXMEMCF, NORTH, MYCOL )
                  END IF
               ELSEIF( TRANSAD .AND. TRANSBE ) THEN
                  CALL DGESD2D( ICTXT, ADROWS, BECOLS,
     $                 C((LJCS-1)*LLDC+LICS), LLDC, NORTH,
     $                 MYCOL )
                  DESCSC(RSRC_) = MOD(DESCSC(RSRC_) + NPROW-1,NPROW)
                  CALL DGERV2D( ICTXT, ADROWS, BECOLS,
     $                 C((LJCS-1)*LLDC+LICS), LLDC, SOUTH,
     $                 MYCOL )
                  CALL DGESD2D( ICTXT, ADROWS, BECOLS,
     $                 F((LJFS-1)*LLDF+LIFS), LLDF, NORTH,
     $                 MYCOL )
                  DESCSF(RSRC_) = MOD(DESCSF(RSRC_) + NPROW-1,NPROW)
                  CALL DGERV2D( ICTXT, ADROWS, BECOLS,
     $                 F((LJFS-1)*LLDF+LIFS), LLDF, SOUTH,
     $                 MYCOL )
                  IF( ADEXT .OR. BEEXT ) THEN
                     CALL DGESD2D( ICTXT, EXMEMCF, 1, DWORK(EXC),
     $                    EXMEMCF, NORTH, MYCOL )
                     CALL DGERV2D( ICTXT, EXMEMCF, 1, DWORK(EXC),
     $                    EXMEMCF, SOUTH, MYCOL )
                     CALL DGESD2D( ICTXT, EXMEMCF, 1, DWORK(EXF),
     $                    EXMEMCF, NORTH, MYCOL )
                     CALL DGERV2D( ICTXT, EXMEMCF, 1, DWORK(EXF),
     $                    EXMEMCF, SOUTH, MYCOL )
                  END IF
               END IF
            END IF
         ELSEIF (SHIFT .AND. NPROCS.GT.1) THEN 
C     
C     Shift (B,E) NorthWest and (C,F) West if NPCOL > 1 
C     or shift (B,E)**T SouthEast and (C,F) East if NPCOL > 1 
C     
            IF( TRANSZ ) THEN
               TRANSAD = .NOT. TRANSAD
               TRANSBE = .NOT. TRANSBE
            END IF
            IF( (.NOT.TRANSAD) .AND. (.NOT.TRANSBE) ) THEN
               CALL DGESD2D( ICTXT, BEROWS, BECOLS,
     $              B((LJBS-1)*LLDB+LIBS), LLDB, NORTH,
     $              WEST )
               DESCSB(RSRC_) = MOD(DESCSB(RSRC_) + NPROW - 1, NPROW)
               DESCSB(CSRC_) = MOD(DESCSB(CSRC_) + NPCOL - 1, NPCOL)
               CALL DGERV2D( ICTXT, BEROWS, BECOLS,
     $              B((LJBS-1)*LLDB+LIBS), LLDB, SOUTH,
     $              EAST )
               CALL DGESD2D( ICTXT, BEROWS, BECOLS,
     $              E((LJES-1)*LLDE+LIES), LLDE, NORTH,
     $              WEST )
               DESCSE(RSRC_) = MOD(DESCSE(RSRC_) + NPROW - 1, NPROW)
               DESCSE(CSRC_) = MOD(DESCSE(CSRC_) + NPCOL - 1, NPCOL)
               CALL DGERV2D( ICTXT, BEROWS, BECOLS,
     $              E((LJES-1)*LLDE+LIES), LLDE, SOUTH,
     $              EAST )
               IF( BEEXT ) THEN
                  CALL DGESD2D( ICTXT, EXMEMBE, 1, DWORK(EXB),
     $                 EXMEMBE, NORTH, WEST )
                  CALL DGERV2D( ICTXT, EXMEMBE, 1, DWORK(EXB),
     $                 EXMEMBE, SOUTH, EAST )
                  CALL DGESD2D( ICTXT, EXMEMBE, 1, DWORK(EXE),
     $                 EXMEMBE, NORTH, WEST )
                  CALL DGERV2D( ICTXT, EXMEMBE, 1, DWORK(EXE),
     $                 EXMEMBE, SOUTH, EAST )
               END IF
            ELSEIF( TRANSAD .AND. TRANSBE ) THEN
               CALL DGESD2D( ICTXT, BEROWS, BECOLS, 
     $              B((LJBS-1)*LLDB+LIBS), LLDB, SOUTH,
     $              EAST )
               DESCSB(RSRC_) = MOD(DESCSB(RSRC_) + 1, NPROW)
               DESCSB(CSRC_) = MOD(DESCSB(CSRC_) + 1, NPCOL)
               CALL DGERV2D( ICTXT, BEROWS, BECOLS,
     $              B((LJBS-1)*LLDB+LIBS), LLDB, NORTH, 
     $              WEST )
               CALL DGESD2D( ICTXT, BEROWS, BECOLS,
     $              E((LJES-1)*LLDE+LIES), LLDE, SOUTH,
     $              EAST )
               DESCSE(RSRC_) = MOD(DESCSE(RSRC_) + 1, NPROW)
               DESCSE(CSRC_) = MOD(DESCSE(CSRC_) + 1, NPCOL)
               CALL DGERV2D( ICTXT, BEROWS, BECOLS,
     $              E((LJES-1)*LLDE+LIES), LLDE, NORTH,
     $              WEST )
               IF( BEEXT ) THEN
                  CALL DGESD2D( ICTXT, EXMEMBE, 1, DWORK(EXB),
     $                 EXMEMBE, SOUTH, EAST )
                  CALL DGERV2D( ICTXT, EXMEMBE, 1, DWORK(EXB),
     $                 EXMEMBE, NORTH, WEST )
                  CALL DGESD2D( ICTXT, EXMEMBE, 1, DWORK(EXE),
     $                 EXMEMBE, SOUTH, EAST )
                  CALL DGERV2D( ICTXT, EXMEMBE, 1, DWORK(EXE),
     $                 EXMEMBE, NORTH, WEST )
               END IF
            END IF
            IF ( NPCOL.GT.1 ) THEN
               IF( (.NOT.TRANSAD) .AND. (.NOT.TRANSBE) ) THEN
                  CALL DGESD2D( ICTXT, ADROWS, BECOLS,
     $                 C((LJCS-1)*LLDC+LICS), LLDC, MYROW,
     $                 WEST )
                  DESCSC(CSRC_) = MOD(DESCSC(CSRC_) + NPCOL-1,NPCOL)
                  CALL DGERV2D( ICTXT, ADROWS, BECOLS,
     $                 C((LJCS-1)*LLDC+LICS), LLDC, MYROW, 
     $                 EAST )
                  CALL DGESD2D( ICTXT, ADROWS, BECOLS,
     $                 F((LJFS-1)*LLDF+LIFS), LLDF, MYROW,
     $                 WEST )
                  DESCSF(CSRC_) = MOD(DESCSF(CSRC_) + NPCOL-1,NPCOL)
                  CALL DGERV2D( ICTXT, ADROWS, BECOLS,
     $                 F((LJFS-1)*LLDF+LIFS), LLDF, MYROW, 
     $                 EAST )
                  IF( ADEXT .OR. BEEXT ) THEN
                     CALL DGESD2D( ICTXT, EXMEMCF, 1, DWORK(EXC),
     $                    EXMEMCF, MYROW, WEST)
                     CALL DGERV2D( ICTXT, EXMEMCF, 1, DWORK(EXC),
     $                    EXMEMCF, MYROW, EAST )
                     CALL DGESD2D( ICTXT, EXMEMCF, 1, DWORK(EXF),
     $                    EXMEMCF, MYROW, WEST )
                     CALL DGERV2D( ICTXT, EXMEMCF, 1, DWORK(EXF),
     $                    EXMEMCF, MYROW, EAST )
                  END IF
               ELSEIF( TRANSAD .AND. TRANSBE ) THEN
                  CALL DGESD2D( ICTXT, ADROWS, BECOLS,
     $                 C((LJCS-1)*LLDC+LICS), LLDC, MYROW,
     $                 EAST )
                  DESCSC(CSRC_) = MOD(DESCSC(CSRC_) + 1, NPCOL)
                  CALL DGERV2D( ICTXT, ADROWS, BECOLS,
     $                 C((LJCS-1)*LLDC+LICS), LLDC, MYROW, 
     $                 WEST )
                  CALL DGESD2D( ICTXT, ADROWS, BECOLS,
     $                 F((LJFS-1)*LLDF+LIFS), LLDF, MYROW,
     $                 EAST )
                  DESCSF(CSRC_) = MOD(DESCSF(CSRC_) + 1, NPCOL)
                  CALL DGERV2D( ICTXT, ADROWS, BECOLS,
     $                 F((LJFS-1)*LLDF+LIFS), LLDF, MYROW, 
     $                 WEST )
                  IF( ADEXT .OR. BEEXT ) THEN
                     CALL DGESD2D( ICTXT, EXMEMCF, 1, DWORK(EXC),
     $                    EXMEMCF, MYROW, EAST )
                     CALL DGERV2D( ICTXT, EXMEMCF, 1, DWORK(EXC),
     $                    EXMEMCF, MYROW, WEST )
                     CALL DGESD2D( ICTXT, EXMEMCF, 1, DWORK(EXF),
     $                    EXMEMCF, MYROW, EAST )
                     CALL DGERV2D( ICTXT, EXMEMCF, 1, DWORK(EXF),
     $                    EXMEMCF, MYROW, WEST )
                  END IF
               END IF
            END IF
         END IF
         IF( SHIFT .AND. TRANSZ ) THEN
            TRANSAD = .NOT. TRANSAD
            TRANSBE = .NOT. TRANSBE
         END IF
C     
C     Recompute the number of block columns of (A,D) and (B,E) held by 
C     my node  
C
         IF( SHIFT ) THEN
            ADCOLS2 = NUMROC( M+IROFFAD, MB, MYCOL, DESCSA(CSRC_), 
     $                       NPCOL )
            BECOLS2 = NUMROC( N+IROFFBE, NB, MYCOL, DESCSB(CSRC_), 
     $                       NPCOL )
            NBCAD = ICEIL( ADCOLS2, MB )
            NBCBE = ICEIL( BECOLS2, NB )
         END IF
C      
         JJS = JS        
C     
C     Solve subsystems on current diagonal in parallel
C       
         PHASES = ICEIL( IS-IEND+1, MNPDIM )
         DO 20 PHASE = 1, PHASES, 1
            IIS = IS-(PHASE-1)*MNPDIM
            IIE = MAX(IS-PHASE*MNPDIM+1,IEND)
            IF( (.NOT.TRANSAD).AND.(.NOT.TRANSBE) ) THEN
               JJS = JS - (PHASE-1)*MNPDIM
            ELSEIF( TRANSAD.AND.(.NOT.TRANSBE) ) THEN
               JJS = JS + (PHASE-1)*MNPDIM
            ELSEIF( (.NOT.TRANSAD).AND.TRANSBE ) THEN
               JJS = JS + (PHASE-1)*MNPDIM
            ELSEIF( TRANSAD.AND.TRANSBE ) THEN
               JJS = JS - (PHASE-1)*MNPDIM
            END IF
            DO 22 JDEEP = JDEEPS, JDEEPE, JDEEPU
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
C     Check if (Aii,Dii) is extended and set some variables describing 
C     (Cij,Fij) for this particular solve. To handle submatrices we
C     distinguish between I = 1 and I > 1.
C     
               IF( I.EQ.1 ) THEN
                  IF( IWORK(EXADINF).EQ.0 ) THEN
                     ROWS = MIN( MB - IROFFAD, M )
                     GSI = 1 + IROFFAD
                     EXROW = .FALSE.
                  ELSEIF( IWORK(EXADINF).EQ.1 ) THEN
                     ROWS = MB - IROFFAD + 1
                     GSI = 1 + IROFFAD
                     EXROW = .TRUE.
                  END IF 
               ELSE
                  IF( IWORK(EXADINF+(I-1)).EQ.0 ) THEN
                     ROWS = MIN(MB, (M - MB + IROFFAD) - (I - 2) * MB)
                     GSI = (I - 1) * MB + 1
                     EXROW = .FALSE.
                  ELSEIF( IWORK(EXADINF+(I-1)).EQ.1 ) THEN
                     ROWS = MB + 1
                     GSI = (I - 1) * MB + 1
                     EXROW = .TRUE.
                  ELSEIF( IWORK(EXADINF+(I-1)).EQ.2 ) THEN
                     ROWS = MIN(MB, (M - MB + IROFFAD) - 
     $                    (I - 2) * MB) - 1
                     GSI = (I - 1) * MB + 2
                     EXROW = .FALSE.
                  ELSEIF( IWORK(EXADINF+(I-1)).EQ.3 ) THEN
                     ROWS = MB
                     GSI = (I - 1) * MB + 2
                     EXROW = .TRUE.
                  END IF
               END IF
               IWORK( IROWS + IIS - I ) = ROWS
               IWORK( IGSI + IIS - I ) = GSI
               IWORK( IEXRW + IIS - I ) = ILG2NT(EXROW)
C     
C     Check if (Bjj,Ejj) is extended and set some variables describing 
C     (Cij,Fij)
C     
               IF( J.EQ.1 ) THEN
                  IF( IWORK(EXBEINF).EQ.0 ) THEN
                     COLS = MIN( NB - IROFFBE, N )
                     GSJ = 1 + IROFFBE
                     EXCOL = .FALSE.
                  ELSEIF( IWORK(EXBEINF).EQ.1 ) THEN
                     COLS = NB - IROFFBE + 1
                     GSJ = 1 + IROFFBE
                     EXCOL = .TRUE.
                  END IF 
               ELSE
                  IF( IWORK(EXBEINF+(J-1)).EQ.0 ) THEN
                     COLS = MIN(NB, (N - NB + IROFFBE) - (J - 2) * NB)
                     GSJ = (J - 1) * NB + 1
                     EXCOL = .FALSE.
                  ELSEIF( IWORK(EXBEINF+(J-1)).EQ.1 ) THEN
                     COLS = NB + 1
                     GSJ = (J - 1) * NB + 1
                     EXCOL = .TRUE.
                  ELSEIF( IWORK(EXBEINF+(J-1)).EQ.2 ) THEN
                     COLS = MIN(NB, (N - NB + IROFFBE) - 
     $                    (J - 2) * NB) - 1
                     GSJ = (J - 1) * NB + 2
                     EXCOL = .FALSE.
                  ELSEIF( IWORK(EXBEINF+(J-1)).EQ.3 ) THEN
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
     $              MYCOL, LIAD, LJAD, ADRSRC, ADCSRC )
               CALL INFOG2L( GSJ, GSJ, DESCSB, NPROW, NPCOL, MYROW, 
     $              MYCOL, LIBE, LJBE, BERSRC, BECSRC )
               CALL INFOG2L( GSI, GSJ, DESCSC, NPROW, NPCOL, MYROW, 
     $              MYCOL, IX, JX, CFRSRC, CFCSRC )
C
C     Submatrix construction is omitted most of the times for
C     deep pipelining
C
               IF( IDEEP.NE.IDEEPS .OR. JDEEP.NE.JDEEPS ) GO TO 32
C     
C     Build the extended matrix (Aii,Dii) and send it to the process 
C     holding (Cij,Fij)
C     
               IF( MYROW.EQ.ADRSRC .AND. MYCOL.EQ.ADCSRC ) THEN
                  CALL DBEXMAT( EXROW, EXROW, ROWS, ROWS, LIAD, LJAD, 
     $                 A((LJAS-1)*LLDA+LIAS), LLDA, NBCAD, MB, MB, 
     $                 DWORK( EXA ), DWORK( ADII ), MB + 1 )
                  CALL DBEXMAT( EXROW, EXROW, ROWS, ROWS, LIAD, LJAD, 
     $                 D((LJDS-1)*LLDD+LIDS), LLDD, NBCAD, MB, MB, 
     $                 DWORK( EXD ), DWORK( ADII+(MB+1)*ROWS ), 
     $                 MB + 1 )
                  IF( (ADRSRC.NE.CFRSRC).OR.(ADCSRC.NE.CFCSRC) ) THEN
                     CALL DGESD2D( ICTXT, ROWS, 2*ROWS, DWORK( ADII ),
     $                    MB + 1, CFRSRC, CFCSRC )
                  END IF
               END IF
C     
C     Build the extended matrix (Bjj,Ejj) and send it to the process 
C     holding (Cij,Fij)
C     
               IF( MYROW.EQ.BERSRC .AND. MYCOL.EQ.BECSRC ) THEN
                  CALL DBEXMAT( EXCOL, EXCOL, COLS, COLS, LIBE, LJBE, 
     $                 B((LJBS-1)*LLDB+LIBS), LLDB, NBCBE, NB, NB, 
     $                 DWORK( EXB ), DWORK( BEJJ ), NB + 1 )
                  CALL DBEXMAT( EXCOL, EXCOL, COLS, COLS, LIBE, LJBE, 
     $                 E((LJES-1)*LLDE+LIES), LLDE, NBCBE, NB, NB, 
     $                 DWORK( EXE ), DWORK( BEJJ+(NB+1)*COLS ), 
     $                 NB + 1 )
                  IF( (BERSRC.NE.CFRSRC).OR.(BECSRC.NE.CFCSRC) ) THEN
                     CALL DGESD2D( ICTXT, COLS, 2*COLS, DWORK( BEJJ ),
     $                    NB + 1, CFRSRC, CFCSRC )
                  END IF
               END IF
C
C     Skipped submatrix construction in deep pipelining? 
C
 32            CONTINUE
C     
C     Set some solution variables
C     
               SNODE = MYROW.EQ.CFRSRC.AND.MYCOL.EQ.CFCSRC
               SRSRC = CFRSRC
               SCSRC = CFCSRC
C     
C     The solving process: build the extended matrices Cij and Fij
C     
               IF( SNODE ) THEN
C     
                  IF( IDEEP.EQ.IDEEPS .AND. JDEEP.EQ.JDEEPS ) THEN
                     CALL DBEXMAT( EXROW, EXCOL, ROWS, COLS, IX, JX, 
     $                    C((LJCS-1)*LLDC+LICS), LLDC, NBCBE, MB, NB, 
     $                    DWORK( EXC ), DWORK( CIJ ), MB + 1 )
                     CALL DBEXMAT( EXROW, EXCOL, ROWS, COLS, IX, JX, 
     $                    F((LJFS-1)*LLDF+LIFS), LLDF, NBCBE, MB, NB, 
     $                    DWORK( EXF ), DWORK( FIJ ), MB + 1)   
C     
C     Receive from sends above
C     
                     IF( ADRSRC.NE.CFRSRC .OR. ADCSRC.NE.CFCSRC ) THEN
                        CALL DGERV2D( ICTXT, ROWS, 2*ROWS, DWORK(ADII),
     $                       MB+1, ADRSRC, ADCSRC )           
                     END IF
                     IF( BERSRC.NE.CFRSRC .OR. BECSRC.NE.CFCSRC ) THEN
                        CALL DGERV2D( ICTXT, COLS, 2*COLS, DWORK(BEJJ),
     $                       NB+1, BERSRC, BECSRC )
                     END IF
                  END IF  
               END IF
C     
C     Solve sub-system 
C     
               IF (SNODE) THEN
C     
                  CALL DTRGCSY( 'S', TRANZ, TRANAD, TRANBE, ISGN, IDEEP,
     $                 JDEEP, ROWS, MB2, COLS, NB2, DWORK(ADII), MB+1, 
     $                 DWORK(BEJJ), NB+1, DWORK(CIJ), MB+1, 
     $                 DWORK( ADII+(MB+1)*ROWS ), MB+1, 
     $                 DWORK( BEJJ+(NB+1)*COLS ), NB+1, DWORK(FIJ), 
     $                 MB+1, DWORK(IPW), 2*(MB2+1)*(NB2+1), IWORK(IW), 
     $                 8, XYRWS, XYCLS, XYRIND, XYCIND, SCALOC, INFO )
                  IWORK( IXYRWS + IIS - I ) = XYRWS
                  IWORK( IXYCLS + IIS - I ) = XYCLS
                  IWORK( IXYRIND + IIS - I ) = XYRIND
                  IWORK( IXYCIND + IIS - I ) = XYCIND 
C     
                  IF (LINFO.NE.0) INFO = LINFO
C     
C     Copy the given solutions back to the global matrices and update 
C     the extension elements
C     
                  IF( IDEEP.EQ.IDEEPE .AND. JDEEP.EQ.JDEEPE .AND.
     $                SCALOC.EQ.ONE ) THEN
                     CALL DUBEXMA( EXROW, EXCOL, ROWS, COLS, IX, JX, 
     $                    C((LJCS-1)*LLDC+LICS), LLDC, NBCBE, MB, NB, 
     $                    DWORK( EXC ), DWORK( CIJ ), MB + 1 )
                     CALL DUBEXMA( EXROW, EXCOL, ROWS, COLS, IX, JX, 
     $                    F((LJFS-1)*LLDF+LIFS), LLDF, NBCBE, MB, NB, 
     $                    DWORK( EXF ), DWORK( FIJ ), MB + 1 )
                  END IF
C     
C     Warn for owerflow 
C     
                  IF (SCALOC.NE.ONE) THEN
                     INFO = 3
                  END IF
C     
C     Let the solving process copy the solution into the DWORK( XIJ ) 
C     and DWORK( YIJ ) space as a preparation for the updates of the 
C     global solution
C     
                  CALL DLACPY( 'All', XYRWS, XYCLS, 
     $                 DWORK( CIJ + (XYCIND-1)*(MB+1) + XYRIND - 1), 
     $                 MB + 1,
     $                 DWORK( XIJ1 + (XYCIND-1)*(MB+1) + XYRIND - 1 ), 
     $                 MB + 1 )
                  CALL DLACPY( 'All', XYRWS, XYCLS, 
     $                 DWORK( FIJ + (XYCIND-1)*(MB+1) + XYRIND - 1), 
     $                 MB + 1,
     $                 DWORK( YIJ1 + (XYCIND-1)*(MB+1) + XYRIND - 1 ), 
     $                 MB + 1 )
                  IF( TRANSZ ) THEN
                     CALL DLACPY( 'All', XYRWS, XYCLS, 
     $                    DWORK( CIJ + (XYCIND-1)*(MB+1) + XYRIND - 1), 
     $                    MB + 1,
     $                    DWORK( XIJ2 + (XYCIND-1)*(MB+1) + XYRIND - 1), 
     $                    MB + 1 )
                     CALL DLACPY( 'All', XYRWS, XYCLS, 
     $                    DWORK( FIJ + (XYCIND-1)*(MB+1) + XYRIND - 1), 
     $                    MB + 1,
     $                    DWORK( YIJ2 + (XYCIND-1)*(MB+1) + XYRIND - 1), 
     $                    MB + 1 )
                  END IF
               END IF
 35            CONTINUE
C     
C     Update inner loop variable j (column index)
C     
               IF( .NOT. TRANSZ ) THEN
                  IF( (.NOT.TRANSAD).AND.(.NOT.TRANSBE) ) THEN
                     J = J - 1
                  ELSEIF( TRANSAD.AND.(.NOT.TRANSBE) ) THEN
                     J = J + 1
                  ELSEIF( (.NOT.TRANSAD).AND.TRANSBE ) THEN
                     J = J + 1
                  ELSEIF( TRANSAD.AND.TRANSBE ) THEN
                     J = J - 1
                  END IF
               ELSE
                  IF( (.NOT.TRANSAD).AND.(.NOT.TRANSBE) ) THEN
                     J = J - 1
                  ELSEIF( TRANSAD.AND.TRANSBE ) THEN
                     J = J - 1
                  END IF
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
               IF( .NOT. TRANSZ ) THEN
                  IF( (.NOT.TRANSAD).AND.(.NOT.TRANSBE) ) THEN
                     J = J - 1
                  ELSEIF( TRANSAD.AND.(.NOT.TRANSBE) ) THEN
                     J = J + 1
                  ELSEIF( (.NOT.TRANSAD).AND.TRANSBE ) THEN
                     J = J + 1
                  ELSEIF( TRANSAD.AND.TRANSBE ) THEN
                     J = J - 1
                  END IF
               ELSE
                  IF( (.NOT.TRANSAD).AND.(.NOT.TRANSBE) ) THEN
                     J = J - 1
                  ELSEIF( TRANSAD.AND.TRANSBE ) THEN
                     J = J - 1
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
     $                    0, MYCOL )
                  END IF
               END IF
            ELSE
               IF( NPCOL.GT.1 ) THEN
                  IF( MYCOL.EQ.0 ) THEN
                     CALL DGEBS2D( ICTXT, 'Row', ' ', 1, 1, SCALOC, 1 )
                  ELSE
                     CALL DGEBR2D( ICTXT, 'Row', ' ', 1, 1, SCALOC, 1,
     $                    MYROW, 0 )
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
                  CALL PDSCAL( M, SCALOC, F, IF, JF+KKK-1, DESCF, 1 )
 123           CONTINUE
               IF( ADEXT.OR.BEEXT ) THEN
                  CALL DSCAL( EXMEMCF, SCALOC, DWORK( EXC ), 1 )
                  CALL DSCAL( EXMEMCF, SCALOC, DWORK( EXF ), 1 )
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
                     XYRWS = IWORK( IXYRWS + IIS - I ) 
                     XYCLS = IWORK( IXYCLS + IIS - I ) 
                     XYRIND = IWORK( IXYRIND + IIS - I )
                     XYCIND = IWORK( IXYCIND + IIS - I )
                     IF( SCALOC.NE.LSCALC ) THEN
                        DO 333 KKK = 1, XYCLS
                           CALL DSCAL( XYRWS, SCALOC / LSCALC, 
     $                          DWORK( XIJ1 + (XYCIND-1)*(MB+1) + 
     $                          XYRIND - 1 + (KKK-1)*(MB+1) ), 1 )
                           IF( TRANSZ ) 
     $                          CALL DSCAL( XYRWS, SCALOC / LSCALC, 
     $                          DWORK( XIJ2 + (XYCIND-1)*(MB+1) + 
     $                          XYRIND - 1 + (KKK-1)*(MB+1) ), 1 )
 333                    CONTINUE
                        DO 444 KKK = 1, XYCLS
                           CALL DSCAL( XYRWS, SCALOC / LSCALC, 
     $                          DWORK( YIJ1 + (XYCIND-1)*(MB+1) + 
     $                          XYRIND - 1 + (KKK-1)*(MB+1) ), 1 )
                           IF( TRANSZ ) 
     $                          CALL DSCAL( XYRWS, SCALOC / LSCALC, 
     $                          DWORK( YIJ2 + (XYCIND-1)*(MB+1) + 
     $                          XYRIND - 1 + (KKK-1)*(MB+1) ), 1 )
 444                    CONTINUE
                     END IF
                     DO 555 KKK = 1, COLS
                        CALL DSCAL( ROWS, SCALOC, 
     $                       DWORK( CIJ + (KKK-1)*(MB+1) ), 1 )
                        CALL DSCAL( ROWS, SCALOC, 
     $                       DWORK( FIJ + (KKK-1)*(MB+1) ), 1 )
 555                 CONTINUE
                     CALL DLACPY( 'All', XYRWS, XYCLS, 
     $                    DWORK( XIJ1 + (XYCIND-1)*(MB+1) +
     $                    XYRIND - 1), MB+1,  DWORK( CIJ + 
     $                    (XYCIND-1)*(MB+1) + XYRIND - 1), MB+1 )
                     CALL DLACPY( 'All', XYRWS, XYCLS, 
     $                    DWORK( YIJ1 + (XYCIND-1)*(MB+1) +
     $                    XYRIND - 1), MB+1,  DWORK( FIJ + 
     $                    (XYCIND-1)*(MB+1) + XYRIND - 1), MB+1 )
                     IF( IDEEP.EQ.IDEEPE .AND. JDEEP.EQ.JDEEPE ) THEN
                        EXROW = INT2LG(IWORK( IEXRW + IIS - I )) 
                        EXCOL = INT2LG(IWORK( IEXCL + IIS - I ))
                        CALL DUBEXMA( EXROW, EXCOL, ROWS, COLS, IX, JX, 
     $                    C((LJCS-1)*LLDC+LICS), LLDC, NBCBE, MB, NB, 
     $                    DWORK( EXC ), DWORK( CIJ ), MB + 1 )
                        CALL DUBEXMA( EXROW, EXCOL, ROWS, COLS, IX, JX, 
     $                    F((LJFS-1)*LLDF+LIFS), LLDF, NBCBE, MB, NB, 
     $                    DWORK( EXF ), DWORK( FIJ ), MB + 1 )
                     END IF
                  END IF
C
 37               CONTINUE
C     
                  IF( .NOT. TRANSZ ) THEN
                     IF( (.NOT.TRANSAD).AND.(.NOT.TRANSBE) ) THEN
                        J = J - 1
                     ELSEIF( TRANSAD.AND.(.NOT.TRANSBE) ) THEN
                        J = J + 1
                     ELSEIF( (.NOT.TRANSAD).AND.TRANSBE ) THEN
                        J = J + 1
                     ELSEIF( TRANSAD.AND.TRANSBE ) THEN
                        J = J - 1
                     END IF
                  ELSE
                     IF( (.NOT.TRANSAD).AND.(.NOT.TRANSBE) ) THEN
                        J = J - 1
                     ELSEIF( TRANSAD.AND.TRANSBE ) THEN
                        J = J - 1
                     END IF
                  END IF
C     
 222           CONTINUE
C     
C     Update value of SCALE according to SCALOC
C     
               SCALE = SCALE * SCALOC
            END IF
C     
C     Broadcast the local solutions in block row i and block column j
C
            IF( K.LT.NROLL ) THEN
               DO 39 RCDIR = 1, 2
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
                     CALL INFOG2L( GSI, GSJ, DESCSC, NPROW, NPCOL, 
     $                    MYROW, MYCOL, IX, JX, CFRSRC, CFCSRC )
C     
                     SNODE = MYROW.EQ.CFRSRC.AND.MYCOL.EQ.CFCSRC
                     SRSRC = CFRSRC
                     SCSRC = CFCSRC
C     
                     IF( SNODE ) THEN  
                        XYRWS = IWORK( IXYRWS + IIS - I ) 
                        XYCLS = IWORK( IXYCLS + IIS - I ) 
                        XYRIND = IWORK( IXYRIND + IIS - I )
                        XYCIND = IWORK( IXYCIND + IIS - I )
                        IBUFF(1) = XYRWS
                        IBUFF(2) = XYCLS
                        IBUFF(3) = XYRIND
                        IBUFF(4) = XYCIND
                        
                        IF( .NOT. TRANSZ ) THEN
                           IF ( NPCOL.GT.1 .AND. RCDIR.EQ.1 ) THEN
                              CALL IGEBS2D( ICTXT, 'Row', ' ', 4, 1,
     $                             IBUFF, 4 )
                              CALL DGEBS2D( ICTXT, 'Row', ' ', XYRWS,
     $                             XYCLS, 
     $                             DWORK( YIJ1+(XYCIND-1)*(MB+1)+
     $                             XYRIND - 1 ), MB + 1 )
                           END IF
                           IF ( NPROW.GT.1 .AND. RCDIR.EQ.2 ) THEN
                              CALL IGEBS2D( ICTXT, 'Col', ' ', 4, 1,
     $                             IBUFF, 4 )
                              CALL DGEBS2D( ICTXT, 'Col', ' ', XYRWS,
     $                             XYCLS, 
     $                             DWORK( XIJ1+(XYCIND-1)*(MB+1)+
     $                             XYRIND - 1 ), MB + 1 )
                           END IF
                        ELSE
                           IF ( NPCOL.GT.1 .AND. RCDIR.EQ.1 ) THEN
                              CALL DLACPY( 'All', XYRWS, XYCLS, 
     $                             DWORK( XIJ1+(XYCIND-1)*(MB+1)+
     $                             XYRIND - 1 ), MB + 1, DWORK(IPW), 
     $                             XYRWS )
                              CALL DLACPY( 'All', XYRWS, XYCLS, 
     $                             DWORK( YIJ1+(XYCIND-1)*(MB+1)+
     $                             XYRIND - 1 ), MB + 1, 
     $                             DWORK(IPW+XYRWS*XYCLS), XYRWS )
                              CALL IGEBS2D( ICTXT, 'Row', ' ', 4, 1,
     $                             IBUFF, 4 )
                              CALL DGEBS2D( ICTXT, 'Row', ' ', XYRWS, 
     $                             2*XYCLS, DWORK( IPW ), XYRWS )
                           END IF
                           IF ( NPROW.GT.1 .AND. RCDIR.EQ.2 ) THEN
                              CALL DLACPY( 'All', XYRWS, XYCLS, 
     $                             DWORK( XIJ2+(XYCIND-1)*(MB+1)+
     $                             XYRIND - 1 ), MB + 1, DWORK(IPW), 
     $                             XYRWS )
                              CALL DLACPY( 'All', XYRWS, XYCLS, 
     $                             DWORK( YIJ2+(XYCIND-1)*(MB+1)+
     $                             XYRIND - 1 ), MB + 1, 
     $                             DWORK(IPW+XYRWS*XYCLS), XYRWS )
                              CALL IGEBS2D( ICTXT, 'Col', ' ', 4, 1,
     $                             IBUFF, 4 )
                              CALL DGEBS2D( ICTXT, 'Col', ' ', XYRWS, 
     $                             2*XYCLS, DWORK( IPW ), XYRWS )
                              
                           END IF
                        END IF
                     ELSE
C     
C     Participate in bcast of solution
C     
                        IF( .NOT. TRANSZ ) THEN
                           IF ( NPCOL.GT.1.AND.MYROW.EQ.SRSRC 
     $                          .AND. RCDIR.EQ.1 ) THEN
                              CALL IGEBR2D( ICTXT, 'Row', ' ', 4, 1,
     $                             IBUFF, 4, SRSRC, SCSRC )
                              IWORK( IXYRWS + IIS - I ) = IBUFF(1)
                              IWORK( IXYCLS + IIS - I ) = IBUFF(2)
                              IWORK( IXYRIND + IIS - I ) = IBUFF(3)
                              IWORK( IXYCIND + IIS - I ) = IBUFF(4)
                              XYRWS = IWORK( IXYRWS + IIS - I ) 
                              XYCLS = IWORK( IXYCLS + IIS - I ) 
                              XYRIND = IWORK( IXYRIND + IIS - I )
                              XYCIND = IWORK( IXYCIND + IIS - I )
                              CALL DGEBR2D( ICTXT, 'Row', ' ', XYRWS,
     $                             XYCLS, DWORK(YIJ1 +(XYCIND-1)*(MB+1)+
     $                             XYRIND - 1 ), MB+1, SRSRC, SCSRC )
                           END IF
                           IF ( NPROW.GT.1.AND.MYCOL.EQ.SCSRC 
     $                          .AND. RCDIR.EQ.2 ) THEN
                              CALL IGEBR2D( ICTXT, 'Col', ' ', 4, 1,
     $                             IBUFF, 4, SRSRC, SCSRC )
                              IWORK( IXYRWS + IIS - I ) = IBUFF(1)
                              IWORK( IXYCLS + IIS - I ) = IBUFF(2)
                              IWORK( IXYRIND + IIS - I ) = IBUFF(3)
                              IWORK( IXYCIND + IIS - I ) = IBUFF(4)
                              XYRWS = IWORK( IXYRWS + IIS - I ) 
                              XYCLS = IWORK( IXYCLS + IIS - I ) 
                              XYRIND = IWORK( IXYRIND + IIS - I )
                              XYCIND = IWORK( IXYCIND + IIS - I )
                              CALL DGEBR2D( ICTXT, 'Col', ' ', XYRWS,
     $                             XYCLS, DWORK(XIJ1+(XYCIND-1)*(MB+1)+
     $                             XYRIND - 1 ), MB+1, SRSRC, SCSRC )
                           END IF
                        ELSE
                           IF ( NPCOL.GT.1.AND.MYROW.EQ.SRSRC 
     $                          .AND. RCDIR.EQ.1 ) THEN   
                              CALL IGEBR2D( ICTXT, 'Row', ' ', 4, 1,
     $                             IBUFF, 4, SRSRC, SCSRC )
                              IWORK( IXYRWS + IIS - I ) = IBUFF(1)
                              IWORK( IXYCLS + IIS - I ) = IBUFF(2)
                              IWORK( IXYRIND + IIS - I ) = IBUFF(3)
                              IWORK( IXYCIND + IIS - I ) = IBUFF(4)
                              XYRWS = IWORK( IXYRWS + IIS - I ) 
                              XYCLS = IWORK( IXYCLS + IIS - I ) 
                              XYRIND = IWORK( IXYRIND + IIS - I )
                              XYCIND = IWORK( IXYCIND + IIS - I )
                              CALL DGEBR2D( ICTXT, 'Row', ' ', XYRWS, 
     $                             2*XYCLS, DWORK( IPW ), XYRWS, SRSRC, 
     $                             SCSRC )
                              CALL DLACPY( 'All', XYRWS, XYCLS, 
     $                             DWORK(IPW), XYRWS,
     $                             DWORK( XIJ1+(XYCIND-1)*(MB+1)+
     $                             XYRIND - 1 ), MB + 1 )
                              CALL DLACPY( 'All', XYRWS, XYCLS, 
     $                             DWORK(IPW+XYRWS*XYCLS), XYRWS,
     $                             DWORK( YIJ1+(XYCIND-1)*(MB+1)+
     $                             XYRIND - 1 ), MB + 1 )
                           END IF
                           IF ( NPROW.GT.1.AND.MYCOL.EQ.SCSRC
     $                          .AND. RCDIR.EQ.2 ) THEN
                              CALL IGEBR2D( ICTXT, 'Col', ' ', 4, 1,
     $                             IBUFF, 4, SRSRC, SCSRC ) 
                              IWORK( IXYRWS + IIS - I ) = IBUFF(1)
                              IWORK( IXYCLS + IIS - I ) = IBUFF(2)
                              IWORK( IXYRIND + IIS - I ) = IBUFF(3)
                              IWORK( IXYCIND + IIS - I ) = IBUFF(4)
                              XYRWS = IWORK( IXYRWS + IIS - I ) 
                              XYCLS = IWORK( IXYCLS + IIS - I ) 
                              XYRIND = IWORK( IXYRIND + IIS - I )
                              XYCIND = IWORK( IXYCIND + IIS - I )
                              CALL DGEBR2D( ICTXT, 'Col', ' ', XYRWS, 
     $                             2*XYCLS, DWORK( IPW ), XYRWS, SRSRC, 
     $                             SCSRC )
                              CALL DLACPY( 'All', XYRWS, XYCLS, 
     $                             DWORK(IPW), XYRWS,
     $                             DWORK( XIJ2+(XYCIND-1)*(MB+1)+
     $                             XYRIND - 1 ), MB + 1 )
                              CALL DLACPY( 'All', XYRWS, XYCLS, 
     $                             DWORK(IPW+XYRWS*XYCLS), XYRWS,
     $                             DWORK( YIJ2+(XYCIND-1)*(MB+1)+
     $                             XYRIND - 1 ), MB + 1 )
                           END IF
                        END IF
                     END IF 
 45                  CONTINUE
C     
C     Update inner loop variable j (column index)
C     
                     IF( .NOT. TRANSZ ) THEN
                        IF( (.NOT.TRANSAD).AND.(.NOT.TRANSBE) ) THEN
                           J = J - 1
                        ELSEIF( TRANSAD.AND.(.NOT.TRANSBE) ) THEN
                           J = J + 1
                        ELSEIF( (.NOT.TRANSAD).AND.TRANSBE ) THEN
                           J = J + 1
                        ELSEIF( TRANSAD.AND.TRANSBE ) THEN
                           J = J - 1
                        END IF
                     ELSE
                        IF( (.NOT.TRANSAD).AND.(.NOT.TRANSBE) ) THEN
                           J = J - 1
                        ELSEIF( TRANSAD.AND.TRANSBE ) THEN
                           J = J - 1
                        END IF
                     END IF
 40               CONTINUE
 39            CONTINUE
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
     $              MYROW, MYCOL, IX, JX, CFRSRC, CFCSRC )
C     
               SNODE = MYROW.EQ.CFRSRC.AND.MYCOL.EQ.CFCSRC
C     
               IF( SNODE ) THEN
                  XYRWS = IWORK( IXYRWS + IIS - I ) 
                  XYCLS = IWORK( IXYCLS + IIS - I ) 
                  XYRIND = IWORK( IXYRIND + IIS - I )
                  XYCIND = IWORK( IXYCIND + IIS - I )
                  CALL DTRGCSY( 'U', TRANZ, TRANAD, TRANBE, ISGN, IDEEP, 
     $                 JDEEP, ROWS, MB2, COLS, NB2, DWORK(ADII), MB+1, 
     $                 DWORK(BEJJ), NB+1, DWORK(CIJ), MB+1, 
     $                 DWORK( ADII+(MB+1)*ROWS ), MB+1, 
     $                 DWORK( BEJJ+(NB+1)*COLS ), NB+1, DWORK(FIJ), 
     $                 MB+1, DWORK(IPW), 2*(MB2+1)*(NB2+1), IWORK(IW), 
     $                 8, XYRWS, XYCLS, XYRIND, XYCIND, SCALOC, INFO )
               END IF
 48            CONTINUE
C     
C     Update inner loop variable j (column index)
C     
               IF( (.NOT.TRANSAD).AND.(.NOT.TRANSBE) ) THEN
                  J = J - 1
               ELSEIF( TRANSAD.AND.(.NOT.TRANSBE) ) THEN
                  J = J + 1
               ELSEIF( (.NOT.TRANSAD).AND.TRANSBE ) THEN
                  J = J + 1
               ELSEIF( TRANSAD.AND.TRANSBE ) THEN
                  J = J - 1
               END IF
 47         CONTINUE
C     
C     Update rest of global system wrt to current solution, that is,
C     update block column j of (C,F) in parallel with respect to Yij and 
C     after that update block row i of (C,F) in parallel with respect to
C     Xij. 
C     
C     For every update we do we also need to find out the dimensions and
C     possible extensions of the submatrices involved.
C     
C     Choose branch depending on the transp. - start with block column j,
C     that is, do one of the operations
C 
C       Ckj = Ckj - Aki * Xij
C       Fkj = Fkj - Dki * Xij
C
C     or 
C
C       Ckj = Ckj - Aik**T * Xij
C       Fkj = Fkj - Fik**T * Xij
C 
C     depending on the transposes. For TRANZ = 'T' things become similar.
C
C
            J = JJS
            DO 50 I = IIS, IIE, -1
               ROWS = IWORK( IROWS + IIS - I )
               GSI = IWORK( IGSI + IIS - I ) 
               EXROW = INT2LG(IWORK( IEXRW + IIS - I )) 
               COLS = IWORK( ICOLS + IIS - I ) 
               GSJ = IWORK( IGSJ + IIS - I )
               EXCOL = INT2LG(IWORK( IEXCL + IIS - I ))
C
               IF( ROWS.EQ.0 .OR. COLS.EQ.0 ) GO TO 55
C
               IF( TRANSZ ) GOTO 100
C     
C     (A,D) not transposed
C     
               IF( .NOT.TRANSAD ) THEN
                  IF( MOD( J, 2 ).EQ.0 ) THEN
                     INDXS = 1
                     INDXE = I - 1
                     INDXU = 1
                  ELSE
                     INDXS = I - 1
                     INDXE = 1
                     INDXU = -1
                  END IF
                  KADUP = 1
                  DO 60 INDX =  INDXS, INDXE, INDXU
                     IF( INDX.EQ.1 ) THEN
                        IF( IWORK(EXADINF).EQ.0 ) THEN
                           ADROWS2 = MB - IROFFAD
                           GSIND = 1 + IROFFAD
                           CEXROW = .FALSE.
                        ELSEIF( IWORK(EXADINF).EQ.1 ) THEN
                           ADROWS2 = MB - IROFFAD + 1
                           GSIND = 1 + IROFFAD
                           CEXROW = .TRUE.
                        END IF
                     ELSE
                        IF( IWORK(EXADINF+(INDX-1)).EQ.0 ) THEN
                           ADROWS2 = MIN(MB, (M - MB + IROFFAD) - 
     $                          (INDX - 2) * MB)
                           GSIND = (INDX - 1) * MB + 1
                           CEXROW = .FALSE.
                        ELSEIF( IWORK(EXADINF+(INDX-1)).EQ.1 ) THEN
                           ADROWS2 = MB + 1
                           GSIND = (INDX - 1) * MB + 1
                           CEXROW = .TRUE.
                        ELSEIF( IWORK(EXADINF+(INDX-1)).EQ.2 ) THEN
                           ADROWS2 = MIN(MB, (M - MB + IROFFAD) - 
     $                          (INDX - 2) * MB) - 1
                           GSIND = (INDX - 1) * MB + 2
                           CEXROW = .FALSE.
                        ELSEIF( IWORK(EXADINF+(INDX-1)).EQ.3 ) THEN
                           ADROWS2 = MB
                           GSIND = (INDX - 1) * MB + 2
                           CEXROW = .TRUE.
                        END IF
                     END IF
C     
                     CALL INFOG2L( GSIND, GSI, DESCSA, NPROW, NPCOL, 
     $                    MYROW, MYCOL, LIAD, LJAD, RSRC, CSRC )
                     CALL INFOG2L( GSIND, GSJ, DESCSC, NPROW, NPCOL, 
     $                    MYROW, MYCOL, IX, JX, RRSRC, RCSRC )
C
C     Submatrix constuction is omitted most of the times for 
C     deep pipelining
C
                     IF( IDEEP.NE.IDEEPS .OR. JDEEP.NE.JDEEPS )
     $                    GO TO 62
C
                     IF( MYROW.EQ.RSRC .AND. MYCOL.EQ.CSRC ) THEN
                        IF( (RSRC.NE.RRSRC) .OR. (CSRC.NE.RCSRC) ) THEN
                           CALL DBEXMAT( CEXROW, EXROW, ADROWS2, ROWS, 
     $                          LIAD, LJAD, A((LJAS-1)*LLDA+LIAS), LLDA, 
     $                          NBCAD, MB, MB, DWORK( EXA ), 
     $                          DWORK(ADUP+(ADUPBL-2)*(MB+1)**2), MB+1)
                           CALL DBEXMAT( CEXROW, EXROW, ADROWS2, ROWS, 
     $                          LIAD, LJAD, D((LJDS-1)*LLDD+LIDS), LLDD, 
     $                          NBCAD, MB, MB, DWORK( EXD ), 
     $                          DWORK(ADUP+(ADUPBL-2)*(MB+1)**2+
     $                          (MB+1)*ROWS ), MB + 1 )
                           CALL DGESD2D( ICTXT, ADROWS2, 2*ROWS, 
     $                          DWORK(ADUP+(ADUPBL-2)*(MB+1)**2), 
     $                          MB + 1, RRSRC,RCSRC )
                        ELSE
                           CALL DBEXMAT( CEXROW, EXROW, ADROWS2, ROWS, 
     $                          LIAD, LJAD, A((LJAS-1)*LLDA+LIAS), LLDA, 
     $                          NBCAD, MB, MB, DWORK( EXA ), 
     $                          DWORK(ADUP+(KADUP-1)*(MB+1)**2), MB+1)
                           CALL DBEXMAT( CEXROW, EXROW, ADROWS2, ROWS, 
     $                          LIAD, LJAD, D((LJDS-1)*LLDD+LIDS), LLDD, 
     $                          NBCAD, MB, MB, DWORK( EXD ), 
     $                          DWORK(ADUP+(KADUP-1)*(MB+1)**2+
     $                          (MB+1)*ROWS ), MB + 1 )
                        END IF
                     END IF
C     
C    Skipped submatrix construction?
C     
 62                  CONTINUE
C     
                     IF(MYROW.EQ.RRSRC.AND.MYCOL.EQ.RCSRC) THEN
                        CALL DBEXMAT( CEXROW, EXCOL, ADROWS2, COLS, IX, 
     $                       JX, C((LJCS-1)*LLDC+LICS), LLDC, NBCBE, 
     $                       MB, NB, DWORK( EXC ), DWORK(CKJ), 
     $                       MB + 1)
                        CALL DBEXMAT( CEXROW, EXCOL, ADROWS2, COLS, IX, 
     $                       JX, F((LJFS-1)*LLDF+LIFS), LLDF, NBCBE, 
     $                       MB, NB, DWORK( EXF ), DWORK(FKJ),
     $                       MB + 1)
                        IF( IDEEP.NE.IDEEPS .OR. JDEEP.NE.JDEEPS )
     $                       GO TO 64
                        IF( (RSRC.NE.RRSRC).OR.(CSRC.NE.RCSRC) ) THEN
                           CALL DGERV2D( ICTXT, ADROWS2, 2*ROWS, 
     $                          DWORK(ADUP+(KADUP-1)*(MB+1)**2), MB+1, 
     $                          RSRC, CSRC )
                        END IF
 64                     CONTINUE
                        XYRWS = IWORK( IXYRWS + IIS - I ) 
                        XYCLS = IWORK( IXYCLS + IIS - I ) 
                        XYRIND = IWORK( IXYRIND + IIS - I )
                        XYCIND = IWORK( IXYCIND + IIS - I ) 
                        CALL DGEMM( 'N', 'N', ADROWS2, XYCLS, XYRWS, 
     $                       -ONE, DWORK( ADUP+(KADUP-1)*(MB+1)**2 +
     $                       (XYRIND-1)*(MB+1)), MB + 1, 
     $                       DWORK(XIJ1+(XYCIND-1)*(MB+1)+XYRIND-1), 
     $                       MB + 1, ONE, DWORK(CKJ+(XYCIND-1)*(MB+1)),
     $                       MB + 1 )
                        CALL DGEMM( 'N', 'N', ADROWS2, XYCLS, XYRWS, 
     $                       -ONE, DWORK( ADUP+(KADUP-1)*(MB+1)**2 +
     $                       (XYRIND-1)*(MB+1) + (MB+1)*ROWS ), MB + 1,
     $                       DWORK(XIJ1+(XYCIND-1)*(MB+1)+XYRIND-1), 
     $                       MB + 1, ONE, DWORK(FKJ+(XYCIND-1)*(MB+1)), 
     $                       MB + 1 )
                        CALL DUBEXMA( CEXROW, EXCOL, ADROWS2, COLS, IX, 
     $                       JX, C((LJCS-1)*LLDC+LICS), LLDC, NBCBE, 
     $                       MB, NB, DWORK( EXC ), DWORK(CKJ), 
     $                       MB + 1 )
                        CALL DUBEXMA( CEXROW, EXCOL, ADROWS2, COLS, IX, 
     $                       JX, F((LJFS-1)*LLDF+LIFS), LLDF, NBCBE, 
     $                       MB, NB, DWORK( EXF ), DWORK(FKJ), 
     $                       MB + 1 )
                        IF( ADUPBL.NE.1 ) KADUP = KADUP + 2
                     END IF 
 60               CONTINUE   
C     
C     (A,D) transposed  
C     
               ELSEIF( TRANSAD ) THEN
                  IF( MOD( J, 2 ).EQ.0 ) THEN
                     INDXS = I + 1
                     INDXE = DBAD
                     INDXU = 1
                  ELSE
                     INDXS = DBAD
                     INDXE = I + 1
                     INDXU = -1
                  END IF
                  KADUP = 1
                  DO 70 INDX =  INDXS, INDXE, INDXU
                     IF( IWORK(EXADINF+(INDX-1)).EQ.0 ) THEN
                        ADCOLS2 = MIN(MB, (M - MB + IROFFAD) - 
     $                       (INDX - 2) * MB)
                        GSIND = (INDX - 1) * MB + 1
                        CEXROW = .FALSE.
                     ELSEIF( IWORK(EXADINF+(INDX-1)).EQ.1 ) THEN
                        ADCOLS2 = MB + 1
                        GSIND = (INDX - 1) * MB + 1
                        CEXROW = .TRUE.
                     ELSEIF( IWORK(EXADINF+(INDX-1)).EQ.2 ) THEN
                        ADCOLS2 = MIN(MB, (M - MB + IROFFAD) - 
     $                       (INDX - 2) * MB) - 1
                        GSIND = (INDX - 1) * MB + 2
                        CEXROW = .FALSE.
                     ELSEIF( IWORK(EXADINF+(INDX-1)).EQ.3 ) THEN
                        ADCOLS2 = MB
                        GSIND = (INDX - 1) * MB + 2
                        CEXROW = .TRUE.
                     END IF
                     CALL INFOG2L( GSI, GSIND, DESCSA, NPROW, NPCOL,
     $                    MYROW, MYCOL, LIAD, LJAD, RSRC, CSRC )
                     CALL INFOG2L( GSIND, GSJ, DESCSC, NPROW, NPCOL, 
     $                    MYROW, MYCOL, IX, JX, RRSRC, RCSRC )
C
C     Submatrix constuction is omitted most of the times for 
C     deep pipelining
C
                     IF( IDEEP.NE.IDEEPS .OR. JDEEP.NE.JDEEPS )
     $                    GO TO 72
C
                     IF( MYROW.EQ.RSRC .AND. MYCOL.EQ.CSRC ) THEN
                        IF( (RSRC.NE.RRSRC) .OR. (CSRC.NE.RCSRC) ) THEN
                           CALL DBEXMAT( EXROW, CEXROW, ROWS, ADCOLS2, 
     $                          LIAD, LJAD, A((LJAS-1)*LLDA+LIAS), LLDA, 
     $                          NBCAD, MB, MB, DWORK( EXA ), 
     $                          DWORK(ADUP+(ADUPBL-2)*(MB+1)**2), 
     $                          MB + 1 )
                           CALL DBEXMAT( EXROW, CEXROW, ROWS, ADCOLS2, 
     $                          LIAD, LJAD, D((LJDS-1)*LLDD+LIDS), LLDD, 
     $                          NBCAD, MB, MB, DWORK( EXD ), 
     $                          DWORK(ADUP+(ADUPBL-2)*(MB+1)**2  + 
     $                          (MB+1)*ADCOLS2 ), MB + 1 ) 
                           CALL DGESD2D( ICTXT, ROWS, 2*ADCOLS2, 
     $                          DWORK(ADUP+(ADUPBL-2)*(MB+1)**2), 
     $                          MB + 1, RRSRC, RCSRC )
                        ELSE
                           CALL DBEXMAT( EXROW, CEXROW, ROWS, ADCOLS2, 
     $                          LIAD, LJAD, A((LJAS-1)*LLDA+LIAS), LLDA, 
     $                          NBCAD, MB, MB, DWORK( EXA ), 
     $                          DWORK(ADUP+(KADUP-1)*(MB+1)**2), 
     $                          MB + 1 )
                           CALL DBEXMAT( EXROW, CEXROW, ROWS, ADCOLS2, 
     $                          LIAD, LJAD, D((LJDS-1)*LLDD+LIDS), LLDD, 
     $                          NBCAD, MB, MB, DWORK( EXD ), 
     $                          DWORK(ADUP+(KADUP-1)*(MB+1)**2  + 
     $                          (MB+1)*ADCOLS2 ), MB + 1 )
                        END IF
                     END IF
C     
C    Skipped submatrix construction?
C     
 72                  CONTINUE
C     
                     IF(MYROW.EQ.RRSRC .AND. MYCOL.EQ.RCSRC) THEN
                        CALL DBEXMAT( CEXROW, EXCOL, ADCOLS2, COLS, IX, 
     $                       JX, C((LJCS-1)*LLDC+LICS), LLDC, NBCBE, 
     $                       MB, NB, DWORK( EXC ), DWORK(CKJ),
     $                       MB + 1 )
                        CALL DBEXMAT( CEXROW, EXCOL, ADCOLS2, COLS, IX, 
     $                       JX, F((LJFS-1)*LLDF+LIFS), LLDF, NBCBE, 
     $                       MB, NB, DWORK( EXF ), DWORK(FKJ),
     $                       MB + 1 )
                        IF( IDEEP.NE.IDEEPS .OR. JDEEP.NE.JDEEPS )
     $                       GO TO 74
                        IF( (RSRC.NE.RRSRC).OR.(CSRC.NE.RCSRC) ) THEN
                           CALL DGERV2D( ICTXT, ROWS, 2*ADCOLS2, 
     $                          DWORK( ADUP+(KADUP-1)*(MB+1)**2 ), 
     $                          MB + 1, RSRC, CSRC )
                        END IF
 74                     CONTINUE
                        XYRWS = IWORK( IXYRWS + IIS - I ) 
                        XYCLS = IWORK( IXYCLS + IIS - I ) 
                        XYRIND = IWORK( IXYRIND + IIS - I )
                        XYCIND = IWORK( IXYCIND + IIS - I ) 
                        CALL DGEMM( 'T', 'N', ADCOLS2, XYCLS, XYRWS, 
     $                       -ONE, DWORK(ADUP+(KADUP-1)*(MB+1)**2+
     $                       XYRIND-1), MB + 1, DWORK(XIJ1+
     $                       (XYCIND-1)*(MB+1)+XYRIND-1), MB + 1, ONE, 
     $                       DWORK(CKJ+(XYCIND-1)*(MB+1)), MB + 1 )
                        CALL DGEMM( 'T', 'N', ADCOLS2, XYCLS, XYRWS, 
     $                       -ONE, DWORK( ADUP+(KADUP-1)*(MB+1)**2+
     $                       XYRIND-1+(MB+1)*ADCOLS2 ), MB + 1, 
     $                       DWORK( XIJ1+(XYCIND-1)*(MB+1)+
     $                       XYRIND-1 ), MB + 1, ONE, DWORK(FKJ+
     $                       (XYCIND-1)*(MB+1)), MB + 1 )
                        CALL DUBEXMA( CEXROW, EXCOL, ADCOLS2, COLS, IX, 
     $                       JX, C((LJCS-1)*LLDC+LICS), LLDC, NBCBE, 
     $                       MB, NB, DWORK( EXC ), DWORK(CKJ), 
     $                       MB + 1 )
                        CALL DUBEXMA( CEXROW, EXCOL, ADCOLS2, COLS, IX, 
     $                       JX, F((LJFS-1)*LLDF+LIFS), LLDF, NBCBE, 
     $                       MB, NB, DWORK( EXF ), DWORK(FKJ), 
     $                       MB + 1 )
                        IF( ADUPBL.NE.1 ) KADUP = KADUP + 2
                     END IF
 70               CONTINUE   
               END IF
C     
C     Now update block row i of (C,F), that is do one of the operations
C     
C     Cik = Cik +/- Yij * Bjk 
C     Fik = Fik +/- Yij * Ejk
C     
C     or
C     
C     Cik = Cik +/- Yij * Bkj**T
C     Fik = Fik +/- Yij * Ekj**T
C     
C     depending on the transposes.
C     
C     (B,E) not transposed
C     
               IF( .NOT.TRANSBE ) THEN
                  IF( MOD( I, 2 ).EQ.0 ) THEN
                     INDXS = J + 1
                     INDXE = DBBE
                     INDXU = 1
                  ELSE
                     INDXS = DBBE
                     INDXE = J + 1
                     INDXU = -1
                  END IF
                  KBEUP = 1
                  DO 80 INDX =  INDXS, INDXE, INDXU
                     IF( IWORK(EXBEINF+(INDX-1)).EQ.0 ) THEN
                        BECOLS2 = MIN(NB, (N - NB + IROFFBE) - 
     $                       (INDX - 2) * NB)
                        GSIND = (INDX - 1) * NB + 1
                        CEXCOL = .FALSE.
                     ELSEIF( IWORK(EXBEINF+(INDX-1)).EQ.1 ) THEN
                        BECOLS2 = NB + 1
                        GSIND = (INDX - 1) * NB + 1
                        CEXCOL = .TRUE.
                     ELSEIF( IWORK(EXBEINF+(INDX-1)).EQ.2 ) THEN
                        BECOLS2 = MIN(NB, (N - NB + IROFFBE) - 
     $                       (INDX - 2) * NB) - 1
                        GSIND = (INDX - 1) * NB + 2
                        CEXCOL = .FALSE.
                     ELSEIF( IWORK(EXBEINF+(INDX-1)).EQ.3 ) THEN
                        BECOLS2 = NB 
                        GSIND = (INDX - 1) * NB + 2
                        CEXCOL = .TRUE.
                     END IF
C     
                     CALL INFOG2L( GSJ, GSIND, DESCSB, NPROW, NPCOL,
     $                    MYROW, MYCOL, LIBE, LJBE, RSRC, CSRC )
                     CALL INFOG2L( GSI, GSIND, DESCSC, NPROW, NPCOL,
     $                    MYROW, MYCOL, IX, JX, RRSRC, RCSRC )
C
C     Submatrix constuction is omitted most of the times for 
C     deep pipelining
C
                     IF( IDEEP.NE.IDEEPS .OR. JDEEP.NE.JDEEPS )
     $                    GO TO 82      
                     IF( MYROW.EQ.RSRC .AND. MYCOL.EQ.CSRC ) THEN
                        IF( (RSRC.NE.RRSRC) .OR. (CSRC.NE.RCSRC) ) THEN
                           CALL DBEXMAT( EXCOL, CEXCOL, COLS, BECOLS2, 
     $                          LIBE, LJBE, B((LJBS-1)*LLDB+LIBS), LLDB, 
     $                          NBCBE, NB, NB, DWORK( EXB ), 
     $                          DWORK( BEUP+(BEUPBL-2)*(NB+1)**2 ), 
     $                          NB + 1 )
                           CALL DBEXMAT( EXCOL, CEXCOL, COLS, BECOLS2, 
     $                          LIBE, LJBE, E((LJES-1)*LLDE+LIES), LLDE, 
     $                          NBCBE, NB, NB, DWORK( EXE ), 
     $                          DWORK( BEUP+(BEUPBL-2)*(NB+1)**2+
     $                          (NB+1)*BECOLS2 ), NB + 1 )
                           CALL DGESD2D( ICTXT, COLS, 2*BECOLS2, 
     $                          DWORK( BEUP+(BEUPBL-2)*(NB+1)**2 ), 
     $                          NB + 1, RRSRC, RCSRC)
                        ELSE
                           CALL DBEXMAT( EXCOL, CEXCOL, COLS, BECOLS2, 
     $                          LIBE, LJBE, B((LJBS-1)*LLDB+LIBS), LLDB, 
     $                          NBCBE, NB, NB, DWORK( EXB ), 
     $                          DWORK( BEUP+(KBEUP-1)*(NB+1)**2 ), 
     $                          NB + 1 )
                           CALL DBEXMAT( EXCOL, CEXCOL, COLS, BECOLS2, 
     $                          LIBE, LJBE, E((LJES-1)*LLDE+LIES), LLDE, 
     $                          NBCBE, NB, NB, DWORK( EXE ), 
     $                          DWORK( BEUP+(KBEUP-1)*(NB+1)**2+
     $                          (NB+1)*BECOLS2 ), NB + 1 )
                        END IF
                     END IF
C     
C    Skipped submatrix construction?
C     
 82                  CONTINUE
C
                     IF(MYROW.EQ.RRSRC .AND. MYCOL.EQ.RCSRC) THEN
                        CALL DBEXMAT( EXROW, CEXCOL, ROWS, BECOLS2,
     $                       IX, JX, C((LJCS-1)*LLDC+LICS), LLDC, 
     $                       NBCBE, MB, NB, DWORK( EXC ), 
     $                       DWORK(CKJ), MB + 1 )
                        CALL DBEXMAT( EXROW, CEXCOL, ROWS, BECOLS2,
     $                       IX, JX, F((LJFS-1)*LLDF+LIFS), LLDF, 
     $                       NBCBE, MB, NB, DWORK( EXF ), 
     $                       DWORK(FKJ), MB + 1 )
                        IF( IDEEP.NE.IDEEPS .OR. JDEEP.NE.JDEEPS )
     $                       GO TO 84
                        IF( (RSRC.NE.RRSRC).OR.(CSRC.NE.RCSRC) ) THEN
                           CALL DGERV2D( ICTXT, COLS, 2*BECOLS2,
     $                          DWORK( BEUP+(KBEUP-1)*(NB+1)**2 ), 
     $                          NB + 1, RSRC, CSRC )
                        END IF
 84                     CONTINUE
                        XYRWS = IWORK( IXYRWS + IIS - I ) 
                        XYCLS = IWORK( IXYCLS + IIS - I ) 
                        XYRIND = IWORK( IXYRIND + IIS - I )
                        XYCIND = IWORK( IXYCIND + IIS - I )
                        CALL DGEMM( 'N', 'N', XYRWS, BECOLS2, XYCLS, 
     $                       USIGN, DWORK( YIJ1+(XYCIND-1)*(MB+1)+
     $                       XYRIND-1 ), MB + 1, DWORK( BEUP+
     $                       (KBEUP-1)*(NB+1)**2+XYCIND-1 ), NB + 1, 
     $                       ONE, DWORK(CKJ+XYRIND-1), MB + 1 )
                        CALL DGEMM( 'N', 'N', XYRWS, BECOLS2, XYCLS, 
     $                       USIGN, DWORK( YIJ1+(XYCIND-1)*(MB+1)+
     $                       XYRIND-1 ), MB + 1, DWORK( BEUP+
     $                       (KBEUP-1)*(NB+1)**2+XYCIND-1+
     $                       (NB+1)*BECOLS2 ), NB + 1, ONE, 
     $                       DWORK(FKJ+XYRIND-1), MB + 1 )
                        CALL DUBEXMA( EXROW, CEXCOL, ROWS, BECOLS2,
     $                       IX, JX, C((LJCS-1)*LLDC+LICS), LLDC, 
     $                       NBCBE, MB, NB, DWORK( EXC ), 
     $                       DWORK(CKJ), MB + 1 )
                        CALL DUBEXMA( EXROW, CEXCOL, ROWS, BECOLS2, IX, 
     $                       JX, F((LJFS-1)*LLDF+LIFS), LLDF, NBCBE, 
     $                       MB, NB, DWORK( EXF ), DWORK(FKJ), 
     $                       MB + 1 )
                        IF( BEUPBL.NE.1 ) KBEUP = KBEUP + 2
                     END IF 
 80               CONTINUE
C     
C     (B,E) transposed
C     
               ELSEIF( TRANSBE ) THEN
                  IF( MOD( I, 2 ).EQ.0 ) THEN
                     INDXS = 1
                     INDXE = J - 1
                     INDXU = 1
                  ELSE
                     INDXS = J - 1
                     INDXE = 1
                     INDXU = -1
                  END IF
                  KBEUP = 1
                  DO 90 INDX =  INDXS, INDXE, INDXU
                     IF( INDX.EQ.1 ) THEN
                        IF( IWORK(EXBEINF+(INDX-1)).EQ.0 ) THEN
                           BEROWS2 = NB - IROFFBE
                           GSIND = 1 + IROFFBE
                           CEXCOL = .FALSE.
                        ELSEIF( IWORK(EXBEINF+(INDX-1)).EQ.1 ) THEN
                           BEROWS2 = NB - IROFFBE + 1
                           GSIND = 1 + IROFFBE
                           CEXCOL = .TRUE.
                        END IF
                     ELSE
                        IF( IWORK(EXBEINF+(INDX-1)).EQ.0 ) THEN
                           BEROWS2 = MIN(NB, (N - NB + IROFFBE) - 
     $                          (INDX - 2) * NB)
                           GSIND = (INDX - 1) * NB + 1
                           CEXCOL = .FALSE.
                        ELSEIF( IWORK(EXBEINF+(INDX-1)).EQ.1 ) THEN
                           BEROWS2 = NB + 1
                           GSIND = (INDX - 1) * NB + 1
                           CEXCOL = .TRUE.
                        ELSEIF( IWORK(EXBEINF+(INDX-1)).EQ.2 ) THEN
                           BEROWS2 = MIN(NB, (N - NB + IROFFBE) - 
     $                          (INDX - 2) * NB) - 1
                           GSIND = (INDX - 1) * NB + 2
                           CEXCOL = .FALSE.
                        ELSEIF( IWORK(EXBEINF+(INDX-1)).EQ.3 ) THEN
                           BEROWS2 = NB 
                           GSIND = (INDX - 1) * NB + 2
                           CEXCOL = .TRUE.
                        END IF
                     END IF
C     
                     CALL INFOG2L( GSIND, GSJ, DESCSB, NPROW, NPCOL,
     $                    MYROW, MYCOL, LIBE, LJBE, RSRC, CSRC )
                     CALL INFOG2L( GSI, GSIND, DESCSC, NPROW, NPCOL, 
     $                    MYROW, MYCOL, IX, JX, RRSRC, RCSRC )
C
C     Submatrix constuction is omitted most of the times for 
C     deep pipelining
C
                     IF( IDEEP.NE.IDEEPS .OR. JDEEP.NE.JDEEPS )
     $                    GO TO 92                     
                     IF( MYROW.EQ.RSRC .AND. MYCOL.EQ.CSRC ) THEN 
                        IF( (RSRC.NE.RRSRC) .OR. (CSRC.NE.RCSRC) ) THEN
                           CALL DBEXMAT( CEXCOL, EXCOL, BEROWS2, COLS, 
     $                          LIBE, LJBE, B((LJBS-1)*LLDB+LIBS), LLDB, 
     $                          NBCBE, NB, NB, DWORK( EXB ), 
     $                          DWORK( BEUP+(BEUPBL-2)*(NB+1)**2 ), 
     $                          NB + 1 )
                           CALL DBEXMAT( CEXCOL, EXCOL, BEROWS2, COLS, 
     $                          LIBE, LJBE, E((LJES-1)*LLDE+LIES), 
     $                          LLDE, NBCBE, NB, NB, DWORK( EXE ), 
     $                          DWORK( BEUP+(BEUPBL-2)*(NB+1)**2+
     $                          (NB+1)*COLS ), NB + 1 )
                           CALL DGESD2D( ICTXT, BEROWS2, 2*COLS, 
     $                          DWORK( BEUP+(BEUPBL-2)*(NB+1)**2 ), 
     $                          NB + 1, RRSRC, RCSRC )
                        ELSE
                           CALL DBEXMAT( CEXCOL, EXCOL, BEROWS2, COLS, 
     $                          LIBE, LJBE, B((LJBS-1)*LLDB+LIBS), LLDB, 
     $                          NBCBE, NB, NB, DWORK( EXB ), 
     $                          DWORK( BEUP+(KBEUP-1)*(NB+1)**2 ), 
     $                          NB + 1 )
                           CALL DBEXMAT( CEXCOL, EXCOL, BEROWS2, COLS, 
     $                          LIBE, LJBE, E((LJES-1)*LLDE+LIES), 
     $                          LLDE, NBCBE, NB, NB, DWORK( EXE ), 
     $                          DWORK( BEUP+(KBEUP-1)*(NB+1)**2+
     $                          (NB+1)*COLS ), NB + 1 )
                        END IF
                     END IF
C     
C    Skipped submatrix construction?
C     
 92                  CONTINUE
C 
                     IF(MYROW.EQ.RRSRC .AND. MYCOL.EQ.RCSRC) THEN
                        CALL DBEXMAT( EXROW, CEXCOL, ROWS, BEROWS2, IX, 
     $                       JX, C((LJCS-1)*LLDC+LICS), LLDC, NBCBE, 
     $                       MB, NB, DWORK( EXC ), DWORK(CKJ),
     $                       MB + 1 )
                        CALL DBEXMAT( EXROW, CEXCOL, ROWS, BEROWS2, IX, 
     $                       JX, F((LJFS-1)*LLDF+LIFS), LLDF, NBCBE, 
     $                       MB, NB, DWORK( EXF ), DWORK(FKJ),
     $                       MB + 1 )
                        IF( IDEEP.NE.IDEEPS .OR. JDEEP.NE.JDEEPS )
     $                       GO TO 94
                        IF( (RSRC.NE.RRSRC) .OR. (CSRC.NE.RCSRC) ) THEN
                           CALL DGERV2D( ICTXT, BEROWS2, 2*COLS, 
     $                          DWORK( BEUP+(KBEUP-1)*(NB+1)**2 ), 
     $                          NB + 1, RSRC, CSRC )
                        END IF
 94                     CONTINUE
                        XYRWS = IWORK( IXYRWS + IIS - I ) 
                        XYCLS = IWORK( IXYCLS + IIS - I ) 
                        XYRIND = IWORK( IXYRIND + IIS - I )
                        XYCIND = IWORK( IXYCIND + IIS - I )
                        CALL DGEMM( 'N', 'T', XYRWS, BEROWS2, XYCLS, 
     $                       USIGN, DWORK( YIJ1+(XYCIND-1)*(MB+1)+
     $                       XYRIND-1), MB + 1, 
     $                       DWORK( BEUP+(KBEUP-1)*(NB+1)**2+
     $                       (XYCIND-1)*(NB+1)), NB + 1, ONE, 
     $                       DWORK( CKJ+XYRIND-1 ), MB + 1 )
                        CALL DGEMM( 'N', 'T', XYRWS, BEROWS2, XYCLS, 
     $                       USIGN, DWORK( YIJ1+(XYCIND-1)*(MB+1)+
     $                       XYRIND-1 ), MB + 1, 
     $                       DWORK( BEUP+(KBEUP-1)*(NB+1)**2+
     $                       (XYCIND-1)*(NB+1)+(NB+1)*COLS ), NB + 1,
     $                       ONE, DWORK( FKJ+XYRIND-1 ), MB + 1 )
                        CALL DUBEXMA( EXROW, CEXCOL, ROWS, BEROWS2, IX, 
     $                       JX, C((LJCS-1)*LLDC+LICS), LLDC, NBCBE, 
     $                       MB, NB, DWORK( EXC ), DWORK(CKJ), 
     $                       MB + 1 )
                        CALL DUBEXMA( EXROW, CEXCOL, ROWS, BEROWS2, IX, 
     $                       JX, F((LJFS-1)*LLDF+LIFS), LLDF, NBCBE, 
     $                       MB, NB, DWORK( EXF ), DWORK(FKJ), 
     $                       MB + 1 )
                        IF( BEUPBL.NE.1 ) KBEUP = KBEUP + 2
                     END IF
 90               CONTINUE
               END IF
C     
               GOTO 55            
C     
C     Updates connected to TRANZ = 'T'
C         
 100           CONTINUE
C     
C     (A,D) transposed
C     
               IF( TRANSAD ) THEN
                  IF( MOD( J, 2 ).EQ.0 ) THEN
                     INDXS = 1
                     INDXE = I - 1
                     INDXU = 1
                  ELSE
                     INDXS = I - 1
                     INDXE = 1
                     INDXU = -1
                  END IF
                  KADUP = 1
                  DO 110 INDX =  INDXS, INDXE, INDXU
                     IF( INDX.EQ.1 ) THEN
                        IF( IWORK(EXADINF).EQ.0 ) THEN
                           ADROWS2 = MB - IROFFAD
                           GSIND = 1 + IROFFAD
                           CEXROW = .FALSE.
                        ELSEIF( IWORK(EXADINF).EQ.1 ) THEN
                           ADROWS2 = MB - IROFFAD + 1
                           GSIND = 1 + IROFFAD
                           CEXROW = .TRUE.
                        END IF
                     ELSE
                        IF( IWORK(EXADINF+(INDX-1)).EQ.0 ) THEN
                           ADROWS2 = MIN(MB, (M - MB + IROFFAD) - 
     $                          (INDX - 2) * MB)
                           GSIND = (INDX - 1) * MB + 1
                           CEXROW = .FALSE.
                        ELSEIF( IWORK(EXADINF+(INDX-1)).EQ.1 ) THEN
                           ADROWS2 = MB + 1
                           GSIND = (INDX - 1) * MB + 1
                           CEXROW = .TRUE.
                        ELSEIF( IWORK(EXADINF+(INDX-1)).EQ.2 ) THEN
                           ADROWS2 = MIN(MB, (M - MB + IROFFAD) - 
     $                          (INDX - 2) * MB) - 1
                           GSIND = (INDX - 1) * MB + 2
                           CEXROW = .FALSE.
                        ELSEIF( IWORK(EXADINF+(INDX-1)).EQ.3 ) THEN
                           ADROWS2 = MB
                           GSIND = (INDX - 1) * MB + 2
                           CEXROW = .TRUE.
                        END IF
                     END IF
C      
                     CALL INFOG2L( GSIND, GSI, DESCSA, NPROW, NPCOL, 
     $                    MYROW, MYCOL, LIAD, LJAD, RSRC, CSRC )
                     CALL INFOG2L( GSIND, GSJ, DESCSC, NPROW, NPCOL, 
     $                    MYROW, MYCOL, IX, JX, RRSRC, RCSRC )
C
C     Submatrix constuction is omitted most of the times for 
C     deep pipelining
C
                     IF( IDEEP.NE.IDEEPS .OR. JDEEP.NE.JDEEPS )
     $                    GO TO 112
C
                     IF( MYROW.EQ.RSRC .AND. MYCOL.EQ.CSRC ) THEN
                        IF( (RSRC.NE.RRSRC) .OR. (CSRC.NE.RCSRC) ) THEN
                           CALL DBEXMAT( CEXROW, EXROW, ADROWS2, ROWS, 
     $                          LIAD, LJAD, A((LJAS-1)*LLDA+LIAS), LLDA, 
     $                          NBCAD, MB, MB, DWORK( EXA ), 
     $                          DWORK(ADUP+(ADUPBL-2)*(MB+1)**2), MB+1)
                           CALL DBEXMAT( CEXROW, EXROW, ADROWS2, ROWS, 
     $                          LIAD, LJAD, D((LJDS-1)*LLDD+LIDS), LLDD, 
     $                          NBCAD, MB, MB, DWORK( EXD ), 
     $                          DWORK(ADUP+(ADUPBL-2)*(MB+1)**2+
     $                          (MB+1)*ROWS ), MB + 1 )
                           CALL DGESD2D( ICTXT, ADROWS2, 2*ROWS, 
     $                          DWORK(ADUP+(ADUPBL-2)*(MB+1)**2), 
     $                          MB + 1, RRSRC,RCSRC )
                        ELSE
                           CALL DBEXMAT( CEXROW, EXROW, ADROWS2, ROWS, 
     $                          LIAD, LJAD, A((LJAS-1)*LLDA+LIAS), LLDA, 
     $                          NBCAD, MB, MB, DWORK( EXA ), 
     $                          DWORK(ADUP+(KADUP-1)*(MB+1)**2), MB+1)
                           CALL DBEXMAT( CEXROW, EXROW, ADROWS2, ROWS, 
     $                          LIAD, LJAD, D((LJDS-1)*LLDD+LIDS), LLDD, 
     $                          NBCAD, MB, MB, DWORK( EXD ), 
     $                          DWORK(ADUP+(KADUP-1)*(MB+1)**2+
     $                          (MB+1)*ROWS ), MB + 1 )
                        END IF  
                     END IF
C     
C    Skipped submatrix construction?
C     
 112                 CONTINUE
C     
                     IF(MYROW.EQ.RRSRC.AND.MYCOL.EQ.RCSRC) THEN
                        CALL DBEXMAT( CEXROW, EXCOL, ADROWS2, COLS, IX, 
     $                       JX, C((LJCS-1)*LLDC+LICS), LLDC, NBCBE, 
     $                       MB, NB, DWORK( EXC ), DWORK(CKJ), 
     $                       MB + 1)
                        IF( IDEEP.NE.IDEEPS .OR. JDEEP.NE.JDEEPS )
     $                       GO TO 114
                        IF( (RSRC.NE.RRSRC).OR.(CSRC.NE.RCSRC) ) THEN
                           CALL DGERV2D( ICTXT, ADROWS2, 2*ROWS, 
     $                          DWORK(ADUP+(KADUP-1)*(MB+1)**2), MB+1, 
     $                          RSRC, CSRC )
                        END IF
 114                    CONTINUE
                        XYRWS = IWORK( IXYRWS + IIS - I ) 
                        XYCLS = IWORK( IXYCLS + IIS - I ) 
                        XYRIND = IWORK( IXYRIND + IIS - I )
                        XYCIND = IWORK( IXYCIND + IIS - I ) 
                        CALL DGEMM( 'N', 'N', ADROWS2, XYCLS, XYRWS, 
     $                       -ONE, DWORK( ADUP+(KADUP-1)*(MB+1)**2 +
     $                       (XYRIND-1)*(MB+1) ), MB + 1, 
     $                       DWORK(XIJ2+(XYCIND-1)*(MB+1)+XYRIND-1), 
     $                       MB + 1, ONE, DWORK(CKJ+(XYCIND-1)*(MB+1)), 
     $                       MB + 1 )
                        CALL DGEMM( 'N', 'N', ADROWS2, XYCLS, XYRWS, 
     $                       -ONE, DWORK( ADUP+(KADUP-1)*(MB+1)**2 +
     $                       (XYRIND-1)*(MB+1)+(MB+1)*ROWS ), MB + 1,
     $                       DWORK(YIJ2+(XYCIND-1)*(MB+1)+XYRIND-1), 
     $                       MB + 1, ONE, DWORK(CKJ+(XYCIND-1)*(MB+1)), 
     $                       MB + 1 )
                        CALL DUBEXMA( CEXROW, EXCOL, ADROWS2, COLS, IX, 
     $                       JX, C((LJCS-1)*LLDC+LICS), LLDC, NBCBE, 
     $                       MB, NB, DWORK( EXC ), DWORK(CKJ), 
     $                       MB + 1 )
                        IF( ADUPBL.NE.1 ) KADUP = KADUP + 2
                     END IF 
 110              CONTINUE   
C     
C     (A,D) not transposed  
C     
               ELSEIF( .NOT. TRANSAD ) THEN
                  IF( MOD( J, 2 ).EQ.0 ) THEN
                     INDXS = I + 1
                     INDXE = DBAD
                     INDXU = 1
                  ELSE
                     INDXS = DBAD
                     INDXE = I + 1
                     INDXU = -1
                  END IF
                  KADUP = 1
                  DO 120 INDX =  INDXS, INDXE, INDXU
                     IF( IWORK(EXADINF+(INDX-1)).EQ.0 ) THEN
                        ADCOLS2 = MIN(MB, (M - MB + IROFFAD) - 
     $                       (INDX - 2) * MB)
                        GSIND = (INDX - 1) * MB + 1
                        CEXROW = .FALSE.
                     ELSEIF( IWORK(EXADINF+(INDX-1)).EQ.1 ) THEN
                        ADCOLS2 = MB + 1
                        GSIND = (INDX - 1) * MB + 1
                        CEXROW = .TRUE.
                     ELSEIF( IWORK(EXADINF+(INDX-1)).EQ.2 ) THEN
                        ADCOLS2 = MIN(MB, (M - MB + IROFFAD) - 
     $                       (INDX - 2) * MB) - 1
                        GSIND = (INDX - 1) * MB + 2
                        CEXROW = .FALSE.
                     ELSEIF( IWORK(EXADINF+(INDX-1)).EQ.3 ) THEN
                        ADCOLS2 = MB
                        GSIND = (INDX - 1) * MB + 2
                        CEXROW = .TRUE.
                     END IF
                     CALL INFOG2L( GSI, GSIND, DESCSA, NPROW, NPCOL,
     $                    MYROW, MYCOL, LIAD, LJAD, RSRC, CSRC )
                     CALL INFOG2L( GSIND, GSJ, DESCSC, NPROW, NPCOL, 
     $                    MYROW, MYCOL, IX, JX, RRSRC, RCSRC )
C
C     Submatrix constuction is omitted most of the times for 
C     deep pipelining
C
                     IF( IDEEP.NE.IDEEPS .OR. JDEEP.NE.JDEEPS )
     $                    GO TO 122
C
                     IF( MYROW.EQ.RSRC .AND. MYCOL.EQ.CSRC ) THEN
                        IF( (RSRC.NE.RRSRC) .OR. (CSRC.NE.RCSRC) ) THEN
                           CALL DBEXMAT( EXROW, CEXROW, ROWS, ADCOLS2, 
     $                          LIAD, LJAD, A((LJAS-1)*LLDA+LIAS), LLDA, 
     $                          NBCAD, MB, MB, DWORK( EXA ), 
     $                          DWORK(ADUP+(ADUPBL-2)*(MB+1)**2), 
     $                          MB + 1 )
                           CALL DBEXMAT( EXROW, CEXROW, ROWS, ADCOLS2, 
     $                          LIAD, LJAD, D((LJDS-1)*LLDD+LIDS), LLDD, 
     $                          NBCAD, MB, MB, DWORK( EXD ), 
     $                          DWORK(ADUP+(ADUPBL-2)*(MB+1)**2  + 
     $                          (MB+1)*ADCOLS2 ), MB + 1 )
                           CALL DGESD2D( ICTXT, ROWS, 2*ADCOLS2, 
     $                          DWORK(ADUP+(ADUPBL-2)*(MB+1)**2), 
     $                          MB + 1, RRSRC, RCSRC )
                        ELSE
                           CALL DBEXMAT( EXROW, CEXROW, ROWS, ADCOLS2, 
     $                          LIAD, LJAD, A((LJAS-1)*LLDA+LIAS), LLDA, 
     $                          NBCAD, MB, MB, DWORK( EXA ), 
     $                          DWORK(ADUP+(KADUP-1)*(MB+1)**2), 
     $                          MB + 1 )
                           CALL DBEXMAT( EXROW, CEXROW, ROWS, ADCOLS2, 
     $                          LIAD, LJAD, D((LJDS-1)*LLDD+LIDS), LLDD, 
     $                          NBCAD, MB, MB, DWORK( EXD ), 
     $                          DWORK(ADUP+(KADUP-1)*(MB+1)**2  + 
     $                          (MB+1)*ADCOLS2 ), MB + 1 )
                        END IF  
                     END IF
C     
C    Skipped submatrix construction?
C     
 122                  CONTINUE
C     
                     IF(MYROW.EQ.RRSRC .AND. MYCOL.EQ.RCSRC) THEN
                        CALL DBEXMAT( CEXROW, EXCOL, ADCOLS2, COLS, IX, 
     $                       JX, C((LJCS-1)*LLDC+LICS), LLDC, NBCBE, 
     $                       MB, NB, DWORK( EXC ), DWORK(CKJ),
     $                       MB + 1 )
                        IF( IDEEP.NE.IDEEPS .OR. JDEEP.NE.JDEEPS )
     $                       GO TO 124
                        IF( (RSRC.NE.RRSRC).OR.(CSRC.NE.RCSRC) ) THEN
                           CALL DGERV2D( ICTXT, ROWS, 2*ADCOLS2, 
     $                          DWORK( ADUP+(KADUP-1)*(MB+1)**2 ), 
     $                          MB + 1, RSRC, CSRC )
                        END IF
 124                    CONTINUE
                        XYRWS = IWORK( IXYRWS + IIS - I ) 
                        XYCLS = IWORK( IXYCLS + IIS - I ) 
                        XYRIND = IWORK( IXYRIND + IIS - I )
                        XYCIND = IWORK( IXYCIND + IIS - I ) 
                        CALL DGEMM( 'T', 'N', ADCOLS2, XYCLS, XYRWS, 
     $                       -ONE, DWORK( ADUP+(KADUP-1)*(MB+1)**2+
     $                       XYRIND-1 ), MB + 1, DWORK( XIJ2+
     $                       (XYCIND-1)*(MB+1)+XYRIND-1 ), MB + 1, ONE, 
     $                       DWORK(CKJ+(XYCIND-1)*(MB+1)), MB + 1 )
                        CALL DGEMM( 'T', 'N', ADCOLS2, XYCLS, XYRWS, 
     $                       -ONE, DWORK( ADUP+(KADUP-1)*(MB+1)**2+
     $                       XYRIND-1+(MB+1)*ADCOLS2 ), MB + 1, 
     $                       DWORK( YIJ2+(XYCIND-1)*(MB+1)+
     $                       XYRIND-1 ), MB + 1, ONE, DWORK(CKJ+
     $                       (XYCIND-1)*(MB+1)), MB + 1 )
                        CALL DUBEXMA( CEXROW, EXCOL, ADCOLS2, COLS, IX, 
     $                       JX, C((LJCS-1)*LLDC+LICS), LLDC, NBCBE, 
     $                       MB, NB, DWORK( EXC ), DWORK(CKJ), 
     $                       MB + 1 )
                        IF( ADUPBL.NE.1 ) KADUP = KADUP + 2
                     END IF
 120              CONTINUE   
               END IF
C     
C     (B,E) transposed
C     
               IF( TRANSBE ) THEN
                  IF( MOD( I, 2 ).EQ.0 ) THEN
                     INDXS = J + 1
                     INDXE = DBBE
                     INDXU = 1
                  ELSE
                     INDXS = DBBE
                     INDXE = J + 1
                     INDXU = -1
                  END IF
                  KBEUP = 1
                  DO 130 INDX =  INDXS, INDXE, INDXU
                     IF( IWORK(EXBEINF+(INDX-1)).EQ.0 ) THEN
                        BECOLS2 = MIN(NB, (N - NB + IROFFBE) - 
     $                       (INDX - 2) * NB)
                        GSIND = (INDX - 1) * NB + 1
                        CEXCOL = .FALSE.
                     ELSEIF( IWORK(EXBEINF+(INDX-1)).EQ.1 ) THEN
                        BECOLS2 = NB + 1
                        GSIND = (INDX - 1) * NB + 1
                        CEXCOL = .TRUE.
                     ELSEIF( IWORK(EXBEINF+(INDX-1)).EQ.2 ) THEN
                        BECOLS2 = MIN(NB, (N - NB + IROFFBE) - 
     $                       (INDX - 2) * NB) - 1
                        GSIND = (INDX - 1) * NB + 2
                        CEXCOL = .FALSE.
                     ELSEIF( IWORK(EXBEINF+(INDX-1)).EQ.3 ) THEN
                        BECOLS2 = NB 
                        GSIND = (INDX - 1) * NB + 2
                        CEXCOL = .TRUE.
                     END IF 
C
                     CALL INFOG2L( GSJ, GSIND, DESCSB, NPROW, NPCOL,
     $                    MYROW, MYCOL, LIBE, LJBE, RSRC, CSRC )
                     CALL INFOG2L( GSI, GSIND, DESCSC, NPROW, NPCOL,
     $                    MYROW, MYCOL, IX, JX, RRSRC, RCSRC )
C     
C     Submatrix constuction is omitted most of the times for 
C     deep pipelining
C
                     IF( IDEEP.NE.IDEEPS .OR. JDEEP.NE.JDEEPS )
     $                    GO TO 132  
C
                     IF( MYROW.EQ.RSRC .AND. MYCOL.EQ.CSRC ) THEN
                        IF( (RSRC.NE.RRSRC) .OR. (CSRC.NE.RCSRC) ) THEN
                           CALL DBEXMAT( EXCOL, CEXCOL, COLS, BECOLS2, 
     $                          LIBE, LJBE, B((LJBS-1)*LLDB+LIBS), LLDB, 
     $                          NBCBE, NB, NB, DWORK( EXB ), 
     $                          DWORK( BEUP+(BEUPBL-2)*(NB+1)**2 ), 
     $                          NB + 1 )
                           CALL DBEXMAT( EXCOL, CEXCOL, COLS, BECOLS2, 
     $                          LIBE, LJBE, E((LJES-1)*LLDE+LIES), LLDE, 
     $                          NBCBE, NB, NB, DWORK( EXE ), 
     $                          DWORK( BEUP+(BEUPBL-2)*(NB+1)**2+
     $                          (NB+1)*BECOLS2 ), NB + 1 )
                           CALL DGESD2D( ICTXT, COLS, 2*BECOLS2, 
     $                          DWORK( BEUP+(BEUPBL-2)*(NB+1)**2 ), 
     $                          NB + 1, RRSRC, RCSRC)
                        ELSE
                           CALL DBEXMAT( EXCOL, CEXCOL, COLS, BECOLS2, 
     $                          LIBE, LJBE, B((LJBS-1)*LLDB+LIBS), LLDB, 
     $                          NBCBE, NB, NB, DWORK( EXB ), 
     $                          DWORK( BEUP+(KBEUP-1)*(NB+1)**2 ), 
     $                          NB + 1 )
                           CALL DBEXMAT( EXCOL, CEXCOL, COLS, BECOLS2, 
     $                          LIBE, LJBE, E((LJES-1)*LLDE+LIES), LLDE, 
     $                          NBCBE, NB, NB, DWORK( EXE ), 
     $                          DWORK( BEUP+(KBEUP-1)*(NB+1)**2+
     $                          (NB+1)*BECOLS2 ), NB + 1 )
                        END IF
                     END IF
C     
C    Skipped submatrix construction?
C     
 132                  CONTINUE
C
                     IF(MYROW.EQ.RRSRC .AND. MYCOL.EQ.RCSRC) THEN
                        CALL DBEXMAT( EXROW, CEXCOL, ROWS, BECOLS2,
     $                       IX, JX, F((LJFS-1)*LLDF+LIFS), LLDF, 
     $                       NBCBE, MB, NB, DWORK( EXF ), 
     $                       DWORK(FKJ), MB + 1 )
                        IF( IDEEP.NE.IDEEPS .OR. JDEEP.NE.JDEEPS )
     $                       GO TO 134
                        IF( (RSRC.NE.RRSRC).OR.(CSRC.NE.RCSRC) ) THEN
                           CALL DGERV2D( ICTXT, COLS, 2*BECOLS2,
     $                          DWORK( BEUP+(KBEUP-1)*(NB+1)**2 ), 
     $                          NB + 1, RSRC, CSRC )
                        END IF
 134                    CONTINUE
                        XYRWS = IWORK( IXYRWS + IIS - I ) 
                        XYCLS = IWORK( IXYCLS + IIS - I ) 
                        XYRIND = IWORK( IXYRIND + IIS - I )
                        XYCIND = IWORK( IXYCIND + IIS - I )
                        CALL DGEMM( 'N', 'N', XYRWS, BECOLS2, XYCLS, 
     $                       USIGN, DWORK( XIJ1 +(XYCIND-1)*(MB+1)+
     $                       XYRIND-1 ), MB + 1, DWORK( BEUP+
     $                       (KBEUP-1)*(NB+1)**2+XYCIND-1 ), NB + 1, 
     $                       ONE, DWORK(FKJ+XYRIND-1), MB + 1 )
                        CALL DGEMM( 'N', 'N', XYRWS, BECOLS2, XYCLS, 
     $                       USIGN, DWORK( YIJ1+(XYCIND-1)*(MB+1)+
     $                       XYRIND-1 ), MB + 1, DWORK( BEUP+
     $                       (KBEUP-1)*(NB+1)**2+XYCIND-1+
     $                       (NB+1)*BECOLS2 ), NB + 1, ONE, 
     $                       DWORK(FKJ+XYRIND-1), MB + 1 )
                        CALL DUBEXMA( EXROW, CEXCOL, ROWS, BECOLS2, IX, 
     $                       JX, F((LJFS-1)*LLDF+LIFS), LLDF, NBCBE, 
     $                       MB, NB, DWORK( EXF ), DWORK(FKJ), 
     $                       MB + 1 )
                        IF( BEUPBL.NE.1 ) KBEUP = KBEUP + 2
                     END IF 
 130              CONTINUE
C     
C     (B,E) not transposed
C     
               ELSEIF( .NOT. TRANSBE ) THEN
                  IF( MOD( I, 2 ).EQ.0 ) THEN
                     INDXS = 1
                     INDXE = J - 1
                     INDXU = 1
                  ELSE
                     INDXS = J - 1
                     INDXE = 1
                     INDXU = -1
                  END IF
                  KBEUP = 1
                  DO 140 INDX =  INDXS, INDXE, INDXU
                     IF( INDX.EQ.1 ) THEN
                        IF( IWORK(EXBEINF+(INDX-1)).EQ.0 ) THEN
                           BEROWS2 = NB - IROFFBE
                           GSIND = 1 + IROFFBE
                           CEXCOL = .FALSE.
                        ELSEIF( IWORK(EXBEINF+(INDX-1)).EQ.1 ) THEN
                           BEROWS2 = NB - IROFFBE + 1
                           GSIND = 1 + IROFFBE
                           CEXCOL = .TRUE.
                        END IF
                     ELSE
                        IF( IWORK(EXBEINF+(INDX-1)).EQ.0 ) THEN
                           BEROWS2 = MIN(NB, (N - NB + IROFFBE) - 
     $                          (INDX - 2) * NB)
                           GSIND = (INDX - 1) * NB + 1
                           CEXCOL = .FALSE.
                        ELSEIF( IWORK(EXBEINF+(INDX-1)).EQ.1 ) THEN
                           BEROWS2 = NB + 1
                           GSIND = (INDX - 1) * NB + 1
                           CEXCOL = .TRUE.
                        ELSEIF( IWORK(EXBEINF+(INDX-1)).EQ.2 ) THEN
                           BEROWS2 = MIN(NB, (N - NB + IROFFBE) - 
     $                          (INDX - 2) * NB) - 1
                           GSIND = (INDX - 1) * NB + 2
                           CEXCOL = .FALSE.
                        ELSEIF( IWORK(EXBEINF+(INDX-1)).EQ.3 ) THEN
                           BEROWS2 = NB 
                           GSIND = (INDX - 1) * NB + 2
                           CEXCOL = .TRUE.
                        END IF
                     END IF
C     
                     CALL INFOG2L( GSIND, GSJ, DESCSB, NPROW, NPCOL,
     $                    MYROW, MYCOL, LIBE, LJBE, RSRC, CSRC )
                     CALL INFOG2L( GSI, GSIND, DESCSC, NPROW, NPCOL, 
     $                    MYROW, MYCOL, IX, JX, RRSRC, RCSRC )
C
C     Submatrix constuction is omitted most of the times for 
C     deep pipelining
C
                     IF( IDEEP.NE.IDEEPS .OR. JDEEP.NE.JDEEPS )
     $                    GO TO 142  
C
                     IF( MYROW.EQ.RSRC .AND. MYCOL.EQ.CSRC ) THEN
                        IF( (RSRC.NE.RRSRC) .OR. (CSRC.NE.RCSRC) ) THEN
                           CALL DBEXMAT( CEXCOL, EXCOL, BEROWS2, COLS, 
     $                          LIBE, LJBE, B((LJBS-1)*LLDB+LIBS), LLDB, 
     $                          NBCBE, NB, NB, DWORK( EXB ), 
     $                          DWORK( BEUP+(BEUPBL-2)*(NB+1)**2 ), 
     $                          NB + 1 )
                           CALL DBEXMAT( CEXCOL, EXCOL, BEROWS2, COLS, 
     $                          LIBE, LJBE, E((LJES-1)*LLDE+LIES), 
     $                          LLDE, NBCBE, NB, NB, DWORK( EXE ), 
     $                          DWORK( BEUP+(BEUPBL-2)*(NB+1)**2+
     $                          (NB+1)*COLS ), NB + 1 )
                           CALL DGESD2D( ICTXT, BEROWS2, 2*COLS, 
     $                          DWORK( BEUP+(BEUPBL-2)*(NB+1)**2 ), 
     $                          NB + 1, RRSRC, RCSRC )
                        ELSE
                           CALL DBEXMAT( CEXCOL, EXCOL, BEROWS2, COLS, 
     $                          LIBE, LJBE, B((LJBS-1)*LLDB+LIBS), LLDB, 
     $                          NBCBE, NB, NB, DWORK( EXB ), 
     $                          DWORK( BEUP+(KBEUP-1)*(NB+1)**2 ), 
     $                          NB + 1 )
                           CALL DBEXMAT( CEXCOL, EXCOL, BEROWS2, COLS, 
     $                          LIBE, LJBE, E((LJES-1)*LLDE+LIES), 
     $                          LLDE, NBCBE, NB, NB, DWORK( EXE ), 
     $                          DWORK( BEUP+(KBEUP-1)*(NB+1)**2+
     $                          (NB+1)*COLS ), NB + 1 )
                        END IF
                     END IF
C     
C    Skipped submatrix construction?
C     
 142                  CONTINUE
C 
                     IF(MYROW.EQ.RRSRC .AND. MYCOL.EQ.RCSRC) THEN
                        CALL DBEXMAT( EXROW, CEXCOL, ROWS, BEROWS2, IX, 
     $                       JX, F((LJFS-1)*LLDF+LIFS), LLDF, NBCBE, 
     $                       MB, NB, DWORK( EXF ), DWORK(FKJ),
     $                       MB + 1 )
                        IF( IDEEP.NE.IDEEPS .OR. JDEEP.NE.JDEEPS )
     $                       GO TO 144
                        IF( (RSRC.NE.RRSRC) .OR. (CSRC.NE.RCSRC) ) THEN
                           CALL DGERV2D( ICTXT, BEROWS2, 2*COLS, 
     $                          DWORK( BEUP+(KBEUP-1)*(NB+1)**2 ), 
     $                          NB + 1, RSRC, CSRC )
                        END IF
 144                    CONTINUE
                        XYRWS = IWORK( IXYRWS + IIS - I ) 
                        XYCLS = IWORK( IXYCLS + IIS - I ) 
                        XYRIND = IWORK( IXYRIND + IIS - I )
                        XYCIND = IWORK( IXYCIND + IIS - I )
                        CALL DGEMM( 'N', 'T', XYRWS, BEROWS2, XYCLS, 
     $                       USIGN, DWORK( XIJ1+(XYCIND-1)*(MB+1)+
     $                       XYRIND-1 ), MB + 1, DWORK( BEUP+
     $                       (KBEUP-1)*(NB+1)**2+(XYCIND-1)*(NB+1) ), 
     $                       NB + 1, ONE, DWORK( FKJ+XYRIND-1 ), 
     $                       MB + 1 )
                        CALL DGEMM( 'N', 'T', XYRWS, BEROWS2, XYCLS, 
     $                       USIGN, DWORK( YIJ1+(XYCIND-1)*(MB+1)+
     $                       XYRIND-1 ), MB+1, DWORK( BEUP+
     $                       (KBEUP-1)*(NB+1)**2+(XYCIND-1)*(NB+1)+
     $                       (NB+1)*COLS ), NB + 1, ONE, 
     $                       DWORK( FKJ+XYRIND-1 ), MB + 1 )
                        CALL DUBEXMA( EXROW, CEXCOL, ROWS, BEROWS2, IX, 
     $                       JX, F((LJFS-1)*LLDF+LIFS), LLDF, NBCBE, 
     $                       MB, NB, DWORK( EXF ), DWORK(FKJ), 
     $                       MB + 1 )
                         IF( BEUPBL.NE.1 ) KBEUP = KBEUP + 2
                     END IF
 140              CONTINUE
               END IF
C     
 55            CONTINUE
C     
C     Update inner loop variable j (column index)
C     
               IF( .NOT. TRANSZ ) THEN
                  IF( (.NOT.TRANSAD).AND.(.NOT.TRANSBE) ) THEN
                     J = J - 1
                  ELSEIF( TRANSAD.AND.(.NOT.TRANSBE) ) THEN
                     J = J + 1
                  ELSEIF( (.NOT.TRANSAD).AND.TRANSBE ) THEN
                     J = J + 1
                  ELSEIF( TRANSAD.AND.TRANSBE ) THEN
                     J = J - 1
                  END IF
               ELSE
                  IF( (.NOT.TRANSAD).AND.(.NOT.TRANSBE) ) THEN
                     J = J - 1
                  ELSEIF( TRANSAD.AND.TRANSBE ) THEN
                     J = J - 1
                  END IF
               END IF
 50         CONTINUE
C     
 24      CONTINUE
 22      CONTINUE
 20      CONTINUE
C     
C     Update outer loop variables
C     
         IF( .NOT. TRANSZ ) THEN
            IF( (.NOT.TRANSAD).AND.(.NOT.TRANSBE) ) THEN
               JS = MIN( JS + 1, DBBE )
               IEND = MAX( 1, IEND - 1 )
            ELSEIF( TRANSAD.AND.(.NOT.TRANSBE) ) THEN
               IF( K .GE. DBAD ) JS = JS + 1
               IS = MIN( IS + 1, DBAD )
            ELSEIF( (.NOT.TRANSAD).AND.TRANSBE ) THEN
               JS = MAX( JS - 1, 1 )
               IEND = MAX( 1, IEND - 1 )
            ELSEIF( TRANSAD.AND.TRANSBE ) THEN
               IF( K .GE. DBAD ) JS = JS - 1
               IS = MIN( IS + 1, DBAD )
            END IF
         ELSE
            IF( (.NOT.TRANSAD).AND.(.NOT.TRANSBE) ) THEN
               IF( K .GE. DBAD ) JS = JS - 1
               IS = MIN( IS + 1, DBAD )
            ELSEIF( TRANSAD.AND.TRANSBE ) THEN
               JS = MIN( JS + 1, DBBE )
               IEND = MAX( 1, IEND - 1 )
            END IF
         END IF
C
 10   CONTINUE
C
C     Before we go on we must back redistributed the elements in (C,F)
C
      IF( ADEXT.OR.BEEXT) THEN
         CALL PDBCKRD( M, N, C, IC, JC, DESCC, IWORK( EXADINF ), 
     $                 IWORK( EXBEINF ), DWORK( EXC ), EXMEMCF,
     $                 DWORK( SND ), MAX( MB, NB ), INFO )
         CALL PDBCKRD( M, N, F, IF, JF, DESCF, IWORK( EXADINF ), 
     $                 IWORK( EXBEINF ), DWORK( EXF ), EXMEMCF,
     $                 DWORK( SND ), MAX( MB, NB ), INFO )
      END IF
C
C     Restore Data in Original Position. Now we don't need to
C     shift the extended elements, we shift just the original
C     matrices. 
C
      IF( SHIFT ) THEN
         ADRSRC = DESCSA( RSRC_ )
         ADCSRC = DESCSA( CSRC_ )
         BERSRC = DESCSB( RSRC_ )
         BECSRC = DESCSB( CSRC_ )
         CFRSRC = DESCSC( RSRC_ )
         CFCSRC = DESCSC( CSRC_ )
         ASI = IA + MB * MOD( NPROW + MYROW - ADRSRC, NPROW ) - 
     $         IROFFAD 
         ASJ = JA + MB * MOD( NPCOL + MYCOL - ADCSRC, NPCOL ) - 
     $         IROFFAD
         BSI = IB + NB * MOD( NPROW + MYROW - BERSRC, NPROW ) - 
     $         IROFFBE
         BSJ = JB + NB * MOD( NPCOL + MYCOL - BECSRC, NPCOL ) - 
     $         IROFFBE
         CSI = IC + MB * MOD( NPROW + MYROW - CFRSRC, NPROW ) - 
     $         IROFFAD
         CSJ = JC + NB * MOD( NPCOL + MYCOL - CFCSRC, NPCOL ) - 
     $         IROFFBE
         DSI = ID + MB * MOD( NPROW + MYROW - ADRSRC, NPROW ) - 
     $         IROFFAD 
         DSJ = JD + MB * MOD( NPCOL + MYCOL - ADCSRC, NPCOL ) - 
     $         IROFFAD
         ESI = IE + NB * MOD( NPROW + MYROW - BERSRC, NPROW ) - 
     $         IROFFBE
         ESJ = JE + NB * MOD( NPCOL + MYCOL - BECSRC, NPCOL ) - 
     $         IROFFBE
         FSI = IF + MB * MOD( NPROW + MYROW - CFRSRC, NPROW ) - 
     $         IROFFAD
         FSJ = JF + NB * MOD( NPCOL + MYCOL - CFCSRC, NPCOL ) - 
     $         IROFFBE
         CALL INFOG2L( ASI, ASJ, DESCA, NPROW, NPCOL, MYROW, MYCOL, 
     $                 LIAS, LJAS, IDUM1, IDUM2 )
         CALL INFOG2L( BSI, BSJ, DESCB, NPROW, NPCOL, MYROW, MYCOL, 
     $                 LIBS, LJBS, IDUM1, IDUM2 )
         CALL INFOG2L( CSI, CSJ, DESCC, NPROW, NPCOL, MYROW, MYCOL, 
     $                 LICS, LJCS, IDUM1, IDUM2 )
         CALL INFOG2L( DSI, DSJ, DESCD, NPROW, NPCOL, MYROW, MYCOL, 
     $                 LIDS, LJDS, IDUM1, IDUM2 )
         CALL INFOG2L( ESI, ESJ, DESCE, NPROW, NPCOL, MYROW, MYCOL, 
     $                 LIES, LJES, IDUM1, IDUM2 )
         CALL INFOG2L( FSI, FSJ, DESCF, NPROW, NPCOL, MYROW, MYCOL, 
     $                 LIFS, LJFS, IDUM1, IDUM2 )
      END IF
C     
      IF (SHIFT .AND. NPROCS.GT.1 .AND. 
     $     (M.EQ.N.OR.N.EQ.NB.OR.M.EQ.MB)) THEN
C     
C     Shift (A,D) East if NPCOL > 1 and (A,D)**T West if NPCOL > 1
C     Shift (B,E) North if NPROW > 1 and (B,E)**T South if NPROW > 1
C     When shifting, also send/receive the extension elements
C     
         IF( TRANSZ ) THEN
            TRANSAD = .NOT. TRANSAD
            TRANSBE = .NOT. TRANSBE
         END IF
         IF ( NPCOL.GT.1.AND.M.NE.MB ) THEN
            IF( (.NOT.TRANSAD) .AND. (.NOT.TRANSBE) ) THEN
               CALL DGESD2D( ICTXT, ADROWS, ADCOLS, 
     $              A((LJAS-1)*LLDA+LIAS), LLDA, MYROW,
     $              EAST )
               DESCSA(CSRC_) = MOD(DESCSA(CSRC_) + 1, NPCOL)
               CALL DGERV2D( ICTXT, ADROWS, ADCOLS,
     $              A((LJAS-1)*LLDA+LIAS), LLDA, MYROW,
     $              WEST )
               CALL DGESD2D( ICTXT, ADROWS, ADCOLS,
     $              D((LJDS-1)*LLDD+LIDS), LLDD, MYROW,
     $              EAST )
               DESCSD(CSRC_) = MOD(DESCSD(CSRC_) + 1, NPCOL)
               CALL DGERV2D( ICTXT, ADROWS, ADCOLS,
     $              D((LJDS-1)*LLDD+LIDS), LLDD, MYROW,
     $              WEST )
            ELSEIF( TRANSAD .AND. TRANSBE ) THEN
               CALL DGESD2D( ICTXT, ADROWS, ADCOLS,
     $              A((LJAS-1)*LLDA+LIAS), LLDA, MYROW,
     $              WEST )  
               DESCSA(CSRC_) = MOD(DESCSA(CSRC_) + NPCOL - 1, NPCOL)
               CALL DGERV2D( ICTXT, ADROWS, ADCOLS,
     $              A((LJAS-1)*LLDA+LIAS), LLDA, MYROW,
     $              EAST )
               CALL DGESD2D( ICTXT, ADROWS, ADCOLS,
     $              D((LJDS-1)*LLDD+LIDS), LLDD, MYROW,
     $              WEST )  
               DESCSD(CSRC_) = MOD(DESCSD(CSRC_) + NPCOL - 1, NPCOL)
               CALL DGERV2D( ICTXT, ADROWS, ADCOLS,
     $              D((LJDS-1)*LLDD+LIDS), LLDD, MYROW,
     $              EAST )
            END IF
         END IF
         IF ( NPROW.GT.1.AND.N.NE.NB ) THEN
            IF( (.NOT.TRANSAD) .AND. (.NOT.TRANSBE) ) THEN
               CALL DGESD2D( ICTXT, BEROWS, BECOLS, 
     $              B((LJBS-1)*LLDB+LIBS), LLDB, NORTH,
     $              MYCOL )
               DESCSB(RSRC_) = MOD(DESCSB(RSRC_) + NPROW - 1,NPROW)
               CALL DGERV2D( ICTXT, BEROWS, BECOLS,
     $              B((LJBS-1)*LLDB+LIBS), LLDB, SOUTH,
     $              MYCOL )
               CALL DGESD2D( ICTXT, BEROWS, BECOLS,
     $              E((LJES-1)*LLDE+LIES), LLDE, NORTH,
     $              MYCOL )
               DESCSE(RSRC_) = MOD(DESCSE(RSRC_) + NPROW - 1,NPROW)
               CALL DGERV2D( ICTXT, BEROWS, BECOLS,
     $              E((LJES-1)*LLDE+LIES), LLDE, SOUTH,
     $              MYCOL )
            ELSEIF( TRANSAD .AND. TRANSBE ) THEN
               CALL DGESD2D( ICTXT, BEROWS, BECOLS,
     $              B((LJBS-1)*LLDB+LIBS), LLDB, SOUTH,
     $              MYCOL )
               DESCSB(RSRC_) = MOD(DESCSB(RSRC_) + 1, NPROW)
               CALL DGERV2D( ICTXT, BEROWS, BECOLS,
     $              B((LJBS-1)*LLDB+LIBS), LLDB, NORTH,
     $              MYCOL )
               CALL DGESD2D( ICTXT, BEROWS, BECOLS, 
     $              E((LJES-1)*LLDE+LIES), LLDE, SOUTH,
     $              MYCOL )
               DESCSE(RSRC_) = MOD(DESCSE(RSRC_) + 1, NPROW)
               CALL DGERV2D( ICTXT, BEROWS, BECOLS,
     $              E((LJES-1)*LLDE+LIES), LLDE, NORTH,
     $              MYCOL )
            END IF
         END IF
      ELSEIF (SHIFT .AND. NPROCS.GT.1 .AND. M.LT.N) THEN 
C     
C     Shift (A,D) SouthEast and (C,F) South  if NPROW > 1 
C     or Shift (A,D)**T NorthWest and (C,F) North if NPROW > 1 
C     
         IF( TRANSZ ) THEN
            TRANSAD = .NOT. TRANSAD
            TRANSBE = .NOT. TRANSBE
         END IF
         IF( (.NOT.TRANSAD) .AND. (.NOT.TRANSBE) ) THEN
            CALL DGESD2D( ICTXT, ADROWS, ADCOLS,
     $           A((LJAS-1)*LLDA+LIAS), LLDA, SOUTH, 
     $           EAST )
            DESCSA(RSRC_) = MOD(DESCSA(RSRC_) + 1, NPROW)
            DESCSA(CSRC_) = MOD(DESCSA(CSRC_) + 1, NPCOL)
            CALL DGERV2D( ICTXT, ADROWS, ADCOLS, 
     $           A((LJAS-1)*LLDA+LIAS), LLDA, NORTH,
     $           WEST )
            CALL DGESD2D( ICTXT, ADROWS, ADCOLS,
     $           D((LJDS-1)*LLDD+LIDS), LLDD, SOUTH,
     $           EAST )
            DESCSD(RSRC_) = MOD(DESCSD(RSRC_) + 1, NPROW)
            DESCSD(CSRC_) = MOD(DESCSD(CSRC_) + 1, NPCOL)
            CALL DGERV2D( ICTXT, ADROWS, ADCOLS,
     $           D((LJDS-1)*LLDD+LIDS), LLDD, NORTH,
     $           WEST )
         ELSEIF( TRANSAD .AND. TRANSBE ) THEN
            CALL DGESD2D( ICTXT, ADROWS, ADCOLS,
     $           A((LJAS-1)*LLDA+LIAS), LLDA, NORTH,
     $           WEST )
            DESCSA(RSRC_) = MOD(DESCSA(RSRC_) + NPROW - 1, NPROW)
            DESCSA(CSRC_) = MOD(DESCSA(CSRC_) + NPCOL - 1, NPCOL)
            CALL DGERV2D( ICTXT, ADROWS, ADCOLS, 
     $           A((LJAS-1)*LLDA+LIAS), LLDA, SOUTH,
     $           EAST )
            CALL DGESD2D( ICTXT, ADROWS, ADCOLS,
     $           D((LJDS-1)*LLDD+LIDS), LLDD, NORTH,
     $           WEST )
            DESCSD(RSRC_) = MOD(DESCSD(RSRC_) + NPROW - 1, NPROW)
            DESCSD(CSRC_) = MOD(DESCSD(CSRC_) + NPCOL - 1, NPCOL)
            CALL DGERV2D( ICTXT, ADROWS, ADCOLS,
     $           D((LJDS-1)*LLDD+LIDS), LLDD, SOUTH,
     $           EAST )
         END IF
         IF ( NPROW.GT.1 ) THEN
            IF( (.NOT.TRANSAD) .AND. (.NOT.TRANSBE) ) THEN
               CALL DGESD2D( ICTXT, ADROWS, BECOLS, 
     $              C((LJCS-1)*LLDC+LICS), LLDC, SOUTH,
     $              MYCOL )
               DESCSC(RSRC_) = MOD(DESCSC(RSRC_) + 1, NPROW)
               CALL DGERV2D( ICTXT, ADROWS, BECOLS, 
     $              C((LJCS-1)*LLDC+LICS), LLDC, NORTH,
     $              MYCOL )
               CALL DGESD2D( ICTXT, ADROWS, BECOLS,
     $              F((LJFS-1)*LLDF+LIFS), LLDF, SOUTH,
     $              MYCOL )
               DESCSF(RSRC_) = MOD(DESCSF(RSRC_) + 1, NPROW)
               CALL DGERV2D( ICTXT, ADROWS, BECOLS,
     $              F((LJFS-1)*LLDF+LIFS), LLDF, NORTH,
     $              MYCOL )
            ELSEIF( TRANSAD .AND. TRANSBE ) THEN
               CALL DGESD2D( ICTXT, ADROWS, BECOLS,
     $              C((LJCS-1)*LLDC+LICS), LLDC, NORTH,
     $              MYCOL )
               DESCSC(RSRC_) = MOD(DESCSC(RSRC_) + NPROW-1,NPROW)
               CALL DGERV2D( ICTXT, ADROWS, BECOLS,
     $              C((LJCS-1)*LLDC+LICS), LLDC, SOUTH,
     $              MYCOL )
               CALL DGESD2D( ICTXT, ADROWS, BECOLS,
     $              F((LJFS-1)*LLDF+LIFS), LLDF, NORTH,
     $              MYCOL )
               DESCSF(RSRC_) = MOD(DESCSF(RSRC_) + NPROW-1,NPROW)
               CALL DGERV2D( ICTXT, ADROWS, BECOLS,
     $              F((LJFS-1)*LLDF+LIFS), LLDF, SOUTH,
     $              MYCOL )
            END IF
         END IF
      ELSEIF (SHIFT .AND. NPROCS.GT.1) THEN 
C     
C     Shift (B,E) NorthWest and (C,F) West if NPCOL > 1 
C     or shift (B,E)**T SouthEast and (C,F) East if NPCOL > 1 
C     
         IF( TRANSZ ) THEN
            TRANSAD = .NOT. TRANSAD
            TRANSBE = .NOT. TRANSBE
         END IF
         IF( (.NOT.TRANSAD) .AND. (.NOT.TRANSBE) ) THEN
            CALL DGESD2D( ICTXT, BEROWS, BECOLS,
     $           B((LJBS-1)*LLDB+LIBS), LLDB, NORTH,
     $           WEST )
            DESCSB(RSRC_) = MOD(DESCSB(RSRC_) + NPROW - 1, NPROW)
            DESCSB(CSRC_) = MOD(DESCSB(CSRC_) + NPCOL - 1, NPCOL)
            CALL DGERV2D( ICTXT, BEROWS, BECOLS,
     $           B((LJBS-1)*LLDB+LIBS), LLDB, SOUTH,
     $           EAST )
            CALL DGESD2D( ICTXT, BEROWS, BECOLS,
     $           E((LJES-1)*LLDE+LIES), LLDE, NORTH,
     $           WEST )
            DESCSE(RSRC_) = MOD(DESCSE(RSRC_) + NPROW - 1, NPROW)
            DESCSE(CSRC_) = MOD(DESCSE(CSRC_) + NPCOL - 1, NPCOL)
            CALL DGERV2D( ICTXT, BEROWS, BECOLS,
     $           E((LJES-1)*LLDE+LIES), LLDE, SOUTH,
     $           EAST )
         ELSEIF( TRANSAD .AND. TRANSBE ) THEN
            CALL DGESD2D( ICTXT, BEROWS, BECOLS, 
     $           B((LJBS-1)*LLDB+LIBS), LLDB, SOUTH,
     $           EAST )
            DESCSB(RSRC_) = MOD(DESCSB(RSRC_) + 1, NPROW)
            DESCSB(CSRC_) = MOD(DESCSB(CSRC_) + 1, NPCOL)
            CALL DGERV2D( ICTXT, BEROWS, BECOLS,
     $           B((LJBS-1)*LLDB+LIBS), LLDB, NORTH, 
     $           WEST )
            CALL DGESD2D( ICTXT, BEROWS, BECOLS,
     $           E((LJES-1)*LLDE+LIES), LLDE, SOUTH,
     $           EAST )
            DESCSE(RSRC_) = MOD(DESCSE(RSRC_) + 1, NPROW)
            DESCSE(CSRC_) = MOD(DESCSE(CSRC_) + 1, NPCOL)
            CALL DGERV2D( ICTXT, BEROWS, BECOLS,
     $           E((LJES-1)*LLDE+LIES), LLDE, NORTH,
     $           WEST )
         END IF
         IF ( NPCOL.GT.1 ) THEN
            IF( (.NOT.TRANSAD) .AND. (.NOT.TRANSBE) ) THEN
               CALL DGESD2D( ICTXT, ADROWS, BECOLS,
     $              C((LJCS-1)*LLDC+LICS), LLDC, MYROW,
     $              WEST )
               DESCSC(CSRC_) = MOD(DESCSC(CSRC_) + NPCOL-1,NPCOL)
               CALL DGERV2D( ICTXT, ADROWS, BECOLS,
     $              C((LJCS-1)*LLDC+LICS), LLDC, MYROW, 
     $              EAST )
               CALL DGESD2D( ICTXT, ADROWS, BECOLS,
     $              F((LJFS-1)*LLDF+LIFS), LLDF, MYROW,
     $              WEST )
               DESCSF(CSRC_) = MOD(DESCSF(CSRC_) + NPCOL-1,NPCOL)
               CALL DGERV2D( ICTXT, ADROWS, BECOLS,
     $              F((LJFS-1)*LLDF+LIFS), LLDF, MYROW, 
     $              EAST )
            ELSEIF( TRANSAD .AND. TRANSBE ) THEN
               CALL DGESD2D( ICTXT, ADROWS, BECOLS,
     $              C((LJCS-1)*LLDC+LICS), LLDC, MYROW,
     $              EAST )
               DESCSC(CSRC_) = MOD(DESCSC(CSRC_) + 1, NPCOL)
               CALL DGERV2D( ICTXT, ADROWS, BECOLS,
     $              C((LJCS-1)*LLDC+LICS), LLDC, MYROW, 
     $              WEST )
               CALL DGESD2D( ICTXT, ADROWS, BECOLS,
     $              F((LJFS-1)*LLDF+LIFS), LLDF, MYROW,
     $              EAST )
               DESCSF(CSRC_) = MOD(DESCSF(CSRC_) + 1, NPCOL)
               CALL DGERV2D( ICTXT, ADROWS, BECOLS,
     $              F((LJFS-1)*LLDF+LIFS), LLDF, MYROW, 
     $              WEST )
            END IF
         END IF
      END IF
      IF( SHIFT .AND. TRANSZ ) THEN
         TRANSAD = .NOT. TRANSAD
         TRANSBE = .NOT. TRANSBE
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
C     End of PTRGCSYD
C     
C *** Last line of PTRGCSYD ***
