CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C         This file is a part of the SCASY software library.           C
C         See http://www.cs.umu.se/parallel/research/scasy             C
C         Contributors: Robert Granat and Bo Kågström.                 C
C         SCASY Version 0.10,    March 31, 2006.                       C
C         Copyright 2006, Umeå University, Sweden.                     C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE PDIMPRED( MODE, M, N, A, IA, JA, DESCA, B, IB, JB, 
     $                     DESCB, C, IC, JC, DESCC, EXTINFA, EXTINFB, 
     $                     EXA, LEXA, EXB, LEXB, EXC, LEXC, DWORK, 
     $                     LDWORK, INFO )
C
C  -- ScaLAPACK-style routine (preliminary version ) --
C     HPC2N and Department of Computing Science,
C     University of Umeå, Sweden.
C     Written by Robert Granat, (granat@cs.umu.se)
C     January 26, 2006.
C
      IMPLICIT NONE
C
C     .. Scalar arguments ..
      INTEGER              M, N, IA, JA, IB, JB, IC, JC, LEXA, LEXB, 
     $                     LEXC, LDWORK, INFO
C     ..
C     .. Array arguments ..
      CHARACTER            MODE( 3 )
      DOUBLE PRECISION     A( * ), B( * ), C( * ), DWORK( * ), EXA( * ),
     $                     EXB( * ), EXC( * )
      INTEGER              DESCA( * ), DESCB( * ), DESCC( * ), 
     $                     EXTINFA( * ), EXTINFB( * ) 
C     ..
C     
C  Purpose and description
C  =======================
C  This subroutine does an implicit redistribution of the elements
C  in the three matrices sub(A), sub(B) and sub(C) involved in the 
C  Sylvester Equation sub(A)X +/- Xsub(B) = sub(C), to avoid 
C  computational errors caused by splittings of any 2-by-2 diagonal 
C  block of sub(A) or sub(B). This is done based on the information 
C  in EXTINFA, EXTINFB, which have been constructed by the
C  subroutine PDEXTCHK.
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
C  MODE     (global input) CHARACTER array of length 3
C           This array decides what matrices to be considered
C           in the redistribution, as follows:
C             If MODE( 1 ) = 'A', A will be considered.
C             If MODE( 2 ) = 'B', B will be considered.
C             If MODE( 3 ) = 'C', C will be considered.
C           The values of MODE also decide what input/output
C           parameters that are referenced (see below).
C
C  Input/Output parameters
C
C  M         (global input) INTEGER
C            With MODE( 1 ) = 'A' or MODE( 3 ) = 'C':
C            The number of rows and columns of the global distributed 
C            matrix A. This is also the number of rows of the
C            global distributed matrix C. M >= 0.
C            Otherwise, M is not referenced.
C
C  N         (global input) INTEGER
C            With MODE( 2 ) = 'B' or MODE( 3 ) = C:
C            The number of rows and columns of the global distributed
C            matrix B. This is also the number of columns of the
C            global distributed matrix C. N >= 0.
C            Otherwise, N is not referenced.      
C
C  A         (local input/output) DOUBLE PRECISION array
C            With MODE( 1 ) = 'A': 
C            Array of dimension (LLD_A,LOCc(M)). Contains the local
C            pieces of the global distributed matrix A. On output,
C            the local parts of the distributed matrix A in real
C            Schur form.
C            Otherwise, A is not referenced.
C
C  IA        (global input) INTEGER
C            With MODE( 1 ) = 'A':       
C            Row start index for sub(A), i.e., the submatrix to operate
C            on. IA >= 1.
C            Otherwise, IA is not referenced.
C
C  JA        (global input) INTEGER
C            With MODE( 1 ) = 'A': 
C            Column start index for sub(A), i.e., the submatrix to 
C            operate on. JA = IA must hold.
C            Otherwise, JA is not referenced.
C
C  DESCA     (global and local input) INTEGER array of dimension DLEN_.
C            With MODE( 1 ) = 'A':      
C            The array descriptor for the global distributed matrix A.
C            Otherwise, DESCA is not referenced.
C
C  B         (local input) DOUBLE PRECISION array
C            With MODE( 2 ) = 'B':
C            Array of dimension (LLD_B,LOCc(N)). Contains the local
C            pieces of the global distributed matrix B. On output,
C            the local parts of the distributed matrix B in real
C            Schur form.
C            Otherwise, B is not referenced.
C
C  IB        (global input) INTEGER
C            With MODE( 2 ) = 'B':
C            Row start index for sub(B), i.e., the submatrix to operate
C            on. IB >= 1.
C            Otherwise, IB is not referenced.
C
C  JB        (global input) INTEGER
C            With MODE( 2 ) = 'B':
C            Column start index for sub(B), i.e., the submatrix to 
C            operate on. JB = IB must hold.
C            Otherwise, JB is not referenced.
C
C  DESCB     (global and local input) INTEGER array of dimension DLEN_.
C            With MODE( 2 ) = 'B':
C            The array descriptor for the global distributed matrix B.
C            Otherwise, DESCB is not referenced.
C
C  C         (local input) DOUBLE PRECISION array 
C            Array of dimension (LLD_C,LOCc(N)).
C            With MODE( 3 ) = 'C': 
C            On entry C contains the local pieces of the global 
C            distributed matrix C.
C            Otherwise, C is not referenced.
C
C  IC        (global input) INTEGER
C            With MODE( 3 ) = 'C': 
C            Row start index for sub(C), i.e., the submatrix to operate 
C            on. MOD(IC,MB_A) = MOD(IA,MB_A) must hold. 
C            Otherwise, IC is not referenced.
C
C  JC        (global input) INTEGER
C            With MODE( 3 ) = 'C':
C            Column start index for sub(C), i.e., the submatrix to 
C            operate on. MOD(JC,MB_B) = MOD(JB,MB_B) must hold.
C            Otherwise, JC is not referenced.
C
C  DESCC     (global and local input) INTEGER array of dimension DLEN_.
C            With MODE( 3 ) = 'C':
C            The array descriptor for the global distributed matrix C.
C            Otherwise, DESCC is not referenced.
C
C  EXTINFA   (global input) INTEGER array of dimension
C            ICEIL( M + MOD( IA - 1, MB_A ), MB_A ). 
C            If MODE( 1 ) = 'A' or MODE( 3 ) = 'C': Carries the 
C            extension information for the matrix A computed by
C            PDEXTCHK.
C            Otherwise, EXTINFA is not referenced.
C  
C  EXTINFB   (global input) INTEGER array of dimension
C            ICEIL( N + MOD( IJ - 1, MB_B ), MB_B ). 
C            If MODE( 2 ) = 'B' or MODE( 3 ) = 'C': Carries the 
C            extension information for the matrix B computed by
C            PDEXTCHK.
C            Otherwise, EXTINFB is not referenced.
C            
C  EXA       (local output) DOUBLE PRECISION array, dimension LEXA.
C            With MODE( 1 ) = 'A':
C            Stores the extension elements associated with the matrix
C            A for usage in calls to the extended block building 
C            routines DBEXMAT and DUBEXMA.
C            Otherwise, EXA is not referenced.
C
C  LEXA      (local input) INTEGER, the length of the array EXA.
C            With MODE( 1 ) = 'A':
C            LEXA >= ICEIL( LOCr( M + MOD( IA - 1, MB_A ), MB_A ) ) * 
C                    ICEIL( LOCc( M + MOD( JA - 1, MB_A ), MB_A ) ) * 
C                    ( 2 * MB_A + 1 )
C            Otherwise, LEXA is not referenced.
C
C  EXB       (local output) DOUBLE PRECISION array, dimension LEXB.
C            With MODE( 2 ) = 'B':
C            Stores the extension elements associated with the matrix
C            B for usage in calls to the extended block building 
C            routines DBEXMAT and DUBEXMA.
C            Otherwise, EXB is not referenced.
C
C  LEXB      (local input) INTEGER, the length of the array EXB.
C            With MODE( 2 ) = 'B':
C            LEXB >= ICEIL( LOCr( N + MOD( IB - 1, NB_A ), NB_A ) ) * 
C                    ICEIL( LOCc( N + MOD( JB - 1, NB_A ), NB_A ) ) * 
C                    ( 2 * NB_A + 1 )
C            Otherwise, LEXB is not referenced.
C
C  EXC       (local output) DOUBLE PRECISION array, dimension LEXC.
C            With MODE( 3 ) = 'C':
C            Stores the extension elements associated with the matrix
C            C for usage in calls to the extended block building 
C            routines DBEXMAT and DUBEXMA.
C            Otherwise, EXC is not referenced.
C
C  LEXC      (local input) INTEGER, the length of the array EXC.
C            With MODE( 3 ) = 'C':
C            LEXC >= ICEIL( LOCr( M + MOD( IC - 1, MB_C ), MB_C ) ) * 
C                    ICEIL( LOCc( N + MOD( JC - 1, NB_C ), NB_C ) ) * 
C                    ( MB_C + NB_C + 1 )
C            Otherwise, LEXC is not referenced.
C
C  Workspace
C
C  DWORK     (local workspace) DOUBLE PRECISION array, dimension
C            LDWORK. With MODE( 1 ) = 'A' or MODE( 2 ) = 'B' or
C            MODE( 3 ) this array is referenced, otherwise not.
C
C  LDWORK    (local or global input) INTEGER
C            The dimension of the array DWORK.
C            With MODE( 1 ) = 'A' or MODE( 2 ) = 'B' or MODE( 3 )
C            LDWORK >= MAX( MB_A, NB_B ).
C            Otherwise, LDWORK is not referenced.
C 
C            If LDWORK = -1, LDWORK is global input and a workspace 
C            query is assumed. The routine will then calculate the 
C            optimal workspace needed, store it in DWORK(1) and return
C            immediately. No error will then be signaled by PXERBLA.
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
C
C  Method
C  ======
C  See the references.
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
C  Limitations
C  ===========
C  None.
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
C  None.
C
C  Revisions
C  =========
C  Please report bugs to <granat@cs.umu.se>.
C
C  Keywords
C  ========
C  Implicit redistribution, real Schur form.
C
C  ======================================================================
C  
C     .. Local Parameters ..
      INTEGER          CTXT_, MB_, NB_, RSRC_, CSRC_, LLD_, DLEN_ 
      DOUBLE PRECISION ZERO
      PARAMETER        ( CTXT_ = 2, MB_ = 5, NB_ = 6, RSRC_ = 7, 
     $                   CSRC_ = 8, LLD_ = 9, DLEN_ = 9, ZERO = 0.0D+0 )
C     ..
C     .. Local Scalars ..
      LOGICAL          LQUERY, A_CARE, B_CARE, C_CARE
      INTEGER          ICTXT, NPROW, NPCOL, MYROW, MYCOL, MB, NB, LLDA,
     $                 LLDB, LLDC, I, J, K, ISTART, JSTART, DBA, DBB, 
     $                 IX, JX, RSRC, CSRC, LBR, LBC, ACOLS, BCOLS, IS, 
     $                 JS, SRSRC, SCSRC, NBC, POS, NPROCS, LEN, IROFFA,
     $                 IROFFB, LIA, LJA, ARSRC, ACSRC, LIB, LJB, BRSRC, 
     $                 BCSRC, SAROWS, SACOLS, SBROWS, SBCOLS, ISX, JSX, 
     $                 ICOFFA, ICOFFB, EXMEMA, EXMEMB, EXMEMC, LIC, LJC,
     $                 RSRC1, CSRC1
      DOUBLE PRECISION CORNER, LOST_ELEM
C     ..
C     .. Local Arrays ..
      INTEGER          DESCSA( DLEN_ ), DESCSB( DLEN_ ), DESCSC( DLEN_ )
C     ..
C     .. External Subroutines ..
      EXTERNAL         INFOG2L, DGESD2D, DGERV2D, DLACPY, BLACS_GRIDINFO
C     ..
C     .. External functions ..
      LOGICAL          LSAME
      INTEGER          ICEIL, NUMROC
      EXTERNAL         LSAME, ICEIL, NUMROC
C     ..
C     .. Executable Statements ..
C 
      ICTXT = DESCA( CTXT_ )
      CALL BLACS_GRIDINFO( ICTXT, NPROW, NPCOL, MYROW, MYCOL )
C
C     Check if the context is valid
C
      IF( NPROW.EQ.-1 ) THEN
         INFO = 1
      END IF
C
C     Check input arguments
C
      LQUERY = LDWORK.EQ.-1
      A_CARE = LSAME( 'A', MODE( 1 ) )
      B_CARE = LSAME( 'B', MODE( 2 ) )
      C_CARE = LSAME( 'C', MODE( 3 ) )
      IF( INFO.EQ.0 ) THEN
         IF( A_CARE )
     $        CALL CHK1MAT( M, 2, M, 2, IA, JA, DESCA, 7, INFO )
         IF( INFO.EQ.0 ) THEN
            IF( A_CARE .AND. IA.LT.1 ) THEN
               INFO = -5
            ELSEIF( A_CARE .AND. JA.NE.IA ) THEN
               INFO = -6
            END IF
         END IF
      END IF
      IF( INFO.EQ.0 ) THEN
         IF( B_CARE )
     $        CALL CHK1MAT( N, 3, N, 3, IB, JB, DESCB, 11, INFO )
         IF( INFO.EQ.0 ) THEN
            IF( B_CARE .AND. IB.LT.1 ) THEN
               INFO = -9
            ELSEIF( B_CARE .AND. JB.NE.IB ) THEN
               INFO = -10
            END IF
         END IF
      END IF
      IF( INFO.EQ.0 ) THEN
         IF( C_CARE ) 
     $        CALL CHK1MAT( M, 2, N, 3, IC, JC, DESCC, 15, INFO )
         IF( INFO.EQ.0 ) THEN
            IF( C_CARE .AND. IC.LT.1 ) THEN
               INFO = -13
            ELSEIF( C_CARE .AND. JC.LT.1 ) THEN
               INFO = -14
            END IF
         END IF
      END IF
C
C     Compute the number of diagonal blocks in A and B and some other
C     values
C  
      IF( INFO.EQ.0 ) THEN
         IF( A_CARE ) THEN
            MB = DESCA( MB_ )
         ELSEIF( C_CARE ) THEN
            MB = DESCC( MB_ )
         ELSE
            MB = 1
         END IF
         IF( B_CARE ) THEN
            NB = DESCB( MB_ )
         ELSEIF( C_CARE ) THEN
            NB = DESCC( NB_ )
         ELSE
            NB = 1
         END IF
         IF( A_CARE ) THEN
            IROFFA = MOD( IA - 1, MB ) 
            DBA = ICEIL( M + IROFFA, MB )
         ELSEIF( C_CARE ) THEN
            IA = IC
            IROFFA = MOD( IC - 1, MB )
            DBA = ICEIL( M + IROFFA, MB )
         ELSE
            IA = 1
            IROFFA = 0
            DBA = 1
         END IF
         IF( B_CARE ) THEN
            IROFFB = MOD( IB - 1, NB )
            DBB = ICEIL( N + IROFFB, NB )
         ELSEIF( C_CARE ) THEN
            IB = JC
            IROFFB = MOD( JC - 1, NB )
            DBB = ICEIL( N + IROFFB, NB )
         ELSE
            IB = 1
            IROFFB = 0
            DBB = 1
         END IF
         IF( A_CARE ) LLDA = DESCA( LLD_ )
         IF( B_CARE ) LLDB = DESCB( LLD_ )
         IF( C_CARE ) LLDC = DESCC( LLD_ )
         NPROCS = NPROW * NPCOL
C     
C     Check who owns the first row and column of sub(A) and sub(B)
C     
         IF( A_CARE )
     $        CALL INFOG2L( IA, JA, DESCA, NPROW, NPCOL, MYROW, MYCOL, 
     $                      LIA, LJA, ARSRC, ACSRC )
         IF( B_CARE ) 
     $        CALL INFOG2L( IB, JB, DESCB, NPROW, NPCOL, MYROW, MYCOL, 
     $                      LIB, LJB, BRSRC, BCSRC )
         IF( C_CARE )
     $        CALL INFOG2L( IC, JC, DESCC, NPROW, NPCOL, MYROW, MYCOL, 
     $                      LIC, LJC, ARSRC, BCSRC ) 
C     
C     Init matrix descriptors for smallest submatrices of A, B and C
C     including sub(A), sub(B) and sub(C) which conform with ScaLAPACK
C     conventions. We need those descriptors for simplifying the
C     computations of indicies in EXA, EXB and EXC later on.
C     
         IF( A_CARE ) THEN
            SAROWS = NUMROC( M + IROFFA, MB, MYROW, ARSRC, NPROW )
            SACOLS = NUMROC( M + IROFFA, MB, MYCOL, ACSRC, NPCOL )
            CALL DESCINIT ( DESCSA, M + IROFFA, M + IROFFA, MB, MB, 
     $                      ARSRC, ACSRC, ICTXT, LLDA, INFO )
            EXMEMA = ICEIL(SAROWS,MB)*ICEIL(SACOLS,MB)*(2*MB+1)
         END IF
         IF( B_CARE ) THEN
            SBROWS = NUMROC( N + IROFFB, NB, MYROW, BRSRC, NPROW )
            SBCOLS = NUMROC( N + IROFFB, NB, MYCOL, BCSRC, 
     $                       NPCOL )
            CALL DESCINIT ( DESCSB, N + IROFFB, N + IROFFB, NB, NB, 
     $                      BRSRC, BCSRC, ICTXT, LLDB, INFO )
            EXMEMB = ICEIL(SBROWS,NB)*ICEIL(SBCOLS,NB)*(2*NB+1)
         END IF
         IF( C_CARE ) THEN
            SAROWS = NUMROC( M + IROFFA, MB, MYROW, ARSRC, NPROW )
            SBCOLS = NUMROC( N + IROFFB, NB, MYCOL, BCSRC, NPCOL )
            CALL DESCINIT ( DESCSC, M + IROFFA, N + IROFFB, DESCC(MB_), 
     $                      DESCC(NB_), ARSRC, BCSRC, ICTXT, LLDC, 
     $                      INFO )
            EXMEMC = ICEIL(SAROWS,MB)*ICEIL(SBCOLS,NB)*(MB+NB+1)
         END IF
C     
C     Check length of output arrays
C     
         IF( A_CARE .AND. LEXA.LT.EXMEMA ) INFO = -19
         IF( INFO.EQ.0 ) THEN
            IF( B_CARE .AND. LEXB.LT.EXMEMB ) INFO = -21
         END IF
         IF( INFO.EQ.0 ) THEN
            IF( C_CARE .AND. LEXC.LT.EXMEMC ) INFO = -23
         END IF
      END IF
C     
C     Check workspace
C
      IF( INFO.EQ.0 .AND. (A_CARE .OR. B_CARE .OR. C_CARE )) THEN
         IF( .NOT. LQUERY .AND. LDWORK.LT. MAX( MB, NB ) ) THEN
            INFO = -19
         ELSEIF( LQUERY ) THEN
            INFO = 0
            DWORK( 1 ) = MAX( MB, NB )
            RETURN
         END IF
      END IF
C     
C     Error messages
C
      IF( INFO.NE.0 ) THEN
         CALL PXERBLA( ICTXT, 'PDIMPRED', -INFO )
      END IF 
C
C     Init the extension arrays with zeros so that we have zeros
C     where the extensions aren't done
C
      IF( A_CARE ) THEN
         DO 5 I = 1, LEXA
            EXA( I ) = ZERO
 5       CONTINUE
      END IF
C
      IF( B_CARE ) THEN
         DO 6 I = 1, LEXB
            EXB( I ) = ZERO
 6       CONTINUE
      END IF
C
      IF( C_CARE ) THEN
         DO 7 I = 1, LEXC
            EXC( I ) = ZERO
 7       CONTINUE
      END IF
C
C     Should we are about A? If not, consider B instead!
C
      IF( .NOT. A_CARE ) GO TO 35
C
C     Extend A on basis of the information in EXTINFA
C
C     Loop through all blocks in A: treat diagonal blocks and other
C     blocks differently. Since A and B are upper (quasi-)triangular
C     we only consider blocks Aij and Bij where j >= i.
C
      NBC = ICEIL( SACOLS, MB )
C
C     Now do the 4 sweeps over A
C
      DO 10 K = 1, 4
C
C     Select which MBxNB block in each 2x2 block to operate on
C     in this sweep
C
         IF( K.EQ.1.OR.K.EQ.2 ) ISTART = 1
         IF( K.EQ.3.OR.K.EQ.4 ) ISTART = 2
         IF( K.EQ.1.OR.K.EQ.3 ) JSTART = 1
         IF( K.EQ.2.OR.K.EQ.4 ) JSTART = 2
C     
C     Loop through all those blocks in A
C     
         DO 20 I = ISTART, DBA, 2
            DO 30 J = JSTART, DBA, 2
C     
C     Check if diagonal block
C     
               IF( I.EQ.J ) THEN
C     
C     Check if we shall extend Aii with one row from the South, one
C     column from the East and one corner element from South-East
C     
                  IF( EXTINFA( I ).EQ.1 .OR. EXTINFA( I ).EQ.3 ) THEN
                     CALL INFOG2L( (I-1) * MB + IA - IROFFA, 
     $                             (J-1) * MB + IA - IROFFA, 
     $                             DESCA, NPROW, NPCOL, MYROW, MYCOL, 
     $                             IX, JX, RSRC, CSRC )
                     CALL INFOG2L( (I-1) * MB + 1, 
     $                             (J-1) * MB + 1, 
     $                             DESCSA, NPROW, NPCOL, MYROW, MYCOL, 
     $                             ISX, JSX, RSRC1, CSRC1 )
                     CALL INFOG2L( I * MB + IA - IROFFA, 
     $                             (J-1) * MB + IA - IROFFA, 
     $                             DESCA, NPROW, NPCOL, MYROW, MYCOL, 
     $                             IS, JS, SRSRC, SCSRC)
C     
C     If I own the extra row send it to the holder of the diagonal block
C     
                     IF( (MYROW.EQ.RSRC.AND.MYCOL.EQ.CSRC) .OR.
     $                    (MYROW.EQ.SRSRC.AND.MYCOL.EQ.SCSRC) ) THEN
                        IF( MYROW.EQ.SRSRC.AND.MYCOL.EQ.SCSRC ) THEN
                           LOST_ELEM = A((JS-1)*LLDA+IS+(MB-1)*LLDA)
                           IF( NPROW.GT.1 )
     $                          CALL DGESD2D( ICTXT, 1, 1, LOST_ELEM, 
     $                                        1, RSRC, CSRC)
                        END IF
C     
C     Receive extra row from process holding the block in the South
C     
                        IF( MYROW.EQ.RSRC.AND.MYCOL.EQ.CSRC ) THEN
                           IF( NPROW.GT.1 )
     $                          CALL DGERV2D( ICTXT, 1, 1, LOST_ELEM,
     $                                        1, SRSRC, SCSRC )
                           LBR = ICEIL( ISX, MB )
                           LBC = ICEIL( JSX, MB )
                           POS = ( (LBR-1)*NBC + (LBC-1) )*(2*MB+1) +
     $                          1 + MB - 1
                           EXA( POS ) = LOST_ELEM
                        END IF
                     END IF
C     
C     If I own extra corner element send it to the holder of the diag.
C     block
C     
                     CALL INFOG2L( I * MB + IA - IROFFA, 
     $                             J * MB + IA - IROFFA,
     $                             DESCA, NPROW, NPCOL, MYROW, MYCOL, 
     $                             IS, JS, SRSRC, SCSRC )
                     IF( (MYROW.EQ.RSRC.AND.MYCOL.EQ.CSRC) .OR.
     $                    (MYROW.EQ.SRSRC.AND.MYCOL.EQ.SCSRC) ) THEN
                        IF( MYROW.EQ.SRSRC.AND.MYCOL.EQ.SCSRC ) THEN
                           CORNER = A((JS-1)*LLDA + IS)
                           IF( NPROCS.GT.1 )
     $                          CALL DGESD2D( ICTXT, 1, 1, CORNER, 1, 
     $                                        RSRC, CSRC )
                        END IF
C     
C     Receive extra corner from process in the South-East
C     
                        IF( MYROW.EQ.RSRC.AND.MYCOL.EQ.CSRC ) THEN
                           IF( NPROCS.GT.1 )
     $                          CALL DGERV2D( ICTXT, 1, 1, CORNER, 1,
     $                                        SRSRC, SCSRC )
                           POS = POS + 1
                           EXA( POS ) = CORNER
                        END IF
                     END IF
C     
C     If I own the extra column send it to the holder of the diagonal
C     block
C     
                     CALL INFOG2L( (I-1) * MB + IA - IROFFA, 
     $                             J * MB + IA - IROFFA, 
     $                             DESCA, NPROW, NPCOL, MYROW, MYCOL,
     $                             IS, JS, SRSRC, SCSRC )
                     IF( (MYROW.EQ.RSRC.AND.MYCOL.EQ.CSRC) .OR.
     $                    (MYROW.EQ.SRSRC.AND.MYCOL.EQ.SCSRC) ) THEN
                        IF( MYROW.EQ.SRSRC.AND.MYCOL.EQ.SCSRC ) THEN
                           CALL DLACPY( 'All', MB, 1, A((JS-1)*LLDA+IS),
     $                                  LLDA, DWORK, MB )
                           IF( NPCOL.GT.1 )
     $                          CALL DGESD2D( ICTXT, MB, 1, DWORK, MB,
     $                                        RSRC, CSRC )
                        END IF
C     
C     Receive extra column from process holding the block in the East
C     
                        IF( MYROW.EQ.RSRC.AND.MYCOL.EQ.CSRC ) THEN
                           IF( NPCOL.GT.1 )
     $                          CALL DGERV2D( ICTXT, MB, 1, DWORK, MB,
     $                                        SRSRC, SCSRC) 
                           POS = POS + 1
                           CALL DLACPY( 'All', MB, 1, DWORK, MB, 
     $                                  EXA(POS), 2 * MB + 1 )
                        END IF
                     END IF
                  END IF
C     
C     Not diagonal block and a block in the upper triangular part of A
C     
               ELSEIF( J.GT.I ) THEN
C     
C     Check if Aii extended, if so extend Aij with one row from
C     the South
C     
                  IF( EXTINFA( I ).EQ.1 .OR. EXTINFA( I ).EQ.3 ) THEN
                     CALL INFOG2L( (I-1) * MB + IA - IROFFA, 
     $                             (J-1) * MB + IA - IROFFA, 
     $                             DESCA, NPROW,NPCOL, MYROW, MYCOL, 
     $                             IX, JX, RSRC, CSRC )
                     CALL INFOG2L( (I-1) * MB + 1, 
     $                             (J-1) * MB + 1, 
     $                             DESCSA, NPROW, NPCOL, MYROW, MYCOL, 
     $                             ISX, JSX, RSRC1, CSRC1 )
                     CALL INFOG2L( I * MB + IA - IROFFA, 
     $                             (J-1) * MB + IA - IROFFA, DESCA,
     $                             NPROW, NPCOL, MYROW, MYCOL, IS, JS, 
     $                             SRSRC, SCSRC )
C     
                     IF( (MYROW.EQ.RSRC.AND.MYCOL.EQ.CSRC) .OR.
     $                    (MYROW.EQ.SRSRC.AND.MYCOL.EQ.SCSRC) ) THEN
C     
C     If I own the extra row send it to the holder of the current block
C     
                        LEN = MIN( MB, M + IROFFA - (J-1)*MB )
                        IF( MYROW.EQ.SRSRC.AND.MYCOL.EQ.SCSRC ) THEN
                           CALL DLACPY( 'All', 1, LEN, 
     $                                  A((JS-1)*LLDA+IS), LLDA, DWORK, 
     $                                  1 )
                           IF( NPROW.GT.1 )
     $                          CALL DGESD2D( ICTXT, LEN, 1, DWORK, MB,
     $                                        RSRC, CSRC )
                        END IF
C     
C     Receive extra row from process holding the block in the South
C     
                        IF( MYROW.EQ.RSRC.AND.MYCOL.EQ.CSRC ) THEN
                           IF( NPROW.GT.1 )
     $                          CALL DGERV2D( ICTXT, LEN, 1, DWORK, MB,
     $                                        SRSRC, SCSRC )
                           LBR = ICEIL( ISX, MB )
                           LBC = ICEIL( JSX, MB )
                           POS = ( (LBR-1)*NBC + (LBC-1) )*( 2*MB+1)+ 1 
                           CALL DLACPY( 'All', LEN, 1, DWORK, MB,
     $                                  EXA(POS), 2 * MB + 1 )
                        END IF
                     END IF
                  END IF
C     
C     Check if Aii and Ajj is extended, if so extend with one extra 
C     corner element from South-East
C     
                  IF(( EXTINFA( I ).EQ.1 .OR. EXTINFA( I ).EQ.3 ).AND.
     $                 ( EXTINFA( J ).EQ.1 .OR. EXTINFA( J ).EQ.3 )) 
     $                 THEN
C     
C     If I own extra corner element send it to the holder of the diag.
C     block
C     
                     CALL INFOG2L( I * MB + IA - IROFFA, 
     $                             J * MB + IA - IROFFA,
     $                             DESCA, NPROW, NPCOL, MYROW, MYCOL, 
     $                             IS, JS, SRSRC, SCSRC )
                     IF( (MYROW.EQ.RSRC.AND.MYCOL.EQ.CSRC) .OR.
     $                    (MYROW.EQ.SRSRC.AND.MYCOL.EQ.SCSRC) ) THEN
                        IF( MYROW.EQ.SRSRC.AND.MYCOL.EQ.SCSRC ) THEN
                           CORNER = A( (JS-1)*LLDA + IS )
                           IF( NPROCS.GT.1 )
     $                          CALL DGESD2D( ICTXT, 1, 1, CORNER, 1, 
     $                                        RSRC, CSRC )
                        END IF
C     
C     Receive extra corner from process in the South-East
C     
                        IF( MYROW.EQ.RSRC.AND.MYCOL.EQ.CSRC ) THEN
                           IF( NPROCS.GT.1 )
     $                          CALL DGERV2D( ICTXT, 1, 1, CORNER, 1,
     $                                        SRSRC, SCSRC )
                           POS = POS + MB 
                           EXA( POS ) = CORNER
                        END IF 
                     END IF
                  END IF
C     
C     Check if Ajj is extended, if so extend with one column from
C     the East
C     
                  IF( EXTINFA( J ).EQ.1 .OR. EXTINFA( J ).EQ.3 ) THEN
                     CALL INFOG2L( (I-1) * MB + IA - IROFFA, 
     $                             (J-1) * MB + IA - IROFFA, 
     $                             DESCA, NPROW, NPCOL, MYROW, MYCOL, 
     $                             IX, JX, RSRC, CSRC )
                     CALL INFOG2L( (I-1) * MB + 1, 
     $                             (J-1) * MB + 1, 
     $                             DESCSA, NPROW, NPCOL, MYROW, MYCOL, 
     $                             ISX, JSX, RSRC1, CSRC1 )
                     CALL INFOG2L( (I-1) * MB + IA - IROFFA, 
     $                             J * MB + IA - IROFFA,
     $                             DESCA, NPROW, NPCOL, MYROW, MYCOL,
     $                             IS, JS, SRSRC, SCSRC )
C
                     IF( (MYROW.EQ.RSRC.AND.MYCOL.EQ.CSRC) .OR.
     $                    (MYROW.EQ.SRSRC.AND.MYCOL.EQ.SCSRC) ) THEN
                        
C     If I own the extra column send it to the holder of the diagonal
C     block
C     
                        IF( MYROW.EQ.SRSRC.AND.MYCOL.EQ.SCSRC ) THEN
                           CALL DLACPY( 'All', MB, 1, A((JS-1)*LLDA+IS),
     $                                  LLDA, DWORK, MB )
                           IF( NPCOL.GT.1 )
     $                          CALL DGESD2D( ICTXT, MB, 1, DWORK, MB,
     $                                        RSRC, CSRC)
                        END IF
C     
C     Receive extra column from process holding the block in the East
C     
                        IF( MYROW.EQ.RSRC.AND.MYCOL.EQ.CSRC ) THEN
                           IF( NPCOL.GT.1 )
     $                          CALL DGERV2D( ICTXT, MB, 1, DWORK, MB,
     $                                        SRSRC, SCSRC )
                           LBR = ICEIL( ISX, MB )
                           LBC = ICEIL( JSX, MB )
                           POS = ( (LBR-1)*NBC + (LBC-1) )*( 2 *MB + 1)+
     $                          1 + MB + 1
                           CALL DLACPY( 'All', MB, 1, DWORK, MB, 
     $                                  EXA(POS), 2 * MB + 1 )
                        END IF
                     END IF
                  END IF
               END IF
 30         CONTINUE
 20      CONTINUE
 10   CONTINUE
C
 35   CONTINUE
C
C     Should we care about B at all? If not, consider C instead!
C
      IF( .NOT. B_CARE ) GO TO 65
C     
C     Now extend B in the same way as we did with A right above.
C     
C     Loop through all blocks in B: threat diagonal blocks and other
C     blocks differently. Only consider Bij where j >= i.
C     
      NBC = ICEIL( SBCOLS, NB )
C     
C     Now do the 4 sweeps over B
C     
      DO 40 K = 1, 4
C
C     Select which MBxNB block in each 2x2 block to operate on
C     in this sweep
C
        IF( K.EQ.1.OR.K.EQ.2 ) ISTART = 1
        IF( K.EQ.3.OR.K.EQ.4 ) ISTART = 2
        IF( K.EQ.1.OR.K.EQ.3 ) JSTART = 1
        IF( K.EQ.2.OR.K.EQ.4 ) JSTART = 2
C
C     Loop through all those blocks in B
C
        DO 50 I = ISTART, DBB, 2
           DO 60 J = JSTART, DBB, 2
C     
C     Check if diagonal block
C     
              IF( I.EQ.J ) THEN
C     
C     Check if we shall extend Bii with one row from the South, one
C     column from the East and one corner element from South-East
C     
                 IF( EXTINFB( I ).EQ.1 .OR. EXTINFB( I ).EQ.3 ) THEN
                    CALL INFOG2L( (I-1) * NB + IB - IROFFB, 
     $                            (J-1) * NB + IB - IROFFB,
     $                            DESCB, NPROW, NPCOL, MYROW, MYCOL, 
     $                            IX, JX, RSRC, CSRC )
                    CALL INFOG2L( (I-1) * NB + 1, 
     $                            (J-1) * NB + 1,
     $                            DESCSB, NPROW, NPCOL, MYROW, MYCOL, 
     $                            ISX, JSX, RSRC1, CSRC1 )
                    CALL INFOG2L( I * NB + IB - IROFFB, 
     $                            (J-1) * NB + IB - IROFFB ,
     $                            DESCB, NPROW, NPCOL, MYROW, MYCOL, IS,
     $                            JS, SRSRC, SCSRC )
C
                    IF( (MYROW.EQ.RSRC.AND.MYCOL.EQ.CSRC) .OR.
     $                   (MYROW.EQ.SRSRC.AND.MYCOL.EQ.SCSRC) ) THEN
C     
C     If I own the extra row send it to the holder of the diagonal block
C     
                       IF( MYROW.EQ.SRSRC.AND.MYCOL.EQ.SCSRC ) THEN
                          LOST_ELEM = B((JS-1)*LLDB+IS+(NB-1)*LLDB)
                          IF( NPROW.GT.1 )
     $                         CALL DGESD2D( ICTXT, 1, 1, LOST_ELEM, 1,
     $                                       RSRC, CSRC )
                       END IF
C     
C     Receive extra row from process holding the block in the South
C     
                       IF( MYROW.EQ.RSRC.AND.MYCOL.EQ.CSRC ) THEN
                          IF( NPROW.GT.1 )
     $                         CALL DGERV2D( ICTXT, 1, 1, LOST_ELEM, 1,
     $                                       SRSRC, SCSRC)
                          LBR = ICEIL( ISX, NB )
                          LBC = ICEIL( JSX, NB )
                          POS = ( (LBR-1)*NBC + (LBC-1) )*( 2 *NB + 1) +
     $                          1 + NB - 1 
                          EXB( POS ) = LOST_ELEM
                       END IF
                    END IF
C     
C     If I own extra corner element send it to the holder of the diag.
C     block
C     
                    CALL INFOG2L( I * NB + IB - IROFFB, 
     $                            J * NB + IB - IROFFB,
     $                            DESCB, NPROW, NPCOL, MYROW, MYCOL, 
     $                            IS, JS, SRSRC, SCSRC )
C
                    IF( (MYROW.EQ.RSRC.AND.MYCOL.EQ.CSRC) .OR.
     $                   (MYROW.EQ.SRSRC.AND.MYCOL.EQ.SCSRC) ) THEN
                       IF( MYROW.EQ.SRSRC.AND.MYCOL.EQ.SCSRC ) THEN
                          CORNER = B((JS-1)*LLDB + IS)
                          IF( NPROCS.GT.1 )
     $                         CALL DGESD2D( ICTXT, 1, 1, CORNER, 1, 
     $                                       RSRC, CSRC )
                       END IF
C     
C     Receive extra corner from process in the South-East
C     
                       IF( MYROW.EQ.RSRC.AND.MYCOL.EQ.CSRC ) THEN
                          IF( NPROCS.GT.1 )
     $                         CALL DGERV2D( ICTXT, 1, 1, CORNER, 1,
     $                                       SRSRC, SCSRC )
                          POS = POS + 1
                          EXB( POS ) = CORNER
                       END IF      
                    END IF
C     
C     If I own the extra column send it to the holder of the diagonal
C     block
C     
                    CALL INFOG2L( (I-1) * NB + IB - IROFFB , 
     $                            J * NB + IB - IROFFB, 
     $                            DESCB, NPROW, NPCOL, MYROW, MYCOL, IS,
     $                            JS, SRSRC, SCSRC )
C
                    IF( (MYROW.EQ.RSRC.AND.MYCOL.EQ.CSRC) .OR.
     $                   (MYROW.EQ.SRSRC.AND.MYCOL.EQ.SCSRC) ) THEN
                       IF( MYROW.EQ.SRSRC.AND.MYCOL.EQ.SCSRC ) THEN
                          CALL DLACPY( 'All', NB, 1, B((JS-1)*LLDB+IS),
     $                                 LLDB, DWORK, NB )
                          IF( NPCOL.GT.1 )
     $                         CALL DGESD2D( ICTXT, NB, 1, DWORK, NB, 
     $                                       RSRC, CSRC )
                       END IF
C     
C     Receive extra column from process holding the block in the East
C     
                       IF( MYROW.EQ.RSRC.AND.MYCOL.EQ.CSRC ) THEN
                          IF( NPCOL.GT.1 )
     $                         CALL DGERV2D( ICTXT, NB, 1, DWORK, NB, 
     $                                       SRSRC, SCSRC ) 
                          POS = POS + 1
                          CALL DLACPY( 'All', NB, 1, DWORK, NB, 
     $                                 EXB(POS), 2 * NB + 1 )
                       END IF
                    END IF
                 END IF
C     
C     Not diagonal block but a block in the upper diagonal part of B
C     
              ELSEIF( J.GT.I ) THEN
C     
C     Check if Bii extended, if so extend with one row from the South
C     
                 IF( EXTINFB( I ).EQ.1 .OR. EXTINFB( I ).EQ.3 ) THEN
                    CALL INFOG2L( (I-1) * NB + IB - IROFFB, 
     $                            (J-1) * NB + IB - IROFFB,
     $                            DESCB, NPROW, NPCOL, MYROW, MYCOL, 
     $                            IX, JX, RSRC, CSRC )
                    CALL INFOG2L( (I-1) * NB + 1, 
     $                            (J-1) * NB + 1,
     $                            DESCSB, NPROW, NPCOL, MYROW, MYCOL, 
     $                            ISX, JSX, RSRC1, CSRC1 )
                    CALL INFOG2L( I * NB + IB - IROFFB, 
     $                            (J-1) * NB + IB - IROFFB, DESCB,
     $                            NPROW, NPCOL, MYROW, MYCOL, IS, JS, 
     $                            SRSRC, SCSRC )
C     
                    IF( (MYROW.EQ.RSRC.AND.MYCOL.EQ.CSRC) .OR.
     $                   (MYROW.EQ.SRSRC.AND.MYCOL.EQ.SCSRC) ) THEN
C     
C     If I own the extra row send it to the holder of Bij
C     
                       LEN = MIN( NB, N + IROFFB - (J-1)*NB )
                       IF( MYROW.EQ.SRSRC.AND.MYCOL.EQ.SCSRC ) THEN
                          CALL DLACPY( 'All', 1, LEN, B((JS-1)*LLDB+IS),
     $                                 LLDB, DWORK, 1 )
                          IF( NPROW.GT.1 )
     $                         CALL DGESD2D( ICTXT, LEN, 1, DWORK, NB,
     $                                       RSRC, CSRC )
                       END IF
C     
C     Receive extra row from process holding the block in the South
C     
                       IF( MYROW.EQ.RSRC.AND.MYCOL.EQ.CSRC ) THEN
                          IF( NPROW.GT.1 )
     $                         CALL DGERV2D( ICTXT, LEN, 1, DWORK, NB,
     $                                       SRSRC, SCSRC )
                          LBR = ICEIL( ISX, NB )
                          LBC = ICEIL( JSX, NB )
                          POS = ( (LBR-1)*NBC + (LBC-1) )*(2 * NB +1) +1 
                          CALL DLACPY( 'All', LEN, 1, DWORK, NB, 
     $                                 EXB(POS), 2 * NB + 1 )
                       END IF
                    END IF
                 END IF
C     
C     Check if Bii and Bjj is extended, if so extend with one extra 
C     corner element from the South-East
C     
                 IF(( EXTINFB( I ).EQ.1 .OR. EXTINFB( I ).EQ.3 ).AND.
     $                ( EXTINFB( J ).EQ.1 .OR. EXTINFB( J ).EQ.3 )) THEN
C     
C     If I own extra corner element send it to the holder of Bij
C     
                    CALL INFOG2L( I * NB + IB - IROFFB, 
     $                            J * NB + IB - IROFFB, 
     $                            DESCB, NPROW, NPCOL, MYROW, MYCOL, IS, 
     $                            JS, SRSRC, SCSRC )
C
                    IF( (MYROW.EQ.RSRC.AND.MYCOL.EQ.CSRC) .OR.
     $                   (MYROW.EQ.SRSRC.AND.MYCOL.EQ.SCSRC) ) THEN
                       IF( MYROW.EQ.SRSRC.AND.MYCOL.EQ.SCSRC ) THEN
                          CORNER = B((JS-1)*LLDB + IS)
                          IF( NPROCS.GT.1 )
     $                         CALL DGESD2D( ICTXT, 1, 1, CORNER, 1, 
     $                                       RSRC, CSRC )
                       END IF
C     
C     Receive extra corner from process in the South-East
C     
                       IF( MYROW.EQ.RSRC.AND.MYCOL.EQ.CSRC ) THEN
                          IF( NPROCS.GT.1 )
     $                         CALL DGERV2D( ICTXT, 1, 1, CORNER, 1,
     $                                       SRSRC, SCSRC )
                          POS = POS + NB 
                          EXB( POS ) = CORNER
                       END IF 
                    END IF
                 END IF
C     
C     Check if Bjj is extended, if so extend with one column from
C     the East
C     
                 IF( EXTINFB( J ).EQ.1 .OR. EXTINFB( J ).EQ.3 ) THEN
                    CALL INFOG2L( (I-1) * NB + IB - IROFFB, 
     $                            (J-1) * NB + IB - IROFFB,
     $                            DESCB, NPROW, NPCOL, MYROW, MYCOL, 
     $                            IX, JX, RSRC, CSRC )
                    CALL INFOG2L( (I-1) * NB + 1, 
     $                            (J-1) * NB + 1,
     $                            DESCSB, NPROW, NPCOL, MYROW, MYCOL, 
     $                            ISX, JSX, RSRC1, CSRC1 )
                    CALL INFOG2L( (I-1) * NB + IB - IROFFB, 
     $                            J * NB + IB - IROFFB, DESCB,
     $                            NPROW, NPCOL, MYROW, MYCOL, IS, JS, 
     $                            SRSRC, SCSRC )
C     
                    IF( (MYROW.EQ.RSRC.AND.MYCOL.EQ.CSRC) .OR.
     $                   (MYROW.EQ.SRSRC.AND.MYCOL.EQ.SCSRC) ) THEN
C     
C     If I own the extra column send it to the holder of the diagonal
C     block
C     
                       IF( MYROW.EQ.SRSRC.AND.MYCOL.EQ.SCSRC ) THEN
                          CALL DLACPY( 'All', NB, 1, B((JS-1)*LLDB+IS),
     $                                 LLDB, DWORK, NB )
                          IF( NPCOL.GT.1 )
     $                         CALL DGESD2D( ICTXT, NB, 1, DWORK, NB, 
     $                                       RSRC,CSRC )
                       END IF
C     
C     Receive extra column from process holding the block in the East
C     
                       IF( MYROW.EQ.RSRC.AND.MYCOL.EQ.CSRC ) THEN
                          IF( NPCOL.GT.1 )
     $                         CALL DGERV2D( ICTXT, NB, 1, DWORK, NB, 
     $                                       SRSRC, SCSRC)
                          LBR = ICEIL( ISX, NB )
                          LBC = ICEIL( JSX, NB )
                          POS = ( (LBR-1)*NBC + (LBC-1) )*( 2 *NB + 1) + 
     $                          1 + NB + 1
                          CALL DLACPY( 'All', NB, 1, DWORK, NB, 
     $                                 EXB(POS), 2 * NB + 1 )
                       END IF
                    END IF
                 END IF
              END IF
 60        CONTINUE
 50     CONTINUE
 40   CONTINUE
C
 65   CONTINUE
C
C     Should we care about C at all? If not, go to end of subroutine
C
      IF( .NOT. C_CARE ) GO TO 95
C     
C     Now extend C in a slightly different way
C     
C     Loop through all blocks in C: treat all blocks equally
C     
C     Now do the 4 sweeps over C
C     
      DO 70 K = 1, 4
C     
C     Select which MBxNB block in each 2x2 block to operate on
C     in this sweep
C     
         IF( K.EQ.1.OR.K.EQ.2 ) ISTART = 1
         IF( K.EQ.3.OR.K.EQ.4 ) ISTART = 2
         IF( K.EQ.1.OR.K.EQ.3 ) JSTART = 1
         IF( K.EQ.2.OR.K.EQ.4 ) JSTART = 2
C     
C     Loop through all those blocks in C
C     
         DO 80 I = ISTART, DBA, 2
            DO 90 J = JSTART, DBB, 2
C     
C     Check if Aii extended, if so extend Cij with one row from
C     the South
C     
               IF( EXTINFA( I ).EQ.1 .OR. EXTINFA( I ).EQ.3 ) THEN
                  CALL INFOG2L( (I-1) * MB + IA - IROFFA, 
     $                          (J-1) * NB + IB - IROFFB, 
     $                          DESCC, NPROW, NPCOL, MYROW, MYCOL, IX, 
     $                          JX, RSRC, CSRC )
                  CALL INFOG2L( (I-1) * MB + 1, 
     $                          (J-1) * NB + 1, 
     $                          DESCSC, NPROW, NPCOL, MYROW, MYCOL, ISX, 
     $                          JSX, RSRC1, CSRC1 )
                  CALL INFOG2L( I * MB + IA - IROFFA, 
     $                          (J-1) * NB + IB - IROFFB, 
     $                          DESCC, NPROW, NPCOL, MYROW, MYCOL, IS, 
     $                          JS, SRSRC, SCSRC )
C
                  IF( (MYROW.EQ.RSRC.AND.MYCOL.EQ.CSRC) .OR.
     $                 (MYROW.EQ.SRSRC.AND.MYCOL.EQ.SCSRC) ) THEN
C     
C     If I own the extra row, send it to the process holding Cij
C     
                     LEN = MIN( NB, N + IROFFB - (J-1)*NB )
                     IF( MYROW.EQ.SRSRC.AND.MYCOL.EQ.SCSRC ) THEN
                        CALL DLACPY( 'All', 1, LEN, C((JS-1)*LLDC+IS),
     $                               LLDC, DWORK, 1 )
                        IF( NPROW.GT.1 ) 
     $                       CALL DGESD2D( ICTXT, LEN, 1, DWORK, NB, 
     $                                     RSRC, CSRC )
                     END IF
C     
C     Receive extra row from the South
C     
                     IF( MYROW.EQ.RSRC.AND.MYCOL.EQ.CSRC ) THEN
                        IF( NPROW.GT.1 )
     $                       CALL DGERV2D( ICTXT, LEN, 1, DWORK, NB, 
     $                                     SRSRC, SCSRC)
                        LBR = ICEIL( ISX, MB )
                        LBC = ICEIL( JSX, NB )
                        POS = ( (LBR-1)*NBC + (LBC-1) )*( MB + NB + 1)+1  
                        CALL DLACPY( 'All', LEN, 1, DWORK, NB, EXC(POS),
     $                               MB + NB + 1 )
                     END IF
                  END IF
               END IF
C     
C     Check if Aii and Bjj is extended, if so extend Cij with an extra
C     corner element from the South-East
C     
               IF( ( EXTINFA( I ).EQ.1 .OR. EXTINFA( I ).EQ.3 ).AND.
     $              ( EXTINFB( J ).EQ.1 .OR. EXTINFB( J ).EQ.3 ) ) THEN   
                  CALL INFOG2L( I * MB + IA - IROFFA, 
     $                          J * NB + IB - IROFFB,
     $                          DESCC, NPROW, NPCOL, MYROW, MYCOL, IS, 
     $                          JS, SRSRC, SCSRC )
C     
C     If I own the extra corner element, send it to the process holding
C     Cij
C     
                  IF( (MYROW.EQ.RSRC.AND.MYCOL.EQ.CSRC) .OR.
     $                 (MYROW.EQ.SRSRC.AND.MYCOL.EQ.SCSRC) ) THEN
                     IF( MYROW.EQ.SRSRC.AND.MYCOL.EQ.SCSRC ) THEN
                        CORNER = C((JS-1)*LLDC + IS) 
                        IF( NPROCS.GT.1 ) 
     $                       CALL DGESD2D( ICTXT, 1, 1, CORNER, 1, RSRC,
     $                                     CSRC )
                     END IF
C     
C     Receive the extra corner element from the South-East
C     
                     IF( MYROW.EQ.RSRC.AND.MYCOL.EQ.CSRC ) THEN
                        IF( NPROCS.GT.1 ) 
     $                       CALL DGERV2D( ICTXT, 1, 1, CORNER, 1, 
     $                                     SRSRC, SCSRC )
                        POS = POS + NB
                        EXC( POS ) = CORNER
                     END IF
                  END IF
               END IF
C     
C     Check if Bjj if extended, then extend Cij with an extra column 
C     from the East
C     
               IF( EXTINFB( J ).EQ.1 .OR. EXTINFB( J ).EQ.3 ) THEN
                  CALL INFOG2L( (I-1) * MB + IA - IROFFA, 
     $                          (J-1) * NB + IB - IROFFB,
     $                          DESCC, NPROW, NPCOL, MYROW, MYCOL, IX, 
     $                          JX, RSRC, CSRC )
                  CALL INFOG2L( (I-1) * MB + 1, 
     $                          (J-1) * NB + 1,
     $                          DESCSC, NPROW, NPCOL, MYROW, MYCOL, ISX, 
     $                          JSX, RSRC1, CSRC1 )
                  CALL INFOG2L( (I-1) * MB + IA - IROFFA, 
     $                          J * NB + IB - IROFFB,
     $                          DESCC, NPROW, NPCOL, MYROW, MYCOL, IS, 
     $                          JS, SRSRC, SCSRC )
C     
                  IF( (MYROW.EQ.RSRC.AND.MYCOL.EQ.CSRC) .OR.
     $                 (MYROW.EQ.SRSRC.AND.MYCOL.EQ.SCSRC) ) THEN
C     
C     If I own the extra column, send it to the process holding Cij
C     
                     LEN = MIN( MB, M + IROFFA - (I-1)*MB )
                     IF( MYROW.EQ.SRSRC.AND.MYCOL.EQ.SCSRC ) THEN
                        CALL DLACPY( 'All', LEN, 1, C((JS-1)*LLDC+IS),
     $                               LLDC, DWORK, MB )
                        IF( NPCOL.GT.1 ) 
     $                       CALL DGESD2D( ICTXT, LEN, 1, DWORK, MB, 
     $                                     RSRC, CSRC )
                     END IF
C     
C     Receive the extra column from the East
C     
                     IF( MYROW.EQ.RSRC.AND.MYCOL.EQ.CSRC ) THEN
                        IF( NPCOL.GT.1 )
     $                       CALL DGERV2D( ICTXT, LEN, 1, DWORK, MB, 
     $                                     SRSRC, SCSRC )
                        LBR = ICEIL( ISX, MB )
                        LBC = ICEIL( JSX, NB )
                        POS =  ( (LBR-1)*NBC + (LBC-1) )*( MB + NB + 1)+ 
     $                         1 + NB + 1 
                        CALL DLACPY( 'All', LEN, 1, DWORK, MB, EXC(POS),
     $                               MB+NB+1)
                     END IF
                  END IF
               END IF
 90         CONTINUE
 80      CONTINUE
 70   CONTINUE
C
 95   CONTINUE
C     
C     Now A, B and C has been implicitely redistributed 
C     
      END
C     
C     End of PDIMPRED
C
C *** Last line of PDIMPRED ***
