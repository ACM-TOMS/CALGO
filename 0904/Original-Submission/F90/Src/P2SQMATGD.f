CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C         This file is a part of the SCASY software library.           C
C         See http://www.cs.umu.se/parallel/research/scasy             C
C         Contributors: Robert Granat and Bo Kågström.                 C
C         SCASY Version 0.10,    March 31, 2006.                       C
C         Copyright 2006, Umeå University, Sweden.                     C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE P2SQMATGD( ADIAG, ASDIAG, BDIAG, UPPER, N, A, DESCA,
     $                      ADW, AMW, ADDIAG, ASUBDIAG, ANQTRBL, ASEED,
     $                      B, DESCB, BDW, BMW, BDDIAG, BSEED, DWORK, 
     $                      LDWORK, INFO )
C
C  -- ScaLAPACK-style routine (preliminary version)--
C     HPC2N and Department of Computing Science,
C     University of Umeå, Sweden.
C     Written by Robert Granat, (granat@cs.umu.se)
C     January 20, 2006.
C
      IMPLICIT NONE
C
C     .. Scalar Arguments ..
      CHARACTER*1        ADIAG, ASDIAG, BDIAG, UPPER
      INTEGER            M, N, LDWORK, INFO, ANQTRBL, ASEED, BSEED
      DOUBLE PRECISION   ADW, AMW, BDW, BMW
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION   A( * ), ADDIAG( * ), ASUBDIAG( * ), B( * ), 
     $                   BDDIAG( * ), DWORK( * )
      INTEGER            DESCA( * ), DESCB( * )
C
C  Purpose
C  =======
C  This subroutine generates the square NxN matrices A and B with real 
C  entries of double precision according to the scheme
C
C     A = P * ( D1 + M1 ) * Q, B = P * ( D2 + M2 ) * Q 
C
C  where P and Q are invertible matrices of order N, D1 and D2 are
C  diagonal matrices of order N, and M1 and M2 are strictly upper 
C  triangular matrices of order N, i.e., with zeros on and below the 
C  diagonal. The values of the entries generated can be scaled by 
C  setting the values of the input parameters _DW (diagonal weight) and
C  _MW (upper triangular weight).
C  If wanted, the diagonals D1 and D2 can be set explicitly by the user 
C  in ADDIAG and BDDIAG. It is also possible to make the upper triangular 
C  matrix ( D1 + M1 ) quasitriangular. In such a case, the user may specify
C  the elements in the subdiagonal ASUBDIAG as well.
C  The matrices P and Q may be chosen as the identity.
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
C  Mode Parameters
C
C  ADIAG    (global input) CHARACTER*1
C           If ADIAG = 'D', the user is supplying the wanted elements
C           in the diagonal matrix D1 in the array ADDIAG.
C           If ADIAG = 'N', the diagonal matrix D1 is generated randomly
C           and the array ADDIAG is not referenced.
C
C  ASDIAG   (global input) CHARACTER*1
C           If ASDIAG = 'D', the upper triangular matrix ( D1 + M1 ) is
C           set as quasitriangular before it is multiplied with the
C           orthogonal matrices P and Q^T. The values on the subdiagonal
C           will be set using ASUBDIAG.
C           If ASDIAG = 'N', and ANQTRBL = 0, the upper triangular matrix 
C           ( D1 + M1 ) will not be set as quasitriangular. 
C           If ASDIAG = 'N' and ANQTRBL >= 1, the matrix ( D1 + M1 ) will 
C           be set as quasitriangular using values from the matrix M1. 
C
C  BDIAG    (global input) CHARACTER*1
C           If BDIAG = 'D', the user is supplying the wanted elements
C           in the diagonal matrix D2 in the array BDDIAG.
C           If BDIAG = 'N', the diagonal matrix D2 is generated randomly
C           and the array BDDIAG is not referenced.
C
C  UPPER    (global input) CHARACTER*1
C           If UPPER = 'U', the othogonal matrices P and Q are chosen 
C           as the identity, and the returned matrices A and B will be
C           upper (quasi)triangular and triangular, respectively.
C           If UPPER = 'N', the orthogonal matrices P and Q will be 
C           generated from QR-factorizing square random matrices.
C
C  Input/Output Arguments
C
C  N        (global input) INTEGER
C           This integer corresponds to the order of the matrices A 
C           and B.
C
C  A        (local output) DOUBLE PRECISION array
C           Array of dimension (LLD_A,LOCc(N)). On output, it contains 
C           the local part of the global distributed generated matrix A.
C           All information contained in A on input will be overwritten.
C
C  DESCA    (global/local input) INTEGER array
C           Array of dimension 9. The array descriptor for the global 
C           distributed matrix A, as defined in ScaLAPACK.
C
C  ADW      (global input) DOUBLE PRECISION
C           ADW sets the number with which each entry on the diagonal
C           matrix D1 shall be weighted ( A(I,I) = A(I,I)*ADW ).
C
C  AMW      (global input) DOUBLE PRECISION
C           AMW sets the number with which each entry in the strictly 
C           upper triangular matrix M1 shall be weighted
C           ( A(I,J) = A(I,J)*AMW, I < J ). 
C
C  ADDIAG   (global input) DOUBLE PRECISION array
C           Array of dimension N. Holds all the wanted diagonal 
C           elements for the matrix D1. Only referenced when ADIAG = 'D'.
C           Notice that this array must be global, i.e., all processes
C           must hold the elements in the whole array. 
C           If ASDIAG = 'D' and ASUBDIAG(I) is nonzero, then ADDIAG(I)
C           and ADDIAG(I+1) must be equal.
C 
C  ASUBDIAG (global input) DOUBLE PRECISION array
C           Array of dimension N-1. Holds all the wanted subdiagonal
C           elements for the matrix (D1 + M1). Only referenced if 
C           ASDIAG = 'D'. If a single entry in ASUBDIAG is nonzero,
C           the preceeding entry and the following entry must be
C           zero, i.e., no neighbouring entries in ASUBDIAG can be
C           nonzero simultaneously.
C
C  ANQTRBL  (global input) INTEGER
C           If SDIAG = 'N', this integer specifies how many 2x2 diagonal
C           blocks that is supposed to be randomly generated on the main 
C           diagonal of the matrix ( D1 + M1 ). May not exceed floor(N/2).
C           The 2x2 blocks generated will be spread out uniformly over
C           the main diagonal. 
C           If ASDIAG = 'N' and ANQTRBL = 0, the matrix ( D1 + M1 ) will be 
C           upper triangular, with no 2x2 blocks on the main diagonal.
C           If ASDIAG = 'D', ANQTRBL will not be referenced.
C
C  ASEED    (global input) INTEGER
C           A seed to the random number generator in PDMATGEN2 for 
C           generating the matrix A.
C
C  B        (local output) DOUBLE PRECISION array
C           Array of dimension (LLD_B,LOCc(N)). On output, it contains 
C           the local part of the global distributed generated matrix B.
C           All information contained in B on input will be overwritten.
C
C  DESCB    (global/local input) INTEGER array
C           Array of dimension 9. The array descriptor for the global 
C           distributed matrix B, as defined in ScaLAPACK.
C
C  BDW      (global input) DOUBLE PRECISION
C           BDW sets the number with which each entry on the diagonal
C           matrix D2 shall be weighted ( B(I,I) = B(I,I)*BDW ).
C
C  BMW      (global input) DOUBLE PRECISION
C           BMW sets the number with which each entry in the strictly 
C           upper triangular matrix M2 shall be weighted 
C           ( B(I,J) = B(I,J)*BMW, I < J).
C
C  BDDIAG   (global input) DOUBLE PRECISION array
C           Array of dimension N. Holds all the wanted diagonal 
C           elements for the matrix D2. Only referenced when BDIAG = 'D'.
C           Notice that this array must be global, i.e., all processes
C           must hold the elements in the whole array.
C           If ASDIAG = 'D' and ASUBDIAG(I) is nonzero, then BDDIAG(I) =
C           BDDIAG(I+1). Notice that BDDIAG(I) = 0 signals an infinite
C           eigenvalue.
C
C  BSEED    (global input) INTEGER
C           A seed to the random number generator in PDMATGEN2 for 
C           generating the matrix B.
C
C  Workspace
C        
C  DWORK    (local input) DOUBLE PRECISION array
C           Array of dimension LDWORK. Used as workspace througout
C           the generation process.
C
C  LDWORK   (local input) INTEGER
C           The dimension of the array DWORK. The optimal value
C           can be found by setting LDWORK = -1. Then this routine
C           will store the optimal value of LDWORK in DWORK(1) and
C           return immediately. No error will then be signaled by
C           PXERBLA.
C
C  Error indicator
C
C  INFO     (local output) INTEGER
C            = 0:  successful exit
C            < 0:  unsuccessful exit. 
C            If the i-th argument is an array and the j-th entry had
C            an illegal value, then INFO = -(i*100+j), if the i-th
C            argument is a scalar and had an illegal value, then
C            INFO = -i.
C            = 1: There was no valid BLACS context.
C
C  Method
C  ======
C  The two matrices are formed by building upper triangular matrices,
C  where each diagonal is weighted with the same factor, and the off-
C  diagonal entries are weighted using another factor. Then the two 
C  matrices are multiplied with two orthogonal transformation matrix P
C  and Q which are formed by QR-factorizing two random matrices. For 
C  details, see the references.
C
C  Additional requirements
C  =======================
C  The global distributed matrices A and B must fulfill the following
C  alignment requirements: DESCA( MB_ ) = DESCB( MB_ ), DESCA( NB_ ) =
C  DESCB( NB_ ),  DESCA( RSRC_ ) = DESCB( RSRC_ ), DESCA( CRSC_ ) =
C  DESCB( CSRC_ ). 
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
C  Specified eigenvalues, random matrices, quasitriangular matrices,
C  matrix pairs.
C
C  *********************************************************************
C
C     .. Parameters ..
      INTEGER            BLOCK_CYCLIC_2D, CSRC_, CTXT_, DLEN_, DTYPE_,
     $                   LLD_, MB_, M_, NB_, N_, RSRC_, COND
      DOUBLE PRECISION   ZERO, MINONE, ONE
      PARAMETER          ( BLOCK_CYCLIC_2D = 1, DLEN_ = 9, DTYPE_ = 1,
     $                     CTXT_ = 2, M_ = 3, N_ = 4, MB_ = 5, NB_ = 6,
     $                     RSRC_ = 7, CSRC_ = 8, LLD_ = 9, 
     $                     ZERO = 0.0D+0, MINONE = -1.0D+0,
     $                     ONE = 1.0D+00, COND = 10 )         
C     .. 
C     .. Local Scalars ..
      INTEGER            ICTXT, NPROW, NPCOL, MYROW, MYCOL, IDUM1,
     $                   IDUM2, AROWS, ACOLS, ASIZE, WORK, SWORK,
     $                   TAUP, TAUQ, LLDA, LLDB, MB, NB, I, J, GI, GJ, 
     $                   STEP, LI1, LJ1, RSRC1, CSRC1, LI2, LJ2, RSRC2, 
     $                   CSRC2, LI3, LJ3, RSRC3, CSRC3, LI4, LJ4, K,
     $                   RSRC4, CSRC4, P, Q, QRWORK, BROWS, BCOLS,
     $                   NPROCS
      LOGICAL            GVADIAG, GVBDIAG, AQUASITRI, GVASDIAG, UPPERTRI
      DOUBLE PRECISION   ELEM
C     ..      
C     .. Local Arrays ..
      INTEGER            DESCP(DLEN_), DESCQ( DLEN_ ), ISEED(4) 
      DOUBLE PRECISION   DPDUM1(1), DPDUM2(1)
C     ..     
C     .. External Subroutines ..
      EXTERNAL           PDMATGEN2, PXERBLA, BLACS_GRIDINFO, BLACS_PNUM,
     $                   CHK1MAT, PCHK1MAT, PDGEQRF, PDORMQR, DGESD2D,
     $                   DGERV2D, DESCINIT, INFOG2L
C     ..
C     .. External Functions ..
      INTEGER            NUMROC, INDXL2G, ICEIL
      LOGICAL            LSAME
      DOUBLE PRECISION   DLARAN
      EXTERNAL           NUMROC, LSAME, INDXL2G, ICEIL
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
C     Check if NPROW and NPCOL are O.K. etc.
C     
      IF( INFO.EQ.0 ) THEN
         CALL CHK1MAT( N, 5, N, 5, 1, 1, DESCA, 7, INFO )
      END IF
C     
C     Check the blocking sizes for being not to large
C     
      IF( INFO.EQ.0 ) THEN
         IF( DESCA( MB_ ).GT.N ) INFO = -(100*7 + MB_)
      END IF
C     
C     Make sure all global variables associated with the global matrix A
C     are indeed global and that INFO is really set to the same value at 
C     all the processes
C     
      IF( INFO.EQ.0 ) THEN
         CALL PCHK1MAT( N, 5, N, 5, 1, 1, DESCA, 7, 0, IDUM1, 
     $        IDUM2, INFO )
      END IF
C
C     Check alignment requirements
C
      IF( INFO.EQ. 0 ) THEN
         IF( DESCA( MB_ ).NE.DESCB( NB_ ) ) INFO = -14*100 + MB_
         IF( DESCA( NB_ ).NE.DESCB( NB_ ) ) INFO = -14*100 + NB_
         IF( DESCA( RSRC_ ).NE.DESCB( RSRC_ ) ) INFO = -14*100 + MB_
         IF( DESCA( CSRC_ ).NE.DESCB( CSRC_ ) ) INFO = -14*100 + NB_
      END IF
C     
C     Check the values of the mode parameters
C     
      IF( INFO.EQ.0 ) THEN
         IF( .NOT.( LSAME( ADIAG, 'D' ) .OR. LSAME( ADIAG, 'N' ))) THEN
            INFO = -1
         END IF
      END IF
C     
      IF( INFO.EQ.0 ) THEN
         IF( .NOT.( LSAME( ASDIAG, 'D' ) .OR. LSAME( ASDIAG, 'N' )))THEN
            INFO = -2
         END IF
      END IF
CC     
      IF( INFO.EQ.0 ) THEN
         IF( .NOT.( LSAME( BDIAG, 'D' ) .OR. LSAME( BDIAG, 'N' ))) THEN
            INFO = -3
         END IF
      END IF
C
      IF( INFO.EQ.0 ) THEN
         IF( .NOT.( LSAME( UPPER, 'U' ) .OR. LSAME( UPPER, 'N' ))) THEN
            INFO = -4
         END IF
      END IF
C     
C     Check that the values of ADDIAG are legal
C     
      IF( INFO.EQ.0 .AND. LSAME( ADIAG, 'D' ) ) THEN
         DO 5 I = 1, N-1
            IF( LSAME(ASDIAG, 'D') .AND. ASUBDIAG(I).NE.ZERO ) THEN
               IF( ADDIAG(I).NE.ADDIAG(I+1) ) INFO = -(10*100 + I + 1)
            END IF  
 5       CONTINUE
      END IF
C     
C     Check that the values of ASUBDIAG are legal
C     
      IF( INFO.EQ.0 .AND. LSAME(ASDIAG, 'D') ) THEN
         IF( ASUBDIAG(1).NE.ZERO .AND. ASUBDIAG(2).NE.ZERO) THEN
            INFO = -(11*100 + 2)
         END IF
         DO 7 I = 2, N-2
            IF( ASUBDIAG(I).NE.ZERO .AND. ( ASUBDIAG(I-1).NE.ZERO .OR.
     $           ASUBDIAG(I+1).NE.ZERO )) INFO = -(11*100 + I)
 7       CONTINUE
         IF( ASUBDIAG(N-1).NE.ZERO .AND. ASUBDIAG(N-2).NE.ZERO ) THEN
            INFO = -(11*100+N-2)
         END IF
      END IF
C     
C     Check the value of ANQTRBL
C     
      IF( INFO.EQ.0 .AND. LSAME( ASDIAG, 'N') ) THEN
         IF( ANQTRBL.GT.(N/2) ) INFO = -12
      END IF
C     
C     Set some logical values to decide how to proceed
C     
      IF( INFO.EQ.0 ) THEN
         GVADIAG = LSAME( ADIAG, 'D' )
         GVASDIAG = LSAME( ASDIAG, 'D' )
         GVBDIAG = LSAME( BDIAG, 'D' )
         IF( .NOT. GVASDIAG ) THEN
            AQUASITRI = ANQTRBL.GT.0
         ELSE
            AQUASITRI = .FALSE.
         END IF
         UPPERTRI = LSAME( UPPER, 'U' )
      END IF
C     
C     Check workspace
C     
      IF( INFO.EQ.0 ) THEN
         WORK = 0
         IF( .NOT. UPPERTRI ) THEN
            AROWS = NUMROC( N, DESCA(MB_), MYROW, DESCA(RSRC_), NPROW )
            ACOLS = NUMROC( N, DESCA(NB_), MYCOL, DESCA(CSRC_), NPCOL )
            BROWS = AROWS
            BCOLS = ACOLS
            ASIZE = AROWS*ACOLS
            WORK = WORK + 2*ASIZE
C     
C     Set a descriptor for the matrices P and Q
C     
            CALL DESCINIT( DESCP, N, N, DESCA(MB_), DESCA(NB_), 
     $           DESCA(RSRC_), DESCA(CSRC_), ICTXT, MAX(1,AROWS),
     $           INFO )
C     
            CALL DESCINIT( DESCQ, N, N, DESCB(MB_), DESCB(NB_), 
     $           DESCB(RSRC_), DESCB(CSRC_), ICTXT, MAX(1,AROWS),
     $           INFO )
C     
C     Compute the needed workspace for the QR-factorization using
C     workspace queries
C     
            CALL PDGEQRF( N, N, A, 1, 1, DESCA, DPDUM1, DPDUM2, -1, 
     $           INFO )
            QRWORK = INT(DPDUM2(1))
            CALL PDORMQR( 'Left', 'NoTranspose', N, N, N, DWORK(1), 1, 
     $           1, DESCQ, DPDUM1, A, 1, 1, DESCA, DPDUM2, -1, INFO )
            QRWORK = MAX( QRWORK, INT(DPDUM2(1)))
            CALL PDORMQR( 'Right', 'Transpose', N, N, N, DWORK(1), 1, 
     $           1, DESCQ, DPDUM1, A, 1, 1, DESCA, DPDUM2, -1, INFO )
            QRWORK = MAX( QRWORK, INT(DPDUM2(1)))
            WORK = WORK + 2*ACOLS + QRWORK 
         END IF
      END IF
C     
      IF( INFO.EQ.0 ) THEN
         IF( WORK.GT.LDWORK.AND.LDWORK.NE.-1 ) INFO = -21
         IF( LDWORK.EQ.-1 ) THEN
            DWORK(1) = WORK
            RETURN
         END IF
      END IF
C     
C     Checking if we may continue or should interrupt
C     
      IF( INFO.NE.0 ) THEN
         CALL PXERBLA( ICTXT, 'P2SQMATGD', -INFO )
         RETURN
      END IF
C     
C     Quick return if possible
C     
      IF( N.EQ.0 ) RETURN
C     
C     Extract the blocking sizes and the leading local dimension of A
C     
      MB = DESCA( MB_ )
      NB = DESCA( NB_ )
      LLDA = DESCA( LLD_ )
      LLDB = DESCB( LLD_ )
C
C     Set seeds for random generator
C
      ISEED( 1 ) = NPCOL * 7
      ISEED( 2 ) = NPROW * 13
      ISEED( 3 ) = MOD( ASEED, 1234 )
      ISEED( 4 ) = MOD( BSEED, 4321 )
C     
C     Divide the workspace between the different integer pointers
C     for ScaLAPACK arrays. 
C     
      IF( .NOT. UPPERTRI ) THEN
         P      = 1
         Q      = P + ASIZE
         TAUP   = Q + ASIZE
         TAUQ   = TAUP + ACOLS
         SWORK  = TAUQ + ACOLS
      ELSE
         SWORK = 1
      END IF
C     
C     Start the real matrix generating work
C
      CALL PDLASET( 'All', N, N, 0D0, 0D0, A, 1, 1, DESCA )
      CALL PDLASET( 'All', N, N, 0D0, 0D0, B, 1, 1, DESCB )
C     
C     Now init a random upper triangular matrix in A using PDMATGEN2:
C     now we have a preliminary version of ( D1 + M1 ). 
C     
      IF( UPPERTRI ) THEN
         CALL P1SQMATGD( ADIAG, ASDIAG, 'Upper', N, A, DESCA, ADW, AMW, 
     $                   ADDIAG, ASUBDIAG, ANQTRBL, ASEED, DWORK, 
     $                   LDWORK, INFO )
         CALL P1SQMATGD( BDIAG, 'No subdiagonal', 'Upper', N, B, DESCB, 
     $                   BDW, BMW, BDDIAG, ASUBDIAG, 0, BSEED, DWORK, 
     $                   LDWORK, INFO )
         RETURN
      END IF
C
C     If a user-specified diagonal matrix D1 is present:
C     Set the new diagonal in the matrix D1 and weight it with ADW.
C     At the same time, weight the matrix M1 with AMW. If no user 
C     specified diagonal matrix D1 is present, just scale the upper 
C     triangular matrix ( D1 + M1 ) using the given values of ADW, AMW.
C
      DO 10 J = 1, ACOLS
         DO 20 I = 1, AROWS
            GJ = INDXL2G( J, NB, MYCOL, DESCA( CSRC_ ), NPCOL )
            GI = INDXL2G( I, MB, MYROW, DESCA( RSRC_ ), NPROW )
            IF( GI.EQ.GJ ) THEN
               IF( GVADIAG ) THEN
                  A( (J-1)*LLDA + I ) = ADW * ADDIAG( GI )
               ELSE
                  A( (J-1)*LLDA + I ) = ADW * A( (J-1)*LLDA + I ) 
               END IF
            ELSEIF( GI.LT.GJ ) THEN
               A( (J-1)*LLDA + I ) = AMW * A( (J-1)*LLDA + I ) 
            END IF
 20      CONTINUE
 10   CONTINUE
C
C     Do the same thing for B...
C
       DO 30 J = 1, BCOLS
         DO 40 I = 1, BROWS
            GJ = INDXL2G( J, NB, MYCOL, DESCB( CSRC_ ), NPCOL )
            GI = INDXL2G( I, MB, MYROW, DESCB( RSRC_ ), NPROW )
            IF( GI.EQ.GJ ) THEN
               IF( GVBDIAG ) THEN
                  B( (J-1)*LLDB + I ) = BDW * BDDIAG( GI )
               ELSE
                  B( (J-1)*LLDB + I ) = BDW * B( (J-1)*LLDB + I ) 
               END IF
            ELSEIF( GI.LT.GJ ) THEN
               B( (J-1)*LLDB + I ) = BMW * B( (J-1)*LLDB + I ) 
            END IF
 40      CONTINUE
 30   CONTINUE
C     
C     Set the subdiagonal elements in A - form a quasitriangular matrix
C     out of the matrix ( D1 + M1 ). Imagine the entries in the 2x2 blocks
C     being labeled as 
C     
C                   |-------------|
C                   |  a   |  b   |
C                   |      |      |
C                   ---------------
C                   | -b   |  a   |
C                   |      |      |
C                   ---------------
C
C     If no user-specified subdiagonal is present and we still want
C     ( D1 + M1 ) to be quasitriangular we set the 2x2 blocks using
C     the corresponing AMW-scaled (b)-elements in M1. 
C
      IF ( AQUASITRI ) THEN
C     
C     Loop through the matrix and create 2x2 blocks on the diagonal
C     
         STEP = ICEIL( N, ANQTRBL)
         DO 50 GI = 1, N-1, STEP
C     
C     Check who owns the four entries we are going to operate on
     
            CALL INFOG2L( GI, GI, DESCA, NPROW, NPCOL, MYROW, MYCOL,
     $                    LI1, LJ1, RSRC1, CSRC1 )
            CALL INFOG2L( GI, GI+1, DESCA, NPROW, NPCOL, MYROW, MYCOL,
     $                    LI2, LJ2, RSRC2, CSRC2 )
            CALL INFOG2L( GI+1, GI, DESCA, NPROW, NPCOL, MYROW, MYCOL,
     $                    LI3, LJ3, RSRC3, CSRC3 )
            CALL INFOG2L( GI+1, GI+1, DESCA, NPROW, NPCOL, MYROW, MYCOL,
     $                    LI4, LJ4, RSRC4, CSRC4 )
C     
C     Setting the bottom-right element in the 2x2 block to the same
C     value as the upper-left (a). If the processors are not the same
C     we have to communicate the value.
C     
            IF( RSRC1.NE.RSRC4 .OR. CSRC1.NE.CSRC4 ) THEN
               IF( MYROW.EQ.RSRC1 .AND. MYCOL.EQ.CSRC1 ) THEN
                  CALL DGESD2D( ICTXT, 1, 1, A((LJ1-1)*LLDA + LI1),
     $                          1, RSRC4, CSRC4 )
               ELSEIF( MYROW.EQ.RSRC4 .AND. MYCOL.EQ.CSRC4 ) THEN
                  CALL DGERV2D( ICTXT, 1, 1, A((LJ4-1)*LLDA + LI4),
     $                          1, RSRC1, CSRC1 )
               END IF
            ELSEIF( MYROW.EQ.RSRC1 .AND. MYCOL.EQ.RSRC1 ) THEN
               A((LJ4-1)*LLDA + LI4) = A((LJ1-1)*LLDA + LI1)
            END IF
C     
C     Setting the bottom-left value in the 2x2 block to the negated
C     value of the upper-right entry (b). If the processors are not
C     the same we have to communicate the value.
C     
            IF( RSRC2.NE.RSRC3 .OR. CSRC2.NE.CSRC3 ) THEN
               IF( MYROW.EQ.RSRC2 .AND. MYCOL.EQ.CSRC2 ) THEN
                  IF( .NOT. UPPERTRI ) THEN
                     A((LJ2-1)*LLDA + LI2) = DLARAN( ISEED )  
                  END IF
                  CALL DGESD2D( ICTXT, 1, 1, A((LJ2-1)*LLDA + LI2),
     $                          1, RSRC3, CSRC3 )
               ELSEIF( MYROW.EQ.RSRC3 .AND. MYCOL.EQ.CSRC3 ) THEN
                  CALL DGERV2D( ICTXT, 1, 1, ELEM, 1, RSRC2, CSRC2 )
                 A((LJ3-1)*LLDA + LI3) = MINONE * ELEM
               END IF
            ELSEIF( MYROW.EQ.RSRC2 .AND. MYCOL.EQ.RSRC2 ) THEN
               IF( .NOT. UPPERTRI ) THEN
                  A((LJ2-1)*LLDA + LI2) = DLARAN( ISEED )  
               END IF
               A((LJ3-1)*LLDA + LI3) = MINONE*A((LJ2-1)*LLDA + LI2)
            END IF 
C
C     Put zeros in the superdiagonal of B corresponding to the 2-by-2
C     block in A to produce the correct eigenvalues
C
            CALL PDELSET( B, GI, GI+1, DESCB, ZERO )
C
 50      CONTINUE
C
C     If the user has specified a subdiagonal, set the corresponding
C     entries in ( D1 + M1 ). At the same time, scale the new b-elements.
C     using ADW.
C
      ELSEIF( GVASDIAG ) THEN 
         DO 60 GI = 1, N-1, 1
            IF( ASUBDIAG(GI).NE.ZERO ) THEN
               CALL INFOG2L( GI, GI+1, DESCA, NPROW, NPCOL, MYROW, 
     $                       MYCOL, LI2, LJ2, RSRC2, CSRC2 )
               CALL INFOG2L( GI+1, GI, DESCA, NPROW, NPCOL, MYROW, 
     $                       MYCOL, LI3, LJ3, RSRC3, CSRC3 )
               IF( MYROW.EQ.RSRC3 .AND. MYCOL.EQ.CSRC3 ) THEN
                  A((LJ3-1)*LLDA + LI3) = ADW * ASUBDIAG(GI)
               END IF
               IF( MYROW.EQ.RSRC2 .AND. MYCOL.EQ.CSRC2 ) THEN
                  A((LJ2-1)*LLDA + LI2) = MINONE * ADW * ASUBDIAG(GI)
               END IF
C
C     Put zeros in the superdiagonal of B corresponding to the 2-by-2
C     block in A to produce the correct eigenvalues
C
               CALL PDELSET( B, GI, GI+1, DESCB, ZERO )
C
            END IF
 60      CONTINUE
      END IF
C
C     If (A,B) should be general (not triangular) perform an equivalence
C     transformation with invertible matrices with condition number 
C     equal to N. Generate these invertible matrices implicitly via
C     a factorized SVD-form: X = U * S * V', where S = diag(1:n) and
C     U and V are different orthogonal matrices. 
C
      IF( .NOT. UPPERTRI ) THEN
C     
         CALL PDMATGEN2( ICTXT, 'Random', 'NoDiagDominant', N, N, MB, 
     $        NB, DWORK(P), DESCP(LLD_), 0, 0, ASEED+7, 0, AROWS, 0, 
     $        ACOLS, MYROW, MYCOL, NPROW, NPCOL )
C         
         CALL PDMATGEN2( ICTXT, 'Random', 'NoDiagDominant', N, N, MB, 
     $        NB, DWORK(Q), DESCQ(LLD_), 0, 0, BSEED+7, 0, AROWS, 0, 
     $        ACOLS, MYROW, MYCOL, NPROW, NPCOL )
C     
         CALL PDGEQRF( N, N, DWORK(P), 1, 1, DESCP, DWORK(TAUP), 
     $        DWORK(SWORK), QRWORK, INFO )
C
         CALL PDGEQRF( N, N, DWORK(Q), 1, 1, DESCQ, DWORK(TAUQ), 
     $        DWORK(SWORK), QRWORK, INFO )
C     
         CALL PDORMQR( 'Left', 'Transpose', N, N, N, DWORK(P), 1, 1,
     $        DESCP, DWORK(TAUP), A, 1, 1, DESCA, DWORK(SWORK),
     $        QRWORK, INFO )
C
         CALL PDORMQR( 'Left', 'Transpose', N, N, N, DWORK(P), 1, 1,
     $        DESCP, DWORK(TAUP), B, 1, 1, DESCB, DWORK(SWORK),
     $        QRWORK, INFO )
C         
         CALL PDORMQR( 'Right', 'NoTranspose', N, N, N, DWORK(Q), 1, 1,
     $        DESCQ, DWORK(TAUQ), A, 1, 1, DESCA, DWORK(SWORK),
     $        QRWORK, INFO )
C
         CALL PDORMQR( 'Right', 'NoTranspose', N, N, N, DWORK(Q), 1, 1,
     $        DESCQ, DWORK(TAUQ), B, 1, 1, DESCB, DWORK(SWORK),
     $        QRWORK, INFO )
C     
         DO 70 K = 1, N
            CALL PDSCAL( N, DBLE(MAX(1,MOD(K,COND))), A, 1, K, DESCA, 1)
            CALL PDSCAL( N, DBLE(MAX(1,MOD(K,COND))), B, 1, K, DESCB, 1)
            CALL PDSCAL( N, ONE / DBLE(MAX(1,MOD(K,COND))), A, K, 1, 
     $                   DESCA, DESCA(M_) )
            CALL PDSCAL( N, ONE / DBLE(MAX(1,MOD(K,COND))), B, K, 1, 
     $                   DESCB, DESCB(N_) )
 70      CONTINUE
C     
         CALL PDMATGEN2( ICTXT, 'Random', 'NoDiagDominant', N, N, MB, 
     $        NB, DWORK(P), DESCP(LLD_), 0, 0, ASEED+13, 0, AROWS, 0, 
     $        ACOLS, MYROW, MYCOL, NPROW, NPCOL )
C     
         CALL PDMATGEN2( ICTXT, 'Random', 'NoDiagDominant', N, N, MB, 
     $        NB, DWORK(Q), DESCQ(LLD_), 0, 0, BSEED+13, 0, AROWS, 0, 
     $        ACOLS, MYROW, MYCOL, NPROW, NPCOL )
C     
         CALL PDGEQRF( N, N, DWORK(P), 1, 1, DESCP, DWORK(TAUP), 
     $        DWORK(SWORK), QRWORK, INFO )
C     
         CALL PDGEQRF( N, N, DWORK(Q), 1, 1, DESCQ, DWORK(TAUQ), 
     $        DWORK(SWORK), QRWORK, INFO )
C     
         CALL PDORMQR( 'Left', 'NoTranspose', N, N, N, DWORK(P), 1, 1,
     $        DESCP, DWORK(TAUP), A, 1, 1, DESCA, DWORK(SWORK),
     $        QRWORK, INFO )
C
         CALL PDORMQR( 'Left', 'NoTranspose', N, N, N, DWORK(P), 1, 1,
     $        DESCP, DWORK(TAUP), B, 1, 1, DESCB, DWORK(SWORK),
     $        QRWORK, INFO )
C         
         CALL PDORMQR( 'Right', 'Transpose', N, N, N, DWORK(Q), 1, 1,
     $        DESCQ, DWORK(TAUQ), A, 1, 1, DESCA, DWORK(SWORK),
     $        QRWORK, INFO )
C
         CALL PDORMQR( 'Right', 'Transpose', N, N, N, DWORK(Q), 1, 1,
     $        DESCQ, DWORK(TAUQ), B, 1, 1, DESCB, DWORK(SWORK),
     $        QRWORK, INFO )
      END IF
C     
      END
C     
C     End of P2SQMATGD
C     
C *** Last line of P2SQMATGD ***
