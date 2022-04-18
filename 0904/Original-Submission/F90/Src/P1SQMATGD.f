CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C         This file is a part of the SCASY software library.           C
C         See http://www.cs.umu.se/parallel/research/scasy             C
C         Contributors: Robert Granat and Bo Kågström.                 C
C         SCASY Version 0.10,    March 31, 2006.                       C
C         Copyright 2006, Umeå University, Sweden.                     C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE P1SQMATGD( DIAG, SDIAG, UPPER, N, A, DESCA, DW, MW, 
     $                      DDIAG, SUBDIAG, NQTRBL, ASEED, DWORK, 
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
      CHARACTER*1        DIAG, SDIAG, UPPER
      INTEGER            N, LDWORK, INFO, NQTRBL, ASEED
      DOUBLE PRECISION   DW, MW
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION   A( * ), DDIAG( * ), SUBDIAG( * ), DWORK( * )
      INTEGER            DESCA( * )
C
C  Purpose
C  =======
C  This subroutine generates a square NxN matrix A with real 
C  entries of double precision according to the scheme
C
C     A = Q * ( D + M ) * Q^T,
C
C  where Q is an orthogonal matrix of order N, D is a diagonal matrix
C  of order N, and M is a stricly triangular matrix of order N, i.e.,
C  with zeros on and below the diagonal. The values of the entries 
C  generated can be scaled by setting the values of the input parameters 
C  DW (diagonal weight) and MW (upper triangular weight).
C  If wanted, the diagonal D can be set explicitely by the user in 
C  DDIAG. It is also possible to make the upper triangular matrix
C  ( D + M ) quasitriangular. In such a case the user may specify
C  the elements in the subdiagonal as well.
C  The orthogonal matrix Q can be chosen as the identity matrix.
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
C  DIAG     (global input) CHARACTER*1
C           If DIAG = 'D', the user is supplying the wanted elements
C           in the diagonal matrix D in the array DDIAG.
C           If DIAG = 'N', the diagonal matrix D is generated randomly
C           and the array DDIAG is not referenced.
C
C  SDIAG    (global input) CHARACTER*1
C           If SDIAG = 'D', the upper triangular matrix ( D + M ) is
C           set as quasitriangular before it is multiplied with the
C           orthogonal matrices Q and Q^T. The values on the subdiagonal
C           will be set using provided values from the array SUBDIAG.
C           If SDIAG = 'N', and NQTRBL = 0, the upper triangular matrix 
C           ( D + M ) will not be set as quasitriangular. 
C           If SDIAG = 'N' and NQTRBL >= 1, the matrix ( D + M ) will be 
C           set as quasitriangular using values from the matrix M. 
C
C  UPPER    (global input) CHARACTER*1
C           If UPPER = 'U', the othogonal matrix Q is chosen as
C           the identity matrix, and the returned matrix A will be
C           upper triangular.
C           If UPPER = 'N', the orthogonal matrix Q will be generated
C           from QR-factorizing a square random matrix.
C
C  Input/Output Arguments
C
C  N        (global input) INTEGER
C           This integer corresponds to the order of the matrix A.
C
C  DW       (global input) DOUBLE PRECISION
C           DW sets the number with which each entry on the diagonal
C           matrix D shall be weighted.
C
C  MW       (global input) DOUBLE PRECISION
C           MW sets the number with which each entry in the striclty 
C           upper triangular matrix M shall be weighted.
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
C  DDIAG    (global input) DOUBLE PRECISION array
C           Array of dimension N. Holds all the wanted diagonal 
C           elements for the matrix D. Only referenced when DIAG = 'D'.
C           Notice that this array must be global, i.e., all processes
C           must hold the elements in the whole array. If SDIAG = N and
C           NQTRBL > 0, some of the values in DDIAG will overwrite 
C           others on the diagonal of D when generating 2x2 diagonal 
c           blocks.
C           If SDIAG = 'D' and SUBDIAG(I) is nonzero, then DDIAG(I)
C           and DDIAG(I+1) must be equal.
C 
C  SUBDIAG  (global input) DOUBLE PRECISION array
C           Array of dimension N-1. Holds all the wanted subdiagonal
C           elements for the matrix (D + M). Only referenced if 
C           SDIAG = 'D'. If a single entry in SUBDIAG is nonzero,
C           the preceeding entry and the following entry must be
C           zero, i.e., no neighbouring entries in SUBDIAG may be
C           nonzero simultaneoulsy.
C
C  NQTRBL   (global input) INTEGER
C           If SDIAG = 'N', this integer specifies how many 2x2 diagonal
C           blocks that is supposed to be randomly generated on the main 
C           diagonal of the matrix ( D + M ). May not exceed floor(N/2).
C           The 2x2 blocks generated will be spread out uniformly over
C           the main diagonal. 
C           If SDIAG = 'N' and NQTRBL = 0, the matrix ( D + M ) will be 
C           upper triangular, with no 2x2 blocks on the main diagonal.
C           If SDIAG = 'D', NQTRBL will not be referenced.
C
C  ASEED    (global input) INTEGER
C           A seed to the random number generator in PDMATGEN2.
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
C            = 1: There was no vaild BLACS context.
C
C  Method
C  ======
C  The matrix is formed by building an upper (quasi-)triangular matrix,
C  where the diagonal is weighted with the same factor, and where the off-
C  diagonal entries are weighted using another factor. Then this matrix
C  is multiplied with an orthogonal similarity transformation matrix Q,
C  which is formed by QR-factorizing a third random matrix. For details, 
C  see the references.
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
C  Specified eigenvalues, random matrices, quasitriangular matrices.
C
C  *********************************************************************
C
C     .. Parameters ..
      INTEGER            BLOCK_CYCLIC_2D, CSRC_, CTXT_, DLEN_, DTYPE_,
     $                   LLD_, MB_, M_, NB_, N_, RSRC_
      DOUBLE PRECISION   ZERO, MINONE
      PARAMETER          ( BLOCK_CYCLIC_2D = 1, DLEN_ = 9, DTYPE_ = 1,
     $                     CTXT_ = 2, M_ = 3, N_ = 4, MB_ = 5, NB_ = 6,
     $                     RSRC_ = 7, CSRC_ = 8, LLD_ = 9, 
     $                     ZERO = 0.0D+0, MINONE = -1.0D+0 )         
C     .. 
C     .. Local Scalars ..
      INTEGER            ICTXT, NPROW, NPCOL, MYROW, MYCOL, IDUM1,
     $                   IDUM2, AROWS, ACOLS, ASIZE, WORK, SWORK,
     $                   TAU, LLDA, MB, NB, I, J, GI, GJ, STEP,
     $                   LI1, LJ1, RSRC1, CSRC1, LI2, LJ2, RSRC2, 
     $                   CSRC2, LI3, LJ3, RSRC3, CSRC3, LI4, LJ4, 
     $                   RSRC4, CSRC4, Q, QRWORK
      LOGICAL            GVDIAG, QUASITRI, GVSDIAG, UPPERTRI
      DOUBLE PRECISION   ELEM
C     ..      
C     .. Local Arrays ..
      INTEGER            DESCQ( DLEN_ ) 
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
         CALL CHK1MAT( N, 4, N, 4, 1, 1, DESCA, 6, INFO )
      END IF
C     
C     Check the blocking sizes for being not to large
C     
      IF( INFO.EQ.0 ) THEN
         IF( DESCA( MB_ ).GT.N ) INFO = -(100*6 + MB_)
      END IF
C     
C     Make sure all global variables associated with the global matrix A
C     are indeed global and that INFO is really set to the same value at 
C     all the processes
C     
      IF( INFO.EQ.0 ) THEN
         CALL PCHK1MAT( N, 4, N, 4, 1, 1, DESCA, 6, 0, IDUM1, 
     $                  IDUM2, INFO )
      END IF
C     
C     Check the values of the mode parameters
C     
      IF( INFO.EQ.0 ) THEN
         IF( .NOT.( LSAME( DIAG, 'D' ) .OR. LSAME( DIAG, 'N' ))) THEN
            INFO = -1
         END IF
      END IF
C     
      IF( INFO.EQ.0 ) THEN
         IF( .NOT.( LSAME( SDIAG, 'D' ) .OR. LSAME( SDIAG, 'N' ))) THEN
            INFO = -2
         END IF
      END IF
C
      IF( INFO.EQ.0 ) THEN
         IF( .NOT.( LSAME( UPPER, 'U' ) .OR. LSAME( UPPER, 'N' ))) THEN
            INFO = -3
         END IF
      END IF
C     
C     Check that the values of DDIAG are legal
C     
      IF( INFO.EQ.0 .AND. LSAME( DIAG, 'D' ) ) THEN
         DO 5 I = 1, N-1
            IF( LSAME( SDIAG, 'D' ) .AND. SUBDIAG( I ).NE.ZERO ) THEN
               IF( DDIAG(I).NE.DDIAG(I+1) ) INFO = -(9*100 + I + 1)
            END IF  
 5       CONTINUE
      END IF
C     
C     Check that the values of SUBDIAG are legal
C     
      IF( INFO.EQ.0 .AND. LSAME(SDIAG, 'D') ) THEN
         IF( SUBDIAG(1).NE.ZERO .AND. SUBDIAG(2).NE.ZERO) THEN
            INFO = -(10*100 + 2)
         END IF
         DO 7 I = 2, N-2
            IF( SUBDIAG(I).NE.ZERO .AND. ( SUBDIAG(I-1).NE.ZERO .OR.
     $           SUBDIAG(I+1).NE.ZERO )) INFO = -(10*100 + I)
 7       CONTINUE
         IF( SUBDIAG(N-1).NE.ZERO .AND. SUBDIAG(N-2).NE.ZERO ) THEN
            INFO = -(10*100+N-2)
         END IF
      END IF
C     
C     Check the value of NQTRBL
C     
      IF( INFO.EQ.0 .AND. LSAME( SDIAG, 'N') ) THEN
         IF( NQTRBL.GT.(N/2) ) INFO = -11
      END IF
C     
C     Check workspace
C      
      IF( INFO.EQ.0 ) THEN
         WORK = 0
C     
         AROWS = NUMROC( N, DESCA(MB_), MYROW, DESCA(RSRC_), NPROW )
         ACOLS = NUMROC( N, DESCA(NB_), MYCOL, DESCA(CSRC_), NPCOL )
         ASIZE = AROWS*ACOLS
         IF( .NOT. LSAME(UPPER, 'U')) WORK = WORK + ASIZE
C     
C     Set a descriptor for the matrix Q
C     
         CALL DESCINIT( DESCQ, N, N, DESCA(MB_), DESCA(NB_), 
     $                  DESCA(RSRC_), DESCA(CSRC_), ICTXT, 
     $                  MAX(1,AROWS), INFO )
C     
C     Compute the needed workspace for the QR-factorization using
C     workspace queries
C     
         CALL PDGEQRF( N, N, A, 1, 1, DESCA, DPDUM1, DPDUM2, -1, INFO )
         QRWORK = INT(DPDUM2(1))
         CALL PDORMQR( 'Left', 'Notranspose', N, N, N, DWORK(1), 1, 1,
     $                 DESCQ, DPDUM1, A, 1, 1, DESCA, DPDUM2, -1, INFO )
         QRWORK = MAX( QRWORK, INT(DPDUM2(1)))
         CALL PDORMQR( 'Right', 'Transpose', N, N, N, DWORK(1), 1, 1,
     $                 DESCQ, DPDUM1, A, 1, 1, DESCA, DPDUM2, -1, INFO )
         QRWORK = MAX( QRWORK, INT(DPDUM2(1)))
         IF( .NOT. LSAME(UPPER, 'U')) WORK = WORK + ACOLS + QRWORK 
      END IF
C     
      IF( INFO.EQ.0 ) THEN
         IF( WORK.GT.LDWORK.AND.LDWORK.NE.-1 ) INFO = -14
         IF( LDWORK.EQ.-1 ) THEN
            DWORK(1) = WORK
            RETURN
         END IF
      END IF
C     
C     Checking if we may continue or should interrupt
C     
      IF( INFO.NE.0 ) THEN
         CALL PXERBLA( ICTXT, 'P1SQMATGD', -INFO )
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
C     
C     Set some logical value to decide how to proceed
C     
      GVDIAG = LSAME( DIAG, 'D' )
      GVSDIAG = LSAME( SDIAG, 'D' )
      IF( .NOT. GVSDIAG ) THEN
         QUASITRI = NQTRBL.GT.0
      ELSE
         QUASITRI = .FALSE.
      END IF
      UPPERTRI = LSAME( UPPER, 'U' )
C     
C     Divide the workspace between the different integer pointers
C     for ScaLAPACK arrays. 
C     
      IF( .NOT. UPPERTRI ) THEN
         Q     = 1
         TAU   = Q + ASIZE
         SWORK = TAU + ACOLS
      ELSE
         SWORK = 1
      END IF
C     
C     Start the real matrix generating work
C
      CALL PDLASET( 'All', N, N, 0D0, 0D0, A, 1, 1, DESCA )
      IF( .NOT. UPPERTRI ) THEN
C     
C     First init a random matrix in Q using PDMATGEN2. PDMATGEN2
C     generates uniformly distributed numbers in the interval [0,1]
C     
         CALL PDMATGEN2( ICTXT, 'Random', 'NoDiagDominant', N, N, MB, 
     $                   NB, DWORK(Q), DESCQ(LLD_), 0, 0, ASEED+7, 0, 
     $                   AROWS, 0, ACOLS, MYROW, MYCOL, NPROW, NPCOL )
C     
C     Now QR-factorize the matrix Q
C     
         CALL PDGEQRF( N, N, DWORK(Q), 1, 1, DESCQ, DWORK(TAU), 
     $                 DWORK(SWORK), QRWORK, INFO )
      END IF
C     
C     Now init a random upper triangular matrix in A using PDMATGEN2:
C     now we have a preliminary version of ( D + M ).
C     
      CALL PDMATGEN2( ICTXT, 'U', 'N', N, N, MB, NB, A, LLDA, 0, 0, 
     $                ASEED, 0, AROWS, 0, ACOLS, MYROW, MYCOL, NPROW, 
     $                NPCOL )
C
C     If a user specified diagonal matrix D is present:
C     Set the new diagonal in the matrix D and weight it with DW.
C     At the same time, weight the matrix M with MW. If no user 
C     specified diagonal matrix D is present, just scale the upper 
C     triangular matrix ( D + M ) using the given values of DW, MW.
C
      DO 10 J = 1, ACOLS
         DO 20 I = 1, AROWS
            GJ = INDXL2G( J, NB, MYCOL, DESCA( CSRC_ ), NPCOL )
            GI = INDXL2G( I, MB, MYROW, DESCA( RSRC_ ), NPROW )
            IF( GI.EQ.GJ ) THEN
               IF( GVDIAG ) THEN
                  A( (J-1)*LLDA + I ) = DW * DDIAG( GI )
               ELSE
                  A( (J-1)*LLDA + I ) = DW * A( (J-1)*LLDA + I ) 
               END IF
            ELSEIF( GI.LT.GJ ) THEN
               A( (J-1)*LLDA + I ) = MW * A( (J-1)*LLDA + I ) 
            END IF
 20      CONTINUE
 10   CONTINUE
C     
C     Set the subdiagonal elements - form a quasitriangular matrix
C     out of the matrix ( D + M ). Imagine the entries in the 2x2 blocks
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
C     If no user specified subdiagonal is present and we still want
C     ( D + M ) to be quasitriangular we set the 2x2 blocks using
C     the corresponing MW-scaled (b)-elements in M. 
C
      IF ( QUASITRI ) THEN
C     
C     Loop thorugh the matrix and create 2x2 blocks on the diagonal
C     
         STEP = ICEIL( N, NQTRBL)
         DO 30 GI = 1, N-1, STEP
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
                  CALL DGESD2D( ICTXT, 1, 1, A((LJ2-1)*LLDA + LI2),
     $                          1, RSRC3, CSRC3 )
               ELSEIF( MYROW.EQ.RSRC3 .AND. MYCOL.EQ.CSRC3 ) THEN
                  CALL DGERV2D( ICTXT, 1, 1, ELEM, 1, RSRC2, CSRC2 )
                  A((LJ3-1)*LLDA + LI3) = MINONE * ELEM
               END IF
            ELSEIF( MYROW.EQ.RSRC2 .AND. MYCOL.EQ.RSRC2 ) THEN
               A((LJ3-1)*LLDA + LI3) = MINONE*A((LJ2-1)*LLDA + LI2)
            END IF  
 30      CONTINUE
C     
C     If the user has specified a subdiagonal, set the corresponding
C     entries in ( D + M ). At the same time, scale the new b-elements.
C     using DW.
C     
      ELSEIF( GVSDIAG ) THEN 
         DO 40 GI = 1, N-1, 1
            IF( SUBDIAG(GI).NE.ZERO ) THEN
               CALL INFOG2L( GI, GI+1, DESCA, NPROW, NPCOL, MYROW, 
     $                       MYCOL, LI2, LJ2, RSRC2, CSRC2 )
               CALL INFOG2L( GI+1, GI, DESCA, NPROW, NPCOL, MYROW, 
     $                       MYCOL, LI3, LJ3, RSRC3, CSRC3 )
               IF( MYROW.EQ.RSRC3 .AND. MYCOL.EQ.RSRC3 ) THEN
                  A((LJ3-1)*LLDA + LI3) = DW * SUBDIAG(GI)
               END IF
               IF( MYROW.EQ.RSRC2 .AND. MYCOL.EQ.RSRC2 ) THEN
                  A((LJ2-1)*LLDA + LI2) = MINONE * DW * SUBDIAG(GI)
               END IF
            END IF
 40      CONTINUE
      END IF
C
      IF( .NOT. UPPERTRI ) THEN
C     
C     Perform the similarity transformation A = Q*A*Q^T. Start by
C     multiplying with Q from the left, and then with Q^T from the
C     right
C     
         CALL PDORMQR( 'Left', 'NoTranspose', N, N, N, DWORK(Q), 1, 1,
     $                 DESCQ, DWORK(TAU), A, 1, 1, DESCA, DWORK(SWORK),
     $                 QRWORK, INFO )
         
         CALL PDORMQR( 'Right', 'Transpose', N, N, N, DWORK(Q), 1, 1,
     $                 DESCQ, DWORK(TAU), A, 1, 1, DESCA, DWORK(SWORK),
     $                 QRWORK, INFO )
      END IF
C     
C     ...and we are done!
C     
      END
C     
C     End of P1SQMATGD
C     
C *** Last line of P1SQMATGD ***
