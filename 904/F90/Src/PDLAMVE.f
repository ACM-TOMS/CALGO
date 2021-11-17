      SUBROUTINE PDLAMVE( UPLO, M, N, A, IA, JA, DESCA, B, IB, JB,
     $                    DESCB, DWORK )

C  -- ScaLAPACK/PSLICOT-style routine --
C     Preliminary version.
C     Dept. Computing Science and HPC2N, Univ. of Umeå, Sweden
C     December 7, 2007.
C
      IMPLICIT NONE
C
C     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            IA, IB, JA, JB, M, N
C     ..
C     .. Array Arguments ..
      INTEGER            DESCA( * ), DESCB( * )
      DOUBLE PRECISION   A( * ), B( * ), DWORK( * )
C     ..
C
C  Purpose
C  =======
C
C  PDLAMVE copies all or part of a distributed matrix A to another
C  distributed matrix B. There is no alignment assumptions at all
C  except that A and B are of the same size. Get it?
C
C  Notes
C  =====
C
C  Each global data object is described by an associated description
C  vector.  This vector stores the information required to establish
C  the mapping between an object element and its corresponding process
C  and memory location.
C
C  Let A be a generic term for any 2D block cyclicly distributed array.
C  Such a global array has an associated description vector DESCA.
C  In the following comments, the character _ should be read as
C  "of the global array".
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
C  Let K be the number of rows or columns of a distributed matrix,
C  and assume that its process grid has dimension p x q.
C  LOCr( K ) denotes the number of elements of K that a process
C  would receive if K were distributed over the p processes of its
C  process column.
C  Similarly, LOCc( K ) denotes the number of elements of K that a
C  process would receive if K were distributed over the q processes of
C  its process row.
C  The values of LOCr() and LOCc() may be determined via a call to the
C  ScaLAPACK tool function, NUMROC:
C          LOCr( M ) = NUMROC( M, MB_A, MYROW, RSRC_A, NPROW ),
C          LOCc( N ) = NUMROC( N, NB_A, MYCOL, CSRC_A, NPCOL ).
C  An upper bound for these quantities may be computed by:
C          LOCr( M ) <= ceil( ceil(M/MB_A)/NPROW )*MB_A
C          LOCc( N ) <= ceil( ceil(N/NB_A)/NPCOL )*NB_A
C
C  Arguments
C  =========
C
C  UPLO    (global input) CHARACTER
C          Specifies the part of the distributed matrix sub( A ) to be
C          copied:
C          = 'U':   Upper triangular part is copied; the strictly
C                   lower triangular part of sub( A ) is not referenced;
C          = 'L':   Lower triangular part is copied; the strictly
C                   upper triangular part of sub( A ) is not referenced;
C          Otherwise:  All of the matrix sub( A ) is copied.
C
C  M       (global input) INTEGER
C          The number of rows to be operated on i.e the number of rows
C          of the distributed submatrix sub( A ). M >= 0.
C
C  N       (global input) INTEGER
C          The number of columns to be operated on i.e the number of
C          columns of the distributed submatrix sub( A ). N >= 0.
C
C  A       (local input) DOUBLE PRECISION pointer into the local memory
C          to an array of dimension (LLD_A, LOCc(JA+N-1) ). This array
C          contains the local pieces of the distributed matrix sub( A )
C          to be copied from.
C
C  IA      (global input) INTEGER
C          The row index in the global array A indicating the first
C          row of sub( A ).
C
C  JA      (global input) INTEGER
C          The column index in the global array A indicating the
C          first column of sub( A ).
C
C  DESCA   (global and local input) INTEGER array of dimension DLEN_.
C          The array descriptor for the distributed matrix A.
C
C  B       (local output) DOUBLE PRECISION pointer into the local memory
C          to an array of dimension (LLD_B, LOCc(JB+N-1) ). This array
C          contains on exit the local pieces of the distributed matrix
C          sub( B ).
C
C  IB      (global input) INTEGER
C          The row index in the global array B indicating the first
C          row of sub( B ).
C
C  JB      (global input) INTEGER
C          The column index in the global array B indicating the
C          first column of sub( B ).
C
C  DESCB   (global and local input) INTEGER array of dimension DLEN_.
C          The array descriptor for the distributed matrix B.
C
C  DWORK   (local workspace) DOUBLE PRECISION array of length ???
C
C  =====================================================================
C
C     .. Parameters ..
      LOGICAL            DEBUG, PRINT
      INTEGER            BLOCK_CYCLIC_2D, CSRC_, CTXT_, DLEN_, DTYPE_,
     $                   LLD_, MB_, M_, NB_, N_, RSRC_
      PARAMETER          ( BLOCK_CYCLIC_2D = 1, DLEN_ = 9, DTYPE_ = 1,
     $                     CTXT_ = 2, M_ = 3, N_ = 4, MB_ = 5, NB_ = 6,
     $                     RSRC_ = 7, CSRC_ = 8, LLD_ = 9,
     $                     DEBUG = .FALSE., PRINT = .FALSE. )
C     ..
C     .. Local Scalars ..
      LOGICAL            UPPER, LOWER, FULL
      INTEGER            ICTXT, NPROW, NPCOL, MYROW, MYCOL, MYPROC,
     $                   NPROCS, AROWS, ACOLS, K, SPROC, SRSRC, SCSRC,
     $                   RPROC, RRSRC, RCSRC, COUNT, J, I, IIA, JJA,
     $                   IIB, JJB, BRSRC, BCSRC, RAROWS, RACOLS,
     $                   INDEX, IDUM, NUMREC, NUMSND
C     ..
C     .. External Subroutines ..
      EXTERNAL           DLACPY, INFOG2L
C     ..
C     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            ICEIL, NUMROC, INDXL2G
      EXTERNAL           ICEIL, LSAME, NUMROC, INDXL2G
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC          MIN, MOD
C     ..
C     .. Executable Statements ..
C
      IF(DEBUG) WRITE(*,*) MYROW,MYCOL, 'PDLAMVE'
      IF(DEBUG) WRITE(*,*) MYROW,MYCOL, M,N,IA,JA,IB,JB
      IF(PRINT) 
     $     CALL PDLAPRNT( M, N, A, IA, JA, DESCA, 0, 0, 'A', 6, DWORK )
C
C     Find underlying mesh properties
C
      ICTXT = DESCA( CTXT_ )
      CALL BLACS_GRIDINFO( ICTXT, NPROW, NPCOL, MYROW, MYCOL )
C
C     Decode input parameters
C
      UPPER = LSAME( UPLO, 'U' )
      IF( .NOT. UPPER ) LOWER = LSAME( UPLO, 'L' )
      FULL = (.NOT. UPPER) .AND. (.NOT. LOWER)
C
C     Assign indiviual numbers based on column major ordering
C
      NPROCS = NPROW*NPCOL
      IF(DEBUG) WRITE(*,*) MYROW,MYCOL, 'NPROCS =',NPROCS
C
C     Do redistribution operation
C
      IF( NPROCS.EQ.1 ) THEN
         CALL DLACPY( UPLO, M, N, A((JA-1)*DESCA(LLD_)+IA),
     $                DESCA(LLD_), B((JB-1)*DESCB(LLD_)+IB),
     $                DESCB(LLD_) )
      ELSEIF( FULL ) THEN
         CALL PDGEMR2D( M, N, A, IA, JA, DESCA, B, IB, JB, DESCB, 
     $                  ICTXT )
      ELSE
         CALL PDGEMR2D( M, N, A, IA, JA, DESCA, DWORK, IB, JB, DESCB, 
     $                  ICTXT )
         CALL PDLACPY( UPLO, M, N, DWORK, IB, JB, DESCB, B, IB, JB,
     $                 DESCB )
      END IF
C
      IF( PRINT ) 
     $     CALL PDLAPRNT( M, N, B, IB, JB, DESCB, 0, 0, 'B', 6, DWORK )
C
      RETURN
C
C *** Last line of PDLAMVE ***
      END
