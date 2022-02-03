      SUBROUTINE PQRRMMM( FILENAME, M, N, A, DESCA, INFO )
C
C     To read in a matrix stored in Matrix Market format into
C     the distributed M-by-N matrix A. 
C
      IMPLICIT NONE
      INTEGER M, N, INFO
      CHARACTER*(*) FILENAME
      INTEGER DESCA(*)
      DOUBLE PRECISION A(*)
C
      INTEGER            BLOCK_CYCLIC_2D, CSRC_, CTXT_, DLEN_, DTYPE_,
     $                   LLD_, MB_, M_, NB_, N_, RSRC_
      DOUBLE PRECISION   ZERO
      PARAMETER          ( BLOCK_CYCLIC_2D = 1, DLEN_ = 9, DTYPE_ = 1,
     $                     CTXT_ = 2, M_ = 3, N_ = 4, MB_ = 5, NB_ = 6,
     $                     RSRC_ = 7, CSRC_ = 8, LLD_ = 9,
     $                     ZERO = 0.0D+00 )
      INTEGER ICTXT, NIN, NNZ, NPROW, NPCOL, MYROW, MYCOL, II, JJ,
     $        LLDA, AROWS, ACOLS, M0, N0, I, J
      DOUBLE PRECISION ALPHA
      INTEGER NUMROC
      EXTERNAL NUMROC
C
      INFO = 0
      ICTXT = DESCA( CTXT_ )
      CALL BLACS_GRIDINFO( ICTXT, NPROW, NPCOL, MYROW, MYCOL )
C     
      OPEN( NIN, FILE=FILENAME, STATUS = 'OLD' )
      READ( NIN, FMT = * )
      READ( NIN, FMT = * ) M0, N0, NNZ
      IF( M.NE.M0 .OR. N.NE.N0 ) THEN
         CLOSE( NIN )
         WRITE(*,*) '*** ERROR: Not correct dimensions on input! ***'
         INFO = 1
         STOP
      END IF
C
      AROWS = NUMROC( M, DESCA(MB_), MYROW, DESCA(RSRC_), NPROW )
      ACOLS = NUMROC( N, DESCA(NB_), MYCOL, DESCA(CSRC_), NPCOL )
      LLDA = DESCA(LLD_)
      DO I = 1, AROWS
         DO J = 1, ACOLS
            A( (J-1)*LLDA + I ) = ZERO
         END DO
      END DO
C     
      DO I = 1, NNZ
         READ( NIN, FMT = * ) II, JJ, ALPHA
         CALL PDELSET( A, II, JJ, DESCA, ALPHA )
      END DO 
C      
      CLOSE( NIN )
C
      END 
      
