***********************************************************************
*                                                                     *
*     pqrrmmm.f: Matrix Market matrix reader                          *
*                                                                     *
***********************************************************************
      SUBROUTINE PQRRMMM( FILENAME, M, N, A, DESCA, INFO )
*
*     To read in a matrix stored in Matrix Market format into
*     the distributed M-by-N matrix A.
*
      IMPLICIT NONE
      INTEGER M, N, INFO
      CHARACTER*(*) FILENAME
      INTEGER DESCA(*)
      DOUBLE PRECISION A(*)
*
      LOGICAL            DEBUG
      PARAMETER          ( DEBUG = .FALSE. )
      INTEGER            BLOCK_CYCLIC_2D, CSRC_, CTXT_, DLEN_, DTYPE_,
     $                   LLD_, MB_, M_, NB_, N_, RSRC_
      INTEGER            LEN, LENGTH
      DOUBLE PRECISION   ZERO
      PARAMETER          ( BLOCK_CYCLIC_2D = 1, DLEN_ = 9, DTYPE_ = 1,
     $                     CTXT_ = 2, M_ = 3, N_ = 4, MB_ = 5, NB_ = 6,
     $                     RSRC_ = 7, CSRC_ = 8, LLD_ = 9,
     $                     ZERO = 0.0D+00, LENGTH = 1000000 )
      INTEGER ICTXT, NIN, NNZ, NPROW, NPCOL, MYROW, MYCOL, II, JJ,
     $        LLDA, AROWS, ACOLS, M0, N0, I, J
      DOUBLE PRECISION ALPHA
      CHARACTER*(128) LINE
      INTEGER IBUF( 2*LENGTH )
      DOUBLE PRECISION DBUF( LENGTH )
      INTEGER NUMROC
      EXTERNAL NUMROC, DGEBS2D, DGEBR2D, IGEBS2D, IGEBR2D
*
      INFO = 0
      ICTXT = DESCA( CTXT_ )
      NIN = 10
      CALL BLACS_GRIDINFO( ICTXT, NPROW, NPCOL, MYROW, MYCOL )
*
      IF(DEBUG .AND. MYROW+MYCOL.EQ.0)
     $   WRITE(*,*), '% ', FILENAME
      IF( MYROW+MYCOL.EQ.0 ) THEN
         OPEN( UNIT=NIN, FILE=FILENAME, STATUS = 'OLD' )
 1       CONTINUE
*           Skip comments
            READ( NIN, '(A)' ) LINE
         IF( INDEX(LINE, '%') .NE.0 ) GO TO 1
         READ( LINE, FMT = * ) M0, N0, NNZ
         IF( M.NE.M0 .OR. N.NE.N0 ) THEN
            CLOSE( NIN )
            WRITE(*,*) '%%% ERROR: Not correct dimensions on input! %%%'
            INFO = 1
            STOP
         END IF
         IF( NPROW*NPCOL.GT.1 )
     $      CALL IGEBS2D( ICTXT, 'All', ' ', 1, 1, NNZ, 1 )
      ELSE
         CALL IGEBR2D( ICTXT, 'All', ' ', 1, 1, NNZ, 1, 0, 0 )
      END IF
*
      AROWS = NUMROC( M, DESCA(MB_), MYROW, DESCA(RSRC_), NPROW )
      ACOLS = NUMROC( N, DESCA(NB_), MYCOL, DESCA(CSRC_), NPCOL )
      LLDA = DESCA(LLD_)
      DO J = 1, ACOLS
         DO I = 1, AROWS
            A( (J-1)*LLDA + I ) = ZERO
         END DO
      END DO
*
      I = 0
 10   CONTINUE
         LEN = MIN( NNZ-I, LENGTH-1 )
         IF( MYROW+MYCOL.EQ.0 ) THEN
            DO J = 1, LEN
               READ( NIN, FMT = * ) II, JJ, ALPHA
               IBUF( J*2-1 ) = II
               IBUF( J*2 ) = JJ
               DBUF( J ) = ALPHA
            END DO
            IBUF( LENGTH ) = LEN
            IF( NPROW*NPCOL.GT.1 ) THEN
               CALL IGEBS2D( ICTXT, 'All', ' ', LENGTH*2, 1, IBUF, 1 )
               CALL DGEBS2D( ICTXT, 'All', ' ', LEN, 1, DBUF, 1 )
            END IF
         ELSE
            CALL IGEBR2D( ICTXT, 'All', ' ', LENGTH*2, 1, IBUF, 1, 0, 0)
            LEN = IBUF( LENGTH )
            CALL DGEBR2D( ICTXT, 'All', ' ', LEN, 1, DBUF, 1, 0, 0 )
         END IF
         DO J = 1, LEN
            CALL PDELSET( A, IBUF( J*2-1 ), IBUF( J*2 ), DESCA,
     $           DBUF( J ) )
         END DO
         I = I + LEN
      IF( I.LT. NNZ ) GO TO 10
*
      IF( MYROW+MYCOL.EQ.0 )
     $   CLOSE( NIN )
*
      END
