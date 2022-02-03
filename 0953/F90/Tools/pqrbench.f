***********************************************************************
*                                                                     *
*     pdmatgen2.f: Benchmark matrix generator                         *
*                                                                     *
***********************************************************************
      SUBROUTINE PQRBENCH( WHICH, N, A, DESCA, WORK, LWORK,
     $                     IWORK, LIWORK, INFO )
      IMPLICIT NONE
      INTEGER            BLOCK_CYCLIC_2D, CSRC_, CTXT_, DLEN_, DTYPE_,
     $                   LLD_, MB_, M_, NB_, N_, RSRC_
      DOUBLE PRECISION   ZERO, ONE, TEN
      PARAMETER          ( BLOCK_CYCLIC_2D = 1, DLEN_ = 9, DTYPE_ = 1,
     $                     CTXT_ = 2, M_ = 3, N_ = 4, MB_ = 5, NB_ = 6,
     $                     RSRC_ = 7, CSRC_ = 8, LLD_ = 9,
     $                     ZERO = 0.0D+00, ONE = 1.0D+00,
     $                     TEN = 10.0D+00 )
      INTEGER WHICH, N, INFO, LWORK, LIWORK, NZR, CSTART, CEND, REND,
     $        RSTART, NX, NY, NZ
      INTEGER DESCA(*), IWORK(*)
      DOUBLE PRECISION A(*), WORK(*)
*
      INTEGER I, J
      DOUBLE PRECISION ALPHA
      INTEGER ISEED(4)
*
      IF( WHICH.EQ.0 ) THEN
         CALL PQRRMMM( 'rg24.mtx', N, N, A, DESCA, INFO )
      ELSEIF( WHICH.EQ.1 ) THEN
         DO I = 1, N
            DO J = 1, N
               IF( I.EQ.1 ) THEN
                  ALPHA = DBLE( N-J+1 )
               ELSEIF( I.EQ.J ) THEN
                  ALPHA = DBLE( I-1 )
               ELSEIF( I.EQ.J+1 ) THEN
                  ALPHA = TEN ** (-3)
               ELSE
                  ALPHA = ZERO
               END IF
               CALL PDELSET( A, I, J, DESCA, ALPHA )
            END DO
         END DO
      ELSEIF( WHICH.EQ.2 ) THEN
         CALL PQRRMMM( 'af23560.mtx', N, N, A, DESCA, INFO )
      ELSEIF( WHICH.EQ.3 ) THEN
         CALL PQRRMMM( 'cry10000.mtx', N, N, A, DESCA, INFO )
      ELSEIF( WHICH.EQ.4 ) THEN
         CALL PQRRMMM( 'olm5000.mtx', N, N, A, DESCA, INFO )
      ELSEIF( WHICH.EQ.5 ) THEN
         CALL PQRRMMM( 'dw8192.mtx', N, N, A, DESCA, INFO )
      ELSEIF( WHICH.EQ.6 ) THEN
         NZR = N / 100
         ISEED(1) = MIN(DESCA(MB_),4095)
         ISEED(2) = MIN(N/DESCA(MB_),4095)
         ISEED(3) = MOD(N,4095)
         ISEED(4) = 7
         CALL MATRAN( N, NZR, ISEED, WORK, IWORK, IWORK(N+2) )
         CALL PDLASET( 'All', N, N, ZERO, ZERO, A, 1, 1, DESCA )
         DO J = 1, N
            CSTART = IWORK(J)
            CEND = IWORK(J+1)-1
            DO I = CSTART, CEND
              CALL PDELSET( A, IWORK(N+2+I-1), J, DESCA, WORK(I) )
            END DO
         END DO
      ELSEIF( WHICH.EQ.7 ) THEN
         NX = INT( SQRT(REAL(N)) )
         NY = N / NX
         IF( N.NE.NX*NY ) THEN
            WRITE(*,*)
     $           '*** ERROR: Not correct dimensions on input! ***'
            INFO = 1
            STOP
         END IF
         NZ = NX*NY + (NX-1)*NY*2 + NX*(NY-1)*2
         CALL MATPDE( I, NX, NY, WORK, IWORK, IWORK(N+2),
     $                WORK(NZ+1) )
         CALL PDLASET( 'All', N, N, ZERO, ZERO, A, 1, 1, DESCA )
         DO I = 1, N
            RSTART = IWORK(I)
            REND = IWORK(I+1)-1
            DO J = RSTART, REND
              CALL PDELSET( A, I, IWORK(N+2+J-1), DESCA, WORK(J) )
            END DO
         END DO
      ELSEIF( WHICH.EQ.8 ) THEN
         DO I = 1, N
            DO J = 1, N
               IF( I.EQ.J ) THEN
                  ALPHA = ONE
               ELSEIF( I.EQ.J+1 ) THEN
                  ALPHA = -ONE
               ELSEIF( J.GE.I+1 .AND. J.LE.I+3 ) THEN
                  ALPHA = ONE
               ELSE
                  ALPHA = ZERO
               END IF
               CALL PDELSET( A, I, J, DESCA, ALPHA )
            END DO
         END DO
      ELSE
         INFO = 1
      END IF
*
      RETURN
*
      END
