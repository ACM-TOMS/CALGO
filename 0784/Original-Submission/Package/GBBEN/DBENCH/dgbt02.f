      SUBROUTINE DGBT02( TABSUB, SIDE, NSIDE, NUPLO, TRNS, NTRNS, NDIAG,
     $                         DIM1, DIM2, NDIM, LDA, NLDA, ALPHA, BETA,
     $                            A, B, C, LD, NMAX, NERR, MXSUB, MXOPT,
     $                                         MXDIM, MXLDA, RUNS, RES )
*     .. Scalar Arguments ..
      INTEGER            LD, NMAX, NERR,
     $                   NSIDE, NUPLO, NTRNS, NDIAG, NDIM, NLDA,
     $                   MXSUB, MXOPT, MXDIM, MXLDA, RUNS
      DOUBLE PRECISION   ALPHA, BETA
*     .. Array Arguments ..
      LOGICAL            TABSUB( MXSUB )
      CHARACTER          SIDE( MXOPT ), TRNS( MXOPT )
      INTEGER            DIM1( MXDIM ), DIM2( MXDIM ), LDA( MXLDA )
      DOUBLE PRECISION   A( LD, NMAX ), B( LD, NMAX ), C( LD, NMAX ),
     $                   RES( MXSUB, MXOPT, MXOPT, MXOPT, MXOPT, MXDIM,
     $                   MXLDA )
*
*
*  Determine problem configurations for DGEMM that, for timing purposes,
*  "correspond" to problem configurations for the remaining Level 3 BLAS
*  routines. Time DGEMM for problems that correspond to the Level 3 BLAS
*  problems timed in DGBT01. Return the performance of DGEMM in Mflops.
*
*  DGBT02 calls a DOUBLE PRECISION function DSECND with no arguments,
*  which is assumed to return the user time for a process in seconds
*  from some fixed starting-time.
*
*
*  -- Written in January-1994.
*     GEMM-Based Level 3 BLAS Benchmark.
*     Per Ling, Institute of Information Processing,
*     University of Umea, Sweden.
*
*  -- Modified in February-1996.
*     Per Ling, Department of Computing Science,
*     Umea University, Sweden.
*
*
*     .. Local Scalars ..
      INTEGER            I, J, M, N, K, NOPS,
     $                   D, L, R, OP1, OP2, OP3, OP4
      DOUBLE PRECISION   TIME, SPEED, TM0, TM1, TM2, TM3, TM4, TM5, TM6,
     $                   TM7, TM8, TM9
*     .. Intrinsic Functions ..
      INTRINSIC          DBLE, MIN
*     .. External Functions ..
      LOGICAL            LSAME
      DOUBLE PRECISION   DSECND
      EXTERNAL           LSAME, DSECND
*     .. External Subroutines ..
      EXTERNAL           DGEMM
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE, SCALE
      PARAMETER        ( ZERO = 0.0D+0, ONE = 1.0D+0, SCALE = 1.0D+6 )
*     ..
*     .. Executable Statements ..
      TM0 = DSECND( )
      TM0 = DSECND( )
      TM0 = DSECND( )
      TM1 = DSECND( )
*
*     ------ Stop indentation ------
*
      DO 180, L = 1, NLDA
      DO 170, OP1 = 1, NSIDE
      DO 160, OP3 = 1, NTRNS
      DO 150, D = 1, NDIM
*
*     ------ Continue indentation ------
*
      RES( 1, OP1, 1, OP3, 1, D, L ) = ZERO
      RES( 2, OP1, 1, OP3, 1, D, L ) = ZERO
      RES( 3, OP1, 1, OP3, 1, D, L ) = ZERO
      RES( 4, OP1, 1, OP3, 1, D, L ) = ZERO
      RES( 5, OP1, 1, OP3, 1, D, L ) = ZERO
      DO 140, R = 1, RUNS
*
*        Time the user-supplied library. Re-initialize C between the
*        timings to avoid overflow and to flush the cache.
*
         IF( TABSUB( 1 ).AND.OP3.EQ.1 )THEN
            DO 20, J = 1, NMAX
               DO 10, I = 1, LD
                  C( I, J ) = ONE + 0.01D+0*DBLE( I+( J-1 )*NMAX )/
     $                                                 DBLE( LD*NMAX+1 )
   10          CONTINUE
   20       CONTINUE
*
*           Time DGEMM for a problem that corresponds to the following
*           problem for DSYMM:
*           DSYMM( SIDE( OP1 ), UPLO( OP2 ), DIM1( D ), DIM2( D ),
*                   ALPHA, A, LDA( L ), B, LDA( L ), BETA, C, LDA( L ) )
*
            IF( LSAME( SIDE( OP1 ), 'L' ) )THEN
*
*              Use K = M.
*
               TM2 = DSECND( )
               CALL DGEMM( 'N', 'N', DIM1( D ), DIM2( D ), DIM1( D ),
     $              ALPHA, A, LDA( L ), B, LDA( L ), BETA, C, LDA( L ) )
               TM3 = DSECND( )
            ELSE IF( LSAME( SIDE( OP1 ), 'R' ) )THEN
*
*              Use K = N.
*
               TM2 = DSECND( )
               CALL DGEMM( 'N', 'N', DIM1( D ), DIM2( D ), DIM2( D ),
     $              ALPHA, B, LDA( L ), A, LDA( L ), BETA, C, LDA( L ) )
               TM3 = DSECND( )
            ELSE
               WRITE( NERR, FMT = * )
     $               'Error: Unknown value for SIDE: ', SIDE( OP1 ), '.'
               STOP
            END IF
         END IF
         IF( TABSUB( 2 ).AND.OP1.EQ.1 )THEN
            DO 40, J = 1, NMAX
               DO 30, I = 1, LD
                  C( I, J ) = ONE + 0.01D+0*DBLE( I+( J-1 )*NMAX )/
     $                                                 DBLE( LD*NMAX+1 )
   30          CONTINUE
   40       CONTINUE
*
*           Time DGEMM for a problem that corresponds to the following
*           problem for DSYRK:
*           DSYRK( UPLO( OP2 ), TRNS( OP3 ), DIM1( D ), DIM2( D ),
*                                ALPHA, A, LDA( L ), BETA, C, LDA( L ) )
*           Use M = N and B = A in the call to DGEMM.
*
            IF( LSAME( TRNS( OP3 ), 'N' ) )THEN
               TM4 = DSECND( )
               CALL DGEMM( 'N', 'T', DIM1( D ), DIM1( D ), DIM2( D ),
     $              ALPHA, A, LDA( L ), A, LDA( L ), BETA, C, LDA( L ) )
               TM5 = DSECND( )
            ELSE IF( LSAME( TRNS( OP3 ), 'T' ) )THEN
               TM4 = DSECND( )
               CALL DGEMM( 'T', 'N', DIM1( D ), DIM1( D ), DIM2( D ),
     $              ALPHA, A, LDA( L ), A, LDA( L ), BETA, C, LDA( L ) )
               TM5 = DSECND( )
            ELSE
               WRITE( NERR, FMT = * )
     $              'Error: Unknown value for TRANS: ', TRNS( OP3 ), '.'
               STOP
            END IF
         END IF
         IF( TABSUB( 3 ).AND.OP1.EQ.1 )THEN
            DO 60, J = 1, NMAX
               DO 50, I = 1, LD
                  C( I, J ) = ONE + 0.01D+0*DBLE( I+( J-1 )*NMAX )/
     $                                                 DBLE( LD*NMAX+1 )
   50          CONTINUE
   60       CONTINUE
*
*           Time DGEMM for a problem that corresponds to the following
*           problem for DSYR2K:
*           DSYR2K( UPLO( OP2 ), TRNS( OP3 ), DIM1( D ), DIM2( D ),
*                   ALPHA, A, LDA( L ), B, LDA( L ), BETA, C, LDA( L ) )
*
            IF( LSAME( TRNS( OP3 ), 'N' ) )THEN
               TM6 = DSECND( )
               CALL DGEMM( 'N', 'T', DIM1( D ), DIM1( D ), DIM2( D ),
     $              ALPHA, A, LDA( L ), B, LDA( L ), BETA, C, LDA( L ) )
               TM7 = DSECND( )
            ELSE IF( LSAME( TRNS( OP3 ), 'T' ) )THEN
               TM6 = DSECND( )
               CALL DGEMM( 'T', 'N', DIM1( D ), DIM1( D ), DIM2( D ),
     $              ALPHA, A, LDA( L ), B, LDA( L ), BETA, C, LDA( L ) )
               TM7 = DSECND( )
            ELSE
               WRITE( NERR, FMT = * )
     $              'Error: Unknown value for TRANS: ', TRNS( OP3 ), '.'
               STOP
            END IF
         END IF
         IF( TABSUB( 4 ).OR.TABSUB( 5 ) )THEN
            DO 80, J = 1, NMAX
               DO 70, I = 1, LD
                  C( I, J ) = ONE + 0.01D+0*DBLE( I+( J-1 )*NMAX )/
     $                                                 DBLE( LD*NMAX+1 )
   70          CONTINUE
   80       CONTINUE
*
*           Time DGEMM for a problem that corresponds to the following
*           problems for DTRMM and DTRSM:
*           DTRMM( SIDE( OP1 ), UPLO( OP2 ), TRNS( OP3 ),
*                              DIAG( OP4 ), DIM1( D ), DIM2( D ), ALPHA,
*                                             A, LDA( L ), C, LDA( L ) )
*           DTRSM( SIDE( OP1 ), UPLO( OP2 ), TRNS( OP3 ),
*                              DIAG( OP4 ), DIM1( D ), DIM2( D ), ALPHA,
*                                             A, LDA( L ), C, LDA( L ) )
*
            IF( LSAME( SIDE( OP1 ), 'L' ) )THEN
*
*              C := alpha*A*C + C or C := alpha*A'*C + C. Use K = M.
*
               TM8 = DSECND( )
               CALL DGEMM( TRNS( OP3 ), 'N', DIM1( D ), DIM2( D ),
     $                       DIM1( D ), ALPHA, A, LDA( L ), B, LDA( L ),
     $                                                ONE, C, LDA( L ) )
               TM9 = DSECND( )
            ELSE IF( LSAME( SIDE( OP1 ), 'R' ) )THEN
*
*              C := alpha*C*A + C or C := alpha*C*A' + C. Use K = N.
*
               TM8 = DSECND( )
               CALL DGEMM( 'N', TRNS( OP3 ), DIM1( D ), DIM2( D ),
     $                       DIM2( D ), ALPHA, B, LDA( L ), A, LDA( L ),
     $                                                ONE, C, LDA( L ) )
               TM9 = DSECND( )
            ELSE
               WRITE( NERR, FMT = * )
     $               'Error: Unknown value for SIDE: ', SIDE( OP1 ), '.'
               STOP
            END IF
         END IF
*
*        Compute the performance of DGEMM in Mflops for problem
*        configurations that corresponds to DSYMM.
*
         IF( TABSUB( 1 ).AND.OP3.EQ.1 )THEN
            TIME = ( TM3 - TM2 ) - ( TM1 - TM0 )
            M = DIM1( D )
            N = DIM2( D )
            IF( LSAME( SIDE( OP1 ), 'L' ) )THEN
               NOPS = ( 2*M + 1 )*M*N + MIN( M*N, M*M )
            ELSE IF( LSAME( SIDE( OP1 ), 'R' ) )THEN
               NOPS = ( 2*N + 1 )*M*N + MIN( M*N, N*N )
            END IF
            IF( TIME.LE.ZERO )THEN
               SPEED = ZERO
            ELSE
               SPEED = DBLE( NOPS )/( TIME*SCALE )
            END IF
            IF( RES( 1, OP1, 1, OP3, 1, D, L ).LT.SPEED )THEN
               DO 90, OP2 = 1, NUPLO
                  RES( 1, OP1, OP2, OP3, 1, D, L ) = SPEED
   90          CONTINUE
            END IF
         END IF
*
*        Compute the performance of DGEMM in Mflops for problem
*        configurations that corresponds to DSYRK.
*
         IF( TABSUB( 2 ).AND.OP1.EQ.1 )THEN
            TIME = ( TM5 - TM4 ) - ( TM1 - TM0 )
            N = DIM1( D )
            K = DIM2( D )
            NOPS = ( 2*K + 1 )*N*N + MIN( N*K, N*N )
            IF( TIME.LE.ZERO )THEN
               SPEED = ZERO
            ELSE
               SPEED = DBLE( NOPS )/( TIME*SCALE )
            END IF
            IF( RES( 2, OP1, 1, OP3, 1, D, L ).LT.SPEED )THEN
               DO 100, OP2 = 1, NUPLO
                  RES( 2, OP1, OP2, OP3, 1, D, L ) = SPEED
  100          CONTINUE
            END IF
         END IF
*
*        Compute the performance of DGEMM in Mflops for problem
*        configurations that corresponds to DSYR2K.
*
         IF( TABSUB( 3 ).AND.OP1.EQ.1 )THEN
            TIME = ( TM7 - TM6 ) - ( TM1 - TM0 )
            N = DIM1( D )
            K = DIM2( D )
            NOPS = ( 2*K + 1 )*N*N + MIN( N*K, N*N )
            IF( TIME.LE.ZERO )THEN
               SPEED = ZERO
            ELSE
               SPEED = DBLE( NOPS )/( TIME*SCALE )
            END IF
            IF( RES( 3, OP1, 1, OP3, 1, D, L ).LT.SPEED )THEN
               DO 110, OP2 = 1, NUPLO
                  RES( 3, OP1, OP2, OP3, 1, D, L ) = SPEED
  110          CONTINUE
            END IF
         END IF
*
*        Compute the performance of DGEMM in Mflops for problem
*        configurations that corresponds to DTRMM and DTRSM.
*
         IF( TABSUB( 4 ).OR.TABSUB( 5 ) )THEN
            TIME = ( TM9 - TM8 ) - ( TM1 - TM0 )
            M = DIM1( D )
            N = DIM2( D )
            IF( LSAME( SIDE( OP1 ), 'L' ) )THEN
               NOPS = ( 2*M - 1 )*M*N + MIN( M*N, M*M )
            ELSE IF( LSAME( SIDE( OP1 ), 'R' ) )THEN
               NOPS = ( 2*N - 1 )*M*N + MIN( M*N, N*N )
            END IF
            IF( TIME.LE.ZERO )THEN
               SPEED = ZERO
            ELSE
               SPEED = DBLE( NOPS )/( TIME*SCALE )
            END IF
            IF( RES( 4, OP1, 1, OP3, 1, D, L ).LT.SPEED )THEN
               DO 130, OP2 = 1, NUPLO
                  DO 120, OP4 = 1, NDIAG
                     RES( 4, OP1, OP2, OP3, OP4, D, L ) = SPEED
                     RES( 5, OP1, OP2, OP3, OP4, D, L ) = SPEED
  120             CONTINUE
  130          CONTINUE
            END IF
         END IF
  140 CONTINUE
*
*     ------ Stop indentation ------
*
  150 CONTINUE
  160 CONTINUE
  170 CONTINUE
  180 CONTINUE
*
*     ------ Continue indentation ------
*
      RETURN
*
*     End of DGBT02.
*
      END
