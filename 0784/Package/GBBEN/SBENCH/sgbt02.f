      SUBROUTINE SGBT02( TABSUB, SIDE, NSIDE, NUPLO, TRNS, NTRNS, NDIAG,
     $                         DIM1, DIM2, NDIM, LDA, NLDA, ALPHA, BETA,
     $                            A, B, C, LD, NMAX, NERR, MXSUB, MXOPT,
     $                                         MXDIM, MXLDA, RUNS, RES )
*     .. Scalar Arguments ..
      INTEGER            LD, NMAX, NERR,
     $                   NSIDE, NUPLO, NTRNS, NDIAG, NDIM, NLDA,
     $                   MXSUB, MXOPT, MXDIM, MXLDA, RUNS
      REAL               ALPHA, BETA
*     .. Array Arguments ..
      LOGICAL            TABSUB( MXSUB )
      CHARACTER          SIDE( MXOPT ), TRNS( MXOPT )
      INTEGER            DIM1( MXDIM ), DIM2( MXDIM ), LDA( MXLDA )
      REAL               A( LD, NMAX ), B( LD, NMAX ), C( LD, NMAX ),
     $                   RES( MXSUB, MXOPT, MXOPT, MXOPT, MXOPT, MXDIM,
     $                   MXLDA )
*
*
*  Determine problem configurations for SGEMM that, for timing purposes,
*  "correspond" to problem configurations for the remaining Level 3 BLAS
*  routines. Time SGEMM for problems that correspond to the Level 3 BLAS
*  problems timed in SGBT01. Return the performance of SGEMM in Mflops.
*
*  SGBT02 calls a REAL function SECOND with no arguments,
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
      REAL               TIME, SPEED, TM0, TM1, TM2, TM3, TM4, TM5, TM6,
     $                   TM7, TM8, TM9
*     .. Intrinsic Functions ..
      INTRINSIC          REAL, MIN
*     .. External Functions ..
      LOGICAL            LSAME
      REAL               SECOND
      EXTERNAL           LSAME, SECOND
*     .. External Subroutines ..
      EXTERNAL           SGEMM
*     .. Parameters ..
      REAL               ZERO, ONE, SCALE
      PARAMETER        ( ZERO = 0.0E+0, ONE = 1.0E+0, SCALE = 1.0E+6 )
*     ..
*     .. Executable Statements ..
      TM0 = SECOND( )
      TM0 = SECOND( )
      TM0 = SECOND( )
      TM1 = SECOND( )
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
                  C( I, J ) = ONE + 0.01E+0*REAL( I+( J-1 )*NMAX )/
     $                                                 REAL( LD*NMAX+1 )
   10          CONTINUE
   20       CONTINUE
*
*           Time SGEMM for a problem that corresponds to the following
*           problem for SSYMM:
*           SSYMM( SIDE( OP1 ), UPLO( OP2 ), DIM1( D ), DIM2( D ),
*                   ALPHA, A, LDA( L ), B, LDA( L ), BETA, C, LDA( L ) )
*
            IF( LSAME( SIDE( OP1 ), 'L' ) )THEN
*
*              Use K = M.
*
               TM2 = SECOND( )
               CALL SGEMM( 'N', 'N', DIM1( D ), DIM2( D ), DIM1( D ),
     $              ALPHA, A, LDA( L ), B, LDA( L ), BETA, C, LDA( L ) )
               TM3 = SECOND( )
            ELSE IF( LSAME( SIDE( OP1 ), 'R' ) )THEN
*
*              Use K = N.
*
               TM2 = SECOND( )
               CALL SGEMM( 'N', 'N', DIM1( D ), DIM2( D ), DIM2( D ),
     $              ALPHA, B, LDA( L ), A, LDA( L ), BETA, C, LDA( L ) )
               TM3 = SECOND( )
            ELSE
               WRITE( NERR, FMT = * )
     $               'Error: Unknown value for SIDE: ', SIDE( OP1 ), '.'
               STOP
            END IF
         END IF
         IF( TABSUB( 2 ).AND.OP1.EQ.1 )THEN
            DO 40, J = 1, NMAX
               DO 30, I = 1, LD
                  C( I, J ) = ONE + 0.01E+0*REAL( I+( J-1 )*NMAX )/
     $                                                 REAL( LD*NMAX+1 )
   30          CONTINUE
   40       CONTINUE
*
*           Time SGEMM for a problem that corresponds to the following
*           problem for SSYRK:
*           SSYRK( UPLO( OP2 ), TRNS( OP3 ), DIM1( D ), DIM2( D ),
*                                ALPHA, A, LDA( L ), BETA, C, LDA( L ) )
*           Use M = N and B = A in the call to SGEMM.
*
            IF( LSAME( TRNS( OP3 ), 'N' ) )THEN
               TM4 = SECOND( )
               CALL SGEMM( 'N', 'T', DIM1( D ), DIM1( D ), DIM2( D ),
     $              ALPHA, A, LDA( L ), A, LDA( L ), BETA, C, LDA( L ) )
               TM5 = SECOND( )
            ELSE IF( LSAME( TRNS( OP3 ), 'T' ) )THEN
               TM4 = SECOND( )
               CALL SGEMM( 'T', 'N', DIM1( D ), DIM1( D ), DIM2( D ),
     $              ALPHA, A, LDA( L ), A, LDA( L ), BETA, C, LDA( L ) )
               TM5 = SECOND( )
            ELSE
               WRITE( NERR, FMT = * )
     $              'Error: Unknown value for TRANS: ', TRNS( OP3 ), '.'
               STOP
            END IF
         END IF
         IF( TABSUB( 3 ).AND.OP1.EQ.1 )THEN
            DO 60, J = 1, NMAX
               DO 50, I = 1, LD
                  C( I, J ) = ONE + 0.01E+0*REAL( I+( J-1 )*NMAX )/
     $                                                 REAL( LD*NMAX+1 )
   50          CONTINUE
   60       CONTINUE
*
*           Time SGEMM for a problem that corresponds to the following
*           problem for SSYR2K:
*           SSYR2K( UPLO( OP2 ), TRNS( OP3 ), DIM1( D ), DIM2( D ),
*                   ALPHA, A, LDA( L ), B, LDA( L ), BETA, C, LDA( L ) )
*
            IF( LSAME( TRNS( OP3 ), 'N' ) )THEN
               TM6 = SECOND( )
               CALL SGEMM( 'N', 'T', DIM1( D ), DIM1( D ), DIM2( D ),
     $              ALPHA, A, LDA( L ), B, LDA( L ), BETA, C, LDA( L ) )
               TM7 = SECOND( )
            ELSE IF( LSAME( TRNS( OP3 ), 'T' ) )THEN
               TM6 = SECOND( )
               CALL SGEMM( 'T', 'N', DIM1( D ), DIM1( D ), DIM2( D ),
     $              ALPHA, A, LDA( L ), B, LDA( L ), BETA, C, LDA( L ) )
               TM7 = SECOND( )
            ELSE
               WRITE( NERR, FMT = * )
     $              'Error: Unknown value for TRANS: ', TRNS( OP3 ), '.'
               STOP
            END IF
         END IF
         IF( TABSUB( 4 ).OR.TABSUB( 5 ) )THEN
            DO 80, J = 1, NMAX
               DO 70, I = 1, LD
                  C( I, J ) = ONE + 0.01E+0*REAL( I+( J-1 )*NMAX )/
     $                                                 REAL( LD*NMAX+1 )
   70          CONTINUE
   80       CONTINUE
*
*           Time SGEMM for a problem that corresponds to the following
*           problems for STRMM and STRSM:
*           STRMM( SIDE( OP1 ), UPLO( OP2 ), TRNS( OP3 ),
*                              DIAG( OP4 ), DIM1( D ), DIM2( D ), ALPHA,
*                                             A, LDA( L ), C, LDA( L ) )
*           STRSM( SIDE( OP1 ), UPLO( OP2 ), TRNS( OP3 ),
*                              DIAG( OP4 ), DIM1( D ), DIM2( D ), ALPHA,
*                                             A, LDA( L ), C, LDA( L ) )
*
            IF( LSAME( SIDE( OP1 ), 'L' ) )THEN
*
*              C := alpha*A*C + C or C := alpha*A'*C + C. Use K = M.
*
               TM8 = SECOND( )
               CALL SGEMM( TRNS( OP3 ), 'N', DIM1( D ), DIM2( D ),
     $                       DIM1( D ), ALPHA, A, LDA( L ), B, LDA( L ),
     $                                                ONE, C, LDA( L ) )
               TM9 = SECOND( )
            ELSE IF( LSAME( SIDE( OP1 ), 'R' ) )THEN
*
*              C := alpha*C*A + C or C := alpha*C*A' + C. Use K = N.
*
               TM8 = SECOND( )
               CALL SGEMM( 'N', TRNS( OP3 ), DIM1( D ), DIM2( D ),
     $                       DIM2( D ), ALPHA, B, LDA( L ), A, LDA( L ),
     $                                                ONE, C, LDA( L ) )
               TM9 = SECOND( )
            ELSE
               WRITE( NERR, FMT = * )
     $               'Error: Unknown value for SIDE: ', SIDE( OP1 ), '.'
               STOP
            END IF
         END IF
*
*        Compute the performance of SGEMM in Mflops for problem
*        configurations that corresponds to SSYMM.
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
               SPEED = REAL( NOPS )/( TIME*SCALE )
            END IF
            IF( RES( 1, OP1, 1, OP3, 1, D, L ).LT.SPEED )THEN
               DO 90, OP2 = 1, NUPLO
                  RES( 1, OP1, OP2, OP3, 1, D, L ) = SPEED
   90          CONTINUE
            END IF
         END IF
*
*        Compute the performance of SGEMM in Mflops for problem
*        configurations that corresponds to SSYRK.
*
         IF( TABSUB( 2 ).AND.OP1.EQ.1 )THEN
            TIME = ( TM5 - TM4 ) - ( TM1 - TM0 )
            N = DIM1( D )
            K = DIM2( D )
            NOPS = ( 2*K + 1 )*N*N + MIN( N*K, N*N )
            IF( TIME.LE.ZERO )THEN
               SPEED = ZERO
            ELSE
               SPEED = REAL( NOPS )/( TIME*SCALE )
            END IF
            IF( RES( 2, OP1, 1, OP3, 1, D, L ).LT.SPEED )THEN
               DO 100, OP2 = 1, NUPLO
                  RES( 2, OP1, OP2, OP3, 1, D, L ) = SPEED
  100          CONTINUE
            END IF
         END IF
*
*        Compute the performance of SGEMM in Mflops for problem
*        configurations that corresponds to SSYR2K.
*
         IF( TABSUB( 3 ).AND.OP1.EQ.1 )THEN
            TIME = ( TM7 - TM6 ) - ( TM1 - TM0 )
            N = DIM1( D )
            K = DIM2( D )
            NOPS = ( 2*K + 1 )*N*N + MIN( N*K, N*N )
            IF( TIME.LE.ZERO )THEN
               SPEED = ZERO
            ELSE
               SPEED = REAL( NOPS )/( TIME*SCALE )
            END IF
            IF( RES( 3, OP1, 1, OP3, 1, D, L ).LT.SPEED )THEN
               DO 110, OP2 = 1, NUPLO
                  RES( 3, OP1, OP2, OP3, 1, D, L ) = SPEED
  110          CONTINUE
            END IF
         END IF
*
*        Compute the performance of SGEMM in Mflops for problem
*        configurations that corresponds to STRMM and STRSM.
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
               SPEED = REAL( NOPS )/( TIME*SCALE )
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
*     End of SGBT02.
*
      END
