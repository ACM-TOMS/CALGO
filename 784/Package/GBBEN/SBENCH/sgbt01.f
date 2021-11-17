      SUBROUTINE SGBT01( SB3LIB, TABSUB, SIDE, NSIDE, UPLO, NUPLO, TRNS,
     $                  NTRNS, DIAG, NDIAG, DIM1, DIM2, NDIM, LDA, NLDA,
     $                      ALPHA, BETA, A, B, C, LD, NMAX, NERR, MXSUB,
     $                                  MXOPT, MXDIM, MXLDA, RUNS, RES )
*     .. Scalar Arguments ..
      CHARACTER          SB3LIB
      INTEGER            LD, NMAX, NERR,
     $                   NSIDE, NUPLO, NTRNS, NDIAG, NDIM, NLDA,
     $                   MXSUB, MXOPT, MXDIM, MXLDA, RUNS
      REAL               ALPHA, BETA
*     .. Array Arguments ..
      LOGICAL            TABSUB( MXSUB )
      CHARACTER          SIDE( MXOPT ), UPLO( MXOPT ), TRNS( MXOPT ),
     $                   DIAG( MXOPT )
      INTEGER            DIM1( MXDIM ), DIM2( MXDIM ), LDA( MXLDA )
      REAL               A( LD, NMAX ), B( LD, NMAX ), C( LD, NMAX ),
     $                   RES( MXSUB, MXOPT, MXOPT, MXOPT, MXOPT, MXDIM,
     $                   MXLDA )
*
*
*  Time all routines except SGEMM in the Level 3 BLAS library specified
*  by the input parameters. The library is either a user-supplied
*  Level 3 BLAS library or the GEMM-Based Level 3 BLAS library included
*  in the benchmark (SGB02, SGB04, SGB06, SGB08, and SGB09). Return the
*  performance in Mflops for each problem configuration.
*
*  SGBT01 calls a REAL function SECOND with no arguments,
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
     $                   TM7, TM8, TM9, TM10, TM11
*     .. Intrinsic Functions ..
      INTRINSIC          REAL, MIN
*     .. External Functions ..
      LOGICAL            LSAME
      REAL               SECOND
      EXTERNAL           LSAME, SECOND
*     .. External Subroutines ..
      EXTERNAL           SSYMM, SSYRK, SSYR2K, STRMM, STRSM,
     $                   SGB02, SGB04, SGB06, SGB08, SGB09
*     .. Parameters ..
      REAL               ZERO, ONE, SCALE
*     .. Parameter Values ..
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
      DO 270, L = 1, NLDA
      DO 260, OP1 = 1, NSIDE
      DO 250, OP2 = 1, NUPLO
      DO 240, OP3 = 1, NTRNS
      DO 230, OP4 = 1, NDIAG
      DO 220, D = 1, NDIM
*
*     ------ Continue indentation ------
*
      RES( 1, OP1, OP2, OP3, OP4, D, L ) = ZERO
      RES( 2, OP1, OP2, OP3, OP4, D, L ) = ZERO
      RES( 3, OP1, OP2, OP3, OP4, D, L ) = ZERO
      RES( 4, OP1, OP2, OP3, OP4, D, L ) = ZERO
      RES( 5, OP1, OP2, OP3, OP4, D, L ) = ZERO
      DO 210, R = 1, RUNS
         IF( LSAME( SB3LIB, 'U' ) )THEN
*
*           Time the user-supplied library. Re-initialize C between the
*           timings to avoid overflow and to flush the cache.
*
            IF( TABSUB( 1 ).AND.OP3.EQ.1.AND.OP4.EQ.1 )THEN
               DO 20, J = 1, NMAX
                  DO 10, I = 1, LD
                     C( I, J ) = ONE + 0.01E+0*REAL( I+( J-1 )*NMAX )/
     $                                                 REAL( LD*NMAX+1 )
   10             CONTINUE
   20          CONTINUE
               TM2 = SECOND( )
               CALL SSYMM( SIDE( OP1 ), UPLO( OP2 ), DIM1( D ),
     $                       DIM2( D ), ALPHA, A, LDA( L ), B, LDA( L ),
     $                                               BETA, C, LDA( L ) )
               TM3 = SECOND( )
            END IF
            IF( TABSUB( 2 ).AND.OP1.EQ.1.AND.OP4.EQ.1 )THEN
               DO 40, J = 1, NMAX
                  DO 30, I = 1, LD
                     C( I, J ) = ONE + 0.01E+0*REAL( I+( J-1 )*NMAX )/
     $                                                 REAL( LD*NMAX+1 )
   30             CONTINUE
   40          CONTINUE
               TM4 = SECOND( )
               CALL SSYRK( UPLO( OP2 ), TRNS( OP3 ), DIM1( D ),
     $                                    DIM2( D ), ALPHA, A, LDA( L ),
     $                                               BETA, C, LDA( L ) )
               TM5 = SECOND( )
            END IF
            IF( TABSUB( 3 ).AND.OP1.EQ.1.AND.OP4.EQ.1 )THEN
               DO 60, J = 1, NMAX
                  DO 50, I = 1, LD
                     C( I, J ) = ONE + 0.01E+0*REAL( I+( J-1 )*NMAX )/
     $                                                 REAL( LD*NMAX+1 )
   50             CONTINUE
   60          CONTINUE
               TM6 = SECOND( )
               CALL SSYR2K( UPLO( OP2 ), TRNS( OP3 ), DIM1( D ),
     $                       DIM2( D ), ALPHA, A, LDA( L ), B, LDA( L ),
     $                                               BETA, C, LDA( L ) )
               TM7 = SECOND( )
            END IF
            IF( TABSUB( 4 ) )THEN
               DO 80, J = 1, NMAX
                  DO 70, I = 1, LD
                     C( I, J ) = ONE + 0.01E+0*REAL( I+( J-1 )*NMAX )/
     $                                                 REAL( LD*NMAX+1 )
   70             CONTINUE
   80          CONTINUE
               TM8 = SECOND( )
               CALL STRMM( SIDE( OP1 ), UPLO( OP2 ), TRNS( OP3 ),
     $                         DIAG( OP4 ), DIM1( D ), DIM2( D ), ALPHA,
     $                                        A, LDA( L ), C, LDA( L ) )
               TM9 = SECOND( )
            END IF
            IF( TABSUB( 5 ) )THEN
               DO 100, J = 1, NMAX
                  DO 90, I = 1, LD
                     C( I, J ) = ONE + 0.01E+0*REAL( I+( J-1 )*NMAX )/
     $                                                 REAL( LD*NMAX+1 )
   90             CONTINUE
  100          CONTINUE
               TM10 = SECOND( )
               CALL STRSM( SIDE( OP1 ), UPLO( OP2 ), TRNS( OP3 ),
     $                         DIAG( OP4 ), DIM1( D ), DIM2( D ), ALPHA,
     $                                        A, LDA( L ), C, LDA( L ) )
               TM11 = SECOND( )
            END IF
         ELSE IF( LSAME( SB3LIB, 'G' ) )THEN
*
*           Time the built-in GEMM-Based Level 3 BLAS library (SGB02,
*           SGB04, SGB06, SGB08, and SGB09). Re-initialize C between the
*           timings to avoid overflow and to flush the cache.
*
            IF( TABSUB( 1 ).AND.OP3.EQ.1.AND.OP4.EQ.1 )THEN
               DO 120, J = 1, NMAX
                  DO 110, I = 1, LD
                     C( I, J ) = ONE + 0.01E+0*REAL( I+( J-1 )*NMAX )/
     $                                                 REAL( LD*NMAX+1 )
  110             CONTINUE
  120          CONTINUE
               TM2 = SECOND( )
               CALL SGB02( SIDE( OP1 ), UPLO( OP2 ), DIM1( D ),
     $                       DIM2( D ), ALPHA, A, LDA( L ), B, LDA( L ),
     $                                               BETA, C, LDA( L ) )
               TM3 = SECOND( )
            END IF
            IF( TABSUB( 2 ).AND.OP1.EQ.1.AND.OP4.EQ.1 )THEN
               DO 140, J = 1, NMAX
                  DO 130, I = 1, LD
                     C( I, J ) = ONE + 0.01E+0*REAL( I+( J-1 )*NMAX )/
     $                                                 REAL( LD*NMAX+1 )
  130             CONTINUE
  140          CONTINUE
               TM4 = SECOND( )
               CALL SGB04( UPLO( OP2 ), TRNS( OP3 ), DIM1( D ),
     $                                    DIM2( D ), ALPHA, A, LDA( L ),
     $                                               BETA, C, LDA( L ) )
               TM5 = SECOND( )
            END IF
            IF( TABSUB( 3 ).AND.OP1.EQ.1.AND.OP4.EQ.1 )THEN
               DO 160, J = 1, NMAX
                  DO 150, I = 1, LD
                     C( I, J ) = ONE + 0.01E+0*REAL( I+( J-1 )*NMAX )/
     $                                                 REAL( LD*NMAX+1 )
  150             CONTINUE
  160          CONTINUE
               TM6 = SECOND( )
               CALL SGB06( UPLO( OP2 ), TRNS( OP3 ), DIM1( D ),
     $                       DIM2( D ), ALPHA, A, LDA( L ), B, LDA( L ),
     $                                               BETA, C, LDA( L ) )
               TM7 = SECOND( )
            END IF
            IF( TABSUB( 4 ) )THEN
               DO 180, J = 1, NMAX
                  DO 170, I = 1, LD
                     C( I, J ) = ONE + 0.01E+0*REAL( I+( J-1 )*NMAX )/
     $                                                 REAL( LD*NMAX+1 )
  170             CONTINUE
  180          CONTINUE
               TM8 = SECOND( )
               CALL SGB08( SIDE( OP1 ), UPLO( OP2 ), TRNS( OP3 ),
     $                         DIAG( OP4 ), DIM1( D ), DIM2( D ), ALPHA,
     $                                        A, LDA( L ), C, LDA( L ) )
               TM9 = SECOND( )
            END IF
            IF( TABSUB( 5 ) )THEN
               DO 200, J = 1, NMAX
                  DO 190, I = 1, LD
                     C( I, J ) = ONE + 0.01E+0*REAL( I+( J-1 )*NMAX )/
     $                                                 REAL( LD*NMAX+1 )
  190             CONTINUE
  200          CONTINUE
               TM10 = SECOND( )
               CALL SGB09( SIDE( OP1 ), UPLO( OP2 ), TRNS( OP3 ),
     $                         DIAG( OP4 ), DIM1( D ), DIM2( D ), ALPHA,
     $                                        A, LDA( L ), C, LDA( L ) )
               TM11 = SECOND( )
            END IF
         ELSE
            WRITE( NERR, FMT = * )
     $      'Error: Unknown Level 3 BLAS library choosen: ', SB3LIB, '.'
         END IF
*
*        Compute the performance of SSYMM in Mflops.
*
         IF( TABSUB( 1 ).AND.OP3.EQ.1.AND.OP4.EQ.1 )THEN
            TIME = ( TM3 - TM2 ) - ( TM1 - TM0 )
            M = DIM1( D )
            N = DIM2( D )
            IF( LSAME( SIDE( OP1 ), 'L' ) )THEN
               NOPS = ( 2*M + 1 )*M*N + MIN( M*N, ( M*( M+1 ) )/2 )
            ELSE IF( LSAME( SIDE( OP1 ), 'R' ) )THEN
               NOPS = ( 2*N + 1 )*M*N + MIN( M*N, ( N*( N+1 ) )/2 )
            ELSE
               WRITE( NERR, FMT = * )
     $               'Error: Unknown value for SIDE: ', SIDE( OP1 ), '.'
            END IF
            IF( TIME.LE.ZERO )THEN
               SPEED = ZERO
            ELSE
               SPEED = REAL( NOPS )/( TIME*SCALE )
            END IF
            IF( RES( 1, OP1, OP2, 1, 1, D, L ).LT.SPEED )THEN
               RES( 1, OP1, OP2, 1, 1, D, L ) = SPEED
            END IF
         END IF
*
*        Compute the performance of SSYRK in Mflops.
*
         IF( TABSUB( 2 ).AND.OP1.EQ.1.AND.OP4.EQ.1 )THEN
            TIME = ( TM5 - TM4 ) - ( TM1 - TM0 )
            N = DIM1( D )
            K = DIM2( D )
            NOPS = ( 2*K + 1 )*( N*( N+1 )/2 ) + MIN( N*K, N*( N+1 )/2 )
            IF( TIME.LE.ZERO )THEN
               SPEED = ZERO
            ELSE
               SPEED = REAL( NOPS )/( TIME*SCALE )
            END IF
            IF( RES( 2, 1, OP2, OP3, 1, D, L ).LT.SPEED )THEN
               RES( 2, 1, OP2, OP3, 1, D, L ) = SPEED
            END IF
         END IF
*
*        Compute the performance of SSYR2K in Mflops.
*
         IF( TABSUB( 3 ).AND.OP1.EQ.1.AND.OP4.EQ.1 )THEN
            TIME = ( TM7 - TM6 ) - ( TM1 - TM0 )
            N = DIM1( D )
            K = DIM2( D )
            NOPS = ( 4*K + 1 )*( N*( N+1 )/2 ) + MIN( 2*N*K, N*( N+1 ) )
            IF( TIME.LE.ZERO )THEN
               SPEED = ZERO
            ELSE
               SPEED = REAL( NOPS )/( TIME*SCALE )
            END IF
            IF( RES( 3, 1, OP2, OP3, 1, D, L ).LT.SPEED )THEN
               RES( 3, 1, OP2, OP3, 1, D, L ) = SPEED
            END IF
         END IF
*
*        Compute the performance of STRMM in Mflops.
*
         IF( TABSUB( 4 ) )THEN
            TIME = ( TM9 - TM8 ) - ( TM1 - TM0 )
            M = DIM1( D )
            N = DIM2( D )
            IF( LSAME( SIDE( OP1 ), 'L' ) )THEN
               IF( LSAME( DIAG( OP4 ), 'N' ) )THEN
                  NOPS = M*M*N + MIN( M*N, ( M*( M + 1 ) )/2 )
               ELSE IF( LSAME( DIAG( OP4 ), 'U' ) )THEN
                  NOPS = M*M*N - M*N + MIN( M*N, ( M*( M + 1 ) )/2 )
               ELSE
                  WRITE( NERR, FMT = * )
     $               'Error: Unknown value for DIAG: ', DIAG( OP4 ), '.'
               END IF
            ELSE IF( LSAME( SIDE( OP1 ), 'R' ) )THEN
               IF( LSAME( DIAG( OP4 ), 'N' ) )THEN
                  NOPS = M*N*N + MIN( M*N, ( N*( N + 1 ) )/2 )
               ELSE IF( LSAME( DIAG( OP4 ), 'U' ) )THEN
                  NOPS = M*N*N - M*N + MIN( M*N, ( N*( N + 1 ) )/2 )
               ELSE
                  WRITE( NERR, FMT = * )
     $               'Error: Unknown value for DIAG: ', DIAG( OP4 ), '.'
               END IF
            ELSE
               WRITE( NERR, FMT = * )
     $               'Error: Unknown value for SIDE: ', SIDE( OP1 ), '.'
            END IF
            IF( TIME.LE.ZERO )THEN
               SPEED = ZERO
            ELSE
               SPEED = REAL( NOPS )/( TIME*SCALE )
            END IF
            IF( RES( 4, OP1, OP2, OP3, OP4, D, L ).LT.SPEED )THEN
               RES( 4, OP1, OP2, OP3, OP4, D, L ) = SPEED
            END IF
         END IF
*
*        Compute the performance of STRSM in Mflops.
*
         IF( TABSUB( 5 ) )THEN
            TIME = ( TM11 - TM10 ) - ( TM1 - TM0 )
            M = DIM1( D )
            N = DIM2( D )
            IF( LSAME( SIDE( OP1 ), 'L' ) )THEN
               IF( LSAME( DIAG( OP4 ), 'N' ) )THEN
                  NOPS = M*M*N + MIN( M*N, ( M*( M + 1 ) )/2 )
               ELSE IF( LSAME( DIAG( OP4 ), 'U' ) )THEN
                  NOPS = M*M*N - M*N + MIN( M*N, ( M*( M + 1 ) )/2 )
               ELSE
                  WRITE( NERR, FMT = * )
     $               'Error: Unknown value for DIAG: ', DIAG( OP4 ), '.'
               END IF
            ELSE IF( LSAME( SIDE( OP1 ), 'R' ) )THEN
               IF( LSAME( DIAG( OP4 ), 'N' ) )THEN
                  NOPS = M*N*N + MIN( M*N, ( N*( N + 1 ) )/2 )
               ELSE IF( LSAME( DIAG( OP4 ), 'U' ) )THEN
                  NOPS = M*N*N - M*N + MIN( M*N, ( N*( N + 1 ) )/2 )
               ELSE
                  WRITE( NERR, FMT = * )
     $               'Error: Unknown value for DIAG: ', DIAG( OP4 ), '.'
               END IF
            ELSE
               WRITE( NERR, FMT = * )
     $               'Error: Unknown value for SIDE: ', SIDE( OP1 ), '.'
            END IF
            IF( TIME.LE.ZERO )THEN
               SPEED = ZERO
            ELSE
               SPEED = REAL( NOPS )/( TIME*SCALE )
            END IF
            IF( RES( 5, OP1, OP2, OP3, OP4, D, L ).LT.SPEED )THEN
               RES( 5, OP1, OP2, OP3, OP4, D, L ) = SPEED
            END IF
         END IF
  210 CONTINUE
*
*     ------ Stop indentation ------
*
  220 CONTINUE
  230 CONTINUE
  240 CONTINUE
  250 CONTINUE
  260 CONTINUE
  270 CONTINUE
*
*     ------ Continue indentation ------
*
      RETURN
*
*     End of SGBT01.
*
      END
