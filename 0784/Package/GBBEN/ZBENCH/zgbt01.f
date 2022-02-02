      SUBROUTINE ZGBT01( ZB3LIB, TABSUB, SIDE, NSIDE, UPLO, NUPLO, TRNS,
     $                  NTRNS, DIAG, NDIAG, DIM1, DIM2, NDIM, LDA, NLDA,
     $                      ALPHA, BETA, A, B, C, LD, NMAX, NERR, MXSUB,
     $                          MXOPT, MXTRNS, MXDIM, MXLDA, RUNS, RES )
*     .. Scalar Arguments ..
      CHARACTER          ZB3LIB
      INTEGER            LD, NMAX, NERR,
     $                   NSIDE, NUPLO, NTRNS, NDIAG, NDIM, NLDA,
     $                   MXSUB, MXOPT, MXTRNS, MXDIM, MXLDA, RUNS
      COMPLEX*16         ALPHA, BETA
*     .. Array Arguments ..
      LOGICAL            TABSUB( MXSUB )
      CHARACTER          SIDE( MXOPT ), UPLO( MXOPT ), TRNS( MXTRNS ),
     $                   DIAG( MXOPT )
      INTEGER            DIM1( MXDIM ), DIM2( MXDIM ), LDA( MXLDA )
      COMPLEX*16         A( LD, NMAX ), B( LD, NMAX ), C( LD, NMAX )
      DOUBLE PRECISION   RES( MXSUB, MXOPT, MXOPT, MXTRNS, MXOPT, MXDIM,
     $                   MXLDA )
*
*
*  Time all routines except ZGEMM in the Level 3 BLAS library specified
*  by the input parameters. The library is either a user-supplied
*  Level 3 BLAS library or the GEMM-Based Level 3 BLAS library included
*  in the benchmark (ZGB02, ZGB04, ZGB04, ZGB05, ZGB06, ZGB07, ZGB08,
*  and ZGB09). Return the performance in Mflops for each problem
*  configuration.
*
*  ZGBT01 calls a DOUBLE PRECISION function DSECND with no arguments,
*  which is assumed to return the user time for a process in seconds
*  from some fixed starting-time.
*
*
*  -- Written in August-1994.
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
      INTEGER            I, J, M, N, K, ADDS, MULTS, NOPS,
     $                   D, L, R, OP1, OP2, OP3, OP4
      DOUBLE PRECISION   TIME, SPEED, TM0, TM1, TM2, TM3, TM4, TM5, TM6,
     $                   TM7, TM8, TM9, TM10, TM11, TM12, TM13, TM14,
     $                   TM15, TM16, TM17
*     .. Intrinsic Functions ..
      INTRINSIC          DBLE, DCMPLX, MIN
*     .. External Functions ..
      LOGICAL            LSAME
      DOUBLE PRECISION   DSECND
      EXTERNAL           LSAME, DSECND
*     .. External Subroutines ..
      EXTERNAL           ZSYMM, ZHEMM, ZSYRK, ZHERK, ZSYR2K, ZHER2K,
     $                   ZTRMM, ZTRSM,
     $                   ZGB02, ZGB03, ZGB04, ZGB05, ZGB06, ZGB07,
     $                   ZGB08, ZGB09
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, SCALE
      COMPLEX*16         Z11
*     .. Parameter Values ..
      PARAMETER        ( ZERO = 0.0D+0, SCALE = 1.0D+6,
     $                   Z11 = ( 1.0D+0, 1.0D+0 ) )
*     ..
*     .. Executable Statements ..
      TM0 = DSECND( )
      TM0 = DSECND( )
      TM0 = DSECND( )
      TM1 = DSECND( )
*
*     ------ Stop indentation ------
*
      DO 390, L = 1, NLDA
      DO 380, OP1 = 1, NSIDE
      DO 370, OP2 = 1, NUPLO
      DO 360, OP3 = 1, NTRNS
      DO 350, OP4 = 1, NDIAG
      DO 340, D = 1, NDIM
*
*     ------ Continue indentation ------
*
      RES( 1, OP1, OP2, OP3, OP4, D, L ) = ZERO
      RES( 2, OP1, OP2, OP3, OP4, D, L ) = ZERO
      RES( 3, OP1, OP2, OP3, OP4, D, L ) = ZERO
      RES( 4, OP1, OP2, OP3, OP4, D, L ) = ZERO
      RES( 5, OP1, OP2, OP3, OP4, D, L ) = ZERO
      DO 330, R = 1, RUNS
         IF( LSAME( ZB3LIB, 'U' ) )THEN
*
*           Time the user-supplied library. Re-initialize C between the
*           timings to avoid overflow and to flush the cache.
*
            IF( TABSUB( 1 ).AND.OP3.EQ.1.AND.OP4.EQ.1 )THEN
               DO 20, J = 1, NMAX
                  DO 10, I = 1, LD
                     C( I, J ) = Z11 + DCMPLX( 0.01D+0*
     $                      DBLE( I+( J-1 )*NMAX )/ DBLE( NMAX*NMAX+1 ),
     $                                   0.06D+0*DBLE( I+( J-1 )*NMAX )/
     $                                             DBLE( NMAX*NMAX+1 ) )
   10             CONTINUE
   20          CONTINUE
               TM2 = DSECND( )
               CALL ZSYMM( SIDE( OP1 ), UPLO( OP2 ), DIM1( D ),
     $                       DIM2( D ), ALPHA, A, LDA( L ), B, LDA( L ),
     $                                               BETA, C, LDA( L ) )
               TM3 = DSECND( )
            END IF
            IF( TABSUB( 2 ).AND.OP3.EQ.1.AND.OP4.EQ.1 )THEN
               DO 40, J = 1, NMAX
                  DO 30, I = 1, LD
                     C( I, J ) = Z11 + DCMPLX( 0.01D+0*
     $                      DBLE( I+( J-1 )*NMAX )/ DBLE( NMAX*NMAX+1 ),
     $                                   0.06D+0*DBLE( I+( J-1 )*NMAX )/
     $                                             DBLE( NMAX*NMAX+1 ) )
   30             CONTINUE
   40          CONTINUE
               TM4 = DSECND( )
               CALL ZHEMM( SIDE( OP1 ), UPLO( OP2 ), DIM1( D ),
     $                       DIM2( D ), ALPHA, A, LDA( L ), B, LDA( L ),
     $                                               BETA, C, LDA( L ) )
               TM5 = DSECND( )
            END IF
            IF( TABSUB( 3 ).AND.OP1.EQ.1.AND.OP4.EQ.1.AND.
     $                              .NOT.LSAME( TRNS( OP3 ), 'C' ) )THEN
               DO 60, J = 1, NMAX
                  DO 50, I = 1, LD
                     C( I, J ) = Z11 + DCMPLX( 0.01D+0*
     $                      DBLE( I+( J-1 )*NMAX )/ DBLE( NMAX*NMAX+1 ),
     $                                   0.06D+0*DBLE( I+( J-1 )*NMAX )/
     $                                             DBLE( NMAX*NMAX+1 ) )
   50             CONTINUE
   60          CONTINUE
               TM6 = DSECND( )
               CALL ZSYRK( UPLO( OP2 ), TRNS( OP3 ), DIM1( D ),
     $                                    DIM2( D ), ALPHA, A, LDA( L ),
     $                                               BETA, C, LDA( L ) )
               TM7 = DSECND( )
            END IF
            IF( TABSUB( 4 ).AND.OP1.EQ.1.AND.OP4.EQ.1.AND.
     $                              .NOT.LSAME( TRNS( OP3 ), 'T' ) )THEN
               DO 80, J = 1, NMAX
                  DO 70, I = 1, LD
                     C( I, J ) = Z11 + DCMPLX( 0.01D+0*
     $                      DBLE( I+( J-1 )*NMAX )/ DBLE( NMAX*NMAX+1 ),
     $                                   0.06D+0*DBLE( I+( J-1 )*NMAX )/
     $                                             DBLE( NMAX*NMAX+1 ) )
   70             CONTINUE
   80          CONTINUE
               TM8 = DSECND( )
               CALL ZHERK( UPLO( OP2 ), TRNS( OP3 ), DIM1( D ),
     $                            DIM2( D ), DBLE( ALPHA ), A, LDA( L ),
     $                                       DBLE( BETA ), C, LDA( L ) )
               TM9 = DSECND( )
            END IF
            IF( TABSUB( 5 ).AND.OP1.EQ.1.AND.OP4.EQ.1.AND.
     $                              .NOT.LSAME( TRNS( OP3 ), 'C' ) )THEN
               DO 100, J = 1, NMAX
                  DO 90, I = 1, LD
                     C( I, J ) = Z11 + DCMPLX( 0.01D+0*
     $                      DBLE( I+( J-1 )*NMAX )/ DBLE( NMAX*NMAX+1 ),
     $                                   0.06D+0*DBLE( I+( J-1 )*NMAX )/
     $                                             DBLE( NMAX*NMAX+1 ) )
   90             CONTINUE
  100          CONTINUE
               TM10 = DSECND( )
               CALL ZSYR2K( UPLO( OP2 ), TRNS( OP3 ), DIM1( D ),
     $                       DIM2( D ), ALPHA, A, LDA( L ), B, LDA( L ),
     $                                               BETA, C, LDA( L ) )
               TM11 = DSECND( )
            END IF
            IF( TABSUB( 6 ).AND.OP1.EQ.1.AND.OP4.EQ.1.AND.
     $                              .NOT.LSAME( TRNS( OP3 ), 'T' ) )THEN
               DO 120, J = 1, NMAX
                  DO 110, I = 1, LD
                     C( I, J ) = Z11 + DCMPLX( 0.01D+0*
     $                      DBLE( I+( J-1 )*NMAX )/ DBLE( NMAX*NMAX+1 ),
     $                                   0.06D+0*DBLE( I+( J-1 )*NMAX )/
     $                                             DBLE( NMAX*NMAX+1 ) )
  110             CONTINUE
  120          CONTINUE
               TM12 = DSECND( )
               CALL ZHER2K( UPLO( OP2 ), TRNS( OP3 ), DIM1( D ),
     $                       DIM2( D ), ALPHA, A, LDA( L ), B, LDA( L ),
     $                                       DBLE( BETA ), C, LDA( L ) )
               TM13 = DSECND( )
            END IF
            IF( TABSUB( 7 ) )THEN
               DO 140, J = 1, NMAX
                  DO 130, I = 1, LD
                     C( I, J ) = Z11 + DCMPLX( 0.01D+0*
     $                      DBLE( I+( J-1 )*NMAX )/ DBLE( NMAX*NMAX+1 ),
     $                                   0.06D+0*DBLE( I+( J-1 )*NMAX )/
     $                                             DBLE( NMAX*NMAX+1 ) )
  130             CONTINUE
  140          CONTINUE
               TM14 = DSECND( )
               CALL ZTRMM( SIDE( OP1 ), UPLO( OP2 ), TRNS( OP3 ),
     $                         DIAG( OP4 ), DIM1( D ), DIM2( D ), ALPHA,
     $                                        A, LDA( L ), C, LDA( L ) )
               TM15 = DSECND( )
            END IF
            IF( TABSUB( 8 ) )THEN
               DO 160, J = 1, NMAX
                  DO 150, I = 1, LD
                     C( I, J ) = Z11 + DCMPLX( 0.01D+0*
     $                      DBLE( I+( J-1 )*NMAX )/ DBLE( NMAX*NMAX+1 ),
     $                                   0.06D+0*DBLE( I+( J-1 )*NMAX )/
     $                                             DBLE( NMAX*NMAX+1 ) )
  150             CONTINUE
  160          CONTINUE
               TM16 = DSECND( )
               CALL ZTRSM( SIDE( OP1 ), UPLO( OP2 ), TRNS( OP3 ),
     $                         DIAG( OP4 ), DIM1( D ), DIM2( D ), ALPHA,
     $                                        A, LDA( L ), C, LDA( L ) )
               TM17 = DSECND( )
            END IF
         ELSE IF( LSAME( ZB3LIB, 'G' ) )THEN
*
*           Time the built-in GEMM-Based Level 3 BLAS library (DGB02,
*           DGB04, DGB06, DGB08, and DGB09). Re-initialize C between the
*           timings to avoid overflow and to flush the cache.
*
            IF( TABSUB( 1 ).AND.OP3.EQ.1.AND.OP4.EQ.1 )THEN
               DO 180, J = 1, NMAX
                  DO 170, I = 1, LD
                     C( I, J ) = Z11 + DCMPLX( 0.01D+0*
     $                      DBLE( I+( J-1 )*NMAX )/ DBLE( NMAX*NMAX+1 ),
     $                                   0.06D+0*DBLE( I+( J-1 )*NMAX )/
     $                                             DBLE( NMAX*NMAX+1 ) )
  170             CONTINUE
  180          CONTINUE
               TM2 = DSECND( )
               CALL ZGB02( SIDE( OP1 ), UPLO( OP2 ), DIM1( D ),
     $                       DIM2( D ), ALPHA, A, LDA( L ), B, LDA( L ),
     $                                               BETA, C, LDA( L ) )
               TM3 = DSECND( )
            END IF
            IF( TABSUB( 2 ).AND.OP3.EQ.1.AND.OP4.EQ.1 )THEN
               DO 200, J = 1, NMAX
                  DO 190, I = 1, LD
                     C( I, J ) = Z11 + DCMPLX( 0.01D+0*
     $                      DBLE( I+( J-1 )*NMAX )/ DBLE( NMAX*NMAX+1 ),
     $                                   0.06D+0*DBLE( I+( J-1 )*NMAX )/
     $                                             DBLE( NMAX*NMAX+1 ) )
  190             CONTINUE
  200          CONTINUE
               TM4 = DSECND( )
               CALL ZGB03( SIDE( OP1 ), UPLO( OP2 ), DIM1( D ),
     $                       DIM2( D ), ALPHA, A, LDA( L ), B, LDA( L ),
     $                                               BETA, C, LDA( L ) )
               TM5 = DSECND( )
            END IF
            IF( TABSUB( 3 ).AND.OP1.EQ.1.AND.OP4.EQ.1.AND.
     $                              .NOT.LSAME( TRNS( OP3 ), 'C' ) )THEN
               DO 220, J = 1, NMAX
                  DO 210, I = 1, LD
                     C( I, J ) = Z11 + DCMPLX( 0.01D+0*
     $                      DBLE( I+( J-1 )*NMAX )/ DBLE( NMAX*NMAX+1 ),
     $                                   0.06D+0*DBLE( I+( J-1 )*NMAX )/
     $                                             DBLE( NMAX*NMAX+1 ) )
  210             CONTINUE
  220          CONTINUE
               TM6 = DSECND( )
               CALL ZGB04( UPLO( OP2 ), TRNS( OP3 ), DIM1( D ),
     $                                    DIM2( D ), ALPHA, A, LDA( L ),
     $                                               BETA, C, LDA( L ) )
               TM7 = DSECND( )
            END IF
            IF( TABSUB( 4 ).AND.OP1.EQ.1.AND.OP4.EQ.1.AND.
     $                              .NOT.LSAME( TRNS( OP3 ), 'T' ) )THEN
               DO 240, J = 1, NMAX
                  DO 230, I = 1, LD
                     C( I, J ) = Z11 + DCMPLX( 0.01D+0*
     $                      DBLE( I+( J-1 )*NMAX )/ DBLE( NMAX*NMAX+1 ),
     $                                   0.06D+0*DBLE( I+( J-1 )*NMAX )/
     $                                             DBLE( NMAX*NMAX+1 ) )
  230             CONTINUE
  240          CONTINUE
               TM8 = DSECND( )
               CALL ZGB05( UPLO( OP2 ), TRNS( OP3 ), DIM1( D ),
     $                            DIM2( D ), DBLE( ALPHA ), A, LDA( L ),
     $                                       DBLE( BETA ), C, LDA( L ) )
               TM9 = DSECND( )
            END IF
            IF( TABSUB( 5 ).AND.OP1.EQ.1.AND.OP4.EQ.1.AND.
     $                              .NOT.LSAME( TRNS( OP3 ), 'C' ) )THEN
               DO 260, J = 1, NMAX
                  DO 250, I = 1, LD
                     C( I, J ) = Z11 + DCMPLX( 0.01D+0*
     $                      DBLE( I+( J-1 )*NMAX )/ DBLE( NMAX*NMAX+1 ),
     $                                   0.06D+0*DBLE( I+( J-1 )*NMAX )/
     $                                             DBLE( NMAX*NMAX+1 ) )
  250             CONTINUE
  260          CONTINUE
               TM10 = DSECND( )
               CALL ZGB06( UPLO( OP2 ), TRNS( OP3 ), DIM1( D ),
     $                       DIM2( D ), ALPHA, A, LDA( L ), B, LDA( L ),
     $                                               BETA, C, LDA( L ) )
               TM11 = DSECND( )
            END IF
            IF( TABSUB( 6 ).AND.OP1.EQ.1.AND.OP4.EQ.1.AND.
     $                              .NOT.LSAME( TRNS( OP3 ), 'T' ) )THEN
               DO 280, J = 1, NMAX
                  DO 270, I = 1, LD
                     C( I, J ) = Z11 + DCMPLX( 0.01D+0*
     $                      DBLE( I+( J-1 )*NMAX )/ DBLE( NMAX*NMAX+1 ),
     $                                   0.06D+0*DBLE( I+( J-1 )*NMAX )/
     $                                             DBLE( NMAX*NMAX+1 ) )
  270             CONTINUE
  280          CONTINUE
               TM12 = DSECND( )
               CALL ZGB07( UPLO( OP2 ), TRNS( OP3 ), DIM1( D ),
     $                       DIM2( D ), ALPHA, A, LDA( L ), B, LDA( L ),
     $                                       DBLE( BETA ), C, LDA( L ) )
               TM13 = DSECND( )
            END IF
            IF( TABSUB( 7 ) )THEN
               DO 300, J = 1, NMAX
                  DO 290, I = 1, LD
                     C( I, J ) = Z11 + DCMPLX( 0.01D+0*
     $                      DBLE( I+( J-1 )*NMAX )/ DBLE( NMAX*NMAX+1 ),
     $                                   0.06D+0*DBLE( I+( J-1 )*NMAX )/
     $                                             DBLE( NMAX*NMAX+1 ) )
  290             CONTINUE
  300          CONTINUE
               TM14 = DSECND( )
               CALL ZGB08( SIDE( OP1 ), UPLO( OP2 ), TRNS( OP3 ),
     $                         DIAG( OP4 ), DIM1( D ), DIM2( D ), ALPHA,
     $                                        A, LDA( L ), C, LDA( L ) )
               TM15 = DSECND( )
            END IF
            IF( TABSUB( 8 ) )THEN
               DO 320, J = 1, NMAX
                  DO 310, I = 1, LD
                     C( I, J ) = Z11 + DCMPLX( 0.01D+0*
     $                      DBLE( I+( J-1 )*NMAX )/ DBLE( NMAX*NMAX+1 ),
     $                                   0.06D+0*DBLE( I+( J-1 )*NMAX )/
     $                                             DBLE( NMAX*NMAX+1 ) )
  310             CONTINUE
  320          CONTINUE
               TM16 = DSECND( )
               CALL ZGB09( SIDE( OP1 ), UPLO( OP2 ), TRNS( OP3 ),
     $                         DIAG( OP4 ), DIM1( D ), DIM2( D ), ALPHA,
     $                                        A, LDA( L ), C, LDA( L ) )
               TM17 = DSECND( )
            END IF
         ELSE
            WRITE( NERR, FMT = * )
     $      'Error: Unknown Level 3 BLAS library choosen: ', ZB3LIB, '.'
         END IF
*
*        Compute the performance of ZSYMM in Mflops.
*
         IF( TABSUB( 1 ).AND.OP3.EQ.1.AND.OP4.EQ.1 )THEN
            TIME = ( TM3 - TM2 ) - ( TM1 - TM0 )
            M = DIM1( D )
            N = DIM2( D )
            IF( LSAME( SIDE( OP1 ), 'L' ) )THEN
               MULTS = ( M + 1 )*M*N + MIN( M*N, ( M*( M+1 ) )/2 )
               ADDS = M*M*N
               NOPS = 6*MULTS + 2*ADDS
            ELSE IF( LSAME( SIDE( OP1 ), 'R' ) )THEN
               MULTS = ( N + 1 )*M*N + MIN( M*N, ( N*( N+1 ) )/2 )
               ADDS = M*N*N
               NOPS = 6*MULTS + 2*ADDS
            ELSE
               WRITE( NERR, FMT = * )
     $               'Error: Unknown value for SIDE: ', SIDE( OP1 ), '.'
            END IF
            IF( TIME.LE.ZERO )THEN
               SPEED = ZERO
            ELSE
               SPEED = DBLE( NOPS )/( TIME*SCALE )
            END IF
            IF( RES( 1, OP1, OP2, 1, 1, D, L ).LT.SPEED )THEN
               RES( 1, OP1, OP2, 1, 1, D, L ) = SPEED
            END IF
         END IF
*
*        Compute the performance of ZHEMM in Mflops.
*
         IF( TABSUB( 2 ).AND.OP3.EQ.1.AND.OP4.EQ.1 )THEN
            TIME = ( TM5 - TM4 ) - ( TM1 - TM0 )
            M = DIM1( D )
            N = DIM2( D )
            IF( LSAME( SIDE( OP1 ), 'L' ) )THEN
               MULTS = ( 6*M + 2 )*M*N + MIN( 6*M*N, 3*M*M - M )
               ADDS = 2*M*M*N
               NOPS = MULTS + ADDS
            ELSE IF( LSAME( SIDE( OP1 ), 'R' ) )THEN
               MULTS = ( 6*N + 2 )*M*N + MIN( 6*M*N, 3*N*N - N )
               ADDS = 2*M*N*N
               NOPS = MULTS + ADDS
            ELSE
               WRITE( NERR, FMT = * )
     $               'Error: Unknown value for SIDE: ', SIDE( OP1 ), '.'
            END IF
            IF( TIME.LE.ZERO )THEN
               SPEED = ZERO
            ELSE
               SPEED = DBLE( NOPS )/( TIME*SCALE )
            END IF
            IF( RES( 2, OP1, OP2, 1, 1, D, L ).LT.SPEED )THEN
               RES( 2, OP1, OP2, 1, 1, D, L ) = SPEED
            END IF
         END IF
*
*        Compute the performance of ZSYRK in Mflops.
*
         IF( TABSUB( 3 ).AND.OP1.EQ.1.AND.OP4.EQ.1.AND.
     $                              .NOT.LSAME( TRNS( OP3 ), 'C' ) )THEN
            TIME = ( TM7 - TM6 ) - ( TM1 - TM0 )
            N = DIM1( D )
            K = DIM2( D )
            MULTS = ( K + 1 )*( N*( N+1 )/2 ) + MIN( N*K, N*( N+1 )/2 )
            ADDS = K*( N*( N+1 )/2 )
            NOPS = 6*MULTS + 2*ADDS
            IF( TIME.LE.ZERO )THEN
               SPEED = ZERO
            ELSE
               SPEED = DBLE( NOPS )/( TIME*SCALE )
            END IF
            IF( RES( 3, 1, OP2, OP3, 1, D, L ).LT.SPEED )THEN
               RES( 3, 1, OP2, OP3, 1, D, L ) = SPEED
            END IF
         END IF
*
*        Compute the performance of ZHERK in Mflops.
*
         IF( TABSUB( 4 ).AND.OP1.EQ.1.AND.OP4.EQ.1.AND.
     $                              .NOT.LSAME( TRNS( OP3 ), 'T' ) )THEN
            TIME = ( TM9 - TM8 ) - ( TM1 - TM0 )
            N = DIM1( D )
            K = DIM2( D )
            MULTS = ( 3*K + 1 )*N*N + MIN( 2*N*K, N*N )
            ADDS = K*N*N
            NOPS = MULTS + ADDS
            IF( TIME.LE.ZERO )THEN
               SPEED = ZERO
            ELSE
               SPEED = DBLE( NOPS )/( TIME*SCALE )
            END IF
            IF( RES( 4, 1, OP2, OP3, 1, D, L ).LT.SPEED )THEN
               RES( 4, 1, OP2, OP3, 1, D, L ) = SPEED
            END IF
         END IF
*
*        Compute the performance of ZSYR2K in Mflops.
*
         IF( TABSUB( 5 ).AND.OP1.EQ.1.AND.OP4.EQ.1.AND.
     $                              .NOT.LSAME( TRNS( OP3 ), 'C' ) )THEN
            TIME = ( TM11 - TM10 ) - ( TM1 - TM0 )
            N = DIM1( D )
            K = DIM2( D )
            MULTS = ( 2*K + 1 )*( N*( N+1 )/2 ) +
     $                                           MIN( 2*N*K, N*( N+1 ) )
            ADDS = K*N*( N+1 )
            NOPS = 6*MULTS + 2*ADDS
            IF( TIME.LE.ZERO )THEN
               SPEED = ZERO
            ELSE
               SPEED = DBLE( NOPS )/( TIME*SCALE )
            END IF
            IF( RES( 5, 1, OP2, OP3, 1, D, L ).LT.SPEED )THEN
               RES( 5, 1, OP2, OP3, 1, D, L ) = SPEED
            END IF
         END IF
*
*        Compute the performance of ZHER2K in Mflops.
*
         IF( TABSUB( 6 ).AND.OP1.EQ.1.AND.OP4.EQ.1.AND.
     $                              .NOT.LSAME( TRNS( OP3 ), 'T' ) )THEN
            TIME = ( TM13 - TM12 ) - ( TM1 - TM0 )
            N = DIM1( D )
            K = DIM2( D )
            MULTS = ( 6*K + 1 )*N*N + MIN( 12*N*K, 6*N*N - 2*N )
            ADDS = 2*K*N*N
            NOPS = MULTS + ADDS
            IF( TIME.LE.ZERO )THEN
               SPEED = ZERO
            ELSE
               SPEED = DBLE( NOPS )/( TIME*SCALE )
            END IF
            IF( RES( 6, 1, OP2, OP3, 1, D, L ).LT.SPEED )THEN
               RES( 6, 1, OP2, OP3, 1, D, L ) = SPEED
            END IF
         END IF
*
*        Compute the performance of ZTRMM in Mflops.
*
         IF( TABSUB( 7 ) )THEN
            TIME = ( TM15 - TM14 ) - ( TM1 - TM0 )
            M = DIM1( D )
            N = DIM2( D )
            IF( LSAME( SIDE( OP1 ), 'L' ) )THEN
               IF( LSAME( DIAG( OP4 ), 'N' ) )THEN
                  MULTS = ( ( M*( M + 1 ) )/2 )*N +
     $                                     MIN( M*N, ( M*( M + 1 ) )/2 )
                  ADDS = ( ( M*( M - 1 ) )/2 )*N
                  NOPS = 6*MULTS + 2*ADDS
               ELSE IF( LSAME( DIAG( OP4 ), 'U' ) )THEN
                  MULTS = ( ( M*( M - 1 ) )/2 )*N +
     $                                     MIN( M*N, ( M*( M + 1 ) )/2 )
                  ADDS = ( ( M*( M - 1 ) )/2 )*N
                  NOPS = 6*MULTS + 2*ADDS
               ELSE
                  WRITE( NERR, FMT = * )
     $               'Error: Unknown value for DIAG: ', DIAG( OP4 ), '.'
               END IF
            ELSE IF( LSAME( SIDE( OP1 ), 'R' ) )THEN
               IF( LSAME( DIAG( OP4 ), 'N' ) )THEN
                  MULTS = ( M*( N*( N + 1 ) )/2 ) +
     $                                     MIN( M*N, ( N*( N + 1 ) )/2 )
                  ADDS = M*( N*( N - 1 ) )/2
                  NOPS = 6*MULTS + 2*ADDS
               ELSE IF( LSAME( DIAG( OP4 ), 'U' ) )THEN
                  MULTS = ( M*( N*( N - 1 ) )/2 ) +
     $                                     MIN( M*N, ( N*( N + 1 ) )/2 )
                  ADDS = M*( N*( N - 1 ) )/2
                  NOPS = 6*MULTS + 2*ADDS
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
               SPEED = DBLE( NOPS )/( TIME*SCALE )
            END IF
            IF( RES( 7, OP1, OP2, OP3, OP4, D, L ).LT.SPEED )THEN
               RES( 7, OP1, OP2, OP3, OP4, D, L ) = SPEED
            END IF
         END IF
*
*        Compute the performance of ZTRSM in Mflops.
*
         IF( TABSUB( 8 ) )THEN
            TIME = ( TM17 - TM16 ) - ( TM1 - TM0 )
            M = DIM1( D )
            N = DIM2( D )
            IF( LSAME( SIDE( OP1 ), 'L' ) )THEN
               IF( LSAME( DIAG( OP4 ), 'N' ) )THEN
                  MULTS = ( ( M*( M + 1 ) )/2 )*N +
     $                                     MIN( M*N, ( M*( M + 1 ) )/2 )
                  ADDS = ( ( M*( M - 1 ) )/2 )*N
                  NOPS = 6*MULTS + 2*ADDS
               ELSE IF( LSAME( DIAG( OP4 ), 'U' ) )THEN
                  MULTS = ( ( M*( M - 1 ) )/2 )*N +
     $                                     MIN( M*N, ( M*( M + 1 ) )/2 )
                  ADDS = ( ( M*( M - 1 ) )/2 )*N
                  NOPS = 6*MULTS + 2*ADDS
               ELSE
                  WRITE( NERR, FMT = * )
     $               'Error: Unknown value for DIAG: ', DIAG( OP4 ), '.'
               END IF
            ELSE IF( LSAME( SIDE( OP1 ), 'R' ) )THEN
               IF( LSAME( DIAG( OP4 ), 'N' ) )THEN
                  MULTS = ( M*( N*( N + 1 ) )/2 ) +
     $                                     MIN( M*N, ( N*( N + 1 ) )/2 )
                  ADDS = M*( N*( N - 1 ) )/2
                  NOPS = 6*MULTS + 2*ADDS
               ELSE IF( LSAME( DIAG( OP4 ), 'U' ) )THEN
                  MULTS = ( M*( N*( N - 1 ) )/2 ) +
     $                                     MIN( M*N, ( N*( N + 1 ) )/2 )
                  ADDS = M*( N*( N - 1 ) )/2
                  NOPS = 6*MULTS + 2*ADDS
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
               SPEED = DBLE( NOPS )/( TIME*SCALE )
            END IF
            IF( RES( 8, OP1, OP2, OP3, OP4, D, L ).LT.SPEED )THEN
               RES( 8, OP1, OP2, OP3, OP4, D, L ) = SPEED
            END IF
         END IF
  330 CONTINUE
*
*     ------ Stop indentation ------
*
  340 CONTINUE
  350 CONTINUE
  360 CONTINUE
  370 CONTINUE
  380 CONTINUE
  390 CONTINUE
*
*     ------ Continue indentation ------
*
      RETURN
*
*     End of ZGBT01.
*
      END
