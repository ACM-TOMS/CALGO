      SUBROUTINE CGBT02( TABSUB, SIDE, NSIDE, NUPLO, TRNS, NTRNS, NDIAG,
     $                         DIM1, DIM2, NDIM, LDA, NLDA, ALPHA, BETA,
     $                            A, B, C, LD, NMAX, NERR, MXSUB, MXOPT,
     $                                 MXTRNS, MXDIM, MXLDA, RUNS, RES )
*     .. Scalar Arguments ..
      INTEGER            LD, NMAX, NERR,
     $                   NSIDE, NUPLO, NTRNS, NDIAG, NDIM, NLDA,
     $                   MXSUB, MXOPT, MXTRNS, MXDIM, MXLDA, RUNS
      COMPLEX            ALPHA, BETA
*     .. Array Arguments ..
      LOGICAL            TABSUB( MXSUB )
      CHARACTER          SIDE( MXOPT ), TRNS( MXTRNS )
      INTEGER            DIM1( MXDIM ), DIM2( MXDIM ), LDA( MXLDA )
      COMPLEX            A( LD, NMAX ), B( LD, NMAX ), C( LD, NMAX )
      REAL               RES( MXSUB, MXOPT, MXOPT, MXTRNS, MXOPT, MXDIM,
     $                   MXLDA )
*
*
*  Determine problem configurations for CGEMM that, for timing purposes,
*  "correspond" to problem configurations for the remaining Level 3 BLAS
*  routines. Time CGEMM for problems that correspond to the Level 3 BLAS
*  problems timed in CGBT01. Return the performance of CGEMM in Mflops.
*
*  CGBT02 calls a REAL function SECOND with no arguments,
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
      REAL               TIME, SPEED, TM0, TM1, TM2, TM3, TM4, TM5, TM6,
     $                   TM7, TM8, TM9, TM10, TM11, TM12, TM13
*     .. Intrinsic Functions ..
      INTRINSIC          REAL, CMPLX, MIN
*     .. External Functions ..
      LOGICAL            LSAME
      REAL               SECOND
      EXTERNAL           LSAME, SECOND
*     .. External Subroutines ..
      EXTERNAL           CGEMM
*     .. Parameters ..
      REAL               ZERO, SCALE
      COMPLEX            ONE, C11
      PARAMETER        ( ZERO = 0.0E+0, SCALE = 1.0E+6,
     $                   ONE = ( 1.0E+0, 0.0E+0 ),
     $                   C11 = ( 1.0E+0, 1.0E+0 ) )
*     ..
*     .. Executable Statements ..
      TM0 = SECOND( )
      TM0 = SECOND( )
      TM0 = SECOND( )
      TM1 = SECOND( )
*
*     ------ Stop indentation ------
*
      DO 240, L = 1, NLDA
      DO 230, OP1 = 1, NSIDE
      DO 220, OP3 = 1, NTRNS
      DO 210, D = 1, NDIM
*
*     ------ Continue indentation ------
*
      RES( 1, OP1, 1, OP3, 1, D, L ) = ZERO
      RES( 2, OP1, 1, OP3, 1, D, L ) = ZERO
      RES( 3, OP1, 1, OP3, 1, D, L ) = ZERO
      RES( 4, OP1, 1, OP3, 1, D, L ) = ZERO
      RES( 5, OP1, 1, OP3, 1, D, L ) = ZERO
      DO 200, R = 1, RUNS
*
*        Time the user-supplied library. Re-initialize C between the
*        timings to avoid overflow and to flush the cache.
*
         IF( ( TABSUB( 1 ).OR.TABSUB( 2 ) ).AND.OP3.EQ.1 )THEN
            DO 20, J = 1, NMAX
               DO 10, I = 1, LD
                  C( I, J ) = C11 + CMPLX( 0.01E+0*
     $                      REAL( I+( J-1 )*NMAX )/ REAL( NMAX*NMAX+1 ),
     $                                   0.06E+0*REAL( I+( J-1 )*NMAX )/
     $                                             REAL( NMAX*NMAX+1 ) )
   10          CONTINUE
   20       CONTINUE
*
*           Time CGEMM for a problem that corresponds to the following
*           problem for CSYMM or CHEMM:
*           CSYMM( SIDE( OP1 ), UPLO( OP2 ), DIM1( D ), DIM2( D ),
*                   ALPHA, A, LDA( L ), B, LDA( L ), BETA, C, LDA( L ) )
*
            IF( LSAME( SIDE( OP1 ), 'L' ) )THEN
*
*              Use K = M.
*
               TM2 = SECOND( )
               CALL CGEMM( 'N', 'N', DIM1( D ), DIM2( D ), DIM1( D ),
     $              ALPHA, A, LDA( L ), B, LDA( L ), BETA, C, LDA( L ) )
               TM3 = SECOND( )
            ELSE IF( LSAME( SIDE( OP1 ), 'R' ) )THEN
*
*              Use K = N.
*
               TM2 = SECOND( )
               CALL CGEMM( 'N', 'N', DIM1( D ), DIM2( D ), DIM2( D ),
     $              ALPHA, B, LDA( L ), A, LDA( L ), BETA, C, LDA( L ) )
               TM3 = SECOND( )
            ELSE
               WRITE( NERR, FMT = * )
     $               'Error: Unknown value for SIDE: ', SIDE( OP1 ), '.'
               STOP
            END IF
         END IF
         IF( TABSUB( 3 ).AND.OP1.EQ.1.AND.
     $                              .NOT.LSAME( TRNS( OP3 ), 'C' ) )THEN
            DO 40, J = 1, NMAX
               DO 30, I = 1, LD
                  C( I, J ) = C11 + CMPLX( 0.01E+0*
     $                      REAL( I+( J-1 )*NMAX )/ REAL( NMAX*NMAX+1 ),
     $                                   0.06E+0*REAL( I+( J-1 )*NMAX )/
     $                                             REAL( NMAX*NMAX+1 ) )
   30          CONTINUE
   40       CONTINUE
*
*           Time CGEMM for a problem that corresponds to the following
*           problem for CSYRK:
*           CSYRK( UPLO( OP2 ), TRNS( OP3 ), DIM1( D ), DIM2( D ),
*                                ALPHA, A, LDA( L ), BETA, C, LDA( L ) )
*           Use M = N and B = A in the call to CGEMM.
*
            IF( LSAME( TRNS( OP3 ), 'N' ) )THEN
               TM4 = SECOND( )
               CALL CGEMM( 'N', 'T', DIM1( D ), DIM1( D ), DIM2( D ),
     $              ALPHA, A, LDA( L ), A, LDA( L ), BETA, C, LDA( L ) )
               TM5 = SECOND( )
            ELSE IF( LSAME( TRNS( OP3 ), 'T' ) )THEN
               TM4 = SECOND( )
               CALL CGEMM( 'T', 'N', DIM1( D ), DIM1( D ), DIM2( D ),
     $              ALPHA, A, LDA( L ), A, LDA( L ), BETA, C, LDA( L ) )
               TM5 = SECOND( )
            ELSE
               WRITE( NERR, FMT = * )
     $              'Error: Unknown value for TRANS: ', TRNS( OP3 ), '.'
               STOP
            END IF
         END IF
         IF( TABSUB( 4 ).AND.OP1.EQ.1.AND.
     $                              .NOT.LSAME( TRNS( OP3 ), 'T' ) )THEN
            DO 60, J = 1, NMAX
               DO 50, I = 1, LD
                  C( I, J ) = C11 + CMPLX( 0.01E+0*
     $                      REAL( I+( J-1 )*NMAX )/ REAL( NMAX*NMAX+1 ),
     $                                   0.06E+0*REAL( I+( J-1 )*NMAX )/
     $                                             REAL( NMAX*NMAX+1 ) )
   50          CONTINUE
   60       CONTINUE
*
*           Time CGEMM for a problem that corresponds to the following
*           problem for CHERK:
*           CHERK( UPLO( OP2 ), TRNS( OP3 ), DIM1( D ), DIM2( D ),
*                REAL( ALPHA ), A, LDA( L ), REAL( BETA ), C, LDA( L ) )
*           Use M = N and B = A in the call to CGEMM.
*
            IF( LSAME( TRNS( OP3 ), 'N' ) )THEN
               TM6 = SECOND( )
               CALL CGEMM( 'N', 'C', DIM1( D ), DIM1( D ), DIM2( D ),
     $              ALPHA, A, LDA( L ), A, LDA( L ), BETA, C, LDA( L ) )
               TM7 = SECOND( )
            ELSE IF( LSAME( TRNS( OP3 ), 'C' ) )THEN
               TM6 = SECOND( )
               CALL CGEMM( 'C', 'N', DIM1( D ), DIM1( D ), DIM2( D ),
     $              ALPHA, A, LDA( L ), A, LDA( L ), BETA, C, LDA( L ) )
               TM7 = SECOND( )
            ELSE
               WRITE( NERR, FMT = * )
     $              'Error: Unknown value for TRANS: ', TRNS( OP3 ), '.'
               STOP
            END IF
         END IF
         IF( TABSUB( 5 ).AND.OP1.EQ.1.AND.
     $                              .NOT.LSAME( TRNS( OP3 ), 'C' ) )THEN
            DO 80, J = 1, NMAX
               DO 70, I = 1, LD
                  C( I, J ) = C11 + CMPLX( 0.01E+0*
     $                      REAL( I+( J-1 )*NMAX )/ REAL( NMAX*NMAX+1 ),
     $                                   0.06E+0*REAL( I+( J-1 )*NMAX )/
     $                                             REAL( NMAX*NMAX+1 ) )
   70          CONTINUE
   80       CONTINUE
*
*           Time CGEMM for a problem that corresponds to the following
*           problem for CSYR2K:
*           CSYR2K( UPLO( OP2 ), TRNS( OP3 ), DIM1( D ), DIM2( D ),
*                   ALPHA, A, LDA( L ), B, LDA( L ), BETA, C, LDA( L ) )
*
            IF( LSAME( TRNS( OP3 ), 'N' ) )THEN
               TM8 = SECOND( )
               CALL CGEMM( 'N', 'T', DIM1( D ), DIM1( D ), DIM2( D ),
     $              ALPHA, A, LDA( L ), B, LDA( L ), BETA, C, LDA( L ) )
               TM9 = SECOND( )
            ELSE IF( LSAME( TRNS( OP3 ), 'T' ) )THEN
               TM8 = SECOND( )
               CALL CGEMM( 'T', 'N', DIM1( D ), DIM1( D ), DIM2( D ),
     $              ALPHA, A, LDA( L ), B, LDA( L ), BETA, C, LDA( L ) )
               TM9 = SECOND( )
            ELSE
               WRITE( NERR, FMT = * )
     $              'Error: Unknown value for TRANS: ', TRNS( OP3 ), '.'
               STOP
            END IF
         END IF
         IF( TABSUB( 6 ).AND.OP1.EQ.1.AND.
     $                              .NOT.LSAME( TRNS( OP3 ), 'T' ) )THEN
            DO 100, J = 1, NMAX
               DO 90, I = 1, LD
                  C( I, J ) = C11 + CMPLX( 0.01E+0*
     $                      REAL( I+( J-1 )*NMAX )/ REAL( NMAX*NMAX+1 ),
     $                                   0.06E+0*REAL( I+( J-1 )*NMAX )/
     $                                             REAL( NMAX*NMAX+1 ) )
   90          CONTINUE
  100       CONTINUE
*
*           Time CGEMM for a problem that corresponds to the following
*           problem for CHER2K:
*           CHER2K( UPLO( OP2 ), TRNS( OP3 ), DIM1( D ), DIM2( D ),
*                                       ALPHA, A, LDA( L ), B, LDA( L ), 
*                                            REAL( BETA ), C, LDA( L ) )
*
            IF( LSAME( TRNS( OP3 ), 'N' ) )THEN
               TM10 = SECOND( )
               CALL CGEMM( 'N', 'C', DIM1( D ), DIM1( D ), DIM2( D ),
     $              ALPHA, A, LDA( L ), B, LDA( L ), BETA, C, LDA( L ) )
               TM11 = SECOND( )
            ELSE IF( LSAME( TRNS( OP3 ), 'C' ) )THEN
               TM10 = SECOND( )
               CALL CGEMM( 'C', 'N', DIM1( D ), DIM1( D ), DIM2( D ),
     $              ALPHA, A, LDA( L ), B, LDA( L ), BETA, C, LDA( L ) )
               TM11 = SECOND( )
            ELSE
               WRITE( NERR, FMT = * )
     $              'Error: Unknown value for TRANS: ', TRNS( OP3 ), '.'
               STOP
            END IF
         END IF
         IF( TABSUB( 7 ).OR.TABSUB( 8 ) )THEN
            DO 120, J = 1, NMAX
               DO 110, I = 1, LD
                  C( I, J ) = C11 + CMPLX( 0.01E+0*
     $                      REAL( I+( J-1 )*NMAX )/ REAL( NMAX*NMAX+1 ),
     $                                   0.06E+0*REAL( I+( J-1 )*NMAX )/
     $                                             REAL( NMAX*NMAX+1 ) )
  110          CONTINUE
  120       CONTINUE
*
*           Time CGEMM for a problem that corresponds to the following
*           problems for CTRMM and CTRSM:
*           CTRMM( SIDE( OP1 ), UPLO( OP2 ), TRNS( OP3 ),
*                              DIAG( OP4 ), DIM1( D ), DIM2( D ), ALPHA,
*                                             A, LDA( L ), C, LDA( L ) )
*           CTRSM( SIDE( OP1 ), UPLO( OP2 ), TRNS( OP3 ),
*                              DIAG( OP4 ), DIM1( D ), DIM2( D ), ALPHA,
*                                             A, LDA( L ), C, LDA( L ) )
*
            IF( LSAME( SIDE( OP1 ), 'L' ) )THEN
*
*              C := alpha*A*C + C or C := alpha*A'*C + C. Use K = M.
*
               TM12 = SECOND( )
               CALL CGEMM( TRNS( OP3 ), 'N', DIM1( D ), DIM2( D ),
     $                       DIM1( D ), ALPHA, A, LDA( L ), B, LDA( L ),
     $                                                ONE, C, LDA( L ) )
               TM13 = SECOND( )
            ELSE IF( LSAME( SIDE( OP1 ), 'R' ) )THEN
*
*              C := alpha*C*A + C or C := alpha*C*A' + C. Use K = N.
*
               TM12 = SECOND( )
               CALL CGEMM( 'N', TRNS( OP3 ), DIM1( D ), DIM2( D ),
     $                       DIM2( D ), ALPHA, B, LDA( L ), A, LDA( L ),
     $                                                ONE, C, LDA( L ) )
               TM13 = SECOND( )
            ELSE
               WRITE( NERR, FMT = * )
     $               'Error: Unknown value for SIDE: ', SIDE( OP1 ), '.'
               STOP
            END IF
         END IF
*
*        Compute the performance of CGEMM in Mflops for problem
*        configurations that corresponds to CSYMM or CHEMM.
*
         IF( ( TABSUB( 1 ).OR.TABSUB( 2 ) ).AND.OP3.EQ.1 )THEN
            TIME = ( TM3 - TM2 ) - ( TM1 - TM0 )
            M = DIM1( D )
            N = DIM2( D )
            IF( LSAME( SIDE( OP1 ), 'L' ) )THEN
               MULTS = ( M + 1 )*M*N + MIN( M*N, M*M )
               ADDS = M*M*N
               NOPS = 6*MULTS + 2*ADDS
            ELSE IF( LSAME( SIDE( OP1 ), 'R' ) )THEN
               MULTS = ( N + 1 )*M*N + MIN( M*N, N*N )
               ADDS = M*N*N
               NOPS = 6*MULTS + 2*ADDS
            END IF
            IF( TIME.LE.ZERO )THEN
               SPEED = ZERO
            ELSE
               SPEED = REAL( NOPS )/( TIME*SCALE )
            END IF
            IF( RES( 1, OP1, 1, OP3, 1, D, L ).LT.SPEED )THEN
               DO 130, OP2 = 1, NUPLO
                  RES( 1, OP1, OP2, OP3, 1, D, L ) = SPEED
                  RES( 2, OP1, OP2, OP3, 1, D, L ) = SPEED
  130          CONTINUE
            END IF
         END IF
*
*        Compute the performance of CGEMM in Mflops for problem
*        configurations that corresponds to CSYRK.
*
         IF( TABSUB( 3 ).AND.OP1.EQ.1 )THEN
            TIME = ( TM5 - TM4 ) - ( TM1 - TM0 )
            N = DIM1( D )
            K = DIM2( D )
            MULTS = ( K + 1 )*N*N + MIN( N*K, N*N )
            ADDS = K*N*N
            NOPS = 6*MULTS + 2*ADDS
            IF( TIME.LE.ZERO )THEN
               SPEED = ZERO
            ELSE
               SPEED = REAL( NOPS )/( TIME*SCALE )
            END IF
            IF( RES( 3, OP1, 1, OP3, 1, D, L ).LT.SPEED )THEN
               DO 140, OP2 = 1, NUPLO
                  RES( 3, OP1, OP2, OP3, 1, D, L ) = SPEED
  140          CONTINUE
            END IF
         END IF
*
*        Compute the performance of CGEMM in Mflops for problem
*        configurations that corresponds to CHERK.
*
         IF( TABSUB( 4 ).AND.OP1.EQ.1 )THEN
            TIME = ( TM7 - TM6 ) - ( TM1 - TM0 )
            N = DIM1( D )
            K = DIM2( D )
            MULTS = ( K + 1 )*N*N + MIN( N*K, N*N )
            ADDS = K*N*N
            NOPS = 6*MULTS + 2*ADDS
            IF( TIME.LE.ZERO )THEN
               SPEED = ZERO
            ELSE
               SPEED = REAL( NOPS )/( TIME*SCALE )
            END IF
            IF( RES( 4, OP1, 1, OP3, 1, D, L ).LT.SPEED )THEN
               DO 150, OP2 = 1, NUPLO
                  RES( 4, OP1, OP2, OP3, 1, D, L ) = SPEED
  150          CONTINUE
            END IF
         END IF
*
*        Compute the performance of CGEMM in Mflops for problem
*        configurations that corresponds to CSYR2K.
*
         IF( TABSUB( 5 ).AND.OP1.EQ.1 )THEN
            TIME = ( TM9 - TM8 ) - ( TM1 - TM0 )
            N = DIM1( D )
            K = DIM2( D )
            NOPS = ( K + 1 )*N*N + MIN( N*K, N*N )
            NOPS = K*N*N
            NOPS = 6*MULTS + 2*ADDS
            IF( TIME.LE.ZERO )THEN
               SPEED = ZERO
            ELSE
               SPEED = REAL( NOPS )/( TIME*SCALE )
            END IF
            IF( RES( 5, OP1, 1, OP3, 1, D, L ).LT.SPEED )THEN
               DO 160, OP2 = 1, NUPLO
                  RES( 5, OP1, OP2, OP3, 1, D, L ) = SPEED
  160          CONTINUE
            END IF
         END IF
*
*        Compute the performance of CGEMM in Mflops for problem
*        configurations that corresponds to CHER2K.
*
         IF( TABSUB( 6 ).AND.OP1.EQ.1 )THEN
            TIME = ( TM11 - TM10 ) - ( TM1 - TM0 )
            N = DIM1( D )
            K = DIM2( D )
            NOPS = ( K + 1 )*N*N + MIN( N*K, N*N )
            NOPS = K*N*N
            NOPS = 6*MULTS + 2*ADDS
            IF( TIME.LE.ZERO )THEN
               SPEED = ZERO
            ELSE
               SPEED = REAL( NOPS )/( TIME*SCALE )
            END IF
            IF( RES( 6, OP1, 1, OP3, 1, D, L ).LT.SPEED )THEN
               DO 170, OP2 = 1, NUPLO
                  RES( 6, OP1, OP2, OP3, 1, D, L ) = SPEED
  170          CONTINUE
            END IF
         END IF
*
*        Compute the performance of CGEMM in Mflops for problem
*        configurations that corresponds to CTRMM and CTRSM.
*
         IF( TABSUB( 7 ).OR.TABSUB( 8 ) )THEN
            TIME = ( TM13 - TM12 ) - ( TM1 - TM0 )
            M = DIM1( D )
            N = DIM2( D )
            IF( LSAME( SIDE( OP1 ), 'L' ) )THEN
               MULTS = M*M*N + MIN( M*N, M*M )
               ADDS = ( M - 1 )*M*N
               NOPS = 6*MULTS + 2*ADDS
            ELSE IF( LSAME( SIDE( OP1 ), 'R' ) )THEN
               MULTS = M*N*N + MIN( M*N, N*N )
               ADDS = ( N - 1 )*M*N
               NOPS = 6*MULTS + 2*ADDS
            END IF
            IF( TIME.LE.ZERO )THEN
               SPEED = ZERO
            ELSE
               SPEED = REAL( NOPS )/( TIME*SCALE )
            END IF
            IF( RES( 7, OP1, 1, OP3, 1, D, L ).LT.SPEED )THEN
               DO 190, OP2 = 1, NUPLO
                  DO 180, OP4 = 1, NDIAG
                     RES( 7, OP1, OP2, OP3, OP4, D, L ) = SPEED
                     RES( 8, OP1, OP2, OP3, OP4, D, L ) = SPEED
  180             CONTINUE
  190          CONTINUE
            END IF
         END IF
  200 CONTINUE
*
*     ------ Stop indentation ------
*
  210 CONTINUE
  220 CONTINUE
  230 CONTINUE
  240 CONTINUE
*
*     ------ Continue indentation ------
*
      RETURN
*
*     End of SGBT02.
*
      END
