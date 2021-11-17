      SUBROUTINE DGBTP2( TAB, TABSUB, SIDE, NSIDE, UPLO, NUPLO, TRNS,
     $                  NTRNS, DIAG, NDIAG, DIM1, DIM2, NDIM, LDA, NLDA,
     $                    NOUT, MXTAB, MXSUB, MXOPT, MXDIM, MXLDA, RUNS,
     $                           ALPHA, BETA, LBL, USRES, MMRES, GBRES )
*     .. Scalar Arguments ..
      INTEGER            NOUT,
     $                   NSIDE, NUPLO, NTRNS, NDIAG, NDIM, NLDA,
     $                   MXTAB, MXSUB, MXOPT, MXDIM, MXLDA, RUNS
      DOUBLE PRECISION   ALPHA, BETA
*     .. Parameters ..
      INTEGER            LST
      PARAMETER        ( LST = 50 )
*     .. Array Arguments ..
      LOGICAL            TABSUB( MXSUB ), TAB( MXTAB )
      CHARACTER          LBL*( LST )
      CHARACTER          SIDE( MXOPT ), UPLO( MXOPT ), TRNS( MXOPT ),
     $                   DIAG( MXOPT )
      INTEGER            DIM1( MXDIM ), DIM2( MXDIM ), LDA( MXLDA )
      DOUBLE PRECISION   USRES( MXSUB, MXOPT, MXOPT, MXOPT, MXOPT,
     $                   MXDIM, MXLDA ),
     $                   GBRES( MXSUB, MXOPT, MXOPT, MXOPT, MXOPT,
     $                   MXDIM, MXLDA ),
     $                   MMRES( MXSUB, MXOPT, MXOPT, MXOPT, MXOPT,
     $                   MXDIM, MXLDA )
*
*
*  DGBTP2 prints tables showing detailed performance results and
*  comparisons between the user-supplied and the built-in GEMM-Based
*  Level 3 BLAS routines. The table results are intended for program
*  developers and others who are interested in detailed performance
*  presentations. Performance of the user-supplied and the built-in
*  GEMM-Based Level 3 BLAS routines are shown. The tables also show
*  GEMM-Efficiency and GEMM-Ratio. See the README and INSTALL files
*  for further information.
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
      INTEGER            D, I, L, NTIM, OP1, OP2, OP3, OP4
      DOUBLE PRECISION   MM, GE, GB, GR, US
*     .. Parameters ..
      INTEGER            MXTOTS, LLN
      DOUBLE PRECISION   ZERO, HUGE
      PARAMETER        ( MXTOTS = 6, LLN = 256, ZERO = 0.0D+0,
     $                   HUGE = 1.0D+10 )
      INTEGER            B1, B2, B3, B4, B5, B6, B7, B8, B9, B10, B11,
     $                   E1, E2, E3, E4, E5, E6, E7, E8, E9, E10, E11
      PARAMETER        ( B1 =  1, B2 =  3, B3 =  5, B4 = 7, B5 = 9,
     $                   B6 = 16, B7 = 23, B8 = 34, B9 = 45, B10 = 56,
     $                   B11 = 66,
     $                   E1 = 2, E2 = 4, E3 = 6, E4 = 8, E5 = 15,
     $                   E6 = 22, E7 = 33, E8 = 44, E9 = 55, E10 = 65,
     $                   E11 = 74 )
*     .. Local Arrays ..
      CHARACTER          OUTLN*( LLN ), OUTLN2*( LLN ), OUTLN3*( LLN )
      DOUBLE PRECISION   MI( MXTOTS ), MA( MXTOTS ), SU( MXTOTS )
*     ..
*     .. Executable Statements ..
*
*     Print an introduction.
*
      WRITE( NOUT, FMT = * )
      WRITE( NOUT, FMT = 9000 )
      WRITE( NOUT, FMT = 9010 )
      WRITE( NOUT, FMT = * )
      WRITE( NOUT, FMT = * )
      WRITE( NOUT, FMT = 9020 )
      WRITE( NOUT, FMT = 9030 ) 'SIDE   ', ( SIDE( I ), I = 1, NSIDE )
      WRITE( NOUT, FMT = 9030 ) 'UPLO   ', ( UPLO( I ), I = 1, NUPLO )
      WRITE( NOUT, FMT = 9030 ) 'TRANS  ', ( TRNS( I ), I = 1, NTRNS )
      WRITE( NOUT, FMT = 9030 ) 'DIAG   ', ( DIAG( I ), I = 1, NDIAG )
      WRITE( NOUT, FMT = 9040 ) 'DIM1   ', ( DIM1( I ), I = 1, NDIM )
      WRITE( NOUT, FMT = 9040 ) 'DIM2   ', ( DIM2( I ), I = 1, NDIM )
      WRITE( NOUT, FMT = 9040 ) 'LDA    ', ( LDA( I ), I = 1, NLDA )
      WRITE( NOUT, FMT = 9050 ) 'ALPHA  ', ALPHA
      WRITE( NOUT, FMT = 9050 ) 'BETA   ', BETA
      WRITE( NOUT, FMT = * )
      WRITE( NOUT, FMT = 9060 ) RUNS
      WRITE( NOUT, FMT = * )
      WRITE( NOUT, FMT = * )
      WRITE( NOUT, FMT = 9070 )
      WRITE( NOUT, FMT = * )
      WRITE( NOUT, FMT = * )
      WRITE( NOUT, FMT = 9080 )
      WRITE( NOUT, FMT = * )
      WRITE( NOUT, FMT = * )
      WRITE( NOUT, FMT = 9090 ) LBL
      WRITE( NOUT, FMT = * )
      WRITE( NOUT, FMT = * )
*
*     Print result tables for DSYMM.
*
      IF( TABSUB( 1 ) )THEN
         WRITE( NOUT, FMT = * )
         WRITE( NOUT, FMT = * )
         WRITE( NOUT, FMT = * )
         WRITE( NOUT, FMT = 9100 ) 'DSYMM ',
     $                                 '           OPTIONS  = SIDE,UPLO'
         WRITE( NOUT, FMT = * )
         WRITE( NOUT, FMT = 9110 ) 'M', 'N'
         DO 50, L = 1, NLDA
            WRITE( NOUT, FMT = * )
            WRITE( NOUT, FMT = 9120 ) LDA( L )
            DO 10, I = 1, MXTOTS
               MI( I ) = HUGE
               MA( I ) = ZERO
               SU( I ) = ZERO
   10       CONTINUE
            NTIM = 0
            DO 40, OP1 = 1, NSIDE
               WRITE( OUTLN( B1:E1 ), FMT = 9130 ) SIDE( OP1 )
               DO 30, OP2 = 1, NUPLO
                  WRITE( OUTLN( B2:E2 ), FMT = 9130 ) UPLO( OP2 )
                  WRITE( OUTLN( B3:E3 ), FMT = 9130 ) '  '
                  WRITE( OUTLN( B4:E4 ), FMT = 9130 ) '  '
                  DO 20, D = 1, NDIM
                     WRITE( OUTLN( B5:E5 ), FMT = 9140 ) DIM1( D )
                     WRITE( OUTLN( B6:E6 ), FMT = 9140 ) DIM2( D )
                     US = USRES( 1, OP1, OP2, 1, 1, D, L )
                     MM = MMRES( 1, OP1, OP2, 1, 1, D, L )
                     GB = GBRES( 1, OP1, OP2, 1, 1, D, L )
                     NTIM = NTIM + 1
                     IF( TAB( 2 ) )THEN
                        WRITE( OUTLN( B7:E7 ), FMT = 9150 ) GB
                        IF( MI( 2 ).GT.GB ) MI( 2 ) = GB
                        IF( MA( 2 ).LT.GB ) MA( 2 ) = GB
                        SU( 2 ) = SU( 2 ) + GB
                     ELSE
                        WRITE( OUTLN( B7:E7 ), FMT = 9160 )
                     END IF
                     IF( TAB( 3 ) )THEN
                        WRITE( OUTLN( B8:E8 ), FMT = 9150 ) US
                        IF( MI( 3 ).GT.US ) MI( 3 ) = US
                        IF( MA( 3 ).LT.US ) MA( 3 ) = US
                        SU( 3 ) = SU( 3 ) + US
                     ELSE
                        WRITE( OUTLN( B7:E8 ), FMT = 9160 )
                     END IF
                     IF( TAB( 4 ) )THEN
                        WRITE( OUTLN( B9:E9 ), FMT = 9150 ) MM
                        IF( MI( 4 ).GT.MM ) MI( 4 ) = MM
                        IF( MA( 4 ).LT.MM ) MA( 4 ) = MM
                        SU( 4 ) = SU( 4 ) + MM
                     ELSE
                        WRITE( OUTLN( B9:E9 ), FMT = 9160 )
                     END IF
                     IF( TAB( 5 ) )THEN
                        IF( MM.GT.ZERO )THEN
                           GE = US/MM
                        ELSE
                           GE = ZERO
                        END IF
                        WRITE( OUTLN( B10:E10 ), FMT = 9170 ) GE
                        IF( MI( 5 ).GT.GE ) MI( 5 ) = GE
                        IF( MA( 5 ).LT.GE ) MA( 5 ) = GE
                        SU( 5 ) = SU( 5 ) + GE
                     ELSE
                        WRITE( OUTLN( B10:E10 ), FMT = 9180 )
                     END IF
                     IF( TAB( 6 ) )THEN
                        IF( US.GT.ZERO )THEN
                           GR = GB/US
                        ELSE
                           GR = ZERO
                        END IF
                        WRITE( OUTLN( B11:E11 ), FMT = 9190 ) GR
                        IF( MI( 6 ).GT.GR ) MI( 6 ) = GR
                        IF( MA( 6 ).LT.GR ) MA( 6 ) = GR
                        SU( 6 ) = SU( 6 ) + GR
                     ELSE
                        WRITE( OUTLN( B11:E11 ), FMT = 9200 )
                     END IF
                     WRITE( NOUT, FMT = 9210 ) OUTLN( B1:E11 )
   20             CONTINUE
   30          CONTINUE
   40       CONTINUE
            WRITE( NOUT, FMT = 9220 )
*
*           Print the min, max, and mean values.
*
            WRITE( OUTLN( B1:E6 ), FMT = 9230 )
            WRITE( OUTLN2( B1:E6 ), FMT = 9240 )
            WRITE( OUTLN3( B1:E6 ), FMT = 9250 )
            IF( TAB( 2 ) )THEN
               WRITE( OUTLN( B7:E7 ), FMT = 9150 ) MI( 2 )
               WRITE( OUTLN2( B7:E7 ), FMT = 9150 ) MA( 2 )
               IF( NTIM.NE.0 )THEN
                  WRITE( OUTLN3( B7:E7 ), FMT = 9150 ) SU( 2 )/NTIM
               END IF
            ELSE
               WRITE( OUTLN( B7:E7 ), FMT = 9160 )
               WRITE( OUTLN2( B7:E7 ), FMT = 9160 )
               WRITE( OUTLN3( B7:E7 ), FMT = 9160 )
            END IF
            IF( TAB( 3 ) )THEN
               WRITE( OUTLN( B8:E8 ), FMT = 9150 ) MI( 3 )
               WRITE( OUTLN2( B8:E8 ), FMT = 9150 ) MA( 3 )
               IF( NTIM.NE.0 )THEN
                  WRITE( OUTLN3( B8:E8 ), FMT = 9150 ) SU( 3 )/NTIM
               END IF
            ELSE
               WRITE( OUTLN( B8:E8 ), FMT = 9160 )
               WRITE( OUTLN2( B8:E8 ), FMT = 9160 )
               WRITE( OUTLN3( B8:E8 ), FMT = 9160 )
            END IF
            IF( TAB( 4 ) )THEN
               WRITE( OUTLN( B9:E9 ), FMT = 9150 ) MI( 4 )
               WRITE( OUTLN2( B9:E9 ), FMT = 9150 ) MA( 4 )
               IF( NTIM.NE.0 )THEN
                  WRITE( OUTLN3( B9:E9 ), FMT = 9150 ) SU( 4 )/NTIM
               END IF
            ELSE
               WRITE( OUTLN( B9:E9 ), FMT = 9160 )
               WRITE( OUTLN2( B9:E9 ), FMT = 9160 )
               WRITE( OUTLN3( B9:E9 ), FMT = 9160 )
            END IF
            IF( TAB( 5 ) )THEN
               WRITE( OUTLN( B10:E10 ), FMT = 9170 ) MI( 5 )
               WRITE( OUTLN2( B10:E10 ), FMT = 9170 ) MA( 5 )
               IF( NTIM.NE.0 )THEN
                  WRITE( OUTLN3( B10:E10 ), FMT = 9170 ) SU( 5 )/NTIM
               END IF
            ELSE
               WRITE( OUTLN( B10:E10 ), FMT = 9180 )
               WRITE( OUTLN2( B10:E10 ), FMT = 9180 )
               WRITE( OUTLN3( B10:E10 ), FMT = 9180 )
            END IF
            IF( TAB( 6 ) )THEN
               WRITE( OUTLN( B11:E11 ), FMT = 9190 ) MI( 6 )
               WRITE( OUTLN2( B11:E11 ), FMT = 9190 ) MA( 6 )
               IF( NTIM.NE.0 )THEN
                  WRITE( OUTLN3( B11:E11 ), FMT = 9190 ) SU( 6 )/NTIM
               END IF
            ELSE
               WRITE( OUTLN( B11:E11 ), FMT = 9200 )
               WRITE( OUTLN2( B11:E11 ), FMT = 9200 )
               WRITE( OUTLN3( B11:E11 ), FMT = 9200 )
            END IF
            WRITE( NOUT, FMT = 9210 ) OUTLN( B1:E11 )
            WRITE( NOUT, FMT = 9210 ) OUTLN2( B1:E11 )
            IF( NTIM.NE.0 )THEN
               WRITE( NOUT, FMT = 9210 ) OUTLN3( B1:E11 )
            END IF
            WRITE( NOUT, FMT = * )
   50    CONTINUE
         WRITE( NOUT, FMT = * )
      END IF
*
*     Print result tables for DSYRK.
*
      IF( TABSUB( 2 ) )THEN
         WRITE( NOUT, FMT = * )
         WRITE( NOUT, FMT = * )
         WRITE( NOUT, FMT = * )
         WRITE( NOUT, FMT = 9100 ) 'DSYRK ',
     $                                 '          OPTIONS  = UPLO,TRANS'
         WRITE( NOUT, FMT = * )
         WRITE( NOUT, FMT = 9110 ) 'N', 'K'
         DO 100, L = 1, NLDA
            WRITE( NOUT, FMT = * )
            WRITE( NOUT, FMT = 9120 ) LDA( L )
            DO 60, I = 1, MXTOTS
               MI( I ) = HUGE
               MA( I ) = ZERO
               SU( I ) = ZERO
   60       CONTINUE
            NTIM = 0
            DO 90, OP2 = 1, NUPLO
               WRITE( OUTLN( B1:E1 ), FMT = 9130 ) UPLO( OP2 )
               DO 80, OP3 = 1, NTRNS
                  WRITE( OUTLN( B2:E2 ), FMT = 9130 ) TRNS( OP3 )
                  WRITE( OUTLN( B3:E3 ), FMT = 9130 ) '  '
                  WRITE( OUTLN( B4:E4 ), FMT = 9130 ) '  '
                  DO 70, D = 1, NDIM
                     WRITE( OUTLN( B5:E5 ), FMT = 9140 ) DIM1( D )
                     WRITE( OUTLN( B6:E6 ), FMT = 9140 ) DIM2( D )
                     US = USRES( 2, 1, OP2, OP3, 1, D, L )
                     MM = MMRES( 2, 1, OP2, OP3, 1, D, L )
                     GB = GBRES( 2, 1, OP2, OP3, 1, D, L )
                     NTIM = NTIM + 1
                     IF( TAB( 2 ) )THEN
                        WRITE( OUTLN( B7:E7 ), FMT = 9150 ) GB
                        IF( MI( 2 ).GT.GB ) MI( 2 ) = GB
                        IF( MA( 2 ).LT.GB ) MA( 2 ) = GB
                        SU( 2 ) = SU( 2 ) + GB
                     ELSE
                        WRITE( OUTLN( B7:E7 ), FMT = 9160 )
                     END IF
                     IF( TAB( 3 ) )THEN
                        WRITE( OUTLN( B8:E8 ), FMT = 9150 ) US
                        IF( MI( 3 ).GT.US ) MI( 3 ) = US
                        IF( MA( 3 ).LT.US ) MA( 3 ) = US
                        SU( 3 ) = SU( 3 ) + US
                     ELSE
                        WRITE( OUTLN( B8:E8 ), FMT = 9160 )
                     END IF
                     IF( TAB( 4 ) )THEN
                        WRITE( OUTLN( B9:E9 ), FMT = 9150 ) MM
                        IF( MI( 4 ).GT.MM ) MI( 4 ) = MM
                        IF( MA( 4 ).LT.MM ) MA( 4 ) = MM
                        SU( 4 ) = SU( 4 ) + MM
                     ELSE
                        WRITE( OUTLN( B9:E9 ), FMT = 9160 )
                     END IF
                     IF( TAB( 5 ) )THEN
                        IF( MM.GT.ZERO )THEN
                           GE = US/MM
                        ELSE
                           GE = ZERO
                        END IF
                        WRITE( OUTLN( B10:E10 ), FMT = 9170 ) GE
                        IF( MI( 5 ).GT.GE ) MI( 5 ) = GE
                        IF( MA( 5 ).LT.GE ) MA( 5 ) = GE
                        SU( 5 ) = SU( 5 ) + GE
                     ELSE
                        WRITE( OUTLN( B10:E10 ), FMT = 9180 )
                     END IF
                     IF( TAB( 6 ) )THEN
                        IF( US.GT.ZERO )THEN
                           GR = GB/US
                        ELSE
                           GR = ZERO
                        END IF
                        WRITE( OUTLN( B11:E11 ), FMT = 9190 ) GR
                        IF( MI( 6 ).GT.GR ) MI( 6 ) = GR
                        IF( MA( 6 ).LT.GR ) MA( 6 ) = GR
                        SU( 6 ) = SU( 6 ) + GR
                     ELSE
                        WRITE( OUTLN( B11:E11 ), FMT = 9200 )
                     END IF
                     WRITE( NOUT, FMT = 9210 ) OUTLN( B1:E11 )
   70             CONTINUE
   80          CONTINUE
   90       CONTINUE
            WRITE( NOUT, FMT = 9220 )
*
*           Print the min, max, and mean values.
*
            WRITE( OUTLN( B1:E6 ), FMT = 9230 )
            WRITE( OUTLN2( B1:E6 ), FMT = 9240 )
            WRITE( OUTLN3( B1:E6 ), FMT = 9250 )
            IF( TAB( 2 ) )THEN
               WRITE( OUTLN( B7:E7 ), FMT = 9150 ) MI( 2 )
               WRITE( OUTLN2( B7:E7 ), FMT = 9150 ) MA( 2 )
               IF( NTIM.NE.0 )THEN
                  WRITE( OUTLN3( B7:E7 ), FMT = 9150 ) SU( 2 )/NTIM
               END IF
            ELSE
               WRITE( OUTLN( B7:E7 ), FMT = 9160 )
               WRITE( OUTLN2( B7:E7 ), FMT = 9160 )
               WRITE( OUTLN3( B7:E7 ), FMT = 9160 )
            END IF
            IF( TAB( 3 ) )THEN
               WRITE( OUTLN( B8:E8 ), FMT = 9150 ) MI( 3 )
               WRITE( OUTLN2( B8:E8 ), FMT = 9150 ) MA( 3 )
               IF( NTIM.NE.0 )THEN
                  WRITE( OUTLN3( B8:E8 ), FMT = 9150 ) SU( 3 )/NTIM
               END IF
            ELSE
               WRITE( OUTLN( B8:E8 ), FMT = 9160 )
               WRITE( OUTLN2( B8:E8 ), FMT = 9160 )
               WRITE( OUTLN3( B8:E8 ), FMT = 9160 )
            END IF
            IF( TAB( 4 ) )THEN
               WRITE( OUTLN( B9:E9 ), FMT = 9150 ) MI( 4 )
               WRITE( OUTLN2( B9:E9 ), FMT = 9150 ) MA( 4 )
               IF( NTIM.NE.0 )THEN
                  WRITE( OUTLN3( B9:E9 ), FMT = 9150 ) SU( 4 )/NTIM
               END IF
            ELSE
               WRITE( OUTLN( B9:E9 ), FMT = 9160 )
               WRITE( OUTLN2( B9:E9 ), FMT = 9160 )
               WRITE( OUTLN3( B9:E9 ), FMT = 9160 )
            END IF
            IF( TAB( 5 ) )THEN
               WRITE( OUTLN( B10:E10 ), FMT = 9170 ) MI( 5 )
               WRITE( OUTLN2( B10:E10 ), FMT = 9170 ) MA( 5 )
               IF( NTIM.NE.0 )THEN
                  WRITE( OUTLN3( B10:E10 ), FMT = 9170 ) SU( 5 )/NTIM
               END IF
            ELSE
               WRITE( OUTLN( B10:E10 ), FMT = 9180 )
               WRITE( OUTLN2( B10:E10 ), FMT = 9180 )
               WRITE( OUTLN3( B10:E10 ), FMT = 9180 )
            END IF
            IF( TAB( 6 ) )THEN
               WRITE( OUTLN( B11:E11 ), FMT = 9190 ) MI( 6 )
               WRITE( OUTLN2( B11:E11 ), FMT = 9190 ) MA( 6 )
               IF( NTIM.NE.0 )THEN
                  WRITE( OUTLN3( B11:E11 ), FMT = 9190 ) SU( 6 )/NTIM
               END IF
            ELSE
               WRITE( OUTLN( B11:E11 ), FMT = 9200 )
               WRITE( OUTLN2( B11:E11 ), FMT = 9200 )
               WRITE( OUTLN3( B11:E11 ), FMT = 9200 )
            END IF
            WRITE( NOUT, FMT = 9210 ) OUTLN( B1:E11 )
            WRITE( NOUT, FMT = 9210 ) OUTLN2( B1:E11 )
            IF( NTIM.NE.0 )THEN
               WRITE( NOUT, FMT = 9210 ) OUTLN3( B1:E11 )
            END IF
            WRITE( NOUT, FMT = * )
  100    CONTINUE
         WRITE( NOUT, FMT = * )
      END IF
*
*     Print result tables for DSYR2K.
*
      IF( TABSUB( 3 ) )THEN
         WRITE( NOUT, FMT = * )
         WRITE( NOUT, FMT = * )
         WRITE( NOUT, FMT = * )
         WRITE( NOUT, FMT = 9100 ) 'DSYR2K',
     $                                 '          OPTIONS  = UPLO,TRANS'
         WRITE( NOUT, FMT = * )
         WRITE( NOUT, FMT = 9110 ) 'N', 'K'
         DO 150, L = 1, NLDA
            WRITE( NOUT, FMT = * )
            WRITE( NOUT, FMT = 9120 ) LDA( L )
            DO 110, I = 1, MXTOTS
               MI( I ) = HUGE
               MA( I ) = ZERO
               SU( I ) = ZERO
  110       CONTINUE
            NTIM = 0
            DO 140, OP2 = 1, NUPLO
               WRITE( OUTLN( B1:E1 ), FMT = 9130 ) UPLO( OP2 )
               DO 130, OP3 = 1, NTRNS
                  WRITE( OUTLN( B2:E2 ), FMT = 9130 ) TRNS( OP3 )
                  WRITE( OUTLN( B3:E3 ), FMT = 9130 ) '  '
                  WRITE( OUTLN( B4:E4 ), FMT = 9130 ) '  '
                  DO 120, D = 1, NDIM
                     WRITE( OUTLN( B5:E5 ), FMT = 9140 ) DIM1( D )
                     WRITE( OUTLN( B6:E6 ), FMT = 9140 ) DIM2( D )
                     US = USRES( 3, 1, OP2, OP3, 1, D, L )
                     MM = MMRES( 3, 1, OP2, OP3, 1, D, L )
                     GB = GBRES( 3, 1, OP2, OP3, 1, D, L )
                     NTIM = NTIM + 1
                     IF( TAB( 2 ) )THEN
                        WRITE( OUTLN( B7:E7 ), FMT = 9150 ) GB
                        IF( MI( 2 ).GT.GB ) MI( 2 ) = GB
                        IF( MA( 2 ).LT.GB ) MA( 2 ) = GB
                        SU( 2 ) = SU( 2 ) + GB
                     ELSE
                        WRITE( OUTLN( B7:E7 ), FMT = 9160 )
                     END IF
                     IF( TAB( 3 ) )THEN
                        WRITE( OUTLN( B8:E8 ), FMT = 9150 ) US
                        IF( MI( 3 ).GT.US ) MI( 3 ) = US
                        IF( MA( 3 ).LT.US ) MA( 3 ) = US
                        SU( 3 ) = SU( 3 ) + US
                     ELSE
                        WRITE( OUTLN( B8:E8 ), FMT = 9160 )
                     END IF
                     IF( TAB( 4 ) )THEN
                        WRITE( OUTLN( B9:E9 ), FMT = 9150 ) MM
                        IF( MI( 4 ).GT.MM ) MI( 4 ) = MM
                        IF( MA( 4 ).LT.MM ) MA( 4 ) = MM
                        SU( 4 ) = SU( 4 ) + MM
                     ELSE
                        WRITE( OUTLN( B9:E9 ), FMT = 9160 )
                     END IF
                     IF( TAB( 5 ) )THEN
                        IF( MM.GT.ZERO )THEN
                           GE = US/MM
                        ELSE
                           GE = ZERO
                        END IF
                        WRITE( OUTLN( B10:E10 ), FMT = 9170 ) GE
                        IF( MI( 5 ).GT.GE ) MI( 5 ) = GE
                        IF( MA( 5 ).LT.GE ) MA( 5 ) = GE
                        SU( 5 ) = SU( 5 ) + GE
                     ELSE
                        WRITE( OUTLN( B10:E10 ), FMT = 9180 )
                     END IF
                     IF( TAB( 6 ) )THEN
                        IF( US.GT.ZERO )THEN
                           GR = GB/US
                        ELSE
                           GR = ZERO
                        END IF
                        WRITE( OUTLN( B11:E11 ), FMT = 9190 ) GR
                        IF( MI( 6 ).GT.GR ) MI( 6 ) = GR
                        IF( MA( 6 ).LT.GR ) MA( 6 ) = GR
                        SU( 6 ) = SU( 6 ) + GR
                     ELSE
                        WRITE( OUTLN( B11:E11 ), FMT = 9200 )
                     END IF
                     WRITE( NOUT, FMT = 9210 ) OUTLN( B1:E11 )
  120             CONTINUE
  130          CONTINUE
  140       CONTINUE
            WRITE( NOUT, FMT = 9220 )
*
*           Print the min, max, and mean values.
*
            WRITE( OUTLN( B1:E6 ), FMT = 9230 )
            WRITE( OUTLN2( B1:E6 ), FMT = 9240 )
            WRITE( OUTLN3( B1:E6 ), FMT = 9250 )
            IF( TAB( 2 ) )THEN
               WRITE( OUTLN( B7:E7 ), FMT = 9150 ) MI( 2 )
               WRITE( OUTLN2( B7:E7 ), FMT = 9150 ) MA( 2 )
               IF( NTIM.NE.0 )THEN
                  WRITE( OUTLN3( B7:E7 ), FMT = 9150 ) SU( 2 )/NTIM
               END IF
            ELSE
               WRITE( OUTLN( B7:E7 ), FMT = 9160 )
               WRITE( OUTLN2( B7:E7 ), FMT = 9160 )
               WRITE( OUTLN3( B7:E7 ), FMT = 9160 )
            END IF
            IF( TAB( 3 ) )THEN
               WRITE( OUTLN( B8:E8 ), FMT = 9150 ) MI( 3 )
               WRITE( OUTLN2( B8:E8 ), FMT = 9150 ) MA( 3 )
               IF( NTIM.NE.0 )THEN
                  WRITE( OUTLN3( B8:E8 ), FMT = 9150 ) SU( 3 )/NTIM
               END IF
            ELSE
               WRITE( OUTLN( B8:E8 ), FMT = 9160 )
               WRITE( OUTLN2( B8:E8 ), FMT = 9160 )
               WRITE( OUTLN3( B8:E8 ), FMT = 9160 )
            END IF
            IF( TAB( 4 ) )THEN
               WRITE( OUTLN( B9:E9 ), FMT = 9150 ) MI( 4 )
               WRITE( OUTLN2( B9:E9 ), FMT = 9150 ) MA( 4 )
               IF( NTIM.NE.0 )THEN
                  WRITE( OUTLN3( B9:E9 ), FMT = 9150 ) SU( 4 )/NTIM
               END IF
            ELSE
               WRITE( OUTLN( B9:E9 ), FMT = 9160 )
               WRITE( OUTLN2( B9:E9 ), FMT = 9160 )
               WRITE( OUTLN3( B9:E9 ), FMT = 9160 )
            END IF
            IF( TAB( 5 ) )THEN
               WRITE( OUTLN( B10:E10 ), FMT = 9170 ) MI( 5 )
               WRITE( OUTLN2( B10:E10 ), FMT = 9170 ) MA( 5 )
               IF( NTIM.NE.0 )THEN
                  WRITE( OUTLN3( B10:E10 ), FMT = 9170 ) SU( 5 )/NTIM
               END IF
            ELSE
               WRITE( OUTLN( B10:E10 ), FMT = 9180 )
               WRITE( OUTLN2( B10:E10 ), FMT = 9180 )
               WRITE( OUTLN3( B10:E10 ), FMT = 9180 )
            END IF
            IF( TAB( 6 ) )THEN
               WRITE( OUTLN( B11:E11 ), FMT = 9190 ) MI( 6 )
               WRITE( OUTLN2( B11:E11 ), FMT = 9190 ) MA( 6 )
               IF( NTIM.NE.0 )THEN
                  WRITE( OUTLN3( B11:E11 ), FMT = 9190 ) SU( 6 )/NTIM
               END IF
            ELSE
               WRITE( OUTLN( B11:E11 ), FMT = 9200 )
               WRITE( OUTLN2( B11:E11 ), FMT = 9200 )
               WRITE( OUTLN3( B11:E11 ), FMT = 9200 )
            END IF
            WRITE( NOUT, FMT = 9210 ) OUTLN( B1:E11 )
            WRITE( NOUT, FMT = 9210 ) OUTLN2( B1:E11 )
            IF( NTIM.NE.0 )THEN
               WRITE( NOUT, FMT = 9210 ) OUTLN3( B1:E11 )
            END IF
            WRITE( NOUT, FMT = * )
  150    CONTINUE
         WRITE( NOUT, FMT = * )
      END IF
*
*     Print result tables for DTRMM.
*
      IF( TABSUB( 4 ) )THEN
         WRITE( NOUT, FMT = * )
         WRITE( NOUT, FMT = * )
         WRITE( NOUT, FMT = * )
         WRITE( NOUT, FMT = 9100 ) 'DTRMM ',
     $                                 'OPTIONS  = SIDE,UPLO,TRANS,DIAG'
         WRITE( NOUT, FMT = * )
         WRITE( NOUT, FMT = 9110 ) 'M', 'N'
         DO 220, L = 1, NLDA
            WRITE( NOUT, FMT = * )
            WRITE( NOUT, FMT = 9120 ) LDA( L )
            DO 160, I = 1, MXTOTS
               MI( I ) = HUGE
               MA( I ) = ZERO
               SU( I ) = ZERO
  160       CONTINUE
            NTIM = 0
            DO 210, OP1 = 1, NSIDE
               WRITE( OUTLN( B1:E1 ), FMT = 9130 ) SIDE( OP1 )
               DO 200, OP2 = 1, NUPLO
                  WRITE( OUTLN( B2:E2 ), FMT = 9130 ) UPLO( OP2 )
                  DO 190, OP3 = 1, NTRNS
                     WRITE( OUTLN( B3:E3 ), FMT = 9130 ) TRNS( OP3 )
                     DO 180, OP4 = 1, NDIAG
                        WRITE( OUTLN( B4:E4 ), FMT = 9130 ) DIAG( OP4 )
                        DO 170, D = 1, NDIM
                           WRITE( OUTLN( B5:E5 ), FMT = 9140 ) DIM1( D )
                           WRITE( OUTLN( B6:E6 ), FMT = 9140 ) DIM2( D )
                           US = USRES( 4, OP1, OP2, OP3, OP4, D, L )
                           MM = MMRES( 4, OP1, OP2, OP3, OP4, D, L )
                           GB = GBRES( 4, OP1, OP2, OP3, OP4, D, L )
                           NTIM = NTIM + 1
                           IF( TAB( 2 ) )THEN
                              WRITE( OUTLN( B7:E7 ), FMT = 9150 ) GB
                              IF( MI( 2 ).GT.GB ) MI( 2 ) = GB
                              IF( MA( 2 ).LT.GB ) MA( 2 ) = GB
                              SU( 2 ) = SU( 2 ) + GB
                           ELSE
                              WRITE( OUTLN( B7:E7 ), FMT = 9160 )
                           END IF
                           IF( TAB( 3 ) )THEN
                              WRITE( OUTLN( B8:E8 ), FMT = 9150 ) US
                              IF( MI( 3 ).GT.US ) MI( 3 ) = US
                              IF( MA( 3 ).LT.US ) MA( 3 ) = US
                              SU( 3 ) = SU( 3 ) + US
                           ELSE
                              WRITE( OUTLN( B8:E8 ), FMT = 9160 )
                           END IF
                           IF( TAB( 4 ) )THEN
                              WRITE( OUTLN( B9:E9 ), FMT = 9150 ) MM
                              IF( MI( 4 ).GT.MM ) MI( 4 ) = MM
                              IF( MA( 4 ).LT.MM ) MA( 4 ) = MM
                              SU( 4 ) = SU( 4 ) + MM
                           ELSE
                              WRITE( OUTLN( B9:E9 ), FMT = 9160 )
                           END IF
                           IF( TAB( 5 ) )THEN
                              IF( MM.GT.ZERO )THEN
                                 GE = US/MM
                              ELSE
                                 GE = ZERO
                              END IF
                              WRITE( OUTLN( B10:E10 ), FMT = 9170 ) GE
                              IF( MI( 5 ).GT.GE ) MI( 5 ) = GE
                              IF( MA( 5 ).LT.GE ) MA( 5 ) = GE
                              SU( 5 ) = SU( 5 ) + GE
                           ELSE
                              WRITE( OUTLN( B10:E10 ), FMT = 9180 )
                           END IF
                           IF( TAB( 6 ) )THEN
                              IF( US.GT.ZERO )THEN
                                 GR = GB/US
                              ELSE
                                 GR = ZERO
                              END IF
                              WRITE( OUTLN( B11:E11 ), FMT = 9190 ) GR
                              IF( MI( 6 ).GT.GR ) MI( 6 ) = GR
                              IF( MA( 6 ).LT.GR ) MA( 6 ) = GR
                              SU( 6 ) = SU( 6 ) + GR
                           ELSE
                              WRITE( OUTLN( B11:E11 ), FMT = 9200 )
                           END IF
                           WRITE( NOUT, FMT = 9210 ) OUTLN( B1:E11 )
  170                   CONTINUE
  180                CONTINUE
  190             CONTINUE
  200          CONTINUE
  210       CONTINUE
            WRITE( NOUT, FMT = 9220 )
*
*           Print the min, max, and mean values.
*
            WRITE( OUTLN( B1:E6 ), FMT = 9230 )
            WRITE( OUTLN2( B1:E6 ), FMT = 9240 )
            WRITE( OUTLN3( B1:E6 ), FMT = 9250 )
            IF( TAB( 2 ) )THEN
               WRITE( OUTLN( B7:E7 ), FMT = 9150 ) MI( 2 )
               WRITE( OUTLN2( B7:E7 ), FMT = 9150 ) MA( 2 )
               IF( NTIM.NE.0 )THEN
                  WRITE( OUTLN3( B7:E7 ), FMT = 9150 ) SU( 2 )/NTIM
               END IF
            ELSE
               WRITE( OUTLN( B7:E7 ), FMT = 9160 )
               WRITE( OUTLN2( B7:E7 ), FMT = 9160 )
               WRITE( OUTLN3( B7:E7 ), FMT = 9160 )
            END IF
            IF( TAB( 3 ) )THEN
               WRITE( OUTLN( B8:E8 ), FMT = 9150 ) MI( 3 )
               WRITE( OUTLN2( B8:E8 ), FMT = 9150 ) MA( 3 )
               IF( NTIM.NE.0 )THEN
                  WRITE( OUTLN3( B8:E8 ), FMT = 9150 ) SU( 3 )/NTIM
               END IF
            ELSE
               WRITE( OUTLN( B8:E8 ), FMT = 9160 )
               WRITE( OUTLN2( B8:E8 ), FMT = 9160 )
               WRITE( OUTLN3( B8:E8 ), FMT = 9160 )
            END IF
            IF( TAB( 4 ) )THEN
               WRITE( OUTLN( B9:E9 ), FMT = 9150 ) MI( 4 )
               WRITE( OUTLN2( B9:E9 ), FMT = 9150 ) MA( 4 )
               IF( NTIM.NE.0 )THEN
                  WRITE( OUTLN3( B9:E9 ), FMT = 9150 ) SU( 4 )/NTIM
               END IF
            ELSE
               WRITE( OUTLN( B9:E9 ), FMT = 9160 )
               WRITE( OUTLN2( B9:E9 ), FMT = 9160 )
               WRITE( OUTLN3( B9:E9 ), FMT = 9160 )
            END IF
            IF( TAB( 5 ) )THEN
               WRITE( OUTLN( B10:E10 ), FMT = 9170 ) MI( 5 )
               WRITE( OUTLN2( B10:E10 ), FMT = 9170 ) MA( 5 )
               IF( NTIM.NE.0 )THEN
                  WRITE( OUTLN3( B10:E10 ), FMT = 9170 ) SU( 5 )/NTIM
               END IF
            ELSE
               WRITE( OUTLN( B10:E10 ), FMT = 9180 )
               WRITE( OUTLN2( B10:E10 ), FMT = 9180 )
               WRITE( OUTLN3( B10:E10 ), FMT = 9180 )
            END IF
            IF( TAB( 6 ) )THEN
               WRITE( OUTLN( B11:E11 ), FMT = 9190 ) MI( 6 )
               WRITE( OUTLN2( B11:E11 ), FMT = 9190 ) MA( 6 )
               IF( NTIM.NE.0 )THEN
                  WRITE( OUTLN3( B11:E11 ), FMT = 9190 ) SU( 6 )/NTIM
               END IF
            ELSE
               WRITE( OUTLN( B11:E11 ), FMT = 9200 )
               WRITE( OUTLN2( B11:E11 ), FMT = 9200 )
               WRITE( OUTLN3( B11:E11 ), FMT = 9200 )
            END IF
            WRITE( NOUT, FMT = 9210 ) OUTLN( B1:E11 )
            WRITE( NOUT, FMT = 9210 ) OUTLN2( B1:E11 )
            IF( NTIM.NE.0 )THEN
               WRITE( NOUT, FMT = 9210 ) OUTLN3( B1:E11 )
            END IF
            WRITE( NOUT, FMT = * )
  220    CONTINUE
         WRITE( NOUT, FMT = * )
      END IF
*
*     Print result tables for DTRSM.
*
      IF( TABSUB( 5 ) )THEN
         WRITE( NOUT, FMT = * )
         WRITE( NOUT, FMT = * )
         WRITE( NOUT, FMT = * )
         WRITE( NOUT, FMT = 9100 ) 'DTRSM ',
     $                                 'OPTIONS  = SIDE,UPLO,TRANS,DIAG'
         WRITE( NOUT, FMT = * )
         WRITE( NOUT, FMT = 9110 ) 'M', 'N'
         DO 290, L = 1, NLDA
            WRITE( NOUT, FMT = * )
            WRITE( NOUT, FMT = 9120 ) LDA( L )
            DO 230, I = 1, MXTOTS
               MI( I ) = HUGE
               MA( I ) = ZERO
               SU( I ) = ZERO
  230       CONTINUE
            NTIM = 0
            DO 280, OP1 = 1, NSIDE
               WRITE( OUTLN( B1:E1 ), FMT = 9130 ) SIDE( OP1 )
               DO 270, OP2 = 1, NUPLO
                  WRITE( OUTLN( B2:E2 ), FMT = 9130 ) UPLO( OP2 )
                  DO 260, OP3 = 1, NTRNS
                     WRITE( OUTLN( B3:E3 ), FMT = 9130 ) TRNS( OP3 )
                     DO 250, OP4 = 1, NDIAG
                        WRITE( OUTLN( B4:E4 ), FMT = 9130 ) DIAG( OP4 )
                        DO 240, D = 1, NDIM
                           WRITE( OUTLN( B5:E5 ), FMT = 9140 ) DIM1( D )
                           WRITE( OUTLN( B6:E6 ), FMT = 9140 ) DIM2( D )
                           US = USRES( 5, OP1, OP2, OP3, OP4, D, L )
                           MM = MMRES( 5, OP1, OP2, OP3, OP4, D, L )
                           GB = GBRES( 5, OP1, OP2, OP3, OP4, D, L )
                           NTIM = NTIM + 1
                           IF( TAB( 2 ) )THEN
                              WRITE( OUTLN( B7:E7 ), FMT = 9150 ) GB
                              IF( MI( 2 ).GT.GB ) MI( 2 ) = GB
                              IF( MA( 2 ).LT.GB ) MA( 2 ) = GB
                              SU( 2 ) = SU( 2 ) + GB
                           ELSE
                              WRITE( OUTLN( B7:E7 ), FMT = 9160 )
                           END IF
                           IF( TAB( 3 ) )THEN
                              WRITE( OUTLN( B8:E8 ), FMT = 9150 ) US
                              IF( MI( 3 ).GT.US ) MI( 3 ) = US
                              IF( MA( 3 ).LT.US ) MA( 3 ) = US
                              SU( 3 ) = SU( 3 ) + US
                           ELSE
                              WRITE( OUTLN( B8:E8 ), FMT = 9160 )
                           END IF
                           IF( TAB( 4 ) )THEN
                              WRITE( OUTLN( B9:E9 ), FMT = 9150 ) MM
                              IF( MI( 4 ).GT.MM ) MI( 4 ) = MM
                              IF( MA( 4 ).LT.MM ) MA( 4 ) = MM
                              SU( 4 ) = SU( 4 ) + MM
                           ELSE
                              WRITE( OUTLN( B9:E9 ), FMT = 9160 )
                           END IF
                           IF( TAB( 5 ) )THEN
                              IF( MM.GT.ZERO )THEN
                                 GE = US/MM
                              ELSE
                                 GE = ZERO
                              END IF
                              WRITE( OUTLN( B10:E10 ), FMT = 9170 ) GE
                              IF( MI( 5 ).GT.GE ) MI( 5 ) = GE
                              IF( MA( 5 ).LT.GE ) MA( 5 ) = GE
                              SU( 5 ) = SU( 5 ) + GE
                           ELSE
                              WRITE( OUTLN( B10:E10 ), FMT = 9190 )
                           END IF
                           IF( TAB( 6 ) )THEN
                              IF( US.GT.ZERO )THEN
                                 GR = GB/US
                              ELSE
                                 GR = ZERO
                              END IF
                              WRITE( OUTLN( B11:E11 ), FMT = 9190 ) GR
                              IF( MI( 6 ).GT.GR ) MI( 6 ) = GR
                              IF( MA( 6 ).LT.GR ) MA( 6 ) = GR
                              SU( 6 ) = SU( 6 ) + GR
                           ELSE
                              WRITE( OUTLN( B11:E11 ), FMT = 9200 )
                           END IF
                           WRITE( NOUT, FMT = 9210 ) OUTLN( B1:E11 )
  240                   CONTINUE
  250                CONTINUE
  260             CONTINUE
  270          CONTINUE
  280       CONTINUE
            WRITE( NOUT, FMT = 9220 )
*
*           Print the min, max, and mean values.
*
            WRITE( OUTLN( B1:E6 ), FMT = 9230 )
            WRITE( OUTLN2( B1:E6 ), FMT = 9240 )
            WRITE( OUTLN3( B1:E6 ), FMT = 9250 )
            IF( TAB( 2 ) )THEN
               WRITE( OUTLN( B7:E7 ), FMT = 9150 ) MI( 2 )
               WRITE( OUTLN2( B7:E7 ), FMT = 9150 ) MA( 2 )
               IF( NTIM.NE.0 )THEN
                  WRITE( OUTLN3( B7:E7 ), FMT = 9150 ) SU( 2 )/NTIM
               END IF
            ELSE
               WRITE( OUTLN( B7:E7 ), FMT = 9160 )
               WRITE( OUTLN2( B7:E7 ), FMT = 9160 )
               WRITE( OUTLN3( B7:E7 ), FMT = 9160 )
            END IF
            IF( TAB( 3 ) )THEN
               WRITE( OUTLN( B8:E8 ), FMT = 9150 ) MI( 3 )
               WRITE( OUTLN2( B8:E8 ), FMT = 9150 ) MA( 3 )
               IF( NTIM.NE.0 )THEN
                  WRITE( OUTLN3( B8:E8 ), FMT = 9150 ) SU( 3 )/NTIM
               END IF
            ELSE
               WRITE( OUTLN( B8:E8 ), FMT = 9160 )
               WRITE( OUTLN2( B8:E8 ), FMT = 9160 )
               WRITE( OUTLN3( B8:E8 ), FMT = 9160 )
            END IF
            IF( TAB( 4 ) )THEN
               WRITE( OUTLN( B9:E9 ), FMT = 9150 ) MI( 4 )
               WRITE( OUTLN2( B9:E9 ), FMT = 9150 ) MA( 4 )
               IF( NTIM.NE.0 )THEN
                  WRITE( OUTLN3( B9:E9 ), FMT = 9150 ) SU( 4 )/NTIM
               END IF
            ELSE
               WRITE( OUTLN( B9:E9 ), FMT = 9160 )
               WRITE( OUTLN2( B9:E9 ), FMT = 9160 )
               WRITE( OUTLN3( B9:E9 ), FMT = 9160 )
            END IF
            IF( TAB( 5 ) )THEN
               WRITE( OUTLN( B10:E10 ), FMT = 9170 ) MI( 5 )
               WRITE( OUTLN2( B10:E10 ), FMT = 9170 ) MA( 5 )
               IF( NTIM.NE.0 )THEN
                  WRITE( OUTLN3( B10:E10 ), FMT = 9170 ) SU( 5 )/NTIM
               END IF
            ELSE
               WRITE( OUTLN( B10:E10 ), FMT = 9180 )
               WRITE( OUTLN2( B10:E10 ), FMT = 9180 )
               WRITE( OUTLN3( B10:E10 ), FMT = 9180 )
            END IF
            IF( TAB( 6 ) )THEN
               WRITE( OUTLN( B11:E11 ), FMT = 9190 ) MI( 6 )
               WRITE( OUTLN2( B11:E11 ), FMT = 9190 ) MA( 6 )
               IF( NTIM.NE.0 )THEN
                  WRITE( OUTLN3( B11:E11 ), FMT = 9190 ) SU( 6 )/NTIM
               END IF
            ELSE
               WRITE( OUTLN( B11:E11 ), FMT = 9200 )
               WRITE( OUTLN2( B11:E11 ), FMT = 9200 )
               WRITE( OUTLN3( B11:E11 ), FMT = 9200 )
            END IF
            WRITE( NOUT, FMT = 9210 ) OUTLN( B1:E11 )
            WRITE( NOUT, FMT = 9210 ) OUTLN2( B1:E11 )
            IF( NTIM.NE.0 )THEN
               WRITE( NOUT, FMT = 9210 ) OUTLN3( B1:E11 )
            END IF
            WRITE( NOUT, FMT = * )
  290    CONTINUE
      END IF
      WRITE( NOUT, FMT = * )
      WRITE( NOUT, FMT = * )
      WRITE( NOUT, FMT = 9260 )
*
      RETURN
*
*     Print formats.
*
 9000 FORMAT( 17X, '****  GEMM-Based Level 3 BLAS Benchmark  ****' )
 9010 FORMAT( 33X, 'Table Results',/,
     $        32X, 'Double Precision' )
 9020 FORMAT(  8X, 'Input parameters.' )
 9030 FORMAT(  8X, A, 3X, 10( A, ' ' ) )
 9040 FORMAT(  8X, A, 1X, 12( I5 ), 2( /, 16X, 12( I5 ) ) )
 9050 FORMAT(  8X, A, F6.1 )
 9060 FORMAT(  8X, 'Results are based on the shortest execution time ',
     $             'of ', I2, ' runs for ',/,
     $         8X, 'each problem configuration.' )
 9070 FORMAT( 27X, 'Performance of a user-supplied',/,
     $        27X, 'Level 3 BLAS routine (megaflops).',/,
     $         8X, 'GEMM-Efficiency = -------------------------------',
     $             '----',/,
     $        27X, 'Performance of the user-supplied',/,
     $        27X, 'DGEMM routine (megaflops).' )
 9080 FORMAT( 22X, 'Performance for the internal GEMM-Based',/,
     $        22X, 'Level 3 BLAS routine Dxxxx (megaflops).',/,
     $         8X, 'GEMM-Ratio = ------------------------------------',
     $             '-----',/,
     $        22X, 'Performance of the user-supplied',/,
     $        22X, 'Level 3 BLAS routine Dxxxx (megaflops).' )
 9090 FORMAT(  8X, 'Test label:  ', A )
 9100 FORMAT(  2X, A, 38X, A )
 9110 FORMAT( 31X, 'GEMM-      User-', /,
     $        29X,'Based lib  suppl lib   DGEMM      GEMM-    GEMM-', /,
     $         2X, 'OPTIONS       ', A,'      ', A,'      ',
     $             'Mflops     Mflops     Mflops     Eff.     Ratio', /,
     $         2X, '==================================================',
     $             '=========================' )
 9120 FORMAT(  2X, '( LDA = ', I4, ' )' )
 9130 FORMAT(  A )
 9140 FORMAT(  I7 )
 9150 FORMAT(  F11.1 )
 9160 FORMAT(  '           ' )
 9170 FORMAT(  F10.2 )
 9180 FORMAT(  '          ' )
 9190 FORMAT(  F9.2 )
 9200 FORMAT(  '         ' )
 9210 FORMAT(  2X, A )
 9220 FORMAT(  2X, '--------------------------------------------------',
     $             '-------------------------' )
 9230 FORMAT(  'Min    ', 15X )
 9240 FORMAT(  'Max    ', 15X )
 9250 FORMAT(  'Mean   ', 15X )
 9260 FORMAT(  1X, '**************************************************',
     $             '****************************' )
*
*     End of DGBTP2.
*
      END
