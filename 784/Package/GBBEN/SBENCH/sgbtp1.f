      SUBROUTINE SGBTP1( TABSUB, SIDE, NSIDE, UPLO, NUPLO, TRNS, NTRNS,
     $                   DIAG, NDIAG, DIM1, DIM2, NDIM, LDA, NLDA, NOUT,
     $                           NERR, MXSUB, MXOPT, MXDIM, MXLDA, RUNS,
     $                                  ALPHA, BETA, LBL, USRES, MMRES )
*     .. Scalar Arguments ..
      INTEGER            NOUT, NERR,
     $                   NSIDE, NUPLO, NTRNS, NDIAG, NDIM, NLDA,
     $                   MXSUB, MXOPT, MXDIM, MXLDA, RUNS
      REAL               ALPHA, BETA
*     .. Parameters ..
      INTEGER            LST
      PARAMETER        ( LST = 50 )
*     .. Array Arguments ..
      LOGICAL            TABSUB( MXSUB )
      CHARACTER          LBL*( LST )
      CHARACTER          SIDE( MXOPT ), UPLO( MXOPT ), TRNS( MXOPT ),
     $                   DIAG( MXOPT )
      INTEGER            DIM1( MXDIM ), DIM2( MXDIM ), LDA( MXLDA )
      REAL               USRES( MXSUB, MXOPT, MXOPT, MXOPT, MXOPT,
     $                   MXDIM, MXLDA ),
     $                   MMRES( MXSUB, MXOPT, MXOPT, MXOPT, MXOPT,
     $                   MXDIM, MXLDA )
*
*
*  SGBTP1 prints the collected benchmark result which is calculated from
*  performance results of the user-supplied Level 3 routines for
*  problems specified in the input file. The result consists of a tuple
*  ( x, y ), where x is the mean value of the GEMM-Efficiency and y is
*  the mean value of the performance of SGEMM in megaflops. SGEMM is
*  timed for problems corresponding to those specified for the remaining
*  Level 3 routines.
*
*  The purpose of the collected benchmark result is to provide an
*  overall performance measure of the user-supplied Level 3 BLAS
*  routines. The intention is to expose the capacity of the target
*  machine for these kinds of problems and to show how well the routines
*  utilize the machine. Furthermore, the collected result is intended to
*  be easy to compare between different target machines. See the README
*  and INSTALL files for further information.
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
      INTEGER            I, D, L, NTIM, OP1, OP2, OP3
      REAL               SPEED, EFF, MM, MMSUM, EFSUM
*     .. Intrinsic Functions ..
      INTRINSIC          REAL
*     .. Parameters ..
      REAL               ZERO
      INTEGER            MXBSUB
      PARAMETER        ( ZERO = 0.0E+0, MXBSUB = 5 )
*     ..
*     .. Executable Statements ..
      IF( MXSUB.GT.MXBSUB )THEN
         WRITE( NERR, FMT = 9000 )
         STOP
      END IF
*
      MMSUM = ZERO
      EFSUM = ZERO
      NTIM = 0
*
*     ------ Stop indentation ------
*
      DO 50, L = 1, NLDA
      DO 40, OP1 = 1, NSIDE
      DO 30, OP2 = 1, NUPLO
      DO 20, OP3 = 1, NTRNS
      DO 10, D = 1, NDIM
*
*     ------ Continue indentation ------
*
*
*     Compute the sum of the performance of SGEMM in megaflops (MMSUM)
*     and the sum of the GEMM-Efficiency (EFSUM).
*
      IF( TABSUB( 1 ).AND.OP3.EQ.1 )THEN
         MM = MMRES( 1, OP1, OP2, 1, 1, D, L )
         MMSUM = MMSUM + MM
         IF( MM.GT.ZERO )THEN
            EFSUM = EFSUM + USRES( 1, OP1, OP2, 1, 1, D, L )/MM
         ELSE
            WRITE( NERR, FMT = 9010 )
            STOP
         END IF
         NTIM = NTIM + 1
      END IF
      IF( TABSUB( 2 ).AND.OP1.EQ.1 )THEN
         MM = MMRES( 2, 1, OP2, OP3, 1, D, L )
         MMSUM = MMSUM + MM
         IF( MM.GT.ZERO )THEN
            EFSUM = EFSUM + USRES( 2, 1, OP2, OP3, 1, D, L )/MM
         ELSE
            WRITE( NERR, FMT = 9010 )
            STOP
         END IF
         NTIM = NTIM + 1
      END IF
      IF( TABSUB( 3 ).AND.OP1.EQ.1 )THEN
         MM = MMRES( 3, 1, OP2, OP3, 1, D, L )
         MMSUM = MMSUM + MM
         IF( MM.GT.ZERO )THEN
            EFSUM = EFSUM + USRES( 3, 1, OP2, OP3, 1, D, L )/MM
         ELSE
            WRITE( NERR, FMT = 9010 )
            STOP
         END IF
         NTIM = NTIM + 1
      END IF
      IF( TABSUB( 4 ) )THEN
         MM = MMRES( 4, OP1, OP2, OP3, 1, D, L )
         MMSUM = MMSUM + MM
         IF( MM.GT.ZERO )THEN
            EFSUM = EFSUM + USRES( 4, OP1, OP2, OP3, 1, D, L )/MM
         ELSE
            WRITE( NERR, FMT = 9010 )
            STOP
         END IF
         NTIM = NTIM + 1
      END IF
      IF( TABSUB( 5 ) )THEN
         MM = MMRES( 5, OP1, OP2, OP3, 1, D, L )
         MMSUM = MMSUM + MM
         IF( MM.GT.ZERO )THEN
            EFSUM = EFSUM + USRES( 5, OP1, OP2, OP3, 1, D, L )/MM
         ELSE
            WRITE( NERR, FMT = 9010 )
            STOP
         END IF
         NTIM = NTIM + 1
      END IF
*
*     ------ Stop indentation ------
*
   10 CONTINUE
   20 CONTINUE
   30 CONTINUE
   40 CONTINUE
   50 CONTINUE
*
*     ------ Continue indentation ------
*
*
*     Compute the collected benchmark result ( x, y ) as the mean value
*     of the GEMM-Efficiency ( x ) and the mean value of the performance
*     of SGEMM in megaflops ( y ).
*
      SPEED = MMSUM/REAL( NTIM )
      EFF = EFSUM/REAL( NTIM )
*
*     Print an introduction and the collected benchmark result.
*
      WRITE( NOUT, FMT = * )
      WRITE( NOUT, FMT = 9020 )
      WRITE( NOUT, FMT = 9030 )
      WRITE( NOUT, FMT = * )
      WRITE( NOUT, FMT = 9040 ) RUNS
      WRITE( NOUT, FMT = * )
      WRITE( NOUT, FMT = 9050 )
      WRITE( NOUT, FMT = 9060 ) 'SIDE   ', ( SIDE( I ), I = 1, NSIDE )
      WRITE( NOUT, FMT = 9060 ) 'UPLO   ', ( UPLO( I ), I = 1, NUPLO )
      WRITE( NOUT, FMT = 9060 ) 'TRANS  ', ( TRNS( I ), I = 1, NTRNS )
      WRITE( NOUT, FMT = 9060 ) 'DIAG   ', ( DIAG( I ), I = 1, NDIAG )
      WRITE( NOUT, FMT = 9070 ) 'DIM1   ', ( DIM1( I ), I = 1, NDIM )
      WRITE( NOUT, FMT = 9070 ) 'DIM2   ', ( DIM2( I ), I = 1, NDIM )
      WRITE( NOUT, FMT = 9070 ) 'LDA    ', ( LDA( I ), I = 1, NLDA )
      WRITE( NOUT, FMT = 9080 ) 'ALPHA  ', ALPHA
      WRITE( NOUT, FMT = 9080 ) 'BETA   ', BETA
      WRITE( NOUT, FMT = * )
      WRITE( NOUT, FMT = 9090 ) LBL
      WRITE( NOUT, FMT = 9100 ) EFF, SPEED
      WRITE( NOUT, FMT = * )
      WRITE( NOUT, FMT = 9110 )
*
      RETURN
*
*     Print formats.
*
 9000 FORMAT(  1X, 'Error: The collected benchmark result could not ',
     $             'be obtained.',/,
     $         1X, 'The value for the input parameter MXSUB is too ',
     $             'large.' )
 9010 FORMAT(  1X, 'Error: The collected benchmark result could not ',
     $             'be obtained.',/,
     $         1X, 'Execution time for SGEMM is zero.' )
 9020 FORMAT( 17X, '****  GEMM-Based Level 3 BLAS Benchmark  ****' )
 9030 FORMAT( 27X, 'Collected Benchmark Result',/,
     $        32X, '      Real      ' )
 9040 FORMAT(  2X, 'The collected benchmark result is a tuple ',
     $             '( x, y ) where x is the mean',/,
     $         2X, 'value of the GEMM-Efficiency and y is the mean ',
     $             'value of the performance',/,
     $         2X, 'of SGEMM in megaflops (see the README file). The ',
     $             'benchmark result is',/,
     $         2X, 'based on the shortest of', I3,' runs for each ',
     $             'problem configuration.' )
 9050 FORMAT(  8X, 'Input parameters.' )
 9060 FORMAT(  8X, A, '   ', 10( A, ' ' ) )
 9070 FORMAT(  8X, A, 1X, 12( I5 ), 2( /, 16X, 12( I5 ) ) )
 9080 FORMAT(  8X, A, F6.1 )
 9090 FORMAT(  8X, 'Test label:         ', A )
 9100 FORMAT(  8X, 'Collected result:   (', F7.2,',', F9.1,'   )' )
 9110 FORMAT(  1X, '**************************************************',
     $             '****************************' )
*
*     End of SGBTP1.
*
      END
