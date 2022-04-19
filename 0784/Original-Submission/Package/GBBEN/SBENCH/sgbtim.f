*
*                  GEMM-Based Level 3 BLAS Benchmark
*                          REAL
*
*  The GEMM-Based Level 3 BLAS Benchmark is a tool for performance
*  evaluation of Level 3 BLAS kernel programs. With the announcement of
*  LAPACK, the need for high performance Level 3 BLAS kernels became
*  apparent. LAPACK is based on calls to the Level 3 BLAS kernels. This
*  benchmark measures and compares performance of a set of user supplied
*  Level 3 BLAS implementations and of the GEMM-Based Level 3 BLAS
*  implementations permanently included in the benchmark. The purpose of
*  the benchmark is to facilitate the user in determining the quality of
*  different Level 3 BLAS implementations. The included GEMM-Based
*  Level 3 BLAS routines provide a lower limit on the performance to be
*  expected from a highly optimized Level 3 BLAS library.
*
*  All routines are written in Fortran 77 for portability. No changes to
*  the code should be necessary in order to run the programs correctly
*  on different target machines. In fact, we strongly recommend the user
*  to avoided changes, except to the user specified parameters and to
*  UNIT numbers for input and output communication. This will ensure
*  that performance results from different target machines are
*  comparable.
*
*  The program calls a REAL function SECOND with no
*  arguments, which is assumed to return the user time for a process in
*  seconds from some fixed starting-time.
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
      PROGRAM SGBTIM
*     .. Parameters ..
      INTEGER            NIN, NOUT, NERR, IERR
      PARAMETER        ( NIN = 5, NOUT = 6, NERR = 6 )
      INTEGER            LD, NMAX
      PARAMETER        ( LD = 530, NMAX = LD )
      INTEGER            LLN, LST, LNM
      PARAMETER        ( LLN = 256, LST = 50, LNM = 6 )
      INTEGER            MXTAB, MXOPT, MXDIM, MXLDA, MXSUB, MXRUNS
      PARAMETER        ( MXTAB = 6, MXSUB = 5, MXOPT = 2, MXDIM = 36,
     $                   MXLDA = 24, MXRUNS = 20 )
      REAL               ONE, ALPHA, BETA
      PARAMETER        ( ONE = 1.0E+0, ALPHA = 0.9E+0, BETA = 1.1E+0 )
*     .. Local Scalars ..
      INTEGER            I, IB, IE, IX, J, JB, JE, KB, KE,
     $                   NTAB, NSIDE, NUPLO, NTRNS, NDIAG, NDIM1, NDIM2,
     $                   NLDA, NRUNS, RUNS, MATCH
      LOGICAL            ERR1, ERR2, ERR3, ERR4, SUB
*     .. Intrinsic Functions ..
      INTRINSIC          REAL
*     .. External Functions ..
      INTEGER            EOLN
      LOGICAL            LSAME, GETWRD
      EXTERNAL           LSAME, GETWRD, EOLN
*     .. External Subroutines ..
      EXTERNAL           SGBT01, SGBT02, SGBTP1, SGBTP2
*     .. Local Arrays ..
      INTEGER            DIM1( MXDIM ), DIM2( MXDIM ), LDA( MXLDA )
      LOGICAL            SUBCHK( MXSUB ), TABSUB( MXSUB ), TAB( MXTAB )
      REAL               A( LD, NMAX ), B( LD, NMAX ), C( LD, NMAX ),
     $                   USRES( MXSUB, MXOPT, MXOPT, MXOPT, MXOPT,
     $                   MXDIM, MXLDA ),
     $                   GBRES( MXSUB, MXOPT, MXOPT, MXOPT, MXOPT,
     $                   MXDIM, MXLDA ),
     $                   MMRES( MXSUB, MXOPT, MXOPT, MXOPT, MXOPT,
     $                   MXDIM, MXLDA )
      COMMON           / SBKCMN / A, B, C, USRES, GBRES, MMRES
      CHARACTER          SIDE( MXOPT ), UPLO( MXOPT ), TRNS( MXOPT ),
     $                   DIAG( MXOPT )
      CHARACTER          INLN*( LLN ), INSTR*( LST ), BLANK*( LST ), 
     $                   LBL*( LST ), NAME( MXSUB )*( LNM )
      CHARACTER          INLNA( LLN )
      EQUIVALENCE      ( INLN, INLNA )
*     .. Data statements ..
      DATA               NTAB/ 0 /, NRUNS/ 0 /, NSIDE/ 0 /, NUPLO/ 0 /,
     $                   NTRNS/ 0 /, NDIAG/ 0 /, NDIM1/ 0 /, NDIM2/ 0 /,
     $                   NLDA/ 0 /
      DATA               TAB/ MXTAB*.FALSE. /, TABSUB/ MXSUB*.FALSE. /,
     $                   SUBCHK/ MXSUB*.FALSE. /,
     $                   SIDE/ MXOPT*' ' /, UPLO/ MXOPT*' '/,
     $                   TRNS/ MXOPT*' ' /, DIAG/ MXOPT*' '/,
     $                   NAME/ 'SSYMM ', 'SSYRK ', 'SSYR2K', 'STRMM ',
     $                   'STRSM '/, SUB/ .FALSE. /
      DATA  BLANK/'                                                  '/,
     $      LBL  /'                                                  '/
*     ..
*     .. Executable Statements ..
*
*     Read the next non-blank/non-comment line from the input parameter
*     file. Store the line in the variable INLN. The first word (token)
*     of the line is stored in INLN( IB:IE ).
*
   10 READ( NIN, FMT = 9000, END = 200 ) INLN
      IF( .NOT.GETWRD( INLNA, LLN, IB, IE ).OR.
     $                                     ( INLN( 1:1 ).EQ.'*' ) )THEN
         GO TO 10
      END IF
*
*     If INLN( IB:IE ) contains the key word for a parameter, then read
*     and store the parameter values given on the same line of the input
*     file, after the key word.
*
      JB = IB
      JE = IE
      I = 0
      ERR1 = .FALSE.
      ERR2 = .FALSE.
      ERR3 = .FALSE.
      ERR4 = .FALSE.
*
*     Read the parameters from the line INLN.
*
      IF( INLN( JB:JE ).EQ.'LBL' )THEN
*
*        Read the label of this test.
*
         IF( LBL.NE.BLANK )THEN
            ERR3 = .TRUE.
         END IF
         IF( GETWRD( INLNA( JE+1 ), LLN-JE-1, KB, KE ) )THEN
            KE = EOLN( INLNA( JE+1 ), LLN-JE-1 )
            JB = JE + KB
            JE = JE + KE
            IF( JE-JB+1.GT.LST )THEN
               ERR4 = .TRUE.
            ELSE
               LBL = INLN( JB:JE )
            END IF
         END IF
         I = 1
      ELSE IF( INLN( JB:JE ).EQ.'TAB' )THEN
*
*        Read which tests to be made.
*
         IF( NTAB.NE.0 )THEN
            ERR3 = .TRUE.
         END IF
   20    IF( GETWRD( INLNA( JE+1 ), LLN-JE-1, KB, KE ) )THEN
            JB = JE + KB
            JE = JE + KE
            I = I + 1
            IF( I.LE.MXTAB )THEN
               INSTR = BLANK( 1:LST-JE+JB-1 )//INLN( JB:JE )
               READ( INSTR, FMT = 9020, IOSTAT = IERR ) IX
               IF( IERR.GT.0.OR.IX.LT.1.OR.IX.GT.MXTAB )THEN
                  ERR1 = .TRUE.
               END IF
               IF( TAB( IX ) )THEN
                   ERR1 = .TRUE.
               END IF
               TAB( IX ) = .TRUE.
            ELSE
               ERR1 = .TRUE.
            END IF
            IF( .NOT.ERR1.AND.JE.LT.LLN )THEN
               GO TO 20
            END IF
         END IF
         NTAB = I
      ELSE IF( INLN( JB:JE ).EQ.'RUNS' )THEN
*
*        Read the number of times each problem is to be executed. The
*        final performance results are computed using the best timing
*        result for each problem.
*
         IF( NRUNS.NE.0 )THEN
            ERR3 = .TRUE.
         END IF
   30    IF( GETWRD( INLNA( JE+1 ), LLN-JE-1, KB, KE ) )THEN
            JB = JE + KB
            JE = JE + KE
            I = I + 1
            IF( I.LE.1 )THEN
               INSTR = BLANK( 1:LST-JE+JB-1 )//INLN( JB:JE )
               READ( INSTR, FMT = 9020, IOSTAT = IERR ) RUNS
               IF( IERR.GT.0.OR.RUNS.LT.1.OR.RUNS.GT.MXRUNS )THEN
                  ERR1 = .TRUE.
               END IF
            ELSE
               ERR1 = .TRUE.
            END IF
            IF( .NOT.ERR1.AND.JE.LT.LLN )THEN
               GO TO 30
            END IF
         END IF
         NRUNS = I
      ELSE IF( INLN( IB:IE ).EQ.'SIDE' )THEN
*
*        Read the values for SIDE.
*
         IF( NSIDE.NE.0 )THEN
            ERR3 = .TRUE.
         END IF
   40    IF( GETWRD( INLNA( JE+1 ), LLN-JE-1, KB, KE ) )THEN
            JB = JE + KB
            JE = JE + KE
            IF( I.LT.MXOPT )THEN
               IF( LSAME( INLN( JB:JB ), 'L' ) )THEN
                  DO 50, J = 1, I
                     IF( LSAME( SIDE( J ), 'L' ) ) ERR1 = .TRUE.
   50             CONTINUE
               ELSE IF( LSAME( INLN( JB:JB ), 'R' ) )THEN
                  DO 60, J = 1, I
                     IF( LSAME( SIDE( J ), 'R' ) ) ERR1 = .TRUE.
   60             CONTINUE
               ELSE
                  ERR1 = .TRUE.
               END IF
            ELSE
               ERR1 = .TRUE.
            END IF
            IF( .NOT.ERR1 )THEN
               I = I + 1
               SIDE( I ) = INLN( JB:JB )
            END IF
            IF( .NOT.ERR1.AND.JE.LT.LLN )THEN
               GO TO 40
            END IF
         END IF
         NSIDE = I
      ELSE IF( INLN( IB:IE ).EQ.'UPLO' )THEN
*
*        Read the values for UPLO.
*
         IF( NUPLO.NE.0 )THEN
            ERR3 = .TRUE.
         END IF
   70    IF( GETWRD( INLNA( JE+1 ), LLN-JE-1, KB, KE ) )THEN
            JB = JE + KB
            JE = JE + KE
            IF( I.LT.MXOPT )THEN
               IF( LSAME( INLN( JB:JB ), 'U' ) )THEN
                  DO 80, J = 1, I
                     IF( LSAME( UPLO( J ), 'U' ) ) ERR1 = .TRUE.
   80             CONTINUE
               ELSE IF( LSAME( INLN( JB:JB ), 'L' ) )THEN
                  DO 90, J = 1, I
                     IF( LSAME( UPLO( J ), 'L' ) ) ERR1 = .TRUE.
   90             CONTINUE
               ELSE
                  ERR1 = .TRUE.
               END IF
            ELSE
               ERR1 = .TRUE.
            END IF
            IF( .NOT.ERR1 )THEN
               I = I + 1
               UPLO( I ) = INLN( JB:JB )
            END IF
            IF( .NOT.ERR1.AND.JE.LT.LLN )THEN
               GO TO 70
            END IF
         END IF
         NUPLO = I
      ELSE IF( INLN( IB:IE ).EQ.'TRANS' )THEN
*
*        Read the values for TRANS.
*
         IF( NTRNS.NE.0 )THEN
            ERR3 = .TRUE.
         END IF
  100    IF( GETWRD( INLNA( JE+1 ), LLN-JE-1, KB, KE ) )THEN
            JB = JE + KB
            JE = JE + KE
            IF( I.LT.MXOPT )THEN
               IF( LSAME( INLN( JB:JB ), 'N' ) )THEN
                  DO 110, J = 1, I
                     IF( LSAME( TRNS( J ), 'N' ) ) ERR1 = .TRUE.
  110             CONTINUE
               ELSE IF( LSAME( INLN( JB:JB ), 'T' ) )THEN
                  DO 120, J = 1, I
                     IF( LSAME( TRNS( J ), 'T' ) ) ERR1 = .TRUE.
  120             CONTINUE
               ELSE
                  ERR1 = .TRUE.
               END IF
            ELSE
               ERR1 = .TRUE.
            END IF
            IF( .NOT.ERR1 )THEN
               I = I + 1
               TRNS( I ) = INLN( JB:JB )
            END IF
            IF( .NOT.ERR1.AND.JE.LT.LLN )THEN
               GO TO 100
            END IF
         END IF
         NTRNS = I
      ELSE IF( INLN( IB:IE ).EQ.'DIAG' )THEN
*
*        Read the values for DIAG.
*
         IF( NDIAG.NE.0 )THEN
            ERR3 = .TRUE.
         END IF
  130    IF( GETWRD( INLNA( JE+1 ), LLN-JE-1, KB, KE ) )THEN
            JB = JE + KB
            JE = JE + KE
            IF( I.LT.MXOPT )THEN
               IF( LSAME( INLN( JB:JB ), 'N' ) )THEN
                  DO 140, J = 1, I
                     IF( LSAME( DIAG( J ), 'N' ) ) ERR1 = .TRUE.
  140             CONTINUE
               ELSE IF( LSAME( INLN( JB:JB ), 'U' ) )THEN
                  DO 150, J = 1, I
                     IF( LSAME( DIAG( J ), 'U' ) ) ERR1 = .TRUE.
  150             CONTINUE
               ELSE
                  ERR1 = .TRUE.
               END IF
            ELSE
               ERR1 = .TRUE.
            END IF
            IF( .NOT.ERR1 )THEN
               I = I + 1
               DIAG( I ) = INLN( JB:JB )
            END IF
            IF( .NOT.ERR1.AND.JE.LT.LLN )THEN
               GO TO 130
            END IF
         END IF
         NDIAG = I
      ELSE IF( INLN( IB:IE ).EQ.'DIM1' )THEN
*
*        Read the values for the first matrix dimension (DIM1).
*
         IF( NDIM1.NE.0 )THEN
            ERR3 = .TRUE.
         END IF
  160    IF( GETWRD( INLNA( JE+1 ), LLN-JE-1, KB, KE ) )THEN
            JB = JE + KB
            JE = JE + KE
            I = I + 1
            IF( I.LE.MXDIM )THEN
               INSTR = BLANK( 1:LST-JE+JB-1 )//INLN( JB:JE )
               READ( INSTR, FMT = 9020, IOSTAT = IERR ) DIM1( I )
               IF( IERR.GT.0.OR.DIM1( I ).LT.0 )THEN
                  ERR1 = .TRUE.
               END IF
               IF( DIM1( I ).GT.NMAX )THEN
                  ERR2 = .TRUE.
               END IF
            ELSE
               ERR1 = .TRUE.
            END IF
            IF( .NOT.( ERR1.OR.ERR2 ).AND.JE.LT.LLN )THEN
               GO TO 160
            END IF
         END IF
         NDIM1 = I
      ELSE IF( INLN( IB:IE ).EQ.'DIM2' )THEN
*
*        Read the values for the second matrix dimension (DIM2).
*
         IF( NDIM2.NE.0 )THEN
            ERR3 = .TRUE.
         END IF
  170    IF( GETWRD( INLNA( JE+1 ), LLN-JE-1, KB, KE ) )THEN
            JB = JE + KB
            JE = JE + KE
            I = I + 1
            IF( I.LE.MXDIM )THEN
               INSTR = BLANK( 1:LST-JE+JB-1 )//INLN( JB:JE )
               READ( INSTR, FMT = 9020, IOSTAT = IERR ) DIM2( I )
               IF( IERR.GT.0.OR.DIM2( I ).LT.0 )THEN
                  ERR1 = .TRUE.
               END IF
               IF( DIM2( I ).GT.NMAX )THEN
                  ERR2 = .TRUE.
               END IF
            ELSE
               ERR1 = .TRUE.
            END IF
            IF( .NOT.( ERR1.OR.ERR2 ).AND.JE.LT.LLN )THEN
               GO TO 170
            END IF
         END IF
         NDIM2 = I
      ELSE IF( INLN( IB:IE ).EQ.'LDA' )THEN
*
*        Read the values for leading dimension (LDA).
*
         IF( NLDA.NE.0 )THEN
            ERR3 = .TRUE.
         END IF
  180    IF( GETWRD( INLNA( JE+1 ), LLN-JE-1, KB, KE ) )THEN
            JB = JE + KB
            JE = JE + KE
            I = I + 1
            IF( I.LE.MXLDA )THEN
               INSTR = BLANK( 1:LST-JE+JB-1 )//INLN( JB:JE )
               READ( INSTR, FMT = 9020, IOSTAT = IERR ) LDA( I )
               IF( IERR.GT.0.OR.LDA( I ).LT.0 )THEN
                  ERR1 = .TRUE.
               END IF
               IF( LDA( I ).GT.NMAX )THEN
                  ERR2 = .TRUE.
               END IF
            ELSE
               ERR1 = .TRUE.
            END IF
            IF( .NOT.( ERR1.OR.ERR2 ).AND.JE.LT.LLN )THEN
               GO TO 180
            END IF
         END IF
         NLDA = I
      ELSE IF( INLN( IB:IE ).EQ.'SSYMM'.OR.INLN( IB:IE ).EQ.'SSYRK'.OR.
     $         INLN( IB:IE ).EQ.'SSYR2K'.OR.INLN( IB:IE ).EQ.'STRMM'.OR.
     $                                    INLN( IB:IE ).EQ.'STRSM' )THEN
*
*        Read which routines to time.
*
         MATCH = 0
         DO 190, I = 1, MXSUB
            IF( NAME( I ).EQ.INLN( IB:IB+5 ) )THEN
               MATCH = I
               IF( SUBCHK( MATCH ) )THEN
                  ERR3 = .TRUE.
               END IF
               SUBCHK( MATCH ) = .TRUE.
            END IF
  190    CONTINUE
         IF( GETWRD( INLNA( JE+1 ), LLN-JE-1, KB, KE ) )THEN
            JB = JE + KB
            JE = JE + KE
*
*           Time the routine if the first non-blank character
*           INLN( JB:JB ) is 'T' or 't'.
*
            TABSUB( MATCH ) = LSAME( INLN( JB:JB ), 'T' )
            IF( .NOT.( TABSUB( MATCH ).OR.
     $                               LSAME( INLN( JB:JB ), 'F' ) ) )THEN
               ERR1 = .TRUE.
            END IF
            SUB = SUB.OR.TABSUB( MATCH )
            I = 1
         ELSE
            I = 0
         END IF
      ELSE
         WRITE( NERR, FMT = * )
         WRITE( NERR, FMT = * )
     $                   'Error: Unknown parameter ', INLN( IB:IE ), '.'
         WRITE( NERR, FMT = * )
         STOP
      END IF
*
      IF( I.EQ.0 )THEN
         WRITE( NERR, FMT = * )
         WRITE( NERR, FMT = * )
     $                    'Error: No values or erroneous values given ',
     $                          'for the parameter ', INLN( IB:IE ), '.'
         WRITE( NERR, FMT = * )
         STOP
      ELSE IF( ERR1 )THEN
         WRITE( NERR, FMT = * )
         WRITE( NERR, FMT = * )
     $           'Erroneus value or too many values for the parameter ',
     $                                                INLN( IB:IE ), '.'
         WRITE( NERR, FMT = * )
         STOP
      ELSE IF( ERR2 )THEN
         WRITE( NERR, FMT = * )
         WRITE( NERR, FMT = * )
     $        'Value too large for ', INLN( IB:IE ), '. Max ', NMAX, '.'
         WRITE( NERR, FMT = * )
         STOP
      ELSE IF( ERR3 )THEN
         WRITE( NERR, FMT = * )
         WRITE( NERR, FMT = * )
     $                'Multiple specifications of the input parameter ',
     $                                                INLN( IB:IE ), '.'
         WRITE( NERR, FMT = * )
         STOP
      ELSE IF( ERR4 )THEN
         WRITE( NERR, FMT = * )
         WRITE( NERR, FMT = 9010 )
     $      'The label of this test is too long. Max ', LST,
     $                                                    ' characters.'
         WRITE( NERR, FMT = * )
         STOP
      END IF
      GO TO 10
*
  200 CONTINUE
      IF( NTAB.LE.0 )THEN
         WRITE( NERR, FMT = * )
         WRITE( NERR, FMT = * )
     $      'Error: No results are chosen to be presented'
         WRITE( NERR, FMT = * )
     $      '       (see the parameter TAB).'
         WRITE( NERR, FMT = * )
         STOP
      END IF
      IF( ( TAB( 2 ).OR.TAB( 3 ).OR.TAB( 4 ).OR.TAB( 5 ).OR.TAB( 6 ) )
     $                 .AND.( NRUNS.LE.0.OR.NSIDE.LE.0.OR.NUPLO.LE.0.OR.
     $                        NTRNS.LE.0.OR.NDIAG.LE.0.OR.NDIM1.LE.0.OR.
     $                   NDIM2.LE.0.OR.NLDA.LE.0.OR.( .NOT.SUB ) ) )THEN
         WRITE( NERR, FMT = * )
         WRITE( NERR, FMT = * )
     $      'Error: A parameter, or values for a parameter, is missing.'
         WRITE( NERR, FMT = * )
     $      'One (or more) of the input parameters RUNS, SIDE, UPLO,'
         WRITE( NERR, FMT = * )
     $      'TRANS, DIAG, DIM1, DIM2, LDA are missing, or none of the'
         WRITE( NERR, FMT = * )
     $      'routines SSYMM, SSYRK, SSYR2K, STRMM, and STRSM are marked'
         WRITE( NERR, FMT = * )
     $      'to be timed', '.'
         WRITE( NERR, FMT = * )
         STOP
      END IF
      IF( NDIM1.NE.NDIM2 )THEN
         WRITE( NERR, FMT = * )
         WRITE( NERR, FMT = * )
     $                         'Error: Different number of dimensions ',
     $                                          'for DIM1 and DIM2', '.'
         WRITE( NERR, FMT = * )
         STOP
      END IF
*
*     Initialize the matrices A and B.
*
      DO 220, J = 1, NMAX
         DO 210, I = 1, NMAX
            A( I, J ) = ONE + 0.08E+0*REAL( I+( J-1 )*NMAX )/
     $                                               REAL( NMAX*NMAX+1 )
  210    CONTINUE
  220 CONTINUE
      DO 240, J = 1, NMAX
         DO 230, I = 1, NMAX
            B( I, J ) = ONE + 0.04E+0*REAL( I+( J-1 )*NMAX )/
     $                                               REAL( NMAX*NMAX+1 )
  230    CONTINUE
  240 CONTINUE
*
*     Time the routines and calculate the results.
*
      IF( TAB( 2 ).OR.TAB( 6 ) )THEN
*
*        Time the internal GEMM-Based Level 3 BLAS routines (SGB02,
*        SGB04, SGB08, SGB08, and SGB09).
*
         CALL SGBT01( 'G', TABSUB, SIDE, NSIDE, UPLO, NUPLO, TRNS,
     $                 NTRNS, DIAG, NDIAG, DIM1, DIM2, NDIM1, LDA, NLDA,
     $                      ALPHA, BETA, A, B, C, LD, NMAX, NERR, MXSUB,
     $                                MXOPT, MXDIM, MXLDA, RUNS, GBRES )
      END IF
      IF( TAB( 1 ).OR.TAB( 3 ).OR.TAB( 5 ).OR.TAB( 6 ) )THEN
*
*        Time the user-supplied Level 3 BLAS library.
*
         CALL SGBT01( 'U', TABSUB, SIDE, NSIDE, UPLO, NUPLO, TRNS,
     $                 NTRNS, DIAG, NDIAG, DIM1, DIM2, NDIM1, LDA, NLDA,
     $                      ALPHA, BETA, A, B, C, LD, NMAX, NERR, MXSUB,
     $                                MXOPT, MXDIM, MXLDA, RUNS, USRES )
      END IF
      IF( TAB( 1 ).OR.TAB( 4 ).OR.TAB( 5 ) )THEN
*
*        Time SGEMM using user specified parameters.
*
         CALL SGBT02( TABSUB, SIDE, NSIDE, NUPLO, TRNS, NTRNS, NDIAG,
     $                        DIM1, DIM2, NDIM1, LDA, NLDA, ALPHA, BETA,
     $                            A, B, C, LD, NMAX, NERR, MXSUB, MXOPT,
     $                                       MXDIM, MXLDA, RUNS, MMRES )
      END IF
      IF( TAB( 1 ) )THEN
*
*        Calculate and print the collected benchmark result.
*
         CALL SGBTP1( TABSUB, SIDE, NSIDE, UPLO, NUPLO, TRNS, NTRNS,
     $                  DIAG, NDIAG, DIM1, DIM2, NDIM1, LDA, NLDA, NOUT,
     $                           NERR, MXSUB, MXOPT, MXDIM, MXLDA, RUNS,
     $                                  ALPHA, BETA, LBL, USRES, MMRES )
      END IF
      IF( TAB( 2 ).OR.TAB( 3 ).OR.TAB( 4 ).OR.TAB( 5 ).OR.TAB( 6 ) )THEN
*
*        Calculate and print the results of TAB choice 2 - 6.
*
         CALL SGBTP2( TAB, TABSUB, SIDE, NSIDE, UPLO, NUPLO, TRNS,
     $                 NTRNS, DIAG, NDIAG, DIM1, DIM2, NDIM1, LDA, NLDA,
     $                    NOUT, MXTAB, MXSUB, MXOPT, MXDIM, MXLDA, RUNS,
     $                           ALPHA, BETA, LBL, USRES, MMRES, GBRES )
      END IF
*
      STOP
*
 9000 FORMAT( A )
 9010 FORMAT( 1X, A, I3, A )
 9020 FORMAT( I50 )
*
*     End of SGBTIM.
*
      END
