*
*                  GEMM-Based Level 3 BLAS Benchmark
*                            Double Complex
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
*  The program calls a DOUBLE PRECISION function DSECND with no
*  arguments, which is assumed to return the user time for a process in
*  seconds from some fixed starting-time.
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
      PROGRAM ZGBTIM
*     .. Parameters ..
      INTEGER            NIN, NOUT, NERR, IERR
      PARAMETER        ( NIN = 5, NOUT = 6, NERR = 6 )
      INTEGER            LD, NMAX
      PARAMETER        ( LD = 530, NMAX = LD )
      INTEGER            LLN, LST, LNM
      PARAMETER        ( LLN = 256, LST = 50, LNM = 6 )
      INTEGER            MXTAB, MXOPT, MXTRNS, MXDIM, MXLDA, MXSUB,
     $                   MXRUNS
      PARAMETER        ( MXTAB = 6, MXSUB = 8, MXOPT = 2, MXTRNS = 3,
     $                   MXDIM = 36, MXLDA = 24, MXRUNS = 20 )
      COMPLEX*16         Z11, ALPHA, BETA
      PARAMETER        ( Z11 = ( 1.0D+0, 1.0D+0 ),
     $                   ALPHA = ( 0.9D+0, 0.05D+0 ),
     $                   BETA = ( 1.1D+0, 0.03D+0 ) )
*     .. Local Scalars ..
      INTEGER            I, IB, IE, IX, J, JB, JE, KB, KE,
     $                   NTAB, NSIDE, NUPLO, NTRNS, NDIAG, NDIM1, NDIM2,
     $                   NLDA, NRUNS, RUNS, MATCH
      LOGICAL            ERR1, ERR2, ERR3, ERR4, SUB
*     .. Intrinsic Functions ..
      INTRINSIC          DBLE, DCMPLX
*     .. External Functions ..
      INTEGER            EOLN
      LOGICAL            LSAME, GETWRD
      EXTERNAL           LSAME, GETWRD, EOLN
*     .. External Subroutines ..
      EXTERNAL           ZGBT01, ZGBT02, ZGBTP1, ZGBTP2
*     .. Local Arrays ..
      INTEGER            DIM1( MXDIM ), DIM2( MXDIM ), LDA( MXLDA )
      LOGICAL            SUBCHK( MXSUB ), TABSUB( MXSUB ), TAB( MXTAB )
      COMPLEX*16         A( LD, NMAX ), B( LD, NMAX ), C( LD, NMAX )
      DOUBLE PRECISION   USRES( MXSUB, MXOPT, MXOPT, MXTRNS, MXOPT,
     $                   MXDIM, MXLDA ),
     $                   GBRES( MXSUB, MXOPT, MXOPT, MXTRNS, MXOPT,
     $                   MXDIM, MXLDA ),
     $                   MMRES( MXSUB, MXOPT, MXOPT, MXTRNS, MXOPT,
     $                   MXDIM, MXLDA )
      COMMON           / ZBKCMN / A, B, C, USRES, GBRES, MMRES
      CHARACTER          SIDE( MXOPT ), UPLO( MXOPT ), TRNS( MXTRNS ),
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
     $                   TRNS/ MXTRNS*' ' /, DIAG/ MXOPT*' '/,
     $                   NAME/ 'ZSYMM ', 'ZHEMM ', 'ZSYRK ', 'ZHERK ',
     $                   'ZSYR2K', 'ZHER2K', 'ZTRMM ', 'ZTRSM '/,
     $                   SUB/ .FALSE. /
      DATA  BLANK/'                                                  '/,
     $      LBL  /'                                                  '/
*     ..
*     .. Executable Statements ..
*
*     Read the next non-blank/non-comment line from the input parameter
*     file. Store the line in the variable INLN. The first word (token)
*     of the line is stored in INLN( IB:IE ).
*
   10 READ( NIN, FMT = 9000, END = 210 ) INLN
      IF( .NOT.GETWRD( INLN, LLN, IB, IE ).OR.
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
         IF( GETWRD( INLN( JE+1:LLN ), LLN-JE, KB, KE ) )THEN
            KE = EOLN( INLN( JE+1:LLN ), LLN-JE )
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
   20    IF( GETWRD( INLN( JE+1:LLN ), LLN-JE, KB, KE ) )THEN
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
   30    IF( GETWRD( INLN( JE+1:LLN ), LLN-JE, KB, KE ) )THEN
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
   40    IF( GETWRD( INLN( JE+1:LLN ), LLN-JE, KB, KE ) )THEN
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
   70    IF( GETWRD( INLN( JE+1:LLN ), LLN-JE, KB, KE ) )THEN
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
  100    IF( GETWRD( INLN( JE+1:LLN ), LLN-JE, KB, KE ) )THEN
            JB = JE + KB
            JE = JE + KE
            IF( I.LT.MXTRNS )THEN
               IF( LSAME( INLN( JB:JB ), 'N' ) )THEN
                  DO 110, J = 1, I
                     IF( LSAME( TRNS( J ), 'N' ) ) ERR1 = .TRUE.
  110             CONTINUE
               ELSE IF( LSAME( INLN( JB:JB ), 'T' ) )THEN
                  DO 120, J = 1, I
                     IF( LSAME( TRNS( J ), 'T' ) ) ERR1 = .TRUE.
  120             CONTINUE
               ELSE IF( LSAME( INLN( JB:JB ), 'C' ) )THEN
                  DO 130, J = 1, I
                     IF( LSAME( TRNS( J ), 'C' ) ) ERR1 = .TRUE.
  130             CONTINUE
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
  140    IF( GETWRD( INLN( JE+1:LLN ), LLN-JE, KB, KE ) )THEN
            JB = JE + KB
            JE = JE + KE
            IF( I.LT.MXOPT )THEN
               IF( LSAME( INLN( JB:JB ), 'N' ) )THEN
                  DO 150, J = 1, I
                     IF( LSAME( DIAG( J ), 'N' ) ) ERR1 = .TRUE.
  150             CONTINUE
               ELSE IF( LSAME( INLN( JB:JB ), 'U' ) )THEN
                  DO 160, J = 1, I
                     IF( LSAME( DIAG( J ), 'U' ) ) ERR1 = .TRUE.
  160             CONTINUE
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
               GO TO 140
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
  170    IF( GETWRD( INLN( JE+1:LLN ), LLN-JE, KB, KE ) )THEN
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
               GO TO 170
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
  180    IF( GETWRD( INLN( JE+1:LLN ), LLN-JE, KB, KE ) )THEN
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
               GO TO 180
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
  190    IF( GETWRD( INLN( JE+1:LLN ), LLN-JE, KB, KE ) )THEN
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
               GO TO 190
            END IF
         END IF
         NLDA = I
      ELSE IF( INLN( IB:IE ).EQ.'ZSYMM'.OR.INLN( IB:IE ).EQ.'ZHEMM'.OR.
     $          INLN( IB:IE ).EQ.'ZSYRK'.OR.INLN( IB:IE ).EQ.'ZHERK'.OR.
     $        INLN( IB:IE ).EQ.'ZSYR2K'.OR.INLN( IB:IE ).EQ.'ZHER2K'.OR.
     $        INLN( IB:IE ).EQ.'ZTRMM'.OR.INLN( IB:IE ).EQ.'ZTRSM' )THEN
*
*        Read which routines to time.
*
         MATCH = 0
         DO 200, I = 1, MXSUB
            IF( NAME( I ).EQ.INLN( IB:IB+5 ) )THEN
               MATCH = I
               IF( SUBCHK( MATCH ) )THEN
                  ERR3 = .TRUE.
               END IF
               SUBCHK( MATCH ) = .TRUE.
            END IF
  200    CONTINUE
         IF( GETWRD( INLN( JE+1:LLN ), LLN-JE, KB, KE ) )THEN
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
  210 CONTINUE
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
     $      'routines ZSYMM, ZHEMM, ZSYRK, ZHERK, ZSYR2K, ZHER2K,'
         WRITE( NERR, FMT = * )
     $      'ZTRMM, and ZTRSM are marked to be timed', '.'
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
      DO 230, J = 1, NMAX
         DO 220, I = 1, NMAX
            A( I, J ) = Z11 + DCMPLX( 0.08D+0*DBLE( I+( J-1 )*NMAX )/
     $              DBLE( NMAX*NMAX+1 ), 0.06D+0*DBLE( I+( J-1 )*NMAX )/
     $                                             DBLE( NMAX*NMAX+1 ) )
  220    CONTINUE
  230 CONTINUE
      DO 250, J = 1, NMAX
         DO 240, I = 1, NMAX
            B( I, J ) = Z11 + DCMPLX( 0.04D+0*DBLE( I+( J-1 )*NMAX )/
     $              DBLE( NMAX*NMAX+1 ), 0.02D+0*DBLE( I+( J-1 )*NMAX )/
     $                                             DBLE( NMAX*NMAX+1 ) )
  240    CONTINUE
  250 CONTINUE
*
*     Time the routines and calculate the results.
*
      IF( TAB( 2 ).OR.TAB( 6 ) )THEN
*
*        Time the internal GEMM-Based Level 3 BLAS routines (ZGB02,
*        ZGB03, ZGB04, ZGB05, ZGB06, ZGB07, ZGB08, and ZGB09).
*
         CALL ZGBT01( 'G', TABSUB, SIDE, NSIDE, UPLO, NUPLO, TRNS,
     $                 NTRNS, DIAG, NDIAG, DIM1, DIM2, NDIM1, LDA, NLDA,
     $                      ALPHA, BETA, A, B, C, LD, NMAX, NERR, MXSUB,
     $                        MXOPT, MXTRNS, MXDIM, MXLDA, RUNS, GBRES )
      END IF
      IF( TAB( 1 ).OR.TAB( 3 ).OR.TAB( 5 ).OR.TAB( 6 ) )THEN
*
*        Time the user-supplied Level 3 BLAS library.
*
         CALL ZGBT01( 'U', TABSUB, SIDE, NSIDE, UPLO, NUPLO, TRNS,
     $                 NTRNS, DIAG, NDIAG, DIM1, DIM2, NDIM1, LDA, NLDA,
     $                      ALPHA, BETA, A, B, C, LD, NMAX, NERR, MXSUB,
     $                        MXOPT, MXTRNS, MXDIM, MXLDA, RUNS, USRES )
      END IF
      IF( TAB( 1 ).OR.TAB( 4 ).OR.TAB( 5 ) )THEN
*
*        Time ZGEMM using user specified parameters.
*
         CALL ZGBT02( TABSUB, SIDE, NSIDE, NUPLO, TRNS, NTRNS, NDIAG,
     $                        DIM1, DIM2, NDIM1, LDA, NLDA, ALPHA, BETA,
     $                            A, B, C, LD, NMAX, NERR, MXSUB, MXOPT,
     $                               MXTRNS, MXDIM, MXLDA, RUNS, MMRES )
      END IF
      IF( TAB( 1 ) )THEN
*
*        Calculate and print the collected benchmark result.
*
         CALL ZGBTP1( TABSUB, SIDE, NSIDE, UPLO, NUPLO, TRNS, NTRNS,
     $                  DIAG, NDIAG, DIM1, DIM2, NDIM1, LDA, NLDA, NOUT,
     $                   NERR, MXSUB, MXOPT, MXTRNS, MXDIM, MXLDA, RUNS,
     $                                  ALPHA, BETA, LBL, USRES, MMRES )
      END IF
      IF( TAB( 2 ).OR.TAB( 3 ).OR.TAB( 4 ).OR.TAB( 5 ).OR.TAB( 6 ) )THEN
*
*        Calculate and print the results of TAB choice 2 - 6.
*
         CALL ZGBTP2( TAB, TABSUB, SIDE, NSIDE, UPLO, NUPLO, TRNS,
     $                 NTRNS, DIAG, NDIAG, DIM1, DIM2, NDIM1, LDA, NLDA,
     $                  NOUT, MXTAB, MXSUB, MXOPT, MXTRNS, MXDIM, MXLDA,
     $                     RUNS, ALPHA, BETA, LBL, USRES, MMRES, GBRES )
      END IF
*
      STOP
*
 9000 FORMAT( A )
 9010 FORMAT( 1X, A, I3, A )
 9020 FORMAT( I50 )
*
*     End of ZGBTIM.
*
      END
