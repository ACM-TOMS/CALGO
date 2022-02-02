      PROGRAM STIMAA
*
*  -- LAPACK timing routine (version 2.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     September 30, 1994
*
*     Rewritten to include the timing of rrqr code.
*
*  Purpose
*  =======
*
*  An annotated example of a data file can be obtained by deleting the
*  first 3 characters from the following lines:
*  LAPACK timing, REAL square matrices
*  5                                Number of values of M
*  100 200 300 400 500              Values of M (row dimension)
*  5                                Number of values of N
*  100 200 300 400 500              Values of N (column dimension)
*  2                                Number of values of K
*  100 400                          Values of K
*  5                                Number of values of NB
*  1 16  32  48  64                 Values of NB (blocksize)
*  0 48 128 128 128                 Values of NX (crossover point)
*  2                                Number of values of LDA
*  512 513                          Values of LDA (leading dimension)
*  0.0                              Minimum time in seconds
*  SQR    T T F
*  SQP    T
*  SRR    T
*
*  The routines are timed for all combinations of applicable values of
*  M, N, K, NB, NX, and LDA, and for all combinations of options such as
*  UPLO and TRANS.  For Level 2 BLAS timings, values of NB are used for
*  INCX.  Certain subroutines, such as the QR factorization, treat the
*  values of M and N as ordered pairs and operate on M x N matrices.
*
*  Internal Parameters
*  ===================
*
*  NMAX    INTEGER
*          The maximum value of M or N for square matrices.
*
*  LDAMAX  INTEGER
*          The maximum value of LDA.
*
*  NMAXB   INTEGER
*          The maximum value of N for band matrices.
*
*  MAXVAL  INTEGER
*          The maximum number of values that can be read in for M, N,
*          K, NB, or NX.
*
*  MXNLDA  INTEGER
*          The maximum number of values that can be read in for LDA.
*
*  NIN     INTEGER
*          The unit number for input.  Currently set to 5 (std input).
*
*  NOUT    INTEGER
*          The unit number for output.  Currently set to 6 (std output).
*
*  =====================================================================
*
*     .. Parameters ..
      INTEGER            NMAX, LDAMAX, NMAXB
      PARAMETER          ( NMAX = 1001, LDAMAX = NMAX+4, NMAXB = 5000 )
      INTEGER            LA
      PARAMETER          ( LA = NMAX*LDAMAX )
      INTEGER            MAXVAL, MXNLDA
      PARAMETER          ( MAXVAL = 12, MXNLDA = 4 )
      INTEGER            NIN, NOUT
      PARAMETER          ( NIN = 5, NOUT = 6 )
*     ..
*     .. Local Scalars ..
      LOGICAL            BLAS, LDAMOK, LDANOK, LDAOK, MOK, NOK, NXNBOK
      CHARACTER          C1
      CHARACTER*2        C2
      CHARACTER*3        C3
      CHARACTER*80       LINE
      INTEGER            I, L, LDR1, LDR2, LDR3, MAXK, MAXLDA,
     $                   MAXM, MAXN, MAXNB, MKMAX, NEED, NK, NLDA, NM,
     $                   NN, NNB
      REAL               S1, S2, TIMMIN
*     ..
*     .. Local Arrays ..
      INTEGER            IWORK( 2*NMAXB ), KVAL( MAXVAL ),
     $                   LDAVAL( MXNLDA ), MVAL( MAXVAL ),
     $                   NBVAL( MAXVAL ), NVAL( MAXVAL ),
     $                   NXVAL( MAXVAL )
      REAL               A( LA, 4 ), D( 2*NMAX, 2 ),
     $                   RESLTS( MAXVAL, MAXVAL, 2*MXNLDA, 4*MAXVAL )
*     ..
*     .. External Functions ..
      LOGICAL            LSAME, LSAMEN
      REAL               SECOND
      EXTERNAL           LSAME, LSAMEN, SECOND
*     ..
*     .. External Subroutines ..
      EXTERNAL           STIMMM, STIMMV, STIMQP, STIMQR,
     $                   STIMRR
*     ..
*     .. Scalars in Common ..
      INTEGER            NB, NEISPK, NPROC, NSHIFT
*     ..
*     .. Common blocks ..
      COMMON             / CENVIR / NB, NPROC, NSHIFT, NEISPK
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX
*     ..
*     .. Executable Statements ..
*
      S1 = SECOND( )
      LDR1 = MAXVAL
      LDR2 = MAXVAL
      LDR3 = 2*MXNLDA
*
*     Read the first line.  The first four characters must be 'BLAS'
*     for the BLAS data file format to be used.  Otherwise, the LAPACK
*     data file format is assumed.
*
      READ( NIN, FMT = '( A80 )' )LINE
      BLAS = LSAMEN( 4, LINE, 'BLAS' )
*
*     Find the last non-blank and print the first line of input as the
*     first line of output.
*
      DO 10 L = 80, 1, -1
         IF( LINE( L: L ).NE.' ' )
     $      GO TO 20
   10 CONTINUE
      L = 1
   20 CONTINUE
      WRITE( NOUT, FMT = '( 1X, A, / )' )LINE( 1: L )
      WRITE( NOUT, FMT = 9992 )
*
*     Read in NM and the values for M.
*
      READ( NIN, FMT = * )NM
      IF( NM.GT.MAXVAL ) THEN
         WRITE( NOUT, FMT = 9999 )'M', 'NM', MAXVAL
         NM = MAXVAL
      END IF
      READ( NIN, FMT = * )( MVAL( I ), I = 1, NM )
      WRITE( NOUT, FMT = 9991 )'M:     ', ( MVAL( I ), I = 1, NM )
*
*     Check that  M <= NMAXB for all values of M.
*
      MOK = .TRUE.
      MAXM = 0
      DO 30 I = 1, NM
         MAXM = MAX( MVAL( I ), MAXM )
         IF( MVAL( I ).GT.NMAXB ) THEN
            WRITE( NOUT, FMT = 9997 )'M', MVAL( I ), NMAXB
            MOK = .FALSE.
         END IF
   30 CONTINUE
      IF( .NOT.MOK )
     $   WRITE( NOUT, FMT = * )
*
*     Read in NN and the values for N.
*
      READ( NIN, FMT = * )NN
      IF( NN.GT.MAXVAL ) THEN
         WRITE( NOUT, FMT = 9999 )'N', 'NN', MAXVAL
         NN = MAXVAL
      END IF
      READ( NIN, FMT = * )( NVAL( I ), I = 1, NN )
      WRITE( NOUT, FMT = 9991 )'N:     ', ( NVAL( I ), I = 1, NN )
*
*     Check that  N <= NMAXB for all values of N.
*
      NOK = .TRUE.
      MAXN = 0
      DO 40 I = 1, NN
         MAXN = MAX( NVAL( I ), MAXN )
         IF( NVAL( I ).GT.NMAXB ) THEN
            WRITE( NOUT, FMT = 9997 )'N', NVAL( I ), NMAXB
            NOK = .FALSE.
         END IF
   40 CONTINUE
      IF( .NOT.NOK )
     $   WRITE( NOUT, FMT = * )
*
*     Read in NK and the values for K.
*
      READ( NIN, FMT = * )NK
      IF( NK.GT.MAXVAL ) THEN
         WRITE( NOUT, FMT = 9999 )'K', 'NK', MAXVAL
         NK = MAXVAL
      END IF
      READ( NIN, FMT = * )( KVAL( I ), I = 1, NK )
      WRITE( NOUT, FMT = 9991 )'K:     ', ( KVAL( I ), I = 1, NK )
*
*     Find the maximum value of K (= NRHS).
*
      MAXK = 0
      DO 50 I = 1, NK
         MAXK = MAX( KVAL( I ), MAXK )
   50 CONTINUE
      MKMAX = MAXM*MAX( 2, MAXK )
*
*     Read in NNB and the values for NB.  For the BLAS input files,
*     NBVAL is used to store values for INCX and INCY.
*
      READ( NIN, FMT = * )NNB
      IF( NNB.GT.MAXVAL ) THEN
         WRITE( NOUT, FMT = 9999 )'NB', 'NNB', MAXVAL
         NNB = MAXVAL
      END IF
      READ( NIN, FMT = * )( NBVAL( I ), I = 1, NNB )
*
*     Find the maximum value of NB.
*
      MAXNB = 0
      DO 60 I = 1, NNB
         MAXNB = MAX( NBVAL( I ), MAXNB )
   60 CONTINUE
*
      IF( BLAS ) THEN
         WRITE( NOUT, FMT = 9991 )'INCX:  ', ( NBVAL( I ), I = 1, NNB )
         DO 70 I = 1, NNB
            NXVAL( I ) = 0
   70    CONTINUE
      ELSE
*
*        LAPACK data files:  Read in the values for NX.
*
         READ( NIN, FMT = * )( NXVAL( I ), I = 1, NNB )
*
         WRITE( NOUT, FMT = 9991 )'NB:    ', ( NBVAL( I ), I = 1, NNB )
         WRITE( NOUT, FMT = 9991 )'NX:    ', ( NXVAL( I ), I = 1, NNB )
      END IF
*
*     Read in NLDA and the values for LDA.
*
      READ( NIN, FMT = * )NLDA
      IF( NLDA.GT.MXNLDA ) THEN
         WRITE( NOUT, FMT = 9999 )'LDA', 'NLDA', MXNLDA
         NLDA = MXNLDA
      END IF
      READ( NIN, FMT = * )( LDAVAL( I ), I = 1, NLDA )
      WRITE( NOUT, FMT = 9991 )'LDA:   ', ( LDAVAL( I ), I = 1, NLDA )
*
*     Check that LDA >= 1 for all values of LDA.
*
      LDAOK = .TRUE.
      MAXLDA = 0
      DO 80 I = 1, NLDA
         MAXLDA = MAX( LDAVAL( I ), MAXLDA )
         IF( LDAVAL( I ).LE.0 ) THEN
            WRITE( NOUT, FMT = 9998 )LDAVAL( I )
            LDAOK = .FALSE.
         END IF
   80 CONTINUE
      IF( .NOT.LDAOK )
     $   WRITE( NOUT, FMT = * )
*
*     Check that MAXLDA*MAXN <= LA (for the dense routines).
*
      LDANOK = .TRUE.
      NEED = MAXLDA*MAXN
      IF( NEED.GT.LA ) THEN
         WRITE( NOUT, FMT = 9995 )MAXLDA, MAXN, NEED
         LDANOK = .FALSE.
      END IF
*
*     Check that MAXLDA*MAXM + MAXM*MAXK <= 3*LA (for band routines).
*
      LDAMOK = .TRUE.
      NEED = MAXLDA*MAXM + MAXM*MAXK
      IF( NEED.GT.3*LA ) THEN
         NEED = ( NEED+2 ) / 3
         WRITE( NOUT, FMT = 9994 )MAXLDA, MAXM, MAXK, NEED
         LDAMOK = .FALSE.
      END IF
*
*     Check that MAXN*MAXNB (or MAXN*INCX) <= LA.
*
      NXNBOK = .TRUE.
      NEED = MAXN*MAXNB
      IF( NEED.GT.LA ) THEN
         WRITE( NOUT, FMT = 9996 )MAXN, MAXNB, NEED
         NXNBOK = .FALSE.
      END IF
*
      IF( .NOT.( MOK .AND. NOK .AND. LDAOK .AND. LDANOK .AND. NXNBOK ) )
     $     THEN
         WRITE( NOUT, FMT = 9984 )
         GO TO 110
      END IF
      IF( .NOT.LDAMOK )
     $   WRITE( NOUT, FMT = * )
*
*     Read the minimum time to time a subroutine.
*
      WRITE( NOUT, FMT = * )
      READ( NIN, FMT = * )TIMMIN
      WRITE( NOUT, FMT = 9993 )TIMMIN
      WRITE( NOUT, FMT = * )
*
*     Read the first input line.
*
      READ( NIN, FMT = '(A)', END = 100 )LINE
*
*     If the first record is the special signal 'NONE', then get the
*     next line but don't time SGEMV.
*
      IF( LSAMEN( 4, LINE, 'NONE' ) ) THEN
         READ( NIN, FMT = '(A)', END = 100 )LINE
      ELSE
         WRITE( NOUT, FMT = 9990 )
*
*        Time SGEMV and SGEMM.
*
         CALL STIMMV( 'SGEMV ', NN, NVAL, NNB, NBVAL, NLDA,
     $                LDAVAL, TIMMIN, A( 1, 1 ), LA, A( 1, 2 ),
     $                A( 1, 3 ), RESLTS, LDR1, LDR2, NOUT )
         CALL STIMMM( 'SGEMM ', 'N', NN, NVAL, NLDA, LDAVAL,
     $                TIMMIN, A( 1, 1 ), A( 1, 2 ), A( 1, 3 ),
     $                RESLTS, LDR1, LDR2, NOUT )
      END IF
*
*     Call the appropriate timing routine for each input line.
*
      WRITE( NOUT, FMT = 9988 )
   90 CONTINUE
      C1 = LINE( 1: 1 )
      C2 = LINE( 2: 3 )
      C3 = LINE( 4: 6 )
*
*     Check first character for correct precision.
*
      IF( .NOT.LSAME( C1, 'Sprecision' ) ) THEN
         WRITE( NOUT, FMT = 9987 )LINE( 1: 6 )
*
      ELSE IF( LSAMEN( 2, C2, 'QR' ) .OR. LSAMEN( 2, C3, 'QR' ) .OR.
     $         LSAMEN( 2, C3( 2: 3 ), 'QR' ) ) THEN
*
*        QR routines
*
         CALL STIMQR( LINE, NN, MVAL, NVAL, NK, KVAL, NNB, NBVAL,
     $                NXVAL, NLDA, LDAVAL, TIMMIN, A( 1, 1 ), D,
     $                A( 1, 2 ), A( 1, 3 ), RESLTS, LDR1, LDR2,
     $                LDR3, NOUT )
*
      ELSE IF( LSAMEN( 2, C2, 'QP' ) .OR. LSAMEN( 3, C3, 'QPF' ) ) THEN
*
*        QR with column pivoting
*
         CALL STIMQP( LINE, NM, MVAL, NVAL, NLDA, LDAVAL, TIMMIN,
     $                A( 1, 1 ), A( 1, 2 ), D( 1, 1 ), A( 1, 3 ), IWORK,
     $                RESLTS, LDR1, LDR2, NOUT )
*
      ELSE IF( LSAMEN( 2, C2, 'RR' ) .OR. LSAMEN( 3, C3, 'RRF' ) ) THEN
*
*        Rank-Revealing QR
*
         CALL STIMRR( LINE, NM, MVAL, NVAL, NK, KVAL,
     $                NNB, NBVAL, NXVAL, NLDA, LDAVAL, TIMMIN,
     $                A( 1, 1 ), A( 1, 2 ), A( 1, 3 ), A( 1, 4 ),
     $                IWORK, RESLTS, LDR1, LDR2, NOUT )
*
      ELSE
*
         WRITE( NOUT, FMT = 9987 )LINE( 1: 6 )
      END IF
*
*     Read the next line of the input file.
*
      READ( NIN, FMT = '(A)', END = 100 )LINE
      GO TO 90
*
*     Branch to this line when the last record is read.
*
  100 CONTINUE
      S2 = SECOND( )
      WRITE( NOUT, FMT = 9986 )
      WRITE( NOUT, FMT = 9985 )S2 - S1
  110 CONTINUE
*
 9999 FORMAT( ' Too many values of ', A, ' using ', A, ' = ', I2 )
 9998 FORMAT( ' *** LDA = ', I7, ' is too small, must have ',
     $      'LDA > 0.' )
 9997 FORMAT( ' *** ', A1, ' = ', I7, ' is too big:  ',
     $      'maximum allowed is', I7 )
 9996 FORMAT( ' *** N*NB is too big for N =', I6, ', NB =', I6,
     $      / ' --> Increase LA to at least ', I8 )
 9995 FORMAT( ' *** LDA*N is too big for the dense routines ', '(LDA =',
     $      I6, ', N =', I6, ')', / ' --> Increase LA to at least ',
     $      I8 )
 9994 FORMAT( ' *** (LDA+K)*M is too big for the band routines ',
     $      '(LDA=', I6, ', M=', I6, ', K=', I6, ')',
     $      / ' --> Increase LA to at least ', I8 )
 9993 FORMAT( ' The minimum time a subroutine will be timed = ', F6.3,
     $      ' seconds' )
 9992 FORMAT( ' The following parameter values will be used:' )
 9991 FORMAT( 4X, A7, 1X, 10I6, / 12X, 10I6 )
 9990 FORMAT( / ' ------------------------------',
     $      / ' >>>>>    Sample BLAS     <<<<<',
     $      / ' ------------------------------' )
 9988 FORMAT( / ' ------------------------------',
     $      / ' >>>>>    Timing data     <<<<<',
     $      / ' ------------------------------' )
 9987 FORMAT( 1X, A6, ':  Unrecognized path or subroutine name', / )
 9986 FORMAT( ' End of tests' )
 9985 FORMAT( ' Total time used = ', F12.2, ' seconds' )
 9984 FORMAT( / ' Tests not done due to input errors' )
*
*     End of STIMAA
*
      END
