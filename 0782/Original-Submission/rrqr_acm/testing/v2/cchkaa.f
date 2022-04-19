      PROGRAM CCHKAA
*
*  -- LAPACK test routine (version 2.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     September 30, 1994
*
*  Purpose
*  =======
*
*  CCHKAA is the main test program for the COMPLEX linear
*  equation routines.
*
*  An annotated example of a data file can be obtained by deleting the
*  first 3 characters from the following lines:
*  Data file for testing COMPLEX LAPACK linear equation routines
*  7                      Number of values of M
*  0 1 2 3 5 10 16        Values of M (row dimension)
*  7                      Number of values of N
*  0 1 2 3 5 10 16        Values of N (column dimension)
*  5                      Number of values of NB
*  1 3 3 3 20             Values of NB (the blocksize)
*  1 0 5 9 1              Values of NX (crossover point)
*  30.0                   Threshold value of test ratio
*  T                      Put T to test the LAPACK routines
*  T                      Put T to test the driver routines
*  T                      Put T to test the error exits
*  CQP    6               List types on next line if 0 < NTYPES <  6
*  CRR    3               List types on next line if 0 < NTYPES <  3
*
*  Internal Parameters
*  ===================
*
*  NMAX    INTEGER
*          The maximum allowable value for N.
*
*  MAXIN   INTEGER
*          The number of different values that can be used for each of
*          M, N, or NB
*
*  MAXRHS  INTEGER
*          The maximum number of right hand sides
*
*  NIN     INTEGER
*          The unit number for input
*
*  NOUT    INTEGER
*          The unit number for output
*
*  LWORK   INTEGER
*          The workspace length
*
*  =====================================================================
*
*     .. Parameters ..
      INTEGER            NMAX
      PARAMETER          ( NMAX = 501 )
      INTEGER            MAXIN
      PARAMETER          ( MAXIN = 12 )
      INTEGER            MAXRHS
      PARAMETER          ( MAXRHS = 10 )
      INTEGER            MATMAX
      PARAMETER          ( MATMAX = 30 )
      INTEGER            NIN, NOUT
      PARAMETER          ( NIN = 5, NOUT = 6 )
      INTEGER            KDMAX
      PARAMETER          ( KDMAX = NMAX+( NMAX+1 ) / 4 )
*
      INTEGER            LWORK
      PARAMETER          ( LWORK = NMAX*( NMAX+MAXRHS ) )
*     ..
*     .. Local Scalars ..
      LOGICAL            FATAL, TSTCHK, TSTERR
      CHARACTER          C1
      CHARACTER*2        C2
      CHARACTER*3        PATH
      CHARACTER*10       INTSTR
      CHARACTER*72       ALINE
      INTEGER            I, IC, J, K, LDA, NB, NM, NMATS, NN,
     $                   NNB, NNB2, NTYPES
      REAL               EPS, S1, S2, THREQ, THRESH
*     ..
*     .. Local Arrays ..
      LOGICAL            DOTYPE( MATMAX )
      INTEGER            IWORK( 3*NMAX ), MVAL( MAXIN ), NBVAL( MAXIN ),
     $                   NBVAL2( MAXIN ), NVAL( MAXIN ), NXVAL( MAXIN )
      REAL               RWORK( 5*NMAX+2*MAXRHS ), S( 2*NMAX )
      COMPLEX            A( ( KDMAX+1 )*NMAX, 7 ), B( NMAX*MAXRHS, 4 ),
     $                   WORK( NMAX, NMAX+MAXRHS )
*     ..
*     .. External Functions ..
      LOGICAL            LSAME, LSAMEN
      REAL               SLAMCH, SECOND
      EXTERNAL           LSAME, LSAMEN, SLAMCH, SECOND
*     ..
*     .. External Subroutines ..
      EXTERNAL           ALAREQ, CCHKQP, CCHKRR
*     ..
*     .. Scalars in Common ..
      LOGICAL            LERR, OK
      CHARACTER*6        SRNAMT
      INTEGER            INFOT, NUNIT
*     ..
*     .. Arrays in Common ..
      INTEGER            IPARMS( 100 )
*     ..
*     .. Common blocks ..
      COMMON             / INFOC / INFOT, NUNIT, OK, LERR
      COMMON             / SRNAMC / SRNAMT
      COMMON             / CLAENV / IPARMS
*     ..
*     .. Data statements ..
      DATA               THREQ / 2.0D0 / , INTSTR / '0123456789' /
*     ..
*     .. Executable Statements ..
*
      S1 = SECOND( )
      LDA = NMAX
      FATAL = .FALSE.
*
*     Read a dummy line.
*
      READ( NIN, FMT = * )
*
*     Report values of parameters.
*
      WRITE( NOUT, FMT = 9993 )
*
*     Read the values of M
*
      READ( NIN, FMT = * )NM
      IF( NM.LT.1 ) THEN
         WRITE( NOUT, FMT = 9996 )' NM ', NM, 1
         NM = 0
         FATAL = .TRUE.
      ELSE IF( NM.GT.MAXIN ) THEN
         WRITE( NOUT, FMT = 9995 )' NM ', NM, MAXIN
         NM = 0
         FATAL = .TRUE.
      END IF
      READ( NIN, FMT = * )( MVAL( I ), I = 1, NM )
      DO 10 I = 1, NM
         IF( MVAL( I ).LT.0 ) THEN
            WRITE( NOUT, FMT = 9996 )' M  ', MVAL( I ), 0
            FATAL = .TRUE.
         ELSE IF( MVAL( I ).GT.NMAX ) THEN
            WRITE( NOUT, FMT = 9995 )' M  ', MVAL( I ), NMAX
            FATAL = .TRUE.
         END IF
   10 CONTINUE
      IF( NM.GT.0 )
     $   WRITE( NOUT, FMT = 9992 )'M   ', ( MVAL( I ), I = 1, NM )
*
*     Read the values of N
*
      READ( NIN, FMT = * )NN
      IF( NN.LT.1 ) THEN
         WRITE( NOUT, FMT = 9996 )' NN ', NN, 1
         NN = 0
         FATAL = .TRUE.
      ELSE IF( NN.GT.MAXIN ) THEN
         WRITE( NOUT, FMT = 9995 )' NN ', NN, MAXIN
         NN = 0
         FATAL = .TRUE.
      END IF
      READ( NIN, FMT = * )( NVAL( I ), I = 1, NN )
      DO 20 I = 1, NN
         IF( NVAL( I ).LT.0 ) THEN
            WRITE( NOUT, FMT = 9996 )' N  ', NVAL( I ), 0
            FATAL = .TRUE.
         ELSE IF( NVAL( I ).GT.NMAX ) THEN
            WRITE( NOUT, FMT = 9995 )' N  ', NVAL( I ), NMAX
            FATAL = .TRUE.
         END IF
   20 CONTINUE
      IF( NN.GT.0 )
     $   WRITE( NOUT, FMT = 9992 )'N   ', ( NVAL( I ), I = 1, NN )
*
*     Read the values of NB
*
      READ( NIN, FMT = * )NNB
      IF( NNB.LT.1 ) THEN
         WRITE( NOUT, FMT = 9996 )'NNB ', NNB, 1
         NNB = 0
         FATAL = .TRUE.
      ELSE IF( NNB.GT.MAXIN ) THEN
         WRITE( NOUT, FMT = 9995 )'NNB ', NNB, MAXIN
         NNB = 0
         FATAL = .TRUE.
      END IF
      READ( NIN, FMT = * )( NBVAL( I ), I = 1, NNB )
      DO 40 I = 1, NNB
         IF( NBVAL( I ).LT.0 ) THEN
            WRITE( NOUT, FMT = 9996 )' NB ', NBVAL( I ), 0
            FATAL = .TRUE.
         END IF
   40 CONTINUE
      IF( NNB.GT.0 )
     $   WRITE( NOUT, FMT = 9992 )'NB  ', ( NBVAL( I ), I = 1, NNB )
*
*     Set NBVAL2 to be the set of unique values of NB
*
      NNB2 = 0
      DO 60 I = 1, NNB
         NB = NBVAL( I )
         DO 50 J = 1, NNB2
            IF( NB.EQ.NBVAL2( J ) )
     $         GO TO 60
   50    CONTINUE
         NNB2 = NNB2 + 1
         NBVAL2( NNB2 ) = NB
   60 CONTINUE
*
*     Read the values of NX
*
      READ( NIN, FMT = * )( NXVAL( I ), I = 1, NNB )
      DO 70 I = 1, NNB
         IF( NXVAL( I ).LT.0 ) THEN
            WRITE( NOUT, FMT = 9996 )' NX ', NXVAL( I ), 0
            FATAL = .TRUE.
         END IF
   70 CONTINUE
      IF( NNB.GT.0 )
     $   WRITE( NOUT, FMT = 9992 )'NX  ', ( NXVAL( I ), I = 1, NNB )
*
*     Read the threshold value for the test ratios.
*
      READ( NIN, FMT = * )THRESH
      WRITE( NOUT, FMT = 9991 )THRESH
*
*     Read the flag that indicates whether to test the LAPACK routines.
*
      READ( NIN, FMT = * )TSTCHK
*
*     Read the flag that indicates whether to test the error exits.
*
      READ( NIN, FMT = * )TSTERR
*
      IF( FATAL ) THEN
         WRITE( NOUT, FMT = 9999 )
         STOP
      END IF
*
*     Calculate and print the machine dependent constants.
*
      EPS = SLAMCH( 'Underflow threshold' )
      WRITE( NOUT, FMT = 9990 )'underflow', EPS
      EPS = SLAMCH( 'Overflow threshold' )
      WRITE( NOUT, FMT = 9990 )'overflow ', EPS
      EPS = SLAMCH( 'Epsilon' )
      WRITE( NOUT, FMT = 9990 )'precision', EPS
      WRITE( NOUT, FMT = * )
*
   80 CONTINUE
*
*     Read a test path and the number of matrix types to use.
*
      READ( NIN, FMT = '(A72)', END = 140 )ALINE
      PATH = ALINE( 1: 3 )
      NMATS = MATMAX
      I = 3
   90 CONTINUE
      I = I + 1
      IF( I.GT.72 )
     $   GO TO 130
      IF( ALINE( I: I ).EQ.' ' )
     $   GO TO 90
      NMATS = 0
  100 CONTINUE
      C1 = ALINE( I: I )
      DO 110 K = 1, 10
         IF( C1.EQ.INTSTR( K: K ) ) THEN
            IC = K - 1
            GO TO 120
         END IF
  110 CONTINUE
      GO TO 130
  120 CONTINUE
      NMATS = NMATS*10 + IC
      I = I + 1
      IF( I.GT.72 )
     $   GO TO 130
      GO TO 100
  130 CONTINUE
      C1 = PATH( 1: 1 )
      C2 = PATH( 2: 3 )
*
*     Check first character for correct precision.
*
      IF( .NOT.LSAME( C1, 'C precision' ) ) THEN
         WRITE( NOUT, FMT = 9989 )PATH
*
      ELSE IF( NMATS.LE.0 ) THEN
*
*        Check for a positive number of tests requested.
*
         WRITE( NOUT, FMT = 9988 )PATH
*
      ELSE IF( LSAMEN( 2, C2, 'QP' ) ) THEN
*
*        QP:  QR factorization with pivoting
*
         NTYPES = 6
         CALL ALAREQ( PATH, NMATS, DOTYPE, NTYPES, NIN, NOUT )
*
         IF( TSTCHK ) THEN
            CALL CCHKQP( DOTYPE, NM, MVAL, NN, NVAL, THRESH,
     $                   TSTERR, A( 1, 1 ), A( 1, 2 ), S( 1 ),
     $                   S( NMAX+1 ), B( 1, 1 ), WORK, RWORK, IWORK,
     $                   NOUT )
         ELSE
            WRITE( NOUT, FMT = 9988 )PATH
         END IF
*
      ELSE IF( LSAMEN( 2, C2, 'RR' ) ) THEN
*
*        RR:  Rank-Revealing QR factorization
*
         NTYPES = 3
         CALL ALAREQ( PATH, NMATS, DOTYPE, NTYPES, NIN, NOUT )
*
         IF( TSTCHK ) THEN
            CALL CCHKRR( DOTYPE, NM, MVAL, NN, NVAL, NNB, NBVAL,
     $                   NXVAL, THRESH, TSTERR, A( 1, 1 ), A( 1, 2 ),
     $                   A( 1, 3 ), S( 1 ), S( NMAX+1 ),
     $                   WORK, LWORK, RWORK, IWORK, NOUT )
         ELSE
            WRITE( NOUT, FMT = 9988 )PATH
         END IF
*
      ELSE
*
         WRITE( NOUT, FMT = 9989 )PATH
      END IF
*
*     Go back to get another input line.
*
      GO TO 80
*
*     Branch to this line when the last record is read.
*
  140 CONTINUE
      CLOSE ( NIN )
      S2 = SECOND( )
      WRITE( NOUT, FMT = 9998 )
      WRITE( NOUT, FMT = 9997 )S2 - S1
*
 9999 FORMAT( / ' Execution not attempted due to input errors' )
 9998 FORMAT( / ' End of tests' )
 9997 FORMAT( ' Total time used = ', F12.2, ' seconds', / )
 9996 FORMAT( ' Invalid input value: ', A4, '=', I6, '; must be >=',
     $      I6 )
 9995 FORMAT( ' Invalid input value: ', A4, '=', I6, '; must be <=',
     $      I6 )
 9993 FORMAT( ' RRQR Tests. COMPLEX LAPACK routines ',
     $      / / ' The following parameter values will be used:' )
 9992 FORMAT( 4X, A4, ':  ', 10I6, / 11X, 10I6 )
 9991 FORMAT( / ' Routines pass computational tests if test ratio is ',
     $      'less than', F8.2, / )
 9990 FORMAT( ' Relative machine ', A, ' is taken to be', D16.6 )
 9989 FORMAT( / 1X, A3, ':  Unrecognized path name' )
 9988 FORMAT( / 1X, A3, ' routines were not tested' )
*
*     End of CCHKAA
*
      END
