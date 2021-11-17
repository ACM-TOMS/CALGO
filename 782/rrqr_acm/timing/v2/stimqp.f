      SUBROUTINE STIMQP( LINE, NM, MVAL, NVAL, NLDA, LDAVAL, TIMMIN, A,
     $                   COPYA, TAU, WORK, IWORK, RESLTS, LDR1, LDR2,
     $                   NOUT )
*
*  -- LAPACK timing routine (version 2.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     March 31, 1993
*
*     .. Scalar Arguments ..
      CHARACTER*80       LINE
      INTEGER            LDR1, LDR2, NLDA, NM, NOUT
      REAL               TIMMIN
*     ..
*     .. Array Arguments ..
      INTEGER            IWORK( * ), LDAVAL( * ), MVAL( * ), NVAL( * )
      REAL               A( * ), COPYA( * ), RESLTS( LDR1, LDR2, * ),
     $                   TAU( * ), WORK( * )
*     ..
*
*  Purpose
*  =======
*
*  STIMQP times the LAPACK routines to perform the QR factorization with
*  column pivoting of a REAL general matrix.
*
*  Two matrix types may be used for timing.  The number of types is
*  set in the parameter NMODE and the matrix types are set in the vector
*  MODES, using the following key:
*     2.  BREAK1    D(1:N-1)=1 and D(N)=1.0/COND in SLATMS
*     3.  GEOM      D(I)=COND**(-(I-1)/(N-1)) in SLATMS
*  These numbers are chosen to correspond with the matrix types in the
*  test code.
*
*  Arguments
*  =========
*
*  LINE    (input) CHARACTER*80
*          The input line that requested this routine.  The first six
*          characters contain either the name of a subroutine or a
*          generic path name.  The remaining characters may be used to
*          specify the individual routines to be timed.  See ATIMIN for
*          a full description of the format of the input line.
*
*  NM      (input) INTEGER
*          The number of values of M and N contained in the vectors
*          MVAL and NVAL.  The matrix sizes are used in pairs (M,N).
*
*  MVAL    (input) INTEGER array, dimension (NM)
*          The values of the matrix row dimension M.
*
*  NVAL    (input) INTEGER array, dimension (NM)
*          The values of the matrix column dimension N.
*
*  NLDA    (input) INTEGER
*          The number of values of LDA contained in the vector LDAVAL.
*
*  LDAVAL  (input) INTEGER array, dimension (NLDA)
*          The values of the leading dimension of the array A.
*
*  TIMMIN  (input) REAL
*          The minimum time a subroutine will be timed.
*
*  A       (workspace) REAL array, dimension (LDAMAX*NMAX)
*          where LDAMAX and NMAX are the maximum values of LDA and N.
*
*  COPYA   (workspace) REAL array, dimension (LDAMAX*NMAX)
*
*  TAU     (workspace) REAL array, dimension (min(M,N))
*
*  WORK    (workspace) REAL array, dimension (3*NMAX)
*
*  IWORK   (workspace) INTEGER array, dimension (2*NMAX)
*
*  RESLTS  (workspace) REAL array, dimension
*                      (LDR1,LDR2,NLDA)
*          The timing results for each subroutine over the relevant
*          values of MODE, (M,N), and LDA.
*
*  LDR1    (input) INTEGER
*          The first dimension of RESLTS.  LDR1 >= max(1,NM).
*
*  LDR2    (input) INTEGER
*          The second dimension of RESLTS.  LDR2 >= max(1,NM).
*
*  NOUT    (input) INTEGER
*          The unit number for output.
*
*  =====================================================================
*
*     .. Parameters ..
      INTEGER            NSUBS, NMODE
      PARAMETER          ( NSUBS = 1, NMODE = 2 )
      REAL               ONE
      PARAMETER          ( ONE = 1.0E+0 )
*     ..
*     .. Local Scalars ..
      CHARACTER*3        PATH
      CHARACTER*6        CNAME
      INTEGER            I, IC, ICL, ILDA, IM, IMODE, INFO, LDA, M,
     $                   MINMN, MODE, N, NB
      REAL               COND, DMAX, OPS, S1, S2, TIME, UNTIME
*     ..
*     .. Local Arrays ..
      LOGICAL            TIMSUB( NSUBS )
      CHARACTER*6        SUBNAM( NSUBS )
      INTEGER            ISEED( 4 ), MODES( NMODE )
*     ..
*     .. External Functions ..
      REAL               SECOND, SLAMCH, SMFLOP, SOPLA
      EXTERNAL           SECOND, SLAMCH, SMFLOP, SOPLA
*     ..
*     .. External Subroutines ..
      EXTERNAL           ATIMCK, ATIMIN, ICOPY, SGEQPF, SLACPY, SLATMS,
     $                   SPRTB5
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MIN, REAL
*     ..
*     .. Data statements ..
      DATA               SUBNAM / 'SGEQPF' /
      DATA               MODES / 2, 3 /
      DATA               ISEED / 0, 0, 0, 1 /
*     ..
*     .. Executable Statements ..
*
*     Extract the timing request from the input line.
*
      PATH( 1: 1 ) = 'Single precision'
      PATH( 2: 3 ) = 'QP'
      CALL ATIMIN( PATH, LINE, NSUBS, SUBNAM, TIMSUB, NOUT, INFO )
      IF( .NOT.TIMSUB( 1 ) .OR. INFO.NE.0 )
     $   GO TO 90
*
*     Check that M <= LDA for the input values.
*
      CNAME = LINE( 1: 6 )
      CALL ATIMCK( 1, CNAME, NM, MVAL, NLDA, LDAVAL, NOUT, INFO )
      IF( INFO.GT.0 ) THEN
         WRITE( NOUT, FMT = 9999 )CNAME
         GO TO 90
      END IF
*
*     Set the condition number and scaling factor for the matrices
*     to be generated.
*
      DMAX = ONE
      COND = ONE / SLAMCH( 'Precision' )
*
*     Do for each pair of values (M,N):
*
      DO 60 IM = 1, NM
         M = MVAL( IM )
         N = NVAL( IM )
         MINMN = MIN( M, N )
*
*        Do for each value of LDA:
*
         DO 50 ILDA = 1, NLDA
            LDA = LDAVAL( ILDA )
            DO 40 IMODE = 1, NMODE
               MODE = MODES( IMODE )
*
*              Generate a test matrix of size m by n using the
*              singular value distribution indicated by MODE.
*
               DO 10 I = 1, N
                  IWORK( N+I ) = 0
   10          CONTINUE
               CALL SLATMS( M, N, 'Uniform', ISEED, 'Nonsymm', TAU,
     $                      MODE, COND, DMAX, M, N, 'No packing', COPYA,
     $                      LDA, WORK, INFO )
*
*              SGEQPF:  QR factorization with column pivoting
*
               CALL SLACPY( 'All', M, N, COPYA, LDA, A, LDA )
               CALL ICOPY( N, IWORK( N+1 ), 1, IWORK, 1 )
               IC = 0
               S1 = SECOND( )
   20          CONTINUE
               CALL SGEQPF( M, N, A, LDA, IWORK, TAU, WORK, INFO )
               S2 = SECOND( )
               TIME = S2 - S1
               IC = IC + 1
               IF( TIME.LT.TIMMIN ) THEN
                  CALL SLACPY( 'All', M, N, COPYA, LDA, A, LDA )
                  CALL ICOPY( N, IWORK( N+1 ), 1, IWORK, 1 )
                  GO TO 20
               END IF
*
*              Subtract the time used in SLACPY and ICOPY.
*
               ICL = 1
               S1 = SECOND( )
   30          CONTINUE
               S2 = SECOND( )
               UNTIME = S2 - S1
               ICL = ICL + 1
               IF( ICL.LE.IC ) THEN
                  CALL SLACPY( 'All', M, N, COPYA, LDA, A, LDA )
                  CALL ICOPY( N, IWORK( N+1 ), 1, IWORK, 1 )
                  GO TO 30
               END IF
*
               TIME = ( TIME-UNTIME ) / REAL( IC )
               OPS = SOPLA( 'SGEQPF', M, N, 0, 0, NB )
               RESLTS( IMODE, IM, ILDA ) = SMFLOP( OPS, TIME, INFO )
*
   40       CONTINUE
   50    CONTINUE
   60 CONTINUE
*
*     Print tables of results
*
      WRITE( NOUT, FMT = 9998 )SUBNAM( 1 )
      IF( NLDA.GT.1 ) THEN
         DO 70 I = 1, NLDA
            WRITE( NOUT, FMT = 9997 )I, LDAVAL( I )
   70    CONTINUE
      END IF
      WRITE( NOUT, FMT = * )
      CALL SPRTB5( 'Type', 'M', 'N', NMODE, MODES, NM, MVAL, NVAL, NLDA,
     $             RESLTS, LDR1, LDR2, NOUT )
   90 CONTINUE
 9999 FORMAT( 1X, A6, ' timing run not attempted', / )
 9998 FORMAT( / ' *** Speed of ', A6, ' in megaflops ***' )
 9997 FORMAT( 5X, 'line ', I2, ' with LDA = ', I5 )
      RETURN
*
*     End of STIMQP
*
      END
