      SUBROUTINE ZTIMRR( LINE, NM, MVAL, NVAL, NK, KVAL,
     $                   NNB, NBVAL, NXVAL, NLDA, LDAVAL, TIMMIN,
     $                   A, COPYA, B, WORK, RWORK, IWORK, RESLTS,
     $                   LDR1, LDR2, NOUT )
*
*  -- LAPACK timing routine (version 2.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     March 31, 1993
*
*     Rewritten for timing rrqr code.
*
*     .. Scalar Arguments ..
      CHARACTER*80       LINE
      INTEGER            LDR1, LDR2, NK, NLDA, NM, NNB, NOUT
      DOUBLE PRECISION   TIMMIN
*     ..
*     .. Array Arguments ..
      INTEGER            IWORK( * ), KVAL( * ), LDAVAL( * ), MVAL( * ),
     $                   NBVAL( * ), NVAL( * ), NXVAL( * )
      DOUBLE PRECISION   RESLTS( LDR1, LDR2, * ), RWORK( * )
      COMPLEX*16         A( * ), COPYA( * ), B( * ), WORK( * )
*     ..
*
*  Purpose
*  =======
*
*  ZTIMRR times the Rank-Revealing QR factorization of a
*  COMPLEX*16 general matrix.
*
*  Two matrix types may be used for timing.  The number of types is
*  set in the parameter NMODE and the matrix types are set in the vector
*  MODES, using the following key:
*     2.  BREAK1    D(1:N-1)=1 and D(N)=1.0/COND in ZLATMS
*     3.  GEOM      D(I)=COND**(-(I-1)/(N-1)) in ZLATMS
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
*  NK      (input) INTEGER
*          The number of values of K in the vector KVAL.
*
*  KVAL    (input) INTEGER array, dimension (NK)
*          The values of the matrix dimension K, used in SORMQR.
*
*  NNB     (input) INTEGER
*          The number of values of NB and NX contained in the
*          vectors NBVAL and NXVAL.  The blocking parameters are used
*          in pairs (NB,NX).
*
*  NBVAL   (input) INTEGER array, dimension (NNB)
*          The values of the blocksize NB.
*
*  NXVAL   (input) INTEGER array, dimension (NNB)
*          The values of the crossover point NX.
*
*  NLDA    (input) INTEGER
*          The number of values of LDA contained in the vector LDAVAL.
*
*  LDAVAL  (input) INTEGER array, dimension (NLDA)
*          The values of the leading dimension of the array A.
*
*  TIMMIN  (input) DOUBLE PRECISION
*          The minimum time a subroutine will be timed.
*
*  A       (workspace) COMPLEX*16 array, dimension (LDAMAX*NMAX)
*          where LDAMAX and NMAX are the maximum values of LDA and N.
*
*  COPYA   (workspace) COMPLEX*16 array, dimension (LDAMAX*NMAX)
*
*  B       (workspace) COMPLEX*16 array, dimension (LDAMAX*NMAX)
*
*
*  WORK    (workspace) COMPLEX*16 array, dimension (3*max(MMAX,NMAX))
*
*  RWORK   (workspace) DOUBLE PRECISION array, dimension (2*NMAX)
*
*  IWORK   (workspace) INTEGER array, dimension (NMAX)
*
*  RESLTS  (workspace) DOUBLE PRECISION array, dimension
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
      PARAMETER          ( NSUBS = 2, NMODE = 2 )
      DOUBLE PRECISION   ONE



      PARAMETER          ( ONE = 1.0D+0 )

*     ..
*     .. Local Scalars ..
      CHARACTER*3        PATH
      CHARACTER*6        CNAME
      INTEGER            I, IC, ICL, ILDA, IM, IMODE, INB, INFO, IK,
     $                   JOB, K, LDA, LW, M, MINMN, MODE, N, NX, NB,
     $                   RANK
      DOUBLE PRECISION   COND, DMAX, OPS, IRCOND, ORCOND, S1, S2,
     $                   TIME, UNTIME
*     ..
*     .. Local Arrays ..
      LOGICAL            TIMSUB( NSUBS )
      CHARACTER*6        SUBNAM( NSUBS )
      INTEGER            ISEED( 4 ), MODES( NMODE )
      DOUBLE PRECISION   SVLUES( 4 )
*     ..
*     .. External Functions ..
      DOUBLE PRECISION   DLAMCH, DMFLOP, DOPLA, DSECND
      EXTERNAL           DLAMCH, DMFLOP, DOPLA, DSECND
*     ..
*     .. External Subroutines ..
      EXTERNAL           ATIMCK, ATIMIN, ZGEQPX, ZGEQPY,
     $                   ZLACPY, ZLATMS
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          DBLE, MIN
*     ..
*     .. Data statements ..
      DATA               SUBNAM / 'ZGEQPX', 'ZGEQPY' /
      DATA               MODES / 2, 3 /
      DATA               ISEED / 0, 0, 0, 1 /
*     ..
*     .. Executable Statements ..
*
*     Extract the timing request from the input line.
*
      PATH( 1: 1 ) = 'Zomplex precision'
      PATH( 2: 3 ) = 'RR'
      CALL ATIMIN( PATH, LINE, NSUBS, SUBNAM, TIMSUB, NOUT, INFO )
      IF( .NOT.TIMSUB( 1 ) .OR. INFO.NE.0 )
     $   GO TO 1000
*
*     Check that M <= LDA for the input values.
*
      CNAME = LINE( 1: 6 )
      CALL ATIMCK( 1, CNAME, NM, MVAL, NLDA, LDAVAL, NOUT, INFO )
      IF( INFO.GT.0 ) THEN
         WRITE( NOUT, FMT = 9999 )CNAME
         GO TO 1000
      END IF
*
*     Set the condition number and scaling factor for the matrices
*     to be generated.
*
      DMAX = ONE
      IRCOND = DLAMCH( 'Precision' )
      COND = ONE / IRCOND
*
*     Do for each value of K:
*
      DO 10 IK = 1, NK
         K = KVAL( IK )
         IF( K.EQ.0 ) THEN
            JOB = 1
         ELSE
            JOB = 2
         END IF
*
*        Do for each type of matrix:
*
         DO 20 IMODE = 1, NMODE
            MODE = MODES( IMODE )
*
*
*           *****************
*           * Timing xGEQPX *
*           *****************
*
*           Do for each value of LDA:
*
            DO 30 ILDA = 1, NLDA
               LDA = LDAVAL( ILDA )
*
*              Do for each pair of values (M,N):
*
               DO 40 IM = 1, NM
                  M = MVAL( IM )
                  N = NVAL( IM )
                  MINMN = MIN( M, N )
*
*                 Generate a test matrix of size m by n using the
*                 singular value distribution indicated by MODE.
*
                  CALL ZLATMS( M, N, 'Uniform', ISEED, 'Nonsymm',
     $                         RWORK, MODE, COND, DMAX, M, N,
     $                         'No packing', COPYA, LDA, WORK, INFO )
*
*                 Do for each pair of values ( NB, NX ) in NBVAL and NXVAL:
*
                  DO 50 INB = 1, NNB
                     NB = NBVAL( INB )
                     CALL XLAENV( 1, NB )
                     NX = NXVAL( INB )
                     CALL XLAENV( 3, NX )
*
                     IF( JOB.EQ.1 ) THEN
                        IF( NB.LT.3 ) THEN
                           LW = MAX( 1, 2*MINMN+3*N )
                        ELSE
                           LW = MAX( 1, 2*MINMN+NB*N )
                        END IF
                     ELSE
                        IF( NB.LT.3 ) THEN
                           LW = MAX( 1, 2*MINMN+2*N+MAX(K,N) )
                        ELSE
                           LW = MAX( 1, 2*MINMN+NB*NB+NB*MAX(K,N) )
                        END IF
                     END IF
*
*                    ZGEQPX:  RRQR factorization
*
                     CALL ZLACPY( 'All', M, N, COPYA, LDA, A, LDA )
                     IC = 0
                     S1 = DSECND( )
   60                CONTINUE
*
                     CALL ZGEQPX( JOB, M, N, K, A, LDA, B, LDA,
     $                            IWORK, IRCOND, ORCOND, RANK, SVLUES,
     $                            WORK, LW, RWORK, INFO )
                     S2 = DSECND( )
*
                     IF( INFO.NE.0 ) THEN
                        WRITE( *,* ) '>>>Warning: INFO returned by ',
     $                               'ZGEQPX is:', INFO
                        INFO = 0
                     END IF
*
                     TIME = S2 - S1
                     IC = IC + 1
                     IF( TIME.LT.TIMMIN ) THEN
                        CALL ZLACPY( 'All', M, N, COPYA, LDA,
     $                                A, LDA )
                        GO TO 60
                     END IF
*
*                    Subtract the time used in ZLACPY.
*
                     ICL = 1
                     S1 = DSECND( )
   70                CONTINUE
                     S2 = DSECND( )
                     UNTIME = S2 - S1
                     ICL = ICL + 1
                     IF( ICL.LE.IC ) THEN
                        CALL ZLACPY( 'All', M, N, COPYA, LDA,
     $                                A, LDA )
                        GO TO 70
                     END IF
*
                     TIME = ( TIME-UNTIME ) / DBLE( IC )
*
*                    The number of flops of yGEQPX is approximately the
*                    the number of flops of yGEQPF plus the number of
*                    flops required by yUNMQR to update matrix C.
*
                     OPS = DOPLA( 'ZGEQPF', M, N, 0, 0, NB )
     $                + DOPLA( 'ZUNMQR', M, K, MINMN, 0, NB )
                     RESLTS( INB, IM, ILDA ) =
     $                              DMFLOP( OPS, TIME, INFO )
*
   50             CONTINUE
   40          CONTINUE
   30       CONTINUE
*
*           Print the results for each value of K and type of matrix.
*
            WRITE( NOUT, FMT = 9995 )SUBNAM( 1 )
            WRITE( NOUT, FMT = 9996 )K, IMODE
            DO 90 I = 1, NLDA
               WRITE( NOUT, FMT = 9997 )I,LDAVAL( I )
   90       CONTINUE
            WRITE( NOUT, FMT = * )
            CALL DPRTB4( '(  NB,  NX)', 'M', 'N', NNB, NBVAL,
     $                NXVAL, NM, MVAL, NVAL, NLDA, RESLTS( 1, 1, 1 ),
     $                LDR1, LDR2, NOUT )
*
*
*           *****************
*           * Timing xGEQPY *
*           *****************
*
*           Do for each value of LDA:
*
            DO 200 ILDA = 1, NLDA
               LDA = LDAVAL( ILDA )
*
*              Do for each pair of values (M,N):
*
               DO 210 IM = 1, NM
                  M = MVAL( IM )
                  N = NVAL( IM )
                  MINMN = MIN( M, N )
*
*                 Generate a test matrix of size m by n using the
*                 singular value distribution indicated by MODE.
*
                  CALL ZLATMS( M, N, 'Uniform', ISEED, 'Nonsymm',
     $                         RWORK, MODE, COND, DMAX, M, N,
     $                         'No packing', COPYA, LDA, WORK, INFO )
*
*                 Do for each pair of values ( NB, NX ) in NBVAL and NXVAL:
*
                  DO 220 INB = 1, NNB
                     NB = NBVAL( INB )
                     CALL XLAENV( 1, NB )
                     NX = NXVAL( INB )
                     CALL XLAENV( 3, NX )
*
                     IF( JOB.EQ.1 ) THEN
                        IF( NB.LT.3 ) THEN
                           LW = MAX( 1, 2*MINMN+3*N )
                        ELSE
                           LW = MAX( 1, 2*MINMN+NB*N )
                        END IF
                     ELSE
                        IF( NB.LT.3 ) THEN
                           LW = MAX( 1, 2*MINMN+2*N+MAX(K,N) )
                        ELSE
                           LW = MAX( 1, 2*MINMN+NB*NB+NB*MAX(K,N) )
                        END IF
                     END IF
*
*                    ZGEQPY:  RRQR factorization
*
                     CALL ZLACPY( 'All', M, N, COPYA, LDA, A, LDA )
                     IC = 0
                     S1 = DSECND( )
  230                CONTINUE
*
                     CALL ZGEQPY( JOB, M, N, K, A, LDA, B, LDA,
     $                            IWORK, IRCOND, ORCOND, RANK, SVLUES,
     $                            WORK, LW, RWORK, INFO )
                     S2 = DSECND( )
*
                     IF( INFO.NE.0 ) THEN
                        WRITE( *,* ) '>>>Warning: INFO returned by ',
     $                               'ZGEQPY is:', INFO
                        INFO = 0
                     END IF
*
                     TIME = S2 - S1
                     IC = IC + 1
                     IF( TIME.LT.TIMMIN ) THEN
                        CALL ZLACPY( 'All', M, N, COPYA, LDA,
     $                                A, LDA )
                        GO TO 230
                     END IF
*
*                    Subtract the time used in ZLACPY.
*
                     ICL = 1
                     S1 = DSECND( )
  240                CONTINUE
                     S2 = DSECND( )
                     UNTIME = S2 - S1
                     ICL = ICL + 1
                     IF( ICL.LE.IC ) THEN
                        CALL ZLACPY( 'All', M, N, COPYA, LDA,
     $                                A, LDA )
                        GO TO 240
                     END IF
*
*                    The number of flops of yGEQPY is approximately the
*                    the number of flops of yGEQPF plus the number of
*                    flops required by yUNMQR to update matrix C.
*
                     TIME = ( TIME-UNTIME ) / DBLE( IC )
*
                     OPS = DOPLA( 'ZGEQPF', M, N, 0, 0, NB )
     $                + DOPLA( 'ZUNMQR', M, K, MINMN, 0, NB )
                     RESLTS( INB, IM, ILDA ) =
     $                              DMFLOP( OPS, TIME, INFO )
*
  220             CONTINUE
  210          CONTINUE
  200       CONTINUE
*
*           Print the results for each value of K and type of matrix.
*
            WRITE( NOUT, FMT = 9995 )SUBNAM( 2 )
            WRITE( NOUT, FMT = 9996 )K, IMODE
            DO 250 I = 1, NLDA
               WRITE( NOUT, FMT = 9997 )I,LDAVAL( I )
  250       CONTINUE
            WRITE( NOUT, FMT = * )
            CALL DPRTB4( '(  NB,  NX)', 'M', 'N', NNB, NBVAL,
     $                NXVAL, NM, MVAL, NVAL, NLDA, RESLTS( 1, 1, 1 ),
     $                LDR1, LDR2, NOUT )
*
*
   20    CONTINUE
   10 CONTINUE
*
 9995 FORMAT( / ' *** Speed of ', A6, ' in megaflops ***' )
 9996 FORMAT( 5X, 'with K = ', I4, '  and type of matrix:', I4 )
 9997 FORMAT( 5X, 'line ', I4, ' with LDA = ', I4 )
 9999 FORMAT( 1X, A6, ' timing run not attempted', / )
*
 1000 CONTINUE
      RETURN
*
*     End of ZTIMRR
*
      END
