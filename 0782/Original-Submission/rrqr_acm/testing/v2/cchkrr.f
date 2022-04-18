      SUBROUTINE CCHKRR( DOTYPE, NM, MVAL, NN, NVAL, NNB, NBVAL,
     $                   NXVAL, THRESH, TSTERR, A, AQ, COPYA, S, COPYS,
     $                   WORK, LNWK, RWORK, IWORK, NOUT )
*
*  -- LAPACK test routine (version 2.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     September 30, 1994
*
*     Rewritten to test RRQR Subroutines.
*
*     .. Scalar Arguments ..
      LOGICAL            TSTERR
      INTEGER            LNWK, NM, NN, NNB, NOUT
      REAL               THRESH
*     ..
*     .. Array Arguments ..
      LOGICAL            DOTYPE( * )
      INTEGER            IWORK( * ), MVAL( * ), NVAL( * ),
     $                   NBVAL( * ), NXVAL( * )
      REAL               COPYS( * ), RWORK( * ), S( * )
      COMPLEX            A( * ), AQ( * ), COPYA( * ), WORK( * )
*     ..
*
*  Purpose
*  =======
*
*  CCHKRR tests CGEQPX and CGEQPY.
*
*  Arguments
*  =========
*
*  DOTYPE  (input) LOGICAL array, dimension (NTYPES)
*          The matrix types to be used for testing.  Matrices of type j
*          (for 1 <= j <= NTYPES) are used for testing if DOTYPE(j) =
*          .TRUE.; if DOTYPE(j) = .FALSE., then type j is not used.
*
*  NM      (input) INTEGER
*          The number of values of M contained in the vector MVAL.
*
*  MVAL    (input) INTEGER array, dimension (NM)
*          The values of the matrix row dimension M.
*
*  NN      (input) INTEGER
*          The number of values of N contained in the vector NVAL.
*
*  NVAL    (input) INTEGER array, dimension (NN)
*          The values of the matrix column dimension N.
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
*  THRESH  (input) REAL
*          The threshold value for the test ratios.  A result is
*          included in the output file if RESULT >= THRESH.  To have
*          every test ratio printed, use THRESH = 0.
*
*  TSTERR  (input) LOGICAL
*          Flag that indicates whether error exits are to be tested.
*
*  A       (workspace) COMPLEX array, dimension (MMAX*NMAX)
*          where MMAX is the maximum value of M in MVAL and NMAX is the
*          maximum value of N in NVAL.
*
*  AQ      (workspace) COMPLEX array, dimension (MMAX*NMAX)
*          where MMAX is the maximum value of M in MVAL and NMAX is the
*          maximum value of N in NVAL.
*
*  COPYA   (workspace) COMPLEX array, dimension (MMAX*NMAX)
*
*  S       (workspace) REAL array, dimension
*                      (min(MMAX,NMAX))
*
*  COPYS   (workspace) REAL array, dimension
*                      (min(MMAX,NMAX))
*
*  WORK    (workspace) COMPLEX array, dimension (LNWK)
*
*  LNWK    (input) INTEGER
*          Workspace length. At least
*            ((max(M*max(M,N) + 4*min(M,N) + max(M,N))).
*
*  RWORK   (workspace) REAL array, dimension (4*NMAX)
*
*  IWORK   (workspace) INTEGER array, dimension (NMAX)
*
*  NOUT    (input) INTEGER
*          The unit number for output.
*
*  =====================================================================
*
*     .. Parameters ..
      INTEGER            NTYPES
      PARAMETER          ( NTYPES = 3 )
      INTEGER            NTESTS
      PARAMETER          ( NTESTS = 14 )
      REAL               ONE, ZERO
      COMPLEX            CONE, CZERO

      PARAMETER          ( ONE = 1.0E+0, ZERO = 0.0E+0,
     $                   CONE = ( 1.0E+0, 0.0E+0 ),
     $                   CZERO = ( 0.0E+0, 0.0E+0 ) )





*     ..
*     .. Local Scalars ..
      CHARACTER*3        PATH
      INTEGER            I, IM, IMODE, IN, INB, INFO, K,
     $                   LDA, LW, LWORK, M, MNMIN, N, NB, NERRS, NFAIL,
     $                   NRUN, NX, RANK
      REAL               EPS, IRCOND, ORCOND
*     ..
*     .. Local Arrays ..
      INTEGER            ISEED( 4 ), ISEEDY( 4 )
      REAL               RESULT( NTESTS ), SVLUES( 4 )
*     ..
*     .. External Functions ..
      REAL               SLAMCH, CQRT12, CRRT01,
     $                   CRRT02
      EXTERNAL           SLAMCH, CQRT12, CRRT01,
     $                   CRRT02
*     ..
*     .. External Subroutines ..
      EXTERNAL           ALAHD, ALASUM, SLAORD, CERRRR,
     $                   CGEQPX, CGEQPY,
     $                   CLACPY, CLASET, CLATMS
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
*     ..
*     .. Scalars in Common ..
      LOGICAL            LERR, OK
      CHARACTER*6        SRNAMT
      INTEGER            INFOT, IOUNIT
*     ..
*     .. Common blocks ..
      COMMON             / INFOC / INFOT, IOUNIT, OK, LERR
      COMMON             / SRNAMC / SRNAMT
*     ..
*     .. Data statements ..
      DATA               ISEEDY / 1988, 1989, 1990, 1991 /
*     ..
*     .. Executable Statements ..
*
*     Initialize constants and the random number seed.
*
      PATH( 1: 1 ) = 'C'
      PATH( 2: 3 ) = 'RR'
      NRUN = 0
      NFAIL = 0
      NERRS = 0
      DO 10 I = 1, 4
         ISEED( I ) = ISEEDY( I )
   10 CONTINUE
      EPS = SLAMCH( 'Epsilon' )
      IRCOND = EPS
*
*     Test the error exits
*
      IF( TSTERR )
     $   CALL CERRRR( PATH, NOUT )
      INFOT = 0
*
      DO 80 IM = 1, NM
*
*        Do for each value of M in MVAL.
*
         M = MVAL( IM )
         LDA = MAX( 1, M )
*
         DO 70 IN = 1, NN
*
*           Do for each value of N in NVAL.
*
            N = NVAL( IN )
            MNMIN = MIN( M, N )
            LWORK = MAX( 1, M*MAX( M, N )+4*MNMIN+MAX( M, N ) )
*
*           Check if there is enough workspace.
*
            IF( LWORK.GT.LNWK ) THEN
               WRITE(*,*) ' Error in CCHKRR.', 'Code 1:',
     $                    ' Actual Workspace:', LNWK,
     $                    '   Needed Workspace :', LWORK
               STOP
            ENDIF
*
            DO 60 IMODE = 1, NTYPES
               IF( .NOT.DOTYPE( IMODE ) )
     $            GO TO 60
*
*              Do for each type of matrix
*                 1:  zero matrix
*                 2:  one small singular value
*                 3:  geometric distribution of singular values
*
*              Generate test matrix of size m by n using
*              singular value distribution indicated by 'mode'.
*
               IF( IMODE.EQ.1 ) THEN
                  CALL CLASET( 'Full', M, N, CZERO, CZERO, COPYA,
     $                          LDA )
                  DO 30 I = 1, MNMIN
                     COPYS( I ) = ZERO
   30             CONTINUE
               ELSE
                  CALL CLATMS( M, N, 'Uniform', ISEED, 'Nonsymm',
     $                         COPYS, IMODE, ONE / EPS, ONE, M, N,
     $                         'No packing', COPYA, LDA, WORK, INFO )
                  CALL SLAORD( 'Decreasing', MNMIN, COPYS, 1 )
               END IF
*
               DO 90 INB = 1, NNB
*
*                 Do for each pair of values (NB,NX) in NBVAL and NXVAL.
*
                  NB = NBVAL( INB )
                  CALL XLAENV( 1, NB )
                  NX = NXVAL( INB )
                  CALL XLAENV( 3, NX )
*
*                 ******************
*                 * Testing yGEQPX *
*                 ******************
*
*                 *-*-*-*-*-*-*-*-*-*-*-*-*-*
*                 * Test yGEQPX when JOB=1  *
*                 *-*-*-*-*-*-*-*-*-*-*-*-*-*
*
*                 Get a work copy of matrix A from COPYA.
*
                  CALL CLACPY( 'All', M, N, COPYA, LDA, A, LDA )
*
*                 Compute the Rank-Revealing QR factorization of A
*
                  IF ( NB.LT.3 ) THEN
                     LW = MAX( 1, 2*MNMIN+N )
                  ELSE
                     LW = MAX( 1, 2*MNMIN+N*NB )
                  END IF
*
*                 Check if there is enough workspace for current
*                 block size.
*
                  IF( LW.GT.LNWK ) THEN
                     WRITE(*,*) ' Error in PREC_CHKRR.', 'Code 2',
     $                          ' Workspace too short for blocksize',nb
                     WRITE(*,*) ' Actual Workspace:', LNWK,
     $                          '   Needed Workspace :', LW
                     STOP
                  ENDIF
*
                  SRNAMT = 'CGEQPX'
                  CALL CGEQPX( 1, M, N, 0, A, LDA, AQ, LDA, IWORK,
     $                         IRCOND, ORCOND, RANK, SVLUES, WORK, LW,
     $                         RWORK, INFO )
*
*                 Compute norm(svd(a) - svd(r))
*
                  RESULT( 1 ) = CQRT12( M, N, A, LDA, COPYS, WORK,
     $                          LWORK, RWORK )

*
*                 *-*-*-*-*-*-*-*-*-*-*-*-*-*
*                 * Test yGEQPX when JOB=2  *
*                 *-*-*-*-*-*-*-*-*-*-*-*-*-*
*
*                 Get a work copy of matrix A from COPYA.
*
                  CALL CLACPY( 'All', M, N, COPYA, LDA, A, LDA )
*
*                 Set WORK to identity.
*
                  CALL CLASET( 'All', M, M, CZERO, CONE, AQ, LDA )
*
*                 Compute the Rank-Revealing QR factorization of A
*
                  IF ( NB.LT.3 ) THEN
                     LW = MAX( 1, 2*MNMIN+MAX(N,M) )
                  ELSE
                     LW = MAX( 1, 2*MNMIN+NB*NB+NB*MAX(N,M) )
                  END IF

                  SRNAMT = 'CGEQPX'
                  CALL CGEQPX( 2, M, N, M, A, LDA, AQ, LDA, IWORK,
     $                         IRCOND, ORCOND, RANK, SVLUES, WORK, LW,
     $                         RWORK, INFO )
*
*                 Compute norm(svd(a) - svd(r))
*
                  RESULT( 2 ) = CQRT12( M, N, A, LDA, COPYS, WORK,
     $                          LWORK, RWORK )
*
*                 Compute norm( A*P - Q*R )
*
                  RESULT( 3 ) = CRRT01( 'Conjugate Transpose',
     $                          M, N, COPYA, LDA, IWORK,
     $                          AQ, LDA, A, LDA, WORK, LWORK )
*
*                 Compute norm( Q'*Q - Identity )
*
                  RESULT( 4 ) = CRRT02( M, AQ, LDA, WORK, LWORK )
*
*                 *-*-*-*-*-*-*-*-*-*-*-*-*-*
*                 * Test yGEQPX when JOB=3  *
*                 *-*-*-*-*-*-*-*-*-*-*-*-*-*
*
*                 Get a work copy of matrix A from COPYA.
*
                  CALL CLACPY( 'All', M, N, COPYA, LDA, A, LDA )
*
*                 Set WORK to identity.
*
                  CALL CLASET( 'All', M, M, CZERO, CONE, AQ, LDA )
*
*                 Compute the Rank-Revealing QR factorization of A
*
                  IF ( NB.LT.3 ) THEN
                     LW = MAX( 1, 2*MNMIN+MAX(N,M) )
                  ELSE
                     LW = MAX( 1, 2*MNMIN+NB*NB+NB*MAX(N,M) )
                  END IF

                  SRNAMT = 'CGEQPX'
                  CALL CGEQPX( 3, M, N, M, A, LDA, AQ, LDA, IWORK,
     $                         IRCOND, ORCOND, RANK, SVLUES, WORK, LW,
     $                         RWORK, INFO )
*
*                 Compute norm(svd(a) - svd(r))
*
                  RESULT( 5 ) = CQRT12( M, N, A, LDA, COPYS, WORK,
     $                          LWORK, RWORK )
*
*                 Compute norm( A*P - Q*R )
*
                  RESULT( 6 ) = CRRT01( 'No transpose', M, N,
     $                          COPYA, LDA, IWORK,
     $                          AQ, LDA, A, LDA, WORK, LWORK )
*
*                 Compute norm( Q'*Q - Identity )
*
                  RESULT( 7 ) = CRRT02( M, AQ, LDA, WORK, LWORK )
*
*
*                 *-*-*-*-*-*-*-*-*-*-*
*                 * Printing results  *
*                 *-*-*-*-*-*-*-*-*-*-*
*
*                 Print information about the tests that did not pass
*                 the threshold.
*
                  DO 100 K = 1, 7
                     IF( RESULT( K ).GE.THRESH ) THEN
                        IF( NFAIL.EQ.0 .AND. NERRS.EQ.0 )
     $                     CALL ALAHD( NOUT, PATH )
                        WRITE( NOUT, FMT = 9999 ) 'CGEQPX', M, N,
     $                                            NB, IMODE, K,
     $                                            RESULT( K )
                        NFAIL = NFAIL + 1
                     END IF
  100             CONTINUE
                  NRUN = NRUN + 7
*
*
*                 ******************
*                 * Testing yGEQPY *
*                 ******************
*
*                 *-*-*-*-*-*-*-*-*-*-*-*-*-*
*                 * Test yGEQPY when JOB=1  *
*                 *-*-*-*-*-*-*-*-*-*-*-*-*-*
*
*                 Get a work copy of matrix A from COPYA.
*
                  CALL CLACPY( 'All', M, N, COPYA, LDA, A, LDA )
*
*                 Compute the Rank-Revealing QR factorization of A
*
                  IF ( NB.LT.3 ) THEN
                     LW = MAX( 1, 2*MNMIN+N )
                  ELSE
                     LW = MAX( 1, 2*MNMIN+N*NB )
                  END IF
*
*                 Check if there is enough workspace for current
*                 block size.
*
                  IF( LW.GT.LNWK ) THEN
                     WRITE(*,*) ' Error in PREC_CHKRR.', 'Code 2',
     $                          ' Workspace too short for blocksize',nb
                     WRITE(*,*) ' Actual Workspace:', LNWK,
     $                          '   Needed Workspace :', LW
                     STOP
                  ENDIF
*
                  SRNAMT = 'CGEQPY'
                  CALL CGEQPY( 1, M, N, 0, A, LDA, AQ, LDA, IWORK,
     $                         IRCOND, ORCOND, RANK, SVLUES, WORK, LW,
     $                         RWORK, INFO )
*
*                 Compute norm(svd(a) - svd(r))
*
                  RESULT( 8 ) = CQRT12( M, N, A, LDA, COPYS, WORK,
     $                          LWORK, RWORK )

*
*                 *-*-*-*-*-*-*-*-*-*-*-*-*-*
*                 * Test yGEQPY when JOB=2  *
*                 *-*-*-*-*-*-*-*-*-*-*-*-*-*
*
*                 Get a work copy of matrix A from COPYA.
*
                  CALL CLACPY( 'All', M, N, COPYA, LDA, A, LDA )
*
*                 Set WORK to identity.
*
                  CALL CLASET( 'All', M, M, CZERO, CONE, AQ, LDA )
*
*                 Compute the Rank-Revealing QR factorization of A
*
                  IF ( NB.LT.3 ) THEN
                     LW = MAX( 1, 2*MNMIN+MAX(N,M) )
                  ELSE
                     LW = MAX( 1, 2*MNMIN+NB*NB+NB*MAX(N,M) )
                  END IF

                  SRNAMT = 'CGEQPY'
                  CALL CGEQPY( 2, M, N, M, A, LDA, AQ, LDA, IWORK,
     $                         IRCOND, ORCOND, RANK, SVLUES, WORK, LW,
     $                         RWORK, INFO )
*
*                 Compute norm(svd(a) - svd(r))
*
                  RESULT( 9 ) = CQRT12( M, N, A, LDA, COPYS, WORK,
     $                          LWORK, RWORK )
*
*                 Compute norm( A*P - Q*R )
*
                  RESULT( 10 ) = CRRT01( 'Conjugate Transpose',
     $                          M, N, COPYA, LDA, IWORK,
     $                          AQ, LDA, A, LDA, WORK, LWORK )
*
*                 Compute norm( Q'*Q - Identity )
*
                  RESULT( 11 ) = CRRT02( M, AQ, LDA, WORK, LWORK )
*
*                 *-*-*-*-*-*-*-*-*-*-*-*-*-*
*                 * Test yGEQPY when JOB=3  *
*                 *-*-*-*-*-*-*-*-*-*-*-*-*-*
*
*                 Get a work copy of matrix A from COPYA.
*
                  CALL CLACPY( 'All', M, N, COPYA, LDA, A, LDA )
*
*                 Set WORK to identity.
*
                  CALL CLASET( 'All', M, M, CZERO, CONE, AQ, LDA )
*
*                 Compute the Rank-Revealing QR factorization of A
*
                  IF ( NB.LT.3 ) THEN
                     LW = MAX( 1, 2*MNMIN+MAX(N,M) )
                  ELSE
                     LW = MAX( 1, 2*MNMIN+NB*NB+NB*MAX(N,M) )
                  END IF

                  SRNAMT = 'CGEQPY'
                  CALL CGEQPY( 3, M, N, M, A, LDA, AQ, LDA, IWORK,
     $                         IRCOND, ORCOND, RANK, SVLUES, WORK, LW,
     $                         RWORK, INFO )
*
*                 Compute norm(svd(a) - svd(r))
*
                  RESULT( 12 ) = CQRT12( M, N, A, LDA, COPYS, WORK,
     $                          LWORK, RWORK )
*
*                 Compute norm( A*P - Q*R )
*
                  RESULT( 13 ) = CRRT01( 'No transpose', M, N,
     $                          COPYA, LDA, IWORK,
     $                          AQ, LDA, A, LDA, WORK, LWORK )
*
*                 Compute norm( Q'*Q - Identity )
*
                  RESULT( 14 ) = CRRT02( M, AQ, LDA, WORK, LWORK )
*
*
*                 *-*-*-*-*-*-*-*-*-*-*
*                 * Printing results  *
*                 *-*-*-*-*-*-*-*-*-*-*
*
*                 Print information about the tests that did not pass
*                 the threshold.
*
                  DO 110 K = 8, 14
                     IF( RESULT( K ).GE.THRESH ) THEN
                        IF( NFAIL.EQ.0 .AND. NERRS.EQ.0 )
     $                     CALL ALAHD( NOUT, PATH )
                        WRITE( NOUT, FMT = 9999 ) 'CGEQPY', M, N,
     $                                            NB, IMODE, K,
     $                                            RESULT( K )
                        NFAIL = NFAIL + 1
                     END IF
  110             CONTINUE
                  NRUN = NRUN + 7
*
   90          CONTINUE
   60       CONTINUE
   70    CONTINUE
   80 CONTINUE
*
*     Print a summary of the results.
*
      CALL ALASUM( PATH, NOUT, NFAIL, NRUN, NERRS )
*
 9999 FORMAT( 1x, a6, ' M =', I5, ', N =', I5, ', NB =', I4,
     $        ', type ', I2, ', test ', I2, ', ratio =', G12.5 )
*
*     End of CCHKRR
*
      END
