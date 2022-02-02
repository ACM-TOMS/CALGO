      SUBROUTINE DTIMQR( LINE, NM, MVAL, NVAL, NK, KVAL, NNB, NBVAL,
     $                   NXVAL, NLDA, LDAVAL, TIMMIN, A, TAU, B, WORK,
     $                   RESLTS, LDR1, LDR2, LDR3, NOUT )
*
*  -- LAPACK timing routine (version 2.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     March 31, 1993
*
*     .. Scalar Arguments ..
      CHARACTER*80       LINE
      INTEGER            LDR1, LDR2, LDR3, NK, NLDA, NM, NNB, NOUT
      DOUBLE PRECISION   TIMMIN
*     ..
*     .. Array Arguments ..
      INTEGER            KVAL( * ), LDAVAL( * ), MVAL( * ), NBVAL( * ),
     $                   NVAL( * ), NXVAL( * )
      DOUBLE PRECISION   A( * ), B( * ), RESLTS( LDR1, LDR2, LDR3, * ),
     $                   TAU( * ), WORK( * )
*     ..
*
*  Purpose
*  =======
*
*  DTIMQR times the LAPACK routines to perform the QR factorization of
*  a DOUBLE PRECISION general matrix.
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
*          The values of the matrix dimension K, used in DORMQR.
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
*  A       (workspace) DOUBLE PRECISION array, dimension (LDAMAX*NMAX)
*          where LDAMAX and NMAX are the maximum values of LDA and N.
*
*  TAU     (workspace) DOUBLE PRECISION array, dimension (min(M,N))
*
*  B       (workspace) DOUBLE PRECISION array, dimension (LDAMAX*NMAX)
*
*  WORK    (workspace) DOUBLE PRECISION array, dimension (LDAMAX*NBMAX)
*          where NBMAX is the maximum value of NB.
*
*  RESLTS  (workspace) DOUBLE PRECISION array, dimension
*                      (LDR1,LDR2,LDR3,2*NK)
*          The timing results for each subroutine over the relevant
*          values of (M,N), (NB,NX), and LDA.
*
*  LDR1    (input) INTEGER
*          The first dimension of RESLTS.  LDR1 >= max(1,NNB).
*
*  LDR2    (input) INTEGER
*          The second dimension of RESLTS.  LDR2 >= max(1,NM).
*
*  LDR3    (input) INTEGER
*          The third dimension of RESLTS.  LDR3 >= max(1,NLDA).
*
*  NOUT    (input) INTEGER
*          The unit number for output.
*
*  Internal Parameters
*  ===================
*
*  MODE    INTEGER
*          The matrix type.  MODE = 3 is a geometric distribution of
*          eigenvalues.  See DLATMS for further details.
*
*  COND    DOUBLE PRECISION
*          The condition number of the matrix.  The singular values are
*          set to values from DMAX to DMAX/COND.
*
*  DMAX    DOUBLE PRECISION
*          The magnitude of the largest singular value.
*
*  =====================================================================
*
*     .. Parameters ..
      INTEGER            NSUBS
      PARAMETER          ( NSUBS = 3 )
      INTEGER            MODE
      DOUBLE PRECISION   COND, DMAX
      PARAMETER          ( MODE = 3, COND = 100.0D0, DMAX = 1.0D0 )
*     ..
*     .. Local Scalars ..
      CHARACTER          LABM, SIDE, TRANS
      CHARACTER*3        PATH
      CHARACTER*6        CNAME
      INTEGER            I, I4, IC, ICL, IK, ILDA, IM, IMX, INB, INFO,
     $                   ISIDE, ISUB, ITOFF, ITRAN, K, K1, LDA, LW, M,
     $                   M1, MINMN, N, N1, NB, NX
      DOUBLE PRECISION   OPS, S1, S2, TIME, UNTIME
*     ..
*     .. Local Arrays ..
      LOGICAL            TIMSUB( NSUBS )
      CHARACTER          SIDES( 2 ), TRANSS( 2 )
      CHARACTER*6        SUBNAM( NSUBS )
      INTEGER            ISEED( 4 ), MUSE( 12 ), NUSE( 12 ), RESEED( 4 )
*     ..
*     .. External Functions ..
      DOUBLE PRECISION   DMFLOP, DOPLA, DSECND
      EXTERNAL           DMFLOP, DOPLA, DSECND
*     ..
*     .. External Subroutines ..
      EXTERNAL           ATIMCK, ATIMIN, DGEQRF, DLACPY, DLATMS, DORGQR,
     $                   DORMQR, DPRTB4, DPRTB5, DTIMMG, ICOPY, XLAENV
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          DBLE, MAX, MIN
*     ..
*     .. Data statements ..
      DATA               SUBNAM / 'DGEQRF', 'DORGQR', 'DORMQR' /
      DATA               SIDES / 'L', 'R' / , TRANSS / 'N', 'T' /
      DATA               ISEED / 0, 0, 0, 1 /
*     ..
*     .. Executable Statements ..
*
*     Extract the timing request from the input line.
*
      PATH( 1: 1 ) = 'Double precision'
      PATH( 2: 3 ) = 'QR'
      CALL ATIMIN( PATH, LINE, NSUBS, SUBNAM, TIMSUB, NOUT, INFO )
      IF( INFO.NE.0 )
     $   GO TO 230
*
*     Check that M <= LDA for the input values.
*
      CNAME = LINE( 1: 6 )
      CALL ATIMCK( 1, CNAME, NM, MVAL, NLDA, LDAVAL, NOUT, INFO )
      IF( INFO.GT.0 ) THEN
         WRITE( NOUT, FMT = 9999 )CNAME
         GO TO 230
      END IF
*
*     Do for each pair of values (M,N):
*
      DO 70 IM = 1, NM
         M = MVAL( IM )
         N = NVAL( IM )
         MINMN = MIN( M, N )
         CALL ICOPY( 4, ISEED, 1, RESEED, 1 )
*
*        Do for each value of LDA:
*
         DO 60 ILDA = 1, NLDA
            LDA = LDAVAL( ILDA )
*
*           Do for each pair of values (NB, NX) in NBVAL and NXVAL.
*
            DO 50 INB = 1, NNB
               NB = NBVAL( INB )
               CALL XLAENV( 1, NB )
               NX = NXVAL( INB )
               CALL XLAENV( 3, NX )
               LW = MAX( 1, N*MAX( 1, NB ) )
*
*              Generate a test matrix of size M by N.
*
               CALL ICOPY( 4, RESEED, 1, ISEED, 1 )
               CALL DLATMS( M, N, 'Uniform', ISEED, 'Nonsymm', TAU,
     $                      MODE, COND, DMAX, M, N, 'No packing', B,
     $                      LDA, WORK, INFO )
*
               IF( TIMSUB( 1 ) ) THEN
*
*                 DGEQRF:  QR factorization
*
                  CALL DLACPY( 'Full', M, N, B, LDA, A, LDA )
                  IC = 0
                  S1 = DSECND( )
   10             CONTINUE
                  CALL DGEQRF( M, N, A, LDA, TAU, WORK, LW, INFO )
                  S2 = DSECND( )
                  TIME = S2 - S1
                  IC = IC + 1
                  IF( TIME.LT.TIMMIN ) THEN
                     CALL DLACPY( 'Full', M, N, B, LDA, A, LDA )
                     GO TO 10
                  END IF
*
*                 Subtract the time used in DLACPY.
*
                  ICL = 1
                  S1 = DSECND( )
   20             CONTINUE
                  S2 = DSECND( )
                  UNTIME = S2 - S1
                  ICL = ICL + 1
                  IF( ICL.LE.IC ) THEN
                     CALL DLACPY( 'Full', M, N, A, LDA, B, LDA )
                     GO TO 20
                  END IF
*
                  TIME = ( TIME-UNTIME ) / DBLE( IC )
                  OPS = DOPLA( 'DGEQRF', M, N, 0, 0, NB )
                  RESLTS( INB, IM, ILDA, 1 ) = DMFLOP( OPS, TIME, INFO )
               ELSE
*
*                 If DGEQRF was not timed, generate a matrix and factor
*                 it using DGEQRF anyway so that the factored form of
*                 the matrix can be used in timing the other routines.
*
                  CALL DLACPY( 'Full', M, N, B, LDA, A, LDA )
                  CALL DGEQRF( M, N, A, LDA, TAU, WORK, LW, INFO )
               END IF
*
               IF( TIMSUB( 2 ) ) THEN
*
*                 DORGQR:  Generate orthogonal matrix Q from the QR
*                 factorization
*
                  CALL DLACPY( 'Full', M, MINMN, A, LDA, B, LDA )
                  IC = 0
                  S1 = DSECND( )
   30             CONTINUE
                  CALL DORGQR( M, MINMN, MINMN, B, LDA, TAU, WORK, LW,
     $                         INFO )
                  S2 = DSECND( )
                  TIME = S2 - S1
                  IC = IC + 1
                  IF( TIME.LT.TIMMIN ) THEN
                     CALL DLACPY( 'Full', M, MINMN, A, LDA, B, LDA )
                     GO TO 30
                  END IF
*
*                 Subtract the time used in DLACPY.
*
                  ICL = 1
                  S1 = DSECND( )
   40             CONTINUE
                  S2 = DSECND( )
                  UNTIME = S2 - S1
                  ICL = ICL + 1
                  IF( ICL.LE.IC ) THEN
                     CALL DLACPY( 'Full', M, MINMN, A, LDA, B, LDA )
                     GO TO 40
                  END IF
*
                  TIME = ( TIME-UNTIME ) / DBLE( IC )
                  OPS = DOPLA( 'DORGQR', M, MINMN, MINMN, 0, NB )
                  RESLTS( INB, IM, ILDA, 2 ) = DMFLOP( OPS, TIME, INFO )
               END IF
*
   50       CONTINUE
   60    CONTINUE
   70 CONTINUE
*
*     Print tables of results
*
      DO 90 ISUB = 1, NSUBS - 1
         IF( .NOT.TIMSUB( ISUB ) )
     $      GO TO 90
         WRITE( NOUT, FMT = 9998 )SUBNAM( ISUB )
         IF( NLDA.GT.1 ) THEN
            DO 80 I = 1, NLDA
               WRITE( NOUT, FMT = 9997 )I, LDAVAL( I )
   80       CONTINUE
         END IF
         WRITE( NOUT, FMT = * )
         IF( ISUB.EQ.2 )
     $      WRITE( NOUT, FMT = 9996 )
         CALL DPRTB4( '(  NB,  NX)', 'M', 'N', NNB, NBVAL, NXVAL, NM,
     $                MVAL, NVAL, NLDA, RESLTS( 1, 1, 1, ISUB ), LDR1,
     $                LDR2, NOUT )
   90 CONTINUE
*
*     Time DORMQR separately.  Here the starting matrix is M by N, and
*     K is the free dimension of the matrix multiplied by Q.
*
      IF( TIMSUB( 3 ) ) THEN
*
*        Check that K <= LDA for the input values.
*
         CALL ATIMCK( 3, CNAME, NK, KVAL, NLDA, LDAVAL, NOUT, INFO )
         IF( INFO.GT.0 ) THEN
            WRITE( NOUT, FMT = 9999 )SUBNAM( 3 )
            GO TO 230
         END IF
*
*        Use only the pairs (M,N) where M >= N.
*
         IMX = 0
         DO 100 IM = 1, NM
            IF( MVAL( IM ).GE.NVAL( IM ) ) THEN
               IMX = IMX + 1
               MUSE( IMX ) = MVAL( IM )
               NUSE( IMX ) = NVAL( IM )
            END IF
  100    CONTINUE
*
*        DORMQR:  Multiply by Q stored as a product of elementary
*        transformations
*
*        Do for each pair of values (M,N):
*
         DO 180 IM = 1, IMX
            M = MUSE( IM )
            N = NUSE( IM )
*
*           Do for each value of LDA:
*
            DO 170 ILDA = 1, NLDA
               LDA = LDAVAL( ILDA )
*
*              Generate an M by N matrix and form its QR decomposition.
*
               CALL DLATMS( M, N, 'Uniform', ISEED, 'Nonsymm', TAU,
     $                      MODE, COND, DMAX, M, N, 'No packing', A,
     $                      LDA, WORK, INFO )
               LW = MAX( 1, N*MAX( 1, NB ) )
               CALL DGEQRF( M, N, A, LDA, TAU, WORK, LW, INFO )
*
*              Do first for SIDE = 'L', then for SIDE = 'R'
*
               I4 = 0
               DO 160 ISIDE = 1, 2
                  SIDE = SIDES( ISIDE )
*
*                 Do for each pair of values (NB, NX) in NBVAL and
*                 NXVAL.
*
                  DO 150 INB = 1, NNB
                     NB = NBVAL( INB )
                     CALL XLAENV( 1, NB )
                     NX = NXVAL( INB )
                     CALL XLAENV( 3, NX )
*
*                    Do for each value of K in KVAL
*
                     DO 140 IK = 1, NK
                        K = KVAL( IK )
*
*                       Sort out which variable is which
*
                        IF( ISIDE.EQ.1 ) THEN
                           M1 = M
                           K1 = N
                           N1 = K
                           LW = MAX( 1, N1*MAX( 1, NB ) )
                        ELSE
                           N1 = M
                           K1 = N
                           M1 = K
                           LW = MAX( 1, M1*MAX( 1, NB ) )
                        END IF
*
*                       Do first for TRANS = 'N', then for TRANS = 'T'
*
                        ITOFF = 0
                        DO 130 ITRAN = 1, 2
                           TRANS = TRANSS( ITRAN )
                           CALL DTIMMG( 0, M1, N1, B, LDA, 0, 0 )
                           IC = 0
                           S1 = DSECND( )
  110                      CONTINUE
                           CALL DORMQR( SIDE, TRANS, M1, N1, K1, A, LDA,
     $                                  TAU, B, LDA, WORK, LW, INFO )
                           S2 = DSECND( )
                           TIME = S2 - S1
                           IC = IC + 1
                           IF( TIME.LT.TIMMIN ) THEN
                              CALL DTIMMG( 0, M1, N1, B, LDA, 0, 0 )
                              GO TO 110
                           END IF
*
*                          Subtract the time used in DTIMMG.
*
                           ICL = 1
                           S1 = DSECND( )
  120                      CONTINUE
                           S2 = DSECND( )
                           UNTIME = S2 - S1
                           ICL = ICL + 1
                           IF( ICL.LE.IC ) THEN
                              CALL DTIMMG( 0, M1, N1, B, LDA, 0, 0 )
                              GO TO 120
                           END IF
*
                           TIME = ( TIME-UNTIME ) / DBLE( IC )
                           OPS = DOPLA( 'DORMQR', M1, N1, K1, ISIDE-1,
     $                           NB )
                           RESLTS( INB, IM, ILDA,
     $                        I4+ITOFF+IK ) = DMFLOP( OPS, TIME, INFO )
                           ITOFF = NK
  130                   CONTINUE
  140                CONTINUE
  150             CONTINUE
                  I4 = 2*NK
  160          CONTINUE
  170       CONTINUE
  180    CONTINUE
*
*        Print tables of results
*
         ISUB = 3
         I4 = 1
         IF( IMX.GE.1 ) THEN
            DO 220 ISIDE = 1, 2
               SIDE = SIDES( ISIDE )
               IF( ISIDE.EQ.1 ) THEN
                  WRITE( NOUT, FMT = 9998 )SUBNAM( ISUB )
                  IF( NLDA.GT.1 ) THEN
                     DO 190 I = 1, NLDA
                        WRITE( NOUT, FMT = 9997 )I, LDAVAL( I )
  190                CONTINUE
                  END IF
               END IF
               DO 210 ITRAN = 1, 2
                  TRANS = TRANSS( ITRAN )
                  DO 200 IK = 1, NK
                     IF( ISIDE.EQ.1 ) THEN
                        N = KVAL( IK )
                        WRITE( NOUT, FMT = 9995 )SUBNAM( ISUB ), SIDE,
     $                     TRANS, 'N', N
                        LABM = 'M'
                     ELSE
                        M = KVAL( IK )
                        WRITE( NOUT, FMT = 9995 )SUBNAM( ISUB ), SIDE,
     $                     TRANS, 'M', M
                        LABM = 'N'
                     END IF
                     CALL DPRTB5( 'NB', LABM, 'K', NNB, NBVAL, IMX,
     $                            MUSE, NUSE, NLDA,
     $                            RESLTS( 1, 1, 1, I4 ), LDR1, LDR2,
     $                            NOUT )
                     I4 = I4 + 1
  200             CONTINUE
  210          CONTINUE
  220       CONTINUE
         ELSE
            WRITE( NOUT, FMT = 9994 )SUBNAM( ISUB )
         END IF
      END IF
  230 CONTINUE
 9999 FORMAT( 1X, A6, ' timing run not attempted', / )
 9998 FORMAT( / ' *** Speed of ', A6, ' in megaflops ***' )
 9997 FORMAT( 5X, 'line ', I2, ' with LDA = ', I5 )
 9996 FORMAT( 5X, 'K = min(M,N)', / )
 9995 FORMAT( / 5X, A6, ' with SIDE = ''', A1, ''', TRANS = ''', A1,
     $      ''', ', A1, ' =', I6, / )
 9994 FORMAT( ' *** No pairs (M,N) found with M >= N:  ', A6,
     $      ' not timed' )
      RETURN
*
*     End of DTIMQR
*
      END
