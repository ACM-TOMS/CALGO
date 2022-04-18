      SUBROUTINE STIMMV( VNAME, NN, NVAL, NK, KVAL, NLDA, LDAVAL,
     $                   TIMMIN, A, LB, B, C, RESLTS, LDR1, LDR2, NOUT )
*
*  -- LAPACK timing routine (version 2.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     March 31, 1993
*
*     .. Scalar Arguments ..
      CHARACTER*( * )    VNAME
      INTEGER            LB, LDR1, LDR2, NK, NLDA, NN, NOUT
      REAL               TIMMIN
*     ..
*     .. Array Arguments ..
      INTEGER            KVAL( * ), LDAVAL( * ), NVAL( * )
      REAL               A( * ), B( * ), C( * ), RESLTS( LDR1, LDR2, * )
*     ..
*
*  Purpose
*  =======
*
*  STIMMV times individual BLAS 2 routines.
*
*  Arguments
*  =========
*
*  VNAME   (input) CHARACTER*(*)
*          The name of the Level 2 BLAS routine to be timed.
*
*  NN      (input) INTEGER
*          The number of values of N contained in the vector NVAL.
*
*  NVAL    (input) INTEGER array, dimension (NN)
*          The values of the matrix dimension N.
*
*  NK      (input) INTEGER
*          The number of values of K contained in the vector KVAL.
*
*  KVAL    (input) INTEGER array, dimension (NK)
*          The values of the bandwidth K.
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
*             where LDAMAX and NMAX are the maximum values permitted
*             for LDA and N.
*
*  LB      (input) INTEGER
*          The length of B and C, needed when timing SGBMV.  If timing
*          SGEMV, LB >= LDAMAX*NMAX.
*
*  B       (workspace) REAL array, dimension (LB)
*
*  C       (workspace) REAL array, dimension (LB)
*
*  RESLTS  (output) REAL array, dimension (LDR1,LDR2,NLDA)
*          The timing results for each subroutine over the relevant
*          values of N and LDA.
*
*  LDR1    (input) INTEGER
*          The first dimension of RESLTS.  LDR1 >= max(1,NK).
*
*  LDR2    (input) INTEGER
*          The second dimension of RESLTS.  LDR2 >= max(1,NN).
*
*  NOUT    (input) INTEGER
*          The unit number for output.
*
*  =====================================================================
*
*     .. Parameters ..
      INTEGER            NSUBS
      REAL               ONE
      PARAMETER          ( NSUBS = 2, ONE = 1.0E+0 )
*     ..
*     .. Local Scalars ..
      CHARACTER          LAB1, LAB2
      CHARACTER*6        CNAME
      INTEGER            I, IB, IC, ICL, IK, ILDA, IN, INFO, ISUB, K,
     $                   KL, KU, LDA, LDB, N, NRHS
      REAL               OPS, S1, S2, TIME, UNTIME
*     ..
*     .. Local Arrays ..
      LOGICAL            TIMSUB( NSUBS )
      CHARACTER*6        SUBNAM( NSUBS )
*     ..
*     .. External Functions ..
      LOGICAL            LSAME, LSAMEN
      REAL               SECOND, SMFLOP, SOPBL2
      EXTERNAL           LSAME, LSAMEN, SECOND, SMFLOP, SOPBL2
*     ..
*     .. External Subroutines ..
      EXTERNAL           ATIMCK, SGBMV, SGEMV, SPRTBL, STIMMG
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN, REAL
*     ..
*     .. Data statements ..
      DATA               SUBNAM / 'SGEMV ', 'SGBMV ' /
*     ..
*     .. Executable Statements ..
*
      CNAME = VNAME
      DO 10 ISUB = 1, NSUBS
         TIMSUB( ISUB ) = LSAMEN( 6, CNAME, SUBNAM( ISUB ) )
         IF( TIMSUB( ISUB ) )
     $      GO TO 20
   10 CONTINUE
      WRITE( NOUT, FMT = 9999 )CNAME
      GO TO 150
   20 CONTINUE
*
*     Check that N or K <= LDA for the input values.
*
      IF( LSAME( CNAME( 3: 3 ), 'B' ) ) THEN
         CALL ATIMCK( 0, CNAME, NK, KVAL, NLDA, LDAVAL, NOUT, INFO )
         LAB1 = 'M'
         LAB2 = 'K'
      ELSE
         CALL ATIMCK( 2, CNAME, NN, NVAL, NLDA, LDAVAL, NOUT, INFO )
         LAB1 = ' '
         LAB2 = 'N'
      END IF
      IF( INFO.GT.0 ) THEN
         WRITE( NOUT, FMT = 9998 )CNAME
         GO TO 150
      END IF
*
*     Print the table header on unit NOUT.
*
      WRITE( NOUT, FMT = 9997 )VNAME
      IF( NLDA.EQ.1 ) THEN
         WRITE( NOUT, FMT = 9996 )LDAVAL( 1 )
      ELSE
         DO 30 I = 1, NLDA
            WRITE( NOUT, FMT = 9995 )I, LDAVAL( I )
   30    CONTINUE
      END IF
      WRITE( NOUT, FMT = * )
*
*     Time SGEMV
*
      IF( TIMSUB( 1 ) ) THEN
         DO 80 ILDA = 1, NLDA
            LDA = LDAVAL( ILDA )
            DO 70 IN = 1, NN
               N = NVAL( IN )
               NRHS = N
               LDB = LDA
               CALL STIMMG( 1, N, N, A, LDA, 0, 0 )
               CALL STIMMG( 0, N, NRHS, B, LDB, 0, 0 )
               CALL STIMMG( 1, N, NRHS, C, LDB, 0, 0 )
               IC = 0
               S1 = SECOND( )
   40          CONTINUE
               IB = 1
               DO 50 I = 1, NRHS
                  CALL SGEMV( 'No transpose', N, N, ONE, A, LDA,
     $                        B( IB ), 1, ONE, C( IB ), 1 )
                  IB = IB + LDB
   50          CONTINUE
               S2 = SECOND( )
               TIME = S2 - S1
               IC = IC + 1
               IF( TIME.LT.TIMMIN ) THEN
                  CALL STIMMG( 1, N, NRHS, C, LDB, 0, 0 )
                  GO TO 40
               END IF
*
*              Subtract the time used in STIMMG.
*
               ICL = 1
               S1 = SECOND( )
   60          CONTINUE
               S2 = SECOND( )
               UNTIME = S2 - S1
               ICL = ICL + 1
               IF( ICL.LE.IC ) THEN
                  CALL STIMMG( 1, N, NRHS, C, LDB, 0, 0 )
                  GO TO 60
               END IF
*
               TIME = ( TIME-UNTIME ) / REAL( IC )
               OPS = NRHS*SOPBL2( 'SGEMV ', N, N, 0, 0 )
               RESLTS( 1, IN, ILDA ) = SMFLOP( OPS, TIME, 0 )
   70       CONTINUE
   80    CONTINUE
*
         CALL SPRTBL( LAB1, LAB2, 1, NVAL, NN, NVAL, NLDA, RESLTS, LDR1,
     $                LDR2, NOUT )
*
      ELSE IF( TIMSUB( 2 ) ) THEN
*
*        Time SGBMV
*
         DO 140 ILDA = 1, NLDA
            LDA = LDAVAL( ILDA )
            DO 130 IN = 1, NN
               N = NVAL( IN )
               DO 120 IK = 1, NK
                  K = MIN( N-1, MAX( 0, KVAL( IK ) ) )
                  KL = K
                  KU = K
                  LDB = N
                  CALL STIMMG( 2, N, N, A, LDA, KL, KU )
                  NRHS = MIN( K, LB / LDB )
                  CALL STIMMG( 0, N, NRHS, B, LDB, 0, 0 )
                  CALL STIMMG( 1, N, NRHS, C, LDB, 0, 0 )
                  IC = 0
                  S1 = SECOND( )
   90             CONTINUE
                  IB = 1
                  DO 100 I = 1, NRHS
                     CALL SGBMV( 'No transpose', N, N, KL, KU, ONE,
     $                           A( KU+1 ), LDA, B( IB ), 1, ONE,
     $                           C( IB ), 1 )
                     IB = IB + LDB
  100             CONTINUE
                  S2 = SECOND( )
                  TIME = S2 - S1
                  IC = IC + 1
                  IF( TIME.LT.TIMMIN ) THEN
                     CALL STIMMG( 1, N, NRHS, C, LDB, 0, 0 )
                     GO TO 90
                  END IF
*
*                 Subtract the time used in STIMMG.
*
                  ICL = 1
                  S1 = SECOND( )
  110             CONTINUE
                  S2 = SECOND( )
                  UNTIME = S2 - S1
                  ICL = ICL + 1
                  IF( ICL.LE.IC ) THEN
                     CALL STIMMG( 1, N, NRHS, C, LDB, 0, 0 )
                     GO TO 110
                  END IF
*
                  TIME = ( TIME-UNTIME ) / REAL( IC )
                  OPS = NRHS*SOPBL2( 'SGBMV ', N, N, KL, KU )
                  RESLTS( IN, IK, ILDA ) = SMFLOP( OPS, TIME, 0 )
  120          CONTINUE
  130       CONTINUE
  140    CONTINUE
*
         CALL SPRTBL( LAB1, LAB2, NN, NVAL, NK, KVAL, NLDA, RESLTS,
     $                LDR1, LDR2, NOUT )
      END IF
*
  150 CONTINUE
 9999 FORMAT( 1X, A6, ':  Unrecognized path or subroutine name', / )
 9998 FORMAT( 1X, A6, ' timing run not attempted', / )
 9997 FORMAT( / ' *** Speed of ', A6, ' in megaflops ***' )
 9996 FORMAT( 5X, 'with LDA = ', I5 )
 9995 FORMAT( 5X, 'line ', I2, ' with LDA = ', I5 )
      RETURN
*
*     End of STIMMV
*
      END
