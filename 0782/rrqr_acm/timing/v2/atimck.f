      SUBROUTINE ATIMCK( ICHK, SUBNAM, NN, NVAL, NLDA, LDAVAL, NOUT,
     $                   INFO )
*
*  -- LAPACK timing routine (version 2.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     March 31, 1993
*
*     .. Scalar Arguments ..
      CHARACTER*6        SUBNAM
      INTEGER            ICHK, INFO, NLDA, NN, NOUT
*     ..
*     .. Array Arguments ..
      INTEGER            LDAVAL( * ), NVAL( * )
*     ..
*
*  Purpose
*  =======
*
*  ATIMCK checks the input values of M, N, or K and LDA to determine
*  if they are valid for type TYPE.  The tests to be performed are
*  specified in the option variable ICHK.
*
*  On exit, INFO contains a count of the number of pairs (N,LDA) that
*  were invalid.
*
*  Arguments
*  =========
*
*  ICHK    (input) INTEGER
*          Specifies the type of comparison
*          = 1:  M <= LDA
*          = 2:  N <= LDA
*          = 3:  K <= LDA
*          = 4:  N*(N+1)/2 <= LA
*          = 0 or other value:  Determined from name passed in SUBNAM
*
*  SUBNAM  (input) CHARACTER*6
*          The name of the subroutine or path for which the input
*          values are to be tested.
*
*  NN      (input) INTEGER
*          The number of values of N contained in the vector NVAL.
*
*  NVAL    (input) INTEGER array, dimension( NN )
*          The values of the matrix size N.
*
*  NLDA    (input) INTEGER
*          The number of values of LDA contained in the vector LDAVAL.
*
*  LDAVAL  (input) INTEGER array, dimension( NLDA )
*          The values of the leading dimension of the array A.
*
*  NOUT    (input) INTEGER
*          The unit number for output.
*
*  INFO    (output) INTEGER
*          The number of pairs (N, LDA) that were invalid.
*
*  =====================================================================
*
*     .. Local Scalars ..
      CHARACTER*2        TYPE
      INTEGER            I, J, LDA, N
*     ..
*     .. External Functions ..
      LOGICAL            LSAMEN
      EXTERNAL           LSAMEN
*     ..
*     .. Executable Statements ..
*
      TYPE = SUBNAM( 2: 3 )
      INFO = 0
*
*     M, N, or K must be less than or equal to LDA.
*
      IF( ICHK.EQ.1 .OR. ICHK.EQ.2 .OR. ICHK.EQ.3 ) THEN
         DO 20 J = 1, NLDA
            LDA = LDAVAL( J )
            DO 10 I = 1, NN
               IF( NVAL( I ).GT.LDA ) THEN
                  INFO = INFO + 1
                  IF( NOUT.GT.0 ) THEN
                     IF( ICHK.EQ.1 ) THEN
                        WRITE( NOUT, FMT = 9999 )SUBNAM, NVAL( I ), LDA
                     ELSE IF( ICHK.EQ.2 ) THEN
                        WRITE( NOUT, FMT = 9998 )SUBNAM, NVAL( I ), LDA
                     ELSE
                        WRITE( NOUT, FMT = 9997 )SUBNAM, NVAL( I ), LDA
                     END IF
                  END IF
               END IF
   10       CONTINUE
   20    CONTINUE
*
*     IF TYPE = 'PP', 'SP', or 'HP',
*     then N*(N+1)/2 must be less than or equal to LA = LDAVAL(1).
*
      ELSE IF( ICHK.EQ.4 ) THEN
         LDA = LDAVAL( 1 )
         DO 30 I = 1, NN
            N = NVAL( I )
            IF( N*( N+1 ) / 2.GT.LDA ) THEN
               INFO = INFO + 1
               IF( NOUT.GT.0 )
     $            WRITE( NOUT, FMT = 9996 )SUBNAM, N, LDA
            END IF
   30    CONTINUE
*
*     IF TYPE = 'GB', then K must satisfy
*        2*K+1 <= LDA,  if SUBNAM = 'xGBMV'
*        3*K+1 <= LDA,  otherwise.
*
      ELSE IF( LSAMEN( 2, TYPE, 'GB' ) ) THEN
         IF( LSAMEN( 3, SUBNAM( 4: 6 ), 'MV ' ) ) THEN
            DO 50 J = 1, NLDA
               LDA = LDAVAL( J )
               DO 40 I = 1, NN
                  IF( 2*NVAL( I )+1.GT.LDA ) THEN
                     INFO = INFO + 1
                     IF( NOUT.GT.0 )
     $                  WRITE( NOUT, FMT = 9994 )SUBNAM, NVAL( I ),
     $                  LDA, 2*NVAL( I ) + 1
                  END IF
   40          CONTINUE
   50       CONTINUE
         ELSE
            DO 70 J = 1, NLDA
               LDA = LDAVAL( J )
               DO 60 I = 1, NN
                  IF( 3*NVAL( I )+1.GT.LDA ) THEN
                     INFO = INFO + 1
                     IF( NOUT.GT.0 )
     $                  WRITE( NOUT, FMT = 9995 )SUBNAM, NVAL( I ),
     $                  LDA, 3*NVAL( I ) + 1
                  END IF
   60          CONTINUE
   70       CONTINUE
         END IF
*
*     IF TYPE = 'PB' or 'TB', then K must satisfy
*        K+1 <= LDA.
*
      ELSE IF( LSAMEN( 2, TYPE, 'PB' ) .OR. LSAMEN( 2, TYPE, 'TB' ) )
     $          THEN
         DO 90 J = 1, NLDA
            LDA = LDAVAL( J )
            DO 80 I = 1, NN
               IF( NVAL( I )+1.GT.LDA ) THEN
                  INFO = INFO + 1
                  IF( NOUT.GT.0 )
     $               WRITE( NOUT, FMT = 9993 )SUBNAM, NVAL( I ), LDA
               END IF
   80       CONTINUE
   90    CONTINUE
*
*     IF TYPE = 'SB' or 'HB', then K must satisfy
*        K+1   <= LDA,  if SUBNAM = 'xxxMV '
*
      ELSE IF( LSAMEN( 2, TYPE, 'SB' ) .OR. LSAMEN( 2, TYPE, 'HB' ) )
     $          THEN
         IF( LSAMEN( 3, SUBNAM( 4: 6 ), 'MV ' ) ) THEN
            DO 110 J = 1, NLDA
               LDA = LDAVAL( J )
               DO 100 I = 1, NN
                  IF( NVAL( I )+1.GT.LDA ) THEN
                     INFO = INFO + 1
                     IF( NOUT.GT.0 )
     $                  WRITE( NOUT, FMT = 9992 )SUBNAM, NVAL( I ), LDA
                  END IF
  100          CONTINUE
  110       CONTINUE
         END IF
*
      END IF
 9999 FORMAT( ' *** Error for ', A6, ':  M > LDA for M =', I6,
     $      ', LDA =', I7 )
 9998 FORMAT( ' *** Error for ', A6, ':  N > LDA for N =', I6,
     $      ', LDA =', I7 )
 9997 FORMAT( ' *** Error for ', A6, ':  K > LDA for K =', I6,
     $      ', LDA =', I7 )
 9996 FORMAT( ' *** Error for ', A6, ':  N*(N+1)/2 > LA for N =', I6,
     $      ', LA =', I7 )
 9995 FORMAT( ' *** Error for ', A6, ':  3*K+1 > LDA for K =', I6,
     $      ', LDA =', I7, / ' --> Increase LDA to at least ', I7 )
 9994 FORMAT( ' *** Error for ', A6, ':  2*K+1 > LDA for K =', I6,
     $      ', LDA =', I7, / ' --> Increase LDA to at least ', I7 )
 9993 FORMAT( ' *** Error for ', A6, ':  K+1 > LDA for K =', I6, ', LD',
     $      'A =', I7 )
 9992 FORMAT( ' *** Error for ', A6, ':  2*K+2 > LDA for K =', I6, ', ',
     $      'LDA =', I7 )
*
      RETURN
*
*     End of ATIMCK
*
      END
