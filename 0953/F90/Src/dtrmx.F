***********************************************************************
*                                                                     *
*     dtrmx.F:                                                        *
*         Auxiliary routine in the package PDHSEQR.                   *
*                                                                     *
*     Contributors: Robert Granat                                     *
*                   Bo Kagstrom                                       *
*                   Daniel Kressner                                   *
*                   Meiyue Shao                                       *
*                                                                     *
*     Department of Computing Science and HPC2N, Umea University      *
*     MATHICSE ANCHP, EPF Lausanne                                    *
*                                                                     *
***********************************************************************
      SUBROUTINE DTRMX ( SIDE, UPLO, TRANSA, DIAG, M, N, ALPHA, A, LDA,
     $                   B, LDB )
*     .. Scalar Arguments ..
      CHARACTER*1        SIDE, UPLO, TRANSA, DIAG
      INTEGER            M, N, LDA, LDB
      DOUBLE PRECISION   ALPHA
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), B( LDB, * )
*     ..
*
*  Purpose
*  =======
*
*  DTRMX  performs one of the matrix-matrix operations
*
*     B := alpha*op( A )*B,   or   B := alpha*B*op( A ),
*
*  where  alpha  is a scalar,  B  is an m by n matrix,  A  is a unit, or
*  non-unit,  upper or lower triangular matrix  and  op( A )  is one  of
*
*     op( A ) = A   or   op( A ) = A'.
*
*
*     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            OMP_GET_NUM_THREADS, OMP_GET_THREAD_NUM
      EXTERNAL           LSAME, OMP_GET_NUM_THREADS, OMP_GET_THREAD_NUM
*     .. External Subroutines ..
      EXTERNAL           XERBLA
*     .. Intrinsic Functions ..
      INTRINSIC          MAX
*     .. Local Scalars ..
      LOGICAL            LSIDE, NOUNIT, UPPER
      INTEGER            I, INFO, J, K, NROWA
      DOUBLE PRECISION   TEMP
      INTEGER            INDX, THREADS, MCHUNK, M2, N2, MYTHREAD, NCHUNK
      LOGICAL            NOTA
*     .. Parameters ..
      DOUBLE PRECISION   ONE         , ZERO
      PARAMETER        ( ONE = 1.0D+0, ZERO = 0.0D+0 )
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      NOTA = LSAME( TRANSA, 'N' )
      LSIDE  = LSAME( SIDE  , 'L' )
      IF( LSIDE )THEN
         NROWA = M
      ELSE
         NROWA = N
      END IF
      NOUNIT = LSAME( DIAG  , 'N' )
      UPPER  = LSAME( UPLO  , 'U' )
*
      INFO   = 0
      IF(      ( .NOT.LSIDE                ).AND.
     $         ( .NOT.LSAME( SIDE  , 'R' ) )      )THEN
         INFO = 1
      ELSE IF( ( .NOT.UPPER                ).AND.
     $         ( .NOT.LSAME( UPLO  , 'L' ) )      )THEN
         INFO = 2
      ELSE IF( ( .NOT.LSAME( TRANSA, 'N' ) ).AND.
     $         ( .NOT.LSAME( TRANSA, 'T' ) ).AND.
     $         ( .NOT.LSAME( TRANSA, 'C' ) )      )THEN
         INFO = 3
      ELSE IF( ( .NOT.LSAME( DIAG  , 'U' ) ).AND.
     $         ( .NOT.LSAME( DIAG  , 'N' ) )      )THEN
         INFO = 4
      ELSE IF( M  .LT.0               )THEN
         INFO = 5
      ELSE IF( N  .LT.0               )THEN
         INFO = 6
      ELSE IF( LDA.LT.MAX( 1, NROWA ) )THEN
         INFO = 9
      ELSE IF( LDB.LT.MAX( 1, M     ) )THEN
         INFO = 11
      END IF
      IF( INFO.NE.0 )THEN
         CALL XERBLA( 'DTRMX ', INFO )
         RETURN
      END IF
*
*
*
*
#ifdef USE_OMP
*$OMP PARALLEL SHARED(THREADS)
*$OMP MASTER
      THREADS = OMP_GET_NUM_THREADS()
*$OMP END MASTER
*$OMP END PARALLEL
      MCHUNK = MAX( 1, M / THREADS )
      NCHUNK = MAX( 1, N / THREADS )
#else
      THREADS = 1
      MCHUNK = M
      NCHUNK = N
#endif
*
*     Do the real work
*
#ifdef USE_OMP
      IF( LSIDE ) THEN
*$OMP PARALLEL DEFAULT( NONE ),
*$OMP& SHARED(SIDE,UPLO,TRANSA,DIAG,M,N,ALPHA,A,LDA,B,LDB,THREADS,
*$OMP& NCHUNK,NOTA), PRIVATE( INDX, N2, MYTHREAD )
         MYTHREAD = OMP_GET_THREAD_NUM()
*$OMP DO
         DO 10 INDX = 1, N, NCHUNK
            N2 = MIN( NCHUNK, N-INDX+1 )
            CALL DTRMM( SIDE, UPLO, TRANSA, DIAG, M, N2, ALPHA, A, LDA,
     $                  B(1,INDX), LDB )
 10      CONTINUE
*$OMP END DO
*$OMP END PARALLEL
      ELSE
*$OMP PARALLEL DEFAULT( NONE ),
*$OMP& SHARED(SIDE,UPLO,TRANSA,DIAG,M,N,ALPHA,A,LDA,B,LDB,THREADS,
*$OMP& MCHUNK,NOTA), PRIVATE( INDX, M2, MYTHREAD )
         MYTHREAD = OMP_GET_THREAD_NUM()
*$OMP DO
         DO 20 INDX = 1, M, MCHUNK
            M2 = MIN( MCHUNK, M-INDX+1 )
            CALL DTRMM( SIDE, UPLO, TRANSA, DIAG, M2, N, ALPHA, A, LDA,
     $                  B(INDX,1), LDB )
 20      CONTINUE
*$OMP END DO
*$OMP END PARALLEL
      END IF
#else
      CALL DTRMM( SIDE, UPLO, TRANSA, DIAG, M, N, ALPHA, A, LDA, B,
     $            LDB )
#endif
*
      RETURN
*
*     End of DTRMX
*
      END

