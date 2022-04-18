      SUBROUTINE DGEMX(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
*
      IMPLICIT NONE
*
*     .. Scalar Arguments ..
      DOUBLE PRECISION ALPHA,BETA
      INTEGER K,LDA,LDB,LDC,M,N
      CHARACTER TRANSA,TRANSB
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION A(LDA,*),B(LDB,*),C(LDC,*)
*     ..
*
*  Purpose
*  =======
*
*  DGEMX  performs one of the matrix-matrix operations
*
*     C := alpha*op( A )*op( B ) + beta*C,
*
*  where  op( X ) is one of
*
*     op( X ) = X   or   op( X ) = X',
*
*  alpha and beta are scalars, and A, B and C are matrices, with op( A )
*  an m by k matrix,  op( B )  a  k by n matrix and  C an m by n matrix.
*
*     .. Parameters ..
*
*     .. External Functions ..
      LOGICAL LSAME
      INTEGER OMP_GET_NUM_THREADS, OMP_GET_THREAD_NUM
      EXTERNAL LSAME, OMP_GET_NUM_THREADS, OMP_GET_THREAD_NUM
*     ..
*     .. External Subroutines ..
      EXTERNAL XERBLA, DGEMM
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC MAX
*     ..
*     .. Local Scalars ..
      DOUBLE PRECISION TEMP
      INTEGER I,INFO,J,L,NCOLA,NROWA,NROWB, THREADS, MCHUNK, INDX, M2, 
     $     N2, MYTHREAD, NCHUNK
      LOGICAL NOTA,NOTB
*     ..
*     .. Parameters ..
      DOUBLE PRECISION ONE,ZERO
      PARAMETER (ONE=1.0D+0,ZERO=0.0D+0)
*     ..
*
*     Set  NOTA  and  NOTB  as  true if  A  and  B  respectively are not
*     transposed and set  NROWA, NCOLA and  NROWB  as the number of rows
*     and  columns of  A  and the  number of  rows  of  B  respectively.
*
      NOTA = LSAME(TRANSA,'N')
      NOTB = LSAME(TRANSB,'N')
      IF (NOTA) THEN
          NROWA = M
          NCOLA = K
      ELSE
          NROWA = K
          NCOLA = M
      END IF
      IF (NOTB) THEN
          NROWB = K
      ELSE
          NROWB = N
      END IF
*
*     Test the input parameters.
*
      INFO = 0
      IF ((.NOT.NOTA) .AND. (.NOT.LSAME(TRANSA,'C')) .AND.
     $    (.NOT.LSAME(TRANSA,'T'))) THEN
          INFO = 1
      ELSE IF ((.NOT.NOTB) .AND. (.NOT.LSAME(TRANSB,'C')) .AND.
     $         (.NOT.LSAME(TRANSB,'T'))) THEN
          INFO = 2
      ELSE IF (M.LT.0) THEN
          INFO = 3
      ELSE IF (N.LT.0) THEN
          INFO = 4
      ELSE IF (K.LT.0) THEN
          INFO = 5
      ELSE IF (LDA.LT.MAX(1,NROWA)) THEN
          INFO = 8
      ELSE IF (LDB.LT.MAX(1,NROWB)) THEN
          INFO = 10
      ELSE IF (LDC.LT.MAX(1,M)) THEN
          INFO = 13
      END IF
      IF (INFO.NE.0) THEN
          CALL XERBLA('DGEMX ',INFO)
          STOP
      END IF
*
*     Quick return if possible.
*
      IF ((M.EQ.0) .OR. (N.EQ.0) .OR.
     $    (((ALPHA.EQ.ZERO).OR. (K.EQ.0)).AND. (BETA.EQ.ONE))) RETURN
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
      IF( M.GT.N ) THEN
*$OMP PARALLEL DEFAULT( NONE ),
*$OMP& SHARED(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC,THREADS,
*$OMP& MCHUNK,NOTA,NOTB), PRIVATE( INDX, M2, MYTHREAD )
         MYTHREAD = OMP_GET_THREAD_NUM()
*$OMP DO
         DO 10 INDX = 1, M, MCHUNK
            M2 = MIN( MCHUNK, M-INDX+1 )
            IF( NOTA ) THEN
               CALL DGEMM(TRANSA,TRANSB,M2,N,K,ALPHA,A(INDX,1),LDA,
     $              B,LDB,BETA,C(INDX,1),LDC)
            ELSE
               CALL DGEMM(TRANSA,TRANSB,M2,N,K,ALPHA,A(1,INDX),LDA,
     $              B,LDB,BETA,C(INDX,1),LDC)
            END IF 
 10      CONTINUE
*$OMP END DO  
*$OMP END PARALLEL   
      ELSE
*$OMP PARALLEL DEFAULT( NONE ),
*$OMP& SHARED(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC,THREADS,
*$OMP& NCHUNK,NOTA,NOTB), PRIVATE( INDX, N2, MYTHREAD )
         MYTHREAD = OMP_GET_THREAD_NUM()   
*$OMP DO
         DO 20 INDX = 1, N, NCHUNK
            N2 = MIN( NCHUNK, N-INDX+1 )
            IF( NOTB ) THEN
               CALL DGEMM(TRANSA,TRANSB,M,N2,K,ALPHA,A,LDA,
     $              B(1,INDX),LDB,BETA,C(1,INDX),LDC)
            ELSE
               CALL DGEMM(TRANSA,TRANSB,M,N2,K,ALPHA,A,LDA,
     $              B(INDX,1),LDB,BETA,C(1,INDX),LDC)
            END IF 
 20      CONTINUE
*$OMP END DO       
*$OMP END PARALLEL
      END IF 
#else
      CALL DGEMM(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
#endif
*
      RETURN
*
*     End of DGEMX
*
      END
