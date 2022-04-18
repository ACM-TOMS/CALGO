      SUBROUTINE ZGB05 ( UPLO, TRANS, N, K, ALPHA, A, LDA,
     $                   BETA, C, LDC )
*     .. Scalar Arguments ..
      CHARACTER*1        UPLO, TRANS
      INTEGER            N, K, LDA, LDC
      DOUBLE PRECISION   ALPHA, BETA
*     .. Array Arguments ..
      COMPLEX*16         A( LDA, * ), C( LDC, * )
*     ..
*
*  Purpose
*  =======
*
*  ZGB05 (ZHERK)  performs one of the hermitian rank k operations
*
*     C := alpha*A*conjg( A' ) + beta*C,
*
*  or
*
*     C := alpha*conjg( A' )*A + beta*C,
*
*  where  alpha and beta  are  real scalars,  C is an  n by n  hermitian
*  matrix and  A  is an  n by k  matrix in the  first case and a  k by n
*  matrix in the second case.
*
*  Parameters
*  ==========
*
*  UPLO   - CHARACTER*1.
*           On  entry,   UPLO  specifies  whether  the  upper  or  lower
*           triangular  part  of the  array  C  is to be  referenced  as
*           follows:
*
*              UPLO = 'U' or 'u'   Only the  upper triangular part of  C
*                                  is to be referenced.
*
*              UPLO = 'L' or 'l'   Only the  lower triangular part of  C
*                                  is to be referenced.
*
*           Unchanged on exit.
*
*  TRANS  - CHARACTER*1.
*           On entry,  TRANS  specifies the operation to be performed as
*           follows:
*
*              TRANS = 'N' or 'n'   C := alpha*A*conjg( A' ) + beta*C.
*
*              TRANS = 'C' or 'c'   C := alpha*conjg( A' )*A + beta*C.
*
*           Unchanged on exit.
*
*  N      - INTEGER.
*           On entry,  N specifies the order of the matrix C.  N must be
*           at least zero.
*           Unchanged on exit.
*
*  K      - INTEGER.
*           On entry with  TRANS = 'N' or 'n',  K  specifies  the number
*           of  columns   of  the   matrix   A,   and  on   entry   with
*           TRANS = 'C' or 'c',  K  specifies  the number of rows of the
*           matrix A.  K must be at least zero.
*           Unchanged on exit.
*
*  ALPHA  - DOUBLE PRECISION.
*           On entry, ALPHA specifies the scalar alpha.
*           Unchanged on exit.
*
*  A      - COMPLEX*16       array of DIMENSION ( LDA, ka ), where ka is
*           k  when  TRANS = 'N' or 'n',  and is  n  otherwise.
*           Before entry with  TRANS = 'N' or 'n',  the  leading  n by k
*           part of the array  A  must contain the matrix  A,  otherwise
*           the leading  k by n  part of the array  A  must contain  the
*           matrix A.
*           Unchanged on exit.
*
*  LDA    - INTEGER.
*           On entry, LDA specifies the first dimension of A as declared
*           in  the  calling  (sub)  program.   When  TRANS = 'N' or 'n'
*           then  LDA must be at least  max( 1, n ), otherwise  LDA must
*           be at least  max( 1, k ).
*           Unchanged on exit.
*
*  BETA   - DOUBLE PRECISION.
*           On entry, BETA specifies the scalar beta.
*           Unchanged on exit.
*
*  C      - COMPLEX*16       array of DIMENSION ( LDC, n ).
*           Before entry  with  UPLO = 'U' or 'u',  the leading  n by n
*           upper triangular part of the array C must contain the upper
*           triangular part  of the  hermitian matrix  and the strictly
*           lower triangular part of C is not referenced.  On exit, the
*           upper triangular part of the array  C is overwritten by the
*           upper triangular part of the updated matrix.
*           Before entry  with  UPLO = 'L' or 'l',  the leading  n by n
*           lower triangular part of the array C must contain the lower
*           triangular part  of the  hermitian matrix  and the strictly
*           upper triangular part of C is not referenced.  On exit, the
*           lower triangular part of the array  C is overwritten by the
*           lower triangular part of the updated matrix.
*           Note that the imaginary parts of the diagonal elements need
*           not be set,  they are assumed to be zero,  and on exit they
*           are set to zero.
*
*  LDC    - INTEGER.
*           On entry, LDC specifies the first dimension of C as declared
*           in  the  calling  (sub)  program.   LDC  must  be  at  least
*           max( 1, n ).
*           Unchanged on exit.
*
*
*  Level 3 Blas routine.
*
*  -- Written on 8-February-1989.
*     Jack Dongarra, Argonne National Laboratory.
*     Iain Duff, AERE Harwell.
*     Jeremy Du Croz, Numerical Algorithms Group Ltd.
*     Sven Hammarling, Numerical Algorithms Group Ltd.
*
*  -- Rewritten in May-1994.
*     GEMM-Based Level 3 BLAS.
*     Per Ling, Institute of Information Processing,
*     University of Umea, Sweden.
*
*
*     .. Local Scalars ..
      INTEGER            INFO, NROWA
      INTEGER            I, II, IX, ISEC, L, LL, LSEC
      LOGICAL            UPPER, NOTR, CLDA, SMALLN, TINYK
      COMPLEX*16         ZALPHA, ZBETA, ZDELTA
*     .. Intrinsic Functions ..
      INTRINSIC          MIN, MAX, DBLE, DCMPLX, DCONJG
*     .. External Functions ..
      LOGICAL            LSAME, ZGB90, ZGB91
      EXTERNAL           LSAME, ZGB90, ZGB91
*     .. External Subroutines ..
      EXTERNAL           XERBLA
      EXTERNAL           ZGEMM, ZGEMV, ZHER, ZCOPY, ZSCAL
*     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      COMPLEX*16         ZONE
      INTEGER            ZIP51, ZIP52
      PARAMETER        ( ONE = 1.0D+0, ZERO = 0.0D+0,
     $                   ZONE = ( 1.0D+0, 0.0D+0 ),
     $                   ZIP51 = 51, ZIP52 = 52 )
*     .. User specified parameters for ZGB05 ..
      INTEGER            RB, CB, RCB
      PARAMETER        ( RCB = 44, RB = 44, CB = 44 )
*     .. Local Arrays ..
      COMPLEX*16         T1( RB, CB ), T2( RCB, RCB ), T3( RCB, RCB ),
     $                   T4( CB )
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      UPPER = LSAME( UPLO, 'U' )
      NOTR = LSAME( TRANS, 'N' )
      IF( NOTR )THEN
         NROWA = N
      ELSE
         NROWA = K
      END IF
      INFO = 0
      IF( ( .NOT.UPPER ).AND.( .NOT.LSAME( UPLO , 'L' ) ) )THEN
         INFO = 1
      ELSE IF( ( .NOT.NOTR ).AND.( .NOT.LSAME( TRANS, 'C' ) ) )THEN
         INFO = 2
      ELSE IF( N.LT.0 )THEN
         INFO = 3
      ELSE IF( K.LT.0 )THEN
         INFO = 4
      ELSE IF( LDA.LT.MAX( 1, NROWA ) )THEN
         INFO = 7
      ELSE IF( LDC.LT.MAX( 1, N ) )THEN
         INFO = 10
      END IF
      IF( INFO.NE.0 )THEN
         CALL XERBLA( 'ZGB05 ', INFO )
         RETURN
      END IF
*
*     Quick return if possible.
*
      IF( ( N.EQ.0 ).OR.
     $    ( ( ( ALPHA.EQ.ZERO ).OR.( K.EQ.0 ) ).AND.( BETA.EQ.ONE ) ) )
     $   RETURN
*
      ZALPHA = DCMPLX( ALPHA, ZERO )
      ZBETA = DCMPLX( BETA, ZERO )
*
*     And when alpha.eq.zero or k.eq.0.
*
      IF( ( ALPHA.EQ.ZERO ).OR.( K.EQ.0 ) )THEN
         IF( UPPER )THEN
            C( 1, 1 ) = DCMPLX( BETA*DBLE( C( 1, 1 ) ), ZERO )
            DO 10, I = 2, N
               CALL ZSCAL ( I-1, ZBETA, C( 1, I ), 1 )
               C( I, I ) = DCMPLX( BETA*DBLE( C( I, I ) ), ZERO )
   10       CONTINUE
         ELSE
            DO 20, I = 1, N-1
               C( I, I ) = DCMPLX( BETA*DBLE( C( I, I ) ), ZERO )
               CALL ZSCAL ( N-I, ZBETA, C( I+1, I ), 1 )
   20       CONTINUE
            C( N, N ) = DCMPLX( BETA*DBLE( C( N, N ) ), ZERO )
         END IF
         RETURN
      END IF
*
*     Start the operations.
*
      IF( UPPER )THEN
         IF( NOTR )THEN
*
*           Form  C := alpha*A*conjg( A' ) + beta*C. Upper, Notr.
*
            SMALLN = .NOT.ZGB90( ZIP51 , N, K )
            IF( SMALLN )THEN
               TINYK = .NOT.ZGB90( ZIP52 , N, K )
               DO 90, II = 1, N, RCB
                  ISEC = MIN( RCB, N-II+1 )
*
*                 C := alpha*A*conjg( A' ) + beta*C, matrix multiply
*                 updating upper vertical blocks of C.
*
                  IF( II.GT.1 )THEN
                     CALL ZGEMM ( 'N', 'C', II-1, ISEC, K, ZALPHA,
     $                                 A( 1, 1 ), LDA, A( II, 1 ), LDA,
     $                                         ZBETA, C( 1, II ), LDC )
                  END IF
                  IF( TINYK )THEN
*
*                    C :=  beta*C, a upper triangular diagonal block
*                    of C is updated with beta. The imaginary part of
*                    the diagonal elements of C are set to ZERO.
*
                     IF( BETA.NE.ONE )THEN
                        C( II, II ) =
     $                         DCMPLX( BETA*DBLE( C( II, II ) ), ZERO )
                        DO 30, I = II+1, II+ISEC-1
                           CALL ZSCAL ( I-II, ZBETA, C( II, I ), 1 )
                           C( I, I ) =
     $                           DCMPLX( BETA*DBLE( C( I, I ) ), ZERO )
   30                   CONTINUE
                     END IF
*
*                    C := alpha*A*conjg( A' ) + C, hermitian matrix
*                    multiply. C is a hermitian diagonal block having
*                    upper triangular storage format.
*
                     DO 40, L = 1, K
                        CALL ZHER  ( 'U', ISEC, ALPHA, A( II, L ),
     $                                            1, C( II, II ), LDC )
   40                CONTINUE
                  ELSE
*
*                    T2 := C, a upper triangular diagonal block of the
*                    hermitian matrix C is copied to the upper
*                    triangular part of T2.
*
                     DO 50, I = II, II+ISEC-1
                        CALL ZCOPY ( I-II+1, C( II, I ), 1,
     $                                             T2( 1, I-II+1 ), 1 )
   50                CONTINUE
*
*                    T2 :=  beta*T2, the upper triangular part of T2 is
*                    updated with beta. The imaginary part of the
*                    diagonal elements of T2 are set to ZERO.
*
                     IF( BETA.NE.ONE )THEN
                        T2( 1, 1 ) =
     $                          DCMPLX( BETA*DBLE( T2( 1, 1 ) ), ZERO )
                        DO 60, I = 2, ISEC
                           CALL ZSCAL ( I-1, ZBETA, T2( 1, I ), 1 )
                           T2( I, I ) =
     $                          DCMPLX( BETA*DBLE( T2( I, I ) ), ZERO )
   60                   CONTINUE
                     END IF
*
*                    T2 := alpha*A*conjg( A' ) + T2, hermitian matrix
*                    multiply. T2 contains a hermitian block having
*                    upper triangular storage format.
*
                     DO 70, L = 1, K
                        CALL ZHER  ( 'U', ISEC, ALPHA, A( II, L ),
     $                                             1, T2( 1, 1 ), RCB )
   70                CONTINUE
*
*                    C := T2, the upper triangular part of T2 is copied
*                    back to C.
*
                     DO 80, I = II, II+ISEC-1
                        CALL ZCOPY ( I-II+1, T2( 1, I-II+1 ), 1,
     $                                                  C( II, I ), 1 )
   80                CONTINUE
                  END IF
   90          CONTINUE
            ELSE
               DO 140, II = 1, N, RB
                  ISEC = MIN( RB, N-II+1 )
*
*                 C := alpha*A*conjg( A' ) + beta*C, matrix multiply
*                 updating upper vertical blocks of C.
*
                  IF( II.GT.1 )THEN
                     CALL ZGEMM ( 'N', 'C', II-1, ISEC, K, ZALPHA,
     $                                 A( 1, 1 ), LDA, A( II, 1 ), LDA,
     $                                         ZBETA, C( 1, II ), LDC )
                  END IF
                  ZDELTA = ZBETA
                  DO 130, LL = 1, K, CB
                     LSEC = MIN( CB, K-LL+1 )
*
*                    T1 := A, a rectangular block of A is copied to T1.
*
                     DO 100, L = LL, LL+LSEC-1
                        CALL ZCOPY ( ISEC, A( II, L ), 1,
     $                                             T1( 1, L-LL+1 ), 1 )
  100                CONTINUE
*
*                    C := alpha*T1*conjg( T1' ) + delta*C, C is
*                    hermitian having triangular storage format. Delta
*                    is used instead of beta to avoid updating the
*                    block of C with beta multiple times. The local
*                    array T4 is used for the conjugated transpose
*                    of vectors of T1.
*
                     DO 120, I = II, II+ISEC-1
                        DO 110, L = LL, LL+LSEC-1
                           T4( L-LL+1 ) =
     $                                   DCONJG( T1( I-II+1, L-LL+1 ) )
  110                   CONTINUE
                        CALL ZGEMV ( 'N', I-II+1, LSEC, ZALPHA,
     $                                      T1( 1, 1 ), RB, T4( 1 ), 1,
     $                                          ZDELTA, C( II, I ), 1 )
                        C( I, I ) = DCMPLX( DBLE( C( I, I ) ), ZERO )
  120                CONTINUE
                     ZDELTA = ZONE
  130             CONTINUE
  140          CONTINUE
            END IF
         ELSE
*
*           Form  C := alpha*conjg( A' )*A + beta*C. Upper, Trans.
*
            SMALLN = .NOT.ZGB90( ZIP51 , N, K )
            IF( SMALLN )THEN
               TINYK = .NOT.ZGB90( ZIP52 , N, K )
               DO 250, II = 1, N, RCB
                  ISEC = MIN( RCB, N-II+1 )
*
*                 C := alpha*conjg( A' )*A + beta*C, matrix multiply
*                 updating upper vertical blocks of C.
*
                  IF( II.GT.1 )THEN
                     CALL ZGEMM ( 'C', 'N', II-1, ISEC, K, ZALPHA,
     $                                 A( 1, 1 ), LDA, A( 1, II ), LDA,
     $                                         ZBETA, C( 1, II ), LDC )
                  END IF
                  IF( TINYK )THEN
*
*                    C :=  beta*C, a upper triangular diagonal block
*                    of C is updated with beta. The imaginary part of
*                    the diagonal elements of C are set to ZERO.
*
                     IF( BETA.NE.ONE )THEN
                        C( II, II ) =
     $                         DCMPLX( BETA*DBLE( C( II, II ) ), ZERO )
                        DO 150, I = II+1, II+ISEC-1
                           CALL ZSCAL ( I-II, ZBETA, C( II, I ), 1 )
                           C( I, I ) =
     $                           DCMPLX( BETA*DBLE( C( I, I ) ), ZERO )
  150                   CONTINUE
                     END IF
*
*                    C := alpha*conjg( A' )*A + C, hermitian matrix
*                    multiply. C is a hermitian diagonal block having
*                    upper triangular storage format. The local array
*                    T3 is used for temporary storage of the conjugate
*                    transposed vectors of A.
*
                     DO 170, L = 1, K
                        DO 160, I = II, II+ISEC-1
                           T3( I-II+1, 1 ) = DCONJG( A( L, I ) )
  160                   CONTINUE
                        CALL ZHER  ( 'U', ISEC, ALPHA, T3( 1, 1 ),
     $                                            1, C( II, II ), LDC )
  170                CONTINUE
                  ELSE
*
*                    T2 := C, a upper triangular diagonal block of the
*                    hermitian matrix C is copied to the upper
*                    triangular part of T2.
*
                     DO 180, I = II, II+ISEC-1
                        CALL ZCOPY ( I-II+1, C( II, I ), 1,
     $                                             T2( 1, I-II+1 ), 1 )
  180                CONTINUE
*
*                    T2 :=  beta*T2, the upper triangular part of T2 is
*                    updated with beta.
*
                     IF( BETA.NE.ONE )THEN
                        DO 190, I = II, II+ISEC-1
                           CALL ZSCAL ( I-II+1, ZBETA,
     $                                             T2( 1, I-II+1 ), 1 )
  190                   CONTINUE
                     END IF
                     DO 230, LL = 1, K, RCB
                        LSEC = MIN( RCB, K-LL+1 )
*
*                       T3 :=  A', the transpose of a square block of A
*                       is copied to T3.
*
                        DO 200, I = II, II+ISEC-1
                           CALL ZCOPY ( LSEC, A( LL, I ), 1,
     $                                           T3( I-II+1, 1 ), RCB )
  200                   CONTINUE
*
*                       T2 := alpha*conjg( T3' )*T3 + T2, hermitian
*                       matrix multiply. T2 contains a hermitian block
*                       having upper triangular storage format. The
*                       local array T3 is used for temporary storage of
*                       the conjugate transposed vectors of A.
*
                        DO 220, L = LL, LL+LSEC-1
                           DO 210, I = 1, ISEC
                              T3( I, L-LL+1 ) =
     $                                        DCONJG( T3( I, L-LL+1 ) )
  210                      CONTINUE
                           CALL ZHER  ( 'U', ISEC, ALPHA,
     $                            T3( 1, L-LL+1 ), 1, T2( 1, 1 ), RCB )
  220                   CONTINUE
  230                CONTINUE
*
*                    C := T2, the upper triangular part of T2 is copied
*                    back to C.
*
                     DO 240, I = II, II+ISEC-1
                        CALL ZCOPY ( I-II+1, T2( 1, I-II+1 ), 1,
     $                                                  C( II, I ), 1 )
  240                CONTINUE
                  END IF
  250          CONTINUE
            ELSE
               CLDA = ZGB91( LDA )
               DO 330, II = 1, N, RB
                  ISEC = MIN( RB, N-II+1 )
*
*                 C := alpha*conjg( A' )*A + beta*C, matrix multiply
*                 updating upper vertical blocks of C.
*
                  IF( II.GT.1 )THEN
                     CALL ZGEMM ( 'C', 'N', II-1, ISEC, K, ZALPHA,
     $                                 A( 1, 1 ), LDA, A( 1, II ), LDA,
     $                                         ZBETA, C( 1, II ), LDC )
                  END IF
                  ZDELTA = ZBETA
                  DO 320, LL = 1, K, CB
                     LSEC = MIN( CB, K-LL+1 )
*
*                    T1 := conjg( A' ), the conjugated transpose of a
*                    rectangular block of A is copied to T1.
*
                     IF( CLDA )THEN
                        DO 270, I = II, II+ISEC-1
                           DO 260, L = LL, LL+LSEC-1
                              T1( I-II+1, L-LL+1 ) =
     $                                              DCONJG( A( L, I ) )
  260                      CONTINUE
  270                   CONTINUE
                     ELSE
                        DO 290, L = LL, LL+LSEC-1
                           DO 280, I = II, II+ISEC-1
                              T1( I-II+1, L-LL+1 ) =
     $                                              DCONJG( A( L, I ) )
  280                      CONTINUE
  290                   CONTINUE
                     END IF
*
*                    C := alpha*T1*conjg( T1' ) + delta*C, C is
*                    hermitian having triangular storage format. Delta
*                    is used instead of beta to avoid updating the
*                    block of C with beta multiple times. The local
*                    array T4 is used for the conjugated transpose
*                    of vectors of T1.
*
                     DO 310, I = II, II+ISEC-1
                        DO 300, L = LL, LL+LSEC-1
                           T4( L-LL+1 ) =
     $                                   DCONJG( T1( I-II+1, L-LL+1 ) )
  300                   CONTINUE
                        CALL ZGEMV ( 'N', I-II+1, LSEC, ZALPHA,
     $                                      T1( 1, 1 ), RB, T4( 1 ), 1,
     $                                          ZDELTA, C( II, I ), 1 )
                        C( I, I ) = DCMPLX( DBLE( C( I, I ) ), ZERO )
  310                CONTINUE
                     ZDELTA = ZONE
  320             CONTINUE
  330          CONTINUE
            END IF
         END IF
      ELSE
         IF( NOTR )THEN
*
*           Form  C := alpha*A*conjg( A' ) + beta*C. Lower, Notr.
*
            SMALLN = .NOT.ZGB90( ZIP51 , N, K )
            IF( SMALLN )THEN
               TINYK = .NOT.ZGB90( ZIP52 , N, K )
               DO 400, IX = N, 1, -RCB
                  II = MAX( 1, IX-RCB+1 )
                  ISEC = IX-II+1
                  IF( TINYK )THEN
*
*                    C :=  beta*C, a lower triangular diagonal block
*                    of C is updated with beta.
*
                     IF( BETA.NE.ONE )THEN
                        DO 340, I = II, II+ISEC-2
                           C( I, I ) =
     $                           DCMPLX( BETA*DBLE( C( I, I ) ), ZERO )
                           CALL ZSCAL ( II+ISEC-I-1, ZBETA,
     $                                                 C( I+1, I ), 1 )
  340                   CONTINUE
                        C( II+ISEC-1, II+ISEC-1 ) =
     $                                 DCMPLX( BETA*DBLE( C( II+ISEC-1,
     $                                            II+ISEC-1 ) ), ZERO )
                     END IF
*
*                    C := alpha*A*conjg( A' ) + C, hermitian matrix
*                    multiply. C is a hermitian diagonal block having
*                    lower triangular storage format.
*
                     DO 350, L = 1, K
                        CALL ZHER  ( 'L', ISEC, ALPHA, A( II, L ),
     $                                            1, C( II, II ), LDC )
  350                CONTINUE
                  ELSE
*
*                    T2 := C, a lower triangular diagonal block of the
*                    hermitian matrix C is copied to the lower
*                    triangular part of T2.
*
                     DO 360, I = II, II+ISEC-1
                        CALL ZCOPY ( II+ISEC-I, C( I, I ), 1,
     $                                        T2( I-II+1, I-II+1 ), 1 )
  360                CONTINUE
*
*                    T2 :=  beta*T2, the lower triangular part of T2 is
*                    updated with beta. The imaginary part of the
*                    diagonal elements of T2 are set to ZERO.
*
                     IF( BETA.NE.ONE )THEN
                        DO 370, I = 1, ISEC-1
                           T2( I, I ) =
     $                          DCMPLX( BETA*DBLE( T2( I, I ) ), ZERO )
                           CALL ZSCAL ( ISEC-I, ZBETA,
     $                                                T2( I+1, I ), 1 )
  370                   CONTINUE
                        T2( ISEC, ISEC ) =
     $                           DCMPLX( BETA*DBLE( T2( ISEC, ISEC ) ),
     $                                                           ZERO )
                     END IF
*
*                    T2 := alpha*A*conjg( A' ) + T2, symmetric matrix
*                    multiply. T2 contains a hermitian block having
*                    lower triangular storage format.
*
                     DO 380, L = 1, K
                        CALL ZHER  ( 'L', ISEC, ALPHA, A( II, L ),
     $                                             1, T2( 1, 1 ), RCB )
  380                CONTINUE
*
*                    C := T2, the lower triangular part of T2 is copied
*                    back to C.
*
                     DO 390, I = II, II+ISEC-1
                        CALL ZCOPY ( II+ISEC-I, T2( I-II+1, I-II+1 ),
     $                                                1, C( I, I ), 1 )
  390                CONTINUE
                  END IF
*
*                 C := alpha*A*conjg( A' ) + beta*C, matrix multiply
*                 on lower vertical blocks of C.
*
                  IF( II+ISEC.LE.N )THEN
                     CALL ZGEMM ( 'N', 'C', N-II-ISEC+1, ISEC, K,
     $                        ZALPHA, A( II+ISEC, 1 ), LDA, A( II, 1 ),
     $                              LDA, ZBETA, C( II+ISEC, II ), LDC )
                  END IF
  400          CONTINUE
            ELSE
               DO 450, IX = N, 1, -RB
                  II = MAX( 1, IX-RB+1 )
                  ISEC = IX-II+1
                  ZDELTA = ZBETA
                  DO 440, LL = 1, K, CB
                     LSEC = MIN( CB, K-LL+1 )
*
*                    T1 := A, a rectangular block of A is copied to T1.
*
                     DO 410, L = LL, LL+LSEC-1
                        CALL ZCOPY ( ISEC, A( II, L ), 1,
     $                                             T1( 1, L-LL+1 ), 1 )
  410                CONTINUE
*
*                    C := alpha*T1*conjg( T1' ) + delta*C, C is
*                    hermitian having triangular storage format. Delta
*                    is used instead of beta to avoid updating the
*                    block of C with beta multiple times. The local
*                    array T4 is used for the conjugated transpose
*                    of vectors of T1.
*
                     DO 430, I = II, II+ISEC-1
                        DO 420, L = LL, LL+LSEC-1
                           T4( L-LL+1 ) =
     $                                   DCONJG( T1( I-II+1, L-LL+1 ) )
  420                   CONTINUE
                        CALL ZGEMV ( 'N', II+ISEC-I, LSEC, ZALPHA,
     $                                 T1( I-II+1, 1 ), RB, T4( 1 ), 1,
     $                                           ZDELTA, C( I, I ), 1 )
                        C( I, I ) = DCMPLX( DBLE( C( I, I ) ), ZERO )
  430                CONTINUE
                     ZDELTA = ZONE
  440             CONTINUE
*
*                 C := alpha*A*conjg( A' ) + beta*C, matrix multiply
*                 updating lower vertical blocks of C.
*
                  IF( II+ISEC.LE.N )THEN
                     CALL ZGEMM ( 'N', 'C', N-II-ISEC+1, ISEC, K,
     $                        ZALPHA, A( II+ISEC, 1 ), LDA, A( II, 1 ),
     $                              LDA, ZBETA, C( II+ISEC, II ), LDC )
                  END IF
  450          CONTINUE
            END IF
         ELSE
*
*           Form  C := alpha*conjg( A' )*A + beta*C. Lower, Trans.
*
            SMALLN = .NOT.ZGB90( ZIP51 , N, K )
            IF( SMALLN )THEN
               TINYK = .NOT.ZGB90( ZIP52 , N, K )
               DO 560, IX = N, 1, -RCB
                  II = MAX( 1, IX-RCB+1 )
                  ISEC = IX-II+1
                  IF( TINYK )THEN
*
*                    C :=  beta*C, a lower triangular diagonal block
*                    of C is updated with beta. The imaginary part of
*                    the diagonal elements of C are set to ZERO.
*
                     IF( BETA.NE.ONE )THEN
                        DO 460, I = II, II+ISEC-2
                           C( I, I ) =
     $                           DCMPLX( BETA*DBLE( C( I, I ) ), ZERO )
                           CALL ZSCAL ( II+ISEC-I-1, ZBETA,
     $                                                 C( I+1, I ), 1 )
  460                   CONTINUE
                        C( II+ISEC-1, II+ISEC-1 ) =
     $                                 DCMPLX( BETA*DBLE( C( II+ISEC-1,
     $                                            II+ISEC-1 ) ), ZERO )
                     END IF
*
*                    C := alpha*conjg( A' )*A + C, hermitian matrix
*                    multiply. C is a hermitian diagonal block having
*                    lower triangular storage format. The local array
*                    T3 is used for temporary storage of the conjugate
*                    transposed vectors of A.
*
                     DO 480, L = 1, K
                        DO 470, I = II, II+ISEC-1
                           T3( I-II+1, 1 ) = DCONJG( A( L, I ) )
  470                   CONTINUE
                        CALL ZHER  ( 'L', ISEC, ALPHA, T3( 1, 1 ),
     $                                            1, C( II, II ), LDC )
  480                CONTINUE
                  ELSE
*
*                    T2 := C, a lower triangular diagonal block of the
*                    symmetric matrix C is copied to the lower
*                    triangular part of T2.
*
                     DO 490, I = II, II+ISEC-1
                        CALL ZCOPY ( II+ISEC-I, C( I, I ), 1,
     $                                        T2( I-II+1, I-II+1 ), 1 )
  490                CONTINUE
*
*                    T2 :=  beta*T2, the lower triangular part of T2 is
*                    updated with beta.
*
                     IF( BETA.NE.ONE )THEN
                        DO 500, I = II, II+ISEC-1
                           CALL ZSCAL ( II+ISEC-I, ZBETA,
     $                                        T2( I-II+1, I-II+1 ), 1 )
  500                   CONTINUE
                     END IF
                     DO 540, LL = 1, K, RCB
                        LSEC = MIN( RCB, K-LL+1 )
*
*                       T3 :=  A', the transpose of a square block of A
*                       is copied to T3.
*
                        DO 510, I = II, II+ISEC-1
                           CALL ZCOPY ( LSEC, A( LL, I ), 1,
     $                                           T3( I-II+1, 1 ), RCB )
  510                   CONTINUE
*
*                       T2 := alpha*conjg( T3' )*T3 + T2, hermitian
*                       matrix multiply. T2 contains a hermitian block
*                       having lower triangular storage format. The
*                       local array T3 is used for temporary storage of
*                       the conjugate transposed vectors of A.
*
                        DO 530, L = LL, LL+LSEC-1
                           DO 520, I = 1, ISEC
                              T3( I, L-LL+1 ) =
     $                                        DCONJG( T3( I, L-LL+1 ) )
  520                      CONTINUE
                           CALL ZHER  ( 'L', ISEC, ALPHA,
     $                            T3( 1, L-LL+1 ), 1, T2( 1, 1 ), RCB )
  530                   CONTINUE
  540                CONTINUE
*
*                    C := T2, the lower triangular part of T2 is copied
*                    back to C.
*
                     DO 550, I = II, II+ISEC-1
                        CALL ZCOPY ( II+ISEC-I, T2( I-II+1, I-II+1 ),
     $                                                1, C( I, I ), 1 )
  550                CONTINUE
                  END IF
*
*                 C := alpha*conjg( A' )*A + beta*C, matrix multiply
*                 updating lower vertical blocks of C.
*
                  IF( II+ISEC.LE.N )THEN
                     CALL ZGEMM ( 'C', 'N', N-II-ISEC+1, ISEC, K,
     $                        ZALPHA, A( 1, II+ISEC ), LDA, A( 1, II ),
     $                              LDA, ZBETA, C( II+ISEC, II ), LDC )
                  END IF
  560          CONTINUE
            ELSE
               CLDA = ZGB91( LDA )
               DO 650, IX = N, 1, -RB
                  II = MAX( 1, IX-RB+1 )
                  ISEC = IX-II+1
                  ZDELTA = ZBETA
                  DO 640, LL = 1, K, CB
                     LSEC = MIN( CB, K-LL+1 )
*
*                    T1 := conjg( A' ), the conjugated transpose of a
*                    rectangular block of A is copied to T1.
*
                     IF( CLDA )THEN
                        DO 580, I = II, II+ISEC-1
                           DO 570, L = LL, LL+LSEC-1
                              T1( I-II+1, L-LL+1 ) =
     $                                              DCONJG( A( L, I ) )
  570                      CONTINUE
  580                   CONTINUE
                     ELSE
                        DO 600, L = LL, LL+LSEC-1
                           DO 590, I = II, II+ISEC-1
                              T1( I-II+1, L-LL+1 ) =
     $                                              DCONJG( A( L, I ) )
  590                      CONTINUE
  600                   CONTINUE
                     END IF
*
*                    C := alpha*T1*conjg( T1' ) + delta*C, C is
*                    hermitian having triangular storage format. Delta
*                    is used instead of beta to avoid updating the
*                    block of C with beta multiple times. The local
*                    array T4 is used for the conjugated transpose
*                    of vectors of T1.
*
                     DO 630, I = II, II+ISEC-1
                        DO 620, L = LL, LL+LSEC-1
                           T4( L-LL+1 ) =
     $                                   DCONJG( T1( I-II+1, L-LL+1 ) )
  620                   CONTINUE
                        CALL ZGEMV ( 'N', II+ISEC-I, LSEC, ZALPHA,
     $                                 T1( I-II+1, 1 ), RB, T4( 1 ), 1,
     $                                           ZDELTA, C( I, I ), 1 )
                        C( I, I ) = DCMPLX( DBLE( C( I, I ) ), ZERO )
  630                CONTINUE
                     ZDELTA = ZONE
  640             CONTINUE
*
*                 C := alpha*conjg( A' )*A + beta*C, matrix multiply
*                 updating lower vertical blocks of C.
*
                  IF( II+ISEC.LE.N )THEN
                     CALL ZGEMM ( 'C', 'N', N-II-ISEC+1, ISEC, K,
     $                        ZALPHA, A( 1, II+ISEC ), LDA, A( 1, II ),
     $                              LDA, ZBETA, C( II+ISEC, II ), LDC )
                  END IF
  650          CONTINUE
            END IF
         END IF
      END IF
*
      RETURN
*
*     End of ZGB05.
*
      END
