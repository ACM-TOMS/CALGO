      SUBROUTINE CLAQPP( M, N, OFFSET, NB, KB, A, LDA, JPVT, RCOND,
     $                   RANK, TAU, SMIN, XMIN, SMAX, XMAX, VN1, VN2,
     $                   AUXV, F, LDF )
*
*     This code is part of a package for solving rank deficient least
*     squares problems, written by:
*     ==================================================================
*     L. Foster                   and   R. Kommu
*     Department of Mathematics         Department of Physics
*     San Jose State University         San Jose State University
*     San Jose, CA 95192                San Jose, CA 95192
*     foster@math.sjsu.edu              rkommu@email.sjsu.edu
*     ==================================================================
*     03/05/2004
*
*     .. Scalar Arguments ..
      INTEGER            KB, LDA, LDF, M, N, NB, OFFSET, RANK
      REAL               RCOND, SMAX, SMIN
*     ..
*     .. Array Arguments ..
      INTEGER            JPVT( * )
      REAL               VN1( * ), VN2( * )
      COMPLEX            A( LDA, * ), AUXV( * ), F( LDF, * ), TAU( * ),
     $                   XMAX( * ), XMIN( * )
*     ..
*
*  Purpose
*  =======
*
*
*  CLAQPP uses a partial QR factorization with column pivoting applied
*  to A(OFFSET+1:M,1:N) to compute P, Q, R01, R02, R11, R12 and, in
*  some cases R22, in the factorization:
*      A * P = [I 0] * [ R01 R02 ]
*              [0 Q]   [ R11 R12 ]
*                      [  0  R22 ].
*  Here R01 and R02 have OFFSET rows, I is an OFFSET by OFFSET
*  identity matrix, Q is an (M - OFFSET) * (M - OFFSET) unitary
*  matrix represented by a product of elementary reflectors and R11
*  is upper triangular.  R11 is defined as the largest order matrix
*  with NB or fewer rows and columns so that
*  o  catastrophic cancellation does not keep CLAQPP from factoring
*     additional columns of A and
*  o  the estimated condition number of
*          R1 = [ R00 R01 ]
*               [  0  R11 ]
*     is less than 1 / RCOND.
*  The OFFSET by OFFSET matrix R00 is not passed into CLAQPP;
*  however estimates of largest and smallest singular values
*  and corresponding approximate singular vectors of R00 are
*  supplied to CLAQPP and these are used to estimate the condition
*  number of R1. On output RANK is the order of R1.
*
*  CLAQPP uses Level 3 BLAS.
*
*
*  Arguments
*  =========
*
*  M       (input) INTEGER
*          The number of rows of the matrix A. M >= 0.
*
*  N       (input) INTEGER
*          The number of columns of the matrix A. N >= 0
*
*  OFFSET  (input) INTEGER
*          The number of rows of A that have been factorized in
*          previous steps.
*
*  NB      (input) INTEGER
*          Factor at most NB columns of A.
*
*  KB      (output) INTEGER
*          KB is NB if catastrophic cancellation was not an issue
*          in factoring columns of A, otherwise KB is the number
*          of columns that can be safely factored.
*
*  A       (input/output) COMPLEX array, dimension (LDA,N)
*          On entry, the M-by-N matrix A.
*          On exit,
*          o A(1:OFFSET,:) contains [R01 R02],
*          o the upper triangular part of
*            A(OFFSET+1:RANK,1:RANK-OFFSET) contains the upper
*            triangular matrix R11,
*          o A(OFFSET+1:RANK,RANK-OFFSET+1:N) contains R12,
*          o the elements of A below the diagonal in
*            A(OFFSET+1:M,1:RANK-OFFSET), together with TAU,
*            represent Q as a product of RANK-OFFSET elementary
*            reflections, and
*          o A(RANK+1:M,RANK-OFFSET+1:N) contains R22 when RANK =
*            OFFSET + KB.  If RANK < OFFSET + KB
*            A(RANK+1:M,RANK-OFFSET+1:N) contains information related
*            to R22 and used in the calculation of P, Q, R01, R02,
*            R11 and R12.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A. LDA >= max(1,M).
*
*  JPVT    (input/output) INTEGER array, dimension (N)
*          JPVT(I) = K <==> Column K of the full matrix A has been
*          permuted into position I in AP.
*
*  RCOND   (input) REAL
*          Stop the factorization of A so that
*                   R1 = [ R00 R01 ]
*                        [  0  R11 ]
*          is the largest dimension matrix of this form with an
*          estimated condition number less than 1/RCOND.
*
*  RANK    (input/output) INTEGER
*          On entry, RANK = OFFSET.
*          On exit, the order of the matrix R1. Note that RANK <
*          OFFSET + KB if the factorization was stopped due to
*          the RCOND test.
*
*
*  TAU     (output) COMPLEX array, dimension (KB)
*          The scalar factors of the elementary reflectors.
*
*  SMIN    (input/output) REAL
*          On entry, an estimate of the smallest singular value of R00.
*          On exit, an estimate of the smallest singular value of R1.
*
*  XMIN    (input/output) COMPLEX array, dimension
*          MIN( (OFFSET+N), M )
*          On entry, an approximate left singular vector of R00
*          corresponding to SMIN.
*          On exit, an approximate left singular vector of R1
*          corresponding to SMIN.
*
*  SMAX    (input/output) REAL
*          On entry, an estimate of the largest singular value of R00.
*          On exit, an estimate of the largest singular value of R1.
*
*  XMAX    (input/output) COMPLEX array, dimension
*          MIN( (OFFSET+N), M )
*          On entry, an approximate left singular vector of R00
*          corresponding to SMAX.
*          On exit, an approximate left singular vector of R1
*          corresponding to SMAX.
*
*  VN1     (input/output) REAL array, dimension (N)
*          The vector with the partial column norms.
*
*  VN2     (input/output) REAL array, dimension (N)
*          The vector with the exact column norms.
*
*  AUXV    (input/output) COMPLEX array, dimension (NB)
*          Auxiliar vector.
*
*  F       (input/output) COMPLEX array, dimension (LDF,NB)
*          Matrix F' = L*Y'*A.
*
*  LDF     (input) INTEGER
*          The leading dimension of the array F. LDF >= max(1,N).
*
*  Further Details
*  ===============
*
*  This is a modification of LAPACK routine xLAQPS written by
*    G. Quintana-Orti, Depto. de Informatica, Universidad Jaime I, Spain
*    X. Sun, Computer Science Dept., Duke University, USA
*  Modified by L. Foster, Department of Mathematics, San Jose State
*    University, San Jose, CA.
*
*  =====================================================================
*
*     .. Parameters ..
      INTEGER            IMAX, IMIN
      PARAMETER          ( IMAX = 1, IMIN = 2 )
      REAL               ZERO, ONE
      COMPLEX            CZERO, CONE
      PARAMETER          ( ZERO = 0.0E+0, ONE = 1.0E+0,
     $                   CZERO = ( 0.0E+0, 0.0E+0 ),
     $                   CONE = ( 1.0E+0, 0.0E+0 ) )
*     ..
*     .. Local Scalars ..
      INTEGER            ITEMP, J, K, LASTRK, LSTICC, PVT, RK
      REAL               SMAXPR, SMINPR, TEMP, TEMP2
      COMPLEX            AKK, C1, C2, S1, S2
*     ..
*     .. External Subroutines ..
      EXTERNAL           CGEMM, CGEMV, CLAIC1, CLARFG, CSWAP
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, REAL, CONJG, MAX, MIN, NINT, SQRT
*     ..
*     .. External Functions ..
      INTEGER            ISAMAX
      REAL               SCNRM2
      EXTERNAL           ISAMAX, SCNRM2
*     ..
*     .. Executable Statements ..
*
      LASTRK = MIN( M, N+OFFSET )
      LSTICC = 0
      K = 0
*
*     Beginning of while loop.
*
   10 CONTINUE
      IF( ( K.LT.NB ) .AND. ( LSTICC.EQ.0 ) ) THEN
         K = K + 1
         RK = OFFSET + K
*
*        Determine ith pivot column and swap if necessary
*
         PVT = ( K-1 ) + ISAMAX( N-K+1, VN1( K ), 1 )
         IF( PVT.NE.K ) THEN
            CALL CSWAP( M, A( 1, PVT ), 1, A( 1, K ), 1 )
            CALL CSWAP( K-1, F( PVT, 1 ), LDF, F( K, 1 ), LDF )
            ITEMP = JPVT( PVT )
            JPVT( PVT ) = JPVT( K )
            JPVT( K ) = ITEMP
            VN1( PVT ) = VN1( K )
            VN2( PVT ) = VN2( K )
         END IF
*
*        Apply previous Householder reflectors to column K:
*        A(RK:M,K) := A(RK:M,K) - A(RK:M,1:K-1)*F(K,1:K-1)'.
*
         IF( K.GT.1 ) THEN
*CC            CALL CGEMM( 'No transpose', 'Conjugate transpose',
*CC     $                  M-RK+1, 1, K-1, -CONE, A( RK, 1 ), LDA,
*CC     $                  F( K, 1 ), LDF, CONE, A( RK, K ), LDA )
            DO 20 J = 1, K - 1
               F( K, J ) = CONJG( F( K, J ) )
   20       CONTINUE
            CALL CGEMV( 'No transpose', M-RK+1, K-1, -CONE, A( RK, 1 ),
     $                  LDA, F( K, 1 ), LDF, CONE, A( RK, K ), 1 )
            DO 30 J = 1, K - 1
               F( K, J ) = CONJG( F( K, J ) )
   30       CONTINUE
         END IF
*
*        Generate elementary reflector H(k).
*
         IF( RK.LT.M ) THEN
            CALL CLARFG( M-RK+1, A( RK, K ), A( RK+1, K ), 1, TAU( K ) )
         ELSE
            CALL CLARFG( 1, A( RK, K ), A( RK, K ), 1, TAU( K ) )
         END IF
*
*        Use incremental condition estimation to determine if the
*        factorization should be stopped
*
         IF( RANK.EQ.0 ) THEN
            XMIN( 1 ) = CONE
            XMAX( 1 ) = CONE
            SMAX = ABS( A( 1, 1 ) )
            SMIN = SMAX
            IF( ABS( A( 1, 1 ) ).EQ.ZERO ) THEN
*              Can exit early, partial factorization done since
*              the calculated rank is zero
               KB = NB
*              Note RANK < OFFSET+KB flags termination due to RCOND test
*              or when the RANK of A is zero
               RETURN
            ELSE
               RANK = 1
            END IF
         ELSE
*
            CALL CLAIC1( IMIN, RANK, XMIN, SMIN, A( 1, K ),
     $                   A( RANK+1, K ), SMINPR, S1, C1 )
            CALL CLAIC1( IMAX, RANK, XMAX, SMAX, A( 1, K ),
     $                   A( RANK+1, K ), SMAXPR, S2, C2 )
*
            IF( SMAXPR*RCOND.LE.SMINPR ) THEN
               DO 40 J = 1, RANK
                  XMIN( J ) = S1*XMIN( J )
                  XMAX( J ) = S2*XMAX( J )
   40          CONTINUE
               XMIN( RANK+1 ) = C1
               XMAX( RANK+1 ) = C2
               SMIN = SMINPR
               SMAX = SMAXPR
               RANK = RANK + 1
            ELSE
*              Can exit early, due to the RCOND test the partial
*              factorization is done
               KB = NB
*              Note RANK < OFFSET+KB flags termination due to RCOND test
               RETURN
            END IF
*
         END IF
*
         AKK = A( RK, K )
         A( RK, K ) = CONE
*
*        Compute Kth column of F:
*
*        Compute  F(K+1:N,K) := tau(K)*A(RK:M,K+1:N)'*A(RK:M,K).
*
         IF( K.LT.N ) THEN
            CALL CGEMV( 'Conjugate transpose', M-RK+1, N-K, TAU( K ),
     $                  A( RK, K+1 ), LDA, A( RK, K ), 1, CZERO,
     $                  F( K+1, K ), 1 )
         END IF
*
*        Padding F(1:K,K) with zeros.
*
         DO 50 J = 1, K
            F( J, K ) = ZERO
   50    CONTINUE
*
*        Incremental updating of F:
*        F(1:N,K) := F(1:N,K) - tau(K)*F(1:N,1:K-1)*A(RK:M,1:K-1)'
*                    *A(RK:M,K).
*
         IF( K.GT.1 ) THEN
            CALL CGEMV( 'Conjugate transpose', M-RK+1, K-1, -TAU( K ),
     $                  A( RK, 1 ), LDA, A( RK, K ), 1, CZERO,
     $                  AUXV( 1 ), 1 )
*
            CALL CGEMV( 'No transpose', N, K-1, CONE, F( 1, 1 ), LDF,
     $                  AUXV( 1 ), 1, CONE, F( 1, K ), 1 )
         END IF
*
*        Update the current row of A:
*        A(RK,K+1:N) := A(RK,K+1:N) - A(RK,1:K)*F(K+1:N,1:K)'.
*
         IF( K.LT.N ) THEN
            CALL CGEMM( 'No transpose', 'Conjugate transpose', 1, N-K,
     $                  K, -CONE, A( RK, 1 ), LDA, F( K+1, 1 ), LDF,
     $                  CONE, A( RK, K+1 ), LDA )
         END IF
*
*        Update partial column norms.
*
         IF( RK.LT.LASTRK ) THEN
            DO 60 J = K + 1, N
               IF( VN1( J ).NE.ZERO ) THEN
                  TEMP = ABS( A( RK, J ) ) / VN1( J )
                  TEMP = MAX( ZERO, ( ONE+TEMP )*( ONE-TEMP ) )
                  TEMP2 = ONE + 0.05*TEMP*( VN1( J ) / VN2( J ) )**2
                  IF( TEMP2.EQ.ONE ) THEN
                     VN2( J ) = REAL( LSTICC )
                     LSTICC = J
                  ELSE
                     VN1( J ) = VN1( J )*SQRT( TEMP )
                  END IF
               END IF
   60       CONTINUE
         END IF
*
         A( RK, K ) = AKK
*
*        End of while loop.
*
         GO TO 10
      END IF
      KB = K
      RK = OFFSET + KB
*
*     Apply the block reflector to the rest of the matrix:
*     A(OFFSET+KB+1:M,KB+1:N) := A(OFFSET+KB+1:M,KB+1:N) -
*                         A(OFFSET+KB+1:M,1:KB)*F(KB+1:N,1:KB)'.
*
      IF( KB.LT.MIN( N, M-OFFSET ) ) THEN
         CALL CGEMM( 'No transpose', 'Conjugate transpose', M-RK, N-KB,
     $               KB, -CONE, A( RK+1, 1 ), LDA, F( KB+1, 1 ), LDF,
     $               CONE, A( RK+1, KB+1 ), LDA )
      END IF
*
*     Recomputation of difficult columns.
*
   70 CONTINUE
      IF( LSTICC.GT.0 ) THEN
         ITEMP = NINT( VN2( LSTICC ) )
         VN1( LSTICC ) = SCNRM2( M-RK, A( RK+1, LSTICC ), 1 )
         VN2( LSTICC ) = VN1( LSTICC )
         LSTICC = ITEMP
         GO TO 70
      END IF
*
      RETURN
*
*     End of CLAQPP
*
      END
