      SUBROUTINE CLAQPQ( M, N, OFFSET, A, LDA, JPVT, RCOND, RANK, TAU,
     $                   SMIN, XMIN, SMAX, XMAX, VN1, VN2, WORK )
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
      INTEGER            LDA, M, N, OFFSET, RANK
      REAL               RCOND, SMAX, SMIN
*     ..
*     .. Array Arguments ..
      INTEGER            JPVT( * )
      REAL               VN1( * ), VN2( * )
      COMPLEX            A( LDA, * ), TAU( * ), WORK( * ), XMAX( * ),
     $                   XMIN( * )
*     ..
*
*  Purpose
*  =======
*
*  CLAQPQ uses a partial QR factorization with column pivoting applied
*  to A(OFFSET+1:M,1:N) to compute P, Q, R01, R02, R11, and R12 in the
*  factorization:
*      A * P = [I 0] * [ R01 R02 ]
*              [0 Q]   [ R11 R12 ]
*                      [  0  R22 ].
*  Here R01 and R02 have OFFSET rows, I is an OFFSET by OFFSET
*  identity matrix, Q is an (M - OFFSET) * (M - OFFSET) unitary
*  matrix represented by a product of elementary reflectors and R11
*  is upper triangular.  R11 is defined as the largest matrix so
*  that the estimated condition number of
*       R1 = [ R00 R01 ]
*            [  0  R11 ]
*  is less than 1/RCOND.  The OFFSET by OFFSET matrix R00 is not
*  passed into CLAQPQ; however estimates of largest and smallest
*  singular values and corresponding approximate singular vectors
*  of R00 are supplied to CLAQPQ and these are used to estimate
*  the condition number of R1. On output RANK is the order of R1.
*
*  CLAQPQ uses Level 2 BLAS.
*
*  Arguments
*  =========
*
*  M       (input) INTEGER
*          The number of rows of the matrix A. M >= 0.
*
*  N       (input) INTEGER
*          The number of columns of the matrix A. N >= 0.
*
*  OFFSET  (input) INTEGER
*          The number of rows of the matrix A that must be pivoted
*          but not factored. OFFSET >= 0.
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
*          o A(RANK+1:M,RANK-OFFSET+1:N) contains information related
*            to R22 and used in the calculation of P, Q, R01,
*            R02, R11 and R12.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A. LDA >= max(1,M).
*
*  JPVT    (input/output) INTEGER array, dimension (N)
*          JPVT(I) = K <==> Column K of the full matrix A has been
*          permuted into position I in AP.
*
*  RCOND   (input) REAL
*          RCOND is used to determine how much of A to factor using
*          the QR factorization with pivoting. Stop the factorization
*          of A so that
*                   R1 = [ R00 R01 ]
*                        [  0  R11 ]
*          is the largest dimension matrix of this form with an
*          estimated condition number less than 1/RCOND.
*
*  RANK    (input/output) INTEGER
*          On entry, RANK = OFFSET.
*          On exit, the order of the matrix R1.
*
*  TAU     (output) COMPLEX array, dimension (min(M,N))
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
*  WORK    (workspace) COMPLEX array, dimension (N)
*
*  Further Details
*  ===============
*
*  This is a modification of LAPACK routine xLAQP2 written by
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
      COMPLEX            CONE
      PARAMETER          ( ZERO = 0.0E+0, ONE = 1.0E+0,
     $                   CONE = ( 1.0E+0, 0.0E+0 ) )
*     ..
*     .. Local Scalars ..
      INTEGER            I, ITEMP, J, MN, OFFPI, PVT
      REAL               SMAXPR, SMINPR, TEMP, TEMP2
      COMPLEX            AII, C1, C2, S1, S2
*     ..
*     .. External Subroutines ..
      EXTERNAL           CLAIC1, CLARF, CLARFG, CSWAP
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, CONJG, MAX, MIN, SQRT
*     ..
*     .. External Functions ..
      INTEGER            ISAMAX
      REAL               SCNRM2
      EXTERNAL           ISAMAX, SCNRM2
*     ..
*     .. Executable Statements ..
*
      MN = MIN( M-OFFSET, N )
      RANK = OFFSET
*
*     Compute factorization.
*
      DO 30 I = 1, MN
*
         OFFPI = OFFSET + I
*
*        Determine ith pivot column and swap if necessary.
*
         PVT = ( I-1 ) + ISAMAX( N-I+1, VN1( I ), 1 )
*
         IF( PVT.NE.I ) THEN
            CALL CSWAP( M, A( 1, PVT ), 1, A( 1, I ), 1 )
            ITEMP = JPVT( PVT )
            JPVT( PVT ) = JPVT( I )
            JPVT( I ) = ITEMP
            VN1( PVT ) = VN1( I )
            VN2( PVT ) = VN2( I )
         END IF
*
*        Generate elementary reflector H(i).
*
         IF( OFFPI.LT.M ) THEN
            CALL CLARFG( M-OFFPI+1, A( OFFPI, I ), A( OFFPI+1, I ), 1,
     $                   TAU( I ) )
         ELSE
            CALL CLARFG( 1, A( M, I ), A( M, I ), 1, TAU( I ) )
         END IF
*
*
*        Determine if the factorization should be stopped using
*        incremental condition estimation
*
         IF( RANK.EQ.0 ) THEN
            XMIN( 1 ) = CONE
            XMAX( 1 ) = CONE
            SMAX = ABS( A( 1, 1 ) )
            SMIN = SMAX
            IF( ABS( A( 1, 1 ) ).EQ.ZERO ) THEN
*              Can exit early, partial factorization done since
*              the rank is zero
               RETURN
            ELSE
               RANK = 1
            END IF
         ELSE
*
            CALL CLAIC1( IMIN, RANK, XMIN, SMIN, A( 1, I ),
     $                   A( RANK+1, I ), SMINPR, S1, C1 )
            CALL CLAIC1( IMAX, RANK, XMAX, SMAX, A( 1, I ),
     $                   A( RANK+1, I ), SMAXPR, S2, C2 )
*
            IF( SMAXPR*RCOND.LE.SMINPR ) THEN
               DO 10 J = 1, RANK
                  XMIN( J ) = S1*XMIN( J )
                  XMAX( J ) = S2*XMAX( J )
   10          CONTINUE
               XMIN( RANK+1 ) = C1
               XMAX( RANK+1 ) = C2
               SMIN = SMINPR
               SMAX = SMAXPR
               RANK = RANK + 1
            ELSE
*               Can exit early, due to the RCOND test the partial
*               factorization is done
               RETURN
            END IF
*
         END IF
*
*
         IF( I.LT.N ) THEN
*
*           Apply H(i)' to A(offset+i:m,i+1:n) from the left.
*
            AII = A( OFFPI, I )
            A( OFFPI, I ) = CONE
            CALL CLARF( 'Left', M-OFFPI+1, N-I, A( OFFPI, I ), 1,
     $                  CONJG( TAU( I ) ), A( OFFPI, I+1 ), LDA,
     $                  WORK( 1 ) )
            A( OFFPI, I ) = AII
         END IF
*
*        Update partial column norms.
*
         DO 20 J = I + 1, N
            IF( VN1( J ).NE.ZERO ) THEN
               TEMP = ONE - ( ABS( A( OFFPI, J ) ) / VN1( J ) )**2
               TEMP = MAX( TEMP, ZERO )
               TEMP2 = ONE + 0.05*TEMP*( VN1( J ) / VN2( J ) )**2
               IF( TEMP2.EQ.ONE ) THEN
                  IF( OFFPI.LT.M ) THEN
                     VN1( J ) = SCNRM2( M-OFFPI, A( OFFPI+1, J ), 1 )
                     VN2( J ) = VN1( J )
                  ELSE
                     VN1( J ) = ZERO
                     VN2( J ) = ZERO
                  END IF
               ELSE
                  VN1( J ) = VN1( J )*SQRT( TEMP )
               END IF
            END IF
   20    CONTINUE
*
   30 CONTINUE
*
      RETURN
*
*     End of CLAQPQ
*
      END
