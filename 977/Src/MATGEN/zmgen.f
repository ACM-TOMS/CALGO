*{SIGMA/MATGEN/zmgen.f}
*
*  Definition:
*  ===========
*      SUBROUTINE ZMGEN( M, N, A, LDA, SCOND, COND, SMODE,
*     $                    CMODE, DISTR, ISEED, D, X, INFO )
*
*
* SIGMA library, MATGEN section updated February 2016. 
* Developed and coded by Zlatko Drmac, Department of Mathematics
* University of Zagreb, Croatia, drmac@math.hr
* Submitted to ACM TOMS.      
*
*     .. Scalar Arguments
*      DOUBLE PRECISION   COND, SCOND
*      INTEGER            CMODE, DISTR, INFO, LDA, M, N, SMODE
*     ..
*     .. Array Arguments
*      DOUBLE COMPLEX    A(LDA,*),  X(*)
*      DOUBLE PRECISION  D(*) 
*      INTEGER           ISEED(4)
*     ..
*
*  Purpose
*  =======
*
*  ZMGEN  generates a random non-singular matrix  A, with
*  the following structure:
*  A = B * D, where B has unit columns (in Euclidean norm)
*  and the spectral condition number of B is SCOND; D is
*  diagonal with condition number COND.
*
*  Arguments
*  =========
*  M       (input) INTEGER
*          Row dimension of the matrix A. M > 0.
*...............................................................................
*  N       (input) INTEGER
*          Column dimension of the matrix A, M >= N > 0.
*...............................................................................
*  A       (output) DOUBLE COMPLEX array, dimension (LDA,N).
*          The generated matrix. A has the following structure:
*          A = B * D, where B has unit columns (in Euclidean norm)
*          and the spectral condition number of B is SCOND; D is
*          diagonal with condition number COND.
*...............................................................................
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= M.
*...............................................................................
*  SCOND   (input) DOUBLE PRECISION
*          See the description of A.
*...............................................................................
*  COND    (input) DOUBLE PRECISION
*          See the description of A.
*...............................................................................
*  SMODE   (input) INTEGER
*          Distribution of the singular values of the matrix B. (See the
*          description of A and HOME/lapack/testing/matgen/*.*)
*...............................................................................
*  CMODE   (input) INTEGER
*          Distribution of the diagonal of the diagonal column scaling.
*          (See the descriptions of A and SMODE.)
*...............................................................................
*  DISTR   (input) INTEGER from {1,2,3}
*          Distribution for random number generator.
*...............................................................................
*  ISEED   (input/output) INTEGER array, dimension (4)
*          On entry ISEED specifies the seed of the random number
*          generator. The array elements should be between 0 and 4095;
*          if not they will be reduced mod 4096.  Also, ISEED(4) must
*          be odd.  The random number generator uses a linear
*          congruential sequence limited to small integers, and so
*          should produce machine independent random numbers. The
*          values of ISEED are changed on exit, and can be used in the
*          next call to continue the same random number sequence.
*...............................................................................
*  D       (workspace) DOUBLE PRECISION array, dimension (N)
*...............................................................................
*  X       (workspace) DOBLE COMPLEX dimension 3*MAX(M,N).
*...............................................................................
*  INFO    (output) INTEGER
*          = 0: successful exit;
*          < 0: if INFO = -i, the i-th argument had an illegal value.
*...............................................................................
*
      SUBROUTINE ZMGEN( M, N, A, LDA, SCOND, COND, SMODE,
     &                    CMODE, DISTR, ISEED, D, X, INFO )
*
*
*     .. Scalar Arguments
      DOUBLE PRECISION   COND, SCOND
      INTEGER            CMODE, DISTR, INFO, LDA, M, N, SMODE
*     ..
*     .. Array Arguments
      COMPLEX*16        A( LDA, * ),  X( * )
      DOUBLE PRECISION  D(*) 
      INTEGER           ISEED(4)
*     .. Parameters ..
      DOUBLE PRECISION    ZERO, ONE        , TWO        , HALF
      PARAMETER ( ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0, HALF =0.5D0 )
*
*     ..
*     .. Local Scalars ..
      INTEGER            i, it, j, k
      DOUBLE PRECISION   Aik, CS, SN, T, TEMP, TOL
      COMPLEX*16         OM
*     ..
*     .. External Subroutines ..
      EXTERNAL           ZLAROR, ZLASET, DLATM1, ZROT, DSCAL, XERBLA
*     ..
*     .. External Functions ..
      DOUBLE PRECISION           DLAMCH, DNRM2
      COMPLEX*16         ZDOTC 
      EXTERNAL           ZDOTC,  DLAMCH, DNRM2
      INTEGER            IDAMAX, IDAMIN
      EXTERNAL           IDAMAX, IDAMIN 
*
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, CONJG, DBLE, DCMPLX, SIGN, SQRT
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
          INFO = 0
      IF( M .LE. 0 ) THEN
          INFO = -1
      ELSE IF ( (N .LE. 0).OR.(N .GT. M) ) THEN
          INFO = -2
      ELSE IF( LDA .LT. M ) THEN
          INFO = -4
      ELSE IF ( SCOND .LT. ONE ) THEN
          INFO = -5
      ELSE IF ( (N.EQ.1) .AND. (SCOND.NE.ONE) ) THEN
          INFO = -5
      ELSE IF ( COND .LT. ONE ) THEN
          INFO = -6
      ELSE IF ( IABS(SMODE) .GT. 7 ) THEN
          INFO = -7
      ELSE IF ( IABS(CMODE) .GT. 7 ) THEN
          INFO = -8
      ELSE IF ( (DISTR .LT. 1) .OR. (DISTR .GT. 3) ) THEN
          INFO = -9
      END IF
      IF( INFO .NE. 0 ) THEN
         CALL XERBLA( 'ZMGEN', -INFO )
         RETURN
      END IF
*
*     Set  A := diag( DLATM1() ) = diag(random singular values).
*     The vector of singular values is normalized to N in 2-norm.
*
      CALL DLATM1( SMODE, SCOND, 0, DISTR, ISEED, D, N, INFO )
      DO 10 i = 1, N
         D(i) = ABS(D(i))
  10   CONTINUE
      i = IDAMIN( N, D, 1 )
      k = IDAMAX( N, D, 1 )
      SN = D(i)
      CS = D(k)
      IF ((CS-SN.NE.ZERO) .AND. ((D(k) .GT. TWO*SCOND*D(i)) .OR.
     &     (D(k) .LT. HALF*SCOND*D(i))) ) THEN
*     .. max ABS(D(:)) / min ABS(D(:)) not close to the desired values SCOND ;
*     transform D(:) to get SCOND. Here we use linear map from the interval
*     [min_i D(i), max_i D(i)] to [1/SCOND,1], for the sake of simplicity.  
*     One can use other choices, e.g. quadratic function.           
      DO 101 j = 1, N
          D(j) = ONE / SCOND +
     &          ( ( ONE-( ONE / SCOND ) ) / ( CS-SN ) ) * ( D(j)-SN )
 101   CONTINUE
      END IF
*
      CALL DSCAL(N, SQRT(DBLE(N)) / DNRM2(N, D, 1), D, 1 )
      CALL ZLASET( 'A', M, N, DCMPLX(ZERO), DCMPLX(ZERO), A, LDA )
      DO 1 j = 1, N
         A(j,j) = DCMPLX(D(j))
   1  CONTINUE
*
*     Set A = U * A * W. Matrices U and W are orthogonal, generated using
*     the method of G.W. Stewart (SIAM J. Numer. Anal. 17, 1980, 403-409).
*
      CALL ZLAROR( 'L', 'N', M, N, A, LDA, ISEED, X, INFO )
      CALL ZLAROR( 'R', 'N', M, N, A, LDA, ISEED, X, INFO )
*
      IF ( N .EQ. 1 ) RETURN 
*     
*     At this point Trace(A^T * A) = N 
*      
*
*     Equilibrate the columns of A: use rotations to make them unit in
*     Euclidean norm; this is implicit complex version of the transformation
*     described in G. Golub, CH. Van Loan: Matrix Compuations, Second edition,
*     1989, problems [P 8.5.3, P 8.5.4]. Theoretically, equilibrated columns
*     are obtained in N-1 steps. Numerically, it does a good job. One could
*     repeat the process, but it does not seem necessary.
*
      TOL = DBLE(M+N)*DLAMCH('Epsilon')
      DO 12 it = 1, 3 
      DO 98 i = 1, N
         D(i) = ABS( ZDOTC(M, A(1,i), 1, A(1,i), 1) )
 98   CONTINUE
      DO 99 i = 1, N - 1
         j = IDAMIN(N-i+1,D(i),1) + i - 1
         k = IDAMAX(N-i+1,D(i),1) + i - 1
         IF ( i .NE. j ) THEN
            CALL ZSWAP( M, A(1,i), 1, A(1,j), 1 )
            TEMP = D(i)
            D(i) = D(j)
            D(j) = TEMP
            IF ( k .EQ. i ) THEN
               k=j
            ELSE IF ( k .EQ. j ) THEN
               k=i
            END IF
         END IF
         IF ((D(i).GE.(ONE-TOL)).OR.(D(k).LE.(ONE+TOL))) THEN
            GO TO 990
         END IF
         Aik = ZDOTC( M, A(1,i), 1, A(1,k), 1 )
         OM  = Aik / ABS(Aik) 
         Aik = ABS(Aik) 
         T   = SIGN( ABS(Aik)
     &         + SQRT( Aik*Aik + ABS((D(i)-ONE)*(D(k)-ONE))), Aik ) /
     &           (D(k) - ONE)
         CS  = ONE / SQRT(ONE+T*T)
         SN  = T * CS
         CALL ZROT( M, A(1,i), 1, A(1,k), 1, CS, CONJG(OM)*SN )
         D(k) = ABS( ZDOTC( M, A(1,k), 1, A(1,k), 1 ) )
 99   CONTINUE
 12   CONTINUE     
*
*     Multiply the columns of A by random diagonal matrix
*     with spectral condition number COND.
*
 990  CONTINUE

      IF ( ( N .GT. 1 ) .AND. ( COND .GT. ONE ) ) THEN

         CALL DLATM1( CMODE, COND, 0, DISTR, ISEED, D, N, INFO )
         DO 20 i = 1, N
            D(i) = DABS(D(i))
 20      CONTINUE
         i = IDAMIN( N, D, 1 )
         k = IDAMAX( N, D, 1 )
         SN = D(i)
         CS = D(k)
        IF ( ( D(k).NE.D(i) ) .AND. ((D(k) .GT. TWO*COND*D(i)) .OR.
     &        (D(k) .LT. HALF*COND*D(i))) ) THEN
            DO 201 j = 1, N
               D(j) = ONE / COND +
     &          ( ( ONE-( ONE / COND ) ) / ( CS-SN ) ) * ( D(j)-SN )
 201        CONTINUE
         END IF
*
         DO 90 j = 1, N
            CALL ZDSCAL( M, ONE/D(j), A(1,j), 1 )
   90    CONTINUE
       END IF
*
      RETURN
*
*     End of  ZMGEN
*
      END
