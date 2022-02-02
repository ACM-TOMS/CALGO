      SUBROUTINE DGESQB( M, N, A, LDA, B, LDB, CS, TAU, DWORK, LDWORK,
     $                   INFO )
      IMPLICIT NONE
C
C     PURPOSE
C
C     Computes a symplectic QR decomposition of a real 2M-by-N matrix
C     [A; B],
C
C               [ A ]             [ R11  R12 ]
C               [   ] = Q * R = Q [          ],
C               [ B ]             [ R21  R22 ]
C
C     where Q is a symplectic orthogonal matrix, R11 is upper triangular
C     and R21 is strictly upper triangular.
C     If [A; B] is symplectic then, theoretically, R21 = 0 and
C     R22 = inv(R11)^T. Blocked version.
C
C     ARGUMENTS
C
C     Input/Output Parameters
C
C     M       (input) INTEGER
C             The number of rows of A. M >= 0.
C
C     N       (input) INTEGER
C             The number of columns of A. N >= 0.
C
C     A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
C             On entry, the leading M-by-N part of this array must
C             contain the matrix A.
C             On exit, the leading M-by-N part of this array contains
C             the matrix [ R11  R12 ] and, in the zero parts of R,
C             information about the elementary reflectors used to
C             compute the symplectic QR decomposition.
C
C     LDA     INTEGER
C             The leading dimension of the array A.  LDA >= MAX(1,M).
C
C     B       (input/output) DOUBLE PRECISION array, dimension (LDB,N)
C             On entry, the leading M-by-N part of this array must
C             contain the matrix B.
C             On exit, the leading M-by-N part of this array contains
C             the matrix [ R21  R22 ] and, in the zero parts of B,
C             information about the elementary reflectors used to
C             compute the symplectic QR decomposition.
C
C     LDB     INTEGER
C             The leading dimension of the array B.  LDB >= MAX(1,M).
C
C     CS      (output) DOUBLE PRECISION array, dimension (2 * min(M,N))
C             On exit, the first 2*min(M,N) elements of this array
C             contain the cosines and sines of the symplectic Givens
C             rotations used to compute the symplectic QR decomposition.
C
C     TAU     (output) DOUBLE PRECISION array, dimension (min(M,N))
C             On exit, the first min(M,N) elements of this array
C             contain the scalar factors of some of the elementary
C             reflectors.
C
C     Workspace
C
C     DWORK   DOUBLE PRECISION array, dimension (LDWORK)
C             On exit, if INFO = 0,  DWORK(1)  returns the optimal
C             value of LDWORK, 9*N*NB + 15*NB*NB, where NB is the
C             optimal block size determined by the function ILAHAP.
C             On exit, if  INFO = -16,  DWORK(1)  returns the minimum
C             value of LDWORK.
C
C     LDWORK  INTEGER
C             The length of the array DWORK.  LDWORK >= MAX(1,N).
C
C     Error Indicator
C
C     INFO    INTEGER
C             = 0:  successful exit;
C             < 0:  if INFO = -i, the i-th argument had an illegal 
C                   value.
C
C     METHOD
C
C     The matrix Q is represented as a product of symplectic reflectors
C     and Givens rotators
C
C     Q = diag( H(1),H(1) ) G(1) diag( F(1),F(1) ) 
C         diag( H(2),H(2) ) G(2) diag( F(2),F(2) )
C                           ....
C         diag( H(k),H(k) ) G(k) diag( F(k),F(k) ),
C
C     where k = min(m,n).
C
C     Each H(i) has the form
C
C           H(i) = I - tau * w * w'
C
C     where tau is a real scalar, and w is a real vector with
C     w(1:i-1) = 0 and w(i) = 1; w(i+1:m) is stored on exit in
C     B(i+1:m,i), and tau in B(i,i).
C
C     Each F(i) has the form
C
C           F(i) = I - nu * v * v'
C
C     where nu is a real scalar, and v is a real vector with
C     v(1:i-1) = 0 and v(i) = 1; v(i+1:m) is stored on exit in
C     A(i+1:m,i), and nu in TAU(i).
C
C     Each G(i) is a Givens rotator acting on rows i of A and B,
C     where the cosine is stored in CS(2*i-1) and the sine in
C     CS(2*i).
C
C     REFERENCES
C
C     [1] A. BUNSE-GERSTNER:
C     Matrix factorizations for symplectic QR-like methods.
C     Linear Algebra Appl., 83:49--77, 1986.
C
C     [2] R. BYERS:
C     Hamiltonian and Symplectic Algorithms for the Algebraic
C     Riccati Equation.
C     Ph.D Dissertation, Center for Applied Mathematics,
C     Cornell University, Ithaca, NY, 1983.
C
C     [3] D. KRESSNER:
C     Block algorithms for orthogonal symplectic factorizations.
C     BIT, 43(4):775-790, 2003.
C
C     NUMERICAL ASPECTS
C
C     The algorithm requires O(M*N*N) floating point operations and is
C     numerically backward stable.
C
C     CONTRIBUTORS
C
C     D. Kressner (Technical Univ. Berlin, Germany) and
C     P. Benner (Technical Univ. Chemnitz, Germany), December 2003.
C
C     KEYWORDS
C
C     Elementary matrix operations, orthogonal symplectic matrix.
C
C     ******************************************************************
C
C     .. Parameters ..
      DOUBLE PRECISION  ONE
      PARAMETER         ( ONE = 1.0D0 )
C     .. Scalar Arguments ..
      INTEGER           INFO, LDA, LDB, LDWORK, M, N
C     .. Array Arguments ..
      DOUBLE PRECISION  A(LDA,*), B(LDB,*), CS(*), DWORK(*), TAU(*)
C     .. Local Scalars ..
      INTEGER           I, IB, IERR, K, NB, NBMIN, NX, PDRS, PDT, PDW,
     $                  WRKOPT
C     .. External Functions ..
      INTEGER           ILAHAP
      EXTERNAL          ILAHAP
C     .. External Subroutines ..
      EXTERNAL          DGESQR, DLAESB, DLAEST, XERBLA
C     .. Intrinsic Functions ..
      INTRINSIC         DBLE, MAX, MIN, SQRT
C
C     .. Executable Statements ..
C
C     Check the scalar input parameters.
C
      INFO  = 0
      NB = ILAHAP( 1, 'DGESQB', ' ', M, N, -1 )
      IF ( M.LT.0 ) THEN
         INFO = -1
      ELSE IF ( N.LT.0 ) THEN
         INFO = -2
      ELSE IF ( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -4
      ELSE IF ( LDB.LT.MAX( 1, M ) ) THEN
         INFO = -6
      ELSE IF ( LDWORK.LT.MAX( 1, N ) ) THEN
         DWORK(1) = DBLE(MAX( 1, N ))
         INFO = -10
      END IF
C
C     Return if there were illegal values.
C
      IF ( INFO.NE.0 ) THEN
         CALL XERBLA( 'DGESQB', -INFO )
         RETURN
      END IF
C
C     Quick return if possible.
C
      K = MIN( M, N )
      IF ( K.EQ.0 ) THEN
         DWORK(1) = ONE
         RETURN
      END IF
C
      NBMIN = 2
      NX = 0
      WRKOPT = N
      IF ( NB.GT.1 .AND. NB.LT.K ) THEN
C
C        Determine when to cross over from blocked to unblocked code.
C
         NX = MAX( 0, ILAHAP( 3, 'DGESQB', ' ', M, N, -1 ) )
         IF ( NX.LT.K ) THEN
C
C           Determine if workspace is large enough for blocked code.
C
            WRKOPT = MAX( WRKOPT, 9*N*NB + 15*NB*NB )
            IF( LDWORK.LT.WRKOPT ) THEN
C
C              Not enough workspace to use optimal NB:  reduce NB and
C              determine the minimum value of NB.
C
               NB = INT( ( SQRT( DBLE( 81*N*N + 60*LDWORK ) )
     $                     - DBLE( 9*N ) ) / 30.0D0 )
               NBMIN = MAX( 2, ILAHAP( 2, 'DGESQB',  ' ', M, N, -1 ) )
            END IF
         END IF
      END IF
C
      IF( NB.GE.NBMIN .AND. NB.LT.K .AND. NX.LT.K ) THEN
C
C        Use blocked code initially
C
         PDRS = 1
         PDT = PDRS + 6*NB*NB
         PDW = PDT + 9*NB*NB
C
         DO 10 I = 1, K - NX, NB
            IB = MIN( K-I+1, NB )
C
C           Compute the symplectic QR decomposition of the current
C           blocks [ A(i:m,i:i+ib-1); B(i:m,i:i+ib-1) ].
C
            CALL DGESQR( M-I+1, IB, A(I,I), LDA, B(I,I), LDB, CS(2*I-1),
     $                   TAU(I), DWORK, LDWORK, IERR )
C
            IF( I+IB.LE.N ) THEN
C
C              Form the triangular factors of the symplectic block
C              reflector SH.
C
               CALL DLAEST( 'Forward', 'Columnwise', 'Columnwise',
     $                      M-I+1, IB, A(I,I), LDA, B(I,I), LDB,
     $                      DWORK(PDRS), NB, DWORK(PDT), NB, CS(2*I-1),
     $                      TAU(I), DWORK(PDW) )
C
C              Apply SH' to [ A(i:m,i+ib:n); B(i:m,i+ib:n) ] from the
C              left.
C
               CALL DLAESB( 'No Structure', 'No Transpose',
     $                      'No Transpose', 'Transpose', 'Forward',
     $                      'Columnwise', 'Columnwise', M-I+1, N-I-IB+1,
     $                      IB, A(I,I), LDA, B(I,I), LDB, DWORK(PDRS),
     $                      NB, DWORK(PDT), NB, A(I,I+IB), LDA,
     $                      B(I,I+IB), LDB,  DWORK(PDW) )
            END IF
   10    CONTINUE
      ELSE
         I = 1
      END IF
C
C     Use unblocked code to factor the last or only block.
C
      IF( I.LE.K )
     $   CALL DGESQR( M-I+1, N-I+1, A(I,I), LDA, B(I,I), LDB, CS(2*I-1),
     $                TAU(I), DWORK, LDWORK, IERR )
C
      DWORK(1) = DBLE( WRKOPT )
C
      RETURN
C *** Last line of DGESQB ***
      END
