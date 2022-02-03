      SUBROUTINE DOSGPV( N, ILO, U1, LDU1, U2, LDU2, CS, TAU, DWORK,
     $                   LDWORK, INFO )
      IMPLICIT NONE
C
C     PURPOSE
C
C     Generates an orthogonal symplectic matrix U, which is defined as
C     a product of symplectic reflectors and Givens rotators
C
C     U = diag( H(1),H(1) )      G(1)  diag( F(1),F(1) ) 
C         diag( H(2),H(2) )      G(2)  diag( F(2),F(2) )
C                                ....
C         diag( H(n-1),H(n-1) ) G(n-1) diag( F(n-1),F(n-1) ).
C
C     as returned by DHAPVL or DSHPVL. The matrix U is returned in terms
C     of its first N/2 rows
C
C                      [  U1   U2 ]
C                  U = [          ].
C                      [ -U2   U1 ]
C
C     Input/Output Parameters
C
C     N       (input) INTEGER
C             The order of the matrices U1 and U2. N >= 0.
C
C     ILO     (input) INTEGER
C             ILO must have the same value as in the previous call of
C             DHAPVL or DSHPVL. U is equal to the unit matrix except in
C             the submatrix
C             U([ilo+1:n n+ilo+1:2*n], [ilo+1:n n+ilo+1:2*n]).
C             1 <= ILO <= N, if N > 0; ILO=1, if N=0.
C
C     U1      (input/output) DOUBLE PRECISION array, dimension (LDU1,N)
C             On entry, the leading N-by-N part of this array must
C             contain in its i-th column the vector which defines the
C             elementary reflector F(i).
C             On exit, the leading N-by-N part of this array contains
C             the matrix Q1.
C
C     LDU1    INTEGER
C             The leading dimension of the array U1.  LDU1 >= N.
C
C     U2      (input/output) DOUBLE PRECISION array, dimension (LDU2,N)
C             On entry, the leading N-by-N part of this array must
C             contain in its i-th column the vector which defines the
C             elementary reflector H(i) and, on the subdiagonal, the
C             scalar factor of H(i).
C             On exit, the leading N-by-N part of this array contains
C             the matrix U2.
C
C     LDU2    INTEGER
C             The leading dimension of the array U2.  LDU2 >= M.
C
C     CS      (input) DOUBLE PRECISION array, dimension (2N-2)
C             On entry, the first 2N-2 elements of this array must
C             contain the cosines and sines of the symplectic Givens
C             rotators G(i).
C
C     TAU     (input) DOUBLE PRECISION array, dimension (N-1)
C             On entry, the first N-1 elements of this array must
C             contain the scalar factors of the elementary reflectors
C             F(i).
C
C     Workspace
C
C     DWORK   DOUBLE PRECISION array, dimension (LDWORK)
C             On exit, if INFO = 0,  DWORK(1)  returns the optimal
C             value of LDWORK.
C             On exit, if  INFO = -10,  DWORK(1)  returns the minimum
C             value of LDWORK.
C
C     LDWORK  INTEGER
C             The length of the array DWORK. LDWORK >= MAX(1,2*(N-ILO)).
C
C     Error Indicator
C
C     INFO    INTEGER
C             = 0:  successful exit;
C             < 0:  if INFO = -i, the i-th argument had an illegal 
C                   value.
C
C     NUMERICAL ASPECTS
C
C     The algorithm requires O(N**3) floating point operations and is
C     strongly backward stable.
C
C     REFERENCES
C
C     [1] C. F. VAN LOAN:
C     A symplectic method for approximating all the eigenvalues of
C     a Hamiltonian matrix.
C     Linear Algebra and its Applications 61 (1984), pp. 233-251.
C
C     [2] D. KRESSNER:
C     Block algorithms for orthogonal symplectic factorizations.
C     To appear in BIT, 2002.
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
      DOUBLE PRECISION  ZERO, ONE
      PARAMETER         ( ZERO = 0.0D+0, ONE = 1.0D+0 )
C     .. Scalar Arguments ..
      INTEGER           ILO, INFO, LDU1, LDU2, LDWORK, N
C     .. Array Arguments ..
      DOUBLE PRECISION  CS(*), DWORK(*), U1(LDU1,*), U2(LDU2,*), TAU(*)
C     .. Local Scalars ..
      INTEGER           I, IERR, J, NH
C     .. External Subroutines ..
      EXTERNAL          DLASET, DOSGSB, XERBLA
C     .. Intrinsic Functions ..
      INTRINSIC         DBLE, MAX
C
C     .. Executable Statements ..
C
C     Check the scalar input parameters.
C
      INFO  = 0
      IF ( N.LT.0 ) THEN
         INFO = -1
      ELSE IF ( ILO.LT.1 .OR. ILO.GT.MAX( 1, N ) ) THEN
         INFO = -2
      ELSE IF ( LDU1.LT.MAX( 1, N ) ) THEN
         INFO = -4
      ELSE IF ( LDU2.LT.MAX( 1, N ) ) THEN
         INFO = -6
      ELSE IF ( LDWORK.LT.MAX( 1, 2*(N-ILO) ) ) THEN
         DWORK(1) = DBLE( MAX( 1, 2*(N-ILO) ) )
         INFO = -10
      END IF
C
C     Return if there were illegal values.
C
      IF ( INFO.NE.0 ) THEN
         CALL XERBLA( 'DOSGPV', -INFO )
         RETURN
      END IF
C
C     Quick return if possible
C
      IF( N.EQ.0 ) THEN
         DWORK(1) = ONE
         RETURN
      END IF
C
C     Shift the vectors which define the elementary reflectors one
C     column to the right, and set the first ilo rows and columns to
C     those of the unit matrix.
C
      DO 30 J = N, ILO + 1, -1
         DO 10 I = 1, J-1
            U1(I,J) = ZERO
   10    CONTINUE     
         DO 20 I = J+1, N
            U1(I,J) = U1(I,J-1)
   20    CONTINUE
   30 CONTINUE
      CALL DLASET( 'All', N, ILO, ZERO, ONE, U1, LDU1 )
      DO 60 J = N, ILO + 1, -1
         DO 40 I = 1, J-1
            U2(I,J) = ZERO
   40    CONTINUE     
         DO 50 I = J, N
            U2(I,J) = U2(I,J-1)
   50    CONTINUE
   60 CONTINUE
      CALL DLASET( 'All', N, ILO, ZERO, ZERO, U2, LDU2 )
      NH = N - ILO
      IF ( NH.GT.0 ) THEN
         CALL DOSGSB( 'No Transpose', 'No Transpose', NH, NH, NH,
     $                U1(ILO+1,ILO+1), LDU1, U2(ILO+1,ILO+1), LDU2,
     $                CS(ILO), TAU(ILO), DWORK, LDWORK, IERR )
      END IF
      RETURN
C *** Last line of DOSGPV ***
      END