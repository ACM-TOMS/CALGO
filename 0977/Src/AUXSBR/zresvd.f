*{SIGMA/AUXSBR/zresvd.f}
*
*  Definition:
*  ===========
*      DOUBLE PRECISION FUNCTION ZRESVD( M, N, K, A, LDA, 
*     $        U, LDU, V, LDV, VVT, S, NRM, DWORK )
* SIGMA library, AUXSBR section updated February 2016. 
* Developed and coded by Zlatko Drmac, Department of Mathematics
* University of Zagreb, Croatia, drmac@math.hr
* Submitted to ACM TOMS.  
*
*      IMPLICIT NONE
*      .. Scalar Arguments .. 
*      INTEGER        M, N, K, LDA, LDU, LDV
*      CHARACTER     VVT, NRM
*      .. Array Arguments ..
*      COMPLEX*16       A( LDA, * ), U( LDU, * ), V( LDV, * )
*      DOUBLE PRECISION S( * ), DWORK( * )      
*
*  Purpose
*  ~~~~~~~
*  ZRESVD computes the residual of a computed SVD of A.
*  ZRESDV = NORM( A - U * DIAG(S) * V^H ) / NORM(A).
*  Here NORM is the matrix norm as computed by DLANGE() with
*  the norm specified in the character NRM. 
*
*  Arguments
*  ~~~~~~~~~
*...............................................................................
*  M (input)
*  M is INTEGER
*  The number of rows of the matrix A.  M >= 0.   
*...............................................................................
*  N (input)
*  N is INTEGER
*  The number of columns of the matrix A.  N>= 0. 
*............................................................................... 
*  K (input)
*  K is INTEGER
*  The rank of the approximate SVD (number of nonzero singular values). K>=0 
*...............................................................................      
*  A (input/output)
*  A DOUBLE COMPLEX array of dimensions LDA x N
*  On entry, the input matrix A. On exit, the residual A - U*DIAG(S)*V^H  
*...............................................................................
*  LDA (input)
*  DA is INTEGER
*  The leading dimension of [A].
*............................................................................... 
*  [U] (input)
*  [U] is DOUBLE COMPLEX ARRAY, LDU-by-N
*  U(:,i) contains the computed i-th left singular vector, i=1:N.
*...............................................................................
*  LDU (input)
*  LDU is INTEGER
*  The leading dimension of [U].
*...............................................................................
*  VVT (input)
*  VVT is CHARACTER*1
*  Specifies whether the right singular vectord are stored as columns or rows of
*  the array [V]. See the description of [V].
*...............................................................................
*  [V] (input)
*  [V] is DOUBLE COMPLEX ARRAY, LDU-by-N
*  If VVT ='N', V(:,i) contains the computed i-th right singular vector, i=1:N.
*  If VVT ='C', V(i,:) is the computed conjugated i-th right singular vector.      
*...............................................................................
*  LDV (input)
*  LDV is INTEGER
*  The leading dimension of [V]. 
*  If VVT ='N', then LDV >= N. If VVT = 'T', then LDV >= K.
*...............................................................................   
*  [S] (input)
*  [S] is DOUBLE PRECISION ARRAY, N-by-1
*  S(i)= i-th computed  singular value;
*...............................................................................  
*  NRM (input)
*  NRM is CHARACTER*1
*  Specifies the norm to be used by DLANGE, 'F', 'M', 'I', '1'. For instance, 
*  NRM='F' defines the norm to be used as the Frobenius norm. 
*  See the description of DLANGE (LAPACK).
*...............................................................................
*  [DWORK] (work space)
*  [DWORK] is DOUBLE PRECISION array of length at least M, that is used as work 
*  space if NRM='I'; otherwise, DWORK is not referenced.
*...............................................................................
*"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""    
      DOUBLE PRECISION FUNCTION ZRESVD( M, N, K, A, LDA, 
     $        U, LDU, V, LDV, VVT, S, NRM, DWORK )
      IMPLICIT NONE
*     .. Scalar Arguments .. 
      INTEGER        M, N, K, LDA, LDU, LDV
      CHARACTER     VVT, NRM
*     .. Array Arguments ..
      COMPLEX*16       A( LDA, * ), U( LDU, * ), V( LDV, * )
      DOUBLE PRECISION S( * ), DWORK( * )
*
*=======================================================================
*     .. Local Parameters ..      
      COMPLEX*16         ONE
      PARAMETER        ( ONE = (1.0D+0,0.0D0) )
      DOUBLE PRECISION   ZERO
      PARAMETER        ( ZERO = 0.0D0 ) 
*     .. Local Variables      
      DOUBLE PRECISION   AFN, TMP 
      INTEGER            i, INFO      
*     .. External Functions (BLAS, LAPACK) 
      EXTERNAL           LSAME,  ZLANGE 
      LOGICAL            LSAME
      DOUBLE PRECISION           ZLANGE
*     .. External subroutines (BLAS, LAPACK)     
      EXTERNAL           ZDSCAL, ZGEMM   
*     .. Intrinsic Functions
      INTRINSIC          MIN
*......................................................................
      INFO = 0 
      IF ( M .LT. 0 ) THEN 
          INFO = -1
      ELSE IF ( N .LT. 0 ) THEN
          INFO = -2
      ELSE IF ( K .LT. 0 ) THEN
          INFO = -3
      ELSE IF ( LDA .LT. M ) THEN
          INFO = -5
      ELSE IF ( LDU .LT. M ) THEN
          INFO = -7
      ELSE IF ( .NOT.( LSAME(VVT,'N').OR.LSAME(VVT,'C') ) ) THEN
          INFO = -8
      ELSE IF ( ( LSAME(VVT,'N') .AND. (LDV.LT.N) ) .OR.
     $          ( LSAME(VVT,'C') .AND. (LDV.LT. K) ) ) THEN
          INFO = -10
      ELSE IF ( .NOT.(LSAME(NRM,'F').OR.LSAME(NRM,'M')
     $    .OR.LSAME(NRM,'I').OR.LSAME(NRM,'1')) ) THEN
          INFO = -12
      END IF
      IF( INFO.NE.0 ) THEN
         ZRESVD = -ONE 
         CALL XERBLA( 'ZRESVD', -INFO )
         RETURN
      END IF
*
      IF ( MIN( M, N, K ) .EQ. 0 ) THEN
          ZRESVD = -ABS(ONE)
          RETURN
      END IF
*......................................................................
*..       
      DO 1 i = 1 , N 
          CALL ZDSCAL( M, S(i), U(1,i), 1 )
 1    CONTINUE
*      
      AFN = ZLANGE( NRM, M, N, A, LDA, DWORK )
      
      IF ( LSAME( VVT, 'N' ) ) THEN 
         CALL ZGEMM( 'N', 'C', M, N, K, -ONE, U, LDU, V, LDV, 
     $            ONE, A, LDA )
      ELSE
         CALL ZGEMM( 'N', 'N', M, N, K, -ONE, U, LDU, V, LDV, 
     $            ONE, A, LDA )
      END IF
*           
      TMP   = ZLANGE( NRM, M, N, A, LDA, DWORK )
      IF ( (AFN .EQ. ZERO) .AND. (TMP .EQ. ZERO) ) THEN 
          ZRESVD = ZERO
      ELSE
          ZRESVD =  TMP / AFN
      END IF
*      
      RETURN
      END
  