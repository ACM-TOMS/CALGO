*{SIGMA/AUXSBR/dchkor.f}
*
      SUBROUTINE DCHKOR( JOB, M, N, V, LDV, DMIN, DMAX, OFFMAX, 
     $           X, LDX, INFO )
* SIGMA library, AUXSBR section updated February 2016. 
* Developed and coded by Zlatko Drmac, Department of Mathematics
* University of Zagreb, Croatia, drmac@math.hr
* Submitted to ACM TOMS.   
*
*     Purpose
*     ~~~~~~~
*     DCHKOR checks whether V is left or right numerically orthogonal.
*     If JOBU = 'T' ('N'), it computes X = V^T * V (X=V * V^T) and then 
*     DMIN   = min_i ABS(X(i,i))     -> should be close to one
*     DMAX   = max_i ABS(X(i,i))     -> should be close to one
*     OFFMAX = max_{i<j} ABS(X(i,j)) -> should be at roundoff level 
*
      IMPLICIT NONE
*     .. Scalar Arguments
      CHARACTER           JOB
      INTEGER             M, N, LDV,  LDX, INFO
      DOUBLE PRECISION    DMIN, DMAX, OFFMAX 
*     .. Array Arguments       
      DOUBLE PRECISION      V( LDV, * )
      DOUBLE PRECISION      X( LDX, * )
*.......................................................................
*=======================================================================
*     .. Local Parameters	
      DOUBLE PRECISION ZERO,         ONE
      PARAMETER      ( ZERO = 0.0D0, ONE = 1.0D0 )
*     .. Local Scalars
      DOUBLE PRECISION TEMP, TEMP1, TEMP2
      INTEGER          i, j, K 
*     .. External Subroutines (BLAS, LAPACK)
      EXTERNAL         DSYRK, XERBLA 
*     .. External Function (BLAS, LAPACK)
      LOGICAL          LSAME
      EXTERNAL         LSAME
*     .. Intrinsic Functions
      INTRINSIC        ABS, MIN, MAX 
*......................................................................	
      INFO = 0 
      IF ( .NOT.( LSAME(JOB,'T') .OR. LSAME(JOB,'N') ) ) THEN 
          INFO = -1
      ELSE IF ( LSAME(JOB,'T') .AND. ( M.LT.N) ) THEN
          INFO = -1
      ELSE IF ( LSAME(JOB,'N') .AND. ( N.LT.M ) ) THEN
          INFO = -1
      ELSE IF ( M.LT.0 ) THEN
          INFO = -2
      ELSE IF ( N.LT.0 ) THEN
          INFO = -3
      ELSE IF ( LDV.LT.M ) THEN
          INFO = -5
      ELSE IF ( LSAME(JOB,'T') .AND. (LDX.LT.N ) ) THEN
          INFO = -10
      ELSE IF ( LSAME(JOB,'N') .AND. (LDX.LT.M ) ) THEN
          INFO = -10 
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DCHKOR', -INFO )
         RETURN
      END IF
*......................................................................	   
      IF ( MIN( M, N ) .EQ. 0 ) THEN
          DMAX   = -ONE
          DMIN   = -ONE
          OFFMAX = -ONE
          RETURN
      END IF
*......................................................................
      IF ( LSAME( JOB, 'T' ) ) THEN 
*        .. test V^T * V  against the N x N identity         
         CALL DSYRK( 'U', 'T', N, M, ONE, V, LDV, ZERO, X, LDX )
         K = N 
      ELSE
*        .. test V * V^T against the M x M identity
         CALL DSYRK( 'U', 'N', M, N, ONE, V, LDV, ZERO, X, LDX )
         K = M 
      END IF
*     
      TEMP = -ONE
      TEMP1 = ABS(X(K,K))
      TEMP2 = ABS(X(K,K))
      DO 99 i = 1, K - 1
          TEMP1 = MAX(TEMP1,ABS(X(i,i)))
          TEMP2 = MIN(TEMP2,ABS(X(i,i)))
          DO 100 j = i + 1, K
              TEMP = MAX(TEMP,ABS(X(i,j)))
 100      CONTINUE
 99	CONTINUE	
      DMAX   = TEMP1
      DMIN   = TEMP2
      OFFMAX = TEMP	
*
      RETURN
      END
