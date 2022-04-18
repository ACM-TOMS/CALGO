*{SIGMA/AUXSBR/zcmpsv.f}
* 
*      SUBROUTINE ZCMPSV( JOB, M, N, U, LDU, S, X, LDX, E, E1, 
*     $             F, ZW, INFO )
* SIGMA library, AUXSBR section updated February 2016. 
* Developed and coded by Zlatko Drmac, Department of Mathematics
* University of Zagreb, Croatia, drmac@math.hr
* Submitted to ACM TOMS.     
*
*       Purpose
*       ~~~~~~~
*       ZCMPSV compares the refence values of the singular vectors 
*       stored in the columns of U with the aproksimate vectors 
*       stored as the columns of X. The reference singular values 
*       are in S(1:N) and are used to compute the relative gaps in
*       the spectrum and the corresponding error bounds as given
*       in the perturbation theory. 
*       The computed output is
*       E  = max_i { || X(:,i) - U(:,i)*(U(:,i)^H * X(:,i))||_2 *  rgap_i }
*       E1 = max_i { || X(:,i) - U(:,i)*(U(:,i)^H * X(:,i))||_2 *  agap_i }
*       Here rgap_i = MIN( (S(i-1)-S(i))/S(i), (S(i)-S(i+1))/S(i) )
*       where S(1)>=S(2)>=..>=S(N)>=0 are the reference singular values.
*       Further, agap_i = MIN( S(i-1)-S(i), S(i)-S(i+1) )/S(1).
*       Optionally, all errros  
*       || X(:,i) - U(:,i)*(U(:,i)^T * X(:,i))||_2 *  rgap_i
*       are returned in the array [F]. 
*       Only positive singular values are considered.
*
*       Arguments
*       ~~~~~~~~~
*       JOB (input)
*       JOB is CHARACTER
*       = 'A' the errors for each particular vector is computed and all
*             computed relative errors are returned in the array F. The
*             largest errors are returned in E and E1. 
*       = 'M' only the maximal errors are computed and returned in the array F.
*...............................................................................
*       M (input)
*       M is INTEGER
*       The number of rows of the input matrix U.  M >= 0.   
*...............................................................................
*       N (input)
*       N is INTEGER
*       The number of columns of the input matrix U.  N>= 0. 
*...............................................................................
*       [U] (input)
*       [U] is DOUBLE COMPLEX ARRAY, LDU-by-N
*       U(:,i) contains the reference i-th singular vector, i=1:N.
*...............................................................................
*       LDU (input)
*       LDU is INTEGER
*       The leading dimension of [U].
*...............................................................................
*       [S] (input)
*       [S] is DOUBLE PRECISION ARRAY, N-by-1
*       S(i) >= S(i+1), ordered reference singular values; at least one positive
*...............................................................................
*       [X] (input)
*       [X] s DOUBLE COMPLEX ARRAY, LDX-by-N
*       X(:,i) contains the i-th singular vector to be compared to the reference
*       value in U(:,i), i=1:N.
*...............................................................................
*       LDX (input)
*       LDX is INTEGER
*       The leading dimension of [X].
*...............................................................................      
*       E (output)
*       E is DOUBLE PRECISION
*       E  = max_i { || X(:,i) - U(:,i)*(U(:,i)^T * X(:,i))||_2 *  rgap_i }
*...............................................................................
*       E1 (output)
*       E1 is DOUBLE PRECISION
*       E1 = max_i { || X(:,i) - U(:,i)*(U(:,i)^T * X(:,i))||_2 *  agap_i }
*...............................................................................
*       [F] (output)
*       F is DOUBLE PRECISION ARRAY
*       If JOB='A', then 
*       F(i) = || X(:,i) - U(:,i)*(U(:,i)^T * X(:,i))||_2 *  rgap_i, i=1:N.
*       Otherwise, not referenced.
*.......................................................................
*       [ZW] (workspace)
*       ZW is DOUBLE COMPLEX ARRAY 
*       Used as a workspace.
*
*"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
      SUBROUTINE ZCMPSV( JOB, M, N, U, LDU, S, X, LDX, E, E1, 
     $             F, ZW, INFO )      
*       .. Scalar Arguments    
      CHARACTER        JOB 
      INTEGER          M, N, LDU, LDX, INFO 
      DOUBLE PRECISION E, E1
*     .. Array Arguments 	  
      COMPLEX*16       U(LDU,*),  X(LDX, *), ZW(*)
      DOUBLE PRECISION S(*), F(*) 
*......................................................................
*     .. Local Parameters	
      DOUBLE PRECISION ZERO,         ONE
      PARAMETER      ( ZERO = 0.0D0, ONE = 1.0D0 )           
*     .. Local Scalars
      LOGICAL           ALL
      INTEGER           j 
      COMPLEX*16        A 
      DOUBLE PRECISION  B, GAP, GAP1	 
*     .. External Subroutines (BLAS, LAPACK)
      EXTERNAL         ZAXPY, ZCOPY, XERBLA
*     .. External Functions (BLAS,LAPACK)	  
      COMPLEX*16       ZDOTC
      EXTERNAL         ZDOTC 
      DOUBLE PRECISION       DZNRM2
      EXTERNAL               DZNRM2
      LOGICAL          LSAME
      EXTERNAL         LSAME
*     .. Intrinsic Functions 
      INTRINSIC        MAX, MIN
*......................................................................	  
 
      IF ( .NOT.(LSAME(JOB,'A').OR.LSAME(JOB,'M'))) THEN 
        INFO = -1 
      ELSE IF ( M .LT. 0 ) THEN
          INFO = -2 
      ELSE IF ( N .LT. 0 ) THEN
          INFO = -3
      ELSE IF ( LDU .LT. M ) THEN
          INFO = -5
      ELSE IF ( S(1) .LE. ZERO ) THEN
          INFO = -6
      ELSE IF ( LDX .LT. M ) THEN
          INFO = -8
      END IF     
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZCMPSV', -INFO )
         RETURN
      END IF
*
      IF ( LSAME( JOB, 'A') ) THEN 
          ALL = .TRUE. 
      ELSE IF ( LSAME( JOB, 'M') ) THEN
        ALL = .FALSE.
      END IF
      IF ( MIN( M, N ) .EQ. 0 ) THEN
          E  = -ONE
        E1 = -ONE
          RETURN
      END IF
*......................................................................	
      IF ( N .GT. 1 ) THEN 
        E  = ZERO
          E1 = ZERO
        DO 1 j = 1, N
*
           IF ( S(j) .GT. ZERO ) THEN
*
           CALL ZCOPY( M,     X(1,j), 1, ZW, 1 )
           A  = ZDOTC( M,     U(1,j), 1, ZW, 1 )
           CALL ZAXPY( M, -A, U(1,j), 1, ZW, 1 )
           B = DZNRM2( M, ZW, 1 )
           IF ( j .EQ. 1 ) THEN 
              GAP  = (S(1) - S(2)) / (S(1))
              GAP1 = S(1) - S(2)
           ELSE IF ( j .EQ. N ) THEN
              GAP  = (S(N-1) - S(N)) / (S(N))
              GAP1 =  S(N-1) - S(N)
           ELSE
             GAP  = MIN( (S(j-1)-S(j))/(S(j)), (S(j)-S(j+1))/(S(j)) ) 
             GAP1 = MIN( S(j-1)-S(j), S(j)-S(j+1) )
           END IF 
           IF ( GAP1 .LT. ZERO ) THEN
               INFO = - 6 
               CALL XERBLA( 'ZCMPSV', -INFO )
               RETURN
           END IF
           GAP  = MIN( 2.0D0, GAP )
           GAP1 = GAP1 / S(1) 
           E  = MAX( E,  B * GAP  ) 
           E1 = MAX( E1, B * GAP1 )
           IF ( ALL ) F(j) = B * GAP
*
             ELSE IF ( S(j) .LT. ZERO ) THEN
             INFO = - 6 
             CALL XERBLA( 'ZCMPSV', -INFO )
             RETURN
*             
             END IF
*
 1      CONTINUE 
      ELSE
        CALL ZCOPY( M,     X(1,1), 1, ZW, 1 )  
        A  = ZDOTC( M,     U(1,1), 1, ZW, 1 )
      CALL ZAXPY( M, -A, U(1,1), 1, ZW, 1 )
      E  = DZNRM2( M, ZW, 1 )   
      E1 = E 
      IF ( ALL ) F(1) = E 
      END IF
*      
      RETURN
      END