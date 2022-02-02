*{SIGMA/AUXSBR/xges2d.f}
*
      SUBROUTINE XGES2D( M, N, SA, LDSA, DA, LDDA )
*
* SIGMA library, AUXSBR section updated October 2015. 
* Developed and coded by Zlatko Drmac, Department of Mathematics
* University of Zagreb, Croatia, drmac@math.hr
* Submitted to ACM TOMS.       
*
*     Purpose
*     ~~~~~~~
*     XGEC2Z copies COMPLEX M x N matrix CA(LDCA,*) into the 
*            DOUBLE COMPLEX M x N matrix ZA(LDZA,*).
*
*     .. Scalar Arguments ..     
      INTEGER        M, N, LDSA, LDDA  
*     .. Array Arguments ..      
      REAL             SA( LDSA, * )
      DOUBLE PRECISION DA( LDDA, * )
*.......................................................................
*     .. Local Scalars .. 
      INTEGER i, INFO, j    
*     .. External Subroutines (BLAS, LAPACK)
      EXTERNAL  XERBLA
*     .. Intrinsic Functions
      INTRINSIC DBLE 
*.......................................................................
      INFO = 0 
      IF ( M .LT. 0 ) THEN
          INFO = 1
      ELSE IF ( N .LT. 0 ) THEN
          INFO = 2
      ELSE IF ( LDSA .LT. M ) THEN 
          INFO = 4 
      ELSE IF ( LDDA .LT. M ) THEN 
          INFO = 6
      END IF 
      IF (INFO.NE.0) THEN
          CALL XERBLA('XGES2D',INFO)
          RETURN
      END IF
      DO 1 j = 1, N 
         DO 2 i = 1, M 
            DA(i,j) = DBLE(SA(i,j))
 2       CONTINUE 
 1    CONTINUE
      RETURN
      END 