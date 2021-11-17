*{SIGMA/AUXSBR/xged2s.f}
*
      SUBROUTINE XGED2S( M, N, DA, LDDA, SA, LDSA )
*
* SIGMA library, AUXSBR section updated February 2016. 
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
      INTRINSIC REAL 
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
          CALL XERBLA('XGED2S',INFO)
          RETURN
      END IF
      DO 1 j = 1, N 
         DO 2 i = 1, M 
            SA(i,j) = REAL(DA(i,j))
 2       CONTINUE 
 1    CONTINUE
      RETURN
      END 