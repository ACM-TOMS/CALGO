*{SIGMA/AUXSBR/xgec2z.f}
*
      SUBROUTINE XGEC2Z( M, N, CA, LDCA, ZA, LDZA )
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
      INTEGER        M, N, LDCA, LDZA  
*     .. Array Arguments ..      
      COMPLEX        CA( LDCA, * )
      COMPLEX*16 ZA( LDZA, * )
*.......................................................................
*     .. Local Scalars .. 
      INTEGER i, INFO, j    
*     .. External Subroutines (BLAS, LAPACK)
      EXTERNAL  XERBLA
*     .. Intrinsic Functions
      INTRINSIC CMPLX 
*.......................................................................
      INFO = 0 
      IF ( M .LT. 0 ) THEN
          INFO = 1
      ELSE IF ( N .LT. 0 ) THEN
          INFO = 2
      ELSE IF ( LDCA .LT. M ) THEN 
          INFO = 4 
      ELSE IF ( LDZA .LT. M ) THEN 
          INFO = 6
      END IF 
      IF (INFO.NE.0) THEN
          CALL XERBLA('XGEC2Z',INFO)
          RETURN
      END IF
      DO 1 j = 1, N 
         DO 2 i = 1, M 
            ZA(i,j) = CMPLX(CA(i,j),KIND=KIND(0.0D0))
 2       CONTINUE 
 1    CONTINUE
      RETURN
      END 