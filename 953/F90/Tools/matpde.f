***********************************************************************
*                                                                     *
*     matpde.f: Random matrix generator                               *
*                                                                     *
***********************************************************************
      SUBROUTINE MATPDE( N, NX, NY, A, ROWPTR, COLIND, RHS ) 
      IMPLICIT REAL*8    (A-H,O-Z) 
*     ..
*     .. Scalar Arguments ..
      INTEGER            N, NX, NY
*     .. 
*     .. Array Arguments .. 
      INTEGER            ROWPTR( * ), COLIND( * ) 
      DOUBLE PRECISION   A( * ), RHS( * )
*
*  Purpose
*  =======
*
*  Forms the five-point central differences operator for the elliptic 
*  PDE                    
*
*   -(P Ux)x -(Q Uy)y  +  R Ux + (R U)x  +  S Uy + (S U)y + T U = F
* 
*  where P, Q, R, S and T are the functions of x and y.  The domain 
*  is the unit square (0,1)x(0,1), and the boundary condition are 
*  Dirichlet.                                            
*
*  The matrix is a block tridiagonal matrix, there are NY blocks of
*  size NX by NX on the diagonal (each block is a tridiagonal matrix),
*  and then NY-1 blocks of size NX by NX on the sub- and super-block
*  diagonal positions.  
*
*  Important note: the matrix A is stored in compressed ROW format. 
*
*  Arguments
*  =========
*
*  N       (output) INTEGER 
*          The order of matrix A.  N = NX*NY 
*
*  NX      (input) INTEGER 
*          The number of mesh points in the X-axis
* 
*  NY      (input) INTEGER 
*          The number of mesh points in the Y-axis
*
*  A       (output) DOUBLE PRECISION array, dimension (NZ)  
*          The numerical values of the nonzero elements of matrix. 
*          where NZ is the number of nonzero elements in the generated
*          matrix. Specifically, NZ = NX*NY + (NX-1)*NY*2 + NX*(NY-1)*2. 
*
*  ROWPTR  (output) INTEGER array, dimension (N+1)  
*          The row start pointers.  
*
*  COLIND  (output) INTEGER array, dimension (NZ)  
*          the column indices.
*
*  RHS     (output) DOUBLE PRECISION array, dimension (N)  
*          The right hand side b of the linear system Ax = b. 
*
*
*  Internal parameters 
*  ===================
*
*  BETA, GAMMA are the coefficients of the elliptic operator,
*  BNDRY( 4 ) are the boundary Dirichlet conditions. These parameters
*  may be changed, which will result in different matrices and spectral
*  distribution. 
*                                                                       
*  ==================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE, HALF
      PARAMETER	         ( ZERO = 0.0D+0, ONE = 1.0D+0, HALF = 0.5D+0 ) 
      DOUBLE PRECISION   BETA, GAMMA
      PARAMETER          ( BETA = 2.0D+1, GAMMA = 0.0D+0 )
*     .. 
*     .. Scalar variables .. 
      INTEGER		 I, J, IX, JY, INDEX, JD, IAI 
*
*     .. Array variables .. 
      DOUBLE PRECISION   BNDRY( 4 ) 
*     ..
*     .. External functions .. 
      DOUBLE PRECISION   F, PC, QC, RC, SC, U
      EXTERNAL           F, PC, QC, RC, SC, U
*     ..
*     .. Intrinsic functions .. 
      INTRINSIC	         MOD, DBLE
*                                                                       
*     set up initial boundary conditions 
*
      BNDRY( 1 ) = ZERO  
      BNDRY( 2 ) = ZERO   
      BNDRY( 3 ) = ZERO   
      BNDRY( 4 ) = ZERO   
*
*     set up other parameters 
*
      N = NX*NY 
      HX = ONE / (DBLE(NX) + ONE) 
      HY = ONE / (DBLE(NY) + ONE)
      HY2 = HY*HY                                                     
      RA = (HY/HX)**2                                                 
      RB  = HY2/HX                                                    
*                                                                       
*     FORM ROWPTR.                                                             
*
      ROWPTR(1) = 1                                                       
      DO 10 I = 1,N                                                     
         IAI = 1                                                      
         IF (I.GT.NX) 
     $      IAI = IAI+1                                      
         IF(MOD(I-1,NX).NE.0) 
     $     IAI = IAI+1                             
         IF(MOD(I,NX).NE.0) 
     $     IAI = IAI+1                               
         IF(I.LE.N-NX) 
     $     IAI = IAI+1                                    
         ROWPTR(I+1) = ROWPTR(I) + IAI                                        
 10   CONTINUE                                                        
*                                                                       
      J = 1                                                           
      DO 110 JY=1,NY                                                  
         YJ = HY*DBLE(JY)
         DO 100 IX=1,NX                                               
*
            XI = HX*DBLE(IX) 
            INDEX = IX + (JY-1)*NX                                    
            RHS(INDEX) = HY2*F(BETA,GAMMA,XI,YJ) 
            P12  = PC (XI + HALF*HX,YJ) 
            PM12 = PC (XI - HALF*HX,YJ)  
            Q12  = QC (XI,YJ + HALF*HY)
            QM12 = QC (XI,YJ - HALF*HY) 
            RIJ  = RC (BETA,XI,YJ)                                         
            R1   = RC (BETA,XI + HX,YJ)                                    
            RM1  = RC (BETA,XI - HX,YJ)                                    
            SIJ  = SC (GAMMA,XI,YJ)                                         
            S1   = SC (GAMMA,XI,YJ + HY)                                    
            SM1  = SC (GAMMA,XI,YJ - HY)                                    
*                                                                       
*           DIAGONAL.   
*
            A(J) = RA*(P12 + PM12)  +  Q12 + QM12  + HY2*TC(XI,YJ)    
            COLIND(J) = INDEX                                             
            JD = J                                                    
            J = J+1                                                   
*                                                                       
*           LOWEST BAND. 
*
            IF (JY.EQ.1) 
     $         GOTO 20                                      
            A(J) = - ( QM12 + HALF*HY*(SIJ+SM1) )
            COLIND(J) = INDEX - NX                                     
            J = J+1                                                
*                                                                       
*           SUB-DIAGONAL.  
*
 20         IF (IX.EQ.1) 
     $         GOTO 30                                      
            A(J) = - ( RA*PM12 + HALF*RB*(RIJ+RM1) ) 
            COLIND(J) = INDEX - 1                                      
            J = J+1                                                
*                                                                       
*           SUPER-DIAGONAL.  
*
 30         IF (IX.EQ.NX) 
     $         GOTO 40                                     
            A(J) = -RA*P12 + HALF*RB*(RIJ+R1) 
            COLIND(J) = INDEX + 1                                      
            J = J+1                                                
*                                                                       
*           HIGHEST BAND.  
*
 40         IF (JY.EQ.NY) 
     $         GOTO 50                                     
            A(J) = -Q12 + HALF*HY*(SIJ+S1) 
            COLIND(J) = INDEX + NX                                     
            J = J+1                                                
*                                                                       
 50         CONTINUE                                                  
*
*           BOUNDARY CONDITIONS.     
*
            IF (JY.NE.1) 
     $         GOTO 60                                      
            COEF = - ( QM12 + HALF*HY*(SIJ+SM1) )
            IF (BNDRY(1).EQ.0)                                     
     $         RHS(INDEX) = RHS(INDEX) - COEF*U(XI,YJ-HY)          
            IF (BNDRY(1).EQ.1)  
     $         A(JD) = A(JD) + COEF               
*                                                                       
 60         IF (IX.NE.1) GOTO 70                                     
               COEF = - ( RA*PM12 + HALF*RB*(RIJ+RM1) ) 
            IF (BNDRY(2).EQ.0)                                     
     $         RHS(INDEX) = RHS(INDEX) - COEF*U(XI-HX,YJ)          
            IF (BNDRY(2).EQ.1)  
     $         A(JD) = A(JD) + COEF               
*                                                                       
 70         IF (IX.NE.NX) GOTO 80                                     
               COEF = -RA*P12 + HALF*RB*(RIJ+R1) 
            IF (BNDRY(3).EQ.0)                                     
     $         RHS(INDEX) = RHS(INDEX) - COEF*U(XI+HX,YJ)          
            IF (BNDRY(3).EQ.1)  
     $         A(JD) = A(JD) + COEF               
*                                                                       
 80         IF (JY.NE.NY) 
     $         GOTO 100                                    
            COEF = -Q12 + HALF*HY*(SIJ+S1) 
            IF (BNDRY(4).EQ.0)                                     
     $         RHS(INDEX) = RHS(INDEX) - COEF*U(XI,YJ+HY)          
            IF (BNDRY(4).EQ.1)  
     $         A(JD) = A(JD) + COEF               
*                                                                       
 100     CONTINUE                                                     
*
 110  CONTINUE                                                        
*                                                                       
      RETURN                                                          
*
*     End of MATPDE
*
      END                                                             
*                                                                       
      DOUBLE PRECISION FUNCTION  PC (X,Y)
      DOUBLE PRECISION   X, Y 
      INTRINSIC          EXP
      PC = EXP(-X*Y)                                                  
      RETURN                                                          
      END                                                             
*                                                                       
      DOUBLE PRECISION FUNCTION  QC (X,Y)
      DOUBLE PRECISION   X, Y 
      INTRINSIC          EXP
      QC = EXP(X*Y)                                                   
      RETURN                                                          
      END                                                             
*
      DOUBLE PRECISION FUNCTION TC (X,Y) 
      DOUBLE PRECISION   X, Y 
      DOUBLE PRECISION   ONE
      PARAMETER          ( ONE = 1.0D+0 ) 
      TC = ONE / (ONE + X + Y)
      RETURN                                                          
      END                                                             
*
      DOUBLE PRECISION FUNCTION  U (X,Y) 
      DOUBLE PRECISION   X, Y
      DOUBLE PRECISION   PI
      PARAMETER          (PI = 3.14159265358979D+0)
      INTRINSIC          SIN, EXP
      U = X * EXP(X*Y) * SIN(PI*X) * SIN(PI*Y)                        
      RETURN                                                          
      END                                                             
*                                                                       
      DOUBLE PRECISION FUNCTION  RC (BETA,X,Y) 
      DOUBLE PRECISION   BETA, X, Y 
      RC = BETA * (X+Y)                                           
      RETURN                                                          
      END                                                             
*                                                                       
      DOUBLE PRECISION FUNCTION  SC (GAMMA,X,Y)
      DOUBLE PRECISION   GAMMA, X, Y 
      SC = GAMMA * (X+Y)                                           
      RETURN                                                          
      END                                                             
*                                                                       
      DOUBLE PRECISION FUNCTION F (BETA, GAMMA, X,Y) 
      DOUBLE PRECISION   BETA, GAMMA, X, Y
*
*
      DOUBLE PRECISION   PI
      PARAMETER          (PI = 3.14159265358979D+0) 
      DOUBLE PRECISION   ONE, TWO
      PARAMETER          (ONE = 1.0D+0, TWO = 2.0D+0) 
*
      DOUBLE PRECISION   PXY, PXXY, QXY, QYXY, CXY, RXY, RXXY,
     $                   SXY, SYXY, EXY, SX, CX, SY, CY, A, AX,
     $                   AXX, AY, AYY
*
      DOUBLE PRECISION   PC, QC, RC, SC
      EXTERNAL	         PC, QC, RC, SC
      INTRINSIC          COS, SIN, EXP
*                                                                       
      PXY  = PC(X,Y)                                                  
      PXXY = -Y * EXP(-X*Y)                                           
      QXY  = QC(X,Y)                                                  
      QYXY = X * EXP(X*Y)                                             
      CXY  = ONE / ( ONE + X + Y)                                        
      RXY  = RC(BETA,X,Y)                                                  
      RXXY = BETA 
      SXY  = SC(GAMMA,X,Y)                                                  
      SYXY = GAMMA
*                                                                       
      EXY = EXP(X*Y)                                                  
      SX  = SIN(PI*X)                                                 
      CX  = COS(PI*X)                                                 
      SY  = SIN(PI*Y)                                                 
      CY  = COS(PI*Y)                                                
*                                                                       
      A   = X * EXY * SX * SY                                         
      AX  = EXY * (PI * X * CX + (ONE + X * Y) * SX) * SY 
      AXX = EXY * (PI * (-PI * X * SX + (X * Y + TWO) * CX) + Y * SX)
     $        * SY + Y * AX                                             
      AY  = EXY * (PI * CY + X * SY) * X * SX                         
      AYY = X * (EXY * PI * (-PI * SY + X * CY) * SX + AY)            
*                                                                       
      F = - (PXY*AXX + QXY*AYY)  +                                    
     $      (TWO*RXY - PXXY) * AX  +  (TWO*SXY - QYXY) * AY  + 
     *      (RXXY + SYXY + CXY) * A                                     
*
      RETURN                                                          
      END                                                             
