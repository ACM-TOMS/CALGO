      SUBROUTINE PADUA2(DEG,DEGMAX,NPD,WPD,FPD,RAUX1,RAUX2,C0,ESTERR)
C
C     ******************************************************************
C
C                                                           From Padua2D
C
C                                    Marco Caliari and Stefano De Marchi
C                                         Department of Computer Science
C                                           University of Verona (Italy)
C                              {marco.caliari,stefano.demarchi}@univr.it
C
C                                                         Marco Vianello
C                             Department of Pure and Applied Mathematics
C                                            University of Padua (Italy)
C                                                   marcov@math.unipd.it
C                                                    November 14th, 2007
C
C       This subroutine computes the coefficient matrix C0, in the 
C     orthonormal Chebyshev basis T_j(x)T_{k-j}(y), 0 <= j <= k <= DEG, 
C     T_0(x)=1, T_j(x) = sqrt(2) * cos(j * acos(x)), of the 
C     interpolation polynomial of degree DEG of the function values FPD 
C     at the set of NPD Padua points (PD1,PD2) in the square [-1,1]^2. 
C     The interpolant may be evaluated at an arbitrary point by the 
C     function PD2VAL. PD1, PD2 and WPD are the Padua points and weights 
C     computed by PDPTS.
C
C     ******************************************************************
C
C     PARAMS type                   i=input, o=output, a=auxiliary           
C     
C     DEG    integer                i: degree of approximation
C     DEGMAX integer                i: maximum degree allowed
C     NPD    integer                i: number of Padua points
C     WPD    double(NPD)            i: array of the weights
C     FPD    double(NPD)            i: array of the values of the function  
C                                      at the Padua points
C     RAUX1  double(0:DEGMAX,DEG+2) a: auxiliary matrix
C     RAUX2  double(0:DEGMAX,DEG+2) a: auxiliary matrix
C     C0     double(0:DEG,0:DEG)    o: coefficient matrix C0
C     ESTERR double                 o: estimated error
C
C     ******************************************************************
C
      INTEGER DEG,DEGMAX,NPD
      DOUBLE PRECISION ESTERR
      DOUBLE PRECISION WPD(NPD),FPD(NPD),
     1     RAUX1(0:DEGMAX,DEG + 2),RAUX2(0:DEGMAX,DEG + 2),
     2     C0(0:DEGMAX + 1,0:DEG + 1)
C
C     LOCALS
C
      INTEGER I,J,K,DEG1
      DOUBLE PRECISION PI
      PARAMETER (PI = 3.1415926535897931D0)
      DEG1 = DEG + 1
C
C     Build the matrix P_2 and store it in RAUX2
C
      DO 10 I = 0,DEG1
         CALL CHEB(DEG,-COS(I * PI / DEG1),RAUX2(0,I + 1))
 10   CONTINUE
C
C     Build the matrix G(f) and store it in C0
C
      K = 0
      DO 20 I = 0,DEG1
         DO 30 J = 0,DEG
            IF (MOD(I+J,2) .EQ. 0) THEN
               K = K + 1
               C0(J,I) = FPD(K) * WPD(K)
            ELSE
               C0(J,I) = 0.0D0
            END IF
 30      CONTINUE
 20   CONTINUE
C
C     Compute the matrix-matrix product G(f)*P_2' and store it in RAUX1
C
      CALL DGEMM('N','T',DEG1,DEG1,DEG1 + 1,1.0D0,
     1        C0,DEGMAX + 2,RAUX2,DEGMAX + 1,0.0D0,RAUX1,DEGMAX + 1)
C
C     Build the matrix P_1 and store it in RAUX2
C
      DO 40 I = 0,DEG
         CALL CHEB(DEG,-COS(I * PI / DEG),RAUX2(0,I + 1))
 40   CONTINUE
C
C     Compute the matrix-matrix product C(f) = P_1*(G(f)*P_2') 
C     and store it in C0
C
      CALL DGEMM('N','N',DEG1,DEG1,DEG1,1.0D0,
     1        RAUX2,DEGMAX + 1,RAUX1,DEGMAX + 1,0.0D0,C0,DEGMAX + 2)
      C0(DEG,0) = C0(DEG,0) / 2
C
C     Compute the estimated error
C
      ESTERR = 0.0D0
      DO 200 J = 0,2
         DO 300 I = 0,DEG - J
            ESTERR = ESTERR + ABS(C0(I,DEG - I - J))
 300      CONTINUE
 200   CONTINUE
      ESTERR = 2 * ESTERR
      END
