      DOUBLE PRECISION FUNCTION PD2VAL(DEG,DEGMAX,C0,TG1,TG2,TTG1,TTG2)
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
C       This function returns the value of the interpolant at (TG1,TG2).
C     C0 is the matrix of the coefficients computed by PADUA2.
C
C     *********************************************************************
C
C     PARAMS type                i=input, o=output, a=auxiliary           
C     
C     DEG    integer             i: degree of approximation
C     DEGMAX integer             i: maximum degree allowed
C     C0     double(0:DEG,0:DEG) i: coefficient matrix C0
C     TG1    double              i: first coordinate of the target point
C     TG2    double              i: second coordinate of the target point
C     TTG1   double(0:DEG)       a: Normalized Chebyshev polynomials at
C                                   the first coordinate of the target 
C                                   point
C     TTG2   double(0:DEG)       a: Normalized Chebyshev polynomials at
C                                   the second coordinate of the target 
C                                   point
C
C     **********************************************************************
C
      INTEGER DEG,DEGMAX
      DOUBLE PRECISION TG1,TG2,C0(0:DEGMAX + 1,0:DEG),
     1     TTG1(0:DEG),TTG2(0:DEG)
C
C     LOCALS
C
      INTEGER I
      DOUBLE PRECISION DDOT
      EXTERNAL DDOT
C     
C     Compute the normalized Chebyshev polynomials at the target point
C     
      CALL CHEB(DEG,TG1,TTG1)
      CALL CHEB(DEG,TG2,TTG2)
C
C     Evaluate the interpolant
C
      PD2VAL = 0.0D0
      DO 10 I = DEG,0,-1
         PD2VAL = PD2VAL + DDOT(I + 1,TTG1(0),1,C0(0,DEG - I),1) * 
     1        TTG2(DEG - I)
 10   CONTINUE
      RETURN
      END
