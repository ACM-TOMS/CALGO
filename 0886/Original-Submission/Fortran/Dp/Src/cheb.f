      SUBROUTINE CHEB(DEG,PT,TCHEB)
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
C       This subroutine computes the array TCHEB of normalized Chebyshev 
C     polynomials from degree 0 to DEG, T_0(x)=1, 
C     T_j(x) = sqrt(2) * cos(j * acos(x)) at the point x = PT.
C
C     ******************************************************************
C
C     PARAMS type                i=input, o=output, a=auxiliary           
C     
C     DEG    integer             i: degree
C     PT     double              i: point
C     TCHEB  double(0:DEG)       o: Normalized Chebyshev polynomials at
C                                   the point PT
C
C     **********************************************************************
      INTEGER DEG
      DOUBLE PRECISION PT,TCHEB(0:DEG)
C
C     LOCALS
C
      DOUBLE PRECISION SQRT2
      PARAMETER (SQRT2 = 1.4142135623730951D0)
      INTEGER J
      IF (DEG .GE. 0) THEN
            TCHEB(0) = 1.0D0
      END IF      
      IF (DEG .GE. 1) THEN
         TCHEB(1) = SQRT2 * PT
      ELSE
         RETURN
      END IF
      IF (DEG .GE. 2) THEN
         TCHEB(2) = 2 * PT * TCHEB(1) - SQRT2 * TCHEB(0)
      ELSE
         RETURN
      END IF
C
C     Chebyshev recurrence
C
      DO 10 J = 3,DEG
         TCHEB(J) = 2 * PT * TCHEB(J - 1) - TCHEB(J - 2)
 10   CONTINUE
      RETURN
      END
