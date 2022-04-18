      SUBROUTINE PDPTS(DEG,PD1,PD2,WPD,NPD)
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
C       This subroutine computes the first family of Padua points and 
C     weights corresponding to degree DEG.
C
C     ******************************************************************
C
C     PARAMS type                i=input, o=output, a=auxiliary           
C     
C     DEG    integer             i: degree of approximation
C     PD1    double(*)           o: array of the first coordinates of 
C                                   the Padua points
C     PD2    double(*)           o: array of the second coordinates of 
C                                   the Padua points
C     WPD    double(*)           o: array of the weights
C     NPD    integer             o: number of Padua points = 
C                                   (DEG + 1) * (DEG + 2) / 2
C
C     ******************************************************************
C
      INTEGER DEG,NPD
      DOUBLE PRECISION PD1(*),PD2(*),WPD(*)
C
C     LOCALS
C
      INTEGER J,K,DEG1,ITEMP0
      DOUBLE PRECISION PI,RTEMP0
      PARAMETER (PI = 3.1415926535897931D0)
C
C     Compute the Padua points of the first family at degree DEG
C
      IF (DEG .EQ. 0) THEN
C     
C     DEG .EQ. 0
C     
         PD1(1) = -1.D0
         PD2(1) = -1.D0
         WPD(1) = 2.0D0
         NPD = 1
         RETURN
      END IF   
      NPD = 0
      DEG1 = DEG + 1
      ITEMP0 = DEG * DEG1
      RTEMP0 = PI / ITEMP0
      DO 10 J = 0,DEG1
         DO 20 K = MOD(J,2),DEG,2
            NPD = NPD + 1
            PD1(NPD) = -COS(DEG1 * K * RTEMP0)
            PD2(NPD) = -COS(DEG * J * RTEMP0)
            WPD(NPD) = 2.0D0 / ITEMP0
            IF ((K .EQ. 0) .OR. (K .EQ. DEG)) THEN
               WPD(NPD) = WPD(NPD) / 2
            END IF
            IF ((J .EQ. 0) .OR. (J .EQ. DEG1)) THEN
               WPD(NPD) = WPD(NPD) / 2
            END IF
 20      CONTINUE
 10   CONTINUE
      RETURN
      END
