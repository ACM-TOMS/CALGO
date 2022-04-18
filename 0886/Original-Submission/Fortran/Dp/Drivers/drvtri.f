      PROGRAM DRVTRI
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
C       This driver computes the interpolation of the Franke function
C     on the triangle T(U,V,W) with vertices U=(U1,U2)=(0,0), 
C     V=(V1,V2)=(1,0) and W=(W1,W2)=(0,1) (unit triangle) 
C     at the first family of Padua points. 
C     The degree of interpolation is DEG = 60 and the number of target 
C     points is NTG = NTG1 ** 2 - NTG1 + 1, NTG1 = 100.
C     The maps from the reference square [-1,1]^2 to the triangle
C     are SIGMA1 and SIGMA2 with inverses ISIGM1 and ISIGM2.
C
C     ******************************************************************
C
C     LOCALS type                
C     
C     DEGMAX integer                   maximum degree of interpolation
C     NPDMAX integer                   maximum number of Padua points
C                                      = (DEGMAX + 1) * (DEGMAX + 2) / 2
C     NTG1MX integer                   maximum value of the parameter
C                                      determining the number of target 
C                                      points
C     NTGMAX integer                   maximum number of target points,
C                                      dependent on NTG1MX
C     DEG    integer                   degree of interpolation
C     NTG1   integer                   parameter determining the number 
C                                      of target points
C     NPD    integer                   number of Padua points
C                                      = (DEG + 1) * (DEG + 2) / 2
C     NTG    integer                   number of target points,
C                                      dependent on NTG1
C     PD1    double(NPDMAX)            array of the first coordinates 
C                                      of the Padua points
C     PD2    double(NPDMAX)            array of the second coordinates
C                                      of the Padua points
C     WPD    double(NPDMAX)            array of the weights
C     FPD    double(NPDMAX)            values of the function at 
C                                      the Padua points
C     RAUX1  double((DEGMAX+1)*(DEGMAX+2)) auxiliary array
C     RAUX2  double((DEGMAX+1)*(DEGMAX+2)) auxiliary array
C     C0     double(0:DEGMAX+1,0:DEGMAX+1) coefficient matrix 
C     TG1    double(NTGMAX)            array of the first coordinates 
C                                      of the target points
C     TG2    double(NTGMAX)            array of the second coordinates
C                                      of the target points
C     INTFTG double(NTGMAX)            values of the interpolated 
C                                      function
C     F      double function           Franke function at the target points
C
C     SIGMA1 double function           first component of the map
C                                      [-1,1]^2 -> T(U,V,W)
C     SIGMA2 double function           second component of the map
C                                      [-1,1]^2 -> T(U,V,W)
C     ISIGM1 double function           first component of the map
C                                      T(U,V,W) -> [-1,1]^2
C     ISIGM2 double function           second component of the map
C                                      T(U,V,W) -> [-1,1]^2
C     MAXERR double                    maximum norm of the error at
C                                      target points
C     ESTERR double                    estimated error
C     
C     ******************************************************************
C
      INTEGER DEGMAX,NPDMAX,NTGMAX,NTG1MX
      PARAMETER (DEGMAX = 60,NPDMAX = (DEGMAX + 1) * (DEGMAX + 2) / 2,
     1     NTG1MX = 100,NTGMAX = NTG1MX ** 2 - NTG1MX + 1)
      INTEGER DEG,NPD,NTG,FAMILY,FTYPE
      DOUBLE PRECISION PD1(NPDMAX),PD2(NPDMAX),WPD(NPDMAX),FPD(NPDMAX)
      DOUBLE PRECISION RAUX1((DEGMAX + 1) * (DEGMAX + 2)),
     1     RAUX2((DEGMAX + 1)* (DEGMAX + 2))
      DOUBLE PRECISION C0(0:DEGMAX + 1,0:DEGMAX + 1)
      DOUBLE PRECISION TG1(NTGMAX),TG2(NTGMAX),INTFTG(NTGMAX)
      DOUBLE PRECISION F,MAXERR,SIGMA1,SIGMA2,PD2VAL,ESTERR,
     1     MINVAL,MAXVAL,MAXDEV,MEAN,ISIGM1,ISIGM2
      DOUBLE PRECISION U1,V1,W1,U2,V2,W2
      DATA U1/0.0D0/,U2/0.0D0/,V1/1.0D0/,V2/0.0D0/,W1/0.0D0/,W2/1.0D0/
      EXTERNAL F,SIGMA1,SIGMA2,ISIGM1,ISIGM2,PD2VAL
      INTEGER I,NTG1
      FAMILY = 1
      DEG = 60
      NTG1 = 100
      FTYPE = 1
      WRITE(6,*)
      WRITE(6,*) 'Interpolation of the Franke function'
      WRITE(6,*) 'on the unit triangle T((0,0),(1,0),(0,1))'
      WRITE(6,*) 'at degree =',DEG
      IF (DEG .GT. DEGMAX) THEN
         WRITE(6,*) 'IN _DRVTRI: DEG = ',DEG,' > DEGMAX = ',DEGMAX
         WRITE(6,*) 'STOP'
         STOP
      END IF
C     
C     Build the first family of Padua points in the square [-1,1]^2
C     
      CALL PDPTS(DEG,PD1,PD2,WPD,NPD)
      OPEN(10,FILE='fpd.dat',STATUS='UNKNOWN')
      WRITE(10,*) '# Function values at ',NPD,' Padua points'
      WRITE(6,*)
      WRITE(6,3) 'Function values at ',NPD,'Padua points       '//
     1     'in fpd.dat   '
      DO 10 I = 1,NPD
C     
C     Compute the Franke function at Padua points mapped to 
C     T(U,V,W)
C     
         FPD(I) = F(FTYPE,SIGMA1(PD1(I),PD2(I),U1,U2,V1,V2,W1,W2),
     1        SIGMA2(PD1(I),PD2(I),U1,U2,V1,V2,W1,W2))
         WRITE(10,2) SIGMA1(PD1(I),PD2(I),U1,U2,V1,V2,W1,W2),
     1        SIGMA2(PD1(I),PD2(I),U1,U2,V1,V2,W1,W2),FPD(I)
 10    CONTINUE
      WRITE(10,*) '# End'
      CLOSE(10)
C     
C     Compute the matrix C0 of the coefficients in the bivariate
C     orthonormal Chebyshev basis
C     
      CALL PADUA2(DEG,DEGMAX,NPD,WPD,FPD,RAUX1,RAUX2,C0,ESTERR)
C     
C     Build the set of target points on T(U,V,W)
C     
      CALL TARGET(U1,U2,V1,V2,W1,W2,NTG1,NTGMAX,TG1,TG2,NTG)
      DO 20 I = 1,NTG
C     
C     Evaluate the interpolant at target points
C     
         INTFTG(I) = PD2VAL(DEG,DEGMAX,C0,
     1        ISIGM1(TG1(I),TG2(I),U1,U2,V1,V2,W1,W2),
     2        ISIGM2(TG1(I),TG2(I),U1,U2,V1,V2,W1,W2),RAUX1,RAUX2)
 20   CONTINUE
C     
C     Compute the error relative to the max deviation from the mean
C     
      MAXERR = 0.0D0
      MEAN = 0.0D0
      MAXVAL = F(FTYPE,TG1(1),TG2(1))
      MINVAL = MAXVAL
      OPEN(10,FILE='ftg.dat',STATUS='UNKNOWN')
      WRITE(10,*) '# Function values at ',NTG,' target points'
      WRITE(6,3) 'Function values at ',NTG, 'target points      '//
     1     'in ftg.dat   '
      OPEN(11,FILE='intftg.dat',STATUS='UNKNOWN')
      WRITE(11,*) '# Interpolated function values at ',NTG,
     1     ' target points'
      WRITE(6,*) 'Interpolated function values at target points'//
     1     ' in intftg.dat'
      DO 30 I = 1,NTG
         MAXERR = MAX(MAXERR,ABS(F(FTYPE,TG1(I),TG2(I)) - INTFTG(I)))
         MEAN = MEAN + F(FTYPE,TG1(I),TG2(I))
         MAXVAL = MAX(MAXVAL,F(FTYPE,TG1(I),TG2(I)))
         MINVAL = MIN(MINVAL,F(FTYPE,TG1(I),TG2(I)))
         WRITE(10,2) TG1(I),TG2(I),F(FTYPE,TG1(I),TG2(I))
         WRITE(11,2) TG1(I),TG2(I),INTFTG(I)
         IF (MOD(I,NTG1).EQ.0) THEN
            WRITE(10,*)
            WRITE(11,*)
         END IF   
 30   CONTINUE   
      IF (MAXVAL .EQ. MINVAL) THEN
         MAXDEV = 1.0D0
      ELSE
         MEAN = MEAN / NTG
         MAXDEV = MAX(MAXVAL - MEAN,MEAN - MINVAL)
      END IF   
      WRITE(10,*) '# End'
      CLOSE(10)
      WRITE(11,*) '# End'
      CLOSE(11)
      WRITE(6,*)
      WRITE(6,1) 'Estimated error:     ',ESTERR / MAXDEV
      WRITE(6,1) 'Actual error:        ',MAXERR / MAXDEV
      WRITE(6,1) 'Expected error:      ',0.1226D-09
      WRITE(6,*)
 1    FORMAT(1X,A21,E10.4)
 2    FORMAT(3(1X,E14.8))
 3    FORMAT(1X,A18,1X,I7,1X,A32)
      STOP
      END
      DOUBLE PRECISION FUNCTION SIGMA1(T1,T2,U1,U2,V1,V2,W1,W2)
C     
C     This functions returns the first component of the map
C     from the square [-1,1]^2 to the triangle T(U,V,W). 
C
      DOUBLE PRECISION T1,T2
      DOUBLE PRECISION U1,U2,V1,V2,W1,W2
      SIGMA1 = (V1 - U1) * (1 + T1) * (1 - T2) / 4 + 
     1     (W1 - U1) * (1 + T2) / 2 + U1
      RETURN
      END
      DOUBLE PRECISION FUNCTION ISIGM1(SIGMA1,SIGMA2,U1,U2,V1,V2,W1,W2)
C     
C     This functions returns the first component of the map
C     from the the triangle T(U,V,W) to the square [-1,1]^2.
C
      DOUBLE PRECISION SIGMA1,SIGMA2,U1,U2,V1,V2,W1,W2
      DOUBLE PRECISION RHO1,RHO2
      RHO1 = (SIGMA1 * (W2 - U2) - 
     1     SIGMA2 * (W1 - U1) + (W1 - U1) * U2 - (W2 - U2) * U1) /
     2     ((V1 - U1) * (W2 - U2) - (V2 - U2) * (W1 - U1))
      RHO2 = (SIGMA1 * (V2 - U2) - 
     1     SIGMA2 * (V1 - U1) + (V1 - U1) * U2 - (V2 - U2) * U1) /
     2     ((W1 - U1) * (V2 - U2) - (W2 - U2) * (V1 - U1))
      IF (RHO2 .EQ. 1.0D0) THEN
         ISIGM1 = 0.0D0
      ELSE
         ISIGM1 = 2 * RHO1 / (1 - RHO2) - 1
      END IF   
      RETURN
      END
      DOUBLE PRECISION FUNCTION SIGMA2(T1,T2,U1,U2,V1,V2,W1,W2)
C     
C     This functions returns the second component of the map
C     from the square [-1,1]^2 to the triangle T(U,V,W).
C
      DOUBLE PRECISION T1,T2
      DOUBLE PRECISION U1,U2,V1,V2,W1,W2
      SIGMA2 = (V2 - U2) * (1 + T1) * (1 - T2) / 4 + 
     1     (W2 - U2) * (1 + T2) / 2 + U2
      RETURN
      END
      DOUBLE PRECISION FUNCTION ISIGM2(SIGMA1,SIGMA2,U1,U2,V1,V2,W1,W2)
C     
C     This functions returns the second component of the map
C     from the the triangle T(U,V,W) to the square [-1,1]^2.
C
      DOUBLE PRECISION SIGMA1,SIGMA2,U1,U2,V1,V2,W1,W2
      DOUBLE PRECISION RHO2
      RHO2 = (SIGMA1 * (V2 - U2) - 
     1     SIGMA2 * (V1 - U1) + (V1 - U1) * U2 - (V2 - U2) * U1) /
     2     ((W1 - U1) * (V2 - U2) - (W2 - U2) * (V1 - U1))
      IF (RHO2 .EQ. 1.0D0) THEN
         ISIGM2 = 1.0D0
      ELSE
         ISIGM2 = 2 * RHO2 - 1
      END IF   
      RETURN
      END
      SUBROUTINE TARGET(U1,U2,V1,V2,W1,W2,NTG1,NTGMAX,TG1,TG2,NTG)
C
C     Target points on the triangle T(U,V,W).
C     The number of target points is NTG = NTG1 ** 2 - NGT1 + 1.
C
      INTEGER NTG1,NTGMAX,NTG
      DOUBLE PRECISION U1,U2,V1,V2,W1,W2,TG1(NTGMAX),TG2(NTGMAX)
      INTEGER I,J
      IF (NTG1 .LT. 2) THEN
         WRITE(6,*) 'IN _TARGET: NTG1 = ',NTG1,' < 2'
         WRITE(6,*) 'STOP'
         STOP
      END IF      
      IF (NTG1 ** 2 - NTG1 + 1 .GT. NTGMAX) THEN
         WRITE(6,*) 'IN _TARGET: NTG1 ** 2 - NTG1 + 1 = ',
     1        NTG1 ** 2 - NTG1 + 1,' > NTGMAX = ',NTGMAX
         WRITE(6,*) 'STOP'
         STOP
      END IF      
      NTG = 0
      DO 10 I = 1,NTG1 - 1
         DO 20 J = 1,NTG1
            NTG = NTG + 1
            TG1(NTG) = 
     1           (V1 - U1) * (I - 1.D0) / (NTG1 - 1) + 
     2           (W1 - U1) * ((J - 1.D0) / (NTG1 - 1)) * 
     3           (1 - (I - 1.D0) / (NTG1 - 1)) + U1
            TG2(NTG) = 
     1           (V2 - U2) * (I - 1.D0) / (NTG1 - 1) + 
     2           (W2 - U2) * ((J - 1.D0) / (NTG1 - 1)) * 
     3           (1 - (I - 1.D0) / (NTG1 - 1)) + U2
 20      CONTINUE
 10   CONTINUE
      I = NTG1
      J = 1
      NTG = NTG + 1
      TG1(NTG) = 
     1     (V1 - U1) * (I - 1.D0) / (NTG1 - 1) + 
     2     (W1 - U1) * ((J - 1.D0) / (NTG1 - 1)) * 
     3     (1 - (I - 1.D0) / (NTG1 - 1)) + U1
      TG2(NTG) = 
     1     (V2 - U2) * (I - 1.D0) / (NTG1 - 1) + 
     2     (W2 - U2) * ((J - 1.D0) / (NTG1 - 1)) * 
     3     (1 - (I - 1.D0) / (NTG1 - 1)) + U2
      RETURN
      END
