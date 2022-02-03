      PROGRAM DRVELL
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
C        This driver computes the interpolation of the Franke function
C     on the disk E((C1,C2),ALPHA,BETA) = E((0.5,0.5),0.5,0.5)  
C     at the first family of Padua points. 
C     The degree of interpolation DEG = 60 and the number of target 
C     points is NTG = NTG1 ** 2 - 2 * NTG1 + 2, NTG1 = 100.  
C     The maps from the reference square [-1,1]^2 to the current domain 
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
C                                      [-1,1]^2 -> E((C1,C2),ALPHA,BETA)
C     SIGMA2 double function           second component of the map
C                                      [-1,1]^2 -> E((C1,C2),ALPHA,BETA)
C     ISIGM1 double function           first component of the map
C                                      E((C1,C2),ALPHA,BETA) -> [-1,1]^2
C     ISIGM2 double function           second component of the map
C                                      E((C1,C2),ALPHA,BETA) -> [-1,1]^2
C     MAXERR double                    maximum norm of the error at
C                                      target points
C     ESTERR double                    estimated error
C     
C     ******************************************************************
C
      INTEGER DEGMAX,NPDMAX,NTGMAX,NTG1MX
      PARAMETER (DEGMAX = 60,NPDMAX = (DEGMAX + 1) * (DEGMAX + 2) / 2,
     1     NTG1MX = 100,NTGMAX = NTG1MX ** 2 - 2 * NTG1MX + 2)
      INTEGER DEG,NPD,NTG,FAMILY,FTYPE
      DOUBLE PRECISION PD1(NPDMAX),PD2(NPDMAX),WPD(NPDMAX),FPD(NPDMAX)
      DOUBLE PRECISION RAUX1((DEGMAX + 1) * (DEGMAX + 2)),
     1     RAUX2((DEGMAX + 1)* (DEGMAX + 2))
      DOUBLE PRECISION C0(0:DEGMAX + 1,0:DEGMAX + 1)
      DOUBLE PRECISION TG1(NTGMAX),TG2(NTGMAX),INTFTG(NTGMAX)
      DOUBLE PRECISION F,MAXERR,SIGMA1,SIGMA2,PD2VAL,ESTERR,
     1     MINVAL,MAXVAL,MAXDEV,MEAN,ISIGM1,ISIGM2
      DOUBLE PRECISION ALPHA,BETA,C1,C2
      DATA C1/0.5D0/,C2/0.5D0/,ALPHA/0.5D0/,BETA/0.5D0/
      EXTERNAL F,SIGMA1,SIGMA2,ISIGM1,ISIGM2,PD2VAL
      INTEGER I,NTG1
      FAMILY = 1
      DEG = 60
      NTG1 = 100
      FTYPE = 1
      WRITE(6,*)
      WRITE(6,*) 'Interpolation of the Franke function'
      WRITE(6,*) 'on the disk with center = (0.5,0.5) and radius = 0.5'
      WRITE(6,*) 'at degree =',DEG
      IF (DEG .GT. DEGMAX) THEN
         WRITE(6,*) 'IN _DRVELL: DEG = ',DEG,' > DEGMAX = ',DEGMAX
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
C     E((C1,C2),ALPHA,BETA)
C     
         FPD(I) = F(1,SIGMA1(PD1(I),PD2(I),C1,C2,ALPHA,BETA),
     1        SIGMA2(PD1(I),PD2(I),C1,C2,ALPHA,BETA))
         WRITE(10,2) SIGMA1(PD1(I),PD2(I),C1,C2,ALPHA,BETA),
     1        SIGMA2(PD1(I),PD2(I),C1,C2,ALPHA,BETA),FPD(I)
 10    CONTINUE
      WRITE(10,*) '# End'
      CLOSE(10)
C     
C     Compute the matrix C0 of the coefficients in the bivariate
C     orthonormal Chebyshev basis
C     
      CALL PADUA2(DEG,DEGMAX,NPD,WPD,FPD,RAUX1,RAUX2,C0,ESTERR)
C     
C     Build the set of target points on E((C1,C2),ALPHA,BETA)
C     
      CALL TARGET(C1,C2,ALPHA,BETA,NTG1,NTGMAX,TG1,TG2,NTG)
      DO 20 I = 1,NTG
C     
C     Evaluate the interpolant at target points
C     
         INTFTG(I) = PD2VAL(DEG,DEGMAX,C0,
     1        ISIGM1(TG1(I),TG2(I),C1,C2,ALPHA,BETA),
     2        ISIGM2(TG1(I),TG2(I),C1,C2,ALPHA,BETA),RAUX1,RAUX2)
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
      WRITE(6,1) 'Expected error:      ',0.1769D-09
      WRITE(6,*)
 1    FORMAT(1X,A21,E10.4)
 2    FORMAT(3(1X,E14.8))
 3    FORMAT(1X,A18,1X,I7,1X,A32)
      STOP
      END
      DOUBLE PRECISION FUNCTION SIGMA1(T1,T2,C1,C2,ALPHA,BETA)
C
C     This function returns the first component of the map 
C     from the square [-1,1]^2 to the ellipse E((C1,C2),ALPHA,BETA).
C     
      DOUBLE PRECISION T1,T2,C1,C2,ALPHA,BETA
      DOUBLE PRECISION PHI
      EXTERNAL PHI
      SIGMA1 = C1 - ALPHA * T2 * SIN(PHI(T1))
      RETURN
      END
      DOUBLE PRECISION FUNCTION ISIGM1(SIGMA1,SIGMA2,C1,C2,ALPHA,BETA)
C
C     This function returns the first component of the map 
C     from the ellipse E((C1,C2),ALPHA,BETA) to the square [-1,1]^2. 
C 
      DOUBLE PRECISION SIGMA1,SIGMA2,C1,C2,ALPHA,BETA
      DOUBLE PRECISION IPHI
      EXTERNAL IPHI
      IF (SIGMA2 .EQ. C2) THEN
            ISIGM1 = 1.0D0
      ELSE
         ISIGM1 = IPHI(ATAN(BETA * (C1 - SIGMA1) / 
     1        (ALPHA * (SIGMA2 - C2))))
      END IF   
      RETURN
      END
      DOUBLE PRECISION FUNCTION SIGMA2(T1,T2,C1,C2,ALPHA,BETA)
C
C     This function returns the second component of the map 
C     from the square [-1,1]^2 to the ellipse E((C1,C2),ALPHA,BETA).
C 
      DOUBLE PRECISION T1,T2,C1,C2,ALPHA,BETA
      DOUBLE PRECISION PHI
      EXTERNAL PHI
      SIGMA2 = C2 + BETA * T2 * COS(PHI(T1))
      RETURN
      END
      DOUBLE PRECISION FUNCTION ISIGM2(SIGMA1,SIGMA2,C1,C2,ALPHA,BETA)
C
C     This function returns the second component of the map 
C     from the ellipse E((C1,C2),ALPHA,BETA) to the square [-1,1]^2. 
C 
      DOUBLE PRECISION SIGMA1,SIGMA2,C1,C2,ALPHA,BETA
      IF (SIGMA2 .EQ. C2) THEN
            ISIGM2 = (C1 - SIGMA1) / ALPHA
      ELSE
         ISIGM2 = SQRT(BETA ** 2 * (C1 - SIGMA1) ** 2 + 
     1        ALPHA ** 2 * (C2 - SIGMA2) ** 2) / (ALPHA * BETA) *
     2        SIGN(1.0D0,SIGMA2 - C2)
      END IF   
      RETURN
      END
      DOUBLE PRECISION FUNCTION PHI(X)
C     
C     This function returns the map 
C     from the segment [-1,1] to the segment [-pi/2,pi/2].
C
      DOUBLE PRECISION X
      DOUBLE PRECISION PI
      PARAMETER (PI = 3.1415926535897931D0)
      PHI = PI * X / 2
      RETURN
      END
      DOUBLE PRECISION FUNCTION IPHI(X)
C     
C     This function returns the map 
C     from the segment [-pi/2,pi/2] to the segment [-1,1].
C
      DOUBLE PRECISION X
      DOUBLE PRECISION PI
      PARAMETER (PI = 3.1415926535897931D0)
      IPHI = 2 * X / PI
      RETURN
      END
      SUBROUTINE TARGET(C1,C2,ALPHA,BETA,NTG1,NTGMAX,TG1,TG2,NTG)
C
C     Target points on the ellipse E((C1,C2),ALPHA,BETA).
C     The number of target points is NTG = NTG1 ** 2 - 2 * NTG1 + 2.
C
      INTEGER NTG1,NTGMAX,NTG
      DOUBLE PRECISION C1,C2,ALPHA,BETA,TG1(NTGMAX),TG2(NTGMAX)
      INTEGER I,J
      IF (NTG1 .LT. 2) THEN
         WRITE(6,*) 'IN _TARGET: NTG1 = ',NTG1,' < 2'
         WRITE(6,*) 'STOP'
         STOP
      END IF      
      IF (NTG1 ** 2 - 2 * NTG1 + 2 .GT. NTGMAX) THEN
         WRITE(6,*) 'IN _TARGET: NTG1 ** 2 - 2 * NTG1 + 2 = ',
     1        NTG1 ** 2 - 2 * NTG1 + 2,' > NTGMAX = ',NTGMAX
         WRITE(6,*) 'STOP'
         STOP
      END IF      
      I = 1
      J = 1
      NTG = 0
      NTG = NTG + 1
      TG1(NTG) = ALPHA * 
     1     (-1.0D0 + (I - 1) * 2.0D0 / (NTG1 - 1)) + C1
      TG2(NTG) =  BETA * 
     1     (-1.0D0 + (J - 1) * 2.0D0 / (NTG1 - 1)) * 
     2     SQRT(1 - (-1.0D0 + (I - 1) * 2.0D0 / (NTG1 - 1)) ** 2) 
     3     + C2
      DO 10 I = 2,NTG1-1
         DO 20 J = 1,NTG1
            NTG = NTG + 1
            TG1(NTG) = ALPHA * 
     1           (-1.0D0 + (I - 1) * 2.0D0 / (NTG1 - 1)) + C1
            TG2(NTG) =  BETA * 
     1           (-1.0D0 + (J - 1) * 2.0D0 / (NTG1 - 1)) * 
     2           SQRT(1 - (-1.0D0 + (I - 1) * 2.0D0 / (NTG1 - 1)) ** 2) 
     3           + C2
 20      CONTINUE
 10   CONTINUE
      I = NTG1
      J = 1
      NTG = NTG + 1
      TG1(NTG) = ALPHA * 
     1     (-1.0D0 + (I - 1) * 2.0D0 / (NTG1 - 1)) + C1
      TG2(NTG) =  BETA * 
     1     (-1.0D0 + (J - 1) * 2.0D0 / (NTG1 - 1)) * 
     2     SQRT(1 - (-1.0D0 + (I - 1) * 2.0D0 / (NTG1 - 1)) ** 2) 
     3     + C2
      RETURN
      END
