      SUBROUTINE TSTFN2 (K,X,Y,IFLAG, F,FX,FY,FXX,FXY,FYY)
      INTEGER K, IFLAG
      DOUBLE PRECISION X, Y, F, FX, FY, FXX, FXY, FYY
C
C***********************************************************
C
C                                            Robert J. Renka
C                                  Dept. of Computer Science
C                                       Univ. of North Texas
C                                           renka@cs.unt.edu
C                                                   10/14/98
C
C   This subroutine computes the value and, optionally, the
C first and/or second partial derivatives of one of ten 
C bivariate test functions.  The first six functions were
C chosen by Richard Franke to test interpolation software.
C (See the reference below.)  The last four functions repre-
C sent more challenging surface fitting problems.
C
C On input:
C
C       K = Integer in the range 1 to 10 which determines
C           the choice of function as follows:
C
C               K = 1 - Exponential
C               K = 2 - Cliff
C               K = 3 - Saddle
C               K = 4 - Gentle
C               K = 5 - Steep
C               K = 6 - Sphere
C               K = 7 - Trig
C               K = 8 - Synergistic Gaussian
C               K = 9 - Cloverleaf Asymmetric Peak/Valley
C               K = 10 - Cosine Peak
C
C   Note that function 6 is only defined inside a circle of
C radius 8/9 centered at (.5,.5).  Thus, if (X-.5)**2 +
C (Y-.5)**2 .GE. 64/81, the value (and partials if IFLAG=1)
C are set to 0 for this function.  Also, the first partial
C derivatives of function 10 are not defined at (.5,.5) --
C again, zeros are returned.
C
C       X,Y = Coordinates of the point at which the selected
C             function is to be evaluated.
C
C       IFLAG = Derivative option indicator:
C               IFLAG = 0 if only a function value is
C                         required.
C               IFLAG = 1 if both the function and its first
C                         partial derivatives are to be
C                         evaluated.
C               IFLAG = 2 if the function, its first partial
C                         derivatives, and its second partial
C                         derivatives are to be evaluated.
C
C Input parameters are not altered by this routine.
C
C On output:
C
C       F = Value of function K at (X,Y).
C
C       FX,FY = First partial derivatives of function K at
C               (X,Y) if IFLAG >= 1, unaltered otherwise.
C
C       FXX,FXY,FYY = Second partial derivatives of function
C                     K at (X,Y) if IFLAG >= 2, unaltered
C                     otherwise.
C
C Intrinsic functions called by TSTFN2:  COS, EXP, SIN,
C                                          SQRT, TANH
C
C Reference:  R. Franke, A Critical Comparison of Some
C               Methods for Interpolation of Scattered Data,
C               Naval Postgraduate School Technical Report,
C               NPS-53-79-003, 1979.
C
C Changed the GO TO statement into the IF ELSE IF END IF statement.
C Marco Caliari <marco.caliari@univr.it>
C
C***********************************************************
C
      DOUBLE PRECISION T1, T2, T3, T4, T5, T6
      IF (K .LT. 1  .OR.  K .GT. 10) RETURN
      IF (K .EQ. 1) THEN
C
C Exponential:
C
         F = .75*EXP(-((9.*X-2.)**2 + (9.*Y-2.)**2)/4.) +
     .        .75*EXP(-((9.*X+1.)**2)/49. - (9.*Y+1.)/10.) +
     .        .5*EXP(-((9.*X-7.)**2 + (9.*Y-3.)**2)/4.) -
     .        .2*EXP(-(9.*X-4.)**2 - (9.*Y-7.)**2)
         IF (IFLAG .LT. 1) RETURN
         T1 = EXP(-((9.*X-2.)**2 + (9.*Y-2.)**2)/4.)
         T2 = EXP(-((9.*X+1.)**2)/49. - (9.*Y+1.)/10.)
         T3 = EXP(-((9.*X-7.)**2 + (9.*Y-3.)**2)/4.)
         T4 = EXP(-(9.*X-4.)**2 - (9.*Y-7.)**2)
         FX = -3.375*(9.*X-2.)*T1 - (27./98.)*(9.*X+1.)*T2
     .        -2.25*(9.*X-7.)*T3 + 3.6*(9.*X-4.)*T4
         FY = -3.375*(9.*Y-2.)*T1 - .675*T2
     .        -2.25*(9.*Y-3.)*T3 + 3.6*(9.*Y-7.)*T4
         IF (IFLAG .LT. 2) RETURN
         FXX = 15.1875*((9.*X-2.)**2 - 2.)*T1 +
     .        60.75*((9.*X+1.)**2 - 24.5)*T2 +
     .        10.125*((9.*X-7.)**2 - 2.)*T3 -
     .        64.8*((9.*X-4.)**2 - .5)*T4
         FXY = 15.1875*(9.*X-2.)*(9.*Y-2.)*T1 +
     .        (243./980.)*(9.*X+1.)*T2 +
     .        10.125*(9.*X-7.)*(9.*Y-3.)*T3 -
     .        64.8*(9.*X-4.)*(9.*Y-7.)*T4
         FYY = 15.1875*((9.*Y-2.)**2 - 2.)*T1 +
     .        .6075*T2 +
     .        10.125*((9.*Y-3.)**2 - 2.)*T3 -
     .        64.8*((9.*Y-7.)**2 - .5)*T4
      ELSE IF (K .EQ. 2) THEN
C     
C Cliff:
C
         F = (TANH(9.0*(Y-X)) + 1.0)/9.0
         IF (IFLAG .LT. 1) RETURN
         T1 = 18.0*(Y-X)
         FX = -4.0/(EXP(T1) + 2.0 + EXP(-T1))
         FY = -FX
         IF (IFLAG .LT. 2) RETURN
         FXX = 18.0*TANH(0.5*T1)*FX
         FXY = -FXX
         FYY = FXX
      ELSE IF (K .EQ. 3) THEN  
C     
C Saddle:
C     
         F = (1.25 + COS(5.4*Y))/(6.0 + 6.0*(3.0*X-1.0)**2)
         IF (IFLAG .LT. 1) RETURN
         T1 = 5.4*Y
         T2 = 1.0 + (3.0*X-1.)**2
         FX = -(3.0*X-1.0)*(1.25 + COS(T1))/(T2**2)
         FY = -.9*SIN(T1)/T2
         IF (IFLAG .LT. 2) RETURN
         FXX = 3.0*(1.25 + COS(T1))*(3.0*T2-4.0)/(T2**3)
         FXY = 5.4*(3.0*X-1.0)*SIN(T1)/(T2**2)
         FYY = -4.86*COS(T1)/T2
      ELSE IF (K .EQ. 4) THEN  
C     
C Gentle:
C
         F = EXP(-5.0625*((X-.5)**2 + (Y-.5)**2))/3.0
         IF (IFLAG .LT. 1) RETURN
         T1 = X - .5
         T2 = Y - .5
         T3 = -3.375*EXP(-5.0625*(T1**2 + T2**2))
         FX = T1*T3
         FY = T2*T3
         IF (IFLAG .LT. 2) RETURN
         FXX = (1.0 - 10.125*T1*T1)*T3
         FXY = -10.125*T1*T2*T3
         FYY = (1.0 - 10.125*T2*T2)*T3
      ELSE IF (K .EQ. 5) THEN  
C     
C Steep:
C
         F = EXP(-20.25*((X-.5)**2 + (Y-.5)**2))/3.0
         IF (IFLAG .LT. 1) RETURN
         T1 = X - .5
         T2 = Y - .5
         T3 = -13.5*EXP(-20.25*(T1**2 + T2**2))
         FX = T1*T3
         FY = T2*T3
         IF (IFLAG .LT. 2) RETURN
         FXX = (1.0 - 40.5*T1*T1)*T3
         FXY = -40.5*T1*T2*T3
         FYY = (1.0 - 40.5*T2*T2)*T3
      ELSE IF (K .EQ. 6) THEN  
C     
C Sphere:
C     
         T4 = 64.0 - 81.0*((X-.5)**2 + (Y-.5)**2)
         F = 0.
         IF (T4 .GE. 0.) F = SQRT(T4)/9.0 - .5
         IF (IFLAG .LT. 1) RETURN
         T1 = X - .5
         T2 = Y - .5
         T3 = 0.
         IF (T4 .GT. 0.) T3 = -9.0/SQRT(T4)
         FX = T1*T3
         FY = T2*T3
         IF (IFLAG .LT. 2) RETURN
         FXX = (1.0 + FX*FX)*T3
         FXY = FX*FY
         FYY = (1.0 + FY*FY)*T3
      ELSE IF (K .EQ. 7) THEN  
C
C Trig:
C
         F = 2.0*COS(10.0*X)*SIN(10.0*Y) + SIN(10.0*X*Y)
         IF (IFLAG .LT. 1) RETURN
         T1 = 10.0*X
         T2 = 10.0*Y
         T3 = 10.0*COS(10.0*X*Y)
         FX = -20.0*SIN(T1)*SIN(T2) + T3*Y
         FY = 20.0*COS(T1)*COS(T2) + T3*X
         IF (IFLAG .LT. 2) RETURN
         T4 = 100.0*SIN(10.0*X*Y)
         FXX = -200.0*COS(T1)*SIN(T2) - T4*Y*Y
         FXY = -200.0*SIN(T1)*COS(T2) + T3 - T4*X*Y
         FYY = -200.0*COS(T1)*SIN(T2) - T4*X*X
      ELSE IF (K .EQ. 8) THEN  
C     
C   Gaussx(1,.5,.1) + Gaussy(.75,.5,.1) + Gaussx(1,.5,.1)*
C   Gaussy(.75,.5,.1), where Gaussx(a,b,c) is the Gaussian
C   function of x with amplitude a, center (mean) b, and
C   width (standard deviation) c.
C
         T1 = 5.0 - 10.0*X
         T2 = 5.0 - 10.0*Y
         T3 = EXP(-.5*T1*T1)
         T4 = EXP(-.5*T2*T2)
         F = T3 + .75*T4*(1.0+T3)
         IF (IFLAG .LT. 1) RETURN
         FX = T1*T3*(10.0 + 7.5*T4)
         FY = T2*T4*(7.5 + 7.5*T3)
         IF (IFLAG .LT. 2) RETURN
         FXX = T3*(T1*T1-1.0)*(100.0 + 75.0*T4)
         FXY = 75.0*T1*T2*T3*T4
         FYY = T4*(T2*T2-1.0)*(75.0 + 75.0*T3)
      ELSE IF (K .EQ. 9) THEN
C     
C Cloverleaf Asymmetric Hill/Valley:
C
         T1 = EXP((10.0 - 20.0*X)/3.0)
         T2 = EXP((10.0 - 20.0*Y)/3.0)
         T3 = 1.0/(1.0 + T1)
         T4 = 1.0/(1.0 + T2)
         T5 = 20.0/3.0
         F = (T5**3 * T1*T2)**2 * (T3*T4)**5 *
     .        (T1-2.0*T3)*(T2-2.0*T4)
         IF (IFLAG .LT. 1) RETURN
         T6 = (T5*T1*T2)**2 * (T5*T3*T4)**5
         FX = T6 * (T2-2.0*T4) *
     .        ((12.0*T3 - 3.0)*T3 + 2.0*T1 - 5.0)
         FY = T6 * (T1-2.0*T3) *
     .        ((12.0*T4 - 3.0)*T4 + 2.0*T2 - 5.0)
         IF (IFLAG .LT. 2) RETURN
         FXX = T5*T6 * (T2-2.0*T4) *
     .        (((-84.0*T3 + 78.0)*T3 + 23.0)*T3 + 4.0*T1-25.0)
         FXY = T5*T6 *
     .        ((12.0*T4 - 3.0)*T4 + 2.0*T2 - 5.0) *
     .        ((12.0*T3 - 3.0)*T3 + 2.0*T1 - 5.0)
         FYY = T5*T6 * (T1-2.0*T3) *
     .        (((-84.0*T4 + 78.0)*T4 + 23.0)*T4 + 4.0*T2-25.0)
      ELSE IF (K .EQ. 10) THEN  
C     
C Cosine Peak:
C
         T1 = SQRT( (80.0*X - 40.0)**2 + (90.0*Y - 45.0)**2 )
         T2 = EXP(-.04*T1)
         T3 = COS(.15*T1)
         F = T2*T3
         IF (IFLAG .LT. 1) RETURN
         T4 = SIN(.15*T1)
         FX = 0.
         FY = 0.
         IF (T1 .EQ. 0.) RETURN
         T4 = SIN(.15*T1)
         FX = -T2*(12.0*T4 + 3.2*T3)*(80.0*X - 40.0)/T1
         FY = -T2*(13.5*T4 + 3.6*T3)*(90.0*Y - 45.0)/T1
         IF (IFLAG .LT. 2) RETURN
         FXX = 0.
         FXY = 0.
         FYY = 0.
         IF (T1 .EQ. 0.) RETURN
         T5 = T2/(T1**3)
         FXX = T5*(T1*(76.8*T4 - 133.76*T3)*(80.0*X-40.0)**2 -
     .        (960.0*T4 + 256.0*T3)*(90.0*Y-45.0)**2 )
         FXY = T5*(T1*(86.4*T4 - 150.48*T3) + 1080.0*T4 +
     .        288.0*T3)*(80.0*X-40.0)*(90.0*Y-45.0)
         FYY = T5*(T1*(97.2*T4 - 169.29*T3)*(90.0*Y-45.0)**2 -
     .        (1215.0*T4 + 324.0*T3)*(80.0*X-40.0)**2)
      END IF   
      RETURN
      END
