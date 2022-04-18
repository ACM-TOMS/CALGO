
      SUBROUTINE QUAT_STEP(PIOUT,QOUT,H,AINERT,PIIN,QIN,NP)

c$$$%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c$$$%
c$$$% Computes a single step (angular momentum Pi, attitude quaternion Q) 
c$$$% in the time interval [0,h] for a free rigid body with inertia 
c$$$% momenta aInert = (A,B,C) where A < B < C.
c$$$%
c$$$% Dimension Piin(3),Qin(4),Piout(3),Qout(4),h(1),aInert(3),NP(1)
c$$$%
c$$$% Solves Pi'= Pi x Omega
c$$$%        Q' = Q hatmap(Omega)
c$$$% (Omega = (Pi(1)/A, Pi(2)/B, Pi(3)/c), Omega angular velocity, 
c$$$% Pi angular momentum)
c$$$%
c$$$% Input: h (time step), aInert (vector of inertia moments, 0<A<B<C),
c$$$% Piin angular momentum at time t=0, Qin attitude quaternion at 
c$$$% time t=0, NP= integer quadrature order. Use NP=0 for exact elliptic 
c$$$% integrals, while 1<= NP <= 10 gives Gaussian quadrature of 
c$$$% order 2*NP.
c$$$%
c$$$% Output: Piout angular momentum at time t=h, Qout attitude quaternion 
c$$$% at time t=h.
c$$$%
c$$$% OBS: Piin normalized so that ||Piin||=1 (euclidean norm)
c$$$%
c$$$%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


c     Define parameters
C     .. Scalar Arguments ..
      DOUBLE PRECISION H
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION AINERT(3),PIIN(3),PIOUT(3),QIN(4),QOUT(4)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION A1DIVA2,A2DIVA1,AAA0,AAA1,AAA2,ADUE,
     +                 AIBETA,AK,ALPHA2,ANI,ANN,AP01,AP11,AP21,AP31,APP,
     +                 APPP,C,CI,CN,CONST0,CONST1,CONST2,DN,
     +                 DNPREC,ELLPIVAL,EM10,EM11,EM20,EM21,EM30,EM31,F1,
     +                 PHINOLL,PHIOUT,PI,S,SGN,SI,SN,
     +                 SNPREC,T,THETA
     +                 ,C1,C2,D1,D2,D3,AK2,ALAM
      INTEGER I,NP
C     ..
C     .. Local Arrays ..
      DOUBLE PRECISION Q1(4),Q2(4),Q3(4)
C     ..
C     .. External Subroutines ..
      EXTERNAL EGELLIPX_STEP,ELLINT_PI,QUAT_PROD,EINT_GAUSS
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ATAN,ATAN2,COS,SIN,SQRT
C     ..
      PI = 3.141592653589793238462643383279502884197169399375105820974D0

      T = PIIN(1)**2/AINERT(1)+PIIN(2)**2/AINERT(2)+PIIN(3)**2/AINERT(3)

      C1 = AINERT(1)/AINERT(2)
     +     *(AINERT(3)-AINERT(2))/(AINERT(3)-AINERT(1))
      ALAM = SQRT((AINERT(3)-AINERT(2))*(AINERT(3)-AINERT(1))
     +     /(AINERT(1)*AINERT(2)*AINERT(3)**2))
      
      C2 = 1.0D0 - C1
      C1 = 1.0D0 - C2

      D1 = SQRT(PIIN(1)**2+C1*PIIN(2)**2)
      D3 = SQRT(C2*PIIN(2)**2+PIIN(3)**2)

c     End parameters

      IF (ABS(D3).LE.1.0D-14) THEN
c     Begin{special cases} in which the is only rotation around an axis
         DO I = 1,3
            PIOUT(I) = PIIN(I)
         END DO
         
         THETA = PIIN(1)*H/AINERT(1)
         S = SIN(THETA*0.50D0)
         C = COS(THETA*0.50D0)
         
         Q3(1) = C
         Q3(2) = S
         Q3(3) = 0.0D0
         Q3(4) = 0.0D0
         
         CALL QUAT_PROD(QOUT,QIN,Q3)
         
         RETURN

      ELSE IF (ABS(D1).LE.1.0D-14) THEN
         DO I = 1,3
            PIOUT(I) = PIIN(I)
         END DO
         
         THETA = PIIN(3)*H/AINERT(3)
         S = SIN(THETA*0.50D0)
         C = COS(THETA*0.50D0)
         
         Q3(1) = C
         Q3(2) = 0.0D0
         Q3(3) = 0.0D0
         Q3(4) = S
         
         CALL QUAT_PROD(QOUT,QIN,Q3)
         
         RETURN
         
      ELSEIF (ABS((AINERT(3)-AINERT(1))).LT.1.0D-15) THEN
c         spherical body, AINERT(3)=AINERT(2)=AINERT(1)
         DO I=1,3
            PIOUT(I) = PIIN(I)
         ENDDO
         THETA=H/AINERT(1)
         S=SIN(THETA/2.0D0)
         C=COS(THETA/2.0D0)       
         Q3(1)=C
         Q3(2)=S*PIIN(1)
         Q3(3)=S*PIIN(2)
         Q3(4)=S*PIIN(3)
         CALL QUAT_PROD(QOUT, QIN,Q3)
         
         RETURN
c     End{special cases}

       ELSE

         AK2 = (C2/C1)*(D1/D3)**2
         AK = SQRT(C2/C1)*(D1/D3)
         SGN = 1.0D0
         
         
c     Generic motion
         
c$$$  %*****************************************************
c$$$  % First case: k<=1 
c$$$  %*****************************************************
          IF (AK2.LE.1.0D0) THEN
             ANI = ALAM*D3
             
             D2 = SQRT(PIIN(1)**2/C1+PIIN(2)**2)
             PHINOLL = ATAN2(-PIIN(2)/D2, -PIIN(1)/D1)+PI
             IF (PIIN(3).LT.0.0D0) THEN
                ANI = -ANI
                D3=-D3
                SGN = -SGN
             END IF
                        
            CALL EGELLIPX_STEP(SN,CN,DN,PHIOUT,H,PHINOLL,AK2,ANI)
            
            PIOUT(1) = D1*CN
            PIOUT(2) = D2*SN
            PIOUT(3) =SGN*SQRT(C2*(PIIN(2)**2-PIOUT(2)**2)+PIIN(3)**2)
            
            EM10 = PIIN(1)
            EM20 = PIIN(2)
            EM30 = PIIN(3)
            
            EM11 = PIOUT(1)
            EM21 = PIOUT(2)
            EM31 = PIOUT(3)
            
           SNPREC = PIIN(2)/D2
           DNPREC = PIIN(3)/D3

            ALPHA2=D1**2

            AAA0 = 1.0D0 - ALPHA2
            AAA2 = SQRT(AK2*AAA0+ALPHA2)
            AAA1 = SQRT(AAA0)
            A1DIVA2 = AAA1/AAA2
            A2DIVA1 = AAA2/AAA1
            
            APPP = 1.0D0 + EM10
            APP = SQRT(APPP)
            ADUE=SQRT(2.0D0)
            Q2(1) = APP/ADUE
            Q2(2) = 0.0D0
            Q2(3) = -Q2(1)*EM30/APPP
            Q2(4) = Q2(1)*EM20/APPP
            
            APPP = 1.0D0 + EM11
            APP = SQRT(APPP)
            ADUE=SQRT(2.0D0)
            AP01 = APP/ADUE
            AP21 = AP01*EM31/APPP
            AP31 = -AP01*EM21/APPP

            ANN = ALPHA2/ (ALPHA2-1.0D0)

c     Compute elliptic integral of the third kind
c     NP=0 uses exact solution 
c     1<= NP <= 10 uses Gaussian quadrature of order 2*NP
            
            IF (NP.EQ.0) THEN 
               CALL ELLINT_PI(ELLPIVAL,PHINOLL,PHIOUT,AK,ANN)
            ELSE
               CALL EINT_GAUSS(ELLPIVAL,PHINOLL,PHIOUT,AK,ANN,NP)
            ENDIF
            
            F1 = A1DIVA2* (ATAN(A2DIVA1*SN/DN)-
     +           ATAN(A2DIVA1*SNPREC/DNPREC))
            
            CONST0 = 0.5D0/AINERT(1)
            CONST1 = (T/2.0D0-CONST0)/ (ANI*AAA0)
            CONST2 = D1


            AIBETA = H*CONST0 + CONST1* (ELLPIVAL-CONST2*F1)
            
            CI = COS(AIBETA)
            SI = SIN(AIBETA)

C     Two simplified quaternion products
            
            Q1(1) = CI*AP01
            Q1(2) = SI*AP01
            Q1(3) = CI*AP21 - SI*AP31
            Q1(4) = CI*AP31 + SI*AP21

            Q3(1) = Q2(1)*Q1(1) - Q2(3)*Q1(3) - Q2(4)*Q1(4) 
            Q3(2) = Q2(1)*Q1(2) - Q2(4)*Q1(3) + Q2(3)*Q1(4)
            Q3(3) = Q2(3)*Q1(1) + Q2(4)*Q1(2) + Q2(1)*Q1(3) 
            Q3(4) = Q2(4)*Q1(1) - Q2(3)*Q1(2) + Q2(1)*Q1(4) 
            
         ELSE
c$$$  %*****************************************************
c$$$  % Second case: k>1
c$$$  %*****************************************************
             AK2 = (C1/C2)*(D3/D1)**2
             AK = SQRT(C1/C2)*(D3/D1)

             D2= SQRT(PIIN(3)**2/C2+PIIN(2)**2)
             ALAM = SQRT((AINERT(1)-AINERT(2))*(AINERT(1)-AINERT(3))
     +     /(AINERT(1)**2*AINERT(2)*AINERT(3)))

             ANI = ALAM*D1
             PHINOLL = ATAN2(-PIIN(2)/D2, -PIIN(3)/D3)+PI

             IF (PIIN(1).LT.0.0D0) THEN
                ANI = -ANI
                D1 = -D1
                SGN = -SGN
             END IF

             CALL EGELLIPX_STEP(SN,CN,DN,PHIOUT,H,PHINOLL,AK2,ANI)

            PIOUT(3) = D3*CN
            PIOUT(2) = D2*SN
            PIOUT(1) = SGN*SQRT(PIIN(1)**2+C1*(PIIN(2)**2-PIOUT(2)**2))
            
            
            EM10 = PIIN(1)
            EM20 = PIIN(2)
            EM30 = PIIN(3)
            
            EM11 = PIOUT(1)
            EM21 = PIOUT(2)
            EM31 = PIOUT(3)
            
            SNPREC = PIIN(2)/D2
            DNPREC = PIIN(1)/D1
            
            ALPHA2 = D3**2


            AAA0 = 1.0D0 - ALPHA2
            AAA2 = SQRT(AK2*AAA0+ALPHA2)
            AAA1 = SQRT(AAA0)
            A1DIVA2 = AAA1/AAA2
            A2DIVA1 = AAA2/AAA1

            APPP = 1.0D0 + EM30
            APP = SQRT(APPP)
            ADUE=SQRT(2.0D0)
            Q2(1) = APP/ADUE
            Q2(4) = -0.0D0
            Q2(2) = -Q2(1)*EM20/APPP
            Q2(3) = Q2(1)*EM10/APPP
            
            APPP = 1.0D0 + EM31
            APP = SQRT(APPP)
            ADUE = SQRT(2.0D0)
            AP01 = APP/ADUE
            AP11 = AP01*EM21/APPP
            AP21 = -AP01*EM11/APPP
            
            CONST0 = 0.5D0/AINERT(3)
            CONST1 = (T/2.0D0-CONST0)/ (ANI*AAA0)
            CONST2 = D3

            
            ANN = ALPHA2/ (ALPHA2-1.0D0)
            
c     Compute elliptic integral of the third kind
c     NP=0 uses exact solution 
c     1<= NP <= 10 uses Gaussian quadrature of order 2*NP
            
            IF (NP.EQ.0) THEN 
               CALL ELLINT_PI(ELLPIVAL,PHINOLL,PHIOUT,AK,ANN)
            ELSE
               CALL EINT_GAUSS(ELLPIVAL,PHINOLL,PHIOUT,AK,ANN,NP)
            ENDIF
            
            F1 = A1DIVA2* (ATAN(A2DIVA1*SN/DN)-
     +           ATAN(A2DIVA1*SNPREC/DNPREC))
            AIBETA = H*CONST0 + CONST1* (ELLPIVAL-CONST2*F1)
            
            CI = COS(AIBETA)
            SI = SIN(AIBETA)
            
C     Two simplified quaternion products
            
            Q1(1) = CI*AP01
            Q1(2) = CI*AP11 - SI*AP21
            Q1(3) = CI*AP21 + SI*AP11
            Q1(4) = SI*AP01
            
            Q3(1) = Q2(1)*Q1(1) - Q2(2)*Q1(2) - Q2(3)*Q1(3) 
            Q3(2) = Q2(2)*Q1(1) + Q2(1)*Q1(2) + Q2(3)*Q1(4)
            Q3(3) = Q2(3)*Q1(1) + Q2(1)*Q1(3) - Q2(2)*Q1(4)
            Q3(4) = -Q2(3)*Q1(2) + Q2(2)*Q1(3) + Q2(1)*Q1(4)
            
         END IF
         
         
         CALL QUAT_PROD(QOUT,QIN,Q3)
         
      ENDIF
      
C     Normalization avoids propagation of rounding errors      
      
      APP=QOUT(1)**2+QOUT(2)**2+QOUT(3)**2+QOUT(4)**2
      APP=SQRT(APP)
      
      QOUT(1) = QOUT(1)/APP
      QOUT(2) = QOUT(2)/APP
      QOUT(3) = QOUT(3)/APP
      QOUT(4) = QOUT(4)/APP

      RETURN

      END

C
C     --------------------------------
C

      SUBROUTINE ELLINT_PI(AOUT,PHI1,PHI2,AK,AM)

c     The Legendre elliptic integral of the 3rd kind  Pi(phi,n, k),
c     between phi1 and phi2 is computed as AOUT = EINT2 - EINT1,
c     where EINTj := Pi(phij, n, k), j=1,2. EINTj is
c     computed using the subroutine ELLPI.

C     .. Scalar Arguments ..
      DOUBLE PRECISION AK,AM,AOUT,PHI1,PHI2
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION EINT1,EINT2,HALV,PH_1X,PH_2X,PI_HALF,SGN,VAL
      INTEGER N0,N1
C     ..
C     .. External Subroutines ..
      EXTERNAL ELLPI
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC DBLE,MOD
C     ..
      PI_HALF = 1.5707963267948966192313216916397514420985846996875529D0

C     call full elliptic integral
      CALL ELLPI(HALV,PI_HALF,-AM,AK)

      EINT1 = 0.0D0
      EINT2 = 0.0D0
      SGN = 1.0D0

c     compute eint1
      IF (PHI1.LT.0.0D0) THEN
          PHI1 = -PHI1
          SGN = -SGN
      END IF

      N0 = INT(PHI1/PI_HALF)
      PH_1X = MOD(PHI1,PI_HALF)
      IF (MOD(N0,2).EQ.0) THEN
          CALL ELLPI(VAL,PH_1X,-AM,AK)
          EINT1 = DBLE(N0)*HALV + VAL

      ELSE
          CALL ELLPI(VAL,PI_HALF-PH_1X,-AM,AK)
          EINT1 = DBLE(N0+1)*HALV - VAL
      END IF

      EINT1 = SGN*EINT1

      SGN = 1.0D0
c     compute eint2
      IF (PHI2.LT.0.0D0) THEN
          PHI2 = -PHI2
          SGN = -SGN
      END IF

      N1 = INT(PHI2/PI_HALF)
      PH_2X = MOD(PHI2,PI_HALF)

      IF (MOD(N1,2).EQ.0) THEN
          CALL ELLPI(VAL,PH_2X,-AM,AK)
          EINT2 = DBLE(N1)*HALV + VAL

      ELSE
          CALL ELLPI(VAL,PI_HALF-PH_2X,-AM,AK)
          EINT2 = DBLE(N1+1)*HALV - VAL
      END IF

      EINT2 = SGN*EINT2

      AOUT = EINT2 - EINT1

      RETURN

      END

C
C     -----------------------
C

      SUBROUTINE ELLPI(ELLPIVAL,PHI,EN,AK)
C     Legendre elliptic integral of the 3rd kind Pi(phi, n, k),
C     evaluated using Carlson’s functions RJ and
C     RF . (Note that the sign convention on n is opposite that of
C     Abramowitz and Stegun.) The
C     ranges of phi and k are 0 ≤ phi ≤ pi/2, 0 ≤ k sin phi ≤ 1.



C     .. Scalar Arguments ..
      DOUBLE PRECISION AK,ELLPIVAL,EN,PHI
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION CC,ENSS,Q,RFXYZ,RJXYZP,S
C     ..
C     .. External Subroutines ..
      EXTERNAL RF,RJ
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC COS,SIN
C     ..
      S = SIN(PHI)
      ENSS = EN*S*S
      CC = (COS(PHI))**2
      Q = (1.0D0-S*AK)* (1.0D0+S*AK)

      CALL RF(RFXYZ,CC,Q,1.0D0)
      CALL RJ(RJXYZP,CC,Q,1.0D0,1.0D0+ENSS)

      ELLPIVAL = S* (RFXYZ-ENSS*RJXYZP/3.0D0)

      RETURN

      END

C     --------------------------------------

      SUBROUTINE RF(RFXYZ,X,Y,Z)
C     Computes Carlson’s elliptic integral of the first kind, 
C     Rf(x, y, z). x and y must be nonnegative,
C     and at most one can be zero. z must be positive. 
C     TINY must be at least 5 times
C     the machine underflow limit.
C     BIG must be at most 1/5th of the machine overflow limit.


C     .. Scalar Arguments ..
      DOUBLE PRECISION RFXYZ,X,Y,Z
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION ALAMB,AVE,BIG,DELX,DELY,DELZ,E2,E3,
     +                 ERRTOL,SQRTX,SQRTY,SQRTZ,TINY,XT,XT1,YT,
     +                 YT1,ZT,ZT1
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS,MAX,MIN,SQRT
C     ..
      ERRTOL = 0.0025D0
      TINY = 1.5D-38
      BIG = 3.0D37

      IF ((MIN(X,Y,Z).LT.0.0D0) .OR. (MIN(X+Y,X+Z,Y+Z).LT.TINY) .OR.
     +    (MAX(X,Y,Z).GT.BIG)) THEN
          WRITE (6,FMT=*) 'invalid arguments in rf'

      ELSE

          XT = X
          YT = Y
          ZT = Z
          AVE = (XT+YT+ZT)/3.0D0
          DELX = (AVE-XT)/AVE
          DELY = (AVE-YT)/AVE
          DELZ = (AVE-ZT)/AVE


          DO WHILE (MAX(ABS(DELX),ABS(DELY),ABS(DELZ)).GT.ERRTOL)
              SQRTX = SQRT(XT)
              SQRTY = SQRT(YT)
              SQRTZ = SQRT(ZT)
              ALAMB = SQRTX* (SQRTY+SQRTZ) + SQRTY*SQRTZ
              XT1 = XT
              YT1 = YT
              ZT1 = ZT
              XT = 0.25D0* (XT1+ALAMB)
              YT = 0.25D0* (YT1+ALAMB)
              ZT = 0.25D0* (ZT1+ALAMB)
              AVE = (XT+YT+ZT)/3.0D0
              DELX = (AVE-XT)/AVE
              DELY = (AVE-YT)/AVE
              DELZ = (AVE-ZT)/AVE
          END DO

          E2 = DELX*DELY - DELZ*DELZ
          E3 = DELX*DELY*DELZ
          RFXYZ = (1.0D0+ (E2/24.0D0-0.1D0-3.0D0*E3/44.0D0)*E2
     +         +E3/14.0D0)/SQRT(AVE)
      END IF

      RETURN

      END


C     -------------------------------------


      SUBROUTINE RJ(RJXYZP,X,Y,Z,P)
C     Computes Carlson’s elliptic integral of the third kind,
C     RJ (x, y, z, p). x, y, and z must be
C     nonnegative, and at most one can be zero. p must be nonzero.
C     If p < 0, the Cauchy principal
C     value is returned. TINY must be at least twice the cube root
C     of the machine underflow limit,
C     BIG at most one fifth the cube root of the machine overflow
C     limit.


C     .. Scalar Arguments ..
      DOUBLE PRECISION P,RJXYZP,X,Y,Z
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION A,ALAMB,ALPHA,AVE,B,BETA,BIG,C1,C2,C3,C4,C5,C6,
     +                 C7,C8,DELP,DELX,DELY,DELZ,EA,EB,EC,ED,EE,ERRTOL,
     +                 FAC,FAC1,PT,PT1,RCAB,RCX,RFXYZ,RHO,RJXYZP1,SQRTX,
     +                 SQRTY,SQRTZ,SUM,SUM1,TAU,TINY,XT,XT1,YT,YT1,ZT,
     +                 ZT1
C     ..
C     .. External Subroutines ..
      EXTERNAL RC,RF
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS,MAX,MIN,SQRT
C     ..
      ERRTOL = 0.0015D0
      TINY = 2.5D-13
      BIG = 9.0D11
      C1 = (3.0D0/14.0D0)
      C2 = (1.0D0/3.0D0)
      C3 = (3.0D0/22.0D0)
      C4 = (3.0D0/26.0D0)
      C5 = (0.75D0*C3)
      C6 = (1.5D0*C4)
      C7 = (0.5D0*C2)
      C8 = (C3+C3)


      IF ((MIN(X,Y,Z).LT.0.0D0) .OR. (MIN(X+Y,X+Z,Y+Z,
     +    ABS(P)).LT.TINY) .OR. (MAX(X,Y,Z,ABS(P)).GT.BIG)) THEN
          WRITE (6,FMT=*) 'invalid arguments in rj'

      ELSE
          SUM = 0.0D0
          FAC = 1.0D0
          IF (P.GT.0.0D0) THEN
              XT = X
              YT = Y
              ZT = Z
              PT = P
              DELX = 1.0D0
              DELY = 1.0D0
              DELZ = 1.0D0
              DELP = 1.0D0

          ELSE
              XT = MIN(X,Y,Z)
              ZT = MAX(X,Y,Z)
              YT = X + Y + Z - XT - ZT
              A = 1.0D0/ (YT-P)
              B = A* (ZT-YT)* (YT-XT)
              PT = YT + B
              RHO = XT*ZT/YT
              TAU = P*PT/YT
              CALL RC(RCX,RHO,TAU)
              DELX = 1.0D0
              DELY = 1.0D0
              DELZ = 1.0D0
              DELP = 1.0D0
          END IF

          DO WHILE (MAX(ABS(DELX),ABS(DELY),ABS(DELZ),ABS(DELP)).GT.
     +       ERRTOL)
              SQRTX = SQRT(XT)
              SQRTY = SQRT(YT)
              SQRTZ = SQRT(ZT)
              ALAMB = SQRTX* (SQRTY+SQRTZ) + SQRTY*SQRTZ
              ALPHA = (PT* (SQRTX+SQRTY+SQRTZ)+SQRTX*SQRTY*SQRTZ)**2
              BETA = PT* (PT+ALAMB)**2
              CALL RC(RCAB,ALPHA,BETA)
              SUM1 = SUM
              SUM = SUM1 + FAC*RCAB
              FAC1 = FAC
              FAC = 0.25D0*FAC1
              XT1 = XT
              YT1 = YT
              ZT1 = ZT
              PT1 = PT
              XT = 0.25D0* (XT1+ALAMB)
              YT = 0.25D0* (YT1+ALAMB)
              ZT = 0.25D0* (ZT1+ALAMB)
              PT = 0.25D0* (PT1+ALAMB)
              AVE = 0.2D0* (XT+YT+ZT+PT+PT)
              DELX = (AVE-XT)/AVE
              DELY = (AVE-YT)/AVE
              DELZ = (AVE-ZT)/AVE
              DELP = (AVE-PT)/AVE

          END DO

          EA = DELX* (DELY+DELZ) + DELY*DELZ
          EB = DELX*DELY*DELZ
          EC = DELP*DELP
          ED = EA - 3.0D0*EC
          EE = EB + 2.0D0*DELP* (EA-EC)
          RJXYZP = 3.0D0*SUM + FAC*(1.0D0+ED* (-C1+C5*ED-C6*EE)+
     +             EB* (C7+DELP* (-C8+DELP*C4))+DELP*EA* (C2-DELP*C3)-
     +             C2*DELP*EC)/ (AVE*SQRT(AVE))

          IF (P.LE.0.0D0) THEN
              CALL RF(RFXYZ,XT,YT,ZT)
              RJXYZP1 = RJXYZP
              RJXYZP = A* (B*RJXYZP1+3.0D0*(RCX-RFXYZ))
          END IF

      END IF

      RETURN

      END

C     -----------------------------

      SUBROUTINE RC(RCXY,X,Y)
C     Computes Carlson's degenerate elliptic integral, RC(x, y). 
C     x must be nonnegative and y must be nonzero. If y < 0, the Cauchy 
c     principal value is returned. TINY must be at least 5 times
C     the machine underflow limit, BIG at most one fifth the machine 
C     maximum overflow limit.



C     .. Scalar Arguments ..
      DOUBLE PRECISION RCXY,X,Y
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION ALAMB,AVE,BIG,C1,C2,C3,C4,COMP1,COMP2,ERRTOL,S,
     +                 SQRTNY,THIRD,TINY,TNBG,W,XT,XT1,YT,YT1
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS,SQRT
C     ..
      ERRTOL = 0.0012D0
      TINY = 1.69D-38
      SQRTNY = 1.3D-19
      BIG = 3.0D37
      TNBG = (TINY*BIG)
      COMP1 = (2.236D0/SQRTNY)
      COMP2 = (TNBG*TNBG/25.0D0)
      THIRD = (1.0D0/3.0D0)
      C1 = 0.3D0
      C2 = (1.0D0/7.0D0)
      C3 = 0.375D0
      C4 = (9.0D0/22.0D0)


      IF ((X.LT.0.0D0) .OR. (Y.EQ.0.0D0) .OR. ((X+ABS(Y)).LT.TINY) .OR.
     +    ((X+ABS(Y)).GT.BIG) .OR. ((Y.LT.-COMP1).AND. (X.GT.0.0D0).AND.
     +    (X.LT.COMP2))) THEN
          WRITE (6,FMT=*) 'invalid arguments in rc'
      ELSE
          IF (Y.GT.0.0D0) THEN
              XT = X
              YT = Y
              W = 1.0D0
          ELSE
              XT = X - Y
              YT = -Y
              W = SQRT(X)/SQRT(XT)
          END IF

          S = 1.0D0
          DO WHILE (ABS(S).GT.ERRTOL)
              ALAMB = 2.0D0*SQRT(XT)*SQRT(YT) + YT
              XT1 = XT
              YT1 = YT
              XT = 0.25D0* (XT1+ALAMB)
              YT = 0.25D0* (YT1+ALAMB)
              AVE = THIRD* (XT+YT+YT)
              S = (YT-AVE)/AVE
          END DO

          RCXY = W* (1.0D0+S*S* (C1+S* (C2+S* (C3+S*C4))))/SQRT(AVE)

      END IF

      RETURN
      END


      SUBROUTINE QUAT_PROD(QOUT,Q1,Q2)
C     .. Array Arguments ..
      DOUBLE PRECISION Q1(4),Q2(4),QOUT(4)
C     ..
      QOUT(1) = Q1(1)*Q2(1) - Q1(2)*Q2(2) - Q1(3)*Q2(3) - Q1(4)*Q2(4)
      QOUT(2) = Q1(2)*Q2(1) + Q1(1)*Q2(2) - Q1(4)*Q2(3) + Q1(3)*Q2(4)
      QOUT(3) = Q1(3)*Q2(1) + Q1(4)*Q2(2) + Q1(1)*Q2(3) - Q1(2)*Q2(4)
      QOUT(4) = Q1(4)*Q2(1) - Q1(3)*Q2(2) + Q1(2)*Q2(3) + Q1(1)*Q2(4)

      RETURN

      END

C
C ---------------------------------------------------------------
C
      SUBROUTINE EGELLIPX_STEP(SN,CN,DN,PHIOUT,H,PHIIN,AM,ANY)
c$$$%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c$$$% Solve the elliptic integral (of the first kind) from 0 to phi.
c$$$%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


C     .. Scalar Arguments ..
      DOUBLE PRECISION AK,ANY,CN,DN,H,PHIIN,PHIOUT,SN
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION AM,ANTE,CH,EPS,PHIJF,PI,TEMPF,TOPI,VECPH
      INTEGER I,J,N
C     ..
C     .. Local Arrays ..
      DOUBLE PRECISION A(20),B(20),C(20),PHIIT(20),PHIN(20)
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS,ASIN,ATAN,COS,COSH,DBLE,MOD,SIN,SQRT,TAN,TANH

C     ..
      EPS = 2.220446049250313D-16
      PI = 3.141592653589793238462643383279502884197169399375105820974D0
      TOPI = 2.0D0*PI

      A(1) = 1.0D0
      AK = SQRT(AM)
      C(1) = ABS(AK)
      B(1) = SQRT(1.0D0-AM)

      I = 1
      DO WHILE (ABS(C(I)).GT.EPS)
          I = I + 1
          A(I) = 0.5D0* (A(I-1)+B(I-1))
          B(I) = SQRT(A(I-1)*B(I-1))
          C(I) = 0.5D0* (A(I-1)-B(I-1))
      END DO
      N = I - 1
      PHIIT(1) = PHIIN
      DO J = 1,N
          PHIJF = MOD(PHIIT(J),PI)
          ANTE = (PHIIT(J)-PHIJF)/PI
          TEMPF = ATAN((B(J)/A(J))*TAN(PHIJF))

          IF (TEMPF.LT.0.0D0) THEN
              TEMPF = TEMPF + PI
          END IF

          PHIIT(J+1) = TEMPF + ANTE*PI + PHIIT(J)

      END DO

      PHIN(I) = DBLE(2**N)*A(I)*ANY*H + PHIIT(I)

      VECPH = PHIN(I)
      DO WHILE (I.GT.1)
          I = I - 1
          PHIN(I) = PHIN(I+1)
          PHIN(I) = 0.5D0* (ASIN(C(I+1)*SIN(MOD(PHIN(I+1),
     +              TOPI))/A(I+1))+PHIN(I+1))
      END DO
      PHIOUT = PHIN(1)


c$$$%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c$$$% Return the Jacobi elliptic funtions
c$$$%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SN = SIN(MOD(PHIOUT,TOPI))
      CN = COS(MOD(PHIOUT,TOPI))
      DN = SQRT(1.0D0-AM*SN**2)

C      % special case m = 1
      IF (ABS(AM-1.0D0).LE.EPS) THEN
          CH = COSH(VECPH)
          SN = TANH(VECPH)
          CN = 1.0D0/CH
          DN = CN
      END IF

      RETURN

      END

C     
C     --------------------------------
C
      SUBROUTINE EINT_GAUSS(ELLPIVAL,PHI0,PHI1,AK,EN,NORDER)

C
C     Gauss-Legendre quadrature nodes (A) and weights (B).
C     The elements are in columns, i.e. for NORDER=p (order 2p), 
C     the corresponding nodes and
C     weights are A(:,p), and B(:,p)
C

C     .. Scalar Arguments ..
      DOUBLE PRECISION AK,ELLPIVAL,EN,PHI0,PHI1
      INTEGER NORDER
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION AK2,DELTA_PHI
      INTEGER K,NMAX
C     ..
C     .. Local Arrays ..
      DOUBLE PRECISION A(10,10),B(5,10),F(10),S(10)
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC SIN,SQRT
C     ..

      NMAX = 10

      IF ((NORDER.GT.NMAX) .OR. (NORDER.LT.1)) THEN
          WRITE (6,FMT=*) 'Error: only values between p=1 and p=',
     +      NMAX,'supported'
          RETURN

      END IF


C     Order 2 coefficients
      A(1,1) = 0.5D0
      B(1,1) = 1.0D0

C     Order 4 coefficients
      A(1,2) = 0.5D0 - SQRT(3.0D0)/6.0D0
      A(2,2) = 0.5D0 + SQRT(3.0D0)/6.0D0
      B(1,2) = 0.5D0

C     Order 6 coefficients
      A(1,3) = 0.5D0 - sqrt(3.0D0/5.0D0)/2.0D0
      A(2,3) = 0.5D0
      A(3,3) = 0.5D0 + sqrt(3.0D0/5.0D0)/2.0D0

      B(1,3) = 5.0D0/18.0D0
      B(2,3) = 8.0D0/18.0D0

C     Order 8 coefficients
      A(1,4) =0.5D0-sqrt((15.0D0+2.0D0*sqrt(30.0D0))/35.0D0)/2.0D0
      A(2,4) =0.5D0-sqrt((15.0D0-2.0D0*sqrt(30.0D0))/35.0D0)/2.0D0
      A(3,4) =0.5D0+sqrt((15.0D0-2.0D0*sqrt(30.0D0))/35.0D0)/2.0D0
      A(4,4) =0.5D0+sqrt((15.0D0+2.0D0*sqrt(30.0D0))/35.0D0)/2.0D0

      B(1,4) = (18.0D0-sqrt(30.0D0))/72.0D0
      B(2,4) = (18.0D0+sqrt(30.0D0))/72.0D0

C     Order 10 coefficients
      A(1,5) = 0.5D0-sqrt((35.0D0+2.0D0*sqrt(70.0D0))/7.0D0)/6.0D0
      A(2,5) = 0.5D0-sqrt((35.0D0-2.0D0*sqrt(70.0D0))/7.0D0)/6.0D0
      A(3,5) = 0.5D0
      A(4,5) = 0.5D0+sqrt((35.0D0-2.0D0*sqrt(70.0D0))/7.0D0)/6.0D0
      A(5,5) = 0.5D0+sqrt((35.0D0+2.0D0*sqrt(70.0D0))/7.0D0)/6.0D0

      B(1,5) = (322.0D0-13.0D0*sqrt(70.0D0))/1800.0D0
      B(2,5) = (322.0D0+13.0D0*sqrt(70.0D0))/1800.0D0
      B(3,5) = 128.0D0/450.0D0

C     Order 12 coefficients
      A(1,6)=0.0337652428984239860938492227530026954326171311438550875D0
      A(2,6)=0.1693953067668677431693002024900473264967757178024149645D0
      A(3,6)=0.3806904069584015456847491391596440322906946849299893249D0
      A(4,6)=0.6193095930415984543152508608403559677093053150700106750D0
      A(5,6)=0.8306046932331322568306997975099526735032242821975850354D0
      A(6,6)=0.9662347571015760139061507772469973045673828688561449124D0

      B(1,6)=0.0856622461895851725201480710863664467634112507420219911D0
      B(2,6)=0.1803807865240693037849167569188580558307609463733727411D0
      B(3,6)=0.2339569672863455236949351719947754974058278028846052676D0

C     Order 14 coefficients
      A(1,7)=0.0254460438286207377369051579760743687996145311646911082D0
      A(2,7)=0.1292344072003027800680676133596057964629261764293048699D0
      A(3,7)=0.2970774243113014165466967939615192683263089929503149368D0
      A(4,7)=0.5D0
      A(5,7)=0.7029225756886985834533032060384807316736910070496850631D0
      A(6,7)=0.8707655927996972199319323866403942035370738235706951300D0
      A(7,7)=0.9745539561713792622630948420239256312003854688353088917D0

      B(1,7)=0.0647424830844348466353057163395410091642937011299733319D0
      B(2,7)=0.1398526957446383339507338857118897912434625326132993822D0
      B(3,7)=0.1909150252525594724751848877444875669391825417669313673D0
      B(4,7)=0.2089795918367346938775510204081632653061224489795918367D0

C     Order 16 coefficients
      A(1,8)=0.0198550717512318841582195657152635047858823828492739808D0
      A(2,8)=0.1016667612931866302042230317620847815814141341920175839D0
      A(3,8)=0.2372337950418355070911304754053768254790178784398035711D0
      A(4,8)=0.4082826787521750975302619288199080096666210935435131088D0
      A(5,8)=0.5917173212478249024697380711800919903333789064564868911D0
      A(6,8)=0.7627662049581644929088695245946231745209821215601964288D0
      A(7,8)=0.8983332387068133697957769682379152184185858658079824160D0
      A(8,8)=0.9801449282487681158417804342847364952141176171507260191D0

      B(1,8)=0.0506142681451881295762656771549810950576970455258424785D0
      B(2,8)=0.1111905172266872352721779972131204422150654350256247823D0
      B(3,8)=0.1568533229389436436689811009933006566301644995013674688D0
      B(4,8)=0.1813418916891809914825752246385978060970730199471652702D0


C     Order 18 coefficients
      A(1,9)=0.0159198802461869550822118985481635649752975997540373352D0
      A(2,9)=0.0819844463366821028502851059651325617279466409376620019D0
      A(3,9)=0.1933142836497048013456489803292629076071396975297176535D0
      A(4,9)=0.3378732882980955354807309926783316957140218696315134555D0
      A(5,9)=0.5D0
      A(6,9)=0.6621267117019044645192690073216683042859781303684865444D0
      A(7,9)=0.8066857163502951986543510196707370923928603024702823464D0
      A(8,9)=0.9180155536633178971497148940348674382720533590623379980D0
      A(9,9)=0.9840801197538130449177881014518364350247024002459626647D0

      B(1,9)=0.0406371941807872059859460790552618253378308603912053753D0
      B(2,9)=0.0903240803474287020292360156214564047571689108660202422D0
      B(3,9)=0.1303053482014677311593714347093164248859201022186499759D0
      B(4,9)=0.1561735385200014200343152032922218327993774306309523227D0
      B(5,9)=0.1651196775006298815822625346434870244394053917863441672D0

C     Order 20 coefficients
      A(1,10)=0.013046735741414139961017993957773973285865026653808940D0
      A(2,10)=0.067468316655507744633951655788253475736228492517334773D0
      A(3,10)=0.160295215850487796882836317442563212115352644082595266D0
      A(4,10)=0.283302302935376404600367028417107918899964081171876751D0
      A(5,10)=0.425562830509184394557586999435140007691217570289654152D0
      A(6,10)=0.574437169490815605442413000564859992308782429710345847D0
      A(7,10)=0.716697697064623595399632971582892081100035918828123248D0
      A(8,10)=0.839704784149512203117163682557436787884647355917404733D0
      A(9,10)=0.932531683344492255366048344211746524263771507482665226D0
      A(10,10)=0.98695326425858586003898200604222602671413497334619105D0

      B(1,10)=0.033335672154344068796784404946665896428932417160079072D0
      B(2,10)=0.074725674575290296572888169828848666201278319834713683D0
      B(3,10)=0.109543181257991021997767467114081596229385935261338544D0
      B(4,10)=0.134633359654998177545613460784734676429879969230441897D0
      B(5,10)=0.147762112357376435086946497325669164710523358513426800D0
C

      DELTA_PHI = PHI1 - PHI0
      AK2 = AK**2

      DO K = 1,NORDER
          S(K) = (SIN(PHI0+DELTA_PHI*A(K,NORDER)))**2
      END DO


      DO K = 1,NORDER
          F(K) = 1.0D0/ ((1.0D0-EN*S(K))*SQRT(1.0D0-AK2*S(K)))
      END DO

      ELLPIVAL = 0.0D0
      DO K = 1,NORDER/2
          ELLPIVAL = ELLPIVAL + B(K,NORDER)*(F(K)+F(NORDER-(K-1)))
      END DO
      IF (NORDER.NE.2*(NORDER/2)) THEN
         ELLPIVAL = ELLPIVAL + B(NORDER/2+1,NORDER)*F(NORDER/2+1)
      ENDIF

      ELLPIVAL = DELTA_PHI*ELLPIVAL

      RETURN

      END
