Module Mathieu_Zhang_Jin
!
! by Danilo Erricolo
! University of Illinois at Chicago
! July 19, 2002 
!
!
! This module contains some of the functions originally developed and copyrighted by
! Shanjie Zhang and Jianming Jin and described in
! COMPUTATION OF SPECIAL FUNCTIONS
! JOHN WILEY & SONS, 1996
! ISBN 0-471-11963-6 
! The subroutines used in this module were modified with permission from the authors by Danilo Erricolo and
!     downloaded from the URL:
!	  http://iris-lee3.ece.uiuc.edu/~jjin/jin_home.html
!     where it is written that "..All the programs and subroutines contained
!     in this archive are copyrighted. However, we give permission to the user
!     who downloads these routines to incorporate any of these routines into his
!     or her programs provided that the copyright is acknowledged. "
!
!
! The modifications consist of:
! 1) Removal of most "GOTO" statements
! 2) Introduction of the statement "implicit none"
! 3) Explicit declaration of the type of all variables
! 4) Use of the intrinsic function "KIND" to improve the portability on different platforms.

use Constants
implicit none
contains



SUBROUTINE CVA2(KD,M,Q,A)
!======================================================
!       Purpose: Calculate a specific characteristic value of
!                Mathieu functions
!       Input :  m  --- Order of Mathieu functions
!                q  --- Parameter of Mathieu functions
!                KD --- Case code
!                       KD=1 for cem(x,q)  ( m = 0,2,4,...)
!                       KD=2 for cem(x,q)  ( m = 1,3,5,...)
!                       KD=3 for sem(x,q)  ( m = 1,3,5,...)
!                       KD=4 for sem(x,q)  ( m = 2,4,6,...)
!       Output:  A  --- Characteristic value
!       Routines called:
!             (1) REFINE for finding accurate characteristic
!                 value using an iteration method
!             (2) CV0 for finding initial characteristic
!                 values using polynomial approximation
!             (3) CVQM for computing initial characteristic
!                 values for q <= 3*m
!             (3) CVQL for computing initial characteristic
!                 values for q >= m*m
!======================================================
!Modified on July 19, 2002 by
!Danilo Erricolo, University of Illinois at Chicago
!Replaced all implicit with "implicit none"
!Added the kind of the real variables to make it more portable     
!=================================================================
implicit none

	integer :: I,IFLAG,KD,M,NDIV,NN
	real(kind=double):: a,a1,a2,delta,q,q1,q2,qq

        IF (M.LE.12.OR.Q.LE.3.0*M.OR.Q.GT.M*M) THEN
            CALL CV0(KD,M,Q,A)
            IF (Q.NE.0.0D0) THEN
				Iflag=1
				CALL REFINE(KD,M,Q,A,IFLAG)
			END IF
        ELSE
           NDIV=10
           DELTA=(M-3.0)*M/NDIV
           IF ((Q-3.0*M).LE.(M*M-Q)) THEN
5             NN=INT((Q-3.0*M)/DELTA)+1
              DELTA=(Q-3.0*M)/NN
              Q1=2.0*M
              CALL CVQM(M,Q1,A1)
              Q2=3.0*M
              CALL CVQM(M,Q2,A2)
              QQ=3.0*M
              DO 10 I=1,NN
                 QQ=QQ+DELTA
                 A=(A1*Q2-A2*Q1+(A2-A1)*QQ)/(Q2-Q1)
                 IFLAG=1
                 IF (I.EQ.NN) IFLAG=-1
                 CALL REFINE(KD,M,QQ,A,IFLAG)
                 Q1=Q2
                 Q2=QQ
                 A1=A2
                 A2=A
10            CONTINUE
              IF (IFLAG.EQ.-10) THEN
                 NDIV=NDIV*2
                 DELTA=(M-3.0)*M/NDIV
                 GO TO 5
              ENDIF
           ELSE
15            NN=INT((M*M-Q)/DELTA)+1
              DELTA=(M*M-Q)/NN
              Q1=M*(M-1.0)
              CALL CVQL(KD,M,Q1,A1)
              Q2=M*M
              CALL CVQL(KD,M,Q2,A2)
              QQ=M*M
              DO 20 I=1,NN
                 QQ=QQ-DELTA
                 A=(A1*Q2-A2*Q1+(A2-A1)*QQ)/(Q2-Q1)
                 IFLAG=1
                 IF (I.EQ.NN) IFLAG=-1
                 CALL REFINE(KD,M,QQ,A,IFLAG)
                 Q1=Q2
                 Q2=QQ
                 A1=A2
                 A2=A
20            CONTINUE
              IF (IFLAG.EQ.-10) THEN
                 NDIV=NDIV*2
                 DELTA=(M-3.0)*M/NDIV
                 GO TO 15
              ENDIF
           ENDIF
        ENDIF
        RETURN
        END subroutine CVA2
!=============================================================



        SUBROUTINE REFINE(KD,M,Q,A,IFLAG) 
!
!       =====================================================
!       Purpose: calculate the accurate characteristic value
!                by the secant method
!       Input :  m --- Order of Mathieu functions
!                q --- Parameter of Mathieu functions
!                A --- Initial characteristic value
!       Output:  A --- Refineed characteristic value
!       Routine called:  CVF for computing the value of F for
!                        characteristic equation
!======================================================
!Modified on July 19, 2002 by
!Danilo Erricolo, University of Illinois at Chicago
!Replace all implicit with "implicit none"
!Added the kind of the real variables to make it more portable     
!Removed some GOTO statements and applied corrections indicated in the ERRATA
!=================================================================
        IMPLICIT none
		integer :: IFLAG,IT,KD,M,MJ
		real(kind=double)::A,EPS,F,F0,F1,Q,X,X0,X1
		

        EPS=1.0D-14
        MJ=10+M
        X0=A
        CALL CVF(KD,M,Q,X0,MJ,F0)
        X1=1.002*A
        CALL CVF(KD,M,Q,X1,MJ,F1)
        DO IT=1,100
           MJ=MJ+1
           X=X1-(X1-X0)/(1.0D0-F0/F1)
           CALL CVF(KD,M,Q,X,MJ,F)
           IF (ABS(1.0-X1/X).LT.EPS.OR.F.EQ.0.0) exit
           X0=X1
           F0=F1
           X1=X
           F1=F
		end do
        A=X  
        END subroutine REFINE
!       ======================================================



        SUBROUTINE CVF(KD,M,Q,A,MJ,F)
!
!       ======================================================
!       Purpose: Compute the value of F for characteristic
!                equation of Mathieu functions
!       Input :  m --- Order of Mathieu functions
!                q --- Parameter of Mathieu functions
!                A --- Characteristic value
!       Output:  F --- Value of F for characteristic equation
!======================================================
!Modified on July 19, 2002 by
!Danilo Erricolo, University of Illinois at Chicago
!Replace all implicit with "implicit none"
!Added the kind of the real variables to make it more portable     
!Removed some GOTO statements and applied corrections indicated in the ERRATA
!=================================================================

		implicit none
		integer :: IC,J,J0,JF,KD,L,L0,M,MJ
		real(kind=double):: A,B,F,Q,T0,T1,T2


        B=A
        IC=INT(M/2)
        L=0
        L0=0
        J0=2
        JF=IC
        IF (KD.EQ.1) L0=2
        IF (KD.EQ.1) J0=3
        IF (KD.EQ.2.OR.KD.EQ.3) L=1
        IF (KD.EQ.4) JF=IC-1
        T1=0.0D0
        DO  J=MJ,IC+1,-1
			 T1=-Q*Q/((2.0D0*J+L)**2-B+T1)
		end do
        IF (M.LE.2) THEN
           T2=0.0D0
           IF (KD.EQ.1.AND.M.EQ.0) T1=T1+T1
           IF (KD.EQ.1.AND.M.EQ.2) T1=-2.0*Q*Q/(4.0-B+T1)-4.0
           IF (KD.EQ.2.AND.M.EQ.1) T1=T1+Q
           IF (KD.EQ.3.AND.M.EQ.1) T1=T1-Q
        ELSE
           IF (KD.EQ.1) T0=4.0D0-B+2.0D0*Q*Q/B
           IF (KD.EQ.2) T0=1.0D0-B+Q
           IF (KD.EQ.3) T0=1.0D0-B-Q
           IF (KD.EQ.4) T0=4.0D0-B
           T2=-Q*Q/T0
           DO  J=J0,JF
            T2=-Q*Q/((2.0D0*J-L-L0)**2-B+T2)
           end do
        ENDIF
        F=(2.0D0*IC+L)**2+T1+T2-B

        END subroutine CVF
!=============================================================


        SUBROUTINE CV0(KD,M,Q,A0)
!
!       =====================================================
!       Purpose: Compute the initial characteristic value of
!                Mathieu functions for m <= 12  or q >= 300 or
!                q <= m*m
!       Input :  m  --- Order of Mathieu functions
!                q  --- Parameter of Mathieu functions
!       Output:  A0 --- Characteristic value
!       Routines called:
!             (1) CVQM for computing initial characteristic
!                 value for q <= 3*m
!             (2) CVQL for computing initial characteristic
!                 value for q >= m*m
!======================================================
!Modified on July 19, 2002 by
!Danilo Erricolo, University of Illinois at Chicago
!Replace all implicit with "implicit none"
!Added the kind of the real variables to make it more portable     
!=================================================================

        IMPLICIT NONE
		integer::KD, M
		real(kind=double) :: A0, Q,Q2

        Q2=Q*Q
        IF (M.EQ.0) THEN
           IF (Q.LE.1.0) THEN
              A0=(((.0036392*Q2-.0125868)*Q2+.0546875)*Q2-.5)*Q2
           ELSE IF (Q.LE.10.0) THEN
              A0=((3.999267D-3*Q-9.638957D-2)*Q-.88297)*Q+.5542818
           ELSE
              CALL CVQL(KD,M,Q,A0)
           ENDIF
        ELSE IF (M.EQ.1) THEN
           IF (Q.LE.1.0.AND.KD.EQ.2) THEN
              A0=(((-6.51E-4*Q-.015625)*Q-.125)*Q+1.0)*Q+1.0
           ELSE IF (Q.LE.1.0.AND.KD.EQ.3) THEN
              A0=(((-6.51E-4*Q+.015625)*Q-.125)*Q-1.0)*Q+1.0
           ELSE IF (Q.LE.10.0.AND. KD.EQ.2) THEN
              A0=(((-4.94603D-4*Q+1.92917D-2)*Q-.3089229)*Q+1.33372)*Q+.811752
           ELSE IF (Q.LE.10.0.AND.KD.EQ.3) THEN
              A0=((1.971096D-3*Q-5.482465D-2)*Q-1.152218)*Q+1.10427
           ELSE
              CALL CVQL(KD,M,Q,A0)
           ENDIF
        ELSE IF (M.EQ.2) THEN
           IF (Q.LE.1.0.AND.KD.EQ.1) THEN
              A0=(((-.0036391*Q2+.0125888)*Q2-.0551939)*Q2+.416667)*Q2+4.0
           ELSE IF (Q.LE.1.0.AND.KD.EQ.4) THEN
              A0=(.0003617*Q2-.0833333)*Q2+4.0
           ELSE IF (Q.LE.15.AND.KD.EQ.1) THEN
              A0=(((3.200972D-4*Q-8.667445D-3)*Q-1.829032D-4)*Q+.9919999)*Q+3.3290504
           ELSE IF (Q.LE.10.0.AND.KD.EQ.4) THEN
              A0=((2.38446D-3*Q-.08725329)*Q-4.732542D-3)*Q+4.00909
           ELSE
              CALL CVQL(KD,M,Q,A0)
           ENDIF
        ELSE IF (M.EQ.3) THEN
           IF (Q.LE.1.0.AND.KD.EQ.2) THEN
              A0=((6.348E-4*Q+.015625)*Q+.0625)*Q2+9.0
           ELSE IF (Q.LE.1.0.AND.KD.EQ.3) THEN
              A0=((6.348E-4*Q-.015625)*Q+.0625)*Q2+9.0
           ELSE IF (Q.LE.20.0.AND.KD.EQ.2) THEN
              A0=(((3.035731D-4*Q-1.453021D-2)*Q+.19069602)*Q-.1039356)*Q+8.9449274
           ELSE IF (Q.LE.15.0.AND.KD.EQ.3) THEN
              A0=((9.369364D-5*Q-.03569325)*Q+.2689874)*Q+8.771735
           ELSE
              CALL CVQL(KD,M,Q,A0)
           ENDIF
        ELSE IF (M.EQ.4) THEN
           IF (Q.LE.1.0.AND.KD.EQ.1) THEN
              A0=((-2.1E-6*Q2+5.012E-4)*Q2+.0333333)*Q2+16.0
           ELSE IF (Q.LE.1.0.AND.KD.EQ.4) THEN
              A0=((3.7E-6*Q2-3.669E-4)*Q2+.0333333)*Q2+16.0
           ELSE IF (Q.LE.25.0.AND.KD.EQ.1) THEN
              A0=(((1.076676D-4*Q-7.9684875D-3)*Q+.17344854)*Q-.5924058)*Q+16.620847
           ELSE IF (Q.LE.20.0.AND.KD.EQ.4) THEN
              A0=((-7.08719D-4*Q+3.8216144D-3)*Q+.1907493)*Q+15.744
           ELSE
              CALL CVQL(KD,M,Q,A0)
           ENDIF
        ELSE IF (M.EQ.5) THEN
           IF (Q.LE.1.0.AND.KD.EQ.2) THEN
              A0=((6.8E-6*Q+1.42E-5)*Q2+.0208333)*Q2+25.0
           ELSE IF (Q.LE.1.0.AND.KD.EQ.3) THEN
              A0=((-6.8E-6*Q+1.42E-5)*Q2+.0208333)*Q2+25.0
           ELSE IF (Q.LE.35.0.AND.KD.EQ.2) THEN
              A0=(((2.238231D-5*Q-2.983416D-3)*Q+.10706975)*Q-.600205)*Q+25.93515
           ELSE IF (Q.LE.25.0.AND.KD.EQ.3) THEN
              A0=((-7.425364D-4*Q+2.18225D-2)*Q+4.16399D-2)*Q+24.897
           ELSE
              CALL CVQL(KD,M,Q,A0)
           ENDIF
        ELSE IF (M.EQ.6) THEN
           IF (Q.LE.1.0) THEN
              A0=(.4D-6*Q2+.0142857)*Q2+36.0
           ELSE IF (Q.LE.40.0.AND.KD.EQ.1) THEN
              A0=(((-1.66846D-5*Q+4.80263D-4)*Q+2.53998D-2)*Q-.181233)*Q+36.423
           ELSE IF (Q.LE.35.0.AND.KD.EQ.4) THEN
              A0=((-4.57146D-4*Q+2.16609D-2)*Q-2.349616D-2)*Q+35.99251
           ELSE
              CALL CVQL(KD,M,Q,A0)
           ENDIF
        ELSE IF (M.EQ.7) THEN
           IF (Q.LE.10.0) THEN
              CALL CVQM(M,Q,A0)
           ELSE IF (Q.LE.50.0.AND.KD.EQ.2) THEN
              A0=(((-1.411114D-5*Q+9.730514D-4)*Q-3.097887D-3)*Q+3.533597D-2)*Q+49.0547
           ELSE IF (Q.LE.40.0.AND.KD.EQ.3) THEN
              A0=((-3.043872D-4*Q+2.05511D-2)*Q-9.16292D-2)*Q+49.19035
           ELSE
              CALL CVQL(KD,M,Q,A0)
           ENDIF
        ELSE IF (M.GE.8) THEN
           IF (Q.LE.3.*M) THEN
              CALL CVQM(M,Q,A0)
           ELSE IF (Q.GT.M*M) THEN
              CALL CVQL(KD,M,Q,A0)
           ELSE
              IF (M.EQ.8.AND.KD.EQ.1) THEN
                 A0=(((8.634308D-6*Q-2.100289D-3)*Q+.169072)*Q-4.64336)*Q+109.4211
              ELSE IF (M.EQ.8.AND.KD.EQ.4) THEN
                 A0=((-6.7842D-5*Q+2.2057D-3)*Q+.48296)*Q+56.59
              ELSE IF (M.EQ.9.AND.KD.EQ.2) THEN
                 A0=(((2.906435D-6*Q-1.019893D-3)*Q+.1101965)*Q-3.821851)*Q+127.6098
              ELSE IF (M.EQ.9.AND.KD.EQ.3) THEN
                 A0=((-9.577289D-5*Q+.01043839)*Q+.06588934)*Q+78.0198
              ELSE IF (M.EQ.10.AND.KD.EQ.1) THEN
                 A0=(((5.44927D-7*Q-3.926119D-4)*Q+.0612099)*Q-2.600805)*Q+138.1923
              ELSE IF (M.EQ.10.AND.KD.EQ.4) THEN
                 A0=((-7.660143D-5*Q+.01132506)*Q-.09746023)*Q+99.29494
              ELSE IF (M.EQ.11.AND.KD.EQ.2) THEN
                 A0=(((-5.67615D-7*Q+7.152722D-6)*Q+.01920291)*Q-1.081583)*Q+140.88
              ELSE IF (M.EQ.11.AND.KD.EQ.3) THEN
                 A0=((-6.310551D-5*Q+.0119247)*Q-.2681195)*Q+123.667
              ELSE IF (M.EQ.12.AND.KD.EQ.1) THEN
                 A0=(((-2.38351D-7*Q-2.90139D-5)*Q+.02023088)*Q-1.289)*Q+171.2723
              ELSE IF (M.EQ.12.AND.KD.EQ.4) THEN
                 A0=(((3.08902D-7*Q-1.577869D-4)*Q+.0247911)*Q-1.05454)*Q+161.471
              ENDIF
           ENDIF
        ENDIF
        END subroutine CV0
!==========================================================================================



        SUBROUTINE CVQL(KD,M,Q,A0)
!
!       ========================================================
!       Purpose: Compute the characteristic value of Mathieu
!                functions  for q >= 3m
!       Input :  m  --- Order of Mathieu functions
!                q  --- Parameter of Mathieu functions
!       Output:  A0 --- Initial characteristic value
!       ========================================================
!======================================================
!Modified on July 19, 2002 by
!Danilo Erricolo, University of Illinois at Chicago
!Replace all implicit with "implicit none"
!Added the kind of the real variables to make it more portable     
!=================================================================

        IMPLICIT NONE
		
		integer::  KD,M
		real(kind=double):: A0,C1,CV1,CV2,D1,D2,D3,D4,P1,P2,Q,W,W2,W3,W4,W6

        IF (KD.EQ.1.OR.KD.EQ.2) W=2.0D0*M+1.0D0
        IF (KD.EQ.3.OR.KD.EQ.4) W=2.0D0*M-1.0D0
        W2=W*W
        W3=W*W2
        W4=W2*W2
        W6=W2*W4
        D1=5.0+34.0/W2+9.0/W4
        D2=(33.0+410.0/W2+405.0/W4)/W
        D3=(63.0+1260.0/W2+2943.0/W4+486.0/W6)/W2
        D4=(527.0+15617.0/W2+69001.0/W4+41607.0/W6)/W3
        C1=128.0
        P2=Q/W4
        P1=SQRT(P2)
        CV1=-2.0*Q+2.0*W*SQRT(Q)-(W2+1.0)/8.0
        CV2=(W+3.0/W)+D1/(32.0*P1)+D2/(8.0*C1*P2)
        CV2=CV2+D3/(64.0*C1*P1*P2)+D4/(16.0*C1*C1*P2*P2)
        A0=CV1-CV2/(C1*P1)        
        END subroutine CVQL
!===========================================================================


        SUBROUTINE CVQM(M,Q,A0)
!
!       =====================================================
!       Purpose: Compute the characteristic value of Mathieu
!                functions for q <= m*m
!       Input :  m  --- Order of Mathieu functions
!                q  --- Parameter of Mathieu functions
!       Output:  A0 --- Initial characteristic value
!       =====================================================
!======================================================
!Modified on July 19, 2002 by
!Danilo Erricolo, University of Illinois at Chicago
!Replace all implicit with "implicit none"
!Added the kind of the real variables to make it more portable     
!=================================================================

        IMPLICIT NONE
		integer:: m
		real(kind=double) :: A0, HM1,HM3,HM5,Q

        HM1=.5*Q/(M*M-1.0)
        HM3=.25*HM1**3/(M*M-4.0)
        HM5=HM1*HM3*Q/((M*M-1.0)*(M*M-9.0))
        A0=M*M+Q*(HM1+(5.0*M*M+7.0)*HM3+(9.0*M**4+58.0*M*M+29.0)*HM5)
        RETURN
        END subroutine CVQM



!=================================================================


	SUBROUTINE JYNB(N,X,NM,BJ,DJ,BY,DY)
!
!       =====================================================
!       Purpose: Compute Bessel functions Jn(x), Yn(x) and
!                their derivatives
!       Input :  x --- Argument of Jn(x) and Yn(x) ( x ø 0 )
!                n --- Order of Jn(x) and Yn(x)
!       Output:  BJ(n) --- Jn(x)
!                DJ(n) --- Jn'(x)
!                BY(n) --- Yn(x)
!                DY(n) --- Yn'(x)
!                NM --- Highest order computed
!       Routines called:
!                MSTA1 and MSTA2 to calculate the starting 
!                point for backward recurrence
!======================================================
!Modified on July 19, 2002 by
!Danilo Erricolo, University of Illinois at Chicago
!Replaced all implicit with "implicit none"
!Added the kind of the real variables to make it more portable     
!Moved older DATA statements into the declaration part of the subroutine
!=================================================================
	IMPLICIT NONE
	
	integer :: K,N,NM,M
	real(kind=double),parameter::PI=3.141592653589793D0
	real(kind=double),parameter::R2P=.63661977236758D0
	real(kind=double):: BJ0,BJ1,BJk,BS,BY0,BY1,BYK,CU,EC,F,F1,F2,P0,P1,Q0,Q1,S0,SU,SV,T1,T2,X
	real(kind=double),dimension(0:n)::BJ,BY,DJ,DY
real(kind=double),dimension(4),parameter:: A=(/-.7031250000000000D-01,.1121520996093750D+00,-.5725014209747314D+00,.6074042001273483D+01/)
real(kind=double),dimension(4),parameter:: B=(/ .7324218750000000D-01,-.2271080017089844D+00,.1727727502584457D+01,-.2438052969955606D+02/)
real(kind=double),dimension(4),parameter:: A1=(/.1171875000000000D+00,-.1441955566406250D+00,.6765925884246826D+00,-.6883914268109947D+01/)
real(kind=double),dimension(4),parameter:: B1=(/-.1025390625000000D+00,.2775764465332031D+00,-.1993531733751297D+01,.2724882731126854D+02/)
	


	 
	   
	NM=N
	IF (X.LT.1.0D-100) THEN
	   DO  K=0,N
	      BJ(K)=0.0D0
	      DJ(K)=0.0D0
	      BY(K)=-1.0D+300
          DY(K)=1.0D+300
	   end do
	   BJ(0)=1.0D0
	   DJ(1)=0.5D0
	   RETURN
	ENDIF
	IF (X.LE.300.0.OR.N.GT.INT(0.9*X)) THEN
	   IF (N.EQ.0) NM=1
	   M=MSTA1(X,400)
	   IF (M.LT.NM) THEN
	      NM=M
	   ELSE
	      M=MSTA2(X,NM,20)
	   ENDIF
	   BS=0.0D0
	   SU=0.0D0
	   SV=0.0D0
	   F2=0.0D0
	   F1=1.0D-100
	   DO K=M,0,-1
	      F=2.0D0*(K+1.0D0)/X*F1-F2
	      IF (K.LE.NM) BJ(K)=F
	      IF (K.EQ.2*INT(K/2).AND.K.NE.0) THEN
		 BS=BS+2.0D0*F
		 SU=SU+(-1)**(K/2)*F/K
	      ELSE IF (K.GT.1) THEN
		 SV=SV+(-1)**(K/2)*K/(K*K-1.0)*F
	      ENDIF
	      F2=F1
          F1=F
       END DO
	   S0=BS+F
	   DO  K=0,NM 
			BJ(K)=BJ(K)/S0
	   end do
	   EC=LOG(X/2.0D0)+0.5772156649015329D0
	   BY0=R2P*(EC*BJ(0)-4.0D0*SU/S0)
	   BY(0)=BY0
	   BY1=R2P*((EC-1.0D0)*BJ(1)-BJ(0)/X-4.0D0*SV/S0)
	   BY(1)=BY1
	ELSE
	   T1=X-0.25D0*PI
	   P0=1.0D0
	   Q0=-0.125D0/X
	   DO K=1,4
	      P0=P0+A(K)*X**(-2*K)
          Q0=Q0+B(K)*X**(-2*K-1)
	   end do
	   CU=SQRT(R2P/X)
	   BJ0=CU*(P0*COS(T1)-Q0*SIN(T1))
	   BY0=CU*(P0*SIN(T1)+Q0*COS(T1))
	   BJ(0)=BJ0
	   BY(0)=BY0
	   T2=X-0.75D0*PI
	   P1=1.0D0
	   Q1=0.375D0/X
	   DO K=1,4
	      P1=P1+A1(K)*X**(-2*K)
          Q1=Q1+B1(K)*X**(-2*K-1)
	   end do
	   BJ1=CU*(P1*COS(T2)-Q1*SIN(T2))
	   BY1=CU*(P1*SIN(T2)+Q1*COS(T2))
	   BJ(1)=BJ1
	   BY(1)=BY1
	   DO K=2,NM
	      BJK=2.0D0*(K-1.0D0)/X*BJ1-BJ0
	      BJ(K)=BJK
	      BJ0=BJ1
       BJ1=BJK
	   end do
	ENDIF
	DJ(0)=-BJ(1)
	DO K=1,NM
        DJ(K)=BJ(K-1)-K/X*BJ(K)
    end do
	DO K=2,NM
	   BYK=2.0D0*(K-1.0D0)*BY1/X-BY0
	   BY(K)=BYK
	   BY0=BY1
       BY1=BYK
    end do
	DY(0)=-BY(1)
	DO K=1,NM
         DY(K)=BY(K-1)-K*BY(K)/X
	end do
	END subroutine JYNB
!==========================================================


!==========================================================
FUNCTION MSTA1(X,MP)
!
!       ===================================================
!       Purpose: Determine the starting point for backward  
!                recurrence such that the magnitude of    
!                Jn(x) at that point is about 10^(-MP)
!       Input :  x     --- Argument of Jn(x)
!                MP    --- Value of magnitude
!       Output:  MSTA1 --- Starting point   
!======================================================
!Modified on July 19, 2002 by
!Danilo Erricolo, University of Illinois at Chicago
!Replaced all implicit with "implicit none"
!Added the kind of the real variables to make it more portable     
!=================================================================

	IMPLICIT NONE
	integer:: IT,MP,N0,N1,NN,MSTA1
	real(kind=double):: A0,F,F0,F1,X

	A0=ABS(X)
	N0=INT(1.1*A0)+1
	F0=ENVJ(N0,A0)-MP
	N1=N0+5
	F1=ENVJ(N1,A0)-MP
	DO  IT=1,20             
	   NN=N1-(N1-N0)/(1.0D0-F0/F1)                  
	   F=ENVJ(NN,A0)-MP
	   IF(ABS(NN-N1).LT.1) exit
	   N0=N1
	   F0=F1
	   N1=NN
       F1=F
    end do
    MSTA1=NN

	END function MSTA1
!======================================================================




!========================================================================
FUNCTION MSTA2(X,N,MP)
!
!       ===================================================
!       Purpose: Determine the starting point for backward
!                recurrence such that all Jn(x) has MP
!                significant digits
!       Input :  x  --- Argument of Jn(x)
!                n  --- Order of Jn(x)
!                MP --- Significant digit
!       Output:  MSTA2 --- Starting point
!======================================================
!Modified on July 19, 2002 by
!Danilo Erricolo, University of Illinois at Chicago
!Replaced all implicit with "implicit none"
!Added the kind of the real variables to make it more portable     
!=================================================================

	IMPLICIT NONE
	integer:: IT,MP,MSTA2,N,NN,N0,N1
	real(kind=double):: A0,EJN,F,F0,F1,HMP,OBJ,X

	A0=ABS(X)
	HMP=0.5D0*MP
	EJN=ENVJ(N,A0)
	IF (EJN.LE.HMP) THEN
	   OBJ=MP
	   N0=INT(1.1*A0)
	ELSE
	   OBJ=HMP+EJN
	   N0=N
	ENDIF
	F0=ENVJ(N0,A0)-OBJ
	N1=N0+5
	F1=ENVJ(N1,A0)-OBJ
	DO IT=1,20
	   NN=N1-(N1-N0)/(1.0D0-F0/F1)
	   F=ENVJ(NN,A0)-OBJ
	   IF (ABS(NN-N1).LT.1) exit
	   N0=N1
	   F0=F1
	   N1=NN
	   F1=F
	end do
	MSTA2=NN+10
END function MSTA2

!=================================================================



!=====================================================
FUNCTION ENVJ(N,X)
!======================================================
!Modified on July 19, 2002 by
!Danilo Erricolo, University of Illinois at Chicago
!Replaced all implicit with "implicit none"
!Added the kind of the real variables to make it more portable     
!=================================================================
implicit none

integer :: N
real(kind=double):: X,ENVJ 
	
	ENVJ=0.5D0*LOG10(6.28D0*N)-N*LOG10(1.36D0*X/N)
	
END function ENVJ
!======================================================





end module Mathieu_Zhang_Jin