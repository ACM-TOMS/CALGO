Module Blanch
!
!	  copyright by Danilo Erricolo
!	  University of Illinois at Chicago
!	  Oct 1, 2002 
!
! This module contains 
!
! function MathieuAngular:   for the evaluation of Mathieu angular functions
! function MathieuRadial:    for the evaluation of Mathieu radial functions
! subroutine MTU12_Blanch:   to compute the Mathieu radial functions
! subroutine MTU0_Blanch:    to compute the Mathieu angular functions

! subroutine Blanch_Coefficients :to compute the expansion coefficients using Blanch's algorithm

! subroutine normalization:  to normalize the expansion coefficients

! function cm: used by the subroutine Blanch_Coefficients
! function V:  used by the subroutine Blanch_Coefficients
! function tail: used by the subroutine Blanch_Coefficients
! function Fm: used by the subroutine Blanch_Coefficients
!
!
!     The subroutines MTU0_Blanch and MTU12_Blanch were originally developed and copyrighted by
!	  Shanjie Zhang and Jianming Jin and described in
!     COMPUTATION OF SPECIAL FUNCTIONS
!	  JOHN WILEY & SONS, 1996
!	  ISBN 0-471-11963-6 
!	  The subroutines MTU0_Blanch and MTU12_Blanch were modified with permission from the authors by Danilo Erricolo and
!     downloaded from the URL
!	  http://iris-lee3.ece.uiuc.edu/~jjin/jin_home.html
!     where it is written that "..All the programs and subroutines contained
!     in this archive are copyrighted. However, we give permission to the user
!     who downloads these routines to incorporate any of these routines into his
!     or her programs provided that the copyright is acknowledged. "
! 
!
use constants
use Mathieu_Zhang_Jin
implicit none

private
public:: MathieuAngular, MathieuRadial, Blanch_Coefficients


contains

!*********************************************************************
!
     function MathieuAngular (parity,der,n,s,x,F,Norm,choice)
!     
!	  copyright by Danilo Erricolo
!	  University of Illinois at Chicago
!	  Oct 1, 2002 
!     ========================================================
!     Purpose
!     Compute the Angular Mathieu function Se, So 
!     parity: 0  Even, such as Se
!     INPUT
!             1  Odd, such as So
!     der     0 No derivatives required
!             1 Derivatives are required
!     n       order of the Mathieu function  
!     s       parameter of the Mathieu function, according to the notation of Blanch
!     x       argument of the Mathieu function
!	  OUTPUT
!     F       array of expansion coefficients 
!	  Norm    normalization coefficient
!	  choice  normalization requested: 1 (Stratton), (2) Ince, (3) Neutral
!     ========================================================
!
!	
	use constants        
	IMPLICIT NONE
        
    integer,intent(in) ::  parity, der,n
	real(kind=double),intent(in)::  s, x
	real(kind=double),intent(out),dimension(Max_+1)::F
	real(kind=double),intent(out)::  Norm
	integer,intent(in)::choice
	

	real(kind=double):: MathieuAngular,q

!   variables required by the subroutines called by this function
    integer :: kf
    real(kind=double)::  csf, csd

	
!	The routines of Zhang and Jin use the parameter q=s/4
	q=s/4

       
!	Determine case function (KF) and compute the appropriate normalization coefficient
       if (parity.eq.0) then
         KF=1
      else
       KF=2
      end if
   
       
     if ((parity==1) .and. (n==0)) then
		MathieuAngular=zero
		Norm=one
		print *, "Warning: angular odd functions are ZERO for M=0 !!!"
     else
	    
		call MTU0_Blanch(KF,n,Q,X,CSF,CSD,F,Norm,choice) 
             
	     if (der.eq.0) then
			MathieuAngular=CSF 
		 else
            MathieuAngular=CSD 
        end if
	  end if
       		 
      
      return
      end function MathieuAngular
!**************************************************************************************************


!     ==========================================================			

      function MathieuRadial(parity,mkind,der,n,s,x,km)

!     
!	  copyright by Danilo Erricolo
!	  University of Illinois at Chicago
!	  Oct 1, 2002

!     ==========================================================
!     Purpose
!     Compute the Mathieu radila funcitons Re1,Re2,Re3,Re4,Ro1,Ro2,Ro3,Ro4
!     parity: 0  Even, Re, ..
!             1  Odd,  Ro, ...
!     mkind    1,2,3,4 
!     der     0 No derivatives required
!             1 Derivatives are required
!     n       order of the Mathieu function  
!     s       parameter of the Mathieu function, according to the notation of Blanch
!     x       argument of the Mathieu function
!	  km	  number of terms actually computed in the series expansion of the Mathieu
!             radial function
!     =============================================================

	
	use constants
      IMPLICIT NONE
!     variables required by this function
      integer,intent(in) :: parity, mkind, der, n
	  integer, intent(out):: KM
	  real(kind=double), intent(in):: s,x
      real(kind=double)::  pi,q
      complex(kind=double) :: MathieuRadial
      

!     variables required by the subroutines called by this function
       integer :: kf,KC
       real(kind=double):: F1R, D1R, F2R, D2R
       

	   pi=acos(-one)	   
!	The routines of Zhang and Jin use a parameter q=s/4
	   q=s/4

       if (parity.eq.0) then
        KF=1
       else
         KF=2
       end if;

	  !Force computation of radial functions of first and second kind
	  KC=3 

      call MTU12_Blanch(KF,KC,N,Q,X,F1R,D1R,F2R,D2R,KM)

      if (der.eq.0) then
!     NO DERIVATIVES IN THE RESULT
         select case (mkind)
           case (1)
            MathieuRadial=CMPLX(f1r,zero,double)
           case (2)
            MathieuRadial=CMPLX(f2r,zero,double)
           case (3)
            MathieuRadial=CMPLX(f1r,f2r,double)
           case (4)
            MathieuRadial=CMPLX(f1r,-f2r,double)
         end select
      else
!     DERIVATIVES requested in the result
        select case (mkind)
          case (1)
            MathieuRadial=CMPLX(d1r,zero,double)
         case (2)
            MathieuRadial=CMPLX(d2r,zero,double)
         case (3)
            MathieuRadial=CMPLX(d1r,d2r,double)
         case (4)
           MathieuRadial=CMPLX(d1r,-d2r,double)
        end select
      end if

	MathieuRadial=sqrt(pi/2)*MathieuRadial ! Stratton's normalization
           
      
     end function MathieuRadial
!***************************************************************************



SUBROUTINE MTU12_Blanch(KF,KC,N,Q,X,F1R,D1R,F2R,D2R,KM)
!
!	  by Danilo Erricolo
!	  University of Illinois at Chicago
!	  Oct 1, 2002
!
!     this is a modification with permission from the authors  of the subroutine
!     MTU12
!	  originally developed and copyrighted by 
!	  Shanjie Zhang and Jianming Jin and described in
!     COMPUTATION OF SPECIAL FUNCTIONS
!	  JOHN WILEY & SONS, 1996
!	  ISBN 0-471-11963-6 
!     
!     This version was downloaded from the URL
!	  http://iris-lee3.ece.uiuc.edu/~jjin/jin_home.html
!     where it is written that "..All the programs and subroutines contained
!     in this archive are copyrighted. However, we give permission to the user
!     who downloads these routines to incorporate any of these routines into his
!     or her programs provided that the copyright is acknowledged. "
!
!
!	  The modifications consist of:
!	  1) introduced IMPLICIT NONE 
!     2) explicitly declared all variables
!	  3) introduced function KIND to allow for change of precision
!	  4) Removed GOTO statements
!     5) Replaced call to the subroutine FCOEF that computes the Mathieu coefficients 
!        with a call to the subroutine Blanch_Coefficients and made
!        related changes
!
!       ==============================================================
!       Purpose: Compute modified Mathieu functions of the first and
!                second kinds, Mcm(1)(2)(x,q) and Msm(1)(2)(x,q),
!                and their derivatives
!       Input:   KF --- Function code
!                       KF=1 for computing Mcm(x,q)
!                       KF=2 for computing Msm(x,q)
!                KC --- Function Code
!                       KC=1 for computing the first kind
!                       KC=2 for computing the second kind
!                            or Msm(2)(x,q) and Msm(2)'(x,q)
!                       KC=3 for computing both the first
!                            and second kinds
!                n  --- Order of Mathieu functions
!                q  --- Parameter of Mathieu functions ( q ø 0 ), according to the notation of McLachlan
!                x  --- Argument of Mathieu functions

!       Output:  F1R --- Mcm(1)(x,q) or Msm(1)(x,q)
!                D1R --- Derivative of Mcm(1)(x,q) or Msm(1)(x,q)
!                F2R --- Mcm(2)(x,q) or Msm(2)(x,q)
!                D2R --- Derivative of Mcm(2)(x,q) or Msm(2)(x,q)
!				 km	 --- number of terms actually computed in the series expansion of the Mathieu
!					     radial function
!				 
!
!       Routines called:
!            (1) CVA2 for computing the characteristic values
!            (2) FCOEF for computing expansion coefficients
!            (3) JYNB for computing Jn(x), Yn(x) and their
!                derivatives
!       ==============================================================

    use constants
	IMPLICIT NONE
	
	integer:: IC,K,KD,MKIND,NM
	integer,intent(out):: KM
	integer, intent(in)::KC, KF, N
	integer:: choice
	real(kind=double), intent(in):: X
	real(kind=double),parameter:: EPS=1.0D-14
	real(kind=double):: A,C1,C2,D1R,D2R,F1R,F2R,Q,U1,U2,W1,W2
	real(kind=double),dimension(0:Max_+1)::BJ1,BJ2,BY1,BY2,DJ1,DJ2,DY1,DY2
	real(kind=double),dimension(Max_+1)::FG
	real(kind=double):: Norm


	


!	Compute Mathieu eigenvalue
	IF (KF.EQ.1.AND.N.EQ.2*INT(N/2)) KD=1
	IF (KF.EQ.1.AND.N.NE.2*INT(N/2)) KD=2
	IF (KF.EQ.2.AND.N.NE.2*INT(N/2)) KD=3
	IF (KF.EQ.2.AND.N.EQ.2*INT(N/2)) KD=4

	CALL CVA2(KD,N,Q,A)
	
!   Compute expansion coefficients	
	select case(KD)
		case(1)
			mkind=0
		case(2)
			mkind=1
		case(3)
			mkind=3
		case(4)
			mkind=2
	end select
    
!	The parameter choice is not requested for radial functions. However, since the subroutine
!	Blanch_coefficients requires it, then choice is set to 1
	
	choice=1
	CALL Blanch_Coefficients(A,Q,N,mkind,FG,KM,norm,choice)

	
	IC=INT(N/2)+1
	IF (KD.EQ.4) IC=N/2
	C1=EXP(-X)
	C2=EXP(X)
	U1=SQRT(Q)*C1
	U2=SQRT(Q)*C2
	CALL JYNB(KM,U1,NM,BJ1,DJ1,BY1,DY1)
	CALL JYNB(KM,U2,NM,BJ2,DJ2,BY2,DY2)

	
!	Evaluate Mathieu radial function	
    W1=zero
	W2=zero
	IF ( (KC==1) .OR. (KC==3)) THEN
		F1R=zero
		DO  K=1,KM
			IF (KD.EQ.1) THEN
				F1R=F1R+(-1)**(IC+K)*FG(K)*BJ1(K-1)*BJ2(K-1)
			ELSE IF (KD.EQ.2.OR.KD.EQ.3) THEN
				F1R=F1R+(-1)**(IC+K)*FG(K)*(BJ1(K-1)*BJ2(K)+(-1)**KD*BJ1(K)*BJ2(K-1))
			ELSE
				F1R=F1R+(-1)**(IC+K)*FG(K)*(BJ1(K-1)*BJ2(K+1)-BJ1(K+1)*BJ2(K-1))
			ENDIF
			IF (K.GE.5.AND.ABS(F1R-W1).LT.ABS(F1R)*EPS) exit
			W1=F1R
		end do    
		F1R=F1R/FG(1)

		D1R=zero
		DO K=1,KM
			IF (KD.EQ.1) THEN
				D1R=D1R+(-1)**(IC+K)*FG(K)*(C2*BJ1(K-1)*DJ2(K-1)-C1*DJ1(K-1)*BJ2(K-1))
			ELSE IF (KD.EQ.2.OR.KD.EQ.3) THEN
				D1R=D1R+(-1)**(IC+K)*FG(K)*(C2*(BJ1(K-1)*DJ2(K)+(-1)**KD*BJ1(K)*DJ2(K-1))-C1*(DJ1(K-1)*BJ2(K)+(-1)**KD*DJ1(K)*BJ2(K-1)))
			ELSE
				D1R=D1R+(-1)**(IC+K)*FG(K)*(C2*(BJ1(K-1)*DJ2(K+1)-BJ1(K+1)*DJ2(K-1))-C1*(DJ1(K-1)*BJ2(K+1)-DJ1(K+1)*BJ2(K-1)))
			ENDIF
			IF (K.GE.5.AND.ABS(D1R-W2).LT.ABS(D1R)*EPS) exit
				W2=D1R
		end do
		D1R=D1R*SQRT(Q)/FG(1)
	END IF

	IF ( (KC==2) .OR. (KC==3)) THEN
		F2R=zero
		DO K=1,KM
			IF (KD.EQ.1) THEN
				F2R=F2R+(-1)**(IC+K)*FG(K)*BJ1(K-1)*BY2(K-1)
			ELSE IF (KD.EQ.2.OR.KD.EQ.3) THEN
				F2R=F2R+(-1)**(IC+K)*FG(K)*(BJ1(K-1)*BY2(K)+(-1)**KD*BJ1(K)*BY2(K-1))
			ELSE
				F2R=F2R+(-1)**(IC+K)*FG(K)*(BJ1(K-1)*BY2(K+1)-BJ1(K+1)*BY2(K-1))
			ENDIF
			IF (K.GE.5.AND.ABS(F2R-W1).LT.ABS(F2R)*EPS) EXIT
				W1=F2R
		END DO
		F2R=F2R/FG(1)

		D2R=zero
		DO K=1,KM
			IF (KD.EQ.1) THEN
				D2R=D2R+(-1)**(IC+K)*FG(K)*(C2*BJ1(K-1)*DY2(K-1)-C1*DJ1(K-1)*BY2(K-1))
			ELSE IF (KD.EQ.2.OR.KD.EQ.3) THEN
				D2R=D2R+(-1)**(IC+K)*FG(K)*(C2*(BJ1(K-1)*DY2(K)+(-1)**KD*BJ1(K)*DY2(K-1))-C1*(DJ1(K-1)*BY2(K)+(-1)**KD*DJ1(K)*BY2(K-1)))
			ELSE
				D2R=D2R+(-1)**(IC+K)*FG(K)*(C2*(BJ1(K-1)*DY2(K+1)-BJ1(K+1)*DY2(K-1))-C1*(DJ1(K-1)*BY2(K+1)-DJ1(K+1)*BY2(K-1)))
			ENDIF
			IF (K.GE.5.AND.ABS(D2R-W2).LT.ABS(D2R)*EPS) EXIT
				W2=D2R
		END DO
		D2R=D2R*SQRT(Q)/FG(1)

	END IF
END subroutine MTU12_Blanch
!===============================================================




!===============================================================
SUBROUTINE MTU0_Blanch(KF,N,Q,XR,CSF,CSD,FG,NORM,CHOICE)
!
!
!	  by Danilo Erricolo
!	  University of Illinois at Chicago
!	  July 19, 2002
!
!     this is a modification with permission from the authors of the subroutine
!     MTU0
!	  originally developed and copyrighted by 
!	  Shanjie Zhang and Jianming Jin and described in
!     COMPUTATION OF SPECIAL FUNCTIONS
!	  JOHN WILEY & SONS, 1996
!	  ISBN 0-471-11963-6 
!
!     This version was downloaded from the URL
!	  http://iris-lee3.ece.uiuc.edu/~jjin/jin_home.html
!     where it is written that "..All the programs and subroutines contained
!     in this archive are copyrighted. However, we give permission to the user
!     who downloads these routines to incorporate any of these routines into his
!     or her programs provided that the copyright is acknowledged. "
!
!	  The modifications consist of:
!	  1) introduced IMPLICIT NONE 
!     2) explicitly declared all variable
!	  3) introduced function KIND to allow for change of precision
!	  4) Removed GOTO statements
!     5) Replaced call to the subroutine FCOEF that computes the Mathieu coefficients 
!        with a call to the subroutine Blanch_Coefficients and made
!        related changes
!     6) Modified to accept only arguments in radians
!
!       ===============================================================
!       Purpose: Compute Mathieu functions cem(x,q) and sem(x,q)
!                and their derivatives ( q >= 0 )
!       Input :  KF  --- Function code
!                        KF=1 for computing cem(x,q) and cem'(x,q)
!                        KF=2 for computing sem(x,q) and sem'(x,q)
!               n   --- Order of Mathieu functions
!               q   --- Parameter of Mathieu functions
!               xr   --- Argument of Mathieu functions (IN RADIANS - MODIFIED BY DANILO ERRICOLO)
!      Output:  CSF --- cem(x,q) or sem(x,q)
!               CSD --- cem'x,q) or sem'x,q)
!               FG  --- Expansion coefficient of the Mathieu angular function
!      Routines called:
!           (1) CVA2 to compute the characteristic values
!           (2) Blanch_Coefficients to compute the expansion coefficients
!======================================================


		use constants
        IMPLICIT NONE
	
		integer,intent(in):: KF,N,choice
		real(kind=double), intent(in):: Q, XR
		real(kind=double), intent(OUT):: csf,csd,NORM

		integer :: IC,KD,K,KM,mkind
		real(kind=double),parameter::EPS=1.0D-14
		real(kind=double):: A
		real(kind=double),dimension(MAX_+1)::FG 
						
        IF (KF.EQ.1.AND.N.EQ.2*INT(N/2)) KD=1
        IF (KF.EQ.1.AND.N.NE.2*INT(N/2)) KD=2
        IF (KF.EQ.2.AND.N.NE.2*INT(N/2)) KD=3
        IF (KF.EQ.2.AND.N.EQ.2*INT(N/2)) KD=4

!		Compute Mathieu eigenvalue

        select case(KD)
		case(1)
			mkind=0
		case(2)
			mkind=1
		case(3)
			mkind=3
		case(4)
			mkind=2
		end select

		CALL CVA2(KD,N,Q,A)

!		Compute expansion coefficients
		CALL Blanch_Coefficients(A,Q,N,mkind,FG,KM,Norm,choice)

!		Evaluate Mathieu angular function		
		IC=INT(N/2)+1
		        
        CSF=zero
        DO K=1,KM
           IF (KD.EQ.1) THEN
              CSF=CSF+FG(K)*COS((2*K-2)*XR)
           ELSE IF (KD.EQ.2) THEN
              CSF=CSF+FG(K)*COS((2*K-1)*XR)
           ELSE IF (KD.EQ.3) THEN
              CSF=CSF+FG(K)*SIN((2*K-1)*XR)
           ELSE IF (KD.EQ.4) THEN
              CSF=CSF+FG(K)*SIN(2*K*XR)
           ENDIF
           IF (K.GE.IC.AND.ABS(FG(K)).LT.ABS(CSF)*EPS) exit
        end do

	    CSD=zero
        DO K=1,KM
           IF (KD.EQ.1) THEN
              CSD=CSD-(2*K-2)*FG(K)*SIN((2*K-2)*XR)
           ELSE IF (KD.EQ.2) THEN
              CSD=CSD-(2*K-1)*FG(K)*SIN((2*K-1)*XR)
           ELSE IF (KD.EQ.3) THEN
              CSD=CSD+(2*K-1)*FG(K)*COS((2*K-1)*XR)
           ELSE IF (KD.EQ.4) THEN
              CSD=CSD+2.0D0*K*FG(K)*COS(2*K*XR)
           ENDIF
           IF (K.GE.IC.AND.ABS(FG(K)).LT.ABS(CSD)*EPS) exit
        end do

        END subroutine
!=============================================================



!=======================================================
subroutine Blanch_Coefficients(A,Q,order,mkind,D_m,k_max,norm,choice)
!=======================================================
! copyright by Danilo Erricolo
! University of Illinois at Chicago
! Oct 14, 2002
! 
! Purpose 
!
! Compute the expansion coefficients of the series that represents the 
! Mathieu function. The subroutine implements the algorithm developed by
! Gertrude Blanch and described in: "Numerical Aspects of Mathieu eigenvalues",
! in Rend. Circ. Mat. Paler. Series II, vol. 15. 51--97, 1966.
!
! INPUT
! A: eigenvalue of the Mathieu equation that corresponds to a solution of order "m" and type given by "kind"
! Q: parameter 
! order: order of the Mathieu function 
! mkind: assumes four possible values
!        0 - even solution of period pi
!		 1 - even solution with period 2*pi 
!		 2 - odd solution of period pi
!        3 - odd solution of period 2*pi
! 
! choice:parameter that selects the type of nomrlaization according to
!		 1 - Morse-Stratton
!        2 - Ince
!        3 - Neutral

! OUTPUT
! D_m:   array of expansion coefficients	
! k_max: Number of cofficients contained in D_m
! norm:  normalization coefficient, only provided when "choice==2"

! This subroutine uses the integer parameter "exp_coeff_factor" defined in the module constants.f90
! Its purpose is to force the computation of a number of expansion coeffiencients that is at least exp_coeff_factor*order/2
! When this subroutine is run at double precision, a good number for exp_coeff_factor is 4. In case of quadruple precision,
! a good number is 7.



use constants
implicit none

real(kind=double), intent(in) :: A,Q
integer, intent(in):: order,mkind,choice
real(kind=double), dimension(0:MAX_),intent(out) :: D_m
integer, intent(out):: k_max 
real(kind=double),intent(out):: norm


integer :: k,m,m_max,m_add,m_start,M2S,M1
real(kind=double):: F_m,pi,Tm,Tm_p,U,FL
real(kind=double),dimension(1:MAX_) :: G_1,G_2,H_m
real(kind=double), dimension(0:MAX_):: A_m
logical:: loop, additional_terms

pi=acos(-one)

A_m(0)=zero
do m=1,Max
	G_2(m)=zero
	G_1(m)=zero 
	A_m(m)=zero
end do


!*********
! STAGE 1
!*********

FL=2.0D0**126



!determine the initial value of m and G_1,m according to the value of kind
select case(mkind)
	case(0) 
		m_start=2
		G_1(1)=V(0,A,Q)
	case(1)
		m_start=3
		G_1(1)=V(1,A,Q)-one
	case(2)
		m_start=4
		G_1(1)=V(2,A,Q)
	case(3)
		m_start=3
		G_1(1)=V(1,A,Q)+one
end select
m=m_start

if (Q>zero) then

!	Computation of M2S=2*s according to 1.91
	U=0.5*sqrt(3.0D0*Q+A) 
	M2s=2*floor(U)+m+4

!	Computation of coefficients Gm,1 according to 1.10 and 1.11 and Appendix 2
	H_m(1)=1/G_1(1)	
	m=m_start
	do while (m<M2S)
		M=M+2
		k=(m-m_start)/2+1
		G_1(k)=V(M-2,A,Q)-cm(m-4)*H_m(k-1)
		H_m(k)=1/G_1(k)	
		IF (abs(G_1(k))<one) exit
	END DO

	M1=M !set starting value of chaining coefficient

	IF ((M==M2S) .and. (abs(G_1(k))>one)) then
!	execute backtracking routine
		G_2(k)=tail(M2S,A,Q)
		Tm_p=G_2(k)-G_1(k) !previous value of Tm function
		m1=m		 
		
		do while (m-m_start>0) !statement (2.14)
			m=m-2
			k=k-1
			G_2(k)=cm(m-2)/(V(m,a,Q)-G_2(k+1))
			TM=G_2(k)-G_1(k)
			if (abs(TM)>abs(Tm_p)) then
				m1=m+2
				exit
			else
				Tm_p=Tm
				m1=m
			end if
		end do		
	end if	

!	***************
!	STAGE 2 and 3 
!	Note that we don't need to track the values of Tm since we already have available a good approximation for 
!	the eigenvalue a
!	***************

	if (m1<M2S) then
		k=(m2s-m_start)/2+1
		G_2(k)=tail(M2S,a,q)
		m=M2S-2
		k=k-1
		do while (m>=m1) !outer while
			H_m(k)=Hm2(m,a,q,G_2,m_start)
			if (abs(H_m(k))>=one) then
				G_2(k)=1/H_m(k)			
				H_m(k)=FL
			else
				loop=.true.
				do while (loop)
					if (m==m1) then
						if (abs(H_m(k))==zero) then
							print *,"Fatal error"
						else						
							G_2(k)=1/H_m(k)						
							H_m(k)=FL
						end if
						loop=.false.
					else									
						G_2(k)=FL
						m=m-2
						k=k-1
						F_m=Fm(M,A,Q,H_m,m_start)
						if (abs(F_m)>=abs(H_m(k+1))) then
							G_2(k)=H_m(k+1)/F_m
							loop=.false.
							H_m(k)=FL
						else
							H_m(k)=F_m/H_m(k+1)
							G_2(k)=FL
						end if
					end if
				end do
			end if
			m=m-2		
			k=k-1
		end do
	end if

!	Start computing coefficients according to normalization 2
	m=m1-2
	k=(m-m_start)/2+1
	A_m(k)=one
	k=k-1
	do while (k>=0)
		A_m(k)=A_m(k+1)/G_1(k+1)
		k=k-1
	end do

	m=m1
	k=(m-m_start)/2+1
	do while (m<=M2S)
		if (G_2(k)>=FL) then
			A_m(k)=A_m(k-2)/Fm(M,A,Q,H_m,m_start)
		else
			A_m(k)=G_2(k)*A_m(k-1)
		end if 
		m=m+2 
		k=k+1
	end do
	
!	Now compute the additional number of coefficients m_add necessary to obtain a truncation error
!	below some specified tolerance TOL. Hence evaluate 6.04 and 6.05

	m=m2s
	k=(m-m_start)/2+1
	m_add=ceiling(log(abs(A_m(k)/TOL))/log(1/abs(G_2(k))))+1


    additional_terms=.false.
	if (m_add>0) then
		m_max=m2s+m_add+2
		additional_terms=.true.
	end if

!   Modify Blanch's algorithm to force that the number of coefficients computed is at least as large as the order
	if (m_max<exp_coeff_factor*order) then
		m_max=exp_coeff_factor*order
		additional_terms=.true.
	end if

	if (additional_terms) then
!		Compute the additional G_{m,2} that are required
		m=m_max
		k=(m-m_start)/2+1
		G_2(k)=tail(m,A,Q)
		m=m-2
		k=k-1
		do while (m>M2S)	
			G_2(k)=cm(m-2)/(V(m,A,Q)-G_2(k+1))
			m=m-2
			k=k-1
		end do

!		Compute the additional remaining coefficients A_m according to the normalization (2)
		do while (m<=m_max)
			A_m(k)=G_2(k)*A_m(k-1)
			m=m+2
			k=k+1
		end do
	end if


	call normalization(mkind,choice,A_m,m1,m_max,D_m,norm)

	k_max=m_max/2+1


else
!	when q=0 only one coefficient is different from zero
	k=(order-m_start)/2+1
	m_max=order

!	normalize coefficients
    if (choice==1) then
		if (mkind<=1) then
			D_m(k)=one
		else
			D_m(k)=one/order
		end if
		if ((mkind==1) .and. (order==0)) then
			Norm=2*pi*D_m(0)**2
        else
			Norm=pi*D_m(k)**2
		end if
	else if (choice==2) then
		if ((mkind==0) .and. (order==0)) then
			D_m(k)=1/sqrt(2.0D0)
		else
			D_m(k)=one
		end if
	else
		D_m(k)=one
	end if
	
	k_max=k		   	
end if



end subroutine Blanch_Coefficients
!************************************************************************************************


!************************************************************************************************
subroutine normalization(mkind,choice,A_m,m1,m_max,B_m,norm)
!	  copyright by Danilo Erricolo
!	  University of Illinois at Chicago
!	  July 19, 2002

! choice==1 Stratton
! choice==2 Ince
! choice==3 Neutral
! norm returns a normalization coefficient only for Stratton's normalization

use constants
implicit none



integer, intent(in):: mkind,m1,m_max,choice
real(kind=double), dimension(0:MAX),intent(in) :: A_m
real(kind=double), dimension(0:MAX),intent(out) :: B_m
real(kind=double),intent(out):: norm

integer:: k,k1,k2,k3,m_start,p
real(kind=double):: A,pi,sum,sum1,sum2,sum_sq,sum_sq1,sum_sq2

pi=acos(-one)

select case(mkind)
	case(0) 
		m_start=2
	case(1)
		m_start=3
	case(2)
		m_start=4
		p=2
	case(3)
		m_start=3
		p=1
end select


if (choice==1) then
! Stratton's normalization
	sum1=zero
	sum2=zero
	sum_sq1=zero
	sum_sq2=zero
	k1=((m1-2)-m_start)/2+1
	k2=(m1-m_start)/2+1
	k3=(m_max-m_start)/2+1

!	compute sum of all coefficients starting from the numbers smaller in magnitude	
	if (mkind<=1) then
		do k=0,k1
			sum1=sum1+A_m(k)
		end do
		do k=k3,k2,-1
			sum2=sum2+A_m(k)
		end do
	
	else
		do k=0,k1
			sum1=sum1+(2*k+p)*A_m(k)
		end do
		do k=k3,k2,-1
			sum2=sum2+(2*k+p)*A_m(k)
		end do
	end if
	
!	Compute the sum of the squares of the coefficients starting from the numbers smaller in magnitude
	do k=0,k1		
		sum_sq1=sum_sq1+A_m(k)**2
	end do
	do k=k3,k2,-1
		sum_sq2=sum_sq2+A_m(k)**2
	end do

!	rinormalize all coefficients
	sum=sum1+sum2
	do k=0,k3
		B_m(k)=A_m(k)/sum
	end do

!	compute normalization coefficient of Stratton
	if (mkind==0) then
		Norm=pi*(sum_sq1+sum_sq2+A_m(0)**2)/sum**2
	else
		Norm=pi*(sum_sq1+sum_sq2)/sum**2
	end if
else if (choice==2) then
! Ince's normalization
!	compute sum of all coefficients starting from the numbers smaller in magnitude
	sum_sq1=zero
	k1=((m1-2)-m_start)/2+1
	do k=0,k1
		sum_sq1=sum_sq1+A_m(k)**2
	end do

	sum_sq2=zero
    k3=(m_max-m_start)/2+1
	k2=(m1-m_start)/2+1
	do k=k3,k2,-1
		sum_sq2=sum_sq2+A_m(k)**2
	end do
	
	if (mkind==0) then
		sum_sq=sqrt(sum_sq1+sum_sq2+A_m(0)**2)*A_m(0)/abs(A_m(0))
	else
		sum_sq=sqrt(sum_sq1+sum_sq2)*A_m(0)/abs(A_m(0))
	end if

!	rinormalize all coefficients
	do k=0,k3
		B_m(k)=A_m(k)/sum_sq
	end do
	Norm=one !this coefficient is not useful for Ince's normalization
else
!Neutral normalization
!Find the largest element in magnitude
    A=A_m(0)
	k1=(m_max-m_start)/2+1
	do k=1,k1
		if (abs(A_m(k))>abs(A)) then
		  A=A_m(k)
		end if
	end do
!	rinormalize all coefficients
	do k=0,k1
		B_m(k)=A_m(k)/A
	end do
	Norm=-one
end if
end subroutine normalization
!************************************************************************************************


!=========================================
function cm(m)
!=========================================
! returns cm as defined by 1.051
!	  copyright by Danilo Erricolo
!	  University of Illinois at Chicago
!	  July 19, 2002

implicit none

integer m,cm
 
 if (m==0) then
   cm=2
 else
   cm=1
 end if
end function cm
!=========================================


!=========================================
function V(m,a,q)
!=========================================
use constants
implicit none

real(kind=double):: V
real(kind=double),INTENT(IN):: A,Q
integer, intent(IN):: m

V=(a-m*m)/Q

end function V
!=========================================


!===========================================
function tail(M2S,alpha,q)
!===========================================

use constants
implicit none
real(kind=double):: tail
integer, intent(IN) ::M2S
real(kind=double),intent(IN):: alpha,q

integer:: k,mf
real(kind=double):: A,Ap,App,B,Bp,Bpp,bk,C,Q_k,Q_p


	!COMPUTATION OF THE TAIL OF THE CONTINUED FRACTION according to Method 1 of Appendix 1
	
	App= one
	Ap=  zero
	Bpp= zero
	Bp=  one

	k=1
	bk=V(M2S,alpha,Q)

	A=bk*Ap+App
	B=bk*Bp+Bpp
	
	Q_k=A/B
	Q_p=zero
	
	App=Ap
	Ap=A
	Bpp=Bp
	Bp=B

	do while (abs(Q_k-Q_p)>TOL)
		k=k+1
		mf=M2S+2*(k-1)
		bk=(alpha-(M2S+2*(k-1))**2)/Q
 
		A=bk*Ap-App
		B=bk*Bp-Bpp

		Q_p=Q_k
		Q_k=A/B
		App=Ap
		Ap=A
		Bpp=Bp
		Bp=B
	end do

!Computation of the tail of the continued fraction according to Method 3 of Appendix I
	C=zero
	k=mf-2
	do 
		C=1/(V(k,alpha,Q)-C)
		if (k==M2S) exit
		k=k-2
	end do

    tail=0.5D0*(Q_k+C)

end function tail


!=========================================
function Hm2(m,a,q,G_2,m_start)
!
! This function accepts as input the array G_2. The positions of the elements
! inside G_2 are related to the index m according to
! k=(m-m_start)/2+1

!	  copyright by Danilo Erricolo
!	  University of Illinois at Chicago
!	  July 19, 2002


!=========================================
use constants
implicit none


integer,intent(in):: m,m_start
real(kind=double),intent(in):: A,Q
real(kind=double),dimension(1:Max), intent(in):: G_2

integer:: k
real(kind=double):: Hm2,Vm
!compute Hm,2

k=((m+2)-m_start)/2+1
if (m>=m_start+2) then
	Vm=(a-m*m)/Q
	Hm2=Vm-G_2(k)
else
	Hm2=(Vm-G_2(k))/2.0D0
end if


end function Hm2
!=========================================

!=========================================
function Fm(m,a,q,H_m,m_start)
!
! This function accepts as input the array G_2. The positions of the elements
! inside G_2 are related to the index m according to
! k=(m-m_start)/2+1
!=========================================
use constants
implicit none


integer,INTENT(in):: m,m_start
real(kind=double),INTENT(in):: A,q
real(kind=double),DIMENSION(1:Max),intent(IN)::H_M

integer :: k
real(kind=double):: Fm,Vm

k=((m+2)-m_start)/2+1
if (m>=m_start+2) then
	Vm=(a-m*m)/Q
	Fm=Vm*H_m(k)-one
ELSE
	Fm=(Vm*H_m(k)-one)/2.0D0
end if

end function
!=========================================

end module  Blanch