Program Driver
! copyright by Danilo Erricolo 
! University of Illinois at Chicago 
! 
!
! This is a driver program that shows the use of the subroutine Blanch_coefficients described in the 
! associated paper submitted to the ACM Transaction on Mathematical software entitled:
! Algorithm XXX: Expansion coefficients of Mathieu functions using Blanch's algorithm
!
!
! This Fortran 90 program requires 3 other modules
! constants.f90
! Mathieu_Zhang_Jin.f90 
! Blanch.f90
! and is usually compiled with the command:
! f90 -o Driver constants.f90 Mathieu_Zhang_Jin.f90 Blanch.f90 Driver.f90
!
! The output is a text file that contains the data of Table V described in the associated paper.

use Constants
use BLANCH


implicit none 
		call TestWronskianProperty()
end program Driver


!******************************



!*******************************************************************************************
subroutine TestWronskianProperty()


use constants
use Blanch 
implicit none

integer :: n,km1,km2,ii

real(kind=double) :: s,x
complex(kind=double):: Re_1p, Re_1,Re_2p,Re_2,Ro_1p, Ro_1,Ro_2p,Ro_2
complex(kind=double),dimension(0:100):: W1,W2


x=3.0D0
s=8.0D0 


open(unit=12,file="TestWronskianProperty.dat",status="replace",action="write")


do ii=1,49
    n=1+2*(ii-1)

!   Wronskian computed using Blanch subroutines
	Re_1p=MathieuRadial(0,1,1,n,s,x,km1)		 		
	Re_2 =MathieuRadial(0,2,0,n,s,x,km1)				
	Re_1 =MathieuRadial(0,1,0,n,s,x,km1)				
	Re_2p=MathieuRadial(0,2,1,n,s,x,km1)				

	W1(ii)=Re_1*Re_2p-Re_2*Re_1p-one
 
	Ro_1p=MathieuRadial(1,1,1,n,s,x,km2)				
	Ro_2 =MathieuRadial(1,2,0,n,s,x,km2)				
	Ro_1 =MathieuRadial(1,1,0,n,s,x,km2)				
	Ro_2p=MathieuRadial(1,2,1,n,s,x,km2)				


	W2(ii)=Ro_1*Ro_2p-Ro_2*Ro_1p-one

	print *, n

!	Prepare an output file that can be read by latex
	write(12,'(I4, A2, E20.8E3, A2, I4, A2, E20.8E3, A2, I4, A2)') n,'&',real(W1(ii)),'&',Km1,'&',real(W2(ii)),'&',Km2,'\\'
	
	
	
	
end do



close(12)



end subroutine TestWronskianProperty



