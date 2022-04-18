Module Constants
!
!	  copyright by Danilo Erricolo
!	  University of Illinois at Chicago
!	  Oct 1, 2002 
implicit none

public :: double,MAX,MAX_,TOL,j 
integer, parameter:: MAX=400,MAX_=400,exp_coeff_factor=5


! MAX  : maximum dimension of arrays used by the subroutines contained in the module Blanch
! MAX_ : maximum dimension of arrays used by the subroutines in module Blanch that refer to subroutines inside Mathieu_Zhang_Jin
! The purpose of exp_coeff factor is to force the computation of a number of expansion coeffiencients that is at least exp_coeff_factor*order/2 inside
! the subroutine Blanch_coefficients. When this subroutine is run at double precision, a good number for exp_coeff_factor is 4. In case of quadruple precision,
! a good number is 7.


! The following line is to be used with double precision
integer, parameter:: double=kind(0.0D0)

!The following line is to be used with quadruple precision
!integer, parameter:: double=kind(0.0Q0)

real(kind=double), parameter::zero=0.0_double,one=1.0_double

real(kind=double),parameter:: TOL=1.0E-40_double
complex(kind=double), parameter:: j=(0.0_double,1.0_double)


end module  