module quadrature_inner_f

use rksuite_90_prec

implicit none

contains 

function hp(x,h)
real(kind=wp), intent(in) :: x
real(kind=wp), intent(in) :: h
real(kind=wp) :: hp

hp = (1.0_wp + x)**3

end function hp

end module quadrature_inner_f


module quadrature_f

! PROBLEM
!
! Compute to five figures the integral
!
!     1  z
!    I  I  (x+1)**3*(z+1)**2  dxdz
!   0  0 
!
! SOLUTION
! 
! This simple nested quadrature problem can be solved by recursive use
! of the integrators provided in rksuite_90. First, we define the 
! h'(x)
!
!    h'(x) = (x+1)**3     with h(0) = 0.0
!
! which is an initial value problem. The required integral can then be
! expressed as
!
!     1
!    I  h(z)*(z+1)**2 dz
!   0
!
! Similarly, we define g'(z)
!
!    g'(z) = h(z)*(z+1)**2     with g(0) = 0.0,
!
! another initial value problem. Then, the answer to the problem is 
! given by g(1)
!
! Hence we use RANGE_INTEGRATE to integrate g'(z) to compute the value
! g(1) and in the definition of g'(z) (the procedure GP) we use 
! RANGE_INTEGRATE again (recursively) to integrate h'(x) (the procedure
! HP) to compute the required values of h(z).
! 
! Given the nature of the integrand and the area of integration we expect
! the answer to be of order of magnitude 1. Hence, we use a relative
! TOLERANCE = 5.0e-5 and a value for THRESHOLD = 1.0e-10.
 
use rksuite_90
use quadrature_inner_f

implicit none

real(kind=wp) :: x_lower=0.0_wp, x_inner_got, h_of_z_got, &
   tol=5.0e-6, thresh = 1.0e-10_wp

type(rk_comm_real_0d) :: comm_inner

contains

function gp(z,g)
real(kind=wp), intent(in) :: z
real(kind=wp), intent(in) :: g
real(kind=wp) :: gp
if (z == 0.0_wp) then
   h_of_z_got = 0.0_wp
else
   call setup(comm_inner,x_lower,0.0_wp,z,tol,thresh)
   call range_integrate(comm_inner,hp,z,x_inner_got,y_got=h_of_z_got)
   call collect_garbage(comm_inner)
end if

gp = h_of_z_got * (1.0_wp + z)**2

end function gp

end module quadrature_f

program quadrature

use quadrature_f

implicit none

integer :: totf
real(kind=wp) :: z_lower=0.0_wp, z_upper=1.0_wp, z_got, g_got

type(rk_comm_real_0d) :: comm

call setup(comm,z_lower,0.0_wp,z_upper,tol,thresh)
call range_integrate(comm,gp,z_upper,z_got,y_got=g_got)

write (*,'(/a,f9.4)') ' ans = ',g_got

call statistics(comm,total_f_calls=totf)

write (*,'(/a,i10)') &
      ' The cost of the outer integration in calls to gp was', totf

call collect_garbage(comm)

end program quadrature
