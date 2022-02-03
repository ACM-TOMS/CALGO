module orr_somm_f

! PROBLEM
! 
! Compute a trial solution for the Orr-Sommerfeld equation for plane
! Poiseuille flow
!    
!  y'''' - 2 k**2 y'' + k**4 y - i k R [(u(t) - c)(y'' - k**2 y) - u''(t)y] = 0
!  
! y(t) is the amplitude of the stream function, R is the Reynolds number,
! k is the wave number of the disturbance, i is sqrt(-1), c is the eigenvalue,
! and, in this instance, the laminar velocity profile is u(t) = 1 - t**2
! 
! The boundary conditions are y'(0) = y'''(0) = 0 = y(1) = y'(1).
! 
! Use c = 0.06659252 - 0.01398327 i, k = 1, R = 1.0e+06,  and the 
! arbitrary initial conditions y(0) = y''(0) = 1.
! 
! SOLUTION
! 
! The real objective in solving this equation is to compute eigenvalues
! c such that the imaginary part of y'(1) = 0. See "Computational
! solution of two point BVPs via orthonormalization", MR Scott and HA Watts
! SIAM J. Num. Anal. 14, pp 40-40, 1977 and the references therein.  
! 
! We recast the problem in to a system of complex first order ODEs using 
! y1 = y, y2 = y', y3 = y'', y4 = y''', to obtain
! 
!    y1' = y2
!    y2' = y3
!    y3' = y4
!    y4' = 2 k**2 y3 - k**4 y1 + i k R [(u(t) - c)(y3 - k**2 y1) - u''(t)y1]
! 
!    given y1(0) = 1, y2(0) = 0, y3(0) = 1, y4(0) = 0,
! 
! and use the complex dependent variable facility of rksuite_90.
! 
! For such an arbitrary selection of initial conditions it is likely that
! the computed solution will "blow up". We use STEP_INTEGRATE to monitor
! the progress of the computation and if any solution component exceeds
! 1.0e+10 in magnitude the computation is abandoned.
! 
! NOTE
! 
! This example demonstrates the use of the complex dependent variable
! facility  offerred by rksuite_90. The example eg_orr_somm_r.f90
! illustrates how the same system can be handled in the more traditional
! (and more complicated) manner of treating separately the real imaginary
! parts of y as a larger system of real dependent variables.  
 
use rksuite_90_prec
 
real(kind=wp), parameter :: k=1.0_wp, k2=k*k, k4=k2*k2, r=1.0e+06_wp, kr=k*r
complex(kind=wp) :: c = (0.06659252_wp, -0.01398327_wp)

contains

function f(t,y)
real(kind=wp), intent(in) :: t
complex(kind=wp), dimension(:), intent(in) :: y
complex(kind=wp), dimension(size(y)) :: f

f(1:3) = y(2:4)
f(4) = 2.0_wp*k2*y(3) - k4*y(1) +  cmplx(0.0_wp,1.0_wp)*kr* &
       ( (u(t)-c)*(y(3)-k**2*y(1)) - upp(t)*y(1) )

end function f

function u(t)
real(kind=wp), intent(in) :: t
real(kind=wp) :: u
u = 1.0_wp - t**2
end function u

function upp(t)
real(kind=wp), intent(in) :: t
real(kind=wp) :: upp
upp = -2.0_wp
end function upp

end module orr_somm_f

program orr_somm

use orr_somm_f
use rksuite_90

implicit none

integer :: steps
real(kind=wp) :: t_start=0.0_wp, t_end=1.0_wp, tol=5.0e-5_wp, t_now

real(kind=wp), dimension(4) :: thres = 1.0e-10_wp
complex(kind=wp), dimension(4) :: y_start, y_now, yderiv_now

type(rk_comm_complex_1d) :: comm

y_start(1) = 1.0_wp
y_start(2) = 0.0_wp
y_start(3) = 1.0_wp
y_start(4) = 0.0_wp

call setup(comm,t_start,y_start,t_end,tol,thres,task='step')

do

   call step_integrate(comm,f,t_now,y_now=y_now,yderiv_now=yderiv_now)

   if (maxval(abs(y_now(:))) > 1.0e+10_wp .or. t_now == t_end) exit

end do

call statistics(comm,num_succ_steps=steps)

write (*,'(a,e11.4,a,i4,a,e11.4)') &
   " At t = ",t_now," in ",steps," steps with imag(y'(t)) = ",aimag(y_now(2))

call collect_garbage(comm)

end program orr_somm







