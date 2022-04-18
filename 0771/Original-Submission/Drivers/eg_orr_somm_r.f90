module orr_somm_r_f
 
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
! We recast the problem in to a system of real first order ODEs using 
! y1 = real(y), y3 = real(y'), y5 = real(y''), y7 = real(y'''),
! y2 = imag(y), y4 = imag(y'), y6 = imag(y''), y8 = imag(y'''), to obtain
! 
!    y1' = y3
!    y2' = y4
!    y3' = y5
!    y4' = y6
!    y5' = y7
!    y6' = y8
!    y7' = 2 k**2 y5 - k**4 y1 - k R * 
!              [(u(t) - cr)(y6 - k**2 y2) + ci (y5 - k**2 y1) - u''(t) y2 ]
!    y8' = 2 k**2 y6 - k**4 y2 + k R *
!              [(u(t) - cr)(y5 - k**2 y1) - ci (y6 - k**2 y2) - u''(t) y1 ]
! 
!    given y1(0) = 1, y3(0) = 0, y5(0) = 1, y7(0) = 0,
!          y2(0) = 0, y4(0) = 0, y6(0) = 0, y8(0) = 0,
! 
!    where cr = real(c) and ci = imag(c).
! 
! For such an arbitrary selection of initial conditions it is likely that
! the computed solution will "blow up". We use STEP_INTEGRATE to monitor
! the progress of the computation and if any solution component exceeds
! 1.0e+10 in magnitude the computation is abandoned.
! 
! NOTE
! 
! This example demonstrates the use of the real dependent variable
! facility  offerred by rksuite_90 when in fact the complex dependent
! variable facility would have been more appropriate. The example 
! eg_orr_somm.f90 illustrates how the same system can be handled more 
! conveniently.
 
use rksuite_90_prec

real(kind=wp), parameter :: k=1.0_wp, k2=k*k, k4=k2*k2, r=1.0e+06_wp, kr=k*r, &
   cr=0.06659252_wp, ci=-0.01398327_wp

contains

function f(t,y)
real(kind=wp), intent(in) :: t
real(kind=wp), dimension(:), intent(in) :: y
real(kind=wp), dimension(size(y)) :: f

real(kind=wp) :: ucr, y2rk2yr, y2ik2yi

f(1:6) = y(3:8)

ucr = u(t) - cr
y2rk2yr = y(5) - k2*y(1)
y2ik2yi = y(6) - k2*y(2)

f(7) = 2.0_wp*k2*y(5) - k4*y(1) - k*r* & 
       ( ucr*y2ik2yi + ci*y2rk2yr - upp(t)*y(2) )

f(8) = 2.0_wp*k2*y(6) - k4*y(2) + k*r* &
       ( ucr*y2rk2yr - ci*y2ik2yi - upp(t)*y(1) )
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

end module orr_somm_r_f

program orr_somm_r

use orr_somm_r_f
use rksuite_90

implicit none

integer :: steps
real(kind=wp) :: t_start=0.0_wp, t_end=1.0_wp, tol=5.0e-5_wp, t_now
   
real(kind=wp), dimension(8) :: thres = 1.0e-10_wp
real(kind=wp), dimension(8) :: y_start, y_now, yderiv_now

type(rk_comm_real_1d) :: comm

y_start(:) = 0.0_wp
y_start(1) = 1.0_wp
y_start(5) = 1.0_wp

call setup(comm,t_start,y_start,t_end,tol,thres,task='step')

do

   call step_integrate(comm,f,t_now,y_now=y_now,yderiv_now=yderiv_now)

   if (maxval(abs(y_now(:))) > 1.0e+10_wp .or. t_now == t_end) exit

end do

call statistics(comm,num_succ_steps=steps)

write (*,'(a,e11.4,a,i4,a,e11.4)') &
   " At t = ",t_now," in ",steps," steps and imag(y'(t)) = ",y_now(4)

call collect_garbage(comm)

end program orr_somm_r
