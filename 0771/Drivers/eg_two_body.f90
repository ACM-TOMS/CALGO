module two_body_f

! PROBLEM
! 
! Integrate a two body problem.  The equations for the coordinates
! (x(t),y(t)) of one body as functions of time t in a suitable frame of
! reference are
! 
!     x'' = - x/r**3,
!     y'' = - y/r**3,   where   r = SQRT(x**2 + y**2).
! 
! The initial conditions lead to elliptic motion with eccentricity ECC.
! This parameter will be taken to be 0.9.
! 
!     x(0) = 1-ECC,     x'(0) = 0,
!     y(0) = 0,         y'(0) = SQRT((1+ECC)/(1-ECC)).
! 
! An accurate solution that shows the general behavior of the orbit is
! desired.  The coordinates will be returned at every time step in [0,20].
! This is a standard test problem for which there is a semi-analytical
! solution.  It will be compared to the computed solution at the end of
! the interval.
! 
! SOLUTION
! 
! For illustrative purposes we shall pose the system using a two dimensional
! array of unkowns
!        
!      ( x   x' )
!      ( y   y' )
! 
! to obtain
! 
!      ( x   x' )'  =  ( x'  -x/r**3 )
!      ( y   y' )      ( y'  -y/r**3 )
! 
! Since it is the general behavior of the solution that is desired, it is
! best to let the integrator choose where to provide answers.  It will
! produce answers more frequently where the solution changes more rapidly.
! Because the solution is inspected at every step, the task is to use
! STEP_INTEGRATE.
! 
! On physical grounds the solution is expected to be somewhat unstable
! when one body approaches the other. To obtain an accurate solution, a
! stringent relative error tolerance should be imposed -- TOL = 1.0D-10
! will be used.  At a tolerance this stringent the highest order pair,
! METHOD = 'H', is likely to be the most efficient choice. This method is
! inefficient when it is to produce answers at a great many specific
! points. It is most effective when used as in this template.  The
! solution components are expected to be of order 1, so threshold values
! of THRESHOLDS = 1.0D-13 are reasonable.  When a solution component is smaller
! in magnitude than this threshold, the code will control the local error
! to be no more than TOL*THRESHOLDS = 1.0D-23.  The reasonableness of this
! choice will be monitored by printing out the maximum value seen for each
! solution component in the course of the integration.  Error and warning
! messages, if any, will be printed out.
! 
! This is the standard test problem D5 of T.E. Hull, W.H. Enright,
! B.M. Fellen, and A.E. Sedgwick, "Comparing Numerical Methods
! for Ordinary Differential Equations," SIAM Journal on Numerical
! Analysis, Vol. 9, pp. 603-637, 1972.  The analytical solution
! in terms of the numerical solution of Kepler's equation can be
! found there as well as in most discussions of the two body
! problem.  The results for the particular choice of eccentricity,
! initial conditions, and interval  are providedin TRUE_Y_AT_TEND.
! 
! Global error assessment has been selected as a check on the reliability
! of the results. From the results generated it will be seen that the
! global error at the end of the run is about 2.3E-9, rather bigger than
! the local error tolerance TOL.  This illustrates the point that at best
! one can anticipate global errors comparable to the tolerance.  In point
! of fact, this problem is unstable at some points of the integration and
! the global error assessment reveals that the worst global error is
! considerably worse than the error at the end -- an example of the value
! of the global error assessment capability.

use rksuite_90_prec

contains

function f(t,y)
real(kind=wp), intent(in) :: t
real(kind=wp), dimension(:,:), intent(in) :: y
real(kind=wp), dimension(size(y,1),size(y,2)) :: f

real(kind=wp) :: r_cubed

r_cubed = sqrt( sum ( y(:,1)**2 ) )**3
f(:,1) =  y(:,2);
f(:,2) = -y(:,1)/r_cubed

end function f

end module two_body_f

program two_body

use two_body_f
use rksuite_90

implicit none

integer :: flag, i, totf
real(kind=wp) :: t_start=0.0_wp, t_end=20.0_wp, tol=1.0e-10_wp, t_now, &
     ecc=0.9_wp, true_error, waste, max_error, t_max_error
real(kind=wp), dimension(2,2) :: y_start, thres = 1.0e-13_wp, y_now, &
     y_maxvals, true_y_at_t_end, assessed_error

type(rk_comm_real_2d) :: comm

true_y_at_t_end(:,1) = (/-1.29526625098758_wp,  0.400393896379232_wp /)
true_y_at_t_end(:,2) = (/-0.67753909247075_wp, -0.127083815427869_wp /)

y_start(:,1) = (/ 1.0_wp - ecc, 0.0_wp /)
y_start(:,2) = (/ 0.0_wp,  sqrt((1.0_wp+ecc)/(1.0_wp-ecc)) /)

call setup(comm,t_start,y_start,t_end,tol,thres,method='high',task='step', &
           error_assess=.true.)

do
   call step_integrate(comm,f,t_now,y_now=y_now,flag=flag)
!  write (*,'(5e14.5') t_now, y_now
   if (t_now==t_end .or. flag > 3) exit
end do

call statistics(comm,y_maxvals=y_maxvals,total_f_calls=totf,waste=waste)

write (*,'(/a,1p2e9.2/21x,1p2e9.2)') &
      '             ymax    ', (y_maxvals(i,1:2),i=1,2)
write (*,'(/a,1pe10.2,0p/a,i10/a,f10.2)') &
      ' The integration reached                      ', t_now, &
      ' The cost of the integration in calls to F was', totf, &
      ' The fraction of failed steps was             ', waste

if (t_now==t_end) then
   true_error = maxval(abs(true_y_at_t_end - y_now)/max(abs(y_now),thres))
   call global_error(comm,rms_error=assessed_error,max_error=max_error, &
                     t_max_error=t_max_error)
   write (*,'(/(a,1p,e9.2/a,e9.2/a,e9.2,a,e9.2))') &
      ' At t = 20, the true error is',true_error, &
      '            assessed error is',maxval(assessed_error), &
      '            maximum error was',max_error,' at ',t_max_error
end if

call collect_garbage(comm)

end program two_body
