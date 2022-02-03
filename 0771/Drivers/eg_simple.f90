module simple_f

! PROBLEM
!
! Compute about four correct figures in the solution of
!
!     y'' = -y
!
! on the range [0,2*pi] at intervals of length pi/4, given that y(0)=0,
! y'(0)=1.
!
! SOLUTION
! 
! Let y1 = y and y2 = y' to obtain the first order system
! 
!    y1' =   y2     with initial values   y1 = 0
!    y2' = - y1                           y2 = 1
! 
! This is a "usual task" that is appropriately solved with
! RANGE_INTEGRATE.  Although the code controls the local error rather than
! the true error, it is "tuned" so that the true error will be comparable
! to the local error tolerance for typical problems. Thus a relative error
! tolerance TOL of 5.0D-5 is appropriate.  In this range of tolerances,
! METHOD = 'M' (the default) is likely to be the most efficient choice.
! The solution components are expected to get as large as 1.0D0.  With
! this in mind, solution components smaller than, say, 1.0D-10 are not
! very interesting, and requiring five correct figures then is not worth
! the cost. For this reason, the threshold values are specified as
! THRESHOLDS = 1.0D-10.  When a solution component is smaller than this
! threshold, the code will control the local error to be no more than
! TOL*THRESHOLDS = 5.0D-15.  Answers will be computed at equally spaced
! points, and the true values of y and y' will be printed out for
! comparison.
 
use rksuite_90_prec

contains

function f(t,y)
real(kind=wp), intent(in) :: t
real(kind=wp), dimension(:), intent(in) :: y
real(kind=wp), dimension(size(y)) :: f

f(:) = (/ y(2), -y(1) /)

end function f

end module simple_f

program simple

use simple_f
use rksuite_90

implicit none

integer :: nout, totf, l
real(kind=wp) :: t_start=0.0_wp, t_end, tol=5.0e-5_wp, t_got, &
   t_want, pi, t_inc
real(kind=wp), dimension(2) :: y_start = (/ 0.0_wp,1.0_wp /), &
   thres = 1.0e-10_wp, y_got, y_maxvals

type(rk_comm_real_1d) :: comm

pi = 4.0_wp*atan(1.0_wp)
t_end = 2.0_wp*pi

call setup(comm,t_start,y_start,t_end,tol,thres)

nout = 8; t_inc = (t_end-t_start)/nout

do l = 1, nout

   t_want = t_end + (l-nout)*t_inc
   call range_integrate(comm,f,t_want,t_got,y_got=y_got)
   write (*,'(1x,f6.3,3x,f9.4,3x,f9.4)') t_got, y_got(1), sin(t_got)
   write (*,'(1x,9x,f9.4,3x,f9.4/)')            y_got(2), cos(t_got)

end do

call statistics(comm,y_maxvals=y_maxvals,total_f_calls=totf)

write (*,'(/a,1p2e9.2)') &
       '             ymax    ', y_maxvals
write (*,'(/a,i10)') &
      ' The cost of the integration in calls to F was', totf

call collect_garbage(comm)

end program simple
