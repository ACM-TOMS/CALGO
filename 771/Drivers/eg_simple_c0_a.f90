module simple_c0_a_f

! PROBLEM
!
! Compute about seven correct figures in the solution of
!
!     y' = i y
!
! on the range [0,2*pi] at intervals of length pi/4, given that y(0)=(1,0).
!
! SOLUTION
!
! This is a variant of "eg_simple_c0.f90", where here somewhat more
! accuracy is required. As for "eg_simple_c0.f90", this is a "usual task"
! that is more appropriately solved with RANGE_INTEGRATE.  However, for
! illustrative purposes, we shall use STEP_INTEGRATE in conjunction with
! RESET_T_END.
! Although the code controls the local error rather than
! the true error, it is "tuned" so that the true error will be comparable
! to the local error tolerance for typical problems. Thus a relative error
! tolerance TOL of 5.0D-8 is appropriate.  In this range of tolerances,
! METHOD = 'H' is likely to be the most efficient choice.
! The solution components are expected to get as large as 1.0D0.  With
! this in mind, solution components smaller than, say, 1.0D-15 are not
! very interesting, and requiring seven correct figures then is not worth
! the cost. For this reason, the threshold values are specified as
! THRESHOLDS = 1.0D-15.  When a solution component is smaller than this
! threshold, the code will control the local error to be no more than
! TOL*THRESHOLDS = 5.0D-23. The true value of y will be printed out for
! comparison.
!
! NOTE
!
! This problem is similar to that in simple.f90. Here, for illustrative
! purposes, we use the complex scalar dependent variable facility of
! rksuite_90.
 
use rksuite_90_prec

contains

function f(t,y)
real(kind=wp), intent(in) :: t
complex(kind=wp), intent(in) :: y
complex(kind=wp) :: f

f = cmplx(0.0_wp,1.0_wp) * y

end function f

end module simple_c0_a_f

program simple_c0_a

use simple_c0_a_f
use rksuite_90

implicit none

integer :: nout, totf, l
real(kind=wp) :: t_start=0.0_wp, t_end_final, tol=5.0e-5_wp, t_now, &
   t_want, pi, t_inc, thres = 1.0e-10_wp, y_maxvals
complex(kind=wp) :: y_start, y_now

type(rk_comm_complex_0d) :: comm

y_start = cmplx(1.0_wp,0.0_wp)
pi = 4.0_wp*atan(1.0_wp)
t_end_final = 2.0_wp*pi

nout = 8; t_inc = (t_end_final-t_start)/nout
l = -7; t_want = t_end_final + l*t_inc

call setup(comm,t_start,y_start,t_want,tol,thres,task='step',method='high')

do

   call step_integrate(comm,f,t_now,y_now=y_now)

   if (t_want == t_now) then
      write (*,'(1x,f6.3,3x,f12.7,3x,f12.7)') t_want, real(y_now), cos(t_want)
      write (*,'(1x,9x,f12.7,3x,f12.7/)')            aimag(y_now), sin(t_want)
      l = l + 1; t_want = t_end_final + l*t_inc
      call reset_t_end(comm,t_want)
   end if

   if (t_now == t_end_final) exit
end do

call statistics(comm,y_maxvals=y_maxvals,total_f_calls=totf)

write (*,'(/a,1pe9.2)') &
       '             ymax    ', y_maxvals
write (*,'(/a,i10)') &
      ' The cost of the integration in calls to F was', totf

call collect_garbage(comm)

end program simple_c0_a
