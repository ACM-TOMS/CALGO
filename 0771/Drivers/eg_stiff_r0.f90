module stiff_r0_f

! PROBLEM
!
! Demonstrate that a scalar ODE can be stiff
!
! SOLUTION
! 
! This problem is considered in detail in 
!    L.F. Shampine,"numerical solution of Ordinary Differential Equations",
!       Chapman and Hall, 1994, p384
! We use the example there
!
!    y' = J ( y - p(x) ) + p'(x),    y(0) = A
!
! which has solution y(x) = ( A - p(0) ) exp(Jx) + p(x). For p(x)
! slowly varying and J large and negative this problem is stiff. We choose
! p(x) = cos(x), A = 0 and J = -10**(n) with n=3.
!
! We use RANGE_INTEGRATE and switch of the printing of error messages by
! setting MESSAGE=.FALSE. in the call to SETUP. 
! 
! NOTE
!
! This problem is for illustrative purposes. We use the real scalar 
! dependent variable facility of rksuite_90
 
use rksuite_90_prec

real(kind=wp) :: lam=-1000.0_wp, a=0.0_wp

contains

function f(t,y)
real(kind=wp), intent(in) :: t
real(kind=wp), intent(in) :: y
real(kind=wp) :: f

f = lam*(y-cos(t)) - sin(t)

end function f

function exact(t)
real(kind=wp), intent(in) :: t
real(kind=wp)  :: exact

exact = (a-cos(0.0_wp))*exp(lam*t) + cos(t)

end function exact

end module stiff_r0_f

program stiff_r0

use stiff_r0_f
use rksuite_90

implicit none

integer :: flag
real(kind=wp) :: t_start=0.0_wp, t_end=10.0_wp, tol=1.0e-4_wp, t_got, &
   thres = 1.0_wp, y_start, y_got

type(rk_comm_real_0d) :: comm

y_start = a

call setup(comm,t_start,y_start,t_end,tol,thres,message=.false.)

call range_integrate(comm,f,t_end,t_got,y_got=y_got,flag=flag)

write (*,'(a,f6.3)') ' Have integrated to t = ',t_got

if (flag==4) then
    write (*,'(a)') ' The problem has been diagnosed as stiff'
else if (flag/=1) then
    write (*,'(a/a)') ' The problem was not diagnosed as stiff !!! ', &
                      ' did you change the program...????'
end if

call collect_garbage(comm)

end program stiff_r0
