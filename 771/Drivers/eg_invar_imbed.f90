module invar_imbed_fw

use rksuite_90_prec


contains

function g(t,y)
real(kind=wp), intent(in) :: t
real(kind=wp), dimension(:), intent(in) :: y
real(kind=wp), dimension(size(y)) :: g

g(:) = (/ 1.0_wp - y(1)**2, -y(1)*y(2) /)

end function g

end module invar_imbed_fw

module invar_imbed_bw

use rksuite_90
use invar_imbed_fw

real(kind=wp) :: u, v

real(kind=wp) :: t_fw_start=0.0_wp, t_fw_end, t_fw_got, tol=1.0e-5_wp
real(kind=wp), dimension(2) :: thres(2)=1.0e-20_wp, &
    y_fw_start, y_fw_got
type(rk_comm_real_1d) :: comm_fw

contains

function f(t,y)

real(kind=wp), intent(in) :: t
real(kind=wp), intent(in) :: y
real(kind=wp) :: f

t_fw_end = t
u = 0.0_wp; v = 1.0_wp; y_fw_start = (/u, v/)

if (t_fw_end /= t_fw_start) then
   call setup(comm_fw,t_fw_start,y_fw_start,t_fw_end,tol,thres)
   call range_integrate(comm_fw,g,t_fw_end,t_fw_got,y_got=y_fw_got)
   call collect_garbage(comm_fw)
   u = y_fw_got(1); v = y_fw_got(2)
end if

f = u*y + v

end function f

end module invar_imbed_bw

program invar_imbed

use invar_imbed_bw

implicit none

real(kind=wp) :: t_bw_start=5.0_wp, t_bw_end, y_bw_start, &
   y_bw_now, t_bw_now

type(rk_comm_real_0d) :: comm_bw

t_fw_end = t_bw_start; y_fw_start = (/0.0_wp, 1.0_wp/)

call setup(comm_fw,t_fw_start,y_fw_start,t_fw_end,tol,thres)
call range_integrate(comm_fw,g,t_fw_end,t_fw_got,y_got=y_fw_got)
call collect_garbage(comm_fw)

u = y_fw_got(1); v = y_fw_got(2)
y_bw_start = (exp(-t_fw_end) - v)/u
    
t_bw_end = t_fw_start

call setup(comm_bw,t_bw_start,y_bw_start,t_bw_end,tol,thres(1),task='s')

write (*,'(a/)') '      t         approx          exact'
do 
   call step_integrate(comm_bw,f,t_bw_now,y_now=y_bw_now)
   write (*,'(1x,f6.3,3x,e12.5,3x,e12.5)') t_bw_now, u*y_bw_now+v, &
      exp(-t_bw_now)
   if (t_bw_now==t_bw_end) exit
end do   

call collect_garbage(comm_bw)

end program invar_imbed
