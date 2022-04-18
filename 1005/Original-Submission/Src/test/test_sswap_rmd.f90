PROGRAM TEST_SSWAP_RMD
  use rmd_stesttools
  implicit none
  procedure(F_interf)     :: F
  procedure(F_rmd_interf) :: F_rmd
  ! External functions
  !
  ! Changing values
  integer :: n, incx, incy
  type(fdata) :: testcase
  ! Constants
  real    :: tol = rmd_stolerance
  integer :: rep = 1
  !    
  do n = 1, 5
    do incx = 1, 3, 2
    do incy = 1, 3, 2
      testcase = fdata([real::], [n, incx, incy], ['0'])
      call rmd_stestrandom(F, F_rmd, 2*n, 2*n, tol, rep, testcase)
    end do
    end do
  end do
END PROGRAM TEST_SSWAP_RMD

SUBROUTINE F(testcase, u, v)
  use rmd_stesttools
  implicit none
  type(fdata), intent(inout) :: testcase
  real, intent(in)  :: u(*)
  real, intent(out) :: v(*)
  !
  real, allocatable :: x(:), y(:)
  integer :: n, incx, incy
  !  
  n = testcase % integers(1)
  incx = testcase % integers(2)
  incy = testcase % integers(3)
  !
  allocate(x(n*incx), y(n*incy))
  !  
  call scopy(n, u, 1, x, incx)
  call scopy(n, u(n+1), 1, y, incy)
  !  
  call sswap(n, x, incx, y, incy)
  !  
  call scopy(n, x, incx, v, 1)
  call scopy(n, y, incy, v(n+1), 1)
  !
END SUBROUTINE F

SUBROUTINE F_RMD(testcase, u, v, ua, va)
  use rmd_stesttools
  implicit none
  type(fdata), intent(inout) :: testcase
  real, intent(in)        :: u(*), v(*), va(*)
  real, intent(inout)     :: ua(*)
  !
  integer :: n, incx, incy
  real, allocatable ::  xa(:), ya(:)
  !
  n = testcase % integers(1)
  incx = testcase % integers(2)
  incy = testcase % integers(3)
  !  
  allocate(xa(n*incx), ya(n*incy))
  !
  call scopy(n, va, 1, xa, incx) ! xa = va(1:n)
  call scopy(n, va(n+1), 1, ya, incy) ! ya = va(n+1:2*n)
  !
  call sswap_rmd(n, incx, incy, xa, ya)
  !
  call scopy(n, xa, incx, ua, 1) ! ua(1:n)     = xa
  call scopy(n, ya, incy, ua(n+1), 1) ! ua(n+1:2*n) = ya
  !
END SUBROUTINE F_RMD
