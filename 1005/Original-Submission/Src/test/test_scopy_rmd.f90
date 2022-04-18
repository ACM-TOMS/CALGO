PROGRAM TEST_SCOPY_RMD
  use rmd_stesttools
  implicit none
  procedure(F_interf)     :: F
  procedure(F_rmd_interf) :: F_rmd
  ! External functions
  ! Changing values
  integer :: n
  type(fdata) :: testcase
  ! Constants
  real    :: tol = rmd_stolerance
  integer :: rep = 1, incx, incy
  !    
  do n = 1, 5
    do incx = 1, 3, 2
    do incy = 1, 3, 2
      testcase = fdata([real::], [n, incx, incy], ['0'])
      call rmd_stestrandom(F, F_rmd, n, n, tol, rep, testcase)
    end do
    end do
  end do
END PROGRAM TEST_SCOPY_RMD

SUBROUTINE F(testcase, u, v)
  use rmd_stesttools
  implicit none
  type(fdata), intent(inout) :: testcase
  real, intent(in)  :: u(*)
  real, intent(out) :: v(*)
  !
  integer :: n
  !  
  n = testcase % integers(1)
  call scopy(n, u, 1, v, 1)
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
  real, allocatable :: x(:), y(:), xa(:), ya(:)
  !
  n = testcase % integers(1)
  incx = testcase % integers(2)
  incy = testcase % integers(3)
  !  
  allocate(x(n*incx), y(n*incy), xa(n*incx), ya(n*incy))
  !
  call scopy(n, u, 1, x, incx)   ! x = u(1:n)
  call scopy(n, ua, 1, xa, incx) ! xa = ua(1:n)
  call scopy(n, va, 1, ya, incy) ! ya = va
  !
  call scopy_rmd(n, incx, incy, xa, ya) ! Note: Only xa gets updated
  !
  call scopy(n, xa, incx, ua, 1) ! ua(1:n) = xa
  !
END SUBROUTINE F_RMD
