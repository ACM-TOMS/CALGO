PROGRAM TEST_SAXPY_RMD
  use rmd_stesttools
  implicit none
  procedure(F_interf)     :: F
  procedure(F_rmd_interf) :: F_rmd
  ! External functions
  ! Changing values
  integer :: m, incx, incy
  type(fdata) :: testcase
  ! Constants
  real    :: tol = rmd_stolerance
  integer :: rep = 1
  !    
  do m = 1, 5
    do incx = 1, 3, 2
      do incy = 1, 3, 2
        testcase = fdata([real::], [m, incx, incy], [character::])
        call rmd_stestrandom(F, F_rmd, 2*m+1, m, tol, rep, testcase)
      end do
    end do
  end do
END PROGRAM TEST_SAXPY_RMD

SUBROUTINE F(testcase, u, v)
  use rmd_stesttools
  implicit none
  type(fdata), intent(inout) :: testcase
  real, intent(in)  :: u(*)
  real, intent(out) :: v(*)
  !
  integer :: m, incx, incy
  real    :: alpha
  real, allocatable :: x(:), y(:)
  !
  m     = testcase % integers(1)
  incx  = testcase % integers(2)
  incy  = testcase % integers(3)
  !  
  allocate(x(m*incx), y(m*incy))
  !
  call scopy(m, u(1), 1, x, incx)   ! x = u(1:m)
  call scopy(m, u(m+1), 1, y, incy)   ! y = u(m+1:2*m)
  alpha = u(2*m+1)
  call saxpy(m, alpha, x, incx, y, incy) ! y = alpha*x + y
  call scopy(m, y, incy, v, 1)        ! v = y
  !
END SUBROUTINE F

SUBROUTINE F_RMD(testcase, u, v, ua, va)
  use rmd_stesttools
  implicit none
  type(fdata), intent(inout) :: testcase
  real, intent(in)        :: u(*), v(*), va(*)
  real, intent(inout)     :: ua(*)
  !
  integer :: m, incx, incy
  real    :: alpha, alphaa
  real, allocatable :: x(:), y(:), xa(:), ya(:)
  !
  m = testcase % integers(1)
  incx = testcase % integers(2)
  incy = testcase % integers(3)
  !  
  allocate(x(m*incx), y(m*incy), xa(m*incx), ya(m*incy))
  !
  call scopy(m, u(1), 1, x, incx) ! x = u(1:m)
  call scopy(m, ua, 1, xa, incx)  ! xa = ua(1:m)
  call scopy(m, va, 1, ya, incy)  ! ya = va
  alpha = u(2*m+1)
  alphaa = ua(2*m+1)
  !
  call saxpy_rmds(m, x, incx, incy, alphaa, ya) ! Note: Only xa gets updated
  call saxpy_rmd(m, alpha, incx, incy, xa, ya) ! Note: Only xa gets updated
  !
  call scopy(m, xa, incx, ua, 1)  ! ua(1:m) = xa
  call scopy(m, ya, incy, ua(m+1), 1)
  ua(2*m+1) = alphaa
  !
END SUBROUTINE F_RMD
