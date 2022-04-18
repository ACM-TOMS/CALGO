PROGRAM TEST_SGER_RMD
  use rmd_stesttools
  implicit none
  procedure(F_interf)     :: F
  procedure(F_rmd_interf) :: F_rmd
  !  
  integer :: m, n
  type(fdata) :: testcase
  
  real    :: tol = rmd_stolerance
  integer :: rep = 1
  !
  !m = 3
  !n = 3
  !testcase = fdata([real::], [m, n], [character::])
  !call rmd_stest123(F, F_rmd, m*n + m + n + 1, m*n, tol, testcase)
  do m = 1, 5
    do n = 1, 5
      testcase = fdata([real::], [m, n], [character::])
      call rmd_stestrandom(F, F_rmd, m*n+m+n+1, m*n, tol, rep, testcase)
    end do
  end do
END PROGRAM TEST_SGER_RMD

SUBROUTINE F(testcase, u, v)
  use rmd_stesttools
  !use dispmodule
  type(fdata), intent(inout) :: testcase
  real, intent(in)  :: u(*)
  real, intent(out) :: v(*)
  !
  integer :: m, n
  real :: alpha
  real, allocatable :: x(:), y(:), A(:, :)
  !
  m     = testcase % integers(1)
  n     = testcase % integers(2)
  !
  allocate(x(m), y(n), A(m, n))
  !  
  call scopy(m*n, u, 1, A, 1)
  call scopy(m, u(m*n+1), 1, x, 1)
  call scopy(n, u(m*n+m+1), 1, y, 1)
  alpha = u(m*n+m+n+1)                    
  !
  call sger(m, n, alpha, x, 1, y, 1, A, m) ! A = alpha*x*y' + A
  call scopy(m*n, A, 1, v, 1)
  !
END SUBROUTINE F

SUBROUTINE F_RMD(testcase, u, v, ua, va)
  use rmd_stesttools
  !use dispmodule
  type(fdata), intent(inout) :: testcase
  real, intent(in)        :: u(*), v(*), va(*)
  real, intent(inout)     :: ua(*)
  !
  integer :: m, n
  real :: alpha, alphaa
  real, allocatable :: x(:), y(:), Aa(:, :), xa(:), ya(:)
  !
  m     = testcase % integers(1)
  n     = testcase % integers(2)
  !
  allocate(x(m), y(n), Aa(m, n), xa(m), ya(n))
  !
  call scopy(m, u(m*n+1), 1, x, 1)
  call scopy(n, u(m*n+m+1), 1, y, 1)
  alpha = u(m*n+m+n+1)
  call scopy(m*n, va, 1, Aa, 1)
  call scopy(m, ua(m*n+1), 1, xa, 1)
  call scopy(n, ua(m*n+m+1), 1, ya, 1)
  alphaa = ua(m*n + m + n + 1)
  !
  call sger_rmds(m, n, x, 1, y, 1, m, alphaa, Aa)
  call sger_rmd(m, n, alpha, x, 1, y, 1, m, xa, ya, Aa, '111')
  !
  ! First m*n elements of ua remain constant
  call scopy(m*n, Aa, 1, ua, 1)        ! copy back from Aa
  call scopy(m, xa, 1, ua(m*n+1), 1)   ! copy back from xa
  call scopy(n, ya, 1, ua(m*n+m+1), 1) ! copy back from ya
  ua(m*n+m+n+1) = alphaa
  !
END SUBROUTINE F_RMD
