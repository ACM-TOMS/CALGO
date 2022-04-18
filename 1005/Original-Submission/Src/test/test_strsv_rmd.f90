PROGRAM TEST_STRSV_RMD
  use rmd_stesttools
  implicit none
  procedure(F_interf)     :: F
  procedure(F_rmd_interf) :: F_rmd
  type(fdata) :: testcase
  integer     :: n, sizeIn, sizeOut, ul, t, d
  real    :: tol = rmd_stolerance
  integer :: esa, incx, i, nn, jsel
  character :: uplo(2), trans(2), diag(2), sel*2, wrksel*2 = '01'
  real, allocatable :: A(:, :), x(:), u(:), wrk(:), Aa(:, :), xa(:)
  !
  uplo  = ['U', 'L']
  trans = ['N', 'T']
  diag  = ['N', 'U']
  !
  do n = 1, 5
    nn = n*(n+1)/2
    do ul = 1, 2
      do t = 1, 2
        do d = 1, 2
          do esa = 0, 3, 3
            do incx = 1, 3, 2
              sizeIn  = n + nn
              sizeOut = n
              !
              allocate(A(n+esa, n), x(n), u(sizeIn), Aa(n+esa, n), xa(n), wrk(n))
              call random_number(A(1:n, 1:n))
              call random_number(x)
              call random_number(Aa(1:n, 1:n))
              call random_number(xa)
              do i = 1, n
                A(i, i) = rmd_srandom() + 0.5
                if(rmd_srandom() < 0.5) A(i, i) = -A(i, i)
              end do
              ! test that workspsace is not used when not needed:
              do jsel = 1, 2
                sel = wrksel(jsel:jsel) // '1'
                wrk = 37.0
                call strsv_rmd(uplo(ul), trans(t), diag(d), n, A, n+esa, x, 1, Aa, xa, wrk, sel)
                if (any(wrk /= 37.0)) then
                  print *, "TEST_STRSV_RMD: Workspace changed when it shouldn't"
                  stop 1
                endif
              enddo
              ! main test of correct adjoint:
              u(nn + 1 : nn + n) = x
              call rmd_strishape_m2v(uplo(ul), n, u, A, n+esa)
              testcase = fdata([real::], [n, esa, incx], [uplo(ul), trans(t), diag(d)])
              call rmd_stestf(F, F_rmd, u, sizeIn, sizeOut, tol, testcase)
              deallocate(A, x, u, Aa, xa, wrk)
            end do
          end do
        end do
      end do
    end do
  end do
END PROGRAM TEST_STRSV_RMD

SUBROUTINE F(testcase, u, v)
  use rmd_stesttools
  implicit none
  ! Arguments
  type(fdata), intent(inout) :: testcase
  real, intent(in)  :: u(*)
  real, intent(out) :: v(*)
  !
  ! Local variables
  integer :: n, sizeA, esa, incx
  character uplo, trans, diag
  !  
  real, allocatable :: A(:, :), x(:)
  !
  uplo  = testcase % characters(1)
  trans = testcase % characters(2)
  diag  = testcase % characters(3)
  n     = testcase % integers(1)
  esa = testcase % integers(2)
  incx = testcase % integers(3)
  !
  sizeA = (n*(n+1))/2
  !  
  allocate(A(n+esa, n), x(n*incx))
  !
  call rmd_strishape_v2m(uplo, n, u, A, n+esa) ! A
  call scopy(n, u(sizeA+1), 1, x, incx)           ! x
  !
  ! BLAS operation
  call strsv(uplo, trans, diag, n, A, n+esa, x, incx)
  !  
  ! Copy back from x to v
  call scopy(n, x, incx, v, 1)
  !
END SUBROUTINE F

SUBROUTINE F_RMD(testcase, u, v, ua, va)
  use rmd_stesttools
  implicit none
  ! Arguments
  type(fdata), intent(inout) :: testcase
  real, intent(in)     :: u(*), va(*), v(*)
  real, intent(inout)  :: ua(*)
  !  
  ! Local variables
  real, allocatable :: A(:, :), x(:), Aa(:, :), xa(:)
  integer :: n, sizeA, esa, incx
  character uplo, trans, diag
  real :: dummy(1)
  !
  uplo   = testcase % characters(1)
  trans  = testcase % characters(2)
  diag   = testcase % characters(3)
  n      = testcase % integers(1)
  esa = testcase % integers(2)
  incx = testcase % integers(3)
  !  
  allocate(A(n+esa, n), x(n*incx), Aa(n+esa, n), xa(n*incx))
  !
  sizeA = (n*(n+1))/2
  !
  ! Copy inital values to local variables
  call rmd_strishape_v2m(uplo, n, u, A, n+esa)
  call rmd_strishape_v2m(uplo, n, ua, Aa, n+esa)
  call scopy(n, v, 1, x, incx)                   
  call scopy(n, va, 1, xa, incx)                   
  !
  ! RMD operation
  call strsv_rmd(uplo, trans, diag, n, A, n+esa, x, incx, Aa, xa, dummy, '11')
  !  
  ! Copy result to output variables
  call rmd_strishape_m2v(uplo, n, ua, Aa, n+esa) ! Aa to ua
  call scopy(n, xa, incx, ua(sizeA+1), 1)           ! xa to ua(sizeA+1)
  !  
END SUBROUTINE F_RMD
