PROGRAM TEST_STPSV_RMD
  use rmd_stesttools
  implicit none
  procedure(F_interf)     :: F
  procedure(F_rmd_interf) :: F_rmd
  !
  type(fdata) :: testcase
  integer     :: n, sizeIn, sizeOut, u, t, d
  real    :: tol = rmd_stolerance
  integer :: incx, i, j
  character :: uplo(2), trans(2), diag(2)
  real, allocatable :: arg(:)
  !
  uplo  = ['L', 'U']
  trans = ['N', 'T']
  diag  = ['N', 'U']
  !
  do n = 1, 5
    do u = 1, 2
      do t = 1, 2
        do d = 1, 2
          do incx = 1, 3, 2
            sizeIn  = n + (n*(n+1))/2
            sizeOut = n
            allocate(arg(sizeIn))
            !
            do i = 1, sizeIn
              arg(i) = rmd_srandom()*4.0 - 2.0
            end do
            !
            if(u == 2) then 
              j = 0
              do i = 1, n
                j = j + i
                arg(j) = rmd_srandom() + 1
                if(rmd_srandom() < 0.5) arg(j) = -arg(j)
              end do
            else 
              j = 0
              do i = 1, n
                j = j + i
                arg((n*(n+1))/2 + 1 - j) = rmd_srandom() + 1
                if(rmd_srandom() < 0.5) arg(j) = -arg(j)
              end do
            end if
            !
            testcase = fdata([real::], [n, incx], [uplo(u), trans(t), diag(d)])
            call rmd_stestf(F, F_rmd, arg, sizeIn, sizeOut, tol, testcase)
            !
            deallocate(arg)
          end do
        end do
      end do
    end do
  end do
END PROGRAM TEST_STPSV_RMD

SUBROUTINE F(testcase, u, v)
  use rmd_stesttools
  implicit none
  ! Arguments
  type(fdata), intent(inout) :: testcase
  real, intent(in)  :: u(*)
  real, intent(out) :: v(*)
  !
  ! Local variables
  integer :: n, sizeAP, incx
  character uplo, trans, diag
  !  
  real, allocatable :: AP(:), x(:)
  !
  uplo  = testcase % characters(1)
  trans = testcase % characters(2)
  diag  = testcase % characters(3)
  n     = testcase % integers(1)
  incx = testcase % integers(2)
  !
  sizeAP = (n*(n+1))/2
  !  
  allocate(AP(sizeAP), x(n*incx))
  !
  call scopy(sizeAP, u, 1, AP, 1) ! AP
  call scopy(     n, u(sizeAP+1), 1, x, incx) ! x
  !
  ! BLAS operation
  call stpsv(uplo, trans, diag, n, AP, x, incx)
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
  real, allocatable :: AP(:), x(:), APa(:), xa(:)
  integer :: n, sizeAP, incx
  character uplo, trans, diag
  real dummy(1)
  !
  uplo   = testcase % characters(1)
  trans  = testcase % characters(2)
  diag   = testcase % characters(3)
  n      = testcase % integers(1)
  incx = testcase % integers(2)
  sizeAP = (n*(n+1))/2
  !  
  allocate(AP(sizeAP), x(n*incx), APa(sizeAP), xa(n*incx))
  !
  ! Copy inital values to local variables
  call scopy(sizeAP, u, 1, AP, 1) ! AP
  call scopy(sizeAP, ua, 1, APa, 1) ! APa
  call scopy(n, v, 1, x, incx)       ! x
  call scopy(n, va, 1, xa, incx)       ! xa
  !
  ! RMD operation
  call stpsv_rmd(uplo, trans, diag, n, AP, x, incx, APa, xa, dummy, '11')
  !  
  ! Copy result to output variables
  call scopy(sizeAP, APa, 1, ua, 1) ! APa to ua
  call scopy(     n, xa, incx, ua(sizeAP+1), 1) ! xa to ua(sizeAP+1)
  !  
END SUBROUTINE F_RMD
