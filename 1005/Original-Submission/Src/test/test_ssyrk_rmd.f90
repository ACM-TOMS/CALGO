PROGRAM TEST_SSYRK_RMD
  use rmd_stesttools
  implicit none
  procedure(F_interf)     :: F
  procedure(F_rmd_interf) :: F_rmd
  real :: tol
  character :: uplo(2), trans(2)
  integer :: n, k, rep, ul, tr, esa, esc, dimu
  type(fdata) :: testcase
  !
  tol = rmd_stolerance
  n = 1
  k = 1
  dimu = n*(n+1)/2 + n*k + 2
  testcase = fdata([real::], [n, k, 0, 0, dimu], ['L', 'N'])
  call rmd_stest123(F, F_rmd, dimu, n*(n+1)/2, tol, testcase)
  rep = 1
  uplo = ['u', 'l']
  trans = ['n', 't']
  !
  do n = 1, 5
    do k = 1, 5
      do ul = 1, 2
        do tr = 1, 2
          do esa = 0, 3, 3
            do esc = 0, 3, 3
              dimu = n*(n+1)/2+n*k + 2
              testcase = fdata([real::], [n, k, esa, esc, dimu], [uplo(ul), trans(tr)])
              call rmd_stestrandom(F, F_rmd, dimu, n*(n+1)/2, tol, rep, testcase)
            end do
          end do
        end do
      end do
    end do
  end do
  !
END PROGRAM TEST_SSYRK_RMD

SUBROUTINE F(testcase, u, v)
  use rmd_stesttools
  implicit none
  type(fdata), intent(inout) :: testcase
  real, intent(in) :: u(*)
  real, intent(out) :: v(*)
  !
  real, allocatable :: A(:, :), C(:, :)
  integer :: n, k, ma, na, esa, esc, dimu
  real :: alpha, beta
  character uplo, trans
  !
  n = testcase % integers(1)
  k = testcase % integers(2)
  uplo = testcase % characters(1)
  trans = testcase % characters(2)
  esa = testcase % integers(3)
  esc = testcase % integers(4)
  dimu = testcase % integers(5)
  alpha = u(dimu-1)
  beta = u(dimu)
  !
  if(trans == 'n' .or. trans == 'N') then
    ma = n
    na = k
  else
    ma = k
    na = n
  end if
  !
  allocate(A(ma+esa, na), C(n+esc, n))
  !
  A(1:ma, 1:na) = reshape(u(1:ma*na), [ma, na])
  call rmd_strishape_v2m(uplo, n, u(ma*na+1:ma*na+(n*(n+1))/2), C, n+esc)
  !
  call ssyrk(uplo, trans, n, k, alpha, A, ma+esa, beta, C, n+esc)
  !
  call rmd_strishape_m2v(uplo, n, v, C, n+esc)
  !
END SUBROUTINE F

SUBROUTINE F_RMD(testcase, u, v, ua, va)
  use rmd_stesttools
  !  
  implicit none
  type(fdata), intent(inout) :: testcase
  real, intent(in) :: u(*)
  real, intent(inout) :: ua(*)
  real, intent(in) :: v(*), va(*)
  !
  real, allocatable :: A(:, :), Aa(:, :), Ca(:, :), C0(:, :)
  integer :: n, k, ma, na, esa, esc, dimu, lda, ldc
  real :: alpha, beta, alphaa, betaa
  character uplo, trans
  !
  n = testcase % integers(1)
  k = testcase % integers(2)
  uplo = testcase % characters(1)
  trans = testcase % characters(2)
  esa = testcase % integers(3)
  esc = testcase % integers(4)
  dimu = testcase % integers(5)
  !  
  if(trans == 'n' .or. trans == 'N') then
    ma = n
    na = k
  else
    ma = k
    na = n
  end if
  lda = ma + esa
  ldc = n + esc
  !
  allocate(A(lda, na))
  allocate(Aa(lda, na), Ca(ldc, n), C0(ldc, n))
  !  
  A(1:ma, 1:na) = reshape(u(1:ma*na), [ma, na])
  Aa(1:ma, 1:na) = reshape(ua(1:ma*na), [ma, na])
  call rmd_strishape_v2m(uplo, n, va, Ca, ldc)
  call rmd_strishape_v2m(uplo, n, u(ma*na+1:ma*na+(n*(n+1))/2), C0, n+esc)
  !
  dimu   = testcase % integers(5)
  alpha = u(dimu-1)
  beta = u(dimu)
  alphaa = ua(dimu-1)
  betaa = ua(dimu)

  !
  call ssyrk_rmds(uplo, trans, n, k, A, lda, C0, ldc, alphaa, betaa, Ca, '11')
  call ssyrk_rmd(uplo, trans, n, k, alpha, A, lda, beta, ldc, Aa, Ca, '11')
  !
  ua(1:ma*na) = reshape(Aa(1:ma, 1:na), [ma*na])
  call rmd_strishape_m2v(uplo, n, ua(ma*na+1:ma*na+n*(n+1)/2), Ca, ldc)
  ua(dimu-1) = alphaa
  ua(dimu) = betaa
  !
END SUBROUTINE F_RMD
