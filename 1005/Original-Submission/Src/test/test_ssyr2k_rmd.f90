PROGRAM TEST_SSYR2K_RMD
  use rmd_stesttools
  implicit none
  procedure(F_interf)     :: F
  procedure(F_rmd_interf) :: F_rmd
  real :: tol
  character :: uplo(2), trans(2)
  integer :: n, k, rep, ul, tr, esa, esb, esc, dimu, dimv
  type(fdata) :: testcase
  !
  tol = rmd_stolerance
  n = 1
  k = 2
  dimu = n*(n+1)/2 + n*k*2 + 2
  dimv = n*(n+1)/2
  testcase = fdata([real::], [n, k, 1, 0, 0, dimu], ['L', 'N'])
  call rmd_stest123(F, F_rmd, dimu, dimv, tol, testcase)
  rep = 1
  uplo = ['u', 'l']
  trans = ['n', 't']
  !
  do n = 1, 5
    do k = 1, 5
      do ul = 1, 2
        do tr = 1, 2
          do esa = 0, 3, 3
            do esb = 0, 3, 3
              do esc = 0, 3, 3
                dimu = n*k*2 + n*(n+1)/2 + 2
                dimv = n*(n+1)/2
                testcase = fdata([real::], [n, k, esa, esb, esc, dimu], [uplo(ul), trans(tr)])
                call rmd_stestrandom(F, F_rmd, dimu, dimv, tol, rep, testcase)
              end do
            end do
          end do
        end do
      end do
    end do
  end do
  !
END PROGRAM TEST_SSYR2K_RMD

SUBROUTINE F(testcase, u, v)
  use rmd_stesttools
  implicit none
  type(fdata), intent(inout) :: testcase
  real, intent(in) :: u(*)
  real, intent(out) :: v(*)
  !
  real, allocatable :: A(:, :), B(:, :), C(:, :)
  integer :: n, k, ma, na, esa, esb, esc, lda, ldb, ldc, dimu
  real :: alpha, beta
  character uplo, trans
  !
  uplo  = testcase % characters(1)
  trans = testcase % characters(2)
  n     = testcase % integers(1)
  k     = testcase % integers(2)
  esa   = testcase % integers(3)
  esb   = testcase % integers(4)
  esc   = testcase % integers(5)
  dimu  = testcase % integers(6)
  alpha = u(dimu-1)
  beta  = u(dimu)
  !
  if(trans == 'n' .or. trans == 'N') then
    ma = n
    na = k
  else
    ma = k
    na = n
  end if
  lda = esa + ma
  ldb = esb + ma
  ldc = esc + n
  !
  allocate(A(lda, na), B(ldb, na), C(ldc, n))
  !
  A(1:ma, 1:na) = reshape(u(1:ma*na), [ma, na])
  B(1:ma, 1:na) = reshape(u(ma*na+1:ma*na*2), [ma, na])
  call rmd_strishape_v2m(uplo, n, u(ma*na*2+1:ma*na*2+(n*(n+1))/2), C, ldc)
  !
  call ssyr2k(uplo, trans, n, k, alpha, A, lda, B, ldb, beta, C, ldc)
  !  
  call rmd_strishape_m2v(uplo, n, v, C, ldc)
  !  
END SUBROUTINE F

SUBROUTINE F_RMD(testcase, u, v, ua, va)
  use rmd_stesttools
  implicit none
  type(fdata), intent(inout) :: testcase
  real, intent(in) :: u(*)
  real, intent(inout) :: ua(*)
  real, intent(in) :: v(*), va(*)
  !
  real, allocatable :: A(:, :), B(:, :), Aa(:, :), Ba(:, :), Ca(:, :), C0(:, :)
  integer :: n, k, ma, na, esa, esb, esc, lda, ldb, ldc, dimu
  real :: alpha, beta, alphaa, betaa
  character uplo, trans
  !  
  n = testcase % integers(1)
  k = testcase % integers(2)
  uplo = testcase % characters(1)
  trans = testcase % characters(2)
  esa = testcase % integers(3)
  esb = testcase % integers(4)
  esc = testcase % integers(5)
  !
  if(trans == 'n' .or. trans == 'N') then
    ma = n
    na = k
  else
    ma = k
    na = n
  end if
  lda = ma + esa
  ldb = ma + esb
  ldc = n + esc
  !
  allocate(A(lda, na), B(ldb, na))
  allocate(Aa(lda, na), Ba(ldb, na), Ca(ldc, n), C0(ldc, n))
  !
  A(1:ma, 1:na) = reshape(u(1:ma*na), [ma, na])
  Aa(1:ma, 1:na) = reshape(ua(1:ma*na), [ma, na])
  B(1:ma, 1:na) = reshape(u(ma*na+1:ma*na*2), [ma, na])
  Ba(1:ma, 1:na) = reshape(ua(ma*na+1:ma*na*2), [ma, na])
  call rmd_strishape_v2m(uplo, n, va, Ca, ldc)
  call rmd_strishape_v2m(uplo, n, u(ma*na*2+1:ma*na*2+(n*(n+1))/2), C0, ldc)
  dimu   = testcase % integers(6)
  alpha = u(dimu-1)
  beta = u(dimu)
  alphaa = ua(dimu-1)
  betaa = ua(dimu)
  !
  call ssyr2k_rmds(uplo, trans, n, k, A, lda, B, ldb, C0, ldc, alphaa, betaa,&
    & Ca, '11')
  call ssyr2k_rmd(uplo, trans, n, k, alpha, A, lda, B, ldb, beta, ldc, &
    Aa, Ba, Ca, '111')
  !
  ua(1:ma*na) = reshape(Aa(1:ma, 1:na), [ma*na])
  ua(ma*na+1:ma*na*2) = reshape(Ba(1:ma, 1:na), [ma*na])
  call rmd_strishape_m2v(uplo, n, ua(ma*na*2+1:ma*na*2+n*(n+1)/2), Ca, ldc)
  ua(dimu-1) = alphaa
  ua(dimu) = betaa
  !  
END SUBROUTINE F_RMD
