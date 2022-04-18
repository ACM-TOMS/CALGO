PROGRAM TEST_SSYMM_RMD
  use rmd_stesttools
  implicit none
  procedure(F_interf)     :: F
  procedure(F_rmd_interf) :: F_rmd
  real :: tol
  character :: uplo(2)
  integer :: m, n, ul, rep, esa, esb, esc, dimu
  type(fdata) :: testcase
  !  
  tol = rmd_stolerance
  rep = 1
  uplo = ['u', 'l']
  !
  do m = 1, 5
    do n = 1, 5
      do ul = 1, 2
        do esa = 0, 3, 3
          do esb = 0, 3, 3
            do esc = 0, 3, 3
              dimu = m*n*2+m*(m+1)/2 + 2
              ! side == 'l': A is mxm
              testcase = fdata([real::], [m, n, esa, esb, esc, dimu], ['l', uplo(ul)])
              call rmd_stestrandom(F, F_rmd, dimu, m*n, tol, rep, testcase)
              ! side == 'r': A is nxn
              dimu = m*n*2+(n*(n+1))/2 + 2
              testcase = fdata([real::], [m, n, esa, esb, esc, dimu], ['r', uplo(ul)])
              call rmd_stestrandom(F, F_rmd, dimu, m*n, tol, rep, testcase)
            end do
          end do
        end do
      end do
    end do
  end do
  !  
END PROGRAM TEST_SSYMM_RMD

SUBROUTINE F(testcase, u, v)
  use rmd_stesttools
  implicit none
  ! Arguments
  type(fdata), intent(inout) :: testcase
  real, intent(in)        :: u(*) ! u only contains the variable part of A
  real, intent(out)       :: v(*)
  !
  ! Local variables
  real, allocatable       :: A(:, :), B(:, :), C(:, :)
  integer :: m, n, ma, la, esa, esb, esc, lda, ldb, ldc, dimu
  real :: alpha, beta
  character side, uplo
  !
  m = testcase % integers(1)
  n = testcase % integers(2)
  side = testcase % characters(1)
  uplo = testcase % characters(2)
  esa = testcase % integers(3)
  esb = testcase % integers(4)
  esc = testcase % integers(5)
  dimu = testcase % integers(6)
  !
  if(side == 'l' .or. side == 'L') then
    ma = m
  else
    ma = n
  end if
  la = ma*(ma+1)/2
  lda = ma + esa
  ldb = m + esb
  ldc = m + esc
  !
  allocate (A(lda, ma), B(ldb, n), C(ldc, n))
  !
  ! Reshape the first part of u into A, one column at a time.
  call rmd_strishape_v2m(uplo, ma, u, A, lda)
  !  
  B(1:m, 1:n) = reshape(u(la+1:la+m*n), [m, n])
  C(1:m, 1:n) = reshape(u(la+m*n+1:la+2*m*n), [m, n])
  alpha = u(dimu-1)
  beta = u(dimu)
  !
  call ssymm(side, uplo, m, n, alpha, A, lda, B, ldb, beta, C, ldc)
  !
  ! Shape C into v.
  v(1:m*n) = reshape(C(1:m, 1:n), [m*n])
  !
END SUBROUTINE F

SUBROUTINE F_RMD(testcase, u, v, ua, va)
  use rmd_stesttools
  !  
  implicit none
  type(fdata), intent(inout) :: testcase
  real, intent(in)    :: u(*), v(*), va(*)
  real, intent(inout) :: ua(*)
  !  
  real, allocatable :: A(:, :), B(:, :), Aa(:, :), Ba(:, :), Ca(:, :), C0(:, :)
  integer :: m, n, ma, la, esa, esb, esc, lda, ldb, ldc, dimu
  real :: alpha, beta, alphaa, betaa
  character side, uplo
  !
  m = testcase % integers(1)
  n = testcase % integers(2)
  side = testcase % characters(1)
  uplo = testcase % characters(2)
  esa = testcase % integers(3)
  esb = testcase % integers(4)
  esc = testcase % integers(5)
  dimu = testcase % integers(6)
  !
  if(side == 'l' .or. side == 'L') then
    ma = m
  else
    ma = n
  end if
  la = ma*(ma+1)/2
  lda = ma + esa
  ldb = m + esb
  ldc = m + esc
  !
  allocate (A(lda, ma), B(ldb, n))
  allocate (Aa(lda, ma), Ba(ldb, n), Ca(ldc, n), C0(ldc, n))
  !  
  ! Reshape the first part of u into A, one column at a time.
  call rmd_strishape_v2m(uplo, ma, u, A, lda)
  !
  ! Reshape the first part of ua into Aa, one column at a time.
  call rmd_strishape_v2m(uplo, ma, ua, Aa, lda)
  !  
  B(1:m, 1:n) = reshape(u(la+1:la+m*n), [m, n])
  Ba(1:m, 1:n) = reshape(ua(la+1:la+m*n), [m, n])
  C0(1:m, 1:n) = reshape(u(la+m*n+1:la+2*m*n), [m, n])
  Ca(1:m, 1:n) = reshape(va(1:m*n), [m, n])
  alpha = u(dimu-1)
  beta = u(dimu)
  alphaa = ua(dimu-1)
  betaa = ua(dimu)
  !
  call ssymm_rmds(side, uplo, m, n, A, lda, B, ldb, C0, ldc, alphaa, betaa, Ca, '11')
  call ssymm_rmd(side, uplo, m, n, alpha, A, lda, B, ldb, beta, ldc, &
    Aa, Ba, Ca, '111')
  !
  ! Reshape Aa into the first part of ua
  call rmd_strishape_m2v(uplo, ma, ua, Aa, lda)
  !
  ua(la+1:la+m*n) = reshape(Ba(1:m, 1:n), [m*n])
  ua(la+m*n+1:la+m*n*2) = reshape(Ca(1:m, 1:n), [m*n])
  ua(dimu-1) = alphaa
  ua(dimu) = betaa
  !
END SUBROUTINE F_RMD
