PROGRAM TEST_SGEMM_RMD
  use rmd_stesttools
  implicit none
  procedure(F_interf)     :: F
  procedure(F_rmd_interf) :: F_rmd
  real :: tol
  integer :: m, n, k, rep, at, bt, esa, esb, esc, dimu
  character trans(2)
  type(fdata) :: testcase
  !
  tol = rmd_stolerance
  rep = 1
  trans = ['n', 't']
  !
  do m = 1, 5 ! Dimensions of A, B, and C
    do n = 1, 5
      do k = 1, 5
        do at = 1, 2 ! Transposition of A and B
          do bt = 1, 2
            do esa = 0, 3, 3 ! Extra size of leading dimensions of A, B, and C
              do esb = 0, 3, 3
                do esc = 0, 3, 3
                  dimu = m*k + k*n + m*n + 2
                  testcase = fdata([real::], [m, n, k, esa, esb, esc, dimu], [trans(at), trans(bt)])
                  call rmd_stestrandom(F, F_rmd, dimu, m*n, tol, rep, testcase)
                end do
              end do
            end do
          end do
        end do
      end do
    end do
  end do
END PROGRAM TEST_SGEMM_RMD

SUBROUTINE F(testcase, u, v)
  use rmd_stesttools
  type(fdata), intent(inout) :: testcase
  real, intent(in) :: u(*)
  real, intent(out) :: v(*)
  !
  real, allocatable :: A(:, :), B(:, :), C(:, :)
  integer :: m, n, k, ma, na, mb, nb, esa, esb, esc, dimu
  real :: alpha, beta
  character transa, transb
  !
  m = testcase % integers(1)
  n = testcase % integers(2)
  k = testcase % integers(3)
  transa = testcase % characters(1)
  transb = testcase % characters(2)
  esa = testcase % integers(4)
  esb = testcase % integers(5)
  esc = testcase % integers(6)
  dimu = testcase % integers(7)
  alpha = u(dimu-1)
  beta = u(dimu)
  !  
  if(transa == 'n' .or. transa == 'N') then
    ma = m
    na = k
  else
    ma = k
    na = m
  end if
  !
  if(transb == 'n' .or. transb == 'N') then
    mb = k
    nb = n
  else
    mb = n
    nb = k
  end if
  !
  allocate(A(ma+esa, na), B(mb+esb, nb), C(m+esc, n))
  !
  A(1:ma, 1:na) = reshape(u(1:ma*na), [ma, na])
  B(1:mb, 1:nb) = reshape(u(ma*na+1:ma*na+mb*nb), [mb, nb])
  C(1:m, 1:n) = reshape(u(ma*na+mb*nb+1:ma*na+mb*nb+m*n), [m, n])
  !
  call sgemm(transa, transb, m, n, k, alpha, A, ma+esa, B, mb+esb, beta, C, m+esc)
  !
  v(1:m*n) = reshape(C(1:m, 1:n), [m*n])
  !
END SUBROUTINE F

SUBROUTINE F_RMD(testcase, u, v, ua, va)
  use rmd_stesttools
  type(fdata), intent(inout) :: testcase
  real, intent(in) :: u(*)
  real, intent(inout) :: ua(*)
  real, intent(in) :: v(*), va(*)
  !
  real, allocatable :: A(:, :), B(:, :), Aa(:, :), Ba(:, :), Ca(:, :), C0(:, :)
  ! Note: C is not used.
  integer :: m, n, k, ma, na, mb, nb, esa, esb, esc, dimu, lda, ldb, ldc
  real :: alpha, beta, alphaa, betaa
  character transa, transb
  !
  m = testcase % integers(1)
  n = testcase % integers(2)
  k = testcase % integers(3)
  transa = testcase % characters(1)
  transb = testcase % characters(2)
  esa = testcase % integers(4)
  esb = testcase % integers(5)
  esc = testcase % integers(6)
  dimu = testcase % integers(7)
  alpha = u(dimu-1)
  beta = u(dimu)
  alphaa = ua(dimu-1)
  betaa = ua(dimu)
  !
  if(transa == 'n' .or. transa == 'N') then
    ma = m
    na = k
  else
    ma = k
    na = m
  end if
  !
  if(transb == 'n' .or. transb == 'N') then
    mb = k
    nb = n
  else
    mb = n
    nb = k
  end if
  !
  lda = ma + esa
  ldb = mb + esb
  ldc = m + esc
  allocate(A(ma+esa, na), B(mb+esb, nb))
  allocate(Aa(ma+esa, na), Ba(mb+esb, nb), Ca(m+esc, n), C0(m + esc, n))
  !
  ! Initialise matrices to NaN to hopefully catch indexing errors.
  ! (Is testcase redundant?)
  !  
  A(1:ma, 1:na) = reshape(u(1:ma*na), [ma, na])
  B(1:mb, 1:nb) = reshape(u(ma*na+1:ma*na+mb*nb), [mb, nb])
  C0(1:m, 1:n)  = reshape(u(ma*na+mb*nb+1:ma*na+mb*nb+m*n), [m, n])
  !  
  Aa(1:ma, 1:na) = reshape(ua(1:ma*na), [ma, na])
  Ba(1:mb, 1:nb) = reshape(ua(ma*na+1:ma*na+mb*nb), [mb, nb])
  Ca(1:m, 1:n)   = reshape(va(1:m*n), [m, n])
  !  
  call sgemm_rmds(transa, transb, m, n, k, A, lda, B, ldb, C0, ldc, alphaa,&
    & betaa, Ca, '11')
  call sgemm_rmd(transa, transb, m, n, k, alpha, A, lda, B, ldb, &
      beta, ldc, Aa, Ba, Ca, '111')
      !  
  ua(1:ma*na) = reshape(Aa(1:ma, 1:na), [ma*na])
  ua(ma*na+1:ma*na+mb*nb) = reshape(Ba(1:mb, 1:nb), [mb*nb])
  ua(ma*na+mb*nb+1:ma*na+mb*nb+m*n) = reshape(Ca(1:m, 1:n), [m*n])
  ua(dimu-1) = alphaa
  ua(dimu) = betaa
  !
END SUBROUTINE F_RMD
