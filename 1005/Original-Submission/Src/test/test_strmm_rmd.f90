PROGRAM TEST_STRMM_RMD
  use rmd_stesttools
  implicit none
  procedure(F_interf)     :: F
  procedure(F_rmd_interf) :: F_rmd
  !
  type(fdata) :: testcase
  integer     :: m, n, k, dimu, dimv, s, u, t, d, esa, esb
  !
  real    :: tol = rmd_stolerance
  integer :: rep = 1
  !
  character :: side(2), uplo(2), transA(2), diag(2)
  !
  side   = ['L', 'R']
  uplo   = ['L', 'U']
  transA = ['N', 'T']
  diag   = ['N', 'U']
  m = 1
  n = 1
  k = 1
  dimu = m*n + k*(k+1)/2 + 1
  testcase = fdata([real::], [m, n, 0, 0, dimu], ['L', 'L', 'N', 'N'])
  call rmd_stest123(F, F_rmd, dimu, m*n, tol, testcase)
  !
  do m = 1, 5
    do n = 1, 5
      do s = 1, 2
        do u = 1, 2
          do t = 1, 2
            do d = 1, 2
              do esa = 0, 3, 3
                do esb = 0, 3, 3
                  if( s == 1 ) then
                    k = m
                  else
                    k = n
                  end if
                  !
                  dimu  = m*n + k*(k+1)/2 + 1
                  dimv = m*n
                  testcase = fdata([real::], [m, n, esa, esb, dimu], [side(s), uplo(u), transA(t), diag(d)])
                  call rmd_stestrandom(F, F_rmd, dimu, dimv, tol, rep, testcase)
                end do
              end do
            end do
          end do
        end do
      end do
    end do
  end do
END PROGRAM TEST_STRMM_RMD

SUBROUTINE F(testcase, u, v)
  use rmd_stesttools
  implicit none
  ! Arguments
  type(fdata), intent(inout) :: testcase
  real, intent(in)  :: u(*)
  real, intent(out) :: v(*)
  !
  ! Local variables
  integer :: m, n, k, sizeA, esa, esb, lda, ldb, dimu
  real :: alpha
  character side, uplo, transA, diag
  real, allocatable :: A(:, :), B(:, :)
  !
  side   = testcase % characters(1)
  uplo   = testcase % characters(2)
  transA = testcase % characters(3)
  diag   = testcase % characters(4)
  m      = testcase % integers(1)
  n      = testcase % integers(2)
  esa  = testcase % integers(3)
  esb  = testcase % integers(4)
  dimu = testcase % integers(5)
  !
  if((side == 'L').or.(side == 'l')) then
    k = m
  else
    k = n
  end if
  lda = k + esa
  ldb = m + esb
  !
  sizeA = (k*(k+1))/2
  !  
  allocate(A(lda, k), B(ldb, n))
  ! Initialise A, B, and C to NaN to hopefully catch indexing errors.
  !
  call rmd_strishape_v2m(uplo, k, u, A, lda) ! A
  B(1:m, 1:n) = reshape(u(sizeA+1:sizeA+m*n), [m, n])   ! B
  alpha = u(dimu)
  !
  ! BLAS operation
  call strmm(side, uplo, transA, diag, m, n, alpha, A, lda, B, ldb)
  !  
  ! Copy back from B to v
  v(1:m*n) = reshape(B(1:m, 1:n), [m*n])
  !
END SUBROUTINE F

SUBROUTINE F_RMD(testcase, u, v, ua, va)
  use rmd_stesttools
  implicit none
  ! Arguments
  type(fdata), intent(inout) :: testcase
  real, intent(in) :: u(*), v(*), va(*)
  real, intent(inout) :: ua(*)
  !  
  ! Local variables
  real, allocatable :: A(:, :), B0(:, :), Aa(:, :), Ba(:, :), wrk(:)
  integer :: m, n, k, sizeA, esa, esb, lda, ldb, dimu
  real :: alpha, alphaa
  character side, uplo, transA, diag
  !
  side   = testcase % characters(1)
  uplo   = testcase % characters(2)
  transA = testcase % characters(3)
  diag   = testcase % characters(4)
  m      = testcase % integers(1)
  n      = testcase % integers(2)
  esa  = testcase % integers(3)
  esb  = testcase % integers(4)
  dimu = testcase % integers(5)
  !
  if((side == 'L').or.(side == 'l')) then
    k = m
  else
    k = n
  end if
  lda = k + esa
  ldb = m + esb
  !  
  allocate(A(lda, k), B0(ldb, n), Aa(lda, k), Ba(ldb, n), wrk(max(n,m)))
  ! Initialise A, B0, and C to NaN to hopefully catch indexing errors.
  !
  sizeA = k*(k+1)/2
  !
  ! Copy inital values to local variables
  call rmd_strishape_v2m(uplo, k, u, A, lda)
  call rmd_strishape_v2m(uplo, k, ua, Aa, lda)
  B0(1:m, 1:n) = reshape(u(sizeA+1:sizeA+m*n), [m, n])
  Ba(1:m, 1:n) = reshape(va(1:m*n), [m, n])
  alpha = u(dimu)
  alphaa = ua(dimu)
  !
  ! RMD operation
  call strmm_rmds(side, uplo, transa, diag, m, n, A, lda, B0, ldb, alphaa, Ba, wrk)
  call strmm_rmd(side, uplo, transA, diag, m, n, alpha, &
    A, lda, B0, ldb, Aa, Ba, '11')
  !  
  ! Copy result to output variables
  call rmd_strishape_m2v(uplo, k, ua, Aa, lda) ! Aa to ua
  ua(sizeA+1:sizeA+m*n) = reshape(Ba(1:m, 1:n), [m*n]) ! Ba to ua(sizeA+1)
  ua(dimu) = alphaa
  !
END SUBROUTINE F_RMD
