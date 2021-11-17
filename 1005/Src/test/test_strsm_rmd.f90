PROGRAM TEST_STRSM_RMD
  use rmd_stesttools
  implicit none
  procedure(F_interf)     :: F
  procedure(F_rmd_interf) :: F_rmd
  !
  type(fdata) :: testcase
  integer     :: m, n, k, dimu, dimv, s, ul, t, d, esa, esb
  real        :: tol = rmd_stolerance
  character   :: side(2), uplo(2), transA(2), diag(2)
  real, allocatable :: u(:)
  !
  side   = ['L', 'R']
  uplo   = ['U', 'L']
  transA = ['N', 'T']
  diag   = ['N', 'U']
  !
  do m = 1, 5
    do n = 1, 5
      do s = 1, 2
        do ul = 1, 2
          do t = 1, 2
            do d = 1, 2
              do esa = 0, 3, 3
                do esb = 0, 3, 3
                  !
                  if( s == 1 ) then
                    k = m
                  else
                    k = n
                  end if
                  if (d == 1) dimu  = k*(k+1)/2 + m*n + 1
                  if (d == 2) dimu  = k*(k-1)/2 + m*n + 1
                  dimv = m*n
                  allocate(u(dimu))
                  call random_number(u)
                  u = 2.0*u - 1.0
                  if (d == 1) call fix_diagonal(uplo(ul), u, k)
                  testcase = fdata([real::], [m, n, esa, esb, dimu], &
                    &              [side(s), uplo(ul), transA(t), diag(d)])
                  call rmd_stestf(F, F_rmd, u, dimu, dimv, tol, testcase)
                  deallocate(u)
                end do
              end do
            end do
          end do
        end do
      end do
    end do
  end do
CONTAINS
  subroutine fix_diagonal(uplo, u, k)
    ! Make diag(A) bounded away from zero
    real :: u(*)
    integer :: i, k, jj
    character :: uplo
    jj = 1
    do i=1,k
      if (abs(u(jj)) < 0.5) u(jj) = sign(0.5, u(jj))
      if (uplo == 'u' .or. uplo == 'U') then
        jj = jj + i+1
      else
        jj = jj + k-i+1
      endif
    enddo
  end subroutine fix_diagonal
END PROGRAM TEST_STRSM_RMD

SUBROUTINE F(testcase, u, v)
  use rmd_stesttools
  implicit none
  ! Arguments
  type(fdata), intent(inout) :: testcase
  real, intent(in)  :: u(*)
  real, intent(out) :: v(*)
  !
  ! Local variables
  integer :: m, n, k, lda, ldb, esa, esb, dimu, dimA
  real :: alpha
  character side, uplo, transA, diag
  !  
  real, allocatable :: A(:,:), B(:,:)
  !
  side   = testcase % characters(1)
  uplo   = testcase % characters(2)
  transA = testcase % characters(3)
  diag   = testcase % characters(4)
  m      = testcase % integers(1)
  n      = testcase % integers(2)
  esa    = testcase % integers(3)
  esb    = testcase % integers(4)
  dimu   = testcase % integers(5)
  !
  if((side == 'L').or.(side == 'l')) then
    k = m
  else
    k = n
  end if
  lda = k+esa
  ldb = m+esb
  allocate(A(lda, k), B(ldb, n))
  !
  dimA = k*(k+1)/2
  if (diag == 'N') then
      call rmd_strishape_v2m(uplo, k, u(1:dimA), A, lda) ! u to A
  else
    call rmd_scopy_unit_v2m(uplo, k, u(1:dimA), A, lda)
  endif
  B(1:m, 1:n) = reshape(u(dimu-m*n : dimu-1), [m, n]) ! u to B
  alpha = u(dimu)
  !
  ! BLAS operation
  call strsm(side, uplo, transA, diag, m, n, alpha, A, lda, B, ldb)
  !  
  ! Copy back from B to v
  v(1:m*n) = reshape(B(1:m, 1:n), [m*n])
END SUBROUTINE F

SUBROUTINE F_RMD(testcase, u, v, ua, va)
  use rmd_stesttools
  implicit none
  ! Arguments
  type(fdata), intent(inout) :: testcase
  real, intent(in)    :: u(*), va(*), v(*)
  real, intent(inout) :: ua(*)
  !  
  ! Local variables
  real, allocatable :: A(:, :), B(:, :), Aa(:, :), Ba(:, :), B0(:, :)
  integer :: m, n, k, lda, ldb, esa, esb, dimu
  real :: alpha, dummy(1), alphaa
  character side, uplo, transA, diag
  !
  side   = testcase % characters(1)
  uplo   = testcase % characters(2)
  transA = testcase % characters(3)
  diag   = testcase % characters(4)
  m      = testcase % integers(1)
  n      = testcase % integers(2)
  esa    = testcase % integers(3)
  esb    = testcase % integers(4)
  dimu   = testcase % integers(5)
  !
  if((side == 'L').or.(side == 'l')) then
    k = m
  else
    k = n
  end if
  lda = k+esa
  ldb = m+esb
  !  
  allocate(A(lda, k), B(ldb, n), Aa(lda, k), Ba(ldb, n), B0(ldb, n))
  !
  ! Copy inital values to local variables
  A = 77
  Aa = 77
  if (diag == 'N') then
    call rmd_strishape_v2m(uplo, k, u, A, lda) ! A
    call rmd_strishape_v2m(uplo, k, ua, Aa, lda) ! Aa
  else
    call rmd_scopy_unit_v2m(uplo, k, u, A, lda)
    call rmd_scopy_unit_v2m(uplo, k, ua, Aa, lda)
  endif
  B0(1:m, 1:n) = reshape(u(dimu-m*n : dimu-1), [m, n]) ! u to B0
  B(1:m, 1:n)  = reshape(v(1:m*n), [m, n]) ! B
  Ba(1:m, 1:n) = reshape(va(1:m*n), [m, n]) ! Ba
  alpha = u(dimu)
  alphaa = ua(dimu)
  !
  ! RMD operation
  call strsm_rmds(side, uplo, transA, diag, m, n, A, lda, B0, ldb, alphaa, Ba)
  call strsm_rmd(side, uplo, transA, diag, m, n, alpha, A, lda, B, ldb, Aa, Ba, dummy, '11')
  !  
  ! Copy result to output variables
  if (diag == 'N') then
    call rmd_strishape_m2v(uplo, k, ua, Aa, lda) ! Aa to ua
  else
    call rmd_scopy_unit_m2v(uplo, k, ua, Aa, lda) ! Aa to ua
  endif
  ua(dimu-m*n : dimu-1) = reshape(Ba(1:m, 1:n), [m*n]) ! Ba to ua
  ua(dimu) = alphaa
  !  
END SUBROUTINE F_RMD
