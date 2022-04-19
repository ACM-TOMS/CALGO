SUBROUTINE F(testcase, u, v)
  use rmd_stesttools
  type(fdata), intent(inout) :: testcase
  real, intent(in) :: u(*)
  real, intent(out) :: v(*)
  !
  real, allocatable :: A(:, :)
  integer n, esa, lda, info
  character uplo
  !
  n = testcase % integers(1)
  uplo = testcase % characters(1)
  esa = testcase % integers(2)
  lda = n + esa
  !
  allocate(A(lda, n))
  call rmd_strishape_v2m(uplo, n, u, A, lda)
  !
  call spotrf(uplo, n, A, lda, info)
  if( info /= 0 ) print *, "Cholesky factorization impossible", info
  !
  call rmd_strishape_m2v(uplo, n, v, A, lda)
END SUBROUTINE F

SUBROUTINE F_RMD(testcase, u, v, ua, va)
  use rmd_stesttools
  type(fdata), intent(inout) :: testcase
  real, intent(in) :: u(*)
  real, intent(inout) :: ua(*)
  real, intent(in) :: v(*), va(*)
  !
  real, allocatable :: A(:, :), Aa(:, :) 
  integer n, esa, lda
  character uplo
  !
  n = testcase % integers(1)
  uplo = testcase % characters(1)
  esa = testcase % integers(2)
  lda = n + esa
  !
  allocate(A(lda, n), Aa(lda, n))
  !
  call rmd_strishape_v2m(uplo, n, v, A, lda)
  call rmd_strishape_v2m(uplo, n, va, Aa, lda)
  !
  call spotrf_rmd(uplo, n, A, lda, Aa)
  !
  call rmd_strishape_m2v(uplo, n, ua, Aa, lda)

END SUBROUTINE F_RMD

PROGRAM TEST_SPOTRF_RMD
  use rmd_stesttools
  implicit none
  real :: tol = rmd_stolerance
  integer :: n, j, esc, iuplo, info
  type(fdata) :: testcase
  procedure(F_interf) :: F
  procedure(F_rmd_interf) :: F_rmd
  real, allocatable :: A(:, :), Aa(:, :), u(:)
  character uplo(2), ul
  !
  uplo = ['L', 'U']
  ! First use the same call as the Matlab program test_chol_rmd uses first
  ! (for use in debugging)
  call hilb(3, A)
  call lehmer(3, Aa)
  call spotrf('l', 3, A, 3, info)
  call spotrf_rmd('l', 3, A, 3, Aa)
  call spotrf('u', 3, A, 3, info)
  call spotrf_rmd('u', 3, A, 3, Aa)
  ! Test with Hilbert matrices + I (known to be pos.def.)
  do n = 3, 10
    do esc = 0, 3, 3
      allocate(u(n*(n+1)/2))
      call hilb(n, A)
      do j = 1, n
        A(j, j) = A(j, j) + 1.0
      end do
      do iuplo = 1, 2
        ul = uplo(iuplo)
        call rmd_strishape_m2v(ul, n, u, A, n)
        testcase = fdata([real::], [n, esc], [ul])
        call rmd_stestf(F, F_rmd, u, n*(n+1)/2, n*(n+1)/2, tol, testcase)
      end do
      deallocate(u)
    end do
  end do

CONTAINS
  SUBROUTINE HILB(n, H)
    integer, intent(in) :: n
    real, allocatable, intent(out) :: H(:, :)
    integer i, j
    allocate(H(n, n))
    do j = 1, n
      do i = 1, n
        H(i, j) = 1.0/(i + j - 1)
      enddo
    end do
  END SUBROUTINE HILB

  SUBROUTINE LEHMER(n, H)
    integer, intent(in) :: n
    real, allocatable, intent(out) :: H(:, :)
    integer i, j
    allocate(H(n, n))
    do j = 1, n
      do i = 1, n
        H(i, j) = (min(i, j) + 0.0)/max(i, j)
      enddo
    end do
  END SUBROUTINE LEHMER
END PROGRAM
