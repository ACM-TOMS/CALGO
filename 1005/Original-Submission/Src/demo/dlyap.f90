! DLYAP  Solve discrete Lyapunov equation
! 
! CALL DLYAP(transp, A, Sig, S) solves the discrete Lyapunov equation
! 
!             S - op(A)*S*op(A)' = Sig,
!             
! for S, where op(A) is A or A' depending on whether transp is 'n' or 't'. It
! uses Gaussian elimination to solve the equation directly, and is adapted from
! vyw_factorize / vyw_solve in [1]. It has the advantage over SLICOT that it is
! not GPL-licensed but it is however considerably slower (4 times for n = 10,
! 100 times for n = 100). Note that it is possible to call the routine with
! Sig=S: CALL DLYAP(transp, A, S, S) to let the solution overwrite the right
! hand side.
!
! [1] K Jonasson: Algorithm 878: Exact VARMA likelihood and its gradient for
!     complete and incomplete data with Matlab, ACM TOMS 2008.

subroutine dlyap(transp, n, A, Sig, S)
  !use dispmodule
  implicit none
  integer, parameter :: dble = kind(0d0)
  character, intent(in)     :: transp   ! work with A or A'
  integer, intent(in)       :: n        ! matrix size
  real(dble), intent(in)    :: A(n,n)   ! n by n, coefficients
  real(dble), intent(inout) :: Sig(n,n) ! n by n, rhs, symmetric, stored in lower triangle
  real(dble), intent(inout) :: S(n,n)   ! n by n, solution, symmetric, both triangles stored
  
  real(dble), allocatable :: F(:,:), y(:)
  integer :: nn, i, j, k, l, i1, i2, j1, j2, info
  integer, allocatable :: ipiv(:)
  nn = n*(n+1)/2
  allocate (F(nn, nn), y(nn), ipiv(nn))
  j1 = 1
  j2 = n
  l = n
  do j=1,n
    i1 = 1
    i2 = n
    k = n
    if (transp == 'n' .or. transp == 'N') then 
      do i=1,n
        F(i1:i2, j1:j2) = - A(i,j)*A(i:n,j:n)
        !- A(i:n,j)*A(i,j:n):
        call dger(k, l, -1d0, A(i,j), 1, A(i,j), n, F(i1, j1), nn)
        k = k - 1
        i1 = i2 + 1
        i2 = i2 + k
      end do
    else
      do i=1,n
        F(i1:i2, j1:j2) = - A(j,i)*transpose(A(j:n,i:n))
        !- A(i:n,j)*A(i,j:n):
        call dger(k, l, -1d0, A(j,i), n, A(j,i), 1, F(i1, j1), nn)
        k = k - 1
        i1 = i2 + 1
        i2 = i2 + k
      end do      
    endif
    F(j1,j1) = F(j1,j1) + 1d0
    y(j1:j2) = Sig(j:n,j)
    l = l - 1
    j1 = j2 + 1
    j2 = j2 + l
  end do
  do i=1,nn
    F(i,i) = F(i,i) + 1d0
  enddo
  call dgesv(nn, 1, F, nn, ipiv, y, nn, info)
  if (info /= 0) stop 'Singular matrix in dlyap'
  j1 = 1
  j2 = n
  do i=1,n ! copy y to lower triangle of S
    S(i:n,i) = y(j1:j2)
    S(i,i:n) = y(j1:j2) ! .. and to upper triangle
    S(i,i) = 2*S(i,i)   ! double diagonal
    j1 = j2 + 1
    j2 = j2 + n - i
  enddo
end subroutine dlyap
