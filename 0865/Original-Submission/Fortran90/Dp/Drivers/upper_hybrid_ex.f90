program upper_hybrid_ex
! Example code to illustrate factorization in blocked hybrid format when
! the factor is required in packed format for compatibility with other code. 

! The code reads a matrix in upper packed format, converts to
! upper blocked hybrid format, performs Cholesky factorization, transforms
! back to upper packed format and solves sets of equations using Lapack. 

  use block_hybrid_Cholesky
  implicit none
  integer, parameter :: wp = kind(1.0d0)
  integer :: i, info, n, nrhs
  real(wp), allocatable :: ap(:), b(:,:)

! Read the matrix order and number of rhs
  read(*,*) n, nrhs

! Allocate the arrays
  allocate(ap(n*(n+1)/2), b(n,nrhs))

! Read the upper-triangular matrix in the upper packed format
  read(*,*) ap(1:n*(n+1)/2)

! Transform the matrix to the upper blocked hybrid format
  call PpHpp('u', n, ap, info )
  if (info /= 0) call terminate("PpHpp")

! Factorize the matrix
  call HppTf( 'u', n, ap, info )
  if (info /= 0) call terminate("HppTf")

! Transform to the upper packed format
  call HppPp('u', n, ap, info )
  if (info /= 0) call terminate("HppPp")

! Read the right-hand sides and solve the equations
  read(*,*) b(1:n,1:nrhs)
  call HppTs( 'u', n, nrhs, ap, b, n, info )
  if (info /= 0) call terminate("dpptrs")
  do i = 1,nrhs
    write(*,'(8f10.3)') b(1:n,i)
  end do

contains

  subroutine terminate(name)
    character(*) name
    write(*,*) "Stopping after failure in ",name," with info = ", info
    stop
  end subroutine terminate

end program upper_hybrid_ex
