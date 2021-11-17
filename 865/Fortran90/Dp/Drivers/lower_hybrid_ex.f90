program lower_hybrid_ex

! Example code to illustrate factorization in blocked hybrid format and
! solution of a set of equations using the factor in blocked hybrid format.

! The code reads a matrix in lower packed format, converts to
! lower blocked hybrid format, performs Cholesky factorization, and
! solves a set of equations. 

  use block_hybrid_Cholesky
  implicit none
  integer, parameter :: wp = kind(1.0d0)
  integer :: info, n
  real(wp), allocatable :: ap(:), b(:)

! Read the matrix order
  read(*,*) n

! Allocate the arrays
  allocate(ap(n*(n+1)/2), b(n))

! Read the lower-triangular matrix in the lower packed format
  read(*,*) ap(1:n*(n+1)/2)

! Transform the matrix to lower blocked hybrid format
  call PpHpp('l', n, ap, info )
  if (info /= 0) call terminate("PpHpp")

! Factorize the matrix
  call HppTf( 'l', n, ap, info )
  if (info /= 0) call terminate("HppTf")

! Read the right-hand side and solve the equation
  read(*,*) b(1:n)
  call HppTs1( 'l', n, ap, b, info )
  if (info /= 0) call terminate("HppTs1")
  write(*,'(8f10.3)') b(1:n)

contains

  subroutine terminate(name)
    character(*) name
    write(*,*) "Stopping after failure in ",name," with info = ", info
    stop
  end subroutine terminate

end program lower_hybrid_ex
