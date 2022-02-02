program kernel_ex
! Example code calling dkcf
  use block_hybrid_Cholesky
  implicit none
  integer, parameter :: wp = kind(1.0d0)
  integer :: i, info, n
  real(wp), allocatable :: a(:,:), b(:)
 
! Read the matrix order and allocate the arrays
  read(*,*)n
  allocate(a(n,n),b(n))

! Read the rows of the upper-triangular part of the matrix
  do i = 1,n
    read(*,*) a(i,i:n)
  end do

! Factorize the matrix
  call dkcf(n, a, n, info)
  if (info /= 0)then
    write(*,*) "Stopping after failure with info = ", info
    stop 
  end if

! Read the right-hand side and solve the equation
  read(*,*) b(1:n)
  call dtrsv('u','t','n',n,a,n,b,1)
  call dtrsv('u','n','n',n,a,n,b,1)
  write(*,'(8f10.3)') b(1:n)
end program kernel_ex
