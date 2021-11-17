module minresqlpMexModule
  use  minresqlpDataModule,    only : dp, sp, ip

  real(dp), allocatable, public :: A(:,:)
  real(dp), allocatable, public :: M(:,:)
contains

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine allocMem(n)
    integer(ip), intent(in)  :: n
    allocate( A(n,n) )
    allocate( M(n,n) )
  end subroutine allocMem
    
  subroutine freeMem(n)
    deallocate( A )
    deallocate( M )
  end subroutine freeMem

  subroutine Aprod (n,x,y)
    use  minresqlpDataModule,    only : dp, sp, ip
    integer(ip), intent(in)  :: n
    real(dp),    intent(in)  :: x(n)
    real(dp),    intent(out) :: y(n)

    !-------------------------------------------------------------------
    ! Aprod  computes y = A*x for some matrix A.
    ! This is a simple example for testing MINRESQLP.
    !-------------------------------------------------------------------

    integer(ip) :: i, j

    do i = 1, n
       y(i) = 0
       do j = 1,n
           y(i) = y(i) + A(i,j) * x(j)
       end do
    end do

  end subroutine Aprod
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine Msolve(n,x,y)
    use  minresqlpDataModule,    only : dp, sp, ip      
    integer(ip), intent(in)  :: n
    real(dp),    intent(in)  :: x(n)
    real(dp),    intent(out) :: y(n)

    !-------------------------------------------------------------------
    ! Msolve solves M*y = x for some symmetric positive-definite matrix M.
    ! This is a simple example for testing MINRESQLP.
    ! Ashift will be the same as shift in MINRESQLP.
    !
    ! If Mpert = 0, the preconditioner will be exact, so
    ! MINRESQLP should require either one or two iterations,
    ! depending on whether (A - shift*I) is positive definite or not.
    !
    ! If Mpert is nonzero, somewhat more iterations will be required.
    !-------------------------------------------------------------------

    do i = 1, n
       y(i) = 0
       do j = 1,n
           y(i) = y(i) + M(i,j) * x(i)
       end do
    end do

  end subroutine Msolve
end module minresqlpMexModule
