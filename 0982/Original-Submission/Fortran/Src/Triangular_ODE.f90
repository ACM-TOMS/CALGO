module Triangular_ODE

  implicit NONE
  public

  interface Compute_Z_Matrix
    module procedure Compute_Z_Matrix_D, Compute_Z_Matrix_S
  end interface

  interface Compute_Solution
    module procedure Compute_Solution_D, Compute_Solution_S
  end interface

contains

  subroutine Compute_Z_Matrix_D ( A, Z )
    integer, parameter :: RK = kind(0.0d0)
!   real(rk), intent(in) :: A(:,:)
!   real(rk), intent(inout) :: Z(:,:)
    include "Compute_Z_Matrix.f9h"
  end subroutine Compute_Z_Matrix_D

  subroutine Compute_Z_Matrix_S ( A, Z )
    integer, parameter :: RK = kind(0.0e0)
!   real(rk), intent(in) :: A(:,:)
!   real(rk), intent(inout) :: Z(:,:)
    include "Compute_Z_Matrix.f9h"
  end subroutine Compute_Z_Matrix_S

  subroutine Compute_Solution_D ( T, A, Z, N0, Amount )
    integer, parameter :: RK = kind(0.0d0)
!   real(rk), intent(in) :: T
!   real(rk), intent(in) :: A(:,:) ! Only the diagonal is needed
!   real(rk), intent(in) :: Z(:,:)
!   real(rk), intent(in) :: N0(:)
!   real(rk), intent(out) :: Amount(:)
    include "Compute_Solution.f9h"
  end subroutine Compute_Solution_D

  subroutine Compute_Solution_S ( T, A, Z, N0, Amount )
    integer, parameter :: RK = kind(0.0e0)
!   real(rk), intent(in) :: T
!   real(rk), intent(in) :: A(:,:) ! Only the diagonal is needed
!   real(rk), intent(in) :: Z(:,:)
!   real(rk), intent(in) :: N0(:)
!   real(rk), intent(out) :: Amount(:)
    include "Compute_Solution.f9h"
  end subroutine Compute_Solution_S

end module Triangular_ODE
