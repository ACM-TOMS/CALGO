!*********************************************************************************************************************************
!*
!* MODULE: DOMAIN
!*
!* PURPOSE: builds the domain
!*
!* DEPENDENCIES:
!*               - PRECISION
!*
!* GLOBALS:
!*          - LOW (integer): number of bytes for simple precision variables (from PRECISION)
!*          - HIGH (integer): number of bytes for double precision variables (from PRECISION)
!*          - PLO (integer): precision for logical variables (from PRECISION)
!*          - PCH (integer): precision for character variables (from PRECISION)
!*          - PIN (integer): precision for integer variables (from PRECISION)
!*          - PRE (integer): precision for real variables (from PRECISION)
!*
!*          - NT (integer): number of time steps
!*          - A (real): initial time
!*          - DELTAT (real): time step
!*          - T (real array): vector of time
!*
!* SCOPE: this module belongs to the example
!*
!* AUTHOR: Germán Ariel Torres
!*         Facultad de Matemática, Astronomía y Física
!*         Universidad Nacional de Córdoba
!*         Argentina
!*
!* EMAIL: torres@famaf.unc.edu.ar
!*
!*********************************************************************************************************************************
!*
module domain
  !*
  use precision
  !*
  implicit none
  !*
  integer( kind = pin ) , public :: nt !* number of time steps
  real( kind = pre ) , public :: a !* initial time
  real( kind = pre ) , public :: deltat !* time step
  real( kind = pre ) , allocatable , dimension( : ) , public :: t !* vector of time
  !*
contains
  !*
  !* set the domain parameters
  subroutine domain_parametersup()
    implicit none
    nt = 250_pin
    a = 0.0_pre
    deltat = 0.1_pre
  end subroutine domain_parametersup
  !*
  !* set the time domain
  subroutine domain_domain()
    implicit none
    integer( kind = pin ) :: l
    allocate( t( nt ) )
    do l = 1_pin , nt
       t( l ) = a + ( l - 1_pin ) * deltat
    end do
  end subroutine domain_domain
  !*
  !* free allocated memory
  subroutine domain_destructor()
    implicit none
    if ( allocated( t ) ) deallocate( t )
  end subroutine domain_destructor
  !*
end module domain

