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
!*          - NT (integer): number of nodes in t-direction
!*          - NX (integer): number of nodes in x-direction
!*          - NY (integer): number of nodes in y-direction
!*          - AT (real): left side of the interval in t-direction
!*          - AX (real): left side of the interval in x-direction
!*          - AY (real): left side of the interval in y-direction
!*          - BX (real): right side of the interval in x-direction
!*          - BY (real): right side of the interval in y-direction
!*          - DELTAT (real): step in t-direction
!*          - DELTAX (real): step in x-direction
!*          - DELTAY (real): step in y-direction
!*          - T (real array): t-domain
!*          - X (real array): x-domain
!*          - Y (real array): y-domain
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
module domain
  !*
  use precision
  !*
  implicit none
  !*
  integer( kind = pin ) , public :: nt !* number of nodes in t-direction
  integer( kind = pin ) , public :: nx !* number of nodes in x-direction
  integer( kind = pin ) , public :: ny !* number of nodes in y-direction
  real( kind = pre ) , public :: at !* left side of the interval in t-direction
  real( kind = pre ) , public :: ax !* left side of the interval in x-direction
  real( kind = pre ) , public :: ay !* left side of the interval in y-direction
  real( kind = pre ) , public :: bx !* right side of the interval in x-direction
  real( kind = pre ) , public :: by !* right side of the interval in y-direction
  real( kind = pre ) , public :: deltat !* step in t-direction
  real( kind = pre ) , public :: deltax !* step in x-direction
  real( kind = pre ) , public :: deltay !* step in y-direction
  real( kind = pre ) , allocatable , dimension( : ) , public :: t !* t-domain
  real( kind = pre ) , allocatable , dimension( : ) , public :: x !* x-domain
  real( kind = pre ) , allocatable , dimension( : ) , public :: y !* y-domain
  !*
contains
  !*
  !* setting parameters
  subroutine domain_parametersup()
    implicit none
    at = 0.0_pre
    ax = -1.0_pre
    ay = -1.0_pre
    bx = 1.0_pre
    by = 1.0_pre
    nt = 1110_pin
    nx = 21_pin
    ny = 21_pin
    deltat = 0.0001_pre
    deltax = ( bx - ax ) / ( nx - 1.0_pre )
    deltay = ( by - ay ) / ( ny - 1.0_pre )
  end subroutine domain_parametersup
  !*
  !* setting domain in all directions
  subroutine domain_domain()
    implicit none
    integer( kind = pin ) :: l
    integer( kind = pin ) :: i
    integer( kind = pin ) :: j
    allocate( t( nt ) )
    allocate( x( nx ) )
    allocate( y( ny ) )
    do l = 1_pin , nt
       t( l ) = at + ( l - 1_pin ) * deltat        
    end do
    do i = 1_pin , nx
       x( i ) = ax + ( i - 1_pin ) * deltax
    end do
    do j = 1_pin , ny
       y( j ) = ay + ( j - 1_pin ) * deltay
    end do
  end subroutine domain_domain
  !*
  !* destructing
  subroutine domain_destructor()
    implicit none
    if ( allocated( t ) ) deallocate( t )
    if ( allocated( x ) ) deallocate( x )
    if ( allocated( y ) ) deallocate( y )
  end subroutine domain_destructor
  !*
end module domain
