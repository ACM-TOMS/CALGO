!*********************************************************************************************************************************
!*
!* MODULE: INITIALIZE
!*
!* PURPOSE: initialize filter
!*
!* DEPENDENCIES:
!*               - PRECISION
!*               - TOOLS
!*               - RANDOM
!*               - SOLUTION
!*               - DOMAIN
!* 
!* GLOBALS:
!*          - LOW (integer): number of bytes for simple precision variables (from PRECISION)
!*          - HIGH (integer): number of bytes for double precision variables (from PRECISION)
!*          - PLO (integer): precision for logical variables (from PRECISION)
!*          - PCH (integer): precision for character variables (from PRECISION)
!*          - PIN (integer): precision for integer variables (from PRECISION)
!*          - PRE (integer): precision for real variables (from PRECISION)
!*
!*          - IDUM (integer): seed for random generators, it must be set to a negative value (from RANDOM)
!* 
!*          - BK (real): parameter for solution (from SOLUTION)
!*          - TF (real): parameter for solution (from SOLUTION)
!*          - SF (real): parameter for solution (from SOLUTION)
!*          - DF (real): parameter for solution (from SOLUTION)
!*
!*          - NT (integer): number of nodes in t-direction (from DOMAIN)
!*          - NX (integer): number of nodes in x-direction (from DOMAIN)
!*          - NY (integer): number of nodes in y-direction (from DOMAIN)
!*          - AT (real): left side of the interval in t-direction (from DOMAIN)
!*          - AX (real): left side of the interval in x-direction (from DOMAIN)
!*          - AY (real): left side of the interval in y-direction (from DOMAIN)
!*          - BX (real): right side of the interval in x-direction (from DOMAIN)
!*          - BY (real): right side of the interval in y-direction (from DOMAIN)
!*          - DELTAT (real): step in t-direction (from DOMAIN)
!*          - DELTAX (real): step in x-direction (from DOMAIN)
!*          - DELTAY (real): step in y-direction (from DOMAIN)
!*          - T (real array): t-domain (from DOMAIN)
!*          - X (real array): x-domain (from DOMAIN)
!*          - Y (real array): y-domain (from DOMAIN)
!*
!*          - DIMSPACESTATE (integer): dimension of the space state
!*          - NUMBERSAMPLES (integer): number of samples
!*          - MODESANALYSIS (integer): number of analysis modes
!*          - FIRST (logical): flag to determine the first step
!*          - ERRORIN (real): a measure of the initial error
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
module initialize
  !*
  use precision
  use tools
  use random
  use solution
  use domain
  !*
  implicit none
  !*
  integer( kind = pin ) , public :: dimspacestate !* dimension of the space state
  integer( kind = pin ) , public :: numbersamples !* number of samples
  integer( kind = pin ) , public :: modesanalysis !* number of analysis modes
  logical( kind = plo ) , public :: first !* flag to determine the first step
  real( kind = pre ) , public :: errorin !* a measure of the initial error
  !*
contains
  !*
  !* set initialization parameters
  subroutine initialize_parametersup()
    implicit none
    errorin = 10.0_pre
    dimspacestate = nx * ny
    numbersamples = 50_pin
    modesanalysis = 200_pin
    first = .true.
  end subroutine initialize_parametersup
  !*
  !* set initial state and initial covariance matrix
  subroutine initialize_initialize1( state , covariance )
    implicit none
    real( kind = pre ) , dimension( dimspacestate ) , intent( out ) :: state
    real( kind = pre ) , dimension( dimspacestate , dimspacestate ) , intent( out ) :: covariance
    integer( kind = pin ) :: i
    integer( kind = pin ) :: j
    integer( kind = pin ) :: ij
    real( kind = pre ) :: aux
    real( kind = pre ) :: sigmain
    covariance = 0.0_pre
    ij = 1_pin
    do j = 1_pin , ny
       do i = 1_pin , nx
          aux = solution_solution( t( 1_pin ) , x( i ) , y( j ) )
          sigmain = errorin * abs( aux )
          state( ij ) = random_normal0d( aux , sigmain )
          covariance( ij , ij ) = sigmain**2_pin
          ij = ij + 1_pin
       end do
    end do
  end subroutine initialize_initialize1
  !*
  !* set initial state and initial covariance matrix
  subroutine initialize_initialize2( state , sqrtcov )
    implicit none
    real( kind = pre ) , dimension( dimspacestate ) , intent( out ) :: state
    real( kind = pre ) , dimension( dimspacestate , modesanalysis ) , intent( out ) :: sqrtcov
    integer( kind = pin ) :: i
    integer( kind = pin ) :: j
    integer( kind = pin ) :: ij
    integer( kind = pin ) :: m
    integer( kind = pin ) , dimension( 2_pin ) :: ijmax
    real( kind = pre ) :: aux
    real( kind = pre ) :: smax
    real( kind = pre ) , allocatable , dimension( : , : ) :: sigmain
    sqrtcov = 0.0_pre
    allocate( sigmain( nx , ny ) )
    ij = 1_pin
    do j = 1_pin , ny
       do i = 1_pin , nx
          aux = solution_solution( t( 1_pin ) , x( i ) , y( j ) )
          sigmain( i , j ) = errorin * abs( aux )
          state( ij ) = random_normal0d( aux , sigmain( i , j ) )
          ij = ij + 1_pin
       end do
    end do
    do m = 1_pin , modesanalysis
       ijmax = maxloc( sigmain( 1_pin : nx , 1_pin : ny ) )
       smax = maxval( sigmain( 1_pin : nx , 1_pin : ny ) )
       call tools_ij2s( nx , ijmax( 1_pin ) , ijmax( 2_pin ) , ij )
       sqrtcov( ij , m ) = smax
       sigmain( ijmax( 1_pin ) , ijmax( 2_pin ) ) = 0.0_pre
    end do
    deallocate( sigmain )
  end subroutine initialize_initialize2
  !*
end module initialize
