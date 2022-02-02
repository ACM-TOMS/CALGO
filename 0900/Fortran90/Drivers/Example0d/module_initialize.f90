!*********************************************************************************************************************************
!*
!* MODULE: INITIALIZE
!*
!* PURPOSE: initialize filter
!*
!* DEPENDENCIES:
!*               - PRECISION
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
!*          - NT (integer): number of time steps (from DOMAIN)
!*          - A (real): initial time (from DOMAIN)
!*          - DELTAT (real): time step (from DOMAIN)
!*          - T (real array): vector of time (from DOMAIN)
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
!*
module initialize
  !*
  use precision
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
    dimspacestate = 1_pin
    numbersamples = 100_pin
    modesanalysis = 1_pin
    first = .true.
    errorin = 2.5_pre
  end subroutine initialize_parametersup
  !*
  !* set initial state and initial covariance matrix
  subroutine initialize_initialize1( state , covariance )
    implicit none
    real( kind = pre ) , dimension( dimspacestate ) , intent( out ) :: state
    real( kind = pre ) , dimension( dimspacestate , dimspacestate ) , intent( out ) :: covariance
!!$    real( kind = pre ) :: aux
!!$    real( kind = pre ) :: sigmain
    !*----------------------------------------------------------------------------------------------------------------------------
    !* The state and covariance arrays are set in order you get the same picture plotted in the paper. In order to play a little    
    !* with the initial conditions, you should remove the following two lines and uncomment the commented lines.
    state( 1_pin ) = 4.0_pre
    covariance( 1_pin , 1_pin ) = 2.0**2_pin
!!$    aux = solution_solution( t( 1_pin ) )
!!$    sigmain = errorin * abs( aux )
!!$    state( 1_pin ) = random_normal0d( aux , sigmain )
!!$    covariance( 1_pin , 1_pin ) = sigmain**2_pin
  end subroutine initialize_initialize1
  !*
  !* set initial state and initial square root covariance matrix
  subroutine initialize_initialize2( state , sqrtcov )
    implicit none
    real( kind = pre ) , dimension( dimspacestate ) , intent( out ) :: state
    real( kind = pre ) , dimension( dimspacestate , modesanalysis ) , intent( out ) :: sqrtcov
!!$    real( kind = pre ) :: aux
!!$    real( kind = pre ) :: sigmain
    !*----------------------------------------------------------------------------------------------------------------------------
    !* The state and covariance arrays are set in order you get the same picture plotted in the paper. In order to play a little    
    !* with the initial conditions, you should remove the following two lines and uncomment the commented lines.
    state( 1_pin ) = 4.0_pre
    sqrtcov( 1_pin , 1_pin ) = abs( 2.0_pre )
!!$    aux = solution_solution( t( 1_pin ) )
!!$    sigmain = errorin * abs( aux )
!!$    state( 1_pin ) = random_normal0D( aux , sigmain )
!!$    sqrtcov( 1_pin , 1_pin ) = abs( sigmain )
  end subroutine initialize_initialize2  
  !*
end module initialize
