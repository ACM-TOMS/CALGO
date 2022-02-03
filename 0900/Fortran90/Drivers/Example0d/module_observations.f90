!*********************************************************************************************************************************
!*
!* MODULE: OBSERVATIONS
!*
!* PURPOSE: sets the observations and covariance matrices
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
!*          - OBSSTEP (integer): index step for observations
!*          - NO (integer): number of steps in which there are observations
!*          - NUMBEROBS (integer): number of observations
!*          - MODESOBS (integer): number of observation modes
!*          - IFOBS (logical): flag to determine if there are observations
!*          - ERROROBS (real): percentage for observation errors
!*          - TO (real array): time of observations
!*          - OBSVALUE (real array): observation values
!*          - SIGMAOBS (real array): variance of observation errors
!*          - COVOBS (real array): covariance matrix of observation errors
!*          - SQRTCOVOBS (real array): square root of the covariance matrix of observation errors
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
module observations
  !*
  use precision
  use random
  use solution
  use domain
  !*
  implicit none
  !*
  integer( kind = pin ) , public :: obsstep !* number of time steps to have an observation
  integer( kind = pin ) , public :: no !* number of steps where we have observations
  integer( kind = pin ) , public :: numberobs !* number of observations
  integer( kind = pin ) , public :: modesobs !* number of observation modes
  logical( kind = plo ) , public :: ifobs !* flag to determine if there are observations
  real( kind = pre ) , public :: errorobs !* a measure of the observation error
  real( kind = pre ) , allocatable , dimension( : ) , public :: to !* observation time
  real( kind = pre ) , allocatable , dimension( : ) , public :: obsvalue !* observation value
  real( kind = pre ) , allocatable , dimension( : ) , public :: sigmaobs !* variance to generate noise to create observations
  real( kind = pre ) , allocatable , dimension( : , : ) , public :: covobs !* covariance matrix of observation errors
  real( kind = pre ) , allocatable , dimension( : , : ) , public :: sqrtcovobs !* square root observation errors covariance matrix
  !*
contains
  !*
  !* set parameters
  subroutine observations_parametersup()
    implicit none
    obsstep = 5_pin
    no = nt / obsstep
    modesobs = 1_pin 
    errorobs = 0.30_pre
  end subroutine observations_parametersup
  !*
  !* decide if there are observations
  subroutine observations_ifobservations( l )
    implicit none
    integer( kind = pin ) , intent( in ) :: l
    if ( mod( l , obsstep ) == 0_pin ) then
       ifobs = .true.
    else
       ifobs = .false.
    end if
  end subroutine observations_ifobservations
  !*
  !* set observation time
  subroutine observations_obstimeup()
    implicit none
    integer( kind = pin ) :: l
    integer( kind = pin ) :: lo
    allocate( to( no ) )
    lo = 1_pin
    do l = 2_pin , nt
       call observations_ifobservations( l )
       if ( ifobs ) then
          to( lo ) = t( l )
          lo = lo + 1_pin
       end if
    end do
  end subroutine observations_obstimeup
  !*
  !* free memory
  subroutine observations_destructor()
    implicit none
    if ( allocated( to ) ) deallocate( to )
  end subroutine observations_destructor
  !*
  !* set number of observations
  subroutine observations_numberobs( l )
    implicit none
    integer( kind = pin ) , intent( in ) :: l
    numberobs = 1_pin
  end subroutine observations_numberobs
  !*
  !* set observation values
  subroutine observations_obsvalue( l )
    implicit none
    integer( kind = pin ) , intent( in ) :: l
    real( kind = pre ) :: aux
    aux = solution_solution( t( l ) )
    sigmaobs( 1_pin ) = errorobs * abs( aux )
    obsvalue( 1_pin ) = random_normal0d( aux , sigmaobs( 1_pin ) )
  end subroutine observations_obsvalue
  !*
  !* observation operator
  subroutine observations_obsop( l , vectorin , vectorout )
    implicit none
    integer( kind = pin ) , intent( in ) :: l
    real( kind = pre ) , dimension( 1_pin ) , intent( in ) :: vectorin
    real( kind = pre ) , dimension( numberobs ) , intent( out ) :: vectorout
    vectorout = vectorin
  end subroutine observations_obsop
  !*
  !* tangent observation operator
  subroutine observations_tangobsop( l , vectorin , vectorout )
    implicit none
    integer( kind = pin ) , intent( in ) :: l
    real( kind = pre ) , dimension( 1_pin ) , intent( in ) :: vectorin
    real( kind = pre ) , dimension( numberobs ) , intent( out ) :: vectorout
    vectorout = vectorin
  end subroutine observations_tangobsop
  !*
  !* set covariance matrix of observation errors
  subroutine observations_covobs( l )
    implicit none
    integer( kind = pin ) , intent( in ) :: l
    covobs( 1_pin , 1_pin ) = sigmaobs( 1_pin )**2_pin
  end subroutine observations_covobs
  !*
  !* set the square root observation errors covariance matrix
  subroutine observations_sqrtcovobs( l )
    implicit none
    integer( kind = pin ) , intent( in ) :: l
    sqrtcovobs( 1_pin , 1_pin ) = abs( sigmaobs( 1_pin ) )
  end subroutine observations_sqrtcovobs
  !*
end module observations
