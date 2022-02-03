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
!*          - T (real array of dimension NT): t-domain (from DOMAIN)
!*          - X (real array of dimension NX): x-domain (from DOMAIN)
!*          - Y (real array of dimension NY): y-domain (from DOMAIN)
!*
!*          - OBSSTEP (integer): index step for observations
!*          - NO (integer): number of steps in which there are observations
!*          - MODESOBS (integer): number of observation modes
!*          - NSTAT (integer): number of stations
!*          - NUMBEROBS (integer): number of observations 
!*          - IFOBS (logical): flag to determine if there are observations
!*          - ERROROBS (real): percentage for observation errors
!*          - INDEXSTATIONSX (integer array): indices for stations in x-direction
!*          - INDEXSTATIONSY (integer array): indices for stations in y-direction 
!*          - TO (real array): time of observations
!*          - STATIONSX (real array): x-coordinates of stations
!*          - STATIONSY (real array): y-coordinates of stations
!*          - SIGMAOBS (real array): variance of observation errors
!*          - OBSVALUE (real array): observation values
!*          - SQRTCOVOBS (real array): square root of the covariance matrix of observation errors
!*          - COVOBS (real array): covariance matrix of observation errors
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
module observations
  !*
  use precision
  use random
  use solution
  use domain
  !*
  implicit none
  !*
  integer( kind = pin ) , public :: obsstep !* index step for observations
  integer( kind = pin ) , public :: no !* number of steps in which there are observations
  integer( kind = pin ) , public :: modesobs !* number of observation modes
  integer( kind = pin ) , public :: nstat !* number of stations
  integer( kind = pin ) , public :: numberobs !* number of observations 
  logical( kind = plo ) , public :: ifobs !* flag to determine if there are observations
  real( kind = pre ) , public :: errorobs !* percentage for observation errors
  integer( kind = pin ) , allocatable , dimension( : ) :: indexstationsx !* indices for stations in x-direction
  integer( kind = pin ) , allocatable , dimension( : ) :: indexstationsy !* indices for stations in y-direction 
  real( kind = pre ) , allocatable , dimension( : ) , public :: to !* time of observations
  real( kind = pre ) , allocatable , dimension( : ) , public :: stationsx !* x-coordinates of stations
  real( kind = pre ) , allocatable , dimension( : ) , public :: stationsy !* y-coordinates of stations
  real( kind = pre ) , allocatable , dimension( : ) , public :: sigmaobs !* variance of observation errors
  real( kind = pre ) , allocatable , dimension( : ) , public :: obsvalue !* observation values
  real( kind = pre ) , allocatable , dimension( : , : ) , public :: sqrtcovobs !* square root of the covariance matrix of observa&
                                                                               !* &tion errors
  real( kind = pre ) , allocatable , dimension( : , : ) , public :: covobs !* covariance matrix of observation errors
  !*
contains
  !*
  !* setting parameters
  subroutine observations_parametersup()
    implicit none
    obsstep = 5_pin
    no = nt / obsstep
    errorobs = 0.30_pre
    modesobs = 12_pin
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
  !* setting observation time
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
  !* setting stations
  subroutine observations_stationsup()
    implicit none
    integer( kind = pin ) :: s 
    nstat = 12_pin
    allocate( indexstationsx( nstat ) )
    indexstationsx(  1_pin )=  6_pin
    indexstationsx(  2_pin )= 14_pin
    indexstationsx(  3_pin )=  3_pin
    indexstationsx(  4_pin )= 10_pin
    indexstationsx(  5_pin )= 18_pin
    indexstationsx(  6_pin )=  6_pin
    indexstationsx(  7_pin )= 14_pin
    indexstationsx(  8_pin )=  3_pin
    indexstationsx(  9_pin )= 10_pin
    indexstationsx( 10_pin )= 18_pin
    indexstationsx( 11_pin )=  6_pin
    indexstationsx( 12_pin )= 14_pin
    allocate( indexstationsy( nstat ) )
    indexstationsy(  1_pin )=  3_pin
    indexstationsy(  2_pin )=  3_pin
    indexstationsy(  3_pin )=  6_pin
    indexstationsy(  4_pin )=  6_pin
    indexstationsy(  5_pin )=  6_pin
    indexstationsy(  6_pin )= 10_pin
    indexstationsy(  7_pin )= 10_pin
    indexstationsy(  8_pin )= 14_pin
    indexstationsy(  9_pin )= 14_pin
    indexstationsy( 10_pin )= 14_pin
    indexstationsy( 11_pin )= 18_pin
    indexstationsy( 12_pin )= 18_pin
    allocate( stationsx( nstat ) )
    do s = 1_pin , nstat
       stationsx( s ) = x( indexstationsx( s ) )
    end do
    allocate( stationsy( nstat ) )
    do s = 1_pin , nstat
       stationsy( s ) = y( indexstationsy( s ) )
    end do
  end subroutine observations_stationsup
  !*
  !* destructing
  subroutine observations_destructor()
    implicit none
    if ( allocated( to ) ) deallocate( to )
    if ( allocated( indexstationsx ) ) deallocate( indexstationsx )
    if ( allocated( indexstationsy ) ) deallocate( indexstationsy )
    if ( allocated( stationsx ) ) deallocate( stationsx )
    if ( allocated( stationsy ) ) deallocate( stationsy )
  end subroutine observations_destructor
  !*
  !* setting number of observations
  subroutine observations_numberobs( l )
    implicit none
    integer( kind = pin ) , intent( in ) :: l
    numberobs = nstat
  end subroutine observations_numberobs
  !*
  !* setting observation value
  subroutine observations_obsvalue( l )
    implicit none
    integer( kind = pin ) , intent( in ) :: l
    integer( kind = pin ) :: s
    real( kind = pre ) :: aux
    do s = 1_pin , numberobs
       aux = solution_solution( t( l ) , stationsx( s ) , stationsy( s ) )
       sigmaobs( s ) = errorobs * abs( aux )
       obsvalue( s ) = random_normal0d( aux , sigmaobs( s ) )
    end do
  end subroutine observations_obsvalue
  !*
  !* tangent of the observation operator
  subroutine observations_tangobsop( l , vectorin , vectorout )
    implicit none
    integer( kind = pin ) , intent( in ) :: l
    real( kind = pre ) , dimension( nx , ny ) , intent( in ) :: vectorin
    real( kind = pre ) , dimension( numberobs ) , intent( out ) :: vectorout
    integer( kind = pin ) :: s
    do s = 1_pin , numberobs
       vectorout( s ) = vectorin( indexstationsx( s ) , indexstationsy( s ) )
    end do
  end subroutine observations_tangobsop
  !*
  !* observation operator
  subroutine observations_obsop( l , vectorin , vectorout )
    implicit none
    integer( kind = pin ) , intent( in ) :: l
    real( kind = pre ) , dimension( nx , ny ) , intent( in ) :: vectorin
    real( kind = pre ) , dimension( numberobs ) , intent( out ) :: vectorout
    integer( kind = pin ) :: s
    do s = 1_pin , numberobs
       vectorout( s ) = vectorin( indexstationsx( s ) , indexstationsy( s ) )
    end do
  end subroutine observations_obsop
  !*
  !* setting square root of the covariance matrix of observation errors
  subroutine observations_sqrtcovobs( l )
    implicit none
    integer( kind = pin ) , intent( in ) :: l
    integer( kind = pin ) :: m
    integer( kind = pin ) , dimension( 1_pin ) :: imax
    real( kind = pre ) :: smax
    real( kind = pre ) , allocatable , dimension( : ) :: sigmaobsaux
    allocate( sigmaobsaux( numberobs ) )
    sigmaobsaux = sigmaobs
    sqrtcovobs = 0.0_pre
    do m = 1_pin , modesobs
       imax = maxloc( sigmaobsaux( 1_pin : numberobs ) )
       smax = maxval( sigmaobsaux( 1_pin : numberobs ) )
       sqrtcovobs( imax( 1_pin ) , m ) = smax
       sigmaobsaux( imax( 1_pin ) ) = 0.0_pre
    end do
    deallocate( sigmaobsaux )
  end subroutine observations_sqrtcovobs
  !*
  !* setting covariance matrix of observation errors
  subroutine observations_covobs( l )
    implicit none
    integer( kind = pin ) , intent( in ) :: l
    integer( kind = pin ) :: s
    covobs = 0.0_pre
    do s = 1_pin , numberobs
       covobs( s , s ) = sigmaobs( s )**2_pin
    end do
  end subroutine observations_covobs
  !*
end module observations
