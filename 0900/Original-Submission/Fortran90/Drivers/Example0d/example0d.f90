!*********************************************************************************************************************************
!*
!* PROGRAM: EXAMPLE2D_EKF 
!*
!* PURPOSE: test the Extended Kalman Filter
!*
!* DEPENDENCIES:
!*                - PRECISION
!*                - TOOLS
!*                - RANDOM
!*                - PARAMETERS
!*                - SOLUTION
!*                - DOMAIN
!*                - MODEL
!*                - OBSERVATIONS
!*                - PARALLEL
!*                - INITIALIZE
!*                - EKF
!*                - ENKF
!*                - RRSQRTKF
!*                - RRSQRTENKF
!*
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
!*          - UNITTRUTH (integer): file unit for truth (from PARAMETERS)
!*          - UNITMODEL (integer): file unit for model (from PARAMETERS)
!*          - UNITKF (integer): file unit for ekf (from PARAMETERS)
!*          - UNITRRSQRTKF (integer): file unit for rrsqrtkf (from PARAMETERS)
!*          - UNITENKF (integer): file unit for enkf (from PARAMETERS)
!*          - UNITRRSQRTENKF (integer): file unit for rrsqrtenkf (from PARAMETERS)
!*          - UNITOBSER (integer): file unit for observations (from PARAMETERS)
!*          - FILETRUTH (character*100): file name for truth (from PARAMETERS)
!*          - FILEMODEL (character*100): file name for model (from PARAMETERS)
!*          - FILEEKF (character*100): file name for ekf (from PARAMETERS)
!*          - FILERRSQRTKF (character*100): file name for rrsqrtkf (from PARAMETERS)
!*          - FILEENKF (character*100): file name for enkf (from PARAMETERS)
!*          - FILERRSQRTENKF (character*100): file name for rrsqrtenkf (from PARAMETERS)
!*          - FILEOBSER (character*100): file name for observations (from PARAMETERS)
!*
!*          - NT (integer): number of time steps (from DOMAIN)
!*          - A (real): initial time (from DOMAIN)
!*          - DELTAT (real): time step (from DOMAIN)
!*          - T (real array): vector of time (from DOMAIN)
!*
!*          - MODESMODEL (integer): number of model modes (from MODEL)
!*          - MODELERROR (real array): covariance matrix of model errors (from MODEL)
!*          - SQRTMODELERROR (real array): square root of the covariance matrix of model errors (from MODEL)
!*
!*          - OBSSTEP (integer): index step for observations (from OBSERVATIONS)
!*          - NO (integer): number of steps in which there are observations (from OBSERVATIONS)
!*          - NUMBEROBS (integer): number of observations (from OBSERVATIONS)
!*          - MODESOBS (integer): number of observation modes (from OBSERVATIONS)
!*          - IFOBS (logical): flag to determine if there are observations (from OBSERVATIONS)
!*          - ERROROBS (real): percentage for observation errors (from OBSERVATIONS)
!*          - TO (real array): time of observations (from OBSERVATIONS)
!*          - OBSVALUE (real array): observation values (from OBSERVATIONS)
!*          - SIGMAOBS (real array): variance of observation errors (from OBSERVATIONS)
!*          - COVOBS (real array): covariance matrix of observation errors (from OBSERVATIONS)
!*          - SQRTCOVOBS (real array): square root of the covariance matrix of observation errors (from OBSERVATIONS)
!*
!*          - IERROR (integer): variable to detect errors (from PARALLEL)
!*          - RANK (integer): which machine (from PARALLEL)
!*          - NPROC (integer): number of processors (from PARALLEL)
!*          - BLOCK_SIZE (integer): block size (from PARALLEL)
!*          - NPROW (integer): number of rows in the processor grid (from PARALLEL)
!*          - NPCOL (integer): number of columns in the processor grid (from PARALLEL)
!*          - CONTEXT (integer): it is a universe where messages exist and do not interact with other context's messages (from PA&
!*                               &RALLEL)
!*          - MYROW (integer): calling processor's row number in the processor grid (from PARALLEL)
!*          - MYCOL (integer): calling processor's column number in the processor grid (from PARALLEL)
!*
!*          - DIMSPACESTATE (integer): dimension of the space state (from INITIALIZE)
!*          - NUMBERSAMPLES (integer): number of samples (from INITIALIZE)
!*          - MODESANALYSIS (integer): number of analysis modes (from INITIALIZE)
!*          - FIRST (logical): flag to determine the first step (from INITIALIZE)
!*          - ERRORIN (real): a measure of the initial error (from INITIALIZE)
!*
!* SCOPE: this program belongs to the example
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
program example0d
  !*
  use precision
  use tools
  use random
  use parameters
  use solution
  use domain
  use model
  use observations
  use initialize
  use ekf
  use rrsqrtkf
  use enkf
  use rrsqrtenkf
  use precision
  use parallel
  !*
  implicit none
  !*
  integer( kind = pin ) :: l
  integer( kind = pin ) :: lo
  real( kind = pre ) , allocatable , dimension( : ) :: state
  real( kind = pre ) , allocatable , dimension( : ) :: stateekf
  real( kind = pre ) , allocatable , dimension( : ) :: staterrsqrtkf
  real( kind = pre ) , allocatable , dimension( : ) :: stateenkf
  real( kind = pre ) , allocatable , dimension( : ) :: staterrsqrtenkf
  real( kind = pre ) , allocatable , dimension( : ) :: truth
  real( kind = pre ) , allocatable , dimension( : ) :: modl
  real( kind = pre ) , allocatable , dimension( : ) :: assimekf
  real( kind = pre ) , allocatable , dimension( : ) :: assimrrsqrtkf
  real( kind = pre ) , allocatable , dimension( : ) :: assimenkf
  real( kind = pre ) , allocatable , dimension( : ) :: assimrrsqrtenkf
  real( kind = pre ) , allocatable , dimension( : ) :: obser
  real( kind = pre ) , allocatable , dimension( : , : ) :: covariance
  real( kind = pre ) , allocatable , dimension( : , : ) :: sqrtcov
  real( kind = pre ) , allocatable , dimension( : , : ) :: sqrtcovaux
  real( kind = pre ) , allocatable , dimension( : , : ) :: sqrtcovfor
  real( kind = pre ) , allocatable , dimension( : , : ) :: covarianceekf
  real( kind = pre ) , allocatable , dimension( : , : ) :: covarianceenkf
  real( kind = pre ) , allocatable , dimension( : , : ) :: sqrtcovrrsqrtkf
  real( kind = pre ) , allocatable , dimension( : , : ) :: sqrtcovrrsqrtenkf
  real( kind = pre ) , allocatable , dimension( : , : ) :: samplesenkf
  real( kind = pre ) , allocatable , dimension( : , : ) :: samplesrrsqrtenkf
  !*
  !* initialization of parallelization
  call parallel_init()
  !*
  !* getting parallelization parameters 
  call parallel_ranksize()
  !*
  !* precision parameters
  print*,''
  print*,'precision parameters'
  print*,'plo = ' , plo
  print*,'pch = ' , pch
  print*,'pin = ' , pin
  print*,'pre = ' , pre
  !*
  !* random parameters
  call random_parametersup()
  print*,''
  print*,'random parameters'
  print*,'idum = ' , idum
  !*
  !* parameters
  call parameters_parametersup()
  print*,''
  print*,'local parameters'
  print*,'unittruth       = ' , unittruth
  print*,'unitmodel       = ' , unitmodel
  print*,'unitekf         = ' , unitekf
  print*,'unitrrsqrtkf    = ' , unitrrsqrtkf
  print*,'unitenkf        = ' , unitenkf
  print*,'unitrrsqrtenkf  = ' , unitrrsqrtenkf
  print*,'unitobser       = ' , unitobser
  print*,'filetruth       = ' , trim( filetruth )
  print*,'filemodel       = ' , trim( filemodel )
  print*,'fileekf         = ' , trim( fileekf )
  print*,'filerrsqrtkf    = ' , trim( filerrsqrtkf )
  print*,'fileenkf        = ' , trim( fileenkf )
  print*,'filerrsqrtenkf  = ' , trim( filerrsqrtenkf )
  print*,'fileobser       = ' , trim( fileobser )
  !*
  !* domain parameters
  call domain_parametersup()
  print*,''
  print*,'domain parameters'
  print*,'nt     = ' , nt
  print*,'a      = ' , a
  print*,'deltat = ' , deltat
  !*
  !* domain
  call domain_domain()
  print*,''
  print*,'domain'
  do l = 1_pin , nt
     print*,'t(',l,') = ' , t( l )
  end do
  !*
  !* model
  call model_parametersup()
  print*,''
  print*,'model parameters'
  print*,'modesmodel = ' , modesmodel
  !*
  !* observations parameters
  call observations_parametersup()
  print*,''
  print*,'observations'
  print*,'obsstep  = ' , obsstep
  print*,'no       = ' , no
  print*,'modesobs = ' , modesobs
  print*,'errorobs = ' , errorobs
  !*
  !* observations time
  call observations_obstimeup()
  print*,''
  print*,'observations time'
  do lo = 1_pin , no
     print*,'to(',lo,') = ' , to( lo )
  end do
  !*
  !* assimilation parameters
  call initialize_parametersup()
  print*,''
  print*,'initialization parameters'
  print*,'dimspacestate = ' , dimspacestate
  print*,'numbersamples = ' , numbersamples
  print*,'modesanalysis = ' , modesanalysis
  print*,'first         = ' , first
  print*,'errorin       = ' , errorin
  !*
  !* opening files
  open( unit = unittruth , file = filetruth )
  open( unit = unitmodel , file = filemodel )
  open( unit = unitekf , file = fileekf )
  open( unit = unitrrsqrtkf , file = filerrsqrtkf )
  open( unit = unitenkf , file = fileenkf )
  open( unit = unitrrsqrtenkf , file = filerrsqrtenkf )
  open( unit = unitobser , file = fileobser )
  !*
  !* allocations
  allocate( state( dimspacestate ) )
  allocate( stateekf( dimspacestate ) )
  allocate( staterrsqrtkf( dimspacestate ) )
  allocate( stateenkf( dimspacestate ) )
  allocate( staterrsqrtenkf( dimspacestate ) )
  allocate( covariance( dimspacestate , dimspacestate ) )
  allocate( sqrtcov( dimspacestate , modesanalysis ) )
  allocate( sqrtcovaux( dimspacestate , numbersamples ) )
  allocate( sqrtcovfor( dimspacestate , modesanalysis + modesmodel ) )
  allocate( sqrtcovrrsqrtkf( dimspacestate , modesanalysis ) )
  allocate( sqrtcovrrsqrtenkf( dimspacestate , modesanalysis ) )
  allocate( covarianceekf( dimspacestate , dimspacestate ) )
  allocate( covarianceenkf( dimspacestate , dimspacestate ) )
  allocate( truth( nt ) )
  allocate( modl( nt ) )
  allocate( assimekf( nt ) )
  allocate( assimrrsqrtkf( nt ) )
  allocate( assimenkf( nt ) )
  allocate( assimrrsqrtenkf( nt ) )
  allocate( obser( no ) )
  allocate( samplesenkf( dimspacestate , numbersamples ) )
  allocate( samplesrrsqrtenkf( dimspacestate , numbersamples ) )
  !*
  !************************************************ step 1 ***********************************************************************
  !*
  print*,''
  print*,'resolution'
  !*
  !* truth
  truth( 1_pin ) = solution_solution( t( 1_pin ) )
  !*
  !* initialization
  call initialize_initialize1( state , covariance )
  call initialize_initialize2( state , sqrtcov )
  !*
  stateekf = state
  covarianceekf = covariance
  !*
  staterrsqrtkf = state
  sqrtcovrrsqrtkf = sqrtcov
  !*
  stateenkf = state
  covarianceenkf = covariance
  !*
  staterrsqrtenkf = state
  sqrtcovrrsqrtenkf = sqrtcov
  !*
  !* model
  modl( 1_pin ) = state( 1_pin )
  !*
  !* assimilation
  assimekf( 1_pin ) = stateekf( 1_pin )
  assimrrsqrtkf( 1_pin ) = staterrsqrtkf( 1_pin )
  assimenkf( 1_pin ) = stateenkf( 1_pin )
  assimrrsqrtenkf( 1_pin ) = staterrsqrtenkf( 1_pin )
  !*
  !* print results
  write( unittruth , * ) t( 1_pin ) , truth( 1_pin )
  write( unitmodel , * ) t( 1_pin ) , modl( 1_pin )
  write( unitekf , * ) t( 1_pin ) , assimekf( 1_pin )
  write( unitrrsqrtkf , * ) t( 1_pin ) , assimrrsqrtkf( 1_pin )
  write( unitenkf , * ) t( 1_pin ) , assimenkf( 1_pin )
  write( unitrrsqrtenkf , * ) t( 1_pin ) , assimrrsqrtenkf( 1_pin )
  !*
  !* printing on screen
  write( * , '("step: ",I5)' ) 1_pin 
  !*
  !****************************************************** step l + 1 *************************************************************
  !*
  lo = 1_pin
  !*
  do l = 1_pin , nt - 1_pin
     !*
     !* truth
     truth( l + 1_pin ) = solution_solution( t( l + 1_pin ) )
     !*
     !* model
     call model_model( l , modl( l ) , modl( l + 1_pin ) )
     !*
     !* assim
     call ekf_predictor( l , stateekf , covarianceekf )
     call rrsqrtkf_predictor( l , staterrsqrtkf , sqrtcovrrsqrtkf , sqrtcovfor )
     if ( l == 1_pin ) then
        call enkf_predictor( l , stateenkf , covarianceenkf , samplesenkf )
        first = .true.
        call rrsqrtenkf_predictor( l , staterrsqrtenkf , sqrtcovrrsqrtenkf , samplesrrsqrtenkf , sqrtcovaux )
     else
        call enkf_predictor( l , stateenkf , covarianceenkf , samplesenkf )
        call rrsqrtenkf_predictor( l , staterrsqrtenkf , sqrtcovrrsqrtenkf , samplesrrsqrtenkf , sqrtcovaux )
     end if
     !*
     !* decide if there are observations
     call observations_ifobservations( l + 1_pin )
     !*
     if ( ifobs ) then
        !*
        !* set number of observations
        call observations_numberobs( l + 1_pin )
        !*
        allocate( obsvalue( numberobs ) )
        allocate( sigmaobs( numberobs ) )
        !*
        !* set observation values
        call observations_obsvalue( l + 1_pin )
        !*
        obser( lo ) = obsvalue( 1_pin ) 
        !*
        write( unitobser , * ) to( lo ) , obser( lo )
        !*
        !* correction step
        call ekf_corrector( l , stateekf , covarianceekf )
        call rrsqrtkf_corrector( l , staterrsqrtkf , sqrtcovrrsqrtkf , sqrtcovfor )
        call enkf_corrector( l , stateenkf , covarianceenkf , samplesenkf )
        call rrsqrtenkf_corrector( l , staterrsqrtenkf , sqrtcovrrsqrtenkf , samplesrrsqrtenkf , sqrtcovaux )
        !*
        deallocate( obsvalue )
        deallocate( sigmaobs )
        !*
        lo = lo + 1_pin
        !*
     end if
     !*
     assimekf( l + 1_pin ) = stateekf( 1_pin )
     assimrrsqrtkf( l + 1_pin ) = staterrsqrtkf( 1_pin )
     assimenkf( l + 1_pin ) = stateenkf( 1_pin )
     assimrrsqrtenkf( l + 1_pin ) = staterrsqrtenkf( 1_pin )
     !*
     !* print results
     write( unittruth , * ) t( l + 1_pin ) , truth( l + 1_pin )
     write( unitmodel , * ) t( l + 1_pin ) , modl( l + 1_pin )
     write( unitekf , * ) t( l + 1_pin ) , assimekf( l + 1_pin )
     write( unitrrsqrtkf , * ) t( l + 1_pin ) , assimrrsqrtkf( l + 1_pin )
     write( unitenkf , * ) t( l + 1_pin ) , assimenkf( l + 1_pin )
     write( unitrrsqrtenkf , * ) t( l + 1_pin ) , assimrrsqrtenkf( l + 1_pin )
     !*
     !* printing on screen
     if ( ifobs ) then
        write( * , '("step: " ,I5," --> assimilation")' ) l + 1_pin
     else
        write( * , '("step: ",I5)' ) l + 1_pin 
     end if
     !*
  end do
  !*
  !*******************************************************************************************************************************
  !*
  !* deallocation
  call domain_destructor()
  call observations_destructor()
  deallocate( state )
  deallocate( stateekf )
  deallocate( staterrsqrtkf )
  deallocate( stateenkf )
  deallocate( staterrsqrtenkf )
  deallocate( covariance )
  deallocate( sqrtcov )
  deallocate( sqrtcovaux )
  deallocate( sqrtcovfor )
  deallocate( covarianceekf )
  deallocate( covarianceenkf )
  deallocate( sqrtcovrrsqrtkf )
  deallocate( sqrtcovrrsqrtenkf )
  deallocate( truth )
  deallocate( modl )
  deallocate( assimekf )
  deallocate( assimrrsqrtkf )
  deallocate( assimenkf )
  deallocate( assimrrsqrtenkf )
  deallocate( obser )
  deallocate( samplesenkf )
  deallocate( samplesrrsqrtenkf )
  !*
  !* closing files
  close( unittruth )
  close( unitmodel )
  close( unitekf )
  close( unitrrsqrtkf )
  close( unitenkf )
  close( unitrrsqrtenkf )
  close( unitobser )
  !*
  !* end parallelization
  call parallel_finalize()
  !*
end program example0d

