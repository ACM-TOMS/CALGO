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
!*          - UNITASSIMEKF (integer): file unit for ekf (from PARAMETERS)     
!*          - UNITASSIMRRSQRTKF (integer): file unit for rrsqrtkf (from PARAMETERS)     
!*          - UNITASSIMENKF (integer): file unit for enkf (from PARAMETERS)     
!*          - UNITASSIMRRSQRTENKF (integer): file unit for rrsqrtenkf (from PARAMETERS)     
!*          - UNITOBSER (integer): file unit for observations (from PARAMETERS)    
!*          - FILETRUTH (character*100): file name for truth (from PARAMETERS)          
!*          - FILEMODEL (character*100): file name for model (from PARAMETERS)     
!*          - FILEASSIMEKF (character*100): file name for ekf (from PARAMETERS)   
!*          - FILEASSIMRRSQRTKF (character*100): file name for rrsqrtkf (from PARAMETERS)   
!*          - FILEASSIMENKF (character*100): file name for enkf (from PARAMETERS)     
!*          - FILEASSIMRRSQRTENKF (character*100): file name for rrsqrtenkf (from PARAMETERS)     
!*          - FILEOBSER (character*100): file name for observations (from PARAMETERS)     
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
!*          - MODESMODEL (integer): number of model modes (from MODEL)
!*          - ERRORBC (real): percentage of boundary condition errors (from MODEL)
!*          - BC (real array): boundary condition (from MODEL)
!*          - SIGMABC (real array): standard deviation of boundary condition (from MODEL)
!*          - SQRTMODELERROR (real array): square root of covariance matrix of model errors (from MODEL)
!*          - MODELERROR (real array): covariance matrix of model errors (from MODEL)
!*
!*          - OBSSTEP (integer): index step for observations (from OBSERVATIONS)
!*          - NO (integer): number of steps in which there are observations (from OBSERVATIONS)
!*          - MODESOBS (integer): number of observation modes (from OBSERVATIONS)
!*          - NSTAT (integer): number of stations (from OBSERVATIONS)
!*          - NUMBEROBS (integer): number of observations (from OBSERVATIONS)
!*          - IFOBS (logical): flag to determine if there are observations (from OBSERVATIONS)
!*          - ERROROBS (real): percentage for observation errors (from OBSERVATIONS)
!*          - INDEXSTATIONSX (integer array): indices for stations in x-direction (from OBSERVATIONS)
!*          - INDEXSTATIONSY (integer array): indices for stations in y-direction (from OBSERVATIONS) 
!*          - TO (real array): time of observations (from OBSERVATIONS)
!*          - STATIONSX (real array): x-coordinates of stations (from OBSERVATIONS)
!*          - STATIONSY (real array): y-coordinates of stations (from OBSERVATIONS)
!*          - SIGMAOBS (real array): variance of observation errors (from OBSERVATIONS)
!*          - OBSVALUE (real array): observation values (from OBSERVATIONS)
!*          - SQRTCOVOBS (real array): square root of the covariance matrix of observation errors (from OBSERVATIONS)
!*          - COVOBS (real array): covariance matrix of observation errors (from OBSERVATIONS)
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
program example2d_ekf
  !*
  use precision
  use tools
  use random
  use parameters
  use solution
  use domain
  use model
  use observations
  use parallel
  use initialize
  use ekf
  !*
  implicit none
  !*
  integer( kind = pin ) :: l
  integer( kind = pin ) :: i
  integer( kind = pin ) :: j
  integer( kind = pin ) :: lo
  integer( kind = pin ) :: s
  real( kind = pre ) , allocatable , dimension( : , : ) :: state
  real( kind = pre ) , allocatable , dimension( : , : ) :: stateekf
  real( kind = pre ) , allocatable , dimension( : , : ) :: covariance
  real( kind = pre ) , allocatable , dimension( : , : ) :: covarianceekf
  real( kind = pre ) , allocatable , dimension( : , : , : ) :: truth
  real( kind = pre ) , allocatable , dimension( : , : , : ) :: modl
  real( kind = pre ) , allocatable , dimension( : , : ) :: obser
  real( kind = pre ) , allocatable , dimension( : , : , : ) :: assimekf
  !*
  !* initialization of parallelization
  call parallel_init()
  !*
  !* getting parallelization parameters 
  call parallel_ranksize()
  !*
  !* precision parameters
  if ( rank == 0_pin ) then
     print*,''
     print*,'precision parameters'
     print*,'plo = ' , plo
     print*,'pch = ' , pch
     print*,'pin = ' , pin
     print*,'pre = ' , pre
  end if
  !*
  !* random parameters
  if ( rank == 0_pin ) then
     call random_parametersup()
     print*,''
     print*,'random parameters'
     print*,'idum = ' , idum
  end if
  call mpi_bcast( idum , 1_pin , mpi_integer , 0_pin , mpi_comm_world , ierror )
  !*
  !* parameters
  if ( rank == 0_pin ) then
     call parameters_parametersup()
     print*,''
     print*,'parameters'
     print*,'unittruth    = ' , unittruth     
     print*,'unitmodel    = ' , unitmodel   
     print*,'unitassimekf = ' , unitassimekf   
     print*,'unitobser    = ' , unitobser   
     print*,'filetruth    = ' , trim( filetruth )   
     print*,'filemodel    = ' , trim( filemodel )  
     print*,'fileassimekf = ' , trim( fileassimekf )  
     print*,'fileobser    = ' , trim( fileobser )  
  end if
  call mpi_bcast( unittruth    , 1_pin , mpi_integer , 0_pin , mpi_comm_world , ierror )
  call mpi_bcast( unitmodel    , 1_pin , mpi_integer , 0_pin , mpi_comm_world , ierror )
  call mpi_bcast( unitassimekf , 1_pin , mpi_integer , 0_pin , mpi_comm_world , ierror )
  call mpi_bcast( unitobser    , 1_pin , mpi_integer , 0_pin , mpi_comm_world , ierror )
  call mpi_bcast( filetruth    , 100_pin , mpi_character , 0_pin , mpi_comm_world , ierror )
  call mpi_bcast( filemodel    , 100_pin , mpi_character , 0_pin , mpi_comm_world , ierror )
  call mpi_bcast( fileassimekf , 100_pin , mpi_character , 0_pin , mpi_comm_world , ierror )
  call mpi_bcast( fileobser    , 100_pin , mpi_character , 0_pin , mpi_comm_world , ierror )
  !*
  !* solution parameters
  if ( rank == 0_pin ) then
     call solution_parametersup()
     print*,''
     print*,'solution parameters'
     print*,'bk = ' , bk
     print*,'tf = ' , tf
     print*,'sf = ' , sf
     print*,'df = ' , df
  end if
  if ( pre == low ) then
     call mpi_bcast( bk , 1_pin , mpi_real , 0_pin , mpi_comm_world , ierror )
     call mpi_bcast( tf , 1_pin , mpi_real , 0_pin , mpi_comm_world , ierror )
     call mpi_bcast( sf , 1_pin , mpi_real , 0_pin , mpi_comm_world , ierror )
     call mpi_bcast( df , 1_pin , mpi_real , 0_pin , mpi_comm_world , ierror )
  else if ( pre == high ) then
     call mpi_bcast( bk , 1_pin , mpi_double_precision , 0_pin , mpi_comm_world , ierror )
     call mpi_bcast( tf , 1_pin , mpi_double_precision , 0_pin , mpi_comm_world , ierror )
     call mpi_bcast( sf , 1_pin , mpi_double_precision , 0_pin , mpi_comm_world , ierror )
     call mpi_bcast( df , 1_pin , mpi_double_precision , 0_pin , mpi_comm_world , ierror )
  end if
  !*
  !* domain parameters
  if ( rank == 0_pin ) then
     call domain_parametersup()
     print*,''
     print*,'domain parameters'
     print*,'at     = ' , at     
     print*,'ax     = ' , ax    
     print*,'ay     = ' , ay    
     print*,'bx     = ' , bx    
     print*,'by     = ' , by    
     print*,'nt     = ' , nt    
     print*,'nx     = ' , nx    
     print*,'ny     = ' , ny    
     print*,'deltat = ' , deltat
     print*,'deltax = ' , deltax
     print*,'deltay = ' , deltay
  end if
  call mpi_bcast( nt , 1_pin , mpi_integer , 0_pin , mpi_comm_world , ierror )
  call mpi_bcast( nx , 1_pin , mpi_integer , 0_pin , mpi_comm_world , ierror )
  call mpi_bcast( ny , 1_pin , mpi_integer , 0_pin , mpi_comm_world , ierror )
  if ( pre == low ) then
     call mpi_bcast( at     , 1_pin , mpi_real , 0_pin , mpi_comm_world , ierror )
     call mpi_bcast( ax     , 1_pin , mpi_real , 0_pin , mpi_comm_world , ierror )
     call mpi_bcast( ay     , 1_pin , mpi_real , 0_pin , mpi_comm_world , ierror )
     call mpi_bcast( bx     , 1_pin , mpi_real , 0_pin , mpi_comm_world , ierror )
     call mpi_bcast( by     , 1_pin , mpi_real , 0_pin , mpi_comm_world , ierror )
     call mpi_bcast( deltat , 1_pin , mpi_real , 0_pin , mpi_comm_world , ierror )
     call mpi_bcast( deltax , 1_pin , mpi_real , 0_pin , mpi_comm_world , ierror )
     call mpi_bcast( deltay , 1_pin , mpi_real , 0_pin , mpi_comm_world , ierror )
  else if ( pre == high ) then
     call mpi_bcast( at     , 1_pin , mpi_double_precision , 0_pin , mpi_comm_world , ierror )
     call mpi_bcast( ax     , 1_pin , mpi_double_precision , 0_pin , mpi_comm_world , ierror )
     call mpi_bcast( ay     , 1_pin , mpi_double_precision , 0_pin , mpi_comm_world , ierror )
     call mpi_bcast( bx     , 1_pin , mpi_double_precision , 0_pin , mpi_comm_world , ierror )
     call mpi_bcast( by     , 1_pin , mpi_double_precision , 0_pin , mpi_comm_world , ierror )
     call mpi_bcast( deltat , 1_pin , mpi_double_precision , 0_pin , mpi_comm_world , ierror )
     call mpi_bcast( deltax , 1_pin , mpi_double_precision , 0_pin , mpi_comm_world , ierror )
     call mpi_bcast( deltay , 1_pin , mpi_double_precision , 0_pin , mpi_comm_world , ierror )
  end if
  !*
  !* domain
  call domain_domain() 
  if ( rank == 0_pin ) then
     print*,''
     print*,'domain'
     do l = 1_pin , nt
        print*,'t(',l,') = ' , t( l )
     end do
     print*,''
     do i = 1_pin , nx
        print*,'x(',i,') = ' , x( i )
     end do
     print*,''
     do j = 1_pin , ny
        print*,'y(',j,') = ' , y( j )
     end do
  end if
  if ( pre == low ) then
     call mpi_bcast( t( 1 : nt ) , nt , mpi_real , 0_pin , mpi_comm_world , ierror )
     call mpi_bcast( x , nx , mpi_real , 0_pin , mpi_comm_world , ierror )
     call mpi_bcast( y , ny , mpi_real , 0_pin , mpi_comm_world , ierror )
  else if ( pre == high ) then
     call mpi_bcast( t( 1 : nt ) , nt , mpi_double_precision , 0_pin , mpi_comm_world , ierror )
     call mpi_bcast( x , nx , mpi_double_precision , 0_pin , mpi_comm_world , ierror )
     call mpi_bcast( y , ny , mpi_double_precision , 0_pin , mpi_comm_world , ierror )
  end if
  !*
  !* model parameters
  call model_parametersup()
  if ( rank == 0_pin ) then
     print*,''
     print*,'model parameters'
     print*,'errorbc = ' , errorbc
     print*,'modesmodel = ' , modesmodel
  end if
  call mpi_bcast( modesmodel , 1_pin , mpi_integer , 0_pin , mpi_comm_world , ierror )
  if ( pre == low ) then
     call mpi_bcast( errorbc , 1_pin , mpi_real , 0_pin , mpi_comm_world , ierror )
  else if ( pre == high ) then
     call mpi_bcast( errorbc , 1_pin , mpi_double_precision , 0_pin , mpi_comm_world , ierror )
  end if
  !*
  !* observations parameters
  if ( rank == 0_pin ) then
     call observations_parametersup()
     print*,''
     print*,'observations parameters'
     print*,'errorobs = ' , errorobs
     print*,'obsstep  = ' , obsstep
     print*,'no       = ' , no
     print*,'modesobs = ' , modesobs
  end if
  call mpi_bcast( obsstep  , 1_pin , mpi_integer , 0_pin , mpi_comm_world , ierror )
  call mpi_bcast( no       , 1_pin , mpi_integer , 0_pin , mpi_comm_world , ierror )
  call mpi_bcast( modesobs , 1_pin , mpi_integer , 0_pin , mpi_comm_world , ierror )
  if ( pre == low ) then
     call mpi_bcast( errorobs , 1_pin , mpi_real , 0_pin , mpi_comm_world , ierror )
  else if ( pre == high ) then
     call mpi_bcast( errorobs , 1_pin , mpi_double_precision , 0_pin , mpi_comm_world , ierror )
  end if
  !*
  !* observations time
  call observations_obstimeup()
  if ( rank == 0_pin ) then
     print*,''
     print*,'observations time'
     do lo = 1_pin , no
        print*,'to(',lo,') = ' , to( lo )
     end do
  end if
  if ( pre == low ) then
     call mpi_bcast( to , no , mpi_real , 0_pin , mpi_comm_world , ierror )
  else if ( pre == high ) then
     call mpi_bcast( to , no , mpi_double_precision , 0_pin , mpi_comm_world , ierror )
  end if
  !*
  !* stations
  call observations_stationsup() 
  if ( rank == 0_pin ) then
     print*,''
     print*,'nstat    = ' , nstat
     print*,'stations'
     do s = 1_pin , nstat
        print*,''
        print*,'station ' , s
        print*,'indexstationsx(',s,') = ' , indexstationsx( s ) 
        print*,'indexstationsy(',s,') = ' , indexstationsy( s ) 
        print*,'stationsx(',s,') = ' , stationsx( s )
        print*,'stationsy(',s,') = ' , stationsy( s )
     end do
  end if
  call mpi_bcast( nstat , 1_pin , mpi_integer , 0_pin , mpi_comm_world , ierror )
  call mpi_bcast( indexstationsx , nstat , mpi_integer , 0_pin , mpi_comm_world , ierror )
  call mpi_bcast( indexstationsy , nstat , mpi_integer , 0_pin , mpi_comm_world , ierror )
  if ( pre == low ) then
     call mpi_bcast( stationsx , nstat , mpi_real , 0_pin , mpi_comm_world , ierror )
     call mpi_bcast( stationsy , nstat , mpi_real , 0_pin , mpi_comm_world , ierror )
  else if ( pre == high ) then
     call mpi_bcast( stationsx , nstat , mpi_double_precision , 0_pin , mpi_comm_world , ierror )
     call mpi_bcast( stationsy , nstat , mpi_double_precision , 0_pin , mpi_comm_world , ierror )
  end if
  !*
  !* initialization parameters
  if ( rank == 0_pin ) then
     call initialize_parametersup()
     print*,''
     print*,'initialization parameters'
     print*,'errorin       = ' , errorin
     print*,'dimspacestate = ' , dimspacestate
     print*,'numbersamples = ' , numbersamples
     print*,'modesanalysis = ' , modesanalysis
     print*,'first         = ' , first
  end if
  call mpi_bcast( dimspacestate , 1_pin , mpi_integer , 0_pin , mpi_comm_world , ierror )
  call mpi_bcast( numbersamples , 1_pin , mpi_integer , 0_pin , mpi_comm_world , ierror )
  call mpi_bcast( modesanalysis , 1_pin , mpi_integer , 0_pin , mpi_comm_world , ierror )
  call mpi_bcast( first         , 1_pin , mpi_logical , 0_pin , mpi_comm_world , ierror )
  if ( pre == low ) then
     call mpi_bcast( errorin , 1_pin , mpi_real , 0_pin , mpi_comm_world , ierror )
  else if ( pre == high ) then
     call mpi_bcast( errorin , 1_pin , mpi_double_precision , 0_pin , mpi_comm_world , ierror )
  end if
  !*
  !* open files
  if ( rank == 0_pin ) then
     open( unit = unittruth , file = filetruth )
     open( unit = unitmodel , file = filemodel )
     open( unit = unitobser , file = fileobser )
     open( unit = unitassimekf , file = fileassimekf )
  end if
  !*
  !* allocations
  allocate( state( nx , ny ) )
  allocate( stateekf( nx , ny ) )
  allocate( covariance( dimspacestate , dimspacestate ) )
  allocate( covarianceekf( dimspacestate , dimspacestate ) )
  allocate( truth( nt , nx , ny ) )
  allocate( modl( nt , nx , ny ) )
  allocate( obser( no , nstat ) )
  allocate( assimekf( nt , nx , ny ) )
  !*
  !***************************************************** step 1 ******************************************************************
  !*
  if ( rank == 0_pin ) then
     print*,''
     print*,'resolution'
  end if
  !*
  !* truth
  do i = 1_pin , nx
     do j = 1_pin , ny
        truth( 1_pin , i , j ) = solution_solution( t( 1_pin ) , x( i ) , y( j ) )
     end do
  end do
  !*
  !* initialization
  if ( rank == 0_pin ) then
     call initialize_initialize1( state , covariance )
  end if
  if ( pre == low ) then
     call mpi_bcast( state , dimspacestate , mpi_real , 0_pin , mpi_comm_world , ierror )
     call mpi_bcast( covariance , dimspacestate * dimspacestate , mpi_real , 0_pin , mpi_comm_world , ierror )
  else if ( pre == high ) then
     call mpi_bcast( state , dimspacestate , mpi_double_precision , 0_pin , mpi_comm_world , ierror )
     call mpi_bcast( covariance , dimspacestate * dimspacestate , mpi_double_precision , 0_pin , mpi_comm_world , ierror )
  end if
  !*
  stateekf = state
  covarianceekf = covariance
  !*
  !* model
  modl( 1_pin , 1_pin : nx , 1_pin : ny ) = state
  !*
  !* assimilation
  assimekf( 1_pin , 1_pin : nx , 1_pin : ny ) = stateekf( 1_pin : nx , 1_pin : ny )
  !*
  !* writing to file
  if ( rank == 0_pin ) then
     do i = 1_pin , nx
        do j = 1_pin , ny
           write( unittruth , '( f15.8 , f15.8 , f15.8 , f15.8 )' ) t( 1_pin ) , x( i ) , y( j ) , truth( 1_pin , i , j ) 
           write( unitmodel , '( f15.8 , f15.8 , f15.8 , f15.8 )' ) t( 1_pin ) , x( i ) , y( j ) , modl( 1_pin , i , j ) 
           write( unitassimekf , '( f15.8 , f15.8 , f15.8 , f15.8 )' ) t( 1_pin ) , x( i ) , y( j ) , assimekf( 1_pin , i , j )
        end do
     end do
  end if
  !*
  !* printing on screen
  if ( rank == 0_pin ) then
     write( * , '("step: ",i5)' ) 1_pin 
  end if
  !*
  !***************************************************** step l+1 ****************************************************************
  !*
  lo = 1_pin
  !*
  do l = 1_pin , nt - 1_pin
     !*
     !* truth
     if ( rank == 0_pin ) then
        do i = 1_pin , nx
           do j = 1_pin , ny
              truth( l + 1_pin , i , j ) = solution_solution( t( l + 1_pin ) , x( i ) , y( j ) )
           end do
        end do
     end if
     if ( pre == low ) then
        call mpi_bcast( truth( l + 1_pin , 1 : nx , 1 : ny ) , nx * ny , mpi_real , mpi_comm_world , ierror )
     else if ( pre == high ) then
        call mpi_bcast( truth( l + 1_pin , 1 : nx , 1 : ny ) , nx * ny , mpi_double_precision , 0_pin , mpi_comm_world , ierror )
     end if
     !*
     !* boundary conditions at time l + 1
     if ( rank == 0_pin ) then
        call model_bc( l + 1_pin )
     end if
     if ( pre == low ) then
        call mpi_bcast( sigmabc( 1 : nx , 1 : ny ) , nx * ny , mpi_real , mpi_comm_world , ierror )
        call mpi_bcast( bc( 1 : nx , 1 : ny ) , nx * ny , mpi_real , mpi_comm_world , ierror )
     else if ( pre == high ) then
        call mpi_bcast( sigmabc( 1 : nx , 1 : ny ) , nx * ny , mpi_double_precision , 0_pin , mpi_comm_world , ierror )
        call mpi_bcast( bc( 1 : nx , 1 : ny ) , nx * ny , mpi_double_precision , 0_pin , mpi_comm_world , ierror )
     end if
     !*
     !* model
     if ( rank == 0_pin ) then
        call model_model( l , modl( l , 1_pin : nx , 1_pin : ny ) , modl( l + 1_pin , 1_pin : nx , 1_pin : ny ) )
     end if
     if ( pre == low ) then
        call mpi_bcast( modl( l + 1_pin , 1 : nx , 1 : ny ) , nx * ny , mpi_real , mpi_comm_world , ierror )
     else if ( pre == high ) then
        call mpi_bcast( modl( l + 1_pin , 1 : nx , 1 : ny ) , nx * ny , mpi_double_precision , 0_pin , mpi_comm_world , ierror )
     end if
     !*
     !* assimilation - prediction
     call ekf_predictor( l , stateekf , covarianceekf )
     !*
     !* assimilation - correction
     if ( rank == 0_pin ) then
        call observations_ifobservations( l + 1_pin )
     end if
     call mpi_bcast( ifobs , 1_pin , mpi_logical , 0_pin , mpi_comm_world , ierror )
     !*
     if ( ifobs ) then
        !*
        if ( rank == 0_pin ) then
           call observations_numberobs( l + 1_pin )
        end if
        call mpi_bcast( numberobs , 1_pin , mpi_integer , 0_pin , mpi_comm_world , ierror )
        !*
        allocate( sigmaobs( numberobs ) )
        allocate( obsvalue( numberobs ) )
        !*
        if ( rank == 0_pin ) then
           call observations_obsvalue( l + 1_pin )
        end if
        if ( pre == low ) then
           call mpi_bcast( sigmaobs , numberobs , mpi_real , 0_pin , mpi_comm_world , ierror )
           call mpi_bcast( obsvalue , numberobs , mpi_real , 0_pin , mpi_comm_world , ierror )
        else if ( pre == high ) then
           call mpi_bcast( sigmaobs , numberobs , mpi_double_precision , 0_pin , mpi_comm_world , ierror )
           call mpi_bcast( obsvalue , numberobs , mpi_double_precision , 0_pin , mpi_comm_world , ierror )
        end if
        !*
        obser( lo , 1_pin : nstat ) = obsvalue
        !*
        if ( rank == 0_pin ) then
           do s = 1_pin , nstat
              write( unitobser , '( f15.8 , f15.8 , f15.8 , f15.8 )' ) to( lo ) , stationsx( s ) , stationsy( s ) , obser( lo , s&
                   & )
           end do
        end if
        !*
        call ekf_corrector( l , stateekf , covarianceekf )
        !*
        deallocate( obsvalue )
        deallocate( sigmaobs )
        !*
        lo = lo + 1_pin
        !*
     end if
     !*
     assimekf( l + 1_pin , 1_pin : nx , 1_pin : ny ) = stateekf
     !*
     !* writing to file
     if ( rank == 0_pin ) then
        do i = 1_pin , nx
           do j = 1_pin , ny
              write( unittruth , '( f15.8 , f15.8 , f15.8 , f15.8 )' ) t( l + 1_pin ) , x( i ) , y( j ) , truth( l + 1_pin , i , &
                   &j ) 
              write( unitmodel , '( f15.8 , f15.8 , f15.8 , f15.8 )' ) t( l + 1_pin ) , x( i ) , y( j ) , modl( l + 1_pin , i , j&
                   & ) 
              write( unitassimekf , '( f15.8 , f15.8 , f15.8 , f15.8 )' ) t( l + 1_pin ) , x( i ) , y( j ) , assimekf( l + 1_pin &
                   &, i , j ) 
           end do
        end do
     end if
     !*
     !* printing on screen
     if ( rank == 0_pin ) then
        if ( ifobs ) then
           write( * , '("step: ",i5," --> assimilation")' ) l + 1_pin
        else
           write( * , '("step: ",i5)' ) l + 1_pin 
        end if
     end if
     !*
  end do
  !*
  !*******************************************************************************************************************************
  !*
  !* deallocations
  call domain_destructor()
  call model_destructor()
  call observations_destructor()
  deallocate( state )
  deallocate( stateekf )
  deallocate( covariance )
  deallocate( covarianceekf )
  deallocate( truth )
  deallocate( modl )
  deallocate( obser )
  deallocate( assimekf )
  !*
  !* closing files
  if ( rank == 0_pin ) then
     close( unittruth )
     close( unitmodel )
     close( unitobser )
     close( unitassimekf )
  end if
  !*
  !* end parallelization
  call parallel_finalize()
  !*
end program example2d_ekf
