!*********************************************************************************************************************************
!*
!* MODULE: RRSQRTENKF
!*
!* PURPOSE: implementation of the Reduced Rank Square Root Ensemble Kalman Filter
!*
!* DEPENDENCIES:
!*               - PRECISION
!*               - RANDOM
!*               - MODEL
!*               - OBSERVATIONS
!*               - PARALLEL
!*               - INITIALIZE
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
!*          - MODESMODEL (integer): number of model modes (from MODEL)
!*          - SQRTMODELERROR (real array): square root of covariance matrix of model errors (from MODEL)
!*          - MODELERROR (real array): covariance matrix of model errors (from MODEL)
!*
!*          - MODESOBS (integer): number of observation modes (from OBSERVATIONS)
!*          - NUMBEROBS (integer): number of observations (from OBSERVATIONS)
!*          - IFOBS (logical): flag to determine if there are observations (from OBSERVATIONS)
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
!*
!* SCOPE: this module belongs to the library
!* 
!* AUTHOR: Germán Ariel Torres
!*         Facultad de Matemática, Astronomía y Física
!*         Universidad Nacional de Córdoba
!*         Argentina
!*
!* EMAIL: torres@famaf.unc.edu.ar
!*
!*********************************************************************************************************************************
module rrsqrtenkf
  !*
  use precision
  use random
  use model
  use observations
  use parallel
  use initialize
  !*
  implicit none
  !*
contains
  !*******************************************************************************************************************************
  !* 
  !* FUNCTION/SUBROUTINE: RRSQRTENKF_PREDICTOR
  !*
  !* PURPOSE: makes the prediction step of the Reduced Rank Square Root Ensemble Kalman Filter, that is:
  !*
  !*          ONLY THE FIRST TIME: Generate x_1^a , ... , x_numbersamples^a random vectors such that belongs to N( analysis(l) , &
  !*                               &sqrtcovanalysis(l) * sqrtcovanalysis(l)^T ) 
  !*          Propagate: x_i^f = model_{l -> l+1} ( x_i^a ) + eta_i, eta_i belongs to N( 0 , sqrtmodelerror(l) * sqrtmodelerror(l&
  !*                             &)^T ),  i=1,...,numbersamples
  !*          forecast(l+1) = ( x_1^f + ... + x_numbersamples^f ) / numbersamples
  !*          sqrtcovforecast(l+1) = [ ( x_1^f - forecast(l+1) ) ... ( x_numbersamples^f - forecast(l+1) ) ] / sqrt( numbersample&
  !*                                 &s - 1 )
  !*
  !* INPUTS:
  !*         - L (integer): current time step
  !*         - SQRTCOV (real array of size DIMSPACESTATE x MODESANALYSIS): reduced rank square root of the covariance matrix of a&
  !*                                                                       &nalysis errors at time L
  !*
  !* INPUTS/OUTPUTS:
  !*                 - STATE (real array of size DIMSPACESTATE): on input the analysis at time L, on output the forecast at time &
  !*                                                             &L+1
  !*                 - SAMPLES (real array of size DIMSPACESTATE x NUMBERSAMPLES): on input the members of the ensemble represent&
  !*                                                                               &ing posible states of analysis at time L, on &
  !*                                                                               &output the members of the ensemble representi&
  !*                                                                               &ng posible states of forecast at time L+1
  !*
  !* OUTPUTS:
  !*          - SQRTCOVAUX (real array of size DIMSPACESTATE x NUMBERSAMPLES): on output the reduced rank square root of the cova&
  !*                                                                           &riance matrix of forecast errors at time L+1
  !*
  !* CALLS:
  !*        - RANDOM_NORMAL1DSQRTMATRIX (from RANDOM)
  !*        - RANDOM_MEANESTIMATOR1D (from RANDOM)
  !*        - RANDOM_SQRTCOVESTIMATOR1D (from RANDOM)
  !*        - PARALLEL_OPERATOR (from PARALLEL)
  !*        - MODEL_SQRTMODELERROR (from MODEL)
  !*  
  !* COMMENTS: none
  !*
  !*******************************************************************************************************************************
  subroutine rrsqrtenkf_predictor( l , state , sqrtcov , samples , sqrtcovaux )
    implicit none
    integer( kind = pin ) , intent( in ) :: l
    real( kind = pre ) , dimension( dimspacestate ) , intent( inout ) :: state
    real( kind = pre ) , dimension( dimspacestate , modesanalysis ) , intent( in ) :: sqrtcov
    real( kind = pre ) , dimension( dimspacestate , numbersamples ) , intent( inout ) :: samples
    real( kind = pre ) , dimension( dimspacestate , numbersamples ) , intent( out ) :: sqrtcovaux
    character( kind = pch , len = 100_pin ) :: option
    real( kind = pre ) , allocatable , dimension( : ) :: vectoraux1
    real( kind = pre ) , allocatable , dimension( : , : ) :: samplesmodel
    real( kind = pre ) , allocatable , dimension( : , : ) :: matrixaux1
    if ( first ) then
       call random_normal1dsqrtmatrix( state , sqrtcov , dimspacestate , modesanalysis , numbersamples , samples )
       first = .false.
    end if
    allocate( sqrtmodelerror( dimspacestate , modesmodel ) )
    call model_sqrtmodelerror( l )
    allocate( samplesmodel( dimspacestate , numbersamples ) )
    allocate( vectoraux1( dimspacestate ) )
    vectoraux1 = 0.0_pre
    call random_normal1dsqrtmatrix( vectoraux1 , sqrtmodelerror , dimspacestate , modesmodel , numbersamples , samplesmodel )
    deallocate( vectoraux1 )
    deallocate( sqrtmodelerror )
    allocate( matrixaux1( dimspacestate , numbersamples ) )
    option = 'modelN'
    call parallel_operator( l , samples , matrixaux1 , option )
    samples = matrixaux1 + samplesmodel
    deallocate( matrixaux1 )
    deallocate( samplesmodel )
    call random_meanestimator1d( dimspacestate , numbersamples , samples , state )
    call random_sqrtcovestimator1d( dimspacestate , numbersamples , samples , state , sqrtcovaux )
  end subroutine rrsqrtenkf_predictor
  !*******************************************************************************************************************************
  !* 
  !* FUNCTION/SUBROUTINE: RRSQRTENKF_CORRECTOR
  !*
  !* PURPOSE: makes the correction step of the Reduced Rank Square Root Ensemble Kalman Filter, that is:
  !*
  !*          gainmatrix(l+1) = sqrtcovforecast(l+1) * ( tangobsop(l+1) * sqrtcovforecast(l+1)  )^T * ( [ tangobsop(l+1) * sqrtco&
  !*                            &vforecast(l+1) | sqrtcovobs(l+1) ] * [ tangobsop(l+1) * sqrtcovforecast(l+1) | sqrtcovobs(l+1) ]&
  !*                            &^T  )^{-1}
  !*          x_i^a = x_i^f + gainmatrix(l+1) * [ y_i - obsop(l+1) ( x_i^f ) ], y_i belongs to N( obsvalue(l+1) , sqrtcovobs(l+1)&
  !*                  & * sqrtcovobs(l+1)^T ), i=1,...,numbersamples  
  !*          analysis(l+1) = ( x_1^a + ... + x_numbersamples^a ) / numbersamples
  !*          sqrtcovanalysis(l+1) = [ ( x_1^a - analysis(l+1) ) ... ( x_numbersamples^a - analysis(l+1) ) ] / sqrt( numbersample&
  !*                                 &s - 1 )
  !*          sqrtcovanalysis(l+1) = cut sqrtcovanalysis(l+1) to modesanalysis columns using SVD if necessary
  !*
  !* INPUTS:
  !*         - L (integer): current time step
  !*
  !* INPUTS/OUTPUTS:
  !*                 - STATE (real array of size DIMSPACESTATE): on input the forecast at time L+1, on exit the analysis at time &
  !*                                                             &L+1
  !*                 - SAMPLES (real array of size DIMSPACESTATE x NUMBERSAMPLES): on input the members of the ensemble represent&
  !*                                                                               &ing posible states of forecast at time L+1, o&
  !*                                                                               &n output the members of the ensemble represen&
  !*                                                                               &ting posible states of analysis at time L+1
  !*                 - SQRTCOVAUX (real array of size DIMSPACESTATE x NUMBERSAMPLES): on input the reduced rank square root of th&
  !*                                                                                  &e covariance matrix of forecast errors at &
  !*                                                                                  &time L+1, on output the reduced rank squar&
  !*                                                                                  &e root of the covariance matrix of analysi&
  !*                                                                                  &s errors at time L+1 without reducing rank
  !*
  !* OUTPUTS:
  !*          - SQRTCOV (real array of size DIMSPACESTATE x MODESANALYSIS): reduced rank square root of the covariance matrix of &
  !*                                                                        &analysis errors at time L+1
  !*
  !* CALLS:
  !*        - PARALLEL_OPERATOR (from PARALLEL)
  !*        - OBSERVATIONS_SQRTCOVOBS: (from OBSERVATIONS)
  !*        - RANDOM_NORMAL1DSQRTMATRIX: (from RANDOM)
  !*        - PARALLEL_GEMM (from PARALLEL)
  !*        - PARALLEL_GESV (from PARALLEL)
  !*        - RANDOM_MEANESTIMATOR1D: (from RANDOM)
  !*        - RANDOM_SQRTCOVESTIMATOR1D: (from RANDOM)
  !*        - PARALLEL_SVDCUT (from PARALLEL)
  !*
  !* COMMENTS:
  !*           - NUMBEROBS: it is the number of observations at time L+1
  !*           - OBSVALUE: vector of observations at time L+1
  !*           - SQRTCOVOBS: it sets the reduced rank square root of the covariance matrix of observations errors at time L+1 dur&
  !*                         &ing execution, and after that it is destroyed
  !*
  !*******************************************************************************************************************************
  subroutine rrsqrtenkf_corrector( l , state , sqrtcov , samples , sqrtcovaux )
    implicit none
    integer( kind = pin ) , intent( in ) :: l
    real( kind = pre ) , dimension( dimspacestate ) , intent( inout ) :: state
    real( kind = pre ) , dimension( dimspacestate , modesanalysis ) , intent( out ) :: sqrtcov
    real( kind = pre ) , dimension( dimspacestate , numbersamples ) , intent( inout ) :: samples
    real( kind = pre ) , dimension( dimspacestate , numbersamples ) , intent( inout ) :: sqrtcovaux
    character( kind = pch , len = 100_pin ) :: option
    real( kind = pre ) , allocatable , dimension( : , : ) :: matrixaux2
    real( kind = pre ) , allocatable , dimension( : , : ) :: matrixaux3
    real( kind = pre ) , allocatable , dimension( : , : ) :: matrixaux4
    real( kind = pre ) , allocatable , dimension( : , : ) :: matrixaux5
    real( kind = pre ) , allocatable , dimension( : , : ) :: matrixaux6
    real( kind = pre ) , allocatable , dimension( : , : ) :: samplesobs
    allocate( matrixaux2( numberobs , numbersamples ) )
    option = 'tangobsopN'
    call parallel_operator( l + 1_pin , sqrtcovaux , matrixaux2 , option )
    allocate( sqrtcovobs( numberobs , modesobs ) )
    call observations_sqrtcovobs( l + 1_pin )
    allocate( samplesobs( numberobs , numbersamples ) )
    call random_normal1dsqrtmatrix( obsvalue , sqrtcovobs , numberobs , modesobs , numbersamples , samplesobs )
    allocate( matrixaux3( numberobs , numbersamples ) )
    option = 'obsopN'
    call parallel_operator( l + 1_pin , samples , matrixaux3 , option )
    matrixaux3 = samplesobs - matrixaux3
    deallocate( samplesobs )
    allocate( matrixaux4( numberobs , numbersamples + modesobs ) )
    matrixaux4( 1_pin : numberobs , 1_pin : numbersamples ) = matrixaux2( 1_pin : numberobs , 1_pin : numbersamples )
    matrixaux4( 1_pin : numberobs , numbersamples + 1_pin : numbersamples + modesobs ) = sqrtcovobs( 1_pin : numberobs , 1_pin : &
         &modesobs )
    deallocate( sqrtcovobs )
    allocate( matrixaux5( numberobs , numberobs ) )
    call parallel_gemm( 'N' , 'T' , numberobs , numbersamples + modesobs , matrixaux4 , numberobs , numbersamples + modesobs , ma&
         &trixaux4 , numberobs , numberobs , matrixaux5 , 1.0_pre , 0.0_pre )
    deallocate( matrixaux4 )
    call parallel_gesv( numberobs , matrixaux5 , numbersamples , matrixaux3 )
    deallocate( matrixaux5 )
    allocate( matrixaux6( numbersamples , numbersamples ) )
    call parallel_gemm( 'T' , 'N' , numberobs , numbersamples , matrixaux2 , numberobs , numbersamples , matrixaux3 , numbersampl&
         &es , numbersamples , matrixaux6 , 1.0_pre , 0.0_pre )
    deallocate( matrixaux2 )
    deallocate( matrixaux3 )
    call parallel_gemm( 'N' , 'N' , dimspacestate , numbersamples , sqrtcovaux , numbersamples , numbersamples , matrixaux6 , dim&
         &spacestate , numbersamples , samples , 1.0_pre , 1.0_pre )
    deallocate( matrixaux6 )
    call random_meanestimator1d( dimspacestate , numbersamples , samples , state )
    call random_sqrtcovestimator1d( dimspacestate , numbersamples , samples , state , sqrtcovaux )
    if ( numbersamples > modesanalysis ) then
       call parallel_svdcut( dimspacestate , numbersamples , modesanalysis , sqrtcovaux , sqrtcov )
    end if
  end subroutine rrsqrtenkf_corrector
  !*
end module rrsqrtenkf
