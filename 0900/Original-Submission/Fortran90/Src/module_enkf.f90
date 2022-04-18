!*********************************************************************************************************************************
!*
!* MODULE: ENKF
!*
!* PURPOSE: implementation of the Ensemble Kalman Filter
!*
!* DEPENDENCIES:
!*               -  PRECISION
!*               -  RANDOM
!*               -  MODEL
!*               -  OBSERVATIONS
!*               -  PARALLEL
!*               -  INITIALIZE
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
module enkf
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
  !* FUNCTION/SUBROUTINE: ENKF_PREDICTOR
  !*
  !* PURPOSE: makes the prediction step of the Ensemble Kalman Filter, that is:
  !*
  !*          ONLY THE FIRST TIME: Generate x_1^a , ... , x_numbersamples^a random vectors such that belongs to N( analysis(l) , &
  !*                               &covanalysis(l) ) 
  !*          Propagate: x_i^f = model_{l -> l+1} ( x_i^a ) + eta_i, eta_i belongs to N( 0 , modelerror(l) ),  i=1,...,numbersamp&
  !*                                                                                                             &les
  !*          forecast(l+1) = ( x_1^f + ... + x_numbersamples^f ) / numbersamples
  !*          covforecast(l+1) = [ ( x_1^f - forecast(l+1) ) * ( x_1^f - forecast(l+1) )^T + ... + ( x_numbersamples^f - forecast&
  !*                             &(l+1) ) * ( x_numbersamples^f - forecast(l+1) )^T ] / ( numbersamples - 1 )
  !*
  !* INPUTS:
  !*         - L (integer): current time step
  !*
  !* INPUTS/OUTPUTS:
  !*                 - STATE (real array of size DIMSPACESTATE): on input the analysis at time L, on output the forecast at time &
  !*                                                             &L+1
  !*                 - COVARIANCE (real array of size DIMSPACESTATE x DIMSPACESTATE): on input the covariance matrix of analysis &
  !*                                                                                  &errors at time L, on output the covariance&
  !*                                                                                  & matrix of forecast errors at time L+1
  !*                 - SAMPLES (real array of size DIMSPACESTATE x NUMBERSAMPLES): on input the members of the ensemble represent&
  !*                                                                               &ing posible states of analysis at time L, on &
  !*                                                                               &output the members of the ensemble representi&
  !*                                                                               &ng posible states of forecast at time L+1
  !*
  !* OUTPUTS: none
  !*
  !* CALLS:
  !*        - RANDOM_NORMAL1DMATRIX (from RANDOM)
  !*        - MODEL_MODELERROR (from MODEL)
  !*        - PARALLEL_OPERATOR (from OPERATOR)
  !*        - RANDOM_MEANESTIMATOR1D (from RANDOM)
  !*        - RANDOM_COVARIANCEESTIMATOR1D (from RANDOM)
  !* 
  !* COMMENTS: none
  !*
  !*******************************************************************************************************************************
  subroutine enkf_predictor( l , state , covariance , samples )
    implicit none
    integer( kind = pin ) , intent( in ) :: l
    real( kind = pre ) , dimension( dimspacestate ) , intent( inout ) :: state
    real( kind = pre ) , dimension( dimspacestate , dimspacestate ) , intent( inout ) :: covariance
    real( kind = pre ) , dimension( dimspacestate , numbersamples ) , intent( inout ) :: samples
    character( kind = pch , len = 100_pin ) :: option
    real( kind = pre ) , allocatable , dimension( : ) :: vectoraux1
    real( kind = pre ) , allocatable , dimension( : , : ) :: samplesmodel
    real( kind = pre ) , allocatable , dimension( : , : ) :: matrixaux1
    if ( first ) then
       call random_normal1dmatrix( state , covariance , dimspacestate , numbersamples , samples )
       first = .false.  
    end if
    allocate( modelerror( dimspacestate , dimspacestate ) )
    call model_modelerror( l )
    allocate( samplesmodel( dimspacestate , numbersamples ) )
    allocate( vectoraux1( dimspacestate ) )
    vectoraux1 = 0.0_pre
    call random_normal1dmatrix( vectoraux1 , modelerror , dimspacestate , numbersamples , samplesmodel )
    deallocate( vectoraux1 )
    deallocate( modelerror )
    allocate( matrixaux1( dimspacestate , numbersamples ) )
    option = 'modelN'
    call parallel_operator( l , samples , matrixaux1 , option )
    samples = matrixaux1 + samplesmodel
    deallocate( matrixaux1 )
    deallocate( samplesmodel )
    call random_meanestimator1d( dimspacestate , numbersamples , samples , state )
    call random_covarianceestimator1d( dimspacestate , numbersamples , samples , state , covariance )
  end subroutine enkf_predictor
  !*******************************************************************************************************************************
  !* 
  !* FUNCTION/SUBROUTINE: ENKF_CORRECTOR
  !*
  !* PURPOSE: makes the correction step of the Ensemble Kalman Filter, that is:
  !* 
  !*            gainmatrix(l+1) = covforecast(l+1) * tangobsop(l+1) * [ tangobsop(l+1) * covforecast(l+1) * tangobsop(l+1)^T + co&
  !*                              &vobs(l+1) ]^{-1}
  !*            x_i^a = x_i^f + gainmatrix(l+1) * [ y_i - obsop(l+1) ( x_i^f ) ], y_i belongs to N( obsvalue(l+1) , covobs(l+1) )&
  !*                    &, i=1,...,numbersamples  
  !*            analysis(l+1) = ( x_1^a + ... + x_numbersamples^a ) / numbersamples
  !*            covanalysis(l+1) = [ ( x_1^a - analysis(l+1) ) * ( x_1^a - analysis(l+1) )^T + ... + ( x_numbersamples^a - analys&
  !*                               &is(l+1) ) * ( x_numbersamples^a - analysis(l+1) )^T ] / ( numbersamples - 1 )
  !*
  !* INPUTS:
  !*         - L (integer): current time step
  !*
  !* INPUTS/OUTPUTS:
  !*                 - STATE (real array of size DIMSPACESTATE): on input the forecast at time L+1, on exit the analysis at time &
  !*                                                             &L+1
  !*                 - COVARIANCE (real array of size DIMSPACESTATE x DIMSPACESTATE): on input the covariance matrix of forecast &
  !*                                                                                  &errors at time L+1, on output the covarian&
  !*                                                                                  &ce matrix of analysis errors at time L+1
  !*                 - SAMPLES (real array of size DIMSPACESTATE x NUMBERSAMPLES): on input the members of the ensemble represent&
  !*                                                                               &ing posible states of forecast at time L+1, o&
  !*                                                                               &n output the members of the ensemble represen&
  !*                                                                               &ting posible states of analysis at time L+1
  !*
  !* OUTPUTS: none
  !*
  !* CALLS:
  !*        - OBSERVATIONS_COVOBS: (from OBSERVATIONS)
  !*        - PARALLEL_OPERATOR: (from PARALLEL)
  !*        - PARALLEL_GESV: (from PARALLEL)
  !*        - PARALLEL_GEMM: (from PARALLEL)
  !*        - RANDOM_NORMAL1DMATRIX: (from RANDOM)
  !*        - RANDOM_MEANESTIMATOR1D: (from RANDOM)
  !*        - RANDOM_COVARIANCEESTIMATOR1D: (from RANDOM)
  !*
  !* COMMENTS:
  !*           - NUMBEROBS: it is the number of observations at time L+1
  !*           - OBSVALUE: vector of observations at time L+1
  !*           - COVOBS: it sets the covariance matrix of observations errors at time L+1 during execution, and after that it is &
  !*                     &destroyed
  !*
  !*******************************************************************************************************************************
  subroutine enkf_corrector( l , state , covariance , samples )
    implicit none
    integer( kind = pin ) , intent( in ) :: l
    real( kind = pre ) , dimension( dimspacestate ) , intent( inout ) :: state
    real( kind = pre ) , dimension( dimspacestate , dimspacestate ) , intent( inout ) :: covariance
    real( kind = pre ) , dimension( dimspacestate , numbersamples ) , intent( inout ) :: samples
    character( kind = pch , len = 100_pin ) :: option
    real( kind = pre ) , allocatable , dimension( : , : ) :: matrixaux2
    real( kind = pre ) , allocatable , dimension( : , : ) :: matrixaux3
    real( kind = pre ) , allocatable , dimension( : , : ) :: matrixaux4
    real( kind = pre ) , allocatable , dimension( : , : ) :: samplesobs
    allocate( covobs( numberobs , numberobs ) )
    call observations_covobs( l + 1_pin )
    allocate( matrixaux2( numberobs , dimspacestate ) )
    option = 'tangobsopT'
    call parallel_operator( l + 1_pin , covariance , matrixaux2 , option )
    allocate( matrixaux3( numberobs , numberobs ) )
    option = 'tangobsopT'
    call parallel_operator( l + 1_pin , matrixaux2 , matrixaux3 , option )
    matrixaux3 = matrixaux3 + covobs
    allocate( samplesobs( numberobs , numbersamples ) )
    call random_normal1dmatrix( obsvalue , covobs , numberobs , numbersamples , samplesobs )
    deallocate( covobs )
    allocate( matrixaux4( numberobs , numbersamples ) )
    option = 'obsopN'
    call parallel_operator( l + 1_pin , samples , matrixaux4 , option )
    matrixaux4 = samplesobs - matrixaux4
    deallocate( samplesobs )
    call parallel_gesv( numberobs , matrixaux3 , numbersamples , matrixaux4 )
    deallocate( matrixaux3 )
    call parallel_gemm( 'T' , 'N' , numberobs , dimspacestate , matrixaux2 , numberobs , numbersamples , matrixaux4 , dimspacesta&
         &te , numbersamples , samples , 1.0_pre , 1.0_pre )
    deallocate( matrixaux4 )
    deallocate( matrixaux2 )
    call random_meanestimator1d( dimspacestate , numbersamples , samples , state )
    call random_covarianceestimator1d( dimspacestate , numbersamples , samples , state , covariance )
  end subroutine enkf_corrector
  !*
end module enkf
