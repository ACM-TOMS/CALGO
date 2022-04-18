!*********************************************************************************************************************************
!*
!* MODULE: RRSQRTKF
!*
!* PURPOSE: implementation of the Reduced Rank Square Root Kalman Filter
!*
!* DEPENDENCIES:
!*               - PRECISION
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
module rrsqrtkf
  !*
  use precision
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
  !* FUNCTION/SUBROUTINE: RRSQRTKF_PREDICTOR
  !*
  !* PURPOSE: makes the prediction step of the Reduced Rank Square Root Kalman Filter, that is:
  !*
  !*          forecast(l+1) = model_{l -> l+1} ( analysis(l) )
  !*          sqrtcovforecast(l+1) = [ tangmodel_{l -> l+1} sqrtcovanalysis(l) | sqrtmodelerror(l) ]
  !*
  !* INPUTS:
  !*         - L (integer): current time step
  !*         - SQRTCOV (real array of size DIMSPACESTATE x MODESANALYSIS): reduced rank square root of the covariance matrix of a&
  !*                                                                       &nalysis errors at time L
  !*
  !* INPUTS/OUTPUTS:
  !*                 - STATE (real array of size DIMSPACESTATE): on input the analysis at time L, on output the forecast at time &
  !*                                                             &L+1
  !*
  !* OUTPUTS:
  !*          - SQRTCOVFOR (real array of size DIMSPACESTATE x ( MODESANALYSIS + MODESMODEL ) ): on output the reduced rank squar&
  !*                                                                                             &e root of the covariance matrix&
  !*                                                                                             & of forecast errors at time L+1
  !*
  !* CALLS:
  !*        - MODEL_MODEL (from MODEL)
  !*        - MODEL_SQRTMODELERROR (from MODEL)
  !*        - PARALLEL_OPERATOR (from PARALLEL)
  !*  
  !* COMMENTS: none
  !*
  !*******************************************************************************************************************************
  subroutine rrsqrtkf_predictor( l , state , sqrtcov , sqrtcovfor )
    implicit none
    integer( kind = pin ) , intent( in ) :: l
    real( kind = pre ) , dimension( dimspacestate ) , intent( inout ) :: state
    real( kind = pre ) , dimension( dimspacestate , modesanalysis ) , intent( in ) :: sqrtcov
    real( kind = pre ) , dimension( dimspacestate , modesanalysis + modesmodel ) , intent( out ) :: sqrtcovfor
    character( kind = pch , len = 100_pin ) :: option
    real( kind = pre ) , allocatable , dimension( : ) :: stateaux
    allocate( stateaux( dimspacestate ) )
    call model_model( l , state , stateaux )
    state = stateaux
    deallocate( stateaux )
    allocate( sqrtmodelerror( dimspacestate , modesmodel ) )
    call model_sqrtmodelerror( l )
    option = 'tangmodelN'
    call parallel_operator( l , sqrtcov( 1_pin : dimspacestate , 1_pin : modesanalysis ) , sqrtcovfor( 1_pin : dimspacestate , 1_&
         &pin : modesanalysis ) , option )
    sqrtcovfor( 1_pin : dimspacestate , modesanalysis + 1_pin : modesanalysis + modesmodel ) = sqrtmodelerror( 1_pin : dimspacest&
         &ate , 1_pin : modesmodel )
    deallocate( sqrtmodelerror )
  end subroutine rrsqrtkf_predictor
  !*******************************************************************************************************************************
  !* 
  !* FUNCTION/SUBROUTINE: RRSQRTKF_CORRECTOR
  !*
  !* PURPOSE: makes the correction step of the Reduced Rank Square Root Kalman Filter, that is:
  !*
  !*          gainmatrix(l+1) = sqrtcovforecast(l+1) * ( tangobsop(l+1) * sqrtcovforecast(l+1)  )^T * ( [ tangobsop(l+1) * sqrtco&
  !*                            &vforecast(l+1) | sqrtcovobs(l+1) ] * [ tangobsop(l+1) * sqrtcovforecast(l+1) | sqrtcovobs(l+1) ]&
  !*                            &^T  )^{-1}
  !*          analysis(l+1) = forecast(l+1) + gainmatrix(l+1) * [ obsvalue(l+1) - obsop(l+1)( forecast(l+1) ) ]
  !*          sqrtcovanalysis(l+1) = [ ( I - gainmatrix(l+1) * tangobsop(l+1) ) * sqrtcovforecast(l+1) | gainmatrix(l+1) * sqrtco&
  !*                                 &vobs(l+1) ]
  !*          sqrtcovanalysis(l+1) = reduce sqrtcovanalysis(l+1) to MODESANALYSIS columns
  !*
  !* INPUTS:
  !*         - L (integer): current time step
  !*         - SQRTCOVFOR (real array of size DIMSPACESTATE x (MODESANALYSIS + MODESMODEL)): on input the reduced rank square roo&
  !*                                                                                         &t of the covariance matrix of forec&
  !*                                                                                         &ast errors at time L+1, on output t&
  !*                                                                                         &he reduced rank square root of the &
  !*                                                                                         &covariance matrix of analysis error&
  !*                                                                                         &s at time L+1 without reducing rank
  !* INPUTS/OUTPUTS:
  !*                 - STATE (real array of size DIMSPACESTATE): on input the forecast at time L+1, on exit the analysis at time &
  !*                                                             &L+1
  !*
  !* OUTPUTS:
  !*          - SQRTCOV (real array of size DIMSPACESTATE x MODESANALYSIS): reduced rank square root of the covariance matrix of &
  !*                                                                        &analysis errors at time L+1
  !*
  !* CALLS:
  !*        - OBSERVATIONS_SQRTCOVOBS: (from OBSERVATIONS)
  !*        - OBSERVATIONS_OBSOP: (from OBSERVATIONS)
  !*        - PARALLEL_OPERATOR: (from PARALLEL)
  !*        - PARALLEL_GEMM: (from PARALLEL)
  !*        - PARALLEL_GESV: (from PARALLEL)
  !*        - PARALLEL_SVDCUT: (from PARALLEL)
  !*
  !* COMMENTS:
  !*           - NUMBEROBS: it is the number of observations at time L+1
  !*           - OBSVALUE: vector of observations at time L+1
  !*           - SQRTCOVOBS: it sets the reduced rank square root of the covariance matrix of observations errors at time L+1 dur&
  !*                         &ing execution, and after that it is destroyed
  !*
  !*******************************************************************************************************************************
  subroutine rrsqrtkf_corrector( l , state , sqrtcov , sqrtcovfor )
    implicit none
    integer( kind = pin ) , intent( in ) :: l
    real( kind = pre ) , dimension( dimspacestate ) , intent( inout ) :: state
    real( kind = pre ) , dimension( dimspacestate , modesanalysis ) , intent( out ) :: sqrtcov
    real( kind = pre ) , dimension( dimspacestate , modesanalysis + modesmodel ) , intent( in ) :: sqrtcovfor
    character( kind = pch , len = 100_pin ) :: option
    real( kind = pre ) , allocatable , dimension( : ) :: vectoraux1
    real( kind = pre ) , allocatable , dimension( : , : ) :: matrixaux1
    real( kind = pre ) , allocatable , dimension( : , : ) :: matrixaux2
    real( kind = pre ) , allocatable , dimension( : , : ) :: matrixaux3
    real( kind = pre ) , allocatable , dimension( : , : ) :: matrixaux4
    real( kind = pre ) , allocatable , dimension( : , : ) :: matrixaux5
    real( kind = pre ) , allocatable , dimension( : , : ) :: matrixaux6
    allocate( sqrtcovobs( numberobs , modesobs ) )
    call observations_sqrtcovobs( l + 1_pin )
    allocate( matrixaux1( numberobs , modesanalysis + modesmodel ) )
    option = 'tangobsopN'
    call parallel_operator( l + 1_pin , sqrtcovfor , matrixaux1 , option )
    allocate( matrixaux2( numberobs , modesanalysis + modesmodel + modesobs ) )
    matrixaux2( 1_pin : numberobs , 1_pin : modesanalysis + modesmodel ) = matrixaux1( 1_pin : numberobs , 1_pin : modesanalysis &
         &+ modesmodel )
    matrixaux2( 1_pin : numberobs , modesanalysis + modesmodel + 1_pin : modesanalysis + modesmodel + modesobs ) = sqrtcovobs( 1_&
         &pin : numberobs , 1_pin : modesobs )
    allocate( matrixaux3( numberobs , numberobs ) )
    call parallel_gemm( 'N' , 'T' , numberobs , modesanalysis + modesmodel + modesobs , matrixaux2 , numberobs , modesanalysis + &
         &modesmodel + modesobs , matrixaux2 , numberobs , numberobs , matrixaux3 , 1.0_pre , 0.0_pre )
    deallocate( matrixaux2 )
    allocate( vectoraux1( numberobs ) )
    call observations_obsop( l + 1_pin , state , vectoraux1 )
    vectoraux1 = obsvalue - vectoraux1
    allocate( matrixaux4( numberobs , 1_pin + modesanalysis + modesmodel + modesobs ) )
    matrixaux4( 1_pin : numberobs , 1_pin ) = vectoraux1( 1_pin : numberobs )
    matrixaux4( 1_pin : numberobs , 1_pin + 1_pin : 1_pin + modesanalysis + modesmodel ) = matrixaux1( 1_pin : numberobs , 1_pin &
         &: modesanalysis + modesmodel )
    matrixaux4( 1_pin : numberobs , 1_pin + modesanalysis + modesmodel + 1_pin : 1_pin + modesanalysis + modesmodel + modesobs ) &
         &= sqrtcovobs( 1_pin : numberobs , 1_pin : modesobs ) 
    deallocate( sqrtcovobs )
    deallocate( vectoraux1 )
    call parallel_gesv( numberobs , matrixaux3 , 1_pin + modesanalysis + modesmodel + modesobs , matrixaux4 )
    deallocate( matrixaux3 )
    allocate( matrixaux5( modesanalysis + modesmodel , 1_pin + modesanalysis + modesmodel + modesobs ) )
    call parallel_gemm( 'T' , 'N' , numberobs , modesanalysis + modesmodel , matrixaux1 , numberobs , 1_pin + modesanalysis + mod&
         &esmodel + modesobs , matrixaux4 , modesanalysis + modesmodel , 1_pin + modesanalysis + modesmodel + modesobs , matrixau&
         &x5 , 1.0_pre , 0.0_pre )
    deallocate( matrixaux4 )
    deallocate( matrixaux1 )
    call parallel_gemm( 'N' , 'N' , dimspacestate , modesanalysis + modesmodel , sqrtcovfor , modesanalysis + modesmodel , 1_pin &
         &, matrixaux5( 1_pin : modesanalysis + modesmodel , 1_pin ) , dimspacestate , 1_pin , state , 1.0_pre , 1.0_pre )
    allocate( matrixaux6( dimspacestate , modesanalysis + modesmodel + modesobs ) )
    matrixaux6( 1_pin : dimspacestate , 1_pin : modesanalysis + modesmodel ) = sqrtcovfor( 1_pin : dimspacestate , 1_pin : modesa&
         &nalysis + modesmodel )
    call parallel_gemm( 'N' , 'N' , dimspacestate , modesanalysis + modesmodel , sqrtcovfor , modesanalysis + modesmodel , modesa&
         &nalysis + modesmodel , matrixaux5( 1_pin : modesanalysis + modesmodel , 1_pin + 1_pin : 1_pin + modesanalysis + modesmo&
         &del ) , dimspacestate , modesanalysis + modesmodel , matrixaux6( 1_pin : dimspacestate , 1_pin : modesanalysis + modesm&
         &odel ) , - 1.0_pre , 1.0_pre )
    call parallel_gemm( 'N' , 'N' , dimspacestate , modesanalysis + modesmodel , sqrtcovfor , modesanalysis + modesmodel , modeso&
         &bs , matrixaux5( 1_pin : modesanalysis + modesmodel , 1_pin + modesanalysis + modesmodel + 1_pin : 1_pin + modesanalysi&
         &s + modesmodel + modesobs ) , dimspacestate , modesobs , matrixaux6( 1_pin : dimspacestate , modesanalysis + modesmodel&
         & + 1_pin : modesanalysis + modesmodel + modesobs ) , 1.0_pre , 0.0_pre )
    deallocate( matrixaux5 )
    call parallel_svdcut( dimspacestate , modesanalysis + modesmodel + modesobs , modesanalysis , matrixaux6 , sqrtcov )
    deallocate( matrixaux6 )
  end subroutine rrsqrtkf_corrector
  !*
end module rrsqrtkf
