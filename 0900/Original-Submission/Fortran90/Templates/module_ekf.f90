!*********************************************************************************************************************************
!*
!* MODULE: EKF
!*
!* PURPOSE: implementation of the Extended Kalman Filter
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
!*          - MODELERROR (real array): covariance matrix of model errors (from MODEL)
!*          - SQRTMODELERROR (real array): square root of covariance matrix of model errors (from MODEL)
!*
!*          - MODESOBS (integer): number of observation modes (from OBSERVATIONS)
!*          - NUMBEROBS (integer): number of observations (from OBSERVATIONS)
!*          - IFOBS (logical): flag to determine if there are observations (from OBSERVATIONS)
!*          - OBSVALUE (real array): observation values (from OBSERVATIONS)
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
!*
!* AUTHOR: Germán Ariel Torres
!*         Facultad de Matemática, Astronomía y Física
!*         Universidad Nacional de Córdoba
!*         Argentina
!*
!* EMAIL: torres@famaf.unc.edu.ar
!*
!*********************************************************************************************************************************
module ekf
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
  !* FUNCTION/SUBROUTINE: EKF_PREDICTOR
  !*
  !* PURPOSE: makes the prediction step of the Extended Kalman Filter, that is:
  !*
  !*            forecast(l+1) = model_{l -> l+1} ( analysis(l) )
  !*            covforecast(l+1) = tangmodel_{l -> l+1} covanalysis(l) tangmodel_{l -> l+1}^T + coverrormodel(l)
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
  !*
  !* OUTPUTS: none
  !*
  !* CALLS:
  !*        - MODEL_MODEL (from MODEL)
  !*        - PARALLEL_OPERATOR (from PARALLEL)
  !*        - MODEL_MODELERROR (from MODEL)
  !* 
  !* COMMENTS: none
  !*
  !*******************************************************************************************************************************
  subroutine ekf_predictor( l , state , covariance )
    implicit none
    integer( kind = pin ) , intent( in ) :: l
    real( kind = pre ) , dimension( dimspacestate ) , intent( inout ) :: state
    real( kind = pre ) , dimension( dimspacestate , dimspacestate ) , intent( inout ) :: covariance
    character( kind = pch , len = 1000_pin ) :: option
    real( kind = pre ) , allocatable , dimension( : ) :: stateaux
    real( kind = pre ) , allocatable , dimension( : , : ) :: matrixaux1
    real( kind = pre ) , allocatable , dimension( : , : ) :: matrixaux3
    allocate( stateaux( dimspacestate ) )
    call model_model( l , state , stateaux )
    state = stateaux
    deallocate( stateaux )
    allocate( matrixaux1( dimspacestate , dimspacestate ) )
    option = 'tangmodelT'
    call parallel_operator( l , covariance , matrixaux1 , option )
    allocate( matrixaux3( dimspacestate , dimspacestate ) )
    option = 'tangmodelT'
    call parallel_operator( l , matrixaux1 , matrixaux3 , option )
    deallocate( matrixaux1 )
    allocate( modelerror( dimspacestate , dimspacestate ) )
    call model_modelerror( l )
    covariance = matrixaux3 + modelerror
    deallocate( matrixaux3 )
    deallocate( modelerror )
  end subroutine ekf_predictor
  !*******************************************************************************************************************************
  !* 
  !* FUNCTION/SUBROUTINE: EKF_CORRECTOR
  !*
  !* PURPOSE: makes the correction step of the Extended Kalman Filter, that is:
  !*            
  !*            gainmatrix(l+1) = covforecast(l+1) * tangobsop(l+1) * [ tangobsop(l+1) * covforecast(l+1) * tangobsop(l+1)^T + co&
  !*                              &vobs(l+1) ]^{-1}
  !*            analysis(l+1) = forecast(l+1) + gainmatrix(l+1) * [ obsvalue(l+1) - obsop(l+1) ( forecast(l+1) ) ]
  !*            covanalysis(l+1) = [ I - gainmatrix(l+1) * tangobsop(l+1) ] * covforecast(l+1) * [ I - gainmatrix(l+1) * tangobso&
  !*                               &p(l+1) ]^T + gainmatrix(l+1) * covobs(l+1) * gainmatrix(l+1)^T
  !*
  !* INPUTS:
  !*         - L (integer): current time step
  !*
  !* INPUTS/OUTPUTS:
  !*                 - STATE (real array of size DIMSPACESTATE): on input the forecast at time L+1, on exit the analysis at time L&
  !*                                                             &+1
  !*                 - COVARIANCE (real array of size DIMSPACESTATE x DIMSPACESTATE): on input the covariance matrix of forecast &
  !*                                                                                  &errors at time L+1, on output the covarian&
  !*                                                                                  &ce matrix of analysis errors at time L+1
  !*
  !* OUTPUTS: none
  !*
  !* CALLS:
  !*        - PARALLEL_OPERATOR (from PARALLEL)        
  !*        - PARALLEL_GETRF (from PARALLEL)
  !*        - PARALLEL_GEMM (from PARALLEL)
  !*        - OBSERVATIONS_COVOBS (from OBSERVATIONS)
  !*        - OBSERVATIONS_OBSOP (from OBSERVATIONS)
  !*        - DGETRS (from LAPACK)
  !*        - SGETRS (from LAPACK)
  !* 
  !* COMMENTS:
  !*           - NUMBEROBS: it is the number of observations at time L+1
  !*           - OBSVALUE: vector of observations at time L+1
  !*           - COVOBS: it sets the covariance matrix of observation errors at time L+1 during execution, and after that it is d&
  !*                     &estroyed
  !*
  !*******************************************************************************************************************************
  subroutine ekf_corrector( l , state , covariance )
    implicit none
    integer( kind = pin ) , intent( in ) :: l
    real( kind = pre ) , dimension( dimspacestate ) , intent( inout ) :: state
    real( kind = pre ) , dimension( dimspacestate , dimspacestate ) , intent( inout ) :: covariance
    character( kind = pch , len = 1000_pin ) :: option
    integer( kind = pin ) :: info
    integer( kind = pin ) , allocatable , dimension( : ) :: ipiv
    real( kind = pre ) , allocatable , dimension( : ) :: vectoraux1
    real( kind = pre ) , allocatable , dimension( : , : ) :: matrixaux2
    real( kind = pre ) , allocatable , dimension( : , : ) :: matrixaux4
    real( kind = pre ) , allocatable , dimension( : , : ) :: matrixaux5
    real( kind = pre ) , allocatable , dimension( : , : ) :: matrixaux6
    real( kind = pre ) , allocatable , dimension( : , : ) :: matrixaux7
    real( kind = pre ) , allocatable , dimension( : , : ) :: matrixaux8
    allocate( matrixaux2( numberobs , dimspacestate ) )
    option = 'tangobsopT'
    call parallel_operator( l + 1_pin , covariance , matrixaux2 , option )
    allocate( matrixaux4( numberobs , numberobs ) )
    option = 'tangobsopT'
    call parallel_operator( l + 1_pin , matrixaux2 , matrixaux4 , option )
    allocate( covobs( numberobs , numberobs ) )
    call observations_covobs( l + 1_pin )
    matrixaux4 = matrixaux4 + covobs
    allocate( vectoraux1( numberobs ) )
    call observations_obsop( l + 1_pin , state , vectoraux1 )
    vectoraux1 = obsvalue - vectoraux1
    allocate( matrixaux5( numberobs , 1_pin + dimspacestate + numberobs ) )
    matrixaux5( 1_pin : numberobs , 1_pin ) = vectoraux1
    matrixaux5( 1_pin : numberobs , 2_pin : 1_pin + dimspacestate ) = matrixaux2
    matrixaux5( 1_pin : numberobs , 2_pin + dimspacestate : 1_pin + dimspacestate + numberobs ) = covobs
    deallocate( vectoraux1 )
    deallocate( covobs )
    allocate( ipiv( numberobs ) )
    call parallel_getrf( numberobs , numberobs , matrixaux4 , ipiv )
    !* trans = 'N'
    !* n     = numberobs 
    !* nrhs  = 1_pin + dimspacestate + numberobs
    !* a     = matrixaux4
    !* lda   = numberobs
    !* ipiv  = ipiv
    !* b     = matrixaux5
    !* ldb   = numberobs
    !* info  = info
    if ( pre == high ) then
       call dgetrs( 'N' , numberobs , 1_pin + dimspacestate + numberobs , matrixaux4 , numberobs , ipiv , matrixaux5 , numberobs &
            &, info )
    else if ( pre == low ) then
       call sgetrs( 'N' , numberobs , 1_pin + dimspacestate + numberobs , matrixaux4 , numberobs , ipiv , matrixaux5 , numberobs &
            &, info )
    end if
    if ( info /= 0_pin ) then
       print*,'error: info /= 0'
       print*,'info = ' , info
       stop
    end if
    allocate( matrixaux6( 1_pin + dimspacestate + numberobs , dimspacestate ) )
    call parallel_gemm( 'T' , 'N' , numberobs , 1_pin + dimspacestate + numberobs , matrixaux5 , numberobs , dimspacestate , matr&
         &ixaux2 , 1_pin + dimspacestate + numberobs , dimspacestate , matrixaux6 , 1.0_pre , 0.0_pre )
    deallocate( matrixaux5 )
    state = state + matrixaux6( 1_pin , 1_pin : dimspacestate )
    covariance = covariance - matrixaux6( 2_pin : dimspacestate + 1_pin , 1_pin : dimspacestate )
    allocate( matrixaux7( numberobs , dimspacestate + dimspacestate ) )
    option = 'tangobsopN'
    call parallel_operator( l + 1_pin , covariance , matrixaux7( 1_pin : numberobs , 1_pin : dimspacestate ) , option )
    matrixaux7( 1_pin : numberobs , dimspacestate + 1_pin : dimspacestate + dimspacestate ) = matrixaux6( 2_pin + dimspacestate :&
         & 1_pin + dimspacestate + numberobs , 1_pin : dimspacestate ) 
    deallocate( matrixaux6 )
    if ( pre == high ) then
       call dgetrs( 'N' , numberobs , dimspacestate + dimspacestate , matrixaux4 , numberobs , ipiv , matrixaux7 , numberobs , in&
            &fo )
    else if ( pre == low ) then
       call sgetrs( 'N' , numberobs , dimspacestate + dimspacestate , matrixaux4 , numberobs , ipiv , matrixaux7 , numberobs , in&
            &fo )
    end if
    if ( info /= 0_pin ) then
       print*,'error: info /= 0'
       print*,'info = ' , info
       stop
    end if
    deallocate( ipiv )
    deallocate( matrixaux4 )
    allocate( matrixaux8( dimspacestate , dimspacestate + dimspacestate ) )
    call parallel_gemm( 'T' , 'N' , numberobs , dimspacestate , matrixaux2 , numberobs , dimspacestate + dimspacestate , matrixau&
         &x7 , dimspacestate , dimspacestate + dimspacestate , matrixaux8 , 1.0_pre , 0.0_pre )
    deallocate( matrixaux2 )
    deallocate( matrixaux7 )
    covariance = covariance - matrixaux8( 1_pin : dimspacestate , 1_pin : dimspacestate ) + matrixaux8( 1_pin : dimspacestate , d&
         &imspacestate + 1_pin : dimspacestate + dimspacestate )
    deallocate( matrixaux8 )
  end subroutine ekf_corrector
  !*
end module ekf
