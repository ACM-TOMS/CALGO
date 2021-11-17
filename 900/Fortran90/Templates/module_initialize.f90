!*********************************************************************************************************************************
!*
!* MODULE: INITIALIZE
!*
!* PURPOSE: initialize filter
!*
!* DEPENDENCIES:
!*               - PRECISION
!* 
!* GLOBALS:
!*          - LOW (integer): number of bytes for simple precision variables (from PRECISION)
!*          - HIGH (integer): number of bytes for double precision variables (from PRECISION)
!*          - PLO (integer): precision for logical variables (from PRECISION)
!*          - PCH (integer): precision for character variables (from PRECISION)
!*          - PIN (integer): precision for integer variables (from PRECISION)
!*          - PRE (integer): precision for real variables (from PRECISION)
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
  !*
  implicit none
  !*
  integer( kind = pin ) , public :: dimspacestate !* dimension of the space state
  integer( kind = pin ) , public :: numbersamples !* number of samples
  integer( kind = pin ) , public :: modesanalysis !* number of analysis modes
  logical( kind = plo ) , public :: first !* flag to determine the first step
  !*
contains
  !*
  !* set initialization parameters
  subroutine initialize_parametersup()
    implicit none
    !*****************************************************************************************************************************
    !* EDIT HERE          
    !*****************************************************************************************************************************
    !* dimspacestate =    
    !* numbersamples =    
    !* modesanalysis =    
    !* first =            
    !*****************************************************************************************************************************
  end subroutine initialize_parametersup
  !*
  !* set initial state and initial covariance matrix
  subroutine initialize_initialize1( state , covariance )
    implicit none
    real( kind = pre ) , dimension( dimspacestate ) , intent( out ) :: state
    real( kind = pre ) , dimension( dimspacestate , dimspacestate ) , intent( out ) :: covariance
    !*****************************************************************************************************************************
    !* SET HERE STATE AND COVARIANCE 
    !*****************************************************************************************************************************
  end subroutine initialize_initialize1
  !*
  !* set initial state and initial covariance matrix
  subroutine initialize_initialize2( state , sqrtcov )
    implicit none
    real( kind = pre ) , dimension( dimspacestate ) , intent( out ) :: state
    real( kind = pre ) , dimension( dimspacestate , modesanalysis ) , intent( out ) :: sqrtcov
    !*****************************************************************************************************************************
    !* SET HERE STATE AND SQRTCOV (IN CASE YOU USE A SQUARE ROOT FILTER)
    !*****************************************************************************************************************************
  end subroutine initialize_initialize2
  !*
end module initialize
