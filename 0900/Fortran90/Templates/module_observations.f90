!*********************************************************************************************************************************
!*
!* MODULE: OBSERVATIONS
!* 
!* PURPOSE: sets the observations and covariance matrices
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
!*          - MODESOBS (integer): number of observation modes
!*          - NUMBEROBS (integer): number of observations 
!*          - IFOBS (logical): flag to determine if there are observations
!*          - OBSVALUE (real array): observation values
!*          - COVOBS (real array): covariance matrix of observation errors
!*          - SQRTCOVOBS (real array): square root of the covariance matrix of observation errors
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
  !*
  implicit none
  !*
  integer( kind = pin ) , public :: modesobs !* number of observation modes
  integer( kind = pin ) , public :: numberobs !* number of observations 
  logical( kind = plo ) , public :: ifobs !* flag to determine if there are observations
  real( kind = pre ) , allocatable , dimension( : ) , public :: obsvalue !* observation values
  real( kind = pre ) , allocatable , dimension( : , : ) , public :: covobs !* covariance matrix of observation errors
  real( kind = pre ) , allocatable , dimension( : , : ) , public :: sqrtcovobs !* square root of the covariance matrix of observa&
                                                                               !* &tion errors
  !*
contains
  !*
  !* setting parameters
  subroutine observations_parametersup()
    implicit none
    !*****************************************************************************************************************************
    !* SET MODESOBS
    !*****************************************************************************************************************************
    !* modesobs = 
    !*****************************************************************************************************************************
  end subroutine observations_parametersup
  !*
  !* decide if there are observations
  subroutine observations_ifobservations( l )
    implicit none
    integer( kind = pin ) , intent( in ) :: l
    !*****************************************************************************************************************************
    !* SET IFOBS
    !*****************************************************************************************************************************
  end subroutine observations_ifobservations
  !*
  !* setting number of observations
  subroutine observations_numberobs( l )
    implicit none
    integer( kind = pin ) , intent( in ) :: l
    !*****************************************************************************************************************************
    !* SET NUMBEROBS
    !*****************************************************************************************************************************
  end subroutine observations_numberobs
  !*
  !* setting observation value
  subroutine observations_obsvalue( l )
    implicit none
    integer( kind = pin ) , intent( in ) :: l
    !*****************************************************************************************************************************
    !* SET OBSVALUE
    !*****************************************************************************************************************************
  end subroutine observations_obsvalue
  !*
  !* tangent of the observation operator
  subroutine observations_tangobsop( l , vectorin , vectorout )
    implicit none
    integer( kind = pin ) , intent( in ) :: l
    !*****************************************************************************************************************************
    !* SET RANK AND SIZE OF VECTORIN AND VECTOROUT
    !*****************************************************************************************************************************
    !* real( kind = pre ) , dimension( ... ) , intent( in ) :: vectorin
    !* real( kind = pre ) , dimension( ... ) , intent( out ) :: vectorout
    !*****************************************************************************************************************************
    !* SET STATEOUT
    !*****************************************************************************************************************************
  end subroutine observations_tangobsop
  !*
  !* observation operator
  subroutine observations_obsop( l , vectorin , vectorout )
    implicit none
    integer( kind = pin ) , intent( in ) :: l
    !*****************************************************************************************************************************
    !* SET RANK AND SIZE OF VECTORIN AND VECTOROUT
    !*****************************************************************************************************************************
    !* real( kind = pre ) , dimension( ... ) , intent( in ) :: vectorin
    !* real( kind = pre ) , dimension( ... ) , intent( out ) :: vectorout
    !*****************************************************************************************************************************
    !* SET VECTOROUT
    !*****************************************************************************************************************************
  end subroutine observations_obsop
  !*
  !* setting square root of the covariance matrix of observation errors
  subroutine observations_sqrtcovobs( l )
    implicit none
    integer( kind = pin ) , intent( in ) :: l
    !*****************************************************************************************************************************
    !* SET SQRTCOVOBS
    !*****************************************************************************************************************************
  end subroutine observations_sqrtcovobs
  !*
  !* setting covariance matrix of observation errors
  subroutine observations_covobs( l )
    implicit none
    integer( kind = pin ) , intent( in ) :: l
    !*****************************************************************************************************************************
    !* SET COVOBS
    !*****************************************************************************************************************************
  end subroutine observations_covobs
  !*
end module observations
