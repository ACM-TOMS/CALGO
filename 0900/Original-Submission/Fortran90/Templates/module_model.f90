!*********************************************************************************************************************************
!*
!* MODULE: MODEL
!*
!* PURPOSE: it sets the propagation model
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
!*          - MODESMODEL (integer): number of model modes
!*          - MODELERROR (real array): covariance matrix of model errors
!*          - SQRTMODELERROR (real array): square root of the covariance matrix of model errors
!*
!* AUTHOR: Germán Ariel Torres
!*         Facultad de Matemática, Astronomía y Física
!*         Universidad Nacional de Córdoba
!*         Argentina
!*
!* EMAIL: torres@famaf.unc.edu.ar
!*
!*********************************************************************************************************************************
module model
  !*
  use precision
  !*
  implicit none
  !*
  integer( kind = pin ) , public :: modesmodel !* number of model modes
  real( kind = pre ) , allocatable , dimension( : , : ) , public :: modelerror !* covariance matrix of model errors
  real( kind = pre ) , allocatable , dimension( : , : ) , public :: sqrtmodelerror !* square root of covariance matrix of model e&
                                                                                   !* &rrors
  !*
contains
  !*
  !* model parameters
  subroutine model_parametersup()
    implicit none
    !*****************************************************************************************************************************
    !* EDIT HERE
    !*****************************************************************************************************************************
    !* modesmodel = 
    !*****************************************************************************************************************************
  end subroutine model_parametersup
  !*
  !* model
  subroutine model_model( l , statein , stateout )
    implicit none
    integer( kind = pin ) , intent( in ) :: l
    !*****************************************************************************************************************************
    !* SET RANK AND SIZE OF STATEIN AND STATEOUT
    !*****************************************************************************************************************************
    !* real( kind = pre ) , dimension( ... ) , intent( in ) :: statein
    !* real( kind = pre ) , dimension( ... ) , intent( out ) :: stateout
    !*****************************************************************************************************************************
    !* SET STATEOUT
    !*****************************************************************************************************************************
  end subroutine model_model
  !*
  !* tangent model
  subroutine model_tangmodel( l , statein , stateout )
    implicit none
    integer( kind = pin ) , intent( in ) :: l
    !*****************************************************************************************************************************
    !* SET RANK AND SIZE OF STATEIN AND STATEOUT
    !*****************************************************************************************************************************
    !* real( kind = pre ) , dimension( ... ) , intent( in ) :: statein
    !* real( kind = pre ) , dimension( ... ) , intent( out ) :: stateout
    !*****************************************************************************************************************************
    !* SET STATEOUT
    !*****************************************************************************************************************************
  end subroutine model_tangmodel
  !*
  !* square root of the covariance matrix of model errors 
  subroutine model_sqrtmodelerror( l )
    implicit none
    integer( kind = pin ) , intent( in ) :: l
    !*****************************************************************************************************************************
    !* SET SQRTMODELERROR
    !*****************************************************************************************************************************
  end subroutine model_sqrtmodelerror
  !*
  !* covariance matrix of model errors
  subroutine model_modelerror( l )
    implicit none
    integer( kind = pin ) , intent( in ) :: l
    !*****************************************************************************************************************************
    !* SET MODELERROR
    !*****************************************************************************************************************************
  end subroutine model_modelerror
  !*
end module model
