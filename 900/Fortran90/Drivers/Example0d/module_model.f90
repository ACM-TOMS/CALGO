!*********************************************************************************************************************************
!*
!* MODULE: MODEL
!*
!* PURPOSE: sets the propagation model
!*
!* DEPENDENCIES:
!*               - PRECISION
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
!*          - NT (integer): number of time steps (from DOMAIN)
!*          - A (real): initial time (from DOMAIN)
!*          - DELTAT (real): time step (from DOMAIN)
!*          - T (real array): vector of time (from DOMAIN)
!*
!*          - MODESMODEL (integer): number of model modes
!*          - MODELERROR (real array): covariance matrix of model errors
!*          - SQRTMODELERROR (real array): square root of the covariance matrix of model errors
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
!*
module model
  !*
  use precision
  use solution
  use domain
  !*
  implicit none
  !*
  integer( kind = pin ) , public :: modesmodel !* number of model modes
  real( kind = pre ) , allocatable , dimension( : , : ) , public :: modelerror !* model errors covariance matrix
  real( kind = pre ) , allocatable , dimension( : , : ) , public :: sqrtmodelerror !* square root of the model errors covariance &
                                                                                   !* &matrix
  !*
contains
  !*
  !* set model parameters
  subroutine model_parametersup()
    implicit none
    modesmodel = 1_pin
  end subroutine model_parametersup
  !*
  !* function f of the ordinary differential equation
  function model_f( t , x )
    implicit none
    real( kind = pre ) , intent( in ) :: t 
    real( kind = pre ) , intent( in ) :: x 
    real( kind = pre ) :: model_f
    model_f = cos( t )
  end function model_f
  !*
  !* time derivative of f
  function model_dfdt( t , x )
    implicit none
    real( kind = pre ) , intent( in ) :: t
    real( kind = pre ) , intent( in ) :: x
    real( kind = pre ) :: model_dfdt
    model_dfdt = - sin( t )
  end function model_dfdt
  !*
  !* x-derivative of f
  function model_dfdx( t , x )
    implicit none
    real( kind = pre ) , intent( in ) :: t
    real( kind = pre ) , intent( in ) :: x
    real( kind = pre ) :: model_dfdx
    model_dfdx = 0.0_pre
  end function model_dfdx
  !*
  !* model ( euler's method )
  subroutine model_model( l , statein , stateout )
    implicit none
    integer( kind = pin ) , intent( in ) :: l
    real( kind = pre ) , dimension( 1_pin ) , intent( in ) :: statein
    real( kind = pre ) , dimension( 1_pin ) , intent( out ) :: stateout
    stateout( 1_pin ) = statein( 1_pin ) + deltat * model_f( t( l ) , statein( 1_pin ) )
  end subroutine model_model
  !*
  !* tangent linear model
  subroutine model_tangmodel( l , statein , stateout )
    implicit none
    integer( kind = pin ) , intent( in ) :: l
    real( kind = pre ) , dimension( 1_pin ) , intent( in ) :: statein
    real( kind = pre ) , dimension( 1_pin ) , intent( out ) :: stateout
    stateout( 1_pin ) = 1.0_pre + deltat * model_dfdx( t( l ) , statein( 1_pin ) )
  end subroutine model_tangmodel
  !*
  !* set model errors covariance matrix
  subroutine model_modelerror( l )
    implicit none
    integer( kind = pin ) , intent( in ) :: l
    real( kind = pre ) :: sol
    real( kind = pre ) :: dfdt
    real( kind = pre ) :: dfdx
    real( kind = pre ) :: f
    real( kind = pre ) :: q
    sol = solution_solution( t( l ) )
    dfdt = model_dfdt( t( l ) , sol )
    dfdx = model_dfdx( t( l ) , sol )
    f = model_f( t( l ) , sol )
    q = - deltat**2_pin / 2.0_pre * ( dfdt + dfdx * f )
    modelerror( 1_pin , 1_pin ) = q**2_pin
  end subroutine model_modelerror
  !*
  !* set the square root of the model errors covariance matrix
  subroutine model_sqrtmodelerror( l )
    implicit none
    integer( kind = pin ) , intent( in ) :: l
    real( kind = pre ) :: sol
    real( kind = pre ) :: dfdt
    real( kind = pre ) :: dfdx
    real( kind = pre ) :: f
    real( kind = pre ) :: q
    sol = solution_solution( t( l ) )
    dfdt = model_dfdt( t( l ) , sol )
    dfdx = model_dfdx( t( l ) , sol )
    f = model_f( t( l ) , sol )
    q = - deltat**2_pin / 2.0_pre * ( dfdt + dfdx * f )
    sqrtmodelerror( 1_pin , 1_pin ) = abs( q )
  end subroutine model_sqrtmodelerror
  !*
end module model
