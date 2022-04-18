!*********************************************************************************************************************************
!*
!* MODULE: SOLUTION
!*
!* PURPOSE: sets the exact solution (useful to compute errors)
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
!* SCOPE: This module belongs to the example
!* 
!* AUTHOR: Germán Ariel Torres
!*         Facultad de Matemática, Astronomía y Física
!*         Universidad Nacional de Córdoba
!*         Argentina
!*
!* EMAIL: torres@famaf.unc.edu.ar
!*
!*********************************************************************************************************************************
module solution
  !*
  use precision
  !*
  implicit none
  !*
  real( kind = pre ) , public :: bk !* parameter for solution
  real( kind = pre ) , public :: tf !* parameter for solution
  real( kind = pre ) , public :: sf !* parameter for solution
  real( kind = pre ) , public :: df !* parameter for solution
  !*
contains
  !*
  !* setting parameters 
  subroutine solution_parametersup()
    implicit none
    bk = 3.0_pre
    tf = 1000.0_pre
    sf = 5.0_pre
    df = 0.5_pre
  end subroutine solution_parametersup
  !*
  !* exact solution
  function solution_solution( t , x , y )
    implicit none
    real( kind = pre ) , intent( in ) :: t
    real( kind = pre ) , intent( in ) :: x
    real( kind = pre ) , intent( in ) :: y
    real( kind = pre ) :: solution_solution
    solution_solution = bk + sin( tf * t ) * sin( sf * x * y ) * exp( - df * x )
  end function solution_solution
  !*
  !* derivative of the solution
  function solution_ddt( t , x, y )
    implicit none
    real( kind = pre ) , intent( in ) :: t
    real( kind = pre ) , intent( in ) :: x
    real( kind = pre ) , intent( in ) :: y
    real( kind = pre ) :: solution_ddt
    solution_ddt = cos( tf * t ) * tf * sin( sf * x * y ) * exp( - df * x )
  end function solution_ddt
  !*
  !* derivative of the solution
  function solution_d2dt2( t , x , y )
    implicit none
    real( kind = pre ) , intent( in ) :: t
    real( kind = pre ) , intent( in ) :: x
    real( kind = pre ) , intent( in ) :: y
    real( kind = pre ) :: solution_d2dt2
    solution_d2dt2 = - sin( tf * t ) * tf**2_pin * sin( sf * x * y ) * exp( - df * x )
  end function solution_d2dt2
  !*
  !* derivative of the solution
  function solution_ddx( t , x , y )
    implicit none
    real( kind = pre ) , intent( in ) :: t
    real( kind = pre ) , intent( in ) :: x
    real( kind = pre ) , intent( in ) :: y
    real( kind = pre ) :: solution_ddx
    solution_ddx = sin( tf * t ) * cos( sf * x * y ) * sf * y * exp( - df * x ) - sin( tf * t ) * sin( sf * x * y ) * df * exp( -&
         & df * x )
  end function solution_ddx
  !*
  !* derivative of the solution
  function solution_d2dx2( t , x , y )
    implicit none
    real( kind = pre ) , intent( in ) :: t
    real( kind = pre ) , intent( in ) :: x
    real( kind = pre ) , intent( in ) :: y
    real( kind = pre ) :: solution_d2dx2
    solution_d2dx2 = - sin( tf * t ) * sin( sf * x * y ) * sf**2_pin * y**2_pin * exp( - df * x ) - 2.0_pre * sin( tf * t ) * cos&
         &( sf * x * y ) * sf * y * df * exp( - df * x ) + sin( tf * t ) * sin( sf * x * y ) * df**2_pin * exp( - df * x )
  end function solution_d2dx2
  !*
  !* derivative of the solution
  function solution_d3dx3( t , x , y )
    implicit none
    real( kind = pre ) , intent( in ) :: t
    real( kind = pre ) , intent( in ) :: x
    real( kind = pre ) , intent( in ) :: y
    real( kind = pre ) :: solution_d3dx3
    solution_d3dx3 = - sin( tf * t) * cos( sf * x * y ) * sf**3_pin * y**3_pin * exp( - df * x ) + 3.0_pre * sin( tf * t ) * sin(&
         & sf * x * y ) * sf**2_pin * y**2_pin * df * exp( - df * x ) + 3.0_pre * sin( tf * t ) * cos( sf * x * y ) * sf * y * df&
         &**2_pin * exp( - df * x ) - sin( tf * t ) * sin( sf * x * y ) * df**3_pin * exp( - df * x )
  end function solution_d3dx3
  !*
  !* derivative of the solution
  function solution_d4dx4( t , x , y )
    implicit none
    real( kind = pre ) , intent( in ) :: t
    real( kind = pre ) , intent( in ) :: x
    real( kind = pre ) , intent( in ) :: y
    real( kind = pre ) :: solution_d4dx4
    solution_d4dx4 = sin( tf * t ) * sin( sf * x * y ) * sf**4_pin * y**4_pin * exp( - df * x ) + 4.0_pre * sin( tf * t ) * cos( &
         &sf * x * y ) * sf**3_pin * y**3_pin * df * exp( - df * x ) - 6.0_pre * sin( tf * t ) * sin( sf * x * y ) * sf**2_pin * &
         &y**2_pin * df**2_pin * exp( - df * x ) - 4.0_pre * sin( tf * t ) * cos( sf * x * y ) * sf * y * df**3_pin * exp( - df *&
         & x ) + sin( tf * t ) * sin( sf * x * y ) * df**4_pin * exp( - df * x )
  end function solution_d4dx4
  !*
  !* derivative of the solution
  function solution_ddy( t , x , y )
    implicit none
    real( kind = pre ) , intent( in ) :: t
    real( kind = pre ) , intent( in ) :: x
    real( kind = pre ) , intent( in ) :: y
    real( kind = pre ) :: solution_ddy
    solution_ddy = sin( tf * t ) * cos( sf * x * y ) * sf * x * exp( - df * x )
  end function solution_ddy
  !*
  !* derivative of the solution
  function solution_d2dy2( t , x , y )
    implicit none
    real( kind = pre ) , intent( in ) :: t
    real( kind = pre ) , intent( in ) :: x
    real( kind = pre ) , intent( in ) :: y
    real( kind = pre ) :: solution_d2dy2
    solution_d2dy2 = - sin( tf * t ) * sin( sf * x * y ) * sf**2_pin * x**2_pin * exp( - df * x )
  end function solution_d2dy2
  !*
  !* derivative of the solution
  function solution_d3dy3( t , x , y )
    implicit none
    real( kind = pre ) , intent( in ) :: t
    real( kind = pre ) , intent( in ) :: x
    real( kind = pre ) , intent( in ) :: y
    real( kind = pre ) :: solution_d3dy3
    solution_d3dy3 = - sin( tf * t ) * cos( sf * x * y ) * sf**3_pin * x**3_pin * exp( - df * x )
  end function solution_d3dy3
  !*
  !* derivative of the solution
  function solution_d4dy4( t , x , y )
    implicit none
    real( kind = pre ) , intent( in ) :: t
    real( kind = pre ) , intent( in ) :: x
    real( kind = pre ) , intent( in ) :: y
    real( kind = pre ) :: solution_d4dy4
    solution_d4dy4 = sin( tf * t ) * sin( sf * x * y ) * sf**4_pin * x**4_pin * exp( - df * x )
  end function solution_d4dy4
  !*
end module solution
