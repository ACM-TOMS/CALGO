!*********************************************************************************************************************************
!*
!* MODULE: MODEL
!*
!* PURPOSE: sets the propagation model
!*
!* DEPENDENCIES:
!*               - PRECISION
!*               - TOOLS
!*               - RANDOM
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
!*          - IDUM (integer): seed for random generators, it must be set to a negative value (from RANDOM)
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
!*          - MODESMODEL (integer): number of model modes
!*          - ERRORBC (real): percentage of boundary condition errors
!*          - BC (real array): boundary condition
!*          - SIGMABC (real array): standard deviation of boundary condition
!*          - SQRTMODELERROR (real array) x MODESMODEL ): square root of the covariance matrix of model errors
!*          - MODELERROR (real array): covariance matrix of model errors
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
module model
  !*
  use precision
  use tools
  use random
  use solution
  use domain
  !*
  implicit none
  !*
  integer( kind = pin ) , public :: modesmodel !* number of model modes
  real( kind = pre ) , public :: errorbc !* percentage of boundary condition errors
  real( kind = pre ) , allocatable , dimension( : , : ) , public :: bc !* boundary condition
  real( kind = pre ) , allocatable , dimension( : , : ) , public :: sigmabc !* standard deviation of boundary condition
  real( kind = pre ) , allocatable , dimension( : , : ) , public :: sqrtmodelerror !* square root of covariance matrix of model e&
                                                                                   !* &rrors
  real( kind = pre ) , allocatable , dimension( : , : ) , public :: modelerror !* covariance matrix of model errors
  !*
contains
  !*
  !* model parameters
  subroutine model_parametersup()
    implicit none
    modesmodel = 200_pin
    errorbc = 0.50_pre
    allocate( bc( nx , ny ) )
    allocate( sigmabc( nx , ny ) )
  end subroutine model_parametersup
  !*
  !* destructing
  subroutine model_destructor()
    implicit none
    if ( allocated( bc ) ) deallocate( bc )
    if ( allocated( sigmabc ) ) deallocate( sigmabc )
  end subroutine model_destructor
  !*
  !* setting boundary conditions
  subroutine model_bc( l )
    implicit none
    integer( kind = pin ) , intent( in ) :: l
    integer( kind = pin ) :: i
    integer( kind = pin ) :: j
    real( kind = pre ) :: aux
    do i = 1_pin , nx
       do j = 1_pin , ny
          aux = solution_solution( t( l ) , x( i ) , y( j ) )
          sigmabc( i , j ) = errorbc * abs( aux )
          bc( i , j ) = random_normal0d( aux , sigmabc( i , j ) )
       end do
    end do
  end subroutine model_bc
  !*
  !* setting part of the model
  function model_u( t , x , y )
    implicit none
    real( kind = pre ) , intent( in ) :: t
    real( kind = pre ) , intent( in ) :: x
    real( kind = pre ) , intent( in ) :: y
    real( kind = pre ) :: model_u
    model_u = 1.0_pre
  end function model_u
  !*
  !* setting part of the model
  function model_v( t , x , y )
    implicit none
    real( kind = pre ) , intent( in ) :: t
    real( kind = pre ) , intent( in ) :: x
    real( kind = pre ) , intent( in ) :: y
    real( kind = pre ) :: model_v
    model_v = 1.0_pre
  end function model_v
  !*
  !* setting part of the model
  function model_alfa( t , x , y )
    implicit none
    real( kind = pre ) , intent( in ) :: t
    real( kind = pre ) , intent( in ) :: x
    real( kind = pre ) , intent( in ) :: y
    real( kind = pre ) :: model_alfa
    model_alfa = 1.0_pre
  end function model_alfa
  !*
  !* setting part of the model
  function model_beta( t , x , y )
    implicit none
    real( kind = pre ) , intent( in ) :: t
    real( kind = pre ) , intent( in ) :: x
    real( kind = pre ) , intent( in ) :: y
    real( kind = pre ) :: model_beta
    model_beta = 1.0_pre
  end function model_beta
  !*
  !* setting part of the model
  function model_source( t , x , y )
    implicit none
    real( kind = pre ) , intent( in ) :: t
    real( kind = pre ) , intent( in ) :: x
    real( kind = pre ) , intent( in ) :: y
    real( kind = pre ) :: model_source
    real( kind = pre ) :: ddt
    real( kind = pre ) :: ddx
    real( kind = pre ) :: ddy
    real( kind = pre ) :: d2dx2
    real( kind = pre ) :: d2dy2
    real( kind = pre ) :: u
    real( kind = pre ) :: v
    real( kind = pre ) :: alfa
    real( kind = pre ) :: beta
    ddt = solution_ddt( t , x , y ) 
    ddx = solution_ddx( t , x , y )
    ddy = solution_ddy( t , x , y )
    d2dx2 = solution_d2dx2( t , x , y )
    d2dy2 = solution_d2dy2( t , x , y )
    u = model_u( t , x , y )
    v = model_v( t , x , y )
    alfa = model_alfa( t , x , y )
    beta = model_beta( t , x , y )
    model_source = ddt + u * ddx + v * ddy - alfa * d2dx2 - beta * d2dy2
  end function model_source
  !*
  !* model
  subroutine model_model( l , statein , stateout )
    implicit none
    integer( kind = pin ) , intent( in ) :: l
    real( kind = pre ) , dimension( nx , ny ) , intent( in ) :: statein
    real( kind = pre ) , dimension( nx , ny ) , intent( out ) :: stateout
    integer( kind = pin ) :: i
    integer( kind = pin ) :: j
    real( kind = pre ) :: coefu
    real( kind = pre ) :: coefv
    real( kind = pre ) :: coefalfa
    real( kind = pre ) :: coefbeta
    real( kind = pre ) :: sourcedt
    stateout = bc
    do i = 2_pin , nx - 1_pin
       do j = 2_pin , ny - 1_pin
          coefu = model_u( t( l ) , x( i ) , y( j ) ) * deltat / ( 2.0_pre * deltax )
          coefv = model_v( t( l ) , x( i ) , y( j ) ) * deltat / ( 2.0_pre * deltay )
          coefalfa = model_alfa( t( l ) , x( i ) , y( j ) ) * deltat / deltax**2_pin
          coefbeta = model_beta( t( l ) , x( i ) , y( j ) ) * deltat / deltay**2_pin
          sourcedt = model_source( t( l ) , x( i ) , y( j ) ) * deltat
          stateout( i , j ) = statein( i , j ) - coefu * ( statein( i + 1_pin , j ) - statein( i - 1_pin , j ) ) - coefv * ( stat&
               &ein( i , j + 1_pin ) - statein( i , j - 1_pin ) ) + coefalfa * ( statein( i + 1_pin , j ) - 2.0_pre * statein( i &
               &, j ) + statein( i - 1_pin , j ) ) + coefbeta * ( statein( i , j + 1_pin ) - 2.0_pre * statein( i , j ) + statein&
               &( i , j - 1_pin ) ) + sourcedt
       end do
    end do
  end subroutine model_model
  !*
  !* tangent model
  subroutine model_tangmodel( l , statein , stateout )
    implicit none
    integer( kind = pin ) , intent( in ) :: l
    real( kind = pre ) , dimension( nx , ny ) , intent( in ) :: statein
    real( kind = pre ) , dimension( nx , ny ) , intent( out ) :: stateout
    integer( kind = pin ) :: i
    integer( kind = pin ) :: j
    real( kind = pre ) :: coefu
    real( kind = pre ) :: coefv
    real( kind = pre ) :: coefalfa
    real( kind = pre ) :: coefbeta
    stateout = 0.0_pre
    do i = 2_pin , nx - 1_pin
       do j = 2_pin , ny - 1_pin
          coefu = model_u( t( l ) , x( i ) , y( j ) ) * deltat / ( 2.0_pre * deltax )
          coefv = model_v( t( l ) , x( i ) , y( j ) ) * deltat / ( 2.0_pre * deltay )
          coefalfa = model_alfa( t( l ) , x( i ) , y( j ) ) * deltat / deltax**2_pin
          coefbeta = model_beta( t( l ) , x( i ) , y( j ) ) * deltat / deltay**2_pin
          stateout( i , j ) = statein( i , j ) - coefu * ( statein( i + 1_pin , j ) - statein( i - 1_pin , j ) ) - coefv * ( stat&
               &ein( i , j + 1_pin ) - statein( i , j - 1_pin ) ) + coefalfa * ( statein( i + 1_pin , j ) - 2.0_pre * statein( i &
               &, j ) + statein( i - 1_pin , j ) ) + coefbeta * ( statein( i , j + 1_pin ) - 2.0_pre * statein( i , j ) + statein&
               &( i , j - 1_pin ) )
       end do
    end do
  end subroutine model_tangmodel
  !*
  !* square root of the covariance matrix of model errors 
  subroutine model_sqrtmodelerror( l )
    implicit none
    integer( kind = pin ) , intent( in ) :: l
    integer( kind = pin ) :: i
    integer( kind = pin ) :: j
    integer( kind = pin ) :: m
    integer( kind = pin ) :: ij
    integer( kind = pin ) , dimension( 2_pin ) :: ijmax
    real( kind = pre ) :: d2dt2
    real( kind = pre ) :: d3dx3
    real( kind = pre ) :: d4dx4
    real( kind = pre ) :: d3dy3
    real( kind = pre ) :: d4dy4
    real( kind = pre ) :: coef1
    real( kind = pre ) :: coef2
    real( kind = pre ) :: coef3
    real( kind = pre ) :: coef4
    real( kind = pre ) :: coef5
    real( kind = pre ) :: aux
    real( kind = pre ) :: smax
    real( kind = pre ) , allocatable , dimension( : , : ) :: sigmamodelerror
    sqrtmodelerror = 0.0_pre
    allocate( sigmamodelerror( nx , ny ) )
    do j = 1_pin , ny
       do i = 1_pin , nx
          if ( i == 1_pin .or. i == nx .or. j == 1_pin .or. j == ny ) then
             sigmamodelerror( i , j ) = sigmabc( i , j )
          else
             d2dt2 = solution_d2dt2( t( l ) , x( i ) , y( j ) )
             d3dx3 = solution_d3dx3( t( l ) , x( i ) , y( j ) )
             d4dx4 = solution_d4dx4( t( l ) , x( i ) , y( j ) )
             d3dy3 = solution_d3dy3( t( l ) , x( i ) , y( j ) )
             d4dy4 = solution_d4dy4( t( l ) , x( i ) , y( j ) )
             coef1 = deltat**2_pin / 2.0_pre
             coef2 = model_u( t( l ) , x( i ) , y( j ) ) * deltat * deltax**2_pin / 6.0_pre
             coef3 = model_v( t( l ) , x( i ) , y( j ) ) * deltat * deltay**2_pin / 6.0_pre
             coef4 = model_alfa( t( l ) , x( i ) , y( j ) ) * deltat * deltax**2_pin / 12.0_pre
             coef5 = model_beta( t( l ) , x( i ) , y( j ) ) * deltat * deltay**2_pin / 12.0_pre
             aux = coef1 * d2dt2 + coef2 * d3dx3 + coef3 * d3dy3 - coef4 * d4dx4 - coef5 * d4dy4
             sigmamodelerror( i , j ) = abs( aux )
          end if
       end do
    end do
    do m = 1_pin , modesmodel
       ijmax = maxloc( sigmamodelerror( 1_pin : nx , 1_pin : ny ) )
       smax = maxval( sigmamodelerror( 1_pin : nx , 1_pin : ny ) )
       call tools_ij2s( nx , ijmax( 1_pin ) , ijmax( 2_pin ) , ij )
       sqrtmodelerror( ij , m ) = smax
       sigmamodelerror( ijmax( 1_pin ) , ijmax( 2_pin ) ) = 0.0_pre
    end do
    deallocate( sigmamodelerror )
  end subroutine model_sqrtmodelerror
  !*
  !* covariance matrix of model errors
  subroutine model_modelerror( l )
    implicit none
    integer( kind = pin ) , intent( in ) :: l
    integer( kind = pin ) :: i
    integer( kind = pin ) :: j
    integer( kind = pin ) :: ij
    real( kind = pre ) :: d2dt2
    real( kind = pre ) :: d3dx3
    real( kind = pre ) :: d4dx4
    real( kind = pre ) :: d3dy3
    real( kind = pre ) :: d4dy4
    real( kind = pre ) :: coef1
    real( kind = pre ) :: coef2
    real( kind = pre ) :: coef3
    real( kind = pre ) :: coef4
    real( kind = pre ) :: coef5
    real( kind = pre ) :: aux
    modelerror = 0.0_pre
    ij = 1_pin
    do j = 1_pin , ny
       do i = 1_pin , nx
          if ( i == 1_pin .or. i == nx .or. j == 1_pin .or. j == ny ) then
             modelerror( ij , ij ) = sigmabc( i , j )**2_pin
          else
             d2dt2 = solution_d2dt2( t( l ) , x( i ) , y( j ) )
             d3dx3 = solution_d3dx3( t( l ) , x( i ) , y( j ) )
             d4dx4 = solution_d4dx4( t( l ) , x( i ) , y( j ) )
             d3dy3 = solution_d3dy3( t( l ) , x( i ) , y( j ) )
             d4dy4 = solution_d4dy4( t( l ) , x( i ) , y( j ) )
             coef1 = deltat**2_pin / 2.0_pre
             coef2 = model_u( t( l ) , x( i ) , y( j ) ) * deltat * deltax**2_pin / 6.0_pre
             coef3 = model_v( t( l ) , x( i ) , y( j ) ) * deltat * deltay**2_pin / 6.0_pre
             coef4 = model_alfa( t( l ) , x( i ) , y( j ) ) * deltat * deltax**2_pin / 12.0_pre
             coef5 = model_beta( t( l ) , x( i ) , y( j ) ) * deltat * deltay**2_pin / 12.0_pre
             aux = coef1 * d2dt2 + coef2 * d3dx3 + coef3 * d3dy3 - coef4 * d4dx4 - coef5 * d4dy4
             modelerror( ij , ij ) = aux**2_pin
          end if
          ij = ij + 1_pin
       end do
    end do
  end subroutine model_modelerror
  !*
end module model
