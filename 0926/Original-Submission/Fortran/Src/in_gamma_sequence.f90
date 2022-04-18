module in_gamma_sequence

implicit none
save

!-----

private

public :: scaled_in_gamma

!-----

!64 bit integer and real types. In Fortran 2008,
!these could be replaced with int64 and real64.
integer   , parameter :: li  = selected_int_kind( 18 )
integer   , parameter :: dp  = selected_real_kind( 15 ) 

integer   , parameter :: dg  = digits( 1.0_dp )

real (dp) , parameter :: eps = epsilon( 1.0_dp )

!-----

integer :: j

real (dp) , dimension (2*dg+2) , parameter :: recips &
                  & = [ ( 1.0_dp / real( j , kind = dp  ) , j = 1 , 2 * dg + 2 ) ]

!-----

!Don't like to have save variables inside procedures, so we will put them here.

!After reading precomputed data, mu_min will be set to lbound( igz , 1 ) which will
!always be negative. Setting an initial value of zero allows us to simultaneously
!check whether data has been read, and, if so, whether that data is sufficient.
integer :: mu_min = 0

!Store for precomputed data. igz(j,1) is the machine number closest to x_j.
!                            igz(j,2) is the machine number closest to S_j( igz(j,1) ).
real (dp) , allocatable , dimension (:,:) :: igz

!-----

contains

!-----

  subroutine error( proc_name , msg )

    implicit none

    character (*) , intent (in) :: proc_name , msg

    !-----

    write (*,fmt='( a , " : error --- " , a )') proc_name , msg
    stop

    !-----

  end subroutine

!-----

  subroutine scaled_in_gamma( n0 , n1 , x , res )
  !Computes S( n , x ) for n = n0 , ... , n1 and x >= 0.

    implicit none

    integer   , intent (in) :: n0 , n1

    real (dp) , intent (in) :: x

    real (dp) , intent (out) , allocatable , dimension (:) :: res

    !-----

    integer       :: n
    integer       :: lambda , mu
    integer       :: up_rec_stop

    real     (dp) :: x_recip

    character (6) :: x_max_char

    !This temporary array allows us to avoid a number of conditional statements.
    real (dp) , allocatable , dimension (:) :: res_tmp

    !-----

    !Check for bad input values
    if ( x < 0.0_dp .or. n1 < n0 ) call error( "scaled_in_gamma" , "invalid input values" )

    !Recurrence in the direction of decreasing n (>0) must cease at lambda + 1.
    lambda     = ceiling( x - 1.5_dp )

#   ifdef DIAGNOSTIC
      write (*,*) " lambda = " , lambda
#   endif

    !mu is the largest integer such that x_mu >= x. We start with
    !a lower bound. In most cases this is the actual value of mu,
    !but if x is close to x_mu, it may evaluate to mu - 1.
     mu       = -ceiling( x + 0.18_dp ) 

    !Allocate a temporary storage array that always includes indices n0,
    !n1, mu and mu + 1 (for the actual value of mu, not the lower bound).
    allocate( res_tmp(min( n0 , mu ):max( n1 , mu + 2 )) )

    !-----

    if ( n1 > lambda ) then 

      res_tmp(n1) = s_series2( n1 , x )

      !Downwards recurrence
      do n = n1 - 1 , max( lambda + 1 , n0 ) , -1
        res_tmp(n) = ( 1.0_dp - x * res_tmp(n+1) ) / ( real( n , kind = dp ) + 0.5_dp )
      end do

      !Stopping point for upwards recurrence
      up_rec_stop = lambda
    else
      up_rec_stop = n1
    end if

    !-----

    need_s_mu : if ( n0 <= lambda ) then

      !Using the lower bound, we can check the value of mu_min without looking at the stored data.
      if ( mu < mu_min ) then

        if ( mu_min == 0 ) call read_precomp_data()

        if ( mu < mu_min ) then
          write (x_max_char,fmt='( f6.2 )') real( abs( mu_min + 0.18 ) )
          call error( "scaled_in_gamma" , "x is too large! " // new_line( "." )       // &
                    & "Precomputed data is available for x <= " // trim( x_max_char ) // new_line( "." ) &
                    & // "More data can be computed using in_gamma_precomp.mpl." )
        end if

      end if

      !Correct the approximate value of mu, if necessary. Note that igz(0,1) 
      !(i.e. x0), which is not defined in the paper, is set to -1 by 
      !read_precomp_data(), so it is always safe to perform this test.
      if ( igz(mu+1,1) >= x ) mu = mu + 1

#     ifdef DIAGNOSTIC
        write (*,*) " mu = " , mu
#     endif

      !-----

      special_case : if ( mu == -1 ) then

        !Note: if mu = -1 then x <= 1, so the possible values for
        !lambda are -1 and 0. We already have S_0 if lambda = -1. 
        if ( lambda == 0 ) res_tmp(0) = 2.0_dp * ( 1.0_dp - x * res_tmp(1) )

        if ( n0 < 0 ) then

          !We need to compute S_{-1}
          if ( x * res_tmp(0) < 0.5_dp ) then
            !Recurrence relation
            res_tmp(-1) = 2.0_dp * ( x * res_tmp(0) - 1.0_dp ) 
          else if ( x < 0.461_dp ) then
            !Series expansion
            res_tmp(-1) = sm1_series1( x )
          else
            !Analytic continuation
            res_tmp(-1) = s_by_continuation( -1 , x )
          end if

          !Fill in the remaining values (if any) by recurrence
          do n = -2 , n0 , -1
            res_tmp(n) = ( 1.0_dp - x * res_tmp(n+1) ) / ( real( n , kind = dp ) + 0.5_dp )
          end do

        end if

        !-----

      else special_case

        x_recip = 1.0_dp / x

        !In this situation, we always need to calculate S( n , x ) at the
        !value of n that minimises | x - xn |, which is mu or mu + 1.
        nearest : if ( igz(mu,1) - x <= x - igz(mu+1,1) ) then

          res_tmp(mu) = s_by_continuation( mu , x )

          if ( n1 > mu ) then

            !In this case we also need S_{mu+1}

            !Note: both 2 * mu + 1 and s_mu are negative
            if ( real( 2 * mu + 1 , kind = dp ) * res_tmp(mu) < 1.0_dp ) then
              !Recurrence relation
              res_tmp(mu+1) = ( 1.0_dp - ( 0.5_dp + real( mu , kind = dp ) ) * res_tmp(mu) ) * x_recip
            else
              !Analytic continuation
              res_tmp(mu+1) = s_by_continuation( mu + 1 , x )
            end if

          end if

        else nearest

          res_tmp(mu+1) = s_by_continuation( mu + 1 , x )

          if ( n0 <= mu ) then

            !In this case we also need S_mu
          
            !Note: both x and S_{mu+1} are positive
            if ( x * res_tmp(mu+1) < 0.5_dp ) then
              !Recurrence relation
              res_tmp(mu) = ( 1.0_dp - x * res_tmp(mu+1) ) / ( real( mu , kind = dp ) + 0.5_dp )
            else
              !Analytic continuation
              res_tmp(mu) = s_by_continuation( mu , x )
            end if

          end if

        end if nearest

        !-----

        !Fill in remaining values by recurrence if necessary. 
        do n = mu - 1, n0 , -1
          res_tmp(n) = ( 1.0_dp - x * res_tmp(n+1) ) / ( 0.5_dp + real( n , kind = dp ) )
        end do

        do n = mu + 2 , up_rec_stop
          res_tmp(n) = ( 1.0_dp + ( 0.5_dp - real( n , kind = dp ) ) * res_tmp(n-1) ) * x_recip
        end do

        !-----

      end if special_case

    end if need_s_mu

    !-----

    allocate( res(n0:n1) , source = res_tmp(n0:n1) )

    !-----

  end subroutine

!-----

  real (dp) function sm1_series1( x ) result( s )
  !First series expansion [equation (12) in the revised paper] 
  !for S( -1 , x ). Note: this routine is only called when x <= 1, 
  !so there is no possibility of overstepping the recips array.

    implicit none

    real (dp) , intent (in) :: x
    
    !-----

    integer   :: j 

    real (dp) :: t , xjj

    !-----

    xjj =  2.0_dp ! 2 * x**j / j!
    s   = -2.0_dp

    j = 1
    do 

      xjj = xjj * x * recips(j)
      t   = xjj * recips(2*j-1)

      if ( t < -s * eps ) exit

      s = s + t
      j = j + 1

    end do

    s = s * exp( -x )

    !-----

  end function

!-----

  real (dp) function s_series2( n , x ) result( s )
  !Computes S_n(x) using the method in section 4 of the revised paper.

    implicit none

    integer   , intent (in) :: n

    real (dp) , intent (in) :: x

    !-----

    integer   :: j

    real (dp) :: t

    !-----

    !We will use a stack to facilitate backward summation.
    type :: stack_node 

      real (dp) :: t

      type (stack_node) , pointer :: prev

    end type

    type (stack_node) , target  :: t0

    type (stack_node) , pointer :: top , tmp

    !-----

#   ifdef DIAGNOSTIC
      write (*,*)
      write (*,100) n , x
      100 format( " Call to function s_series2 with  n = " , i5 , " and  x = " , f15.10 , " ..." )
#   endif

    !First term in the series
    t0 % t = 1.0_dp / ( real( n , kind = dp ) + 0.5_dp )

    !Initialise stack
    nullify( t0 % prev )
    top => t0

    !Note: provided that n >= 0, a bound on the relative error can be obtained 
    !by multiplying the current term by x and taking the modulus. The next term 
    !can then be obtained by dividing the result by n + j + 1/2. Therefore in 
    !the loop below we will update t in two steps.
    t = x * t0 % t

    j = 1

    do

      if ( abs( t ) < eps ) exit

      t = -t / ( real( n + j , kind = dp ) + 0.5_dp )

      !Create new node
      allocate( tmp )
      tmp % prev => top
      tmp % t    = t

      !New node is now stack top
      top => tmp

      t = t * x
      j = j + 1

    end do

    !-----

    !Sum the values on the stack, and clean house
    s = top % t

    do

      if ( .not. associated( top % prev ) ) exit

      tmp => top
      top => top % prev

      s = s + top % t

      deallocate( tmp )

    end do

    !-----

#   ifdef DIAGNOSTIC
      write (*,300) s , j + 1
      300 format ( " Approximate value of S is " , f25.10 , ";" , i5 , " terms summed" )
#   endif

    !-----

  end function

!-----

  real (dp) function s_by_continuation( n , x ) result( s )
  !Uses analytic continuation and precomputed data to calculate S( n , x ).
  !For any given x, there are, at most, two appropriate values of n. These 
  !are calculated by a higher level routine.

    implicit none

    integer   , intent (in) :: n

    real (dp) , intent (in) :: x

    !-----

    integer   :: j , p , p_max
    integer   :: nsteps
    integer   :: an 

    real (dp) :: t , dt , mdt , dtt
    real (dp) :: xp
    real (dp) :: t_recip
    real (dp) :: max_step_recip
    real (dp) :: mdt_pow
    real (dp) :: w , g 
    real (dp) :: f , fst
    real (dp) :: upd
    real (dp) :: y , tmp , corr

    real (dp) , dimension (dg+2) :: gamma_vals

    !-----

#   ifdef DIAGNOSTIC
      write (*,*)
      write (*,100) n , x
      100 format( " Call to function s_by_continuation with  n = " , i5 , " and x = " , f15.10 )
#   endif

    !-----
 
    t  = igz(n,1) !Starting point for stepping algorithm
    s  = igz(n,2) !Value of s at starting point

    t_recip = 1.0_dp / t

#   ifdef DIAGNOSTIC
      write(*,fmt='( " Initial value: t = " , f20.15 , " s = " , e30.20 )') t , s
#   endif

    !-----

    !Calculate delta t.

    !Are we increasing or decreasing t with each step?
    if ( x > t ) then

      !In this case, the maximum allowable step size increases with each 
      !successive step, but for efficiency we will use steps of equal
      !size. Therefore we must consider the size of the first step.
      max_step_recip = real( 1 - 2 * n , kind = dp ) * t_recip
      nsteps         = ceiling( ( x - t ) * max_step_recip )

    else

      !Here, the maximum allowable step size decreases with each successive
      !step. The maximum is calculated by solving for the largest safe
      !step that reaches x.
      max_step_recip = -real( 2 * n , kind = dp ) / x
      nsteps         = ceiling( ( t - x ) * max_step_recip )

    end if

    dt  = ( x - t ) * recips(nsteps)
    mdt = -dt
    xp  = exp( dt )

    !-----

#   ifdef DIAGNOSTIC
      write (*,fmt='( " dt = " , f15.12 , "  steps required: " , i2 )') dt , nsteps
#   endif

    !-----

    !This is an upper bound on the number of terms required by the update series
    !(0.69315 is slightly greater than log( 2 )). Note that we must calculate one
    !more term than we actually use.
    p_max = ceiling( -real( dg , kind = dp ) * 0.69315_dp &
        & / log( abs( dt ) * t_recip * ( 0.5_dp - real( n , kind = dp ) ) ) ) + 1 

    !Calculate gamma( p_max , -dt ) / (-dt)**p_max using equation (7) in the revised
    !paper (good for large p_max). Note that we will always have |dt| < 1, so
    !there is no danger of overstepping the recips array.

    w = recips(p_max) ! (-dt )**j / ( p_max )_{j+1}

    g = w

    p = 1
    gmx : do

      w = w * mdt * recips(p_max+p)

      if ( abs( w ) < abs( g ) * eps ) exit gmx

      g = g + w
      p = p + 1

    end do gmx

    gamma_vals(p_max) = g * xp

    !-----

    !Obtain lower values by recurrence
    do p = p_max - 1 , 1 , -1
      gamma_vals(p) = ( xp - dt * gamma_vals(p+1) ) * recips(p)
    end do

    !Multiply by appropriate power of (-dt)
    mdt_pow       = mdt
    gamma_vals(1) = gamma_vals(1) * mdt_pow 

    do p = 2 , p_max
      mdt_pow = mdt * mdt_pow
      gamma_vals(p) = gamma_vals(p) * mdt_pow
    end do

    !-----

#   ifdef SUPERDIAGNOSTIC
      do p = 1 , p_max
        write(*,fmt='( " gamma( " , i4 , " " , f20.15 , " ) = " , f20.15 )') p , mdt , gamma_vals(p)
      end do
#   endif

    !-----

    step : do j = 1 , nsteps

      !Calculate update series; each term has magnitude less than 1/2 that of the
      !previous term, so there is no possibility of overstepping the recips array.
      upd = gamma_vals(1)
      f   = 1.0_dp    ! (1/2-n)_p / ( p! * t**p )
      p   = 1

#     ifdef SUPERDIAGNOSTIC

        write (*,fmt='( "Calculating update series: " )')
        write (*,fmt='( i3 , 2( " " , f25.20 ) )') 1 , upd , upd

#     endif

      update : do 

       f   = f * ( 0.5_dp - real( n + 1 - p , kind = dp ) ) * recips(p) * t_recip
       w   = f * gamma_vals(p+1)

        if ( abs( w ) < eps * abs( upd ) ) exit update

        upd = upd + w

#       ifdef SUPERDIAGNOSTIC
          write (*,fmt='( i3 , 2( " " , f25.20 ) )') p , w , upd
#       endif

        p   = p + 1

      end do update

      !-----

      !Calculate ( 1 + dt / t )**(-n) 
      !The binomial expansion + Kahan summation seems to give the best accuracy.
      dtt  = dt * t_recip

#     ifdef DIAGNOSTIC
        write (*,fmt='( " dt / t = " , f20.15 )') dtt
#     endif

      f    = 1.0_dp
      fst  = 1.0_dp
      corr = 0.0_dp

      an   = -n
       
      do p = 1 , an

        !Calculate next term in series
        fst  = fst * dtt * real( an - p + 1 , kind = dp ) / real( p , kind = dp )

        !Calculate correction
        y    = fst - corr
        tmp  = f + y
        corr = ( tmp - f ) - y

        !Update series
        f    = tmp

        !Watch out for floating underflow. Note that the maximum value for ( an - p + 1 ) / p 
        !is -n, and ( -n dt / t ) <= 0.5, so the series is monotonic.
        if ( abs( fst ) < f * eps ) exit

      end do

      f = f / ( sqrt( 1.0_dp + dtt ) )

      !Update s
      s       = f / xp * ( s - upd * t_recip )

      !Update t
      t       = t + dt
      t_recip = 1.0_dp / t

      !-----

#     ifdef DIAGNOSTIC
        write (*,fmt='( "Update  : t = " , f20.15 , " s = " , e30.20 )') t , s
#     endif

      !-----

    end do step

    !-----

#   ifdef DIAGNOSTIC
      write (*,*)
      write (*,fmt='( " function s_by_continuation executed ok. Return value is " , e30.20 )') s
#   endif

    !-----

  end function

!-----

  subroutine read_precomp_data()

    implicit none

    !-----

    logical :: in_use

    integer :: ios , u
    integer :: j
    integer :: d

    integer (li) :: x , gx

    !-----

#   ifdef DIAGNOSTIC
      write (*,*)
      write (*,fmt='( " Call to read_precomp_data..." )')
#   endif

    !-----
 
    !Find a safe unit number. In Fortran 2008
    !we could use newunit in the open statement.
    u = 10

    do
      inquire( unit = u , opened = in_use )
      if ( .not. in_use ) exit
      u = u + 1
    end do
 
    !-----
 
    open( unit = u , file = "in_gamma_precomp.dat" , action = 'read' , iostat = ios )
    if ( ios /= 0 ) call error( "read_precomp_data" , "unable to open file" )

    !Skip three header lines
    do j = 1 , 3
      read (u,*,iostat=ios)
      if ( ios /= 0 ) call error( "read_precomp_data" , "invalid data file structure" )
    end do

    read (u,*,iostat=ios) mu_min
    if ( ios /= 0 ) call error( "read_precomp_data" , "invalid data file structure" )

    backspace( u )

    allocate( igz(mu_min:0,2) )

    !Note: the paper does not define x0, but it's convenient to set it
    !to an arbitrary negative value here. (See the need_s_mu block in
    !subroutine scaled_in_gamma for details.)
    igz(0,1) = -1.0_dp

    !-----

    do j = mu_min , -1

      read (u,*,iostat=ios) d , x , gx
      if ( ios /= 0 ) call error( "read_precomp_data" , "incomplete data file" )

      igz(j,:) = [ transfer( x , 1.0_dp ) , transfer( gx , 1.0_dp ) ]

    end do

    !d doesn't actually do anything; this just prevents a compiler warning
    d = d - 1 

    !-----

#   ifdef DIAGNOSTIC
      write (*,fmt='( " Data obtained for mu = " , i5 , " ... -1" )') mu_min
#   endif

    !-----

  end subroutine

!-----

end module
