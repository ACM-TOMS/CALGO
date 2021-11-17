module rk_erfcx_Cody
    ! This module provides generic interfaces for functions that evaluate
    ! erf(x), erfc(x), and exp(x*x)*erfc(x) in either single or double precision.
   use set_rk, only: rk
    implicit none
    private
    public :: rk_Erf, rk_Erfc, rk_Erfcx

    interface rk_erf
        module procedure rk_erf
    end interface

    interface rk_erfc
        module procedure rk_erfc
    end interface

    interface rk_erfcx
        module procedure rk_erfcx
    end interface


contains

    elemental function rk_erf ( x )

        !*****************************************************************************80
        !
        !! rk_ERF evaluates the error function.
        !
        !  Discussion:
        !
        !    This routine computes approximate values for erf(x).
        !
        !    This routine was renamed from "DERF" to "rk_ERF" to avoid being
        !    overshadowed by corresponding functions supplied by some compilers.
        !
        !    See comments heading rk_CALERF.
        !
        !  Modified:
        !
        !    23 January 2008
        !
        !  Author:
        !
        !    William Cody
        !
        !  Parameters:
        !
        !    Input, real ( kind = rk ) X, the argument of the function.
        !
        !    Output, real ( kind = rk ) rk_ERF, the value of the function.
        !
        implicit none
   !     integer, parameter :: rk = rpr  ! double
        real ( rk ), intent(in) :: x
        integer :: jint
        real    ( rk ) rk_erf

        jint = 0
        call rk_calerf ( x, rk_erf, jint )

    end function rk_erf



    elemental function rk_erfc ( x )

        !*****************************************************************************80
        !
        !! rk_ERFC evaluates the complementary error function.
        !
        !  Discussion:
        !
        !    This routine computes approximate values for erfc(x).
        !
        !    This routine was renamed from "DERFC" to "rk_ERFC" to avoid being
        !    overshadowed by corresponding functions supplied by some compilers.
        !
        !    See comments heading CALERF.
        !
        !  Modified:
        !
        !    23 January 2008
        !
        !  Author:
        !
        !    William Cody
        !
        !  Parameters:
        !
        !    Input, real ( kind = rk ) X, the argument of the function.
        !
        !    Output, real ( kind = rk ) rk_ERFC, the value of the function.
        !
        implicit none
  !      integer, parameter :: rk = rpr  ! double
        integer::jint
        real    ( rk ), intent(in) :: x
        real    ( rk ) rk_erfc
        jint = 1
        call rk_calerf ( x, rk_erfc, jint )

    end function rk_erfc



    elemental function rk_erfcx ( x )

        !*****************************************************************************80
        !
        !! rk_ERFCX evaluates the exponentially scaled complementary error function.
        !
        !  Discussion:
        !
        !    This routine computes approximate values for exp(x*x) * erfc(x).
        !
        !    This routine was renamed from "DERFCX" to "rk_ERFCX" to match the
        !    renamings of DERF and DERFC.
        !
        !    See comments heading CALERF.
        !
        !  Modified:
        !
        !    23 January 2008
        !
        !  Author:
        !
        !    William Cody
        !
        !  Parameters:
        !13 significant digits
        !    Input, real ( kind = rk ) X, the argument of the function.
        !
        !    Output, real ( kind = rk ) rk_ERFCX, the value of the function.
        !
        implicit none
 !       integer, parameter :: rk = rpr  ! double
        real ( rk ), intent(in) :: x
        real ( rk ) rk_erfcx
        integer :: jint
        jint = 2
        call rk_calerf ( x, rk_erfcx, jint )

    end function rk_erfcx


    elemental subroutine rk_calerf ( arg, result, jint )

        !*****************************************************************************80
        !
        !! CALERF computes various forms of the error function.
        !
        !  Discussion:
        !
        !    This routine evaluates erf(x), erfc(x), and exp(x*x)*erfc(x)
        !    for a real argument x.
        !
        !  Modified:
        !
        !    03 April 2007
        !
        !  Author:ifort rk_erfcx_Cody.f90  Faddeyeva_v2_mod_rk.f90 Faddeyeva_driver_rk.f90 -o Faddeyeva_driver_rk
        !
        !    William Cody
        !
        !  Reference:
        !
        !    William Cody,
        !    Rational Chebyshev Approximations for the Error Function,
        !    Mathematics of Computation,
        !    Volume 23, Number 107, July 1969, pages 631-638.
        !
        !  Parameters:
        !
        !    Input, real ( kind = rk ) ARG, the argument.  If JINT is 1, the
        !    argument must be less than XBIG.  If JINT is 2, the argument
        !    must lie between XNEG and XMAX.
        !
        !    Output, real ( kind = rk ) RESULT, the value of the function,
        !    which depends on the input value of JINT:
        !    0, RESULT = erf(x);
        !    1, RESULT = erfc(x) = 1 - erf(x);
        !    2, RESULT = exp(x*x)*erfc(x) = exp(x*x) - erf(x*x)*erf(x).
        !
        !    Input, integer ( kind = i4 ) JINT, chooses the function to be computed.
        !    0, erf(x);
        !    1, erfc(x);
        !    2, exp(x*x)*erfc(x).
        !
        implicit none
        integer, intent(in)::jint
        real    ( rk ), intent(in) :: arg
        real    ( rk ), intent(out) :: result

        real    ( rk ) del
        real    ( rk ), parameter :: four = 4.0e0_rk
        real    ( rk ), parameter :: half = 0.5e0_rk
        integer :: i
        real    ( rk ), parameter :: one = 1.0e0_rk
        real    ( rk ), parameter :: sixten = 16.0e0_rk
        real( rk ), parameter :: sqrpi = 5.641895835477563e-001_rk
        real    ( rk ), parameter :: two = 2.0e0_rk
        real    ( rk ), parameter :: thresh = 0.46875e0_rk
        real    ( rk ) x
        real    ( rk ) xden
        real    ( rk ) xnum
        real    ( rk ) y
        real    ( rk ) ysq
        real    ( rk ), parameter :: zero = 0.0e0_rk

        !  Machine-dependent constants

        real    ( rk ), parameter :: xinf = huge(1.0_rk)
        real    ( rk ), parameter :: xneg = -26.628_rk
        real    ( rk ), parameter :: xsmall = 1.11e-16_rk
        real    ( rk ), parameter :: xbig = 26.543_rk
        real    ( rk), parameter :: xhuge = 6.71e7_rk
        real    ( rk ), parameter :: xmax = sqrpi/tiny(1.0_rk)!

        !  Coefficients for approximation to  erf  in first interval

        real    ( rk ), parameter :: a(5) = &
            (/ 3.16112374387056560e00_rk,1.13864154151050156e02_rk, &
            3.77485237685302021e02_rk,3.20937758913846947e03_rk, &
            1.85777706184603153e-1_rk /)
        real    ( rk ), parameter :: b(4) = &
            (/ 2.36012909523441209e01_rk,2.44024637934444173e02_rk, &
            1.28261652607737228e03_rk,2.84423683343917062e03_rk /)

        !  Coefficients for approximation to  erfc  in second interval

        real    ( rk ), parameter :: c(9) = &
            (/ 5.64188496988670089e-1_rk,8.88314979438837594e0_rk, &
            6.61191906371416295e01_rk,2.98635138197400131e02_rk, &
            8.81952221241769090e02_rk,1.71204761263407058e03_rk, &
            2.05107837782607147e03_rk,1.23033935479799725e03_rk, &
            2.15311535474403846e-8_rk /)
        real    ( rk ), parameter :: d(8) = &
            (/ 1.57449261107098347e01_rk,1.17693950891312499e02_rk, &
            5.37181101862009858e02_rk,1.62138957456669019e03_rk, &
            3.29079923573345963e03_rk,4.36261909014324716e03_rk, &
            3.43936767414372164e03_rk,1.23033935480374942e03_rk /)

        !  Coefficients for approximation to  erfc  in third interval

        real    ( rk ), parameter :: p(6) = &
            (/ 3.05326634961232344e-1_rk,3.60344899949804439e-1_rk, &
            1.25781726111229246e-1_rk,1.60837851487422766e-2_rk, &
            6.58749161529837803e-4_rk,1.63153871373020978e-2_rk /)
        real    ( rk ), parameter :: q(5) = &
            (/ 2.56852019228982242e00_rk,1.87295284992346047e00_rk, &
            5.27905102951428412e-1_rk,6.05183413124413191e-2_rk, &
            2.33520497626869185e-3_rk /)

        x = arg
        y = abs ( x )

        !  Evaluate erf for |X| <= 0.46875.

        if ( y <= thresh ) then

            ysq = zero
            if ( xsmall < y ) then
                ysq = y * y
            end if

            xnum = a(5) * ysq
            xden = ysq

            do i = 1, 3
                xnum = ( xnum + a(i) ) * ysq
                xden = ( xden + b(i) ) * ysq
            end do

            result = x * ( xnum + a(4) ) / ( xden + b(4) )

            if ( jint /= 0 ) then
                result = one - result
            end if

            if ( jint == 2 ) then
                result = exp ( ysq ) * result
            end if

            return

            !  Evaluate erfc for 0.46875 <= |X| <= 4.0.

        else if ( y <= four ) then

            xnum = c(9) * y
            xden = y

            do i = 1, 7
                xnum = ( xnum + c(i) ) * y
                xden = ( xden + d(i) ) * y
            end do

            result = ( xnum + c(8) ) / ( xden + d(8) )

            if ( jint /= 2 ) then
                ysq = aint ( y * sixten ) / sixten
                del = ( y - ysq ) * ( y + ysq )
                result = exp ( -ysq * ysq ) * exp ( -del ) * result
            end if

            !  Evaluate erfc for 4.0 < |X|.

        else

            result = zero

            if ( xbig <= y ) then

                if ( jint /= 2 .or. xmax <= y ) then
                    go to 300
                end if

                if ( xhuge <= y ) then
                    result = sqrpi / y
                    go to 300
                end if

            end if

            ysq = one / ( y * y )
            xnum = p(6) * ysq
            xden = ysq
            do i = 1, 4
                xnum = ( xnum + p(i) ) * ysq
                xden = ( xden + q(i) ) * ysq
            end do

            result = ysq * ( xnum + p(5) ) / ( xden + q(5) )
            result = ( sqrpi -  result ) / y

            if ( jint /= 2 ) then
                ysq = aint ( y * sixten ) / sixten
                del = ( y - ysq ) * ( y + ysq )
                result = exp ( -ysq * ysq ) * exp ( -del ) * result
            end if

        end if

        !  Fix up for negative argument, erf, etc.

        300 continue

        if ( jint == 0 ) then

            result = ( half - result ) + half
            if ( x < zero ) then
                result = -result
            end if

        else if ( jint == 1 ) then

            if ( x < zero ) then
                result = two - result
            end if

        else

            if ( x < zero ) then

                if ( x < xneg ) then
                    result = xinf
                else
                    ysq = aint ( x * sixten ) / sixten
                    del = ( x - ysq ) * ( x + ysq )
                    y = exp ( ysq * ysq ) * exp ( del )
                    result = ( y + y ) - result
                end if

            end if

        end if

    end subroutine rk_calerf

end module rk_erfcx_Cody
