function SLAPY2( x, y )
   use LA_CONSTANTS32, only: wp, zero, one, rtmin, rtmax, safmin, safmax
   real(wp) :: SLAPY2
!
!  LAPACK auxiliary routine
!  E. Anderson
!  July 30, 2016
!
!  .. Scalar Arguments ..
   real(wp) :: x, y
!  ..
!
!  Purpose
!  =======
!
!  SLAPY2 returns sqrt(x**2+y**2), taking care not to cause unnecessary
!  overflow.
!
!  Arguments
!  =========
!
!  X       (input) REAL
!  Y       (input) REAL
!          X and Y specify the values x and y.
!
! =====================================================================
!
!  .. Local Scalars ..
   real(wp) :: w, ww, xabs, yabs
!  ..
!
   xabs = abs( x )
   yabs = abs( y )
   if( xabs == zero ) then
      SLAPY2 = yabs
   else if( yabs == zero ) then
      SLAPY2 = xabs
   else if( xabs > rtmin .and. xabs < rtmax .and. &
            yabs > rtmin .and. yabs < rtmax ) then
      SLAPY2 = sqrt( x*x + y*y )
   else
      w = min( safmax, max( safmin, xabs, yabs ) )
      ww = one / w
      SLAPY2 = w*sqrt( (x*ww)**2 + (y*ww)**2 )
   end if
   return
end function
