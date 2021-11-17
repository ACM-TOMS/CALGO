subroutine DTIMSWAP( n, nrep, x, y, xs, ys, thresh, rtim, nout )
   use LA_CONSTANTS, only: wp, zero, one, ulp
!
!  .. Scalar Arguments ..
   integer :: n, nrep, nout
   real(wp) :: rtim, thresh
!  ..
!  .. Array Arguments ..
   real(wp) :: x(*), xs(*), y(*), ys(*)
!  ..
!
!  DTIMSWAP times DSWAP from the Level 1 BLAS, a subroutine to swap
!  two vectors.  This is the in-cache test with the same vectors used
!  for each trial.
!
!  E. Anderson
!  July 8, 2016
!
!  .. Local Scalars ..
   character*6 :: subnam
   integer :: i, imax, incx, incxs, incy, incys, ix, iy, &
      j, jmax, nerrs, nfail, ns, ntests
   real(wp) :: t1, t2, tottim, tmax, trat, umax, urat, &
      xexp, xmax, xnrm, yexp, ymax, ynrm
!  ..
!  .. External Functions ..
   real(wp) :: DSECND
   external :: DSECND
!  ..
!
   nfail = 0
   ntests = 0
   subnam = 'DSWAP '
   ns = n
   incxs = 1
   incys = 1
!
!  Set x(i) = real(i) and y(i) = real(i)/1000
!
   do i = 1, n
     xs(i) = real(i)
     ys(i) = real(i) / 1000._wp
   end do
!
!  Do the test nrep times but only test the final one
!
   n = ns
   x(1:n) = xs(1:n)
   incx = incxs
   y(1:n) = ys(1:n)
   incy = incys
   if( nrep == 1 ) then
     jmax = 2
   else
     jmax = nrep - mod(nrep,2)
   end if
   t1 = DSECND()
   do j = 1, jmax-1, 2
     call DSWAP( n, x, incx, y, incy )
     call DSWAP( n, x, incx, y, incy )
   end do
   t2 = DSECND()
   tottim = t2 - t1
   rtim = tottim / jmax
!
   n = ns
   x(1:n) = xs(1:n)
   incx = incxs
   y(1:n) = ys(1:n)
   incy = incys
   call DSWAP( n, x, incx, y, incy )
!  write(nout,99) 'DSWAP ', ns, rtim
!
!  See what values changed inside the subroutine
!
   nerrs = 0
   call IARGCKS( nerrs, subnam, 'N', n, ns, nout )
   call IARGCKS( nerrs, subnam, 'INCX', incx, incxs, nout )
   call IARGCKS( nerrs, subnam, 'INCY', incy, incys, nout )
!
!  Check that the expected values are in x
!
   imax = 1
   xmax = zero
   tmax = -one
   ix = 1
   if( incxs < 0 ) ix = 1 - (ns-1)*incxs
   do i = 1, ns
     xexp = ys(i)
     xnrm = abs( xexp )
     if( xnrm == zero ) then
       trat = abs( x(ix) ) / ulp
     else
       trat = ( abs( x(ix) - xexp ) / xnrm ) / ulp
     end if
     if( trat > tmax ) then
       imax = i
       xmax = x(ix)
       tmax = trat
     end if
     ix = ix + incxs
   end do
!
!  Check that the expected values are in y
!
   jmax = 1
   ymax = zero
   umax = -one
   iy = 1
   if( incys < 0 ) iy = 1 - (ns-1)*incys
   do i = 1, ns
     yexp = xs(i)
     ynrm = abs( yexp )
     if( ynrm == zero ) then
       urat = abs( y(iy) ) / ulp
     else
       urat = ( abs( y(iy) - yexp ) / ynrm ) / ulp
     end if
     if( urat > umax ) then
       jmax = i
       ymax = y(iy)
       umax = urat
     end if
     iy = iy + incys
   end do
   if( max(tmax,umax) >= thresh ) then
     write(nout,98) subnam, ns, incxs, incys, max(tmax,umax)
     xexp = ys(imax)
     yexp = xs(jmax)
     write(nout,97) imax, xmax, xexp, jmax, ymax, yexp
     nerrs = nerrs + 1
   end if
   if( nerrs > 0 ) nfail = nfail + 1
   ntests = ntests + 1
!
   if( nfail /= 0 ) then
     write(nout,95) subnam, nfail, ntests
   end if
   return
!
99 format( A6, ': N=', I8, ', time=', E15.8 )
98 format( A6, ': N=', I6, ', INCX=', I4, ', INCY=', I4, ', test=', E15.8 )
97 format( '   occurs at X(', I6, ')=', E15.8, ', expected ', E15.8, /, &
           '   and/or at Y(', I6, ')=', E15.8, ', expected ', E15.8 )
96 format( A6, ' passed the computational tests (', I6, ' tests run)', / )
95 format( A6, ': ', I6, ' of ', I6,' tests failed', / )
end subroutine
