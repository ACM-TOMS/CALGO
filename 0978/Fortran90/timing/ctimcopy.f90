subroutine CTIMCOPY( n, nrep, x, y, xs, ys, thresh, rtim, nout )
   use LA_CONSTANTS32, only: wp, zero, one, ulp, czero
!
!  .. Scalar Arguments ..
   integer :: n, nrep, nout
   real(wp) :: rtim, thresh
!  ..
!  .. Array Arguments ..
   complex(wp) :: x(*), xs(*), y(*), ys(*)
!  ..
!
!  CTIMCOPY times CCOPY from the Level 1 BLAS, a subroutine to copy
!  a vector x to a vector y.  This is the in-cache test with the same
!  vectors used for each trial.
!
!  E. Anderson
!  July 11, 2016
!
!  .. Local Scalars ..
   character*6 :: subnam
   integer :: i, imax, incx, incxs, incy, incys, iy, &
      j, nerrs, nfail, ns, ntests
   real(wp) :: rogue, t1, t2, tottim, tmax, trat, ynrm
   complex(wp) :: yexp, ymax
!  ..
!  .. External Functions ..
   real(wp) :: SECOND
!  ..
!
   nfail = 0
   ntests = 0
   rogue = -999.999_wp
   subnam = 'CCOPY '
   ns = n
   incxs = 1
   incys = 1
!
!  Set x(i) = cmplx( real(i), -real(i) )
!
   do i = 1, n
     xs(i) = cmplx( real(i), -real(i), wp )
   end do
   ys(1:n) = cmplx( rogue, -rogue, wp )
!
!  Do the test nrep times but only test the final one
!
   n = ns
   x(1:n) = xs(1:n)
   incx = incxs
   y(1:n) = ys(1:n)
   incy = incys
   t1 = SECOND()
   do j = 1, nrep
     call CCOPY( n, x, incx, y, incy )
   end do
   t2 = SECOND()
   tottim = t2 - t1
   rtim = tottim / nrep
!  write(nout,99) 'CCOPY ', ns, rtim
!
!  See what values changed inside the subroutine
!
   nerrs = 0
   call IARGCKS( nerrs, subnam, 'N', n, ns, nout )
   call CARGCKV( nerrs, subnam, 'X', ns, x, xs, nout )
   call IARGCKS( nerrs, subnam, 'INCX', incx, incxs, nout )
   call IARGCKS( nerrs, subnam, 'INCY', incy, incys, nout )
!
!  Check that the expected values are in y
!
   iy = 1
   if( incys < 0 ) iy = 1 - (ns-1)*incys
   imax = 1
   ymax = czero
   tmax = -one
   do i = 1, ns
     yexp = xs(i)
     ynrm = abs( yexp )
     if( ynrm == zero ) then
       trat = abs( y(iy) ) / ulp
     else
       trat = ( abs( y(iy) - yexp ) / ynrm ) / ulp
     end if
     if( trat > tmax ) then
       imax = i
       ymax = y(iy)
       tmax = trat
     end if
     iy = iy + incys
   end do
   if( tmax >= thresh ) then
     write(nout,98) subnam, ns, incxs, incys, tmax
     yexp = xs(imax)
     write(nout,97) imax, ymax, yexp
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
98 format( A6, ': N=', I8, ', INCX=', I4, ', INCY=', I4, ', test=', E15.8 )
97 format( '   occurs at Y(', I8, ')=(', E15.8, ',', E15.8, &
           '), expected=(', E15.8, ',', E15.8, ')' )
96 format( A6, ' passed the computational tests (', I6, ' tests run)', / )
95 format( A6, ': ', I6, ' of ', I6,' tests failed', / )
end subroutine
