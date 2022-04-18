subroutine DTIMPY2( n, nrep, x, y, z, xs, ys, thresh, rtim, nout )
   use LA_CONSTANTS, only: wp, zero, ulp, one, two
!
!  .. Scalar Arguments ..
   integer :: n, nrep, nout
   real(wp) :: rtim, thresh
!  ..
!  .. Array Arguments ..
   real(wp) :: x(*), xs(*), y(*), ys(*), z(*)
!  ..
!
!  DTIMPY2 times DLAPY2 from LAPACK, a function to compute the
!  quantity sqrt( x**2 + y**2 ) with scaling to avoid overflow or
!  underflow.  This is the in-cache test with the same vector
!  used for each trial.
!
!  E. Anderson
!  July 8, 2016
!
!  .. Local Scalars ..
   character*6 :: subnam
   integer :: i, imax, iy, j, k, nerrs, nfail, ns, ntests
   real(wp) :: rogue, tmax, trat, zexp, zmax, t1, t2, tottim
!  ..
!  .. External Functions ..
   real(wp) :: DSECND, DLAPY2
!  ..
!
   nfail = 0
   ntests = 0
   rogue = -ONE
   subnam = 'DLAPY2'
   ns = n
!
!  Set xs and ys to random values from (-1,1).
!
   call random_number(xs(1:n))
   do i = 1, n
     xs(i) = one - two*xs(i)
   end do
   call random_number(ys(1:n))
   do i = 1, n
     ys(i) = one - two*ys(i)
   end do
!
!  Do the test nrep times but only test the final one
!
   n = ns
   x(1:n) = xs(1:n)
   y(1:n) = ys(1:n)
   t1 = DSECND()
   do j = 1, nrep
     do k = 1, n
       z(k) = DLAPY2( x(k), y(k) )
     end do
   end do
   t2 = DSECND()
   tottim = t2 - t1
   rtim = tottim / nrep
!  write(nout,99) 'DLAPY2', ns, tottim / nrep
!
!  See what values changed inside the subroutine
!
   nerrs = 0
   call IARGCKS( nerrs, subnam, 'N', n, ns, nout )
   call DARGCKV( nerrs, subnam, 'X', n, x, xs, nout )
   call DARGCKV( nerrs, subnam, 'Y', n, y, ys, nout )
!
!    Check that the expected values are in z.
!    In case of failure, only report the largest one.
!
   imax = 1
   zmax = zero
   tmax = -one
   iy = 1
   do i = 1, ns
     zexp = sqrt( xs(i)**2 + ys(i)**2 )
     if( zexp == zero ) then
       trat = abs( z(i) ) / ulp
     else
       trat = ( abs( z(i) - zexp ) / zexp ) / ulp
     end if
     if( trat > tmax ) then
       imax = i
       zmax = z(i)
       tmax = trat
     end if
   end do
   if( tmax >= thresh ) then
     write(nout,98) subnam, ns, tmax
     zexp = sqrt( xs(imax)**2 + ys(imax)**2 )
     write(nout,97) imax, zmax, zexp
     nerrs = nerrs + 1
   end if
   if( nerrs > 0 ) nfail = nfail + 1
   ntests = ntests + 1
   if( nfail /= 0 ) then
     write(*,95) subnam, nfail, ntests
   end if
   return
!
99 format( A6, ': N=', I8, ', time=', E15.8 )
98 format( A6, ': N=', I8, ', test=', E15.8 )
97 format( '   occurs at Z(', I8, ')=', E15.8, ', expected=', E15.8 )  
96 format( A6, ' passed the computational tests (', I6, ' tests run)', / )
95 format( A6, ': ', I6, ' of ', I6,' tests failed', / )
end subroutine
