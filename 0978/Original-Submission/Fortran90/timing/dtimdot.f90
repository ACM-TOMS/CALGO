subroutine DTIMDOT( n, nrep, x, y, xs, ys, thresh, rtim, nout )
   use LA_CONSTANTS, only: wp, zero, one, two, ulp
!
!  .. Scalar Arguments ..
   integer :: n, nrep, nout
   real(wp) :: rtim, thresh
!  ..
!  .. Array Arguments ..
   real(wp) :: x(*), xs(*), y(*), ys(*)
!  ..
!
!  DTIMDOT times DDOT from the Level 1 BLAS, a function to compute the
!  dot product of a vector x and a vector y.  This is the in-cache test
!  with the same vectors used for each trial.
!
!  E. Anderson
!  July 8, 2016
!
!  .. Local Scalars ..
   character*6 :: subnam
   integer :: i, incx, incxs, incy, incys, j, nerrs, nfail, ns, ntests
   real(wp) :: rogue1, rogue2, sexp, snrm, trat, t1, t2, tottim
!  ..
!  .. Local Arrays ..
   real(wp), allocatable :: sbla(:)
!  ..
!  .. External Functions ..
   real(wp) :: DDOT, DSECND
   external :: DDOT, DSECND
!  ..
!
   nfail = 0
   ntests = 0
   rogue1 = -666.666_wp
   rogue2 = -999.999_wp
   subnam = 'DDOT  '
   ns = n
   incxs = 1
   incys = 1
   allocate(sbla(nrep))
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
!  Calculate the expected value of the dot product
!
   sexp = zero
   do i = 1, n
     sexp = sexp + xs(i)*ys(i)
   end do
!
!  Do the test nrep times but only test the final one
!
   n = ns
   x(1:n) = xs(1:n)
   incx = incxs
   y(1:n) = ys(1:n)
   incy = incys
   t1 = DSECND()
   do j = 1, nrep
     sbla(j) = DDOT( n, x, incx, y, incy )
   end do
   t2 = DSECND()
   tottim = t2 - t1
   rtim = tottim / nrep
!  write(nout,99) 'DDOT  ', ns, rtim
!
!  See what values changed inside the subroutine
!
   nerrs = 0
   call IARGCKS( nerrs, subnam, 'N', n, ns, nout )
   call DARGCKV( nerrs, subnam, 'X', ns, x, xs, nout )
   call IARGCKS( nerrs, subnam, 'INCX', incx, incxs, nout )
   call DARGCKV( nerrs, subnam, 'Y', ns, y, ys, nout )
   call IARGCKS( nerrs, subnam, 'INCY', incy, incys, nout )
!
!  Compare sbla and sexp
!
   snrm = abs( sexp )
   if( snrm == zero ) then
     trat = abs( sbla(nrep) ) / ulp
   else
     trat = ( abs( sbla(nrep) - sexp ) / snrm ) / ( real(max(1,n),wp)*ulp )
   end if
   if (trat >= thresh) then
     write(nout,98) subnam, ns, incxs, incys, trat
     nerrs = nerrs + 1
   end if
   if( nerrs > 0 ) nfail = nfail + 1
   ntests = ntests + 1
!
   if( nfail /= 0 ) then
     write(nout,96) subnam, nfail, ntests
   end if
   deallocate(sbla)
   return
99 format( A6, ': N=', I8, ', time=', E15.8 )
98 format( A6, ': N=', I8,', INCX=', I4, ', INCY=', I4, ', test=', E15.8 )
97 format( A6, ' passed the computational tests (', I6, ' tests run)', / )
96 format( A6, ': ', I6, ' of ', I6,' tests failed to pass threshold', / )
end subroutine
