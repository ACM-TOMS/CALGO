subroutine STIMASUM( n, nrep, x, xs, thresh, rtim, nout )
   use LA_CONSTANTS32, only: wp, zero, one, two, ulp
!
!  .. Scalar Arguments ..
   integer :: n, nrep, nout
   real(wp) :: rtim, thresh
!  ..
!  .. Array Arguments ..
   real(wp) :: x(*), xs(*)
!  ..
!
!  STIMASUM times SASUM from the Level 1 BLAS, a function to compute
!  the 1-norm of a vector x.  This is the in-cache test with the same
!  vector used for each trial.
!
!  E. Anderson
!  July 8, 2016
!
!  .. Local Scalars ..
   character*6 :: subnam
   integer :: i, incx, incxs, j, nerrs, nfail, ns, ntests
   real(wp) :: t1, t2, tottim, trat, znrm
!  ..
!  .. Local Arrays ..
   real(wp), allocatable :: sbla(:)
!  ..
!  .. External Functions ..
   real(wp) :: SASUM, SECOND
   external :: SASUM, SECOND
!  ..
!
   nfail = 0
   ntests = 0
   subnam = 'SASUM '
   ns = n
   incxs = 1
   allocate(sbla(nrep))
!
!  Fill the vector x with random values from (-1,1).
!
   call random_number(xs(1:n))
   do i = 1, n
     xs(i) = one - two*xs(i)
   end do
!
!  Compute the expected value of the 1-norm
!
   znrm = zero
   do i = 1, n
     znrm = znrm + abs(xs(i))
   end do
!
!  Do the test nrep times but only test the final one
!
   n = ns
   x(1:n) = xs(1:n)
   incx = incxs
   t1 = SECOND()
   do j = 1, nrep
     sbla(j) = SASUM( n, x, incx )
   end do
   t2 = SECOND()
   tottim = t2 - t1
   rtim = tottim / nrep
!  write(nout,99) 'SASUM ', ns, rtim
!
!  See what values changed inside the subroutine
!
   nerrs = 0
   call IARGCKS( nerrs, subnam, 'N', n, ns, nout )
   call SARGCKV( nerrs, subnam, 'X', ns, x, xs, nout )
   call IARGCKS( nerrs, subnam, 'INCX', incx, incxs, nout )
!
!  Compare sbla and zrnm
!
   if( znrm == zero ) then
     trat = sbla(nrep) / ulp
   else
     trat = ( abs( sbla(nrep) - znrm ) / znrm ) / ( real(max(1,n),wp)*ulp )
   end if
   if (trat >= thresh) then
     write(nout,98) subnam, ns, incxs, trat
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
!
99 format( A6, ': N=', I8, ', time=', E15.8 )
98 format( 1X, A6, ': N=', I8,', INCX=', I4, ', test=', E15.8 )
97 format( 1X, A6, ' passed the computational tests (', I6, ' tests run)', / )
96 format( 1X, A6, ': ', I6, ' of ', I6,' tests failed to pass threshold', / )
end subroutine
