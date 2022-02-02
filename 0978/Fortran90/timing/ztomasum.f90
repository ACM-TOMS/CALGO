subroutine ZTOMASUM( n, nrep, x, ldx, xs, work, thresh, nout )
   use LA_CONSTANTS, only: wp, zero, one, two, ulp
!
!  .. Scalar Arguments ..
   integer :: ldx, n, nrep, nout
   real(wp) :: thresh
!  ..
!  .. Array Arguments ..
   real(wp) :: work(*)
   complex(wp) :: x(ldx,*), xs(*)
!  ..
!
!  ZTOMASUM times DZASUM from the Level 1 BLAS, a function to compute
!  the 1-norm of a vector x.  This is the out-of-cache test with a
!  different vector used for each trial.
!
!  E. Anderson
!  July 18, 2016
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
   real(wp) :: DZASUM, DSECND
   external :: DZASUM, DSECND
!  ..
!
   nfail = 0
   ntests = 0
   subnam = 'DZASUM'
   ns = n
   incxs = 1
   call random_number(work(1:2*ns))
   do i = 1, ns
     xs(i) = cmplx( one - two*work(2*i-1), one - two*work(2*i), wp )
   end do
   allocate(sbla(nrep))
!
!  Compute the expected value of the 1-norm
!
   znrm = zero
   do i = 1, ns
     znrm = znrm + abs(real(xs(i))) + abs(aimag(xs(i)))
   end do
!
!  Do the test nrep times but only test the final one
!
   tottim = zero
   n = ns
   do j = 1, nrep
     x(1:n,j) = xs(1:n)
   end do
   t1 = DSECND()
   do j = 1, nrep
     n = ns
     incx = incxs
     sbla(j) = DZASUM( n, x(1,j), incx )
   end do
   t2 = DSECND()
   tottim = tottim + t2 - t1
   write(nout,99) 'DZASUM', ns, tottim / nrep
!
!  See what values changed inside the subroutine
!
   nerrs = 0
   call IARGCKS( nerrs, subnam, 'N', n, ns, nout )
   call ZARGCKV( nerrs, subnam, 'X', ns, x(1,nrep), xs, nout )
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
