subroutine CTIMAMAX( n, nrep, x, xs, work, rtim, nout )
   use LA_CONSTANTS32, only: wp, zero, one, two, bignum
!
!  .. Scalar Arguments ..
   integer :: n, nrep, nout
   real(wp) :: rtim
!  ..
!  .. Array Arguments ..
   real(wp) :: work(*)
   complex(wp) :: x(*), xs(*)
!  ..
!
!  CTIMAMAX times ICAMAX from the Level 1 BLAS, a function to find
!  the index of the largest element in absolute value in a vector x.
!  This is the in-cache test with the same vector used for each trial.
!
!  E. Anderson
!  July 11, 2016
!
!  .. Local Scalars ..
   character*6 :: subnam
   integer :: i, imax, incx, incxs, j, nerrs, nfail, ns, ntests
   integer, allocatable :: ibla(:)
   real(wp) :: t1, t2, tottim
!  ..
!  .. External Functions ..
   integer :: ICAMAX
   real(wp) :: SECOND
   external :: ICAMAX, SECOND
!  ..
!
   nfail = 0
   ntests = 0
   subnam = 'ICAMAX'
   ns = n
   incxs = 1
   allocate(ibla(nrep))
!
!  Fill x with random values and put the largest element in x(n)
!
   call random_number(work(1:2*n))
   do i = 1, n
     xs(i) = cmplx( one - two*work(2*i-1), one - two*work(2*i), wp )
   end do
   imax = n
   xs(imax) = cmplx( bignum, bignum, wp )
!
!  Do the test nrep times but only test the final one
!
   x(1:n) = xs(1:n)
   incx = incxs
   t1 = SECOND()
   do j = 1, nrep
     ibla(j) = ICAMAX( n, x, incx )
   end do
   t2 = SECOND()
   tottim = t2 - t1
   rtim = tottim / nrep
!  write(nout,99) 'ICAMAX', ns, rtim
!
!  Check that it found the right maximum
!
   nerrs = 0
   if( ibla(nrep) /= imax ) then
     write(nout,98) subnam, ns, incxs, ibla(nrep), imax
   end if
!
!  See what values changed inside the subroutine
!
   call IARGCKS( nerrs, subnam, 'N', n, ns, nout )
   call CARGCKV( nerrs, subnam, 'X', ns, x, xs, nout )
   call IARGCKS( nerrs, subnam, 'INCX', incx, incxs, nout )
   if( nerrs > 0 ) nfail = nfail + 1
   ntests = ntests + 1
!
   if( nfail /= 0 ) then
     write(nout,96) subnam, nfail, ntests
   end if
   deallocate(ibla)
   return
99 format( A6, ': N=', I8, ', time=', E15.8 )
98 format( A6, ': N=', I8,', INCX=', I4, ', IMAX=', I6, &
           ', expected=', I6 )
97 format( A6, ' passed the computational tests (', I6, ' tests run)', / )
96 format( A6, ': ', I6, ' of ', I6,' tests failed', / )
end subroutine
