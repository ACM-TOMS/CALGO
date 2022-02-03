subroutine STIMSCAL( n, nrep, nal, alvals, x, xs, thresh, rtim, nout )
   use LA_CONSTANTS32, only: wp, zero, one, two, ulp
!
!  .. Scalar Arguments ..
   integer :: n, nrep, nal, nout
   real(wp) :: thresh
!  ..
!  .. Array Arguments ..
   real(wp) :: alvals(*), rtim(*), x(*), xs(*)
!  ..
!
!  STIMSCAL times SSCAL from the Level 1 BLAS, a subroutine to scale a
!  vector x by a scalar alpha.  This is the in-cache test with the same
!  vector used for each trial.
!  
!  E. Anderson
!  July 11, 2016
!
!  .. Local Scalars ..
   character*6 :: subnam
   integer :: i, imax, incx, incxs, ix, j, jmax, k, nerrs, nfail, ns, ntests
   real(wp) :: alpha, alphas, ralpha, t1, t2, tottim, tmax, trat, &
               xexp, xmax, xnrm
!  ..
!  .. External Functions ..
   real(wp) :: SECOND
   external :: SECOND
!  ..
!  .. External Subroutines ..
   external :: SSCAL
!  ..
!
   nfail = 0
   ntests = 0
   subnam = 'SSCAL '
   ns = n
   incxs = 1
!
!  Fill the vector x with random values from (-1,1).
!
   call random_number(xs(1:n))
   do i = 1, n
     xs(i) = one - two*xs(i)
   end do
!
!  Test for each value of alpha
!
   do k = 1, nal
     alphas = alvals(k)
!
!    Do the test nrep times but only test the final one
!
     n = ns
     alpha = alphas
     if( alpha == zero ) then
       ralpha = zero
     else
       ralpha = one / alpha
     end if
     x(1:n) = xs(1:n)
     incx = incxs
     if( nrep == 1 ) then
       jmax = 2
     else
       jmax = nrep - mod(nrep,2)
     end if
     t1 = SECOND()
     do j = 1, jmax-1, 2
       call SSCAL( n, alpha, x, incx )
       call SSCAL( n, ralpha, x, incx )
     end do
     t2 = SECOND()
     tottim = t2 - t1
     rtim(k) = tottim / jmax
!
     n = ns
     alpha = alphas
     x(1:n) = xs(1:n)
     incx = incxs
     call SSCAL( n, alpha, x, incx )
!    write(nout,99) 'SSCAL ', ns, alphas, rtim
!
!    See what values changed inside the subroutine
!
     nerrs = 0
     call IARGCKS( nerrs, subnam, 'N', n, ns, nout )
     call SARGCKS( nerrs, subnam, 'ALPHA', alpha, alphas, nout )
     call IARGCKS( nerrs, subnam, 'INCX', incx, incxs, nout )
!
!    Check that the values of x were scaled correctly
!
     ix = 1
     if( incxs < 0 ) ix = 1 - (ns-1)*incxs
     imax = 1
     xmax = zero
     tmax = -one
     do i = 1, ns
       xexp = alphas*xs(i)
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
     if( tmax >= thresh ) then
       write(nout,98) subnam, ns, alphas, incxs, trat
       xexp = alphas*xs(imax)
       write(nout,97) imax, xmax, xexp
       nerrs = nerrs + 1
     end if
     if( nerrs > 0 ) nfail = nfail + 1
     ntests = ntests + 1
   end do
   if( nfail /= 0 ) then
     write(nout,95) subnam, nfail, ntests
   end if
   return
99 format( A6, ': N=', I8, ', ALPHA=', E15.8, ', time=', E15.8 )
98 format( A6, ': N=', I8, ', ALPHA=', E15.8, ', INCX=', I4, &
           ', test=', E15.8 )
97 format( '   occurs at X(', I8, ')=', E15.8, ', expected=', E15.8 )
96 format( A6, ' passed the computational tests (', I6, ' tests run)', / )
95 format( A6, ': ', I6, ' of ', I6,' tests failed to pass threshold', / )
end subroutine
