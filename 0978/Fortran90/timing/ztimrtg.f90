subroutine ZTIMRTG( n, nrep, x, y, z, xs, ys, work, thresh, rtim, nout )
   use LA_CONSTANTS, only: wp, zero, half, one, two, ulp, smlnum, bignum
!
!  .. Scalar Arguments ..
   integer :: n, nrep, nout
   real(wp) :: thresh
!  ..
!  .. Array Arguments ..
   real(wp) :: rtim(*), z(*), work(*)
   complex(wp) :: x(*), xs(*), y(*), ys(*)
!  ..
!
!  ZTIMRTG times ZLARGV from LAPACK, a subroutine to compute a vector
!  of plane rotations.  This is the in-cache test with the same vectors
!  used for each trial.
!
!  E. Anderson
!  July 18, 2016
!
!  .. Parameters ..
   integer, parameter :: ntests = 4
!  ..
!  .. Local Scalars ..
   character*6 :: subnam
   integer :: i, incx, incxs, incy, incys, incz, inczs, ix, ixtype, &
      iytype, iztype, j, k, nerrs, nrun, ns
   real(wp) :: rogue, sx, sy, t1, t2, tottim, x1
!  ..
!  .. Local Arrays ..
   integer :: iworst(ntests)
   real(wp) :: result(ntests), worst(ntests)
!  ..
!  .. External Functions ..
   real(wp) :: DSECND
   external :: DSECND
!  ..
!
   nerrs = 0
   nrun = 0
   rogue = -666.666_wp
   subnam = 'ZLARGV'
   ns = n
   incxs = 1
   incys = 1
   inczs = 1
!
!  Set x and y to random values in (-1,-ulp] U [ulp,1)
!
   call random_number(work(1:2*n))
   do i = 1, 2*n
     x1 = work(i)
     if( x1 < half ) then
        work(i) = ulp + (one-ulp)*two*x1
     else
        work(i) = -ulp - (one-ulp)*(two*x1 - one)
     end if
   end do
   do i = 1, n
     xs(i) = cmplx( work(2*i-1), work(2*i), wp )
   end do
   call random_number(work(1:2*n))
   do i = 1, 2*n
     x1 = work(i)
     if( x1 < half ) then
        work(i) = ulp + (one-ulp)*two*x1
     else
        work(i) = -ulp - (one-ulp)*(two*x1 - one)
     end if
   end do
   do i = 1, n
     ys(i) = cmplx( work(2*i-1), work(2*i), wp )
   end do
!
!  Do for each of 9 scaling combinations:
!   1. x small, y small
!   2. x small, y medium
!   3. x small, y large
!   4. x medium, y small
!   5. x medium, y medium
!   6. x medium, y large
!   7. x large, y small
!   8. x large, y medium
!   9. x large, y large
!
   do ixtype = 1, 3
     if( ixtype == 1 ) then
       sx = smlnum
     else if( ixtype == 2 ) then
       sx = one
     else
       sx = bignum
     end if
     do iytype = 1, 3
       if( iytype == 1 ) then
         sy = smlnum
       else if( iytype == 2 ) then
         sy = one
       else
         sy = bignum
       end if
       iztype = (ixtype-1)*3+iytype
!
!      Do the test nrep times but only test the final one
!
       t1 = DSECND()
       do j = 1, nrep
         n = ns
         x(1:n) = sx*xs(1:n)
         incx = incxs
         y(1:n) = sy*ys(1:n)
         incy = incys
         z(1:n) = rogue
         incz = inczs
         call ZLARGV( n, x, incx, y, incy, z, incz )
       end do
       t2 = DSECND()
       tottim = t2 - t1
       rtim(iztype) = tottim / nrep
!
!      See what values changed inside the subroutine
!
       call IARGCKS( nerrs, subnam, 'N', n, ns, nout )
       call IARGCKS( nerrs, subnam, 'INCX', incx, incxs, nout )
       call IARGCKS( nerrs, subnam, 'INCY', incy, incys, nout )
       call IARGCKS( nerrs, subnam, 'INCC', incz, inczs, nout )
!
!      Check results for each i
!
       ix = 1
       worst(1:ntests) = zero
       do i = 1, ns
         call ZRTG01( sx*xs(i), sy*ys(i), z(i), y(i), x(i), result )
!
!        Compare tests to the threshold
!
         do k = 1, ntests
           if( result(k) < thresh ) then
           else
              nerrs = nerrs + 1
           end if
           if( worst(k) < result(k) ) then
!
!            Track the worst result for each test k
!
             iworst(k) = i
             worst(k) = result(k)
           end if
         end do
         nrun = nrun + ntests
       end do
!
!      Print the result of any failing tests
!
       do k = 1, ntests
         if( worst(k) < thresh ) then
         else
           i = iworst(k)
           if( k == 1 ) then
             write(nout,97) subnam, iztype, i, z(i)*xs(i)+y(i)*ys(i), &
                            x(i), k, result(k)
           else if( k == 2 ) then
             write(nout,96) subnam, iztype, i, -conjg(y(i))*xs(i)+z(i)*ys(i), &
                            k, result(k)
           else if( k == 3 .or. k == 4 ) then
             write(nout,98) subnam, iztype, i, z(i), y(i), k, result(k)
           end if
         end if
       end do
!
       if( nerrs /= 0 ) then
         write(nout,95) subnam, nerrs, nrun
       end if
     end do
   end do
   return
!
98 format( A6, ': type=', I1, ', I = ', I7, ', C = ', G13.5, ', S = (', &
           G13.5, ',', G13.5, '), test(', I1, ') =', G13.5 )
97 format( A6, ': type=', I1, ', I = ', I7, ', C*F + S*G = (', G13.5, &
           ',', G13.5, '), R = (', G13.5, ',', G13.5, '), test(', I1, &
           ') =', G13.5 )
96 format( A6, ': type=', I1, ', I = ', I7, ', -S*F + C*G = (', G13.5, &
           ',', G13.5, '), test(', I1, ') =', G13.5 )
95 format( A6, ': ', I8, ' of ', I8,' tests failed', / )
end subroutine
