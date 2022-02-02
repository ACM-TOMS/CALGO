subroutine DTOMRTG( n, nrep, x, ldx, y, ldy, z, ldz, xs, ys, &
                    thresh, nout )
   use LA_CONSTANTS, only: wp, zero, half, one, two, ulp, smlnum, bignum
!
!  .. Scalar Arguments ..
   integer :: ldx, ldy, ldz, n, nrep, nout
   real(wp) :: thresh
!  ..
!  .. Array Arguments ..
   real(wp) :: x(ldx,*), xs(*), y(ldy,*), ys(*), z(ldz,*)
!  ..
!
!  DTOMRTG times DLARGV from LAPACK, a subroutine to compute a vector
!  of plane rotations.  This is the out-of-cache test with a different
!  vector used for each trial.
!
!  E. Anderson
!  July 22, 2016
!
!  .. Parameters ..
   integer, parameter :: ntests = 4
!  ..
!  .. Local Scalars ..
   character*6 :: subnam
   integer :: i, incx, incxs, incy, incys, incz, inczs, ixtype, &
      iytype, iztype, j, k, likbla, nerrs, nrun, ns
   real(wp) :: rogue, sx, sy, t1, t2, tottim
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
   rogue = -666.666_wp
   subnam = 'DLARGV'
   ns = n
   incxs = 1
   incys = 1
   inczs = 1
!
!  Set x and y to random values in (-1,ulp] U [ulp,1)
!
   call random_number(xs(1:n))
   do i = 1, n
     if( abs( xs(i) ) < half ) then 
        xs(i) = ulp + (one-ulp)*two*xs(i)
     else
        xs(i) = -ulp - (one-ulp)*(two*xs(i) - one)
     end if
   end do
   call random_number(ys(1:n))
   do i = 1, n
     if( abs( ys(i) ) < half ) then 
        ys(i) = ulp + (one-ulp)*two*ys(i)
     else
        ys(i) = -ulp - (one-ulp)*(two*ys(i) - one)
     end if
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
       x(1:n,1) = smlnum*xs(1:n)
       sx = smlnum
     else if( ixtype == 2 ) then
       x(1:n,1) = xs(1:n)
       sx = one
     else
       x(1:n,1) = bignum*xs(1:n)
       sx = bignum
     end if
     do iytype = 1, 3
       if( iytype == 1 ) then
         y(1:n,1) = smlnum*ys(1:n)
         sy = smlnum
       else if( iytype == 2 ) then
         y(1:n,1) = ys(1:n)
         sy = one
       else
         y(1:n,1) = bignum*ys(1:n)
         sy = bignum
       end if
       iztype = (ixtype-1)*3+iytype
       z(1:n,1) = rogue
       do j = 2, nrep
         x(1:n,j) = x(1:n,1)
         y(1:n,j) = y(1:n,1)
         z(1:n,j) = z(1:n,1)
       end do
!
!      Do the test nrep times but only test the final one
!
       n = ns
       incx = incxs
       incy = incys
       incz = inczs
       t1 = DSECND()
       do j = 1, nrep
         call DLARGV( n, x(1,j), incx, y(1,j), incy, z(1,j), incz )
       end do
       t2 = DSECND()
       tottim = t2 - t1
       write(nout,99) 'DLARGV', ns, iztype, tottim / nrep
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
       nerrs = 0
       nrun = 0
       likbla = 0
       j = nrep
       worst(1:ntests) = zero
       do i = 1, ns
         call DRTG01( likbla, sx*xs(i), sy*ys(i), z(i,j), y(i,j), x(i,j), result )
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
             write(nout,97) subnam, iztype, i, &
                z(i,j)*(sx*xs(i))+y(i,j)*(sy*ys(i)), x(i,j), k, worst(k)
           else if( k == 2 ) then
             write(nout,96) subnam, iztype, i, &
                -y(i,j)*(sx*xs(i))+z(i,j)*(sy*ys(i)), k, worst(k)
           else if( k == 3 .or. k == 4 ) then
             write(nout,98) subnam, iztype, i, z(i,j), y(i,j), k, worst(k)
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
99 format( A6, ': N=', I8, ', type=', I1, ', time=', E15.8 )
98 format( A6, ': type=', I1, ', I = ', I7, ', C = ', G13.5, ', S = ', G13.5, ', test(', I1, &
           ') =', G13.5 )
97 format( A6, ': type=', I1, ', I = ', I7, ', C*F + S*G = ', G13.5, ', R = ', G13.5, &
           ', test(', I1, ') =', G13.5 )
96 format( A6, ': type=', I1, ', I = ', I7, ', -S*F + C*G = ', G13.5, ', test(', I1, ') =', &
           G13.5 )
95 format( A6, ': ', I8, ' of ', I8,' tests failed', / )
end subroutine
