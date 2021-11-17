      FUNCTION REXPU()
!
!     Random-number generator for the exponential distribution
!     Algorithm EA from J. H. Ahrens and U. Dieter,
!     Communications of the ACM, 31 (1988) 1330--1337.
!     Coded by K. G. Hamilton, December 1996, with corrections.
!
      data alog2/0.6931 4718 0559 9453/
      data     a/5.7133 6315 2645 4228/
      data     b/3.4142 1356 2373 0950/
      data    c/-1.6734 0532 4028 4925/
      data     p/0.9802 5814 3468 5472/
      data    aa/5.6005 7075 6973 8080/
      data    bb/3.3468 1064 8056 9850/
      data    hh/0.0026 1067 2360 2095/
      data    dd/0.0857 8643 7626 9050/
!
    1 call random_number(u)
!       Comment out the following line if your RNG can never return 0.0
      if (u.le.0) go to 1               ! Zero-protector
      g = c
    2 u = u+u
      if (u.ge.1.0) go to 4
      g = g + alog2
      go to 2
    4 u = u-1.0
      if (u.gt.p) go to 6
      rexpu = g + aa/(bb-u)
      return
    6 call random_number(u)
      y = a/(b-u)
      call random_number(up)
      if ((up*hh+dd)*(b-u)**2 .gt. exp(-(y+c))) go to 6
      rexpu = g+y
      return
      end
      FUNCTION REXPS()
!
!     Random-number generator for the exponential distribution
!     Algorithm EA from J. H. Ahrens and U. Dieter,
!     Communications of the ACM, 31 (1988) 1330--1337.
!     Coded by K. G. Hamilton, December 1996, with corrections.
!
      data alog2/0.6931 4718 0559 9453/
      data     a/5.7133 6315 2645 4228/
      data     b/3.4142 1356 2373 0950/
      data    c/-1.6734 0532 4028 4925/
      data     p/0.9802 5814 3468 5472/
      data    aa/5.6005 7075 6973 8080/
      data    bb/3.3468 1064 8056 9850/
      data    hh/0.0026 1067 2360 2095/
      data    dd/0.0857 8643 7626 9050/
!
      call random_number(u)
      do while (u.le.0)                 ! Comment out this block 
        call random_number(u)           ! if your RNG can never
      enddo                             ! return exact zero
      g = c
      u = u+u
      do while (u.lt.1.0)
         g = g + alog2
         u = u+u
      enddo
      u = u-1.0
      if (u.le.p) then
        rexps = g + aa/(bb-u)
        return
      endif
      do
        call random_number(u)
        y = a/(b-u)
        call random_number(up)
        if ((up*hh+dd)*(b-u)**2 .le. exp(-(y+c))) then
          rexps = g+y
          return
        endif
      enddo
      end
