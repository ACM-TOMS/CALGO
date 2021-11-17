      double precision function dasum(n,dx,incx)
c***begin prologue  dasum
c***date written   791001   (yymmdd)
c***revision date  820801   (yymmdd)
c***category no.  d1a3a
c***keywords  add,blas,double precision,linear algebra,magnitude,sum,
c             vector
c***author  lawson, c. l., (jpl)
c           hanson, r. j., (snla)
c           kincaid, d. r., (u. of texas)
c           krogh, f. t., (jpl)
c***purpose  sum of magnitudes of d.p. vector components
c***description
c
c                b l a s  subprogram
c    description of parameters
c
c     --input--
c        n  number of elements in input vector(s)
c       dx  double precision vector with n elements
c     incx  storage spacing between elements of dx
c
c     --output--
c    dasum  double precision result (zero if n .le. 0)
c
c     returns sum of magnitudes of double precision dx.
c     dasum = sum from 0 to n-1 of dabs(dx(1+i*incx))
c***references  lawson c.l., hanson r.j., kincaid d.r., krogh f.t.,
c                 *basic linear algebra subprograms for fortran usage*,
c                 algorithm no. 539, transactions on mathematical
c                 software, volume 5, number 3, september 1979, 308-323
c***routines called  (none)
c***end prologue  dasum
c
      double precision dx(1)
c***first executable statement  dasum
      dasum = 0.d0
      if(n.le.0)return
      if(incx.eq.1)goto 20
c
c        code for increments not equal to 1.
c
      ns = n*incx
          do 10 i=1,ns,incx
          dasum = dasum + dabs(dx(i))
   10     continue
      return
c
c        code for increments equal to 1.
c
c
c        clean-up loop so remaining vector length is a multiple of 6.
c
   20 m = mod(n,6)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
         dasum = dasum + dabs(dx(i))
   30 continue
      if( n .lt. 6 ) return
   40 mp1 = m + 1
      do 50 i = mp1,n,6
         dasum = dasum + dabs(dx(i)) + dabs(dx(i+1)) + dabs(dx(i+2))
     1   + dabs(dx(i+3)) + dabs(dx(i+4)) + dabs(dx(i+5))
   50 continue
      return
      end
