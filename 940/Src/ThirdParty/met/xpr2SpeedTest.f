c
c
      subroutine calcabc(size,a,b,x,y,w)
      implicit none
      integer size,size1, r,c,i,k
      double precision a(size+1,size), b(size+1, size),
     $     x(size), y(size), w(size+1)
c
c     For exact measurement of C++ expression template performance,
c     do-loop must not be merged.
c
      size1 = size + 1
      do c=1,size
         do r=1,size1
            a(r,c) = - dble((r-1)*4*7)/dble(size) * 2.0
     $           + 3.0 * dble((c-1)*2*7)/dble(size)
         enddo
      enddo
      do c=1,size
         do r=1,size1
            b(r,c) = 3.0 * dble((r-1)*7)/dble(size)
     $           + dexp( dble((c-1)*3*7)/dble(size) * 0.1 )
         enddo
      enddo
      do i=1,size
         x(i) = 7.0 - dble((i-1)*7)/dble(size)
      enddo
      do i=1,size
         y(i) = dexp( dble((i-1)*7)/dble(size) * 0.3 )
      enddo
      do i=1,size1
         w(i) = 0.0
      enddo
      do k=1,size
         do i=1,size1
            w(i) = w(i) + ( a(i,k) * x(k) - b(i,k) * y(k) ) * 3.0
         enddo
      enddo
      end
