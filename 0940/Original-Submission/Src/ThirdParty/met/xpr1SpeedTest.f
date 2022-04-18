c
c
      subroutine calcabc(size,a,b,sa,sb,c)
      implicit none
      integer size,i
      double precision a(size),b(size),sa, sb, c(size)
c
c     For exact measurement of C++ expression template performance,
c     do-loop must not be merged.
c
      do i=1,size
         a(i) = dble(-i+1+size)
      enddo
      do i=1,size  
         b(i) = dble((i-1)*3 + 2) + dexp( dble(i-1)/dble(size) )
      enddo
      do i=1,size
         c(i) = -a(i)*sa + sb*b(i)
c         c(i) = a(i) + b(i)
      enddo
      end


         
