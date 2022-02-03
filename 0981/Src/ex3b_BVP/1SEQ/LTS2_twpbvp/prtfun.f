      subroutine sprt(n, array)
* sprt prints an array.
      implicit double precision (a-h,o-z)
      dimension array(n)
      
      write(6,900) (array(i), i=1,n)
  900 format(1h ,(7(1pe11.3)))
      return
      end

      subroutine mprt(nrowd, nrow, ncol, array)
* mprt prints a matrix.
      implicit double precision (a-h,o-z)
      dimension array(nrowd, ncol)

      do 400 i = 1, nrow
         write(6,900) i,(array(i,j), j=1,ncol)
  400 continue
  900 format(1h ,i5,(6(1pe11.3)))
      return
      end

