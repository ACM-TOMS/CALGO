      program main

c***********************************************************************
c
cc TOMS434_PRB tests CONP.
c
      implicit none

      integer c
      integer r
      parameter ( c = 3 )
!     parameter ( r = 3 )
      parameter ( r = 2 )

      integer nc
      integer nr
      parameter ( nc = c + 1 )
      parameter ( nr = r + 1 )

      integer i
      integer j
      integer matrix(nr,nc)
      real pc
      real ps
      real pt

      save matrix
c
c  Matrix data is by COLUMNS.
c
      data matrix /
     &   1,  1,  2,
     &   2,  0,  2,
     &   3,  2,  5,
     &   6,  3,  9 /

!    &  14, 10,  5, 29,
!    &   2,  4, 10, 16,
!    &   4,  6,  5, 15,
!    &  20, 20, 20, 60 /

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TOMS434_PRB'
      write ( *, '(a)' ) '  Test TOMS algorithm 434, probability'
      write ( *, '(a)' ) '  of RxC contingency tables.'

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  The contingency table:'
      write ( *, '(a)' ) ' '
      do i = 1, nr
        write ( *, '(6(2x,i4))' ) ( matrix(i,j), j = 1, nc )
      end do

      call conp ( matrix, nr, nc, pt, ps, pc )

      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6)' ) '  PT, prob of this table ', pt
      write ( *, '(a,g14.6)' ) 
     &  '  PS, prob of this or less likely table ', ps
      write ( *, '(a,g14.6)' ) '  PC, prob of some table ', pc

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TOMS434_PRB'
      write ( *, '(a)' ) '  Normal end of execution.'

      stop
      end
      function factorial ( n )
      integer i,n
      real factorial
      factorial = 1.0
      do i = 1, n
        factorial = factorial * i
      end do
      return
      end
