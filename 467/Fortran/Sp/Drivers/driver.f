      program main

c***********************************************************************
c
cc TOMS467_PRB tests XPOSE.
c
      implicit none

      integer a_max
      integer moved_max

      parameter ( a_max = 3000 )
      parameter ( moved_max = 100 )

      real a(a_max)
      logical moved(moved_max)
      integer n1
      integer n12
      integer n2
      integer nwork

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TOMS467_PRB'
      write ( *, '(a)' ) '  Test TOMS algorithm 467, in place'
      write ( *, '(a)' ) '  matrix transposition.'

      n1 = 10
      n2 = 10
      n12 = n1 * n2
      nwork = ( n1 + n2 ) / 2

      write ( *, '(a)' ) ' '
      write ( *, '(a,i6)' ) '  Row dimension N1 =    ', n1
      write ( *, '(a,i6)' ) '  Column dimension N2 = ', n2
      write ( *, '(a,i6)' ) '  Total size N12 =      ', n12 
      write ( *, '(a,i6)' ) '  Workspace NWORK =     ', nwork

      call set_a ( n1, n2, a )

      call print_a ( n1, n2, a, 1, 5, 1, 5 )

      call xpose ( a, n1, n2, n12, moved, nwork )

      call print_a ( n2, n1, a, 1, 5, 1, 5 )

      n1 = 7
      n2 = 30
      n12 = n1 * n2
      nwork = ( n1 + n2 ) / 2

      write ( *, '(a)' ) ' '
      write ( *, '(a,i6)' ) '  Row dimension N1 =    ', n1
      write ( *, '(a,i6)' ) '  Column dimension N2 = ', n2
      write ( *, '(a,i6)' ) '  Total size N12 =      ', n12 
      write ( *, '(a,i6)' ) '  Workspace NWORK =     ', nwork

      call set_a ( n1, n2, a )

      call print_a ( n1, n2, a, 1, 5, 1, 5 )

      call xpose ( a, n1, n2, n12, moved, nwork )

      call print_a ( n2, n1, a, 1, 5, 1, 5 )

      n1 = 24
      n2 = 8
      n12 = n1 * n2
      nwork = ( n1 + n2 ) / 2

      write ( *, '(a)' ) ' '
      write ( *, '(a,i6)' ) '  Row dimension N1 =    ', n1
      write ( *, '(a,i6)' ) '  Column dimension N2 = ', n2
      write ( *, '(a,i6)' ) '  Total size N12 =      ', n12 
      write ( *, '(a,i6)' ) '  Workspace NWORK =     ', nwork

      call set_a ( n1, n2, a )

      call print_a ( n1, n2, a, 1, 5, 1, 5 )

      call xpose ( a, n1, n2, n12, moved, nwork )

      call print_a ( n2, n1, a, 1, 5, 1, 5 )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TOMS467_PRB'
      write ( *, '(a)' ) '  Normal end of execution.'

      stop
      end
      subroutine set_a ( n1, n2, a )

c***********************************************************************
c
cc SET_A sets the matrix A.
c
      implicit none

      integer n1
      integer n2

      real a(n1,n2)
      integer i
      integer j

      do i = 1, n1
        do j = 1, n2
          a(i,j) = 1000 * i + j
        end do
      end do

      return
      end
      subroutine print_a ( m, n, a, i_lo, i_hi, j_lo, j_hi )

c***********************************************************************
c
cc PRINT_A prints the matrix A.
c
      implicit none

      integer m
      integer n

      real a(m,n)
      integer i
      integer i_hi
      integer i_lo
      integer j
      integer j_hi
      integer j_lo

      write ( *, '(a)' ) ' '

      do i = i_lo, i_hi
        write ( *, '(2x,5f8.0)' ) ( a(i,j), j = j_lo, j_hi )
      end do

      return
      end

