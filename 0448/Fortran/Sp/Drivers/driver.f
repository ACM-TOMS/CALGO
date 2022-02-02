      program main

c***********************************************************************
c
cc TOMS448_PRB tests COUNT.
c
      implicit none

      integer k_max
      integer n_max

      parameter ( k_max = 10 )
      parameter ( n_max = 100 )

      integer c(k_max)
      external count
      integer i
      integer k
      integer n
      integer p(n_max)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TOMS448_PRB'
      write ( *, '(a)' ) '  Test TOMS algorithm 448, number of'
      write ( *, '(a)' ) '  partitions of an integer '
      write ( *, '(a)' ) '  restricted to a set C.'

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'Test 1:'
      write ( *, '(a)' ) '  Let C = 1, 2, 3, ..., N'
      write ( *, '(a)' ) ' '

      n = 10
      k = 10

      do i = 1, k
        c(i) = i
      end do

      call count ( c, k, p, n )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '       I    P(I)'
      write ( *, '(a)' ) ' '
      do i = 1, n
        write ( *, '(2x,i6,2x,i6)' ) i, p(i)
      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'Test 2:'
      write ( *, '(a)' ) '  Let C = 1, 3, 5, 7, ( N + 1 ) / 2'
      write ( *, '(a)' ) ' '

      n = 10
      k = ( n + 1 ) / 2

      do i = 1, k
        c(i) = 2 * i - 1
      end do

      call count ( c, k, p, n )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '       I    P(I)'
      write ( *, '(a)' ) ' '
      do i = 1, n
        write ( *, '(2x,i6,2x,i6)' ) i, p(i)
      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'Test 3:'
      write ( *, '(a)' ) '  Let C = 1, 2, 2'
      write ( *, '(a)' ) ' '

      n = 10
      k = 3

      c(1) = 1
      c(2) = 2
      c(3) = 2

      call count ( c, k, p, n )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '       I    P(I)'
      write ( *, '(a)' ) ' '
      do i = 1, n
        write ( *, '(2x,i6,2x,i6)' ) i, p(i)
      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'Test 4:'
      write ( *, '(a)' ) '  Let C = 1, 5, 10, 25, 50, 100'
      write ( *, '(a)' ) ' '

      n = 100
      k = 6

      c(1) = 1
      c(2) = 5
      c(3) = 10
      c(4) = 25
      c(5) = 50
      c(6) = 100

      call count ( c, k, p, n )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '       I    P(I)'
      write ( *, '(a)' ) ' '
      do i = 1, n
        write ( *, '(2x,i6,2x,i6)' ) i, p(i)
      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TOMS448_PRB'
      write ( *, '(a)' ) '  Normal end of execution.'

      stop
      end

