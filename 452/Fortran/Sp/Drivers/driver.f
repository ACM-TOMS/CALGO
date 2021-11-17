      program main

c***********************************************************************
c
cc TOMS452_PRB tests NXCBN.
c
      implicit none

      integer i
      integer ic(10)
      integer j
      integer m
      integer n
      integer number

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TOMS452_PRB'
      write ( *, '(a)' ) '  Test TOMS algorithm 452, which generates'
      write ( *, '(a)' ) '  combinations of M things from a set of N.'

      n = 10
      m = 3
      number = 5

      do i = 1, n
        if ( i <= m ) then
          ic(i) = 1
        else
          ic(i) = 0
        end if
      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a,i2)' ) '  Let N = ', n
      write ( *, '(a,i2)' ) '  and M = ', m
      write ( *, '(a,i2,a)' ) '  Generate ', number, ' combinations.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Initial combination:'
      write ( *, '(a)' ) ' '
      write ( *, '(2x,10i1)' ) ( ic(i), i = 1, n )
      write ( *, '(a)' ) ' '

      do j = 1, number
        call nxcbn ( n, m, ic )
        write ( *, '(2x,10i1)' ) ( ic(i), i = 1, n )
      end do

      n = 5
      m = 2
      number = 15

      do i = 1, n
        if ( i <= m ) then
          ic(i) = 1
        else
          ic(i) = 0
        end if
      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a,i2)' ) '  Let N = ', n
      write ( *, '(a,i2)' ) '  and M = ', m
      write ( *, '(a,i2,a)' ) '  Generate ', number, ' combinations.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Initial combination:'
      write ( *, '(a)' ) ' '
      write ( *, '(2x,10i1)' ) ( ic(i), i = 1, n )
      write ( *, '(a)' ) ' '

      do j = 1, number
        call nxcbn ( n, m, ic )
        write ( *, '(2x,10i1)' ) ( ic(i), i = 1, n )
      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TOMS452_PRB'
      write ( *, '(a)' ) '  Normal end of execution.'

      stop
      end
