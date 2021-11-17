      program main

c***********************************************************************
c
cc TOMS470_PRB tests FAKUB.
c
c  Modified:
c
c    15 January 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TOMS470_PRB'
      write ( *, '(a)' ) '  Test TOMS algorithm 470, for solving'
      write ( *, '(a)' ) '  an "almost tridiagonal" linear system.'

      call test01
      call test02
      call test03

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TOMS470_PRB'
      write ( *, '(a)' ) '  Normal end of execution.'

      write ( *, '(a)' ) ' '

      stop
      end
      subroutine test01

c***********************************************************************
c
cc TEST01 tests FAKUB.
c
c  Discussion:
c
c    Here, we just give a simple tridiagonal system.
c
c    Our linear system is:
c
c     ( 4  -1   0   0   0 )   (1)    ( 2)
c     (-1   4  -1   0   0 )   (2)    ( 4)
c     ( 0  -1   4  -1   0 ) * (3)  = ( 6)
c     ( 0   0  -1   4  -1 )   (4)    ( 8)
c     ( 0   0   0  -1   4 )   (5)    (16)
c
c  Modified:
c
c    15 January 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer m
      integer mm
      integer n
      integer nn

      parameter ( m = 1 )
      parameter ( n = 5 )
      parameter ( mm = m )
      parameter ( nn = n )

      real alfa
      real b(nn,mm)
      real d(n)
      real eps
      real h(n)
      integer i
      integer jprom(1)
      real s(n)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST01'
      write ( *, '(a)' ) '  Test FAKUB, a linear system solver'
      write ( *, '(a)' ) '  for "almost tridiagonal" systems.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Actually, try a true tridiagonal system.'

      s(1) = 0.0E+00
      do i = 2, n
        s(i) = -1.0E+00
      end do

      do i = 1, n
        d(i) = 4.0E+00
      end do

      do i = 1, n-1
        h(i) = -1.0E+00
      end do
      h(n) = 0.0E+00

      do i = 1, n
        b(i,1) = 2.0E+00 * i
      end do
      b(n,1) = b(n,1) + n + 1

      alfa = 0.25E+00
      eps = 0.000001E+00

      write ( *, '(a)' ) ' '
      write ( *, '(a,i6)' ) '  Number of equations, N = ', n

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Tridiagonal elements:'
      write ( *, '(a)' ) ' '

      do i = 1, n
        write ( *, * ) s(i), d(i), h(i)
      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Right hand side of linear system:'
      write ( *, '(a)' ) ' '

      do i = 1, n
        write ( *, * ) b(i,1)
      end do

      call fakub ( n, s, d, h, b, m, nn, mm, jprom, alfa, eps )

      if ( m .le. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a,i6)' ) '  M = ', m
        write ( *, '(a)' ) ' '
      end if

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Solution, which should be (1,2,3,...,n):'
      write ( *, '(a)' ) ' '

      do i = 1, n
        write ( *, * ) b(i,1)
      end do

      return
      end
      subroutine test02

c***********************************************************************
c
cc TEST02 tests FAKUB.
c
c  Discussion:
c
c    Here we give an "almost tridiagonal" system.
c
c    Our linear system is:
c
c     ( 4  -1   0   0  -1 )   (1)    (-3)
c     (-1   4  -1   0   0 )   (2)    ( 4)
c     ( 0  -1   4  -1   0 ) * (3)  = ( 6)
c     ( 0   0  -1   4  -1 )   (4)    ( 8)
c     (-1   0   0  -1   4 )   (5)    (15)
c
c  Modified:
c
c    15 January 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer m
      integer mm
      integer n
      integer nn

      parameter ( m = 3 )
      parameter ( n = 5 )
      parameter ( mm = m )
      parameter ( nn = n )

      real alfa
      real b(nn,mm)
      real d(n)
      real eps
      real h(n)
      integer i
      integer j
      integer jprom(m-1)
      real s(n)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST02'
      write ( *, '(a)' ) '  Test FAKUB, a linear system solver'
      write ( *, '(a)' ) '  for "almost tridiagonal" systems.'

      s(1) = 0.0E+00
      do i = 2, n
        s(i) = -1.0E+00
      end do

      do i = 1, n
        d(i) = 4.0E+00
      end do

      do i = 1, n-1
        h(i) = -1.0E+00
      end do
      h(n) = 0.0E+00

      do i = 1, n
        b(i,1) = 2.0E+00 * i
      end do
      b(1,1) = b(1,1) - n
      b(n,1) = b(n,1) + n

      do i = 1, n
        do j = 2, m
          b(i,j) = 0.0E+00
        end do
      end do

      b(1,2) = -1.0E+00
      b(n,3) = -1.0E+00

      jprom(1) = n 
      jprom(2) = 1

      alfa = 0.25E+00
      eps = 0.000001E+00

      write ( *, '(a)' ) ' '
      write ( *, '(a,i6)' ) '  Number of equations, N = ', n

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Tridiagonal elements:'
      write ( *, '(a)' ) ' '

      do i = 1, n
        write ( *, * ) s(i), d(i), h(i)
      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Indices of unknowns with nonzero'
      write ( *, '(a)' ) '  non-tridiagonal coefficients.'
      write ( *, '(a)' ) ' '

      do j = 1, m-1
        write ( *, * ) jprom(j)
      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Right hand side of linear system:'
      write ( *, '(a)' ) ' '

      do i = 1, n
        write ( *, * ) b(i,1)
      end do

      call fakub ( n, s, d, h, b, m, nn, mm, jprom, alfa, eps )

      if ( m .le. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a,i6)' ) '  M = ', m
        write ( *, '(a)' ) ' '
      end if

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Solution, which should be (1,2,3,...,n):'
      write ( *, '(a)' ) ' '

      do i = 1, n
        write ( *, * ) b(i,1)
      end do

      return
      end
      subroutine test03

c***********************************************************************
c
cc TEST03 tests GAUSD.
c
      implicit none

      integer n
      integer nn
      
      parameter ( n = 5 )
      parameter ( nn = n )

      real a(nn,nn)
      real b(nn)
      integer i
      integer j
      integer m
      real x(nn)

      save a

      data a /
     &   10.0E+00,  1.0E+00,  2.0E+00,  3.0E+00,  4.0E+00,
     &    1.0E+00, 10.0E+00,  2.0E+00,  3.0E+00,  4.0E+00,
     &    1.0E+00,  2.0E+00, 10.0E+00,  3.0E+00,  4.0E+00,
     &    1.0E+00,  2.0E+00,  3.0E+00, 10.0E+00,  4.0E+00,
     &    1.0E+00,  2.0E+00,  3.0E+00,  4.0E+00, 10.0E+00 /

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST03'
      write ( *, '(a)' ) '  Test GAUSD, a linear system solver.'

      do i = 1, n
        x(i) = i
      end do

      do i = 1, n
        b(i) = 0.0E+00
        do j = 1, n
          b(i) = b(i) + a(i,j) * x(j)
        end do
      end do

      call gausd ( n, a, b, m, nn )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Solution, which should be (1,2,3,...,n):'
      write ( *, '(a)' ) ' '

      do i = 1, n
        write ( *, * ) b(i)
      end do

      return
      end
