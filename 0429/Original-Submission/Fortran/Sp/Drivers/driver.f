      program main

c***********************************************************************
c
cc TOMS429_PRB tests POLYAN, ACM TOMS algorithm 429.
c
      implicit none

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TOMS429_PRB'
      write ( *, '(a)' ) '  Test TOMS algorithm 429, which obtains'
      write ( *, '(a)' ) '  information about the roots of '
      write ( *, '(a)' ) '  a polynomial.'

      call test01
      call test02
      call test03

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TOMS429_PRB'
      write ( *, '(a)' ) '  Normal end of execution.'

      stop
      end
      subroutine test01

c***********************************************************************
c
cc TEST01 tests POLYAN with a polynomial suggested by Driessen and Hunt.
c
c  Discussion:
c
c    The original version of POLYAN, when given this polynomial, reported:
c
c    Roots are in annulus of inner radius 0.454E+00 and outer radius 0.836E+01,
c    there are no real positive roots,
c    the negative roots (if any) are between -.454E+00 and -0.836D+01,
c    there are no roots with positive real parts.
c
c    But this is incorrect for the given polynomial.
c
c  Reference:
c
c    H B Driessen and E W Hunt,
c    Remark on Algorithm 429,
c    Communications of the ACM,
c    Volume 16, Number 9, page 579, September 1973.
c
      implicit none

      integer n
      parameter ( n = 4 )

      real c(n)
      real cm(n)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST01'
      write ( *, '(a)' ) '  p(x) = x^4 + 5.6562x^3 + 5.8854x^2'
      write ( *, '(a)' ) '             + 7.3646x   + 6.1354'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Approximate roots:'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  -1.001,'
      write ( *, '(a)' ) '  -4.7741, '
      write ( *, '(a)' ) '   0.0089 + 1.1457 i,'
      write ( *, '(a)' ) '   0.0089 - 1.1457 i.'
      write ( *, '(a)' ) ' '

      c(1) = 5.6562
      c(2) = 5.8854
      c(3) = 7.3646
      c(4) = 6.1354

      call polyan ( c, cm, n )

      return
      end
      subroutine test02

c***********************************************************************
c
cc TEST02 tests POLYAN with a polynomial with a single root.
c
      implicit none

      integer n
      parameter ( n = 4 )

      real c(n)
      real cm(n)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST02'
      write ( *, '(a)' ) '  p(x) = x^4 - 8 x^3 + 24 x^2 - 32 x + 16'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Exact roots:'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  2 (multiplicity 4)'
      write ( *, '(a)' ) ' '

      c(1) =  -8.0
      c(2) =  24.0
      c(3) = -32.0
      c(4) =  16.0

      call polyan ( c, cm, n )

      return
      end
      subroutine test03

c***********************************************************************
c
cc TEST03 tests POLYAN with a polynomial in the roots of unity.
c
      implicit none

      integer n
      parameter ( n = 5 )

      real c(n)
      real cm(n)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST03'
      write ( *, '(a)' ) '  p(x) = x^5 - 1'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Exact roots:'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  0.3090 + 0.9510 i'
      write ( *, '(a)' ) ' -0.8090 + 0.5877 i'
      write ( *, '(a)' ) ' -0.8090 - 0.5877 i'
      write ( *, '(a)' ) '  0.3090 - 0.9510 i'
      write ( *, '(a)' ) '  1'
      write ( *, '(a)' ) ' '

      c(1) =  0.0
      c(2) =  0.0
      c(3) =  0.0
      c(4) =  0.0
      c(5) = -1.0

      call polyan ( c, cm, n )

      return
      end
