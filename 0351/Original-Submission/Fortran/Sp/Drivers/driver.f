      program main

c***********************************************************************
c
cc TOMS351_PRB tests TOMS351.
c
c  Modified:
c
c    05 January 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TOMS351_PRB'
      write ( *, '(a)' ) '  Test TOMS algorithm 351,'
      write ( *, '(a)' ) '  Rombert integration.'
      write ( *, '(a)' ) ' '

      call test01
      call test02

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TOMS351_PRB'
      write ( *, '(a)' ) '  Normal end of execution.'

      stop
      end
      subroutine test01

c***********************************************************************
c
cc TEST01 tests ROMINT.
c
c  Modified:
c
c    05 January 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none

      real a
      real b
      real eps
      real err
      real f01
      external f01
      integer maxe
      integer maxe_in
      integer n
      real pi
      real val

      eps = 0.000001E+00
      pi = 3.14159265E+00

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST01'
      write ( *, '(a)' ) '  Test ROMINT, for Romberg integration.'
      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6)' ) '  Error tolerance EPS = ', eps
      write ( *, '(a)' ) ' '

      write ( *, '(a)' ) 
     &  '      A         B   MAXE  MAXE     N      ERR             VAL'
      write ( *, '(a)' ) 
     &  '                    in    out'
      write ( *, '(a)' ) ' '

      a = 0.0E+00
      b = pi / 2.0E+00

      do 10 maxe_in = 1, 5

        maxe = maxe_in

        call romint ( val, err, eps, a, b, n, maxe, f01 )

        write ( *, 
     &    '(2x,f8.4,2x,f8.4,2x,i2,4x,i2,2x,i4,2x,g14.6,2x,g14.6)' ) 
     &    a, b, maxe_in, maxe, n, err, val

10    continue

      return
      end
      subroutine test02

c***********************************************************************
c
cc TEST02 tests ROMINT.
c
c  Modified:
c
c    05 January 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none

      real a
      real b
      real eps
      real err
      real f02
      external f02
      integer maxe
      integer maxe_in
      integer n
      real val

      eps = 0.000001E+00

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST02'
      write ( *, '(a)' ) '  Test ROMINT, for Romberg integration.'
      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6)' ) '  Error tolerance EPS = ', eps
      write ( *, '(a)' ) ' '

      write ( *, '(a)' ) 
     &  '      A         B   MAXE  MAXE     N      ERR             VAL'
      write ( *, '(a)' ) 
     &  '                    in    out'
      write ( *, '(a)' ) ' '

      a = 0.0E+00
      b = 4.3E+00

      do 10 maxe_in = 1, 5

        maxe = maxe_in

        call romint ( val, err, eps, a, b, n, maxe, f02 )

        write ( *, 
     &    '(2x,f8.4,2x,f8.4,2x,i2,4x,i2,2x,i4,2x,g14.6,2x,g14.6)' ) 
     &    a, b, maxe_in, maxe, n, err, val

10    continue

      return
      end
      function f01 ( x )

c***********************************************************************
c
cc F01 evaluates test function 1.
c
      implicit none

      real f01
      real x

      f01 = cos ( x )

      return
      end
      function f02 ( x )

c***********************************************************************
c
cc F02 evaluates test function 2.
c
      implicit none

      real f02
      real x

      f02 = exp ( - x * x )

      return
      end
