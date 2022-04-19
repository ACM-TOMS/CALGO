      program main

c***********************************************************************
c
cc TOMS450_PRB tests ROMIN.
c
      implicit none

      integer n, ncyc
      parameter ( n = 2 )

      real f
      real f2
      external funct
      integer i
      external monitr
      real step
      real x(n)
      real x2(n)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TOMS450_PRB'
      write ( *, '(a)' ) '  Test TOMS algorithm 450, Rosenbrock'
      write ( *, '(a)' ) '  function minimization.'
      write ( *, '(a)' ) ' '
      write ( *, '(a,i6)' ) '  Dimension of problem N = ', n

      x(1) = -1.2E+00
      x(2) =  1.0E+00

      ncyc = 3

      step = 0.1E+00

      call funct ( n, x, f )

      write ( *, '(a)' ) ' '
      write ( *, '(a,2g14.6)' ) '  Initial X = ', ( x(i), i = 1, n )
      write ( *, '(a,g14.6)' ) '  Initial F(X) = ', f

      call romin ( n, x, funct, step, ncyc, monitr )

      call funct ( n, x, f )

      write ( *, '(a)' ) ' '
      write ( *, '(a,2g14.6)' ) '  Final X = ', ( x(i), i = 1, n )
      write ( *, '(a,g14.6)' ) '  Final F(X) = ', f

      x2(1) = 0.99513E+00
      x2(2) = 0.99053E+00

      call funct ( n, x2, f2 )

      write ( *, '(a)' ) ' '
      write ( *, '(a,2g14.6)' ) '  Expected X = ', ( x2(i), i = 1, n )
      write ( *, '(a,g14.6)' ) '  Expected F(X) = ', f2

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TOMS452_PRB'
      write ( *, '(a)' ) '  Normal end of execution.'

      stop
      end
      subroutine monitr ( n, x, f, r, b, con, nr )

c***********************************************************************
c
cc MONITR monitors the progress of the iteration.
c
c  Discussion:
c
c    This version of the monitor routine is designed simply
c    to stop the iteration after 200 steps have been taken.
c
c  Parameters:
c
c    R IS THE ACTUAL NUMBER OF FUNCTION EVALUATIONS (FOR THE
c      INITIAL ESTIMATE R = 0.)
c
c    B IS THE VALUE OF THE EUCLIDEAN NORM OF THE VECTOR
c      REPRESENTING THE TOTAL PROGRESS MADE SINCE THE
c      AXES WERE LAST ROTATED.
c
c    CON IS A LOGICAL VARIABLE.  AT THE START OF THE
c      SUBROUTINE ROMIN CON=.FALSE.  IF THE CONVERGENCE
c      CRITERIA OF THE ROUTINE MONITOR ARE SATISFIED
c      CON MUST BE SET TO .TRUE. TO STOP THE PROCESS.
c
c    NR IS THE MONITOR INDEX.
c
      implicit none

      integer n

      real b
      logical con
      real f
      integer nr
      integer r
      real x(n)

      if ( 200 .le. r ) then
        con = .true.
      end if

      return
      end
      subroutine funct ( n, x, f )

c***********************************************************************
c
cc FUNCT evaluates the scalar function of N variables to be minimized.
c
      implicit none

      integer n

      real f
      real x(n)

      f = 100.0E+00 * ( x(2) - x(1) * x(1) )**2 + ( 1.0E+00 - x(1) )**2

      return
      end
