      program main

c***********************************************************************
c
cc TOMS460_PRB tests ADIP.
c
      implicit none

      integer n
      parameter ( n = 10 )

      double precision a
      double precision b
      double precision c
      double precision d
      integer i
      integer ier
      integer iopt
      integer itns
      double precision dmu
      double precision omeh(n)
      double precision omev(n)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TOMS460_PRB'
      write ( *, '(a)' ) '  Test TOMS algorithm 460, optimal'
      write ( *, '(a)' ) '  ADI parameters.'
      write ( *, '(a)' ) ' '

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Given ITNS, compute DMU:'
      write ( *, '(a)' ) ' '

      a = 0.5D+00
      b = 2.0D+00
      c = 0.25D+00
      d = 0.75D+00
      iopt = 1

      do itns = 1, 5

        call adip ( a, b, c, d, iopt, n, itns, dmu, omeh, 
     &    omev, ier )

        write ( *, '(2x,i6,2x,g14.6,2x,i6)' ) itns, dmu, ier

      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Given DMU, compute ITNS:'
      write ( *, '(a)' ) ' '

      a = 0.5D+00
      b = 2.0D+00
      c = 0.25D+00
      d = 0.75D+00
      iopt = 2
      dmu = 1.0D+00

      do i = 1, 8

        dmu = dmu / 10.0D+00

        call adip ( a, b, c, d, iopt, n, itns, dmu, omeh, 
     &    omev, ier )

         write ( *, '(2x,g14.6,2x,i6,2x,i6)' ) dmu, itns, ier

      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TOMS460_PRB'
      write ( *, '(a)' ) '  Normal end of execution.'

      stop
      end

