      program main

c***********************************************************************
c
cc TOMS441_PRB tests DIPOLE.
c
      implicit none

      integer sample_num
      integer test_num

      parameter ( sample_num = 1000 )
      parameter ( test_num = 3 )

      real a
      real alpha_test(test_num)
      real alpha
      real b
      real dipole
      integer i
      real mean
      real r
      real r_test(test_num)
      integer seed
      integer test
      real variance
      real x(sample_num)
      real xmax
      real xmin

      save alpha_test
      save r_test

      data alpha_test / 0.0, 0.785398163397448, 1.57079632679490 /
      data r_test / 1.0, 0.5, 0.0 /

      seed = 123456789

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TOMS441_PRB'
      write ( *, '(a)' ) '  Test TOMS algorithm 441, generating'
      write ( *, '(a)' ) '  random deviates from the dipole '
      write ( *, '(a)' ) '  distribution.'

      do test = 1, test_num

        alpha = alpha_test(test)
        r = r_test(test)

        a = r * cos ( alpha )
        b = r * sin ( alpha )

        write ( *, '(a)' ) ' '
        write ( *, '(a,g14.6)' ) ' A = ', a
        write ( *, '(a,g14.6)' ) ' B = ', b

        xmax = -1.0E+30
        xmin = +1.0D+30
        mean = 0.0D+00

        do i = 1, sample_num
          x(i) = dipole ( a, b, seed )
          xmax = max ( xmax, x(i) )
          xmin = min ( xmin, x(i) )
          mean = mean + x(i)
        end do

        mean = mean / real ( sample_num )

        variance = 0.0
        do i = 1, sample_num
          variance = variance + ( x(i) - mean )**2
        end do
        variance = variance / real ( sample_num - 1 )

        write ( *, '(a)' ) ' '
        write ( *, '(a,i6)'    ) '  Sample size =     ', sample_num
        write ( *, '(a,g14.6)' ) '  Sample mean =     ', mean
        write ( *, '(a,g14.6)' ) '  Sample variance = ', variance
        write ( *, '(a,g14.6)' ) '  Sample maximum =  ', xmax
        write ( *, '(a,g14.6)' ) '  Sample minimum =  ', xmin

      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TOMS441_PRB'
      write ( *, '(a)' ) '  Normal end of execution.'

      stop
      end
      function r11 ( seed )

c*******************************************************************************
c
cc R11 returns a pseudorandom number between -1 and +1.
c
c  Modified:
c
c    06 December 2005
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input/output, integer SEED, the "seed" value, which should NOT be 0.
c    On output, SEED has been updated.
c
c    Output, real RD11, a new pseudorandom variate, strictly between -1 and 1.
c
      implicit none

      integer k
      real r11
      integer seed

      k = seed / 127773

      seed = 16807 * ( seed - k * 127773 ) - k * 2836

      if ( seed < 0 ) then
        seed = seed + 2147483647
      end if
c
c  Although SEED can be represented exactly as a 32 bit integer,
c  it generally cannot be represented exactly as a 32 bit real number!
c
      r11 = 2.0 * real ( dble ( seed ) * 4.656612875D-10 ) - 1.0

      return
      end
