      program main

c*******************************************************************************
c
cc TOMS365_PRB tests the CRF routine.
c
c  Modified:
c
c    06 January 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none


      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TOMS365_PRB:'
      write ( *, '(a)' ) '  Test the ACM TOMS 365 Algorithm CRF.'

      call test01
      call test02

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TOMS365_PRB:'
      write ( *, '(a)' ) '  Normal end of execution.'

      write ( *, '(a)' ) ' '

      stop
      end
      subroutine test01

c*******************************************************************************
c
cc TEST01 tests CRF on the function F(Z) = Z + 1.
c
c  Modified:
c
c    08 January 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none

      real de
      real dm
      real ds
      real he
      real hm
      real hs
      complex f01
      external f01
      integer n
      complex ze
      complex zs

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST01:'
      write ( *, '(a)' ) '  CRF uses the downhill method to find' 
      write ( *, '(a)' ) '  a root of a complex analytic function.'
      write ( *, '(a)' ) '  Here, we use F(Z) = Z + 1.'

      zs = ( 2.0E+00, 0.5E+00 )
      hs = 0.25E+00
      hm = 0.0001E+00
      dm = 0.00001E+00

      write ( *, '(a)' ) ' '
      write ( *, '(a,2g14.6)' ) '  Initial estimate ZS =     ', zs
      write ( *, '(a,g14.6)' ) '  Initial stepsize HS =     ', hs
      write ( *, '(a,g14.6)' ) '  Minimum stepsize HM =     ', hm
      write ( *, '(a,g14.6)' ) '  Deviation tolerance DM =  ', dm

      call crf ( zs, hs, hm, dm, f01, ds, ze, he, de, n )

      write ( *, '(a)' ) ' '
      write ( *, '(a,2g14.6)' ) '  Final estimate ZE =       ', ze
      write ( *, '(a,g14.6)' ) '  Final stepsize HE =       ', he
      write ( *, '(a,g14.6)' ) '  Initial deviation DS =    ', ds
      write ( *, '(a,g14.6)' ) '  Final deviation DE =      ', de
      write ( *, '(a,i6)' ) '  Number of iterations, N = ', n

      return
      end
      subroutine test02

c*******************************************************************************
c
cc TEST02 tests CRF on the function F(Z) = Z**5 + 1.
c
c  Modified:
c
c    08 January 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none

      real de
      real dm
      real ds
      real he
      real hm
      real hs
      complex f02
      external f02
      integer n
      complex ze
      complex zs

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST02:'
      write ( *, '(a)' ) '  CRF uses the downhill method to find' 
      write ( *, '(a)' ) '  a root of a complex analytic function.'
      write ( *, '(a)' ) '  Here, we use F(Z) = Z**5 + 1.'

      zs = ( 2.0E+00, 0.5E+00 )
      hs = 0.25E+00
      hm = 0.0001E+00
      dm = 0.00001E+00

      write ( *, '(a)' ) ' '
      write ( *, '(a,2g14.6)' ) '  Initial estimate ZS =     ', zs
      write ( *, '(a,g14.6)' ) '  Initial stepsize HS =     ', hs
      write ( *, '(a,g14.6)' ) '  Minimum stepsize HM =     ', hm
      write ( *, '(a,g14.6)' ) '  Deviation tolerance DM =  ', dm

      call crf ( zs, hs, hm, dm, f02, ds, ze, he, de, n )

      write ( *, '(a)' ) ' '
      write ( *, '(a,2g14.6)' ) '  Final estimate ZE =       ', ze
      write ( *, '(a,g14.6)' ) '  Final stepsize HE =       ', he
      write ( *, '(a,g14.6)' ) '  Initial deviation DS =    ', ds
      write ( *, '(a,g14.6)' ) '  Final deviation DE =      ', de
      write ( *, '(a,i6)' ) '  Number of iterations, N = ', n

      return
      end
      function f01 ( z )

c*******************************************************************************
c
cc F01 evaluates the function F(Z) = Z + 1.
c
c  Modified:
c
c    06 January 2006
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, complex Z, the aergument.
c
c    Output, complex F01, the function value.
c
      implicit none

      complex f01
      complex z


      f01 = z + 1.0E+00

      return
      end
      function f02 ( z )

c*******************************************************************************
c
cc F02 evaluates the function F(Z) = Z**5 + 1.
c
c  Modified:
c
c    06 January 2006
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, complex Z, the aergument.
c
c    Output, complex F02, the function value.
c
      implicit none

      complex f02
      complex z

      f02 = z**5 + 1.0E+00

      return
      end
