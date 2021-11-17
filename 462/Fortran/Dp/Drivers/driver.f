      program main

c***********************************************************************
c
cc TOMS462_PRB tests BIVNOR.
c
      implicit none

      double precision ah
      double precision ak
      double precision bivnor
      double precision cdf
      double precision expect
      double precision r

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TOMS462_PRB'
      write ( *, '(a)' ) '  Test TOMS algorithm 462,'
      write ( *, '(a)' ) '  bivariate normal distribution.'
      write ( *, '(a)' ) ' '

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 
     &  '         H         K         R       BIVNOR          EXPECT'
      write ( *, '(a)' ) ' '

      ah =  0.8D+00
      ak = -1.5D+00
      r =  -0.9D+00
      expect = 0.148D+00

      cdf = bivnor ( ah, ak, r )

      write ( *, '(2x,f9.5,2x,f9.5,2x,f9.5,2x,g14.6,2x,g14.6)' ) 
     &  ah, ak, r, cdf, expect

      ah =  0.6D+00
      ak = -1.4D+00
      r =  -0.7D+00
      expect = 0.208D+00

      cdf = bivnor ( ah, ak, r )

      write ( *, '(2x,f9.5,2x,f9.5,2x,f9.5,2x,g14.6,2x,g14.6)' ) 
     &  ah, ak, r, cdf, expect

      ah =  0.2D+00
      ak = -1.0D+00
      r =  -0.5D+00
      expect = 0.304D+00

      cdf = bivnor ( ah, ak, r )

      write ( *, '(2x,f9.5,2x,f9.5,2x,f9.5,2x,g14.6,2x,g14.6)' ) 
     &  ah, ak, r, cdf, expect

      ah = -1.2D+00
      ak =  0.1D+00
      r =   0.0D+00
      expect = 0.407D+00

      cdf = bivnor ( ah, ak, r )

      write ( *, '(2x,f9.5,2x,f9.5,2x,f9.5,2x,g14.6,2x,g14.6)' ) 
     &  ah, ak, r, cdf, expect

      ah = -1.2D+00
      ak = -0.1D+00
      r =   0.3D+00
      expect = 0.501D+00

      cdf = bivnor ( ah, ak, r )

      write ( *, '(2x,f9.5,2x,f9.5,2x,f9.5,2x,g14.6,2x,g14.6)' ) 
     &  ah, ak, r, cdf, expect

      ah = -0.4D+00
      ak = -0.9D+00
      r =   0.6D+00
      expect = 0.601D+00

      cdf = bivnor ( ah, ak, r )

      write ( *, '(2x,f9.5,2x,f9.5,2x,f9.5,2x,g14.6,2x,g14.6)' ) 
     &  ah, ak, r, cdf, expect
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TOMS462_PRB'
      write ( *, '(a)' ) '  Normal end of execution.'

      stop
      end

