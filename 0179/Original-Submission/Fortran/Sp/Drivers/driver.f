      program main

c***********************************************************************
c
cc TOMS179_PRB tests TOMS179.
c
c  Modified:
c
c    03 January 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TOMS179_PRB'
      write ( *, '(a)' ) '  Test TOMS algorithm 179, to evaluate'
      write ( *, '(a)' ) '  the modified Beta function.'

      call test01
      call test02

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TOMS179_PRB'
      write ( *, '(a)' ) '  Normal end of execution.'

      stop
      end
      subroutine test01

c***********************************************************************
c
cc TEST01 tests DLGAMA.
c
c  Modified:
c
c    03 January 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision dlgama
      double precision fx
      double precision fx2
      integer n_data
      double precision x

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST01'
      write ( *, '(a)' ) '  Test DLGAMA, which estimates the logarithm'
      write ( *, '(a)' ) '  of the Gamma function.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 
     &  '      X         Exact Value      Computed'
      write ( *, '(a)' ) ' '

      n_data = 0

10    continue

      call gamma_log_values ( n_data, x, fx )

      if ( n_data <= 0 ) then
        go to 20
      end if

      fx2 = dlgama ( x )

      write ( *, '(2x,f8.4,2x,g16.8,2x,g16.8)' ) x, fx, fx2

      go to 10

20    continue

      return
      end
      subroutine test02

c***********************************************************************
c
cc TEST02 tests MDBETA.
c
c  Modified:
c
c    03 January 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision fx
      double precision fx2
      integer ier
      integer n_data
      double precision p
      double precision prob
      double precision q
      double precision x

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST02'
      write ( *, '(a)' ) '  Test MDBETA, which estimates the value of'
      write ( *, '(a)' ) '  the modified Beta function.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 
     &  '      X        P        Q         Exact Value      Computed'
      write ( *, '(a)' ) ' '

      n_data = 0

10    continue

      call beta_cdf_values ( n_data, p, q, x, fx )

      if ( n_data <= 0 ) then
        go to 20
      end if

      call mdbeta ( x, p, q, prob, ier )

      write ( *, '(2x,f8.4,2x,f8.4,2x,f8.4,2x,g16.8,2x,g16.8)' ) 
     &  x, p, q, fx, prob

      go to 10

20    continue

      return
      end
      function dlgama ( x )

c*******************************************************************************
c
cc DLGAMA calculates the natural logarithm of GAMMA ( X ).
c
c  Discussion:
c
c    Computation is based on an algorithm outlined in references 1 and 2.
c    The program uses rational functions that theoretically approximate
c    LOG(GAMMA(X)) to at least 18 significant decimal digits.  The
c    approximation for 12 < X is from Hart et al, while approximations
c    for X < 12.0D+00 are similar to those in Cody and Hillstrom,
c    but are unpublished.
c
c    The accuracy achieved depends on the arithmetic system, the compiler,
c    intrinsic functions, and proper selection of the machine dependent
c    constants.
c
c  Modified:
c
c    03 January 2006
c
c  Author:
c
c    W. J. Cody and L. Stoltz
c    Argonne National Laboratory
c
c  Reference:
c
c    W. J. Cody and Kenneth Hillstrom,
c    Chebyshev Approximations for the Natural Logarithm of the Gamma Function,
c    Mathematics of Computation,
c    Volume 21, 1967, pages 198-203.
c
c    Kenneth Hillstrom,
c    ANL/AMD Program ANLC366S, DGAMMA/DLGAMA,
c    May 1969.
c
c    Hart, Cheney, Lawson, Maehly, Mesztenyi, Rice, Thacher, Witzgall,
c    Computer Approximations,
c    Wiley, 1968.
c
c  Parameters:
c
c    Input, double precision X, the argument of the Gamma function.
c    X must be positive.
c
c    Output, double precision DLGAMA, the logarithm of the Gamma
c    function of X.
c
c  Local Parameters:
c
c    Local, double precision BETA, the radix for the floating-point
c    representation.
c
c    Local, integer MAXEXP, the smallest positive power of BETA that overflows.
c
c    Local, double precision XBIG, the largest argument for which
c    LN(GAMMA(X)) is representable in the machine, the solution to the equation
c      LN(GAMMA(XBIG)) = BETA**MAXEXP.
c
c    Local, double precision FRTBIG, a rough estimate of the fourth root
c    of XBIG.
c
c  Approximate values for some important machines are:
c
c                            BETA      MAXEXP         XBIG     FRTBIG
c
c  CRAY-1        (S.P.)        2        8191       9.62D+2461  3.13D+615
c  Cyber 180/855 (S.P.)        2        1070       1.72D+319   6.44D+79
c  IEEE (IBM/XT) (S.P.)        2         128       4.08D+36    1.42D+9
c  IEEE (IBM/XT) (D.P.)        2        1024       2.55D+305   2.25D+76
c  IBM 3033      (D.P.)       16          63       4.29D+73    2.56D+18
c  VAX D-Format  (D.P.)        2         127       2.05D+36    1.20D+9
c  VAX G-Format  (D.P.)        2        1023       1.28D+305   1.89D+76
c
      implicit none

      double precision c(7)
      double precision corr
      double precision d1
      double precision d2
      double precision d4
      double precision dlgama
      integer i
      double precision frtbig
      double precision p1(8)
      double precision p2(8)
      double precision p4(8)
      double precision pnt68
      double precision q1(8)
      double precision q2(8)
      double precision q4(8)
      double precision res
      double precision sqrtpi
      double precision x
      double precision xbig
      double precision xden
      double precision xeps
      double precision xm1
      double precision xm2
      double precision xm4
      double precision xnum
      double precision xsq

      data c /
     &  -1.910444077728D-03, 
     &   8.4171387781295D-04, 
     &  -5.952379913043012D-04, 
     &   7.93650793500350248D-04, 
     &  -2.777777777777681622553D-03, 
     &   8.333333333333333331554247D-02, 
     &   5.7083835261D-03 /
      data d1 / -5.772156649015328605195174D-01 /
      data d2 /  4.227843350984671393993777D-01 /
      data d4 /  1.791759469228055000094023D+00 /
      data frtbig / 1.42D+09 /
      data p1 /
     &  4.945235359296727046734888D+00, 
     &  2.018112620856775083915565D+02, 
     &  2.290838373831346393026739D+03, 
     &  1.131967205903380828685045D+04, 
     &  2.855724635671635335736389D+04, 
     &  3.848496228443793359990269D+04, 
     &  2.637748787624195437963534D+04, 
     &  7.225813979700288197698961D+03 /
      data p2 /
     &  4.974607845568932035012064D+00, 
     &  5.424138599891070494101986D+02, 
     &  1.550693864978364947665077D+04, 
     &  1.847932904445632425417223D+05, 
     &  1.088204769468828767498470D+06, 
     &  3.338152967987029735917223D+06, 
     &  5.106661678927352456275255D+06, 
     &  3.074109054850539556250927D+06 /
      data p4 /
     &  1.474502166059939948905062D+04, 
     &  2.426813369486704502836312D+06, 
     &  1.214755574045093227939592D+08, 
     &  2.663432449630976949898078D+09, 
     &  2.940378956634553899906876D+10, 
     &  1.702665737765398868392998D+11, 
     &  4.926125793377430887588120D+11, 
     &  5.606251856223951465078242D+11 /
      data pnt68 / 0.6796875D+00 /
      data q1 /
     &  6.748212550303777196073036D+01, 
     &  1.113332393857199323513008D+03, 
     &  7.738757056935398733233834D+03, 
     &  2.763987074403340708898585D+04, 
     &  5.499310206226157329794414D+04, 
     &  6.161122180066002127833352D+04, 
     &  3.635127591501940507276287D+04, 
     &  8.785536302431013170870835D+03 /
      data q2 /
     &  1.830328399370592604055942D+02, 
     &  7.765049321445005871323047D+03, 
     &  1.331903827966074194402448D+05, 
     &  1.136705821321969608938755D+06, 
     &  5.267964117437946917577538D+06, 
     &  1.346701454311101692290052D+07, 
     &  1.782736530353274213975932D+07, 
     &  9.533095591844353613395747D+06 /
      data q4 /
     &  2.690530175870899333379843D+03, 
     &  6.393885654300092398984238D+05, 
     &  4.135599930241388052042842D+07, 
     &  1.120872109616147941376570D+09, 
     &  1.488613728678813811542398D+10, 
     &  1.016803586272438228077304D+11, 
     &  3.417476345507377132798597D+11, 
     &  4.463158187419713286462081D+11 /
      data sqrtpi / 0.9189385332046727417803297D+00 /
      data xbig / 4.08D+36 /
      data xeps / 2.23D-16 /
c
c  Return immediately if the argument is out of range.
c
      if ( x <= 0.0D+00 .or. xbig < x ) then
        dlgama = log ( xbig )
        return
      end if

      if ( x <= xeps ) then

        res = -log ( x )

      else if ( x <= 1.5D+00 ) then

        if ( x < pnt68 ) then
          corr = - log ( x )
          xm1 = x
        else
          corr = 0.0D+00
          xm1 = ( x - 0.5D+00 ) - 0.5D+00
        end if

        if ( x <= 0.5D+00 .or. pnt68 <= x ) then

          xden = 1.0D+00
          xnum = 0.0D+00

          do i = 1, 8
            xnum = xnum * xm1 + p1(i)
            xden = xden * xm1 + q1(i)
          end do

          res = corr + ( xm1 * ( d1 + xm1 * ( xnum / xden ) ) )

        else

          xm2 = ( x - 0.5D+00 ) - 0.5D+00
          xden = 1.0D+00
          xnum = 0.0D+00
          do i = 1, 8
            xnum = xnum * xm2 + p2(i)
            xden = xden * xm2 + q2(i)
          end do

          res = corr + xm2 * ( d2 + xm2 * ( xnum / xden ) )

        end if

      else if ( x <= 4.0D+00 ) then

        xm2 = x - 2.0D+00
        xden = 1.0D+00
        xnum = 0.0D+00
        do i = 1, 8
          xnum = xnum * xm2 + p2(i)
          xden = xden * xm2 + q2(i)
        end do

        res = xm2 * ( d2 + xm2 * ( xnum / xden ) )

      else if ( x <= 12.0D+00 ) then

        xm4 = x - 4.0D+00
        xden = - 1.0D+00
        xnum = 0.0D+00
        do i = 1, 8
          xnum = xnum * xm4 + p4(i)
          xden = xden * xm4 + q4(i)
        end do

        res = d4 + xm4 * ( xnum / xden )

      else

        res = 0.0D+00

        if ( x <= frtbig ) then

          res = c(7)
          xsq = x * x

          do i = 1, 6
            res = res / xsq + c(i)
          end do

        end if

        res = res / x
        corr = log ( x )
        res = res + sqrtpi - 0.5D+00 * corr
        res = res + x * ( corr - 1.0D+00 )

      end if

      dlgama = res

      return
      end
      subroutine gamma_log_values ( n_data, x, fx )

c*******************************************************************************
c
cc GAMMA_LOG_VALUES returns some values of the Log Gamma function.
c
c  Discussion:
c
c    In Mathematica, the function can be evaluated by:
c
c      Log[Gamma[x]]
c
c  Modified:
c
c    03 January 2006
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Milton Abramowitz and Irene Stegun,
c    Handbook of Mathematical Functions,
c    US Department of Commerce, 1964.
c
c    Stephen Wolfram,
c    The Mathematica Book,
c    Fourth Edition,
c    Wolfram Media / Cambridge University Press, 1999.
c
c  Parameters:
c
c    Input/output, integer N_DATA.  The user sets N_DATA to 0 before the
c    first call.  On each call, the routine increments N_DATA by 1, and
c    returns the corresponding data; when there is no more data, the
c    output value of N_DATA will be 0 again.
c
c    Output, double precision X, the argument of the function.
c
c    Output, double precision FX, the value of the function.
c
      implicit none

      integer n_max
      parameter ( n_max = 20 )

      double precision fx
      double precision fx_vec(n_max)
      integer n_data
      double precision x
      double precision x_vec(n_max)

      save fx_vec
      save x_vec

      data fx_vec /
     &  0.1524063822430784D+01,  0.7966778177017837D+00, 
     &  0.3982338580692348D+00,  0.1520596783998375D+00, 
     &  0.0000000000000000D+00, -0.4987244125983972D-01, 
     & -0.8537409000331584D-01, -0.1081748095078604D+00, 
     & -0.1196129141723712D+00, -0.1207822376352452D+00, 
     & -0.1125917656967557D+00, -0.9580769740706586D-01, 
     & -0.7108387291437216D-01, -0.3898427592308333D-01, 
     &  0.00000000000000000D+00,  0.69314718055994530D+00, 
     &  0.17917594692280550D+01,  0.12801827480081469D+02, 
     &  0.39339884187199494D+02,  0.71257038967168009D+02 /
      data x_vec /
     &  0.20D+00,  0.40D+00,  0.60D+00,  0.80D+00, 
     &  1.00D+00,  1.10D+00,  1.20D+00,  1.30D+00, 
     &  1.40D+00,  1.50D+00,  1.60D+00,  1.70D+00, 
     &  1.80D+00,  1.90D+00,  2.00D+00,  3.00D+00, 
     &  4.00D+00, 10.00D+00, 20.00D+00, 30.00D+00 /

      if ( n_data < 0 ) then
        n_data = 0
      end if

      n_data = n_data + 1

      if ( n_max < n_data ) then
        n_data = 0
        x = 0.0D+00
        fx = 0.0D+00
      else
        x = x_vec(n_data)
        fx = fx_vec(n_data)
      end if

      return
      end
      subroutine beta_cdf_values ( n_data, a, b, x, fx )

c*******************************************************************************
c
cc BETA_CDF_VALUES returns some values of the Beta CDF.
c
c  Discussion:
c
c    The incomplete Beta function may be written
c
c      BETA_INC(A,B,X) = Integral (0 to X) T**(A-1) * (1-T)**(B-1) dT
c                      / Integral (0 to 1) T**(A-1) * (1-T)**(B-1) dT
c
c    Thus,
c
c      BETA_INC(A,B,0.0) = 0.0
c      BETA_INC(A,B,1.0) = 1.0
c
c    The incomplete Beta function is also sometimes called the
c    "modified" Beta function, or the "normalized" Beta function
c    or the Beta CDF (cumulative density function.
c
c    In Mathematica, the function can be evaluated by:
c
c      BETA[X,A,B] / BETA[A,B]
c
c    The function can also be evaluated by using the Statistics package:
c
c      Needs["Statistics`ContinuousDistributions`"]
c      dist = BetaDistribution [ a, b ]
c      CDF [ dist, x ]
c
c  Modified:
c
c    04 January 2006
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Milton Abramowitz and Irene Stegun,
c    Handbook of Mathematical Functions,
c    US Department of Commerce, 1964.
c
c    Karl Pearson,
c    Tables of the Incomplete Beta Function,
c    Cambridge University Press, 1968.
c
c    Stephen Wolfram,
c    The Mathematica Book,
c    Fourth Edition,
c    Wolfram Media / Cambridge University Press, 1999.
c
c  Parameters:
c
c    Input/output, integer N_DATA.  The user sets N_DATA to 0 before the
c    first call.  On each call, the routine increments N_DATA by 1, and
c    returns the corresponding data; when there is no more data, the
c    output value of N_DATA will be 0 again.
c
c    Output, double precision A, B, the parameters of the function.
c
c    Output, double precision X, the argument of the function.
c
c    Output, double precision FX, the value of the function.
c
      implicit none

      integer n_max, i
      parameter ( n_max = 42 )

      double precision a
      double precision a_vec(n_max)
      double precision b
      double precision b_vec(n_max)
      double precision fx
      double precision fx_vec(n_max)
      integer n_data
      double precision x
      double precision x_vec(n_max)

      save a_vec
      save b_vec
      save fx_vec
      save x_vec

      data a_vec /   0.5D+00,   0.5D+00,   0.5D+00,   1.0D+00, 
     &   1.0D+00,   1.0D+00,   1.0D+00,   1.0D+00,   2.0D+00, 
     &   2.0D+00,   2.0D+00,   2.0D+00,   2.0D+00,   2.0D+00, 
     &   2.0D+00,   2.0D+00,   2.0D+00,   5.5D+00,  10.0D+00, 
     &  10.0D+00,  10.0D+00,  10.0D+00,  20.0D+00,  20.0D+00, 
     &  20.0D+00,  20.0D+00,  20.0D+00,  30.0D+00,  30.0D+00, 
     &  40.0D+00,  0.1D+01,   0.1D+01,   0.1D+01,   0.1D+01, 
     &   0.1D+01,   0.1D+01,   0.1D+01,   0.1D+01,   0.2D+01, 
     &   0.3D+01,   0.4D+01,   0.5D+01 /
      data b_vec /
     &   0.5D+00,   0.5D+00,   0.5D+00,   0.5D+00,    0.5D+00, 
     &   0.5D+00,   0.5D+00,   1.0D+00,   2.0D+00,   2.0D+00, 
     &   2.0D+00,   2.0D+00,   2.0D+00,   2.0D+00,   2.0D+00, 
     &   2.0D+00,   2.0D+00,   5.0D+00,   0.5D+00,   5.0D+00, 
     &   5.0D+00,  10.0D+00,   5.0D+00,  10.0D+00,  10.0D+00, 
     &  20.0D+00,  20.0D+00,  10.0D+00,  10.0D+00,  20.0D+00, 
     &   0.5D+00,   0.5D+00,   0.5D+00,   0.5D+00,   0.2D+01, 
     &   0.3D+01,   0.4D+01,   0.5D+01,   0.2D+01,   0.2D+01, 
     &   0.2D+01,   0.2D+01 /
      data (fx_vec(i), i=1,20) /
     &  0.6376856085851985D-01,  0.2048327646991335D+00, 
     &  0.1000000000000000D+01,  0.0000000000000000D+00, 
     &  0.5012562893380045D-02,  0.5131670194948620D-01, 
     &  0.2928932188134525D+00,  0.5000000000000000D+00, 
     &  0.2800000000000000D-01,  0.1040000000000000D+00, 
     &  0.2160000000000000D+00,  0.3520000000000000D+00, 
     &  0.5000000000000000D+00,  0.6480000000000000D+00, 
     &  0.7840000000000000D+00,  0.8960000000000000D+00, 
     &  0.9720000000000000D+00,  0.4361908850559777D+00, 
     &  0.1516409096347099D+00,  0.8978271484375000D-01/
      data (fx_vec(i), i=21,42) /
     &  0.1000000000000000D+01,  0.5000000000000000D+00, 
     &  0.4598773297575791D+00,  0.2146816102371739D+00, 
     &  0.9507364826957875D+00,  0.5000000000000000D+00, 
     &  0.8979413687105918D+00,  0.2241297491808366D+00, 
     &  0.7586405487192086D+00,  0.7001783247477069D+00, 
     &  0.5131670194948620D-01,  0.1055728090000841D+00, 
     &  0.1633399734659245D+00,  0.2254033307585166D+00, 
     &  0.3600000000000000D+00,  0.4880000000000000D+00, 
     &  0.5904000000000000D+00,  0.6723200000000000D+00, 
     &  0.2160000000000000D+00,  0.8370000000000000D-01, 
     &  0.3078000000000000D-01,  0.1093500000000000D-01 /
      data x_vec /
     &  0.01D+00,  0.10D+00,  1.00D+00,  0.00D+00,  0.01D+00, 
     &  0.10D+00,  0.50D+00,  0.50D+00,  0.10D+00,  0.20D+00, 
     &  0.30D+00,  0.40D+00,  0.50D+00,  0.60D+00,  0.70D+00, 
     &  0.80D+00,  0.90D+00,  0.50D+00,  0.90D+00,  0.50D+00, 
     &  1.00D+00,  0.50D+00,  0.80D+00,  0.60D+00,  0.80D+00, 
     &  0.50D+00,  0.60D+00,  0.70D+00,  0.80D+00,  0.70D+00, 
     &  0.10D+00,  0.20D+00,  0.30D+00,  0.40D+00,  0.20D+00, 
     &  0.20D+00,  0.20D+00,  0.20D+00,  0.30D+00,  0.30D+00, 
     &  0.30D+00,  0.30D+00 /

      if ( n_data < 0 ) then
        n_data = 0
      end if

      n_data = n_data + 1

      if ( n_max < n_data ) then
        n_data = 0
        a = 0.0D+00
        b = 0.0D+00
        x = 0.0D+00
        fx = 0.0D+00
      else
        a = a_vec(n_data)
        b = b_vec(n_data)
        x = x_vec(n_data)
        fx = fx_vec(n_data)
      end if

      return
      end
