C      ALGORITHM 713, COLLECTED ALGORITHMS FROM ACM.
C      THIS WORK PUBLISHED IN TRANSACTIONS ON MATHEMATICAL SOFTWARE,
C      VOL. 19, NO. 1, MARCH, 1993, P. 131.

# for Sun
FC = f77
FFLAGS = -O2
LDFLAGS =

# for Convex
#FC = fc
#FFLAGS = -O2 -fi
#LDFLAGS = -fi

# for Cray
#FC = cf77
#FFLAGS = 
#LDFLAGS =

 OBJ =  testall.o  vfnlib.o  fnlib.o util.o
DOBJ = dtestall.o dvfnlib.o dfnlib.o util.o


all : testall dtestall

testall : $(OBJ)
	$(FC) -o $@ $(LDFLAGS) $(OBJ) 

dtestall : $(DOBJ)
	$(FC) -o $@ $(LDFLAGS) $(DOBJ) 

.f :
	$(FC) -c $(FFLAGS) $<

clean :
	rm -f testall dtestall core *.o

    Portable Vectorized Software for Bessel Function Evaluation

           Ronald F. Boisvert and Bonita V. Saunders
         Computing and Applied Mathematics Laboratory
        National Institute of Standards and Technology
                    Gaithersburg, MD 20899

                    boisvert@cam.nist.gov
                    saunders@cam.nist.gov

                        June 19, 1991

  --------
  Contents
  --------

  Makefile ...... A Unix makefile which generates the test programs.

  dfnlib.f ...... FNLIB (double precision)

  dtestall.f .... Test program (double precision)

  dvfnlib.f ..... VFNLIB (double precision)

  fnlib.f ....... FNLIB (single precision)

  util.f ........ Machine-dependent utilities.

  testall.f ..... Test program (single precision)

  vfnlib.f ...... VFNLIB (single precision)


  -----
  Notes
  -----

  The VFNLIB package is contained in the files vfnlib.f, dvfnlib.f,
  and util.f.  The test programs testall.f and dtestall.f compare
  the values computed by VFNLIB with the original FNLIB package.
  The required FNLIB codes are contained in the files fnlib.f and
  dfnlib.f.

  Parts of the package which are machine-dependent are stored in the
  file util.f.  The routines in this file are as follows:

    I1MACH    These return machine-dependent constants according
    R1MACH    to the PORT framework.  They are set up for machines
    D1MACH    using IEEE arithmetic, like the Sun.  To modify them,
              simply comment out the Sun constants and un-comment
              the lines corresponding to your local machine.

    SECOND    Returns the elapsed CP time in seconds since the
              beginning of the run in seconds.  This is set up to
              make standard Unix system calls.  On the Cray, this
              routine is supplied by the system, so it should be
              deleted from this file.
    

  ------------
  Dependencies
  ------------

   testall.f  requires  vfnlib.f  fnlib.f util.f
  dtestall.f  requires dvfnlib.f dfnlib.f util.f

   vfnlib.f   requires util.f
  dvfnlib.f   requires util.f

   fnlib.f    requires util.f
  dfnlib.f    requires util.f
      function besi0 (x)
c april 1977 version.  w. fullerton, c3, los alamos scientific lab.
      dimension bi0cs(12)
c     external alog, besi0e, csevl, exp, inits, r1mach, sqrt
c
c series for bi0        on the interval  0.          to  9.00000d+00
c                                        with weighted error   2.46e-18
c                                         log weighted error  17.61
c                               significant figures required  17.90
c                                    decimal places required  18.15
c
      data bi0 cs( 1) /   -.0766054725 2839144951e0 /
      data bi0 cs( 2) /   1.9273379539 93808270e0 /
      data bi0 cs( 3) /    .2282644586 920301339e0 /
      data bi0 cs( 4) /    .0130489146 6707290428e0 /
      data bi0 cs( 5) /    .0004344270 9008164874e0 /
      data bi0 cs( 6) /    .0000094226 5768600193e0 /
      data bi0 cs( 7) /    .0000001434 0062895106e0 /
      data bi0 cs( 8) /    .0000000016 1384906966e0 /
      data bi0 cs( 9) /    .0000000000 1396650044e0 /
      data bi0 cs(10) /    .0000000000 0009579451e0 /
      data bi0 cs(11) /    .0000000000 0000053339e0 /
      data bi0 cs(12) /    .0000000000 0000000245e0 /
c
      data nti0, xsml, xmax / 0, 0., 0. /
c
      if (nti0.ne.0) go to 10
      nti0 = inits (bi0cs, 12, 0.1*r1mach(3))
      xsml = sqrt (4.0*r1mach(3))
      xmax = alog (r1mach(2))
c
 10   y = abs(x)
      if (y.gt.3.0) go to 20
c
      besi0 = 1.0
      if (y.gt.xsml) besi0 = 2.75 + csevl (y*y/4.5-1.0, bi0cs, nti0)
      return
c
 20   if (y.gt.xmax) call seteru (34hbesi0   abs(x) so big i0 overflows,
     1  34, 1, 2)
c
      besi0 = exp(y) * besi0e(x)
c
      return
      end
      function besi1 (x)
c oct 1983 version.  w. fullerton, c3, los alamos scientific lab.
      dimension bi1cs(11)
c     external alog, besi1e, csevl, exp, inits, r1mach, sqrt
c
c series for bi1        on the interval  0.          to  9.00000d+00
c                                        with weighted error   2.40e-17
c                                         log weighted error  16.62
c                               significant figures required  16.23
c                                    decimal places required  17.14
c
      data bi1 cs( 1) /   -.0019717132 61099859e0 /
      data bi1 cs( 2) /    .4073488766 7546481e0 /
      data bi1 cs( 3) /    .0348389942 99959456e0 /
      data bi1 cs( 4) /    .0015453945 56300123e0 /
      data bi1 cs( 5) /    .0000418885 21098377e0 /
      data bi1 cs( 6) /    .0000007649 02676483e0 /
      data bi1 cs( 7) /    .0000000100 42493924e0 /
      data bi1 cs( 8) /    .0000000000 99322077e0 /
      data bi1 cs( 9) /    .0000000000 00766380e0 /
      data bi1 cs(10) /    .0000000000 00004741e0 /
      data bi1 cs(11) /    .0000000000 00000024e0 /
c
c
      data nti1, xmin, xsml, xmax / 0, 3*0. /
c
      if (nti1.ne.0) go to 10
      nti1 = inits (bi1cs, 11, 0.1*r1mach(3))
      xmin = 2.0*r1mach(1)
      xsml = sqrt (8.0*r1mach(3))
      xmax = alog (r1mach(2))
c
 10   y = abs(x)
      if (y.gt.3.0) go to 20
c
      besi1 = 0.0
      if (y.eq.0.0) return
c
      if (y.le.xmin) call seteru (
     1  37hbesi1   abs(x) so small i1 underflows, 37, 1, 0)
      if (y.gt.xmin) besi1 = 0.5*x
      if (y.gt.xsml) besi1 = x * (.875 + csevl(y*y/4.5-1., bi1cs, nti1))
      return
c
 20   if (y.gt.xmax) call seteru (
     1  34hbesi1   abs(x) so big i1 overflows, 34, 2, 2)
c
      besi1 = exp(y) * besi1e(x)
c
      return
      end
      function besj0 (x)
c april 1977 version.  w. fullerton, c3, los alamos scientific lab.
      dimension bj0cs(13), bm0cs(21), bth0cs(24)
c     external cos, csevl, inits, r1mach, sqrt
c
c series for bj0        on the interval  0.          to  1.60000d+01
c                                        with weighted error   7.47e-18
c                                         log weighted error  17.13
c                               significant figures required  16.98
c                                    decimal places required  17.68
c
      data bj0 cs( 1) /    .1002541619 68939137e0 /
      data bj0 cs( 2) /   -.6652230077 64405132e0 /
      data bj0 cs( 3) /    .2489837034 98281314e0 /
      data bj0 cs( 4) /   -.0332527231 700357697e0 /
      data bj0 cs( 5) /    .0023114179 304694015e0 /
      data bj0 cs( 6) /   -.0000991127 741995080e0 /
      data bj0 cs( 7) /    .0000028916 708643998e0 /
      data bj0 cs( 8) /   -.0000000612 108586630e0 /
      data bj0 cs( 9) /    .0000000009 838650793e0 /
      data bj0 cs(10) /   -.0000000000 124235515e0 /
      data bj0 cs(11) /    .0000000000 001265433e0 /
      data bj0 cs(12) /   -.0000000000 000010619e0 /
      data bj0 cs(13) /    .0000000000 000000074e0 /
c
c series for bm0        on the interval  0.          to  6.25000d-02
c                                        with weighted error   4.98e-17
c                                         log weighted error  16.30
c                               significant figures required  14.97
c                                    decimal places required  16.96
c
      data bm0 cs( 1) /    .0928496163 7381644e0 /
      data bm0 cs( 2) /   -.0014298770 7403484e0 /
      data bm0 cs( 3) /    .0000283057 9271257e0 /
      data bm0 cs( 4) /   -.0000014330 0611424e0 /
      data bm0 cs( 5) /    .0000001202 8628046e0 /
      data bm0 cs( 6) /   -.0000000139 7113013e0 /
      data bm0 cs( 7) /    .0000000020 4076188e0 /
      data bm0 cs( 8) /   -.0000000003 5399669e0 /
      data bm0 cs( 9) /    .0000000000 7024759e0 /
      data bm0 cs(10) /   -.0000000000 1554107e0 /
      data bm0 cs(11) /    .0000000000 0376226e0 /
      data bm0 cs(12) /   -.0000000000 0098282e0 /
      data bm0 cs(13) /    .0000000000 0027408e0 /
      data bm0 cs(14) /   -.0000000000 0008091e0 /
      data bm0 cs(15) /    .0000000000 0002511e0 /
      data bm0 cs(16) /   -.0000000000 0000814e0 /
      data bm0 cs(17) /    .0000000000 0000275e0 /
      data bm0 cs(18) /   -.0000000000 0000096e0 /
      data bm0 cs(19) /    .0000000000 0000034e0 /
      data bm0 cs(20) /   -.0000000000 0000012e0 /
      data bm0 cs(21) /    .0000000000 0000004e0 /
c
c series for bth0       on the interval  0.          to  6.25000d-02
c                                        with weighted error   3.67e-17
c                                         log weighted error  16.44
c                               significant figures required  15.53
c                                    decimal places required  17.13
c
      data bth0cs( 1) /   -.2463916377 4300119e0 /
      data bth0cs( 2) /    .0017370983 07508963e0 /
      data bth0cs( 3) /   -.0000621836 33402968e0 /
      data bth0cs( 4) /    .0000043680 50165742e0 /
      data bth0cs( 5) /   -.0000004560 93019869e0 /
      data bth0cs( 6) /    .0000000621 97400101e0 /
      data bth0cs( 7) /   -.0000000103 00442889e0 /
      data bth0cs( 8) /    .0000000019 79526776e0 /
      data bth0cs( 9) /   -.0000000004 28198396e0 /
      data bth0cs(10) /    .0000000001 02035840e0 /
      data bth0cs(11) /   -.0000000000 26363898e0 /
      data bth0cs(12) /    .0000000000 07297935e0 /
      data bth0cs(13) /   -.0000000000 02144188e0 /
      data bth0cs(14) /    .0000000000 00663693e0 /
      data bth0cs(15) /   -.0000000000 00215126e0 /
      data bth0cs(16) /    .0000000000 00072659e0 /
      data bth0cs(17) /   -.0000000000 00025465e0 /
      data bth0cs(18) /    .0000000000 00009229e0 /
      data bth0cs(19) /   -.0000000000 00003448e0 /
      data bth0cs(20) /    .0000000000 00001325e0 /
      data bth0cs(21) /   -.0000000000 00000522e0 /
      data bth0cs(22) /    .0000000000 00000210e0 /
      data bth0cs(23) /   -.0000000000 00000087e0 /
      data bth0cs(24) /    .0000000000 00000036e0 /
c
      data pi4 / 0.7853981633 9744831e0 /
      data ntj0, ntm0, ntth0, xsml, xmax / 3*0, 2*0./
c
      if (ntj0.ne.0) go to 10
      ntj0 = inits (bj0cs, 13, 0.1*r1mach(3))
      ntm0 = inits (bm0cs, 21, 0.1*r1mach(3))
      ntth0 = inits (bth0cs, 24, 0.1*r1mach(3))
c
      xsml = sqrt (4.0*r1mach(3))
      xmax = 1.0/r1mach(4)
c
 10   y = abs(x)
      if (y.gt.4.0) go to 20
c
      besj0 = 1.0
      if (y.gt.xsml) besj0 = csevl (.125*y*y-1., bj0cs, ntj0)
      return
c
 20   if (y.gt.xmax) call seteru (
     1  42hbesj0   no precision because abs(x) is big, 42, 1, 2)
c
      z = 32.0/y**2 - 1.0
      ampl = (0.75 + csevl (z, bm0cs, ntm0)) / sqrt(y)
      theta = y - pi4 + csevl (z, bth0cs, ntth0) / y
      besj0 = ampl * cos (theta)
c
      return
      end
      function besj1 (x)
c sept 1983 edition.  w. fullerton, c3, los alamos scientific lab.
      dimension bj1cs(12), bm1cs(21), bth1cs(24)
c     external cos, csevl, inits, r1mach, sqrt
c
c series for bj1        on the interval  0.          to  1.60000d+01
c                                        with weighted error   4.48e-17
c                                         log weighted error  16.35
c                               significant figures required  15.77
c                                    decimal places required  16.89
c
      data bj1 cs( 1) /   -.1172614151 3332787e0 /
      data bj1 cs( 2) /   -.2536152183 0790640e0 /
      data bj1 cs( 3) /    .0501270809 84469569e0 /
      data bj1 cs( 4) /   -.0046315148 09625081e0 /
      data bj1 cs( 5) /    .0002479962 29415914e0 /
      data bj1 cs( 6) /   -.0000086789 48686278e0 /
      data bj1 cs( 7) /    .0000002142 93917143e0 /
      data bj1 cs( 8) /   -.0000000039 36093079e0 /
      data bj1 cs( 9) /    .0000000000 55911823e0 /
      data bj1 cs(10) /   -.0000000000 00632761e0 /
      data bj1 cs(11) /    .0000000000 00005840e0 /
      data bj1 cs(12) /   -.0000000000 00000044e0 /
c
c series for bm1        on the interval  0.          to  6.25000d-02
c                                        with weighted error   5.61e-17
c                                         log weighted error  16.25
c                               significant figures required  14.97
c                                    decimal places required  16.91
c
      data bm1 cs( 1) /    .1047362510 931285e0 /
      data bm1 cs( 2) /    .0044244389 3702345e0 /
      data bm1 cs( 3) /   -.0000566163 9504035e0 /
      data bm1 cs( 4) /    .0000023134 9417339e0 /
      data bm1 cs( 5) /   -.0000001737 7182007e0 /
      data bm1 cs( 6) /    .0000000189 3209930e0 /
      data bm1 cs( 7) /   -.0000000026 5416023e0 /
      data bm1 cs( 8) /    .0000000004 4740209e0 /
      data bm1 cs( 9) /   -.0000000000 8691795e0 /
      data bm1 cs(10) /    .0000000000 1891492e0 /
      data bm1 cs(11) /   -.0000000000 0451884e0 /
      data bm1 cs(12) /    .0000000000 0116765e0 /
      data bm1 cs(13) /   -.0000000000 0032265e0 /
      data bm1 cs(14) /    .0000000000 0009450e0 /
      data bm1 cs(15) /   -.0000000000 0002913e0 /
      data bm1 cs(16) /    .0000000000 0000939e0 /
      data bm1 cs(17) /   -.0000000000 0000315e0 /
      data bm1 cs(18) /    .0000000000 0000109e0 /
      data bm1 cs(19) /   -.0000000000 0000039e0 /
      data bm1 cs(20) /    .0000000000 0000014e0 /
      data bm1 cs(21) /   -.0000000000 0000005e0 /
c
c series for bth1       on the interval  0.          to  6.25000d-02
c                                        with weighted error   4.10e-17
c                                         log weighted error  16.39
c                               significant figures required  15.96
c                                    decimal places required  17.08
c
      data bth1cs( 1) /    .7406014102 6313850e0 /
      data bth1cs( 2) /   -.0045717556 59637690e0 /
      data bth1cs( 3) /    .0001198185 10964326e0 /
      data bth1cs( 4) /   -.0000069645 61891648e0 /
      data bth1cs( 5) /    .0000006554 95621447e0 /
      data bth1cs( 6) /   -.0000000840 66228945e0 /
      data bth1cs( 7) /    .0000000133 76886564e0 /
      data bth1cs( 8) /   -.0000000024 99565654e0 /
      data bth1cs( 9) /    .0000000005 29495100e0 /
      data bth1cs(10) /   -.0000000001 24135944e0 /
      data bth1cs(11) /    .0000000000 31656485e0 /
      data bth1cs(12) /   -.0000000000 08668640e0 /
      data bth1cs(13) /    .0000000000 02523758e0 /
      data bth1cs(14) /   -.0000000000 00775085e0 /
      data bth1cs(15) /    .0000000000 00249527e0 /
      data bth1cs(16) /   -.0000000000 00083773e0 /
      data bth1cs(17) /    .0000000000 00029205e0 /
      data bth1cs(18) /   -.0000000000 00010534e0 /
      data bth1cs(19) /    .0000000000 00003919e0 /
      data bth1cs(20) /   -.0000000000 00001500e0 /
      data bth1cs(21) /    .0000000000 00000589e0 /
      data bth1cs(22) /   -.0000000000 00000237e0 /
      data bth1cs(23) /    .0000000000 00000097e0 /
      data bth1cs(24) /   -.0000000000 00000040e0 /
c
c
      data pi4 / 0.7853981633 9744831e0 /
      data ntj1, ntm1, ntth1, xsml, xmin, xmax / 3*0, 3*0./
c
      if (ntj1.ne.0) go to 10
      ntj1 = inits (bj1cs, 12, 0.1*r1mach(3))
      ntm1 = inits (bm1cs, 21, 0.1*r1mach(3))
      ntth1 = inits (bth1cs, 24, 0.1*r1mach(3))
c
      xsml = sqrt (8.0*r1mach(3))
      xmin = 2.0*r1mach(1)
      xmax = 1.0/r1mach(4)
c
 10   y = abs(x)
      if (y.gt.4.0) go to 20
c
      besj1 = 0.0
      if (y.eq.0.0) return
      if (y.le.xmin) call seteru (
     1  37hbesj1   abs(x) so small j1 underflows, 37, 1, 0)
      if (y.gt.xmin) besj1 = 0.5*x
      if (y.gt.xsml) besj1 = x * (.25 + csevl(.125*y*y-1., bj1cs, ntj1))
      return
c
 20   if (y.gt.xmax) call seteru (
     1  42hbesj1   no precision because abs(x) is big, 42, 2, 2)
      z = 32.0/y**2 - 1.0
      ampl = (0.75 + csevl (z, bm1cs, ntm1)) / sqrt(y)
      theta = y - 3.0*pi4 + csevl (z, bth1cs, ntth1) / y
      besj1 = sign (ampl, x) * cos (theta)
c
      return
      end
      function besk0 (x)
c april 1977 version.  w. fullerton, c3, los alamos scientific lab.
      dimension bk0cs(11)
c     external alog, besi0, besk0e, csevl, exp, inits, r1mach, sqrt
c
c series for bk0        on the interval  0.          to  4.00000d+00
c                                        with weighted error   3.57e-19
c                                         log weighted error  18.45
c                               significant figures required  17.99
c                                    decimal places required  18.97
c
      data bk0 cs( 1) /   -.0353273932 3390276872e0 /
      data bk0 cs( 2) /    .3442898999 246284869e0 /
      data bk0 cs( 3) /    .0359799365 1536150163e0 /
      data bk0 cs( 4) /    .0012646154 1144692592e0 /
      data bk0 cs( 5) /    .0000228621 2103119451e0 /
      data bk0 cs( 6) /    .0000002534 7910790261e0 /
      data bk0 cs( 7) /    .0000000019 0451637722e0 /
      data bk0 cs( 8) /    .0000000000 1034969525e0 /
      data bk0 cs( 9) /    .0000000000 0004259816e0 /
      data bk0 cs(10) /    .0000000000 0000013744e0 /
      data bk0 cs(11) /    .0000000000 0000000035e0 /
c
      data ntk0, xsml, xmax / 0, 0., 0. /
c
      if (ntk0.ne.0) go to 10
      ntk0 = inits (bk0cs, 11, 0.1*r1mach(3))
      xsml = sqrt (4.0*r1mach(3))
      xmax = -alog(r1mach(1))
      xmax = xmax - 0.5*xmax*alog(xmax)/(xmax+0.5) - 0.01
c
 10   if (x.le.0.) call seteru (29hbesk0   x is zero or negative, 29,
     1  2, 2)
      if (x.gt.2.) go to 20
c
      y = 0.
      if (x.gt.xsml) y = x*x
      besk0 = -alog(0.5*x)*besi0(x) - .25 + csevl (.5*y-1., bk0cs, ntk0)
      return
c
 20   besk0 = 0.
      if (x.gt.xmax) call seteru (30hbesk0   x so big k0 underflows, 30,
     1  1, 0)
      if (x.gt.xmax) return
c
      besk0 = exp(-x) * besk0e(x)
c
      return
      end
      function besk1 (x)
c july 1980 version.  w. fullerton, c3, los alamos scientific lab.
      dimension bk1cs(11)
c     external alog, besi1, besk1e, csevl, exp, inits, r1mach, sqrt
c
c series for bk1        on the interval  0.          to  4.00000d+00
c                                        with weighted error   7.02e-18
c                                         log weighted error  17.15
c                               significant figures required  16.73
c                                    decimal places required  17.67
c
      data bk1 cs( 1) /    .0253002273 389477705e0 /
      data bk1 cs( 2) /   -.3531559607 76544876e0 /
      data bk1 cs( 3) /   -.1226111808 22657148e0 /
      data bk1 cs( 4) /   -.0069757238 596398643e0 /
      data bk1 cs( 5) /   -.0001730288 957513052e0 /
      data bk1 cs( 6) /   -.0000024334 061415659e0 /
      data bk1 cs( 7) /   -.0000000221 338763073e0 /
      data bk1 cs( 8) /   -.0000000001 411488392e0 /
      data bk1 cs( 9) /   -.0000000000 006666901e0 /
      data bk1 cs(10) /   -.0000000000 000024274e0 /
      data bk1 cs(11) /   -.0000000000 000000070e0 /
c
      data ntk1, xmin, xsml, xmax /0, 0., 0., 0. /
c
      if (ntk1.ne.0) go to 10
      ntk1 = inits (bk1cs, 11, 0.1*r1mach(3))
      xmin = exp (amax1(alog(r1mach(1)), -alog(r1mach(2))) + .01)
      xsml = sqrt (4.0*r1mach(3))
      xmax = -alog(r1mach(1))
      xmax = xmax - 0.5*xmax*alog(xmax)/(xmax+0.5) - 0.01
c
 10   if (x.le.0.) call seteru (29hbesk1   x is zero or negative, 29,
     1  2, 2)
      if (x.gt.2.0) go to 20
c
      if (x.lt.xmin) call seteru (31hbesk1   x so small k1 overflows,
     1  31, 3, 2)
      y = 0.
      if (x.gt.xsml) y = x*x
      besk1 = alog(0.5*x)*besi1(x) +
     1  (0.75 + csevl (.5*y-1., bk1cs, ntk1))/x
      return
c
 20   besk1 = 0.
      if (x.gt.xmax) call seteru (30hbesk1   x so big k1 underflows,
     1  30, 1, 0)
      if (x.gt.xmax) return
c
      besk1 = exp(-x) * besk1e(x)
c
      return
      end
      function besy0 (x)
c august 1980 version.  w. fullerton, c3, los alamos scientific lab.
      dimension by0cs(13), bm0cs(21), bth0cs(24)
c     external alog, besj0, csevl, inits, r1mach, sin, sqrt
c
c series for by0        on the interval  0.          to  1.60000d+01
c                                        with weighted error   1.20e-17
c                                         log weighted error  16.92
c                               significant figures required  16.15
c                                    decimal places required  17.48
c
      data by0 cs( 1) /   -.0112778393 92865573e0 /
      data by0 cs( 2) /   -.1283452375 6042035e0 /
      data by0 cs( 3) /   -.1043788479 9794249e0 /
      data by0 cs( 4) /    .0236627491 83969695e0 /
      data by0 cs( 5) /   -.0020903916 47700486e0 /
      data by0 cs( 6) /    .0001039754 53939057e0 /
      data by0 cs( 7) /   -.0000033697 47162423e0 /
      data by0 cs( 8) /    .0000000772 93842676e0 /
      data by0 cs( 9) /   -.0000000013 24976772e0 /
      data by0 cs(10) /    .0000000000 17648232e0 /
      data by0 cs(11) /   -.0000000000 00188105e0 /
      data by0 cs(12) /    .0000000000 00001641e0 /
      data by0 cs(13) /   -.0000000000 00000011e0 /
c
c series for bm0        on the interval  0.          to  6.25000d-02
c                                        with weighted error   4.98e-17
c                                         log weighted error  16.30
c                               significant figures required  14.97
c                                    decimal places required  16.96
c
      data bm0 cs( 1) /    .0928496163 7381644e0 /
      data bm0 cs( 2) /   -.0014298770 7403484e0 /
      data bm0 cs( 3) /    .0000283057 9271257e0 /
      data bm0 cs( 4) /   -.0000014330 0611424e0 /
      data bm0 cs( 5) /    .0000001202 8628046e0 /
      data bm0 cs( 6) /   -.0000000139 7113013e0 /
      data bm0 cs( 7) /    .0000000020 4076188e0 /
      data bm0 cs( 8) /   -.0000000003 5399669e0 /
      data bm0 cs( 9) /    .0000000000 7024759e0 /
      data bm0 cs(10) /   -.0000000000 1554107e0 /
      data bm0 cs(11) /    .0000000000 0376226e0 /
      data bm0 cs(12) /   -.0000000000 0098282e0 /
      data bm0 cs(13) /    .0000000000 0027408e0 /
      data bm0 cs(14) /   -.0000000000 0008091e0 /
      data bm0 cs(15) /    .0000000000 0002511e0 /
      data bm0 cs(16) /   -.0000000000 0000814e0 /
      data bm0 cs(17) /    .0000000000 0000275e0 /
      data bm0 cs(18) /   -.0000000000 0000096e0 /
      data bm0 cs(19) /    .0000000000 0000034e0 /
      data bm0 cs(20) /   -.0000000000 0000012e0 /
      data bm0 cs(21) /    .0000000000 0000004e0 /
c
c series for bth0       on the interval  0.          to  6.25000d-02
c                                        with weighted error   3.67e-17
c                                         log weighted error  16.44
c                               significant figures required  15.53
c                                    decimal places required  17.13
c
      data bth0cs( 1) /   -.2463916377 4300119e0 /
      data bth0cs( 2) /    .0017370983 07508963e0 /
      data bth0cs( 3) /   -.0000621836 33402968e0 /
      data bth0cs( 4) /    .0000043680 50165742e0 /
      data bth0cs( 5) /   -.0000004560 93019869e0 /
      data bth0cs( 6) /    .0000000621 97400101e0 /
      data bth0cs( 7) /   -.0000000103 00442889e0 /
      data bth0cs( 8) /    .0000000019 79526776e0 /
      data bth0cs( 9) /   -.0000000004 28198396e0 /
      data bth0cs(10) /    .0000000001 02035840e0 /
      data bth0cs(11) /   -.0000000000 26363898e0 /
      data bth0cs(12) /    .0000000000 07297935e0 /
      data bth0cs(13) /   -.0000000000 02144188e0 /
      data bth0cs(14) /    .0000000000 00663693e0 /
      data bth0cs(15) /   -.0000000000 00215126e0 /
      data bth0cs(16) /    .0000000000 00072659e0 /
      data bth0cs(17) /   -.0000000000 00025465e0 /
      data bth0cs(18) /    .0000000000 00009229e0 /
      data bth0cs(19) /   -.0000000000 00003448e0 /
      data bth0cs(20) /    .0000000000 00001325e0 /
      data bth0cs(21) /   -.0000000000 00000522e0 /
      data bth0cs(22) /    .0000000000 00000210e0 /
      data bth0cs(23) /   -.0000000000 00000087e0 /
      data bth0cs(24) /    .0000000000 00000036e0 /
c
c twodpi = 2.0/pi
      data twodpi / 0.6366197723 6758134e0 /
      data pi4 / 0.7853981633 9744831e0 /
      data alnhaf / -0.6931471805 59945309e0 /
      data nty0, ntm0, ntth0, xsml, xmax / 3*0, 2*0./
c
      if (nty0.ne.0) go to 10
      nty0 = inits (by0cs, 13, 0.1*r1mach(3))
      ntm0 = inits (bm0cs, 21, 0.1*r1mach(3))
      ntth0 = inits (bth0cs, 24, 0.1*r1mach(3))
c
      xsml = sqrt (4.0*r1mach(3))
      xmax = 1.0/r1mach(4)
c
 10   if (x.le.0.) call seteru (29hbesy0   x is zero or negative, 29,
     1  1, 2)
      if (x.gt.4.0) go to 20
c
      y = 0.
      if (x.gt.xsml) y = x*x
      besy0 = twodpi*(alnhaf+alog(x))*besj0(x) + .375 +
     1  csevl (.125*y-1., by0cs, nty0)
      return
c
 20   if (x.gt.xmax) call seteru (
     1  37hbesy0   no precision because x is big, 37, 2, 2)
c
      z = 32.0/x**2 - 1.0
      ampl = (0.75 + csevl (z, bm0cs, ntm0)) / sqrt(x)
      theta = x - pi4 + csevl (z, bth0cs, ntth0) / x
      besy0 = ampl * sin (theta)
c
      return
      end
      function besy1 (x)
c april 1977 version.  w. fullerton, c3, los alamos scientific lab.
      dimension by1cs(14), bm1cs(21), bth1cs(24)
c     external alog, besj1, csevl, exp, inits, r1mach, sin, sqrt
c
c series for by1        on the interval  0.          to  1.60000d+01
c                                        with weighted error   1.87e-18
c                                         log weighted error  17.73
c                               significant figures required  17.83
c                                    decimal places required  18.30
c
      data by1 cs( 1) /    .0320804710 0611908629e0 /
      data by1 cs( 2) /   1.2627078974 33500450e0 /
      data by1 cs( 3) /    .0064999618 9992317500e0 /
      data by1 cs( 4) /   -.0893616452 8860504117e0 /
      data by1 cs( 5) /    .0132508812 2175709545e0 /
      data by1 cs( 6) /   -.0008979059 1196483523e0 /
      data by1 cs( 7) /    .0000364736 1487958306e0 /
      data by1 cs( 8) /   -.0000010013 7438166600e0 /
      data by1 cs( 9) /    .0000000199 4539657390e0 /
      data by1 cs(10) /   -.0000000003 0230656018e0 /
      data by1 cs(11) /    .0000000000 0360987815e0 /
      data by1 cs(12) /   -.0000000000 0003487488e0 /
      data by1 cs(13) /    .0000000000 0000027838e0 /
      data by1 cs(14) /   -.0000000000 0000000186e0 /
c
c series for bm1        on the interval  0.          to  6.25000d-02
c                                        with weighted error   5.61e-17
c                                         log weighted error  16.25
c                               significant figures required  14.97
c                                    decimal places required  16.91
c
      data bm1 cs( 1) /    .1047362510 931285e0 /
      data bm1 cs( 2) /    .0044244389 3702345e0 /
      data bm1 cs( 3) /   -.0000566163 9504035e0 /
      data bm1 cs( 4) /    .0000023134 9417339e0 /
      data bm1 cs( 5) /   -.0000001737 7182007e0 /
      data bm1 cs( 6) /    .0000000189 3209930e0 /
      data bm1 cs( 7) /   -.0000000026 5416023e0 /
      data bm1 cs( 8) /    .0000000004 4740209e0 /
      data bm1 cs( 9) /   -.0000000000 8691795e0 /
      data bm1 cs(10) /    .0000000000 1891492e0 /
      data bm1 cs(11) /   -.0000000000 0451884e0 /
      data bm1 cs(12) /    .0000000000 0116765e0 /
      data bm1 cs(13) /   -.0000000000 0032265e0 /
      data bm1 cs(14) /    .0000000000 0009450e0 /
      data bm1 cs(15) /   -.0000000000 0002913e0 /
      data bm1 cs(16) /    .0000000000 0000939e0 /
      data bm1 cs(17) /   -.0000000000 0000315e0 /
      data bm1 cs(18) /    .0000000000 0000109e0 /
      data bm1 cs(19) /   -.0000000000 0000039e0 /
      data bm1 cs(20) /    .0000000000 0000014e0 /
      data bm1 cs(21) /   -.0000000000 0000005e0 /
c
c series for bth1       on the interval  0.          to  6.25000d-02
c                                        with weighted error   4.10e-17
c                                         log weighted error  16.39
c                               significant figures required  15.96
c                                    decimal places required  17.08
c
      data bth1cs( 1) /    .7406014102 6313850e0 /
      data bth1cs( 2) /   -.0045717556 59637690e0 /
      data bth1cs( 3) /    .0001198185 10964326e0 /
      data bth1cs( 4) /   -.0000069645 61891648e0 /
      data bth1cs( 5) /    .0000006554 95621447e0 /
      data bth1cs( 6) /   -.0000000840 66228945e0 /
      data bth1cs( 7) /    .0000000133 76886564e0 /
      data bth1cs( 8) /   -.0000000024 99565654e0 /
      data bth1cs( 9) /    .0000000005 29495100e0 /
      data bth1cs(10) /   -.0000000001 24135944e0 /
      data bth1cs(11) /    .0000000000 31656485e0 /
      data bth1cs(12) /   -.0000000000 08668640e0 /
      data bth1cs(13) /    .0000000000 02523758e0 /
      data bth1cs(14) /   -.0000000000 00775085e0 /
      data bth1cs(15) /    .0000000000 00249527e0 /
      data bth1cs(16) /   -.0000000000 00083773e0 /
      data bth1cs(17) /    .0000000000 00029205e0 /
      data bth1cs(18) /   -.0000000000 00010534e0 /
      data bth1cs(19) /    .0000000000 00003919e0 /
      data bth1cs(20) /   -.0000000000 00001500e0 /
      data bth1cs(21) /    .0000000000 00000589e0 /
      data bth1cs(22) /   -.0000000000 00000237e0 /
      data bth1cs(23) /    .0000000000 00000097e0 /
      data bth1cs(24) /   -.0000000000 00000040e0 /
c
c twodpi = 2.0/pi
      data twodpi / 0.6366197723 6758134e0 /
      data pi4 / 0.7853981633 9744831e0 /
      data nty1, ntm1, ntth1, xmin, xsml, xmax / 3*0, 3*0./
c
      if (nty1.ne.0) go to 10
      nty1 = inits (by1cs, 14, 0.1*r1mach(3))
      ntm1 = inits (bm1cs, 21, 0.1*r1mach(3))
      ntth1 = inits (bth1cs, 24, 0.1*r1mach(3))
c
      xmin = 1.571*exp ( amax1(alog(r1mach(1)), -alog(r1mach(2)))+.01)
      xsml = sqrt (4.0*r1mach(3))
      xmax = 1.0/r1mach(4)
c
 10   if (x.le.0.) call seteru (29hbesy1   x is zero or negative, 29,
     1  1, 2)
      if (x.gt.4.0) go to 20
c
      if (x.lt.xmin) call seteru (31hbesy1   x so small y1 overflows,
     1  31, 3, 2)
      y = 0.
      if (x.gt.xsml) y = x*x
      besy1 = twodpi*alog(0.5*x)*besj1(x) +
     1  (0.5 + csevl (.125*y-1., by1cs, nty1))/x
      return
c
 20   if (x.gt.xmax) call seteru (
     1  37hbesy1   no precision because x is big, 37, 2, 2)
c
      z = 32.0/x**2 - 1.0
      ampl = (0.75 + csevl (z, bm1cs, ntm1)) / sqrt(x)
      theta = x - 3.0*pi4 + csevl (z, bth1cs, ntth1) / x
      besy1 = ampl * sin (theta)
c
      return
      end
      function inits (os, nos, eta)
c april 1977 version.  w. fullerton, c3, los alamos scientific lab.
c
c initialize the orthogonal series so that inits is the number of terms
c needed to insure the error is no larger than eta.  ordinarily, eta
c will be chosen to be one-tenth machine precision.
c
c             input arguments --
c os     array of nos coefficients in an orthogonal series.
c nos    number of coefficients in os.
c eta    requested accuracy of series.
c
      dimension os(nos)
c
      if (nos.lt.1) call seteru (
     1  35hinits   number of coefficients lt 1, 35, 2, 2)
c
      err = 0.
      do 10 ii=1,nos
        i = nos + 1 - ii
        err = err + abs(os(i))
        if (err.gt.eta) go to 20
 10   continue
c
 20   if (i.eq.nos) call seteru (28hinits   eta may be too small, 28,
     1  1, 2)
      inits = i
c
      return
      end
      function csevl (x, cs, n)
c april 1977 version.  w. fullerton, c3, los alamos scientific lab.
c
c evaluate the n-term chebyshev series cs at x.  adapted from
c r. broucke, algorithm 446, c.a.c.m., 16, 254 (1973).  also see fox
c and parker, chebyshev polys in numerical analysis, oxford press, p.56.
c
c             input arguments --
c x      value at which the series is to be evaluated.
c cs     array of n terms of a chebyshev series.  in eval-
c        uating cs, only half the first coef is summed.
c n      number of terms in array cs.
c
      dimension cs(1)
c
      if (n.lt.1) call seteru (28hcsevl   number of terms le 0, 28, 2,2)
      if (n.gt.1000) call seteru (31hcsevl   number of terms gt 1000,
     1  31, 3, 2)
      if (x.lt.(-1.1) .or. x.gt.1.1) call seteru (
     1  25hcsevl   x outside (-1,+1), 25, 1, 1)
c
      b1 = 0.
      b0 = 0.
      twox = 2.*x
      do 10 i=1,n
        b2 = b1
        b1 = b0
        ni = n + 1 - i
        b0 = twox*b1 - b2 + cs(ni)
 10   continue
c
      csevl = 0.5 * (b0-b2)
c
      return
      end
      subroutine seteru (messg, nmessg, nerr, iopt)
      common /cseter/ iunflo
      integer messg(1)
      data iunflo / 0 /
c
      if (iopt.ne.0) call seterr (messg, nmessg, nerr, iopt)
      if (iopt.ne.0) return
c
      if (iunflo.le.0) return
      call seterr (messg, nmessg, nerr, 1)
c
      return
      end
      function besi0e (x)
c july 1977 version.  w. fullerton, c3, los alamos scientific lab.
      dimension bi0cs(12), ai0cs(21), ai02cs(22)
c     external csevl, exp, inits, r1mach, sqrt
c
c series for bi0        on the interval  0.          to  9.00000d+00
c                                        with weighted error   2.46e-18
c                                         log weighted error  17.61
c                               significant figures required  17.90
c                                    decimal places required  18.15
c
      data bi0 cs( 1) /   -.0766054725 2839144951e0 /
      data bi0 cs( 2) /   1.9273379539 93808270e0 /
      data bi0 cs( 3) /    .2282644586 920301339e0 /
      data bi0 cs( 4) /    .0130489146 6707290428e0 /
      data bi0 cs( 5) /    .0004344270 9008164874e0 /
      data bi0 cs( 6) /    .0000094226 5768600193e0 /
      data bi0 cs( 7) /    .0000001434 0062895106e0 /
      data bi0 cs( 8) /    .0000000016 1384906966e0 /
      data bi0 cs( 9) /    .0000000000 1396650044e0 /
      data bi0 cs(10) /    .0000000000 0009579451e0 /
      data bi0 cs(11) /    .0000000000 0000053339e0 /
      data bi0 cs(12) /    .0000000000 0000000245e0 /
c
c series for ai0        on the interval  1.25000d-01 to  3.33333d-01
c                                        with weighted error   7.87e-17
c                                         log weighted error  16.10
c                               significant figures required  14.69
c                                    decimal places required  16.76
c
      data ai0 cs( 1) /    .0757599449 4023796e0 /
      data ai0 cs( 2) /    .0075913808 1082334e0 /
      data ai0 cs( 3) /    .0004153131 3389237e0 /
      data ai0 cs( 4) /    .0000107007 6463439e0 /
      data ai0 cs( 5) /   -.0000079011 7997921e0 /
      data ai0 cs( 6) /   -.0000007826 1435014e0 /
      data ai0 cs( 7) /    .0000002783 8499429e0 /
      data ai0 cs( 8) /    .0000000082 5247260e0 /
      data ai0 cs( 9) /   -.0000000120 4463945e0 /
      data ai0 cs(10) /    .0000000015 5964859e0 /
      data ai0 cs(11) /    .0000000002 2925563e0 /
      data ai0 cs(12) /   -.0000000001 1916228e0 /
      data ai0 cs(13) /    .0000000000 1757854e0 /
      data ai0 cs(14) /    .0000000000 0112822e0 /
      data ai0 cs(15) /   -.0000000000 0114684e0 /
      data ai0 cs(16) /    .0000000000 0027155e0 /
      data ai0 cs(17) /   -.0000000000 0002415e0 /
      data ai0 cs(18) /   -.0000000000 0000608e0 /
      data ai0 cs(19) /    .0000000000 0000314e0 /
      data ai0 cs(20) /   -.0000000000 0000071e0 /
      data ai0 cs(21) /    .0000000000 0000007e0 /
c
c series for ai02       on the interval  0.          to  1.25000d-01
c                                        with weighted error   3.79e-17
c                                         log weighted error  16.42
c                               significant figures required  14.86
c                                    decimal places required  17.09
c
      data ai02cs( 1) /    .0544904110 1410882e0 /
      data ai02cs( 2) /    .0033691164 7825569e0 /
      data ai02cs( 3) /    .0000688975 8346918e0 /
      data ai02cs( 4) /    .0000028913 7052082e0 /
      data ai02cs( 5) /    .0000002048 9185893e0 /
      data ai02cs( 6) /    .0000000226 6668991e0 /
      data ai02cs( 7) /    .0000000033 9623203e0 /
      data ai02cs( 8) /    .0000000004 9406022e0 /
      data ai02cs( 9) /    .0000000000 1188914e0 /
      data ai02cs(10) /   -.0000000000 3149915e0 /
      data ai02cs(11) /   -.0000000000 1321580e0 /
      data ai02cs(12) /   -.0000000000 0179419e0 /
      data ai02cs(13) /    .0000000000 0071801e0 /
      data ai02cs(14) /    .0000000000 0038529e0 /
      data ai02cs(15) /    .0000000000 0001539e0 /
      data ai02cs(16) /   -.0000000000 0004151e0 /
      data ai02cs(17) /   -.0000000000 0000954e0 /
      data ai02cs(18) /    .0000000000 0000382e0 /
      data ai02cs(19) /    .0000000000 0000176e0 /
      data ai02cs(20) /   -.0000000000 0000034e0 /
      data ai02cs(21) /   -.0000000000 0000027e0 /
      data ai02cs(22) /    .0000000000 0000003e0 /
c
      data nti0, ntai0, ntai02, xsml / 3*0, 0. /
c
      if (nti0.ne.0) go to 10
      nti0 = inits (bi0cs, 12, 0.1*r1mach(3))
      ntai0 = inits (ai0cs, 21, 0.1*r1mach(3))
      ntai02 = inits (ai02cs, 22, 0.1*r1mach(3))
      xsml = sqrt (4.0*r1mach(3))
c
 10   y = abs(x)
      if (y.gt.3.0) go to 20
c
      besi0e = 1.0
      if (y.gt.xsml) besi0e = exp(-y) * ( 2.75 +
     1  csevl (y*y/4.5-1.0, bi0cs, nti0) )
      return
c
 20   if (y.le.8.) besi0e = (.375 + csevl ((48./y-11.)/5., ai0cs, ntai0)
     1  ) / sqrt(y)
      if (y.gt.8.) besi0e = (.375 + csevl (16./y-1., ai02cs, ntai02))
     1  / sqrt(y)
c
      return
      end
      function besi1e (x)
c oct 1983 version.  w. fullerton, c3, los alamos scientific lab.
      dimension bi1cs(11), ai1cs(21), ai12cs(22)
c     external csevl, exp, inits, r1mach, sqrt
c
c series for bi1        on the interval  0.          to  9.00000d+00
c                                        with weighted error   2.40e-17
c                                         log weighted error  16.62
c                               significant figures required  16.23
c                                    decimal places required  17.14
c
      data bi1 cs( 1) /   -.0019717132 61099859e0 /
      data bi1 cs( 2) /    .4073488766 7546481e0 /
      data bi1 cs( 3) /    .0348389942 99959456e0 /
      data bi1 cs( 4) /    .0015453945 56300123e0 /
      data bi1 cs( 5) /    .0000418885 21098377e0 /
      data bi1 cs( 6) /    .0000007649 02676483e0 /
      data bi1 cs( 7) /    .0000000100 42493924e0 /
      data bi1 cs( 8) /    .0000000000 99322077e0 /
      data bi1 cs( 9) /    .0000000000 00766380e0 /
      data bi1 cs(10) /    .0000000000 00004741e0 /
      data bi1 cs(11) /    .0000000000 00000024e0 /
c
c series for ai1        on the interval  1.25000d-01 to  3.33333d-01
c                                        with weighted error   6.98e-17
c                                         log weighted error  16.16
c                               significant figures required  14.53
c                                    decimal places required  16.82
c
      data ai1 cs( 1) /   -.0284674418 1881479e0 /
      data ai1 cs( 2) /   -.0192295323 1443221e0 /
      data ai1 cs( 3) /   -.0006115185 8579437e0 /
      data ai1 cs( 4) /   -.0000206997 1253350e0 /
      data ai1 cs( 5) /    .0000085856 1914581e0 /
      data ai1 cs( 6) /    .0000010494 9824671e0 /
      data ai1 cs( 7) /   -.0000002918 3389184e0 /
      data ai1 cs( 8) /   -.0000000155 9378146e0 /
      data ai1 cs( 9) /    .0000000131 8012367e0 /
      data ai1 cs(10) /   -.0000000014 4842341e0 /
      data ai1 cs(11) /   -.0000000002 9085122e0 /
      data ai1 cs(12) /    .0000000001 2663889e0 /
      data ai1 cs(13) /   -.0000000000 1664947e0 /
      data ai1 cs(14) /   -.0000000000 0166665e0 /
      data ai1 cs(15) /    .0000000000 0124260e0 /
      data ai1 cs(16) /   -.0000000000 0027315e0 /
      data ai1 cs(17) /    .0000000000 0002023e0 /
      data ai1 cs(18) /    .0000000000 0000730e0 /
      data ai1 cs(19) /   -.0000000000 0000333e0 /
      data ai1 cs(20) /    .0000000000 0000071e0 /
      data ai1 cs(21) /   -.0000000000 0000006e0 /
c
c series for ai12       on the interval  0.          to  1.25000d-01
c                                        with weighted error   3.55e-17
c                                         log weighted error  16.45
c                               significant figures required  14.69
c                                    decimal places required  17.12
c
      data ai12cs( 1) /    .0285762350 1828014e0 /
      data ai12cs( 2) /   -.0097610974 9136147e0 /
      data ai12cs( 3) /   -.0001105889 3876263e0 /
      data ai12cs( 4) /   -.0000038825 6480887e0 /
      data ai12cs( 5) /   -.0000002512 2362377e0 /
      data ai12cs( 6) /   -.0000000263 1468847e0 /
      data ai12cs( 7) /   -.0000000038 3538039e0 /
      data ai12cs( 8) /   -.0000000005 5897433e0 /
      data ai12cs( 9) /   -.0000000000 1897495e0 /
      data ai12cs(10) /    .0000000000 3252602e0 /
      data ai12cs(11) /    .0000000000 1412580e0 /
      data ai12cs(12) /    .0000000000 0203564e0 /
      data ai12cs(13) /   -.0000000000 0071985e0 /
      data ai12cs(14) /   -.0000000000 0040836e0 /
      data ai12cs(15) /   -.0000000000 0002101e0 /
      data ai12cs(16) /    .0000000000 0004273e0 /
      data ai12cs(17) /    .0000000000 0001041e0 /
      data ai12cs(18) /   -.0000000000 0000382e0 /
      data ai12cs(19) /   -.0000000000 0000186e0 /
      data ai12cs(20) /    .0000000000 0000033e0 /
      data ai12cs(21) /    .0000000000 0000028e0 /
      data ai12cs(22) /   -.0000000000 0000003e0 /
c
      data nti1, ntai1, ntai12, xmin, xsml / 3*0, 2*0. /
c
      if (nti1.ne.0) go to 10
      nti1 = inits (bi1cs, 11, 0.1*r1mach(3))
      ntai1 = inits (ai1cs, 21, 0.1*r1mach(3))
      ntai12 = inits (ai12cs, 22, 0.1*r1mach(3))
c
      xmin = 2.0*r1mach(1)
      xsml = sqrt (8.0*r1mach(3))
c
 10   y = abs(x)
      if (y.gt.3.0) go to 20
c
      besi1e = 0.0
      if (y.eq.0.0) return
c
      if (y.le.xmin) call seteru (
     1  37hbesi1e  abs(x) so small i1 underflows, 37, 1, 0)
      besi1e = 0.0
      if (y.gt.xmin) besi1e = 0.5*x
      if (y.gt.xsml) besi1e = x * (.875 + csevl(y*y/4.5-1., bi1cs,nti1))
      besi1e = exp(-y) * besi1e
      return
c
 20   if (y.le.8.) besi1e = (.375 + csevl ((48./y-11.)/5., ai1cs, ntai1)
     1  ) / sqrt(y)
      if (y.gt.8.) besi1e = (.375 + csevl (16./y-1.0, ai12cs, ntai12))
     1  / sqrt(y)
      besi1e = sign (besi1e, x)
c
      return
      end
      function besk0e (x)
c july 1980 version.  w. fullerton, c3, los alamos scientific lab.
      dimension bk0cs(11), ak0cs(17), ak02cs(14)
c     external alog, besi0, csevl, exp, inits, r1mach, sqrt
c
c series for bk0        on the interval  0.          to  4.00000d+00
c                                        with weighted error   3.57e-19
c                                         log weighted error  18.45
c                               significant figures required  17.99
c                                    decimal places required  18.97
c
      data bk0 cs( 1) /   -.0353273932 3390276872e0 /
      data bk0 cs( 2) /    .3442898999 246284869e0 /
      data bk0 cs( 3) /    .0359799365 1536150163e0 /
      data bk0 cs( 4) /    .0012646154 1144692592e0 /
      data bk0 cs( 5) /    .0000228621 2103119451e0 /
      data bk0 cs( 6) /    .0000002534 7910790261e0 /
      data bk0 cs( 7) /    .0000000019 0451637722e0 /
      data bk0 cs( 8) /    .0000000000 1034969525e0 /
      data bk0 cs( 9) /    .0000000000 0004259816e0 /
      data bk0 cs(10) /    .0000000000 0000013744e0 /
      data bk0 cs(11) /    .0000000000 0000000035e0 /
c
c series for ak0        on the interval  1.25000d-01 to  5.00000d-01
c                                        with weighted error   5.34e-17
c                                         log weighted error  16.27
c                               significant figures required  14.92
c                                    decimal places required  16.89
c
      data ak0 cs( 1) /   -.0764394790 3327941e0 /
      data ak0 cs( 2) /   -.0223565260 5699819e0 /
      data ak0 cs( 3) /    .0007734181 1546938e0 /
      data ak0 cs( 4) /   -.0000428100 6688886e0 /
      data ak0 cs( 5) /    .0000030817 0017386e0 /
      data ak0 cs( 6) /   -.0000002639 3672220e0 /
      data ak0 cs( 7) /    .0000000256 3713036e0 /
      data ak0 cs( 8) /   -.0000000027 4270554e0 /
      data ak0 cs( 9) /    .0000000003 1694296e0 /
      data ak0 cs(10) /   -.0000000000 3902353e0 /
      data ak0 cs(11) /    .0000000000 0506804e0 /
      data ak0 cs(12) /   -.0000000000 0068895e0 /
      data ak0 cs(13) /    .0000000000 0009744e0 /
      data ak0 cs(14) /   -.0000000000 0001427e0 /
      data ak0 cs(15) /    .0000000000 0000215e0 /
      data ak0 cs(16) /   -.0000000000 0000033e0 /
      data ak0 cs(17) /    .0000000000 0000005e0 /
c
c series for ak02       on the interval  0.          to  1.25000d-01
c                                        with weighted error   2.34e-17
c                                         log weighted error  16.63
c                               significant figures required  14.67
c                                    decimal places required  17.20
c
      data ak02cs( 1) /   -.0120186982 6307592e0 /
      data ak02cs( 2) /   -.0091748526 9102569e0 /
      data ak02cs( 3) /    .0001444550 9317750e0 /
      data ak02cs( 4) /   -.0000040136 1417543e0 /
      data ak02cs( 5) /    .0000001567 8318108e0 /
      data ak02cs( 6) /   -.0000000077 7011043e0 /
      data ak02cs( 7) /    .0000000004 6111825e0 /
      data ak02cs( 8) /   -.0000000000 3158592e0 /
      data ak02cs( 9) /    .0000000000 0243501e0 /
      data ak02cs(10) /   -.0000000000 0020743e0 /
      data ak02cs(11) /    .0000000000 0001925e0 /
      data ak02cs(12) /   -.0000000000 0000192e0 /
      data ak02cs(13) /    .0000000000 0000020e0 /
      data ak02cs(14) /   -.0000000000 0000002e0 /
c
      data ntk0, ntak0, ntak02, xsml / 3*0, 0. /
c
      if (ntk0.ne.0) go to 10
      ntk0 = inits (bk0cs, 11, 0.1*r1mach(3))
      ntak0 = inits (ak0cs, 17, 0.1*r1mach(3))
      ntak02 = inits (ak02cs, 14, 0.1*r1mach(3))
      xsml = sqrt (4.0*r1mach(3))
c
 10   if (x.le.0.) call seteru (29hbesk0e  x is zero or negative, 29,
     1  1, 2)
      if (x.gt.2.) go to 20
c
      y = 0.
      if (x.gt.xsml) y = x*x
      besk0e = exp(x) * (-alog(0.5*x)*besi0(x)
     1  - .25 + csevl (.5*y-1., bk0cs, ntk0) )
      return
c
 20   if (x.le.8.) besk0e = (1.25 + csevl ((16./x-5.)/3., ak0cs, ntak0))
     1  / sqrt(x)
      if (x.gt.8.) besk0e = (1.25 + csevl (16./x-1., ak02cs, ntak02))
     1  / sqrt(x)
c
      return
      end
      function besk1e (x)
c july 1980 version.  w. fullerton, c3, los alamos scientific lab.
      dimension bk1cs(11), ak1cs(17), ak12cs(14)
c     external alog, besi1, csevl, exp, inits, r1mach, sqrt
c
c series for bk1        on the interval  0.          to  4.00000d+00
c                                        with weighted error   7.02e-18
c                                         log weighted error  17.15
c                               significant figures required  16.73
c                                    decimal places required  17.67
c
      data bk1 cs( 1) /    .0253002273 389477705e0 /
      data bk1 cs( 2) /   -.3531559607 76544876e0 /
      data bk1 cs( 3) /   -.1226111808 22657148e0 /
      data bk1 cs( 4) /   -.0069757238 596398643e0 /
      data bk1 cs( 5) /   -.0001730288 957513052e0 /
      data bk1 cs( 6) /   -.0000024334 061415659e0 /
      data bk1 cs( 7) /   -.0000000221 338763073e0 /
      data bk1 cs( 8) /   -.0000000001 411488392e0 /
      data bk1 cs( 9) /   -.0000000000 006666901e0 /
      data bk1 cs(10) /   -.0000000000 000024274e0 /
      data bk1 cs(11) /   -.0000000000 000000070e0 /
c
c series for ak1        on the interval  1.25000d-01 to  5.00000d-01
c                                        with weighted error   6.06e-17
c                                         log weighted error  16.22
c                               significant figures required  15.41
c                                    decimal places required  16.83
c
      data ak1 cs( 1) /    .2744313406 973883e0 /
      data ak1 cs( 2) /    .0757198995 3199368e0 /
      data ak1 cs( 3) /   -.0014410515 5647540e0 /
      data ak1 cs( 4) /    .0000665011 6955125e0 /
      data ak1 cs( 5) /   -.0000043699 8470952e0 /
      data ak1 cs( 6) /    .0000003540 2774997e0 /
      data ak1 cs( 7) /   -.0000000331 1163779e0 /
      data ak1 cs( 8) /    .0000000034 4597758e0 /
      data ak1 cs( 9) /   -.0000000003 8989323e0 /
      data ak1 cs(10) /    .0000000000 4720819e0 /
      data ak1 cs(11) /   -.0000000000 0604783e0 /
      data ak1 cs(12) /    .0000000000 0081284e0 /
      data ak1 cs(13) /   -.0000000000 0011386e0 /
      data ak1 cs(14) /    .0000000000 0001654e0 /
      data ak1 cs(15) /   -.0000000000 0000248e0 /
      data ak1 cs(16) /    .0000000000 0000038e0 /
      data ak1 cs(17) /   -.0000000000 0000006e0 /
c
c series for ak12       on the interval  0.          to  1.25000d-01
c                                        with weighted error   2.58e-17
c                                         log weighted error  16.59
c                               significant figures required  15.22
c                                    decimal places required  17.16
c
      data ak12cs( 1) /    .0637930834 3739001e0 /
      data ak12cs( 2) /    .0283288781 3049721e0 /
      data ak12cs( 3) /   -.0002475370 6739052e0 /
      data ak12cs( 4) /    .0000057719 7245160e0 /
      data ak12cs( 5) /   -.0000002068 9392195e0 /
      data ak12cs( 6) /    .0000000097 3998344e0 /
      data ak12cs( 7) /   -.0000000005 5853361e0 /
      data ak12cs( 8) /    .0000000000 3732996e0 /
      data ak12cs( 9) /   -.0000000000 0282505e0 /
      data ak12cs(10) /    .0000000000 0023720e0 /
      data ak12cs(11) /   -.0000000000 0002176e0 /
      data ak12cs(12) /    .0000000000 0000215e0 /
      data ak12cs(13) /   -.0000000000 0000022e0 /
      data ak12cs(14) /    .0000000000 0000002e0 /
c
      data ntk1, ntak1, ntak12, xmin, xsml / 3*0, 2*0. /
c
      if (ntk1.ne.0) go to 10
      ntk1 = inits (bk1cs, 11, 0.1*r1mach(3))
      ntak1 = inits (ak1cs, 17, 0.1*r1mach(3))
      ntak12 = inits (ak12cs, 14, 0.1*r1mach(3))
c
      xmin = exp (amax1(alog(r1mach(1)), -alog(r1mach(2))) + .01)
      xsml = sqrt (4.0*r1mach(3))
c
 10   if (x.le.0.) call seteru (29hbesk1e  x is zero or negative, 29,
     1  1, 2)
      if (x.gt.2.0) go to 20
c
      if (x.lt.xmin) call seteru (31hbesk1e  x so small k1 overflows,
     1  31, 2, 2)
      y = 0.
      if (x.gt.xsml) y = x*x
      besk1e = exp(x) * (alog(0.5*x)*besi1(x) +
     1  (0.75 + csevl (.5*y-1., bk1cs, ntk1))/x )
      return
c
 20   if (x.le.8.) besk1e = (1.25 + csevl ((16./x-5.)/3., ak1cs, ntak1))
     1  / sqrt(x)
      if (x.gt.8.) besk1e = (1.25 + csevl (16./x-1., ak12cs, ntak12))
     1  / sqrt(x)
c
      return
      end
      integer function i8save(isw,ivalue,set)
c
c  if (isw = 1) i8save returns the current error number and
c               sets it to ivalue if set = .true. .
c
c  if (isw = 2) i8save returns the current recovery switch and
c               sets it to ivalue if set = .true. .
c
      logical set
c
      integer iparam(2)
c  iparam(1) is the error number and iparam(2) is the recovery switch.
c
c  start execution error free and with recovery turned off.
c
      data iparam(1) /0/,  iparam(2) /2/
c
      i8save=iparam(isw)
      if (set) iparam(isw)=ivalue
c
      return
c
      end
      subroutine eprint
c
c  this subroutine prints the last error message, if any.
c
      integer messg(1)
c
      call e9rint(messg,1,1,.false.)
      return
c
      end
      subroutine s88fmt( n, w, ifmt )
c
c  s88fmt  replaces ifmt(1), ... , ifmt(n) with
c  the characters corresponding to the n least significant
c  digits of w.
c
      integer n,w,ifmt(n)
c
      integer nt,wt,digits(10)
c
      data digits( 1) / 1h0 /
      data digits( 2) / 1h1 /
      data digits( 3) / 1h2 /
      data digits( 4) / 1h3 /
      data digits( 5) / 1h4 /
      data digits( 6) / 1h5 /
      data digits( 7) / 1h6 /
      data digits( 8) / 1h7 /
      data digits( 9) / 1h8 /
      data digits(10) / 1h9 /
c
      nt = n
      wt = w
c
 10   if (nt .le. 0) return
        idigit = mod( wt, 10 )
        ifmt(nt) = digits(idigit+1)
        wt = wt/10
        nt = nt - 1
        go to 10
c
      end
      subroutine seterr (messg, nmessg, nerr, iopt)
c
c  this version modified by w. fullerton to dump if iopt = 1 and
c  not recovering.
c  seterr sets lerror = nerr, optionally prints the message and dumps
c  according to the following rules...
c
c    if iopt = 1 and recovering      - just remember the error.
c    if iopt = 1 and not recovering  - print, dump and stop.
c    if iopt = 2                     - print, dump and stop.
c
c  input
c
c    messg  - the error message.
c    nmessg - the length of the message, in characters.
c    nerr   - the error number. must have nerr non-zero.
c    iopt   - the option. must have iopt=1 or 2.
c
c  error states -
c
c    1 - message length not positive.
c    2 - cannot have nerr=0.
c    3 - an unrecovered error followed by another error.
c    4 - bad value for iopt.
c
c  only the first 72 characters of the message are printed.
c
c  the error handler calls a subroutine named fdump to produce a
c  symbolic dump. to complete the package, a dummy version of fdump
c  is supplied, but it should be replaced by a locally written version
c  which at least gives a trace-back.
c
      integer messg(1)
      external i1mach, i8save
c
c  the unit for error messages.
c
      iwunit=i1mach(4)
c
      if (nmessg.ge.1) go to 10
c
c  a message of non-positive length is fatal.
c
        write(iwunit,9000)
 9000   format(52h1error    1 in seterr - message length not positive.)
        go to 60
c
c  nw is the number of words the message occupies.
c
 10   nw=(min0(nmessg,72)-1)/i1mach(6)+1
c
      if (nerr.ne.0) go to 20
c
c  cannot turn the error state off using seterr.
c
        write(iwunit,9001)
 9001   format(42h1error    2 in seterr - cannot have nerr=0//
     1         34h the current error message follows///)
        call e9rint(messg,nw,nerr,.true.)
        itemp=i8save(1,1,.true.)
        go to 50
c
c  set lerror and test for a previous unrecovered error.
c
 20   if (i8save(1,nerr,.true.).eq.0) go to 30
c
        write(iwunit,9002)
 9002   format(23h1error    3 in seterr -,
     1         48h an unrecovered error followed by another error.//
     2         48h the previous and current error messages follow.///)
        call eprint
        call e9rint(messg,nw,nerr,.true.)
        go to 50
c
c  save this message in case it is not recovered from properly.
c
 30   call e9rint(messg,nw,nerr,.true.)
c
      if (iopt.eq.1 .or. iopt.eq.2) go to 40
c
c  must have iopt = 1 or 2.
c
        write(iwunit,9003)
 9003   format(42h1error    4 in seterr - bad value for iopt//
     1         34h the current error message follows///)
        go to 50
c
c  test for recovery.
c
 40   if (iopt.eq.2) go to 50
c
      if (i8save(2,0,.false.).eq.1) return
c
c     call eprint
c     stop
c
 50   call eprint
 60   call fdump
      stop
c
      end
      subroutine e9rint(messg,nw,nerr,save)
c
c  this routine stores the current error message or prints the old one,
c  if any, depending on whether or not save = .true. .
c
      integer messg(nw)
      logical save
      external i1mach, i8save
c
c  messgp stores at least the first 72 characters of the previous
c  message. its length is machine dependent and must be at least
c
c       1 + 71/(the number of characters stored per integer word).
c
      integer messgp(36),fmt(14),ccplus
c
c  start with no previous message.
c
      data messgp(1)/1h1/, nwp/0/, nerrp/0/
c
c  set up the format for printing the error message.
c  the format is simply (a1,14x,72axx) where xx=i1mach(6) is the
c  number of characters stored per integer word.
c
      data ccplus  / 1h+ /
c
      data fmt( 1) / 1h( /
      data fmt( 2) / 1ha /
      data fmt( 3) / 1h1 /
      data fmt( 4) / 1h, /
      data fmt( 5) / 1h1 /
      data fmt( 6) / 1h4 /
      data fmt( 7) / 1hx /
      data fmt( 8) / 1h, /
      data fmt( 9) / 1h7 /
      data fmt(10) / 1h2 /
      data fmt(11) / 1ha /
      data fmt(12) / 1hx /
      data fmt(13) / 1hx /
      data fmt(14) / 1h) /
c
      if (.not.save) go to 20
c
c  save the message.
c
        nwp=nw
        nerrp=nerr
        do 10 i=1,nw
 10     messgp(i)=messg(i)
c
        go to 30
c
 20   if (i8save(1,0,.false.).eq.0) go to 30
c
c  print the message.
c
        iwunit=i1mach(4)
        write(iwunit,9000) nerrp
 9000   format(7h error ,i4,4h in )
c
        call s88fmt(2,i1mach(6),fmt(12))
        write(iwunit,fmt) ccplus,(messgp(i),i=1,nwp)
c
 30   return
c
      end

      subroutine fdump
      return
      end
C
C  ---------------------------------------
C  Compare scalar fnlib routines to vfnlib
C  ---------------------------------------
C
C  Ron Boisvert and Bonita Saunders
C  Computing and Applied Mathematics Laboratory
C  National Institute of Standards and Technology
C  Gaithersburg, MD 20899
C
C  boisvert@cam.nist.gov
C  saunders@cam.nist.gov
C
C  21 Feb 91
C
C----------------------------------------------------------------------
C
      INTEGER M
      PARAMETER ( M=2000 )
C
      REAL      XLO, XHI, X, F, FX, WORK
      INTEGER   I, IWORK, J
C
      DIMENSION XLO(3,8), XHI(3,8), X(M), F(M), FX(M), WORK(7*M),
     +          IWORK(M)
C
      SAVE XLO, XHI
C
      EXTERNAL  VI0, BESI0, VI1, BESI1, VJ0, BESJ0, VJ1, BESJ1,
     +          VK0, BESK0, VK1, BESK1, VY0, BESY0, VY1, BESY1
C
      DATA ((XLO(I,J),XHI(I,J),I=1,3),J=1,8) /
     +     0.0E0, 3.0E0 , -8.0E0, -3.0E0, 10.0E0, 50.0E0,
     +     0.0E0, 3.0E0 , -8.0E0, -3.0E0, 10.0E0, 50.0E0,
     +     0.0E0, 4.0E0 , -8.0E0, -4.0E0,  8.0E0, 50.0E0,
     +     0.0E0, 4.0E0 , -8.0E0, -4.0E0,  8.0E0, 50.0E0,
     +     0.0E0, 2.0E0 ,  2.0E0,  8.0E0, 10.0E0, 50.0E0,
     +     0.0E0, 2.0E0 ,  2.0E0,  8.0E0, 10.0E0, 50.0E0,
     +     0.0E0, 1.0E-3,  1.0E0,  4.0E0,  4.0E0, 50.0E0,
     +     0.0E0, 1.0E-3,  1.0E0,  4.0E0,  4.0E0, 50.0E0  /
C
C----------------------------------------------------------------------
C
C     ... TEST EACH USER-CALLABLE ROUTINE
C
      CALL TEST('I0',BESI0,VI0,XLO(1,1),XHI(1,1),M,X,F,FX,WORK,IWORK)
      CALL TEST('I1',BESI1,VI1,XLO(1,2),XHI(1,2),M,X,F,FX,WORK,IWORK)
      CALL TEST('J0',BESJ0,VJ0,XLO(1,3),XHI(1,3),M,X,F,FX,WORK,IWORK)
      CALL TEST('J1',BESJ1,VJ1,XLO(1,4),XHI(1,4),M,X,F,FX,WORK,IWORK)
      CALL TEST('K0',BESK0,VK0,XLO(1,5),XHI(1,5),M,X,F,FX,WORK,IWORK)
      CALL TEST('K1',BESK1,VK1,XLO(1,6),XHI(1,6),M,X,F,FX,WORK,IWORK)
      CALL TEST('Y0',BESY0,VY0,XLO(1,7),XHI(1,7),M,X,F,FX,WORK,IWORK)
      CALL TEST('Y1',BESY1,VY1,XLO(1,8),XHI(1,8),M,X,F,FX,WORK,IWORK)
C
      END

      SUBROUTINE TEST (NAME, SFUN, VSUB, XLOW, XHIGH, M, X, F, FX, 
     +                 WORK, IWORK)
C
C----------------------------------------------------------------------
C
C  compare scalar function sfun to vectorized subroutine vsub for m
C  arguments distributed in various ways in the intervals
C  (xlo(i),xhi(i)), i=1,2,3.
C
C----------------------------------------------------------------------
C
C  ----------
C  PARAMETERS
C  ----------
C
      REAL    SFUN, XLOW, XHIGH, X, F, FX, WORK
      INTEGER IWORK, M
      CHARACTER*(*) NAME
C
      DIMENSION XLOW(3), XHIGH(3), X(M), F(M), FX(M), WORK(6*M),
     +          IWORK(M)
C
      EXTERNAL SFUN, VSUB
C
C  ---------------
C  LOCAL VARIABLES
C  ---------------
C
      LOGICAL DEBUG
      INTEGER NREP, NTRY
      PARAMETER ( DEBUG=.FALSE., NREP=15, NTRY=19 )
C
      INTEGER I, INFO, IREP, IOFF, IPC(NTRY), JPC(NTRY), K, KPC,
     +        MA, MB, MC, NFAIL, NREPV
      REAL    DIFF, DMAX, DX, EPS, R1MACH
      REAL SECOND, SU, SUMAX, SUMIN, SUAVE, T0, T1, TS, TSAVE, TV,
     +        TVAVE
C
      SAVE IPC, JPC
C
      DATA IPC /  0,  0,  0,  0,  0, 25, 25, 25, 25, 33, 50, 50, 50,
     +           75, 75,100,  1,  1, 98/
      DATA JPC /  0, 25, 50, 75,100,  0, 25, 50, 75, 33,  0, 25, 50,
     +            0, 25,  0,  1, 98,  1/
C
C----------------------------------------------------------------------
C
C  ntry = number of trials (each with different argument distribution)
C  nrep = number of times test repeated to get reliable timings
C
C  the distribution of arguments in trial k is
C
C     ipc(k)%            in (xlow(1),xhigh(1))
C     jpc(k)%            in (xlow(2),xhigh(2))
C     100-ipc(k)-jpc(k)% in (xlow(3),xhigh(3))
C
C  scalar results (fx) and vector results (f) are compared and accepted
C  if   abs((f - fx)/fx) < eps  for fx > 1  (relative difference)
C  or   abs(f - fx)      < eps  for fx < 1  (absolute difference)
C
C----------------------------------------------------------------------
C
      EPS = 5.0E0*R1MACH(4)

C     ... TRIGGER INITIALIZATIONS
C
      X(1) = 1.0E0
      CALL VSUB(1,X,F,WORK,IWORK,INFO)
      FX(1) = SFUN(X(1))
C
      PRINT *
      PRINT *
      PRINT *,' ----------'
      PRINT *,' TEST OF ',NAME
      PRINT *,' ----------'
      PRINT *
      PRINT *,' Vector Length = ',M
      PRINT *,' Repetitions   = ',NREP
      PRINT *,' Eps           = ',EPS
      PRINT *
      PRINT *,' Range A = (',XLOW(1),',',XHIGH(1),')'
      PRINT *,' Range B = (',XLOW(2),',',XHIGH(2),')'
      PRINT *,' Range C = (',XLOW(3),',',XHIGH(3),')'
      PRINT *
      PRINT *,' %A      = percent arguments in range A'
      PRINT *,' %B      = percent arguments in range B'
      PRINT *,' %C      = percent arguments in range C'
      PRINT *,' Scalar  = time for scalar code (FNLIB)'
      PRINT *,' Vector  = time for vector code (VFNLIB)'
      PRINT *,' SU      = speedup (Scalar/Vector)'
      PRINT *,' Nfail   = number of scalar-vector differences .gt. eps'
      PRINT *,' Dmax    = maximum scalar-vector difference'
      PRINT *
      PRINT *,' All times in seconds/argument.'
      PRINT *
      PRINT *,' Scalar-vector differences measured by'
      PRINT *,'     abs(s-v)/max(abs(s),1)'
      PRINT *,' where s,v are scalar,vector results.'
      PRINT *
      PRINT *,'  %A  %B  %C    Scalar    Vector   SU   Nfail   Dmax'
      PRINT *,' -------------------------------------------------------'
C
      SUMAX = 0.0
      SUMIN = 1000.0
      SUAVE = 0.0
      TSAVE = 0.0
      TVAVE = 0.0
C
C     ... PERFORM EACH TEST
C
      DO 100 K=1,NTRY
C
         MA = IPC(K)*M/100
         MB = JPC(K)*M/100
         KPC = 100 - IPC(K) - JPC(K)
         MC = M - MA - MB
C
C        ... SELECT X VECTOR
C
         IOFF = 0
C
         IF (MA .GT. 0) THEN
            DX = (XHIGH(1) - XLOW(1))/MA
            DO 10 I=1,MA
               X(I) = XLOW(1) + I*DX
   10       CONTINUE
            IOFF = IOFF + MA
         ENDIF
C
         IF (MB .GT. 0) THEN
            DX = (XHIGH(2) - XLOW(2))/MB
            DO 20 I=1,MB
               X(I+IOFF) = XLOW(2) + I*DX
   20       CONTINUE
            IOFF = IOFF + MB
         ENDIF
C
         IF (MC .GT. 0) THEN
            DX = (XHIGH(3) - XLOW(3))/MC
            DO 30 I=1,MC
               X(I+IOFF) = XLOW(3) + I*DX
   30       CONTINUE
         ENDIF
C
C        ... RUN VECTOR CODE
C
	 NREPV = 5*NREP
         T0 = SECOND()
         DO 40 IREP=1,NREPV
            CALL VSUB(M,X,F,WORK,IWORK,INFO)
   40    CONTINUE
         T1 = SECOND()
         TV = (T1 - T0)/NREPV
         IF (INFO .LT. 0) PRINT *,'Warning in ',NAME,' : Info = ',INFO
         IF (INFO .GT. 0) PRINT *,'Failure in ',NAME,' : Info = ',INFO,
     +                            '  Position = ',IWORK(1)
C
C        ... RUN SCALAR CODE
C
         T0 = SECOND()
         DO 60 IREP=1,NREP
            DO 50 I=1,M
               FX(I) = SFUN(X(I))
   50       CONTINUE
   60    CONTINUE
         T1 = SECOND()
         TS = (T1 - T0)/NREP
C
C        ... COMPUTE TIMINGS
C
         SU = TS/TV
         SUMIN = MIN(SU,SUMIN)
         SUMAX = MAX(SU,SUMAX)
         SUAVE = SUAVE + SU
         TSAVE = TSAVE + TS
         TVAVE = TVAVE + TV
C
C        ... EVALUATE RESULTS
C
         DMAX = 0.0E0
         NFAIL = 0
         DO 70 I=1,M
            DIFF = ABS(FX(I)-F(I))/MAX(ABS(FX(I)),1.0E0)
            DMAX = MAX(DIFF,DMAX)
            IF (DIFF .GT. EPS) THEN
               NFAIL = NFAIL + 1
               IF (DEBUG) PRINT *,'Failure ',I,X(I),F(I),FX(I),DIFF
            ENDIF
  70     CONTINUE
C
C         ... PRINT RESULTS
C
          PRINT 1000, IPC(K),JPC(K),KPC,TS/M,TV/M,SU,NFAIL,DMAX
C
  100 CONTINUE
C
      SUAVE = SUAVE/NTRY
      TSAVE = TSAVE/M/NTRY
      TVAVE = TVAVE/M/NTRY
C
      PRINT *
      PRINT 1002,'Ave Scalar   = ',TSAVE,' sec/argument'
      PRINT 1002,'Ave Vector   = ',TVAVE,' sec/argument'
      PRINT *
      PRINT 1001,'Max speed up = ',SUMAX
      PRINT 1001,'Min speed up = ',SUMIN
      PRINT *
      PRINT 1001,'Ave speed up = ',SUAVE
C
      RETURN
 1000 FORMAT(3(1X,I3),1P,2(3X,E7.1),0P,1X,F5.1,1X,I4,3X,1P,E7.1)
 1001 FORMAT(1X,A,F5.1,A)
 1002 FORMAT(1X,A,1P,E7.1,A)
      END
      INTEGER FUNCTION I1MACH(I)
C***BEGIN PROLOGUE  I1MACH
C***DATE WRITTEN   750101   (YYMMDD)
C***REVISION DATE  910131   (YYMMDD)
C***CATEGORY NO.  R1
C***KEYWORDS  MACHINE CONSTANTS
C***AUTHOR  FOX, P. A., (BELL LABS)
C           HALL, A. D., (BELL LABS)
C           SCHRYER, N. L., (BELL LABS)
C***PURPOSE  Returns integer machine dependent constants
C***DESCRIPTION
C
C     This is the CMLIB version of I1MACH, the integer machine
C     constants subroutine originally developed for the PORT library.
C
C     I1MACH can be used to obtain machine-dependent parameters
C     for the local machine environment.  It is a function
C     subroutine with one (input) argument, and can be called
C     as follows, for example
C
C          K = I1MACH(I)
C
C     where I=1,...,16.  The (output) value of K above is
C     determined by the (input) value of I.  The results for
C     various values of I are discussed below.
C
C  I/O unit numbers.
C    I1MACH( 1) = the standard input unit.
C    I1MACH( 2) = the standard output unit.
C    I1MACH( 3) = the standard punch unit.
C    I1MACH( 4) = the standard error message unit.
C
C  Words.
C    I1MACH( 5) = the number of bits per integer storage unit.
C    I1MACH( 6) = the number of characters per integer storage unit.
C
C  Integers.
C    assume integers are represented in the S-digit, base-A form
C
C               sign ( X(S-1)*A**(S-1) + ... + X(1)*A + X(0) )
C
C               where 0 .LE. X(I) .LT. A for I=0,...,S-1.
C    I1MACH( 7) = A, the base.
C    I1MACH( 8) = S, the number of base-A digits.
C    I1MACH( 9) = A**S - 1, the largest magnitude.
C
C  Floating-Point Numbers.
C    Assume floating-point numbers are represented in the T-digit,
C    base-B form
C               sign (B**E)*( (X(1)/B) + ... + (X(T)/B**T) )
C
C               where 0 .LE. X(I) .LT. B for I=1,...,T,
C               0 .LT. X(1), and EMIN .LE. E .LE. EMAX.
C    I1MACH(10) = B, the base.
C
C  Single-Precision
C    I1MACH(11) = T, the number of base-B digits.
C    I1MACH(12) = EMIN, the smallest exponent E.
C    I1MACH(13) = EMAX, the largest exponent E.
C
C  Double-Precision
C    I1MACH(14) = T, the number of base-B digits.
C    I1MACH(15) = EMIN, the smallest exponent E.
C    I1MACH(16) = EMAX, the largest exponent E.
C
C  To alter this function for a particular environment,
C  the desired set of DATA statements should be activated by
C  removing the C from column 1.  Also, the values of
C  I1MACH(1) - I1MACH(4) should be checked for consistency
C  with the local operating system.
C***REFERENCES  FOX P.A., HALL A.D., SCHRYER N.L.,*FRAMEWORK FOR A
C                 PORTABLE LIBRARY*, ACM TRANSACTIONS ON MATHEMATICAL
C                 SOFTWARE, VOL. 4, NO. 2, JUNE 1978, PP. 177-188.
C***ROUTINES CALLED  (NONE)
C***END PROLOGUE  I1MACH
C
      INTEGER IMACH(16),OUTPUT
      EQUIVALENCE (IMACH(4),OUTPUT)
C
C     MACHINE CONSTANTS FOR IEEE ARITHMETIC MACHINES, SUCH AS THE AT&T
C     3B SERIES, MOTOROLA 68000 BASED MACHINES (E.G. SUN 3 AND AT&T
C     PC 7300), AND 8087 BASED MICROS (E.G. IBM PC AND AT&T 6300).
C
C === MACHINE = IEEE.MOST-SIG-BYTE-FIRST
C === MACHINE = IEEE.LEAST-SIG-BYTE-FIRST
C === MACHINE = SUN
C === MACHINE = 68000
C === MACHINE = 8087
C === MACHINE = IBM.PC
C === MACHINE = ATT.3B
C === MACHINE = ATT.7300
C === MACHINE = ATT.6300
       DATA IMACH( 1) /    5 /
       DATA IMACH( 2) /    6 /
       DATA IMACH( 3) /    7 /
       DATA IMACH( 4) /    6 /
       DATA IMACH( 5) /   32 /
       DATA IMACH( 6) /    4 /
       DATA IMACH( 7) /    2 /
       DATA IMACH( 8) /   31 /
       DATA IMACH( 9) / 2147483647 /
       DATA IMACH(10) /    2 /
       DATA IMACH(11) /   24 /
       DATA IMACH(12) / -125 /
       DATA IMACH(13) /  128 /
       DATA IMACH(14) /   53 /
       DATA IMACH(15) / -1021 /
       DATA IMACH(16) /  1024 /
C
C     MACHINE CONSTANTS FOR AMDAHL MACHINES.
C
C === MACHINE = AMDAHL
C      DATA IMACH( 1) /   5 /
C      DATA IMACH( 2) /   6 /
C      DATA IMACH( 3) /   7 /
C      DATA IMACH( 4) /   6 /
C      DATA IMACH( 5) /  32 /
C      DATA IMACH( 6) /   4 /
C      DATA IMACH( 7) /   2 /
C      DATA IMACH( 8) /  31 /
C      DATA IMACH( 9) / 2147483647 /
C      DATA IMACH(10) /  16 /
C      DATA IMACH(11) /   6 /
C      DATA IMACH(12) / -64 /
C      DATA IMACH(13) /  63 /
C      DATA IMACH(14) /  14 /
C      DATA IMACH(15) / -64 /
C      DATA IMACH(16) /  63 /
C
C     MACHINE CONSTANTS FOR THE BURROUGHS 1700 SYSTEM.
C
C === MACHINE = BURROUGHS.1700
C      DATA IMACH( 1) /    7 /
C      DATA IMACH( 2) /    2 /
C      DATA IMACH( 3) /    2 /
C      DATA IMACH( 4) /    2 /
C      DATA IMACH( 5) /   36 /
C      DATA IMACH( 6) /    4 /
C      DATA IMACH( 7) /    2 /
C      DATA IMACH( 8) /   33 /
C      DATA IMACH( 9) / Z1FFFFFFFF /
C      DATA IMACH(10) /    2 /
C      DATA IMACH(11) /   24 /
C      DATA IMACH(12) / -256 /
C      DATA IMACH(13) /  255 /
C      DATA IMACH(14) /   60 /
C      DATA IMACH(15) / -256 /
C      DATA IMACH(16) /  255 /
C
C     MACHINE CONSTANTS FOR THE BURROUGHS 5700 SYSTEM.
C
C === MACHINE = BURROUGHS.5700
C      DATA IMACH( 1) /   5 /
C      DATA IMACH( 2) /   6 /
C      DATA IMACH( 3) /   7 /
C      DATA IMACH( 4) /   6 /
C      DATA IMACH( 5) /  48 /
C      DATA IMACH( 6) /   6 /
C      DATA IMACH( 7) /   2 /
C      DATA IMACH( 8) /  39 /
C      DATA IMACH( 9) / O0007777777777777 /
C      DATA IMACH(10) /   8 /
C      DATA IMACH(11) /  13 /
C      DATA IMACH(12) / -50 /
C      DATA IMACH(13) /  76 /
C      DATA IMACH(14) /  26 /
C      DATA IMACH(15) / -50 /
C      DATA IMACH(16) /  76 /
C
C     MACHINE CONSTANTS FOR THE BURROUGHS 6700/7700 SYSTEMS.
C
C === MACHINE = BURROUGHS.6700
C === MACHINE = BURROUGHS.7700
C      DATA IMACH( 1) /   5 /
C      DATA IMACH( 2) /   6 /
C      DATA IMACH( 3) /   7 /
C      DATA IMACH( 4) /   6 /
C      DATA IMACH( 5) /  48 /
C      DATA IMACH( 6) /   6 /
C      DATA IMACH( 7) /   2 /
C      DATA IMACH( 8) /  39 /
C      DATA IMACH( 9) / O0007777777777777 /
C      DATA IMACH(10) /   8 /
C      DATA IMACH(11) /  13 /
C      DATA IMACH(12) / -50 /
C      DATA IMACH(13) /  76 /
C      DATA IMACH(14) /  26 /
C      DATA IMACH(15) / -32754 /
C      DATA IMACH(16) /  32780 /
C
C     MACHINE CONSTANTS FOR THE CONVEX C-120 (NATIVE MODE)
C
C === MACHINE = CONVEX.C1
C      DATA IMACH( 1) /    5 /
C      DATA IMACH( 2) /    6 /
C      DATA IMACH( 3) /    0 /
C      DATA IMACH( 4) /    6 /
C      DATA IMACH( 5) /   32 /
C      DATA IMACH( 6) /    4 /
C      DATA IMACH( 7) /    2 /
C      DATA IMACH( 8) /   31 /
C      DATA IMACH( 9) / 2147483647 /
C      DATA IMACH(10) /    2 /
C      DATA IMACH(11) /   24 /
C      DATA IMACH(12) / -127 /
C      DATA IMACH(13) /  127 /
C      DATA IMACH(14) /   53 /
C      DATA IMACH(15) / -1023 /
C      DATA IMACH(16) /  1023 /
C
C     MACHINE CONSTANTS FOR THE CONVEX (NATIVE MODE)
C     WITH -R8 OPTION
C
C === MACHINE = CONVEX.C1.R8
C      DATA IMACH( 1) /     5 /
C      DATA IMACH( 2) /     6 /
C      DATA IMACH( 3) /     0 /
C      DATA IMACH( 4) /     6 /
C      DATA IMACH( 5) /    32 /
C      DATA IMACH( 6) /     4 /
C      DATA IMACH( 7) /     2 /
C      DATA IMACH( 8) /    31 /
C      DATA IMACH( 9) / 2147483647 /
C      DATA IMACH(10) /     2 /
C      DATA IMACH(11) /    53 /
C      DATA IMACH(12) / -1023 /
C      DATA IMACH(13) /  1023 /
C      DATA IMACH(14) /    53 /
C      DATA IMACH(15) / -1023 /
C      DATA IMACH(16) /  1023 /
C
C     MACHINE CONSTANTS FOR THE CONVEX C-120 (IEEE MODE)
C
C === MACHINE = CONVEX.C1.IEEE
C      DATA IMACH( 1) /    5 /
C      DATA IMACH( 2) /    6 /
C      DATA IMACH( 3) /    0 /
C      DATA IMACH( 4) /    6 /
C      DATA IMACH( 5) /   32 /
C      DATA IMACH( 6) /    4 /
C      DATA IMACH( 7) /    2 /
C      DATA IMACH( 8) /   31 /
C      DATA IMACH( 9) / 2147483647 /
C      DATA IMACH(10) /    2 /
C      DATA IMACH(11) /   24 /
C      DATA IMACH(12) / -125 /
C      DATA IMACH(13) /  128 /
C      DATA IMACH(14) /   53 /
C      DATA IMACH(15) / -1021 /
C      DATA IMACH(16) /  1024 /
C
C     MACHINE CONSTANTS FOR THE CONVEX (IEEE MODE)
C     WITH -R8 OPTION
C
C === MACHINE = CONVEX.C1.IEEE.R8
C      DATA IMACH( 1) /     5 /
C      DATA IMACH( 2) /     6 /
C      DATA IMACH( 3) /     0 /
C      DATA IMACH( 4) /     6 /
C      DATA IMACH( 5) /    32 /
C      DATA IMACH( 6) /     4 /
C      DATA IMACH( 7) /     2 /
C      DATA IMACH( 8) /    31 /
C      DATA IMACH( 9) / 2147483647 /
C      DATA IMACH(10) /     2 /
C      DATA IMACH(11) /    53 /
C      DATA IMACH(12) / -1021 /
C      DATA IMACH(13) /  1024 /
C      DATA IMACH(14) /    53 /
C      DATA IMACH(15) / -1021 /
C      DATA IMACH(16) /  1024 /
C
C     MACHINE CONSTANTS FOR THE CYBER 170/180 SERIES USING NOS (FTN5).
C
C === MACHINE = CYBER.170.NOS
C === MACHINE = CYBER.180.NOS
C      DATA IMACH( 1) /    5 /
C      DATA IMACH( 2) /    6 /
C      DATA IMACH( 3) /    7 /
C      DATA IMACH( 4) /    6 /
C      DATA IMACH( 5) /   60 /
C      DATA IMACH( 6) /   10 /
C      DATA IMACH( 7) /    2 /
C      DATA IMACH( 8) /   48 /
C      DATA IMACH( 9) / O"00007777777777777777" /
C      DATA IMACH(10) /    2 /
C      DATA IMACH(11) /   48 /
C      DATA IMACH(12) / -974 /
C      DATA IMACH(13) / 1070 /
C      DATA IMACH(14) /   96 /
C      DATA IMACH(15) / -927 /
C      DATA IMACH(16) / 1070 /
C
C     MACHINE CONSTANTS FOR THE CDC 180 SERIES USING NOS/VE
C
C === MACHINE = CYBER.180.NOS/VE
C      DATA IMACH( 1) /     5 /
C      DATA IMACH( 2) /     6 /
C      DATA IMACH( 3) /     7 /
C      DATA IMACH( 4) /     6 /
C      DATA IMACH( 5) /    64 /
C      DATA IMACH( 6) /     8 /
C      DATA IMACH( 7) /     2 /
C      DATA IMACH( 8) /    63 /
C      DATA IMACH( 9) / 9223372036854775807 /
C      DATA IMACH(10) /     2 /
C      DATA IMACH(11) /    47 /
C      DATA IMACH(12) / -4095 /
C      DATA IMACH(13) /  4094 /
C      DATA IMACH(14) /    94 /
C      DATA IMACH(15) / -4095 /
C      DATA IMACH(16) /  4094 /
C
C     MACHINE CONSTANTS FOR THE CYBER 205
C
C === MACHINE = CYBER.205
C      DATA IMACH( 1) /      5 /
C      DATA IMACH( 2) /      6 /
C      DATA IMACH( 3) /      7 /
C      DATA IMACH( 4) /      6 /
C      DATA IMACH( 5) /     64 /
C      DATA IMACH( 6) /      8 /
C      DATA IMACH( 7) /      2 /
C      DATA IMACH( 8) /     47 /
C      DATA IMACH( 9) / X'00007FFFFFFFFFFF' /
C      DATA IMACH(10) /      2 /
C      DATA IMACH(11) /     47 /
C      DATA IMACH(12) / -28625 /
C      DATA IMACH(13) /  28718 /
C      DATA IMACH(14) /     94 /
C      DATA IMACH(15) / -28625 /
C      DATA IMACH(16) /  28718 /
C
C     MACHINE CONSTANTS FOR THE CDC 6000/7000 SERIES.
C
C === MACHINE = CDC.6000
C === MACHINE = CDC.7000
C      DATA IMACH( 1) /    5 /
C      DATA IMACH( 2) /    6 /
C      DATA IMACH( 3) /    7 /
C      DATA IMACH( 4) /    6 /
C      DATA IMACH( 5) /   60 /
C      DATA IMACH( 6) /   10 /
C      DATA IMACH( 7) /    2 /
C      DATA IMACH( 8) /   48 /
C      DATA IMACH( 9) / 00007777777777777777B /
C      DATA IMACH(10) /    2 /
C      DATA IMACH(11) /   48 /
C      DATA IMACH(12) / -974 /
C      DATA IMACH(13) / 1070 /
C      DATA IMACH(14) /   96 /
C      DATA IMACH(15) / -927 /
C      DATA IMACH(16) / 1070 /
C
C     MACHINE CONSTANTS FOR THE CRAY 1, XMP, 2, AND 3.
C     USING THE 46 BIT INTEGER COMPILER OPTION
C
C === MACHINE = CRAY.46-BIT-INTEGER
C      DATA IMACH( 1) /     5 /
C      DATA IMACH( 2) /     6 /
C      DATA IMACH( 3) /   102 /
C      DATA IMACH( 4) /     6 /
C      DATA IMACH( 5) /    64 /
C      DATA IMACH( 6) /     8 /
C      DATA IMACH( 7) /     2 /
C      DATA IMACH( 8) /    46 /
C      DATA IMACH( 9) /  777777777777777777777B /
C      DATA IMACH(10) /     2 /
C      DATA IMACH(11) /    47 /
C      DATA IMACH(12) / -8189 /
C      DATA IMACH(13) /  8190 /
C      DATA IMACH(14) /    94 /
C      DATA IMACH(15) / -8099 /
C      DATA IMACH(16) /  8190 /
C
C     MACHINE CONSTANTS FOR THE CRAY 1, XMP, 2, AND 3.
C     USING THE 64 BIT INTEGER COMPILER OPTION
C
C === MACHINE = CRAY.64-BIT-INTEGER
C      DATA IMACH( 1) /     5 /
C      DATA IMACH( 2) /     6 /
C      DATA IMACH( 3) /   102 /
C      DATA IMACH( 4) /     6 /
C      DATA IMACH( 5) /    64 /
C      DATA IMACH( 6) /     8 /
C      DATA IMACH( 7) /     2 /
C      DATA IMACH( 8) /    63 /
C      DATA IMACH( 9) /  777777777777777777777B /
C      DATA IMACH(10) /     2 /
C      DATA IMACH(11) /    47 /
C      DATA IMACH(12) / -8189 /
C      DATA IMACH(13) /  8190 /
C      DATA IMACH(14) /    94 /
C      DATA IMACH(15) / -8099 /
C      DATA IMACH(16) /  8190 /C
C     MACHINE CONSTANTS FOR THE DATA GENERAL ECLIPSE S/200
C
C === MACHINE = DATA_GENERAL.ECLIPSE.S/200
C      DATA IMACH( 1) /   11 /
C      DATA IMACH( 2) /   12 /
C      DATA IMACH( 3) /    8 /
C      DATA IMACH( 4) /   10 /
C      DATA IMACH( 5) /   16 /
C      DATA IMACH( 6) /    2 /
C      DATA IMACH( 7) /    2 /
C      DATA IMACH( 8) /   15 /
C      DATA IMACH( 9) /32767 /
C      DATA IMACH(10) /   16 /
C      DATA IMACH(11) /    6 /
C      DATA IMACH(12) /  -64 /
C      DATA IMACH(13) /   63 /
C      DATA IMACH(14) /   14 /
C      DATA IMACH(15) /  -64 /
C      DATA IMACH(16) /   63 /
C
C     ELXSI  6400
C
C === MACHINE = ELSXI.6400
C      DATA IMACH( 1) /     5 /
C      DATA IMACH( 2) /     6 /
C      DATA IMACH( 3) /     6 /
C      DATA IMACH( 4) /     6 /
C      DATA IMACH( 5) /    32 /
C      DATA IMACH( 6) /     4 /
C      DATA IMACH( 7) /     2 /
C      DATA IMACH( 8) /    32 /
C      DATA IMACH( 9) / 2147483647 /
C      DATA IMACH(10) /     2 /
C      DATA IMACH(11) /    24 /
C      DATA IMACH(12) /  -126 /
C      DATA IMACH(13) /   127 /
C      DATA IMACH(14) /    53 /
C      DATA IMACH(15) / -1022 /
C      DATA IMACH(16) /  1023 /
C
C     MACHINE CONSTANTS FOR THE HARRIS 220
C     MACHINE CONSTANTS FOR THE HARRIS SLASH 6 AND SLASH 7
C
C === MACHINE = HARRIS.220
C === MACHINE = HARRIS.SLASH6
C === MACHINE = HARRIS.SLASH7
C      DATA IMACH( 1) /       5 /
C      DATA IMACH( 2) /       6 /
C      DATA IMACH( 3) /       0 /
C      DATA IMACH( 4) /       6 /
C      DATA IMACH( 5) /      24 /
C      DATA IMACH( 6) /       3 /
C      DATA IMACH( 7) /       2 /
C      DATA IMACH( 8) /      23 /
C      DATA IMACH( 9) / 8388607 /
C      DATA IMACH(10) /       2 /
C      DATA IMACH(11) /      23 /
C      DATA IMACH(12) /    -127 /
C      DATA IMACH(13) /     127 /
C      DATA IMACH(14) /      38 /
C      DATA IMACH(15) /    -127 /
C      DATA IMACH(16) /     127 /
C
C     MACHINE CONSTANTS FOR THE HONEYWELL 600/6000 SERIES.
C     MACHINE CONSTANTS FOR THE HONEYWELL DPS 8/70 SERIES.
C
C === MACHINE = HONEYWELL.600/6000
C === MACHINE = HONEYWELL.DPS.8/70
C      DATA IMACH( 1) /    5 /
C      DATA IMACH( 2) /    6 /
C      DATA IMACH( 3) /   43 /
C      DATA IMACH( 4) /    6 /
C      DATA IMACH( 5) /   36 /
C      DATA IMACH( 6) /    4 /
C      DATA IMACH( 7) /    2 /
C      DATA IMACH( 8) /   35 /
C      DATA IMACH( 9) / O377777777777 /
C      DATA IMACH(10) /    2 /
C      DATA IMACH(11) /   27 /
C      DATA IMACH(12) / -127 /
C      DATA IMACH(13) /  127 /
C      DATA IMACH(14) /   63 /
C      DATA IMACH(15) / -127 /
C      DATA IMACH(16) /  127 /
C
C     MACHINE CONSTANTS FOR THE HP 2100
C     3 WORD DOUBLE PRECISION OPTION WITH FTN4
C
C === MACHINE = HP.2100.3_WORD_DP
C      DATA IMACH(1) /      5/
C      DATA IMACH(2) /      6 /
C      DATA IMACH(3) /      4 /
C      DATA IMACH(4) /      1 /
C      DATA IMACH(5) /     16 /
C      DATA IMACH(6) /      2 /
C      DATA IMACH(7) /      2 /
C      DATA IMACH(8) /     15 /
C      DATA IMACH(9) /  32767 /
C      DATA IMACH(10)/      2 /
C      DATA IMACH(11)/     23 /
C      DATA IMACH(12)/   -128 /
C      DATA IMACH(13)/    127 /
C      DATA IMACH(14)/     39 /
C      DATA IMACH(15)/   -128 /
C      DATA IMACH(16)/    127 /
C
C     MACHINE CONSTANTS FOR THE HP 2100
C     4 WORD DOUBLE PRECISION OPTION WITH FTN4
C
C === MACHINE = HP.2100.4_WORD_DP
C      DATA IMACH(1) /      5 /
C      DATA IMACH(2) /      6 /
C      DATA IMACH(3) /      4 /
C      DATA IMACH(4) /      1 /
C      DATA IMACH(5) /     16 /
C      DATA IMACH(6) /      2 /
C      DATA IMACH(7) /      2 /
C      DATA IMACH(8) /     15 /
C      DATA IMACH(9) /  32767 /
C      DATA IMACH(10)/      2 /
C      DATA IMACH(11)/     23 /
C      DATA IMACH(12)/   -128 /
C      DATA IMACH(13)/    127 /
C      DATA IMACH(14)/     55 /
C      DATA IMACH(15)/   -128 /
C      DATA IMACH(16)/    127 /
C
C     HP 9000
C
C === MACHINE = HP.9000
C      DATA IMACH( 1) /     5 /
C      DATA IMACH( 2) /     6 /
C      DATA IMACH( 3) /     6 /
C      DATA IMACH( 4) /     7 /
C      DATA IMACH( 5) /    32 /
C      DATA IMACH( 6) /     4 /
C      DATA IMACH( 7) /     2 /
C      DATA IMACH( 8) /    32 /
C      DATA IMACH( 9) / 2147483647 /
C      DATA IMACH(10) /     2 /
C      DATA IMACH(11) /    24 /
C      DATA IMACH(12) /  -126 /
C      DATA IMACH(13) /   127 /
C      DATA IMACH(14) /    53 /
C      DATA IMACH(15) / -1015 /
C      DATA IMACH(16) /  1017 /
C
C     MACHINE CONSTANTS FOR THE IBM 360/370 SERIES,
C     THE XEROX SIGMA 5/7/9 AND THE SEL SYSTEMS 85/86 AND
C     THE INTERDATA 3230 AND INTERDATA 7/32.
C
C === MACHINE = IBM.360
C === MACHINE = IBM.370
C === MACHINE = XEROX.SIGMA.5
C === MACHINE = XEROX.SIGMA.7
C === MACHINE = XEROX.SIGMA.9
C === MACHINE = SEL.85
C === MACHINE = SEL.86
C === MACHINE = INTERDATA.3230
C === MACHINE = INTERDATA.7/32
C      DATA IMACH( 1) /   5 /
C      DATA IMACH( 2) /   6 /
C      DATA IMACH( 3) /   7 /
C      DATA IMACH( 4) /   6 /
C      DATA IMACH( 5) /  32 /
C      DATA IMACH( 6) /   4 /
C      DATA IMACH( 7) /   2 /
C      DATA IMACH( 8) /  31 /
C      DATA IMACH( 9) / Z7FFFFFFF /
C      DATA IMACH(10) /  16 /
C      DATA IMACH(11) /   6 /
C      DATA IMACH(12) / -64 /
C      DATA IMACH(13) /  63 /
C      DATA IMACH(14) /  14 /
C      DATA IMACH(15) / -64 /
C      DATA IMACH(16) /  63 /
C
C     MACHINE CONSTANTS FOR THE INTERDATA 8/32
C     WITH THE UNIX SYSTEM FORTRAN 77 COMPILER.
C
C     FOR THE INTERDATA FORTRAN VII COMPILER REPLACE
C     THE Z'S SPECIFYING HEX CONSTANTS WITH Y'S.
C
C === MACHINE = INTERDATA.8/32.UNIX
C      DATA IMACH( 1) /   5 /
C      DATA IMACH( 2) /   6 /
C      DATA IMACH( 3) /   6 /
C      DATA IMACH( 4) /   6 /
C      DATA IMACH( 5) /  32 /
C      DATA IMACH( 6) /   4 /
C      DATA IMACH( 7) /   2 /
C      DATA IMACH( 8) /  31 /
C      DATA IMACH( 9) / Z'7FFFFFFF' /
C      DATA IMACH(10) /  16 /
C      DATA IMACH(11) /   6 /
C      DATA IMACH(12) / -64 /
C      DATA IMACH(13) /  62 /
C      DATA IMACH(14) /  14 /
C      DATA IMACH(15) / -64 /
C      DATA IMACH(16) /  62 /
C
C     MACHINE CONSTANTS FOR THE PDP-10 (KA PROCESSOR).
C
C === MACHINE = PDP-10.KA
C      DATA IMACH( 1) /    5 /
C      DATA IMACH( 2) /    6 /
C      DATA IMACH( 3) /    7 /
C      DATA IMACH( 4) /    6 /
C      DATA IMACH( 5) /   36 /
C      DATA IMACH( 6) /    5 /
C      DATA IMACH( 7) /    2 /
C      DATA IMACH( 8) /   35 /
C      DATA IMACH( 9) / "377777777777 /
C      DATA IMACH(10) /    2 /
C      DATA IMACH(11) /   27 /
C      DATA IMACH(12) / -128 /
C      DATA IMACH(13) /  127 /
C      DATA IMACH(14) /   54 /
C      DATA IMACH(15) / -101 /
C      DATA IMACH(16) /  127 /
C
C     MACHINE CONSTANTS FOR THE PDP-10 (KI PROCESSOR).
C
C === MACHINE = PDP-10.KI
C      DATA IMACH( 1) /    5 /
C      DATA IMACH( 2) /    6 /
C      DATA IMACH( 3) /    7 /
C      DATA IMACH( 4) /    6 /
C      DATA IMACH( 5) /   36 /
C      DATA IMACH( 6) /    5 /
C      DATA IMACH( 7) /    2 /
C      DATA IMACH( 8) /   35 /
C      DATA IMACH( 9) / "377777777777 /
C      DATA IMACH(10) /    2 /
C      DATA IMACH(11) /   27 /
C      DATA IMACH(12) / -128 /
C      DATA IMACH(13) /  127 /
C      DATA IMACH(14) /   62 /
C      DATA IMACH(15) / -128 /
C      DATA IMACH(16) /  127 /
C
C     MACHINE CONSTANTS FOR PDP-11 FORTRANS SUPPORTING
C     32-BIT INTEGER ARITHMETIC.
C
C === MACHINE = PDP-11.32-BIT
C      DATA IMACH( 1) /    5 /
C      DATA IMACH( 2) /    6 /
C      DATA IMACH( 3) /    7 /
C      DATA IMACH( 4) /    6 /
C      DATA IMACH( 5) /   32 /
C      DATA IMACH( 6) /    4 /
C      DATA IMACH( 7) /    2 /
C      DATA IMACH( 8) /   31 /
C      DATA IMACH( 9) / 2147483647 /
C      DATA IMACH(10) /    2 /
C      DATA IMACH(11) /   24 /
C      DATA IMACH(12) / -127 /
C      DATA IMACH(13) /  127 /
C      DATA IMACH(14) /   56 /
C      DATA IMACH(15) / -127 /
C      DATA IMACH(16) /  127 /
C
C     MACHINE CONSTANTS FOR PDP-11 FORTRANS SUPPORTING
C     16-BIT INTEGER ARITHMETIC.
C
C === MACHINE = PDP-11.16-BIT
C      DATA IMACH( 1) /    5 /
C      DATA IMACH( 2) /    6 /
C      DATA IMACH( 3) /    7 /
C      DATA IMACH( 4) /    6 /
C      DATA IMACH( 5) /   16 /
C      DATA IMACH( 6) /    2 /
C      DATA IMACH( 7) /    2 /
C      DATA IMACH( 8) /   15 /
C      DATA IMACH( 9) / 32767 /
C      DATA IMACH(10) /    2 /
C      DATA IMACH(11) /   24 /
C      DATA IMACH(12) / -127 /
C      DATA IMACH(13) /  127 /
C      DATA IMACH(14) /   56 /
C      DATA IMACH(15) / -127 /
C      DATA IMACH(16) /  127 /
C
C     MACHINE CONSTANTS FOR THE SEQUENT BALANCE 8000.
C
C === MACHINE = SEQUENT.BALANCE.8000
C      DATA IMACH( 1) /     0 /
C      DATA IMACH( 2) /     0 /
C      DATA IMACH( 3) /     7 /
C      DATA IMACH( 4) /     0 /
C      DATA IMACH( 5) /    32 /
C      DATA IMACH( 6) /     1 /
C      DATA IMACH( 7) /     2 /
C      DATA IMACH( 8) /    31 /
C      DATA IMACH( 9) /  2147483647 /
C      DATA IMACH(10) /     2 /
C      DATA IMACH(11) /    24 /
C      DATA IMACH(12) /  -125 /
C      DATA IMACH(13) /   128 /
C      DATA IMACH(14) /    53 /
C      DATA IMACH(15) / -1021 /
C      DATA IMACH(16) /  1024 /
C
C     MACHINE CONSTANTS FOR THE UNIVAC 1100 SERIES. FTN COMPILER
C
C === MACHINE = UNIVAC.1100
C      DATA IMACH( 1) /    5 /
C      DATA IMACH( 2) /    6 /
C      DATA IMACH( 3) /    1 /
C      DATA IMACH( 4) /    6 /
C      DATA IMACH( 5) /   36 /
C      DATA IMACH( 6) /    4 /
C      DATA IMACH( 7) /    2 /
C      DATA IMACH( 8) /   35 /
C      DATA IMACH( 9) / O377777777777 /
C      DATA IMACH(10) /    2 /
C      DATA IMACH(11) /   27 /
C      DATA IMACH(12) / -128 /
C      DATA IMACH(13) /  127 /
C      DATA IMACH(14) /   60 /
C      DATA IMACH(15) /-1024 /
C      DATA IMACH(16) / 1023 /
C
C     MACHINE CONSTANTS FOR THE VAX 11/780
C
C === MACHINE = VAX.11/780
C      DATA IMACH(1) /    5 /
C      DATA IMACH(2) /    6 /
C      DATA IMACH(3) /    5 /
C      DATA IMACH(4) /    6 /
C      DATA IMACH(5) /   32 /
C      DATA IMACH(6) /    4 /
C      DATA IMACH(7) /    2 /
C      DATA IMACH(8) /   31 /
C      DATA IMACH(9) /2147483647 /
C      DATA IMACH(10)/    2 /
C      DATA IMACH(11)/   24 /
C      DATA IMACH(12)/ -127 /
C      DATA IMACH(13)/  127 /
C      DATA IMACH(14)/   56 /
C      DATA IMACH(15)/ -127 /
C      DATA IMACH(16)/  127 /
C
C
C***FIRST EXECUTABLE STATEMENT  I1MACH
C
      I1MACH=IMACH(I)
      RETURN
C
      END
      REAL FUNCTION R1MACH(I)
C***BEGIN PROLOGUE  R1MACH
C***DATE WRITTEN   790101   (YYMMDD)
C***REVISION DATE  910131   (YYMMDD)
C***CATEGORY NO.  R1
C***KEYWORDS  MACHINE CONSTANTS
C***AUTHOR  FOX, P. A., (BELL LABS)
C           HALL, A. D., (BELL LABS)
C           SCHRYER, N. L., (BELL LABS)
C***PURPOSE  Returns single precision machine dependent constants
C***DESCRIPTION
C
C     This is the CMLIB version of R1MACH, the real machine
C     constants subroutine originally developed for the PORT library.
C
C     R1MACH can be used to obtain machine-dependent parameters
C     for the local machine environment.  It is a function
C     subroutine with one (input) argument, and can be called
C     as follows, for example
C
C          A = R1MACH(I)
C
C     where I=1,...,5.  The (output) value of A above is
C     determined by the (input) value of I.  The results for
C     various values of I are discussed below.
C
C  Single-Precision Machine Constants
C  R1MACH(1) = B**(EMIN-1), the smallest positive magnitude.
C  R1MACH(2) = B**EMAX*(1 - B**(-T)), the largest magnitude.
C  R1MACH(3) = B**(-T), the smallest relative spacing.
C  R1MACH(4) = B**(1-T), the largest relative spacing.
C  R1MACH(5) = LOG10(B)
C***REFERENCES  FOX, P.A., HALL, A.D., SCHRYER, N.L, *FRAMEWORK FOR
C                 A PORTABLE LIBRARY*, ACM TRANSACTIONS ON MATHE-
C                 MATICAL SOFTWARE, VOL. 4, NO. 2, JUNE 1978,
C                 PP. 177-188.
C***ROUTINES CALLED  (NONE)
C***END PROLOGUE  R1MACH
C
      INTEGER SMALL(2)
      INTEGER LARGE(2)
      INTEGER RIGHT(2)
      INTEGER DIVER(2)
      INTEGER LOG10(2)
C
      REAL RMACH(5)
C
      EQUIVALENCE (RMACH(1),SMALL(1))
      EQUIVALENCE (RMACH(2),LARGE(1))
      EQUIVALENCE (RMACH(3),RIGHT(1))
      EQUIVALENCE (RMACH(4),DIVER(1))
      EQUIVALENCE (RMACH(5),LOG10(1))
C
C
C     MACHINE CONSTANTS FOR IEEE ARITHMETIC MACHINES, SUCH AS THE AT&T
C     3B SERIES, MOTOROLA 68000 BASED MACHINES (E.G. SUN 3 AND AT&T
C     PC 7300), AND 8087 BASED MICROS (E.G. IBM PC AND AT&T 6300).
C
C === MACHINE = IEEE.MOST-SIG-BYTE-FIRST
C === MACHINE = IEEE.LEAST-SIG-BYTE-FIRST
C === MACHINE = SUN
C === MACHINE = 68000
C === MACHINE = 8087
C === MACHINE = IBM.PC
C === MACHINE = ATT.3B
C === MACHINE = ATT.6300
C === MACHINE = ATT.7300
       DATA SMALL(1) /     8388608 /
       DATA LARGE(1) /  2139095039 /
       DATA RIGHT(1) /   864026624 /
       DATA DIVER(1) /   872415232 /
       DATA LOG10(1) /  1050288283 /
C
C     MACHINE CONSTANTS FOR AMDAHL MACHINES.
C
C === MACHINE = AMDAHL
C      DATA SMALL(1) /    1048576 /
C      DATA LARGE(1) / 2147483647 /
C      DATA RIGHT(1) /  990904320 /
C      DATA DIVER(1) / 1007681536 /
C      DATA LOG10(1) / 1091781651 /
C
C     MACHINE CONSTANTS FOR THE BURROUGHS 1700 SYSTEM.
C
C === MACHINE = BURROUGHS.1700
C      DATA RMACH(1) / Z400800000 /
C      DATA RMACH(2) / Z5FFFFFFFF /
C      DATA RMACH(3) / Z4E9800000 /
C      DATA RMACH(4) / Z4EA800000 /
C      DATA RMACH(5) / Z500E730E8 /
C
C     MACHINE CONSTANTS FOR THE BURROUGHS 5700/6700/7700 SYSTEMS.
C
C === MACHINE = BURROUGHS.5700
C === MACHINE = BURROUGHS.6700
C === MACHINE = BURROUGHS.7700
C      DATA RMACH(1) / O1771000000000000 /
C      DATA RMACH(2) / O0777777777777777 /
C      DATA RMACH(3) / O1311000000000000 /
C      DATA RMACH(4) / O1301000000000000 /
C      DATA RMACH(5) / O1157163034761675 /
C
C     MACHINE CONSTANTS FOR THE CONVEX C-120 (NATIVE MODE)
C
C === MACHINE = CONVEX.C1
C      DATA RMACH(1) / 2.9387360E-39 /
C      DATA RMACH(2) / 1.7014117E+38 /
C      DATA RMACH(3) / 5.9604645E-08 /
C      DATA RMACH(4) / 1.1920929E-07 /
C      DATA RMACH(5) / 3.0102999E-01 /
C
C     MACHINE CONSTANTS FOR THE CONVEX C-120 (NATIVE MODE)
C     WITH -R8 OPTION
C
C === MACHINE = CONVEX.C1.R8
C      DATA RMACH(1) / 5.562684646268007D-309 /
C      DATA RMACH(2) / 8.988465674311577D+307 /
C      DATA RMACH(3) / 1.110223024625157D-016 /
C      DATA RMACH(4) / 2.220446049250313D-016 /
C      DATA RMACH(5) / 3.010299956639812D-001 /
C
C     MACHINE CONSTANTS FOR THE CONVEX C-120 (IEEE MODE)
C
C === MACHINE = CONVEX.C1.IEEE
C      DATA RMACH(1) / 1.1754945E-38 /
C      DATA RMACH(2) / 3.4028234E+38 /
C      DATA RMACH(3) / 5.9604645E-08 /
C      DATA RMACH(4) / 1.1920929E-07 /
C      DATA RMACH(5) / 3.0102999E-01 /
C
C     MACHINE CONSTANTS FOR THE CONVEX C-120 (IEEE MODE)
C     WITH -R8 OPTION
C
C === MACHINE = CONVEX.C1.IEEE.R8
C      DATA RMACH(1) / 2.225073858507202D-308 /
C      DATA RMACH(2) / 1.797693134862315D+308 /
C      DATA RMACH(3) / 1.110223024625157D-016 /
C      DATA RMACH(4) / 2.220446049250313D-016 /
C      DATA RMACH(5) / 3.010299956639812D-001 /
C
C     MACHINE CONSTANTS FOR THE CYBER 170/180 SERIES USING NOS (FTN5).
C
C === MACHINE = CYBER.170.NOS
C === MACHINE = CYBER.180.NOS
C      DATA RMACH(1) / O"00014000000000000000" /
C      DATA RMACH(2) / O"37767777777777777777" /
C      DATA RMACH(3) / O"16404000000000000000" /
C      DATA RMACH(4) / O"16414000000000000000" /
C      DATA RMACH(5) / O"17164642023241175720" /
C
C     MACHINE CONSTANTS FOR THE CDC 180 SERIES USING NOS/VE
C
C === MACHINE = CYBER.180.NOS/VE
C      DATA RMACH(1) / Z"3001800000000000" /
C      DATA RMACH(2) / Z"4FFEFFFFFFFFFFFE" /
C      DATA RMACH(3) / Z"3FD2800000000000" /
C      DATA RMACH(4) / Z"3FD3800000000000" /
C      DATA RMACH(5) / Z"3FFF9A209A84FBCF" /
C
C     MACHINE CONSTANTS FOR THE CYBER 205
C
C === MACHINE = CYBER.205
C      DATA RMACH(1) / X'9000400000000000' /
C      DATA RMACH(2) / X'6FFF7FFFFFFFFFFF' /
C      DATA RMACH(3) / X'FFA3400000000000' /
C      DATA RMACH(4) / X'FFA4400000000000' /
C      DATA RMACH(5) / X'FFD04D104D427DE8' /
C
C     MACHINE CONSTANTS FOR THE CDC 6000/7000 SERIES.
C
C === MACHINE = CDC.6000
C === MACHINE = CDC.7000
C      DATA RMACH(1) / 00014000000000000000B /
C      DATA RMACH(2) / 37767777777777777777B /
C      DATA RMACH(3) / 16404000000000000000B /
C      DATA RMACH(4) / 16414000000000000000B /
C      DATA RMACH(5) / 17164642023241175720B /
C
C     MACHINE CONSTANTS FOR THE CRAY 1, XMP, 2, AND 3.
C
C === MACHINE = CRAY.46-BIT-INTEGER
C === MACHINE = CRAY.64-BIT-INTEGER
C      DATA RMACH(1) / 200034000000000000000B /
C      DATA RMACH(2) / 577767777777777777776B /
C      DATA RMACH(3) / 377224000000000000000B /
C      DATA RMACH(4) / 377234000000000000000B /
C      DATA RMACH(5) / 377774642023241175720B /
C
C     MACHINE CONSTANTS FOR THE DATA GENERAL ECLIPSE S/200
C
C     NOTE - IT MAY BE APPROPRIATE TO INCLUDE THE FOLLOWING LINE -
C     STATIC RMACH(5)
C
C === MACHINE = DATA_GENERAL.ECLIPSE.S/200
C      DATA SMALL/20K,0/,LARGE/77777K,177777K/
C      DATA RIGHT/35420K,0/,DIVER/36020K,0/
C      DATA LOG10/40423K,42023K/
C
C     ELXSI 6400
C
C === MACHINE = ELSXI.6400
C      DATA SMALL(1) / '00800000'X /
C      DATA LARGE(1) / '7F7FFFFF'X /
C      DATA RIGHT(1) / '33800000'X /
C      DATA DIVER(1) / '34000000'X /
C      DATA LOG10(1) / '3E9A209B'X /
C
C     MACHINE CONSTANTS FOR THE HARRIS 220
C     MACHINE CONSTANTS FOR THE HARRIS SLASH 6 AND SLASH 7
C
C === MACHINE = HARRIS.220
C === MACHINE = HARRIS.SLASH6
C === MACHINE = HARRIS.SLASH7
C      DATA SMALL(1),SMALL(2) / '20000000, '00000201 /
C      DATA LARGE(1),LARGE(2) / '37777777, '00000177 /
C      DATA RIGHT(1),RIGHT(2) / '20000000, '00000352 /
C      DATA DIVER(1),DIVER(2) / '20000000, '00000353 /
C      DATA LOG10(1),LOG10(2) / '23210115, '00000377 /
C
C     MACHINE CONSTANTS FOR THE HONEYWELL 600/6000 SERIES.
C     MACHINE CONSTANTS FOR THE HONEYWELL DPS 8/70 SERIES.
C
C === MACHINE = HONEYWELL.600/6000
C === MACHINE = HONEYWELL.DPS.8/70
C      DATA RMACH(1) / O402400000000 /
C      DATA RMACH(2) / O376777777777 /
C      DATA RMACH(3) / O714400000000 /
C      DATA RMACH(4) / O716400000000 /
C      DATA RMACH(5) / O776464202324 /
C
C     MACHINE CONSTANTS FOR THE HP 2100
C     3 WORD DOUBLE PRECISION WITH FTN4
C
C === MACHINE = HP.2100.3_WORD_DP
C      DATA SMALL(1), SMALL(2) / 40000B,       1 /
C      DATA LARGE(1), LARGE(2) / 77777B, 177776B /
C      DATA RIGHT(1), RIGHT(2) / 40000B,    325B /
C      DATA DIVER(1), DIVER(2) / 40000B,    327B /
C      DATA LOG10(1), LOG10(2) / 46420B,  46777B /
C
C     MACHINE CONSTANTS FOR THE HP 2100
C     4 WORD DOUBLE PRECISION WITH FTN4
C
C === MACHINE = HP.2100.4_WORD_DP
C      DATA SMALL(1), SMALL(2) / 40000B,       1 /
C      DATA LARGE91), LARGE(2) / 77777B, 177776B /
C      DATA RIGHT(1), RIGHT(2) / 40000B,    325B /
C      DATA DIVER(1), DIVER(2) / 40000B,    327B /
C      DATA LOG10(1), LOG10(2) / 46420B,  46777B /
C
C     HP 9000
C
C      R1MACH(1) = 1.17549435E-38
C      R1MACH(2) = 1.70141163E+38
C      R1MACH(3) = 5.960464478E-8
C      R1MACH(4) = 1.119209290E-7
C      R1MACH(5) = 3.01030010E-1
C
C === MACHINE = HP.9000
C      DATA SMALL(1) / 00040000000B /
C      DATA LARGE(1) / 17677777777B /
C      DATA RIGHT(1) / 06340000000B /
C      DATA DIVER(1) / 06400000000B /
C      DATA LOG10(1) / 07646420233B /
C
C     MACHINE CONSTANTS FOR THE IBM 360/370 SERIES,
C     THE XEROX SIGMA 5/7/9, THE SEL SYSTEMS 85/86 AND
C     THE INTERDATA 3230 AND INTERDATA 7/32.
C
C === MACHINE = IBM.360
C === MACHINE = IBM.370
C === MACHINE = XEROX.SIGMA.5
C === MACHINE = XEROX.SIGMA.7
C === MACHINE = XEROX.SIGMA.9
C === MACHINE = SEL.85
C === MACHINE = SEL.86
C === MACHINE = INTERDATA.3230
C === MACHINE = INTERDATA.7/32
C      DATA RMACH(1) / Z00100000 /
C      DATA RMACH(2) / Z7FFFFFFF /
C      DATA RMACH(3) / Z3B100000 /
C      DATA RMACH(4) / Z3C100000 /
C      DATA RMACH(5) / Z41134413 /
C
C     MACHINE CONSTANTS FOR THE INTERDATA 8/32
C     WITH THE UNIX SYSTEM FORTRAN 77 COMPILER.
C
C     FOR THE INTERDATA FORTRAN VII COMPILER REPLACE
C     THE Z'S SPECIFYING HEX CONSTANTS WITH Y'S.
C
C === MACHINE = INTERDATA.8/32.UNIX
C      DATA RMACH(1) / Z'00100000' /
C      DATA RMACH(2) / Z'7EFFFFFF' /
C      DATA RMACH(3) / Z'3B100000' /
C      DATA RMACH(4) / Z'3C100000' /
C      DATA RMACH(5) / Z'41134413' /
C
C     MACHINE CONSTANTS FOR THE PDP-10 (KA OR KI PROCESSOR).
C
C === MACHINE = PDP-10.KA
C === MACHINE = PDP-10.KI
C      DATA RMACH(1) / "000400000000 /
C      DATA RMACH(2) / "377777777777 /
C      DATA RMACH(3) / "146400000000 /
C      DATA RMACH(4) / "147400000000 /
C      DATA RMACH(5) / "177464202324 /
C
C     MACHINE CONSTANTS FOR PDP-11 FORTRAN SUPPORTING
C     32-BIT INTEGERS (EXPRESSED IN INTEGER AND OCTAL).
C
C === MACHINE = PDP-11.32-BIT
C      DATA SMALL(1) /    8388608 /
C      DATA LARGE(1) / 2147483647 /
C      DATA RIGHT(1) /  880803840 /
C      DATA DIVER(1) /  889192448 /
C      DATA LOG10(1) / 1067065499 /
C
C      DATA RMACH(1) / O00040000000 /
C      DATA RMACH(2) / O17777777777 /
C      DATA RMACH(3) / O06440000000 /
C      DATA RMACH(4) / O06500000000 /
C      DATA RMACH(5) / O07746420233 /
C
C     MACHINE CONSTANTS FOR PDP-11 FORTRAN SUPPORTING
C     16-BIT INTEGERS  (EXPRESSED IN INTEGER AND OCTAL).
C
C === MACHINE = PDP-11.16-BIT
C      DATA SMALL(1),SMALL(2) /   128,     0 /
C      DATA LARGE(1),LARGE(2) / 32767,    -1 /
C      DATA RIGHT(1),RIGHT(2) / 13440,     0 /
C      DATA DIVER(1),DIVER(2) / 13568,     0 /
C      DATA LOG10(1),LOG10(2) / 16282,  8347 /
C
C      DATA SMALL(1),SMALL(2) / O000200, O000000 /
C      DATA LARGE(1),LARGE(2) / O077777, O177777 /
C      DATA RIGHT(1),RIGHT(2) / O032200, O000000 /
C      DATA DIVER(1),DIVER(2) / O032400, O000000 /
C      DATA LOG10(1),LOG10(2) / O037632, O020233 /
C
C     MACHINE CONSTANTS FOR THE SEQUENT BALANCE 8000.
C
C === MACHINE = SEQUENT.BALANCE.8000
C      DATA SMALL(1) / $00800000 /
C      DATA LARGE(1) / $7F7FFFFF /
C      DATA RIGHT(1) / $33800000 /
C      DATA DIVER(1) / $34000000 /
C      DATA LOG10(1) / $3E9A209B /
C
C     MACHINE CONSTANTS FOR THE UNIVAC 1100 SERIES.
C
C === MACHINE = UNIVAC.1100
C      DATA RMACH(1) / O000400000000 /
C      DATA RMACH(2) / O377777777777 /
C      DATA RMACH(3) / O146400000000 /
C      DATA RMACH(4) / O147400000000 /
C      DATA RMACH(5) / O177464202324 /
C
C     MACHINE CONSTANTS FOR THE VAX 11/780
C    (EXPRESSED IN INTEGER AND HEXADECIMAL)
C  *** THE INTEGER FORMAT SHOULD BE OK FOR UNIX SYSTEMS***
C
C === MACHINE = VAX.11/780
C      DATA SMALL(1) /       128 /
C      DATA LARGE(1) /    -32769 /
C      DATA RIGHT(1) /     13440 /
C      DATA DIVER(1) /     13568 /
C      DATA LOG10(1) / 547045274 /
C
C  ***THE HEX FORMAT BELOW MAY NOT BE SUITABLE FOR UNIX SYSTEMS***
C
C      DATA SMALL(1) / Z00000080 /
C      DATA LARGE(1) / ZFFFF7FFF /
C      DATA RIGHT(1) / Z00003480 /
C      DATA DIVER(1) / Z00003500 /
C      DATA LOG10(1) / Z209B3F9A /
C
C
C***FIRST EXECUTABLE STATEMENT  R1MACH
C
      R1MACH = RMACH(I)
      RETURN
C
      END
      DOUBLE PRECISION FUNCTION D1MACH(I)
C***BEGIN PROLOGUE  D1MACH
C***DATE WRITTEN   750101   (YYMMDD)
C***REVISION DATE  910131   (YYMMDD)
C***CATEGORY NO.  R1
C***KEYWORDS  MACHINE CONSTANTS
C***AUTHOR  FOX, P. A., (BELL LABS)
C           HALL, A. D., (BELL LABS)
C           SCHRYER, N. L., (BELL LABS)
C***PURPOSE  Returns double precision machine dependent constants
C***DESCRIPTION
C
C     This is the CMLIB version of D1MACH, the double precision machine
C     constants subroutine originally developed for the PORT library.
C
C     D1MACH can be used to obtain machine-dependent parameters
C     for the local machine environment.  It is a function
C     subprogram with one (input) argument, and can be called
C     as follows, for example
C
C          D = D1MACH(I)
C
C     where I=1,...,5.  The (output) value of D above is
C     determined by the (input) value of I.  The results for
C     various values of I are discussed below.
C
C  Double-precision machine constants
C  D1MACH( 1) = B**(EMIN-1), the smallest positive magnitude.
C  D1MACH( 2) = B**EMAX*(1 - B**(-T)), the largest magnitude.
C  D1MACH( 3) = B**(-T), the smallest relative spacing.
C  D1MACH( 4) = B**(1-T), the largest relative spacing.
C  D1MACH( 5) = LOG10(B)
C***REFERENCES  FOX P.A., HALL A.D., SCHRYER N.L.,*FRAMEWORK FOR A
C                 PORTABLE LIBRARY*, ACM TRANSACTIONS ON MATHEMATICAL
C                 SOFTWARE, VOL. 4, NO. 2, JUNE 1978, PP. 177-188.
C***ROUTINES CALLED  (NONE)
C***END PROLOGUE  D1MACH
C
      INTEGER SMALL(4)
      INTEGER LARGE(4)
      INTEGER RIGHT(4)
      INTEGER DIVER(4)
      INTEGER LOG10(4)
C
      DOUBLE PRECISION DMACH(5)
C
      EQUIVALENCE (DMACH(1),SMALL(1))
      EQUIVALENCE (DMACH(2),LARGE(1))
      EQUIVALENCE (DMACH(3),RIGHT(1))
      EQUIVALENCE (DMACH(4),DIVER(1))
      EQUIVALENCE (DMACH(5),LOG10(1))
C
C     MACHINE CONSTANTS FOR IEEE ARITHMETIC MACHINES, SUCH AS THE AT&T
C     3B SERIES AND MOTOROLA 68000 BASED MACHINES (E.G. SUN 3 AND AT&T
C     PC 7300), IN WHICH THE MOST SIGNIFICANT BYTE IS STORED FIRST.
C
C === MACHINE = IEEE.MOST-SIG-BYTE-FIRST
C === MACHINE = SUN
C === MACHINE = 68000
C === MACHINE = ATT.3B
C === MACHINE = ATT.7300
       DATA SMALL(1),SMALL(2) /    1048576,          0 /
       DATA LARGE(1),LARGE(2) / 2146435071,         -1 /
       DATA RIGHT(1),RIGHT(2) / 1017118720,          0 /
       DATA DIVER(1),DIVER(2) / 1018167296,          0 /
       DATA LOG10(1),LOG10(2) / 1070810131, 1352628735 /
C
C     MACHINE CONSTANTS FOR IEEE ARITHMETIC MACHINES AND 8087-BASED
C     MICROS, SUCH AS THE IBM PC AND AT&T 6300, IN WHICH THE LEAST
C     SIGNIFICANT BYTE IS STORED FIRST.
C
C === MACHINE = IEEE.LEAST-SIG-BYTE-FIRST
C === MACHINE = 8087
C === MACHINE = IBM.PC
C === MACHINE = ATT.6300
C      DATA SMALL(1),SMALL(2) /          0,    1048576 /
C      DATA LARGE(1),LARGE(2) /         -1, 2146435071 /
C      DATA RIGHT(1),RIGHT(2) /          0, 1017118720 /
C      DATA DIVER(1),DIVER(2) /          0, 1018167296 /
C      DATA LOG10(1),LOG10(2) / 1352628735, 1070810131 /
C
C     MACHINE CONSTANTS FOR AMDAHL MACHINES.
C
C === MACHINE = AMDAHL
C      DATA SMALL(1),SMALL(2) /    1048576,          0 /
C      DATA LARGE(1),LARGE(2) / 2147483647,         -1 /
C      DATA RIGHT(1),RIGHT(2) /  856686592,          0 /
C      DATA DIVER(1),DIVER(2) /  873463808,          0 /
C      DATA LOG10(1),LOG10(2) / 1091781651, 1352628735 /
C
C     MACHINE CONSTANTS FOR THE BURROUGHS 1700 SYSTEM.
C
C === MACHINE = BURROUGHS.1700
C      DATA SMALL(1) / ZC00800000 /
C      DATA SMALL(2) / Z000000000 /
C      DATA LARGE(1) / ZDFFFFFFFF /
C      DATA LARGE(2) / ZFFFFFFFFF /
C      DATA RIGHT(1) / ZCC5800000 /
C      DATA RIGHT(2) / Z000000000 /
C      DATA DIVER(1) / ZCC6800000 /
C      DATA DIVER(2) / Z000000000 /
C      DATA LOG10(1) / ZD00E730E7 /
C      DATA LOG10(2) / ZC77800DC0 /
C
C     MACHINE CONSTANTS FOR THE BURROUGHS 5700 SYSTEM.
C
C === MACHINE = BURROUGHS.5700
C      DATA SMALL(1) / O1771000000000000 /
C      DATA SMALL(2) / O0000000000000000 /
C      DATA LARGE(1) / O0777777777777777 /
C      DATA LARGE(2) / O0007777777777777 /
C      DATA RIGHT(1) / O1461000000000000 /
C      DATA RIGHT(2) / O0000000000000000 /
C      DATA DIVER(1) / O1451000000000000 /
C      DATA DIVER(2) / O0000000000000000 /
C      DATA LOG10(1) / O1157163034761674 /
C      DATA LOG10(2) / O0006677466732724 /
C
C     MACHINE CONSTANTS FOR THE BURROUGHS 6700/7700 SYSTEMS.
C
C === MACHINE = BURROUGHS.6700
C === MACHINE = BURROUGHS.7700
C      DATA SMALL(1) / O1771000000000000 /
C      DATA SMALL(2) / O7770000000000000 /
C      DATA LARGE(1) / O0777777777777777 /
C      DATA LARGE(2) / O7777777777777777 /
C      DATA RIGHT(1) / O1461000000000000 /
C      DATA RIGHT(2) / O0000000000000000 /
C      DATA DIVER(1) / O1451000000000000 /
C      DATA DIVER(2) / O0000000000000000 /
C      DATA LOG10(1) / O1157163034761674 /
C      DATA LOG10(2) / O0006677466732724 /
C
C     MACHINE CONSTANTS FOR THE CONVEX C-120 (NATIVE MODE)
C     WITH OR WITHOUT -R8 OPTION
C
C === MACHINE = CONVEX.C1
C === MACHINE = CONVEX.C1.R8
C      DATA DMACH(1) / 5.562684646268007D-309 /
C      DATA DMACH(2) / 8.988465674311577D+307 /
C      DATA DMACH(3) / 1.110223024625157D-016 /
C      DATA DMACH(4) / 2.220446049250313D-016 /
C      DATA DMACH(5) / 3.010299956639812D-001 /
C
C     MACHINE CONSTANTS FOR THE CONVEX C-120 (IEEE MODE)
C     WITH OR WITHOUT -R8 OPTION
C
C === MACHINE = CONVEX.C1.IEEE
C === MACHINE = CONVEX.C1.IEEE.R8
C      DATA DMACH(1) / 2.225073858507202D-308 /
C      DATA DMACH(2) / 1.797693134862315D+308 /
C      DATA DMACH(3) / 1.110223024625157D-016 /
C      DATA DMACH(4) / 2.220446049250313D-016 /
C      DATA DMACH(5) / 3.010299956639812D-001 /
C
C     MACHINE CONSTANTS FOR THE CYBER 170/180 SERIES USING NOS (FTN5).
C
C === MACHINE = CYBER.170.NOS
C === MACHINE = CYBER.180.NOS
C      DATA SMALL(1) / O"00604000000000000000" /
C      DATA SMALL(2) / O"00000000000000000000" /
C      DATA LARGE(1) / O"37767777777777777777" /
C      DATA LARGE(2) / O"37167777777777777777" /
C      DATA RIGHT(1) / O"15604000000000000000" /
C      DATA RIGHT(2) / O"15000000000000000000" /
C      DATA DIVER(1) / O"15614000000000000000" /
C      DATA DIVER(2) / O"15010000000000000000" /
C      DATA LOG10(1) / O"17164642023241175717" /
C      DATA LOG10(2) / O"16367571421742254654" /
C
C     MACHINE CONSTANTS FOR THE CDC 180 SERIES USING NOS/VE
C
C === MACHINE = CYBER.180.NOS/VE
C      DATA SMALL(1) / Z"3001800000000000" /
C      DATA SMALL(2) / Z"3001000000000000" /
C      DATA LARGE(1) / Z"4FFEFFFFFFFFFFFE" /
C      DATA LARGE(2) / Z"4FFE000000000000" /
C      DATA RIGHT(1) / Z"3FD2800000000000" /
C      DATA RIGHT(2) / Z"3FD2000000000000" /
C      DATA DIVER(1) / Z"3FD3800000000000" /
C      DATA DIVER(2) / Z"3FD3000000000000" /
C      DATA LOG10(1) / Z"3FFF9A209A84FBCF" /
C      DATA LOG10(2) / Z"3FFFF7988F8959AC" /
C
C     MACHINE CONSTANTS FOR THE CYBER 205
C
C === MACHINE = CYBER.205
C      DATA SMALL(1) / X'9000400000000000' /
C      DATA SMALL(2) / X'8FD1000000000000' /
C      DATA LARGE(1) / X'6FFF7FFFFFFFFFFF' /
C      DATA LARGE(2) / X'6FD07FFFFFFFFFFF' /
C      DATA RIGHT(1) / X'FF74400000000000' /
C      DATA RIGHT(2) / X'FF45000000000000' /
C      DATA DIVER(1) / X'FF75400000000000' /
C      DATA DIVER(2) / X'FF46000000000000' /
C      DATA LOG10(1) / X'FFD04D104D427DE7' /
C      DATA LOG10(2) / X'FFA17DE623E2566A' /
C
C     MACHINE CONSTANTS FOR THE CDC 6000/7000 SERIES.
C
C === MACHINE = CDC.6000
C === MACHINE = CDC.7000
C      DATA SMALL(1) / 00604000000000000000B /
C      DATA SMALL(2) / 00000000000000000000B /
C      DATA LARGE(1) / 37767777777777777777B /
C      DATA LARGE(2) / 37167777777777777777B /
C      DATA RIGHT(1) / 15604000000000000000B /
C      DATA RIGHT(2) / 15000000000000000000B /
C      DATA DIVER(1) / 15614000000000000000B /
C      DATA DIVER(2) / 15010000000000000000B /
C      DATA LOG10(1) / 17164642023241175717B /
C      DATA LOG10(2) / 16367571421742254654B /
C
C     MACHINE CONSTANTS FOR THE CRAY 1, XMP, 2, AND 3.
C
C === MACHINE = CRAY.46-BIT-INTEGER
C === MACHINE = CRAY.64-BIT-INTEGER
C      DATA SMALL(1) / 201354000000000000000B /
C      DATA SMALL(2) / 000000000000000000000B /
C      DATA LARGE(1) / 577767777777777777777B /
C      DATA LARGE(2) / 000007777777777777776B /
C      DATA RIGHT(1) / 376434000000000000000B /
C      DATA RIGHT(2) / 000000000000000000000B /
C      DATA DIVER(1) / 376444000000000000000B /
C      DATA DIVER(2) / 000000000000000000000B /
C      DATA LOG10(1) / 377774642023241175717B /
C      DATA LOG10(2) / 000007571421742254654B /
C
C     MACHINE CONSTANTS FOR THE DATA GENERAL ECLIPSE S/200
C
C     NOTE - IT MAY BE APPROPRIATE TO INCLUDE THE FOLLOWING LINE -
C     STATIC DMACH(5)
C
C === MACHINE = DATA_GENERAL.ECLIPSE.S/200
C      DATA SMALL/20K,3*0/,LARGE/77777K,3*177777K/
C      DATA RIGHT/31420K,3*0/,DIVER/32020K,3*0/
C      DATA LOG10/40423K,42023K,50237K,74776K/
C
C     ELXSI 6400
C
C === MACHINE = ELSXI.6400
C      DATA SMALL(1), SMALL(2) / '00100000'X,'00000000'X /
C      DATA LARGE(1), LARGE(2) / '7FEFFFFF'X,'FFFFFFFF'X /
C      DATA RIGHT(1), RIGHT(2) / '3CB00000'X,'00000000'X /
C      DATA DIVER(1), DIVER(2) / '3CC00000'X,'00000000'X /
C      DATA LOG10(1), DIVER(2) / '3FD34413'X,'509F79FF'X /
C
C     MACHINE CONSTANTS FOR THE HARRIS 220
C     MACHINE CONSTANTS FOR THE HARRIS SLASH 6 AND SLASH 7
C
C === MACHINE = HARRIS.220
C === MACHINE = HARRIS.SLASH6
C === MACHINE = HARRIS.SLASH7
C      DATA SMALL(1),SMALL(2) / '20000000, '00000201 /
C      DATA LARGE(1),LARGE(2) / '37777777, '37777577 /
C      DATA RIGHT(1),RIGHT(2) / '20000000, '00000333 /
C      DATA DIVER(1),DIVER(2) / '20000000, '00000334 /
C      DATA LOG10(1),LOG10(2) / '23210115, '10237777 /
C
C     MACHINE CONSTANTS FOR THE HONEYWELL 600/6000 SERIES.
C     MACHINE CONSTANTS FOR THE HONEYWELL DPS 8/70 SERIES.
C
C === MACHINE = HONEYWELL.600/6000
C === MACHINE = HONEYWELL.DPS.8/70
C      DATA SMALL(1),SMALL(2) / O402400000000, O000000000000 /
C      DATA LARGE(1),LARGE(2) / O376777777777, O777777777777 /
C      DATA RIGHT(1),RIGHT(2) / O604400000000, O000000000000 /
C      DATA DIVER(1),DIVER(2) / O606400000000, O000000000000 /
C      DATA LOG10(1),LOG10(2) / O776464202324, O117571775714 /
C
C      MACHINE CONSTANTS FOR THE HP 2100
C      3 WORD DOUBLE PRECISION OPTION WITH FTN4
C
C === MACHINE = HP.2100.3_WORD_DP
C      DATA SMALL(1), SMALL(2), SMALL(3) / 40000B,       0,       1 /
C      DATA LARGE(1), LARGE(2), LARGE(3) / 77777B, 177777B, 177776B /
C      DATA RIGHT(1), RIGHT(2), RIGHT(3) / 40000B,       0,    265B /
C      DATA DIVER(1), DIVER(2), DIVER(3) / 40000B,       0,    276B /
C      DATA LOG10(1), LOG10(2), LOG10(3) / 46420B,  46502B,  77777B /
C
C      MACHINE CONSTANTS FOR THE HP 2100
C      4 WORD DOUBLE PRECISION OPTION WITH FTN4
C
C === MACHINE = HP.2100.4_WORD_DP
C      DATA SMALL(1), SMALL(2) /  40000B,       0 /
C      DATA SMALL(3), SMALL(4) /       0,       1 /
C      DATA LARGE(1), LARGE(2) /  77777B, 177777B /
C      DATA LARGE(3), LARGE(4) / 177777B, 177776B /
C      DATA RIGHT(1), RIGHT(2) /  40000B,       0 /
C      DATA RIGHT(3), RIGHT(4) /       0,    225B /
C      DATA DIVER(1), DIVER(2) /  40000B,       0 /
C      DATA DIVER(3), DIVER(4) /       0,    227B /
C      DATA LOG10(1), LOG10(2) /  46420B,  46502B /
C      DATA LOG10(3), LOG10(4) /  76747B, 176377B /
C
C     HP 9000
C
C      D1MACH(1) = 2.8480954D-306
C      D1MACH(2) = 1.40444776D+306
C      D1MACH(3) = 2.22044605D-16
C      D1MACH(4) = 4.44089210D-16
C      D1MACH(5) = 3.01029996D-1
C
C === MACHINE = HP.9000
C      DATA SMALL(1), SMALL(2) / 00040000000B, 00000000000B /
C      DATA LARGE(1), LARGE(2) / 17737777777B, 37777777777B /
C      DATA RIGHT(1), RIGHT(2) / 07454000000B, 00000000000B /
C      DATA DIVER(1), DIVER(2) / 07460000000B, 00000000000B /
C      DATA LOG10(1), LOG10(2) / 07764642023B, 12047674777B /
C
C     MACHINE CONSTANTS FOR THE IBM 360/370 SERIES,
C     THE XEROX SIGMA 5/7/9, THE SEL SYSTEMS 85/86, AND
C     THE INTERDATA 3230 AND INTERDATA 7/32.
C
C === MACHINE = IBM.360
C === MACHINE = IBM.370
C === MACHINE = XEROX.SIGMA.5
C === MACHINE = XEROX.SIGMA.7
C === MACHINE = XEROX.SIGMA.9
C === MACHINE = SEL.85
C === MACHINE = SEL.86
C === MACHINE = INTERDATA.3230
C === MACHINE = INTERDATA.7/32
C      DATA SMALL(1),SMALL(2) / Z00100000, Z00000000 /
C      DATA LARGE(1),LARGE(2) / Z7FFFFFFF, ZFFFFFFFF /
C      DATA RIGHT(1),RIGHT(2) / Z33100000, Z00000000 /
C      DATA DIVER(1),DIVER(2) / Z34100000, Z00000000 /
C      DATA LOG10(1),LOG10(2) / Z41134413, Z509F79FF /
C
C     MACHINE CONSTANTS FOR THE INTERDATA 8/32
C     WITH THE UNIX SYSTEM FORTRAN 77 COMPILER.
C
C     FOR THE INTERDATA FORTRAN VII COMPILER REPLACE
C     THE Z'S SPECIFYING HEX CONSTANTS WITH Y'S.
C
C === MACHINE = INTERDATA.8/32.UNIX
C      DATA SMALL(1),SMALL(2) / Z'00100000', Z'00000000' /
C      DATA LARGE(1),LARGE(2) / Z'7EFFFFFF', Z'FFFFFFFF' /
C      DATA RIGHT(1),RIGHT(2) / Z'33100000', Z'00000000' /
C      DATA DIVER(1),DIVER(2) / Z'34100000', Z'00000000' /
C      DATA LOG10(1),LOG10(2) / Z'41134413', Z'509F79FF' /
C
C     MACHINE CONSTANTS FOR THE PDP-10 (KA PROCESSOR).
C
C === MACHINE = PDP-10.KA
C      DATA SMALL(1),SMALL(2) / "033400000000, "000000000000 /
C      DATA LARGE(1),LARGE(2) / "377777777777, "344777777777 /
C      DATA RIGHT(1),RIGHT(2) / "113400000000, "000000000000 /
C      DATA DIVER(1),DIVER(2) / "114400000000, "000000000000 /
C      DATA LOG10(1),LOG10(2) / "177464202324, "144117571776 /
C
C     MACHINE CONSTANTS FOR THE PDP-10 (KI PROCESSOR).
C
C === MACHINE = PDP-10.KI
C      DATA SMALL(1),SMALL(2) / "000400000000, "000000000000 /
C      DATA LARGE(1),LARGE(2) / "377777777777, "377777777777 /
C      DATA RIGHT(1),RIGHT(2) / "103400000000, "000000000000 /
C      DATA DIVER(1),DIVER(2) / "104400000000, "000000000000 /
C      DATA LOG10(1),LOG10(2) / "177464202324, "047674776746 /
C
C     MACHINE CONSTANTS FOR PDP-11 FORTRAN SUPPORTING
C     32-BIT INTEGERS (EXPRESSED IN INTEGER AND OCTAL).
C
C === MACHINE = PDP-11.32-BIT
C      DATA SMALL(1),SMALL(2) /    8388608,           0 /
C      DATA LARGE(1),LARGE(2) / 2147483647,          -1 /
C      DATA RIGHT(1),RIGHT(2) /  612368384,           0 /
C      DATA DIVER(1),DIVER(2) /  620756992,           0 /
C      DATA LOG10(1),LOG10(2) / 1067065498, -2063872008 /
C
C      DATA SMALL(1),SMALL(2) / O00040000000, O00000000000 /
C      DATA LARGE(1),LARGE(2) / O17777777777, O37777777777 /
C      DATA RIGHT(1),RIGHT(2) / O04440000000, O00000000000 /
C      DATA DIVER(1),DIVER(2) / O04500000000, O00000000000 /
C      DATA LOG10(1),LOG10(2) / O07746420232, O20476747770 /
C
C     MACHINE CONSTANTS FOR PDP-11 FORTRAN SUPPORTING
C     16-BIT INTEGERS (EXPRESSED IN INTEGER AND OCTAL).
C
C === MACHINE = PDP-11.16-BIT
C      DATA SMALL(1),SMALL(2) /    128,      0 /
C      DATA SMALL(3),SMALL(4) /      0,      0 /
C      DATA LARGE(1),LARGE(2) /  32767,     -1 /
C      DATA LARGE(3),LARGE(4) /     -1,     -1 /
C      DATA RIGHT(1),RIGHT(2) /   9344,      0 /
C      DATA RIGHT(3),RIGHT(4) /      0,      0 /
C      DATA DIVER(1),DIVER(2) /   9472,      0 /
C      DATA DIVER(3),DIVER(4) /      0,      0 /
C      DATA LOG10(1),LOG10(2) /  16282,   8346 /
C      DATA LOG10(3),LOG10(4) / -31493, -12296 /
C
C      DATA SMALL(1),SMALL(2) / O000200, O000000 /
C      DATA SMALL(3),SMALL(4) / O000000, O000000 /
C      DATA LARGE(1),LARGE(2) / O077777, O177777 /
C      DATA LARGE(3),LARGE(4) / O177777, O177777 /
C      DATA RIGHT(1),RIGHT(2) / O022200, O000000 /
C      DATA RIGHT(3),RIGHT(4) / O000000, O000000 /
C      DATA DIVER(1),DIVER(2) / O022400, O000000 /
C      DATA DIVER(3),DIVER(4) / O000000, O000000 /
C      DATA LOG10(1),LOG10(2) / O037632, O020232 /
C      DATA LOG10(3),LOG10(4) / O102373, O147770 /
C
C     MACHINE CONSTANTS FOR THE SEQUENT BALANCE 8000
C
C === MACHINE = SEQUENT.BALANCE.8000
C      DATA SMALL(1),SMALL(2) / $00000000,  $00100000 /
C      DATA LARGE(1),LARGE(2) / $FFFFFFFF,  $7FEFFFFF /
C      DATA RIGHT(1),RIGHT(2) / $00000000,  $3CA00000 /
C      DATA DIVER(1),DIVER(2) / $00000000,  $3CB00000 /
C      DATA LOG10(1),LOG10(2) / $509F79FF,  $3FD34413 /
C
C     MACHINE CONSTANTS FOR THE UNIVAC 1100 SERIES. FTN COMPILER
C
C === MACHINE = UNIVAC.1100
C      DATA SMALL(1),SMALL(2) / O000040000000, O000000000000 /
C      DATA LARGE(1),LARGE(2) / O377777777777, O777777777777 /
C      DATA RIGHT(1),RIGHT(2) / O170540000000, O000000000000 /
C      DATA DIVER(1),DIVER(2) / O170640000000, O000000000000 /
C      DATA LOG10(1),LOG10(2) / O177746420232, O411757177572 /
C
C     MACHINE CONSTANTS FOR VAX 11/780
C     (EXPRESSED IN INTEGER AND HEXADECIMAL)
C    *** THE INTEGER FORMAT SHOULD BE OK FOR UNIX SYSTEMS***
C
C === MACHINE = VAX.11/780
C      DATA SMALL(1), SMALL(2) /        128,           0 /
C      DATA LARGE(1), LARGE(2) /     -32769,          -1 /
C      DATA RIGHT(1), RIGHT(2) /       9344,           0 /
C      DATA DIVER(1), DIVER(2) /       9472,           0 /
C      DATA LOG10(1), LOG10(2) /  546979738,  -805796613 /
C
C    ***THE HEX FORMAT BELOW MAY NOT BE SUITABLE FOR UNIX SYSYEMS***
C      DATA SMALL(1), SMALL(2) / Z00000080, Z00000000 /
C      DATA LARGE(1), LARGE(2) / ZFFFF7FFF, ZFFFFFFFF /
C      DATA RIGHT(1), RIGHT(2) / Z00002480, Z00000000 /
C      DATA DIVER(1), DIVER(2) / Z00002500, Z00000000 /
C      DATA LOG10(1), LOG10(2) / Z209A3F9A, ZCFF884FB /
C
C   MACHINE CONSTANTS FOR VAX 11/780 (G-FLOATING)
C     (EXPRESSED IN INTEGER AND HEXADECIMAL)
C    *** THE INTEGER FORMAT SHOULD BE OK FOR UNIX SYSTEMS***
C
C      DATA SMALL(1), SMALL(2) /         16,           0 /
C      DATA LARGE(1), LARGE(2) /     -32769,          -1 /
C      DATA RIGHT(1), RIGHT(2) /      15552,           0 /
C      DATA DIVER(1), DIVER(2) /      15568,           0 /
C      DATA LOG10(1), LOG10(2) /  1142112243, 2046775455 /
C
C    ***THE HEX FORMAT BELOW MAY NOT BE SUITABLE FOR UNIX SYSYEMS***
C      DATA SMALL(1), SMALL(2) / Z00000010, Z00000000 /
C      DATA LARGE(1), LARGE(2) / ZFFFF7FFF, ZFFFFFFFF /
C      DATA RIGHT(1), RIGHT(2) / Z00003CC0, Z00000000 /
C      DATA DIVER(1), DIVER(2) / Z00003CD0, Z00000000 /
C      DATA LOG10(1), LOG10(2) / Z44133FF3, Z79FF509F /
C
C
C***FIRST EXECUTABLE STATEMENT  D1MACH
C
      D1MACH = DMACH(I)
      RETURN
C
      END
      SUBROUTINE WFERR (SNAME, MESSAG, LEVEL)
C***BEGIN PROLOGUE  WFERR
C***PURPOSE  Processes an error message for the vectorized special
C            function package
C***LIBRARY   VFNLIB
C***TYPE      CHARACTER
C***AUTHOR  SAUNDERS, B. V., (NIST)
C             COMPUTING AND APPLIED MATHEMATICS LABORATORY
C             NATIONAL INSTITUTE OF STANDARDS AND TECHNOLOGY
C             GAITHERSBURG, MD 20899
C           BOISVERT, R. F., (NIST)
C             COMPUTING AND APPLIED MATHEMATICS LABORATORY
C             NATIONAL INSTITUTE OF STANDARDS AND TECHNOLOGY
C             GAITHERSBURG, MD 20899
C***DESCRIPTION
C
C   WFERR processes an error message for the vectorized special
C   function package.
C
C   P A R A M E T E R S
C
C   SNAME  (Input) Character*(*)
C          The name of the routine issuing the message.
C
C   MESSAG (Input) Character*(*)
C          The error message.  Must be less than 80 characters.
C
C   LEVEL  (Input) Integer
C          The message level.  1 for recoverable, 2 for fatal.
C
C   This routine issues a STOP for fatal errors.
C
C
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   910226  DATE WRITTEN
C***END PROLOGUE  WFERR
C
C     ... PARAMETERS
C
      CHARACTER*(*)  MESSAG, SNAME
      INTEGER        LEVEL
C
C     ... LOCAL VARIABLES
C
      INTEGER IUNIT, I1MACH, MLEN
C
C***FIRST EXECUTABLE STATEMENT  WFERR
C
      IUNIT = I1MACH(4)
      MLEN = MIN(LEN(MESSAG),80)
C
      WRITE(IUNIT,*)
      IF (LEVEL .EQ. 2) THEN
         WRITE(IUNIT,*) 'FATAL ERROR IN ',SNAME,' ...'
         WRITE(IUNIT,*) MESSAG(1:MLEN)
         STOP
      ELSE
         WRITE(IUNIT,*) 'RECOVERABLE ERROR IN ',SNAME,' ...'
         WRITE(IUNIT,*) MESSAG(1:MLEN)
         RETURN
      ENDIF
C
      END
      REAL FUNCTION SECOND ()
C
C  OBTAIN ELAPSED CPU TIME SINCE START OF RUN.
C
C  ON CRAY -- THIS ROUTINE NOT NEEDED.
C
C  ON CONVEX -- MORE ACCURATE TIMINGS ARE OBTAINED BY CALLING
C               THE VECLIB ROUTINE CPUTIME AS FOLLOWS
C
C     SECOND = CPUTIME(0.0E0)
C
C  ON UNIX -- (SUN, FOR EXAMPLE) THE FOLLOWING MAY BE USED
C
      REAL TARRAY(2), ETIME, TOTAL
      TOTAL = ETIME(TARRAY)
      SECOND = TARRAY(1) + TARRAY(2)
C
      RETURN
      END
      SUBROUTINE VI0 (M, X, F, WORK, IWORK, INFO)
C***BEGIN PROLOGUE  VI0
C***PURPOSE  Computes the hyperbolic Bessel function of the first kind
C            of order zero (I0) for a vector of real arguments
C***LIBRARY   VFNLIB
C***CATEGORY  C10B1
C***TYPE      SINGLE PRECISION (VI0-S, DVI0-D)
C***KEYWORDS  BESSEL FUNCTION,FIRST KIND,HYPERBOLIC BESSEL FUNCTION,
C             MODIFIED BESSEL FUNCTION, ORDER ZERO, VECTORIZED
C***AUTHOR  SAUNDERS, B. V., (NIST)
C             COMPUTING AND APPLIED MATHEMATICS LABORATORY
C             NATIONAL INSTITUTE OF STANDARDS AND TECHNOLOGY
C             GAITHERSBURG, MD 20899
C           BOISVERT, R. F., (NIST)
C             COMPUTING AND APPLIED MATHEMATICS LABORATORY
C             NATIONAL INSTITUTE OF STANDARDS AND TECHNOLOGY
C             GAITHERSBURG, MD 20899
C***DESCRIPTION
C
C   VI0 computes the modified (hyperbolic) Bessel function of the
C   first kind of order zero (I0) for real arguments using uniform
C   approximation by Chebyshev polynomials.
C
C
C   P A R A M E T E R S
C
C   M      (Input) Integer (M .GT. 0)
C          The number of arguments at which the function is to be
C          evaluated.
C
C   X      (Input) Real array of length M
C          The arguments at which the function is to be evaluated are
C          stored in X(1) to X(M) in any order.
C
C   F      (Output) Real array of length M
C          F(i) contains the value of the function at X(i), i=1,..,M.
C
C   WORK   (Work) Real vector of length 7*M
C
C   IWORK  (Work) Integer vector of length M
C
C   INFO   (Output) Integer
C          Indicates status of computed result. The following table
C          lists possible values and their meanings.  If OK=Yes then
C          all F(i) have been set, otherwise none have been set.
C
C          INFO  OK            Description
C          ------------------------------------------------------------
C            0   Yes  Successfull execution.
C            1   No   Error: M .LE. 0
C            2   No   Error: Some abs(X(i)) so big I0 overflows.
C                     The index of the first offending argument is
C                     returned in IWORK(1).
C
C
C *********************************************************************
C   This routine is a modification of the function BESI0 developed by
C   W. Fullerton of LANL.
C *********************************************************************
C
C***REFERENCES  (NONE)
C***ROUTINES CALLED  WI0
C***REVISION HISTORY  (YYMMDD)
C   910226  DATE WRITTEN
C***END PROLOGUE  VI0
C
C     ----------
C     PARAMETERS
C     ----------
C
      INTEGER INFO, IWORK, M
      REAL    F, X, WORK
C
      DIMENSION X(M), F(M), WORK(7*M), IWORK(M)
C
C     ---------------
C     LOCAL VARIABLES
C     ---------------
C
      INTEGER IWB0, IWB1, IWB2, IWY, IWTC, IWYC, IWZC, JIN
C
C***FIRST EXECUTABLE STATEMENT  VI0
C
C     ... PARTITION WORK ARRAYS
C
      IWY   = 1
      IWTC  = IWY  + M
      IWYC  = IWTC + M
      IWZC  = IWYC + M
      IWB0  = IWZC + M
      IWB1  = IWB0 + M
      IWB2  = IWB1 + M
C     Total = IWB2 + M
C
      JIN   = 1
C     Total = JIN  + M
C
C     ... WI0 DOES ALL THE WORK
C
      CALL WI0(M,X,F,WORK(IWY),WORK(IWTC),WORK(IWYC),WORK(IWZC),
     +         IWORK(JIN),WORK(IWB0),WORK(IWB1),WORK(IWB2),INFO)
C
      RETURN
      END
      SUBROUTINE WI0 (M, X, F, Y, TCMP, YCMP, ZCMP, INDX, B0, B1, B2, 
     +   INFO)
C***BEGIN PROLOGUE  WI0
C***SUBSIDIARY
C***PURPOSE  Computes the hyperbolic Bessel function of the first kind
C            of order zero (I0) for a vector of arguments
C***LIBRARY   VFNLIB
C***AUTHOR  SAUNDERS, B. V., (NIST)
C             COMPUTING AND APPLIED MATHEMATICS LABORATORY
C             NATIONAL INSTITUTE OF STANDARDS AND TECHNOLOGY
C             GAITHERSBURG, MD 20899
C           BOISVERT, R. F., (NIST)
C             COMPUTING AND APPLIED MATHEMATICS LABORATORY
C             NATIONAL INSTITUTE OF STANDARDS AND TECHNOLOGY
C             GAITHERSBURG, MD 20899
C***ROUTINES CALLED  IWCS, R1MACH, WNGT, WGTHR, WGTLE, WGT, WLE,
C                    WSCTR, WCS
C***REVISION HISTORY  (YYMMDD)
C   910226  DATE WRITTEN
C***END PROLOGUE  WI0
C
C  ----------
C  PARAMETERS
C  ----------
C
      INTEGER INFO, INDX, M
      REAL    B0, B1, B2, F, X, Y, TCMP, YCMP, ZCMP
C
      DIMENSION B0(M), B1(M), B2(M), F(M), INDX(M), X(M), Y(M),
     +          TCMP(M), YCMP(M), ZCMP(M)
C
C  ---------------
C  LOCAL VARIABLES
C  ---------------
C
      INTEGER LAI0, LAI02, LBI0
      PARAMETER ( LAI0=21, LAI02=22, LBI0=12 )
C
      INTEGER I, IWCS, J, KEY, N, NTI0, NTAI0, NTAI02
      REAL    AI0CS, AI02CS, BI0CS, EPMACH, EPS, R1MACH, XSML,
     +        XMAX
C
      DIMENSION AI0CS(LAI0), AI02CS(LAI02), BI0CS(LBI0)
C
      SAVE AI0CS, AI02CS, BI0CS, N, NTAI0, NTAI02, NTI0, XSML, XMAX
C
C----------------------------------------------------------------------
C
C Series for BI0        on the interval  0.          to  9.00000D+00
C                                        with weighted error   2.46E-18
C                                         log weighted error  17.61
C                               significant figures required  17.90
C                                    decimal places required  18.15
C
      DATA BI0 CS( 1) /   -.0766054725 2839144951E0 /
      DATA BI0 CS( 2) /   1.9273379539 93808270E0 /
      DATA BI0 CS( 3) /    .2282644586 920301339E0 /
      DATA BI0 CS( 4) /    .0130489146 6707290428E0 /
      DATA BI0 CS( 5) /    .0004344270 9008164874E0 /
      DATA BI0 CS( 6) /    .0000094226 5768600193E0 /
      DATA BI0 CS( 7) /    .0000001434 0062895106E0 /
      DATA BI0 CS( 8) /    .0000000016 1384906966E0 /
      DATA BI0 CS( 9) /    .0000000000 1396650044E0 /
      DATA BI0 CS(10) /    .0000000000 0009579451E0 /
      DATA BI0 CS(11) /    .0000000000 0000053339E0 /
      DATA BI0 CS(12) /    .0000000000 0000000245E0 /
C
C----------------------------------------------------------------------
C
C Series for AI0        on the interval  1.25000D-01 to  3.33333D-01
C                                        with weighted error   7.87E-17
C                                         log weighted error  16.10
C                               significant figures required  14.69
C                                    decimal places required  16.76
C
      DATA AI0 CS( 1) /    .0757599449 4023796E0 /
      DATA AI0 CS( 2) /    .0075913808 1082334E0 /
      DATA AI0 CS( 3) /    .0004153131 3389237E0 /
      DATA AI0 CS( 4) /    .0000107007 6463439E0 /
      DATA AI0 CS( 5) /   -.0000079011 7997921E0 /
      DATA AI0 CS( 6) /   -.0000007826 1435014E0 /
      DATA AI0 CS( 7) /    .0000002783 8499429E0 /
      DATA AI0 CS( 8) /    .0000000082 5247260E0 /
      DATA AI0 CS( 9) /   -.0000000120 4463945E0 /
      DATA AI0 CS(10) /    .0000000015 5964859E0 /
      DATA AI0 CS(11) /    .0000000002 2925563E0 /
      DATA AI0 CS(12) /   -.0000000001 1916228E0 /
      DATA AI0 CS(13) /    .0000000000 1757854E0 /
      DATA AI0 CS(14) /    .0000000000 0112822E0 /
      DATA AI0 CS(15) /   -.0000000000 0114684E0 /
      DATA AI0 CS(16) /    .0000000000 0027155E0 /
      DATA AI0 CS(17) /   -.0000000000 0002415E0 /
      DATA AI0 CS(18) /   -.0000000000 0000608E0 /
      DATA AI0 CS(19) /    .0000000000 0000314E0 /
      DATA AI0 CS(20) /   -.0000000000 0000071E0 /
      DATA AI0 CS(21) /    .0000000000 0000007E0 /
C
C----------------------------------------------------------------------
C
C Series for AI02       on the interval  0.          to  1.25000D-01
C                                        with weighted error   3.79E-17
C                                         log weighted error  16.42
C                               significant figures required  14.86
C                                    decimal places required  17.09
C
      DATA AI02CS( 1) /    .0544904110 1410882E0 /
      DATA AI02CS( 2) /    .0033691164 7825569E0 /
      DATA AI02CS( 3) /    .0000688975 8346918E0 /
      DATA AI02CS( 4) /    .0000028913 7052082E0 /
      DATA AI02CS( 5) /    .0000002048 9185893E0 /
      DATA AI02CS( 6) /    .0000000226 6668991E0 /
      DATA AI02CS( 7) /    .0000000033 9623203E0 /
      DATA AI02CS( 8) /    .0000000004 9406022E0 /
      DATA AI02CS( 9) /    .0000000000 1188914E0 /
      DATA AI02CS(10) /   -.0000000000 3149915E0 /
      DATA AI02CS(11) /   -.0000000000 1321580E0 /
      DATA AI02CS(12) /   -.0000000000 0179419E0 /
      DATA AI02CS(13) /    .0000000000 0071801E0 /
      DATA AI02CS(14) /    .0000000000 0038529E0 /
      DATA AI02CS(15) /    .0000000000 0001539E0 /
      DATA AI02CS(16) /   -.0000000000 0004151E0 /
      DATA AI02CS(17) /   -.0000000000 0000954E0 /
      DATA AI02CS(18) /    .0000000000 0000382E0 /
      DATA AI02CS(19) /    .0000000000 0000176E0 /
      DATA AI02CS(20) /   -.0000000000 0000034E0 /
      DATA AI02CS(21) /   -.0000000000 0000027E0 /
      DATA AI02CS(22) /    .0000000000 0000003E0 /
C
C----------------------------------------------------------------------
C
      DATA NTI0 / 0 /
C
C***FIRST EXECUTABLE STATEMENT  WI0
C
      IF (M .LE. 0)  GO TO 910
C
      IF (NTI0.EQ.0) THEN
         EPMACH = R1MACH(3)
         EPS = 0.10E0*EPMACH
         NTI0 = IWCS(BI0CS, LBI0, EPS)
         NTAI0 = IWCS(AI0CS, LAI0, EPS)
         NTAI02 = IWCS(AI02CS, LAI02, EPS)
         XSML = SQRT(4.0E0*EPMACH)
         XMAX = LOG(R1MACH(2))
      ENDIF
C
      DO 10 I=1,M
         Y(I) = ABS(X(I))
  10  CONTINUE
C
      CALL WNGT(M,Y,XMAX,KEY)
      IF (KEY .NE. 0)  GO TO 920
C
C  ----------------
C  CASE Y .LE. XSML
C  ----------------
C
      DO 15 I=1,M
         F(I) = 1.0E0
  15  CONTINUE
C
C  --------------------------
C  CASE  XSML .LT. Y .LE. 3.0
C  --------------------------
C
      CALL WGTLE(M,Y,XSML,3.0E0,N,INDX)
      IF (N .GT. 0) THEN
         CALL WGTHR(N,Y,INDX,YCMP)
         DO 20 J=1,N
            TCMP(J) = YCMP(J)**2/4.50E0 - 1.0E0
  20     CONTINUE
         CALL WCS(N,TCMP,BI0CS,NTI0,ZCMP,B0,B1,B2)
         DO 30 J=1,N
            ZCMP(J) = 2.750E0 + ZCMP(J)
  30     CONTINUE
         CALL WSCTR(N,ZCMP,INDX,F)
      ENDIF
C
C  -------------------------
C  CASE  3.0 .LT. Y .LE. 8.0
C  -------------------------
C
      CALL WGTLE(M,Y,3.0E0,8.0E0,N,INDX)
      IF (N .GT. 0) THEN
         CALL WGTHR(N,Y,INDX,YCMP)
         DO 50 J=1,N
            TCMP(J) = (48.0E0/YCMP(J) - 11.0E0)/5.0E0
  50     CONTINUE
         CALL WCS(N,TCMP,AI0CS,NTAI0,ZCMP,B0,B1,B2)
         DO 60 J=1,N
            ZCMP(J) = EXP(YCMP(J))*(0.3750E0+ZCMP(J))/SQRT(YCMP(J))
  60     CONTINUE
         CALL WSCTR(N,ZCMP,INDX,F)
      ENDIF
C
C  ----------------
C  CASE  Y .GT. 8.0
C  ----------------
C
      CALL WGT(M,Y,8.0E0,N,INDX)
      IF (N .GT. 0) THEN
         CALL WGTHR(N,Y,INDX,YCMP)
        DO 80 J=1,N
            TCMP(J) = 16.0E0/YCMP(J) - 1.0E0
  80     CONTINUE
         CALL WCS(N,TCMP,AI02CS,NTAI02,ZCMP,B0,B1,B2)
         DO 90 J=1,N
            ZCMP(J) = EXP(YCMP(J))*(0.3750E0+ZCMP(J))/SQRT(YCMP(J))
  90     CONTINUE
         CALL WSCTR(N,ZCMP,INDX,F)
      ENDIF
C
C  -----
C  EXITS
C  -----
C
C     ... NORMAL
C
      INFO = 0
      GO TO 999
C
C     ... M .LE. 0
C
  910 CONTINUE
      INFO = 1
      GO TO 999
C
C     ... ABS(X) SO LARGE I0 OVERFLOWS
C
  920 CONTINUE
      INFO = 2
      INDX(1) = KEY
      GO TO 999
C
  999 CONTINUE
      RETURN
      END

       SUBROUTINE VI1 (M, X, F, WORK, IWORK, INFO)
C***BEGIN PROLOGUE  VI1
C***PURPOSE  Computes the hyperbolic Bessel function of the first kind
C            of order one (I1) for a vector of real arguments
C***LIBRARY   VFNLIB
C***CATEGORY  C10B1
C***TYPE      SINGLE PRECISION (VI1-S, DVI1-D)
C***KEYWORDS  BESSEL FUNCTION,FIRST KIND,HYPERBOLIC BESSEL FUNCTION,
C             MODIFIED BESSEL FUNCTION, ORDER ONE, VECTORIZED
C***AUTHOR  SAUNDERS, B. V., (NIST)
C             COMPUTING AND APPLIED MATHEMATICS LABORATORY
C             NATIONAL INSTITUTE OF STANDARDS AND TECHNOLOGY
C             GAITHERSBURG, MD 20899
C           BOISVERT, R. F., (NIST)
C             COMPUTING AND APPLIED MATHEMATICS LABORATORY
C             NATIONAL INSTITUTE OF STANDARDS AND TECHNOLOGY
C             GAITHERSBURG, MD 20899
C***DESCRIPTION
C
C   VI1 computes the modified (hyperbolic) Bessel function of the
C   first kind of order one (I1) for real arguments using uniform
C   approximation by Chebyshev polynomials.
C
C
C   P A R A M E T E R S
C
C   M      (Input) Integer
C          The number of arguments at which the function is to be
C          evaluated.
C
C   X      (Input) Real array of length M
C          The arguments at which the function is to be evaluated are
C          stored in X(1) to X(M) in any order.
C
C   F      (Output) Real array of length M
C          F(i) contains the value of the function at X(i), i=1,..,M.
C
C   WORK   (Work) Real vector of length 7*M
C
C   IWORK  (Work) Integer vector of length M
C
C   INFO   (Output) Integer
C          Indicates status of computed result. The following table
C          lists possible values and their meanings.  If OK=Yes then
C          all F(i) have been set, otherwise none have been set.
C
C          INFO  OK            Description
C          ------------------------------------------------------------
C           -1   Yes  Warning: Some abs(X(i)) so small I1 underflows.
C                     The corresponding F(i) are set to zero.
C            0   Yes  Successfull execution.
C            1   No   Error: M .LE. 0
C            2   No   Error: Some abs(X(i)) so big I1 overflows.
C                     The index of the first offending argument is
C                     returned in IWORK(1).
C
C
C *********************************************************************
C   This routine is a modification of the function BESI1 developed by
C   W. Fullerton of LANL.
C *********************************************************************
C
C***REFERENCES  (NONE)
C***ROUTINES CALLED  WI1
C***REVISION HISTORY  (YYMMDD)
C   910226  DATE WRITTEN
C***END PROLOGUE  VI1
C
C     ----------
C     PARAMETERS
C     ----------
C
      INTEGER INFO, IWORK, M
      REAL    F, X, WORK
C
      DIMENSION X(M), F(M), WORK(7*M), IWORK(M)
C
C     ---------------
C     LOCAL VARIABLES
C     ---------------
C
      INTEGER IWB0, IWB1, IWB2, IWY, IWTC, IWYC, IWZC, JIN
C
C***FIRST EXECUTABLE STATEMENT  VI1
C
C     ... PARTITION WORK ARRAYS
C
      IWY   = 1
      IWTC  = IWY  + M
      IWYC  = IWTC + M
      IWZC  = IWYC + M
      IWB0  = IWZC + M
      IWB1  = IWB0 + M
      IWB2  = IWB1 + M
C     Total = IWB2 + M
C
      JIN   = 1
C     Total = JIN  + M
C
C     ... WI1 DOES ALL THE WORK
C
      CALL WI1(M,X,F,WORK(IWY),WORK(IWTC),WORK(IWYC),WORK(IWZC),
     +         IWORK(JIN),WORK(IWB0),WORK(IWB1),WORK(IWB2),INFO)
C
      RETURN
      END

      SUBROUTINE WI1(M, X, F, Y, TCMP, YCMP, ZCMP, INDX, B0, B1, B2, 
     +   INFO)
C***BEGIN PROLOGUE  WI1
C***SUBSIDIARY
C***PURPOSE  Computes the hyperbolic Bessel function of the first kind
C            of order one (I1) for a vector of arguments
C***LIBRARY   VFNLIB
C***AUTHOR  SAUNDERS, B. V., (NIST)
C             COMPUTING AND APPLIED MATHEMATICS LABORATORY
C             NATIONAL INSTITUTE OF STANDARDS AND TECHNOLOGY
C             GAITHERSBURG, MD 20899
C           BOISVERT, R. F., (NIST)
C             COMPUTING AND APPLIED MATHEMATICS LABORATORY
C             NATIONAL INSTITUTE OF STANDARDS AND TECHNOLOGY
C             GAITHERSBURG, MD 20899
C***ROUTINES CALLED  IWCS, R1MACH, WNGT, WGTHR, WGTLE, WGT, WLE,
C                    WSCTR, WCS
C***REVISION HISTORY  (YYMMDD)
C   910226  DATE WRITTEN
C***END PROLOGUE  WI1
C
C  ----------
C  PARAMETERS
C  ----------
C
      INTEGER LAI1, LAI12, LBI1
      PARAMETER ( LAI1=21, LAI12=22, LBI1=11 )
C
      INTEGER INFO, INDX, M
      REAL    B0, B1, B2, F, X, Y, TCMP, YCMP, ZCMP
C
      DIMENSION B0(M), B1(M), B2(M), F(M), INDX(M), X(M), Y(M),
     +          TCMP(M), YCMP(M), ZCMP(M)
C
C  ---------------
C  LOCAL VARIABLES
C  ---------------
C
      INTEGER I, IWCS, J, KEY, N, NTI1, NTAI1, NTAI12, NTOT
      REAL    AI1CS, AI12CS, BI1CS, EPMACH, EPS, R1MACH,
     +        XMAX, XMIN, XSML
C
      DIMENSION AI1CS(LAI1), AI12CS(LAI12), BI1CS(LBI1)
C
      SAVE AI1CS, AI12CS, BI1CS, N, NTAI1, NTAI12, NTI1,
     +     XMAX, XMIN, XSML
C
C----------------------------------------------------------------------
C
C Series for BI1        ON THE INTERVAL  0.          to  9.00000D+00
C                                        WITH WEIGHTED ERROR   2.40E-17
C                                         LOG WEIGHTED ERROR  16.62
C                               SIGNIFICANT FIGURES REQUIRED  16.23
C                                    DECIMAL PLACES REQUIRED  17.14
C
      DATA BI1 CS( 1) /   -.0019717132 61099859E0 /
      DATA BI1 CS( 2) /    .4073488766 7546481E0 /
      DATA BI1 CS( 3) /    .0348389942 99959456E0 /
      DATA BI1 CS( 4) /    .0015453945 56300123E0 /
      DATA BI1 CS( 5) /    .0000418885 21098377E0 /
      DATA BI1 CS( 6) /    .0000007649 02676483E0 /
      DATA BI1 CS( 7) /    .0000000100 42493924E0 /
      DATA BI1 CS( 8) /    .0000000000 99322077E0 /
      DATA BI1 CS( 9) /    .0000000000 00766380E0 /
      DATA BI1 CS(10) /    .0000000000 00004741E0 /
      DATA BI1 CS(11) /    .0000000000 00000024E0 /
C
C----------------------------------------------------------------------
C
C Series for AI1        on the interval  1.25000D-01 to  3.33333D-01
C                                        with weighted error   6.98E-17
C                                         log weighted error  16.16
C                               significant figures required  14.53
C                                    decimal places required  16.82
C
      DATA AI1 CS( 1) /   -.0284674418 1881479E0 /
      DATA AI1 CS( 2) /   -.0192295323 1443221E0 /
      DATA AI1 CS( 3) /   -.0006115185 8579437E0 /
      DATA AI1 CS( 4) /   -.0000206997 1253350E0 /
      DATA AI1 CS( 5) /    .0000085856 1914581E0 /
      DATA AI1 CS( 6) /    .0000010494 9824671E0 /
      DATA AI1 CS( 7) /   -.0000002918 3389184E0 /
      DATA AI1 CS( 8) /   -.0000000155 9378146E0 /
      DATA AI1 CS( 9) /    .0000000131 8012367E0 /
      DATA AI1 CS(10) /   -.0000000014 4842341E0 /
      DATA AI1 CS(11) /   -.0000000002 9085122E0 /
      DATA AI1 CS(12) /    .0000000001 2663889E0 /
      DATA AI1 CS(13) /   -.0000000000 1664947E0 /
      DATA AI1 CS(14) /   -.0000000000 0166665E0 /
      DATA AI1 CS(15) /    .0000000000 0124260E0 /
      DATA AI1 CS(16) /   -.0000000000 0027315E0 /
      DATA AI1 CS(17) /    .0000000000 0002023E0 /
      DATA AI1 CS(18) /    .0000000000 0000730E0 /
      DATA AI1 CS(19) /   -.0000000000 0000333E0 /
      DATA AI1 CS(20) /    .0000000000 0000071E0 /
      DATA AI1 CS(21) /   -.0000000000 0000006E0 /
C
C----------------------------------------------------------------------
C
C Series for AI12       on the interval  0.          to  1.25000D-01
C                                        with weighted error   3.55E-17
C                                         log weighted error  16.45
C                               significant figures required  14.69
C                                    decimal places required  17.12
C
      DATA AI12CS( 1) /    .0285762350 1828014E0 /
      DATA AI12CS( 2) /   -.0097610974 9136147E0 /
      DATA AI12CS( 3) /   -.0001105889 3876263E0 /
      DATA AI12CS( 4) /   -.0000038825 6480887E0 /
      DATA AI12CS( 5) /   -.0000002512 2362377E0 /
      DATA AI12CS( 6) /   -.0000000263 1468847E0 /
      DATA AI12CS( 7) /   -.0000000038 3538039E0 /
      DATA AI12CS( 8) /   -.0000000005 5897433E0 /
      DATA AI12CS( 9) /   -.0000000000 1897495E0 /
      DATA AI12CS(10) /    .0000000000 3252602E0 /
      DATA AI12CS(11) /    .0000000000 1412580E0 /
      DATA AI12CS(12) /    .0000000000 0203564E0 /
      DATA AI12CS(13) /   -.0000000000 0071985E0 /
      DATA AI12CS(14) /   -.0000000000 0040836E0 /
      DATA AI12CS(15) /   -.0000000000 0002101E0 /
      DATA AI12CS(16) /    .0000000000 0004273E0 /
      DATA AI12CS(17) /    .0000000000 0001041E0 /
      DATA AI12CS(18) /   -.0000000000 0000382E0 /
      DATA AI12CS(19) /   -.0000000000 0000186E0 /
      DATA AI12CS(20) /    .0000000000 0000033E0 /
      DATA AI12CS(21) /    .0000000000 0000028E0 /
      DATA AI12CS(22) /   -.0000000000 0000003E0 /
C
C----------------------------------------------------------------------
C
      DATA NTI1 / 0 /
C
C***FIRST EXECUTABLE STATEMENT  WI1
C
      IF (M .LE. 0)  GO TO 910
C
      IF (NTI1.EQ.0) THEN
         EPMACH = R1MACH(3)
         EPS = 0.10E0*EPMACH
         NTI1 = IWCS (BI1CS, LBI1, EPS)
         NTAI1 = IWCS (AI1CS, LAI1, EPS)
         NTAI12 = IWCS (AI12CS, LAI12, EPS)
         XMIN = 2.0E0*R1MACH(1)
         XSML = SQRT(8.0E0*EPMACH)
         XMAX = LOG(R1MACH(2))
      ENDIF
C
      NTOT = 0
C
      DO 10 I=1,M
         Y(I) = ABS(X(I))
  10  CONTINUE
C
      CALL WNGT(M,Y,XMAX,KEY)
      IF (KEY .NE. 0)  GO TO 920
C
C  ---------------------------
C  CASE   Y=0  OR  Y TOO SMALL
C  ---------------------------
C
C     NOTE -- I0 UNDERFLOWS FOR X .LE. XMIN
C
      DO 15 I=1,M
         F(I) = 0.0E0
  15  CONTINUE
C
C  ----------------------------
C  CASE   XMIN .LT. Y .LE. XSML
C  ----------------------------
C
      CALL WGTLE(M,Y,XMIN,XSML,N,INDX)
      IF (N .GT. 0) THEN
         NTOT = NTOT + N
         CALL WGTHR(N,Y,INDX,YCMP)
         DO 20 J=1,N
            ZCMP(J) = 0.50E0*YCMP(J)
  20     CONTINUE
         CALL WSCTR(N,ZCMP,INDX,F)
      ENDIF
C
C  ---------------------------
C  CASE   XSML .LT. Y .LE. 3.0
C  ---------------------------
C
      CALL WGTLE(M,Y,XSML,3.0E0,N,INDX)
      IF (N .GT. 0) THEN
         NTOT = NTOT + N
         CALL WGTHR(N,Y,INDX,YCMP)
         DO 30 J=1,N
            TCMP(J) = YCMP(J)**2/4.50E0 - 1.0E0
  30     CONTINUE
         CALL WCS(N,TCMP,BI1CS,NTI1,ZCMP,B0,B1,B2)
         DO 40 J=1,N
            ZCMP(J) = YCMP(J)*(0.8750E0 + ZCMP(J))
  40     CONTINUE
         CALL WSCTR(N,ZCMP,INDX,F)
      ENDIF
C
C  -------------------------
C  CASE  3.0 .LT. Y .LE. 8.0
C  -------------------------
C
      CALL WGTLE(M,Y,3.0E0,8.0E0,N,INDX)
      IF (N .GT. 0) THEN
         NTOT = NTOT + N
         CALL WGTHR(N,Y,INDX,YCMP)
         DO 70 J=1,N
            TCMP(J) = (48.0E0/YCMP(J) - 11.0E0)/5.0E0
  70     CONTINUE
         CALL WCS(N,TCMP,AI1CS,NTAI1,ZCMP,B0,B1,B2)
         DO 80 J=1,N
            ZCMP(J) = EXP(YCMP(J))*(0.3750E0 + ZCMP(J))/SQRT(YCMP(J))
  80     CONTINUE
         CALL WSCTR(N,ZCMP,INDX,F)
      ENDIF
C
C  ----------------
C  CASE  8.0 .LT. Y
C  ----------------
C
      CALL WGT(M,Y,8.0E0,N,INDX)
      IF (N .GT. 0) THEN
         NTOT = NTOT + N
         CALL WGTHR(N,Y,INDX,YCMP)
         DO 100 J=1,N
            TCMP(J) = 16.0E0/YCMP(J) - 1.0E0
  100     CONTINUE
         CALL WCS(N,TCMP,AI12CS,NTAI12,ZCMP,B0,B1,B2)
         DO 110 J=1,N
            ZCMP(J) = EXP(YCMP(J))*(0.3750E0 + ZCMP(J))/SQRT(YCMP(J))
  110    CONTINUE
         CALL WSCTR(N,ZCMP,INDX,F)
      ENDIF
C
C  -------------------------------------------
C  REVERSE SIGN FOR NEGATIVE X (RESULT IS ODD)
C  -------------------------------------------
C
      DO 200 I=1,M
         IF (X(I) .LT. 0.0E0)  F(I) = -F(I)
  200 CONTINUE
C
C  -----
C  EXITS
C  -----
C
C     ... NORMAL
C
      INFO = 0
      IF (NTOT .NE. M)  INFO = -1
      GO TO 999
C
C     ... M .LE. 0
C
  910 CONTINUE
      INFO = 1
      GO TO 999
C
C     ... ABS(X) SO LARGE I1 OVERFLOWS
C
  920 CONTINUE
      INFO = 2
      INDX(1) = KEY
      GO TO 999
C
  999 CONTINUE
      RETURN
      END
      SUBROUTINE VK0 (M, X, F, WORK, IWORK, INFO)
C***BEGIN PROLOGUE  VK0
C***PURPOSE  Computes the hyperbolic Bessel function of the third kind
C            of order zero (K0) for a vector of real arguments
C***LIBRARY   VFNLIB
C***CATEGORY  C10B1
C***TYPE      SINGLE PRECISION (VK0-S, DVK0-D)
C***KEYWORDS  BESSEL FUNCTION,THIRD KIND,HYPERBOLIC BESSEL FUNCTION,
C             MODIFIED BESSEL FUNCTION, ORDER ZERO, VECTORIZED
C***AUTHOR  SAUNDERS, B. V., (NIST)
C             COMPUTING AND APPLIED MATHEMATICS LABORATORY
C             NATIONAL INSTITUTE OF STANDARDS AND TECHNOLOGY
C             GAITHERSBURG, MD 20899
C           BOISVERT, R. F., (NIST)
C             COMPUTING AND APPLIED MATHEMATICS LABORATORY
C             NATIONAL INSTITUTE OF STANDARDS AND TECHNOLOGY
C             GAITHERSBURG, MD 20899
C***DESCRIPTION
C
C   VK0 computes the modified (hyperbolic) Bessel function of the
C   third kind of order zero (K0) for real arguments using uniform
C   approximation by Chebyshev polynomials.
C
C
C   P A R A M E T E R S
C
C   M      (Input) Integer
C          The number of arguments at which the function is to be
C          evaluated.
C
C   X      (Input) Real array of length M
C          The arguments at which the function is to be evaluated are
C          stored in X(1) to X(M) in any order.
C
C   F      (Output) Real array of length M
C          F(i) contains the value of the function at X(i), i=1,..,M.
C
C   WORK   (Work) Real vector of length 7*M
C
C   IWORK  (Work) Integer vector of length M
C
C   INFO   (Output) Integer
C          Indicates status of computed result. The following table
C          lists possible values and their meanings.  If OK=Yes then
C          all F(i) have been set, otherwise none have been set.
C
C          INFO  OK            Description
C          ------------------------------------------------------------
C           -1   Yes  Warning: Some X(i) so big K0 underflows.
C                     The corresponding F(i) are set to zero.
C            0   Yes  Successfull execution.
C            1   No   Error: M .LE. 0
C            2   No   Error: Some X(i) is zero or negative.
C                     The index of the first offending argument is
C                     returned in IWORK(1).
C
C
C *********************************************************************
C   This routine is a modification of the function BESK0 developed by
C   W. Fullerton of LANL.
C *********************************************************************
C
C***REFERENCES  (NONE)
C***ROUTINES CALLED  WK0
C***REVISION HISTORY  (YYMMDD)
C   910226  DATE WRITTEN
C***END PROLOGUE  VK0
C
C     ----------
C     PARAMETERS
C     ----------
C
      INTEGER INFO, IWORK, M
      REAL    F, X, WORK
C
      DIMENSION X(M), F(M), WORK(7*M), IWORK(M)
C
C     ---------------
C     LOCAL VARIABLES
C     ---------------
C
      INTEGER IWB0, IWB1, IWB2, IWTC, IWXC, IWYC, IWZC, JIN
C
C***FIRST EXECUTABLE STATEMENT  VK0
C
C     ... PARTITION WORK ARRAYS
C
      IWTC   = 1
      IWXC  = IWTC + M
      IWYC  = IWXC + M
      IWZC  = IWYC + M
      IWB0  = IWZC + M
      IWB1  = IWB0 + M
      IWB2  = IWB1 + M
C     Total = IB2  + M
C
      JIN   = 1
C     Total = JIN  + M
C
C     ... WK0 DOES ALL THE WORK
C
      CALL WK0(M,X,F,WORK(IWTC),WORK(IWXC),WORK(IWYC),WORK(IWZC),
     +         IWORK(JIN),WORK(IWB0),WORK(IWB1),WORK(IWB2),INFO)
C
      RETURN
      END


      SUBROUTINE WK0 (M, X, F, TCMP, XCMP, YCMP, ZCMP, INDX, B0, B1,  
     +   B2, INFO)
C***BEGIN PROLOGUE  WK0
C***SUBSIDIARY
C***PURPOSE  Computes the hyperbolic Bessel function of the third kind
C            of order zero (K0) for a vector of arguments
C***LIBRARY   VFNLIB
C***AUTHOR  SAUNDERS, B. V., (NIST)
C             COMPUTING AND APPLIED MATHEMATICS LABORATORY
C             NATIONAL INSTITUTE OF STANDARDS AND TECHNOLOGY
C             GAITHERSBURG, MD 20899
C           BOISVERT, R. F., (NIST)
C             COMPUTING AND APPLIED MATHEMATICS LABORATORY
C             NATIONAL INSTITUTE OF STANDARDS AND TECHNOLOGY
C             GAITHERSBURG, MD 20899
C***ROUTINES CALLED  IWCS, R1MACH, WNLE, WGTHR, WGTLE, WGT, WLE,
C                    WSCTR, WCS
C***REVISION HISTORY  (YYMMDD)
C   910226  DATE WRITTEN
C***END PROLOGUE  WK0
C
C  ----------
C  PARAMETERS
C  ----------
C
      INTEGER INFO, INDX, M
      REAL    B0, B1, B2, F, X, TCMP, XCMP, YCMP, ZCMP
C
      DIMENSION B0(M), B1(M), B2(M), F(M), INDX(M), X(M), TCMP(M),
     +          XCMP(M), YCMP(M), ZCMP(M)
C
C  ---------------
C  LOCAL VARIABLES
C  ---------------
C
      INTEGER LAK0, LAK02, LBI0, LBK0
      PARAMETER ( LAK0=17, LAK02=14, LBI0=12, LBK0=11 )
C
      INTEGER I, IWCS, J, KEY, N, NTI0, NTK0, NTAK0, NTAK02, NTOT
      REAL    AK0CS, AK02CS, BI0CS, BK0CS, EPMACH, EPS, R1MACH,
     +        XMAX, XSML
C
      DIMENSION  AK0CS(LAK0), AK02CS(LAK02), BI0CS(LBI0), BK0CS(LBK0)
C
      SAVE BK0CS, AK0CS, AK02CS, BI0CS, N, NTAK0, NTAK02, NTI0, NTK0,
     +     XMAX, XSML
C
C----------------------------------------------------------------------
C
C Series for BK0        on the interval  0.          to  4.00000D+00
C                                        with weighted error   3.57E-19
C                                         log weighted error  18.45
C                               significant figures required  17.99
C                                    decimal places required  18.97
C
      DATA BK0 CS( 1) /   -.0353273932 3390276872E0 /
      DATA BK0 CS( 2) /    .3442898999 246284869E0 /
      DATA BK0 CS( 3) /    .0359799365 1536150163E0 /
      DATA BK0 CS( 4) /    .0012646154 1144692592E0 /
      DATA BK0 CS( 5) /    .0000228621 2103119451E0 /
      DATA BK0 CS( 6) /    .0000002534 7910790261E0 /
      DATA BK0 CS( 7) /    .0000000019 0451637722E0 /
      DATA BK0 CS( 8) /    .0000000000 1034969525E0 /
      DATA BK0 CS( 9) /    .0000000000 0004259816E0 /
      DATA BK0 CS(10) /    .0000000000 0000013744E0 /
      DATA BK0 CS(11) /    .0000000000 0000000035E0 /
C
C----------------------------------------------------------------------
C
C Series for AK0        on the interval  1.25000D-01 to  5.00000D-01
C                                        with weighted error   5.34E-17
C                                         log weighted error  16.27
C                               significant figures required  14.92
C                                    decimal places required  16.89
C
      DATA AK0 CS( 1) /   -.0764394790 3327941E0 /
      DATA AK0 CS( 2) /   -.0223565260 5699819E0 /
      DATA AK0 CS( 3) /    .0007734181 1546938E0 /
      DATA AK0 CS( 4) /   -.0000428100 6688886E0 /
      DATA AK0 CS( 5) /    .0000030817 0017386E0 /
      DATA AK0 CS( 6) /   -.0000002639 3672220E0 /
      DATA AK0 CS( 7) /    .0000000256 3713036E0 /
      DATA AK0 CS( 8) /   -.0000000027 4270554E0 /
      DATA AK0 CS( 9) /    .0000000003 1694296E0 /
      DATA AK0 CS(10) /   -.0000000000 3902353E0 /
      DATA AK0 CS(11) /    .0000000000 0506804E0 /
      DATA AK0 CS(12) /   -.0000000000 0068895E0 /
      DATA AK0 CS(13) /    .0000000000 0009744E0 /
      DATA AK0 CS(14) /   -.0000000000 0001427E0 /
      DATA AK0 CS(15) /    .0000000000 0000215E0 /
      DATA AK0 CS(16) /   -.0000000000 0000033E0 /
      DATA AK0 CS(17) /    .0000000000 0000005E0 /
C
C----------------------------------------------------------------------
C
C Series for AK02       on the interval  0.          to  1.25000D-01
C                                        with weighted error   2.34E-17
C                                         log weighted error  16.63
C                               significant figures required  14.67
C                                    decimal places required  17.20
C
      DATA AK02CS( 1) /   -.0120186982 6307592E0 /
      DATA AK02CS( 2) /   -.0091748526 9102569E0 /
      DATA AK02CS( 3) /    .0001444550 9317750E0 /
      DATA AK02CS( 4) /   -.0000040136 1417543E0 /
      DATA AK02CS( 5) /    .0000001567 8318108E0 /
      DATA AK02CS( 6) /   -.0000000077 7011043E0 /
      DATA AK02CS( 7) /    .0000000004 6111825E0 /
      DATA AK02CS( 8) /   -.0000000000 3158592E0 /
      DATA AK02CS( 9) /    .0000000000 0243501E0 /
      DATA AK02CS(10) /   -.0000000000 0020743E0 /
      DATA AK02CS(11) /    .0000000000 0001925E0 /
      DATA AK02CS(12) /   -.0000000000 0000192E0 /
      DATA AK02CS(13) /    .0000000000 0000020E0 /
      DATA AK02CS(14) /   -.0000000000 0000002E0 /
C
C----------------------------------------------------------------------
C
C Series for BI0        on the interval  0.          to  9.00000D+00
C                                        with weighted error   2.46E-18
C                                         log weighted error  17.61
C                               significant figures required  17.90
C                                    decimal places required  18.15
C
      DATA BI0 CS( 1) /   -.0766054725 2839144951E0 /
      DATA BI0 CS( 2) /   1.9273379539 93808270E0 /
      DATA BI0 CS( 3) /    .2282644586 920301339E0 /
      DATA BI0 CS( 4) /    .0130489146 6707290428E0 /
      DATA BI0 CS( 5) /    .0004344270 9008164874E0 /
      DATA BI0 CS( 6) /    .0000094226 5768600193E0 /
      DATA BI0 CS( 7) /    .0000001434 0062895106E0 /
      DATA BI0 CS( 8) /    .0000000016 1384906966E0 /
      DATA BI0 CS( 9) /    .0000000000 1396650044E0 /
      DATA BI0 CS(10) /    .0000000000 0009579451E0 /
      DATA BI0 CS(11) /    .0000000000 0000053339E0 /
      DATA BI0 CS(12) /    .0000000000 0000000245E0 /
C
C----------------------------------------------------------------------
C
      DATA NTK0 / 0 /
C
C***FIRST EXECUTABLE STATEMENT  WK0
C
      IF (M .LE. 0)  GO TO 910
C
      IF (NTK0 .EQ. 0) THEN
         EPMACH = R1MACH(3)
         EPS = 0.10E0*EPMACH
         NTK0 = IWCS(BK0CS, LBK0, EPS)
         NTAK0 = IWCS(AK0CS, LAK0, EPS)
         NTAK02 = IWCS(AK02CS, LAK02, EPS)
         NTI0 = IWCS(BI0CS, LBI0, EPS)
         XSML = SQRT (4.0E0*EPMACH)
         XMAX = -LOG(R1MACH(1))
         XMAX = XMAX - 0.50E0*XMAX*LOG(XMAX)/(XMAX+0.50E0) - 0.010E0
      ENDIF
C
      NTOT = 0
C
      CALL WNLE(M,X,0.0E0,KEY)
      IF (KEY .NE. 0)  GO TO 920
C
C  ------------------
C  CASE   X .GT. XMAX
C  ------------------
C
C     NOTE -- K0 UNDERFLOWS FOR X .GT. XMAX
C
      DO 5 I=1,M
         F(I) = 0.0E0
   5  CONTINUE
C
C  ----------------
C  CASE  X .LE. 2.0
C  ----------------
C
      CALL WLE(M,X,2.0E0,N,INDX)
      IF (N .GT. 0) THEN
         NTOT = NTOT + N
         CALL WGTHR(N,X,INDX,XCMP)
C
C        ... COMPUTE I0(X) ... RESULT IN ZCMP
C
         DO 10 J=1,N
            TCMP(J) = XCMP(J)**2/4.50E0 - 1.0E0
   10    CONTINUE
         CALL WCS(N,TCMP,BI0CS,NTI0,ZCMP,B0,B1,B2)
         DO 15 J=1,N
            ZCMP(J) = 2.750E0 + ZCMP(J)
   15    CONTINUE
C
         DO 20 J=1,N
            TCMP(J) = 0.50E0*XCMP(J)**2 - 1.0E0
   20    CONTINUE
         CALL WCS(N,TCMP,BK0CS,NTK0,YCMP,B0,B1,B2)
         DO 30 J=1,N
            YCMP(J) = -LOG(0.50E0*XCMP(J))*ZCMP(J) - 0.250E0 + YCMP(J)
   30    CONTINUE
         CALL WSCTR(N,YCMP,INDX,F)
      ENDIF
C
C  -------------------------
C  CASE  2.0 .LT. X .LE. 8.0
C  -------------------------
C
      CALL WGTLE(M,X,2.0E0,8.0E0,N,INDX)
      IF (N .GT. 0) THEN
         NTOT = NTOT + N
         CALL WGTHR(N,X,INDX,XCMP)
         DO 40 J=1,N
            TCMP(J) = (16.0E0/XCMP(J) - 5.0E0)/3.0E0
   40    CONTINUE
         CALL WCS(N,TCMP,AK0CS,NTAK0,ZCMP,B0,B1,B2)
         DO 50 J=1,N
            YCMP(J) = EXP(-XCMP(J))*(1.250E0 + ZCMP(J))/SQRT(XCMP(J))
   50    CONTINUE
         CALL WSCTR(N,YCMP,INDX,F)
      ENDIF
C
C  --------------------------
C  CASE  8.0 .LT. X .LE. XMAX
C  --------------------------
C
      CALL WGTLE(M,X,8.0E0,XMAX,N,INDX)
      IF (N .GT. 0) THEN
         NTOT = NTOT + N
         CALL WGTHR(N,X,INDX,XCMP)
         DO 60 J=1,N
            TCMP(J) = 16.0E0/XCMP(J) - 1.0E0
   60    CONTINUE
         CALL WCS(N,TCMP,AK02CS,NTAK02,ZCMP,B0,B1,B2)
         DO 70 J=1,N
            YCMP(J) = EXP(-XCMP(J))*(1.250E0 + ZCMP(J))/SQRT(XCMP(J))
   70    CONTINUE
         CALL WSCTR(N,YCMP,INDX,F)
      ENDIF
C
C  -----
C  EXITS
C  -----
C
C     ... NORMAL
C
      INFO = 0
      GO TO 999
C
C     ... M .LE. 0
C
  910 CONTINUE
      INFO = 1
      IF (NTOT .NE. M)  INFO = -1
      GO TO 999
C
C     ... X IS ZERO OR NEGATIVE
C
  920 CONTINUE
      INFO = 2
      INDX(1) = KEY
      GO TO 999
C
  999 CONTINUE
      RETURN
      END


      SUBROUTINE VK1 (M, X, F, WORK, IWORK, INFO)
C***BEGIN PROLOGUE  VK1
C***PURPOSE  Computes the hyperbolic Bessel function of the third kind
C            of order one (K1) for a vector of real arguments
C***LIBRARY   VFNLIB
C***CATEGORY  C10B1
C***TYPE      SINGLE PRECISION (VK1-S, DVK1-D)
C***KEYWORDS  BESSEL FUNCTION,THIRD KIND,HYPERBOLIC BESSEL FUNCTION,
C             MODIFIED BESSEL FUNCTION, ORDER ONE, VECTORIZED
C***AUTHOR  SAUNDERS, B. V., (NIST)
C             COMPUTING AND APPLIED MATHEMATICS LABORATORY
C             NATIONAL INSTITUTE OF STANDARDS AND TECHNOLOGY
C             GAITHERSBURG, MD 20899
C           BOISVERT, R. F., (NIST)
C             COMPUTING AND APPLIED MATHEMATICS LABORATORY
C             NATIONAL INSTITUTE OF STANDARDS AND TECHNOLOGY
C             GAITHERSBURG, MD 20899
C***DESCRIPTION
C
C   VK1 computes the modified (hyperbolic) Bessel function of the
C   third kind of order one (K1) for real arguments using uniform
C   approximation by Chebyshev polynomials.
C
C
C   P A R A M E T E R S
C
C   M      (Input) Integer
C          The number of arguments at which the function is to be
C          evaluated.
C
C   X      (Input) Real array of length M
C          The arguments at which the function is to be evaluated are
C          stored in X(1) to X(M) in any order.
C
C   F      (Output) Real array of length M
C          F(i) contains the value of the function at X(i), i=1,..,M.
C
C   WORK   (Work) Real vector of length 7*M
C
C   IWORK  (Work) Integer vector of length M
C
C   INFO   (Output) Integer
C          Indicates status of computed result. The following table
C          lists possible values and their meanings.  If OK=Yes then
C          all F(i) have been set, otherwise none have been set.
C
C          INFO  OK            Description
C          ------------------------------------------------------------
C           -1   Yes  Warning: Some X(i) so big K1 underflows.
C                     The corresponding F(i) are set to zero.
C            0   Yes  Successfull execution.
C            1   No   Error: M .LE. 0
C            2   No   Error: Some X(i) is zero or negative.
C                     The index of the first offending argument is
C                     returned in IWORK(1).
C            3   No   Error: Some X(i) so small K1 overflows.
C                     The index of the first offending argument is
C                     returned in IWORK(1).
C
C
C *********************************************************************
C   This routine is a modification of the function BESK1 developed by
C   W. Fullerton of LANL.
C *********************************************************************
C
C***REFERENCES  (NONE)
C***ROUTINES CALLED  WK1
C***REVISION HISTORY  (YYMMDD)
C   910226  DATE WRITTEN
C***END PROLOGUE  VK1
C
C     ----------
C     PARAMETERS
C     ----------
C
      INTEGER INFO, IWORK, M
      REAL    F, X, WORK
C
      DIMENSION X(M), F(M), WORK(7*M), IWORK(M)
C
C     ---------------
C     LOCAL VARIABLES
C     ---------------
C
      INTEGER IWB0, IWB1, IWB2, IWTC, IWXC, IWYC, IWZC, JIN
C
C***FIRST EXECUTABLE STATEMENT  VK1
C
C     ... PARTITION WORK ARRAYS
C
      IWTC   = 1
      IWXC  = IWTC + M
      IWYC  = IWXC + M
      IWZC  = IWYC + M
      IWB0  = IWZC + M
      IWB1  = IWB0 + M
      IWB2  = IWB1 + M
C     Total = IB2  + M
C
      JIN   = 1
C     Total = JIN  + M
C
C     ... WK1 DOES ALL THE WORK
C
      CALL WK1(M,X,F,WORK(IWTC),WORK(IWXC),WORK(IWYC),WORK(IWZC),
     +         IWORK(JIN),WORK(IWB0),WORK(IWB1),WORK(IWB2),INFO)
C
      RETURN
      END

      SUBROUTINE WK1 (M, X, F, TCMP, XCMP, YCMP, ZCMP, INDX, B0, B1,  
     +   B2, INFO)
C***BEGIN PROLOGUE  WK1
C***SUBSIDIARY
C***PURPOSE  Computes the hyperbolic Bessel function of the third kind
C            of order one (K1) for a vector of arguments
C***LIBRARY   VFNLIB
C***AUTHOR  SAUNDERS, B. V., (NIST)
C             COMPUTING AND APPLIED MATHEMATICS LABORATORY
C             NATIONAL INSTITUTE OF STANDARDS AND TECHNOLOGY
C             GAITHERSBURG, MD 20899
C           BOISVERT, R. F., (NIST)
C             COMPUTING AND APPLIED MATHEMATICS LABORATORY
C             NATIONAL INSTITUTE OF STANDARDS AND TECHNOLOGY
C             GAITHERSBURG, MD 20899
C***ROUTINES CALLED  IWCS, R1MACH, WNLE, WGTHR, WGTLE, WGT, WLE,
C                    WSCTR, WCS
C***REVISION HISTORY  (YYMMDD)
C   910226  DATE WRITTEN
C***END PROLOGUE  WK1
C
C  ----------
C  PARAMETERS
C  ----------
C
      INTEGER INFO, INDX, M
      REAL    B0, B1, B2, F, X, TCMP, XCMP, YCMP, ZCMP
C
      DIMENSION B0(M), B1(M), B2(M), F(M), INDX(M), X(M), TCMP(M),
     +          XCMP(M), YCMP(M), ZCMP(M)
C
C  ---------------
C  LOCAL VARIABLES
C  ---------------
C
      INTEGER LAK1, LAK12, LBI1, LBK1
      PARAMETER ( LAK1=17, LAK12=14, LBI1=11, LBK1=11 )
C
      INTEGER I, IWCS, J, KEY, N, NTI1, NTK1, NTAK1, NTAK12, NTOT
      REAL    AK1CS, AK12CS, BI1CS, BK1CS, EPMACH, EPS, R1MACH,
     +        XMAX, XMIN, XSML
C
      DIMENSION AK1CS(LAK1), AK12CS(LAK12), BI1CS(LBI1), BK1CS(LBK1)
C
      SAVE AK1CS, AK12CS, BI1CS, BK1CS, NTAK1, NTAK12, NTI1, NTK1,
     +     XMAX, XMIN, XSML
C
C----------------------------------------------------------------------
C
C Series for BK1        on the interval  0.          to  4.00000D+00
C                                        with weighted error   7.02E-18
C                                         log weighted error  17.15
C                               significant figures required  16.73
C                                    decimal places required  17.67
C
      DATA BK1 CS( 1) /    .0253002273 389477705E0 /
      DATA BK1 CS( 2) /   -.3531559607 76544876E0 /
      DATA BK1 CS( 3) /   -.1226111808 22657148E0 /
      DATA BK1 CS( 4) /   -.0069757238 596398643E0 /
      DATA BK1 CS( 5) /   -.0001730288 957513052E0 /
      DATA BK1 CS( 6) /   -.0000024334 061415659E0 /
      DATA BK1 CS( 7) /   -.0000000221 338763073E0 /
      DATA BK1 CS( 8) /   -.0000000001 411488392E0 /
      DATA BK1 CS( 9) /   -.0000000000 006666901E0 /
      DATA BK1 CS(10) /   -.0000000000 000024274E0 /
      DATA BK1 CS(11) /   -.0000000000 000000070E0 /
C
C----------------------------------------------------------------------
C
C Series for BI1        ON THE INTERVAL  0.          to  9.00000D+00
C                                        WITH WEIGHTED ERROR   2.40E-17
C                                         LOG WEIGHTED ERROR  16.62
C                               SIGNIFICANT FIGURES REQUIRED  16.23
C                                    DECIMAL PLACES REQUIRED  17.14
C
      DATA BI1 CS( 1) /   -.0019717132 61099859E0 /
      DATA BI1 CS( 2) /    .4073488766 7546481E0 /
      DATA BI1 CS( 3) /    .0348389942 99959456E0 /
      DATA BI1 CS( 4) /    .0015453945 56300123E0 /
      DATA BI1 CS( 5) /    .0000418885 21098377E0 /
      DATA BI1 CS( 6) /    .0000007649 02676483E0 /
      DATA BI1 CS( 7) /    .0000000100 42493924E0 /
      DATA BI1 CS( 8) /    .0000000000 99322077E0 /
      DATA BI1 CS( 9) /    .0000000000 00766380E0 /
      DATA BI1 CS(10) /    .0000000000 00004741E0 /
      DATA BI1 CS(11) /    .0000000000 00000024E0 /
C
C----------------------------------------------------------------------
C
C Series for AK1        on the interval  1.25000D-01 to  5.00000D-01
C                                        with weighted error   6.06E-17
C                                         log weighted error  16.22
C                               significant figures required  15.41
C                                    decimal places required  16.83
C
      DATA AK1 CS( 1) /    .2744313406 973883E0 /
      DATA AK1 CS( 2) /    .0757198995 3199368E0 /
      DATA AK1 CS( 3) /   -.0014410515 5647540E0 /
      DATA AK1 CS( 4) /    .0000665011 6955125E0 /
      DATA AK1 CS( 5) /   -.0000043699 8470952E0 /
      DATA AK1 CS( 6) /    .0000003540 2774997E0 /
      DATA AK1 CS( 7) /   -.0000000331 1163779E0 /
      DATA AK1 CS( 8) /    .0000000034 4597758E0 /
      DATA AK1 CS( 9) /   -.0000000003 8989323E0 /
      DATA AK1 CS(10) /    .0000000000 4720819E0 /
      DATA AK1 CS(11) /   -.0000000000 0604783E0 /
      DATA AK1 CS(12) /    .0000000000 0081284E0 /
      DATA AK1 CS(13) /   -.0000000000 0011386E0 /
      DATA AK1 CS(14) /    .0000000000 0001654E0 /
      DATA AK1 CS(15) /   -.0000000000 0000248E0 /
      DATA AK1 CS(16) /    .0000000000 0000038E0 /
      DATA AK1 CS(17) /   -.0000000000 0000006E0 /
C
C----------------------------------------------------------------------
C
C Series for AK12       on the interval  0.          to  1.25000D-01
C                                        with weighted error   2.58E-17
C                                         log weighted error  16.59
C                               significant figures required  15.22
C                                    decimal places required  17.16
C
      DATA AK12CS( 1) /    .0637930834 3739001E0 /
      DATA AK12CS( 2) /    .0283288781 3049721E0 /
      DATA AK12CS( 3) /   -.0002475370 6739052E0 /
      DATA AK12CS( 4) /    .0000057719 7245160E0 /
      DATA AK12CS( 5) /   -.0000002068 9392195E0 /
      DATA AK12CS( 6) /    .0000000097 3998344E0 /
      DATA AK12CS( 7) /   -.0000000005 5853361E0 /
      DATA AK12CS( 8) /    .0000000000 3732996E0 /
      DATA AK12CS( 9) /   -.0000000000 0282505E0 /
      DATA AK12CS(10) /    .0000000000 0023720E0 /
      DATA AK12CS(11) /   -.0000000000 0002176E0 /
      DATA AK12CS(12) /    .0000000000 0000215E0 /
      DATA AK12CS(13) /   -.0000000000 0000022E0 /
      DATA AK12CS(14) /    .0000000000 0000002E0 /
C
C----------------------------------------------------------------------
C
      DATA NTK1 / 0 /
C
C***FIRST EXECUTABLE STATEMENT  WK1
C
      IF (M .LE. 0)  GO TO 910
C
      IF (NTK1 .EQ. 0) THEN
         EPMACH = R1MACH(3)
         EPS = 0.10E0*EPMACH
         NTK1 = IWCS(BK1CS, LBK1, EPS)
         NTI1 = IWCS(BI1CS, LBI1, EPS)
         NTAK1 = IWCS(AK1CS, LAK1, EPS)
         NTAK12 = IWCS(AK12CS, LAK12, EPS)
         XMIN = EXP(MAX( LOG(R1MACH(1)), -LOG(R1MACH(2))) + 0.010E0)
         XSML = SQRT(4.0E0*EPMACH)
         XMAX = -LOG(R1MACH(1))
         XMAX = XMAX - 0.50E0*XMAX*LOG(XMAX)/(XMAX + 0.50E0)
      ENDIF
C
      NTOT = 0
C
      CALL WNLE(M,X,0.0E0,KEY)
      IF (KEY .NE. 0)  GO TO 920
C
      CALL WNLE(M,X,XMIN,KEY)
      IF (KEY .NE. 0)  GO TO 930
C
C  ------------------
C  CASE   X .GT. XMAX
C  ------------------
C
C     NOTE -- K0 UNDERFLOWS FOR X .GT. XMAX
C
      DO 5 I=1,M
         F(I) = 0.0E0
   5  CONTINUE
C
C  ----------------
C  CASE  X .LE. 2.0
C  ----------------
C
      CALL WLE(M,X,2.0E0,N,INDX)
      IF (N .GT. 0) THEN
         NTOT = NTOT + N
         CALL WGTHR(N,X,INDX,XCMP)
C
C        ... COMPUTE I1(X) ... RESULT IN ZCMP
C
         DO 20 J=1,N
            TCMP(J) = XCMP(J)**2/4.50E0 - 1.0E0
   20    CONTINUE
         CALL WCS(N,TCMP,BI1CS,NTI1,ZCMP,B0,B1,B2)
         DO 30 J=1,N
            ZCMP(J) = XCMP(J)*(0.8750E0 + ZCMP(J))
   30    CONTINUE
C
         DO 40 J=1,N
            TCMP(J) = 0.50E0*XCMP(J)**2 - 1.0E0
   40    CONTINUE
         CALL WCS(N,TCMP,BK1CS,NTK1,YCMP,B0,B1,B2)
         DO 50 J=1,N
            ZCMP(J) = LOG(0.50E0*XCMP(J))*ZCMP(J) +
     +                (0.750E0 + YCMP(J))/XCMP(J)
   50    CONTINUE
         CALL WSCTR(N,ZCMP,INDX,F)
      ENDIF
C
C  -------------------------
C  CASE  2.0 .LT. X .LE. 8.0
C  -------------------------
C
      CALL WGTLE(M,X,2.0E0,8.0E0,N,INDX)
      IF (N .GT. 0) THEN
         NTOT = NTOT + N
         CALL WGTHR(N,X,INDX,XCMP)
         DO 60 J=1,N
            TCMP(J) = (16.0E0/XCMP(J) - 5.0E0)/3.0E0
   60    CONTINUE
         CALL WCS(N,TCMP,AK1CS,NTAK1,YCMP,B0,B1,B2)
         DO 70 J=1,N
            ZCMP(J) = EXP(-XCMP(J))*(1.250E0 + YCMP(J))/SQRT(XCMP(J))
   70    CONTINUE
         CALL WSCTR(N,ZCMP,INDX,F)
      ENDIF
C
C  --------------------------
C  CASE  8.0 .LT. X .LE. XMAX
C  --------------------------
C
      CALL WGTLE(M,X,8.0E0,XMAX,N,INDX)
      IF (N .GT. 0) THEN
         NTOT = NTOT + N
         CALL WGTHR(N,X,INDX,XCMP)
         DO 80 J=1,N
            TCMP(J) = 16.0E0/XCMP(J) - 1.0E0
   80    CONTINUE
         CALL WCS(N,TCMP,AK12CS,NTAK12,YCMP,B0,B1,B2)
         DO 90 J=1,N
            ZCMP(J) = EXP(-XCMP(J))*(1.250E0 + YCMP(J))/SQRT(XCMP(J))
   90    CONTINUE
         CALL WSCTR(N,ZCMP,INDX,F)
      ENDIF
C
C  -----
C  EXITS
C  -----
C
C     ... NORMAL
C
      INFO = 0
      IF (NTOT .NE. M)  INFO = -1
      GO TO 999
C
C     ... M .LE. 0
C
  910 CONTINUE
      INFO = 1
      GO TO 999
C
C     ... X IS ZERO OR NEGATIVE
C
  920 CONTINUE
      INFO = 2
      INDX(1) = KEY
      GO TO 999
C
C     ... X SO SMALL K1 OVERFLOWS
C
  930 CONTINUE
      INFO = 3
      INDX(1) = KEY
      GO TO 999
C
  999 CONTINUE
      RETURN
      END
      SUBROUTINE VJ0 (M, X, F, WORK, IWORK, INFO)
C***BEGIN PROLOGUE  VJ0
C***PURPOSE  Computes the Bessel function of the first kind of order
C            zero (J0) for a vector of real arguments
C***LIBRARY   VFNLIB
C***CATEGORY  C10B1
C***TYPE      SINGLE PRECISION (VJ0-S, DVJ0-D)
C***KEYWORDS  BESSEL FUNCTION,FIRST KIND, ORDER ZERO, SPECIAL FUNCTION,
C             VECTORIZED
C***AUTHOR  SAUNDERS, B. V., (NIST)
C             COMPUTING AND APPLIED MATHEMATICS LABORATORY
C             NATIONAL INSTITUTE OF STANDARDS AND TECHNOLOGY
C             GAITHERSBURG, MD 20899
C           BOISVERT, R. F., (NIST)
C             COMPUTING AND APPLIED MATHEMATICS LABORATORY
C             NATIONAL INSTITUTE OF STANDARDS AND TECHNOLOGY
C             GAITHERSBURG, MD 20899
C***DESCRIPTION
C
C   VJ0 computes the Bessel function of the first kind of order zero
C   (J0) for real arguments using uniform approximation by Chebyshev
C   polynomials.
C
C
C   P A R A M E T E R S
C
C   M      (Input) Integer (M .GT. 0)
C          The number of arguments at which the function is to be
C          evaluated.
C
C   X      (Input) Real array of length M
C          The arguments at which the function is to be evaluated are
C          stored in X(1) to X(M) in any order.
C
C   F      (Output) Real array of length M
C          F(i) contains the value of the function at X(i), i=1,..,M.
C
C   WORK   (Work) Real vector of length 7*M
C
C   IWORK  (Work) Integer vector of length M
C
C   INFO   (Output) Integer
C          Indicates status of computed result. The following table
C          lists possible values and their meanings.  If OK=Yes then
C          all F(i) have been set, otherwise none have been set.
C
C          INFO  OK            Description
C          ------------------------------------------------------------
C            0   Yes  Successfull execution.
C            1   No   Error: M .LE. 0
C            2   No   Error: Some abs(X(i)) so big that no precision
C                     possible in computing J0.  The index of the
C                     first offending argument is returned in IWORK(1).
C
C *********************************************************************
C   This routine is a modification of the function BESJ0 developed by
C   W. Fullerton of LANL.
C *********************************************************************
C
C***REFERENCES  (NONE)
C***ROUTINES CALLED  WJ0
C***REVISION HISTORY  (YYMMDD)
C   910226  DATE WRITTEN
C***END PROLOGUE  VJ0
C
C     ----------
C     PARAMETERS
C     ----------
C
      INTEGER INFO, IWORK, M
      REAL    F, X, WORK
C
      DIMENSION X(M), F(M), WORK(7*M), IWORK(M)
C
C     ---------------
C     LOCAL VARIABLES
C     ---------------
C
      INTEGER IWB0, IWB1, IWB2, IWY, IWTC, IWYC, IWZC, JIN
C
C
C***FIRST EXECUTABLE STATEMENT  VJ0
C
C     ... PARTITION WORK ARRAYS
C
      IWY   = 1
      IWTC  = IWY  + M
      IWYC  = IWTC + M
      IWZC  = IWYC + M
      IWB0  = IWZC + M
      IWB1  = IWB0 + M
      IWB2  = IWB1 + M
C     Total = IWB2 + M
C
      JIN   = 1
C     Total = JIN  + M
C
C     ... WJ0 DOES ALL THE WORK
C
      CALL WJ0(M,X,F,WORK(IWY),WORK(IWTC),WORK(IWYC),WORK(IWZC),
     +         IWORK(JIN),WORK(IWB0),WORK(IWB1),WORK(IWB2),INFO)
C
      RETURN
      END
      SUBROUTINE WJ0 (M, X, F, Y, TCMP, YCMP, ZCMP, INDX, B0, B1, B2, 
     +   INFO)
C***BEGIN PROLOGUE  WJ0
C***SUBSIDIARY
C***PURPOSE  Computes the Bessel function of the first kind
C            of order zero (J0) for a vector of arguments
C***LIBRARY   VFNLIB
C***AUTHOR  SAUNDERS, B. V., (NIST)
C             COMPUTING AND APPLIED MATHEMATICS LABORATORY
C             NATIONAL INSTITUTE OF STANDARDS AND TECHNOLOGY
C             GAITHERSBURG, MD 20899
C           BOISVERT, R. F., (NIST)
C             COMPUTING AND APPLIED MATHEMATICS LABORATORY
C             NATIONAL INSTITUTE OF STANDARDS AND TECHNOLOGY
C             GAITHERSBURG, MD 20899
C***ROUTINES CALLED  IWCS, R1MACH, WNGT, WGTHR, WGTLE, WGT, WLE,
C                    WSCTR, WCS
C***REVISION HISTORY  (YYMMDD)
C   910226  DATE WRITTEN
C***END PROLOGUE  WJ0
C
C  ----------
C  PARAMETERS
C  ----------
C
      INTEGER INFO, INDX, M
      REAL    B0, B1, B2, F, X, Y, TCMP, YCMP, ZCMP
C
      DIMENSION B0(M), B1(M), B2(M), F(M), INDX(M), X(M), Y(M),
     +          TCMP(M), YCMP(M), ZCMP(M)
C
C  ---------------
C  LOCAL VARIABLES
C  ---------------
C
      INTEGER LBJ0, LBM0, LBM02, LBTH0, LBT02
      PARAMETER (LBJ0=13, LBM0=21, LBM02=2, LBTH0=24, LBT02=2)
C
      INTEGER I, IWCS, J, JH, KEY, N, NA, NB, NTJ0, NTM0, NTM02,
     +        NTTH0, NTT02
      REAL    BJ0CS, BM0CS, BM02CS, BTH0CS, BT02CS, C1, C2,
     +        EPMACH, EPS, PI4, R1MACH, XSML, XMAX
C
      DIMENSION BJ0CS(LBJ0), BM0CS(LBM0), BM02CS(LBM02), BTH0CS(LBTH0),
     +          BT02CS(LBT02)
C
      SAVE BJ0CS, BM0CS, BM02CS, BTH0CS, BT02CS, NTJ0, NTM0, NTM02,
     +     NTTH0, NTT02, PI4, XSML, XMAX
C
C----------------------------------------------------------------------
C
C Series for BJ0        on the interval  0.          to  1.60000D+01
C                                        with weighted error   7.47E-18
C                                         log weighted error  17.13
C                               significant figures required  16.98
C                                    decimal places required  17.68
C
      DATA BJ0 CS( 1) /    .1002541619 68939137E0 /
      DATA BJ0 CS( 2) /   -.6652230077 64405132E0 /
      DATA BJ0 CS( 3) /    .2489837034 98281314E0 /
      DATA BJ0 CS( 4) /   -.0332527231 700357697E0 /
      DATA BJ0 CS( 5) /    .0023114179 304694015E0 /
      DATA BJ0 CS( 6) /   -.0000991127 741995080E0 /
      DATA BJ0 CS( 7) /    .0000028916 708643998E0 /
      DATA BJ0 CS( 8) /   -.0000000612 108586630E0 /
      DATA BJ0 CS( 9) /    .0000000009 838650793E0 /
      DATA BJ0 CS(10) /   -.0000000000 124235515E0 /
      DATA BJ0 CS(11) /    .0000000000 001265433E0 /
      DATA BJ0 CS(12) /   -.0000000000 000010619E0 /
      DATA BJ0 CS(13) /    .0000000000 000000074E0 /
C
C----------------------------------------------------------------------
C
C Series for BM0        on the interval  0.          to  6.25000D-02
C                                        with weighted error   4.98E-17
C                                         log weighted error  16.30
C                               significant figures required  14.97
C                                    decimal places required  16.96
C
      DATA BM0 CS( 1) /    .0928496163 7381644E0 /
      DATA BM0 CS( 2) /   -.0014298770 7403484E0 /
      DATA BM0 CS( 3) /    .0000283057 9271257E0 /
      DATA BM0 CS( 4) /   -.0000014330 0611424E0 /
      DATA BM0 CS( 5) /    .0000001202 8628046E0 /
      DATA BM0 CS( 6) /   -.0000000139 7113013E0 /
      DATA BM0 CS( 7) /    .0000000020 4076188E0 /
      DATA BM0 CS( 8) /   -.0000000003 5399669E0 /
      DATA BM0 CS( 9) /    .0000000000 7024759E0 /
      DATA BM0 CS(10) /   -.0000000000 1554107E0 /
      DATA BM0 CS(11) /    .0000000000 0376226E0 /
      DATA BM0 CS(12) /   -.0000000000 0098282E0 /
      DATA BM0 CS(13) /    .0000000000 0027408E0 /
      DATA BM0 CS(14) /   -.0000000000 0008091E0 /
      DATA BM0 CS(15) /    .0000000000 0002511E0 /
      DATA BM0 CS(16) /   -.0000000000 0000814E0 /
      DATA BM0 CS(17) /    .0000000000 0000275E0 /
      DATA BM0 CS(18) /   -.0000000000 0000096E0 /
      DATA BM0 CS(19) /    .0000000000 0000034E0 /
      DATA BM0 CS(20) /   -.0000000000 0000012E0 /
      DATA BM0 CS(21) /    .0000000000 0000004E0 /
C
C----------------------------------------------------------------------
C
C Series for BM02       Used in double precision version only
C
      DATA BM02CS( 1) /   1.0E0 /
      DATA BM02CS( 2) /   0.0E0 /
C
C----------------------------------------------------------------------
C
C Series for BTH0       on the interval  0.          to  6.25000D-02
C                                        with weighted error   3.67E-17
C                                         log weighted error  16.44
C                               significant figures required  15.53
C                                    decimal places required  17.13
C
      DATA BTH0CS( 1) /   -.2463916377 4300119E0 /
      DATA BTH0CS( 2) /    .0017370983 07508963E0 /
      DATA BTH0CS( 3) /   -.0000621836 33402968E0 /
      DATA BTH0CS( 4) /    .0000043680 50165742E0 /
      DATA BTH0CS( 5) /   -.0000004560 93019869E0 /
      DATA BTH0CS( 6) /    .0000000621 97400101E0 /
      DATA BTH0CS( 7) /   -.0000000103 00442889E0 /
      DATA BTH0CS( 8) /    .0000000019 79526776E0 /
      DATA BTH0CS( 9) /   -.0000000004 28198396E0 /
      DATA BTH0CS(10) /    .0000000001 02035840E0 /
      DATA BTH0CS(11) /   -.0000000000 26363898E0 /
      DATA BTH0CS(12) /    .0000000000 07297935E0 /
      DATA BTH0CS(13) /   -.0000000000 02144188E0 /
      DATA BTH0CS(14) /    .0000000000 00663693E0 /
      DATA BTH0CS(15) /   -.0000000000 00215126E0 /
      DATA BTH0CS(16) /    .0000000000 00072659E0 /
      DATA BTH0CS(17) /   -.0000000000 00025465E0 /
      DATA BTH0CS(18) /    .0000000000 00009229E0 /
      DATA BTH0CS(19) /   -.0000000000 00003448E0 /
      DATA BTH0CS(20) /    .0000000000 00001325E0 /
      DATA BTH0CS(21) /   -.0000000000 00000522E0 /
      DATA BTH0CS(22) /    .0000000000 00000210E0 /
      DATA BTH0CS(23) /   -.0000000000 00000087E0 /
      DATA BTH0CS(24) /    .0000000000 00000036E0 /
C
C----------------------------------------------------------------------
C
C Series for BT02       Used in double precision version only
C
      DATA BT02CS( 1) /   1.0E0 /
      DATA BT02CS( 2) /   0.0E0 /
C
C----------------------------------------------------------------------
C
      DATA NTJ0 / 0 /
C
C***FIRST EXECUTABLE STATEMENT  WJ0
C
      IF (M .LE. 0) GO TO 910
C
      IF (NTJ0.EQ.0) THEN
         PI4 = ATAN(1.0E0)
         EPMACH = R1MACH(3)
         EPS = 0.10E0*EPMACH
         NTJ0 = IWCS(BJ0CS, LBJ0, EPS)
         NTM0 = IWCS(BM0CS, LBM0, EPS)
         NTM02 = IWCS(BM02CS, LBM02, EPS)
         NTTH0 = IWCS(BTH0CS, LBTH0, EPS)
         NTT02 = IWCS(BT02CS, LBT02, EPS)
         XSML = SQRT(4.0E0*EPMACH)
         XMAX = 1.0E0/R1MACH(4)
      ENDIF
C
      DO 10 I=1,M
         Y(I) = ABS(X(I))
  10  CONTINUE
C
      CALL WNGT(M,Y,XMAX,KEY)
      IF (KEY .NE. 0) GO TO 920
C
C  ----------------
C  CASE Y .LE. XSML
C  ----------------
C
      DO 15 I=1,M
         F(I) = 1.0E0
  15  CONTINUE
C
C  --------------------------
C  CASE  XSML .LT. Y .LE. 4.0
C  --------------------------
C
      CALL WGTLE(M,Y,XSML,4.0E0,N,INDX)
      IF (N .GT. 0) THEN
         CALL WGTHR(N,Y,INDX,YCMP)
         DO 20 J=1,N
            TCMP(J) = 0.1250E0*YCMP(J)**2 - 1.0E0
  20     CONTINUE
         CALL WCS(N,TCMP,BJ0CS,NTJ0,ZCMP,B0,B1,B2)
         CALL WSCTR(N,ZCMP,INDX,F)
      ENDIF
C
C  ----------------
C  CASE  Y .GT. 4.0
C  ----------------
C
      CALL WGT(M,Y,4.0E0,N,INDX)
      IF (N .GT. 0) THEN
         CALL WGTHR(N,Y,INDX,YCMP)
         DO 50 J=1,N
            TCMP(J) = 32.0E0/YCMP(J)**2 - 1.0E0
  50     CONTINUE
         CALL WCS(N,TCMP,BM0CS,NTM0,Y,B0,B1,B2)
         CALL WCS(N,TCMP,BTH0CS,NTTH0,ZCMP,B0,B1,B2)
         DO 60 J=1,N
            Y(J) = (0.750E0 + Y(J)) / SQRT(YCMP(J))
            ZCMP(J) = (YCMP(J) - PI4) + ZCMP(J) / YCMP(J)
            ZCMP(J) = Y(J) * COS(ZCMP(J))
  60     CONTINUE
         CALL WSCTR(N,ZCMP,INDX,F)
      ENDIF
C
C  -----
C  EXITS
C  -----
C
C     ... NORMAL
C
      INFO = 0
      GO TO 999
C
C     ... M .LE. 0
C
  910 CONTINUE
      INFO = 1
      GO TO 999
C
C     ... NO PRECISION BECAUSE X BIG
C
  920 CONTINUE
      INFO = 2
      INDX(1) = KEY
      GO TO 999
C
  999 CONTINUE
      RETURN
      END
      SUBROUTINE VJ1 (M, X, F, WORK, IWORK, INFO)
C***BEGIN PROLOGUE  VJ1
C***PURPOSE  Computes the Bessel function of the first kind
C            of order one (J1) for a vector of real arguments
C***LIBRARY   VFNLIB
C***CATEGORY  C10B1
C***TYPE      SINGLE PRECISION (VJ1-S, DVJ1-D)
C***KEYWORDS  BESSEL FUNCTION, FIRST KIND, SPECIAL FUNCTION,
C             ORDER ONE, VECTORIZED
C***AUTHOR  SAUNDERS, B. V., (NIST)
C             COMPUTING AND APPLIED MATHEMATICS LABORATORY
C             NATIONAL INSTITUTE OF STANDARDS AND TECHNOLOGY
C             GAITHERSBURG, MD 20899
C           BOISVERT, R. F., (NIST)
C             COMPUTING AND APPLIED MATHEMATICS LABORATORY
C             NATIONAL INSTITUTE OF STANDARDS AND TECHNOLOGY
C             GAITHERSBURG, MD 20899
C***DESCRIPTION
C
C   VJ1 computes the  Bessel function of the first kind of order
C   one (J1) for real arguments using uniform approximation by
C   Chebyshev polynomials.
C
C
C   P A R A M E T E R S
C
C   M      (Input) Integer (M .GT. 0)
C          The number of arguments at which the function is to be
C          evaluated.
C
C   X      (Input) Real array of length M
C          The arguments at which the function is to be evaluated are
C          stored in X(1) to X(M) in any order.
C
C   F      (Output) Real array of length M
C          F(i) contains the value of the function at X(i), i=1,..,M.
C
C   WORK   (Work) Real vector of length 7*M
C
C   IWORK  (Work) Integer vector of length M
C
C
C   INFO   (Output) Integer
C          Indicates status of computed result. The following table
C          lists possible values and their meanings.  If OK=Yes then
C          all F(i) have been set, otherwise none have been set.
C
C          INFO  OK            Description
C          ------------------------------------------------------------
C           -1   Yes  Warning: Some abs(X(i)) so small J1 underflows.
C                     The corresponding F(i) are set to zero.
C            0   Yes  Successfull execution.
C            1   No   Error: M .LE. 0
C            2   No   Error: Some abs(X(i)) so big that no precision
C                     possible in computing J1.  The index of the
C                     first offending argument is returned in IWORK(1).
C
C *********************************************************************
C   This routine is a modification of the function BESJ1 developed by
C   W. Fullerton of LANL.
C *********************************************************************
C
C***REFERENCES  (NONE)
C***ROUTINES CALLED  WJ1
C***REVISION HISTORY  (YYMMDD)
C   910226  DATE WRITTEN
C***END PROLOGUE  VJ1
C
C     ----------
C     PARAMETERS
C     ----------
C
      INTEGER INFO, IWORK, M
      REAL    F, X, WORK
C
      DIMENSION X(M), F(M), WORK(7*M), IWORK(M)
C
C     ---------------
C     LOCAL VARIABLES
C     ---------------
C
      INTEGER IWB0, IWB1, IWB2, IWY, IWTC, IWYC, IWZC, JIN
C
C***FIRST EXECUTABLE STATEMENT  VJ1
C
C     ... PARTITION WORK ARRAYS
C
      IWY   = 1
      IWTC  = IWY  + M
      IWYC  = IWTC + M
      IWZC  = IWYC + M
      IWB0  = IWZC + M
      IWB1  = IWB0 + M
      IWB2  = IWB1 + M
C     Total = IWB2 + M
C
      JIN   = 1
C     Total = JIN  + M
C
C     ... WJ1 DOES ALL THE WORK
C
      CALL WJ1(M,X,F,WORK(IWY),WORK(IWTC),WORK(IWYC),WORK(IWZC),
     +         IWORK(JIN),WORK(IWB0),WORK(IWB1),WORK(IWB2),INFO)
C
      RETURN
      END
      SUBROUTINE WJ1(M, X, F, Y, TCMP, YCMP, ZCMP, INDX, B0, B1, B2, 
     +   INFO)
C***BEGIN PROLOGUE  WJ1
C***SUBSIDIARY
C***PURPOSE  Computes the Bessel function of the first kind
C            of order one (J1) for a vector of arguments
C***LIBRARY   VFNLIB
C***AUTHOR  SAUNDERS, B. V., (NIST)
C             COMPUTING AND APPLIED MATHEMATICS LABORATORY
C             NATIONAL INSTITUTE OF STANDARDS AND TECHNOLOGY
C             GAITHERSBURG, MD 20899
C           BOISVERT, R. F., (NIST)
C             COMPUTING AND APPLIED MATHEMATICS LABORATORY
C             NATIONAL INSTITUTE OF STANDARDS AND TECHNOLOGY
C             GAITHERSBURG, MD 20899
C***ROUTINES CALLED  IWCS, R1MACH, WNGT, WGTHR, WGTLE, WGT,
C                    WLE, WSCTR, WCS
C***REVISION HISTORY  (YYMMDD)
C   910226  DATE WRITTEN
C***END PROLOGUE  WJ1
C
C  ----------
C  PARAMETERS
C  ----------
C
      INTEGER INFO, INDX, M
      REAL    B0, B1, B2, F, X, Y, TCMP, YCMP, ZCMP
C
      DIMENSION B0(M), B1(M), B2(M), F(M), INDX(M), X(M), Y(M),
     +          TCMP(M), YCMP(M), ZCMP(M)
C
C  ---------------
C  LOCAL VARIABLES
C  ---------------
C
      INTEGER LBJ1, LBM1, LBM12, LBTH1, LBT12
      PARAMETER (LBJ1=12, LBM1=21, LBM12=2, LBTH1=24, LBT12=2)
C
      INTEGER I, IWCS, J, JH, K, KEY, N, NA, NB, NTJ1, NTM1, NTM12,
     +        NTTH1, NTT12, NTOT
      REAL    BJ1CS, BM1CS, BM12CS, BTH1CS, BT12CS, C1, C2,
     +        EPMACH, EPS, PI4, R1MACH, TPI4, XMAX, XMIN, XSML
C
      DIMENSION BJ1CS(LBJ1), BM1CS(LBM1), BM12CS(LBM12), BTH1CS(LBTH1),
     +          BT12CS(LBT12)
C
      SAVE BJ1CS, BM1CS, BM12CS, BTH1CS, BT12CS, NTJ1, NTM1, NTM12,
     +      NTTH1, NTT12, PI4, TPI4, XMIN, XSML, XMAX
C
C----------------------------------------------------------------------
C
C Series for BJ1        on the interval  0.          to  1.60000D+01
C                                        with weighted error   4.48E-17
C                                         log weighted error  16.35
C                               significant figures required  15.77
C                                    decimal places required  16.89
C
      DATA BJ1 CS( 1) /   -.1172614151 3332787E0 /
      DATA BJ1 CS( 2) /   -.2536152183 0790640E0 /
      DATA BJ1 CS( 3) /    .0501270809 84469569E0 /
      DATA BJ1 CS( 4) /   -.0046315148 09625081E0 /
      DATA BJ1 CS( 5) /    .0002479962 29415914E0 /
      DATA BJ1 CS( 6) /   -.0000086789 48686278E0 /
      DATA BJ1 CS( 7) /    .0000002142 93917143E0 /
      DATA BJ1 CS( 8) /   -.0000000039 36093079E0 /
      DATA BJ1 CS( 9) /    .0000000000 55911823E0 /
      DATA BJ1 CS(10) /   -.0000000000 00632761E0 /
      DATA BJ1 CS(11) /    .0000000000 00005840E0 /
      DATA BJ1 CS(12) /   -.0000000000 00000044E0 /
C
C----------------------------------------------------------------------
C
C Series for BM1        on the interval  0.          to  6.25000D-02
C                                        with weighted error   5.61E-17
C                                         log weighted error  16.25
C                               significant figures required  14.97
C                                    decimal places required  16.91
C
      DATA BM1 CS( 1) /    .1047362510 931285E0 /
      DATA BM1 CS( 2) /    .0044244389 3702345E0 /
      DATA BM1 CS( 3) /   -.0000566163 9504035E0 /
      DATA BM1 CS( 4) /    .0000023134 9417339E0 /
      DATA BM1 CS( 5) /   -.0000001737 7182007E0 /
      DATA BM1 CS( 6) /    .0000000189 3209930E0 /
      DATA BM1 CS( 7) /   -.0000000026 5416023E0 /
      DATA BM1 CS( 8) /    .0000000004 4740209E0 /
      DATA BM1 CS( 9) /   -.0000000000 8691795E0 /
      DATA BM1 CS(10) /    .0000000000 1891492E0 /
      DATA BM1 CS(11) /   -.0000000000 0451884E0 /
      DATA BM1 CS(12) /    .0000000000 0116765E0 /
      DATA BM1 CS(13) /   -.0000000000 0032265E0 /
      DATA BM1 CS(14) /    .0000000000 0009450E0 /
      DATA BM1 CS(15) /   -.0000000000 0002913E0 /
      DATA BM1 CS(16) /    .0000000000 0000939E0 /
      DATA BM1 CS(17) /   -.0000000000 0000315E0 /
      DATA BM1 CS(18) /    .0000000000 0000109E0 /
      DATA BM1 CS(19) /   -.0000000000 0000039E0 /
      DATA BM1 CS(20) /    .0000000000 0000014E0 /
      DATA BM1 CS(21) /   -.0000000000 0000005E0 /
C
C----------------------------------------------------------------------
C
C Series for BM12       Used in double precision version only
C
      DATA BM12CS( 1) /   1.0E0 /
      DATA BM12CS( 2) /   0.0E0 /
C
C----------------------------------------------------------------------
C
C Series for BTH1       on the interval  0.          to  6.25000D-02
C                                        with weighted error   4.10E-17
C                                         log weighted error  16.39
C                               significant figures required  15.96
C
      DATA BTH1CS( 1) /    .7406014102 6313850E0 /
      DATA BTH1CS( 2) /   -.0045717556 59637690E0 /
      DATA BTH1CS( 3) /    .0001198185 10964326E0 /
      DATA BTH1CS( 4) /   -.0000069645 61891648E0 /
      DATA BTH1CS( 5) /    .0000006554 95621447E0 /
      DATA BTH1CS( 6) /   -.0000000840 66228945E0 /
      DATA BTH1CS( 7) /    .0000000133 76886564E0 /
      DATA BTH1CS( 8) /   -.0000000024 99565654E0 /
      DATA BTH1CS( 9) /    .0000000005 29495100E0 /
      DATA BTH1CS(10) /   -.0000000001 24135944E0 /
      DATA BTH1CS(11) /    .0000000000 31656485E0 /
      DATA BTH1CS(12) /   -.0000000000 08668640E0 /
      DATA BTH1CS(13) /    .0000000000 02523758E0 /
      DATA BTH1CS(14) /   -.0000000000 00775085E0 /
      DATA BTH1CS(15) /    .0000000000 00249527E0 /
      DATA BTH1CS(16) /   -.0000000000 00083773E0 /
      DATA BTH1CS(17) /    .0000000000 00029205E0 /
      DATA BTH1CS(18) /   -.0000000000 00010534E0 /
      DATA BTH1CS(19) /    .0000000000 00003919E0 /
      DATA BTH1CS(20) /   -.0000000000 00001500E0 /
      DATA BTH1CS(21) /    .0000000000 00000589E0 /
      DATA BTH1CS(22) /   -.0000000000 00000237E0 /
      DATA BTH1CS(23) /    .0000000000 00000097E0 /
      DATA BTH1CS(24) /   -.0000000000 00000040E0 /
C
C----------------------------------------------------------------------
C
C Series for BT12       Used in double precision version only
C
      DATA BT12CS( 1) /   1.0E0 /
      DATA BT12CS( 2) /   0.0E0 /
C
C----------------------------------------------------------------------
C
      DATA NTJ1 / 0 /
C
C***FIRST EXECUTABLE STATEMENT  WJ1
C
      IF (M .LE. 0) GO TO 910
C
      IF (NTJ1.EQ.0) THEN
         EPMACH = R1MACH(3)
         EPS = 0.10E0*EPMACH
         NTJ1  = IWCS(BJ1CS, LBJ1, EPS)
         NTM1  = IWCS(BM1CS, LBM1, EPS)
         NTM12 = IWCS(BM12CS, LBM12, EPS)
         NTTH1 = IWCS(BTH1CS, LBTH1, EPS)
         NTT12 = IWCS(BT12CS, LBT12, EPS)
         XMIN = 2.0E0*R1MACH(1)
         XSML = SQRT(8.0E0*EPMACH)
         XMAX = 1.0E0/R1MACH(4)
         PI4 = ATAN(1.0E0)
         TPI4 = 3.0E0*PI4
      ENDIF
C
      NTOT = 0
C
      DO 10 I=1,M
         Y(I) = ABS(X(I))
  10  CONTINUE
C
      CALL WNGT(M,Y,XMAX,KEY)
      IF (KEY .NE. 0) GO TO 920
C
C  ---------------------------
C  CASE   Y=0  OR  Y TOO SMALL
C  ---------------------------
C
C     NOTE --- J1 UNDERFLOWS FOR  0 .LT. Y .LE. XMIN
C
      DO 15 I=1,M
         F(I) = 0.0E0
  15  CONTINUE
C
C  ----------------------------
C  CASE   XMIN .LT. Y .LE. XSML
C  ----------------------------
C
      CALL WGTLE(M,Y,XMIN,XSML,N,INDX)
      IF (N .GT. 0) THEN
         NTOT = NTOT + N
         CALL WGTHR(N,Y,INDX,YCMP)
         DO 20 J=1,N
            ZCMP(J) = 0.50E0*YCMP(J)
  20     CONTINUE
         CALL WSCTR(N,ZCMP,INDX,F)
      ENDIF
C
C  ---------------------------
C  CASE   XSML .LT. Y .LE. 4.0
C  ---------------------------
C
      CALL WGTLE(M,Y,XSML,4.0E0,N,INDX)
      IF (N .GT. 0) THEN
         NTOT = NTOT + N
         CALL WGTHR(N,Y,INDX,YCMP)
         DO 30 J=1,N
            TCMP(J) = 0.1250E0*YCMP(J)**2 - 1.0E0
  30     CONTINUE
         CALL WCS(N,TCMP,BJ1CS,NTJ1,ZCMP,B0,B1,B2)
         DO 40 K=1,N
            ZCMP(K) = YCMP(K)*(0.250E0 + ZCMP(K))
  40     CONTINUE
         CALL WSCTR(N,ZCMP,INDX,F)
      ENDIF
C
C  ----------------
C  CASE  Y .GT. 4.0
C  ----------------
C
      CALL WGT(M,Y,4.0E0,N,INDX)
      IF (N .GT. 0) THEN
         NTOT = NTOT + N
         CALL WGTHR(N,Y,INDX,YCMP)
         DO 50 J=1,N
            TCMP(J) = 32.0E0/YCMP(J)**2 - 1.0E0
  50     CONTINUE
         CALL WCS(N,TCMP,BM1CS,NTM1,Y,B0,B1,B2)
         CALL WCS(N,TCMP,BTH1CS,NTTH1,ZCMP,B0,B1,B2)
         DO 60 J=1,N
            Y(J) = (0.750E0 + Y(J)) / SQRT(YCMP(J))
            ZCMP(J) = ( YCMP(J) - TPI4 ) + ZCMP(J)/YCMP(J)
            ZCMP(J) = Y(J) * COS(ZCMP(J))
  60     CONTINUE
         CALL WSCTR(N,ZCMP,INDX,F)
      ENDIF
C
C  -------------------------------------------
C  REVERSE SIGN FOR NEGATIVE X (RESULT IS ODD)
C  -------------------------------------------
C
      DO 200 I=1,M
         IF (X(I) .LT. 0.0E0)  F(I) = -F(I)
  200 CONTINUE
C
C  -----
C  EXITS
C  -----
C
C     ... NORMAL
C
      INFO = 0
      IF (NTOT .NE. M)  INFO = -1
      GO TO 999
C
C     ... M .LE. 0
C
  910 CONTINUE
      INFO = 1
      GO TO 999
C
C     ... NO PRECISION BECAUSE X BIG
C
  920 CONTINUE
      INFO = 2
      INDX(1) = KEY
      GO TO 999
C
  999 CONTINUE
      RETURN
      END
      SUBROUTINE VY0 (M, X, F, WORK, IWORK, INFO)
C***BEGIN PROLOGUE  VY0
C***PURPOSE  Computes the Bessel function of the second kind
C            of order zero (Y0) for a vector of real arguments
C***LIBRARY   VFNLIB
C***CATEGORY  C10A1
C***TYPE      SINGLE PRECISION (VY0-S, DVY0-D)
C***KEYWORDS  BESSEL FUNCTION,SECOND KIND, ORDER ZERO, VECTORIZED
C***AUTHOR  SAUNDERS, B. V., (NIST)
C             COMPUTING AND APPLIED MATHEMATICS LABORATORY
C             NATIONAL INSTITUTE OF STANDARDS AND TECHNOLOGY
C             GAITHERSBURG, MD 20899
C           BOISVERT, R. F., (NIST)
C             COMPUTING AND APPLIED MATHEMATICS LABORATORY
C             NATIONAL INSTITUTE OF STANDARDS AND TECHNOLOGY
C             GAITHERSBURG, MD 20899
C***DESCRIPTION
C
C   VY0 computes the  Bessel function of the second kind of order
C   zero (Y0) for real arguments using uniform approximation by
C   Chebyshev polynomials.
C
C
C   P A R A M E T E R S
C
C   M      (Input) Integer (M .GT. 0)
C          The number of arguments at which the function is to be
C          evaluated.
C
C   X      (Input) Real array of length M
C          The arguments at which the function is to be evaluated are
C          stored in X(1) to X(M) in any order.
C
C   F      (Output) Real array of length M
C          F(i) contains the value of the function at X(i), i=1,..,M.
C
C   WORK   (Work) Real vector of length 7*M
C
C   IWORK  (Work) Integer vector of length M
C
C   INFO   (Output) Integer
C          Indicates status of computed result. The following table
C          lists possible values and their meanings.  If OK=Yes then
C          all F(i) have been set, otherwise none have been set.
C
C          INFO  OK            Description
C          ------------------------------------------------------------
C            0   Yes  Successfull execution.
C            1   No   Error: M .LE. 0
C            2   No   Error: X(i) zero or negative for some i
C            3   No   Error: Some X(i) so big that no precision possible
C                     in computing Y0. The index of the first offending
C                     argument is returned in IWORK(1).
C
C
C *********************************************************************
C   This routine is a modification of the function BESY0 developed by
C   W. Fullerton of LANL.
C *********************************************************************
C
C***REFERENCES  (NONE)
C***ROUTINES CALLED  WY0
C***REVISION HISTORY  (YYMMDD)
C   910226  DATE WRITTEN
C***END PROLOGUE  VY0
C
C     ----------
C     PARAMETERS
C     ----------
C
      INTEGER INFO, IWORK, M
      REAL    F, X, WORK
C
      DIMENSION X(M), F(M), WORK(7*M), IWORK(M)
C
C     ---------------
C     LOCAL VARIABLES
C     ---------------
C
      INTEGER IWB0, IWB1, IWB2, IWTC, IWXC, IWYC, IWZC, JIN
C
C
C***FIRST EXECUTABLE STATEMENT  VY0
C
C     ... PARTITION WORK ARRAYS
C
      IWTC   = 1
      IWXC  = IWTC + M
      IWYC  = IWXC + M
      IWZC  = IWYC + M
      IWB0  = IWZC + M
      IWB1  = IWB0 + M
      IWB2  = IWB1 + M
C     Total = IB2  + M
C
      JIN   = 1
C     Total = JIN  + M
C
C     ... WY0 DOES ALL THE WORK
C
      CALL WY0(M,X,F,WORK(IWTC),WORK(IWXC),WORK(IWYC),WORK(IWZC),
     +         IWORK(JIN),WORK(IWB0),WORK(IWB1),WORK(IWB2),INFO)
C
      RETURN
      END

      SUBROUTINE WY0 (M, X, F, TCMP, XCMP, YCMP, ZCMP, INDX, B0, B1,  
     +   B2, INFO)
C***BEGIN PROLOGUE  WY0
C***SUBSIDIARY
C***PURPOSE  Computes the Bessel function of the second kind
C            of order zero (Y0) for a vector of arguments
C***LIBRARY   VFNLIB
C***AUTHOR  SAUNDERS, B. V., (NIST)
C             COMPUTING AND APPLIED MATHEMATICS LABORATORY
C             NATIONAL INSTITUTE OF STANDARDS AND TECHNOLOGY
C             GAITHERSBURG, MD 20899
C           BOISVERT, R. F., (NIST)
C             COMPUTING AND APPLIED MATHEMATICS LABORATORY
C             NATIONAL INSTITUTE OF STANDARDS AND TECHNOLOGY
C             GAITHERSBURG, MD 20899
C***ROUTINES CALLED  IWCS, R1MACH, WNGT, WGTHR, WGTLE, WGT, WLE,
C                    WSCTR, WCS
C***REVISION HISTORY  (YYMMDD)
C   910226  DATE WRITTEN
C***END PROLOGUE  WY0
C
C  ----------
C  PARAMETERS
C  ----------
C
      INTEGER INFO, INDX, M
      REAL    B0, B1, B2, F, X, TCMP, XCMP, YCMP, ZCMP
C
      DIMENSION B0(M), B1(M), B2(M), F(M), INDX(M), X(M), TCMP(M),
     +          XCMP(M), YCMP(M), ZCMP(M)
C
C  ---------------
C  LOCAL VARIABLES
C  ---------------
C
      INTEGER LBJ0, LBM0, LBM02, LBTH0, LBT02, LBY0
      PARAMETER (LBJ0=13, LBM0=21, LBM02=2, LBTH0=24, LBT02=2,
     +           LBY0=13)
C
      INTEGER IWCS, J, JH, KEY, N, NA, NB, NTJ0, NTM0, NTM02, NTTH0,
     +        NTT02, NTY0
      REAL    BJ0CS, BM0CS, BM02CS, BTH0CS, BT02CS, BY0CS, C1,
     +         C2, EPMACH, EPS, R1MACH, PI4, TWODPI, XSML, XMAX
C
      DIMENSION BY0CS(LBY0), BM0CS(LBM0), BM02CS(LBM02), BTH0CS(LBTH0),
     +          BT02CS(LBT02), BJ0CS(LBJ0)
C
      SAVE BJ0CS, BM0CS, BTH0CS, BY0CS, NTJ0, NTM0, NTTH0, NTY0,
     +     PI4, TWODPI, XSML, XMAX
C
C----------------------------------------------------------------------
C
C Series for BY0        on the interval  0.          to  1.60000D+01
C                                        with weighted error   1.20E-17
C                                         log weighted error  16.92
C                               significant figures required  16.15
C                                    decimal places required  17.48
C
      DATA BY0 CS( 1) /   -.0112778393 92865573E0 /
      DATA BY0 CS( 2) /   -.1283452375 6042035E0 /
      DATA BY0 CS( 3) /   -.1043788479 9794249E0 /
      DATA BY0 CS( 4) /    .0236627491 83969695E0 /
      DATA BY0 CS( 5) /   -.0020903916 47700486E0 /
      DATA BY0 CS( 6) /    .0001039754 53939057E0 /
      DATA BY0 CS( 7) /   -.0000033697 47162423E0 /
      DATA BY0 CS( 8) /    .0000000772 93842676E0 /
      DATA BY0 CS( 9) /   -.0000000013 24976772E0 /
      DATA BY0 CS(10) /    .0000000000 17648232E0 /
      DATA BY0 CS(11) /   -.0000000000 00188105E0 /
      DATA BY0 CS(12) /    .0000000000 00001641E0 /
      DATA BY0 CS(13) /   -.0000000000 00000011E0 /
C
C----------------------------------------------------------------------
C
C Series for BM0        on the interval  0.          to  6.25000D-02
C                                        with weighted error   4.98E-17
C                                         log weighted error  16.30
C                               significant figures required  14.97
C                                    decimal places required  16.96
C
      DATA BM0 CS( 1) /    .0928496163 7381644E0 /
      DATA BM0 CS( 2) /   -.0014298770 7403484E0 /
      DATA BM0 CS( 3) /    .0000283057 9271257E0 /
      DATA BM0 CS( 4) /   -.0000014330 0611424E0 /
      DATA BM0 CS( 5) /    .0000001202 8628046E0 /
      DATA BM0 CS( 6) /   -.0000000139 7113013E0 /
      DATA BM0 CS( 7) /    .0000000020 4076188E0 /
      DATA BM0 CS( 8) /   -.0000000003 5399669E0 /
      DATA BM0 CS( 9) /    .0000000000 7024759E0 /
      DATA BM0 CS(10) /   -.0000000000 1554107E0 /
      DATA BM0 CS(11) /    .0000000000 0376226E0 /
      DATA BM0 CS(12) /   -.0000000000 0098282E0 /
      DATA BM0 CS(13) /    .0000000000 0027408E0 /
      DATA BM0 CS(14) /   -.0000000000 0008091E0 /
      DATA BM0 CS(15) /    .0000000000 0002511E0 /
      DATA BM0 CS(16) /   -.0000000000 0000814E0 /
      DATA BM0 CS(17) /    .0000000000 0000275E0 /
      DATA BM0 CS(18) /   -.0000000000 0000096E0 /
      DATA BM0 CS(19) /    .0000000000 0000034E0 /
      DATA BM0 CS(20) /   -.0000000000 0000012E0 /
      DATA BM0 CS(21) /    .0000000000 0000004E0 /
C
C----------------------------------------------------------------------
C
C Series for BM02       Used in double precision version only
C
      DATA BM02CS( 1) /   1.0E0 /
      DATA BM02CS( 2) /   0.0E0 /
C
C----------------------------------------------------------------------
C
C Series for BTH0       on the interval  0.          to  6.25000D-02
C                                        with weighted error   3.67E-17
C                                         log weighted error  16.44
C                               significant figures required  15.53
C                                    decimal places required  17.13
C
C
      DATA BTH0CS( 1) /   -.2463916377 4300119E0 /
      DATA BTH0CS( 2) /    .0017370983 07508963E0 /
      DATA BTH0CS( 3) /   -.0000621836 33402968E0 /
      DATA BTH0CS( 4) /    .0000043680 50165742E0 /
      DATA BTH0CS( 5) /   -.0000004560 93019869E0 /
      DATA BTH0CS( 6) /    .0000000621 97400101E0 /
      DATA BTH0CS( 7) /   -.0000000103 00442889E0 /
      DATA BTH0CS( 8) /    .0000000019 79526776E0 /
      DATA BTH0CS( 9) /   -.0000000004 28198396E0 /
      DATA BTH0CS(10) /    .0000000001 02035840E0 /
      DATA BTH0CS(11) /   -.0000000000 26363898E0 /
      DATA BTH0CS(12) /    .0000000000 07297935E0 /
      DATA BTH0CS(13) /   -.0000000000 02144188E0 /
      DATA BTH0CS(14) /    .0000000000 00663693E0 /
      DATA BTH0CS(15) /   -.0000000000 00215126E0 /
      DATA BTH0CS(16) /    .0000000000 00072659E0 /
      DATA BTH0CS(17) /   -.0000000000 00025465E0 /
      DATA BTH0CS(18) /    .0000000000 00009229E0 /
      DATA BTH0CS(19) /   -.0000000000 00003448E0 /
      DATA BTH0CS(20) /    .0000000000 00001325E0 /
      DATA BTH0CS(21) /   -.0000000000 00000522E0 /
      DATA BTH0CS(22) /    .0000000000 00000210E0 /
      DATA BTH0CS(23) /   -.0000000000 00000087E0 /
      DATA BTH0CS(24) /    .0000000000 00000036E0 /
C
C----------------------------------------------------------------------
C
C Series for BT02       Used in double precision version only
C
      DATA BT02CS( 1) /   1.0E0 /
      DATA BT02CS( 2) /   0.0E0 /
C
C----------------------------------------------------------------------
C
C Series for BJ0        on the interval  0.          to  1.60000D+01
C                                        with weighted error   7.47E-18
C                                         log weighted error  17.13
C                               significant figures required  16.98
C                                    decimal places required  17.68
C
      DATA BJ0 CS( 1) /    .1002541619 68939137E0 /
      DATA BJ0 CS( 2) /   -.6652230077 64405132E0 /
      DATA BJ0 CS( 3) /    .2489837034 98281314E0 /
      DATA BJ0 CS( 4) /   -.0332527231 700357697E0 /
      DATA BJ0 CS( 5) /    .0023114179 304694015E0 /
      DATA BJ0 CS( 6) /   -.0000991127 741995080E0 /
      DATA BJ0 CS( 7) /    .0000028916 708643998E0 /
      DATA BJ0 CS( 8) /   -.0000000612 108586630E0 /
      DATA BJ0 CS( 9) /    .0000000009 838650793E0 /
      DATA BJ0 CS(10) /   -.0000000000 124235515E0 /
      DATA BJ0 CS(11) /    .0000000000 001265433E0 /
      DATA BJ0 CS(12) /   -.0000000000 000010619E0 /
      DATA BJ0 CS(13) /    .0000000000 000000074E0 /
C
C----------------------------------------------------------------------
C
      DATA NTY0   / 0 /
C
C***FIRST EXECUTABLE STATEMENT  WY0
C
      IF (M .LE. 0)  GOTO 910
C
      IF (NTY0 .EQ. 0) THEN
         PI4 = ATAN(1.0E0)
         TWODPI = 0.50E0/PI4
         EPMACH = R1MACH(3)
         EPS = 0.10E0*EPMACH
         NTY0 = IWCS(BY0CS, LBY0, EPS)
         NTM0 = IWCS(BM0CS, LBM0, EPS)
         NTM02 = IWCS(BM02CS, LBM02, EPS)
         NTTH0 = IWCS(BTH0CS, LBTH0, EPS)
         NTT02 = IWCS(BT02CS, LBT02, EPS)
         NTJ0 = IWCS(BJ0CS, LBJ0, EPS)
         XSML = SQRT (4.0E0*EPMACH)
         XMAX = 1.0E0/R1MACH(4)
      ENDIF
C
      CALL WNLE(M,X,0.0E0,KEY)
      IF (KEY .NE. 0) GO TO 920
C
      CALL WNGT(M,X,XMAX,KEY)
      IF (KEY .NE. 0) GO TO 930
C
C  ----------------
C  CASE  X .LE. 4.0
C  ----------------
C
      CALL WLE(M,X,4.0E0,N,INDX)
      IF (N .GT. 0) THEN
         CALL WGTHR(N,X,INDX,XCMP)
C
C        ... COMPUTE J0(X) ... RESULT IN ZCMP
C
         DO 10 J=1,N
            TCMP(J) = 0.125E0*XCMP(J)**2 - 1.0E0
   10    CONTINUE
         CALL WCS(N,TCMP,BJ0CS,NTJ0,ZCMP,B0,B1,B2)
C
         CALL WCS(N,TCMP,BY0CS,NTY0,YCMP,B0,B1,B2)
         DO 30 J=1,N
            YCMP(J) = TWODPI*LOG(0.50E0*XCMP(J))*ZCMP(J)
     +                 + 0.375E0 + YCMP(J)
   30    CONTINUE
         CALL WSCTR(N,YCMP,INDX,F)
      ENDIF
C
C  ----------------
C  CASE  X .GT. 4.0
C  ----------------
C
      CALL WGT(M,X,4.0E0,N,INDX)
      IF (N .GT. 0) THEN
         CALL WGTHR(N,X,INDX,XCMP)
         DO 50 J=1,N
            TCMP(J) = 32.0E0/XCMP(J)**2 - 1.0E0
  50     CONTINUE
         CALL WCS(N,TCMP,BTH0CS,NTTH0,ZCMP,B0,B1,B2)
         CALL WCS(N,TCMP,BM0CS,NTM0,YCMP,B0,B1,B2)
         DO 60 J=1,N
            YCMP(J) = (0.750E0 + YCMP(J)) / SQRT(XCMP(J))
            ZCMP(J) = (XCMP(J) - PI4) + ZCMP(J) / XCMP(J)
            YCMP(J) = YCMP(J) * SIN(ZCMP(J))
  60     CONTINUE
         CALL WSCTR(N,YCMP,INDX,F)
      ENDIF
C
C  -----
C  EXITS
C  -----
C
C     ... NORMAL
C
      INFO = 0
      GO TO 999
C
C     ... M .LE. 0
C
  910 CONTINUE
      INFO = 1
      GO TO 999
C
C     ... X I ZERO OR NEGATIVE
C
  920 CONTINUE
      INFO = 2
      INDX(1) = KEY
      GO TO 999
C
C     ... NO PRECISION BECAUSE X BIG
C
  930 CONTINUE
      INFO = 3
      INDX(1) = KEY
      GO TO 999
C
  999 CONTINUE
      RETURN
      END
      SUBROUTINE VY1 (M, X, F, WORK, IWORK, INFO)
C***BEGIN PROLOGUE  VY1
C***PURPOSE  Computes the Bessel function of the second kind
C            of order one (Y1) for a vector of real arguments
C***LIBRARY   VFNLIB
C***CATEGORY  C10A1
C***TYPE      SINGLE PRECISION (VY1-S, DVY1-D)
C***KEYWORDS  BESSEL FUNCTION,SECOND KIND, ORDER ONE, VECTORIZED
C***AUTHOR  SAUNDERS, B. V., (NIST)
C             COMPUTING AND APPLIED MATHEMATICS LABORATORY
C             NATIONAL INSTITUTE OF STANDARDS AND TECHNOLOGY
C             GAITHERSBURG, MD 20899
C           BOISVERT, R. F., (NIST)
C             COMPUTING AND APPLIED MATHEMATICS LABORATORY
C             NATIONAL INSTITUTE OF STANDARDS AND TECHNOLOGY
C             GAITHERSBURG, MD 20899
C***DESCRIPTION
C
C   VY1 computes the  Bessel function of the second kind of order
C   one (Y1) for real arguments using uniform approximation by
C   Chebyshev polynomials.
C
C
C   P A R A M E T E R S
C
C   M      (Input) Integer (M .GT. 0)
C          The number of arguments at which the function is to be
C          evaluated.
C
C   X      (Input) Real array of length M
C          The arguments at which the function is to be evaluated are
C          stored in X(1) to X(M) in any order.
C
C   F      (Output) Real array of length M
C          F(i) contains the value of the function at X(i), i=1,..,M.
C
C   WORK   (Work) Real vector of length 7*M
C
C   IWORK  (Work) Integer vector of length M
C
C   INFO   (Output) Integer
C          Indicates status of computed result. The following table
C          lists possible values and their meanings.  If OK=Yes then
C          all F(i) have been set, otherwise none have been set.
C
C          INFO  OK            Description
C          ------------------------------------------------------------
C            0   Yes  Successfull execution.
C            1   No   Error: M .LE. 0
C            2   No   Error: X(i) zero or negative for some i
C            3   No   Error: Some X(i) so small Y1 overflows.
C                     The index of the first offending argument is
C                     returned in IWORK(1).
C
C
C *********************************************************************
C   This routine is a modification of the function BESY1 developed by
C   W. Fullerton of LANL.
C *********************************************************************
C
C***REFERENCES  (NONE)
C***ROUTINES CALLED  WY1
C***REVISION HISTORY  (YYMMDD)
C   910226  DATE WRITTEN
C***END PROLOGUE  VY1
C
C  ----------
C  PARAMETERS
C  ----------
      INTEGER INFO, IWORK, M
      REAL    F, X, WORK
C
      DIMENSION X(M), F(M), WORK(7*M), IWORK(M)
C
C ---------------
C  LOCAL VARIABLES
C  ---------------
C
      INTEGER IWB0, IWB1, IWB2, IWTC, IWXC, IWYC, IWZC, JIN
C
C
C***FIRST EXECUTABLE STATEMENT  VY1
C
C     ... PARTITION WORK ARRAYS
C
      IWTC   = 1
      IWXC  = IWTC + M
      IWYC  = IWXC + M
      IWZC  = IWYC + M
      IWB0  = IWZC + M
      IWB1  = IWB0 + M
      IWB2  = IWB1 + M
C     Total = IB2  + M
C
      JIN   = 1
C     Total = JIN  + M
C
C     ... WY1 DOES ALL THE WORK
C
      CALL WY1(M,X,F,WORK(IWTC),WORK(IWXC),WORK(IWYC),WORK(IWZC),
     +         IWORK(JIN),WORK(IWB0),WORK(IWB1),WORK(IWB2),INFO)
C
      RETURN
      END

      SUBROUTINE WY1 (M, X, F, TCMP, XCMP, YCMP, ZCMP, INDX, B0, B1,  
     +   B2, INFO)
C***BEGIN PROLOGUE  WY1
C***SUBSIDIARY
C***PURPOSE  Computes the Bessel function of the second kind
C            of order one (Y1) for a vector of arguments
C***LIBRARY   VFNLIB
C***AUTHOR  SAUNDERS, B. V., (NIST)
C             COMPUTING AND APPLIED MATHEMATICS LABORATORY
C             NATIONAL INSTITUTE OF STANDARDS AND TECHNOLOGY
C             GAITHERSBURG, MD 20899
C           BOISVERT, R. F., (NIST)
C             COMPUTING AND APPLIED MATHEMATICS LABORATORY
C             NATIONAL INSTITUTE OF STANDARDS AND TECHNOLOGY
C             GAITHERSBURG, MD 20899
C***ROUTINES CALLED  IWCS, R1MACH, WNGT, WGTHR, WGTLE, WGT, WLE,
C                    WSCTR, WCS, WNLE
C***REVISION HISTORY  (YYMMDD)
C   910226  DATE WRITTEN
C***END PROLOGUE  WY1
C
C  ----------
C  PARAMETERS
C  ----------
C
      INTEGER INFO, INDX, M
      REAL    B0, B1, B2, F, X, TCMP, XCMP, YCMP, ZCMP
C
      DIMENSION B0(M), B1(M), B2(M), F(M), INDX(M), X(M), TCMP(M),
     +          XCMP(M), YCMP(M), ZCMP(M)
C
C  ---------------
C  LOCAL VARIABLES
C  ---------------
C
      INTEGER LBJ1, LBM1, LBM12, LBTH1, LBT12, LBY1
      PARAMETER (LBJ1=12, LBM1=21, LBM12=2, LBTH1=24, LBT12=2,
     +           LBY1=14)
C
      INTEGER IWCS, J, JH, KEY, N, NA, NB, NTJ1, NTM1, NTM12, NTTH1,
     +        NTT12, NTY1
      REAL    BJ1CS, BM1CS, BM12CS, BTH1CS, BT12CS, BY1CS,
     +        C1, C2, EPMACH, EPS, PI4, R1MACH, TPI4, TWODPI, XMIN,
     +        XSML, XMAX
C
      DIMENSION BJ1CS(LBJ1), BM1CS(LBM1), BM12CS(LBM12), BTH1CS(LBTH1),
     +          BT12CS(LBT12), BY1CS(LBY1)
C
      SAVE BJ1CS, BM1CS, BTH1CS, BY1CS, NTJ1, NTM1, NTM12, NTTH1, NTT12,
     +      NTY1, PI4, TPI4, TWODPI, XMIN, XSML, XMAX
C
C----------------------------------------------------------------------
C
C
C Series for BY1        on the interval  0.          to  1.60000D+01
C                                        with weighted error   1.87E-18
C                                         log weighted error  17.73
C                               significant figures required  17.83
C                                    decimal places required  18.30
C
      DATA BY1 CS( 1) /    .0320804710 0611908629E0 /
      DATA BY1 CS( 2) /   1.2627078974 33500450E0 /
      DATA BY1 CS( 3) /    .0064999618 9992317500E0 /
      DATA BY1 CS( 4) /   -.0893616452 8860504117E0 /
      DATA BY1 CS( 5) /    .0132508812 2175709545E0 /
      DATA BY1 CS( 6) /   -.0008979059 1196483523E0 /
      DATA BY1 CS( 7) /    .0000364736 1487958306E0 /
      DATA BY1 CS( 8) /   -.0000010013 7438166600E0 /
      DATA BY1 CS( 9) /    .0000000199 4539657390E0 /
      DATA BY1 CS(10) /   -.0000000003 0230656018E0 /
      DATA BY1 CS(11) /    .0000000000 0360987815E0 /
      DATA BY1 CS(12) /   -.0000000000 0003487488E0 /
      DATA BY1 CS(13) /    .0000000000 0000027838E0 /
      DATA BY1 CS(14) /   -.0000000000 0000000186E0 /
C
C----------------------------------------------------------------------
C
C Series for BM1        on the interval  0.          to  6.25000D-02
C                                        with weighted error   5.61E-17
C                                         log weighted error  16.25
C                               significant figures required  14.97
C                                    decimal places required  16.91
C
      DATA BM1 CS( 1) /    .1047362510 931285E0 /
      DATA BM1 CS( 2) /    .0044244389 3702345E0 /
      DATA BM1 CS( 3) /   -.0000566163 9504035E0 /
      DATA BM1 CS( 4) /    .0000023134 9417339E0 /
      DATA BM1 CS( 5) /   -.0000001737 7182007E0 /
      DATA BM1 CS( 6) /    .0000000189 3209930E0 /
      DATA BM1 CS( 7) /   -.0000000026 5416023E0 /
      DATA BM1 CS( 8) /    .0000000004 4740209E0 /
      DATA BM1 CS( 9) /   -.0000000000 8691795E0 /
      DATA BM1 CS(10) /    .0000000000 1891492E0 /
      DATA BM1 CS(11) /   -.0000000000 0451884E0 /
      DATA BM1 CS(12) /    .0000000000 0116765E0 /
      DATA BM1 CS(13) /   -.0000000000 0032265E0 /
      DATA BM1 CS(14) /    .0000000000 0009450E0 /
      DATA BM1 CS(15) /   -.0000000000 0002913E0 /
      DATA BM1 CS(16) /    .0000000000 0000939E0 /
      DATA BM1 CS(17) /   -.0000000000 0000315E0 /
      DATA BM1 CS(18) /    .0000000000 0000109E0 /
      DATA BM1 CS(19) /   -.0000000000 0000039E0 /
      DATA BM1 CS(20) /    .0000000000 0000014E0 /
      DATA BM1 CS(21) /   -.0000000000 0000005E0 /
C
C----------------------------------------------------------------------
C
C Series for BM12       Used in double precision version only
C
      DATA BM12CS( 1) /   1.0E0 /
      DATA BM12CS( 2) /   0.0E0 /
C
C----------------------------------------------------------------------
C
C Series for BTH1       on the interval  0.          to  6.25000D-02
C                                        with weighted error   4.10E-17
C                                         log weighted error  16.39
C                               significant figures required  15.96
C                                    decimal places required  17.08
C
      DATA BTH1CS( 1) /    .7406014102 6313850E0 /
      DATA BTH1CS( 2) /   -.0045717556 59637690E0 /
      DATA BTH1CS( 3) /    .0001198185 10964326E0 /
      DATA BTH1CS( 4) /   -.0000069645 61891648E0 /
      DATA BTH1CS( 5) /    .0000006554 95621447E0 /
      DATA BTH1CS( 6) /   -.0000000840 66228945E0 /
      DATA BTH1CS( 7) /    .0000000133 76886564E0 /
      DATA BTH1CS( 8) /   -.0000000024 99565654E0 /
      DATA BTH1CS( 9) /    .0000000005 29495100E0 /
      DATA BTH1CS(10) /   -.0000000001 24135944E0 /
      DATA BTH1CS(11) /    .0000000000 31656485E0 /
      DATA BTH1CS(12) /   -.0000000000 08668640E0 /
      DATA BTH1CS(13) /    .0000000000 02523758E0 /
      DATA BTH1CS(14) /   -.0000000000 00775085E0 /
      DATA BTH1CS(15) /    .0000000000 00249527E0 /
      DATA BTH1CS(16) /   -.0000000000 00083773E0 /
      DATA BTH1CS(17) /    .0000000000 00029205E0 /
      DATA BTH1CS(18) /   -.0000000000 00010534E0 /
      DATA BTH1CS(19) /    .0000000000 00003919E0 /
      DATA BTH1CS(20) /   -.0000000000 00001500E0 /
      DATA BTH1CS(21) /    .0000000000 00000589E0 /
      DATA BTH1CS(22) /   -.0000000000 00000237E0 /
      DATA BTH1CS(23) /    .0000000000 00000097E0 /
      DATA BTH1CS(24) /   -.0000000000 00000040E0 /
C
C----------------------------------------------------------------------
C
C Series for BT12       Used in double precision version only
C
      DATA BT12CS( 1) /   1.0E0 /
      DATA BT12CS( 2) /   0.0E0 /
C
C----------------------------------------------------------------------
C
C Series for BJ1        on the interval  0.          to  1.60000D+01
C                                        with weighted error   4.48E-17
C                                         log weighted error  16.35
C                               significant figures required  15.77
C                                    decimal places required  16.89
C
      DATA BJ1 CS( 1) /   -.1172614151 3332787E0 /
      DATA BJ1 CS( 2) /   -.2536152183 0790640E0 /
      DATA BJ1 CS( 3) /    .0501270809 84469569E0 /
      DATA BJ1 CS( 4) /   -.0046315148 09625081E0 /
      DATA BJ1 CS( 5) /    .0002479962 29415914E0 /
      DATA BJ1 CS( 6) /   -.0000086789 48686278E0 /
      DATA BJ1 CS( 7) /    .0000002142 93917143E0 /
      DATA BJ1 CS( 8) /   -.0000000039 36093079E0 /
      DATA BJ1 CS( 9) /    .0000000000 55911823E0 /
      DATA BJ1 CS(10) /   -.0000000000 00632761E0 /
      DATA BJ1 CS(11) /    .0000000000 00005840E0 /
      DATA BJ1 CS(12) /   -.0000000000 00000044E0 /
C
C----------------------------------------------------------------------
C
      DATA NTY1 / 0 /
C
C***FIRST EXECUTABLE STATEMENT  WY1
C
      IF (M .LE. 0)  GO TO 910
C
      IF (NTY1 .EQ. 0) THEN
         EPMACH = R1MACH(3)
         EPS = 0.10E0*EPMACH
         NTY1 = IWCS(BY1CS, LBY1, EPS)
         NTM1 = IWCS(BM1CS, LBM1, EPS)
         NTM12 = IWCS(BM12CS, LBM12, EPS)
         NTTH1 = IWCS(BTH1CS, LBTH1, EPS)
         NTT12 = IWCS(BT12CS, LBT12, EPS)
         NTJ1 = IWCS(BJ1CS, LBJ1, EPS)
         XMIN = EXP(MAX( LOG(R1MACH(1)), -LOG(R1MACH(2))) + 0.010E0)
         XSML = SQRT(4.0E0*EPMACH)
         XMAX = 1.0/R1MACH(4)
         PI4 = ATAN(1.0E0)
         TPI4 = 3.0E0*PI4
         TWODPI = 0.50E0/PI4
      ENDIF
C
      CALL WNLE(M,X,0.0E0,KEY)
      IF (KEY .NE. 0) GO TO 920
C
      CALL WNLE(M,X,XMIN,KEY)
      IF (KEY .NE. 0) GO TO 930
C
C  ----------------
C  CASE  X .LE. 4.0
C  ----------------
C
      CALL WLE(M,X,4.0E0,N,INDX)
      IF (N .GT. 0) THEN
         CALL WGTHR(N,X,INDX,XCMP)
C
C        ... COMPUTE J1(X) ... RESULT IN ZCMP
C
         DO 20 J=1,N
            TCMP(J) = 0.125E0*XCMP(J)**2 - 1.0E0
   20    CONTINUE
         CALL WCS(N,TCMP,BJ1CS,NTJ1,ZCMP,B0,B1,B2)
         DO 30 J=1,N
            ZCMP(J) = XCMP(J)*(0.250E0 + ZCMP(J))
   30    CONTINUE
C
         CALL WCS(N,TCMP,BY1CS,NTY1,YCMP,B0,B1,B2)
         DO 40 J=1,N
            ZCMP(J) = TWODPI*LOG(0.50E0*XCMP(J))*ZCMP(J) +
     +                (0.50E0 + YCMP(J))/XCMP(J)
   40    CONTINUE
         CALL WSCTR(N,ZCMP,INDX,F)
      ENDIF
C
C  ----------------
C  CASE  X .GT. 4.0
C  ----------------
C
      CALL WGT(M,X,4.0E0,N,INDX)
      IF (N .GT. 0) THEN
         CALL WGTHR(N,X,INDX,XCMP)
         DO 50 J=1,N
            TCMP(J) = 32.0E0/XCMP(J)**2 - 1.0E0
  50     CONTINUE
         CALL WCS(N,TCMP,BTH1CS,NTTH1,ZCMP,B0,B1,B2)
         CALL WCS(N,TCMP,BM1CS,NTM1,YCMP,B0,B1,B2)
         DO 60 J=1,N
            YCMP(J) = (0.750E0 + YCMP(J)) / SQRT(XCMP(J))
            ZCMP(J) = (XCMP(J) - TPI4) + ZCMP(J) / XCMP(J)
            YCMP(J) = YCMP(J) * SIN(ZCMP(J))
  60     CONTINUE
         CALL WSCTR(N,YCMP,INDX,F)
      ENDIF
C
C  -----
C  EXITS
C  -----
C
C     ... NORMAL
C
      INFO = 0
      GO TO 999
C
C     ... M .LE. 0
C
  910 CONTINUE
      INFO = 1
      GO TO 999
C
C     ... X IS ZERO OR NEGATIVE
C
  920 CONTINUE
      INFO = 2
      INDX(1) = KEY
      GO TO 999
C
C     ... X SO SMALL THAT Y1 OVERFLOWS
C
  930 CONTINUE
      INFO = 3
      INDX(1) = KEY
      GO TO 999
C
  999 CONTINUE
      RETURN
      END
      SUBROUTINE WCS (M, X, CS, N, F, B0, B1, B2)
C***BEGIN PROLOGUE  WCS
C***PURPOSE  Evaluate a Chebyshev series for a vector of real arguments
C***LIBRARY   VFNLIB
C***CATEGORY  C3A2
C***TYPE      SINGLE PRECISION (WCS-S, DWCS-D)
C***KEYWORDS  CHEBYSHEV SERIES EVALUATION, VECTORIZED
C***AUTHOR  SAUNDERS, B. V,. (NIST)
C             COMPUTING AND APPLIED MATHEMATICS LABORATORY
C             NATIONAL INSTITUTE OF STANDARDS AND TECHNOLOGY
C             GAITHERSBURG, MD 20899
C           BOISVERT, R. F., (NIST)
C             COMPUTING AND APPLIED MATHEMATICS LABORATORY
C             NATIONAL INSTITUTE OF STANDARDS AND TECHNOLOGY
C             GAITHERSBURG, MD 20899
C***DESCRIPTION
C
C   WCS evaluates a given Chebyshev series for a vector of real
C   arguments.
C
C
C   P A R A M E T E R S
C
C   M   (Input) Integer
C       The number of arguments at which the series is to be
C       evaluated.
C
C   X   (Input) Real array of length M
C       The arguments at which the function is to be evaluated are
C       stored in X(1) to X(M) in any order.
C
C   CS  (Input) Real array of length .GE. N
C       The N coefficients of the Chebyshev series. (Note -- only half
C       the first coefficient is used in summing the series.)
C
C   N   (Input) Integer (0 .LE. N .LE. 1000)
C       The number of terms in the Chebyshev series.
C
C   F   (Output) Real array of length M
C       F(i) contains the value of the series at X(i), i=1,..,M.
C
C   B0  (Work) Real array of length M
C
C   B1  (Work) Real array of length M
C
C   B2  (Work) Real array of length M
C
C
C *********************************************************************
C   This routine is a modification of the function CSEVL  developed by
C   W. Fullerton of LANL.
C *********************************************************************
C
C***REFERENCES  (NONE)
C***ROUTINES CALLED  WFERR
C***REVISION HISTORY  (YYMMDD)
C   910226  DATE WRITTEN
C***END PROLOGUE  WCS
C
C     ----------
C     PARAMETERS
C     ----------
C
      INTEGER N, M
      REAL    X, CS, F, B0, B1, B2
C
      DIMENSION B0(M), B1(M), B2(M), CS(N), F(M), X(M)
C
C     ---------------
C     LOCAL VARIABLES
C     ---------------
C
      INTEGER I, IMOD, J, NI
      REAL    CS1, CS2, CSNI, CSNI1, CSNI2
C
C***FIRST EXECUTABLE STATEMENT  WCS
C
      IF (N .LT.    0)  CALL WFERR('WCS','NUMBER OF TERMS LT 0',2)
      IF (N .GT. 1000)  CALL WFERR('WCS','NUMBER OF TERMS GT 1000',2)
C
C     ... INITIALIZATION
C
      DO 10 I=1,M
         B1(I) = 0.0E0
         B2(I) = 0.0E0
         F(I)  = 2.0E0*X(I)
  10  CONTINUE
C
C     ... THREE-TERM RECURSION
C         (DO THREE STEPS AT A TIME TO AVOID VECTOR COPIES)
C
      IMOD = MOD(N,3)
      DO 30 I=1,N-IMOD,3
         NI = N + 1 - I
         CSNI = CS(NI)
         CSNI1 = CS(NI-1)
         CSNI2 = CS(NI-2)
         DO 20 J=1,M
            B0(J) = ( F(J)*B1(J)-B2(J) ) + CSNI
            B2(J) = ( F(J)*B0(J)-B1(J) ) + CSNI1
            B1(J) = ( F(J)*B2(J)-B0(J) ) + CSNI2
  20     CONTINUE
  30  CONTINUE
C
C     ... LAST STEP
C         (CLEANUP FOR CASE N NOT DIVISIBLE BY 3)
C
      IF (IMOD .EQ. 0) THEN
         DO 40 J=1,M
            F(J) = 0.50E0*(B1(J) - B0(J))
  40     CONTINUE
      ELSEIF (IMOD .EQ. 1) THEN
         CS1 = CS(1)
         DO 50 J=1,M
            B0(J) = ( F(J)*B1(J) - B2(J) ) + CS1
            F(J)  = 0.50E0*(B0(J) - B2(J))
  50     CONTINUE
      ELSE
         CS1 = CS(1)
         CS2 = CS(2)
         DO 60 J=1,M
            B0(J) = ( F(J)*B1(J) - B2(J) ) + CS2
            B2(J) = ( F(J)*B0(J) - B1(J) ) + CS1
            F(J)  = 0.50E0*(B2(J) - B1(J))
  60     CONTINUE
      ENDIF
C
      RETURN
      END
      INTEGER FUNCTION IWCS (OS, NOS, ETA)
C***BEGIN PROLOGUE  IWCS
C***PURPOSE  Determines the number of terms of an orthogonal series
C            required to meet a specified error tolerance.
C***LIBRARY   VFNLIB
C***CATEGORY  C3A2
C***TYPE      SINGLE PRECISION (IWCS-S, IDWCS-D)
C***KEYWORDS  CHEBYSHEV SERIES, INITIALIZATION
C***AUTHOR  BOISVERT, R. F., (NIST)
C             COMPUTING AND APPLIED MATHEMATICS LABORATORY
C             NATIONAL INSTITUTE OF STANDARDS AND TECHNOLOGY
C             GAITHERSBURG, MD 20899
C***DESCRIPTION
C
C   IWCS returns the number of terms of the Chebyshev series OS
C   required to insure the error in evaluating it is no larger than
C   ETA.
C
C   Ordinarily, ETA is chosen to be one-tenth machine precision.
C
C
C   P A R A M E T E R S
C
C   OS     (Input) Real array of length .GE. NOS
C          Coefficients of the NOS-term orthogonal series.
C
C   NOS    (Input) Integer (.GE. 1)
C          Number of terms in the orthogonal series.
C
C   ETA    (Input) Real
C          Requested accuracy of the series.
C
C
C *********************************************************************
C   This routine is a modification of the function INITS developed by
C   W. Fullerton of LANL.
C *********************************************************************
C
C***REFERENCES  (NONE)
C***ROUTINES CALLED  WFERR
C***REVISION HISTORY  (YYMMDD)
C   910226  DATE WRITTEN
C***END PROLOGUE  IWCS
C
C     ----------
C     PARAMETERS
C     ----------
C
      INTEGER NOS
      REAL    ETA, OS
C
      DIMENSION OS(NOS)
C
C     ---------------
C     LOCAL VARIABLES
C     ---------------
C
      INTEGER I, II
      REAL    ERR
C
C***FIRST EXECUTABLE STATEMENT  IWCS
C
      IF (NOS .LT. 1)
     +   CALL WFERR('IWCS','NUMBER OF COEFFICIENTS LT 1',2)
C
      ERR = 0.0E0
      DO 10 II=NOS,1,-1
         I = II
         ERR = ERR + ABS(OS(I))
         IF (ERR .GT. ETA) GO TO 20
 10   CONTINUE
C
 20   CONTINUE
      IF (I .EQ. NOS)
     +   CALL WFERR('IWCS','TOO MUCH ACCURACY REQUESTED',2)
C
      IWCS = I
C
      RETURN
      END
      SUBROUTINE WNLE (M, X, SCALR, KEY)
C***BEGIN PROLOGUE  WNLE
C***PURPOSE  Determines whether elements of a given vector are less
C            than or equal to a given scalar
C***LIBRARY   VFNLIB
C***TYPE      SINGLE PRECISION (WNLE-S, DWNLE-D)
C***AUTHOR  SAUNDERS, B. V., (NIST)
C             COMPUTING AND APPLIED MATHEMATICS LABORATORY
C             NATIONAL INSTITUTE OF STANDARDS AND TECHNOLOGY
C             GAITHERSBURG, MD 20899
C           BOISVERT, R. F., (NIST)
C             COMPUTING AND APPLIED MATHEMATICS LABORATORY
C             NATIONAL INSTITUTE OF STANDARDS AND TECHNOLOGY
C             GAITHERSBURG, MD 20899
C***DESCRIPTION
C
C   WNLE determines whether elements of a given vector are less than
C   or equal to a given scalar.
C
C
C   P A R A M E T E R S
C
C   M      (Input) Integer (M .GT. 0)
C          The number of elements of the input vector.
C
C   X      (Input) Real array of length M
C          The input vector.
C
C   SCALR  (Input) Real
C          The scalar used for comparison with vector elements.
C
C   KEY    (Output) Integer
C          If KEY=0 then no vector elements satisfy X(i).LE.SCALR.
C          If KEY>0 then some vector elements satisfy X(i).LE.SCALR;
C          the first to do so is X(KEY).
C
C
C   CAUTION : constraints on the input variables are not checked by
C             this routine.
C
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   910226  DATE WRITTEN
C***END PROLOGUE  WNLE
C
C     ----------
C     PARAMETERS
C     ----------
C
      INTEGER KEY, M
      REAL    SCALR, X
C
      DIMENSION X(M)
C
C     ---------------
C     LOCAL VARIABLES
C     ---------------
C
      INTEGER I
C
C***FIRST EXECUTABLE STATEMENT  WNLE
C
C     ... QUICK CHECK (VECTORIZABLE)
C
      KEY = 0
      DO 10 I=1,M
         IF (X(I) .LE. SCALR) KEY = 1
   10 CONTINUE
C
C     ... IF CHECK FAILED THEN FIND INDEX OF FIRST FAILURE
C
      IF (KEY .NE. 0) THEN
         KEY = 0
         DO 20 I=1,M
            IF (X(I) .LE. SCALR) THEN
               KEY = I
               GO TO 30
            ENDIF
   20    CONTINUE
      ENDIF
C
   30 CONTINUE
      RETURN
      END

      SUBROUTINE WNGT (M, X, SCALR, KEY)
C***BEGIN PROLOGUE  WNGT
C***PURPOSE  Determines whether elements of a given vector are greater
C            than a given scalar
C***LIBRARY   VFNLIB
C***TYPE      SINGLE PRECISION (WNGT-S, DWNGT-D)
C***AUTHOR  SAUNDERS, B. V., (NIST)
C             COMPUTING AND APPLIED MATHEMATICS LABORATORY
C             NATIONAL INSTITUTE OF STANDARDS AND TECHNOLOGY
C             GAITHERSBURG, MD 20899
C           BOISVERT, R. F., (NIST)
C             COMPUTING AND APPLIED MATHEMATICS LABORATORY
C             NATIONAL INSTITUTE OF STANDARDS AND TECHNOLOGY
C             GAITHERSBURG, MD 20899
C***DESCRIPTION
C
C   WNGT determines whether elements of a given vector are greater
C   than a given scalar.
C
C
C   P A R A M E T E R S
C
C   M      (Input) Integer (M .GT. 0)
C          The number of elements of the input vector.
C
C   X      (Input) Real array of length M
C          The input vector.
C
C   SCALR  (Input) Real
C          The scalar used for comparison with vector elements.
C
C   KEY    (Output) Integer
C          If KEY=0 then no vector elements satisfy X(i).GT.SCALR.
C          If KEY>0 then some vector elements satisfy X(i).GT.SCALR;
C          the first to do so is X(KEY).
C
C
C   CAUTION : constraints on the input variables are not checked by
C             this routine.
C
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   910226  DATE WRITTEN
C***END PROLOGUE  WNGT
C
C     ----------
C     PARAMETERS
C     ----------
C
      INTEGER KEY, M
      REAL    SCALR, X
C
      DIMENSION X(M)
C
C     ---------------
C     LOCAL VARIABLES
C     ---------------
C
      INTEGER I
C
C***FIRST EXECUTABLE STATEMENT  WNGT
C
C     ... QUICK CHECK (VECTORIZABLE)
C
      KEY = 0
      DO 10 I=1,M
         IF (X(I) .GT. SCALR) KEY = 1
   10 CONTINUE
C
C     ... IF CHECK FAILED THEN FIND INDEX OF FIRST FAILURE
C
      IF (KEY .NE. 0) THEN
         KEY = 0
         DO 20 I=1,M
            IF (X(I) .GT. SCALR) THEN
               KEY = I
               GO TO 30
            ENDIF
   20    CONTINUE
      ENDIF
C
   30 CONTINUE
      RETURN
      END

      SUBROUTINE WGTLE (M, X, SCALR1, SCALR2, N, INDX)
C***BEGIN PROLOGUE  WGTLE
C***PURPOSE  Builds an array of indices of vector elements that are
C            greater than a given scalar and less than or equal to a
C            second given scalar
C***LIBRARY   VFNLIB
C***TYPE      SINGLE PRECISION (WGTLE-S, DWGTLE-D)
C***AUTHOR  SAUNDERS, B. V., (NIST)
C             COMPUTING AND APPLIED MATHEMATICS LABORATORY
C             NATIONAL INSTITUTE OF STANDARDS AND TECHNOLOGY
C             GAITHERSBURG, MD 20899
C           BOISVERT, R. F., (NIST)
C             COMPUTING AND APPLIED MATHEMATICS LABORATORY
C             NATIONAL INSTITUTE OF STANDARDS AND TECHNOLOGY
C             GAITHERSBURG, MD 20899
C***DESCRIPTION
C
C   WGTLE builds an array of indices of vector elements that are
C   greater than a given scalar and less than or equal to a second
C   given scalar.
C
C
C   P A R A M E T E R S
C
C   M      (Input) Integer (M .GT. 0)
C          The number of elements of the input vector.
C
C   X      (Input) Real array of length M
C          The input vector.
C
C   SCALR1,
C   SCALR2 (Input) Real
C          The scalars used for comparison with vector elements.
C          Indices i are selected if SCALR1 .LT. X(i) .LE. SCALR2.
C
C   N      (Output) Integer
C          The number of elements of X that satisfy relationship with
C          scalars.
C
C   INDX   (Output) Integer array of length N
C          Array containing indices of vector elements that satisfy
C          relationship with scalars.
C
C
C   CAUTION : constraints on the input variables are not checked by
C             this routine.
C
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   910226  DATE WRITTEN
C***END PROLOGUE  WGTLE
C
C     ----------
C     PARAMETERS
C     ----------
C
      INTEGER INDX, M, N
      REAL    SCALR1, SCALR2, X
C
      DIMENSION X(M), INDX(*)
C
C     ---------------
C     LOCAL VARIABLES
C     ---------------
C
      INTEGER I
      REAL    ELEMT
C
C***FIRST EXECUTABLE STATEMENT WGTLE
C
      N = 0
C
      DO 10 I=1,M
         ELEMT = X(I)
         IF (ELEMT .GT. SCALR1 .AND. ELEMT .LE. SCALR2) THEN
            N = N + 1
            INDX(N) = I
         ENDIF
  10  CONTINUE
C
      RETURN
      END

      SUBROUTINE WGT (M, X, SCALR, N, INDX)
C***BEGIN PROLOGUE  WGT
C***PURPOSE  Builds an array of indices of vector elements that are
C            greater than a specified scalar
C***LIBRARY   VFNLIB
C***TYPE      SINGLE PRECISION (WGT-S, DWGT-D)
C***AUTHOR  SAUNDERS, B. V., (NIST)
C             COMPUTING AND APPLIED MATHEMATICS LABORATORY
C             NATIONAL INSTITUTE OF STANDARDS AND TECHNOLOGY
C             GAITHERSBURG, MD 20899
C           BOISVERT, R. F., (NIST)
C             COMPUTING AND APPLIED MATHEMATICS LABORATORY
C             NATIONAL INSTITUTE OF STANDARDS AND TECHNOLOGY
C             GAITHERSBURG, MD 20899
C***DESCRIPTION
C
C   WGT builds an array of indices of vector elements that are
C   greater than a specified scalar.
C
C
C   P A R A M E T E R S
C
C   M      (Input) Integer (M .GT. 0)
C          The number of elements of the input vector.
C
C   X      (Input) Real array of length M
C          The input vector.
C
C   SCALR  (Input) Real
C          The scalar used for comparison with vector elements.
C          Indices i are selected if X(i) .GT. SCALR.
C
C   N      (Output) Integer
C          The number of elements of X that satisfy relationship with
C          scalar.
C
C   INDX   (Output) Integer array of length N
C          Array containing indices of vector elements that satisfy
C          relationship with scalar.
C
C
C   CAUTION : constraints on the input variables are not checked by
C             this routine.
C
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   910226  DATE WRITTEN
C***END PROLOGUE  WGT
C
C     ----------
C     PARAMETERS
C     ----------
C
      INTEGER INDX, M, N
      REAL    SCALR, X
C
      DIMENSION X(M), INDX(*)
C
C     ---------------
C     LOCAL VARIABLES
C     ---------------
C
      INTEGER I
C
C***FIRST EXECUTABLE STATEMENT WGT
C
      N = 0
C
      DO 10 I=1,M
         IF (X(I) .GT. SCALR) THEN
            N = N + 1
            INDX(N) = I
         ENDIF
  10  CONTINUE
C
      RETURN
      END

      SUBROUTINE WLE (M, X, SCALR, N, INDX)
C***BEGIN PROLOGUE  WLE
C***PURPOSE  Builds an array of indices of vector elements that are
C            less than  or equal to a given scalar
C***LIBRARY   VFNLIB
C***TYPE      SINGLE PRECISION (WLE-S, DWLE-D)
C***AUTHOR  SAUNDERS, B. V., (NIST)
C             COMPUTING AND APPLIED MATHEMATICS LABORATORY
C             NATIONAL INSTITUTE OF STANDARDS AND TECHNOLOGY
C             GAITHERSBURG, MD 20899
C           BOISVERT, R. F., (NIST)
C             COMPUTING AND APPLIED MATHEMATICS LABORATORY
C             NATIONAL INSTITUTE OF STANDARDS AND TECHNOLOGY
C             GAITHERSBURG, MD 20899
C***DESCRIPTION
C
C   WLE builds an array of indices of vector elements that are
C   less than or equal to a given scalar.
C
C
C   P A R A M E T E R S
C
C   M      (Input) Integer (M .GT. 0)
C          The number of elements of the input vector.
C
C   X      (Input) Real array of length M
C          The input vector.
C
C   SCALR  (Input) Real
C          The scalar used for comparison with vector elements.
C          Indices i are selected if X(i) .LE. SCALR.
C
C   N      (Output) Integer
C          The number of elements of X that satisfy relationship with
C          scalars.
C
C   INDX   (Output) Integer array of length N
C          Array containing indices of vector elements that satisfy
C          relationship with scalars.
C
C
C   CAUTION : constraints on the input variables are not checked by
C             this routine.
C
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   910226  DATE WRITTEN
C***END PROLOGUE  WLE
C
C     ----------
C     PARAMETERS
C     ----------
C
      INTEGER INDX, M, N
      REAL    SCALR, X
C
      DIMENSION X(M), INDX(*)
C
C     ---------------
C     LOCAL VARIABLES
C     ---------------
C
      INTEGER I
C
C***FIRST EXECUTABLE STATEMENT WGE
C
      N = 0
C
      DO 10 I = 1,M
         IF (X(I) .LE. SCALR) THEN
            N = N + 1
            INDX(N) = I
         ENDIF
  10  CONTINUE
C
      RETURN
      END

      SUBROUTINE WGTHR (N, X, INDX, Y)
C***BEGIN PROLOGUE  WGTHR
C***PURPOSE  Performs a gather by index on a given vector
C***LIBRARY   VFNLIB
C***TYPE      SINGLE PRECISION (WGTHR-S, DWGTHR-D)
C***AUTHOR  SAUNDERS, B. V., (NIST)
C             COMPUTING AND APPLIED MATHEMATICS LABORATORY
C             NATIONAL INSTITUTE OF STANDARDS AND TECHNOLOGY
C             GAITHERSBURG, MD 20899
C           BOISVERT, R. F., (NIST)
C             COMPUTING AND APPLIED MATHEMATICS LABORATORY
C             NATIONAL INSTITUTE OF STANDARDS AND TECHNOLOGY
C             GAITHERSBURG, MD 20899
C***DESCRIPTION
C
C   WGTHR performs a gather by index on a given vector based on the
C   indices provided in an array: Y(i) = X(INDX(i)), i=1,...,N.
C
C
C   P A R A M E T E R S
C
C   N      (Input) Integer (N .GT. 0)
C          The number of elements to be gathered from the input vector
C          X.
C
C   X      (Input) Real array of length M
C          The input vector.
C
C   INDX   (Input) Integer array of length N
C          Array specifying indices of elements to be gathered from
C          the input vector.
C
C   Y      (Output) Real array of length N
C          The array containing the compressed vector once the gather is
C          completed.  Y(i) = X(INDX(i)), i=1,...,N.
C
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   910226  DATE WRITTEN
C***END PROLOGUE  WGTHR
C
C     ----------
C     PARAMETERS
C     ----------
C
      INTEGER INDX, N
      REAL    X, Y
C
      DIMENSION X(*), Y(N), INDX(N)
C
C     ---------------
C     LOCAL VARIABLES
C     ---------------
C
      INTEGER J
C
C***FIRST EXECUTABLE STATEMENT WGTHR
C
      IF (N .LE. 0) RETURN
C
      DO 10 J=1,N
         Y(J) = X(INDX(J))
  10  CONTINUE
C
      RETURN
      END

      SUBROUTINE WSCTR (N, Y, INDX, X)
C***BEGIN PROLOGUE  WSCTR
C***PURPOSE  Performs a scatter by index on a given vector
C***LIBRARY   VFNLIB
C***TYPE      SINGLE PRECISION (WSCTR-S, DWSCTR-D)
C***AUTHOR  SAUNDERS, B. V., (NIST)
C             COMPUTING AND APPLIED MATHEMATICS LABORATORY
C             NATIONAL INSTITUTE OF STANDARDS AND TECHNOLOGY
C             GAITHERSBURG, MD 20899
C           BOISVERT, R. F., (NIST)
C             COMPUTING AND APPLIED MATHEMATICS LABORATORY
C             NATIONAL INSTITUTE OF STANDARDS AND TECHNOLOGY
C             GAITHERSBURG, MD 20899
C***DESCRIPTION
C
C   WSCTR performs a scatter by index on a given vector based on the
C   indices provided in an array: X(INDX(i)) = Y(i), i=1,...,N.
C
C
C   P A R A M E T E R S
C
C   N      (Input) Integer (N .GT. 0)
C          The number of elements in the input vector Y.
C
C   Y      (Input) Real array of length N
C          The compressed vector.
C
C   INDX   (Input) Integer array of length N
C          Array of indices of positions elements of compressed
C          vector will occupy when scattered into vector X.
C
C   X      (Output) Real array
C          The vector to receive the scattered elements of Y.
C          X(INDX(i)) = Y(i), i=1,...,N.  Only elements of X
C          whose indices are in INDX are changed.
C
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   910226  DATE WRITTEN
C***END PROLOGUE  WSCTR
C
C
C     ----------
C     PARAMETERS
C     ----------
C
      INTEGER INDX, N
      REAL    Y, X
C
      DIMENSION Y(N), X(*), INDX(N)
C
C     ---------------
C     LOCAL VARIABLES
C     ---------------
C
      INTEGER J
C
C***FIRST EXECUTABLE STATEMENT WGTHR
C
      IF (N .LE. 0) RETURN
C
      DO 10 J=1,N
         X(INDX(J)) = Y(J)
  10  CONTINUE
C
      RETURN
      END
     


      double precision function dbesi0 (x)
c july 1977 edition.  w. fullerton, c3, los alamos scientific lab.
      double precision x, bi0cs(18), xmax, xsml, y, d1mach,
     1  dcsevl, dbsi0e, dexp, dlog, dsqrt
c     external d1mach, dbsi0e, dcsevl, dexp, dlog, dsqrt, initds
c
c series for bi0        on the interval  0.          to  9.00000e+00
c                                        with weighted error   9.51e-34
c                                         log weighted error  33.02
c                               significant figures required  33.31
c                                    decimal places required  33.65
c
      data bi0 cs(  1) / -.7660547252 8391449510 8189497624 3285 d-1   /
      data bi0 cs(  2) / +.1927337953 9938082699 5240875088 1196 d+1   /
      data bi0 cs(  3) / +.2282644586 9203013389 3702929233 0415 d+0   /
      data bi0 cs(  4) / +.1304891466 7072904280 7933421069 1888 d-1   /
      data bi0 cs(  5) / +.4344270900 8164874513 7868268102 6107 d-3   /
      data bi0 cs(  6) / +.9422657686 0019346639 2317174411 8766 d-5   /
      data bi0 cs(  7) / +.1434006289 5106910799 6209187817 9957 d-6   /
      data bi0 cs(  8) / +.1613849069 6617490699 1541971999 4611 d-8   /
      data bi0 cs(  9) / +.1396650044 5356696994 9509270814 2522 d-10  /
      data bi0 cs( 10) / +.9579451725 5054453446 2752317189 3333 d-13  /
      data bi0 cs( 11) / +.5333981859 8625021310 1510774400 0000 d-15  /
      data bi0 cs( 12) / +.2458716088 4374707746 9678591999 9999 d-17  /
      data bi0 cs( 13) / +.9535680890 2487700269 4434133333 3333 d-20  /
      data bi0 cs( 14) / +.3154382039 7214273367 8933333333 3333 d-22  /
      data bi0 cs( 15) / +.9004564101 0946374314 6666666666 6666 d-25  /
      data bi0 cs( 16) / +.2240647369 1236700160 0000000000 0000 d-27  /
      data bi0 cs( 17) / +.4903034603 2428373333 3333333333 3333 d-30  /
      data bi0 cs( 18) / +.9508172606 1226666666 6666666666 6666 d-33  /
c
      data nti0, xsml, xmax / 0, 2*0.d0 /
c
      if (nti0.ne.0) go to 10
      nti0 = initds (bi0cs, 18, 0.1*sngl(d1mach(3)))
      xsml = dsqrt (8.0d0*d1mach(3))
      xmax = dlog (d1mach(2))
c
 10   y = dabs(x)
      if (y.gt.3.0d0) go to 20
c
      dbesi0 = 1.0d0
      if (y.gt.xsml) dbesi0 = 2.75d0 + dcsevl (y*y/4.5d0-1.d0, bi0cs,
     1  nti0)
      return
c
 20   if (y.gt.xmax) call seteru (
     1  35hdbesi0  dabs(x) so big i0 overflows, 35, 2, 2)
c
      dbesi0 = dexp(y) * dbsi0e(x)
c
      return
      end
      double precision function dbesi1 (x)
c oct 1983 edition.  w. fullerton, c3, los alamos scientific lab.
      double precision x, bi1cs(17), xmax, xmin, xsml, y, d1mach,
     1  dcsevl, dbsi1e, dexp, dlog, dsqrt
c     external d1mach, dbsi1e, dcsevl, dexp, dlog, dsqrt, initds
c
c series for bi1        on the interval  0.          to  9.00000e+00
c                                        with weighted error   1.44e-32
c                                         log weighted error  31.84
c                               significant figures required  31.45
c                                    decimal places required  32.46
c
      data bi1 cs(  1) / -.1971713261 0998597316 1385032181 49 d-2     /
      data bi1 cs(  2) / +.4073488766 7546480608 1553936520 14 d+0     /
      data bi1 cs(  3) / +.3483899429 9959455866 2450377837 87 d-1     /
      data bi1 cs(  4) / +.1545394556 3001236038 5984010584 89 d-2     /
      data bi1 cs(  5) / +.4188852109 8377784129 4588320041 20 d-4     /
      data bi1 cs(  6) / +.7649026764 8362114741 9597039660 69 d-6     /
      data bi1 cs(  7) / +.1004249392 4741178689 1798080372 38 d-7     /
      data bi1 cs(  8) / +.9932207791 9238106481 3712980548 63 d-10    /
      data bi1 cs(  9) / +.7663801791 8447637275 2001716813 49 d-12    /
      data bi1 cs( 10) / +.4741418923 8167394980 3880919481 60 d-14    /
      data bi1 cs( 11) / +.2404114404 0745181799 8631720320 00 d-16    /
      data bi1 cs( 12) / +.1017150500 7093713649 1211007999 99 d-18    /
      data bi1 cs( 13) / +.3645093565 7866949458 4917333333 33 d-21    /
      data bi1 cs( 14) / +.1120574950 2562039344 8106666666 66 d-23    /
      data bi1 cs( 15) / +.2987544193 4468088832 0000000000 00 d-26    /
      data bi1 cs( 16) / +.6973231093 9194709333 3333333333 33 d-29    /
      data bi1 cs( 17) / +.1436794822 0620800000 0000000000 00 d-31    /
c
      data nti1, xmin, xsml, xmax / 0, 3*0.d0 /
c
      if (nti1.ne.0) go to 10
      nti1 = initds (bi1cs, 17, 0.1*sngl(d1mach(3)))
      xmin = 2.0d0*d1mach(1)
      xsml = dsqrt (8.0d0*d1mach(3))
      xmax = dlog (d1mach(2))
c
 10   y = dabs(x)
      if (y.gt.3.0d0) go to 20
c
      dbesi1 = 0.0d0
      if (y.eq.0.0d0) return
c
      if (y.le.xmin) call seteru (
     1  38hdbesi1  dabs(x) so small i1 underflows, 38, 1, 0)
      if (y.gt.xmin) dbesi1 = 0.5d0*x
      if (y.gt.xsml) dbesi1 = x*(0.875d0 + dcsevl (y*y/4.5d0-1.d0,
     1  bi1cs, nti1))
      return
c
 20   if (y.gt.xmax) call seteru (
     1  35hdbesi1  dabs(x) so big i1 overflows, 35, 2, 2)
c
      dbesi1 = dexp(y) * dbsi1e(x)
c
      return
      end
      double precision function dbesj0 (x)
c july 1977 edition.  w. fullerton, c3, los alamos scientific lab.
      double precision x, bj0cs(19), ampl, theta, xsml, y,
     1  d1mach, dcsevl, dcos, dsqrt
c     external d1mach, dcos, dcsevl, dsqrt, initds
c
c series for bj0        on the interval  0.          to  1.60000e+01
c                                        with weighted error   4.39e-32
c                                         log weighted error  31.36
c                               significant figures required  31.21
c                                    decimal places required  32.00
c
      data bj0 cs(  1) / +.1002541619 6893913701 0731272640 74 d+0     /
      data bj0 cs(  2) / -.6652230077 6440513177 6787578311 24 d+0     /
      data bj0 cs(  3) / +.2489837034 9828131370 4604687266 80 d+0     /
      data bj0 cs(  4) / -.3325272317 0035769653 8843415038 54 d-1     /
      data bj0 cs(  5) / +.2311417930 4694015462 9049241177 29 d-2     /
      data bj0 cs(  6) / -.9911277419 9508092339 0485193365 49 d-4     /
      data bj0 cs(  7) / +.2891670864 3998808884 7339037470 78 d-5     /
      data bj0 cs(  8) / -.6121085866 3032635057 8184074815 16 d-7     /
      data bj0 cs(  9) / +.9838650793 8567841324 7687486364 15 d-9     /
      data bj0 cs( 10) / -.1242355159 7301765145 5158970068 36 d-10    /
      data bj0 cs( 11) / +.1265433630 2559045797 9158272103 63 d-12    /
      data bj0 cs( 12) / -.1061945649 5287244546 9148175129 59 d-14    /
      data bj0 cs( 13) / +.7470621075 8024567437 0989155840 00 d-17    /
      data bj0 cs( 14) / -.4469703227 4412780547 6270079999 99 d-19    /
      data bj0 cs( 15) / +.2302428158 4337436200 5230933333 33 d-21    /
      data bj0 cs( 16) / -.1031914479 4166698148 5226666666 66 d-23    /
      data bj0 cs( 17) / +.4060817827 4873322700 8000000000 00 d-26    /
      data bj0 cs( 18) / -.1414383600 5240913919 9999999999 99 d-28    /
      data bj0 cs( 19) / +.4391090549 6698880000 0000000000 00 d-31    /
c
      data ntj0, xsml / 0, 0.d0 /
c
      if (ntj0.ne.0) go to 10
      ntj0 = initds (bj0cs, 19, 0.1*sngl(d1mach(3)))
      xsml = dsqrt (4.0d0*d1mach(3))
c
 10   y = dabs(x)
      if (y.gt.4.0d0) go to 20
c
      dbesj0 = 1.0d0
      if (y.gt.xsml) dbesj0 = dcsevl (.125d0*y*y-1.d0, bj0cs, ntj0)
      return
c
 20   call d9b0mp (y, ampl, theta)
      dbesj0 = ampl * dcos(theta)
c
      return
      end
      double precision function dbesj1 (x)
c sept 1983 edition.  w. fullerton, c3, los alamos scientific lab.
      double precision x, bj1cs(19), ampl, theta, xsml, xmin, y,
     1  d1mach, dcsevl, dcos, dsqrt
c     external d1mach, dcos, dcsevl, dsqrt, initds
c
c series for bj1        on the interval  0.          to  1.60000e+01
c                                        with weighted error   1.16e-33
c                                         log weighted error  32.93
c                               significant figures required  32.36
c                                    decimal places required  33.57
c
      data bj1 cs(  1) / -.1172614151 3332786560 6240574524 003 d+0    /
      data bj1 cs(  2) / -.2536152183 0790639562 3030884554 698 d+0    /
      data bj1 cs(  3) / +.5012708098 4469568505 3656363203 743 d-1    /
      data bj1 cs(  4) / -.4631514809 6250819184 2619728789 772 d-2    /
      data bj1 cs(  5) / +.2479962294 1591402453 9124064592 364 d-3    /
      data bj1 cs(  6) / -.8678948686 2788258452 1246435176 416 d-5    /
      data bj1 cs(  7) / +.2142939171 4379369150 2766250991 292 d-6    /
      data bj1 cs(  8) / -.3936093079 1831797922 9322764073 061 d-8    /
      data bj1 cs(  9) / +.5591182317 9468800401 8248059864 032 d-10   /
      data bj1 cs( 10) / -.6327616404 6613930247 7695274014 880 d-12   /
      data bj1 cs( 11) / +.5840991610 8572470032 6945563268 266 d-14   /
      data bj1 cs( 12) / -.4482533818 7012581903 9135059199 999 d-16   /
      data bj1 cs( 13) / +.2905384492 6250246630 6018688000 000 d-18   /
      data bj1 cs( 14) / -.1611732197 8414416541 2118186666 666 d-20   /
      data bj1 cs( 15) / +.7739478819 3927463729 8346666666 666 d-23   /
      data bj1 cs( 16) / -.3248693782 1119984114 3466666666 666 d-25   /
      data bj1 cs( 17) / +.1202237677 2274102272 0000000000 000 d-27   /
      data bj1 cs( 18) / -.3952012212 6513493333 3333333333 333 d-30   /
      data bj1 cs( 19) / +.1161678082 2664533333 3333333333 333 d-32   /
c
      data ntj1, xsml, xmin / 0, 2*0.d0 /
c
      if (ntj1.ne.0) go to 10
      ntj1 = initds (bj1cs, 19, 0.1*sngl(d1mach(3)))
c
      xsml = dsqrt (4.0d0*d1mach(3))
      xmin = 2.0d0*d1mach(1)
c
 10   y = dabs(x)
      if (y.gt.4.0d0) go to 20
c
      dbesj1 = 0.0d0
      if (y.eq.0.0d0) return
c
      if (y.le.xmin) call seteru (
     1  38hdbesj1  dabs(x) so small j1 underflows, 38, 1, 0)
      if (y.gt.xmin) dbesj1 = 0.5d0*x
      if (y.gt.xsml) dbesj1 = x*(.25d0 + dcsevl (.125d0*y*y-1.d0,
     1  bj1cs, ntj1) )
      return
c
 20   call d9b1mp (y, ampl, theta)
C**Next line changed from*****************
C     dbesj1 = ampl * dcos(theta)
C  to
      dbesj1 = sign(ampl,x) * dcos(theta)
C**by Ron Boisvert, 19 Mar 91 ************
c
      return
      end
      double precision function dbesk0 (x)
c july 1977 edition.  w. fullerton, c3, los alamos scientific lab.
      double precision x, bk0cs(16), xmax, xsml, y,
     1  d1mach, dcsevl, dbesi0, dbsk0e, dexp, dlog, dsqrt
c     external d1mach, dbesi0, dbsk0e, dcsevl, dexp, dlog, dsqrt, initds
c
c series for bk0        on the interval  0.          to  4.00000e+00
c                                        with weighted error   3.08e-33
c                                         log weighted error  32.51
c                               significant figures required  32.05
c                                    decimal places required  33.11
c
      data bk0 cs(  1) / -.3532739323 3902768720 1140060063 153 d-1    /
      data bk0 cs(  2) / +.3442898999 2462848688 6344927529 213 d+0    /
      data bk0 cs(  3) / +.3597993651 5361501626 5721303687 231 d-1    /
      data bk0 cs(  4) / +.1264615411 4469259233 8479508673 447 d-2    /
      data bk0 cs(  5) / +.2286212103 1194517860 8269830297 585 d-4    /
      data bk0 cs(  6) / +.2534791079 0261494573 0790013428 354 d-6    /
      data bk0 cs(  7) / +.1904516377 2202088589 7214059381 366 d-8    /
      data bk0 cs(  8) / +.1034969525 7633624585 1008317853 089 d-10   /
      data bk0 cs(  9) / +.4259816142 7910825765 2445327170 133 d-13   /
      data bk0 cs( 10) / +.1374465435 8807508969 4238325440 000 d-15   /
      data bk0 cs( 11) / +.3570896528 5083735909 9688597333 333 d-18   /
      data bk0 cs( 12) / +.7631643660 1164373766 7498666666 666 d-21   /
      data bk0 cs( 13) / +.1365424988 4407818590 8053333333 333 d-23   /
      data bk0 cs( 14) / +.2075275266 9066680831 9999999999 999 d-26   /
      data bk0 cs( 15) / +.2712814218 0729856000 0000000000 000 d-29   /
      data bk0 cs( 16) / +.3082593887 9146666666 6666666666 666 d-32   /
c
      data ntk0, xsml, xmax / 0, 2*0.d0 /
c
      if (ntk0.ne.0) go to 10
      ntk0 = initds (bk0cs, 16, 0.1*sngl(d1mach(3)))
      xsml = dsqrt (4.0d0*d1mach(3))
      xmax = -dlog(d1mach(1))
      xmax = xmax - 0.5d0*xmax*dlog(xmax)/(xmax+0.5d0)
c
 10   if (x.le.0.d0) call seteru (29hdbesk0  x is zero or negative,
     1  29, 2, 2)
      if (x.gt.2.0d0) go to 20
c
      y = 0.d0
      if (x.gt.xsml) y = x*x
      dbesk0 = -dlog(0.5d0*x)*dbesi0(x) - 0.25d0 + dcsevl (.5d0*y-1.d0,
     1  bk0cs, ntk0)
      return
c
 20   dbesk0 = 0.d0
      if (x.gt.xmax) call seteru (30hdbesk0  x so big k0 underflows,
     1  30, 1, 0)
      if (x.gt.xmax) return
c
      dbesk0 = dexp(-x) * dbsk0e(x)
c
      return
      end
      double precision function dbesk1 (x)
c july 1980 edition.  w. fullerton, c3, los alamos scientific lab.
      double precision x, bk1cs(16), xmax, xmin, xsml, y,
     1  d1mach, dcsevl, dbesi1, dbsk1e, dexp, dlog, dsqrt
c     external d1mach, dbesi1, dbsk1e, dcsevl, dexp, dlog, dsqrt, initds
c
c series for bk1        on the interval  0.          to  4.00000e+00
c                                        with weighted error   9.16e-32
c                                         log weighted error  31.04
c                               significant figures required  30.61
c                                    decimal places required  31.64
c
      data bk1 cs(  1) / +.2530022733 8947770532 5311208685 33 d-1     /
      data bk1 cs(  2) / -.3531559607 7654487566 7238316918 01 d+0     /
      data bk1 cs(  3) / -.1226111808 2265714823 4790679300 42 d+0     /
      data bk1 cs(  4) / -.6975723859 6398643501 8129202960 83 d-2     /
      data bk1 cs(  5) / -.1730288957 5130520630 1765073689 79 d-3     /
      data bk1 cs(  6) / -.2433406141 5659682349 6007350301 64 d-5     /
      data bk1 cs(  7) / -.2213387630 7347258558 3152525451 26 d-7     /
      data bk1 cs(  8) / -.1411488392 6335277610 9583302126 08 d-9     /
      data bk1 cs(  9) / -.6666901694 1993290060 8537512643 73 d-12    /
      data bk1 cs( 10) / -.2427449850 5193659339 2631968648 53 d-14    /
      data bk1 cs( 11) / -.7023863479 3862875971 7837971200 00 d-17    /
      data bk1 cs( 12) / -.1654327515 5100994675 4910293333 33 d-19    /
      data bk1 cs( 13) / -.3233834745 9944491991 8933333333 33 d-22    /
      data bk1 cs( 14) / -.5331275052 9265274999 4666666666 66 d-25    /
      data bk1 cs( 15) / -.7513040716 2157226666 6666666666 66 d-28    /
      data bk1 cs( 16) / -.9155085717 6541866666 6666666666 66 d-31    /
c
      data ntk1, xmin, xsml, xmax / 0, 3*0.d0 /
c
      if (ntk1.ne.0) go to 10
      ntk1 = initds (bk1cs, 16, 0.1*sngl(d1mach(3)))
      xmin = dexp (dmax1(dlog(d1mach(1)), -dlog(d1mach(2))) + 0.01d0)
      xsml = dsqrt (4.0d0*d1mach(3))
      xmax = -dlog(d1mach(1))
      xmax = xmax - 0.5d0*xmax*dlog(xmax)/(xmax+0.5d0) - 0.01d0
c
 10   if (x.le.0.d0) call seteru (29hdbesk1  x is zero or negative,
     1  29, 2, 2)
      if (x.gt.2.0d0) go to 20
c
      if (x.lt.xmin) call seteru (31hdbesk1  x so small k1 overflows,
     1  31, 3, 2)
      y = 0.d0
      if (x.gt.xsml) y = x*x
      dbesk1 = dlog(0.5d0*x)*dbesi1(x) + (0.75d0 + dcsevl (.5d0*y-1.d0,
     1  bk1cs, ntk1))/x
      return
c
 20   dbesk1 = 0.d0
      if (x.gt.xmax) call seteru (30hdbesk1  x so big k1 underflows,
     1  30, 1, 0)
      if (x.gt.xmax) return
c
      dbesk1 = dexp(-x) * dbsk1e(x)
c
      return
      end
      double precision function dbesy0 (x)
c august 1980 edition.  w. fullerton, c3, los alamos scientific lab.
      double precision x, by0cs(19), ampl, theta, twodpi, xsml,
     1  y, alnhaf, d1mach, dcsevl, dbesj0, dlog, dsin, dsqrt
c     external d1mach, dbesj0, dcsevl, dlog, dsin, dsqrt, initds
c
c series for by0        on the interval  0.          to  1.60000e+01
c                                        with weighted error   8.14e-32
c                                         log weighted error  31.09
c                               significant figures required  30.31
c                                    decimal places required  31.73
c
      data by0 cs(  1) / -.1127783939 2865573217 9398054602 8 d-1      /
      data by0 cs(  2) / -.1283452375 6042034604 8088453183 8 d+0      /
      data by0 cs(  3) / -.1043788479 9794249365 8176227661 8 d+0      /
      data by0 cs(  4) / +.2366274918 3969695409 2415926461 3 d-1      /
      data by0 cs(  5) / -.2090391647 7004862391 9622395034 2 d-2      /
      data by0 cs(  6) / +.1039754539 3905725209 9924657638 1 d-3      /
      data by0 cs(  7) / -.3369747162 4239720967 1877534503 7 d-5      /
      data by0 cs(  8) / +.7729384267 6706671585 2136721637 1 d-7      /
      data by0 cs(  9) / -.1324976772 6642595914 4347606896 4 d-8      /
      data by0 cs( 10) / +.1764823261 5404527921 0038936315 8 d-10     /
      data by0 cs( 11) / -.1881055071 5801962006 0282301206 9 d-12     /
      data by0 cs( 12) / +.1641865485 3661495027 9223718574 9 d-14     /
      data by0 cs( 13) / -.1195659438 6046060857 4599100672 0 d-16     /
      data by0 cs( 14) / +.7377296297 4401858424 9411242666 6 d-19     /
      data by0 cs( 15) / -.3906843476 7104373307 4090666666 6 d-21     /
      data by0 cs( 16) / +.1795503664 4361579498 2912000000 0 d-23     /
      data by0 cs( 17) / -.7229627125 4480104789 3333333333 3 d-26     /
      data by0 cs( 18) / +.2571727931 6351685973 3333333333 3 d-28     /
      data by0 cs( 19) / -.8141268814 1636949333 3333333333 3 d-31     /
c
      data twodpi / 0.6366197723 6758134307 5535053490 057 d0 /
      data alnhaf /-0.6931471805 5994530941 7232121458 18d0 /
      data nty0, xsml / 0, 0.d0 /
c
      if (nty0.ne.0) go to 10
      nty0 = initds (by0cs, 19, 0.1*sngl(d1mach(3)))
      xsml = dsqrt (4.0d0*d1mach(3))
c
 10   if (x.le.0.d0) call seteru (29hdbesy0  x is zero or negative, 29,
     1  1, 2)
      if (x.gt.4.0d0) go to 20
c
      y = 0.d0
      if (x.gt.xsml) y = x*x
      dbesy0 = twodpi*(alnhaf+dlog(x))*dbesj0(x) + .375d0
     1  + dcsevl (.125d0*y-1.d0, by0cs, nty0)
      return
c
 20   call d9b0mp (x, ampl, theta)
      dbesy0 = ampl * dsin(theta)
      return
c
      end
      double precision function dbesy1 (x)
c july 1977 edition.  w. fullerton, c3, los alamos scientific lab.
      double precision x, by1cs(20), ampl, theta, twodpi, xmin, xsml,
     1  y, d1mach, dcsevl, dbesj1, dexp, dlog, dsin, dsqrt
c     external d1mach, dbesj1, dcsevl, dexp, dlog, dsin, dsqrt, initds
c
c series for by1        on the interval  0.          to  1.60000e+01
c                                        with weighted error   8.65e-33
c                                         log weighted error  32.06
c                               significant figures required  32.17
c                                    decimal places required  32.71
c
      data by1 cs(  1) / +.3208047100 6119086293 2352018628 015 d-1    /
      data by1 cs(  2) / +.1262707897 4335004495 3431725999 727 d+1    /
      data by1 cs(  3) / +.6499961899 9231750009 7490637314 144 d-2    /
      data by1 cs(  4) / -.8936164528 8605041165 3144160009 712 d-1    /
      data by1 cs(  5) / +.1325088122 1757095451 2375510370 043 d-1    /
      data by1 cs(  6) / -.8979059119 6483523775 3039508298 105 d-3    /
      data by1 cs(  7) / +.3647361487 9583067824 2287368165 349 d-4    /
      data by1 cs(  8) / -.1001374381 6660005554 9075523845 295 d-5    /
      data by1 cs(  9) / +.1994539657 3901739703 1159372421 243 d-7    /
      data by1 cs( 10) / -.3023065601 8033816728 4799332520 743 d-9    /
      data by1 cs( 11) / +.3609878156 9478119611 6252914242 474 d-11   /
      data by1 cs( 12) / -.3487488297 2875824241 4552947409 066 d-13   /
      data by1 cs( 13) / +.2783878971 5591766581 3507698517 333 d-15   /
      data by1 cs( 14) / -.1867870968 6194876876 6825352533 333 d-17   /
      data by1 cs( 15) / +.1068531533 9116825975 7070336000 000 d-19   /
      data by1 cs( 16) / -.5274721956 6844822894 3872000000 000 d-22   /
      data by1 cs( 17) / +.2270199403 1556641437 0133333333 333 d-24   /
      data by1 cs( 18) / -.8595390353 9452310869 3333333333 333 d-27   /
      data by1 cs( 19) / +.2885404379 8337945600 0000000000 000 d-29   /
      data by1 cs( 20) / -.8647541138 9371733333 3333333333 333 d-32   /
c
      data twodpi / 0.6366197723 6758134307 5535053490 057 d0 /
      data nty1, xmin, xsml / 0, 2*0.d0 /
c
      if (nty1.ne.0) go to 10
      nty1 = initds (by1cs, 20, 0.1*sngl(d1mach(3)))
c
      xmin = 1.571d0 * dexp (dmax1(dlog(d1mach(1)), -dlog(d1mach(2))) +
     1  0.01d0)
      xsml = dsqrt (4.0d0*d1mach(3))
c
 10   if (x.le.0.d0) call seteru (29hdbesy1  x is zero or negative, 29,
     1  1, 2)
      if (x.gt.4.0d0) go to 20
c
      if (x.lt.xmin) call seteru (31hdbesy1  x so small y1 overflows,
     1  31, 3, 2)
      y = 0.d0
      if (x.gt.xsml) y = x*x
      dbesy1 = twodpi * dlog(0.5d0*x)*dbesj1(x) + (0.5d0 +
     1  dcsevl (.125d0*y-1.d0, by1cs, nty1))/x
      return
c
 20   call d9b1mp (x, ampl, theta)
      dbesy1 = ampl * dsin(theta)
      return
c
      end
      function initds (dos, nos, eta)
c june 1977 edition.   w. fullerton, c3, los alamos scientific lab.
c
c initialize the double precision orthogonal series dos so that initds
c is the number of terms needed to insure the error is no larger than
c eta.  ordinarily eta will be chosen to be one-tenth machine precision.
c
c             input arguments --
c dos    dble prec array of nos coefficients in an orthogonal series.
c nos    number of coefficients in dos.
c eta    requested accuracy of series.
c
      double precision dos(nos)
c
      if (nos.lt.1) call seteru (
     1  35hinitds  number of coefficients lt 1, 35, 2, 2)
c
      err = 0.
      do 10 ii=1,nos
        i = nos + 1 - ii
        err = err + abs(sngl(dos(i)))
        if (err.gt.eta) go to 20
 10   continue
c
 20   if (i.eq.nos) call seteru (28hinitds  eta may be too small, 28,
     1  1, 2)
      initds = i
c
      return
      end
      double precision function dcsevl (x, a, n)
c
c evaluate the n-term chebyshev series a at x.  adapted from
c r. broucke, algorithm 446, c.a.c.m., 16, 254 (1973).
c
c             input arguments --
c x      dble prec value at which the series is to be evaluated.
c a      dble prec array of n terms of a chebyshev series.  in eval-
c        uating a, only half the first coef is summed.
c n      number of terms in array a.
c
      double precision a(n), x, twox, b0, b1, b2
c
      if (n.lt.1) call seteru (28hdcsevl  number of terms le 0, 28, 2,2)
      if (n.gt.1000) call seteru (31hdcsevl  number of terms gt 1000,
     1  31, 3, 2)
      if (x.lt.(-1.1d0) .or. x.gt.1.1d0) call seteru (
     1  25hdcsevl  x outside (-1,+1), 25, 1, 1)
c
      twox = 2.0d0*x
      b1 = 0.d0
      b0 = 0.d0
      do 10 i=1,n
        b2 = b1
        b1 = b0
        ni = n - i + 1
        b0 = twox*b1 - b2 + a(ni)
 10   continue
c
      dcsevl = 0.5d0 * (b0-b2)
c
      return
      end
      subroutine seteru (messg, nmessg, nerr, iopt)
      common /cseter/ iunflo
      integer messg(1)
      data iunflo / 0 /
c
      if (iopt.ne.0) call seterr (messg, nmessg, nerr, iopt)
      if (iopt.ne.0) return
c
      if (iunflo.le.0) return
      call seterr (messg, nmessg, nerr, 1)
c
      return
      end
      double precision function dbsi0e (x)
c july 1977 edition.  w. fullerton, c3, los alamos scientific lab.
      double precision x, bi0cs(18), ai0cs(46), ai02cs(69),
     1  xsml, y, d1mach, dcsevl, dexp, dsqrt
c     external d1mach, dcsevl, dexp, dsqrt, initds
c
c series for bi0        on the interval  0.          to  9.00000e+00
c                                        with weighted error   9.51e-34
c                                         log weighted error  33.02
c                               significant figures required  33.31
c                                    decimal places required  33.65
c
      data bi0 cs(  1) / -.7660547252 8391449510 8189497624 3285 d-1   /
      data bi0 cs(  2) / +.1927337953 9938082699 5240875088 1196 d+1   /
      data bi0 cs(  3) / +.2282644586 9203013389 3702929233 0415 d+0   /
      data bi0 cs(  4) / +.1304891466 7072904280 7933421069 1888 d-1   /
      data bi0 cs(  5) / +.4344270900 8164874513 7868268102 6107 d-3   /
      data bi0 cs(  6) / +.9422657686 0019346639 2317174411 8766 d-5   /
      data bi0 cs(  7) / +.1434006289 5106910799 6209187817 9957 d-6   /
      data bi0 cs(  8) / +.1613849069 6617490699 1541971999 4611 d-8   /
      data bi0 cs(  9) / +.1396650044 5356696994 9509270814 2522 d-10  /
      data bi0 cs( 10) / +.9579451725 5054453446 2752317189 3333 d-13  /
      data bi0 cs( 11) / +.5333981859 8625021310 1510774400 0000 d-15  /
      data bi0 cs( 12) / +.2458716088 4374707746 9678591999 9999 d-17  /
      data bi0 cs( 13) / +.9535680890 2487700269 4434133333 3333 d-20  /
      data bi0 cs( 14) / +.3154382039 7214273367 8933333333 3333 d-22  /
      data bi0 cs( 15) / +.9004564101 0946374314 6666666666 6666 d-25  /
      data bi0 cs( 16) / +.2240647369 1236700160 0000000000 0000 d-27  /
      data bi0 cs( 17) / +.4903034603 2428373333 3333333333 3333 d-30  /
      data bi0 cs( 18) / +.9508172606 1226666666 6666666666 6666 d-33  /
c
c series for ai0        on the interval  1.25000e-01 to  3.33333e-01
c                                        with weighted error   2.74e-32
c                                         log weighted error  31.56
c                               significant figures required  30.15
c                                    decimal places required  32.39
c
      data ai0 cs(  1) / +.7575994494 0237959427 2987203743 8 d-1      /
      data ai0 cs(  2) / +.7591380810 8233455072 9297873320 4 d-2      /
      data ai0 cs(  3) / +.4153131338 9237505018 6319749138 2 d-3      /
      data ai0 cs(  4) / +.1070076463 4390730735 8242970217 0 d-4      /
      data ai0 cs(  5) / -.7901179979 2128946607 5031948573 0 d-5      /
      data ai0 cs(  6) / -.7826143501 4387522697 8898980690 9 d-6      /
      data ai0 cs(  7) / +.2783849942 9488708063 8118538985 7 d-6      /
      data ai0 cs(  8) / +.8252472600 6120271919 6682913319 8 d-8      /
      data ai0 cs(  9) / -.1204463945 5201991790 5496089110 3 d-7      /
      data ai0 cs( 10) / +.1559648598 5060764436 1228752792 8 d-8      /
      data ai0 cs( 11) / +.2292556367 1033165434 7725480285 7 d-9      /
      data ai0 cs( 12) / -.1191622884 2790646036 7777423447 8 d-9      /
      data ai0 cs( 13) / +.1757854916 0324098302 1833124774 3 d-10     /
      data ai0 cs( 14) / +.1128224463 2189005171 4441135682 4 d-11     /
      data ai0 cs( 15) / -.1146848625 9272988777 2963387698 2 d-11     /
      data ai0 cs( 16) / +.2715592054 8036628726 4365192160 6 d-12     /
      data ai0 cs( 17) / -.2415874666 5626878384 4247572028 1 d-13     /
      data ai0 cs( 18) / -.6084469888 2551250646 0609963922 4 d-14     /
      data ai0 cs( 19) / +.3145705077 1754772937 0836026730 3 d-14     /
      data ai0 cs( 20) / -.7172212924 8711877179 6217505917 6 d-15     /
      data ai0 cs( 21) / +.7874493403 4541033960 8390960332 7 d-16     /
      data ai0 cs( 22) / +.1004802753 0094624023 4524457183 9 d-16     /
      data ai0 cs( 23) / -.7566895365 3505348534 2843588881 0 d-17     /
      data ai0 cs( 24) / +.2150380106 8761198878 1205128784 5 d-17     /
      data ai0 cs( 25) / -.3754858341 8308744291 5158445260 8 d-18     /
      data ai0 cs( 26) / +.2354065842 2269925769 0075710532 2 d-19     /
      data ai0 cs( 27) / +.1114667612 0479285302 2637335511 0 d-19     /
      data ai0 cs( 28) / -.5398891884 3969903786 9677932270 9 d-20     /
      data ai0 cs( 29) / +.1439598792 2407526770 4285840452 2 d-20     /
      data ai0 cs( 30) / -.2591916360 1110934064 6081840196 2 d-21     /
      data ai0 cs( 31) / +.2238133183 9985839074 3409229824 0 d-22     /
      data ai0 cs( 32) / +.5250672575 3647711727 7221683199 9 d-23     /
      data ai0 cs( 33) / -.3249904138 5332307841 7343228586 6 d-23     /
      data ai0 cs( 34) / +.9924214103 2050379278 5728471040 0 d-24     /
      data ai0 cs( 35) / -.2164992254 2446695231 4655429973 3 d-24     /
      data ai0 cs( 36) / +.3233609471 9435940839 7333299199 9 d-25     /
      data ai0 cs( 37) / -.1184620207 3967424898 2473386666 6 d-26     /
      data ai0 cs( 38) / -.1281671853 9504986505 4833868799 9 d-26     /
      data ai0 cs( 39) / +.5827015182 2793905116 0556885333 3 d-27     /
      data ai0 cs( 40) / -.1668222326 0261097193 6450150399 9 d-27     /
      data ai0 cs( 41) / +.3625309510 5415699757 0068480000 0 d-28     /
      data ai0 cs( 42) / -.5733627999 0557135899 4595839999 9 d-29     /
      data ai0 cs( 43) / +.3736796722 0630982296 4258133333 3 d-30     /
      data ai0 cs( 44) / +.1602073983 1568519633 6551253333 3 d-30     /
      data ai0 cs( 45) / -.8700424864 0572298845 2249599999 9 d-31     /
      data ai0 cs( 46) / +.2741320937 9374811456 0341333333 3 d-31     /
c
c series for ai02       on the interval  0.          to  1.25000e-01
c                                        with weighted error   1.97e-32
c                                         log weighted error  31.71
c                               significant figures required  30.15
c                                    decimal places required  32.63
c
      data ai02cs(  1) / +.5449041101 4108831607 8960962268 0 d-1      /
      data ai02cs(  2) / +.3369116478 2556940898 9785662979 9 d-2      /
      data ai02cs(  3) / +.6889758346 9168239842 6263914301 1 d-4      /
      data ai02cs(  4) / +.2891370520 8347564829 6692402323 2 d-5      /
      data ai02cs(  5) / +.2048918589 4690637418 2760534093 1 d-6      /
      data ai02cs(  6) / +.2266668990 4981780645 9327743136 1 d-7      /
      data ai02cs(  7) / +.3396232025 7083863451 5084396952 3 d-8      /
      data ai02cs(  8) / +.4940602388 2249695891 0482449783 5 d-9      /
      data ai02cs(  9) / +.1188914710 7846438342 4084525196 3 d-10     /
      data ai02cs( 10) / -.3149916527 9632413645 3864862961 9 d-10     /
      data ai02cs( 11) / -.1321581184 0447713118 7540739926 7 d-10     /
      data ai02cs( 12) / -.1794178531 5068061177 7943574026 9 d-11     /
      data ai02cs( 13) / +.7180124451 3836662336 7106429346 9 d-12     /
      data ai02cs( 14) / +.3852778382 7421427011 4089801777 6 d-12     /
      data ai02cs( 15) / +.1540086217 5214098269 1325823339 7 d-13     /
      data ai02cs( 16) / -.4150569347 2872220866 2689972015 6 d-13     /
      data ai02cs( 17) / -.9554846698 8283076487 0214494312 5 d-14     /
      data ai02cs( 18) / +.3811680669 3526224207 4605535511 8 d-14     /
      data ai02cs( 19) / +.1772560133 0565263836 0493266675 8 d-14     /
      data ai02cs( 20) / -.3425485619 6772191346 1924790328 2 d-15     /
      data ai02cs( 21) / -.2827623980 5165834849 4205593759 4 d-15     /
      data ai02cs( 22) / +.3461222867 6974610930 9706250813 4 d-16     /
      data ai02cs( 23) / +.4465621420 2967599990 1042054284 3 d-16     /
      data ai02cs( 24) / -.4830504485 9441820712 5525403795 4 d-17     /
      data ai02cs( 25) / -.7233180487 8747539545 6227240924 5 d-17     /
      data ai02cs( 26) / +.9921475412 1736985988 8046093981 0 d-18     /
      data ai02cs( 27) / +.1193650890 8459820855 0439949924 2 d-17     /
      data ai02cs( 28) / -.2488709837 1508072357 2054491660 2 d-18     /
      data ai02cs( 29) / -.1938426454 1609059289 8469781132 6 d-18     /
      data ai02cs( 30) / +.6444656697 3734438687 8301949394 9 d-19     /
      data ai02cs( 31) / +.2886051596 2892243264 8171383073 4 d-19     /
      data ai02cs( 32) / -.1601954907 1749718070 6167156200 7 d-19     /
      data ai02cs( 33) / -.3270815010 5923147208 9193567485 9 d-20     /
      data ai02cs( 34) / +.3686932283 8264091811 4600723939 3 d-20     /
      data ai02cs( 35) / +.1268297648 0309501530 1359529710 9 d-22     /
      data ai02cs( 36) / -.7549825019 3772739076 9636664410 1 d-21     /
      data ai02cs( 37) / +.1502133571 3778353496 3712789053 4 d-21     /
      data ai02cs( 38) / +.1265195883 5096485349 3208799248 3 d-21     /
      data ai02cs( 39) / -.6100998370 0836807086 2940891600 2 d-22     /
      data ai02cs( 40) / -.1268809629 2601282643 6872095924 2 d-22     /
      data ai02cs( 41) / +.1661016099 8907414578 4038487490 5 d-22     /
      data ai02cs( 42) / -.1585194335 7658855793 7970504881 4 d-23     /
      data ai02cs( 43) / -.3302645405 9682178009 5381766755 6 d-23     /
      data ai02cs( 44) / +.1313580902 8392397817 4039623117 4 d-23     /
      data ai02cs( 45) / +.3689040246 6711567933 1425637280 4 d-24     /
      data ai02cs( 46) / -.4210141910 4616891492 1978247249 9 d-24     /
      data ai02cs( 47) / +.4791954591 0828657806 3171401373 0 d-25     /
      data ai02cs( 48) / +.8459470390 2218217952 9971707412 4 d-25     /
      data ai02cs( 49) / -.4039800940 8728324931 4607937181 0 d-25     /
      data ai02cs( 50) / -.6434714653 6504313473 0100850469 5 d-26     /
      data ai02cs( 51) / +.1225743398 8756659903 4464736990 5 d-25     /
      data ai02cs( 52) / -.2934391316 0257089231 9879821175 4 d-26     /
      data ai02cs( 53) / -.1961311309 1949829262 0371205728 9 d-26     /
      data ai02cs( 54) / +.1503520374 8221934241 6229900309 8 d-26     /
      data ai02cs( 55) / -.9588720515 7448265520 3386388206 9 d-28     /
      data ai02cs( 56) / -.3483339380 8170454863 9441108511 4 d-27     /
      data ai02cs( 57) / +.1690903610 2630436730 6244960725 6 d-27     /
      data ai02cs( 58) / +.1982866538 7356030438 9400115718 8 d-28     /
      data ai02cs( 59) / -.5317498081 4918162145 7583002528 4 d-28     /
      data ai02cs( 60) / +.1803306629 8883929462 3501450390 1 d-28     /
      data ai02cs( 61) / +.6213093341 4548931758 8405311242 2 d-29     /
      data ai02cs( 62) / -.7692189292 7721618632 0072806673 0 d-29     /
      data ai02cs( 63) / +.1858252826 1117025426 2556016596 3 d-29     /
      data ai02cs( 64) / +.1237585142 2813957248 9927154554 1 d-29     /
      data ai02cs( 65) / -.1102259120 4092238032 1779478779 2 d-29     /
      data ai02cs( 66) / +.1886287118 0397044900 7787447943 1 d-30     /
      data ai02cs( 67) / +.2160196872 2436589131 4903141406 0 d-30     /
      data ai02cs( 68) / -.1605454124 9197432005 8446594965 5 d-30     /
      data ai02cs( 69) / +.1965352984 5942906039 3884807331 8 d-31     /
c
      data nti0, ntai0, ntai02, xsml / 3*0, 0.d0 /
c
      if (nti0.ne.0) go to 10
      eta = 0.1*sngl(d1mach(3))
      nti0 = initds (bi0cs, 18, eta)
      ntai0 = initds (ai0cs, 46, eta)
      ntai02 = initds (ai02cs, 69, eta)
      xsml = dsqrt (8.0d0*d1mach(3))
c
 10   y = dabs(x)
      if (y.gt.3.0d0) go to 20
c
      dbsi0e = 1.0d0
      if (y.gt.xsml) dbsi0e = dexp(-y) * (2.75d0 +
     1  dcsevl (y*y/4.5d0-1.d0, bi0cs, nti0) )
      return
c
 20   if (y.le.8.d0) dbsi0e = (0.375d0 + dcsevl ((48.d0/y-11.d0)/5.d0,
     1  ai0cs, ntai0))/dsqrt(y)
      if (y.gt.8.d0) dbsi0e = (0.375d0 + dcsevl (16.d0/y-1.d0, ai02cs,
     1  ntai02))/dsqrt(y)
c
      return
      end
      double precision function dbsi1e (x)
c oct 1983 edition.  w. fullerton, c3, los alamos scientific lab.
      double precision x, bi1cs(17), ai1cs(46), ai12cs(69), xmin,
     1  xsml, y, d1mach, dcsevl, dexp, dsqrt
c     external d1mach, dcsevl, dexp, dsqrt, initds
c
c series for bi1        on the interval  0.          to  9.00000e+00
c                                        with weighted error   1.44e-32
c                                         log weighted error  31.84
c                               significant figures required  31.45
c                                    decimal places required  32.46
c
      data bi1 cs(  1) / -.1971713261 0998597316 1385032181 49 d-2     /
      data bi1 cs(  2) / +.4073488766 7546480608 1553936520 14 d+0     /
      data bi1 cs(  3) / +.3483899429 9959455866 2450377837 87 d-1     /
      data bi1 cs(  4) / +.1545394556 3001236038 5984010584 89 d-2     /
      data bi1 cs(  5) / +.4188852109 8377784129 4588320041 20 d-4     /
      data bi1 cs(  6) / +.7649026764 8362114741 9597039660 69 d-6     /
      data bi1 cs(  7) / +.1004249392 4741178689 1798080372 38 d-7     /
      data bi1 cs(  8) / +.9932207791 9238106481 3712980548 63 d-10    /
      data bi1 cs(  9) / +.7663801791 8447637275 2001716813 49 d-12    /
      data bi1 cs( 10) / +.4741418923 8167394980 3880919481 60 d-14    /
      data bi1 cs( 11) / +.2404114404 0745181799 8631720320 00 d-16    /
      data bi1 cs( 12) / +.1017150500 7093713649 1211007999 99 d-18    /
      data bi1 cs( 13) / +.3645093565 7866949458 4917333333 33 d-21    /
      data bi1 cs( 14) / +.1120574950 2562039344 8106666666 66 d-23    /
      data bi1 cs( 15) / +.2987544193 4468088832 0000000000 00 d-26    /
      data bi1 cs( 16) / +.6973231093 9194709333 3333333333 33 d-29    /
      data bi1 cs( 17) / +.1436794822 0620800000 0000000000 00 d-31    /
c
c series for ai1        on the interval  1.25000e-01 to  3.33333e-01
c                                        with weighted error   2.81e-32
c                                         log weighted error  31.55
c                               significant figures required  29.93
c                                    decimal places required  32.38
c
      data ai1 cs(  1) / -.2846744181 8814786741 0037246830 7 d-1      /
      data ai1 cs(  2) / -.1922953231 4432206510 4444877497 9 d-1      /
      data ai1 cs(  3) / -.6115185857 9437889822 5624991778 5 d-3      /
      data ai1 cs(  4) / -.2069971253 3502277088 8282377797 9 d-4      /
      data ai1 cs(  5) / +.8585619145 8107255655 3694467313 8 d-5      /
      data ai1 cs(  6) / +.1049498246 7115908625 1745399786 0 d-5      /
      data ai1 cs(  7) / -.2918338918 4479022020 9343232669 7 d-6      /
      data ai1 cs(  8) / -.1559378146 6317390001 6068096907 7 d-7      /
      data ai1 cs(  9) / +.1318012367 1449447055 2530287390 9 d-7      /
      data ai1 cs( 10) / -.1448423418 1830783176 3913446781 5 d-8      /
      data ai1 cs( 11) / -.2908512243 9931420948 2504099301 0 d-9      /
      data ai1 cs( 12) / +.1266388917 8753823873 1115969040 3 d-9      /
      data ai1 cs( 13) / -.1664947772 9192206706 2417839858 0 d-10     /
      data ai1 cs( 14) / -.1666653644 6094329760 9593715499 9 d-11     /
      data ai1 cs( 15) / +.1242602414 2907682652 3216847201 7 d-11     /
      data ai1 cs( 16) / -.2731549379 6724323972 5146142863 3 d-12     /
      data ai1 cs( 17) / +.2023947881 6458037807 0026268898 1 d-13     /
      data ai1 cs( 18) / +.7307950018 1168836361 9869812612 3 d-14     /
      data ai1 cs( 19) / -.3332905634 4046749438 1377861713 3 d-14     /
      data ai1 cs( 20) / +.7175346558 5129537435 4225466567 0 d-15     /
      data ai1 cs( 21) / -.6982530324 7962563558 5062922365 6 d-16     /
      data ai1 cs( 22) / -.1299944201 5627607600 6044608058 7 d-16     /
      data ai1 cs( 23) / +.8120942864 2427988920 5467834286 0 d-17     /
      data ai1 cs( 24) / -.2194016207 4107368981 5626664378 3 d-17     /
      data ai1 cs( 25) / +.3630516170 0296548482 7986093233 4 d-18     /
      data ai1 cs( 26) / -.1695139772 4391041663 0686679039 9 d-19     /
      data ai1 cs( 27) / -.1288184829 8979078071 1688253822 2 d-19     /
      data ai1 cs( 28) / +.5694428604 9670527801 0999107310 9 d-20     /
      data ai1 cs( 29) / -.1459597009 0904800565 4550990028 7 d-20     /
      data ai1 cs( 30) / +.2514546010 6757173140 8469133448 5 d-21     /
      data ai1 cs( 31) / -.1844758883 1391248181 6040002901 3 d-22     /
      data ai1 cs( 32) / -.6339760596 2279486419 2860979199 9 d-23     /
      data ai1 cs( 33) / +.3461441102 0310111111 0814662656 0 d-23     /
      data ai1 cs( 34) / -.1017062335 3713935475 9654102357 3 d-23     /
      data ai1 cs( 35) / +.2149877147 0904314459 6250077866 6 d-24     /
      data ai1 cs( 36) / -.3045252425 2386764017 4620617386 6 d-25     /
      data ai1 cs( 37) / +.5238082144 7212859821 7763498666 6 d-27     /
      data ai1 cs( 38) / +.1443583107 0893824464 1678950399 9 d-26     /
      data ai1 cs( 39) / -.6121302074 8900427332 0067071999 9 d-27     /
      data ai1 cs( 40) / +.1700011117 4678184183 4918980266 6 d-27     /
      data ai1 cs( 41) / -.3596589107 9842441585 3521578666 6 d-28     /
      data ai1 cs( 42) / +.5448178578 9484185766 5051306666 6 d-29     /
      data ai1 cs( 43) / -.2731831789 6890849891 6256426666 6 d-30     /
      data ai1 cs( 44) / -.1858905021 7086007157 7190399999 9 d-30     /
      data ai1 cs( 45) / +.9212682974 5139334411 2776533333 3 d-31     /
      data ai1 cs( 46) / -.2813835155 6535611063 7083306666 6 d-31     /
c
c series for ai12       on the interval  0.          to  1.25000e-01
c                                        with weighted error   1.83e-32
c                                         log weighted error  31.74
c                               significant figures required  29.97
c                                    decimal places required  32.66
c
      data ai12cs(  1) / +.2857623501 8280120474 4984594846 9 d-1      /
      data ai12cs(  2) / -.9761097491 3614684077 6516445730 2 d-2      /
      data ai12cs(  3) / -.1105889387 6262371629 1256921277 5 d-3      /
      data ai12cs(  4) / -.3882564808 8776903934 5654477627 4 d-5      /
      data ai12cs(  5) / -.2512236237 8702089252 9452002212 1 d-6      /
      data ai12cs(  6) / -.2631468846 8895195068 3705236523 2 d-7      /
      data ai12cs(  7) / -.3835380385 9642370220 4500678796 8 d-8      /
      data ai12cs(  8) / -.5589743462 1965838068 6811252222 9 d-9      /
      data ai12cs(  9) / -.1897495812 3505412344 9892503323 8 d-10     /
      data ai12cs( 10) / +.3252603583 0154882385 5508067994 9 d-10     /
      data ai12cs( 11) / +.1412580743 6613781331 6336633284 6 d-10     /
      data ai12cs( 12) / +.2035628544 1470895072 2452613684 0 d-11     /
      data ai12cs( 13) / -.7198551776 2459085120 9258989044 6 d-12     /
      data ai12cs( 14) / -.4083551111 0921973182 2849963969 1 d-12     /
      data ai12cs( 15) / -.2101541842 7726643130 1984572746 2 d-13     /
      data ai12cs( 16) / +.4272440016 7119513542 9778833699 7 d-13     /
      data ai12cs( 17) / +.1042027698 4128802764 1741449994 8 d-13     /
      data ai12cs( 18) / -.3814403072 4370078047 6707253539 6 d-14     /
      data ai12cs( 19) / -.1880354775 5107824485 1273453396 3 d-14     /
      data ai12cs( 20) / +.3308202310 9209282827 3190335240 5 d-15     /
      data ai12cs( 21) / +.2962628997 6459501390 6854654205 2 d-15     /
      data ai12cs( 22) / -.3209525921 9934239587 7837353288 7 d-16     /
      data ai12cs( 23) / -.4650305368 4893583255 7128281897 9 d-16     /
      data ai12cs( 24) / +.4414348323 0717079499 4611375964 1 d-17     /
      data ai12cs( 25) / +.7517296310 8421048054 2545808029 5 d-17     /
      data ai12cs( 26) / -.9314178867 3268833756 8484784515 7 d-18     /
      data ai12cs( 27) / -.1242193275 1948909561 1678448869 7 d-17     /
      data ai12cs( 28) / +.2414276719 4548484690 0515390217 6 d-18     /
      data ai12cs( 29) / +.2026944384 0532851789 7192286069 2 d-18     /
      data ai12cs( 30) / -.6394267188 2690977870 4391988681 1 d-19     /
      data ai12cs( 31) / -.3049812452 3730958960 8488450357 1 d-19     /
      data ai12cs( 32) / +.1612841851 6514802251 3462230769 1 d-19     /
      data ai12cs( 33) / +.3560913964 3099250545 1027090462 0 d-20     /
      data ai12cs( 34) / -.3752017947 9364390796 6682800324 6 d-20     /
      data ai12cs( 35) / -.5787037427 0747993459 5198231074 1 d-22     /
      data ai12cs( 36) / +.7759997511 6481619619 8236963209 2 d-21     /
      data ai12cs( 37) / -.1452790897 2022333940 6445987408 5 d-21     /
      data ai12cs( 38) / -.1318225286 7390367021 2192275337 4 d-21     /
      data ai12cs( 39) / +.6116654862 9030707018 7999133171 7 d-22     /
      data ai12cs( 40) / +.1376279762 4271264277 3024338363 4 d-22     /
      data ai12cs( 41) / -.1690837689 9593478849 1983938230 6 d-22     /
      data ai12cs( 42) / +.1430596088 5954331539 8720108538 5 d-23     /
      data ai12cs( 43) / +.3409557828 0905940204 0536772990 2 d-23     /
      data ai12cs( 44) / -.1309457666 2707602278 4573872642 4 d-23     /
      data ai12cs( 45) / -.3940706411 2402574360 9352141755 7 d-24     /
      data ai12cs( 46) / +.4277137426 9808765808 0616679735 2 d-24     /
      data ai12cs( 47) / -.4424634830 9826068819 0028312302 9 d-25     /
      data ai12cs( 48) / -.8734113196 2307149721 1530978874 7 d-25     /
      data ai12cs( 49) / +.4045401335 6835333921 4340414242 8 d-25     /
      data ai12cs( 50) / +.7067100658 0946894656 5160771780 6 d-26     /
      data ai12cs( 51) / -.1249463344 5651052230 0286451860 5 d-25     /
      data ai12cs( 52) / +.2867392244 4034370329 7948339142 6 d-26     /
      data ai12cs( 53) / +.2044292892 5042926702 8177957421 0 d-26     /
      data ai12cs( 54) / -.1518636633 8204625683 7134680291 1 d-26     /
      data ai12cs( 55) / +.8110181098 1875758861 3227910703 7 d-28     /
      data ai12cs( 56) / +.3580379354 7735860911 2717370327 0 d-27     /
      data ai12cs( 57) / -.1692929018 9279025095 9305717544 8 d-27     /
      data ai12cs( 58) / -.2222902499 7024276390 6775852777 4 d-28     /
      data ai12cs( 59) / +.5424535127 1459696550 4860040112 8 d-28     /
      data ai12cs( 60) / -.1787068401 5780186887 6491299330 4 d-28     /
      data ai12cs( 61) / -.6565479068 7228149388 2392943788 0 d-29     /
      data ai12cs( 62) / +.7807013165 0611452809 2206770683 9 d-29     /
      data ai12cs( 63) / -.1816595260 6689797173 7933315222 1 d-29     /
      data ai12cs( 64) / -.1287704952 6600848203 7687559895 9 d-29     /
      data ai12cs( 65) / +.1114548172 9881645474 1370927369 4 d-29     /
      data ai12cs( 66) / -.1808343145 0393369391 5936887668 7 d-30     /
      data ai12cs( 67) / -.2231677718 2037719522 3244822893 9 d-30     /
      data ai12cs( 68) / +.1619029596 0803415106 1790980361 4 d-30     /
      data ai12cs( 69) / -.1834079908 8049414139 0130843921 0 d-31     /
c
      data nti1, ntai1, ntai12, xmin, xsml / 3*0, 2*0.d0 /
c
      if (nti1.ne.0) go to 10
      eta = 0.1*sngl(d1mach(3))
      nti1 = initds (bi1cs, 17, eta)
      ntai1 = initds (ai1cs, 46, eta)
      ntai12 = initds (ai12cs, 69, eta)
c
      xmin = 2.0d0*d1mach(1)
      xsml = dsqrt (8.0d0*d1mach(3))
c
 10   y = dabs(x)
      if (y.gt.3.0d0) go to 20
c
      dbsi1e = 0.0d0
      if (y.eq.0.0d0) return
c
      if (y.le.xmin) call seteru (
     1  38hdbsi1e  dabs(x) so small i1 underflows, 37, 1, 0)
      if (y.gt.xmin) dbsi1e = 0.5d0*x
      if (y.gt.xsml) dbsi1e = x*(0.875d0 + dcsevl (y*y/4.5d0-1.d0,
     1  bi1cs, nti1) )
      dbsi1e = dexp(-y) * dbsi1e
      return
c
 20   if (y.le.8.d0) dbsi1e = (0.375d0 + dcsevl ((48.d0/y-11.d0)/5.d0,
     1  ai1cs, ntai1))/dsqrt(y)
      if (y.gt.8.d0) dbsi1e = (0.375d0 + dcsevl (16.d0/y-1.d0, ai12cs,
     1  ntai12))/dsqrt(y)
      dbsi1e = dsign (dbsi1e, x)
c
      return
      end
      subroutine d9b0mp (x, ampl, theta)
c july 1977 edition.  w. fullerton, c3, los alamos scientific lab.
c
c evaluate the modulus and phase for the bessel j0 and y0 functions.
c
      double precision x, ampl, theta, bm0cs(37), bt02cs(39),
     1  bm02cs(40), bth0cs(44), xmax, pi4, z, d1mach, dcsevl,
     2  dsqrt
c     external d1mach, dcsevl, dsqrt, initds
c
c series for bm0        on the interval  1.56250e-02 to  6.25000e-02
c                                        with weighted error   4.40e-32
c                                         log weighted error  31.36
c                               significant figures required  30.02
c                                    decimal places required  32.14
c
      data bm0 cs(  1) / +.9211656246 8277427125 7376773018 2 d-1      /
      data bm0 cs(  2) / -.1050590997 2719051024 8071637175 5 d-2      /
      data bm0 cs(  3) / +.1470159840 7687597540 5639285095 2 d-4      /
      data bm0 cs(  4) / -.5058557606 0385542233 4792932770 2 d-6      /
      data bm0 cs(  5) / +.2787254538 6324441766 3035613788 1 d-7      /
      data bm0 cs(  6) / -.2062363611 7809148026 1884101897 3 d-8      /
      data bm0 cs(  7) / +.1870214313 1388796751 3817259626 1 d-9      /
      data bm0 cs(  8) / -.1969330971 1356362002 4173077782 5 d-10     /
      data bm0 cs(  9) / +.2325973793 9992754440 1250881805 2 d-11     /
      data bm0 cs( 10) / -.3009520344 9382502728 5122473448 2 d-12     /
      data bm0 cs( 11) / +.4194521333 8506691814 7120676864 6 d-13     /
      data bm0 cs( 12) / -.6219449312 1884458259 7326742956 4 d-14     /
      data bm0 cs( 13) / +.9718260411 3360684696 0176588526 9 d-15     /
      data bm0 cs( 14) / -.1588478585 7010752073 6663596693 7 d-15     /
      data bm0 cs( 15) / +.2700072193 6713088900 8621732445 8 d-16     /
      data bm0 cs( 16) / -.4750092365 2340089924 7750478677 3 d-17     /
      data bm0 cs( 17) / +.8615128162 6043708731 9170374656 0 d-18     /
      data bm0 cs( 18) / -.1605608686 9561448157 4560270335 9 d-18     /
      data bm0 cs( 19) / +.3066513987 3144829751 8853980159 9 d-19     /
      data bm0 cs( 20) / -.5987764223 1939564306 9650561706 6 d-20     /
      data bm0 cs( 21) / +.1192971253 7482483064 8906984106 6 d-20     /
      data bm0 cs( 22) / -.2420969142 0448054894 8468258133 3 d-21     /
      data bm0 cs( 23) / +.4996751760 5106164533 7100287999 9 d-22     /
      data bm0 cs( 24) / -.1047493639 3511585100 9504051199 9 d-22     /
      data bm0 cs( 25) / +.2227786843 7974681010 4818346666 6 d-23     /
      data bm0 cs( 26) / -.4801813239 3981628623 7054293333 3 d-24     /
      data bm0 cs( 27) / +.1047962723 4709599564 7699626666 6 d-24     /
      data bm0 cs( 28) / -.2313858165 6786153251 0126080000 0 d-25     /
      data bm0 cs( 29) / +.5164823088 4626742116 3519999999 9 d-26     /
      data bm0 cs( 30) / -.1164691191 8500653895 2540159999 9 d-26     /
      data bm0 cs( 31) / +.2651788486 0433192829 5833600000 0 d-27     /
      data bm0 cs( 32) / -.6092559503 8257284976 9130666666 6 d-28     /
      data bm0 cs( 33) / +.1411804686 1442593080 3882666666 6 d-28     /
      data bm0 cs( 34) / -.3298094961 2317372457 5061333333 3 d-29     /
      data bm0 cs( 35) / +.7763931143 0740650317 1413333333 3 d-30     /
      data bm0 cs( 36) / -.1841031343 6614584784 2133333333 3 d-30     /
      data bm0 cs( 37) / +.4395880138 5943107371 0079999999 9 d-31     /
c
c series for bth0       on the interval  0.          to  1.56250e-02
c                                        with weighted error   2.66e-32
c                                         log weighted error  31.57
c                               significant figures required  30.67
c                                    decimal places required  32.40
c
      data bth0cs(  1) / -.2490178086 2128936717 7097937899 67 d+0     /
      data bth0cs(  2) / +.4855029960 9623749241 0486155354 85 d-3     /
      data bth0cs(  3) / -.5451183734 5017204950 6562735635 05 d-5     /
      data bth0cs(  4) / +.1355867305 9405964054 3774459299 03 d-6     /
      data bth0cs(  5) / -.5569139890 2227626227 5832184149 20 d-8     /
      data bth0cs(  6) / +.3260903182 4994335304 0042057194 68 d-9     /
      data bth0cs(  7) / -.2491880786 2461341125 2379038779 93 d-10    /
      data bth0cs(  8) / +.2344937742 0882520554 3524135648 91 d-11    /
      data bth0cs(  9) / -.2609653444 4310387762 1775747661 36 d-12    /
      data bth0cs( 10) / +.3335314042 0097395105 8699550149 23 d-13    /
      data bth0cs( 11) / -.4789000044 0572684646 7507705574 09 d-14    /
      data bth0cs( 12) / +.7595617843 6192215972 6425685452 48 d-15    /
      data bth0cs( 13) / -.1313155601 6891440382 7733974876 33 d-15    /
      data bth0cs( 14) / +.2448361834 5240857495 4268207383 55 d-16    /
      data bth0cs( 15) / -.4880572981 0618777683 2567619183 31 d-17    /
      data bth0cs( 16) / +.1032728502 9786316149 2237563612 04 d-17    /
      data bth0cs( 17) / -.2305763381 5057217157 0047445270 25 d-18    /
      data bth0cs( 18) / +.5404444300 1892693993 0171084837 65 d-19    /
      data bth0cs( 19) / -.1324069519 4366572724 1550328823 85 d-19    /
      data bth0cs( 20) / +.3378079562 1371970203 4247921247 22 d-20    /
      data bth0cs( 21) / -.8945762915 7111779003 0269262922 99 d-21    /
      data bth0cs( 22) / +.2451990688 9219317090 8999086514 05 d-21    /
      data bth0cs( 23) / -.6938842287 6866318680 1399331576 57 d-22    /
      data bth0cs( 24) / +.2022827871 4890138392 9463033377 91 d-22    /
      data bth0cs( 25) / -.6062850000 2335483105 7941953717 64 d-23    /
      data bth0cs( 26) / +.1864974896 4037635381 8237883962 70 d-23    /
      data bth0cs( 27) / -.5878373238 4849894560 2450365308 67 d-24    /
      data bth0cs( 28) / +.1895859144 7999563485 5311795035 13 d-24    /
      data bth0cs( 29) / -.6248197937 2258858959 2916207285 65 d-25    /
      data bth0cs( 30) / +.2101790168 4551024686 6386335290 74 d-25    /
      data bth0cs( 31) / -.7208430093 5209253690 8139339924 46 d-26    /
      data bth0cs( 32) / +.2518136389 2474240867 1564059767 46 d-26    /
      data bth0cs( 33) / -.8951804225 8785778806 1439459536 43 d-27    /
      data bth0cs( 34) / +.3235723747 9762298533 2562358685 87 d-27    /
      data bth0cs( 35) / -.1188301051 9855353657 0471441137 96 d-27    /
      data bth0cs( 36) / +.4430628690 7358104820 5792319417 31 d-28    /
      data bth0cs( 37) / -.1676100964 8834829495 7920101356 81 d-28    /
      data bth0cs( 38) / +.6429294692 1207466972 5323939660 88 d-29    /
      data bth0cs( 39) / -.2499226116 6978652421 2072136827 63 d-29    /
      data bth0cs( 40) / +.9839979429 9521955672 8282603553 18 d-30    /
      data bth0cs( 41) / -.3922037524 2408016397 9891316261 58 d-30    /
      data bth0cs( 42) / +.1581810703 0056522138 5906188456 92 d-30    /
      data bth0cs( 43) / -.6452550614 4890715944 3440983654 26 d-31    /
      data bth0cs( 44) / +.2661111136 9199356137 1770183463 67 d-31    /
c
c series for bm02       on the interval  0.          to  1.56250e-02
c                                        with weighted error   4.72e-32
c                                         log weighted error  31.33
c                               significant figures required  30.00
c                                    decimal places required  32.13
c
      data bm02cs(  1) / +.9500415145 2283813693 3086133556 0 d-1      /
      data bm02cs(  2) / -.3801864682 3656709917 4808156685 1 d-3      /
      data bm02cs(  3) / +.2258339301 0314811929 5182992722 4 d-5      /
      data bm02cs(  4) / -.3895725802 3722287647 3062141260 5 d-7      /
      data bm02cs(  5) / +.1246886416 5120816979 3099052972 5 d-8      /
      data bm02cs(  6) / -.6065949022 1025037798 0383505838 7 d-10     /
      data bm02cs(  7) / +.4008461651 4217469910 1527597104 5 d-11     /
      data bm02cs(  8) / -.3350998183 3980942184 6729879457 4 d-12     /
      data bm02cs(  9) / +.3377119716 5174173670 6326434199 6 d-13     /
      data bm02cs( 10) / -.3964585901 6350127005 6935629582 3 d-14     /
      data bm02cs( 11) / +.5286111503 8838572173 8793974473 5 d-15     /
      data bm02cs( 12) / -.7852519083 4508523136 5464024349 3 d-16     /
      data bm02cs( 13) / +.1280300573 3866822010 1163407344 9 d-16     /
      data bm02cs( 14) / -.2263996296 3914297762 8709924488 4 d-17     /
      data bm02cs( 15) / +.4300496929 6567903886 4641029047 7 d-18     /
      data bm02cs( 16) / -.8705749805 1325870797 4753545145 5 d-19     /
      data bm02cs( 17) / +.1865862713 9620951411 8144277205 0 d-19     /
      data bm02cs( 18) / -.4210482486 0930654573 4508697230 1 d-20     /
      data bm02cs( 19) / +.9956676964 2284009915 8162741784 2 d-21     /
      data bm02cs( 20) / -.2457357442 8053133596 0592147854 7 d-21     /
      data bm02cs( 21) / +.6307692160 7620315680 8735370705 9 d-22     /
      data bm02cs( 22) / -.1678773691 4407401426 9333117238 8 d-22     /
      data bm02cs( 23) / +.4620259064 6739044337 7087813608 7 d-23     /
      data bm02cs( 24) / -.1311782266 8603087322 3769340249 6 d-23     /
      data bm02cs( 25) / +.3834087564 1163028277 4792244027 6 d-24     /
      data bm02cs( 26) / -.1151459324 0777412710 7261329357 6 d-24     /
      data bm02cs( 27) / +.3547210007 5233385230 7697134521 3 d-25     /
      data bm02cs( 28) / -.1119218385 8150046462 6435594217 6 d-25     /
      data bm02cs( 29) / +.3611879427 6298378316 9840499425 7 d-26     /
      data bm02cs( 30) / -.1190687765 9133331500 9264176246 3 d-26     /
      data bm02cs( 31) / +.4005094059 4039681318 0247644953 6 d-27     /
      data bm02cs( 32) / -.1373169422 4522123905 9519391601 7 d-27     /
      data bm02cs( 33) / +.4794199088 7425315859 9649152643 7 d-28     /
      data bm02cs( 34) / -.1702965627 6241095840 0699447645 2 d-28     /
      data bm02cs( 35) / +.6149512428 9363300715 0357516132 4 d-29     /
      data bm02cs( 36) / -.2255766896 5818283499 4430023724 2 d-29     /
      data bm02cs( 37) / +.8399707509 2942994860 6165835320 0 d-30     /
      data bm02cs( 38) / -.3172997595 5626023555 6742393615 2 d-30     /
      data bm02cs( 39) / +.1215205298 8812985545 8333302651 4 d-30     /
      data bm02cs( 40) / -.4715852749 7544386930 1321056804 5 d-31     /
c
c series for bt02       on the interval  1.56250e-02 to  6.25000e-02
c                                        with weighted error   2.99e-32
c                                         log weighted error  31.52
c                               significant figures required  30.61
c                                    decimal places required  32.32
c
      data bt02cs(  1) / -.2454829521 3424597462 0504672493 24 d+0     /
      data bt02cs(  2) / +.1254412103 9084615780 7853317782 99 d-2     /
      data bt02cs(  3) / -.3125395041 4871522854 9734467095 71 d-4     /
      data bt02cs(  4) / +.1470977824 9940831164 4534269693 14 d-5     /
      data bt02cs(  5) / -.9954348893 7950033643 4688503511 58 d-7     /
      data bt02cs(  6) / +.8549316673 3203041247 5787113977 51 d-8     /
      data bt02cs(  7) / -.8698975952 6554334557 9855121791 92 d-9     /
      data bt02cs(  8) / +.1005209953 3559791084 5401010821 53 d-9     /
      data bt02cs(  9) / -.1282823060 1708892903 4836236855 44 d-10    /
      data bt02cs( 10) / +.1773170078 1805131705 6557504510 23 d-11    /
      data bt02cs( 11) / -.2617457456 9485577488 6362841809 25 d-12    /
      data bt02cs( 12) / +.4082835138 9972059621 9664812211 03 d-13    /
      data bt02cs( 13) / -.6675166823 9742720054 6067495542 61 d-14    /
      data bt02cs( 14) / +.1136576139 3071629448 3924695499 51 d-14    /
      data bt02cs( 15) / -.2005118962 0647160250 5592664121 17 d-15    /
      data bt02cs( 16) / +.3649797879 4766269635 7205914641 06 d-16    /
      data bt02cs( 17) / -.6830963756 4582303169 3558437888 00 d-17    /
      data bt02cs( 18) / +.1310758314 5670756620 0571042679 46 d-17    /
      data bt02cs( 19) / -.2572336310 1850607778 7571306495 99 d-18    /
      data bt02cs( 20) / +.5152165744 1863959925 2677809493 33 d-19    /
      data bt02cs( 21) / -.1051301756 3758802637 9407414613 33 d-19    /
      data bt02cs( 22) / +.2182038199 1194813847 3010845013 33 d-20    /
      data bt02cs( 23) / -.4600470121 0362160577 2259054933 33 d-21    /
      data bt02cs( 24) / +.9840700692 5466818520 9536511999 99 d-22    /
      data bt02cs( 25) / -.2133403803 5728375844 7359863466 66 d-22    /
      data bt02cs( 26) / +.4683103642 3973365296 0662869333 33 d-23    /
      data bt02cs( 27) / -.1040021369 1985747236 5133823999 99 d-23    /
      data bt02cs( 28) / +.2334910567 7301510051 7777408000 00 d-24    /
      data bt02cs( 29) / -.5295682532 3318615788 0497493333 33 d-25    /
      data bt02cs( 30) / +.1212634195 2959756829 1962879999 99 d-25    /
      data bt02cs( 31) / -.2801889708 2289428760 2756266666 66 d-26    /
      data bt02cs( 32) / +.6529267898 7012873342 5937066666 66 d-27    /
      data bt02cs( 33) / -.1533798006 1873346427 8357333333 33 d-27    /
      data bt02cs( 34) / +.3630588430 6364536682 3594666666 66 d-28    /
      data bt02cs( 35) / -.8656075571 3629122479 1722666666 66 d-29    /
      data bt02cs( 36) / +.2077990997 2536284571 2383999999 99 d-29    /
      data bt02cs( 37) / -.5021117022 1417221674 3253333333 33 d-30    /
      data bt02cs( 38) / +.1220836027 9441714184 1919999999 99 d-30    /
      data bt02cs( 39) / -.2986005626 7039913454 2506666666 66 d-31    /
c
      data pi4 / 0.7853981633 9744830961 5660845819 876 d0 /
      data nbm0, nbt02, nbm02, nbth0, xmax / 4*0, 0.d0 /
c
      if (nbm0.ne.0) go to 10
      eta = 0.1*sngl(d1mach(3))
      nbm0 = initds (bm0cs, 37, eta)
      nbt02 = initds (bt02cs, 39, eta)
      nbm02 = initds (bm02cs, 40, eta)
      nbth0 = initds (bth0cs, 44, eta)
c
      xmax = 1.0d0/d1mach(4)
c
 10   if (x.lt.4.d0) call seteru (22hd9b0mp  x must be ge 4, 22, 1, 2)
c
      if (x.gt.8.d0) go to 20
      z = (128.d0/(x*x) - 5.d0)/3.d0
      ampl = (.75d0 + dcsevl (z, bm0cs, nbm0))/dsqrt(x)
      theta = x - pi4 + dcsevl (z, bt02cs, nbt02)/x
      return
c
 20   if (x.gt.xmax) call seteru (
     1  37hd9b0mp  no precision because x is big, 37, 2, 2)
c
      z = 128.d0/(x*x) - 1.d0
      ampl = (.75d0 + dcsevl (z, bm02cs, nbm02))/dsqrt(x)
      theta = x - pi4 + dcsevl (z, bth0cs, nbth0)/x
      return
c
      end
      subroutine d9b1mp (x, ampl, theta)
c july 1977 edition.  w. fullerton, c3, los alamos scientific lab.
c
c evaluate the modulus and phase for the bessel j1 and y1 functions.
c
      double precision x, ampl, theta, bm1cs(37), bt12cs(39),
     1  bm12cs(40), bth1cs(44), xmax, pi4, z, d1mach, dcsevl,
     2  dsqrt
c     external d1mach, dcsevl, dsqrt, initds
c
c series for bm1        on the interval  1.56250e-02 to  6.25000e-02
c                                        with weighted error   4.91e-32
c                                         log weighted error  31.31
c                               significant figures required  30.04
c                                    decimal places required  32.09
c
      data bm1 cs(  1) / +.1069845452 6180630149 6998530853 8 d+0      /
      data bm1 cs(  2) / +.3274915039 7159649007 2905514344 5 d-2      /
      data bm1 cs(  3) / -.2987783266 8316985920 3044577793 8 d-4      /
      data bm1 cs(  4) / +.8331237177 9919745313 9322266902 3 d-6      /
      data bm1 cs(  5) / -.4112665690 3020073048 9638172549 8 d-7      /
      data bm1 cs(  6) / +.2855344228 7892152207 1975766316 1 d-8      /
      data bm1 cs(  7) / -.2485408305 4156238780 6002659605 5 d-9      /
      data bm1 cs(  8) / +.2543393338 0725824427 4248439717 4 d-10     /
      data bm1 cs(  9) / -.2941045772 8229675234 8975082790 9 d-11     /
      data bm1 cs( 10) / +.3743392025 4939033092 6505615362 6 d-12     /
      data bm1 cs( 11) / -.5149118293 8211672187 2054824352 7 d-13     /
      data bm1 cs( 12) / +.7552535949 8651439080 3404076419 9 d-14     /
      data bm1 cs( 13) / -.1169409706 8288464441 6629062246 4 d-14     /
      data bm1 cs( 14) / +.1896562449 4347915717 2182460506 0 d-15     /
      data bm1 cs( 15) / -.3201955368 6932864206 6477531639 4 d-16     /
      data bm1 cs( 16) / +.5599548399 3162041144 8416990549 3 d-17     /
      data bm1 cs( 17) / -.1010215894 7304324431 1939044454 4 d-17     /
      data bm1 cs( 18) / +.1873844985 7275629833 0204271957 3 d-18     /
      data bm1 cs( 19) / -.3563537470 3285802192 7430143999 9 d-19     /
      data bm1 cs( 20) / +.6931283819 9712383304 2276351999 9 d-20     /
      data bm1 cs( 21) / -.1376059453 4065001522 5140893013 3 d-20     /
      data bm1 cs( 22) / +.2783430784 1070802205 9977932799 9 d-21     /
      data bm1 cs( 23) / -.5727595364 3205616893 4866943999 9 d-22     /
      data bm1 cs( 24) / +.1197361445 9188926725 3575679999 9 d-22     /
      data bm1 cs( 25) / -.2539928509 8918719766 4144042666 6 d-23     /
      data bm1 cs( 26) / +.5461378289 6572959730 6961919999 9 d-24     /
      data bm1 cs( 27) / -.1189211341 7733202889 8628949333 3 d-24     /
      data bm1 cs( 28) / +.2620150977 3400815949 5782400000 0 d-25     /
      data bm1 cs( 29) / -.5836810774 2556859019 2093866666 6 d-26     /
      data bm1 cs( 30) / +.1313743500 0805957734 2361599999 9 d-26     /
      data bm1 cs( 31) / -.2985814622 5103803553 3277866666 6 d-27     /
      data bm1 cs( 32) / +.6848390471 3346049376 2559999999 9 d-28     /
      data bm1 cs( 33) / -.1584401568 2224767211 9296000000 0 d-28     /
      data bm1 cs( 34) / +.3695641006 5709380543 0101333333 3 d-29     /
      data bm1 cs( 35) / -.8687115921 1446682430 1226666666 6 d-30     /
      data bm1 cs( 36) / +.2057080846 1587634629 2906666666 6 d-30     /
      data bm1 cs( 37) / -.4905225761 1162255185 2373333333 3 d-31     /
c
c series for bt12       on the interval  1.56250e-02 to  6.25000e-02
c                                        with weighted error   3.33e-32
c                                         log weighted error  31.48
c                               significant figures required  31.05
c                                    decimal places required  32.27
c
      data bt12cs(  1) / +.7382386012 8742974662 6208397927 64 d+0     /
      data bt12cs(  2) / -.3336111317 4483906384 4701476811 89 d-2     /
      data bt12cs(  3) / +.6146345488 8046964698 5148994201 86 d-4     /
      data bt12cs(  4) / -.2402458516 1602374264 9776354695 68 d-5     /
      data bt12cs(  5) / +.1466355557 7509746153 2105919972 04 d-6     /
      data bt12cs(  6) / -.1184191730 5589180567 0051475049 83 d-7     /
      data bt12cs(  7) / +.1157419896 3919197052 1254663030 55 d-8     /
      data bt12cs(  8) / -.1300116112 9439187449 3660077945 71 d-9     /
      data bt12cs(  9) / +.1624539114 1361731937 7421662736 67 d-10    /
      data bt12cs( 10) / -.2208963682 1403188752 1554417701 28 d-11    /
      data bt12cs( 11) / +.3218030425 8553177090 4743586537 78 d-12    /
      data bt12cs( 12) / -.4965314793 2768480785 5520211353 81 d-13    /
      data bt12cs( 13) / +.8043890043 2847825985 5588826393 17 d-14    /
      data bt12cs( 14) / -.1358912131 0161291384 6947126822 82 d-14    /
      data bt12cs( 15) / +.2381050439 7147214869 6765296059 73 d-15    /
      data bt12cs( 16) / -.4308146636 3849106724 4712414207 99 d-16    /
      data bt12cs( 17) / +.8020254403 2771002434 9935125504 00 d-17    /
      data bt12cs( 18) / -.1531631064 2462311864 2300274687 99 d-17    /
      data bt12cs( 19) / +.2992860635 2715568924 0730405546 66 d-18    /
      data bt12cs( 20) / -.5970996465 8085443393 8156366506 66 d-19    /
      data bt12cs( 21) / +.1214028966 9415185024 1608526506 66 d-19    /
      data bt12cs( 22) / -.2511511469 6612948901 0069777066 66 d-20    /
      data bt12cs( 23) / +.5279056717 0328744850 7383807999 99 d-21    /
      data bt12cs( 24) / -.1126050922 7550498324 3611613866 66 d-21    /
      data bt12cs( 25) / +.2434827735 9576326659 6634624000 00 d-22    /
      data bt12cs( 26) / -.5331726123 6931800130 0384426666 66 d-23    /
      data bt12cs( 27) / +.1181361505 9707121039 2059903999 99 d-23    /
      data bt12cs( 28) / -.2646536828 3353523514 8567893333 33 d-24    /
      data bt12cs( 29) / +.5990339404 1361503945 5778133333 33 d-25    /
      data bt12cs( 30) / -.1369085463 0829503109 1363839999 99 d-25    /
      data bt12cs( 31) / +.3157679015 4380228326 4136533333 33 d-26    /
      data bt12cs( 32) / -.7345791508 2084356491 4005333333 33 d-27    /
      data bt12cs( 33) / +.1722808148 0722747930 7059200000 00 d-27    /
      data bt12cs( 34) / -.4071690796 1286507941 0688000000 00 d-28    /
      data bt12cs( 35) / +.9693474513 6779622700 3733333333 33 d-29    /
      data bt12cs( 36) / -.2323763633 7765716765 3546666666 66 d-29    /
      data bt12cs( 37) / +.5607451067 3522029406 8906666666 66 d-30    /
      data bt12cs( 38) / -.1361646539 1539005860 5226666666 66 d-30    /
      data bt12cs( 39) / +.3326310923 3894654388 9066666666 66 d-31    /
c
c series for bm12       on the interval  0.          to  1.56250e-02
c                                        with weighted error   5.01e-32
c                                         log weighted error  31.30
c                               significant figures required  29.99
c                                    decimal places required  32.10
c
      data bm12cs(  1) / +.9807979156 2330500272 7209354693 7 d-1      /
      data bm12cs(  2) / +.1150961189 5046853061 7548348460 2 d-2      /
      data bm12cs(  3) / -.4312482164 3382054098 8935809773 2 d-5      /
      data bm12cs(  4) / +.5951839610 0888163078 1302980183 2 d-7      /
      data bm12cs(  5) / -.1704844019 8269098574 0070158647 8 d-8      /
      data bm12cs(  6) / +.7798265413 6111095086 5817382740 1 d-10     /
      data bm12cs(  7) / -.4958986126 7664158094 9175495186 5 d-11     /
      data bm12cs(  8) / +.4038432416 4211415168 3820226514 4 d-12     /
      data bm12cs(  9) / -.3993046163 7251754457 6548384664 5 d-13     /
      data bm12cs( 10) / +.4619886183 1189664943 1334243277 5 d-14     /
      data bm12cs( 11) / -.6089208019 0953833013 4547261933 3 d-15     /
      data bm12cs( 12) / +.8960930916 4338764821 5704804124 9 d-16     /
      data bm12cs( 13) / -.1449629423 9420231229 1651891892 5 d-16     /
      data bm12cs( 14) / +.2546463158 5377760561 6514964806 8 d-17     /
      data bm12cs( 15) / -.4809472874 6478364442 5926371862 0 d-18     /
      data bm12cs( 16) / +.9687684668 2925990490 8727583912 4 d-19     /
      data bm12cs( 17) / -.2067213372 2779660232 4503811755 1 d-19     /
      data bm12cs( 18) / +.4646651559 1503847318 0276780959 0 d-20     /
      data bm12cs( 19) / -.1094966128 8483341382 4135132833 9 d-20     /
      data bm12cs( 20) / +.2693892797 2886828609 0570761278 5 d-21     /
      data bm12cs( 21) / -.6894992910 9303744778 1897002685 7 d-22     /
      data bm12cs( 22) / +.1830268262 7520629098 9066855474 0 d-22     /
      data bm12cs( 23) / -.5025064246 3519164281 5611355322 4 d-23     /
      data bm12cs( 24) / +.1423545194 4548060396 3169363419 4 d-23     /
      data bm12cs( 25) / -.4152191203 6164503880 6888676980 1 d-24     /
      data bm12cs( 26) / +.1244609201 5039793258 8233007654 7 d-24     /
      data bm12cs( 27) / -.3827336370 5693042994 3191866128 6 d-25     /
      data bm12cs( 28) / +.1205591357 8156175353 7472398183 5 d-25     /
      data bm12cs( 29) / -.3884536246 3764880764 3185936112 4 d-26     /
      data bm12cs( 30) / +.1278689528 7204097219 0489528346 1 d-26     /
      data bm12cs( 31) / -.4295146689 4479462720 6193691591 2 d-27     /
      data bm12cs( 32) / +.1470689117 8290708864 5680270798 3 d-27     /
      data bm12cs( 33) / -.5128315665 1060731281 8037401779 6 d-28     /
      data bm12cs( 34) / +.1819509585 4711693854 8143737328 6 d-28     /
      data bm12cs( 35) / -.6563031314 8419808676 1863505037 3 d-29     /
      data bm12cs( 36) / +.2404898976 9199606531 9891487583 4 d-29     /
      data bm12cs( 37) / -.8945966744 6906124732 3495824297 9 d-30     /
      data bm12cs( 38) / +.3376085160 6572310266 3714897824 0 d-30     /
      data bm12cs( 39) / -.1291791454 6206563609 1309991696 6 d-30     /
      data bm12cs( 40) / +.5008634462 9588105206 8495150125 4 d-31     /
c
c series for bth1       on the interval  0.          to  1.56250e-02
c                                        with weighted error   2.82e-32
c                                         log weighted error  31.55
c                               significant figures required  31.12
c                                    decimal places required  32.37
c
      data bth1cs(  1) / +.7474995720 3587276055 4434839696 95 d+0     /
      data bth1cs(  2) / -.1240077714 4651711252 5457775413 84 d-2     /
      data bth1cs(  3) / +.9925244240 4424527376 6414976895 92 d-5     /
      data bth1cs(  4) / -.2030369073 7159711052 4193753756 08 d-6     /
      data bth1cs(  5) / +.7535961770 5690885712 1840175836 29 d-8     /
      data bth1cs(  6) / -.4166161271 5343550107 6300238562 28 d-9     /
      data bth1cs(  7) / +.3070161807 0834890481 2451020912 16 d-10    /
      data bth1cs(  8) / -.2817849963 7605213992 3240088839 24 d-11    /
      data bth1cs(  9) / +.3079069673 9040295476 0281468216 47 d-12    /
      data bth1cs( 10) / -.3880330026 2803434112 7873475547 81 d-13    /
      data bth1cs( 11) / +.5509603960 8630904934 5617262085 62 d-14    /
      data bth1cs( 12) / -.8659006076 8383779940 1033989539 94 d-15    /
      data bth1cs( 13) / +.1485604914 1536749003 4236890606 83 d-15    /
      data bth1cs( 14) / -.2751952981 5904085805 3712121250 09 d-16    /
      data bth1cs( 15) / +.5455079609 0481089625 0362236409 23 d-17    /
      data bth1cs( 16) / -.1148653450 1983642749 5436310271 77 d-17    /
      data bth1cs( 17) / +.2553521337 7973900223 1990525335 22 d-18    /
      data bth1cs( 18) / -.5962149019 7413450395 7682879078 49 d-19    /
      data bth1cs( 19) / +.1455662290 2372718620 2883020058 33 d-19    /
      data bth1cs( 20) / -.3702218542 2450538201 5797760195 93 d-20    /
      data bth1cs( 21) / +.9776307412 5345357664 1684345179 24 d-21    /
      data bth1cs( 22) / -.2672682163 9668488468 7237753930 52 d-21    /
      data bth1cs( 23) / +.7545330038 4983271794 0381906557 64 d-22    /
      data bth1cs( 24) / -.2194789991 9802744897 8923833716 47 d-22    /
      data bth1cs( 25) / +.6564839462 3955262178 9069998174 93 d-23    /
      data bth1cs( 26) / -.2015560429 8370207570 7840768695 19 d-23    /
      data bth1cs( 27) / +.6341776855 6776143492 1446671856 70 d-24    /
      data bth1cs( 28) / -.2041927788 5337895634 8137699555 91 d-24    /
      data bth1cs( 29) / +.6719146422 0720567486 6589800185 51 d-25    /
      data bth1cs( 30) / -.2256907911 0207573595 7090036873 36 d-25    /
      data bth1cs( 31) / +.7729771989 2989706370 9269598719 29 d-26    /
      data bth1cs( 32) / -.2696744451 2294640913 2114240809 20 d-26    /
      data bth1cs( 33) / +.9574934451 8502698072 2955219336 27 d-27    /
      data bth1cs( 34) / -.3456916844 8890113000 1756808276 27 d-27    /
      data bth1cs( 35) / +.1268123481 7398436504 2119862383 74 d-27    /
      data bth1cs( 36) / -.4723253663 0722639860 4649937134 45 d-28    /
      data bth1cs( 37) / +.1785000847 8186376177 8586197964 17 d-28    /
      data bth1cs( 38) / -.6840436100 4510395406 2152235667 46 d-29    /
      data bth1cs( 39) / +.2656602867 1720419358 2934226722 12 d-29    /
      data bth1cs( 40) / -.1045040252 7914452917 7141614846 70 d-29    /
      data bth1cs( 41) / +.4161829082 5377144306 8619171970 64 d-30    /
      data bth1cs( 42) / -.1677163920 3643714856 5013478828 87 d-30    /
      data bth1cs( 43) / +.6836199777 6664389173 5359280285 28 d-31    /
      data bth1cs( 44) / -.2817224786 1233641166 7395746228 10 d-31    /
c
      data pi4 / 0.7853981633 9744830961 5660845819 876 d0 /
      data nbm1, nbt12, nbm12, nbth1, xmax / 4*0, 0.d0 /
c
      if (nbm1.ne.0) go to 10
      eta = 0.1*sngl(d1mach(3))
      nbm1 = initds (bm1cs, 37, eta)
      nbt12 = initds (bt12cs, 39, eta)
      nbm12 = initds (bm12cs, 40, eta)
      nbth1 = initds (bth1cs, 44, eta)
c
      xmax = 1.0d0/d1mach(4)
c
 10   if (x.lt.4.d0) call seteru (22hd9b1mp  x must be ge 4, 22, 1, 2)
c
      if (x.gt.8.d0) go to 20
      z = (128.d0/(x*x) - 5.d0)/3.d0
      ampl = (0.75d0 + dcsevl (z, bm1cs, nbm1))/dsqrt(x)
      theta = x - 3.d0*pi4 + dcsevl (z, bt12cs, nbt12)/x
      return
c
 20   if (x.gt.xmax) call seteru (
     1  37hd9b1mp  no precision because x is big, 37, 2, 2)
c
      z = 128.d0/(x*x) - 1.d0
      ampl = (0.75d0 + dcsevl (z, bm12cs, nbm12))/dsqrt(x)
      theta = x - 3.d0*pi4 + dcsevl (z, bth1cs, nbth1)/x
      return
c
      end
      double precision function dbsk0e (x)
c july 1980 edition.  w. fullerton, c3, los alamos scientific lab.
      double precision x, bk0cs(16), ak0cs(38), ak02cs(33),
     1  xsml, y, d1mach, dcsevl, dbesi0, dexp, dlog, dsqrt
c     external d1mach, dbesi0, dcsevl, dexp, dlog, dsqrt, initds
c
c series for bk0        on the interval  0.          to  4.00000e+00
c                                        with weighted error   3.08e-33
c                                         log weighted error  32.51
c                               significant figures required  32.05
c                                    decimal places required  33.11
c
      data bk0 cs(  1) / -.3532739323 3902768720 1140060063 153 d-1    /
      data bk0 cs(  2) / +.3442898999 2462848688 6344927529 213 d+0    /
      data bk0 cs(  3) / +.3597993651 5361501626 5721303687 231 d-1    /
      data bk0 cs(  4) / +.1264615411 4469259233 8479508673 447 d-2    /
      data bk0 cs(  5) / +.2286212103 1194517860 8269830297 585 d-4    /
      data bk0 cs(  6) / +.2534791079 0261494573 0790013428 354 d-6    /
      data bk0 cs(  7) / +.1904516377 2202088589 7214059381 366 d-8    /
      data bk0 cs(  8) / +.1034969525 7633624585 1008317853 089 d-10   /
      data bk0 cs(  9) / +.4259816142 7910825765 2445327170 133 d-13   /
      data bk0 cs( 10) / +.1374465435 8807508969 4238325440 000 d-15   /
      data bk0 cs( 11) / +.3570896528 5083735909 9688597333 333 d-18   /
      data bk0 cs( 12) / +.7631643660 1164373766 7498666666 666 d-21   /
      data bk0 cs( 13) / +.1365424988 4407818590 8053333333 333 d-23   /
      data bk0 cs( 14) / +.2075275266 9066680831 9999999999 999 d-26   /
      data bk0 cs( 15) / +.2712814218 0729856000 0000000000 000 d-29   /
      data bk0 cs( 16) / +.3082593887 9146666666 6666666666 666 d-32   /
c
c series for ak0        on the interval  1.25000e-01 to  5.00000e-01
c                                        with weighted error   2.85e-32
c                                         log weighted error  31.54
c                               significant figures required  30.19
c                                    decimal places required  32.33
c
      data ak0 cs(  1) / -.7643947903 3279414240 8297827008 8 d-1      /
      data ak0 cs(  2) / -.2235652605 6998190520 2309555079 1 d-1      /
      data ak0 cs(  3) / +.7734181154 6938582353 0061817404 7 d-3      /
      data ak0 cs(  4) / -.4281006688 8860994644 5214643541 6 d-4      /
      data ak0 cs(  5) / +.3081700173 8629747436 5001482666 0 d-5      /
      data ak0 cs(  6) / -.2639367222 0096649740 6744889272 3 d-6      /
      data ak0 cs(  7) / +.2563713036 4034692062 9408826574 2 d-7      /
      data ak0 cs(  8) / -.2742705549 9002012638 5721191524 4 d-8      /
      data ak0 cs(  9) / +.3169429658 0974995920 8083287340 3 d-9      /
      data ak0 cs( 10) / -.3902353286 9621841416 0106571796 2 d-10     /
      data ak0 cs( 11) / +.5068040698 1885754020 5009212728 6 d-11     /
      data ak0 cs( 12) / -.6889574741 0078706795 4171355798 4 d-12     /
      data ak0 cs( 13) / +.9744978497 8259176913 8820133683 1 d-13     /
      data ak0 cs( 14) / -.1427332841 8845485053 8985534012 2 d-13     /
      data ak0 cs( 15) / +.2156412571 0214630395 5806297652 7 d-14     /
      data ak0 cs( 16) / -.3349654255 1495627721 8878205853 0 d-15     /
      data ak0 cs( 17) / +.5335260216 9529116921 4528039260 1 d-16     /
      data ak0 cs( 18) / -.8693669980 8907538076 3962237883 7 d-17     /
      data ak0 cs( 19) / +.1446404347 8622122278 8776344234 6 d-17     /
      data ak0 cs( 20) / -.2452889825 5001296824 0467875157 3 d-18     /
      data ak0 cs( 21) / +.4233754526 2321715728 2170634240 0 d-19     /
      data ak0 cs( 22) / -.7427946526 4544641956 9534129493 3 d-20     /
      data ak0 cs( 23) / +.1323150529 3926668662 7796746240 0 d-20     /
      data ak0 cs( 24) / -.2390587164 7396494513 3598146559 9 d-21     /
      data ak0 cs( 25) / +.4376827585 9232261401 6571255466 6 d-22     /
      data ak0 cs( 26) / -.8113700607 3451180593 3901141333 3 d-23     /
      data ak0 cs( 27) / +.1521819913 8321729583 1037815466 6 d-23     /
      data ak0 cs( 28) / -.2886041941 4833977702 3595861333 3 d-24     /
      data ak0 cs( 29) / +.5530620667 0547179799 9261013333 3 d-25     /
      data ak0 cs( 30) / -.1070377329 2498987285 9163306666 6 d-25     /
      data ak0 cs( 31) / +.2091086893 1423843002 9632853333 3 d-26     /
      data ak0 cs( 32) / -.4121713723 6462038274 1026133333 3 d-27     /
      data ak0 cs( 33) / +.8193483971 1213076401 3568000000 0 d-28     /
      data ak0 cs( 34) / -.1642000275 4592977267 8075733333 3 d-28     /
      data ak0 cs( 35) / +.3316143281 4802271958 9034666666 6 d-29     /
      data ak0 cs( 36) / -.6746863644 1452959410 8586666666 6 d-30     /
      data ak0 cs( 37) / +.1382429146 3184246776 3541333333 3 d-30     /
      data ak0 cs( 38) / -.2851874167 3598325708 1173333333 3 d-31     /
c
c series for ak02       on the interval  0.          to  1.25000e-01
c                                        with weighted error   2.30e-32
c                                         log weighted error  31.64
c                               significant figures required  29.68
c                                    decimal places required  32.40
c
      data ak02cs(  1) / -.1201869826 3075922398 3934621245 2 d-1      /
      data ak02cs(  2) / -.9174852691 0256953106 5256107571 3 d-2      /
      data ak02cs(  3) / +.1444550931 7750058210 4884387805 7 d-3      /
      data ak02cs(  4) / -.4013614175 4357097286 7102107787 9 d-5      /
      data ak02cs(  5) / +.1567831810 8523106725 9034899033 3 d-6      /
      data ak02cs(  6) / -.7770110438 5217377103 1579975446 0 d-8      /
      data ak02cs(  7) / +.4611182576 1797178825 3313052958 6 d-9      /
      data ak02cs(  8) / -.3158592997 8605657705 2666580330 9 d-10     /
      data ak02cs(  9) / +.2435018039 3650411278 3588781432 9 d-11     /
      data ak02cs( 10) / -.2074331387 3983478977 0985337350 6 d-12     /
      data ak02cs( 11) / +.1925787280 5899170847 4273650469 3 d-13     /
      data ak02cs( 12) / -.1927554805 8389561036 0034718221 8 d-14     /
      data ak02cs( 13) / +.2062198029 1978182782 8523786964 4 d-15     /
      data ak02cs( 14) / -.2341685117 5792424026 0364019507 1 d-16     /
      data ak02cs( 15) / +.2805902810 6430422468 1517882845 8 d-17     /
      data ak02cs( 16) / -.3530507631 1618079458 1548246357 3 d-18     /
      data ak02cs( 17) / +.4645295422 9351082674 2421633706 6 d-19     /
      data ak02cs( 18) / -.6368625941 3442664739 2205346133 3 d-20     /
      data ak02cs( 19) / +.9069521310 9865155676 2234880000 0 d-21     /
      data ak02cs( 20) / -.1337974785 4236907398 4500531199 9 d-21     /
      data ak02cs( 21) / +.2039836021 8599523155 2208896000 0 d-22     /
      data ak02cs( 22) / -.3207027481 3678405000 6086997333 3 d-23     /
      data ak02cs( 23) / +.5189744413 6623099636 2635946666 6 d-24     /
      data ak02cs( 24) / -.8629501497 5405721929 6460799999 9 d-25     /
      data ak02cs( 25) / +.1472161183 1025598552 0803840000 0 d-25     /
      data ak02cs( 26) / -.2573069023 8670112838 1235199999 9 d-26     /
      data ak02cs( 27) / +.4601774086 6435165873 7664000000 0 d-27     /
      data ak02cs( 28) / -.8411555324 2010937371 3066666666 6 d-28     /
      data ak02cs( 29) / +.1569806306 6353689393 0154666666 6 d-28     /
      data ak02cs( 30) / -.2988226453 0057577889 7919999999 9 d-29     /
      data ak02cs( 31) / +.5796831375 2168365206 1866666666 6 d-30     /
      data ak02cs( 32) / -.1145035994 3476813321 5573333333 3 d-30     /
      data ak02cs( 33) / +.2301266594 2496828020 0533333333 3 d-31     /
      data ntk0, ntak0, ntak02, xsml / 3*0, 0.0d0 /
c
      if (ntk0.ne.0) go to 10
      eta = 0.1*sngl(d1mach(3))
      ntk0 = initds (bk0cs, 16, eta)
      ntak0 = initds (ak0cs, 38, eta)
      ntak02 = initds (ak02cs, 33, eta)
      xsml = dsqrt (4.0d0*d1mach(3))
c
 10   if (x.le.0.d0) call seteru (29hdbsk0e  x is zero or negative, 29,
     1  1, 2)
      if (x.gt.2.0d0) go to 20
c
      y = 0.d0
      if (x.gt.xsml) y = x*x
      dbsk0e = dexp(x)*(-dlog(0.5d0*x)*dbesi0(x) - 0.25d0 +
     1  dcsevl (.5d0*y-1.d0, bk0cs, ntk0))
      return
c
 20   if (x.le.8.d0) dbsk0e = (1.25d0 + dcsevl ((16.d0/x-5.d0)/3.d0,
     1  ak0cs, ntak0))/dsqrt(x)
      if (x.gt.8.d0) dbsk0e = (1.25d0 +
     1  dcsevl (16.d0/x-1.d0, ak02cs, ntak02))/dsqrt(x)
c
      return
      end
      double precision function dbsk1e (x)
c july 1980 edition.  w. fullerton, c3, los alamos scientific lab.
      double precision x, bk1cs(16), ak1cs(38), ak12cs(33), xmin,
     1  xsml, y, d1mach, dcsevl, dbesi1, dexp, dlog, dsqrt
c     external d1mach, dbesi1, dcsevl, dexp, dlog, dsqrt, initds
c
c series for bk1        on the interval  0.          to  4.00000e+00
c                                        with weighted error   9.16e-32
c                                         log weighted error  31.04
c                               significant figures required  30.61
c                                    decimal places required  31.64
c
      data bk1 cs(  1) / +.2530022733 8947770532 5311208685 33 d-1     /
      data bk1 cs(  2) / -.3531559607 7654487566 7238316918 01 d+0     /
      data bk1 cs(  3) / -.1226111808 2265714823 4790679300 42 d+0     /
      data bk1 cs(  4) / -.6975723859 6398643501 8129202960 83 d-2     /
      data bk1 cs(  5) / -.1730288957 5130520630 1765073689 79 d-3     /
      data bk1 cs(  6) / -.2433406141 5659682349 6007350301 64 d-5     /
      data bk1 cs(  7) / -.2213387630 7347258558 3152525451 26 d-7     /
      data bk1 cs(  8) / -.1411488392 6335277610 9583302126 08 d-9     /
      data bk1 cs(  9) / -.6666901694 1993290060 8537512643 73 d-12    /
      data bk1 cs( 10) / -.2427449850 5193659339 2631968648 53 d-14    /
      data bk1 cs( 11) / -.7023863479 3862875971 7837971200 00 d-17    /
      data bk1 cs( 12) / -.1654327515 5100994675 4910293333 33 d-19    /
      data bk1 cs( 13) / -.3233834745 9944491991 8933333333 33 d-22    /
      data bk1 cs( 14) / -.5331275052 9265274999 4666666666 66 d-25    /
      data bk1 cs( 15) / -.7513040716 2157226666 6666666666 66 d-28    /
      data bk1 cs( 16) / -.9155085717 6541866666 6666666666 66 d-31    /
c
c series for ak1        on the interval  1.25000e-01 to  5.00000e-01
c                                        with weighted error   3.07e-32
c                                         log weighted error  31.51
c                               significant figures required  30.71
c                                    decimal places required  32.30
c
      data ak1 cs(  1) / +.2744313406 9738829695 2576662272 66 d+0     /
      data ak1 cs(  2) / +.7571989953 1993678170 8923781492 90 d-1     /
      data ak1 cs(  3) / -.1441051556 4754061229 8531161756 25 d-2     /
      data ak1 cs(  4) / +.6650116955 1257479394 2513854770 36 d-4     /
      data ak1 cs(  5) / -.4369984709 5201407660 5808450891 67 d-5     /
      data ak1 cs(  6) / +.3540277499 7630526799 4171390085 34 d-6     /
      data ak1 cs(  7) / -.3311163779 2932920208 9826882457 04 d-7     /
      data ak1 cs(  8) / +.3445977581 9010534532 3114997709 92 d-8     /
      data ak1 cs(  9) / -.3898932347 4754271048 9819374927 58 d-9     /
      data ak1 cs( 10) / +.4720819750 4658356400 9474493390 05 d-10    /
      data ak1 cs( 11) / -.6047835662 8753562345 3735915628 90 d-11    /
      data ak1 cs( 12) / +.8128494874 8658747888 1938379856 63 d-12    /
      data ak1 cs( 13) / -.1138694574 7147891428 9239159510 42 d-12    /
      data ak1 cs( 14) / +.1654035840 8462282325 9729482050 90 d-13    /
      data ak1 cs( 15) / -.2480902567 7068848221 5160104405 33 d-14    /
      data ak1 cs( 16) / +.3829237890 7024096948 4292272991 57 d-15    /
      data ak1 cs( 17) / -.6064734104 0012418187 7682103773 86 d-16    /
      data ak1 cs( 18) / +.9832425623 2648616038 1940046506 66 d-17    /
      data ak1 cs( 19) / -.1628416873 8284380035 6666201156 26 d-17    /
      data ak1 cs( 20) / +.2750153649 6752623718 2841203370 66 d-18    /
      data ak1 cs( 21) / -.4728966646 3953250924 2810695680 00 d-19    /
      data ak1 cs( 22) / +.8268150002 8109932722 3920503466 66 d-20    /
      data ak1 cs( 23) / -.1468140513 6624956337 1939648853 33 d-20    /
      data ak1 cs( 24) / +.2644763926 9208245978 0858948266 66 d-21    /
      data ak1 cs( 25) / -.4829015756 4856387897 9698688000 00 d-22    /
      data ak1 cs( 26) / +.8929302074 3610130180 6563327999 99 d-23    /
      data ak1 cs( 27) / -.1670839716 8972517176 9977514666 66 d-23    /
      data ak1 cs( 28) / +.3161645603 4040694931 3686186666 66 d-24    /
      data ak1 cs( 29) / -.6046205531 2274989106 5064106666 66 d-25    /
      data ak1 cs( 30) / +.1167879894 2042732700 7184213333 33 d-25    /
      data ak1 cs( 31) / -.2277374158 2653996232 8678400000 00 d-26    /
      data ak1 cs( 32) / +.4481109730 0773675795 3058133333 33 d-27    /
      data ak1 cs( 33) / -.8893288476 9020194062 3360000000 00 d-28    /
      data ak1 cs( 34) / +.1779468001 8850275131 3920000000 00 d-28    /
      data ak1 cs( 35) / -.3588455596 7329095821 9946666666 66 d-29    /
      data ak1 cs( 36) / +.7290629049 2694257991 6799999999 99 d-30    /
      data ak1 cs( 37) / -.1491844984 5546227073 0240000000 00 d-30    /
      data ak1 cs( 38) / +.3073657387 2934276300 7999999999 99 d-31    /
c
c series for ak12       on the interval  0.          to  1.25000e-01
c                                        with weighted error   2.41e-32
c                                         log weighted error  31.62
c                               significant figures required  30.25
c                                    decimal places required  32.38
c
      data ak12cs(  1) / +.6379308343 7390010366 0048853410 2 d-1      /
      data ak12cs(  2) / +.2832887813 0497209358 3503028470 8 d-1      /
      data ak12cs(  3) / -.2475370673 9052503454 1454556673 2 d-3      /
      data ak12cs(  4) / +.5771972451 6072488204 7097662576 3 d-5      /
      data ak12cs(  5) / -.2068939219 5365483027 4553319655 2 d-6      /
      data ak12cs(  6) / +.9739983441 3818041803 0921309788 7 d-8      /
      data ak12cs(  7) / -.5585336140 3806249846 8889551112 9 d-9      /
      data ak12cs(  8) / +.3732996634 0461852402 2121285473 1 d-10     /
      data ak12cs(  9) / -.2825051961 0232254451 3506575492 8 d-11     /
      data ak12cs( 10) / +.2372019002 4841441736 4349695548 6 d-12     /
      data ak12cs( 11) / -.2176677387 9917539792 6830166793 8 d-13     /
      data ak12cs( 12) / +.2157914161 6160324539 3956268970 6 d-14     /
      data ak12cs( 13) / -.2290196930 7182692759 9155133815 4 d-15     /
      data ak12cs( 14) / +.2582885729 8232749619 1993956522 6 d-16     /
      data ak12cs( 15) / -.3076752641 2684631876 2109817344 0 d-17     /
      data ak12cs( 16) / +.3851487721 2804915970 9489684479 9 d-18     /
      data ak12cs( 17) / -.5044794897 6415289771 1728250880 0 d-19     /
      data ak12cs( 18) / +.6888673850 4185442370 1829222399 9 d-20     /
      data ak12cs( 19) / -.9775041541 9501183030 0213248000 0 d-21     /
      data ak12cs( 20) / +.1437416218 5238364610 0165973333 3 d-21     /
      data ak12cs( 21) / -.2185059497 3443473734 9973333333 3 d-22     /
      data ak12cs( 22) / +.3426245621 8092206316 4538880000 0 d-23     /
      data ak12cs( 23) / -.5531064394 2464082325 0124800000 0 d-24     /
      data ak12cs( 24) / +.9176601505 6859954037 8282666666 6 d-25     /
      data ak12cs( 25) / -.1562287203 6180249114 4874666666 6 d-25     /
      data ak12cs( 26) / +.2725419375 4843331323 4943999999 9 d-26     /
      data ak12cs( 27) / -.4865674910 0748279923 7802666666 6 d-27     /
      data ak12cs( 28) / +.8879388552 7235025873 5786666666 6 d-28     /
      data ak12cs( 29) / -.1654585918 0392575489 3653333333 3 d-28     /
      data ak12cs( 30) / +.3145111321 3578486743 0399999999 9 d-29     /
      data ak12cs( 31) / -.6092998312 1931276124 1600000000 0 d-30     /
      data ak12cs( 32) / +.1202021939 3698158346 2399999999 9 d-30     /
      data ak12cs( 33) / -.2412930801 4594088413 8666666666 6 d-31     /
c
      data ntk1, ntak1, ntak12, xmin, xsml / 3*0, 2*0.d0 /
c
      if (ntk1.ne.0) go to 10
      eta = 0.1*sngl(d1mach(3))
      ntk1 = initds (bk1cs, 16, eta)
      ntak1 = initds (ak1cs, 38, eta)
      ntak12 = initds (ak12cs, 33, eta)
c
      xmin = dexp (dmax1(dlog(d1mach(1)), -dlog(d1mach(2))) + 0.01d0)
      xsml = dsqrt (4.0d0*d1mach(3))
c
 10   if (x.le.0.d0) call seteru (29hdbsk1e  x is zero or negative, 29,
     1  1, 2)
      if (x.gt.2.0d0) go to 20
c
      if (x.lt.xmin) call seteru (31hdbsk1e  x so small k1 overflows,
     1  31, 2, 2)
      y = 0.d0
      if (x.gt.xsml) y = x*x
      dbsk1e = dexp(x)*(dlog(0.5d0*x)*dbesi1(x) + (0.75d0 +
     1  dcsevl (0.5d0*y-1.d0, bk1cs, ntk1))/x )
      return
c
 20   if (x.le.8.d0) dbsk1e = (1.25d0 + dcsevl ((16.d0/x-5.d0)/3.d0,
     1  ak1cs, ntak1))/dsqrt(x)
      if (x.gt.8.d0) dbsk1e = (1.25d0 +
     1  dcsevl (16.d0/x-1.d0, ak12cs, ntak12))/dsqrt(x)
c
      return
      end
      subroutine seterr (messg, nmessg, nerr, iopt)
c
c  this version modified by w. fullerton to dump if iopt = 1 and
c  not recovering.
c  seterr sets lerror = nerr, optionally prints the message and dumps
c  according to the following rules...
c
c    if iopt = 1 and recovering      - just remember the error.
c    if iopt = 1 and not recovering  - print, dump and stop.
c    if iopt = 2                     - print, dump and stop.
c
c  input
c
c    messg  - the error message.
c    nmessg - the length of the message, in characters.
c    nerr   - the error number. must have nerr non-zero.
c    iopt   - the option. must have iopt=1 or 2.
c
c  error states -
c
c    1 - message length not positive.
c    2 - cannot have nerr=0.
c    3 - an unrecovered error followed by another error.
c    4 - bad value for iopt.
c
c  only the first 72 characters of the message are printed.
c
c  the error handler calls a subroutine named fdump to produce a
c  symbolic dump. to complete the package, a dummy version of fdump
c  is supplied, but it should be replaced by a locally written version
c  which at least gives a trace-back.
c
      integer messg(1)
c     external i1mach, i8save
c
c  the unit for error messages.
c
      iwunit=i1mach(4)
c
      if (nmessg.ge.1) go to 10
c
c  a message of non-positive length is fatal.
c
        write(iwunit,9000)
 9000   format(52h1error    1 in seterr - message length not positive.)
        go to 60
c
c  nw is the number of words the message occupies.
c
 10   nw=(min0(nmessg,72)-1)/i1mach(6)+1
c
      if (nerr.ne.0) go to 20
c
c  cannot turn the error state off using seterr.
c
        write(iwunit,9001)
 9001   format(42h1error    2 in seterr - cannot have nerr=0//
     1         34h the current error message follows///)
        call e9rint(messg,nw,nerr,.true.)
        itemp=i8save(1,1,.true.)
        go to 50
c
c  set lerror and test for a previous unrecovered error.
c
 20   if (i8save(1,nerr,.true.).eq.0) go to 30
c
        write(iwunit,9002)
 9002   format(23h1error    3 in seterr -,
     1         48h an unrecovered error followed by another error.//
     2         48h the previous and current error messages follow.///)
        call eprint
        call e9rint(messg,nw,nerr,.true.)
        go to 50
c
c  save this message in case it is not recovered from properly.
c
 30   call e9rint(messg,nw,nerr,.true.)
c
      if (iopt.eq.1 .or. iopt.eq.2) go to 40
c
c  must have iopt = 1 or 2.
c
        write(iwunit,9003)
 9003   format(42h1error    4 in seterr - bad value for iopt//
     1         34h the current error message follows///)
        go to 50
c
c  test for recovery.
c
 40   if (iopt.eq.2) go to 50
c
      if (i8save(2,0,.false.).eq.1) return
c
c     call eprint
c     stop
c
 50   call eprint
 60   call fdump
      stop
c
      end
      subroutine e9rint(messg,nw,nerr,save)
c
c  this routine stores the current error message or prints the old one,
c  if any, depending on whether or not save = .true. .
c
      integer messg(nw)
      logical save
c     external i1mach, i8save
c
c  messgp stores at least the first 72 characters of the previous
c  message. its length is machine dependent and must be at least
c
c       1 + 71/(the number of characters stored per integer word).
c
      integer messgp(36),fmt(14),ccplus
c
c  start with no previous message.
c
      data messgp(1)/1h1/, nwp/0/, nerrp/0/
c
c  set up the format for printing the error message.
c  the format is simply (a1,14x,72axx) where xx=i1mach(6) is the
c  number of characters stored per integer word.
c
      data ccplus  / 1h+ /
c
      data fmt( 1) / 1h( /
      data fmt( 2) / 1ha /
      data fmt( 3) / 1h1 /
      data fmt( 4) / 1h, /
      data fmt( 5) / 1h1 /
      data fmt( 6) / 1h4 /
      data fmt( 7) / 1hx /
      data fmt( 8) / 1h, /
      data fmt( 9) / 1h7 /
      data fmt(10) / 1h2 /
      data fmt(11) / 1ha /
      data fmt(12) / 1hx /
      data fmt(13) / 1hx /
      data fmt(14) / 1h) /
c
      if (.not.save) go to 20
c
c  save the message.
c
        nwp=nw
        nerrp=nerr
        do 10 i=1,nw
 10     messgp(i)=messg(i)
c
        go to 30
c
 20   if (i8save(1,0,.false.).eq.0) go to 30
c
c  print the message.
c
        iwunit=i1mach(4)
        write(iwunit,9000) nerrp
 9000   format(7h error ,i4,4h in )
c
        call s88fmt(2,i1mach(6),fmt(12))
        write(iwunit,fmt) ccplus,(messgp(i),i=1,nwp)
c
 30   return
c
      end
      integer function i8save(isw,ivalue,set)
c
c  if (isw = 1) i8save returns the current error number and
c               sets it to ivalue if set = .true. .
c
c  if (isw = 2) i8save returns the current recovery switch and
c               sets it to ivalue if set = .true. .
c
      logical set
c
      integer iparam(2)
c  iparam(1) is the error number and iparam(2) is the recovery switch.
c
c  start execution error free and with recovery turned off.
c
      data iparam(1) /0/,  iparam(2) /2/
c
      i8save=iparam(isw)
      if (set) iparam(isw)=ivalue
c
      return
c
      end
      subroutine eprint
c
c  this subroutine prints the last error message, if any.
c
      integer messg(1)
c
      call e9rint(messg,1,1,.false.)
      return
c
      end
      subroutine s88fmt( n, w, ifmt )
c
c  s88fmt  replaces ifmt(1), ... , ifmt(n) with
c  the characters corresponding to the n least significant
c  digits of w.
c
      integer n,w,ifmt(n)
c
      integer nt,wt,digits(10)
c
      data digits( 1) / 1h0 /
      data digits( 2) / 1h1 /
      data digits( 3) / 1h2 /
      data digits( 4) / 1h3 /
      data digits( 5) / 1h4 /
      data digits( 6) / 1h5 /
      data digits( 7) / 1h6 /
      data digits( 8) / 1h7 /
      data digits( 9) / 1h8 /
      data digits(10) / 1h9 /
c
      nt = n
      wt = w
c
 10   if (nt .le. 0) return
        idigit = mod( wt, 10 )
        ifmt(nt) = digits(idigit+1)
        wt = wt/10
        nt = nt - 1
        go to 10
c
      end
      subroutine fdump
      return
      end

C
C  ---------------------------------------
C  Compare scalar fnlib routines to vfnlib
C  ---------------------------------------
C
C  Ron Boisvert and Bonita Saunders
C  Computing and Applied Mathematics Laboratory
C  National Institute of Standards and Technology
C  Gaithersburg, MD 20899
C
C  boisvert@cam.nist.gov
C  saunders@cam.nist.gov
C
C  21 Feb 91
C
C----------------------------------------------------------------------
C
      INTEGER M
      PARAMETER ( M=2000 )
C
      DOUBLE PRECISION   XLO, XHI, X, F, FX, WORK
      INTEGER   I, IWORK, J
C
      DIMENSION XLO(3,8), XHI(3,8), X(M), F(M), FX(M), WORK(7*M),
     +          IWORK(M)
C
      SAVE XLO, XHI
C
      EXTERNAL  DVI0, DBESI0, DVI1, DBESI1, DVJ0, DBESJ0, DVJ1, DBESJ1,
     +          DVK0, DBESK0, DVK1, DBESK1, DVY0, DBESY0, DVY1, DBESY1
C
      DATA ((XLO(I,J),XHI(I,J),I=1,3),J=1,8) /
     +     0.0D0, 3.0D0 , -8.0D0, -3.0D0, 10.0D0, 50.0D0,
     +     0.0D0, 3.0D0 , -8.0D0, -3.0D0, 10.0D0, 50.0D0,
     +     0.0D0, 4.0D0 , -8.0D0, -4.0D0,  8.0D0, 50.0D0,
     +     0.0D0, 4.0D0 , -8.0D0, -4.0D0,  8.0D0, 50.0D0,
     +     0.0D0, 2.0D0 ,  2.0D0,  8.0D0, 10.0D0, 50.0D0,
     +     0.0D0, 2.0D0 ,  2.0D0,  8.0D0, 10.0D0, 50.0D0,
     +     0.0D0, 1.0D-3,  1.0D0,  4.0D0,  4.0D0, 50.0D0,
     +     0.0D0, 1.0D-3,  1.0D0,  4.0D0,  4.0D0, 50.0D0  /
C
C----------------------------------------------------------------------
C
C     ... TEST EACH USER-CALLABLE ROUTINE
C
      CALL TEST('I0',DBESI0,DVI0,XLO(1,1),XHI(1,1),M,X,F,FX,WORK,IWORK)
      CALL TEST('I1',DBESI1,DVI1,XLO(1,2),XHI(1,2),M,X,F,FX,WORK,IWORK)
      CALL TEST('J0',DBESJ0,DVJ0,XLO(1,3),XHI(1,3),M,X,F,FX,WORK,IWORK)
      CALL TEST('J1',DBESJ1,DVJ1,XLO(1,4),XHI(1,4),M,X,F,FX,WORK,IWORK)
      CALL TEST('K0',DBESK0,DVK0,XLO(1,5),XHI(1,5),M,X,F,FX,WORK,IWORK)
      CALL TEST('K1',DBESK1,DVK1,XLO(1,6),XHI(1,6),M,X,F,FX,WORK,IWORK)
      CALL TEST('Y0',DBESY0,DVY0,XLO(1,7),XHI(1,7),M,X,F,FX,WORK,IWORK)
      CALL TEST('Y1',DBESY1,DVY1,XLO(1,8),XHI(1,8),M,X,F,FX,WORK,IWORK)
C
      END

      SUBROUTINE TEST (NAME, SFUN, VSUB, XLOW, XHIGH, M, X, F, FX, 
     +                 WORK, IWORK)
C
C----------------------------------------------------------------------
C
C  compare scalar function sfun to vectorized subroutine vsub for m
C  arguments distributed in various ways in the intervals
C  (xlo(i),xhi(i)), i=1,2,3.
C
C----------------------------------------------------------------------
C
C  ----------
C  PARAMETERS
C  ----------
C
      DOUBLE PRECISION SFUN, XLOW, XHIGH, X, F, FX, WORK
      INTEGER IWORK, M
      CHARACTER*(*) NAME
C
      DIMENSION XLOW(3), XHIGH(3), X(M), F(M), FX(M), WORK(6*M),
     +          IWORK(M)
C
      EXTERNAL SFUN, VSUB
C
C  ---------------
C  LOCAL VARIABLES
C  ---------------
C
      LOGICAL DEBUG
      INTEGER NREP, NTRY
      PARAMETER ( DEBUG=.FALSE., NREP=15, NTRY=19 )
C
      INTEGER I, INFO, IREP, IOFF, IPC(NTRY), JPC(NTRY), K, KPC,
     +        MA, MB, MC, NFAIL, NREPV
      DOUBLE PRECISION DIFF, DMAX, DX, EPS, D1MACH
      REAL SECOND, SU, SUMAX, SUMIN, SUAVE, T0, T1, TS, TSAVE, TV,
     +        TVAVE
C
      SAVE IPC, JPC
C
      DATA IPC /  0,  0,  0,  0,  0, 25, 25, 25, 25, 33, 50, 50, 50,
     +           75, 75,100,  1,  1, 98/
      DATA JPC /  0, 25, 50, 75,100,  0, 25, 50, 75, 33,  0, 25, 50,
     +            0, 25,  0,  1, 98,  1/
C
C----------------------------------------------------------------------
C
C  ntry = number of trials (each with different argument distribution)
C  nrep = number of times test repeated to get reliable timings
C
C  the distribution of arguments in trial k is
C
C     ipc(k)%            in (xlow(1),xhigh(1))
C     jpc(k)%            in (xlow(2),xhigh(2))
C     100-ipc(k)-jpc(k)% in (xlow(3),xhigh(3))
C
C  scalar results (fx) and vector results (f) are compared and accepted
C  if   abs((f - fx)/fx) < eps  for fx > 1  (relative difference)
C  or   abs(f - fx)      < eps  for fx < 1  (absolute difference)
C
C----------------------------------------------------------------------
C
      EPS = 5.0D0*D1MACH(4)

C     ... TRIGGER INITIALIZATIONS
C
      X(1) = 1.0D0
      CALL VSUB(1,X,F,WORK,IWORK,INFO)
      FX(1) = SFUN(X(1))
C
      PRINT *
      PRINT *
      PRINT *,' ----------'
      PRINT *,' TEST OF ',NAME
      PRINT *,' ----------'
      PRINT *
      PRINT *,' Vector Length = ',M
      PRINT *,' Repetitions   = ',NREP
      PRINT *,' Eps           = ',EPS
      PRINT *
      PRINT *,' Range A = (',XLOW(1),',',XHIGH(1),')'
      PRINT *,' Range B = (',XLOW(2),',',XHIGH(2),')'
      PRINT *,' Range C = (',XLOW(3),',',XHIGH(3),')'
      PRINT *
      PRINT *,' %A      = percent arguments in range A'
      PRINT *,' %B      = percent arguments in range B'
      PRINT *,' %C      = percent arguments in range C'
      PRINT *,' Scalar  = time for scalar code (FNLIB)'
      PRINT *,' Vector  = time for vector code (VFNLIB)'
      PRINT *,' SU      = speedup (Scalar/Vector)'
      PRINT *,' Nfail   = number of differences .gt. eps'
      PRINT *,' Dmax    = maximum difference'
      PRINT *
      PRINT *,' Times given in seconds/argument.'
      PRINT *
      PRINT *,' Differences are abs(s-v)/max(abs(s),1))'
      PRINT *,' where s = FNLIB result and v = VFNLIB result.'
      PRINT *
      PRINT *,'  %A  %B  %C    Scalar    Vector   SU   Nfail   Dmax'
      PRINT *,' -------------------------------------------------------'
C
      SUMAX = 0.0
      SUMIN = 1000.0
      SUAVE = 0.0
      TSAVE = 0.0
      TVAVE = 0.0
C
C     ... PERFORM EACH TEST
C
      DO 100 K=1,NTRY
C
         MA = IPC(K)*M/100
         MB = JPC(K)*M/100
         KPC = 100 - IPC(K) - JPC(K)
         MC = M - MA - MB
C
C        ... SELECT X VECTOR
C
         IOFF = 0
C
         IF (MA .GT. 0) THEN
            DX = (XHIGH(1) - XLOW(1))/MA
            DO 10 I=1,MA
               X(I) = XLOW(1) + I*DX
   10       CONTINUE
            IOFF = IOFF + MA
         ENDIF
C
         IF (MB .GT. 0) THEN
            DX = (XHIGH(2) - XLOW(2))/MB
            DO 20 I=1,MB
               X(I+IOFF) = XLOW(2) + I*DX
   20       CONTINUE
            IOFF = IOFF + MB
         ENDIF
C
         IF (MC .GT. 0) THEN
            DX = (XHIGH(3) - XLOW(3))/MC
            DO 30 I=1,MC
               X(I+IOFF) = XLOW(3) + I*DX
   30       CONTINUE
         ENDIF
C
C        ... RUN VECTOR CODE
C
	 NREPV = 5*NREP
         T0 = SECOND()
         DO 40 IREP=1,NREPV
            CALL VSUB(M,X,F,WORK,IWORK,INFO)
   40    CONTINUE
         T1 = SECOND()
         TV = (T1 - T0)/NREPV
         IF (INFO .LT. 0) PRINT *,'Warning in ',NAME,' : Info = ',INFO
         IF (INFO .GT. 0) PRINT *,'Failure in ',NAME,' : Info = ',INFO,
     +                            '  Position = ',IWORK(1)
C
C        ... RUN SCALAR CODE
C
         T0 = SECOND()
         DO 60 IREP=1,NREP
            DO 50 I=1,M
               FX(I) = SFUN(X(I))
   50       CONTINUE
   60    CONTINUE
         T1 = SECOND()
         TS = (T1 - T0)/NREP
C
C        ... COMPUTE TIMINGS
C
         SU = TS/TV
         SUMIN = MIN(SU,SUMIN)
         SUMAX = MAX(SU,SUMAX)
         SUAVE = SUAVE + SU
         TSAVE = TSAVE + TS
         TVAVE = TVAVE + TV
C
C        ... EVALUATE RESULTS
C
         DMAX = 0.0D0
         NFAIL = 0
         DO 70 I=1,M
            DIFF = ABS(FX(I)-F(I))/MAX(ABS(FX(I)),1.0D0)
            DMAX = MAX(DIFF,DMAX)
            IF (DIFF .GT. EPS) THEN
               NFAIL = NFAIL + 1
               IF (DEBUG) PRINT *,'Failure ',I,X(I),F(I),FX(I),DIFF
            ENDIF
  70     CONTINUE
C
C         ... PRINT RESULTS
C
          PRINT 1000, IPC(K),JPC(K),KPC,TS/M,TV/M,SU,NFAIL,DMAX
C
  100 CONTINUE
C
      SUAVE = SUAVE/NTRY
      TSAVE = TSAVE/M/NTRY
      TVAVE = TVAVE/M/NTRY
C
      PRINT *
      PRINT 1002,'Ave Scalar   = ',TSAVE,' sec/argument'
      PRINT 1002,'Ave Vector   = ',TVAVE,' sec/argument'
      PRINT *
      PRINT 1001,'Max speed up = ',SUMAX
      PRINT 1001,'Min speed up = ',SUMIN
      PRINT *
      PRINT 1001,'Ave speed up = ',SUAVE
C
      RETURN
 1000 FORMAT(3(1X,I3),1P,2(3X,E7.1),0P,1X,F5.1,1X,I4,3X,1P,E7.1)
 1001 FORMAT(1X,A,F5.1,A)
 1002 FORMAT(1X,A,1P,E7.1,A)
      END
      SUBROUTINE DVI0 (M, X, F, WORK, IWORK, INFO)
C***BEGIN PROLOGUE  DVI0
C***PURPOSE  Computes the hyperbolic Bessel function of the first kind
C            of order zero (I0) for a vector of real arguments
C***LIBRARY   VFNLIB
C***CATEGORY  C10B1
C***TYPE      DOUBLE PRECISION (VI0-S, DVI0-D)
C***KEYWORDS  BESSEL FUNCTION,FIRST KIND,HYPERBOLIC BESSEL FUNCTION,
C             MODIFIED BESSEL FUNCTION, ORDER ZERO, VECTORIZED
C***AUTHOR  SAUNDERS, B. V., (NIST)
C             COMPUTING AND APPLIED MATHEMATICS LABORATORY
C             NATIONAL INSTITUTE OF STANDARDS AND TECHNOLOGY
C             GAITHERSBURG, MD 20899
C           BOISVERT, R. F., (NIST)
C             COMPUTING AND APPLIED MATHEMATICS LABORATORY
C             NATIONAL INSTITUTE OF STANDARDS AND TECHNOLOGY
C             GAITHERSBURG, MD 20899
C***DESCRIPTION
C
C   DVI0 computes the modified (hyperbolic) Bessel function of the
C   first kind of order zero (I0) for real arguments using uniform
C   approximation by Chebyshev polynomials.
C
C
C   P A R A M E T E R S
C
C   M      (Input) Integer (M .GT. 0)
C          The number of arguments at which the function is to be
C          evaluated.
C
C   X      (Input) Double precision array of length M
C          The arguments at which the function is to be evaluated are
C          stored in X(1) to X(M) in any order.
C
C   F      (Output) Double precision array of length M
C          F(i) contains the value of the function at X(i), i=1,..,M.
C
C   WORK   (Work) Double precision vector of length 7*M
C
C   IWORK  (Work) Integer vector of length M
C
C   INFO   (Output) Integer
C          Indicates status of computed result. The following table
C          lists possible values and their meanings.  If OK=Yes then
C          all F(i) have been set, otherwise none have been set.
C
C          INFO  OK            Description
C          ------------------------------------------------------------
C            0   Yes  Successfull execution.
C            1   No   Error: M .LE. 0
C            2   No   Error: Some abs(X(i)) so big I0 overflows.
C                     The index of the first offending argument is
C                     returned in IWORK(1).
C
C
C *********************************************************************
C   This routine is a modification of the function DBESI0 developed by
C   W. Fullerton of LANL.
C *********************************************************************
C
C***REFERENCES  (NONE)
C***ROUTINES CALLED  DWI0
C***REVISION HISTORY  (YYMMDD)
C   910226  DATE WRITTEN
C***END PROLOGUE  DVI0
C
C     ----------
C     PARAMETERS
C     ----------
C
      INTEGER INFO, IWORK, M
      DOUBLE PRECISION F, X, WORK
C
      DIMENSION X(M), F(M), WORK(7*M), IWORK(M)
C
C     ---------------
C     LOCAL VARIABLES
C     ---------------
C
      INTEGER IWB0, IWB1, IWB2, IWY, IWTC, IWYC, IWZC, JIN
C
C***FIRST EXECUTABLE STATEMENT  DVI0
C
C     ... PARTITION WORK ARRAYS
C
      IWY   = 1
      IWTC  = IWY  + M
      IWYC  = IWTC + M
      IWZC  = IWYC + M
      IWB0  = IWZC + M
      IWB1  = IWB0 + M
      IWB2  = IWB1 + M
C     Total = IWB2 + M
C
      JIN   = 1
C     Total = JIN  + M
C
C     ... DWI0 DOES ALL THE WORK
C
      CALL DWI0(M,X,F,WORK(IWY),WORK(IWTC),WORK(IWYC),WORK(IWZC),
     +         IWORK(JIN),WORK(IWB0),WORK(IWB1),WORK(IWB2),INFO)
C
      RETURN
      END
      SUBROUTINE DWI0 (M, X, F, Y, TCMP, YCMP, ZCMP, INDX, B0, B1, B2, 
     +   INFO)
C***BEGIN PROLOGUE  DWI0
C***SUBSIDIARY
C***PURPOSE  Computes the hyperbolic Bessel function of the first kind
C            of order zero (I0) for a vector of arguments
C***LIBRARY   VFNLIB
C***AUTHOR  SAUNDERS, B. V., (NIST)
C             COMPUTING AND APPLIED MATHEMATICS LABORATORY
C             NATIONAL INSTITUTE OF STANDARDS AND TECHNOLOGY
C             GAITHERSBURG, MD 20899
C           BOISVERT, R. F., (NIST)
C             COMPUTING AND APPLIED MATHEMATICS LABORATORY
C             NATIONAL INSTITUTE OF STANDARDS AND TECHNOLOGY
C             GAITHERSBURG, MD 20899
C***ROUTINES CALLED  IDWCS, D1MACH, DWNGT, DWGTHR, DWGTLE, DWGT, DWLE,
C                    DWSCTR, DWCS
C***REVISION HISTORY  (YYMMDD)
C   910226  DATE WRITTEN
C***END PROLOGUE  DWI0
C
C  ----------
C  PARAMETERS
C  ----------
C
      INTEGER INFO, INDX, M
      DOUBLE PRECISION B0, B1, B2, F, X, Y, TCMP, YCMP, ZCMP
C
      DIMENSION B0(M), B1(M), B2(M), F(M), INDX(M), X(M), Y(M),
     +          TCMP(M), YCMP(M), ZCMP(M)
C
C  ---------------
C  LOCAL VARIABLES
C  ---------------
C
      INTEGER LAI0, LAI02, LBI0
      PARAMETER ( LAI0=46, LAI02=69, LBI0=18 )
C
      INTEGER I, IDWCS, J, KEY, N, NTI0, NTAI0, NTAI02
      DOUBLE PRECISION AI0CS, AI02CS, BI0CS, EPMACH, EPS, D1MACH, XSML,
     +        XMAX
C
      DIMENSION AI0CS(LAI0), AI02CS(LAI02), BI0CS(LBI0)
C
      SAVE AI0CS, AI02CS, BI0CS, N, NTAI0, NTAI02, NTI0, XSML, XMAX
C
C----------------------------------------------------------------------
C
C Series for BI0        on the interval  0.          to  9.00000D+00
C                                        with weighted error   9.51E-34
C                                         log weighted error  33.02
C                               significant figures required  33.31
C                                    decimal places required  33.65
C
      DATA BI0 CS(  1) / -.7660547252 8391449510 8189497624 3285 D-1   /
      DATA BI0 CS(  2) / +.1927337953 9938082699 5240875088 1196 D+1   /
      DATA BI0 CS(  3) / +.2282644586 9203013389 3702929233 0415 D+0   /
      DATA BI0 CS(  4) / +.1304891466 7072904280 7933421069 1888 D-1   /
      DATA BI0 CS(  5) / +.4344270900 8164874513 7868268102 6107 D-3   /
      DATA BI0 CS(  6) / +.9422657686 0019346639 2317174411 8766 D-5   /
      DATA BI0 CS(  7) / +.1434006289 5106910799 6209187817 9957 D-6   /
      DATA BI0 CS(  8) / +.1613849069 6617490699 1541971999 4611 D-8   /
      DATA BI0 CS(  9) / +.1396650044 5356696994 9509270814 2522 D-10  /
      DATA BI0 CS( 10) / +.9579451725 5054453446 2752317189 3333 D-13  /
      DATA BI0 CS( 11) / +.5333981859 8625021310 1510774400 0000 D-15  /
      DATA BI0 CS( 12) / +.2458716088 4374707746 9678591999 9999 D-17  /
      DATA BI0 CS( 13) / +.9535680890 2487700269 4434133333 3333 D-20  /
      DATA BI0 CS( 14) / +.3154382039 7214273367 8933333333 3333 D-22  /
      DATA BI0 CS( 15) / +.9004564101 0946374314 6666666666 6666 D-25  /
      DATA BI0 CS( 16) / +.2240647369 1236700160 0000000000 0000 D-27  /
      DATA BI0 CS( 17) / +.4903034603 2428373333 3333333333 3333 D-30  /
      DATA BI0 CS( 18) / +.9508172606 1226666666 6666666666 6666 D-33  /

C
C-------------------------------------------------------------------
C
C Series for AI0        on the interval  1.25000D-01 to  3.33333D-01
C                                        with weighted error   2.74E-32
C                                         log weighted error  31.56
C                               significant figures required  30.15
C                                    decimal places required  32.39
C
      DATA AI0 CS(  1) / +.7575994494 0237959427 2987203743 8 D-1      /
      DATA AI0 CS(  2) / +.7591380810 8233455072 9297873320 4 D-2      /
      DATA AI0 CS(  3) / +.4153131338 9237505018 6319749138 2 D-3      /
      DATA AI0 CS(  4) / +.1070076463 4390730735 8242970217 0 D-4      /
      DATA AI0 CS(  5) / -.7901179979 2128946607 5031948573 0 D-5      /
      DATA AI0 CS(  6) / -.7826143501 4387522697 8898980690 9 D-6      /
      DATA AI0 CS(  7) / +.2783849942 9488708063 8118538985 7 D-6      /
      DATA AI0 CS(  8) / +.8252472600 6120271919 6682913319 8 D-8      /
      DATA AI0 CS(  9) / -.1204463945 5201991790 5496089110 3 D-7      /
      DATA AI0 CS( 10) / +.1559648598 5060764436 1228752792 8 D-8      /
      DATA AI0 CS( 11) / +.2292556367 1033165434 7725480285 7 D-9      /
      DATA AI0 CS( 12) / -.1191622884 2790646036 7777423447 8 D-9      /
      DATA AI0 CS( 13) / +.1757854916 0324098302 1833124774 3 D-10     /
      DATA AI0 CS( 14) / +.1128224463 2189005171 4441135682 4 D-11     /
      DATA AI0 CS( 15) / -.1146848625 9272988777 2963387698 2 D-11     /
      DATA AI0 CS( 16) / +.2715592054 8036628726 4365192160 6 D-12     /
      DATA AI0 CS( 17) / -.2415874666 5626878384 4247572028 1 D-13     /
      DATA AI0 CS( 18) / -.6084469888 2551250646 0609963922 4 D-14     /
      DATA AI0 CS( 19) / +.3145705077 1754772937 0836026730 3 D-14     /
      DATA AI0 CS( 20) / -.7172212924 8711877179 6217505917 6 D-15     /
      DATA AI0 CS( 21) / +.7874493403 4541033960 8390960332 7 D-16     /
      DATA AI0 CS( 22) / +.1004802753 0094624023 4524457183 9 D-16     /
      DATA AI0 CS( 23) / -.7566895365 3505348534 2843588881 0 D-17     /
      DATA AI0 CS( 24) / +.2150380106 8761198878 1205128784 5 D-17     /
      DATA AI0 CS( 25) / -.3754858341 8308744291 5158445260 8 D-18     /
      DATA AI0 CS( 26) / +.2354065842 2269925769 0075710532 2 D-19     /
      DATA AI0 CS( 27) / +.1114667612 0479285302 2637335511 0 D-19     /
      DATA AI0 CS( 28) / -.5398891884 3969903786 9677932270 9 D-20     /
      DATA AI0 CS( 29) / +.1439598792 2407526770 4285840452 2 D-20     /
      DATA AI0 CS( 30) / -.2591916360 1110934064 6081840196 2 D-21     /
      DATA AI0 CS( 31) / +.2238133183 9985839074 3409229824 0 D-22     /
      DATA AI0 CS( 32) / +.5250672575 3647711727 7221683199 9 D-23     /
      DATA AI0 CS( 33) / -.3249904138 5332307841 7343228586 6 D-23     /
      DATA AI0 CS( 34) / +.9924214103 2050379278 5728471040 0 D-24     /
      DATA AI0 CS( 35) / -.2164992254 2446695231 4655429973 3 D-24     /
      DATA AI0 CS( 36) / +.3233609471 9435940839 7333299199 9 D-25     /
      DATA AI0 CS( 37) / -.1184620207 3967424898 2473386666 6 D-26     /
      DATA AI0 CS( 38) / -.1281671853 9504986505 4833868799 9 D-26     /
      DATA AI0 CS( 39) / +.5827015182 2793905116 0556885333 3 D-27     /
      DATA AI0 CS( 40) / -.1668222326 0261097193 6450150399 9 D-27     /
      DATA AI0 CS( 41) / +.3625309510 5415699757 0068480000 0 D-28     /
      DATA AI0 CS( 42) / -.5733627999 0557135899 4595839999 9 D-29     /
      DATA AI0 CS( 43) / +.3736796722 0630982296 4258133333 3 D-30     /
      DATA AI0 CS( 44) / +.1602073983 1568519633 6551253333 3 D-30     /
      DATA AI0 CS( 45) / -.8700424864 0572298845 2249599999 9 D-31     /
      DATA AI0 CS( 46) / +.2741320937 9374811456 0341333333 3 D-31     /
C
C-------------------------------------------------------------------
C
C Series for AI02       on the interval  0.          to  1.25000D-01
C                                        with weighted error   1.97E-32
C                                         log weighted error  31.71
C                               significant figures required  30.15
C                                    decimal places required  32.63
C
      DATA AI02CS(  1) / +.5449041101 4108831607 8960962268 0 D-1      /
      DATA AI02CS(  2) / +.3369116478 2556940898 9785662979 9 D-2      /
      DATA AI02CS(  3) / +.6889758346 9168239842 6263914301 1 D-4      /
      DATA AI02CS(  4) / +.2891370520 8347564829 6692402323 2 D-5      /
      DATA AI02CS(  5) / +.2048918589 4690637418 2760534093 1 D-6      /
      DATA AI02CS(  6) / +.2266668990 4981780645 9327743136 1 D-7      /
      DATA AI02CS(  7) / +.3396232025 7083863451 5084396952 3 D-8      /
      DATA AI02CS(  8) / +.4940602388 2249695891 0482449783 5 D-9      /
      DATA AI02CS(  9) / +.1188914710 7846438342 4084525196 3 D-10     /
      DATA AI02CS( 10) / -.3149916527 9632413645 3864862961 9 D-10     /
      DATA AI02CS( 11) / -.1321581184 0447713118 7540739926 7 D-10     /
      DATA AI02CS( 12) / -.1794178531 5068061177 7943574026 9 D-11     /
      DATA AI02CS( 13) / +.7180124451 3836662336 7106429346 9 D-12     /
      DATA AI02CS( 14) / +.3852778382 7421427011 4089801777 6 D-12     /
      DATA AI02CS( 15) / +.1540086217 5214098269 1325823339 7 D-13     /
      DATA AI02CS( 16) / -.4150569347 2872220866 2689972015 6 D-13     /
      DATA AI02CS( 17) / -.9554846698 8283076487 0214494312 5 D-14     /
      DATA AI02CS( 18) / +.3811680669 3526224207 4605535511 8 D-14     /
      DATA AI02CS( 19) / +.1772560133 0565263836 0493266675 8 D-14     /
      DATA AI02CS( 20) / -.3425485619 6772191346 1924790328 2 D-15     /
      DATA AI02CS( 21) / -.2827623980 5165834849 4205593759 4 D-15     /
      DATA AI02CS( 22) / +.3461222867 6974610930 9706250813 4 D-16     /
      DATA AI02CS( 23) / +.4465621420 2967599990 1042054284 3 D-16     /
      DATA AI02CS( 24) / -.4830504485 9441820712 5525403795 4 D-17     /
      DATA AI02CS( 25) / -.7233180487 8747539545 6227240924 5 D-17     /
      DATA AI02CS( 26) / +.9921475412 1736985988 8046093981 0 D-18     /
      DATA AI02CS( 27) / +.1193650890 8459820855 0439949924 2 D-17     /
      DATA AI02CS( 28) / -.2488709837 1508072357 2054491660 2 D-18     /
      DATA AI02CS( 29) / -.1938426454 1609059289 8469781132 6 D-18     /
      DATA AI02CS( 30) / +.6444656697 3734438687 8301949394 9 D-19     /
      DATA AI02CS( 31) / +.2886051596 2892243264 8171383073 4 D-19     /
      DATA AI02CS( 32) / -.1601954907 1749718070 6167156200 7 D-19     /
      DATA AI02CS( 33) / -.3270815010 5923147208 9193567485 9 D-20     /
      DATA AI02CS( 34) / +.3686932283 8264091811 4600723939 3 D-20     /
      DATA AI02CS( 35) / +.1268297648 0309501530 1359529710 9 D-22     /
      DATA AI02CS( 36) / -.7549825019 3772739076 9636664410 1 D-21     /
      DATA AI02CS( 37) / +.1502133571 3778353496 3712789053 4 D-21     /
      DATA AI02CS( 38) / +.1265195883 5096485349 3208799248 3 D-21     /
      DATA AI02CS( 39) / -.6100998370 0836807086 2940891600 2 D-22     /
      DATA AI02CS( 40) / -.1268809629 2601282643 6872095924 2 D-22     /
      DATA AI02CS( 41) / +.1661016099 8907414578 4038487490 5 D-22     /
      DATA AI02CS( 42) / -.1585194335 7658855793 7970504881 4 D-23     /
      DATA AI02CS( 43) / -.3302645405 9682178009 5381766755 6 D-23     /
      DATA AI02CS( 44) / +.1313580902 8392397817 4039623117 4 D-23     /
      DATA AI02CS( 45) / +.3689040246 6711567933 1425637280 4 D-24     /
      DATA AI02CS( 46) / -.4210141910 4616891492 1978247249 9 D-24     /
      DATA AI02CS( 47) / +.4791954591 0828657806 3171401373 0 D-25     /
      DATA AI02CS( 48) / +.8459470390 2218217952 9971707412 4 D-25     /
      DATA AI02CS( 49) / -.4039800940 8728324931 4607937181 0 D-25     /
      DATA AI02CS( 50) / -.6434714653 6504313473 0100850469 5 D-26     /
      DATA AI02CS( 51) / +.1225743398 8756659903 4464736990 5 D-25     /
      DATA AI02CS( 52) / -.2934391316 0257089231 9879821175 4 D-26     /
      DATA AI02CS( 53) / -.1961311309 1949829262 0371205728 9 D-26     /
      DATA AI02CS( 54) / +.1503520374 8221934241 6229900309 8 D-26     /
      DATA AI02CS( 55) / -.9588720515 7448265520 3386388206 9 D-28     /
      DATA AI02CS( 56) / -.3483339380 8170454863 9441108511 4 D-27     /
      DATA AI02CS( 57) / +.1690903610 2630436730 6244960725 6 D-27     /
      DATA AI02CS( 58) / +.1982866538 7356030438 9400115718 8 D-28     /
      DATA AI02CS( 59) / -.5317498081 4918162145 7583002528 4 D-28     /
      DATA AI02CS( 60) / +.1803306629 8883929462 3501450390 1 D-28     /
      DATA AI02CS( 61) / +.6213093341 4548931758 8405311242 2 D-29     /
      DATA AI02CS( 62) / -.7692189292 7721618632 0072806673 0 D-29     /
      DATA AI02CS( 63) / +.1858252826 1117025426 2556016596 3 D-29     /
      DATA AI02CS( 64) / +.1237585142 2813957248 9927154554 1 D-29     /
      DATA AI02CS( 65) / -.1102259120 4092238032 1779478779 2 D-29     /
      DATA AI02CS( 66) / +.1886287118 0397044900 7787447943 1 D-30     /
      DATA AI02CS( 67) / +.2160196872 2436589131 4903141406 0 D-30     /
      DATA AI02CS( 68) / -.1605454124 9197432005 8446594965 5 D-30     /
      DATA AI02CS( 69) / +.1965352984 5942906039 3884807331 8 D-31     /
C
C-------------------------------------------------------------------
C
      DATA NTI0 / 0 /
C
C***FIRST EXECUTABLE STATEMENT  DWI0
C
      IF (M .LE. 0)  GO TO 910
C
      IF (NTI0.EQ.0) THEN
         EPMACH = D1MACH(3)
         EPS = 0.10D0*EPMACH
         NTI0 = IDWCS(BI0CS, LBI0, EPS)
         NTAI0 = IDWCS(AI0CS, LAI0, EPS)
         NTAI02 = IDWCS(AI02CS, LAI02, EPS)
         XSML = SQRT(4.0D0*EPMACH)
         XMAX = LOG(D1MACH(2))
      ENDIF
C
      DO 10 I=1,M
         Y(I) = ABS(X(I))
  10  CONTINUE
C
      CALL DWNGT(M,Y,XMAX,KEY)
      IF (KEY .NE. 0)  GO TO 920
C
C  ----------------
C  CASE Y .LE. XSML
C  ----------------
C
      DO 15 I=1,M
         F(I) = 1.0D0
  15  CONTINUE
C
C  --------------------------
C  CASE  XSML .LT. Y .LE. 3.0
C  --------------------------
C
      CALL DWGTLE(M,Y,XSML,3.0D0,N,INDX)
      IF (N .GT. 0) THEN
         CALL DWGTHR(N,Y,INDX,YCMP)
         DO 20 J=1,N
            TCMP(J) = YCMP(J)**2/4.50D0 - 1.0D0
  20     CONTINUE
         CALL DWCS(N,TCMP,BI0CS,NTI0,ZCMP,B0,B1,B2)
         DO 30 J=1,N
            ZCMP(J) = 2.750D0 + ZCMP(J)
  30     CONTINUE
         CALL DWSCTR(N,ZCMP,INDX,F)
      ENDIF
C
C  -------------------------
C  CASE  3.0 .LT. Y .LE. 8.0
C  -------------------------
C
      CALL DWGTLE(M,Y,3.0D0,8.0D0,N,INDX)
      IF (N .GT. 0) THEN
         CALL DWGTHR(N,Y,INDX,YCMP)
         DO 50 J=1,N
            TCMP(J) = (48.0D0/YCMP(J) - 11.0D0)/5.0D0
  50     CONTINUE
         CALL DWCS(N,TCMP,AI0CS,NTAI0,ZCMP,B0,B1,B2)
         DO 60 J=1,N
            ZCMP(J) = EXP(YCMP(J))*(0.3750D0+ZCMP(J))/SQRT(YCMP(J))
  60     CONTINUE
         CALL DWSCTR(N,ZCMP,INDX,F)
      ENDIF
C
C  ----------------
C  CASE  Y .GT. 8.0
C  ----------------
C
      CALL DWGT(M,Y,8.0D0,N,INDX)
      IF (N .GT. 0) THEN
         CALL DWGTHR(N,Y,INDX,YCMP)
        DO 80 J=1,N
            TCMP(J) = 16.0D0/YCMP(J) - 1.0D0
  80     CONTINUE
         CALL DWCS(N,TCMP,AI02CS,NTAI02,ZCMP,B0,B1,B2)
         DO 90 J=1,N
            ZCMP(J) = EXP(YCMP(J))*(0.3750D0+ZCMP(J))/SQRT(YCMP(J))
  90     CONTINUE
         CALL DWSCTR(N,ZCMP,INDX,F)
      ENDIF
C
C  -----
C  EXITS
C  -----
C
C     ... NORMAL
C
      INFO = 0
      GO TO 999
C
C     ... M .LE. 0
C
  910 CONTINUE
      INFO = 1
      GO TO 999
C
C     ... ABS(X) SO LARGE I0 OVERFLOWS
C
  920 CONTINUE
      INFO = 2
      INDX(1) = KEY
      GO TO 999
C
  999 CONTINUE
      RETURN
      END

       SUBROUTINE DVI1 (M, X, F, WORK, IWORK, INFO)
C***BEGIN PROLOGUE  DVI1
C***PURPOSE  Computes the hyperbolic Bessel function of the first kind
C            of order one (I1) for a vector of real arguments
C***LIBRARY   VFNLIB
C***CATEGORY  C10B1
C***TYPE      DOUBLE PRECISION (VI1-S, DVI1-D)
C***KEYWORDS  BESSEL FUNCTION,FIRST KIND,HYPERBOLIC BESSEL FUNCTION,
C             MODIFIED BESSEL FUNCTION, ORDER ONE, VECTORIZED
C***AUTHOR  SAUNDERS, B. V., (NIST)
C             COMPUTING AND APPLIED MATHEMATICS LABORATORY
C             NATIONAL INSTITUTE OF STANDARDS AND TECHNOLOGY
C             GAITHERSBURG, MD 20899
C           BOISVERT, R. F., (NIST)
C             COMPUTING AND APPLIED MATHEMATICS LABORATORY
C             NATIONAL INSTITUTE OF STANDARDS AND TECHNOLOGY
C             GAITHERSBURG, MD 20899
C***DESCRIPTION
C
C   DVI1 computes the modified (hyperbolic) Bessel function of the
C   first kind of order one (I1) for real arguments using uniform
C   approximation by Chebyshev polynomials.
C
C
C   P A R A M E T E R S
C
C   M      (Input) Integer
C          The number of arguments at which the function is to be
C          evaluated.
C
C   X      (Input) Double precision array of length M
C          The arguments at which the function is to be evaluated are
C          stored in X(1) to X(M) in any order.
C
C   F      (Output) Double precision array of length M
C          F(i) contains the value of the function at X(i), i=1,..,M.
C
C   WORK   (Work) Double precision vector of length 7*M
C
C   IWORK  (Work) Integer vector of length M
C
C   INFO   (Output) Integer
C          Indicates status of computed result. The following table
C          lists possible values and their meanings.  If OK=Yes then
C          all F(i) have been set, otherwise none have been set.
C
C          INFO  OK            Description
C          ------------------------------------------------------------
C           -1   Yes  Warning: Some abs(X(i)) so small I1 underflows.
C                     The corresponding F(i) are set to zero.
C            0   Yes  Successfull execution.
C            1   No   Error: M .LE. 0
C            2   No   Error: Some abs(X(i)) so big I1 overflows.
C                     The index of the first offending argument is
C                     returned in IWORK(1).
C
C
C *********************************************************************
C   This routine is a modification of the function DBESI1 developed by
C   W. Fullerton of LANL.
C *********************************************************************
C
C***REFERENCES  (NONE)
C***ROUTINES CALLED  DWI1
C***REVISION HISTORY  (YYMMDD)
C   910226  DATE WRITTEN
C***END PROLOGUE  DVI1
C
C     ----------
C     PARAMETERS
C     ----------
C
      INTEGER INFO, IWORK, M
      DOUBLE PRECISION F, X, WORK
C
      DIMENSION X(M), F(M), WORK(7*M), IWORK(M)
C
C     ---------------
C     LOCAL VARIABLES
C     ---------------
C
      INTEGER IWB0, IWB1, IWB2, IWY, IWTC, IWYC, IWZC, JIN
C
C***FIRST EXECUTABLE STATEMENT  DVI1
C
C     ... PARTITION WORK ARRAYS
C
      IWY   = 1
      IWTC  = IWY  + M
      IWYC  = IWTC + M
      IWZC  = IWYC + M
      IWB0  = IWZC + M
      IWB1  = IWB0 + M
      IWB2  = IWB1 + M
C     Total = IWB2 + M
C
      JIN   = 1
C     Total = JIN  + M
C
C     ... DWI1 DOES ALL THE WORK
C
      CALL DWI1(M,X,F,WORK(IWY),WORK(IWTC),WORK(IWYC),WORK(IWZC),
     +         IWORK(JIN),WORK(IWB0),WORK(IWB1),WORK(IWB2),INFO)
C
      RETURN
      END

      SUBROUTINE DWI1(M, X, F, Y, TCMP, YCMP, ZCMP, INDX, B0, B1, B2, 
     +   INFO)
C***BEGIN PROLOGUE  DWI1
C***SUBSIDIARY
C***PURPOSE  Computes the hyperbolic Bessel function of the first kind
C            of order one (I1) for a vector of arguments
C***LIBRARY   VFNLIB
C***AUTHOR  SAUNDERS, B. V., (NIST)
C             COMPUTING AND APPLIED MATHEMATICS LABORATORY
C             NATIONAL INSTITUTE OF STANDARDS AND TECHNOLOGY
C             GAITHERSBURG, MD 20899
C           BOISVERT, R. F., (NIST)
C             COMPUTING AND APPLIED MATHEMATICS LABORATORY
C             NATIONAL INSTITUTE OF STANDARDS AND TECHNOLOGY
C             GAITHERSBURG, MD 20899
C***ROUTINES CALLED  IDWCS, D1MACH, DWNGT, DWGTHR, DWGTLE, DWGT, DWLE,
C                    DWSCTR, DWCS
C***REVISION HISTORY  (YYMMDD)
C   910226  DATE WRITTEN
C***END PROLOGUE  DWI1
C
C  ----------
C  PARAMETERS
C  ----------
C
      INTEGER LAI1, LAI12, LBI1
      PARAMETER ( LAI1=46, LAI12=69, LBI1=17 )
C
      INTEGER INFO, INDX, M
      DOUBLE PRECISION B0, B1, B2, F, X, Y, TCMP, YCMP, ZCMP
C
      DIMENSION B0(M), B1(M), B2(M), F(M), INDX(M), X(M), Y(M),
     +          TCMP(M), YCMP(M), ZCMP(M)
C
C  ---------------
C  LOCAL VARIABLES
C  ---------------
C
      INTEGER I, IDWCS, J, KEY, N, NTI1, NTAI1, NTAI12, NTOT
      DOUBLE PRECISION AI1CS, AI12CS, BI1CS, EPMACH, EPS, D1MACH,
     +        XMAX, XMIN, XSML
C
      DIMENSION AI1CS(LAI1), AI12CS(LAI12), BI1CS(LBI1)
C
      SAVE AI1CS, AI12CS, BI1CS, N, NTAI1, NTAI12, NTI1,
     +     XMAX, XMIN, XSML
C
C----------------------------------------------------------------------
C
C Series for BI1        ON THE INTERVAL  0.          to  9.00000D+00
C                                        with weighted error   1.44E-32
C                                         log weighted error  31.84
C                               significant figures required  31.45
C                                    decimal places required  32.46
C
      DATA BI1 CS(  1) / -.1971713261 0998597316 1385032181 49 D-2     /
      DATA BI1 CS(  2) / +.4073488766 7546480608 1553936520 14 D+0     /
      DATA BI1 CS(  3) / +.3483899429 9959455866 2450377837 87 D-1     /
      DATA BI1 CS(  4) / +.1545394556 3001236038 5984010584 89 D-2     /
      DATA BI1 CS(  5) / +.4188852109 8377784129 4588320041 20 D-4     /
      DATA BI1 CS(  6) / +.7649026764 8362114741 9597039660 69 D-6     /
      DATA BI1 CS(  7) / +.1004249392 4741178689 1798080372 38 D-7     /
      DATA BI1 CS(  8) / +.9932207791 9238106481 3712980548 63 D-10    /
      DATA BI1 CS(  9) / +.7663801791 8447637275 2001716813 49 D-12    /
      DATA BI1 CS( 10) / +.4741418923 8167394980 3880919481 60 D-14    /
      DATA BI1 CS( 11) / +.2404114404 0745181799 8631720320 00 D-16    /
      DATA BI1 CS( 12) / +.1017150500 7093713649 1211007999 99 D-18    /
      DATA BI1 CS( 13) / +.3645093565 7866949458 4917333333 33 D-21    /
      DATA BI1 CS( 14) / +.1120574950 2562039344 8106666666 66 D-23    /
      DATA BI1 CS( 15) / +.2987544193 4468088832 0000000000 00 D-26    /
      DATA BI1 CS( 16) / +.6973231093 9194709333 3333333333 33 D-29    /
      DATA BI1 CS( 17) / +.1436794822 0620800000 0000000000 00 D-31    /
C
C-------------------------------------------------------------------
C
C Series for AI1        on the interval  1.25000D-01 to  3.33333D-01
C                                        with weighted error   2.81E-32
C                                         log weighted error  31.55
C                               significant figures required  29.93
C                                    decimal places required  32.38
C
      DATA AI1 CS(  1) / -.2846744181 8814786741 0037246830 7 D-1      /
      DATA AI1 CS(  2) / -.1922953231 4432206510 4444877497 9 D-1      /
      DATA AI1 CS(  3) / -.6115185857 9437889822 5624991778 5 D-3      /
      DATA AI1 CS(  4) / -.2069971253 3502277088 8282377797 9 D-4      /
      DATA AI1 CS(  5) / +.8585619145 8107255655 3694467313 8 D-5      /
      DATA AI1 CS(  6) / +.1049498246 7115908625 1745399786 0 D-5      /
      DATA AI1 CS(  7) / -.2918338918 4479022020 9343232669 7 D-6      /
      DATA AI1 CS(  8) / -.1559378146 6317390001 6068096907 7 D-7      /
      DATA AI1 CS(  9) / +.1318012367 1449447055 2530287390 9 D-7      /
      DATA AI1 CS( 10) / -.1448423418 1830783176 3913446781 5 D-8      /
      DATA AI1 CS( 11) / -.2908512243 9931420948 2504099301 0 D-9      /
      DATA AI1 CS( 12) / +.1266388917 8753823873 1115969040 3 D-9      /
      DATA AI1 CS( 13) / -.1664947772 9192206706 2417839858 0 D-10     /
      DATA AI1 CS( 14) / -.1666653644 6094329760 9593715499 9 D-11     /
      DATA AI1 CS( 15) / +.1242602414 2907682652 3216847201 7 D-11     /
      DATA AI1 CS( 16) / -.2731549379 6724323972 5146142863 3 D-12     /
      DATA AI1 CS( 17) / +.2023947881 6458037807 0026268898 1 D-13     /
      DATA AI1 CS( 18) / +.7307950018 1168836361 9869812612 3 D-14     /
      DATA AI1 CS( 19) / -.3332905634 4046749438 1377861713 3 D-14     /
      DATA AI1 CS( 20) / +.7175346558 5129537435 4225466567 0 D-15     /
      DATA AI1 CS( 21) / -.6982530324 7962563558 5062922365 6 D-16     /
      DATA AI1 CS( 22) / -.1299944201 5627607600 6044608058 7 D-16     /
      DATA AI1 CS( 23) / +.8120942864 2427988920 5467834286 0 D-17     /
      DATA AI1 CS( 24) / -.2194016207 4107368981 5626664378 3 D-17     /
      DATA AI1 CS( 25) / +.3630516170 0296548482 7986093233 4 D-18     /
      DATA AI1 CS( 26) / -.1695139772 4391041663 0686679039 9 D-19     /
      DATA AI1 CS( 27) / -.1288184829 8979078071 1688253822 2 D-19     /
      DATA AI1 CS( 28) / +.5694428604 9670527801 0999107310 9 D-20     /
      DATA AI1 CS( 29) / -.1459597009 0904800565 4550990028 7 D-20     /
      DATA AI1 CS( 30) / +.2514546010 6757173140 8469133448 5 D-21     /
      DATA AI1 CS( 31) / -.1844758883 1391248181 6040002901 3 D-22     /
      DATA AI1 CS( 32) / -.6339760596 2279486419 2860979199 9 D-23     /
      DATA AI1 CS( 33) / +.3461441102 0310111111 0814662656 0 D-23     /
      DATA AI1 CS( 34) / -.1017062335 3713935475 9654102357 3 D-23     /
      DATA AI1 CS( 35) / +.2149877147 0904314459 6250077866 6 D-24     /
      DATA AI1 CS( 36) / -.3045252425 2386764017 4620617386 6 D-25     /
      DATA AI1 CS( 37) / +.5238082144 7212859821 7763498666 6 D-27     /
      DATA AI1 CS( 38) / +.1443583107 0893824464 1678950399 9 D-26     /
      DATA AI1 CS( 39) / -.6121302074 8900427332 0067071999 9 D-27     /
      DATA AI1 CS( 40) / +.1700011117 4678184183 4918980266 6 D-27     /
      DATA AI1 CS( 41) / -.3596589107 9842441585 3521578666 6 D-28     /
      DATA AI1 CS( 42) / +.5448178578 9484185766 5051306666 6 D-29     /
      DATA AI1 CS( 43) / -.2731831789 6890849891 6256426666 6 D-30     /
      DATA AI1 CS( 44) / -.1858905021 7086007157 7190399999 9 D-30     /
      DATA AI1 CS( 45) / +.9212682974 5139334411 2776533333 3 D-31     /
      DATA AI1 CS( 46) / -.2813835155 6535611063 7083306666 6 D-31     /
C
C-------------------------------------------------------------------
C
C Series for AI12       on the interval  0.          to  1.25000D-01
C                                        with weighted error   1.83E-32
C                                         log weighted error  31.74
C                               significant figures required  29.97
C                                    decimal places required  32.66
C
      DATA AI12CS(  1) / +.2857623501 8280120474 4984594846 9 D-1      /
      DATA AI12CS(  2) / -.9761097491 3614684077 6516445730 2 D-2      /
      DATA AI12CS(  3) / -.1105889387 6262371629 1256921277 5 D-3      /
      DATA AI12CS(  4) / -.3882564808 8776903934 5654477627 4 D-5      /
      DATA AI12CS(  5) / -.2512236237 8702089252 9452002212 1 D-6      /
      DATA AI12CS(  6) / -.2631468846 8895195068 3705236523 2 D-7      /
      DATA AI12CS(  7) / -.3835380385 9642370220 4500678796 8 D-8      /
      DATA AI12CS(  8) / -.5589743462 1965838068 6811252222 9 D-9      /
      DATA AI12CS(  9) / -.1897495812 3505412344 9892503323 8 D-10     /
      DATA AI12CS( 10) / +.3252603583 0154882385 5508067994 9 D-10     /
      DATA AI12CS( 11) / +.1412580743 6613781331 6336633284 6 D-10     /
      DATA AI12CS( 12) / +.2035628544 1470895072 2452613684 0 D-11     /
      DATA AI12CS( 13) / -.7198551776 2459085120 9258989044 6 D-12     /
      DATA AI12CS( 14) / -.4083551111 0921973182 2849963969 1 D-12     /
      DATA AI12CS( 15) / -.2101541842 7726643130 1984572746 2 D-13     /
      DATA AI12CS( 16) / +.4272440016 7119513542 9778833699 7 D-13     /
      DATA AI12CS( 17) / +.1042027698 4128802764 1741449994 8 D-13     /
      DATA AI12CS( 18) / -.3814403072 4370078047 6707253539 6 D-14     /
      DATA AI12CS( 19) / -.1880354775 5107824485 1273453396 3 D-14     /
      DATA AI12CS( 20) / +.3308202310 9209282827 3190335240 5 D-15     /
      DATA AI12CS( 21) / +.2962628997 6459501390 6854654205 2 D-15     /
      DATA AI12CS( 22) / -.3209525921 9934239587 7837353288 7 D-16     /
      DATA AI12CS( 23) / -.4650305368 4893583255 7128281897 9 D-16     /
      DATA AI12CS( 24) / +.4414348323 0717079499 4611375964 1 D-17     /
      DATA AI12CS( 25) / +.7517296310 8421048054 2545808029 5 D-17     /
      DATA AI12CS( 26) / -.9314178867 3268833756 8484784515 7 D-18     /
      DATA AI12CS( 27) / -.1242193275 1948909561 1678448869 7 D-17     /
      DATA AI12CS( 28) / +.2414276719 4548484690 0515390217 6 D-18     /
      DATA AI12CS( 29) / +.2026944384 0532851789 7192286069 2 D-18     /
      DATA AI12CS( 30) / -.6394267188 2690977870 4391988681 1 D-19     /
      DATA AI12CS( 31) / -.3049812452 3730958960 8488450357 1 D-19     /
      DATA AI12CS( 32) / +.1612841851 6514802251 3462230769 1 D-19     /
      DATA AI12CS( 33) / +.3560913964 3099250545 1027090462 0 D-20     /
      DATA AI12CS( 34) / -.3752017947 9364390796 6682800324 6 D-20     /
      DATA AI12CS( 35) / -.5787037427 0747993459 5198231074 1 D-22     /
      DATA AI12CS( 36) / +.7759997511 6481619619 8236963209 2 D-21     /
      DATA AI12CS( 37) / -.1452790897 2022333940 6445987408 5 D-21     /
      DATA AI12CS( 38) / -.1318225286 7390367021 2192275337 4 D-21     /
      DATA AI12CS( 39) / +.6116654862 9030707018 7999133171 7 D-22     /
      DATA AI12CS( 40) / +.1376279762 4271264277 3024338363 4 D-22     /
      DATA AI12CS( 41) / -.1690837689 9593478849 1983938230 6 D-22     /
      DATA AI12CS( 42) / +.1430596088 5954331539 8720108538 5 D-23     /
      DATA AI12CS( 43) / +.3409557828 0905940204 0536772990 2 D-23     /
      DATA AI12CS( 44) / -.1309457666 2707602278 4573872642 4 D-23     /
      DATA AI12CS( 45) / -.3940706411 2402574360 9352141755 7 D-24     /
      DATA AI12CS( 46) / +.4277137426 9808765808 0616679735 2 D-24     /
      DATA AI12CS( 47) / -.4424634830 9826068819 0028312302 9 D-25     /
      DATA AI12CS( 48) / -.8734113196 2307149721 1530978874 7 D-25     /
      DATA AI12CS( 49) / +.4045401335 6835333921 4340414242 8 D-25     /
      DATA AI12CS( 50) / +.7067100658 0946894656 5160771780 6 D-26     /
      DATA AI12CS( 51) / -.1249463344 5651052230 0286451860 5 D-25     /
      DATA AI12CS( 52) / +.2867392244 4034370329 7948339142 6 D-26     /
      DATA AI12CS( 53) / +.2044292892 5042926702 8177957421 0 D-26     /
      DATA AI12CS( 54) / -.1518636633 8204625683 7134680291 1 D-26     /
      DATA AI12CS( 55) / +.8110181098 1875758861 3227910703 7 D-28     /
      DATA AI12CS( 56) / +.3580379354 7735860911 2717370327 0 D-27     /
      DATA AI12CS( 57) / -.1692929018 9279025095 9305717544 8 D-27     /
      DATA AI12CS( 58) / -.2222902499 7024276390 6775852777 4 D-28     /
      DATA AI12CS( 59) / +.5424535127 1459696550 4860040112 8 D-28     /
      DATA AI12CS( 60) / -.1787068401 5780186887 6491299330 4 D-28     /
      DATA AI12CS( 61) / -.6565479068 7228149388 2392943788 0 D-29     /
      DATA AI12CS( 62) / +.7807013165 0611452809 2206770683 9 D-29     /
      DATA AI12CS( 63) / -.1816595260 6689797173 7933315222 1 D-29     /
      DATA AI12CS( 64) / -.1287704952 6600848203 7687559895 9 D-29     /
      DATA AI12CS( 65) / +.1114548172 9881645474 1370927369 4 D-29     /
      DATA AI12CS( 66) / -.1808343145 0393369391 5936887668 7 D-30     /
      DATA AI12CS( 67) / -.2231677718 2037719522 3244822893 9 D-30     /
      DATA AI12CS( 68) / +.1619029596 0803415106 1790980361 4 D-30     /
      DATA AI12CS( 69) / -.1834079908 8049414139 0130843921 0 D-31     /
C
C-------------------------------------------------------------------
C
      DATA NTI1 / 0 /
C
C***FIRST EXECUTABLE STATEMENT  DWI1
C
      IF (M .LE. 0)  GO TO 910
C
      IF (NTI1.EQ.0) THEN
         EPMACH = D1MACH(3)
         EPS = 0.10D0*EPMACH
         NTI1 = IDWCS (BI1CS, LBI1, EPS)
         NTAI1 = IDWCS (AI1CS, LAI1, EPS)
         NTAI12 = IDWCS (AI12CS, LAI12, EPS)
         XMIN = 2.0D0*D1MACH(1)
         XSML = SQRT(8.0D0*EPMACH)
         XMAX = LOG(D1MACH(2))
      ENDIF
C
      NTOT = 0
C
      DO 10 I=1,M
         Y(I) = ABS(X(I))
  10  CONTINUE
C
      CALL DWNGT(M,Y,XMAX,KEY)
      IF (KEY .NE. 0)  GO TO 920
C
C  ---------------------------
C  CASE   Y=0  OR  Y TOO SMALL
C  ---------------------------
C
C     NOTE -- I0 UNDERFLOWS FOR X .LE. XMIN
C
      DO 15 I=1,M
         F(I) = 0.0D0
  15  CONTINUE
C
C  ----------------------------
C  CASE   XMIN .LT. Y .LE. XSML
C  ----------------------------
C
      CALL DWGTLE(M,Y,XMIN,XSML,N,INDX)
      IF (N .GT. 0) THEN
         NTOT = NTOT + N
         CALL DWGTHR(N,Y,INDX,YCMP)
         DO 20 J=1,N
            ZCMP(J) = 0.50D0*YCMP(J)
  20     CONTINUE
         CALL DWSCTR(N,ZCMP,INDX,F)
      ENDIF
C
C  ---------------------------
C  CASE   XSML .LT. Y .LE. 3.0
C  ---------------------------
C
      CALL DWGTLE(M,Y,XSML,3.0D0,N,INDX)
      IF (N .GT. 0) THEN
         NTOT = NTOT + N
         CALL DWGTHR(N,Y,INDX,YCMP)
         DO 30 J=1,N
            TCMP(J) = YCMP(J)**2/4.50D0 - 1.0D0
  30     CONTINUE
         CALL DWCS(N,TCMP,BI1CS,NTI1,ZCMP,B0,B1,B2)
         DO 40 J=1,N
            ZCMP(J) = YCMP(J)*(0.8750D0 + ZCMP(J))
  40     CONTINUE
         CALL DWSCTR(N,ZCMP,INDX,F)
      ENDIF
C
C  -------------------------
C  CASE  3.0 .LT. Y .LE. 8.0
C  -------------------------
C
      CALL DWGTLE(M,Y,3.0D0,8.0D0,N,INDX)
      IF (N .GT. 0) THEN
         NTOT = NTOT + N
         CALL DWGTHR(N,Y,INDX,YCMP)
         DO 70 J=1,N
            TCMP(J) = (48.0D0/YCMP(J) - 11.0D0)/5.0D0
  70     CONTINUE
         CALL DWCS(N,TCMP,AI1CS,NTAI1,ZCMP,B0,B1,B2)
         DO 80 J=1,N
            ZCMP(J) = EXP(YCMP(J))*(0.3750D0 + ZCMP(J))/SQRT(YCMP(J))
  80     CONTINUE
         CALL DWSCTR(N,ZCMP,INDX,F)
      ENDIF
C
C  ----------------
C  CASE  8.0 .LT. Y
C  ----------------
C
      CALL DWGT(M,Y,8.0D0,N,INDX)
      IF (N .GT. 0) THEN
         NTOT = NTOT + N
         CALL DWGTHR(N,Y,INDX,YCMP)
         DO 100 J=1,N
            TCMP(J) = 16.0D0/YCMP(J) - 1.0D0
  100     CONTINUE
         CALL DWCS(N,TCMP,AI12CS,NTAI12,ZCMP,B0,B1,B2)
         DO 110 J=1,N
            ZCMP(J) = EXP(YCMP(J))*(0.3750D0 + ZCMP(J))/SQRT(YCMP(J))
  110    CONTINUE
         CALL DWSCTR(N,ZCMP,INDX,F)
      ENDIF
C
C  -------------------------------------------
C  REVERSE SIGN FOR NEGATIVE X (RESULT IS ODD)
C  -------------------------------------------
C
      DO 200 I=1,M
         IF (X(I) .LT. 0.0D0)  F(I) = -F(I)
  200 CONTINUE
C
C  -----
C  EXITS
C  -----
C
C     ... NORMAL
C
      INFO = 0
      IF (NTOT .NE. M)  INFO = -1
      GO TO 999
C
C     ... M .LE. 0
C
  910 CONTINUE
      INFO = 1
      GO TO 999
C
C     ... ABS(X) SO LARGE I1 OVERFLOWS
C
  920 CONTINUE
      INFO = 2
      INDX(1) = KEY
      GO TO 999
C
  999 CONTINUE
      RETURN
      END
      SUBROUTINE DVK0 (M, X, F, WORK, IWORK, INFO)
C***BEGIN PROLOGUE  DVK0
C***PURPOSE  Computes the hyperbolic Bessel function of the third kind
C            of order zero (K0) for a vector of real arguments
C***LIBRARY   VFNLIB
C***CATEGORY  C10B1
C***TYPE      DOUBLE PRECISION (VK0-S, DVK0-D)
C***KEYWORDS  BESSEL FUNCTION,THIRD KIND,HYPERBOLIC BESSEL FUNCTION,
C             MODIFIED BESSEL FUNCTION, ORDER ZERO, VECTORIZED
C***AUTHOR  SAUNDERS, B. V., (NIST)
C             COMPUTING AND APPLIED MATHEMATICS LABORATORY
C             NATIONAL INSTITUTE OF STANDARDS AND TECHNOLOGY
C             GAITHERSBURG, MD 20899
C           BOISVERT, R. F., (NIST)
C             COMPUTING AND APPLIED MATHEMATICS LABORATORY
C             NATIONAL INSTITUTE OF STANDARDS AND TECHNOLOGY
C             GAITHERSBURG, MD 20899
C***DESCRIPTION
C
C   DVK0 computes the modified (hyperbolic) Bessel function of the
C   third kind of order zero (K0) for real arguments using uniform
C   approximation by Chebyshev polynomials.
C
C
C   P A R A M E T E R S
C
C   M      (Input) Integer
C          The number of arguments at which the function is to be
C          evaluated.
C
C   X      (Input) Double precision array of length M
C          The arguments at which the function is to be evaluated are
C          stored in X(1) to X(M) in any order.
C
C   F      (Output) Double precision array of length M
C          F(i) contains the value of the function at X(i), i=1,..,M.
C
C   WORK   (Work) Double precision vector of length 7*M
C
C   IWORK  (Work) Integer vector of length M
C
C   INFO   (Output) Integer
C          Indicates status of computed result. The following table
C          lists possible values and their meanings.  If OK=Yes then
C          all F(i) have been set, otherwise none have been set.
C
C          INFO  OK            Description
C          ------------------------------------------------------------
C           -1   Yes  Warning: Some X(i) so big K0 underflows.
C                     The corresponding F(i) are set to zero.
C            0   Yes  Successfull execution.
C            1   No   Error: M .LE. 0
C            2   No   Error: Some X(i) is zero or negative.
C                     The index of the first offending argument is
C                     returned in IWORK(1).
C
C
C *********************************************************************
C   This routine is a modification of the function DBESK0 developed by
C   W. Fullerton of LANL.
C *********************************************************************
C
C***REFERENCES  (NONE)
C***ROUTINES CALLED  DWK0
C***REVISION HISTORY  (YYMMDD)
C   910226  DATE WRITTEN
C***END PROLOGUE  DVK0
C
C     ----------
C     PARAMETERS
C     ----------
C
      INTEGER INFO, IWORK, M
      DOUBLE PRECISION F, X, WORK
C
      DIMENSION X(M), F(M), WORK(7*M), IWORK(M)
C
C     ---------------
C     LOCAL VARIABLES
C     ---------------
C
      INTEGER IWB0, IWB1, IWB2, IWTC, IWXC, IWYC, IWZC, JIN
C
C***FIRST EXECUTABLE STATEMENT  DVK0
C
C     ... PARTITION WORK ARRAYS
C
      IWTC   = 1
      IWXC  = IWTC + M
      IWYC  = IWXC + M
      IWZC  = IWYC + M
      IWB0  = IWZC + M
      IWB1  = IWB0 + M
      IWB2  = IWB1 + M
C     Total = IB2  + M
C
      JIN   = 1
C     Total = JIN  + M
C
C     ... DWK0 DOES ALL THE WORK
C
      CALL DWK0(M,X,F,WORK(IWTC),WORK(IWXC),WORK(IWYC),WORK(IWZC),
     +         IWORK(JIN),WORK(IWB0),WORK(IWB1),WORK(IWB2),INFO)
C
      RETURN
      END


      SUBROUTINE DWK0 (M, X, F, TCMP, XCMP, YCMP, ZCMP, INDX, B0, B1,  
     +   B2, INFO)
C***BEGIN PROLOGUE  DWK0
C***SUBSIDIARY
C***PURPOSE  Computes the hyperbolic Bessel function of the third kind
C            of order zero (K0) for a vector of arguments
C***LIBRARY   VFNLIB
C***AUTHOR  SAUNDERS, B. V., (NIST)
C             COMPUTING AND APPLIED MATHEMATICS LABORATORY
C             NATIONAL INSTITUTE OF STANDARDS AND TECHNOLOGY
C             GAITHERSBURG, MD 20899
C           BOISVERT, R. F., (NIST)
C             COMPUTING AND APPLIED MATHEMATICS LABORATORY
C             NATIONAL INSTITUTE OF STANDARDS AND TECHNOLOGY
C             GAITHERSBURG, MD 20899
C***ROUTINES CALLED  IDWCS, D1MACH, DWNLE, DWGTHR, DWGTLE, DWGT, DWLE,
C                    DWSCTR, DWCS
C***REVISION HISTORY  (YYMMDD)
C   910226  DATE WRITTEN
C***END PROLOGUE  DWK0
C
C  ----------
C  PARAMETERS
C  ----------
C
      INTEGER INFO, INDX, M
      DOUBLE PRECISION B0, B1, B2, F, X, TCMP, XCMP, YCMP, ZCMP
C
      DIMENSION B0(M), B1(M), B2(M), F(M), INDX(M), X(M), TCMP(M),
     +          XCMP(M), YCMP(M), ZCMP(M)
C
C  ---------------
C  LOCAL VARIABLES
C  ---------------
C
      INTEGER LAK0, LAK02, LBI0, LBK0
      PARAMETER ( LAK0=38, LAK02=33, LBI0=18, LBK0=16 )
C
      INTEGER I, IDWCS, J, KEY, N, NTI0, NTK0, NTAK0, NTAK02, NTOT
      DOUBLE PRECISION AK0CS, AK02CS, BI0CS, BK0CS, EPMACH, EPS, D1MACH,
     +        XMAX, XSML
C
      DIMENSION  AK0CS(LAK0), AK02CS(LAK02), BI0CS(LBI0), BK0CS(LBK0)
C
      SAVE BK0CS, AK0CS, AK02CS, BI0CS, N, NTAK0, NTAK02, NTI0, NTK0,
     +     XMAX, XSML
C
C----------------------------------------------------------------------
C
C Series for BK0        on the interval  0.          to  4.00000D+00
C                                        with weighted error   3.08E-33
C                                         log weighted error  32.51
C                               significant figures required  32.05
C                                    decimal places required  33.11
C
      DATA BK0 CS(  1) / -.3532739323 3902768720 1140060063 153 D-1    /
      DATA BK0 CS(  2) / +.3442898999 2462848688 6344927529 213 D+0    /
      DATA BK0 CS(  3) / +.3597993651 5361501626 5721303687 231 D-1    /
      DATA BK0 CS(  4) / +.1264615411 4469259233 8479508673 447 D-2    /
      DATA BK0 CS(  5) / +.2286212103 1194517860 8269830297 585 D-4    /
      DATA BK0 CS(  6) / +.2534791079 0261494573 0790013428 354 D-6    /
      DATA BK0 CS(  7) / +.1904516377 2202088589 7214059381 366 D-8    /
      DATA BK0 CS(  8) / +.1034969525 7633624585 1008317853 089 D-10   /
      DATA BK0 CS(  9) / +.4259816142 7910825765 2445327170 133 D-13   /
      DATA BK0 CS( 10) / +.1374465435 8807508969 4238325440 000 D-15   /
      DATA BK0 CS( 11) / +.3570896528 5083735909 9688597333 333 D-18   /
      DATA BK0 CS( 12) / +.7631643660 1164373766 7498666666 666 D-21   /
      DATA BK0 CS( 13) / +.1365424988 4407818590 8053333333 333 D-23   /
      DATA BK0 CS( 14) / +.2075275266 9066680831 9999999999 999 D-26   /
      DATA BK0 CS( 15) / +.2712814218 0729856000 0000000000 000 D-29   /
      DATA BK0 CS( 16) / +.3082593887 9146666666 6666666666 666 D-32   /
C
C-------------------------------------------------------------------
C
C Series for AK0        on the interval  1.25000D-01 to  5.00000D-01
C                                        with weighted error   2.85E-32
C                                         log weighted error  31.54
C                               significant figures required  30.19
C                                    decimal places required  32.33
C
      DATA AK0 CS(  1) / -.7643947903 3279414240 8297827008 8 D-1      /
      DATA AK0 CS(  2) / -.2235652605 6998190520 2309555079 1 D-1      /
      DATA AK0 CS(  3) / +.7734181154 6938582353 0061817404 7 D-3      /
      DATA AK0 CS(  4) / -.4281006688 8860994644 5214643541 6 D-4      /
      DATA AK0 CS(  5) / +.3081700173 8629747436 5001482666 0 D-5      /
      DATA AK0 CS(  6) / -.2639367222 0096649740 6744889272 3 D-6      /
      DATA AK0 CS(  7) / +.2563713036 4034692062 9408826574 2 D-7      /
      DATA AK0 CS(  8) / -.2742705549 9002012638 5721191524 4 D-8      /
      DATA AK0 CS(  9) / +.3169429658 0974995920 8083287340 3 D-9      /
      DATA AK0 CS( 10) / -.3902353286 9621841416 0106571796 2 D-10     /
      DATA AK0 CS( 11) / +.5068040698 1885754020 5009212728 6 D-11     /
      DATA AK0 CS( 12) / -.6889574741 0078706795 4171355798 4 D-12     /
      DATA AK0 CS( 13) / +.9744978497 8259176913 8820133683 1 D-13     /
      DATA AK0 CS( 14) / -.1427332841 8845485053 8985534012 2 D-13     /
      DATA AK0 CS( 15) / +.2156412571 0214630395 5806297652 7 D-14     /
      DATA AK0 CS( 16) / -.3349654255 1495627721 8878205853 0 D-15     /
      DATA AK0 CS( 17) / +.5335260216 9529116921 4528039260 1 D-16     /
      DATA AK0 CS( 18) / -.8693669980 8907538076 3962237883 7 D-17     /
      DATA AK0 CS( 19) / +.1446404347 8622122278 8776344234 6 D-17     /
      DATA AK0 CS( 20) / -.2452889825 5001296824 0467875157 3 D-18     /
      DATA AK0 CS( 21) / +.4233754526 2321715728 2170634240 0 D-19     /
      DATA AK0 CS( 22) / -.7427946526 4544641956 9534129493 3 D-20     /
      DATA AK0 CS( 23) / +.1323150529 3926668662 7796746240 0 D-20     /
      DATA AK0 CS( 24) / -.2390587164 7396494513 3598146559 9 D-21     /
      DATA AK0 CS( 25) / +.4376827585 9232261401 6571255466 6 D-22     /
      DATA AK0 CS( 26) / -.8113700607 3451180593 3901141333 3 D-23     /
      DATA AK0 CS( 27) / +.1521819913 8321729583 1037815466 6 D-23     /
      DATA AK0 CS( 28) / -.2886041941 4833977702 3595861333 3 D-24     /
      DATA AK0 CS( 29) / +.5530620667 0547179799 9261013333 3 D-25     /
      DATA AK0 CS( 30) / -.1070377329 2498987285 9163306666 6 D-25     /
      DATA AK0 CS( 31) / +.2091086893 1423843002 9632853333 3 D-26     /
      DATA AK0 CS( 32) / -.4121713723 6462038274 1026133333 3 D-27     /
      DATA AK0 CS( 33) / +.8193483971 1213076401 3568000000 0 D-28     /
      DATA AK0 CS( 34) / -.1642000275 4592977267 8075733333 3 D-28     /
      DATA AK0 CS( 35) / +.3316143281 4802271958 9034666666 6 D-29     /
      DATA AK0 CS( 36) / -.6746863644 1452959410 8586666666 6 D-30     /
      DATA AK0 CS( 37) / +.1382429146 3184246776 3541333333 3 D-30     /
      DATA AK0 CS( 38) / -.2851874167 3598325708 1173333333 3 D-31     /
C
C-------------------------------------------------------------------
C
C Series for AK02       on the interval  0.          to  1.25000D-01
C                                        with weighted error   2.30E-32
C                                         log weighted error  31.64
C                               significant figures required  29.68
C                                    decimal places required  32.40
C
      DATA AK02CS(  1) / -.1201869826 3075922398 3934621245 2 D-1      /
      DATA AK02CS(  2) / -.9174852691 0256953106 5256107571 3 D-2      /
      DATA AK02CS(  3) / +.1444550931 7750058210 4884387805 7 D-3      /
      DATA AK02CS(  4) / -.4013614175 4357097286 7102107787 9 D-5      /
      DATA AK02CS(  5) / +.1567831810 8523106725 9034899033 3 D-6      /
      DATA AK02CS(  6) / -.7770110438 5217377103 1579975446 0 D-8      /
      DATA AK02CS(  7) / +.4611182576 1797178825 3313052958 6 D-9      /
      DATA AK02CS(  8) / -.3158592997 8605657705 2666580330 9 D-10     /
      DATA AK02CS(  9) / +.2435018039 3650411278 3588781432 9 D-11     /
      DATA AK02CS( 10) / -.2074331387 3983478977 0985337350 6 D-12     /
      DATA AK02CS( 11) / +.1925787280 5899170847 4273650469 3 D-13     /
      DATA AK02CS( 12) / -.1927554805 8389561036 0034718221 8 D-14     /
      DATA AK02CS( 13) / +.2062198029 1978182782 8523786964 4 D-15     /
      DATA AK02CS( 14) / -.2341685117 5792424026 0364019507 1 D-16     /
      DATA AK02CS( 15) / +.2805902810 6430422468 1517882845 8 D-17     /
      DATA AK02CS( 16) / -.3530507631 1618079458 1548246357 3 D-18     /
      DATA AK02CS( 17) / +.4645295422 9351082674 2421633706 6 D-19     /
      DATA AK02CS( 18) / -.6368625941 3442664739 2205346133 3 D-20     /
      DATA AK02CS( 19) / +.9069521310 9865155676 2234880000 0 D-21     /
      DATA AK02CS( 20) / -.1337974785 4236907398 4500531199 9 D-21     /
      DATA AK02CS( 21) / +.2039836021 8599523155 2208896000 0 D-22     /
      DATA AK02CS( 22) / -.3207027481 3678405000 6086997333 3 D-23     /
      DATA AK02CS( 23) / +.5189744413 6623099636 2635946666 6 D-24     /
      DATA AK02CS( 24) / -.8629501497 5405721929 6460799999 9 D-25     /
      DATA AK02CS( 25) / +.1472161183 1025598552 0803840000 0 D-25     /
      DATA AK02CS( 26) / -.2573069023 8670112838 1235199999 9 D-26     /
      DATA AK02CS( 27) / +.4601774086 6435165873 7664000000 0 D-27     /
      DATA AK02CS( 28) / -.8411555324 2010937371 3066666666 6 D-28     /
      DATA AK02CS( 29) / +.1569806306 6353689393 0154666666 6 D-28     /
      DATA AK02CS( 30) / -.2988226453 0057577889 7919999999 9 D-29     /
      DATA AK02CS( 31) / +.5796831375 2168365206 1866666666 6 D-30     /
      DATA AK02CS( 32) / -.1145035994 3476813321 5573333333 3 D-30     /
      DATA AK02CS( 33) / +.2301266594 2496828020 0533333333 3 D-31     /
C
C-------------------------------------------------------------------
C
C Series for BI0        on the interval  0.          to  9.00000D+00
C                                        with weighted error   9.51E-34
C                                         log weighted error  33.02
C                               significant figures required  33.31
C                                    decimal places required  33.65
C
      DATA BI0 CS(  1) / -.7660547252 8391449510 8189497624 3285 D-1   /
      DATA BI0 CS(  2) / +.1927337953 9938082699 5240875088 1196 D+1   /
      DATA BI0 CS(  3) / +.2282644586 9203013389 3702929233 0415 D+0   /
      DATA BI0 CS(  4) / +.1304891466 7072904280 7933421069 1888 D-1   /
      DATA BI0 CS(  5) / +.4344270900 8164874513 7868268102 6107 D-3   /
      DATA BI0 CS(  6) / +.9422657686 0019346639 2317174411 8766 D-5   /
      DATA BI0 CS(  7) / +.1434006289 5106910799 6209187817 9957 D-6   /
      DATA BI0 CS(  8) / +.1613849069 6617490699 1541971999 4611 D-8   /
      DATA BI0 CS(  9) / +.1396650044 5356696994 9509270814 2522 D-10  /
      DATA BI0 CS( 10) / +.9579451725 5054453446 2752317189 3333 D-13  /
      DATA BI0 CS( 11) / +.5333981859 8625021310 1510774400 0000 D-15  /
      DATA BI0 CS( 12) / +.2458716088 4374707746 9678591999 9999 D-17  /
      DATA BI0 CS( 13) / +.9535680890 2487700269 4434133333 3333 D-20  /
      DATA BI0 CS( 14) / +.3154382039 7214273367 8933333333 3333 D-22  /
      DATA BI0 CS( 15) / +.9004564101 0946374314 6666666666 6666 D-25  /
      DATA BI0 CS( 16) / +.2240647369 1236700160 0000000000 0000 D-27  /
      DATA BI0 CS( 17) / +.4903034603 2428373333 3333333333 3333 D-30  /
      DATA BI0 CS( 18) / +.9508172606 1226666666 6666666666 6666 D-33  /

C
C-------------------------------------------------------------------
C
      DATA NTK0 / 0 /
C
C***FIRST EXECUTABLE STATEMENT  DWK0
C
      IF (M .LE. 0)  GO TO 910
C
      IF (NTK0 .EQ. 0) THEN
         EPMACH = D1MACH(3)
         EPS = 0.10D0*EPMACH
         NTK0 = IDWCS(BK0CS, LBK0, EPS)
         NTAK0 = IDWCS(AK0CS, LAK0, EPS)
         NTAK02 = IDWCS(AK02CS, LAK02, EPS)
         NTI0 = IDWCS(BI0CS, LBI0, EPS)
         XSML = SQRT (4.0D0*EPMACH)
         XMAX = -LOG(D1MACH(1))
         XMAX = XMAX - 0.50D0*XMAX*LOG(XMAX)/(XMAX+0.50D0) - 0.010D0
      ENDIF
C
      NTOT = 0
C
      CALL DWNLE(M,X,0.0D0,KEY)
      IF (KEY .NE. 0)  GO TO 920
C
C  ------------------
C  CASE   X .GT. XMAX
C  ------------------
C
C     NOTE -- K0 UNDERFLOWS FOR X .GT. XMAX
C
      DO 5 I=1,M
         F(I) = 0.0D0
   5  CONTINUE
C
C  ----------------
C  CASE  X .LE. 2.0
C  ----------------
C
      CALL DWLE(M,X,2.0D0,N,INDX)
      IF (N .GT. 0) THEN
         NTOT = NTOT + N
         CALL DWGTHR(N,X,INDX,XCMP)
C
C        ... COMPUTE I0(X) ... RESULT IN ZCMP
C
         DO 10 J=1,N
            TCMP(J) = XCMP(J)**2/4.50D0 - 1.0D0
   10    CONTINUE
         CALL DWCS(N,TCMP,BI0CS,NTI0,ZCMP,B0,B1,B2)
         DO 15 J=1,N
            ZCMP(J) = 2.750D0 + ZCMP(J)
   15    CONTINUE
C
         DO 20 J=1,N
            TCMP(J) = 0.50D0*XCMP(J)**2 - 1.0D0
   20    CONTINUE
         CALL DWCS(N,TCMP,BK0CS,NTK0,YCMP,B0,B1,B2)
         DO 30 J=1,N
            YCMP(J) = -LOG(0.50D0*XCMP(J))*ZCMP(J) - 0.250D0 + YCMP(J)
   30    CONTINUE
         CALL DWSCTR(N,YCMP,INDX,F)
      ENDIF
C
C  -------------------------
C  CASE  2.0 .LT. X .LE. 8.0
C  -------------------------
C
      CALL DWGTLE(M,X,2.0D0,8.0D0,N,INDX)
      IF (N .GT. 0) THEN
         NTOT = NTOT + N
         CALL DWGTHR(N,X,INDX,XCMP)
         DO 40 J=1,N
            TCMP(J) = (16.0D0/XCMP(J) - 5.0D0)/3.0D0
   40    CONTINUE
         CALL DWCS(N,TCMP,AK0CS,NTAK0,ZCMP,B0,B1,B2)
         DO 50 J=1,N
            YCMP(J) = EXP(-XCMP(J))*(1.250D0 + ZCMP(J))/SQRT(XCMP(J))
   50    CONTINUE
         CALL DWSCTR(N,YCMP,INDX,F)
      ENDIF
C
C  --------------------------
C  CASE  8.0 .LT. X .LE. XMAX
C  --------------------------
C
      CALL DWGTLE(M,X,8.0D0,XMAX,N,INDX)
      IF (N .GT. 0) THEN
         NTOT = NTOT + N
         CALL DWGTHR(N,X,INDX,XCMP)
         DO 60 J=1,N
            TCMP(J) = 16.0D0/XCMP(J) - 1.0D0
   60    CONTINUE
         CALL DWCS(N,TCMP,AK02CS,NTAK02,ZCMP,B0,B1,B2)
         DO 70 J=1,N
            YCMP(J) = EXP(-XCMP(J))*(1.250D0 + ZCMP(J))/SQRT(XCMP(J))
   70    CONTINUE
         CALL DWSCTR(N,YCMP,INDX,F)
      ENDIF
C
C  -----
C  EXITS
C  -----
C
C     ... NORMAL
C
      INFO = 0
      GO TO 999
C
C     ... M .LE. 0
C
  910 CONTINUE
      INFO = 1
      IF (NTOT .NE. M)  INFO = -1
      GO TO 999
C
C     ... X IS ZERO OR NEGATIVE
C
  920 CONTINUE
      INFO = 2
      INDX(1) = KEY
      GO TO 999
C
  999 CONTINUE
      RETURN
      END


      SUBROUTINE DVK1 (M, X, F, WORK, IWORK, INFO)
C***BEGIN PROLOGUE  DVK1
C***PURPOSE  Computes the hyperbolic Bessel function of the third kind
C            of order one (K1) for a vector of real arguments
C***LIBRARY   VFNLIB
C***CATEGORY  C10B1
C***TYPE      DOUBLE PRECISION (VK1-S, DVK1-D)
C***KEYWORDS  BESSEL FUNCTION,THIRD KIND,HYPERBOLIC BESSEL FUNCTION,
C             MODIFIED BESSEL FUNCTION, ORDER ONE, VECTORIZED
C***AUTHOR  SAUNDERS, B. V., (NIST)
C             COMPUTING AND APPLIED MATHEMATICS LABORATORY
C             NATIONAL INSTITUTE OF STANDARDS AND TECHNOLOGY
C             GAITHERSBURG, MD 20899
C           BOISVERT, R. F., (NIST)
C             COMPUTING AND APPLIED MATHEMATICS LABORATORY
C             NATIONAL INSTITUTE OF STANDARDS AND TECHNOLOGY
C             GAITHERSBURG, MD 20899
C***DESCRIPTION
C
C   DVK1 computes the modified (hyperbolic) Bessel function of the
C   third kind of order one (K1) for real arguments using uniform
C   approximation by Chebyshev polynomials.
C
C
C   P A R A M E T E R S
C
C   M      (Input) Integer
C          The number of arguments at which the function is to be
C          evaluated.
C
C   X      (Input) Double precision array of length M
C          The arguments at which the function is to be evaluated are
C          stored in X(1) to X(M) in any order.
C
C   F      (Output) Double precision array of length M
C          F(i) contains the value of the function at X(i), i=1,..,M.
C
C   WORK   (Work) Double precision vector of length 7*M
C
C   IWORK  (Work) Integer vector of length M
C
C   INFO   (Output) Integer
C          Indicates status of computed result. The following table
C          lists possible values and their meanings.  If OK=Yes then
C          all F(i) have been set, otherwise none have been set.
C
C          INFO  OK            Description
C          ------------------------------------------------------------
C           -1   Yes  Warning: Some X(i) so big K1 underflows.
C                     The corresponding F(i) are set to zero.
C            0   Yes  Successfull execution.
C            1   No   Error: M .LE. 0
C            2   No   Error: Some X(i) is zero or negative.
C                     The index of the first offending argument is
C                     returned in IWORK(1).
C            3   No   Error: Some X(i) so small K1 overflows.
C                     The index of the first offending argument is
C                     returned in IWORK(1).
C
C
C *********************************************************************
C   This routine is a modification of the function DBESK1 developed by
C   W. Fullerton of LANL.
C *********************************************************************
C
C***REFERENCES  (NONE)
C***ROUTINES CALLED  DWK1
C***REVISION HISTORY  (YYMMDD)
C   910226  DATE WRITTEN
C***END PROLOGUE  DVK1
C
C     ----------
C     PARAMETERS
C     ----------
C
      INTEGER INFO, IWORK, M
      DOUBLE PRECISION F, X, WORK
C
      DIMENSION X(M), F(M), WORK(7*M), IWORK(M)
C
C     ---------------
C     LOCAL VARIABLES
C     ---------------
C
      INTEGER IWB0, IWB1, IWB2, IWTC, IWXC, IWYC, IWZC, JIN
C
C***FIRST EXECUTABLE STATEMENT  DVK1
C
C     ... PARTITION WORK ARRAYS
C
      IWTC   = 1
      IWXC  = IWTC + M
      IWYC  = IWXC + M
      IWZC  = IWYC + M
      IWB0  = IWZC + M
      IWB1  = IWB0 + M
      IWB2  = IWB1 + M
C     Total = IB2  + M
C
      JIN   = 1
C     Total = JIN  + M
C
C     ... DWK1 DOES ALL THE WORK
C
      CALL DWK1(M,X,F,WORK(IWTC),WORK(IWXC),WORK(IWYC),WORK(IWZC),
     +         IWORK(JIN),WORK(IWB0),WORK(IWB1),WORK(IWB2),INFO)
C
      RETURN
      END

      SUBROUTINE DWK1 (M, X, F, TCMP, XCMP, YCMP, ZCMP, INDX, B0, B1,  
     +   B2, INFO)
C***BEGIN PROLOGUE  DWK1
C***SUBSIDIARY
C***PURPOSE  Computes the hyperbolic Bessel function of the third kind
C            of order one (K1) for a vector of arguments
C***LIBRARY   VFNLIB
C***AUTHOR  SAUNDERS, B. V., (NIST)
C             COMPUTING AND APPLIED MATHEMATICS LABORATORY
C             NATIONAL INSTITUTE OF STANDARDS AND TECHNOLOGY
C             GAITHERSBURG, MD 20899
C           BOISVERT, R. F., (NIST)
C             COMPUTING AND APPLIED MATHEMATICS LABORATORY
C             NATIONAL INSTITUTE OF STANDARDS AND TECHNOLOGY
C             GAITHERSBURG, MD 20899
C***ROUTINES CALLED  IDWCS, D1MACH, DWNLE, DWGTHR, DWGTLE, DWGT, DWLE,
C                    DWSCTR, DWCS
C***REVISION HISTORY  (YYMMDD)
C   910226  DATE WRITTEN
C***END PROLOGUE  DWK1
C
C  ----------
C  PARAMETERS
C  ----------
C
      INTEGER INFO, INDX, M
      DOUBLE PRECISION B0, B1, B2, F, X, TCMP, XCMP, YCMP, ZCMP
C
      DIMENSION B0(M), B1(M), B2(M), F(M), INDX(M), X(M), TCMP(M),
     +          XCMP(M), YCMP(M), ZCMP(M)
C
C  ---------------
C  LOCAL VARIABLES
C  ---------------
C
      INTEGER LAK1, LAK12, LBI1, LBK1
      PARAMETER ( LAK1=38, LAK12=33, LBI1=17, LBK1=16 )
C
      INTEGER I, IDWCS, J, KEY, N, NTI1, NTK1, NTAK1, NTAK12, NTOT
      DOUBLE PRECISION AK1CS, AK12CS, BI1CS, BK1CS, EPMACH, EPS, D1MACH,
     +        XMAX, XMIN, XSML
C
      DIMENSION AK1CS(LAK1), AK12CS(LAK12), BI1CS(LBI1), BK1CS(LBK1)
C
      SAVE AK1CS, AK12CS, BI1CS, BK1CS, NTAK1, NTAK12, NTI1, NTK1,
     +     XMAX, XMIN, XSML
C
C----------------------------------------------------------------------
C
C Series for BK1        on the interval  0.          to  4.00000D+00
C                                        with weighted error   9.16E-32
C                                         log weighted error  31.04
C                               significant figures required  30.61
C                                    decimal places required  31.64
C
      DATA BK1 CS(  1) / +.2530022733 8947770532 5311208685 33 D-1     /
      DATA BK1 CS(  2) / -.3531559607 7654487566 7238316918 01 D+0     /
      DATA BK1 CS(  3) / -.1226111808 2265714823 4790679300 42 D+0     /
      DATA BK1 CS(  4) / -.6975723859 6398643501 8129202960 83 D-2     /
      DATA BK1 CS(  5) / -.1730288957 5130520630 1765073689 79 D-3     /
      DATA BK1 CS(  6) / -.2433406141 5659682349 6007350301 64 D-5     /
      DATA BK1 CS(  7) / -.2213387630 7347258558 3152525451 26 D-7     /
      DATA BK1 CS(  8) / -.1411488392 6335277610 9583302126 08 D-9     /
      DATA BK1 CS(  9) / -.6666901694 1993290060 8537512643 73 D-12    /
      DATA BK1 CS( 10) / -.2427449850 5193659339 2631968648 53 D-14    /
      DATA BK1 CS( 11) / -.7023863479 3862875971 7837971200 00 D-17    /
      DATA BK1 CS( 12) / -.1654327515 5100994675 4910293333 33 D-19    /
      DATA BK1 CS( 13) / -.3233834745 9944491991 8933333333 33 D-22    /
      DATA BK1 CS( 14) / -.5331275052 9265274999 4666666666 66 D-25    /
      DATA BK1 CS( 15) / -.7513040716 2157226666 6666666666 66 D-28    /
      DATA BK1 CS( 16) / -.9155085717 6541866666 6666666666 66 D-31    /
C
C-------------------------------------------------------------------
C
C Series for BI1        ON THE INTERVAL  0.          to  9.00000D+00
C                                        with weighted error   1.44E-32
C                                         log weighted error  31.84
C                               significant figures required  31.45
C                                    decimal places required  32.46
C
      DATA BI1 CS(  1) / -.1971713261 0998597316 1385032181 49 D-2     /
      DATA BI1 CS(  2) / +.4073488766 7546480608 1553936520 14 D+0     /
      DATA BI1 CS(  3) / +.3483899429 9959455866 2450377837 87 D-1     /
      DATA BI1 CS(  4) / +.1545394556 3001236038 5984010584 89 D-2     /
      DATA BI1 CS(  5) / +.4188852109 8377784129 4588320041 20 D-4     /
      DATA BI1 CS(  6) / +.7649026764 8362114741 9597039660 69 D-6     /
      DATA BI1 CS(  7) / +.1004249392 4741178689 1798080372 38 D-7     /
      DATA BI1 CS(  8) / +.9932207791 9238106481 3712980548 63 D-10    /
      DATA BI1 CS(  9) / +.7663801791 8447637275 2001716813 49 D-12    /
      DATA BI1 CS( 10) / +.4741418923 8167394980 3880919481 60 D-14    /
      DATA BI1 CS( 11) / +.2404114404 0745181799 8631720320 00 D-16    /
      DATA BI1 CS( 12) / +.1017150500 7093713649 1211007999 99 D-18    /
      DATA BI1 CS( 13) / +.3645093565 7866949458 4917333333 33 D-21    /
      DATA BI1 CS( 14) / +.1120574950 2562039344 8106666666 66 D-23    /
      DATA BI1 CS( 15) / +.2987544193 4468088832 0000000000 00 D-26    /
      DATA BI1 CS( 16) / +.6973231093 9194709333 3333333333 33 D-29    /
      DATA BI1 CS( 17) / +.1436794822 0620800000 0000000000 00 D-31    /
C
C-------------------------------------------------------------------
C
C Series for AK1        on the interval  1.25000D-01 to  5.00000D-01
C                                        with weighted error   3.07E-32
C                                         log weighted error  31.51
C                               significant figures required  30.71
C                                    decimal places required  32.30
C
      DATA AK1 CS(  1) / +.2744313406 9738829695 2576662272 66 D+0     /
      DATA AK1 CS(  2) / +.7571989953 1993678170 8923781492 90 D-1     /
      DATA AK1 CS(  3) / -.1441051556 4754061229 8531161756 25 D-2     /
      DATA AK1 CS(  4) / +.6650116955 1257479394 2513854770 36 D-4     /
      DATA AK1 CS(  5) / -.4369984709 5201407660 5808450891 67 D-5     /
      DATA AK1 CS(  6) / +.3540277499 7630526799 4171390085 34 D-6     /
      DATA AK1 CS(  7) / -.3311163779 2932920208 9826882457 04 D-7     /
      DATA AK1 CS(  8) / +.3445977581 9010534532 3114997709 92 D-8     /
      DATA AK1 CS(  9) / -.3898932347 4754271048 9819374927 58 D-9     /
      DATA AK1 CS( 10) / +.4720819750 4658356400 9474493390 05 D-10    /
      DATA AK1 CS( 11) / -.6047835662 8753562345 3735915628 90 D-11    /
      DATA AK1 CS( 12) / +.8128494874 8658747888 1938379856 63 D-12    /
      DATA AK1 CS( 13) / -.1138694574 7147891428 9239159510 42 D-12    /
      DATA AK1 CS( 14) / +.1654035840 8462282325 9729482050 90 D-13    /
      DATA AK1 CS( 15) / -.2480902567 7068848221 5160104405 33 D-14    /
      DATA AK1 CS( 16) / +.3829237890 7024096948 4292272991 57 D-15    /
      DATA AK1 CS( 17) / -.6064734104 0012418187 7682103773 86 D-16    /
      DATA AK1 CS( 18) / +.9832425623 2648616038 1940046506 66 D-17    /
      DATA AK1 CS( 19) / -.1628416873 8284380035 6666201156 26 D-17    /
      DATA AK1 CS( 20) / +.2750153649 6752623718 2841203370 66 D-18    /
      DATA AK1 CS( 21) / -.4728966646 3953250924 2810695680 00 D-19    /
      DATA AK1 CS( 22) / +.8268150002 8109932722 3920503466 66 D-20    /
      DATA AK1 CS( 23) / -.1468140513 6624956337 1939648853 33 D-20    /
      DATA AK1 CS( 24) / +.2644763926 9208245978 0858948266 66 D-21    /
      DATA AK1 CS( 25) / -.4829015756 4856387897 9698688000 00 D-22    /
      DATA AK1 CS( 26) / +.8929302074 3610130180 6563327999 99 D-23    /
      DATA AK1 CS( 27) / -.1670839716 8972517176 9977514666 66 D-23    /
      DATA AK1 CS( 28) / +.3161645603 4040694931 3686186666 66 D-24    /
      DATA AK1 CS( 29) / -.6046205531 2274989106 5064106666 66 D-25    /
      DATA AK1 CS( 30) / +.1167879894 2042732700 7184213333 33 D-25    /
      DATA AK1 CS( 31) / -.2277374158 2653996232 8678400000 00 D-26    /
      DATA AK1 CS( 32) / +.4481109730 0773675795 3058133333 33 D-27    /
      DATA AK1 CS( 33) / -.8893288476 9020194062 3360000000 00 D-28    /
      DATA AK1 CS( 34) / +.1779468001 8850275131 3920000000 00 D-28    /
      DATA AK1 CS( 35) / -.3588455596 7329095821 9946666666 66 D-29    /
      DATA AK1 CS( 36) / +.7290629049 2694257991 6799999999 99 D-30    /
      DATA AK1 CS( 37) / -.1491844984 5546227073 0240000000 00 D-30    /
      DATA AK1 CS( 38) / +.3073657387 2934276300 7999999999 99 D-31    /
C
C-------------------------------------------------------------------
C
C Series for AK12       on the interval  0.          to  1.25000D-01
C                                        with weighted error   2.41E-32
C                                         log weighted error  31.62
C                               significant figures required  30.25
C                                    decimal places required  32.38
C
      DATA AK12CS(  1) / +.6379308343 7390010366 0048853410 2 D-1      /
      DATA AK12CS(  2) / +.2832887813 0497209358 3503028470 8 D-1      /
      DATA AK12CS(  3) / -.2475370673 9052503454 1454556673 2 D-3      /
      DATA AK12CS(  4) / +.5771972451 6072488204 7097662576 3 D-5      /
      DATA AK12CS(  5) / -.2068939219 5365483027 4553319655 2 D-6      /
      DATA AK12CS(  6) / +.9739983441 3818041803 0921309788 7 D-8      /
      DATA AK12CS(  7) / -.5585336140 3806249846 8889551112 9 D-9      /
      DATA AK12CS(  8) / +.3732996634 0461852402 2121285473 1 D-10     /
      DATA AK12CS(  9) / -.2825051961 0232254451 3506575492 8 D-11     /
      DATA AK12CS( 10) / +.2372019002 4841441736 4349695548 6 D-12     /
      DATA AK12CS( 11) / -.2176677387 9917539792 6830166793 8 D-13     /
      DATA AK12CS( 12) / +.2157914161 6160324539 3956268970 6 D-14     /
      DATA AK12CS( 13) / -.2290196930 7182692759 9155133815 4 D-15     /
      DATA AK12CS( 14) / +.2582885729 8232749619 1993956522 6 D-16     /
      DATA AK12CS( 15) / -.3076752641 2684631876 2109817344 0 D-17     /
      DATA AK12CS( 16) / +.3851487721 2804915970 9489684479 9 D-18     /
      DATA AK12CS( 17) / -.5044794897 6415289771 1728250880 0 D-19     /
      DATA AK12CS( 18) / +.6888673850 4185442370 1829222399 9 D-20     /
      DATA AK12CS( 19) / -.9775041541 9501183030 0213248000 0 D-21     /
      DATA AK12CS( 20) / +.1437416218 5238364610 0165973333 3 D-21     /
      DATA AK12CS( 21) / -.2185059497 3443473734 9973333333 3 D-22     /
      DATA AK12CS( 22) / +.3426245621 8092206316 4538880000 0 D-23     /
      DATA AK12CS( 23) / -.5531064394 2464082325 0124800000 0 D-24     /
      DATA AK12CS( 24) / +.9176601505 6859954037 8282666666 6 D-25     /
      DATA AK12CS( 25) / -.1562287203 6180249114 4874666666 6 D-25     /
      DATA AK12CS( 26) / +.2725419375 4843331323 4943999999 9 D-26     /
      DATA AK12CS( 27) / -.4865674910 0748279923 7802666666 6 D-27     /
      DATA AK12CS( 28) / +.8879388552 7235025873 5786666666 6 D-28     /
      DATA AK12CS( 29) / -.1654585918 0392575489 3653333333 3 D-28     /
      DATA AK12CS( 30) / +.3145111321 3578486743 0399999999 9 D-29     /
      DATA AK12CS( 31) / -.6092998312 1931276124 1600000000 0 D-30     /
      DATA AK12CS( 32) / +.1202021939 3698158346 2399999999 9 D-30     /
      DATA AK12CS( 33) / -.2412930801 4594088413 8666666666 6 D-31     /
C
C-------------------------------------------------------------------
C
      DATA NTK1 / 0 /
C
C***FIRST EXECUTABLE STATEMENT  DWK1
C
      IF (M .LE. 0)  GO TO 910
C
      IF (NTK1 .EQ. 0) THEN
         EPMACH = D1MACH(3)
         EPS = 0.10D0*EPMACH
         NTK1 = IDWCS(BK1CS, LBK1, EPS)
         NTI1 = IDWCS(BI1CS, LBI1, EPS)
         NTAK1 = IDWCS(AK1CS, LAK1, EPS)
         NTAK12 = IDWCS(AK12CS, LAK12, EPS)
         XMIN = EXP(MAX( LOG(D1MACH(1)), -LOG(D1MACH(2))) + 0.010D0)
         XSML = SQRT(4.0D0*EPMACH)
         XMAX = -LOG(D1MACH(1))
         XMAX = XMAX - 0.50D0*XMAX*LOG(XMAX)/(XMAX + 0.50D0)
      ENDIF
C
      NTOT = 0
C
      CALL DWNLE(M,X,0.0D0,KEY)
      IF (KEY .NE. 0)  GO TO 920
C
      CALL DWNLE(M,X,XMIN,KEY)
      IF (KEY .NE. 0)  GO TO 930
C
C  ------------------
C  CASE   X .GT. XMAX
C  ------------------
C
C     NOTE -- K0 UNDERFLOWS FOR X .GT. XMAX
C
      DO 5 I=1,M
         F(I) = 0.0D0
   5  CONTINUE
C
C  ----------------
C  CASE  X .LE. 2.0
C  ----------------
C
      CALL DWLE(M,X,2.0D0,N,INDX)
      IF (N .GT. 0) THEN
         NTOT = NTOT + N
         CALL DWGTHR(N,X,INDX,XCMP)
C
C        ... COMPUTE I1(X) ... RESULT IN ZCMP
C
         DO 20 J=1,N
            TCMP(J) = XCMP(J)**2/4.50D0 - 1.0D0
   20    CONTINUE
         CALL DWCS(N,TCMP,BI1CS,NTI1,ZCMP,B0,B1,B2)
         DO 30 J=1,N
            ZCMP(J) = XCMP(J)*(0.8750D0 + ZCMP(J))
   30    CONTINUE
C
         DO 40 J=1,N
            TCMP(J) = 0.50D0*XCMP(J)**2 - 1.0D0
   40    CONTINUE
         CALL DWCS(N,TCMP,BK1CS,NTK1,YCMP,B0,B1,B2)
         DO 50 J=1,N
            ZCMP(J) = LOG(0.50D0*XCMP(J))*ZCMP(J) +
     +                (0.750D0 + YCMP(J))/XCMP(J)
   50    CONTINUE
         CALL DWSCTR(N,ZCMP,INDX,F)
      ENDIF
C
C  -------------------------
C  CASE  2.0 .LT. X .LE. 8.0
C  -------------------------
C
      CALL DWGTLE(M,X,2.0D0,8.0D0,N,INDX)
      IF (N .GT. 0) THEN
         NTOT = NTOT + N
         CALL DWGTHR(N,X,INDX,XCMP)
         DO 60 J=1,N
            TCMP(J) = (16.0D0/XCMP(J) - 5.0D0)/3.0D0
   60    CONTINUE
         CALL DWCS(N,TCMP,AK1CS,NTAK1,YCMP,B0,B1,B2)
         DO 70 J=1,N
            ZCMP(J) = EXP(-XCMP(J))*(1.250D0 + YCMP(J))/SQRT(XCMP(J))
   70    CONTINUE
         CALL DWSCTR(N,ZCMP,INDX,F)
      ENDIF
C
C  --------------------------
C  CASE  8.0 .LT. X .LE. XMAX
C  --------------------------
C
      CALL DWGTLE(M,X,8.0D0,XMAX,N,INDX)
      IF (N .GT. 0) THEN
         NTOT = NTOT + N
         CALL DWGTHR(N,X,INDX,XCMP)
         DO 80 J=1,N
            TCMP(J) = 16.0D0/XCMP(J) - 1.0D0
   80    CONTINUE
         CALL DWCS(N,TCMP,AK12CS,NTAK12,YCMP,B0,B1,B2)
         DO 90 J=1,N
            ZCMP(J) = EXP(-XCMP(J))*(1.250D0 + YCMP(J))/SQRT(XCMP(J))
   90    CONTINUE
         CALL DWSCTR(N,ZCMP,INDX,F)
      ENDIF
C
C  -----
C  EXITS
C  -----
C
C     ... NORMAL
C
      INFO = 0
      IF (NTOT .NE. M)  INFO = -1
      GO TO 999
C
C     ... M .LE. 0
C
  910 CONTINUE
      INFO = 1
      GO TO 999
C
C     ... X IS ZERO OR NEGATIVE
C
  920 CONTINUE
      INFO = 2
      INDX(1) = KEY
      GO TO 999
C
C     ... X SO SMALL K1 OVERFLOWS
C
  930 CONTINUE
      INFO = 3
      INDX(1) = KEY
      GO TO 999
C
  999 CONTINUE
      RETURN
      END
      SUBROUTINE DVJ0 (M, X, F, WORK, IWORK, INFO)
C***BEGIN PROLOGUE  DVJ0
C***PURPOSE  Computes the Bessel function of the first kind of order
C            zero (J0) for a vector of real arguments
C***LIBRARY   VFNLIB
C***CATEGORY  C10B1
C***TYPE      DOUBLE PRECISION (VJ0-S, DVJ0-D)
C***KEYWORDS  BESSEL FUNCTION,FIRST KIND, ORDER ZERO, SPECIAL FUNCTION,
C             VECTORIZED
C***AUTHOR  SAUNDERS, B. V., (NIST)
C             COMPUTING AND APPLIED MATHEMATICS LABORATORY
C             NATIONAL INSTITUTE OF STANDARDS AND TECHNOLOGY
C             GAITHERSBURG, MD 20899
C           BOISVERT, R. F., (NIST)
C             COMPUTING AND APPLIED MATHEMATICS LABORATORY
C             NATIONAL INSTITUTE OF STANDARDS AND TECHNOLOGY
C             GAITHERSBURG, MD 20899
C***DESCRIPTION
C
C   DVJ0 computes the Bessel function of the first kind of order zero
C   (J0) for real arguments using uniform approximation by Chebyshev
C   polynomials.
C
C
C   P A R A M E T E R S
C
C   M      (Input) Integer (M .GT. 0)
C          The number of arguments at which the function is to be
C          evaluated.
C
C   X      (Input) Double precision array of length M
C          The arguments at which the function is to be evaluated are
C          stored in X(1) to X(M) in any order.
C
C   F      (Output) Double precision array of length M
C          F(i) contains the value of the function at X(i), i=1,..,M.
C
C   WORK   (Work) Double precision vector of length 7*M
C
C   IWORK  (Work) Integer vector of length M
C
C   INFO   (Output) Integer
C          Indicates status of computed result. The following table
C          lists possible values and their meanings.  If OK=Yes then
C          all F(i) have been set, otherwise none have been set.
C
C          INFO  OK            Description
C          ------------------------------------------------------------
C            0   Yes  Successfull execution.
C            1   No   Error: M .LE. 0
C            2   No   Error: Some abs(X(i)) so big that no precision
C                     possible in computing J0.  The index of the
C                     first offending argument is returned in IWORK(1).
C
C *********************************************************************
C   This routine is a modification of the function DBESJ0 developed by
C   W. Fullerton of LANL.
C *********************************************************************
C
C***REFERENCES  (NONE)
C***ROUTINES CALLED  DWJ0
C***REVISION HISTORY  (YYMMDD)
C   910226  DATE WRITTEN
C***END PROLOGUE  DVJ0
C
C     ----------
C     PARAMETERS
C     ----------
C
      INTEGER INFO, IWORK, M
      DOUBLE PRECISION F, X, WORK
C
      DIMENSION X(M), F(M), WORK(7*M), IWORK(M)
C
C     ---------------
C     LOCAL VARIABLES
C     ---------------
C
      INTEGER IWB0, IWB1, IWB2, IWY, IWTC, IWYC, IWZC, JIN
C
C
C***FIRST EXECUTABLE STATEMENT  DVJ0
C
C     ... PARTITION WORK ARRAYS
C
      IWY   = 1
      IWTC  = IWY  + M
      IWYC  = IWTC + M
      IWZC  = IWYC + M
      IWB0  = IWZC + M
      IWB1  = IWB0 + M
      IWB2  = IWB1 + M
C     Total = IWB2 + M
C
      JIN   = 1
C     Total = JIN  + M
C
C     ... DWJ0 DOES ALL THE WORK
C
      CALL DWJ0(M,X,F,WORK(IWY),WORK(IWTC),WORK(IWYC),WORK(IWZC),
     +         IWORK(JIN),WORK(IWB0),WORK(IWB1),WORK(IWB2),INFO)
C
      RETURN
      END
      SUBROUTINE DWJ0 (M, X, F, Y, TCMP, YCMP, ZCMP, INDX, B0, B1, B2, 
     +   INFO)
C***BEGIN PROLOGUE  DWJ0
C***SUBSIDIARY
C***PURPOSE  Computes the Bessel function of the first kind
C            of order zero (J0) for a vector of arguments
C***LIBRARY   VFNLIB
C***AUTHOR  SAUNDERS, B. V., (NIST)
C             COMPUTING AND APPLIED MATHEMATICS LABORATORY
C             NATIONAL INSTITUTE OF STANDARDS AND TECHNOLOGY
C             GAITHERSBURG, MD 20899
C           BOISVERT, R. F., (NIST)
C             COMPUTING AND APPLIED MATHEMATICS LABORATORY
C             NATIONAL INSTITUTE OF STANDARDS AND TECHNOLOGY
C             GAITHERSBURG, MD 20899
C***ROUTINES CALLED  IDWCS, D1MACH, DWNGT, DWGTHR, DWGTLE, DWGT, DWLE,
C                    DWSCTR, DWCS
C***REVISION HISTORY  (YYMMDD)
C   910226  DATE WRITTEN
C***END PROLOGUE  DWJ0
C
C  ----------
C  PARAMETERS
C  ----------
C
      INTEGER INFO, INDX, M
      DOUBLE PRECISION B0, B1, B2, F, X, Y, TCMP, YCMP, ZCMP
C
      DIMENSION B0(M), B1(M), B2(M), F(M), INDX(M), X(M), Y(M),
     +          TCMP(M), YCMP(M), ZCMP(M)
C
C  ---------------
C  LOCAL VARIABLES
C  ---------------
C
      INTEGER LBJ0, LBM0, LBM02, LBTH0, LBT02
      PARAMETER (LBJ0=19, LBM0=37, LBM02=40, LBTH0=44, LBT02=39)
C
      INTEGER I, IDWCS, J, JH, KEY, N, NA, NB, NTJ0, NTM0, NTM02,
     +        NTTH0, NTT02
      DOUBLE PRECISION BJ0CS, BM0CS, BM02CS, BTH0CS, BT02CS, C1, C2,
     +        EPMACH, EPS, PI4, D1MACH, XSML, XMAX
C
      DIMENSION BJ0CS(LBJ0), BM0CS(LBM0), BM02CS(LBM02), BTH0CS(LBTH0),
     +          BT02CS(LBT02)
C
      SAVE BJ0CS, BM0CS, BM02CS, BTH0CS, BT02CS, NTJ0, NTM0, NTM02,
     +     NTTH0, NTT02, PI4, XSML, XMAX
C
C----------------------------------------------------------------------
C
C Series for BJ0        on the interval  0.          to  1.60000D+01
C                                        with weighted error   4.39E-32
C                                         log weighted error  31.36
C                               significant figures required  31.21
C                                    decimal places required  32.00
C
      DATA BJ0 CS(  1) / +.1002541619 6893913701 0731272640 74 D+0     /
      DATA BJ0 CS(  2) / -.6652230077 6440513177 6787578311 24 D+0     /
      DATA BJ0 CS(  3) / +.2489837034 9828131370 4604687266 80 D+0     /
      DATA BJ0 CS(  4) / -.3325272317 0035769653 8843415038 54 D-1     /
      DATA BJ0 CS(  5) / +.2311417930 4694015462 9049241177 29 D-2     /
      DATA BJ0 CS(  6) / -.9911277419 9508092339 0485193365 49 D-4     /
      DATA BJ0 CS(  7) / +.2891670864 3998808884 7339037470 78 D-5     /
      DATA BJ0 CS(  8) / -.6121085866 3032635057 8184074815 16 D-7     /
      DATA BJ0 CS(  9) / +.9838650793 8567841324 7687486364 15 D-9     /
      DATA BJ0 CS( 10) / -.1242355159 7301765145 5158970068 36 D-10    /
      DATA BJ0 CS( 11) / +.1265433630 2559045797 9158272103 63 D-12    /
      DATA BJ0 CS( 12) / -.1061945649 5287244546 9148175129 59 D-14    /
      DATA BJ0 CS( 13) / +.7470621075 8024567437 0989155840 00 D-17    /
      DATA BJ0 CS( 14) / -.4469703227 4412780547 6270079999 99 D-19    /
      DATA BJ0 CS( 15) / +.2302428158 4337436200 5230933333 33 D-21    /
      DATA BJ0 CS( 16) / -.1031914479 4166698148 5226666666 66 D-23    /
      DATA BJ0 CS( 17) / +.4060817827 4873322700 8000000000 00 D-26    /
      DATA BJ0 CS( 18) / -.1414383600 5240913919 9999999999 99 D-28    /
      DATA BJ0 CS( 19) / +.4391090549 6698880000 0000000000 00 D-31    /
C
C-------------------------------------------------------------------
C
C Series for BM0        on the interval  0.          to  6.25000D-02
C                                        with weighted error   4.40E-32
C                                         log weighted error  31.36
C                               significant figures required  30.02
C                                    decimal places required  32.14
C
      DATA BM0 CS(  1) / +.9211656246 8277427125 7376773018 2 D-1      /
      DATA BM0 CS(  2) / -.1050590997 2719051024 8071637175 5 D-2      /
      DATA BM0 CS(  3) / +.1470159840 7687597540 5639285095 2 D-4      /
      DATA BM0 CS(  4) / -.5058557606 0385542233 4792932770 2 D-6      /
      DATA BM0 CS(  5) / +.2787254538 6324441766 3035613788 1 D-7      /
      DATA BM0 CS(  6) / -.2062363611 7809148026 1884101897 3 D-8      /
      DATA BM0 CS(  7) / +.1870214313 1388796751 3817259626 1 D-9      /
      DATA BM0 CS(  8) / -.1969330971 1356362002 4173077782 5 D-10     /
      DATA BM0 CS(  9) / +.2325973793 9992754440 1250881805 2 D-11     /
      DATA BM0 CS( 10) / -.3009520344 9382502728 5122473448 2 D-12     /
      DATA BM0 CS( 11) / +.4194521333 8506691814 7120676864 6 D-13     /
      DATA BM0 CS( 12) / -.6219449312 1884458259 7326742956 4 D-14     /
      DATA BM0 CS( 13) / +.9718260411 3360684696 0176588526 9 D-15     /
      DATA BM0 CS( 14) / -.1588478585 7010752073 6663596693 7 D-15     /
      DATA BM0 CS( 15) / +.2700072193 6713088900 8621732445 8 D-16     /
      DATA BM0 CS( 16) / -.4750092365 2340089924 7750478677 3 D-17     /
      DATA BM0 CS( 17) / +.8615128162 6043708731 9170374656 0 D-18     /
      DATA BM0 CS( 18) / -.1605608686 9561448157 4560270335 9 D-18     /
      DATA BM0 CS( 19) / +.3066513987 3144829751 8853980159 9 D-19     /
      DATA BM0 CS( 20) / -.5987764223 1939564306 9650561706 6 D-20     /
      DATA BM0 CS( 21) / +.1192971253 7482483064 8906984106 6 D-20     /
      DATA BM0 CS( 22) / -.2420969142 0448054894 8468258133 3 D-21     /
      DATA BM0 CS( 23) / +.4996751760 5106164533 7100287999 9 D-22     /
      DATA BM0 CS( 24) / -.1047493639 3511585100 9504051199 9 D-22     /
      DATA BM0 CS( 25) / +.2227786843 7974681010 4818346666 6 D-23     /
      DATA BM0 CS( 26) / -.4801813239 3981628623 7054293333 3 D-24     /
      DATA BM0 CS( 27) / +.1047962723 4709599564 7699626666 6 D-24     /
      DATA BM0 CS( 28) / -.2313858165 6786153251 0126080000 0 D-25     /
      DATA BM0 CS( 29) / +.5164823088 4626742116 3519999999 9 D-26     /
      DATA BM0 CS( 30) / -.1164691191 8500653895 2540159999 9 D-26     /
      DATA BM0 CS( 31) / +.2651788486 0433192829 5833600000 0 D-27     /
      DATA BM0 CS( 32) / -.6092559503 8257284976 9130666666 6 D-28     /
      DATA BM0 CS( 33) / +.1411804686 1442593080 3882666666 6 D-28     /
      DATA BM0 CS( 34) / -.3298094961 2317372457 5061333333 3 D-29     /
      DATA BM0 CS( 35) / +.7763931143 0740650317 1413333333 3 D-30     /
      DATA BM0 CS( 36) / -.1841031343 6614584784 2133333333 3 D-30     /
      DATA BM0 CS( 37) / +.4395880138 5943107371 0079999999 9 D-31     /
C
C-------------------------------------------------------------------
C
C Series for BM02       Used in double precision version only
C                                        with weighted error   4.72E-32
C                                         log weighted error  31.33
C                               significant figures required  30.00
C                                    decimal places required  32.13
C
      DATA BM02CS(  1) / +.9500415145 2283813693 3086133556 0 D-1      /
      DATA BM02CS(  2) / -.3801864682 3656709917 4808156685 1 D-3      /
      DATA BM02CS(  3) / +.2258339301 0314811929 5182992722 4 D-5      /
      DATA BM02CS(  4) / -.3895725802 3722287647 3062141260 5 D-7      /
      DATA BM02CS(  5) / +.1246886416 5120816979 3099052972 5 D-8      /
      DATA BM02CS(  6) / -.6065949022 1025037798 0383505838 7 D-10     /
      DATA BM02CS(  7) / +.4008461651 4217469910 1527597104 5 D-11     /
      DATA BM02CS(  8) / -.3350998183 3980942184 6729879457 4 D-12     /
      DATA BM02CS(  9) / +.3377119716 5174173670 6326434199 6 D-13     /
      DATA BM02CS( 10) / -.3964585901 6350127005 6935629582 3 D-14     /
      DATA BM02CS( 11) / +.5286111503 8838572173 8793974473 5 D-15     /
      DATA BM02CS( 12) / -.7852519083 4508523136 5464024349 3 D-16     /
      DATA BM02CS( 13) / +.1280300573 3866822010 1163407344 9 D-16     /
      DATA BM02CS( 14) / -.2263996296 3914297762 8709924488 4 D-17     /
      DATA BM02CS( 15) / +.4300496929 6567903886 4641029047 7 D-18     /
      DATA BM02CS( 16) / -.8705749805 1325870797 4753545145 5 D-19     /
      DATA BM02CS( 17) / +.1865862713 9620951411 8144277205 0 D-19     /
      DATA BM02CS( 18) / -.4210482486 0930654573 4508697230 1 D-20     /
      DATA BM02CS( 19) / +.9956676964 2284009915 8162741784 2 D-21     /
      DATA BM02CS( 20) / -.2457357442 8053133596 0592147854 7 D-21     /
      DATA BM02CS( 21) / +.6307692160 7620315680 8735370705 9 D-22     /
      DATA BM02CS( 22) / -.1678773691 4407401426 9333117238 8 D-22     /
      DATA BM02CS( 23) / +.4620259064 6739044337 7087813608 7 D-23     /
      DATA BM02CS( 24) / -.1311782266 8603087322 3769340249 6 D-23     /
      DATA BM02CS( 25) / +.3834087564 1163028277 4792244027 6 D-24     /
      DATA BM02CS( 26) / -.1151459324 0777412710 7261329357 6 D-24     /
      DATA BM02CS( 27) / +.3547210007 5233385230 7697134521 3 D-25     /
      DATA BM02CS( 28) / -.1119218385 8150046462 6435594217 6 D-25     /
      DATA BM02CS( 29) / +.3611879427 6298378316 9840499425 7 D-26     /
      DATA BM02CS( 30) / -.1190687765 9133331500 9264176246 3 D-26     /
      DATA BM02CS( 31) / +.4005094059 4039681318 0247644953 6 D-27     /
      DATA BM02CS( 32) / -.1373169422 4522123905 9519391601 7 D-27     /
      DATA BM02CS( 33) / +.4794199088 7425315859 9649152643 7 D-28     /
      DATA BM02CS( 34) / -.1702965627 6241095840 0699447645 2 D-28     /
      DATA BM02CS( 35) / +.6149512428 9363300715 0357516132 4 D-29     /
      DATA BM02CS( 36) / -.2255766896 5818283499 4430023724 2 D-29     /
      DATA BM02CS( 37) / +.8399707509 2942994860 6165835320 0 D-30     /
      DATA BM02CS( 38) / -.3172997595 5626023555 6742393615 2 D-30     /
      DATA BM02CS( 39) / +.1215205298 8812985545 8333302651 4 D-30     /
      DATA BM02CS( 40) / -.4715852749 7544386930 1321056804 5 D-31     /
C
C-------------------------------------------------------------------
C
C Series for BTH0       on the interval  0.          to  6.25000D-02
C                                        with weighted error   2.66E-32
C                                         log weighted error  31.57
C                               significant figures required  30.67
C                                    decimal places required  32.40
C
      DATA BTH0CS(  1) / -.2490178086 2128936717 7097937899 67 D+0     /
      DATA BTH0CS(  2) / +.4855029960 9623749241 0486155354 85 D-3     /
      DATA BTH0CS(  3) / -.5451183734 5017204950 6562735635 05 D-5     /
      DATA BTH0CS(  4) / +.1355867305 9405964054 3774459299 03 D-6     /
      DATA BTH0CS(  5) / -.5569139890 2227626227 5832184149 20 D-8     /
      DATA BTH0CS(  6) / +.3260903182 4994335304 0042057194 68 D-9     /
      DATA BTH0CS(  7) / -.2491880786 2461341125 2379038779 93 D-10    /
      DATA BTH0CS(  8) / +.2344937742 0882520554 3524135648 91 D-11    /
      DATA BTH0CS(  9) / -.2609653444 4310387762 1775747661 36 D-12    /
      DATA BTH0CS( 10) / +.3335314042 0097395105 8699550149 23 D-13    /
      DATA BTH0CS( 11) / -.4789000044 0572684646 7507705574 09 D-14    /
      DATA BTH0CS( 12) / +.7595617843 6192215972 6425685452 48 D-15    /
      DATA BTH0CS( 13) / -.1313155601 6891440382 7733974876 33 D-15    /
      DATA BTH0CS( 14) / +.2448361834 5240857495 4268207383 55 D-16    /
      DATA BTH0CS( 15) / -.4880572981 0618777683 2567619183 31 D-17    /
      DATA BTH0CS( 16) / +.1032728502 9786316149 2237563612 04 D-17    /
      DATA BTH0CS( 17) / -.2305763381 5057217157 0047445270 25 D-18    /
      DATA BTH0CS( 18) / +.5404444300 1892693993 0171084837 65 D-19    /
      DATA BTH0CS( 19) / -.1324069519 4366572724 1550328823 85 D-19    /
      DATA BTH0CS( 20) / +.3378079562 1371970203 4247921247 22 D-20    /
      DATA BTH0CS( 21) / -.8945762915 7111779003 0269262922 99 D-21    /
      DATA BTH0CS( 22) / +.2451990688 9219317090 8999086514 05 D-21    /
      DATA BTH0CS( 23) / -.6938842287 6866318680 1399331576 57 D-22    /
      DATA BTH0CS( 24) / +.2022827871 4890138392 9463033377 91 D-22    /
      DATA BTH0CS( 25) / -.6062850000 2335483105 7941953717 64 D-23    /
      DATA BTH0CS( 26) / +.1864974896 4037635381 8237883962 70 D-23    /
      DATA BTH0CS( 27) / -.5878373238 4849894560 2450365308 67 D-24    /
      DATA BTH0CS( 28) / +.1895859144 7999563485 5311795035 13 D-24    /
      DATA BTH0CS( 29) / -.6248197937 2258858959 2916207285 65 D-25    /
      DATA BTH0CS( 30) / +.2101790168 4551024686 6386335290 74 D-25    /
      DATA BTH0CS( 31) / -.7208430093 5209253690 8139339924 46 D-26    /
      DATA BTH0CS( 32) / +.2518136389 2474240867 1564059767 46 D-26    /
      DATA BTH0CS( 33) / -.8951804225 8785778806 1439459536 43 D-27    /
      DATA BTH0CS( 34) / +.3235723747 9762298533 2562358685 87 D-27    /
      DATA BTH0CS( 35) / -.1188301051 9855353657 0471441137 96 D-27    /
      DATA BTH0CS( 36) / +.4430628690 7358104820 5792319417 31 D-28    /
      DATA BTH0CS( 37) / -.1676100964 8834829495 7920101356 81 D-28    /
      DATA BTH0CS( 38) / +.6429294692 1207466972 5323939660 88 D-29    /
      DATA BTH0CS( 39) / -.2499226116 6978652421 2072136827 63 D-29    /
      DATA BTH0CS( 40) / +.9839979429 9521955672 8282603553 18 D-30    /
      DATA BTH0CS( 41) / -.3922037524 2408016397 9891316261 58 D-30    /
      DATA BTH0CS( 42) / +.1581810703 0056522138 5906188456 92 D-30    /
      DATA BTH0CS( 43) / -.6452550614 4890715944 3440983654 26 D-31    /
      DATA BTH0CS( 44) / +.2661111136 9199356137 1770183463 67 D-31    /
C
C-------------------------------------------------------------------
C
C Series for BT02       Used in double precision version only
C                                        with weighted error   2.99E-32
C                                         log weighted error  31.52
C
      DATA BT02CS(  1) / -.2454829521 3424597462 0504672493 24 D+0     /
      DATA BT02CS(  2) / +.1254412103 9084615780 7853317782 99 D-2     /
      DATA BT02CS(  3) / -.3125395041 4871522854 9734467095 71 D-4     /
      DATA BT02CS(  4) / +.1470977824 9940831164 4534269693 14 D-5     /
      DATA BT02CS(  5) / -.9954348893 7950033643 4688503511 58 D-7     /
      DATA BT02CS(  6) / +.8549316673 3203041247 5787113977 51 D-8     /
      DATA BT02CS(  7) / -.8698975952 6554334557 9855121791 92 D-9     /
      DATA BT02CS(  8) / +.1005209953 3559791084 5401010821 53 D-9     /
      DATA BT02CS(  9) / -.1282823060 1708892903 4836236855 44 D-10    /
      DATA BT02CS( 10) / +.1773170078 1805131705 6557504510 23 D-11    /
      DATA BT02CS( 11) / -.2617457456 9485577488 6362841809 25 D-12    /
      DATA BT02CS( 12) / +.4082835138 9972059621 9664812211 03 D-13    /
      DATA BT02CS( 13) / -.6675166823 9742720054 6067495542 61 D-14    /
      DATA BT02CS( 14) / +.1136576139 3071629448 3924695499 51 D-14    /
      DATA BT02CS( 15) / -.2005118962 0647160250 5592664121 17 D-15    /
      DATA BT02CS( 16) / +.3649797879 4766269635 7205914641 06 D-16    /
      DATA BT02CS( 17) / -.6830963756 4582303169 3558437888 00 D-17    /
      DATA BT02CS( 18) / +.1310758314 5670756620 0571042679 46 D-17    /
      DATA BT02CS( 19) / -.2572336310 1850607778 7571306495 99 D-18    /
      DATA BT02CS( 20) / +.5152165744 1863959925 2677809493 33 D-19    /
      DATA BT02CS( 21) / -.1051301756 3758802637 9407414613 33 D-19    /
      DATA BT02CS( 22) / +.2182038199 1194813847 3010845013 33 D-20    /
      DATA BT02CS( 23) / -.4600470121 0362160577 2259054933 33 D-21    /
      DATA BT02CS( 24) / +.9840700692 5466818520 9536511999 99 D-22    /
      DATA BT02CS( 25) / -.2133403803 5728375844 7359863466 66 D-22    /
      DATA BT02CS( 26) / +.4683103642 3973365296 0662869333 33 D-23    /
      DATA BT02CS( 27) / -.1040021369 1985747236 5133823999 99 D-23    /
      DATA BT02CS( 28) / +.2334910567 7301510051 7777408000 00 D-24    /
      DATA BT02CS( 29) / -.5295682532 3318615788 0497493333 33 D-25    /
      DATA BT02CS( 30) / +.1212634195 2959756829 1962879999 99 D-25    /
      DATA BT02CS( 31) / -.2801889708 2289428760 2756266666 66 D-26    /
      DATA BT02CS( 32) / +.6529267898 7012873342 5937066666 66 D-27    /
      DATA BT02CS( 33) / -.1533798006 1873346427 8357333333 33 D-27    /
      DATA BT02CS( 34) / +.3630588430 6364536682 3594666666 66 D-28    /
      DATA BT02CS( 35) / -.8656075571 3629122479 1722666666 66 D-29    /
      DATA BT02CS( 36) / +.2077990997 2536284571 2383999999 99 D-29    /
      DATA BT02CS( 37) / -.5021117022 1417221674 3253333333 33 D-30    /
      DATA BT02CS( 38) / +.1220836027 9441714184 1919999999 99 D-30    /
      DATA BT02CS( 39) / -.2986005626 7039913454 2506666666 66 D-31    /
C
C-------------------------------------------------------------------
C
      DATA NTJ0 / 0 /
C
C***FIRST EXECUTABLE STATEMENT  DWJ0
C
      IF (M .LE. 0) GO TO 910
C
      IF (NTJ0.EQ.0) THEN
         PI4 = ATAN(1.0D0)
         EPMACH = D1MACH(3)
         EPS = 0.10D0*EPMACH
         NTJ0 = IDWCS(BJ0CS, LBJ0, EPS)
         NTM0 = IDWCS(BM0CS, LBM0, EPS)
         NTM02 = IDWCS(BM02CS, LBM02, EPS)
         NTTH0 = IDWCS(BTH0CS, LBTH0, EPS)
         NTT02 = IDWCS(BT02CS, LBT02, EPS)
         XSML = SQRT(4.0D0*EPMACH)
         XMAX = 1.0D0/D1MACH(4)
      ENDIF
C
      DO 10 I=1,M
         Y(I) = ABS(X(I))
  10  CONTINUE
C
      CALL DWNGT(M,Y,XMAX,KEY)
      IF (KEY .NE. 0) GO TO 920
C
C  ----------------
C  CASE Y .LE. XSML
C  ----------------
C
      DO 15 I=1,M
         F(I) = 1.0D0
  15  CONTINUE
C
C  --------------------------
C  CASE  XSML .LT. Y .LE. 4.0
C  --------------------------
C
      CALL DWGTLE(M,Y,XSML,4.0D0,N,INDX)
      IF (N .GT. 0) THEN
         CALL DWGTHR(N,Y,INDX,YCMP)
         DO 20 J=1,N
            TCMP(J) = 0.1250D0*YCMP(J)**2 - 1.0D0
  20     CONTINUE
         CALL DWCS(N,TCMP,BJ0CS,NTJ0,ZCMP,B0,B1,B2)
         CALL DWSCTR(N,ZCMP,INDX,F)
      ENDIF
C
C  ----------------
C  CASE  Y .GT. 4.0
C  ----------------
C
      CALL DWGTLE(M,Y,4.0D0,8.0D0,NA,INDX)
      JH = NA + 1
      CALL DWGT(M,Y,8.0D0,NB,INDX(JH))
      N = NA + NB
      IF (N .GT. 0) THEN
         CALL DWGTHR(N,Y,INDX,YCMP)
         C1 = 128.0D0/3.0D0
         C2 = 5.0D0/3.0D0
         DO 50 J=1,NA
            ZCMP(J) = C1/YCMP(J)**2 - C2
  50     CONTINUE
         CALL DWCS(NA,ZCMP,BM0CS,NTM0,Y,B0,B1,B2)
         CALL DWCS(NA,ZCMP,BT02CS,NTT02,ZCMP,B0,B1,B2)
         DO 60 J=JH,N
            ZCMP(J) = 128.0D0/YCMP(J)**2 - 1.0D0
  60     CONTINUE
         CALL DWCS(NB,ZCMP(JH),BM02CS,NTM02,Y(JH),B0,B1,B2)
         CALL DWCS(NB,ZCMP(JH),BTH0CS,NTTH0,ZCMP(JH),B0,B1,B2)
         DO 70 J=1,N
            Y(J) = (0.750D0 + Y(J)) / SQRT(YCMP(J))
            ZCMP(J) = (YCMP(J) - PI4) + ZCMP(J) / YCMP(J)
            ZCMP(J) = Y(J) * COS(ZCMP(J))
  70     CONTINUE
         CALL DWSCTR(N,ZCMP,INDX,F)
      ENDIF
C
C  -----
C  EXITS
C  -----
C
C     ... NORMAL
C
      INFO = 0
      GO TO 999
C
C     ... M .LE. 0
C
  910 CONTINUE
      INFO = 1
      GO TO 999
C
C     ... NO PRECISION BECAUSE X BIG
C
  920 CONTINUE
      INFO = 2
      INDX(1) = KEY
      GO TO 999
C
  999 CONTINUE
      RETURN
      END
      SUBROUTINE DVJ1 (M, X, F, WORK, IWORK, INFO)
C***BEGIN PROLOGUE  DVJ1
C***PURPOSE  Computes the Bessel function of the first kind
C            of order one (J1) for a vector of real arguments
C***LIBRARY   VFNLIB
C***CATEGORY  C10B1
C***TYPE      DOUBLE PRECISION (VJ1-S, DVJ1-D)
C***KEYWORDS  BESSEL FUNCTION, FIRST KIND, SPECIAL FUNCTION,
C             ORDER ONE, VECTORIZED
C***AUTHOR  SAUNDERS, B. V., (NIST)
C             COMPUTING AND APPLIED MATHEMATICS LABORATORY
C             NATIONAL INSTITUTE OF STANDARDS AND TECHNOLOGY
C             GAITHERSBURG, MD 20899
C           BOISVERT, R. F., (NIST)
C             COMPUTING AND APPLIED MATHEMATICS LABORATORY
C             NATIONAL INSTITUTE OF STANDARDS AND TECHNOLOGY
C             GAITHERSBURG, MD 20899
C***DESCRIPTION
C
C   DVJ1 computes the  Bessel function of the first kind of order
C   one (J1) for real arguments using uniform approximation by
C   Chebyshev polynomials.
C
C
C   P A R A M E T E R S
C
C   M      (Input) Integer (M .GT. 0)
C          The number of arguments at which the function is to be
C          evaluated.
C
C   X      (Input) Double precision array of length M
C          The arguments at which the function is to be evaluated are
C          stored in X(1) to X(M) in any order.
C
C   F      (Output) Double precision array of length M
C          F(i) contains the value of the function at X(i), i=1,..,M.
C
C   WORK   (Work) Double precision vector of length 7*M
C
C   IWORK  (Work) Integer vector of length M
C
C
C   INFO   (Output) Integer
C          Indicates status of computed result. The following table
C          lists possible values and their meanings.  If OK=Yes then
C          all F(i) have been set, otherwise none have been set.
C
C          INFO  OK            Description
C          ------------------------------------------------------------
C           -1   Yes  Warning: Some abs(X(i)) so small J1 underflows.
C                     The corresponding F(i) are set to zero.
C            0   Yes  Successfull execution.
C            1   No   Error: M .LE. 0
C            2   No   Error: Some abs(X(i)) so big that no precision
C                     possible in computing J1.  The index of the
C                     first offending argument is returned in IWORK(1).
C
C *********************************************************************
C   This routine is a modification of the function DBESJ1 developed by
C   W. Fullerton of LANL.
C *********************************************************************
C
C***REFERENCES  (NONE)
C***ROUTINES CALLED  DWJ1
C***REVISION HISTORY  (YYMMDD)
C   910226  DATE WRITTEN
C***END PROLOGUE  DVJ1
C
C     ----------
C     PARAMETERS
C     ----------
C
      INTEGER INFO, IWORK, M
      DOUBLE PRECISION F, X, WORK
C
      DIMENSION X(M), F(M), WORK(7*M), IWORK(M)
C
C     ---------------
C     LOCAL VARIABLES
C     ---------------
C
      INTEGER IWB0, IWB1, IWB2, IWY, IWTC, IWYC, IWZC, JIN
C
C***FIRST EXECUTABLE STATEMENT  DVJ1
C
C     ... PARTITION WORK ARRAYS
C
      IWY   = 1
      IWTC  = IWY  + M
      IWYC  = IWTC + M
      IWZC  = IWYC + M
      IWB0  = IWZC + M
      IWB1  = IWB0 + M
      IWB2  = IWB1 + M
C     Total = IWB2 + M
C
      JIN   = 1
C     Total = JIN  + M
C
C     ... DWJ1 DOES ALL THE WORK
C
      CALL DWJ1(M,X,F,WORK(IWY),WORK(IWTC),WORK(IWYC),WORK(IWZC),
     +         IWORK(JIN),WORK(IWB0),WORK(IWB1),WORK(IWB2),INFO)
C
      RETURN
      END
      SUBROUTINE DWJ1(M, X, F, Y, TCMP, YCMP, ZCMP, INDX, B0, B1, B2, 
     +   INFO)
C***BEGIN PROLOGUE  DWJ1
C***SUBSIDIARY
C***PURPOSE  Computes the Bessel function of the first kind
C            of order one (J1) for a vector of arguments
C***LIBRARY   VFNLIB
C***AUTHOR  SAUNDERS, B. V., (NIST)
C             COMPUTING AND APPLIED MATHEMATICS LABORATORY
C             NATIONAL INSTITUTE OF STANDARDS AND TECHNOLOGY
C             GAITHERSBURG, MD 20899
C           BOISVERT, R. F., (NIST)
C             COMPUTING AND APPLIED MATHEMATICS LABORATORY
C             NATIONAL INSTITUTE OF STANDARDS AND TECHNOLOGY
C             GAITHERSBURG, MD 20899
C***ROUTINES CALLED  IDWCS, D1MACH, DWNGT, DWGTHR, DWGTLE, DWGT,
C                    DWLE, DWSCTR, DWCS
C***REVISION HISTORY  (YYMMDD)
C   910226  DATE WRITTEN
C***END PROLOGUE  DWJ1
C
C  ----------
C  PARAMETERS
C  ----------
C
      INTEGER INFO, INDX, M
      DOUBLE PRECISION B0, B1, B2, F, X, Y, TCMP, YCMP, ZCMP
C
      DIMENSION B0(M), B1(M), B2(M), F(M), INDX(M), X(M), Y(M),
     +          TCMP(M), YCMP(M), ZCMP(M)
C
C  ---------------
C  LOCAL VARIABLES
C  ---------------
C
      INTEGER LBJ1, LBM1, LBM12, LBTH1, LBT12
      PARAMETER (LBJ1=19, LBM1=37, LBM12=40, LBTH1=44, LBT12=39)
C
      INTEGER I, IDWCS, J, JH, K, KEY, N, NA, NB, NTJ1, NTM1, NTM12,
     +        NTTH1, NTT12, NTOT
      DOUBLE PRECISION BJ1CS, BM1CS, BM12CS, BTH1CS, BT12CS, C1, C2,
     +        EPMACH, EPS, PI4, D1MACH, TPI4, XMAX, XMIN, XSML
C
      DIMENSION BJ1CS(LBJ1), BM1CS(LBM1), BM12CS(LBM12), BTH1CS(LBTH1),
     +          BT12CS(LBT12)
C
      SAVE BJ1CS, BM1CS, BM12CS, BTH1CS, BT12CS, NTJ1, NTM1, NTM12,
     +      NTTH1, NTT12, PI4, TPI4, XMIN, XSML, XMAX
C
C----------------------------------------------------------------------
C
C Series for BJ1        on the interval  0.          to  1.60000D+01
C                                        with weighted error   1.16E-33
C                                         log weighted error  32.93
C                               significant figures required  32.36
C                                    decimal places required  33.57
C
      DATA BJ1 CS(  1) / -.1172614151 3332786560 6240574524 003 D+0    /
      DATA BJ1 CS(  2) / -.2536152183 0790639562 3030884554 698 D+0    /
      DATA BJ1 CS(  3) / +.5012708098 4469568505 3656363203 743 D-1    /
      DATA BJ1 CS(  4) / -.4631514809 6250819184 2619728789 772 D-2    /
      DATA BJ1 CS(  5) / +.2479962294 1591402453 9124064592 364 D-3    /
      DATA BJ1 CS(  6) / -.8678948686 2788258452 1246435176 416 D-5    /
      DATA BJ1 CS(  7) / +.2142939171 4379369150 2766250991 292 D-6    /
      DATA BJ1 CS(  8) / -.3936093079 1831797922 9322764073 061 D-8    /
      DATA BJ1 CS(  9) / +.5591182317 9468800401 8248059864 032 D-10   /
      DATA BJ1 CS( 10) / -.6327616404 6613930247 7695274014 880 D-12   /
      DATA BJ1 CS( 11) / +.5840991610 8572470032 6945563268 266 D-14   /
      DATA BJ1 CS( 12) / -.4482533818 7012581903 9135059199 999 D-16   /
      DATA BJ1 CS( 13) / +.2905384492 6250246630 6018688000 000 D-18   /
      DATA BJ1 CS( 14) / -.1611732197 8414416541 2118186666 666 D-20   /
      DATA BJ1 CS( 15) / +.7739478819 3927463729 8346666666 666 D-23   /
      DATA BJ1 CS( 16) / -.3248693782 1119984114 3466666666 666 D-25   /
      DATA BJ1 CS( 17) / +.1202237677 2274102272 0000000000 000 D-27   /
      DATA BJ1 CS( 18) / -.3952012212 6513493333 3333333333 333 D-30   /
      DATA BJ1 CS( 19) / +.1161678082 2664533333 3333333333 333 D-32   /
C
C-------------------------------------------------------------------
C
C Series for BM1        on the interval  0.          to  6.25000D-02
C                                        with weighted error   4.91E-32
C                                         log weighted error  31.31
C                               significant figures required  30.04
C                                    decimal places required  32.09
C
      DATA BM1 CS(  1) / +.1069845452 6180630149 6998530853 8 D+0      /
      DATA BM1 CS(  2) / +.3274915039 7159649007 2905514344 5 D-2      /
      DATA BM1 CS(  3) / -.2987783266 8316985920 3044577793 8 D-4      /
      DATA BM1 CS(  4) / +.8331237177 9919745313 9322266902 3 D-6      /
      DATA BM1 CS(  5) / -.4112665690 3020073048 9638172549 8 D-7      /
      DATA BM1 CS(  6) / +.2855344228 7892152207 1975766316 1 D-8      /
      DATA BM1 CS(  7) / -.2485408305 4156238780 6002659605 5 D-9      /
      DATA BM1 CS(  8) / +.2543393338 0725824427 4248439717 4 D-10     /
      DATA BM1 CS(  9) / -.2941045772 8229675234 8975082790 9 D-11     /
      DATA BM1 CS( 10) / +.3743392025 4939033092 6505615362 6 D-12     /
      DATA BM1 CS( 11) / -.5149118293 8211672187 2054824352 7 D-13     /
      DATA BM1 CS( 12) / +.7552535949 8651439080 3404076419 9 D-14     /
      DATA BM1 CS( 13) / -.1169409706 8288464441 6629062246 4 D-14     /
      DATA BM1 CS( 14) / +.1896562449 4347915717 2182460506 0 D-15     /
      DATA BM1 CS( 15) / -.3201955368 6932864206 6477531639 4 D-16     /
      DATA BM1 CS( 16) / +.5599548399 3162041144 8416990549 3 D-17     /
      DATA BM1 CS( 17) / -.1010215894 7304324431 1939044454 4 D-17     /
      DATA BM1 CS( 18) / +.1873844985 7275629833 0204271957 3 D-18     /
      DATA BM1 CS( 19) / -.3563537470 3285802192 7430143999 9 D-19     /
      DATA BM1 CS( 20) / +.6931283819 9712383304 2276351999 9 D-20     /
      DATA BM1 CS( 21) / -.1376059453 4065001522 5140893013 3 D-20     /
      DATA BM1 CS( 22) / +.2783430784 1070802205 9977932799 9 D-21     /
      DATA BM1 CS( 23) / -.5727595364 3205616893 4866943999 9 D-22     /
      DATA BM1 CS( 24) / +.1197361445 9188926725 3575679999 9 D-22     /
      DATA BM1 CS( 25) / -.2539928509 8918719766 4144042666 6 D-23     /
      DATA BM1 CS( 26) / +.5461378289 6572959730 6961919999 9 D-24     /
      DATA BM1 CS( 27) / -.1189211341 7733202889 8628949333 3 D-24     /
      DATA BM1 CS( 28) / +.2620150977 3400815949 5782400000 0 D-25     /
      DATA BM1 CS( 29) / -.5836810774 2556859019 2093866666 6 D-26     /
      DATA BM1 CS( 30) / +.1313743500 0805957734 2361599999 9 D-26     /
      DATA BM1 CS( 31) / -.2985814622 5103803553 3277866666 6 D-27     /
      DATA BM1 CS( 32) / +.6848390471 3346049376 2559999999 9 D-28     /
      DATA BM1 CS( 33) / -.1584401568 2224767211 9296000000 0 D-28     /
      DATA BM1 CS( 34) / +.3695641006 5709380543 0101333333 3 D-29     /
      DATA BM1 CS( 35) / -.8687115921 1446682430 1226666666 6 D-30     /
      DATA BM1 CS( 36) / +.2057080846 1587634629 2906666666 6 D-30     /
      DATA BM1 CS( 37) / -.4905225761 1162255185 2373333333 3 D-31     /
C
C-------------------------------------------------------------------
C
C Series for BM12       Used in double precision version only
C                                        with weighted error   5.01E-32
C                                         log weighted error  31.30
C                               significant figures required  29.99
C                                    decimal places required  32.10
C
      DATA BM12CS(  1) / +.9807979156 2330500272 7209354693 7 D-1      /
      DATA BM12CS(  2) / +.1150961189 5046853061 7548348460 2 D-2      /
      DATA BM12CS(  3) / -.4312482164 3382054098 8935809773 2 D-5      /
      DATA BM12CS(  4) / +.5951839610 0888163078 1302980183 2 D-7      /
      DATA BM12CS(  5) / -.1704844019 8269098574 0070158647 8 D-8      /
      DATA BM12CS(  6) / +.7798265413 6111095086 5817382740 1 D-10     /
      DATA BM12CS(  7) / -.4958986126 7664158094 9175495186 5 D-11     /
      DATA BM12CS(  8) / +.4038432416 4211415168 3820226514 4 D-12     /
      DATA BM12CS(  9) / -.3993046163 7251754457 6548384664 5 D-13     /
      DATA BM12CS( 10) / +.4619886183 1189664943 1334243277 5 D-14     /
      DATA BM12CS( 11) / -.6089208019 0953833013 4547261933 3 D-15     /
      DATA BM12CS( 12) / +.8960930916 4338764821 5704804124 9 D-16     /
      DATA BM12CS( 13) / -.1449629423 9420231229 1651891892 5 D-16     /
      DATA BM12CS( 14) / +.2546463158 5377760561 6514964806 8 D-17     /
      DATA BM12CS( 15) / -.4809472874 6478364442 5926371862 0 D-18     /
      DATA BM12CS( 16) / +.9687684668 2925990490 8727583912 4 D-19     /
      DATA BM12CS( 17) / -.2067213372 2779660232 4503811755 1 D-19     /
      DATA BM12CS( 18) / +.4646651559 1503847318 0276780959 0 D-20     /
      DATA BM12CS( 19) / -.1094966128 8483341382 4135132833 9 D-20     /
      DATA BM12CS( 20) / +.2693892797 2886828609 0570761278 5 D-21     /
      DATA BM12CS( 21) / -.6894992910 9303744778 1897002685 7 D-22     /
      DATA BM12CS( 22) / +.1830268262 7520629098 9066855474 0 D-22     /
      DATA BM12CS( 23) / -.5025064246 3519164281 5611355322 4 D-23     /
      DATA BM12CS( 24) / +.1423545194 4548060396 3169363419 4 D-23     /
      DATA BM12CS( 25) / -.4152191203 6164503880 6888676980 1 D-24     /
      DATA BM12CS( 26) / +.1244609201 5039793258 8233007654 7 D-24     /
      DATA BM12CS( 27) / -.3827336370 5693042994 3191866128 6 D-25     /
      DATA BM12CS( 28) / +.1205591357 8156175353 7472398183 5 D-25     /
      DATA BM12CS( 29) / -.3884536246 3764880764 3185936112 4 D-26     /
      DATA BM12CS( 30) / +.1278689528 7204097219 0489528346 1 D-26     /
      DATA BM12CS( 31) / -.4295146689 4479462720 6193691591 2 D-27     /
      DATA BM12CS( 32) / +.1470689117 8290708864 5680270798 3 D-27     /
      DATA BM12CS( 33) / -.5128315665 1060731281 8037401779 6 D-28     /
      DATA BM12CS( 34) / +.1819509585 4711693854 8143737328 6 D-28     /
      DATA BM12CS( 35) / -.6563031314 8419808676 1863505037 3 D-29     /
      DATA BM12CS( 36) / +.2404898976 9199606531 9891487583 4 D-29     /
      DATA BM12CS( 37) / -.8945966744 6906124732 3495824297 9 D-30     /
      DATA BM12CS( 38) / +.3376085160 6572310266 3714897824 0 D-30     /
      DATA BM12CS( 39) / -.1291791454 6206563609 1309991696 6 D-30     /
      DATA BM12CS( 40) / +.5008634462 9588105206 8495150125 4 D-31     /
C
C-------------------------------------------------------------------
C
C Series for BTH1       on the interval  0.          to  6.25000D-02
C                                        with weighted error   2.82E-32
C                                         log weighted error  31.55
C                               significant figures required  31.12
C
      DATA BTH1CS(  1) / +.7474995720 3587276055 4434839696 95 D+0     /
      DATA BTH1CS(  2) / -.1240077714 4651711252 5457775413 84 D-2     /
      DATA BTH1CS(  3) / +.9925244240 4424527376 6414976895 92 D-5     /
      DATA BTH1CS(  4) / -.2030369073 7159711052 4193753756 08 D-6     /
      DATA BTH1CS(  5) / +.7535961770 5690885712 1840175836 29 D-8     /
      DATA BTH1CS(  6) / -.4166161271 5343550107 6300238562 28 D-9     /
      DATA BTH1CS(  7) / +.3070161807 0834890481 2451020912 16 D-10    /
      DATA BTH1CS(  8) / -.2817849963 7605213992 3240088839 24 D-11    /
      DATA BTH1CS(  9) / +.3079069673 9040295476 0281468216 47 D-12    /
      DATA BTH1CS( 10) / -.3880330026 2803434112 7873475547 81 D-13    /
      DATA BTH1CS( 11) / +.5509603960 8630904934 5617262085 62 D-14    /
      DATA BTH1CS( 12) / -.8659006076 8383779940 1033989539 94 D-15    /
      DATA BTH1CS( 13) / +.1485604914 1536749003 4236890606 83 D-15    /
      DATA BTH1CS( 14) / -.2751952981 5904085805 3712121250 09 D-16    /
      DATA BTH1CS( 15) / +.5455079609 0481089625 0362236409 23 D-17    /
      DATA BTH1CS( 16) / -.1148653450 1983642749 5436310271 77 D-17    /
      DATA BTH1CS( 17) / +.2553521337 7973900223 1990525335 22 D-18    /
      DATA BTH1CS( 18) / -.5962149019 7413450395 7682879078 49 D-19    /
      DATA BTH1CS( 19) / +.1455662290 2372718620 2883020058 33 D-19    /
      DATA BTH1CS( 20) / -.3702218542 2450538201 5797760195 93 D-20    /
      DATA BTH1CS( 21) / +.9776307412 5345357664 1684345179 24 D-21    /
      DATA BTH1CS( 22) / -.2672682163 9668488468 7237753930 52 D-21    /
      DATA BTH1CS( 23) / +.7545330038 4983271794 0381906557 64 D-22    /
      DATA BTH1CS( 24) / -.2194789991 9802744897 8923833716 47 D-22    /
      DATA BTH1CS( 25) / +.6564839462 3955262178 9069998174 93 D-23    /
      DATA BTH1CS( 26) / -.2015560429 8370207570 7840768695 19 D-23    /
      DATA BTH1CS( 27) / +.6341776855 6776143492 1446671856 70 D-24    /
      DATA BTH1CS( 28) / -.2041927788 5337895634 8137699555 91 D-24    /
      DATA BTH1CS( 29) / +.6719146422 0720567486 6589800185 51 D-25    /
      DATA BTH1CS( 30) / -.2256907911 0207573595 7090036873 36 D-25    /
      DATA BTH1CS( 31) / +.7729771989 2989706370 9269598719 29 D-26    /
      DATA BTH1CS( 32) / -.2696744451 2294640913 2114240809 20 D-26    /
      DATA BTH1CS( 33) / +.9574934451 8502698072 2955219336 27 D-27    /
      DATA BTH1CS( 34) / -.3456916844 8890113000 1756808276 27 D-27    /
      DATA BTH1CS( 35) / +.1268123481 7398436504 2119862383 74 D-27    /
      DATA BTH1CS( 36) / -.4723253663 0722639860 4649937134 45 D-28    /
      DATA BTH1CS( 37) / +.1785000847 8186376177 8586197964 17 D-28    /
      DATA BTH1CS( 38) / -.6840436100 4510395406 2152235667 46 D-29    /
      DATA BTH1CS( 39) / +.2656602867 1720419358 2934226722 12 D-29    /
      DATA BTH1CS( 40) / -.1045040252 7914452917 7141614846 70 D-29    /
      DATA BTH1CS( 41) / +.4161829082 5377144306 8619171970 64 D-30    /
      DATA BTH1CS( 42) / -.1677163920 3643714856 5013478828 87 D-30    /
      DATA BTH1CS( 43) / +.6836199777 6664389173 5359280285 28 D-31    /
      DATA BTH1CS( 44) / -.2817224786 1233641166 7395746228 10 D-31    /
C
C-------------------------------------------------------------------
C
C Series for BT12       Used in double precision version only
C                                        with weighted error   3.33E-32
C                                         log weighted error  31.48
C                               significant figures required  31.05
C                                    decimal places required  32.27
C
      DATA BT12CS(  1) / +.7382386012 8742974662 6208397927 64 D+0     /
      DATA BT12CS(  2) / -.3336111317 4483906384 4701476811 89 D-2     /
      DATA BT12CS(  3) / +.6146345488 8046964698 5148994201 86 D-4     /
      DATA BT12CS(  4) / -.2402458516 1602374264 9776354695 68 D-5     /
      DATA BT12CS(  5) / +.1466355557 7509746153 2105919972 04 D-6     /
      DATA BT12CS(  6) / -.1184191730 5589180567 0051475049 83 D-7     /
      DATA BT12CS(  7) / +.1157419896 3919197052 1254663030 55 D-8     /
      DATA BT12CS(  8) / -.1300116112 9439187449 3660077945 71 D-9     /
      DATA BT12CS(  9) / +.1624539114 1361731937 7421662736 67 D-10    /
      DATA BT12CS( 10) / -.2208963682 1403188752 1554417701 28 D-11    /
      DATA BT12CS( 11) / +.3218030425 8553177090 4743586537 78 D-12    /
      DATA BT12CS( 12) / -.4965314793 2768480785 5520211353 81 D-13    /
      DATA BT12CS( 13) / +.8043890043 2847825985 5588826393 17 D-14    /
      DATA BT12CS( 14) / -.1358912131 0161291384 6947126822 82 D-14    /
      DATA BT12CS( 15) / +.2381050439 7147214869 6765296059 73 D-15    /
      DATA BT12CS( 16) / -.4308146636 3849106724 4712414207 99 D-16    /
      DATA BT12CS( 17) / +.8020254403 2771002434 9935125504 00 D-17    /
      DATA BT12CS( 18) / -.1531631064 2462311864 2300274687 99 D-17    /
      DATA BT12CS( 19) / +.2992860635 2715568924 0730405546 66 D-18    /
      DATA BT12CS( 20) / -.5970996465 8085443393 8156366506 66 D-19    /
      DATA BT12CS( 21) / +.1214028966 9415185024 1608526506 66 D-19    /
      DATA BT12CS( 22) / -.2511511469 6612948901 0069777066 66 D-20    /
      DATA BT12CS( 23) / +.5279056717 0328744850 7383807999 99 D-21    /
      DATA BT12CS( 24) / -.1126050922 7550498324 3611613866 66 D-21    /
      DATA BT12CS( 25) / +.2434827735 9576326659 6634624000 00 D-22    /
      DATA BT12CS( 26) / -.5331726123 6931800130 0384426666 66 D-23    /
      DATA BT12CS( 27) / +.1181361505 9707121039 2059903999 99 D-23    /
      DATA BT12CS( 28) / -.2646536828 3353523514 8567893333 33 D-24    /
      DATA BT12CS( 29) / +.5990339404 1361503945 5778133333 33 D-25    /
      DATA BT12CS( 30) / -.1369085463 0829503109 1363839999 99 D-25    /
      DATA BT12CS( 31) / +.3157679015 4380228326 4136533333 33 D-26    /
      DATA BT12CS( 32) / -.7345791508 2084356491 4005333333 33 D-27    /
      DATA BT12CS( 33) / +.1722808148 0722747930 7059200000 00 D-27    /
      DATA BT12CS( 34) / -.4071690796 1286507941 0688000000 00 D-28    /
      DATA BT12CS( 35) / +.9693474513 6779622700 3733333333 33 D-29    /
      DATA BT12CS( 36) / -.2323763633 7765716765 3546666666 66 D-29    /
      DATA BT12CS( 37) / +.5607451067 3522029406 8906666666 66 D-30    /
      DATA BT12CS( 38) / -.1361646539 1539005860 5226666666 66 D-30    /
      DATA BT12CS( 39) / +.3326310923 3894654388 9066666666 66 D-31    /
C
C-------------------------------------------------------------------
C
      DATA NTJ1 / 0 /
C
C***FIRST EXECUTABLE STATEMENT  DWJ1
C
      IF (M .LE. 0) GO TO 910
C
      IF (NTJ1.EQ.0) THEN
         EPMACH = D1MACH(3)
         EPS = 0.10D0*EPMACH
         NTJ1  = IDWCS(BJ1CS, LBJ1, EPS)
         NTM1  = IDWCS(BM1CS, LBM1, EPS)
         NTM12 = IDWCS(BM12CS, LBM12, EPS)
         NTTH1 = IDWCS(BTH1CS, LBTH1, EPS)
         NTT12 = IDWCS(BT12CS, LBT12, EPS)
         XMIN = 2.0D0*D1MACH(1)
         XSML = SQRT(8.0D0*EPMACH)
         XMAX = 1.0D0/D1MACH(4)
         PI4 = ATAN(1.0D0)
         TPI4 = 3.0D0*PI4
      ENDIF
C
      NTOT = 0
C
      DO 10 I=1,M
         Y(I) = ABS(X(I))
  10  CONTINUE
C
      CALL DWNGT(M,Y,XMAX,KEY)
      IF (KEY .NE. 0) GO TO 920
C
C  ---------------------------
C  CASE   Y=0  OR  Y TOO SMALL
C  ---------------------------
C
C     NOTE --- J1 UNDERFLOWS FOR  0 .LT. Y .LE. XMIN
C
      DO 15 I=1,M
         F(I) = 0.0D0
  15  CONTINUE
C
C  ----------------------------
C  CASE   XMIN .LT. Y .LE. XSML
C  ----------------------------
C
      CALL DWGTLE(M,Y,XMIN,XSML,N,INDX)
      IF (N .GT. 0) THEN
         NTOT = NTOT + N
         CALL DWGTHR(N,Y,INDX,YCMP)
         DO 20 J=1,N
            ZCMP(J) = 0.50D0*YCMP(J)
  20     CONTINUE
         CALL DWSCTR(N,ZCMP,INDX,F)
      ENDIF
C
C  ---------------------------
C  CASE   XSML .LT. Y .LE. 4.0
C  ---------------------------
C
      CALL DWGTLE(M,Y,XSML,4.0D0,N,INDX)
      IF (N .GT. 0) THEN
         NTOT = NTOT + N
         CALL DWGTHR(N,Y,INDX,YCMP)
         DO 30 J=1,N
            TCMP(J) = 0.1250D0*YCMP(J)**2 - 1.0D0
  30     CONTINUE
         CALL DWCS(N,TCMP,BJ1CS,NTJ1,ZCMP,B0,B1,B2)
         DO 40 K=1,N
            ZCMP(K) = YCMP(K)*(0.250D0 + ZCMP(K))
  40     CONTINUE
         CALL DWSCTR(N,ZCMP,INDX,F)
      ENDIF
C
C  ----------------
C  CASE  Y .GT. 4.0
C  ----------------
C
      CALL DWGTLE(M,Y,4.0D0,8.0D0,NA,INDX)
      JH = NA + 1
      CALL DWGT(M,Y,8.0D0,NB,INDX(JH))
      N = NA + NB
      IF (N .GT. 0) THEN
         NTOT = NTOT + N
         CALL DWGTHR(N,Y,INDX,YCMP)
         C1 = 128.0D0/3.0D0
         C2 = 5.0D0/3.0D0
         DO 50 J=1,NA
            ZCMP(J) = C1/YCMP(J)**2 - C2
  50     CONTINUE
         CALL DWCS(NA,ZCMP,BM1CS,NTM1,Y,B0,B1,B2)
         CALL DWCS(NA,ZCMP,BT12CS,NTT12,ZCMP,B0,B1,B2)
         DO 60 J=JH,N
            ZCMP(J) = 128.0D0/YCMP(J)**2 - 1.0D0
  60     CONTINUE
         CALL DWCS(NB,ZCMP(JH),BM12CS,NTM12,Y(JH),B0,B1,B2)
         CALL DWCS(NB,ZCMP(JH),BTH1CS,NTTH1,ZCMP(JH),B0,B1,B2)
         DO 70 J=1,N
            Y(J) = (0.750D0 + Y(J)) / SQRT(YCMP(J))
            ZCMP(J) = (YCMP(J) - TPI4) + ZCMP(J) / YCMP(J)
            ZCMP(J) = Y(J) * COS(ZCMP(J))
  70     CONTINUE
         CALL DWSCTR(N,ZCMP,INDX,F)
      ENDIF
C
C  -------------------------------------------
C  REVERSE SIGN FOR NEGATIVE X (RESULT IS ODD)
C  -------------------------------------------
C
      DO 200 I=1,M
         IF (X(I) .LT. 0.0D0)  F(I) = -F(I)
  200 CONTINUE
C
C  -----
C  EXITS
C  -----
C
C     ... NORMAL
C
      INFO = 0
      IF (NTOT .NE. M)  INFO = -1
      GO TO 999
C
C     ... M .LE. 0
C
  910 CONTINUE
      INFO = 1
      GO TO 999
C
C     ... NO PRECISION BECAUSE X BIG
C
  920 CONTINUE
      INFO = 2
      INDX(1) = KEY
      GO TO 999
C
  999 CONTINUE
      RETURN
      END
      SUBROUTINE DVY0 (M, X, F, WORK, IWORK, INFO)
C***BEGIN PROLOGUE  DVY0
C***PURPOSE  Computes the Bessel function of the second kind
C            of order zero (Y0) for a vector of real arguments
C***LIBRARY   VFNLIB
C***CATEGORY  C10A1
C***TYPE      DOUBLE PRECISION (VY0-S, DVY0-D)
C***KEYWORDS  BESSEL FUNCTION,SECOND KIND, ORDER ZERO, VECTORIZED
C***AUTHOR  SAUNDERS, B. V., (NIST)
C             COMPUTING AND APPLIED MATHEMATICS LABORATORY
C             NATIONAL INSTITUTE OF STANDARDS AND TECHNOLOGY
C             GAITHERSBURG, MD 20899
C           BOISVERT, R. F., (NIST)
C             COMPUTING AND APPLIED MATHEMATICS LABORATORY
C             NATIONAL INSTITUTE OF STANDARDS AND TECHNOLOGY
C             GAITHERSBURG, MD 20899
C***DESCRIPTION
C
C   DVY0 computes the  Bessel function of the second kind of order
C   zero (Y0) for real arguments using uniform approximation by
C   Chebyshev polynomials.
C
C
C   P A R A M E T E R S
C
C   M      (Input) Integer (M .GT. 0)
C          The number of arguments at which the function is to be
C          evaluated.
C
C   X      (Input) Double precision array of length M
C          The arguments at which the function is to be evaluated are
C          stored in X(1) to X(M) in any order.
C
C   F      (Output) Double precision array of length M
C          F(i) contains the value of the function at X(i), i=1,..,M.
C
C   WORK   (Work) Double precision vector of length 7*M
C
C   IWORK  (Work) Integer vector of length M
C
C   INFO   (Output) Integer
C          Indicates status of computed result. The following table
C          lists possible values and their meanings.  If OK=Yes then
C          all F(i) have been set, otherwise none have been set.
C
C          INFO  OK            Description
C          ------------------------------------------------------------
C            0   Yes  Successfull execution.
C            1   No   Error: M .LE. 0
C            2   No   Error: X(i) zero or negative for some i
C            3   No   Error: Some X(i) so big that no precision possible
C                     in computing Y0. The index of the first offending
C                     argument is returned in IWORK(1).
C
C
C *********************************************************************
C   This routine is a modification of the function DBESY0 developed by
C   W. Fullerton of LANL.
C *********************************************************************
C
C***REFERENCES  (NONE)
C***ROUTINES CALLED  DWY0
C***REVISION HISTORY  (YYMMDD)
C   910226  DATE WRITTEN
C***END PROLOGUE  DVY0
C
C     ----------
C     PARAMETERS
C     ----------
C
      INTEGER INFO, IWORK, M
      DOUBLE PRECISION F, X, WORK
C
      DIMENSION X(M), F(M), WORK(7*M), IWORK(M)
C
C     ---------------
C     LOCAL VARIABLES
C     ---------------
C
      INTEGER IWB0, IWB1, IWB2, IWTC, IWXC, IWYC, IWZC, JIN
C
C
C***FIRST EXECUTABLE STATEMENT  DVY0
C
C     ... PARTITION WORK ARRAYS
C
      IWTC   = 1
      IWXC  = IWTC + M
      IWYC  = IWXC + M
      IWZC  = IWYC + M
      IWB0  = IWZC + M
      IWB1  = IWB0 + M
      IWB2  = IWB1 + M
C     Total = IB2  + M
C
      JIN   = 1
C     Total = JIN  + M
C
C     ... DWY0 DOES ALL THE WORK
C
      CALL DWY0(M,X,F,WORK(IWTC),WORK(IWXC),WORK(IWYC),WORK(IWZC),
     +         IWORK(JIN),WORK(IWB0),WORK(IWB1),WORK(IWB2),INFO)
C
      RETURN
      END

      SUBROUTINE DWY0 (M, X, F, TCMP, XCMP, YCMP, ZCMP, INDX, B0, B1,  
     +   B2, INFO)
C***BEGIN PROLOGUE  DWY0
C***SUBSIDIARY
C***PURPOSE  Computes the Bessel function of the second kind
C            of order zero (Y0) for a vector of arguments
C***LIBRARY   VFNLIB
C***AUTHOR  SAUNDERS, B. V., (NIST)
C             COMPUTING AND APPLIED MATHEMATICS LABORATORY
C             NATIONAL INSTITUTE OF STANDARDS AND TECHNOLOGY
C             GAITHERSBURG, MD 20899
C           BOISVERT, R. F., (NIST)
C             COMPUTING AND APPLIED MATHEMATICS LABORATORY
C             NATIONAL INSTITUTE OF STANDARDS AND TECHNOLOGY
C             GAITHERSBURG, MD 20899
C***ROUTINES CALLED  IDWCS, D1MACH, DWNGT, DWGTHR, DWGTLE, DWGT, DWLE,
C                    DWSCTR, DWCS
C***REVISION HISTORY  (YYMMDD)
C   910226  DATE WRITTEN
C***END PROLOGUE  DWY0
C
C  ----------
C  PARAMETERS
C  ----------
C
      INTEGER INFO, INDX, M
      DOUBLE PRECISION B0, B1, B2, F, X, TCMP, XCMP, YCMP, ZCMP
C
      DIMENSION B0(M), B1(M), B2(M), F(M), INDX(M), X(M), TCMP(M),
     +          XCMP(M), YCMP(M), ZCMP(M)
C
C  ---------------
C  LOCAL VARIABLES
C  ---------------
C
      INTEGER LBJ0, LBM0, LBM02, LBTH0, LBT02, LBY0
      PARAMETER (LBJ0=19, LBM0=37, LBM02=40, LBTH0=44, LBT02=39,
     +           LBY0=19)
C
      INTEGER IDWCS, J, JH, KEY, N, NA, NB, NTJ0, NTM0, NTM02, NTTH0,
     +        NTT02, NTY0
      DOUBLE PRECISION BJ0CS, BM0CS, BM02CS, BTH0CS, BT02CS, BY0CS, C1,
     +         C2, EPMACH, EPS, D1MACH, PI4, TWODPI, XSML, XMAX
C
      DIMENSION BY0CS(LBY0), BM0CS(LBM0), BM02CS(LBM02), BTH0CS(LBTH0),
     +          BT02CS(LBT02), BJ0CS(LBJ0)
C
      SAVE BJ0CS, BM0CS, BTH0CS, BY0CS, NTJ0, NTM0, NTTH0, NTY0,
     +     PI4, TWODPI, XSML, XMAX
C
C----------------------------------------------------------------------
C
C Series for BY0        on the interval  0.          to  1.60000D+01
C                                        with weighted error   8.14E-32
C                                         log weighted error  31.09
C                               significant figures required  30.31
C                                    decimal places required  31.73
C
      DATA BY0 CS(  1) / -.1127783939 2865573217 9398054602 8 D-1      /
      DATA BY0 CS(  2) / -.1283452375 6042034604 8088453183 8 D+0      /
      DATA BY0 CS(  3) / -.1043788479 9794249365 8176227661 8 D+0      /
      DATA BY0 CS(  4) / +.2366274918 3969695409 2415926461 3 D-1      /
      DATA BY0 CS(  5) / -.2090391647 7004862391 9622395034 2 D-2      /
      DATA BY0 CS(  6) / +.1039754539 3905725209 9924657638 1 D-3      /
      DATA BY0 CS(  7) / -.3369747162 4239720967 1877534503 7 D-5      /
      DATA BY0 CS(  8) / +.7729384267 6706671585 2136721637 1 D-7      /
      DATA BY0 CS(  9) / -.1324976772 6642595914 4347606896 4 D-8      /
      DATA BY0 CS( 10) / +.1764823261 5404527921 0038936315 8 D-10     /
      DATA BY0 CS( 11) / -.1881055071 5801962006 0282301206 9 D-12     /
      DATA BY0 CS( 12) / +.1641865485 3661495027 9223718574 9 D-14     /
      DATA BY0 CS( 13) / -.1195659438 6046060857 4599100672 0 D-16     /
      DATA BY0 CS( 14) / +.7377296297 4401858424 9411242666 6 D-19     /
      DATA BY0 CS( 15) / -.3906843476 7104373307 4090666666 6 D-21     /
      DATA BY0 CS( 16) / +.1795503664 4361579498 2912000000 0 D-23     /
      DATA BY0 CS( 17) / -.7229627125 4480104789 3333333333 3 D-26     /
      DATA BY0 CS( 18) / +.2571727931 6351685973 3333333333 3 D-28     /
      DATA BY0 CS( 19) / -.8141268814 1636949333 3333333333 3 D-31     /
C
C-------------------------------------------------------------------
C
C Series for BM0        on the interval  0.          to  6.25000D-02
C                                        with weighted error   4.40E-32
C                                         log weighted error  31.36
C                               significant figures required  30.02
C                                    decimal places required  32.14
C
      DATA BM0 CS(  1) / +.9211656246 8277427125 7376773018 2 D-1      /
      DATA BM0 CS(  2) / -.1050590997 2719051024 8071637175 5 D-2      /
      DATA BM0 CS(  3) / +.1470159840 7687597540 5639285095 2 D-4      /
      DATA BM0 CS(  4) / -.5058557606 0385542233 4792932770 2 D-6      /
      DATA BM0 CS(  5) / +.2787254538 6324441766 3035613788 1 D-7      /
      DATA BM0 CS(  6) / -.2062363611 7809148026 1884101897 3 D-8      /
      DATA BM0 CS(  7) / +.1870214313 1388796751 3817259626 1 D-9      /
      DATA BM0 CS(  8) / -.1969330971 1356362002 4173077782 5 D-10     /
      DATA BM0 CS(  9) / +.2325973793 9992754440 1250881805 2 D-11     /
      DATA BM0 CS( 10) / -.3009520344 9382502728 5122473448 2 D-12     /
      DATA BM0 CS( 11) / +.4194521333 8506691814 7120676864 6 D-13     /
      DATA BM0 CS( 12) / -.6219449312 1884458259 7326742956 4 D-14     /
      DATA BM0 CS( 13) / +.9718260411 3360684696 0176588526 9 D-15     /
      DATA BM0 CS( 14) / -.1588478585 7010752073 6663596693 7 D-15     /
      DATA BM0 CS( 15) / +.2700072193 6713088900 8621732445 8 D-16     /
      DATA BM0 CS( 16) / -.4750092365 2340089924 7750478677 3 D-17     /
      DATA BM0 CS( 17) / +.8615128162 6043708731 9170374656 0 D-18     /
      DATA BM0 CS( 18) / -.1605608686 9561448157 4560270335 9 D-18     /
      DATA BM0 CS( 19) / +.3066513987 3144829751 8853980159 9 D-19     /
      DATA BM0 CS( 20) / -.5987764223 1939564306 9650561706 6 D-20     /
      DATA BM0 CS( 21) / +.1192971253 7482483064 8906984106 6 D-20     /
      DATA BM0 CS( 22) / -.2420969142 0448054894 8468258133 3 D-21     /
      DATA BM0 CS( 23) / +.4996751760 5106164533 7100287999 9 D-22     /
      DATA BM0 CS( 24) / -.1047493639 3511585100 9504051199 9 D-22     /
      DATA BM0 CS( 25) / +.2227786843 7974681010 4818346666 6 D-23     /
      DATA BM0 CS( 26) / -.4801813239 3981628623 7054293333 3 D-24     /
      DATA BM0 CS( 27) / +.1047962723 4709599564 7699626666 6 D-24     /
      DATA BM0 CS( 28) / -.2313858165 6786153251 0126080000 0 D-25     /
      DATA BM0 CS( 29) / +.5164823088 4626742116 3519999999 9 D-26     /
      DATA BM0 CS( 30) / -.1164691191 8500653895 2540159999 9 D-26     /
      DATA BM0 CS( 31) / +.2651788486 0433192829 5833600000 0 D-27     /
      DATA BM0 CS( 32) / -.6092559503 8257284976 9130666666 6 D-28     /
      DATA BM0 CS( 33) / +.1411804686 1442593080 3882666666 6 D-28     /
      DATA BM0 CS( 34) / -.3298094961 2317372457 5061333333 3 D-29     /
      DATA BM0 CS( 35) / +.7763931143 0740650317 1413333333 3 D-30     /
      DATA BM0 CS( 36) / -.1841031343 6614584784 2133333333 3 D-30     /
      DATA BM0 CS( 37) / +.4395880138 5943107371 0079999999 9 D-31     /
C
C-------------------------------------------------------------------
C
C Series for BM02       Used in double precision version only
C                                        with weighted error   4.72E-32
C                                         log weighted error  31.33
C                               significant figures required  30.00
C                                    decimal places required  32.13
C
      DATA BM02CS(  1) / +.9500415145 2283813693 3086133556 0 D-1      /
      DATA BM02CS(  2) / -.3801864682 3656709917 4808156685 1 D-3      /
      DATA BM02CS(  3) / +.2258339301 0314811929 5182992722 4 D-5      /
      DATA BM02CS(  4) / -.3895725802 3722287647 3062141260 5 D-7      /
      DATA BM02CS(  5) / +.1246886416 5120816979 3099052972 5 D-8      /
      DATA BM02CS(  6) / -.6065949022 1025037798 0383505838 7 D-10     /
      DATA BM02CS(  7) / +.4008461651 4217469910 1527597104 5 D-11     /
      DATA BM02CS(  8) / -.3350998183 3980942184 6729879457 4 D-12     /
      DATA BM02CS(  9) / +.3377119716 5174173670 6326434199 6 D-13     /
      DATA BM02CS( 10) / -.3964585901 6350127005 6935629582 3 D-14     /
      DATA BM02CS( 11) / +.5286111503 8838572173 8793974473 5 D-15     /
      DATA BM02CS( 12) / -.7852519083 4508523136 5464024349 3 D-16     /
      DATA BM02CS( 13) / +.1280300573 3866822010 1163407344 9 D-16     /
      DATA BM02CS( 14) / -.2263996296 3914297762 8709924488 4 D-17     /
      DATA BM02CS( 15) / +.4300496929 6567903886 4641029047 7 D-18     /
      DATA BM02CS( 16) / -.8705749805 1325870797 4753545145 5 D-19     /
      DATA BM02CS( 17) / +.1865862713 9620951411 8144277205 0 D-19     /
      DATA BM02CS( 18) / -.4210482486 0930654573 4508697230 1 D-20     /
      DATA BM02CS( 19) / +.9956676964 2284009915 8162741784 2 D-21     /
      DATA BM02CS( 20) / -.2457357442 8053133596 0592147854 7 D-21     /
      DATA BM02CS( 21) / +.6307692160 7620315680 8735370705 9 D-22     /
      DATA BM02CS( 22) / -.1678773691 4407401426 9333117238 8 D-22     /
      DATA BM02CS( 23) / +.4620259064 6739044337 7087813608 7 D-23     /
      DATA BM02CS( 24) / -.1311782266 8603087322 3769340249 6 D-23     /
      DATA BM02CS( 25) / +.3834087564 1163028277 4792244027 6 D-24     /
      DATA BM02CS( 26) / -.1151459324 0777412710 7261329357 6 D-24     /
      DATA BM02CS( 27) / +.3547210007 5233385230 7697134521 3 D-25     /
      DATA BM02CS( 28) / -.1119218385 8150046462 6435594217 6 D-25     /
      DATA BM02CS( 29) / +.3611879427 6298378316 9840499425 7 D-26     /
      DATA BM02CS( 30) / -.1190687765 9133331500 9264176246 3 D-26     /
      DATA BM02CS( 31) / +.4005094059 4039681318 0247644953 6 D-27     /
      DATA BM02CS( 32) / -.1373169422 4522123905 9519391601 7 D-27     /
      DATA BM02CS( 33) / +.4794199088 7425315859 9649152643 7 D-28     /
      DATA BM02CS( 34) / -.1702965627 6241095840 0699447645 2 D-28     /
      DATA BM02CS( 35) / +.6149512428 9363300715 0357516132 4 D-29     /
      DATA BM02CS( 36) / -.2255766896 5818283499 4430023724 2 D-29     /
      DATA BM02CS( 37) / +.8399707509 2942994860 6165835320 0 D-30     /
      DATA BM02CS( 38) / -.3172997595 5626023555 6742393615 2 D-30     /
      DATA BM02CS( 39) / +.1215205298 8812985545 8333302651 4 D-30     /
      DATA BM02CS( 40) / -.4715852749 7544386930 1321056804 5 D-31     /
C
C-------------------------------------------------------------------
C
C Series for BTH0       on the interval  0.          to  6.25000D-02
C                                        with weighted error   2.66E-32
C                                         log weighted error  31.57
C                               significant figures required  30.67
C                                    decimal places required  32.40
C
      DATA BTH0CS(  1) / -.2490178086 2128936717 7097937899 67 D+0     /
      DATA BTH0CS(  2) / +.4855029960 9623749241 0486155354 85 D-3     /
      DATA BTH0CS(  3) / -.5451183734 5017204950 6562735635 05 D-5     /
      DATA BTH0CS(  4) / +.1355867305 9405964054 3774459299 03 D-6     /
      DATA BTH0CS(  5) / -.5569139890 2227626227 5832184149 20 D-8     /
      DATA BTH0CS(  6) / +.3260903182 4994335304 0042057194 68 D-9     /
      DATA BTH0CS(  7) / -.2491880786 2461341125 2379038779 93 D-10    /
      DATA BTH0CS(  8) / +.2344937742 0882520554 3524135648 91 D-11    /
      DATA BTH0CS(  9) / -.2609653444 4310387762 1775747661 36 D-12    /
      DATA BTH0CS( 10) / +.3335314042 0097395105 8699550149 23 D-13    /
      DATA BTH0CS( 11) / -.4789000044 0572684646 7507705574 09 D-14    /
      DATA BTH0CS( 12) / +.7595617843 6192215972 6425685452 48 D-15    /
      DATA BTH0CS( 13) / -.1313155601 6891440382 7733974876 33 D-15    /
      DATA BTH0CS( 14) / +.2448361834 5240857495 4268207383 55 D-16    /
      DATA BTH0CS( 15) / -.4880572981 0618777683 2567619183 31 D-17    /
      DATA BTH0CS( 16) / +.1032728502 9786316149 2237563612 04 D-17    /
      DATA BTH0CS( 17) / -.2305763381 5057217157 0047445270 25 D-18    /
      DATA BTH0CS( 18) / +.5404444300 1892693993 0171084837 65 D-19    /
      DATA BTH0CS( 19) / -.1324069519 4366572724 1550328823 85 D-19    /
      DATA BTH0CS( 20) / +.3378079562 1371970203 4247921247 22 D-20    /
      DATA BTH0CS( 21) / -.8945762915 7111779003 0269262922 99 D-21    /
      DATA BTH0CS( 22) / +.2451990688 9219317090 8999086514 05 D-21    /
      DATA BTH0CS( 23) / -.6938842287 6866318680 1399331576 57 D-22    /
      DATA BTH0CS( 24) / +.2022827871 4890138392 9463033377 91 D-22    /
      DATA BTH0CS( 25) / -.6062850000 2335483105 7941953717 64 D-23    /
      DATA BTH0CS( 26) / +.1864974896 4037635381 8237883962 70 D-23    /
      DATA BTH0CS( 27) / -.5878373238 4849894560 2450365308 67 D-24    /
      DATA BTH0CS( 28) / +.1895859144 7999563485 5311795035 13 D-24    /
      DATA BTH0CS( 29) / -.6248197937 2258858959 2916207285 65 D-25    /
      DATA BTH0CS( 30) / +.2101790168 4551024686 6386335290 74 D-25    /
      DATA BTH0CS( 31) / -.7208430093 5209253690 8139339924 46 D-26    /
      DATA BTH0CS( 32) / +.2518136389 2474240867 1564059767 46 D-26    /
      DATA BTH0CS( 33) / -.8951804225 8785778806 1439459536 43 D-27    /
      DATA BTH0CS( 34) / +.3235723747 9762298533 2562358685 87 D-27    /
      DATA BTH0CS( 35) / -.1188301051 9855353657 0471441137 96 D-27    /
      DATA BTH0CS( 36) / +.4430628690 7358104820 5792319417 31 D-28    /
      DATA BTH0CS( 37) / -.1676100964 8834829495 7920101356 81 D-28    /
      DATA BTH0CS( 38) / +.6429294692 1207466972 5323939660 88 D-29    /
      DATA BTH0CS( 39) / -.2499226116 6978652421 2072136827 63 D-29    /
      DATA BTH0CS( 40) / +.9839979429 9521955672 8282603553 18 D-30    /
      DATA BTH0CS( 41) / -.3922037524 2408016397 9891316261 58 D-30    /
      DATA BTH0CS( 42) / +.1581810703 0056522138 5906188456 92 D-30    /
      DATA BTH0CS( 43) / -.6452550614 4890715944 3440983654 26 D-31    /
      DATA BTH0CS( 44) / +.2661111136 9199356137 1770183463 67 D-31    /
C
C-------------------------------------------------------------------
C
C Series for BT02       Used in double precision version only
C                                        with weighted error   2.99E-32
C                                         log weighted error  31.52
C
      DATA BT02CS(  1) / -.2454829521 3424597462 0504672493 24 D+0     /
      DATA BT02CS(  2) / +.1254412103 9084615780 7853317782 99 D-2     /
      DATA BT02CS(  3) / -.3125395041 4871522854 9734467095 71 D-4     /
      DATA BT02CS(  4) / +.1470977824 9940831164 4534269693 14 D-5     /
      DATA BT02CS(  5) / -.9954348893 7950033643 4688503511 58 D-7     /
      DATA BT02CS(  6) / +.8549316673 3203041247 5787113977 51 D-8     /
      DATA BT02CS(  7) / -.8698975952 6554334557 9855121791 92 D-9     /
      DATA BT02CS(  8) / +.1005209953 3559791084 5401010821 53 D-9     /
      DATA BT02CS(  9) / -.1282823060 1708892903 4836236855 44 D-10    /
      DATA BT02CS( 10) / +.1773170078 1805131705 6557504510 23 D-11    /
      DATA BT02CS( 11) / -.2617457456 9485577488 6362841809 25 D-12    /
      DATA BT02CS( 12) / +.4082835138 9972059621 9664812211 03 D-13    /
      DATA BT02CS( 13) / -.6675166823 9742720054 6067495542 61 D-14    /
      DATA BT02CS( 14) / +.1136576139 3071629448 3924695499 51 D-14    /
      DATA BT02CS( 15) / -.2005118962 0647160250 5592664121 17 D-15    /
      DATA BT02CS( 16) / +.3649797879 4766269635 7205914641 06 D-16    /
      DATA BT02CS( 17) / -.6830963756 4582303169 3558437888 00 D-17    /
      DATA BT02CS( 18) / +.1310758314 5670756620 0571042679 46 D-17    /
      DATA BT02CS( 19) / -.2572336310 1850607778 7571306495 99 D-18    /
      DATA BT02CS( 20) / +.5152165744 1863959925 2677809493 33 D-19    /
      DATA BT02CS( 21) / -.1051301756 3758802637 9407414613 33 D-19    /
      DATA BT02CS( 22) / +.2182038199 1194813847 3010845013 33 D-20    /
      DATA BT02CS( 23) / -.4600470121 0362160577 2259054933 33 D-21    /
      DATA BT02CS( 24) / +.9840700692 5466818520 9536511999 99 D-22    /
      DATA BT02CS( 25) / -.2133403803 5728375844 7359863466 66 D-22    /
      DATA BT02CS( 26) / +.4683103642 3973365296 0662869333 33 D-23    /
      DATA BT02CS( 27) / -.1040021369 1985747236 5133823999 99 D-23    /
      DATA BT02CS( 28) / +.2334910567 7301510051 7777408000 00 D-24    /
      DATA BT02CS( 29) / -.5295682532 3318615788 0497493333 33 D-25    /
      DATA BT02CS( 30) / +.1212634195 2959756829 1962879999 99 D-25    /
      DATA BT02CS( 31) / -.2801889708 2289428760 2756266666 66 D-26    /
      DATA BT02CS( 32) / +.6529267898 7012873342 5937066666 66 D-27    /
      DATA BT02CS( 33) / -.1533798006 1873346427 8357333333 33 D-27    /
      DATA BT02CS( 34) / +.3630588430 6364536682 3594666666 66 D-28    /
      DATA BT02CS( 35) / -.8656075571 3629122479 1722666666 66 D-29    /
      DATA BT02CS( 36) / +.2077990997 2536284571 2383999999 99 D-29    /
      DATA BT02CS( 37) / -.5021117022 1417221674 3253333333 33 D-30    /
      DATA BT02CS( 38) / +.1220836027 9441714184 1919999999 99 D-30    /
      DATA BT02CS( 39) / -.2986005626 7039913454 2506666666 66 D-31    /
C
C-------------------------------------------------------------------
C
C Series for BJ0        on the interval  0.          to  1.60000D+01
C                                        with weighted error   4.39E-32
C                                         log weighted error  31.36
C                               significant figures required  31.21
C                                    decimal places required  32.00
C
      DATA BJ0 CS(  1) / +.1002541619 6893913701 0731272640 74 D+0     /
      DATA BJ0 CS(  2) / -.6652230077 6440513177 6787578311 24 D+0     /
      DATA BJ0 CS(  3) / +.2489837034 9828131370 4604687266 80 D+0     /
      DATA BJ0 CS(  4) / -.3325272317 0035769653 8843415038 54 D-1     /
      DATA BJ0 CS(  5) / +.2311417930 4694015462 9049241177 29 D-2     /
      DATA BJ0 CS(  6) / -.9911277419 9508092339 0485193365 49 D-4     /
      DATA BJ0 CS(  7) / +.2891670864 3998808884 7339037470 78 D-5     /
      DATA BJ0 CS(  8) / -.6121085866 3032635057 8184074815 16 D-7     /
      DATA BJ0 CS(  9) / +.9838650793 8567841324 7687486364 15 D-9     /
      DATA BJ0 CS( 10) / -.1242355159 7301765145 5158970068 36 D-10    /
      DATA BJ0 CS( 11) / +.1265433630 2559045797 9158272103 63 D-12    /
      DATA BJ0 CS( 12) / -.1061945649 5287244546 9148175129 59 D-14    /
      DATA BJ0 CS( 13) / +.7470621075 8024567437 0989155840 00 D-17    /
      DATA BJ0 CS( 14) / -.4469703227 4412780547 6270079999 99 D-19    /
      DATA BJ0 CS( 15) / +.2302428158 4337436200 5230933333 33 D-21    /
      DATA BJ0 CS( 16) / -.1031914479 4166698148 5226666666 66 D-23    /
      DATA BJ0 CS( 17) / +.4060817827 4873322700 8000000000 00 D-26    /
      DATA BJ0 CS( 18) / -.1414383600 5240913919 9999999999 99 D-28    /
      DATA BJ0 CS( 19) / +.4391090549 6698880000 0000000000 00 D-31    /
C
C-------------------------------------------------------------------
C
      DATA NTY0   / 0 /
C
C***FIRST EXECUTABLE STATEMENT  DWY0
C
      IF (M .LE. 0)  GOTO 910
C
      IF (NTY0 .EQ. 0) THEN
         PI4 = ATAN(1.0D0)
         TWODPI = 0.50D0/PI4
         EPMACH = D1MACH(3)
         EPS = 0.10D0*EPMACH
         NTY0 = IDWCS(BY0CS, LBY0, EPS)
         NTM0 = IDWCS(BM0CS, LBM0, EPS)
         NTM02 = IDWCS(BM02CS, LBM02, EPS)
         NTTH0 = IDWCS(BTH0CS, LBTH0, EPS)
         NTT02 = IDWCS(BT02CS, LBT02, EPS)
         NTJ0 = IDWCS(BJ0CS, LBJ0, EPS)
         XSML = SQRT (4.0D0*EPMACH)
         XMAX = 1.0D0/D1MACH(4)
      ENDIF
C
      CALL DWNLE(M,X,0.0D0,KEY)
      IF (KEY .NE. 0) GO TO 920
C
      CALL DWNGT(M,X,XMAX,KEY)
      IF (KEY .NE. 0) GO TO 930
C
C  ----------------
C  CASE  X .LE. 4.0
C  ----------------
C
      CALL DWLE(M,X,4.0D0,N,INDX)
      IF (N .GT. 0) THEN
         CALL DWGTHR(N,X,INDX,XCMP)
C
C        ... COMPUTE J0(X) ... RESULT IN ZCMP
C
         DO 10 J=1,N
            TCMP(J) = 0.125E0*XCMP(J)**2 - 1.0D0
   10    CONTINUE
         CALL DWCS(N,TCMP,BJ0CS,NTJ0,ZCMP,B0,B1,B2)
C
         CALL DWCS(N,TCMP,BY0CS,NTY0,YCMP,B0,B1,B2)
         DO 30 J=1,N
            YCMP(J) = TWODPI*LOG(0.50D0*XCMP(J))*ZCMP(J)
     +                 + 0.375E0 + YCMP(J)
   30    CONTINUE
         CALL DWSCTR(N,YCMP,INDX,F)
      ENDIF
C
C  ----------------
C  CASE  X .GT. 4.0
C  ----------------
C
      CALL DWGTLE(M,X,4.0D0,8.0D0,NA,INDX)
      JH = NA + 1
      CALL DWGT(M,X,8.0D0,NB,INDX(JH))
      N = NA + NB
      IF (N .GT. 0) THEN
         CALL DWGTHR(N,X,INDX,XCMP)
         C1 = 128.0D0/3.0D0
         C2 = 5.0D0/3.0D0
         DO 50 J=1,NA
            ZCMP(J) = C1/XCMP(J)**2 - C2
  50     CONTINUE
         CALL DWCS(NA,ZCMP,BM0CS,NTM0,YCMP,B0,B1,B2)
         CALL DWCS(NA,ZCMP,BT02CS,NTT02,ZCMP,B0,B1,B2)
         DO 60 J=JH,N
            ZCMP(J) = 128.0D0/XCMP(J)**2 - 1.0D0
  60     CONTINUE
         CALL DWCS(NB,ZCMP(JH),BM02CS,NTM02,YCMP(JH),B0,B1,B2)
         CALL DWCS(NB,ZCMP(JH),BTH0CS,NTTH0,ZCMP(JH),B0,B1,B2)
         DO 70 J=1,N
            YCMP(J) = (0.750D0 + YCMP(J)) / SQRT(XCMP(J))
            ZCMP(J) = (XCMP(J) - PI4) + ZCMP(J) / XCMP(J)
            ZCMP(J) = YCMP(J) * SIN(ZCMP(J))
  70     CONTINUE
         CALL DWSCTR(N,ZCMP,INDX,F)
      ENDIF

C
C  -----
C  EXITS
C  -----
C
C     ... NORMAL
C
      INFO = 0
      GO TO 999
C
C     ... M .LE. 0
C
  910 CONTINUE
      INFO = 1
      GO TO 999
C
C     ... X I ZERO OR NEGATIVE
C
  920 CONTINUE
      INFO = 2
      INDX(1) = KEY
      GO TO 999
C
C     ... NO PRECISION BECAUSE X BIG
C
  930 CONTINUE
      INFO = 3
      INDX(1) = KEY
      GO TO 999
C
  999 CONTINUE
      RETURN
      END
      SUBROUTINE DVY1 (M, X, F, WORK, IWORK, INFO)
C***BEGIN PROLOGUE  DVY1
C***PURPOSE  Computes the Bessel function of the second kind
C            of order one (Y1) for a vector of real arguments
C***LIBRARY   VFNLIB
C***CATEGORY  C10A1
C***TYPE      DOUBLE PRECISION (VY1-S, DVY1-D)
C***KEYWORDS  BESSEL FUNCTION,SECOND KIND, ORDER ONE, VECTORIZED
C***AUTHOR  SAUNDERS, B. V., (NIST)
C             COMPUTING AND APPLIED MATHEMATICS LABORATORY
C             NATIONAL INSTITUTE OF STANDARDS AND TECHNOLOGY
C             GAITHERSBURG, MD 20899
C           BOISVERT, R. F., (NIST)
C             COMPUTING AND APPLIED MATHEMATICS LABORATORY
C             NATIONAL INSTITUTE OF STANDARDS AND TECHNOLOGY
C             GAITHERSBURG, MD 20899
C***DESCRIPTION
C
C   DVY1 computes the  Bessel function of the second kind of order
C   one (Y1) for real arguments using uniform approximation by
C   Chebyshev polynomials.
C
C
C   P A R A M E T E R S
C
C   M      (Input) Integer (M .GT. 0)
C          The number of arguments at which the function is to be
C          evaluated.
C
C   X      (Input) Double precision array of length M
C          The arguments at which the function is to be evaluated are
C          stored in X(1) to X(M) in any order.
C
C   F      (Output) Double precision array of length M
C          F(i) contains the value of the function at X(i), i=1,..,M.
C
C   WORK   (Work) Double precision vector of length 7*M
C
C   IWORK  (Work) Integer vector of length M
C
C   INFO   (Output) Integer
C          Indicates status of computed result. The following table
C          lists possible values and their meanings.  If OK=Yes then
C          all F(i) have been set, otherwise none have been set.
C
C          INFO  OK            Description
C          ------------------------------------------------------------
C            0   Yes  Successfull execution.
C            1   No   Error: M .LE. 0
C            2   No   Error: X(i) zero or negative for some i
C            3   No   Error: Some X(i) so small Y1 overflows.
C                     The index of the first offending argument is
C                     returned in IWORK(1).
C
C
C *********************************************************************
C   This routine is a modification of the function DBESY1 developed by
C   W. Fullerton of LANL.
C *********************************************************************
C
C***REFERENCES  (NONE)
C***ROUTINES CALLED  DWY1
C***REVISION HISTORY  (YYMMDD)
C   910226  DATE WRITTEN
C***END PROLOGUE  DVY1
C
C  ----------
C  PARAMETERS
C  ----------
      INTEGER INFO, IWORK, M
      DOUBLE PRECISION F, X, WORK
C
      DIMENSION X(M), F(M), WORK(7*M), IWORK(M)
C
C ---------------
C  LOCAL VARIABLES
C  ---------------
C
      INTEGER IWB0, IWB1, IWB2, IWTC, IWXC, IWYC, IWZC, JIN
C
C
C***FIRST EXECUTABLE STATEMENT  DVY1
C
C     ... PARTITION WORK ARRAYS
C
      IWTC   = 1
      IWXC  = IWTC + M
      IWYC  = IWXC + M
      IWZC  = IWYC + M
      IWB0  = IWZC + M
      IWB1  = IWB0 + M
      IWB2  = IWB1 + M
C     Total = IB2  + M
C
      JIN   = 1
C     Total = JIN  + M
C
C     ... DWY1 DOES ALL THE WORK
C
      CALL DWY1(M,X,F,WORK(IWTC),WORK(IWXC),WORK(IWYC),WORK(IWZC),
     +         IWORK(JIN),WORK(IWB0),WORK(IWB1),WORK(IWB2),INFO)
C
      RETURN
      END

      SUBROUTINE DWY1 (M, X, F, TCMP, XCMP, YCMP, ZCMP, INDX, B0, B1,  
     +   B2, INFO)
C***BEGIN PROLOGUE  DWY1
C***SUBSIDIARY
C***PURPOSE  Computes the Bessel function of the second kind
C            of order one (Y1) for a vector of arguments
C***LIBRARY   VFNLIB
C***AUTHOR  SAUNDERS, B. V., (NIST)
C             COMPUTING AND APPLIED MATHEMATICS LABORATORY
C             NATIONAL INSTITUTE OF STANDARDS AND TECHNOLOGY
C             GAITHERSBURG, MD 20899
C           BOISVERT, R. F., (NIST)
C             COMPUTING AND APPLIED MATHEMATICS LABORATORY
C             NATIONAL INSTITUTE OF STANDARDS AND TECHNOLOGY
C             GAITHERSBURG, MD 20899
C***ROUTINES CALLED  IDWCS, D1MACH, DWNGT, DWGTHR, DWGTLE, DWGT, DWLE,
C                    DWSCTR, DWCS, DWNLE
C***REVISION HISTORY  (YYMMDD)
C   910226  DATE WRITTEN
C***END PROLOGUE  DWY1
C
C  ----------
C  PARAMETERS
C  ----------
C
      INTEGER INFO, INDX, M
      DOUBLE PRECISION B0, B1, B2, F, X, TCMP, XCMP, YCMP, ZCMP
C
      DIMENSION B0(M), B1(M), B2(M), F(M), INDX(M), X(M), TCMP(M),
     +          XCMP(M), YCMP(M), ZCMP(M)
C
C  ---------------
C  LOCAL VARIABLES
C  ---------------
C
      INTEGER LBJ1, LBM1, LBM12, LBTH1, LBT12, LBY1
      PARAMETER (LBJ1=19, LBM1=37, LBM12=40, LBTH1=44, LBT12=39,
     +           LBY1=20)
C
      INTEGER IDWCS, J, JH, KEY, N, NA, NB, NTJ1, NTM1, NTM12, NTTH1,
     +        NTT12, NTY1
      DOUBLE PRECISION BJ1CS, BM1CS, BM12CS, BTH1CS, BT12CS, BY1CS,
     +        C1, C2, EPMACH, EPS, PI4, D1MACH, TPI4, TWODPI, XMIN,
     +        XSML, XMAX
C
      DIMENSION BJ1CS(LBJ1), BM1CS(LBM1), BM12CS(LBM12), BTH1CS(LBTH1),
     +          BT12CS(LBT12), BY1CS(LBY1)
C
      SAVE BJ1CS, BM1CS, BTH1CS, BY1CS, NTJ1, NTM1, NTM12, NTTH1, NTT12,
     +      NTY1, PI4, TPI4, TWODPI, XMIN, XSML, XMAX
C
C----------------------------------------------------------------------
C
C
C Series for BY1        on the interval  0.          to  1.60000D+01
C                                        with weighted error   8.65E-33
C                                         log weighted error  32.06
C                               significant figures required  32.17
C                                    decimal places required  32.71
C
      DATA BY1 CS(  1) / +.3208047100 6119086293 2352018628 015 D-1    /
      DATA BY1 CS(  2) / +.1262707897 4335004495 3431725999 727 D+1    /
      DATA BY1 CS(  3) / +.6499961899 9231750009 7490637314 144 D-2    /
      DATA BY1 CS(  4) / -.8936164528 8605041165 3144160009 712 D-1    /
      DATA BY1 CS(  5) / +.1325088122 1757095451 2375510370 043 D-1    /
      DATA BY1 CS(  6) / -.8979059119 6483523775 3039508298 105 D-3    /
      DATA BY1 CS(  7) / +.3647361487 9583067824 2287368165 349 D-4    /
      DATA BY1 CS(  8) / -.1001374381 6660005554 9075523845 295 D-5    /
      DATA BY1 CS(  9) / +.1994539657 3901739703 1159372421 243 D-7    /
      DATA BY1 CS( 10) / -.3023065601 8033816728 4799332520 743 D-9    /
      DATA BY1 CS( 11) / +.3609878156 9478119611 6252914242 474 D-11   /
      DATA BY1 CS( 12) / -.3487488297 2875824241 4552947409 066 D-13   /
      DATA BY1 CS( 13) / +.2783878971 5591766581 3507698517 333 D-15   /
      DATA BY1 CS( 14) / -.1867870968 6194876876 6825352533 333 D-17   /
      DATA BY1 CS( 15) / +.1068531533 9116825975 7070336000 000 D-19   /
      DATA BY1 CS( 16) / -.5274721956 6844822894 3872000000 000 D-22   /
      DATA BY1 CS( 17) / +.2270199403 1556641437 0133333333 333 D-24   /
      DATA BY1 CS( 18) / -.8595390353 9452310869 3333333333 333 D-27   /
      DATA BY1 CS( 19) / +.2885404379 8337945600 0000000000 000 D-29   /
      DATA BY1 CS( 20) / -.8647541138 9371733333 3333333333 333 D-32   /
C
C-------------------------------------------------------------------
C
C Series for BM1        on the interval  0.          to  6.25000D-02
C                                        with weighted error   4.91E-32
C                                         log weighted error  31.31
C                               significant figures required  30.04
C                                    decimal places required  32.09
C
      DATA BM1 CS(  1) / +.1069845452 6180630149 6998530853 8 D+0      /
      DATA BM1 CS(  2) / +.3274915039 7159649007 2905514344 5 D-2      /
      DATA BM1 CS(  3) / -.2987783266 8316985920 3044577793 8 D-4      /
      DATA BM1 CS(  4) / +.8331237177 9919745313 9322266902 3 D-6      /
      DATA BM1 CS(  5) / -.4112665690 3020073048 9638172549 8 D-7      /
      DATA BM1 CS(  6) / +.2855344228 7892152207 1975766316 1 D-8      /
      DATA BM1 CS(  7) / -.2485408305 4156238780 6002659605 5 D-9      /
      DATA BM1 CS(  8) / +.2543393338 0725824427 4248439717 4 D-10     /
      DATA BM1 CS(  9) / -.2941045772 8229675234 8975082790 9 D-11     /
      DATA BM1 CS( 10) / +.3743392025 4939033092 6505615362 6 D-12     /
      DATA BM1 CS( 11) / -.5149118293 8211672187 2054824352 7 D-13     /
      DATA BM1 CS( 12) / +.7552535949 8651439080 3404076419 9 D-14     /
      DATA BM1 CS( 13) / -.1169409706 8288464441 6629062246 4 D-14     /
      DATA BM1 CS( 14) / +.1896562449 4347915717 2182460506 0 D-15     /
      DATA BM1 CS( 15) / -.3201955368 6932864206 6477531639 4 D-16     /
      DATA BM1 CS( 16) / +.5599548399 3162041144 8416990549 3 D-17     /
      DATA BM1 CS( 17) / -.1010215894 7304324431 1939044454 4 D-17     /
      DATA BM1 CS( 18) / +.1873844985 7275629833 0204271957 3 D-18     /
      DATA BM1 CS( 19) / -.3563537470 3285802192 7430143999 9 D-19     /
      DATA BM1 CS( 20) / +.6931283819 9712383304 2276351999 9 D-20     /
      DATA BM1 CS( 21) / -.1376059453 4065001522 5140893013 3 D-20     /
      DATA BM1 CS( 22) / +.2783430784 1070802205 9977932799 9 D-21     /
      DATA BM1 CS( 23) / -.5727595364 3205616893 4866943999 9 D-22     /
      DATA BM1 CS( 24) / +.1197361445 9188926725 3575679999 9 D-22     /
      DATA BM1 CS( 25) / -.2539928509 8918719766 4144042666 6 D-23     /
      DATA BM1 CS( 26) / +.5461378289 6572959730 6961919999 9 D-24     /
      DATA BM1 CS( 27) / -.1189211341 7733202889 8628949333 3 D-24     /
      DATA BM1 CS( 28) / +.2620150977 3400815949 5782400000 0 D-25     /
      DATA BM1 CS( 29) / -.5836810774 2556859019 2093866666 6 D-26     /
      DATA BM1 CS( 30) / +.1313743500 0805957734 2361599999 9 D-26     /
      DATA BM1 CS( 31) / -.2985814622 5103803553 3277866666 6 D-27     /
      DATA BM1 CS( 32) / +.6848390471 3346049376 2559999999 9 D-28     /
      DATA BM1 CS( 33) / -.1584401568 2224767211 9296000000 0 D-28     /
      DATA BM1 CS( 34) / +.3695641006 5709380543 0101333333 3 D-29     /
      DATA BM1 CS( 35) / -.8687115921 1446682430 1226666666 6 D-30     /
      DATA BM1 CS( 36) / +.2057080846 1587634629 2906666666 6 D-30     /
      DATA BM1 CS( 37) / -.4905225761 1162255185 2373333333 3 D-31     /
C
C-------------------------------------------------------------------
C
C Series for BM12       Used in double precision version only
C                                        with weighted error   5.01E-32
C                                         log weighted error  31.30
C                               significant figures required  29.99
C                                    decimal places required  32.10
C
      DATA BM12CS(  1) / +.9807979156 2330500272 7209354693 7 D-1      /
      DATA BM12CS(  2) / +.1150961189 5046853061 7548348460 2 D-2      /
      DATA BM12CS(  3) / -.4312482164 3382054098 8935809773 2 D-5      /
      DATA BM12CS(  4) / +.5951839610 0888163078 1302980183 2 D-7      /
      DATA BM12CS(  5) / -.1704844019 8269098574 0070158647 8 D-8      /
      DATA BM12CS(  6) / +.7798265413 6111095086 5817382740 1 D-10     /
      DATA BM12CS(  7) / -.4958986126 7664158094 9175495186 5 D-11     /
      DATA BM12CS(  8) / +.4038432416 4211415168 3820226514 4 D-12     /
      DATA BM12CS(  9) / -.3993046163 7251754457 6548384664 5 D-13     /
      DATA BM12CS( 10) / +.4619886183 1189664943 1334243277 5 D-14     /
      DATA BM12CS( 11) / -.6089208019 0953833013 4547261933 3 D-15     /
      DATA BM12CS( 12) / +.8960930916 4338764821 5704804124 9 D-16     /
      DATA BM12CS( 13) / -.1449629423 9420231229 1651891892 5 D-16     /
      DATA BM12CS( 14) / +.2546463158 5377760561 6514964806 8 D-17     /
      DATA BM12CS( 15) / -.4809472874 6478364442 5926371862 0 D-18     /
      DATA BM12CS( 16) / +.9687684668 2925990490 8727583912 4 D-19     /
      DATA BM12CS( 17) / -.2067213372 2779660232 4503811755 1 D-19     /
      DATA BM12CS( 18) / +.4646651559 1503847318 0276780959 0 D-20     /
      DATA BM12CS( 19) / -.1094966128 8483341382 4135132833 9 D-20     /
      DATA BM12CS( 20) / +.2693892797 2886828609 0570761278 5 D-21     /
      DATA BM12CS( 21) / -.6894992910 9303744778 1897002685 7 D-22     /
      DATA BM12CS( 22) / +.1830268262 7520629098 9066855474 0 D-22     /
      DATA BM12CS( 23) / -.5025064246 3519164281 5611355322 4 D-23     /
      DATA BM12CS( 24) / +.1423545194 4548060396 3169363419 4 D-23     /
      DATA BM12CS( 25) / -.4152191203 6164503880 6888676980 1 D-24     /
      DATA BM12CS( 26) / +.1244609201 5039793258 8233007654 7 D-24     /
      DATA BM12CS( 27) / -.3827336370 5693042994 3191866128 6 D-25     /
      DATA BM12CS( 28) / +.1205591357 8156175353 7472398183 5 D-25     /
      DATA BM12CS( 29) / -.3884536246 3764880764 3185936112 4 D-26     /
      DATA BM12CS( 30) / +.1278689528 7204097219 0489528346 1 D-26     /
      DATA BM12CS( 31) / -.4295146689 4479462720 6193691591 2 D-27     /
      DATA BM12CS( 32) / +.1470689117 8290708864 5680270798 3 D-27     /
      DATA BM12CS( 33) / -.5128315665 1060731281 8037401779 6 D-28     /
      DATA BM12CS( 34) / +.1819509585 4711693854 8143737328 6 D-28     /
      DATA BM12CS( 35) / -.6563031314 8419808676 1863505037 3 D-29     /
      DATA BM12CS( 36) / +.2404898976 9199606531 9891487583 4 D-29     /
      DATA BM12CS( 37) / -.8945966744 6906124732 3495824297 9 D-30     /
      DATA BM12CS( 38) / +.3376085160 6572310266 3714897824 0 D-30     /
      DATA BM12CS( 39) / -.1291791454 6206563609 1309991696 6 D-30     /
      DATA BM12CS( 40) / +.5008634462 9588105206 8495150125 4 D-31     /
C
C-------------------------------------------------------------------
C
C Series for BTH1       on the interval  0.          to  6.25000D-02
C                                        with weighted error   2.82E-32
C                                         log weighted error  31.55
C                               significant figures required  31.12
C
      DATA BTH1CS(  1) / +.7474995720 3587276055 4434839696 95 D+0     /
      DATA BTH1CS(  2) / -.1240077714 4651711252 5457775413 84 D-2     /
      DATA BTH1CS(  3) / +.9925244240 4424527376 6414976895 92 D-5     /
      DATA BTH1CS(  4) / -.2030369073 7159711052 4193753756 08 D-6     /
      DATA BTH1CS(  5) / +.7535961770 5690885712 1840175836 29 D-8     /
      DATA BTH1CS(  6) / -.4166161271 5343550107 6300238562 28 D-9     /
      DATA BTH1CS(  7) / +.3070161807 0834890481 2451020912 16 D-10    /
      DATA BTH1CS(  8) / -.2817849963 7605213992 3240088839 24 D-11    /
      DATA BTH1CS(  9) / +.3079069673 9040295476 0281468216 47 D-12    /
      DATA BTH1CS( 10) / -.3880330026 2803434112 7873475547 81 D-13    /
      DATA BTH1CS( 11) / +.5509603960 8630904934 5617262085 62 D-14    /
      DATA BTH1CS( 12) / -.8659006076 8383779940 1033989539 94 D-15    /
      DATA BTH1CS( 13) / +.1485604914 1536749003 4236890606 83 D-15    /
      DATA BTH1CS( 14) / -.2751952981 5904085805 3712121250 09 D-16    /
      DATA BTH1CS( 15) / +.5455079609 0481089625 0362236409 23 D-17    /
      DATA BTH1CS( 16) / -.1148653450 1983642749 5436310271 77 D-17    /
      DATA BTH1CS( 17) / +.2553521337 7973900223 1990525335 22 D-18    /
      DATA BTH1CS( 18) / -.5962149019 7413450395 7682879078 49 D-19    /
      DATA BTH1CS( 19) / +.1455662290 2372718620 2883020058 33 D-19    /
      DATA BTH1CS( 20) / -.3702218542 2450538201 5797760195 93 D-20    /
      DATA BTH1CS( 21) / +.9776307412 5345357664 1684345179 24 D-21    /
      DATA BTH1CS( 22) / -.2672682163 9668488468 7237753930 52 D-21    /
      DATA BTH1CS( 23) / +.7545330038 4983271794 0381906557 64 D-22    /
      DATA BTH1CS( 24) / -.2194789991 9802744897 8923833716 47 D-22    /
      DATA BTH1CS( 25) / +.6564839462 3955262178 9069998174 93 D-23    /
      DATA BTH1CS( 26) / -.2015560429 8370207570 7840768695 19 D-23    /
      DATA BTH1CS( 27) / +.6341776855 6776143492 1446671856 70 D-24    /
      DATA BTH1CS( 28) / -.2041927788 5337895634 8137699555 91 D-24    /
      DATA BTH1CS( 29) / +.6719146422 0720567486 6589800185 51 D-25    /
      DATA BTH1CS( 30) / -.2256907911 0207573595 7090036873 36 D-25    /
      DATA BTH1CS( 31) / +.7729771989 2989706370 9269598719 29 D-26    /
      DATA BTH1CS( 32) / -.2696744451 2294640913 2114240809 20 D-26    /
      DATA BTH1CS( 33) / +.9574934451 8502698072 2955219336 27 D-27    /
      DATA BTH1CS( 34) / -.3456916844 8890113000 1756808276 27 D-27    /
      DATA BTH1CS( 35) / +.1268123481 7398436504 2119862383 74 D-27    /
      DATA BTH1CS( 36) / -.4723253663 0722639860 4649937134 45 D-28    /
      DATA BTH1CS( 37) / +.1785000847 8186376177 8586197964 17 D-28    /
      DATA BTH1CS( 38) / -.6840436100 4510395406 2152235667 46 D-29    /
      DATA BTH1CS( 39) / +.2656602867 1720419358 2934226722 12 D-29    /
      DATA BTH1CS( 40) / -.1045040252 7914452917 7141614846 70 D-29    /
      DATA BTH1CS( 41) / +.4161829082 5377144306 8619171970 64 D-30    /
      DATA BTH1CS( 42) / -.1677163920 3643714856 5013478828 87 D-30    /
      DATA BTH1CS( 43) / +.6836199777 6664389173 5359280285 28 D-31    /
      DATA BTH1CS( 44) / -.2817224786 1233641166 7395746228 10 D-31    /
C
C-------------------------------------------------------------------
C
C Series for BT12       Used in double precision version only
C                                        with weighted error   3.33E-32
C                                         log weighted error  31.48
C                               significant figures required  31.05
C                                    decimal places required  32.27
C
      DATA BT12CS(  1) / +.7382386012 8742974662 6208397927 64 D+0     /
      DATA BT12CS(  2) / -.3336111317 4483906384 4701476811 89 D-2     /
      DATA BT12CS(  3) / +.6146345488 8046964698 5148994201 86 D-4     /
      DATA BT12CS(  4) / -.2402458516 1602374264 9776354695 68 D-5     /
      DATA BT12CS(  5) / +.1466355557 7509746153 2105919972 04 D-6     /
      DATA BT12CS(  6) / -.1184191730 5589180567 0051475049 83 D-7     /
      DATA BT12CS(  7) / +.1157419896 3919197052 1254663030 55 D-8     /
      DATA BT12CS(  8) / -.1300116112 9439187449 3660077945 71 D-9     /
      DATA BT12CS(  9) / +.1624539114 1361731937 7421662736 67 D-10    /
      DATA BT12CS( 10) / -.2208963682 1403188752 1554417701 28 D-11    /
      DATA BT12CS( 11) / +.3218030425 8553177090 4743586537 78 D-12    /
      DATA BT12CS( 12) / -.4965314793 2768480785 5520211353 81 D-13    /
      DATA BT12CS( 13) / +.8043890043 2847825985 5588826393 17 D-14    /
      DATA BT12CS( 14) / -.1358912131 0161291384 6947126822 82 D-14    /
      DATA BT12CS( 15) / +.2381050439 7147214869 6765296059 73 D-15    /
      DATA BT12CS( 16) / -.4308146636 3849106724 4712414207 99 D-16    /
      DATA BT12CS( 17) / +.8020254403 2771002434 9935125504 00 D-17    /
      DATA BT12CS( 18) / -.1531631064 2462311864 2300274687 99 D-17    /
      DATA BT12CS( 19) / +.2992860635 2715568924 0730405546 66 D-18    /
      DATA BT12CS( 20) / -.5970996465 8085443393 8156366506 66 D-19    /
      DATA BT12CS( 21) / +.1214028966 9415185024 1608526506 66 D-19    /
      DATA BT12CS( 22) / -.2511511469 6612948901 0069777066 66 D-20    /
      DATA BT12CS( 23) / +.5279056717 0328744850 7383807999 99 D-21    /
      DATA BT12CS( 24) / -.1126050922 7550498324 3611613866 66 D-21    /
      DATA BT12CS( 25) / +.2434827735 9576326659 6634624000 00 D-22    /
      DATA BT12CS( 26) / -.5331726123 6931800130 0384426666 66 D-23    /
      DATA BT12CS( 27) / +.1181361505 9707121039 2059903999 99 D-23    /
      DATA BT12CS( 28) / -.2646536828 3353523514 8567893333 33 D-24    /
      DATA BT12CS( 29) / +.5990339404 1361503945 5778133333 33 D-25    /
      DATA BT12CS( 30) / -.1369085463 0829503109 1363839999 99 D-25    /
      DATA BT12CS( 31) / +.3157679015 4380228326 4136533333 33 D-26    /
      DATA BT12CS( 32) / -.7345791508 2084356491 4005333333 33 D-27    /
      DATA BT12CS( 33) / +.1722808148 0722747930 7059200000 00 D-27    /
      DATA BT12CS( 34) / -.4071690796 1286507941 0688000000 00 D-28    /
      DATA BT12CS( 35) / +.9693474513 6779622700 3733333333 33 D-29    /
      DATA BT12CS( 36) / -.2323763633 7765716765 3546666666 66 D-29    /
      DATA BT12CS( 37) / +.5607451067 3522029406 8906666666 66 D-30    /
      DATA BT12CS( 38) / -.1361646539 1539005860 5226666666 66 D-30    /
      DATA BT12CS( 39) / +.3326310923 3894654388 9066666666 66 D-31    /
C
C-------------------------------------------------------------------
C
C Series for BJ1        on the interval  0.          to  1.60000D+01
C                                        with weighted error   1.16E-33
C                                         log weighted error  32.93
C                               significant figures required  32.36
C                                    decimal places required  33.57
C
      DATA BJ1 CS(  1) / -.1172614151 3332786560 6240574524 003 D+0    /
      DATA BJ1 CS(  2) / -.2536152183 0790639562 3030884554 698 D+0    /
      DATA BJ1 CS(  3) / +.5012708098 4469568505 3656363203 743 D-1    /
      DATA BJ1 CS(  4) / -.4631514809 6250819184 2619728789 772 D-2    /
      DATA BJ1 CS(  5) / +.2479962294 1591402453 9124064592 364 D-3    /
      DATA BJ1 CS(  6) / -.8678948686 2788258452 1246435176 416 D-5    /
      DATA BJ1 CS(  7) / +.2142939171 4379369150 2766250991 292 D-6    /
      DATA BJ1 CS(  8) / -.3936093079 1831797922 9322764073 061 D-8    /
      DATA BJ1 CS(  9) / +.5591182317 9468800401 8248059864 032 D-10   /
      DATA BJ1 CS( 10) / -.6327616404 6613930247 7695274014 880 D-12   /
      DATA BJ1 CS( 11) / +.5840991610 8572470032 6945563268 266 D-14   /
      DATA BJ1 CS( 12) / -.4482533818 7012581903 9135059199 999 D-16   /
      DATA BJ1 CS( 13) / +.2905384492 6250246630 6018688000 000 D-18   /
      DATA BJ1 CS( 14) / -.1611732197 8414416541 2118186666 666 D-20   /
      DATA BJ1 CS( 15) / +.7739478819 3927463729 8346666666 666 D-23   /
      DATA BJ1 CS( 16) / -.3248693782 1119984114 3466666666 666 D-25   /
      DATA BJ1 CS( 17) / +.1202237677 2274102272 0000000000 000 D-27   /
      DATA BJ1 CS( 18) / -.3952012212 6513493333 3333333333 333 D-30   /
      DATA BJ1 CS( 19) / +.1161678082 2664533333 3333333333 333 D-32   /
C
C-------------------------------------------------------------------
C
      DATA NTY1 / 0 /
C
C***FIRST EXECUTABLE STATEMENT  DWY1
C
      IF (M .LE. 0)  GO TO 910
C
      IF (NTY1 .EQ. 0) THEN
         EPMACH = D1MACH(3)
         EPS = 0.10D0*EPMACH
         NTY1 = IDWCS(BY1CS, LBY1, EPS)
         NTM1 = IDWCS(BM1CS, LBM1, EPS)
         NTM12 = IDWCS(BM12CS, LBM12, EPS)
         NTTH1 = IDWCS(BTH1CS, LBTH1, EPS)
         NTT12 = IDWCS(BT12CS, LBT12, EPS)
         NTJ1 = IDWCS(BJ1CS, LBJ1, EPS)
         XMIN = EXP(MAX( LOG(D1MACH(1)), -LOG(D1MACH(2))) + 0.010D0)
         XSML = SQRT(4.0D0*EPMACH)
         XMAX = 1.0/D1MACH(4)
         PI4 = ATAN(1.0D0)
         TPI4 = 3.0D0*PI4
         TWODPI = 0.50D0/PI4
      ENDIF
C
      CALL DWNLE(M,X,0.0D0,KEY)
      IF (KEY .NE. 0) GO TO 920
C
      CALL DWNLE(M,X,XMIN,KEY)
      IF (KEY .NE. 0) GO TO 930
C
C  ----------------
C  CASE  X .LE. 4.0
C  ----------------
C
      CALL DWLE(M,X,4.0D0,N,INDX)
      IF (N .GT. 0) THEN
         CALL DWGTHR(N,X,INDX,XCMP)
C
C        ... COMPUTE J1(X) ... RESULT IN ZCMP
C
         DO 20 J=1,N
            TCMP(J) = 0.125E0*XCMP(J)**2 - 1.0D0
   20    CONTINUE
         CALL DWCS(N,TCMP,BJ1CS,NTJ1,ZCMP,B0,B1,B2)
         DO 30 J=1,N
            ZCMP(J) = XCMP(J)*(0.250D0 + ZCMP(J))
   30    CONTINUE
C
         CALL DWCS(N,TCMP,BY1CS,NTY1,YCMP,B0,B1,B2)
         DO 40 J=1,N
            ZCMP(J) = TWODPI*LOG(0.50D0*XCMP(J))*ZCMP(J) +
     +                (0.50D0 + YCMP(J))/XCMP(J)
   40    CONTINUE
         CALL DWSCTR(N,ZCMP,INDX,F)
      ENDIF
C
C  ----------------
C  CASE  X .GT. 4.0
C  ----------------
C
      CALL DWGTLE(M,X,4.0D0,8.0D0,NA,INDX)
      JH = NA + 1
      CALL DWGT(M,X,8.0D0,NB,INDX(JH))
      N = NA + NB
      IF (N .GT. 0) THEN
         CALL DWGTHR(N,X,INDX,XCMP)
         C1 = 128.0D0/3.0D0
         C2 = 5.0D0/3.0D0
         DO 50 J=1,NA
            ZCMP(J) = C1/XCMP(J)**2 - C2
  50     CONTINUE
         CALL DWCS(NA,ZCMP,BM1CS,NTM1,YCMP,B0,B1,B2)
         CALL DWCS(NA,ZCMP,BT12CS,NTT12,ZCMP,B0,B1,B2)
         DO 60 J=JH,N
            ZCMP(J) = 128.0D0/XCMP(J)**2 - 1.0D0
  60     CONTINUE
         CALL DWCS(NB,ZCMP(JH),BM12CS,NTM12,YCMP(JH),B0,B1,B2)
         CALL DWCS(NB,ZCMP(JH),BTH1CS,NTTH1,ZCMP(JH),B0,B1,B2)
         DO 70 J=1,N
            YCMP(J) = (0.750D0 + YCMP(J)) / SQRT(XCMP(J))
            ZCMP(J) = (XCMP(J) - TPI4) + ZCMP(J) / XCMP(J)
            ZCMP(J) = YCMP(J) * SIN(ZCMP(J))
  70     CONTINUE
         CALL DWSCTR(N,ZCMP,INDX,F)
      ENDIF
C
C  -----
C  EXITS
C  -----
C
C     ... NORMAL
C
      INFO = 0
      GO TO 999
C
C     ... M .LE. 0
C
  910 CONTINUE
      INFO = 1
      GO TO 999
C
C     ... X IS ZERO OR NEGATIVE
C
  920 CONTINUE
      INFO = 2
      INDX(1) = KEY
      GO TO 999
C
C     ... X SO SMALL THAT Y1 OVERFLOWS
C
  930 CONTINUE
      INFO = 3
      INDX(1) = KEY
      GO TO 999
C
  999 CONTINUE
      RETURN
      END
      SUBROUTINE DWCS (M, X, CS, N, F, B0, B1, B2)
C***BEGIN PROLOGUE  DWCS
C***PURPOSE  Evaluate a Chebyshev series for a vector of real arguments
C***LIBRARY   VFNLIB
C***CATEGORY  C3A2
C***TYPE      DOUBLE PRECISION (DWCS-S, DWCS-D)
C***KEYWORDS  CHEBYSHEV SERIES EVALUATION, VECTORIZED
C***AUTHOR  SAUNDERS, B. V,. (NIST)
C             COMPUTING AND APPLIED MATHEMATICS LABORATORY
C             NATIONAL INSTITUTE OF STANDARDS AND TECHNOLOGY
C             GAITHERSBURG, MD 20899
C           BOISVERT, R. F., (NIST)
C             COMPUTING AND APPLIED MATHEMATICS LABORATORY
C             NATIONAL INSTITUTE OF STANDARDS AND TECHNOLOGY
C             GAITHERSBURG, MD 20899
C***DESCRIPTION
C
C   DWCS evaluates a given Chebyshev series for a vector of real
C   arguments.
C
C
C   P A R A M E T E R S
C
C   M   (Input) Integer
C       The number of arguments at which the series is to be
C       evaluated.
C
C   X   (Input) Double precision array of length M
C       The arguments at which the function is to be evaluated are
C       stored in X(1) to X(M) in any order.
C
C   CS  (Input) Double precision array of length .GE. N
C       The N coefficients of the Chebyshev series. (Note -- only half
C       the first coefficient is used in summing the series.)
C
C   N   (Input) Integer (0 .LE. N .LE. 1000)
C       The number of terms in the Chebyshev series.
C
C   F   (Output) Double precision array of length M
C       F(i) contains the value of the series at X(i), i=1,..,M.
C
C   B0  (Work) Double precision array of length M
C
C   B1  (Work) Double precision array of length M
C
C   B2  (Work) Double precision array of length M
C
C
C *********************************************************************
C   This routine is a modification of the function CSEVL  developed by
C   W. Fullerton of LANL.
C *********************************************************************
C
C***REFERENCES  (NONE)
C***ROUTINES CALLED  WFERR
C***REVISION HISTORY  (YYMMDD)
C   910226  DATE WRITTEN
C***END PROLOGUE  DWCS
C
C     ----------
C     PARAMETERS
C     ----------
C
      INTEGER N, M
      DOUBLE PRECISION X, CS, F, B0, B1, B2
C
      DIMENSION B0(M), B1(M), B2(M), CS(N), F(M), X(M)
C
C     ---------------
C     LOCAL VARIABLES
C     ---------------
C
      INTEGER I, IMOD, J, NI
      DOUBLE PRECISION CS1, CS2, CSNI, CSNI1, CSNI2
C
C***FIRST EXECUTABLE STATEMENT  DWCS
C
      IF (N .LT.    0)  CALL WFERR('DWCS','NUMBER OF TERMS LT 0',2)
      IF (N .GT. 1000)  CALL WFERR('DWCS','NUMBER OF TERMS GT 1000',2)
C
C     ... INITIALIZATION
C
      DO 10 I=1,M
         B1(I) = 0.0D0
         B2(I) = 0.0D0
         F(I)  = 2.0D0*X(I)
  10  CONTINUE
C
C     ... THREE-TERM RECURSION
C         (DO THREE STEPS AT A TIME TO AVOID VECTOR COPIES)
C
      IMOD = MOD(N,3)
      DO 30 I=1,N-IMOD,3
         NI = N + 1 - I
         CSNI = CS(NI)
         CSNI1 = CS(NI-1)
         CSNI2 = CS(NI-2)
         DO 20 J=1,M
            B0(J) = ( F(J)*B1(J)-B2(J) ) + CSNI
            B2(J) = ( F(J)*B0(J)-B1(J) ) + CSNI1
            B1(J) = ( F(J)*B2(J)-B0(J) ) + CSNI2
  20     CONTINUE
  30  CONTINUE
C
C     ... LAST STEP
C         (CLEANUP FOR CASE N NOT DIVISIBLE BY 3)
C
      IF (IMOD .EQ. 0) THEN
         DO 40 J=1,M
            F(J) = 0.50D0*(B1(J) - B0(J))
  40     CONTINUE
      ELSEIF (IMOD .EQ. 1) THEN
         CS1 = CS(1)
         DO 50 J=1,M
            B0(J) = ( F(J)*B1(J) - B2(J) ) + CS1
            F(J)  = 0.50D0*(B0(J) - B2(J))
  50     CONTINUE
      ELSE
         CS1 = CS(1)
         CS2 = CS(2)
         DO 60 J=1,M
            B0(J) = ( F(J)*B1(J) - B2(J) ) + CS2
            B2(J) = ( F(J)*B0(J) - B1(J) ) + CS1
            F(J)  = 0.50D0*(B2(J) - B1(J))
  60     CONTINUE
      ENDIF
C
      RETURN
      END
      INTEGER FUNCTION IDWCS (OS, NOS, ETA)
C***BEGIN PROLOGUE  IDWCS
C***PURPOSE  Determines the number of terms of an orthogonal series
C            required to meet a specified error tolerance.
C***LIBRARY   VFNLIB
C***CATEGORY  C3A2
C***TYPE      DOUBLE PRECISION (IDWCS-S, IDWCS-D)
C***KEYWORDS  CHEBYSHEV SERIES, INITIALIZATION
C***AUTHOR  BOISVERT, R. F., (NIST)
C             COMPUTING AND APPLIED MATHEMATICS LABORATORY
C             NATIONAL INSTITUTE OF STANDARDS AND TECHNOLOGY
C             GAITHERSBURG, MD 20899
C***DESCRIPTION
C
C   IDWCS returns the number of terms of the Chebyshev series OS
C   required to insure the error in evaluating it is no larger than
C   ETA.
C
C   Ordinarily, ETA is chosen to be one-tenth machine precision.
C
C
C   P A R A M E T E R S
C
C   OS     (Input) Double precision array of length .GE. NOS
C          Coefficients of the NOS-term orthogonal series.
C
C   NOS    (Input) Integer (.GE. 1)
C          Number of terms in the orthogonal series.
C
C   ETA    (Input) Double precision
C          Requested accuracy of the series.
C
C
C *********************************************************************
C   This routine is a modification of the function INITS developed by
C   W. Fullerton of LANL.
C *********************************************************************
C
C***REFERENCES  (NONE)
C***ROUTINES CALLED  WFERR
C***REVISION HISTORY  (YYMMDD)
C   910226  DATE WRITTEN
C***END PROLOGUE  IDWCS
C
C     ----------
C     PARAMETERS
C     ----------
C
      INTEGER NOS
      DOUBLE PRECISION ETA, OS
C
      DIMENSION OS(NOS)
C
C     ---------------
C     LOCAL VARIABLES
C     ---------------
C
      INTEGER I, II
      DOUBLE PRECISION ERR
C
C***FIRST EXECUTABLE STATEMENT  IDWCS
C
      IF (NOS .LT. 1)
     +   CALL WFERR('IDWCS','NUMBER OF COEFFICIENTS LT 1',2)
C
      ERR = 0.0D0
      DO 10 II=NOS,1,-1
         I = II
         ERR = ERR + ABS(OS(I))
         IF (ERR .GT. ETA) GO TO 20
 10   CONTINUE
C
 20   CONTINUE
      IF (I .EQ. NOS)
     +   CALL WFERR('IDWCS','TOO MUCH ACCURACY REQUESTED',2)
C
      IDWCS = I
C
      RETURN
      END
      SUBROUTINE DWNLE (M, X, SCALR, KEY)
C***BEGIN PROLOGUE  DWNLE
C***PURPOSE  Determines whether elements of a given vector are less
C            than or equal to a given scalar
C***LIBRARY   VFNLIB
C***TYPE      DOUBLE PRECISION (DWNLE-S, DDWNLE-D)
C***AUTHOR  SAUNDERS, B. V., (NIST)
C             COMPUTING AND APPLIED MATHEMATICS LABORATORY
C             NATIONAL INSTITUTE OF STANDARDS AND TECHNOLOGY
C             GAITHERSBURG, MD 20899
C           BOISVERT, R. F., (NIST)
C             COMPUTING AND APPLIED MATHEMATICS LABORATORY
C             NATIONAL INSTITUTE OF STANDARDS AND TECHNOLOGY
C             GAITHERSBURG, MD 20899
C***DESCRIPTION
C
C   DWNLE determines whether elements of a given vector are less than
C   or equal to a given scalar.
C
C
C   P A R A M E T E R S
C
C   M      (Input) Integer (M .GT. 0)
C          The number of elements of the input vector.
C
C   X      (Input) Double precision array of length M
C          The input vector.
C
C   SCALR  (Input) Double precision
C          The scalar used for comparison with vector elements.
C
C   KEY    (Output) Integer
C          If KEY=0 then no vector elements satisfy X(i).LE.SCALR.
C          If KEY>0 then some vector elements satisfy X(i).LE.SCALR;
C          the first to do so is X(KEY).
C
C
C   CAUTION : constraints on the input variables are not checked by
C             this routine.
C
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   910226  DATE WRITTEN
C***END PROLOGUE  DWNLE
C
C     ----------
C     PARAMETERS
C     ----------
C
      INTEGER KEY, M
      DOUBLE PRECISION SCALR, X
C
      DIMENSION X(M)
C
C     ---------------
C     LOCAL VARIABLES
C     ---------------
C
      INTEGER I
C
C***FIRST EXECUTABLE STATEMENT  DWNLE
C
C     ... QUICK CHECK (VECTORIZABLE)
C
      KEY = 0
      DO 10 I=1,M
         IF (X(I) .LE. SCALR) KEY = 1
   10 CONTINUE
C
C     ... IF CHECK FAILED THEN FIND INDEX OF FIRST FAILURE
C
      IF (KEY .NE. 0) THEN
         KEY = 0
         DO 20 I=1,M
            IF (X(I) .LE. SCALR) THEN
               KEY = I
               GO TO 30
            ENDIF
   20    CONTINUE
      ENDIF
C
   30 CONTINUE
      RETURN
      END

      SUBROUTINE DWNGT (M, X, SCALR, KEY)
C***BEGIN PROLOGUE  DWNGT
C***PURPOSE  Determines whether elements of a given vector are greater
C            than a given scalar
C***LIBRARY   VFNLIB
C***TYPE      DOUBLE PRECISION (DWNGT-S, DWNGT-D)
C***AUTHOR  SAUNDERS, B. V., (NIST)
C             COMPUTING AND APPLIED MATHEMATICS LABORATORY
C             NATIONAL INSTITUTE OF STANDARDS AND TECHNOLOGY
C             GAITHERSBURG, MD 20899
C           BOISVERT, R. F., (NIST)
C             COMPUTING AND APPLIED MATHEMATICS LABORATORY
C             NATIONAL INSTITUTE OF STANDARDS AND TECHNOLOGY
C             GAITHERSBURG, MD 20899
C***DESCRIPTION
C
C   DWNGT determines whether elements of a given vector are greater
C   than a given scalar.
C
C
C   P A R A M E T E R S
C
C   M      (Input) Integer (M .GT. 0)
C          The number of elements of the input vector.
C
C   X      (Input) Double precision array of length M
C          The input vector.
C
C   SCALR  (Input) Double precision
C          The scalar used for comparison with vector elements.
C
C   KEY    (Output) Integer
C          If KEY=0 then no vector elements satisfy X(i).GT.SCALR.
C          If KEY>0 then some vector elements satisfy X(i).GT.SCALR;
C          the first to do so is X(KEY).
C
C
C   CAUTION : constraints on the input variables are not checked by
C             this routine.
C
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   910226  DATE WRITTEN
C***END PROLOGUE  DWNGT
C
C     ----------
C     PARAMETERS
C     ----------
C
      INTEGER KEY, M
      DOUBLE PRECISION SCALR, X
C
      DIMENSION X(M)
C
C     ---------------
C     LOCAL VARIABLES
C     ---------------
C
      INTEGER I
C
C***FIRST EXECUTABLE STATEMENT  DWNGT
C
C     ... QUICK CHECK (VECTORIZABLE)
C
      KEY = 0
      DO 10 I=1,M
         IF (X(I) .GT. SCALR) KEY = 1
   10 CONTINUE
C
C     ... IF CHECK FAILED THEN FIND INDEX OF FIRST FAILURE
C
      IF (KEY .NE. 0) THEN
         KEY = 0
         DO 20 I=1,M
            IF (X(I) .GT. SCALR) THEN
               KEY = I
               GO TO 30
            ENDIF
   20    CONTINUE
      ENDIF
C
   30 CONTINUE
      RETURN
      END

      SUBROUTINE DWGTLE (M, X, SCALR1, SCALR2, N, INDX)
C***BEGIN PROLOGUE  DWGTLE
C***PURPOSE  Builds an array of indices of vector elements that are
C            greater than a given scalar and less than or equal to a
C            second given scalar
C***LIBRARY   VFNLIB
C***TYPE      DOUBLE PRECISION (DWGTLE-S, DDWGTLE-D)
C***AUTHOR  SAUNDERS, B. V., (NIST)
C             COMPUTING AND APPLIED MATHEMATICS LABORATORY
C             NATIONAL INSTITUTE OF STANDARDS AND TECHNOLOGY
C             GAITHERSBURG, MD 20899
C           BOISVERT, R. F., (NIST)
C             COMPUTING AND APPLIED MATHEMATICS LABORATORY
C             NATIONAL INSTITUTE OF STANDARDS AND TECHNOLOGY
C             GAITHERSBURG, MD 20899
C***DESCRIPTION
C
C   DWGTLE builds an array of indices of vector elements that are
C   greater than a given scalar and less than or equal to a second
C   given scalar.
C
C
C   P A R A M E T E R S
C
C   M      (Input) Integer (M .GT. 0)
C          The number of elements of the input vector.
C
C   X      (Input) Double precision array of length M
C          The input vector.
C
C   SCALR1,
C   SCALR2 (Input) Double precision
C          The scalars used for comparison with vector elements.
C          Indices i are selected if SCALR1 .LT. X(i) .LE. SCALR2.
C
C   N      (Output) Integer
C          The number of elements of X that satisfy relationship with
C          scalars.
C
C   INDX   (Output) Integer array of length N
C          Array containing indices of vector elements that satisfy
C          relationship with scalars.
C
C
C   CAUTION : constraints on the input variables are not checked by
C             this routine.
C
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   910226  DATE WRITTEN
C***END PROLOGUE  DWGTLE
C
C     ----------
C     PARAMETERS
C     ----------
C
      INTEGER INDX, M, N
      DOUBLE PRECISION SCALR1, SCALR2, X
C
      DIMENSION X(M), INDX(*)
C
C     ---------------
C     LOCAL VARIABLES
C     ---------------
C
      INTEGER I
      DOUBLE PRECISION ELEMT
C
C***FIRST EXECUTABLE STATEMENT DWGTLE
C
      N = 0
C
      DO 10 I=1,M
         ELEMT = X(I)
         IF (ELEMT .GT. SCALR1 .AND. ELEMT .LE. SCALR2) THEN
            N = N + 1
            INDX(N) = I
         ENDIF
  10  CONTINUE
C
      RETURN
      END

      SUBROUTINE DWGT (M, X, SCALR, N, INDX)
C***BEGIN PROLOGUE  DWGT
C***PURPOSE  Builds an array of indices of vector elements that are
C            greater than a specified scalar
C***LIBRARY   VFNLIB
C***TYPE      DOUBLE PRECISION (DWGT-S, DWGT-D)
C***AUTHOR  SAUNDERS, B. V., (NIST)
C             COMPUTING AND APPLIED MATHEMATICS LABORATORY
C             NATIONAL INSTITUTE OF STANDARDS AND TECHNOLOGY
C             GAITHERSBURG, MD 20899
C           BOISVERT, R. F., (NIST)
C             COMPUTING AND APPLIED MATHEMATICS LABORATORY
C             NATIONAL INSTITUTE OF STANDARDS AND TECHNOLOGY
C             GAITHERSBURG, MD 20899
C***DESCRIPTION
C
C   DWGT builds an array of indices of vector elements that are
C   greater than a specified scalar.
C
C
C   P A R A M E T E R S
C
C   M      (Input) Integer (M .GT. 0)
C          The number of elements of the input vector.
C
C   X      (Input) Double precision array of length M
C          The input vector.
C
C   SCALR  (Input) Double precision
C          The scalar used for comparison with vector elements.
C          Indices i are selected if X(i) .GT. SCALR.
C
C   N      (Output) Integer
C          The number of elements of X that satisfy relationship with
C          scalar.
C
C   INDX   (Output) Integer array of length N
C          Array containing indices of vector elements that satisfy
C          relationship with scalar.
C
C
C   CAUTION : constraints on the input variables are not checked by
C             this routine.
C
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   910226  DATE WRITTEN
C***END PROLOGUE  DWGT
C
C     ----------
C     PARAMETERS
C     ----------
C
      INTEGER INDX, M, N
      DOUBLE PRECISION SCALR, X
C
      DIMENSION X(M), INDX(*)
C
C     ---------------
C     LOCAL VARIABLES
C     ---------------
C
      INTEGER I
C
C***FIRST EXECUTABLE STATEMENT DWGT
C
      N = 0
C
      DO 10 I=1,M
         IF (X(I) .GT. SCALR) THEN
            N = N + 1
            INDX(N) = I
         ENDIF
  10  CONTINUE
C
      RETURN
      END

      SUBROUTINE DWLE (M, X, SCALR, N, INDX)
C***BEGIN PROLOGUE  DWLE
C***PURPOSE  Builds an array of indices of vector elements that are
C            less than  or equal to a given scalar
C***LIBRARY   VFNLIB
C***TYPE      DOUBLE PRECISION (DWLE-S, DWLE-D)
C***AUTHOR  SAUNDERS, B. V., (NIST)
C             COMPUTING AND APPLIED MATHEMATICS LABORATORY
C             NATIONAL INSTITUTE OF STANDARDS AND TECHNOLOGY
C             GAITHERSBURG, MD 20899
C           BOISVERT, R. F., (NIST)
C             COMPUTING AND APPLIED MATHEMATICS LABORATORY
C             NATIONAL INSTITUTE OF STANDARDS AND TECHNOLOGY
C             GAITHERSBURG, MD 20899
C***DESCRIPTION
C
C   DWLE builds an array of indices of vector elements that are
C   less than or equal to a given scalar.
C
C
C   P A R A M E T E R S
C
C   M      (Input) Integer (M .GT. 0)
C          The number of elements of the input vector.
C
C   X      (Input) Double precision array of length M
C          The input vector.
C
C   SCALR  (Input) Double precision
C          The scalar used for comparison with vector elements.
C          Indices i are selected if X(i) .LE. SCALR.
C
C   N      (Output) Integer
C          The number of elements of X that satisfy relationship with
C          scalars.
C
C   INDX   (Output) Integer array of length N
C          Array containing indices of vector elements that satisfy
C          relationship with scalars.
C
C
C   CAUTION : constraints on the input variables are not checked by
C             this routine.
C
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   910226  DATE WRITTEN
C***END PROLOGUE  DWLE
C
C     ----------
C     PARAMETERS
C     ----------
C
      INTEGER INDX, M, N
      DOUBLE PRECISION SCALR, X
C
      DIMENSION X(M), INDX(*)
C
C     ---------------
C     LOCAL VARIABLES
C     ---------------
C
      INTEGER I
C
C***FIRST EXECUTABLE STATEMENT WGE
C
      N = 0
C
      DO 10 I = 1,M
         IF (X(I) .LE. SCALR) THEN
            N = N + 1
            INDX(N) = I
         ENDIF
  10  CONTINUE
C
      RETURN
      END

      SUBROUTINE DWGTHR (N, X, INDX, Y)
C***BEGIN PROLOGUE  DWGTHR
C***PURPOSE  Performs a gather by index on a given vector
C***LIBRARY   VFNLIB
C***TYPE      DOUBLE PRECISION (WGTHR-S, DWGTHR-D)
C***AUTHOR  SAUNDERS, B. V., (NIST)
C             COMPUTING AND APPLIED MATHEMATICS LABORATORY
C             NATIONAL INSTITUTE OF STANDARDS AND TECHNOLOGY
C             GAITHERSBURG, MD 20899
C           BOISVERT, R. F., (NIST)
C             COMPUTING AND APPLIED MATHEMATICS LABORATORY
C             NATIONAL INSTITUTE OF STANDARDS AND TECHNOLOGY
C             GAITHERSBURG, MD 20899
C***DESCRIPTION
C
C   DWGTHR performs a gather by index on a given vector based on the
C   indices provided in an array: Y(i) = X(INDX(i)), i=1,...,N.
C
C
C   P A R A M E T E R S
C
C   N      (Input) Integer (N .GT. 0)
C          The number of elements to be gathered from the input vector
C          X.
C
C   X      (Input) Double precision array of length M
C          The input vector.
C
C   INDX   (Input) Integer array of length N
C          Array specifying indices of elements to be gathered from
C          the input vector.
C
C   Y      (Output) Double precision array of length N
C          The array containing the compressed vector once the gather is
C          completed.  Y(i) = X(INDX(i)), i=1,...,N.
C
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   910226  DATE WRITTEN
C***END PROLOGUE  DWGTHR
C
C     ----------
C     PARAMETERS
C     ----------
C
      INTEGER INDX, N
      DOUBLE PRECISION X, Y
C
      DIMENSION X(*), Y(N), INDX(N)
C
C     ---------------
C     LOCAL VARIABLES
C     ---------------
C
      INTEGER J
C
C***FIRST EXECUTABLE STATEMENT DWGTHR
C
      IF (N .LE. 0) RETURN
C
      DO 10 J=1,N
         Y(J) = X(INDX(J))
  10  CONTINUE
C
      RETURN
      END

      SUBROUTINE DWSCTR (N, Y, INDX, X)
C***BEGIN PROLOGUE  DWSCTR
C***PURPOSE  Performs a scatter by index on a given vector
C***LIBRARY   VFNLIB
C***TYPE      DOUBLE PRECISION (WSCTR-S, DWSCTR-D)
C***AUTHOR  SAUNDERS, B. V., (NIST)
C             COMPUTING AND APPLIED MATHEMATICS LABORATORY
C             NATIONAL INSTITUTE OF STANDARDS AND TECHNOLOGY
C             GAITHERSBURG, MD 20899
C           BOISVERT, R. F., (NIST)
C             COMPUTING AND APPLIED MATHEMATICS LABORATORY
C             NATIONAL INSTITUTE OF STANDARDS AND TECHNOLOGY
C             GAITHERSBURG, MD 20899
C***DESCRIPTION
C
C   DWSCTR performs a scatter by index on a given vector based on the
C   indices provided in an array: X(INDX(i)) = Y(i), i=1,...,N.
C
C
C   P A R A M E T E R S
C
C   N      (Input) Integer (N .GT. 0)
C          The number of elements in the input vector Y.
C
C   Y      (Input) Double precision array of length N
C          The compressed vector.
C
C   INDX   (Input) Integer array of length N
C          Array of indices of positions elements of compressed
C          vector will occupy when scattered into vector X.
C
C   X      (Output) Double precision array
C          The vector to receive the scattered elements of Y.
C          X(INDX(i)) = Y(i), i=1,...,N.  Only elements of X
C          whose indices are in INDX are changed.
C
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   910226  DATE WRITTEN
C***END PROLOGUE  DWSCTR
C
C
C     ----------
C     PARAMETERS
C     ----------
C
      INTEGER INDX, N
      DOUBLE PRECISION Y, X
C
      DIMENSION Y(N), X(*), INDX(N)
C
C     ---------------
C     LOCAL VARIABLES
C     ---------------
C
      INTEGER J
C
C***FIRST EXECUTABLE STATEMENT DWGTHR
C
      IF (N .LE. 0) RETURN
C
      DO 10 J=1,N
         X(INDX(J)) = Y(J)
  10  CONTINUE
C
      RETURN
      END
