c  File: sfrenl.[for|f|c]
c  Contains procedures: SFRENC(), SFRENF(), SFRENG(), SFRENS()
c  and private low-level procedure: SFREN1().
c     .  Copyright (C) 1992, California Institute of Technology.
c     .  U. S. Government sponsorship under
c     .  NASA contract NAS7-918 is acknowledged.
c>> 1996-01-08 SFRENL WV Snyder Use SCSPXX for cos(Pi/2 x**2), etc.
C>> 1995-11-03 SFRENL Krogh  Removed blanks in numbers for C conversion.
c>> 1994-11-02 SFRENL Krogh  Changes to use M77CON
c>> 1994-10-18 SFRENL WV Snyder More specializing instructions
c>> 1993-02-25 SFRENL CLL. Edited to eliminate ENTRY and EQUIVALENCE.
c>> 1992-09-15 SFRENL WV Snyder Specializing instructions
c>> 1992-04-13 SFRENL WV Snyder Declare SFRENF, SFRENG, SFRENS
c>> 1992-03-18 SFRENL WV Snyder Move declarations for coefficient arrays
c>> 1992-01-24 SFRENL WV Snyder Original code
c--S replaces "?": ?FRENC, ?FREN1, ?FRENF, ?FRENG, ?FRENS, ?FRENL
c--&               ?CSPXX, ?SNPXX
c Subprograms in this file compute the Fresnel Cosine and Sine
c integrals C(x) and S(x), and the auxiliary functions f(x) and g(x),
c for any X:
c     SFRENC(X) for Fresnel integral C(X)
c     SFRENS(X) for Fresnel integral S(x)
c     SFRENF(X) for Fresnel integral auxiliary function f(x)
c     SFRENG(X) for Fresnel integral auxiliary function g(x).
c
c Developed by W. V. Snyder, Jet Propulsion Laboratory, 24 January 1992.
c
c Ref: W. J. Cody, "Chebyshev Approximations for the Fresnel Integrals",
c Mathematics of Computation, 1968, pp 450-453 plus Microfiche Suppl.
c W. V. Snyder, "Algorithm 723: Fresnel Integrals," ACM Trans. Math.
c Softw. 19, 4 (December 1993) 452-456.
c Accuracies of highest order formulae, where E is relative error:
c
c Range           Function   -log10(E)   Function   -log10(E)
c |X|<=1.2          C(x)       16.24       S(x)       17.26
c 1.2<|X|<=1.6      C(x)       17.47       S(x)       18.66
c 1.6<|X|<=1.9      f(x)       17.13       g(x)       16.25
c 1.9<|X|<=2.4      f(x)       16.64       g(x)       15.65
c 2.4<|X|           f(x)       16.89       g(x)       15.58
c
c Refer to Cody for accuracy of other approximations.
c
c     ==================================================================
      real             function SFRENC(X)
      external SFREN1
      real             X, SFREN1
      SFRENC = SFREN1(1,X)
      return
      end
c     ==================================================================
      real             function SFRENF(X)
      external SFREN1
      real             X, SFREN1
      SFRENF = SFREN1(3,X)
      return
      end
c     ==================================================================
      real             function SFRENG(X)
      external SFREN1
      real             X, SFREN1
      SFRENG = SFREN1(4,X)
      return
      end
c     ==================================================================
      real             function SFRENS(X)
      external SFREN1
      real             X, SFREN1
      SFRENS = SFREN1(2,X)
      return
      end
c     ==================================================================
      real             function SFREN1(MODE, X)
c     MODE = 1 means compute C.
c     MODE = 2 means compute S.
c     MODE = 3 means compute F.
c     MODE = 4 means compute G.
c     ------------------------------------------------------------------
c                        Internal variables.
c
c PID2 is pi / 2.
c RPI is the reciprocal of PI.
c RPISQ is the reciprocal of PI squared.
c AX is abs(x).
c BIGX is 1/sqrt(round-off).  If X > BIGX then to the working
c         precision x**2 is an integer (which we assume to be a multiple
c         of four), so cos(pi/2 * x**2) = 1, and sin(pi/2 * x**2) = 0.
c C and S are values of C(x) and S(x), respectively.
c CX and SX are cos(pi/2 * ax**2) and sin(pi/2 * ax**2), respectively.
c F and G are used to compute f(x) and g(x) when X > 1.6.
c HAVEC, HAVEF, HAVEG, HAVES are logical variables that indicate
c         whether the values stored in C, F, G and S correspond to the
c         value stored in X.  HAVEF indicates we have both F and G when
c         XSAVE .le. 1.6, and HAVEC indicates we have both C and S when
c         XSAVE .gt. 1.6.
c LARGEF is 1/(pi * underflow).  If X > LARGEF then f ~ 0.
c LARGEG is cbrt(1/(pi**2 * underflow)).  If X > LARGEG then g ~ 0.
c LARGEX is 1/sqrt(sqrt(underflow)).  If X > LARGEX then f ~ 1/(pi * x)
c         and g ~ 1/(pi**2 * x**3).
c MODE indicates the function to be computed: 1 = C(x), 2 = S(x),
c         3 = f(x), 4 = g(x).
c NEEDC, NEEDF, NEEDG, NEEDS are arrays indexed by MODE (MODE+4 when
c         X .gt. 1.6) that indicate what functions are needed.
c WANTC indicates whether C and S must be computed from F and G.
c WANTF and WANTG indicate we computed F and G on the present call.
c XSAVE is the most recently provided value of X.
c X4 is either X ** 4 or (1.0/X) ** 4.
c     If you change the order of approximation, you must change the
c     declarations and DATA statements for the coefficient arrays,
c     and the executable statements that evaluate the approximations.
c     ------------------------------------------------------------------
      external R1MACH, SCSPXX, SSNPXX
      real             R1MACH, SCSPXX, SSNPXX
      real             PID2, RPI, RPISQ
      parameter (PID2 = 1.570796326794896619231321691639751442099e0)
      parameter (RPI = 0.3183098861837906715377675267450287240689e0)
      parameter (RPISQ = RPI * RPI)
      real             AX, BIGX, C, CX, F, G, LARGEF, LARGEG, LARGEX
      real             S, SX, X, XSAVE, X4
      integer MODE
      logical HAVEC, HAVEF, HAVEG, HAVES, WANTC, WANTF, WANTG
      logical NEEDC(8), NEEDF(8), NEEDG(8), NEEDS(8)
      save BIGX, C, F, G, LARGEF, LARGEG, LARGEX, S, XSAVE
      save HAVEC, HAVEF, HAVEG, HAVES
c++   default digits=16
c++   code for digits <= 3 is inactive
c      real             PC1(0:1), QC1(1:1)
c++   code for digits > 3 & digits <= 6 is inactive
c      real             PC1(0:2), QC1(1:2)
c++   code for digits > 6 & digits <= 11 is inactive
c      real             PC1(0:3), QC1(1:3)
c++   code for digits > 11 is active
      real             PC1(0:4), QC1(1:4)
c++   end
c++   code for digits <= 2 is inactive
c      real             PC2(0:1), QC2(1:1)
c++   code for digits > 2 & digits <= 5 is inactive
c      real             PC2(0:2), QC2(1:2)
c++   code for digits > 5 & digits <= 8 is inactive
c      real             PC2(0:3), QC2(1:3)
c++   code for digits > 8 & digits <= 12 is inactive
c      real             PC2(0:4), QC2(1:4)
c++   code for digits > 12 is active
      real             PC2(0:5), QC2(1:5)
c++   end
c++   code for digits <= 3 is inactive
c      real             PS1(0:1), QS1(1:1)
c++   code for digits > 3 & digits <= 7 is inactive
c      real             PS1(0:2), QS1(1:2)
c++   code for digits > 7 & digits <= 12 is inactive
c      real             PS1(0:3), QS1(1:3)
c++   code for digits > 12 is active
      real             PS1(0:4), QS1(1:4)
c++   end
c++   code for digits <= 2 is inactive
c      real             PS2(0:1), QS2(1:1)
c++   code for digits > 2 & digits <= 5 is inactive
c      real             PS2(0:2), QS2(1:2)
c++   code for digits > 5 & digits <= 9 is inactive
c      real             PS2(0:3), QS2(1:3)
c++   code for digits > 9 & digits <= 14 is inactive
c      real             PS2(0:4), QS2(1:4)
c++   code for digits > 14 is active
      real             PS2(0:5), QS2(1:5)
c++   end
c++   code for digits <= 5 is inactive
c      real             PF1(0:1), QF1(1:1)
c++   code for digits > 5 & digits <= 8 is inactive
c      real             PF1(0:2), QF1(1:2)
c++   code for digits > 8 & digits <= 11 is inactive
c      real             PF1(0:3), QF1(1:3)
c++   code for digits > 11 & digits <= 14 is inactive
c      real             PF1(0:4), QF1(1:4)
c++   code for digits > 14 is active
      real             PF1(0:5), QF1(1:5)
c++   end
c++   code for digits <= 5 is inactive
c      real             PF2(0:1), QF2(1:1)
c++   code for digits > 5 & digits <= 8 is inactive
c      real             PF2(0:2), QF2(1:2)
c++   code for digits > 8 & digits <= 11 is inactive
c      real             PF2(0:3), QF2(1:3)
c++   code for digits > 11 & digits <= 14 is inactive
c      real             PF2(0:4), QF2(1:4)
c++   code for digits > 14 is active
      real             PF2(0:5), QF2(1:5)
c++   end
c++   code for digits <= 3 is inactive
c      real             PF3(0:0)
c++   code for digits > 3 & digits <= 6 is inactive
c      real             PF3(0:1), QF3(1:1)
c++   code for digits > 6 & digits <= 9 is inactive
c      real             PF3(0:2), QF3(1:2)
c++   code for digits > 9 & digits <= 11 is inactive
c      real             PF3(0:3), QF3(1:3)
c++   code for digits > 11 & digits <= 13 is inactive
c      real             PF3(0:4), QF3(1:4)
c++   code for digits > 13 & digits <= 15 is inactive
c      real             PF3(0:5), QF3(1:5)
c++   code for digits > 15 is active
      real             PF3(0:6), QF3(1:6)
c++   end
c++   code for digits <= 4 is inactive
c      real             PG1(0:1), QG1(1:1)
c++   code for digits > 4 & digits <= 7 is inactive
c      real             PG1(0:2), QG1(1:2)
c++   code for digits > 7 & digits <= 10 is inactive
c      real             PG1(0:3), QG1(1:3)
c++   code for digits > 10 & digits <= 13 is inactive
c     real             PG1(0:4), QG1(1:4)
c++   code for digits > 13 is active
      real             PG1(0:5), QG1(1:5)
c++   end
c++   code for digits <= 4 is inactive
c      real             PG2(0:1), QG2(1:1)
c++   code for digits > 4 & digits <= 7 is inactive
c      real             PG2(0:2), QG2(1:2)
c++   code for digits > 7 & digits <= 10 is inactive
c      real             PG2(0:3), QG2(1:3)
c++   code for digits > 10 & digits <= 13 is inactive
c      real             PG2(0:4), QG2(1:4)
c++   code for digits > 13 is active
      real             PG2(0:5), QG2(1:5)
c++   end
c++   code for digits <= 3 is inactive
c      real             PG3(0:0)
c++   code for digits > 3 & digits <= 5 is inactive
c      real             PG3(0:1), QG3(r:1)
c++   code for digits > 5 & digits <= 8 is inactive
c      real             PG3(0:2), QG3(1:2)
c++   code for digits > 8 & digits <= 10 is inactive
c      real             PG3(0:3), QG3(1:3)
c++   code for digits > 10 & digits <= 12 is inactive
c      real             PG3(0:4), QG3(1:4)
c++   code for digits = 13 is inactive
c      real             PG3(0:5), QG3(1:5)
c++   code for digits > 13 is active
      real             PG3(0:6), QG3(1:6)
c++   end
c
      data BIGX /-1.0e0/
      data C /0.0e0/, F /0.5e0/, G /0.5e0/, S /0.0e0/, XSAVE /0.0e0/
      data HAVEC/.TRUE./, HAVEF/.TRUE./, HAVEG/.TRUE./, HAVES/.TRUE./
c        C(x)    S(x)    f(x)    g(x)    C(x)    S(x)    f(x)    g(x)
      data NEEDC
     1 /.TRUE., .FALSE.,.TRUE., .TRUE., .TRUE., .FALSE.,.FALSE.,.FALSE./
      data NEEDS
     1 /.FALSE.,.TRUE., .TRUE., .TRUE., .FALSE.,.TRUE., .FALSE.,.FALSE./
      data NEEDF
     1 /.FALSE.,.FALSE.,.TRUE., .FALSE.,.TRUE., .TRUE., .TRUE., .FALSE./
      data NEEDG
     1 /.FALSE.,.FALSE.,.FALSE.,.TRUE. ,.TRUE., .TRUE., .FALSE.,.TRUE. /
c
c     Coefficients for C(x), |X| <= 1.2
c++   code for digits <= 3 is inactive
c      data pc1 / 1.00053e0, -1.12353e-1/, qc1 /1.38937e-1/
c++   code for digits > 3 & digits <= 6 is inactive
c      data pc1 / 9.99999896e-1, -1.63090954e-1, 1.06388604e-2/
c      data qc1 / 8.36467414e-2,  3.10155884e-3/
c++   code for digits > 6 & digits <= 11 is inactive
c      data pc1 / 1.0000000000042e0,  -1.8651127631106e-1,
c     *           1.5065663274457e-2, -3.1058074693185e-4/
c      data qc1 / 6.0228833908128e-2,  1.7410300558198e-3,
c     *            .6308617899210e-5/
c++   code for digits > 11 is active
      data pc1 / 9.999999999999999421e-1 ,
     *          -1.994608988261842706e-1 ,
     *           1.761939525434914045e-2 ,
     *          -5.280796513726226960e-4 ,
     *           5.477113856826871660e-6/
      data qc1 / 4.727921120104532689e-2 ,
     *           1.099572150256418851e-3 ,
     *           1.552378852769941331e-5 ,
     *           1.189389014228757184e-7/
c++   end
c
c     Coefficients for C(x), 1.2 < |X| <= 1.6
c++   code for digits <= 2 is inactive
c      data pc2 / 1.3139e0, -5.8827e-2/, qc2 / 4.7324e-1/
c++   code for digits > 2 & digits <= 5 is inactive
c      data pc2 / 9.9790890e-1, -1.5391895e-1, 9.9817933e-3/
c      data qc2 / 8.9878647e-2,  5.6003071e-3/
c++   code for digits > 5 & digits <= 8 is inactive
c      data pc2 / 1.00000440109e0, -1.83413851908e-1,
c     *           1.45776170828e-2,-2.78980705270e-4/
c      data qc2 / 6.33347445637e-2, 2.01245738369e-3,
c     *           4.03949004646e-5/
c++   code for digits > 8 & digits <= 12 is inactive
c      data pc2 / 9.999999966496876e-1, -1.980300987022688e-1,
c     *           1.735518748450023e-2, -5.069442297935788e-4,
c     *           5.059962254678234e-6/
c      data qc2 / 4.871000308997918e-2,  1.188406974004084e-3,
c     *           1.824527635843850e-5,  1.678314257874103e-7/
c++   code for digits > 12 is active
      data pc2 / 1.00000000000111043640e0  ,
     *          -2.07073360335323894245e-1 ,
     *           1.91870279431746926505e-2 ,
     *          -6.71376034694922109230e-4 ,
     *           1.02365435056105864908e-5 ,
     *          -5.68293310121870728343e-8/
      data qc2 / 3.96667496952323433510e-2 ,
     *           7.88905245052359907842e-4 ,
     *           1.01344630866749406081e-5 ,
     *           8.77945377892369265356e-8 ,
     *           4.41701374065009620393e-10/
c++   end
c
c     Coefficients for S(x), |X| <= 1.2
c++   code for digits <= 3 is inactive
c      data ps1 / 5.23677e-1, -4.67900e-2/, qs1 / 8.81622e-2/
c++   code for digits > 3 & digits <= 7 is inactive
c      data ps1 / 5.235987665e-1,-5.837763961e-2, 2.165920196e-3/
c      data qs1 / 6.474944765e-2, 1.713287588e-3/
c++   code for digits > 7 & digits <= 12 is inactive
c      data ps1 / 5.23598775598566e-1, -6.59149581139046e-2,
c     *           3.21501649828293e-3, -4.88704436240178e-5/
c      data qs1 / 5.03546388670085e-2,  1.17835980356588e-3,
c     *           1.37089875826980e-5/
c++   code for digits > 12 is active
      data ps1 / 5.2359877559829887021e-1 ,
     *          -7.0748991514452302596e-2 ,
     *           3.8778212346368287939e-3 ,
     *          -8.4555728435277680591e-5 ,
     *           6.7174846662514086196e-7/
      data qs1 / 4.1122315114238422205e-2 ,
     *           8.1709194215213447204e-4 ,
     *           9.6269087593903403370e-6 ,
     *           5.9528122767840998345e-8/
c++   end
c
c     coefficients for S(x), 1.2 < |X| <= 1.6
c++   code for digits <= 2 is inactive
c      data ps2 / 5.4766e-1, -3.7151e-2/, qs2 / 1.4559e-1/
c++   code for digits > 2 & digits <= 5 is inactive
c      data ps2 / 5.2343994e-1, -5.4828347e-2, 1.9020881e-3/
c      data qs2 / 7.1104704e-2,  2.5596274e-3/
c++   code for digits > 5 & digits <= 9 is inactive
c      data ps2 / 5.23599040498e-1, -6.46593392426e-2,
c     *           3.08030794361e-3, -4.40800854418e-5/
c      data qs2 / 5.27536681685e-2,  1.34311026821e-3,
c     *           1.90476612849e-5/
c++   code for digits > 9 & digits <= 14 is inactive
c      data ps2 / 5.2359877543509178e-1,-7.0149076634833662e-2,
c     *           3.8031581605987038e-3,-8.0964948714408156e-5,
c     *           6.1908080210052772e-7/
c      data qs2 / 4.2268067370395487e-2, 8.7642753831073237e-4,
c     *           1.1088542889789282e-5, 7.8725829545478464e-8/
c++   code for digits > 14 is active
      data ps2 / 5.23598775598344165913e-1 ,
     *          -7.37766914010191323867e-2 ,
     *           4.30730526504366510217e-3 ,
     *          -1.09540023911434994566e-4 ,
     *           1.28531043742724820610e-6 ,
     *          -5.76765815593088804567e-9/
      data qs2 / 3.53398342767472162540e-2 ,
     *           6.18224620195473216538e-4 ,
     *           6.87086265718620117905e-6 ,
     *           5.03090581246612375866e-8 ,
     *           2.05539124458579596075e-10/
c++   end
c
c     coefficients for f(x), 1.6 < |X| <= 1.9
c++   code for digits <= 5 is inactive
c      data pf1 / 3.1803519e-1, 4.2352332e-1/, qf1 / 1.6015403e0/
c++   code for digits > 5 & digits <= 8 is inactive
c      data pf1 / 3.18285252094e-1,  2.02860275713e0,
c     *           1.08108136048e0 /
c      data qf1 / 6.67171548135e0,   4.52997553972e0 /
c++   code for digits > 8 & digits <= 11 is inactive
c      data pf1 / 3.1830646311448e-1,  4.6077249266684e0,
c     *           1.1914578438074e1,   3.5555073665830e0 /
c      data qf1 / 1.4778464938936e1,   4.0902971163705e1,
c     *           1.6130432784819e1 /
c++   code for digits > 11 & digits <= 14 is inactive
c      data pf1 / 3.1830926850490599e-1,  8.0358812280394156e0,
c     *           4.8034065557792487e1,   6.9853426160102065e1,
c     *           1.3530423554038784e1 /
c      data qf1 / 2.5549161843579529e1,   1.5761100558012250e2,
c     *           2.4956199380517229e2,   6.5563064008391554e1 /
c++   code for digits > 14 is active
      data pf1 / 3.1830975293580985290e-1 ,
     *           1.2226000551672961219e1  ,
     *           1.2924886131901657025e2  ,
     *           4.3886367156695547655e2  ,
     *           4.1466722177958961672e2  ,
     *           5.6771463664185116454e1 /
      data qf1 / 3.8713003365583442831e1  ,
     *           4.1674359830705629745e2  ,
     *           1.4740030733966610568e3  ,
     *           1.5371675584895759916e3  ,
     *           2.9113088788847831515e2 /
c++   end
c
c     coefficients for f(x), 1.9 < |X| <= 2.4
c++   code for digits <= 5 is inactive
c      data pf2 / 3.1825112e-1, 5.8395951e-1/, qf2 / 2.1243944e0 /
c++   code for digits > 5 & digits <= 8 is inactive
c      data pf2 / 3.1830699932e-1, 2.9993457087e0, 2.3956340010e0/
c      data qf2 / 9.7254517756e0,  9.4814077696e0 /
c++   code for digits > 8 & digits <= 11 is inactive
c      data pf2 / 3.1830963944521e-1,  7.0770431878327e0,
c     *           2.8756083262903e1,   1.3583685742326e1 /
c      data qf2 / 2.2536994405207e1,   9.6127568469278e1,
c     *           5.7407028004031e1 /
c++   code for digits > 11 & digits <= 14 is inactive
c      data pf2 / 3.183098568640159e-1,  1.265129469683175e1,
c     *           1.219467616498339e2,   2.910033655512762e2,
c     *           9.278397828631516e1/
c      data qf2 / 4.004915302781009e1,   3.942059697951583e2,
c     *           1.001368403495691e3,   4.142676224222433e2/
c++   code for digits > 14 is active
      data pf2 / 3.183098818220169217e-1 ,
     *           1.958839410219691002e1  ,
     *           3.398371349269842400e2  ,
     *           1.930076407867157531e3  ,
     *           3.091451615744296552e3  ,
     *           7.177032493651399590e2 /
      data qf2 / 6.184271381728873709e1  ,
     *           1.085350675006501251e3  ,
     *           6.337471558511437898e3  ,
     *           1.093342489888087888e4  ,
     *           3.361216991805511494e3 /
c++   end
c
c     coefficients for f(x), 2.4 < |X|
c++   code for digits <= 3 is inactive
c      data pf3 /-8.97969e-2/
c++   code for digits > 3 & digits <= 6 is inactive
c      data pf3 / -9.67122165e-2, -3.73565920e-1/
c      data qf3 / 7.29466572e0 /
c++   code for digits > 6 & digits <= 9 is inactive
c      data pf3 /-9.67541443648e-2, -2.26566293818e0,
c      *         -3.53429821084e0 /
c      data qf3 / 2.69596939977e1,   9.72883957469e1 /
c++   code for digits > 9 & digits <= 11 is inactive
c      data pf3 /-9.6754596389025e-2, -5.5254943840897e0,
c     *          -5.8606282086171e1,  -5.2235560918394e1 /
c      data qf3 / 6.0654483020979e1,   7.8528841711294e2,
c     *           1.8592498578831e3 /
c++   code for digits > 11 & digits <= 13 is inactive
c      data pf3 /-9.675460316952504e-2, -1.023876428129288e1,
c     *          -2.712579634037998e2,  -1.766119345282127e3,
c     *          -1.043464426656267e3 /
c      data qf3 / 1.093682244053459e2,   3.155843461920557e3,
c     *           2.625634316204420e4,   4.578252057246393e4 /
c++   code for digits > 13 & digits <= 15 is inactive
c      data pf3 /-9.67546032967090380e-2,-1.64797712841245767e1,
c     *          -8.16343401784374598e2, -1.34922028171857248e4,
c     *          -6.13547113614699772e4, -2.61294753225141779e4/
c      data qf3 / 1.73871690673649114e2,  9.01827596231524147e3,
c     *           1.65946462621853184e5,  1.00105478900791339e6,
c     *           1.37012364817225972e6 /
c++   code for digits > 15 is active
      data pf3 /-9.675460329952532343e-2 ,
     *          -2.431275407194161683e1  ,
     *          -1.947621998306889176e3  ,
     *          -6.059852197160773639e4  ,
     *          -7.076806952837779823e5  ,
     *          -2.417656749061154155e6  ,
     *          -7.834914590078317336e5 /
      data qf3 / 2.548289012949732752e2  ,
     *           2.099761536857815105e4  ,
     *           6.924122509827708985e5  ,
     *           9.178823229918143780e6  ,
     *           4.292733255630186679e7  ,
     *           4.803294784260528342e7 /
c++   end
c
c     coefficients for g(x), 1.6 < |X| <= 1.9
c++   code for digits <= 4 is inactive
c      data pg1 / 1.007011e-1, 1.457780e-1/, qg1 / 2.700253e0 /
c++   code for digits > 4 & digits <= 7 is inactive
c      data pg1 / 1.012500483e-1, 7.735207446e-1, 3.878282770e-1/
c      data qg1 / 9.101956178e0,  1.002234185e1 /
c++   code for digits > 7 & digits <= 10 is inactive
c      data pg1 / 1.013095436817e-1,  1.731937984173e0,
c     *           5.036588245265e0,   1.290707246507e0 /
c      data qg1 / 1.860077430076e1,   6.901951935725e1,
c     *           4.308918659989e1 /
c++   code for digits > 10 & digits <= 13 is inactive
c      data pg1 / 1.013188096509180e-1,  2.966220554725899e0,
c     *           2.015435299505393e1,   3.155605679387908e1,
c     *           4.916012830646366e0 /
c      data qg1 / 3.079178736724045e1,   2.362874318980473e2,
c     *           4.953861258248338e2,   2.026144493403599e2 /
c++   code for digits > 13 is active
      data pg1 / 1.013206188102747985e-1 ,
     *           4.445338275505123778e0  ,
     *           5.311228134809894481e1  ,
     *           1.991828186789025318e2  ,
     *           1.962320379716626191e2  ,
     *           2.054214324985006303e1 /
      data qg1 / 4.539250196736893605e1  ,
     *           5.835905757164290666e2  ,
     *           2.544731331818221034e3  ,
     *           3.481121478565452837e3  ,
     *           1.013794833960028555e3 /
c++   end
c
c     coefficients for g(x), 1.9 < |X| <= 2.4
c++   code for digits <= 4 is inactive
c      data pg2 / 1.011711e-1,  2.239630e-1/, qg2 / 3.609237e0 /
c++   code for digits > 4 & digits <= 7 is inactive
c      data pg2 / 1.013115463e-1,  1.182021191e0,  9.886315969e-1/
c      data qg2 / 1.317219285e1,   2.101104994e1 /
c++   code for digits > 7 & digits <= 10 is inactive
c      data pg2 / 1.013202025653e-1,  2.690767013770e0,
c     *           1.283624749271e1,   5.791337587723e0 /
c      data qg2 / 2.807457005500e1,   1.598704313522e2,
c     *           1.525630385045e2 /
c++   code for digits > 10 & digits <= 13 is inactive
c      data pg2 / 1.013210509409046e-1,  4.682769769757399e0,
c     *           5.241351613346472e1,   1.411735944550041e2,
c     *           4.019756012788710e1 /
c      data qg2 / 4.773652966704634e1,   5.802036967947208e2,
c     *           1.947179327244406e3,   1.266041236960445e3 /
c++   code for digits > 13 is active
      data pg2 / 1.01321161761804586e-1 ,
     *           7.11205001789782823e0  ,
     *           1.40959617911315524e2  ,
     *           9.08311749529593938e2  ,
     *           1.59268006085353864e3  ,
     *           3.13330163068755950e2 /
      data qg2 / 7.17128596939302198e1  ,
     *           1.49051922797329229e3  ,
     *           1.06729678030580897e4  ,
     *           2.41315567213369742e4  ,
     *           1.15149832376260604e4 /
c++   end
c
c     coefficients for g(x), 2.4 < |X|
c++   code for digits <= 3 is inactive
c      data pg3 /-1.3549e-1/
c++   code for digits > 3 & digits <= 5 is inactive
c      data pg3 /-1.5382824e-1, -6.4547464e-1/, qg3 / 1.0290328e1/
c++   code for digits > 5 & digits <= 8 is inactive
c      data pg3 /-1.5398760302e-1, -4.2772857997e0,
c     *          -6.5940181141e0 /
c      data qg3 / 3.4149986477e1,   1.7071708816e2 /
c++   code for digits > 8 & digits <= 10 is inactive
c      data pg3 /-1.539896971616e-1, -1.024638314446e1,
c     *          -1.252507402120e2,  -1.029189676144e2 /
c      data qg3 / 7.292229303754e1,   1.186528673248e3,
c     *           3.844430908476e3 /
c++   code for digits > 10 & digits <= 12 is inactive
c      data pg3 /-1.53989733057971e-1, -1.86483223831639e1,
c     *          -5.66882778026550e2,  -4.16714647017489e3,
c     *          -2.14678074364341e3 /
c      data qg3 / 1.27484298075273e2,   4.40258858159927e3,
c     *           4.57583830694684e4,   1.08207883281275e5 /
c++   code for digits = 13 is inactive
c      data pg3 /-1.539897338019247e-1, -2.959078231855258e1,
c     *          -1.661867087183632e3,  -3.112796454920657e4,
c     *          -1.573536286819766e5,  -5.574884113746041e4 /
c      data qg3 / 1.985439813544440e2,   1.196693134508007e4,
c     *           2.625573864512749e5,   1.969055829311119e6,
c     *           3.629081333313312e6 /
c++   code for digits > 13 is active
      data pg3 /-1.53989733819769316e-1 ,
     *          -4.31710157823357568e1  ,
     *          -3.87754141746378493e3  ,
     *          -1.35678867813756347e5  ,
     *          -1.77758950838029676e6  ,
     *          -6.66907061668636416e6  ,
     *          -1.72590224654836845e6 /
      data qg3 / 2.86733194975899483e2  ,
     *           2.69183180396242536e4  ,
     *           1.02878693056687506e6  ,
     *           1.62095600500231646e7  ,
     *           9.38695862531635179e7  ,
     *           1.40622441123580005e8 /
c++   end
c     ------------------------------------------------------------------
      if (bigx .lt. 0.0e0) then
         bigx = 1.0e0 / R1MACH(4)
         largef = rpi / R1MACH(1)
         largeg = (rpi * largef) ** (1.0e0 / 3.0e0)
         largex = 1.0e0/sqrt(sqrt(R1MACH(1)))
      end if
      if (x .ne. xsave) then
         havec = .false.
         havef = .false.
         haveg = .false.
         haves = .false.
      end if
      ax = abs(x)
      if (ax .le. 1.6e0) then
         x4 = ax**4
         if (needc(mode) .and. .not. havec) then
            if (ax .le. 1.2e0) then
c++             code for digits <= 3 is inactive
c               c = x * (pc1(1)*x4+pc1(0)) / (qc1(1)*x4+1.0e0)
c++             code for digits > 3 & digits <= 6 is inactive
c               c = x * ((pc1(2)*x4+pc1(1))*x4+pc1(0))
c     1               / ((qc1(2)*x4+qc1(1))*x4+1.0e0)
c++             code for digits > 6 & digits <= 11 is inactive
c               c = x * (((pc1(3)*x4+pc1(2))*x4+pc1(1))*x4+pc1(0))
c     1               / (((qc1(3)*x4+qc1(2))*x4+qc1(1))*x4+1.0e0)
c++             code for digits > 11 is active
               c = x * ((((pc1(4)*x4+pc1(3))*x4+pc1(2))*x4+pc1(1))*x4+
     1                     pc1(0))
     2           / ((((qc1(4)*x4+qc1(3))*x4+qc1(2))*x4+qc1(1))*x4+1.0e0)
c++             end
            else
c++             code for digits <= 2 is inactive
c               c = x * (pc2(1)*x4+pc2(0)) / (qc2(1)*x4+1.0e0)
c++             code for digits > 2 & digits <= 5 is inactive
c               c = x * ((pc2(2)*x4+pc2(1))*x4+pc2(0))
c     1               / ((qc2(2)*x4+qc2(1))*x4+1.0e0)
c++             code for digits > 5 & digits <= 8 is inactive
c               c = x * (((pc2(3)*x4+pc2(2))*x4+pc2(1))*x4+pc2(0))
c     1               / (((qc2(3)*x4+qc2(2))*x4+qc2(1))*x4+1.0e0)
c++             code for digits > 8 & digits <= 12 is inactive
c               c = x * ((((pc2(4)*x4+pc2(3))*x4+pc2(2))*x4+pc2(1))*x4+
c     1                   pc2(0))
c     2               / ((((qc2(4)*x4+qc2(3))*x4+qc2(2))*x4+qc2(1))*x4+
c     3                   1.0e0)
c++             code for digits > 12 is active
               c = x * (((((pc2(5)*x4+pc2(4))*x4+pc2(3))*x4+pc2(2))*x4+
     1                     pc2(1))*x4+pc2(0))
     2               / (((((qc2(5)*x4+qc2(4))*x4+qc2(3))*x4+qc2(2))*x4+
     3                   qc2(1))*x4+1.0e0)
c++             end
            end if
            havec = .true.
         end if
         if (needs(mode) .and. .not. haves) then
            if (ax .le. 1.2e0) then
c++             code for digits <= 3 is inactive
c               s = x**3*(ps1(1)*x4+ps1(0)) / (qs1(1)*x4+1.0e0)
c++             code for digits > 3 & digits <= 7 is inactive
c               s = x**3*((ps1(2)*x4+ps1(1))*x4+ps1(0))
c     1                / ((qs1(2)*x4+qs1(1))*x4+1.0e0)
c++             code for digits > 7 & digits <= 12 is inactive
c               s = x**3*(((ps1(3)*x4+ps1(2))*x4+ps1(1))*x4+ps1(0))
c     1                / (((qs1(3)*x4+qs1(2))*x4+qs1(1))*x4+1.0e0)
c++             code for digits > 12 is active
               s = x**3*((((ps1(4)*x4+ps1(3))*x4+ps1(2))*x4+ps1(1))*x4+
     1                      ps1(0))
     2           / ((((qs1(4)*x4+qs1(3))*x4+qs1(2))*x4+qs1(1))*x4+1.0e0)
c++             end
            else
c++             code for digits <= 2 is inactive
c               s = x**3*(ps2(1)*x4+ps2(0)) / (qs2(1)*x4+1.0e0)
c++             code for digits > 2 & digits <= 5 is inactive
c               s = x**3*((ps2(2)*x4+ps2(1))*x4+ps2(0))
c     1              /   ((qs2(2)*x4+qs2(1))*x4+1.0e0)
c++             code for digits > 5 & digits <= 9 is inactive
c               s = x**3*(((ps2(3)*x4+ps2(2))*x4+ps2(1))*x4+ps2(0))
c     1              /   (((qs2(3)*x4+qs2(2))*x4+qs2(1))*x4+1.0e0)
c++             code for digits > 9 & digits <= 14 is inactive
c               s = x**3*((((ps2(4)*x4+ps2(3))*x4+ps2(2))*x4+ps2(1))*x4+
c     1                      ps2(0))
c     2              /   ((((qs2(4)*x4+qs2(3))*x4+qs2(2))*x4+qs2(1))*x4+
c     3                      1.0e0)
c++             code for digits > 14 is active
               s = x**3*(((((ps2(5)*x4+ps2(4))*x4+ps2(3))*x4+ps2(2))*x4+
     1                      ps2(1))*x4+ps2(0))
     2              /   (((((qs2(5)*x4+qs2(4))*x4+qs2(3))*x4+qs2(2))*x4+
     3                      qs2(1))*x4+1.0e0)
c++             end
            end if
            haves = .true.
         end if
         if ((needf(mode) .or. needg(mode)) .and. .not. havef) then
            cx = scspxx (ax)
            sx = ssnpxx (ax)
            f = (0.5e0 - s) * cx - (0.5e0 - c) * sx
            g = (0.5e0 - c) * cx + (0.5e0 - s) * sx
            havef = .true.
         end if
      else
         if (ax .le. largex) then
            x4 = (1.0e0 / ax) ** 4
            wantf = needf(mode+4) .and. .not. havef
            if (wantf) then
               if (ax .le. 1.9e0) then
c++                code for digits <= 5 is inactive
c                  f = (pf1(1)*x4+pf1(0)) / ((qf1(1)*x4+1.0e0) * ax)
c++                code for digits > 5 & digits <= 8 is inactive
c                  f = ((pf1(2)*x4+pf1(1))*x4+pf1(0))
c     1             / (((qf1(2)*x4+qf1(1))*x4+1.0e0) * ax)
c++                code for digits > 8 & digits <= 11 is inactive
c                  f = (((pf1(3)*x4+pf1(2))*x4+pf1(1))*x4+pf1(0))
c     1             / ((((qf1(3)*x4+qf1(2))*x4+qf1(1))*x4+1.0e0) * ax)
c++                code for digits > 11 & digits <= 14 is inactive
c                  f = ((((pf1(4)*x4+pf1(3))*x4+pf1(2))*x4+pf1(1))*x4+
c     1                    pf1(0))
c     2             / (((((qf1(4)*x4+qf1(3))*x4+qf1(2))*x4+qf1(1))*x4+
c     3                    1.0e0) * ax)
c++                code for digits > 14 is active
                  f = (((((pf1(5)*x4+pf1(4))*x4+pf1(3))*x4+pf1(2))*x4+
     1                    pf1(1))*x4+pf1(0))
     2             / ((((((qf1(5)*x4+qf1(4))*x4+qf1(3))*x4+qf1(2))*x4+
     3                    qf1(1))*x4+1.0e0) * ax)
c++                end
               else if (ax .le. 2.4) then
c++                code for digits <= 5 is inactive
c                  f = (pf2(1)*x4+pf2(0)) / ((qf2(1))*x4+1.0e0) * ax)
c++                code for digits > 5 & digits <= 8 is inactive
c                  f = ((pf2(2)*x4+pf2(1))*x4+pf2(0))
c     1             / (((qf2(2)*x4+qf2(1))*x4+1.0e0) * ax)
c++                code for digits > 8 & digits <= 11 is inactive
c                  f = (((pf2(3)*x4+pf2(2))*x4+pf2(1))*x4+pf2(0))
c     1             / ((((qf2(3)*x4+qf2(2))*x4+qf2(1))*x4+1.0e0) * ax)
c++                code for digits > 11 & digits <= 14 is inactive
c                  f = ((((pf2(4)*x4+pf2(3))*x4+pf2(2))*x4+pf2(1))*x4+
c     1                    pf2(0))
c     2             / (((((qf2(4)*x4+qf2(3))*x4+qf2(2))*x4+qf2(1))*x4+
c     3                    1.0e0) * ax)
c++                code for digits > 14 is active
                  f = (((((pf2(5)*x4+pf2(4))*x4+pf2(3))*x4+pf2(2))*x4+
     1                    pf2(1))*x4+pf2(0))
     2             / ((((((qf2(5)*x4+qf2(4))*x4+qf2(3))*x4+qf2(2))*x4+
     3                    qf2(1))*x4+1.0e0) * ax)
c++                end
               else
c++                code for digits <= 3 is inactive
c                  f = (rpi + x4*pf3(0)) / ax
c++                code for digits > 3 & digits <= 6 is inactive
c                  f = (rpi + x4*(pf3(1)*x4+pf3(0)) / (qf3(1)*x4+1.0e0))
c                      / ax
c++                code for digits > 6 & digits <= 9 is inactive
c                  f = (rpi +
c     1              x4*((pf3(2)*x4+pf3(1))*x4+pf3(0))
c     2            /    ((qf3(2)*x4+qf3(1))*x4+1.0e0)) / ax
c++                code for digits > 9 & digits <= 11 is inactive
c                  f = (rpi +
c     1              x4*(((pf3(3)*x4+pf3(2))*x4+pf3(1))*x4+pf3(0))
c     2            /    (((qf3(3)*x4+qf3(2))*x4+qf3(1))*x4+1.0e0)) / ax
c++                code for digits > 11 & digits <= 13 is inactive
c                  f = (rpi +
c     1              x4*((((pf3(4)*x4+pf3(3))*x4+pf3(2))*x4+pf3(1))*x4+
c     2                   pf3(0))
c     3            /    ((((qf3(4)*x4+qf3(3))*x4+qf3(2))*x4+qf3(1))*x4+
c     4                   1.0e0)) / ax
c++                code for digits > 13 & digits <= 15 is inactive
c                  f = (rpi +
c     1              x4*(((((pf3(5)*x4+pf3(4))*x4+pf3(3))*x4+pf3(2))*x4+
c     2                   pf3(1))*x4+pf3(0))
c     3            /    (((((qf3(5)*x4+qf3(4))*x4+qf3(3))*x4+qf3(2))*x4+
c     4                   qf3(1))*x4+1.0e0)) / ax
c++                code for digits > 15 is active
                  f = (rpi +
     1              x4*((((((pf3(6)*x4+pf3(5))*x4+pf3(4))*x4+pf3(3))*x4+
     2                   pf3(2))*x4+pf3(1))*x4+pf3(0))
     3            /    ((((((qf3(6)*x4+qf3(5))*x4+qf3(4))*x4+qf3(3))*x4+
     4                   qf3(2))*x4+qf3(1))*x4+1.0e0)) / ax
c++                end
               end if
               havef = .true.
            end if
            wantg = needg(mode+4) .and. .not. haveg
            if (wantg) then
               if (x .le. 1.9e0) then
c++                code for digits <=4 is inactive
c                  g = (pg1(1)*x4+pg1(0)) / ((qg1(1)*x4+1.0e0) * ax**3)
c++                code for digits > 4 & digits <= 7 is inactive
c                  g = ((pg1(2)*x4+pg1(1))*x4+pg1(0))
c     1             / (((qg1(2)*x4+qg1(1))*x4+1.0e0) * ax**3)
c++                code for digits > 7 & digits <= 10 is inactive
c                  g = (((pg1(3)*x4+pg1(2))*x4+pg1(1))*x4+pg1(0))
c     1             / ((((qg1(3)*x4+qg1(2))*x4+qg1(1))*x4+1.0e0) *ax**3)
c++                code for digits > 10 & digits <= 13 is inactive
c                  g = ((((pg1(4)*x4+pg1(3))*x4+pg1(2))*x4+pg1(1))*x4+
c     1                    pg1(0))
c     2             / (((((qg1(4)*x4+qg1(3))*x4+qg1(2))*x4+qg1(1))*x4+
c     3                    1.0e0) * ax**3)
c++                code for digits > 13 is active
                  g = (((((pg1(5)*x4+pg1(4))*x4+pg1(3))*x4+pg1(2))*x4+
     1                    pg1(1))*x4+pg1(0))
     2             / ((((((qg1(5)*x4+qg1(4))*x4+qg1(3))*x4+qg1(2))*x4+
     3                    qg1(1))*x4+1.0e0) * ax**3)
c++                end
               else if (ax .le. 2.4e0) then
c++                code for digits <=4 is inactive
c                  g = (pg2(1)*x4+pg2(0)) / ((qg2(1)*x4+1.0e0) * ax**3)
c++                code for digits > 4 & digits <= 7 is inactive
c                  g = ((pg2(2)*x4+pg2(1))*x4+pg2(0))
c     1             / (((qg2(2)*x4+qg2(1))*x4+1.0e0) * ax**3)
c++                code for digits > 7 & digits <= 10 is inactive
c                  g = (((pg2(3)*x4+pg2(2))*x4+pg2(1))*x4+pg2(0))
c     1             / ((((qg2(3)*x4+qg2(2))*x4+qg2(1))*x4+1.0e0) *ax**3)
c++                code for digits > 10 & digits <= 13 is inactive
c                  g = ((((pg2(4)*x4+pg2(3))*x4+pg2(2))*x4+pg2(1))*x4+
c     1                    pg2(0))
c     2             / (((((qg2(4)*x4+qg2(3))*x4+qg2(2))*x4+qg2(1))*x4+
c     3                    1.0e0) * ax**3)
c++                code for digits > 13 is active
                  g = (((((pg2(5)*x4+pg2(4))*x4+pg2(3))*x4+pg2(2))*x4+
     1                     pg2(1))*x4+pg2(0))
     2             / ((((((qg2(5)*x4+qg2(4))*x4+qg2(3))*x4+qg2(2))*x4+
     3                    qg2(1))*x4+1.0e0) * ax**3)
c++                end
               else
c++                code for digits <=3 is inactive
c                  g = (rpisq + x4*pg3(0)) / ax**3
c++                code for digits > 3 & digits <= 5 is inactive
c                  g = (rpisq + x4*(pg3(1)*x4+pg3(0))
c     1            /    (qg3(1)*x4+1.0e0)) / ax**3
c++                code for digits > 5 & digits <= 8 is inactive
c                  g = (rpisq + x4*((pg3(2)*x4+pg3(1))*x4+pg3(0))
c     1            /    ((qg3(2)*x4+qg3(1))*x4+1.0e0)) / ax**3
c++                code for digits > 8 & digits <= 10 is inactive
c                  g = (rpisq +
c     1              x4*(((pg3(3)*x4+pg3(2))*x4+pg3(1))*x4+pg3(0))
c     2            /    (((qg3(3)*x4+qg3(2))*x4+qg3(1))*x4+1.0e0))
c     3            / ax**3
c++                code for digits > 10 & digits <= 12 is inactive
c                  g = (rpisq +
c     1              x4*((((pg3(4)*x4+pg3(3))*x4+pg3(2))*x4+pg3(1))*x4+
c     2                     pg3(0))
c     3            /    ((((qg3(4)*x4+qg3(3))*x4+qg3(2))*x4+qg3(1))*x4+
c     4                     1.0e0)) / ax**3
c++                code for digits = 13 is inactive
c                  g = (rpisq +
c     1              x4*(((((pg3(5)*x4+pg3(4))*x4+pg3(3))*x4+pg3(2))*x4+
c     2                      pg3(1))*x4+pg3(0))
c     3            /    (((((qg3(5)*x4+qg3(4))*x4+qg3(3))*x4+qg3(2))*x4+
c     4                      qg3(1))*x4+1.0e0)) / ax**3
c++                code for digits > 13 is active
                  g = (rpisq +
     1              x4*((((((pg3(6)*x4+pg3(5))*x4+pg3(4))*x4+pg3(3))*x4+
     2                   pg3(2))*x4+pg3(1))*x4+pg3(0))
     3            /    ((((((qg3(6)*x4+qg3(5))*x4+qg3(4))*x4+qg3(3))*x4+
     4                   qg3(2))*x4+qg3(1))*x4+1.0e0)) / ax**3
c++                end
               end if
               haveg = .true.
            end if
         else
            wantf = needf(mode)
            if (wantf) then
               if (x .le. largef) then
                  f = rpi / ax
               else
                  f = 0.0e0
               end if
            end if
            wantg = needg(mode)
            if (wantg) then
               if (x .le. largeg) then
                  g = rpisq / ax**3
               else
                  g = 0.0e0
               end if
            end if
         end if
         wantc = (needc(mode+4) .or. needs(mode+4)) .and. .not. havec
         if (wantc .or. x.lt.0.0e0) then
            cx = scspxx (ax)
            sx = ssnpxx (ax)
            if (wantc) then
               c = 0.5e0 + f*sx - g*cx
               s = 0.5e0 - f*cx - g*sx
               if (x .lt. 0.0) then
                  c = -c
                  s = -s
               end if
               havec = .true.
            end if
            if (x .lt. 0.0) then
c              We COULD do the following before the preceeding, and then
c              not put in a test in the preceeding for x .lt. 0, but
c              even though the results are mathematically identical, we
c              would have some cancellation above if we did so.
               if (wantg) g = cx + sx - g
               if (wantf) f = cx - sx - f
            end if
          end if
      end if
      xsave = x
      go to (1, 2, 3, 4), MODE
    1    SFREN1 = c
         return
    2    SFREN1 = s
         return
    3    SFREN1 = f
         return
    4    SFREN1 = g
         return
      end
