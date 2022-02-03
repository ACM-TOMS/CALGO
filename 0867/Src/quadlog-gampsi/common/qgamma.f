      REAL*16 FUNCTION  qgamma(x)
************************************************************************
*     (Quadruple-precision Gamma(x))
*     Compute and return the value of the Gamma(x) function for
*     quadruple-precision x.
*
*     This code correctly handles the case where x is NaN, for which
*     Gamma(NaN) is a NaN. and the case where x is sufficiently large
*     and positive, or takes one of the integer values 0, -1, -2, -3,
*     ..., for which Gamma(x) = +Infinity, a special value in IEEE 754
*     arithmetic.
*
*     This code is derived from code given by
*
*         W. J. Cody
*         Algorithm 715: SPECFUN --- A Portable FORTRAN
*         Package of Special Function Routines and Test Drivers
*         ACM Trans. Math. Software 19(1) 22--32, March 1993.
*
*     but has been augmented for support of IEEE 754 arithmetic, and
*     extended with a new rational Pade approximation for x*Gamma(x)
*     in [1,2], and new algorithms for argument reduction.
*     (04-Aug-2000)
************************************************************************
*     FOREWORD: This function is derived from code available on the
*     World-Wide Web at ftp://ftp.netlib.org/toms/715, corresponding to
*     the article
*
*         W. J. Cody
*         Algorithm 715: SPECFUN --- A Portable FORTRAN
*         Package of Special Function Routines and Test Drivers
*         ACM Trans. Math. Software 19(1) 22--32 March 1993.
*
*     The code has been prettyprinted to match the conventions of the
*     quadrature package in which it is distributed, and the
*     single-precision code embedded in comments has been discarded.
*
*     Additional code has been added to handle the cases of x = NaN, for
*     which Gamma(NaN) = NaN, and large x, and x = 0, -1, -2, -3, ...,
*     for which Gamma(x) = +Infinity, a special value in IEEE 754
*     arithmetic.
*
*     Where relevant, constants are given to about 50 decimal digits,
*     and the Pad{\'e} polynomial approximation for x in 1..2 has been
*     chosen to give agreement to about 35 decimal digits, which is
*     sufficient for all current quadruple-precision arithmetic
*     implementations (Compaq/DEC VAX, IBM S/360, IEEE 754, and
*     approximations to IEEE 754).
*
*     Numerical tests of the relative errors on several platforms show
*     it to be less than one ulp for x in 0+ .. 1.0q-30, and up to four
*     ulps for x in 1.0q-30 .. xinf, on systems (Compaq/DEC, IBM, Sun)
*     with careful implementations of quadruple-precision arithmetic.
*     On HP systems, worst case relative errors of 1 to 2 decimal digits
*     are seen, and on SGI systems, relative errors of up to 4 decimal
*     digits.
*
*     The Gamma() function has asymptotes at x = 0, -1, -2, -3, ...,
*     with the value approaching +Infinity from one side, and
*     -Infinity from the other.  For such arguments, this function
*     returns qinf(), a run-time trappable +Infinity on IEEE 754 systems.
*
*     Tests with trap handlers enabled for underflow, overflow, invalid
*     operand (NaN), and zero-divide for arguments away from the
*     underflow limit reported no floating-point exceptions.
*
*     To facilitate modifications for other floating-point systems,
*     short Maple programs given in comments can be used to generate
*     needed constants to arbitrary precision.
*
*     The constants eps, xbig, xminin, and xinf below have been given
*     values suitable for IEEE 754 quadruple-precision arithmetic.  This
*     is applicable on all current Intel, Motorola, and most RISC CPUs,
*     and even on IBM S/390 mainframes since 1999.  Older Compaq/DEC
*     VAX, Convex, Cray, and IBM S/3x0 systems will require changes in
*     those values.  It has not yet been possible to carry out numerical
*     tests on other than IEEE 754 arithmetic systems.
*
*     IBM RS/6000 AIX quadruple-precision has a different format,
*     consisting of two consecutive double-precision values, with NO
*     constraint relation between the two exponents.  This reduces the
*     precision from p = 113 bits to p = 106 bits, and the power-of-ten
*     exponent range from -4932 .. +4932 to -308 .. +308.  Comments
*     below give values suitable for this unusual quadruple-precision
*     representation.
*
*     Silicon Graphics (SGI) IRIX quadruple-precision also uses two
*     consecutive double-precision values, but the exponents are
*     somewhat constrained, eliminating the bogus machine-epsilon
*     problem of the IBM RS/6000 format.  The exponent range and
*     precision are the same as on the latter.
*     (04-Aug-2000)
************************************************************************
*
*     NB: The comments in the remainder of this comment header block are
*     from the original ACM Algorithm 715, and do not reflect the higher
*     precision actually implemented in this quadruple-precision version
*     of the algorithms.
*
*     This routine calculates the gamma function for a real argument x.
*
*     Computation is based on an algorithm outlined in W. J. Cody, ``An
*     overview of software development for special functions'', Lecture
*     Notes in Mathematics, 506, Numerical Analysis Dundee, 1975,
*     G. A. Watson (Ed.), Springer Verlag, Berlin, 1976.
*
*     The program uses rational functions that approximate the gamma
*     function to at least 20 significant decimal digits.  Coefficients
*     for the approximation over the interval (1,2) are unpublished.
*     Those for the approximation for x .GE. 12 are from Hart, et. al.,
*     Computer Approximations, Wiley and Sons, New York, 1968.  Lower
*     order approximations can be substituted for these on machines with
*     less precise arithmetic.
*
************************************************************************
*
*     Explanation of machine-dependent constants:
*
*     eps    - the smallest positive floating-point number such that
*              1.0 + eps .gt. 1.0
*     xbig   - the largest argument for which Gamma(x) is
*              representable in the machine, i.e., the solution to the
*              equation
*                      Gamma(xbig) = xinf.
*     xminin - the smallest positive floating-point number such that
*              1/xminin is machine representable.
*     xinf   - the largest machine representable floating-point number.
*
*     Approximate values for some important machines are:
*
*             IBM/195    CDC/7600  UNIVAC/1108   VAX 11/780 (UNIX)
*              (d.p.)  (s.p.,rndg)    (d.p.)     (s.p.)     (d.p.)
*
*     eps     2.221d-16  3.553e-15   1.735d-18   5.961e-08  1.388d-17
*     xbig    57.574     177.802     171.489     34.844     34.844
*     xminin  1.382d-76  3.132e-294  1.113d-308  5.883e-39  5.883d-39
*     xinf    7.23d+75   1.26e+322   8.98d+307   1.70e+38   1.70d+38
*
************************************************************************
*
*     Error returns:
*
*      The program returns the value xinf for singularities or
*      when overflow would occur.  The computation is believed
*      to be free of underflow and overflow.
*
*     Other subprograms required (single precision version):
*
*         alog,exp,float,ifix,sin
*
*     Other subprograms required (double precision version):
*
*         dble,dexp,dlog,dsin,float,ifix,sngl
*
*      Author: W. J. Cody
*              Applied Mathematics Division
*              Argonne National Laboratory
*              Argonne, IL 60439
*
*      Latest modification: May 18, 1982
*
************************************************************************
*
*     Intrinsic functions
*
      INTRINSIC           iqint,       mod,         qcos,        qexp
      INTRINSIC           qint,        qlog,        qsin
*
*     Built-in functions
*
      INTEGER             iqint,       mod
*
      REAL*16             qcos,        qexp,        qint,        qlog
      REAL*16             qsin
*
*     External functions
*
      EXTERNAL            isqnan,      qfloat,      qinf,        qnan
*
      LOGICAL             isqnan
*
      REAL*16             qfloat,      qinf,        qnan
*
*     Parameter variables
*
      INTEGER             npq
      PARAMETER           (npq = 16)
*
      REAL*16             half
      PARAMETER           (half = 0.5q+00)
*
      REAL*16             one
      PARAMETER           (one = 1.0q+00)
*
      REAL*16             two
      PARAMETER           (two = 2.0q+00)
*
      REAL*16             four
      PARAMETER           (four = 4.0q+00)
*
      REAL*16             zero
      PARAMETER           (zero = 0.0q+00)
*
      INCLUDE 'qeps.inc'
*
      INCLUDE 'ql2pi2.inc'
*
*     REAL*16             pi
*     PARAMETER           (pi =
*    X    3.14159265358979323846264338327950288419716939937511q+00)
*
*     piov4 = pi / 4
*
      REAL*16             piov4
      PARAMETER           (piov4 =
     X    7.853981633974483096156608458198757210492923498437764552q-01)
*
      INCLUDE 'qxinf.inc'
*
      INCLUDE 'qxbig.inc'
*
      INCLUDE 'qxmin.inc'
*
*     Maple V6:
*         Digits := 0;
*         # DEC, HP, IBM S/390 IEEE 754, Sun:
*         printf("%.49e\n", 2^(-16382)):
*             3.3621031431120935062626778173217526025980793448465e-4932
*             (smallest normalized floating-point number)
*         printf("%.49e\n", 2^(16382)):
*             2.9743287383930794127143983165700178269086117177413e+4931
*             (representable)
*         printf("%.49e\n", 1/((1 - 2^(-113)) * 2^(16384))):
*             8.4052578577802337656566945433043823158920882918693e-4933
*             (NOT representable as normalized)
*         # IBM RS/6000 and SGI:
*         printf("%.49e\n", 2^(-1022)):
*             2.2250738585072013830902327173324040642192159804623e-308
*         printf("%.49e\n", 2^(1022)):
*             4.4942328371557897693232629769725618340449424473558e+307
*             (representable)
*         printf("%.49e\n", 1/((1 - 2^(-53)) * 2^(1024))):
*             5.5626846462680040753076390948892589466409921927067e-309
*             (NOT representable as normalized)
*
*     Thus, the reciprocal is representable, since it is less than xinf.
*     There are smaller denormalized values whose reciprocals are
*     representable, the smallest being
*
*         # DEC, HP, IBM S/390 IEEE 754, Sun:
*         printf("%.49e\n", 2^(-16382 - 112)):
*             6.4751751194380251109244389582276465524995693380347e-4966
*         # IBM RS/6000 and SGI:
*         printf("%.49e\n", 2^(-1022 - 53)):
*             2.4703282292062327208828439643411068618252990130716e-324
*
*     However, some architectures carry severe performance penalties for
*     denormalized values, so we avoid using them.
*
      REAL*16             xminin
      PARAMETER           (xminin = xmin)
*
*     The asymptotic series for ln(Gamma(x)) is used for x in
*     (cutoff,INF).  The number of terms computed is nc (see above),
*     which determines the size of the array c(*) below, and the
*     number of precomputed initial values stored in it.
*
*     Subject to that limit, the cutoff can be adjusted upwards to
*     reduce the accuracy loss from the use of qlog() and qexp()
*     functions.  In practice, this is what limits accuracy of qgamma(x)
*     for large x.  About two decimal figures are lost by switching to
*     the asymptotic series.
*
*     For IEEE 754 quadruple precision, p = 113, and ulp = 2**(p-1) =
*     1.93e-34.
*
*     For IBM RS/6000 and SGI quadruple precision, p = 105, and ulp =
*     2**(p-1) = 2.47e-32.
*
*     The numbers of terms tabulated below were determined by a Maple
*     program that exited the summation loop when the relative size of
*     the last-added term dropped below the ulp value. The starred
*     entries correspond to integer cutoff values just below the limit
*     at which Gamma(x) overflows.
*
*      =======     ======================
*       cutoff         number of terms
*                  IEEE 754   IBM RS/6000
*                               and SGI
*      =======     ======================
*           10          201         201
*           11          201          29
*           12           29          23
*           13           24          20
*           14           22          19
*           15           20          17
*           16           19          17
*           17           18          16
*           18           17          15
*           19           16          15
*           20           16          14
*           21           15          14
*           22           15          13
*           23           14          13
*           25           14          12
*           26           13          12
*           30           13          11
*           31           12          11
*           37           11          10
*           47           10          10
*           49           10           9
*           64            9           9
*           71            9           8
*           95            8           8
*          115            8           7
*          162            7           7
*          171            7       *** 7
*          226            7           6
*          336            6           6
*          603            6           5
*          971            5           5
*         1755        *** 5           5
*         5000            5           4
*         6000            4           4
*        10000            4           4
*       100000            3           3
*      1000000            3           3
*      =======     ======================
*
*     Although the asymptotic series could be truncated earlier for
*     large x, that would preclude the use of the fast Horner nested
*     multiplication, so the code below always uses a fixed number of
*     terms.
*
*     In order to avoid the loss of about two decimal digits for large
*     arguments, this version of the code sets the cutoff value larger
*     than xbig, so that the asymptotic series will NEVER be used.
*
*     It must be emphasized that the error in the asymptotic series is
*     NOT from inaccuracies in qexp() and qlog(), but rather from the
*     unavoidable relative error magnification in exp(x), which scales
*     like the argument.  From elementary calculus, relerr(f(x)) = x *
*     (fprime(x)/f(x)) * relerr(x), so relerr(exp(x)) = x * relerr(x).
*     For x .EQ. 100, two digits are unavoidably lost, and this effect
*     is very evident in plots of the log of the absolute value of the
*     relative error against x: at the cutoff, there is a dramatic drop
*     in accuracy.
*
*     The penalty for avoiding the asymptotic series is up to
*     x floating-point multiplies from the recursion Gamma(x+1) =
*     x*Gamma(x), making the code somewhat slower for large arguments,
*     but the benefit is significantly enhanced accuracy.
*
      REAL*16             cutoff
*     PARAMETER           (cutoff = 12.0q+00)
*     PARAMETER           (cutoff = 18.0q+00)
*     PARAMETER           (cutoff = 25.0q+00)
*     PARAMETER           (cutoff = 100.0q+00)
*     PARAMETER           (cutoff = 1000.0q+00)
      PARAMETER           (cutoff = xbig + one)
*
*     ncuse is the number of terms actually used in the summation of the
*     asymptotic series; it must be set to match cutoff values in the
*     table above, and must also be in the range 1..nc.
*
      INTEGER             ncuse
*     PARAMETER           (ncuse = 29)
*     PARAMETER           (ncuse = 17)
*     PARAMETER           (ncuse = 14)
*     PARAMETER           (ncuse = 8)
*     PARAMETER           (ncuse = 5)
      PARAMETER           (ncuse = 5)
*
*     Argument variables
*
      REAL*16             x
*
*     Local variables
*
      INTEGER             i,           n,           umod8
*
      LOGICAL             parity
*
      REAL*16             fact,        p(npq),      q(npq),      result
      REAL*16             sinval,      sum,         u,           v
      REAL*16             w,           xden,        xnum,        y
      REAL*16             y1,          yfrac,       yint,        ysqinv
      REAL*16             z
*
*-----------------------------------------------------------------------
*     Numerator and denominator coefficients for rational Pad{\'e}
*     approximation of x*Gamma(x) over [1,2], computed in Maple from
*     xgampade.map, the critical kernel of which is
*
*         Digits := 250;
*         the_approx := pade('x'*GAMMA('x'), 'x' = 1.5, [15,15]):
*
*     Maximum relative error in Gamma(x) = 1.84e-36
*     at x = 1.0000 for x on [1,2] with Pade degree [15,15]
*-----------------------------------------------------------------------
      DATA p /
     X   3.6135523544306598716552588412695060493014867824318004515q-01,
     X   2.8818029022995208875842831666077908230519147436798285210q-01,
     X   1.3503819200503020719414581752972474066235674547490945473q-01,
     X   4.4375365153128757820728276535292785241838397698368086060q-02,
     X   1.1421462791142745982654316776504349501851191038534645809q-02,
     X   2.3904572277690682235111327212058189148210417297532561119q-03,
     X   4.2166356067825107442394987772842408012194655441541586203q-04,
     X   6.3357745620341267995600028088153345545272688239610034514q-05,
     X   8.2411401587982997113706993384408325145665835983531226944q-06,
     X   9.2503162271282023453917129714955744235163393710655384060q-07,
     X   9.0235003446813139727911126553585258987256515429453721083q-08,
     X   7.5006893823952877333966241443823690708044271063516432811q-09,
     X   5.3067143220820585748515032453178395409317546081165585093q-10,
     X   3.0174329665515119279666429573594868907032721735080469620q-11,
     X   1.3270356888012001336209797431388382462571725343016693689q-12,
     X   3.7197137589151210347452188276357493996725738165846843864q-14 /
      DATA q /
     X   3.6135523544306598716322080629544401703731661826821707767q-01,
     X   4.9676019272187137588968132132574393227835922770496574306q-01,
     X   6.4375394885483302025039816790289219754351284786676800168q-02,
     X  -8.1867480139729561423544667922545858518464347180460796538q-02,
     X  -3.4577026887258298576651469940520862753335840263462379455q-03,
     X   6.9510887099734225848145981013457732742924066192890034558q-03,
     X  -7.0192853565588217241275948792454052308195337506535207766q-04,
     X  -2.1318515145476545192318792191550191423634710172969933255q-04,
     X   5.7405623408060753599202564969260608435685761385398808853q-05,
     X  -3.1429140956847597747670694873706069906380556216380536609q-06,
     X  -7.4345368999022378486135311168119148154467799579100931497q-07,
     X   1.7012094423833662452290119949173386670892329150502228806q-07,
     X  -1.6929974237066533451436087535720618742418690925567075938q-08,
     X   9.6510934451089183645265124480142619453059072223597097147q-10,
     X  -3.0894496234882073515506394582340492027251781345226329568q-11,
     X   4.3649844254680481243346166967009621318574333862749601292q-13 /
*-----------------------------------------------------------------------
      INCLUDE 'qcgam.inc'
*-----------------------------------------------------------------------
      parity = .FALSE.
      fact = one
      n = 0
      y = x
      IF (isqnan(x)) THEN
*-----------------------------------------------------------------------
*          Argument is a NaN, so Gamma(NaN) is too.
*-----------------------------------------------------------------------
           qgamma = qnan()
           RETURN
      ELSE IF (y .LE. zero) THEN
*-----------------------------------------------------------------------
*          Argument is zero or negative.  Gamma(x) is finite, except at
*          the poles where x = 0, -1, -2, -3, ...
*-----------------------------------------------------------------------
           y = -x
           yint = qint(y)
           yfrac = y - yint
           IF (yfrac .NE. zero) THEN
                parity = (yint .NE. qint(yint*half)*two)
*-----------------------------------------------------------------------
*               We want to set fact = -pi/qsin(pi*yfrac).  However, pi
*               has two leading zero bits in IBM S/360 arithmetic, so
*               compute this more accurately as fact = -4*(pi/4) /
*               sin((pi/4)*(4*yfrac)), and reduce this further by
*               4*yfrac = w = u (integer) + v (fraction), where w, u,
*               and v are all nonnegative and small, producing eight
*               cases involving arguments with v and 1-v (which does not
*               suffer from subtraction loss).
*-----------------------------------------------------------------------
                w = four * yfrac
                u = qint(w)
                v = w - u
                umod8 = mod(iqint(u),8)
                IF (umod8 .eq. 0) THEN
                     sinval = qsin(piov4 * v)
                ELSE IF (umod8 .eq. 1) THEN
                     sinval = qcos(piov4 * (one -v))
                ELSE IF (umod8 .eq. 2) THEN
                     sinval = qcos(piov4 * v)
                ELSE IF (umod8 .eq. 3) THEN
                     sinval = qsin(piov4 * (one - v))
                ELSE IF (umod8 .eq. 4) THEN
                     sinval = -qsin(piov4 * v)
                ELSE IF (umod8 .eq. 5) THEN
                     sinval = -qcos(piov4 * (one - v))
                ELSE IF (umod8 .eq. 6) THEN
                     sinval = -qcos(piov4 * v)
                ELSE
                     sinval = -qsin(piov4 * (one - v))
                END IF
                fact = -four * (piov4 / sinval)
*OLD:           fact = - pi/qsin(pi*yfrac)
                y = y + one
           ELSE
*----------------------------------------------------------------------
*               Argument is negative integer, so Gamma(x) = +/-Infinity
*----------------------------------------------------------------------
                result = qinf()
                GO TO 400
           END IF
      END IF
*-----------------------------------------------------------------------
*     Argument is now positive, with fact and parity recording history
*     of negative argument.
*-----------------------------------------------------------------------
      IF (y .LT. eps) THEN
*-----------------------------------------------------------------------
*          Argument .LT. eps
*-----------------------------------------------------------------------
           IF (y .GE. xminin) THEN
                result = one/y
           ELSE
                result = qinf()
                GO TO 400
           END IF
      ELSE IF (y .LT. cutoff) THEN
           y1 = y
           IF (y .LT. one) THEN
*-----------------------------------------------------------------------
*               0.0 .LT. argument .LT. 1.0
*-----------------------------------------------------------------------
                z = y + one
                y = y + one
           ELSE
*-----------------------------------------------------------------------
*               1.0 .LT. argument .LT. cutoff, reduce argument if
*               necessary
*-----------------------------------------------------------------------
                n = iqint(y) - 1
                y = y - qfloat(n)
                z = y
           END IF
*-----------------------------------------------------------------------
*          Evaluate approximation for 1.0 .LT. argument .LT. 2.0
*          Compute the Pade approximation
*              Gamma(x) = (sum(k=1,npq)p(k)x^(k-1)/
*                          sum(k=1,npq)q(k)x^(k-1)) / x
*-----------------------------------------------------------------------
           xnum = p(npq)
           xden = q(npq)
           DO 100 i = npq - 1, 1, -1
                xnum = xnum*z + p(i)
                xden = xden*z + q(i)
  100      CONTINUE
           result = xnum/(xden * z)
           IF (y1 .LT. y) THEN
*-----------------------------------------------------------------------
*               Adjust result for case 0.0 .LT. argument .LT. 1.0
*-----------------------------------------------------------------------
                result = result/y1
           ELSE IF (y1 .GT. y) THEN
*-----------------------------------------------------------------------
*               Adjust result for case 2.0 .LT. argument .LT. cutoff
*-----------------------------------------------------------------------
                DO 200 i = 1, n
                     result = result*y
                     y = y + one
  200           CONTINUE
           END IF
      ELSE
*-----------------------------------------------------------------------
*          Evaluate for argument .GE. cutoff
*-----------------------------------------------------------------------
           IF (y .LE. xbig) THEN
*-----------------------------------------------------------------------
*               Use asymptotic series for ln(Gamma(x)).  The four
*               parts are summed in order from smallest to largest,
*               positive parts first.
*-----------------------------------------------------------------------
                ysqinv = one/(y*y)
                sum = c(ncuse)
                DO 300 i = (ncuse - 1), 1, -1
                     sum = sum*ysqinv + c(i)
  300           CONTINUE
                sum = ((y - half)*qlog(y) +(ln2pi2 + sum/y)) - y
                result = qexp(sum)
           ELSE
                result = qinf()
                GO TO 400
           END IF
      END IF
*-----------------------------------------------------------------------
*     Final adjustments and return
*-----------------------------------------------------------------------
  400 IF (parity) result = - result
      IF (fact .NE. one) result = fact/result
      qgamma = result
      END
