************************************************************************
*     There are currently two implementations of IEEE 754
*     quadruple-precision arithmetic:
*
*     (1) DEC, HP, and Sun (all in software), and IBM S/390 G5 (in
*         hardware):
*            1-bit sign, 15-bit exponent, 1 hidden bit, 112-bit
*            fraction, exponent bias = 16383
*
*     (2) IBM RS/6000 and SGI (in software): two double-precision words
*         of
*            1-bit sign, 11-bit exponent, 1 hidden bit, 52-bit fraction,
*            1-bit sign, 11-bit exponent, 1 hidden bit, 52-bit fraction,
*            exponent bias = 1023
*         This corresponds roughly to
*            1-bit sign, 11-bit exponent, 1 hidden bit, 105-bit
*            fraction, exponent bias = 1023
*         although the exponents are independent, instead of differing
*         by a constant value.
*
*     For a p-bit fraction (EXCLUDING the hidden bit), and q-bit
*     exponent, the largest biased exponent is 2^q - 1, but that largest
*     value is reserved for Infinity and NaN, so the largest biased
*     exponent for normal numbers is 2^q - 2.  Thus, the largest
*     representable normal floating-point number is
*
*         xinf = (2 - 2^(-p))*2^(2^q - 2 - bias)
*
*     Maple V5.1:
*        Digits := 50;
*        xinf := proc(p,q,bias)
*                     evalf((2 - 2^(-p))*2^(2^q - 2  - bias))
*                 end;
*        printf("%60.45e\n", xinf(23,8,127)); # IEEE 754 single
*            3.402823466385288598117041834845169254400000000e+38
*        printf("%60.45e\n", xinf(52,11,1023)); # IEEE 754 double
*            1.797693134862315708145274237317043567980705675e+308
*        printf("%60.45e\n", xinf(112,15,16383)); # IEEE 754 quadruple
*            1.189731495357231765085759326628007016196469053e+4932
*        printf("%60.45e\n", xinf(105,11,1023)); # IBM RS/6000 quadruple
*            1.797693134862315907729305190789002575339325774e+308
*
*     It is not possible on paired-double architectures (IBM RS/6000 and
*     SGI) to produce reliably via addition, or input, a value like that
*     last shown.  IBM incorrectly inputs it as 7FEFFFFFFFFFFFFF
*     7C9FFFFFFFFFFFFE (off by 1 in the last bit), and SGI inputs it as
*     Infinity.  On the IBM RS/6000, careful manipulation of sums can
*     produce the even larger value
*
*         3.59538626972463101712248656769011712E+0308 =
*         7FEFFFFF FFFFFFFE 7FEFFFFF FFFFFFFE
*
*     but this still has errors in two bit positions.  We therefore fall
*     back to the largest double-precision number instead for these
*     aberrant systems.
*
*     Maple V5.1:
*         Digits := 75:
*         # DEC, HP, IBM S/390 IEEE 754, Sun:
*         xinf := (1 - 2^(-113)) * 2^(16384):
*         printf("%.50e\n", xinf);
*             1.18973149535723176508575932662800701619646905264169e+4932
*         # IBM RS/6000 and SGI:
*         xinf := (2 - 2^(-52)) * 2^(1023):
*         printf("%.50e\n", xinf):
*             1.79769313486231570814527423731704356798070567525845e+308
*
*     Tests of 14 Fortran compilers on 10 UNIX architectures showed that
*     all successfully converted values with 33 or more digits to the
*     largest representable quadruple-precision floating-point value,
*     and none erroneously produced Infinity instead.
*
*     Tests of 41 Fortran compilers on 10 UNIX architectures showed that
*     all successfully converted values with 17 or more digits to the
*     largest representable double-precision floating-point value, and
*     none erroneously produced Infinity instead.
************************************************************************
      REAL*16             xinf
*
*     DEC, HP, IBM S/390 IEEE 754, Sun:
*
      PARAMETER           (xinf =
     X    1.18973149535723176508575932662800701619646905264169q+4932)
*
*     IBM RS/6000 and SGI:
*
*     PARAMETER           (xinf =
*    X    1.79769313486231570814527423731704356798070567525845q+308)
