************************************************************************
*     IEEE 754 double-precision arithmetic:
*
*         1-bit sign, 11-bit exponent, 1 hidden bit, 52-bit fraction,
*         exponent bias = 1023
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
*
*     Tests of 41 Fortran compilers on 10 UNIX architectures showed that
*     all successfully converted values with 17 or more digits to the
*     largest representable double-precision floating-point value, and
*     none erroneously produced Infinity instead.
************************************************************************
      DOUBLE PRECISION    xinf
      PARAMETER           (xinf =
     X    1.797693134862315708145274237317043567980705675d+308)
