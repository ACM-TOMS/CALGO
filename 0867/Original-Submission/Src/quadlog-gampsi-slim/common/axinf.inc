************************************************************************
*     IEEE 754 single-precision arithmetic:
*
*            1-bit sign, 8-bit exponent, 1 hidden bit, 23-bit
*            fraction, exponent bias = 127
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
************************************************************************
      REAL             xinf
*
*     IEEE 754:
*
      PARAMETER           (xinf =
     X    3.4028234663852885981170418348451692544e+38)
