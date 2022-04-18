************************************************************************
*     Machine epsilon: the smallest number that can be added to 1.0
*     and still differ from it.
*
*     IEEE 754 single-precision normalized numbers have the value
*         1.f * 2^p
*     where f is a 23-bit binary fraction, and p ranges from -126 to
*     127.  Denormalized values have the form
*         0.f * 2^p
*     where f may have leading zero bits.
*
*     Maple V5.1:
*         Digits := 50;
*         # IEEE 754:
*         printf("%.35e\n", 2^(-23)):
*             1.19209289550781250000000000000000000e-07
************************************************************************
      REAL                eps
*     IEEE 754:
      PARAMETER           (eps = 1.1920928955078125e-07)
