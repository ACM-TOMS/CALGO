************************************************************************
*     Machine epsilon: the smallest number that can be added to 1.0
*     and still differ from it.
*
*     IEEE 754 double-precision normalized numbers have the value
*         1.f * 2^p
*     where f is a 52-bit binary fraction, and p ranges from -1022 to
*     1023.  Denormalized values have the form
*         0.f * 2^p
*     where f may have leading zero bits.
*
*     Maple V5.1:
*         Digits := 25;
*         printf("%.25e\n", 2^(-52))
*             2.2204460492503130808472630e-16
************************************************************************
      DOUBLE PRECISION    eps
      PARAMETER           (eps = 2.2204460492503130808472630d-16)
