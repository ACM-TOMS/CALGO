************************************************************************
*     Machine epsilon: the smallest number that can be added to 1.0
*     and still differ from it.
*
*     IEEE 754 quadruple-precision normalized numbers have the value
*         1.f * 2^p
*     where f is a 112-bit binary fraction, and p ranges from -16382 to
*     16383.  Denormalized values have the form
*         0.f * 2^p
*     where f may have leading zero bits.  The pair-of-doubles format
*     used by IBM RS/6000 and SGI has a 105-bit binary fraction, and
*     p ranges from -1022 to 1023.
*
*     Maple V5.1:
*         Digits := 50;
*         # DEC, HP, IBM S/390 IEEE 754, Sun:
*         printf("%.35e\n", 2^(-112)):
*             1.92592994438723585305597794258492732e-34
*         # IBM RS/6000 and SGI:
*         printf("%.35e\n", 2^(-105)):
*             2.46519032881566189191165176650870697e-32
************************************************************************
      REAL*16             eps
*     DEC, HP, IBM S/390 IEEE 754, Sun:
      PARAMETER           (eps =
     X    1.92592994438723585305597794258492732q-34)
*     IBM RS/6000 and SGI:
*     PARAMETER           (eps =
*    X    2.46519032881566189191165176650870697q-32)
