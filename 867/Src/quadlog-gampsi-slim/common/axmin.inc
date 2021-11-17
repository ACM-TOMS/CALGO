************************************************************************
*     Smallest representable nonzero normalized single-precision
*     floating-point number (because some architectures carry severe
*     performance penalties for denormalized values, and others do not
*     even implement them, we avoid using them here):
*
*     Maple V5.1:
*         Digits := 75;
*
*         # IEEE 754:
*         printf("%.50e\n", 2^(-126));
*             1.17549435082228750796873653722224567781866555677209e-38
*
************************************************************************
      REAL             xmin
*     IEEE 754:
      PARAMETER           (xmin =
     X    1.17549435082228750796873653722224567781866555677209e-38)
