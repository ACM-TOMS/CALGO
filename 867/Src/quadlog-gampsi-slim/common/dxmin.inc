************************************************************************
*     Smallest representable nonzero normalized double-precision
*     floating-point number (because some architectures carry severe
*     performance penalties for denormalized values, and others do not
*     even implement them, we avoid using them here):
*
*     Maple V5.1:
*         Digits := 50;
*
*         printf("%.50e\n", 2^(-1022));
*             2.22507385850720138309023271733240406421921598046230e-308
*
*     Tests of 41 Fortran compilers on 10 UNIX architectures showed
*     that all successfully converted values with 17 or more digits to
*     the correct value, and none erroneously produced zero instead.
*     The nagf90 and nagf95 compilers on some systems refused to
*     compile this constant with more than 36 digits, so we reduce
*     the precision, being careful to round up in the last digit.
************************************************************************
      DOUBLE PRECISION    xmin
      PARAMETER           (xmin =
     X    2.22507385850720138309023271733240407d-308)
