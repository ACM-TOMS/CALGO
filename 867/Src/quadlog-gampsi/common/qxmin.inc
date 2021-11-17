************************************************************************
*     Smallest representable nonzero normalized quadruple-precision
*     floating-point number (because some architectures carry severe
*     performance penalties for denormalized values, and others do not
*     even implement them, we avoid using them here):
*
*     Maple V5.1:
*         Digits := 50;
*
*         # DEC, HP, IBM S/390 IEEE 754, Sun:
*         printf("%.50e\n", 2^(-16382));
*             3.36210314311209350626267781732175260259807934484650e-4932
*         # IBM RS/6000 and SGI:
*         printf("%.50e\n", 2^(-1022));
*             2.22507385850720138309023271733240406421921598046230e-308
*
*     Tests of 11 Fortran compilers on 10 UNIX architectures showed
*     that all successfully converted rounded-up values with 33 or
*     more digits to the smallest normalized quadruple-precision
*     floating-point value, and none erroneously produced zero
*     instead.
*
*     Tests of 41 Fortran compilers on 10 UNIX architectures showed that
*     all successfully converted values with 17 or more digits to the
*     largest representable double-precision floating-point value, and
*     none erroneously produced Infinity instead.
************************************************************************
      REAL*16             xmin
*     DEC, HP, IBM S/390 IEEE 754, Sun:
      PARAMETER           (xmin =
     X    3.36210314311209350626267781732175260259807934484650q-4932)
*     IBM RS/6000 and SGI:
*     PARAMETER           (xmin =
*    X    2.22507385850720138309023271733240406421921598046230q-308)
