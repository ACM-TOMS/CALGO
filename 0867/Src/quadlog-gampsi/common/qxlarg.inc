*
*     Maple V5.1:
*         Digits := 50;
*         # DEC, HP, IBM S/390 IEEE 754, Sun:
*         printf("%60.45e\n", fsolve(x*ln(x) = 2^113, x));
*             1.402920063718526600227830569191951958902892674e+32
*         # IBM RS/6000 and SGI:
*         printf("%60.45e\n", fsolve(x*ln(x) = 2^106, x));
*             1.171782682388372063216052705125821796211923258e+30
*
      REAL*16             xlarge
*     DEC, HP, IBM S/390 IEEE 754, Sun:
      PARAMETER           (xlarge = 1.40q+32)
*     IBM RS/6000 and SGI:
*     PARAMETER           (xlarge = 1.17q+30)
