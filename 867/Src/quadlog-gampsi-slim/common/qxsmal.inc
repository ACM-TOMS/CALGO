*
*     Maple V5.1:
*         Digits := 50;
*         # DEC, HP, IBM S/390 IEEE 754, Sun:
*         eps := 2^(-112);
*         printf("%60.45e\n", evalf(sqrt(3*eps)/Pi));
*             7.651225422980357370063654563932107352494472746e-18
*         # IBM RS/6000 and SGI:
*         eps := 2^(-105);
*         printf("%60.45e\n", evalf(sqrt(3*eps)/Pi));
*             8.656373409562113739944200777301572379756799586e-17
*
      REAL*16             xsmall
*     DEC, HP, IBM S/390 IEEE 754, Sun:
      PARAMETER           (xsmall =
     X    7.651225422980357370063654563932107352494472746q-18)
*     IBM RS/6000 and SGI:
*     PARAMETER           (xsmall =
*    X    8.656373409562113739944200777301572379756799586q-17)
*
