*
*     Maple V6:
*         Digits := 75;
*         # DEC, HP, IBM S/390 IEEE 754, Sun:
*         xinf := (1 - 2^(-113)) * 2^(16384):
*         printf("%.50e\n", xinf):
*             1.18973149535723176508575932662800701619646905264169e+4932
*         xbig := fsolve(GAMMA(x) = xinf, x, 1755..1756):
*         printf("%.47f\n", xbig):
*             1755.54834290446291700388921607020336178078553303361
*         # IBM RS/6000 and SGI:
*         xinf := (2 - 2^(-52)) * 2^(1023):
*         printf("%.50e\n", xinf):
*             1.79769313486231570814527423731704356798070567525845e+308
*         xbig := fsolve(GAMMA(x) = xinf, x, 170..175):
*         printf("%.48f\n", xbig);
*             171.62437695630272079085475647646959367985015842820538
*
      REAL*16             xbig
*     DEC, HP, IBM S/390 IEEE 754, Sun:
      PARAMETER           (xbig =
     X    1755.54834290446291700388921607020336178078553303361q+00)
*     IBM RS/6000 and SGI:
*     PARAMETER           (xbig =
*    X    171.62437695630272079085475647646959367985015842820538q+00)
