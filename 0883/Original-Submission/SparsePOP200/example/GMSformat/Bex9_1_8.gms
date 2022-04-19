*  NLP written by GAMS Convert at 07/19/01 13:40:19
*  
*  Equation counts
*     Total       E       G       L       N       X
*        13      12       0       1       0       0
*  
*  Variable counts
*                 x       b       i     s1s     s2s      sc      si
*     Total    cont  binary integer    sos1    sos2   scont    sint
*        15      15       0       0       0       0       0       0
*  FX     0       0       0       0       0       0       0       0
*  
*  Nonzero counts
*     Total   const      NL     DLL
*        33      23      10       0
*
*  Solve m using NLP minimizing objvar;

Variables  objvar,x2,x3,x4,x5,x6,x7,x8,x9,x11,x12,x13,x14;

Positive Variables x2,x3,x4,x5,x6,x7,x8,x9,x11,x12,x13,x14;

Equations  e1,e2,e3,e4,e5,e6,e7,e8,e9,e10,e11,e12;


e1..  - objvar - 2*x2 + x3 + 0.5*x4 =E= 0;

e2..    x2 + x3 =L= 2;

e3..  - 2*x2 + x4 - x5 + x6 =E= -2.5;

e4..    x2 - 3*x3 + x5 + x7 =E= 2;

e5..  - x4 + x8 =E= 0;

e6..  - x5 + x9 =E= 0;

e7.. x11*x6 =E= 0;

e8.. x12*x7 =E= 0;

e9.. x13*x8 =E= 0;

e10.. x14*x9 =E= 0;

*e11.. x15*x10 =E= 0;

e11..    x11 - x13 =E= 4;

e12..    x11 + x12 - x14 =E= -1;

* set non default bounds

x2.up = 2;
x3.up = 2; 
x4.up = 10;
x5.up = 8;
x6.up = 10; 
x7.up = 8; 
x8.up = 10;
x9.up = 8;
x11.up = 10;
x12.up = 10;
x13.up = 10;
x14.up = 10;

* set non default levels


* set non default marginals


Model m / all /;

m.limrow=0; m.limcol=0;

$if NOT '%gams.u1%' == '' $include '%gams.u1%'

Solve m using NLP minimizing objvar;
