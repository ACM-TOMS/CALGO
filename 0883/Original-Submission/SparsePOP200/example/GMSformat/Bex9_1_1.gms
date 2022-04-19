*  NLP written by GAMS Convert at 07/19/01 13:40:17
*  
*  Equation counts
*     Total       E       G       L       N       X
*        13      13       0       0       0       0
*  
*  Variable counts
*                 x       b       i     s1s     s2s      sc      si
*     Total    cont  binary integer    sos1    sos2   scont    sint
*        14      14       0       0       0       0       0       0
*  FX     0       0       0       0       0       0       0       0
*  
*  Nonzero counts
*     Total   const      NL     DLL
*        37      27      10       0
*
*  Solve m using NLP minimizing objvar;


Variables  x1,x2,objvar,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14;

Positive Variables x1,x2,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14;

Equations  e1,e2,e3,e4,e5,e6,e7,e8,e9,e10,e11,e12,e13;


e1..  - 3*x1 + 2*x2 - objvar - x4 =E= 0;

e2..    x1 + 4*x2 - 2*x4 + x5 =E= 16;

e3..    3*x1 - 2*x2 + 8*x4 + x6 =E= 48;

e4..    x1 - 3*x2 - 2*x4 + x7 =E= -12;

e5..  - x1 + x8 =E= 0;

e6..    x1 + x9 =E= 4;

e7.. x10*x5 =E= 0;

e8.. x11*x6 =E= 0;

e9.. x12*x7 =E= 0;

e10.. x13*x8 =E= 0;

e11.. x14*x9 =E= 0;

e12..    x10 + 3*x11 + x12 - x13 + x14 =E= 1;

e13..    2*x11 - 3*x12 =E= 0;

* set non default bounds

x1.up = 20; 
x2.up = 20; 
x4.up = 20; 
x5.up = 20; 
x6.up = 20; 
x7.up = 20; 
x8.up = 20; 
x9.up = 20; 
x10.up = 20; 
x11.up = 20; 
x12.up = 20; 
x13.up = 20; 
x14.up = 20; 

* set non default levels


* set non default marginals


Model m / all /;

m.limrow=0; m.limcol=0;

$if NOT '%gams.u1%' == '' $include '%gams.u1%'

Solve m using NLP minimizing objvar;
