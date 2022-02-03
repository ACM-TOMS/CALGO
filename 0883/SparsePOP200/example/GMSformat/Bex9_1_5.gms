*  NLP written by GAMS Convert at 07/19/01 13:40:18
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
*        33      23      10       0
*
*  Solve m using NLP minimizing objvar;


Variables  objvar,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14;

Positive Variables x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14;

Equations  e1,e2,e3,e4,e5,e6,e7,e8,e9,e10,e11,e12,e13;


e1..  - objvar - x2 + 10*x3 - x4 =E= 0;

e2..    x2 + x3 + x5 =E= 1;

e3..    x2 + x4 + x6 =E= 1;

e4..    x3 + x4 + x7 =E= 1;

e5..  - x3 + x8 =E= 0;

e6..  - x4 + x9 =E= 0;

e7.. x10*x5 =E= 0;

e8.. x11*x6 =E= 0;

e9.. x12*x7 =E= 0;

e10.. x13*x8 =E= 0;

e11.. x14*x9 =E= 0;

e12..    x10 + x12 - x13 =E= 1;

e13..    x11 + x12 - x14 =E= 1;

* set non default bounds

x2.up = 5.0; 
x3.up = 5.0; 
x4.up = 5.0; 
x5.up = 5.0; 
x6.up = 5.0; 
x7.up = 5.0; 
x8.up = 5.0; 
x9.up = 5.0; 
x10.up = 5.0; 
x11.up = 5.0; 
x12.up = 5.0; 
x13.up = 5.0; 
x14.up = 5.0; 

* set non default levels


* set non default marginals


Model m / all /;

m.limrow=0; m.limcol=0;

$if NOT '%gams.u1%' == '' $include '%gams.u1%'

Solve m using NLP minimizing objvar;
