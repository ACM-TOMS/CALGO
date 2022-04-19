*  NLP written by GAMS Convert at 07/30/01 17:04:16
*  
*  Equation counts
*     Total       E       G       L       N       X
*        10       8       0       2       0       0
*  
*  Variable counts
*                 x       b       i     s1s     s2s      sc      si
*     Total    cont  binary integer    sos1    sos2   scont    sint
*        13      13       0       0       0       0       0       0
*  FX     0       0       0       0       0       0       0       0
*  
*  Nonzero counts
*     Total   const      NL     DLL
*        34      27       7       0
*
*  Solve m using NLP minimizing objvar;


Variables  x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,objvar;

Positive Variables x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12;

Equations  e1,e2,e3,e4,e5,e6,e7,e8,e9,e10;


e1..    x1 - 6*x3 - 16*x4 - 10*x5 =E= 0;

e2..    x2 - 9*x6 - 15*x7 =E= 0;

e3..    x6 - x8 - x10 =E= 0;

e4..    x7 - x9 - x11 =E= 0;

e5..    x3 + x4 - x10 - x11 =E= 0;

e6..    x5 - x8 - x9 =E= 0;

e7.. x12*(x10 + x11) - 3*x3 - x4 =E= 0;

e8.. x12*x10 - 2.5*x10 - 0.5*x8 =L= 0;

e9.. x12*x11 - 1.5*x11 + 0.5*x9 =L= 0;

e10..    x1 - x2 - objvar =E= 0;

* set non default bounds

* upper bounds added ---> 
x1.up = 5000;
x2.up = 5000;
x3.up = 300;
x4.up = 300;
x5.up = 300;
* <--- upper bounds added
x6.up = 100; 
x7.up = 200; 
* upper bounds added ---> 
x8.up = 100;
x9.up = 200;
x10.up = 100;
x11.up = 200;
x12.up = 1000;
* <--- upper bounds added

* set non default levels

x8.l = 1; 
x9.l = 1; 
x10.l = 1; 
x11.l = 1; 
x12.l = 1; 

* set non default marginals


Model m / all /;

m.limrow=0; m.limcol=0;

$if NOT '%gams.u1%' == '' $include '%gams.u1%'

Solve m using NLP minimizing objvar;
