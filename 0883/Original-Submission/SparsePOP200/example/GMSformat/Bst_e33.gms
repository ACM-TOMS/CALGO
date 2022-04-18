*  NLP written by GAMS Convert at 08/30/02 09:08:49
*  
*  Equation counts
*     Total       E       G       L       N       X       C
*         7       5       0       2       0       0       0
*  
*  Variable counts
*                 x       b       i     s1s     s2s      sc      si
*     Total    cont  binary integer    sos1    sos2   scont    sint
*        10      10       0       0       0       0       0       0
*  FX     0       0       0       0       0       0       0       0
*  
*  Nonzero counts
*     Total   const      NL     DLL
*        30      23       7       0
*
*  Solve m using NLP minimizing objvar;


Variables  x1,x2,x3,x4,x5,x6,x7,x8,x9,objvar;

Positive Variables x1,x2,x3,x4,x5,x6,x7,x8;

Equations  e1,e2,e3,e4,e5,e6,e7;


e1..    x1 + x2 - x3 - x4 =E= 0;

e2..  - x9*(x3 + x4) + 0.03*x1 + 0.01*x2 =E= 0;

e3..    x3 - x5 + x6 =E= 0;

e4..    x4 + x7 - x8 =E= 0;

e5.. x9*x3 - 0.025*x5 + 0.02*x6 =L= 0;

e6.. x9*x4 + 0.02*x7 - 0.015*x8 =L= 0;

e7..  - 6*x1 - 16*x2 + 9*x5 - 10*x6 - 10*x7 + 15*x8 + objvar =E= 0;

* set non default bounds

x1.up = 300; 
x2.up = 300; 
x3.up = 100; 
x4.up = 200; 
x5.up = 100; 
x6.up = 300; 
x7.up = 100; 
x8.up = 200; 
x9.lo = 0.01; 
x9.up = 0.03; 

* set non default levels


* set non default marginals


Model m / all /;

m.limrow=0; m.limcol=0;

$if NOT '%gams.u1%' == '' $include '%gams.u1%'

Solve m using NLP minimizing objvar;
