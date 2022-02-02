*  NLP written by GAMS Convert at 07/19/01 13:40:21
*  
*  Equation counts
*     Total       E       G       L       N       X
*         8       8       0       0       0       0
*  
*  Variable counts
*                 x       b       i     s1s     s2s      sc      si
*     Total    cont  binary integer    sos1    sos2   scont    sint
*         9       9       0       0       0       0       0       0
*  FX     0       0       0       0       0       0       0       0
*  
*  Nonzero counts
*     Total   const      NL     DLL
*        19      13       6       0
*
*  Solve m using NLP minimizing objvar;


Variables  objvar,x2,x3,x4,x5,x6,x7,x8,x9;

Positive Variables x3,x4,x5,x6,x7,x8,x9;

Equations  e1,e2,e3,e4,e5,e6,e7,e8;


* e1.. (0.5*x4 - 1)*(x4 - 2) + (0.5*x5 - 1)*(x5 - 2) - objvar =E= 0;
e1.. 0.5*x4^2 - 2*x4 + 0.5*x5^2 -2*x5 - objvar =E= -4;

e2..  - x3 + x4 + x5 =E= 0;

e3..  - x4 + x6 =E= 0;

e4..  - x5 + x7 =E= 0;

e5.. x6*x8 =E= 0;

e6.. x7*x9 =E= 0;

e7..    x2 + x4 - x8 =E= 0;

e8..    x2 - x9 =E= -1;

* set non default bounds
x2.lo = -1; 
x2.up = 199;
x3.up = 400; 
x4.up = 200; 
x5.up = 200; 
x6.up = 200; 
x7.up = 200; 
x8.up = 200; 
x9.up = 200; 

* set non default levels


* set non default marginals


Model m / all /;

m.limrow=0; m.limcol=0;

$if NOT '%gams.u1%' == '' $include '%gams.u1%'

Solve m using NLP minimizing objvar;
