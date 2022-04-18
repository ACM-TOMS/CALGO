*  NLP written by GAMS Convert at 07/19/01 13:39:37
*  
*  Equation counts
*     Total       E       G       L       N       X
*         7       5       0       2       0       0
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

Positive Variables x1,x2,x3,x4,x5,x6,x7,x8,x9;

Equations  e1,e2,e3,e4,e5,e6,e7;


e1..  - 9*x1 - 15*x2 + 6*x3 + 16*x4 + 10*x5 + 10*x6 - objvar =E= 0;
*e1.. x7 + objvar =E= 0; 

e2..  - x3 - x4 + x8 + x9 =E= 0;

e3..    x1 - x5 - x8 =E= 0;

e4..    x2 - x6 - x9 =E= 0;

e5.. x7*x8 - 2.5*x1 + 2*x5 =L= 0;

e6.. x7*x9 - 1.5*x2 + 2*x6 =L= 0;

e7.. x7*x8 + x7*x9 - 3*x3 - x4 =E= 0;

* set non default bounds

*
* Lower and upper bounds of variables are updated by SparsePOP
*

x1.up = 100; 
x2.up = 200; 
x3.up = 110; 
x4.up = 300; 
x5.up = 100; 
x6.up = 124; 
x7.up = 500; 
x8.up = 100; 
x9.up = 200; 

* x1.up = 1; 
* x2.lo = 199;
* x2.up = 200; 
* x3.up = 1; 
* x4.up = 101; 
* x5.up = 1; 
* x6.up = 101; 
* x7.up = 2; 
* x8.up = 1.0; 
* x9.up = 101; 

* set non default levels


* set non default marginals


Model m / all /;

m.limrow=0; m.limcol=0;

$if NOT '%gams.u1%' == '' $include '%gams.u1%'

Solve m using NLP minimizing objvar;
