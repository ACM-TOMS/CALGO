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
*        22      14       8       0
*
*  Solve m using NLP minimizing objvar;


Variables  x1,objvar,x3,x4,x5,x6,x7,x8,x9;

Positive Variables x1,x3,x4,x5,x6,x7,x8,x9;

Equations  e1,e2,e3,e4,e5,e6,e7,e8;


* e1.. (x3 - 3)*(x3 - 3) + (x1 - 2)*(x1 - 2) - objvar =E= 0;
e1.. x1*x1 - 4*x1 + x3*x3 - 6*x3 -objvar =E= -13; 

e2..    x1 - 2*x3 + x4 =E= 1;

e3..  - 2*x1 + x3 + x5 =E= 2;

e4..    2*x1 + x3 + x6 =E= 14;

e5.. x4*x7 =E= 0;

e6.. x5*x8 =E= 0;

e7.. x6*x9 =E= 0;

e8..    2*x1 + x7 - 2*x8 + 2*x9 =E= 10;

* set non default bounds

x1.up = 10; 
x3.up = 8; 
x4.up = 10; 
x5.up = 10; 
x6.up = 10; 
x7.up = 10; 
x8.up = 10; 
x9.up = 10; 

* set non default levels


* set non default marginals


Model m / all /;

m.limrow=0; m.limcol=0;

$if NOT '%gams.u1%' == '' $include '%gams.u1%'

Solve m using NLP minimizing objvar;
