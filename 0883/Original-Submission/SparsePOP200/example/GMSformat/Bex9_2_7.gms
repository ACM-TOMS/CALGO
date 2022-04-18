*  NLP written by GAMS Convert at 07/19/01 13:40:22
*  
*  Equation counts
*     Total       E       G       L       N       X
*        10      10       0       0       0       0
*  
*  Variable counts
*                 x       b       i     s1s     s2s      sc      si
*     Total    cont  binary integer    sos1    sos2   scont    sint
*        11      11       0       0       0       0       0       0
*  FX     0       0       0       0       0       0       0       0
*  
*  Nonzero counts
*     Total   const      NL     DLL
*        28      18      10       0
*
*  Solve m using NLP minimizing objvar;


Variables  objvar,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11;

Positive Variables x2,x3,x4,x5,x6,x7,x8,x9,x10,x11;

Equations  e1,e2,e3,e4,e5,e6,e7,e8,e9,e10;


*e1.. (x2 - 5)*(x2 - 5) + (1 + 2*x3)*(1 + 2*x3) - objvar =E= 0;
e1.. x2^2 - 10*x2 + 4*x3^2 +4*x3 - objvar =E= -26;

e2..  - 3*x2 + x3 + x4 =E= -3;

e3..    x2 - 0.5*x3 + x5 =E= 4;

e4..    x2 + x3 + x6 =E= 7;

e5..  - x3 + x7 =E= 0;

e6.. x4*x8 =E= 0;

e7.. x5*x9 =E= 0;

e8.. x6*x10 =E= 0;

e9.. x7*x11 =E= 0;

e10..  - 1.5*x2 + 2*x3 + x8 - 0.5*x9 + x10 - x11 =E= 2;

* set non default bounds

x2.up = 7;
x3.up = 7;
x4.up = 7.5; 
x5.up = 7; 
x6.up = 7; 
x7.up = 7; 
x8.up = 20; 
x9.up = 20; 
x10.up = 20; 
x11.up = 20; 

* set non default levels


* set non default marginals


Model m / all /;

m.limrow=0; m.limcol=0;

$if NOT '%gams.u1%' == '' $include '%gams.u1%'

Solve m using NLP minimizing objvar;
