*  NLP written by GAMS Convert at 08/31/02 18:40:01
*  
*  Equation counts
*     Total       E       G       L       N       X       C
*        11       1       0      10       0       0       0
*  
*  Variable counts
*                 x       b       i     s1s     s2s      sc      si
*     Total    cont  binary integer    sos1    sos2   scont    sint
*        11      11       0       0       0       0       0       0
*  FX     0       0       0       0       0       0       0       0
*  
*  Nonzero counts
*     Total   const      NL     DLL
*        56      46      10       0
*
*  Solve m using NLP minimizing objvar;


Variables  x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,objvar;

Positive Variables x1,x2,x3,x4,x5,x6,x7,x8,x9,x10;

Equations  e1,e2,e3,e4,e5,e6,e7,e8,e9,e10,e11;


e1..  - 8*x1 - 6*x3 + 7*x4 - 7*x5 =L= 1;

e2..  - 6*x1 + 2*x2 - 3*x3 + 9*x4 - 3*x5 =L= 3;

e3..    6*x1 - 7*x3 - 8*x4 + 2*x5 =L= 5;

e4..  - x1 + x2 - 8*x3 - 7*x4 - 5*x5 =L= 4;

e5..    4*x1 - 7*x2 + 4*x3 + 5*x4 + x5 =L= 0;

e6..    5*x7 - 4*x8 + 9*x9 - 7*x10 =L= 0;

e7..    7*x6 + 4*x7 + 3*x8 + 7*x9 + 5*x10 =L= 7;

e8..    6*x6 + x7 - 8*x8 + 8*x9 =L= 3;

e9..  - 3*x6 + 2*x7 + 7*x8 + x10 =L= 6;

e10..  - 2*x6 - 3*x7 + 8*x8 + 5*x9 - 2*x10 =L= 2;

e11..  x1*x6 + 2*x1 + 3*x6 + x2*x7 - 4*x2 - x7 + x3*x8 + 8*x3 - 2*x8 + x4*x9 + 4*x4 - 4*x9 + x5*x10 + 9*x5 + 5*x10 - objvar =E= 0;

* set non default bounds

x1.up = 20; 
x2.up = 20; 
x3.up = 20; 
x4.up = 20; 
x5.up = 20; 
x6.up = 20; 
x7.up = 20; 
x8.up = 20; 
x9.up = 20; 
x10.up = 20; 

* set non default levels


* set non default marginals


Model m / all /;

m.limrow=0; m.limcol=0;

$if NOT '%gams.u1%' == '' $include '%gams.u1%'

Solve m using NLP minimizing objvar;
