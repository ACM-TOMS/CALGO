*  NLP written by GAMS Convert at 07/19/01 13:40:22
*  
*  Equation counts
*     Total       E       G       L       N       X
*        13      13       0       0       0       0
*  
*  Variable counts
*                 x       b       i     s1s     s2s      sc      si
*     Total    cont  binary integer    sos1    sos2   scont    sint
*        17      17       0       0       0       0       0       0
*  FX     0       0       0       0       0       0       0       0
*  
*  Nonzero counts
*     Total   const      NL     DLL
*        33      17      16       0
*
*  Solve m using NLP minimizing objvar;


Variables  objvar,x2,x3,x4,x5,x6,x7,x8,x9,x12,x13,x14,x15;

Positive Variables x2,x3,x4,x5,x6,x7,x8,x9,x12,x13,x14,x15;

Equations  e1,e2,e3,e4,e5,e6,e7,e8,e9,e10,e11;


e1.. x2^2 - 2*x2 + x3^2 - 2*x3 + x4^2 + x5^2 - objvar =E= 0;

e2..  - x4 + x6 =E= -0.5;

e3..  - x5 + x7 =E= -0.5;

e4..    x4 + x8 =E= 1.5;

e5..    x5 + x9 =E= 1.5;

e6.. x6*x12 =E= 0;

e7.. x7*x13 =E= 0;

e8.. x8*x14 =E= 0;

e9.. x9*x15 =E= 0;

e10..  - 2*x2 + 2*x4 - x12 + x14 =E= 0;

e11..  - 2*x3 + 2*x5 - x13 + x15 =E= 0;

* set non default bounds

x2.up = 203;
x3.up = 203;
x4.up = 1.5;
x5.up = 1.5; 
x6.up = 1.0; 
x7.up = 1.0; 
x8.up = 1.0; 
x9.up = 1.0; 
x12.up = 200; 
x13.up = 200; 
x14.up = 200; 
x15.up = 200; 

* set non default levels


* set non default marginals


Model m / all /;

m.limrow=0; m.limcol=0;

$if NOT '%gams.u1%' == '' $include '%gams.u1%'

Solve m using NLP minimizing objvar;
