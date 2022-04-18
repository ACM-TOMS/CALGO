*  NLP written by GAMS Convert at 07/19/01 13:39:39
*  
*  Equation counts
*     Total       E       G       L       N       X
*        17      17       0       0       0       0
*  
*  Variable counts
*                 x       b       i     s1s     s2s      sc      si
*     Total    cont  binary integer    sos1    sos2   scont    sint
*        23      23       0       0       0       0       0       0
*  FX     0       0       0       0       0       0       0       0
*  
*  Nonzero counts
*     Total   const      NL     DLL
*        64      40      24       0
*
*  Solve m using NLP minimizing objvar;


Variables  x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22,objvar;

Positive Variables x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22;

Equations  e1,e2,e3,e4,e5,e6,e7,e8,e9,e10,e11,e12,e13,e14,e15,e16,e17;


e1..    x1 + x2 + x3 + x4 =E= 300;

e2..    x5 - x6 - x7 =E= 0;

e3..    x8 - x9 - x10 - x11 =E= 0;

e4..    x12 - x13 - x14 - x15 =E= 0;

e5..    x16 - x17 - x18 =E= 0;

e6.. x13*x21 + 0.333*x1 - x5 =E= 0;

e7.. x13*x22 - x8*x20 + 0.333*x1 =E= 0;

e8..  - x8*x19 + 0.333*x1 =E= 0;

e9..  - x12*x21 - 0.333*x2 =E= 0;

e10.. x9*x20 - x12*x22 + 0.333*x2 =E= 0;

e11.. x9*x19 + 0.333*x2 - x16 =E= 0;

e12.. x14*x21 + 0.333*x3 + x6 =E= 30;

e13.. x10*x20 + x14*x22 + 0.333*x3 =E= 50;

e14.. x10*x19 + 0.333*x3 + x17 =E= 30;

e15..    x19 + x20 =E= 1;

e16..    x21 + x22 =E= 1;

e17..  - 0.00432*x1 - 0.01517*x2 - 0.01517*x9 - 0.00432*x13 + objvar =E= 0.9979;

* set non default bounds

x1.up = 300; 
x2.up = 300; 
x3.up = 151; 
x4.up = 300; 
x5.up = 30; 
x6.up = 30; 
x7.up = 300; 
x8.up = 300; 
x9.up = 300; 
x10.up = 300; 
x11.up = 300; 
x12.up = 300; 
x13.up = 300; 
x14.up = 300; 
x15.up = 300; 
x16.up = 300; 
x17.up = 30; 
x18.up = 300; 
x19.up = 1; 
x20.up = 1; 
x21.up = 1; 
x22.up = 1; 

* set non default levels


* set non default marginals


Model m / all /;

m.limrow=0; m.limcol=0;

$if NOT '%gams.u1%' == '' $include '%gams.u1%'

Solve m using NLP minimizing objvar;
