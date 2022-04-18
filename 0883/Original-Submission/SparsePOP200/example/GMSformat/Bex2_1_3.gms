*  NLP written by GAMS Convert at 07/19/01 13:39:28
*  
*  Equation counts
*     Total       E       G       L       N       X
*        10       1       0       9       0       0
*  
*  Variable counts
*                 x       b       i     s1s     s2s      sc      si
*     Total    cont  binary integer    sos1    sos2   scont    sint
*        14      14       0       0       0       0       0       0
*  FX     0       0       0       0       0       0       0       0
*  
*  Nonzero counts
*     Total   const      NL     DLL
*        41      37       4       0
*
*  Solve m using NLP minimizing objvar;


Variables  x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,objvar;

Positive Variables x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13;

Equations  e1,e2,e3,e4,e5,e6,e7,e8,e9,e10;


* e1..  - (5*x1 - 0.5*(10*x1*x1 + 10*x2*x2 + 10*x3*x3 + 10*x4*x4) + 5*x2 + 5*x3
*       + 5*x4) + x5 + x6 + x7 + x8 + x9 + x10 + x11 + x12 + x13 + objvar =E= 0;
e1..   5*x1 - 5*x1*x1 - 5*x2*x2 - 5*x3*x3 - 5*x4*x4 + 5*x2 + 5*x3 + 5*x4 - x5 - x6 - x7 - x8 - x9 - x10 - x11 - x12 - x13 - objvar =E= 0;

e2..    2*x1 + 2*x2 + x10 + x11 =L= 10;

e3..    2*x1 + 2*x3 + x10 + x12 =L= 10;

e4..    2*x2 + 2*x3 + x11 + x12 =L= 10;

e5..  - 8*x1 + x10 =L= 0;

e6..  - 8*x2 + x11 =L= 0;

e7..  - 8*x3 + x12 =L= 0;

e8..  - 2*x4 - x5 + x10 =L= 0;

e9..  - 2*x6 - x7 + x11 =L= 0;

e10..  - 2*x8 - x9 + x12 =L= 0;

* set non default bounds

x1.up = 1; 
x2.up = 1; 
x3.up = 1; 
x4.up = 1; 
x5.up = 1; 
x6.up = 1; 
x7.up = 1; 
x8.up = 1; 
x9.up = 1; 
x10.up = 8;
x11.up = 8; 
x12.up = 8; 
x13.up = 1; 

* set non default levels

x1.l = 1; 
x2.l = 1; 
x3.l = 1; 
x4.l = 1; 
x5.l = 1; 
x6.l = 1; 
x7.l = 1; 
x8.l = 1; 
x9.l = 1; 
x10.l = 3; 
x11.l = 3; 
x12.l = 3; 
x13.l = 1; 

* set non default marginals


Model m / all /;

m.limrow=0; m.limcol=0;

$if NOT '%gams.u1%' == '' $include '%gams.u1%'

Solve m using NLP minimizing objvar;
