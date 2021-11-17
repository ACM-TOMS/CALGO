*  NLP written by GAMS Convert at 07/19/01 13:39:29
*  
*  Equation counts
*     Total       E       G       L       N       X
*         6       1       0       5       0       0
*  
*  Variable counts
*                 x       b       i     s1s     s2s      sc      si
*     Total    cont  binary integer    sos1    sos2   scont    sint
*         7       7       0       0       0       0       0       0
*  FX     0       0       0       0       0       0       0       0
*  
*  Nonzero counts
*     Total   const      NL     DLL
*        37      36       1       0
*
*  Solve m using NLP minimizing objvar;


Variables  x1,x2,x3,x4,x5,x6,objvar;

Positive Variables x1,x2,x3,x4,x5,x6;

Equations  e1,e2,e3,e4,e5,e6;


e1..  - 6.5*x1 + 0.5*x1*x1 + x2 + 2*x3 + 3*x4 + 2*x5 + x6 + objvar =E= 0;

e2..    x1 + 2*x2 + 8*x3 + x4 + 3*x5 + 5*x6 =L= 16;

e3..  - 8*x1 - 4*x2 - 2*x3 + 2*x4 + 4*x5 - x6 =L= -1;

e4..    2*x1 + 0.5*x2 + 0.2*x3 - 3*x4 - x5 - 4*x6 =L= 24;

e5..    0.2*x1 + 2*x2 + 0.1*x3 - 4*x4 + 2*x5 + 2*x6 =L= 12;

e6..  - 0.1*x1 - 0.5*x2 + 2*x3 + 5*x4 - 5*x5 + 3*x6 =L= 3;

* set non default bounds

x1.up = 1; 
x2.up = 8;
x3.up = 2;
x4.up = 1; 
x5.up = 1; 
x6.up = 2; 

* set non default levels

x2.l = 6; 
x4.l = 1; 
x5.l = 1; 

* set non default marginals


Model m / all /;

m.limrow=0; m.limcol=0;

$if NOT '%gams.u1%' == '' $include '%gams.u1%'

Solve m using NLP minimizing objvar;
