*  NLP written by GAMS Convert at 07/19/01 13:39:28
*  
*  Equation counts
*     Total       E       G       L       N       X
*         3       1       0       2       0       0
*  
*  Variable counts
*                 x       b       i     s1s     s2s      sc      si
*     Total    cont  binary integer    sos1    sos2   scont    sint
*         7       7       0       0       0       0       0       0
*  FX     0       0       0       0       0       0       0       0
*  
*  Nonzero counts
*     Total   const      NL     DLL
*        15      10       5       0
*
*  Solve m using NLP minimizing objvar;


Variables  x1,x2,x3,x4,x5,x6,objvar;

Positive Variables x1,x2,x3,x4,x5,x6;

Equations  e1,e2,e3;


* e1..  - (-0.5*(x1*x1 + x2*x2 + x3*x3 + x4*x4 + x5*x5) - 10.5*x1 - 7.5*x2 - 3.5*
*     x3 - 2.5*x4  -1.5*x5) + 10*x6 + objvar =E= 0;
e1..  +0.5*x1*x1 +0.5*x2*x2 +0.5*x3*x3 +0.5*x4*x4 +0.5*x5*x5 +10.5*x1 +7.5*x2 +3.5*x3 +2.5*x4 +1.5*x5 +10*x6 + objvar =E= 0;

e2..    6*x1 + 3*x2 + 3*x3 + 2*x4 + x5 =L= 6.5;

e3..    10*x1 + 10*x3 + x6 =L= 20;

* set non default bounds

x1.up = 1; 
x2.up = 1; 
x3.up = 1; 
x4.up = 1; 
x5.up = 1; 
x6.up = 20; 

* set non default levels

x2.l = 1; 
x4.l = 1; 
x5.l = 1; 
x6.l = 20; 

* set non default marginals


Model m / all /;

m.limrow=0; m.limcol=0;

$if NOT '%gams.u1%' == '' $include '%gams.u1%'

Solve m using NLP minimizing objvar;
