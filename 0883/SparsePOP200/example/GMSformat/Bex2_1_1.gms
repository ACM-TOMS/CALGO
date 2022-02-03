*  NLP written by GAMS Convert at 07/19/01 13:39:27
*  
*  Equation counts
*     Total       E       G       L       N       X
*         2       1       0       1       0       0
*  
*  Variable counts
*                 x       b       i     s1s     s2s      sc      si
*     Total    cont  binary integer    sos1    sos2   scont    sint
*         6       6       0       0       0       0       0       0
*  FX     0       0       0       0       0       0       0       0
*  
*  Nonzero counts
*     Total   const      NL     DLL
*        11       6       5       0
*
*  Solve m using NLP minimizing objvar;


Variables  x1,x2,x3,x4,x5,objvar;

Positive Variables x1,x2,x3,x4,x5;

Equations  e1,e2;


e1..  42*x1 - 50*x1*x1 - 50*x2*x2 - 50*x3* x3 - 50*x4*x4 - 50*x5*x5 + 44*x2 + 45*x3 + 47*x4 + 47.5*x5 - objvar =E= 0;

e2..    20*x1 + 12*x2 + 11*x3 + 7*x4 + 4*x5 =L= 40;

* set non default bounds

x1.up = 1; 
x2.up = 1; 
x3.up = 1; 
x4.up = 1; 
x5.up = 1; 

* set non default levels

x1.l = 1; 
x2.l = 1; 
x4.l = 1; 

* set non default marginals


Model m / all /;

m.limrow=0; m.limcol=0;

$if NOT '%gams.u1%' == '' $include '%gams.u1%'

Solve m using NLP minimizing objvar;
