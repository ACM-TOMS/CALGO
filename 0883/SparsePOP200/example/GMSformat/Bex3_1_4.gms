*  NLP written by GAMS Convert at 07/19/01 13:39:32
*  
*  Equation counts
*     Total       E       G       L       N       X
*         4       1       1       2       0       0
*  
*  Variable counts
*                 x       b       i     s1s     s2s      sc      si
*     Total    cont  binary integer    sos1    sos2   scont    sint
*         4       4       0       0       0       0       0       0
*  FX     0       0       0       0       0       0       0       0
*  
*  Nonzero counts
*     Total   const      NL     DLL
*        12       9       3       0
*
*  Solve m using NLP minimizing objvar;


Variables  x1,x2,x3,objvar;

Positive Variables x1,x2,x3;

Equations  e1,e2,e3,e4;


e1..    2*x1 - x2 + x3 + objvar =E= 0;

e2.. 4*x1^2 - 4*x1*x2 + 4*x1*x3 + 2*x2^2 - 2*x2*x3 + 2*x3^2 - 20*x1 + 9*x2 - 13*x3 =G= -24;

e3..    x1 + x2 + x3 =L= 4;

e4..    3*x2 + x3 =L= 6;

* set non default bounds

x1.up = 2;
x2.up = 2;  
x3.up = 3; 

* set non default levels

x1.l = 0.5; 
x3.l = 3; 

* set non default marginals


Model m / all /;

m.limrow=0; m.limcol=0;

$if NOT '%gams.u1%' == '' $include '%gams.u1%'

Solve m using NLP minimizing objvar;

