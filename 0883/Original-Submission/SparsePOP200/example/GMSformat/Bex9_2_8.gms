
*  NLP written by GAMS Convert at 07/19/01 13:40:23
*  
*  Equation counts
*     Total       E       G       L       N       X
*         6       6       0       0       0       0
*  
*  Variable counts
*                 x       b       i     s1s     s2s      sc      si
*     Total    cont  binary integer    sos1    sos2   scont    sint
*         7       7       0       0       0       0       0       0
*  FX     2       2       0       0       0       0       0       0
*  
*  Nonzero counts
*     Total   const      NL     DLL
*        14       8       6       0
*
*  Solve m using NLP minimizing objvar;


Variables  objvar,x2,x3,x4,x5,x6,x7;

Positive Variables x2,x3,x4,x5,x6,x7;

Equations  e1,e2,e3,e4,e5,e6;


e1.. 3*x3 - 4*x2*x3 + 2*x2 - objvar =E= -1;

e2..  - x3 + x4 =E= 0;

e3..    x3 + x5 =E= 1;

e4.. x6*x4 =E= 0;

e5.. x7*x5 =E= 0;

e6..    4*x2 - x6 + x7 =E= 1;

* set non default bounds

x2.up = 1; 
x3.up = 20; 
x4.up = 20; 
x5.up = 20; 
x6.fx = 0; 
x7.fx = 0; 

* set non default levels


* set non default marginals


Model m / all /;

m.limrow=0; m.limcol=0;

$if NOT '%gams.u1%' == '' $include '%gams.u1%'

Solve m using NLP minimizing objvar;

