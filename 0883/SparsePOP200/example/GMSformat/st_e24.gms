*  NLP written by GAMS Convert at 08/29/02 12:49:54
*  
*  Equation counts
*     Total       E       G       L       N       X       C
*         5       1       1       3       0       0       0
*  
*  Variable counts
*                 x       b       i     s1s     s2s      sc      si
*     Total    cont  binary integer    sos1    sos2   scont    sint
*         3       3       0       0       0       0       0       0
*  FX     0       0       0       0       0       0       0       0
*  
*  Nonzero counts
*     Total   const      NL     DLL
*        11       9       2       0
*
*  Solve m using NLP minimizing objvar;


Variables  x1,x2,objvar;

Positive Variables x1;

Equations  e1,e2,e3,e4,e5;


e1..    2*x1 + 3*x2 =G= 9;

e2..    3*x1 - x2 =L= 8;

e3..  - x1 + 2*x2 =L= 8;

e4..    x1 + 2*x2 =L= 12;

e5..  (5 + x1 - x2)*(x1 + x2 - 1) + x1 - objvar =E= 0;

* set non default bounds

x1.up = 4; 
x2.lo = 1; 
x2.up = 5; 

* set non default levels


* set non default marginals


Model m / all /;

m.limrow=0; m.limcol=0;

$if NOT '%gams.u1%' == '' $include '%gams.u1%'

Solve m using NLP minimizing objvar;