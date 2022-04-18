*  NLP written by GAMS Convert at 08/31/02 19:07:50
*  
*  Equation counts
*     Total       E       G       L       N       X       C
*         3       1       0       2       0       0       0
*  
*  Variable counts
*                 x       b       i     s1s     s2s      sc      si
*     Total    cont  binary integer    sos1    sos2   scont    sint
*         3       3       0       0       0       0       0       0
*  FX     0       0       0       0       0       0       0       0
*  
*  Nonzero counts
*     Total   const      NL     DLL
*         7       5       2       0
*
*  Solve m using NLP minimizing objvar;


Variables  x1,x2,objvar;

Positive Variables x1,x2;

Equations  e1,e2,e3;


e1..  - 6*x1 + 8*x2 =L= 3;

e2..    3*x1 - x2 =L= 3;

e3..  x1*x2 - x1 - x2 - objvar =E= 0;

* set non default bounds

x1.up = 5; 
x2.up = 5; 

* set non default levels


* set non default marginals


Model m / all /;

m.limrow=0; m.limcol=0;

$if NOT '%gams.u1%' == '' $include '%gams.u1%'

Solve m using NLP minimizing objvar;
