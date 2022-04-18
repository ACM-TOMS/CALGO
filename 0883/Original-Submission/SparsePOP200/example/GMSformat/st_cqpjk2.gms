*  NLP written by GAMS Convert at 08/31/02 18:59:46
*  
*  Equation counts
*     Total       E       G       L       N       X       C
*         2       1       0       1       0       0       0
*  
*  Variable counts
*                 x       b       i     s1s     s2s      sc      si
*     Total    cont  binary integer    sos1    sos2   scont    sint
*         4       4       0       0       0       0       0       0
*  FX     0       0       0       0       0       0       0       0
*  
*  Nonzero counts
*     Total   const      NL     DLL
*         7       4       3       0
*
*  Solve m using NLP minimizing objvar;


Variables  x1,x2,x3,objvar;

Positive Variables x1,x2,x3;

Equations  e1,e2;


e1..    x1 + x2 + x3 =L= 10000000000;

e2..  9*x1*x1 - 15*x1 + 9*x2*x2 - 12*x2 + 9*x3*x3 - 9*x3 - objvar =E= 0;

* set non default bounds

x1.up = 1; 
x2.up = 1; 
x3.up = 1; 

* set non default levels


* set non default marginals


Model m / all /;

m.limrow=0; m.limcol=0;

$if NOT '%gams.u1%' == '' $include '%gams.u1%'

Solve m using NLP minimizing objvar;
