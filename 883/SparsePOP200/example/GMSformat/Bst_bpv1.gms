*  NLP written by GAMS Convert at 08/31/02 18:48:44
*  
*  Equation counts
*     Total       E       G       L       N       X       C
*         5       1       3       1       0       0       0
*  
*  Variable counts
*                 x       b       i     s1s     s2s      sc      si
*     Total    cont  binary integer    sos1    sos2   scont    sint
*         5       5       0       0       0       0       0       0
*  FX     0       0       0       0       0       0       0       0
*  
*  Nonzero counts
*     Total   const      NL     DLL
*        13       9       4       0
*
*  Solve m using NLP minimizing objvar;


Variables  x1,x2,x3,x4,objvar;

Positive Variables x1,x2,x3,x4;

Equations  e1,e2,e3,e4,e5;


e1..    x1 + 3*x2 =G= 30;

e2..    2*x1 + x2 =G= 20;

e3..  - 1.6667*x3 + x4 =G= 10;

e4..    x3 + x4 =L= 15;

e5..  - (x1*x3 + x2*x4) + objvar =E= 0;

* set non default bounds

x1.up = 27; 
x2.up = 16; 
x3.up = 10; 
x4.up = 10; 

* set non default levels


* set non default marginals


Model m / all /;

m.limrow=0; m.limcol=0;

$if NOT '%gams.u1%' == '' $include '%gams.u1%'

Solve m using NLP minimizing objvar;
