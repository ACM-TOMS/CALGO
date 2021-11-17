*  NLP written by GAMS Convert at 08/31/02 18:52:41
*  
*  Equation counts
*     Total       E       G       L       N       X       C
*         6       1       4       1       0       0       0
*  
*  Variable counts
*                 x       b       i     s1s     s2s      sc      si
*     Total    cont  binary integer    sos1    sos2   scont    sint
*         5       5       0       0       0       0       0       0
*  FX     0       0       0       0       0       0       0       0
*  
*  Nonzero counts
*     Total   const      NL     DLL
*        14      11       3       0
*
*  Solve m using NLP minimizing objvar;


Variables  x1,x2,x3,x4,objvar;

Positive Variables x1,x2,x3,x4;

Equations  e1,e2,e3,e4,e5,e6;


e1..  - 4*x1 - x2 =G= -12;

e2..    3*x1 - x2 =G= -1;

e3..    4*x3 - x4 =L= 12;

e4..  - x3 - x4 =G= -8;

e5..    4*x3 - x4 =G= -1;

e6..  - (x2*x3 + x2 + x3 - x1*x3) + objvar =E= 0;

* set non default bounds

x1.up = 4; 
x2.up = 4; 
x3.up = 5; 
x4.up = 5; 

* set non default levels


* set non default marginals


Model m / all /;

m.limrow=0; m.limcol=0;

$if NOT '%gams.u1%' == '' $include '%gams.u1%'

Solve m using NLP minimizing objvar;
