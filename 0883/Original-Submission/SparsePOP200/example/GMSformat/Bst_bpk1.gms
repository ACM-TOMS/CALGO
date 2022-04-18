*  NLP written by GAMS Convert at 08/31/02 18:44:46
*  
*  Equation counts
*     Total       E       G       L       N       X       C
*         7       1       0       6       0       0       0
*  
*  Variable counts
*                 x       b       i     s1s     s2s      sc      si
*     Total    cont  binary integer    sos1    sos2   scont    sint
*         5       5       0       0       0       0       0       0
*  FX     0       0       0       0       0       0       0       0
*  
*  Nonzero counts
*     Total   const      NL     DLL
*        17      13       4       0
*
*  Solve m using NLP minimizing objvar;


Variables  x1,x2,x3,x4,objvar;

Positive Variables x1,x2,x3,x4;

Equations  e1,e2,e3,e4,e5,e6,e7;


e1..    x1 + 4*x2 =L= 8;

e2..    4*x1 + x2 =L= 12;

e3..    3*x1 + 4*x2 =L= 12;

e4..    2*x3 + x4 =L= 8;

e5..    x3 + 2*x4 =L= 8;

e6..    x3 + x4 =L= 5;

e7..  - (x1 - x1*x3 - x3 + x1*x4 + x2*x3 - x2 - x2*x4) + objvar =E= 0;

* set non default bounds


* set non default levels


* set non default marginals


Model m / all /;

m.limrow=0; m.limcol=0;

$if NOT '%gams.u1%' == '' $include '%gams.u1%'

Solve m using NLP minimizing objvar;
